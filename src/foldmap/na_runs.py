# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Nucleic-acid strands from stem to stem: single-stranded runs and the joints between pieces of a bent duplex,
routed around the elements like loops (straight, round an elbow, bowed or along the grid, in stages that keep off
the strands drawn before them), hairpin loops drawn as arcs over a stem's tip, internal-loop sides as bows."""

from __future__ import annotations

import math

import numpy as np

from .layout import NT_STEP, PlacedStem
from .na_ends import _from_beside
from .na_paths import (
    _BEAD,
    TEXT_TOLL,
    NAPath,
    _bow,
    _dot,
    _face,
    _length,
    _na_clear,
    _rounded,
    _spaced,
    _xy,
    graze,
)
from .route import _INFLATE, _clear, _curve, _Router, clip_len

_BOW_MAX = 8  # internal-loop sides up to this many nucleotides are drawn as bows, not routed
_JOINT_STRAIGHT = 1.5  # a joint between two pieces of a bent duplex shorter than this may be drawn straight
_JOINT_NEAR = 2.5  # ... and, when no route reaches it, one shorter than this (9C75 P2b -> P2a, 1.83: random-100 C24),
# or between pieces whose boxes touch
# ... or, farther apart, one shorter than this that leaves its face at most 60° off its axis (_JOINT_FORWARD), keeps off
# the stems' and domains' names and the protein loops, and runs into neither piece: the duplex reads as carrying on
# (random-100 C2b-4: 8VWQ's P26a -> P26b, 4.41 apart but 76° off its face, stays a curve)
_JOINT_REACH = 5.0
_JOINT_FORWARD = 0.5
# a joint no route reaches and no line or elbow may take is drawn as a smooth curve out of its face along the axis and
# into the partner's face along its axis, its control points this far out (page units), the first clear taken: the
# duplex bends round; random-100 C2b-6 (8VWQ's P26a <-> P26b, 73-76° off their faces; P124a's pieces turned back)
_JOINT_BENDS = (1.0, 1.5, 2.0, 0.5, 3.0)
# ... and when none of those is clear, the same swung sideways by this much at each end (page units): outward on its
# rail's side first (pieces turned back from each other, 8VWQ's P124a and 7Y7C's P103, meet past the corners of their
# bodies, which reach beyond the rails), then the other way (pieces side by side, 7Y7C's P138c1 and P138c2)
_JOINT_SWINGS = (0.0, 1.0, 2.0, -1.0, -2.0)
# the extra cost of a grid cell under a strand drawn before, for a route from the free cells beside a face (the last
# pass's, _from_beside): it runs on top of another strand only where nothing else routes, as a loop keeps off a stack
# gap's cells (route._GAP_TOLL); C2b-6 review minor 3: 7Y7C's a:520-521 ran 1.2 along a strand, through its nucleotide
_BESIDE_TOLL = 50.0
# a run keeps off every strand drawn before it (and every protein loop) by a detour at most this many times as
# long as its plain drawing, which keeps off the hairpin loops and free ends only. Measured on 24 campaign entries,
# crossing a strand only a cost: run-run crossings 250 -> 239, protein loops crossed 178 -> 175, no entry worse in
# either (bound 2: 245; unbounded: 250, with 4YBB, 2GIS and 4P9R worse)
_DETOUR = 3.0
# either stage's drawing is also at most this many nucleotide steps longer than the plain one per nucleotide on it
# (and per end): no perimeter loops (4YBB's
# P72->P71 was drawn 5.7 times its direct distance, 1GID's P15->P16 3.6). Measured on all 39 entries: 6 keeps every
# entry's protein-loop crossings at or below step 2's (total 158 -> 165; 4 and 5 put 5B43 and 4YBB above it), the
# longest run at 3.1 times its direct distance (was 5.7)
_DETOUR_STEPS = 6.0
# the A* budget of that first stage (route._MAX_NODES for any other route): a detour longer than the caps above is not
# taken, so there is no need to search far. 5000 (was 20000) keeps 4YBB's routing within 1.5x of step 2's time and
# changes no entry's crossings but 1Y26's run-run (3 -> 5)
_FIRST_NODES = 5_000
# failing that, a strand pays for the cells beside the strands (a run with nucleotides) and protein loops it passes
# as well, where its nucleotides would touch theirs, and that makes crossing them cost more than a cell; taken when
# at most this many times as long as its plain drawing. Measured on all 39 campaign entries with the first stage:
# run-run crossings 255 -> 237, protein loops crossed 188 -> 158, no entry's protein-loop crossings worse (1.25 and
# 1.5 with the beside cells in the first stage too: 4YBB and 5B43 worse; 3: 3CRX's strand exchange C:18-19 goes
# round a whole arm; without the cells beside protein loops: 4YBB 51 -> 58)
_BESIDE_DETOUR = 2.0


def _arc(d: PlacedStem, ka, kb, m: int) -> list[tuple[float, float]]:
    """A hairpin loop: a near-circle beyond the stem's tip from nucleotide ka to kb, room for m nucleotides."""
    pa, pb = d.point(ka), d.point(kb)
    out = d.direction if d.axial[ka] >= 0 else -d.direction
    centre, r = d.hairpin_circle(ka, kb, m)
    a0 = np.arctan2(*(pa - centre)[::-1])
    a1 = np.arctan2(*(pb - centre)[::-1])
    cands = [a1 + 2 * np.pi * k for k in (-1, 0, 1) if 0 < abs(a1 + 2 * np.pi * k - a0) < 2 * np.pi]
    a1 = max(cands, key=lambda c: _dot((np.cos((a0 + c) / 2), np.sin((a0 + c) / 2)), out))
    t = np.linspace(a0, a1, max(24, 4 * (m + 2)))
    return _xy(centre + r * np.c_[np.cos(t), np.sin(t)])


def _grazes(pts, own) -> bool:
    """Does the polyline run more than a bead's radius (_BEAD / 2) into either of its own pieces?"""
    return any(graze(p, q, d) > _BEAD / 2 for d in own for p, q in zip(pts, pts[1:]))


def _elbow(router: _Router, da: PlacedStem, db: PlacedStem, fa, ea, fb, eb, strands=None):
    """A joint the straight line cannot take (it runs along one of its own pieces): out of face fa along ea, then
    into face fb against its outward eb, a short way each (half, one, two nucleotide steps), joined straight or
    round one corner (where the two arms' lines meet ahead of both, or either axis-aligned corner). The first clear
    by _INFLATE of every element and stem but its own two, into neither of those beyond a bead's radius, and
    crossing no strand but its partner; None if none is."""
    fa, ea, fb, eb = (np.asarray(v, float) for v in (fa, ea, fb, eb))
    for h in (0.5 * NT_STEP, NT_STEP, 2 * NT_STEP):
        a, b = fa + ea * h, fb + eb * h
        corners = [[]]
        det = ea[0] * -eb[1] + eb[0] * ea[1]
        if abs(det) > 1e-9:  # a + ea t = b + eb u, both ahead
            t = ((b - a)[0] * -eb[1] + eb[0] * (b - a)[1]) / det
            u = (ea[0] * (b - a)[1] - ea[1] * (b - a)[0]) / det
            if t > 0 and u > 0:
                corners.append([a + ea * t])
        corners += [[np.array([a[0], b[1]])], [np.array([b[0], a[1]])]]
        for c in corners:
            pts = _xy([fa, a, *c, b, fb])
            if _grazes(pts, (da, db)) or not _clear(pts, (da, db), *router.elements, grow=_INFLATE, kink=True):
                continue
            if strands is None or not strands.crosses(pts, skip=frozenset((da.id, db.id))):
                return pts
    return None


def _bend(fa, ea, fb, eb, h: float, n: int = 24, sa=(0.0, 0.0), sb=(0.0, 0.0)) -> list[tuple[float, float]]:
    """The cubic from face fa (leaving along ea) to face fb (entered against its outward eb), its control points h
    out along each, and swung sideways by sa and sb."""
    p0, p3 = np.asarray(fa, float), np.asarray(fb, float)
    p1 = p0 + np.asarray(ea, float) * h + np.asarray(sa, float)
    p2 = p3 + np.asarray(eb, float) * h + np.asarray(sb, float)
    t = np.linspace(0.0, 1.0, n)[:, None]
    pts = (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3
    pts[0], pts[-1] = p0, p3
    return _xy(pts)


def _joint_bend(
    router: _Router, da: PlacedStem, db: PlacedStem, fa, ea, fb, eb, strands=None, names=(), loops=None, keys=()
):
    """A joint no route reaches and no straight line or elbow may take, drawn as a smooth curve (_bend, its control
    points _JOINT_BENDS out, then swung sideways _JOINT_SWINGS, outward on its rails' sides first: `keys`, the
    joined nucleotides; the first that is clear): clear by _INFLATE of every element and stem but its own two, into
    neither of those at all, crossing no strand drawn so far but its partner, off every name (`names`, as _off_names)
    and protein loop (`loops`); None if none is."""
    own = (da, db)
    side_a = da.normal * da.side[keys[0]] if keys else np.zeros(2)
    side_b = db.normal * db.side[keys[1]] if keys else np.zeros(2)
    for k, h in ((k, h) for k in _JOINT_SWINGS for h in _JOINT_BENDS):
        if k and not keys:
            break
        pts = _bend(fa, ea, fb, eb, h, sa=side_a * k, sb=side_b * k)
        if not _clear(pts, own, *router.elements, grow=_INFLATE, own_body=True):
            continue
        if strands is not None and strands.crosses(pts, skip=frozenset((da.id, db.id))):
            continue
        if not _off_names([pts[0], *pts, pts[-1]], names) or (loops is not None and loops.crosses(pts)):
            continue
        return pts
    return None


def _straight_joint(
    router: _Router, da: PlacedStem, db: PlacedStem, fa, ea, fb, strands=None, names=(), loops=None
) -> bool:
    """May a joint no route reaches be drawn straight from face fa to face fb, as a kink joins two pieces of one
    helix? Only going on forward (along fa's outward direction ea), clear by _INFLATE of every element and stem but its
    own two (a line that runs into one of those beyond a bead's radius is drawn as an _elbow instead where one is
    clear), and crossing no nucleic-acid strand drawn so far (`strands`, a _Strands) other than the partner strand of
    this same joint (two strands may cross where the duplex carries on, as in a helix). Between pieces that are near
    (faces closer than _JOINT_NEAR) or whose boxes touch; or farther apart, up to _JOINT_REACH, only when the line
    also leaves its face within 60° of its axis (_JOINT_FORWARD), runs into neither piece, and keeps off every name
    (`names`, as _off_names) and protein loop (`loops`, a _Strands). The near and touching joints keep the old test:
    held to the names and loops too, two more of 4YBB's became plain curves, which cross more than a line does
    (random-100 C2b-4 fix round 1)."""
    a, b = da.rect, db.rect
    touching = min(a[2], b[2]) - max(a[0], b[0]) >= -1e-9 and min(a[3], b[3]) - max(a[1], b[1]) >= -1e-9
    gap, forward = math.dist(fa, fb), _dot(np.subtract(fb, fa), ea)
    line = _xy([fa, fb])
    if forward <= 0:
        return False
    if not (gap < _JOINT_NEAR or touching):
        if gap >= _JOINT_REACH or forward < _JOINT_FORWARD * gap or _grazes(line, (da, db)):
            return False
        # _off_names judges a route between its faces (it drops the first and last point): doubled, the line itself
        if not _off_names([fa, fa, fb, fb], names) or (loops is not None and loops.crosses(line)):
            return False
    if not _clear(line, (da, db), *router.elements, grow=_INFLATE, kink=True):
        return False
    return strands is None or not strands.crosses(line, skip=frozenset((da.id, db.id)))


def _strand_between(
    router: _Router,
    s: int,
    i: int,
    j: int,
    da: PlacedStem,
    db: PlacedStem,
    kind,
    names=(),
    strands=None,
    loops=None,
    extra=None,
) -> NAPath:
    """The single strand from stem nucleotide i to stem nucleotide j (the m nucleotides between drawn on it): along
    its rail from i to its stem's end face, from that face to the next stem's end face, along the rail to j. Between
    the faces it is checked against every element and stem, its own two by their true bodies (_clear, own_body).
    Its plain drawing keeps off the hairpin loops and free ends drawn so far (`strands`, a _Strands), or, when none
    does, pays to cross them (_drawing levels 1, 0). When that crosses any strand (but a joint's partner) or protein
    loop (`loops`, a _Strands), two other drawings are tried, the first taken that is at most so many times as long
    as the plain one and keeps off the stem names if that does: for a run with nucleotides, one off every strand and
    protein loop, their grid cells blocked and no drawing crossing them taken (level 2; _DETOUR); then one that also
    pays for the cells in `extra`, beside strands and protein loops (_BESIDE_DETOUR). In each pass `names` (the
    stem-label boxes) are hard obstacles for the first choice of drawing: the first one clear that also keeps off
    every name; only when none does, the first one clear."""
    m, nts, between = j - i - 1, list(range(i + 1, j)), ((s, i), (s, j))
    own, need = (da, db), (j - i) * NT_STEP
    ra, rb = da.point((s, i)), db.point((s, j))
    if da is db:  # over the stem's tip: a hairpin loop
        pts = _arc(da, (s, i), (s, j), m)
        return NAPath(
            "hairpin", s, nts, pts, _spaced(pts, m), fallback=not _na_clear(pts, own, router.elements), between=between
        )
    (fa, ea), (fb, eb) = _face(da, (s, i), True), _face(db, (s, j), False)
    joint = m == 0 and da.parent and da.parent == db.parent
    args = (router, s, nts, between, da, db, (ra, fa, ea, fb, eb, rb), need, joint, kind, names, strands, loops)
    if strands is None:
        return _drawing(*args, 0)
    grid = router.grid
    start, counts = grid.used, (router.tight, router.rejected)

    def attempt(levels, cost):  # the first of these passes that draws, with this extra cost: (path, cells used, counts)
        router.tight, router.rejected = counts  # the debug counts of the attempt taken only
        for level in levels:
            grid.used = start + cost if cost is not None else start.copy()
            path = _drawing(*args, level)
            if path is not None:
                return path, (grid.used - cost if cost is not None else grid.used), (router.tight, router.rejected)
        return None, None, None

    plain, used, plain_counts = attempt((1, 0), None)  # off the hairpin loops and free ends, else crossing them
    skip = frozenset((da.id, db.id)) if joint else None
    grid.used = used
    if not strands.crosses(plain.points[1:-1], skip) and (loops is None or not loops.crosses(plain.points[1:-1])):
        return plain  # it crosses nothing: no other drawing can do better
    stages = [((2,), None, _DETOUR)] if nts else []  # off every strand and loop: a run with nucleotides only
    stages += [((1, 0), extra, _BESIDE_DETOUR)] if extra is not None else []  # paying for the cells beside them
    for levels, cost, bound in stages:
        path, marked, taken = attempt(levels, cost)
        if path is not None and _takes_over(path, plain, bound, len(nts), names):
            grid.used = marked
            router.tight, router.rejected = taken
            return path
    grid.used = used
    router.tight, router.rejected = plain_counts
    return plain


def _takes_over(path: NAPath, plain: NAPath, bound: float, m: int, names) -> bool:
    """May a later stage's drawing of a strand with m nucleotides replace its plain drawing? Only if it is no fallback
    where the plain one is clean, keeps off the stem names if the plain one does, and is no long detour: at most
    `bound` times as long as the plain drawing and at most _DETOUR_STEPS nucleotide steps per nucleotide (and end)
    longer."""
    longest = min(bound * _length(plain.points), _length(plain.points) + _DETOUR_STEPS * NT_STEP * (m + 2))
    return (
        _length(path.points) <= longest + 1e-9
        and (not path.fallback or plain.fallback)
        and (_off_names(path.points, names) or not _off_names(plain.points, names))
    )


def _off_names(pts, names) -> bool:
    """Does a run's drawing keep off every stem-label box (`names`) between its faces?"""
    mid = pts[1:-1]
    if len(mid) < 2:
        return True
    xs, ys = [p[0] for p in mid], [p[1] for p in mid]
    x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
    near = [b for b in names if b[0] < x1 and x0 < b[2] and b[1] < y1 and y0 < b[3]]  # the boxes it can meet
    for p, q in zip(mid, mid[1:]):
        lo_x, hi_x, lo_y, hi_y = min(p[0], q[0]), max(p[0], q[0]), min(p[1], q[1]), max(p[1], q[1])
        for b in near:
            if b[0] < hi_x and lo_x < b[2] and b[1] < hi_y and lo_y < b[3] and clip_len(p, q, b) > 1e-9:
                return False
    return True


def _drawing(router: _Router, s, nts, between, da, db, ends, need, joint, kind, names, strands, loops, level):
    """One pass of _strand_between: the strand as drawn, or None when no drawing clear of the elements keeps off
    what the pass makes hard: (level 2) every strand drawn so far and every protein loop; (1) the hairpin loops and
    free ends; (0) nothing, crossing them is only a cost."""
    hard = level > 0
    ra, fa, ea, fb, eb, rb = ends
    m, own = len(nts), (da, db)
    grid = router.grid
    on_rails = lambda mid: _xy([ra, *mid, rb])  # noqa: E731
    skip = frozenset((da.id, db.id)) if joint else None  # a joint's partner strand may cross it

    def ok(pts) -> bool:  # clear between the faces, and off what this pass makes hard
        if not _clear(pts[1:-1], own, *router.elements, own_body=True):
            return False
        if level == 2 and (strands.crosses(pts[1:-1], skip) or (loops is not None and loops.crosses(pts[1:-1]))):
            return False
        return not (level == 1 and strands.crosses(pts[1:-1], skip, fixed_only=True))

    def off_names(pts) -> bool:
        return _off_names(pts, names)

    reach = max(need - math.dist(ra, fa) - math.dist(fb, rb), 0.0)  # what the part between the faces must hold
    if joint and math.dist(fa, fb) < _JOINT_STRAIGHT and _dot(fb - fa, ea) > 0:
        pts = on_rails([fa, fb])  # the duplex just carries on
        if ok(pts):
            for c in grid.cells_along(fa, fb):
                grid.used[c] += 1
            return NAPath("joint", s, [], pts, [], between=between)
    bowed = None
    if kind == "internal" and m <= _BOW_MAX:
        # one side of an internal loop: bowed out from face to face on its own side, so the two sides never meet
        bowed = on_rails(_side_bow(da, between[0], fa, fb, reach))
        if not ok(bowed):
            bowed = None
        elif off_names(bowed):
            return NAPath("run", s, nts, bowed, _spaced(bowed, m), between=between)

    def search(find):  # a route search on the grid as this pass makes it (hard: the strands it keeps off blocked)
        if not hard:
            return find()
        saved, budget = grid.blocked.copy(), (router.max_nodes, router.wide_nodes)
        grid.blocked |= (strands.cells | (grid.used > 0)) if level == 2 else strands.fixed_cells
        if level == 2:  # a detour longer than _DETOUR times the second stage's is not taken: no need to search far
            router.max_nodes, router.wide_nodes = _FIRST_NODES, 0
        try:
            return find()
        finally:
            grid.blocked, (router.max_nodes, router.wide_nodes) = saved, budget

    route = search(lambda: _faces_route(router, fa, ea, da, fb, eb, db))
    bend = (router, da, db, fa, ea, fb, eb, strands, names, loops, between)
    line = grazing = None
    if route is None and joint and _straight_joint(router, da, db, fa, ea, fb, strands, names, loops):
        # pieces that touch (no margin is kept between them): the duplex carries on straight across their meeting
        # corner, as a kink joins two pieces of one helix; past one of them, not through it (a shallow graze of a
        # corner by the line leaving its face, as 4OO8 P6b -> P6a, is kept), round an elbow or a bend where one is
        # clear, else routed from beside its faces (random-100 C2b-6: 7Y7C's a:520-521 ran into P50a's corner)
        line = [fa, fb]
        if _grazes(line, (da, db)):
            line, grazing = _elbow(router, da, db, fa, ea, fb, eb, strands) or _joint_bend(*bend), line
    elif route is None and joint:  # no straight line may join them: a bend, when one is clear
        line = _joint_bend(*bend)
    if route is None and line is None and not hard:  # no way straight out of a face, even loose (the last pass:
        # not where the strands drawn before block the way): from the free cells beside it, off those strands
        route = _beside_off_strands(router, (fa, ea, da), (fb, eb, db), strands)
    # a grazing line only in the last pass, once the routes from beside its faces are tried
    line = line if line is not None or route is not None or hard else grazing
    if line is not None:
        for a, b in zip(line, line[1:]):
            for c in grid.cells_along(a, b):
                grid.used[c] += 1
        return NAPath("joint", s, [], on_rails(line), [], between=between)
    fallback = route is None
    if fallback and hard and bowed is None:
        return None
    if route is None:
        route = _curve(fa, ea, fb, eb)
    pts = on_rails(route)
    shapes = [pts]  # the drawings clear of the elements, in order of preference
    if not fallback:
        rounded = on_rails(_rounded(route))
        shapes = [rounded, pts] if ok(rounded) else [pts]
    if not fallback and _length(shapes[0]) < need:  # too short for its nucleotides: bowed out, if that is clear
        left = np.array([fa[1] - fb[1], fb[0] - fa[0]])
        bows = (on_rails(_bow(fa, fb, reach, out=left)), on_rails(_bow(fa, fb, reach, out=-left)))
        shapes = [b for b in bows if ok(b)] + shapes
    if hard:
        shapes = [q for q in shapes if not fallback and ok(q)]
        if not shapes and bowed is None:
            return None
    soft = [bowed] if bowed is not None else [shapes[0]]  # the first choice when every drawing meets a name
    pts = next((q for q in shapes if not fallback and ok(q) and off_names(q)), soft[0])
    kind_drawn = "joint" if joint else "run"
    return NAPath(kind_drawn, s, nts, pts, _spaced(pts, m), fallback=fallback and pts is not bowed, between=between)


def _beside_off_strands(router: _Router, a, b, strands=None):
    """_from_beside with every grid cell under a strand drawn so far (`strands`, a _Strands) costing _BESIDE_TOLL more
    to enter, so the route keeps off those strands unless nothing else routes."""
    grid = router.grid
    saved = grid.toll
    if strands is not None and strands.cells is not None:
        toll = strands.cells * _BESIDE_TOLL
        grid.toll = toll if saved is None else saved + toll
    try:
        return _from_beside(router, a, b)
    finally:
        grid.toll = saved


def _faces_route(router: _Router, fa, ea, da: PlacedStem, fb, eb, db: PlacedStem):
    """The route from end face fa to end face fb (router.route, its own two stems checked by their true bodies); when
    none starts or ends clear of the other obstacles' margins, again with its ends loose (_Grid.free_cell_along), the
    drawn route's own check deciding, as for a cropped gap (na_ends._gap): a face may look into the margin of a stem
    beside it, as two pieces of one bent duplex kept end to end do into each other's (random-100 C2b-4: 7Y7C's
    a:843-844 and 8VWQ's 1:836 were plain curves)."""
    route, _ = router.route(fa, ea, da, fb, eb, db, own_body=True, text_toll=TEXT_TOLL)
    grid = router.grid
    if route is None and (grid.free_cell_along(fa, ea, da) is None or grid.free_cell_along(fb, eb, db) is None):
        route, _ = router.route(fa, ea, da, fb, eb, db, own_body=True, loose=True, text_toll=TEXT_TOLL)
    return route


def _side_bow(d: PlacedStem, key, fa, fb, reach: float) -> list[tuple[float, float]]:
    """One side of an internal loop from end face fa to end face fb, bowed out on the side of its rail on stem d
    (at nucleotide `key`) until the part between the faces is `reach` long."""
    return _bow(fa, fb, reach, out=d.normal * d.side[key])
