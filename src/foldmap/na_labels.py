# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Text and marks beside the drawn nucleic-acid strands: contact beads kept apart along runs, a letter for each
single-stranded nucleotide, the "... N nt ..." labels of cropped gaps, isolated pairs (ties) as dotted lines or
labelled marks, and the page bounds that hold them all."""

from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np

from .layout import (
    _GLYPH_W,
    _LETTER_EM,
    HEAD,
    HEAD_LEN,
    NT_ROW,
    PITCH,
    SHAFT,
    Layout,
    _box_meets_body,
    _element_label_box,
    _footprint,
    _ink_points,
    _reach,
    glyph,
    label_boxes,
    line_clearance,
    stem_bodies,
    stem_name_keeps,
    stem_points,
    text_em,
)
from .na_paths import _BEAD, NAPath, NATie, _dot, _near, _seg_dist, _seg_dists, _spaced, crosses_body, in_body, midway
from .route import clip_len, clip_lengths

_TIE_SHORT = 3.0  # stem widths: a longer tie is drawn as two labelled marks, not a line
_TICK_GAP, _TICK = 0.15, 0.35  # a tie mark: a tick this far out from its nucleotide, this long
_TIE_TURNS = (0.0, 45.0, -45.0, 90.0, -90.0, 135.0, -135.0, 180.0)  # degrees from the partner a mark may point
_TICK_REACH = (1.0, 2.0, 3.0)  # times _TICK a mark's tick is drawn, the shortest that leaves its label clear


_SLIDES = tuple(f * k for k in (0.1, 0.2, 0.3, 0.4, 0.48) for f in (1, -1))  # spacings a run's nucleotides may slide


def _unbunch(layout: Layout, paths: list[NAPath]) -> None:
    """Slide each run's nucleotides along it, as little as keeps them a bead's width (_BEAD) from every nucleotide of
    another strand or stem, where its own spacing puts one on top of another's (two runs crossing or passing close);
    if no slide does, the one that keeps them furthest, and then each nucleotide still on another's spot slides alone
    within the room its neighbours leave (_slide_each). Never into its own stems' bodies."""
    at = {key: d for d in layout.stems for key in d.axial}
    taken = _Points(_BEAD)
    for d in layout.stems:
        for key in d.axial:
            taken.add(d.point(key))
    runs = [p for p in paths if p.kind == "run" and p.nts]
    for p in paths:
        if p.kind != "run":
            for xy in p.nt_xy:
                taken.add(xy)

    def room(xy, own) -> float:  # how near the nearest other nucleotide comes; -1 inside its own stems
        if any(in_body(q, d) for q in xy for d in own):
            return -1.0
        return min(taken.nearest(q) for q in xy)

    def room_one(xy, i, mine, own) -> float:  # room(): one nucleotide of a run, its run's others (mine) counted too
        if any(in_body(xy, d) for d in own):
            return -1.0
        near = [math.dist(xy, q) for k, q in enumerate(mine) if k != i]
        return min([taken.nearest(xy), *near])

    for p in runs:
        own = [at[key] for key in p.between]

        best, space, best_shift = p.nt_xy, room(p.nt_xy, own), 0.0
        for shift in _SLIDES if space < _BEAD else ():
            xy = _spaced(p.points, len(p.nts), shift)
            r = room(xy, own)
            if r > space:
                best, space, best_shift = xy, r, shift
            if r >= _BEAD:
                break
        if space < _BEAD:  # no slide of the whole run clears it: each nucleotide still on another's spot alone
            best = _slide_each(
                p.points, len(p.nts), best_shift, lambda xy, i, mine, own=own: room_one(xy, i, mine, own)
            )
        p.nt_xy = best
        for xy in best:
            taken.add(xy)
    _clear_junctions(layout, paths, runs, at)


def _clear_junctions(layout: Layout, paths: list[NAPath], runs: list[NAPath], at: dict) -> None:
    """Plan 28 F6 (M7): a run nucleotide _unbunch left under a bead's width from another run's (each run is slid
    against the runs placed before it, which stay put: 6QIK's A:G1447 and y:G85, 0.28 apart where three runs meet)
    has the other run slide instead: every run holding such a nucleotide, in turn, slides its nucleotides alone from
    where they are (_slide_from), against every other nucleotide as now placed. A slide is kept only when the run's
    nucleotides then fall short of a bead's width by less in all (a nucleotide's room is its nearest neighbour's
    distance, so a slide onto a third nucleotide counts against itself)."""
    fixed = [d.point(key) for d in layout.stems for key in d.axial]
    fixed += [xy for p in paths if p.kind != "run" for xy in p.nt_xy]
    every = _Points(_BEAD)  # every nucleotide as placed: the runs short of room are found once
    for xy in fixed:
        every.add(xy)
    whose = [-1] * len(fixed)
    for r, p in enumerate(runs):
        for xy in p.nt_xy:
            whose.append(r)
            every.add(xy)
    tight = sorted({whose[k] for k, xy in enumerate(every.xy) if whose[k] >= 0 and every.nearest(xy, skip=k) < _BEAD})
    for r in tight:
        p = runs[r]
        others = _Points(_BEAD)
        for xy in fixed:
            others.add(xy)
        for k, q in enumerate(runs):
            if k != r:
                for xy in q.nt_xy:
                    others.add(xy)
        own = [at[key] for key in p.between]

        def room(xy, i, mine, own=own, others=others) -> float:
            if any(in_body(xy, d) for d in own):
                return -1.0
            return min([others.nearest(xy), *(math.dist(xy, q) for k, q in enumerate(mine) if k != i)])

        def short(xy, room=room) -> float:  # how far the run's nucleotides fall short of a bead's width, together
            return sum(max(0.0, _BEAD - room(q, i, xy)) for i, q in enumerate(xy))

        was = short(p.nt_xy)
        moved = _slide_from(p.points, list(p.nt_xy), room)
        if short(moved) < was - 1e-9:
            p.nt_xy = moved


def _slide_from(pts, xy: list, room) -> list[tuple[float, float]]:
    """_slide_each from the nucleotides' present places along the run `pts`: each one whose room(xy, i, all) is under
    a bead's width slides alone the least distance that clears it (else to the place that leaves it most room), never
    nearer than a bead's width (along the run) to its neighbours or the run's two ends."""
    p, seg, cum, at = _run_line(pts)

    def where(q) -> float:  # the place along the run nearest q
        a, b = p[:-1], p[1:]
        d = b - a
        f = np.clip(((q - a) * d).sum(axis=1) / np.maximum((d * d).sum(axis=1), 1e-12), 0.0, 1.0)
        near = a + f[:, None] * d
        k = int(np.argmin(np.hypot(*(near - q).T)))
        return float(cum[k] + f[k] * seg[k])

    return _slide_loop([where(np.asarray(q, float)) for q in xy], xy, float(cum[-1]), at, room)


_EACH_STEP = 0.05  # page units: the steps a nucleotide still on another's spot slides alone along its run (_slide_each)


def _run_line(pts):
    """A run's polyline as (points, segment lengths, cumulative lengths, at): at(t) is the point t along it."""
    p = np.asarray(pts, float)
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])

    def at(t: float) -> tuple[float, float]:
        k = min(int(np.searchsorted(cum, t, side="right")) - 1, len(seg) - 1)
        f = 0.0 if seg[k] == 0 else (t - cum[k]) / seg[k]
        q = p[k] + f * (p[k + 1] - p[k])
        return float(q[0]), float(q[1])

    return p, seg, cum, at


def _slide_loop(ts: list[float], xy: list, total: float, at, room) -> list[tuple[float, float]]:
    """Each nucleotide (place ts[i] along the run, xy[i] on the page) whose room(xy, i, all) is under a bead's width
    slid alone the least distance that clears it (else the one that leaves it most room), never nearer than a bead's
    width (along the run) to its neighbours or the run's two ends (0 and `total`); _slide_each's and _slide_from's."""
    m = len(xy)
    for i in range(m):
        space = room(xy[i], i, xy)
        if space >= _BEAD:
            continue
        t0 = ts[i]
        lo = (ts[i - 1] if i > 0 else 0.0) + _BEAD
        hi = (ts[i + 1] if i < m - 1 else total) - _BEAD
        for k in range(1, int(max(t0 - lo, hi - t0, 0.0) / _EACH_STEP) + 1):
            for t in (t0 + k * _EACH_STEP, t0 - k * _EACH_STEP):
                if lo <= t <= hi:
                    q = at(t)
                    r = room(q, i, xy)
                    if r > space:
                        ts[i], xy[i], space = t, q, r
                    if r >= _BEAD:
                        break
            if space >= _BEAD:
                break
    return xy


def _slide_each(pts, m: int, shift: float, room) -> list[tuple[float, float]]:
    """The m nucleotides of a run spaced along `pts` (as _spaced with `shift`), each one whose room(xy, i, all) is under
    a bead's width slid alone along the run, the least distance that clears it (else the one that leaves it most
    room), never nearer than a bead's width (along the run) to its neighbours or the run's two ends: two runs side by
    side, spaced differently, whose nucleotides no slide of a whole run keeps apart (6S13 a:692-694 beside
    a:713-716)."""
    _, _, cum, at = _run_line(pts)
    total = float(cum[-1])
    ts = [float(t) for t in (np.arange(1, m + 1) + shift) / (m + 1) * total]
    return _slide_loop(ts, [at(t) for t in ts], total, at, room)


_LETTER_GAP = 0.38  # page units: closest two base letters (or a letter and a nucleotide) may sit
_LETTER_RINGS = tuple(NT_ROW * f for f in (3.8, 4.6, 5.5, 6.5))  # how far out a crowded letter looks last, nearest
_RING_TURNS = tuple(np.radians(a) for a in range(0, 360, 15))  # first, every 15° round its nucleotide


def _outline(e) -> list:
    """The parts of an element's drawn outline as bodies (cx, cy, length, width, direction, normal): a strand's shaft
    and its arrowhead (as render's _arrow draws them, the head taken as wide at its base all along), not its box, which
    is the strand pitch wide; any other element or stem, its body."""
    u = np.array([math.cos(e.angle), math.sin(e.angle)])
    n = np.array([-u[1], u[0]])
    if getattr(getattr(e, "sse", None), "kind", None) != "E":
        return [SimpleNamespace(cx=e.cx, cy=e.cy, length=e.length, width=e.width, angle=e.angle, direction=u, normal=n)]
    k = e.width / PITCH  # Style.strand_scale, carried by the strand's width
    head = min(HEAD_LEN, e.length * 0.45)
    shaft = SimpleNamespace(cx=e.cx - u[0] * head / 2, cy=e.cy - u[1] * head / 2, length=e.length - head)
    tip = SimpleNamespace(cx=e.cx + u[0] * (e.length - head) / 2, cy=e.cy + u[1] * (e.length - head) / 2, length=head)
    shaft.width, tip.width = SHAFT * k, HEAD * k
    for part in (shaft, tip):
        part.direction, part.normal, part.angle = u, n, e.angle
    return [shaft, tip]


def _in_sight(label, anchor, bodies, rects: np.ndarray) -> bool:
    """Does the straight line from a label at `label` to what it names at `anchor` (its nucleotide, its gap's middle,
    its end's nucleotide) cross no body's drawn outline (_outline) but the one the anchor lies under? A label set
    further out than beside its anchor (rings, outlines, far reaches) must not stand across an element from what it
    names: it would read as another's (review of held-out 2FD8: a letter left of strand L naming a nucleotide right
    of it). `rects`: the bodies' boxes (_boxes_of)."""
    (ax, ay), (bx, by) = label, anchor
    lo_x, hi_x, lo_y, hi_y = min(ax, bx), max(ax, bx), min(ay, by), max(ay, by)
    near = (rects[:, 0] <= hi_x) & (lo_x <= rects[:, 2]) & (rects[:, 1] <= hi_y) & (lo_y <= rects[:, 3])
    for k in np.flatnonzero(near).tolist():
        parts = _outline(bodies[k])
        if any(in_body(anchor, part) for part in parts):
            continue  # the body its anchor lies under
        if any(crosses_body(label, anchor, part) for part in parts):
            return False
    return True


def _boxes_of(bodies) -> np.ndarray:
    """Each body's (element's, stem's, turned) box (x0, y0, x1, y1), from its centre, length, width and angle."""
    if not bodies:
        return np.empty((0, 4))
    b = np.array([(e.cx, e.cy, e.length, e.width, e.angle) for e in bodies], float)
    c, s = np.abs(np.cos(b[:, 4])), np.abs(np.sin(b[:, 4]))
    hx, hy = c * b[:, 2] / 2 + s * b[:, 3] / 2, s * b[:, 2] / 2 + c * b[:, 3] / 2
    return np.stack([b[:, 0] - hx, b[:, 1] - hy, b[:, 0] + hx, b[:, 1] + hy], axis=1)


def _strand_lines(paths: list[NAPath]) -> tuple[np.ndarray, np.ndarray]:
    """Every strand's drawn line as segments (x0, y0, x1, y1), and the index in `paths` of the strand each is of."""
    lines = [(k, np.asarray(p.points, float).reshape(-1, 2)) for k, p in enumerate(paths) if len(p.points) > 1]
    if not lines:
        return np.empty((0, 4)), np.empty(0, int)
    segs = np.vstack([np.hstack([q[:-1], q[1:]]) for _, q in lines])
    return segs, np.concatenate([np.full(len(q) - 1, k) for k, q in lines])


def _own_nearest(q, k: int, segs: np.ndarray, of: np.ndarray) -> bool:
    """Is strand k's line (segments `segs` of strands `of`, _strand_lines) the nearest strand to the spot q, no other
    strand's line nearer it (a tie is its own)? A letter nearer another strand reads as that strand's. A strand with
    no line (one point) is nobody's nearest and takes no part: such a letter is judged as before."""
    a, b = segs[:, :2], segs[:, 2:]
    ab = b - a
    along = ((np.asarray(q, float) - a) * ab).sum(axis=1)
    length2 = (ab * ab).sum(axis=1)
    t = np.clip(np.divide(along, length2, out=np.zeros_like(along), where=length2 > 0), 0.0, 1.0)
    d = np.hypot(*(a + t[:, None] * ab - np.asarray(q, float)).T)
    mine = d[of == k]
    return mine.size == 0 or not bool((d[of != k] < mine.min()).any())


def _place_letters(
    layout: Layout, paths: list[NAPath], loops=(), only=None
) -> dict[tuple[int, int], tuple[float, float]]:
    """A letter for each single-stranded nucleotide beside its strand: on the side away from the stems it leaves
    and enters (outside a hairpin loop, in the row of the stem's letters along a free end), else the other, else
    further out on either side, or a little along the strand, whichever first keeps clear of other letters,
    nucleotides, end labels, elements, stem names, the strands drawn (runs, ends, gaps: a spot a line passes
    through is not clear, random-100 C14), the protein loops (`loops`, polylines: random-100 held-out 13FZ's B:2
    and B:8 on loops) and the elements' names. When none is clear (a strand lying over the protein that holds it:
    held-out 2FD8's three nucleotides across a jelly roll), the same spots judged by the elements' drawn outlines (a
    strand's shaft is narrower than its box: a letter fits between two shafts), then rings further out
    (_LETTER_RINGS), nearest first, each in sight of its nucleotide (_in_sight); when none of those is clear either,
    the one furthest from the other letters. A spot beside or along the strand is clear only where no other strand
    drawn is nearer it than its own (_own_nearest; plan 29 Task 8b: 22MZ's U57, beside run 6 but nearer run 2's line,
    on run 2's far side, read as run 2's).
    `only` (a set of (strand, nucleotide)): place those letters again, every other letter staying where
    layout.na_letters has it, clear of the ties' and gaps' labels placed since too (_letters_off_loops)."""
    where = stem_points(layout)  # on its curve where a bent duplex is drawn as one (na_curve)
    lines = [np.asarray(p.points, float).reshape(-1, 2) for p in paths if len(p.points) > 1]
    lines += [np.asarray(pts, float).reshape(-1, 2) for pts in loops if len(pts) > 1]
    strand_segs = np.vstack([np.hstack([q[:-1], q[1:]]) for q in lines]) if lines else np.empty((0, 4))
    seg_lo = np.minimum(strand_segs[:, :2], strand_segs[:, 2:])
    seg_hi = np.maximum(strand_segs[:, :2], strand_segs[:, 2:])
    em = _em(layout)
    line_clear = line_clearance(layout)  # figcheck.LINE_TEXT_CLEARANCE_PX at the page's text size
    ink = {b: _reach(_ink_points(b, False, _LETTER_EM * em, 0.0)) for b in "ACGUTIN"}

    def struck(q, key, clear: float | None = None) -> bool:  # a strand through the letter's ink, or near it
        base = layout.nucleic.strands[key[0]].bases[key[1]]
        x0, y0, x1, y1 = ink.get(base, ink["N"]) + np.array([q[0], q[1], q[0], q[1]])
        c = line_clear if clear is None else clear
        box = (x0 - c, y0 - c, x1 + c, y1 + c)
        near = (seg_lo[:, 0] < box[2]) & (seg_hi[:, 0] > box[0]) & (seg_lo[:, 1] < box[3]) & (seg_hi[:, 1] > box[1])
        return bool(near.any()) and bool((clip_lengths(strand_segs[near], [box]) > 0).any())

    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    taken = _Points(1.0)  # every letter, nucleotide and end label so far: nearer than 1 is measured exactly
    for xy in stem_points(layout, NT_ROW).values():  # each stem letter (on its curve where one is drawn)
        taken.add(xy)
    index = {key: taken.add(xy) for key, xy in where.items()}  # a letter's own nucleotide does not count
    for _, _, xy in layout.na_ends:
        taken.add(xy)
    stay = {} if only is None else {k: xy for k, xy in layout.na_letters.items() if k not in only}
    for xy in stay.values():
        taken.add(xy)
    bodies = [*layout.placed.values(), *stem_bodies(layout)]  # by their true (tilted) outlines, not their boxes
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    g = 1.2 * gh  # a stem's name stands a letter and a half clear of any base letter
    boxes = [np.add(b, (-g, -g, g, g)) for b in stem_name_keeps(layout, gw / _GLYPH_W, keep=False)]
    # and the elements' names, a letter's half size round them (held-out 13FZ: B:1's letter on A:987-990's name)
    half = (-gw / 2, -gh / 2, gw / 2, gh / 2)
    named = [p for p in layout.placed.values() if getattr(getattr(p, "sse", None), "kind", None) in ("H", "G")]
    boxes += [np.add(_element_label_box(p, gw / _GLYPH_W), half) for p in named]  # as label_boxes has them
    if only is not None:  # the ties' and gaps' labels, placed after the letters
        tie = (gw / 2 + 0.05, gh / 2)
        boxes += [
            (x - tie[0] - half[2], y - tie[1] - half[3], x + tie[0] + half[2], y + tie[1] + half[3])
            for t in layout.na_ties
            for _, _, (x, y) in t.marks
        ]
        boxes += [np.add(_gap_box(p, _em(layout)), half) for p in paths if p.kind == "gap" and p.label]
    out: dict[tuple[int, int], tuple[float, float]] = {}

    def room(q, own: int) -> float:
        if any(x0 < q[0] < x1 and y0 < q[1] < y1 for x0, y0, x1, y1 in boxes) or any(_near(q, e) for e in bodies):
            return 0.0
        return taken.nearest(q, own)

    rects = _boxes_of(bodies)

    def seen(q, xy) -> bool:  # a spot further out than beside it names its own nucleotide (_in_sight)
        return _in_sight(q, xy, bodies, rects)

    def outline_room(q, key, own: int) -> float:  # room, judged by the letter's ink (pad 0.05) on the drawn outlines
        if any(x0 < q[0] < x1 and y0 < q[1] < y1 for x0, y0, x1, y1 in boxes):
            return 0.0
        base = layout.nucleic.strands[key[0]].bases[key[1]]
        box = ink.get(base, ink["N"]) + np.array([q[0] - 0.05, q[1] - 0.05, q[0] + 0.05, q[1] + 0.05])
        near = (rects[:, 0] < box[2]) & (box[0] < rects[:, 2]) & (rects[:, 1] < box[3]) & (box[1] < rects[:, 3])
        if near.any() and any(_box_meets_body(box, part) for k in np.flatnonzero(near) for part in _outline(bodies[k])):
            return 0.0
        return taken.nearest(q, own)

    own_segs, own_of = _strand_lines(paths)

    def own(q, k: int) -> bool:  # no other strand's line nearer the spot than its own strand's (paths[k])
        return _own_nearest(q, k, own_segs, own_of)

    at = {key: d for d in layout.stems for key in d.axial}
    for k_path, p in enumerate(paths):
        pts = np.asarray(p.points, float)
        stems = [at[key] for key in p.between]
        centre = np.mean([[d.cx, d.cy] for d in stems], axis=0) if stems else None
        for m, (i, xy) in enumerate(zip(p.nts, p.nt_xy)):
            xy = np.asarray(xy)
            ahead = np.asarray(p.nt_xy[min(m + 1, len(p.nt_xy) - 1)]) - np.asarray(p.nt_xy[max(m - 1, 0)])
            if len(p.nt_xy) < 2:
                ahead = pts[-1] - pts[0]
            ahead = ahead / max(float(np.hypot(*ahead)), 1e-9)
            side = np.array([-ahead[1], ahead[0]])
            if centre is not None and _dot(side, xy - centre) < 0:
                side = -side
            options = [xy + side * NT_ROW * f for f in (1.0, -1.0, 1.7, -1.7, 2.4, -2.4)]
            options += [xy + side * NT_ROW * f + ahead * t for f in (1.0, -1.0) for t in (0.4, -0.4)]
            key = (p.strand, i)
            if key in stay:
                out[key] = stay[key]
                continue
            spaces = [room(q, index[key]) for q in options]
            free = [r >= _LETTER_GAP and not struck(q, key) and own(q, k_path) for q, r in zip(options, spaces)]
            if not any(free):  # crowded: further out, and further along the strand
                more = [xy + side * NT_ROW * f for f in (3.1, -3.1)]
                more += [xy + side * NT_ROW * f + ahead * t for f in (1.0, -1.0, 1.7, -1.7) for t in (0.8, -0.8)]
                options += more
                spaces += [room(q, index[key]) for q in more]
                free += [
                    r >= _LETTER_GAP and not struck(q, key) and own(q, k_path)
                    for q, r in zip(more, spaces[-len(more) :])
                ]
            if not any(free) and p.kind != "hairpin":  # by the drawn outlines, then on rings further out
                ok = [
                    outline_room(q, key, index[key]) >= _LETTER_GAP and not struck(q, key) and seen(q, xy)
                    for q in options
                ]
                if not any(ok):
                    ring = [xy + r * np.array([math.cos(t), math.sin(t)]) for r in _LETTER_RINGS for t in _RING_TURNS]
                    found = next(
                        (q for q in ring if room(q, index[key]) >= _LETTER_GAP and not struck(q, key) and seen(q, xy)),
                        None,
                    )
                    if found is not None:
                        options, spaces, ok = [*options, found], [*spaces, _LETTER_GAP], [*ok, True]
                free = ok if any(ok) else free
            if p.kind == "hairpin":  # a hairpin loop's letters stay in its stem's footprint (layout keeps it clear):
                # a spot inside clear with the line clearance, else one no strand runs through; else outside (8SA2)
                d = at[p.between[0]]
                fx0, fy0, fx1, fy1 = _footprint(d, gw / _GLYPH_W)
                g = gh / 2
                held = [
                    fx0 <= q[0] - d.cx - g
                    and q[0] - d.cx + g <= fx1
                    and fy0 <= q[1] - d.cy - g
                    and q[1] - d.cy + g <= fy1
                    for q in options
                ]
                inked = [h and r >= _LETTER_GAP and not struck(q, key, 0.0) for q, r, h in zip(options, spaces, held)]
                if any(ok and h for ok, h in zip(free, held)):  # a clear spot inside: only those
                    free = [ok and h for ok, h in zip(free, held)]
                elif any(inked):  # inside, no strand through its ink (within the clearance only)
                    free = inked
                elif not any(free) and any(held):  # none clear anywhere: the roomiest inside
                    spaces = [r if h else -1.0 for r, h in zip(spaces, held)]
            clear = [q for q, ok in zip(options, free) if ok]  # clear of letters and of every strand drawn
            spot = clear[0] if clear else options[int(np.argmax(spaces))]
            out[(p.strand, i)] = (float(spot[0]), float(spot[1]))
            taken.add(spot)
    return out


def _letters_off_loops(layout: Layout, loops) -> set[tuple[int, int]]:
    """In place: the base letters a protein loop (`loops`, polylines) still runs through, or nearer than the line
    clearance, once route._off_na_text has routed such loops again where it could, placed again clear of the loops
    (_place_letters, `only`); every other letter stays. Random-100 held-out 13FZ: B:2's and B:8's letters on loops
    no route could take off them. Letters are placed blind to the protein loops first so that a loop routed off a
    letter, as before, keeps that route (a letter kept off the loop from the start left the loop crossing the strand
    it would have been routed off: one more protein-loop crossing in 1IHF, 1QBJ, 1QTQ and 1URN). Returns the
    letters placed again."""
    lines = [np.asarray(pts, float).reshape(-1, 2) for pts in loops if len(pts) > 1]
    if not lines or not layout.na_letters:
        return set()
    segs = np.vstack([np.hstack([q[:-1], q[1:]]) for q in lines])
    em = _em(layout)
    c = line_clearance(layout)
    ink = {b: _reach(_ink_points(b, False, _LETTER_EM * em, 0.0)) for b in "ACGUTIN"}
    strands = layout.nucleic.strands
    hit = set()
    for (s, i), (x, y) in sorted(layout.na_letters.items()):
        x0, y0, x1, y1 = ink.get(strands[s].bases[i], ink["N"]) + np.array([x - c, y - c, x + c, y + c])
        near = (np.minimum(segs[:, 0], segs[:, 2]) < x1) & (np.maximum(segs[:, 0], segs[:, 2]) > x0)
        near &= (np.minimum(segs[:, 1], segs[:, 3]) < y1) & (np.maximum(segs[:, 1], segs[:, 3]) > y0)
        if near.any() and (clip_lengths(segs[near], [(x0, y0, x1, y1)]) > 0).any():
            hit.add((s, i))
    if hit:
        again = _place_letters(layout, layout.na_paths, [q.tolist() for q in lines], only=hit)
        layout.na_letters.update({k: again[k] for k in sorted(hit)})
    return hit


_GAP_ALONG = tuple(0.5 + s * k / 20 for k in range(9) for s in ((1,) if k == 0 else (-1, 1)))  # middle out
_GAP_OUT = (1.0, 2.0, 3.0, 4.0, 6.0)  # how far beside its gap a label may sit, in multiples of the least that clears
# it (4 and 6, random-100 C2b-6: 7Y7C's "... 62 nt ..." beside a gap 0.57 across, between stems)


_GAP_RINGS = tuple(1.0 + k / 4 for k in range(9))  # then round its midpoint, in multiples of its half width and half
_GAP_TURNS = tuple(math.radians(a) for a in range(0, 360, 5))  # height, every 5°, nearest first (random-100 held-out
# 8ZTU: "... 5 nt ..." by a gap 0.54 long among runs, where every spot beside it met one; the one clear spot, between
# two runs 0.7 apart, 2.7 such steps out)


_GAP_TEXT = 0.75  # a gap label's text size, in ems of the figure's text (render draws it at 0.75 of the font)
_GAP_CHAR = 0.5  # its mean character width, in ems of its own text ("…", digits and spaces are narrow)
_GAP_ROOM = 1.15  # room for the page to grow (its bounds widen after routing) and for a larger font


def _gap_half(label: str, em: float) -> tuple[float, float]:
    """Half the width and height of a gap label's box, one em of the figure's text being `em` page units
    (measured: 0.48 and 0.43 em per character and 1.0 and 0.94 em high on 2NZ4 and the 4YBB fixture)."""
    size = _GAP_TEXT * em * _GAP_ROOM
    return _GAP_CHAR * size * len(label) / 2 + 0.05, size / 2


def _em(layout: Layout) -> float:
    """Page units of one em of the figure's text on this layout's page (layout.text_em, 1 unit of pad each side)."""
    x0, _, x1, _ = layout.bounds
    return text_em(x1 - x0 + 2.0)


def _gap_box(p: NAPath, em: float) -> tuple[float, float, float, float]:
    """The box a gap's label takes (one em of text `em` page units): centred at label_xy, else standing on the
    gap's midpoint."""
    hw, hh = _gap_half(p.label, em)
    if p.label_xy is None:
        x, y = midway(p.points)
        return (x - hw, y, x + hw, y + 2 * hh)
    x, y = p.label_xy
    return (x - hw, y - hh, x + hw, y + hh)


def _along(pts, f: float):
    """(the point a fraction f of the way along the polyline, the unit direction there)."""
    p = np.asarray(pts, float)
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    t = f * cum[-1]
    k = min(int(np.searchsorted(cum, t, side="right")) - 1, len(seg) - 1)
    u = (p[k + 1] - p[k]) / seg[k] if seg[k] > 0 else np.array([1.0, 0.0])
    return p[k] + (t - cum[k]) * u, u


def _label_gaps(layout: Layout, paths: list[NAPath], loops=()) -> None:
    """Put each labelled gap's "... N nt ..." beside its dashed path (NAPath.label_xy), as base letters are placed:
    at a fraction of the way along it (_GAP_ALONG, the middle first), on either side, then further out, the
    first spot whose box meets nothing: no element's or stem's body, strand, tie or protein loop (`loops`, polylines),
    base letter, nucleotide (where its bead goes), end label, name or gap label placed before it, nor its own
    dashed line; when every spot meets something, the first spot round the gap's midpoint (_GAP_RINGS, nearest
    first) that meets nothing and is in sight of the gap (_in_sight) at its midpoint or either end (_gap_ring); when
    none does, the spot beside it meeting the fewest other things, then the least
    of its own line (a label beside its own gap reads as its own)."""
    gaps = [p for p in paths if p.kind == "gap" and p.label and len(p.points) > 1]
    if not gaps:
        return
    em = _em(layout)
    bodies = [*layout.placed.values(), *layout.ghosts, *stem_bodies(layout)]
    segs = [(a, b) for p in paths for a, b in zip(p.points, p.points[1:])]
    owner = [k for k, p in enumerate(paths) for _ in p.points[1:]]  # the path each segment is part of
    segs += [(a, b) for pts in loops for a, b in zip(pts, pts[1:])]
    # and the ties as drawn: a dotted line, or a tick at each end (random-100 C2b-6: 7Y7C's "... 62 nt ..." on a tick)
    segs += [(a, b) for t in layout.na_ties for a, b, _ in t.marks] + [
        tuple(t.ends) for t in layout.na_ties if not t.label
    ]
    owner += [-1] * (len(segs) - len(owner))
    seg_box = np.array([(min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1])) for a, b in segs])
    seg_box = seg_box.reshape(-1, 4)
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    g = gw / 2 + 0.05, gh / 2
    glyphs = list(stem_points(layout, NT_ROW).values()) + list(layout.na_letters.values())
    glyphs += [xy for _, _, xy in layout.na_ends]
    taken = [(x - g[0], y - g[1], x + g[0], y + g[1]) for x, y in glyphs] + [b for _, b in label_boxes(layout, [])]
    beads = list(stem_points(layout).values()) + [q for p in paths for q in p.nt_xy]
    beads = np.asarray(beads, float).reshape(-1, 2)

    body_box = _boxes_of(bodies)

    def meets(box, own: int) -> tuple[int, int]:  # (other things met, segments of its own gap met)
        x0, y0, x1, y1 = box
        edges = [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]
        under = (body_box[:, 0] < x1) & (x0 < body_box[:, 2]) & (body_box[:, 1] < y1) & (y0 < body_box[:, 3])
        n = sum(  # only a body whose box meets the box can meet it
            1
            for e in (bodies[k] for k in np.flatnonzero(under).tolist())
            if _near(((x0 + x1) / 2, (y0 + y1) / 2), e, 0.0)
            or x0 < e.cx < x1
            and y0 < e.cy < y1
            or any(crosses_body(a, b, e) for a, b in edges)
        )
        near = np.flatnonzero((seg_box[:, 0] < x1) & (x0 < seg_box[:, 2]) & (seg_box[:, 1] < y1) & (y0 < seg_box[:, 3]))
        hit = [k for k in near.tolist() if clip_len(*segs[k], box) > 1e-9]
        mine = sum(1 for k in hit if owner[k] == own)
        n += len(hit) - mine
        n += sum(1 for b in taken if b[0] < x1 and x0 < b[2] and b[1] < y1 and y0 < b[3])
        if len(beads):
            dx = beads[:, 0] - np.clip(beads[:, 0], x0, x1)
            dy = beads[:, 1] - np.clip(beads[:, 1], y0, y1)
            n += int(np.count_nonzero(dx * dx + dy * dy < (_BEAD / 2) ** 2))
        return n, mine

    for p in gaps:
        own = next(k for k, q in enumerate(paths) if q is p)
        hw, hh = _gap_half(p.label, em)
        best = None
        for out in _GAP_OUT:
            for f in _GAP_ALONG:
                q, u = _along(p.points, f)
                normal = np.array([-u[1], u[0]])
                reach = out * (abs(normal[0]) * hw + abs(normal[1]) * hh + 0.1)  # the box just clear of the gap
                for side in (1.0, -1.0):
                    x, y = q + side * reach * normal
                    box = (x - hw, y - hh, x + hw, y + hh)
                    n = meets(box, own)
                    if best is None or n < best[0]:
                        best = (n, (float(x), float(y)), box)
                    if n == (0, 0):
                        break
                if best[0] == (0, 0):
                    break
            if best[0] == (0, 0):
                break
        if best[0] != (0, 0):  # nowhere beside it: round its midpoint, nearest first, a spot that meets nothing
            ring = _gap_ring(p.points, hw, hh, lambda box, own=own: meets(box, own) == (0, 0), bodies, body_box)
            if ring is not None:
                best = ((0, 0), *ring)
        p.label_xy = best[1]
        taken.append(best[2])


def _gap_ring(points, hw: float, hh: float, clear, bodies, body_box):
    """(centre, box) of the first spot round a gap's midpoint (_GAP_RINGS times hw + hh out, every _GAP_TURNS step,
    nearest first) whose box (hw, hh: half extents) is clear(box) and that is in sight (_in_sight) of the gap's
    midpoint or of either end: a dashed line bent round a stem is read from any part of it (held-out 8ZTU's gap 499,
    whose one clear spot lay across a stem from its midpoint, beside its end); None if there is none."""
    mid = np.asarray(midway(points), float)
    anchors = (mid, np.asarray(points[0], float), np.asarray(points[-1], float))
    for f in _GAP_RINGS:
        for t in _GAP_TURNS:
            x, y = mid + f * (hw + hh) * np.array([math.cos(t), math.sin(t)])
            box = (x - hw, y - hh, x + hw, y + hh)
            if clear(box) and any(_in_sight((x, y), q, bodies, body_box) for q in anchors):
                return (float(x), float(y)), box
    return None


def _tie_name(n: int) -> str:
    """a, b, ..., z, aa, ab, ...: the label of the n-th tie drawn as marks."""
    return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n)


def _ties(layout: Layout, paths: list[NAPath], loops=()) -> list[NATie]:
    """Each isolated pair as drawn: a straight dotted line when it is short (at most _TIE_SHORT stem widths) and
    crosses no stem, element or letter, otherwise a tick at each nucleotide with the same label (a, b, ... in
    the order of the ties), each tick turned from its partner as little as keeps it and its label clear of
    bodies, letters, nucleotides, labels, strands and protein loops (`loops`, polylines); where no turn does, the
    tick is drawn longer (_TICK_REACH) and the turns tried again (random-100 C2b-4: 4P9R's and 7Y7C's labels took
    the first turn, on a strand or loop); where none does at any length, the one that meets the fewest lines, bodies
    and texts (random-100 C2b-6), the first such on a tie. A label is clear of a body when its box is (held-out 6S13:
    a label whose centre stood just off a stem had half its box on it)."""
    na = layout.nucleic
    where = stem_points(layout)  # on its curve where a bent duplex is drawn as one (na_curve)
    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    letters = stem_points(layout, NT_ROW)  # on its curve where a bent duplex is drawn as one (na_curve)
    letters.update({key: np.asarray(xy) for key, xy in layout.na_letters.items()})
    bodies = [*layout.placed.values(), *stem_bodies(layout)]
    outlines = [_outline(e) for e in bodies]  # as drawn: a label's box is clear of a strand's shaft and head
    rects = _boxes_of(bodies)

    def _under(box) -> list[int]:  # the bodies whose boxes the box meets: only those can be under it
        x0, y0, x1, y1 = box
        return np.flatnonzero(
            (rects[:, 0] < x1) & (x0 < rects[:, 2]) & (rects[:, 1] < y1) & (y0 < rects[:, 3])
        ).tolist()

    width = max((d.width for d in layout.stems), default=1.0)
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    boxes = [b for _, b in label_boxes(layout, [])] + [
        (x - gw / 2, y - gh / 2, x + gw / 2, y + gh / 2) for _, _, (x, y) in layout.na_ends
    ]
    segments = [(a, b) for p in paths for a, b in zip(p.points, p.points[1:])]
    segments += [(a, b) for pts in loops for a, b in zip(pts, pts[1:])]
    seg = np.asarray(segments, float).reshape(-1, 4)
    seg_lo, seg_hi = np.minimum(seg[:, :2], seg[:, 2:]), np.maximum(seg[:, :2], seg[:, 2:])
    taken = [*letters.values(), *where.values()]
    hw, hh = gw / 2 + 0.05, gh / 2

    def near_segments(box, b) -> list:  # only these can meet the box or come within 0.15 of b (protein loops: many)
        lo = np.minimum(box[:2], np.subtract(b, 0.15))
        hi = np.maximum(box[2:], np.add(b, 0.15))
        rows = np.flatnonzero(((seg_lo <= hi) & (seg_hi >= lo)).all(axis=1))
        return [segments[k] for k in rows.tolist()]

    def clear(tick, at, own) -> bool:  # own: the tie's nucleotide, which the tick starts beside
        (a, b), (x, y) = tick, at
        box = (x - hw, y - hh, x + hw, y + hh)
        if any(_box_meets_body(box, part) for k in _under(box) for part in outlines[k]):  # its box, not its centre
            return False
        if any(crosses_body(a, b, e) for e in bodies):
            return False
        pts = np.asarray(taken, float).reshape(-1, 2)  # vectorised (batch B: 37 s of 7Y7C's 140 s was here)
        near = pts[np.hypot(*(pts - np.asarray(own, float)).T) > 1e-9]
        if (np.hypot(*(near - np.asarray(at, float)).T) < _LETTER_GAP + hh).any():
            return False
        if (_seg_dists(near, a, b) < 0.2).any():
            return False
        if any(min(box[2], o[2]) > max(box[0], o[0]) and min(box[3], o[3]) > max(box[1], o[1]) for o in boxes):
            return False
        return not any(clip_len(p, q, box) > 1e-9 or _seg_dist(b, p, q) < 0.15 for p, q in near_segments(box, b))

    def misses(tick, at, own) -> tuple:  # how badly a mark that is not clear sits: what a reader sees first
        (a, b), (x, y) = tick, at
        box = (x - hw, y - hh, x + hw, y + hh)
        through = sum(1 for p, q in near_segments(box, b) if clip_len(p, q, box) > 1e-9)  # lines across its label
        on = sum(1 for k in _under(box) if any(_box_meets_body(box, part) for part in outlines[k])) + sum(  # text
            1 for o in boxes if min(box[2], o[2]) > max(box[0], o[0]) and min(box[3], o[3]) > max(box[1], o[1])
        )
        tick_on = sum(1 for e in bodies if crosses_body(a, b, e))
        pts = np.asarray(taken, float).reshape(-1, 2)
        near = pts[np.hypot(*(pts - np.asarray(own, float)).T) > 1e-9]
        crowd = int((np.hypot(*(near - np.asarray(at, float)).T) < _LETTER_GAP + hh).sum())
        return through + on, tick_on, crowd

    out: list[NATie] = []
    for p in na.structure.ties:
        ka, kb = (p.a, p.i), (p.b, p.j)
        if ka not in where or kb not in where:
            continue  # cropped away
        pa, pb = where[ka], where[kb]
        tie = NATie(ka, kb, [tuple(map(float, pa)), tuple(map(float, pb))])
        long = math.dist(pa, pb) > _TIE_SHORT * width
        hit = any(crosses_body(pa, pb, e) for e in bodies) or any(
            _seg_dist(xy, pa, pb) < gh / 2 for key, xy in letters.items() if key not in (ka, kb)
        )
        if long or hit:
            tie.label = _tie_name(sum(1 for t in out if t.label))
            for xy, other in ((pa, pb), (pb, pa)):
                toward = np.subtract(other, xy)
                toward = toward / max(float(np.hypot(*toward)), 1e-9)
                options = []
                for reach in _TICK_REACH:
                    for deg in _TIE_TURNS:
                        c, s_ = math.cos(math.radians(deg)), math.sin(math.radians(deg))
                        u = np.array([toward[0] * c - toward[1] * s_, toward[0] * s_ + toward[1] * c])
                        a, b = xy + u * _TICK_GAP, xy + u * (_TICK_GAP + _TICK * reach)
                        options.append(((a, b), b + u * (hh + 0.1)))
                found = next((o for o in options if clear(*o, xy)), None)
                if found is None:  # none clear: the one whose label the fewest lines, bodies and texts meet, then
                    # whose tick crosses fewest bodies, then nearest fewest nucleotides; the first on a tie
                    # (random-100 C2b-6: 7Y7C's and 8VWQ's crowded ties took the first turn, across runs)
                    found = min(options, key=lambda o: misses(*o, xy))
                (a, b), at = found
                tie.marks.append(((float(a[0]), float(a[1])), (float(b[0]), float(b[1])), (float(at[0]), float(at[1]))))
                boxes.append((at[0] - hw, at[1] - hh, at[0] + hw, at[1] + hh))
                taken.append(np.asarray(at))
        out.append(tie)
    return out


class _Points:
    """Points in a cell hash of cell size `reach`: the nearest one within reach is found among the 3 x 3 cells
    around a query, so placing n letters is linear in n, not quadratic."""

    def __init__(self, reach: float):
        self.reach = reach
        self.xy: list[tuple[float, float]] = []
        self.cells: dict[tuple[int, int], list[int]] = {}

    def _cell(self, x: float, y: float) -> tuple[int, int]:
        return math.floor(x / self.reach), math.floor(y / self.reach)

    def add(self, q) -> int:
        """Add a point; its index."""
        x, y = float(q[0]), float(q[1])
        self.xy.append((x, y))
        self.cells.setdefault(self._cell(x, y), []).append(len(self.xy) - 1)
        return len(self.xy) - 1

    def nearest(self, q, skip: int = -1) -> float:
        """Distance from q to the nearest point but `skip`, or `reach` when none is nearer."""
        x, y = float(q[0]), float(q[1])
        cx, cy = self._cell(x, y)
        best = self.reach
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for k in self.cells.get((cx + dx, cy + dy), ()):
                    if k != skip:
                        tx, ty = self.xy[k]
                        best = min(best, float(np.hypot(x - tx, y - ty)))
        return best


def _widen(layout: Layout, paths: list[NAPath]) -> None:
    """Grow layout.bounds to hold the strands, every base letter, the end labels, the gaps' labels and the
    ties' labels."""
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    pad = NT_ROW + max(gw, gh)
    boxes = [(x - pad, y - pad, x + pad, y + pad) for p in paths for x, y in p.points]
    g = max(gw, gh) / 2 + 0.05  # a letter turned any way
    letters = list(stem_points(layout, NT_ROW).values()) + list(layout.na_letters.values())
    letters += [xy for _, _, xy in layout.na_ends]
    boxes += [(x - g, y - g, x + g, y + g) for x, y in letters]
    boxes += [_gap_box(p, _em(layout)) for p in paths if p.label]
    boxes += [b for name, b in label_boxes(layout, []) if name.startswith("tie-label:")]
    if not boxes:
        return
    b = np.asarray(boxes, float)
    x0, y0, x1, y1 = layout.bounds
    layout.bounds = (
        min(x0, float(b[:, 0].min())),
        min(y0, float(b[:, 1].min())),
        max(x1, float(b[:, 2].max())),
        max(y1, float(b[:, 3].max())),
    )
