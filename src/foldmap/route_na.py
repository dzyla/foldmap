# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Nucleic-acid strands on the page: single strands routed like loops, hairpin arcs, free ends, cropped gaps."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np

from .layout import _GLYPH_H, _GLYPH_W, NT_ROW, NT_STEP, Layout, PlacedStem, _stem_label_box, label_boxes
from .route import (
    _DIRS,
    _INFLATE,
    _REUSE,
    _TURN,
    _clear,
    _corners,
    _curve,
    _element_obstacles,
    _Grid,
    _Router,
    _simplify,
    _snap_to_port,
    clip_len,
)

log = logging.getLogger(__name__)


_BOW_MAX = 8  # internal-loop sides up to this many nucleotides are drawn as bows, not routed
SHORT_GAP = 3  # a cropped gap of fewer nucleotides is drawn as a bare break, without its "… N nt …" label
_JOINT_STRAIGHT = 1.5  # a joint between two pieces of a bent duplex shorter than this may be drawn straight
_TAIL_ROUTED = 6  # a free end of more nucleotides that no straight drawing keeps clear is routed out to free space
# a run keeps off every strand drawn before it (and every protein loop) by a detour at most this many times as
# long as its plain drawing, which keeps off the hairpin loops and free ends only. Measured on 24 campaign entries,
# crossing a strand only a cost: run-run crossings 250 -> 239, protein loops crossed 178 -> 175, no entry worse in
# either (bound 2: 245; unbounded: 250, with 4YBB, 2GIS and 4P9R worse)
_DETOUR = 3.0
_FIRST_NODES = 20_000  # the A* budget of that first stage (route._MAX_NODES for any other route)
# failing that, a strand pays for the cells beside the strands (a run with nucleotides) and protein loops it passes
# as well, where its nucleotides would touch theirs, and that makes crossing them cost more than a cell; taken when
# at most this many times as long as its plain drawing. Measured on all 39 campaign entries with the first stage:
# run-run crossings 255 -> 237, protein loops crossed 188 -> 158, no entry's protein-loop crossings worse (1.25 and
# 1.5 with the beside cells in the first stage too: 4YBB and 5B43 worse; 3: 3CRX's strand exchange C:18-19 goes
# round a whole arm; without the cells beside protein loops: 4YBB 51 -> 58)
_BESIDE_DETOUR = 2.0


@dataclass(eq=False)
class NAPath:
    """A drawn stretch of single strand. kind: run (stem to stem), hairpin (an arc over a stem's tip), joint
    (between two pieces of a bent duplex), tail (a strand's free end), trace (a strand with no pairs) or gap
    (nucleotides cropped away, dashed; labelled "… N nt …" unless fewer than SHORT_GAP). nts: the nucleotide
    indices drawn on it, 5'->3', at nt_xy along `points`. fallback: no drawing clear of the elements was found.
    between: the (strand, nucleotide) of the stem nucleotides it leaves and enters (one for a tail)."""

    kind: str
    strand: int
    nts: list[int]
    points: list[tuple[float, float]]
    nt_xy: list[tuple[float, float]]
    label: str = ""
    fallback: bool = False
    between: tuple = ()


_TIE_SHORT = 3.0  # stem widths: a longer tie is drawn as two labelled marks, not a line
_TICK_GAP, _TICK = 0.15, 0.35  # a tie mark: a tick this far out from its nucleotide, this long
_TIE_TURNS = (0.0, 45.0, -45.0, 90.0, -90.0, 135.0, -135.0, 180.0)  # degrees from the partner a mark may point


@dataclass(eq=False)
class NATie:
    """An isolated base pair (a tie) as drawn, between nucleotides `a` and `b` at `ends`. Unlabelled: a thin
    dotted line (short, crossing no stem, element or letter). Labelled: at each end a tick (from, to) and the
    label's centre, `marks`, both ends carrying the same lowercase label."""

    a: tuple[int, int]
    b: tuple[int, int]
    ends: list[tuple[float, float]]
    label: str = ""
    marks: list = field(default_factory=list)


def path_name(path: NAPath, na) -> str:
    """How a warning names a drawn strand: "strand <chain>:<first>-<last>" by the nucleotides on it (or, with
    none, the stem nucleotides it joins; a gap by the strand piece before it)."""
    labels = na.strands[path.strand].labels
    idx = path.nts or [i for _, i in path.between] or [len(labels) - 1]
    first, last = labels[idx[0]], labels[idx[-1]]
    return f"strand {first.chain}:{first.tag}-{last.tag}"


def _use(grid: _Grid, pts, beside: np.ndarray | None = None) -> None:
    """Mark the cells under a drawn strand as used, so later routes avoid crossing it (at the reuse cost); and in
    `beside` the cells beside them, which a run with nucleotides on it pays to pass through as well (there its
    nucleotides would touch the strand's, and crossing the strand costs more than crossing a cell)."""
    under = {c for a, b in zip(pts, pts[1:]) for c in grid.cells_along(a, b)}
    for c in sorted(under):
        grid.used[c] += 1
    if beside is not None:
        near = {(x + dx, y + dy) for x, y in under for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))} - under
        for c in sorted(near):
            if 0 <= c[0] < grid.nx and 0 <= c[1] < grid.ny:
                beside[c] += 1


def _dot(a, b) -> float:
    """a . b of two 2-vectors, written out: numpy's dot goes through BLAS, whose kernels differ in the last bit."""
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1])


def _xy(pts) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in pts]


def _length(pts) -> float:
    p = np.asarray(pts, float)
    return float(np.hypot(*np.diff(p, axis=0).T).sum()) if len(p) > 1 else 0.0


def _spaced(pts, m: int, shift: float = 0.0) -> list[tuple[float, float]]:
    """m points evenly spaced along the polyline, ends excluded; `shift` (a fraction of their spacing, under a half)
    slides them all along it."""
    p = np.asarray(pts, float)
    if m <= 0 or len(p) < 2:
        return []
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    out = []
    for t in (np.arange(1, m + 1) + shift) / (m + 1) * cum[-1]:
        k = min(int(np.searchsorted(cum, t, side="right")) - 1, len(seg) - 1)
        f = 0.0 if seg[k] == 0 else (t - cum[k]) / seg[k]
        q = p[k] + f * (p[k + 1] - p[k])
        out.append((float(q[0]), float(q[1])))
    return out


_NA_CORNER = 0.6  # radius of a routed strand's rounded corners


def _rounded(pts, radius: float = _NA_CORNER, n: int = 8) -> list[tuple[float, float]]:
    """The polyline with each corner replaced by a quadratic curve (as a protein loop's corners are drawn),
    sampled, so the nucleotides spaced along it sit on the drawn line."""
    p = [np.asarray(q, float) for q in pts]
    out = [p[0]]
    t = np.linspace(0.0, 1.0, n)[:, None]
    for a, b, c in zip(p, p[1:], p[2:]):
        la, lc = float(np.hypot(*(b - a))), float(np.hypot(*(c - b)))
        r = min(radius, la / 2, lc / 2)
        if r < 1e-9:
            out.append(b)
            continue
        q0, q2 = b - (b - a) / la * r, b + (c - b) / lc * r
        out += list((1 - t) ** 2 * q0 + 2 * t * (1 - t) * b + t**2 * q2)
    out.append(p[-1])
    return _xy(_simplify(out))


def midway(pts) -> tuple[float, float]:
    """The point half way along a polyline (where a gap's label goes)."""
    return _spaced(pts, 1)[0] if len(pts) > 1 else tuple(pts[0])


def _bow(pa, pb, need: float, n: int = 32, out=None) -> list[tuple[float, float]]:
    """A curve from pa to pb bowed out (to the left, or towards `out`) until it is about `need` long."""
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    d = pb - pa
    chord = max(float(np.hypot(*d)), 1e-9)
    left = np.array([-d[1], d[0]]) / chord
    if out is not None and _dot(left, out) < 0:
        left = -left
    h = 0.0
    while True:
        ctrl = (pa + pb) / 2 + left * h
        t = np.linspace(0, 1, n)[:, None]
        pts = (1 - t) ** 2 * pa + 2 * t * (1 - t) * ctrl + t**2 * pb
        if _length(pts) >= need or h > 4 * need:
            return _xy(pts)
        h += 0.25


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


_TURNS = (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4)  # a free end's bends


def _free_end(d: PlacedStem, key, m: int, five: bool, turn: float = 0.0):
    """A strand's free end of m nucleotides, straight on from the stem and away from it (after its first
    nucleotide turned by `turn` radians): (points 5'->3', nucleotide positions 5'->3', where the end label goes)."""
    _, out = d.port(key, leaving=not five)
    rail = d.point(key)
    c, s = math.cos(turn), math.sin(turn)
    bent = np.array([out[0] * c - out[1] * s, out[0] * s + out[1] * c])
    first = rail + out * (NT_STEP / 2 + NT_STEP)
    xy = [first + bent * (k - 1) * NT_STEP for k in range(1, m + 1)]  # nearest the stem first
    tip = first + bent * (m * NT_STEP + NT_ROW)
    pts = [rail] + xy
    if five:
        pts, xy = pts[::-1], xy[::-1]
    return _xy(pts), _xy(xy), (float(tip[0]), float(tip[1]))


def _tail(
    router: _Router, layout: Layout, d: PlacedStem, s: int, key, nts: list[int], five: bool, labels, strands=None
):
    """A free end drawn straight on, or bent after its first nucleotide when that is what keeps it (and its end
    label) clear of the elements, and its label clear of loops and other text; one of more than _TAIL_ROUTED
    nucleotides that no such drawing keeps clear is routed out to free space (_routed_end), keeping off the strands
    drawn so far (`strands`); else it is drawn straight, a fallback. A drawing that keeps out of its partner's
    doorway (_doorway) comes first: (path, end label position, outward direction at its tip)."""
    found = None
    tries = (True, False) if labels is not None else (False,)
    door = _doorway(d, key, five)
    for keep_door in (True, False):
        for need_label in tries:
            for turn in _TURNS:
                pts, xy, tip = _free_end(d, key, len(nts), five, turn)
                drawn = [*(pts[::-1] if five else pts), tip]
                if (
                    _na_clear(drawn, (d,), router.elements)
                    and not (keep_door and door is not None and _crosses_segment(drawn[:-1], *door))
                    and (not need_label or _label_room(router, layout, tip, labels))
                ):
                    found = pts, xy, tip
                    break
            if found:
                break
        # routed: first off every cell a loop or strand uses; crossing them only as a cost when the doorway is given up
        for strict, need_label in [(True, n) for n in tries] + ([] if keep_door else [(False, n) for n in tries]):
            if found is not None or len(nts) <= _TAIL_ROUTED:
                break
            room = labels if need_label else None
            found = _routed_end(
                router, layout, d, key, len(nts), five, room, strands, door if keep_door else None, strict
            )
        if found:
            break
    pts, xy, tip = found or _free_end(d, key, len(nts), five)
    path = NAPath("tail", s, nts, pts, xy, fallback=found is None, between=(key,))
    far = xy[0] if five else xy[-1]
    v = np.subtract(tip, far)
    return path, tip, v / max(float(np.hypot(*v)), 1e-9)


def _doorway(d: PlacedStem, key, five: bool):
    """Where the partner of a free end's stem nucleotide leaves (or enters) the stem at the same end: the segment from
    its end face 3 nucleotide steps straight out, the way a strand from there to another stem starts. A free end
    turned across it would close it off. None if the nucleotide is unpaired."""
    partner = next(((b, j) for a, i, b, j in d.pairs if (a, i) == key), None)
    partner = partner or next(((a, i) for a, i, b, j in d.pairs if (b, j) == key), None)
    if partner is None or partner not in d.axial:
        return None
    face, out = _face(d, partner, five)
    return face, face + out * 3 * NT_STEP


def _crosses_segment(pts, a, b) -> bool:
    """Does the polyline cross segment a-b (properly, not merely touch)?"""

    def turn(o, u, v) -> float:
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])

    for p, q in zip(pts, pts[1:]):
        if turn(p, q, a) * turn(p, q, b) < -1e-12 and turn(a, b, p) * turn(a, b, q) < -1e-12:
            return True
    return False


def _along(pts, at: float) -> np.ndarray:
    """The point `at` along the polyline (its last point when it is shorter)."""
    p = np.asarray(pts, float)
    for a, b in zip(p, p[1:]):
        step = float(np.hypot(*(b - a)))
        if at <= step:
            return a + (b - a) * (at / step if step else 0.0)
        at -= step
    return p[-1]


def _upto(pts, at: float) -> list[tuple[float, float]]:
    """The polyline cut `at` along it."""
    p = np.asarray(pts, float)
    out = [p[0]]
    for a, b in zip(p, p[1:]):
        step = float(np.hypot(*(b - a)))
        if at <= step:
            out.append(a + (b - a) * (at / step if step else 0.0))
            return _xy(out)
        out.append(b)
        at -= step
    return _xy(out)


def _strung(line, gaps):
    """Marks strung along the polyline from its start, each the first point at exactly gaps[k] (measured straight)
    from the one before: (the marks, how far along the polyline the last but one lies), or None if it is too
    short."""
    p = np.asarray(line, float)
    marks: list[tuple[float, float]] = []
    prev, done, along, at = p[0], 0, 0.0, None
    for a, b in zip(p, p[1:]):
        v = b - a
        vv = _dot(v, v)
        t0 = 0.0
        while vv > 0:
            w = a - prev
            # the furthest root of |w + v t| = gap: where the segment leaves the circle around the last mark
            bq, cq = 2 * _dot(w, v), _dot(w, w) - gaps[done] ** 2
            disc = bq * bq - 4 * vv * cq
            t = (-bq + math.sqrt(disc)) / (2 * vv) if disc >= 0 else -1.0
            if not t0 <= t <= 1.0:
                break
            q = a + v * t
            marks.append((float(q[0]), float(q[1])))
            prev, done, t0 = q, done + 1, t
            if done == len(gaps) - 1:
                at = along + math.sqrt(vv) * t
            if done == len(gaps):
                return marks, at
        along += math.sqrt(vv)
    return None


def _reach_out(grid: _Grid, start, heading: int, steps: int, max_nodes: int = 50_000):
    """The cheapest orthogonal path of `steps` free cells from `start`, leaving along `heading` (an index into
    route._DIRS): a cell costs 1, a turn _TURN, a cell an earlier route uses _REUSE more, so it runs out into free
    space. Cells from start; None if no path is that long."""
    import heapq

    best = {(start, heading): 0.0}
    heap = [(0.0, 0, start, heading, 1)]
    parent: dict = {}
    tie, nodes = 0, 0
    while heap:
        g, _, cell, d, n = heapq.heappop(heap)
        if n >= steps:
            cells = [cell]
            key = (cell, d)
            while key in parent:
                key = parent[key]
                cells.append(key[0])
            return cells[::-1]
        if g > best.get((cell, d), math.inf):
            continue
        nodes += 1
        if nodes > max_nodes:
            return None
        for nd in range(4):
            if nd == (d + 2) % 4:
                continue
            nxt = (cell[0] + _DIRS[nd][0], cell[1] + _DIRS[nd][1])
            if not grid.free(nxt):
                continue
            cost = g + 1 + (_TURN if nd != d else 0) + (_REUSE if grid.used[nxt] else 0)
            if cost < best.get((nxt, nd), math.inf):
                best[(nxt, nd)] = cost
                parent[(nxt, nd)] = (cell, d)
                tie += 1
                heapq.heappush(heap, (cost, tie, nxt, nd, n + 1))
    return None


def _routed_end(
    router: _Router,
    layout: Layout,
    d: PlacedStem,
    key,
    m: int,
    five: bool,
    labels,
    strands=None,
    door=None,
    strict=False,
):
    """A long free end routed from its stem's end face out to free space along the grid (around the elements, off
    the cells of the strands drawn so far, `strands`), its m nucleotides spaced NT_STEP apart along it as on a
    straight end, checked like one: (points, nucleotide positions, end label position) as _free_end gives them, or
    None when no such drawing is clear (or, with `labels`, leaves its end label no room, or crosses `door`).
    `strict`: the cells any loop or strand drawn so far uses are blocked, not only paid for."""
    grid = router.grid
    rail = d.point(key)
    face, out = _face(d, key, not five)
    start = grid.free_cell_along(face, out, d)
    if start is None:
        return None
    # the first nucleotide 1.5 steps from the rail, the next ones a step apart, the end label beyond, as on a
    # straight end (measured straight, so a corner cannot bunch them up)
    gaps = [1.5 * NT_STEP] + [NT_STEP] * (m - 1) + [NT_STEP + NT_ROW]
    heading = min(range(4), key=lambda k: -_dot(_DIRS[k], out))
    saved = grid.blocked.copy()
    if strands is not None:
        grid.blocked |= strands.cells
    for c in grid.cells_along(*door) if door is not None else ():
        grid.blocked[c] = True
    if strict:
        grid.blocked |= grid.used > 0
    grid.blocked[start] = saved[start]
    try:
        cells = _reach_out(grid, start, heading, int(math.ceil(sum(gaps) / grid.res)) + 4)
    finally:
        grid.blocked = saved
    if cells is None:
        return None
    centres = [list(grid.center(c)) for c in _corners(cells)] if len(cells) > 1 else [list(grid.center(cells[0]))]
    _snap_to_port(centres, face, out, True)
    line = [_xy([rail])[0], *_rounded(_simplify([tuple(face), *[tuple(c) for c in centres]]))]
    strung = _strung(line, gaps)
    if strung is None:
        return None
    marks, last = strung
    if min(math.dist(p, q) for k, p in enumerate(marks) for q in marks[k + 1 :]) < 0.9 * NT_STEP:
        return None  # it turns back on itself: its nucleotides would bunch up
    xy, tip = marks[:-1], marks[-1]
    pts = _upto(line, last)
    if not _na_clear([*pts, tip], (d,), router.elements) or (strands is not None and strands.crosses(pts)):
        return None
    if door is not None and _crosses_segment(pts, *door):
        return None
    if labels is not None and not _label_room(router, layout, tip, labels):
        return None
    for c in cells:
        grid.used[c] += 1
    if five:
        pts, xy = pts[::-1], xy[::-1]
    return pts, xy, tip


_END_SPOTS = ((1.0, 0.0), (1.0, 0.6), (2.0, 0.0), (2.0, 0.6), (1.0, -0.6), (3.0, 0.0), (0.0, 2.4), (3.0, 1.2))


def _face_label(router: _Router, layout: Layout, d: PlacedStem, key, leaving: bool, labels) -> tuple:
    """Where the end label of a strand ending at a stem goes: just beyond its rail's end, else a little out to
    the side or further on, the first spot clear of loops, elements, stems and other text."""
    port, along = d.port(key, leaving)
    side = d.normal * d.side[key]  # away from the stem's axis
    spots = [port + along * NT_ROW * a + side * _GLYPH_H * b for a, b in _END_SPOTS]
    clear = (q for q in spots if labels is not None and _label_room(router, layout, q, labels))
    spot = next(clear, spots[0])
    return float(spot[0]), float(spot[1])


def _label_room(router: _Router, layout: Layout, xy, labels: list) -> bool:
    """An end label at xy keeps clear of the loops routed so far (their grid cells), every element's and stem's
    body, the text and stub boxes, and the end labels already placed."""
    hw, hh = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2
    x0, y0, x1, y1 = xy[0] - hw, xy[1] - hh, xy[0] + hw, xy[1] + hh
    grid = router.grid
    (i0, j0), (i1, j1) = grid.cell((x0, y0)), grid.cell((x1, y1))
    i0, j0, i1, j1 = max(i0, 0), max(j0, 0), min(i1, grid.nx - 1), min(j1, grid.ny - 1)
    if i0 <= i1 and j0 <= j1 and grid.used[i0 : i1 + 1, j0 : j1 + 1].any():
        return False
    if any(bx0 < x1 and x0 < bx1 and by0 < y1 and y0 < by1 for bx0, by0, bx1, by1 in router.keep_out):
        return False
    if any(abs(x - xy[0]) < 2 * hw and abs(y - xy[1]) < 2 * hh for x, y in labels):
        return False
    return not any(_near(xy, e, hh + 0.05) for e in [*layout.placed.values(), *layout.stems])


def _face(d: PlacedStem, key, leaving: bool) -> tuple[np.ndarray, np.ndarray]:
    """Where the strand crosses the stem's end face at this (end) nucleotide, along its rail, and the outward
    direction there: the port a strand is routed from, on the stem's outline as an element's port is."""
    _, out = d.port(key, leaving)
    reach = max(d.length / 2 - d.axial[key] * _dot(out, d.direction), 0.0)
    return d.point(key) + out * reach, out


def _local(pt, d) -> tuple[float, float]:
    c, s = math.cos(d.angle), math.sin(d.angle)
    x, y = pt[0] - d.cx, pt[1] - d.cy
    return x * c + y * s, y * c - x * s


def in_body(pt, d) -> bool:
    """Is pt strictly inside the element's (or stem's) tilted body?"""
    u, v = _local(pt, d)
    return abs(u) < d.length / 2 - 1e-9 and abs(v) < d.width / 2 - 1e-9


def _near(pt, e, margin: float = _GLYPH_W / 2) -> bool:
    """A letter at pt would touch the element's (or stem's) body."""
    u, v = _local(pt, e)
    return abs(u) < e.length / 2 + margin and abs(v) < e.width / 2 + margin


def _leaves_at(p, q, d) -> float:
    """The fraction of the way along p-q (p inside d's body) where it leaves the body."""
    (u0, v0), (u1, v1) = _local(p, d), _local(q, d)
    t = 1.0
    for a, b, h in ((u0, u1, d.length / 2), (v0, v1, d.width / 2)):
        if b > a:
            t = min(t, (h - a) / (b - a))
        elif b < a:
            t = min(t, (-h - a) / (b - a))
    return max(0.0, t)


def outside_own(pts, own) -> list[tuple[float, float]]:
    """The polyline without its first and last stretches inside its own stems' bodies (a strand runs along its
    rail to the stem's end face before it is drawn among the other elements)."""

    def trim(p):
        for i in range(len(p) - 1):
            mine = [d for d in own if in_body(p[i], d)]
            if not mine:
                return p[i:]
            if any(in_body(p[i + 1], d) for d in own):
                continue
            t = max(_leaves_at(p[i], p[i + 1], d) for d in mine)
            q = np.add(p[i], t * np.subtract(p[i + 1], p[i]))
            return [(float(q[0]), float(q[1]))] + list(p[i + 1 :])
        return [p[-1], p[-1]]

    out = trim(trim(list(pts))[::-1])[::-1]
    return out if len(out) > 1 else [out[0], out[0]]


def _na_clear(pts, own, elements) -> bool:
    """A strand drawn clear of every element and stem (by _clear), from where it leaves its own stems."""
    return _clear(outside_own(pts, own), own, *elements, own_body=True)


def route_nucleic(layout: Layout, router: _Router | None = None, loops=()) -> list[NAPath]:
    """Every single-stranded stretch as drawn: runs routed around the elements like loops (as are the joints
    between pieces of a bent duplex), hairpin loops as arcs, free ends straight on, strands without pairs along
    their projection, cropped gaps dashed. Every drawing is checked by _clear; one that fails is a fallback.
    `loops`: the protein loops already drawn (polylines), whose crossings a run weighs against crossing strands.
    Also fills layout.na_ends (where each strand's 5' and 3' labels go) and layout.na_letters (where each
    single-stranded nucleotide's letter goes), and widens layout.bounds to hold it all."""
    na = layout.nucleic
    if router is None:
        keep_out = [b for _, b in label_boxes(layout, [])]
        router = _Router(_Grid(layout, keep_out), _element_obstacles(layout), keep_out)
    elements = router.elements
    at = {key: d for d in layout.stems for key in d.axial}
    kind_of = {(r.strand, r.lo): r.kind for r in na.structure.runs}
    ends: list[tuple[int, str, tuple[float, float]]] = []
    exits: dict[tuple[int, int], tuple] = {}  # strand end -> (port, outward direction, owner) a gap leaves from

    def checked(path: NAPath, own) -> NAPath:
        path.fallback = path.fallback or not _na_clear(path.points, own, elements)
        return path

    gapped = {(g.before, "3′") for g in na.gaps} | {(g.after, "5′") for g in na.gaps}
    labels: list[tuple[float, float]] = []  # the end labels placed so far (a cropped end into a gap has none)

    def room(s: int, e: str):  # the labels an end's own label must keep clear of; None: it has no label
        return None if (s, e) in gapped else labels

    def end(s: int, e: str, xy) -> None:
        ends.append((s, e, xy))
        if (s, e) not in gapped:
            labels.append(xy)

    strands = _Strands(router.grid)  # the strands drawn so far: hard obstacles for a run's first routing pass
    crossed = _Strands(None)  # the protein loops, which a run keeps off the strands only if it crosses no more
    for pts in loops:
        crossed.add(pts)
    fixed: dict[tuple[int, int], NAPath] = {}  # (strand, nucleotide before it) -> a path of fixed shape

    def fix(key, path: NAPath) -> None:
        fixed[key] = path
        strands.add(among(path, at), fixed=path.kind in ("hairpin", "tail"))

    between: list[tuple] = []  # strands from stem to stem still to route, once the fixed shapes are known
    for s, strand in enumerate(na.strands):
        n = len(strand)
        if s in layout.na_traces:
            pts = _xy(layout.na_traces[s])
            fix((s, -1), checked(NAPath("trace", s, list(range(n)), pts, pts), ()))
            for label, a, b in (("5′", pts[0], pts[1]), ("3′", pts[-1], pts[-2])):
                v = np.subtract(a, b)
                v = v / max(float(np.hypot(*v)), 1e-9)
                end(s, label, tuple(float(c) for c in np.add(a, v * (NT_ROW + _GLYPH_W))))
                exits[(s, 0 if label == "5′" else n - 1)] = (np.add(a, v * NT_STEP / 2), v, None)
            continue
        held = [i for i in range(n) if (s, i) in at]
        if not held:
            continue
        first, last = held[0], held[-1]
        d0 = at[(s, first)]
        if first > 0:
            path, tip, v = _tail(router, layout, d0, s, (s, first), list(range(first)), True, room(s, "5′"), strands)
            fix((s, -1), path)
            exits[(s, 0)] = (np.add(path.nt_xy[0], v * NT_STEP / 2), v, d0)
        else:
            tip = _face_label(router, layout, d0, (s, 0), False, room(s, "5′"))
            exits[(s, 0)] = (*_face(d0, (s, 0), False), d0)
        end(s, "5′", tip)
        for i, j in zip(held, held[1:]):
            da, db = at[(s, i)], at[(s, j)]
            if da is db and da.side[(s, i)] == da.side[(s, j)]:
                continue  # along one rail of a stem (a bulge is drawn on the rail)
            if da is db:  # a hairpin loop: an arc over the tip
                fix((s, i), _strand_between(router, s, i, j, da, db, None))
            else:
                between.append((s, i, j, da, db, kind_of.get((s, i + 1))))
        dn = at[(s, last)]
        if last < n - 1:
            path, tip, v = _tail(
                router, layout, dn, s, (s, last), list(range(last + 1, n)), False, room(s, "3′"), strands
            )
            fix((s, last), path)
            exits[(s, n - 1)] = (np.add(path.nt_xy[-1], v * NT_STEP / 2), v, dn)
        else:
            tip = _face_label(router, layout, dn, (s, last), True, room(s, "3′"))
            exits[(s, n - 1)] = (*_face(dn, (s, last), True), dn)
        end(s, "3′", tip)
    # routes keep off the strands' end labels (as a protein loop keeps off its chain's N and C)
    g = _GLYPH_H / 2 + 0.05
    router.grid.keep_off([(x - g, y - g, x + g, y + g) for s, e, (x, y) in ends if (s, e) not in gapped])
    # and pay to cross a hairpin loop or a free end, as a loop pays to share an earlier loop's cells, or to run
    # beside one; beside a protein loop costs the same, so crossing one costs as much as crossing a strand
    beside, beside_loops = np.zeros_like(router.grid.used), np.zeros_like(router.grid.used)
    for p in fixed.values():
        if p.kind in ("hairpin", "tail"):
            _use(router.grid, p.points, beside)
    for pts in loops:
        _use(router.grid, pts, beside_loops)
    names = [_stem_label_box(d) for d in layout.stems]
    routed = {}
    for s, i, j, da, db, kind in between:
        # a run with nucleotides may pay for the cells beside strands and protein loops too, one without
        # nucleotides for those beside protein loops
        extra = beside + beside_loops if j - i > 1 else beside_loops
        p = routed[(s, i)] = _strand_between(router, s, i, j, da, db, kind, names, strands, crossed, extra)
        strands.add(among(p, at), frozenset((da.id, db.id)) if p.kind == "joint" else None)
        _use(router.grid, among(p, at), beside)
    out = [p for _, p in sorted({**fixed, **routed}.items())]
    for g in na.gaps:  # cropped away: the 3' end of one piece to the 5' end of the next, dashed
        a, b = exits.get((g.before, len(na.strands[g.before]) - 1)), exits.get((g.after, 0))
        if a is not None and b is not None:
            out.append(_gap(router, g, a, b))
    _unbunch(layout, out)
    layout.na_ends = [e for e in ends if (e[0], e[1]) not in gapped]
    layout.na_letters = _place_letters(layout, out)
    layout.na_ties = _ties(layout, out)
    _widen(layout, out)
    if log.isEnabledFor(logging.DEBUG):
        log.debug("%d nucleic-acid paths: %d fallbacks", len(out), sum(p.fallback for p in out))
    return out


class _Strands:
    """The nucleic-acid strands drawn so far, each as it lies among the elements (among): their segments, for an
    exact crossing test, and (given the grid) the grid cells under them, which a strand's first routing pass may not
    enter. A joint's segments carry the ids of the two pieces it joins (its partner strand may cross it). Without a
    grid it holds other lines to count crossings with (the protein loops)."""

    def __init__(self, grid: _Grid | None):
        self.grid = grid
        self.cells = np.zeros(grid.blocked.shape, bool) if grid is not None else None
        self.fixed_cells = np.zeros(grid.blocked.shape, bool) if grid is not None else None
        self._parts: list[tuple[np.ndarray, int, bool]] = []  # (segments, joint number (0: not a joint), fixed)
        self._joints: dict[frozenset, int] = {}  # the ids of two pieces -> their joint number
        self._all = None  # (segments, joint numbers) of every part, concatenated when first needed

    def add(self, pts, ids=None, fixed: bool = False) -> None:
        """Add a polyline: `ids`, the two pieces a joint joins; `fixed`, a hairpin loop or free end."""
        p = np.asarray(pts, float).reshape(-1, 2)
        if len(p) < 2:
            return
        joint = self._joints.setdefault(ids, len(self._joints) + 1) if ids is not None else 0
        self._parts.append((np.hstack([p[:-1], p[1:]]), joint, fixed))
        self._all = None
        for a, b in zip(pts, pts[1:]) if self.grid is not None else ():
            for c in self.grid.cells_along(a, b):
                self.cells[c] = True
                self.fixed_cells[c] |= fixed

    def crosses(self, pts, skip=None, fixed_only: bool = False) -> bool:
        """Does the polyline cross (properly, not merely touch) a strand drawn so far (a hairpin loop or free end
        only, with `fixed_only`), but a joint between the pieces `skip` (ids)?"""
        return self.crossings(pts, skip, fixed_only) > 0

    def crossings(self, pts, skip=None, fixed_only: bool = False) -> int:
        """How many (segment, strand segment) pairs cross properly, a joint between the pieces `skip` (ids) left
        out."""
        q = np.asarray(pts, float).reshape(-1, 2)
        if len(q) < 2 or not self._parts:
            return 0
        if self._all is None:
            segs = np.vstack([seg for seg, _, _ in self._parts])
            joints = np.concatenate([np.full(len(seg), j) for seg, j, _ in self._parts])
            self._all = segs, joints, np.concatenate([np.full(len(seg), f) for seg, _, f in self._parts])
        segs, joints, fixed = self._all
        keep = fixed.copy() if fixed_only else np.ones(len(segs), bool)
        if skip is not None and skip in self._joints:
            keep &= joints != self._joints[skip]
        x0, y0 = q.min(axis=0)
        x1, y1 = q.max(axis=0)
        keep &= (np.minimum(segs[:, 0], segs[:, 2]) <= x1) & (np.maximum(segs[:, 0], segs[:, 2]) >= x0)
        keep &= (np.minimum(segs[:, 1], segs[:, 3]) <= y1) & (np.maximum(segs[:, 1], segs[:, 3]) >= y0)
        other = segs[keep]
        if not len(other):
            return 0
        a = np.hstack([q[:-1], q[1:]])
        p0, p1 = a[:, None, 0:2], a[:, None, 2:4]
        u, v = other[None, :, 0:2], other[None, :, 2:4]

        def turn(o, e, w):  # written out elementwise, never through BLAS
            return (e[..., 0] - o[..., 0]) * (w[..., 1] - o[..., 1]) - (e[..., 1] - o[..., 1]) * (w[..., 0] - o[..., 0])

        crossed = (turn(p0, p1, u) * turn(p0, p1, v) < -1e-12) & (turn(u, v, p0) * turn(u, v, p1) < -1e-12)
        return int(crossed.sum())


_BEAD = 0.45  # page units: a contact bead's width (render: 1.05 x the backbone tube); nearer nucleotides' beads meet
_SLIDES = tuple(f * k for k in (0.1, 0.2, 0.3, 0.4, 0.48) for f in (1, -1))  # spacings a run's nucleotides may slide


def _unbunch(layout: Layout, paths: list[NAPath]) -> None:
    """Slide each run's nucleotides along it, as little as keeps them a bead's width (_BEAD) from every nucleotide of
    another strand or stem, where its own spacing puts one on top of another's (two runs crossing or passing close);
    if no slide does, the one that keeps them furthest. Never into its own stems' bodies."""
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

    for p in runs:
        own = [at[key] for key in p.between]

        best, space = p.nt_xy, room(p.nt_xy, own)
        for shift in _SLIDES if space < _BEAD else ():
            xy = _spaced(p.points, len(p.nts), shift)
            r = room(xy, own)
            if r > space:
                best, space = xy, r
            if r >= _BEAD:
                break
        p.nt_xy = best
        for xy in best:
            taken.add(xy)


def _straight_joint(router: _Router, da: PlacedStem, db: PlacedStem, fa, ea, fb, strands=None) -> bool:
    """May a joint no route reaches be drawn straight from face fa to face fb, as a kink joins two pieces of one
    helix? Only between pieces that are near (faces closer than _JOINT_STRAIGHT) or whose boxes touch, going on
    forward (along fa's outward direction ea), clear by _INFLATE of every element and stem but its own two, and
    crossing no nucleic-acid strand drawn so far (`strands`, a _Strands) other than the partner strand of this
    same joint (two strands may cross where the duplex carries on, as in a helix)."""
    a, b = da.rect, db.rect
    touching = min(a[2], b[2]) - max(a[0], b[0]) >= -1e-9 and min(a[3], b[3]) - max(a[1], b[1]) >= -1e-9
    if not (math.dist(fa, fb) < _JOINT_STRAIGHT or touching) or _dot(np.subtract(fb, fa), ea) <= 0:
        return False
    if not _clear(_xy([fa, fb]), (da, db), *router.elements, grow=_INFLATE, kink=True):
        return False
    return strands is None or not strands.crosses(_xy([fa, fb]), skip=frozenset((da.id, db.id)))


def among(path: NAPath, at: dict) -> list[tuple[float, float]]:
    """The part of a drawn strand that lies among the other elements. A run or a joint (stem to stem) is drawn as
    [rail point, end face, ..., end face, rail point]: without its two rail stretches, its first and last segments.
    Any other strand without its first and last stretches inside its own stems (outside_own). `at`: stem nucleotide
    -> its PlacedStem."""
    if path.kind in ("run", "joint") and len(path.between) == 2:
        return list(path.points[1:-1])
    return outside_own(path.points, tuple(at[k] for k in path.between if k in at))


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
    start = grid.used

    def attempt(levels, cost):  # the first of these passes that draws, with this extra cost: (path, cells used)
        for level in levels:
            grid.used = start + cost if cost is not None else start.copy()
            path = _drawing(*args, level)
            if path is not None:
                return path, (grid.used - cost if cost is not None else grid.used)
        return None, None

    plain, used = attempt((1, 0), None)  # off the hairpin loops and free ends, else crossing them only a cost
    skip = frozenset((da.id, db.id)) if joint else None
    grid.used = used
    if not strands.crosses(plain.points[1:-1], skip) and (loops is None or not loops.crosses(plain.points[1:-1])):
        return plain  # it crosses nothing: no other drawing can do better
    stages = [((2,), None, _DETOUR)] if nts else []  # off every strand and loop: a run with nucleotides only
    stages += [((1, 0), extra, _BESIDE_DETOUR)] if extra is not None else []  # paying for the cells beside them
    for levels, cost, bound in stages:
        path, marked = attempt(levels, cost)
        if (
            path is not None
            and _length(path.points) <= bound * _length(plain.points) + 1e-9
            and (_off_names(path.points, names) or not _off_names(plain.points, names))
        ):
            grid.used = marked
            return path
    grid.used = used
    return plain


def _off_names(pts, names) -> bool:
    """Does a run's drawing keep off every stem-label box (`names`) between its faces?"""
    return not any(clip_len(p, q, b) > 1e-9 for b in names for p, q in zip(pts[1:-1], pts[2:-1]))


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
    if hard:
        saved, budget = grid.blocked.copy(), router.max_nodes
        grid.blocked |= (strands.cells | (grid.used > 0)) if level == 2 else strands.fixed_cells
        if level == 2:  # a detour longer than _DETOUR times the second stage's is not taken: no need to search far
            router.max_nodes = _FIRST_NODES
        try:
            route, _ = router.route(fa, ea, da, fb, eb, db, own_body=True)
        finally:
            grid.blocked, router.max_nodes = saved, budget
    else:
        route, _ = router.route(fa, ea, da, fb, eb, db, own_body=True)
    if route is None and joint and _straight_joint(router, da, db, fa, ea, fb, strands):
        # pieces that touch (no margin is kept between them): the duplex carries on straight across their meeting
        # corner, as a kink joins two pieces of one helix
        for c in grid.cells_along(fa, fb):
            grid.used[c] += 1
        return NAPath("joint", s, [], on_rails([fa, fb]), [], between=between)
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


def _side_bow(d: PlacedStem, key, fa, fb, reach: float) -> list[tuple[float, float]]:
    """One side of an internal loop from end face fa to end face fb, bowed out on the side of its rail on stem d
    (at nucleotide `key`) until the part between the faces is `reach` long."""
    return _bow(fa, fb, reach, out=d.normal * d.side[key])


def _gap(router: _Router, g, a, b) -> NAPath:
    """A cropped gap from one strand piece's 3' end to the next one's 5' end: routed, dashed and (unless short)
    labelled with the number of nucleotides left out."""
    (pa, ea, own_a), (pb, eb, own_b) = a, b
    own = tuple(d for d in (own_a, own_b) if d is not None)
    route = None
    if own_a is not None and own_b is not None:
        route, _ = router.route(pa, ea, own_a, pb, eb, own_b, own_body=True)
        if route is None:  # an end in another element's clearance margin
            route, _ = router.route(pa, ea, own_a, pb, eb, own_b, own_body=True, loose=True)
    fallback = route is None
    pts = _xy(route if route is not None else [pa, pb])
    if not fallback:
        pts = next((q for q in (_rounded(pts),) if _na_clear(q, own, router.elements)), pts)
    label = f"… {g.count} nt …" if g.count >= SHORT_GAP else ""
    fallback = fallback or not _na_clear(pts, own, router.elements)
    return NAPath("gap", g.before, [], pts, [], label, fallback)


_LETTER_GAP = 0.38  # page units: closest two base letters (or a letter and a nucleotide) may sit


def _place_letters(layout: Layout, paths: list[NAPath]) -> dict[tuple[int, int], tuple[float, float]]:
    """A letter for each single-stranded nucleotide beside its strand: on the side away from the stems it leaves
    and enters (outside a hairpin loop, in the row of the stem's letters along a free end), else the other, else
    further out on either side, or a little along the strand, whichever first keeps clear of other letters,
    nucleotides, end labels, elements and stem names; when none is clear, the one furthest from the others."""
    where = {key: d.point(key) for d in layout.stems for key in d.axial}
    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    taken = _Points(1.0)  # every letter, nucleotide and end label so far: nearer than 1 is measured exactly
    for d in layout.stems:
        for key in d.axial:
            taken.add(d.point(key, NT_ROW))
    index = {key: taken.add(xy) for key, xy in where.items()}  # a letter's own nucleotide does not count
    for _, _, xy in layout.na_ends:
        taken.add(xy)
    bodies = [*layout.placed.values(), *layout.stems]  # by their true (tilted) outlines, not their boxes
    g = 1.2 * _GLYPH_H  # a stem's name stands a letter and a half clear of any base letter
    boxes = [np.add(_stem_label_box(d), (-g, -g, g, g)) for d in layout.stems]
    out: dict[tuple[int, int], tuple[float, float]] = {}

    def room(q, own: int) -> float:
        if any(x0 < q[0] < x1 and y0 < q[1] < y1 for x0, y0, x1, y1 in boxes) or any(_near(q, e) for e in bodies):
            return 0.0
        return taken.nearest(q, own)

    at = {key: d for d in layout.stems for key in d.axial}
    for p in paths:
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
            spaces = [room(q, index[(p.strand, i)]) for q in options]
            clear = [q for q, r in zip(options, spaces) if r >= _LETTER_GAP]
            spot = clear[0] if clear else options[int(np.argmax(spaces))]
            out[(p.strand, i)] = (float(spot[0]), float(spot[1]))
            taken.add(spot)
    return out


def _tie_name(n: int) -> str:
    """a, b, ..., z, aa, ab, ...: the label of the n-th tie drawn as marks."""
    return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n)


def _seg_dist(q, a, b) -> float:
    """Distance from point q to segment a-b."""
    q, a, b = (np.asarray(x, float) for x in (q, a, b))
    ab = b - a
    t = 0.0 if not ab.any() else float(np.clip(_dot(q - a, ab) / _dot(ab, ab), 0.0, 1.0))
    return float(np.hypot(*(q - a - t * ab)))


def crosses_body(p, q, e) -> bool:
    """Does segment p-q pass through the element's (or stem's) true, tilted outline?"""
    (u0, v0), (u1, v1) = _local(p, e), _local(q, e)
    return clip_len((u0, v0), (u1, v1), (-e.length / 2, -e.width / 2, e.length / 2, e.width / 2)) > 1e-9


def _ties(layout: Layout, paths: list[NAPath]) -> list[NATie]:
    """Each isolated pair as drawn: a straight dotted line when it is short (at most _TIE_SHORT stem widths) and
    crosses no stem, element or letter, otherwise a tick at each nucleotide with the same label (a, b, ... in
    the order of the ties), each tick turned from its partner as little as keeps it and its label clear of
    bodies, letters, nucleotides, labels and strands."""
    na = layout.nucleic
    where = {key: d.point(key) for d in layout.stems for key in d.axial}
    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    letters = {key: d.point(key, NT_ROW) for d in layout.stems for key in d.axial}
    letters.update({key: np.asarray(xy) for key, xy in layout.na_letters.items()})
    bodies = [*layout.placed.values(), *layout.stems]
    width = max((d.width for d in layout.stems), default=1.0)
    boxes = [b for _, b in label_boxes(layout, [])] + [
        (x - _GLYPH_W / 2, y - _GLYPH_H / 2, x + _GLYPH_W / 2, y + _GLYPH_H / 2) for _, _, (x, y) in layout.na_ends
    ]
    segments = [(a, b) for p in paths for a, b in zip(p.points, p.points[1:])]
    taken = [*letters.values(), *where.values()]
    hw, hh = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2

    def clear(tick, at, own) -> bool:  # own: the tie's nucleotide, which the tick starts beside
        (a, b), (x, y) = tick, at
        box = (x - hw, y - hh, x + hw, y + hh)
        if any(_near(at, e, 0.0) or crosses_body(a, b, e) for e in bodies):
            return False
        near = [t for t in taken if math.dist(t, own) > 1e-9]
        if any(math.dist(at, t) < _LETTER_GAP + hh for t in near) or any(_seg_dist(t, a, b) < 0.2 for t in near):
            return False
        if any(min(box[2], o[2]) > max(box[0], o[0]) and min(box[3], o[3]) > max(box[1], o[1]) for o in boxes):
            return False
        return not any(clip_len(p, q, box) > 1e-9 or _seg_dist(b, p, q) < 0.15 for p, q in segments)

    out: list[NATie] = []
    for p in na.structure.ties:
        ka, kb = (p.a, p.i), (p.b, p.j)
        if ka not in where or kb not in where:
            continue  # cropped away
        pa, pb = where[ka], where[kb]
        tie = NATie(ka, kb, [tuple(map(float, pa)), tuple(map(float, pb))])
        long = math.dist(pa, pb) > _TIE_SHORT * width
        hit = any(crosses_body(pa, pb, e) for e in bodies) or any(
            _seg_dist(xy, pa, pb) < _GLYPH_H / 2 for key, xy in letters.items() if key not in (ka, kb)
        )
        if long or hit:
            tie.label = _tie_name(sum(1 for t in out if t.label))
            for xy, other in ((pa, pb), (pb, pa)):
                toward = np.subtract(other, xy)
                toward = toward / max(float(np.hypot(*toward)), 1e-9)
                options = []
                for deg in _TIE_TURNS:
                    c, s_ = math.cos(math.radians(deg)), math.sin(math.radians(deg))
                    u = np.array([toward[0] * c - toward[1] * s_, toward[0] * s_ + toward[1] * c])
                    a, b = xy + u * _TICK_GAP, xy + u * (_TICK_GAP + _TICK)
                    options.append(((a, b), b + u * (hh + 0.1)))
                (a, b), at = next((o for o in options if clear(*o, xy)), options[0])
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
    pad = NT_ROW + max(_GLYPH_W, _GLYPH_H)
    boxes = [(x - pad, y - pad, x + pad, y + pad) for p in paths for x, y in p.points]
    g = max(_GLYPH_W, _GLYPH_H) / 2 + 0.05  # a letter turned any way
    letters = [d.point(key, NT_ROW) for d in layout.stems for key in d.axial] + list(layout.na_letters.values())
    letters += [xy for _, _, xy in layout.na_ends]
    boxes += [(x - g, y - g, x + g, y + g) for x, y in letters]
    for p in paths:
        if p.label:
            x, y = midway(p.points)
            half = _GLYPH_W * len(p.label) / 2 + 0.05
            boxes.append((x - half, y, x + half, y + _GLYPH_H))
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
