# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A nucleic-acid strand's free parts: free ends (tails) drawn straight on, bent or routed out to free space, where
the strands' 5' and 3' labels go, the clearance of strands without pairs drawn along their projection (traces), and
the dashed gaps where nucleotides are cropped away."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .layout import NT_ROW, NT_STEP, Layout, PlacedStem, glyph
from .na_paths import (
    _BEAD,
    TEXT_TOLL,
    NAPath,
    _bow,
    _crosses_segment,
    _dot,
    _face,
    _na_clear,
    _near,
    _rounded,
    _strung,
    _upto,
    _xy,
)
from .route import (
    _DIRS,
    _PUSH,
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
)

SHORT_GAP = 3  # a cropped gap of fewer nucleotides is drawn as a bare break, without its "… N nt …" label
_TAIL_ROUTED = 6  # a free end of more nucleotides that no straight drawing keeps clear is routed out to free space


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


def _enters(pts, circles, d: PlacedStem) -> bool:
    """Does the polyline come inside any of `circles` ((centre, radius, stem) of hairpin loops) but stem d's own?"""
    p = np.asarray(pts, float).reshape(-1, 2)
    a, b = (p[:-1], p[1:]) if len(p) > 1 else (p, p)
    ab = b - a
    for centre, r, owner in circles:
        if owner is d:
            continue
        t = np.clip(((centre - a) * ab).sum(axis=1) / np.maximum((ab * ab).sum(axis=1), 1e-12), 0.0, 1.0)
        if (np.hypot(*(a + ab * t[:, None] - centre).T) < r).any():
            return True
    return False


def _tail(
    router: _Router,
    layout: Layout,
    d: PlacedStem,
    s: int,
    key,
    nts: list[int],
    five: bool,
    labels,
    strands=None,
    circles=(),
    target=None,
):
    """A free end drawn straight on, or bent after its first nucleotide when that is what keeps it (and its end
    label) clear of the elements, and its label clear of loops and other text; one of more than _TAIL_ROUTED
    nucleotides that no such drawing keeps clear is routed out to free space (_routed_end), keeping off the strands
    drawn so far (`strands`); else it is drawn straight, a fallback. A drawing that keeps out of its partner's
    doorway (_doorway) comes first: (path, end label position, outward direction at its tip). No drawing comes into
    another stem's hairpin loop (`circles`: (centre, radius, stem) of every hairpin loop, _enters): a gap leaving
    the end from there would cross the loop (4YBB). A free end the protein touches tries its turns nearest its
    touching elements first (_toward), so it lies along what holds it (1URN's loop on U1A's RRM). Every drawing
    whose nucleotides keep a bead's width (_BEAD) off the strands drawn so far comes before any that does not (9CJI's
    3' end was turned onto its 5' end's line, nucleotides on top of each other): first the straight or turned ones
    off the strands, then the routed ones kept off their cells (strict); only then the old order, a routed end
    crossing strands as a cost last. `target`: where the molecule goes on
    past a gap at this 3' end, or comes from before a gap at this 5' end (_toward)."""
    found = None
    tries = (True, False) if labels is not None else (False,)
    door = _doorway(d, key, five)
    turns = _toward(layout, d, s, key, nts, five, target)
    for apart in (True, False) if strands is not None else (False,):  # off the strands drawn so far first
        for keep_door in (True, False):
            for need_label in tries:
                for turn in turns:
                    pts, xy, tip = _free_end(d, key, len(nts), five, turn)
                    drawn = [*(pts[::-1] if five else pts), tip]
                    if (
                        _na_clear(drawn, (d,), router.elements, both=False)  # its tip never back in its stem
                        and not _enters(drawn, circles, d)
                        and not (keep_door and door is not None and _crosses_segment(drawn[:-1], *door))
                        and (not need_label or _label_room(router, layout, tip, labels))
                        and not (apart and strands.near(xy, _BEAD))  # never along another strand (9CJI)
                    ):
                        found = pts, xy, tip
                        break
                if found:
                    break
            # routed: first off every cell a loop or strand uses; crossing them only as a cost when the doorway is
            # given up, and only after every drawing off the strands failed (a crossing is no better than running
            # by a strand: batch C2a review Q10)
            crossing = [] if keep_door or apart else [(False, n) for n in tries]
            for strict, need_label in [(True, n) for n in tries] + crossing:
                if found is not None or len(nts) <= _TAIL_ROUTED:
                    break
                room = labels if need_label else None
                found = _routed_end(
                    router, layout, d, key, len(nts), five, room, strands, door if keep_door else None, strict
                )
                if found is not None and _enters(found[3], circles, d):  # along its line out to its label
                    found = None
                found = found[:3] if found is not None else None
            if found:
                break
        if found:
            break
    pts, xy, tip = found or _free_end(d, key, len(nts), five)
    path = NAPath("tail", s, nts, pts, xy, fallback=found is None, between=(key,))
    far = xy[0] if five else xy[-1]
    v = np.subtract(tip, far)
    return path, tip, v / max(float(np.hypot(*v)), 1e-9)


def _touchers(layout: Layout, s: int, nts) -> dict[int, list[tuple[float, float, float, float]]]:
    """Nucleotide of strand s (of `nts`) -> the boxes of the placed elements touching it."""
    na, want = layout.nucleic, set(nts)
    out: dict[int, list] = {}
    for p in layout.placed.values():
        held = {k for r in range(p.sse.start, p.sse.end + 1) for t, k in na.contacts.get(r, ()) if t == s and k in want}
        for k in sorted(held):
            out.setdefault(k, []).append(p.rect)
    return out


def _toward(layout: Layout, d: PlacedStem, s: int, key, nts: list[int], five: bool, target=None) -> tuple[float, ...]:
    """_TURNS, of those at most a right angle off straight the ones that bring a free end's touched nucleotides
    nearest their touching elements' boxes (mean gap) first, then the rest. A free end at a gap (`target`) turns its
    tip nearest that first: a 3' end the protein does not touch, toward the stem point where the molecule goes on past
    the gap (3OWW: P5's loop is broken by a missing nucleotide; its 5' piece pointed away and the dashed gap went round
    the stem across a run); a 5' end, toward where the molecule comes from (the previous piece's 3' exit), weighed with
    its touching elements when the protein touches it (their mean gap plus the tip's distance; plan 28 F4: 4AS1's 5'
    end B:36-40 after the gap from B:31 turned to its protein, away from B:31, and the gap crossed it to reach its
    tip). Else _TURNS as they are."""
    held = _touchers(layout, s, nts)
    if held and not five:
        target = None  # a 3' end the protein touches lies along it (as before F4: 6QIK's gap after strand 4 crossed it)
    if not held and target is None:
        return _TURNS

    def gap(turn: float) -> float:
        tip = _free_end(d, key, len(nts), five, turn)
        cost = 0.0 if target is None else math.dist(tip[2], target)
        if held:
            xy = dict(zip(nts, tip[1]))
            cost += float(np.mean([_box_gap(xy[k], b) for k in sorted(held) for b in held[k]]))
        return cost

    ahead = sorted((t for t in _TURNS if abs(t) <= np.pi / 2 + 1e-9), key=lambda t: _snap_gap(gap(t)))
    return (*ahead, *(t for t in _TURNS if abs(t) > np.pi / 2 + 1e-9))  # folding back along the stem only after


def _box_gap(xy, b) -> float:
    """Distance from a point to a box (x0, y0, x1, y1), 0 inside it."""
    return float(math.hypot(max(b[0] - xy[0], xy[0] - b[2], 0.0), max(b[1] - xy[1], xy[1] - b[3], 0.0)))


def _snap_gap(v: float) -> float:
    return round(v, 6)  # ties (and their order, _TURNS') kept the same on every machine


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
    straight end, checked like one: (points, nucleotide positions, end label position) as _free_end gives them, and
    the line from the stem out to the label (what it is judged along), or None when no such drawing is clear (or,
    with `labels`, leaves its end label no room, or crosses `door`).
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
    marks, last, end = strung
    if min(math.dist(p, q) for k, p in enumerate(marks) for q in marks[k + 1 :]) < 0.9 * NT_STEP:
        return None  # it turns back on itself: its nucleotides would bunch up
    xy, tip = marks[:-1], marks[-1]
    pts = _upto(line, last)
    # judged along its line out to its label, as a straight end is: a chord from its last nucleotide to the label
    # cut the corner of an element the line keeps clear of (5LJ3's E:-16--5)
    out_to_label = _upto(line, end)
    if not _na_clear(out_to_label, (d,), router.elements, both=False) or (strands is not None and strands.crosses(pts)):
        return None
    if door is not None and _crosses_segment(pts, *door):
        return None
    if labels is not None and not _label_room(router, layout, tip, labels):
        return None
    for c in cells:
        grid.used[c] += 1
    if five:
        pts, xy = pts[::-1], xy[::-1]
    return pts, xy, tip, out_to_label


_END_SPOTS = ((1.0, 0.0), (1.0, 0.6), (2.0, 0.0), (2.0, 0.6), (1.0, -0.6), (3.0, 0.0), (0.0, 2.4), (3.0, 1.2))


def _face_label(router: _Router, layout: Layout, d: PlacedStem, key, leaving: bool, labels) -> tuple:
    """Where the end label of a strand ending at a stem goes: just beyond its rail's end, else a little out to
    the side or further on, the first spot clear of loops, elements, stems and other text."""
    port, along = d.port(key, leaving)
    side = d.normal * d.side[key]  # away from the stem's axis
    spots = [port + along * NT_ROW * a + side * glyph(layout)[1] * b for a, b in _END_SPOTS]
    clear = (q for q in spots if labels is not None and _label_room(router, layout, q, labels))
    spot = next(clear, spots[0])
    return float(spot[0]), float(spot[1])


def _label_room(router: _Router, layout: Layout, xy, labels: list) -> bool:
    """An end label at xy keeps clear of the loops routed so far (their grid cells), every element's and stem's
    body (an element by its drawn outline, na_labels._outline: a strand's shaft and head, not the pitch-wide box
    round it, as base letters are judged; held-out 2FD8's 5′ between two strands of its sheet), the text and stub
    boxes, and the end labels already placed."""
    from .na_labels import _outline

    gw, gh = glyph(layout)  # sized from the text drawn (C4)
    hw, hh = gw / 2 + 0.05, gh / 2
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
    parts = [q for e in layout.placed.values() for q in _outline(e)] + list(layout.stems)
    return not any(_near(xy, e, hh + 0.05) for e in parts)


_TRACE_NEAR = 2  # nucleotides: an element touching a trace lies under it this far either side of what it holds


def _trace_segments(layout: Layout, s: int):
    """ok(i, a, b): may segment i of strand s's trace (from nucleotide i to i + 1) be drawn from a to b? Clear of
    every element but of the elements touching its own nucleotides or a neighbour's (_TRACE_NEAR) and the sheets
    holding them. A single strand lies on the protein it wraps where it holds it (1EYG's ssDNA over the SSB tetramer;
    layout keeps the others off it), not along the whole trace. A ghost is the strand it copies, so it is exempt
    where that is."""
    na = layout.nucleic
    holds: dict[str, set[int]] = {}  # element -> the nucleotides of strand s it touches
    for p in [*layout.placed.values(), *layout.ghosts]:
        got = {x[1] for r in range(p.sse.start, p.sse.end + 1) for x in na.contacts.get(r, ()) if x[0] == s}
        if got:
            holds[p.sse.id] = got
    blocks = [set(b) for b in layout.sheets]  # barrels of one pairing network count as one (Layout.sheet_joins: 1EYG's
    for a, b in getattr(layout, "sheet_joins", []):  # OB folds, one sheet until plan 12 drew each barrel apart)
        ka, kb = (next((k for k, x in enumerate(blocks) if e in x), None) for e in (a, b))
        if ka is not None and kb is not None and ka != kb:
            blocks[ka] |= blocks[kb]
            blocks[kb] = set()
    for block in blocks:  # a sheet is placed whole: all its strands lie where its touching ones do
        got = set().union(*(holds.get(e, set()) for e in block))
        for e in sorted(block) if got else ():
            holds[e] = got
    rects, owners = _element_obstacles(layout, ghost_owner=True)  # the router's obstacles; a ghost is its strand's
    free = [k for k, o in enumerate(owners) if o not in holds]  # never exempt
    free_rects, free_owners = rects[:, free], [owners[k] for k in free]
    held = [k for k, o in enumerate(owners) if o in holds]

    def ok(i: int, a, b) -> bool:
        seg = _xy([a, b])
        if not _clear(seg, (), free_rects, free_owners, own_body=True):
            return False
        lo, hi = i - _TRACE_NEAR, i + 1 + _TRACE_NEAR
        keep = [k for k in held if not any(lo <= t <= hi for t in holds[owners[k]])]
        return not keep or _clear(seg, (), rects[:, keep], [owners[k] for k in keep], own_body=True)

    return ok


def _trace_clear(layout: Layout, s: int, pts) -> bool:
    """Is strand s's trace (drawn at its projected place) clear of every element, each segment as _trace_segments
    judges it?"""
    ok = _trace_segments(layout, s)
    xy = _xy(pts)
    return all(ok(i, xy[i], xy[i + 1]) for i in range(len(xy) - 1))


_OFF_STEPS = tuple(1.0 + k / 4 for k in range(7))  # beads (_TRACE_HALF) a trace nucleotide may be nudged off an
# element not holding it, nearest first: up to 2.5 (random-100 held-out 13FZ: B:3 lay 1.07 inside A:1027-1031's end)
_OFF_TURNS = tuple(np.radians(a) for a in range(0, 360, 5))  # in any direction, the first that clears


def off_unheld(layout: Layout, s: int, pts) -> np.ndarray:
    """Strand s's trace with each nucleotide whose two segments are not clear (_trace_segments: they cross an element
    that does not hold them) nudged, as little as clears them both (_OFF_STEPS, _OFF_TURNS), to a spot at least a
    bead's width (_TRACE_HALF) from every other trace nucleotide; one no nudge clears stays (the trace is then a
    fallback). Run once the elements are placed (layout.layout): random-100 held-out 13FZ's RNA B, held by its sheet at
    nucleotides 5-6, projected 1.07 inside the end of A:1027-1031 of that sheet at B:3, which nothing holds there."""
    from .layout import _TRACE_HALF, _snap

    ok = _trace_segments(layout, s)
    xy = np.array(pts, float).reshape(-1, 2)
    n = len(xy)
    others = [np.asarray(q, float).reshape(-1, 2) for t, q in sorted(layout.na_traces.items()) if t != s]

    def fine(i: int, q) -> bool:  # both of nucleotide i's segments clear with it at q
        return (i == 0 or ok(i - 1, xy[i - 1], q)) and (i == n - 1 or ok(i, q, xy[i + 1]))

    for i in range(n):
        if fine(i, xy[i]):
            continue
        rest = np.vstack([np.delete(xy, i, axis=0), *others])
        spot = next(
            (
                q
                for step in _OFF_STEPS
                for t in _OFF_TURNS
                for q in [_snap(xy[i] + step * _TRACE_HALF * np.array([math.cos(t), math.sin(t)]))]
                if (not len(rest) or float(np.hypot(*(rest - q).T).min()) >= _TRACE_HALF) and fine(i, q)
            ),
            None,
        )
        if spot is not None:
            xy[i] = spot
    return xy


_GAP_BOWS = (1.2, 1.6, 2.2)  # a gap with no way on the grid may bow out to this many times its chord


@dataclass(eq=False)
class _Point:
    """A point a strand's route leaves from or enters, as the owner of that port: the end of a strand drawn along its
    projection (a trace, which no stem holds) or the free cell beside a stem's end a route starts from (_side_start).
    Its rect is the point and its body is empty, so no route must keep out of it."""

    cx: float
    cy: float
    id: str = "point"
    length: float = 0.0
    width: float = 0.0
    angle: float = 0.0

    @property
    def rect(self) -> tuple[float, float, float, float]:
        return self.cx, self.cy, self.cx, self.cy


def _from_beside(router: _Router, a, b):
    """A strand's route (a gap's, a run's or a joint's) when the way straight out of an end (a = (port, outward
    direction, owner), b the same) has no free cell, even loose (_Grid.free_cell_along): from the nearest free cell
    ahead of that end (_side_start), reached by a straight stub, to the other end; None if neither end needs it or no
    route is found (random-100 C2b-6: 7Y7C's domains view, the gap after a:1439, whose end faced P80a 0.4 to 3.15 out,
    was drawn straight across three stems)."""
    grid = router.grid
    ends = []
    for port, out, own in (a, b):
        side = _side_start(router, port, out, own) if grid.free_cell_along(port, out, own, True) is None else None
        if side is None:
            ends.append((port, out, own, False))
        else:
            v = side - np.asarray(port, float)
            ends.append((side, v / float(np.hypot(*v)), _Point(float(side[0]), float(side[1])), True))
    if not (ends[0][3] or ends[1][3]):
        return None
    (qa, ua, oa, sa), (qb, ub, ob, sb) = ends
    mid, _ = router.route(qa, ua, oa, qb, ub, ob, own_body=True, loose=True, text_toll=TEXT_TOLL)
    if mid is None:
        return None
    return ([tuple(map(float, a[0]))] if sa else []) + list(mid) + ([tuple(map(float, b[0]))] if sb else [])


def _side_start(router: _Router, port, out, own):
    """The nearest free grid cell ahead of a port (within _PUSH's reach), whose centre a straight stub from the port
    reaches clear of every element and stem (its own by its body): where a gap starts when the way straight out of its
    end runs through a neighbouring stem's margin for the whole of _PUSH. None if there is none. Nearest first, the
    lower cell on a tie."""
    grid = router.grid
    port, out = np.asarray(port, float), np.asarray(out, float)
    reach = _PUSH[-1]
    c0 = grid.cell(port)
    n = int(math.ceil(reach / grid.res)) + 1
    found = []
    for dx in range(-n, n + 1):
        for dy in range(-n, n + 1):
            c = (c0[0] + dx, c0[1] + dy)
            if not grid.free(c):
                continue
            v = np.subtract(grid.center(c), port)
            d = float(np.hypot(*v))
            if d <= reach and float(v @ out) > 0:
                found.append((round(d, 6), c))
    for _, c in sorted(found):
        q = tuple(float(x) for x in grid.center(c))
        if _clear([tuple(port), q], (own,), *router.elements, own_body=True):
            return np.array(q)
    return None


def _gap(router: _Router, g, a, b, strands=None) -> NAPath:
    """A cropped gap from one strand piece's 3' end to the next one's 5' end: routed, dashed and (unless short)
    labelled with the number of nucleotides left out. An end no stem holds (a trace's: a strand drawn along its
    projection) is routed from as a point (_Point): random-100 C2b-6, 7Y7C's gaps after Z:37, Z:61 and V:21 were
    never routed and, with no curve or bow clear, drawn straight across six stems and a terminus. As a run's first
    pass does, it keeps off the hairpin loops and free ends drawn (`strands`, a _Strands: their cells blocked) when a
    route does, and only then crosses them at the reuse cost (6S13's gaps after a:36 and a:462 ran across the free
    5' end of a:37-72 they lead to). Plan 28 F4: with the free ends blocked it may also arrive along the chord from
    where it comes, and two ends nearer than two grid cells (turned to meet: 4AS1) are joined by the chord itself when
    it is clear and crosses no strand drawn so far; a longer chord is never taken."""
    (pa, ea, own_a), (pb, eb, own_b) = a, b
    own_a = own_a if own_a is not None else _Point(float(pa[0]), float(pa[1]))
    own_b = own_b if own_b is not None else _Point(float(pb[0]), float(pb[1]))
    own = (own_a, own_b)
    route = None
    if strands is not None and strands.fixed_cells is not None and strands.fixed_cells.any():
        grid = router.grid
        chord = np.subtract(pa, pb)
        chord = chord / max(float(np.hypot(*chord)), 1e-9)
        # arriving along the end's own outward heading first; then, with the free ends still blocked, along the chord
        # from where the gap comes (plan 28 F4: 4AS1's gap to B:36's tip, whose heading pointed away from B:31, found
        # no blocked route that way and fell to one crossing the tail it joins)
        for heading in (eb, chord):
            saved, used = grid.blocked.copy(), grid.used.copy()
            grid.blocked |= strands.fixed_cells
            try:
                route, _ = router.route(pa, ea, own_a, pb, heading, own_b, own_body=True, text_toll=TEXT_TOLL)
            finally:
                grid.blocked = saved
            if route is not None and strands.crosses(_xy(route), fixed_only=True):
                route, grid.used = None, used  # the cells it marked used are not drawn
            if route is not None:
                break
        straight = [tuple(map(float, pa)), tuple(map(float, pb))]
        if (  # two ends turned to meet (4AS1, F4), nearer than a grid route can join: the chord itself, when it is
            # that short (two grid cells), clear of the elements and crossing no strand drawn so far
            route is None
            and math.dist(pa, pb) <= 2 * grid.res
            and _na_clear(straight, own, router.elements)
            and not strands.crosses(straight)
        ):
            route = straight
    if route is None:
        route, _ = router.route(pa, ea, own_a, pb, eb, own_b, own_body=True, text_toll=TEXT_TOLL)
    if route is None:  # an end in another element's clearance margin
        route, _ = router.route(pa, ea, own_a, pb, eb, own_b, own_body=True, loose=True, text_toll=TEXT_TOLL)
    if route is None:  # an end whose way straight out lies in a neighbour's margin all along: from beside it
        route = _from_beside(router, (pa, ea, own_a), (pb, eb, own_b))
    if route is None:  # no way on the grid: a curve leaving and entering along the ends, or a bow either side
        chord = math.dist(pa, pb)
        shapes = [_xy(_curve(pa, ea, pb, eb))]
        shapes += [
            _bow(pa, pb, chord * f, out=s * np.array([pa[1] - pb[1], pb[0] - pa[0]]))
            for f in _GAP_BOWS
            for s in (1, -1)
        ]
        route = next((q for q in shapes if _na_clear(q, own, router.elements)), None)
    fallback = route is None
    pts = _xy(route if route is not None else [pa, pb])
    if not fallback:
        pts = next((q for q in (_rounded(pts),) if _na_clear(q, own, router.elements)), pts)
    label = f"… {g.count} nt …" if g.count >= SHORT_GAP else ""
    fallback = fallback or not _na_clear(pts, own, router.elements)
    return NAPath("gap", g.before, [], pts, [], label, fallback)
