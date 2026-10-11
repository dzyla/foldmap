# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Orthogonal loop routing: A* on a coarse grid around the placed elements."""

from __future__ import annotations

import heapq
import logging
import math
from dataclasses import dataclass

import numpy as np

from .features import KINK_GAP
from .layout import _GLYPH_H, _GLYPH_W, NT_ROW, NT_STEP, Layout, Placed, PlacedStem, _stem_label_box, label_boxes
from .model import SSE

RES = 0.25  # grid cell size, page units
_INFLATE = 0.15  # keep loops this far from arrows
# A drawn loop nearer than this to an element counts as crossing it. A free cell centre reaches RES/2 into the
# _INFLATE clearance margin (at any res, the reach makes up the rest of a coarser cell), and the port snap
# moves a run anywhere inside its cells, so a snapped route can come as close as _INFLATE - RES/2 = 0.025:
# 0.02 leaves a margin. A drawing that still breaches it is redrawn unsnapped, then becomes the flagged fallback.
_TOUCH = 0.02
_TURN = 2.5
_REUSE = 6.0  # extra cost for a cell an earlier loop already uses
_PAD = 3.0
_MAX_NODES = 200_000
_MAX_CELLS = 1_500_000
_DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
_PUSH = [0.4 + 0.25 * k for k in range(12)]  # how far out along the exit the path may start

log = logging.getLogger(__name__)


@dataclass(eq=False)
class Loop:
    a_id: str  # source element (loop leaves its C-terminal end)
    b_id: str  # target element (loop enters its N-terminal end)
    points: list[tuple[float, float]]
    dashed: bool  # residues missing between the two elements
    fallback: bool  # no orthogonal route found; points trace a smooth curve
    crosses_labels: bool = False  # only a route over a label or chain-end stub was found


def _oid(e) -> str:
    """The owner id of an element's or a stem's obstacle rect."""
    return e.sse.id if isinstance(e, Placed) else f"stem:{e.id}"


def _element_obstacles(layout: Layout) -> tuple[np.ndarray, list]:
    """The rects every drawn connector must keep clear of (elements, ghosts, stems) as a (4, n) array of x0, y0,
    x1, y1 rows, and the owner id of each (None for a ghost)."""
    own = [(_oid(p), p.rect) for p in layout.placed.values()] + [(None, g.rect) for g in layout.ghosts]
    own += [(_oid(d), d.rect) for d in layout.stems]
    return np.array([r for _, r in own], float).reshape(-1, 4).T.copy(), [k for k, _ in own]


class _Grid:
    def __init__(self, layout: Layout, keep_out: list[tuple[float, float, float, float]] = ()):
        x0, y0, x1, y1 = layout.bounds
        w, h = x1 - x0 + 2 * _PAD, y1 - y0 + 2 * _PAD
        self.res = max(RES, math.sqrt(w * h / _MAX_CELLS))
        self.x0, self.y0 = x0 - _PAD, y0 - _PAD
        self.nx, self.ny = int(math.ceil(w / self.res)), int(math.ceil(h / self.res))
        self.used = np.zeros((self.nx, self.ny), np.int32)
        self.xs = self.x0 + (np.arange(self.nx) + 0.5) * self.res
        self.ys = self.y0 + (np.arange(self.ny) + 0.5) * self.res

        # On a coarse grid a cell is blocked if any part of it meets the obstacle: the open interval of blocked
        # centres, the rect grown by _INFLATE and then by reach on each side, is w + res + 0.05 > res long on
        # each axis, so even a zero-width obstacle always blocks at least one cell.
        self.reach = (self.res - RES) / 2
        self.blocked = self._mask(p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems])
        self.soft = self._mask(keep_out) & ~self.blocked  # text, stubs: avoided, given up only if nothing else routes
        self.blocked |= self.soft

    def _mask(self, rects) -> np.ndarray:
        xs, ys, reach = self.xs, self.ys, self.reach
        m = np.zeros((self.nx, self.ny), bool)
        for rx0, ry0, rx1, ry1 in rects:
            m |= np.outer(
                (xs + reach > rx0 - _INFLATE) & (xs - reach < rx1 + _INFLATE),
                (ys + reach > ry0 - _INFLATE) & (ys - reach < ry1 + _INFLATE),
            )
        return m

    def keep_off(self, rects) -> None:
        """More text (or drawn strand) for routes to avoid like the keep-out boxes: given up only if nothing
        else routes."""
        more = self._mask(rects) & ~self.blocked
        self.soft |= more
        self.blocked |= more

    def cell(self, pt) -> tuple[int, int]:
        return int(math.floor((pt[0] - self.x0) / self.res)), int(math.floor((pt[1] - self.y0) / self.res))

    def center(self, c) -> tuple[float, float]:
        return self.x0 + (c[0] + 0.5) * self.res, self.y0 + (c[1] + 0.5) * self.res

    def free(self, c) -> bool:
        return 0 <= c[0] < self.nx and 0 <= c[1] < self.ny and not self.blocked[c]

    def cells_along(self, p, q) -> set[tuple[int, int]]:
        """Free cells under the straight segment p-q."""
        n = max(2, int(math.dist(p, q) / (self.res / 4)))
        out = set()
        for t in np.linspace(0.0, 1.0, n):
            c = self.cell((p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t))
            if self.free(c):
                out.add(c)
        return out

    def stub_cells(self, port, end) -> set[tuple[int, int]]:
        """Free cells under the straight stub from a port to its first routing cell (that cell excluded),
        so no other part of the same loop can pass between the two and be crossed by the stub."""
        return self.cells_along(port, self.center(end)) - {end}

    def free_cell_along(self, port, exit_vec, own: Placed, loose: bool = False):
        """First free cell on the way out of a port. Only the cells the element's own bounding box blocks
        (a tilted helix end can sit inside it) may be skipped; any other obstacle in the way means no route.
        `loose` (a strand's loose end, which may sit in another element's clearance margin): any blocked cell
        may be skipped, and the drawn route's final check decides."""
        g = _INFLATE + self.reach
        rx0, ry0, rx1, ry1 = own.rect
        for t in _PUSH:
            pt = (port[0] + exit_vec[0] * t, port[1] + exit_vec[1] * t)
            c = self.cell(pt)
            if self.free(c):
                return c
            cx, cy = self.center(c)
            if not (loose or (rx0 - g < cx < rx1 + g and ry0 - g < cy < ry1 + g)):
                return None
        return None


def clip_len(p, q, rect) -> float:
    """Length of the part of segment p-q strictly inside rect (Liang-Barsky clipping); 0.0 when it misses
    or only touches the edge."""
    x0, y0, x1, y1 = rect
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for pk, qk in ((-dx, p[0] - x0), (dx, x1 - p[0]), (-dy, p[1] - y0), (dy, y1 - p[1])):
        if pk == 0:
            if qk <= 0:
                return 0.0
            continue
        if pk < 0:
            t0 = max(t0, qk / pk)
        else:
            t1 = min(t1, qk / pk)
        if t0 >= t1:
            return 0.0
    return (t1 - t0) * math.hypot(dx, dy)


def _crosses_body(p, q, e: Placed) -> bool:
    """True when segment p-q runs through the element's own (possibly tilted) body, not just its end."""
    c, s = math.cos(e.angle), math.sin(e.angle)
    u0, v0 = (p[0] - e.cx) * c + (p[1] - e.cy) * s, (p[1] - e.cy) * c - (p[0] - e.cx) * s  # along, across it
    u1, v1 = (q[0] - e.cx) * c + (q[1] - e.cy) * s, (q[1] - e.cy) * c - (q[0] - e.cx) * s
    hl, hw = e.length / 2 - 1e-6, e.width / 2 - 1e-6
    return clip_len((u0, v0), (u1, v1), (-hl, -hw, hl, hw)) > 1e-9


def _clear(path, own: tuple, rects: np.ndarray, owners: list, grow=_TOUCH, kink=False, own_body=False) -> bool:
    """True when the connector `path` between its own two elements comes no nearer than `grow` to any of
    `rects` (a (4, n) array of x0, y0, x1, y1 rows; owners[k] is the element id owning rect k, or None), by
    exact clipping. The end segments leave from ports on the own elements' boxes (inside them when tilted),
    so the own boxes are not checked there. A route with turns must still leave and enter its ports from
    outside: its end segments may touch the own elements' true outlines but not enter them. A kink link
    (`kink`: the straight join of two pieces of one element) skips its own two elements and their labels.
    `own_body` (a nucleic-acid strand, drawn from where it leaves its own stems' bodies, curved or not): no
    segment may enter an own element's (possibly tilted) body, and its bounding box does not count."""
    skip = tuple(_oid(e) for e in own)
    xs, ys = [p[0] for p in path], [p[1] for p in path]
    x0, y0, x1, y1 = min(xs) - grow, min(ys) - grow, max(xs) + grow, max(ys) + grow
    rx0, ry0, rx1, ry1 = rects
    segs = list(zip(path, path[1:]))
    last = len(segs) - 1
    for k in np.flatnonzero((rx0 <= x1) & (rx1 >= x0) & (ry0 <= y1) & (ry1 >= y0)).tolist():
        box = (rx0[k] - grow, ry0[k] - grow, rx1[k] + grow, ry1[k] + grow)
        mine = owners[k] in skip
        for i, (p, q) in enumerate(segs):
            if mine and (own_body or i == 0 or i == last):
                continue
            if (
                max(p[0], q[0]) > box[0]
                and min(p[0], q[0]) < box[2]
                and max(p[1], q[1]) > box[1]
                and min(p[1], q[1]) < box[3]  # else it misses the box: no need to clip
                and clip_len(p, q, box) > 1e-9
            ):
                return False
    ends = range(len(segs)) if own_body else sorted({0, last})
    return kink or not any(_crosses_body(*segs[i], e) for i in ends for e in own)


def _axis(vec) -> int:
    x, y = vec
    if abs(x) >= abs(y):
        return 0 if x > 0 else 2
    return 1 if y > 0 else 3


def _astar(grid: _Grid, start, start_dir, goal, goal_dir):
    """Cheapest orthogonal path. Either end may be taken in any heading but one pointing back into its
    element, paying a turn when it is not the port's own axis; the goal pushes a terminal node. (Insisting
    on the port's heading made paths circle the end cell to take it: a visible curl.)"""
    tie = 0
    h0 = abs(goal[0] - start[0]) + abs(goal[1] - start[1])
    best, heap = {}, []
    for d in range(4):  # leave in any heading but back into the element; a sideways start pays its turn
        if d != (start_dir + 2) % 4:
            g0 = 0.0 if d == start_dir else _TURN
            best[(start[0], start[1], d)] = g0
            heap.append((g0 + h0, d, g0, start[0], start[1], d, False))
            tie = max(tie, d)
    heapq.heapify(heap)
    parent: dict = {}
    nodes = 0
    while heap:
        _, _, g, x, y, d, done = heapq.heappop(heap)
        if done:
            cells = [(x, y)]
            key = (x, y, d)
            while key in parent:
                key = parent[key]
                cells.append((key[0], key[1]))
            return cells[::-1]
        if g > best.get((x, y, d), math.inf):
            continue
        if (x, y) == goal and d != (goal_dir + 2) % 4:
            tie += 1
            heapq.heappush(heap, (g + (_TURN if d != goal_dir else 0), tie, g, x, y, d, True))
            continue
        nodes += 1
        if nodes > _MAX_NODES:
            return None
        for nd in range(4):
            if nd == (d + 2) % 4:
                continue
            nx_, ny_ = x + _DIRS[nd][0], y + _DIRS[nd][1]
            if not grid.free((nx_, ny_)):
                continue
            cost = g + 1 + (_TURN if nd != d else 0) + (_REUSE if grid.used[nx_, ny_] else 0)
            key = (nx_, ny_, nd)
            if cost < best.get(key, math.inf):
                best[key] = cost
                parent[key] = (x, y, d)
                tie += 1
                h = abs(goal[0] - nx_) + abs(goal[1] - ny_)
                heapq.heappush(heap, (cost + h, tie, cost, nx_, ny_, nd, False))
    return None


def _corners(cells):
    out = [cells[0]]
    for a, b, c in zip(cells, cells[1:], cells[2:]):
        if (b[0] - a[0], b[1] - a[1]) != (c[0] - b[0], c[1] - b[1]):
            out.append(b)
    out.append(cells[-1])
    return out


def _snap_to_port(pts: list[list[float]], port, exit_vec, at_start: bool) -> None:
    """Move the first (or last) straight run onto the port's axis line so the stub is exactly axial."""
    if min(abs(exit_vec[0]), abs(exit_vec[1])) > 1e-6:
        return  # tilted helix end: leave the short stub diagonal
    perp = 0 if abs(exit_vec[0]) < 1e-6 else 1
    seq = range(len(pts)) if at_start else range(len(pts) - 1, -1, -1)
    ref = pts[seq[0]][perp]
    moved = []
    for i in seq:
        if abs(pts[i][perp] - ref) > 1e-9:
            break
        moved.append(i)
    if len(moved) == len(pts):  # one straight run joins both ports; only the start snaps
        if not at_start:
            return
    for i in moved:
        pts[i][perp] = port[perp]


def _simplify(points):
    out = [points[0]]
    for a, b, c in zip(points, points[1:], points[2:]):
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if abs(cross) > 1e-9:
            out.append(b)
    out.append(points[-1])
    return out


def _curve(pa, ea, pb, eb, n=24):
    p0, p3 = np.array(pa), np.array(pb)
    p1, p2 = p0 + np.array(ea) * 2.0, p3 + np.array(eb) * 2.0
    t = np.linspace(0, 1, n)[:, None]
    pts = (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3
    pts[0], pts[-1] = p0, p3
    return [(float(x), float(y)) for x, y in pts]


def _route_pair(grid: _Grid, a: Placed, b: Placed):
    """Grid cells of the cheapest orthogonal route from a's C port to b's N port, or None."""
    return _route_ports(grid, a.c_port, a.exit_c, a, b.n_port, b.exit_n, b)


def _route_ports(grid: _Grid, pa, ea, own_a, pb, eb, own_b, loose=False):
    """Grid cells of the cheapest orthogonal route from port pa (leaving along ea) to port pb (entered against
    eb), or None; own_a and own_b (anything with a rect) are what the ports belong to. `loose`: see
    _Grid.free_cell_along."""
    start = grid.free_cell_along(pa, ea, own_a, loose)
    goal = grid.free_cell_along(pb, eb, own_b, loose)
    if start is None or goal is None:
        return None
    stubs = grid.stub_cells(pa, start) | grid.stub_cells(pb, goal)  # the straight ends drawn to each port
    for c in stubs:
        grid.blocked[c] = True
    cells = _astar(grid, start, _axis(ea), goal, _axis(-np.asarray(eb)))
    for c in stubs:
        grid.blocked[c] = False
    return cells


def _polyline(grid: _Grid, cells, pa, ea, pb, eb, snap: bool):
    """The route through `cells` from port to port, its end runs moved onto the ports' axes if `snap`."""
    pts = [list(grid.center(c)) for c in _corners(cells)]
    if snap:
        _snap_to_port(pts, pa, ea, True)
        _snap_to_port(pts, pb, eb, False)
    return _simplify([tuple(pa)] + [(x, y) for x, y in pts] + [tuple(pb)])


class _Router:
    """Draws routes on one grid and validates each drawing by _clear against the elements (and stems)."""

    def __init__(self, grid: _Grid, elements: tuple[np.ndarray, list]):
        self.grid, self.elements = grid, elements
        self.tight = self.rejected = 0

    def drawn(self, pa, ea, own_a, pb, eb, own_b, own_body=False, loose=False):
        """The A* route from port pa to port pb as drawn, snapped to the ports or else unsnapped; None if there
        is no route or every drawing of it crosses an element."""
        grid = self.grid
        cells = _route_ports(grid, pa, ea, own_a, pb, eb, own_b, loose)
        if cells is None:
            return None
        for snap in (True, False):
            pts = _polyline(grid, cells, pa, ea, pb, eb, snap)
            if _clear(pts, (own_a, own_b), *self.elements, own_body=own_body):
                for c in cells:
                    grid.used[c] += 1
                if log.isEnabledFor(logging.DEBUG):
                    self.tight += not _clear(pts, (own_a, own_b), *self.elements, grow=_INFLATE - 1e-6)
                return pts
            self.rejected += 1
        return None

    def route(self, pa, ea, own_a, pb, eb, own_b, own_body=False, loose=False):
        """(points, over labels): a drawn route, kept off text and stubs unless nothing else routes; None
        points if there is no clear route at all."""
        grid = self.grid
        pts = self.drawn(pa, ea, own_a, pb, eb, own_b, own_body, loose)
        if pts is not None or not grid.soft.any():
            return pts, False
        grid.blocked &= ~grid.soft
        pts = self.drawn(pa, ea, own_a, pb, eb, own_b, own_body, loose)
        grid.blocked |= grid.soft
        return pts, pts is not None


def route_loops(layout: Layout, sses: list[SSE], bb) -> list[Loop]:
    unjoined = getattr(layout, "unjoined", set())  # neighbours' fragments with unshown elements between
    pairs = [
        (a, b)
        for a, b in zip(sses, sses[1:])
        if a.chain == b.chain
        and a.id in layout.placed
        and b.id in layout.placed
        and (a.id, b.id) not in unjoined
        and a.chain not in getattr(layout, "partial", ())  # neighbours' fragments: shown, not traced
    ]
    if not pairs and layout.nucleic is None:
        return []
    keep_out = label_boxes(layout, sses)
    grid = _Grid(layout, [b for _, b in keep_out])
    # the elements every drawn loop must not cross, then the text and stubs a kink link also keeps clear of
    elements = _element_obstacles(layout)
    owners = elements[1] + [n.removeprefix("label:") if n.startswith("label:") else None for n, _ in keep_out]
    rects = np.hstack([elements[0], np.array([r for _, r in keep_out], float).reshape(-1, 4).T])
    router = _Router(grid, elements)

    dist = lambda ab: math.dist(layout.placed[ab[0].id].c_port, layout.placed[ab[1].id].n_port)  # noqa: E731
    routed: dict[tuple[str, str], tuple[list, bool, bool]] = {}  # points, fallback curve, over labels
    for a, b in sorted(pairs, key=dist):
        pa, pb = layout.placed[a.id], layout.placed[b.id]
        gap = np.subtract(pb.n_port, pa.c_port)
        # a kink is drawn straight with no grid clearance, so it clears by the full _INFLATE (drawn routes use _TOUCH)
        if (
            b.start - a.end - 1 <= KINK_GAP
            and np.hypot(*gap) < 1.5
            and gap @ pa.exit_c > 0
            and _clear([pa.c_port, pb.n_port], (pa, pb), rects, owners, grow=_INFLATE, kink=True)  # and labels
        ):
            for c in grid.cells_along(pa.c_port, pb.n_port):
                grid.used[c] += 1
            routed[(a.id, b.id)] = ([pa.c_port, pb.n_port], False, False)  # a kink: the helix just carries on
            continue
        pts, over = router.route(pa.c_port, pa.exit_c, pa, pb.n_port, pb.exit_n, pb)
        if pts is None:
            routed[(a.id, b.id)] = (_curve(pa.c_port, pa.exit_c, pb.n_port, pb.exit_n), True, False)
        else:
            routed[(a.id, b.id)] = (pts, False, over)
    log.debug(
        "%d loops: %d fallback curves, %d drawn nearer than %.2f to an element, %d drawings rejected as crossing",
        len(pairs),
        sum(fb for _, fb, _ in routed.values()),
        router.tight,
        _INFLATE,
        router.rejected,
    )
    out = []
    for a, b in pairs:
        pts, fb, over = routed[(a.id, b.id)]
        dashed = any(not bb.prev[k] for k in range(a.end + 1, b.start + 1))
        out.append(Loop(a.id, b.id, [(float(x), float(y)) for x, y in pts], dashed, fb, over))
    if layout.nucleic is not None:  # after the loops, so strands see the cells they use
        layout.na_paths = route_nucleic(layout, router)
    return out


_BOW_MAX = 8  # internal-loop sides up to this many nucleotides are drawn as bows, not routed
SHORT_GAP = 3  # a cropped gap of fewer nucleotides is drawn as a bare break, without its "… N nt …" label
_JOINT_STRAIGHT = 1.5  # a joint between two pieces of a bent duplex shorter than this may be drawn straight


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


def _dot(a, b) -> float:
    """a . b of two 2-vectors, written out: numpy's dot goes through BLAS, whose kernels differ in the last bit."""
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1])


def _xy(pts) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in pts]


def _length(pts) -> float:
    p = np.asarray(pts, float)
    return float(np.hypot(*np.diff(p, axis=0).T).sum()) if len(p) > 1 else 0.0


def _spaced(pts, m: int) -> list[tuple[float, float]]:
    """m points evenly spaced along the polyline, ends excluded."""
    p = np.asarray(pts, float)
    if m <= 0 or len(p) < 2:
        return []
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    out = []
    for t in np.arange(1, m + 1) / (m + 1) * cum[-1]:
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
    half = float(np.hypot(*(pb - pa))) / 2
    r = max(half * 1.02, (m + 1) * NT_STEP / (1.6 * np.pi))
    centre = (pa + pb) / 2 + out * np.sqrt(max(r * r - half * half, 0.0))
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


def _tail(d: PlacedStem, s: int, key, nts: list[int], five: bool, elements) -> tuple[NAPath, tuple, np.ndarray]:
    """A free end drawn straight on, or bent after its first nucleotide when that is what keeps it (and its end
    label) clear of the elements: (path, end label position, outward direction at its tip)."""
    for turn in _TURNS:
        pts, xy, tip = _free_end(d, key, len(nts), five, turn)
        drawn = pts[::-1] if five else pts
        if _na_clear([*drawn, tip], (d,), elements):
            break
    else:
        pts, xy, tip = _free_end(d, key, len(nts), five)
        turn = None
    path = NAPath("tail", s, nts, pts, xy, fallback=turn is None, between=(key,))
    far = xy[0] if five else xy[-1]
    v = np.subtract(tip, far)
    return path, tip, v / max(float(np.hypot(*v)), 1e-9)


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


def _in_body(pt, d) -> bool:
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


def _outside(pts, own) -> list[tuple[float, float]]:
    """The polyline without its first and last stretches inside its own stems' bodies (a strand runs along its
    rail to the stem's end face before it is drawn among the other elements)."""

    def trim(p):
        for i in range(len(p) - 1):
            mine = [d for d in own if _in_body(p[i], d)]
            if not mine:
                return p[i:]
            if any(_in_body(p[i + 1], d) for d in own):
                continue
            t = max(_leaves_at(p[i], p[i + 1], d) for d in mine)
            q = np.add(p[i], t * np.subtract(p[i + 1], p[i]))
            return [(float(q[0]), float(q[1]))] + list(p[i + 1 :])
        return [p[-1], p[-1]]

    out = trim(trim(list(pts))[::-1])[::-1]
    return out if len(out) > 1 else [out[0], out[0]]


def _na_clear(pts, own, elements) -> bool:
    """A strand drawn clear of every element and stem (by _clear), from where it leaves its own stems."""
    return _clear(_outside(pts, own), own, *elements, own_body=True)


def route_nucleic(layout: Layout, router: _Router | None = None) -> list[NAPath]:
    """Every single-stranded stretch as drawn: runs routed around the elements like loops (as are the joints
    between pieces of a bent duplex), hairpin loops as arcs, free ends straight on, strands without pairs along
    their projection, cropped gaps dashed. Every drawing is checked by _clear; one that fails is a fallback.
    Also fills layout.na_ends (where each strand's 5' and 3' labels go) and layout.na_letters (where each
    single-stranded nucleotide's letter goes), and widens layout.bounds to hold it all."""
    na = layout.nucleic
    if router is None:
        router = _Router(_Grid(layout, [b for _, b in label_boxes(layout, [])]), _element_obstacles(layout))
    elements = router.elements
    at = {key: d for d in layout.stems for key in d.axial}
    kind_of = {(r.strand, r.lo): r.kind for r in na.structure.runs}
    ends: list[tuple[int, str, tuple[float, float]]] = []
    exits: dict[tuple[int, int], tuple] = {}  # strand end -> (port, outward direction, owner) a gap leaves from

    def checked(path: NAPath, own) -> NAPath:
        path.fallback = path.fallback or not _na_clear(path.points, own, elements)
        return path

    gapped = {(g.before, "3′") for g in na.gaps} | {(g.after, "5′") for g in na.gaps}
    fixed: dict[tuple[int, int], NAPath] = {}  # (strand, nucleotide before it) -> a path of fixed shape
    between: list[tuple] = []  # strands from stem to stem still to route, once the fixed shapes are known
    for s, strand in enumerate(na.strands):
        n = len(strand)
        if s in layout.na_traces:
            pts = _xy(layout.na_traces[s])
            fixed[(s, -1)] = checked(NAPath("trace", s, list(range(n)), pts, pts), ())
            for label, a, b in (("5′", pts[0], pts[1]), ("3′", pts[-1], pts[-2])):
                v = np.subtract(a, b)
                v = v / max(float(np.hypot(*v)), 1e-9)
                ends.append((s, label, tuple(float(c) for c in np.add(a, v * (NT_ROW + _GLYPH_W)))))
                exits[(s, 0 if label == "5′" else n - 1)] = (np.add(a, v * NT_STEP / 2), v, None)
            continue
        held = [i for i in range(n) if (s, i) in at]
        if not held:
            continue
        first, last = held[0], held[-1]
        d0 = at[(s, first)]
        if first > 0:
            path, tip, v = _tail(d0, s, (s, first), list(range(first)), True, elements)
            fixed[(s, -1)] = path
            exits[(s, 0)] = (np.add(path.nt_xy[0], v * NT_STEP / 2), v, d0)
        else:
            port, along = d0.port((s, 0), False)
            tip = tuple(float(c) for c in port + along * NT_ROW)
            exits[(s, 0)] = (*_face(d0, (s, 0), False), d0)
        ends.append((s, "5′", tip))
        for i, j in zip(held, held[1:]):
            da, db = at[(s, i)], at[(s, j)]
            if da is db and da.side[(s, i)] == da.side[(s, j)]:
                continue  # along one rail of a stem (a bulge is drawn on the rail)
            if da is db:  # a hairpin loop: an arc over the tip
                fixed[(s, i)] = _strand_between(router, s, i, j, da, db, None)
            else:
                between.append((s, i, j, da, db, kind_of.get((s, i + 1))))
        dn = at[(s, last)]
        if last < n - 1:
            path, tip, v = _tail(dn, s, (s, last), list(range(last + 1, n)), False, elements)
            fixed[(s, last)] = path
            exits[(s, n - 1)] = (np.add(path.nt_xy[-1], v * NT_STEP / 2), v, dn)
        else:
            port, along = dn.port((s, last), True)
            tip = tuple(float(c) for c in port + along * NT_ROW)
            exits[(s, n - 1)] = (*_face(dn, (s, last), True), dn)
        ends.append((s, "3′", tip))
    # routes keep off the strands' end labels (as a protein loop keeps off its chain's N and C)
    g = _GLYPH_H / 2 + 0.05
    router.grid.keep_off([(x - g, y - g, x + g, y + g) for s, e, (x, y) in ends if (s, e) not in gapped])
    routed = {(s, i): _strand_between(router, s, i, j, da, db, kind) for s, i, j, da, db, kind in between}
    out = [p for _, p in sorted({**fixed, **routed}.items())]
    for g in na.gaps:  # cropped away: the 3' end of one piece to the 5' end of the next, dashed
        a, b = exits.get((g.before, len(na.strands[g.before]) - 1)), exits.get((g.after, 0))
        if a is not None and b is not None:
            out.append(_gap(router, g, a, b))
    layout.na_ends = [e for e in ends if (e[0], e[1]) not in gapped]
    layout.na_letters = _place_letters(layout, out)
    _widen(layout, out)
    if log.isEnabledFor(logging.DEBUG):
        log.debug("%d nucleic-acid paths: %d fallbacks", len(out), sum(p.fallback for p in out))
    return out


def _strand_between(router: _Router, s: int, i: int, j: int, da: PlacedStem, db: PlacedStem, kind) -> NAPath:
    """The single strand from stem nucleotide i to stem nucleotide j (the m nucleotides between drawn on it)."""
    m, nts, between = j - i - 1, list(range(i + 1, j)), ((s, i), (s, j))
    own, need = (da, db), (j - i) * NT_STEP
    ra, rb = da.point((s, i)), db.point((s, j))
    ok = lambda pts: _na_clear(pts, own, router.elements)  # noqa: E731
    if da is db:  # over the stem's tip: a hairpin loop
        pts = _arc(da, (s, i), (s, j), m)
        return NAPath("hairpin", s, nts, pts, _spaced(pts, m), fallback=not ok(pts), between=between)
    (fa, ea), (fb, eb) = _face(da, (s, i), True), _face(db, (s, j), False)
    joint = m == 0 and da.parent and da.parent == db.parent
    if joint and math.dist(fa, fb) < _JOINT_STRAIGHT and _dot(fb - fa, ea) > 0:
        pts = _xy([ra, fa, fb, rb])  # the duplex just carries on
        if ok(pts):
            for c in router.grid.cells_along(fa, fb):
                router.grid.used[c] += 1
            return NAPath("joint", s, [], pts, [], between=between)
    if kind == "internal" and m <= _BOW_MAX:
        # one side of an internal loop: bowed out on its own side, so the two sides never meet
        pts = _bow(ra, rb, need, out=da.normal * da.side[(s, i)])
        if ok(pts):
            return NAPath("run", s, nts, pts, _spaced(pts, m), between=between)
    route, _ = router.route(fa, ea, da, fb, eb, db, own_body=True)
    fallback = route is None
    if route is None:
        route = _curve(fa, ea, fb, eb)
    pts = _xy([ra, *route, rb])
    if not fallback:
        pts = next((q for q in (_rounded(pts),) if ok(q)), pts)
    if not fallback and _length(pts) < need:  # too short for its nucleotides: bowed out, if that is clear
        left = np.array([ra[1] - rb[1], rb[0] - ra[0]])
        pts = next((b for b in (_bow(ra, rb, need, out=left), _bow(ra, rb, need, out=-left)) if ok(b)), pts)
    return NAPath("joint" if joint else "run", s, nts, pts, _spaced(pts, m), fallback=fallback, between=between)


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
    taken = [d.point(key, NT_ROW) for d in layout.stems for key in d.axial] + list(where.values())
    taken += [np.asarray(xy) for _, _, xy in layout.na_ends]
    bodies = [*layout.placed.values(), *layout.stems]  # by their true (tilted) outlines, not their boxes
    g = 1.2 * _GLYPH_H  # a stem's name stands a letter and a half clear of any base letter
    boxes = [np.add(_stem_label_box(d), (-g, -g, g, g)) for d in layout.stems]
    out: dict[tuple[int, int], tuple[float, float]] = {}

    def room(q, own) -> float:
        if any(x0 < q[0] < x1 and y0 < q[1] < y1 for x0, y0, x1, y1 in boxes) or any(_near(q, e) for e in bodies):
            return 0.0
        return min((float(np.hypot(*(q - t))) for t in taken if t is not own), default=np.inf)

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
            spaces = [room(q, where[(p.strand, i)]) for q in options]
            clear = [q for q, r in zip(options, spaces) if r >= _LETTER_GAP]
            spot = clear[0] if clear else options[int(np.argmax(spaces))]
            out[(p.strand, i)] = (float(spot[0]), float(spot[1]))
            taken.append(spot)
    return out


def _widen(layout: Layout, paths: list[NAPath]) -> None:
    """Grow layout.bounds to hold the strands, every base letter, the end labels and the gaps' labels."""
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
