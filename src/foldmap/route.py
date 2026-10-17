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
from .layout import (
    _GLYPH_W,
    Layout,
    Placed,
    _element_label_box,
    _ink_em,
    _snap,
    _text_width,
    _trace_hits,
    glyph,
    helix_label_pos,
    label_boxes,
    overlaps_any,
)
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
_SLAB_COST = 3.0  # extra cost for a cell inside the membrane slab: loops between TM helices go round it
_SLAB_RUN = 12.0  # ... a step along the slab, not across it: a loop leaves the slab at a face (random-100 C26)
_PAD = 3.0
_MAX_NODES = 200_000
# a search that ran out of _MAX_NODES with its two ends joined by free cells is run once more with this budget
# (random-100 C2b-2): on a large page a route round a wall of long elements, or between ends far apart, takes more.
# Measured on the batch's pages: 7JU4's, 7QUP's and 8QEM's stuck loops (and their loops over labels) need at most
# 800 000; 8VWQ's and 7Y7C's longest nucleic-acid runs 1 600 000, at about 5 s a search
_MAX_NODES_WIDE = 1_600_000
# ... and every second search of one figure's routing (loops, exits and passes, nucleic-acid strands) together
# spends at most this many nodes, counted as A* expands them; once spent, a search that runs out of _MAX_NODES
# falls back as before (fix round 1: 7QUP's focus view took 158 s of its 180). About 25 s at most on a large page
_MAX_NODES_PAGE = 8_000_000
_MAX_CELLS = 1_500_000
_DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
_PUSH = [0.4 + 0.25 * k for k in range(12)]  # how far out along the exit the path may start
_PUSH_NEAR = [0.15 + 0.25 * k for k in range(13)]  # ... out of a port facing a stack gap (Grid.near_ports): nearer
_PUSH_OVER = [0.4 + 0.25 * k for k in range(60)]  # from a helix drawn over a sheet: out past the sheet's edge

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


def _element_obstacles(layout: Layout, ghost_owner: bool = False) -> tuple[np.ndarray, list]:
    """The rects every drawn connector must keep clear of (elements, ghosts, stems) as a (4, n) array of x0, y0,
    x1, y1 rows, and the owner id of each: None for a ghost, or with ghost_owner the id of the element it copies. A
    helix drawn over a sheet (Placed.over) is one too: only its own loops and its sheet's strands' may pass under it
    (_clear's `exempt`, from _over_exempt)."""
    ghost = lambda g: g.sse.id if ghost_owner else None  # noqa: E731
    own = [(_oid(p), p.rect) for p in layout.placed.values()]
    own += [(ghost(g), g.rect) for g in layout.ghosts]
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
        self.blocked = self._mask(
            p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems] if not getattr(p, "over", False)
        )
        # a helix over a sheet: the cells it covers beyond the strands are soft, opened to its own loops and its
        # sheet's (route_loops) and to any loop only when nothing else routes
        self.over = {k: self._mask([p.rect]) & ~self.blocked for k, p in sorted(layout.placed.items()) if p.over}
        self.soft = self._mask(keep_out) & ~self.blocked  # text, stubs: avoided, given up only if nothing else routes
        for cells in self.over.values():
            self.soft |= cells
        self.blocked |= self.soft
        self.toll: np.ndarray | None = None  # extra cost of entering a cell (text placed after routing), when set
        # ports a route may start nearer (_PUSH_NEAR), when set; so may any port whose first step is in the slab
        self.near_ports: set[tuple[float, float]] = set()
        self.slab = np.zeros((self.nx, self.ny), bool)  # cells strictly inside the membrane band (one per row)
        if layout.membrane:
            for lo, hi in layout.membrane.get("bands", [layout.membrane["y"]]):
                self.slab[:, (self.ys > lo + 0.1) & (self.ys < hi - 0.1)] = True

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

    def _into_slab(self, port, exit_vec) -> bool:
        """Does a route's usual first step out of this port (_PUSH[0] along its exit) land in the membrane slab? Such
        a route starts nearer (_PUSH_NEAR), so it turns round its own element's end at the port's own row, not 0.4
        into the slab: 6JO6's helix ending 0.03 below a face, 5LDW's in-slab helix facing a 0.81 gap under the one
        stacked on it (random-100 C2b-5)."""
        c = self.cell((port[0] + exit_vec[0] * _PUSH[0], port[1] + exit_vec[1] * _PUSH[0]))
        return bool(0 <= c[0] < self.nx and 0 <= c[1] < self.ny and self.slab[c])

    def free_cell_along(self, port, exit_vec, own: Placed, loose: bool = False):
        """First free cell on the way out of a port. Only the cells the element's own bounding box blocks
        (a tilted helix end can sit inside it) may be skipped; any other obstacle in the way means no route.
        `loose` (a strand's loose end, which may sit in another element's clearance margin; any end of a helix
        drawn over a sheet, which lies on the strands): any blocked cell may be skipped, and the drawn route's
        final check decides. From a helix over a sheet the way out runs over the strands to the sheet's edge
        (_PUSH_OVER) and stops there: a soft cell (text, another helix) next means no way out this side."""
        g = _INFLATE + self.reach
        rx0, ry0, rx1, ry1 = own.rect
        over = loose and getattr(own, "over", False)
        near = (round(float(port[0]), 6), round(float(port[1]), 6)) in self.near_ports
        near = near or self._into_slab(port, exit_vec)
        for t in _PUSH_OVER if over else _PUSH_NEAR if near else _PUSH:
            pt = (port[0] + exit_vec[0] * t, port[1] + exit_vec[1] * t)
            c = self.cell(pt)
            if self.free(c):
                return c
            if over and 0 <= c[0] < self.nx and 0 <= c[1] < self.ny and self.soft[c]:
                return None  # past the strands, text or another helix: the walk stops at the sheet's edge
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


def clip_lengths(segments, boxes) -> np.ndarray:
    """clip_len for every (segment, box) pair at once: segments (x0, y0, x1, y1) by boxes (x0, y0, x1, y1), the length
    of each segment strictly inside each box (Liang-Barsky, as clip_len). The one segment-against-box test of the
    text and tether checks (route._text_hits, na_labels letter spots, render._over_bodies; fix round 1)."""
    seg = np.asarray(segments, float).reshape(-1, 4)
    box = np.asarray(boxes, float).reshape(-1, 4)
    px, py = seg[:, 0:1], seg[:, 1:2]
    dx, dy = seg[:, 2:3] - px, seg[:, 3:4] - py
    t0 = np.zeros((len(seg), len(box)))
    t1 = np.ones((len(seg), len(box)))
    empty = np.zeros((len(seg), len(box)), bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        for pk, qk in ((-dx, px - box[:, 0]), (dx, box[:, 2] - px), (-dy, py - box[:, 1]), (dy, box[:, 3] - py)):
            pk = np.broadcast_to(pk, qk.shape)
            flat = pk == 0
            empty |= flat & (qk <= 0)
            r = qk / pk
            t0 = np.where(~flat & (pk < 0), np.maximum(t0, r), t0)
            t1 = np.where(~flat & (pk > 0), np.minimum(t1, r), t1)
    inside = ~empty & (t0 < t1)
    return np.where(inside, (t1 - t0) * np.hypot(dx, dy), 0.0)


def boxes_meeting(boxes, box) -> np.ndarray:
    """Which of `boxes` share area with `box` (open overlap), as a bool array."""
    b = np.asarray(boxes, float).reshape(-1, 4)
    return (b[:, 0] < box[2]) & (b[:, 2] > box[0]) & (b[:, 1] < box[3]) & (b[:, 3] > box[1])


def _crosses_body(p, q, e: Placed) -> bool:
    """True when segment p-q runs through the element's own (possibly tilted) body, not just its end."""
    c, s = math.cos(e.angle), math.sin(e.angle)
    u0, v0 = (p[0] - e.cx) * c + (p[1] - e.cy) * s, (p[1] - e.cy) * c - (p[0] - e.cx) * s  # along, across it
    u1, v1 = (q[0] - e.cx) * c + (q[1] - e.cy) * s, (q[1] - e.cy) * c - (q[0] - e.cx) * s
    hl, hw = e.length / 2 - 1e-6, e.width / 2 - 1e-6
    return clip_len((u0, v0), (u1, v1), (-hl, -hw, hl, hw)) > 1e-9


def _clear(
    path,
    own: tuple,
    rects: np.ndarray,
    owners: list,
    grow=_TOUCH,
    kink=False,
    own_body=False,
    under=((), ()),
    exempt=(),
) -> bool:
    """True when the connector `path` between its own two elements comes no nearer than `grow` to any of
    `rects` (a (4, n) array of x0, y0, x1, y1 rows; owners[k] is the element id owning rect k, or None), by
    exact clipping. The end segments leave from ports on the own elements' boxes (inside them when tilted),
    so the own boxes are not checked there. A route with turns must still leave and enter its ports from
    outside: its end segments may touch the own elements' true outlines but not enter them. A kink link
    (`kink`: the straight join of two pieces of one element) skips its own two elements and their labels.
    `own_body` (a nucleic-acid strand, drawn from where it leaves its own stems' bodies, curved or not): no
    segment may enter an own element's (possibly tilted) body, and its bounding box does not count. `under`: the
    owners the first and the last segment may cross (from a helix drawn over a sheet: its sheet's strands, which its
    ports lie on; into the other end: that helix, whose body may cover the end's port). `exempt`: owners the whole
    path may cross (the helices over a sheet the path's elements belong to: _over_exempt)."""
    skip = tuple(_oid(e) for e in own)
    xs, ys = [p[0] for p in path], [p[1] for p in path]
    x0, y0, x1, y1 = min(xs) - grow, min(ys) - grow, max(xs) + grow, max(ys) + grow
    rx0, ry0, rx1, ry1 = rects
    segs = list(zip(path, path[1:]))
    last = len(segs) - 1
    for k in np.flatnonzero((rx0 <= x1) & (rx1 >= x0) & (ry0 <= y1) & (ry1 >= y0)).tolist():
        box = (rx0[k] - grow, ry0[k] - grow, rx1[k] + grow, ry1[k] + grow)
        if owners[k] in exempt:
            continue
        mine = owners[k] in skip
        for i, (p, q) in enumerate(segs):
            if mine and (own_body or i == 0 or i == last):
                continue
            if (i == 0 and owners[k] in under[0]) or (i == last and owners[k] in under[1]):
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
    passes = {0: set(under[0])}  # an end over a sheet passes over its strands, its own among them
    passes[last] = passes.get(last, set()) | set(under[1])
    return kink or not any(_crosses_body(*segs[i], e) for i in ends for e in own if _oid(e) not in passes.get(i, ()))


def _axis(vec) -> int:
    x, y = vec
    if abs(x) >= abs(y):
        return 0 if x > 0 else 2
    return 1 if y > 0 else 3


def _astar(grid: _Grid, start, start_dir, goal, goal_dir, max_nodes: int = _MAX_NODES, spent: list | None = None):
    """Cheapest orthogonal path. Either end may be taken in any heading but one pointing back into its
    element, paying a turn when it is not the port's own axis; the goal pushes a terminal node. (Insisting
    on the port's heading made paths circle the end cell to take it: a visible curl.) `spent` (a one-item list):
    the nodes expanded are added to it."""
    spent = spent if spent is not None else [0]
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
    toll = grid.toll
    while heap:
        _, _, g, x, y, d, done = heapq.heappop(heap)
        if done:
            cells = [(x, y)]
            key = (x, y, d)
            while key in parent:
                key = parent[key]
                cells.append((key[0], key[1]))
            spent[0] += nodes
            return cells[::-1]
        if g > best.get((x, y, d), math.inf):
            continue
        if (x, y) == goal and d != (goal_dir + 2) % 4:
            tie += 1
            heapq.heappush(heap, (g + (_TURN if d != goal_dir else 0), tie, g, x, y, d, True))
            continue
        nodes += 1
        if nodes > max_nodes:  # this node is not expanded: the budget is spent exactly
            spent[0] += max_nodes
            return None
        for nd in range(4):
            if nd == (d + 2) % 4:
                continue
            nx_, ny_ = x + _DIRS[nd][0], y + _DIRS[nd][1]
            if not grid.free((nx_, ny_)):
                continue
            cost = g + 1 + (_TURN if nd != d else 0) + (_REUSE if grid.used[nx_, ny_] else 0)
            cost += (_SLAB_RUN if nd in (0, 2) else _SLAB_COST) if grid.slab[nx_, ny_] else 0
            if toll is not None:
                cost += toll[nx_, ny_]
            key = (nx_, ny_, nd)
            if cost < best.get(key, math.inf):
                best[key] = cost
                parent[key] = (x, y, d)
                tie += 1
                h = abs(goal[0] - nx_) + abs(goal[1] - ny_)
                heapq.heappush(heap, (cost + h, tie, cost, nx_, ny_, nd, False))
    spent[0] += nodes
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


def _route_ports(
    grid: _Grid, pa, ea, own_a, pb, eb, own_b, loose=False, max_nodes: int = _MAX_NODES, wide: "_Router | None" = None
):
    """Grid cells of the cheapest orthogonal route from port pa (leaving along ea) to port pb (entered against
    eb), or None; own_a and own_b (anything with a rect) are what the ports belong to. `loose`: see
    _Grid.free_cell_along; a pair (for pa, for pb) or one value for both. A search that finds no route within
    max_nodes is run again (_Router.second) when `wide` (the router) allows: free cells join its two ends, and the
    page's budget is not spent."""
    la, lb = loose if isinstance(loose, tuple) else (loose, loose)
    start = grid.free_cell_along(pa, ea, own_a, la)
    goal = grid.free_cell_along(pb, eb, own_b, lb)
    if start is None or goal is None:
        return None
    stubs = grid.stub_cells(pa, start) | grid.stub_cells(pb, goal)  # the straight ends drawn to each port
    for c in stubs:
        grid.blocked[c] = True
    cells = _astar(grid, start, _axis(ea), goal, _axis(-np.asarray(eb)), max_nodes)
    if cells is None and wide is not None:
        cells = wide.second(start, _axis(ea), goal, _axis(-np.asarray(eb)), max_nodes, stubs)
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

    def __init__(self, grid: _Grid, elements: tuple[np.ndarray, list], keep_out=()):
        self.grid, self.elements = grid, elements
        self.keep_out = list(keep_out)  # the text and stub boxes the grid's soft cells were made from
        self.tight = self.rejected = 0  # debug counts: routes drawn tight, drawings rejected as crossing
        self.max_nodes = _MAX_NODES  # the A* search's budget
        self.wide_nodes = _MAX_NODES_WIDE  # ... and its second search's, when a route exists (_route_ports)
        self.wide_left = _MAX_NODES_PAGE  # what the second searches of this page may still spend
        self.wide_searches = self.wide_refused = 0  # debug counts: second searches run, refused for the budget
        self._free: tuple[bytes, np.ndarray] | None = None  # the last free-cell labelling (_joined) and its key

    def _joined(self, a, b, stubs) -> bool:
        """Are cells a and b joined by free cells (side by side, as A* steps)? The grid as it stands without the
        port stubs `route` blocked for this search, labelled once per grid state (successive loops' searches see
        the same one)."""
        from hashlib import blake2b

        from scipy.ndimage import label

        free = ~self.grid.blocked
        for c in stubs:
            free[c] = True
        key = blake2b(np.packbits(free).tobytes(), digest_size=16).digest()
        if self._free is None or self._free[0] != key:
            self._free = (key, label(free)[0])
        marks = self._free[1]
        return bool(marks[a] and marks[a] == marks[b])

    def second(self, start, start_dir, goal, goal_dir, max_nodes: int, stubs):
        """The second search of a route the first (within max_nodes) did not find: within wide_nodes, or what the
        page has left of _MAX_NODES_PAGE if less, only when free cells join its ends; None if not run or no
        route."""
        budget = min(self.wide_nodes, self.wide_left)
        if budget <= max_nodes:
            self.wide_refused += self.wide_nodes > max_nodes  # a wider search would have run but for the budget
            return None
        if not self._joined(start, goal, stubs):
            return None
        spent = [0]
        cells = _astar(self.grid, start, start_dir, goal, goal_dir, budget, spent)
        self.wide_left -= spent[0]
        self.wide_searches += 1
        return cells

    def drawn(self, pa, ea, own_a, pb, eb, own_b, own_body=False, loose=False, under=((), ()), exempt=()):
        """The A* route from port pa to port pb as drawn, snapped to the ports or else unsnapped; None if there
        is no route or every drawing of it crosses an element (`under`: see _clear)."""
        grid = self.grid
        cells = _route_ports(grid, pa, ea, own_a, pb, eb, own_b, loose, self.max_nodes, self)
        if cells is None:
            return None
        for snap in (True, False):
            pts = _polyline(grid, cells, pa, ea, pb, eb, snap)
            if _clear(pts, (own_a, own_b), *self.elements, own_body=own_body, under=under, exempt=exempt):
                for c in cells:
                    grid.used[c] += 1
                if log.isEnabledFor(logging.DEBUG):
                    self.tight += not _clear(
                        pts, (own_a, own_b), *self.elements, grow=_INFLATE - 1e-6, under=under, exempt=exempt
                    )
                return pts
            self.rejected += 1
        return None

    def route(self, pa, ea, own_a, pb, eb, own_b, own_body=False, loose=False, under=((), ()), exempt=()):
        """(points, over labels): a drawn route, kept off text and stubs unless nothing else routes; None
        points if there is no clear route at all. Over labels: found with the keep-outs open, it runs through one of
        their boxes (keep_out). A label's margin cells by a port can end the strict search though the route drawn
        passes beside the label (random-100 C2b-2: 8RRX's D and G loops)."""
        grid = self.grid
        pts = self.drawn(pa, ea, own_a, pb, eb, own_b, own_body, loose, under, exempt)
        if pts is not None or not grid.soft.any():
            return pts, False
        grid.blocked &= ~grid.soft
        pts = self.drawn(pa, ea, own_a, pb, eb, own_b, own_body, loose, under, exempt)
        grid.blocked |= grid.soft
        return pts, pts is not None and _text_hits(pts, self.keep_out) > 0


def _cardinal(v) -> np.ndarray:
    """The page direction (right, up, left or down) nearest v."""
    x, y = float(v[0]), float(v[1])
    return np.array([np.sign(x) or 1.0, 0.0]) if abs(x) >= abs(y) else np.array([0.0, np.sign(y)])


def _over_exits(pa: Placed, pb: Placed) -> list[tuple[np.ndarray, np.ndarray]]:
    """The exits a loop between pa and pb may take. A helix drawn over a sheet is left straight in the page
    direction nearest its own (a stub across the strand its port lies on, out past the sheet's edge); when that
    direction is across the page, the stub up or down (its side of the helix's) is tried next: across, it would
    cross the next strands. Any other element is left along its axis."""

    def ways(p: Placed, out: np.ndarray) -> list[np.ndarray]:
        if not p.over:
            return [out]
        first = _cardinal(out)
        return [first] if first[0] == 0 else [first, np.array([0.0, np.sign(out[1]) or 1.0])]

    return [(ea, eb) for ea in ways(pa, pa.exit_c) for eb in ways(pb, pb.exit_n)]


def _over_exempt(layout: Layout, a: str, b: str) -> list[str]:
    """The helices drawn over a sheet that a loop between elements a and b may pass under: the one it joins, and
    any over the sheet of a strand it joins."""
    return [h for h, ids in sorted(layout.on_sheet.items()) if {a, b} & {h, *ids}]


def _strands_under(layout: Layout, sid: str, port) -> set[str]:
    """The strands of the sheet a helix is drawn over (Layout.on_sheet) whose box holds the helix's port: the only
    ones its loop's stub may cross."""
    out = set()
    for k in layout.on_sheet.get(sid, ()):
        x0, y0, x1, y1 = layout.placed[k].rect
        if x0 <= port[0] <= x1 and y0 <= port[1] <= y1:
            out.add(k)
    return out


_NAME_INK_UP, _NAME_INK_DOWN = 0.46, 0.30  # em: a helix name's ink above and below its centre (Liberation Sans)
_NAME_ALONG = (0.0, -0.25, 0.25, -0.4, 0.4)  # where along its helix a name that must move may sit: helix lengths


def _name_ink(layout: Layout, p: Placed, em: float) -> tuple[float, float, float, float]:
    """A helix name's drawn ink as a box: its width from the text's own metrics (_text_width, regular), its height
    from the digits' and Greek letters' ascent and descent. The keep-out box (_element_label_box) is a rough glyph
    size, smaller than the text on a wide page, where the text is held at 4 pt or more."""
    cx, cy = helix_label_pos(p)
    hw = _text_width(p.label, bold=False) * em / 2
    return (cx - hw, cy - _NAME_INK_DOWN * em, cx + hw, cy + _NAME_INK_UP * em)


def _clear_helix_names(layout: Layout, sses: list[SSE], loops: list[Loop]) -> None:
    """A helix name that a line runs through or nearer than figcheck.LINE_TEXT_CLEARANCE_PX of its ink (a loop's,
    a nucleic-acid strand's, a nucleotide trace's; a loop out of the helix's own end included) moves to the spot that
    lines run nearest least. Routed loops keep clear of the names already (random-100 C6: their keep-outs are sized
    from the drawn text, label_boxes and layout.glyph, at the 6 pt floor); this is the fallback for lines no route
    holds: a fallback curve where no route exists, a nucleic-acid path or trace. It moves a name to the spot lines
    run nearest least: the other side of the helix, or along it (_NAME_ALONG), the earlier spot (the current side
    first) on a tie. Lines through the name's ink (within layout.line_band: as figcheck's A6 reads a line by a text)
    count before lines near it: on a page of small text the room is wider than the name and every spot may have a
    line within it, and a spot with fewer such lines, one of them through the name, is no better (batch C2b-1, 6YF9
    focus). A spot must be on the page and meet no element, other name, stub or letter. A name no spot
    improves stays. A mirror copy's name (Layout.mirror_pairs) follows its mate's, as _mirror_sides sets it, and both
    must be clear. Run after the strands are routed (route_nucleic): their paths and traces count. Routes are
    unchanged."""
    from .layout import _mirror_sides, line_band, line_clearance, lines_meeting, segments_of

    names = sorted((p for p in layout.placed.values() if p.sse.kind in ("H", "G") and p.label), key=lambda p: p.sse.id)
    if not names:
        return
    em = _ink_em(layout, layout.font_scale)
    grow = line_clearance(layout)  # clearance and a line's half width
    band = line_band(layout)  # a line through the name, as A6 reads it
    segs = segments_of([lp.points for lp in loops] + [path.points for path in layout.na_paths])
    traces = [np.asarray(pts, float).reshape(-1, 2) for _, pts in sorted(layout.na_traces.items())]
    if not len(segs) and not traces:
        return
    solid = [q.rect for q in [*layout.placed.values(), *layout.ghosts, *layout.stems]]
    x0, y0, x1, y1 = layout.bounds
    texts = dict(label_boxes(layout, sses))
    mate_of = dict(layout.mirror_pairs)  # copy id -> mate id
    copy_of = {m: c for c, m in mate_of.items()}

    def struck(box) -> tuple[int, int]:
        """(lines through the ink, lines within the clearance of it); a nucleotide's bead and letter count as both."""
        on = _trace_hits(box, traces) if traces else 0
        through, near = (int(lines_meeting(box, segs, r).sum()) + on for r in (band, grow))
        return through, near

    def strikes(group) -> tuple[int, int]:
        got = [struck(_name_ink(layout, q, em)) for q in group]
        return sum(t for t, _ in got), sum(n for _, n in got)

    def place(group, side, along, taken):
        """The group's names at (side, along): their strikes, or None when a name is off the page or on something."""
        for q in group:
            q.label_side, q.label_along = side, float(_snap(along * q.length))
        if len(group) > 1:
            _mirror_sides([(group[1], group[0])])
        for q in group:
            ink = _name_ink(layout, q, em)
            if not (ink[0] >= x0 and ink[1] >= y0 and ink[2] <= x1 and ink[3] <= y1) or overlaps_any(ink, taken):
                return None
        return strikes(group)

    moved = 0  # names moved, a mirror copy's counted with its mate's
    k = glyph(layout)[0] / _GLYPH_W
    for p in names:
        if p.sse.id in mate_of:  # a mirror copy goes with its mate
            continue
        group = [p] + [layout.placed[copy_of[p.sse.id]]] if p.sse.id in copy_of else [p]
        now = strikes(group)
        if not now[1]:
            continue
        keys = {f"label:{q.sse.id}" for q in group}
        taken = solid + [box for name, box in texts.items() if name not in keys]
        start = [(q.label_side, q.label_along) for q in group]
        best = (now, start)  # the fewest lines struck, the earliest spot on a tie
        for along in _NAME_ALONG:
            for side in (start[0][0], -start[0][0]):
                got = place(group, side, along, taken)
                if got is not None and got < best[0]:
                    best = (got, [(q.label_side, q.label_along) for q in group])
        moved += len(group) if best[1] != start else 0
        for q, (side, along) in zip(group, best[1]):
            q.label_side, q.label_along = side, along
            texts[f"label:{q.sse.id}"] = _element_label_box(q, k)  # as label_boxes sized it (layout.glyph)
    log.debug("%d helix names moved off lines", moved)  # random-100 C6: what the campaign reports of the mover


def route_loops(layout: Layout, sses: list[SSE], bb) -> list[Loop]:
    """The loop between each pair of sequence-adjacent elements, routed around everything placed. When the layout
    holds nucleic acid it also routes the single strands (route_na.route_nucleic): that sets layout.na_paths,
    na_ends and na_letters and widens layout.bounds to hold them, so call it once per layout (a second call
    routes again on top of the widened page)."""
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
    keep_out = label_boxes(layout, sses, route=True)
    grid = _Grid(layout, [b for _, b in keep_out])
    # the elements every drawn loop must not cross, then the text and stubs a kink link also keeps clear of
    elements = _element_obstacles(layout)
    owners = elements[1] + [n.removeprefix("label:") if n.startswith("label:") else None for n, _ in keep_out]
    rects = np.hstack([elements[0], np.array([r for _, r in keep_out], float).reshape(-1, 4).T])
    router = _Router(grid, elements, [b for _, b in keep_out])

    dist = lambda ab: math.dist(layout.placed[ab[0].id].c_port, layout.placed[ab[1].id].n_port)  # noqa: E731
    routed: dict[tuple[str, str], tuple[list, bool, bool]] = {}  # points, fallback curve, over labels
    kinks_rejected = 0  # debug count: kink links not drawn straight for something (a third element, text) in the way
    for a, b in sorted(pairs, key=dist):
        pa, pb = layout.placed[a.id], layout.placed[b.id]
        gap = np.subtract(pb.n_port, pa.c_port)
        # a kink is drawn straight with no grid clearance, so it clears by the full _INFLATE (drawn routes use _TOUCH)
        kinked = b.start - a.end - 1 <= KINK_GAP or (a.id, b.id) in getattr(layout, "kink_links", ())  # a copy's
        # pieces drawn as their mates', a kink there (layout._regroup_copies)
        if kinked and np.hypot(*gap) < 1.5 and gap @ pa.exit_c > 0:
            if _clear(  # and labels
                [pa.c_port, pb.n_port],
                (pa, pb),
                rects,
                owners,
                grow=_INFLATE,
                kink=True,
                exempt=_over_exempt(layout, a.id, b.id),
            ):
                for c in grid.cells_along(pa.c_port, pb.n_port):
                    grid.used[c] += 1
                routed[(a.id, b.id)] = ([pa.c_port, pb.n_port], False, False)  # a kink: the helix just carries on
                continue
            kinks_rejected += 1
        # a helix over a sheet: its loops leave it straight across the strand its port lies on, and reach the other
        # end under it
        under = (
            _strands_under(layout, a.id, pa.c_port) | ({b.id} if pb.over else set()),
            _strands_under(layout, b.id, pb.n_port) | ({a.id} if pa.over else set()),
        )
        exempt = _over_exempt(layout, a.id, b.id)
        for h in exempt:  # its cells open to this loop
            grid.blocked &= ~grid.over[h]
        pts, over = None, False
        for ea, eb in _over_exits(pa, pb):
            pts, over = router.route(
                pa.c_port, ea, pa, pb.n_port, eb, pb, loose=(pa.over, pb.over), under=under, exempt=exempt
            )
            if pts is not None:
                break
        for h in exempt:
            grid.blocked |= grid.over[h]
        if pts is None:
            routed[(a.id, b.id)] = (_curve(pa.c_port, pa.exit_c, pb.n_port, pb.exit_n), True, False)
        else:
            routed[(a.id, b.id)] = (pts, False, over)
    if log.isEnabledFor(logging.DEBUG):
        crossing = sum(  # fallback curves through an element (other than their own two)
            fb and not _clear(pts, (layout.placed[a], layout.placed[b]), *elements)
            for (a, b), (pts, fb, _) in sorted(routed.items())
        )
        log.debug(
            "%d loops: %d fallback curves, %d of them crossing an element, %d kink links rejected, "
            "%d drawn nearer than %.2f to an element, %d drawings rejected as crossing",
            len(pairs),
            sum(fb for _, fb, _ in routed.values()),
            crossing,
            kinks_rejected,
            router.tight,
            _INFLATE,
            router.rejected,
        )
    out = []
    for a, b in pairs:
        pts, fb, over = routed[(a.id, b.id)]
        dashed = any(not bb.prev[k] for k in range(a.end + 1, b.start + 1))
        out.append(Loop(a.id, b.id, [(float(x), float(y)) for x, y in pts], dashed, fb, over))
    if layout.stack_gaps:
        _apart_in_stack_gaps(layout, router, out)
    if layout.nucleic is not None:  # after the loops, so strands see the cells they use
        from .route_na import route_nucleic

        layout.na_paths = route_nucleic(layout, router, [pts for pts, _, _ in routed.values()])
        _off_na_text(layout, router, out)
    _clear_helix_names(layout, sses, out)
    log.debug(  # fix round 1: what the second searches spent of the page's budget
        "%d second-search nodes spent of %d: %d searches, %d refused for the budget",
        _MAX_NODES_PAGE - router.wide_left,
        _MAX_NODES_PAGE,
        router.wide_searches,
        router.wide_refused,
    )
    return out


_NA_TEXT_TOLL = 25.0  # the extra cost of a cell of that text, while such a loop is routed again
_GAP_TOLL = 50.0  # the extra cost of a cell a loop must not take in a stack gap (_apart_in_stack_gaps)


def _cells(grid: _Grid, points) -> list:
    return sorted({c for p, q in zip(points, points[1:]) for c in grid.cells_along(p, q)})


def _route_again(layout: Layout, router: _Router, lp: Loop):
    """Loop `lp` routed again on the router's grid as it stands (its toll set by the caller), its own cells freed
    while it looks; (points, over a keep-out) or None, its cells given back either way (the caller marks the route
    it takes)."""
    grid = router.grid
    pa, pb = layout.placed[lp.a_id], layout.placed[lp.b_id]
    cells = _cells(grid, lp.points)
    for c in cells:
        grid.used[c] = max(grid.used[c] - 1, 0)
    under = (
        _strands_under(layout, lp.a_id, pa.c_port) | ({lp.b_id} if pb.over else set()),
        _strands_under(layout, lp.b_id, pb.n_port) | ({lp.a_id} if pa.over else set()),
    )
    exempt = _over_exempt(layout, lp.a_id, lp.b_id)
    for h in exempt:
        grid.blocked &= ~grid.over[h]
    best, over = None, False
    for ea, eb in _over_exits(pa, pb):
        best, over = router.route(
            pa.c_port, ea, pa, pb.n_port, eb, pb, loose=(pa.over, pb.over), under=under, exempt=exempt
        )
        if best is not None:
            break
    for h in exempt:
        grid.blocked |= grid.over[h]
    if best is not None:  # router.route marked its cells: give them back
        for c in _cells(grid, best):
            grid.used[c] = max(grid.used[c] - 1, 0)
    for c in cells:
        grid.used[c] += 1
    return None if best is None else (best, over)


def _in_box(points, box) -> list[tuple[float, float]]:
    """The points of a polyline, sampled every 0.05 along it, that lie inside box (x0, y0, x1, y1), rounded."""
    out = []
    for p, q in zip(points, points[1:]):
        n = max(1, int(math.dist(p, q) / 0.05))
        for t in np.linspace(0.0, 1.0, n + 1):
            x, y = p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t
            if box[0] - 1e-9 <= x <= box[2] + 1e-9 and box[1] - 1e-9 <= y <= box[3] + 1e-9:
                out.append((round(float(x), 3), round(float(y), 3)))
    return out


def _meet_in(a, b, box) -> bool:
    """Do two polylines come within 0.05 of each other inside the box?"""
    pa, pb = np.array(_in_box(a, box)).reshape(-1, 2), np.array(_in_box(b, box)).reshape(-1, 2)
    if not len(pa) or not len(pb):
        return False
    d = pa[:, None, :] - pb[None, :, :]
    return bool((np.hypot(d[..., 0], d[..., 1]) < 0.05).any())


def _apart_in_stack_gaps(layout: Layout, router: _Router, loops: list[Loop]) -> None:
    """In place: the two loops leaving the facing ends of strands drawn one above the other in a column
    (Layout.stack_gaps) are kept apart in the gap between them, each in its own half (5KQJ: B's loop up and F's loop
    down met on one line, so the column read as one strand). When they meet there, the lower strand's loop is routed
    again to turn out of the gap as low as it can (a toll on the gap's cells, and those beside it, rising with the
    height), then the upper's as high as it can and off the lower loop; the new routes are taken only when both are
    found and no longer meet in the gap (fix round 1). Each leaves the gap on the side its route takes: where a
    neighbouring strand spans the gap's height on one side (5KQJ's E, left of B and F), both leave by the other."""
    grid = router.grid
    for x0, y0, x1, y1 in layout.stack_gaps:
        cx = (x0 + x1) / 2
        ends = {}
        for k, lp in enumerate(loops):
            for p in (lp.points[0], lp.points[-1]):
                if abs(p[0] - cx) < 1e-6 and (abs(p[1] - y0) < 1e-6 or abs(p[1] - y1) < 1e-6):
                    ends["lower" if abs(p[1] - y0) < 1e-6 else "upper"] = k
        if len(ends) < 2 or loops[ends["lower"]].fallback or loops[ends["upper"]].fallback:
            continue
        lo, up = loops[ends["lower"]], loops[ends["upper"]]
        if not _meet_in(lo.points, up.points, (x0, y0, x1, y1)):
            continue
        # the lower loop turns out of the gap as low as it can, the upper as high: a toll on the gap's cells (and
        # beside it, where the loops turn out) rising with the height, or with the depth
        # their routes may start nearer the facing ends, so a tight gap holds both turns
        grid.near_ports = {(round(cx, 6), round(y0, 6)), (round(cx, 6), round(y1, 6))}
        near = grid._mask([(x0 - (x1 - x0), y0, x1 + (x1 - x0), y1)])
        rows = np.broadcast_to(grid.ys[None, :], near.shape)
        # clipped at 0: cells beside the gap reach past its bottom and top, where a negative toll paid A* to curl
        grid.toll = np.where(near, _GAP_TOLL * np.clip((rows - y0) / max(y1 - y0, 1e-9), 0.0, None), 0.0)
        new_lo = _route_again(layout, router, lo)
        if new_lo is None:
            grid.toll, grid.near_ports = None, set()
            continue
        grid.toll = np.where(near, _GAP_TOLL * np.clip((y1 - rows) / max(y1 - y0, 1e-9), 0.0, None), 0.0)
        for c in _cells(grid, new_lo[0]):
            grid.toll[c] += _GAP_TOLL
        new_up = _route_again(layout, router, up)
        grid.toll, grid.near_ports = None, set()
        if new_up is None or _meet_in(new_lo[0], new_up[0], (x0, y0, x1, y1)):
            continue
        for old, (pts, over), k in ((lo, new_lo, ends["lower"]), (up, new_up, ends["upper"])):
            for c in _cells(grid, old.points):
                grid.used[c] = max(grid.used[c] - 1, 0)
            for c in _cells(grid, pts):
                grid.used[c] += 1
            loops[k] = Loop(old.a_id, old.b_id, [(float(x), float(y)) for x, y in pts], old.dashed, False, over)


def _text_hits(points, boxes) -> int:
    """How many of `boxes` the polyline runs through."""
    pts = np.asarray(points, float).reshape(-1, 2)
    if len(pts) < 2 or not len(boxes):
        return 0
    return int((clip_lengths(np.hstack([pts[:-1], pts[1:]]), boxes) > 0).any(axis=0).sum())


_REROUTE_TURNS = 4  # corners a route off NA text may add to the one it replaces when it only grazed the clearance
_REROUTE_TURNS_INK = 12  # ... and when it ran through the text's ink (7XSP's B:355-374 loop: 6 more to free 2 letters)


def _strand_crossings(points, strands) -> int:
    """How many times a polyline properly crosses the nucleic-acid strands (`strands`, an na_paths._Strands)."""
    return strands.crossings(points) if strands is not None else 0


def _better(new, old, boxes, text, strands=None, ink=()) -> bool:
    """Does route `new` cross less nucleic-acid text (`boxes`) than `old`, and no more text in all (`text`: the
    names and stubs routes keep out of), with at most _REROUTE_TURNS more corners (_REROUTE_TURNS_INK when `old` runs
    through the text's ink itself, `ink`, not only its clearance), crossing the strands (`strands`:
    routed round the old route) no more often than it crosses less text (fix round 1: one crossing more per letter or
    label freed at most)?"""
    turns = _REROUTE_TURNS_INK if _text_hits(old, ink) else _REROUTE_TURNS
    if len(new) > len(old) + turns:
        return False
    freed = _text_hits(old, boxes) - _text_hits(new, boxes)
    if _strand_crossings(new, strands) - _strand_crossings(old, strands) > freed:
        return False
    fewer = _text_hits(new, boxes) < _text_hits(old, boxes)
    return fewer and _text_hits(new, [*boxes, *text]) < _text_hits(old, [*boxes, *text])


def _off_na_text(layout: Layout, router: _Router, loops: list[Loop]) -> None:
    """In place: a protein loop crossing the nucleic acid's text (its letters and 5'/3' labels, placed only after the
    loops were routed: layout.na_text_boxes) is routed again with that text as a keep-out, and takes the new route
    when it crosses less of it (random-100 C11: 6SO9's RRM loops ran over the bound RNA's letters). Fallback
    curves and kinks stay as they are; a loop with no better route keeps its own."""
    from .layout import line_clearance, na_text_boxes

    g = line_clearance(layout)  # a protein loop keeps this far from the nucleic acid's text, beyond its ink
    ink = na_text_boxes(layout)
    boxes = [(x0 - g, y0 - g, x1 + g, y1 + g) for x0, y0, x1, y1 in ink]
    hit = [k for k, lp in enumerate(loops) if not lp.fallback and len(lp.points) > 2 and _text_hits(lp.points, boxes)]
    if not hit:
        return
    from .na_paths import _Strands

    grid = router.grid
    text = grid._mask(boxes)
    grid.toll = np.where(text, _NA_TEXT_TOLL, 0.0)
    strands = _Strands(None)  # the nucleic acid's strands as drawn: a reroute crosses no more of them
    for path in layout.na_paths:
        strands.add(path.points)
    for k in hit:
        lp = loops[k]
        pa, pb = layout.placed[lp.a_id], layout.placed[lp.b_id]
        cells = sorted({c for p, q in zip(lp.points, lp.points[1:]) for c in grid.cells_along(p, q)})
        for c in cells:  # its own cells free while it looks again
            grid.used[c] = max(grid.used[c] - 1, 0)
        under = (  # as route_loops routed it: a helix over a sheet leaves it across its strands
            _strands_under(layout, lp.a_id, pa.c_port) | ({lp.b_id} if pb.over else set()),
            _strands_under(layout, lp.b_id, pb.n_port) | ({lp.a_id} if pa.over else set()),
        )
        exempt = _over_exempt(layout, lp.a_id, lp.b_id)
        for h in exempt:
            grid.blocked &= ~grid.over[h]
        best, over = None, False
        for ea, eb in _over_exits(pa, pb):
            best, over = router.route(
                pa.c_port, ea, pa, pb.n_port, eb, pb, loose=(pa.over, pb.over), under=under, exempt=exempt
            )
            if best is not None:
                break
        for h in exempt:
            grid.blocked |= grid.over[h]
        better = best is not None and _better(best, lp.points, boxes, router.keep_out, strands, ink)
        if better and over <= lp.crosses_labels:
            loops[k] = Loop(lp.a_id, lp.b_id, [(float(x), float(y)) for x, y in best], lp.dashed, False, over)
            continue
        if best is not None:  # not taken: give back the cells drawn() marked
            for c in sorted({c for p, q in zip(best, best[1:]) for c in grid.cells_along(p, q)}):
                grid.used[c] = max(grid.used[c] - 1, 0)
        for c in cells:
            grid.used[c] += 1
    grid.toll = None
