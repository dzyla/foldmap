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
from .layout import Layout, Placed, label_boxes
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
_PAD = 3.0
_MAX_NODES = 200_000
_MAX_CELLS = 1_500_000
_DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
_PUSH = [0.4 + 0.25 * k for k in range(12)]  # how far out along the exit the path may start
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
        for t in _PUSH_OVER if over else _PUSH:
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


def _astar(grid: _Grid, start, start_dir, goal, goal_dir, max_nodes: int = _MAX_NODES):
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
        if nodes > max_nodes:
            return None
        for nd in range(4):
            if nd == (d + 2) % 4:
                continue
            nx_, ny_ = x + _DIRS[nd][0], y + _DIRS[nd][1]
            if not grid.free((nx_, ny_)):
                continue
            cost = g + 1 + (_TURN if nd != d else 0) + (_REUSE if grid.used[nx_, ny_] else 0)
            cost += _SLAB_COST if grid.slab[nx_, ny_] else 0
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


def _route_ports(grid: _Grid, pa, ea, own_a, pb, eb, own_b, loose=False, max_nodes: int = _MAX_NODES):
    """Grid cells of the cheapest orthogonal route from port pa (leaving along ea) to port pb (entered against
    eb), or None; own_a and own_b (anything with a rect) are what the ports belong to. `loose`: see
    _Grid.free_cell_along; a pair (for pa, for pb) or one value for both."""
    la, lb = loose if isinstance(loose, tuple) else (loose, loose)
    start = grid.free_cell_along(pa, ea, own_a, la)
    goal = grid.free_cell_along(pb, eb, own_b, lb)
    if start is None or goal is None:
        return None
    stubs = grid.stub_cells(pa, start) | grid.stub_cells(pb, goal)  # the straight ends drawn to each port
    for c in stubs:
        grid.blocked[c] = True
    cells = _astar(grid, start, _axis(ea), goal, _axis(-np.asarray(eb)), max_nodes)
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

    def drawn(self, pa, ea, own_a, pb, eb, own_b, own_body=False, loose=False, under=((), ()), exempt=()):
        """The A* route from port pa to port pb as drawn, snapped to the ports or else unsnapped; None if there
        is no route or every drawing of it crosses an element (`under`: see _clear)."""
        grid = self.grid
        cells = _route_ports(grid, pa, ea, own_a, pb, eb, own_b, loose, self.max_nodes)
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
        points if there is no clear route at all."""
        grid = self.grid
        pts = self.drawn(pa, ea, own_a, pb, eb, own_b, own_body, loose, under, exempt)
        if pts is not None or not grid.soft.any():
            return pts, False
        grid.blocked &= ~grid.soft
        pts = self.drawn(pa, ea, own_a, pb, eb, own_b, own_body, loose, under, exempt)
        grid.blocked |= grid.soft
        return pts, pts is not None


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
        if b.start - a.end - 1 <= KINK_GAP and np.hypot(*gap) < 1.5 and gap @ pa.exit_c > 0:
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
    if layout.nucleic is not None:  # after the loops, so strands see the cells they use
        from .route_na import route_nucleic

        layout.na_paths = route_nucleic(layout, router, [pts for pts, _, _ in routed.values()])
    return out
