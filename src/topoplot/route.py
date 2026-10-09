"""Orthogonal loop routing: A* on a coarse grid around the placed elements."""

from __future__ import annotations

import heapq
import math
from dataclasses import dataclass

import numpy as np

from .features import KINK_GAP
from .layout import Layout, Placed, label_boxes
from .model import SSE

RES = 0.25  # grid cell size, page units
_INFLATE = 0.15  # keep loops this far from arrows
_TURN = 2.5
_REUSE = 6.0  # extra cost for a cell an earlier loop already uses
_PAD = 3.0
_MAX_NODES = 200_000
_MAX_CELLS = 1_500_000
_DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
_PUSH = [0.4 + 0.25 * k for k in range(12)]  # how far out along the exit the path may start


@dataclass(eq=False)
class Loop:
    a_id: str  # source element (loop leaves its C-terminal end)
    b_id: str  # target element (loop enters its N-terminal end)
    points: list[tuple[float, float]]
    dashed: bool  # residues missing between the two elements
    fallback: bool  # no orthogonal route found; points trace a smooth curve
    crosses_labels: bool = False  # only a route over a label or chain-end stub was found


class _Grid:
    def __init__(self, layout: Layout, keep_out: list[tuple[float, float, float, float]] = ()):
        x0, y0, x1, y1 = layout.bounds
        w, h = x1 - x0 + 2 * _PAD, y1 - y0 + 2 * _PAD
        self.res = max(RES, math.sqrt(w * h / _MAX_CELLS))
        self.x0, self.y0 = x0 - _PAD, y0 - _PAD
        self.nx, self.ny = int(math.ceil(w / self.res)), int(math.ceil(h / self.res))
        self.used = np.zeros((self.nx, self.ny), np.int32)
        xs = self.x0 + (np.arange(self.nx) + 0.5) * self.res
        ys = self.y0 + (np.arange(self.ny) + 0.5) * self.res
        def mask(rects):
            m = np.zeros((self.nx, self.ny), bool)
            for rx0, ry0, rx1, ry1 in rects:
                m |= np.outer((xs > rx0 - _INFLATE) & (xs < rx1 + _INFLATE), (ys > ry0 - _INFLATE) & (ys < ry1 + _INFLATE))
            return m

        self.blocked = mask(p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.dna])
        self.soft = mask(keep_out) & ~self.blocked  # text and stubs: avoided, given up only if nothing else routes
        self.blocked |= self.soft

    def cell(self, pt) -> tuple[int, int]:
        return int(math.floor((pt[0] - self.x0) / self.res)), int(math.floor((pt[1] - self.y0) / self.res))

    def center(self, c) -> tuple[float, float]:
        return self.x0 + (c[0] + 0.5) * self.res, self.y0 + (c[1] + 0.5) * self.res

    def free(self, c) -> bool:
        return 0 <= c[0] < self.nx and 0 <= c[1] < self.ny and not self.blocked[c]

    def stub_cells(self, port, end) -> set[tuple[int, int]]:
        """Free cells under the straight stub from a port to its first routing cell (that cell excluded),
        so no other part of the same loop can pass between the two and be crossed by the stub."""
        ex, ey = self.center(end)
        n = max(2, int(math.dist(port, (ex, ey)) / (self.res / 4)))
        out = set()
        for t in np.linspace(0.0, 1.0, n):
            c = self.cell((port[0] + (ex - port[0]) * t, port[1] + (ey - port[1]) * t))
            if c != end and self.free(c):
                out.add(c)
        return out

    def free_cell_along(self, port, exit_vec, own: Placed):
        """First free cell on the way out of a port. Only the element's own bounding box (a tilted
        helix end can sit inside it) may be skipped; any other obstacle in the way means no route."""
        rx0, ry0, rx1, ry1 = own.rect
        for t in _PUSH:
            pt = (port[0] + exit_vec[0] * t, port[1] + exit_vec[1] * t)
            c = self.cell(pt)
            if self.free(c):
                return c
            cx, cy = self.center(c)
            if not (rx0 - _INFLATE < cx < rx1 + _INFLATE and ry0 - _INFLATE < cy < ry1 + _INFLATE):
                return None
        return None


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
    pa, pb = a.c_port, b.n_port
    start = grid.free_cell_along(pa, a.exit_c, a)
    goal = grid.free_cell_along(pb, b.exit_n, b)
    if start is None or goal is None:
        return None
    stubs = grid.stub_cells(pa, start) | grid.stub_cells(pb, goal)  # the straight ends drawn to each port
    for c in stubs:
        grid.blocked[c] = True
    cells = _astar(grid, start, _axis(a.exit_c), goal, _axis(-b.exit_n))
    for c in stubs:
        grid.blocked[c] = False
    if cells is None:
        return None
    for c in cells:
        grid.used[c] += 1
    pts = [list(grid.center(c)) for c in _corners(cells)]
    _snap_to_port(pts, pa, a.exit_c, True)
    _snap_to_port(pts, pb, b.exit_n, False)
    return _simplify([pa] + [(x, y) for x, y in pts] + [pb])


def route_loops(layout: Layout, sses: list[SSE], bb) -> list[Loop]:
    pairs = [(a, b) for a, b in zip(sses, sses[1:]) if a.chain == b.chain and a.id in layout.placed and b.id in layout.placed]
    if not pairs:
        return []
    grid = _Grid(layout, [b for _, b in label_boxes(layout, sses)])
    dist = lambda ab: math.dist(layout.placed[ab[0].id].c_port, layout.placed[ab[1].id].n_port)  # noqa: E731
    routed: dict[tuple[str, str], tuple[list, bool]] = {}
    for a, b in sorted(pairs, key=dist):
        pa, pb = layout.placed[a.id], layout.placed[b.id]
        gap = np.subtract(pb.n_port, pa.c_port)
        if b.start - a.end - 1 <= KINK_GAP and np.hypot(*gap) < 1.5 and gap @ pa.exit_c > 0:
            routed[(a.id, b.id)] = ([pa.c_port, pb.n_port], False, False)  # a kink: the helix just carries on
            continue
        pts, over = _route_pair(grid, pa, pb), False
        if pts is None and grid.soft.any():
            grid.blocked &= ~grid.soft
            pts, over = _route_pair(grid, pa, pb), True
            grid.blocked |= grid.soft
        if pts is None:
            routed[(a.id, b.id)] = (_curve(pa.c_port, pa.exit_c, pb.n_port, pb.exit_n), True, False)
        else:
            routed[(a.id, b.id)] = (pts, False, over)
    out = []
    for a, b in pairs:
        pts, fb, over = routed[(a.id, b.id)]
        dashed = any(not bb.prev[k] for k in range(a.end + 1, b.start + 1))
        out.append(Loop(a.id, b.id, [(float(x), float(y)) for x, y in pts], dashed, fb, over))
    return out
