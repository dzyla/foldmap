# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Drawn nucleic-acid strands: the NAPath and NATie records, the plane geometry of polylines and stem bodies that
routing (na_runs), free ends and gaps (na_ends) and the text beside them (na_labels) share, the clearance check of a
strand among the elements, and the registry of the strands drawn so far (_Strands)."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .layout import _GLYPH_W, PlacedStem
from .route import _clear, _Grid, _simplify, clip_len


@dataclass(eq=False)
class NAPath:
    """A drawn stretch of single strand. kind: run (stem to stem), hairpin (an arc over a stem's tip), joint
    (between two pieces of a bent duplex), tail (a strand's free end), trace (a strand with no pairs) or gap
    (nucleotides cropped away, dashed; labelled "… N nt …" unless fewer than na_ends.SHORT_GAP). nts: the nucleotide
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


def _crosses_segment(pts, a, b) -> bool:
    """Does the polyline cross segment a-b (properly, not merely touch)?"""

    def turn(o, u, v) -> float:
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])

    for p, q in zip(pts, pts[1:]):
        if turn(p, q, a) * turn(p, q, b) < -1e-12 and turn(a, b, p) * turn(a, b, q) < -1e-12:
            return True
    return False


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


def graze(p, q, d) -> float:
    """How deep segment p-q runs into the stem's (or element's) true, tilted body, page units: 0 when it misses."""
    (u0, v0), (u1, v1) = _local(p, d), _local(q, d)
    t = np.linspace(0.0, 1.0, 33)
    u, v = u0 + t * (u1 - u0), v0 + t * (v1 - v0)
    return float(max(0.0, np.max(np.minimum(d.length / 2 - np.abs(u), d.width / 2 - np.abs(v)))))


def among(path: NAPath, at: dict) -> list[tuple[float, float]]:
    """The part of a drawn strand that lies among the other elements. A run or a joint (stem to stem) is drawn as
    [rail point, end face, ..., end face, rail point]: without its two rail stretches, its first and last segments.
    Any other strand without its first and last stretches inside its own stems (outside_own). `at`: stem nucleotide
    -> its PlacedStem."""
    if path.kind in ("run", "joint") and len(path.between) == 2:
        return list(path.points[1:-1])
    return outside_own(path.points, tuple(at[k] for k in path.between if k in at))


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
