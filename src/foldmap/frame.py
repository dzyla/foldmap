# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Global page frame: which 3D directions become page-right (u), page-up (v) and the view axis (w)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import SSE


@dataclass(eq=False)
class Frame:
    origin: np.ndarray
    u: np.ndarray  # page x
    v: np.ndarray  # page y (up)
    w: np.ndarray  # towards the viewer
    sector: dict | None = None  # chain -> protomer index, for a dimer drawn as a mirror image (mirror_x)
    classes: dict | None = None  # chain -> chain class (its position within a protomer)
    mirror_x: float | None = None  # page x of the line protomer 1 mirrors protomer 0 about (a dyad in the page)

    def project(self, points: np.ndarray, chains=None) -> np.ndarray:
        rel = np.atleast_2d(points) - self.origin
        return np.stack([rel @ self.u, rel @ self.v], axis=-1)

    def depth(self, points: np.ndarray) -> np.ndarray:
        return (np.atleast_2d(points) - self.origin) @ self.w

    def basis(self, point) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Page-right, page-up and toward-viewer directions at a point (the same everywhere here)."""
        return self.u, self.v, self.w

    def direction(self, point, vec) -> np.ndarray:
        """A 3D direction at `point`, as seen on the page (page-right, page-up components)."""
        u, v, _ = self.basis(point)
        return np.array([float(vec @ u), float(vec @ v)])


@dataclass(eq=False)
class CylinderFrame:
    """Unrolled view of a cyclic assembly: the symmetry axis points up the page and every protomer is seen
    from outside, square-on, at its own angle around the axis. Protomer k sits `period` Å to the right of
    protomer k-1 and is drawn exactly like it. Locally a proper rotation (never a mirror)."""

    origin: np.ndarray  # a point on the axis
    axis: np.ndarray  # page up
    ref: np.ndarray  # unit vector perpendicular to the axis, toward protomer 0
    n: int
    radius: float  # Å; arc length per radian
    sector: dict  # chain -> protomer index
    classes: dict = None  # chain -> chain class (its position within a protomer)
    angles: list = None  # each protomer's turn about the axis (radians, unwrapped); default 2*pi*k/n
    rise: float = 0.0  # Å along the axis from one protomer to the next (helical filaments)

    def __post_init__(self):
        if self.angles is None:
            self.angles = [2 * np.pi * k / self.n for k in range(self.n)]

    @property
    def period(self) -> float:
        """Page distance (Å) from one protomer to the next."""
        steps = np.diff(self.angles)
        return self.radius * float(np.mean(steps)) if len(steps) else self.radius * 2 * np.pi

    def _theta(self, p: np.ndarray) -> float:
        rel = p - self.origin
        rho = rel - (rel @ self.axis) * self.axis
        side = np.cross(self.axis, self.ref)
        return float(np.arctan2(rho @ side, rho @ self.ref))

    def _place(self, p: np.ndarray, k: int | None) -> np.ndarray:
        th = self._theta(p)
        wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi  # noqa: E731
        if k is None:  # chains outside any protomer: the protomer they face
            k = int(np.argmin([abs(wrap(th - a)) for a in self.angles]))
        local = wrap(th - self.angles[k])
        x = self.radius * (self.angles[k] - float(np.mean(self.angles))) + self.radius * local
        return np.array([x, float((p - self.origin) @ self.axis)])

    def project(self, points: np.ndarray, chains=None) -> np.ndarray:
        pts = np.atleast_2d(points)
        ks = [None] * len(pts) if chains is None else [self.sector.get(c) for c in chains]
        return np.array([self._place(p, k) for p, k in zip(pts, ks)])

    def basis(self, point) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rel = np.asarray(point, float) - self.origin
        rho = rel - (rel @ self.axis) * self.axis
        out = rho / np.linalg.norm(rho) if np.linalg.norm(rho) > 1e-6 else self.ref
        right = np.cross(self.axis, out)  # e_theta: page-right; (right, up, out) is right-handed
        return right, self.axis, out

    def direction(self, point, vec) -> np.ndarray:
        u, v, _ = self.basis(point)
        return np.array([float(vec @ u), float(vec @ v)])


def _turn3(vec: np.ndarray, axis: np.ndarray, a: float) -> np.ndarray:
    """`vec` turned by `a` radians about the unit `axis` (right-handed: anticlockwise seen from the axis' tip)."""
    vec = np.asarray(vec, float)
    return (
        vec * np.cos(a)
        + np.cross(axis, vec) * np.sin(a)
        + np.outer(vec @ axis, axis).reshape(vec.shape) * (1 - np.cos(a))
    )


def turn2(points: np.ndarray, a: float) -> np.ndarray:
    """Page points (or one point) turned by `a` radians anticlockwise about the page origin."""
    c, s = np.cos(a), np.sin(a)
    p = np.asarray(points, float)
    return np.stack([c * p[..., 0] - s * p[..., 1], s * p[..., 0] + c * p[..., 1]], axis=-1)


@dataclass(eq=False)
class RosetteFrame:
    """A cyclic assembly (Cn, n >= 3) seen down its n-fold axis (Style.cyclic rosette, plan 21): a point of protomer
    k is turned by -angles[k] about the axis (onto protomer 0), projected down the axis with page-up toward protomer
    0, then turned by +angles[k] about the page centre, where the axis projects (page origin). angles[k] is exactly
    k * 2 pi / n for the copy at that place round the ring, so copy k is drawn as protomer 0 turned. Points of chains
    outside every protomer are projected as they lie. Locally a proper rotation (never a mirror).

    A D2 (dihedral_rosette_frame, plan 21 Task 3) has twice the copies: copies 0..order-1 are ring A, the halves of
    the protomers one two-fold (`dyad`) takes onto the other halves, seen down the main axis; copy order + j is ring
    B, the image of ring A's copy j under that two-fold, drawn mirrored in the page's dyad line (page y -> -y), turned
    by angles[order + j] and moved below ring A (layout._rosette). Projecting down the axis, a half-turn about a
    page-right two-fold is that mirror: angles[order + j] is -angles[j] with the dyad page-right, minus a further
    half-turn with the other two-fold across the axis page-right (dihedral_rosette_frame)."""

    origin: np.ndarray  # a point on the axis: the page centre
    axis: np.ndarray  # toward the viewer (w)
    ref: np.ndarray  # page up: perpendicular to the axis, toward protomer 0
    n: int  # copies drawn
    sector: dict  # chain -> copy index (its protomer's; a D2's half-protomer's)
    classes: dict = None  # chain -> chain class (its position within a protomer)
    angles: list = None  # each copy's exact page turn about the axis, radians (k * 2 pi / n)
    mirrors: list = None  # each copy drawn mirrored in the dyad line before its turn (a D2's ring B)
    dyad: np.ndarray | None = None  # a D2's two-fold taking ring A's copy j onto ring B's (local frames)

    def __post_init__(self):
        if self.angles is None:
            self.angles = [2 * np.pi * k / self.n for k in range(self.n)]
        if self.mirrors is None:
            self.mirrors = [False] * self.n

    @property
    def order(self) -> int:
        """The main axis' order: the copies in one ring (n, or half of a D2's)."""
        return self.n // 2 if self.dyad is not None else self.n

    @property
    def u(self) -> np.ndarray:  # page right: u x v = w
        return np.cross(self.ref, self.axis)

    @property
    def v(self) -> np.ndarray:
        return self.ref

    @property
    def w(self) -> np.ndarray:
        return self.axis

    def local(self, k: int = 0) -> Frame:
        """The plain frame protomer k is seen in as protomer 0 is: the page axes turned with it about the axis (its
        elements' places and angles in that frame are protomer 0's, before the page turn)."""
        from .layout import _snap

        if self.mirrors[k]:  # ring B's copy j: ring A's copy j taken by the half-turn about the dyad
            a = self.angles[k - self.order]
            d = self.dyad

            def move(x):
                t = _turn3(x, self.axis, a)
                return 2 * (t @ d) * d - t

            u, v, w = (_snap(move(x)) for x in (self.u, self.v, self.axis))
            return Frame(self.origin, u, v, w)
        a = self.angles[k]
        u, v = (_snap(_turn3(x, self.axis, a)) for x in (self.u, self.v))  # cos, sin: the last bit differs by libm
        return Frame(self.origin, u, v, self.axis)

    def project(self, points: np.ndarray, chains=None) -> np.ndarray:
        pts = np.atleast_2d(points)
        ks = [None] * len(pts) if chains is None else [self.sector.get(c) for c in chains]
        out = []
        for p, k in zip(pts, ks):
            a = 0.0 if k is None else self.angles[k]
            rel = _turn3(p - self.origin, self.axis, -a)
            out.append(turn2(np.array([rel @ self.u, rel @ self.v]), a))
        return np.array(out).reshape(-1, 2)

    def depth(self, points: np.ndarray) -> np.ndarray:
        return (np.atleast_2d(points) - self.origin) @ self.axis

    def basis(self, point) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Page-right, page-up and toward-viewer directions (one turn turns both a point and the page back)."""
        return self.u, self.v, self.axis

    def direction(self, point, vec) -> np.ndarray:
        u, v, _ = self.basis(point)
        return np.array([float(vec @ u), float(vec @ v)])


def rosette_frame(sses: list[SSE], sym) -> RosetteFrame:
    """The rosette of a Cn assembly (n >= 3): seen from the tip of its axis, page-up toward protomer 0's elements.
    Each protomer's turn is the nearest whole k * 2 pi / n to the turn symmetry measured (Symmetry.angles). The axis
    and page-up are snapped before they choose anything (they differ in the last bit between BLAS kernels)."""
    from .layout import _snap

    axis = _snap(_unit(sym.axis))
    own = [s.centroid for s in sses if s.chain in sym.protomers[0]]
    first = np.mean(own, axis=0) - sym.centre if own else np.zeros(3)
    ref = first - (first @ axis) * axis
    ref = _snap(_unit(ref)) if np.linalg.norm(ref) > 1e-6 else _any_perpendicular(axis)
    ref = _snap(_unit(ref - (ref @ axis) * axis))
    n = len(sym.protomers)
    step = 2 * np.pi / n
    measured = list(sym.angles) or [step * k for k in range(n)]
    angles = [step * (round(float(a) / step) % n) for a in measured]
    sector, classes = protomer_maps(sym)
    return RosetteFrame(_snap(sym.centre), axis, ref, n, sector, classes, angles)


def dihedral_rings(sym, perm: dict) -> tuple[list[list[str]], list[list[str]]] | None:
    """A D2's two rings (plan 21 Task 3), as lists of copies' chains: ring A's copy k is the part of protomer k whose
    chains the two-fold `perm` (symmetry.protomer_dyad) takes onto the rest of it (protomer 0's in order, a chain
    joining ring A unless its image already has: the same places in every protomer), ring B's copy k those chains'
    images. None unless the two rings share no chain and hold every protomer's."""
    first = sym.protomers[0]
    places, claimed = [], set()
    for i, c in enumerate(first):
        if c in claimed or perm.get(c) not in first:
            continue
        places.append(i)
        claimed |= {c, perm[c]}
    ring_a = [[p[i] for i in places] for p in sym.protomers]
    if any(c not in perm for copy in ring_a for c in copy):
        return None
    ring_b = [[perm[c] for c in copy] for copy in ring_a]
    a = [c for copy in ring_a for c in copy]
    b = [c for copy in ring_b for c in copy]
    every = sorted(c for p in sym.protomers for c in p)
    if len(set(a) | set(b)) != len(a) + len(b) or sorted(a + b) != every:
        return None
    return ring_a, ring_b


def dihedral_rosette_frame(sses: list[SSE], sym, dyad: np.ndarray, rings) -> RosetteFrame:
    """The mirrored rosette of a D2 (Style.cyclic rosette, plan 21 Task 3): ring A (dihedral_rings) seen down the main
    axis from its own side (it lies nearer the viewer), the two-fold `dyad` that takes it onto ring B along
    page-right, so ring B, seen down the same axis, is ring A mirrored in the dyad line (RosetteFrame). Ring A's
    copies turn by the nearest whole multiples of 2 pi / order to Symmetry.angles; ring B's copy j by minus ring
    A's copy j's turn, after its mirror. The axis and the dyad are snapped before they choose anything."""
    from .layout import _snap

    ring_a, ring_b = rings
    order = len(ring_a)
    axis = _snap(_unit(sym.axis))
    own = [s.centroid for s in sses if any(s.chain in copy for copy in ring_a)]
    side = 1.0 if not own or float((np.mean(own, axis=0) - sym.centre) @ axis) >= 0 else -1.0
    axis = axis * side
    d = np.asarray(dyad, float)
    d = _snap(_unit(d - (d @ axis) * axis))
    # page-right: the dyad, or the other two-fold across the axis (the dyad turned a quarter about it: the half-turn
    # about the axis, then this one, is the dyad's), whichever lies nearer ring A's first copy's direction from the
    # axis, so ring A's two copies stand side by side and the rings make two rows (layout._rosette pushes them apart
    # along page-right). Seen down the axis, ring B is ring A mirrored in either line: with the other, each copy of
    # ring B is the mirror of ring A's other copy, turned by a further half-turn.
    other = _snap(_unit(np.cross(axis, d)))
    mid = np.mean([s.centroid for s in sses if s.chain in ring_a[0]] or [sym.centre], axis=0) - sym.centre
    across = abs(float(mid @ other)) > abs(float(mid @ d)) + 1e-6
    right = other if across else d
    ref = _snap(_unit(np.cross(axis, right)))  # page-right = ref x axis
    step = 2 * np.pi / order
    measured = list(sym.angles) or [step * k for k in range(order)]
    turns = [step * (round(side * float(a) / step) % order) for a in measured]
    angles = turns + [(-(a + (np.pi if across else 0.0))) % (2 * np.pi) for a in turns]
    sector = {c: k for k, copy in enumerate(ring_a + ring_b) for c in copy}
    classes = {c: i for copy in ring_a + ring_b for i, c in enumerate(copy)}
    mirrors = [False] * order + [True] * order
    return RosetteFrame(_snap(sym.centre), axis, ref, 2 * order, sector, classes, angles, mirrors, d)


def protomer_maps(sym) -> tuple[dict, dict]:
    """(sector, classes): each protomer chain's protomer index, and its chain class (its position within its
    protomer), which _symmetry_ties matches copies by."""
    sector = {c: k for k, p in enumerate(sym.protomers) for c in p}
    classes = {c: i for p in sym.protomers for i, c in enumerate(p)}
    return sector, classes


def symmetric_frame(sses: list[SSE], sym) -> CylinderFrame:
    cents = np.array([s.centroid for s in sses])
    rel = cents - sym.centre
    rho = rel - np.outer(rel @ sym.axis, sym.axis)
    radius = float(np.median(np.linalg.norm(rho, axis=1)))
    own = [s.centroid for s in sses if s.chain in sym.protomers[0]]  # none: protomer 0 has no helix or strand
    first = np.mean(own, axis=0) - sym.centre if own else np.zeros(3)
    ref = first - (first @ sym.axis) * sym.axis
    ref = ref / np.linalg.norm(ref) if np.linalg.norm(ref) > 1e-6 else _any_perpendicular(sym.axis)
    sector, classes = protomer_maps(sym)
    return CylinderFrame(
        sym.centre,
        sym.axis,
        ref,
        len(sym.protomers),
        max(radius, 5.0),
        sector,
        classes,
        list(sym.angles) or None,
        float(sym.rise),
    )


def _unit(x: np.ndarray) -> np.ndarray:
    return x / np.linalg.norm(x)


def _any_perpendicular(v: np.ndarray) -> np.ndarray:
    e = np.eye(3)[int(np.argmin(np.abs(v)))]
    return _unit(e - (e @ v) * v)


def view_frame(
    sses: list[SSE],
    rotate: float = 0.0,
    flip_v: bool = False,
    up: list[float] | np.ndarray | None = None,
    view: list[float] | np.ndarray | None = None,
) -> Frame:
    """Page axes chosen so helices and strands lie along the page.

    Page-up is the length-weighted dominant element axis, signed so the first element points up, or
    `up` exactly as given (e.g. a membrane normal pointing to the side that belongs on top). Page-right is
    the widest spread of the element centroids perpendicular to it, signed so the last element lies to the
    right of the first. `rotate` turns
    the picture counter-clockwise (degrees); `flip_v` mirrors it (on request only).
    """
    if not sses:
        return Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    cents = np.array([s.centroid for s in sses])
    origin = cents.mean(axis=0)
    if up is not None:
        up = np.asarray(up, float)
        if up.shape != (3,) or not np.isfinite(up).all() or np.linalg.norm(up) < 1e-9:
            raise ValueError(f"up must be a non-zero 3D vector, got {up.tolist()}")
        v = _unit(up)  # taken as given: it says which side goes on top
    else:
        scatter = sum(len(s) * np.outer(s.axis, s.axis) for s in sses)
        v = np.linalg.eigh(scatter)[1][:, -1]
        if sses[0].axis @ v < 0:
            v = -v
    if view is not None:  # the viewer looks along -view: page-up is `v` made perpendicular to it
        w = np.asarray(view, float)
        if w.shape != (3,) or not np.isfinite(w).all() or np.linalg.norm(w) < 1e-9:
            raise ValueError(f"view must be a non-zero 3D vector, got {w.tolist()}")
        w = _unit(w)
        v = v - (v @ w) * w
        v = _unit(v) if np.linalg.norm(v) > 1e-6 else _any_perpendicular(w)
        u = np.cross(v, w)  # so that u x v = w: a rotation, never a mirror
        theta = np.radians(rotate)
        u, v = np.cos(theta) * u - np.sin(theta) * v, np.sin(theta) * u + np.cos(theta) * v
        return Frame(origin, u, -v if flip_v else v, w)
    rel = cents - origin
    flat = rel - np.outer(rel @ v, v)
    if np.linalg.norm(flat) < 1e-9:
        u = _any_perpendicular(v)
    else:
        u = np.linalg.svd(flat, full_matrices=False)[2][0]
        u = _unit(u - (u @ v) * v)
    if (cents[-1] - cents[0]) @ u < 0:
        u = -u
    w = np.cross(u, v)
    theta = np.radians(rotate)
    u, v = np.cos(theta) * u - np.sin(theta) * v, np.sin(theta) * u + np.cos(theta) * v
    if flip_v:
        v = -v
    return Frame(origin, u, v, w)


def _duplex_axis(nucleic) -> tuple[np.ndarray, np.ndarray]:
    """(centre, axis) of the longest stem: the duplex a protein-DNA picture is drawn around."""
    from .nucleic import stem_axis

    biggest = max(nucleic.structure.stems, key=len)
    return stem_axis(biggest, nucleic.strands)


def dna_view_frame(sses: list[SSE], nucleic, rotate: float = 0.0, flip_v: bool = False) -> Frame:
    """The usual picture of a protein-DNA complex: the longest duplex runs across the page (5' of its first
    strand on the left) and the protein sits below it. A rotation, never a mirror (unless flip_v)."""
    dna_c, u = _duplex_axis(nucleic)  # the same axis the layout draws the duplex along
    pts = np.vstack([s.p for s in nucleic.strands])
    prot = np.array([s.centroid for s in sses]).mean(axis=0) if sses else dna_c
    up = dna_c - prot  # from the protein toward the DNA
    v = up - (up @ u) * u
    v = _unit(v) if np.linalg.norm(v) > 1e-6 else _any_perpendicular(u)
    w = np.cross(u, v)
    origin = np.vstack([pts, *([np.array([s.centroid for s in sses])] if sses else [])]).mean(axis=0)
    theta = np.radians(rotate)
    u, v = np.cos(theta) * u - np.sin(theta) * v, np.sin(theta) * u + np.cos(theta) * v
    if flip_v:
        v = -v
    return Frame(origin, u, v, w)


@dataclass(eq=False)
class RadialFrame:
    """Protein-DNA view unrolled around the duplex axis: page x runs along the DNA, page y is minus the
    distance from the DNA axis. Whatever wraps around the DNA (zinc fingers, clamps) lands below it, in a
    row under the base pairs it binds, at its true distance. Locally a proper rotation."""

    origin: np.ndarray  # duplex centre
    axis: np.ndarray  # page right
    down: np.ndarray  # fallback outward direction on the axis itself (toward the protein)
    sector: dict | None = None  # as Frame: set with mirror_x for a dimer whose dyad lies in the page
    classes: dict | None = None
    mirror_x: float | None = None

    def _outward(self, p: np.ndarray) -> np.ndarray:
        rel = np.asarray(p, float) - self.origin
        rho = rel - (rel @ self.axis) * self.axis
        return rho / np.linalg.norm(rho) if np.linalg.norm(rho) > 1e-6 else self.down

    def project(self, points: np.ndarray, chains=None) -> np.ndarray:
        rel = np.atleast_2d(points) - self.origin
        along = rel @ self.axis
        rho = rel - np.outer(along, self.axis)
        return np.stack([along, -np.linalg.norm(rho, axis=1)], axis=-1)

    def basis(self, point) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        up = -self._outward(point)  # toward the DNA
        return self.axis, up, np.cross(self.axis, up)

    def direction(self, point, vec) -> np.ndarray:
        u, v, _ = self.basis(point)
        return np.array([float(vec @ u), float(vec @ v)])


def dna_radial_frame(sses: list[SSE], nucleic, rotate: float = 0.0, flip_v: bool = False) -> RadialFrame:
    centre, axis = _duplex_axis(nucleic)
    prot = np.array([s.centroid for s in sses]).mean(axis=0) if sses else centre + _any_perpendicular(axis)
    down = prot - centre
    down = down - (down @ axis) * axis
    down = down / np.linalg.norm(down) if np.linalg.norm(down) > 1e-6 else _any_perpendicular(axis)
    return RadialFrame(centre, axis, down)
