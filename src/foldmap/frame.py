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
    outside every protomer are projected as they lie. Locally a proper rotation (never a mirror)."""

    origin: np.ndarray  # a point on the axis: the page centre
    axis: np.ndarray  # toward the viewer (w)
    ref: np.ndarray  # page up: perpendicular to the axis, toward protomer 0
    n: int
    sector: dict  # chain -> protomer index
    classes: dict = None  # chain -> chain class (its position within a protomer)
    angles: list = None  # each protomer's exact turn about the axis, radians (k * 2 pi / n)

    def __post_init__(self):
        if self.angles is None:
            self.angles = [2 * np.pi * k / self.n for k in range(self.n)]

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


def frame_like(frame_a: Frame, rotation: np.ndarray, translation: np.ndarray) -> Frame:
    """Structure B's frame drawn like A's (plan 13a, --like): A's frame carried through the superposition A ≈ R·B + t
    (align.ChainAlignment's rotation and translation), so a point b of B lands on the page where its place R·b + t in
    A lands. B's page axes are Rᵀu, Rᵀv, Rᵀw and its origin Rᵀ(origin_A − t); A's --up, --view, --rotate, --flip-v and
    membrane turn come with it. R is a proper rotation, so B is seen as A is: in a mirror only where A's frame is one
    (--flip-v negates v alone, so A's axes are then improper and B's are too, the two drawn alike). A plain
    Frame only: a frame drawing copies as mirror images, an unrolled cylinder or a rosette places each protomer its own
    way, and cannot be carried by one rotation (TypeError)."""
    if type(frame_a) is not Frame or frame_a.mirror_x is not None or frame_a.sector:
        raise TypeError(f"frame_like carries a plain Frame only, not {type(frame_a).__name__} with copies")
    r = np.asarray(rotation, float)
    t = np.asarray(translation, float)
    return Frame(r.T @ (np.asarray(frame_a.origin, float) - t), r.T @ frame_a.u, r.T @ frame_a.v, r.T @ frame_a.w)


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
