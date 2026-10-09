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


def symmetric_frame(sses: list[SSE], sym) -> CylinderFrame:
    cents = np.array([s.centroid for s in sses])
    rel = cents - sym.centre
    rho = rel - np.outer(rel @ sym.axis, sym.axis)
    radius = float(np.median(np.linalg.norm(rho, axis=1)))
    first = np.array([s.centroid for s in sses if s.chain in sym.protomers[0]]).mean(axis=0) - sym.centre
    ref = first - (first @ sym.axis) * sym.axis
    ref = ref / np.linalg.norm(ref) if np.linalg.norm(ref) > 1e-6 else _any_perpendicular(sym.axis)
    sector = {c: k for k, p in enumerate(sym.protomers) for c in p}
    classes = {c: i for p in sym.protomers for i, c in enumerate(p)}
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


def dna_view_frame(sses: list[SSE], nucleic, rotate: float = 0.0, flip_v: bool = False) -> Frame:
    """The usual picture of a protein-DNA complex: the longest duplex runs across the page (5' of its first
    strand on the left) and the protein sits below it. A rotation, never a mirror (unless flip_v)."""
    from .model import duplex_axis, duplex_groups

    biggest = max(duplex_groups(nucleic), key=lambda g: sum(len(nucleic.strands[k]) for k in g))
    dna_c, u = duplex_axis(nucleic, biggest)  # the same axis the layout draws the duplex along
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
    from .model import duplex_axis, duplex_groups

    biggest = max(duplex_groups(nucleic), key=lambda g: sum(len(nucleic.strands[k]) for k in g))
    centre, axis = duplex_axis(nucleic, biggest)
    prot = np.array([s.centroid for s in sses]).mean(axis=0) if sses else centre + _any_perpendicular(axis)
    down = prot - centre
    down = down - (down @ axis) * axis
    down = down / np.linalg.norm(down) if np.linalg.norm(down) > 1e-6 else _any_perpendicular(axis)
    return RadialFrame(centre, axis, down)
