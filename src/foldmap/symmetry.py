"""Global symmetry of an assembly, in the spirit of the RCSB/BioJava global-symmetry search.

1. Subunits (protein chains) are clustered by aligned sequence, so copies match whatever their numbering,
   gaps or missing loops.
2. Every same-cluster pair superposition proposes an axis. For each axis the highest order n is kept for
   which a rotation by 2*pi/n maps every clustered subunit onto a copy within `tolerance` Å (CA RMSD).
3. A perpendicular two-fold as well makes the group dihedral (Dn); with no closed axis, a screw that walks
   the copies along an axis makes a helical filament.
4. Protomers are the orbits under the main axis, one copy from each orbit, grouped by largest interface.
   Chains that have no copy are left out (`others`) and drawn as usual.

Exact in cryo-EM models refined with symmetry; approximate (non-crystallographic) in X-ray models, hence
the tolerance. The user can ask for a symmetry (C2, D3, helical), give the protomers, or switch it off.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from difflib import SequenceMatcher

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .model import SSE, Backbone

TOLERANCE = 3.0  # Å, CA RMSD between a copy and its image
SAME = 0.9  # aligned identity for two chains to be copies
MIN_RES = 10  # chains shorter than this are left out
MIN_TURN = np.radians(10)  # smaller rotations are lattice-like translations, not symmetry
MIN_RISE = 2.0  # Å along the axis per copy for a helical filament
MAX_ORDER = 24


@dataclass(eq=False)
class Symmetry:
    kind: str  # "C" cyclic, "D" dihedral, "T"/"O"/"I" cubic (tetrahedral, octahedral, icosahedral), "H" helical
    n: int  # order of the main axis (C, D); number of copies along the filament (H)
    axis: np.ndarray  # unit vector of the main axis
    centre: np.ndarray  # a point on it
    protomers: list[list[str]]  # chains of each protomer, in order around (or along) the axis
    rmsd: float
    others: list[str] = field(default_factory=list)  # chains with no symmetry copy
    angles: list[float] = field(default_factory=list)  # each protomer's turn about the axis (radians, unwrapped)
    twist: float = 0.0  # helical: turn per copy (radians)
    rise: float = 0.0  # helical: Å along the axis per copy

    @property
    def label(self) -> str:
        if self.kind == "H":
            return "helical"
        return self.kind if self.kind in "TOI" else f"{self.kind}{self.n}"


# ---------------------------------------------------------------------------------------------- subunits
def _subunits(bb: Backbone) -> dict[str, tuple[str, np.ndarray]]:
    seqs: dict[str, list[str]] = {}
    xyz: dict[str, list] = {}
    for k, l in enumerate(bb.labels):
        info = gemmi.find_tabulated_residue(l.name)
        seqs.setdefault(l.chain, []).append(
            info.one_letter_code.upper() if info and info.one_letter_code != " " else "X"
        )
        xyz.setdefault(l.chain, []).append(bb.ca[k])
    return {c: ("".join(seqs[c]), np.array(xyz[c])) for c in seqs if len(seqs[c]) >= MIN_RES}


def _align(a: str, b: str) -> list[tuple[int, int]]:
    blocks = SequenceMatcher(None, a, b, autojunk=False).get_matching_blocks()
    return [(m.a + i, m.b + i) for m in blocks for i in range(m.size)]


def _clusters(subs) -> tuple[list[list[str]], dict[str, dict[int, int]]]:
    """Chains grouped as copies, and for each chain: residue index -> index in its cluster's first chain."""
    clusters: list[list[str]] = []
    to_rep: dict[str, dict[int, int]] = {}
    for c, (seq, _) in subs.items():
        for cl in clusters:
            pairs = _align(subs[cl[0]][0], seq)
            if len(pairs) >= SAME * min(len(seq), len(subs[cl[0]][0])):
                cl.append(c)
                to_rep[c] = {j: i for i, j in pairs}
                break
        else:
            clusters.append([c])
            to_rep[c] = {i: i for i in range(len(seq))}
    return clusters, to_rep


class _Assembly:
    """Clustered subunits with their CA atoms indexed by position in the cluster's first chain."""

    def __init__(self, bb: Backbone, chains: set[str] | None = None):
        subs = _subunits(bb)
        if chains is not None:
            subs = {c: v for c, v in subs.items() if c in chains}
        self.clusters, to_rep = _clusters(subs)
        self.cluster_of = {c: k for k, cl in enumerate(self.clusters) for c in cl}
        self.atoms = {c: {to_rep[c][i]: subs[c][1][i] for i in to_rep[c]} for c in subs}
        self.ca = {c: subs[c][1] for c in subs}

    def pair(self, a: str, b: str) -> tuple[np.ndarray, np.ndarray]:
        keys = sorted(self.atoms[a].keys() & self.atoms[b].keys())
        return np.array([self.atoms[a][k] for k in keys]), np.array([self.atoms[b][k] for k in keys])

    def centroid(self, c: str) -> np.ndarray:
        return self.ca[c].mean(axis=0)


# ---------------------------------------------------------------------------------------------- geometry
def _kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """R, t with q ~ p @ R.T + t, and the RMSD."""
    pc, qc = p.mean(axis=0), q.mean(axis=0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    t = qc - pc @ r.T
    return r, t, float(np.sqrt(((p @ r.T + t - q) ** 2).sum(axis=1).mean()))


def _rotation(r: np.ndarray) -> tuple[float, np.ndarray]:
    """(angle, unit axis) of a rotation matrix, the angle in [0, pi]."""
    angle = float(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1)))
    axis = np.array([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]])
    if np.linalg.norm(axis) < 1e-6:  # 0 or 180 degrees: the eigenvector of eigenvalue +1
        w, v = np.linalg.eigh((r + r.T) / 2)
        axis = v[:, np.argmax(w)]
    return angle, axis / np.linalg.norm(axis)


def _turn(axis: np.ndarray, angle: float) -> np.ndarray:
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _screw(r: np.ndarray, t: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, float]:
    """A rigid motion x -> R x + t as a screw: (angle, axis, a point on the axis, rise along it)."""
    angle, axis = _rotation(r)
    rise = float(t @ axis)
    perp = t - rise * axis
    centre = np.linalg.lstsq(np.eye(3) - r, perp, rcond=None)[0]
    centre -= (centre @ axis) * axis
    return angle, axis, centre, rise


def _candidate_axes(asm: _Assembly, tol: float) -> list[np.ndarray]:
    """Directions proposed by superposing copies onto each other, merged when within ~8 degrees."""
    axes: list[np.ndarray] = []
    for cl in asm.clusters:
        for i, a in enumerate(cl):
            for b in cl[i + 1 :]:
                p, q = asm.pair(a, b)
                if len(p) < MIN_RES:
                    continue
                r, _, rmsd = _kabsch(p, q)
                angle, ax = _rotation(r)
                if rmsd <= tol and angle >= MIN_TURN and not any(abs(ax @ x) > 0.99 for x in axes):
                    axes.append(ax)
    return axes


# ---------------------------------------------------------------------------------------------- testing
def _maps(asm: _Assembly, motion, tol: float) -> tuple[dict[str, str], float] | None:
    """Apply a motion to every clustered subunit; each must land on a distinct copy within tol Å."""
    perm, worst = {}, 0.0
    for cl in asm.clusters:
        if len(cl) < 2:
            continue
        for a in cl:
            best, score = None, np.inf
            for b in cl:
                p, q = asm.pair(a, b)
                if len(p) < MIN_RES:
                    continue
                d = float(np.sqrt(((motion(p) - q) ** 2).sum(axis=1).mean()))
                if d < score:
                    best, score = b, d
            if best is None or score > tol:
                return None
            perm[a] = best
            worst = max(worst, score)
    if len(set(perm.values())) != len(perm):
        return None
    return perm, worst


def _orbits(perm: dict[str, str]) -> list[list[str]]:
    seen, out = set(), []
    for c in sorted(perm):
        if c in seen:
            continue
        orbit = [c]
        while perm[orbit[-1]] != c:
            orbit.append(perm[orbit[-1]])
        seen |= set(orbit)
        out.append(orbit)
    return out


def _cyclic(asm: _Assembly, axis: np.ndarray, centre: np.ndarray, n: int, tol: float):
    r = _turn(axis, 2 * np.pi / n)
    found = _maps(asm, lambda x: (x - centre) @ r.T + centre, tol)
    if found is None or any(len(o) != n for o in _orbits(found[0])):
        return None
    return found


def _protomers(asm: _Assembly, perm: dict[str, str], n: int) -> list[list[str]]:
    """Protomer 0 takes one subunit from each orbit, packed together; protomer k is its k-th image."""
    orbits = sorted(_orbits(perm), key=lambda o: (asm.cluster_of[o[0]], o[0]))
    tree = {c: cKDTree(asm.ca[c]) for c in perm}
    first = [orbits[0][0]]
    for orbit in orbits[1:]:
        first.append(max(orbit, key=lambda c: (sum(tree[c].count_neighbors(tree[x], 8.0) for x in first), c)))
    out = [first]
    for _ in range(n - 1):
        out.append([perm[c] for c in out[-1]])
    return out


def _within(asm: _Assembly) -> list[str]:
    return [c for cl in asm.clusters if len(cl) >= 2 for c in cl]


def _cyclic_search(asm: _Assembly, tol: float, want: int | None = None):
    """Best (n, axis, centre, perm, rmsd) over candidate axes: highest order, then lowest RMSD."""
    sym = _within(asm)
    if len(sym) < 2:
        return None
    centre = np.mean([asm.centroid(c) for c in sym], axis=0)
    best = None
    for ax in _candidate_axes(asm, tol):
        orders = [want] if want else range(min(MAX_ORDER, len(sym)), 1, -1)
        for n in orders:
            if len(sym) % n:
                continue
            found = _cyclic(asm, ax, centre, n, tol)
            if found is not None:
                if best is None or (n, -found[1]) > (best[0], -best[4]):
                    best = (n, ax, centre, found[0], found[1])
                break
    return best


def _refine(asm: _Assembly, perm: dict[str, str], n: int, axis: np.ndarray, centre: np.ndarray, tol: float):
    """Fit the axis to all copies at once; keep the guess if the fit is no better."""
    p = np.vstack([asm.pair(a, perm[a])[0] for a in perm])
    q = np.vstack([asm.pair(a, perm[a])[1] for a in perm])
    r, t, _ = _kabsch(p, q)
    _, ax, cen, _ = _screw(r, t)
    if ax @ axis < 0:
        ax = -ax
    cen = cen + ((centre - cen) @ ax) * ax  # same height as before
    rr = _turn(ax, 2 * np.pi / n)
    found = _maps(asm, lambda x: (x - cen) @ rr.T + cen, tol)
    if found is not None and found[0] == perm:
        return ax, cen, found[1]
    return None


def _cubic(asm: _Assembly, centre: np.ndarray, tol: float) -> str | None:
    """T, O or I when the copies carry several non-parallel 3-fold axes (with a 4-fold: O; with a 5-fold: I)."""
    orders: list[tuple[np.ndarray, int]] = []
    for ax in _candidate_axes(asm, tol):
        for n in (5, 4, 3):
            if len(_within(asm)) % n == 0 and _cyclic(asm, ax, centre, n, tol) is not None:
                orders.append((ax, n))
                break
    threes = []
    for ax, n in orders:
        if n == 3 and not any(abs(ax @ t) > 0.99 for t in threes):
            threes.append(ax)
    if len(threes) < 2:
        return None
    found = {n for _, n in orders}
    return "I" if 5 in found else "O" if 4 in found else "T"


def _dihedral(asm: _Assembly, axis: np.ndarray, centre: np.ndarray, tol: float) -> bool:
    for ax in _candidate_axes(asm, tol):
        if abs(ax @ axis) < 0.15 and _cyclic(asm, ax, centre, 2, tol) is not None:
            return True
    return False


def _helical(asm: _Assembly, tol: float):
    """A screw that walks each copy onto the next, chaining every clustered subunit along one axis."""
    best = None
    for cl in asm.clusters:
        for i, a in enumerate(cl):
            for b in cl[i + 1 :]:
                p, q = asm.pair(a, b)
                if len(p) < MIN_RES:
                    continue
                r, t, rmsd = _kabsch(p, q)
                angle, ax, cen, rise = _screw(r, t)
                if rmsd > tol or abs(rise) < MIN_RISE:
                    continue
                if rise < 0:
                    continue  # the same screw is found the other way round from (b, a)
                step = {}
                for c in _within(asm):
                    for d in asm.clusters[asm.cluster_of[c]]:
                        pp, qq = asm.pair(c, d)
                        if len(pp) >= MIN_RES and np.sqrt(((pp @ r.T + t - qq) ** 2).sum(axis=1).mean()) <= tol:
                            step[c] = d
                            break
                chains = _helix_chains(step, _within(asm))
                if chains is None:
                    continue
                score = (len(chains[0]), -abs(rise))
                if best is None or score > best[0]:
                    best = (score, angle, ax, cen, rise, chains, rmsd)
    return best


def _helix_chains(step: dict[str, str], members: list[str]) -> list[list[str]] | None:
    """Follow the screw from every chain that nothing steps onto; each cluster must form one path, all equal."""
    starts = [c for c in members if c not in step.values()]
    paths = []
    for s in starts:
        path = [s]
        while path[-1] in step and len(path) <= len(members):
            path.append(step[path[-1]])
        paths.append(path)
    if not paths or sorted(c for p in paths for c in p) != sorted(members) or len({len(p) for p in paths}) != 1:
        return None
    if len(paths[0]) < 3:
        return None
    return paths


# ---------------------------------------------------------------------------------------------- public
def _parse(request: str | None) -> tuple[str, int | None]:
    if request in (None, "", "auto"):
        return "auto", None
    if request == "off":
        return "off", None
    if request in ("helical", "H"):
        return "H", None
    m = re.fullmatch(r"([CD])(\d+)", request)
    if not m or int(m.group(2)) < 2:
        raise ValueError(f"symmetry must be auto, off, helical, Cn or Dn (n >= 2); got {request!r}")
    return m.group(1), int(m.group(2))


def _around(asm: _Assembly, protomers, axis, centre) -> list[float]:
    """Each protomer's turn about the axis relative to protomer 0 (radians, increasing in order)."""
    ref = np.mean([asm.centroid(c) for c in protomers[0]], axis=0) - centre
    ref -= (ref @ axis) * axis
    side = np.cross(axis, ref)
    out, last = [], 0.0
    for p in protomers:
        rel = np.mean([asm.centroid(c) for c in p], axis=0) - centre
        a = float(np.arctan2(rel @ side, rel @ ref))
        while out and a < last - 1e-6:
            a += 2 * np.pi
        out.append(a)
        last = a
    return out


def _from_user(asm: _Assembly, protomers: list[list[str]], tol: float) -> Symmetry:
    known = set(asm.ca)
    missing = sorted({c for p in protomers for c in p} - known)
    if missing:
        raise ValueError(f"protomers name chains that are not protein chains here: {', '.join(missing)}")
    if len(protomers) < 2 or len({len(p) for p in protomers}) != 1:
        raise ValueError("give at least two protomers with the same number of chains, e.g. A,B;C,D")

    p0 = np.vstack([asm.pair(protomers[0][i], protomers[1][i])[0] for i in range(len(protomers[0]))])
    p1 = np.vstack([asm.pair(protomers[0][i], protomers[1][i])[1] for i in range(len(protomers[0]))])
    r, t, rmsd = _kabsch(p0, p1)
    angle, axis, centre, rise = _screw(r, t)
    n = len(protomers)
    others = sorted(known - {c for p in protomers for c in p})
    if abs(rise) >= MIN_RISE and abs(angle - 2 * np.pi / n) > np.radians(15):
        if rise < 0:
            axis, rise = -axis, -rise
        sym = Symmetry("H", n, axis, centre, protomers, rmsd, others, twist=angle, rise=rise)
    else:
        centre = np.mean([asm.centroid(c) for p in protomers for c in p], axis=0)
        sym = Symmetry("C", n, axis, centre, protomers, rmsd, others)
    sym.angles = _around(asm, protomers, sym.axis, sym.centre) if sym.kind != "H" else [k * angle for k in range(n)]
    return sym


def detect_symmetry(
    bb: Backbone,
    sses: list[SSE] | None = None,
    request: str | None = "auto",
    protomers: list[list[str]] | None = None,
    tolerance: float = TOLERANCE,
) -> Symmetry | None:
    """The assembly's symmetry, or None. `request`: auto, off, Cn, Dn or helical (an error if absent);
    `protomers`: the user's own grouping, in order; `tolerance`: CA RMSD in Å between copies."""
    kind, n = _parse(request)
    if kind == "off":
        return None
    asm = _Assembly(bb)
    if protomers:
        return _from_user(asm, protomers, tolerance)
    found = None
    if kind in ("auto", "C", "D"):
        best = _cyclic_search(asm, tolerance, n)
        if best is not None:
            order, axis, centre, perm, rmsd = best
            refined = _refine(asm, perm, order, axis, centre, tolerance)
            if refined is not None:
                axis, centre, rmsd = refined
            dihedral = _dihedral(asm, axis, centre, tolerance)
            cubic = _cubic(asm, centre, tolerance) if dihedral or order >= 3 else None
            if kind != "D" or dihedral:
                groups = _protomers(asm, perm, order)
                others = sorted(set(asm.ca) - set(perm))
                label = cubic if cubic and kind == "auto" else "D" if dihedral else "C"
                found = Symmetry(label, order, axis, centre, groups, rmsd, others)
                found.angles = [2 * np.pi * k / order for k in range(order)]
                if kind == "C" and found.kind == "D":
                    found.kind = "C"  # asked for the cyclic part only
    if found is None and kind in ("auto", "H"):
        best = _helical(asm, tolerance)
        if best is not None:
            _, angle, axis, centre, rise, paths, rmsd = best
            groups = [list(g) for g in zip(*paths)]
            others = sorted(set(asm.ca) - {c for g in groups for c in g})
            found = Symmetry("H", len(groups), axis, centre, groups, rmsd, others, twist=angle, rise=rise)
            found.angles = [k * angle for k in range(len(groups))]
    if found is None and kind != "auto":
        best = detect_symmetry(bb, sses, "auto", tolerance=tolerance) if kind != "auto" else None
        raise ValueError(
            f"no {kind}{n or ''} symmetry within {tolerance:g} Å (found: {best.label if best else 'none'})".replace(
                "Hhelical", "helical"
            )
        )
    return found
