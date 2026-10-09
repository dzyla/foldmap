"""Where the lipid bilayer is: exactly from OPM's dummy atoms when the file has them, otherwise estimated from
the structure (the long helices and strands give the normal; a 30 Å slab is slid to the most hydrophobic,
least charged position). The cytoplasmic side follows the positive-inside rule (more Lys/Arg)."""

from __future__ import annotations

from dataclasses import dataclass

import gemmi
import numpy as np

from .dssp import assign_dssp
from .model import Backbone
from .ss import build_sses

HALF = 15.0  # Å, half the hydrophobic thickness of a typical bilayer
_HYDROPHOBIC = {"ALA", "ILE", "LEU", "MET", "PHE", "VAL", "TRP", "CYS"}
_POLAR = {"LYS", "ARG", "ASP", "GLU", "ASN", "GLN", "HIS"}
_MIN_HYDROPHOBIC, _MAX_POLAR = 0.55, 0.15  # a slab this oily and this uncharged is membrane
_SIDED_RATIO, _SIDED_DIFF = 2.0, 3  # Lys/Arg counts must differ this clearly to call the sides
_MIN_CROSSING = 10.0  # Å an element must reach beyond the centre on both sides to cross the membrane


@dataclass
class Membrane:
    normal: np.ndarray  # unit vector pointing to the outside (non-cytoplasmic side)
    centre: float  # position of the mid-plane along the normal (Å)
    half: float  # half thickness (Å)
    source: str  # 'opm' (dummy atoms in the file) or 'estimate'
    sided: bool = False  # True when the positive-inside rule clearly says which side is the cytoplasm

    def depth(self, xyz: np.ndarray) -> np.ndarray:
        """Signed distance from the mid-plane, positive outside."""
        return np.asarray(xyz) @ self.normal - self.centre


def _opm(path) -> Membrane | None:
    st = gemmi.read_structure(str(path))
    dummies = [a.pos for ch in st[0] for r in ch if r.name == "DUM" for a in r]
    if not dummies:
        return None
    z = np.array([p.z for p in dummies])
    return Membrane(np.array([0.0, 0.0, 1.0]), 0.0, float(np.abs(z).max()), "opm")


def find_membrane(path, bb: Backbone) -> Membrane | None:
    """The membrane, or None for a soluble protein."""
    found = _opm(path)
    if found is None:
        found = _estimate(bb)
    if found is None:
        return None
    return _orient(found, bb)


def _estimate(bb: Backbone) -> Membrane | None:
    sses = build_sses(bb, assign_dssp(bb).ss)
    long = [s for s in sses if (s.kind == "H" and len(s) >= 16) or (s.kind == "E" and len(s) >= 8)]
    if not long:
        return None
    normal = np.linalg.eigh(sum(len(s) * np.outer(s.axis, s.axis) for s in long))[1][:, -1]
    z = bb.ca @ normal
    hyd = np.array([l.name in _HYDROPHOBIC for l in bb.labels])
    pol = np.array([l.name in _POLAR for l in bb.labels])
    best = None
    for c in np.arange(z.min(), z.max() + 1.0, 1.0):
        inside = np.abs(z - c) < HALF
        if inside.sum() < 20:
            continue
        score = hyd[inside].sum() - 1.5 * pol[inside].sum()
        if best is None or score > best[0]:
            best = (score, c, inside)
    if best is None:
        return None
    _, centre, inside = best
    crossing = [
        s
        for s in long
        if z[s.start : s.end + 1].min() < centre - _MIN_CROSSING
        and z[s.start : s.end + 1].max() > centre + _MIN_CROSSING
    ]
    if hyd[inside].mean() < _MIN_HYDROPHOBIC or pol[inside].mean() > _MAX_POLAR or not crossing:
        return None
    return Membrane(normal, float(centre), HALF, "estimate")


def _orient(m: Membrane, bb: Backbone) -> Membrane:
    """Point the normal outward: the side with fewer Lys/Arg just beyond the membrane (positive-inside rule)."""
    d = m.depth(bb.ca)
    basic = np.array([l.name in ("LYS", "ARG") for l in bb.labels])
    band = lambda lo, hi: int((basic & (d > lo) & (d < hi)).sum())  # noqa: E731
    up, down = band(m.half - 5, m.half + 15), band(-m.half - 15, -m.half + 5)
    sided = max(up, down) >= _SIDED_RATIO * min(up, down) and abs(up - down) >= _SIDED_DIFF
    if up > down:
        return Membrane(-m.normal, -m.centre, m.half, m.source, sided)
    return Membrane(m.normal, m.centre, m.half, m.source, sided)
