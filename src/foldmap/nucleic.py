# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Nucleic-acid secondary structure from coordinates: base pairs, stems, single-stranded runs and bends.

Strands are model.NAStrand (labels, P and C1' coordinates); pairing reads the gemmi residues behind them."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import NamedTuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

HBOND_MAX = 3.8  # Å, donor-acceptor distance of a base-pair hydrogen bond (stretched pairs at 2-3 Å resolution)
C1_PAIR = (9.0, 12.0)  # Å, C1'-C1' across a Watson-Crick or wobble pair (10.5 +- 1.5)
PLANE_MAX = 45.0  # degrees between the two base planes (TBP and nucleosome DNA reach 40 of propeller twist)
MIN_LOOP = 3  # a hairpin closes at least this many unpaired nucleotides

# (purine, pyrimidine) parent bases -> the hydrogen bonds checked, how many must hold, the pair kind
_PAIRS = {
    ("G", "C"): ((("N1", "N3"), ("O6", "N4"), ("N2", "O2")), 2, "WC"),
    ("I", "C"): ((("N1", "N3"), ("O6", "N4")), 2, "WC"),
    ("A", "U"): ((("N1", "N3"), ("N6", "O4")), 2, "WC"),
    ("A", "T"): ((("N1", "N3"), ("N6", "O4")), 2, "WC"),
    ("G", "U"): ((("N1", "O2"), ("O6", "N3")), 2, "GU"),
    ("G", "T"): ((("N1", "O2"), ("O6", "N3")), 2, "GU"),
}
_PURINE = {"A", "G", "I"}
_RING = ("N1", "C2", "N3", "C4", "C5", "C6")


class BasePair(NamedTuple):
    """Nucleotide i of strand a paired with nucleotide j of strand b; (a, i) < (b, j)."""

    a: int
    i: int
    b: int
    j: int
    kind: str  # "WC" or "GU"


def parent_base(res) -> str | None:
    """A, C, G, U, T or I: the base a (possibly modified) nucleotide pairs as; None if it has no base."""
    info = gemmi.find_tabulated_residue(res.name)
    code = info.one_letter_code.upper() if info is not None and info.is_nucleic_acid() else " "
    if code in ("A", "C", "G", "U", "T", "I"):
        return code
    names = {a.name for a in res}
    if "N9" in names:  # a purine
        return "G" if "O6" in names else "A" if "N6" in names else None
    if "N4" in names:
        return "C"
    if "O4" in names and "N3" in names:  # uracil or thymine (pseudouridine has the same atom names)
        return "T" if "C7" in names or "C5M" in names else "U"
    return None


def _glycosidic(res) -> np.ndarray | None:
    """C1' -> the base atom bonded to it (N9 of a purine, N1 of a pyrimidine, C5 of pseudouridine)."""
    c1 = res.find_atom("C1'", "*")
    for name in ("N9", "N1", "C5"):
        a = res.find_atom(name, "*")
        if c1 is not None and a is not None and c1.pos.dist(a.pos) < 1.7:
            return np.array(a.pos.tolist()) - np.array(c1.pos.tolist())
    return None


def _cis(x: np.ndarray, gx: np.ndarray | None, y: np.ndarray, gy: np.ndarray | None) -> bool:
    """Are both glycosidic bonds on the same side of the C1'-C1' line (a cis pair, as in a double helix)?
    Trans pairs such as the tRNA Levitt pair G15·C48 use the same atoms but put them on opposite sides."""
    if gx is None or gy is None:
        return True
    u = y - x
    n, m = np.cross(gx, u), np.cross(u, gy)
    return float(n @ m) < -0.5 * float(np.linalg.norm(n) * np.linalg.norm(m))


def _normal(res) -> np.ndarray | None:
    pts = [a.pos.tolist() for a in res if a.name in _RING]
    if len(pts) < 4:
        return None
    pts = np.array(pts)
    return np.linalg.svd(pts - pts.mean(axis=0))[2][2]


def pair_nucleotides(runs: list[list]) -> list[BasePair]:
    """Watson-Crick and G·U wobble pairs among the gemmi residues of each strand (runs[strand][nucleotide]).

    A pair needs its hydrogen bonds (two of the three for G-C), the C1'-C1' span of a paired duplex, roughly
    coplanar bases in the cis orientation of a double helix and, on one strand, more than MIN_LOOP nucleotides
    between the partners. Each nucleotide
    keeps its best pair (more bonds, then shorter ones); only mutual best pairs are kept."""
    bases = []  # (strand, index, parent, residue, C1' position)
    for s, run in enumerate(runs):
        for i, r in enumerate(run):
            c1 = r.find_atom("C1'", "*")
            base = parent_base(r)
            if c1 is not None and base is not None:
                bases.append((s, i, base, r, np.array(c1.pos.tolist())))
    if len(bases) < 2:
        return []
    tree = cKDTree(np.array([b[4] for b in bases]))
    normals: dict[int, np.ndarray | None] = {}
    best: dict[tuple[int, int], tuple[tuple[int, float], tuple[int, int], str]] = {}
    for x, y in sorted(tree.query_pairs(C1_PAIR[1])):
        bx, by = bases[x], bases[y]
        if bx[0] == by[0] and abs(bx[1] - by[1]) <= MIN_LOOP:
            continue
        if float(np.linalg.norm(bx[4] - by[4])) < C1_PAIR[0]:
            continue
        pur, pyr = (bx, by) if bx[2] in _PURINE else (by, bx)
        rule = _PAIRS.get((pur[2], pyr[2]))
        if rule is None:
            continue
        bonds, need, kind = rule
        found = []
        for p_atom, q_atom in bonds:
            p, q = pur[3].find_atom(p_atom, "*"), pyr[3].find_atom(q_atom, "*")
            if p is not None and q is not None and p.pos.dist(q.pos) <= HBOND_MAX:
                found.append(p.pos.dist(q.pos))
        if len(found) < need:
            continue
        for k in (x, y):
            if k not in normals:
                normals[k] = _normal(bases[k][3])
        nx, ny = normals[x], normals[y]
        if nx is not None and ny is not None and abs(float(nx @ ny)) < np.cos(np.radians(PLANE_MAX)):
            continue
        if not _cis(bx[4], _glycosidic(bx[3]), by[4], _glycosidic(by[3])):
            continue
        score = (len(found), -float(np.mean(found)))
        for me, other in ((bx, by), (by, bx)):
            key = (me[0], me[1])
            if key not in best or score > best[key][0]:
                best[key] = (score, (other[0], other[1]), kind)
    out = set()
    for key, (_, other, kind) in best.items():
        if best.get(other, (None, None))[1] == key:
            (a, i), (b, j) = sorted((key, other))
            out.add(BasePair(a, i, b, j, kind))
    return sorted(out)


@dataclass(eq=False)
class NAStructure:  # filled in by Task 3
    stems: list = field(default_factory=list)
    runs: list = field(default_factory=list)
    ties: list = field(default_factory=list)


class NAGap(NamedTuple):  # Task 3 gives it its docstring
    before: int
    after: int
    count: int


def secondary_structure(strands, pairs):  # replaced in Task 3
    return NAStructure()
