# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Non-covalent interactions other than hydrogen bonds, from the heavy atoms of a model: one finder per type.

Partners are chosen as for hydrogen bonds (`hbonds.parse_selection`: off, na, interchain, chains:A+B,C+D), residues
are read as `hbonds` reads them (polymer residues only, so waters and ligands are left out; of residues sharing a
number the first; of each atom its first conformer), and two atoms of one residue never interact.

Salt bridge: an anionic O within 4.0 Å of a cationic N, atom to atom (Barlow & Thornton, J. Mol. Biol. 168:867,
1983; PLIP, Salentin et al. 2015 / Adasme et al. 2021, measures 5.5 Å between charge centres instead).

- anionic O: Asp OD1/OD2, Glu OE1/OE2, the O and OXT of a chain's last protein residue (or of any residue with an
  OXT), nucleic-acid phosphate OP1/OP2;
- cationic N: Lys NZ, Arg NE/NH1/NH2, the N of a chain's first protein residue; His ND1/NE2 with the note
  "His if protonated" (its charge depends on pH and surroundings).

A chain's first and last protein residues are its termini; a gap inside the chain (missing residues) is not, since
the chain goes on. Every cation-anion atom pair within the limit is one Interaction, so an Arg-Asp pair may give
several.

Cation-π: Lys NZ, or the Arg guanidinium centre (the mean of NE, NH1 and NH2, as PLIP takes the centre of the
charged N atoms), at most 6.0 Å from an aromatic ring's centroid, with an offset of at most 2.0 Å: the distance
from the centroid to the cation's projection on the ring plane (PLIP's PICATION_DIST_MAX and PISTACK_OFFSET_MAX,
Salentin et al., Nucleic Acids Res. 43:W443, 2015; Adasme et al., Nucleic Acids Res. 49:W530, 2021; the 6.0 Å
from Gallivan & Dougherty, PNAS 96:9459, 1999). The rings: Phe, Tyr, Trp (its five- and six-membered rings), His,
and the nucleobases (the six-membered ring of every base, the five-membered ring of purines); a ring counts only
when all its atoms are in the model. The ring plane is the least-squares plane through its atoms. The angle
reported is between the ring normal and the centroid-cation line (0° straight above the ring).

Each Interaction names its partners as (chain, residue tag, atom): the residue tag as `ResLabel.tag` (number and
insertion code), the atom its name, "centre" for the Arg guanidinium centre, "ring:6" or "ring:5" for a ring (its
size; one ring of each size per residue). `a` is the cationic partner. Output is sorted by (a, b, kind); distances
and angles are rounded to 1e-3.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from scipy.spatial import cKDTree

from .hbonds import NUCLEIC, PROTEIN, _heavy, _kind, parse_selection
from .model import ResLabel

KINDS = ("hbond", "salt", "stack", "cation-pi", "hydrophobic")

SALT_MAX = 4.0  # Å, anionic O to cationic N
CATION_PI_MAX = 6.0  # Å, cation to ring centroid
CATION_PI_OFFSET = 2.0  # Å, ring centroid to the cation's projection on the ring plane
HIS_NOTE = "His if protonated"

_ANION_SIDE = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
_CATION_SIDE = {"LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2"), "HIS": ("ND1", "NE2")}
_PHOSPHATE = ("OP1", "OP2")
_C_TERMINAL = ("O", "OXT")
_GUANIDINIUM = ("NE", "NH1", "NH2")

_SIX = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
_PROTEIN_RINGS = {
    "PHE": (_SIX,),
    "TYR": (_SIX,),
    "TRP": (("CG", "CD1", "NE1", "CE2", "CD2"), ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")),
    "HIS": (("CG", "ND1", "CD2", "CE1", "NE2"),),
}
_BASE_RINGS = (("N1", "C2", "N3", "C4", "C5", "C6"), ("C4", "C5", "N7", "C8", "N9"))


@dataclass(frozen=True)
class Interaction:
    """One interaction between two partners, each (chain, residue tag, atom or "ring:K" / "centre").
    kind: one of KINDS. distance in Å; angle in degrees (stacking; cation-π: ring normal to the cation), else
    None; note: a caveat to show with it ("His if protonated")."""

    kind: str
    a: tuple[str, str, str]
    b: tuple[str, str, str]
    distance: float
    angle: float | None
    note: str = ""

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"interactions: no kind {self.kind!r}; use one of {', '.join(KINDS)}")


class _Residue(NamedTuple):
    label: ResLabel
    chain: int  # index into the chain names
    kind: int  # PROTEIN or NUCLEIC
    atoms: dict  # heavy atom name -> gemmi.Position
    first: bool  # the chain's first protein residue
    last: bool  # the chain's last protein residue


def _residues(model, select) -> tuple[object, list[_Residue], list[str]]:
    """The selection `select` names, the polymer residues of the chains it can involve, and those chains' names."""
    if isinstance(model, tuple):
        model = model[0]
    sel = parse_selection(select)
    if sel.kind == "off":
        return sel, [], []
    wanted = None
    if sel.kind == "chains":
        present = {c.name for c in model}
        missing = sorted({n for p in sel.pairs for n in p} - present)
        if missing:
            raise ValueError(
                f"interactions: no chain {', '.join(missing)} in the model; it has {', '.join(sorted(present))}"
            )
        wanted = {n for p in sel.pairs for n in p}
    out: list[_Residue] = []
    names: list[str] = []
    for chain in model:
        if wanted is not None and chain.name not in wanted:
            continue
        if chain.name not in names:
            names.append(chain.name)
        c = names.index(chain.name)
        kept, seen = [], set()
        for res in chain:  # one residue per number, the first (microheterogeneity), as hbonds does
            kind = _kind(res)
            key = (res.seqid.num, res.seqid.icode.strip())
            if kind is not None and key not in seen:
                seen.add(key)
                kept.append((res, kind))
        protein = [n for n, (_, kind) in enumerate(kept) if kind == PROTEIN]
        for n, (res, kind) in enumerate(kept):
            label = ResLabel(chain.name, res.seqid.num, res.seqid.icode.strip(), res.name)
            first = bool(protein) and n == protein[0]
            last = bool(protein) and n == protein[-1]
            out.append(_Residue(label, c, kind, _heavy(res), first, last))
    return sel, out, names


def _allowed(sel, names: list[str], residues: list[_Residue], i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Which residue pairs (index arrays into residues) the selection lets interact; never a residue with itself."""
    if not len(i):
        return np.zeros(0, bool)
    chain = np.array([r.chain for r in residues], int)
    kind = np.array([r.kind for r in residues], int)
    if sel.kind == "na":
        keep = kind[i] != kind[j]
    elif sel.kind == "interchain":
        keep = chain[i] != chain[j]
    else:
        allowed = np.zeros((len(names),) * 2, bool)
        for a, b in sel.pairs:
            allowed[names.index(a), names.index(b)] = True
            allowed[names.index(b), names.index(a)] = True
        keep = allowed[chain[i], chain[j]]
    return keep & (i != j)


def _key(res: _Residue, atom: str) -> tuple[str, str, str]:
    return (res.label.chain, res.label.tag, atom)


def _sorted(found: list[Interaction]) -> list[Interaction]:
    return sorted(found, key=lambda x: (x.a, x.b, x.kind))


def _charged(residues: list[_Residue]) -> list[tuple[list[float], int, int, str]]:
    """The charged atoms: (xyz, +1 cation / -1 anion, residue index, atom name), in model order."""
    rows = []
    for r, res in enumerate(residues):
        atoms, name = res.atoms, res.label.name
        if res.kind == PROTEIN:
            cations = list(_CATION_SIDE.get(name, ()))
            anions = list(_ANION_SIDE.get(name, ()))
            if res.first:
                cations.insert(0, "N")
            if res.last or "OXT" in atoms:
                anions.extend(_C_TERMINAL)
        else:
            cations, anions = [], list(_PHOSPHATE)
        for sign, wanted in ((1, cations), (-1, anions)):
            for atom in wanted:
                if atom in atoms:
                    rows.append((atoms[atom].tolist(), sign, r, atom))
    return rows


def find_salt_bridges(model, select) -> list[Interaction]:
    """The salt bridges of `model` (a gemmi Model, or what io.read_model returns) between the partners `select`
    names (a hbonds Selection or its text): every cationic N - anionic O atom pair within SALT_MAX, the cation as
    `a`, sorted by (a, b, kind)."""
    sel, residues, names = _residues(model, select)
    rows = _charged(residues)
    if len(rows) < 2:
        return []
    xyz = np.array([x[0] for x in rows], float)
    sign = np.array([x[1] for x in rows], int)
    owner = np.array([x[2] for x in rows], int)
    ij = cKDTree(xyz).query_pairs(SALT_MAX, output_type="ndarray")
    if not len(ij):
        return []
    i, j = ij[:, 0], ij[:, 1]
    keep = (sign[i] != sign[j]) & _allowed(sel, names, residues, owner[i], owner[j])
    i, j = i[keep], j[keep]
    cation, anion = np.where(sign[i] > 0, i, j), np.where(sign[i] > 0, j, i)
    dist = np.linalg.norm(xyz[cation] - xyz[anion], axis=1)
    out = []
    for c, a, d in zip(cation, anion, dist):
        if d > SALT_MAX:
            continue
        res = residues[owner[c]]
        note = HIS_NOTE if res.label.name == "HIS" and rows[c][3] != "N" else ""
        out.append(
            Interaction(
                "salt", _key(res, rows[c][3]), _key(residues[owner[a]], rows[a][3]), round(float(d), 3), None, note
            )
        )
    return _sorted(out)


def _rings(residues: list[_Residue]) -> list[tuple[np.ndarray, np.ndarray, int, str]]:
    """The aromatic rings with all their atoms in the model: (centroid, unit normal, residue index, "ring:K")."""
    out = []
    for r, res in enumerate(residues):
        if res.kind == NUCLEIC:
            rings = _BASE_RINGS
        else:
            rings = _PROTEIN_RINGS.get(res.label.name, ())
        for ring in rings:
            if not all(a in res.atoms for a in ring):
                continue
            xyz = np.array([res.atoms[a].tolist() for a in ring], float)
            centre = xyz.mean(axis=0)
            normal = np.linalg.svd(xyz - centre)[2][2]
            out.append((centre, normal / np.linalg.norm(normal), r, f"ring:{len(ring)}"))
    return out


def _cations(residues: list[_Residue]) -> list[tuple[np.ndarray, int, str]]:
    """Lys NZ and the Arg guanidinium centre: (xyz, residue index, atom name or "centre")."""
    out = []
    for r, res in enumerate(residues):
        if res.kind != PROTEIN:
            continue
        atoms = res.atoms
        if res.label.name == "LYS" and "NZ" in atoms:
            out.append((np.array(atoms["NZ"].tolist(), float), r, "NZ"))
        elif res.label.name == "ARG" and all(a in atoms for a in _GUANIDINIUM):
            out.append((np.mean([atoms[a].tolist() for a in _GUANIDINIUM], axis=0), r, "centre"))
    return out


def find_cation_pi(model, select) -> list[Interaction]:
    """The cation-π interactions of `model` between the partners `select` names: Lys NZ or the Arg guanidinium
    centre (`a`) within CATION_PI_MAX of a ring centroid (`b`) and at most CATION_PI_OFFSET off the ring's axis;
    angle: between the ring normal and the centroid-cation line. Sorted by (a, b, kind)."""
    sel, residues, names = _residues(model, select)
    rings, cations = _rings(residues), _cations(residues)
    if not rings or not cations:
        return []
    centres = np.array([x[0] for x in rings], float)
    tree = cKDTree(centres)
    near = tree.query_ball_point(np.array([x[0] for x in cations], float), CATION_PI_MAX, return_sorted=True)
    pairs = np.array([(c, k) for c, ks in enumerate(near) for k in ks], int).reshape(-1, 2)
    if not len(pairs):
        return []
    owner_c = np.array([cations[c][1] for c in pairs[:, 0]], int)
    owner_r = np.array([rings[k][2] for k in pairs[:, 1]], int)
    pairs = pairs[_allowed(sel, names, residues, owner_c, owner_r)]
    out = []
    for c, k in pairs:
        xyz, r, atom = cations[c]
        centre, normal, q, ring = rings[k]
        v = xyz - centre
        d = float(np.linalg.norm(v))
        height = float(np.dot(v, normal))
        offset = float(np.linalg.norm(v - height * normal))
        if d == 0.0 or d > CATION_PI_MAX or offset > CATION_PI_OFFSET:
            continue
        angle = float(np.degrees(np.arccos(min(1.0, abs(height) / d))))
        out.append(
            Interaction("cation-pi", _key(residues[r], atom), _key(residues[q], ring), round(d, 3), round(angle, 3))
        )
    return _sorted(out)
