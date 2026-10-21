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

π-stacking: two of those rings with centroids at most 5.5 Å apart, parallel (the angle between their normals at
most 30°) or T-shaped (60-90°), with an offset of at most 2.0 Å: the smaller of the two distances from one centroid
to the other's projection on its ring plane (PLIP's PISTACK_DIST_MAX, PISTACK_ANG_DEV and PISTACK_OFFSET_MAX). The
angle reported is the one between the normals; the note says "parallel" or "T-shaped". Two nucleobases never stack
here, whatever the selection: base-base stacking within nucleic acid is the ladder's to show. Every ring pair is one
Interaction, so a Trp or a purine may stack twice with one ring.

Hydrophobic contact: two carbons, each bonded only to carbon, at most 4.0 Å apart (PLIP's HYDROPH_DIST_MAX and its
hydrophobic atom, a carbon with only carbon or hydrogen neighbours). Bonds are read from the coordinates: heavy
atoms of one residue at most 2.0 Å apart (C-S is 1.8 Å, C-Se 1.95 Å; atoms two bonds apart sit 2.2 Å or more). That
leaves the side-chain carbons (CA is bonded to N, C to O, and Ser CB, Thr CB, Cys CB, Met CG and CE, Pro CD to a
heteroatom), and in nucleic acid the DNA C2' and the base carbons with no heteroatom neighbour (thymine C5 and C7,
the C5 of U and C). One Interaction per residue pair, the shortest (ties to the first atoms by name).

Each Interaction names its partners as (chain, residue tag, atom): the residue tag as `ResLabel.tag` (number and
insertion code), the atom its name, "centre" for the Arg guanidinium centre, "ring:6" or "ring:5" for a ring (its
size; one ring of each size per residue). `a` is the cationic partner (salt bridge, cation-π); of two partners of
like role (stacking, hydrophobic) the one that sorts first. Output is sorted by `sort_key`: a, then b, each by chain,
residue number, insertion code and atom, then kind.

`verify_interactions` re-checks drawn lines from the coordinates by a code path of its own (residues, atoms, rings,
planes and distances read and measured again, with no finder's helpers): the lines that break their rule.

Every measured quantity (distance, offset, angle) is snapped to a 1e-3 grid (Å or degrees) before it meets a limit,
and the limits sit on that grid, so a contact at a limit is judged the same on every CPU: the last bits a BLAS kernel
or a SIMD sum leaves in a centroid or a ring normal never decide it (A12/A13). The reported values are the snapped
ones, so a value shown at the limit always passed it.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import NamedTuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .hbonds import MAX_DIST, MIN_DIST, NUCLEIC, PROTEIN, _heavy, _kind, parse_selection
from .model import ResLabel

KINDS = ("hbond", "salt", "stack", "cation-pi", "hydrophobic")

SALT_MAX = 4.0  # Å, anionic O to cationic N
CATION_PI_MAX = 6.0  # Å, cation to ring centroid
CATION_PI_OFFSET = 2.0  # Å, ring centroid to the cation's projection on the ring plane
HIS_NOTE = "His if protonated"
STACK_MAX = 5.5  # Å, ring centroid to ring centroid
STACK_PARALLEL = 30.0  # degrees between the ring normals, at most, for a parallel stack
STACK_T_MIN = 60.0  # ... at least, for a T-shaped one (up to 90)
STACK_OFFSET = 2.0  # Å, centroid to the other centroid's projection on its ring plane (the smaller of the two)
PARALLEL, T_SHAPED = "parallel", "T-shaped"
HYDROPHOBIC_MAX = 4.0  # Å, carbon to carbon
_COVALENT = 2.0  # Å, two heavy atoms of one residue this close are bonded
_GRID = 1e3  # measured values are snapped to 1/_GRID (Å, degrees) before any limit sees them
_REACH = 0.5 / _GRID  # a tree search reaches this far past a limit, so a value that snaps onto it is not lost

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
    elements: dict  # heavy atom name -> element name ("C", "N", ...)
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
            out.append(_Residue(label, c, kind, _heavy(res), _elements(res), first, last))
    return sel, out, names


def _elements(res) -> dict[str, str]:
    """The element of each heavy atom of a residue, by name, as `_heavy` reads the atoms."""
    out: dict[str, str] = {}
    for atom in res:
        if atom.element.name not in ("H", "D"):
            out.setdefault(atom.name, atom.element.name)
    return out


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


def _snap(v: float) -> float:
    """v on the 1/_GRID grid; never -0.0 (as layout._snap)."""
    return float(np.round(float(v) * _GRID) / _GRID) + 0.0


def _snaps(v: np.ndarray) -> np.ndarray:
    """_snap of each value."""
    return np.round(v * _GRID) / _GRID + 0.0


_TAG = re.compile(r"(-?\d+)(.*)")


def _partner_key(p: tuple[str, str, str]) -> tuple[str, int, str, str]:
    """(chain, residue number, insertion code, atom): residue numbers compare as numbers, so 5 comes before 19."""
    chain, tag, atom = p
    m = _TAG.fullmatch(tag)
    return (chain, int(m.group(1)), m.group(2), atom) if m else (chain, 0, tag, atom)


def _ordered(p: tuple[str, str, str], q: tuple[str, str, str]) -> tuple[tuple[str, str, str], tuple[str, str, str]]:
    """Two partners of like role, the one that sorts first first."""
    return (p, q) if _partner_key(p) <= _partner_key(q) else (q, p)


def sort_key(x: Interaction) -> tuple:
    """The order of found interactions: a, then b (each by chain, residue number, insertion code, atom), then kind."""
    return (_partner_key(x.a), _partner_key(x.b), x.kind)


def _sorted(found: list[Interaction]) -> list[Interaction]:
    return sorted(found, key=sort_key)


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
    `a`, sorted by `sort_key`."""
    sel, residues, names = _residues(model, select)
    rows = _charged(residues)
    if len(rows) < 2:
        return []
    xyz = np.array([x[0] for x in rows], float)
    sign = np.array([x[1] for x in rows], int)
    owner = np.array([x[2] for x in rows], int)
    ij = cKDTree(xyz).query_pairs(SALT_MAX + _REACH, output_type="ndarray")
    if not len(ij):
        return []
    i, j = ij[:, 0], ij[:, 1]
    keep = (sign[i] != sign[j]) & _allowed(sel, names, residues, owner[i], owner[j])
    i, j = i[keep], j[keep]
    cation, anion = np.where(sign[i] > 0, i, j), np.where(sign[i] > 0, j, i)
    dist = np.linalg.norm(xyz[cation] - xyz[anion], axis=1)
    out = []
    for c, a, d in zip(cation, anion, map(_snap, dist)):
        if d > SALT_MAX:
            continue
        res = residues[owner[c]]
        note = HIS_NOTE if res.label.name == "HIS" and rows[c][3] != "N" else ""
        out.append(Interaction("salt", _key(res, rows[c][3]), _key(residues[owner[a]], rows[a][3]), d, None, note))
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
    angle: between the ring normal and the centroid-cation line. Sorted by `sort_key`."""
    sel, residues, names = _residues(model, select)
    rings, cations = _rings(residues), _cations(residues)
    if not rings or not cations:
        return []
    centres = np.array([x[0] for x in rings], float)
    tree = cKDTree(centres)
    near = tree.query_ball_point(np.array([x[0] for x in cations], float), CATION_PI_MAX + _REACH, return_sorted=True)
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
        length = float(np.linalg.norm(v))
        height = float(np.dot(v, normal))
        d, offset = _snap(length), _snap(np.linalg.norm(v - height * normal))
        if d == 0.0 or d > CATION_PI_MAX or offset > CATION_PI_OFFSET:
            continue
        angle = _snap(np.degrees(np.arccos(min(1.0, abs(height) / length))))
        out.append(Interaction("cation-pi", _key(residues[r], atom), _key(residues[q], ring), d, angle))
    return _sorted(out)


def find_pi_stacking(model, select) -> list[Interaction]:
    """The π-stacking ring pairs of `model` between the partners `select` names: centroids within STACK_MAX,
    parallel (normals within STACK_PARALLEL) or T-shaped (from STACK_T_MIN), offset at most STACK_OFFSET; never two
    nucleobases. angle: between the ring normals; note: "parallel" or "T-shaped". Sorted by `sort_key`."""
    sel, residues, names = _residues(model, select)
    rings = _rings(residues)
    if len(rings) < 2:
        return []
    ij = cKDTree(np.array([x[0] for x in rings], float)).query_pairs(STACK_MAX + _REACH, output_type="ndarray")
    if not len(ij):
        return []
    i, j = ij[:, 0], ij[:, 1]
    owner = np.array([x[2] for x in rings], int)
    nucleic = np.array([r.kind == NUCLEIC for r in residues], bool)
    keep = _allowed(sel, names, residues, owner[i], owner[j]) & ~(nucleic[owner[i]] & nucleic[owner[j]])
    out = []
    for p, q in zip(i[keep], j[keep]):
        c1, n1, r1, ring1 = rings[p]
        c2, n2, r2, ring2 = rings[q]
        v = c2 - c1
        d = _snap(np.linalg.norm(v))
        if d == 0.0 or d > STACK_MAX:
            continue
        angle = _snap(np.degrees(np.arccos(min(1.0, abs(float(np.dot(n1, n2)))))))
        if angle <= STACK_PARALLEL:
            note = PARALLEL
        elif angle >= STACK_T_MIN:
            note = T_SHAPED
        else:
            continue
        offset = _snap(min(np.linalg.norm(v - np.dot(v, n1) * n1), np.linalg.norm(v - np.dot(v, n2) * n2)))
        if offset > STACK_OFFSET:
            continue
        a, b = _ordered(_key(residues[r1], ring1), _key(residues[r2], ring2))
        out.append(Interaction("stack", a, b, d, angle, note))
    return _sorted(out)


def _hydrophobic_carbons(residues: list[_Residue]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """The carbons bonded only to carbon (bonds: heavy atoms of one residue within _COVALENT): their xyz, residue
    index and name, in model order."""
    rows = [(pos.tolist(), r, name) for r, res in enumerate(residues) for name, pos in res.atoms.items()]
    if not rows:
        return np.zeros((0, 3)), np.zeros(0, int), []
    xyz = np.array([x[0] for x in rows], float)
    owner = np.array([x[1] for x in rows], int)
    carbon = np.array([residues[r].elements[name] == "C" for _, r, name in rows], bool)
    hetero_bonded = np.zeros(len(rows), bool)
    ij = cKDTree(xyz).query_pairs(_COVALENT + _REACH, output_type="ndarray")
    if len(ij):
        i, j = ij[:, 0], ij[:, 1]
        same = owner[i] == owner[j]
        i, j = i[same], j[same]
        bonded = _snaps(np.linalg.norm(xyz[i] - xyz[j], axis=1)) <= _COVALENT
        i, j = i[bonded], j[bonded]
        hetero_bonded[i[~carbon[j]]] = True
        hetero_bonded[j[~carbon[i]]] = True
    keep = np.flatnonzero(carbon & ~hetero_bonded)
    return xyz[keep], owner[keep], [rows[k][2] for k in keep]


def find_hydrophobic(model, select) -> list[Interaction]:
    """The hydrophobic contacts of `model` between the partners `select` names: of each residue pair with two
    carbons bonded only to carbon within HYDROPHOBIC_MAX, the shortest such pair (ties: the first atoms by name);
    `a` the residue that sorts first. Sorted by `sort_key`."""
    sel, residues, names = _residues(model, select)
    xyz, owner, atom = _hydrophobic_carbons(residues)
    if len(xyz) < 2:
        return []
    ij = cKDTree(xyz).query_pairs(HYDROPHOBIC_MAX + _REACH, output_type="ndarray")
    if not len(ij):
        return []
    i, j = ij[:, 0], ij[:, 1]
    keep = _allowed(sel, names, residues, owner[i], owner[j])
    i, j = i[keep], j[keep]
    dist = _snaps(np.linalg.norm(xyz[i] - xyz[j], axis=1))
    best: dict[tuple, tuple] = {}
    for p, q, d in zip(i, j, dist.tolist()):
        if d > HYDROPHOBIC_MAX:
            continue
        a, b = _ordered(_key(residues[owner[p]], atom[p]), _key(residues[owner[q]], atom[q]))
        pair = (a[:2], b[:2])
        rank = (d, a[2], b[2])
        if pair not in best or rank < best[pair][0]:
            best[pair] = (rank, a, b)
    return _sorted([Interaction("hydrophobic", a, b, rank[0], None) for rank, a, b in best.values()])


# The independent re-check. Nothing below calls a finder or its helpers (_residues, _rings, _cations, _charged,
# _hydrophobic_carbons, _heavy, _kind): residues are read from the model again, rings are planes fitted by an
# eigen-decomposition (the finders take an SVD), angles come from atan2 (the finders take arccos), and distances from
# math.dist. Only the limits and the names of the rule's atoms are shared.

_V_SLACK = 0.5 / _GRID + 1e-5  # a measured value the finders' grid could place on a limit passes it
_V_BONDED = _COVALENT + 0.5 / _GRID - 1e-5  # surely bonded, as the finders judge it
_V_CATIONS = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2"}}
_V_ANIONS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
_V_RINGS = {
    ("PHE", 6): {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    ("TYR", 6): {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    ("TRP", 5): {"CG", "CD1", "NE1", "CE2", "CD2"},
    ("TRP", 6): {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    ("HIS", 5): {"CG", "ND1", "CD2", "CE1", "NE2"},
}
_V_BASE_RINGS = {6: {"N1", "C2", "N3", "C4", "C5", "C6"}, 5: {"C4", "C5", "N7", "C8", "N9"}}
_V_HBOND_ELEMENTS = {"N", "O", "S"}


class _VRes(NamedTuple):
    name: str
    protein: bool
    atoms: dict  # name -> (element, (x, y, z)), the first conformer
    n_term: bool
    c_term: bool


def _v_polymer(res) -> str | None:
    """ "protein", "nucleic" or None (water, ligand, ion)."""
    if res.is_water() or res.entity_type not in (gemmi.EntityType.Polymer, gemmi.EntityType.Unknown):
        return None
    info = gemmi.find_tabulated_residue(res.name)
    if info is not None and info.is_amino_acid():
        return "protein"
    if info is not None and info.is_nucleic_acid():
        return "nucleic"
    if res.entity_type == gemmi.EntityType.Polymer:
        names = {a.name for a in res}
        if {"C1'", "C4'", "O4'"} <= names:
            return "nucleic"
        if {"N", "CA", "C"} <= names:
            return "protein"
    return None


def _v_book(model) -> dict[tuple[str, str], _VRes]:
    """(chain, residue tag) -> its polymer residue, the first of a number."""
    book: dict[tuple[str, str], _VRes] = {}
    for chain in model:
        order: list[tuple[tuple[str, str], str, object]] = []
        seen: set[tuple[str, str]] = set()
        for res in chain:
            kind = _v_polymer(res)
            key = (chain.name, f"{res.seqid.num}{res.seqid.icode.strip()}")
            if kind is None or key in book or key in seen:
                continue
            seen.add(key)
            order.append((key, kind, res))
        proteins = [k for k, kind, _ in order if kind == "protein"]
        for key, kind, res in order:
            atoms: dict[str, tuple[str, tuple[float, float, float]]] = {}
            for at in res:
                if at.element.name not in ("H", "D") and at.name not in atoms:
                    atoms[at.name] = (at.element.name, (at.pos.x, at.pos.y, at.pos.z))
            protein = kind == "protein"
            n_term = protein and key == proteins[0]
            c_term = protein and (key == proteins[-1] or "OXT" in atoms)
            book[key] = _VRes(res.name, protein, atoms, n_term, c_term)
    return book


def _v_sub(p, q):
    return (p[0] - q[0], p[1] - q[1], p[2] - q[2])


def _v_dot(p, q):
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]


def _v_mean(points):
    n = len(points)
    return (sum(p[0] for p in points) / n, sum(p[1] for p in points) / n, sum(p[2] for p in points) / n)


def _v_ring(res: _VRes, label: str):
    """(centroid, unit normal) of the ring `label` ("ring:K") names in `res`, or None."""
    if not label.startswith("ring:") or not label[5:].isdigit():
        return None
    size = int(label[5:])
    atoms = _V_RINGS.get((res.name, size)) if res.protein else _V_BASE_RINGS.get(size)
    if not atoms or not atoms <= res.atoms.keys():
        return None
    points = [res.atoms[a][1] for a in sorted(atoms)]
    centre = _v_mean(points)
    spread = np.array([_v_sub(p, centre) for p in points], float)
    _, vectors = np.linalg.eigh(spread.T @ spread)
    normal = vectors[:, 0]
    return centre, tuple(float(x) for x in normal / math.sqrt(float(normal @ normal)))


def _v_offset(v, normal) -> float:
    """The distance from the end of v to the plane's normal axis through its start."""
    d2, h = _v_dot(v, v), _v_dot(v, normal)
    return math.sqrt(max(0.0, d2 - h * h))


def _v_where(p) -> str:
    return f"{p[0]}:{p[1]}:{p[2]}"


def _v_atom(res: _VRes, p, bad: list[str]):
    got = res.atoms.get(p[2])
    if got is None:
        bad.append(f"no atom {p[2]} in {res.name} {p[0]}:{p[1]}")
    return got


def _v_over(value: float, limit: float, what: str, bad: list[str]) -> None:
    if value > limit + _V_SLACK:
        bad.append(f"{what} {value:.3f} over {limit}")


def _v_salt(ra: _VRes, rb: _VRes, x: Interaction, bad: list[str]):
    na, nb = x.a[2], x.b[2]
    if not (ra.protein and (na in _V_CATIONS.get(ra.name, ()) or (na == "N" and ra.n_term))):
        bad.append(f"{na} of {ra.name} {x.a[0]}:{x.a[1]} is not a cationic N")
    if rb.protein:
        anion = nb in _V_ANIONS.get(rb.name, ()) or (nb in ("O", "OXT") and rb.c_term)
    else:
        anion = nb in ("OP1", "OP2")
    if not anion:
        bad.append(f"{nb} of {rb.name} {x.b[0]}:{x.b[1]} is not an anionic O")
    note = HIS_NOTE if ra.name == "HIS" and na != "N" else ""
    if x.note != note:
        bad.append(f"note {x.note!r}, the rule's {note!r}")
    pa, pb = _v_atom(ra, x.a, bad), _v_atom(rb, x.b, bad)
    if pa is None or pb is None:
        return None
    d = math.dist(pa[1], pb[1])
    _v_over(d, SALT_MAX, "distance", bad)
    return d, None


def _v_hbond(ra: _VRes, rb: _VRes, x: Interaction, bad: list[str]):
    pa, pb = _v_atom(ra, x.a, bad), _v_atom(rb, x.b, bad)
    if pa is None or pb is None:
        return None
    for p, got in ((x.a, pa), (x.b, pb)):
        if got[0] not in _V_HBOND_ELEMENTS:
            bad.append(f"{_v_where(p)} ({got[0]}) is not an H-bond atom (N, O or S)")
    d = math.dist(pa[1], pb[1])
    if not MIN_DIST - _V_SLACK <= d <= MAX_DIST + _V_SLACK:
        bad.append(f"H-bond distance {d:.3f} outside {MIN_DIST}-{MAX_DIST}")
    return d, None


def _v_cation_pi(ra: _VRes, rb: _VRes, x: Interaction, bad: list[str]):
    na = x.a[2]
    cation = None
    if ra.protein and ra.name == "LYS" and na == "NZ":
        got = _v_atom(ra, x.a, bad)
        cation = got[1] if got else None
    elif ra.protein and ra.name == "ARG" and na == "centre":
        if {"NE", "NH1", "NH2"} <= ra.atoms.keys():
            cation = _v_mean([ra.atoms[n][1] for n in ("NE", "NH1", "NH2")])
        else:
            bad.append(f"ARG {x.a[0]}:{x.a[1]} lacks a guanidinium atom")
    else:
        bad.append(f"{na} of {ra.name} {x.a[0]}:{x.a[1]} is not a cation (Lys NZ or the Arg centre)")
    ring = _v_ring(rb, x.b[2])
    if ring is None:
        bad.append(f"{rb.name} {x.b[0]}:{x.b[1]} has no whole {x.b[2]}")
    if cation is None or ring is None:
        return None
    v = _v_sub(cation, ring[0])
    h, offset = abs(_v_dot(v, ring[1])), _v_offset(v, ring[1])
    d = math.dist(cation, ring[0])
    _v_over(d, CATION_PI_MAX, "distance", bad)
    _v_over(offset, CATION_PI_OFFSET, "offset", bad)
    return d, math.degrees(math.atan2(offset, h))


def _v_stack(ra: _VRes, rb: _VRes, x: Interaction, bad: list[str]):
    if not ra.protein and not rb.protein:
        bad.append("two bases: base-base stacking within nucleic acid is left out")
    r1, r2 = _v_ring(ra, x.a[2]), _v_ring(rb, x.b[2])
    for res, p, ring in ((ra, x.a, r1), (rb, x.b, r2)):
        if ring is None:
            bad.append(f"{res.name} {p[0]}:{p[1]} has no whole {p[2]}")
    if r1 is None or r2 is None:
        return None
    (c1, n1), (c2, n2) = r1, r2
    v = _v_sub(c2, c1)
    cross = (n1[1] * n2[2] - n1[2] * n2[1], n1[2] * n2[0] - n1[0] * n2[2], n1[0] * n2[1] - n1[1] * n2[0])
    angle = math.degrees(math.atan2(math.sqrt(_v_dot(cross, cross)), abs(_v_dot(n1, n2))))
    d = math.dist(c1, c2)
    _v_over(d, STACK_MAX, "distance", bad)
    _v_over(min(_v_offset(v, n1), _v_offset(v, n2)), STACK_OFFSET, "offset", bad)
    parallel = angle <= STACK_PARALLEL + _V_SLACK
    t_shaped = angle >= STACK_T_MIN - _V_SLACK
    if not parallel and not t_shaped:
        bad.append(f"angle {angle:.3f} between the ring normals: neither parallel nor T-shaped")
    elif x.note != (PARALLEL if parallel else T_SHAPED):
        bad.append(f"marked {x.note!r} at an angle of {angle:.3f}")
    return d, angle


def _v_hydrophobic(ra: _VRes, rb: _VRes, x: Interaction, bad: list[str]):
    pa, pb = _v_atom(ra, x.a, bad), _v_atom(rb, x.b, bad)
    for res, p, got in ((ra, x.a, pa), (rb, x.b, pb)):
        if got is None:
            continue
        if got[0] != "C":
            bad.append(f"{_v_where(p)} ({got[0]}) is not a carbon")
            continue
        hetero = sorted(
            n for n, (e, xyz) in res.atoms.items() if n != p[2] and e != "C" and math.dist(xyz, got[1]) <= _V_BONDED
        )
        if hetero:
            bad.append(f"{_v_where(p)} is bonded to {', '.join(hetero)}")
    if pa is None or pb is None:
        return None
    d = math.dist(pa[1], pb[1])
    _v_over(d, HYDROPHOBIC_MAX, "distance", bad)
    return d, None


_V_RULES = {
    "hbond": _v_hbond,
    "salt": _v_salt,
    "stack": _v_stack,
    "cation-pi": _v_cation_pi,
    "hydrophobic": _v_hydrophobic,
}


def verify_interactions(model, drawn) -> list[str]:
    """Re-check the lines a figure draws against `model` (a gemmi Model, or what io.read_model returns): `drawn` is
    a list of Interaction (anything with kind, a, b, distance, angle, note). Each line's partners are looked up in
    the model again and measured by this function's own code (see the comment above), and must meet their kind's
    rule: the right atoms for the kind (cationic N, anionic O, a whole ring, a carbon bonded only to carbon, an N, O
    or S for an H-bond), two different residues, the distance (H-bonds: within MIN_DIST-MAX_DIST) and offset within
    their limits, a stack's angle in its window and named by its note, and the drawn distance and angle the measured
    ones (to the finders' 1e-3 grid). Returns one text per line that fails, "<kind> <a> - <b>: <reasons>" with each
    partner chain:tag:atom, in the order drawn; empty when every line holds."""
    if isinstance(model, tuple):
        model = model[0]
    book = _v_book(model)
    out = []
    for x in drawn:
        bad: list[str] = []
        ra, rb = book.get(tuple(x.a[:2])), book.get(tuple(x.b[:2]))
        for p, res in ((x.a, ra), (x.b, rb)):
            if res is None:
                bad.append(f"no residue {p[0]}:{p[1]} in the model")
        if tuple(x.a[:2]) == tuple(x.b[:2]):
            bad.append("a residue with itself")
        rule = _V_RULES.get(x.kind)
        if rule is None:
            bad.append(f"no rule for kind {x.kind!r}")
        elif ra is not None and rb is not None:
            measured = rule(ra, rb, x, bad)
            if measured is not None:
                d, angle = measured
                if abs(x.distance - d) > _V_SLACK:
                    bad.append(f"distance drawn {x.distance}, measured {d:.3f}")
                if angle is None and x.angle is not None:
                    bad.append(f"an angle ({x.angle}) where the kind has none")
                elif angle is not None and (x.angle is None or abs(x.angle - angle) > _V_SLACK):
                    bad.append(f"angle drawn {x.angle}, measured {angle:.3f}")
        if bad:
            out.append(f"{x.kind} {_v_where(x.a)} - {_v_where(x.b)}: {'; '.join(bad)}")
    return out
