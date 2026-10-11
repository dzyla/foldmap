# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Hydrogen bonds between protein and nucleic acid, or between chains, from the heavy atoms of a model.

Most entries carry no hydrogens, so a bond is a donor and an acceptor heavy atom 2.5-3.5 Å apart, and, when the
donor's antecedent (the heavy atom(s) it is bonded to) is known, an angle antecedent-donor-acceptor of at least 90°.
Asn/Gln amide and His ring atoms may be flipped in the model, so they count as donor and acceptor both. Waters and
ligands are left out.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .model import ResLabel

MIN_DIST = 2.5  # Å, donor-acceptor
MAX_DIST = 3.5
MIN_ANGLE = 90.0  # degrees, antecedent-donor-acceptor

D, A, DA = "D", "A", "DA"  # donor, acceptor, either

# atom -> (role, antecedent atoms; a donor's antecedent is their mean position)
_BACKBONE = {"N": (D, ("CA",)), "O": (A, ()), "OXT": (A, ())}
_SIDE_CHAIN = {
    "ARG": {"NE": (D, ("CD", "CZ")), "NH1": (D, ("CZ",)), "NH2": (D, ("CZ",))},
    "LYS": {"NZ": (D, ("CE",))},
    "ASN": {"OD1": (DA, ("CG",)), "ND2": (DA, ("CG",))},
    "GLN": {"OE1": (DA, ("CD",)), "NE2": (DA, ("CD",))},
    "HIS": {"ND1": (DA, ("CG", "CE1")), "NE2": (DA, ("CD2", "CE1"))},
    "SER": {"OG": (DA, ("CB",))},
    "THR": {"OG1": (DA, ("CB",))},
    "TYR": {"OH": (DA, ("CZ",))},
    "ASP": {"OD1": (A, ()), "OD2": (A, ())},
    "GLU": {"OE1": (A, ()), "OE2": (A, ())},
    "TRP": {"NE1": (D, ("CD1", "CE2"))},
    "CYS": {"SG": (DA, ("CB",))},  # weak either way
    "MET": {"SD": (A, ())},
}
_NA_BACKBONE = {
    "OP1": (A, ()),
    "OP2": (A, ()),
    "OP3": (A, ()),
    "O5'": (A, ()),
    "O4'": (A, ()),
    "O3'": (A, ()),
    "O2'": (DA, ("C2'",)),
}
# bases by their exocyclic atoms, so modified nucleotides use their parent's table
_GUANINE = {"N1": (D, ("C2", "C6")), "N2": (D, ("C2",)), "O6": (A, ()), "N3": (A, ()), "N7": (A, ())}
_ADENINE = {"N6": (D, ("C6",)), "N1": (A, ()), "N3": (A, ()), "N7": (A, ())}
_CYTOSINE = {"N4": (D, ("C4",)), "N3": (A, ()), "O2": (A, ())}
# N1 is the glycosidic nitrogen of uracil and thymine, a donor only in pseudouridine (dropped when bonded to C1')
_URACIL = {"N3": (D, ("C2", "C4")), "N1": (D, ("C2", "C6")), "O2": (A, ()), "O4": (A, ())}
_GLYCOSIDIC = 1.7  # Å, an N this close to C1' carries no hydrogen

PROTEIN, NUCLEIC = 0, 1


@dataclass(frozen=True)
class Selection:
    """Which bonds to find. kind: "off", "na" (protein with nucleic acid), "interchain" (any two different chains)
    or "chains" (only between the chains of each pair in `pairs`, any modality)."""

    kind: str
    pairs: tuple[tuple[str, str], ...] = ()

    def __str__(self) -> str:
        if self.kind == "chains":
            return "chains:" + ",".join(f"{a}+{b}" for a, b in self.pairs)
        return self.kind


_FORMS = "off, na, interchain or chains:A+B,C+D"


def parse_selection(text: str | Selection) -> Selection:
    """The selection a `hbonds` value names: off | na | interchain | chains:A+B,C+D (chain names as in the
    model)."""
    if isinstance(text, Selection):
        return text
    value = str(text).strip()
    if value in ("off", "na", "interchain"):
        return Selection(value)
    if value == "ligands":
        raise ValueError("hbonds: 'ligands' is reserved for later; use " + _FORMS)
    if not value.startswith("chains:"):
        raise ValueError(f"hbonds: {value!r} is not a selection; use {_FORMS}")
    pairs = []
    for item in value[len("chains:") :].split(","):
        names = [n.strip() for n in item.split("+")]
        if len(names) != 2 or not all(names):
            raise ValueError(f"hbonds: {item.strip()!r} in {value!r} is not a chain pair like A+B")
        if names[0] == names[1]:
            raise ValueError(f"hbonds: {item.strip()!r} in {value!r} must name two different chains")
        pairs.append((names[0], names[1]))
    return Selection("chains", tuple(pairs))


class Atom(NamedTuple):
    res: ResLabel
    name: str

    def __str__(self) -> str:
        return f"{self.res.chain}:{self.res.name}{self.res.tag} {self.name}"


class HBond(NamedTuple):
    donor: Atom
    acceptor: Atom
    distance: float  # Å

    def __str__(self) -> str:
        return f"{self.donor} → {self.acceptor} {self.distance:.2f} Å"


def _kind(res) -> int | None:
    """PROTEIN, NUCLEIC, or None (water, ligand, ion)."""
    if res.is_water():
        return None
    info = gemmi.find_tabulated_residue(res.name)
    if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
        return PROTEIN if info.is_amino_acid() else NUCLEIC
    if res.entity_type != gemmi.EntityType.Polymer:
        return None
    names = {a.name for a in res}
    if {"C1'", "C4'", "O4'"} <= names:
        return NUCLEIC
    if {"N", "CA", "C"} <= names:
        return PROTEIN
    return None


def _table(res, kind: int, atoms: dict) -> dict:
    if kind == PROTEIN:
        table = dict(_BACKBONE)
        if res.name == "PRO":
            del table["N"]
        table.update(_SIDE_CHAIN.get(res.name, {}))
        return table
    table = dict(_NA_BACKBONE)
    for marker, base in (("O6", _GUANINE), ("N6", _ADENINE), ("N4", _CYTOSINE), ("O4", _URACIL)):
        if marker in atoms:
            table.update(base)
            break
    c1 = atoms.get("C1'")
    if c1 is not None:  # the glycosidic N of any base is no donor
        for name in list(table):
            if name.startswith("N") and name in atoms and atoms[name].dist(c1) < _GLYCOSIDIC:
                del table[name]
    return table


def _polar_atoms(model, chains: set[str] | None):
    """Arrays of the donor and acceptor atoms: xyz, antecedent xyz (nan when unknown), donor?, acceptor?,
    residue index, chain index, kind; with the residue labels, atom names and chain names."""
    xyz, ante, don, acc, owner, chain_of, kind_of, names = [], [], [], [], [], [], [], []
    labels: list[ResLabel] = []
    chain_names: list[str] = []
    for chain in model:
        if chains is not None and chain.name not in chains:
            continue
        if chain.name not in chain_names:
            chain_names.append(chain.name)
        c = chain_names.index(chain.name)
        for res in chain:
            kind = _kind(res)
            if kind is None:
                continue
            atoms: dict[str, gemmi.Position] = {}
            for atom in res:  # first conformer of each atom
                if atom.element.name not in ("H", "D"):
                    atoms.setdefault(atom.name, atom.pos)
            table = _table(res, kind, atoms)
            r = len(labels)
            labels.append(ResLabel(chain.name, res.seqid.num, res.seqid.icode.strip(), res.name))
            for name, (role, antecedents) in table.items():
                pos = atoms.get(name)
                if pos is None:
                    continue
                xyz.append(pos.tolist())
                found = [atoms[a].tolist() for a in antecedents if a in atoms]
                ante.append(np.mean(found, axis=0) if antecedents and len(found) == len(antecedents) else [np.nan] * 3)
                don.append(role in (D, DA))
                acc.append(role in (A, DA))
                owner.append(r)
                chain_of.append(c)
                kind_of.append(kind)
                names.append(name)
    return (
        np.array(xyz, float).reshape(-1, 3),
        np.array(ante, float).reshape(-1, 3),
        np.array(don, bool),
        np.array(acc, bool),
        np.array(owner, int),
        np.array(chain_of, int),
        np.array(kind_of, int),
        labels,
        names,
        chain_names,
    )


def _donates(xyz, ante, don, acc, d, a):
    """Atoms d donate to atoms a: d a donor, a an acceptor, and the antecedent-donor-acceptor angle >= MIN_ANGLE
    (true when the antecedent is unknown)."""
    to_ante = ante[d] - xyz[d]
    to_acc = xyz[a] - xyz[d]
    with np.errstate(divide="ignore", invalid="ignore"):
        cos = np.einsum("ij,ij->i", to_ante, to_acc) / (
            np.linalg.norm(to_ante, axis=1) * np.linalg.norm(to_acc, axis=1)
        )
        angle_ok = np.isnan(cos) | (cos <= np.cos(np.radians(MIN_ANGLE)) + 1e-12)
    return don[d] & acc[a] & angle_ok


def find_hbonds(model, select) -> list[HBond]:
    """The hydrogen bonds of `model` (a gemmi Model, or what io.read_model returns) that `select` (a Selection or
    its text, see parse_selection) asks for, sorted by donor then acceptor. Each donor-acceptor atom pair is
    reported once; when either atom could donate, the first in model order is the donor."""
    if isinstance(model, tuple):
        model = model[0]
    sel = parse_selection(select)
    if sel.kind == "off":
        return []
    chains = None
    if sel.kind == "chains":
        present = {c.name for c in model}
        missing = sorted({n for p in sel.pairs for n in p} - present)
        if missing:
            raise ValueError(f"hbonds: no chain {', '.join(missing)} in the model; it has {', '.join(sorted(present))}")
        chains = {n for p in sel.pairs for n in p}
    xyz, ante, don, acc, owner, chain_of, kind_of, labels, names, chain_names = _polar_atoms(model, chains)
    if len(xyz) < 2:
        return []
    ij = cKDTree(xyz).query_pairs(MAX_DIST, output_type="ndarray")
    if not len(ij):
        return []
    i, j = ij[:, 0], ij[:, 1]
    if sel.kind == "na":
        keep = kind_of[i] != kind_of[j]
    elif sel.kind == "interchain":
        keep = chain_of[i] != chain_of[j]
    else:
        allowed = np.zeros((len(chain_names),) * 2, bool)
        for a, b in sel.pairs:
            allowed[chain_names.index(a), chain_names.index(b)] = True
            allowed[chain_names.index(b), chain_names.index(a)] = True
        keep = allowed[chain_of[i], chain_of[j]]
    keep &= owner[i] != owner[j]
    i, j = i[keep], j[keep]
    dist = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    keep = dist >= MIN_DIST
    i, j, dist = i[keep], j[keep], dist[keep]
    lo, hi = np.minimum(i, j), np.maximum(i, j)  # model order
    forward = _donates(xyz, ante, don, acc, lo, hi)
    backward = ~forward & _donates(xyz, ante, don, acc, hi, lo)
    out = [
        HBond(Atom(labels[owner[d]], names[d]), Atom(labels[owner[a]], names[a]), float(x))
        for d, a, x in zip(
            np.where(forward, lo, hi)[forward | backward],
            np.where(forward, hi, lo)[forward | backward],
            dist[forward | backward],
        )
    ]
    return sorted(out)
