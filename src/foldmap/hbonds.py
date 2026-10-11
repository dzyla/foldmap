# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Hydrogen bonds between protein and nucleic acid, or between chains, from the heavy atoms of a model.

A heavy-atom approximation: most entries carry no hydrogens, so a bond is a donor and an acceptor heavy atom
2.5-3.5 Å apart where the donor's hydrogen could point at the acceptor. That is judged by the angle
antecedent-donor-acceptor, the antecedent being the mean position of the heavy atoms the donor is bonded to in its
residue (CB for Ser OG; C2 and C6 for G N1, so the ring bisector):

- donors whose hydrogen rotates (Ser OG, Thr OG1, Tyr OH, Cys SG, Lys NZ, O2', a terminal O3'/O5' hydroxyl, the
  N-terminal amine): at least 60°;
- sp2 donors (backbone and amide N, ring N-H, amino N, Arg, His, Trp): at least 90°;
- a donor whose antecedent atoms are not all in the model passes.

Asn/Gln amide and His ring atoms may be flipped in the model, so they count as donor and acceptor both. Roles follow
the residue's own bonds (heavy atoms within 1.7 Å), so modified residues keep only what their chemistry allows: an
imino or hydroxyl donor (ring N-H, O2', Ser OG) bonded to a heavy atom outside its antecedents has no hydrogen left
(U N1 is glycosidic, PSU N1 donates; 2'-O-methyl O2' only accepts); an amino donor (G N2, A N6, C N4) keeps one
until two substituents (2MG N2 donates, M2G N2 does not), or until closed into a ring, where it is an imine that
only accepts (wybutosine N2); a ring N with a third heavy substituent (7MG N7, 1MA N1) accepts nothing. A
terminal O3' or O5' (no phosphate on it) is a hydroxyl, donor and acceptor. Bases are typed by their exocyclic
atoms (O6 guanine, N6 adenine, N4 cytosine, O4 uracil and thymine), so a modified nucleotide uses its parent's
table. Of residues sharing a number (microheterogeneity) the first is kept, and of each atom its first conformer.
Waters and ligands (any non-polymer entity) are left out.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .model import ResLabel

MIN_DIST = 2.5  # Å, donor-acceptor
MAX_DIST = 3.5  # sulfur H-bonds (Cys SG, Met SD) usually run 3.3-3.9 Å; capped here at 3.5 Å by design
SP2_ANGLE = 90.0  # degrees, antecedent-donor-acceptor, for a donor whose H lies in a plane
ROTATING_ANGLE = 60.0  # ... for a donor whose H rotates about its bond to the antecedent
_BOND = 1.7  # Å, heavy atoms bonded
_PEPTIDE_BOND = 2.5  # Å, C(i-1)-N(i), as io

D, A, DA = "D", "A", "DA"  # donor, acceptor, either


class _Site(NamedTuple):
    """A donor or acceptor atom. ante: the heavy atoms a donor is bonded to (their mean is its antecedent);
    h: hydrogens it carries, each extra heavy substituent taking one; rotating: its H rotates (60° rule)."""

    role: str
    ante: tuple[str, ...] = ()
    h: int = 1
    rotating: bool = False


_ACC = _Site(A)
_BACKBONE = {"N": _Site(D, ("CA",)), "O": _ACC, "OXT": _ACC}  # Pro N: CD takes its H
_N_TERMINAL = _Site(D, ("CA",), 3, True)
_SIDE_CHAIN = {
    "ARG": {"NE": _Site(D, ("CD", "CZ")), "NH1": _Site(D, ("CZ",), 2), "NH2": _Site(D, ("CZ",), 2)},
    "LYS": {"NZ": _Site(D, ("CE",), 3, True)},
    "ASN": {"OD1": _Site(DA, ("CG",), 2), "ND2": _Site(DA, ("CG",), 2)},
    "GLN": {"OE1": _Site(DA, ("CD",), 2), "NE2": _Site(DA, ("CD",), 2)},
    "HIS": {"ND1": _Site(DA, ("CG", "CE1")), "NE2": _Site(DA, ("CD2", "CE1"))},
    "SER": {"OG": _Site(DA, ("CB",), 1, True)},
    "THR": {"OG1": _Site(DA, ("CB",), 1, True)},
    "TYR": {"OH": _Site(DA, ("CZ",), 1, True)},
    "ASP": {"OD1": _ACC, "OD2": _ACC},
    "GLU": {"OE1": _ACC, "OE2": _ACC},
    "TRP": {"NE1": _Site(D, ("CD1", "CE2"))},
    "CYS": {"SG": _Site(DA, ("CB",), 1, True)},  # weak either way
    "MET": {"SD": _ACC},
}
_NA_BACKBONE = {"OP1": _ACC, "OP2": _ACC, "OP3": _ACC, "O5'": _ACC, "O4'": _ACC, "O3'": _ACC}
_NA_BACKBONE["O2'"] = _Site(DA, ("C2'",), 1, True)
_HYDROXYL = {"O5'": (_Site(DA, ("C5'",), 1, True), "P"), "O3'": (_Site(DA, ("C3'",), 1, True), "P")}
_GUANINE = {"N1": _Site(D, ("C2", "C6")), "N2": _Site(D, ("C2",), 2), "O6": _ACC, "N3": _ACC, "N7": _ACC}
_ADENINE = {"N6": _Site(D, ("C6",), 2), "N1": _ACC, "N3": _ACC, "N7": _ACC}
_CYTOSINE = {"N4": _Site(D, ("C4",), 2), "N3": _ACC, "O2": _ACC}
# N1 of uracil and thymine is glycosidic (C1' takes its H); in pseudouridine it carries one
_URACIL = {"N3": _Site(D, ("C2", "C4")), "N1": _Site(D, ("C2", "C6")), "O2": _ACC, "O4": _ACC}
_BASES = (("O6", _GUANINE), ("N6", _ADENINE), ("N4", _CYTOSINE), ("O4", _URACIL))

PROTEIN, NUCLEIC = 0, 1


_KINDS = ("off", "na", "interchain", "chains")
_FORMS = "off, na, interchain or chains:A+B,C+D"


@dataclass(frozen=True)
class Selection:
    """Which bonds to find. kind: "off", "na" (protein with nucleic acid), "interchain" (any two different chains)
    or "chains" (only between the chains of each pair in `pairs`, any modality). Each pair is kept in sorted order
    and once, so A+B and B+A are one selection."""

    kind: str
    pairs: tuple[tuple[str, str], ...] = ()

    def __post_init__(self):
        if self.kind not in _KINDS:
            raise ValueError(f"hbonds: no selection kind {self.kind!r}; use {_FORMS}")
        if (self.kind == "chains") != bool(self.pairs):
            raise ValueError(f"hbonds: chain pairs go with kind 'chains' only, and it needs at least one; use {_FORMS}")
        pairs: list[tuple[str, str]] = []
        for pair in self.pairs:
            if len(pair) != 2 or not all(pair) or pair[0] == pair[1]:
                raise ValueError(f"hbonds: {pair!r} is not a pair of two different chains")
            key = tuple(sorted(pair))
            if key not in pairs:
                pairs.append(key)
        object.__setattr__(self, "pairs", tuple(pairs))

    def __str__(self) -> str:
        if self.kind == "chains":
            return "chains:" + ",".join(f"{a}+{b}" for a, b in self.pairs)
        return self.kind


def parse_selection(text: str | Selection) -> Selection:
    """The selection a `hbonds` value names: off | na | interchain | chains:A+B,C+D. Keywords in any case; chain
    names as in the model (case kept)."""
    if isinstance(text, Selection):
        return text
    value = str(text).strip()
    word = value.lower()
    if word in ("off", "na", "interchain"):
        return Selection(word)
    if word == "ligands":
        raise ValueError("hbonds: 'ligands' is reserved for later; use " + _FORMS)
    if not word.startswith("chains:"):
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
    """PROTEIN, NUCLEIC, or None (water, ligand, ion: any residue of a non-polymer entity, even of a polymer
    residue type; with no entity types in the file, the tabulated type decides)."""
    if res.entity_type not in (gemmi.EntityType.Polymer, gemmi.EntityType.Unknown) or res.is_water():
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


def _heavy(res) -> dict[str, gemmi.Position]:
    """Heavy atoms of a residue by name, the first conformer of each."""
    atoms: dict[str, gemmi.Position] = {}
    for atom in res:
        if atom.element.name not in ("H", "D"):
            atoms.setdefault(atom.name, atom.pos)
    return atoms


def _linked(a: gemmi.Position | None, b: gemmi.Position | None, cutoff: float = _BOND) -> bool:
    return a is not None and b is not None and a.dist(b) <= cutoff


def _table(kind: int, name: str, atoms: dict, prev: dict | None, nxt: dict | None) -> dict[str, _Site]:
    """The donor and acceptor sites a residue's table gives, before its own bonds are looked at."""
    if kind == PROTEIN:
        table = dict(_BACKBONE)
        table.update(_SIDE_CHAIN.get(name, {}))
        if not _linked((prev or {}).get("C"), atoms.get("N"), _PEPTIDE_BOND):
            table["N"] = _N_TERMINAL
        return table
    table = dict(_NA_BACKBONE)
    for marker, base in _BASES:
        if marker in atoms:
            table.update(base)
            break
    for oxygen, (site, phosphorus) in _HYDROXYL.items():  # no phosphate on it: a hydroxyl
        holder = atoms if oxygen == "O5'" else (nxt or {})
        if not _linked(atoms.get(oxygen), holder.get(phosphorus)):
            table[oxygen] = site
    return table


def _in_ring(bonded: np.ndarray, k: int) -> bool:
    """Atom k lies on a ring: two of its neighbours are joined without it."""
    first, *others = [m for m in np.flatnonzero(bonded[k]) if m != k]
    if not others:
        return False
    seen, todo = {k, first}, [first]
    while todo:
        for m in np.flatnonzero(bonded[todo.pop()]):
            if m in others:
                return True
            if m not in seen:
                seen.add(m)
                todo.append(m)
    return False


def _roles(table: dict[str, _Site], atoms: dict) -> dict[str, _Site]:
    """The sites as the residue's own bonds leave them: each extra heavy substituent on a donor takes one of its
    hydrogens; a ring N with three heavy neighbours has no lone pair to accept with; an amino N closed into a
    ring (wybutosine N2) is a ring imine, an acceptor only."""
    names = list(atoms)
    xyz = np.array([atoms[n].tolist() for n in names])
    bonded = np.linalg.norm(xyz[:, None] - xyz[None], axis=2) <= _BOND
    out = {}
    for name, site in table.items():
        if name not in atoms:
            continue
        k = names.index(name)
        neighbours = {names[m] for m in np.flatnonzero(bonded[k]) if m != k}
        donor, acceptor = site.role in (D, DA), site.role in (A, DA)
        if site.h == 2 and name.startswith("N") and _in_ring(bonded, k):
            donor, acceptor = False, True
        if donor and len(neighbours - set(site.ante)) >= site.h:
            donor = False
        if acceptor and name.startswith("N") and len(neighbours) >= 3:
            acceptor = False
        if donor or acceptor:
            out[name] = site._replace(role=DA if donor and acceptor else D if donor else A)
    return out


@dataclass
class _Sites:
    """The donor and acceptor atoms of a model, as arrays: xyz, antecedent xyz (nan when unknown), the largest
    cosine of the antecedent-donor-acceptor angle allowed, donor?, acceptor?, residue index (into labels), chain
    index (into chains), PROTEIN or NUCLEIC; and the atom names."""

    xyz: np.ndarray
    ante: np.ndarray
    max_cos: np.ndarray
    donor: np.ndarray
    acceptor: np.ndarray
    owner: np.ndarray
    chain: np.ndarray
    kind: np.ndarray
    labels: list[ResLabel]
    names: list[str]
    chains: list[str]


def _sites(model, chains: set[str] | None) -> _Sites:
    rows: list[tuple] = []
    labels: list[ResLabel] = []
    chain_names: list[str] = []
    for chain in model:
        if chains is not None and chain.name not in chains:
            continue
        if chain.name not in chain_names:
            chain_names.append(chain.name)
        c = chain_names.index(chain.name)
        kept, seen = [], set()
        for res in chain:  # one residue per number, the first (microheterogeneity), as io does
            kind = _kind(res)
            key = (res.seqid.num, res.seqid.icode.strip())
            if kind is not None and key not in seen:
                seen.add(key)
                kept.append((res, kind, _heavy(res)))
        for n, (res, kind, atoms) in enumerate(kept):
            prev = kept[n - 1][2] if n else None
            nxt = kept[n + 1][2] if n + 1 < len(kept) else None
            r = len(labels)
            labels.append(ResLabel(chain.name, res.seqid.num, res.seqid.icode.strip(), res.name))
            for name, site in _roles(_table(kind, res.name, atoms, prev, nxt), atoms).items():
                found = [atoms[a].tolist() for a in site.ante if a in atoms]
                ante = np.mean(found, axis=0) if site.ante and len(found) == len(site.ante) else [np.nan] * 3
                angle = ROTATING_ANGLE if site.rotating else SP2_ANGLE
                rows.append((atoms[name].tolist(), ante, np.cos(np.radians(angle)), site.role, r, c, kind, name))
    return _Sites(
        np.array([x[0] for x in rows], float).reshape(-1, 3),
        np.array([x[1] for x in rows], float).reshape(-1, 3),
        np.array([x[2] for x in rows], float),
        np.array([x[3] in (D, DA) for x in rows], bool),
        np.array([x[3] in (A, DA) for x in rows], bool),
        np.array([x[4] for x in rows], int),
        np.array([x[5] for x in rows], int),
        np.array([x[6] for x in rows], int),
        labels,
        [x[7] for x in rows],
        chain_names,
    )


def _donates(sites: _Sites, d, a):
    """Atoms d donate to atoms a: d a donor, a an acceptor, and the antecedent-donor-acceptor angle wide enough
    for d (true when its antecedent is unknown)."""
    to_ante = sites.ante[d] - sites.xyz[d]
    to_acc = sites.xyz[a] - sites.xyz[d]
    with np.errstate(divide="ignore", invalid="ignore"):
        cos = np.einsum("ij,ij->i", to_ante, to_acc) / (
            np.linalg.norm(to_ante, axis=1) * np.linalg.norm(to_acc, axis=1)
        )
        angle_ok = np.isnan(cos) | (cos <= sites.max_cos[d] + 1e-12)
    return sites.donor[d] & sites.acceptor[a] & angle_ok


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
    sites = _sites(model, chains)
    if len(sites.xyz) < 2:
        return []
    ij = cKDTree(sites.xyz).query_pairs(MAX_DIST, output_type="ndarray")
    if not len(ij):
        return []
    i, j = ij[:, 0], ij[:, 1]
    if sel.kind == "na":
        keep = sites.kind[i] != sites.kind[j]
    elif sel.kind == "interchain":
        keep = sites.chain[i] != sites.chain[j]
    else:
        allowed = np.zeros((len(sites.chains),) * 2, bool)
        for a, b in sel.pairs:
            allowed[sites.chains.index(a), sites.chains.index(b)] = True
            allowed[sites.chains.index(b), sites.chains.index(a)] = True
        keep = allowed[sites.chain[i], sites.chain[j]]
    keep &= sites.owner[i] != sites.owner[j]
    i, j = i[keep], j[keep]
    dist = np.linalg.norm(sites.xyz[i] - sites.xyz[j], axis=1)
    keep = dist >= MIN_DIST
    i, j, dist = i[keep], j[keep], dist[keep]
    lo, hi = np.minimum(i, j), np.maximum(i, j)  # model order
    forward = _donates(sites, lo, hi)
    backward = ~forward & _donates(sites, hi, lo)
    found = forward | backward
    labels, names = sites.labels, sites.names
    out = [
        HBond(Atom(labels[sites.owner[d]], names[d]), Atom(labels[sites.owner[a]], names[a]), float(x))
        for d, a, x in zip(np.where(forward, lo, hi)[found], np.where(forward, hi, lo)[found], dist[found])
    ]
    return sorted(out)
