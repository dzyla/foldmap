"""Plain data containers shared by all pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np


class ResLabel(NamedTuple):
    chain: str
    seq: int
    icode: str
    name: str


class Bridge(NamedTuple):
    """Backbone H-bond bridge between residues i < j. kind: 'P' parallel, 'A' antiparallel."""

    i: int
    j: int
    kind: str


@dataclass(eq=False)
class Backbone:
    """Protein backbone of model 1. xyz[i] is rows N, CA, C, O of residue i.

    prev[i] is True when residue i is peptide-bonded to residue i-1.
    """

    labels: list[ResLabel]
    xyz: np.ndarray  # (N, 4, 3)
    prev: np.ndarray  # (N,) bool
    asu_chains: set[str] = field(default_factory=set)  # chains as deposited (the rest were built by symmetry)
    b: np.ndarray | None = None  # (N,) CA B-factor (or whatever the model stores there: local resolution, ADP)

    def __len__(self) -> int:
        return len(self.labels)

    def only(self, chains) -> Backbone:
        """The same backbone restricted to some chains."""
        keep = [k for k, l in enumerate(self.labels) if l.chain in set(chains)]
        return Backbone(
            [self.labels[k] for k in keep],
            self.xyz[keep],
            self.prev[keep],
            set(self.asu_chains),
            self.b[keep] if self.b is not None else None,
        )

    @property
    def nxt(self) -> np.ndarray:
        out = np.zeros(len(self), bool)
        out[:-1] = self.prev[1:]
        return out

    @property
    def ca(self) -> np.ndarray:
        return self.xyz[:, 1, :]


@dataclass(eq=False)
class SSE:
    """A helix ('H') or strand ('E'). start/end are inclusive residue indices."""

    kind: str
    chain: str
    start: int
    end: int
    first: ResLabel
    last: ResLabel
    centroid: np.ndarray
    axis: np.ndarray  # unit vector, points N-terminal -> C-terminal end

    @property
    def id(self) -> str:
        return f"{self.chain}:{self.first.seq}-{self.last.seq}"

    def __len__(self) -> int:
        return self.end - self.start + 1


@dataclass(eq=False)
class Sheet:
    """Strands in sheet order. directions[k] is +1 if strands[k] runs the same way as
    strands[0], else -1. pair_kinds[k] relates strands[k] and strands[k+1] ('P', 'A' or '?')."""

    strands: list[SSE]
    directions: list[int]
    pair_kinds: list[str]
    closed: bool = False  # barrel: a cycle was cut at its weakest pairing
    ambiguous: bool = False  # branched: ordering is a heuristic


@dataclass(eq=False)
class NAStrand:
    """One continuous nucleic-acid strand, 5' -> 3'. p: backbone point per nucleotide (P, else C4'), c1: C1'."""

    chain: str
    labels: list[ResLabel]
    p: np.ndarray  # (n, 3)
    c1: np.ndarray  # (n, 3)
    rna: bool

    @property
    def id(self) -> str:
        return f"{self.chain}:{self.labels[0].seq}-{self.labels[-1].seq}"

    def __len__(self) -> int:
        return len(self.labels)


@dataclass(eq=False)
class Nucleic:
    """Nucleic acid of model 1. pairs: (strand a, nucleotide i, strand b, nucleotide j) with a <= b.
    contacts: protein residue index (Backbone order) -> {(strand, nucleotide)} within contact distance."""

    strands: list[NAStrand] = field(default_factory=list)
    pairs: list[tuple[int, int, int, int]] = field(default_factory=list)
    contacts: dict[int, set[tuple[int, int]]] = field(default_factory=dict)
    contact_chain: dict[int, str] = field(default_factory=dict)  # protein chain of each contacting residue


def duplex_groups(na: Nucleic) -> list[list[int]]:
    """Strands joined by base pairs, as sorted groups (a lone strand is its own group)."""
    group = list(range(len(na.strands)))

    def root(x: int) -> int:
        while group[x] != x:
            x = group[x]
        return x

    for a, _, b, _ in na.pairs:
        group[root(a)] = root(b)
    out: dict[int, list[int]] = {}
    for k in range(len(na.strands)):
        out.setdefault(root(k), []).append(k)
    return sorted(out.values())


def duplex_axis(na: Nucleic, strands: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """(centre, unit axis) of a duplex: principal axis of all its backbone points, signed so the first
    strand runs 5' -> 3' along it."""
    pts = np.vstack([na.strands[s].p for s in strands])
    centre = pts.mean(axis=0)
    axis = np.linalg.svd(pts - centre)[2][0]
    first = na.strands[strands[0]].p
    if (first[-1] - first[0]) @ axis < 0:
        axis = -axis
    return centre, axis


@dataclass
class Ligand:
    """A small molecule or ion bound to the protein. contacts: Backbone indices of the residues holding it
    (coordinating a metal, or within reach of a ligand); additive: a common buffer, cryoprotectant, detergent or
    counter-ion, hidden unless asked for."""

    name: str
    chain: str
    seq: int
    metal: bool
    symbol: str  # element symbol of a metal ion (Zn, Fe, K...), else the residue name
    contacts: list[int] = field(default_factory=list)
    additive: bool = False


@dataclass(eq=False)
class Links:
    """Covalent extras of model 1. disulfides: (residue i, residue j) in Backbone order;
    glycans: (residue the glycan hangs on, its sugars from the root outward)."""

    disulfides: list[tuple[int, int]] = field(default_factory=list)
    glycans: list[tuple[int, list[str]]] = field(default_factory=list)
    ligands: list[Ligand] = field(default_factory=list)
