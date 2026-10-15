# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Plain data containers shared by all pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import NamedTuple

import numpy as np

from .nucleic import BasePair, NAGap, NAStructure, secondary_structure


class ResLabel(NamedTuple):
    chain: str
    seq: int
    icode: str
    name: str

    @property
    def tag(self) -> str:
        """Residue number with its insertion code, as written in a structure file: 100, 100A, -7."""
        return f"{self.seq}{self.icode}"


def res_ref(chain: str, first: ResLabel, last: ResLabel) -> str:
    """A residue-range reference (see layout.resolve) that names exactly first..last of a chain."""
    return f"res:{chain}:{first.tag}-{last.tag}"


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
    linked: frozenset[int] = frozenset()  # residues whose N/CA/C/O are stand-ins (a chromophore fused into the chain)

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
            frozenset(new for new, old in enumerate(keep) if old in self.linked),
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
    kinked: bool = False  # a helix piece continuing the one before it past a 3-10/pi kink (ss.build_sses split them)

    @property
    def id(self) -> str:
        return f"{self.chain}:{self.first.tag}-{self.last.tag}"

    def __len__(self) -> int:
        return self.end - self.start + 1


def group_by_chain(sses, keep=None) -> dict[str, list[SSE]]:
    """Elements grouped by chain, chains in order of first appearance; `keep` (ids) filters the elements."""
    out: dict[str, list[SSE]] = {}
    for s in sses:
        if keep is None or s.id in keep:
            out.setdefault(s.chain, []).append(s)
    return out


def homolog_chains(bb: Backbone) -> dict[str, tuple[str, int]]:
    """Copies of one protein, however much of each is modelled: chain -> (the first chain of its group of copies,
    the offset that takes its residue numbers to that chain's). Two chains are copies when at least 20 residues
    sit at the same numbers and at least 90 % of those have the same name; the numbers may also be shifted by the
    difference of the chains' first numbers (1EYG numbers chain k from 1000 k + 1)."""
    residues: dict[str, dict[tuple[int, str], str]] = {}
    for l in bb.labels:
        residues.setdefault(l.chain, {})[(l.seq, l.icode)] = l.name

    def alike(a: str, b: str, shift: int) -> bool:  # b's residue n is a's residue n - shift
        ra, rb = residues[a], residues[b]
        common = [k for k in rb if (k[0] - shift, k[1]) in ra]
        same = sum(rb[k] == ra[(k[0] - shift, k[1])] for k in common)
        return len(common) >= 20 and same >= 0.9 * len(common)

    first = {c: next(iter(r))[0] for c, r in residues.items()}
    out: dict[str, tuple[str, int]] = {}
    for b in residues:
        tries = ((a, shift) for a in out for shift in dict.fromkeys((0, first[b] - first[a])))
        a, shift = next(((a, shift) for a, shift in tries if alike(a, b, shift)), (b, 0))
        out[b] = (b, 0) if a == b else (out[a][0], out[a][1] + shift)
    return out


# ---------------------------------------------------------------------------------------------- copies alike


def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Residues two ranges (first, last residue number, in one numbering) share; 0 or less when they do not meet."""
    return min(a[1], b[1]) - max(a[0], b[0]) + 1


def same_place(a: tuple[int, int], b: tuple[int, int]) -> bool:
    """Are two ranges one element's place: sharing at least half of the shorter one's residues?"""
    o = overlap(a, b)
    return o > 0 and 2 * o >= min(a[1] - a[0], b[1] - b[0]) + 1


def pair_places(xs: list[tuple[str, int, int]], ys: list[tuple[str, int, int]]) -> dict[int, int]:
    """One-to-one pairs (index in xs -> index in ys) of places (kind, first, last residue number, one numbering):
    the same kind at the same place (same_place), the largest overlaps paired first, then the earliest."""
    found = sorted(
        (-overlap(x[1:], y[1:]), i, j)
        for i, x in enumerate(xs)
        for j, y in enumerate(ys)
        if x[0] == y[0] and same_place(x[1:], y[1:])
    )
    out: dict[int, int] = {}
    taken: set[int] = set()
    for _, i, j in found:
        if i not in out and j not in taken:
            out[i] = j
            taken.add(j)
    return out


@dataclass
class CopyClasses:
    """Elements of copies of one protein that are one element (copy_classes)."""

    reference: dict[str, str]  # first chain of each group of copies (model.homolog_chains) -> its reference copy
    mate: dict[str, str]  # element id -> the element it is: the reference copy's, else the first copy holding it
    groups: dict[str, list[str]]  # first chain -> the copies holding elements, in chain order


def copy_classes(sses: list[SSE], copies: dict[str, tuple[str, int]]) -> CopyClasses:
    """Which elements of copies of one protein (`copies`, homolog_chains) are one element. Two copies' elements
    pair one-to-one (pair_places, residue numbers taken to the first copy's). The reference is the copy pairing
    most elements with the others (more elements, then the earlier chain, break ties): every element pairing
    one of its elements is that element. The others, in chain order, pair one-to-one with those found so far
    that the reference lacks, else start their own. Chains with no copy holding elements are left out."""
    chains = list(dict.fromkeys([*copies, *(s.chain for s in sses)]))
    mine = group_by_chain(sses)

    def places(c: str) -> list[tuple[str, int, int]]:
        shift = copies.get(c, (c, 0))[1]
        return [(s.kind, s.first.seq - shift, s.last.seq - shift) for s in mine[c]]

    groups: dict[str, list[str]] = {}
    for c in chains:
        if c in mine:
            groups.setdefault(copies.get(c, (c, 0))[0], []).append(c)
    groups = {root: members for root, members in groups.items() if len(members) > 1}
    at = {c: places(c) for members in groups.values() for c in members}
    out = CopyClasses({}, {}, groups)
    for root, members in groups.items():
        paired = {m: sum(len(pair_places(at[c], at[m])) for c in members if c != m) for m in members}
        ref = max(members, key=lambda m: (paired[m], len(at[m]), -chains.index(m)))
        out.reference[root] = ref
        out.mate.update({s.id: s.id for s in mine[ref]})
        extra: list[tuple[str, tuple[str, int, int]]] = []  # elements the reference lacks: (first holder's id, place)
        for c in (c for c in members if c != ref):
            found = pair_places(at[c], at[ref])
            out.mate.update({mine[c][i].id: mine[ref][j].id for i, j in found.items()})
            rest = [i for i in range(len(at[c])) if i not in found]
            known = pair_places([at[c][i] for i in rest], [p for _, p in extra])
            for k, i in enumerate(rest):
                out.mate[mine[c][i].id] = extra[known[k]][0] if k in known else mine[c][i].id
            extra += [(mine[c][i].id, at[c][i]) for k, i in enumerate(rest) if k not in known]
    return out


def generic_label(kind: str, n: int) -> str:
    """The label of a chain's n-th (from 0) element of a kind: strands A…Z then S27…, 3-10 helices η1…, α1…"""
    if kind == "E":
        return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[n] if n < 26 else f"S{n + 1}"
    return f"η{n + 1}" if kind == "G" else f"α{n + 1}"


def shared_labels(sses: list[SSE], copies: dict[str, tuple[str, int]], labels: dict[str, str]) -> dict[str, str]:
    """id -> label for the elements of copies of one protein (copy_classes) whose label (`labels`, id -> label)
    changes so that one element has one label in every copy: the reference copy's. An element the reference lacks
    keeps its label unless a copy before it (the reference among them) or its own chain uses that label; then it
    takes the first free label of its kind (generic_label). The same element in later copies takes that label too."""
    classes = copy_classes(sses, copies)
    mine = group_by_chain(sses)
    out: dict[str, str] = {}
    for root, ref in classes.reference.items():
        given = {s.id: labels[s.id] for s in mine[ref]}  # mate -> the label its copies take
        taken, ours = set(given.values()), set(given)
        for c in (c for c in classes.groups[root] if c != ref):
            new = {s.id: given[classes.mate[s.id]] for s in mine[c] if classes.mate[s.id] in ours}
            used = set(new.values())
            for s in (s for s in mine[c] if s.id not in new):
                key = classes.mate[s.id]
                label = given.get(key)
                if label is None or label in used:
                    label = labels[s.id]
                    n = 0
                    while label in used or label in taken:
                        label, n = generic_label(s.kind, n), n + 1
                    given.setdefault(key, label)
                new[s.id] = label
                used.add(label)
            taken |= used
            out.update({k: v for k, v in new.items() if v != labels[k]})
    return out


@dataclass(eq=False)
class Sheet:
    """Strands in sheet order. directions[k] is +1 if strands[k] runs the same way as
    strands[0], else -1. pair_kinds[k] relates strands[k] and strands[k+1] ('P', 'A' or '?')."""

    strands: list[SSE]
    directions: list[int]
    pair_kinds: list[str]
    closed: bool = False  # barrel: a cycle was cut at its weakest pairing
    ambiguous: bool = False  # branched: ordering is a heuristic
    # column of each strand side by side (sheets._columns): two strands end to end along one edge of a partner (both
    # pairing with it, on its same side; not necessarily one strand) share one; empty: each strand its own
    columns: list[int] = field(default_factory=list)


@dataclass(eq=False)
class NAStrand:
    """One continuous nucleic-acid strand, 5' -> 3'. p: backbone point per nucleotide (P, else C4'), c1: C1'."""

    chain: str
    labels: list[ResLabel]
    p: np.ndarray  # (n, 3)
    c1: np.ndarray  # (n, 3)
    rna: bool
    cut5: bool = False  # cropped, or a chain break: the molecule continues beyond its 5' end
    cut3: bool = False  # ... beyond its 3' end
    bases: str = ""  # the parent base letter per nucleotide (A C G U T I; N when it has none); "N" throughout if unset

    def __post_init__(self):
        if not self.bases:
            self.bases = "N" * len(self.labels)

    @property
    def id(self) -> str:
        return f"{self.chain}:{self.labels[0].tag}-{self.labels[-1].tag}"

    def __len__(self) -> int:
        return len(self.labels)


@dataclass(eq=False)
class Nucleic:
    """Nucleic acid of model 1. pairs: (strand a, nucleotide i, strand b, nucleotide j) with a <= b.
    contacts: protein residue index (Backbone order) -> {(strand, nucleotide)} within contact distance.
    base_pairs: the same pairs with their kind; structure: stems, runs and ties built from them; gaps: where
    cropping left nucleotides out between two pieces of one strand; tetrads: G-tetrads of a G-quadruplex, (strand,
    nucleotide) x 4 (read, not drawn: the scope note says so)."""

    strands: list[NAStrand] = field(default_factory=list)
    pairs: list[tuple[int, int, int, int]] = field(default_factory=list)
    contacts: dict[int, set[tuple[int, int]]] = field(default_factory=dict)
    contact_chain: dict[int, str] = field(default_factory=dict)  # protein chain of each contacting residue
    base_pairs: list[BasePair] = field(default_factory=list)
    structure: NAStructure = field(default_factory=NAStructure)
    gaps: list[NAGap] = field(default_factory=list)
    tetrads: list[tuple[tuple[int, int], ...]] = field(default_factory=list)  # G-tetrads (nucleic.g_tetrads), not drawn

    @property
    def size(self) -> int:
        return sum(len(s) for s in self.strands)


def nucleic_from_pairs(strands, base_pairs, contacts, contact_chain, gaps=(), tetrads=()) -> Nucleic:
    """A Nucleic with its derived fields (pairs, structure) filled in."""
    base_pairs = sorted(base_pairs)
    return Nucleic(
        list(strands),
        [tuple(p[:4]) for p in base_pairs],
        contacts,
        contact_chain,
        base_pairs,
        secondary_structure(strands, base_pairs),
        list(gaps),
        sorted(tetrads),
    )


AUTO_NT = 300  # na_extent auto: larger nucleic acids are cropped to the stems the protein touches
_MERGE_GAP = 12  # crop by stems: kept stretches of one strand closer than this many nucleotides are joined


def _touched(na: Nucleic) -> set[tuple[int, int]]:
    return {x for v in na.contacts.values() for x in v}


def crop_nucleic(na: Nucleic, flank: int = 6, *, by: str = "contacts") -> Nucleic:
    """Only what the protein touches. by="contacts": the touched stretch of each strand plus `flank`
    nucleotides each side (and their base-pair partners, so duplexes stay whole); strands never touched are
    dropped. by="stems": the stems with a touched nucleotide plus one stem either side in sequence, and the
    short runs between them; each strand may fall into several pieces, with a gap recorded between them. A stem
    or isolated pair only partly inside goes into the gap whole, so no half stem and no paired nucleotide drawn
    as single-stranded; a piece left holding no stem nor an end of an isolated pair goes into the gap too, as does
    a strand of one nucleotide (ignored when drawn)
    (unless no piece holds a stem), its touchers then cropped like the rest."""
    if by == "contacts":
        window: dict[int, list[int]] = {}
        for k, i in _touched(na):
            lo, hi = window.get(k, [i, i])
            window[k] = [min(lo, i), max(hi, i)]
        for k, w in window.items():
            window[k] = [max(0, w[0] - flank), min(len(na.strands[k]) - 1, w[1] + flank)]
        for _ in range(3):  # partners of every kept nucleotide are kept too
            for a, i, b, j in na.pairs:
                for x, y, u, v in ((a, i, b, j), (b, j, a, i)):
                    if x in window and window[x][0] <= y <= window[x][1]:
                        lo, hi = window.get(u, [v, v])
                        window[u] = [min(lo, v), max(hi, v)]
        return _cut(na, {k: [tuple(w)] for k, w in window.items()})
    if by != "stems":
        raise ValueError(f"crop by contacts or stems; got {by!r}")
    stems = na.structure.stems
    touched = _touched(na)
    hit = [n for n, st in enumerate(stems) if touched & set(st.nucleotides())]
    # with the stem before and after each on its own strand (stems are in 5'->3' order of their first arm)
    keep = {m for n in hit for m in (n - 1, n, n + 1) if 0 <= m < len(stems) and stems[m].a == stems[n].a}
    spans: dict[int, list[list[int]]] = {}
    for n in sorted(keep):
        for s, lo, hi in (stems[n].arm5, stems[n].arm3):
            spans.setdefault(s, []).append([lo, hi])
    tie_of = {x: p for p in na.structure.ties for x in ((p.a, p.i), (p.b, p.j))}
    for k, i in sorted(touched):  # touched single-stranded nucleotides stay too; a touched tie with its partner
        tie = tie_of.get((k, i))
        for s, j in ((tie.a, tie.i), (tie.b, tie.j)) if tie else ((k, i),):
            spans.setdefault(s, []).append([j, j])
    windows: dict[int, list[tuple[int, int]]] = {}
    for k, found in spans.items():
        merged: list[list[int]] = []
        for lo, hi in sorted(found):
            if merged and lo - merged[-1][1] <= _MERGE_GAP:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi])
        windows[k] = [(lo, hi) for lo, hi in merged]
    inside = {(k, i) for k, w in windows.items() for lo, hi in w for i in range(lo, hi + 1)}
    for group in [st.nucleotides() for st in stems] + [[(p.a, p.i), (p.b, p.j)] for p in na.structure.ties]:
        if not all(x in inside for x in group):  # partly inside: its part inside goes into the gap
            inside -= set(group)
    pieces = {k: _stretches(sorted(i for s, i in inside if s == k)) for k in sorted({s for s, _ in inside})}
    held = {x for st in stems for x in st.nucleotides() if x in inside}  # stems, and isolated pairs kept whole
    held |= {x for p in na.structure.ties for x in ((p.a, p.i), (p.b, p.j)) if x in inside}
    stemmed = {
        k: [(lo, hi) for lo, hi in w if any((k, i) in held for i in range(lo, hi + 1)) and len(na.strands[k]) > 1]
        for k, w in pieces.items()
    }
    if any(stemmed.values()):  # a piece holding no stem (a touched single nucleotide far from them) goes into a gap
        pieces = {k: w for k, w in stemmed.items() if w}
    return _cut(na, pieces)


def _stretches(idx: list[int]) -> list[tuple[int, int]]:
    """Sorted indices as runs of consecutive ones, (first, last)."""
    out: list[tuple[int, int]] = []
    for i in idx:
        if out and i == out[-1][1] + 1:
            out[-1] = (out[-1][0], i)
        else:
            out.append((i, i))
    return out


def nucleic_touching(na: Nucleic, chains) -> Nucleic:
    """Only the strands a protein chain in `chains` touches, with every strand base-paired to them."""
    chosen = set(chains)
    contacts = {r: v for r, v in na.contacts.items() if na.contact_chain.get(r) in chosen}  # theirs only
    keep = {s for v in contacts.values() for s, _ in v}
    for _ in range(len(na.strands)):
        grown = keep | {p.b for p in na.base_pairs if p.a in keep} | {p.a for p in na.base_pairs if p.b in keep}
        if grown == keep:
            break
        keep = grown
    held = replace(na, contacts=contacts, contact_chain={r: na.contact_chain[r] for r in contacts})
    return _cut(held, {s: [(0, len(na.strands[s]) - 1)] for s in sorted(keep)})


def nucleic_of_chains(na: Nucleic, chains) -> Nucleic:
    """Only the strands of the nucleic-acid chains named (a chain's pieces and the gaps between them kept)."""
    chosen = set(chains)
    return _cut(na, {s: [(0, len(st) - 1)] for s, st in enumerate(na.strands) if st.chain in chosen})


def _cut(na: Nucleic, windows: dict[int, list[tuple[int, int]]]) -> Nucleic:
    """The nucleotides inside the windows, each window its own strand piece, indices renumbered. A gap joins
    consecutive pieces of one molecule (one strand, or strands joined by na.gaps at a chain break), counting every
    nucleotide left out between them; a piece's end is cut where the molecule goes on beyond it."""
    after = {g.before: g for g in na.gaps}  # strand k's 3' end goes on, after g.count nucleotides, in strand k + 1
    joined = {g.after for g in na.gaps}
    strands, where, gaps = [], {}, []
    last, left_out = None, 0  # the molecule's last piece kept so far, and the nucleotides left out since
    for k, s in enumerate(na.strands):
        if k not in joined:  # a new molecule
            last, left_out = None, 0
        pos = 0
        for lo, hi in windows.get(k, []):
            left_out += lo - pos
            if last is not None:
                gaps.append(NAGap(last, len(strands), left_out))
            for i in range(lo, hi + 1):
                where[(k, i)] = (len(strands), i - lo)
            strands.append(
                NAStrand(
                    s.chain,
                    s.labels[lo : hi + 1],
                    s.p[lo : hi + 1],
                    s.c1[lo : hi + 1],
                    s.rna,
                    s.cut5 or lo > 0 or k in joined,
                    s.cut3 or hi < len(s) - 1 or k in after,
                    s.bases[lo : hi + 1],
                )
            )
            last, left_out, pos = len(strands) - 1, 0, hi + 1
        left_out += len(s) - pos + (after[k].count if k in after else 0)
    pairs = []
    for p in na.base_pairs:
        if (p.a, p.i) in where and (p.b, p.j) in where:
            (a, i), (b, j) = sorted((where[(p.a, p.i)], where[(p.b, p.j)]))
            pairs.append(BasePair(a, i, b, j, p.kind))
    contacts = {r: {where[x] for x in v if x in where} for r, v in na.contacts.items()}
    contacts = {r: v for r, v in contacts.items() if v}
    chains = {r: c for r, c in na.contact_chain.items() if r in contacts}
    tetrads = [tuple(sorted(where[x] for x in t)) for t in na.tetrads if all(x in where for x in t)]
    return nucleic_from_pairs(strands, pairs, contacts, chains, gaps, tetrads)


@dataclass
class Ligand:
    """A small molecule or ion bound to the protein, or a small molecule bound to the nucleic acid alone. contacts:
    Backbone indices of the residues holding it (coordinating a metal, or within reach of a ligand); nucleotides:
    for one touching no protein residue, the nucleotides holding it; additive: a common buffer, cryoprotectant,
    detergent or counter-ion, hidden unless asked for."""

    name: str
    chain: str
    seq: int
    metal: bool
    symbol: str  # element symbol of a metal ion (Zn, Fe, K...), else the residue name
    contacts: list[int] = field(default_factory=list)
    additive: bool = False
    members: list[tuple[str, int]] = field(default_factory=list)  # (chain, residue number) of each part
    codes: list[str] = field(default_factory=list)  # residue codes of a covalently joined ligand, in bond order
    # (chain, number, icode) of the nucleotides holding a ligand that touches no protein residue (contacts empty)
    nucleotides: list[tuple[str, int, str]] = field(default_factory=list)


@dataclass(eq=False)
class Links:
    """Covalent extras of model 1. disulfides: (residue i, residue j) in Backbone order;
    glycans: (residue the glycan hangs on, its sugars from the root outward)."""

    disulfides: list[tuple[int, int]] = field(default_factory=list)
    glycans: list[tuple[int, list[str]]] = field(default_factory=list)
    ligands: list[Ligand] = field(default_factory=list)
    glycan_members: list[list[tuple[str, int]]] = field(default_factory=list)  # sugar residues of each glycan
    # metal ions held by nucleotides alone, by name: not marked (a ribosome holds hundreds), named by the scope note
    na_ions: dict[str, int] = field(default_factory=dict)
