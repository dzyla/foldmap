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
HBOND_STRETCH = 4.0  # Å, the one longer bond a pair stacked in register on a strict pair may have (1HYS, 1AOI)
C1_PAIR = (9.0, 12.0)  # Å, C1'-C1' across a Watson-Crick or wobble pair (10.5 +- 1.5)
PLANE_MAX = 45.0  # degrees between the two base planes (TBP and nucleosome DNA reach 40 of propeller twist)
MIN_LOOP = 3  # a hairpin closes at least this many unpaired nucleotides
SEGMENT_TURN = 25.0  # degrees the duplex axis may turn before a new straight piece starts
PIECE_MIN = 5  # base pairs in the shortest straight piece
SEGMENT_MIN = 2 * PIECE_MIN  # base pairs in the shortest stem that can split (into two shortest pieces)
AXIS_WINDOW = 11  # base pairs (about one helical turn) whose C1' atoms give the local duplex axis

# (purine, pyrimidine) parent bases -> the hydrogen bonds checked, how many must hold, the pair kind.
# 4-thiouridine (4SU) has S4 in place of O4, so by design it never makes the A-U pair's second bond.
_PAIRS = {
    ("G", "C"): ((("N1", "N3"), ("O6", "N4"), ("N2", "O2")), 2, "WC"),
    ("I", "C"): ((("N1", "N3"), ("O6", "N4")), 2, "WC"),
    ("A", "U"): ((("N1", "N3"), ("N6", "O4")), 2, "WC"),
    ("A", "T"): ((("N1", "N3"), ("N6", "O4")), 2, "WC"),
    ("G", "U"): ((("N1", "O2"), ("O6", "N3")), 2, "GU"),
    ("G", "T"): ((("N1", "O2"), ("O6", "N3")), 2, "GU"),
}
_PURINE = {"A", "G", "I"}
# Pseudouridine hangs from C5, not N1: its ring turns over, O2 takes uracil's O4 place and O4 takes O2's
_PSEUDO = {"O2": "O4", "O4": "O2"}
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


def _base_atom(res) -> gemmi.Atom | None:
    """The base atom bonded to C1' (N9 of a purine, N1 of a pyrimidine, C5 of pseudouridine), or None."""
    c1 = res.find_atom("C1'", "*")
    for name in ("N9", "N1", "C5"):
        a = res.find_atom(name, "*")
        if c1 is not None and a is not None and c1.pos.dist(a.pos) < 1.7:
            return a
    return None


def _glycosidic(res) -> np.ndarray | None:
    """C1' -> the base atom bonded to it."""
    a = _base_atom(res)
    if a is None:
        return None
    return np.array(a.pos.tolist()) - np.array(res.find_atom("C1'", "*").pos.tolist())


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

    A pair needs its hydrogen bonds within HBOND_MAX (two of the three for G-C), the C1'-C1' span of a paired
    duplex, roughly coplanar bases in the cis orientation of a double helix and, on one strand, partners more
    than MIN_LOOP apart (at least MIN_LOOP unpaired nucleotides between them). Each nucleotide keeps its best
    pair (more bonds, then shorter ones); only mutual best pairs are kept. A second pass, among the nucleotides
    still unpaired, takes a pair whose one bond is within HBOND_MAX and whose other is stretched to at most
    HBOND_STRETCH, but only where it stacks in register on a pair of the first pass (the next step of the same
    helix): stretched A·T and A·U pairs inside a duplex at 2-3 Å resolution, never an isolated one."""
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
    strict: list[tuple] = []  # (key, other key, score, kind) per candidate, both orders
    stretched: list[tuple] = []
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
        glycosidic = _base_atom(pyr[3])
        turned = glycosidic is not None and glycosidic.name == "C5"  # pseudouridine
        lengths = []
        for p_atom, q_atom in bonds:
            q_atom = _PSEUDO.get(q_atom, q_atom) if turned else q_atom
            p, q = pur[3].find_atom(p_atom, "*"), pyr[3].find_atom(q_atom, "*")
            if p is not None and q is not None:
                lengths.append(p.pos.dist(q.pos))
        lengths.sort()
        found = [d for d in lengths if d <= HBOND_MAX]
        tight = len(found) >= need
        loose = len(lengths) >= need and lengths[need - 2] <= HBOND_MAX and lengths[need - 1] <= HBOND_STRETCH
        if not tight and not loose:  # all but one bond within HBOND_MAX, that one within HBOND_STRETCH
            continue
        for k in (x, y):
            if k not in normals:
                normals[k] = _normal(bases[k][3])
        nx, ny = normals[x], normals[y]
        if nx is not None and ny is not None and abs(float(nx @ ny)) < np.cos(np.radians(PLANE_MAX)):
            continue
        if not _cis(bx[4], _glycosidic(bx[3]), by[4], _glycosidic(by[3])):
            continue
        held = found if tight else [d for d in lengths if d <= HBOND_STRETCH]
        score = (len(held), -float(np.mean(held)))
        for me, other in ((bx, by), (by, bx)):
            (strict if tight else stretched).append(((me[0], me[1]), (other[0], other[1]), score, kind))
    out = _mutual_best(strict)
    paired = {key for p in out for key in ((p.a, p.i), (p.b, p.j))}
    steps = {(p.a, p.i, p.b, p.j) for p in out}

    def in_register(key, other) -> bool:  # stacks on a strict pair: the helix's next step either way
        (a, i), (b, j) = sorted((key, other))
        return (a, i - 1, b, j + 1) in steps or (a, i + 1, b, j - 1) in steps

    out |= _mutual_best([c for c in stretched if not {c[0], c[1]} & paired and in_register(c[0], c[1])])
    return sorted(out)


def _mutual_best(candidates: list[tuple]) -> set[BasePair]:
    """The pairs of (key, other key, score, kind) candidates (both orders listed) that are each nucleotide's best
    (highest score; the first listed on a tie) for both of its nucleotides."""
    best: dict[tuple[int, int], tuple] = {}
    for key, other, score, kind in candidates:
        if key not in best or score > best[key][0]:
            best[key] = (score, other, kind)
    out = set()
    for key, (_, other, kind) in best.items():
        if best.get(other, (None, None))[1] == key:
            (a, i), (b, j) = sorted((key, other))
            out.add(BasePair(a, i, b, j, kind))
    return out


@dataclass(eq=False)
class Stem:
    """Stacked base pairs. pairs[0] is the outermost: the first arm (strand a) runs 5'->3' from it, the second
    arm (strand b) runs 5'->3' back to it. bulges: (strand, nucleotide) left unpaired inside the stem."""

    pairs: list[BasePair]
    bulges: list[tuple[int, int]] = field(default_factory=list)
    id: str = ""
    name: str = ""
    parent: str = ""  # a straight piece of a bent stem: the whole stem's id

    @property
    def a(self) -> int:
        return self.pairs[0].a

    @property
    def b(self) -> int:
        return self.pairs[0].b

    @property
    def arm5(self) -> tuple[int, int, int]:
        """(strand, first, last) nucleotide of the first arm."""
        return self.a, self.pairs[0].i, self.pairs[-1].i

    @property
    def arm3(self) -> tuple[int, int, int]:
        """(strand, first, last) nucleotide of the second arm, 5'->3'."""
        return self.b, self.pairs[-1].j, self.pairs[0].j

    def __len__(self) -> int:
        return len(self.pairs)

    def nucleotides(self) -> list[tuple[int, int]]:
        """Every (strand, nucleotide) of both arms, bulges included: first arm 5'->3', then the second."""
        (a, i0, i1), (b, j0, j1) = self.arm5, self.arm3
        return [(a, i) for i in range(i0, i1 + 1)] + [(b, j) for j in range(j0, j1 + 1)]


class NARun(NamedTuple):
    """Nucleotides lo..hi of one strand that are in no stem (a tie's nucleotides can be among them).
    left/right: the stem nucleotide (strand, nucleotide) just before and after it (None at a strand end).
    kind: hairpin, internal, junction, linker or tail."""

    strand: int
    lo: int
    hi: int
    left: tuple[int, int] | None
    right: tuple[int, int] | None
    kind: str


class NAGap(NamedTuple):
    """Nucleotides left out by cropping between strand piece `before` (its 3' end) and `after` (its 5' end)."""

    before: int
    after: int
    count: int


@dataclass(eq=False)
class NAStructure:
    stems: list[Stem] = field(default_factory=list)
    runs: list[NARun] = field(default_factory=list)
    ties: list[BasePair] = field(default_factory=list)  # isolated pairs: drawn as a bar, not a stem


def _stacks(p: BasePair, q: BasePair, paired: set[tuple[int, int]], strands=None) -> bool:
    """q continues p's stem: the next step, a bulge of up to 2 nt on one side, a 1x1 mismatch, or (given the
    strands) a 2x2 loop in register whose two unpaired steps still sit across the helix (their C1' atoms within
    C1_PAIR: non-canonical or stretched pairs of the same duplex, as in the nucleosome's 1AOI I20-21). A
    nucleotide the step skips must be unpaired (`paired`: every paired (strand, nucleotide)), else it would sit
    in the stem and in another pair."""
    if (p.a, p.b) != (q.a, q.b):
        return False
    gi, gj = q.i - p.i - 1, p.j - q.j - 1
    if gi < 0 or gj < 0:
        return False
    across = (
        gi == gj == 2
        and strands is not None
        and all(
            C1_PAIR[0] <= float(np.linalg.norm(strands[p.a].c1[p.i + k] - strands[p.b].c1[p.j - k])) <= C1_PAIR[1]
            for k in (1, 2)
        )
    )
    if not ((gi == 0 and gj <= 2) or (gj == 0 and gi <= 2) or (gi == gj == 1) or across):
        return False
    skipped = [(p.a, i) for i in range(p.i + 1, q.i)] + [(p.b, j) for j in range(q.j + 1, p.j)]
    return not any(n in paired for n in skipped)


def _group_stems(pairs: list[BasePair], strands=None) -> list[list[BasePair]]:
    paired = {(p.a, p.i) for p in pairs} | {(p.b, p.j) for p in pairs}
    stems: list[list[BasePair]] = []
    for p in sorted(pairs, key=lambda p: (p.a, p.b, p.i, -p.j)):
        for st in stems:
            if _stacks(st[-1], p, paired, strands):
                st.append(p)
                break
        else:
            stems.append([p])
    return stems


def _nested(groups: list[list[BasePair]]) -> list[BasePair]:
    """The pairs of the intra-strand stems (`groups`, each a stem's pairs, in 5'->3' order) that cross no other,
    keeping the pairs of the longest stems first."""
    keep: list[BasePair] = []
    for st in sorted([g for g in groups if g[0].a == g[0].b], key=len, reverse=True):
        if all(not (p.a == q.a and (p.i < q.i < p.j < q.j or q.i < p.i < q.j < p.j)) for p in st for q in keep):
            keep += st
    return keep


def _loop_kind(strand: int, lo: int, hi: int, nested: list[BasePair]) -> str:
    """hairpin / internal / junction by the loop the run lies in (the innermost pair enclosing it on its strand
    and how many stems branch off that loop); linker when nothing encloses it."""
    mine = [p for p in nested if p.a == strand and p.i < lo and hi < p.j]
    if not mine:
        return "linker"
    close = max(mine, key=lambda p: p.i)
    inner = [p for p in nested if p.a == strand and close.i < p.i and p.j < close.j]
    branches = [p for p in inner if not any(q.i < p.i and p.j < q.j for q in inner)]
    return "hairpin" if not branches else "internal" if len(branches) == 1 else "junction"


def secondary_structure(strands, pairs: list[BasePair]) -> NAStructure:
    """Stems (two or more stacked pairs), ties (an isolated pair) and the single-stranded runs between them.
    Stems are named P1, P2 ... (RNA) or DNA 1, DNA 2 ... in 5'->3' order of their first arm."""
    out = NAStructure()
    for group in _group_stems(pairs, strands):
        if len(group) < 2:
            out.ties += group
            continue
        stem = Stem(group)
        for p, q in zip(group, group[1:]):
            stem.bulges += [(p.a, i) for i in range(p.i + 1, q.i)] + [(p.b, j) for j in range(q.j + 1, p.j)]
        stem.bulges.sort()
        out.stems.append(stem)
    out.stems.sort(key=lambda s: (s.a, s.pairs[0].i))
    rna = dna = 0
    for st in out.stems:
        la, lb = strands[st.a].labels, strands[st.b].labels
        (_, i0, i1), (_, j0, j1) = st.arm5, st.arm3
        st.id = f"{la[i0].chain}:{la[i0].tag}-{la[i1].tag}/{lb[j0].chain}:{lb[j0].tag}-{lb[j1].tag}"
        if strands[st.a].rna or strands[st.b].rna:
            rna += 1
            st.name = f"P{rna}"
        else:
            dna += 1
            st.name = f"DNA {dna}"
    inside = {n for st in out.stems for n in st.nucleotides()}
    nested = _nested([st.pairs for st in out.stems])
    for s, strand in enumerate(strands):
        i, n = 0, len(strand)
        while i < n:
            if (s, i) in inside:
                i += 1
                continue
            j = i
            while j + 1 < n and (s, j + 1) not in inside:
                j += 1
            left = (s, i - 1) if i > 0 else None
            right = (s, j + 1) if j + 1 < n else None
            kind = "tail" if left is None or right is None else _loop_kind(s, i, j, nested)
            out.runs.append(NARun(s, i, j, left, right, kind))
            i = j + 1
    return out


def stem_axis(stem: Stem, strands) -> tuple[np.ndarray, np.ndarray]:
    """(centre, unit axis) of a stem (two or more pairs), pointing the way its first arm runs 5'->3'.

    A stem of a helical turn or more: the principal axis of its backbone points (longer than it is wide).
    A shorter one: the direction perpendicular to every pair's C1'-C1' line, since the pairs turn about the
    axis; this holds down to two pairs, where a principal axis would point across the stem instead. With only
    two or three pairs whose C1'-C1' lines are nearly parallel (little twist between them), that direction is
    poorly defined and the axis can tilt far from the true one."""
    a_pts = np.array([strands[p.a].c1[p.i] for p in stem.pairs])
    b_pts = np.array([strands[p.b].c1[p.j] for p in stem.pairs])
    centre = np.vstack([a_pts, b_pts]).mean(axis=0)
    if len(stem) >= AXIS_WINDOW:
        backbone = np.array([strands[s].p[i] for s, i in stem.nucleotides()])
        axis = np.linalg.svd(backbone - backbone.mean(axis=0), full_matrices=False)[2][0]
    else:
        axis = np.linalg.svd(b_pts - a_pts)[2][2]
    run = (a_pts[-1] + b_pts[-1]) / 2 - (a_pts[0] + b_pts[0]) / 2  # from the outermost pair inward
    if axis @ run < 0:
        axis = -axis
    return centre, axis


def _local_axes(stem: Stem, strands) -> np.ndarray:
    """Unit axis at each base pair: the direction perpendicular to the C1'-C1' lines of the pairs around it
    (about a helical turn of them, at least PIECE_MIN), pointing the way the first arm runs."""
    a_pts = np.array([strands[p.a].c1[p.i] for p in stem.pairs])
    b_pts = np.array([strands[p.b].c1[p.j] for p in stem.pairs])
    mids = (a_pts + b_pts) / 2
    n = len(stem)
    w = min(AXIS_WINDOW, max(PIECE_MIN, n // 2))
    out = np.zeros((n, 3))
    for k in range(n):
        lo = min(max(0, k - w // 2), n - w)
        axis = np.linalg.svd(b_pts[lo : lo + w] - a_pts[lo : lo + w])[2][2]
        out[k] = axis if (mids[lo + w - 1] - mids[lo]) @ axis >= 0 else -axis
    return out


def segment_stem(stem: Stem, strands, max_turn: float = SEGMENT_TURN) -> list[Stem]:
    """A long stem cut into straight pieces wherever its axis has turned more than max_turn degrees from the
    start of the current piece (pieces at least PIECE_MIN pairs). Stems under SEGMENT_MIN pairs stay whole.
    A cut falls only between two pairs that stack directly, so a bulge or mismatch stays inside a piece."""
    if len(stem) < SEGMENT_MIN:
        return [stem]
    axes = _local_axes(stem, strands)
    cuts, start = [0], 0
    for k in range(1, len(stem)):
        p, q = stem.pairs[k - 1], stem.pairs[k]
        if q.i != p.i + 1 or q.j != p.j - 1:  # a bulge or mismatch between them: the next plain step takes the cut
            continue
        turned = np.degrees(np.arccos(np.clip(axes[k] @ axes[start], -1.0, 1.0)))
        if turned > max_turn and k - start >= PIECE_MIN:
            cuts.append(k)
            start = k
    if len(cuts) > 1 and len(stem) - cuts[-1] < PIECE_MIN:  # a short last piece joins the one before
        cuts.pop()
    if len(cuts) == 1:
        return [stem]
    letters = "abcdefghijklmnopqrstuvwxyz"
    pieces = []
    for n, (lo, hi) in enumerate(zip(cuts, cuts[1:] + [len(stem)])):
        part = stem.pairs[lo:hi]
        inner = [x for x in stem.bulges if _between(x, part)]
        suffix = letters[n] if n < len(letters) else str(n + 1)
        pieces.append(Stem(part, inner, f"{stem.id}#{n + 1}", f"{stem.name}{suffix}", stem.id))
    return pieces


def _between(x: tuple[int, int], part: list[BasePair]) -> bool:
    s, i = x
    return (s == part[0].a and part[0].i < i < part[-1].i) or (s == part[0].b and part[-1].j < i < part[0].j)


def bend(pieces: list[Stem], strands) -> float:
    """Angle (degrees) between the axes of the first and the last piece."""
    a, b = stem_axis(pieces[0], strands)[1], stem_axis(pieces[-1], strands)[1]
    return float(np.degrees(np.arccos(np.clip(a @ b, -1.0, 1.0))))


def straight_duplex(structure: NAStructure, strands) -> list[int] | None:
    """The two strands of the one long straight duplex that dominates the nucleic acid, or None. A protein-DNA
    complex built around such a duplex keeps the classic picture: the duplex across the page, protein below."""
    if not structure.stems:
        return None
    biggest = max(structure.stems, key=len)
    total = sum(len(s) for s in structure.stems)
    if biggest.a == biggest.b or len(biggest) < 0.8 * total or len(segment_stem(biggest, strands)) > 1:
        return None
    return [biggest.a, biggest.b]


@dataclass(eq=False)
class StemView:
    """A stem as frame.view_frame sees an element: a centroid, an axis and a size."""

    centroid: np.ndarray
    axis: np.ndarray
    n: int

    def __len__(self) -> int:
        return self.n


def stem_views(structure: NAStructure, strands) -> list[StemView]:
    """Each straight stem piece, weighted by its nucleotides, for choosing the page axes."""
    return [
        StemView(*stem_axis(piece, strands), len(piece.nucleotides()))
        for st in structure.stems
        for piece in segment_stem(st, strands)
    ]
