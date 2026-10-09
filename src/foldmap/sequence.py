"""Sequence view data: a chain's sequence with its secondary structure (labelled as in the topology figure),
and multiple sequence alignments mapped onto it."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

from .dssp import assign_dssp
from .io import load_backbone
from .layout import provisional
from .ss import build_sses

GAP = "-"
SIMILAR_GROUPS = ("STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW")  # Clustal's strong groups
AMINO = "ACDEFGHIKLMNPQRSTVWY"


@dataclass
class Element:
    kind: str  # H helix, E strand, G 3-10 helix
    label: str  # as in the topology figure: α1, A, η1
    start: int  # first and last position in the track
    end: int
    id: str = ""  # the element's id in the topology figure


@dataclass
class Track:
    """One chain: its sequence (deposited, so unmodelled residues are included), residue numbers, which
    positions were modelled, a per-position secondary structure code (H E G T, '-' coil, ' ' unmodelled) and
    the elements."""

    chain: str
    letters: str
    numbers: list[int]
    observed: list[bool]
    ss: str
    elements: list[Element] = field(default_factory=list)
    residue: list[int | None] = field(default_factory=list)  # backbone index of each position (None: not modelled)


def _one_letter(names: list[str]) -> str:
    return "".join(
        (gemmi.find_tabulated_residue(n).one_letter_code.upper() or "X").replace(" ", "X")
        if gemmi.find_tabulated_residue(n)
        else "X"
        for n in names
    )


def chain_track(path, chain: str | None = None, full_sequence: bool = False) -> Track:
    """The chain's track. By default it covers the modelled span (first to last observed residue, with
    unmodelled stretches inside it); full_sequence covers the whole deposited sequence."""
    bb = load_backbone(path, "asu")
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, short_helices=True)
    labels = {k: p.label for k, p in provisional(sses).placed.items()}
    chains = list(dict.fromkeys(l.chain for l in bb.labels))
    if not chains:
        raise ValueError(f"{path}: no protein chain found")
    chain = chain or chains[0]
    if chain not in chains:
        raise ValueError(f"{path}: no protein chain {chain!r}; the file has {', '.join(chains)}")
    mine = [k for k, l in enumerate(bb.labels) if l.chain == chain]

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    polymer = st[0][chain].get_polymer()
    entity = st.get_entity_of(polymer) if len(polymer) else None
    full = list(entity.full_sequence) if entity is not None else []
    by_id = {(r.seqid.num, r.seqid.icode.strip()): r.label_seq for r in polymer}
    seats = [by_id.get((bb.labels[k].seq, bb.labels[k].icode)) for k in mine]
    if full and all(s is not None and 1 <= s <= len(full) for s in seats):
        names = [gemmi.Entity.first_mon(m) for m in full]
        position = {k: s - 1 for k, s in zip(mine, seats)}
    else:  # no deposited sequence (e.g. a PDB file without SEQRES): the modelled residues only
        names = [bb.labels[k].name for k in mine]
        position = {k: n for n, k in enumerate(mine)}
    n = len(names)
    residue: list[int | None] = [None] * n
    for k, p in position.items():
        residue[p] = k
    observed = [r is not None for r in residue]

    code = {"H": "H", "G": "G", "E": "E", "T": "T"}
    ss = [code.get(dssp.ss[r], "-") if r is not None else " " for r in residue]
    elements = []
    for s in sses:
        if s.chain != chain:
            continue
        a, b = position[s.start], position[s.end]
        for p in range(a, b + 1):
            ss[p] = s.kind
        elements.append(Element(s.kind, labels[s.id], a, b, s.id))

    numbers = [bb.labels[r].seq if r is not None else 0 for r in residue]
    seen = [p for p in range(n) if observed[p]]
    for p in range(n):  # unmodelled positions numbered from their nearest modelled neighbour
        if not observed[p]:
            q = min(seen, key=lambda q: abs(q - p))
            numbers[p] = numbers[q] + (p - q)

    lo, hi = (0, n - 1) if full_sequence else (seen[0], seen[-1])
    shift = lambda e: Element(e.kind, e.label, e.start - lo, e.end - lo, e.id)  # noqa: E731
    return Track(
        chain,
        _one_letter(names)[lo : hi + 1],
        numbers[lo : hi + 1],
        observed[lo : hi + 1],
        "".join(ss[lo : hi + 1]),
        [shift(e) for e in elements],
        residue[lo : hi + 1],
    )


@dataclass
class Alignment:
    names: list[str]
    rows: list[str]

    @property
    def width(self) -> int:
        return len(self.rows[0]) if self.rows else 0


def read_alignment(path) -> Alignment:
    """An alignment from FASTA (aligned, '-' gaps), Clustal or Stockholm; '.' gaps become '-'."""
    text = Path(path).read_text()
    lines = [l.rstrip() for l in text.splitlines()]
    first = next((l for l in lines if l.strip()), "")
    order: list[str] = []
    seqs: dict[str, str] = {}

    def add(name: str, chunk: str) -> None:
        if name not in seqs:
            order.append(name)
            seqs[name] = ""
        seqs[name] += chunk

    if first.startswith(">"):
        name = None
        for l in lines:
            if l.startswith(">"):
                name = l[1:].split()[0] if l[1:].split() else f"seq{len(order) + 1}"
                add(name, "")
            elif name is not None:
                add(name, l.strip())
    elif first.startswith("# STOCKHOLM"):
        for l in lines:
            if l and not l.startswith("#") and l != "//":
                parts = l.split()
                if len(parts) >= 2:
                    add(parts[0], parts[1])
    elif first.upper().startswith(("CLUSTAL", "MUSCLE", "PROBCONS")):
        for l in lines[1:]:
            if l and not l[0].isspace():
                parts = l.split()
                if len(parts) >= 2:
                    add(parts[0], parts[1])
    else:
        raise ValueError(f"{path}: unknown alignment format; use FASTA, Clustal or Stockholm")
    rows = ["".join(GAP if c in ".-~" else c.upper() for c in seqs[n]) for n in order]
    if len(rows) < 1 or len({len(r) for r in rows}) != 1:
        raise ValueError(f"{path}: alignment rows must all have the same length (got {sorted({len(r) for r in rows})})")
    return Alignment(order, rows)


def _global(a: str, b: str, match: float = 2.0, mismatch: float = -1.0, gap: float = 2.0):
    """Needleman-Wunsch (linear gaps): the aligned index pairs (i in a, j in b)."""
    n, m = len(a), len(b)
    A, B = np.frombuffer(a.encode(), np.uint8), np.frombuffer(b.encode(), np.uint8)
    H = np.zeros((n + 1, m + 1))
    H[0] = -gap * np.arange(m + 1)
    H[:, 0] = -gap * np.arange(n + 1)
    js = np.arange(m + 1)
    for i in range(1, n + 1):
        diag = H[i - 1, :-1] + np.where(A[i - 1] == B, match, mismatch)
        up = H[i - 1, 1:] - gap
        t = np.concatenate([[H[i, 0]], np.maximum(diag, up)])
        H[i] = np.maximum.accumulate(t + gap * js) - gap * js  # best of: here, or a run of gaps from the left
    pairs, i, j = [], n, m
    while i > 0 and j > 0:
        s = match if a[i - 1] == b[j - 1] else mismatch
        if np.isclose(H[i, j], H[i - 1, j - 1] + s):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif np.isclose(H[i, j], H[i - 1, j] - gap):
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


@dataclass
class Mapping:
    reference: int  # row of the alignment that is the structure's chain
    position: list[int | None]  # per alignment column: the track position it holds (None: none)


def map_alignment(aln: Alignment, track: Track, reference: str | None = None) -> Mapping:
    """Which alignment row is the structure (by name, or the row most identical to the chain), and where each
    column falls on the chain's track."""
    if reference is not None:
        hits = [k for k, n in enumerate(aln.names) if reference.lower() in n.lower()]
        if not hits:
            raise ValueError(f"no alignment row named {reference!r}; rows: {', '.join(aln.names)}")
        best = hits[0]
    else:

        def identity(row: str) -> float:
            plain = row.replace(GAP, "")
            pairs = _global(plain, track.letters)
            return sum(plain[i] == track.letters[j] for i, j in pairs) / max(len(track.letters), 1)

        best = max(range(len(aln.rows)), key=lambda k: identity(aln.rows[k]))
    row = aln.rows[best]
    plain = row.replace(GAP, "")
    to_track = dict(_global(plain, track.letters))
    position, k = [], 0
    for c in row:
        if c == GAP:
            position.append(None)
        else:
            position.append(to_track.get(k))
            k += 1
    return Mapping(best, position)


def column_classes(aln: Alignment, threshold: float = 0.7) -> tuple[np.ndarray, np.ndarray]:
    """(strict, similar) per column. Strict: every sequence has the same residue. Similar: not strict, but at
    least `threshold` of all sequences share one residue or one similarity group (gaps count against)."""
    cols = np.array([list(r) for r in aln.rows]).T if aln.rows else np.empty((0, 0))
    strict = np.zeros(aln.width, bool)
    similar = np.zeros(aln.width, bool)
    n = len(aln.rows)
    for c, col in enumerate(cols):
        letters = [x for x in col if x != GAP]
        if len(letters) == n and len(set(letters)) == 1:
            strict[c] = True
            continue
        if not letters:
            continue
        top = max(letters.count(x) for x in set(letters))
        group = max(sum(x in g for x in letters) for g in SIMILAR_GROUPS)
        similar[c] = max(top, group) / n >= threshold
    return strict, similar


def conservation(aln: Alignment) -> np.ndarray:
    """Per column, 0 (variable) .. 1 (invariant): 1 - Shannon entropy / log 20 over the residues present,
    scaled down by the fraction of gaps."""
    out = np.zeros(aln.width)
    n = len(aln.rows)
    for c in range(aln.width):
        letters = [r[c] for r in aln.rows if r[c] != GAP]
        if not letters:
            continue
        _, counts = np.unique(letters, return_counts=True)
        p = counts / counts.sum()
        h = float(-(p * np.log(p)).sum())
        out[c] = (1.0 - h / np.log(20)) * len(letters) / n
    return np.clip(out, 0.0, 1.0)


MIN_IDENTITY = 0.5  # a chain takes an alignment's scores only if it is at least this identical to the reference


def residue_conservation(source, aln: Alignment, reference: str | None = None) -> dict[int, float]:
    """Conservation (0 variable .. 1 invariant) per backbone residue, for every chain the alignment describes
    (homo-oligomers and assembly copies: all of them); chains unrelated to it get nothing. source: a
    structure file (its deposited chains) or a Backbone."""
    from .io import load_backbone
    from .model import Backbone

    bb = source if isinstance(source, Backbone) else load_backbone(source, "asu")
    scores = conservation(aln)
    if reference is not None:
        hits = [k for k, n in enumerate(aln.names) if reference.lower() in n.lower()]
        if not hits:
            raise ValueError(f"no alignment row named {reference!r}; rows: {', '.join(aln.names)}")
        rows = hits[:1]
    else:
        rows = list(range(len(aln.rows)))
    out: dict[int, float] = {}
    for chain in dict.fromkeys(l.chain for l in bb.labels):
        mine = [k for k, l in enumerate(bb.labels) if l.chain == chain]
        seq = _one_letter([bb.labels[k].name for k in mine])
        best = None
        for r in rows:  # the alignment row this chain matches best
            plain = aln.rows[r].replace(GAP, "")
            pairs = _global(plain, seq)
            same = sum(plain[i] == seq[j] for i, j in pairs)
            if best is None or same > best[0]:
                best = (same, r, dict(pairs))
        same, r, to_chain = best
        if same < MIN_IDENTITY * len(seq):
            continue
        k = 0
        for c, ch in enumerate(aln.rows[r]):
            if ch == GAP:
                continue
            if k in to_chain:
                out[mine[to_chain[k]]] = float(scores[c])
            k += 1
    return out
