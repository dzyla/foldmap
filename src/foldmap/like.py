# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A structure B drawn like a reference A (plan 13a, --like): what B's layout takes from A's finished one.

The reference enters B's layout at its frame (frame.frame_like: A's frame carried through the superposition) and at
its homes and choices (build_layout's `like`): each of B's elements with a mate in A (align.element_mates) starts at
its mate's place on A's page, held there harder than a projected home (layout.LIKE_WEIGHT), with its mate's angle,
label and, for identity colourings, colour; B's sheets are drawn as βα rows where their mates' are, and B's helices
over their sheet where their mates are. The optimiser itself is not changed: B's own pairing decides each sheet's
strand order (science S1-S3, S6), read left to right as A reads its mates (Like.order), and _separate clears every
overlap as always.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .align import ChainAlignment, Mates
from .layout import LIKE_WEIGHT, PIECE_GAP, Layout, can_be_row, helix_length, snap, template_choices
from .model import SSE, Sheet, generic_label, res_ref
from .style import Style

IDENTITY_COLOURINGS = ("chain", "sequence", "shade", "sstype", "domain")  # a mate takes its mate's colour
EXTRA_TINT = 0.35  # an extra element's colour: its N-terminal neighbour's, mixed this much with white
_LETTERS = "abcdefghijklmnopqrstuvwxyz"
_PRIMES = ("′", "″", "‴")
_STRAND_LETTER = re.compile(r"[A-Z][′″‴]*")  # a strand named by a capital letter (generic labels, Ig's A-G)
LOW_COVERAGE = 0.7  # aligned residues under this fraction of the shorter chain: the note says how few aligned
_METHODS = {"sequence": "sequence alignment", "structure": "structure alignment"}
_HELIX_SIGN = {"H": "α", "G": "η"}  # the generic label's sign of each helix type (model.generic_label)


@dataclass(frozen=True)
class Like:
    """What B's layout takes from the reference's finished layout (build_layout's `like`)."""

    target: dict[str, tuple[float, float]]  # B element id -> page centre of its mate in A's finished layout
    angle: dict[str, float]  # B element id -> its mate's drawn angle
    label: dict[str, str]  # B element ref (res:chain:a-b) -> label (mates, pieces, extras)
    colour: dict[str, str]  # B element id -> colour (identity colourings only; else empty)
    template: list[bool] | None  # A's βα-row choices, per sheet of B (its own for a sheet with no mate); None: B's own
    # A's helix-over-sheet choices, by B helix id: False takes a helix off its sheet; True keeps it over only where
    # B's own packing (features.helices_on_sheets) puts it on a sheet (one way: an over-sheet helix needs a sheet)
    over: dict[str, bool]
    weight: float  # _K_LIKE / _K_HOME (layout.LIKE_WEIGHT)
    # B element id -> (cx, cy) of the unit its mate was laid out in and the mate's offset there (Layout.units): a
    # unit laid out like its mate's lands exactly on that unit's centre (build_layout's homes, layout._like_homes).
    # Read only where it gives `target` (to 1e-9): a target edited without its unit is followed, on the grid
    unit: dict[str, tuple[float, float, float, float]] = field(default_factory=dict)
    # pairs of B element ids the reference relates as they are related here (their mates neighbours in sequence or
    # bridged by a disulfide there; parts of one helix): no disulfide spring between them when both are held at the
    # reference's places (layout._attract); a disulfide only B has still pulls (a loop between two held items never
    # does: the reference places them). None: every pair
    settled: frozenset[frozenset[str]] | None = None
    # B element id -> its mate's label, for a mate of the other helix type (an α-helix whose mate is a 3-10 helix, or
    # the reverse): named by its own type (_labels), the note says which mate's name it stands for (like_note)
    crossed: dict[str, str] = field(default_factory=dict)
    # per sheet of B (as `template`): True when A draws B's mated track strands left to right in B's track order,
    # False when right to left, None when it does not say (fewer than two mated, mates in two of A's sheet blocks,
    # not in one order there, or pointing otherwise than B's pairing has them: a pairing difference). None: B's own
    # (layout._sheet_order). A sheet seen nearly edge-on has its left-right order decided by a fraction of an
    # ångström; drawn like A it reads as A's does, never as A's mirror image (3WZ0-2AV5, plan 13a-12)
    order: list[bool | None] | None = None


def _family(s: SSE) -> str:
    return "E" if s.kind == "E" else "H"


def like_from(
    layout_a: Layout,
    sses_a: list[SSE],
    sses_b: list[SSE],
    mates: Mates,
    look: Style,
    sheets_b: list[Sheet] | None = None,
) -> Like:
    """B's Like from A's finished layout and the element mates (align.element_mates, with B's elements as `sses_b`
    and A's as `sses_a`, the elements A's layout drew).

    - Places and angles: a mate takes its mate's page centre and drawn angle. A helix B splits where A does not (the
      mate and its pieces, Mates.pieces) is laid end to end along its mate as drawn, in sequence, PIECE_GAP apart,
      each at its own drawn length and the mate's angle, centred on the mate (those places on the grid).
    - Labels: a mate takes its mate's label; of a split helix's parts the first in sequence takes it and the next
      ones the same primed (α3, α3′, α3″). A mate of the other helix type (an α-helix whose mate is the 3-10 helix
      η1) keeps its own type's sign with its mate's number (α1), so no name calls an α-helix 3-10; where the
      reference names another element so (its α1), it is named as an extra is (Like.crossed: the note says which
      mate's name it stands for). An extra is named after the matched element of its own type (strand, α- or
      3-10 helix) before it, lettered (α3a, η1a), the first name not taken; a strand after one named by a capital
      letter is primed instead (C′, C″: Ig's extra strands); an extra with none before it keeps its own label
      unless that is taken, then lettered (or primed) too.
    - Colours (identity colourings only, IDENTITY_COLOURINGS): a mate or piece takes its mate's colour as A draws it
      (render.element_colours), an extra a tint of its N-terminal neighbour's (the first coloured element after it
      when none comes before; coloured by type, sstype, the neighbour is one of its own kind).
    - Sheet order: `sheets_b` gives `order` (_order): each sheet's strands left to right as A draws their mates.
    - Choices: `sheets_b` (B's sheets, as build_layout gets them) gives `template`: each sheet a βα row when most of
      its mated strands stand in one of A's rows (Layout.ba_rows), not when most do not, else (a tie, or no strand
      mated) as B's own choice (layout.template_choices); never a sheet that cannot be a row. `over`: each B helix
      with a mate or piece drawn over its sheet exactly when the mate is in A (Layout.on_sheet).
    """
    placed = layout_a.placed
    by_b = {s.id: s for s in sses_b}
    groups: dict[str, list[str]] = {}  # A element -> its mate and pieces in B
    for b, a in [*mates.mate.items(), *mates.pieces.items()]:
        if a in placed and b in by_b:
            groups.setdefault(a, []).append(b)
    target: dict[str, tuple[float, float]] = {}
    angle: dict[str, float] = {}
    unit: dict[str, tuple[float, float, float, float]] = {}
    for a, bs in groups.items():
        p = placed[a]
        if len(bs) == 1:
            (b,) = bs
            target[b], angle[b] = (p.cx, p.cy), p.angle
            if a in layout_a.units:
                unit[b] = layout_a.units[a]
            continue
        bs = sorted(bs, key=lambda b: (by_b[b].start, b))
        lengths = [helix_length(by_b[b], look) for b in bs]
        along = -(sum(lengths) + PIECE_GAP * (len(bs) - 1)) / 2
        d = np.array([np.cos(p.angle), np.sin(p.angle)])
        for b, n in zip(bs, lengths):
            t = along + n / 2
            # on the grid: cos and sin differ in the last bit between CPUs (A12)
            target[b], angle[b] = tuple(float(v) for v in snap((p.cx + d[0] * t, p.cy + d[1] * t))), p.angle
            along += n + PIECE_GAP
    label, crossed = _labels(placed, sses_b, mates)
    return Like(
        target=target,
        angle=angle,
        label=label,
        colour=_colours(layout_a, sses_a, sses_b, mates, look),
        template=_template(layout_a, sses_b, mates, sheets_b) if sheets_b is not None else None,
        over={
            b: a in layout_a.on_sheet
            for b, a in [*mates.mate.items(), *mates.pieces.items()]
            if b in by_b and by_b[b].kind == "H"
        },
        weight=LIKE_WEIGHT,
        unit=unit,
        settled=_settled(layout_a, sses_a, groups),
        crossed=crossed,
        order=[_order(layout_a, sheet, mates) for sheet in sheets_b] if sheets_b is not None else None,
    )


def _order(layout_a: Layout, sheet: Sheet, mates: Mates) -> bool | None:
    """Does A draw this B sheet's mated track strands left to right in B's track order (True), right to left
    (False), or not say (None: fewer than two mated strands drawn in one of A's sheet blocks, their mates not in
    one order there, or pointing otherwise than B's pairing has them)? Like.order."""
    placed = layout_a.placed
    block = {k: n for n, ids in enumerate(layout_a.sheets) for k in ids}
    mated = [
        (k, mates.mate[s.id])
        for k, s in enumerate(sheet.strands)
        if s.id in mates.mate and mates.mate[s.id] in placed and mates.mate[s.id] in block
    ]
    if len(mated) < 2 or len({block[a] for _, a in mated}) != 1:
        return None
    # each strand drawn up or down as its mate, relative to the first: as B's own pairing has them (directions)
    ups = [(1 if np.sin(placed[a].angle) > 0 else -1) * sheet.directions[k] for k, a in mated]
    if len(set(ups)) != 1:
        return None
    xs = [placed[a].cx for _, a in mated]
    steps = [y - x for x, y in zip(xs, xs[1:])]
    if all(d > 0 for d in steps):
        return True
    if all(d < 0 for d in steps):
        return False
    return None


def _settled(layout_a: Layout, sses_a: list[SSE], groups: dict[str, list[str]]) -> frozenset[frozenset[str]]:
    """Like.settled: the pairs of B elements (mates and pieces, `groups`: A element -> them) whose mates the
    reference relates: one A element (a split helix's parts), A elements next to each other in sequence in one chain
    (a loop), or two joined by a disulfide (each cysteine's element, or the nearest element end, as cli's bridge
    springs place it)."""
    seq = sorted(sses_a, key=lambda s: s.start)
    related = {frozenset((a.id, b.id)) for a, b in zip(seq, seq[1:]) if a.chain == b.chain}
    links = getattr(layout_a, "links", None)
    for i, j in getattr(links, "disulfides", None) or []:
        ends = [_cys_element(layout_a, sses_a, r) for r in (i, j)]
        if None not in ends and ends[0] != ends[1]:
            related.add(frozenset(ends))
    out = set()
    for bs in groups.values():
        out |= {frozenset((x, y)) for x in bs for y in bs if x != y}
    for pair in related:
        x, y = sorted(pair)
        out |= {frozenset((p, q)) for p in groups.get(x, []) for q in groups.get(y, [])}
    return frozenset(out)


def _cys_element(layout_a: Layout, sses_a: list[SSE], r: int) -> str | None:
    """The element a residue (Backbone index of A's layout) is drawn on for a disulfide: its own, or the nearest."""
    if r >= len(layout_a.res_chain):
        return None
    mine = [s for s in sses_a if s.chain == layout_a.res_chain[r]]
    if not mine:
        return None
    return min(mine, key=lambda s: 0 if s.start <= r <= s.end else min(abs(r - s.start), abs(r - s.end))).id


def _labels(placed: dict, sses_b: list[SSE], mates: Mates) -> tuple[dict[str, str], dict[str, str]]:
    """(B element ref -> label, Like.crossed) (like_from: mates, pieces, extras)."""
    names: dict[str, str] = {}  # B element id -> label
    start = {s.id: s.start for s in sses_b}
    kind = {s.id: s.kind for s in sses_b}
    parts: dict[str, list[str]] = {}  # A element -> its mate and pieces in B
    for b, a in [*mates.mate.items(), *mates.pieces.items()]:
        if a in placed and b in start:
            parts.setdefault(a, []).append(b)
    crossed: dict[str, str] = {}  # B element id -> its mate's name, of a mate of the other helix type
    for a, bs in parts.items():  # the first part in sequence takes the name, the next ones primed (α6, α6′)
        for k, b in enumerate(sorted(bs, key=lambda b: (start[b], b))):
            label = placed[a].label
            label = label if k == 0 else label + (_PRIMES[k - 1] if k <= len(_PRIMES) else "′" * k)
            names[b] = _own_type(label, kind[b])
            if names[b] != label:
                crossed[b] = label
    reference = {p.label for p in placed.values()}
    unnamed = {b for b in crossed if names[b] in reference}  # the reference's name of another element: as an extra
    for b in unnamed:
        del names[b]
    taken = reference | set(names.values())
    extras = set(mates.extras) | unnamed
    last: dict[tuple[str, str], str] = {}  # (chain, type) -> the name of the last matched element so far
    counts: dict[tuple[str, str], int] = {}  # B's own generic labels, as layout._assign_labels gives them
    for s in sorted(sses_b, key=lambda s: (s.chain, s.start)):
        key = (s.chain, s.kind)
        n = counts.get(key, 0)
        counts[key] = n + 1
        if s.id in names and s.id in mates.mate:
            last[key] = _own_type(placed[mates.mate[s.id]].label, s.kind)  # its mate's name, unprimed
        elif s.id in extras:
            base = last.get(key)
            if base is None:
                own = generic_label(s.kind, n)
                if own not in taken:
                    taken.add(own)
                    continue  # its own label, free
                base = own
            extra = _primed if s.kind == "E" and _STRAND_LETTER.fullmatch(base) else _lettered
            k = 0
            while extra(base, k) in taken:
                k += 1
            names[s.id] = extra(base, k)
            taken.add(names[s.id])
    labels = {res_ref(s.chain, s.first, s.last): names[s.id] for s in sses_b if s.id in names}
    return labels, crossed


def _own_type(label: str, kind: str) -> str:
    """A generic helix name (α3, η1, primed or not) with the sign of the helix type `kind` (H α, G η): η1 for a 3-10
    helix whose mate is α1; any other name as it is (strands, conventional names such as globins' A-H)."""
    sign = _HELIX_SIGN.get(kind)
    if sign and label[:1] in _HELIX_SIGN.values() and label[:1] != sign and label[1:2].isdigit():
        return sign + label[1:]
    return label


def _lettered(name: str, k: int) -> str:
    return name + (_LETTERS[k] if k < len(_LETTERS) else _LETTERS[-1] + "′" * (k - len(_LETTERS) + 1))


def _primed(name: str, k: int) -> str:
    """C′, C″, C‴, then C′′′′ …: an extra strand after strand C (Ig's C′ and C″), never a lower-case letter on an
    upper-case one (Ca)."""
    return name + (_PRIMES[k] if k < len(_PRIMES) else "′" * (k + 1))


def _colours(layout_a: Layout, sses_a: list[SSE], sses_b: list[SSE], mates: Mates, look: Style) -> dict[str, str]:
    """B element id -> colour (like_from), or nothing when `look` colours by data (each structure its own)."""
    if look.color_by not in IDENTITY_COLOURINGS:
        return {}
    from .palette import tint
    from .render import element_colours

    fill, _ = element_colours(layout_a, sses_a, look)
    out = {b: fill[a] for b, a in [*mates.mate.items(), *mates.pieces.items()] if a in fill}
    order = sorted(sses_b, key=lambda s: (s.chain, s.start))
    for s in order:
        if s.id not in mates.extras:
            continue
        kin = [x for x in order if x.chain == s.chain and (look.color_by != "sstype" or x.kind == s.kind)]
        at = kin.index(s)
        before = [x for x in kin[:at][::-1] if x.id in out]
        after = [x for x in kin[at + 1 :] if x.id in out and x.id not in mates.extras]
        nearest = before or after
        if nearest:
            out[s.id] = tint(out[nearest[0].id], EXTRA_TINT)
    return out


def _template(layout_a: Layout, sses_b: list[SSE], mates: Mates, sheets_b: list[Sheet]) -> list[bool]:
    """Each of B's sheets drawn as a βα row as its mates' sheet in A is (like_from)."""
    rows = {k for row in layout_a.ba_rows for k in row}
    own = template_choices(sheets_b, sses_b, {})
    out = []
    for sheet, mine in zip(sheets_b, own):
        mated = [mates.mate[s.id] for s in sheet.strands if s.id in mates.mate]
        inrow = sum(a in rows for a in mated)
        choice = mine if 2 * inrow == len(mated) else 2 * inrow > len(mated)
        out.append(bool(choice) and can_be_row(sheet))
    return out


def counted(mates: Mates) -> tuple[int, int]:
    """(B's elements matched, B's elements): mates and the pieces of a helix A has whole count as matched."""
    matched = len(mates.mate) + len(mates.pieces)
    return matched, matched + len(mates.extras)


def _coverage(aln: ChainAlignment) -> str:
    """'; only N% of the shorter chain aligned' when the aligned residues are under LOW_COVERAGE of the shorter chain
    (a hinge between domains, a circular permutation, or a chain modelled only in part: the note does not guess
    which), else ''. N is rounded down, so the clause never says more aligned than did."""
    shorter = min(len(aln.residues_a), len(aln.residues_b))
    if not shorter or len(aln.pairs) >= LOW_COVERAGE * shorter:
        return ""
    return f"; only {100 * len(aln.pairs) // shorter}% of the shorter chain aligned"


def reference_name(structure) -> str:
    """The reference as the notes name it: a file by its name (A.cif), a PDB id as given."""
    return Path(str(structure)).name


def like_note(
    name: str,
    chain: str,
    aln: ChainAlignment,
    mates: Mates,
    extras: list[str],
    match: str = "auto",
    crossed: list[tuple[str, str]] = (),
) -> str:
    """The note of a figure drawn like the reference `name` (its chain `chain`): how many of this structure's
    elements were matched, by which alignment, the names of those only this structure has (`extras`, as drawn) and
    of the mates named by their own helix type (`crossed`: (name as drawn, the mate's name), Like.crossed). E.g.
    "drawn like A.cif chain A: 31 of 34 elements matched (structure alignment, TM-score 0.61 over 388 residues);
    only in this structure: α3a, β9a; named by their own helix type: α2a for η1"."""
    matched, total = counted(mates)
    method = _METHODS.get(aln.method) or f"alignment {Path(str(match)).name}"
    text = (
        f"drawn like {name} chain {chain}: {matched} of {total} elements matched ({method}, TM-score "
        f"{aln.tm_score:.2f} over {len(aln.pairs)} residues)"
    )
    if extras:
        text += f"; only in this structure: {', '.join(extras)}"
    if crossed:
        text += f"; named by their own helix type: {', '.join(f'{b} for {a}' for b, a in crossed)}"
    return text + _coverage(aln)


_ADVICE = {  # what draws a figure whose frame cannot be carried like the reference (uncarried_note)
    "copies": "draw one chain of each",
    "stack": "draw it with --mode projected",
    "duplex": "draw it with --no-dna",
}


def uncarried_note(name: str, why: tuple[str, str], mine: bool) -> str:
    """The note of a figure drawn on its own because a page frame cannot be carried (cli._uncarried's (kind, what)):
    this figure's own (`mine`) or the reference's. E.g. "not drawn like A.cif: this structure's frame (a C5 rosette)
    cannot be carried; draw one chain of each", "not drawn like A.cif: its frame (--mode stack) cannot be carried;
    draw it with --mode projected"."""
    kind, what = why
    whose = "this structure's frame" if mine else "its frame"
    return f"not drawn like {name}: {whose} ({what}) cannot be carried; {_ADVICE.get(kind, _ADVICE['copies'])}"


def refused_note(name: str, aln: ChainAlignment, mates: Mates) -> str:
    """The note of a figure drawn on its own because it is not the reference's fold (Mates.same_fold), e.g. "not
    drawn like A.cif: not the same fold (TM-score 0.31; 9 of 30 elements matched)"."""
    matched, total = counted(mates)
    return (
        f"not drawn like {name}: not the same fold (TM-score {aln.tm_score:.2f}; {matched} of {total} elements "
        f"matched)" + _coverage(aln)
    )
