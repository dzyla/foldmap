# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A structure B drawn like a reference A (plan 13a, --like): what B's layout takes from A's finished one.

The reference enters B's layout at its frame (frame.frame_like: A's frame carried through the superposition) and at
its homes and choices (build_layout's `like`): each of B's elements with a mate in A (align.element_mates) starts at
its mate's place on A's page, held there harder than a projected home (layout._K_LIKE), with its mate's drawn angle,
label and, for identity colourings, colour; B's sheets are drawn as βα rows where their mates' are, and B's helices
over their sheet where their mates are. The optimiser itself is not changed: B's own pairing decides each sheet's
strand order (science S1-S3, S6), and _separate clears every overlap as always.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .align import Mates
from .layout import (
    _K_HOME,
    _K_LIKE,
    _PIECE_GAP,
    ETA_W,
    HELIX_RISE,
    HELIX_W,
    SCALE,
    TEN_RISE,
    Layout,
    _row_shape,
    template_choices,
)
from .model import SSE, Sheet, generic_label, res_ref
from .style import Style

IDENTITY_COLOURINGS = ("chain", "sequence", "shade", "sstype", "domain")  # a mate takes its mate's colour
EXTRA_TINT = 0.35  # an extra element's colour: its N-terminal neighbour's, mixed this much with white
_LETTERS = "abcdefghijklmnopqrstuvwxyz"
_PRIMES = ("′", "″", "‴")


@dataclass(frozen=True)
class Like:
    """What B's layout takes from the reference's finished layout (build_layout's `like`)."""

    target: dict[str, tuple[float, float]]  # B element id -> page centre of its mate in A's finished layout
    angle: dict[str, float]  # B element id -> its mate's drawn angle
    label: dict[str, str]  # B element ref (res:chain:a-b) -> label (mates, pieces, extras)
    colour: dict[str, str]  # B element id -> colour (identity colourings only; else empty)
    template: list[bool] | None  # A's βα-row choices, per sheet of B (its own for a sheet with no mate); None: B's own
    over: dict[str, bool]  # A's helix-over-sheet choices, by B helix id
    weight: float  # _K_LIKE / _K_HOME
    # B element id -> (cx, cy) of the unit its mate was laid out in and the mate's offset there (Layout.units): a
    # unit laid out like its mate's lands exactly on that unit's centre (build_layout's homes, layout._like_homes)
    unit: dict[str, tuple[float, float, float, float]] = field(default_factory=dict)


def _drawn_length(s: SSE, look: Style) -> float:
    """A helix's drawn length, as build_layout draws it (_helix_item: true length, at least its width)."""
    if s.kind == "G":
        return max(HELIX_W * ETA_W * look.helix_scale, (len(s) - 1) * TEN_RISE * SCALE)
    return max(HELIX_W * look.helix_scale, (len(s) - 1) * HELIX_RISE * SCALE)


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
      mate and its pieces, Mates.pieces) is laid end to end along its mate as drawn, in sequence, _PIECE_GAP apart,
      each at its own drawn length and the mate's angle, centred on the mate.
    - Labels: a mate takes its mate's label, a piece the same with a prime (α3′); an extra is named after the
      matched element of its family (strand or helix) before it, lettered (α3a, β9a), the first name not taken; an
      extra with none before it keeps its own label unless that is taken, then lettered too.
    - Colours (identity colourings only, IDENTITY_COLOURINGS): a mate or piece takes its mate's colour as A draws it
      (render.element_colours), an extra a tint of its N-terminal neighbour's (the first coloured element after it
      when none comes before).
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
        lengths = [_drawn_length(by_b[b], look) for b in bs]
        along = -(sum(lengths) + _PIECE_GAP * (len(bs) - 1)) / 2
        d = np.array([np.cos(p.angle), np.sin(p.angle)])
        for b, n in zip(bs, lengths):
            t = along + n / 2
            target[b], angle[b] = (float(p.cx + d[0] * t), float(p.cy + d[1] * t)), p.angle
            along += n + _PIECE_GAP
    return Like(
        target=target,
        angle=angle,
        label=_labels(placed, sses_b, mates),
        colour=_colours(layout_a, sses_a, sses_b, mates, look),
        template=_template(layout_a, sses_b, mates, sheets_b) if sheets_b is not None else None,
        over={
            b: a in layout_a.on_sheet
            for b, a in [*mates.mate.items(), *mates.pieces.items()]
            if b in by_b and by_b[b].kind == "H"
        },
        weight=_K_LIKE / _K_HOME,
        unit=unit,
    )


def _labels(placed: dict, sses_b: list[SSE], mates: Mates) -> dict[str, str]:
    """B element ref -> label (like_from: mates, pieces, extras)."""
    names: dict[str, str] = {}  # B element id -> label
    for b, a in mates.mate.items():
        if a in placed:
            names[b] = placed[a].label
    pieces: dict[str, int] = {}
    for b in sorted(mates.pieces, key=lambda b: next(s.start for s in sses_b if s.id == b)):
        a = mates.pieces[b]
        if a in placed:
            k = pieces.get(a, 0)
            pieces[a] = k + 1
            names[b] = placed[a].label + (_PRIMES[k] if k < len(_PRIMES) else "′" * (k + 1))
    taken = {p.label for p in placed.values()} | set(names.values())
    extras = set(mates.extras)
    last: dict[tuple[str, str], str] = {}  # (chain, family) -> the label of the last matched element so far
    counts: dict[tuple[str, str], int] = {}  # B's own generic labels, as layout._assign_labels gives them
    for s in sorted(sses_b, key=lambda s: (s.chain, s.start)):
        key = (s.chain, _family(s))
        n = counts.get((s.chain, s.kind), 0)
        counts[(s.chain, s.kind)] = n + 1
        if s.id in names and s.id in mates.mate:
            last[key] = names[s.id]
        elif s.id in extras:
            base = last.get(key)
            if base is None:
                own = generic_label(s.kind, n)
                if own not in taken:
                    taken.add(own)
                    continue  # its own label, free
                base = own
            k = 0
            while _lettered(base, k) in taken:
                k += 1
            names[s.id] = _lettered(base, k)
            taken.add(names[s.id])
    return {res_ref(s.chain, s.first, s.last): names[s.id] for s in sses_b if s.id in names}


def _lettered(name: str, k: int) -> str:
    return name + (_LETTERS[k] if k < len(_LETTERS) else _LETTERS[-1] + "′" * (k - len(_LETTERS) + 1))


def _colours(layout_a: Layout, sses_a: list[SSE], sses_b: list[SSE], mates: Mates, look: Style) -> dict[str, str]:
    """B element id -> colour (like_from), or nothing when `look` colours by data (each structure its own)."""
    if look.color_by not in IDENTITY_COLOURINGS:
        return {}
    from .palette import tint
    from .render import element_colours

    fill, _ = element_colours(layout_a, sses_a, look)
    out = {b: fill[a] for b, a in [*mates.mate.items(), *mates.pieces.items()] if a in fill}
    order = sorted(sses_b, key=lambda s: (s.chain, s.start))
    for k, s in enumerate(order):
        if s.id not in mates.extras:
            continue
        before = [x for x in order[:k][::-1] if x.chain == s.chain and x.id in out]
        after = [x for x in order[k + 1 :] if x.chain == s.chain and x.id in out and x.id not in mates.extras]
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
        out.append(bool(choice) and _row_shape(sheet))
    return out
