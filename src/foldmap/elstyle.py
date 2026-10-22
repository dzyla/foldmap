# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Each element's own colour and width, set by hand (plan 23 Task 1): the style keys `element_colors` and
`element_widths`, comma-separated REF=VALUE entries, as --color REF=COLOUR and --width REF=FACTOR give them.

A reference names an element as --mark and the edits do (layout.resolve: a label, chain:label, #k, res:A:1-60,
dep:/beta:/fig: names), or one of the other things drawn:

    loop:REF        the loop leaving element REF (its C end), or with REF a residue range the loops it meets;
                    a residue range that meets no element names the loops it meets too
    na:CHAIN        a nucleic-acid chain's backbone (stem rails, ribbons, blocks, single strands)
    lig:CODE        every bound ligand or ion drawn with that code (HEM, ZN), lig:CODE:A142 one of them

A colour wins over the figure's colouring (color_by, a mark) for what it names; a width scales the drawn body (a
helix's diameter, an arrow's width, a loop's stroke, a backbone's weight, a ligand box's edge) of that only, and the
layout keeps the room (build_layout widens the element before it is cleared of its neighbours, line_band keeps text
off a wider loop or backbone)."""

from __future__ import annotations

from typing import NamedTuple

COLORS, WIDTHS = "element_colors", "element_widths"
WIDTH_MIN, WIDTH_MAX = 0.25, 4.0  # a width factor's range: thinner than a quarter vanishes, wider than 4 swamps
PREFIXES = ("loop:", "na:", "lig:")
FORMS = (
    "an element (α1, A:G, #3, res:A:1-60, dep:NAME), loop:REF (the loop after REF), na:CHAIN (a nucleic-acid chain) "
    "or lig:CODE (a ligand, lig:HEM or lig:HEM:A142)"
)


def entries(text: str, key: str) -> list[tuple[str, str]]:
    """(reference, value) of each comma-separated REF=VALUE entry of `text`, in order; ValueError (naming `key`) for
    an entry without both. The value is what follows the last '=' (a reference never holds one)."""
    out = []
    for entry in filter(None, (e.strip() for e in str(text or "").split(","))):
        ref, eq, value = (part.strip() for part in entry.rpartition("="))
        if not eq or not ref or not value:
            what = "COLOUR (a name or #rrggbb)" if key == COLORS else f"FACTOR ({WIDTH_MIN:g} to {WIDTH_MAX:g})"
            raise ValueError(f"{key} entries are REF={what}; got {entry!r}")
        out.append((ref, value))
    return out


def colours(text: str) -> list[tuple[str, str]]:
    """element_colors as (reference, #rrggbb) pairs; ValueError names an entry whose colour is none."""
    from matplotlib.colors import to_hex

    out = []
    for ref, value in entries(text, COLORS):
        try:
            out.append((ref, to_hex(value)))
        except ValueError:
            raise ValueError(f"{COLORS}: {ref} has {value!r}, which is not a colour (use a name or #rrggbb)") from None
    return out


def widths(text: str) -> list[tuple[str, float]]:
    """element_widths as (reference, factor) pairs; ValueError names an entry whose factor is no number in range."""
    out = []
    for ref, value in entries(text, WIDTHS):
        try:
            factor = float(value)
        except ValueError:
            factor = float("nan")
        if not WIDTH_MIN <= factor <= WIDTH_MAX:
            raise ValueError(f"{WIDTHS}: {ref} has {value!r}; a width is a factor from {WIDTH_MIN:g} to {WIDTH_MAX:g}")
        out.append((ref, factor))
    return out


def joined(text: str, items) -> str:
    """`text` (a style value of REF=VALUE entries) with `items` (more of them, as --color/--width give them)
    appended: later entries win for what they name."""
    parts = [p for p in (str(text or "").strip(), *(str(i).strip() for i in items or ())) if p]
    return ",".join(parts)


def references(text: str) -> list[str]:
    """The references of a style value's entries, as written (for spec.check_references); [] when it parses not
    (validate says why)."""
    try:
        return [ref for ref, _ in entries(text, COLORS)]
    except ValueError:
        return []


class Targets(NamedTuple):
    """What one reference names: element ids, loops (each by the id of the element it leaves), nucleic-acid chains
    and ligand selectors ((code, chain+number or ''))."""

    elements: list
    loops: list
    chains: list
    ligands: list


def _loops_after(layout, sses) -> dict[str, tuple]:
    """Each placed element that a loop leaves -> (the element, the next one along its chain), in sequence."""
    placed = sorted((s for s in sses if s.id in layout.placed), key=lambda s: s.start)
    out = {}
    for a, b in zip(placed, placed[1:]):
        if a.chain == b.chain:
            out[a.id] = (a, b)
    return out


def _key(r) -> tuple[int, str]:
    return (r.seq, r.icode)


def loops_meeting(layout, sses, ref: str) -> list[str]:
    """The loops (by the id of the element each leaves) whose residues a residue range meets: those strictly
    between the two elements it joins."""
    from .layout import residue_span

    span = residue_span(ref)
    if span is None:
        return []
    chain, lo, hi = span
    lo, hi = min(lo, hi), max(lo, hi)
    out = []
    for a_id, (a, b) in _loops_after(layout, sses).items():
        if chain not in (None, a.chain):
            continue
        first, last = sorted((_key(a.last), _key(b.first)))
        if lo < last and hi > first:
            out.append(a_id)
    return out


def targets(layout, sses, ref: str) -> Targets:
    """What `ref` names in this figure (module doc); ValueError says why it names nothing."""
    from .layout import residue_span, resolve

    ref = ref.strip()
    if ref.startswith("na:"):
        chain = ref[3:].strip()
        held = list(dict.fromkeys(s.chain for s in layout.nucleic.strands)) if layout.nucleic is not None else []
        if chain not in held:
            have = f"its nucleic-acid chains are {', '.join(held)}" if held else "it draws no nucleic acid"
            raise ValueError(f"no nucleic-acid chain {chain!r} drawn; {have}")
        return Targets([], [], [chain], [])
    if ref.startswith("lig:"):
        code, _, which = ref[4:].strip().partition(":")
        if not code.strip():
            raise ValueError("it names no ligand; write lig:CODE (lig:HEM) or lig:CODE:CHAIN+NUMBER (lig:HEM:A142)")
        return Targets([], [], [], [(code.strip().upper(), which.strip())])
    if ref.startswith("loop:"):
        inner = ref[5:].strip()
        after = _loops_after(layout, sses)
        try:  # the loop leaving each element it names (a layout file keeps loop:α1 as loop:res:A:3-11)
            found = [k for k in resolve(layout, sses, inner) if k in after]
        except ValueError:
            if not inner.startswith("res:"):
                raise
            found = loops_meeting(layout, sses, inner)  # a range in a loop: the loops it meets
        if not found:
            raise ValueError("no loop drawn there (a loop leaves an element for the next along its chain)")
        return Targets([], found, [], [])
    try:
        return Targets(resolve(layout, sses, ref), [], [], [])
    except ValueError as err:
        body = ref if ref.startswith("res:") else f"res:{ref}"
        found = loops_meeting(layout, sses, body) if residue_span(body) is not None else []
        if not found or type(err).__name__ == "AmbiguousName":
            raise
        return Targets([], found, [], [])


def ligand_named(ligand, selectors) -> bool:
    """Whether a bound ligand (model.Ligand) is one of `selectors` ((code, chain+number or '')): its name, symbol or
    one of its codes, and its chain and number when given."""
    names = {str(ligand.name).upper(), str(ligand.symbol).upper(), *(str(c).upper() for c in ligand.codes or ())}
    return any(code in names and which in ("", f"{ligand.chain}{ligand.seq}") for code, which in selectors)


class Hand(NamedTuple):
    """Everything set by hand, resolved for one figure: element id / loop (by the element it leaves) / nucleic-acid
    chain -> colour or factor; ligands as [(selectors, value)], later entries winning."""

    elements: dict
    loops: dict
    chains: dict
    ligands: list

    def ligand(self, ligand, default=None):
        """The value set for a bound ligand, else `default`."""
        out = default
        for selectors, value in self.ligands:
            if ligand_named(ligand, selectors):
                out = value
        return out


def resolve_all(layout, sses, pairs, key: str) -> Hand:
    """`pairs` ((reference, value), in order) resolved against the figure: later entries win. ValueError names the
    key and the reference that names nothing."""
    out = Hand({}, {}, {}, [])
    for ref, value in pairs:
        try:
            found = targets(layout, sses, ref)
        except ValueError as err:
            raise ValueError(f"{key} {ref}: {err}") from None
        out.elements.update(dict.fromkeys(found.elements, value))
        out.loops.update(dict.fromkeys(found.loops, value))
        out.chains.update(dict.fromkeys(found.chains, value))
        if found.ligands:
            out.ligands.append((found.ligands, value))
    return out


def canonical(layout, sses, ref: str) -> list[str]:
    """A reference as a layout file keeps it: elements as residue ranges (layoutfile.canonical), a loop as loop:
    and the range of the element it leaves; nucleic-acid chains, ligands and residue ranges as written."""
    from .layoutfile import canonical as element_ranges

    ref = ref.strip()
    if ref.startswith(("na:", "lig:")):
        return [ref]
    found = targets(layout, sses, ref)
    by_id = {s.id: s for s in sses}
    if found.loops and not ref.startswith("loop:"):  # a residue range in a loop: a range already
        return [ref]
    if found.loops:  # loop:α1 -> loop: and the range of the element it leaves; a range in a loop as written
        from .layout import resolve
        from .model import res_ref

        try:
            resolve(layout, sses, ref[5:])
        except ValueError:
            return [ref]
        return [f"loop:{res_ref(by_id[k].chain, by_id[k].first, by_id[k].last)}" for k in found.loops]
    return element_ranges(layout, sses, ref)
