# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""How a figure was made (plan 12 R7): one ordered list of what in it is inferred, ambiguous, or moved for
readability, beyond the standard method. The command line prints it, the app shows it under the figure, the MCP
server returns it, and the SVG (<desc>), PDF (Subject) and PNG (Description) files carry it.

Each item is a dict: kind, text (one line, as every interface says it), count (how many things it concerns) and
elements (their ids: element ids, or domain, stem and ligand names). Kinds, in the order listed:

    domains            where the domain panels come from: inferred from contacts, CATH's, UniProt's or as given
    ambiguous_pairing  pairings drawn off a sheet's run, named only, or joining barrels drawn apart (sheet_note_items)
    bulge_joined       strands drawn whole across a beta-bulge DSSP splits them at (ss._join_bulged_pieces)
    helix_matched      helices cut or joined so copies of a protein are drawn alike (ss._copy_breaks)
    ghost_strands      strands drawn again as dashed copies: a barrel's closure, or in another copy's sheet
    names_moved        names set off their usual spot to stay readable
    additives_hidden   crystallisation and purification additives not marked (Style.ligands auto)
    crowded_marks      ligand marks with no clear spot, drawn over other ink
    loops_unrouted     loops or nucleic-acid strands with no clear route, drawn as plain curves
    loops_over_text    loops routed over a label or chain end, the only route found
    page_scaled        a page held at its largest with text under 6 pt; a PNG under 300 dpi

An empty list means nothing was inferred or moved beyond the standard method. The facts are those already kept on
the layout and the figure; this module only words them."""

from __future__ import annotations

from xml.sax.saxutils import escape

HEADING = "How this figure was made:"
WARNINGS = ("crowded_marks", "loops_unrouted", "loops_over_text")  # said as warnings on the command line


def _item(kind: str, text: str, elements, count: int | None = None) -> dict:
    elements = list(elements)
    return {"kind": kind, "text": text, "count": len(elements) if count is None else count, "elements": elements}


def stuck(loops, layout) -> list[str]:
    """Loops with no clear route, drawn as plain curves (a>b), then nucleic-acid strands drawn with no clear way."""
    from .na_paths import path_name

    out = [f"{l.a_id}>{l.b_id}" for l in loops if l.fallback]
    na = getattr(layout, "nucleic", None)
    return out + [path_name(p, na) for p in getattr(layout, "na_paths", []) if p.fallback]


def over_labels(loops) -> list[str]:
    """Loops routed over a label or a chain-end stub, the only route found (route.Loop.crosses_labels)."""
    return [f"{l.a_id}>{l.b_id}" for l in loops if l.crosses_labels and not l.fallback]


def _domains(layout) -> list[dict]:
    named = [name for name, ids in layout.domains if ids]
    stems = ", ".join(dict.fromkeys(name.split(" · ")[0] for name in named))
    if layout.domain_note:  # inferred (or none found), or CATH's (or why there are none): cli.domain_note
        return [_item("domains", layout.domain_note, named)]
    if not named:
        return []
    if layout.domain_source == "UniProt":
        return [_item("domains", f"domains {stems} from UniProt's annotation, not inferred", named)]
    return [_item("domains", f"domains {stems} as given, not inferred", named)]


def _bulges(layout) -> list[dict]:
    joined = sorted(
        (p.sse for p in layout.placed.values() if p.sse.kind == "E" and getattr(p.sse, "bulges", ())),
        key=lambda s: s.start,
    )
    if not joined:
        return []
    said = "; ".join(f"{s.id} (bulge at {', '.join(s.bulges)})" for s in joined)
    text = f"{len(joined)} strand(s) drawn whole across a β-bulge, where DSSP splits them: {said}"
    return [_item("bulge_joined", text, [s.id for s in joined])]


def _helix_matched(layout) -> list[dict]:
    ids = [i for i in getattr(layout, "helix_matched_ids", []) if i in layout.placed]
    if not ids:
        return []
    n = layout.helix_matched
    text = (
        f"helix pieces cut or joined so copies of a protein are drawn alike, where DSSP breaks a helix in some copies "
        f"only ({n.get('cut', 0)} cut, {n.get('joined', 0)} joined): {', '.join(ids)}"
    )
    return [_item("helix_matched", text, ids)]


def _ghosts(layout) -> list[dict]:
    from .layout import block_ghosts

    closing, lent = {}, {}
    for k, block in enumerate(layout.sheets):
        for g in block_ghosts(layout, k):
            (closing if g.sse.id in block else lent)[g.sse.id] = None
    for g in layout.ghosts:  # one no block records: it stands for a strand drawn elsewhere
        if g.sse.id not in closing:
            lent[g.sse.id] = None
    out = []
    if closing:
        text = f"{len(closing)} strand(s) drawn again, dashed, after a barrel's last strand to show it closes"
        out.append(_item("ghost_strands", f"{text}: {', '.join(closing)}", closing))
    if lent:
        text = f"{len(lent)} strand(s) drawn again, dashed, in the sheet of another copy they pair with"
        out.append(_item("ghost_strands", f"{text}: {', '.join(lent)}", lent))
    return out


def _names_moved(fig, layout, labels: bool) -> list[dict]:
    """Helix names on the other side of their helix or moved along it (layout._element_label_sides,
    _over_label_sides, route._clear_helix_names; a mirror copy's name moves when its mate's does), strand letters
    moved along their arrow off a disulfide bar (render.strand_label_spot), domain names above their panel with a
    leader (Layout.domain_leaders) and stem names set further out from their letters (PlacedStem.label_out)."""
    parts, elements = [], []

    def moved(p) -> bool:
        return p.label_along != 0 or p.label_side == -1

    if labels:
        mate_of = dict(getattr(layout, "mirror_pairs", []))  # copy id -> mate id
        helices = []
        for i, p in layout.placed.items():
            if p.sse.kind not in ("H", "G") or not p.label:
                continue
            mate = layout.placed.get(mate_of.get(i, ""))
            if p.label_along != 0 or (moved(mate) if mate is not None else p.label_side == -1):
                helices.append(p)
        if helices:
            parts.append(f"helix names on the other side of their helix or along it: {_names(helices)}")
            elements += [p.sse.id for p in helices]
        letters = [layout.placed[i] for i in getattr(fig, "foldmap_moved_letters", []) if i in layout.placed]
        if letters:
            parts.append(f"strand letters along their arrow, off a disulfide bar: {_names(letters)}")
            elements += [p.sse.id for p in letters]
    leaders = list(layout.domain_leaders)
    if leaders:
        parts.append(f"domain names above their panel, with a leader: {', '.join(leaders)}")
        elements += leaders
    stems = [d.name for d in layout.stems if d.label_out > 0 and d.name]
    if stems:
        parts.append(f"stem names set further out from their letters: {', '.join(stems)}")
        elements += stems
    if not parts:
        return []
    text = f"{len(elements)} name(s) moved from their usual spot to stay readable: {'; '.join(parts)}"
    return [_item("names_moved", text, elements)]


def _names(placed) -> str:
    return ", ".join(f"{p.label} ({p.sse.id})" if p.label else p.sse.id for p in placed)


def figure_provenance(fig, layout, loops=(), labels: bool = True) -> list[dict]:
    """How the figure `fig` was drawn from `layout` and its routed `loops` (see the module docstring); `labels`:
    whether element names are drawn (Style.labels)."""
    from .render import page_notes
    from .scope import _hidden_additives

    out = _domains(layout)
    out += [dict(i) for i in getattr(layout, "sheet_items", [])]
    out += _bulges(layout)
    out += _helix_matched(layout)
    out += _ghosts(layout)
    out += _names_moved(fig, layout, labels)
    additives = _hidden_additives(layout)
    if additives:
        text = (
            f"additives {', '.join(additives)} not marked: crystallisation and purification leftovers "
            "(the style ligands=all marks them)"
        )
        out.append(_item("additives_hidden", text, additives))
    crowded = list(getattr(fig, "foldmap_crowded_marks", []))
    if crowded:
        text = f"{len(crowded)} ligand mark(s) with no clear spot, drawn over other ink: {', '.join(crowded)}"
        out.append(_item("crowded_marks", text, crowded))
    if lost := stuck(loops, layout):
        text = (
            f"no clear route for {len(lost)} loop(s) or strand(s), drawn as plain curves that may cross elements: "
            f"{', '.join(lost)}"
        )
        out.append(_item("loops_unrouted", text, lost))
    if over := over_labels(loops):
        text = f"{len(over)} loop(s) drawn over a label or chain end, the only route found: {', '.join(over)}"
        out.append(_item("loops_over_text", text, over))
    out += [_item("page_scaled", note, [], 1) for note in page_notes(fig)]
    return out


def lines(items: list[dict], warn: bool = False) -> list[str]:
    """The items as the list says them, one line each; `warn`: those in WARNINGS start "warning: "."""
    return [("warning: " if warn and i["kind"] in WARNINGS else "") + i["text"] for i in items]


def as_text(items: list[dict]) -> str:
    """The list as plain text, as the figure files carry it: the heading, then "- " and one item a line ('' when
    the list is empty)."""
    return "\n".join([HEADING, *(f"- {line}" for line in lines(items))]) if items else ""


def with_desc(svg: str, items: list[dict]) -> str:
    """The SVG text with the list (as_text) as its <desc>, the first child of its root; unchanged when empty."""
    text = as_text(items)
    if not text:
        return svg
    at = svg.index(">", svg.index("<svg")) + 1
    return f"{svg[:at]}\n <desc>{escape(text)}</desc>{svg[at:]}"
