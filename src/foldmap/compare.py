# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Two structures of one fold on one page (plan 13a, `foldmap compare`): A on the left, B drawn like A on the right.

Both panels are drawn at one page scale, the inches per page unit of the larger panel's page (render.page_size and
layout.page_scale), so a helix as long in both is as long on the page. The page has one legend (the panels' entries,
each once) and, under both, the note of B drawn like A (or why B was drawn on its own: not the same fold, or a frame
that cannot be carried). Each panel keeps its own gids, prefixed `a/` and `b/`, so a check that reads one figure
reads either panel by its prefix. `mates_json` says what was matched: the alignment (method, TM-score, RMSD,
identity, residue pairs), the mates, the pieces of a helix A has whole, B's extras and A's elements B lacks.
"""

from __future__ import annotations

import json
import re
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from matplotlib.figure import Figure

from .layout import font_pt, page_fit, page_scale
from .model import res_ref
from .render import _legible, draw, draw_page_legend, page, page_legend, page_notes, page_size, pin_fonts
from .route import route_loops
from .style import Style

PANEL_GAP = 0.5  # page units between the two panels' pages (each keeps its own pad)
PREFIXES = ("a/", "b/")  # each panel's gids: A's, B's
SIDES = ("left", "right")  # each panel's provenance items start with its side
NOTE_SIZE = 0.85  # the note's text, as a fraction of the figure's (the context caption's)
NOTE_LINE = 1.25  # its line step, in its own ems
NOTE_PAD = 0.3  # page units above and below it
_PAIR = re.compile(r'\[\n\s*("[^"\n]*"),\n\s*("[^"\n]*")\n\s*\]')  # two residue numbers, as json.dumps indents them


@dataclass
class Compared:
    """A compare page and what it was drawn from."""

    figure: Any  # matplotlib Figure
    data: dict  # mates_json
    note: str  # the note under both panels: how B was drawn like A, or why not
    provenance: list[dict] = field(default_factory=list)  # both panels', each item's text led by its side
    panels: tuple = ()  # (layout, sses, bb) of A, then of B


def compare(
    source_a,
    source_b,
    chain_a: str | None = None,
    chain_b: str | None = None,
    match: str = "auto",
    look: Style | None = None,
    assembly: str = "auto",
) -> Compared:
    """The page of `source_b` (a file or PDB id) drawn like `source_a`, one protein chain of each (`chain_a`,
    `chain_b`, else each its first: ruling 7), aligned by `match` (auto, sequence, structure or a two-row alignment
    file, the reference's chain first), drawn with `look`; both read as `assembly` reads them (auto, asu: as
    deposited, under the deposited chain names; or an assembly id)."""
    from .align import align_chains, element_mates
    from .cli import _one_chain, _resolve_structure, make_layout
    from .io import load_backbone, read_model
    from .like import reference_name

    look = (look or Style()).validate()
    path_b = _resolve_structure(str(source_b))
    model_b = read_model(path_b, assembly)
    like = {"structure": str(source_a), "chain": chain_a, "match": match, "layout": None}
    got: dict = {}
    lay_b, sses_b, bb_b = make_layout(
        path_b,
        look=look,
        chains=[chain_b] if chain_b else None,
        assembly=assembly,
        like=like,
        model=model_b,
        like_out=got,
    )
    chain_b = bb_b.labels[0].chain if len(bb_b) else (chain_b or "")
    if got.get("reference") is None:  # B's own frame cannot be carried: drawn on its own; A laid out here
        path_a = _resolve_structure(str(source_a))
        model_a = read_model(path_a, assembly)
        present = list(dict.fromkeys(l.chain for l in load_backbone(path_a, assembly, model=model_a).labels))
        one = _one_chain([chain_a] if chain_a else None, present, "--chain-a", str(source_a))
        lay_a, sses_a, bb_a = make_layout(path_a, look=look, chains=[one], assembly=assembly, model=model_a)
        got["reference"] = (lay_a, sses_a, bb_a, model_a, one)
    lay_a, sses_a, bb_a, model_a, chain_a = got["reference"]
    aln, mates = got.get("alignment"), got.get("mates")
    if aln is None:  # a frame that cannot be carried stopped before the alignment
        aln = align_chains(bb_a, chain_a, bb_b, chain_b, match)
        mates = element_mates(sses_a, sses_b, aln)
    note = lay_b.like_note
    panels = (
        (lay_a, sses_a, bb_a, model_a, f"{reference_name(source_a)} chain {chain_a}"),
        (lay_b, sses_b, bb_b, model_b, f"{reference_name(source_b)} chain {chain_b}"),
    )
    fig, items = _page(panels, note, look)
    return Compared(
        fig,
        mates_json(lay_a, sses_a, lay_b, sses_b, aln, mates),
        note,
        items,
        tuple(p[:3] for p in panels),
    )


def _page(panels, note: str, look: Style) -> tuple[Figure, list[dict]]:
    """The page: each panel (layout, sses, bb, model, title) drawn at the larger page's scale, side by side, top
    aligned, with one legend and the note under both. Returns the figure and the panels' provenance."""
    routed = [route_loops(lay, sses, bb) for lay, sses, bb, _, _ in panels]
    sizes = []
    for (lay, sses, _, model, title), loops in zip(panels, routed):
        size = page_size(lay, loops, sses, title, look=look, model=model)
        if size is None:
            raise ValueError(f"{title} has no helices or strands to draw")
        sizes.append(size)
    (wa, ha, entries_a), (wb, hb, entries_b) = sizes
    larger = (wa, ha) if wa * ha >= wb * hb else (wb, hb)
    scale = page_scale(*larger, look.font_scale)  # the larger page's inches per page unit, for both
    width = wa + PANEL_GAP + wb
    legend = page_legend(list(dict.fromkeys([*entries_a, *entries_b])), width, look)
    font = font_pt(scale) * look.font_scale * NOTE_SIZE
    per_line = max(30, int((width - 2) * scale * 72 / (font * 0.52)))  # characters per line, as the context caption
    lines = textwrap.wrap(note, per_line) if note else []
    note_h = (2 * NOTE_PAD + len(lines) * font * NOTE_LINE / (scale * 72)) if lines else 0.0
    top = max(ha, hb)
    height = top + legend.height + note_h
    fig = Figure(figsize=(width * scale, height * scale), facecolor=look.background)
    base = fig.add_axes([0, 0, 1, 1])  # the page's own page units: the legend and the note
    base.set_xlim(0, width)
    base.set_ylim(0, height)
    base.set_aspect("equal")
    base.axis("off")
    items: list[dict] = []
    left = 0.0
    for k, ((lay, sses, _, model, title), loops, (w, h, _)) in enumerate(zip(panels, routed, sizes)):
        ax = fig.add_axes([left / width, (height - h) / height, w / width, h / height])

        def at(w_u, h_u, _entries, ax=ax, w=w, h=h, title=title):
            if abs(w_u - w) > 1e-9 or abs(h_u - h) > 1e-9:  # drawn as measured (render.page_size)
                raise RuntimeError(f"{title}: page {w_u} x {h_u} drawn, {w} x {h} measured")
            return scale, fig, ax

        draw(lay, loops, sses, title, look=look, model=model, at=at)
        for artist in ax.findobj(lambda a, ax=ax: a is not ax and a.get_gid() is not None):
            artist.set_gid(PREFIXES[k] + artist.get_gid())
        for item in getattr(fig, "foldmap_provenance", []):
            if item["kind"] != "page_scaled":  # the page's own, said once below
                items.append({**item, "text": f"{SIDES[k]}: {item['text']}", "panel": PREFIXES[k][0]})
        left += w + PANEL_GAP
    draw_page_legend(fig, base, legend, (0.0, width, note_h, height), look, scale)
    if lines:
        with page(look):
            t = base.text(
                1.0,
                note_h - NOTE_PAD,
                "\n".join(lines),
                ha="left",
                va="top",
                fontsize=font,
                style="italic",
                color=_legible("#5d6877"),
                linespacing=NOTE_LINE,
            )
        t.set_gid("compare-note")
    fig.foldmap_page_fit = page_fit(width * scale, height * scale, font_pt(scale))
    items += [{"kind": "page_scaled", "text": text, "count": 1, "elements": []} for text in page_notes(fig)]
    fig.foldmap_provenance = items
    return pin_fonts(fig), items


def _round(x: float) -> float:
    return float(round(float(x), 3)) + 0.0  # to 1e-3, and never -0.0


def mates_json(lay_a, sses_a, lay_b, sses_b, aln, mates) -> dict:
    """What was matched, for --json: {alignment: {method, tm_score, rmsd, identity, pairs}, mates, pieces, extras,
    only_a}. Each element is named as drawn (its label on its panel) with its residues (res:CHAIN:FIRST-LAST); pairs
    are residue numbers [A's, B's] as the files write them (insertion codes included). Lists in sequence order
    (mates, pieces and extras by B's, only_a by A's), numbers rounded to 1e-3; json.dumps(sort_keys=True) writes the
    keys sorted (write_json)."""
    a_of = {s.id: s for s in sses_a}
    b_of = {s.id: s for s in sses_b}

    def named(side: str, lay, s) -> dict:
        p = lay.placed.get(s.id)
        label = p.label if p is not None and p.label else s.id
        return {side: label, f"{side}_residues": res_ref(s.chain, s.first, s.last)}

    def order(s):
        return (s.chain, s.start, s.id)

    def pairs_of(table: dict[str, str]) -> list[dict]:
        out = [(b_of[b], a_of[a]) for b, a in table.items() if b in b_of and a in a_of]
        return [{**named("a", lay_a, a), **named("b", lay_b, b)} for b, a in sorted(out, key=lambda ba: order(ba[0]))]

    ra, rb = aln.residues_a, aln.residues_b
    return {
        "alignment": {
            "method": aln.method,
            "tm_score": _round(aln.tm_score),
            "rmsd": _round(aln.rmsd),
            "identity": _round(aln.identity),
            "pairs": [[ra[i].tag, rb[j].tag] for i, j in aln.pairs],
        },
        "mates": pairs_of(mates.mate),
        "pieces": pairs_of(mates.pieces),
        "extras": [named("b", lay_b, s) for s in sorted((b_of[b] for b in mates.extras if b in b_of), key=order)],
        "only_a": [named("a", lay_a, s) for s in sorted((a_of[a] for a in mates.only_a if a in a_of), key=order)],
    }


def write_json(data: dict, path) -> Path:
    """`data` (mates_json) written to `path`, keys sorted, LF line ends."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        text = json.dumps(data, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write(_PAIR.sub(r"[\1, \2]", text) + "\n")  # a residue pair on one line
    return path
