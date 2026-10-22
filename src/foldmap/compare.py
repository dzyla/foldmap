# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Two structures of one fold on one page (plan 13a, `foldmap compare`): A on the left, B drawn like A on the right;
and the same page with what changed drawn on it (plan 13e, `foldmap diff`).

Both panels are drawn at one page scale, the inches per page unit of the larger panel's page (render.page_size and
layout.page_scale), so a helix as long in both is as long on the page. The page has one legend (the panels' entries,
each once) and, under both, the note of B drawn like A (or why B was drawn on its own: not the same fold, or a frame
that cannot be carried). Each panel keeps its own gids, prefixed `a/` and `b/`, so a check that reads one figure
reads either panel by its prefix. `mates_json` says what was matched: the alignment (method, TM-score, RMSD,
identity, residue pairs), the mates, the pieces of a helix A has whole, B's extras and A's elements B lacks.

The diff page (diff_page) is the compare page of the two chains diff.diff_structures compares, aligned as it aligns
them, with its changes drawn: each element one structure has and the other lacks (Diff.only_a, only_b) as a dashed
outline in the other's panel (Layout.absent, render._draw_absent), carried there by the rigid 2D fit of the mates'
page centres (figcheck.rigid_fit, H2's fit); B's panel coloured by each residue's Cα displacement after the core
fit (Diff.displacement, through color_by=data: stripes along each element, a colour bar) when a residue moved
MOVED Å or more; and the rest named in the note under both (diff_note): secondary-structure changes, pairing,
domain turns, interface and contact changes. A structure against itself changes nothing: its diff page is its
compare page, byte for byte. `panel_checks` judges each panel as its own figure (H5).
"""

from __future__ import annotations

import copy
import json
import re
import textwrap
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from matplotlib.figure import Figure

from .layout import font_pt, page_fit, page_scale, snap
from .model import res_ref
from .render import draw, draw_page_legend, legible, page, page_legend, page_notes, page_size, pin_fonts
from .route import route_loops
from .style import Style

PANEL_GAP = 0.5  # page units between the two panels' pages (each keeps its own pad)
PREFIXES = ("a/", "b/")  # each panel's gids: A's, B's
SIDES = ("left", "right")  # each panel's provenance items start with its side
NOTE_SIZE = 0.85  # the note's text, as a fraction of the figure's (the context caption's)
NOTE_LINE = 1.25  # its line step, in its own ems
NOTE_PAD = 0.3  # page units above and below it
MOVED = 0.1  # Å: B's panel is coloured by displacement when a residue moved this far (coordinates' precision, 1e-3)
DISPLACEMENT = "Cα displacement (Å)"  # the colour bar's name
_PAIR = re.compile(r'\[\n\s*("[^"\n]*"),\n\s*("[^"\n]*")\n\s*\]')  # two residue numbers, as json.dumps indents them
_SCALAR = r'(?:"[^"\\\n]*"|-?\d[\d.eE+-]*|true|false|null)'
_ROW = re.compile(rf"\[\n\s*((?:{_SCALAR},\n\s*)+{_SCALAR})\n\s*\]")  # a short list of strings and numbers


@dataclass
class Compared:
    """A compare page and what it was drawn from."""

    figure: Any  # matplotlib Figure
    data: dict  # mates_json (diff_json on a diff page)
    note: str  # the note under both panels: how B was drawn like A, or why not (and what changed, on a diff page)
    provenance: list[dict] = field(default_factory=list)  # both panels', each item's text led by its side
    panels: tuple = ()  # (layout, sses, bb) of A, then of B, as drawn
    # each panel as drawn, for checks that judge it as a figure of its own (panel_checks): layout, loops, sses, bb,
    # model, look (as drawn), plain (the page's look), path, chain
    drawn: tuple = ()
    diff: Any = None  # diff.Diff on a diff page


@dataclass
class _Laid:
    """Both panels laid out (compare's work before drawing)."""

    panels: tuple  # (layout, sses, bb, model, title) of A, then of B
    note: str
    aln: Any
    mates: Any
    alike: bool  # B was drawn like A (its places taken from A's)
    paths: tuple[str, str]
    chains: tuple[str, str]


def _lay_out(
    source_a,
    source_b,
    chain_a,
    chain_b,
    match: str,
    look: Style,
    assembly: str,
    same_protein: bool = False,
    names: tuple = (None, None),
    models: tuple = (None, None),
) -> _Laid:
    """Both panels laid out: B like A (cli.make_layout's like; `same_protein`: drawn alike even under the fold's
    TM-score, two states of one protein), else each on its own. `names`: what the titles and notes call A and B
    (None: reference_name of the source; the app passes an upload's name as uploaded). `models`: A's and B's
    io.read_model as `assembly` reads them, when already read (None: read here)."""
    from .align import align_chains, element_mates
    from .cli import make_layout, one_chain, resolve_structure
    from .io import load_backbone, read_model
    from .like import reference_name

    path_b = resolve_structure(str(source_b))
    model_b = models[1] if models[1] is not None else read_model(path_b, assembly)
    like = {"structure": str(source_a), "chain": chain_a, "match": match, "layout": None}
    if models[0] is not None:
        like["model"] = models[0]
    if same_protein:
        like["same_protein"] = True
    if names[0]:
        like["name"] = names[0]
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
    path_a = resolve_structure(str(source_a))
    if got.get("reference") is None:  # B's own frame cannot be carried: drawn on its own; A laid out here
        model_a = models[0] if models[0] is not None else read_model(path_a, assembly)
        present = list(dict.fromkeys(l.chain for l in load_backbone(path_a, assembly, model=model_a).labels))
        one = one_chain([chain_a] if chain_a else None, present, "--chain-a", str(source_a))
        lay_a, sses_a, bb_a = make_layout(path_a, look=look, chains=[one], assembly=assembly, model=model_a)
        got["reference"] = (lay_a, sses_a, bb_a, model_a, one)
    lay_a, sses_a, bb_a, model_a, chain_a = got["reference"]
    aln, mates = got.get("alignment"), got.get("mates")
    if aln is None:  # a frame that cannot be carried stopped before the alignment
        aln = align_chains(bb_a, chain_a, bb_b, chain_b, match)
        mates = element_mates(sses_a, sses_b, aln)
    panels = (
        (lay_a, sses_a, bb_a, model_a, f"{names[0] or reference_name(source_a)} chain {chain_a}"),
        (lay_b, sses_b, bb_b, model_b, f"{names[1] or reference_name(source_b)} chain {chain_b}"),
    )
    alike = got.get("like") is not None
    return _Laid(panels, lay_b.like_note, aln, mates, alike, (str(path_a), str(path_b)), (chain_a, chain_b))


def _sources(laid: _Laid) -> tuple:
    return tuple(zip(laid.paths, laid.chains))


def compare(
    source_a,
    source_b,
    chain_a: str | None = None,
    chain_b: str | None = None,
    match: str = "auto",
    look: Style | None = None,
    assembly: str = "auto",
    names: tuple = (None, None),
) -> Compared:
    """The page of `source_b` (a file or PDB id) drawn like `source_a`, one protein chain of each (`chain_a`,
    `chain_b`, else each its first: ruling 7), aligned by `match` (auto, sequence, structure or a two-row alignment
    file, the reference's chain first), drawn with `look`; both read as `assembly` reads them (auto, asu: as
    deposited, under the deposited chain names; or an assembly id). `names`: what the titles and the note call A and B
    (None: each by its source)."""
    look = (look or Style()).validate()
    laid = _lay_out(source_a, source_b, chain_a, chain_b, match, look, assembly, names=names)
    (lay_a, sses_a, *_), (lay_b, sses_b, *_) = laid.panels
    drawn: list = []
    fig, items = _page(laid.panels, laid.note, look, sources=_sources(laid), drawn=drawn)
    return Compared(
        fig,
        mates_json(lay_a, sses_a, lay_b, sses_b, laid.aln, laid.mates),
        laid.note,
        items,
        tuple(p[:3] for p in laid.panels),
        tuple(drawn),
    )


def diff_page(
    source_a,
    source_b,
    chain_a: str | None = None,
    chain_b: str | None = None,
    match: str = "auto",
    look: Style | None = None,
    assembly: str = "auto",
    domains: dict | None = None,
    names: tuple = (None, None),
) -> Compared:
    """What changed from `source_a` to `source_b` (files or PDB ids; one protein chain of each, `chain_a`,
    `chain_b`, else each its first) drawn on their compare page (module docstring): diff.diff_structures' Diff, the
    page drawn as compare draws it with the alignment route diff takes (auto: by sequence when the chains are one
    protein, A19's rule, else by structure). Each structure is read once and its chain's DSSP run once: the Diff is
    made from the panels as laid out (diff.drawn_side), with the page's alignment and mates. `domains` ({name:
    residue range or ranges of A's chain}) adds each domain's turn to the note. `names`: as compare's.
    Compared.data is diff_json; Compared.diff the Diff."""
    from .diff import diff_sides, drawn_side, read_structure, route
    from .userdata import ResidueData

    look = (look or Style()).validate()
    read_a = read_structure(source_a, chain_a, "A", assembly)
    read_b = read_structure(source_b, chain_b, "B", assembly)
    way = route(read_a.full.only([read_a.chain]), read_b.full.only([read_b.chain]), match)
    one = match == "auto" and way == "sequence"  # one protein by A19's rule (diff.route): its states drawn alike
    laid = _lay_out(source_a, source_b, chain_a, chain_b, way, look, assembly, one, names, (read_a.model, read_b.model))
    (lay_a, sses_a, *rest_a), (lay_b, sses_b, *rest_b) = laid.panels
    sides = [
        drawn_side(source, got, bb, sses, lay.bridges, assembly)
        for source, got, (lay, sses, bb, *_) in zip((source_a, source_b), (read_a, read_b), laid.panels)
    ]
    d = diff_sides(*sides, way, domains, alignment=laid.aln, mates=laid.mates)
    lay_a, lay_b = copy.copy(lay_a), copy.copy(lay_b)  # the marks go on copies: make_layout's own stay as laid out
    carried = _carried(lay_a, lay_b, d) if laid.alike else None
    if carried is not None:
        lay_b.absent, lay_a.absent = carried
    moved = max(d.displacement.values(), default=0.0)
    look_b = look
    if moved >= MOVED:  # B by how far each residue moved: the per-residue data colouring (plan 13d-5)
        lay_b.residue_data = ResidueData(DISPLACEMENT, "number", dict(d.displacement))
        look_b = replace(look, color_by="data")
    laid.panels = ((lay_a, sses_a, *rest_a), (lay_b, sses_b, *rest_b))
    titles = laid.panels[0][4], laid.panels[1][4]
    said = diff_note(d, lay_a, lay_b, *titles, carried is not None, moved >= MOVED, alike=laid.alike)
    note = f"{laid.note}\n{said}" if said else laid.note
    extra = [{"kind": "diff", "text": said, "count": 1, "elements": []}] if said else []
    drawn: list = []
    fig, items = _page(laid.panels, note, look, (look, look_b), extra, _sources(laid), drawn)
    return Compared(
        fig,
        diff_json(lay_a, sses_a, lay_b, sses_b, d),
        note,
        items,
        tuple(p[:3] for p in laid.panels),
        tuple(drawn),
        d,
    )


def _motion(onto, points):
    """Placed -> Placed: an element of the page `points` come from, moved by the rigid 2D fit of `points` onto
    their matched `onto` (figcheck.rigid_fit), its centre and angle snapped (layout.snap)."""
    from .figcheck import rigid_fit

    turn, c_from, c_to = rigid_fit(onto, points)
    spin = float(np.arctan2(turn[0, 1], turn[0, 0]))

    def moved(p):
        cx, cy = snap((np.array([p.cx, p.cy]) - c_from) @ turn + c_to)
        return replace(p, cx=float(cx), cy=float(cy), angle=float(snap(p.angle + spin)))

    return moved


def _carried(lay_a, lay_b, d):
    """(B's absent, A's absent): A's elements B has no mate for, carried onto B's page, and B's onto A's, by the
    rigid fit of the mates' page centres (criteria.like_pairs; B's centre of a helix drawn in pieces their run's,
    as H2 measures it); None under criteria.LIKE_MIN_MATES mates drawn (no fit to carry them by)."""
    from .criteria import LIKE_MIN_MATES, element_centres, like_pairs, part_centres

    pairs = like_pairs(lay_a, lay_b, d.mates)
    if len(pairs) < LIKE_MIN_MATES:
        return None
    at_a, at_b = element_centres(lay_a, [a for a, *_ in pairs]), part_centres(lay_b, pairs)
    to_b, to_a = _motion(at_b, at_a), _motion(at_a, at_b)

    def drawn(lay, ids):
        return [lay.placed[k] for k in ids if k in lay.placed and not lay.placed[k].ghost]

    return [to_b(p) for p in drawn(lay_a, d.only_a)], [to_a(p) for p in drawn(lay_b, d.only_b)]


def _label(lay, k: str) -> str:
    p = lay.placed.get(k)
    return p.label if p is not None and p.label else k


def _counted(rows, word: str) -> str:
    gained = sum(1 for r in rows if r[-1] == "gained")
    lost = len(rows) - gained
    return f"{word}: " + ", ".join(
        t for t in (f"{gained} gained" if gained else "", f"{lost} lost" if lost else "") if t
    )


def diff_note(d, lay_a, lay_b, title_a: str, title_b: str, carried: bool, coloured: bool, alike: bool = True) -> str:
    """What changed, said under both panels ('' when nothing did): A's elements B lacks (outlined on the right when
    `carried`), B's own (outlined on the left when `carried`; the like note names them when B is drawn `alike`, else
    they are named here), and, when either has elements the other lacks and they are not `carried`, why there are no
    outlines (B not drawn `alike`, or fewer than criteria.LIKE_MIN_MATES mates drawn to carry them by); the residues
    whose secondary structure changes, the pairings gained, lost or shifted (by A's names), the domains' turns, the
    interfaces, ligands and contacts gained or lost, and, when `coloured`, the displacement colouring and its largest
    value."""
    name_a, name_b = title_a.rsplit(" chain ", 1)[0], title_b.rsplit(" chain ", 1)[0]
    parts = []
    if d.only_a:
        where = " (dashed outlines on the right)" if carried else ""
        parts.append(f"only in {name_a}: {', '.join(_label(lay_a, k) for k in d.only_a)}{where}")
    if d.only_b and carried:
        parts.append(f"{name_b}'s own elements outlined dashed on the left")
    elif d.only_b and not alike:  # drawn alike, the like note names them ("only in this structure")
        parts.append(f"only in {name_b}: {', '.join(_label(lay_b, k) for k in d.only_b)}")
    if (d.only_a or d.only_b) and not carried:
        from .criteria import LIKE_MIN_MATES

        why = f"{name_b} drawn on its own" if not alike else f"fewer than {LIKE_MIN_MATES} mates drawn to carry them by"
        parts.append(f"no dashed outlines: {why}")
    if d.ss_changes:
        n = len(d.ss_changes)
        parts.append(f"{n} residue{'' if n == 1 else 's'} change{'s' if n == 1 else ''} secondary structure")
    if d.pairing:
        parts.append("pairing: " + ", ".join(f"{_label(lay_a, x)}–{_label(lay_a, y)} {c}" for x, y, c in d.pairing))
    if d.domain_turns:
        own = d.domain_angles
        said = "domains as given turn " + ", ".join(f"{k} {v:.1f}°" for k, v in d.domain_turns.items())
        said += " after the core fit"
        if own:  # core-free: each domain fitted on its own, against the first given
            first = next(iter(own)).split("/", 1)[0]
            each = ", ".join(f"{k.split('/', 1)[1]} {v:.1f}°" for k, v in own.items())
            said += f"; each fitted on its own, against {first} (the first given): {each}"
        parts.append(said)
    c = d.contacts
    if c.get("interfaces"):
        parts.append("interfaces: " + ", ".join(f"chain {p} {ch}" for p, _, ch in c["interfaces"]))
    if c.get("ligands"):
        parts.append("ligands: " + ", ".join(f"{n} {ch}" for n, ch in c["ligands"]))
    for key, word in (
        ("elements", "element contacts"),
        ("interface_elements", "elements at an interface"),
        ("ligand_elements", "elements at a ligand"),
    ):
        if c.get(key):
            parts.append(_counted(c[key], word))
    if coloured:
        (chain, tag), most = max(d.displacement.items(), key=lambda kv: kv[1])  # the first of the largest
        parts.append(
            f"the right panel coloured by Cα displacement after a fit on {len(d.core)} core residues (largest "
            f"{most:.1f} Å, residue {chain}:{tag})"
        )
    return f"changes from {name_a} to {name_b}: " + "; ".join(parts) if parts else ""


def _page(
    panels, note: str, look: Style, looks: tuple | None = None, extra=(), sources=None, drawn: list | None = None
) -> tuple[Figure, list[dict]]:
    """The page: each panel (layout, sses, bb, model, title) drawn in its look (`looks`, else `look`) at the larger
    page's scale, side by side, top aligned, with one legend (laid out in `look`) and the note under both (each of
    its lines a paragraph). Returns the figure and the panels' provenance (then `extra`, then the page's own); each
    panel as drawn is added to `drawn` (Compared.drawn), `sources` giving its (structure path, chain)."""
    looks = looks or (look, look)
    drawn = [] if drawn is None else drawn
    routed = [route_loops(lay, sses, bb) for lay, sses, bb, _, _ in panels]
    sizes = []
    for (lay, sses, _, model, title), loops, mine in zip(panels, routed, looks):
        size = page_size(lay, loops, sses, title, look=mine, model=model)
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
    lines = [line for para in note.split("\n") for line in textwrap.wrap(para, per_line)] if note else []
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
    for k, ((lay, sses, bb, model, title), loops, (w, h, _), mine) in enumerate(zip(panels, routed, sizes, looks)):
        ax = fig.add_axes([left / width, (height - h) / height, w / width, h / height])

        def at(w_u, h_u, _entries, ax=ax, w=w, h=h, title=title):
            if abs(w_u - w) > 1e-9 or abs(h_u - h) > 1e-9:  # drawn as measured (render.page_size)
                raise RuntimeError(f"{title}: page {w_u} x {h_u} drawn, {w} x {h} measured")
            return scale, fig, ax

        draw(lay, loops, sses, title, look=mine, model=model, at=at)
        for artist in ax.findobj(lambda a, ax=ax: a is not ax and a.get_gid() is not None):
            artist.set_gid(PREFIXES[k] + artist.get_gid())
        for item in getattr(fig, "foldmap_provenance", []):
            if item["kind"] != "page_scaled":  # the page's own, said once below
                items.append({**item, "text": f"{SIDES[k]}: {item['text']}", "panel": PREFIXES[k][0]})
        drawn.append(
            SimpleNamespace(
                layout=lay,
                loops=loops,
                sses=sses,
                bb=bb,
                model=model,
                look=mine,
                plain=look,
                path=sources[k][0] if sources else None,
                chain=sources[k][1] if sources else None,
            )
        )
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
                color=legible("#5d6877"),
                linespacing=NOTE_LINE,
            )
        t.set_gid("compare-note")
    fig.foldmap_page_fit = page_fit(width * scale, height * scale, font_pt(scale))
    items += list(extra)
    items += [{"kind": "page_scaled", "text": text, "count": 1, "elements": []} for text in page_notes(fig)]
    fig.foldmap_provenance = items
    return pin_fonts(fig), items


def panel_checks(made: Compared, assembly: str = "auto", found: dict | None = None) -> dict[str, dict]:
    """H5 for each panel of a page (criteria.h5_good_figure): the panel drawn as a figure of its own, its marks and
    colouring as on the page (no title, as alone), passes every S- and A-check (criteria.check_figure) that its
    structure's chain drawn alone in the page's look passes. {"a": A's panel, "b": B's}. `found`, when given, gets
    each side's two check_figure results: {side: (the panel's, the chain drawn alone's)}."""
    from .criteria import check_figure, h5_good_figure
    from .spec import Drawn, FigureSpec, draw_spec

    out = {}
    for side, p in zip("ab", made.drawn):
        spec = FigureSpec(structure=Path(p.path), look=p.look, chains=[p.chain], assembly=assembly)
        figure = draw(p.layout, p.loops, p.sses, None, look=p.look, model=p.model)
        panel = Drawn(figure, p.loops, p.layout, p.sses, p.bb, p.model)
        alone_spec = replace(spec, look=p.plain)
        alone = draw_spec(alone_spec)
        judged = check_figure(panel, spec, p.path), check_figure(alone, alone_spec, p.path)
        if found is not None:
            found[side] = judged
        out[side] = h5_good_figure(*judged)
    return out


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


DIFF_KEYS = (  # diff_json's own
    "ss_changes",
    "pairing",
    "displacement",
    "core",
    "domain_turns",
    "domain_angles",
    "contacts",
)


def diff_json(lay_a, sses_a, lay_b, sses_b, d) -> dict:
    """What changed, for foldmap diff --json: mates_json of the Diff's alignment and mates (B's elements A lacks are
    its extras, A's that B lacks only_a), then
    ss_changes [{chain, residue, a, b}] (each residue whose element kind differs: H, G, E, C none, - absent; A's
    chain and number, B's for a residue only B has), pairing [{strands: [A's names], change}], displacement
    [[B's chain, B's residue, Å]] (B's order), core (the number of residue pairs the core fit kept), domain_turns
    {name: degrees} (domains given: each domain's turn after the core fit), domain_angles {"FIRST/NAME": degrees}
    (two domains or more: each after the first against it, each fitted on its own, no core) and contacts {key:
    [[...parts, gained | lost]]} (diff.CONTACT_KEYS; elements by A's names); numbers to 1e-3, every list in the
    Diff's order."""
    from .diff import CONTACT_KEYS

    out = mates_json(lay_a, sses_a, lay_b, sses_b, d.alignment, d.mates)
    named = {"elements": (0, 1), "interface_elements": (0,), "ligand_elements": (0,)}  # the parts naming elements

    def row(key: str, parts) -> list:
        return [_label(lay_a, x) if k in named.get(key, ()) else x for k, x in enumerate(parts)]

    out.update(
        {
            "ss_changes": [{"chain": c, "residue": r, "a": x, "b": y} for c, r, x, y in d.ss_changes],
            "pairing": [{"strands": [_label(lay_a, x), _label(lay_a, y)], "change": c} for x, y, c in d.pairing],
            "displacement": [[c, r, _round(v)] for (c, r), v in d.displacement.items()],
            "core": len(d.core),
            "domain_turns": {k: _round(v) for k, v in d.domain_turns.items()},
            "domain_angles": {k: _round(v) for k, v in d.domain_angles.items()},
            "contacts": {key: [row(key, parts) for parts in d.contacts.get(key, ())] for key in CONTACT_KEYS},
        }
    )
    return out


def json_text(data: dict) -> str:
    """`data` (mates_json or diff_json) as write_json writes it: keys sorted, LF line ends; a residue pair, and any
    other short list of strings and numbers, on one line."""
    text = json.dumps(data, indent=1, sort_keys=True, ensure_ascii=False)
    text = _PAIR.sub(r"[\1, \2]", text)  # a residue pair on one line
    text = _ROW.sub(lambda m: "[" + re.sub(r",\n\s*", ", ", m[1]) + "]", text)
    return text + "\n"


def write_json(data: dict, path) -> Path:
    """`data` (mates_json or diff_json) written to `path` (json_text)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json_text(data))
    return path
