# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The sequence figure: secondary structure drawn over a chain's sequence, or over a multiple sequence
alignment through its reference row, in rows of N columns (ESPript-style conservation boxes)."""

from __future__ import annotations

from pathlib import Path as FsPath

import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Circle, PathPatch, Polygon, Rectangle

from .layout import ETA_W, HEAD, HEAD_LEN
from .palette import darken, text_color_on, tint
from .render import (
    DATA_NONE,
    _mix,
    add_stripes,
    body_paint,
    data_symbol,
    eta_box,
    helix_ribbon,
    legible,
    page,
    pin_fonts,
    strand_arrow,
    stripe_paint,
)
from .sequence import GAP, Alignment, chain_track, column_classes, conservation, map_alignment
from .style import Style

STRICT = "#c81e1e"  # identical columns: white letters on red
FRAME = "#2b5fb4"  # runs of conserved columns are framed in blue
UNMODELLED = "#9aa0a6"
BAR = "#8a97a8"
INK = "#1d1d1f"
_MONO = ("Liberation Mono", "DejaVu Sans Mono", "Menlo", "Consolas", "Courier New")  # residue letters


def _mono() -> str:
    """The first installed monospaced font (asking matplotlib for missing ones floods the terminal)."""
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    return next((name for name in _MONO if name in installed), "monospace")


_LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5  # in column widths


def draw_sequence(
    path,
    *,
    chain: str | None = None,
    alignment: Alignment | None = None,
    reference: str | None = None,
    columns: int = 60,
    look: Style | None = None,
    full_sequence: bool = False,
    ss_colour: str = "figure",
    threshold: float = 0.7,
    conservation_bar: bool = True,
    title: str | None = None,
    uniprot: str | None = None,
    msa: str | None = None,
    data: str | None = None,
    data_column: str | None = None,
) -> Figure:
    """ss_colour: 'figure' paints each element as in the topology figure drawn with `look`; 'black' is the
    classic plain style. msa: the alignment or ConSurf table file the figure's colours come from (theme
    'conservation'); `alignment` is what the sequence view shows. data, data_column: the researcher's per-residue
    values (--data, --data-column), drawn under color_by=data as the topology figure draws them: stripes (or each
    element's mean colour) for numbers, a mark on the residue for categories."""
    if columns < 10:
        raise ValueError("columns must be at least 10")
    if ss_colour not in ("figure", "black"):
        raise ValueError("ss_colour is 'figure' or 'black'")
    look = (look or Style()).validate()
    with page(look):
        return pin_fonts(
            _draw(
                path,
                chain,
                alignment,
                reference,
                columns,
                look or Style(),
                full_sequence,
                ss_colour,
                threshold,
                conservation_bar,
                title,
                uniprot,
                msa,
                data,
                data_column,
            )
        )


def _element_colours(path, look: Style, ss_colour: str, msa=None, reference=None, data=None, data_column=None):
    """(element id -> colour, element id -> its residues' stripe colours, residue index -> (colour, shape, key) of
    its data mark), as the topology figure drawn with `look` paints them (render.element_colours, data_stripes,
    data_marks' categories)."""
    if ss_colour == "black":
        return {}, {}, {}
    from .cli import make_layout
    from .render import data_categories, data_stripes, data_values, element_colours

    layout, sses, _ = make_layout(
        path, look=look, assembly="asu", msa=msa, msa_reference=reference, data=data, data_column=data_column
    )
    colours = element_colours(layout, sses, look)[0]
    kinds = data_categories(layout) if look.color_by == "data" else {}
    marks = {
        k: (*kinds[v], f"{layout.res_labels[k].chain}:{layout.res_labels[k].tag}")
        for k, v in enumerate(data_values(layout))
        if v is not None and kinds
    }
    return colours, data_stripes(layout, sses, look, colours), marks


REGION_COLOURS = {
    "Domain": "#4c78a8",
    "Region": "#9c9c9c",
    "Repeat": "#72b7b2",
    "Zinc finger": "#b279a2",
    "DNA binding": "#f58518",
    "Motif": "#e45756",
    "Transmembrane": "#c49c3b",
    "Coiled coil": "#54a24b",
    "Signal": "#ff9da6",
    "Propeptide": "#bab0ac",
    "Topological domain": "#d8d2c4",
}
_FEAT_H, _SITE_H, _MAX_LANES = 1.25, 1.0, 4


def _uniprot_track(path, track, accession):
    """(regions as (feature, first, last track position, lane), sites as (feature, track position)) for the chain."""
    from .io import load_backbone
    from .uniprot import annotate

    ann = annotate(path, load_backbone(path, "asu"), None if accession == "auto" else accession)
    where = {r: p for p, r in enumerate(track.residue) if r is not None}
    regions, lanes = [], []
    spans = sorted(
        (
            (f, min(where[k] for k in ks if k in where), max(where[k] for k in ks if k in where))
            for f, c, ks in ann.regions
            if c == track.chain and any(k in where for k in ks)
        ),
        key=lambda t: (t[1], -(t[2] - t[1])),
    )
    for f, a, b in spans:
        lane = next((i for i, end in enumerate(lanes) if end < a), None)
        if lane is None:
            if len(lanes) >= _MAX_LANES:
                continue
            lanes.append(-1)
            lane = len(lanes) - 1
        lanes[lane] = b
        regions.append((f, a, b, lane))
    sites = sorted(
        {
            (f.type, where[k]): (f, where[k])
            for f, c, ks in ann.sites
            if c == track.chain and f.type not in ("Glycosylation", "Disulfide bond")
            for k in ks
            if k in where
        }.values(),
        key=lambda t: t[1],
    )
    return regions, len(lanes), sites


def _draw(
    path,
    chain,
    alignment,
    reference,
    columns,
    look,
    full_sequence,
    ss_colour,
    threshold,
    cons_bar,
    title,
    uniprot=None,
    msa=None,
    data=None,
    data_column=None,
):
    track = chain_track(path, chain, full_sequence)
    regions, lanes, sites = _uniprot_track(path, track, uniprot) if uniprot else ([], 0, [])
    feat_h = lanes * _FEAT_H + (_SITE_H if sites else 0.0) + (0.4 if lanes or sites else 0.0)
    colour_of, stripes, marks = _element_colours(path, look, ss_colour, msa, reference, data, data_column)
    if alignment is not None:
        mapping = map_alignment(alignment, track, reference)
        rows, names, position, ref = alignment.rows, alignment.names, mapping.position, mapping.reference
        strict, similar = column_classes(alignment, threshold)
        cons = conservation(alignment) if cons_bar else None
    else:
        rows, names = [track.letters], [f"{FsPath(path).stem}:{track.chain}"]
        position, ref = list(range(len(track.letters))), 0
        strict = similar = np.zeros(len(position), bool)
        cons = None
    width, n_rows = len(position), len(rows)
    mono = _mono()
    col_of = {p: c for c, p in enumerate(position) if p is not None}
    blocks = -(-width // columns)
    block_h = _LABEL_H + _GLYPH_H + _NUM_H + n_rows * _ROW_H + feat_h + (_CONS_H if cons is not None else 0.0) + _GAP_H
    margin = max(len(n) for n in names) * 0.62 + 1.6
    title_h = 2.4 if title else 0.6
    total_h = title_h + blocks * block_h

    font = 8.0 * look.font_scale
    unit = font / 72 * 0.8  # inches per column
    fig = Figure(figsize=((margin + columns + 1.0) * unit, total_h * unit), facecolor=look.background)
    ink, frame_c = look.ink, legible(FRAME)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-margin, columns + 1.0)
    ax.set_ylim(-total_h, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.text(-margin + 0.4, -1.1, title, ha="left", va="center", fontsize=font * 1.25, fontweight="bold", color=ink)

    def top(b: int) -> float:
        return -title_h - b * block_h

    def glyph_y(b: int) -> float:
        return top(b) - _LABEL_H - _GLYPH_H / 2

    def row_y(b: int, k: int) -> float:
        return top(b) - _LABEL_H - _GLYPH_H - _NUM_H - (k + 0.5) * _ROW_H

    for b in range(blocks):
        marker = Rectangle((-margin, top(b) - block_h), margin + columns, block_h, fill=False, lw=0, visible=False)
        marker.set_gid(f"block:{b}")
        ax.add_patch(marker)
        lo, hi = b * columns, min(width, (b + 1) * columns)
        for k, name in enumerate(names):
            t = ax.text(
                -margin + 0.3,
                row_y(b, k),
                name,
                ha="left",
                va="center",
                fontsize=font * 0.85,
                fontweight="bold" if k == ref and n_rows > 1 else "normal",
                color=ink,
            )
            t.set_gid(f"name:{k}:{b}")
        y_top, y_bot = row_y(b, 0) + _ROW_H / 2, row_y(b, n_rows - 1) - _ROW_H / 2
        run = None
        for c in range(lo, hi + 1):  # framed runs of conserved columns
            kept = c < hi and (strict[c] or similar[c])
            if kept and run is None:
                run = c
            elif not kept and run is not None:
                frame = Rectangle(
                    (run - lo - 0.5, y_bot), c - run, y_top - y_bot, fill=False, ec=frame_c, lw=0.8, zorder=3
                )
                frame.set_gid(f"frame:{b}:{run}")
                ax.add_patch(frame)
                run = None
        for c in range(lo, hi):
            x, p = c - lo, position[c]
            if strict[c]:
                box = Rectangle((x - 0.5, y_bot), 1.0, y_top - y_bot, fc=STRICT, ec="none", zorder=1)
                box.set_gid(f"strict:{c}")
                ax.add_patch(box)
            for k, row in enumerate(rows):
                ch = row[c]
                colour, weight = ink, "normal"
                if alignment is None and not track.observed[p]:
                    ch, colour = ch.lower(), UNMODELLED
                elif ch == GAP:
                    ch, colour = "·", UNMODELLED
                elif strict[c]:
                    colour, weight = "white", "bold"
                elif similar[c]:
                    colour, weight = STRICT, "bold"
                t = ax.text(
                    x,
                    row_y(b, k),
                    ch,
                    ha="center",
                    va="center",
                    fontsize=font,
                    color=colour,
                    fontweight=weight,
                    fontfamily=mono,
                    zorder=2,
                )
                t.set_gid(f"res:{k}:{c}")
            if p is not None and track.numbers[p] % 10 == 0:
                t = ax.text(
                    x,
                    top(b) - _LABEL_H - _GLYPH_H - _NUM_H / 2,
                    str(track.numbers[p]),
                    ha="center",
                    va="center",
                    fontsize=font * 0.72,
                    color=UNMODELLED,
                )
                t.set_gid(f"num:{c}")
            if cons is not None:
                base = y_bot - feat_h - _CONS_H + 0.15
                bar = Rectangle((x - 0.42, base), 0.84, max(cons[c], 0.02) * (_CONS_H - 0.4), fc=BAR, ec="none")
                bar.set_gid(f"cons:{c}")
                ax.add_patch(bar)

    from .render import SITE_COLOURS

    def feat_y(b: int) -> float:  # top of the UniProt track in block b
        return row_y(b, n_rows - 1) - _ROW_H / 2 - 0.2

    for n, (f, a, z, lane) in enumerate(regions):
        if a not in col_of or z not in col_of:
            continue
        c0, c1 = col_of[a], col_of[z]
        colour = REGION_COLOURS.get(f.type, "#9c9c9c")
        for piece, b in enumerate(range(c0 // columns, c1 // columns + 1)):
            lo, hi = max(c0, b * columns), min(c1, (b + 1) * columns - 1)
            y = feat_y(b) - lane * _FEAT_H - 0.75
            bar = Rectangle(
                (lo - b * columns - 0.45, y), hi - lo + 0.9, 0.7, fc=colour, ec="none", alpha=0.85, zorder=2
            )
            bar.set_gid(f"uniprot-region:{n}:{piece}")
            ax.add_patch(bar)
            if piece == 0 or lo == b * columns:
                room = hi - lo + 0.9
                text = f.label if len(f.label) * 0.55 < room else f.label[: max(0, int(room / 0.55) - 1)] + "…"
                if room > 2.5:
                    ax.text(
                        lo - b * columns,
                        y + 0.35,
                        text,
                        ha="left",
                        va="center",
                        fontsize=font * 0.68,
                        color=text_color_on(colour),
                        zorder=3,
                        clip_on=True,
                    )
    for n, (f, p) in enumerate(sites):
        if p not in col_of:
            continue
        c = col_of[p]
        b = c // columns
        y = feat_y(b) - lanes * _FEAT_H - _SITE_H / 2
        dot = Circle((c - b * columns, y), 0.3, fc=SITE_COLOURS.get(f.type, "#555555"), ec="none", zorder=3)
        dot.set_gid(f"uniprot-site:{n}")
        ax.add_patch(dot)

    for e in track.elements:
        if e.start not in col_of or e.end not in col_of:
            continue
        c0, c1 = col_of[e.start], col_of[e.end]
        colour = colour_of.get(e.id, "#000000")
        for piece, b in enumerate(range(c0 // columns, c1 // columns + 1)):
            lo, hi = max(c0, b * columns), min(c1, (b + 1) * columns - 1)
            x0, x1 = lo - b * columns - 0.5, hi - b * columns + 0.5
            painted = None
            if e.id in stripes:  # each column its residue's stripe, as the figure's element
                first = track.residue[e.start]
                painted = [
                    stripes[e.id][track.residue[p] - first]
                    if p is not None and track.residue[p] is not None
                    else DATA_NONE
                    for p in (position[c] for c in range(lo, hi + 1))
                ]
            _body(ax, e.kind, x0, x1, glyph_y(b), colour, hi == c1, f"{e.id}:{piece}", look, font, painted)
            if lo == c0:
                t = ax.text(
                    x0,
                    top(b) - _LABEL_H * 0.45,
                    e.label,
                    ha="left",
                    va="center",
                    fontsize=font * 0.9,
                    color=legible(darken(colour, 0.85)) if colour != "#000000" else ink,
                )
                t.set_gid(f"sslabel:{e.id}")

    for c, p in enumerate(position):  # data marks (categories) on their residues' columns, over the track
        k = track.residue[p] if p is not None else None
        if k in marks:
            colour, shape, key = marks[k]
            b = c // columns
            at = (c - b * columns, glyph_y(b))
            ring = Circle(at, _MARK_RING, fc=look.background, ec=legible(INK), lw=0.6 * font / 8 * 0.6, zorder=3)
            ring.set_gid(f"data-ring:{key}")
            ax.add_patch(ring)
            mark = data_symbol(shape, at, _MARK)
            mark.set(facecolor=colour, edgecolor="none", zorder=3.1)
            mark.set_gid(f"data:{key}")
            ax.add_patch(mark)

    ss = track.ss
    k = 0
    while k < len(ss):  # turns (TT), ESPript-style, over runs of DSSP turn residues outside elements
        if ss[k] == "T":
            j = k
            while j + 1 < len(ss) and ss[j + 1] == "T":
                j += 1
            if j > k and k in col_of and j in col_of:
                c = (col_of[k] + col_of[j]) / 2
                b = col_of[k] // columns
                t = ax.text(
                    c - b * columns, glyph_y(b), "TT", ha="center", va="center", fontsize=font * 0.75, color=ink
                )
                t.set_gid(f"turn:{k}")
            k = j + 1
        else:
            k += 1
    return fig


_BODY = 1.25  # a helix ribbon's diameter and an arrowhead's width, column widths (the figure: HELIX_W = HEAD)
_RIBBON = 0.28  # the ribbon's width along its axis, a fraction of a turn: a turn here is 3.6 columns long, longer
# for its diameter than the figure's (render._BAND, 0.36), so a narrower band keeps its faces as wide across


def _body(
    ax, kind: str, x0: float, x1: float, y: float, colour: str, head: bool, key: str, look, font: float, painted=None
):
    """An element's piece in one row, drawn from x0 to x1 as the figure draws it (render's ribbon, arrow and 3-10
    box, painted by look.fill and look.helix_shading): a helix ribbon one turn per 3.6 residues, starting and ending
    on its axis and spanning exactly its columns; a strand's arrow, its head (head) ending at its last residue; a
    3-10 helix a rounded box. gids ss:<key> (the ribbon's front faces, the arrow, the box), ss-back:<key> and
    ss-shine:<key>. painted (color_by=data stripes): each column's colour, a stripe each over the body (the
    figure's render.add_stripes), gids ss-stripe:<key>:<k> and ss-stripe-back:<key>:<k>, the edge drawn again over
    them (ss-edge:<key>, ss-edge-back:<key>)."""
    lw = 0.6 * font / 8
    fill, shade, edge, edge_w = body_paint(kind, colour, look, lw)
    if colour == "#000000" and shade == fill and kind == "H" and look.helix_shading != "none":
        shade = _mix(colour, look.background, 0.55)  # the plain black style: a grey back face keeps the twist
    cx, length = (x0 + x1) / 2, x1 - x0
    if kind == "E":
        k = _BODY / HEAD
        pts = strand_arrow(cx, y, length, 0.0, k, head, HEAD_LEN * k)
        patch = Polygon(pts, closed=True, fc=fill, ec=edge, lw=edge_w, joinstyle="round", zorder=2)
    elif kind == "G":
        patch = eta_box(cx, y, length, _BODY * ETA_W, fc=fill, ec=edge, lw=edge_w, zorder=2)
    else:
        turns = max(1, round(length / 3.6))
        front, back, shine = helix_ribbon(cx, y, length, _BODY, 0.0, turns, _RIBBON)
        span = np.ptp(front.vertices[:, 0])  # its tapered ends reach a little past its length: drawn a little
        front, back, shine = helix_ribbon(cx, y, length * length / span, _BODY, 0.0, turns, _RIBBON)  # shorter
        behind = PathPatch(back, fc=shade, ec=edge, lw=edge_w * 0.8, joinstyle="round", zorder=1.9)
        behind.set_gid(f"ss-back:{key}")
        ax.add_patch(behind)
        if painted:
            _stripes(
                ax, behind, np.vstack([front.vertices, back.vertices]), x0, x1, y, painted, kind, look, lw, key, True
            )
        if shine is not None and look.helix_shading == "gloss" and look.fill != "outline":
            glint, see = ("#ffffff", 0.4) if painted else (tint(fill, 0.4), 0.85)  # over stripes: a white glint
            gloss = PathPatch(shine, fc=glint, ec="none", alpha=see, zorder=2.05)
            gloss.set_gid(f"ss-shine:{key}")
            ax.add_patch(gloss)
        patch = PathPatch(front, fc=fill, ec=edge, lw=edge_w, joinstyle="round", zorder=2)
    patch.set_gid(f"ss:{key}")
    ax.add_patch(patch)
    if painted:
        if kind == "E":
            extent = pts
        elif kind == "G":
            h = _BODY * ETA_W / 2
            extent = np.array([(x0, y - h), (x1, y - h), (x1, y + h), (x0, y + h)])
        else:
            extent = np.vstack([front.vertices, back.vertices])
        _stripes(ax, patch, extent, x0, x1, y, painted, kind, look, lw, key)


_MARK, _MARK_RING = 0.24, 0.36  # a data mark's symbol and ring, column widths


def _stripes(ax, body, extent, x0, x1, y, painted, kind, look, lw, key, back=False) -> None:
    """One stripe per column over a body (clipped to it), and its edge drawn again over them."""
    face = "-back" if back else ""
    centres = [(i + 0.5) / len(painted) for i in range(len(painted))]
    paint = stripe_paint(kind, look, lw, back)
    z = 1.91 if back else 2.01
    add_stripes(ax, body, extent, (x0, y), (x1, y), centres, painted, paint, z, f"ss-stripe{face}:{key}")
    edge = PathPatch(body.get_path(), fc="none", ec=body.get_edgecolor(), lw=body.get_linewidth(), zorder=z + 0.01)
    edge.set(transform=body.get_transform(), joinstyle="round")
    edge.set_gid(f"ss-edge{face}:{key}")
    ax.add_patch(edge)
