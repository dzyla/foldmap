"""The sequence figure: secondary structure drawn over a chain's sequence, or over a multiple sequence
alignment through its reference row, in rows of N columns (ESPript-style conservation boxes)."""

from __future__ import annotations

from pathlib import Path as FsPath

import numpy as np
from matplotlib import rc_context
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.path import Path

from .palette import darken
from .render import _RC
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
) -> Figure:
    """ss_colour: 'figure' paints each element as in the topology figure drawn with `look`; 'black' is the
    classic plain style."""
    if columns < 10:
        raise ValueError("columns must be at least 10")
    if ss_colour not in ("figure", "black"):
        raise ValueError("ss_colour is 'figure' or 'black'")
    with rc_context(_RC):
        return _draw(
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
        )


def _element_colours(path, look: Style, ss_colour: str) -> dict[str, str]:
    if ss_colour == "black":
        return {}
    from .cli import make_layout
    from .render import element_colours

    layout, sses, _ = make_layout(path, look=look, assembly="asu")
    return element_colours(layout, sses, look)[0]


def _draw(path, chain, alignment, reference, columns, look, full_sequence, ss_colour, threshold, cons_bar, title):
    track = chain_track(path, chain, full_sequence)
    colour_of = _element_colours(path, look, ss_colour)
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
    block_h = _LABEL_H + _GLYPH_H + _NUM_H + n_rows * _ROW_H + (_CONS_H if cons is not None else 0.0) + _GAP_H
    margin = max(len(n) for n in names) * 0.62 + 1.6
    title_h = 2.4 if title else 0.6
    total_h = title_h + blocks * block_h

    font = 8.0 * look.font_scale
    unit = font / 72 * 0.8  # inches per column
    fig = Figure(figsize=((margin + columns + 1.0) * unit, total_h * unit))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-margin, columns + 1.0)
    ax.set_ylim(-total_h, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.text(-margin + 0.4, -1.1, title, ha="left", va="center", fontsize=font * 1.25, fontweight="bold", color=INK)

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
                color=INK,
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
                    (run - lo - 0.5, y_bot), c - run, y_top - y_bot, fill=False, ec=FRAME, lw=0.8, zorder=3
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
                colour, weight = INK, "normal"
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
                base = y_bot - _CONS_H + 0.15
                bar = Rectangle((x - 0.42, base), 0.84, max(cons[c], 0.02) * (_CONS_H - 0.4), fc=BAR, ec="none")
                bar.set_gid(f"cons:{c}")
                ax.add_patch(bar)

    for e in track.elements:
        if e.start not in col_of or e.end not in col_of:
            continue
        c0, c1 = col_of[e.start], col_of[e.end]
        colour = colour_of.get(e.id, "#000000")
        for piece, b in enumerate(range(c0 // columns, c1 // columns + 1)):
            lo, hi = max(c0, b * columns), min(c1, (b + 1) * columns - 1)
            x0, x1 = lo - b * columns - 0.5, hi - b * columns + 0.5
            gid = f"ss:{e.id}:{piece}"
            if e.kind == "E":
                _strand(ax, x0, x1, glyph_y(b), colour, hi == c1, gid, font)
            else:
                _coil(ax, x0, x1, glyph_y(b), colour, lo - c0, e.kind, gid, font)
            if lo == c0:
                t = ax.text(
                    x0,
                    top(b) - _LABEL_H * 0.45,
                    e.label,
                    ha="left",
                    va="center",
                    fontsize=font * 0.9,
                    color=darken(colour, 0.85) if colour != "#000000" else INK,
                )
                t.set_gid(f"sslabel:{e.id}")

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
                    c - b * columns, glyph_y(b), "TT", ha="center", va="center", fontsize=font * 0.75, color=INK
                )
                t.set_gid(f"turn:{k}")
            k = j + 1
        else:
            k += 1
    return fig


def _coil(ax, x0: float, x1: float, y: float, colour: str, offset: int, kind: str, gid: str, font: float) -> None:
    """A helix as a coiled line, one turn per 3.6 residues (3.0 for 3-10), phase continuous across rows."""
    pitch, amp, lw = (3.6, 0.42, 1.9) if kind == "H" else (3.0, 0.3, 1.4)
    xs = np.linspace(x0, x1, max(int((x1 - x0) * 14), 2))
    ys = y + amp * np.sin(2 * np.pi * (xs - x0 + offset) / pitch)
    line = PathPatch(
        Path(np.column_stack([xs, ys])),
        fc="none",
        ec=colour,
        lw=lw * font / 8,
        capstyle="round",
        joinstyle="round",
        zorder=2,
    )
    line.set_gid(gid)
    ax.add_patch(line)


def _strand(ax, x0: float, x1: float, y: float, colour: str, head: bool, gid: str, font: float) -> None:
    shaft, wing = 0.24, 0.52
    if head:
        h = min(1.1, (x1 - x0) * 0.6)
        pts = [
            (x0, y - shaft),
            (x1 - h, y - shaft),
            (x1 - h, y - wing),
            (x1, y),
            (x1 - h, y + wing),
            (x1 - h, y + shaft),
            (x0, y + shaft),
        ]
    else:
        pts = [(x0, y - shaft), (x1, y - shaft), (x1, y + shaft), (x0, y + shaft)]
    arrow = Polygon(pts, closed=True, fc=colour, ec=darken(colour, 0.7), lw=0.5 * font / 8, joinstyle="round", zorder=2)
    arrow.set_gid(gid)
    ax.add_patch(arrow)
