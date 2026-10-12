# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Draw a Layout and its loops with matplotlib; every artist carries a gid so SVG layers are editable."""

from __future__ import annotations

import logging
import re
import threading
from collections import Counter
from collections.abc import Callable
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, replace
from pathlib import Path as FsPath
from typing import NamedTuple

import numpy as np
from matplotlib import colormaps, rc_context
from matplotlib.axes import Axes
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path
from matplotlib.text import Text
from matplotlib.textpath import text_to_path
from matplotlib.transforms import Affine2D

from .hbonds import find_hbonds, parse_selection
from .layout import (
    DNA_MINOR,
    DNA_PITCH,
    DNA_W,
    END_LABEL,
    END_STUB,
    HEAD,
    IN_PER_UNIT,
    MAX_WIDTH_IN,
    NT_ROW,
    PITCH,
    SHAFT,
    Layout,
    Placed,
    PlacedStem,
    centred_box,
    dna_to_page,
    dna_y,
    domain_panel,
    font_pt,
    helix_label_pos,
    label_boxes,
    overlaps_any,
    resolve,
    stem_label_pos,
    termini,
)
from .model import SSE, group_by_chain
from .na_paths import midway
from .palette import chain_colors, darken, delta_e, text_color_on, tint
from .route import Loop
from .style import _ALIASES, CHOICES, MODALITIES, Style, backbone_overrides

_IN_PER_UNIT = IN_PER_UNIT
_MAX_WIDTH_IN = MAX_WIDTH_IN  # double-column width
_HEAD_HALF = HEAD / 2
_HEAD_LEN = 0.7
_CORNER = 0.3
_EXT = {".svg": "svg", ".pdf": "pdf", ".png": "png"}
_FONTS = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]  # journal sans first
_MEASURE = FontProperties(family=_FONTS, size=10)  # text widths are measured at this size and scaled
_RC = {"font.family": "sans-serif", "font.sans-serif": _FONTS}


class _QuietMissingFonts(logging.Filter):
    """The figure names the whole journal font list, so fonts a machine lacks are expected: matplotlib
    would warn about each one at every first draw. Other font warnings are left alone."""

    def filter(self, record: logging.LogRecord) -> bool:
        return not any(record.getMessage() == f"findfont: Font family {name!r} not found." for name in _FONTS)


logging.getLogger("matplotlib.font_manager").addFilter(_QuietMissingFonts())
_PANEL = "#eef1f4"
_LIPID, _HEADS, _HEAD_EDGE = "#f7f0de", "#e6d3a3", "#c2a765"  # membrane band, lipid head groups
_RUN_EDGE, _RUN_EDGE_W = "#3a3a3a", 1.8  # dark edge under residue-coloured loops (pale colours stay legible)
_SHORT_HELIX = 1.6  # helices/3-10 boxes shorter than this hold their residue numbers past the ends
LONG_LOOP, LONG_TAIL = 25, 10  # loops / chain ends at least this many residues get a residue-count label
_DOMAIN_TONES = ("#4c78a8", "#e45756", "#54a24b", "#b279a2", "#f58518", "#72b7b2")  # one hue per domain
MARK_DEFAULT = "#d1495b"  # a marked element without its own colour
MATE_GREY = "#c4c9cf"  # symmetry mates (or other chains) when one part is highlighted
_MATE_LINE = "#aab0b7"
_LEGEND_ROW = 1.0
_LOOP_ROOM = 0.25  # page units kept between a loop point and the page edge: the casing at the largest line width
STYLES = CHOICES["fill"]
LOOP_STYLES = CHOICES["loops"]


_font_pt = font_pt


_RES_PER_TURN = 3.6
_MIN_PITCH = 0.5  # page units; an end-on helix keeps legible turns instead of a dark smear
_BAND = 0.36  # ribbon width along the axis, as a fraction of one turn
_SHINE = 0.4  # highlight width, as a fraction of the ribbon


def _band(t0: float, t1: float, n: int = 12, taper: int = 0, width: float = _BAND) -> np.ndarray:
    """Ribbon band of the coil y = 0.5 sin(2 pi t) between t0 and t1 (axis units: one turn = 1, diameter = 1).
    The ribbon is a tape wound on the cylinder, so it is offset along the axis and turns edge-on at the top
    and bottom of each turn. taper=+1 narrows it to a point at t0, -1 at t1 (the helix ends on its axis),
    2 at both ends."""
    t = np.linspace(t0, t1, n)
    y = 0.5 * np.sin(2 * np.pi * t)
    half = np.full(n, width / 2)
    ramp = np.linspace(0.0, 1.0, n)
    if taper == 2:  # lens: narrow at both ends
        half *= np.sin(np.pi * ramp)
    elif taper:
        half *= np.sqrt(ramp if taper > 0 else ramp[::-1])
    return np.vstack([np.column_stack([t + half, y]), np.column_stack([t - half, y])[::-1]])


# One helical turn, drawn once and reused: the front band rises across the near face (y -1/2 -> 1/2),
# the back band falls across the far face. A whole helix is these two shapes shifted turn by turn, plus
# half front bands so the ribbon starts and ends on the axis, where the loops attach.
_FRONT = _band(-0.25, 0.25)
_BACK = _band(0.25, 0.75)
_FRONT_HEAD = _band(0.0, 0.25, 9, taper=+1)
_FRONT_TAIL = _band(-0.25, 0.0, 9, taper=-1)
_SHINE_BAND = _band(-0.16, 0.16, 14, taper=2, width=_BAND * _SHINE) + (_BAND * 0.1, 0)  # glint on the near face


def _coil(p: Placed, residues: int) -> tuple[Path, Path, Path]:
    """Front faces, back faces and front highlights of a helix as compound paths in page coordinates."""
    turns = int(max(1, min(round(residues / _RES_PER_TURN), p.length // _MIN_PITCH)))
    front = [_FRONT_HEAD] + [_FRONT + (k, 0) for k in range(1, turns)] + [_FRONT_TAIL + (turns, 0)]
    back = [_BACK + (k, 0) for k in range(turns)]
    shine = [_SHINE_BAND + (k, 0) for k in range(1, turns)]
    to_page = Affine2D().translate(-turns / 2, 0).scale(p.length / turns, p.width).rotate(p.angle).translate(p.cx, p.cy)

    def compound(polys):
        verts = np.vstack([np.vstack([q, q[:1]]) for q in polys])
        codes = np.concatenate([[Path.MOVETO] + [Path.LINETO] * (len(q) - 1) + [Path.CLOSEPOLY] for q in polys])
        return to_page.transform_path(Path(verts, codes))

    return compound(front), compound(back), compound(shine) if shine else None


_DNA = {
    "bold": ("#3d4852", "#aab4bd", "#c3cad1"),
    "pale": ("#6b7781", "#c5ccd2", "#d5dade"),
    "outline": ("#3d4852", "#d5dade", "#d5dade"),
}  # front backbone, back backbone, base-pair rungs
_DNA_DARK = ("#d3dce5", "#7d8b99", "#5a6b7c")  # the same on a dark page


def _dna_colours(style: str) -> tuple[str, str, str]:
    return _DNA_DARK if _luminance(_BACKGROUND.get()) < 0.5 else _DNA[style]


_CHAIN_POOL = (  # na_backbone=chain: colours for nucleic-acid chains, the one furthest from those already used first
    "#8c564b",
    "#6a3d9a",
    "#1b7f79",
    "#b15928",
    "#4d6b1f",
    "#a6275a",
    "#2b5c9e",
    "#7a6a00",
)


def backbone_colours(
    layout: Layout, protein: list[str], look: Style, colors: dict[str, str]
) -> tuple[dict[str, str], dict[str, str]]:
    """(nucleic-acid chain -> backbone colour, protein chain -> loop colour override) for this figure: the
    na_backbone mode, then `backbone_colors` (a modality key, then a chain key). Shared by every nucleic-acid
    style and the protein loops. Raises ValueError naming backbone_colors for a chain the figure lacks.
    Every colour goes through `_legible`: on a dark page a user-given colour that is too dark is lightened (on a
    light page one that is too pale is darkened); colours that already stand out are used exactly as given."""
    over = backbone_overrides(look.backbone_colors)
    strands = layout.nucleic.strands if layout.nucleic is not None else []
    rna = {s.chain: s.rna for s in reversed(strands)}  # a chain's first strand decides its modality
    known = {*protein, *rna, *layout.res_chain}
    for key in over:
        if key not in MODALITIES and key not in known:
            raise ValueError(f"backbone_colors: unknown chain {key!r}; this figure has {', '.join(sorted(known))}")
    mode = look.na_backbone
    taken = [look.ink, *(colors[c] for c in protein if c in colors)]
    na: dict[str, str] = {}
    for chain, is_rna in rna.items():
        if mode == "ink":
            colour = _dna_colours(look.fill)[0]
        elif mode == "chain":
            colour = max(_CHAIN_POOL, key=lambda c: min(delta_e(c, t) for t in taken))
            taken.append(colour)
        elif mode == "modality":
            colour = look.rna_color if is_rna else look.dna_color
        else:
            colour = mode
        colour = over.get("rna" if is_rna else "dna", colour)
        na[chain] = _legible(over.get(chain, colour))
    loops = {c: _legible(over[c] if c in over else over["protein"]) for c in protein if c in over or "protein" in over}
    return na, loops


_DNA_TUBE = 0.42  # backbone tube width, page units
_DNA_SEG = 16  # points along each half-turn


def _tube(pts: np.ndarray, width: float) -> np.ndarray:
    """A centre line thickened to a band of `width` (perpendicular offsets)."""
    t = np.gradient(pts, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    n = np.column_stack([-t[:, 1], t[:, 0]]) * width / 2
    return np.vstack([pts + n, (pts - n)[::-1]])


def _dna_turns(d: PlacedStem, slot: int) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Front and back half-turns of one backbone in duplex-local coordinates: the same half-turn fragment
    repeated every turn (shifted by the minor groove for the second backbone), clipped to the duplex ends."""
    x0, x1 = -d.length / 2, d.length / 2
    off = slot * DNA_MINOR * DNA_PITCH
    front, back = [], []
    for k in range(int(np.floor((x0 - off) / DNA_PITCH)) - 1, int(np.ceil((x1 - off) / DNA_PITCH)) + 1):
        for lo, hi, dest in ((-0.25, 0.25, front), (0.25, 0.75, back)):
            a, b = max(x0, off + (k + lo) * DNA_PITCH), min(x1, off + (k + hi) * DNA_PITCH)
            if b - a > 1e-6:
                xs = np.linspace(a, b, _DNA_SEG)
                dest.append(_tube(np.column_stack([xs, dna_y(slot, xs, d.width)]), _DNA_TUBE * d.width / DNA_W))
    return front, back


def _compound(polys: list[np.ndarray], closed: bool = True) -> Path:
    if not polys:
        return Path(np.empty((0, 2)))
    if closed:
        verts = np.vstack([np.vstack([q, q[:1]]) for q in polys])
        codes = np.concatenate([[Path.MOVETO] + [Path.LINETO] * (len(q) - 1) + [Path.CLOSEPOLY] for q in polys])
    else:
        verts = np.vstack(polys)
        codes = np.concatenate([[Path.MOVETO] + [Path.LINETO] * (len(q) - 1) for q in polys])
    return Path(verts, codes)


def _na_style(layout: Layout, look: Style) -> str:
    """The stem style drawn: a large nucleic acid cropped by na_extent auto goes from ladder to block."""
    return "block" if look.na_style == "ladder" and layout.na_cropped else look.na_style


class _NAInk(NamedTuple):
    """What every nucleic-acid drawing step shares."""

    ax: object
    front: str  # backbone and text colour
    back: str  # the far backbone of a double helix
    rung: str  # base pairs
    fill: str
    lw: float
    font: float
    gu: set  # (a, i, b, j) of the G·U pairs, drawn dotted
    chains: list  # chain of each strand
    backbone: dict  # chain -> backbone colour (the shared nucleic-acid backbone helper's answer)
    seqs: list  # per strand, the residue number of each nucleotide as the structure file writes it

    def colour(self, strand: int) -> str:
        """The backbone colour of a strand: stem rails, ribbons, blocks, paths, hairpins, ties."""
        return self.backbone[self.chains[strand]]


def _upright(angle: float) -> float:
    rot = float(np.degrees(angle))
    return rot - 180 if rot > 90 else rot + 180 if rot < -90 else rot


def _rung_lines(ink: _NAInk, k: int, d: PlacedStem, inset: float = 0.0) -> None:
    for pair in d.pairs:
        a, i, b, j = pair
        p, q = d.point((a, i)), d.point((b, j))
        v = (q - p) / max(float(np.hypot(*(q - p))), 1e-9)
        (line,) = ink.ax.plot(
            *np.c_[p + v * inset, q - v * inset],
            color=ink.rung,
            lw=ink.lw * 1.4,
            ls=(0, (1, 1.5)) if pair in ink.gu else "-",
            solid_capstyle="round",
            zorder=2.1,
        )
        line.set_gid(f"na:stem:{k}:rung")


def _ladder(ink: _NAInk, k: int, d: PlacedStem, where: dict, letters: dict) -> None:
    """Two rails (a bulged nucleotide pushed out a little) and a rung per pair; letters outside the rails."""
    rails: dict[str, list] = {}  # backbone colour -> its rails (one colour unless a duplex mixes chains)
    for sign in (1, -1):
        keys = sorted((key for key in d.axial if d.side[key] == sign), key=lambda key: d.axial[key])
        for colour in dict.fromkeys(ink.colour(key[0]) for key in keys):
            mine = [key for key in keys if ink.colour(key[0]) == colour]
            rails.setdefault(colour, []).append(
                np.array([d.point(key, 0.3 if key in d.bulges else 0.0) for key in mine])
            )
    for n, (colour, parts) in enumerate(rails.items()):
        rail = PathPatch(
            _compound([r for r in parts if len(r) > 1], closed=False),
            fc="none",
            ec=colour,
            lw=ink.lw * 2.2,
            capstyle="round",
            joinstyle="round",
            zorder=2.2,
        )
        rail.set_gid(f"na:stem:{k}:rails" if n == 0 else f"na:stem:{k}:rails:{n + 1}")
        ink.ax.add_patch(rail)
    _rung_lines(ink, k, d)
    for key in d.axial:
        letters[key] = (d.point(key, NT_ROW), _upright(d.angle))


def _draw_stem_label(ink: _NAInk, k: int, d: PlacedStem) -> None:
    t = ink.ax.text(
        *stem_label_pos(d),
        d.name,
        ha="center",
        va="center",
        rotation=_upright(d.angle),
        fontsize=ink.font * 0.8,
        fontweight="bold",
        color=ink.front,
        zorder=4,
    )
    t.set_gid(f"na:label:{k}")


def _ribbon_side(d: PlacedStem, key: tuple[int, int]) -> int:
    """Which backbone of a drawn double helix a nucleotide sits on: 0 for the first arm's, else 1."""
    return 0 if d.side[key] > 0 else 1


def _ribbon_point(d: PlacedStem, key: tuple[int, int]) -> np.ndarray:
    """Where a nucleotide sits on the helix backbone it belongs to."""
    x = d.axial[key]
    return dna_to_page(d, x, dna_y(_ribbon_side(d, key), x, d.width))


def _ribbon(ink: _NAInk, k: int, d: PlacedStem, where: dict, letters: dict) -> None:
    """The double helix, front and back backbones from one repeated turn, a rung per pair; no stem letters
    (the helix itself is the drawing). Beads sit on the backbone, in the strand's backbone colour."""
    colour = {}  # backbone slot -> colour of the strand drawn there
    for key in sorted(d.axial):
        colour.setdefault(_ribbon_side(d, key), ink.colour(key[0]))
    slots: dict[str, list[int]] = {}  # colour -> its backbones (one colour unless a duplex mixes chains)
    for slot, c in sorted(colour.items()):
        slots.setdefault(c, []).append(slot)

    def page(q: np.ndarray) -> np.ndarray:
        return dna_to_page(d, q[:, 0], q[:, 1])

    outline = ink.fill == "outline"
    for n, (c, group) in enumerate(slots.items()):
        fronts, backs = [], []
        for slot in group:
            f, b = _dna_turns(d, slot)
            fronts += f
            backs += b
        suffix = "" if n == 0 else f":{n + 1}"
        back = PathPatch(_compound([page(q) for q in backs]), fc=tint(c, 0.5), ec="none", zorder=2.0)
        back.set_gid(f"na:stem:{k}:back{suffix}")
        ink.ax.add_patch(back)
        front = PathPatch(
            _compound([page(q) for q in fronts]),
            fc="white" if outline else c,
            ec=c if outline else "none",
            lw=ink.lw * 0.8,
            zorder=2.2,
        )
        front.set_gid(f"na:stem:{k}:front{suffix}")
        ink.ax.add_patch(front)
    rungs = [
        dna_to_page(d, [x, x], [dna_y(0, x, d.width), dna_y(1, x, d.width)])
        for x in (d.axial[(a, i)] for a, i, _, _ in d.pairs)
    ]
    rung = PathPatch(
        _compound(rungs, closed=False), fc="none", ec=ink.rung, lw=ink.lw * 2.2, capstyle="round", zorder=2.1
    )
    rung.set_gid(f"na:stem:{k}:rungs")
    ink.ax.add_patch(rung)


def _block(ink: _NAInk, k: int, d: PlacedStem, where: dict, letters: dict) -> None:
    """A bar between the two rails, named inside (P1, DNA 1a ...); no letters. A duplex of two colours is two
    half bars."""
    colour = {}  # side -> colour of the strand on it
    for key in sorted(d.axial):
        colour.setdefault(d.side[key], ink.colour(key[0]))
    halves = [colour.get(1, colour.get(-1)), colour.get(-1, colour.get(1))]
    r = d.rail
    shown = [halves[0]] if halves[0] == halves[1] else halves
    for n, c in enumerate(shown):
        lo, hi = (-r, r) if len(shown) == 1 else ((0.0, r) if n == 0 else (-r, 0.0))
        corners = [
            dna_to_page(d, x, y)
            for x, y in ((-d.length / 2, lo), (d.length / 2, lo), (d.length / 2, hi), (-d.length / 2, hi))
        ]
        body = Polygon(
            corners, closed=True, fc="white" if ink.fill == "outline" else c, ec=c, lw=ink.lw * 1.2, zorder=2.2
        )
        body.set_gid(f"na:stem:{k}:block" if n == 0 else f"na:stem:{k}:block:{n + 1}")
        ink.ax.add_patch(body)
    mixed = _mix(shown[0], shown[-1], 0.5)
    t = ink.ax.text(
        d.cx,
        d.cy,
        d.name,
        ha="center",
        va="center",
        rotation=_upright(d.angle),
        fontsize=ink.font * 0.9,
        fontweight="bold",
        color=_legible(mixed) if ink.fill == "outline" else text_color_on(mixed),
        zorder=4,
    )
    t.set_gid(f"na:label:{k}")
    a, i, b, j = d.pairs[0]  # the first pair (outermost): its two residue numbers, in the slot beside the stem
    partner = f"{ink.chains[b]}:" if ink.chains[b] != ink.chains[a] else ""  # across chains, the partner's chain too
    t = ink.ax.text(
        *stem_label_pos(d),
        f"{ink.seqs[a][i]}\u00b7{partner}{ink.seqs[b][j]}",
        ha="center",
        va="center",
        rotation=_upright(d.angle),
        fontsize=ink.font * 0.75,  # a base letter's size: the stem's name slot, empty in this style
        color=ink.front,
        zorder=4,
    )
    t.set_gid(f"na:stem:{k}:first")


def _letters(ink: _NAInk, k: int, d: PlacedStem, where: dict, letters: dict) -> None:
    """The classic secondary-structure look: each nucleotide's letter on its rail, a dash per pair."""
    _rung_lines(ink, k, d, inset=0.3)
    for key in d.axial:
        letters[key] = (d.point(key), _upright(d.angle))


_STEM_STYLES = {"ladder": _ladder, "ribbon": _ribbon, "block": _block, "letters": _letters}  # na_style -> drawing


def _draw_strands(ink: _NAInk, layout: Layout, style: str, where: dict, letters: dict) -> None:
    """Single strands along their drawn paths (cropped gaps dashed, labelled unless short); their letters beside
    them, where routing placed them (on the strand itself in the letters style)."""
    for n, path in enumerate(layout.na_paths):
        pts = np.asarray(path.points, float)
        gap = path.kind == "gap"
        (line,) = ink.ax.plot(
            pts[:, 0],
            pts[:, 1],
            color=ink.colour(path.strand),
            lw=ink.lw * (1.0 if gap else 1.4),
            ls=((0, (2, 2)) if path.label else (0, (1, 1.5))) if gap else "-",
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=1.9,
        )
        line.set_gid(f"na:run:{n}")
        if path.label:  # where routing put it (_label_gaps), else standing on the gap's midpoint
            xy, va = (path.label_xy, "center") if path.label_xy is not None else (midway(path.points), "bottom")
            t = ink.ax.text(*xy, path.label, ha="center", va=va, fontsize=ink.font * 0.75, color=ink.front)
            t.set_gid(f"na:gap:{n}")
        for i, xy in zip(path.nts, path.nt_xy):
            key = (path.strand, i)
            letters[key] = (np.asarray(xy if style == "letters" else layout.na_letters[key]), 0.0)


_TIE_LABEL = 0.75  # a tie label's size: a base letter's (0.75 of the text size)
TIE_LEGEND = "a…a: tertiary base pair"  # the legend's words for a tie drawn as two labelled marks


@dataclass
class ContactBeads:
    """The contact marks of a drawn figure (left on the figure as `foldmap_beads`): chain -> [(nucleotide
    (strand, index), toucher residue (Backbone index), colour)] in drawing order, and the function that
    coloured them."""

    marks: dict
    colour: Callable[[int, str], str | None]


def _draw_na_marks(ink: _NAInk, layout: Layout, style: str, residue_colour, pt_per_unit, look, where, letters):
    """Isolated pairs as bars, base letters (touched ones in the toucher's colour), beads on the touched
    nucleotides in the protein chain's colour, and the strands' 5' and 3' ends."""
    na = layout.nucleic
    for tie in layout.na_ties:  # a short clear tie as a dotted line; a long or crossing one as two labelled marks
        gid = f"na:tie:{na.strands[tie.a[0]].labels[tie.a[1]].tag}-{na.strands[tie.b[0]].labels[tie.b[1]].tag}"
        if not tie.label:  # as the rungs: Watson-Crick solid, G·U dotted
            wobble = (*tie.a, *tie.b) in ink.gu
            (line,) = ink.ax.plot(
                *np.c_[tie.ends[0], tie.ends[1]],
                color=ink.rung,
                lw=ink.lw,
                ls=(0, (1, 1.5)) if wobble else "-",
                zorder=1.8,
            )
            line.set_gid(gid)
            continue
        for a, b, at in tie.marks:
            (tick,) = ink.ax.plot(
                *np.c_[a, b], color=ink.colour(tie.a[0]), lw=ink.lw * 1.2, solid_capstyle="round", zorder=2.4
            )
            tick.set_gid(f"{gid}:mark")
            t = ink.ax.text(
                *at,
                tie.label,
                ha="center",
                va="center",
                fontsize=ink.font * _TIE_LABEL,  # a base letter's size, upright, bold, in the mark's colour
                fontweight="bold",
                color=_legible(ink.colour(tie.a[0])),
                zorder=4,
            )
            t.set_gid(f"{gid}:label")
    touched: dict[tuple[int, int], tuple[str, str]] = {}  # nucleotide -> (chain, colour) of its first toucher
    toucher: dict[tuple[int, int], int] = {}  # nucleotide -> that protein residue (Backbone index)
    for r, nts in sorted(na.contacts.items()):
        chain = na.contact_chain.get(r)
        colour = residue_colour(r, chain)
        for key in sorted(nts):
            if colour is not None and key not in touched:
                touched[key] = (chain, colour)
                toucher[key] = r
    marks: dict[str, list] = {}  # chain -> [(nucleotide, toucher, colour)] as marked, for nacheck and the hover
    if look.nucleotide_labels and style != "block":
        for key, (xy, rot) in sorted(letters.items()):
            hit = touched.get(key)
            lab = na.strands[key[0]].labels[key[1]]
            if hit and style == "letters":  # the bold coloured letter is the mark
                marks.setdefault(hit[0], []).append((key, toucher[key], hit[1]))
            t = ink.ax.text(
                *xy,
                na.strands[key[0]].bases[key[1]],
                ha="center",
                va="center",
                rotation=rot,
                fontsize=ink.font * 0.75,
                fontweight="bold" if hit else "normal",
                color=hit[1] if hit else ink.colour(key[0]) if style == "letters" else ink.front,
                zorder=4,
            )
            t.set_gid(f"na:nt:{lab.chain}:{lab.tag}")
    if style != "letters":
        size = (_DNA_TUBE * 1.05 * pt_per_unit) ** 2
        for chain in sorted({c for c, _ in touched.values()}):
            keys = [key for key, (c, _) in sorted(touched.items()) if c == chain and key in where]
            if keys:
                marks[chain] = [(key, toucher[key], touched[key][1]) for key in keys]  # scatter order
                xy = np.array([where[key] for key in keys])
                sc = ink.ax.scatter(
                    xy[:, 0],
                    xy[:, 1],
                    s=size,
                    c=[touched[key][1] for key in keys],
                    edgecolors="white",
                    linewidths=ink.lw * 0.6,
                    zorder=2.3,
                )
                sc.set_gid(f"na:contact:{chain}")
    for s, end, (x, y) in layout.na_ends:
        strand = na.strands[s]
        cut = strand.cut5 if end == "5′" else strand.cut3
        t = ink.ax.text(
            x, y, "…" if cut else end, ha="center", va="center", fontsize=ink.font, fontweight="bold", color=ink.front
        )
        t.set_gid(f"na:end:{strand.id}:{end}")
    ink.ax.figure.foldmap_beads = ContactBeads(marks, residue_colour)


def _draw_nucleic(
    ax,
    layout: Layout,
    sses: list[SSE],
    residue_colour,
    fill: str,
    lw: float,
    font: float,
    pt_per_unit: float,
    look: Style,
    backbone: dict[str, str],
) -> None:
    """Stems in the chosen na_style, single strands along their drawn paths, base letters, contact beads in
    the touching protein chain's colour, isolated pairs as bars, cropped gaps and the 5'/3' ends."""
    na = layout.nucleic
    gu = {(p.a, p.i, p.b, p.j) for p in na.base_pairs if p.kind == "GU"}
    ink = _NAInk(
        ax,
        *_dna_colours(fill),
        fill,
        lw,
        font,
        gu,
        [s.chain for s in na.strands],
        backbone,
        [[lab.seq for lab in s.labels] for s in na.strands],
    )
    style = _na_style(layout, look)
    where = nucleotide_points(layout, style)  # nucleotide -> where it is drawn (bead and hydrogen bonds go there)
    letters: dict[tuple[int, int], tuple[np.ndarray, float]] = {}  # nucleotide -> letter position, rotation
    for k, d in enumerate(layout.stems):
        _STEM_STYLES.get(style, _ladder)(ink, k, d, where, letters)
        if style != "block" and look.labels:
            _draw_stem_label(ink, k, d)
    _draw_strands(ink, layout, style, where, letters)
    _draw_na_marks(ink, layout, style, residue_colour, pt_per_unit, look, where, letters)


def _arrow(p: Placed) -> np.ndarray:
    half, head = p.length / 2, min(_HEAD_LEN, p.length * 0.45)
    k = p.width / PITCH  # Style.strand_scale, carried by the strand's width
    sh, hh = SHAFT * k / 2, _HEAD_HALF * k
    local = np.array(
        [
            (-half, -sh),
            (half - head, -sh),
            (half - head, -hh),
            (half, 0),
            (half - head, hh),
            (half - head, sh),
            (-half, sh),
        ]
    )
    c, s = np.cos(p.angle), np.sin(p.angle)
    return local @ np.array([[c, s], [-s, c]]) + [p.cx, p.cy]


def _polyline(points: list[tuple[float, float]]) -> Path:
    return Path(np.array(points, float), [Path.MOVETO] + [Path.LINETO] * (len(points) - 1))


def _rounded_path(points: list[tuple[float, float]], radius: float = _CORNER) -> Path:
    pts = [np.array(p) for p in points]
    verts, codes = [pts[0]], [Path.MOVETO]
    for a, b, c in zip(pts, pts[1:], pts[2:]):
        r = min(radius, np.linalg.norm(b - a) / 2, np.linalg.norm(c - b) / 2)
        u1, u2 = (b - a) / np.linalg.norm(b - a), (c - b) / np.linalg.norm(c - b)
        verts += [b - u1 * r, b, b + u2 * r]
        codes += [Path.LINETO, Path.CURVE3, Path.CURVE3]
    verts.append(pts[-1])
    codes.append(Path.LINETO)
    return Path(np.array(verts), codes)


_CURVE_MAX = 3.0  # longest stretch of a leg a curved corner may take
_T = np.linspace(0.0, 1.0, 25)[:, None]


def _corner_clear(ctrl: list[np.ndarray], rects: np.ndarray) -> bool:
    p0, p1, p2, p3 = ctrl
    pts = (1 - _T) ** 3 * p0 + 3 * (1 - _T) ** 2 * _T * p1 + 3 * (1 - _T) * _T**2 * p2 + _T**3 * p3
    x, y = pts[:, 0:1], pts[:, 1:2]
    inside = (x > rects[:, 0]) & (x < rects[:, 2]) & (y > rects[:, 1]) & (y < rects[:, 3])
    return not inside.any()


def _smooth_path(points: list[tuple[float, float]], rects: np.ndarray) -> Path:
    """Cubic curves in place of the routed path's corners, so the loop keeps the route's way around
    obstacles but loses the right angles. Each corner takes as much of its legs as it can (up to half,
    at most _CURVE_MAX) while the curve stays out of every element box in `rects`."""
    pts = [np.array(p, float) for p in points]
    if len(pts) < 3:
        return Path(np.array(pts), [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1))
    verts, codes = [pts[0]], [Path.MOVETO]
    for k in range(1, len(pts) - 1):
        a, c, b = pts[k - 1], pts[k], pts[k + 1]
        la, lb = np.linalg.norm(c - a), np.linalg.norm(b - c)
        ua, ub = (c - a) / la, (b - c) / lb
        r = min(la / 2, lb / 2, _CURVE_MAX)
        while True:
            s, e = c - ua * r, c + ub * r
            ctrl = [s, s + (c - s) * 2 / 3, e + (c - e) * 2 / 3, e]
            if r < 0.02 or _corner_clear(ctrl, rects):
                break
            r *= 0.6
        verts += ctrl
        codes += [Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    verts.append(pts[-1])
    codes.append(Path.LINETO)
    return Path(np.array(verts), codes)


def _chains_in_order(sses: list[SSE]) -> list[str]:
    seen: list[str] = []
    for s in sses:
        if s.chain not in seen:
            seen.append(s.chain)
    return seen


_RAMP = {"turbo": (0.1, 0.92), "viridis": (0.0, 0.92), "plasma": (0.0, 0.9), "cividis": (0.0, 1.0)}  # skip muddy ends


def sequence_colour(t: float, name: str) -> str:
    """Colour at fraction t (0 = N end, 1 = C end) of the trimmed N -> C ramp."""
    lo, hi = _RAMP[name]
    return to_hex(colormaps[name](lo + float(np.clip(t, 0.0, 1.0)) * (hi - lo)))


def _shade(colour: str, t: float) -> str:
    """The chain's hue from light (t=0, N end) through itself (t=0.5) to dark (t=1, C end)."""
    return tint(colour, 0.55 * (1 - 2 * t)) if t < 0.5 else darken(colour, 1 - 0.4 * (2 * t - 1))


# The page being drawn (set by page()), so line colours can adapt to it. A context variable, not a global:
# concurrent app sessions draw in their own threads.
_BACKGROUND: ContextVar[str] = ContextVar("foldmap_background", default="#ffffff")


# rc_context edits the process-wide rcParams, and matplotlib is not thread-safe, so pages are drawn one at a time.
_DRAWING = threading.RLock()


def pin_fonts(fig: Figure) -> Figure:
    """Name foldmap's fonts on every text that follows the default family, so the figure looks the same
    however it is saved (plain fig.savefig, a notebook) and not only inside page()."""
    for t in fig.findobj(Text):
        if t.get_fontfamily() == ["sans-serif"]:
            t.set_fontfamily(_FONTS)
    return fig


@contextmanager
def page(look: Style):
    """Draw on `look`'s page: its background for colour choices and its ink and fonts for text."""
    with _DRAWING:
        token = _BACKGROUND.set(look.background)
        try:
            with rc_context({**_RC, "text.color": look.ink}):
                yield
        finally:
            _BACKGROUND.reset(token)


def _luminance(colour) -> float:
    from matplotlib.colors import to_rgb

    r, g, b = to_rgb(colour)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _mix(a: str, b: str, t: float) -> str:
    """a blended toward b by t (0 = a, 1 = b)."""
    from matplotlib.colors import to_rgb

    return to_hex([x + (y - x) * t for x, y in zip(to_rgb(a), to_rgb(b))])


def _label_ink(colour: str) -> str:
    """Text in an element's own hue: a deeper shade on light pages, the colour itself on dark ones."""
    if _luminance(_BACKGROUND.get()) < 0.5:
        return _legible(colour)
    return _legible(darken(colour, 0.8))


def _legible(colour: str, ceiling: float = 0.5) -> str:
    """Shift a line colour until it stands out on the page: darker on light pages (pale shades are fine for
    fills, not for thin lines), lighter on dark ones."""
    from matplotlib.colors import to_rgb

    if _luminance(_BACKGROUND.get()) < 0.5:
        for _ in range(12):
            if _luminance(colour) >= 0.55:
                break
            colour = _mix(to_hex(to_rgb(colour)), "#ffffff", 0.2)
        return colour
    for _ in range(12):
        r, g, b = to_rgb(colour)
        if 0.2126 * r + 0.7152 * g + 0.0722 * b <= ceiling:
            break
        colour = darken(to_hex((r, g, b)), 0.85)
    return to_hex(colour)


_SULFUR = "#e3b505"
_SNFG = {  # Symbol Nomenclature for Glycans: (shape, colour)
    "NAG": ("square", "#0090bc"),
    "NDG": ("square", "#0090bc"),
    "GCS": ("square", "#0090bc"),
    "A2G": ("square", "#ffd400"),
    "NGA": ("square", "#ffd400"),
    "MAN": ("circle", "#00a651"),
    "BMA": ("circle", "#00a651"),
    "GAL": ("circle", "#ffd400"),
    "GLA": ("circle", "#ffd400"),
    "GLB": ("circle", "#ffd400"),
    "GLC": ("circle", "#0090bc"),
    "BGC": ("circle", "#0090bc"),
    "FUC": ("triangle", "#ed1c24"),
    "FUL": ("triangle", "#ed1c24"),
    "SIA": ("diamond", "#a54399"),
    "SLB": ("diamond", "#a54399"),
    "NGC": ("diamond", "#8fcce9"),
    "KDN": ("diamond", "#00a651"),
    "XYS": ("star", "#f47920"),
    "XYP": ("star", "#f47920"),
}


def residue_point(layout: Layout, loops: list[Loop], sses: list[SSE], bb, r: int) -> tuple[float, float] | None:
    """Where residue r (Backbone index) sits on the page: along its element, or along the loop it is in, in
    proportion to its place in the sequence; past a chain end, at the end of the terminal stub."""
    chains = [l.chain for l in bb.labels] if hasattr(bb, "labels") else bb
    chain = chains[r]
    mine = [s for s in sses if s.chain == chain and s.id in layout.placed]
    for s in mine:
        if s.start <= r <= s.end:
            p = layout.placed[s.id]
            t = (r - s.start) / max(s.end - s.start, 1)
            return tuple(float(v) for v in np.add(p.n_port, t * np.subtract(p.c_port, p.n_port)))
    before = [s for s in mine if s.end < r]
    after = [s for s in mine if s.start > r]
    if before and after:
        a, b = before[-1], after[0]
        loop = next((l for l in loops if l.a_id == a.id and l.b_id == b.id), None)
        if loop is not None:
            pts = np.asarray(loop.points, float)
            seg = np.hypot(*np.diff(pts, axis=0).T)
            target = (r - a.end) / max(b.start - a.end, 1) * seg.sum()
            k = int(np.searchsorted(np.cumsum(seg), target))
            k = min(k, len(seg) - 1)
            into = (target - (np.cumsum(seg)[k] - seg[k])) / max(seg[k], 1e-9)
            return tuple(float(v) for v in pts[k] + into * (pts[k + 1] - pts[k]))
    if before:
        p = layout.placed[before[-1].id]
        return tuple(float(v) for v in np.add(p.c_port, p.exit_c * END_STUB))
    if after:
        p = layout.placed[after[0].id]
        return tuple(float(v) for v in np.add(p.n_port, p.exit_n * END_STUB))
    return None


def loop_of(layout: Layout, loops: list[Loop], sses: list[SSE], chains, r: int) -> np.ndarray | None:
    """The routed points of the loop residue r lies in (None if it sits in an element or a chain end)."""
    chain = chains[r]
    mine = [s for s in sses if s.chain == chain and s.id in layout.placed]
    if any(s.start <= r <= s.end for s in mine):
        return None
    before = [s for s in mine if s.end < r]
    after = [s for s in mine if s.start > r]
    if not (before and after):
        return None
    loop = next((l for l in loops if l.a_id == before[-1].id and l.b_id == after[0].id), None)
    return None if loop is None else np.asarray(loop.points, float)


def _nearest_on(pts: np.ndarray, q) -> tuple[float, float]:
    """The point of polyline pts closest to q, kept clear of the loop's two end stubs."""
    q = np.asarray(q, float)
    best, d_best = None, np.inf
    for a, b in zip(pts, pts[1:]):
        ab = b - a
        t = float(np.clip(np.dot(q - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0))
        c = a + t * ab
        d = float(np.hypot(*(c - q)))
        if d < d_best:
            best, d_best = c, d
    return float(best[0]), float(best[1])


def disulfide_points(layout, loops, sses, chains, i: int, j: int):
    """Where to draw a disulfide's two cysteines. A cysteine in a loop is placed where its loop passes closest
    to its partner: loop routes are schematic, so any point on its own loop is faithful, and the nearest one
    keeps the bond short. Cysteines in elements stay at their place along the element."""
    a, b = residue_point(layout, loops, sses, chains, i), residue_point(layout, loops, sses, chains, j)
    if a is None or b is None:
        return None, None
    la, lb = loop_of(layout, loops, sses, chains, i), loop_of(layout, loops, sses, chains, j)
    for _ in range(3):  # both in loops: settle by alternating
        if la is not None:
            a = _nearest_on(la, b)
        if lb is not None:
            b = _nearest_on(lb, a)
    if np.hypot(a[0] - b[0], a[1] - b[1]) < _SS_MIN and (la is not None or lb is not None):
        b = _nearest_beyond(lb, a) if lb is not None else b  # loops crossing at the bond: keep a visible bar
        if np.hypot(a[0] - b[0], a[1] - b[1]) < _SS_MIN and la is not None:
            a = _nearest_beyond(la, b)
    return a, b


_SS_MIN = 0.5  # shortest disulfide bar drawn (page units)
_SS_PARALLEL = np.cos(np.radians(10))  # bars closer in direction than this can lie on top of each other
_SS_CLEAR = 2.0  # bars nearer than this many bar widths (side by side) cannot be told apart
_SS_STEP = 3.0  # lane spacing of a detouring bar: at least this many bar widths...
_SS_DOT_STEP = 2.5  # ...and this many sulfur-dot radii, so the legs always stand out from a dot


def disulfide_step(lw: float, pt_per_unit: float) -> float:
    """Spacing between lanes of detouring bars, page units (bar and dot sizes as _draw_links draws them)."""
    bar, dot_radius = lw * 2.4 / pt_per_unit, lw * 5.5 / 2 / pt_per_unit
    return max(_SS_STEP * bar, _SS_DOT_STEP * dot_radius)


def _unit(a, b) -> np.ndarray:
    d = np.asarray(b, float) - np.asarray(a, float)
    return d / (np.hypot(*d) or 1.0)


def _bars_collide(s, t, width: float) -> bool:
    """Two bars (point pairs) nearly parallel, closer than _SS_CLEAR bar widths and overlapping along their length."""
    u, v = _unit(*s), _unit(*t)
    if abs(u @ v) < _SS_PARALLEL:
        return False
    n = np.array([-u[1], u[0]])
    if np.mean([abs((np.asarray(p, float) - np.asarray(s[0], float)) @ n) for p in t]) >= _SS_CLEAR * width:
        return False
    along_s = sorted(float((np.asarray(p, float) - np.asarray(s[0], float)) @ u) for p in s)
    along_t = sorted(float((np.asarray(p, float) - np.asarray(s[0], float)) @ u) for p in t)
    return min(along_s[1], along_t[1]) > max(along_s[0], along_t[0])


def disulfide_lanes(segments, width: float) -> list[int]:
    """A lane for each bar: 0 stays straight; bars that would lie on top of one another get different lanes, by
    interval colouring in order of where they start (so the result does not depend on the order given)."""
    n = len(segments)
    hit = [[_bars_collide(segments[i], segments[j], width) if i != j else False for j in range(n)] for i in range(n)]
    lanes = [0] * n
    done: set[int] = set()
    for i in range(n):  # one conflict group at a time, coloured along the group's first bar
        if i in done:
            continue
        group, todo = [], [i]
        while todo:
            k = todo.pop()
            if k not in done:
                done.add(k)
                group.append(k)
                todo.extend(j for j in range(n) if hit[k][j] and j not in done)
        origin, u = np.asarray(segments[min(group)][0], float), _unit(*segments[min(group)])
        start = {k: min(float((np.asarray(p, float) - origin) @ u) for p in segments[k]) for k in group}
        placed: list[int] = []
        for k in sorted(group, key=lambda k: (start[k], k)):
            used = {lanes[j] for j in placed if hit[k][j]}
            lanes[k] = next(lane for lane in range(n) if lane not in used)
            placed.append(k)
    return lanes


def disulfide_route(a, b, lane: int, side: int, step: float) -> list[tuple[float, float]]:
    """The points of a bar: straight in lane 0, else out from each cysteine, along, and back (side is +/-1)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if lane == 0:
        return [tuple(a), tuple(b)]
    u = _unit(a, b)
    off = side * lane * step * np.array([-u[1], u[0]])
    return [tuple(a), tuple(a + off), tuple(b + off), tuple(b)]


def disulfide_side(a, b, lane: int, taken, width: float, step: float) -> int:
    """Which side (+1 or -1 of the bar's left normal) a detour leaves on: the one with less overlap against the
    boxes in taken; a tie goes below, or to the right for a vertical bar."""
    u = _unit(a, b)
    n = np.array([-u[1], u[0]])
    if abs(n[1]) > 1e-9:
        below = 1 if n[1] < 0 else -1  # the sign whose normal points down
    else:
        below = 1 if n[0] > 0 else -1  # vertical bar: the sign whose normal points right

    def covered(side: int) -> float:
        pts = disulfide_route(a, b, lane, side, step)
        total, h = 0.0, width / 2
        for p, q in zip(pts, pts[1:]):
            x0, y0, x1, y1 = min(p[0], q[0]) - h, min(p[1], q[1]) - h, max(p[0], q[0]) + h, max(p[1], q[1]) + h
            total += sum(
                max(0.0, min(x1, bx1) - max(x0, bx0)) * max(0.0, min(y1, by1) - max(y0, by0))
                for bx0, by0, bx1, by1 in taken
            )
        return total

    less = covered(1) - covered(-1)
    return below if abs(less) < 1e-9 else (1 if less < 0 else -1)


def _nearest_beyond(pts: np.ndarray, q) -> tuple[float, float]:
    """The point of polyline pts closest to q but at least _SS_MIN from it."""
    dense = np.vstack([a + t * (b - a) for a, b in zip(pts, pts[1:]) for t in np.linspace(0, 1, 30)])
    d = np.hypot(*(dense - np.asarray(q, float)).T)
    ok = d >= _SS_MIN
    k = int(np.argmin(np.where(ok, d, np.inf))) if ok.any() else int(np.argmax(d))
    return float(dense[k][0]), float(dense[k][1])


def _symbol(ax, shape: str, colour: str, xy, size: float, lw: float, gid: str) -> None:
    x, y = xy
    h = size / 2
    if shape == "square":
        patch = Rectangle((x - h, y - h), size, size)
    elif shape == "circle":
        from matplotlib.patches import Circle

        patch = Circle((x, y), h)
    else:
        n, r0 = {"triangle": (3, 1.15), "diamond": (4, 1.2), "star": (5, 1.15)}[shape]
        ang = np.pi / 2 + np.arange(n) * 2 * np.pi / n
        patch = Polygon(np.column_stack([x + r0 * h * np.cos(ang), y + r0 * h * np.sin(ang)]), closed=True)
    patch.set(facecolor=colour, edgecolor=_legible("#2b2b2b"), linewidth=lw * 0.6, zorder=3.8)
    patch.set_gid(gid)
    ax.add_patch(patch)


METAL_COLOURS = {  # Jmol-like element colours, darkened enough to carry a white symbol
    "Zn": "#6f7bb0",
    "Fe": "#d0602c",
    "Mg": "#4f9a3a",
    "Ca": "#3d8a3d",
    "K": "#8f40d4",
    "Na": "#7b4fc9",
    "Cu": "#b8752e",
    "Mn": "#8a63b8",
    "Co": "#c75c7a",
    "Ni": "#3f9f6a",
    "Cd": "#b08a2e",
    "Hg": "#7a7a92",
}
LIGAND_FILL, LIGAND_EDGE = "#fff1c1", "#b8860b"
_LIGAND_TETHERS = 4  # a ligand is tethered to the elements it touches most, at most this many
_TETHER_REACH = 6.5  # page units: a tether target farther than this from the heaviest one is another pocket
_SPOT_RADIUS = 3.0  # the spot search stops at this ring radius
# A tether to another chain's residue is at most _TETHER_REACH from the home chain's heaviest target, and the marker
# is within _SPOT_RADIUS of the anchor (an average of the targets): TETHER_CAP bounds the cross-chain tethers
# measured from the marker when the home targets are close to the heaviest.
TETHER_CAP = _TETHER_REACH + _SPOT_RADIUS


def _home_chain(layout: Layout, g) -> str:
    """The chain a ligand is drawn with: its own when it touches it, else the one it touches most (ties go to the
    alphabetically first)."""
    count: dict[str, int] = {}
    for r in g.contacts:
        count[layout.res_chain[r]] = count.get(layout.res_chain[r], 0) + 1
    if g.chain in count:
        return g.chain
    return max(sorted(count), key=count.__getitem__)


def _tether_targets(layout: Layout, loops: list[Loop], sses: list[SSE], g, chain: str | None):
    """Page points the ligand is tethered to, their weights, and their residues; only `chain`'s contacts when given."""
    groups: dict[str, list[int]] = {}  # contacts by the element (or loop) they belong to
    for r in g.contacts:
        if chain is not None and layout.res_chain[r] != chain:
            continue
        home = next((s.id for s in sses if s.start <= r <= s.end and s.id in layout.placed), None)
        if home is None:
            home = "loop:" + str(max((s.end for s in sses if s.end < r), default=-1))
        groups.setdefault(home, []).append(r)
    ranked = sorted(groups.values(), key=len, reverse=True)
    if not g.metal:  # a ligand: one tether per element it sits on, the four it touches most
        ranked = [[sorted(rs)[len(rs) // 2]] for rs in ranked[:_LIGAND_TETHERS]]
    targets, weights, residues = [], [], []
    for rs in ranked:
        for r in rs:
            at = residue_point(layout, loops, sses, layout.res_chain, r)
            if at is not None and all(np.hypot(at[0] - t[0], at[1] - t[1]) > 0.3 for t in targets):
                targets.append(at)
                weights.append(len(groups.get(next(k for k, v in groups.items() if r in v), [])))
                residues.append(r)
    return targets, weights, residues


def _local_targets(layout: Layout, loops: list[Loop], sses: list[SSE], g, metal: bool):
    """The tether targets of a ligand: all of its home chain's, then those of other chains that lie within
    _TETHER_REACH of the home chain's heaviest target (the heaviest of the others when the home chain is not drawn,
    and then that chain is the home). Never fewer than the home chain's."""
    home = _home_chain(layout, g)
    targets, weights, residues = _tether_targets(layout, loops, sses, g, home)
    all_t, all_w, all_r = _tether_targets(layout, loops, sses, g, None)
    if not targets:
        if not all_t:
            return [], [], []
        home = layout.res_chain[all_r[int(np.argmax(all_w))]]
        targets, weights, residues = _tether_targets(layout, loops, sses, g, home)
    ref = targets[int(np.argmax(weights))]  # the heaviest; ties go to the first
    for t, w, r in zip(all_t, all_w, all_r):
        if layout.res_chain[r] == home or (not metal and len(targets) >= _LIGAND_TETHERS):
            continue
        if np.hypot(t[0] - ref[0], t[1] - ref[1]) <= _TETHER_REACH and all(
            np.hypot(t[0] - u[0], t[1] - u[1]) > 0.3 for u in targets
        ):
            targets.append(t)
            weights.append(w)
            residues.append(r)
    return targets, weights, residues


def ligand_marks(layout: Layout, loops: list[Loop], sses: list[SSE], look: Style) -> list[dict]:
    """Where each shown ligand/ion goes: as close as possible to the residues holding it, clear of elements,
    labels and other markers. Each mark: ligand, centre, half size, tether targets, tethered residues."""
    links = layout.links
    if links is None or look.ligands == "none" or not getattr(links, "ligands", None):
        return []
    codes = None if look.ligands in ("auto", "all") else {c.upper() for c in look.ligands.split(",")}
    chosen = [
        g
        for g in links.ligands
        if (codes is not None and {c.upper() for c in (g.codes or [g.name])} & codes)
        or (codes is None and (look.ligands == "all" or not g.additive))
    ]
    if layout.focus_chains:  # one subunit of an assembly: only what binds it
        chosen = [g for g in chosen if any(layout.res_chain[r] in layout.focus_chains for r in g.contacts)]
    taken = [p.rect for p in layout.placed.values()] + [g.rect for g in layout.ghosts]
    taken += [box for _, box in label_boxes(layout, sses)]
    out = []
    for g in sorted(chosen, key=lambda g: g.contacts[0]):
        targets, weights, tethered = _local_targets(layout, loops, sses, g, g.metal)
        if not targets:
            continue
        anchor = np.average(np.asarray(targets), axis=0, weights=weights)
        half = (0.3, 0.3) if g.metal else (0.12 * len(g.name) + 0.2, 0.27)
        spot = None
        for ring in np.arange(0.0, _SPOT_RADIUS, 0.35):
            for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False):
                c = anchor + ring * np.array([np.cos(angle), np.sin(angle)])
                box = (c[0] - half[0] - 0.08, c[1] - half[1] - 0.08, c[0] + half[0] + 0.08, c[1] + half[1] + 0.08)
                if not any(box[0] < t[2] and t[0] < box[2] and box[1] < t[3] and t[1] < box[3] for t in taken):
                    spot = c
                    break
            if spot is not None:
                break
        if spot is None:
            continue
        taken.append((spot[0] - half[0], spot[1] - half[1], spot[0] + half[0], spot[1] + half[1]))
        out.append(
            {
                "ligand": g,
                "centre": (float(spot[0]), float(spot[1])),
                "half": half,
                "targets": targets,
                "residues": tethered,
            }
        )
    return out


def _draw_ligands(ax, marks: list[dict], lw: float, font: float) -> None:
    for m in marks:
        g, (cx, cy), (hw, hh) = m["ligand"], m["centre"], m["half"]
        key = f"{g.name}:{g.chain}{g.seq}"
        for k, (tx, ty) in enumerate(m["targets"]):
            tether = PathPatch(
                Path([(cx, cy), (tx, ty)]),
                fc="none",
                ec=_legible("#6b6b6b"),
                lw=lw * 0.7,
                ls=(0, (1.2, 1.4)),
                capstyle="round",
                zorder=3.6,
            )
            tether.set_gid(f"ligand-tether:{key}:{k}")
            ax.add_patch(tether)
        if g.metal:
            fill = METAL_COLOURS.get(g.symbol, "#7f7f7f")
            marker = Circle((cx, cy), hw, fc=fill, ec=darken(fill, 0.7), lw=lw * 0.7, zorder=3.8)
            ink = text_color_on(fill)
        else:
            marker = FancyBboxPatch(
                (cx - hw, cy - hh),
                2 * hw,
                2 * hh,
                boxstyle="round,pad=0,rounding_size=0.12",
                fc=LIGAND_FILL,
                ec=LIGAND_EDGE,
                lw=lw * 0.7,
                zorder=3.8,
            )
            ink = "#5c4300"
        marker.set_gid(f"ligand:{key}")
        ax.add_patch(marker)
        t = ax.text(
            cx,
            cy,
            g.symbol,
            ha="center",
            va="center",
            fontsize=font * (0.62 if g.metal else 0.6),
            fontweight="bold",
            color=ink,
            zorder=3.9,
        )
        t.set_gid(f"ligand-label:{key}")


_LINKED_STAR = 0.38  # radius of the mark of a residue fused into the chain
_LINKED_LABEL = 0.6  # its label, as a fraction of the font
_LINKED_LABEL_MIN = 7.0  # but never smaller than this (pt)


def linked_marks(layout: Layout, loops: list[Loop], sses: list[SSE], marks: list[dict]) -> list[dict]:
    """A star on the chain at each peptide-linked non-standard residue (a chromophore: GFP's CRO66) of a drawn
    chain, and its code and number beside it, in the first free spot on rings around it. Each mark: residue
    label, star centre, label centre, label half size."""
    shown = [k for k in layout.res_linked if layout.res_chain[k] in {s.chain for s in sses if s.id in layout.placed}]
    if layout.focus_chains:
        shown = [k for k in shown if layout.res_chain[k] in layout.focus_chains]
    if not shown:
        return []
    taken = [p.rect for p in layout.placed.values()] + [g.rect for g in layout.ghosts]
    taken += [box for _, box in label_boxes(layout, sses)]
    taken += [
        (
            m["centre"][0] - m["half"][0],
            m["centre"][1] - m["half"][1],
            m["centre"][0] + m["half"][0],
            m["centre"][1] + m["half"][1],
        )
        for m in marks
    ]
    out = []
    for k in shown:
        at = residue_point(layout, loops, sses, layout.res_chain, k)
        if at is None:
            continue
        label = layout.res_labels[k]
        text = f"{label.name}{label.tag}"
        half = (0.07 * len(text) + 0.1, 0.17)
        taken.append((at[0] - _LINKED_STAR, at[1] - _LINKED_STAR, at[0] + _LINKED_STAR, at[1] + _LINKED_STAR))
        spot = None
        for ring in np.arange(_LINKED_STAR + 0.1, 12.0, 0.35):
            for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False):
                c = np.asarray(at) + ring * np.array([np.cos(angle), np.sin(angle)])
                box = (c[0] - half[0] - 0.12, c[1] - half[1] - 0.05, c[0] + half[0] + 0.12, c[1] + half[1] + 0.05)
                if not any(box[0] < t[2] and t[0] < box[2] and box[1] < t[3] and t[1] < box[3] for t in taken):
                    spot = c
                    break
            if spot is not None:
                break
        if spot is not None:
            taken.append((spot[0] - half[0], spot[1] - half[1], spot[0] + half[0], spot[1] + half[1]))
        out.append(
            {
                "label": label,
                "text": text,
                "centre": (float(at[0]), float(at[1])),
                "text_at": None if spot is None else (float(spot[0]), float(spot[1])),
                "half": half,
            }
        )
    return out


def _star(cx: float, cy: float, r: float) -> np.ndarray:
    """Five points out, notches in at 0.4 of the radius."""
    angles = np.pi / 2 + np.arange(10) * np.pi / 5
    radii = np.where(np.arange(10) % 2 == 0, r, 0.4 * r)
    return np.c_[cx + radii * np.cos(angles), cy + radii * np.sin(angles)]


def _draw_linked(ax, linked: list[dict], lw: float, font: float) -> None:
    for m in linked:
        cx, cy = m["centre"]
        key = f"{m['label'].chain}:{m['label'].tag}"
        star = Polygon(
            _star(cx, cy, _LINKED_STAR),
            closed=True,
            fc=LIGAND_FILL,
            ec=LIGAND_EDGE,
            lw=lw * 0.7,
            joinstyle="round",
            zorder=3.8,
        )
        star.set_gid(f"linked:{key}")
        ax.add_patch(star)
        if m["text_at"] is not None:
            t = ax.text(
                *m["text_at"],
                m["text"],
                ha="center",
                va="center",
                fontsize=max(_LINKED_LABEL_MIN, font * _LINKED_LABEL),
                color="#5c4300",
                zorder=3.9,
            )
            t.set_gid(f"linked-label:{key}")


SITE_COLOURS = {  # UniProt site types, Okabe-Ito based
    "Active site": "#d55e00",
    "Binding site": "#0072b2",
    "Site": "#009e73",
    "Modified residue": "#cc79a7",
    "Lipidation": "#e69f00",
    "Cross-link": "#56b4e9",
    "Glycosylation": "#7f7f7f",
    "Disulfide bond": "#f0c419",
}
SITE_GROUPS = {
    "active": "Active site",
    "binding": "Binding site",
    "site": "Site",
    "modified": "Modified residue",
    "lipidation": "Lipidation",
    "cross-link": "Cross-link",
}


def uniprot_marks(layout, loops, sses, look) -> list[dict]:
    """UniProt sites to draw: (feature, chain, residue, page point). Glycosylation and disulfides come from the
    structure itself, so they are not repeated here."""
    ann = getattr(layout, "uniprot", None)
    if ann is None or look.uniprot_sites == "none":
        return []
    if look.uniprot_sites == "auto":
        wanted = {"Active site", "Binding site"}
    elif look.uniprot_sites == "all":
        wanted = set(SITE_GROUPS.values())
    else:
        wanted = {SITE_GROUPS[t] for t in look.uniprot_sites.split(",") if t in SITE_GROUPS}
    shown = lambda c: not layout.partial or c not in layout.partial  # noqa: E731
    out, seen = [], set()
    for f, chain, ks in ann.sites:
        if f.type not in wanted or not shown(chain):
            continue
        for k in ks:
            if (f.type, k) in seen:
                continue
            seen.add((f.type, k))
            at = residue_point(layout, loops, sses, layout.res_chain, k)
            if at is not None:
                out.append({"feature": f, "chain": chain, "residue": k, "at": at})
    return out


def _draw_uniprot(ax, marks: list[dict], layout, lw: float, font: float) -> None:
    for n, m in enumerate(marks):
        f, (x, y) = m["feature"], m["at"]
        colour = SITE_COLOURS.get(f.type, "#555555")
        halo = Circle((x, y), 0.27, fc=_BACKGROUND.get(), ec=_legible("#1d1d1f"), lw=lw * 0.9, zorder=3.84)
        halo.set_gid(f"feature-ring:{n}")
        ax.add_patch(halo)
        dot = Circle((x, y), 0.17, fc=colour, ec="none", zorder=3.85)
        dot.set_gid(f"feature:{n}")
        ax.add_patch(dot)


def _draw_links(
    ax, layout: Layout, loops: list[Loop], sses: list[SSE], look: Style, lw: float, pt_per_unit: float
) -> None:
    links = layout.links
    if links is None:
        return
    chains = layout.res_chain
    shown = lambda r: not layout.partial or chains[r] not in layout.partial  # noqa: E731 - fragments: no links
    if look.disulfides:
        bars = []
        for i, j in links.disulfides:
            if not (shown(i) and shown(j)):
                continue
            a, b = disulfide_points(layout, loops, sses, chains, i, j)
            if a is not None and b is not None:
                bars.append((i, j, a, b))
        width = lw * 2.4 / pt_per_unit  # bar width in page units
        step = disulfide_step(lw, pt_per_unit)
        lanes = disulfide_lanes([(a, b) for _, _, a, b in bars], width)
        keep_out = [p.rect for p in layout.placed.values()] + [g.rect for g in layout.ghosts]
        keep_out += [box for _, box in label_boxes(layout, sses)]
        for (i, j, a, b), lane in zip(bars, lanes):
            side = disulfide_side(a, b, lane, keep_out, width, step) if lane else 1
            route = disulfide_route(a, b, lane, side, step)
            bar = PathPatch(
                Path(route, [Path.MOVETO] + [Path.LINETO] * (len(route) - 1)),
                fc="none",
                ec=_SULFUR,
                lw=lw * 2.4,
                capstyle="round",
                joinstyle="round",
                zorder=3.7,
            )
            bar.set_gid(f"disulfide:{i}-{j}")
            ax.add_patch(bar)
            dots = ax.scatter(
                [a[0], b[0]],
                [a[1], b[1]],
                s=(lw * 5.5) ** 2,
                c=_SULFUR,
                edgecolors="#8a6d00",
                linewidths=lw * 0.5,
                zorder=3.75,
            )
            dots.set_gid(f"ss-dot:{i}-{j}")
    if look.glycans:
        x0, y0, x1, y1 = layout.bounds
        centre = np.array([(x0 + x1) / 2, (y0 + y1) / 2])
        taken = [p.rect for p in layout.placed.values()] + [g.rect for g in layout.ghosts]
        taken += [box for _, box in label_boxes(layout, sses)]
        for loop in loops:  # loop lines too, as thin boxes along each segment
            pts = np.asarray(loop.points, float)
            for a, b in zip(pts, pts[1:]):
                taken.append(
                    (min(a[0], b[0]) - 0.06, min(a[1], b[1]) - 0.06, max(a[0], b[0]) + 0.06, max(a[1], b[1]) + 0.06)
                )
        size, step = 0.34, 0.42
        for r, sugars in links.glycans:
            if not shown(r):
                continue
            at = residue_point(layout, loops, sses, chains, r)
            if at is None:
                continue
            at = np.asarray(at, float)
            count = min(len(sugars), 5)
            away = at - centre
            best = None
            for out in (np.array([1.0, 0.0]), np.array([-1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.0, -1.0])):
                for lead in np.arange(0.35, 6.0, 0.1):  # the shortest stem that clears everything
                    boxes = [
                        centred_box(*(at + out * (lead + step * k)), size / 2 + 0.03, size / 2 + 0.03)
                        for k in range(count)
                    ]
                    if not any(overlaps_any(b, taken) for b in boxes):
                        cost = lead - 0.15 * float(np.dot(out, away) > 0)  # ties: point outward
                        if best is None or cost < best[0]:
                            best = (cost, out, lead, boxes)
                        break
            if best is None:  # nowhere clear: the direction and stem crossing the fewest things
                tries = []
                for out in (np.array([1.0, 0.0]), np.array([-1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.0, -1.0])):
                    for lead in np.arange(0.35, 6.0, 0.25):
                        boxes = [centred_box(*(at + out * (lead + step * k)), size / 2, size / 2) for k in range(count)]
                        tries.append((sum(overlaps_any(b, taken) for b in boxes), lead, out.tolist(), boxes))
                hits, lead, out, boxes = min(tries, key=lambda t: (t[0], t[1]))
                best = (0.0, np.array(out), lead, boxes)
            _, out, lead, boxes = best
            taken += boxes
            tip = at + out * (lead + step * (count - 0.5))
            stem = PathPatch(
                Path([tuple(at), tuple(tip)], [Path.MOVETO, Path.LINETO]),
                fc="none",
                ec=_legible("#2b2b2b"),
                lw=lw * 0.8,
                zorder=3.75,
            )
            stem.set_gid(f"glycan:{r}:stem")
            ax.add_patch(stem)
            for k, name in enumerate(sugars[:5]):
                shape, colour = _SNFG.get(name, ("circle", "#ffffff"))
                _symbol(ax, shape, colour, at + out * (lead + step * k), size, lw, f"glycan:{r}:{k}")


def _chevron(points, size: float) -> Path | None:
    """An open arrowhead ('>') at the middle of the loop's longest segment, pointing on along the loop."""
    pts = np.asarray(points, float)
    if len(pts) < 2:
        return None
    k = max(range(len(pts) - 1), key=lambda j: np.hypot(*(pts[j + 1] - pts[j])))
    a, b = pts[k], pts[k + 1]
    if np.hypot(*(b - a)) < 3 * size:
        return None
    u = (b - a) / np.hypot(*(b - a))
    n, m = np.array([-u[1], u[0]]), (a + b) / 2
    return Path(
        np.array([m - u * size + n * size * 0.8, m + u * size * 0.6, m - u * size - n * size * 0.8]),
        [Path.MOVETO, Path.LINETO, Path.LINETO],
    )


KYTE_DOOLITTLE = {
    "ILE": 4.5,
    "VAL": 4.2,
    "LEU": 3.8,
    "PHE": 2.8,
    "CYS": 2.5,
    "MET": 1.9,
    "ALA": 1.8,
    "GLY": -0.4,
    "THR": -0.7,
    "SER": -0.8,
    "TRP": -0.9,
    "TYR": -1.3,
    "PRO": -1.6,
    "HIS": -3.2,
    "GLU": -3.5,
    "GLN": -3.5,
    "ASP": -3.5,
    "ASN": -3.5,
    "LYS": -3.9,
    "ARG": -4.5,
    "MSE": 1.9,
}
_PROPERTY = {
    "bfactor": ("RdYlBu_r", 0.06, 0.94, "B-factor: rigid → flexible"),
    "hydropathy": ("BrBG_r", 0.08, 0.92, "hydropathy: hydrophilic → hydrophobic"),
}
_SSTYPE = {"H": "#c8414b", "E": "#1f7a8c", "G": "#e8969c"}  # Richardson: helices warm, strands cool


def property_colour(t: float, kind: str) -> str:
    """Colour at fraction t (0 = low, 1 = high) of a property ramp (bfactor, hydropathy)."""
    name, lo, hi, _ = _PROPERTY[kind]
    return to_hex(colormaps[name](lo + float(np.clip(t, 0.0, 1.0)) * (hi - lo)))


PLDDT_BANDS = (  # AlphaFold's own confidence colours: (lower bound, colour, legend text)
    (90.0, "#0053d6", "very high (> 90)"),
    (70.0, "#65cbf3", "confident (70–90)"),
    (50.0, "#ffdb13", "low (50–70)"),
    (-np.inf, "#ff7d45", "very low (< 50)"),
)


def _flatten(path: Path, per_curve: int = 12) -> np.ndarray:
    """A path's points along its length, Bézier pieces sampled."""
    pts = []
    for seg, _ in path.iter_bezier():
        t = np.linspace(0.0, 1.0, per_curve if seg.degree > 1 else 2)
        q = seg(t)
        pts.extend(q if not pts else q[1:])
    return np.asarray(pts, float)


def _cut(pts: np.ndarray, at: np.ndarray, t0: float, t1: float) -> np.ndarray:
    """The stretch of polyline pts between arc lengths t0 and t1 (at: cumulative arc length)."""
    inner = pts[(at > t0) & (at < t1)]
    a = np.array([np.interp(t0, at, pts[:, 0]), np.interp(t0, at, pts[:, 1])])
    b = np.array([np.interp(t1, at, pts[:, 0]), np.interp(t1, at, pts[:, 1])])
    return np.vstack([a, inner, b])


def residue_runs(path: Path, colours: list) -> list[tuple[str, np.ndarray]]:
    """Split a loop path into runs of equal colour, one equal share of its length per residue, in order."""
    pts = _flatten(path)
    at = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    total, n, runs, start = at[-1], len(colours), [], 0
    for k in range(1, n + 1):
        if k == n or colours[k] != colours[start]:
            runs.append((colours[start], _cut(pts, at, total * start / n, total * k / n)))
            start = k
    return runs


def _typical(colours: list) -> str | None:
    """The most common colour; ties go to the one met first, so the result never depends on hashing."""
    found = [c for c in colours if c is not None]
    return Counter(found).most_common(1)[0][0] if found else None


def plddt_colour(v: float) -> str:
    return next(c for lo, c, _ in PLDDT_BANDS if v >= lo)


def plddt_percent(values) -> list[float]:
    """pLDDT on AlphaFold's 0-100 scale. Some predictors (ESMFold, some ColabFold exports) store 0-1 in the
    B-factor column: when no value exceeds 1, the whole model is read as fractions."""
    v = [float(x) for x in values]
    return [x * 100.0 for x in v] if v and max(v) <= 1.0 else v


CONSURF = (  # ConSurf's nine conservation grades, 1 (variable, turquoise) .. 9 (conserved, maroon)
    "#10c8d1",
    "#8cffff",
    "#d7ffff",
    "#eaffff",
    "#ffffff",
    "#fcedf4",
    "#fac9de",
    "#f07dab",
    "#a02560",
)
NO_DATA = "#d9d9d9"  # residues the alignment says nothing about


def grade_colour(grade: int) -> str:
    """ConSurf's colour for a grade 1..9; grade 0 (no score) is grey."""
    return CONSURF[int(np.clip(grade, 1, 9)) - 1] if grade else NO_DATA


def _residue_values(layout: Layout, kind: str) -> list[float]:
    if kind == "bfactor":
        return list(layout.res_b) or [0.0] * len(layout.res_chain)
    return [KYTE_DOOLITTLE.get(n, 0.0) for n in layout.res_name]


def _colouring(sses: list[SSE], chains: list[str], colors: dict[str, str], look: Style, layout: Layout = None):
    """(element id -> colour, residue colour function). By chain: the chain's colour. By sequence: a ramp
    from the chain's first element (N) to its last (C), restarted for every chain."""
    if look.color_by == "sstype":
        return {s.id: _SSTYPE.get(s.kind, "#7f7f7f") for s in sses}, (lambda r, chain=None: "#7f7f7f")
    if look.color_by == "conservation":  # grades relative to this protein (sequence.grades), ConSurf's colours
        g = np.array(layout.res_grade if layout is not None and layout.res_grade else [], int)
        if not g.size or not g.any():
            return {s.id: NO_DATA for s in sses}, (lambda r, chain=None: None)

        def element(s: SSE) -> str:
            mine = g[s.start : s.end + 1]
            mine = mine[mine > 0]
            return grade_colour(int(round(float(mine.mean())))) if mine.size else NO_DATA

        return {s.id: element(s) for s in sses}, (lambda r, chain=None: grade_colour(int(g[r])))
    if look.color_by == "plddt":  # absolute bands, never rescaled: 'confident' means the same in every figure
        b = plddt_percent(layout.res_b) if layout is not None else []
        if not b:
            return {s.id: plddt_colour(0.0) for s in sses}, (lambda r, chain=None: None)
        return (
            {s.id: plddt_colour(float(np.mean(b[s.start : s.end + 1]))) for s in sses},
            lambda r, chain=None: plddt_colour(b[r]),
        )
    if look.color_by in _PROPERTY:
        values = _residue_values(layout, look.color_by) if layout is not None and layout.res_chain else []
        if not values:
            return {s.id: property_colour(0.5, look.color_by) for s in sses}, (lambda r, chain=None: None)
        means = {s.id: float(np.mean(values[s.start : s.end + 1])) for s in sses}
        lo, hi = min(means.values()), max(means.values())

        def scaled(v: float) -> float:
            return (v - lo) / (hi - lo) if hi > lo else 0.5

        return (
            {k: property_colour(scaled(v), look.color_by) for k, v in means.items()},
            lambda r, chain=None: property_colour(scaled(values[r]), look.color_by),
        )
    if look.color_by == "shade":
        out, span = {}, {}
        groups = group_by_chain(sses)
        for chain in chains:
            mine = groups[chain]
            for k, s in enumerate(mine):
                out[s.id] = _shade(colors[chain], k / (len(mine) - 1) if len(mine) > 1 else 0.5)
            span[chain] = (mine[0].start, mine[-1].end)

        def by_shade(r: int, chain: str | None = None) -> str | None:
            if chain not in span:
                return None
            a, b = span[chain]
            return _shade(colors[chain], float(np.clip((r - a) / max(b - a, 1), 0.0, 1.0)))

        return out, by_shade
    if look.color_by == "chain":

        def by_chain(r: int, chain: str | None = None) -> str | None:
            return colors.get(chain)

        return {s.id: colors[s.chain] for s in sses}, by_chain
    out, span = {}, {}
    groups = group_by_chain(sses)
    for chain in chains:
        mine = groups[chain]
        for k, s in enumerate(mine):
            out[s.id] = sequence_colour(k / (len(mine) - 1) if len(mine) > 1 else 0.0, look.sequence_map)
        span[chain] = (mine[0].start, mine[-1].end)

    def by_residue(r: int, chain: str | None = None) -> str | None:
        if chain not in span:
            return None
        a, b = span[chain]
        return sequence_colour((r - a) / max(b - a, 1), look.sequence_map)

    return out, by_residue


def element_colours(layout: Layout, sses: list[SSE], look: Style):
    """(element id -> fill colour, residue colour function) exactly as the figure draws them: the colouring
    mode, then marked elements, then highlight greys."""
    chains = _chains_in_order(sses)
    colors = chain_colors(chains, look.palette)
    element_colour, residue_colour = _colouring(sses, chains, colors, look, layout)
    for entry in filter(None, (e.strip() for e in look.mark.split(","))):  # marked elements: their own colour
        ref, _, colour = entry.partition("=")
        try:
            colour = to_hex(colour.strip()) if colour else MARK_DEFAULT
        except ValueError:
            raise ValueError(f"mark colour {colour!r} is not a colour (use a name or #rrggbb)") from None
        for sid in resolve(layout, sses, ref):
            element_colour[sid] = colour
    if layout.focus is not None:  # highlight: everything outside the focus in one quiet grey
        dim = {s.chain for s in sses} - layout.focus
        element_colour = {k: (MATE_GREY if s.chain in dim else element_colour[k]) for k, s in ((s.id, s) for s in sses)}
        base_residue = residue_colour

        def residue_colour(r, chain=None, _base=base_residue, _dim=dim):  # noqa: F811
            return MATE_GREY if chain in _dim else _base(r, chain)

    return element_colour, residue_colour


def draw(layout: Layout, loops: list[Loop], sses: list[SSE], title: str | None = None, **kw) -> Figure:
    """`look` carries every visual choice; palette/style (fill)/loop_style are shorthands that override it.
    `model` (the io.read_model result the layout was built from) is needed only for look.hbonds."""
    look = (kw.get("look") or Style()).validate()
    with page(look):  # fonts/colours are fixed when text is made: set them here
        return pin_fonts(_draw(layout, loops, sses, title, **kw))


class HBondLine(NamedTuple):
    """One drawn hydrogen-bond line: a residue (or nucleotide) pair carrying at least one bond."""

    gid: str  # hbond:A:42-1:DG16
    ends: tuple  # ((x, y), (x, y)): page points of the first and second residue
    length: float  # page units
    bonds: tuple  # hbonds.HBond, the donor-acceptor atom pairs between the two, shortest first
    residues: tuple  # Backbone indices of the protein ends
    labels: tuple = ()  # model.ResLabel of the two ends, in the order of `ends`


class HBondOverlay(NamedTuple):
    lines: list  # HBondLine drawn, shortest first
    total: int  # residue pairs with a bond and both ends drawn, before the hbond_max cap
    unplaced: int = 0  # residue pairs with a bond that are left out because an end is not drawn

    def legend(self) -> str:
        n = len(self.lines)
        text = str(n) if n == self.total else f"{n} of {self.total}"
        return f"H-bonds ({text}{f', {self.unplaced} not shown' if self.unplaced else ''})"


def nucleotide_points(layout: Layout, style: str = "ladder") -> dict[tuple[int, int], np.ndarray]:
    """Where each drawn nucleotide (strand, index) sits in the stem `style` (as `_na_style` resolved it): its rail
    point on a stem (its helix backbone point in the ribbon style), its run point on a path. Beads, letters and
    hydrogen-bond ends all go here."""
    where: dict[tuple[int, int], np.ndarray] = {}
    for d in layout.stems:
        for key in d.axial:
            where[key] = _ribbon_point(d, key) if style == "ribbon" else d.point(key)
    for path in layout.na_paths:
        for i, xy in zip(path.nts, path.nt_xy):
            where[(path.strand, i)] = np.asarray(xy)
    return where


def hbond_overlay(layout: Layout, loops: list[Loop], sses: list[SSE], look: Style, model=None) -> HBondOverlay:
    """The hydrogen-bond lines `look.hbonds` asks for, from the model the layout was built from: one per residue
    pair (protein before nucleotide, else by tag), between their points on the page. Bonds with an end that is not
    drawn (cropped nucleic acid, chains outside the focus, residues past a chain's elements) are left out. At most
    `look.hbond_max` lines, the shortest on the page kept (not the shortest in Å). `model` is the io.read_model
    result the layout was built from."""
    sel = parse_selection(look.hbonds)
    if sel.kind == "off":
        return HBondOverlay([], 0)
    if model is None:
        raise ValueError(f"hbonds={look.hbonds} needs the structure: pass model= (the io.read_model result)")
    protein = {(l.chain, l.seq, l.icode): (k, l.name) for k, l in enumerate(layout.res_labels)}
    nts = {}
    if layout.nucleic is not None:
        nts = {
            (lab.chain, lab.seq, lab.icode): ((s, i), lab.name)
            for s, strand in enumerate(layout.nucleic.strands)
            for i, lab in enumerate(strand.labels)
        }
    where = nucleotide_points(layout, _na_style(layout, look))
    cache: dict = {}

    def locate(lab):
        """(sort key, page point, Backbone index or None) of a residue label, or None when it is not drawn."""
        ident = (lab.chain, lab.seq, lab.icode)
        if ident not in cache:
            cache[ident] = None
            hit = protein.get(ident)
            shown = lab.chain not in layout.partial and not (
                layout.focus_chains and lab.chain not in layout.focus_chains
            )
            if hit is not None and hit[1] == lab.name:
                at = residue_point(layout, loops, sses, layout.res_chain, hit[0]) if shown else None
                if at is not None:
                    cache[ident] = ((0, *ident), at, hit[0], lab)
            elif ident in nts and nts[ident][1] == lab.name and nts[ident][0] in where:
                cache[ident] = ((1, *ident), tuple(float(v) for v in where[nts[ident][0]]), None, lab)
        return cache[ident]

    pairs: dict[str, list] = {}
    unplaced = set()
    for bond in find_hbonds(model, sel):
        a, b = locate(bond.donor.res), locate(bond.acceptor.res)
        if a is None or b is None:
            unplaced.add(frozenset(((r.chain, r.seq, r.icode) for r in (bond.donor.res, bond.acceptor.res))))
            continue
        a, b = sorted((a, b), key=lambda e: e[0])
        gid = f"hbond:{a[0][1]}:{a[0][2]}{a[0][3]}-{b[0][1]}:{b[0][2]}{b[0][3]}"
        pairs.setdefault(gid, [a, b, []])[2].append(bond)
    lines = []
    for gid, (a, b, bonds) in pairs.items():
        length = float(np.hypot(a[1][0] - b[1][0], a[1][1] - b[1][1]))
        if length > 0:
            residues = tuple(e[2] for e in (a, b) if e[2] is not None)
            lines.append(
                HBondLine(
                    gid, (a[1], b[1]), length, tuple(sorted(bonds, key=lambda x: x.distance)), residues, (a[3], b[3])
                )
            )
    lines.sort(key=lambda l: (l.length, l.gid))
    return HBondOverlay(lines[: look.hbond_max], len(lines), len(unplaced))


def hbond_alpha(length: float, diagonal: float) -> float:
    """Opacity of a hydrogen-bond line: 1 for short lines, falling linearly to 0.4 at half the page diagonal."""
    t = (length / diagonal - 0.15) / 0.35 if diagonal > 0 else 0.0
    return 1.0 - 0.6 * min(max(t, 0.0), 1.0)


def _draw_hbonds(c: _Canvas, overlay: HBondOverlay) -> None:
    """Thin dashed lines between the residues of each hydrogen-bonded pair: over the panels and the nucleic-acid
    rungs, under the backbone rails, the contact beads, the elements and the letters (gid hbond:<a>-<b>). Long
    lines fade (alpha 1 up to 15% of the content's diagonal, 0.4 from 50% on), so a crowd of them stays legible."""
    x0, y0, x1, y1 = c.content
    diag = float(np.hypot(x1 - x0, y1 - y0))
    for line in overlay.lines:
        (artist,) = c.ax.plot(
            *np.c_[line.ends[0], line.ends[1]],
            color=c.look.hbond_color,
            lw=c.loop_lw * 0.85,
            ls=(0, (3, 2)),
            alpha=hbond_alpha(line.length, diag),
            solid_capstyle="butt",
            zorder=2.15,
        )
        artist.set_gid(line.gid)


@dataclass
class _Canvas:
    """What every drawing step shares: the axes, the resolved look, colours and the page geometry."""

    fig: Figure
    ax: Axes
    look: Style
    colors: dict[str, str]
    element_colour: dict[str, str]
    residue_colour: Callable[[int, str], str | None]
    font: float
    lw: float
    loop_lw: float
    pt_per_unit: float
    pad: float
    content: tuple[float, float, float, float]  # x0, y0, x1, y1: elements, marks, at least one legend entry wide
    page: tuple[float, float, float, float]  # xmin, xmax, ymin, ymax
    backbone: tuple[dict[str, str], dict[str, str]]  # backbone_colours: nucleic-acid chains, protein loop overrides


class _Legend(NamedTuple):
    entries: list
    entry: float  # width of one legend entry (page units)
    per_row: int
    height: float


def _resolve_look(look, palette, style, loop_style) -> Style:
    style = _ALIASES["fill"].get(style, style)
    if style is not None and style not in STYLES:
        raise ValueError(f"unknown style {style!r}; choose from {', '.join(STYLES)}")
    if loop_style is not None and loop_style not in LOOP_STYLES:
        raise ValueError(f"unknown loop style {loop_style!r}; choose from {', '.join(LOOP_STYLES)}")
    return replace(
        look or Style(),
        **{k: v for k, v in (("palette", palette), ("fill", style), ("loops", loop_style)) if v is not None},
    ).validate()


def _empty_figure() -> Figure:
    fig = Figure(figsize=(4.0, 1.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(0.5, 0.5, "no helices or strands found", ha="center", va="center", fontsize=9, transform=ax.transAxes)
    return fig


def _legend_entries(
    layout: Layout,
    look: Style,
    chains: list[str],
    colors: dict[str, str],
    umarks,
    na_colours: dict[str, str],
    hb: HBondOverlay | None = None,
    linked: list[dict] = (),
) -> list:
    entries = [(colors[c], f"Chain {c}") for c in chains]
    if len(chains) == 1 and layout.focus is None and look.color_by == "chain":
        entries = []  # one chain in its own colour: a "Chain A" box says nothing
    if layout.focus is not None:
        entries = [(colors[c], f"Chain {c}") for c in chains if c in layout.focus]
        if layout.focus_chains:
            if layout.partial:
                entries.append((MATE_GREY, "neighbouring subunits (′, ″)"))
        else:
            entries.append((MATE_GREY, "symmetry copies" if look.highlight in ("asu", "protomer") else "other chains"))
    if look.color_by == "sequence":
        entries = [(None, "N → C" + (" (each chain)" if len(chains) > 1 else ""))]
    elif look.color_by in _PROPERTY:
        entries = [(("ramp", look.color_by), _PROPERTY[look.color_by][3])]
    elif look.color_by == "plddt":
        entries = [(c, text) for _, c, text in PLDDT_BANDS]
    elif look.color_by == "conservation":
        entries = [(("grades",), (layout.cons_label or "conservation") + ": variable → conserved")]
    elif look.color_by == "sstype":
        kinds = {p.sse.kind for p in layout.placed.values()}
        entries = [
            (_SSTYPE[k], name) for k, name in (("H", "helix"), ("E", "strand"), ("G", "3₁₀ helix")) if k in kinds
        ]
    for kind in dict.fromkeys(m["feature"].type for m in umarks):  # UniProt site types shown
        entries.append((SITE_COLOURS.get(kind, "#555555"), f"{kind.lower()} (UniProt)"))
    if layout.nucleic is not None and layout.nucleic.strands and look.legend:
        strands = layout.nucleic.strands
        ids = sorted({s.chain for s in strands})
        groups: dict[str, list] = {}  # one entry per backbone colour: the chains drawn in it
        for chain in ids:
            groups.setdefault(na_colours[chain], []).append(chain)
        for colour, mine in groups.items():
            rna = [any(s.rna for s in strands if s.chain == c) for c in mine]
            kind = "RNA" if all(rna) else "DNA" if not any(rna) else "DNA/RNA"
            entries.append((colour, f"{kind} ({', '.join(mine)})"))
    if layout.nucleic is not None and look.legend:  # ties drawn as two labelled marks: what the letters mean
        tie = next((t for t in layout.na_ties if t.label), None)
        if tie is not None:
            entries.append((na_colours[layout.nucleic.strands[tie.a[0]].chain], TIE_LEGEND))
    if hb is not None and hb.lines and look.legend:
        entries.append((look.hbond_color, hb.legend()))
    if linked and look.legend:
        entries.append((LIGAND_FILL, "chromophore / linked residue"))
    return entries


def _canvas(
    layout: Layout,
    look: Style,
    title,
    chains,
    colors,
    element_colour,
    residue_colour,
    marks,
    umarks,
    hb,
    linked=(),
    loops: list[Loop] = (),
) -> tuple[_Canvas, _Legend]:
    """The figure, its axes and every size the steps share; also lays the legend out, as it sets the page width."""
    x0, y0, x1, y1 = layout.bounds
    for m in [*marks, *({"centre": m["text_at"], "half": m["half"]} for m in linked if m["text_at"])]:
        (cx, cy), (hw, hh) = m["centre"], m["half"]  # markers may sit just outside the elements' bounds
        x0, y0, x1, y1 = min(x0, cx - hw), min(y0, cy - hh), max(x1, cx + hw), max(y1, cy + hh)
    pad = 1.0
    # a loop routed around the outermost elements may leave the bounds: the page (bounds + pad, the legend below)
    # holds every loop point with _LOOP_ROOM to spare, so no loop is clipped at the edge or runs into the legend
    pts = np.concatenate([np.asarray(lp.points, float).reshape(-1, 2) for lp in loops] or [np.empty((0, 2))])
    if len(pts):
        (lx0, ly0), (lx1, ly1) = pts.min(axis=0), pts.max(axis=0)
        room = _LOOP_ROOM - pad
        x0, y0, x1, y1 = min(x0, lx0 - room), min(y0, ly0 - room), max(x1, lx1 + room), max(y1, ly1 + room)
    top = 1.9 if title else 0.0
    backbone = backbone_colours(layout, chains, look, colors)  # once per figure
    entries = _legend_entries(layout, look, chains, colors, umarks, backbone[0], hb, linked)

    def legend_entry(x1: float) -> float:  # one legend entry's width in page units: swatch, text, gap
        scale = min(_IN_PER_UNIT, _MAX_WIDTH_IN / (x1 - x0 + 2 * pad))  # as the figure below
        font_u = _font_pt(scale) * look.font_scale / (scale * 72)  # text drawn at that size, in page units
        # the unhinted advance width, which is what is drawn (ink extents differ between matplotlib versions);
        # measured at a size matplotlib does not clamp (it floors font sizes at 1 pt) and scaled to font_u
        widest = max(
            (text_to_path.get_text_width_height_descent(text, _MEASURE, ismath=False)[0] for _, text in entries),
            default=0.0,
        )
        return 0.85 + widest * font_u / _MEASURE.get_size_in_points() + 0.9

    entry = legend_entry(x1)
    if x0 + entry - 0.9 > x1:  # never narrower than one legend entry; a wider page may shrink, so measure again
        entry = legend_entry(x0 + entry - 0.9)
        x1 = max(x1, x0 + entry - 0.9)
    per_row = max(1, int((x1 - x0 + 0.9) // entry))
    if not look.legend:
        entries = []
    legend_h = 0.8 + _LEGEND_ROW * -(-len(entries) // per_row) if entries else 0.3
    if layout.context and look.legend:
        legend_h += 1.3
    xmin, xmax, ymin, ymax = x0 - pad, x1 + pad, y0 - pad - legend_h, y1 + pad + top
    w_u, h_u = xmax - xmin, ymax - ymin
    scale = min(_IN_PER_UNIT, _MAX_WIDTH_IN / w_u)
    fig = Figure(figsize=(w_u * scale, h_u * scale), facecolor=look.background)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    pt_per_unit = scale * 72
    font = _font_pt(scale) * look.font_scale
    lw = float(np.clip(0.9 * scale / _IN_PER_UNIT, 0.5, 1.2))
    canvas = _Canvas(
        fig,
        ax,
        look,
        colors,
        element_colour,
        residue_colour,
        font,
        lw,
        lw * look.loop_width,
        pt_per_unit,
        pad,
        (x0, y0, x1, y1),
        (xmin, xmax, ymin, ymax),
        backbone,
    )
    return canvas, _Legend(entries, entry, per_row, legend_h)


def _draw_loops(c: _Canvas, layout: Layout, loops: list[Loop], sses: list[SSE]) -> None:
    """Each loop: a casing in the background colour, its coloured line, arrows, residue runs and length labels."""
    boxes = np.array([p.rect for p in [*layout.placed.values(), *layout.ghosts]]).reshape(-1, 4)  # none: DNA only
    boxes = boxes + [1e-6, 1e-6, -1e-6, -1e-6]
    sse_by_id = {s.id: s for s in sses}
    over = c.backbone[1]
    for i, loop in enumerate(loops):
        if loop.fallback:  # already a smooth curve
            path = _polyline(loop.points)
        elif c.look.loops == "curved":
            path = _smooth_path(loop.points, boxes)
        else:
            path = _rounded_path(loop.points)
        name = f"{loop.a_id}>{loop.b_id}"
        casing = PathPatch(
            path,
            fc="none",
            ec=c.look.background,
            lw=c.loop_lw * 3.0,
            capstyle="round",
            joinstyle="round",
            zorder=1 + i * 1e-3,
        )
        casing.set_gid(f"loopcase:{name}")
        c.ax.add_patch(casing)
        a_sse = sse_by_id.get(loop.a_id)
        loop_c = {
            "black": c.look.ink,
            "chain": c.colors.get(a_sse.chain, c.look.ink) if a_sse else c.look.ink,
            "element": c.element_colour.get(loop.a_id, c.look.ink),
            "residue": c.look.ink,  # replaced by per-residue runs below
        }[c.look.loop_color]
        loop_c = _legible(loop_c)
        mine = over.get(a_sse.chain) if a_sse else None  # backbone_colors for this protein chain
        if mine:
            loop_c = mine
        if layout.focus is not None and a_sse is not None and a_sse.chain not in layout.focus:
            loop_c = _MATE_LINE
        if c.look.loop_arrows:
            arrow = _chevron(loop.points, 0.22 * c.look.loop_width**0.5)
            if arrow is not None:
                head = PathPatch(
                    arrow,
                    fc="none",
                    ec=loop_c,
                    lw=c.loop_lw,
                    capstyle="round",
                    joinstyle="round",
                    zorder=1 + i * 1e-3 + 7e-4,
                )
                head.set_gid(f"loop-arrow:{name}")
                c.ax.add_patch(head)
        line = PathPatch(
            path,
            fc="none",
            ec=loop_c,
            lw=c.loop_lw,
            capstyle="round",
            joinstyle="round",
            ls=(0, (3, 2)) if loop.dashed else "-",
            zorder=1 + i * 1e-3 + 5e-4,
        )
        line.set_gid(f"loop:{name}")
        c.ax.add_patch(line)
        b_sse = sse_by_id.get(loop.b_id)
        if a_sse is not None and b_sse is not None and b_sse.start - a_sse.end - 1 >= LONG_LOOP:
            pts = np.asarray(loop.points, float)  # the residue count, beside the loop's longest straight run
            k = max(range(len(pts) - 1), key=lambda j: np.hypot(*(pts[j + 1] - pts[j])))
            mid, d = (pts[k] + pts[k + 1]) / 2, pts[k + 1] - pts[k]
            vertical = abs(d[1]) > abs(d[0])
            t = c.ax.text(
                mid[0] + (0.15 if vertical else 0.0),
                mid[1] + (0.0 if vertical else 0.22),
                f"{b_sse.start - a_sse.end - 1} aa",
                ha="left" if vertical else "center",
                va="center" if vertical else "bottom",
                fontsize=c.font * 0.72,
                color=_legible("#5d6877"),
                zorder=4,
                rotation=0,
            )
            t.set_gid(f"loop-length:{name}")
        if c.look.loop_color == "residue" and loop_c != _MATE_LINE and a_sse is not None and not mine:
            b_sse = sse_by_id.get(loop.b_id)
            span = range(a_sse.end + 1, b_sse.start) if b_sse is not None else range(0)
            ends = [a_sse.end, b_sse.start] if b_sse is not None else [a_sse.end]
            shades = [c.residue_colour(r, a_sse.chain) for r in (span or ends)]
            shades = [col or c.look.ink for col in shades]  # true colours, matching the legend ...
            line.set_edgecolor(_RUN_EDGE)  # ... on a dark edge (the loop's own line, also the hover target)
            line.set_linewidth(c.loop_lw * _RUN_EDGE_W)
            for k, (shade, piece) in enumerate(residue_runs(path, shades)):
                run = PathPatch(
                    Path(piece),
                    fc="none",
                    ec=shade,
                    lw=c.loop_lw,
                    capstyle="butt",
                    joinstyle="round",
                    ls=(0, (3, 2)) if loop.dashed else "-",
                    zorder=1 + i * 1e-3 + 6e-4,
                )
                run.set_gid(f"loop-seg:{name}:{k}")
                c.ax.add_patch(run)


def _draw_membrane(c: _Canvas, layout: Layout) -> None:
    """The lipid bilayer behind everything."""
    x0, y0, x1, y1 = c.content
    if layout.membrane is None:
        return
    # the lipid bilayer behind everything: a band with head groups on both faces; one per row of a wrapped stack
    for k, (b0, b1) in enumerate(layout.membrane.get("bands", [layout.membrane["y"]])):
        row = f":{k}" if k else ""
        left, right = x0 - c.pad * 0.7, x1 + c.pad * 0.7
        band = Rectangle((left, b0), right - left, b1 - b0, fc=_LIPID, ec="none", zorder=0.1)
        band.set_gid(f"membrane{row}")
        c.ax.add_patch(band)
        xs = np.arange(left + 0.2, right, 0.42)
        for face, y in (("low", b0), ("high", b1)):
            heads = c.ax.scatter(
                xs,
                np.full(len(xs), y),
                s=(c.lw * 4.2) ** 2,
                c=_HEADS,
                edgecolors=_HEAD_EDGE,
                linewidths=c.lw * 0.4,
                zorder=0.15,
            )
            heads.set_gid(f"membrane-heads:{face}{row}")
        below = layout.membrane["inside"] == "below"
        sides = (("out" if below else "in", b1 + 0.35, "bottom"), ("in" if below else "out", b0 - 0.35, "top"))
        for text, y, va in sides if layout.membrane.get("sides", True) else ():
            t = c.ax.text(left + 0.1, y, text, ha="left", va=va, fontsize=c.font * 0.9, style="italic", color="#7a6a45")
            t.set_gid(f"membrane-side:{text}{row}")


def _draw_domains(c: _Canvas, layout: Layout) -> None:
    """Named domains: a soft tinted panel each."""
    for k, (name, ids) in enumerate(layout.domains):  # named domains: a soft tinted panel, name at top left
        if not ids:
            continue
        x0, y0, x1, y1 = domain_panel(layout, name)
        tone = _DOMAIN_TONES[k % len(_DOMAIN_TONES)]
        panel = FancyBboxPatch(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            boxstyle="round,pad=0,rounding_size=0.6",
            fc=_mix(tone, c.look.background, 0.86),
            ec=_mix(tone, c.look.background, 0.55),
            lw=c.lw * 0.6,
            zorder=0.2,
        )
        panel.set_gid(f"domain:{name}")
        c.ax.add_patch(panel)
        t = c.ax.text(
            x0 + 0.35,
            y1 - 0.3,
            name,
            ha="left",
            va="top",
            fontsize=c.font * 1.05,
            fontweight="bold",
            color=darken(tone, 0.75),
            zorder=0.25,
        )
        t.set_gid(f"domain-label:{name}")


def _draw_sheet_panels(c: _Canvas, layout: Layout) -> None:
    """A pale panel behind each sheet."""
    if c.look.sheet_panels:
        ghosts = {g.sse.id: g for g in layout.ghosts}
        for k, ids in enumerate(layout.sheets):
            if len(ids) < 2:
                continue
            rects = np.array([layout.placed[i].rect for i in ids] + [ghosts[i].rect for i in ids if i in ghosts])
            x0, y0 = rects[:, 0].min() - 0.3, rects[:, 1].min() - 0.3
            x1, y1 = rects[:, 2].max() + 0.3, rects[:, 3].max() + 0.3
            panel = FancyBboxPatch(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                boxstyle="round,pad=0,rounding_size=0.35",
                fc=_mix(c.look.background, c.look.ink, 0.07) if c.look.background != "#ffffff" else _PANEL,
                ec="none",
                zorder=0.3,
            )
            panel.set_gid(f"sheet-panel:{k}")
            c.ax.add_patch(panel)


def _draw_elements(c: _Canvas, layout: Layout) -> None:
    """The strands, 3-10 helices and helices, with their labels."""
    for p in layout.placed.values():
        color = c.element_colour[p.sse.id]
        outline = c.look.fill == "outline"
        edge, edge_w = (color, c.lw * 1.2) if outline else (darken(color), c.lw * 0.8)
        if p.sse.kind == "E":
            strand_fill = c.look.background if outline else tint(color, 0.45) if c.look.fill == "pale" else color
            if c.look.fill == "pale":
                edge = darken(color, 0.5)  # pastel body, firm outline (illustration style)
            patch = Polygon(_arrow(p), closed=True, fc=strand_fill, ec=edge, lw=edge_w, joinstyle="round", zorder=3)
            patch.set_gid(f"strand:{p.sse.id}")
            c.ax.add_patch(patch)
            if c.look.labels:
                t = c.ax.text(
                    p.cx,
                    p.cy,
                    p.label,
                    ha="center",
                    va="center",
                    fontsize=c.font,
                    fontweight="bold",
                    color=_label_ink(color) if outline else text_color_on(strand_fill),
                    zorder=4,
                )
                t.set_gid(f"label:{p.sse.id}")
        elif p.sse.kind == "G":  # a short 3-10 helix: a small rounded box
            fill = {"bold": color, "pale": tint(color, 0.35), "outline": c.look.background}[c.look.fill]
            box = FancyBboxPatch(
                (-p.length / 2, -p.width / 2),
                p.length,
                p.width,
                boxstyle=f"round,pad=0,rounding_size={p.width * 0.3}",
                fc=fill,
                ec=edge,
                lw=edge_w,
                zorder=3,
            )
            box.set_transform(Affine2D().rotate(p.angle).translate(p.cx, p.cy) + c.ax.transData)
            box.set_gid(f"eta:{p.sse.id}")
            c.ax.add_patch(box)
            if c.look.labels:
                t = c.ax.text(
                    *helix_label_pos(p),
                    p.label,
                    ha="center",
                    va="center",
                    fontsize=c.font * 0.9,
                    color=_label_ink(color),
                    zorder=4,
                )
                t.set_gid(f"label:{p.sse.id}")
        else:
            fill, shade = {
                "bold": (color, darken(color, 0.72)),
                "pale": (tint(color, 0.15), tint(color, 0.55)),
                "outline": (c.look.background, _mix(color, c.look.background, 0.7)),
            }[c.look.fill]
            if c.look.helix_shading == "none":
                shade = fill
            front, back, shine = _coil(p, len(p.sse))
            behind = PathPatch(back, fc=shade, ec=edge, lw=edge_w * 0.8, joinstyle="round", zorder=2.9)
            behind.set_gid(f"helix-back:{p.sse.id}")
            c.ax.add_patch(behind)
            bar = PathPatch(front, fc=fill, ec=edge, lw=edge_w, joinstyle="round", zorder=3)
            bar.set_gid(f"helix:{p.sse.id}")
            c.ax.add_patch(bar)
            if shine is not None and c.look.helix_shading == "gloss" and c.look.fill != "outline":
                gloss = PathPatch(shine, fc=tint(fill, 0.4), ec="none", alpha=0.85, zorder=3.05)
                gloss.set_gid(f"helix-shine:{p.sse.id}")
                c.ax.add_patch(gloss)
            if c.look.labels:
                t = c.ax.text(
                    *helix_label_pos(p),
                    p.label,
                    ha="center",
                    va="center",
                    fontsize=c.font,
                    color=_label_ink(color),
                    zorder=4,
                )
                t.set_gid(f"label:{p.sse.id}")


def _draw_residue_numbers(c: _Canvas, layout: Layout) -> None:
    """First and last residue numbers at the ends of each element."""
    if c.look.residue_numbers:  # helices: beside each end; strands: inside the arrow near each end (or just past it)
        for p in layout.placed.values():
            nrm = np.array([-np.sin(p.angle), np.cos(p.angle)])
            colour = c.element_colour[p.sse.id]
            for end, port, num in (("N", p.n_port, p.sse.first.tag), ("C", p.c_port, p.sse.last.tag)):
                inward = p.direction * (1.0 if end == "N" else -1.0)
                if p.sse.kind in ("H", "G"):
                    away = -1.0 if nrm[1] >= 0 else 1.0  # the side opposite the helix name
                    along = 0.35 if p.length >= _SHORT_HELIX else -0.3  # a short box: just past each end instead
                    xy, ink = np.asarray(port) + inward * along + nrm * away * (p.width / 2 + 0.38), _legible("#4a4a4a")
                elif p.length >= 3.0:  # room for both numbers and the strand letter between them
                    reach = 0.35 if end == "N" else min(_HEAD_LEN, p.length * 0.45) + 0.3  # clear of the head
                    xy, ink = (
                        np.asarray(port) + inward * reach,
                        darken(colour) if c.look.fill == "outline" else text_color_on(colour),
                    )
                else:  # too short to hold numbers: just past the end, beside the loop stub
                    xy, ink = np.asarray(port) - inward * 0.45 + nrm * 0.42, _legible("#4a4a4a")
                t = c.ax.text(*xy, str(num), ha="center", va="center", fontsize=c.font * 0.62, color=ink, zorder=4.5)
                t.set_gid(f"resnum:{p.sse.id}:{end}")


def _draw_ghosts(c: _Canvas, layout: Layout) -> None:
    """Dashed copies of elements that belong to a neighbouring subunit."""
    for g in layout.ghosts:
        color = c.element_colour[g.sse.id]
        patch = Polygon(
            _arrow(g),
            closed=True,
            fc=tint(color, 0.8),
            ec=darken(color),
            lw=c.lw * 0.8,
            ls=(0, (2, 1.5)),
            joinstyle="round",
            zorder=3,
        )
        patch.set_gid(f"ghost:{g.sse.id}")
        c.ax.add_patch(patch)
        if c.look.labels:
            t = c.ax.text(
                g.cx,
                g.cy,
                g.label,
                ha="center",
                va="center",
                fontsize=c.font,
                fontweight="bold",
                color=darken(color),
                zorder=4,
            )
            t.set_gid(f"ghostlabel:{g.sse.id}")


def _draw_termini(c: _Canvas, layout: Layout, sses: list[SSE]) -> None:
    """The N and C terminal stubs, their letters and long-tail lengths."""
    groups = group_by_chain(sses)
    over = c.backbone[1]
    for end, chain, port, ex in termini(layout, sses):
        tip = (port[0] + ex[0] * END_STUB, port[1] + ex[1] * END_STUB)
        ink = over.get(chain, c.look.ink)
        if (
            c.look.loop_color == "residue" and layout.res_chain and chain not in over
        ):  # a terminal tail: its typical colour
            mine = groups.get(chain, [])
            ours = [r for r, ch in enumerate(layout.res_chain) if ch == chain]
            if mine and ours:
                tail = range(ours[0], mine[0].start) if end == "N" else range(mine[-1].end + 1, ours[-1] + 1)
                ink = _typical([c.residue_colour(r, chain) for r in tail]) or c.look.ink
        if ink != c.look.ink and chain not in over:
            edge = PathPatch(Path([port, tip]), fc="none", ec=_RUN_EDGE, lw=c.loop_lw * _RUN_EDGE_W, zorder=1.99)
            edge.set_gid(f"stub-edge:{end}:{chain}")
            c.ax.add_patch(edge)
        stub = PathPatch(Path([port, tip], [Path.MOVETO, Path.LINETO]), fc="none", ec=ink, lw=c.loop_lw, zorder=2)
        stub.set_gid(f"stub:{end}:{chain}")
        c.ax.add_patch(stub)
        t = c.ax.text(
            port[0] + ex[0] * END_LABEL,
            port[1] + ex[1] * END_LABEL,
            end,
            ha="center",
            va="center",
            fontsize=c.font,
            fontweight="bold",
            zorder=4,
        )
        t.set_gid(f"terminus:{end}:{chain}")
        mine = groups.get(chain, [])
        ours = [r for r, ch in enumerate(layout.res_chain) if ch == chain] if layout.res_chain else []
        if mine and ours:  # a long unstructured end: say how long, so it is not mistaken for a short stub
            n_tail = (mine[0].start - ours[0]) if end == "N" else (ours[-1] - mine[-1].end)
            if n_tail >= LONG_TAIL:
                at = np.add(port, np.asarray(ex) * (END_LABEL + 0.62))
                t = c.ax.text(
                    *at,
                    f"{n_tail} aa",
                    ha="center",
                    va="center",
                    fontsize=c.font * 0.72,
                    color=_legible("#5d6877"),
                    zorder=4,
                )
                t.set_gid(f"tail-length:{end}:{chain}")


_NA_ENTRY = re.compile(r"(DNA|RNA|DNA/RNA) \(")  # the legend text of a nucleic-acid backbone


def _draw_legend(c: _Canvas, layout: Layout, legend: _Legend) -> None:
    """The legend rows and, for a large assembly, the context caption under them."""
    xmin, xmax, ymin, ymax = c.page
    for k, (colour, text) in enumerate(legend.entries):
        row, col = divmod(k, legend.per_row)
        x, y = xmin + c.pad + col * legend.entry, ymin + legend.height - 1.1 - row * _LEGEND_ROW
        if isinstance(colour, tuple) and colour[0] == "grades":  # the nine grades in ConSurf's colours
            for q, shade in enumerate(CONSURF):
                c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none"))
            c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5))
        elif isinstance(colour, tuple):  # a property ramp
            for q in range(6):
                c.ax.add_patch(
                    Rectangle((x + q * 0.1, y - 0.3), 0.1, 0.6, fc=property_colour(q / 5, colour[1]), ec="none")
                )
        elif colour is None:  # the N -> C ramp
            for q in range(6):
                c.ax.add_patch(
                    Rectangle(
                        (x + q * 0.1, y - 0.3), 0.1, 0.6, fc=sequence_colour(q / 5, c.look.sequence_map), ec="none"
                    )
                )
        elif text.startswith("H-bonds"):  # a short dashed line, as drawn
            (dash,) = c.ax.plot([x, x + 0.6], [y, y], color=colour, lw=c.loop_lw * 0.85, ls=(0, (3, 2)))
            dash.set_gid("legend-hbond")
        elif _NA_ENTRY.match(text):  # a nucleic-acid backbone: a short thick stroke, not a chain's square
            (rail,) = c.ax.plot([x, x + 0.6], [y, y], color=colour, lw=c.lw * 2.4, solid_capstyle="butt")
            rail.set_gid("legend-na")
        elif text == TIE_LEGEND:  # two short ticks, each with its letter, as a tie is drawn
            (ticks,) = c.ax.plot(
                [x + 0.1, x + 0.1, np.nan, x + 0.5, x + 0.5],
                [y - 0.3, y, np.nan, y - 0.3, y],
                color=colour,
                lw=c.lw * 1.2,
                solid_capstyle="round",
            )
            ticks.set_gid("legend-tie")
            for n, tx in enumerate((x + 0.1, x + 0.5)):
                a = c.ax.text(
                    tx,
                    y + 0.22,
                    "a",
                    ha="center",
                    va="center",
                    fontsize=c.font * _TIE_LABEL,
                    fontweight="bold",
                    color=_legible(colour),
                )
                a.set_gid(f"legend-tie-label:{n}")
        elif text.startswith("chromophore"):  # the mark itself: a star
            star = Polygon(_star(x + 0.3, y, 0.42), closed=True, fc=colour, ec=LIGAND_EDGE, lw=c.lw * 0.7)
            star.set_gid("legend-linked")
            c.ax.add_patch(star)
        else:
            swatch = Rectangle((x, y - 0.3), 0.6, 0.6, fc=colour, ec=darken(colour), lw=c.lw * 0.8)
            c.ax.add_patch(swatch)
            if text.startswith("Chain "):  # a chain's legend entry: clickable in the interactive page
                swatch.set_gid(f"legend-chain:{text[6:]}")
        t = c.ax.text(x + 0.85, y, text, ha="left", va="center", fontsize=c.font)
        if text.startswith("Chain "):
            t.set_gid(f"legend-chain-label:{text[6:]}")
    if layout.context and c.look.legend:  # what the figure shows of a large assembly, under the legend
        import textwrap

        per_line = max(30, int((xmax - xmin - 2 * c.pad) * c.pt_per_unit / (c.font * 0.85 * 0.52)))
        t = c.ax.text(
            xmin + c.pad,
            ymin + 0.55,
            "\n".join(textwrap.wrap(layout.context, per_line)),
            ha="left",
            va="center",
            fontsize=c.font * 0.85,
            style="italic",
            color=_legible("#5d6877"),
        )
        t.set_gid("context")


def _draw_title(c: _Canvas, title: str) -> None:
    """The title, centred at the top."""
    xmin, xmax, ymin, ymax = c.page
    c.ax.text((xmin + xmax) / 2, ymax - 0.8, title, ha="center", va="center", fontsize=c.font * 1.3, fontweight="bold")


def _draw(
    layout: Layout,
    loops: list[Loop],
    sses: list[SSE],
    title: str | None = None,
    *,
    look: Style | None = None,
    palette: str | None = None,
    style: str | None = None,
    loop_style: str | None = None,
    model=None,
) -> Figure:
    look = _resolve_look(look, palette, style, loop_style)
    chains = _chains_in_order(sses)
    colors = chain_colors(chains, look.palette)
    element_colour, residue_colour = element_colours(layout, sses, look)
    if not layout.placed and not layout.stems and not layout.na_paths:
        return _empty_figure()
    marks = ligand_marks(layout, loops, sses, look)
    linked = linked_marks(layout, loops, sses, marks)
    umarks = uniprot_marks(layout, loops, sses, look)
    hb = hbond_overlay(layout, loops, sses, look, model)
    c, legend = _canvas(
        layout, look, title, chains, colors, element_colour, residue_colour, marks, umarks, hb, linked, loops
    )
    _draw_loops(c, layout, loops, sses)
    _draw_membrane(c, layout)
    _draw_domains(c, layout)
    _draw_sheet_panels(c, layout)
    _draw_hbonds(c, hb)
    c.fig.foldmap_hbonds = hb  # for the interactive page's hover notes
    if layout.stems or layout.na_paths:
        _draw_nucleic(
            c.ax, layout, sses, c.residue_colour, c.look.fill, c.lw, c.font, c.pt_per_unit, c.look, c.backbone[0]
        )
    _draw_elements(c, layout)
    _draw_links(c.ax, layout, loops, sses, look, c.lw, c.pt_per_unit)
    _draw_ligands(c.ax, marks, c.lw, c.font)
    _draw_linked(c.ax, linked, c.lw, c.font)
    _draw_uniprot(c.ax, umarks, layout, c.lw, c.font)
    _draw_residue_numbers(c, layout)
    _draw_ghosts(c, layout)
    _draw_termini(c, layout, sses)
    _draw_legend(c, layout, legend)
    if title:
        _draw_title(c, title)
    return c.fig


def save_svg(fig: Figure) -> str:
    """The figure as an inline SVG element (no XML prolog), with the same settings as save()."""
    import io as _io

    buf = _io.StringIO()
    with rc_context({"svg.fonttype": "none", "svg.hashsalt": "foldmap", **_RC}):
        fig.savefig(buf, format="svg", metadata={"Date": None}, facecolor=fig.get_facecolor())
    text = buf.getvalue()
    return text[text.index("<svg") :]


def save(fig: Figure, path: str | FsPath, dpi: int = 300) -> FsPath:
    path = FsPath(path)
    fmt = _EXT.get(path.suffix.lower())
    if fmt is None:
        raise ValueError(f"unsupported extension {path.suffix!r}; use .svg, .pdf or .png")
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {"svg": {"Date": None}, "pdf": {"CreationDate": None}, "png": {"Software": None}}[fmt]
    with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42, "svg.hashsalt": "foldmap", **_RC}):
        fig.savefig(path, format=fmt, dpi=dpi, metadata=meta, facecolor=fig.get_facecolor())
    return path
