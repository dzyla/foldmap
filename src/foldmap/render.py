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
from types import SimpleNamespace
from typing import NamedTuple

import numpy as np
from matplotlib import colormaps, rc_context
from matplotlib.axes import Axes
from matplotlib.colors import to_hex, to_rgba
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path
from matplotlib.text import Text
from matplotlib.textpath import text_to_path
from matplotlib.transforms import Affine2D

from . import text_metrics
from .hbonds import find_hbonds
from .interactions import (
    Interaction,
    failing_lines,
    find_cation_pi,
    find_hydrophobic,
    find_pi_stacking,
    find_salt_bridges,
    sort_key,
)
from .layout import (
    _GLYPH_W,
    AXIS_SYMBOL_R,
    DNA_MINOR,
    DNA_PITCH,
    DNA_TUBE,
    DNA_W,
    DOMAIN_NAME_SIZE,
    END_STUB,
    HEAD,
    HEAD_LEN,
    IN_PER_UNIT,
    NA_RUN_W,
    NT_ROW,
    PITCH,
    RUN_EDGE_W,
    SHAFT,
    TAIL_COUNT_SIZE,
    Layout,
    Placed,
    PlacedStem,
    _ink_em,
    _ink_line,
    _reach,
    _snap,
    _text_width,
    block_ghosts,
    centred_box,
    clear_inks,
    count_half,
    count_page,
    dna_to_page,
    dna_y,
    domain_fields,
    end_label_pos,
    first_pair_text,
    font_pt,
    glyph,
    glyph_scale,
    helix_label_pos,
    label_boxes,
    line_band,
    line_pt,
    lines_meeting,
    na_strand_lines,
    na_text_boxes,
    overlaps_any,
    page_fit,
    page_scale,
    page_text,
    panel_members,
    resolve,
    segments_of,
    spans_crossed,
    stem_keep_box,
    stem_label_pos,
    tail_counts,
    termini,
    text_anchor,
    text_em,
    thin_boxes,
    upright_deg,
)
from .model import SSE, group_by_chain
from .na_paths import _local, among, midway
from .palette import BACK_SHADE, PALETTES, chain_colors, darken, delta_e, text_color_on, tint
from .provenance import as_text, figure_provenance, with_desc
from .route import Loop, boxes_meeting, clip_len, clip_lengths
from .style import _ALIASES, CHOICES, MODALITIES, Style, backbone_overrides, interactions_of

_IN_PER_UNIT = IN_PER_UNIT
_HEAD_HALF = HEAD / 2
_HEAD_LEN = HEAD_LEN
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
_LIPID_ALPHA = 0.8  # the band over a sheet's plate lets the plate show through
_RUN_EDGE, _RUN_EDGE_W = "#3a3a3a", RUN_EDGE_W  # dark edge under residue-coloured loops (pale colours stay legible)
_SHORT_HELIX = 1.6  # helices/3-10 boxes shorter than this hold their residue numbers past the ends
LONG_LOOP = 25  # loops at least this many residues get a residue-count label (chain ends: layout.LONG_TAIL)
_DOMAIN_TONES = (  # one hue per domain: Paul Tol's muted set less its indigo (too near the RNA purple), ordered so
    # that neighbours in sequence differ most; colour-blind safe (tests/test_domains.py checks it)
    "#44aa99",  # teal
    "#cc6677",  # rose
    "#ddcc77",  # sand
    "#117733",  # green
    "#88ccee",  # cyan
    "#882255",  # wine
    "#999933",  # olive
    "#aa4499",  # purple
)
_DOMAIN_BAR_H = 1.1  # page units each chain's domain bar takes below the figure (at least; more for larger text)

Z = {  # the drawing order, back to front (plan 10 R8: panel < plate < membrane band < loops < bodies < bars/tethers
    # < labels).
    # Every zorder render sets is taken from here; a loop adds i * 1e-3 (its rank), its casing, arrows and runs a
    # few 1e-4 more, and a helix over a sheet lifts its parts by Z["over"] - Z["body"].
    "domain-panel": 0.2,  # a domain's tinted panel
    "sheet-panel": 0.3,  # a sheet's plate
    "sheet-domain-edge": 0.31,  # a domain's edge over a sheet several domains share
    "membrane": 0.312,  # the lipid band, translucent over the plates (an unrolled beta-barrel's runs across it) ...
    "membrane-heads": 0.315,  # ... and its head groups
    "domain-name-plate": 0.32,  # a domain's name: its plate over the panel's edge ...
    "domain-leader": 0.33,  # ... its leader, when it stands outside
    "domain-name": 0.35,  # ... the name
    "loop": 1.0,  # loops
    "na-tie": 1.8,  # an isolated pair drawn as a short tie
    "na-run": 1.9,  # a single nucleic-acid strand
    "stub-edge": 1.99,  # a chain end's stub: its edge ...
    "stub": 2.0,  # ... and line
    "na-back": 2.0,  # a stem: its back rail ...
    "na-rung": 2.1,  # ... its rungs ...
    "na-front": 2.2,  # ... its front rail or block
    "hbond": 2.15,  # H-bond overlay lines
    "na-contact": 2.3,  # contact beads on a nucleic acid
    "na-tie-mark": 2.4,  # a long tie's marks
    "domain-bar": 2.0,  # the domain bar under the figure ...
    "domain-bar-label": 2.1,  # ... and the names in it
    "body-back": 2.9,  # a helix's back face (Z["body"] - 0.1)
    "body": 3.0,  # strands, helices, 3-10 boxes, ghosts
    "gloss": 3.05,  # a helix's shine (Z["body"] + 0.05)
    "over-loop": 3.2,  # the loops of a helix over a sheet: above the strands, below the helix
    "over": 3.5,  # a helix over a sheet (Placed.over): above the strands, below their letters
    "tether": 3.6,  # ligand tethers and disulfide bars: above the bodies (R8), but drawn only where no body or text
    "bar": 3.7,  # is (visible_pieces), so they read as passing under them (plan 10 Task 6, Q-R10)
    "bar-end": 3.75,  # a disulfide's end dots and end marks (on a body's edge), glycan stems
    "mark": 3.8,  # ligand, linked and glycan marks
    "feature-ring": 3.84,  # a UniProt feature: its ring ...
    "feature": 3.85,  # ... and dot
    "mark-label": 3.9,  # ligand and linked names
    "label": 4.0,  # element labels, stem names, letters, chain ends, loop lengths
    "resnum": 4.5,  # residue numbers at element ends (and a helix-over-sheet's name)
}

_BAR_TICK = 0.3  # page units from a narrow segment to its name below it
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
def _turn(band: float) -> tuple[np.ndarray, ...]:
    """Front band, back band, front head, front tail and glint of one turn, for a ribbon `band` of a turn wide."""
    return (
        _band(-0.25, 0.25, width=band),
        _band(0.25, 0.75, width=band),
        _band(0.0, 0.25, 9, taper=+1, width=band),
        _band(-0.25, 0.0, 9, taper=-1, width=band),
        _band(-0.16, 0.16, 14, taper=2, width=band * _SHINE) + (band * 0.1, 0),  # glint on the near face
    )


_TURN = _turn(_BAND)


def _coil(p: Placed, residues: int) -> tuple[Path, Path, Path]:
    """Front faces, back faces and front highlights of a helix as compound paths in page coordinates."""
    turns = int(max(1, min(round(residues / _RES_PER_TURN), p.length // _MIN_PITCH)))
    return helix_ribbon(p.cx, p.cy, p.length, p.width, p.angle, turns)


def helix_ribbon(
    cx: float, cy: float, length: float, width: float, angle: float, turns: int, band: float = _BAND
) -> tuple[Path, Path, Path | None]:
    """Front faces, back faces and front highlights (None for one turn) of a helix ribbon of `turns` turns, centred
    on (cx, cy), `length` along its axis at `angle`, `width` across, as compound paths in page coordinates; it starts
    and ends on its axis. band: the ribbon's width along the axis, as a fraction of one turn. The figure's helices
    and the sequence view's track are both drawn with it."""
    front_band, back_band, head, tail, glint = _TURN if band == _BAND else _turn(band)
    front = [head] + [front_band + (k, 0) for k in range(1, turns)] + [tail + (turns, 0)]
    back = [back_band + (k, 0) for k in range(turns)]
    shine = [glint + (k, 0) for k in range(1, turns)]
    to_page = Affine2D().translate(-turns / 2, 0).scale(length / turns, width).rotate(angle).translate(cx, cy)

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
    Every colour goes through `legible`: on a dark page a user-given colour that is too dark is lightened (on a
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
        na[chain] = legible(over.get(chain, colour))
    loops = {c: legible(over[c] if c in over else over["protein"]) for c in protein if c in over or "protein" in over}
    return na, loops


_DNA_TUBE = DNA_TUBE  # backbone tube width, page units
_DNA_SEG = 16  # points along each half-turn


def _tube(pts: np.ndarray, width: float) -> np.ndarray:
    """A centre line thickened to a band of `width` (perpendicular offsets)."""
    t = np.gradient(pts, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    n = np.column_stack([-t[:, 1], t[:, 0]]) * width / 2
    return np.vstack([pts + n, (pts - n)[::-1]])


def _dna_turns(d: PlacedStem, slot: int, wide: float = 1.0) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Front and back half-turns of one backbone in duplex-local coordinates: the same half-turn fragment
    repeated every turn (shifted by the minor groove for the second backbone), clipped to the duplex ends; the
    tube `wide` times _DNA_TUBE (Style.na_line)."""
    x0, x1 = -d.length / 2, d.length / 2
    off = slot * DNA_MINOR * DNA_PITCH
    front, back = [], []
    for k in range(int(np.floor((x0 - off) / DNA_PITCH)) - 1, int(np.ceil((x1 - off) / DNA_PITCH)) + 1):
        for lo, hi, dest in ((-0.25, 0.25, front), (0.25, 0.75, back)):
            a, b = max(x0, off + (k + lo) * DNA_PITCH), min(x1, off + (k + hi) * DNA_PITCH)
            if b - a > 1e-6:
                xs = np.linspace(a, b, _DNA_SEG)
                dest.append(_tube(np.column_stack([xs, dna_y(slot, xs, d.width)]), _DNA_TUBE * wide * d.width / DNA_W))
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
    em: float  # page units of one em of the text (font / points per page unit)
    wide: float  # the backbone's weight (Layout.na_line): rails, runs, ribbon tubes and block edges

    def colour(self, strand: int) -> str:
        """The backbone colour of a strand: stem rails, ribbons, blocks, paths, hairpins, ties."""
        return self.backbone[self.chains[strand]]


_upright = upright_deg


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
            zorder=Z["na-rung"],
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
            lw=ink.lw * 2.2 * ink.wide,
            capstyle="round",
            joinstyle="round",
            zorder=Z["na-front"],
        )
        rail.set_gid(f"na:stem:{k}:rails" if n == 0 else f"na:stem:{k}:rails:{n + 1}")
        ink.ax.add_patch(rail)
    _rung_lines(ink, k, d)
    for key in d.axial:
        letters[key] = (d.point(key, NT_ROW), _upright(d.angle))


def _draw_stem_label(ink: _NAInk, k: int, d: PlacedStem) -> None:
    t = ink.ax.text(
        *text_anchor(d.name, True, stem_label_pos(d), ink.em * 0.8, _upright(d.angle)),
        d.name,
        ha="left",
        va="baseline",
        rotation=_upright(d.angle),
        rotation_mode="anchor",  # started where the layout centres it (text_anchor), as any matplotlib draws it
        fontsize=ink.font * 0.8,
        fontweight="bold",
        color=ink.front,
        zorder=Z["label"],
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
            f, b = _dna_turns(d, slot, ink.wide)
            fronts += f
            backs += b
        suffix = "" if n == 0 else f":{n + 1}"
        back = PathPatch(_compound([page(q) for q in backs]), fc=tint(c, 0.5), ec="none", zorder=Z["na-back"])
        back.set_gid(f"na:stem:{k}:back{suffix}")
        ink.ax.add_patch(back)
        front = PathPatch(
            _compound([page(q) for q in fronts]),
            fc="white" if outline else c,
            ec=c if outline else "none",
            lw=ink.lw * 0.8,
            zorder=Z["na-front"],
        )
        front.set_gid(f"na:stem:{k}:front{suffix}")
        ink.ax.add_patch(front)
    rungs = [
        dna_to_page(d, [x, x], [dna_y(0, x, d.width), dna_y(1, x, d.width)])
        for x in (d.axial[(a, i)] for a, i, _, _ in d.pairs)
    ]
    rung = PathPatch(
        _compound(rungs, closed=False), fc="none", ec=ink.rung, lw=ink.lw * 2.2, capstyle="round", zorder=Z["na-rung"]
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
            corners,
            closed=True,
            fc="white" if ink.fill == "outline" else c,
            ec=c,
            lw=ink.lw * 1.2 * ink.wide,
            zorder=Z["na-front"],
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
        color=legible(mixed) if ink.fill == "outline" else text_color_on(mixed),
        zorder=Z["label"],
    )
    t.set_gid(f"na:label:{k}")
    first = first_pair_text(d, ink.chains, ink.seqs)  # the outermost pair's residue numbers, in the name's slot
    t = ink.ax.text(
        *text_anchor(first, False, stem_label_pos(d), ink.em * 0.75, _upright(d.angle)),
        first,
        ha="left",
        va="baseline",
        rotation=_upright(d.angle),
        rotation_mode="anchor",
        fontsize=ink.font * 0.75,  # a base letter's size: the stem's name slot, empty in this style
        color=ink.front,
        zorder=Z["label"],
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
            lw=ink.lw * (1.0 if gap else NA_RUN_W * ink.wide),
            ls=((0, (2, 2)) if path.label else (0, (1, 1.5))) if gap else "-",
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=Z["na-run"],
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
                zorder=Z["na-tie"],
            )
            line.set_gid(gid)
            continue
        for a, b, at in tie.marks:
            (tick,) = ink.ax.plot(
                *np.c_[a, b],
                color=ink.colour(tie.a[0]),
                lw=ink.lw * 1.2,
                solid_capstyle="round",
                zorder=Z["na-tie-mark"],
            )
            tick.set_gid(f"{gid}:mark")
            t = ink.ax.text(
                *at,
                tie.label,
                ha="center",
                va="center",
                fontsize=ink.font * _TIE_LABEL,  # a base letter's size, upright, bold, in the mark's colour
                fontweight="bold",
                color=legible(ink.colour(tie.a[0])),
                zorder=Z["label"],
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
            base = na.strands[key[0]].bases[key[1]]
            t = ink.ax.text(
                *text_anchor(base, bool(hit), xy, ink.em * 0.75, rot),  # centred on xy, as the layout measures it
                base,
                ha="left",
                va="baseline",
                rotation=rot,
                rotation_mode="anchor",
                fontsize=ink.font * 0.75,
                fontweight="bold" if hit else "normal",
                color=hit[1] if hit else ink.colour(key[0]) if style == "letters" else ink.front,
                zorder=Z["label"],
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
                    zorder=Z["na-contact"],
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
        font / pt_per_unit,
        layout.na_line,
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
    k = p.width / PITCH  # Style.strand_scale, carried by the strand's width
    return strand_arrow(p.cx, p.cy, p.length, p.angle, k)


def strand_arrow(
    cx: float, cy: float, length: float, angle: float, k: float = 1.0, head: bool = True, tip: float = _HEAD_LEN
) -> np.ndarray:
    """A strand's arrow, centred on (cx, cy), `length` long at `angle`, its shaft and head k times as wide as the
    figure's (Style.strand_scale), the head `tip` long (at most 0.45 of the arrow) and ending at the far end;
    head=False: the shaft alone (a strand that runs on past a row of the sequence view). The figure's strands and
    the sequence view's track are both drawn with it."""
    half, sh = length / 2, SHAFT * k / 2
    if not head:
        local = np.array([(-half, -sh), (half, -sh), (half, sh), (-half, sh)])
    else:
        tip, hh = min(tip, length * 0.45), _HEAD_HALF * k
        local = np.array(
            [
                (-half, -sh),
                (half - tip, -sh),
                (half - tip, -hh),
                (half, 0),
                (half - tip, hh),
                (half - tip, sh),
                (-half, sh),
            ]
        )
    c, s = np.cos(angle), np.sin(angle)
    return local @ np.array([[c, s], [-s, c]]) + [cx, cy]


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


def _relative_luminance(colour) -> float:
    """WCAG 2 relative luminance (linear-light sRGB)."""
    from matplotlib.colors import to_rgb

    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(colour)]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def contrast_ratio(a, b) -> float:
    """WCAG 2 contrast ratio of two colours, 1 to 21 (4.5 is the AA level for normal text)."""
    la, lb = sorted((_relative_luminance(a), _relative_luminance(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


_NAME_CONTRAST = 4.5  # WCAG AA for normal text: a domain's name against what it is drawn on


def ink_on(fill: str) -> str:
    """Black or white, whichever stands out more on `fill` (always 4.5:1 or more)."""
    return "#000000" if contrast_ratio("#000000", fill) >= contrast_ratio("#ffffff", fill) else "#ffffff"


def domain_name_colour(tone: str, on: str) -> str:
    """A domain's name in its own hue as drawn on `on` (its panel's tint, or the page): a deeper shade on a light
    fill, the tone itself on a dark one, moved toward black (or white) until it reaches _NAME_CONTRAST."""
    start = darken(tone, 0.75) if ink_on(on) == "#000000" else tone
    for k in range(21):
        colour = _mix(start, ink_on(on), k / 20)
        if contrast_ratio(colour, on) >= _NAME_CONTRAST:
            return colour
    return ink_on(on)


def _label_ink(colour: str) -> str:
    """Text in an element's own hue: a deeper shade on light pages, the colour itself on dark ones."""
    if _luminance(_BACKGROUND.get()) < 0.5:
        return legible(colour)
    return legible(darken(colour, 0.8))


def legible(colour: str, ceiling: float = 0.5) -> str:
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


_RESIDUE_STEP = 0.5  # page units along a loop's route per residue from the nearer element end (R3)


def residue_point(layout: Layout, loops: list[Loop], sses: list[SSE], bb, r: int) -> tuple[float, float] | None:
    """Where residue r (Backbone index) sits on the page: along its element, or along the loop it is in,
    the number of residues from the nearer element end (a bounded step each, _RESIDUE_STEP); past a chain end, at
    the end of the terminal stub."""
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
            run = np.cumsum(seg)
            d_a, d_b = r - a.end, b.start - r
            # residue distance from the nearer end, a bounded step along the route from that end (R3: never a
            # fraction of the routed length, which a loop round a barrel stretches far from either element); a
            # short route keeps the residues in order inside it
            walk = min(min(d_a, d_b) * _RESIDUE_STEP, run[-1] * min(d_a, d_b) / max(d_a + d_b, 1))
            target = walk if d_a <= d_b else run[-1] - walk
            k = min(int(np.searchsorted(run, target)), len(seg) - 1)
            into = (target - (run[k] - seg[k])) / max(seg[k], 1e-9)
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


def _nearest_on(pts: np.ndarray, q, avoid=()) -> tuple[float, float]:
    """The point of polyline pts closest to q, kept clear of the loop's two end stubs, and outside everything in
    `avoid` (boxes and outlines: _covered) where any of the polyline is (sampled every _HIDE_STEP there)."""
    q = np.asarray(q, float)
    best, d_best = None, np.inf
    for a, b in zip(pts, pts[1:]):
        ab = b - a
        t = float(np.clip(np.dot(q - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0))
        c = a + t * ab
        d = float(np.hypot(*(c - q)))
        if d < d_best:
            best, d_best = c, d
    if avoid and _covered(best[None, :], avoid)[0]:
        dense = _dense(pts)[0]
        free = ~_covered(dense, avoid)
        if free.any():
            d = np.hypot(*(dense - q).T)
            best = dense[int(np.argmin(np.where(free, d, np.inf)))]
    return float(best[0]), float(best[1])


def _covered(pts: np.ndarray, avoid) -> np.ndarray:
    """Which of the points lie inside any of `avoid`: boxes (x0, y0, x1, y1) and outlines (Path; each asked only of
    the points within its vertices' extent, _within_extents)."""
    shapes = [a for a in avoid if isinstance(a, Path)]
    out = _in_boxes(pts, [a for a in avoid if not isinstance(a, Path)])
    inside = _within_extents(pts, shapes)
    for k in np.flatnonzero(inside.any(axis=0)).tolist():
        near = inside[:, k] & ~out
        if near.any():
            out[near] |= shapes[k].contains_points(pts[near])
    return out


def _within_extents(pts: np.ndarray, paths: list[Path]) -> np.ndarray:
    """(points x paths) which points lie within each path's vertices' extent, edges included: only those can lie
    inside it (plan 15 Task 7: 8D3D's focus view, 848 disulfides under about 1,800 bodies, asked every body of
    every point)."""
    out = np.zeros((len(pts), len(paths)), bool)
    sizes = np.array([len(p.vertices) for p in paths], int)
    some = np.flatnonzero(sizes > 0)  # a path of no vertices holds no point (and has no extent to reduce)
    if not len(some):
        return out
    v = np.concatenate([paths[k].vertices for k in some.tolist()])
    starts = np.cumsum(sizes[some]) - sizes[some]
    lo, hi = np.minimum.reduceat(v, starts), np.maximum.reduceat(v, starts)
    x, y = pts[:, :1], pts[:, 1:2]
    out[:, some] = (x >= lo[:, 0]) & (x <= hi[:, 0]) & (y >= lo[:, 1]) & (y <= hi[:, 1])
    return out


def disulfide_points(layout, loops, sses, chains, i: int, j: int, avoid=()):
    """Where to draw a disulfide's two cysteines. A cysteine in a loop is placed where its loop passes closest
    to its partner, outside `avoid` (text boxes, body outlines) where its loop leaves room: loop routes are
    schematic, so any point on its own loop is faithful, and the nearest one keeps the bond short. Cysteines in
    elements stay at their place along the element."""
    a, b = residue_point(layout, loops, sses, chains, i), residue_point(layout, loops, sses, chains, j)
    if a is None or b is None:
        return None, None
    la, lb = loop_of(layout, loops, sses, chains, i), loop_of(layout, loops, sses, chains, j)
    for _ in range(3):  # both in loops: settle by alternating
        if la is not None:
            a = _nearest_on(la, b, avoid)
        if lb is not None:
            b = _nearest_on(lb, a, avoid)
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
    patch.set(facecolor=colour, edgecolor=legible("#2b2b2b"), linewidth=lw * 0.6, zorder=Z["mark"])
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
_SPOT_FAR = 3  # ... and, when nothing is clear within it, at this many times it (then the mark goes on its anchor)
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


_SPOT_MORE_RINGS = 4  # rings past the first clear spot's where a spot whose tethers cross fewer bodies is looked for


def _over_bodies(c, targets, bodies) -> float:
    """How much of the tethers from spot c runs over element bodies (`bodies`, their boxes): each tether counted
    over every body but the one its target lies on (random-100 C9)."""
    boxes = np.asarray(bodies, float).reshape(-1, 4)
    if not len(boxes) or not len(targets):
        return 0.0
    t = np.asarray(targets, float).reshape(-1, 2)
    lengths = clip_lengths(np.hstack([np.broadcast_to(np.asarray(c, float), t.shape), t]), boxes)
    own = (boxes[None, :, 0] <= t[:, 0:1]) & (t[:, 0:1] <= boxes[None, :, 2])  # the body a target lies on
    own &= (boxes[None, :, 1] <= t[:, 1:2]) & (t[:, 1:2] <= boxes[None, :, 3])
    return float(np.where(own, 0.0, lengths).sum())


_ALONG_NEAR, _ALONG_TURN = (
    0.15,
    np.radians(20),
)  # a tether this near a loop line and this near its direction runs along it


def _along_loops(c, targets, runs) -> float:
    """How much of the tethers from spot c runs along a loop's line (`runs`, segments x0 y0 x1 y1): the length of
    each tether beside a segment within _ALONG_TURN of its direction whose two ends lie within _ALONG_NEAR of the
    tether's line (random-100 C9: 5HI7's Zn tether drawn on the loop it reaches)."""
    segs = np.asarray(runs, float).reshape(-1, 4)
    if not len(segs):
        return 0.0
    c = np.asarray(c, float)
    a, b = segs[:, :2], segs[:, 2:]
    e = b - a
    el = np.hypot(e[:, 0], e[:, 1])
    total = 0.0
    for t in np.asarray(targets, float).reshape(-1, 2):
        d = t - c
        length = float(np.hypot(d[0], d[1]))
        if length < 1e-9:
            continue
        u = d / length
        sin = np.abs(u[0] * e[:, 1] - u[1] * e[:, 0]) / np.maximum(el, 1e-12)
        ra, rb = a - c, b - c
        off_a = np.abs(ra[:, 0] * u[1] - ra[:, 1] * u[0])
        off_b = np.abs(rb[:, 0] * u[1] - rb[:, 1] * u[0])
        pa, pb = ra @ u, rb @ u
        lo, hi = np.clip(np.minimum(pa, pb), 0.0, length), np.clip(np.maximum(pa, pb), 0.0, length)
        beside = (el > 1e-9) & (sin < np.sin(_ALONG_TURN)) & (off_a <= _ALONG_NEAR) & (off_b <= _ALONG_NEAR)
        total += float(np.where(beside, hi - lo, 0.0).sum())
    return total


_RAIL_GRAZE = 0.3  # page units a tether may run inside a stem past the entry its target needs and not count as crossing


def _stems_crossed(c, targets, stems) -> int:
    """(tether, stem) pairs where the tether from spot c runs inside the stem's true (tilted) body (`stems`,
    PlacedStem) more than _RAIL_GRAZE beyond the entry it needs: a tether to a nucleotide on a stem's rail enters that
    stem by the rail's inset only; one from the far side crosses the ladder (random-100 C2b-4 fix round 1: 3OWW's GLY
    A89, after the stems were spread, tethered across P3b's rungs)."""
    n = 0
    c = np.asarray(c, float)
    for t in np.asarray(targets, float).reshape(-1, 2):
        lo, hi = np.minimum(c, t), np.maximum(c, t)
        for d in stems:
            x0, y0, x1, y1 = d.rect
            if x0 > hi[0] or lo[0] > x1 or y0 > hi[1] or lo[1] > y1:
                continue
            hl, hw = d.length / 2, d.width / 2
            inside = clip_len(_local(c, d), _local(t, d), (-hl, -hw, hl, hw))
            u, v = _local(t, d)
            entry = max(min(hl - abs(u), hw - abs(v)), 0.0)  # how deep its own target lies in this stem
            n += inside - entry > _RAIL_GRAZE
    return n


def _strands_crossed(c, targets, strands) -> int:
    """(tether, nucleic-acid strand) pairs whose lines cross properly (a tether ending on a strand touches it, it does
    not cross it): `strands`, each drawn strand's segments as an (n, 4) array (hairpin loops, runs, joints, gaps, free
    ends, traces). Products written out elementwise, never through BLAS."""
    c = np.asarray(c, float)
    n = 0
    for t in np.asarray(targets, float).reshape(-1, 2):
        lo, hi = np.minimum(c, t), np.maximum(c, t)
        for seg in strands:
            if not len(seg):
                continue
            p, q = seg[:, 0:2], seg[:, 2:4]
            near = (np.minimum(p, q) <= hi).all(axis=1) & (np.maximum(p, q) >= lo).all(axis=1)
            if not near.any():
                continue
            p, q = p[near], q[near]

            def turn(o, e, w):
                return (e[..., 0] - o[..., 0]) * (w[..., 1] - o[..., 1]) - (e[..., 1] - o[..., 1]) * (
                    w[..., 0] - o[..., 0]
                )

            a, b = c[None, :], t[None, :]
            if ((turn(p, q, a) * turn(p, q, b) < -1e-12) & (turn(a, b, p) * turn(a, b, q) < -1e-12)).any():
                n += 1
    return n


def _beyond_reach(c, bound) -> bool:
    """Does any of `bound` lie farther than _TETHER_REACH from spot c (the distances snapped, layout._snap)?"""
    return any(float(_snap(np.hypot(*np.subtract(t, c)))) > _TETHER_REACH for t in bound)


def _mark_spot(
    anchor,
    half,
    targets,
    taken,
    texts,
    tethers,
    bodies=(),
    radius=_SPOT_RADIUS,
    runs=(),
    na=None,
    retarget=None,
    bound=(),
):
    """The first spot on rings round the anchor whose mark (grown by 0.08) is clear of `taken` and of the tethers
    drawn before, and whose own tethers cross no text, lengthening none past the first clear spot's (the plan-11
    spot); where its tethers run over element bodies (`bodies`), the spot within _SPOT_MORE_RINGS rings further whose
    tethers run over the least of them (C9: straight tethers cross strands and helices) and along the least of the
    loops' lines (`runs`, _along_loops); else that first clear spot; None when there is none within _SPOT_RADIUS.
    `na` (stems, strands: a mark held by nucleotides alone) judges every clear spot within `radius`, its tethers no
    longer held to the first clear spot's reach, in this order: the stems its tethers cross (past the rail each
    enters, _stems_crossed), the strands they cross as drawn among the elements (hairpin loops, runs,
    _strands_crossed), their length over elements and along loop lines (as above), their length (random-100 C2b-4
    fix rounds 2 and 3 rulings: scored on the stems alone, 8SA2's B1Z was blind to the hairpin loop and runs its
    tethers crossed; the two counts come before the continuous lengths, which almost never tie).
    `retarget` (spot -> targets): the tethers a spot gets, when they depend on it (a mark held by nucleotides is
    tethered to those drawn nearest it, _na_targets), judged as drawn.
    `bound` (the home chain's main contact and the targets on other chains, _place_mark): a spot is taken only when
    each of them lies within _TETHER_REACH of it; its other home tethers are not held to the plan-11 reach (the
    longest tether from the first clear spot): they are drawn whatever their length (the plan-11 ruling), and 8ESZ's
    3PE 5604 has no spot keeping its B2 target within _TETHER_REACH and its far home target within that reach
    (18.8 apart > 12.10 + 6.5). None is returned when no spot within `radius` meets the bound (never the first clear
    spot)."""
    spot = free = None
    best, last_ring = None, None
    for n, ring in enumerate(np.arange(0.0, radius, 0.35)):
        if last_ring is not None and n > last_ring:
            break
        for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False):
            c = anchor + ring * np.array([np.cos(angle), np.sin(angle)])
            box = (c[0] - half[0] - 0.08, c[1] - half[1] - 0.08, c[0] + half[0] + 0.08, c[1] + half[1] + 0.08)
            if overlaps_any(box, taken) or any(clip_len(a, b, box) > 1e-9 for a, b in tethers):
                continue
            ts = targets if retarget is None else retarget(c)
            if free is None:
                free = c
                if not bound:  # a bounded search holds the other home tethers to no reach
                    reach = max(_TETHER_REACH, *(np.hypot(*(np.subtract(t, c))) for t in ts))
            if _beyond_reach(c, bound):
                continue
            if (na is not None or bound or all(np.hypot(*(np.subtract(t, c))) <= reach for t in ts)) and not any(
                clip_len(c, t, b) > 1e-9 for t in ts for b in texts
            ):
                over = _over_bodies(c, ts, bodies) + _along_loops(c, ts, runs)
                if na is not None:  # a tuple, compared in order, each part on a fixed grid (determinism)
                    length = sum(float(np.hypot(*(np.subtract(t, c)))) for t in ts)
                    over = (
                        _stems_crossed(c, ts, na[0]),
                        _strands_crossed(c, ts, na[1]),
                        round(over, 4),
                        round(length, 4),
                    )
                    if best is None or over < best[0]:
                        best = (over, c)
                else:
                    if best is None or over < best[0] - 1e-9:
                        best = (over, c)
                    if over <= 1e-9:
                        break
                if last_ring is None and na is None:
                    last_ring = n + _SPOT_MORE_RINGS
        if best is not None and na is None and best[0] <= 1e-9:
            break
    spot = best[1] if best is not None else None
    return spot if spot is not None or bound else free


def _place_mark(
    anchor, half, targets, taken, lines, texts, tethers, bodies=(), runs=(), na=None, retarget=None, bound=()
):
    """A ligand's spot: clear of the loops' lines and the tethers drawn before where there is room, else of the
    elements and text only (the tethers' crossings given up), else of the elements and text with every tether given
    up; else further out, else on its anchor: a chosen ligand is never left without a mark (C5 fix round 1: 6JO6 lost
    24 on a crowded page; batch B: 7QUP lost one with nothing clear within _SPOT_RADIUS).
    `bound`: the home chain's main contact (its heaviest target), then its targets on other chains. Each of those lies
    within _TETHER_REACH of the main contact (_local_targets), but the anchor, the mean of all its targets, follows
    home contacts drawn far from it: a spot that leaves one of them farther than _TETHER_REACH gives way to the first
    spot, kept off as much (the same keep-outs, or stricter), out to _SPOT_FAR times _SPOT_RADIUS, that keeps them
    and the main contact within it, so the mark stays in its pocket (random-100 C2b-5: 6JLN's CLA B609 drawn 10.15
    from H's helix, 3.80 from the main contact; held to the other chains alone, 1SLV's Cu went 7.35 from its anchor,
    7.67 from its main contact, in the beta view); else the spot stays."""
    layers = ((taken + lines, tethers), (taken, tethers), (taken, ()))

    def kept(spot, upto: int):
        if not _beyond_reach(spot, bound[1:]):
            return spot
        for blocked, drawn in layers[: upto + 1]:
            far = _SPOT_FAR * _SPOT_RADIUS
            near = _mark_spot(anchor, half, targets, blocked, texts, drawn, bodies, far, runs, na, retarget, bound)
            if near is not None:
                return near
        return spot

    for k, (blocked, drawn) in enumerate(layers):
        spot = _mark_spot(anchor, half, targets, blocked, texts, drawn, bodies, runs=runs, na=na, retarget=retarget)
        if spot is not None:
            return kept(spot, k)
    # nothing clear within _SPOT_RADIUS (7QUP's crowded page dropped GTP:9E500): further out, else the least bad
    spot = _mark_spot(
        anchor,
        half,
        targets,
        taken,
        texts,
        (),
        bodies,
        radius=_SPOT_FAR * _SPOT_RADIUS,
        runs=runs,
        na=na,
        retarget=retarget,
    )
    return kept(spot, 2) if spot is not None else _least_bad_spot(anchor, half, taken, bodies)


def _least_bad_spot(anchor, half, taken, bodies) -> np.ndarray:
    """A mark with no clear spot at all (ligand_marks records it as crowded): the ring spot within _SPOT_FAR times
    _SPOT_RADIUS overlapping the fewest element bodies, then the fewest other boxes, the nearest on a tie (fix round
    1: it went on its anchor, over whatever was there)."""
    body = np.asarray(bodies, float).reshape(-1, 4)
    other = np.asarray(taken, float).reshape(-1, 4)
    best, score = np.asarray(anchor, float), None
    for ring in np.arange(0.0, _SPOT_FAR * _SPOT_RADIUS, 0.35):
        for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False):
            c = np.asarray(anchor, float) + ring * np.array([np.cos(angle), np.sin(angle)])
            box = np.array([c[0] - half[0], c[1] - half[1], c[0] + half[0], c[1] + half[1]])

            s = (int(boxes_meeting(body, box).sum()), int(boxes_meeting(other, box).sum()))
            if score is None or s < score:
                best, score = c, s
    return best


def _disulfide_lines(layout: Layout, loops: list[Loop], sses: list[SSE], look: Style) -> list[tuple]:
    """Thin boxes along each disulfide bar's straight run between its two cysteines (its routed lane lies near it):
    ligand marks keep off them (C5)."""
    links = layout.links
    if links is None or not look.disulfides:
        return []
    chains = layout.res_chain
    out = []
    for i, j in links.disulfides:
        if layout.partial and (chains[i] in layout.partial or chains[j] in layout.partial):
            continue
        a, b = disulfide_points(layout, loops, sses, chains, i, j)
        if a is not None and b is not None:
            out += thin_boxes([a, b])
    return out


def _na_spots(layout: Layout, look: Style) -> tuple[dict, list, list]:
    """Where each drawn nucleotide sits (nucleotide_points, in the figure's stem style); the boxes a mark held by
    nucleotides keeps off: the stems, the strands' paths and the nucleic acid's text (letters and ends:
    layout.na_text_boxes; stem names with their margin, as routes keep off them: layout.stem_keep_box); and the boxes
    its tethers keep off, the stem names' own (label_boxes). A tether ends on a stem's rail, past the stem's box, its
    letters and the strand it touches: held to those too, no spot's tethers passed and the mark took the first clear
    spot whatever its tethers crossed (random-100 C2b-4 fix round 1: 3OWW's GLY A89 across P3b's rungs)."""
    if layout.nucleic is None:
        return {}, [], []
    paths = [box for p in layout.na_paths for box in thin_boxes(p.points)]
    paths += [stem_keep_box(d, glyph(layout)[0] / _GLYPH_W) for d in layout.stems]
    keep = [d.rect for d in layout.stems] + paths + na_text_boxes(layout)
    names = [box for name, box in label_boxes(layout, []) if name.startswith("stem-label:")]
    return nucleotide_points(layout, _na_style(layout, look)), keep, names


def _na_targets(layout: Layout, g, where: dict, near=None):
    """The tether targets of a ligand held by nucleotides alone (Ligand.nucleotides): one per stem (or strand) it
    sits on, the _LIGAND_TETHERS it touches most, at its middle nucleotide there, or with `near` (the mark's spot)
    the one drawn nearest it, so a tether to a stem ends on the rail facing the mark and never runs along a rung
    (8XZE's THG, batch C2a review Q5); weights: nucleotides held there; no protein residues."""
    na = layout.nucleic
    if na is None or not g.nucleotides:
        return [], [], []
    key = {
        (lab.chain, lab.seq, lab.icode): (s, i) for s, st in enumerate(na.strands) for i, lab in enumerate(st.labels)
    }
    stem = {k: d.id for d in layout.stems for k in d.axial}
    groups: dict = {}
    for lab in g.nucleotides:
        k = key.get(lab)
        if k is not None and k in where:
            groups.setdefault(stem.get(k, f"strand {k[0]}"), []).append(k)
    targets, weights = [], []
    for ks in sorted(groups.values(), key=len, reverse=True)[:_LIGAND_TETHERS]:
        if near is None:
            k = sorted(ks)[len(ks) // 2]
        else:  # the nearest, the first in order on a tie
            k = min(sorted(ks), key=lambda key: float(np.hypot(*(np.asarray(where[key], float) - near))))
        x, y = (float(v) for v in where[k])
        if all(np.hypot(x - t[0], y - t[1]) > 0.3 for t in targets):
            targets.append((x, y))
            weights.append(len(ks))
    return targets, weights, []


def _on_axis(layout: Layout, g) -> bool:
    """Whether a ligand sits on a rosette's axis (Layout.rosette): held by residues of every one of its copies."""
    n = getattr(layout, "rosette", 0)
    sector = getattr(layout.frame, "sector", None) or {}
    return bool(n) and len({sector.get(layout.res_chain[r]) for r in g.contacts} - {None}) == n


_AXIS_MARK_GAP = 0.12  # page units between two marks of a rosette's axis cluster


def _axis_cluster(ligands: list, taken: list, texts: list) -> list[dict]:
    """The marks of the ligands on a rosette's axis (_on_axis: held by every copy), one compact cluster round the axis
    symbol at the page centre, untethered (the cluster stands for the channel or pore every copy lines): on a square
    grid of one mark's size, the spots nearest the centre first (then by angle from page-right), each clear of the
    symbol, the elements, the names and the marks before it. In order of first holding residue; `taken` and `texts`
    gain each mark."""
    halves = [(0.3, 0.3) if g.metal else (0.12 * len(g.name) + 0.2, 0.27) for g in ligands]
    step = 2 * max(max(h) for h in halves) + _AXIS_MARK_GAP
    reach = int(np.ceil(np.sqrt(len(ligands)))) + 6
    spots = sorted(
        ((i * step, j * step) for i in range(-reach, reach + 1) for j in range(-reach, reach + 1)),
        key=lambda s: (round(float(np.hypot(*s)), 6), round(float(np.arctan2(s[1], s[0]) % (2 * np.pi)), 6)),
    )
    out = []
    order = sorted(range(len(ligands)), key=lambda i: ligands[i].contacts[0])
    for i in order:
        g, (hw, hh) = ligands[i], halves[i]
        spot = next(
            (s for s in spots if not overlaps_any((s[0] - hw, s[1] - hh, s[0] + hw, s[1] + hh), taken)), spots[0]
        )
        mark = (spot[0] - hw, spot[1] - hh, spot[0] + hw, spot[1] + hh)
        crowded = overlaps_any(mark, taken)
        taken.append(mark)
        texts.append(mark)
        out.append(
            {
                "ligand": g,
                "centre": (float(spot[0]), float(spot[1])),
                "half": (hw, hh),
                "targets": [],
                "residues": [],
                "crowded": bool(crowded),
            }
        )
    return out


def ligand_marks(layout: Layout, loops: list[Loop], sses: list[SSE], look: Style, counts=None) -> list[dict]:
    """Where each shown ligand/ion goes: as close as possible to the residues holding it, clear of elements,
    labels and other markers, and of the nucleic acid's stems, strands and text (_na_spots) where there is one. Each
    mark: ligand, centre, half size, tether targets, tethered residues."""
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
    if getattr(layout, "rosette", 0):  # a rosette's axis symbol at the page centre
        taken.append(centred_box(0.0, 0.0, AXIS_SYMBOL_R, AXIS_SYMBOL_R))
    bodies = list(taken)  # what a tether should not run over (C9)
    labels = [box for _, box in label_boxes(layout, sses)]
    if counts is None:
        counts = [t["box"] for t in [*tail_counts(layout, sses, loops), *loop_counts(layout, loops, sses)]]
    taken += labels + counts
    texts = (
        labels
        + counts
        + [  # and the strand letters at their centres (they slide off a disulfide bar only later)
            centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25)
            for p in layout.placed.values()
            if p.sse.kind == "E"
        ]
    )
    lines = [box for lp in loops for box in thin_boxes(lp.points)] + _disulfide_lines(layout, loops, sses, look)
    runs = [(*p, *q) for lp in loops for p, q in zip(lp.points, lp.points[1:])]  # a tether keeps off along them (C9)
    # a mark is kept off the loops' lines (C5) ...
    tethers: list[tuple] = []  # ... and the tethers drawn before: (from, to)
    out = []
    na_spots = na_lines = None  # where nucleotides are drawn, the nucleic acid's own boxes, and its strands as drawn
    # among the elements (each a segment array): for ligands held by it alone, built once per page
    # by first holding residue, in the file's order on a tie (as before C33), those held by nucleotides alone last
    axis = [g for g in chosen if g.contacts and _on_axis(layout, g)]  # a rosette's: one cluster at its centre
    if axis:
        out += _axis_cluster(axis, taken, texts)
    for g in sorted(chosen, key=lambda g: (not g.contacts, g.contacts[0] if g.contacts else 0)):
        if any(g is a for a in axis):
            continue
        bound: list = []  # its main contact, then its targets on other chains than its home (the first target's)
        if g.contacts:
            targets, weights, tethered = _local_targets(layout, loops, sses, g, g.metal)
            keep_off, words = taken, texts
            if layout.nucleic is not None:  # kept off the stems, strands and their text too
                if na_spots is None:  # (random-100 held-out 8ZTU: Mg a6201, held by a protein residue, on P stem 188)
                    na_spots = _na_spots(layout, look)
                keep_off, words = taken + na_spots[1], texts + na_spots[2]
            chains = [layout.res_chain[r] for r in tethered]
            others = [t for t, c in zip(targets, chains) if c != chains[0]]
            home = len(targets) - len(others)  # the home chain's targets come first (_local_targets)
            if others:
                bound = [targets[int(np.argmax(weights[:home]))], *others]
        else:  # held by nucleotides alone (C33: 8XZE): tethered to them, kept off stems, letters and names
            if na_spots is None:
                na_spots = _na_spots(layout, look)
            targets, weights, tethered = _na_targets(layout, g, na_spots[0])
            keep_off, words = taken + na_spots[1], texts + na_spots[2]
        if not targets:
            continue
        anchor = np.average(np.asarray(targets), axis=0, weights=weights)
        half = (0.3, 0.3) if g.metal else (0.12 * len(g.name) + 0.2, 0.27)
        na, retarget = None, None
        if not g.contacts:  # held by nucleotides alone: its tethers, as drawn from each spot, off the stems and strands
            if na_lines is None:  # a run's stretches along its own stems' rails are the stems' (na_paths.among)
                at = {key: d for d in layout.stems for key in d.axial}
                strands = [np.asarray(among(p, at), float).reshape(-1, 2) for p in layout.na_paths]
                na_lines = [np.hstack([q[:-1], q[1:]]) for q in strands if len(q) > 1]
            na = (layout.stems, na_lines)
            retarget = lambda c, g=g, where=na_spots[0]: _na_targets(layout, g, where, np.asarray(c, float))[0]  # noqa: E731
        spot = _place_mark(anchor, half, targets, keep_off, lines, words, tethers, bodies, runs, na, retarget, bound)
        if not g.contacts:  # held by nucleotides alone: each tether to the touched nucleotide nearest the mark
            targets = _na_targets(layout, g, na_spots[0], np.asarray(spot, float))[0]
        crowded = overlaps_any((spot[0] - half[0], spot[1] - half[1], spot[0] + half[0], spot[1] + half[1]), keep_off)
        mark = (spot[0] - half[0], spot[1] - half[1], spot[0] + half[0], spot[1] + half[1])
        taken.append(mark)
        texts.append(mark)  # a later ligand's tether keeps off this one's name
        tethers += [(tuple(spot), tuple(t)) for t in targets]
        out.append(
            {
                "ligand": g,
                "centre": (float(spot[0]), float(spot[1])),
                "half": half,
                "targets": targets,
                "residues": tethered,
                "crowded": bool(crowded),  # no clear spot: drawn over an element or text (Drawn.crowded_marks)
            }
        )
    return out


def _draw_ligands(ax, marks: list[dict], lw: float, font: float, layout: Layout, pt_per_unit: float) -> None:
    """Each ligand's mark and name, and its tethers: drawn under the bodies and the text (only what shows of each,
    visible_pieces), so none runs over a letter."""
    bodies = link_bodies(layout, ax)
    keep = link_keep(lw * 0.7, pt_per_unit, ax.figure.dpi)
    halos = text_halos(ax, pt_per_unit, keep)
    for d in layout.stems:  # names' margin, as routes keep it, the tether's own half width beyond
        x0, y0, x1, y1 = stem_keep_box(d, glyph(layout)[0] / _GLYPH_W)
        halos.append((x0 - keep, y0 - keep, x1 + keep, y1 + keep))
    # and the nucleic acid's letters and ends, drawn after the tethers (as the layout measures their ink)
    halos += [(x0 - keep, y0 - keep, x1 + keep, y1 + keep) for x0, y0, x1, y1 in na_text_boxes(layout)]
    for m in marks:
        g, (cx, cy), (hw, hh) = m["ligand"], m["centre"], m["half"]
        key = f"{g.name}:{g.chain}{g.seq}"
        for k, (tx, ty) in enumerate(m["targets"]):
            pieces = visible_pieces([(cx, cy), (tx, ty)], bodies, halos)
            if not pieces:
                continue
            tether = PathPatch(
                _pieces_path(pieces),
                fc="none",
                ec=legible("#6b6b6b"),
                lw=lw * 0.7,
                ls=(0, (1.2, 1.4)),
                capstyle="round",
                zorder=Z["tether"],
            )
            tether.set_gid(f"ligand-tether:{key}:{k}")
            ax.add_patch(tether)
        if g.metal:
            fill = METAL_COLOURS.get(g.symbol, "#7f7f7f")
            marker = Circle((cx, cy), hw, fc=fill, ec=darken(fill, 0.7), lw=lw * 0.7, zorder=Z["mark"])
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
                zorder=Z["mark"],
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
            zorder=Z["mark-label"],
        )
        t.set_gid(f"ligand-label:{key}")


_LINKED_STAR = 0.38  # radius of the mark of a residue fused into the chain
_LINKED_LABEL = 0.6  # its label, as a fraction of the font
_LINKED_LABEL_MIN = 7.0  # but never smaller than this (pt)


def linked_marks(layout: Layout, loops: list[Loop], sses: list[SSE], marks: list[dict]) -> list[dict]:
    """A star on the chain at each peptide-linked non-standard residue (a chromophore: GFP's CRO66) of a drawn
    chain, and its code and number beside it, in the first spot on rings around it clear of elements, ghosts,
    names, other marks and loop paths. Each mark: residue label, star centre, label centre, label half size."""
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
    segments = [(p, q) for lp in loops for p, q in zip(lp.points, lp.points[1:])]
    out = []
    for k in shown:
        at = residue_point(layout, loops, sses, layout.res_chain, k)
        if at is None:
            continue
        label = layout.res_labels[k]
        text = f"{label.name}{label.tag}"
        pt, unit = page_text(layout)  # its ink as drawn (_draw_linked), at least _LINKED_LABEL_MIN pt (C28)
        size = max(_LINKED_LABEL_MIN, pt * _LINKED_LABEL) * unit
        half = count_half(text, size)
        taken.append((at[0] - _LINKED_STAR, at[1] - _LINKED_STAR, at[0] + _LINKED_STAR, at[1] + _LINKED_STAR))
        spot = None
        for ring in np.arange(_LINKED_STAR + 0.1, 12.0, 0.35):
            for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False):
                c = np.asarray(at) + ring * np.array([np.cos(angle), np.sin(angle)])
                box = (c[0] - half[0] - 0.12, c[1] - half[1] - 0.05, c[0] + half[0] + 0.12, c[1] + half[1] + 0.05)
                if not any(
                    box[0] < t[2] and t[0] < box[2] and box[1] < t[3] and t[1] < box[3] for t in taken
                ) and not any(clip_len(p, q, box) > 1e-9 for p, q in segments):
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
            zorder=Z["mark"],
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
                zorder=Z["mark-label"],
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
        halo = Circle((x, y), 0.27, fc=_BACKGROUND.get(), ec=legible("#1d1d1f"), lw=lw * 0.9, zorder=Z["feature-ring"])
        halo.set_gid(f"feature-ring:{n}")
        ax.add_patch(halo)
        dot = Circle((x, y), 0.17, fc=colour, ec="none", zorder=Z["feature"])
        dot.set_gid(f"feature:{n}")
        ax.add_patch(dot)


_HIDE_STEP = 0.02  # page units: a bar or a tether is sampled this finely where a body or a text hides it
_SS_PIECE_MIN = 0.15  # page units: a stretch of a disulfide bar shorter than this between bodies or text is not drawn


def _length_of(pts) -> float:
    return float(np.hypot(*np.diff(np.asarray(pts, float), axis=0).T).sum())


def _dense(points) -> tuple[np.ndarray, np.ndarray]:
    """The polyline's points with others between them, at most _HIDE_STEP apart, and for each the index of the
    segment it lies on (its end point: the segment it ends)."""
    pts = np.asarray(points, float)
    out, seg = [pts[:1]], [np.zeros(1, int)]
    for k, (a, b) in enumerate(zip(pts, pts[1:])):
        n = max(1, int(np.ceil(np.hypot(*(b - a)) / _HIDE_STEP)))
        out.append(a + np.outer(np.arange(1, n + 1) / n, b - a))
        seg.append(np.full(n, k))
    return np.vstack(out), np.concatenate(seg)


def _in_boxes(pts: np.ndarray, boxes) -> np.ndarray:
    """Which of the points lie strictly inside any of the boxes (x0, y0, x1, y1)."""
    if not len(boxes):
        return np.zeros(len(pts), bool)
    b = np.asarray(boxes, float).reshape(-1, 4)
    x, y = pts[:, :1], pts[:, 1:2]
    return ((x > b[:, 0]) & (x < b[:, 2]) & (y > b[:, 1]) & (y < b[:, 3])).any(axis=1)


def link_bodies(layout: Layout, ax=None) -> list[Path]:
    """The outlines (page units) of the bodies a disulfide bar or a ligand tether passes under: each strand's and
    ghost's arrow, each helix's and 3-10 helix's box (in that order: placed elements, then ghosts), and with `ax`
    the helices' and 3-10 helices' drawn patches after them too (a ribbon's end caps reach past its box)."""
    out = []
    for p in layout.placed.values():
        if p.sse.kind == "E":
            out.append(_ring(_arrow(p)))
            continue
        u = np.array([np.cos(p.angle), np.sin(p.angle)]) * p.length / 2
        v = np.array([-np.sin(p.angle), np.cos(p.angle)]) * p.width / 2
        c = np.array([p.cx, p.cy])
        out.append(_ring([c - u - v, c + u - v, c + u + v, c - u + v]))
    out += [_ring(_arrow(g)) for g in layout.ghosts]
    for a in ax.patches if ax is not None else ():
        if (a.get_gid() or "").startswith(("helix:", "helix-back:", "eta:")):
            out.append(a.get_path().transformed(a.get_transform() - ax.transData))
    return out


def _ring(corners) -> Path:
    """A closed outline through every one of `corners`: the first repeated at the end, as Path's closed=True
    takes the last vertex for its CLOSEPOLY code and ignores where it lies."""
    v = np.asarray(corners, float)
    return Path(np.vstack([v, v[:1]]), closed=True)


def text_halos(ax, pt_per_unit: float, keep: float, skip=()) -> list[tuple[float, float, float, float]]:
    """The ink boxes (page units, from text_metrics, as the layout measures text) of the texts drawn on the page so
    far, grown by `keep`: where a disulfide bar or a ligand tether is not drawn, so none runs through or by a
    letter. Texts whose gid starts with one of `skip` are left out."""
    up_down = (text_metrics.ASCENT + text_metrics.DESCENT) / text_metrics.UNITS  # a line box's height, ems
    out = []
    for t in ax.texts:
        words, gid = t.get_text(), t.get_gid() or ""
        if not words.strip() or t.get_transform() != ax.transData or (skip and gid.startswith(tuple(skip))):
            continue
        size = t.get_fontsize() / pt_per_unit  # page units to the em
        ink, start = _ink_line(words, str(t.get_fontweight()) in ("bold", "700", "heavy", "black"))
        x, y = t.get_position()
        ha, va = t.get_horizontalalignment(), t.get_verticalalignment()
        cx = x + {"left": -start[0], "right": start[0]}.get(ha, 0.0) * size
        cy = y + {"bottom": up_down / 2, "top": -up_down / 2, "baseline": -start[1]}.get(va, 0.0) * size
        turn = np.radians(t.get_rotation())
        rot = np.array([[np.cos(turn), -np.sin(turn)], [np.sin(turn), np.cos(turn)]])
        x0, y0, x1, y1 = _reach(ink @ rot.T * size + (cx, cy))
        out.append((x0 - keep, y0 - keep, x1 + keep, y1 + keep))
    return out


def turned_halos(ax, pt_per_unit: float, keep: float) -> list[tuple[float, float, float, float]]:
    """text_halos for the interaction lines, with turned texts placed as matplotlib places them: rotation_mode
    "anchor" turns the alignment offset with the text about its anchor (text_halos leaves the offset unturned, so a
    turned left- or baseline-aligned text, a stem's first-pair numbers, is boxed beside its ink), "default" aligns
    the turned line box by ha and va. Each box holds the turned ink (axis-aligned), grown by `keep`. text_halos
    stays as it is: the disulfide bars and tethers of every default figure are placed by it."""
    up_down = (text_metrics.ASCENT + text_metrics.DESCENT) / text_metrics.UNITS  # a line box's height, ems
    out = []
    for t in ax.texts:
        words = t.get_text()
        if not words.strip() or t.get_transform() != ax.transData:
            continue
        size = t.get_fontsize() / pt_per_unit  # page units to the em
        ink, start = _ink_line(words, str(t.get_fontweight()) in ("bold", "700", "heavy", "black"))
        ha, va = t.get_horizontalalignment(), t.get_verticalalignment()
        turn = np.radians(t.get_rotation())
        rot = np.array([[np.cos(turn), -np.sin(turn)], [np.sin(turn), np.cos(turn)]])
        if t.get_rotation_mode() == "anchor" or not t.get_rotation() % 360:
            shift = np.array(
                [
                    {"left": -start[0], "right": start[0]}.get(ha, 0.0),
                    {"bottom": up_down / 2, "top": -up_down / 2, "baseline": -start[1]}.get(va, 0.0),
                ]
            )
            centre = np.asarray(t.get_position(), float) + rot @ shift * size
        else:  # the turned line box, aligned
            half = np.array([-start[0], up_down / 2])
            corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * half @ rot.T
            lo, hi = corners.min(axis=0), corners.max(axis=0)
            shift = np.array(
                [
                    {"left": -lo[0], "right": -hi[0]}.get(ha, 0.0),
                    {"top": -hi[1]}.get(va, 0.0 if va == "center" else -lo[1]),
                ]
            )
            centre = np.asarray(t.get_position(), float) + shift * size
        x0, y0, x1, y1 = _reach(ink @ rot.T * size + centre)
        out.append((x0 - keep, y0 - keep, x1 + keep, y1 + keep))
    return out


def visible_pieces(points, bodies: list[Path], halos=()) -> list[np.ndarray]:
    """The stretches of polyline `points`, in order along it, outside every body (drawn over it) and every halo
    (around a text): what of a disulfide bar or a ligand tether is drawn, each stretch from where it comes out to
    where it goes under, through the polyline's own corners between (all of it, as given, when nothing hides it). A
    body hides the rest; a text keeps it off its letters."""
    pts = np.asarray(points, float)
    dense, seg = _dense(pts)
    hidden = _in_boxes(dense, halos)
    within = _within_extents(dense, bodies)
    for k in np.flatnonzero(within.any(axis=0)).tolist():
        near = within[:, k] & ~hidden
        if near.any():
            hidden[near] |= bodies[k].contains_points(dense[near])
    pieces, start = [], None
    for k, h in enumerate([*hidden, True]):
        if not h and start is None:
            start = k
        elif h and start is not None:
            if k - start > 1:  # its first and last points, and the corners of the polyline between them
                corners = pts[seg[start] + 1 : seg[k - 1] + 1]
                pieces.append(np.vstack([dense[start], corners, dense[k - 1]]))
            start = None
    return pieces


def _pieces_path(pieces: list[np.ndarray]) -> Path:
    """One path of the visible pieces of a line (each its own stroke)."""
    codes = np.concatenate([[Path.MOVETO] + [Path.LINETO] * (len(p) - 1) for p in pieces])
    return Path(np.vstack(pieces), codes)


def link_keep(lw_pt: float, pt_per_unit: float, dpi: float) -> float:
    """How far (page units) a line of width lw_pt points keeps its centre from a letter's ink: half its width and
    figcheck.LINE_TEXT_CLEARANCE_PX (as figcheck measures a line by a text), with half a pixel to spare."""
    from .figcheck import LINE_TEXT_CLEARANCE_PX

    return lw_pt / 2 / pt_per_unit + (LINE_TEXT_CLEARANCE_PX + 0.5) / (pt_per_unit * dpi / 72)


def disulfide_bars(
    layout: Layout, loops: list[Loop], sses: list[SSE], look: Style, lw: float, pt_per_unit: float, avoid=()
):
    """Each disulfide bar to draw: (i, j, cysteine point a, cysteine point b, the bar's route); a cysteine in a loop
    keeps out of the boxes in `avoid` where its loop leaves room (disulfide_points)."""
    links = layout.links
    if links is None or not look.disulfides:
        return []
    chains = layout.res_chain
    shown = lambda r: not layout.partial or chains[r] not in layout.partial  # noqa: E731 - fragments: no links
    bars = []
    for i, j in links.disulfides:
        if not (shown(i) and shown(j)):
            continue
        a, b = disulfide_points(layout, loops, sses, chains, i, j, avoid)
        if a is not None and b is not None:
            bars.append((i, j, a, b))
    width = lw * 2.4 / pt_per_unit  # bar width in page units
    step = disulfide_step(lw, pt_per_unit)
    lanes = disulfide_lanes([(a, b) for _, _, a, b in bars], width)
    keep_out = [p.rect for p in layout.placed.values()] + [g.rect for g in layout.ghosts]
    keep_out += [box for _, box in label_boxes(layout, sses)]
    out = []
    for (i, j, a, b), lane in zip(bars, lanes):
        side = disulfide_side(a, b, lane, keep_out, width, step) if lane else 1
        out.append((i, j, a, b, disulfide_route(a, b, lane, side, step)))
    return out


def strand_label_spot(p: Placed, label: str, lines, em: float, keep: float) -> tuple[float, float]:
    """Where a strand's letter goes: its centre, else the spot nearest the centre along the strand's shaft whose ink
    box (text_metrics) keeps `keep` page units clear of every line in `lines` (the disulfide bars' routes, which pass
    under the strand and show at its edges); the centre when no spot along the shaft is clear."""
    if not lines:
        return p.cx, p.cy
    ink, _ = _ink_line(label, True)
    x0, y0, x1, y1 = _reach(ink * em)
    axis = np.array([np.cos(p.angle), np.sin(p.angle)])
    room = p.length / 2 - min(HEAD_LEN, p.length * 0.45) - max(abs(y0), abs(y1))  # the shaft, either side
    for k in range(int(max(room, 0.0) / 0.05) + 1):
        for t in (0.0,) if k == 0 else (-0.05 * k, 0.05 * k):
            cx, cy = p.cx + axis[0] * t, p.cy + axis[1] * t
            box = (cx + x0 - keep, cy + y0 - keep, cx + x1 + keep, cy + y1 + keep)
            if not any(clip_len(u, v, box) > 1e-9 for q in lines for u, v in zip(q, q[1:])):
                return float(cx), float(cy)
    return p.cx, p.cy


def _draw_links(
    ax,
    layout: Layout,
    loops: list[Loop],
    sses: list[SSE],
    look: Style,
    lw: float,
    pt_per_unit: float,
    bars=None,
    counts=None,
) -> None:
    """Disulfide bars (`bars`: disulfide_bars' routes, the ones the strand letters slid off; computed here when
    None) and glycans."""
    links = layout.links
    if links is None:
        return
    chains = layout.res_chain
    shown = lambda r: not layout.partial or chains[r] not in layout.partial  # noqa: E731 - fragments: no links
    if look.disulfides:
        keep = link_keep(lw * 2.4, pt_per_unit, ax.figure.dpi)
        halos = text_halos(ax, pt_per_unit, keep)
        bodies = link_bodies(layout, ax)
        if bars is None:
            bars = disulfide_bars(layout, loops, sses, look, lw, pt_per_unit, halos)
        for i, j, a, b, route in bars:
            # drawn under the bodies and the text (only what shows of it); a cysteine in the open keeps its sulfur
            # dot, one under a body gets its end mark on the body's edge, where the bar comes out from under it
            pieces = visible_pieces(route, bodies, halos)
            kept = [q for q in pieces if _length_of(q) >= _SS_PIECE_MIN] or pieces[:1]  # no slivers
            if pieces:
                bar = PathPatch(
                    _pieces_path(kept),
                    fc="none",
                    ec=_SULFUR,
                    lw=lw * 2.4,
                    capstyle="butt",
                    joinstyle="round",
                    zorder=Z["bar"],
                )
                bar.set_gid(f"disulfide:{i}-{j}")
                ax.add_patch(bar)
            hidden = [any(body.contains_point(q) for body in bodies) for q in (a, b)]
            ends = [pieces[0][0], pieces[-1][-1]] if pieces else [a, b]
            for gid, at in (
                (f"ss-dot:{i}-{j}", [q for q, h in zip((a, b), hidden) if not h or not pieces]),
                (f"disulfide-end:{i}-{j}", [e for e, h in zip(ends, hidden) if h and pieces]),
            ):
                if not at:
                    continue
                dots = ax.scatter(
                    [float(q[0]) for q in at],
                    [float(q[1]) for q in at],
                    s=(lw * 5.5) ** 2,
                    c=_SULFUR,
                    edgecolors="#8a6d00",
                    linewidths=lw * 0.5,
                    zorder=Z["bar-end"],
                )
                dots.set_gid(gid)
    if look.glycans:
        x0, y0, x1, y1 = layout.bounds
        centre = np.array([(x0 + x1) / 2, (y0 + y1) / 2])
        taken = [p.rect for p in layout.placed.values()] + [g.rect for g in layout.ghosts]
        taken += [box for _, box in label_boxes(layout, sses)]
        for loop in loops:  # loop lines too, as thin boxes along each segment
            taken += thin_boxes(loop.points)
        size, step = 0.34, 0.42
        stem_halos = text_halos(ax, pt_per_unit, link_keep(lw * 0.8, pt_per_unit, ax.figure.dpi))
        texts = [box for _, box in label_boxes(layout, sses)] + [  # what a stem must not cross (C28)
            centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25)
            for p in layout.placed.values()
            if p.sse.kind == "E"
        ]
        texts += (
            counts
            if counts is not None
            else [t["box"] for t in [*tail_counts(layout, sses, loops), *loop_counts(layout, loops, sses)]]
        )
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
            for clear_stem in (True, False):  # a stem across no text where there is one, else as before
                for out in (np.array([1.0, 0.0]), np.array([-1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.0, -1.0])):
                    for lead in np.arange(0.35, 6.0, 0.1):  # the shortest stem that clears everything
                        boxes = [
                            centred_box(*(at + out * (lead + step * k)), size / 2 + 0.03, size / 2 + 0.03)
                            for k in range(count)
                        ]
                        if not any(overlaps_any(b, taken) for b in boxes):
                            stem_at = at + out * (lead - size / 2)
                            if clear_stem and any(clip_len(at, stem_at, b) > 1e-9 for b in texts):
                                break  # this way the stem crosses a name: the next way
                            cost = lead - 0.15 * float(np.dot(out, away) > 0)  # ties: point outward
                            if best is None or cost < best[0]:
                                best = (cost, out, lead, boxes)
                            break
                if best is not None:
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
            # drawn off the letters it would cross (as a ligand's tether): 7SMR's only way out of its sheet ran down
            # strand G over the letter G
            pieces = visible_pieces([tuple(at), tuple(tip)], [], stem_halos)
            stem = PathPatch(
                _pieces_path(pieces) if pieces else Path(np.empty((0, 2))),
                fc="none",
                ec=legible("#2b2b2b"),
                lw=lw * 0.8,
                zorder=Z["bar-end"],
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
    if domain_coloured(look, layout):
        return _by_domain(sses, layout)
    if look.color_by in ("chain", "domain"):  # by domain with no domain named: by chain

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


def domain_tone(k: int) -> str:
    """The colour of the k-th named domain (Layout.domains order): its panel's tint, its bar segment and, coloured by
    domain, its elements. Eight tones (Cas9 has eight domains), used again from the ninth domain on."""
    return _DOMAIN_TONES[k % len(_DOMAIN_TONES)]


def tone_of(layout: Layout, d: int) -> int:
    """The tone index of the d-th domain (Layout.domains order): copies of one domain share one
    (Layout.domain_tones), else its place."""
    return layout.domain_tones[d] if d < len(layout.domain_tones) else d


_GREY_FILLS = (PALETTES["greys"][2].lower(), PALETTES["greys"][1].lower())  # print's bar: two light greys, alternating
_HATCHES = ("", "////", "\\\\\\\\", "xxx", "....", "||||", "----", "++")  # ... and a hatch per domain


def domain_paint(k: int, look: Style) -> tuple[str, str]:
    """(fill, hatch) of the k-th domain's panel tint and bar segments. The domain tones (domain_tone), unless the
    theme keeps its own palette (domain_colours off, as print and blueprint do, and no color_by=domain): then the
    palette's colours, or for the greys palette two light greys alternating, each domain with a hatch of its own."""
    if look.domain_colours or look.color_by == "domain":
        return domain_tone(k), ""
    if look.palette == "greys":
        return _GREY_FILLS[k % 2], _HATCHES[k % len(_HATCHES)]
    colours = PALETTES[look.palette]
    return colours[k % len(colours)].lower(), ""


def _domain_entries(layout: Layout) -> list[tuple[str, str]]:
    """One legend entry per domain tone: the domain's name, or for copies of one domain in several chains the name
    without its chain (D1 · A1, D1 · A2 -> D1)."""
    names: dict[int, list[str]] = {}
    for d, (name, ids) in enumerate(layout.domains):
        if ids:
            names.setdefault(tone_of(layout, d), []).append(name)
    return [(domain_tone(k), ns[0] if len(ns) == 1 else ns[0].split(" · ")[0]) for k, ns in names.items()]


def domain_coloured(look: Style, layout: Layout | None) -> bool:
    """Are elements coloured by domain: with a domain named and drawn, and color_by=domain, or color_by=chain (the
    default) with domain_colours on. Any other colouring (a theme's or the user's) is kept."""
    if layout is None or not any(ids for _, ids in layout.domains):
        return False
    return look.color_by == "domain" or (look.color_by == "chain" and look.domain_colours)


def _by_domain(sses: list[SSE], layout: Layout):
    """Each element in its domain's tone, an element in no domain grey; a loop residue in the colour of the nearer
    element it joins (a loop is split at its middle, as the domain bar splits it)."""
    tone = {k: domain_tone(tone_of(layout, d)) for d, (_, ids) in enumerate(layout.domains) for k in ids}
    out = {s.id: tone.get(s.id, MATE_GREY) for s in sses}
    runs: dict[str, list[tuple[int, int, str]]] = {}
    for s in sorted(sses, key=lambda s: s.start):
        runs.setdefault(s.chain, []).append((s.start, s.end, out[s.id]))

    def by_domain(r: int, chain: str | None = None) -> str | None:
        mine = runs.get(chain)
        if not mine:
            return None
        for (_, a1, ca), (b0, _, cb) in zip(mine, mine[1:]):
            if r <= a1:
                return ca
            if r < b0:
                return ca if r - a1 <= b0 - r else cb
        return mine[-1][2]

    return out, by_domain


def element_colours(layout: Layout, sses: list[SSE], look: Style):
    """(element id -> fill colour, residue colour function) exactly as the figure draws them: the colouring
    mode (drawn like a reference under an identity colouring: the reference's colours, Layout.like_colour), then
    marked elements, then highlight greys."""
    chains = _chains_in_order(sses)
    colors = chain_colors(chains, look.palette)
    element_colour, residue_colour = _colouring(sses, chains, colors, look, layout)
    if getattr(layout, "like_colour", None):  # drawn like a reference (--like): its mates' colours, by identity
        from .like import IDENTITY_COLOURINGS

        if look.color_by in IDENTITY_COLOURINGS:
            element_colour.update({k: v for k, v in layout.like_colour.items() if k in element_colour})
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
    `model` (the io.read_model result the layout was built from) is needed only for look.hbonds / interactions."""
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


INTERACTION_INK = {  # kind -> (colour, dash) of its lines (plan 13c-4). H-bonds keep look.hbond_color (None here)
    # and their short dashes. The other four: dark enough for a thin line on white (contrast 5.2-9.4 to 1) and
    # chosen by CIEDE2000 under normal vision and simulated protan, deutan and tritan sight (Machado 2009): each
    # stands >= 6.6 from every chain colour of the okabe-ito, tol-bright and tol-muted palettes, the H-bond magenta
    # and the DNA/RNA backbones, and >= 9 from the others, so a line never reads as the element it ends on. The
    # dashes alone tell the five apart in print and greyscale: short dashes, long dashes, dense dots, dash-dot,
    # sparse dots.
    "hbond": (None, (0, (3, 2))),
    "salt": ("#cc2a5c", (0, (7, 2))),  # crimson, long dashes
    "stack": ("#125129", (0, (1, 1.2))),  # dark green, dense dots
    "cation-pi": ("#093fe0", (0, (5, 1.5, 1, 1.5))),  # blue, dash-dot
    "hydrophobic": ("#5d5605", (0, (1, 3.5))),  # dark olive, sparse dots
}
INK_ONLY_PALETTES = ("greys", "blueprint")  # the print and blueprint themes: every kind in hbond_color, by dash only
INTERACTION_NAMES = {  # kind -> its legend's name
    "hbond": "H-bonds",
    "salt": "Salt bridges",
    "stack": "π-stacking",
    "cation-pi": "Cation-π",
    "hydrophobic": "Hydrophobic contacts",
}
INTERACTION_FINDERS = {
    "salt": find_salt_bridges,
    "stack": find_pi_stacking,
    "cation-pi": find_cation_pi,
    "hydrophobic": find_hydrophobic,
}
INTERACTION_GAP = 2.0  # line widths between the parallel lines of two kinds drawn between one residue pair


def interaction_ink(kind: str, look: Style) -> tuple[str, tuple]:
    """(colour, dash) of `kind`'s lines: INTERACTION_INK's, but in look.hbond_color for H-bonds and for every kind
    under the print and blueprint palettes (grey on paper, pale on the dark page), told apart by the dash alone."""
    colour, dash = INTERACTION_INK[kind]
    return (look.hbond_color if colour is None or look.palette in INK_ONLY_PALETTES else colour), dash


class InteractionLine(NamedTuple):
    """One drawn line of a kind other than H-bonds: a residue (or nucleotide) pair with at least one interaction of
    that kind. One line per kind and pair: the pair's salt-bridge atom pairs, or its stacked ring pairs (a Trp, a
    purine), are one line, all of them kept in `found`."""

    kind: str  # interactions.KINDS
    gid: str  # ia:salt:A:19-1:DT7
    ends: tuple  # ((x, y), (x, y)): page points of the first and second residue (the line may lie beside, slots)
    length: float  # page units
    found: tuple  # interactions.Interaction between the two, shortest first
    residues: tuple  # Backbone indices of the protein ends
    labels: tuple = ()  # model.ResLabel of the two ends, in the order of `ends`


class KindOverlay(NamedTuple):
    """The lines of one kind other than H-bonds (as HBondOverlay is the H-bonds')."""

    kind: str
    lines: list  # InteractionLine drawn, shortest first
    total: int  # residue pairs with the interaction and both ends drawn, before the interaction_max cap
    unplaced: int = 0  # residue pairs with the interaction left out because an end is not drawn

    def legend(self) -> str:
        n = len(self.lines)
        text = f"{n:,}" if n == self.total else f"{n:,} of {self.total:,} shown"
        if self.unplaced:
            text += f", {self.unplaced:,} with an end not drawn"
        return f"{INTERACTION_NAMES[self.kind]} ({text})"


class InteractionOverlay(NamedTuple):
    """Every interaction line a figure draws: the H-bonds' (their own overlay, hbond: gids, hbond_max) and those of
    each other kind asked for (ia: gids, interaction_max each), with where two or more kinds share a residue pair."""

    hbonds: HBondOverlay
    kinds: tuple = ()  # KindOverlay of each other kind asked for, in interactions.KINDS order
    slots: dict = {}  # gid -> (place, count): a line's place among the kinds drawn between its pair, when count > 1
    dropped: tuple = ()  # interactions.verify_interactions texts of the found ones it flagged: not drawn

    def lines(self) -> list:
        """Every line drawn: the H-bonds', then each kind's."""
        return [*self.hbonds.lines, *(line for part in self.kinds for line in part.lines)]

    def records(self) -> list[dict]:
        """Every line drawn, as the interactive page's hover and the MCP result give it (interaction_record)."""
        return [interaction_record(line) for line in self.lines()]


INTERACTION_ONE = {  # kind -> what one line of it is called (the hover, the MCP result)
    "hbond": "H-bond",
    "salt": "Salt bridge",
    "stack": "π-stacking",
    "cation-pi": "Cation-π",
    "hydrophobic": "Hydrophobic contact",
}


def residue_name(lab) -> str:
    """A residue or nucleotide as the hover notes name it: A:Gln44, 1:DA4."""
    return f"{lab.chain}:{lab.name.title() if len(lab.name) == 3 else lab.name}{lab.tag}"


def _partner_name(atom: str) -> str:
    """An interaction partner's atom as a reader reads it: NZ, ring:6 -> 6-ring, centre -> guanidinium."""
    if atom.startswith("ring:"):
        return f"{atom[len('ring:') :]}-ring"
    return "guanidinium" if atom == "centre" else atom


def interaction_record(line) -> dict:
    """One drawn line (HBondLine or InteractionLine) as data: `key` (its gid), `kind`, `a` and `b` (its residues,
    as residue_name, in the order of its ends), `residues` (Backbone indices of its protein ends), `het` (its
    nucleotide ends, [chain, number], for a 3D view), `distance` (Å, the shortest of the pair) and `atoms`: one
    {a, b, distance, angle, note} per interaction behind it (atoms in the order of a and b; angle in degrees or
    None; note "" or a caveat), shortest first. `text`: the hover note, "Salt bridge A:Arg19 – 1:DT7: NH1 – OP1
    3.1 Å; ...", H-bond lines as before ("H-bond A:Gln44 – 1:DA4: NE2 – N7 3.3 Å")."""
    first, second = line.labels
    kind = getattr(line, "kind", "hbond")
    atoms = []
    if kind == "hbond":
        for x in line.bonds:
            pair = (x.donor.name, x.acceptor.name) if x.donor.res == first else (x.acceptor.name, x.donor.name)
            atoms.append({"a": pair[0], "b": pair[1], "distance": x.distance, "angle": None, "note": ""})
    else:
        for x in line.found:
            pair = (x.a[2], x.b[2]) if (x.a[0], x.a[1]) == (first.chain, first.tag) else (x.b[2], x.a[2])
            atoms.append({"a": pair[0], "b": pair[1], "distance": x.distance, "angle": x.angle, "note": x.note})

    def said(x) -> str:
        more = [f"{x['angle']:.0f}°"] if x["angle"] is not None else []
        more += [x["note"]] if x["note"] else []
        extra = f" ({', '.join(more)})" if more else ""
        return f"{_partner_name(x['a'])} – {_partner_name(x['b'])} {x['distance']:.1f} Å{extra}"

    return {
        "key": line.gid,
        "kind": kind,
        "a": residue_name(first),
        "b": residue_name(second),
        "residues": list(line.residues),
        "het": [[lab.chain, lab.seq] for lab in line.labels[len(line.residues) :]],  # protein ends sort first
        "distance": min(x["distance"] for x in atoms),
        "atoms": atoms,
        "text": f"{INTERACTION_ONE[kind]} {residue_name(first)} – {residue_name(second)}: "
        + "; ".join(said(x) for x in atoms),
    }


def _locator(layout: Layout, loops: list[Loop], sses: list[SSE], look: Style):
    """locate(ResLabel) -> (sort key, page point, Backbone index or None, ResLabel), or None when that residue (or
    nucleotide) is not drawn: a protein residue at render.residue_point (not in a chain left partial or out of the
    focus), a nucleotide at nucleotide_points. Also the drawn residues by (chain, tag)."""
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

    labels = {(l.chain, l.tag): l for l in layout.res_labels}
    if layout.nucleic is not None:
        labels.update({(lab.chain, lab.tag): lab for strand in layout.nucleic.strands for lab in strand.labels})
    return locate, labels


def _as_interaction(bond) -> Interaction:
    """An hbonds.HBond as the re-check reads it."""
    d, a = bond.donor, bond.acceptor
    return Interaction("hbond", (d.res.chain, d.res.tag, d.name), (a.res.chain, a.res.tag, a.name), bond.distance, None)


def interaction_overlay(
    layout: Layout, loops: list[Loop], sses: list[SSE], look: Style, model=None, kinds: tuple[str, ...] | None = None
) -> InteractionOverlay:
    """The interaction lines the look asks for (style.interactions_of: `interactions`, else the older `hbonds`;
    `kinds` narrows them), from the model the layout was built from (io.read_model's result): one line per kind and
    residue pair (protein before nucleotide, else by tag), between their points on the page. Interactions with an
    end that is not drawn (cropped nucleic acid, chains outside the focus, residues past a chain's elements) are left
    out and counted. Every interaction found is re-checked first by interactions.verify_interactions, its own
    reading of the coordinates; one it flags is not drawn and is named in `dropped` (the figure's notes say so).
    At most `look.hbond_max` H-bond lines and `look.interaction_max` lines of each other kind, the shortest on the
    page kept (not the shortest in Å). H-bond lines keep their hbond: gids and draw as before plan 13c; the others'
    gids are ia:<kind>:<a>-<b>."""
    asked = interactions_of(look, layout.nucleic is not None)  # interactions=..., else the older hbonds=...
    wanted = [k for k in asked.kinds if kinds is None or k in kinds]
    if not wanted:
        return InteractionOverlay(HBondOverlay([], 0))
    sel = asked.between
    if model is None:
        given = f"interactions={look.interactions}" if look.interactions != "off" else f"hbonds={look.hbonds}"
        raise ValueError(f"{given} needs the structure: pass model= (the io.read_model result)")
    locate, labels = _locator(layout, loops, sses, look)
    bonds = find_hbonds(model, sel) if "hbond" in wanted else []
    found = {k: INTERACTION_FINDERS[k](model, sel) for k in wanted if k != "hbond"}
    checked = [*(_as_interaction(b) for b in bonds), *(x for k in found for x in found[k])]
    failing = failing_lines(model, checked)
    if failing:  # the re-check's flagged ones are not drawn
        flagged = {n for n, _ in failing}
        bonds = [b for n, b in enumerate(bonds) if n not in flagged]  # checked starts with the H-bonds
        found = {k: [] for k in found}
        for n, x in enumerate(checked):
            if n not in flagged and x.kind != "hbond":
                found[x.kind].append(x)
    hb = _hbond_lines(bonds, locate, look) if "hbond" in wanted else HBondOverlay([], 0)
    parts = tuple(_kind_lines(k, found[k], locate, labels, look) for k in wanted if k != "hbond")
    pairs: dict[str, list[str]] = {}  # residue pair -> the gids drawn between it, in KINDS order
    for line in [*hb.lines, *(line for part in parts for line in part.lines)]:
        pairs.setdefault(line.gid.split(":", 2 if line.gid.startswith("ia:") else 1)[-1], []).append(line.gid)
    slots = {gid: (k, len(gids)) for gids in pairs.values() if len(gids) > 1 for k, gid in enumerate(gids)}
    return InteractionOverlay(hb, parts, slots, tuple(text for _, text in failing))


def _hbond_lines(bonds, locate, look: Style) -> HBondOverlay:
    """The H-bond lines of `bonds` (hbonds.HBond): one per residue pair, gid hbond:<a>-<b>, at most look.hbond_max."""
    pairs: dict[str, list] = {}
    unplaced = set()
    for bond in bonds:
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


def _kind_lines(kind: str, found: list, locate, labels: dict, look: Style) -> KindOverlay:
    """The lines of one kind other than H-bonds from its interactions.Interaction list: one per residue pair, gid
    ia:<kind>:<a>-<b>, at most look.interaction_max, the shortest on the page kept."""
    pairs: dict[str, list] = {}
    unplaced = set()
    for x in found:
        ends = [labels.get((p[0], p[1])) for p in (x.a, x.b)]
        a, b = (locate(lab) if lab is not None else None for lab in ends)
        if a is None or b is None:
            unplaced.add(frozenset((x.a[:2], x.b[:2])))
            continue
        a, b = sorted((a, b), key=lambda e: e[0])
        gid = f"ia:{kind}:{a[0][1]}:{a[0][2]}{a[0][3]}-{b[0][1]}:{b[0][2]}{b[0][3]}"
        pairs.setdefault(gid, [a, b, []])[2].append(x)
    lines = []
    for gid, (a, b, xs) in pairs.items():
        length = float(np.hypot(a[1][0] - b[1][0], a[1][1] - b[1][1]))
        if length > 0:
            residues = tuple(e[2] for e in (a, b) if e[2] is not None)
            xs = tuple(sorted(xs, key=lambda x: (x.distance, sort_key(x))))
            lines.append(InteractionLine(kind, gid, (a[1], b[1]), length, xs, residues, (a[3], b[3])))
    lines.sort(key=lambda l: (l.length, l.gid))
    return KindOverlay(kind, lines[: look.interaction_max], len(lines), len(unplaced))


def hbond_overlay(layout: Layout, loops: list[Loop], sses: list[SSE], look: Style, model=None) -> HBondOverlay:
    """The hydrogen-bond lines the look asks for (style.interactions_of: `interactions` with hbond among its kinds,
    else the older `hbonds`): interaction_overlay's H-bond part. One per residue pair (protein before nucleotide,
    else by tag), between their points on the page; bonds with an end that is not drawn are left out. At most
    `look.hbond_max` lines, the shortest on the page kept (not the shortest in Å). `model` is the io.read_model
    result the layout was built from."""
    return interaction_overlay(layout, loops, sses, look, model, kinds=("hbond",)).hbonds


def hbond_alpha(length: float, diagonal: float) -> float:
    """Opacity of a hydrogen-bond line: 1 for short lines, falling linearly to 0.4 at half the page diagonal."""
    t = (length / diagonal - 0.15) / 0.35 if diagonal > 0 else 0.0
    return 1.0 - 0.6 * min(max(t, 0.0), 1.0)


def _beside(ends: tuple, place: int, count: int, gap: float) -> tuple:
    """The two ends moved square to the line by `gap` page units per place, the `count` places centred on it."""
    (x0, y0), (x1, y1) = ends
    length = float(np.hypot(x1 - x0, y1 - y0))
    s = (place - (count - 1) / 2) * gap / length
    nx, ny = -(y1 - y0) * s, (x1 - x0) * s
    return ((x0 + nx, y0 + ny), (x1 + nx, y1 + ny))


def _draw_lines(c: _Canvas, kind: str, lines: list, slots: dict, halos=None) -> list[str]:
    """One kind's lines, dashed in its ink (interaction_ink); long ones fade (hbond_alpha). A line sharing its
    residue pair with another kind's lies beside it, INTERACTION_GAP line widths apart, the pair's lines centred
    on the residues. With `halos` (turned_halos: the letters on the page, grown by the line's keep), a line is
    broken where it would run under a letter (visible_pieces), so no text is crossed; H-bond lines are drawn whole,
    as before plan 13c (halos None). Returns the gids of lines wholly under text, not drawn."""
    x0, y0, x1, y1 = c.content
    diag = float(np.hypot(x1 - x0, y1 - y0))
    colour, dash = interaction_ink(kind, c.look)
    lw = c.loop_lw * 0.85
    hidden = []
    for line in lines:
        ends = line.ends
        if line.gid in slots:
            ends = _beside(ends, *slots[line.gid], INTERACTION_GAP * lw / c.pt_per_unit)
        xy = np.c_[ends[0], ends[1]]
        if halos is not None:
            pieces = visible_pieces(np.asarray(ends, float), [], halos)
            if not pieces:
                hidden.append(line.gid)
                continue
            gap = np.full((1, 2), np.nan)
            xy = np.vstack([np.vstack([p, gap]) for p in pieces][:-1] + [pieces[-1]]).T
        (artist,) = c.ax.plot(
            *xy,
            color=colour,
            lw=lw,
            ls=dash,
            alpha=hbond_alpha(line.length, diag),
            solid_capstyle="butt",
            zorder=Z["hbond"],
        )
        artist.set_gid(line.gid)
    return hidden


def _draw_hbonds(c: _Canvas, overlay: HBondOverlay, slots: dict | None = None) -> None:
    """Thin dashed lines between the residues of each hydrogen-bonded pair: over the panels and the nucleic-acid
    rungs, under the backbone rails, the contact beads, the elements and the letters (gid hbond:<a>-<b>). Long
    lines fade (alpha 1 up to 15% of the content's diagonal, 0.4 from 50% on), so a crowd of them stays legible."""
    _draw_lines(c, "hbond", overlay.lines, slots or {})


def _draw_interactions(c: _Canvas, overlay: InteractionOverlay) -> list[str]:
    """The lines of every kind but H-bonds (gid ia:<kind>:<a>-<b>), in interactions.KINDS order: drawn after every
    text in the figure (but the legend's and the title), so each keeps off the letters (turned_halos, link_keep), at
    the H-bonds' depth (Z["hbond"]: under the elements, beads and letters). Returns the gids of those wholly under
    text, not drawn (the figure's notes name them)."""
    if not overlay.kinds:
        return []
    halos = turned_halos(c.ax, c.pt_per_unit, link_keep(c.loop_lw * 0.85, c.pt_per_unit, c.fig.dpi))
    return [gid for part in overlay.kinds for gid in _draw_lines(c, part.kind, part.lines, overlay.slots, halos)]


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
    bar_h: float = 0.0  # page units each domain bar takes, under the content (0: no bar)


class _Legend(NamedTuple):
    entries: list
    entry: float  # width of one legend entry (page units)
    per_row: int
    height: float
    row: float = _LEGEND_ROW  # the step between rows (page units): _LEGEND_ROW, grown with the text (C4)


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
    hb: InteractionOverlay | None = None,
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
    elif domain_coloured(look, layout):  # one entry per domain, in place of the chains; the greys kept
        grey = [e for e in entries if e[0] == MATE_GREY]
        named = {k for _, ids in layout.domains for k in ids}
        dim = set(chains) - layout.focus if layout.focus is not None else set()
        if any(k not in named and p.sse.chain not in dim for k, p in layout.placed.items()):
            grey.append((MATE_GREY, "no named domain"))
        entries = _domain_entries(layout) + grey
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
    if hb is not None and hb.hbonds.lines and look.legend:
        entries.append((look.hbond_color, hb.hbonds.legend()))
    for part in hb.kinds if hb is not None and look.legend else ():  # each other kind drawn: its dash, as drawn
        if part.lines:
            entries.append((("ia", part.kind), part.legend()))
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
    texts=(),
    at=None,
) -> tuple[_Canvas, _Legend]:
    """The figure, its axes and every size the steps share; also lays the legend out, as it sets the page width.
    `at` (a panel of a page of several figures, compare.py): called with the page's width and height (page units)
    and its legend entries, it gives (inches per page unit, the figure, the axes) to draw in; the legend is then
    left to the page (no entries drawn here, no room kept for them)."""
    x0, y0, x1, y1 = layout.bounds
    for bx0, by0, bx1, by1 in texts:  # texts placed outside the elements' bounds (tail counts)
        x0, y0, x1, y1 = min(x0, bx0), min(y0, by0), max(x1, bx1), max(y1, by1)
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
        return _legend_entry(entries, x1 - x0 + 2 * pad, look.font_scale)

    entry = legend_entry(x1)
    if x0 + entry - 0.9 > x1:  # never narrower than one legend entry; a wider page may shrink, so measure again
        entry = legend_entry(x0 + entry - 0.9)
        x1 = max(x1, x0 + entry - 0.9)
    per_row = max(1, int((x1 - x0 + 0.9) // entry))
    if not look.legend:
        entries = []
    shown = entries if at is None else []  # a panel's entries go to its page's one legend
    em = _legend_em(x1 - x0 + 2 * pad, look.font_scale)
    row = _LEGEND_ROW * glyph_scale(em)  # a legend row, sized from the text drawn in it (C4)
    legend_h = 0.8 + row * -(-len(shown) // per_row) if shown else 0.3
    if layout.context and look.legend and at is None:
        legend_h += 1.3
    bars = domain_bars(layout) if domain_bar_shown(look, layout) else []
    bar_h = max(_DOMAIN_BAR_H, 2.4 * em + _BAR_TICK) if bars else 0.0
    xmin, xmax, ymin, ymax = x0 - pad, x1 + pad, y0 - pad - len(bars) * bar_h - legend_h, y1 + pad + top
    w_u, h_u = xmax - xmin, ymax - ymin
    if at is None:
        scale = page_scale(w_u, h_u, look.font_scale)  # text at MIN_TEXT_PT or larger: past MAX_WIDTH_IN it grows (C4)
        fig = Figure(figsize=(w_u * scale, h_u * scale), facecolor=look.background)
        ax = fig.add_axes([0, 0, 1, 1])
    else:
        scale, fig, ax = at(w_u, h_u, list(entries))
    fig.foldmap_page_fit = page_fit(w_u * scale, h_u * scale, _font_pt(scale))  # the notes' page line, if any
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    pt_per_unit = scale * 72
    font = _font_pt(scale) * look.font_scale
    lw = line_pt(scale)
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
        bar_h,
    )
    return canvas, _Legend(shown, entry, per_row, legend_h, row)


def _legend_entry(entries: list, width: float, font_scale: float) -> float:
    """One legend entry's width in page units (swatch, text, gap) for the widest of `entries`, on a page `width`
    page units wide (any page's scale: the text's em is the same page units on every page)."""
    font_u = _legend_em(width, font_scale)  # text drawn at that size, in page units
    # the unhinted advance width, which is what is drawn (ink extents differ between matplotlib versions);
    # measured at a size matplotlib does not clamp (it floors font sizes at 1 pt) and scaled to font_u
    widest = max(
        (text_to_path.get_text_width_height_descent(text, _MEASURE, ismath=False)[0] for _, text in entries),
        default=0.0,
    )
    return 0.85 + widest * font_u / _MEASURE.get_size_in_points() + 0.9


def _legend_em(width: float, font_scale: float) -> float:
    """The text's em, page units, on a page `width` page units wide: the same at any scale (not the figure's own
    scale, as its height is not known yet, and need not be)."""
    scale = page_scale(width, font_scale=font_scale)
    return _font_pt(scale) * font_scale / (scale * 72)


class _Measured(Exception):
    """page_size's way out of draw once the page is known (nothing is drawn)."""


def page_size(layout: Layout, loops: list[Loop], sses: list[SSE], title: str | None = None, **kw):
    """(width, height, legend entries) of the page draw(layout, loops, sses, title, **kw) would make, in page units,
    with no legend drawn on it (draw's `at`): for a page of several figures at one scale with one legend
    (compare.py). None when there is nothing to draw."""

    def measure(width: float, height: float, entries: list):
        raise _Measured(width, height, entries)

    try:
        draw(layout, loops, sses, title, at=measure, **kw)
    except _Measured as m:
        return m.args
    return None


def page_legend(entries: list, width: float, look: Style) -> _Legend:
    """The legend of a page `width` page units wide (a pad of 1 at each side) holding `entries`, laid out as a
    figure's own (_canvas): for a page of several figures with one legend (compare.py). Its height is 0 with no
    entries."""
    entry = _legend_entry(entries, width, look.font_scale)
    per_row = max(1, int((width - 2 + 0.9) // entry))
    row = _LEGEND_ROW * glyph_scale(_legend_em(width, look.font_scale))
    height = 0.8 + row * -(-len(entries) // per_row) if entries else 0.0
    return _Legend(list(entries), entry, per_row, height, row)


def draw_page_legend(fig: Figure, ax: Axes, legend: _Legend, bounds: tuple, look: Style, scale: float) -> None:
    """Draw `legend` (page_legend) on `ax`, whose page units are `bounds` (xmin, xmax, ymin, ymax: the legend's rows
    are laid from ymin up, as a figure's are), at `scale` inches per page unit."""
    lw, font = line_pt(scale), _font_pt(scale) * look.font_scale
    pt = scale * 72
    c = _Canvas(
        fig, ax, look, {}, {}, lambda *_: None, font, lw, lw * look.loop_width, pt, 1.0, bounds, bounds, ({}, {})
    )
    with page(look):
        _draw_legend(c, SimpleNamespace(context=None), legend)


_OVER_Z = Z["over"]  # a helix over a sheet (Placed.over): above the strands, below their letters; its name at 4.5
_OVER_ALPHA = 0.6  # its fill and edge: the strands beneath show through
_OVER_TONE = 0.7  # coloured as its strands: darkened this much, to stand apart from them ...
_OVER_SAME_ALPHA = 0.8  # ... its fill at this alpha, its edge (darker still) opaque
_OVER_LOOP_Z = Z["over-loop"]  # the loops of a helix over a sheet: above the strands, below the helix


_COUNT_GAP = 0.15  # page units from a loop's line to its residue count's ink
_COUNT_RUNS = 4  # the loop's longest straight runs tried first for its count
_COUNT_AT = (0.5, 1 / 3, 2 / 3)  # where along a run its count is tried: its middle and thirds ...
_COUNT_AT_ENDS = (1 / 6, 5 / 6)  # ... then nearer its ends


def _count_spots(pts: np.ndarray, half) -> list[tuple[int, tuple[float, float, float, float]]]:
    """A loop count's candidate places, in order, each with the run it sits beside (index into pts): beside the
    _COUNT_RUNS longest runs at _COUNT_AT, then beside the other runs at _COUNT_AT, then beside every run at
    _COUNT_AT_ENDS; each on either side of its run (right then left of an upright run, above then below a level
    one), its ink _COUNT_GAP from the line. Longer runs first, the earlier on a tie; centres snapped (layout._snap)."""
    runs = sorted(range(len(pts) - 1), key=lambda j: (-np.hypot(*(pts[j + 1] - pts[j])), j))
    out = []
    for group, along in ((runs[:_COUNT_RUNS], _COUNT_AT), (runs[_COUNT_RUNS:], _COUNT_AT), (runs, _COUNT_AT_ENDS)):
        for j in group:
            p, q = pts[j], pts[j + 1]
            d = q - p
            vertical = abs(d[1]) > abs(d[0])
            sides = ((1.0, 0.0), (-1.0, 0.0)) if vertical else ((0.0, 1.0), (0.0, -1.0))
            for t in along:
                for n in map(np.array, sides):
                    reach = _COUNT_GAP + abs(n[0]) * half[0] + abs(n[1]) * half[1]
                    c = _snap(p + d * t + n * reach)
                    out.append((j, centred_box(float(c[0]), float(c[1]), *half)))
    return out


def _foot(p, q, c) -> tuple[float, float]:
    """The point of segment p-q nearest c."""
    p, q, c = (np.asarray(v, float) for v in (p, q, c))
    d = q - p
    t = float(np.clip(np.dot(c - p, d) / max(float(np.dot(d, d)), 1e-12), 0.0, 1.0))
    return float(p[0] + d[0] * t), float(p[1] + d[1] * t)


def loop_counts(layout: Layout, loops: list[Loop], sses: list[SSE], tails=None, hidden=None, drawn=None) -> list[dict]:
    """The residue count ("NN aa") of each loop of at least LONG_LOOP residues: name (a>b), text, centre and ink box
    (page units, sized from the drawn text as tail counts are). Its place is the first of its spots (_count_spots:
    beside its runs, longest first) where its ink is clear (layout.clear_inks, as figcheck's A6 reads it) of
    elements, text (label_boxes, tail counts), the counts placed before, every other loop's line as drawn (`drawn`:
    drawn_lines, else the routes), every nucleic-acid strand's and its own route's but the run it sits beside, on
    the page, with no other line between it and its run (layout.spans_crossed: it reads as its loop's). A count
    with no clear spot is left out, never drawn over ink: {"kind": "loop", "name", "text"} is appended to `hidden`
    (a list), which the figure's provenance lists (counts_hidden). Placed after routing: no keep-out."""
    em = _ink_em(layout, layout.font_scale) or text_em(layout.width or 1.0) * layout.font_scale
    by_id = {s.id: s for s in sses}
    taken = [p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems]]
    taken += [b for _, b in label_boxes(layout, sses)] + [
        t["box"] for t in (tails if tails is not None else tail_counts(layout, sses, loops))
    ]
    lines = [segments_of([lp.points]) for lp in loops]
    seen = lines if drawn is None else [segments_of([pts]) for pts in drawn]  # other loops: as drawn
    strands = segments_of(na_strand_lines(layout))  # nucleic-acid strands' and ties' lines
    page = count_page(layout, np.concatenate([strands, *seen]))
    out = []
    for i, loop in enumerate(loops):
        a, b = by_id.get(loop.a_id), by_id.get(loop.b_id)
        if a is None or b is None or b.start - a.end - 1 < LONG_LOOP or len(loop.points) < 2:
            continue
        text = f"{b.start - a.end - 1} aa"
        size = TAIL_COUNT_SIZE * em
        spots = _count_spots(np.asarray(loop.points, float), count_half(text, size))
        ink = count_half(text, size, 0.0)
        inks = [centred_box((x0 + x1) / 2, (y0 + y1) / 2, *ink) for _, (x0, y0, x1, y1) in spots]
        others = np.concatenate([strands, *(seg for k, seg in enumerate(seen) if k != i)])
        ok = clear_inks(inks, taken, others, layout, page)
        pts = np.asarray(loop.points, float)
        centres = [((x0 + x1) / 2, (y0 + y1) / 2) for _, (x0, y0, x1, y1) in spots]
        spans = [(*c, *_foot(pts[j], pts[j + 1], c)) for (j, _), c in zip(spots, centres)]
        ok &= ~spans_crossed(spans, others)  # no other line between it and its run: it reads as its loop's
        own = lines_meeting(inks, lines[i], line_band(layout))  # its own segments, but the run it sits beside
        ok &= ~np.array([own[np.arange(len(own)) != j, k].any() for k, (j, _) in enumerate(spots)], bool)
        name = f"{loop.a_id}>{loop.b_id}"
        box = next((bx for (_, bx), clear in zip(spots, ok.tolist()) if clear), None)
        if box is None:
            if hidden is not None:
                hidden.append({"kind": "loop", "name": name, "text": text})
            continue
        taken.append(box)
        out.append({"name": name, "text": text, "centre": ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2), "box": box})
    return out


def _loop_boxes(layout: Layout) -> np.ndarray:
    """The element boxes a curved loop's corners keep out of (_smooth_path), a hair inside each."""
    boxes = np.array([p.rect for p in [*layout.placed.values(), *layout.ghosts]]).reshape(-1, 4)  # none: DNA only
    return boxes + [1e-6, 1e-6, -1e-6, -1e-6]


def _loop_path(loop: Loop, look: Style, boxes: np.ndarray) -> Path:
    """A loop's line as drawn: a fallback curve as it is, else the route with smooth (Style.loops curved) or rounded
    corners."""
    if loop.fallback:  # already a smooth curve
        return _polyline(loop.points)
    if look.loops == "curved":
        return _smooth_path(loop.points, boxes)
    return _rounded_path(loop.points)


def drawn_lines(layout: Layout, loops: list[Loop], look: Style) -> list[np.ndarray]:
    """Each loop's line as drawn (_loop_path), its curves flattened as figcheck reads them (Path.to_polygons): the
    lines the counts placed after routing keep from (loop_counts, layout.tail_counts)."""
    boxes = _loop_boxes(layout)
    out = []
    for loop in loops:
        pieces = _loop_path(loop, look, boxes).to_polygons(closed_only=False) if len(loop.points) > 1 else []
        out.append(np.concatenate(pieces) if pieces else np.asarray(loop.points, float).reshape(-1, 2))
    return out


def _draw_loops(c: _Canvas, layout: Layout, loops: list[Loop], sses: list[SSE], lcounts=None) -> None:
    """Each loop: a casing in the background colour, its coloured line, arrows, residue runs and length labels."""
    boxes = _loop_boxes(layout)
    sse_by_id = {s.id: s for s in sses}
    over = c.backbone[1]
    counts = {count["name"]: count for count in (lcounts if lcounts is not None else loop_counts(layout, loops, sses))}
    lifted = {k for k, p in layout.placed.items() if p.over}  # helices over a sheet: their loops lie on it too
    for i, loop in enumerate(loops):
        z = (_OVER_LOOP_Z if lifted & {loop.a_id, loop.b_id} else Z["loop"]) + i * 1e-3
        path = _loop_path(loop, c.look, boxes)
        name = f"{loop.a_id}>{loop.b_id}"
        casing = PathPatch(
            path,
            fc="none",
            ec=c.look.background,
            lw=c.loop_lw * 3.0,
            capstyle="round",
            joinstyle="round",
            zorder=z,
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
        loop_c = legible(loop_c)
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
                    zorder=z + 7e-4,
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
            zorder=z + 5e-4,
        )
        line.set_gid(f"loop:{name}")
        c.ax.add_patch(line)
        count = counts.get(name)
        if count is not None:  # the residue count, beside one of the loop's long straight runs
            t = c.ax.text(
                *count["centre"],
                count["text"],
                ha="center",
                va="center",
                fontsize=c.font * TAIL_COUNT_SIZE,
                color=legible("#5d6877"),
                zorder=Z["label"],
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
                    zorder=z + 6e-4,
                )
                run.set_gid(f"loop-seg:{name}:{k}")
                c.ax.add_patch(run)


def _draw_axis_symbol(c: _Canvas, layout: Layout) -> None:
    """A rosette's n-fold axis (Layout.rosette, Style.cyclic): a filled regular n-gon at the page centre, one corner
    up, as crystallographers mark a rotation axis seen end on."""
    n = getattr(layout, "rosette", 0)
    if not n:
        return
    turns = np.pi / 2 + 2 * np.pi * np.arange(n) / n
    corners = np.stack([AXIS_SYMBOL_R * np.cos(turns), AXIS_SYMBOL_R * np.sin(turns)], axis=1)
    mark = Polygon(corners, closed=True, fc=c.look.ink, ec="none", zorder=Z["mark"])
    mark.set_gid(f"axis:{n}")
    c.ax.add_patch(mark)


def _draw_membrane(c: _Canvas, layout: Layout) -> None:
    """The lipid bilayer: over the panels and plates, behind the loops and elements."""
    x0, y0, x1, y1 = c.content
    if layout.membrane is None:
        return
    # a translucent band with head groups on both faces; one per row of a wrapped stack
    for k, (b0, b1) in enumerate(layout.membrane.get("bands", [layout.membrane["y"]])):
        row = f":{k}" if k else ""
        left, right = x0 - c.pad * 0.7, x1 + c.pad * 0.7
        band = Rectangle(
            (left, b0), right - left, b1 - b0, fc=_LIPID, ec="none", alpha=_LIPID_ALPHA, zorder=Z["membrane"]
        )
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
                zorder=Z["membrane-heads"],
            )
            heads.set_gid(f"membrane-heads:{face}{row}")
        below = layout.membrane["inside"] == "below"
        sides = (("out" if below else "in", b1 + 0.35, "bottom"), ("in" if below else "out", b0 - 0.35, "top"))
        for text, y, va in sides if layout.membrane.get("sides", True) else ():
            t = c.ax.text(left + 0.1, y, text, ha="left", va=va, fontsize=c.font * 0.9, style="italic", color="#7a6a45")
            t.set_gid(f"membrane-side:{text}{row}")


_NAME_PLATE_ALPHA = 0.85  # a domain name's plate: the page colour at this alpha over the panel's edge
_NAME_PLATE_PAD = 0.06  # page units the plate reaches beyond the name's box


def domain_outline(fields, name: str) -> Path | None:
    """A domain's panel outline (layout.domain_fields): the zero contour of its field, as one compound path
    (pieces and holes), None when the panel is empty."""
    from contourpy import FillType, contour_generator

    gen = contour_generator(fields.xs, fields.ys, fields.values[name], fill_type=FillType.OuterCode)
    points, codes = gen.filled(0.0, 2.0)
    if not points:
        return None
    return Path(np.vstack(points), np.concatenate(codes))


_COPY_DASHES = ("-", (0, (5, 2)), (0, (1.5, 1.5)), (0, (5, 2, 1.5, 2)))  # copies of one domain: one edge each


def _draw_domains(c: _Canvas, layout: Layout) -> dict[str, Path]:
    """Named domains: a soft tinted panel each, following its elements (layout.domain_fields), its name on a plate
    over the panel's top edge (or above it, with a leader). Copies of one domain share its tone (tone_of) but each
    has an edge of its own, deeper and dashed differently (_COPY_DASHES), so interlocked copies read apart. The
    outlines drawn, by domain."""
    fields = domain_fields(layout)
    tones = [tone_of(layout, k) for k in range(len(layout.domains))]
    shared = Counter(t for t, (_, ids) in zip(tones, layout.domains) if ids)
    seen: Counter = Counter()
    outlines: dict[str, Path] = {}
    for k, (name, ids) in enumerate(layout.domains):
        if not ids or name not in fields.values:
            continue
        outline = domain_outline(fields, name)
        tone = domain_paint(tones[k], c.look)[0]
        fill = _mix(tone, c.look.background, 0.86)
        copy = seen[tones[k]]
        seen[tones[k]] += 1
        if outline is not None:
            outlines[name] = outline
            many = shared[tones[k]] > 1
            panel = PathPatch(
                outline,
                fc=fill,
                ec=_mix(tone, c.look.background, (0.15 + 0.15 * (copy % 3)) if many else 0.55),
                lw=c.lw * (0.8 if many else 0.6),
                ls=_COPY_DASHES[copy % len(_COPY_DASHES)] if many else "-",
                zorder=Z["domain-panel"],
            )
            panel.set_gid(f"domain:{name}")
            c.ax.add_patch(panel)
        if name not in layout.domain_names:
            continue
        x0, y0, x1, y1 = layout.domain_names[name]
        plate = FancyBboxPatch(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            boxstyle=f"round,pad={_NAME_PLATE_PAD},rounding_size={(y1 - y0) / 2 + _NAME_PLATE_PAD}",
            fc=to_rgba(c.look.background, _NAME_PLATE_ALPHA),
            ec="none",
            zorder=Z["domain-name-plate"],
        )
        plate.set_gid(f"domain-name-plate:{name}")
        c.ax.add_patch(plate)
        if name in layout.domain_leaders:
            (lx, ly0), (_, ly1) = layout.domain_leaders[name]
            (leader,) = c.ax.plot(
                [lx, lx], [ly0, ly1], color=_mix(tone, c.look.background, 0.3), lw=c.lw * 0.8, zorder=Z["domain-leader"]
            )
            leader.set_gid(f"domain-leader:{name}")
        t = c.ax.text(
            (x0 + x1) / 2,
            (y0 + y1) / 2,
            name,
            ha="center",
            va="center",
            fontsize=c.font * DOMAIN_NAME_SIZE,
            fontweight="bold",
            color=domain_name_colour(tone, fill),
            zorder=Z["domain-name"],
        )
        t.set_gid(f"domain-label:{name}")
    return outlines


class BarSegment(NamedTuple):
    """One stretch of a chain's domain bar: residues lo..hi (a boundary may fall half way between two residues)."""

    name: str
    gid: str  # domain-bar:<chain>:<name>, then #2, #3 for the later stretches of a domain in pieces
    lo: float
    hi: float
    tone: str
    hatch: str = ""  # with a theme's greys: each domain's own hatch (domain_paint)


def domain_bar_shown(look: Style, layout: Layout) -> bool:
    """domain_bar on: whenever a named domain is drawn; auto: when two or more are; off: never."""
    drawn = sum(bool(ids) for _, ids in layout.domains)
    if look.domain_bar == "off":
        return False
    return drawn >= (1 if look.domain_bar == "on" else 2)


def domain_bars(layout: Layout, look: Style | None = None) -> list[tuple[str, float, float, list[BarSegment]]]:
    """(chain, first residue, last residue, segments) for each chain holding domain elements (not the neighbours that
    only lend a few, Layout.partial), in chain order. The domain elements in sequence make the segments: a run of
    one domain's elements is one, from its first element's first residue to its last element's last; between two
    runs the loop is split at its middle, and the end segments reach the chain's first and last drawn residues, so
    the bar has no gap. Residue numbers as in the file."""
    owner = {k: d for d, (_, ids) in enumerate(layout.domains) for k in ids}
    elements = sorted((p.sse for k, p in layout.placed.items() if k in owner), key=lambda s: s.start)
    out = []
    for chain in dict.fromkeys(s.chain for s in elements if s.chain not in layout.partial):  # not a neighbour's few
        runs: list[list] = []  # [domain, first residue, last residue]
        for s in (s for s in elements if s.chain == chain):
            if runs and runs[-1][0] == owner[s.id]:
                runs[-1][2] = s.last.seq
            else:
                runs.append([owner[s.id], s.first.seq, s.last.seq])
        seqs = [r.seq for r in layout.res_labels if r.chain == chain]
        first, last = float(min([*seqs, runs[0][1]])), float(max([*seqs, runs[-1][2]]))
        cuts = [first, *((a[2] + b[1]) / 2 for a, b in zip(runs, runs[1:])), last]
        seen: Counter = Counter()
        segments = []
        for (d, _, _), lo, hi in zip(runs, cuts, cuts[1:]):
            name = layout.domains[d][0]
            seen[name] += 1
            gid = f"domain-bar:{chain}:{name}" + (f"#{seen[name]}" if seen[name] > 1 else "")
            segments.append(BarSegment(name, gid, lo, hi, *domain_paint(tone_of(layout, d), look or Style())))
        out.append((chain, first, last, segments))
    return out


def _spread(centres: list[float], widths: list[float], lo: float, hi: float, gap: float) -> list[float]:
    """Centres for texts of these widths, in order, as near the wanted ones as they can be while each keeps `gap`
    from the next and all stay within lo..hi (when they fit)."""
    xs = []
    for x, w in zip(centres, widths):  # left to right: clear of the previous one and of lo
        floor = lo + w / 2 if not xs else xs[-1] + (widths[len(xs) - 1] + w) / 2 + gap
        xs.append(max(x, floor))
    for k in range(len(xs) - 1, -1, -1):  # right to left: clear of the next one and of hi
        ceiling = hi - widths[k] / 2 if k == len(xs) - 1 else xs[k + 1] - (widths[k + 1] + widths[k]) / 2 - gap
        xs[k] = min(xs[k], ceiling)
    return xs


DOMAIN_SOURCE_CAPTIONS = {  # Layout.domain_source -> what the domain bar says of the domains (plan 12 R2)
    "inferred": "domains inferred from contacts",
    "CATH": "domains: CATH",
    "UniProt": "domains: UniProt",
}


def domain_source_caption(layout: Layout) -> str:
    """The domain bar's word on where its domains come from: inferred (automatic ones are a guess from contacts,
    not a published assignment), CATH or UniProt; '' for domains the user named."""
    return DOMAIN_SOURCE_CAPTIONS.get(getattr(layout, "domain_source", ""), "")


def _draw_domain_bar(c: _Canvas, layout: Layout) -> None:
    """Under the figure, one bar per chain holding domain elements (domain_bars), each in a band of c.bar_h, from the
    content's left edge on one residue scale for every bar (the longest chain spans the content width; Zyla 2024 Fig
    3A draws F2 and F1 so), each segment in its domain's tone, named inside when the name fits at 0.7 of the text
    size, else below on a short tick; the chain's end residue numbers over its ends, and between them the chain's
    name when there are several bars, and over the first bar where the domains come from (domain_source_caption)."""
    em = c.font / c.pt_per_unit  # page units of one em of the figure's text
    bars = domain_bars(layout, c.look)
    x0 = c.content[0]
    unit = (c.content[2] - x0) / max(max(last - first for _, first, last, _ in bars), 1.0)  # page units per residue
    for k, (chain, first, last, segments) in enumerate(bars):
        x1 = x0 + (last - first) * unit
        y = c.content[1] - c.pad - k * c.bar_h - 1.65 * em  # the bar's bottom; above it the residue numbers
        small = {"fontsize": c.font * 0.6, "va": "bottom"}
        for text, x, ha, end in ((first, x0, "left", "N"), (last, x1, "right", "C")):
            t = c.ax.text(x, y + em + 0.08, f"{text:g}", ha=ha, color=c.look.ink, **small)
            t.set_gid(f"domain-bar-end:{chain}:{end}")
        source = domain_source_caption(layout) if k == 0 else ""
        if len(bars) > 1 or source:
            text = " · ".join(t for t in (f"chain {chain}" if len(bars) > 1 else "", source) if t)
            t = c.ax.text(
                (x0 + x1) / 2,
                y + em + 0.08,
                text,
                ha="center",
                style="italic",
                color=legible("#5d6877"),
                **small,
            )
            t.set_gid(f"domain-bar-chain:{chain}" if len(bars) > 1 else "domain-bar-source")
        below = []  # segments too narrow for their name: (segment centre, name width, segment)
        for seg in segments:
            xa, xb = x0 + (seg.lo - first) * unit, x0 + (seg.hi - first) * unit
            rect = Rectangle(
                (xa, y),
                xb - xa,
                em,
                fc=seg.tone,
                ec=darken(seg.tone),
                lw=c.lw * 0.6,
                hatch=seg.hatch or None,
                zorder=Z["domain-bar"],
            )
            if seg.hatch and hasattr(rect, "set_hatchcolor"):  # matplotlib >= 3.11; earlier: the rc hatch colour
                rect.set_hatchcolor(darken(seg.tone, 0.8))
            rect.set_gid(seg.gid)
            c.ax.add_patch(rect)
            width = _text_width(seg.name) * 0.7 * em  # bold metrics: a little wider than the regular text drawn
            if width + 0.2 <= xb - xa:
                t = c.ax.text(
                    (xa + xb) / 2,
                    y + em / 2,
                    seg.name,
                    ha="center",
                    va="center",
                    fontsize=c.font * 0.7,
                    color=ink_on(seg.tone),
                    zorder=Z["domain-bar-label"],
                    bbox={"boxstyle": "square,pad=0.15", "fc": seg.tone, "ec": "none"} if seg.hatch else None,
                )
                t.set_gid(seg.gid.replace("domain-bar:", "domain-bar-label:", 1))
            else:
                below.append(((xa + xb) / 2, width, seg))
        xs = _spread([b[0] for b in below], [b[1] for b in below], c.page[0] + 0.1, c.page[1] - 0.1, 0.2)
        for (xc, _, seg), xt in zip(below, xs):
            (tick,) = c.ax.plot([xc, xt], [y, y - _BAR_TICK], color=darken(seg.tone), lw=c.lw * 0.6)
            tick.set_gid(seg.gid.replace("domain-bar:", "domain-bar-tick:", 1))
            t = c.ax.text(
                xt,
                y - _BAR_TICK - 0.04,
                seg.name,
                ha="center",
                va="top",
                fontsize=c.font * 0.7,
                color=domain_name_colour(seg.tone, c.look.background),
            )
            t.set_gid(seg.gid.replace("domain-bar:", "domain-bar-label:", 1))


def _copy_turn(layout: Layout, chain: str) -> float:
    """The turn (radians) a rosette (Layout.rosette) draws a chain's copy at about the page centre; 0 elsewhere."""
    if not getattr(layout, "rosette", 0):
        return 0.0
    k = layout.frame.sector.get(chain)
    return float(layout.frame.angles[k]) if k is not None else 0.0


def _turned_back(p: Placed, turn: float) -> Placed:
    """An element as its copy's own frame draws it, before a rosette's turn (itself when `turn` is 0)."""
    if not turn:
        return p
    c, s = np.cos(-turn), np.sin(-turn)
    return Placed(p.sse, c * p.cx - s * p.cy, s * p.cx + c * p.cy, p.length, p.angle - turn, p.width)


def _draw_sheet_panels(c: _Canvas, layout: Layout, outlines: dict[str, Path] | None = None) -> None:
    """A pale panel behind each sheet. With named domains (`outlines`, _draw_domains), a sheet whose strands are all
    one domain's takes that domain's tint, a little deeper than its panel; a sheet several domains share stays pale,
    with each of those domains' edges drawn over it, clipped to it, so the domains stay visible inside sheets."""
    if not c.look.sheet_panels:
        return
    owner = {i: d for d, (name, _) in enumerate(layout.domains) for i in panel_members(layout, name)}
    for k, ids in enumerate(layout.sheets):
        if len(ids) < 2:
            continue
        turn = _copy_turn(layout, layout.placed[ids[0]].sse.chain)  # a rosette's copy: its panel turned with it
        members = [layout.placed[i] for i in ids] + list(block_ghosts(layout, k))
        rects = np.array([_turned_back(m, turn).rect for m in members])
        x0, y0 = rects[:, 0].min() - 0.3, rects[:, 1].min() - 0.3
        x1, y1 = rects[:, 2].max() + 0.3, rects[:, 3].max() + 0.3
        fill = _mix(c.look.background, c.look.ink, 0.07) if c.look.background != "#ffffff" else _PANEL
        mine = sorted({owner[i] for i in ids if i in owner})
        names = [layout.domains[d][0] for d in mine if outlines and layout.domains[d][0] in outlines]
        if len(mine) == 1 and len(names) == 1 and all(i in owner for i in ids):
            fill = _mix(domain_paint(tone_of(layout, mine[0]), c.look)[0], c.look.background, 0.78)
        panel = FancyBboxPatch(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            boxstyle="round,pad=0,rounding_size=0.35",
            fc=fill,
            ec="none",
            zorder=Z["sheet-panel"],
        )
        if turn:
            panel.set_transform(Affine2D().rotate(turn) + c.ax.transData)
        panel.set_gid(f"sheet-panel:{k}")
        c.ax.add_patch(panel)
        if len(names) > 1:
            for d, name in zip(mine, names):
                tone = domain_paint(tone_of(layout, d), c.look)[0]
                edge = PathPatch(
                    outlines[name],
                    fc="none",
                    ec=_mix(tone, c.look.background, 0.3),
                    lw=c.lw * 0.8,
                    zorder=Z["sheet-domain-edge"],
                )
                edge.set_clip_path(panel)
                edge.set_gid(f"sheet-domain-edge:{k}:{name}")
                c.ax.add_patch(edge)


def body_paint(kind: str, color: str, look: Style, lw: float) -> tuple[str, str, str, float]:
    """(fill, shade, edge, edge width) of an element's body in `color` as look.fill and look.helix_shading paint it:
    shade is a helix's back face (the fill for a strand or a 3-10 helix); lw: the figure's line width. The figure's
    elements and the sequence view's track are both painted with it."""
    outline = look.fill == "outline"
    edge, edge_w = (color, lw * 1.2) if outline else (darken(color), lw * 0.8)
    if kind == "E":
        fill = look.background if outline else tint(color, 0.45) if look.fill == "pale" else color
        if look.fill == "pale":
            edge = darken(color, 0.5)  # pastel body, firm outline (illustration style)
        return fill, fill, edge, edge_w
    if kind == "G":
        fill = {"bold": color, "pale": tint(color, 0.35), "outline": look.background}[look.fill]
        return fill, fill, edge, edge_w
    fill, shade = {
        "bold": (color, darken(color, BACK_SHADE)),
        "pale": (tint(color, 0.15), tint(color, 0.55)),
        "outline": (look.background, _mix(color, look.background, 0.7)),
    }[look.fill]
    return fill, fill if look.helix_shading == "none" else shade, edge, edge_w


def eta_box(cx: float, cy: float, length: float, width: float, **kw) -> FancyBboxPatch:
    """A 3-10 helix's small rounded box, centred on (cx, cy) (the figure draws it at (0, 0) and rotates and moves it
    with its transform); kw: the patch's paint. The figure's 3-10 helices and the sequence view's are both this."""
    return FancyBboxPatch(
        (cx - length / 2, cy - width / 2), length, width, boxstyle=f"round,pad=0,rounding_size={width * 0.3}", **kw
    )


def _draw_elements(c: _Canvas, layout: Layout, bars=()) -> None:
    """The strands, 3-10 helices and helices, with their labels; a strand's letter slides along its shaft off the
    disulfide bars (`bars`: their routes) that pass under it (strand_label_spot)."""
    k = glyph(layout)[0] / _GLYPH_W  # helix names stand as far out as the text drawn needs (helix_label_pos)
    keep = link_keep(c.lw * 2.4, c.pt_per_unit, c.fig.dpi)
    c.fig.foldmap_moved_letters = []  # strand letters off their arrow's centre (provenance: names moved)
    for p in layout.placed.values():
        color = c.element_colour[p.sse.id]
        outline = c.look.fill == "outline"
        fill, shade, edge, edge_w = body_paint(p.sse.kind, color, c.look, c.lw)
        if p.sse.kind == "E":
            strand_fill = fill
            patch = Polygon(
                _arrow(p), closed=True, fc=strand_fill, ec=edge, lw=edge_w, joinstyle="round", zorder=Z["body"]
            )
            patch.set_gid(f"strand:{p.sse.id}")
            c.ax.add_patch(patch)
            if c.look.labels:
                spot = strand_label_spot(p, p.label, bars, c.font / c.pt_per_unit, keep)
                if spot != (p.cx, p.cy):
                    c.fig.foldmap_moved_letters.append(p.sse.id)
                t = c.ax.text(
                    *spot,
                    p.label,
                    ha="center",
                    va="center",
                    fontsize=c.font,
                    fontweight="bold",
                    color=_label_ink(color) if outline else text_color_on(strand_fill),
                    zorder=Z["label"],
                )
                t.set_gid(f"label:{p.sse.id}")
        elif p.sse.kind == "G":  # a short 3-10 helix: a small rounded box
            box = eta_box(0.0, 0.0, p.length, p.width, fc=fill, ec=edge, lw=edge_w, zorder=Z["body"])
            box.set_transform(Affine2D().rotate(p.angle).translate(p.cx, p.cy) + c.ax.transData)
            box.set_gid(f"eta:{p.sse.id}")
            c.ax.add_patch(box)
            if c.look.labels:
                t = c.ax.text(
                    *helix_label_pos(p, k),
                    p.label,
                    ha="center",
                    va="center",
                    fontsize=c.font * 0.9,
                    color=_label_ink(color),
                    zorder=Z["label"],
                )
                t.set_gid(f"label:{p.sse.id}")
        else:
            front, back, shine = _coil(p, len(p.sse))
            lift, alpha = (
                (_OVER_Z - Z["body"], _OVER_ALPHA) if p.over else (0.0, None)
            )  # over a sheet: on top, see-through
            fc_fill, fc_shade = fill, shade
            if p.over and all(c.element_colour.get(k) == color for k in layout.on_sheet.get(p.sse.id, ())):
                # the strands' own colour: a darker tone, fill partly see-through, edge opaque and darker still
                fill, shade, alpha = darken(fill, _OVER_TONE), darken(shade, _OVER_TONE), None
                fc_fill, fc_shade = to_rgba(fill, _OVER_SAME_ALPHA), to_rgba(shade, _OVER_SAME_ALPHA)
                edge = to_rgba(darken(color, 0.45), 1.0)
            behind = PathPatch(
                back,
                fc=fc_shade,
                ec=edge,
                lw=edge_w * 0.8,
                joinstyle="round",
                zorder=Z["body-back"] + lift,
                alpha=alpha,
            )
            behind.set_gid(f"helix-back:{p.sse.id}")
            c.ax.add_patch(behind)
            bar = PathPatch(
                front, fc=fc_fill, ec=edge, lw=edge_w, joinstyle="round", zorder=Z["body"] + lift, alpha=alpha
            )
            bar.set_gid(f"helix:{p.sse.id}")
            c.ax.add_patch(bar)
            if shine is not None and c.look.helix_shading == "gloss" and c.look.fill != "outline":
                gloss = PathPatch(
                    shine, fc=tint(fill, 0.4), ec="none", alpha=0.85 * (alpha or 1.0), zorder=Z["gloss"] + lift
                )
                gloss.set_gid(f"helix-shine:{p.sse.id}")
                c.ax.add_patch(gloss)
            if c.look.labels:
                t = c.ax.text(
                    *helix_label_pos(p, k),
                    p.label,
                    ha="center",
                    va="center",
                    fontsize=c.font,
                    color=_label_ink(color),
                    zorder=Z["label"] + lift,
                )
                if p.over:  # on the strands: a plate of the page colour keeps the name readable
                    t.set_bbox({"boxstyle": "round,pad=0.12", "fc": c.look.background, "ec": "none", "alpha": 0.85})
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
                    xy, ink = np.asarray(port) + inward * along + nrm * away * (p.width / 2 + 0.38), legible("#4a4a4a")
                elif p.length >= 3.0:  # room for both numbers and the strand letter between them
                    reach = 0.35 if end == "N" else min(_HEAD_LEN, p.length * 0.45) + 0.3  # clear of the head
                    xy, ink = (
                        np.asarray(port) + inward * reach,
                        darken(colour) if c.look.fill == "outline" else text_color_on(colour),
                    )
                else:  # too short to hold numbers: just past the end, beside the loop stub
                    xy, ink = np.asarray(port) - inward * 0.45 + nrm * 0.42, legible("#4a4a4a")
                t = c.ax.text(
                    *xy, str(num), ha="center", va="center", fontsize=c.font * 0.62, color=ink, zorder=Z["resnum"]
                )
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
            zorder=Z["body"],
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
                zorder=Z["label"],
            )
            t.set_gid(f"ghostlabel:{g.sse.id}")


def _draw_termini(c: _Canvas, layout: Layout, sses: list[SSE], loops: list[Loop], tails=None) -> None:
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
            edge = PathPatch(
                Path([port, tip]), fc="none", ec=_RUN_EDGE, lw=c.loop_lw * _RUN_EDGE_W, zorder=Z["stub-edge"]
            )
            edge.set_gid(f"stub-edge:{end}:{chain}")
            c.ax.add_patch(edge)
        stub = PathPatch(
            Path([port, tip], [Path.MOVETO, Path.LINETO]), fc="none", ec=ink, lw=c.loop_lw, zorder=Z["stub"]
        )
        stub.set_gid(f"stub:{end}:{chain}")
        c.ax.add_patch(stub)
        t = c.ax.text(
            *end_label_pos(layout, end, chain, port, ex),
            end,
            ha="center",
            va="center",
            fontsize=c.font,
            fontweight="bold",
            zorder=Z["label"],
        )
        t.set_gid(f"terminus:{end}:{chain}")
    for count in (
        tails if tails is not None else tail_counts(layout, sses, loops)
    ):  # a long unstructured end: how long, not taken for a short stub
        t = c.ax.text(
            *count["centre"],
            count["text"],
            ha="center",
            va="center",
            fontsize=c.font * TAIL_COUNT_SIZE,
            color=legible("#5d6877"),
            zorder=Z["label"],
        )
        t.set_gid(f"tail-length:{count['end']}:{count['chain']}")


_NA_ENTRY = re.compile(r"(DNA|RNA|DNA/RNA) \(")  # the legend text of a nucleic-acid backbone


def _draw_legend(c: _Canvas, layout: Layout, legend: _Legend) -> None:
    """The legend rows and, for a large assembly, the context caption under them."""
    xmin, xmax, ymin, ymax = c.page
    for k, (colour, text) in enumerate(legend.entries):
        row, col = divmod(k, legend.per_row)
        x, y = xmin + c.pad + col * legend.entry, ymin + legend.height - 1.1 - row * legend.row
        if isinstance(colour, tuple) and colour[0] == "grades":  # the nine grades in ConSurf's colours
            for q, shade in enumerate(CONSURF):
                c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none"))
            c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5))
        elif isinstance(colour, tuple) and colour[0] == "ia":  # an interaction kind: a short line, as drawn
            ink, dash = interaction_ink(colour[1], c.look)
            (mark,) = c.ax.plot([x, x + 0.6], [y, y], color=ink, lw=c.loop_lw * 0.85, ls=dash)
            mark.set_gid(f"legend-ia:{colour[1]}")
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
            (rail,) = c.ax.plot(  # as thick as the backbone is drawn (Style.na_line)
                [x, x + 0.6], [y, y], color=colour, lw=c.lw * 2.4 * c.look.na_line, solid_capstyle="butt"
            )
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
                    color=legible(colour),
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
            color=legible("#5d6877"),
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
    at=None,
) -> Figure:
    look = _resolve_look(look, palette, style, loop_style)
    chains = _chains_in_order(sses)
    colors = chain_colors(chains, look.palette)
    element_colour, residue_colour = element_colours(layout, sses, look)
    if not layout.placed and not layout.stems and not layout.na_paths:
        return _empty_figure()
    hidden: list[dict] = []  # counts with no clear spot, left out (provenance: counts_hidden)
    drawn = drawn_lines(layout, loops, look)  # the loops' lines as drawn, which the counts keep off
    as_drawn = [SimpleNamespace(points=pts) for pts in drawn]
    tails = tail_counts(layout, sses, as_drawn, hidden)  # placed once: every text that keeps off them gets these boxes
    lcounts = loop_counts(layout, loops, sses, tails, hidden, drawn)
    counts = [t["box"] for t in [*tails, *lcounts]]
    marks = ligand_marks(layout, loops, sses, look, counts)
    linked = linked_marks(layout, loops, sses, marks)
    umarks = uniprot_marks(layout, loops, sses, look)
    hb = interaction_overlay(layout, loops, sses, look, model)
    c, legend = _canvas(
        layout,
        look,
        title,
        chains,
        colors,
        element_colour,
        residue_colour,
        marks,
        umarks,
        hb,
        linked,
        loops,
        counts,
        at,
    )
    _draw_loops(c, layout, loops, sses, lcounts)
    _draw_membrane(c, layout)
    _draw_axis_symbol(c, layout)
    outlines = _draw_domains(c, layout)
    _draw_sheet_panels(c, layout, outlines)
    _draw_hbonds(c, hb.hbonds, hb.slots)  # the other kinds after every text: they keep off the letters
    c.fig.foldmap_hbonds = hb.hbonds  # the H-bond lines (the interactive page reads foldmap_interactions)
    c.fig.foldmap_interactions = hb  # every line drawn, and those the re-check flagged (provenance)
    c.fig.foldmap_counts_hidden = hidden
    c.fig.foldmap_crowded_marks = [  # ligand marks with no clear spot, drawn over other ink (Drawn.crowded_marks)
        f"{m['ligand'].name}:{m['ligand'].chain}{m['ligand'].seq}" for m in marks if m.get("crowded")
    ]
    if layout.stems or layout.na_paths:
        _draw_nucleic(
            c.ax, layout, sses, c.residue_colour, c.look.fill, c.lw, c.font, c.pt_per_unit, c.look, c.backbone[0]
        )
    _draw_residue_numbers(c, layout)
    _draw_termini(c, layout, sses, loops, tails)
    # the bars once: a loop's cysteine off the text drawn so far (loop lengths, chain ends, residue numbers); a
    # strand's letter then slides off the bars that pass under it (not when residue numbers share its shaft)
    keep = link_keep(c.lw * 2.4, c.pt_per_unit, c.fig.dpi)
    avoid = [*text_halos(c.ax, c.pt_per_unit, keep), *link_bodies(layout)]  # a loop's cysteine shows: not under
    bars = disulfide_bars(layout, loops, sses, look, c.lw, c.pt_per_unit, avoid)  # a letter, nor under a body
    _draw_elements(c, layout, [] if look.residue_numbers else [route for *_, route in bars])
    _draw_ghosts(c, layout)
    # after every text but the marks' own: bars and tethers keep off the letters (text_halos)
    _draw_links(c.ax, layout, loops, sses, look, c.lw, c.pt_per_unit, bars=bars, counts=counts)
    _draw_ligands(c.ax, marks, c.lw, c.font, layout, c.pt_per_unit)
    _draw_linked(c.ax, linked, c.lw, c.font)
    _draw_uniprot(c.ax, umarks, layout, c.lw, c.font)
    c.fig.foldmap_interactions_hidden = _draw_interactions(c, hb)  # lines wholly under text (provenance)
    if c.bar_h:
        _draw_domain_bar(c, layout)
    if at is None:  # a panel's legend is its page's (_canvas)
        _draw_legend(c, layout, legend)
    if title:
        _draw_title(c, title)
    # how it was made, last: what is inferred, ambiguous or moved (Drawn.provenance; save writes it into the file)
    c.fig.foldmap_provenance = figure_provenance(c.fig, layout, loops, labels=c.look.labels)
    return c.fig


def save_svg(fig: Figure) -> str:
    """The figure as an inline SVG element (no XML prolog), with the same settings as save()."""
    import io as _io

    buf = _io.StringIO()
    with rc_context({"svg.fonttype": "none", "svg.hashsalt": "foldmap", **_RC}):
        fig.savefig(buf, format="svg", metadata={"Date": None}, facecolor=fig.get_facecolor())
    text = with_desc(buf.getvalue(), getattr(fig, "foldmap_provenance", []))
    return text[text.index("<svg") :]


PNG_DPI = 300  # a PNG's resolution ...
MAX_PNG_PIXELS = 50_000_000  # ... up to this many pixels (200 MB to draw): a larger page's PNG takes a lower dpi


def png_dpi(fig: Figure, dpi: int = PNG_DPI) -> int:
    """The dpi save() writes the figure's PNG at: `dpi`, or less on a page so large (random-100 C4: past
    MAX_WIDTH_IN the page grows to keep its text at 6 pt) that it would take more than MAX_PNG_PIXELS."""
    w, h = fig.get_size_inches()
    return int(min(dpi, (MAX_PNG_PIXELS / (w * h)) ** 0.5))


def page_notes(fig: Figure) -> list[str]:
    """What the notes (CLI, app, MCP) say of the figure's page: a page held at layout.MAX_PAGE_IN with its text under
    6 pt (layout.page_fit), and a PNG written at less than PNG_DPI (png_dpi)."""
    out = [fit] if (fit := getattr(fig, "foldmap_page_fit", "")) else []
    if (dpi := png_dpi(fig)) < PNG_DPI:
        w, h = fig.get_size_inches()
        out.append(f"its PNG is {dpi} dpi, not {PNG_DPI}: the page is {w:.0f} x {h:.0f} in (SVG and PDF keep it whole)")
    return out


def save(fig: Figure, path: str | FsPath, dpi: int = PNG_DPI) -> FsPath:
    """The figure as .svg, .pdf or .png (at `dpi`, at most png_dpi's), with how it was made (Drawn.provenance) in
    the file: the SVG's <desc>, the PDF's Subject, the PNG's Description (none when there is nothing to say)."""
    import io as _io

    path = FsPath(path)
    fmt = _EXT.get(path.suffix.lower())
    if fmt is None:
        raise ValueError(f"unsupported extension {path.suffix!r}; use .svg, .pdf or .png")
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {"svg": {"Date": None}, "pdf": {"CreationDate": None}, "png": {"Software": None}}[fmt]
    items = getattr(fig, "foldmap_provenance", [])
    if items and fmt != "svg":
        meta[{"pdf": "Subject", "png": "Description"}[fmt]] = as_text(items)
    if fmt == "png":
        dpi = min(dpi, png_dpi(fig))
    with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42, "svg.hashsalt": "foldmap", **_RC}):
        if fmt == "svg":  # as matplotlib writes it, with the <desc> added, and LF line ends on every platform
            buf = _io.StringIO()  # (matplotlib, given a path, writes in text mode: CRLF on Windows)
            fig.savefig(buf, format=fmt, dpi=dpi, metadata=meta, facecolor=fig.get_facecolor())
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(with_desc(buf.getvalue(), items) if items else buf.getvalue())
        else:
            fig.savefig(path, format=fmt, dpi=dpi, metadata=meta, facecolor=fig.get_facecolor())
    return path
