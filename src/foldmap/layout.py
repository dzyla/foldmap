# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Place sheet blocks and helices on the page (page units: 1 unit = one strand pitch)."""

from __future__ import annotations

import itertools
import logging
import math
import re
from collections import Counter
from collections.abc import Callable
from contextvars import ContextVar
from dataclasses import KW_ONLY, dataclass, field
from types import SimpleNamespace
from typing import NamedTuple

import numpy as np

from . import panels, text_metrics
from .features import KINK_GAP, helix_bundles
from .frame import Frame, symmetric_frame, view_frame
from .model import SSE, Nucleic, Sheet, generic_label, group_by_chain, pair_places, same_place
from .nucleic import Stem, StemView, segment_stem, stem_axis, stem_views
from .style import CHOICES, Style

log = logging.getLogger(__name__)

PITCH = 1.1  # distance between neighbouring strand centres
SHAFT = 0.62  # arrow shaft width
SCALE = 0.23  # page units per Å (PITCH / 4.8 Å strand spacing)
HEAD = 1.0  # strand arrowhead width
HEAD_LEN = 0.7  # strand arrowhead length (at most 0.45 of a short strand)
HELIX_W = 1.0  # helix ribbon diameter, ~4.4 Å: one arrowhead wide, as in a cartoon (Style.helix_scale)
STRAND_RISE = 3.3  # Å per residue along a strand
HELIX_RISE = 1.5  # Å per residue along a helix
TEN_RISE = 2.0  # Å per residue along a 3-10 helix
ETA_W = 0.7  # 3-10 box width as a fraction of the helix width
MARGIN = 0.8  # clear space kept around each element so loops can pass
_ROW_GAP = 2.5  # between neighbours in a stack row
_ROW_SPACING = 3.2  # between stack rows: room for a chain-end label below one row and above the next
LABEL_GAP = 0.3  # a helix label's centre sits this far beyond the side of the element box
END_STUB = 0.7  # chain-end stub length
END_LABEL = 1.15  # distance from a chain-end port to its N/C letter
LONG_TAIL = 10  # a chain end with this many residues or more past its last element gets a count ("14 aa")
TAIL_COUNT_SIZE = 0.72  # that count's text size, as a fraction of the figure's
_TAIL_COUNT_GAP = 0.08  # page units between the N/C letter's ink and the count's
_TAIL_COUNT_FAR = (0.3, 0.6, 0.9, 1.2)  # page units further out a count looks when its near spots are all taken
_GLYPH_W, _GLYPH_H = 0.32, 0.42  # rough label glyph size in page units, for keep-out boxes (at _EM_TUNED: glyph_scale)
_GLYPH_K: ContextVar[float] = ContextVar("foldmap_glyph_scale", default=1.0)  # glyph_scale of the text the layout
# being built is drawn with (build_layout sets it from Style.font_scale): for the glyph room of things that hold no
# layout (an item's footprint, a stem's name box, a trace's letters); layout.glyph(layout) wherever a layout is held
IN_PER_UNIT = 0.30  # render: inches per page unit, at most
MAX_WIDTH_IN = 7.2  # render: a figure is drawn at most this wide (double column), inches, while its text stays ...
MIN_TEXT_PT = 6.0  # ... this size or larger: a wider page grows past MAX_WIDTH_IN instead (random-100 C4 ruling)
MAX_PAGE_IN = 48.0  # render: the longest side a page grows to (A0's, about); past it the whole figure, its text
# with it, is scaled to fit: text under MIN_TEXT_PT, said by the figure's notes (page_fit)
_MAX_ROW = 30.0
_MIN_STRAND = 1.0  # shortest drawn strand (room for the arrow head)
# items: at or below it the dense energy keeps recorded layouts bit-identical (the pair-list summation order
# moves L-BFGS minima)
_PAIR_ENERGY_ABOVE = 200
_W_SEQ = 1.0  # pull between the two ends of a loop
_LOOP_SLACK = 1.5  # loop ends closer than this are left alone
_W_SS = 2.0  # pull between two bridged cysteines (twice a loop: disulfides should read as short bars)
_SS_SLACK = 0.6
_W_CONTACT = 0.3  # pull between elements that touch in 3D (at _CONTACT_CAP CA pairs)
_CONTACT_CAP = 20
_K_HOME = 0.02  # pull back toward the projected position, keeps the fold recognisable
_STEP = 0.15
_SWEEPS = 400
_SEP_SWEEPS = 5000  # push-apart passes before the guaranteed slide
_SEP_ORBIT_SWEEPS = 1000  # push-apart passes that keep symmetry copies alike, before they may move alone
_SEP_SYM_ROOM = 1.05  # copies drawn alike may take at most this much more area than items moved one by one
# ... and copies drawn as mirror images about a dyad in the page, pushed for as many passes as items moved one by
# one: the figure says they are mirror images (9GEP's nucleosome needed 1.07 and about 3,000 passes)
_SEP_MIRROR_ROOM = 1.10
_SEP_PUSH = 0.005  # each push clears an overlap by this much more, so touching boxes cannot trade it back and forth
_SLIDE_STEP = 0.5  # the guaranteed slide moves an item by this much at a time
_LAMBDAS = (0.3, 3.0, 30.0, 300.0)  # overlap penalty, stage by stage
_MAX_MOVE = 0.5  # largest step an item takes in one sweep
_GRID = 1e4  # optimiser inputs and item geometry are rounded to 1/_GRID page units (radians for angles)


@dataclass(eq=False)
class Placed:
    sse: SSE
    cx: float
    cy: float
    length: float
    angle: float  # direction N -> C, radians (0 = right, pi/2 = up)
    width: float
    label: str = ""
    ghost: bool = False  # repeat of a barrel's first strand, drawn to show the sheet closes
    ends: str = ""  # "N", "C" or "NC" when the chain starts/ends here (its stub and letter need room)
    label_side: int = 1  # a helix's label: 1 on its upper side, -1 on the other (_choose_label_sides)
    over: bool = False  # a helix packed on a sheet face, drawn over the sheet's strands (_over_sheet)
    label_along: float = 0.0  # its name moved this far along its axis, N -> C (_over_label_sides)

    @property
    def direction(self) -> np.ndarray:
        return np.array([np.cos(self.angle), np.sin(self.angle)])

    @property
    def n_port(self) -> tuple[float, float]:
        p = np.array([self.cx, self.cy]) - self.direction * self.length / 2
        return float(p[0]), float(p[1])

    @property
    def c_port(self) -> tuple[float, float]:
        p = np.array([self.cx, self.cy]) + self.direction * self.length / 2
        return float(p[0]), float(p[1])

    @property
    def exit_c(self) -> np.ndarray:
        return self.direction

    @property
    def exit_n(self) -> np.ndarray:
        return -self.direction

    @property
    def rect(self) -> tuple[float, float, float, float]:
        hx = abs(np.cos(self.angle)) * self.length / 2 + abs(np.sin(self.angle)) * self.width / 2
        hy = abs(np.sin(self.angle)) * self.length / 2 + abs(np.cos(self.angle)) * self.width / 2
        return (self.cx - hx, self.cy - hy, self.cx + hx, self.cy + hy)


DNA_RISE = 3.4  # Å per base pair (B-DNA)
DNA_W = 20.0 * SCALE  # B-DNA diameter
DNA_PITCH = 10.5 * DNA_RISE * SCALE  # one turn, page units
DNA_MINOR = 0.37  # the second backbone trails the first by this fraction of a turn (minor groove)


NA_W = DNA_W  # stem width: a double helix is about 20 Å across, A-form RNA as B-DNA
RAIL_IN = 0.6  # ladder rails sit this far inside the stem's long sides, base letters outside them
NT_STEP = 0.55  # spacing of unpaired nucleotides along a drawn single strand (schematic, not to scale)


@dataclass(eq=False)
class PlacedStem:
    """A stem (base-paired double helix: an RNA stem, a DNA duplex, or one straight piece of a bent duplex)
    drawn as a rigid element. axial[(strand, nt)]: the nucleotide's position along the drawn axis from the
    centre, page units, positive along `direction`. side[(strand, nt)]: +1 on the first arm's rail (left of
    `direction`), -1 on the second arm's. The first arm runs 5'->3' along `direction`, the second back."""

    id: str
    name: str
    strands: list[int]  # [first arm's strand, second arm's]; equal for an RNA hairpin stem
    pairs: list[tuple[int, int, int, int]]
    cx: float
    cy: float
    length: float
    angle: float
    width: float
    axial: dict[tuple[int, int], float]
    side: dict[tuple[int, int], int]
    bulges: set = field(default_factory=set)
    parent: str = ""  # a piece of a bent duplex: the whole stem's id
    first: int = 1  # the side the first arm's rail is on
    margin: float = 0.0  # Å the first arm is seen beside the second (_side_signal): how clearly the side is seen
    hairpin: int = -1  # nucleotides of the hairpin loop over the stem's inner tip (-1: no hairpin loop there)
    ghost: bool = False
    label_side: int = 1  # its name beside the first rail (1) or the second (-1): _choose_label_sides
    label_out: float = 0.0  # and this much further out than the first slot (_choose_label_sides)
    label_at: float | None = None  # its name's centre this far from the axis, beside its letters (_choose_label_sides)
    label_along: float = 0.0  # and this far along it from its centre, midway along its row (_choose_label_sides)
    label_keep: tuple[float, float, float, float] | None = None  # its name's ink box grown by _INK_MARGIN em, around
    # the name's centre: what routes keep out of (stem_keep_box); None until _choose_label_sides

    @property
    def direction(self) -> np.ndarray:
        return np.array([np.cos(self.angle), np.sin(self.angle)])

    @property
    def normal(self) -> np.ndarray:
        return np.array([-np.sin(self.angle), np.cos(self.angle)])

    @property
    def rail(self) -> float:
        """Distance of each ladder rail from the axis."""
        return max(self.width / 2 - RAIL_IN, self.width / 4)

    def point(self, key: tuple[int, int], off: float = 0.0) -> np.ndarray:
        """Page position of a nucleotide on its rail, or `off` further out from the axis."""
        r = (self.rail + off) * self.side[key]
        return np.array([self.cx, self.cy]) + self.axial[key] * self.direction + r * self.normal

    def port(self, key: tuple[int, int], leaving: bool) -> tuple[np.ndarray, np.ndarray]:
        """Where the strand leaves (or enters) the stem at this nucleotide, and the outward direction there:
        along the arm's own 5'->3' sense when leaving, against it when entering."""
        sense = 1 if self.side[key] == self.first else -1  # the first arm runs along `direction`, the second back
        along = self.direction * sense * (1 if leaving else -1)
        return self.point(key) + along * NT_STEP / 2, along

    rect = Placed.rect

    def hairpin_circle(self, ka: tuple[int, int], kb: tuple[int, int], m: int) -> tuple[np.ndarray, float]:
        """(centre, radius) of the near-circle a hairpin loop of m nucleotides from ka to kb is drawn on, beyond
        the stem's tip."""
        pa, pb = self.point(ka), self.point(kb)
        out = self.direction if self.axial[ka] >= 0 else -self.direction
        half = float(np.hypot(*(pb - pa))) / 2
        r = max(half * 1.02, (m + 1) * NT_STEP / (1.6 * np.pi))
        return (pa + pb) / 2 + out * np.sqrt(max(r * r - half * half, 0.0)), r


@dataclass(eq=False)
class Layout:
    placed: dict[str, Placed]
    bounds: tuple[float, float, float, float]
    _: KW_ONLY  # the rest by name only
    ghosts: list[Placed] = field(default_factory=list)
    stems: list[PlacedStem] = field(default_factory=list)
    nucleic: Nucleic | None = None
    sheets: list[list[str]] = field(default_factory=list)  # strand ids of each sheet block, in layout order
    sheet_ghosts: list[list[str]] = field(default_factory=list)  # ghost strand ids drawn in each block of `sheets` (a
    # barrel's repeated strand, a strand lent across the cut of an unrolled ring: _unwrap); [] when not recorded
    domains: list[tuple[str, list[str]]] = field(default_factory=list)  # named domains and their element ids
    domain_lent: dict[str, list[str]] = field(default_factory=dict)  # domain -> elements a neighbour lends to it
    # (focus mode, chains in `partial`, drawn grey): inside its panel (panel_members), not its members (random-100 C23)
    domain_names: dict[str, tuple[float, float, float, float]] = field(default_factory=dict)  # each drawn domain's
    # name box (_place_domain_names), on the top edge of its panel or, with a leader, above it
    domain_leaders: dict[str, tuple] = field(default_factory=dict)  # domain -> leader (from the name, to its panel)
    domain_tones: list[int] = field(default_factory=list)  # tone of each domain (features.domain_tones; [] by order)
    na_paths: list = field(default_factory=list)  # route.NAPath: single strands as drawn (filled by route_loops)
    na_traces: dict = field(default_factory=dict)  # strand -> projected backbone (page units), for unpaired strands
    na_ends: list = field(default_factory=list)  # (strand, "5′" or "3′", label position), filled by route_loops
    na_letters: dict = field(default_factory=dict)  # single-stranded (strand, nt) -> its letter's place (route_loops)
    na_ties: list = field(default_factory=list)  # na_paths.NATie: isolated pairs as drawn (route_loops)
    na_cropped: bool = False  # na_extent auto cropped a large nucleic acid to the stems the protein touches
    na_read: int = 0  # nucleotides read before that crop (scope.drawn_scope says how many of them were drawn)
    frame: object = None  # the frame.Frame (or RadialFrame, CylinderFrame) that took 3D to the page
    focus: set[str] | None = None  # chains drawn in colour when highlighting; the rest in grey
    links: object = None  # model.Links: disulfides and glycans, residue indices in Backbone order
    res_chain: list[str] = field(default_factory=list)  # chain of each Backbone residue
    res_name: list[str] = field(default_factory=list)  # residue type of each Backbone residue
    res_b: list[float] = field(default_factory=list)  # CA B-factor of each Backbone residue
    res_labels: list = field(default_factory=list)  # model.ResLabel of each Backbone residue
    res_linked: list[int] = field(default_factory=list)  # residues fused into the chain (a chromophore), sorted
    focus_chains: set | None = None  # large assemblies: the subunit drawn in full (None: everything drawn)
    unjoined: set = field(default_factory=set)  # (a, b) element pairs shown side by side but not adjacent in sequence
    partial: set = field(default_factory=set)  # chains present only as the elements that complete its fold
    context: str | None = None  # what the figure shows of the assembly, e.g. "1 of 6 subunits · helical ..."
    subunits: tuple[int, str] | None = None  # focus on a symmetric assembly: its subunits and symmetry (sym.label,
    # or "helical")
    chains_read: list[str] = field(default_factory=list)  # every chain read and asked for (protein, then nucleic
    # acid), before focus kept one subunit
    chains_asked: list[str] = field(default_factory=list)  # --chains: the protein chains asked for (scope says so)
    page_axes: list | None = None  # 3x3 rows: page-right, page-up, toward-viewer (for the first protomer)
    page_centre: list | None = None  # the 3D point those axes are taken around
    uniprot: object = None  # uniprot.Annotation when UniProt features were requested
    membrane: dict | None = None  # {"y": (bottom, top) of the band, "inside": "below"|"above", "source": ...}
    # with "bands": one (bottom, top) per stack row ("y" is the first)
    stack_rows: list[list[str]] = field(default_factory=list)  # a depth-keeping stack wrapped: each row's ids
    ba_rows: list[list[str]] = field(default_factory=list)  # each βα template row's elements (no ghost), as drawn
    # left to right: the strands and the helices between them (_ba_template)
    on_sheet: dict[str, list[str]] = field(default_factory=dict)  # a helix drawn over a sheet -> the sheet's strands
    # the gaps between two strands sharing a column (Sheet.columns), as boxes (_stack_gaps): the loops leaving
    # their facing ends are kept apart there (route._apart_in_stack_gaps)
    stack_gaps: list[tuple] = field(default_factory=list)
    res_cons: list[float] = field(default_factory=list)  # alignment conservation of each residue (nan: none)
    res_grade: list[int] = field(default_factory=list)  # conservation grade 1..9 of each residue (0: none)
    cons_label: str = ""  # legend text saying where the grades come from
    copy_ties: dict[str, tuple[str, bool]] = field(default_factory=dict)  # copy chain -> (its protomer-0 mate chain,
    # drawn as its mirror image): the symmetry ties asked for (_symmetry_ties), kept when _separate moves one by one
    mirror_pairs: list[tuple[str, str]] = field(
        default_factory=list
    )  # (copy, mate) helix ids of mirror-image protomers
    tied_ids: set[str] = field(default_factory=set)  # ids of the elements of every item a symmetry tie holds, copy
    # and mate alike (a sheet block among them is drawn as its mate is: figcheck.shared_sheet_elements)
    font_scale: float = 1.0  # Style.font_scale the figure is drawn with (texts placed in the layout are sized by it)
    helix_matched: dict = field(default_factory=dict)  # {"cut", "joined"}: helix pieces made alike in copies (ss)
    na_letter_rows: tuple[float | None, bool] = (None, False)  # where a stem's base letters sit (NT_ROW, 0 on the
    # rails, None: not drawn) and whether single-stranded letters are drawn, as render draws them (na_text_boxes)
    ligands: str = "auto"  # Style.ligands the figure is drawn with (auto hides additives: scope.drawn_scope names them)
    coil_chains: list[tuple[str, int]] = field(default_factory=list)  # (chain, residues) of each protein chain read
    # with no helix or strand: nothing of it can be drawn (scope.drawn_scope says so)

    @property
    def dna(self) -> list[PlacedStem]:  # earlier name
        return self.stems

    @property
    def width(self) -> float:
        return self.bounds[2] - self.bounds[0]

    @property
    def height(self) -> float:
        return self.bounds[3] - self.bounds[1]


@dataclass(eq=False)
class _Item:
    """A sheet block or a single helix, moved as a unit. members hold offsets from (cx, cy)."""

    members: list[Placed]
    offsets: list[tuple[float, float]]
    cx: float
    cy: float
    home: tuple[float, float]
    order: int  # first residue index, for stack mode
    chain: str  # chain of the lowest-numbered element, for stack mode

    @property
    def box(self) -> tuple[float, float, float, float]:
        """Footprint relative to (cx, cy): every member with its label room."""
        fp = [_footprint(m) for m in self.members]
        return (
            min(o[0] + f[0] for o, f in zip(self.offsets, fp)),
            min(o[1] + f[1] for o, f in zip(self.offsets, fp)),
            max(o[0] + f[2] for o, f in zip(self.offsets, fp)),
            max(o[1] + f[3] for o, f in zip(self.offsets, fp)),
        )

    @property
    def half(self) -> tuple[float, float]:
        x0, y0, x1, y1 = self.box
        return (x1 - x0) / 2, (y1 - y0) / 2

    @property
    def mid(self) -> tuple[float, float]:
        """Centre of the footprint relative to (cx, cy)."""
        x0, y0, x1, y1 = self.box
        return (x0 + x1) / 2, (y0 + y1) / 2

    def apply(self) -> None:
        for m, (ox, oy) in zip(self.members, self.offsets):
            m.cx, m.cy = self.cx + ox, self.cy + oy


def _footprint(m, k: float | None = None) -> tuple[float, float, float, float]:
    """A member's box relative to its own centre; a helix's or a stem's includes the room its label needs, glyphs
    `k` (glyph_scale; the layout being built's, _GLYPH_K, when None) times _GLYPH_W x _GLYPH_H."""
    k = _GLYPH_K.get() if k is None else k
    gw, gh = _GLYPH_W * k, _GLYPH_H * k
    x0, y0, x1, y1 = m.rect
    box = (x0 - m.cx, y0 - m.cy, x1 - m.cx, y1 - m.cy)
    if isinstance(m, Placed) and m.sse.kind in ("H", "G"):
        lx, ly = helix_label_pos(Placed(m.sse, 0.0, 0.0, m.length, m.angle, m.width))
        hx, hy = gw * 1.5 + 0.05, gh / 2  # up to three glyphs: α10
        box = (min(box[0], lx - hx), min(box[1], ly - hy), max(box[2], lx + hx), max(box[3], ly + hy))
    if isinstance(m, PlacedStem):  # a stem's name beside it
        x0, y0, x1, y1 = _stem_label_box(m, k)
        box = (min(box[0], x0 - m.cx), min(box[1], y0 - m.cy), max(box[2], x1 - m.cx), max(box[3], y1 - m.cy))
    if isinstance(m, PlacedStem) and m.hairpin >= 0:  # and its hairpin loop over the tip, with the loop's letters
        a, i, b, j = m.pairs[-1]
        (x, y), r = m.hairpin_circle((a, i), (b, j), m.hairpin)
        x, y, r = x - m.cx, y - m.cy, r + NT_ROW + gh / 2
        box = (min(box[0], x - r), min(box[1], y - r), max(box[2], x + r), max(box[3], y + r))
    if isinstance(m, Placed) and m.ends:
        d = m.direction
        for end in m.ends:
            out = d if end == "C" else -d
            t = out * (m.length / 2 + END_LABEL)  # the N/C letter beyond the stub
            hx, hy = gw / 2 + 0.05, gh / 2
            box = (min(box[0], t[0] - hx), min(box[1], t[1] - hy), max(box[2], t[0] + hx), max(box[3], t[1] + hy))
    return box


def _strand_length(sse: SSE) -> float:
    return max(_MIN_STRAND, (len(sse) - 1) * STRAND_RISE * SCALE)


def _anchor(strands: list[SSE]) -> tuple[int, str]:
    first = min(strands, key=lambda s: s.start)
    return first.start, first.chain


def _sheet_order(sheet: Sheet, frame: Frame) -> tuple[list[int], int, np.ndarray]:
    """How a sheet is drawn: the strands' indices left to right, s = +1 when strands[0] points up (else -1), and
    the strands' projected centroids (Å)."""
    strands = sheet.strands
    cents = frame.project(np.array([s.centroid for s in strands]), [s.chain for s in strands])
    mean_dir = sum(d * s.axis for d, s in zip(sheet.directions, strands))
    at = np.mean([x.centroid for x in strands], axis=0)
    seen = frame.direction(at, mean_dir)
    s = 1 if seen[1] >= 0 else -1  # +1: strands[0] points up
    # Drawing the strands vertical turns the projected sheet in the page plane; the strand order must turn
    # with it (left to right = across the drawn-up direction), or the sheet shows its back face (mirror image).
    up2d = s * seen
    right = np.array([up2d[1], -up2d[0]]) if np.linalg.norm(up2d) > 1e-6 else np.array([1.0, 0.0])
    order = list(range(len(strands)))
    if len(strands) > 1 and (cents[-1] - cents[0]) @ right < 0:
        order.reverse()
    return order, s, cents


_STACK_PAD = 1.2  # page units between two strands drawn one above the other in one column: room for the loops
# leaving their facing ends (0.2 left 5KQJ's B and F loops no route, 1.0 5JB1's; set on the 5 sample entries)


def _stack_gaps(lay: Layout) -> list[tuple[float, float, float, float]]:
    """The gaps between two strands drawn one above the other in a column (Layout.stack_gaps): the box between their
    facing ends, the column's width wide. route._apart_in_stack_gaps keeps the two loops leaving those ends apart
    there, so the column never reads as one strand (5KQJ's B and F)."""
    out = []
    for ids in lay.sheets:
        strands = [lay.placed[k] for k in ids]
        for p in strands:
            above = [q for q in strands if q is not p and abs(q.cx - p.cx) < 1e-6 and q.cy > p.cy]
            if above:
                q = min(above, key=lambda x: x.cy)
                top, bottom = max(p.n_port[1], p.c_port[1]), min(q.n_port[1], q.c_port[1])
                out.append((p.cx - p.width / 2, top, p.cx + p.width / 2, bottom))
    return out


def _column_heights(sheet: Sheet, order: list[int], s: int) -> dict[int, float]:
    """The heights (page units) of the strands that share a column (Sheet.columns), each as its centroid lies along
    the sheet's drawn up (its strands' mean direction, by their senses, times the flip s) from the column's middle,
    moved apart evenly from that middle until they keep _STACK_PAD between them; snapped."""
    up3 = _unit(s * sum(d * x.axis for d, x in zip(sheet.directions, sheet.strands)))
    out: dict[int, float] = {}
    for c in sorted(set(sheet.columns)):
        ks = [k for k in order if sheet.columns[k] == c]
        if len(ks) < 2:
            continue
        cents = [_snap(sheet.strands[k].centroid) for k in ks]
        mid = _snap(np.mean(cents, axis=0))
        ys = sorted(((float(_snap((x - mid) @ up3 * SCALE)), k) for x, k in zip(cents, ks)))
        placed = [ys[0][0]]
        for (_, a), (y, b) in zip(ys, ys[1:]):  # each above the one below it, at least its reach apart
            need = (_strand_length(sheet.strands[a]) + _strand_length(sheet.strands[b])) / 2 + _STACK_PAD
            placed.append(max(y, placed[-1] + need))
        shift = (ys[0][0] + ys[-1][0]) / 2 - (placed[0] + placed[-1]) / 2  # the column's middle kept where it was
        out.update({k: float(_snap(y + shift)) for y, (_, k) in zip(placed, ys)})
    return out


def _sheet_item(
    sheet: Sheet, frame: Frame, pitch: float = PITCH, over: list[SSE] = (), helix_w: float = HELIX_W
) -> _Item:
    """A sheet as one block of strands at `pitch`, each pointing as it runs relative to strands[0]; a closed
    barrel's first strand repeated as a ghost after the last. `over` (helices packed on its face, in sequence): each
    drawn over the strands where it packs (_over_sheet), unless it would lie on one drawn there before it."""
    strands = sheet.strands
    order, s, cents = _sheet_order(sheet, frame)
    members, offsets = [], []
    column = sheet.columns or list(range(len(strands)))  # strands sharing a column (sheets._columns) stand in it
    slot_of = {c: n for n, c in enumerate(dict.fromkeys(column[k] for k in order))}
    slots = len(slot_of) + (1 if sheet.closed else 0)
    heights = _column_heights(sheet, order, s) if sheet.columns else {}
    for k in order:
        up = sheet.directions[k] * s > 0
        sse = strands[k]
        members.append(Placed(sse, 0, 0, _strand_length(sse), np.pi / 2 if up else -np.pi / 2, pitch))
        offsets.append(((slot_of[column[k]] - (slots - 1) / 2) * pitch, heights.get(k, 0.0)))
    if sheet.closed:  # the barrel wraps: its first strand pairs again with the last
        f = members[0]
        members.append(Placed(f.sse, 0, 0, f.length, f.angle, f.width, ghost=True))
        offsets.append(((slots - 1) / 2 * pitch, 0.0))
    laid: list[tuple[float, float, float, float]] = []
    for h in over:  # in sequence; one that would lie on a helix already over the sheet stays beside it
        member, (ox, oy) = _over_sheet(
            h, sheet, order, s, [o[0] for o in offsets[: len(order)]], members[: len(order)], helix_w
        )
        x0, y0, x1, y1 = member.rect
        box = (x0 + ox, y0 + oy, x1 + ox, y1 + oy)
        if overlaps_any(box, laid):
            continue
        laid.append(box)
        members.append(member)
        offsets.append((ox, oy))
    home = tuple(cents.mean(axis=0) * SCALE)
    return _Item(members, offsets, home[0], home[1], home, *_anchor(strands))


def _unit(v) -> np.ndarray:
    return _snap(np.asarray(v, float) / np.linalg.norm(v))


def _over_sheet(
    h: SSE, sheet: Sheet, order: list[int], s: int, xs: list[float], strands: list[Placed], width: float = HELIX_W
) -> tuple[Placed, tuple[float, float]]:
    """A helix packed on the sheet's face (features.helices_on_sheets), placed over the drawn strands as it packs on
    them, true length: the member and its offset from the sheet item's centre. The sheet's drawn axes in 3D: up3,
    the strands' mean direction (each by its sense) times the item's flip s, and across3, from the first drawn
    strand's centroid to the last's, made orthogonal to up3; n = across3 x up3, the face as drawn.
    x: between the two strands whose axis lines lie nearest the helix centroid in the plane of the face (both seen
    along n), by inverse distance, from their drawn x (`xs`, in drawn order); beyond an end strand, that strand's x.
    y: the centroid's height along up3 above the strands' mean centroid, scaled, within the strands' drawn extent.
    Angle: atan2(axis . up3, axis . across3), not snapped to 90 degrees (Burroughs 2007 draws ubiquitin's α1 diagonal).
    Every 3D vector and the results are snapped."""
    cents = [_snap(sheet.strands[k].centroid) for k in order]
    up3 = _unit(s * sum(d * x.axis for d, x in zip(sheet.directions, sheet.strands)))
    across = cents[-1] - cents[0]
    across3 = _unit(across - (across @ up3) * up3)
    n = _snap(np.cross(across3, up3))
    centre, axis = _snap(h.centroid), _snap(h.axis)
    mid = _snap(np.mean(cents, axis=0))

    def flat(v):  # seen along n: in the plane of the face
        return v - (v @ n) * n

    dist = []
    for k, c in zip(order, cents):
        a = _unit(flat(sheet.strands[k].axis))
        v = flat(centre - c)
        dist.append(float(_snap(np.linalg.norm(v - (v @ a) * a))))
    t = float(_snap((centre - cents[0]) @ across3))
    end = float(_snap((cents[-1] - cents[0]) @ across3))
    ranked = sorted(range(len(order)), key=lambda i: (dist[i], i))
    i, j = ranked[0], ranked[-1] if len(ranked) == 1 else ranked[1]
    if (i == 0 and t <= 0) or (i == len(order) - 1 and t >= end) or dist[i] == 0:
        x = xs[i]
    else:
        x = (xs[i] / dist[i] + xs[j] / dist[j]) / (1 / dist[i] + 1 / dist[j])
    reach = max(m.length for m in strands if not m.ghost) / 2
    y = float(np.clip((centre - mid) @ up3 * SCALE, -reach, reach))
    angle = float(_snap(np.arctan2(axis @ up3, axis @ across3)))
    length = max(width, (len(h) - 1) * HELIX_RISE * SCALE)
    return Placed(h, 0, 0, length, angle, width, over=True), (float(_snap(x)), float(_snap(y)))


_BA_GAP = 0.3  # page units between the footprints of neighbours in a βα template row
_BA_MAX = 10  # strands: a longer parallel sheet (2BNH's 17-strand leucine-rich repeat) is no α/β barrel or Rossmann
# sheet, and as one row it would run across the page
_BA_SHARE = (3, 5)  # the template applies when at least this share of a sheet's βα units are in its row


class BAUnit(NamedTuple):
    """Two strands of one chain consecutive in sequence among a sheet's strands, and the helices between them."""

    p: SSE
    q: SSE
    slots: tuple[int, int]  # their places in sheet.strands; in a barrel's closing unit, the first's is len(strands)
    helices: list[SSE]  # the H elements of the chain between p.end and q.start, in sequence
    inrow: bool  # the strands are neighbours in the sheet (or close the barrel) and a helix joins them


def _ba_units(sheet: Sheet, sses: list[SSE]) -> list[BAUnit]:
    """The βα units of a sheet, chain by chain in sequence. A closed sheet's ghost counts as slot len(strands): the
    unit joining its first and last slots (the closing unit) is in its row, between the last strand and the
    ghost."""
    slot = {s.id: k for k, s in enumerate(sheet.strands)}
    n = len(sheet.strands)
    out = []
    for chain in sorted({s.chain for s in sheet.strands}):
        mine = sorted((s for s in sheet.strands if s.chain == chain), key=lambda s: s.start)
        for p, q in zip(mine, mine[1:]):
            sp, sq = slot[p.id], slot[q.id]
            if sheet.closed and n > 2 and {sp, sq} == {0, n - 1}:
                sp, sq = (n if sp == 0 else sp), (n if sq == 0 else sq)
            hs = [s for s in sses if s.chain == chain and s.kind == "H" and p.end < s.start and s.end < q.start]
            out.append(BAUnit(p, q, (sp, sq), hs, abs(sp - sq) == 1 and bool(hs)))
    return out


def _open_at_ends(sheet: Sheet) -> Sheet:
    """A barrel of one chain opened between its last and first strand in sequence (β8|β1 of a TIM barrel, the cut of
    the textbook figures), when those two pair: the same ring of strands, turned so the row (_ba_template) reads in
    sequence and no βα unit closes the barrel round the row. Anything else is returned as it is."""
    n = len(sheet.strands)
    if not sheet.closed or n < 3 or len({s.chain for s in sheet.strands}) != 1:
        return sheet
    by_seq = sorted(range(n), key=lambda k: sheet.strands[k].start)
    lo, hi = sorted((by_seq[0], by_seq[-1]))
    if hi - lo != 1:  # already the ends (0, n - 1), or not neighbours in the ring
        return sheet
    order = list(range(hi, n)) + list(range(hi))
    kinds = {frozenset((k, k + 1)): kind for k, kind in enumerate(sheet.pair_kinds)}
    d = sheet.directions
    return Sheet(
        [sheet.strands[k] for k in order],
        [d[k] * d[order[0]] for k in order],
        [kinds.get(frozenset((a, b)), "P" if d[a] == d[b] else "A") for a, b in zip(order, order[1:])],
        True,
        sheet.ambiguous,
    )


def _row_shape(sheet: Sheet) -> bool:
    """Can the sheet be a βα row at all: all parallel, four to _BA_MAX strands, its order not a guess?"""
    return 4 <= len(sheet.strands) <= _BA_MAX and not sheet.ambiguous and all(k == "P" for k in sheet.pair_kinds)


def template_applies(sheet: Sheet, sses: list[SSE]) -> bool:
    """Is the sheet, on its own, drawn as a βα row (_ba_template)? An all-parallel sheet of 4 to _BA_MAX strands, its
    order not a guess, with at least 3/5 of its βα units in the row (1TIM: 7 of 7, 1LDM: 4 of 5). Mixed and
    antiparallel sheets never are. Pure: it changes nothing. Copies of a protein share the decision: see
    template_choices."""
    if not _row_shape(sheet):
        return False
    units = _ba_units(sheet, sses)
    return bool(units) and _BA_SHARE[1] * sum(u.inrow for u in units) >= _BA_SHARE[0] * len(units)


def template_choices(sheets: list[Sheet], sses: list[SSE], copies: dict[str, tuple[str, int]]) -> list[bool]:
    """Which sheets are drawn as βα rows, copies of one protein alike. `copies` (model.homolog_chains): chain ->
    (first chain of its copies, residue-number offset). Sheets of one chain each match when they belong to copies of
    one protein and at least half the strands of the smaller pair one-to-one with the other sheet's (model.pair_places,
    residue numbers taken to the first copy's).
    In each group of matching sheets: if any takes the template on its own (template_applies) and every one can be a
    row (_row_shape), all of them take it (a copy missing a helix or two is drawn as its mates are); if any cannot
    be a row, none takes it, so the copies stay alike. Pure."""
    own = [template_applies(sh, sses) for sh in sheets]

    def mapped(sh):
        chains = {s.chain for s in sh.strands}
        if len(chains) != 1:
            return None
        (c,) = chains
        root, shift = copies.get(c, (c, 0))
        return root, c, [("E", s.first.seq - shift, s.last.seq - shift) for s in sh.strands]

    keys = [mapped(sh) for sh in sheets]
    links = []
    for i, j in itertools.combinations(range(len(sheets)), 2):
        a, b = keys[i], keys[j]
        if a is None or b is None or a[0] != b[0] or a[1] == b[1]:
            continue
        if 2 * len(pair_places(a[2], b[2])) >= min(len(a[2]), len(b[2])):
            links.append((i, j))
    group = _components(len(sheets), links)
    out = []
    for i in range(len(sheets)):
        mine = [k for k in range(len(sheets)) if group[k] == group[i]]
        out.append(any(own[k] for k in mine) and all(_row_shape(sheets[k]) for k in mine))
    return out


def template_rows(sheets: list[Sheet], sses: list[SSE], copies: dict[str, tuple[str, int]]) -> list[list[str]]:
    """The elements of each sheet drawn as a βα row (template_choices), in projected mode with fold_template auto:
    its strands and every helix of its in-row βα units (the one in the row and any drawn free beside it). Pure."""
    rows = []
    for sh, row in zip(sheets, template_choices(sheets, sses, copies)):
        if row:
            units = _ba_units(_open_at_ends(sh), sses)
            helices = sorted({h.id for u in units if u.inrow for h in u.helices})
            rows.append([s.id for s in sh.strands] + helices)
    return rows


def _ba_template(
    sheet: Sheet, units: list[BAUnit], frame: Frame, pitch: float = PITCH, helix_w: float = HELIX_W
) -> _Item:
    """An all-parallel α/β sheet (template_applies) as one row, as the hand-drawn TIM barrel and Rossmann figures
    draw it: the strands in sheet order and pointing as _sheet_item draws them (the front face), and after each
    strand, when an in-row unit joins it to the next, that unit's longest helix (the first on a tie), true length,
    running against the strands; a barrel's closing helix between the last strand and the ghost. Each footprint
    (labels included) starts _BA_GAP after the last one ends, all centred on one line, the row centred on the item.
    Other helices are placed as free helices. Home: the strands' projected mean, as _sheet_item."""
    strands = sheet.strands
    order, s, cents = _sheet_order(sheet, frame)
    if strands[order[0]].start > strands[order[-1]].start:  # the row reads left to right in sequence: turned
        order, s = order[::-1], -s  # half round in the page (never mirrored), its strands then point down (K9: under
        # --up 1TIM read β8 ... β1; 1GOX's barrels and two of 1LDM's Rossmann rows did by default)
    n = len(order)
    drawn = {k: d for d, k in enumerate(order)}
    after: dict[int, SSE] = {}  # drawn slot -> the helix drawn right of it
    for u in units:
        if u.inrow:
            gap = n - 1 if n in u.slots else min(drawn[u.slots[0]], drawn[u.slots[1]])
            after[gap] = max(u.helices, key=lambda h: (len(h), -h.start))
    members: list[Placed] = []
    for d, k in enumerate(order):
        up = sheet.directions[k] * s > 0
        sse = strands[k]
        members.append(Placed(sse, 0, 0, _strand_length(sse), np.pi / 2 if up else -np.pi / 2, pitch))
        if d in after:
            h = after[d]
            length = max(helix_w, (len(h) - 1) * HELIX_RISE * SCALE)
            members.append(Placed(h, 0, 0, length, -np.pi / 2 if up else np.pi / 2, helix_w))
    if sheet.closed:
        f = members[0]
        members.append(Placed(f.sse, 0, 0, f.length, f.angle, f.width, ghost=True))
    xs, right = [], None
    for m in members:
        fp = _footprint(m)
        x = -fp[0] if right is None else right + _BA_GAP - fp[0]
        xs.append(x)
        right = x + fp[2]
    mid = (_footprint(members[0])[0] + right) / 2
    offsets = [(float(_snap(x - mid)), 0.0) for x in xs]
    home = tuple(cents.mean(axis=0) * SCALE)
    return _Item(members, offsets, home[0], home[1], home, *_anchor(strands))


def _helix_item(
    sse: SSE, frame: Frame, angle_mode: str = "snap", width: float = HELIX_W, rise: float = HELIX_RISE
) -> _Item:
    """Always true length. angle_mode: 'tilted' as seen in the view; 'snap' to the nearest of up, down, left,
    right (the topology-diagram convention); 'upright' up or down like a strand (stack mode)."""
    au, av = frame.direction(sse.centroid, sse.axis)
    length = max(width, (len(sse) - 1) * rise * SCALE)  # true length even when tilted: equal pitch everywhere
    angle = _angle(float(np.arctan2(av, au)), float(av), angle_mode)
    pos = frame.project(sse.centroid, [sse.chain])[0] * SCALE
    placed = Placed(sse, 0, 0, length, angle, width)
    return _Item(
        [placed], [(0.0, 0.0)], float(pos[0]), float(pos[1]), (float(pos[0]), float(pos[1])), sse.start, sse.chain
    )


_KINK_CORNER = 5  # residues: the longest corner between two pieces of a kinked helix (a π turn)


def _kink_groups(sses: list[SSE], skip: set[str]) -> list[list[SSE]]:
    """The pieces of each helix split at a kink (SSE.kinked: build_sses cut it at a 3-10/pi corner), in order: each
    kinked piece with the helix piece of its chain that ends at most _KINK_CORNER residues before it (a subset of the
    elements may leave that piece out). Groups any of whose pieces is in `skip` (a bundle holds them) are left out."""
    groups: list[list[SSE]] = []
    last: dict[str, SSE] = {}  # each chain's last helix so far
    for b in sorted(sses, key=lambda s: s.start):
        a = last.get(b.chain)
        if b.kind == "H" and b.kinked and a is not None and 0 <= b.start - a.end - 1 <= _KINK_CORNER:
            if groups and groups[-1][-1] is a:
                groups[-1].append(b)
            else:
                groups.append([a, b])
        if b.kind == "H":
            last[b.chain] = b
    return [g for g in groups if not any(s.id in skip for s in g)]


_KINK_ROOM = 0.25  # one routing cell (route.RES): the corner loop's room on each leg between two kinked pieces


def _kink_item(group: list[SSE], frame: Frame, angle_mode: str, width: float = HELIX_W) -> _Item:
    """The pieces of a kinked helix as one unit, end to end, so the corner reads as one helix turning (2HHB α's AB
    corner, which an element seen between the pieces kept five widths apart): each piece as _helix_item draws it
    (true length, its own angle). Each next piece starts where the loop from the last one's C port turns once into
    it: a piece across the last starts just beyond its end and just off its axis, both legs a routing cell clear of
    the two bodies (_KINK_ROOM); a piece running on starts a little further along; a piece turned back runs beside
    the last, MARGIN apart, on the side it is seen. The unit sits at the mean of the pieces' projected places,
    offsets on the grid."""
    parts = [_helix_item(s, frame, angle_mode, width) for s in group]
    members = [p.members[0] for p in parts]
    at = [np.zeros(2)]
    for k, (a, b) in enumerate(zip(members, members[1:])):
        end = at[-1] + a.direction * a.length / 2
        turn = float(a.direction @ b.direction)
        normal = np.array([-a.direction[1], a.direction[0]])
        if turn < -0.5:  # turned back: a hairpin beside the last piece, on the side the piece turns to
            out = float(_snap(b.direction @ normal))
            if out == 0.0:  # straight back (snapped angles): the side it is seen on
                out = float(np.subtract((parts[k + 1].cx, parts[k + 1].cy), (parts[k].cx, parts[k].cy)) @ normal)
            start = end + normal * (1 if out >= 0 else -1) * ((a.width + b.width) / 2 + MARGIN)
        elif turn > 0.5:  # running on
            start = end + a.direction * 4 * _KINK_ROOM
        else:  # across: one turn of the loop between the two
            start = end + a.direction * (b.width / 2 + _KINK_ROOM) + b.direction * (a.width / 2 + _KINK_ROOM)
            # a piece turned more than 90 degrees (drawn tilted) leans back over the last one's end: its near
            # corners kept _KINK_ROOM beyond that end
            side = np.array([-b.direction[1], b.direction[0]]) * b.width / 2
            back = max(float((end - c) @ a.direction) for c in (start + side, start - side)) + _KINK_ROOM
            if back > 1e-9:
                start = start + a.direction * back
        at.append(start + b.direction * b.length / 2)
    mean = np.mean(at, axis=0)
    offsets = [tuple(float(v) for v in _snap(x - mean)) for x in at]
    x, y = (float(v) for v in _snap(np.mean([(p.cx, p.cy) for p in parts], axis=0)))
    return _Item(members, offsets, x, y, (x, y), group[0].start, group[0].chain)


def _angle(seen: float, av: float, mode: str) -> float:
    if mode == "upright":
        return np.pi / 2 if av >= 0 else -np.pi / 2
    if mode == "snap":
        return float(np.round(seen / (np.pi / 2)) * (np.pi / 2))
    return float(_snap(seen))  # as seen, on the grid (radians)


def _bundle_item(group: list[SSE], frame: Frame, angle_mode: str, width: float) -> _Item:
    """Helices of several chains packed along their length, drawn as one unit: one column per helix (kinked
    pieces of one helix share a column, end to end), columns side by side in their order across the bundle at
    their true spacing (at least a label's width apart), staggered along the shared axis as in 3D;
    antiparallel helices point the other way."""
    from .features import KINK_GAP

    ref = group[0].axis
    mean = sum(np.sign(float(s.axis @ ref)) * s.axis * len(s) for s in group)
    m = mean / np.linalg.norm(mean)
    centre = np.mean([s.centroid for s in group], axis=0)
    a2 = frame.direction(centre, m)
    angle = _angle(float(np.arctan2(a2[1], a2[0])), float(a2[1]), angle_mode)
    d2, side2 = np.array([np.cos(angle), np.sin(angle)]), np.array([-np.sin(angle), np.cos(angle)])
    columns: list[list[SSE]] = []
    for s in sorted(group, key=lambda s: (s.chain, s.start)):
        last = columns[-1][-1] if columns else None
        if last is not None and last.chain == s.chain and s.start - last.end - 1 <= KINK_GAP:
            columns[-1].append(s)  # the next piece of a kinked helix
        else:
            columns.append([s])
    lateral = [float(np.mean([frame.project(s.centroid, [s.chain])[0] @ side2 for s in col])) for col in columns]
    columns = [c for _, c in sorted(zip(lateral, columns), key=lambda x: x[0])]
    pts = [np.mean([s.centroid for s in col], axis=0) for col in columns]
    gaps = [np.linalg.norm((a - b) - ((a - b) @ m) * m) for a, b in zip(pts, pts[1:])]
    spacing = max(width + 1.2, float(np.mean(gaps)) * SCALE if gaps else 0.0)
    members, offsets = [], []
    for k, col in enumerate(columns):
        last_end = -np.inf
        for s in sorted(col, key=lambda s: float((s.centroid - centre) @ m)):
            length = max(width, (len(s) - 1) * HELIX_RISE * SCALE)
            along = max(float((s.centroid - centre) @ m) * SCALE, last_end + 0.3 + length / 2)  # kinked pieces in a row
            last_end = along + length / 2
            members.append(Placed(s, 0, 0, length, angle if float(s.axis @ m) >= 0 else angle + np.pi, width))
            lat = (k - (len(columns) - 1) / 2) * spacing
            offsets.append(tuple(float(v) for v in _snap(d2 * along + side2 * lat)))  # on the grid, like positions
    pos = frame.project(centre)[0] * SCALE
    first = min(group, key=lambda s: s.start)
    return _Item(
        members, offsets, float(pos[0]), float(pos[1]), (float(pos[0]), float(pos[1])), first.start, first.chain
    )


_DENSE_PAIRS = 256  # _overlapping tests every pair of up to this many footprints at once; more, only near ones


def _overlapping(xy: np.ndarray, halves: np.ndarray, mids: np.ndarray, apart: np.ndarray | None = None) -> np.ndarray:
    """(i, j) with i < j, in row order, for every pair of footprints (plus MARGIN) that overlap; pairs marked in
    the boolean matrix `apart` (_joint_matrix) are never counted. Past _DENSE_PAIRS footprints only the pairs whose
    centres lie near enough across the page to overlap are measured (random-100 C33: every sweep measured all
    pairs, 20 s of 7JU4's layout), by the same arithmetic, so the pairs are the same."""
    q = xy + mids
    n = len(q)
    if n <= _DENSE_PAIRS:
        depth = halves[:, None, :] + halves[None, :, :] + MARGIN - np.abs(q[None, :, :] - q[:, None, :])
        hit = (depth[..., 0] > 1e-9) & (depth[..., 1] > 1e-9)
        if apart is not None:
            hit &= ~apart
        return np.argwhere(np.triu(hit, 1))
    order = np.argsort(q[:, 0], kind="stable")
    sx = q[order, 0]
    reach = (halves[order, 0] + halves[:, 0].max() + MARGIN) * (1 + 1e-9) + 1e-6  # past any pair that can overlap
    stop = np.searchsorted(sx, sx + reach, side="right")
    start = np.arange(n) + 1
    count = np.maximum(stop - start, 0)
    first = np.repeat(order, count)
    other = order[np.repeat(start - np.cumsum(count) + count, count) + np.arange(count.sum())]
    a, b = np.minimum(first, other), np.maximum(first, other)
    dx = halves[a, 0] + halves[b, 0] + MARGIN - np.abs(q[b, 0] - q[a, 0])  # as the dense matrix's [a, b]
    dy = halves[a, 1] + halves[b, 1] + MARGIN - np.abs(q[b, 1] - q[a, 1])
    hit = (dx > 1e-9) & (dy > 1e-9)
    if apart is not None:
        hit &= ~apart[a, b]
    a, b = a[hit], b[hit]
    rows = np.lexsort((b, a))
    return np.stack([a[rows], b[rows]], axis=1).reshape(-1, 2)


def _joint_matrix(n: int, joints) -> np.ndarray:
    """Boolean n x n matrix, True for both orders of every item pair in `joints` (_joints)."""
    m = np.zeros((n, n), bool)
    for i, j in joints:
        m[i, j] = m[j, i] = True
    return m


class Walls(NamedTuple):
    """Boxes that stay where they are (x0, y0, x1, y1 rows: each nucleotide of an unpaired strand's trace, a bead
    either side), and for each (item, box) whether the item must keep off it: an item holding an element that
    touches that trace is exempt (a single strand lies on the protein it wraps)."""

    boxes: np.ndarray  # (m, 4)
    hard: np.ndarray  # (items, m) bool


def _separate(
    items: list[_Item],
    ties: dict[int, tuple[int, int]] | None = None,
    joints: set[tuple[int, int]] = frozenset(),
    fixed: Walls | None = None,
    mirror: float | None = None,
    room: float = _SEP_SYM_ROOM,
    level: np.ndarray | None = None,
) -> None:
    """Push overlapping items apart, then slide any left over clear. With `ties` (_symmetry_ties) the pushes
    keep every protomer drawn alike (_push_apart), and that result is kept if it clears every overlap within
    _SEP_ORBIT_SWEEPS sweeps and its box takes at most `room` times the area of moving the items one by
    one from the same start (widening the spread can otherwise cost far more room); else the one-by-one
    result is used. Consecutive pieces of one bent duplex (`joints`, _joints) in a short chain (_rigid_joints)
    are moved together and never pushed apart. Items are pushed off the `fixed` boxes they must keep off (Walls),
    which stay where they are. `mirror`: the copies are their mates' mirror images about that page x (_push_apart);
    then the tied pushes run as many passes as the one-by-one ones and may take up to _SEP_MIRROR_ROOM. `level`
    (bool per item): items that keep their height (the transmembrane elements of a membrane view), only ever pushed
    and slid sideways (_push_apart, _slide_clear)."""
    halves = np.array([it.half for it in items], float).reshape(-1, 2)
    mids = np.array([it.mid for it in items], float).reshape(-1, 2)
    start = np.array([(it.cx, it.cy) for it in items], float).reshape(-1, 2)
    home = np.array([it.home for it in items], float).reshape(-1, 2)
    joints = _rigid_joints(items, joints)
    apart = _joint_matrix(len(items), joints)
    xy = start.copy()
    alone = start.copy()
    _push_apart(alone, halves, mids, home, {}, _SEP_SWEEPS, apart, fixed, level=level)
    if mirror is not None:  # mirror images about a dyad in the page: drawn so whenever room allows (9GEP)
        sweeps, room = _SEP_SWEEPS, max(room, _SEP_MIRROR_ROOM)
    else:
        sweeps = _SEP_ORBIT_SWEEPS
    together = ties and _push_apart(xy, halves, mids, home, ties, sweeps, apart, fixed, mirror, level)
    if not together or _box_area(xy, halves, mids) > room * _box_area(alone, halves, mids):
        xy = alone
    for it, (x, y) in zip(items, xy.tolist()):
        it.cx, it.cy = x, y
    _slide_clear(items, joints, fixed, level)


def _components(n: int, links) -> np.ndarray:
    """Component id (its least member) of each of n items, linked pairwise by `links` (union-find)."""
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for i, j in links:
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)
    return np.array([find(i) for i in range(n)], int)


def _body_overlap(ca, a: PlacedStem, cb, b: PlacedStem) -> float:
    """How deep the bodies (length x width rectangles turned by their angles) of stems a and b, centred at ca and
    cb, overlap: the least penetration over the four separating axes; 0 when they do not."""
    d = np.subtract(cb, ca)
    depth = np.inf
    for axis in (a.direction, a.normal, b.direction, b.normal):
        ra = abs(axis @ a.direction) * a.length / 2 + abs(axis @ a.normal) * a.width / 2
        rb = abs(axis @ b.direction) * b.length / 2 + abs(axis @ b.normal) * b.width / 2
        depth = min(depth, ra + rb - abs(axis @ d))
    return max(float(depth), 0.0)


def body_overlap(a: PlacedStem, b: PlacedStem) -> float:
    """How deep two placed stems' bodies overlap (_body_overlap), page units."""
    return _body_overlap((a.cx, a.cy), a, (b.cx, b.cy), b)


_JOINT_GROUP = 4  # at most this many consecutive pieces of a duplex are separated as one group
_JOINT_OVERLAP = 0.25  # widths two joined pieces' bodies may overlap and still be moved together


def _rigid_joints(items: list[_Item], joints: set[tuple[int, int]]) -> set[tuple[int, int]]:
    """The joints (_joints) whose pieces _separate moves together: those in chains of at most _JOINT_GROUP pieces (a
    longer chain moved as one, as the nucleosome's 22 pieces, carries the protein it wraps out of its ring), and of
    those only pairs whose bodies overlap by at most _JOINT_OVERLAP widths. The other joined pieces are pushed
    apart like any two items."""
    if not joints:
        return set()
    cid = _components(len(items), sorted(joints))
    size = np.bincount(cid, minlength=len(items))
    rigid = set()
    for i, j in sorted(joints):
        if size[cid[i]] > _JOINT_GROUP:
            continue
        a, b = items[i].members[0], items[j].members[0]
        if _body_overlap((items[i].cx, items[i].cy), a, (items[j].cx, items[j].cy), b) <= _JOINT_OVERLAP * max(
            a.width, b.width
        ):
            rigid.add((i, j))
    return rigid


def _box_area(xy: np.ndarray, halves: np.ndarray, mids: np.ndarray) -> float:
    """Area of the box around every footprint."""
    q = xy + mids
    lo, hi = (q - halves).min(axis=0), (q + halves).max(axis=0)
    return float((hi[0] - lo[0]) * (hi[1] - lo[1]))


def _push_apart(
    xy: np.ndarray,
    halves: np.ndarray,
    mids: np.ndarray,
    home: np.ndarray,
    ties: dict[int, tuple[int, int]],
    sweeps: int,
    apart: np.ndarray | None = None,
    fixed: Walls | None = None,
    mirror: float | None = None,
    level: np.ndarray | None = None,
) -> bool:
    """Push overlapping footprints apart in place (axis of least penetration); spring back to home early on.
    Each sweep finds the overlapping pairs at once and resolves them in order. With `ties` a push moves an
    item's whole orbit (it and its symmetry copies) by the same amount, and a sideways push between two
    protomers moves protomer k by k times the amount when that widens the spread between protomers (also
    between two copies of one element); otherwise two copies of one element, which no orbit-wide move can
    part, move alone. Items joined by `apart` (consecutive pieces of one bent duplex) are never pushed apart: a push
    moves the whole chain of pieces (with their orbits); two items of one such group that overlap move alone.
    An item overlapping a `fixed` box it must keep off (Walls) is pushed (with its group) the whole way off it.
    With `mirror` (page x of a dyad in the page) each copy is its mate's mirror image: a sideways push moves the
    copies of a group the other way, an item and its own copy are pushed apart sideways only, symmetrically about
    the mirror line, and the spread is never widened. `level` (bool per item): an item (with its group) that keeps
    its height; a push up or down moves only the other one of a pair, and two such items are pushed sideways. True
    if no overlap is left."""
    n = len(xy)
    root = list(range(n))  # the protomer-0 item of each item's orbit
    copy_k = np.zeros(n)  # protomer of each copy (0 for everything else)
    for i, (j, k) in ties.items():
        root[i], copy_k[i] = j, k
    orbits: dict[int, list[int]] = {}
    for i, r in enumerate(root):
        orbits.setdefault(r, []).append(i)
    joined = np.argwhere(np.triu(apart, 1)).tolist() if apart is not None else []
    cid = _components(n, [(i, r) for i, r in enumerate(root)] + joined)  # orbits and joined pieces move as one
    members = {c: np.flatnonzero(cid == c) for c in sorted(set(cid.tolist()))}
    group = [members[c] for c in cid.tolist()]
    protomer = [copy_k[i] if len(orbits[root[i]]) > 1 else None for i in range(n)]  # None: item has no copies
    widen = np.sign(sum((xy[i, 0] - xy[j, 0]) / k for i, (j, k) in ties.items())) or 1.0  # protomer k lies this way
    sx = np.where((copy_k > 0) & (mirror is not None), -1.0, 1.0)  # a mirrored copy moves against its mate in x
    walls = fixed
    held = None if level is None else np.array([bool(level[g].any()) for g in group])  # its group keeps its height
    for sweep in range(sweeps):
        pairs = _overlapping(xy, halves, mids, apart)
        if not len(pairs) and not _off_walls(xy, halves, mids, walls, group, move=False) and sweep >= 100:
            return True
        for i, j in pairs.tolist():
            dx = xy[j, 0] + mids[j, 0] - xy[i, 0] - mids[i, 0]
            dy = xy[j, 1] + mids[j, 1] - xy[i, 1] - mids[i, 1]
            px = halves[i, 0] + halves[j, 0] + MARGIN - abs(dx)
            py = halves[i, 1] + halves[j, 1] + MARGIN - abs(dy)
            if px <= 1e-9 or py <= 1e-9:  # an earlier push this sweep already cleared it
                continue
            gi, gj = (group[i], group[j]) if cid[i] != cid[j] else (i, j)
            if mirror is not None and cid[i] == cid[j] and root[i] == root[j]:  # an item and its own mirror copy
                py = np.inf
            if held is not None and held[i] and held[j]:  # both keep their height: sideways
                py = np.inf
            if px < py:
                s = (1 if dx > 0 else -1) if dx != 0 else 1
                if mirror is None and None not in (protomer[i], protomer[j]) and protomer[i] != protomer[j]:
                    wider = s * (px + 2 * _SEP_PUSH) / (protomer[j] - protomer[i])
                    if wider * widen > 0:  # items without copies (copy_k 0) stay put, whichever protomer
                        xy[:, 0] += copy_k * wider
                        continue
                xy[gi, 0] -= sx[gi] * sx[i] * s * (px / 2 + _SEP_PUSH)
                xy[gj, 0] += sx[gj] * sx[j] * s * (px / 2 + _SEP_PUSH)
            elif held is not None and (held[i] or held[j]):  # the other one moves the whole way
                s = (1 if dy > 0 else -1) if dy != 0 else 1
                if held[i]:
                    xy[gj, 1] += s * (py + _SEP_PUSH)
                else:
                    xy[gi, 1] -= s * (py + _SEP_PUSH)
            else:
                s = (1 if dy > 0 else -1) if dy != 0 else 1
                xy[gi, 1] -= s * (py / 2 + _SEP_PUSH)
                xy[gj, 1] += s * (py / 2 + _SEP_PUSH)
        _off_walls(xy, halves, mids, walls, group, sx=sx, held=held)
        if sweep < 100:
            xy += 0.03 * (home - xy)
    return not len(_overlapping(xy, halves, mids, apart)) and not _off_walls(xy, halves, mids, walls, group, move=False)


def _off_walls(
    xy, halves, mids, walls: Walls | None, group, move: bool = True, sx: np.ndarray | None = None, held=None
) -> bool:
    """Push each footprint overlapping a fixed box it must keep off (Walls), MARGIN kept, the whole way off it
    along the axis of least overlap, moving its group (`group[i]`, the items moved with item i) in place; `sx`
    (-1 for a mirrored copy, _push_apart) turns a sideways move of a group's copies the other way; `held` (bool per
    item, _push_apart's): a group keeping its height is moved sideways only. True if any overlapped (with
    move=False, only that is reported)."""
    if walls is None or not len(walls.boxes):
        return False
    boxes = np.asarray(walls.boxes, float)
    centres, half = (boxes[:, :2] + boxes[:, 2:]) / 2, (boxes[:, 2:] - boxes[:, :2]) / 2
    d = centres[None, :, :] - (xy + mids)[:, None, :]
    pen = halves[:, None, :] + half[None, :, :] + MARGIN - np.abs(d)
    hits = np.argwhere((pen[..., 0] > 1e-9) & (pen[..., 1] > 1e-9) & walls.hard)
    if not move or not len(hits):
        return bool(len(hits))
    for i, w in hits.tolist():  # in order; an earlier push may already have cleared it
        dx, dy = centres[w] - xy[i] - mids[i]
        px = halves[i, 0] + half[w, 0] + MARGIN - abs(dx)
        py = halves[i, 1] + half[w, 1] + MARGIN - abs(dy)
        if px <= 1e-9 or py <= 1e-9:
            continue
        if px < py or (held is not None and held[i]):
            turn = 1.0 if sx is None else sx[group[i]] * sx[i]
            xy[group[i], 0] -= turn * ((1 if dx > 0 else -1) if dx != 0 else 1) * (px + _SEP_PUSH)
        else:
            xy[group[i], 1] -= ((1 if dy > 0 else -1) if dy != 0 else 1) * (py + _SEP_PUSH)
    return True


def _slide_clear(
    items: list[_Item],
    joints: set[tuple[int, int]] = frozenset(),
    fixed: Walls | None = None,
    level: np.ndarray | None = None,
) -> None:
    """Guarantee: each item in turn, if it still overlaps one already settled (other than a piece it is joined
    to, `joints`) or a `fixed` box it must keep off, takes the nearest clear spot straight right, left, up or down
    (in that order on a tie), in _SLIDE_STEP steps; one that keeps its height (`level`) right or left only."""
    halves = [it.half for it in items]
    mids = [it.mid for it in items]
    boxes = [] if fixed is None else np.asarray(fixed.boxes, float).tolist()
    walls = [
        (SimpleNamespace(cx=(x0 + x1) / 2, cy=(y0 + y1) / 2), ((x1 - x0) / 2, (y1 - y0) / 2), (0.0, 0.0))
        for x0, y0, x1, y1 in boxes
    ]
    for i, it in enumerate(items):
        x0, y0 = it.cx, it.cy
        mine = [w for k, w in enumerate(walls) if fixed.hard[i, k]] if walls else []
        settled = mine + [(items[j], halves[j], mids[j]) for j in range(i) if (j, i) not in joints]
        k = 0
        ways = ((1, 0), (-1, 0)) if level is not None and level[i] else ((1, 0), (-1, 0), (0, 1), (0, -1))
        while any(_collide(it, halves[i], mids[i], *other) for other in settled):
            k += 1
            for dx, dy in ways:
                it.cx, it.cy = x0 + k * dx * _SLIDE_STEP, y0 + k * dy * _SLIDE_STEP
                if not any(_collide(it, halves[i], mids[i], *other) for other in settled):
                    break


def _owner(items: list[_Item]) -> dict[str, tuple[int, Placed, tuple[float, float]]]:
    return {
        (m.id if isinstance(m, PlacedStem) else m.sse.id): (k, m, o)
        for k, it in enumerate(items)
        for m, o in zip(it.members, it.offsets)
        if not m.ghost
    }


def _contact_weights(items: list[_Item], contacts: dict[frozenset[str], int]) -> np.ndarray:
    """Symmetric pull strength between items whose elements touch in 3D."""
    owner = _owner(items)
    w = np.zeros((len(items), len(items)))
    for pair, count in sorted((tuple(sorted(k)), v) for k, v in contacts.items()):
        if len(pair) != 2:
            continue  # an element touching itself carries no pull
        a, b = (owner.get(x) for x in pair)
        if a and b and a[0] != b[0]:
            w[a[0], b[0]] += _W_CONTACT * min(count, _CONTACT_CAP) / _CONTACT_CAP
            w[b[0], a[0]] = w[a[0], b[0]]
    return w


def _loop_springs(items: list[_Item], sses: list[SSE]) -> list[tuple[int, np.ndarray, int, np.ndarray]]:
    """(item of a, C-port offset of a, item of b, N-port offset of b) for each loop between items."""
    owner = _owner(items)
    springs = []
    for a, b in zip(sses, sses[1:]):
        if a.chain != b.chain or a.id not in owner or b.id not in owner:
            continue
        (i, ma, oa), (j, mb, ob) = owner[a.id], owner[b.id]
        if i != j:
            springs.append((i, np.add(oa, ma.direction * ma.length / 2), j, np.add(ob, -mb.direction * mb.length / 2)))
    return springs


def _mass(it: _Item) -> float:
    """How firmly an item holds its projected place: one per element; a duplex like a sheet of
    one strand per two base pairs (the protein moves around the DNA, not the DNA around the protein)."""
    return sum(len(m.axial) / 4 if isinstance(m, PlacedStem) else 1 for m in it.members)


def _energy(flat, home, anchor_w, want, pull, loops, lam, mids=None, fixed=None):
    """Layout energy and gradient over only the overlapping boxes and the contact pairs `pull` = (i, j, w), i < j.
    Same terms as _energy_dense."""
    P = flat.reshape(-1, 2)
    grad = 2 * anchor_w * (P - home)
    E = float((anchor_w[:, 0] * ((P - home) ** 2).sum(1)).sum())
    Q = P if mids is None else P + mids  # footprint centres
    E += _fixed_overlap(Q, fixed, lam, grad)
    oi, oj = np.nonzero(np.triu(want[..., 0] > np.abs(Q[None, :, 0] - Q[:, None, 0]), 1))
    keep = want[oi, oj, 1] > np.abs(Q[oj, 1] - Q[oi, 1])
    oi, oj = oi[keep], oj[keep]
    if len(oi):  # overlapping boxes: pushed apart along the axis of least overlap
        d = Q[oj] - Q[oi]
        px, py = want[oi, oj, 0] - np.abs(d[:, 0]), want[oi, oj, 1] - np.abs(d[:, 1])
        on_x = px < py
        E += lam * float((np.where(on_x, px, py) ** 2).sum())
        _pair_grad(grad, oi, oj, d, np.where(on_x, 2 * lam * px, 0.0), np.where(on_x, 0.0, 2 * lam * py))
    ci, cj, cw = pull
    if len(ci):  # touching in 3D but apart on the page: pulled together until the boxes meet
        d = Q[cj] - Q[ci]
        px, py = want[ci, cj, 0] - np.abs(d[:, 0]), want[ci, cj, 1] - np.abs(d[:, 1])
        gap = np.maximum(-px, -py)  # how far apart two boxes are beyond touching
        cg = np.where(gap > 0, 2 * cw * gap, 0.0)
        E += float((np.where(gap > 0, cw * gap**2, 0.0)).sum())
        on_x = -px >= -py
        _pair_grad(grad, ci, cj, d, -np.where(on_x, cg, 0.0), -np.where(on_x, 0.0, cg))
    I, OI, J, OJ, w, slack = loops
    if len(I):
        g = (P[J] + OJ) - (P[I] + OI)
        dist = np.maximum(np.sqrt(g[:, 0] * g[:, 0] + g[:, 1] * g[:, 1]), 1e-12)
        excess = np.clip(dist - slack, 0.0, None)
        E += float((w * excess**2).sum())
        dg = (2 * w * excess / dist)[:, None] * g
        np.add.at(grad, J, dg)
        np.add.at(grad, I, -dg)
    return E, grad.ravel()


def _fixed_overlap(Q, fixed, lam, grad) -> float:
    """lam * weight * squared overlap of each footprint (centres Q) with each fixed box, `fixed` = (box centres
    (m, 2), wanted centre distances (n, m, 2): half extents of both plus MARGIN, weights (n, m)), along the axis of
    least overlap; the gradient (added to `grad` in place) moves the footprints only."""
    if fixed is None or not len(fixed[0]):
        return 0.0
    centres, want, weight = fixed
    d = centres[None, :, :] - Q[:, None, :]
    px, py = want[..., 0] - np.abs(d[..., 0]), want[..., 1] - np.abs(d[..., 1])
    hit, on_x = (px > 0) & (py > 0), px < py
    w = lam * weight
    grad[:, 0] += np.where(hit & on_x, 2 * w * px * np.sign(d[..., 0]), 0.0).sum(1)
    grad[:, 1] += np.where(hit & ~on_x, 2 * w * py * np.sign(d[..., 1]), 0.0).sum(1)
    return float((w * np.where(hit, np.where(on_x, px, py), 0.0) ** 2).sum())


def _pair_grad(grad, i, j, d, gx, gy) -> None:
    """Add dE/d(overlap depth) of pairs (i, j), d = Q_j - Q_i, to the per-item gradient (in place)."""
    n = len(grad)
    fx, fy = gx * np.sign(d[:, 0]), gy * np.sign(d[:, 1])
    grad[:, 0] += np.bincount(i, fx, n) - np.bincount(j, fx, n)
    grad[:, 1] += np.bincount(i, fy, n) - np.bincount(j, fy, n)


def _energy_dense(flat, home, anchor_w, want, W, loops, lam, mids=None, fixed=None):
    """Layout energy and its gradient. Terms: stay near the projected position (fold); loop ends within
    _LOOP_SLACK of each other; touching elements' (and DNA binders') boxes adjacent;
    lam * squared box overlap (raised stage by stage until overlaps are gone), with each other and with the
    `fixed` boxes (_fixed_overlap)."""
    P = flat.reshape(-1, 2)
    n = len(P)
    grad = 2 * anchor_w * (P - home)
    E = float((anchor_w[:, 0] * ((P - home) ** 2).sum(1)).sum())
    Q = P if mids is None else P + mids  # footprint centres
    E += _fixed_overlap(Q, fixed, lam, grad)
    d = Q[None, :, :] - Q[:, None, :]  # d[i, j] = Q_j - Q_i
    sg = np.sign(d)
    px, py = want[..., 0] - np.abs(d[..., 0]), want[..., 1] - np.abs(d[..., 1])  # overlap depth per axis
    upper = np.triu(np.ones((n, n), bool), 1)

    def pair_grad(gx, gy):  # dE/d(overlap depth) per pair (upper triangle) -> per item
        grad[:, 0] += (gx * sg[..., 0]).sum(1) - (gx * sg[..., 0]).sum(0)
        grad[:, 1] += (gy * sg[..., 1]).sum(1) - (gy * sg[..., 1]).sum(0)

    on_x = px < py
    ov = upper & (px > 0) & (py > 0)
    depth = np.where(on_x, px, py)
    E += lam * float((depth[ov] ** 2).sum())
    pair_grad(np.where(ov & on_x, 2 * lam * px, 0.0), np.where(ov & ~on_x, 2 * lam * py, 0.0))
    gap = np.maximum(-px, -py)  # how far apart two boxes are beyond touching
    ct = upper & (W > 0) & (gap > 0)
    E += float((W[ct] * gap[ct] ** 2).sum())
    cg = np.where(ct, 2 * W * gap, 0.0)
    gx_side = -px >= -py
    pair_grad(-np.where(gx_side, cg, 0.0), -np.where(~gx_side, cg, 0.0))
    I, OI, J, OJ, w, slack = loops
    if len(I):
        g = (P[J] + OJ) - (P[I] + OI)
        dist = np.maximum(np.sqrt(g[:, 0] * g[:, 0] + g[:, 1] * g[:, 1]), 1e-12)
        excess = np.clip(dist - slack, 0.0, None)
        E += float((w * excess**2).sum())
        dg = (2 * w * excess / dist)[:, None] * g
        np.add.at(grad, J, dg)
        np.add.at(grad, I, -dg)
    return E, grad.ravel()


def _snap(v) -> np.ndarray:
    """Round to the 1/_GRID grid, so noise of ~1e-16 in a projection (np.linalg.svd and small dot products
    differ between BLAS kernels) never reaches the optimiser, which would amplify it, nor the drawn geometry.
    The rounding removes that noise except when a value sits within noise of a half-grid point (probability
    ~1e-12 per value); exact libm values such as cos 0 = 1 are safe."""
    return np.round(np.asarray(v, float) * _GRID) / _GRID


def _lbfgs(fun, x0, maxiter=300, m=10, pgtol=1e-5, ftol=2.2e-9):
    """Minimise fun (x -> (value, gradient)) from x0 by L-BFGS (two-loop recursion, last m steps) with a
    strong-Wolfe line search; stops like scipy's L-BFGS-B: largest gradient component <= pgtol, relative
    decrease <= ftol, or maxiter iterations. numpy only and every dot product an elementwise (a * b).sum(): no
    BLAS, whose results depend on the CPU kernel, so the same input gives the same minimum everywhere."""
    x = np.array(x0, float)
    f, g = fun(x)
    S, Y, R = [], [], []  # steps, gradient changes, 1 / (s . y)
    for _ in range(maxiter):
        if np.abs(g).max() <= pgtol:
            break
        d = -g
        alpha = []
        for s, y, r in zip(reversed(S), reversed(Y), reversed(R)):
            a = r * (s * d).sum()
            alpha.append(a)
            d -= a * y
        if S:
            d *= (S[-1] * Y[-1]).sum() / (Y[-1] * Y[-1]).sum()
        for s, y, r, a in zip(S, Y, R, reversed(alpha)):
            d += (a - r * (y * d).sum()) * s
        if S and (g * d).sum() >= 0:  # not a descent direction: forget the curvature pairs
            S, Y, R = [], [], []
            d = -g
        step = 1.0 if S else min(1.0 / np.sqrt((d * d).sum()), 1e10)
        found = _wolfe(fun, x, f, g, d, step)
        if found is None:
            if not S:
                break  # no step along steepest descent satisfies the conditions
            S, Y, R = [], [], []  # retry once from steepest descent, as L-BFGS-B does
            continue
        step, f_new, g_new = found
        s, y = step * d, g_new - g
        x, f_old, f, g = x + s, f, f_new, g_new
        if np.abs(g).max() <= pgtol or f_old - f <= ftol * max(abs(f_old), abs(f), 1.0):
            break
        sy = (s * y).sum()
        if sy > np.finfo(float).eps * (y * y).sum():  # keep only pairs with positive curvature
            S.append(s)
            Y.append(y)
            R.append(1.0 / sy)
            if len(S) > m:
                del S[0], Y[0], R[0]
    return x


def _wolfe(fun, x, f0, g0, d, step, c1=1e-3, c2=0.9, evals=20):
    """A step along d meeting the strong Wolfe conditions (sufficient decrease c1, curvature c2): grow the
    step until it brackets one, then shrink the bracket by safeguarded cubic interpolation. Returns
    (step, value, gradient) there, or None after `evals` evaluations."""
    slope0 = (g0 * d).sum()
    lo, f_lo, s_lo = 0.0, f0, slope0  # the best step so far that decreases enough
    hi = None
    for k in range(evals):
        if hi is not None:  # zoom: cubic through both ends, kept in the bracket's inner 80 %, else bisect
            a, fa, sa, b, fb, sb = lo, f_lo, s_lo, *hi
            if a == b:
                return None  # the bracket has shrunk to a point
            d1 = sa + sb - 3 * (fa - fb) / (a - b)
            rad = d1 * d1 - sa * sb
            t = (a + b) / 2
            if rad >= 0:
                d2 = np.copysign(np.sqrt(rad), b - a)
                den = sb - sa + 2 * d2
                if den != 0:
                    c = b - (b - a) * (sb + d2 - d1) / den
                    if min(a, b) + 0.1 * abs(b - a) <= c <= max(a, b) - 0.1 * abs(b - a):
                        t = c
            step = float(t)
        f, g = fun(x + step * d)
        slope = (g * d).sum()
        if f > f0 + c1 * step * slope0 or (f >= f_lo and (hi is not None or k > 0)):
            hi = (step, f, slope)  # too far: the minimum lies between lo and here
        elif abs(slope) <= -c2 * slope0:
            return step, f, g
        else:
            if hi is not None and slope * (hi[0] - lo) >= 0:
                hi = (lo, f_lo, s_lo)
            elif hi is None and slope >= 0:
                hi = (lo, f_lo, s_lo)
            lo, f_lo, s_lo = step, f, slope
            if hi is None:
                step *= 2.0  # still going down: try further
    return None


def _attract(
    items: list[_Item],
    sses: list[SSE],
    contacts: dict[frozenset[str], int],
    ties: dict[int, tuple[int, int]] | None = None,
    spread: float = 0.0,
    lift: float = 0.0,
    bridges: list[tuple[str, float, str, float]] = (),
    extra: list[tuple[int, np.ndarray, int, np.ndarray, float, float]] = (),
    joints: set[tuple[int, int]] = frozenset(),
    fixed: Walls | None = None,
    mirror: float | None = None,
) -> None:
    """Minimise the layout energy (_energy) from the projected positions, raising the overlap penalty in
    stages. The result becomes each item's new home; _separate then guarantees no overlap remains. Every input
    is snapped to the grid (_snap) and _lbfgs runs without BLAS, so the minimum is the same on every machine.
    `extra`: further springs (item i, port offset, item j, port offset, weight, slack), the nucleic-acid strands
    from stem to stem (_stem_springs). `joints` (_joints): pairs of items kept apart without MARGIN, so their
    boxes may touch and no element fits between them for long. `fixed` (Walls): boxes that stay where they are (the
    nucleotides of unpaired strands' traces), the items kept MARGIN off them, firmly or (an element touching that
    trace) not at all. `mirror` (page x of a dyad in the page): each copy of `ties` sits at its mate's mirror image
    about it, at the same height, instead of `spread` along."""
    n = len(items)
    W = _contact_weights(items, contacts)
    ci, cj = np.nonzero(np.triu(W > 0, 1))
    pull = (ci, cj, W[ci, cj])
    springs = _loop_springs(items, sses)
    weight = [_W_SEQ] * len(springs)
    slack = [_LOOP_SLACK] * len(springs)
    for i, oi, j, oj, w, sl in extra:  # nucleic-acid strands from stem to stem
        springs.append((i, oi, j, oj))
        weight.append(float(_snap(w)))
        slack.append(float(_snap(sl)))
    owner = _owner(items)
    bridges_used = []
    for a, ta, b, tb in bridges:  # disulfides: the two cysteines themselves pulled together
        if a in owner and b in owner and owner[a][0] != owner[b][0]:
            bridges_used.append((a, b))
            (i, ma, oa), (j, mb, ob) = owner[a], owner[b]
            springs.append(
                (
                    i,
                    np.add(oa, ma.direction * (ta - 0.5) * ma.length),
                    j,
                    np.add(ob, mb.direction * (tb - 0.5) * mb.length),
                )
            )
            weight.append(_W_SS)
            slack.append(_SS_SLACK)
    springs = [(i, _snap(oi), j, _snap(oj)) for i, oi, j, oj in springs]
    loops = (
        np.array([s[0] for s in springs], int),
        np.array([s[1] for s in springs], float).reshape(-1, 2),
        np.array([s[2] for s in springs], int),
        np.array([s[3] for s in springs], float).reshape(-1, 2),
        np.array(weight, float),
        np.array(slack, float),
    )
    home = _snap([it.home for it in items])
    h = _snap([it.half for it in items])
    want = h[:, None, :] + h[None, :, :] + MARGIN
    joined = _joint_matrix(n, joints)  # consecutive pieces of a bent duplex: they may touch, with no margin between
    want[joined] = (h[:, None, :] + h[None, :, :])[joined]
    anchor_w = _K_HOME * np.array([[_mass(it)] for it in items], float)
    mids = _snap([it.mid for it in items])
    boxes = np.asarray(fixed.boxes if fixed is not None else [], float).reshape(-1, 4)
    centres, halves = _snap((boxes[:, :2] + boxes[:, 2:]) / 2), _snap((boxes[:, 2:] - boxes[:, :2]) / 2)
    weight = fixed.hard.astype(float) if fixed is not None else np.ones((n, 0))
    walls = (centres, h[:, None, :] + halves[None, :, :] + MARGIN, weight)
    spread, lift = float(_snap(spread)), float(_snap(lift))
    mirror = None if mirror is None else float(_snap(mirror))
    ties = ties or {}
    free = [i for i in range(n) if i not in ties]
    slot = {i: k for k, i in enumerate(free)}
    base = np.array([slot[ties[i][0]] if i in ties else slot[i] for i in range(n)])
    copy_k = np.array([ties[i][1] if i in ties else 0 for i in range(n)], float)
    flip = np.stack([np.where(copy_k > 0, -1.0, 1.0), np.ones(n)], axis=1)  # a mirrored copy: x against its mate

    def expand(z):
        P = z[:-1].reshape(-1, 2)[base].copy()
        if mirror is not None:  # protomer 1 is protomer 0 seen in the mirror line
            P[copy_k > 0, 0] = 2 * mirror - P[copy_k > 0, 0]
            return P
        P[:, 0] += copy_k * z[-1]  # protomer k sits k spreads to the right of protomer 0
        P[:, 1] += copy_k * lift  # and, along a filament, k rises higher
        return P

    def tied(z):  # one position per orbit of symmetry copies, plus the spread between protomers
        if n > _PAIR_ENERGY_ABOVE:
            E, g = _energy(expand(z).ravel(), home, anchor_w, want, pull, loops, lam, mids, walls)
        else:
            E, g = _energy_dense(expand(z).ravel(), home, anchor_w, want, W, loops, lam, mids, walls)
        g = g.reshape(-1, 2)
        gz = np.zeros((len(free), 2))
        if mirror is not None:  # the spread is fixed at 0 (no gradient)
            np.add.at(gz, base, g * flip)
            return E, np.append(gz.ravel(), 0.0)
        np.add.at(gz, base, g)
        return E, np.append(gz.ravel(), float((g[:, 0] * copy_k).sum()))

    def solve(z, stages=_LAMBDAS):
        nonlocal lam
        for lam in stages:  # noqa: B007 - read by tied() through nonlocal
            z = _lbfgs(tied, z)
        return z

    lam = _LAMBDAS[0]
    z = solve(np.append(_snap([[items[i].cx, items[i].cy] for i in free]).ravel(), spread))
    lam = _LAMBDAS[-1]
    best = tied(z)[0]
    mass = np.array([_mass(it) for it in items])
    for k in range(len(springs) - len(bridges_used), len(springs)):  # disulfides still stretched: a second start
        i, oi, j, oj = springs[k]
        P = expand(z)
        gx, gy = (P[j] + oj) - (P[i] + oi)
        if np.sqrt(gx * gx + gy * gy) <= 2 * _SS_SLACK + 1.0:
            continue
        mover, anchor, om, oa = (i, j, oi, oj) if mass[i] <= mass[j] else (j, i, oj, oi)
        if mover in ties:
            continue
        hm, ha = h[mover], h[anchor]
        gap_x, gap_y = hm[0] + ha[0] + MARGIN, hm[1] + ha[1] + MARGIN
        for side in ((0.0, -1.0), (-gap_x, 0.0), (gap_x, 0.0), (0.0, gap_y), (0.0, -gap_y)):
            # the lighter element brought beside its partner: just below the cysteine, or clear of its footprint
            trial = _snap(z)  # every start on the grid
            trial[2 * slot[mover] : 2 * slot[mover] + 2] = _snap(P[anchor] + oa - om + np.array(side))
            trial = solve(trial, _LAMBDAS[1:])
            lam = _LAMBDAS[-1]
            energy = tied(trial)[0]
            if energy < best:
                z, best = trial, energy
    x = expand(z).ravel()
    for it, p in zip(items, x.reshape(n, 2)):
        it.cx, it.cy = float(p[0]), float(p[1])
        it.home = (it.cx, it.cy)


def _collide(a: _Item, ha, ma, b: _Item, hb, mb) -> bool:
    """Do the footprints (half extents h, centre offsets m) of a and b, plus MARGIN, overlap?"""
    dx, dy = b.cx + mb[0] - a.cx - ma[0], b.cy + mb[1] - a.cy - ma[1]
    return ha[0] + hb[0] + MARGIN - abs(dx) > 1e-9 and ha[1] + hb[1] + MARGIN - abs(dy) > 1e-9


def _stack(items: list[_Item], depth: bool = False) -> list[list[_Item]]:
    """Rows in sequence order, one row break per chain and when a row is full. With `depth` (a membrane is
    drawn) every item keeps its height along the frame's up: whole chains side by side, packed in chain order
    into rows no wider than _MAX_ROW, each row moved down as a whole. The rows are returned."""
    items.sort(key=lambda it: (it.chain, it.order))
    if depth:

        def width(row: list[_Item]) -> float:
            return sum(it.box[2] - it.box[0] for it in row) + _ROW_GAP * (len(row) - 1)

        rows = []
        for c in dict.fromkeys(it.chain for it in items):
            mine = [it for it in items if it.chain == c]
            # a single chain wider than _MAX_ROW has a row of its own: its elements keep their depth, so it
            # cannot fold back
            if rows and width(rows[-1] + mine) <= _MAX_ROW:
                rows[-1] += mine
            else:
                rows.append(mine)
        bottom = None
        for row in rows:
            ys = [float(_snap(it.home[1])) for it in row]
            top = max(y + it.box[3] for y, it in zip(ys, row))
            shift = 0.0 if bottom is None else float(_snap(bottom - _ROW_SPACING - top))
            x = 0.0
            for y, it in zip(ys, row):
                it.cx = x - it.box[0]
                it.cy = y + shift
                x += it.box[2] - it.box[0] + _ROW_GAP
            bottom = min(it.cy + it.box[1] for it in row)
        return rows
    rows: list[list[_Item]] = []
    width = 0.0
    for it in items:
        hx = it.half[0]
        new_chain = rows and rows[-1][-1].chain != it.chain
        if not rows or new_chain or width + 2 * hx > _MAX_ROW:
            rows.append([])
            width = 0.0
        rows[-1].append(it)
        width += 2 * hx + _ROW_GAP
    y_top = 0.0
    for row in rows:
        height = max(it.box[3] - it.box[1] for it in row)
        x = 0.0
        for it in row:  # footprints (with label room) side by side, centred in the row's height
            it.cx = x - it.box[0]
            it.cy = y_top - height / 2 - it.mid[1]
            x += it.box[2] - it.box[0] + _ROW_GAP
        dna_row = all(isinstance(m, PlacedStem) for it in row for m in it.members)
        y_top -= height + (MARGIN if dna_row else _ROW_SPACING)  # binders hang right under the DNA
    return rows


def _centre_dna_over_binders(items: list[_Item], contacts: dict[frozenset[str], int]) -> None:
    owner = _owner(items)
    for it in items:
        dna = [m for m in it.members if isinstance(m, PlacedStem)]
        if not dna:
            continue
        xs = [
            items[owner[x][0]].cx + owner[x][2][0]
            for key in contacts
            if dna[0].id in key
            for x in key
            if x != dna[0].id and x in owner
        ]
        if xs:
            it.cx = float(np.mean(xs))


def _assign_labels(placed: list[Placed]) -> None:
    counts: dict[tuple[str, str], int] = {}
    for p in sorted(placed, key=lambda p: (p.sse.chain, p.sse.start)):
        n = counts.get((p.sse.chain, p.sse.kind), 0)
        counts[(p.sse.chain, p.sse.kind)] = n + 1
        p.label = generic_label(p.sse.kind, n)


def stem_pieces(na: Nucleic) -> list[Stem]:
    """Every stem of the nucleic acid, long bent ones cut into straight pieces, in 5'->3' order."""
    return [piece for st in na.structure.stems for piece in segment_stem(st, na.strands)]


_MIN_RISE = 2.6  # Å per base pair, the least a stem is drawn with (A-form rises 2.8, B-form 3.4)


def _rungs(piece: Stem, na: Nucleic, centre: np.ndarray, axis: np.ndarray) -> dict[tuple[int, int], float]:
    """Each nucleotide's place along the drawn stem, centred: one straight rung per pair, the pairs spread over
    the stem's length along its axis (so a curved stem's ends never crowd), each step between two pairs at
    least (1 + the most nucleotides bulged out on either arm there) * _MIN_RISE long, bulged nucleotides spaced
    evenly between their neighbours on their own arm."""
    pairs = piece.pairs
    proj = [float(((na.strands[p.a].c1[p.i] + na.strands[p.b].c1[p.j]) / 2 - centre) @ axis) * SCALE for p in pairs]
    steps = [1 + max(q.i - p.i - 1, p.j - q.j - 1) for p, q in zip(pairs, pairs[1:])]
    span = max(abs(proj[-1] - proj[0]), sum(steps) * _MIN_RISE * SCALE)
    at = np.concatenate([[0], np.cumsum(steps)]) / max(sum(steps), 1)  # each pair's fraction of the way along
    out: dict[tuple[int, int], float] = {}
    for k, p in enumerate(pairs):
        out[(p.a, p.i)] = out[(p.b, p.j)] = -span / 2 + span * float(at[k])
    for p, q in zip(pairs, pairs[1:]):  # bulges: evenly between the two pairs around them, on their own arm
        x0, x1 = out[(p.a, p.i)], out[(q.a, q.i)]
        for arm, (s, lo, hi) in enumerate(((p.a, p.i, q.i), (p.b, q.j, p.j))):
            inner = range(lo + 1, hi)  # the second arm runs back along the stem (also in a hairpin, a == b)
            for m, i in enumerate(reversed(inner) if arm else inner):
                out[(s, i)] = x0 + (x1 - x0) * (m + 1) / (len(inner) + 1)
    return out


def _side_signal(piece: Stem, na: Nucleic, frame: Frame, centre: np.ndarray, angle: float) -> float:
    """How far (Å, on the grid) the first arm is seen to the left of the drawn direction `angle` from the second:
    its mean C1' against the second arm's. Its sign is the side that fits the projection best, so the drawing
    turns the stem, never mirrors it; near zero (the arms seen nearly on top of each other) both sides fit about
    equally. Not the outermost pair's C1'-C1' line: it turns with the helical twist, so its side depends on where
    in the turn the stem starts."""
    (a, i0, i1), (b, j0, j1) = piece.arm5, piece.arm3
    across = np.mean([na.strands[a].c1[i] for i in range(i0, i1 + 1)], axis=0) - np.mean(
        [na.strands[b].c1[j] for j in range(j0, j1 + 1)], axis=0
    )
    seen = frame.direction(centre, across)
    return float(_snap(np.cos(angle) * seen[1] - np.sin(angle) * seen[0]))


def _hairpin(piece: Stem, na: Nucleic) -> int:
    """Nucleotides of the hairpin loop closed by the piece's innermost pair, or -1 if that pair closes none."""
    a, i, b, j = piece.pairs[-1][:4]
    return next((r.hi - r.lo + 1 for r in na.structure.runs if r.left == (a, i) and r.right == (b, j)), -1)


def _stem_item(piece: Stem, na: Nucleic, frame: Frame, flat: bool, width: float = NA_W) -> _Item:
    """A stem piece as a rigid item at its projected place, turned onto the page (never mirrored). Its axis
    differs in the last bit between BLAS kernels, so every number the layout takes from it is snapped (_snap)."""
    centre, axis = stem_axis(piece, na.strands)
    (a, i0, i1), (b, j0, j1) = piece.arm5, piece.arm3
    keys = [(a, i) for i in range(i0, i1 + 1)] + [(b, j) for j in range(j0, j1 + 1)]
    # on the grid; + 0.0 turns a -0.0 (a centre pair a hair below zero) into 0.0
    along = {key: float(_snap(v)) + 0.0 for key, v in _rungs(piece, na, centre, axis).items()}
    length = float(_snap(max(along.values()) - min(along.values()) + DNA_RISE * SCALE))  # half a step beyond each end
    a2 = frame.direction(centre, axis)
    angle = 0.0 if flat or np.linalg.norm(a2) < 1e-6 else float(_snap(np.arctan2(a2[1], a2[0])))
    signal = _side_signal(piece, na, frame, centre, angle)
    first = 1 if signal >= 0 else -1
    side = {k: first if n <= i1 - i0 else -first for n, k in enumerate(keys)}  # keys: first arm, then second
    x, y = (float(v) for v in _snap(frame.project(centre)[0] * SCALE))
    d = PlacedStem(
        piece.id,
        piece.name,
        [a, b],
        [tuple(p[:4]) for p in piece.pairs],
        0.0,
        0.0,
        length,
        angle,
        width,
        along,
        side,
        set(piece.bulges),
        piece.parent,
        first,
        abs(signal),
        _hairpin(piece, na),
    )
    return _Item([d], [(0.0, 0.0)], x, y, (x, y), -1, "")


# Å: a stem seen with its arms at least this far apart keeps the side it is seen on. On the test structures every
# stem below it fits its projection about as well turned round; the closest one that does not is 4YBB P5 at 2.81 Å.
_SIDE_CLEAR = 2.75
_CYCLE_MAX = 12  # unclear stems linked round a cycle: up to this many are settled by trying every turning


def _cross2(u: np.ndarray, v: np.ndarray) -> float:
    return float(u[0] * v[1] - u[1] * v[0])


def _continuous_sides(stems: list[PlacedStem], na: Nucleic) -> None:
    """Turn round (swap the rails of) the stems whose side the projection cannot tell (margin < _SIDE_CLEAR), so
    their strands run on uncrossed to the stems they are joined to both ways: the next piece of a bent duplex,
    or the stem across an internal loop or bulge (one strand goes from A to B, another from B back to A).
    Stems seen clearly keep their side. Each group of unclear stems linked to one another is settled exactly
    (_two_colour): fewest crossed links, any crossing that cannot be avoided at the link whose stems are seen
    least clearly, then fewest stems turned; so the result does not depend on the order of the stems."""
    at = {key: k for k, d in enumerate(stems) for key in d.axial}
    steps: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
    for s, strand in enumerate(na.strands):
        held = [i for i in range(len(strand)) if (s, i) in at]
        for i, j in zip(held, held[1:]):
            if at[(s, i)] != at[(s, j)]:
                steps.setdefault((at[(s, i)], at[(s, j)]), []).append((s, i, j))
    links: list[tuple[int, int, int]] = []  # (stem, stem, +1: same turn wanted, -1: opposite turns wanted)
    for (k, m), there in sorted(steps.items()):
        back = steps.get((m, k), [])
        if k > m or len(there) != 1 or len(back) != 1:
            continue
        (s, i, j), (t, q, r) = there[0], back[0]
        a, b = stems[k], stems[m]
        # the two chords A(s,i)->B(s,j) and B(t,q)->A(t,r) cross when, seen along the way from A's end to B's,
        # A's two rail points lie the other way round from B's (turning a stem round swaps its rail points)
        u = a.port((s, i), True)[1] - b.port((s, j), False)[1]
        ca = float(_snap(_cross2(u, a.point((t, r)) - a.point((s, i)))))
        cb = float(_snap(_cross2(u, b.point((t, q)) - b.point((s, j)))))
        if ca and cb:
            links.append((k, m, 1 if (ca > 0) == (cb > 0) else -1))
    fixed = {k for k, d in enumerate(stems) if d.margin >= _SIDE_CLEAR}
    turn = dict.fromkeys(fixed, 1)
    nbr: dict[int, list[int]] = {}
    for k, m, _ in links:
        if k not in fixed and m not in fixed:
            nbr.setdefault(k, []).append(m)
            nbr.setdefault(m, []).append(k)
    for k in range(len(stems)):
        if k in turn:
            continue
        group, todo = {k}, [k]
        while todo:
            for m in nbr.get(todo.pop(), ()):
                if m not in group:
                    group.add(m)
                    todo.append(m)
        mine = [ln for ln in links if ln[0] in group or ln[1] in group]
        turn.update(_two_colour(sorted(group), mine, nbr, turn, stems))
    for k, d in enumerate(stems):
        if turn[k] < 0:
            d.first = -d.first
            d.side = {key: -v for key, v in d.side.items()}


def _two_colour(group, links, nbr, fixed, stems) -> dict[int, int]:
    """+1 (keep) or -1 (turn round) for each stem of `group` (sorted), given the settled stems in `fixed`: the
    least cost over `links` (k, m, wanted sign), counting each crossed link as one, below that the rank of how
    clearly its stems are seen (least clear cheapest), below that each stem turned. Exact by dynamic programming
    on a tree (chains of pieces and loops are trees); a group with a cycle is solved by trying every turning up
    to 2**_CYCLE_MAX, beyond that on a spanning tree."""
    order = sorted(
        links,
        key=lambda ln: (
            max(stems[ln[0]].margin, stems[ln[1]].margin),
            min(stems[ln[0]].margin, stems[ln[1]].margin),
            ln,
        ),
    )
    rank = {ln: r for r, ln in enumerate(order)}
    big = len(links) * len(links) + 1  # one crossing fewer outweighs any ranks
    unit = len(group) + 1  # and any ranks outweigh the number of stems turned

    def cost(ln) -> int:
        return (big + rank[ln]) * unit

    inner = [ln for ln in links if ln[0] in group and ln[1] in group]
    if len(inner) >= len(group) and len(group) <= _CYCLE_MAX:  # a cycle: every turning
        best, pick = None, None
        for signs in itertools.product((1, -1), repeat=len(group)):
            t = {**fixed, **dict(zip(group, signs))}
            c = sum(cost(ln) for ln in links if t[ln[0]] * t[ln[1]] != ln[2]) + signs.count(-1)
            if best is None or c < best:
                best, pick = c, dict(zip(group, signs))
        return pick
    edge: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
    for ln in links:
        edge.setdefault((ln[0], ln[1]), []).append(ln)
        edge.setdefault((ln[1], ln[0]), []).append(ln)
    root, parent, bfs = group[0], {}, [group[0]]
    children: dict[int, list[int]] = {v: [] for v in group}
    for v in bfs:
        for m in sorted(nbr.get(v, ())):
            if m in children and m != root and m not in parent:
                parent[m] = v
                children[v].append(m)
                bfs.append(m)

    def unary(v, s) -> int:  # links to settled stems, and the turn itself
        out = 0 if s > 0 else 1
        for ln in links:
            other = ln[1] if ln[0] == v else ln[0] if ln[1] == v else None
            if other is not None and other in fixed and s * fixed[other] != ln[2]:
                out += cost(ln)
        return out

    def between(v, s, c, t) -> int:
        return sum(cost(ln) for ln in edge.get((v, c), ()) if s * t != ln[2])

    best: dict[int, dict[int, int]] = {}
    pick: dict[tuple[int, int], int] = {}
    for v in reversed(bfs):
        best[v] = {}
        for s in (1, -1):
            total = unary(v, s)
            for c in children[v]:
                options = [(between(v, s, c, t) + best[c][t], -t) for t in (1, -1)]
                value, neg = min(options)
                pick[(c, s)] = -neg
                total += value
            best[v][s] = total
    out = {root: 1 if best[root][1] <= best[root][-1] else -1}
    for v in bfs[1:]:
        out[v] = pick[(v, out[parent[v]])]
    return out


_TRACE_HALF = (
    0.45  # page units: a bead's width, which a trace nucleotide's box reaches either side of it (na_paths._BEAD)
)


_WALL_MIN = 3  # a trace of fewer nucleotides (a dinucleotide ligand-like oligo) is drawn but is no wall


def _trace_walls(items: list[_Item], traces: dict[int, np.ndarray], na: Nucleic) -> Walls:
    """Each nucleotide of a trace of at least _WALL_MIN as a fixed box (_TRACE_HALF either side), which the items
    keep off, but for an item holding an element that touches that trace: a single strand wrapped over its protein
    (1EYG's ssDNA across the SSB tetramer) does not push away the protein it lies on."""
    boxes, owner = [], []
    for s, pts in traces.items():
        if len(pts) < _WALL_MIN:
            continue
        for x, y in np.asarray(pts, float).tolist():
            boxes.append((x - _TRACE_HALF, y - _TRACE_HALF, x + _TRACE_HALF, y + _TRACE_HALF))
            owner.append(s)
    hard = np.ones((len(items), len(boxes)), bool)
    for k, it in enumerate(items):
        touched = {
            x[0]
            for m in it.members
            if isinstance(m, Placed) and not m.ghost
            for r in range(m.sse.start, m.sse.end + 1)
            for x in na.contacts.get(r, ())
        }
        hard[k] = [s not in touched for s in owner]
    return Walls(np.asarray(boxes, float).reshape(-1, 4), hard)


_NUDGES = tuple(f * k for k in (1.05, 1.25, 1.5, 2.0) for f in (1, -1))  # sideways steps, in beads (_TRACE_HALF)


def _room(others: np.ndarray, xy) -> float:
    """How near the nearest of `others` ((n, 2) page points) comes to xy."""
    return float(np.hypot(*(others - xy).T).min()) if len(others) else np.inf


def _traces(na: Nucleic, frame: Frame) -> dict[int, np.ndarray]:
    """Strands with no base pair at all, drawn along their projected backbone (page units, on the grid). A
    nucleotide projected within a bead's width (_TRACE_HALF) of one drawn before it (another strand's, or its own,
    its 5' neighbour seen end-on included) is nudged sideways off its trace, as little as clears it (else as far as
    it can): two strands crossing in projection, or one running toward the viewer (1EYG), would draw their beads
    on each other."""
    held = {n for st in na.structure.stems for n in st.nucleotides()}
    out: dict[int, np.ndarray] = {}
    traced = [
        s
        for s, strand in enumerate(na.strands)
        if len(strand) > 1 and not any((s, i) in held for i in range(len(strand)))
    ]
    taken = np.zeros((sum(len(na.strands[s]) for s in traced), 2))  # every trace nucleotide placed so far
    count = 0
    for s in traced:
        pts = _snap(frame.project(na.strands[s].p) * SCALE)
        for i in range(len(pts)):
            others = taken[:count]
            space = _room(others, pts[i])
            if space < _TRACE_HALF:
                along = pts[min(i + 1, len(pts) - 1)] - pts[max(i - 1, 0)]
                if float(np.hypot(*along)) < 1e-9:  # its neighbours on one spot: its 5' neighbour's way, or across
                    along = pts[i] - pts[i - 1] if i > 0 else np.zeros(2)
                if float(np.hypot(*along)) < 1e-9:
                    along = np.array([0.0, 1.0])
                side = np.array([-along[1], along[0]]) / float(np.hypot(*along))
                for k in _NUDGES:
                    xy = _snap(pts[i] + side * k * _TRACE_HALF)
                    r = _room(others, xy)
                    if r > space:
                        pts[i], space = xy, r
                    if r >= _TRACE_HALF:
                        break
            taken[count] = pts[i]
            count += 1
        out[s] = pts
    return out


def _nearest_piece(na: Nucleic, pieces: list[Stem]) -> dict[tuple[int, int], int]:
    """Each nucleotide's stem piece: its own if paired (or bulged) in one, otherwise the piece of the nearest
    such nucleotide along its strand, so a protein holding a loop or a free end is pulled to the stems beside it."""
    own = {key: k for k, piece in enumerate(pieces) for key in piece.nucleotides()}
    out = dict(own)
    for s, strand in enumerate(na.strands):
        held = [i for i in range(len(strand)) if (s, i) in own]
        if not held:
            continue
        for i in range(len(strand)):
            if (s, i) not in own:
                out[(s, i)] = own[(s, min(held, key=lambda h: (abs(h - i), h)))]
    return out


def _stem_contacts(sses: list[SSE], na: Nucleic, pieces: list[Stem]) -> dict[frozenset[str], int]:
    """Residues of each element touching each stem (or the single strands beside it), as extra adjacency pull."""
    piece_of = _nearest_piece(na, pieces)
    out: dict[frozenset[str], int] = {}
    for s in sses:
        for k, piece in enumerate(pieces):
            n = sum(1 for r in range(s.start, s.end + 1) if any(piece_of.get(x) == k for x in na.contacts.get(r, ())))
            if n:
                out[frozenset((s.id, piece.id))] = 4 * n  # one residue in the groove counts like several CA pairs
    return out


_W_JOINT = 4.0  # pull between consecutive pieces of one bent duplex: the chain follows the duplex
_JOINT_SLACK = 0.3
# a strand passing through at most _SHORT_RUN nucleotides between two stems of more than two strands (a strand
# exchange, as at a Holliday junction) holds the stem ends close, firmer than a joint, so arms no protein holds come
# in too. Not within one RNA chain: there it pulled a tRNA's T arm across its cloverleaf.
_SHORT_RUN = 2
_W_SHORT = 10.0
# two stems of one helix (nucleic._helices) with at most this many nucleotides on either side of the loop between
# them are held close like a strand exchange. Measured over the 39 campaign entries (NA x protein crossings, total and
# 4YBB): none held 174/76, 2: 168/70, 3: 152/54, 4: 164/55 (5F9R 8 -> 20), 6: 148/39 (two joints through their pieces)
_HELIX_LOOP = 3


def _strand_steps(items: list[_Item], na: Nucleic):
    """(strand, nucleotide i, next held nucleotide j, (item, stem) of i, (item, stem) of j) for each step of a
    strand from one stem item to another."""
    at = {key: (k, m) for k, it in enumerate(items) for m in it.members if isinstance(m, PlacedStem) for key in m.axial}
    for s, strand in enumerate(na.strands):
        held = [i for i in range(len(strand)) if (s, i) in at]
        for i, j in zip(held, held[1:]):
            if at[(s, i)][0] != at[(s, j)][0]:
                yield s, i, j, at[(s, i)], at[(s, j)]


def _is_joint(n: int, a: PlacedStem, b: PlacedStem) -> bool:
    """A step through n nucleotides from piece a to piece b is a joint of one bent duplex."""
    return n == 0 and bool(a.parent) and a.parent == b.parent


def lying_over(layout: Layout) -> set[frozenset[str]]:
    """The ids of each helix drawn over a sheet and each strand of that sheet (Layout.on_sheet): their boxes overlap
    by design."""
    return {frozenset((h, k)) for h, ids in layout.on_sheet.items() for k in ids}


def joined_pieces(stems: list[PlacedStem], na: Nucleic) -> set[frozenset[str]]:
    """The ids of every two consecutive pieces of one bent duplex (a strand steps from one to the next with no
    nucleotide between): they may meet or overlap on the page (_joints)."""
    at = {key: d for d in stems for key in d.axial}
    return {
        frozenset((d.id, e.id))
        for s, strand in enumerate(na.strands)
        for i in range(len(strand) - 1)
        if (d := at.get((s, i))) is not None
        and (e := at.get((s, i + 1))) is not None
        and d is not e
        and _is_joint(0, d, e)
    }


def _joints(items: list[_Item], na: Nucleic) -> set[tuple[int, int]]:
    """(i, j), i < j: the items of two consecutive pieces of one bent duplex (joined_pieces). _attract keeps no
    MARGIN between them (their boxes may touch, the joint spring holding them end to end), and _separate moves a
    short chain of them together (_rigid_joints): no margin is kept open for an element to fill."""
    item = {m.id: k for k, it in enumerate(items) for m in it.members if isinstance(m, PlacedStem)}
    stems = [m for it in items for m in it.members if isinstance(m, PlacedStem)]
    return {tuple(sorted(item[x] for x in pair)) for pair in joined_pieces(stems, na)}


def _stem_springs(items: list[_Item], na: Nucleic) -> list[tuple[int, np.ndarray, int, np.ndarray, float, float]]:
    """(item, port offset, item, port offset, weight, slack) along each strand from one stem to the next: the
    joint between two pieces of a bent duplex holds them end to end; a strand exchange (at most _SHORT_RUN
    nucleotides between stems of more than two strands), or a side of a short internal loop or bulge between two
    stems of one helix (_short_loops), holds the two stem ends close (slack: its nucleotides' steps); any other
    run pulls its two stems together like a loop, with slack for its nucleotides or, when longer, the distance
    between its ends at the items' home (projected) places."""
    owner = _owner(items)
    short = _short_loops(na)
    out = []
    for s, i, j, (ki, mi), (kj, mj) in _strand_steps(items, na):
        oi = np.subtract(mi.port((s, i), True)[0], (mi.cx, mi.cy)) + owner[mi.id][2]
        oj = np.subtract(mj.port((s, j), False)[0], (mj.cx, mj.cy)) + owner[mj.id][2]
        n = j - i - 1
        mates = frozenset((mi.parent or mi.id, mj.parent or mj.id)) in short
        if _is_joint(n, mi, mj):  # a joint: the duplex carries on
            w, slack = _W_JOINT, _JOINT_SLACK
        elif (n <= _SHORT_RUN and len({*mi.strands, *mj.strands}) > 2) or mates:  # ends held close
            w, slack = _W_SHORT, max(_JOINT_SLACK, n * NT_STEP)
        else:  # pulls only beyond the distance seen between its ends: the stems keep their 3D arrangement
            seen = np.add(items[kj].home, oj) - np.add(items[ki].home, oi)
            w, slack = _W_SEQ, max(_LOOP_SLACK, 0.5 * n * NT_STEP, float(_snap(np.hypot(*seen))))
        out.append((ki, oi, kj, oj, w, slack))
    return out


def _short_loops(na: Nucleic) -> set[frozenset[str]]:
    """The ids of every two stems of one helix (nucleic._helices: the inner continuing the outer) between which
    neither side of the internal loop or bulge holds more than _HELIX_LOOP nucleotides: one helix the drawing keeps
    together (4YBB's h44 decoding site, 2GIS P3a/P3b)."""
    out = set()
    stems = na.structure.stems
    by_helix: dict[int, list] = {}
    for st in stems:
        if st.helix >= 0:
            by_helix.setdefault(st.helix, []).append(st)
    for helix in by_helix.values():
        helix.sort(key=lambda st: (st.a, st.pairs[0].i))  # outermost first
        for a, b in zip(helix, helix[1:]):
            x, y = a.pairs[-1], b.pairs[0]
            if max(y.i - x.i - 1, x.j - y.j - 1) <= _HELIX_LOOP:
                out.add(frozenset((a.id, b.id)))
    return out


def _within_protomers(bundles: list[list[str]], by_id: dict[str, SSE], sector: dict) -> list[list[str]]:
    """Helix bundles (features.helix_bundles) cut where they cross between two protomers: the part of such a
    bundle in one protomer stays a bundle if it still spans two chains; else the pieces of one kinked helix in it
    (consecutive, at most features.KINK_GAP residues apart) stay together, and its other helices are drawn alone. A
    bundle of a dimer's two copies of one helix (1LMB α5) would hold them in one item, side by side across their
    axis, where no mirror line can pass between them; apart, each copy can sit at its mate's mirror image. A
    bundle of a ring's copies (the measles F stalk, 8UUP) keeps each copy's part with its protomer, drawn alike. A
    bundle within one protomer, or of one protomer and chains outside every protomer (sector -1: a bound
    peptide), is kept whole."""
    out = []
    for g in bundles:
        parts: dict[int, list[str]] = {}
        for k in g:
            parts.setdefault(sector.get(by_id[k].chain, -1), []).append(k)
        if len(parts.keys() - {-1}) < 2:
            out.append(g)
            continue
        for p in sorted(parts):
            if len({by_id[k].chain for k in parts[p]}) > 1:
                out.append(parts[p])
                continue
            runs: list[list[str]] = []
            for k in sorted(parts[p], key=lambda k: by_id[k].start):
                if runs and by_id[k].start - by_id[runs[-1][-1]].end - 1 <= KINK_GAP:
                    runs[-1].append(k)
                else:
                    runs.append([k])
            out += [r for r in runs if len(r) > 1]
    return out


def _mirrored(angle: float) -> float:
    """An element's angle seen in a vertical mirror (pi - angle), in (-pi, pi]."""
    return float(-((angle - 2 * np.pi) % (2 * np.pi) - np.pi)) + 0.0  # + 0.0: never -0.0


def _ring(frame) -> int | None:
    """The number of protomers round a closed ring (a cyclic assembly, or a dimer drawn side by side), counted from
    the frame's sectors; None along a helical filament (its protomers rise and never come round)."""
    sector = getattr(frame, "sector", None) or {}
    if getattr(frame, "rise", 0.0):
        return None
    return len(set(sector.values())) or None


def _host(members: list[Placed], sector: dict, classes: dict) -> int | None:
    """The protomer an item holding elements of two or more protomers belongs to: the one holding most of its
    residues (a copy's own strands with strands lent by its neighbour, as in a measles F sheet; each copy's such sheet
    is then its mates' image); on a tie, that of its first element by chain class, then residue number. None when
    two protomers' copies of one element come first (an item symmetric in itself)."""
    real = sorted(
        (classes[m.sse.chain], m.sse.first.seq, m.sse.last.seq, sector[m.sse.chain]) for m in members if not m.ghost
    )
    held: dict[int, int] = {}
    for _, a, b, k in real:
        held[k] = held.get(k, 0) + b - a + 1
    real = [r for r in real if held[r[3]] == max(held.values(), default=0)]
    if not real:
        return None
    c, a, b, k = real[0]
    if any(x == c and y <= b and a <= z and s != k for x, y, z, s in real[1:]):
        return None
    return k


def _protomers(items: list[_Item], frame, mirror: bool) -> list[int | None]:
    """The protomer of each item (frame.sector): that of all its elements, or of an item holding several protomers'
    elements its _host, round a closed ring (_ring) or along a filament (a FimA subunit's sheet with the next
    one's donor strand; not with a mirror line, where such items stay as they were: untied); None for an item with
    chains outside every protomer."""
    sector, classes = getattr(frame, "sector", None) or {}, getattr(frame, "classes", None) or {}
    out = []
    for it in items:
        members = [m for m in it.members if isinstance(m, Placed)]
        ks = {sector.get(m.sse.chain) for m in members}
        if len(ks) == 1:
            out.append(ks.pop())
        else:
            out.append(None if mirror or None in ks or not members else _host(members, sector, classes))
    return out


def _symmetry_ties(
    items: list[_Item], frame, mirror: bool = True, shift: dict[str, int] | None = None
) -> dict[int, tuple[int, int]]:
    """Copy item -> (its protomer-0 mate, protomer index k): matched by chain class and shared residue numbers.
    The copy also takes the mate's shape and orientation, so every protomer is drawn the same; in a frame with a
    mirror line (`mirror_x`, a dyad in the page) it takes them mirrored, protomer 1 drawn as protomer 0's mirror
    image; `mirror` False takes them unmirrored all the same (build_layout drops a mirror line a βα row crosses).
    An item holding elements of several protomers (_protomers: a sheet with a neighbour's strand) is matched with
    each element's protomer counted from the item's own, round the ring (_ring). A copy whose elements each match
    one of its mate's, but not all of them (DSSP missed a strand in one copy), is tied too when they share at least
    half the smaller one's residues; and when no element of a copy matches two of its mate's, those it shares take
    the mate's shape and place, the rest keep theirs beside them. Residue numbers are taken to the first copy's by
    `shift` (chain -> offset, model.homolog_chains: 2MJW numbers chain D 95-188 where B is 1-94)."""
    shift = shift or {}
    sector, classes = getattr(frame, "sector", None), getattr(frame, "classes", None)
    if not sector or not classes:
        return {}
    mirror = mirror and getattr(frame, "mirror_x", None) is not None
    ring = _ring(frame)
    homes = _protomers(items, frame, mirror)

    def span(m, k):  # its residues by chain class and protomer counted from the item's own, k
        if not isinstance(m, Placed) or m.sse.chain not in classes:
            return set()
        d = sector[m.sse.chain] - k
        place = (classes[m.sse.chain], d % ring if ring else d)
        off = shift.get(m.sse.chain, 0)
        return {(place, r - off) for r in range(m.sse.first.seq, m.sse.last.seq + 1)}

    def residues(i):  # chains outside every protomer (bound peptides, unpaired chains) take no part
        return set().union(*(span(m, homes[i]) for m in items[i].members if isinstance(m, Placed) and not m.ghost))

    sets = [residues(i) if homes[i] is not None else set() for i in range(len(items))]
    zero = [i for i in range(len(items)) if homes[i] == 0 and sets[i]]
    ties = {}
    for i, it in enumerate(items):
        k = homes[i]
        if not k or not sets[i]:
            continue
        best = max(zero, key=lambda j: len(sets[i] & sets[j]) / len(sets[i] | sets[j]), default=None)
        if best is None:
            continue
        if len(sets[i] & sets[best]) < 0.5 * len(sets[i] | sets[best]):
            best = _part_mate(it, k, zero, items, sets[i], sets, span)
            if best is None:
                continue
        mate = items[best]
        ghosts = [q for q, x in enumerate(mate.members) if x.ghost]
        pairs = [
            ghosts[0]  # a barrel's ghost takes the place of the mate's, whichever strand each repeats
            if m.ghost and len(ghosts) == 1
            else max(
                range(len(mate.members)),
                key=lambda q: (len(span(m, k) & span(mate.members[q], 0)), mate.members[q].ghost == m.ghost),
            )
            for m in it.members
        ]
        if len(mate.members) == len(it.members) and sorted(pairs) == list(range(len(pairs))):
            for m, q in zip(it.members, pairs):
                if m.ghost:  # it repeats the strand that pairs with the one the mate's ghost repeats (a barrel seen
                    # from its other face opened its row at the other end)
                    m.sse = next(
                        (
                            x.sse
                            for x, r in zip(it.members, pairs)
                            if not x.ghost and mate.members[r].sse is mate.members[q].sse
                        ),
                        m.sse,
                    )
            for m, q in zip(it.members, pairs):  # one-to-one: take the mate's shape member by member
                a = mate.members[q].angle
                m.angle, m.length = (_mirrored(a) if mirror else a), mate.members[q].length
            it.offsets = [(-mate.offsets[q][0], mate.offsets[q][1]) if mirror else mate.offsets[q] for q in pairs]
        else:  # DSSP split, merged or missed an element in one copy: match the orientation
            shared = [(p, q) for p, q in enumerate(pairs) if span(it.members[p], k) & span(mate.members[q], 0)]
            votes = [
                np.cos(it.members[p].angle - (np.pi - mate.members[q].angle if mirror else mate.members[q].angle))
                for p, q in shared
            ]
            if votes and np.mean(votes) < 0:  # drawn upside down: turn the block 180 degrees in the page
                for m in it.members:
                    m.angle += np.pi
                it.offsets = [(-ox, -oy) for ox, oy in it.offsets]
            if not mirror and shared and len({q for _, q in shared}) == len(shared):
                _take_shared(it, mate, shared)
        ties[i] = (best, k)
    return ties


def _part_mate(it: _Item, k: int, zero: list[int], items: list[_Item], mine: set, sets: list[set], span) -> int | None:
    """The protomer-0 item whose elements hold each of this copy's elements, one each, or the other way round (DSSP
    missed a strand in one copy), sharing at least half the smaller one's residues; the most shared first."""

    def into(a: list, ka: int, b: list, kb: int) -> bool:  # every element of a overlaps its own element of b
        hit = [next((q for q, y in enumerate(b) if span(x, ka) & span(y, kb)), None) for x in a if span(x, ka)]
        return None not in hit and len(set(hit)) == len(hit)

    mine_m = [m for m in it.members if isinstance(m, Placed) and not m.ghost]
    found = []
    for j in zero:
        both = len(mine & sets[j])
        if not both or both < 0.5 * min(len(mine), len(sets[j])):
            continue
        theirs = [m for m in items[j].members if isinstance(m, Placed) and not m.ghost]
        if into(mine_m, k, theirs, 0) or into(theirs, 0, mine_m, k):
            found.append((-both, j))
    return min(found)[1] if found else None


def _take_shared(it: _Item, mate: _Item, shared: list[tuple[int, int]]) -> None:
    """A copy's elements `shared` (member, its mate's member) take the mate's shape and place; the copy's other
    elements move with them, by their mean step, so each keeps its place beside the elements it is drawn with."""
    step = np.mean([np.subtract(mate.offsets[q], it.offsets[p]) for p, q in shared], axis=0)
    offsets = [tuple(float(v) for v in _snap(np.add(o, step))) for o in it.offsets]
    for p, q in shared:
        it.members[p].angle, it.members[p].length = mate.members[q].angle, mate.members[q].length
        offsets[p] = mate.offsets[q]
    it.offsets = offsets


def _unwrap(items: list[_Item], ties: dict[int, tuple[int, int]], frame) -> dict[int, int]:
    """Elements lent across the cut of a ring unrolled into a row. A copy's element drawn in its neighbour's item
    (a strand that pairs in the next protomer's sheet) sits on the same side of its own protomer in every copy but
    one: the copy whose neighbour lies across the cut, at the row's other end. There the element leaves the item
    for an item of its own, tied to the protomer-0 mate's place the same number of protomers along as in the other
    copies (one past the row's end), so every copy is drawn alike; a strand leaves a ghost (dashed, as a barrel's)
    where it pairs, and takes ghosts of the strands it pairs with (_partners) to stand beside. New items are
    appended to `items` and their ties to `ties`. Returns each new item's index -> the index of the item it left."""
    ring = _ring(frame)
    sector = getattr(frame, "sector", None) or {}
    if not ring or not ties:
        return {}
    homes = _protomers(items, frame, False)
    roots = {j for j, _ in ties.values()}
    orbit = sorted([(j, j, 0) for j in roots] + [(i, j, k) for i, (j, k) in ties.items() if homes[i] == k])
    out = {}
    for x, root, h in orbit:
        it = items[x]
        away: dict[int, list[int]] = {}
        for p, m in enumerate(it.members):
            if not isinstance(m, Placed) or m.ghost or m.sse.chain not in sector or sector[m.sse.chain] == h:
                continue
            g = sector[m.sse.chain]
            r = (g - h) % ring
            usual = -r if 2 * r < ring else ring - r  # protomers from the element's own to the item holding it
            if h - g != usual:
                away.setdefault(g + usual, []).append(p)
        for k, ps in sorted(away.items()):
            moved = [it.members[p] for p in ps]
            partners = _partners(it, ps)
            items.append(
                _Item(
                    moved + [Placed(m.sse, 0, 0, m.length, m.angle, m.width, ghost=True) for m in partners],
                    [it.offsets[p] for p in ps] + [it.offsets[it.members.index(m)] for m in partners],
                    it.cx,
                    it.cy,
                    it.home,
                    *_anchor([m.sse for m in moved]),
                )
            )
            ties[len(items) - 1] = (root, k)
            out[len(items) - 1] = x
            for p in ps:  # a strand leaves its ghost; a helix just leaves (and lies over no sheet now)
                m = it.members[p]
                m.over = False
                ghost = Placed(m.sse, 0, 0, m.length, m.angle, m.width, ghost=True)
                it.members[p] = ghost if m.sse.kind == "E" else None
            keep = [p for p, m in enumerate(it.members) if m is not None]
            it.members, it.offsets = [it.members[p] for p in keep], [it.offsets[p] for p in keep]
    return out


def _partners(it: _Item, ps: list[int]) -> list[Placed]:
    """The strands of item `it` beside its strands `ps` in the sheet (one strand pitch from one of them, centre to
    centre), not among them: those they pair with, whose ghosts stand beside them where they leave (_unwrap)."""
    out = []
    for q, m in enumerate(it.members):
        if q in ps or not isinstance(m, Placed) or m.ghost or m.sse.kind != "E":
            continue
        if any(
            it.members[p].sse.kind == "E"
            and np.hypot(*np.subtract(it.offsets[q], it.offsets[p])) <= max(m.width, it.members[p].width) + 1e-6
            for p in ps
        ):
            out.append(m)
    return out


def _copy_ties(
    items: list[_Item], ties: dict[int, tuple[int, int]], frame, mirrored: bool
) -> dict[str, tuple[str, bool]]:
    """Chain -> (the protomer-0 chain of its class, mirrored) for each chain of protomer k holding an element of a
    copy item of `ties` (_symmetry_ties) tied as protomer k: the copies the symmetry asked to draw alike, whether or
    not _separate keeps them tied (an element lent by a neighbouring protomer is that protomer's)."""
    sector, classes = getattr(frame, "sector", None) or {}, getattr(frame, "classes", None) or {}
    first = {classes[c]: c for c in sorted(classes) if sector.get(c) == 0}
    out = {}
    for i in sorted(ties):
        for m in items[i].members:
            c = m.sse.chain if isinstance(m, Placed) and not m.ghost else None
            if c in classes and classes[c] in first and sector.get(c) == ties[i][1]:
                out[c] = (first[classes[c]], mirrored)
    return dict(sorted(out.items()))


_RES_RANGE = re.compile(r"(-?\d+)([A-Za-z]?)(?:-(-?\d+)([A-Za-z]?))?")  # 150, 150-159, -7-4, 100A-100K


def resolve(layout: Layout, sses: list[SSE], ref: str) -> list[str]:
    """Element ids for a reference: a label (G, α2; every chain), chain:label (A:G), #k (the k-th element
    from the N terminus, 1-based), or res:a-b / res:chain:a-b (elements overlapping those residues; negative
    numbers and insertion codes as in the file, e.g. res:A:-7-4, res:H:100A-100K)."""
    ref = ref.strip()
    placed = [s for s in sses if s.id in layout.placed]
    if ref.startswith("#"):
        k = int(ref[1:]) if ref[1:].isdigit() else 0
        if not 1 <= k <= len(placed):
            raise ValueError(f"no element {ref!r}: elements are numbered #1 to #{len(placed)}")
        return [placed[k - 1].id]
    if ref.startswith("res:"):
        body = ref[4:]
        chain = None
        if body.count(":") == 1:
            chain, body = body.split(":")
        m = _RES_RANGE.fullmatch(body)
        if m is None:
            raise ValueError(
                f"residue reference {ref!r} must look like res:150-159, res:A:150-159, res:A:-7-4 or res:H:100A-100K"
            )
        lo = (int(m[1]), m[2])
        hi = (int(m[3]), m[4]) if m[3] is not None else lo
        in_chain = [s for s in placed if chain is None or s.chain == chain]
        # An element's own reference first: descending numbering (1H ... 1A, 1) has first > last.
        out = [s.id for s in in_chain if (s.first.tag, s.last.tag) == (f"{lo[0]}{lo[1]}", f"{hi[0]}{hi[1]}")]
        if not out:
            lo, hi = min(lo, hi), max(lo, hi)
            out = [
                s.id
                for s in in_chain
                if min((s.first.seq, s.first.icode), (s.last.seq, s.last.icode)) <= hi
                and max((s.first.seq, s.first.icode), (s.last.seq, s.last.icode)) >= lo
            ]
        if not out and not lo[1] and not hi[1]:
            # Saved by foldmap 0.1.0, which ignored insertion codes: match by residue number alone.
            out = [
                s.id
                for s in in_chain
                if min(s.first.seq, s.last.seq) <= hi[0] and max(s.first.seq, s.last.seq) >= lo[0]
            ]
    elif ":" in ref:
        chain, label = ref.split(":", 1)
        out = [s.id for s in placed if s.chain == chain and layout.placed[s.id].label == label]
    else:
        out = [s.id for s in placed if layout.placed[s.id].label == ref]
    if not out:
        raise ValueError(f"no element {ref!r}; use a label such as A or α1, chain:label, #k or res:a-b")
    return out


def _one(layout: Layout, sses: list[SSE], ref: str) -> str:
    ids = resolve(layout, sses, ref)
    if len(ids) != 1:
        raise ValueError(f"{ref!r} matches {len(ids)} elements ({', '.join(ids)}); name one, e.g. chain:label")
    return ids[0]


def _adjust(
    items: list[_Item],
    lay: Layout,
    sses: list[SSE],
    rename: dict[str, str],
    swap: list[tuple[str, str]],
    move: list[tuple[str, tuple[float, float]]],
    joints: set[tuple[int, int]] = frozenset(),
    fixed: Walls | None = None,
    level: np.ndarray | None = None,
) -> None:
    """The user's own hand on the layout: new labels, exchanged places, nudges; then overlaps are cleared
    (consecutive pieces of a bent duplex, `joints`, may still meet; items keep off the `fixed` boxes; the items of
    `level`, a membrane view's transmembrane elements, keep their height: _separate)."""
    where = {
        m.sse.id: (it, k) for it in items for k, m in enumerate(it.members) if isinstance(m, Placed) and not m.ghost
    }
    for ref, name in rename.items():
        for sid in resolve(lay, sses, ref):
            it, k = where[sid]
            it.members[k].label = name
            for g in it.members:
                if g.ghost and g.sse.id == sid:
                    g.label = name
    moved = False
    for a_ref, b_ref in swap:
        a, b = _one(lay, sses, a_ref), _one(lay, sses, b_ref)
        (ia, ka), (ib, kb) = where[a], where[b]
        if ia is ib:  # two strands of one sheet: exchange their slots
            ia.offsets[ka], ia.offsets[kb] = ia.offsets[kb], ia.offsets[ka]
        else:
            ia.cx, ib.cx, ia.cy, ib.cy = ib.cx, ia.cx, ib.cy, ia.cy
        moved = True
    for ref, (dx, dy) in move:
        for sid in resolve(lay, sses, ref):
            it, _ = where[sid]
            it.cx += dx
            it.cy += dy
            moved = True
    if moved:
        for it in items:
            it.home = (it.cx, it.cy)
        _separate(items, None, joints, fixed, level=level)
        for it in items:
            it.apply()


def _tm_ends(it: _Item, level: set[str], y: float) -> list[tuple[float, float]]:
    """Page y of the two ends (ports) of each transmembrane element (`level`) of an item centred at height y."""
    out = []
    for m, (_, oy) in zip(it.members, it.offsets):
        if getattr(m, "sse", None) is not None and m.sse.id in level and not m.ghost:
            reach = abs(float(np.sin(m.angle))) * m.length / 2
            out.append((y + oy - reach, y + oy + reach))
    return out


def _to_face(it: _Item, level: set[str], band: tuple[float, float], y: float) -> float:
    """The least move up or down that brings each transmembrane element of the item at height y drawn wholly
    between the faces (a tilted helix, drawn upright at its length, is shorter than the slab it crosses) to a face;
    0 when none is."""
    inside = [(y0, y1) for y0, y1 in _tm_ends(it, level, y) if band[0] < y0 and y1 < band[1]]
    if not inside:
        return 0.0
    up = max(band[1] - y1 for _, y1 in inside)
    down = min(band[0] - y0 for y0, _ in inside)
    return up if up < -down else down


def _text_scaled(build):
    """build_layout with _GLYPH_K set, while it runs, to the glyph_scale of the text its style draws (Style.font_scale):
    the glyph room of footprints and stem names follows the text size, as layout.glyph does where a layout is held."""
    import functools

    @functools.wraps(build)
    def run(*args, **kwargs):
        style = kwargs.get("style") or Style()
        token = _GLYPH_K.set(glyph_scale(_EM_TUNED * style.font_scale))
        try:
            return build(*args, **kwargs)
        finally:
            _GLYPH_K.reset(token)

    return run


@_text_scaled
def build_layout(
    sses: list[SSE],
    sheets: list[Sheet],
    frame: Frame,
    mode: str = "projected",
    *,
    attract: bool = True,
    contacts: dict[frozenset[str], int] | None = None,
    nucleic: Nucleic | None = None,
    style: Style | None = None,
    bridges: list[tuple[str, float, str, float]] | None = None,
    domains: list[tuple[str, list[str]]] | None = None,
    partial: set | None = None,
    rename: dict[str, str] | None = None,
    swap: list[tuple[str, str]] | None = None,
    move: list[tuple[str, tuple[float, float]]] | None = None,
    keep_depth: bool = False,
    na_cropped: bool = False,
    template: str | None = None,
    copies: dict[str, tuple[str, int]] | None = None,
    relabel: dict[str, str] | None = None,
    names: Callable[[Layout, list[SSE]], None] | None = None,
    over: dict[str, int] | None = None,
    level: set[str] | None = None,
    band: tuple[float, float] | None = None,
) -> Layout:
    """Labels: the generic ones (_assign_labels), then `relabel` (ref -> label: copies' shared labels, primes of
    elements lent by neighbours), then `names` (names.apply_names: --names, given the placed layout with its stems
    and ba_rows), then `rename` (the user's own, through _adjust like `swap` and `move`); the label sides and boxes
    are chosen after (_choose_label_sides).
    `template` (Style.fold_template, the style's when None): "auto" draws each sheet that template_applies to as a
    βα row (_ba_template), in projected mode; `copies` (model.homolog_chains) makes copies of a protein share that
    choice (template_choices). `over` (features.helices_on_sheets: helix id -> sheet index): each helix drawn over
    its sheet's strands (_sheet_item, _over_sheet), in projected mode, unless the sheet is a βα row, the helix is a
    row's, a bundle's or a piece of a kinked helix (a rigid unit drawn whole beside the sheet), or a copy of it in
    another chain of the protein (`copies`) is drawn but not over its sheet (_over_alike). `keep_depth` (stack mode
    with a membrane): each element stays at its height along the frame's up. `na_cropped` (na_extent auto cropped
    the nucleic acid): drawn as blocks, without base letters. `attract` (projected mode) pulls loop-joined and
    touching elements together before overlaps are removed; `contacts` comes from features.sse_contacts. `nucleic`
    adds its stems (straight pieces of base-paired DNA/RNA), placed like any element and pulling the elements that
    touch them; strands with no pair at all are kept as projected traces. `level` (projected mode; the
    transmembrane elements of a membrane view, cli): the items holding them keep their projected height, so they
    cross the band where the membrane is; the rest moves round them (_separate). `band` (page y of the membrane's
    two faces as projected): a lone transmembrane helix drawn wholly between them (a tilted helix drawn upright at
    its length is shorter than the slab it crosses) is moved, the least way, until one end meets a face."""
    if mode not in ("projected", "stack"):
        raise ValueError(f"unknown layout mode {mode!r}")
    style = style or Style()
    angle_mode = "upright" if mode == "stack" else style.helix_angle
    template = style.fold_template if template is None else template
    if template not in CHOICES["fold_template"]:
        raise ValueError(f"fold_template must be one of {', '.join(CHOICES['fold_template'])}; got {template!r}")
    as_rows = mode == "projected" and template == "auto"  # a stack row draws elements in sequence instead
    rowed = template_choices(sheets, sses, copies or {}) if as_rows else [False] * len(sheets)
    templates = {
        k: _ba_template(
            ring := _open_at_ends(sh),
            _ba_units(ring, sses),
            frame,
            PITCH * style.strand_scale,
            HELIX_W * style.helix_scale,
        )
        for k, (sh, row) in enumerate(zip(sheets, rowed))
        if row
    }
    in_rows = {m.sse.id for it in templates.values() for m in it.members if m.sse.kind == "H"}  # helices of βα rows
    mirror = getattr(frame, "mirror_x", None)  # a dyad in the page: protomer 1 drawn as protomer 0's mirror image
    rows_x = [(it.cx + it.box[0], it.cx + it.box[2]) for it in templates.values()]
    if mirror is not None and any(x0 < mirror < x1 for x0, x1 in rows_x):
        mirror = None  # a βα row across the mirror line would meet its own image: the copies go side by side instead
    by_id = {s.id: s for s in sses}
    bundles = helix_bundles([s for s in sses if s.id not in in_rows])  # bundles are coiled coils, never a row's
    if getattr(frame, "sector", None):  # each protomer's part of a bundle drawn with its own copy
        bundles = _within_protomers(bundles, by_id, frame.sector)
    bundled = {k for g in bundles for k in g}
    kinked = (
        _kink_groups(sses, bundled | in_rows) if mode == "projected" else []
    )  # a stack row draws each piece on its own
    on_sheet = (
        _over_alike(over or {}, sses, rowed, bundled | in_rows | {s.id for g in kinked for s in g}, copies or {})
        if mode == "projected"
        else {}
    )
    items = [
        templates[k]
        if k in templates
        else _sheet_item(
            sh,
            frame,
            PITCH * style.strand_scale,
            [by_id[h] for h, at in sorted(on_sheet.items(), key=lambda x: by_id[x[0]].start) if at == k],
            HELIX_W * style.helix_scale,
        )
        for k, sh in enumerate(sheets)
    ]
    items += [_bundle_item([by_id[k] for k in g], frame, angle_mode, HELIX_W * style.helix_scale) for g in bundles]
    items += [_kink_item(g, frame, angle_mode, HELIX_W * style.helix_scale) for g in kinked]
    lying = {m.sse.id for it in items for m in it.members if m.over}  # over a sheet: on_sheet less any crowded out
    grouped = bundled | in_rows | {s.id for g in kinked for s in g} | lying
    items += [
        _helix_item(s, frame, angle_mode, HELIX_W * style.helix_scale)
        for s in sses
        if s.kind == "H" and s.id not in grouped
    ]
    items += [
        _helix_item(s, frame, angle_mode, HELIX_W * ETA_W * style.helix_scale, TEN_RISE) for s in sses if s.kind == "G"
    ]
    contacts = dict(contacts or {})
    pieces = stem_pieces(nucleic) if nucleic is not None and nucleic.strands else []
    stem_items = [_stem_item(p, nucleic, frame, mode == "stack", NA_W * style.dna_scale) for p in pieces]
    if stem_items:
        _continuous_sides([it.members[0] for it in stem_items], nucleic)
    items += stem_items
    for key, n in (_stem_contacts(sses, nucleic, pieces) if pieces else {}).items():
        contacts[key] = contacts.get(key, 0) + n
    traces = _traces(nucleic, frame) if nucleic is not None and nucleic.strands else {}
    if not items and not traces:
        return Layout({}, (0.0, 0.0, 0.0, 0.0))
    items.sort(key=lambda it: it.order)
    members = {m.sse.id: m for it in items for m in it.members if isinstance(m, Placed) and not m.ghost}
    for chain, mine in group_by_chain(sses, keep=members).items():
        if chain not in (partial or ()):  # fragments of neighbours have no termini of their own
            members[mine[0].id].ends += "N"
            members[mine[-1].id].ends += "C"
    joints = _joints(items, nucleic) if pieces else set()
    walls = _trace_walls(items, traces, nucleic) if traces else None  # traces stay: items keep off them
    rows: list[list[_Item]] = []
    kept = None  # items keeping their height (projected mode, `level`)
    mirrored: list[tuple] = []  # (copy, mate) members of protomers drawn as mirror images
    copy_ties: dict[str, tuple[str, bool]] = {}
    tied_ids: set[str] = set()
    if mode == "stack":
        rows = _stack(items, keep_depth)
        _centre_dna_over_binders(items, contacts)
    else:
        ties = _symmetry_ties(items, frame, mirror is not None, {c: off for c, (_, off) in (copies or {}).items()})
        left = _unwrap(items, ties, frame) if mirror is None else {}
        if left and walls is not None:  # an element that left its item keeps off the traces it kept off
            walls = Walls(walls.boxes, np.vstack([walls.hard, walls.hard[list(left.values())]]))
        copy_ties = _copy_ties(items, ties, frame, mirror is not None)
        tied_ids = {m.sse.id for i, (j, _) in ties.items() for k in (i, j) for m in items[k].members}
        spread, lift = 0.0, 0.0
        if ties and mirror is not None:  # a dyad in the page: each copy starts at its mate's mirror image
            mirrored = [
                (a, b) for i, (j, _) in sorted(ties.items()) for a, b in zip(items[i].members, items[j].members)
            ]
            for i, (j, _) in ties.items():
                items[i].cx, items[i].cy = 2 * mirror - items[j].cx, items[j].cy
                items[i].home = (2 * mirror - items[j].home[0], items[j].home[1])
        elif ties:  # protomers side by side at least one protomer's width apart; copies start on their mates
            mates = {j for j, _ in ties.values()}
            x0 = min(items[j].cx + items[j].box[0] for j in mates)
            x1 = max(items[j].cx + items[j].box[2] for j in mates)
            # protomer k sits k spreads along: rightwards when the copies turn anticlockwise seen from the top
            # of the axis, leftwards for a left-handed filament (or an assembly turned over), as projected
            period = getattr(frame, "period", 0.0)  # a CylinderFrame's; a Frame's (its mirror dropped) is 0
            spread = float(np.sign(period) or 1.0) * max(abs(period) * SCALE, x1 - x0 + 2 * MARGIN)
            lift = getattr(frame, "rise", 0.0) * SCALE
            for i, (j, k) in ties.items():
                items[i].cx, items[i].cy = items[j].cx + k * spread, items[j].cy + k * lift
                items[i].home = (items[j].home[0] + k * spread, items[j].home[1] + k * lift)
        kept = (
            np.array(
                [any(getattr(m, "sse", None) is not None and m.sse.id in level for m in it.members) for it in items]
            )
            if level
            else None
        )
        heights = [it.home[1] for it in items]
        if attract and len(items) > 1:
            extra = _stem_springs(items, nucleic) if pieces else []
            rise = getattr(frame, "rise", 0.0) * SCALE
            _attract(
                items,
                sses,
                contacts,
                ties,
                spread,
                rise,
                bridges or [],
                extra=extra,
                joints=joints,
                fixed=walls,
                mirror=mirror if ties else None,
            )
        if kept is not None and kept.any():  # transmembrane elements at their projected height, reaching a face
            for k in np.flatnonzero(kept).tolist():
                shift = _to_face(items[k], level, band, heights[k]) if band is not None else 0.0
                y = float(_snap(heights[k] + shift))
                items[k].cy = y
                items[k].home = (items[k].home[0], y)
        _separate(items, ties, joints, walls, mirror=mirror if ties else None, level=kept)
        if domains and _part_domains(items, domains, ties, mirror if ties else None):
            # the domains, parted with their copies, keep them alike whatever room that takes
            _separate(items, ties, joints, walls, mirror=mirror if ties else None, room=np.inf, level=kept)
    for it in items:
        it.apply()
    dna = [m for it in items for m in it.members if isinstance(m, PlacedStem)]
    placed = [m for it in items for m in it.members if isinstance(m, Placed) and not m.ghost]
    ghosts = [m for it in items for m in it.members if isinstance(m, Placed) and m.ghost]
    _assign_labels(placed)
    real = {p.sse.id: p for p in placed}
    named = Layout(real, (0.0, 0.0, 0.0, 0.0), stems=dna, nucleic=nucleic if dna else None)
    named.ba_rows = _ba_rows(items)
    if relabel:  # labels shared by copies, primes of lent elements: before any name
        _adjust(items, named, sses, relabel, [], [])
    if names is not None:  # conventional names (names.apply_names): before the user's own
        names(named, sses)
    if rename or swap or move:
        _adjust(items, named, sses, rename or {}, swap or [], move or [], joints, walls, kept)
    for g in ghosts:
        g.label = real[g.sse.id].label
    rects = np.array([p.rect for p in placed + ghosts + dna]).reshape(-1, 4)
    blocks = [  # the strands of each item holding any: a βα row's helices lie between them, on its panel
        (ids, [m.sse.id for m in it.members if isinstance(m, Placed) and m.ghost])
        for it in items
        if (ids := [m.sse.id for m in it.members if isinstance(m, Placed) and m.sse.kind == "E" and not m.ghost])
    ]
    lay = Layout(
        real,
        (0.0, 0.0, 0.0, 0.0),
        ghosts=ghosts,
        stems=dna,
        nucleic=nucleic if dna or traces else None,
        sheets=[ids for ids, _ in blocks],
        sheet_ghosts=[ghost_ids for _, ghost_ids in blocks],
    )
    lay.na_traces = traces
    lay.frame = frame
    lay.copy_ties = copy_ties
    lay.tied_ids = tied_ids
    lay.mirror_pairs = [(a.sse.id, b.sse.id) for a, b in mirrored if isinstance(a, Placed) and a.sse.kind in ("H", "G")]
    lay.font_scale = style.font_scale
    lay.ligands = style.ligands
    lay.ba_rows = _ba_rows(items)  # after the user's swaps
    lay.stack_gaps = _stack_gaps(lay)
    lay.on_sheet = {
        m.sse.id: [x.sse.id for x in it.members if x.sse.kind == "E" and not x.ghost]
        for it in items
        for m in it.members
        if isinstance(m, Placed) and m.over
    }
    if keep_depth and len(rows) > 1:  # wrapped into rows of whole chains: each row gets its own band
        lay.stack_rows = [[m.sse.id for it in row for m in it.members if isinstance(m, Placed)] for row in rows]
    if traces:
        pts = np.vstack(list(traces.values()))
        rects = np.vstack([rects, [[pts[:, 0].min(), pts[:, 1].min(), pts[:, 0].max(), pts[:, 1].max()]]])
    lay.partial = set(partial or ())
    chain_of = {k: p.sse.chain for k, p in lay.placed.items()}
    lay.domains = [
        (name, [k for k in ids if k in real and chain_of[k] not in lay.partial]) for name, ids in domains or []
    ]
    lent = {name: [k for k in ids if k in real and chain_of[k] in lay.partial] for name, ids in domains or []}
    lay.domain_lent = {name: ids for name, ids in lent.items() if ids}
    rects = np.vstack([rects, *[[domain_panel(lay, name)] for name, _ in lay.domains]]) if lay.domains else rects
    lay.na_cropped = na_cropped
    _choose_label_sides(lay, sses, style, mirrored)  # after naming: a name's length sets its box
    _place_domain_names(lay, sses, style.font_scale)  # after the element and stem names it keeps clear of
    rects = np.vstack([rects, [b for _, b in label_boxes(lay, sses)] or np.empty((0, 4))])
    lay.bounds = (
        float(rects[:, 0].min()),
        float(rects[:, 1].min()),
        float(rects[:, 2].max()),
        float(rects[:, 3].max()),
    )
    return lay


def _ba_rows(items: list[_Item]) -> list[list[str]]:
    """Each βα row's elements (no ghost) as drawn, left to right: the items holding strands and helices in line
    (a helix over a sheet is not in a row)."""
    return [
        [m.sse.id for m in sorted((m for m in it.members if not m.ghost), key=lambda m: m.cx)]
        for it in items
        if any(isinstance(m, Placed) and m.sse.kind == "E" for m in it.members)
        and any(isinstance(m, Placed) and m.sse.kind == "H" and not m.over for m in it.members)
    ]


def _over_alike(
    over: dict[str, int], sses: list[SSE], rowed: list[bool], taken: set[str], copies: dict[str, tuple[str, int]]
) -> dict[str, int]:
    """The helices of `over` (helix id -> sheet index) that are drawn over their sheet: not on a βα row's sheet
    (`rowed`, template_choices), not `taken` by another unit (a row, a bundle, a kinked helix), and alike in every
    copy of the protein: a helix whose copy in another chain (`copies`, model.homolog_chains; a helix at its place,
    model.same_place, once numbered as the first copy's) is drawn but not over its sheet stays beside, as that copy
    does."""
    by_id = {s.id: s for s in sses}
    keep = {h: k for h, k in sorted(over.items()) if h in by_id and not rowed[k] and h not in taken}

    def key(s):
        root, shift = copies.get(s.chain, (s.chain, 0))
        return root, s.first.seq - shift, s.last.seq - shift

    helices = [(s, key(s)) for s in sses if s.kind == "H"]
    out = {}
    for h, k in keep.items():
        root, a, b = key(by_id[h])
        mates = [
            g for g, (r, ga, gb) in helices if g.chain != by_id[h].chain and r == root and same_place((a, b), (ga, gb))
        ]
        if all(g.id in keep for g in mates):
            out[h] = k
    return out


DOMAIN_PAD = 0.7  # panel margin around a domain's elements
_DOMAIN_GAP = 2 * DOMAIN_PAD + 0.6  # clear space between two domains' elements: both panels plus a gutter


_PART_BLOWUP = 2.0  # parting that left domains unparted and the items' width + height this many times: undone


def _part_domains(
    items: list[_Item],
    domains: list[tuple[str, list[str]]],
    ties: dict[int, tuple[int, int]] | None = None,
    mirror: float | None = None,
    sweeps: int = 60,
) -> bool:
    """Move whole domains apart (each as a rigid group) until their boxes no longer overlap. Each item goes with
    the domain most of its elements belong to. With `ties` (_symmetry_ties) copies stay alike: a domain moves with
    its copies' domains and every item of their orbits (a unit, _domain_units), each copy by its tie (a mirrored
    copy, `mirror`, the other way across); two domains of one unit, copies of each other, part by widening the
    spread between protomers (protomer k moves k times the step), or for a mirror image sideways symmetrically.
    When the sweeps do not part them and leave the items' width plus height over _PART_BLOWUP times (3L9G: a sheet
    shared by four copies holds ten domains, so every sweep widened the spread again: random-100 C31), they start
    over with domains that share an item never parted (nor any domain joined to them so), stopping when a sweep no
    longer shrinks the overlap. True if anything moved. A domain moves as a whole, up or down too, its transmembrane
    elements in a membrane view with it (build_layout `level` is not held here: a domain keeps its own shape, and
    holding some of its items' heights would tear it; the _separate after parting holds them again)."""
    ties = ties or {}
    owner = {k: d for d, (_, ids) in enumerate(domains) for k in ids}
    groups: dict[int, list[int]] = {}
    bound = list(range(len(domains)))  # domains joined (sharing an item), on the second try only
    shared: list[list[int]] = []  # the domains of each item holding several

    def root(d: int) -> int:
        while bound[d] != d:
            d = bound[d]
        return d

    for n, it in enumerate(items):
        votes = [owner[m.sse.id] for m in it.members if isinstance(m, Placed) and m.sse.id in owner]
        if votes:
            groups.setdefault(max(set(votes), key=votes.count), []).append(n)
            if len(set(votes)) > 1:
                shared.append(sorted(set(votes)))
    keys = sorted(groups)
    if len(keys) < 2:
        return False
    unit, members = _domain_units(len(items), groups, ties)
    copy_k = np.zeros(len(items))
    for i, (_, k) in ties.items():
        copy_k[i] = k
    sx = np.where((copy_k > 0) & (mirror is not None), -1.0, 1.0)  # a mirrored copy moves against its mate in x

    rel = [it.box for it in items]  # each item's footprint about its centre: only centres move here (C33: the
    # property, measured afresh four times an item for every pair of domains, was 165 s of 7JU4's domains view)

    def box(d):
        b = np.array(
            [
                (
                    items[n].cx + rel[n][0],
                    items[n].cy + rel[n][1],
                    items[n].cx + rel[n][2],
                    items[n].cy + rel[n][3],
                )
                for n in groups[d]
            ]
        )
        return b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max()

    def shift(ns, dx, dy) -> None:
        for n in ns:
            it = items[n]
            step = _snap(np.array([sx[n] * dx, dy]))
            it.cx, it.cy = it.cx + float(step[0]), it.cy + float(step[1])
            it.home = (it.home[0] + float(step[0]), it.home[1] + float(step[1]))  # rest here from now on

    def overlap() -> float:  # summed over the pairs that may part
        total = 0.0
        for i, a in enumerate(keys):
            for b in keys[i + 1 :]:
                if root(a) != root(b):
                    A, B = box(a), box(b)
                    ox = min(A[2], B[2]) - max(A[0], B[0]) + _DOMAIN_GAP
                    oy = min(A[3], B[3]) - max(A[1], B[1]) + _DOMAIN_GAP
                    total += min(ox, oy) if ox > 0 and oy > 0 else 0.0
        return total

    start = [(it.cx, it.cy, it.home) for it in items]

    def extent() -> float:  # the items' width and height, summed
        b = np.array([(it.cx + r[0], it.cy + r[1], it.cx + r[2], it.cy + r[3]) for it, r in zip(items, rel)])
        return float(b[:, 2].max() - b[:, 0].min() + b[:, 3].max() - b[:, 1].min())

    before = extent()
    moved = parted = False
    for second in (False, True):
        if second:  # the sweeps never parted them and blew the page up: again, domains sharing an item joined
            if parted or not shared or extent() <= _PART_BLOWUP * before:
                break
            for it, (cx, cy, home) in zip(items, start):
                it.cx, it.cy, it.home = cx, cy, home
            for ds in shared:
                for d in ds:
                    a, b = root(ds[0]), root(d)
                    bound[max(a, b)] = min(a, b)
        moved, parted = _part_sweeps(
            keys, groups, unit, members, copy_k, mirror, box, shift, overlap, root, sweeps, second
        )
    return moved


def _part_sweeps(keys, groups, unit, members, copy_k, mirror, box, shift, overlap, root, sweeps: int, guard: bool):
    """_part_domains' sweeps: (anything moved, parted: a sweep found no clash); with `guard`, they also stop when a
    sweep no longer shrinks the overlap."""
    moved, last = False, (overlap() if guard else 0.0)  # summed only when it may stop the sweeps (its boxes cost)
    for _ in range(sweeps):
        clash = False
        for i, a in enumerate(keys):
            for b in keys[i + 1 :]:
                if root(a) == root(b):
                    continue
                A, B = box(a), box(b)
                ox = min(A[2], B[2]) - max(A[0], B[0]) + _DOMAIN_GAP
                oy = min(A[3], B[3]) - max(A[1], B[1]) + _DOMAIN_GAP
                if ox <= 0 or oy <= 0:
                    continue
                ka, kb = _domain_level(groups[a], copy_k), _domain_level(groups[b], copy_k)
                sa, sb = (-1.0 if k > 0 and mirror is not None else 1.0 for k in (ka, kb))
                s = 1.0 if (B[0] + B[2]) >= (A[0] + A[2]) else -1.0  # b lies this way of a, in x
                if unit[a] != unit[b]:
                    if ox <= oy:  # the cheaper direction, each side going half way (as the domain's copy moves)
                        shift(members[unit[a]], -s * ox / 2 * sa, 0.0)
                        shift(members[unit[b]], s * ox / 2 * sb, 0.0)
                    else:
                        t = 1.0 if (B[1] + B[3]) >= (A[1] + A[3]) else -1.0
                        shift(members[unit[a]], 0.0, -t * oy / 2)
                        shift(members[unit[b]], 0.0, t * oy / 2)
                elif sa != sb:  # a domain and its mirror image: apart sideways, symmetrically
                    shift(members[unit[a]], s * ox / (sb - sa), 0.0)
                elif kb != ka and mirror is None:  # copies of one domain: the spread widens
                    wider = s * ox / (kb - ka)
                    for n, k in enumerate(copy_k.tolist()):
                        if k:
                            shift([n], k * wider, 0.0)
                else:
                    continue  # no move keeps them alike: left to _separate
                clash = moved = True
        if not clash:
            return moved, True
        if guard:
            now = overlap()
            if now >= last - 1e-6:  # the moves no longer help
                break
            last = now
    return moved, False


def _domain_level(ns: list[int], copy_k: np.ndarray) -> float:
    """The protomer a domain's items (`ns`) belong to: the copy level most of them have, the lower on a tie."""
    levels = Counter(copy_k[n] for n in ns)
    return float(max(sorted(levels), key=lambda k: levels[k]))


def _domain_units(
    n: int, groups: dict[int, list[int]], ties: dict[int, tuple[int, int]]
) -> tuple[dict[int, int], dict[int, list[int]]]:
    """Domains that must move together to keep copies alike, and the items each such unit moves: the domain groups
    (domain -> its items) joined through the symmetry orbits of their items (ties: copy -> its protomer-0 mate),
    with every item of those orbits. (domain -> unit, unit -> items, in order.)"""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def join(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[max(a, b)] = min(a, b)

    for i, (j, _) in sorted(ties.items()):
        join(i, j)
    for d in sorted(groups):
        for m in groups[d][1:]:
            join(groups[d][0], m)
    unit = {d: find(groups[d][0]) for d in sorted(groups)}
    used = set(unit.values())
    members: dict[int, list[int]] = {}
    for m in range(n):
        if find(m) in used:
            members.setdefault(find(m), []).append(m)
    return unit, members


def block_ghosts(layout: Layout, k: int) -> list[Placed]:
    """The ghosts drawn in sheet block k (Layout.sheets): Layout.sheet_ghosts when recorded, else those repeating
    one of its strands (a barrel's)."""
    ghosts = {g.sse.id: g for g in layout.ghosts}
    mine = layout.sheet_ghosts[k] if layout.sheet_ghosts else layout.sheets[k]
    return [ghosts[i] for i in mine if i in ghosts]


def panel_members(layout: Layout, name: str) -> list[str]:
    """The elements a domain's panel is drawn round: its own, then those its neighbours lend to it (focus mode:
    drawn grey and primed, a neighbour's strand completing its sheet; Layout.domain_lent), which are no member."""
    return [*dict(layout.domains)[name], *layout.domain_lent.get(name, [])]


def domain_panel(layout: Layout, name: str) -> tuple[float, float, float, float]:
    """A domain's panel (its elements plus margin) with room above for its name. A ghost counts where it repeats a
    strand of its own block (a barrel's): one lent across an unrolled ring's cut stands in another copy's sheet."""
    ids = panel_members(layout, name)
    mine = set(ids)
    repeats = [
        g.rect
        for k, block in enumerate(layout.sheets)
        for g in block_ghosts(layout, k)
        if g.sse.id in block and g.sse.id in mine
    ]
    r = np.array([layout.placed[k].rect for k in ids] + repeats)
    x0, y0, x1, y1 = r[:, 0].min(), r[:, 1].min(), r[:, 2].max(), r[:, 3].max()
    return (float(x0 - DOMAIN_PAD), float(y0 - DOMAIN_PAD), float(x1 + DOMAIN_PAD), float(y1 + DOMAIN_PAD + 0.6))


def domain_name_box(layout: Layout, name: str) -> tuple[float, float, float, float]:
    """The box a domain's name is drawn in (_place_domain_names), on its plate."""
    return layout.domain_names[name]


_PLATE_PAD = 0.3  # a sheet's plate reaches this far beyond its strands (render._draw_sheet_panels)
DOMAIN_NAME_SIZE = 1.05  # a domain's name, in ems of the figure's text (render draws it bold at this size)
_NAME_HEIGHT = 1.0  # ems of its own size a name's box takes up and down (capitals, a descender, the plate's room)
_NAME_INSET = 0.35  # page units a name keeps from its panel's side
_NAME_SINK = 0.1  # page units a name's middle sits inside its panel's top edge: it lies on the edge, on its plate
_NAME_STEP = 0.25  # page units between the places tried along the top edge
_NAME_CLEAR = 0.08  # page units a name keeps from anything else
_NAME_KEEP = 0.3  # page units loops and nucleic-acid strands keep from a name: room kept round it where there is
# any, and its routing keep-out (label_boxes) grown by _NAME_ROUTE
_NAME_ROUTE = 0.15  # (a wider keep-out forces loops over it where room is short)
_NAME_NEAR_NT = 0.3  # page units a name keeps from a nucleotide of an unpaired strand, beyond its letter
_NAME_STUB = 0.8  # page units: the straight run of a loop out of an element's end, kept clear of a name
_SHORT_LOOP = 6.0  # page units: a loop whose ends are this near takes the short way between them, kept clear too
_LEADER = 0.8  # page units a name moved outside its panel stands above it, joined to it by a short leader


def _domain_boxes(layout: Layout) -> tuple[dict[str, list[tuple]], list[tuple]]:
    """Each drawn domain's element boxes (its elements, and the ghosts repeating a strand of their own block: a
    barrel's, as domain_panel counts them), and the boxes of what is in no domain (elements, ghosts, stems)."""
    owned: dict[str, list[tuple]] = {}
    used: set[str] = set()
    repeats: set[int] = set()
    for name, ids in layout.domains:
        if not ids:
            continue
        ids = panel_members(layout, name)
        mine = set(ids)
        own = [
            g
            for k, block in enumerate(layout.sheets)
            for g in block_ghosts(layout, k)
            if g.sse.id in block and g.sse.id in mine
        ]
        owned[name] = [layout.placed[k].rect for k in ids] + [g.rect for g in own]
        used |= mine
        repeats |= {id(g) for g in own}
    free = [p.rect for k, p in layout.placed.items() if k not in used]
    free += [g.rect for g in layout.ghosts if id(g) not in repeats]
    return owned, free


def domain_fields(layout: Layout, names: bool = True) -> panels.Fields:
    """The panels of the drawn domains as fields (panels.panel_fields): each follows its own elements, padded by
    DOMAIN_PAD, and wraps its name's box when the name sits on its edge (names=True; a name with a leader stands
    outside, on its own plate)."""
    owned, free = _domain_boxes(layout)
    stems = [(d.cx, d.cy, d.length / 2, d.width / 2, d.angle) for d in layout.stems]
    on_edge = {n: b for n, b in layout.domain_names.items() if n not in layout.domain_leaders} if names else None
    return panels.panel_fields(owned, free, DOMAIN_PAD, on_edge, stems)


def _name_obstacles(layout: Layout, sses: list[SSE]) -> tuple[list[tuple], list[tuple]]:
    """What a domain's name must keep clear of: (solid) elements and ghosts, stems with their letters, sheet plates
    and the other names (labels, chain ends and their stubs, stem names, tie letters, domain names placed before),
    and (ways) the straight run of a loop out of each element's end and a short loop's way between its two ends. A
    name's leader keeps off the solid ones only: it may cross a loop."""
    solid = [_grown(p.rect, _NAME_CLEAR) for p in [*layout.placed.values(), *layout.ghosts]]
    for k, ids in enumerate(layout.sheets):
        if len(ids) > 1:
            r = np.array([layout.placed[i].rect for i in ids] + [g.rect for g in block_ghosts(layout, k)])
            solid.append((*(r[:, :2].min(axis=0) - _PLATE_PAD), *(r[:, 2:].max(axis=0) + _PLATE_PAD)))
    solid += [_grown(b, _NAME_CLEAR) for _, b in label_boxes(layout, sses)]
    ways = []
    for q in [*layout.placed.values(), *layout.ghosts]:
        for port, out in ((q.n_port, q.exit_n), (q.c_port, q.exit_c)):
            (ax, ay), (bx, by) = port, np.add(port, out * _NAME_STUB)
            ways.append((min(ax, bx) - 0.2, min(ay, by) - 0.2, max(ax, bx) + 0.2, max(ay, by) + 0.2))
    for mine in group_by_chain(sses, keep=layout.placed).values():
        for a, b in zip(mine, mine[1:]):
            p, q = layout.placed[a.id], layout.placed[b.id]
            if np.hypot(*np.subtract(p.c_port, q.n_port)) <= _SHORT_LOOP:
                ends = np.array(
                    [
                        p.c_port,
                        np.add(p.c_port, p.exit_c * _NAME_STUB),
                        q.n_port,
                        np.add(q.n_port, q.exit_n * _NAME_STUB),
                    ]
                )
                ways.append((*(ends.min(axis=0) - 0.2), *(ends.max(axis=0) + 0.2)))
    return solid, ways


def _place_domain_names(layout: Layout, sses: list[SSE], font_scale: float = 1.0) -> None:
    """Each drawn domain's name (Layout.domain_names), in domain order, at the first clear place (_name_obstacles,
    and the unpaired strands' traces) of _name_spots: its panel's top left and top right on the top edge, then
    along that edge, then outside with a short leader (Layout.domain_leaders). When nothing is clear, the top
    left."""
    layout.domain_names, layout.domain_leaders = {}, {}
    fields = domain_fields(layout, names=False)
    if not fields.values:
        return
    em = _ink_em(layout, font_scale) or text_em(layout.width or 1.0) * font_scale
    traces = [np.asarray(pts, float).reshape(-1, 2) for _, pts in sorted(layout.na_traces.items())]
    solid, ways = _name_obstacles(layout, sses)

    letters = (0.0 if layout.na_cropped else NT_ROW + glyph(layout)[1]) + _NAME_CLEAR

    def clear(box, foot, margin) -> bool:
        if overlaps_any(_grown(box, margin), solid + ways) or _trace_hits(_grown(box, _NAME_NEAR_NT), traces):
            return False
        if any(_box_meets_body(_grown(box, letters), d) for d in layout.stems):  # a stem, with its letters
            return False
        if foot is None:
            return True
        (x, ya), (_, yb) = foot
        return not overlaps_any((x - 0.05, min(ya, yb), x + 0.05, max(ya, yb)), solid)

    for name, value in fields.values.items():
        size = DOMAIN_NAME_SIZE * em
        spots = list(_name_spots(value > 0, fields.xs, fields.ys, size * _text_width(name) + 0.1, size * _NAME_HEIGHT))
        if not spots:
            continue
        box, foot = next(  # a place with room for the margin loops keep (_NAME_KEEP), else one without
            (spot for margin in (_NAME_KEEP, 0.0) for spot in spots if clear(*spot, margin)), spots[0]
        )
        layout.domain_names[name] = box
        if foot is not None:
            layout.domain_leaders[name] = foot
        solid.append(_grown(box, _NAME_CLEAR))


def _name_spots(inside: np.ndarray, xs: np.ndarray, ys: np.ndarray, w: float, h: float):
    """The places a name of w x h may take on a panel (its cells `inside`, rows along ys), in order, as (box,
    leader or None): on the top edge (the name's middle _NAME_SINK inside it, its plate over the edge) at the top
    left, the top right, then along the edge left to right, _NAME_INSET from the panel's sides; then outside, _LEADER
    and twice that above the edge, with a leader down to it; then the same on the bottom edge. The largest piece of
    a panel in pieces first, at each of these steps; on a piece narrower than the name, the name centred over it."""
    from scipy.ndimage import label as pieces

    marks, n = pieces(inside)
    order = sorted(range(1, n + 1), key=lambda m: (-int((marks == m).sum()), m))
    rows = np.arange(len(ys))[:, None]
    half = panels.GRID / 2
    for up in (1, -1):
        for lift in (0.0, _LEADER, 2 * _LEADER):
            for m in order:
                mine = marks == m
                cols = np.flatnonzero(mine.any(axis=0))
                edge = np.full(len(xs), np.nan)
                if up > 0:
                    edge[cols] = ys[np.where(mine, rows, -1).max(axis=0)[cols]] + half
                else:
                    edge[cols] = ys[np.where(mine, rows, len(ys)).min(axis=0)[cols]] - half
                x0, x1 = xs[cols[0]] - half + _NAME_INSET, xs[cols[-1]] + half - _NAME_INSET
                if x1 - x0 < w:  # a narrow piece (a helix on end): the name centred over it, its plate overhanging
                    starts = [(x0 + x1 - w) / 2]
                else:
                    starts = [*np.arange(x0, x1 - w, _NAME_STEP), x1 - w]
                for bx in dict.fromkeys([starts[0], starts[-1], *starts[1:-1]]):
                    span = (xs >= bx) & (xs <= bx + w) & ~np.isnan(edge)
                    if not span.any():
                        continue
                    at = float(edge[span].max() if up > 0 else edge[span].min())
                    mid = at + up * (lift + h / 2 if lift else -_NAME_SINK)
                    box = tuple(float(v) for v in _snap([bx, mid - h / 2, bx + w, mid + h / 2]))
                    foot = None
                    if lift:
                        x = float(_snap(bx + w / 2))
                        below = float(_snap(edge[span][np.argmin(np.abs(xs[span] - x))]))
                        foot = ((x, box[1] if up > 0 else box[3]), (x, below))
                    yield box, foot


def provisional(sses: list[SSE]) -> Layout:
    """A layout holding only element labels, for resolving references before anything is placed."""
    placed = [Placed(s, 0.0, 0.0, 1.0, 0.0, 1.0) for s in sses]
    _assign_labels(placed)
    return Layout({p.sse.id: p for p in placed}, (0.0, 0.0, 0.0, 0.0))


def helix_label_pos(p: Placed) -> tuple[float, float]:
    """Beside the bar, on its upper side (the other with label_side -1), clear of the body; label_along further
    along it."""
    nx, ny = -np.sin(p.angle), np.cos(p.angle)
    side = (p.width / 2 + LABEL_GAP) * (1 if ny >= 0 else -1) * p.label_side
    return float(p.cx + nx * side + ny * p.label_along), float(p.cy + ny * side - nx * p.label_along)


def termini(layout: Layout, sses: list[SSE]) -> list[tuple[str, str, tuple[float, float], np.ndarray]]:
    """(N or C, chain, port, outward direction) for each chain's first and last placed element."""
    out = []
    for chain, mine in group_by_chain(sses, keep=layout.placed).items():
        if chain not in layout.partial:
            first, last = layout.placed[mine[0].id], layout.placed[mine[-1].id]
            out += [("N", chain, first.n_port, first.exit_n), ("C", chain, last.c_port, last.exit_c)]
    return out


def centred_box(cx: float, cy: float, hx: float, hy: float) -> tuple[float, float, float, float]:
    """(x0, y0, x1, y1) of the box with centre (cx, cy) and half-sizes (hx, hy)."""
    return (cx - hx, cy - hy, cx + hx, cy + hy)


def overlaps_any(box, rects) -> bool:
    """True when the box overlaps any of the rects (x0, y0, x1, y1), touching edges excluded."""
    return any(box[0] < t[2] and t[0] < box[2] and box[1] < t[3] and t[1] < box[3] for t in rects)


def dna_y(slot: int, x: np.ndarray | float, width: float = DNA_W) -> np.ndarray | float:
    """Sideways offset of backbone `slot` (0 or 1) at axial position x on a drawn duplex (page units)."""
    return width / 2 * np.sin(2 * np.pi * (np.asarray(x) / DNA_PITCH - slot * DNA_MINOR))


def dna_to_page(d: PlacedStem, x, y) -> np.ndarray:
    """Duplex-local (axial, sideways) to page coordinates."""
    c, s = np.cos(d.angle), np.sin(d.angle)
    x, y = np.asarray(x, float), np.asarray(y, float)
    return np.stack([d.cx + x * c - y * s, d.cy + x * s + y * c], axis=-1)


NT_ROW = 0.55  # a base letter sits this far out from its nucleotide


def stem_label_pos(d: PlacedStem) -> tuple[float, float]:
    """Where a stem's name goes: beside it on its first rail's side (the second's with label_side -1). Once placed
    (label_at, _choose_label_sides), midway along its letter row on that side and label_at from its axis,
    _INK_MARGIN em beyond its letters; before that (the room placement keeps for it) a glyph beyond them, level with
    its centre."""
    if d.label_at is None:
        room = NT_ROW + _GLYPH_H * _GLYPH_K.get() + d.label_out  # a glyph beyond its letters (glyph_scale)
        p = np.array([d.cx, d.cy]) + d.normal * d.label_side * (d.width / 2 + room)
    else:
        p = np.array([d.cx, d.cy]) + d.direction * d.label_along + d.normal * d.label_side * d.label_at
    return float(p[0]), float(p[1])


def _middle_pair(d: PlacedStem, side: int, letters: bool) -> float:
    """Where along the stem, from its centre, a name centred on it sits: midway along its letter row on `side` (the
    first and last letters of that side's nucleotides), or, with no letters drawn, midway along its middle pair
    (pairs)."""
    if not d.pairs:
        return 0.0
    if letters:
        at = [d.axial[key] for key in d.axial if d.side[key] == side]
        return float(_snap((min(at) + max(at)) / 2)) if at else 0.0
    at = sorted((d.axial[(a, i)] + d.axial[(b, j)]) / 2 for a, i, b, j in d.pairs)
    return float(_snap((at[(len(at) - 1) // 2] + at[len(at) // 2]) / 2))


def stem_keep_box(d: PlacedStem, k: float | None = None) -> tuple[float, float, float, float]:
    """The box routes keep out of round a stem's name: its ink grown by _INK_MARGIN em (label_keep) where its place
    is chosen, else its text box (_stem_label_box, glyphs `k`)."""
    if d.label_keep is None:
        return _stem_label_box(d, k)
    x, y = stem_label_pos(d)
    x0, y0, x1, y1 = d.label_keep
    return x + x0, y + y0, x + x1, y + y1


def _stem_label_box(d: PlacedStem, k: float | None = None) -> tuple[float, float, float, float]:
    """The box around a stem's name, written along the stem: glyphs `k` (glyph_scale; the layout being built's,
    _GLYPH_K, when None) times _GLYPH_W x _GLYPH_H."""
    k = _GLYPH_K.get() if k is None else k
    along, across = _GLYPH_W * k * len(d.name) / 2 + 0.05, _GLYPH_H * k / 2
    c, s = abs(np.cos(d.angle)), abs(np.sin(d.angle))
    return centred_box(*stem_label_pos(d), c * along + s * across, s * along + c * across)


def _box_meets_body(box, d) -> bool:
    """Does the box (x0, y0, x1, y1) overlap the body (length x width, turned by its angle) of the stem or
    element d? Separating axes: the page's two and the body's two."""
    x0, y0, x1, y1 = box
    c = np.array([(x0 + x1) / 2 - d.cx, (y0 + y1) / 2 - d.cy])
    hx, hy = (x1 - x0) / 2, (y1 - y0) / 2
    for u in (np.array([1.0, 0.0]), np.array([0.0, 1.0]), d.direction, d.normal):
        reach = abs(u[0]) * hx + abs(u[1]) * hy + abs(u @ d.direction) * d.length / 2 + abs(u @ d.normal) * d.width / 2
        if abs(u @ c) >= reach:
            return False
    return True


def _body_distance(q, d) -> float:
    """How far the point q is from the body (length x width, turned) of the stem d; 0 inside it."""
    v = np.subtract(q, (d.cx, d.cy))
    a, b = abs(v @ d.direction) - d.length / 2, abs(v @ d.normal) - d.width / 2
    return float(np.hypot(max(a, 0.0), max(b, 0.0)))


_SLOTS = (0.25, 0.5, 0.75)  # page units a stem's name may move further out than beside its letters, when neither
# side is clear there
_INK_MARGIN = 1.0  # ems kept clear between a stem's name and every base letter (its own letter row, or its rail
# when no letters are drawn), and by every route: at 0.25 em 1EHZ's D and T read as letters of their rows, at 0.5 em
# T still stacked on a G of the T arm (plan 10 Task 4 review, fix round 1)
_NAME_EM, _LETTER_EM = 0.8, 0.75  # ems: render draws a stem's name, and a base letter, at these sizes
DNA_TUBE = 0.42  # a drawn double helix's backbone tube width, page units (render's ribbon)
_ink_lines: dict[tuple[str, bool], tuple[np.ndarray, np.ndarray]] = {}


def _ink_line(text: str, bold: bool) -> tuple[np.ndarray, np.ndarray]:
    """One line of text as render draws stem names and base letters (text_anchor), in ems, unturned: its ink (its
    glyphs' convex hulls from text_metrics at their pen positions, with pair kerning; inside the drawn ink) around
    its centre, and where its pen starts on its baseline from that centre. The centre is the middle of its line box,
    which spans its advance widths along and the font's typographic ascender and descender (or the ink, where it
    reaches past them) up and down: where matplotlib 3.11 centres a line (ha and va "center"). Blank text: no ink."""
    key = (text, bold)
    if key not in _ink_lines:
        tm = text_metrics
        advance, kern, hull = (
            (tm.ADVANCE, tm.KERN, tm.HULL) if bold else (tm.ADVANCE_NORMAL, tm.KERN_NORMAL, tm.HULL_NORMAL)
        )
        pen, parts = 0.0, []
        for k, c in enumerate(text):
            pen += kern.get(text[k - 1 : k + 1], 0) if k else 0
            if c in hull:
                parts.append(np.array(hull[c].split(), float).reshape(-1, 2) + (pen, 0.0))
            pen += advance.get(c, _FALLBACK_WIDTH * tm.UNITS)
        p = np.vstack(parts) if parts else np.empty((0, 2))
        up, down = max(p[:, 1].max(initial=0.0), tm.ASCENT), max(-p[:, 1].min(initial=0.0), tm.DESCENT)
        start = np.array([-pen / 2, -(up - down) / 2])
        _ink_lines[key] = ((p + start) / tm.UNITS, start / tm.UNITS)
    return _ink_lines[key]


def _turn(points: np.ndarray, rotation: float) -> np.ndarray:
    """Points turned by `rotation` degrees, anticlockwise."""
    t = np.radians(rotation)
    c, s = np.cos(t), np.sin(t)
    return points @ np.array([[c, s], [-s, c]])


def _ink_points(text: str, bold: bool, size: float, angle: float) -> np.ndarray:
    """The text's ink (_ink_line) drawn `size` page units to the em, turned as render turns stem text (upright along
    `angle`), around the point it is centred on."""
    return _turn(_ink_line(text, bold)[0] * size, upright_deg(angle))


def text_anchor(text: str, bold: bool, centre, size: float, rotation: float) -> tuple[float, float]:
    """Where render starts a stem name or a base letter (ha "left", va "baseline", rotation_mode "anchor") so that
    it is centred on `centre` as the layout measures it (_ink_line), `size` page units to the em, turned by
    `rotation` degrees: matplotlib 3.10 and 3.11 place a line's glyphs from its start alike, but centre it on boxes of
    their own (3.10's grow with hinting at small sizes and spaces: 1KX5's "DNA 1d" sat 0.05 page units along)."""
    x, y = np.asarray(centre, float) + _turn(_ink_line(text, bold)[1] * size, rotation)
    return float(x), float(y)


def _reach(points: np.ndarray) -> np.ndarray:
    """(x0, y0, x1, y1) of the points, around the point they are drawn about (zero for no ink)."""
    if not len(points):
        return np.zeros(4)
    return np.concatenate([points.min(axis=0), points.max(axis=0)])


def _box_gaps(box, boxes: np.ndarray) -> np.ndarray:
    """The distance from the box (x0, y0, x1, y1) to each of the boxes (n, 4); 0 where they meet (as figcheck
    measures a stem name's gap)."""
    if not len(boxes):
        return np.empty(0)
    dx = np.maximum(np.maximum(box[0] - boxes[:, 2], boxes[:, 0] - box[2]), 0.0)
    dy = np.maximum(np.maximum(box[1] - boxes[:, 3], boxes[:, 1] - box[3]), 0.0)
    return np.hypot(dx, dy)


def upright_deg(angle: float) -> float:
    """The rotation (degrees) render writes stem text with along a stem at `angle` (radians): turned to read left
    to right."""
    rot = float(np.degrees(angle))
    return rot - 180 if rot > 90 else rot + 180 if rot < -90 else rot


def first_pair_text(d: PlacedStem, chains: list, seqs: list) -> str:
    """The block style's text beside a stem: its first (outermost) pair's two residue numbers (seqs: per strand, each
    nucleotide's number as the file writes it), with the partner's chain (chains: each strand's) across chains."""
    a, i, b, j = d.pairs[0]
    partner = f"{chains[b]}:" if chains[b] != chains[a] else ""
    return f"{seqs[a][i]}\u00b7{partner}{seqs[b][j]}"


_FALLBACK_WIDTH = 0.6  # ems: a character the metric table lacks (about a bold lower-case letter)


def _text_width(text: str, bold: bool = True) -> float:
    """A bold name's width in ems of its own size, from the built-in metrics of Liberation Sans Bold (Arial Bold's
    metrics; text_metrics.py; Regular when not `bold`), never from the fonts this machine happens to have, so the
    layout is the same everywhere: its advances plus pair kerning, as render's TextPath measure adds them up;
    snapped."""
    advance, kern = text_metrics.ADVANCE, text_metrics.KERN
    if not bold:
        advance, kern = text_metrics.ADVANCE_NORMAL, text_metrics.KERN_NORMAL
    fallback = _FALLBACK_WIDTH * text_metrics.UNITS
    pen = sum(advance.get(c, fallback) for c in text) + sum(kern.get(a + b, 0) for a, b in zip(text, text[1:]))
    return float(_snap(pen / text_metrics.UNITS))


def page_text(layout: Layout) -> tuple[float, float]:
    """The figure's text size in points (Layout.font_scale applied) and page units per point, on the page the
    placed elements, stems and traces span so far (as _ink_em; a layout with nothing placed: its bounds)."""
    rects = [p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems] if hasattr(p, "rect")]
    rects += [(*pts.min(axis=0), *pts.max(axis=0)) for pts in map(np.asarray, layout.na_traces.values())]
    if rects:
        r = np.asarray(rects, float)
        width, height = float(r[:, 2].max() - r[:, 0].min()) + 2.0, float(r[:, 3].max() - r[:, 1].min()) + 2.0
    else:
        width, height = layout.width or 1.0, layout.height or 0.0
    scale = page_scale(width, height)
    return font_pt(scale) * layout.font_scale, 1.0 / (scale * 72)


def line_clearance(layout: Layout) -> float:
    """Page units a line keeps from text: figcheck.LINE_TEXT_CLEARANCE_PX and half a pixel of line width, at the
    page's text scale (page_text; 100 dpi, as figcheck measures), so as many pixels on a wide page of small text as
    on a small one (fix round 1: C14 and C11 kept a share of an em, short of it below 7 pt)."""
    from .figcheck import LINE_TEXT_CLEARANCE_PX

    return (LINE_TEXT_CLEARANCE_PX + 0.5) * 0.72 * page_text(layout)[1]


def _ink_em(layout: Layout, font_scale: float) -> float:
    """Page units of one em of the drawn text, from the page the placed elements, stems and traces span so far
    (the drawn page is at most 4% wider on the NA fixtures drawn with letters, whose text then grows with it)."""
    rects = [p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems]]
    rects += [(*pts.min(axis=0), *pts.max(axis=0)) for pts in map(np.asarray, layout.na_traces.values())]
    if not rects:
        return 0.0
    r = np.asarray(rects, float)
    return text_em(float(r[:, 2].max() - r[:, 0].min()) + 2.0) * font_scale


def _trace_hits(box, traces: list[np.ndarray]) -> int:
    """How many nucleotides of the traces (unpaired strands drawn along their projection) a name's box would sit on,
    counting the room their letters take beside them, plus the trace segments it would cross."""
    x0, y0, x1, y1 = box
    reach = NT_ROW + _GLYPH_H * _GLYPH_K.get() / 2  # a nucleotide's bead and its letter beside it (glyph_scale)
    n = 0
    for pts in traces:
        n += int(
            np.count_nonzero(
                (pts[:, 0] > x0 - reach)
                & (pts[:, 0] < x1 + reach)
                & (pts[:, 1] > y0 - reach)
                & (pts[:, 1] < y1 + reach)
            )
        )
        n += sum(1 for a, b in zip(pts, pts[1:]) if _segment_meets_box(a, b, box, _TRACE_HALF / 2))
    return n


def _segment_meets_box(a, b, box, pad: float) -> bool:
    """Does segment a-b pass within `pad` of the box (Liang-Barsky on the box grown by pad)?"""
    x0, y0, x1, y1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    t0, t1 = 0.0, 1.0
    dx, dy = b[0] - a[0], b[1] - a[1]
    for p, q in ((-dx, a[0] - x0), (dx, x1 - a[0]), (-dy, a[1] - y0), (dy, y1 - a[1])):
        if p == 0:
            if q < 0:
                return False
            continue
        t = q / p
        if p < 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return False
    return True


def _element_label_sides(layout: Layout, sses: list[SSE]) -> None:
    """A helix's name beside its upper side (label_side 1) unless it lies on a nucleic-acid trace, its beads or its
    letters (1EYG's ssDNA over chain A's α1): then on the other side, if that lies on fewer of them and meets no
    element, stem, other name or chain end."""
    traces = [np.asarray(pts, float).reshape(-1, 2) for _, pts in sorted(layout.na_traces.items())]
    if not traces:
        return
    k = glyph(layout)[0] / _GLYPH_W
    for p in layout.placed.values():
        if p.sse.kind not in ("H", "G"):
            continue
        p.label_side = 1
        first = _trace_hits(_element_label_box(p, k), traces)
        if not first:
            continue
        p.label_side = -1
        box = _element_label_box(p, k)
        taken = [q.rect for q in [*layout.placed.values(), *layout.ghosts] if q is not p]
        taken += [b for name, b in label_boxes(layout, sses) if name != f"label:{p.sse.id}"]
        if (
            _trace_hits(box, traces) >= first
            or overlaps_any(box, taken)
            or any(_box_meets_body(box, d) for d in layout.stems)
        ):
            p.label_side = 1


_OVER_STUB = 0.6  # page units: the straight run of a loop out of a strand's end, kept clear of a name over the sheet
_OVER_ALONG = (0.0, -1 / 3, 1 / 3)  # where a name over a sheet may sit along its helix: helix lengths from its middle


def _box_overlap(a, b) -> float:
    return max(0.0, min(a[2], b[2]) - max(a[0], b[0])) * max(0.0, min(a[3], b[3]) - max(a[1], b[1]))


def _over_label_sides(layout: Layout) -> None:
    """A helix drawn over a sheet lies on its strands: its name goes beside it where it meets the fewest strand ends
    (a loop leaves each strand straight out of its end) and arrowheads, then covers least of the strands' letters
    and of the other names, then lies nearest the middle (_OVER_ALONG: the middle, or a third of the way to either
    end), upper side first on a tie; so the name stays readable and clear of the loops."""
    strands = [q for q in [*layout.placed.values(), *layout.ghosts] if q.sse.kind == "E"]
    gw, gh = glyph(layout)
    kg = gw / _GLYPH_W
    letters = [centred_box(q.cx, q.cy, gw * len(q.label) / 2 + 0.05, gh / 2) for q in strands]
    ends = []
    for q in strands:
        for port, out in ((q.n_port, q.exit_n), (q.c_port, q.exit_c)):
            (ax, ay), (bx, by) = port, np.add(port, out * _OVER_STUB)
            ends.append((min(ax, bx) - 0.15, min(ay, by) - 0.15, max(ax, bx) + 0.15, max(ay, by) + 0.15))
        tip, side = np.array(q.c_port), np.array([-q.direction[1], q.direction[0]]) * HEAD / 2 * q.width / PITCH
        base = tip - q.direction * min(HEAD_LEN, q.length * 0.45)
        xs, ys = zip(tip, base + side, base - side)  # its arrowhead
        ends.append((min(xs), min(ys), max(xs), max(ys)))
    for p in layout.placed.values():
        if not p.over:
            continue
        names = [_element_label_box(q, kg) for q in layout.placed.values() if q is not p and q.sse.kind in ("H", "G")]
        best = None
        for k, along in enumerate(_OVER_ALONG):
            for side in (1, -1):
                p.label_side, p.label_along = side, float(_snap(along * p.length))
                box = _element_label_box(p, kg)
                cost = (
                    sum(_box_overlap(box, e) > 0 for e in ends),
                    float(_snap(sum(_box_overlap(box, b) for b in letters + names))),
                    k,
                )
                if best is None or cost < best[0]:
                    best = (cost, side, p.label_along)
        _, p.label_side, p.label_along = best


def _mirror_sides(mirrored: list[tuple]) -> None:
    """Each mirror copy's name on the mirror image of its mate's side: its offset from its body points the way the
    mate's does, reflected in the page's (vertical) dyad line."""
    for copy, mate in mirrored:
        if type(copy) is not type(mate) or isinstance(copy, Placed) and copy.sse.kind != mate.sse.kind:
            continue
        if isinstance(copy, Placed) and copy.sse.kind not in ("H", "G"):
            continue
        where = stem_label_pos if isinstance(copy, PlacedStem) else helix_label_pos
        mx, my = np.subtract(where(mate), (mate.cx, mate.cy))
        copy.label_side = 1
        if isinstance(copy, PlacedStem):
            copy.label_out = mate.label_out
        if float(np.dot(np.subtract(where(copy), (copy.cx, copy.cy)), (-mx, my))) < 0:
            copy.label_side = -1


def _choose_label_sides(layout: Layout, sses: list[SSE], style: Style, mirrored: list[tuple] = ()) -> None:
    """Element names first (_element_label_sides). Then each stem's name beside its own letter row (_name_at): its
    ink _INK_MARGIN em beyond the row's ink across the stem, centred along the row; beside its first
    rail (label_side 1) unless its box meets another stem's name decided before it, another stem's body, an element
    or an element's name, or sits nearer another stem than its own, or its ink comes within _INK_MARGIN em of a base
    letter of another stem (where the style draws them); then beside the second rail if that is clear; then further
    out (_SLOTS), either side; else (logged) the slot with the most room, nearer its own stem where any slot is
    (_pick_slot). A mirror copy's name (the (copy, mate) members in `mirrored`, of protomers drawn as mirror images)
    then takes the mirror of its mate's side (_mirror_sides). Each name's ink box grown by _INK_MARGIN em is what
    routes keep out of (label_keep; not in block style, which draws the name on the block). Runs once the names (and
    so the boxes' widths) are final."""
    _element_label_sides(layout, sses)
    _over_label_sides(layout)
    if mirrored:
        _mirror_sides([(a, b) for a, b in mirrored if isinstance(a, Placed)])
    em = _ink_em(layout, style.font_scale)
    rects = [p.rect for p in [*layout.placed.values(), *layout.ghosts]]
    rects += [b for name, b in label_boxes(Layout(layout.placed, layout.bounds), []) if name.startswith("label:")]
    na_style = "block" if style.na_style == "ladder" and layout.na_cropped else style.na_style  # as render draws it
    row = {"ladder": NT_ROW, "letters": 0.0}.get(na_style) if style.nucleotide_labels else None  # stems' letters
    layout.na_letter_rows = (row, bool(style.nucleotide_labels) and na_style != "block")
    bases = [s.bases for s in layout.nucleic.strands] if layout.nucleic is not None else []
    chains = [s.chain for s in layout.nucleic.strands] if layout.nucleic is not None else []
    seqs = [[lab.seq for lab in s.labels] for s in layout.nucleic.strands] if layout.nucleic is not None else []
    letter_ink = {  # stem -> its nucleotides with the ink of each one's base letter on the page, as render draws them
        id(d): [
            (key, _ink_points(bases[key[0]][key[1]], False, _LETTER_EM * em, d.angle) + d.point(key, row))
            for key in sorted(d.axial)
        ]
        if row is not None and bases
        else []
        for d in layout.stems
    }
    letters = {id(d): np.array([_reach(p) for _, p in letter_ink[id(d)]]).reshape(-1, 4) for d in layout.stems}
    block = na_style == "block"  # its slot beside the stem holds the outermost pair's numbers, the name inside
    ink = {  # stem -> the ink of the text drawn in its name's slot, around its centre
        id(d): _ink_points(
            first_pair_text(d, chains, seqs) if block else d.name,
            not block,
            (_LETTER_EM if block else _NAME_EM) * em,
            d.angle,
        )
        for d in layout.stems
    }
    beside = {}  # (stem, side) -> its name's (label_at, label_along) beside its letters (or its edge) there
    for d in layout.stems:
        edge = (d.width + DNA_TUBE * d.width / DNA_W) / 2 if na_style == "ribbon" else d.rail  # drawn, without letters
        for side in (1, -1):
            beside[id(d), side] = _name_at(d, side, ink[id(d)], letter_ink[id(d)], edge, _INK_MARGIN * em)
    others_letters = {
        id(d): np.vstack([letters[id(e)] for e in layout.stems if e is not d] or [np.empty((0, 4))])
        for d in layout.stems
    }
    named: list[tuple[float, float, float, float]] = []
    slots = [(side, out) for out in (0.0, *_SLOTS) for side in (1, -1)]
    for d in layout.stems:
        others = [e for e in layout.stems if e is not d]
        reach = _reach(ink[id(d)])
        rooms = []
        for s in slots:
            rooms.append(_name_room(d, s, beside, others, named, rects, reach, others_letters[id(d)]))
            if _pick_slot(rooms[-1:], _INK_MARGIN * em) == (0, False):
                break
        k, fallback = _pick_slot(rooms, _INK_MARGIN * em)
        if fallback:
            clear, nearer, room = rooms[k]
            log.info(
                "stem %s: no slot is clear, nearer its own stem and %.2f em from every other letter; took side %d "
                "out %.2f (clear %s, nearer %s, room %.3f)",
                d.name,
                _INK_MARGIN,
                *slots[k],
                clear,
                nearer,
                room,
            )
        d.label_side, d.label_out = slots[k]
        d.label_at, d.label_along = beside[id(d), d.label_side][0] + d.label_out, beside[id(d), d.label_side][1]
        named.append(_stem_label_box(d))
    if mirrored:
        _mirror_sides([(a, b) for a, b in mirrored if isinstance(a, PlacedStem)])
        for d in layout.stems:  # a copy's name beside its own letters on the side it took
            d.label_at, d.label_along = beside[id(d), d.label_side][0] + d.label_out, beside[id(d), d.label_side][1]
    for d in layout.stems:  # what routes keep out of: the name's ink and the margin round it; a block's name is
        # drawn on its body, its slot holds only the first pair's numbers (4YBB's 23 % more crossings with a margin)
        d.label_keep = None if block else _grown(tuple(float(v) for v in _reach(ink[id(d)])), _INK_MARGIN * em)


def _pick_slot(rooms: list[tuple[bool, bool, float]], margin: float) -> tuple[int, bool]:
    """Which slot a stem's name takes, from each slot's (clear: its box meets no name placed before it, element,
    element's name or other stem's body; nearer: it is nearer its own stem than any other; room: the gap from its
    ink to the nearest letter of another stem), in the order tried: the first clear and nearer one with `margin` of
    room; else (a fallback, logged) the one with the most room among the clear and nearer ones, else among the clear
    ones (a name must stay nearer its own stem than any other unless no slot allows it), else the first; the first of
    equals. -> (its index, whether it fell back)"""
    for k, (clear, nearer, room) in enumerate(rooms):
        if clear and nearer and room >= margin:
            return k, False
    return max(range(len(rooms)), key=lambda k: (rooms[k][0], rooms[k][1], rooms[k][2], -k)), True


def _name_at(
    d: PlacedStem, side: int, ink: np.ndarray, letters: list, edge: float, margin: float
) -> tuple[float, float]:
    """How far from the stem's axis, and along it from its centre (_middle_pair), its name (ink: its points around
    its centre, as drawn) is centred beside `side`: `margin` beyond its letter row's ink across the stem (letters:
    (nucleotide, its letter's ink on the page) of the stem; the row is that side's letters, measured in the stem's
    own frame, as figcheck.stem_name_gaps measures it), centred along the row; with no letters drawn, `margin`
    beyond the stem's drawn edge (`edge` from the axis). On the 1e-4 grid, outward."""
    normal = d.normal * side
    inner = float(np.max(ink @ -normal, initial=0.0))  # how far the name's ink reaches back toward the stem
    centre = np.array([d.cx, d.cy])
    row = [p for key, p in letters if d.side[key] == side]
    along = _middle_pair(d, side, bool(row))
    if not row:
        return math.ceil((edge + margin + inner) * _GRID) / _GRID, along
    out = max(float(np.max((p - centre) @ normal)) for p in row)  # the row's ink reaches this far out
    return math.ceil((out + margin + inner) * _GRID) / _GRID, along


def _name_room(d: PlacedStem, slot, beside, others, named, rects, ink, letters) -> tuple[bool, bool, float]:
    """The stem's name in this slot (label_side, label_out), as _pick_slot weighs it: (clear: its box meets no
    other stem's name decided before it, element, element's name or other stem's body; nearer: it sits nearer its
    own stem than any other; room: the gap from its ink box to the nearest letter of another stem, inf for none).
    See _choose_label_sides. Leaves d there."""
    d.label_side, d.label_out = slot
    d.label_at, d.label_along = beside[id(d), d.label_side][0] + d.label_out, beside[id(d), d.label_side][1]
    box, at = _stem_label_box(d), stem_label_pos(d)
    own = _body_distance(at, d)
    clear = not (overlaps_any(box, named) or overlaps_any(box, rects) or any(_box_meets_body(box, e) for e in others))
    nearer = all(_body_distance(at, e) > own for e in others)
    room = float(_box_gaps(ink + np.tile(at, 2), letters).min()) if len(letters) else math.inf
    return clear, nearer, room


def _grown(box, pad: float) -> tuple[float, float, float, float]:
    return box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad


def _element_label_box(p: Placed, k: float = 1.0) -> tuple[float, float, float, float]:
    """The box around a helix's name: `k` (glyph_scale of the text drawn) glyphs of _GLYPH_W x _GLYPH_H a letter."""
    return centred_box(*helix_label_pos(p), _GLYPH_W * k * len(p.label) / 2 + 0.05, _GLYPH_H * k / 2)


def na_text_boxes(layout: Layout) -> list[tuple[float, float, float, float]]:
    """The ink boxes of the nucleic acid's text as render draws it, for protein loops to keep off once it is placed
    (route.route_loops, random-100 C11): each single-stranded nucleotide's letter (Layout.na_letters), each stem's
    base letters on their row (Layout.na_letter_rows) and the strands' 5'/3' labels (a full em, bold)."""
    na = layout.nucleic
    row, single = layout.na_letter_rows
    if na is None:
        return []
    em = _ink_em(layout, layout.font_scale)
    out = []

    def add(text: str, bold: bool, size: float, angle: float, at) -> None:
        if text.strip():
            out.append(tuple(_reach(_ink_points(text, bold, size, angle) + np.asarray(at, float))))

    if single:
        for (s, i), xy in sorted(layout.na_letters.items()):
            add(na.strands[s].bases[i], False, _LETTER_EM * em, 0.0, xy)
    if row is not None:
        for d in layout.stems:
            for key in sorted(d.axial):
                add(na.strands[key[0]].bases[key[1]], False, _LETTER_EM * em, d.angle, d.point(key, row))
    for _, end, xy in layout.na_ends:
        add(end, True, em, 0.0, xy)
    return out


def thin_boxes(points, pad: float = 0.06) -> list[tuple[float, float, float, float]]:
    """A line as thin boxes along its segments (for box overlap tests)."""
    pts = np.asarray(points, float).reshape(-1, 2)
    return [(*(np.minimum(a, b) - pad), *(np.maximum(a, b) + pad)) for a, b in zip(pts, pts[1:])]


def count_half(text: str, size: float) -> tuple[float, float]:
    """Half the width and height of a count's or note's ink box: `text` at `size` page units per em."""
    return size * _text_width(text, bold=False) / 2 + 0.04, size * 0.4 + 0.04


def tail_counts(layout: Layout, sses: list[SSE], loops=()) -> list[dict]:
    """The residue count written beyond each long chain end (at least LONG_TAIL residues past its last element):
    end, chain, text ("14 aa"), centre and ink box (page units). Placed after routing, so not a route keep-out: the
    box is sized from the figure's text (TAIL_COUNT_SIZE em of _ink_em) and set from the N/C letter's far edge,
    along the stub, beside the letter on either side, then further along, then further out beside the letter and along
    the stub (_TAIL_COUNT_FAR), the first clear of elements, text (label_boxes), the loops' lines (`loops`,
    route.Loop) and the counts placed before; else the one of them meeting the fewest of those."""
    if not layout.res_chain:
        return []
    em = _ink_em(layout, layout.font_scale) or text_em(layout.width or 1.0) * layout.font_scale
    gw, gh = glyph(layout)
    groups = group_by_chain(sses)
    rows: dict[str, list[int]] = {}
    for r, ch in enumerate(layout.res_chain):
        rows.setdefault(ch, []).append(r)
    taken = [p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems]]
    taken += [b for _, b in label_boxes(layout, sses)]
    for loop in loops:  # loop lines as thin boxes along each segment
        taken += thin_boxes(loop.points)
    out = []
    for end, chain, port, ex in termini(layout, sses):
        mine, ours = groups.get(chain, []), rows.get(chain, [])
        if not mine or not ours:
            continue
        n = (mine[0].start - ours[0]) if end == "N" else (ours[-1] - mine[-1].end)
        if n < LONG_TAIL:
            continue
        text = f"{n} aa"
        size = TAIL_COUNT_SIZE * em
        half = count_half(text, size)
        ex = np.asarray(ex, float)
        letter = np.add(port, ex * END_LABEL)
        own = (max(gw / 2 + 0.05, em * _text_width(end) / 2), max(gh / 2, em * 0.4))  # its keep-out (label_boxes)
        spots = []
        for u in (ex, np.array([-ex[1], ex[0]]), np.array([ex[1], -ex[0]])):
            reach = abs(u[0]) * (own[0] + half[0]) + abs(u[1]) * (own[1] + half[1]) + _TAIL_COUNT_GAP
            spots.append(_snap(letter + u * reach))
        spots += [_snap(spots[0] + ex * k) for k in (0.3, 0.6)]
        # crowded (batch C2a: 7BBH, 6JO6, 5LDW): further out beside the letter, then further along its stub
        side = (np.array([-ex[1], ex[0]]), np.array([ex[1], -ex[0]]))
        spots += [_snap(spots[1 + j] + side[j] * k) for k in _TAIL_COUNT_FAR for j in (0, 1)]
        spots += [_snap(spots[0] + ex * k) for k in _TAIL_COUNT_FAR[2:]]
        boxes = [centred_box(float(c[0]), float(c[1]), *half) for c in spots]
        box = next((b for b in boxes if not overlaps_any(b, taken)), None)
        if box is None:  # nowhere clear: the spot meeting the fewest things, the first of them
            box = min(boxes, key=lambda b: sum(overlaps_any(b, (t,)) for t in taken))
        taken.append(box)
        out.append(
            {
                "end": end,
                "chain": chain,
                "text": text,
                "centre": ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2),
                "box": box,
            }
        )
    return out


def label_boxes(
    layout: Layout, sses: list[SSE], route: bool = False
) -> list[tuple[str, tuple[float, float, float, float]]]:
    """Keep-out boxes for text and chain-end stubs: loops route around them, the page bounds include them. With
    `route`, a stem's name is what routes keep out of: its ink grown by _INK_MARGIN em (stem_keep_box). Names and
    letters are sized from the text drawn (glyph)."""
    out = []
    gw, gh = glyph(layout)
    k = gw / _GLYPH_W
    for p in layout.placed.values():
        if p.sse.kind in ("H", "G"):
            out.append((f"label:{p.sse.id}", _element_label_box(p, k)))
    for end, chain, port, ex in termini(layout, sses):
        a, b = np.add(port, ex * 0.15), np.add(port, ex * END_STUB)
        out.append(
            (
                f"{end}-stub:{chain}",
                (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05),
            )
        )
        t = np.add(port, ex * END_LABEL)
        out.append((f"{end}:{chain}", centred_box(float(t[0]), float(t[1]), gw / 2 + 0.05, gh / 2)))
    for n, d in enumerate(layout.stems):
        out.append((f"stem-label:{n}", stem_keep_box(d, k) if route else _stem_label_box(d, k)))
    for n, tie in enumerate(layout.na_ties):  # a long tie's labels at its two marks
        for end, (_, _, (x, y)) in enumerate(tie.marks):
            out.append((f"tie-label:{n}:{end}", centred_box(x, y, gw / 2 + 0.05, gh / 2)))
    for name, box in layout.domain_names.items():  # a domain's name, on its panel's edge (_place_domain_names), with
        out.append((f"domain-label:{name}", _grown(box, _NAME_ROUTE)))  # room so no line or nucleotide runs by it
    for name, ((x, y0), (_, y1)) in layout.domain_leaders.items():
        out.append((f"domain-leader:{name}", (x - 0.1, min(y0, y1), x + 0.1, max(y0, y1))))
    return out


_EM_TUNED = 0.55 * SHAFT  # page units per em of the figure's text (before Style.font_scale), on every page
_MIN_SCALE = MIN_TEXT_PT / (_EM_TUNED * 72)  # inches per page unit at which the text is MIN_TEXT_PT


def page_scale(width: float, height: float = 0.0) -> float:
    """Inches per page unit of a page `width` x `height` page units, as render draws it: IN_PER_UNIT, less on a
    page wider than MAX_WIDTH_IN, but never so little that the text falls under MIN_TEXT_PT (the page grows past
    MAX_WIDTH_IN instead; random-100 C4: on wide pages text held at 4 pt outgrew every keep-out sized in page
    units), nor so much that a side passes MAX_PAGE_IN (past it the whole figure is scaled to fit, its text with
    it: the text keeps its page units, so nothing it keeps clear of changes; only its print size falls)."""
    scale = max(min(IN_PER_UNIT, MAX_WIDTH_IN / width), _MIN_SCALE)
    return min(scale, MAX_PAGE_IN / max(width, height))


def font_pt(scale: float) -> float:
    """The figure's text size in points (before Style.font_scale) at `scale` inches per page unit (page_scale):
    _EM_TUNED page units per em, at most 9 pt."""
    return float(min(_EM_TUNED * scale * 72, 9.0))


def text_em(width: float, height: float = 0.0) -> float:
    """How many page units one em of the figure's text takes on a page `width` x `height` units, as render sizes it
    (before Style.font_scale): _EM_TUNED on every page (page_scale grows a wide page rather than shrink its text
    under MIN_TEXT_PT, and scales a page past MAX_PAGE_IN as a whole)."""
    scale = page_scale(width, height)
    return font_pt(scale) / (scale * 72)


def glyph_scale(em: float) -> float:
    """How much the keep-outs sized for text in page units (the _GLYPH_W x _GLYPH_H glyph box of names, chain-end
    letters and nucleotide letters; render's legend row) grow for text drawn `em` page units per em: in
    proportion past _EM_TUNED, the em they were sized at, never below it (random-100 C4 ruling: every route
    keep-out and legend row is sized from the drawn text); snapped. The one helper they all go through."""
    return max(1.0, float(_snap(em / _EM_TUNED)))


def glyph(layout: Layout) -> tuple[float, float]:
    """A glyph's keep-out width and height (page units) for the text this layout is drawn with (glyph_scale of
    its em: text_em with Layout.font_scale)."""
    k = glyph_scale(text_em(layout.width or 1.0, layout.height or 0.0) * layout.font_scale)
    return _GLYPH_W * k, _GLYPH_H * k


def page_fit(width_in: float, height_in: float, font: float) -> str:
    """What the figure's notes say of its page (random-100 C4): '' while its text is MIN_TEXT_PT or larger (before
    Style.font_scale: `font`, the text size drawn, is compared at font_scale 1); past MAX_PAGE_IN, that the
    whole figure was scaled to fit and its text is smaller (rounded down: 5.97 pt is "5.9 pt", under 6)."""
    if font >= MIN_TEXT_PT - 1e-6:
        return ""
    return (
        f"the page ({width_in:.0f} x {height_in:.0f} in) is held at {MAX_PAGE_IN:.0f} in on its longer side, so its "
        f"text is {math.floor(font * 10) / 10:.1f} pt, under {MIN_TEXT_PT:.0f} pt; draw fewer chains (--chains, "
        "--focus) for larger text"
    )


NA_TORN_TURN = 90.0  # degrees: linked stems turned against one another by more than this when unrolled are torn


NA_HELD = 3  # residues of one protomer touching a stem piece for that protomer to hold it (as cli.NA_HOLD_MIN)


def _held_by_two(nucleic, sym, at: dict, n: int) -> bool:
    """Are most of the n stem pieces (more than half) each held by two protomers of `sym`, each touching it with at
    least NA_HELD residues? Unrolled, a piece is drawn beside one protomer only, so the other's contacts land a
    protomer's width away: when that is most of the nucleic acid, the unrolling tears it from the protein that holds
    it (3L2V: each viral DNA end of the PFV intasome is gripped by both inner integrase chains). One shared stem
    among several stays unrolled (30QM's P3, by the NA campaign's C4 ruling). `at`: nucleotide -> piece."""
    sector = {c: k for k, chains in enumerate(sym.protomers) for c in chains}
    holders: dict[int, dict[int, set[int]]] = {}  # piece -> protomer -> residues touching it
    for r, touched in sorted(nucleic.contacts.items()):
        k = sector.get(nucleic.contact_chain.get(r))
        if k is None:
            continue
        for key in sorted(touched):
            if key in at:
                holders.setdefault(at[key], {}).setdefault(k, set()).add(r)
    shared = sum(sum(len(rs) >= NA_HELD for rs in by.values()) >= 2 for by in holders.values())
    return 2 * shared > n


def na_torn_by_unrolling(nucleic, sses, sym) -> bool:
    """Would unrolling the assembly (symmetric_frame) tear the nucleic acid apart, against the plain projected frame
    (view_frame) the figure would otherwise use? Each stem is drawn whole where its centre is seen, so a single stem
    cannot be torn; what can: (a) a duplex split into pieces whose consecutive pieces' centres face different
    protomers (they land on either side of a seam: 1RVA's kinked site, drawn as a 180° hairpin), or (b) a group of
    stems linked through at most _SHORT_RUN nucleotides (a junction, a bent duplex, coaxial stacks) whose axes, after
    the best common turn, are turned against one another by more than NA_TORN_TURN between the two frames (4LD0's
    junction arms: 139°; 30QM's linked stems at most 57°), or (c) most stem pieces each held by two protomers
    (_held_by_two: 3L2V's two viral DNA ends, each gripped by both inner integrase chains)."""
    import math

    cyl = symmetric_frame(sses, sym)
    stems = [StemView(_snap(v.centroid), _snap(v.axis), v.n) for v in stem_views(nucleic.structure, nucleic.strands)]
    plain = view_frame(sses + stems, False, False, None, None)
    pieces = stem_pieces(nucleic)
    axes = [tuple(_snap(x) for x in stem_axis(p, nucleic.strands)) for p in pieces]

    def wrap(a: float) -> float:
        return (a + math.pi) % (2 * math.pi) - math.pi

    def facing(point) -> int:  # the protomer whose sector the point is unrolled into (as CylinderFrame places it)
        th = cyl._theta(point)
        return min(range(len(cyl.angles)), key=lambda k: (abs(wrap(th - cyl.angles[k])), k))

    for k in range(len(pieces) - 1):
        a, b = pieces[k], pieces[k + 1]
        if a.parent and a.parent == b.parent and facing(axes[k][0]) != facing(axes[k + 1][0]):
            return True
    at = {key: k for k, p in enumerate(pieces) for key in p.nucleotides()}
    if _held_by_two(nucleic, sym, at, len(pieces)):
        return True
    links = []  # pieces a strand steps between through at most _SHORT_RUN nucleotides
    for s, strand in enumerate(nucleic.strands):
        held = [i for i in range(len(strand)) if (s, i) in at]
        links += [(at[(s, i)], at[(s, j)]) for i, j in zip(held, held[1:]) if j - i - 1 <= _SHORT_RUN]
    group = _components(len(pieces), links)
    for g in sorted(set(group)):
        members = [k for k in range(len(pieces)) if group[k] == g]
        if len(members) < 2:
            continue
        turns = []
        for k in members:
            (cx, cy), (px, py) = cyl.direction(*axes[k]), plain.direction(*axes[k])
            turns.append(wrap(math.atan2(cy, cx) - math.atan2(py, px)))
        common = math.atan2(sum(math.sin(t) for t in turns), sum(math.cos(t) for t in turns))
        if max(abs(wrap(t - common)) for t in turns) > math.radians(NA_TORN_TURN):
            return True
    return False
