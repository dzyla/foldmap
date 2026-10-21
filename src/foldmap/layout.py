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
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from . import panels, text_metrics
from .features import KINK_GAP, helix_bundles
from .frame import Frame, RosetteFrame, symmetric_frame, view_frame
from .model import SSE, Nucleic, Sheet, generic_label, group_by_chain, pair_places, res_ref, same_place
from .nucleic import Stem, StemView, segment_stem, stem_axis, stem_views
from .style import CHOICES, Style

if TYPE_CHECKING:
    from .like import Like

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
LABEL_GAP = 0.3  # a helix label's centre sits this far beyond the side of the element box (times glyph_scale)
_NAME_PAD = 0.05  # a helix name's keep-out box reaches this far past its glyphs (times glyph_scale)
END_STUB = 0.7  # chain-end stub length
END_LABEL = 1.15  # distance from a chain-end port to its N/C letter
LONG_TAIL = 10  # a chain end with this many residues or more past its last element gets a count ("14 aa")
TAIL_COUNT_SIZE = 0.72  # that count's text size, as a fraction of the figure's
_TAIL_COUNT_GAP = 0.08  # page units between the N/C letter's ink and the count's
_TAIL_COUNT_NEAR = 1.5  # em: the most a tail count's ink stands from its letter's (it reads as the letter's)
_TAIL_COUNT_FAR = (0.3, 0.6, 0.9, 1.2)  # page units further out a count looks when its near spots are all taken
_GLYPH_W, _GLYPH_H = 0.32, 0.42  # rough label glyph size in page units, for keep-out boxes (at _EM_TUNED: glyph_scale)
_GLYPH_K: ContextVar[float] = ContextVar("foldmap_glyph_scale", default=1.0)  # glyph_scale of the text the layout
# being built is drawn with (build_layout sets it from Style.font_scale): for the glyph room of things that hold no
# layout (an item's footprint, a stem's name box, a trace's letters); layout.glyph(layout) wherever a layout is held
# The shared sheets whose way up is a symmetry tie that the layout being built draws with strands[0] down, and a dict
# it reports into (build_layout, _both_ways): the turnable sheets and whether the copies were drawn alike
_TIE_TURN: ContextVar[tuple[frozenset[str], dict | None]] = ContextVar("foldmap_tie_turn", default=(frozenset(), None))
IN_PER_UNIT = 0.30  # render: inches per page unit, at most
MAX_WIDTH_IN = 7.2  # render: a figure is drawn at most this wide (double column), inches, while its text stays ...
MIN_TEXT_PT = 6.0  # ... this size or larger: a wider page grows past MAX_WIDTH_IN instead (random-100 C4 ruling)
MAX_PAGE_IN = 48.0  # render: the longest side a page grows to (A0's, about); past it the whole figure, its text
# with it, is scaled to fit: text under MIN_TEXT_PT, said by the figure's notes (page_fit); never so far that text
# falls under matplotlib's 1 pt (page_scale: past that the page passes MAX_PAGE_IN)
_MAX_ROW = 30.0
_MIN_STRAND = 1.0  # shortest drawn strand (room for the arrow head)
# items: at or below it the dense energy keeps recorded layouts bit-identical (the pair-list summation order
# moves L-BFGS minima)
_PAIR_ENERGY_ABOVE = 200
_W_SEQ = 1.0  # pull between the two ends of a loop
_LOOP_SLACK = 1.5  # loop ends closer than this are left alone
_W_SS = 2.0  # pull between two bridged cysteines (twice a loop: disulfides should read as short bars)
_SS_SLACK = 0.6
# the second starts of one layout's stretched disulfides (_attract) spend at most this many pair evaluations of the
# energy (each evaluation counts its items squared, what it costs), counted as they run, never time; once spent, a
# disulfide still stretched keeps the first minimum. A second start solves the whole page again, up to five per
# bridge, so a page of many chains and bridges had no bound (random-100 held-out 8D3D's focus view: 896 items, 848
# bridges, about 2 s a start, killed after 900 s). Measured at plan 15 Task 7 on the default and focus views of both
# random-100 samples and every tests/data structure (419 layouts): the most spent was 3.1e9 (8RRX's focus view, 16 s),
# so none of them reaches it; 8D3D's focus view spends 3.0e9 a bridge. The budget is per _attract, and cli.make_layout
# can lay a figure out up to three times (the let_go second layout, the membrane view's relayout), each with its own
_SS_RESTART_PAIRS = 5_000_000_000
_W_CONTACT = 0.3  # pull between elements that touch in 3D (at _CONTACT_CAP CA pairs)
_CONTACT_CAP = 20
_K_HOME = 0.02  # pull back toward the projected position, keeps the fold recognisable
# an item drawn like a reference (build_layout's `like`, plan 13a): pulled to its mates' place in the reference's
# finished layout this much harder, so it stays there unless it must move to clear an overlap (tuned in 13a-10 only)
_K_LIKE = 8 * _K_HOME
LIKE_WEIGHT = _K_LIKE / _K_HOME  # Like.weight: the reference's home spring in units of a projected home's
_STEP = 0.15
_SWEEPS = 400
_SEP_SWEEPS = 5000  # push-apart passes before the guaranteed slide
_SEP_ORBIT_SWEEPS = 1000  # push-apart passes that keep symmetry copies alike, before they may move alone
_SEP_SYM_ROOM = 1.05  # copies drawn alike may take at most this much more area than items moved one by one
# ... and copies drawn as mirror images about a dyad in the page, pushed for as many passes as items moved one by
# one: the figure says they are mirror images (9GEP's nucleosome needed 1.07 and about 3,000 passes)
_SEP_MIRROR_ROOM = 1.10
# ... and, when neither holds, copies drawn alike with every item that has no copy yielding to them (_separate's
# second chance) may take up to this much (random-100 C2b-3: 8QEM's twelve ring copies drawn whole need 1.61 times
# the room of moving items one by one; 6JLN 1.31, 7QUP 1.10, 9P1X 1.06). Copies drawn alike are the rule.
_SEP_TIE_ROOM = 2.0
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
    label_side_default: int = 1  # the side its name was given before any later move (_element_label_sides, mirrors)

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
_HAIRPIN_STEP = 0.4545  # page units: the least step between neighbours on a hairpin loop, a bead's width (_TRACE_HALF,
# 0.45) and the 1 % the drawn polyline cuts off the arc


def _hairpin_step(half: float, r: float, m: int) -> float:
    """The straight step between neighbours of m nucleotides spaced over a hairpin loop's arc: radius r, its feet
    `half` either side of its chord's middle, the arc the long way round beyond the chord."""
    span = 2 * np.pi - 2 * np.arcsin(min(half / r, 1.0))
    return float(2 * r * np.sin(span / (2 * (m + 1))))


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
        # a bead apart (_TRACE_HALF): sized by the stem's width, the arc is only about half a circle and its
        # nucleotides stood closer (random-100 C2b-4: 7Y7C's 13- to 18-nucleotide loops, 0.40-0.44)
        # (with no nucleotide the step is the chord whatever the radius: nothing to space)
        while m > 0 and _hairpin_step(half, r, m) < _HAIRPIN_STEP:
            r = float(_snap(r + 0.01))
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
    domain_source: str = ""  # where the domains come from: "inferred" (--domains auto), "CATH", "UniProt", "" (given)
    domain_note: str = ""  # what the notes say of them (cli.domain_note): inferred, or CATH's, or none found
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
    # (N or C, chain) -> where that chain end's letter is drawn when a line ran through its usual spot (route_loops:
    # route._clear_termini); end_label_pos reads it
    end_labels: dict = field(default_factory=dict)
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
    # the gaps between two strands sharing a column (Sheet.accessory), as boxes (_stack_gaps): the loops leaving
    # their facing ends are kept apart there (route._apart_in_stack_gaps)
    stack_gaps: list[tuple] = field(default_factory=list)
    accessory: dict[str, str] = field(default_factory=dict)  # a strand drawn off its sheet's track (Sheet.accessory),
    # beside its parent -> that parent's id: no barrel strand (science S3)
    sheet_notes: list[str] = field(default_factory=list)  # what is said of the sheets as drawn (sheet_note_items)
    sheet_items: list[dict] = field(default_factory=list)  # the same as provenance items (sheet_note_items)
    # (accessory, partner) pairings drawn apart: no height level with the partner was free (_column_heights); the
    # notes name them ("also pairs with"), science S6 reads them
    unshown: list[tuple[str, str]] = field(default_factory=list)
    sheet_joins: list[tuple[str, str]] = field(default_factory=list)  # strand pairs joining two sheets built apart
    # (Sheet.joined: two barrels of one pairing network, 1EYG's OB folds)
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
    kink_links: set[tuple[str, str]] = field(default_factory=set)  # consecutive helix pieces of a copy whose mates
    # are a kink (features.KINK_GAP): drawn as their mates, their loop is a kink link too (_regroup_copies, route)
    let_go: int = 0  # copy items build_layout let go (_one_per_mate): make_layout keeps them only if more alike
    two_fold_assembly: bool = False  # a two-fold (C2) assembly off a straight duplex: Style.two_fold chooses how
    # make_layout draws it (the app shows that choice only then)
    cyclic_assembly: bool = False  # a C4 or C7 assembly drawn whole that a rosette can draw (cli.ROSETTE_ORDERS):
    # Style.cyclic chooses how it is drawn (the app shows that choice only then)
    lent_across_cut: bool = False  # a strand lent across an unrolled ring's cut (_unwrap): Style.copy_ghosts chooses
    # whether its dashed copies are drawn (the app shows that choice only then)
    lent_apart: dict[str, int] = field(default_factory=dict)  # copy_ghosts off: each strand lent across the cut, drawn
    # past the row's end with no ghost -> the block of `sheets` it pairs in (figcheck.order_checks, provenance)
    font_scale: float = 1.0  # Style.font_scale the figure is drawn with (texts placed in the layout are sized by it)
    loop_stroke: float = 1.0  # a loop's drawn stroke in line_pt: Style.loop_width, RUN_EDGE_W times it under
    # residue colours (its dark edge); for line_band
    helix_matched: dict = field(default_factory=dict)  # {"cut", "joined"}: helix pieces made alike in copies (ss)
    helix_matched_ids: list[str] = field(default_factory=list)  # the helices those cuts and joins made
    na_letter_rows: tuple[float | None, bool] = (None, False)  # where a stem's base letters sit (NT_ROW, 0 on the
    # rails, None: not drawn) and whether single-stranded letters are drawn, as render draws them (na_text_boxes)
    ligands: str = "auto"  # Style.ligands the figure is drawn with (auto hides additives: scope.drawn_scope names them)
    coil_chains: list[tuple[str, int]] = field(default_factory=list)  # (chain, residues) of each protein chain read
    # with no helix or strand: nothing of it can be drawn (scope.drawn_scope says so)
    rosette: int = 0  # a Cn assembly drawn as a rosette (Style.cyclic, plan 21): n, its n-fold axis at page (0, 0)
    rosette_push: float = 0.0  # page units each copy was pushed out from the axis to clear its neighbours
    rosette_strays: list[str] = field(default_factory=list)  # copy elements with no protomer-0 mate, placed alone
    rosette_dropped: list[str] = field(default_factory=list)  # edits (move, swap, rename) naming a copy's element,
    # not applied: a rosette draws every copy as protomer 0 turned (provenance notes them)
    rosette_mates: dict[str, tuple[str, int]] = field(default_factory=dict)  # a rosette's copy element -> (its
    # protomer-0 mate, its protomer k): drawn as the mate turned by frame.angles[k] (route turns the mate's loops)
    cyclic_note: str = ""  # what the notes say of Style.cyclic: why a rosette asked for was drawn unrolled, or how
    # units: each element drawn -> (cx, cy) of the unit (sheet block, row, bundle, helix) it was laid out in and its
    # offset there (build_layout; a rosette records none): a structure drawn like this one takes its homes from them
    units: dict[str, tuple[float, float, float, float]] = field(default_factory=dict)
    like_note: str = ""  # what the notes say of --like (cli.make_layout): drawn like the reference, or why not
    like_extras: list[str] = field(default_factory=list)  # drawn like a reference: the elements it lacks (ids)
    # element id -> colour its mate in the reference has (Like.colour): render.element_colours draws these for the
    # identity colourings (like.IDENTITY_COLOURINGS); data colourings keep this structure's own values
    like_colour: dict[str, str] = field(default_factory=dict)

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
        lx, ly = helix_label_pos(Placed(m.sse, 0.0, 0.0, m.length, m.angle, m.width), k)
        hx, hy = gw * 1.5 + _NAME_PAD * k, gh / 2  # up to three glyphs: α10
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


def _drawn_strands(sheet: Sheet) -> tuple[list[SSE], list[int]]:
    """Every strand a sheet's block draws, its track (Sheet.strands) then its accessory strands, with their
    directions."""
    return [*sheet.strands, *(a.strand for a in sheet.accessory)], [
        *sheet.directions,
        *(a.direction for a in sheet.accessory),
    ]


def _sheet_seen(sheet: Sheet, frame: Frame) -> np.ndarray:
    """The page direction (snapped) of a sheet's strands, each by its sense: its drawn up, strands[0]'s way."""
    every, senses = _drawn_strands(sheet)
    mean_dir = sum(d * s.axis for d, s in zip(senses, every))
    at = np.mean([x.centroid for x in every], axis=0)
    return _snap(frame.direction(at, mean_dir))


def _sheet_order(
    sheet: Sheet, frame: Frame, ref_angle: dict[str, float] | None = None
) -> tuple[list[int], int, np.ndarray]:
    """How a sheet is drawn: the track strands' indices left to right, s = +1 when strands[0] points up (else -1),
    and the projected centroids (Å) of every strand drawn (_drawn_strands: the track's, then the accessory's).
    `ref_angle` (strand id -> drawn angle, a reference's: Like.angle): the first track strand it names points as
    there, the sheet's way up taken from it; the order follows from s as always (the front face, never a mirror)."""
    strands = sheet.strands
    every, _ = _drawn_strands(sheet)
    cents = frame.project(np.array([s.centroid for s in every]), [s.chain for s in every])
    told = next(((k, ref_angle[x.id]) for k, x in enumerate(strands) if x.id in ref_angle), None) if ref_angle else None
    # Both decisions below are made on the grid (_snap): a sheet straddling a dyad in the page (1SLV's, shared by the
    # two copies) runs straight across the page by symmetry, and its end strands, dyad mates, stand level: its
    # height and their step across it are 0 up to the last bits of the element axes, which differ between BLAS
    # kernels (Windows). Such a sheet draws strands[0] up, and its strands in their own order; unless the layout
    # being built turns it down (_TIE_TURN: _both_ways lays out a shared one both ways).
    seen = _sheet_seen(sheet, frame)
    s = 1 if seen[1] > 0 or (seen[1] == 0 and strands[0].id not in _TIE_TURN.get()[0]) else -1  # +1: strands[0] up
    if told is not None:  # strand k drawn up exactly when directions[k] * s > 0
        s = sheet.directions[told[0]] * (1 if np.sin(told[1]) > 0 else -1)
    # Drawing the strands vertical turns the projected sheet in the page plane; the strand order must turn
    # with it (left to right = across the drawn-up direction), or the sheet shows its back face (mirror image).
    up2d = s * seen
    right = np.array([up2d[1], -up2d[0]]) if np.linalg.norm(up2d) > 1e-6 else np.array([1.0, 0.0])
    order = list(range(len(strands)))
    if len(strands) > 1 and float(_snap((cents[len(strands) - 1] - cents[0]) @ right)) < 0:
        order.reverse()
    return order, s, cents


_PIECE_GAP = 0.3  # page units between two pieces of one helix drawn end to end (_bundle_item, _regroup_copies)
PIECE_GAP = _PIECE_GAP  # (its public name: like.like_from lays a split helix's pieces so)
_STACK_PAD = 1.3  # page units between two strands drawn end to end in one column: room for the loops leaving
# their facing ends (0.2 left 5KQJ's B and F loops no route, 1.0 5JB1's; C34 set 1.2 on 5 sample entries). Plan 12
# fix round 3 (accessories beside their partner): swept 1.2-1.6 on the 37 sample entries with accessory strands (5
# views) and the 30QM and 4YBB fixtures; entries newly failing against t12c2b: 1.2: 7 (and 4YBB's run crossing),
# 1.25: 5, 1.3: 8RRX, 1.35: 5JB1, 8RRX, 1.4: 5JB1, 8RRX (and 30QM's beads), 1.45: 7GLY, 7XSP, 8RRX (and 30QM's
# beads), 1.5: 5, 1.6: 5. The stable plateau is 1.3-1.4 (a step changes at most one entry); within it 1.3 alone also
# keeps 6ZS5's disulfide bars no longer with their springs than without (1.35 and 1.4: 0.008 longer)


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


# a stack gap's way out: the longest run of its height beside which the column on one side has no strand. Less
# than _STACK_SHUT both sides (the columns either side standing past the gap's ends but for a sliver: 9P1X 0.15,
# 8OZJ 0.54), the loops leaving its facing ends have no way out, or squeeze through one grid cell together (8OZJ's
# default view): the strand stacked last stands further off, by at most _STACK_NUDGE, meeting nothing and as level
# with its partner as it was, until a side is open by _STACK_OPEN (two loops' grid cells, route.RES, each _INFLATE
# clear of the strands), else by _STACK_SHUT (one cell, 0.55, and a margin), else it stays (random-100 C2b-2)
_STACK_SHUT = 0.6
_STACK_OPEN = 0.85
_STACK_NUDGE = 1.0
_NUDGE_STEP = 0.05
# an accessory that closes a run of free height in a column all round (stacks beside stacks: no run beside it leads
# out, however tall the run) moves by at most this, level with its partner or not, until it closes none (random-100
# held-out 8BCU: G:111-114, level with G:119-122, closed the run four loops needed; 12BM: C:71-79 closed the gap
# between C:84-89 and B:107-116)
_POCKET_REACH = 2.0
_TOUCH = 1e-6  # page units two drawn spans must share to stand level (science._TOUCH)
_CLEAR = 1e-4  # page units past a block's end a cleared height stands, one snap step (_column_heights)


# more room in a stack gap a chain end faces: its stub and N/C letter. Fix round 3: at most this, so 9FWZ's
# D:662-666, end to end with D:634-641, still stands level with its partner D:645-659 (END_LABEL + 0.6 left it above
# D:645-659's end; the D:618-629 loop into that gap then passes the letter, a route over a label)
_END_ROOM = END_LABEL - 0.15


def _sheet_pairs(sheet: Sheet) -> set[frozenset]:
    """The strand-id pairs of a sheet that pair (may stand side by side): its track's neighbours, a closed sheet's
    ends, its chords, and each accessory with its parent and its other partners (Accessory.others)."""
    s = sheet.strands
    out = {frozenset((a.id, b.id)) for a, b in zip(s, s[1:])}
    if sheet.closed and len(s) > 2:
        out.add(frozenset((s[0].id, s[-1].id)))
    out |= {frozenset((a.id, b.id)) for a, b in sheet.chords}
    for a in sheet.accessory:
        out |= {frozenset((a.strand.id, x.id)) for x in (a.parent, *a.others)}
    return out


def _column_heights(
    sheet: Sheet,
    column: dict[int, int],
    s: int,
    ends: dict[str, str] | None = None,
    ghost: int = 0,
) -> dict[int, float]:
    """The heights (page units) of a sheet's accessory strands (`column`: index into _drawn_strands, the ghost as -1,
    -> drawn column; `ghost`: the index of the strand the ghost repeats). The track strands and the ghost keep the row
    (height 0). Each accessory stands in its column (Accessory.column: beside its partner on its 3D side), level with
    its partner where it can: at its partner's height when nothing in that column meets it there, else end to end
    with what does (one column means collinear), on the side its centroid lies on along the sheet's drawn up (its
    strands' mean direction, by their senses, times the flip s), _STACK_PAD apart, _END_ROOM more for each chain end
    facing the gap (`ends`: element id -> "N", "C" or both, the chain ends drawn; 9FWZ's C end, D:662-666; 7JU4's
    two, chain E's N end under chain f's C end, whose letters met with room for one, random-100 C2b-3). One that
    cannot stand level with its partner stands level with another partner of its own (Accessory.others), in a column
    beside that one, when it can; else as near its partner as it can (_unshown reads that off the drawing, and the
    note says "also pairs with"; fix round 2/3 rulings, 7ARD). Wherever it stands, it is never level with a strand
    it does not pair with (_sheet_pairs) in the column on either side of its own (random-100 C2b-S, 9P1X focus:
    a:106-109 under c:2-8 stood beside c:50-61; I:94-97, level with H:85-97, beside I:2-8 and I:29-30): such a height
    is given up for the free one nearest its partner's, in its own column. One whose stack gap stays shut there
    (_STACK_SHUT: the loops leaving its facing ends have no way out) stands on its partner's other side instead, by
    the same end of its partner, when it stands level with its partner there, beside no stranger, its gap open
    (random-100 held-out 2FR2: A:115-118 under A:103-109, beside its partner A:89-100 and A:125-135, both standing
    past the gap's ends; held level with the partner, no nudge could open it). Last, one that closes a run of free
    height in any column (unpocketed) moves clear of it. Snapped."""
    every, senses = _drawn_strands(sheet)
    ends = ends or {}
    n = len(sheet.strands)
    sense = dict(zip(range(len(senses)), senses)) | {-1: senses[ghost]}
    strand = lambda k: every[k] if k >= 0 else every[ghost]  # noqa: E731
    length = lambda k: _strand_length(strand(k))  # noqa: E731
    up3 = _unit(s * sum(d * x.axis for d, x in zip(senses, every)))
    index = {id(x): k for k, x in enumerate(every)}

    def room(lo: int, hi: int) -> float:  # lo stands below hi: lo's top end and hi's bottom end face the gap
        up_lo, up_hi = sense[lo] * s > 0, sense[hi] * s > 0
        facing = int(("C" if up_lo else "N") in ends.get(strand(lo).id, "")) if lo >= 0 else 0
        facing += int(("N" if up_hi else "C") in ends.get(strand(hi).id, "")) if hi >= 0 else 0
        return _STACK_PAD + _END_ROOM * facing  # each chain end facing the gap: its own stub and letter

    height: dict[int, float] = {k: 0.0 for k in range(n)} | ({-1: 0.0} if -1 in column else {})
    occupied: dict[int, list[int]] = {}
    for k in height:
        occupied.setdefault(column[k], []).append(k)

    def spot(k: int, c: int, want: float, way: int = 0) -> float:  # the height nearest `want` where k meets
        # nothing in column c (`way`: past what it meets that way, else on the side it lies on in 3D)
        half = length(k) / 2
        blocks = sorted(
            (height[q] - length(q) / 2 - room(k, q) - half, height[q] + length(q) / 2 + room(q, k) + half, q)
            for q in occupied.get(c, [])
        )
        y = want
        for _ in range(len(blocks) + 1):
            hit = next(((lo, hi, q) for lo, hi, q in blocks if lo < y < hi), None)
            if hit is None:
                return y
            lo, hi, q = hit
            if not way:  # past the end of what meets it, on the side it lies on in 3D; then on the same way
                way = 1 if float(_snap((_snap(strand(k).centroid) - _snap(strand(q).centroid)) @ up3)) >= 0 else -1
            y = hi if way > 0 else lo
        return y

    def level(k: int, q: int, y: float) -> bool:
        return abs(y - height[q]) < (length(k) + length(q)) / 2

    def runs(x: int, own: tuple = ()) -> list[tuple[float, float]]:
        """The free runs of height of column x (`own`: k's (column, bottom, top), not yet placed), the ends open."""
        spans = [(height[q] - length(q) / 2, height[q] + length(q) / 2) for q in occupied.get(x, [])]
        spans += [own[1:]] if own and own[0] == x else []
        out, y = [], -np.inf
        for lo, hi in sorted(spans):
            if lo > y:
                out.append((y, lo))
            y = max(y, hi)
        return out + [(y, np.inf)]

    def way_out(g0: float, g1: float, c: int, need: float, own: tuple = ()) -> bool:
        """Does gap (g0, g1) in column c open so far: a run of free height at least `need` beside it, in the next
        column, that leads on, run beside run, to the sheet's edge or past a column's last strand? A run between two
        strands of the next column shut on its own far side is no way out (random-100 held-out 12BM, 8BCU: stacks
        beside stacks, each gap open only into another closed one). `own`: k's (column, bottom, top), not yet placed."""
        cols = sorted(set(column.values()))
        todo, seen = [(c, g0, g1)], set()
        while todo:
            x, lo, hi = todo.pop()
            for side in (-1, 1):
                at = cols.index(x) + side
                if not 0 <= at < len(cols):
                    return True  # the sheet's edge
                for a, b in runs(cols[at], own):
                    if min(hi, b) - max(lo, a) < need - 1e-9:
                        continue
                    if np.isinf(a) or np.isinf(b):
                        return True  # past the column's last strand
                    if (cols[at], a, b) not in seen:
                        seen.add((cols[at], a, b))
                        todo.append((cols[at], a, b))
        return False

    def shut(k: int, c: int, y: float, need: float) -> list[int]:  # what k stands end to end with, gap shut
        half, out = length(k) / 2, []
        for q in occupied.get(c, []):
            top, bottom = height[q] + length(q) / 2, height[q] - length(q) / 2
            if 0 < (y - half) - top <= room(q, k) + _STACK_NUDGE + 1e-9:
                gap = (top, y - half)
            elif 0 < bottom - (y + half) <= room(k, q) + _STACK_NUDGE + 1e-9:
                gap = (y + half, bottom)
            else:
                continue
            if not way_out(*gap, c, need, (c, y - half, y + half)):
                out.append(q)
        return out

    def meets(k: int, c: int, y: float) -> bool:  # does k at height y come nearer than its room to anything in c?
        half = length(k) / 2
        return any(
            height[q] - length(q) / 2 - room(k, q) - half + 1e-9
            < y
            < height[q] + length(q) / 2 + room(q, k) + half - 1e-9
            for q in occupied.get(c, [])
        )

    pairs = _sheet_pairs(sheet)

    def beside(c: int) -> list[int]:  # the columns standing next to column c: the nearest taken on either side
        taken = sorted(x for x in occupied if x != c and occupied[x])
        lower, higher = [x for x in taken if x < c], [x for x in taken if x > c]
        return lower[-1:] + higher[:1]

    def strangers(k: int, c: int) -> list[int]:  # the strands beside column c that k does not pair with
        me = strand(k).id
        return [
            q
            for x in beside(c)
            for q in occupied[x]
            if strand(q).id != me and frozenset((me, strand(q).id)) not in pairs
        ]

    def clash(k: int, c: int, y: float) -> bool:  # does k at height y in column c stand level with a stranger?
        half = length(k) / 2
        return any(
            min(y + half, height[q] + length(q) / 2) - max(y - half, height[q] - length(q) / 2) > _TOUCH
            for q in strangers(k, c)
        )

    def clear(k: int, c: int, want: float) -> float:
        """The height nearest `want` where k meets nothing in column c and no stranger beside it (the lower on a
        tie)."""
        half = length(k) / 2
        cands = [want]  # past each block's ends by _CLEAR, so snapping never brings it back in
        for q in occupied.get(c, []):
            cands += [height[q] - length(q) / 2 - room(k, q) - half, height[q] + length(q) / 2 + room(q, k) + half]
        for q in strangers(k, c):
            cands += [height[q] - length(q) / 2 - half, height[q] + length(q) / 2 + half]
        ys = sorted({float(_snap(y + _CLEAR * np.sign(y - want))) for y in cands}, key=lambda y: (abs(y - want), y))
        free = next((y for y in ys if not meets(k, c, y) and not clash(k, c, y)), None)
        # the farthest block end is always free (each block is an open interval, cleared by _CLEAR past its end)
        assert free is not None, f"no free height for {strand(k).id} in column {c}"
        return free

    def opened(k: int, c: int, y: float, p: int) -> float:
        """k's height, moved off the one strand it stands end to end with across a gap shut (_STACK_SHUT) until
        the gap opens to a side by _STACK_OPEN, else by _STACK_SHUT, by at most _STACK_NUDGE, meeting nothing else
        and as level with its partner p as it was; else as it stands."""
        walls = shut(k, c, y, _STACK_SHUT)
        if len(walls) != 1:
            return y
        way, was, some = (1 if y > height[walls[0]] else -1), level(k, p, y), None
        for step in range(1, round(_STACK_NUDGE / _NUDGE_STEP) + 1):
            z = float(_snap(y + way * step * _NUDGE_STEP))
            if meets(k, c, z) or level(k, p, z) != was or clash(k, c, z):
                break
            if not shut(k, c, z, _STACK_OPEN):
                return z
            if some is None and not shut(k, c, z, _STACK_SHUT):
                some = z
        return y if some is None else some

    def closed(own: tuple = ()) -> int:
        """How many runs of free height between two strands of a column (`own`: k's (column, bottom, top) among
        them) have no way out (way_out, by _STACK_SHUT): the loops leaving the ends facing one have no route."""
        return sum(
            not way_out(a, b, x, _STACK_SHUT, own)
            for x in sorted(set(column.values()))
            for a, b in runs(x, own)
            if not (np.isinf(a) or np.isinf(b))
        )

    def unpocketed(k: int, c: int, y: float) -> float:
        """k's height, moved when standing at y it closes a run of free height that had a way out, or faces one
        closed (closed: stacks beside stacks, however tall the run), to the nearest height within _POCKET_REACH,
        _NUDGE_STEP at a time, up first on a tie, that closes none, meeting nothing and level with no stranger,
        level with its partner or not (the notes then name the pairing); else as it stands."""
        half, base = length(k) / 2, closed()
        if closed((c, y - half, y + half)) <= base:
            return y
        for step in range(1, round(_POCKET_REACH / _NUDGE_STEP) + 1):
            for way in (1, -1):
                z = float(_snap(y + way * step * _NUDGE_STEP))
                if not meets(k, c, z) and not clash(k, c, z) and closed((c, z - half, z + half)) <= base:
                    return z
        return y

    for i, a in enumerate(sheet.accessory):
        k, p = n + i, index[id(a.parent)]
        y = spot(k, column[k], height[p])
        if not level(k, p, y):
            for x in a.others:  # level with another partner of its own, in a column beside it, if it can
                q = index.get(id(x))
                if q is None or q not in height:
                    continue
                for c in (column[q] - 1, column[q] + 1):
                    z = spot(k, c, height[q])
                    if level(k, q, z) and not clash(k, c, z):
                        column[k], y = c, z
                        break
                else:
                    continue
                break
        y = opened(k, column[k], y, p)
        if clash(k, column[k], y):
            y = clear(k, column[k], height[p])
        if shut(k, column[k], y, _STACK_SHUT):  # still no way out of its stack gap: beside its partner's other side
            c = 2 * column[p] - column[k]
            z = opened(k, c, spot(k, c, y, 1 if y >= height[p] else -1), p)  # at its partner's same end
            if c != column[k] and level(k, p, z) and not clash(k, c, z) and not shut(k, c, z, _STACK_SHUT):
                column[k], y = c, z
        y = unpocketed(k, column[k], y)
        height[k] = float(_snap(y))
        occupied.setdefault(column[k], []).append(k)
    return {k: h for k, h in height.items() if k >= n}


def _sheet_item(
    sheet: Sheet,
    frame: Frame,
    pitch: float = PITCH,
    over: list[SSE] = (),
    helix_w: float = HELIX_W,
    ends: dict[str, str] | None = None,
    ref_angle: dict[str, float] | None = None,
) -> _Item:
    """A sheet as one block of strands at `pitch`, each pointing as it runs relative to strands[0]; a closed
    barrel's first strand repeated as a ghost after the last. Each accessory strand (Sheet.accessory) stands in its
    column (Accessory.column, turned with the track when the track is drawn right to left), stacked with what stands
    there (_column_heights) or past the track's end: beside its parent, never between two strands that pair. `over`
    (helices packed on its face, in sequence): each drawn over the strands where it packs (_over_sheet), unless it
    would lie on one drawn there before it. `ends`: the chain ends drawn. `ref_angle`: a reference's strand angles
    (_sheet_order)."""
    strands = sheet.strands
    n = len(strands)
    order, s, cents = _sheet_order(sheet, frame, ref_angle)
    forward = order[0] == 0
    drawn = {k: d for d, k in enumerate(order)}  # index into _drawn_strands (the ghost: -1) -> drawn column
    drawn.update({n + i: a.column if forward else n - 1 - a.column for i, a in enumerate(sheet.accessory)})
    if sheet.closed:
        drawn[-1] = n  # the ghost, after the last
    heights = _column_heights(sheet, drawn, s, ends, order[0]) if sheet.accessory else {}  # may move some
    slot_of = {c: i for i, c in enumerate(sorted(set(drawn.values())))}
    slots = len(slot_of)
    members, offsets = [], []
    every, senses = _drawn_strands(sheet)

    def put(k: int) -> None:
        if k < 0:  # the barrel wraps: its first strand pairs again with the last
            f = members[0]
            members.append(Placed(f.sse, 0, 0, f.length, f.angle, f.width, ghost=True))
        else:
            up = senses[k] * s > 0
            members.append(Placed(every[k], 0, 0, _strand_length(every[k]), np.pi / 2 if up else -np.pi / 2, pitch))
        offsets.append(((slot_of[drawn[k]] - (slots - 1) / 2) * pitch, heights.get(k, 0.0)))

    for k in [*order, *([-1] if sheet.closed else []), *range(n, len(every))]:  # track, ghost, accessory strands
        put(k)
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
    return _Item(members, offsets, home[0], home[1], home, *_anchor(every))


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
    angle = _half_turn(_snap(np.arctan2(axis @ up3, axis @ across3)))
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


def can_be_row(sheet: Sheet) -> bool:
    """Can the sheet be drawn as a βα row at all (_row_shape)? Pure; like.like_from asks it of a mated sheet."""
    return _row_shape(sheet)


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
    sheet: Sheet,
    units: list[BAUnit],
    frame: Frame,
    pitch: float = PITCH,
    helix_w: float = HELIX_W,
    ref_angle: dict[str, float] | None = None,
) -> _Item:
    """An all-parallel α/β sheet (template_applies) as one row, as the hand-drawn TIM barrel and Rossmann figures
    draw it: the strands in sheet order and pointing as _sheet_item draws them (the front face), and after each
    strand, when an in-row unit joins it to the next, that unit's longest helix (the first on a tie), true length,
    running against the strands; a barrel's closing helix between the last strand and the ghost. Each footprint
    (labels included) starts _BA_GAP after the last one ends, all centred on one line, the row centred on the item.
    Other helices are placed as free helices. Home: the strands' projected mean, as _sheet_item. `ref_angle`: a
    reference's strand angles (_sheet_order), the row then pointing as the reference's does, never turned again."""
    strands = sheet.strands
    order, s, cents = _sheet_order(sheet, frame, ref_angle)
    told = bool(ref_angle) and any(x.id in ref_angle for x in strands)
    if not told and strands[order[0]].start > strands[order[-1]].start:  # the row reads in sequence: turned
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
    sse: SSE,
    frame: Frame,
    angle_mode: str = "snap",
    width: float = HELIX_W,
    rise: float = HELIX_RISE,
    angle: float | None = None,
) -> _Item:
    """Always true length. angle_mode: 'tilted' as seen in the view; 'snap' to the nearest of up, down, left,
    right (the topology-diagram convention); 'upright' up or down like a strand (stack mode). `angle`: drawn so
    instead (a reference's, Like.angle)."""
    au, av = frame.direction(sse.centroid, sse.axis)
    length = max(width, (len(sse) - 1) * rise * SCALE)  # true length even when tilted: equal pitch everywhere
    if angle is None:
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


def _kink_item(
    group: list[SSE], frame: Frame, angle_mode: str, width: float = HELIX_W, ref_angle: dict[str, float] | None = None
) -> _Item:
    """The pieces of a kinked helix as one unit, end to end, so the corner reads as one helix turning (2HHB α's AB
    corner, which an element seen between the pieces kept five widths apart): each piece as _helix_item draws it
    (true length, its own angle). Each next piece starts where the loop from the last one's C port turns once into
    it: a piece across the last starts just beyond its end and just off its axis, both legs a routing cell clear of
    the two bodies (_KINK_ROOM); a piece running on starts a little further along; a piece turned back runs beside
    the last, MARGIN apart, on the side it is seen. The unit sits at the mean of the pieces' projected places,
    offsets on the grid. `ref_angle`: the pieces a reference draws at its angles (Like.angle)."""
    ref_angle = ref_angle or {}
    parts = [_helix_item(s, frame, angle_mode, width, angle=ref_angle.get(s.id)) for s in group]
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
    """The drawn angle of an element seen at angle `seen` (page-up component av). Wherever the sign of av decides
    (up or down; pointing left as +pi or -pi) it is av on the grid, so noise either side of zero, which depends on
    the CPU, decides nothing: a helix along the page points up, or left as +pi."""
    up = float(_snap(av)) >= 0
    if mode == "upright":
        return np.pi / 2 if up else -np.pi / 2
    if mode == "snap":
        k = float(np.round(seen / (np.pi / 2)))
        if abs(k) == 2:  # left: +pi or -pi as av's side (on the grid)
            k = 2.0 if up else -2.0
        return k * (np.pi / 2) + 0.0
    return _half_turn(_snap(seen))  # as seen, on the grid (radians)


def _bundle_item(
    group: list[SSE], frame: Frame, angle_mode: str, width: float, ref_angle: dict[str, float] | None = None
) -> _Item:
    """Helices of several chains packed along their length, drawn as one unit: one column per helix (kinked
    pieces of one helix share a column, end to end), columns side by side in their order across the bundle at
    their true spacing (at least a label's width apart), staggered along the shared axis as in 3D;
    antiparallel helices point the other way. A piece drawn the other way from the last faces it end to end, both
    ends' loops leaving into the gap between (0.3 left them no way out: random-100 C2b-2, 6JLN's E:42-47, 7ARD's
    H:121-151). One that is no kinked piece of it by features' rule (both at least KINK_SHORT long, their axes
    bent more than KINK_COS: 7ARD's, cos -1.0, a hairpin) is another helix: it takes a column of its own. A kinked
    piece drawn the other way (a short one, whose axis is no guide: 6JLN's six residues) stays in the column,
    _STACK_PAD from it, as two strands end to end in a sheet's column. `ref_angle` (Like.angle): the bundle turned
    so the first helix (N to C) a reference draws, of those running along it if any, points as there."""
    from .features import KINK_COS, KINK_GAP, KINK_SHORT

    ref = group[0].axis
    mean = sum(np.sign(float(s.axis @ ref)) * s.axis * len(s) for s in group)
    m = mean / np.linalg.norm(mean)
    centre = np.mean([s.centroid for s in group], axis=0)
    a2 = frame.direction(centre, m)
    up = lambda s: float(s.axis @ m) >= 0  # noqa: E731 - drawn along the bundle (else the other way), as below
    told = sorted((not up(s), s.start) + (s.id,) for s in group if s.id in ref_angle) if ref_angle else []
    if not told:
        angle = _angle(float(np.arctan2(a2[1], a2[0])), float(a2[1]), angle_mode)
    else:  # each helix is drawn at the bundle's angle, or half a turn round when it runs the other way
        back, _, sid = told[0]
        angle = _wrap_angle(ref_angle[sid] - np.pi) if back else ref_angle[sid]
    d2, side2 = np.array([np.cos(angle), np.sin(angle)]), np.array([-np.sin(angle), np.cos(angle)])
    columns: list[list[SSE]] = []

    def hairpin(a: SSE, b: SSE) -> bool:  # b drawn the other way from a, and no kinked piece of it
        bent = float(_snap(float(a.axis @ b.axis))) <= KINK_COS
        return up(a) != up(b) and min(len(a), len(b)) >= KINK_SHORT and bent

    for s in sorted(group, key=lambda s: (s.chain, s.start)):
        last = columns[-1][-1] if columns else None
        if last is not None and last.chain == s.chain and s.start - last.end - 1 <= KINK_GAP and not hairpin(last, s):
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
        last_end, prev = -np.inf, None
        for s in sorted(col, key=lambda s: float((s.centroid - centre) @ m)):
            length = max(width, (len(s) - 1) * HELIX_RISE * SCALE)
            gap = _PIECE_GAP if prev is None or up(s) == up(prev) else _STACK_PAD  # ends facing: room for loops
            along = max(float((s.centroid - centre) @ m) * SCALE, last_end + gap + length / 2)  # kinked pieces in a row
            last_end, prev = along + length / 2, s
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
    shared: set[int] = frozenset(),
    again: bool = False,
) -> str:
    """Push overlapping items apart, then slide any left over clear. With `ties` (_symmetry_ties) the pushes
    keep every protomer drawn alike (_push_apart), and that result is kept if it clears every overlap within
    _SEP_ORBIT_SWEEPS sweeps and its box takes at most `room` times the area of moving the items one by
    one from the same start (widening the spread can otherwise cost far more room); else the one-by-one
    result is used. Consecutive pieces of one bent duplex (`joints`, _joints) in a short chain (_rigid_joints)
    are moved together and never pushed apart. Items are pushed off the `fixed` boxes they must keep off (Walls),
    which stay where they are. `mirror`: the copies are their mates' mirror images about that page x (_push_apart);
    then the tied pushes run as many passes as the one-by-one ones and may take up to _SEP_MIRROR_ROOM. `level`
    (bool per item): items that keep their height (the transmembrane elements of a membrane view), only ever pushed
    and slid sideways (_push_apart, _slide_clear). `shared`: items holding elements of several protomers with no
    copy of their own (a dimer's continuous sheet, 1SLV): when the tied pushes fail, they are run again with these
    yielding to the tied items, which keep their copies alike (_push_apart), and sliding clear last
    (_slide_clear). With `again` (build_layout's `alike_first`), tied pushes not kept get a second chance
    (_tied_again): every item with no copy yields, and the result is kept if it takes at most _SEP_TIE_ROOM times
    the one-by-one room (both measured after their slides). Returns "alike" when the tied pushes are kept (the
    copies drawn alike), "wide" when the second chance keeps them, else "alone" (moved one by one; also with no
    ties)."""
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
    lone: list[int] = []
    if ties and shared and not together:  # the tied pushes fail: again, the shared items yielding to the copies
        xy, lone = start.copy(), sorted(shared)
        together = _push_apart(xy, halves, mids, home, ties, sweeps, apart, fixed, mirror, level, shared)
        if not together and _lone_left(xy, halves, mids, apart, set(lone)):
            together = True  # what is left overlaps only the shared items: they slide clear, the copies stay alike
    if not together or _box_area(xy, halves, mids) > room * _box_area(alone, halves, mids):
        xy, lone = alone, []
    for it, (x, y) in zip(items, xy.tolist()):
        it.cx, it.cy = x, y
    _slide_clear(items, joints, fixed, level, lone)
    if ties and xy is not alone:
        return "alike"
    if ties and again and _tied_again(items, ties, start, halves, mids, home, joints, apart, fixed, mirror, level):
        return "wide"
    return "alone"


def _one_per_mate(
    items: list[_Item], ties: dict[int, tuple[int, int]], helices: bool = True
) -> dict[int, tuple[int, int]]:
    """`ties` with one copy item per mate and protomer: of several tied to one mate as one protomer (_symmetry_ties'
    partial mates: a copy's sheet DSSP split into two blocks where its mate's is one, a helix cut in two), the one
    with the most elements (ghosts aside), then the larger footprint, then the first, stays tied; the others are
    items of their own. Two items held to one place in one orbit overlap and the pushes, which move an orbit as one,
    can only part them by halves, each off its mate's place (random-100 held-out A15: 2BBJ E~A 1.42; 5K7Z, once
    cli._side_by_side ties its dimers, C~A 1.10 without this; 9GTP's 1K~1l and S~X). With `helices` False
    (_regroup_copies runs next and joins a helix piece let go to its helix: 7QUP 9C) only items holding strands are
    judged; helix pieces stay tied."""

    def sheet(i: int) -> bool:
        return helices or any(isinstance(m, Placed) and m.sse.kind == "E" and not m.ghost for m in items[i].members)

    best: dict[tuple[int, int], tuple[tuple[int, float], int]] = {}
    for i, (j, k) in sorted(ties.items()):
        if not sheet(i):
            continue
        size = (
            sum(isinstance(m, Placed) and not m.ghost for m in items[i].members),
            float(np.prod(items[i].half)),
        )
        if (j, k) not in best or size > best[(j, k)][0]:
            best[(j, k)] = (size, i)
    keep = {i for _, i in best.values()}
    return {i: t for i, t in ties.items() if i in keep or not sheet(i)}


def _tied_again(
    items: list[_Item],
    ties: dict[int, tuple[int, int]],
    start: np.ndarray,
    halves: np.ndarray,
    mids: np.ndarray,
    home: np.ndarray,
    joints: set[tuple[int, int]],
    apart: np.ndarray,
    fixed: Walls | None,
    mirror: float | None,
    level: np.ndarray | None,
) -> bool:
    """_separate's second chance for copies drawn alike, the items already placed one by one (and slid clear). The
    tied pushes run again from `start` for _SEP_SWEEPS passes with every item that has no copy (no tie holds it or
    holds it as a mate: elements one copy lacks, a bundle with a chain outside every protomer, a shared sheet)
    yielding to the orbits (_push_apart's `shared`), and sliding clear last. Pushed against an orbit, such an item
    would otherwise widen the whole spread at each push (random-100 C2b-3: 6JLN's from 80 to 215 page units).
    Nucleic-acid stems never yield: they are placed by the protein they touch, and pushed off it they lose the shape
    the structure has (30QM's Spearman 0.496 -> 0.436, held at 0.446). What the pushes leave overlapping slides
    clear with each orbit moved as one (_slide_clear's `ties`: a jam the pushes trade back and forth between items
    and the walls of a trace, the same in every copy, 13FZ's A and C), the yielding items last. Kept, and True, when
    nothing is left overlapping and the result takes at most _SEP_TIE_ROOM times the one-by-one room; else the items
    are put back. Not where nucleic-acid stems are drawn: they never yield, and orbits slid round them leave the
    protein off the stems it holds (30QM: four more crossings of its RNA over protein, Procrustes 0.794 -> 0.795);
    there the second chance stands only when what the pushes leave overlaps only the yielding items, as before."""
    free = {
        i
        for i in set(range(len(items))) - set(ties) - {j for j, _ in ties.values()}
        if not any(isinstance(m, PlacedStem) for m in items[i].members)
    }
    stems = any(isinstance(m, PlacedStem) for it in items for m in it.members)
    xy = start.copy()
    done = _push_apart(xy, halves, mids, home, ties, _SEP_SWEEPS, apart, fixed, mirror, level, free)
    if stems and not done and not _lone_left(xy, halves, mids, apart, free):
        return False
    alone = [(it.cx, it.cy) for it in items]
    room = _box_area(np.array(alone, float).reshape(-1, 2), halves, mids)
    for it, (x, y) in zip(items, xy.tolist()):
        it.cx, it.cy = x, y
    _slide_clear(items, joints, fixed, level, sorted(free), None if stems else ties, mirror)
    now = np.array([(it.cx, it.cy) for it in items], float).reshape(-1, 2)
    clear = stems or not len(_overlapping(now, halves, mids, apart))
    if clear and _box_area(now, halves, mids) <= _SEP_TIE_ROOM * room:
        return True
    for it, (x, y) in zip(items, alone):
        it.cx, it.cy = x, y
    return False


def _lone_left(xy, halves, mids, apart, lone: set[int]) -> bool:
    """Does every overlap left after the tied pushes (_push_apart) involve one of the `lone` (shared) items?"""
    pairs = _overlapping(xy, halves, mids, apart)
    return bool(len(pairs)) and all(i in lone or j in lone for i, j in pairs.tolist())


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
    shared: set[int] = frozenset(),
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
            lone = [q for q in (i, j) if q in shared and cid[i] != cid[j]]
            if ties and len(lone) == 1 and protomer[j if lone[0] == i else i] is not None:
                # an item of several protomers with no copy (`shared`: a dimer's sheet, 1SLV) yields to an orbit, the
                # whole way: moving the orbit moves its copies too, each against its own neighbours, and they part
                q, sign = lone[0], (-1 if lone[0] == i else 1)
                if px < py:
                    xy[group[q], 0] += sign * ((1 if dx > 0 else -1) if dx != 0 else 1) * (px + _SEP_PUSH)
                else:
                    xy[group[q], 1] += sign * ((1 if dy > 0 else -1) if dy != 0 else 1) * (py + _SEP_PUSH)
                continue
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
    last: list[int] = (),
    ties: dict[int, tuple[int, int]] | None = None,
    mirror: float | None = None,
) -> None:
    """Guarantee: each item in turn, if it still overlaps one already settled (other than a piece it is joined
    to, `joints`) or a `fixed` box it must keep off, takes the nearest clear spot straight right, left, up or down
    (in that order on a tie), in _SLIDE_STEP steps; one that keeps its height (`level`) right or left only. Pieces
    of a bent duplex `joints` holds together (_rigid_joints) slide as one, their offsets kept, each exempt from the
    others (random-100 C2b-4: slid one at a time, 7Y7C's P138c2 landed on P138c1 and P85a2 28 units from P85a1); the
    group settles when its first item comes up in that order (so a group with one item in `last` and one not settles
    early; unreachable today: stems never yield in _tied_again, so none is in `last`). The items in `last` (with no
    symmetry copies, _separate) are settled after every other, so they move, not the copies (1SLV: the dimer's shared
    sheet between the two copies' eta2 helices). With `ties` (_symmetry_ties) an item and its symmetry copies (its
    orbit) slide as one, by the same step, so the copies stay alike (a copy of a `mirror` the other way across);
    they may still overlap one another, which the caller checks."""
    halves = [it.half for it in items]
    mids = [it.mid for it in items]
    boxes = [] if fixed is None else np.asarray(fixed.boxes, float).tolist()
    walls = [
        (SimpleNamespace(cx=(x0 + x1) / 2, cy=(y0 + y1) / 2), ((x1 - x0) / 2, (y1 - y0) / 2), (0.0, 0.0))
        for x0, y0, x1, y1 in boxes
    ]
    later = set(last)
    order = [i for i in range(len(items)) if i not in later] + [i for i in range(len(items)) if i in later]
    links = sorted(joints) + [(i, j) for i, (j, _) in sorted((ties or {}).items())]
    cid = _components(len(items), links).tolist()
    turn = [1.0] * len(items)  # a mirrored copy slides the other way across (as _push_apart's sx)
    if mirror is not None:
        for i, (_, k) in (ties or {}).items():
            turn[i] = -1.0 if k > 0 else 1.0
    done: set[int] = set()  # the groups settled
    settled: list[int] = []
    for i in order:
        if cid[i] in done:
            continue
        group = [j for j in order if cid[j] == cid[i]]  # i alone, or its rigid chain of duplex pieces
        start = [(items[j].cx, items[j].cy) for j in group]
        others = [(items[j], halves[j], mids[j]) for j in settled]
        mine = {j: [w for k, w in enumerate(walls) if fixed.hard[j, k]] if walls else [] for j in group}
        pairs = [(j, o) for j in group for o in mine[j] + others]
        k = 0
        sideways = level is not None and any(level[j] for j in group)
        ways = ((1, 0), (-1, 0)) if sideways else ((1, 0), (-1, 0), (0, 1), (0, -1))
        while any(_collide(items[j], halves[j], mids[j], *o) for j, o in pairs):
            k += 1
            for dx, dy in ways:
                for j, (x0, y0) in zip(group, start):
                    items[j].cx, items[j].cy = x0 + k * dx * turn[j] * _SLIDE_STEP, y0 + k * dy * _SLIDE_STEP
                if not any(_collide(items[j], halves[j], mids[j], *o) for j, o in pairs):
                    break
        settled += group
        done.add(cid[i])


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


def _x_overlaps(Q, want, reach=None) -> tuple[np.ndarray, np.ndarray]:
    """The pairs (i, j), i < j, whose footprints (centres Q) overlap along x: want[i, j, 0] > |Q_j - Q_i| in x, in
    row-major order, as np.nonzero(np.triu(...)) gives them, element for element (the summation order of the terms
    over them moves L-BFGS minima). `reach` (want[..., 0].max(axis=1), once per want): only the items within an
    item's reach in x are compared with it, found by sorting, not all n x n (plan 15 Task 7: 8D3D's 896 items)."""
    n = len(Q)
    if reach is None:
        return np.nonzero(np.triu(want[..., 0] > np.abs(Q[None, :, 0] - Q[:, None, 0]), 1))
    x = Q[:, 0]
    order = np.argsort(x, kind="stable")
    xs = x[order]
    pad = reach * (1 + 1e-9) + 1e-9  # a wider window than the test: the test below decides
    lo, hi = np.searchsorted(xs, x - pad, "left"), np.searchsorted(xs, x + pad, "right")
    counts = hi - lo
    i = np.repeat(np.arange(n), counts)
    j = order[np.arange(int(counts.sum())) - np.repeat(np.cumsum(counts) - counts - lo, counts)]
    hit = (j > i) & (want[i, j, 0] > np.abs(x[j] - x[i]))
    i, j = i[hit], j[hit]
    k = np.lexsort((j, i))
    return i[k], j[k]


def _energy(flat, home, anchor_w, want, pull, loops, lam, mids=None, fixed=None, reach=None):
    """Layout energy and gradient over only the overlapping boxes and the contact pairs `pull` = (i, j, w), i < j.
    Same terms as _energy_dense. `reach`: see _x_overlaps."""
    P = flat.reshape(-1, 2)
    grad = 2 * anchor_w * (P - home)
    E = float((anchor_w[:, 0] * ((P - home) ** 2).sum(1)).sum())
    Q = P if mids is None else P + mids  # footprint centres
    E += _fixed_overlap(Q, fixed, lam, grad)
    oi, oj = _x_overlaps(Q, want, reach)
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
    ~1e-12 per value); exact libm values such as cos 0 = 1 are safe. A value a hair either side of zero snaps to
    0.0, never -0.0 (+ 0.0): the two compare equal but sin, atan2 and the layout signature tell them apart, and
    which one noise gives depends on the CPU (1LMB's duplex drawn along the page, angle 0, under AVX off)."""
    return np.round(np.asarray(v, float) * _GRID) / _GRID + 0.0


def snap(v) -> np.ndarray:
    """Page values on the layout's grid (_snap), for places computed outside layout.py (like.like_from)."""
    return _snap(v)


_HALF_TURN = min(np.pi, round(np.pi * _GRID) / _GRID)  # the lesser of pi and pi on the grid


def _half_turn(angle) -> float:
    """A snapped angle (or a multiple of pi/2) in (-pi, pi]: -0.0 is 0.0 and -pi (or pi on the grid, negated) is
    +pi. atan2 of a direction along the page gives +-0.0 or +-pi by the sign of a component that is zero up to
    noise (the last bits of np.linalg.svd or a dot product), and which one depends on the CPU."""
    a = float(angle) + 0.0
    return -a if a <= -_HALF_TURN else a


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
    held: dict[int, tuple[float, float]] | None = None,
    held_k: float = _K_LIKE,
    settled: frozenset[frozenset[str]] | None = None,
) -> None:
    """Minimise the layout energy (_energy) from the projected positions, raising the overlap penalty in
    stages. The result becomes each item's new home; _separate then guarantees no overlap remains. Every input
    is snapped to the grid (_snap) and _lbfgs runs without BLAS, so the minimum is the same on every machine.
    `extra`: further springs (item i, port offset, item j, port offset, weight, slack), the nucleic-acid strands
    from stem to stem (_stem_springs). `joints` (_joints): pairs of items kept apart without MARGIN, so their
    boxes may touch and no element fits between them for long. `fixed` (Walls): boxes that stay where they are (the
    nucleotides of unpaired strands' traces), the items kept MARGIN off them, firmly or (an element touching that
    trace) not at all. `mirror` (page x of a dyad in the page): each copy of `ties` sits at its mate's mirror image
    about it, at the same height, instead of `spread` along.
    `held` (item index -> home, plan 13a: items drawn like a reference, at their mates' places in its finished
    layout, a copy's tie aside): each starts at its home, held there by a spring `held_k` per element (_K_LIKE) in
    place of _K_HOME, and is moved as a displacement from it (the home itself is not snapped, so an item nothing
    moves stays exactly where the reference has it). No contact pull and no loop spring acts between two held items
    (a loop only this structure has, round elements the reference has and it lacks, is routed as it is: pulled, it
    dragged 2HHB α's E helix off β's place, median deviation 0.25 -> 1.41 page units); no disulfide spring either
    where the reference bridges their mates the same way (`settled`: element id pairs, Like.settled; None: every
    pair), but a disulfide only this structure has still pulls, and when stretched may start again from beside its
    partner, a held item too (a bond drawn across the page is wrong whatever the reference)."""
    n = len(items)
    held = {i: h for i, h in (held or {}).items() if i not in (ties or {})}

    def known(i: int, j: int, a: str, b: str) -> bool:  # a bridge between two held items the reference has too
        return i in held and j in held and (settled is None or frozenset((a, b)) in settled)

    W = _contact_weights(items, contacts)
    if held:  # two items held where a reference has them: no pull between them
        both = np.array(sorted(held), int)
        W[np.ix_(both, both)] = 0.0
    ci, cj = np.nonzero(np.triu(W > 0, 1))
    pull = (ci, cj, W[ci, cj])
    springs = _loop_springs(items, sses)
    if held:  # every loop between two held items is the reference's to place
        springs = [x for x in springs if not (x[0] in held and x[2] in held)]
    weight = [_W_SEQ] * len(springs)
    slack = [_LOOP_SLACK] * len(springs)
    for i, oi, j, oj, w, sl in extra:  # nucleic-acid strands from stem to stem
        if i in held and j in held:
            continue
        springs.append((i, oi, j, oj))
        weight.append(float(_snap(w)))
        slack.append(float(_snap(sl)))
    owner = _owner(items)
    bridges_used = []
    for a, ta, b, tb in bridges:  # disulfides: the two cysteines themselves pulled together
        if a in owner and b in owner and owner[a][0] != owner[b][0]:
            if known(owner[a][0], owner[b][0], a, b):
                continue
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
    if held:  # at the reference's places, unsnapped, and held there harder
        rows = np.array(sorted(held), int)
        home[rows] = np.array([held[i] for i in sorted(held)], float)
        anchor_w[rows] = held_k * np.array([[_mass(items[i])] for i in sorted(held)], float)
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
    at = np.zeros((len(free), 2))  # a held item's z is its displacement from its home
    for i, xy in held.items():
        at[slot[i]] = xy

    def expand(z):
        if held:
            P = (z[:-1].reshape(-1, 2) + at)[base]
        else:
            P = z[:-1].reshape(-1, 2)[base].copy()
        if mirror is not None:  # protomer 1 is protomer 0 seen in the mirror line
            P[copy_k > 0, 0] = 2 * mirror - P[copy_k > 0, 0]
            return P
        P[:, 0] += copy_k * z[-1]  # protomer k sits k spreads to the right of protomer 0
        P[:, 1] += copy_k * lift  # and, along a filament, k rises higher
        return P

    def tied(z):  # one position per orbit of symmetry copies, plus the spread between protomers
        evals[0] += 1
        if n > _PAIR_ENERGY_ABOVE:
            E, g = _energy(expand(z).ravel(), home, anchor_w, want, pull, loops, lam, mids, walls, reach)
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

    evals = [0]  # energy evaluations so far
    reach = want[..., 0].max(axis=1) if n > _PAIR_ENERGY_ABOVE else None  # _x_overlaps
    lam = _LAMBDAS[0]
    start = _snap([[items[i].cx, items[i].cy] for i in free]).reshape(-1, 2)
    for i in held:
        start[slot[i]] = 0.0
    z = solve(np.append(start.ravel(), spread))
    lam = _LAMBDAS[-1]
    best = tied(z)[0]
    mass = np.array([_mass(it) for it in items])
    restart_from, starts, unstarted = evals[0], 0, 0  # the second starts' budget (_SS_RESTART_PAIRS), debug counts
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
        if (evals[0] - restart_from) * n * n >= _SS_RESTART_PAIRS:
            unstarted += 1
            continue
        starts += 1
        for side in ((0.0, -1.0), (-gap_x, 0.0), (gap_x, 0.0), (0.0, gap_y), (0.0, -gap_y)):
            # the lighter element brought beside its partner: just below the cysteine, or clear of its footprint
            trial = _snap(z)  # every start on the grid
            trial[2 * slot[mover] : 2 * slot[mover] + 2] = _snap(P[anchor] + oa - om + np.array(side))
            if mover in held:  # a held item's z is its displacement from its home
                trial[2 * slot[mover] : 2 * slot[mover] + 2] -= at[slot[mover]]
            trial = solve(trial, _LAMBDAS[1:])
            lam = _LAMBDAS[-1]
            energy = tied(trial)[0]
            if energy < best:
                z, best = trial, energy
    if unstarted:
        log.debug(
            "%d stretched disulfides started again, %d not for the budget (%d items, %d pair evaluations)",
            starts,
            unstarted,
            n,
            (evals[0] - restart_from) * n * n,
        )
    x = expand(z).ravel()
    for it, p in zip(items, x.reshape(n, 2)):
        it.cx, it.cy = float(p[0]), float(p[1])
        it.home = (it.cx, it.cy)


_STEM_ROOM = 1.5  # the room a nucleic acid's stems take on the page: their footprints' area (with half a MARGIN round
# each) times this, about what rectangles packed with gaps fill. Measured on the NA-only entries of the NA campaign and
# random-100 (random-100 C2b-4): 3OWW's Spearman 0.484 -> 0.651, 1GID 0.601 -> 0.737, 3IGI 0.83 -> 0.856, the rest
# within 0.001 (1.2: 4LVV 0.961 -> 0.886; 2.0: 4P9R 0.82 -> 0.798)


def _room_for_stems(items: list[_Item]) -> None:
    """A nucleic acid drawn with no protein element (every item a stem; a protein with no helix or strand, drawn as
    coil only, does not count): its projected places spread about their middle, all
    by one factor, until the box round their footprints holds _STEM_ROOM times the footprints' area, before _attract
    pulls and parts them. Packed tighter than that, the overlaps were cleared by pushes each stem took its own way, and
    the stems' order was lost (random-100 C2b-4: 3OWW's ten stems in a box of 0.98 times their footprints' area,
    Spearman 0.792 as projected, 0.493 after _attract). One factor keeps every ratio of the projected distances, so
    their rank. Nothing moves when the room is there, or when any item holds a protein element."""
    if not items or not all(isinstance(m, PlacedStem) for it in items for m in it.members):
        return
    q = _snap([it.home for it in items]) + _snap([it.mid for it in items])
    h = _snap([it.half for it in items]) + MARGIN / 2
    need = float(sum(4.0 * float(x) * float(y) for x, y in h.tolist()))
    lo, hi = (q - h).min(axis=0), (q + h).max(axis=0)
    have = float((hi[0] - lo[0]) * (hi[1] - lo[1]))
    s = float(_snap(math.sqrt(_STEM_ROOM * need / have))) if have > 0 else 1.0
    if s <= 1.0:
        return
    c = _snap(q.sum(axis=0) / len(items))
    for it in items:
        hx, hy = (float(v) for v in _snap(c + s * (np.asarray(it.home, float) - c)))
        it.cx, it.cy = float(_snap(it.cx + hx - it.home[0])), float(_snap(it.cy + hy - it.home[1]))
        it.home = (hx, hy)


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
    along = {key: float(_snap(v)) for key, v in _rungs(piece, na, centre, axis).items()}  # on the grid
    length = float(_snap(max(along.values()) - min(along.values()) + DNA_RISE * SCALE))  # half a step beyond each end
    a2 = frame.direction(centre, axis)
    angle = 0.0 if flat or np.linalg.norm(a2) < 1e-6 else _half_turn(_snap(np.arctan2(a2[1], a2[0])))
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
    items: list[_Item],
    frame,
    mirror: bool = True,
    shift: dict[str, int] | None = None,
    honest: tuple[set, dict[str, str]] | None = None,
    own: set[int] | None = None,
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
    `shift` (chain -> offset, model.homolog_chains: 2MJW numbers chain D 95-188 where B is 1-94).
    `honest` holds the strand pairs that pair (_sheet_pairs of every sheet) and each accessory's partner. A copy
    takes its mate's places only when that draws nothing against its pairing that its own drawing did not
    (_drawn_wrong, worse): no strand side by side with one it does not pair with, and no accessory away from a
    partner it stood beside, unless the mate's drawing names that same pairing. Otherwise the take is undone and the
    copy keeps its own drawing, turned as its mate is (ruling C2b-3/R2: science outranks A15; random-100 C2b-S,
    9P1X focus: L's sheet took V's places the other way round, and O's barrel took E's for two of R's strands with
    its ghost between them). A one-to-one take that is undone falls through to the partial take below: the copy is
    turned 180 degrees if its mate points the other way, then takes the places of the elements it shares
    (_take_shared, reflected when they run the other way), which `worse` judges again; if that passes, the copy is
    drawn as its mate's in part. The id of each copy item that keeps its own drawing goes in `own` (tied in its own
    shape, as _regroup_copies' are: no Layout.tied_ids)."""
    shift = shift or {}
    paired, partner = honest if honest is not None else (None, {})

    def keep(it: _Item):  # what a take changes: each member's strand, angle and length, and the offsets
        return [(m.sse, m.angle, m.length) if isinstance(m, Placed) else None for m in it.members], list(it.offsets)

    def undo(it: _Item, kept) -> None:
        for m, was in zip(it.members, kept[0]):
            if was is not None:
                m.sse, m.angle, m.length = was
        it.offsets = kept[1]

    def wrong(it: _Item) -> set:
        return set() if paired is None else _drawn_wrong(it, paired, partner)

    def worse(it: _Item, before: set, mate: _Item, taken: dict[str, str]) -> bool:
        """Does the take draw anything against the pairing that the copy's own drawing did not? An accessory the
        take draws away from its partner, named in the notes, is allowed only as its mate's own drawing: the take
        gave it its mate's place and its partner its mate's partner's (`taken`: copy strand id -> its mate's), and
        the mate's drawing names that very pairing (4XG1 B:174-175, as A:174-175); never when it only moved with
        them (9P1X R:2-8, which E's sheet lacks, would leave R:52-61), nor where its mate stands level beside its
        partner (fix round 1)."""
        added = wrong(it) - before
        if not added:
            return False
        named = {(w[1], w[2]) for w in wrong(mate) if w[0] == "named"}
        return any(w[0] != "named" or (taken.get(w[1]), taken.get(w[2])) not in named for w in added)

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
        one_to_one = len(mate.members) == len(it.members) and sorted(pairs) == list(range(len(pairs)))
        if one_to_one:
            was, before = keep(it), wrong(it)
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
            taken = {
                m.sse.id: mate.members[q].sse.id
                for m, q in zip(it.members, pairs)
                if isinstance(m, Placed) and not m.ghost and not mate.members[q].ghost
            }
            if worse(it, before, mate, taken):  # against its pairing
                undo(it, was)
                one_to_one = False
                if own is not None:
                    own.add(id(it))
        if not one_to_one:  # DSSP split, merged or missed an element in one copy: match the orientation
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
                was, before = keep(it), wrong(it)
                _take_shared(it, mate, shared)
                taken = {it.members[p].sse.id: mate.members[q].sse.id for p, q in shared if not it.members[p].ghost}
                if worse(it, before, mate, taken):
                    undo(it, was)
                    if own is not None:
                        own.add(id(it))
                elif own is not None:
                    own.discard(id(it))  # drawn as its mate's, in part
        ties[i] = (best, k)
    return ties


def _regroup_copies(
    items: list[_Item], ties: dict[int, tuple[int, int]], frame, shift: dict[str, int], keep: set[int]
) -> tuple[dict[int, tuple[int, int]], set[tuple[str, str]], set[int]]:
    """A copy's helices and loose elements grouped into items as its protomer-0 mate's are, where _symmetry_ties
    could not tie them item to item (random-100 C2b-3 A15: 6JLN e's helices, bundled in E's copy with R, a chain
    with no copy; 7QUP 13D's; 9P1X's). Elements are paired by model.pair_places (the same kind at the same place,
    residue numbers taken to the first copy's by `shift`), in the order the figure holds them, as
    figcheck.copies_alike pairs them, chain to chain of one class (one chain of the class in each protomer).
    - A copy item is clean when it is tied to the item holding its paired elements' mates, holds every element
      paired with that item's, and each sits at its mate's offset with its shape. A copy item holding paired
      elements that is not clean is dissolved (with the other holders of the same mate items), and its paired
      elements regrouped, one item per mate item, at the mate's offsets with its shape, tied to it; its unpaired
      elements join that group by their mean step, slid clear of its members. A helix piece let go (tied, with no
      pair, to a mate item another copy item pairs with) joins its own helix end to end (_on_end).
    - A sheet is never dissolved, split, merged or reshaped, nor are its ghosts dropped: a copy's strands are placed
      by that copy's own pairing (science S1-S3, S6; fix round 1, ruling R2: 6JLN b's sheet had been split as B's,
      drawing b:353-356 beside b:398-402 with no bridge). Nor is a βα row (`keep`, by id). A copy item no tie holds
      but whose elements pair (a sheet DSSP grouped otherwise: 7QUP 7D's) keeps its own shape and is tied to the
      mate item holding most of its pairs' mates, its offsets moved by their mean step: it stands with its copy.
    - Consecutive helix pieces whose mates are a kink (at most features.KINK_GAP apart) while they are not are
      returned as kink links (Layout.kink_links: their loop is drawn as their mates', 6JLN e:10-13>e:17-39).
    Items with elements of several protomers and the copies of a mirror are left as they are. Returns the ties, the
    kink links and the ids of the items tied in their own shape (kept out of Layout.tied_ids: not drawn as their
    mates'); `items` is changed in place (dissolved items dropped, new ones appended)."""
    sector, classes = getattr(frame, "sector", None) or {}, getattr(frame, "classes", None) or {}
    if not sector or not classes:
        return ties, set(), set()

    def plain(it: _Item, ks: set) -> bool:  # every element (ghosts aside) in protomers ks (None: outside every one)
        return (
            all(isinstance(m, Placed) for m in it.members)
            and {sector.get(m.sse.chain) for m in it.members if not m.ghost} <= ks
        )

    def fixed(it: _Item) -> bool:  # a sheet (its strands placed by this copy's own pairing, ghosts and barrel
        # closure too: C2b-3 fix round 1, ruling R2) or a βα row (`keep`): never dissolved
        return id(it) in keep or any(m.ghost or m.sse.kind == "E" for m in it.members)

    def at(members: list, m) -> int:  # the index of member m itself (equal members are not it)
        return next(q for q, x in enumerate(members) if x is m)

    owner = {id(m): i for i, it in enumerate(items) for m in it.members}
    chains: dict[tuple[int, int], list[str]] = {}  # (protomer, class) -> its chains
    for c in sorted(sector):
        if c in classes:
            chains.setdefault((sector[c], classes[c]), []).append(c)
    roots = [
        i
        for i, it in enumerate(items)
        if plain(it, {0, None}) and any(sector.get(m.sse.chain) == 0 for m in it.members)
    ]
    root_set = set(roots)
    ties = dict(ties)
    drop: set[int] = set()
    new: list[tuple[_Item, int, int]] = []  # (item, mate, protomer)
    left: list[_Item] = []  # unpaired elements of a dissolved item with no pair: untied, where they were
    kinks: set[tuple[str, str]] = set()
    own: set[int] = set()  # ids of the copy items tied in their own shape (not drawn as their mates': no tied_ids)

    def step_to(it: _Item, mine: _Item, ps) -> np.ndarray:  # mean move of it's members ps onto their mates in mine
        return np.mean(
            [np.subtract(mine.offsets[at(mine.members, mate[id(it.members[p])])], it.offsets[p]) for p in ps], axis=0
        )

    for k in sorted({v for v in sector.values() if v > 0}):
        copies = {i for i, it in enumerate(items) if plain(it, {k})}
        mate: dict[int, Placed] = {}  # id of a copy element -> its protomer-0 element
        kinked: list[tuple[Placed, Placed]] = []  # consecutive copy helices whose mates are a kink
        for (p, cls), mine in sorted(chains.items()):
            if p != k or len(mine) != 1 or len(chains.get((0, cls), [])) != 1:
                continue
            c0, ck = chains[(0, cls)][0], mine[0]

            def elements(chain: str, pool: set[int]) -> list[Placed]:  # in the order the figure holds them, as
                # figcheck._copy_pairs pairs them (lay.placed: item by item), which breaks ties between overlaps
                return [m for i in sorted(pool) for m in items[i].members if m.sse.chain == chain and not m.ghost]

            xs, ys = elements(c0, root_set), elements(ck, copies)

            def place(m: Placed) -> tuple[str, int, int]:
                off = shift.get(m.sse.chain, 0)
                return m.sse.kind, m.sse.first.seq - off, m.sse.last.seq - off

            for a, b in pair_places([place(m) for m in xs], [place(m) for m in ys]).items():
                mate[id(ys[b])] = xs[a]
            helices = sorted((m for m in ys if m.sse.kind == "H" and id(m) in mate), key=lambda m: m.sse.start)
            for a, b in zip(helices, helices[1:]):  # regrouped as their mates, a kink of the mate's is a kink here
                ma, mb = mate[id(a)].sse, mate[id(b)].sse
                if 0 <= mb.start - ma.end - 1 <= KINK_GAP and b.sse.start - a.sse.end - 1 > KINK_GAP:
                    kinked.append((a, b))

        def at_mates(it: _Item, mine: _Item, mate: dict[int, Placed] = mate) -> bool:  # each paired element at its
            # mate's offset, with its shape
            return all(
                it.offsets[p] == mine.offsets[q] and (m.angle, m.length) == (x.angle, x.length)
                for p, m in enumerate(it.members)
                if id(m) in mate
                for q, x in enumerate(mine.members)
                if x is mate[id(m)]
            )

        held: dict[int, set[int]] = {}  # mate item -> copy items holding its elements' pairs
        for i in sorted(copies):
            for m in items[i].members:
                if id(m) in mate:
                    held.setdefault(owner[id(mate[id(m)])], set()).add(i)
        dirty = set()
        for i in sorted(copies):
            got = {owner[id(mate[id(m)])] for m in items[i].members if id(m) in mate}
            if not got:  # tied to an item whose elements others pair with, it took the place of one of them (9P1X
                # U:227-229 on U:231-237, both under A:226-236)
                if i in ties and ties[i][0] in held and ties[i][1] == k:
                    dirty.add(i)
                continue
            (r,) = got if len(got) == 1 else (None,)
            if r is None or held[r] != {i} or ties.get(i) != (r, k) or not at_mates(items[i], items[r]):
                dirty.add(i)
        grown = True
        while grown:  # a mate item's other holders are regrouped too, and theirs in turn
            more = {i for held_by in held.values() if held_by & dirty for i in held_by}
            grown, dirty = bool(more - dirty), dirty | more
        locked = {i for i in dirty if fixed(items[i])}  # a sheet or a βα row is never dissolved: what shares a
        while locked:  # mate with one stays as it is, and so on
            more = {i for held_by in held.values() if held_by & locked for i in held_by} - locked
            dirty, locked = dirty - more, locked | more if more else set()
        dirty = {i for i in dirty if not fixed(items[i])}
        for i in sorted(copies - dirty - set(ties)):  # a copy item no tie holds (a sheet its mate's groups otherwise:
            # 7QUP 7D's) keeps its own shape, tied to the mate item holding most of its pairs' mates: it stands where
            # its mates' copy stands, not where it was projected (42 page units off), and moves with its copy
            ps = [p for p, m in enumerate(items[i].members) if id(m) in mate]
            if ps:
                to = [owner[id(mate[id(items[i].members[p])])] for p in ps]
                most = min(set(to), key=lambda r: (-to.count(r), r))
                step = step_to(items[i], items[most], [p for p, r in zip(ps, to) if r == most])
                items[i].offsets = [tuple(float(v) for v in _snap(np.add(o, step))) for o in items[i].offsets]
                ties[i] = (most, k)
                own.add(id(items[i]))
        if not dirty:
            continue
        groups: dict[int, list[tuple[Placed, tuple[float, float]]]] = {}
        extras: dict[int, list[tuple[Placed, tuple[float, float]]]] = {}
        for i in sorted(dirty):
            it = items[i]
            to = {p: owner[id(mate[id(m)])] for p, m in enumerate(it.members) if id(m) in mate}
            for p, r in to.items():
                q = at(items[r].members, mate[id(it.members[p])])
                it.members[p].angle, it.members[p].length = items[r].members[q].angle, items[r].members[q].length
                groups.setdefault(r, []).append((it.members[p], items[r].offsets[q]))
            rest = [p for p, m in enumerate(it.members) if p not in to and not m.ghost]
            if not rest:
                continue
            if to:  # what the mate lacks joins the group that took most of its item's elements, by their mean step
                most = min(set(to.values()), key=lambda r: (-sum(v == r for v in to.values()), r))
                step = step_to(it, items[most], [p for p, r in to.items() if r == most])
                extras.setdefault(most, []).extend(
                    (it.members[p], np.add(it.offsets[p], step), it, it.offsets[p]) for p in rest
                )
            elif not _on_end(items, ties, held, dirty, fixed, i, it, rest):  # a helix piece let go joins its helix
                mine_m = [it.members[p] for p in rest]
                left.append(
                    _Item(
                        mine_m, [it.offsets[p] for p in rest], it.cx, it.cy, it.home, *_anchor([m.sse for m in mine_m])
                    )
                )
        drop |= dirty
        regrouped = {id(m) for got in groups.values() for m, _ in got}
        kinks |= {(a.sse.id, b.sse.id) for a, b in kinked if id(a) in regrouped and id(b) in regrouped}
        for r, got in sorted(groups.items()):
            mates = items[r].members
            got.sort(key=lambda x: (at(mates, mate[id(x[0])]), x[0].sse.start))
            refused: dict[int, list] = {}
            for m, old, at_old in _into_group(got, extras.get(r, [])):
                refused.setdefault(id(old), [old]).append((m, at_old))
            for old, *those in refused.values():
                mine_m = [m for m, _ in those]
                left.append(
                    _Item(mine_m, [o for _, o in those], old.cx, old.cy, old.home, *_anchor([m.sse for m in mine_m]))
                )
            mine_m = [m for m, _ in got]
            item = _Item(
                mine_m, [o for _, o in got], items[r].cx, items[r].cy, items[r].home, *_anchor([m.sse for m in mine_m])
            )
            new.append((item, r, k))
    if not drop:
        return ties, kinks, own
    kept = [i for i in range(len(items)) if i not in drop]
    index = {old: n for n, old in enumerate(kept)}
    out = {index[i]: (index[j], k) for i, (j, k) in sorted(ties.items()) if i in index}
    items[:] = [items[i] for i in kept] + [it for it, _, _ in new] + left
    for n, (_, r, k) in enumerate(new):
        out[len(kept) + n] = (index[r], k)
    return dict(sorted(out.items())), kinks, own


def _boxes_meet(a, b) -> bool:
    """Do two boxes (x0, y0, x1, y1), MARGIN round each, overlap?"""
    return a[0] - MARGIN < b[2] and b[0] - MARGIN < a[2] and a[1] - MARGIN < b[3] and b[1] - MARGIN < a[3]


_JOIN_TRIES = 40  # _SLIDE_STEP moves an unpaired member joining its regrouped copy may take to clear its members


def _into_group(got: list, extra: list) -> list:
    """_regroup_copies: the unpaired members of a dissolved copy item (`extra`: (member, offset by the mean step,
    its old item, its old offset)) join the regrouped item `got` ((member, offset) pairs, appended to), each slid
    sideways away from the group's middle, _SLIDE_STEP at a time on snapped offsets, to the first spot clear of every
    member (at most _JOIN_TRIES steps). Those with no clear spot are returned as (member, old item, old offset), to
    stay untied where they were: never drawn over another element of the item (fix round 2)."""
    refused = []
    for m, o, old, at_old in extra:
        o = _snap(np.asarray(o, float))
        way = 1.0 if o[0] >= float(_snap(np.mean([q[0] for _, q in got] or [0.0]))) else -1.0
        for _ in range(_JOIN_TRIES):
            if _clear_of(m, o, got):
                got.append((m, tuple(float(v) for v in o)))
                break
            o = _snap(o + (way * _SLIDE_STEP, 0.0))
        else:
            refused.append((m, old, at_old))
    return refused


def _clear_of(m: Placed, o, placed: list) -> bool:
    """Is member m, at offset o of an item, clear (by MARGIN) of each (member, offset) of `placed`, by footprint?"""
    box = np.add(_footprint(m), (o[0], o[1], o[0], o[1]))
    return not any(_boxes_meet(box, np.add(_footprint(x), (q[0], q[1], q[0], q[1]))) for x, q in placed)


def _on_end(items, ties, held, dirty, fixed, i: int, it: _Item, rest: list[int]) -> bool:
    """A helix piece let go by _regroup_copies (tied to a mate item whose elements another copy item pairs with:
    7QUP 9C:240-243, cut from 9C:224-238 where 1C:224-243 is whole) joins that copy item when it is no sheet: end to
    end with the helix of its own chain next to it in sequence, along its axis, _PIECE_GAP apart (as a kinked
    helix's pieces stand), so it stays with its copy (it drifted 62 page units off as an item of its own). A spot
    another element of the item takes is passed by: the piece moves on along the axis, _SLIDE_STEP at a time, to the
    first clear spot (at most _JOIN_TRIES steps); with none, False and nothing joins (it stays untied)."""
    if i not in ties or len(held.get(ties[i][0], ())) != 1:
        return False
    (h,) = held[ties[i][0]]
    host = items[h]
    pieces = [it.members[p] for p in rest]
    if h in dirty or fixed(host) or any(m.sse.kind not in ("H", "G") for m in pieces):
        return False
    placed = list(zip(host.members, host.offsets))
    joins = []
    for m in sorted(pieces, key=lambda m: m.sse.start):
        near = [x for x in host.members if x.sse.chain == m.sse.chain and x.sse.kind in ("H", "G")]
        if not near:
            return False
        x = min(near, key=lambda x: (min(abs(m.sse.start - x.sse.end), abs(x.sse.start - m.sse.end)), x.sse.start))
        q = next(n for n, y in enumerate(host.members) if y is x)
        way = 1.0 if m.sse.start > x.sse.end else -1.0
        shape = Placed(m.sse, 0.0, 0.0, m.length, x.angle, m.width)  # the piece as it would be drawn
        for n in range(_JOIN_TRIES):
            reach = way * (x.length / 2 + m.length / 2 + _PIECE_GAP + n * _SLIDE_STEP)
            o = _snap(np.add(host.offsets[q], reach * x.direction))
            if _clear_of(shape, o, placed):
                break
        else:
            return False
        placed.append((shape, o))
        joins.append((m, x.angle, tuple(float(v) for v in o)))
    for m, angle, o in joins:
        m.angle = angle
        host.members.append(m)
        host.offsets.append(o)
    return True


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


def _drawn_wrong(it: _Item, pairs: set, partner: dict[str, str]) -> set:
    """What an item's strands, at its offsets, draw against their pairing (figcheck.science_checks' reading of S1 and
    S6, inside one item): two strands side by side (adjacent columns, across the page, spans overlapping) that do not
    pair (`pairs`: frozensets of strand ids), two overlapping in one column, and an accessory (`partner`: its id ->
    its partner's) not level beside its partner or the partner's ghost by columns: "away" when it is within a strand's
    width of it and level, so that _unshown counts it beside and the notes do not name the pairing (9P1X O's ghost
    wedged between R:29-35 and R:52-61), else "named" (the notes name it: S6's ruled case). Columns by x rounded to
    1e-3."""
    spots = []
    for m, (ox, oy) in zip(it.members, it.offsets):
        if isinstance(m, Placed) and m.sse.kind == "E" and not m.over:
            half = abs(np.sin(m.angle)) * m.length / 2
            spots.append((m.sse.id, round(float(ox), 3), (oy - half, oy + half), m.ghost, m.width))
    xs = sorted({x for _, x, _, _, _ in spots})
    col = {x: n for n, x in enumerate(xs)}

    def meet(a, b) -> bool:
        return min(a[1], b[1]) - max(a[0], b[0]) > _TOUCH

    out = set()
    for a, xa, sa, ga, wa in spots:
        for b, xb, sb, _, _ in spots:
            if a < b and meet(sa, sb):
                d = abs(col[xa] - col[xb])
                if d == 0 or (d == 1 and frozenset((a, b)) not in pairs):
                    out.add(("column" if d == 0 else "side", a, b))
        if ga or a not in partner:
            continue
        at = [(x, sp, w) for i, x, sp, _, w in spots if i == partner[a]]
        level = any(abs(col[x] - col[xa]) == 1 and meet(sa, sp) for x, sp, _ in at)
        counted = any(_TOUCH < abs(x - xa) <= max(w, wa) + _TOUCH and meet(sa, sp) for x, sp, w in at)  # _unshown
        if not level and at:
            out.add(("away" if counted else "named", a, partner[a]))
    return out


def _take_shared(it: _Item, mate: _Item, shared: list[tuple[int, int]]) -> None:
    """A copy's elements `shared` (member, its mate's member) take the mate's shape and place; the copy's other
    elements move with them, by their mean step, so each keeps its place beside the elements it is drawn with. When
    the mate's places run the other way across the page (the shared elements' x in reverse order: 6Y7S's FimA
    sheets, 9P1X L's against V's), the others are turned with them, mirrored across the page (x and angle), so a
    strand the mate lacks still stands beside the strand it pairs with, on the same side (random-100 C2b-S: B:122-123
    stood beside B:29-38, 0 bridges)."""
    taken = {p for p, _ in shared}
    crowded = {  # the copy's other elements already on another of its elements before the take: not the take's doing
        p for p, m in enumerate(it.members) if p not in taken and _on_member(m, it.offsets[p], _others(it, p))
    }
    own = np.array([it.offsets[p][0] for p, _ in shared])
    theirs = np.array([mate.offsets[q][0] for _, q in shared])
    turn = float(_snap(np.dot(own - own.mean(), theirs - theirs.mean()))) < 0  # on snapped offsets
    sign = -1.0 if turn else 1.0
    step = np.mean(
        [np.subtract(mate.offsets[q], (sign * it.offsets[p][0], it.offsets[p][1])) for p, q in shared], axis=0
    )
    offsets = [tuple(float(v) for v in _snap(np.add((sign * ox, oy), step))) for ox, oy in it.offsets]
    if turn:
        for p, m in enumerate(it.members):
            if p not in taken and isinstance(m, Placed):
                m.angle = _mirrored(m.angle)
    for p, q in shared:
        it.members[p].angle, it.members[p].length = mate.members[q].angle, mate.members[q].length
        offsets[p] = mate.offsets[q]
    it.offsets = offsets
    for p, m in enumerate(it.members):  # an element the mate lacks, landed on the shared ones: slid off them
        if p in taken or p in crowded or not _on_member(m, it.offsets[p], _others(it, p)):
            continue
        it.offsets[p] = _slid_off(m, it.offsets[p], _others(it, p))


def _others(it: _Item, p: int) -> list:
    """The (member, offset) pairs of item `it` but its member p."""
    return [(x, o) for q, (x, o) in enumerate(zip(it.members, it.offsets)) if q != p]


def _on_member(m, o, placed: list) -> bool:
    """Does member m, at offset o of an item, overlap any (member, offset) of `placed` by footprint (no margin)?"""
    a = np.add(_footprint(m), (o[0], o[1], o[0], o[1]))
    for x, q in placed:
        b = np.add(_footprint(x), (q[0], q[1], q[0], q[1]))
        if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]:
            return True
    return False


def _slid_off(m, o, placed: list) -> tuple[float, float]:
    """_take_shared: the nearest spot to offset o, _SLIDE_STEP at a time up, down, or sideways away from the item's
    middle (in that order at each distance), where member m stands clear (by MARGIN) of every (member, offset) of
    `placed` (random-100 held-out 2NAD: chain B's helix B:267-270, in its βα row, where chain A's row has none, took
    the place between two strands that its mate's row draws side by side, so it lay on both and its two loops had no
    route). At most _JOIN_TRIES steps; none clear: o as it was."""
    o = np.asarray(o, float)
    mid = float(_snap(np.mean([q[0] for _, q in placed] or [0.0])))
    way = 1.0 if o[0] >= mid else -1.0
    for n in range(1, _JOIN_TRIES + 1):
        for d in ((0.0, 1.0), (0.0, -1.0), (way, 0.0)):
            at = _snap(o + np.multiply(d, n * _SLIDE_STEP))
            if _clear_of(m, at, placed):
                return tuple(float(v) for v in at)
    return tuple(float(v) for v in o)


def _unwrap(items: list[_Item], ties: dict[int, tuple[int, int]], frame, ghosts: bool = True) -> dict[int, int]:
    """Elements lent across the cut of a ring unrolled into a row. A copy's element drawn in its neighbour's item
    (a strand that pairs in the next protomer's sheet) sits on the same side of its own protomer in every copy but
    one: the copy whose neighbour lies across the cut, at the row's other end. There the element leaves the item
    for an item of its own, tied to the protomer-0 mate's place the same number of protomers along as in the other
    copies (one past the row's end), so every copy is drawn alike; a strand leaves a ghost (dashed, as a barrel's)
    where it pairs, and takes ghosts of the strands it pairs with (_partners) to stand beside; with `ghosts` False
    (Style.copy_ghosts off) neither: the strand stands alone in its own item past the row's end, and no dashed copy
    of another copy's strand is drawn (a barrel's closing ghost is not made here and stays). New items are
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
            partners = _partners(it, ps) if ghosts else []
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
                it.members[p] = ghost if m.sse.kind == "E" and ghosts else None
            keep = [p for p, m in enumerate(it.members) if m is not None]
            it.members, it.offsets = [it.members[p] for p in keep], [it.offsets[p] for p in keep]
    return out


def _partners(it: _Item, ps: list[int]) -> list[Placed]:
    """The strands of item `it` beside its strands `ps` in the sheet, not among them: those they pair with, whose
    ghosts stand beside them where they leave (_unwrap): the strands in the next column, at most one pitch across,
    whose spans meet its own. A strand in its own column stands end to end with it (one column means collinear, no
    pairing: an accessory beside its partner shares a neighbour's column), so it lends no ghost (fix round 4: 8UTF's
    c:440-442 took the ghost of a:340-343, which paired with nothing it was drawn with)."""

    def strand(q):
        m = it.members[q]
        return isinstance(m, Placed) and not m.ghost and m.sse.kind == "E"

    out = []
    for p in ps:
        if not strand(p):
            continue
        mp, (x, y) = it.members[p], it.offsets[p]
        found = [
            q
            for q in range(len(it.members))
            if q != p
            and strand(q)
            and 1e-6 < abs(it.offsets[q][0] - x) <= max(mp.width, it.members[q].width) + 1e-6
            and abs(it.offsets[q][1] - y) < (mp.length + it.members[q].length) / 2
        ]
        out += [it.members[q] for q in found if q not in ps and it.members[q] not in out]
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


def _both_ways(build):
    """build_layout over a symmetry tie: a sheet shared by tied copies (_separate's `shared`) that runs straight
    across the page may point either way (_sheet_order draws strands[0] up), and one way the tied copies can need
    more room than _SEP_SYM_ROOM allows (2ILU: 1.07). When the copies are not drawn alike and such sheets exist, the
    layout is made again with them turned down (_TIE_TURN), and that one is kept if it draws the copies alike; else
    (both or neither alike) the first, strands[0] up. Only those layouts are made twice (CI ruling, plan 12). Copies
    drawn alike only in _separate's wider room (_tied_again: report "wide") count as not alike here, but beat copies
    not drawn alike at all (random-100 C2b-3)."""
    import functools

    @functools.wraps(build)
    def run(*args, **kwargs):
        def made(down: frozenset[str]):
            report: dict = {}
            token = _TIE_TURN.set((down, report))
            try:
                return build(*args, **kwargs), report
            finally:
                _TIE_TURN.reset(token)

        lay, first = made(frozenset())
        if first.get("alike", True) or not first.get("turnable"):
            return lay
        other, second = made(first["turnable"])
        return other if second.get("alike") or (second.get("wide") and not first.get("wide")) else lay

    return run


AXIS_SYMBOL_R = 0.45  # page units: a rosette's n-fold symbol at the page centre (render), corner to centre
_ROSETTE_STEP = 0.25  # page units: the rosette's radial push grows by this much until neighbouring copies clear
_ROSETTE_GAP = 2 * MARGIN  # clear space kept between two copies' elements (with their names): room for loops
_ROSETTE_SAMPLE = 0.4  # page units between the points an element is sampled at along its axis, for that clearance
_ROSETTE_MAX_PUSH = 200.0  # page units: the push never grows past this (a copy farther out clears by then anyway)


def _ref_chain(ref: str) -> str | None:
    """The chain a reference (resolve) names: res:A:..., A:label; None for one naming every chain."""
    ref = ref.strip()
    if ref.startswith("res:"):
        body = ref[4:]
        return body.split(":")[0] if body.count(":") == 1 else None
    if ref.startswith("#") or ":" not in ref:
        return None
    return ref.split(":", 1)[0]


def _copy_matches(zero: list[SSE], sses: list[SSE], frame, shift: dict[str, int]) -> list[dict[str, SSE]]:
    """For each protomer k, protomer 0's element id -> its copy in protomer k: of the same kind, in the chain of the
    same class (frame.classes), whose residues (numbered as the first copy's, `shift`: model.homolog_chains) overlap
    it most; one to one, the largest overlaps first (ties in sequence order)."""
    sector, classes = frame.sector, frame.classes or {}
    out = []
    for k in range(frame.n):
        if k == 0:
            out.append({s.id: s for s in zero})
            continue
        mine = [s for s in sses if sector.get(s.chain) == k]
        pairs = []
        for i, s in enumerate(zero):
            a0, b0 = s.first.seq - shift.get(s.chain, 0), s.last.seq - shift.get(s.chain, 0)
            for j, t in enumerate(mine):
                if t.kind != s.kind or classes.get(t.chain) != classes.get(s.chain):
                    continue
                a1, b1 = t.first.seq - shift.get(t.chain, 0), t.last.seq - shift.get(t.chain, 0)
                both = min(b0, b1) - max(a0, a1) + 1
                if both > 0:
                    pairs.append((-both, i, j))
        got: dict[str, SSE] = {}
        taken: set[int] = set()
        for _, i, j in sorted(pairs):
            if zero[i].id not in got and j not in taken:
                got[zero[i].id] = mine[j]
                taken.add(j)
        out.append(got)
    return out


def _rosette_samples(members: list[Placed]) -> tuple[np.ndarray, np.ndarray]:
    """Points along each member's axis (page units) and the radius each stands for: half the member's footprint
    across its axis (its name's room included, _footprint), so two copies clear when their points do."""
    pts, rad = [], []
    for m in members:
        x0, y0, x1, y1 = _footprint(m)
        across = max(-x0, x1) if abs(np.cos(m.angle)) < 0.5 else max(-y0, y1)
        along = max(-x0, x1) if abs(np.cos(m.angle)) >= 0.5 else max(-y0, y1)
        steps = max(2, int(np.ceil(2 * along / _ROSETTE_SAMPLE)) + 1)
        for t in np.linspace(-along, along, steps):
            pts.append([m.cx + m.direction[0] * t, m.cy + m.direction[1] * t])
            rad.append(across)
    return np.array(pts, float).reshape(-1, 2), np.array(rad, float)


def _rosette_push(members: list[Placed], angles: list[float], out: np.ndarray) -> float:
    """The smallest push (page units, a multiple of _ROSETTE_STEP) along `out` (each copy along its own turn of it)
    that keeps every copy's members _ROSETTE_GAP clear of every other's: protomer 0's against each turned copy, which
    by the symmetry covers every pair (the plan-21 spike's explode())."""
    from scipy.spatial.distance import cdist

    pts, rad = _rosette_samples(members)
    if not len(pts):
        return 0.0
    reach = rad[:, None] + rad[None, :] + _ROSETTE_GAP
    from .frame import turn2

    delta = 0.0
    while delta < _ROSETTE_MAX_PUSH:
        q = pts + delta * out
        # the turned copy snapped to the grid, as local() snaps it: a cos/sin ulp never decides the push
        if all((cdist(q, _snap(turn2(q, a))) - reach).min() >= -1e-9 for a in angles[1:]):
            break
        delta += _ROSETTE_STEP
    return float(delta)


def _box_on_body(box, e) -> bool:
    """Whether a page box meets an element's body (its length x width rectangle along its own axis, at any angle), by
    separating axes; touching is not meeting."""
    x0, y0, x1, y1 = box
    d = np.array([np.cos(e.angle), np.sin(e.angle)])
    n = np.array([-d[1], d[0]])
    c = np.array([e.cx, e.cy])
    body = np.array([c + a * e.length / 2 * d + b * e.width / 2 * n for a in (-1, 1) for b in (-1, 1)])
    square = np.array([[x0, y0], [x0, y1], [x1, y0], [x1, y1]])
    for axis in (np.array([1.0, 0.0]), np.array([0.0, 1.0]), d, n):
        a, b = body @ axis, square @ axis
        if min(a.max(), b.max()) - max(a.min(), b.min()) <= 1e-6:
            return False
    return True


def _rosette_label_sides(lay: Layout) -> None:
    """A turned copy's helix names (text upright), each side chosen after the turn: the side its mate's name is on,
    unless that box meets an element, a ghost, the axis symbol or a name chosen before it and the other side's does
    not. Protomer 0 keeps its own sides; copies are taken in order of their elements' ids."""
    k = glyph(lay)[0] / _GLYPH_W
    mates = lay.rosette_mates
    bodies = [*lay.placed.values(), *lay.ghosts]
    r = AXIS_SYMBOL_R + 0.15
    names = [
        _element_label_box(p, k) for p in lay.placed.values() if p.sse.kind in ("H", "G") and p.sse.id not in mates
    ]
    names.append((-r, -r, r, r))
    for i in sorted(mates):
        p = lay.placed.get(i)
        if p is None or p.sse.kind not in ("H", "G"):
            continue

        def meets(side: int, p: Placed = p) -> bool:  # its name's box on another (tilted) body, a name, the symbol
            p.label_side = side
            box = _element_label_box(p, k)
            return overlaps_any(box, names) or any(_box_on_body(box, q) for q in bodies if q is not p)

        first = p.label_side
        if meets(first) and not meets(-first):
            p.label_side = -first
        else:
            p.label_side = first
        p.label_side_default = p.label_side
        names.append(_element_label_box(p, k))


def _wrap_angle(a: float) -> float:
    """An angle in (-pi, pi], on the grid."""
    a = float(_snap((a + np.pi) % (2 * np.pi) - np.pi))
    return float(_HALF_TURN) if a <= -_HALF_TURN + 1e-9 else a + 0.0


def _rosette(sses: list[SSE], sheets: list[Sheet], frame, mode: str, kw: dict) -> Layout:
    """A Cn assembly seen down its axis (frame.RosetteFrame; Style.cyclic rosette, plan 21 Task 1): protomer 0 laid
    out alone in its own frame (frame.local(0): snapping, sheets, de-overlap as for any figure), then copy k drawn
    as that layout turned by frame.angles[k] (exactly k * 360 / n degrees) about the axis, which projects to page
    (0, 0), every copy pushed out from it by the smallest common push that clears its neighbours (_rosette_push).
    An element's angle is snapped in protomer 0's frame and turned with its copy, so each copy's elements snap in
    their own frame. Copy elements with no protomer-0 mate (DSSP found them in one copy only) are placed alone,
    seen in their copy's own frame (frame.local(k)) and moved out from the axis until clear. `sheets` are the
    sheets within protomers (cli makes them so: no sheet across two copies)."""
    from .frame import turn2

    style = kw.get("style") or Style()
    sector = frame.sector
    zero = [s for s in sses if sector.get(s.chain) == 0]
    ids0 = {s.id for s in zero}
    chains0 = {c for c, k in sector.items() if k == 0}
    sheets0 = [sh for sh in sheets if all(sector.get(s.chain) == 0 for s in sh.strands)]

    def mine(ref: str) -> bool:  # a reference to protomer 0 (or to every chain)
        return _ref_chain(ref) in (None, *chains0)

    sub = dict(kw)
    sub["relabel"] = {r: v for r, v in (kw.get("relabel") or {}).items() if mine(r)}
    sub["rename"] = {r: v for r, v in (kw.get("rename") or {}).items() if mine(r)}
    sub["swap"] = [(a, b) for a, b in kw.get("swap") or [] if mine(a) and mine(b)]
    sub["move"] = [(r, d) for r, d in kw.get("move") or [] if mine(r)]
    dropped = [  # the user's edits naming a copy's element: not applied (each copy is protomer 0 turned), noted
        *(r for r in kw.get("rename") or {} if not mine(r)),
        *(r for pair in kw.get("swap") or [] for r in pair if not mine(r)),
        *(r for r, _ in kw.get("move") or [] if not mine(r)),
    ]
    sub["bridges"] = [b for b in kw.get("bridges") or [] if b[0] in ids0 and b[2] in ids0]
    sub["domains"] = [(name, [i for i in ids if i in ids0]) for name, ids in kw.get("domains") or []]
    lay = build_layout(zero, sheets0, frame.local(0), mode, **sub)
    shift = {c: off for c, (_, off) in (kw.get("copies") or {}).items()}
    maps = _copy_matches(zero, sses, frame, shift)
    own = list(lay.placed.values())
    mid = np.mean([[p.cx, p.cy] for p in own], axis=0) if own else np.zeros(2)
    out = _snap(mid / np.linalg.norm(mid)) if np.linalg.norm(mid) > 1e-6 else np.array([0.0, 1.0])
    delta = _rosette_push(own + list(lay.ghosts), frame.angles, out)
    shove = out * delta

    def turned(m: Placed, t: SSE, a: float) -> Placed:
        x, y = (float(v) for v in _snap(turn2(np.array([m.cx, m.cy]) + shove, a)))
        q = Placed(t, x, y, m.length, _wrap_angle(m.angle + a), m.width, label=m.label, ghost=m.ghost, over=m.over)
        up = lambda angle: 1 if np.cos(angle) >= 0 else -1  # noqa: E731  the sign of the normal's page-up part
        q.label_side = q.label_side_default = m.label_side * up(m.angle) * up(q.angle)  # the mate's side of it
        return q

    placed: dict[str, Placed] = {}
    ghosts: list[Placed] = []
    for k, a in enumerate(frame.angles):
        for p in own:
            if p.sse.id in maps[k]:
                q = turned(p, maps[k][p.sse.id], a)
                placed[q.sse.id] = q
        ghosts += [turned(g, maps[k][g.sse.id], a) for g in lay.ghosts if g.sse.id in maps[k]]
    relabel = kw.get("relabel") or {}
    strays: list[Placed] = []
    for k, a in enumerate(frame.angles):  # copy elements no protomer-0 element maps to: placed alone, clear
        matched = {t.id for t in maps[k].values()}
        local = frame.local(k)
        for t in sses:
            if sector.get(t.chain) != k or t.id in matched:
                continue
            if t.kind == "E":
                it = _sheet_item(Sheet([t], [1], []), local, PITCH * style.strand_scale)
            else:
                width = HELIX_W * style.helix_scale * (ETA_W if t.kind == "G" else 1.0)
                # projected mode only (cli._no_rosette keeps a stack row unrolled): the style's helix angles
                it = _helix_item(t, local, style.helix_angle, width, TEN_RISE if t.kind == "G" else HELIX_RISE)
            m = it.members[0]
            m.cx, m.cy = it.cx + it.offsets[0][0], it.cy + it.offsets[0][1]
            q = turned(m, t, a)
            q.label = relabel.get(res_ref(t.chain, t.first, t.last))
            if q.label is None:  # the first label of its kind its chain does not use (a copy's lone η: no second η1)
                used = {r.label for r in [*placed.values(), *strays] if r.sse.chain == t.chain}
                free = (generic_label(t.kind, n) for n in range(len(used) + 1))
                q.label = next(label for label in free if label not in used)
            radial = np.array([q.cx, q.cy])
            radial = radial / np.linalg.norm(radial) if np.linalg.norm(radial) > 1e-6 else turn2(out, a)
            taken = [r.rect for r in [*placed.values(), *ghosts, *strays]]
            for _ in range(int(_ROSETTE_MAX_PUSH / _ROSETTE_STEP)):
                x0, y0, x1, y1 = q.rect
                if not overlaps_any((x0 - MARGIN, y0 - MARGIN, x1 + MARGIN, y1 + MARGIN), taken):
                    break
                q.cx, q.cy = (float(v) for v in _snap(np.array([q.cx, q.cy]) + radial * _ROSETTE_STEP))
            strays.append(q)
    order = {s.id: i for i, s in enumerate(sses)}
    every = sorted([*placed.values(), *strays], key=lambda p: order[p.sse.id])
    lay.placed = {p.sse.id: p for p in every}
    lay.ghosts = ghosts
    for p in every:
        p.ends = ""
    for chain, mine_ in group_by_chain(sses, keep=lay.placed).items():
        if chain not in (kw.get("partial") or ()):
            lay.placed[mine_[0].id].ends += "N"
            lay.placed[mine_[-1].id].ends += "C"

    def each(ids) -> list[list[str]]:  # protomer 0's ids as each copy's (those with a copy there)
        return [[maps[k][i].id for i in ids if i in maps[k]] for k in range(frame.n)]

    sheets_, ghosts_ = [], []
    for ids, gids in zip(lay.sheets, lay.sheet_ghosts):
        for got, gone in zip(each(ids), each(gids)):
            if got:
                sheets_.append(got)
                ghosts_.append(gone)
    for q in strays:
        if q.sse.kind == "E":
            sheets_.append([q.sse.id])
            ghosts_.append([])
    lay.sheets, lay.sheet_ghosts = sheets_, ghosts_
    lay.ba_rows = [row for r in lay.ba_rows for row in each(r) if row]
    lay.on_sheet = {
        maps[k][h].id: [maps[k][s].id for s in strands if s in maps[k]]
        for h, strands in lay.on_sheet.items()
        for k in range(frame.n)
        if h in maps[k]
    }
    lay.accessory = {
        maps[k][a].id: maps[k][b].id
        for a, b in lay.accessory.items()
        for k in range(frame.n)
        if a in maps[k] and b in maps[k]
    }
    lay.unshown = sorted(
        (maps[k][a].id, maps[k][b].id) for a, b in lay.unshown for k in range(frame.n) if a in maps[k] and b in maps[k]
    )
    lay.kink_links = {
        (maps[k][a].id, maps[k][b].id)
        for a, b in lay.kink_links
        for k in range(frame.n)
        if a in maps[k] and b in maps[k]
    }
    lay.sheet_items = sheet_note_items(sheets, lay.unshown)
    lay.sheet_notes = [i["text"] for i in lay.sheet_items]
    lay.sheet_joins = sorted({tuple(sorted((a.id, b.id))) for sh in sheets for a, b in sh.joined})
    zero_of = {c: k for c, k in (frame.classes or {}).items() if sector.get(c) == 0}
    first = {i: c for c, i in zero_of.items()}
    lay.copy_ties = {c: (first[frame.classes[c]], False) for c, k in sector.items() if k and frame.classes[c] in first}
    lay.tied_ids = {t.id for k in range(1, frame.n) for i, t in maps[k].items()} | {
        i for k in range(1, frame.n) for i in maps[k]
    }
    lay.mirror_pairs = []
    chain_of = {k: p.sse.chain for k, p in lay.placed.items()}
    lay.domains = [
        (name, [k for k in ids if k in lay.placed and chain_of[k] not in lay.partial])
        for name, ids in kw.get("domains") or []
    ]
    lay.domain_lent = {}
    lay.frame = frame
    lay.rosette = frame.n
    lay.rosette_push = delta
    lay.rosette_strays = [q.sse.id for q in strays]
    lay.rosette_dropped = list(dict.fromkeys(dropped))
    lay.rosette_mates = {t.id: (i, k) for k in range(1, frame.n) for i, t in maps[k].items()}
    lay.stack_gaps = _stack_gaps(lay)
    _rosette_label_sides(lay)
    _choose_label_sides(lay, sses, style, [])
    _place_domain_names(lay, sses, style.font_scale)
    rects = np.array([p.rect for p in [*lay.placed.values(), *lay.ghosts]]).reshape(-1, 4)
    rects = np.vstack([rects, *[[domain_panel(lay, name)] for name, _ in lay.domains if _]]) if lay.domains else rects
    rects = np.vstack([rects, [b for _, b in label_boxes(lay, sses)] or np.empty((0, 4))])
    lay.bounds = (
        float(rects[:, 0].min()),
        float(rects[:, 1].min()),
        float(rects[:, 2].max()),
        float(rects[:, 3].max()),
    )
    return lay


_UNIT_AGREES = 1e-9  # page units: a Like.unit is read only when it gives its target's place to this


def _like_homes(items: list[_Item], like: Like) -> dict[int, tuple[float, float]]:
    """Each item holding an element with a place in the reference (Like.target) starts at the place they give it,
    which becomes its home: the mean over those members of their mate's place less their own offset. A member's place
    is read as its mate's unit centre plus the difference of the two offsets (Like.unit) when that unit agrees with
    its target (within _UNIT_AGREES), else as its target on the grid (_snap: a place computed with sin and cos, which
    differ in the last bit between CPUs, never reaches _attract unsnapped) less its offset; the mean is taken from the
    first one's, so an item laid out as its mate's unit (a structure drawn like itself) lands exactly on that unit's
    centre. Returns item index -> home (_attract's `held`)."""
    out: dict[int, tuple[float, float]] = {}
    for k, it in enumerate(items):
        got = []
        for m, (ox, oy) in zip(it.members, it.offsets):
            if not isinstance(m, Placed) or m.ghost or m.sse.id not in like.target:
                continue
            tx, ty = like.target[m.sse.id]
            unit = like.unit.get(m.sse.id)
            agrees = unit is not None and max(abs(unit[0] + unit[2] - tx), abs(unit[1] + unit[3] - ty)) <= _UNIT_AGREES
            if agrees:
                ux, uy, ax, ay = unit
                got.append((ux + (ax - ox), uy + (ay - oy)))
            else:
                sx, sy = (float(v) for v in _snap((tx, ty)))
                got.append((sx - ox, sy - oy))
        if not got:
            continue
        x0, y0 = got[0]
        x = x0 + math.fsum(x - x0 for x, _ in got) / len(got)
        y = y0 + math.fsum(y - y0 for _, y in got) / len(got)
        it.cx, it.cy = x, y
        it.home = (x, y)
        out[k] = (x, y)
    return out


@_text_scaled
@_both_ways
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
    alike_first: bool = True,
    let_go: bool = True,
    mirror_only: bool = False,
    like: Like | None = None,
) -> Layout | None:
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
    its length is shorter than the slab it crosses) is moved, the least way, until one end meets a face.
    `alike_first`: a copy's elements are grouped as its mate's (_regroup_copies) and the copies the symmetry ties may
    be drawn alike in up to _SEP_TIE_ROOM (_separate's second chance, _tied_again), random-100 C2b-3; cli turns it
    off in a membrane view, where the transmembrane helices' heights are judged instead (C17).
    `let_go`: of several copy items tied to one mate in one protomer, only the largest stays tied (_one_per_mate);
    Layout.let_go counts the others, and cli.make_layout lays out again without letting them go and keeps the more
    alike of the two (random-100 held-out A15). `mirror_only` (cli's dyad view of a two-fold assembly, plan 18): None
    at once where the frame's mirror line is dropped (a βα row across it), for cli to unroll the copies instead.
    `like` (like.Like, plan 13a; `frame` then frame.frame_like's, the reference's carried over): the elements are
    labelled as their mates in the reference (merged over `relabel`); in projected mode a matched item (one holding an
    element with a mate) starts at the reference's place for it, the mean over its matched members of their mates'
    places less their offsets, and is held there (_attract's `held`, Like.weight times _K_HOME) while the others keep
    their projected homes; _separate runs as always. Matched sheets take the reference's βα-row choices
    (Like.template) and way up (_sheet_order: their own strand order, science S1-S3, S6), matched helices its drawn
    angles and over-sheet choices (Like.over), one way only: a helix whose mate the reference draws beside its sheet
    is drawn beside, but one whose mate lies over its sheet is drawn over only when `over` (features.helices_on_sheets)
    finds it packed on a sheet here too (an over-sheet helix needs its own sheet). With like=None nothing of this
    runs."""
    if mode not in ("projected", "stack"):
        raise ValueError(f"unknown layout mode {mode!r}")
    if like is not None and (type(frame) is not Frame or frame.sector or frame.mirror_x is not None):
        raise ValueError(f"a layout drawn like a reference takes a plain frame (frame_like): {type(frame).__name__}")
    if isinstance(frame, RosetteFrame):  # a Cn assembly seen down its axis: protomer 0 laid out, the copies turned
        given = dict(attract=attract, contacts=contacts, nucleic=nucleic, style=style, bridges=bridges)
        given.update(domains=domains, partial=partial, rename=rename, swap=swap, move=move, keep_depth=keep_depth)
        given.update(na_cropped=na_cropped, template=template, copies=copies, relabel=relabel, names=names, over=over)
        given.update(level=level, band=band, alike_first=alike_first, let_go=let_go)
        return _rosette(sses, sheets, frame, mode, given)
    style = style or Style()
    angle_mode = "upright" if mode == "stack" else style.helix_angle
    template = style.fold_template if template is None else template
    if template not in CHOICES["fold_template"]:
        raise ValueError(f"fold_template must be one of {', '.join(CHOICES['fold_template'])}; got {template!r}")
    as_rows = mode == "projected" and template == "auto"  # a stack row draws elements in sequence instead
    rowed = template_choices(sheets, sses, copies or {}) if as_rows else [False] * len(sheets)
    seen_as = like.angle if like is not None and mode == "projected" else None  # the reference's drawn angles
    if like is not None and as_rows and like.template is not None:  # each sheet a βα row as the reference's mate
        if len(like.template) != len(sheets):
            raise ValueError(f"Like.template has {len(like.template)} choices for {len(sheets)} sheets")
        rowed = [bool(x) for x in like.template]
    if like is not None and mode == "projected":  # over its sheet only where the reference's mate is
        over = {h: k for h, k in (over or {}).items() if like.over.get(h, True)}
    templates = {
        k: _ba_template(
            ring := _open_at_ends(sh),
            _ba_units(ring, sses),
            frame,
            PITCH * style.strand_scale,
            HELIX_W * style.helix_scale,
            seen_as,
        )
        for k, (sh, row) in enumerate(zip(sheets, rowed))
        if row
    }
    in_rows = {m.sse.id for it in templates.values() for m in it.members if m.sse.kind == "H"}  # helices of βα rows
    mirror = getattr(frame, "mirror_x", None)  # a dyad in the page: protomer 1 drawn as protomer 0's mirror image
    rows_x = [(it.cx + it.box[0], it.cx + it.box[2]) for it in templates.values()]
    if mirror is not None and any(x0 < mirror < x1 for x0, x1 in rows_x):
        mirror = None  # a βα row across the mirror line would meet its own image: the copies go side by side instead
        if mirror_only:
            return None
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
    chain_ends: dict[str, str] = {}  # the chain ends drawn (as set on the members below), for stacked columns
    for chain, mine in group_by_chain([s for s in sses if s.kind in ("H", "G", "E")]).items():
        if chain not in (partial or ()):
            chain_ends[mine[0].id] = chain_ends.get(mine[0].id, "") + "N"
            chain_ends[mine[-1].id] = chain_ends.get(mine[-1].id, "") + "C"
    items = [
        templates[k]
        if k in templates
        else _sheet_item(
            sh,
            frame,
            PITCH * style.strand_scale,
            [by_id[h] for h, at in sorted(on_sheet.items(), key=lambda x: by_id[x[0]].start) if at == k],
            HELIX_W * style.helix_scale,
            chain_ends,
            seen_as,
        )
        for k, sh in enumerate(sheets)
    ]
    helix_w = HELIX_W * style.helix_scale
    items += [_bundle_item([by_id[k] for k in g], frame, angle_mode, helix_w, seen_as) for g in bundles]
    items += [_kink_item(g, frame, angle_mode, helix_w, seen_as) for g in kinked]
    lying = {m.sse.id for it in items for m in it.members if m.over}  # over a sheet: on_sheet less any crowded out
    grouped = bundled | in_rows | {s.id for g in kinked for s in g} | lying
    seen_as = seen_as or {}
    items += [
        _helix_item(s, frame, angle_mode, helix_w, angle=seen_as.get(s.id))
        for s in sses
        if s.kind == "H" and s.id not in grouped
    ]
    items += [
        _helix_item(s, frame, angle_mode, HELIX_W * ETA_W * style.helix_scale, TEN_RISE, angle=seen_as.get(s.id))
        for s in sses
        if s.kind == "G"
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
    held = _like_homes(items, like) if like is not None and mode == "projected" else {}
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
    kink_links: set[tuple[str, str]] = set()
    own: set[int] = set()
    gone = 0  # copy items _one_per_mate let go
    lends = False  # a strand lent across an unrolled ring's cut (_unwrap): Style.copy_ghosts is a choice
    apart: dict[str, int] = {}  # copy_ghosts off: a strand lent across the cut -> id() of the item it left
    if mode == "stack":
        rows = _stack(items, keep_depth)
        _centre_dna_over_binders(items, contacts)
    else:
        shift = {c: off for c, (_, off) in (copies or {}).items()}
        honest = (
            set().union(*map(_sheet_pairs, sheets)),
            {a.strand.id: a.parent.id for sh in sheets for a in sh.accessory},
        )
        kept_own: set[int] = set()  # copy items whose mates' places would misdraw their pairing: their own shape
        ties = _symmetry_ties(items, frame, mirror is not None, shift, honest, kept_own)
        regroup = alike_first and mirror is None and not pieces and walls is None  # (no stem or trace indexes items)
        # a copy split in two (_one_per_mate): the larger part tied, the other let go once _regroup_copies has run
        # (let go before it, a sheet is tied again beside its copy: 9GTP's whole file spread 35 % wider)
        alone = (
            {id(items[i]) for i in set(ties) - set(_one_per_mate(items, ties, helices=not regroup))}
            if let_go
            else set()
        )
        gone = len(alone)
        own = set(kept_own)
        if regroup:
            ties, kink_links, own = _regroup_copies(items, ties, frame, shift, {id(it) for it in templates.values()})
            own |= kept_own
        ties = {i: t for i, t in ties.items() if id(items[i]) not in alone}
        sector_of = getattr(frame, "sector", None) or {}
        left = _unwrap(items, ties, frame, style.copy_ghosts == "on") if mirror is None else {}
        lends = any(m.sse.kind == "E" for i in left for m in items[i].members if isinstance(m, Placed) and not m.ghost)
        if style.copy_ghosts == "off":  # each lent strand -> the item it left, whose sheet it pairs in
            apart = {m.sse.id: id(items[x]) for i, x in left.items() for m in items[i].members if m.sse.kind == "E"}
        if left and walls is not None:  # an element that left its item keeps off the traces it kept off
            walls = Walls(walls.boxes, np.vstack([walls.hard, walls.hard[list(left.values())]]))
        copy_ties = _copy_ties(items, ties, frame, mirror is not None)
        tied_ids = {
            m.sse.id for i, (j, _) in ties.items() for k in (i, j) if id(items[k]) not in own for m in items[k].members
        }
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
            if pieces and not ties and walls is None:
                _room_for_stems(items)
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
                held=held,
                held_k=_K_HOME * like.weight if like is not None else _K_LIKE,
                settled=like.settled if like is not None else None,
            )
        if kept is not None and kept.any():  # transmembrane elements at their projected height, reaching a face
            for k in np.flatnonzero(kept).tolist():
                shift = _to_face(items[k], level, band, heights[k]) if band is not None else 0.0
                y = float(_snap(heights[k] + shift))
                items[k].cy = y
                items[k].home = (items[k].home[0], y)
        shared = {  # an item of several protomers that no tie holds: a dimer's continuous sheet (1SLV)
            k
            for k, it in enumerate(items)
            if k not in ties
            and k not in {j for j, _ in ties.values()}
            and len({sector_of.get(m.sse.chain) for m in it.members if isinstance(m, Placed)} - {None}) > 1
        }
        state = _separate(  # "alike", "wide" (alike only by the second chance, in more room) or "alone"
            items, ties, joints, walls, mirror=mirror if ties else None, level=kept, shared=shared, again=alike_first
        )
        if domains and _part_domains(items, domains, ties, mirror if ties else None):
            # the domains, parted with their copies, keep them alike whatever room that takes
            state = _separate(
                items, ties, joints, walls, mirror=mirror if ties else None, room=np.inf, level=kept, again=alike_first
            )
        report = _TIE_TURN.get()[1]
        if report is not None:  # the shared sheets whose way up is a symmetry tie: _both_ways may turn them
            held = {m.sse.id for k in shared for m in items[k].members}
            report["turnable"] = frozenset(
                sh.strands[0].id
                for k, sh in enumerate(sheets)
                if k not in templates and sh.strands[0].id in held and _sheet_seen(sh, frame)[1] == 0
            )
            report["alike"], report["wide"] = state == "alike", state == "wide"
    for it in items:
        it.apply()
    dna = [m for it in items for m in it.members if isinstance(m, PlacedStem)]
    placed = [m for it in items for m in it.members if isinstance(m, Placed) and not m.ghost]
    ghosts = [m for it in items for m in it.members if isinstance(m, Placed) and m.ghost]
    _assign_labels(placed)
    real = {p.sse.id: p for p in placed}
    named = Layout(real, (0.0, 0.0, 0.0, 0.0), stems=dna, nucleic=nucleic if dna else None)
    named.ba_rows = _ba_rows(items)
    if like is not None and like.label:  # the reference's labels for its mates (and B's extras named after them)
        relabel = {**(relabel or {}), **like.label}
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
        (ids, [m.sse.id for m in it.members if isinstance(m, Placed) and m.ghost], id(it))
        for it in items
        if (ids := [m.sse.id for m in it.members if isinstance(m, Placed) and m.sse.kind == "E" and not m.ghost])
    ]
    block_of = {item: k for k, (_, _, item) in enumerate(blocks)}
    lay = Layout(
        real,
        (0.0, 0.0, 0.0, 0.0),
        ghosts=ghosts,
        stems=dna,
        nucleic=nucleic if dna or traces else None,
        sheets=[ids for ids, _, _ in blocks],
        sheet_ghosts=[ghost_ids for _, ghost_ids, _ in blocks],
    )
    lay.na_traces = traces
    lay.frame = frame
    lay.units = {
        m.sse.id: (it.cx, it.cy, ox, oy)
        for it in items
        for m, (ox, oy) in zip(it.members, it.offsets)
        if isinstance(m, Placed) and not m.ghost
    }
    lay.copy_ties = copy_ties
    lay.tied_ids = tied_ids
    lay.kink_links = kink_links
    lay.let_go = gone
    lay.lent_across_cut = lends
    lay.lent_apart = {i: block_of[item] for i, item in apart.items() if i in real and item in block_of}
    lay.mirror_pairs = [(a.sse.id, b.sse.id) for a, b in mirrored if isinstance(a, Placed) and a.sse.kind in ("H", "G")]
    lay.font_scale = style.font_scale
    lay.loop_stroke = style.loop_width * (RUN_EDGE_W if style.loop_color == "residue" else 1.0)
    lay.ligands = style.ligands
    lay.ba_rows = _ba_rows(items)  # after the user's swaps
    lay.stack_gaps = _stack_gaps(lay)
    lay.accessory = {a.strand.id: a.parent.id for sh in sheets for a in sh.accessory if a.strand.id in real}
    lay.unshown = _unshown(real, sheets, ghosts)
    lay.sheet_items = sheet_note_items(sheets, lay.unshown)
    lay.sheet_notes = [i["text"] for i in lay.sheet_items]
    lay.sheet_joins = sorted({tuple(sorted((a.id, b.id))) for sh in sheets for a, b in sh.joined})
    lay.on_sheet = {
        m.sse.id: [x.sse.id for x in it.members if x.sse.kind == "E" and not x.ghost]
        for it in items
        for m in it.members
        if isinstance(m, Placed) and m.over
    }
    if traces:  # a trace's nucleotides nudged off the elements that do not hold them (na_ends.off_unheld)
        from .na_ends import off_unheld

        for s in sorted(traces):
            traces[s] = off_unheld(lay, s, traces[s])
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


def _level_beside(x: Placed, p: Placed) -> bool:
    """x stands in a column next to p's, their spans overlapping (level with it)."""
    beside = 1e-6 < abs(x.cx - p.cx) <= max(x.width, p.width) + 1e-6
    (x0, x1), (p0, p1) = sorted((x.n_port[1], x.c_port[1])), sorted((p.n_port[1], p.c_port[1]))
    return beside and min(x1, p1) - max(x0, p0) > 1e-6


def extra_pairs(sheet: Sheet) -> list[tuple]:
    """The pairings of a sheet's accessory strands other than with their parent: with a strand in Accessory.others,
    a strand of the track or another accessory (a pair of two accessories once). Not those of a parent's child,
    which is that child's own (accessory, parent) pairing."""
    own = {frozenset((a.strand.id, a.parent.id)) for a in sheet.accessory}
    out, seen = [], set(own)
    for a in sheet.accessory:
        for o in a.others:
            key = frozenset((a.strand.id, o.id))
            if key not in seen:
                seen.add(key)
                out.append((a.strand, o))
    return out


def _unshown(placed: dict[str, Placed], sheets: list[Sheet], ghosts: list[Placed] = ()) -> list[tuple[str, str]]:
    """The (accessory, partner) pairings the figure does not draw side by side, as drawn (after copies take their
    mates' shapes): an accessory not in a column next to its partner's, or not level with it (spans apart). The
    partner's ghost (a barrel's closure strand) standing level beside it counts as the partner. Also every pairing
    of an accessory with a strand of Accessory.others (the track, or another accessory) not drawn so."""
    ghost_of: dict[str, list[Placed]] = {}
    for g in ghosts:
        ghost_of.setdefault(g.sse.id, []).append(g)
    out = []
    for sh in sheets:
        pairs = [(a.strand, a.parent) for a in sh.accessory] + extra_pairs(sh)
        for xs, ps in pairs:
            x, p = placed.get(xs.id), placed.get(ps.id)
            if x is None or p is None:
                continue
            if not any(_level_beside(x, q) for q in (p, *ghost_of.get(ps.id, ()))):
                out.append((xs.id, ps.id))
    return sorted(out)


def accessory_lines(sheet: Sheet, unshown: list[tuple[str, str]] = (), level: set[str] | None = None) -> list[str]:
    """One line per pairing a sheet draws off its run of paired strands or not at all, worded once for the notes
    (sheet_note_items) and `foldmap summary`: each accessory strand with its partner, "D:662-666 drawn beside its
    partner D:645-659", or, drawn away from it (`unshown`: not beside it and level, _unshown), "x:67-68 also pairs
    with x:79-86"; and each pairing of its track not drawn side by side (Sheet.chords: 9P1X's ring chords), or of an
    accessory with another strand of its sheet likewise, "O:29-36 also pairs with P:28-35". `level`, for `foldmap
    summary` (which reads a figure from outside): the accessories that figure draws level beside their partner; any
    other, drawn away or not drawn (focus), reads "A1:2-6 pairs with A1:9-23 (off the run)"."""
    away = set(unshown)

    def line(a) -> str:
        if level is not None:
            return (
                f"{a.strand.id} drawn beside its partner {a.parent.id}"
                if a.strand.id in level
                else f"{a.strand.id} pairs with {a.parent.id} (off the run)"
            )
        if (a.strand.id, a.parent.id) in away:
            return f"{a.strand.id} also pairs with {a.parent.id}"
        return f"{a.strand.id} drawn beside its partner {a.parent.id}"

    out = [line(a) for a in sheet.accessory]
    extras = [(a, b) for a, b in extra_pairs(sheet) if (a.id, b.id) in away]  # drawn side by side: nothing to say
    return out + [f"{a.id} also pairs with {b.id}" for a, b in [*sheet.chords, *extras]]


def sheet_note_items(sheets: list[Sheet], unshown: list[tuple[str, str]] = ()) -> list[dict]:
    """What the figure says of its sheets, as provenance items (kind ambiguous_pairing: text, count, elements): the
    pairings joining barrels drawn apart (Sheet.joined); then every pairing drawn off a sheet's run of paired strands
    or not drawn (accessory_lines), with how many accessories are named only, not drawn level with their partner
    (`unshown`, plan 12 ruling); then the pairings of a track not drawn side by side (Sheet.chords)."""
    out = []

    def item(text: str, count: int, pairs) -> dict:
        ids = list(dict.fromkeys(i for pair in pairs for i in pair))
        return {"kind": "ambiguous_pairing", "text": text, "count": count, "elements": ids}

    joined = sorted({tuple(sorted((a.id, b.id))) for sh in sheets for a, b in sh.joined})
    if joined:  # two barrels of one pairing network (sheets._other_rings), drawn each as its own
        joins = ", ".join(" with ".join(pair) for pair in joined)
        out.append(item(f"barrels drawn apart, each closed, that also pair: {joins}", len(joined), joined))
    lone, chords, lone_pairs, chord_pairs = [], [], [], []  # each sheet's accessories first, then its chords
    for sh in sheets:
        lines = accessory_lines(sh, unshown)
        lone += lines[: len(sh.accessory)]
        chords += lines[len(sh.accessory) :]
        lone_pairs += [(a.strand.id, a.parent.id) for a in sh.accessory]
        chord_pairs += [(a.id, b.id) for a, b in sh.chords]
        chord_pairs += [(a.id, b.id) for a, b in extra_pairs(sh) if (a.id, b.id) in set(unshown)]
    if lone:
        named = sum(pair in set(unshown) for pair in lone_pairs)
        text = f"{len(lone)} strand(s) off a sheet's run of paired strands: {'; '.join(lone)}"
        if named:  # drawn away from their partner: the note is all that pairs them
            text += f" ({named} not level with their partner, named only)"
        out.append(item(text, len(lone), lone_pairs))
    if chords:
        out.append(item(f"pairings not drawn side by side: {'; '.join(chords)}", len(chords), chord_pairs))
    return out


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
_NAME_OPEN = _NAME_ROUTE + 0.1  # page units a name keeps from a loop's way out of an element's end (_ways_out, 0.2
# either side of it) at the least: its routing keep-out (_NAME_ROUTE out) 0.1 off the way, with a route's clearance
# (route._INFLATE, 0.15), leaves the cells along the way's middle free, half a cell (route.RES) either side of it, so
# a route can leave the end (random-100 held-out 7NZ2's domains view: "D1 · B1" 0.08 from A1:1330-1337's way out,
# with no room to keep, closed it); a leader's keep-out keeps the same 0.1
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


def _ways_out(layout: Layout) -> list[tuple]:
    """The straight run of a loop out of each element's (and ghost's) end, _NAME_STUB long, 0.2 either side."""
    ways = []
    for q in [*layout.placed.values(), *layout.ghosts]:
        for port, out in ((q.n_port, q.exit_n), (q.c_port, q.exit_c)):
            (ax, ay), (bx, by) = port, np.add(port, out * _NAME_STUB)
            ways.append((min(ax, bx) - 0.2, min(ay, by) - 0.2, max(ax, bx) + 0.2, max(ay, by) + 0.2))
    return ways


@dataclass(frozen=True, eq=False)
class _Way:
    """A turned box (length along `direction`, width across it) centred at (cx, cy), as _box_meets_body reads one."""

    cx: float
    cy: float
    direction: np.ndarray
    normal: np.ndarray
    length: float
    width: float


def _stem_ways(layout: Layout) -> list[_Way]:
    """The straight way of a strand out of each end of each stem: _NAME_STUB along the stem's axis past its end
    face, as wide as the stem and 0.2 either side, as a protein loop's (_ways_out). A nucleic-acid run, a cropped
    gap or a free end leaves a stem's end along its axis (na_paths._face): a domain's name there is crossed by it
    (random-100 held-out 8ZTU's domains view: "D2 · H" 0.6 past a stem's end, on the cropped gap leaving it)."""
    out = []
    for d in layout.stems:
        for sense in (1, -1):
            c = np.array([d.cx, d.cy]) + d.direction * sense * (d.length + _NAME_STUB) / 2
            out.append(_Way(float(c[0]), float(c[1]), d.direction, d.normal, _NAME_STUB, d.width + 0.4))
    return out


def _name_obstacles(layout: Layout, sses: list[SSE], plates: bool = True) -> tuple[list[tuple], list[tuple]]:
    """What a domain's name must keep clear of: (solid) elements and ghosts, stems with their letters, sheet plates
    (unless not `plates`) and the other names (labels, chain ends and their stubs, stem names, tie letters, domain
    names placed before), and (ways) the straight run of a loop out of each element's end and a short loop's way
    between its two ends. A name's leader keeps off the solid ones and the ways out of the ends (_ways_out): it may
    cross a loop, not stand on a loop's way out of an element's end, which leaves that loop only a route over it
    (random-100 held-out 6ZKT, 7NZ2, 8BCU's domains view: a leader 0.2-0.6 from a port, across its way out)."""
    solid = [_grown(p.rect, _NAME_CLEAR) for p in [*layout.placed.values(), *layout.ghosts]]
    for k, ids in enumerate(layout.sheets):
        if plates and len(ids) > 1:
            r = np.array([layout.placed[i].rect for i in ids] + [g.rect for g in block_ghosts(layout, k)])
            solid.append((*(r[:, :2].min(axis=0) - _PLATE_PAD), *(r[:, 2:].max(axis=0) + _PLATE_PAD)))
    solid += [_grown(b, _NAME_CLEAR) for _, b in label_boxes(layout, sses)]
    ways = _ways_out(layout)
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
    along that edge, then outside with a short leader (Layout.domain_leaders). Clear: of the obstacles with the room
    loops keep (_NAME_KEEP), then without it, and in both with room for a route out of every element's end
    (_NAME_OPEN: by the name and its leader). When nothing is clear, the first place clear of all but the sheet
    plates (a name reads on a plate; random-100 C2b-2: 7JU4's D2 of chain m lies inside a large sheet's plate, and
    its name took the top left, on four loops' ways out of their ends); else, of the places clear but for that room
    (else of all), the one that leaves the ways out the most (random-100 held-out 7NZ2: a name falling back to the
    top left of a panel round one helix stood on its way out)."""
    layout.domain_names, layout.domain_leaders = {}, {}
    fields = domain_fields(layout, names=False)
    if not fields.values:
        return
    em = _ink_em(layout, font_scale) or text_em(layout.width or 1.0) * font_scale
    traces = [np.asarray(pts, float).reshape(-1, 2) for _, pts in sorted(layout.na_traces.items())]
    solid, ways = _name_obstacles(layout, sses)
    bare, _ = _name_obstacles(layout, sses, plates=False)
    outs = _ways_out(layout)  # a leader's too (_name_obstacles)

    letters = (0.0 if layout.na_cropped else NT_ROW + glyph(layout)[1]) + _NAME_CLEAR

    stem_ways = _stem_ways(layout)

    def clear(box, foot, margin, solid=solid, room=True) -> bool:
        """Is the place clear: with `room`, leaving every way out of an end room for a route (_NAME_OPEN)."""
        if overlaps_any(_grown(box, margin), solid + ways) or _trace_hits(_grown(box, _NAME_NEAR_NT), traces):
            return False
        if room and overlaps_any(_grown(box, max(margin, _NAME_OPEN)), outs):  # room for a route out of every end
            return False
        if any(_box_meets_body(_grown(box, letters), d) for d in layout.stems):  # a stem, with its letters
            return False
        if any(_box_meets_body(_grown(box, margin), w) for w in stem_ways):  # a strand's way out of a stem's end
            return False
        if foot is None:
            return True
        (x, ya), (_, yb) = foot
        if overlaps_any((x - 0.05, min(ya, yb), x + 0.05, max(ya, yb)), solid):
            return False
        if not room:
            return True
        # a leader may cross a loop, not stand on a loop's way out of an element's (or a stem's) end: its routing
        # keep-out (label_boxes), with the room a name keeps, off the ways out, and room for a route out at the least
        keep = _grown((x - 0.1, min(ya, yb), x + 0.1, max(ya, yb)), max(margin, _NAME_OPEN - _NAME_ROUTE))
        return not overlaps_any(keep, outs) and not any(_box_meets_body(keep, w) for w in stem_ways)

    def room(box, foot) -> float:
        """The room a place leaves the ways out of the ends: the least gap between its routing keep-outs (the name's,
        its leader's: label_boxes) and a way out (_ways_out), negative where they overlap; at most what a route
        needs (_NAME_OPEN - _NAME_ROUTE), so places with enough room tie."""
        keep = [_grown(box, _NAME_ROUTE)]
        if foot is not None:
            (x, ya), (_, yb) = foot
            keep.append((x - 0.1, min(ya, yb), x + 0.1, max(ya, yb)))
        gaps = [max(a[0] - b[2], b[0] - a[2], a[1] - b[3], b[1] - a[3]) for a in keep for b in outs]
        return min([_NAME_OPEN - _NAME_ROUTE, *gaps])

    for name, value in fields.values.items():
        size = DOMAIN_NAME_SIZE * em
        spots = list(_name_spots(value > 0, fields.xs, fields.ys, size * _text_width(name) + 0.1, size * _NAME_HEIGHT))
        if not spots:
            continue
        spot = next(  # a place with room for the margin loops keep (_NAME_KEEP), else one without; else on a plate
            (
                spot
                for under in (solid, bare)
                for margin in (_NAME_KEEP, 0.0)
                for spot in spots
                if clear(*spot, margin, under)
            ),
            None,
        )
        if spot is None:  # else of the places clear but for that room (else of all), the one leaving the ways most
            spot = max([s for s in spots if clear(*s, 0.0, bare, room=False)] or spots, key=lambda s: room(*s))
        box, foot = spot
        layout.domain_names[name] = box
        if foot is not None:
            layout.domain_leaders[name] = foot
        solid.append(_grown(box, _NAME_CLEAR))


def _name_spots(inside: np.ndarray, xs: np.ndarray, ys: np.ndarray, w: float, h: float):
    """The places a name of w x h may take on a panel (its cells `inside`, rows along ys), in order, as (box,
    leader or None): on the top edge (the name's middle _NAME_SINK inside it, its plate over the edge) at the top
    left, the top right, then along the edge left to right, _NAME_INSET from the panel's sides; then outside, _LEADER
    and twice that above the edge, with a leader down to it; then the same on the bottom edge. The largest piece of
    a panel in pieces first, at each of these steps; on a piece narrower than the name, the name centred over it.
    Last, each place outside again with its leader off the name's middle, at either end of the edge under the name
    (_NAME_INSET in): a leader from the middle of a name over a helix on end stands on the helix's way out, and
    leaves its loop only a route over the name or the leader (random-100 held-out 7NZ2's domains view, D10 · B1)."""
    from scipy.ndimage import label as pieces

    marks, n = pieces(inside)
    order = sorted(range(1, n + 1), key=lambda m: (-int((marks == m).sum()), m))
    rows = np.arange(len(ys))[:, None]
    half = panels.GRID / 2
    later = []  # (box, up, the edge's columns under it and their heights): the places outside, for leaders off middle
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
                        later.append((box, up, xs[span], edge[span]))
                    yield box, foot
    for box, up, sx, se in later:
        mid = float(_snap((box[0] + box[2]) / 2))
        ends = (max(sx.min() + _NAME_INSET, box[0] + 0.05), min(sx.max() - _NAME_INSET, box[2] - 0.05))
        for x in dict.fromkeys(float(_snap(v)) for v in ends):
            if abs(x - mid) > 1e-6 and ends[0] <= ends[1]:
                below = float(_snap(se[np.argmin(np.abs(sx - x))]))
                yield box, ((x, box[1] if up > 0 else box[3]), (x, below))


def provisional(sses: list[SSE]) -> Layout:
    """A layout holding only element labels, for resolving references before anything is placed."""
    placed = [Placed(s, 0.0, 0.0, 1.0, 0.0, 1.0) for s in sses]
    _assign_labels(placed)
    return Layout({p.sse.id: p for p in placed}, (0.0, 0.0, 0.0, 0.0))


def helix_label_pos(p: Placed, k: float | None = None) -> tuple[float, float]:
    """Beside the bar, on its upper side (the other with label_side -1), clear of the body; label_along further
    along it. Its centre stands LABEL_GAP times `k` (glyph_scale of the text drawn; the layout being built's,
    _GLYPH_K, when None) beyond the bar's side, so a name keeps its gap to its bar, in ems, at any Style.font_scale;
    and, where its keep-out box (up to three glyphs, as _footprint's: α10) reaches further across the bar's side
    than LABEL_GAP (a bar turned from the page's x: a name's width lies across it), further by that excess times
    k - 1: the box then reaches no further into the bar, nor towards its ports, than at k = 1, where loops leave
    the bar beside it (plan 22: at a fixed LABEL_GAP, names drawn at font_scale 2 or more sat on their bars; at
    LABEL_GAP times k alone, η1's box on 1UBQ at 3 covered the port column of its vertical bar)."""
    k = _GLYPH_K.get() if k is None else k
    nx, ny = -np.sin(p.angle), np.cos(p.angle)
    across = abs(nx) * (_GLYPH_W * 1.5 + _NAME_PAD) + abs(ny) * _GLYPH_H / 2  # the box's half-reach across the side
    side = (p.width / 2 + LABEL_GAP * k + max(0.0, across - LABEL_GAP) * (k - 1)) * (1 if ny >= 0 else -1)
    side *= p.label_side
    return float(p.cx + nx * side + ny * p.label_along), float(p.cy + ny * side - nx * p.label_along)


def end_label_pos(layout: Layout, end: str, chain: str, port, ex) -> tuple[float, float]:
    """Where a chain end's N or C letter is drawn: END_LABEL out along its stub, or where route._clear_termini moved it
    (Layout.end_labels)."""
    moved = layout.end_labels.get((end, chain))
    if moved is not None:
        return moved
    return float(port[0] + ex[0] * END_LABEL), float(port[1] + ex[1] * END_LABEL)


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
    scale = page_scale(width, height, layout.font_scale)
    return font_pt(scale) * layout.font_scale, 1.0 / (scale * 72)


def line_clearance(layout: Layout) -> float:
    """Page units a line keeps from text: figcheck.LINE_TEXT_CLEARANCE_PX and half a pixel of line width, at the
    page's text scale (page_text; 100 dpi, as figcheck measures), so as many pixels on a wide page of small text as
    on a small one (fix round 1: C14 and C11 kept a share of an em, short of it below 7 pt)."""
    from .figcheck import LINE_TEXT_CLEARANCE_PX

    return (LINE_TEXT_CLEARANCE_PX + 0.5) * 0.72 * page_text(layout)[1]


RUN_EDGE_W = 1.8  # render: the dark edge under residue-coloured loops, in loop widths
NA_RUN_W = 1.4  # render: a nucleic-acid strand's run, in line_pt


def line_pt(scale: float) -> float:
    """The width (points) render draws an element outline at, on a page of `scale` inches per page unit; a loop is
    drawn this times Style.loop_width."""
    return float(np.clip(0.9 * scale / IN_PER_UNIT, 0.5, 1.2))


def line_band(layout: Layout) -> float:
    """How far (page units) a loop's ink reaches from its centre line as figcheck's A6 reads a line by a text
    (text_collisions at clearance 0: half the drawn stroke, line_pt times Layout.loop_stroke, or a nucleic-acid
    strand's NA_RUN_W if wider), plus one pixel at figcheck's 100 dpi to spare (page_text): what a text placed after
    routing keeps from every line (lines_meeting)."""
    unit = page_text(layout)[1]  # page units per point
    stroke = max(layout.loop_stroke, NA_RUN_W if layout.na_paths else 0.0)
    return (line_pt(1.0 / (72 * unit)) * stroke / 2 + 0.72) * unit


def lines_meeting(boxes, segments, band: float) -> np.ndarray:
    """Which segments meet which boxes: a (segment, box) bool array, true where the segment (x0, y0, x1, y1) runs
    within `band` of the box (x0, y0, x1, y1), exact (route.clip_lengths); segments far from every box are not
    clipped. The one line-by-text test of the texts placed after routing: loop and tail counts, helix names
    (route._clear_helix_names), with band line_band (A6's) or wider."""
    from .route import clip_lengths

    box = np.asarray(boxes, float).reshape(-1, 4) + np.array([-band, -band, band, band])
    seg = np.asarray(segments, float).reshape(-1, 4)
    out = np.zeros((len(seg), len(box)), bool)
    if not len(seg) or not len(box):
        return out
    lo, hi = box[:, :2].min(axis=0), box[:, 2:].max(axis=0)  # the boxes' neighbourhood
    near = (np.minimum(seg[:, 0], seg[:, 2]) < hi[0]) & (np.maximum(seg[:, 0], seg[:, 2]) > lo[0])
    near &= (np.minimum(seg[:, 1], seg[:, 3]) < hi[1]) & (np.maximum(seg[:, 1], seg[:, 3]) > lo[1])
    if near.any():
        out[near] = clip_lengths(seg[near], box) > 1e-9
    return out


def spans_crossed(spans, segments, bodies=()) -> np.ndarray:
    """Which spans (x0, y0, x1, y1: from a text to what it names) a line or a body runs across: each span turned to
    lie flat along x, the segments turned with it, then lines_meeting on the flat span (a hair of band: a crossing,
    not a line beside it); a body (x0, y0, x1, y1) by route.clip_lengths. What says a count sits beside what it
    counts, with nothing drawn between (batch C2b-1 fix round 2, 7SMR)."""
    from .route import clip_lengths

    span = np.asarray(spans, float).reshape(-1, 4)
    seg = np.asarray(segments, float).reshape(-1, 4)
    out = np.zeros(len(span), bool)
    for k, (x0, y0, x1, y1) in enumerate(span):
        length = math.hypot(x1 - x0, y1 - y0)
        if length < 1e-9:
            continue
        c, s = (x1 - x0) / length, (y1 - y0) / length
        turn = np.array([[c, -s], [s, c]])  # page -> the span's frame: (x - x0, y - y0) @ turn
        flat = (seg.reshape(-1, 2, 2) - [x0, y0]) @ turn
        out[k] = lines_meeting((0.0, 0.0, length, 0.0), flat.reshape(-1, 4), 1e-6).any()
    rect = np.asarray(bodies, float).reshape(-1, 4)
    if len(rect) and len(span):
        out |= (clip_lengths(span, rect) > 1e-9).any(axis=1)
    return out


def segments_of(lines) -> np.ndarray:
    """The segments (x0, y0, x1, y1) of point lists (loops' and nucleic-acid strands' points), as one array."""
    segs = [
        (*p[k], *p[k + 1]) for pts in lines for p in [np.asarray(pts, float).reshape(-1, 2)] for k in range(len(p) - 1)
    ]
    return np.asarray(segs, float).reshape(-1, 4)


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
        p.label_side_default = p.label_side  # a deliberate flip off a nucleic-acid trace is no later move


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
        if isinstance(copy, Placed):
            copy.label_side_default = copy.label_side


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
    leaving = _leaving_stretches(layout) if block else np.empty((0, 4))
    band = _run_band(layout) if block else 0.0
    for d in layout.stems:
        others = [e for e in layout.stems if e is not d]
        reach = _reach(ink[id(d)])
        rooms, free = [], []
        for s in slots:
            rooms.append(_name_room(d, s, beside, others, named, rects, reach, others_letters[id(d)]))
            # block style: the first pair's numbers off the strands leaving the stems' ends along their rails, as
            # figcheck's A6 reads a line by a text (its page box)
            box = _reach(ink[id(d)] + np.asarray(stem_label_pos(d)))
            free.append(not len(leaving) or not lines_meeting(box, leaving, band).any())
            if _pick_slot(rooms[-1:], _INK_MARGIN * em) == (0, False) and free[-1]:
                break
        k, fallback = next(
            ((k, False) for k, r in enumerate(rooms) if free[k] and _pick_slot([r], _INK_MARGIN * em) == (0, False)),
            _pick_slot(rooms, _INK_MARGIN * em),
        )
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


def _leaving_stretches(layout: Layout) -> np.ndarray:
    """Where single strands leave the stems' ends (x0, y0, x1, y1 rows): along its rail from the stem's outermost
    nucleotide there through the end face and a route's first step (route._PUSH[0]) beyond, as every run, joint,
    hairpin loop, free end or cropped gap leaving that end is drawn, before routing decides where it goes. A strand
    that carries on along the same rail (a bulge) or ends there (a 5' or 3' end with no gap) leaves no stretch. What
    block style's first-pair numbers keep off (random-100 C2b-6 D: on a short stem the numbers, as long as the stem,
    lay beside the rail stretch of a run leaving its end: 8VWQ's P49b, P53b, P95, P141c and 7Y7C's P141 in every
    view)."""
    from .route import _PUSH  # route imports layout: read when called

    na = layout.nucleic
    if na is None or not layout.stems:
        return np.empty((0, 4))
    at = {key: d for d in layout.stems for key in d.axial}
    held = {s: sorted(i for t, i in at if t == s) for s in sorted({t for t, _ in at})}
    gapped = {(g.before, 1) for g in na.gaps} | {(g.after, -1) for g in na.gaps}
    out = []
    for d in layout.stems:
        for key in sorted(d.axial):
            s, i = key
            row = held[s]
            k = row.index(i)
            for way in (1, -1):
                k2 = k + way
                if 0 <= k2 < len(row):
                    e = at[(s, row[k2])]
                    if e is d and d.side[(s, row[k2])] == d.side[key]:
                        continue  # along its rail (a bulge on it)
                elif not (0 <= i + way < len(na.strands[s]) or (s, way) in gapped):
                    continue  # the strand's own end, no gap beyond: its 5' or 3' label, no line
                _, ex = d.port(key, way == 1)
                reach = max(d.length / 2 - d.axial[key] * float(np.dot(ex, d.direction)), 0.0)
                tip = d.point(key) + ex * (reach + _PUSH[0])
                out.append((*d.point(key), *tip))
    return np.asarray(out, float).reshape(-1, 4)


def _run_band(layout: Layout) -> float:
    """How far a nucleic-acid run's ink reaches from its centre line as figcheck's A6 reads a line by a text (as
    line_band, with a run's NA_RUN_W stroke), on the page as it stands."""
    unit = page_text(layout)[1]
    stroke = max(layout.loop_stroke, NA_RUN_W)
    return (line_pt(1.0 / (72 * unit)) * stroke / 2 + 0.72) * unit


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
    """The box around a helix's name, where it stands at `k` (helix_label_pos): `k` (glyph_scale of the text drawn)
    glyphs of _GLYPH_W x _GLYPH_H a letter, and _NAME_PAD times `k` beyond them along the text."""
    return centred_box(*helix_label_pos(p, k), _GLYPH_W * k * len(p.label) / 2 + _NAME_PAD * k, _GLYPH_H * k / 2)


def na_strand_lines(layout: Layout) -> list:
    """Every nucleic-acid line drawn: the strands (Layout.na_paths) and the isolated pairs (Layout.na_ties), a short
    tie's dotted line or a long one's two ticks; for the counts placed after routing to keep off (held-out 6S13: a
    loop's "NN aa" on the dotted line of the tie 448-506)."""
    out = [path.points for path in layout.na_paths]
    for tie in layout.na_ties:
        out += [[a, b] for a, b, _ in tie.marks] if tie.label else [list(tie.ends)]
    return out


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


def count_half(text: str, size: float, pad: float = 0.04) -> tuple[float, float]:
    """Half the width and height of a count's or note's ink box: `text` at `size` page units per em, grown by `pad`
    (the room it is placed with; 0: its ink, as figcheck measures it)."""
    return size * _text_width(text, bold=False) / 2 + pad, size * 0.4 + pad


COUNT_ROOM = 1.0  # page units a count may stand past the layout's bounds and lines (render's page margin, `pad`)


def count_page(layout: Layout, segments) -> tuple[float, float, float, float]:
    """Where a text placed after routing may go: the layout's bounds and every line's extent, and COUNT_ROOM around
    them (render's page holds the counts, growing to them: this keeps one from wandering off the figure)."""
    x0, y0, x1, y1 = layout.bounds
    if len(segments):
        seg = np.asarray(segments, float).reshape(-1, 4)
        x0, x1 = min(x0, seg[:, [0, 2]].min()), max(x1, seg[:, [0, 2]].max())
        y0, y1 = min(y0, seg[:, [1, 3]].min()), max(y1, seg[:, [1, 3]].max())
    return (float(x0) - COUNT_ROOM, float(y0) - COUNT_ROOM, float(x1) + COUNT_ROOM, float(y1) + COUNT_ROOM)


def clear_inks(inks, boxes, segments, layout: Layout, page) -> np.ndarray:
    """Which ink boxes (x0, y0, x1, y1) a count may be drawn at, as figcheck's A6 reads a text: inside `page`, on no
    box (elements, names, letters, earlier counts) within one pixel at figcheck's 100 dpi, and no line within
    line_band (lines_meeting: the exact segments)."""
    ink = np.asarray(inks, float).reshape(-1, 4)
    px = 0.72 * page_text(layout)[1]
    grown = ink + np.array([-px, -px, px, px])
    ok = (ink[:, 0] >= page[0]) & (ink[:, 1] >= page[1]) & (ink[:, 2] <= page[2]) & (ink[:, 3] <= page[3])
    rect = np.asarray(boxes, float).reshape(-1, 4)
    if len(rect):
        ok &= ~(
            (grown[:, None, 0] < rect[None, :, 2])
            & (rect[None, :, 0] < grown[:, None, 2])
            & (grown[:, None, 1] < rect[None, :, 3])
            & (rect[None, :, 1] < grown[:, None, 3])
        ).any(axis=1)
    return ok & ~lines_meeting(ink, segments, line_band(layout)).any(axis=0)


def tail_counts(layout: Layout, sses: list[SSE], loops=(), hidden=None) -> list[dict]:
    """The residue count written beyond each long chain end (at least LONG_TAIL residues past its last element):
    end, chain, text ("14 aa"), centre and ink box (page units). Placed after routing, so not a route keep-out: the
    box is sized from the figure's text (TAIL_COUNT_SIZE em of _ink_em) and set from the N/C letter's far edge,
    along the stub, beside the letter on either side, then further along, then further out beside the letter and along
    the stub (_TAIL_COUNT_FAR), the first beside its letter (within _TAIL_COUNT_NEAR em, no line or body between them:
    spans_crossed) where its ink is clear (clear_inks, as figcheck's A6 reads it) of elements,
    text (label_boxes), the loops' lines (`loops`, route.Loop), the nucleic-acid strands' lines and the counts placed
    before, on the page. A count with no clear spot is left
    out, never drawn over ink: {"kind": "tail", "name": "<end>:<chain>", "text"} is appended to `hidden` (a list),
    which the figure's provenance lists (counts_hidden)."""
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
    segs = segments_of([line.points for line in loops] + na_strand_lines(layout))  # loop and nucleic-acid lines
    page = count_page(layout, segs)
    bodies = [p.rect for p in [*layout.placed.values(), *layout.ghosts, *layout.stems]]
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
        letter = np.asarray(end_label_pos(layout, end, chain, port, ex), float)
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
        ink = count_half(text, size, 0.0)
        ok = clear_inks([centred_box(float(c[0]), float(c[1]), *ink) for c in spots], taken, segs, layout, page)
        # beside its letter (fix round 2, 7SMR): within _TAIL_COUNT_NEAR em of the letter's ink, nothing between
        gap = [max(abs(c[0] - letter[0]) - own[0] - ink[0], abs(c[1] - letter[1]) - own[1] - ink[1]) for c in spots]
        ok &= np.array(gap) <= _TAIL_COUNT_NEAR * em
        ok &= ~spans_crossed([(*letter, *c) for c in spots], segs, bodies)
        box = next((b for b, clear in zip(boxes, ok.tolist()) if clear), None)
        if box is None:  # nowhere clear: left out (batch C2b-1: 7SMR, 5LDW, 8OZJ)
            if hidden is not None:
                hidden.append({"kind": "tail", "name": f"{end}:{chain}", "text": text})
            continue
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
        t = end_label_pos(layout, end, chain, port, ex)
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
    if getattr(layout, "rosette", 0):  # a rosette's axis symbol at the page centre: no line runs over it
        r = AXIS_SYMBOL_R + 0.15
        out.append(("axis", (-r, -r, r, r)))
    return out


_EM_TUNED = 0.55 * SHAFT  # page units per em of the figure's text (before Style.font_scale), on every page
_MIN_SCALE = MIN_TEXT_PT / (_EM_TUNED * 72)  # inches per page unit at which the text is MIN_TEXT_PT
MPL_MIN_PT = 1.0  # matplotlib draws no text smaller than this (FontProperties.set_size floors a size under 1 pt)
SMALLEST_TEXT = 0.6  # ems of the figure's text: the smallest render draws (a ligand's name, a residue number)
_FLOOR_SCALE = MPL_MIN_PT / (SMALLEST_TEXT * _EM_TUNED * 72)  # inches per page unit at which that is MPL_MIN_PT
MAX_GROWN_IN = 200.0  # the longest side that floor grows a page to: PDF readers open no page past 200 in (14 400 pt;
# matplotlib writes no UserUnit), so past it the page is held there and its smallest text falls under MPL_MIN_PT


def page_scale(width: float, height: float = 0.0, font_scale: float = 1.0) -> float:
    """Inches per page unit of a page `width` x `height` page units, as render draws it: IN_PER_UNIT, less on a
    page wider than MAX_WIDTH_IN, but never so little that the text falls under MIN_TEXT_PT (the page grows past
    MAX_WIDTH_IN instead; random-100 C4: on wide pages text held at 4 pt outgrew every keep-out sized in page
    units), nor so much that a side passes MAX_PAGE_IN (past it the whole figure is scaled to fit, its text with
    it: the text keeps its page units, so nothing it keeps clear of changes; only its print size falls), but never
    so far that its smallest text (SMALLEST_TEXT) falls under MPL_MIN_PT: matplotlib would draw it at 1 pt, larger
    than every keep-out sized for it, and the page grows past MAX_PAGE_IN instead (random-100 held-out 7BKC's
    domains view, 2263 x 110 units: held at 48 in its text was 0.5 pt, drawn at 1 pt, 62 texts on loops). The floor
    is for the text as drawn: a Style.font_scale under 1 (the journal theme's 0.9) raises it; and it grows a page to
    MAX_GROWN_IN at most (a PDF reader's largest page), past which the page is held there, its smallest text under
    MPL_MIN_PT (page_fit says so)."""
    scale = max(min(IN_PER_UNIT, MAX_WIDTH_IN / width), _MIN_SCALE)
    longest = max(width, height)
    floor = min(_FLOOR_SCALE / min(font_scale, 1.0), MAX_GROWN_IN / longest)
    return min(scale, max(MAX_PAGE_IN / longest, floor))


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
    whole figure was scaled to fit and its text is smaller (rounded down: 5.97 pt is "5.9 pt", under 6); a page
    past MAX_PAGE_IN (page_scale's MPL_MIN_PT floor), that it is; a page held at MAX_GROWN_IN, that it is and that
    some text is under MPL_MIN_PT (drawn at it, larger than the room kept for it)."""
    if font >= MIN_TEXT_PT - 1e-6:
        return ""
    longest = max(width_in, height_in)
    if longest <= MAX_PAGE_IN + 1e-6:
        held = f"held at {MAX_PAGE_IN:.0f} in on its longer side"
    elif longest < MAX_GROWN_IN - 1e-6:
        held = f"past {MAX_PAGE_IN:.0f} in, so that no text is under {MPL_MIN_PT:.0f} pt"
    else:
        held = (
            f"held at {MAX_GROWN_IN:.0f} in on its longer side (the largest page a PDF reader opens), so some text is "
            f"under {MPL_MIN_PT:.0f} pt and drawn at {MPL_MIN_PT:.0f} pt, larger than the room kept for it"
        )
    return (
        f"the page ({width_in:.0f} x {height_in:.0f} in) is {held}, so its text is "
        f"{math.floor(font * 10) / 10:.1f} pt, under {MIN_TEXT_PT:.0f} pt; draw fewer chains (--chains, --focus) for "
        "larger text"
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
