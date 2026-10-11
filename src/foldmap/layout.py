# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Place sheet blocks and helices on the page (page units: 1 unit = one strand pitch)."""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import NamedTuple

import numpy as np

from .features import helix_bundles
from .frame import Frame, symmetric_frame, view_frame
from .model import SSE, Nucleic, Sheet, group_by_chain
from .nucleic import Stem, StemView, segment_stem, stem_axis, stem_views
from .style import Style

PITCH = 1.1  # distance between neighbouring strand centres
SHAFT = 0.62  # arrow shaft width
SCALE = 0.23  # page units per Å (PITCH / 4.8 Å strand spacing)
HEAD = 1.0  # strand arrowhead width
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
_GLYPH_W, _GLYPH_H = 0.32, 0.42  # rough label glyph size in page units, for keep-out boxes
IN_PER_UNIT = 0.30  # render: inches per page unit, at most
MAX_WIDTH_IN = 7.2  # render: the widest figure (double-column width), inches
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
    ghosts: list[Placed] = field(default_factory=list)
    stems: list[PlacedStem] = field(default_factory=list)
    nucleic: Nucleic | None = None
    sheets: list[list[str]] = field(default_factory=list)  # strand ids of each sheet block, in layout order
    domains: list[tuple[str, list[str]]] = field(default_factory=list)  # named domains and their element ids
    na_paths: list = field(default_factory=list)  # route.NAPath: single strands as drawn (filled by route_loops)
    na_traces: dict = field(default_factory=dict)  # strand -> projected backbone (page units), for unpaired strands
    na_ends: list = field(default_factory=list)  # (strand, "5′" or "3′", label position), filled by route_loops
    na_letters: dict = field(default_factory=dict)  # single-stranded (strand, nt) -> its letter's place (route_loops)
    na_ties: list = field(default_factory=list)  # na_paths.NATie: isolated pairs as drawn (route_loops)
    na_cropped: bool = False  # na_extent auto cropped a large nucleic acid to the stems the protein touches
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
    page_axes: list | None = None  # 3x3 rows: page-right, page-up, toward-viewer (for the first protomer)
    page_centre: list | None = None  # the 3D point those axes are taken around
    uniprot: object = None  # uniprot.Annotation when UniProt features were requested
    membrane: dict | None = None  # {"y": (bottom, top) of the band, "inside": "below"|"above", "source": ...}
    res_cons: list[float] = field(default_factory=list)  # alignment conservation of each residue (nan: none)
    res_grade: list[int] = field(default_factory=list)  # conservation grade 1..9 of each residue (0: none)
    cons_label: str = ""  # legend text saying where the grades come from

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


def _footprint(m) -> tuple[float, float, float, float]:
    """A member's box relative to its own centre; a helix's or a stem's includes the room its label needs."""
    x0, y0, x1, y1 = m.rect
    box = (x0 - m.cx, y0 - m.cy, x1 - m.cx, y1 - m.cy)
    if isinstance(m, Placed) and m.sse.kind in ("H", "G"):
        lx, ly = helix_label_pos(Placed(m.sse, 0.0, 0.0, m.length, m.angle, m.width))
        hx, hy = _GLYPH_W * 1.5 + 0.05, _GLYPH_H / 2  # up to three glyphs: α10
        box = (min(box[0], lx - hx), min(box[1], ly - hy), max(box[2], lx + hx), max(box[3], ly + hy))
    if isinstance(m, PlacedStem):  # a stem's name beside it
        x0, y0, x1, y1 = _stem_label_box(m)
        box = (min(box[0], x0 - m.cx), min(box[1], y0 - m.cy), max(box[2], x1 - m.cx), max(box[3], y1 - m.cy))
    if isinstance(m, PlacedStem) and m.hairpin >= 0:  # and its hairpin loop over the tip, with the loop's letters
        a, i, b, j = m.pairs[-1]
        (x, y), r = m.hairpin_circle((a, i), (b, j), m.hairpin)
        x, y, r = x - m.cx, y - m.cy, r + NT_ROW + _GLYPH_H / 2
        box = (min(box[0], x - r), min(box[1], y - r), max(box[2], x + r), max(box[3], y + r))
    if isinstance(m, Placed) and m.ends:
        d = m.direction
        for end in m.ends:
            out = d if end == "C" else -d
            t = out * (m.length / 2 + END_LABEL)  # the N/C letter beyond the stub
            hx, hy = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2
            box = (min(box[0], t[0] - hx), min(box[1], t[1] - hy), max(box[2], t[0] + hx), max(box[3], t[1] + hy))
    return box


def _strand_length(sse: SSE) -> float:
    return max(_MIN_STRAND, (len(sse) - 1) * STRAND_RISE * SCALE)


def _anchor(strands: list[SSE]) -> tuple[int, str]:
    first = min(strands, key=lambda s: s.start)
    return first.start, first.chain


def _sheet_item(sheet: Sheet, frame: Frame, pitch: float = PITCH) -> _Item:
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
    members, offsets = [], []
    slots = len(order) + (1 if sheet.closed else 0)
    for slot, k in enumerate(order):
        up = sheet.directions[k] * s > 0
        sse = strands[k]
        members.append(Placed(sse, 0, 0, _strand_length(sse), np.pi / 2 if up else -np.pi / 2, pitch))
        offsets.append(((slot - (slots - 1) / 2) * pitch, 0.0))
    if sheet.closed:  # the barrel wraps: its first strand pairs again with the last
        f = members[0]
        members.append(Placed(f.sse, 0, 0, f.length, f.angle, f.width, ghost=True))
        offsets.append(((slots - 1) / 2 * pitch, 0.0))
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


def _overlapping(xy: np.ndarray, halves: np.ndarray, mids: np.ndarray, apart: np.ndarray | None = None) -> np.ndarray:
    """(i, j) with i < j, in row order, for every pair of footprints (plus MARGIN) that overlap; pairs marked in
    the boolean matrix `apart` (_joint_matrix) are never counted."""
    q = xy + mids
    depth = halves[:, None, :] + halves[None, :, :] + MARGIN - np.abs(q[None, :, :] - q[:, None, :])
    hit = (depth[..., 0] > 1e-9) & (depth[..., 1] > 1e-9)
    if apart is not None:
        hit &= ~apart
    return np.argwhere(np.triu(hit, 1))


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
) -> None:
    """Push overlapping items apart, then slide any left over clear. With `ties` (_symmetry_ties) the pushes
    keep every protomer drawn alike (_push_apart), and that result is kept if it clears every overlap within
    _SEP_ORBIT_SWEEPS sweeps and its box takes at most _SEP_SYM_ROOM times the area of moving the items one by
    one from the same start (widening the spread can otherwise cost far more room); else the one-by-one
    result is used. Consecutive pieces of one bent duplex (`joints`, _joints) in a short chain (_rigid_joints)
    are moved together and never pushed apart. Items are pushed off the `fixed` boxes they must keep off (Walls),
    which stay where they are."""
    halves = np.array([it.half for it in items], float).reshape(-1, 2)
    mids = np.array([it.mid for it in items], float).reshape(-1, 2)
    start = np.array([(it.cx, it.cy) for it in items], float).reshape(-1, 2)
    home = np.array([it.home for it in items], float).reshape(-1, 2)
    joints = _rigid_joints(items, joints)
    apart = _joint_matrix(len(items), joints)
    xy = start.copy()
    alone = start.copy()
    _push_apart(alone, halves, mids, home, {}, _SEP_SWEEPS, apart, fixed)
    together = ties and _push_apart(xy, halves, mids, home, ties, _SEP_ORBIT_SWEEPS, apart, fixed)
    if not together or _box_area(xy, halves, mids) > _SEP_SYM_ROOM * _box_area(alone, halves, mids):
        xy = alone
    for it, (x, y) in zip(items, xy.tolist()):
        it.cx, it.cy = x, y
    _slide_clear(items, joints, fixed)


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
) -> bool:
    """Push overlapping footprints apart in place (axis of least penetration); spring back to home early on.
    Each sweep finds the overlapping pairs at once and resolves them in order. With `ties` a push moves an
    item's whole orbit (it and its symmetry copies) by the same amount, and a sideways push between two
    protomers moves protomer k by k times the amount when that widens the spread between protomers (also
    between two copies of one element); otherwise two copies of one element, which no orbit-wide move can
    part, move alone. Items joined by `apart` (consecutive pieces of one bent duplex) are never pushed apart: a push
    moves the whole chain of pieces (with their orbits); two items of one such group that overlap move alone.
    An item overlapping a `fixed` box it must keep off (Walls) is pushed (with its group) the whole way off it.
    True if no overlap is left."""
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
    walls = fixed
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
            if px < py:
                s = (1 if dx > 0 else -1) if dx != 0 else 1
                if None not in (protomer[i], protomer[j]) and protomer[i] != protomer[j]:
                    wider = s * (px + 2 * _SEP_PUSH) / (protomer[j] - protomer[i])
                    if wider * widen > 0:  # items without copies (copy_k 0) stay put, whichever protomer
                        xy[:, 0] += copy_k * wider
                        continue
                xy[gi, 0] -= s * (px / 2 + _SEP_PUSH)
                xy[gj, 0] += s * (px / 2 + _SEP_PUSH)
            else:
                s = (1 if dy > 0 else -1) if dy != 0 else 1
                xy[gi, 1] -= s * (py / 2 + _SEP_PUSH)
                xy[gj, 1] += s * (py / 2 + _SEP_PUSH)
        _off_walls(xy, halves, mids, walls, group)
        if sweep < 100:
            xy += 0.03 * (home - xy)
    return not len(_overlapping(xy, halves, mids, apart)) and not _off_walls(xy, halves, mids, walls, group, move=False)


def _off_walls(xy, halves, mids, walls: Walls | None, group, move: bool = True) -> bool:
    """Push each footprint overlapping a fixed box it must keep off (Walls), MARGIN kept, the whole way off it
    along the axis of least overlap, moving its group (`group[i]`, the items moved with item i) in place. True if
    any overlapped (with move=False, only that is reported)."""
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
        if px < py:
            xy[group[i], 0] -= ((1 if dx > 0 else -1) if dx != 0 else 1) * (px + _SEP_PUSH)
        else:
            xy[group[i], 1] -= ((1 if dy > 0 else -1) if dy != 0 else 1) * (py + _SEP_PUSH)
    return True


def _slide_clear(items: list[_Item], joints: set[tuple[int, int]] = frozenset(), fixed: Walls | None = None) -> None:
    """Guarantee: each item in turn, if it still overlaps one already settled (other than a piece it is joined
    to, `joints`) or a `fixed` box it must keep off, takes the nearest clear spot straight right, left, up or down
    (in that order on a tie), in _SLIDE_STEP steps."""
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
        while any(_collide(it, halves[i], mids[i], *other) for other in settled):
            k += 1
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
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
) -> None:
    """Minimise the layout energy (_energy) from the projected positions, raising the overlap penalty in
    stages. The result becomes each item's new home; _separate then guarantees no overlap remains. Every input
    is snapped to the grid (_snap) and _lbfgs runs without BLAS, so the minimum is the same on every machine.
    `extra`: further springs (item i, port offset, item j, port offset, weight, slack), the nucleic-acid strands
    from stem to stem (_stem_springs). `joints` (_joints): pairs of items kept apart without MARGIN, so their
    boxes may touch and no element fits between them for long. `fixed` (Walls): boxes that stay where they are (the
    nucleotides of unpaired strands' traces), the items kept MARGIN off them, firmly or (an element touching that
    trace) not at all."""
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
    ties = ties or {}
    free = [i for i in range(n) if i not in ties]
    slot = {i: k for k, i in enumerate(free)}
    base = np.array([slot[ties[i][0]] if i in ties else slot[i] for i in range(n)])
    copy_k = np.array([ties[i][1] if i in ties else 0 for i in range(n)], float)

    def expand(z):
        P = z[:-1].reshape(-1, 2)[base].copy()
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


def _stack(items: list[_Item]) -> None:
    """Rows in sequence order, one row break per chain and when a row is full."""
    items.sort(key=lambda it: (it.chain, it.order))
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
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    counts: dict[tuple[str, str], int] = {}
    for p in sorted(placed, key=lambda p: (p.sse.chain, p.sse.start)):
        n = counts.get((p.sse.chain, p.sse.kind), 0)
        counts[(p.sse.chain, p.sse.kind)] = n + 1
        if p.sse.kind == "E":
            p.label = letters[n] if n < 26 else f"S{n + 1}"
        elif p.sse.kind == "G":
            p.label = f"η{n + 1}"
        else:
            p.label = f"α{n + 1}"


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


def _symmetry_ties(items: list[_Item], frame) -> dict[int, tuple[int, int]]:
    """Copy item -> (its protomer-0 mate, protomer index k): matched by chain class and shared residue numbers.
    The copy also takes the mate's shape and orientation, so every protomer is drawn the same."""
    sector, classes = getattr(frame, "sector", None), getattr(frame, "classes", None)
    if not sector or not classes:
        return {}

    def residues(it):  # chains outside every protomer (bound peptides, unpaired chains) take no part
        return {
            (classes[m.sse.chain], r)
            for m in it.members
            if isinstance(m, Placed) and not m.ghost and m.sse.chain in classes
            for r in range(m.sse.first.seq, m.sse.last.seq + 1)
        }

    def protomer(it):
        ks = {sector.get(m.sse.chain) for m in it.members if isinstance(m, Placed)}
        return ks.pop() if len(ks) == 1 else None

    sets = [residues(it) for it in items]
    zero = [i for i, it in enumerate(items) if protomer(it) == 0 and sets[i]]
    ties = {}
    for i, it in enumerate(items):
        k = protomer(it)
        if not k or not sets[i]:
            continue
        best = max(zero, key=lambda j: len(sets[i] & sets[j]) / len(sets[i] | sets[j]), default=None)
        if best is None or len(sets[i] & sets[best]) < 0.5 * len(sets[i] | sets[best]):
            continue
        mate = items[best]

        def span(m):
            if m.sse.chain not in classes:
                return set()
            return {(classes[m.sse.chain], r) for r in range(m.sse.first.seq, m.sse.last.seq + 1)}

        pairs = [
            max(
                range(len(mate.members)),
                key=lambda q: (len(span(m) & span(mate.members[q])), mate.members[q].ghost == m.ghost),
            )
            for m in it.members
        ]
        if len(mate.members) == len(it.members) and sorted(pairs) == list(range(len(pairs))):
            for m, q in zip(it.members, pairs):  # one-to-one: take the mate's shape member by member
                m.angle, m.length = mate.members[q].angle, mate.members[q].length
            it.offsets = [mate.offsets[q] for q in pairs]
        else:  # DSSP split or merged an element in one copy: keep the shape, match the orientation
            votes = [
                np.cos(m.angle - mate.members[q].angle)
                for m, q in zip(it.members, pairs)
                if span(m) & span(mate.members[q])
            ]
            if votes and np.mean(votes) < 0:  # drawn upside down: turn the block 180 degrees in the page
                for m in it.members:
                    m.angle += np.pi
                it.offsets = [(-ox, -oy) for ox, oy in it.offsets]
        ties[i] = (best, k)
    return ties


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
) -> None:
    """The user's own hand on the layout: new labels, exchanged places, nudges; then overlaps are cleared
    (consecutive pieces of a bent duplex, `joints`, may still meet; items keep off the `fixed` boxes)."""
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
        _separate(items, None, joints, fixed)
        for it in items:
            it.apply()


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
) -> Layout:
    """`attract` (projected mode) pulls loop-joined and touching elements together before overlaps
    are removed; `contacts` comes from features.sse_contacts. `nucleic` adds its stems (straight pieces of
    base-paired DNA/RNA), placed like any element and pulling the elements that touch them; strands with no pair
    at all are kept as projected traces."""
    if mode not in ("projected", "stack"):
        raise ValueError(f"unknown layout mode {mode!r}")
    style = style or Style()
    angle_mode = "upright" if mode == "stack" else style.helix_angle
    items = [_sheet_item(sh, frame, PITCH * style.strand_scale) for sh in sheets]
    bundles = helix_bundles(sses)
    bundled = {k for g in bundles for k in g}
    by_id = {s.id: s for s in sses}
    items += [_bundle_item([by_id[k] for k in g], frame, angle_mode, HELIX_W * style.helix_scale) for g in bundles]
    items += [
        _helix_item(s, frame, angle_mode, HELIX_W * style.helix_scale)
        for s in sses
        if s.kind == "H" and s.id not in bundled
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
    if mode == "stack":
        _stack(items)
        _centre_dna_over_binders(items, contacts)
    else:
        ties = _symmetry_ties(items, frame)
        spread, lift = 0.0, 0.0
        if ties:  # protomers side by side at least one protomer's width apart; copies start on their mates
            mates = {j for j, _ in ties.values()}
            x0 = min(items[j].cx + items[j].box[0] for j in mates)
            x1 = max(items[j].cx + items[j].box[2] for j in mates)
            # protomer k sits k spreads along: rightwards when the copies turn anticlockwise seen from the top
            # of the axis, leftwards for a left-handed filament (or an assembly turned over), as projected
            spread = float(np.sign(frame.period) or 1.0) * max(abs(frame.period) * SCALE, x1 - x0 + 2 * MARGIN)
            lift = getattr(frame, "rise", 0.0) * SCALE
            for i, (j, k) in ties.items():
                items[i].cx, items[i].cy = items[j].cx + k * spread, items[j].cy + k * lift
                items[i].home = (items[j].home[0] + k * spread, items[j].home[1] + k * lift)
        if attract and len(items) > 1:
            extra = _stem_springs(items, nucleic) if pieces else []
            rise = getattr(frame, "rise", 0.0) * SCALE
            _attract(items, sses, contacts, ties, spread, rise, bridges or [], extra=extra, joints=joints, fixed=walls)
        _separate(items, ties, joints, walls)
        if domains and _part_domains(items, domains):
            _separate(items, ties, joints, walls)
    for it in items:
        it.apply()
    dna = [m for it in items for m in it.members if isinstance(m, PlacedStem)]
    placed = [m for it in items for m in it.members if isinstance(m, Placed) and not m.ghost]
    ghosts = [m for it in items for m in it.members if isinstance(m, Placed) and m.ghost]
    _assign_labels(placed)
    real = {p.sse.id: p for p in placed}
    if rename or swap or move:
        _adjust(items, Layout(real, (0.0, 0.0, 0.0, 0.0)), sses, rename or {}, swap or [], move or [], joints, walls)
    for g in ghosts:
        g.label = real[g.sse.id].label
    rects = np.array([p.rect for p in placed + ghosts + dna]).reshape(-1, 4)
    sheet_blocks = [
        [m.sse.id for m in it.members if not m.ghost]
        for it in items
        if it.members and all(isinstance(m, Placed) and m.sse.kind == "E" for m in it.members)
    ]
    lay = Layout(
        real,
        (0.0, 0.0, 0.0, 0.0),
        ghosts,
        dna,
        nucleic if dna or traces else None,
        sheet_blocks,
        bool(dna or traces) and style.nucleotide_labels,
    )
    lay.na_traces = traces
    lay.frame = frame
    if traces:
        pts = np.vstack(list(traces.values()))
        rects = np.vstack([rects, [[pts[:, 0].min(), pts[:, 1].min(), pts[:, 0].max(), pts[:, 1].max()]]])
    lay.partial = set(partial or ())
    lay.domains = [(name, [k for k in ids if k in real]) for name, ids in domains or []]
    rects = np.vstack([rects, *[[domain_panel(lay, name)] for name, _ in lay.domains]]) if lay.domains else rects
    rects = np.vstack([rects, [b for _, b in label_boxes(lay, sses)] or np.empty((0, 4))])
    lay.bounds = (
        float(rects[:, 0].min()),
        float(rects[:, 1].min()),
        float(rects[:, 2].max()),
        float(rects[:, 3].max()),
    )
    return lay


DOMAIN_PAD = 0.7  # panel margin around a domain's elements
_DOMAIN_GAP = 2 * DOMAIN_PAD + 0.6  # clear space between two domains' elements: both panels plus a gutter


def _part_domains(items: list[_Item], domains: list[tuple[str, list[str]]], sweeps: int = 60) -> bool:
    """Move whole domains apart (each as a rigid group) until their panels no longer overlap. Each item goes
    with the domain most of its elements belong to. True if anything moved."""
    owner = {k: d for d, (_, ids) in enumerate(domains) for k in ids}
    groups: dict[int, list[_Item]] = {}
    for it in items:
        votes = [owner[m.sse.id] for m in it.members if isinstance(m, Placed) and m.sse.id in owner]
        if votes:
            groups.setdefault(max(set(votes), key=votes.count), []).append(it)
    keys = sorted(groups)
    if len(keys) < 2:
        return False

    def box(d):
        b = np.array([(it.cx + it.box[0], it.cy + it.box[1], it.cx + it.box[2], it.cy + it.box[3]) for it in groups[d]])
        return b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max()

    moved = False
    for _ in range(sweeps):
        clash = False
        for i, a in enumerate(keys):
            for b in keys[i + 1 :]:
                A, B = box(a), box(b)
                ox = min(A[2], B[2]) - max(A[0], B[0]) + _DOMAIN_GAP
                oy = min(A[3], B[3]) - max(A[1], B[1]) + _DOMAIN_GAP
                if ox <= 0 or oy <= 0:
                    continue
                clash = moved = True
                ca, cb = (
                    np.array([(A[0] + A[2]) / 2, (A[1] + A[3]) / 2]),
                    np.array([(B[0] + B[2]) / 2, (B[1] + B[3]) / 2]),
                )
                if ox <= oy:  # the cheaper direction, each side going half way
                    step = np.array([ox / 2 if cb[0] >= ca[0] else -ox / 2, 0.0])
                else:
                    step = np.array([0.0, oy / 2 if cb[1] >= ca[1] else -oy / 2])
                for it, sign in [(it, -1) for it in groups[a]] + [(it, 1) for it in groups[b]]:
                    it.cx, it.cy = it.cx + sign * step[0], it.cy + sign * step[1]
                    it.home = (it.home[0] + sign * step[0], it.home[1] + sign * step[1])  # rest here from now on
        if not clash:
            break
    return moved


def domain_panel(layout: Layout, name: str) -> tuple[float, float, float, float]:
    """A domain's panel (its elements plus margin) with room above for its name."""
    ids = dict(layout.domains)[name]
    r = np.array([layout.placed[k].rect for k in ids] + [g.rect for g in layout.ghosts if g.sse.id in ids])
    x0, y0, x1, y1 = r[:, 0].min(), r[:, 1].min(), r[:, 2].max(), r[:, 3].max()
    return (float(x0 - DOMAIN_PAD), float(y0 - DOMAIN_PAD), float(x1 + DOMAIN_PAD), float(y1 + DOMAIN_PAD + 0.6))


def provisional(sses: list[SSE]) -> Layout:
    """A layout holding only element labels, for resolving references before anything is placed."""
    placed = [Placed(s, 0.0, 0.0, 1.0, 0.0, 1.0) for s in sses]
    _assign_labels(placed)
    return Layout({p.sse.id: p for p in placed}, (0.0, 0.0, 0.0, 0.0))


def helix_label_pos(p: Placed) -> tuple[float, float]:
    """Beside the bar, on its upper side, clear of the body."""
    nx, ny = -np.sin(p.angle), np.cos(p.angle)
    side = (p.width / 2 + LABEL_GAP) * (1 if ny >= 0 else -1)
    return float(p.cx + nx * side), float(p.cy + ny * side)


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
    """Where a stem's name goes: beside it on its first rail's side, a glyph beyond its base letters."""
    p = np.array([d.cx, d.cy]) + d.normal * (d.width / 2 + NT_ROW + _GLYPH_H)
    return float(p[0]), float(p[1])


def _stem_label_box(d: PlacedStem) -> tuple[float, float, float, float]:
    """The box around a stem's name, written along the stem."""
    along, across = _GLYPH_W * len(d.name) / 2 + 0.05, _GLYPH_H / 2
    c, s = abs(np.cos(d.angle)), abs(np.sin(d.angle))
    return centred_box(*stem_label_pos(d), c * along + s * across, s * along + c * across)


def label_boxes(layout: Layout, sses: list[SSE]) -> list[tuple[str, tuple[float, float, float, float]]]:
    """Keep-out boxes for text and chain-end stubs: loops route around them, the page bounds include them."""
    out = []
    for p in layout.placed.values():
        if p.sse.kind in ("H", "G"):
            out.append(
                (
                    f"label:{p.sse.id}",
                    centred_box(*helix_label_pos(p), _GLYPH_W * len(p.label) / 2 + 0.05, _GLYPH_H / 2),
                )
            )
    for end, chain, port, ex in termini(layout, sses):
        a, b = np.add(port, ex * 0.15), np.add(port, ex * END_STUB)
        out.append(
            (
                f"{end}-stub:{chain}",
                (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05),
            )
        )
        t = np.add(port, ex * END_LABEL)
        out.append((f"{end}:{chain}", centred_box(float(t[0]), float(t[1]), _GLYPH_W / 2 + 0.05, _GLYPH_H / 2)))
    for k, d in enumerate(layout.stems):
        out.append((f"stem-label:{k}", _stem_label_box(d)))
    for n, tie in enumerate(layout.na_ties):  # a long tie's labels at its two marks
        for end, (_, _, (x, y)) in enumerate(tie.marks):
            out.append((f"tie-label:{n}:{end}", centred_box(x, y, _GLYPH_W / 2 + 0.05, _GLYPH_H / 2)))
    for name, ids in layout.domains:  # a domain's name, top left inside its panel
        if ids:
            x0, _, _, y1 = domain_panel(layout, name)
            out.append(
                (
                    f"domain-label:{name}",
                    (x0 + 0.25, y1 - 0.35 - _GLYPH_H * 1.1, x0 + 0.45 + _GLYPH_W * 1.1 * len(name), y1 - 0.2),
                )
            )
    return out


def font_pt(scale: float) -> float:
    """The figure's text size in points (before Style.font_scale) at `scale` inches per page unit."""
    return float(np.clip(0.55 * SHAFT * scale * 72, 4.0, 9.0))


def text_em(width: float) -> float:
    """How many page units one em of the figure's text takes on a page `width` units wide, as render sizes it
    (before Style.font_scale): a wide page is scaled down to MAX_WIDTH_IN and its text, held at 4 pt or more,
    takes more page units."""
    scale = min(IN_PER_UNIT, MAX_WIDTH_IN / width)
    return font_pt(scale) / (scale * 72)


NA_TORN_TURN = 90.0  # degrees: linked stems turned against one another by more than this when unrolled are torn


def na_torn_by_unrolling(nucleic, sses, sym) -> bool:
    """Would unrolling the assembly (symmetric_frame) tear the nucleic acid apart, against the plain projected frame
    (view_frame) the figure would otherwise use? Each stem is drawn whole where its centre is seen, so a single stem
    cannot be torn; what can: (a) a duplex split into pieces whose consecutive pieces' centres face different
    protomers (they land on either side of a seam: 1RVA's kinked site, drawn as a 180° hairpin), or (b) a group of
    stems linked through at most _SHORT_RUN nucleotides (a junction, a bent duplex, coaxial stacks) whose axes, after
    the best common turn, are turned against one another by more than NA_TORN_TURN between the two frames (4LD0's
    junction arms: 139°; 30QM's linked stems at most 57°)."""
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
