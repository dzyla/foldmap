# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Page checks for any figure, protein or nucleic acid, as nacheck is for nucleic acids: element overlaps, text
against text, bodies and lines, text off the page, labels repeated in a chain, copies drawn alike, domain panels and
names, stem names by their own letters, chains or nucleotides left out without a note, and N->C and sheet order.

Every check takes what spec.draw_spec returns (spec.Drawn: figure, loops, layout, sses, bb, model) or its parts.
Page measures are in page units (the axes' data units); pixel measures use the figure as drawn, at its own dpi.
Artists are read by their gid only (render gives every artist one)."""

from __future__ import annotations

import math
import re
from collections import Counter
from itertools import combinations

import matplotlib
import numpy as np
from matplotlib.path import Path as MPath
from matplotlib.transforms import Affine2D, Bbox

from .layout import block_ghosts, crosses_band, joined_pieces, lying_over
from .membrane import (
    BAND_EPS,
    band_ends,
    band_set,
    crosses,
    crossing_excess,
    loop_crossing,
    loop_excess,
    loop_lanes,
    loop_side_bound,
    own_bound,
    slab_helix,
    slab_range,
)
from .membrane import band_eps as _eps
from .membrane import nearest_band as _nearest_band
from .membrane import page_depth as _drawn_depth
from .na_curve import _FAR, on_curves

# Lines a text can be crossed by (render's gids, matched by prefix). A loop's casing (loopcase:) and a stub's edge
# (stub-edge:) are drawn under their line in the page or a pale colour: neither is read. A chain-end stub is its
# terminus letter's and its tail length's line, a domain leader its name's (Page.own_lines): read against every
# other text (final review B, Minor 1).
LINE_GIDS = (
    "loop:",
    "loop-seg:",
    "loop-arrow:",
    "na:run:",
    "na:tie:",
    "ligand-tether:",
    "disulfide:",
    "hbond:",
    "ia:",
    "stub:",
    "domain-leader:",
    "absent:",  # foldmap diff: the dashed outline of an element only the other structure has (render._draw_absent)
)
GLYCAN_STEM = re.compile(r"glycan:[^:]+:stem$")
# Bodies: opaque shapes a text must not sit on unless it is theirs, and that hide a line drawn under them. Element
# bodies are owned by their element id (label:<id>, resnum:<id>:<end>, ghostlabel:<id> are theirs).
ELEMENT_BODIES = ("strand:", "helix:", "helix-back:", "eta:", "ghost:")
OWN_BODIES = {"ligand:": "ligand-label:", "linked:": "linked-label:", "domain-name-plate:": "domain-label:"}
DATA_BODIES = ("data-ring:",)  # a category mark of the researcher's data (render.data_marks): a text on one hides it
BESIDE_BODIES = ("helix:", "helix-back:", "eta:")  # bodies whose name stands beside them, never on them: A6 reads
# their own name over them (plan 22: names over their own bars at font_scale 2-3.5 went unread, owned)
TEXT_OVERLAP_PX2 = 1.0  # text over text or a body: overlaps up to this many square pixels are anti-aliasing
TEXT_DEPTH_PX = 2.0  # ... and a text or body reaching no deeper than this into a text's ink is a near-touch (R13):
# the same scale as LINE_TEXT_CLEARANCE_PX, wider than what matplotlib versions disagree by
OFF_PAGE_PX = 1.0  # text may reach this far past the page edge
PANEL_GRID = 0.05  # page units: the raster domain panel outlines are measured on
_STEP_PX = 0.5  # sampling step (pixels) along a line inside a text box, and across a text box over a body
_MPL = tuple(int(x) for x in re.findall(r"\d+", matplotlib.__version__)[:2])
_AUTO_DOMAIN = re.compile(r"D\d+( · .+)?")  # features.auto_domains names: D1, D2 ... with " · chain" when split


# ------------------------------------------------------------------------------------------- loop crossings


def loop_crossings(loops) -> int:
    """Proper crossings between the polylines of different loops (shared end points do not count)."""

    def cross(p, q, r, s) -> bool:
        def side(a, b, c):
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

        return side(r, s, p) * side(r, s, q) < 0 and side(p, q, r) * side(p, q, s) < 0

    segs = [list(zip(lp.points, lp.points[1:])) for lp in loops]
    flat = [(i, p, q) for i, ss in enumerate(segs) for p, q in ss]
    if not flat:
        return 0
    own = np.array([i for i, _, _ in flat])
    box = np.array([(min(p[0], q[0]), min(p[1], q[1]), max(p[0], q[0]), max(p[1], q[1])) for _, p, q in flat], float)
    n = 0
    for k, (i, p, q) in enumerate(flat):  # only segments of later loops whose boxes meet this one's can cross it
        later = np.flatnonzero(own > i)
        near = later[
            (box[later, 0] <= box[k, 2])
            & (box[k, 0] <= box[later, 2])
            & (box[later, 1] <= box[k, 3])
            & (box[k, 1] <= box[later, 3])
        ]
        n += sum(cross(p, q, flat[m][1], flat[m][2]) for m in near.tolist())
    return n


def honest_crossings(lay, loops) -> int:
    """Every crossing a reader sees: loop x loop (loop_crossings), plus each loop drawn through the body of an element
    it does not join, a helix drawn over a sheet included (a loop passing under a helix not its own), except a loop
    leaving a helix drawn over a sheet across the strand its port lies on (Task 14 review)."""
    from .route import _crosses_body

    n = loop_crossings(loops)
    placed = list(lay.placed.items())
    rects = np.array([e.rect for _, e in placed], float).reshape(-1, 4)
    for lp in loops:
        segs = list(zip(lp.points, lp.points[1:]))
        a, b = lay.placed.get(lp.a_id), lay.placed.get(lp.b_id)
        ports = [p for q, p in ((a, lp.points[0]), (b, lp.points[-1])) if q is not None and q.over]
        pts = np.asarray(lp.points, float).reshape(-1, 2)
        x0, y0 = pts.min(axis=0) - 1e-6
        x1, y1 = pts.max(axis=0) + 1e-6
        near = (rects[:, 0] <= x1) & (x0 <= rects[:, 2]) & (rects[:, 1] <= y1) & (y0 <= rects[:, 3])
        for m in np.flatnonzero(near).tolist():  # an element whose box the loop's extent misses is not crossed
            k, e = placed[m]
            if k in (lp.a_id, lp.b_id) or not any(_crosses_body(p, q, e) for p, q in segs):
                continue
            x0, y0, x1, y1 = e.rect
            if e.sse.kind == "E" and any(x0 <= x <= x1 and y0 <= y <= y1 for x, y in ports):
                continue
            n += 1
    return n


# ------------------------------------------------------------------------------------------- element boxes


COLUMN_TOL = 1e-3  # page units (in strand pitches across a sheet, _across): two strands this close share a column


def sheet_columns(values) -> dict[float, int]:
    """Each value across a sheet (_across) -> its column, counted from 0 in increasing order: a value within
    COLUMN_TOL of its column's first value shares that column (3-decimal rounding split 3KGS's 23.2525 +- 8e-6 in
    two); measured from the first, never the one before, so a run of close values cannot drift a column wide."""
    place, n, first = {}, -1, None
    for x in sorted(set(values)):
        if first is None or x - first > COLUMN_TOL:
            n, first = n + 1, x
        place[x] = n
    return place


ROSETTE_TURN_TOL = 1e-3  # page units: a turned copy back within this of its mate stands exactly where its mate does


def local_layout(lay):
    """A rosette's layout (Layout.rosette) with each element and ghost turned back by its copy's turn about the page
    centre (frame.angles[k] of its chain's protomer): every copy then stands as protomer 0 does, in its own frame,
    where a sheet's columns run across the page and its strands up and down as the layout drew protomer 0's. The
    checks that read a sheet block in page axes (science_checks: columns and spans) read it so; a block is one copy's
    (a rosette's sheets are within protomers). Any other layout as it is."""
    if not getattr(lay, "rosette", 0):
        return lay
    import copy as _copy

    sector, angles = lay.frame.sector, lay.frame.angles
    mates = getattr(lay, "rosette_mates", {}) or {}
    ghost_of = {g.sse.id: g for g in lay.ghosts if sector.get(g.sse.chain) == 0}

    def back(p):
        k = sector.get(p.sse.chain)
        if not k:
            return p
        c, s = math.cos(-angles[k]), math.sin(-angles[k])
        q = _copy.copy(p)
        q.cx, q.cy = c * p.cx - s * p.cy, s * p.cx + c * p.cy
        q.angle = p.angle - angles[k]
        mate = mates.get(p.sse.id, (None,))[0]
        m = (ghost_of if p.ghost else lay.placed).get(mate)
        if (
            m is not None
            and math.hypot(q.cx - m.cx, q.cy - m.cy) < ROSETTE_TURN_TOL
            and math.cos(q.angle - m.angle) > 0.99
        ):
            q.cx, q.cy, q.angle = m.cx, m.cy, m.angle  # its mate's place, which it is turned from (to the grid)
        return q

    out = _copy.copy(lay)
    out.placed = {i: back(p) for i, p in lay.placed.items()}
    out.ghosts = [back(g) for g in lay.ghosts]
    return out


def _body_touch(p, q) -> float:
    """Page units two bodies of different rosette copies may meet by and still only touch: what the grid can move a
    corner by, a copy's angle on the 1/_GRID rad grid (half the longer body's length times that) and its centre on
    the 1/_GRID page-unit grid (layout._rosette snaps both after the turn)."""
    from .layout import _GRID

    return (max(p.length, q.length) / 2 + 1.0) / _GRID


def _bodies_meet(p, q, touch: float = 1e-6) -> bool:
    """Whether two elements' bodies (rectangles along their own axes, length x width) overlap by more than `touch`,
    by separating axes: a rosette's turned copies stand at any angle, where their rects (the boxes round them) meet
    without the bodies meeting."""

    def corners(e):
        d = np.array([math.cos(e.angle), math.sin(e.angle)])
        n = np.array([-d[1], d[0]])
        c = np.array([e.cx, e.cy])
        hl, hw = e.length / 2, e.width / 2
        return np.array([c + a * hl * d + b * hw * n for a in (-1, 1) for b in (-1, 1)]), (d, n)

    (cp, ap), (cq, aq) = corners(p), corners(q)
    for axis in (*ap, *aq):
        a, b = cp @ axis, cq @ axis
        if min(a.max(), b.max()) - max(a.min(), b.min()) <= touch:
            return False
    return True


def _own_rect(p) -> tuple[float, float, float, float]:
    """A rosette element's rect in its copy's own frame (local_layout's `p`), its angle read on the page axis it was
    laid out along when within one 1/_GRID rad step of it: layout._rosette stores every angle on that grid (π/2 as
    1.5708), which tilts an upright strand's box by a few 1e-6 page units, and two sheet neighbours that touch in the
    protomer's own layout would then meet by more than the default rule's 1e-6."""
    from .layout import _GRID

    quarter = round(p.angle / (math.pi / 2)) * (math.pi / 2)
    a = quarter if abs(p.angle - quarter) <= 1.0 / _GRID else p.angle
    hx = abs(math.cos(a)) * p.length / 2 + abs(math.sin(a)) * p.width / 2
    hy = abs(math.sin(a)) * p.length / 2 + abs(math.cos(a)) * p.width / 2
    return (p.cx - hx, p.cy - hy, p.cx + hx, p.cy + hy)


def _rects_meet(a, b) -> bool:
    """Whether two page boxes overlap by more than 1e-6 both ways (the default figure's element-overlap rule)."""
    return min(a[2], b[2]) - max(a[0], b[0]) > 1e-6 and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6


def element_overlaps(lay) -> list[tuple[str, str]]:
    """Pairs of element or stem boxes (their rects) that overlap, except two consecutive pieces of one bent duplex,
    which meet at their joint (layout.joined_pieces), and a helix drawn over a sheet with that sheet's strands
    (layout.lying_over). nacheck.render_checks counts these as its overlaps. In a rosette (Layout.rosette) two
    elements of one copy overlap when their rects do in the copy's own frame (local_layout, _own_rect: the default
    rule); two
    of different copies when their bodies do (_bodies_meet, by what the grid can move a corner: _body_touch), as
    the turned copies' page rects are boxes round tilted bodies. A bent duplex drawn as a curve (Layout.na_curves) is
    its stretches (_stem_boxes), named once by the curve's name, and overlaps itself also where its band (rails and
    letters) comes over itself (NACurve.overlap): (name, name)."""
    joined = joined_pieces(lay.stems, lay.nucleic) if lay.nucleic is not None else set()
    joined |= lying_over(lay)
    if getattr(lay, "rosette", 0):
        sector, own = lay.frame.sector, local_layout(lay).placed
        out = []
        for p, q in combinations(lay.placed.values(), 2):
            if frozenset((p.sse.id, q.sse.id)) in joined:
                continue
            if sector.get(p.sse.chain) == sector.get(q.sse.chain):
                meet = _rects_meet(_own_rect(own[p.sse.id]), _own_rect(own[q.sse.id]))
            else:
                meet = _bodies_meet(p, q, _body_touch(p, q))
            if meet:
                out.append((p.sse.id, q.sse.id))
        return out
    boxes = [(p.sse.id, p.rect) for p in lay.placed.values()] + _stem_boxes(lay, joined)
    out = [(i, j) for (i, a), (j, b) in combinations(boxes, 2) if _rects_meet(a, b) and frozenset((i, j)) not in joined]
    if not getattr(lay, "na_curves", None):
        return out
    whole = lambda i: i.split("@", 1)[0]  # noqa: E731 - a curve's stretch: the curve, one thing however many meet
    out = list(dict.fromkeys((whole(i), whole(j)) for i, j in out))
    over = [(c.name, c.name) for c in lay.na_curves if c.overlap]  # its band (rails and letters) over itself
    return out + [pair for pair in over if pair not in out]


def _stem_boxes(lay, joined: set) -> list[tuple[str, tuple[float, float, float, float]]]:
    """(id, rect) of each stem body as drawn: every stem; with bent duplexes drawn as curves (Layout.na_curves), each
    curve's straight stretches (layout.stem_bodies) in place of its pieces, named "<curve name>@<arc length>" (counted
    as the curve: element_overlaps). Two stretches of one curve less than two band widths apart along it
    (na_curve._FAR: one stretch of the band, not two that come together) are added to `joined`: they meet where the
    curve turns, as two consecutive pieces at their joint; two further apart overlap as any two stems do (their boxes
    meet), as two pieces of one duplex that are not consecutive do."""
    curves = getattr(lay, "na_curves", None)
    if not curves:
        return [(d.id, d.rect) for d in lay.stems]
    from .layout import stem_bodies

    out, along = [], []
    for d in stem_bodies(lay):
        n = getattr(d, "curve", -1)
        if n < 0:
            out.append((d.id, d.rect))
            continue
        mid = (d.s0 + d.s1) / 2
        name = f"{curves[n].name}@{mid:.2f}"
        out.append((name, d.rect))
        along.append((n, mid, name, d))
    for (n, a, x, d), (m, b, y, e) in combinations(along, 2):
        if n != m or not _rects_meet(d.rect, e.rect):  # (only boxes that meet can count)
            continue
        if abs(a - b) < _FAR * curves[n].band:
            joined.add(frozenset((x, y)))
    return out


# ------------------------------------------------------------------------------------------- text on the page

# Ink is measured under one hinting setting, whatever the caller's rcParams: matplotlib 3.10 and 3.11 hint by
# default differently, and the hinted layout moved text collisions (5 of 9 fixtures) and stem-name gaps (up to
# 0.061 page units) between them (plan 10 Task 1 review).
MEASURE_RC = {"text.hinting": "no_hinting"}
if _MPL < (3, 11):  # 3.11 deprecates the factor and hints at 1 always
    MEASURE_RC["text.hinting_factor"] = 1
INK_FONT = "Liberation Sans"  # the font the ink gates are calibrated on (CI installs fonts-liberation)
INK_CHECKS = ("text_collisions", "off_page", "stem_name_gaps", "domain_name_collisions")  # gated on INK_FONT
# Tolerances later tasks gate with (R7 revised, Task 1 fix round 1). Under MEASURE_RC, with Liberation Sans, the
# same figures measured with matplotlib 3.10.9 (the minimum-deps pins) and 3.11.2 agree exactly on protein figures;
# on nucleic-acid figures (text at other sizes) ink boxes still move by up to 0.15 px, stem-name gaps by up to 0.047
# page units (1KX5; 0.015 elsewhere) and 3 of 14 figures gain or lose one near-touch text collision. So a stem name's
# gap to its letters is gated at stem_gap_tolerance(margin), and a line is counted against a text within
# LINE_TEXT_CLEARANCE_PX of its ink (T6), wider than that spread.
STEM_GAP_FLOOR = 0.08  # page units
LINE_TEXT_CLEARANCE_PX = 2.0


def stem_gap_tolerance(margin: float) -> float:
    """The largest gap (page units) a stem name may keep from its own letters when it should sit `margin` from
    them: the margin plus 1e-3, never under STEM_GAP_FLOOR (what matplotlib versions disagree by)."""
    return max(margin + 1e-3, STEM_GAP_FLOOR)


def _renderer(fig):
    from matplotlib import rc_context
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    canvas = FigureCanvasAgg(fig)
    with rc_context(MEASURE_RC):
        canvas.draw()
    return canvas.get_renderer()


def _texts(fig) -> list:
    return [t for ax in fig.axes for t in ax.texts if t.get_visible() and t.get_text().strip()] + [
        t for t in fig.texts if t.get_visible() and t.get_text().strip()
    ]


def _names(texts) -> list[str]:
    """Each text's gid; a text drawn without one (the title, legend entries) by its place among the texts and its
    words, so two such texts are never taken for one."""
    return [t.get_gid() or f"text#{k}:{t.get_text()!r}" for k, t in enumerate(texts)]


def _ink_boxes(fig, renderer=None) -> list[tuple[object, Bbox]]:
    """(text, ink box in pixels) for every drawn text: the extents of its glyph outlines (_ink_outlines), tighter
    than the layout box get_window_extent gives (line height, descent)."""
    return [(t, Bbox.union([p.get_extents() for p in paths])) for t, paths in _ink_outlines(fig, renderer)]


def _ink_outlines(fig, renderer=None) -> list[tuple[object, list[MPath]]]:
    """(text, its glyph outlines in pixels) for every drawn text: matplotlib's TextToPath, placed and turned as the
    text is drawn, under MEASURE_RC. A text that never reaches the glyph recorder is an error, not a text left
    out."""
    from matplotlib import rc_context
    from matplotlib.textpath import TextToPath

    renderer = renderer or _renderer(fig)
    to_path = TextToPath()
    got: list[Bbox] = []
    calls = [0]

    def record(gc, x, y, s, prop, angle, ismath=False, mtext=None):
        calls[0] += 1
        if not s.strip():
            return
        verts, codes = to_path.get_text_path(prop, s, ismath=ismath)
        if len(verts):
            scale = renderer.points_to_pixels(prop.get_size_in_points()) / to_path.FONT_SCALE
            if renderer.flipy():  # as RendererBase._draw_text_as_path places them: y counted from the top
                y = renderer.get_canvas_width_height()[1] - y
            place = Affine2D().scale(scale).rotate_deg(angle).translate(x, y)
            got.append(MPath(verts, codes).transformed(place))

    renderer.draw_text = record  # this instance only: the glyph outlines are measured, nothing is drawn
    out = []
    texts = _texts(fig)
    try:
        with rc_context(MEASURE_RC):
            for t, name in zip(texts, _names(texts)):
                got.clear()
                calls[0] = 0
                t.draw(renderer)
                if not calls[0]:
                    raise RuntimeError(f"figcheck: text {name} drew nothing to measure")
                if got:
                    out.append((t, list(got)))
    finally:
        del renderer.draw_text
    return out


def text_fonts(fig) -> list[str]:
    """The font files the figure's texts resolve to on this machine, sorted."""
    from matplotlib.font_manager import findfont

    return sorted({findfont(t.get_fontproperties()) for t in _texts(fig)})


def ink_font_reason(fig) -> str:
    """'' when every text resolves to INK_FONT, else why the ink-gated checks are skipped."""
    from matplotlib.ft2font import FT2Font

    other = [f for f in text_fonts(fig) if FT2Font(f).family_name != INK_FONT]
    if not other:
        return ""
    return f"ink checks are calibrated on {INK_FONT}; this machine draws text in {', '.join(other)}"


def _overlap(a: Bbox, b: Bbox) -> float:
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return w * h if w > 0 and h > 0 else 0.0


def _owner(gid: str) -> str | None:
    """The element, ligand, stem ... a text or a body belongs to (one key for both), None for neither."""
    for prefix in ("label:", "ghostlabel:", *ELEMENT_BODIES):
        if gid.startswith(prefix):
            return "el:" + gid[len(prefix) :]
    if gid.startswith("resnum:"):
        return "el:" + gid[len("resnum:") :].rsplit(":", 1)[0]
    for body, text in OWN_BODIES.items():
        for prefix in (body, text):
            if gid.startswith(prefix):
                return body + gid[len(prefix) :]
    return None


class Page:
    """A figure measured once (measure()): its texts with their names and ink boxes, its bodies and lines, all in
    pixels, its legend marks' extents, and the fonts its text resolves to."""

    def __init__(self, drawn):
        fig = drawn.figure
        self.fig = fig
        renderer = _renderer(fig)
        self.outlines = dict((id(t), paths) for t, paths in _ink_outlines(fig, renderer))  # text -> pixels
        self._turned_boxes: dict[tuple[int, float], Bbox] = {}
        boxes = {k: Bbox.union([p.get_extents() for p in paths]) for k, paths in self.outlines.items()}
        texts = _texts(fig)
        self.texts = [(t, n, boxes[id(t)]) for t, n in zip(texts, _names(texts)) if id(t) in boxes]
        self.fonts = text_fonts(fig)
        self.font_reason = ink_font_reason(fig)
        self.ax = fig.axes[0] if fig.axes else None
        self.layout = getattr(drawn, "layout", None)
        artists = [a for ax in fig.axes for a in (*ax.patches, *ax.lines)]
        self.bodies = []  # (gid, path in pixels, zorder, opaque)
        self.lines = []  # (gid, [polyline (n, 2) in pixels], zorder, half its stroke width in pixels)
        self.legend = []  # (gid, extent in pixels) of the legend's marks
        for a in artists:
            gid = a.get_gid() or ""
            if gid.startswith(LINE_GIDS) or GLYCAN_STEM.match(gid):
                self.lines.append((gid, _polylines(a), a.get_zorder(), a.get_linewidth() * fig.dpi / 72 / 2))
            elif gid.startswith((*ELEMENT_BODIES, *OWN_BODIES, *DATA_BODIES, "na:stem:", "na:curve:")):
                path = a.get_path().transformed(a.get_transform())
                self.bodies.append((gid, path, a.get_zorder(), _opaque(a)))
            if gid.startswith("legend"):
                self.legend.append((gid, a.get_window_extent(renderer)))
        self.stems = self._stem_bodies()

    def pair_boxes(self, ta, a: Bbox, tb, b: Bbox) -> tuple[Bbox, Bbox]:
        """Two texts' ink boxes to compare: for two texts turned alike (a stem's name and its letters), the boxes
        of their outlines in that turned frame, which hold their inks tightly; else their page boxes (a and b),
        which hold a turned text's ink loosely (its corners)."""
        turn = ta.get_rotation()
        if abs(turn - tb.get_rotation()) > 1e-6 or turn % 90 == 0:
            return a, b
        return self._turned(ta, turn), self._turned(tb, turn)

    def _turned(self, t, turn: float) -> Bbox:
        key = (id(t), turn)
        if key not in self._turned_boxes:
            back = Affine2D().rotate_deg(-turn)
            self._turned_boxes[key] = Bbox.union([p.transformed(back).get_extents() for p in self.outlines[id(t)]])
        return self._turned_boxes[key]

    def _stem_bodies(self) -> list[tuple[str, MPath, set[str]]]:
        """Each stem's body (its length x width box, turned) in pixels, with the texts that are its own: its name,
        its first pair's numbers and its nucleotides' letters. A bent duplex drawn as a curve (Layout.na_curves) is
        its pieces' place: its rails' width along it in stretches a band long (na:curve:<n>), each owning the curve's
        texts within na_curve._FAR bands of it along the curve only (its name, :first, its letters), as a piece owns
        its own: one gyre's letters on another gyre are flagged (A6), as on another piece."""
        lay = self.layout
        if lay is None or not lay.stems or self.ax is None:
            return []
        from .layout import dna_to_page

        na = lay.nucleic
        curves = getattr(lay, "na_curves", None) or []
        hidden = on_curves(lay)
        out = []
        for k, d in enumerate(lay.stems):
            if k in hidden:
                continue
            ends = ((-1, -1), (1, -1), (1, 1), (-1, 1))
            corners = [dna_to_page(d, x * d.length / 2, y * d.width / 2) for x, y in ends]
            path = MPath(self.ax.transData.transform(np.array([*corners, corners[0]])), closed=True)
            own = {f"na:label:{k}", f"na:stem:{k}:first"}
            if na is not None:
                own |= {f"na:nt:{na.strands[s].labels[i].chain}:{na.strands[s].labels[i].tag}" for s, i in d.axial}
            out.append((f"na:stem:{k}", path, own))

        to_data = self.ax.transData.inverted()
        drawn_at = {gid: box for _, gid, box in self.texts}  # each text's ink box (pixels), by its gid

        def arc_of(c, gid: str) -> float | None:  # where along the curve a text of it stands: nearest its ink
            box = drawn_at.get(gid)
            if box is None:
                return None
            x, y = to_data.transform(((box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2))
            return float(c.s[int(np.argmin(np.hypot(c.xy[:, 0] - x, c.xy[:, 1] - y)))])

        for n, c in enumerate(curves):  # one body a band long at a time, owning the texts within _FAR bands along
            lo, hi = c.extent()
            # its name and :first own the stretches near where each is drawn (a block's name: inside, mid-curve)
            texts = {g: arc_of(c, g) for g in (f"na:curve:{n}:label", f"na:curve:{n}:first")}
            edges = np.linspace(lo, hi, max(1, math.ceil((hi - lo) / max(c.band, 0.1))) + 1)
            for s0, s1 in zip(edges[:-1].tolist(), edges[1:].tolist()):
                ring = np.vstack([c.offset_line(s0, s1, c.width / 2), c.offset_line(s0, s1, -c.width / 2)[::-1]])
                path = MPath(self.ax.transData.transform(np.vstack([ring, ring[:1]])), closed=True)
                reach = (s0 - _FAR * c.band, s1 + _FAR * c.band)
                own = {g for g, at in texts.items() if at is not None and reach[0] <= at <= reach[1]}
                if na is not None:
                    own |= {
                        f"na:nt:{na.strands[s].labels[i].chain}:{na.strands[s].labels[i].tag}"
                        for (s, i), at in c.at.items()
                        if reach[0] <= at <= reach[1]
                    }
                out.append((f"na:curve:{n}", path, own))
        return out

    def own_lines(self) -> dict[str, set[str]]:
        """Text gid -> the gids of the lines it belongs to: a loop's length, a gap's label, a tie's letters, a
        ligand's tethers, a strand's letters and ends to its own drawn runs, a chain end's letter and tail length
        to its stub, a domain's name to its leader."""
        out: dict[str, set[str]] = {}
        for gid, *_ in self.lines:
            if gid.startswith("loop"):
                name = gid.split(":", 1)[1]
                name = name.rsplit(":", 1)[0] if gid.startswith("loop-seg:") else name  # loop-seg:<name>:<k>
                out.setdefault(f"loop-length:{name}", set()).add(gid)
            elif gid.startswith("ligand-tether:"):
                out.setdefault("ligand-label:" + gid[len("ligand-tether:") :].rsplit(":", 1)[0], set()).add(gid)
            elif gid.startswith("na:tie:"):
                out.setdefault(gid.removesuffix(":mark") + ":label", set()).add(gid)
            elif gid.startswith("stub:"):  # stub:<end>:<chain>: its letter and its tail's length
                end_chain = gid[len("stub:") :]
                for text in (f"terminus:{end_chain}", f"tail-length:{end_chain}"):
                    out.setdefault(text, set()).add(gid)
            elif gid.startswith("domain-leader:"):
                out.setdefault("domain-label:" + gid[len("domain-leader:") :], set()).add(gid)
        lay = self.layout
        na = getattr(lay, "nucleic", None)
        for n, p in enumerate(getattr(lay, "na_paths", None) or []):
            if na is None:
                break
            strand = na.strands[p.strand]
            mine = {f"na:gap:{n}", f"na:end:{strand.id}:5′", f"na:end:{strand.id}:3′"}
            mine |= {f"na:nt:{strand.labels[i].chain}:{strand.labels[i].tag}" for i in p.nts}
            for t in mine:
                out.setdefault(t, set()).add(f"na:run:{n}")
        return out


def measure(drawn) -> Page:
    """The figure drawn and measured once, for any number of the page checks (their `page` argument)."""
    return Page(drawn)


def _opaque(a) -> bool:
    """A filled patch, not see-through (a helix drawn over a sheet is): it hides what is drawn under it."""
    if not hasattr(a, "get_fill") or not a.get_fill():
        return False
    alpha = a.get_alpha()
    return a.get_facecolor()[3] >= 0.99 and (alpha is None or alpha >= 0.99)


def _polylines(a) -> list[np.ndarray]:
    """An artist's path in pixels as polylines (curves flattened, broken at NaN)."""
    if hasattr(a, "get_xydata"):  # Line2D
        xy = a.get_transform().transform(a.get_xydata())
        cuts = np.flatnonzero(~np.isfinite(xy).all(axis=1))
        pieces = [p[np.isfinite(p).all(axis=1)] for p in np.split(xy, cuts)]
        return [p for p in pieces if len(p) >= 2]
    path = a.get_path()
    return [p for p in path.to_polygons(a.get_transform(), closed_only=False) if len(p) >= 2]


def _clip(p: np.ndarray, q: np.ndarray, box: Bbox) -> tuple[float, float] | None:
    """The stretch (t0, t1) of the segment p->q inside the box (Liang-Barsky), None when it misses."""
    t0, t1 = 0.0, 1.0
    d = q - p
    for dv, lo, hi, v in ((d[0], box.x0, box.x1, p[0]), (d[1], box.y0, box.y1, p[1])):
        if abs(dv) < 1e-12:
            if v < lo or v > hi:
                return None
            continue
        a, b = (lo - v) / dv, (hi - v) / dv
        t0, t1 = max(t0, min(a, b)), min(t1, max(a, b))
        if t0 > t1:
            return None
    return t0, t1


def _visible(points: np.ndarray, covers: list[MPath]) -> bool:
    """Is any of the points (of a line) hidden under none of `covers` (the opaque bodies drawn above the line)?"""
    hidden = np.zeros(len(points), bool)
    for path in covers:
        hidden |= path.contains_points(points)
        if hidden.all():
            return False
    return not hidden.all()


def _crosses(polys: list[np.ndarray], box: Bbox, covers: list[MPath]) -> bool:
    for poly in polys:
        for p, q in zip(poly, poly[1:]):
            span = _clip(p, q, box)
            if span is None:
                continue
            n = max(2, int(np.hypot(*(q - p)) * (span[1] - span[0]) / _STEP_PX) + 1)
            t = np.linspace(span[0], span[1], n)[:, None]
            if _visible(p + t * (q - p), covers):
                return True
    return False


def _area_on(box: Bbox, path: MPath) -> float:
    """Square pixels of the box that lie on the (closed) path, sampled every _STEP_PX."""
    xs = np.arange(box.x0 + _STEP_PX / 2, box.x1, _STEP_PX)
    ys = np.arange(box.y0 + _STEP_PX / 2, box.y1, _STEP_PX)
    if not len(xs) or not len(ys):
        return 0.0
    gx, gy = np.meshgrid(xs, ys)
    return float(path.contains_points(np.c_[gx.ravel(), gy.ravel()]).sum()) * _STEP_PX**2


def _on(box: Bbox, path: MPath, extent: Bbox | None = None) -> bool:
    """Does the body cover more than TEXT_OVERLAP_PX2 of the text's ink lying deeper than TEXT_DEPTH_PX inside it?
    `extent`: the path's extents, when known (text_collisions measures each body's once)."""
    core = _core(box)
    return (
        core is not None
        and _overlap(core, path.get_extents() if extent is None else extent) > TEXT_OVERLAP_PX2
        and _area_on(core, path) > TEXT_OVERLAP_PX2
    )


def _extents(boxes) -> np.ndarray:
    """(n, 4) x0, y0, x1, y1 rows of Bboxes (each normalised as Bbox.overlaps reads it)."""
    e = np.array([b.extents for b in boxes], float).reshape(-1, 4)
    return np.hstack([np.minimum(e[:, :2], e[:, 2:]), np.maximum(e[:, :2], e[:, 2:])])


def _meeting(box, rows: np.ndarray, open_: bool) -> np.ndarray:
    """Indices of the rows (_extents) that `box` (x0, y0, x1, y1) meets: sharing area (open_, as _overlap > 0
    needs) or touching (as Bbox.overlaps counts it). A prefilter only: what it leaves out could not pass the test
    it stands before, so the checks give what they gave one pair at a time (random-100 C33: 7JU4's text checks
    took 90 s)."""
    x0, y0, x1, y1 = box
    if open_:
        m = (rows[:, 0] < x1) & (x0 < rows[:, 2]) & (rows[:, 1] < y1) & (y0 < rows[:, 3])
    else:
        m = (rows[:, 0] <= x1) & (x0 <= rows[:, 2]) & (rows[:, 1] <= y1) & (y0 <= rows[:, 3])
    return np.flatnonzero(m)


def _depth(*boxes: Bbox) -> float:
    """How deep an overlap must reach into the texts' inks to count: TEXT_DEPTH_PX, or a third of the smallest ink
    side for text smaller than that (its whole ink covered still counts)."""
    return min(TEXT_DEPTH_PX, *(min(b.width, b.height) / 3 for b in boxes))


def _core(box: Bbox) -> Bbox | None:
    """A text's ink box less _depth all round (None when nothing is left)."""
    d = _depth(box)
    if box.width <= 2 * d or box.height <= 2 * d:
        return None
    return Bbox.from_extents(box.x0 + d, box.y0 + d, box.x1 - d, box.y1 - d)


def _texts_meet(a: Bbox, b: Bbox) -> bool:
    """Two texts' inks overlap by more than TEXT_OVERLAP_PX2, deeper than _depth both ways."""
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    d = _depth(a, b)
    return w > d and h > d and w * h > TEXT_OVERLAP_PX2


def text_collisions(drawn, *, clearance_px: float = 0.0, page: Page | None = None) -> list[str]:
    """Every "<text> × <thing gid>" a reader sees on the page (a text is named by its gid, one drawn without a gid
    by its place and words: _names):
    - text over text, their ink boxes overlapping by more than TEXT_OVERLAP_PX2 and deeper than TEXT_DEPTH_PX both
      ways (a text's own parts, one gid, never); two texts turned alike are compared in their turned frame
      (Page.pair_boxes: a stem's name beside its letter row);
    - text over a body (an element, a stem, a ligand or linked mark, a data category mark) it does not belong to (its
      label, its residue numbers, its letters, its name) or a helix's or 3-10 helix's bar under its own name
      (BESIDE_BODIES: the name stands beside its bar), by more than TEXT_OVERLAP_PX2 of its ink deeper than
      TEXT_DEPTH_PX inside it, once per element (a helix is drawn as helix: and helix-back:); a helix drawn over a
      sheet and that sheet's strands are each other's (layout.lying_over);
    - text crossed by a line (LINE_GIDS: loops, nucleic-acid runs and ties, tethers, disulfide bars, glycan stems,
      H-bonds and the other interaction lines) whose stroke comes within `clearance_px` of its ink box where no
      opaque body drawn above the line hides it there: a line under the text still runs through its letters, a
      loop under the strand its letter sits on does not show. A text's own line (a loop's length, a gap's label,
      a tie's letter, a ligand's tethers, a strand's letters and ends) is not counted."""
    page = page or measure(drawn)
    lay = page.layout
    over = {frozenset(f"el:{k}" for k in pair) for pair in lying_over(lay)} if lay is not None else set()
    own_lines = page.own_lines()
    hits: list[str] = []
    lines = [  # each line, the extent of its drawing, its reach and the opaque bodies drawn above it
        (gid, polys, Bbox.union([Bbox.from_extents(*p.min(axis=0), *p.max(axis=0)) for p in polys]), pad, covers)
        for gid, polys, z, pad in page.lines
        if polys
        for covers in [[p for _, p, bz, opaque in page.bodies if opaque and bz > z]]
    ]
    inks = _extents([box for *_, box in page.texts])
    for i, (ta, ga, a) in enumerate(page.texts):  # pairs in combinations' order; those whose inks share no area skipped
        for j in _meeting(inks[i], inks[i + 1 :], True).tolist():
            tb, gb, b = page.texts[i + 1 + j]
            if ga != gb and _texts_meet(a, b) and _texts_meet(*page.pair_boxes(ta, a, tb, b)):
                hits.append(f"{ga} × {gb}")
    stems = ("na:stem:", "na:curve:")  # stems and curves: their whole body below (page.stems)
    bodies = [(body, path, _owner(body) or body) for body, path, _, _ in page.bodies if not body.startswith(stems)]
    body_rows = _extents([path.get_extents() for _, path, _ in bodies])  # each body measured once
    stem_rows = _extents([path.get_extents() for _, path, _ in page.stems])
    line_rows = _extents([extent for _, _, extent, _, _ in lines])
    pads = np.array([pad for _, _, _, pad, _ in lines], float)
    for t, name, box in page.texts:
        gid = t.get_gid() or ""
        mine = _owner(gid)
        named = gid.startswith("label:")  # an element's name: a helix's stands beside its bar (BESIDE_BODIES)
        met: set[str] = set()
        core = _core(box)
        if core is not None:  # stems: their whole body below
            for k in _meeting(core.extents, body_rows, True).tolist():
                body, path, theirs = bodies[k]
                owned = theirs == mine and not (named and body.startswith(BESIDE_BODIES))  # a bar: not under its name
                if owned or theirs in met or frozenset((theirs, mine)) in over:
                    continue
                if _on(box, path, Bbox(body_rows[k].reshape(2, 2))):
                    met.add(theirs)
                    hits.append(f"{name} × {body.replace('helix-back:', 'helix:', 1)}")
            stems_met: set[str] = set()  # a curve is several stretches (Page._stem_bodies): met once
            for k in _meeting(core.extents, stem_rows, True).tolist():
                body, path, own = page.stems[k]
                if body not in stems_met and gid not in own and _on(box, path, Bbox(stem_rows[k].reshape(2, 2))):
                    stems_met.add(body)
                    hits.append(f"{name} × {body}")
        r = pads + clearance_px  # each line's stroke meets the ink, or comes within the clearance
        near = (line_rows[:, 0] <= box.x1 + r) & (box.x0 - r <= line_rows[:, 2])
        near &= (line_rows[:, 1] <= box.y1 + r) & (box.y0 - r <= line_rows[:, 3])
        for k in np.flatnonzero(near).tolist():
            line, polys, extent, pad, covers = lines[k]
            r = pad + clearance_px
            reach = Bbox.from_extents(box.x0 - r, box.y0 - r, box.x1 + r, box.y1 + r)
            if line in own_lines.get(gid, ()) or not extent.overlaps(reach):
                continue
            if _crosses(polys, reach, covers):
                hits.append(f"{name} × {line}")
    return hits


def off_page(drawn, *, page: Page | None = None) -> list[str]:
    """Texts (their ink) and legend marks reaching more than OFF_PAGE_PX past the page edge."""
    page = page or measure(drawn)
    edge = page.fig.bbox
    out = [name for _, name, box in page.texts if _beyond(box, edge)]
    return out + [gid for gid, box in page.legend if _beyond(box, edge)]


def _beyond(box: Bbox, page: Bbox) -> bool:
    return (
        box.x0 < page.x0 - OFF_PAGE_PX
        or box.y0 < page.y0 - OFF_PAGE_PX
        or box.x1 > page.x1 + OFF_PAGE_PX
        or box.y1 > page.y1 + OFF_PAGE_PX
    )


# ------------------------------------------------------------------------------------------- labels and copies


def label_repeats(lay) -> dict[str, list[str]]:
    """Chain -> the labels two or more of its drawn elements carry; "<na>" -> stem names (a curve's, not its pieces')
    given to two or more stems (ghosts aside)."""
    count: dict[str, Counter] = {}
    for p in lay.placed.values():
        count.setdefault(p.sse.chain, Counter())[p.label] += 1
    repeated = {c: sorted(k for k, n in got.items() if n > 1) for c, got in count.items()}
    hidden = on_curves(lay)  # drawn as curves: their names
    names = Counter(d.name for k, d in enumerate(lay.stems) if not d.ghost and d.name and k not in hidden)
    names.update(c.name for c in getattr(lay, "na_curves", None) or [] if c.name)
    repeated["<na>"] = sorted(k for k, n in names.items() if n > 1)
    return {c: ks for c, ks in repeated.items() if ks}


def fit_deviation(a, b, mirrored: bool = False) -> float:
    """The largest distance between matched points a[k] and b[k] after the best rigid motion of b onto a (least
    squares: a rotation and a translation; a reflection and a translation when `mirrored`)."""
    return float(fit_deviations(a, b, mirrored).max())


def fit_deviations(a, b, mirrored: bool = False) -> np.ndarray:
    """The distance between each pair of matched points a[k] and b[k] after the best rigid motion of b onto a
    (fit_deviation's: A15's, and H2's for a figure drawn like a reference, criteria.h2_places)."""
    a, b = np.asarray(a, float).reshape(-1, 2), np.asarray(b, float).reshape(-1, 2)
    if mirrored:
        b = b * [-1.0, 1.0]
    turn, cb, ca = rigid_fit(a, b)
    fitted = (b - cb) @ turn + ca
    return np.hypot(*(fitted - a).T)


def rigid_fit(a, b) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(turn, cb, ca): the best rigid motion (a turn and a shift, no mirror) of the 2D points b onto their matched
    points a, a point q of b's moved to (q - cb) @ turn + ca (fit_deviations' Procrustes)."""
    a, b = np.asarray(a, float).reshape(-1, 2), np.asarray(b, float).reshape(-1, 2)
    ca, cb = a.mean(axis=0), b.mean(axis=0)
    u, _, vt = np.linalg.svd((b - cb).T @ (a - ca))
    turn = u @ vt
    if np.linalg.det(turn) < 0:  # the best proper rotation
        u[:, -1] *= -1
        turn = u @ vt
    return turn, cb, ca


COPY_MIN_MATCHED = 3  # a tied copy is compared on at least this many elements, and on half the smaller copy's


def shared_sheet_elements(lay, copy: str, ref: str) -> set[str]:
    """The strands of every sheet block (Layout.sheets) holding strands of both a copy and its reference (a dimer's
    continuous sheet, a swapped strand) that the layout left untied: one item holding both copies, placed by the
    sheet, so its strands cannot sit at their copy's tied place: copies_alike leaves them out (C10 ruling). A sheet a
    symmetry tie holds (Layout.tied_ids: a ring's sheets, with a strand lent by a third copy or not) stays in, and so
    is measured (8UTF b~a, 8UUP D~B: fix round 1)."""
    out: set[str] = set()
    for block in lay.sheets:
        if {copy, ref} <= {lay.placed[k].sse.chain for k in block if k in lay.placed} and not set(block) & lay.tied_ids:
            out |= set(block)
    return out


COPY_IDENTITY = 0.9  # two chains are copies of one protein when this fraction of their aligned residues match ...
COPY_LENGTH = 0.8  # ... and the shorter is at least this fraction of the longer
COPY_ALIGNED = 20  # ... and at least this many do (as model.homolog_chains asks of residues at the same numbers)


def aligned_positions(xs, ys) -> list[tuple[int, int]]:
    """(i, j) for every residue xs[i] the global alignment of two residue-name sequences puts against ys[j]
    (gemmi.align_string_sequences; equal sequences align one to one)."""
    xs, ys = list(xs), list(ys)
    if xs == ys:
        return [(i, i) for i in range(len(xs))]
    import gemmi

    cigar = gemmi.align_string_sequences(xs, ys, [], gemmi.AlignmentScoring()).cigar_str()
    out, i, j = [], 0, 0
    for n, op in re.findall(r"(\d+)([MID])", cigar):
        n = int(n)
        if op == "M":
            out += [(i + k, j + k) for k in range(n)]
        i += n if op in "MI" else 0
        j += n if op in "MD" else 0
    return out


def same_protein(xs, ys, pairs: list[tuple[int, int]]) -> bool:
    """Are two chains (residue names; pairs from aligned_positions) copies of one protein?"""
    same = sum(xs[i] == ys[j] for i, j in pairs)
    longer = max(len(xs), len(ys))  # a short chain scattered over a long one is no copy of it
    return same >= COPY_ALIGNED and same >= COPY_IDENTITY * len(pairs) and min(len(xs), len(ys)) >= COPY_LENGTH * longer


def _untied_copies(lay, residues: dict[str, list[str]], chains: list[str]) -> list[tuple[str, str]]:
    """(copy, reference) for each drawn chain no tie holds that is a copy of one protein (by sequence) with a tied
    reference or an earlier untied chain, when the figure ties any copies (a symmetric assembly): such a copy
    must be drawn alike too (C18 ruling)."""
    if not lay.copy_ties:
        return []
    tied = set(lay.copy_ties) | {ref for ref, _ in lay.copy_ties.values()}
    roots = sorted({ref for ref, _ in lay.copy_ties.values()}, key=chains.index)
    out = []
    for c in (c for c in chains if c not in tied):
        root = next((r for r in roots if same_protein(residues[r], residues[c], _aligned(residues, r, c))), None)
        if root is None:
            roots.append(c)
        else:
            out.append((c, root))
    return out


def copies_untied(drawn) -> list[str]:
    """The "copy~reference" of each drawn chain that is a copy of one protein (by sequence, same_protein) with an
    earlier drawn chain, when the figure ties no copies (Layout.copy_ties empty): copies drawn with no symmetry tie,
    which A15 reports as not judged rather than as no copies (random-100 held-out 5K7Z: two dimers on a duplex frame,
    untied, said "no copies to compare"). Chains lent by neighbours (Layout.partial: a focus view's fragments) are
    left out; [] when the figure ties any copy (_untied_copies judges those)."""
    lay = drawn.layout
    if lay.copy_ties:
        return []
    residues: dict[str, list[str]] = {}
    first: dict[str, int] = {}
    for k, label in enumerate(drawn.bb.labels):
        residues.setdefault(label.chain, []).append(label.name)
        first.setdefault(label.chain, k)
    partial = getattr(lay, "partial", set()) or set()
    chains = sorted(
        {p.sse.chain for p in lay.placed.values() if not p.ghost and p.sse.chain in residues} - set(partial),
        key=first.__getitem__,
    )
    roots: list[str] = []
    out = []
    for c in chains:
        root = next((r for r in roots if same_protein(residues[r], residues[c], _aligned(residues, r, c))), None)
        if root is None:
            roots.append(c)
        else:
            out.append(f"{c}~{root}")
    return out


def _aligned(residues: dict[str, list[str]], a: str, b: str) -> list[tuple[int, int]]:
    return aligned_positions(residues[a], residues[b])


def _copy_pairs(drawn, notes: dict | None = None) -> dict[tuple[str, str], tuple[float | None, int, int, str]]:
    """(copy, reference) -> (deviation, elements matched, elements of the smaller of the two, why not judged) for each
    tie of Layout.copy_ties, and for each copy no tie holds (_untied_copies). The copies' elements pair one-to-one
    by sequence (C18 ruling): the copy's residues aligned to the reference's (aligned_positions), its element ranges
    taken to the reference's residues, then model.pair_places (the same kind, sharing at least half of the shorter).
    Strands of a sheet shared between the two copies are left out (shared_sheet_elements), and so are the strands
    of sheets the two copies group differently (_sheets_differ, C2b-3 ruling R2; a note for each such tie goes in
    `notes`). The deviation is that of fit_deviation (a reflection for a mirror tie; for an untied copy only when
    the figure mirrors a tie). When fewer than COPY_MIN_MATCHED or than half the smaller copy's elements match (two
    points always fit), the tie is not judged (deviation None) if the smaller copy has fewer than COPY_MIN_MATCHED
    elements to match ("n/a: k elements") or a shared or differing sheet took the rest ("n/a: shared sheet",
    "n/a: sheets differ"), C22 and C10 rulings; else it is inf."""
    from .model import overlap, pair_places

    lay = drawn.layout
    differ_note = {} if notes is None else notes
    residues: dict[str, list[str]] = {}
    first: dict[str, int] = {}
    for k, label in enumerate(drawn.bb.labels):
        residues.setdefault(label.chain, []).append(label.name)
        first.setdefault(label.chain, k)
    mine: dict[str, list[tuple[str, object]]] = {}
    for k, p in lay.placed.items():
        if not p.ghost and p.sse.chain in residues:
            mine.setdefault(p.sse.chain, []).append((k, p))
    chains = sorted(mine, key=first.__getitem__)
    ties = [(copy, ref, (mirrored,)) for copy, (ref, mirrored) in sorted(lay.copy_ties.items())]
    hands = (False, True) if any(m for _, m in lay.copy_ties.values()) else (False,)
    ties += [(copy, ref, hands) for copy, ref in _untied_copies(lay, residues, chains)]
    tie_parts = []
    for copy, ref, mirrors in ties:
        shared = shared_sheet_elements(lay, copy, ref)
        a = [p for k, p in mine.get(ref, []) if k not in shared]
        b = [p for k, p in mine.get(copy, []) if k not in shared]
        to_ref = {j: i for i, j in _aligned(residues, ref, copy)} if ref in residues and copy in residues else {}

        def place(p, chain, mapping=None):
            local = range(p.sse.start - first[chain], p.sse.end - first[chain] + 1)
            got = [mapping[j] for j in local if j in mapping] if mapping is not None else list(local)
            return (p.sse.kind, min(got), max(got)) if got else (p.sse.kind, -2, -3)

        xs, ys = [place(p, ref) for p in a], [place(p, copy, to_ref) for p in b]
        found = sorted(pair_places(xs, ys).items())
        judgeable = len(found) >= COPY_MIN_MATCHED and 2 * len(found) >= min(len(a), len(b))

        def lacking(xs, ys, els):  # strands with no element of the other copy at their place, of any kind
            return {
                e.sse.id
                for e, x in zip(els, xs)
                if e.sse.kind == "E" and not any(overlap(x[1:], y[1:]) > 0 for y in ys)
            }

        differ = _sheets_differ(
            lay, copy, ref, [(a[i], b[j]) for i, j in found], lacking(xs, ys, a) | lacking(ys, xs, b)
        )
        tie_parts.append((copy, ref, mirrors, shared, a, b, xs, ys, found, judgeable, differ))
    # a copy's block another of its chains finds grouped apart from its mate's (9P1X O's 9-strand barrel holding R's
    # strands against E's 8 with J's beside: O~E and R~J differ, so P's strands in it differ too: C2b-S)
    apart = {n for *_, differ in tie_parts for _, mine, _ in differ for n in mine}
    where = {k: n for n, (ids, _) in enumerate(_blocks(lay)) for k in ids}
    out = {}
    for copy, ref, mirrors, shared, a, b, xs, ys, found, judgeable, differ in tie_parts:
        own_groups = {n for _, mine, _ in differ for n in mine}
        through = sorted(
            {
                (where[b[j].sse.id], where.get(a[i].sse.id, -1))
                for i, j in found
                if where.get(b[j].sse.id) in apart - own_groups
            }
        )
        if through:
            differ = differ + [
                (
                    {a[i].sse.id for i, j in found if where.get(b[j].sse.id) == c}
                    | {k for k in _blocks(lay)[c][0] if k in lay.placed and lay.placed[k].sse.chain == copy},
                    [c],
                    sorted({r for cc, r in through if cc == c and r >= 0}),
                )
                for c in sorted({c for c, _ in through})
            ]
        if differ:  # their strands are placed by each copy's own pairing: left out, reported (C2b-3 R2)
            gone = {k for strands, _, _ in differ for k in strands}
            kept_a, kept_b = (
                [n for n, p in enumerate(a) if p.sse.id not in gone],
                [n for n, p in enumerate(b) if p.sse.id not in gone],
            )
            a, b = [a[n] for n in kept_a], [b[n] for n in kept_b]
            found = sorted(pair_places([xs[n] for n in kept_a], [ys[n] for n in kept_b]).items())
            differ_note[(copy, ref)] = (
                "copies' sheets differ: "
                + "; ".join(
                    f"{copy} sheet{'s' if len(c) > 1 else ''} {'+'.join(map(str, c))} vs "
                    f"{ref} sheet{'s' if len(r) > 1 else ''} {'+'.join(map(str, r))}"
                    for _, c, r in differ
                )
                + f"; {len(gone)} strands not judged"
            )
        smaller = min(len(a), len(b))
        why = ""
        if len(found) >= COPY_MIN_MATCHED and 2 * len(found) >= smaller:
            pa = [(a[i].cx, a[i].cy) for i, _ in found]
            pb = [(b[j].cx, b[j].cy) for _, j in found]
            dev: float | None = min(fit_deviation(pa, pb, m) for m in mirrors)
        elif differ and judgeable:  # R2 left too few to fit: its own reason, never C22's, never inf
            dev, why = None, f"n/a: sheets differ ({len(found)} matched)"
        elif 0 < smaller < COPY_MIN_MATCHED:
            dev, why = None, f"n/a: {smaller} element{'s' if smaller > 1 else ''}"
        elif shared and smaller == 0:
            dev, why = None, "n/a: shared sheet"
        elif shared and len(found) == COPY_MIN_MATCHED - 1 and 2 * len(found) >= smaller:
            dev, why = None, f"n/a: {len(found)} matched"  # a shared sheet left two to match, too few to fit (6YQD)
        else:
            dev = math.inf
        out[(copy, ref)] = (dev, len(found), smaller, why)
    return out


def _blocks(lay) -> list[tuple[list[str], list[str]]]:
    """Each drawn sheet block (Layout.sheets: its strands) with its ghosts (Layout.sheet_ghosts)."""
    ghosts = list(getattr(lay, "sheet_ghosts", []) or [])
    return [(list(ids), ghosts[n] if n < len(ghosts) else []) for n, ids in enumerate(lay.sheets)]


def _sheets_differ(
    lay, copy: str, ref: str, pairs, alone: set[str] = frozenset()
) -> list[tuple[set[str], list[int], list[int]]]:
    """The strands (ids, both copies) of the sheet blocks two tied copies group differently (C2b-3 ruling R2): the
    paired strands (`pairs`: (reference element, copy element)) of a block of one copy lie in more than one block of
    the other, a block no tie draws as its mate's (Layout.tied_ids) holds a strand the other copy lacks (`alone`: no
    element of the other copy at its place; 7QUP 7D's sheet, which keeps its own shape), or the two blocks close
    differently (one draws a barrel's ghost of its own strand, the other not; a strand lent across a ring's cut
    leaves a ghost that is no closure). Each copy's strands are placed by its own pairing (S1-S3), so these cannot
    sit alike: A15 leaves them out, as C10 leaves out a shared sheet, and reports them (copies_sheets_differ). For
    each group of blocks so joined: the strands of `copy` and `ref` in it, and its copy's and reference's block
    numbers (Layout.sheets indices, as A11 numbers sheets)."""
    blocks = _blocks(lay)
    where = {k: n for n, (ids, _) in enumerate(blocks) for k in ids}
    links = sorted(
        {(where[pa.sse.id], where[pb.sse.id]) for pa, pb in pairs if pa.sse.id in where and pb.sse.id in where}
    )
    if not links:
        return []
    comp = {}  # union-find over ("r", n) and ("c", n)

    def root(x):
        while comp.setdefault(x, x) != x:
            x = comp[x]
        return x

    for r, c in links:
        comp[root(("r", r))] = root(("c", c))
    groups: dict = {}
    for r, c in links:
        groups.setdefault(root(("r", r)), set()).update({("r", r), ("c", c)})
    out = []
    for group in (groups[k] for k in sorted(groups, key=str)):
        refs = sorted(n for side, n in group if side == "r")
        mine = sorted(n for side, n in group if side == "c")
        closed = {any(ghost in blocks[n][0] for ghost in blocks[n][1]) for n in refs + mine}
        tied = any(k in lay.tied_ids for n in mine for k in blocks[n][0])  # drawn as its mate's, by the tie
        lacks = not tied and any(k in alone for n in refs + mine for k in blocks[n][0])
        if len(refs) == 1 and len(mine) == 1 and len(closed) == 1 and not lacks:
            continue
        strands = {
            k for n in refs + mine for k in blocks[n][0] if k in lay.placed and lay.placed[k].sse.chain in (copy, ref)
        }
        out.append((strands, mine, refs))
    return out


def copies_sheets_differ(drawn) -> dict[tuple[str, str], str]:
    """(copy chain, reference chain) -> a note naming the sheet blocks the two copies group differently, whose
    strands A15 leaves out (reported, not judged: C2b-3 ruling R2)."""
    notes: dict[tuple[str, str], str] = {}
    _copy_pairs(drawn, notes)
    return notes


def copies_alike(drawn) -> dict[tuple[str, str], float]:
    """(copy chain, reference chain) -> the largest element-centre deviation (page units) between the two after the
    best rigid 2D motion, for every tie the symmetry asked for (Layout.copy_ties) and judged; inf when too few
    elements match to say (_copy_pairs). Elements one copy lacks are left out; ties not judged are in
    copies_not_judged."""
    return {pair: dev for pair, (dev, *_) in _copy_pairs(drawn).items() if dev is not None}


def copies_not_judged(drawn) -> dict[tuple[str, str], str]:
    """(copy chain, reference chain) -> why, for every tie too small to judge (reported, not failed: _copy_pairs)."""
    return {pair: why for pair, (dev, _, _, why) in _copy_pairs(drawn).items() if dev is None}


def copy_matches(drawn) -> dict[tuple[str, str], tuple[int, int]]:
    """(copy chain, reference chain) -> (elements matched, elements of the smaller copy) for every tie."""
    return {pair: (n, smaller) for pair, (_, n, smaller, _) in _copy_pairs(drawn).items()}


# ------------------------------------------------------------------------------------------- domains


def outline_overlap(a: MPath, b: MPath, grid: float = PANEL_GRID) -> float:
    """The area two closed outlines (page units) share, over the smaller one's area, both rasterised at `grid`."""
    return _cells_overlap(_cells(a, grid), _cells(b, grid))


def _cells(path: MPath, grid: float = PANEL_GRID) -> np.ndarray:
    """The cells of the shared raster (`grid`, anchored at the origin) whose centres lie inside the closed path (by
    the nonzero rule, as the panel is filled: its holes are outside), each as one sorted integer key."""
    box = path.get_extents()
    i0, j0 = math.floor(box.x0 / grid), math.floor(box.y0 / grid)
    i1, j1 = math.ceil(box.x1 / grid), math.ceil(box.y1 / grid)
    if i1 <= i0 or j1 <= j0:
        return np.empty(0, np.int64)
    gi, gj = np.meshgrid(np.arange(i0, i1), np.arange(j0, j1))
    gi, gj = gi.ravel(), gj.ravel()
    centres = np.c_[(gi + 0.5) * grid, (gj + 0.5) * grid]
    winding = np.zeros(len(centres), int)  # the nonzero rule Agg fills by: a ring turning the other way is a hole
    for ring in path.to_polygons(closed_only=True):  # (Path.contains_points tests each ring alone: holes filled)
        turn = np.sum(ring[:-1, 0] * ring[1:, 1] - ring[1:, 0] * ring[:-1, 1])
        winding += np.where(MPath(ring).contains_points(centres), 1 if turn > 0 else -1, 0)
    hit = winding != 0
    return np.sort((gi[hit].astype(np.int64) << 32) + (gj[hit].astype(np.int64) + (1 << 31)))


def _cells_overlap(a: np.ndarray, b: np.ndarray) -> float:
    if not len(a) or not len(b):
        return 0.0
    return len(np.intersect1d(a, b, assume_unique=True)) / min(len(a), len(b))


def domain_checks(drawn, collisions: list[str] | None = None) -> dict:
    """Named domains on the page:
    cross_chain: domains whose elements lie in more than one chain;
    unalike: elements that are one element in copies of a protein (model.copy_classes) but sit in differently named
        domains ("<element> <domain> vs <element> <domain>", names without their " · chain");
    small: automatic domains (named D1, D2 ...) with fewer than features.DOMAIN_MIN_SSE_RESIDUES residues in their
        elements;
    panel_overlaps: (domain, domain, fraction) for every two panels (the drawn domain:<name> outlines) that share any
        area, the fraction of the smaller panel's area (outline_overlap);
    name_collisions: text_collisions (`collisions`, when already measured) involving a domain name (domain-label:)."""
    from .features import DOMAIN_MIN_SSE_RESIDUES
    from .model import copy_classes, homolog_chains

    lay, fig = drawn.layout, drawn.figure
    domains = [(name, [k for k in ids if k in lay.placed]) for name, ids in lay.domains]
    cross = [name for name, ids in domains if len({lay.placed[k].sse.chain for k in ids}) > 1]
    small = [
        name
        for name, ids in domains
        if _AUTO_DOMAIN.fullmatch(name) and sum(len(lay.placed[k].sse) for k in ids) < DOMAIN_MIN_SSE_RESIDUES
    ]
    mate = copy_classes(drawn.sses, homolog_chains(drawn.bb)).mate
    seen: dict[str, tuple[str, str]] = {}  # mate -> (first element of it in a domain, its domain)
    unalike = []
    for name, ids in domains:
        base = name.split(" · ")[0]
        for k in ids:
            first = seen.setdefault(mate.get(k, k), (k, base))
            if first[1] != base:
                unalike.append(f"{first[0]} {first[1]} vs {k} {base}")
    outlines = {}
    if fig.axes:
        to_page = fig.axes[0].transData.inverted()
        for a in fig.axes[0].patches:
            gid = a.get_gid() or ""
            if gid.startswith("domain:"):
                outlines[gid[len("domain:") :]] = a.get_path().transformed(a.get_transform() + to_page)
    overlaps = []
    cells = {name: _cells(path) for name, path in outlines.items()}
    for (na, ca), (nb, cb) in combinations(cells.items(), 2):
        f = _cells_overlap(ca, cb)
        if f > 0:
            overlaps.append((na, nb, round(f, 4)))
    names = [h for h in (text_collisions(drawn) if collisions is None else collisions) if "domain-label:" in h]
    return {
        "cross_chain": cross,
        "unalike": unalike,
        "small": small,
        "panel_overlaps": overlaps,
        "name_collisions": names,
    }


# ------------------------------------------------------------------------------------------- stem names


def _box_gap(a, b) -> float:
    """Distance between two boxes (x0, y0, x1, y1); 0 when they meet."""
    dx = max(a[0] - b[2], b[0] - a[2], 0.0)
    dy = max(a[1] - b[3], b[1] - a[3], 0.0)
    return float(np.hypot(dx, dy))


def stem_name_gaps(drawn, *, page: Page | None = None) -> dict[str, tuple[float, float]]:
    """na:label:<k> -> (gap from the name's ink to stem k's own letters (na:nt: of its nucleotides), or to the stem's
    rect when its letters are not drawn; gap from the name's ink box to the nearest other stem's rect), page units.
    The gap to the letters is measured in the stem's own frame (along it and across it, as names and letters are
    written along it): across the stem from its letter row on the name's side to the name, where the name lies
    beside that row (a name centred between two letters sits that far from the row, not from either letter); else
    the gap from the name to the nearest letter. A bent duplex drawn as a curve (Layout.na_curves), na:curve:<n>:label:
    in the curve's frame where the name stands (its point nearest the name, along its direction there), against its
    letters within the name's reach plus a band along it; the other stems' bodies as drawn (layout.stem_bodies), its
    own stretches aside."""
    lay, fig = drawn.layout, drawn.figure
    if not lay.stems or not fig.axes:
        return {}
    page = page or measure(drawn)
    to_page = fig.axes[0].transData.inverted()
    texts = {name: t for t, name, _ in page.texts}
    ink = {name: tuple(to_page.transform_bbox(box).extents) for _, name, box in page.texts}
    na = lay.nucleic
    curves = getattr(lay, "na_curves", None) or []
    hidden = on_curves(lay)
    from .layout import stem_bodies

    bodies = stem_bodies(lay)
    out = {}
    for k, d in enumerate(lay.stems):
        name = ink.get(f"na:label:{k}")
        if name is None or k in hidden:
            continue
        tags = [f"na:nt:{na.strands[s].labels[i].chain}:{na.strands[s].labels[i].tag}" for s, i in d.axial]
        frame = Affine2D().translate(-d.cx, -d.cy).rotate(-d.angle)  # page -> (along, across) the stem
        own = _own_gap(page, texts, to_page, f"na:label:{k}", tags, frame, name, [d.rect])
        others = [_box_gap(name, e.rect) for e in bodies if e is not d]  # (stem_bodies: layout.stems itself)
        out[f"na:label:{k}"] = (float(own), min(others) if others else math.inf)
    for n, c in enumerate(curves):
        label = f"na:curve:{n}:label"
        name = ink.get(label)
        if name is None:
            continue
        q = np.array([(name[0] + name[2]) / 2, (name[1] + name[3]) / 2])
        i = int(np.argmin(np.hypot(c.xy[:, 0] - q[0], c.xy[:, 1] - q[1])))
        reach = math.hypot(name[2] - name[0], name[3] - name[1]) / 2 + c.band
        keys = [key for key in c.at if abs(c.at[key] - float(c.s[i])) <= reach]
        tags = [f"na:nt:{na.strands[s].labels[j].chain}:{na.strands[s].labels[j].tag}" for s, j in keys]
        turn = math.atan2(float(c.tangent[i, 1]), float(c.tangent[i, 0]))
        frame = Affine2D().translate(-float(c.xy[i, 0]), -float(c.xy[i, 1])).rotate(-turn)
        at = float(c.s[i])

        def far(e, at=at, reach=_FAR * c.band) -> bool:  # another stretch of the band, far along it
            return max(e.s0 - at, at - e.s1, 0.0) > reach

        mine = [e.rect for e in bodies if getattr(e, "curve", -1) == n and not far(e)]
        own = _own_gap(page, texts, to_page, label, tags, frame, name, mine)
        # every other stem's body, and the curve's own stretches more than _FAR bands along it (a piece's siblings)
        others = [_box_gap(name, e.rect) for e in bodies if getattr(e, "curve", -1) != n or far(e)]
        out[label] = (float(own), min(others) if others else math.inf)
    return out


def _own_gap(page: Page, texts: dict, to_page, label: str, tags: list[str], frame, name, rects) -> float:
    """stem_name_gaps' gap from a stem's name (`label`, its ink box `name` in page units) to its own letters (`tags`)
    in the stem's frame (`frame`: page -> along, across), or to its body's rects when no letter is drawn."""
    local = {
        g: Bbox.union([p.transformed(frame).get_extents() for p in _paths_on_page(page, texts[g], to_page)])
        for g in [label, *tags]
        if g in texts
    }
    letters = [local[g] for g in tags if g in local]
    if not letters:
        return min(_box_gap(name, r) for r in rects)
    box = local[label]
    side = 1.0 if box.y0 + box.y1 > 0 else -1.0
    row = [b for b in letters if (b.y0 + b.y1) * side > 0]
    across = (box.y0 - max(b.y1 for b in row) if side > 0 else min(b.y0 for b in row) - box.y1) if row else -math.inf
    beside = row and box.x1 > min(b.x0 for b in row) and box.x0 < max(b.x1 for b in row)
    return across if beside and across > 0 else min(_box_gap(box.extents, b.extents) for b in letters)


def _segment_box_gap(p, q, box) -> float:
    """Distance between the segment p-q and the box (x0, y0, x1, y1); 0 when they meet."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    if _clip(p, q, Bbox.from_extents(*box)) is not None:
        return 0.0
    x0, y0, x1, y1 = box
    ends = min(_box_gap(box, (*v, *v)) for v in (p, q))  # an end nearest the box
    d = q - p
    n = float(d @ d)
    if n < 1e-24:
        return ends
    corners = np.array([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    t = np.clip((corners - p) @ d / n, 0.0, 1.0)  # or a corner nearest the segment
    return float(min(ends, np.hypot(*(p + t[:, None] * d - corners).T).min()))


def stem_name_line_gaps(drawn, *, page: Page | None = None) -> dict[str, tuple[float, str]]:
    """na:label:<k> (a curve's: na:curve:<n>:label) -> (gap from the name's ink box to the nearest drawn line's stroke,
    page units; that line's gid): every line figcheck draws a text against (loops, nucleic-acid runs, ties, tethers,
    bonds), segment by segment. Routes keep a stem's name layout._INK_MARGIN em clear (layout.stem_keep_box)."""
    lay, fig = drawn.layout, drawn.figure
    if not lay.stems or not fig.axes:
        return {}
    page = page or measure(drawn)
    to_page = fig.axes[0].transData.inverted()
    per_px = float(to_page.transform((1, 0))[0] - to_page.transform((0, 0))[0])  # page units per pixel
    ink = {name: tuple(to_page.transform_bbox(box).extents) for _, name, box in page.texts}
    lines = [(gid, to_page.transform(poly), pad * per_px) for gid, polys, _, pad in page.lines for poly in polys]
    out = {}
    labels = [f"na:label:{k}" for k in range(len(lay.stems))]
    labels += [f"na:curve:{n}:label" for n in range(len(getattr(lay, "na_curves", None) or []))]  # curves (na_bend)
    for label in labels:
        box = ink.get(label)
        if box is None:
            continue
        best = (math.inf, "")
        for gid, poly, pad in lines:
            lo, hi = poly.min(axis=0), poly.max(axis=0)
            if _box_gap(box, (*lo, *hi)) - pad >= best[0]:
                continue
            gap = min(_segment_box_gap(p, q, box) for p, q in zip(poly, poly[1:]))
            best = min(best, (max(gap - pad, 0.0), gid))
        out[label] = best
    return out


def _paths_on_page(page: Page, text, to_page) -> list[MPath]:
    """A text's glyph outlines in page units."""
    return [p.transformed(to_page) for p in page.outlines[id(text)]]


# ------------------------------------------------------------------------------------------- scope


def _polymer_chains(model) -> dict[str, list[tuple[int, str]]]:
    """Chain -> its polymer residues (number, insertion code), for protein and nucleic-acid chains (gemmi)."""
    import gemmi

    out = {}
    for chain in model:
        polymer = chain.get_polymer()
        kind = polymer.check_polymer_type() if len(polymer) else gemmi.PolymerType.Unknown
        if kind in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD) or _nucleic(kind):
            out.setdefault(chain.name, []).extend((r.seqid.num, r.seqid.icode.strip()) for r in polymer)
    return out


def _nucleic(kind) -> bool:
    import gemmi

    return kind in (gemmi.PolymerType.Dna, gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid)


def scope_honest(drawn, path, assembly: str = "auto", focus: str = "auto", *, chains=None) -> list[str]:
    """What the figure leaves out without saying so: "chain <id>" for each polymer chain not drawn at all, and
    "nt <chain>:<first>-<last>" for each stretch of at least 4 nucleotides of a drawn nucleic-acid chain not drawn.
    The chains are counted with gemmi from the model read_model read (the asymmetric unit when `assembly` auto drew
    one copy of several), not from foldmap's loaders. Left out on purpose, and so not counted: what the drawn_scope
    reason (the note printed with the figure) covers (assembly copies, the subunits a focus leaves out), chains not
    among `chains` (--chains), and the chains an explicit --focus leaves undrawn (the chains it draws are still
    checked nucleotide by nucleotide), and what na_extent auto crops from a large nucleic acid once the drawn_scope
    reason says so (plan 10 Task 8): its ends and the gaps between the pieces it kept (_cropped_away); a run missing
    inside a kept piece is still counted. A protein chain with no helix or strand is counted as said when the
    reason says it is not drawn (random-100 C15).

    The file's own chains are then counted again with gemmi, apart from the assembly foldmap built (final review B,
    Important 2: 1RM5's assembly 1 is built from A and B alone): a chain of the file that no drawn chain was built
    from is left out, unless the reason says so (io.ModelRead.left_out)."""
    import gemmi

    from .scope import drawn_scope

    lay, model = drawn.layout, drawn.model
    several = getattr(model, "several", None) if assembly == "auto" else None
    in_file = gemmi.read_structure(str(path))[0]
    read = in_file if several else model[0]
    expected = _polymer_chains(read)
    na = lay.nucleic
    strands = na.strands if na is not None else []
    shown: dict[str, set[tuple[int, str]]] = {}  # nucleic-acid chain -> nucleotides drawn

    def show(s: int, nts) -> None:
        lab = strands[s].labels
        shown.setdefault(strands[s].chain, set()).update((lab[i].seq, lab[i].icode.strip()) for i in nts)

    for d in lay.stems:
        for s, i in d.axial:
            show(s, [i])
    for p in getattr(lay, "na_paths", None) or []:
        show(p.strand, p.nts)
    for s in getattr(lay, "na_traces", None) or {}:
        show(s, range(len(strands[s])))
    drawn_chains = {p.sse.chain for p in lay.placed.values()} | {g.sse.chain for g in lay.ghosts} | set(shown)
    scope = drawn_scope(model, lay, assembly, focus)
    covered: set[str] = set(scope.get("undrawn_chains", []))  # said not drawn: no helix or strand (C15)
    if scope["reason"]:
        if several:
            covered |= set(expected) - {c.name for c in model[0]}
        if focus == "auto" and lay.focus_chains and lay.subunits:
            covered |= _focus_covered(set(expected), lay, _nucleic_chains(read), set(scope.get("na_outside", [])))
    chosen = focus not in ("auto", "none")
    cropped = bool(lay.na_cropped and scope.get("crop"))
    out = []
    for chain in sorted(expected):
        if chain in drawn_chains or chain in covered or chosen or (chains and chain not in chains):
            continue
        if cropped and chain in _nucleic_chains(read):
            continue
        out.append(f"chain {chain}")
    file_chains = _polymer_chains(in_file)
    built_from = {_file_chain(c.name, file_chains) for c in model[0]}
    said_out = set(getattr(model, "left_out", None) or []) if scope["chains_reason"] else set()
    asked = {_file_chain(c, file_chains) for c in chains} if chains else None
    for chain in sorted(set(file_chains) - built_from - set(expected)):  # no copy of it in the model drawn from
        if chain in said_out or chosen or (asked is not None and chain not in asked):
            continue
        if cropped and chain in _nucleic_chains(in_file):
            continue
        out.append(f"chain {chain}")
    said = _cropped_away(na, expected) if cropped else {}  # the note says these were cropped, and how to draw them
    for chain in sorted(shown):
        drawn_nts = shown[chain] | said.get(chain, set())
        run: list[tuple[int, str]] = []
        for r in [*expected.get(chain, []), None]:
            if r is not None and r not in drawn_nts:
                run.append(r)
                continue
            if len(run) >= 4:
                out.append(f"nt {chain}:{run[0][0]}{run[0][1]}-{run[-1][0]}{run[-1][1]}")
            run = []
    return out


def _focus_covered(expected: set[str], lay, nucleic: set[str], named: set[str]) -> set[str]:
    """The chains a subunit focus's note covers: every chain outside the focus that is another subunit's, and a
    nucleic-acid chain only when the note names it (F1 review, Important 3: 8QEM's duplex was covered by "one of 12
    subunits")."""
    return {c for c in expected - set(lay.focus_chains) if c not in nucleic or c in named}


def _file_chain(name: str, file_chains) -> str:
    """The chain of the file a drawn chain was built from: itself, or for an assembly copy (io.read_model names
    copies AddNumber: A -> A1, A2) the longest chain of the file its name is followed by digits after."""
    if name in file_chains:
        return name
    found = [c for c in file_chains if name.startswith(c) and name[len(c) :].isdigit()]
    return max(found, key=len) if found else name


def _cropped_away(na, expected: dict) -> dict[str, set[tuple[int, str]]]:
    """Chain -> the nucleotides na_extent auto's crop took away, as its note says: those before the first piece and
    after the last piece the crop kept of the chain, and those of each gap between two of its pieces (NAGap); a
    nucleotide missing inside a kept piece is not among them."""
    out: dict[str, set[tuple[int, str]]] = {}
    pieces: dict[str, list[tuple[int, int]]] = {}  # chain -> (first, last) index into expected, per piece
    for st in na.strands:
        order = {r: k for k, r in enumerate(expected.get(st.chain, []))}
        at = [order[(lab.seq, lab.icode.strip())] for lab in st.labels if (lab.seq, lab.icode.strip()) in order]
        if at:
            pieces.setdefault(st.chain, []).append((min(at), max(at)))
    for g in na.gaps:
        a, b = na.strands[g.before], na.strands[g.after]
        if a.chain == b.chain and a.chain in pieces:
            order = {r: k for k, r in enumerate(expected[a.chain])}
            lo = order.get((a.labels[-1].seq, a.labels[-1].icode.strip()))
            hi = order.get((b.labels[0].seq, b.labels[0].icode.strip()))
            if lo is not None and hi is not None:
                out.setdefault(a.chain, set()).update(expected[a.chain][lo + 1 : hi])
    for chain, spans in pieces.items():
        rows = expected[chain]
        first, last = min(lo for lo, _ in spans), max(hi for _, hi in spans)
        out.setdefault(chain, set()).update(rows[:first] + rows[last + 1 :])
    return out


def _nucleic_chains(model) -> set[str]:
    return {c.name for c in model if len(c.get_polymer()) and _nucleic(c.get_polymer().check_polymer_type())}


# ------------------------------------------------------------------------------------------- order


def order_checks(drawn) -> list[str]:
    """N->C and sheet order as drawn:
    - each chain's loops join its drawn elements one after the next, in sequence (sses order): a loop that joins two
      elements with a drawn element of its chain between them ("<a>><b> skips <c>"), and two consecutive drawn
      elements with no loop and no such loop over them ("no loop <a>><b>"); chains only lent to the figure
      (Layout.partial) and pairs shown side by side but not adjacent (Layout.unjoined) have no loops;
    - each sheet block (Layout.sheets) holds the strands of one sheet of sheets.build_sheets on the backbone's
      bridges, its track side by side in that sheet's order (either way round; a closed barrel from any strand) and
      each accessory strand (Sheet.accessory) in its parent's column or the next (science S6 judges which), each two
      neighbours pointing the same way when they pair in parallel and opposite ways when antiparallel. With
      copy_ghosts off, a strand lent across an unrolled ring's cut stands past the row's end with no ghost
      (Layout.lent_apart): the block it pairs in is read with it, its order and directions judged among the strands
      it draws, and the strand's own block is no sheet to judge."""
    from .dssp import assign_dssp
    from .sheets import build_sheets

    lay, sses = local_layout(drawn.layout), drawn.sses  # a rosette's sheets read in their own copy's frame
    out = []
    order = {s.id: k for k, s in enumerate(sses)}
    chain_of = {s.id: s.chain for s in sses}
    unjoined = {tuple(p) for p in lay.unjoined}
    shown = [s for s in sses if s.id in lay.placed and s.chain not in lay.partial]
    skipping: set[tuple[str, str]] = set()
    for lp in drawn.loops:
        a, b = lp.a_id, lp.b_id
        if a not in order or b not in order or chain_of[a] != chain_of[b]:
            continue
        between = [s.id for s in shown if s.chain == chain_of[a] and order[a] < order[s.id] < order[b]]
        if between:
            out.append(f"{a}>{b} skips {', '.join(between)}")
            ids = [a, *between, b]
            skipping |= set(zip(ids, ids[1:]))
    joined = {(lp.a_id, lp.b_id) for lp in drawn.loops}
    for x, y in zip(shown, shown[1:]):
        pair = (x.id, y.id)
        if x.chain == y.chain and pair not in joined and pair not in skipping and pair not in unjoined:
            out.append(f"no loop {x.id}>{y.id}")
    if not lay.sheets:
        return out
    bridges = assign_dssp(drawn.bb).bridges
    if getattr(lay, "rosette", 0):  # a rosette's sheets are each protomer's own (cli._within_protomer)
        sector, chain = lay.frame.sector, drawn.bb.labels
        bridges = [b for b in bridges if sector.get(chain[b.i].chain, -1) == sector.get(chain[b.j].chain, -2)]
    sheets = {
        frozenset([*(s.id for s in sh.strands), *(a.strand.id for a in sh.accessory)]): sh
        for sh in build_sheets(sses, bridges)
    }
    at = {  # each block's strands as drawn: its own, and the ghosts standing in it for strands drawn elsewhere
        k: {**{i: lay.placed[i] for i in ids}, **{g.sse.id: g for g in block_ghosts(lay, k) if g.sse.id not in ids}}
        for k, ids in enumerate(lay.sheets)
    }
    apart = dict(getattr(lay, "lent_apart", {}) or {})
    whole = [k for k in at if frozenset(at[k]) in sheets]
    for k, ids in enumerate(lay.sheets):
        sheet = sheets.get(frozenset(at[k]))
        if sheet is None and any(set(ids) <= set(at[j]) for j in whole if j != k):
            continue  # strands lent across an unrolled ring's cut (layout._unwrap): checked where their ghosts stand
        if sheet is None and set(ids) <= set(apart):
            continue  # strands lent across the cut with no ghost (copy_ghosts off): read with the block they left
        if sheet is None and (gone := [i for i, j in apart.items() if j == k]):
            sheet = sheets.get(frozenset([*at[k], *gone]))
        if sheet is None:
            out.append(f"sheet {k}: {', '.join(at[k])} is no sheet of the bridges")
            continue
        want = [s.id for s in sheet.strands if s.id in at[k]]  # (one lent apart: the rest in order)
        got = _drawn_order(at[k])
        seen, wanted = got, want
        if sheet.accessory:  # the track in order; each accessory strand in a column next to its parent's
            drawn = {i: float(np.array([at[k][i].cx, at[k][i].cy]) @ _across(at[k])) for i in at[k]}
            place = sheet_columns(drawn.values())
            seen = [i for i in got if i in want]
            for a in sheet.accessory:
                i, p = a.strand.id, a.parent.id
                if i in drawn and p in drawn and abs(place[drawn[i]] - place[drawn[p]]) > 1:
                    out.append(f"sheet {k}: {i} not drawn at {p}, the strand it pairs with")
        if not _same_order(seen, wanted, sheet.closed):
            out.append(f"sheet {k}: drawn {', '.join(got)}, paired {', '.join(want)}")
        track = [s.id for s in sheet.strands]
        for (a, da), (b, db) in zip(zip(track, sheet.directions), zip(track[1:], sheet.directions[1:])):
            if a not in at[k] or b not in at[k]:
                continue
            pa, pb = at[k][a], at[k][b]
            same = float(pa.direction @ pb.direction) > 0
            if same != (da == db):
                out.append(f"sheet {k}: {a} and {b} drawn {'parallel' if same else 'antiparallel'}")
    return out


def _across(at: dict) -> np.ndarray:
    """The page direction across drawn strands (id -> its drawn strand or ghost): normal to their mean direction."""
    ps = list(at.values())
    ref = ps[0].direction
    mean = sum(p.direction * (1.0 if p.direction @ ref >= 0 else -1.0) for p in ps)
    return np.array([mean[1], -mean[0]])


def _drawn_order(at: dict) -> list[str]:
    """Strand ids (id -> its drawn strand, or the ghost standing for it) in the order they stand side by side on
    the page: across their mean direction."""
    across = _across(at)
    return [k for _, k in sorted(((float(np.array([p.cx, p.cy]) @ across), k) for k, p in at.items()))]


def _same_order(got: list[str], want: list[str], closed: bool) -> bool:
    if closed:
        n = len(want)
        turns = [want[k:] + want[:k] for k in range(n)]
        return any(got in (t, t[::-1]) for t in turns)
    return got in (want, want[::-1])


# ------------------------------------------------------------------------------------------- science (S1-S5)


def drawn_columns(layout, k: int) -> list[list[tuple[str, bool]]]:
    """Sheet block k (Layout.sheets) as drawn: its columns left to right (across the strands' mean direction, as
    _drawn_order reads them), each the (strand id, is a ghost) standing in it; strands at one place across share a
    column (stacked one above the other)."""
    at = [(i, layout.placed[i], False) for i in layout.sheets[k]]
    at += [(g.sse.id, g, True) for g in block_ghosts(layout, k)]
    across = _across({n: p for n, (_, p, _) in enumerate(at)})
    cols: dict[float, list[tuple[str, bool]]] = {}
    for i, p, ghost in at:
        cols.setdefault(round(float(np.array([p.cx, p.cy]) @ across), 3), []).append((i, ghost))
    return [sorted(cols[x]) for x in sorted(cols)]


def _span(p) -> tuple[float, float]:
    """A drawn strand's extent along the page's up: (bottom, top) of its two ports."""
    return tuple(sorted((float(p.n_port[1]), float(p.c_port[1]))))


def residue_places(drawn, k: int):
    """science.Place of residue k (Backbone index) on its drawn figure: its element's axis (whole length; expected at
    its residue's fraction along it), the loop route it is in (from its nearer element end, RESIDUE_STEP per residue;
    expected there, or in proportion to its residue distances when the route is too short to hold that), or a chain
    end's stub (expected at its end); None when its chain or element is not drawn. Worked out apart from
    render.residue_point, which S4 judges."""
    from .layout import END_STUB
    from .science import RESIDUE_STEP

    lay, loops = drawn.layout, drawn.loops
    mine = [s for s in drawn.sses if s.chain == lay.res_chain[k] and s.id in lay.placed]
    for s in mine:
        if s.start <= k <= s.end:
            p = lay.placed[s.id]
            axis = [tuple(p.n_port), tuple(p.c_port)]
            length = float(math.hypot(*np.subtract(axis[1], axis[0])))
            return axis, 0, length, length * (k - s.start) / max(s.end - s.start, 1)
    before, after = [s for s in mine if s.end < k], [s for s in mine if s.start > k]
    if before and after:
        a, b = before[-1], after[0]
        loop = next((lp for lp in loops if lp.a_id == a.id and lp.b_id == b.id), None)
        if loop is None:
            return None
        d_a, d_b = k - a.end, b.start - k
        pts = [tuple(q) for q in loop.points]
        run = sum(math.hypot(*np.subtract(q, p)) for p, q in zip(pts, pts[1:]))
        bound = min(d_a, d_b) * RESIDUE_STEP
        return pts, 0 if d_a <= d_b else -1, bound, min(bound, run * min(d_a, d_b) / (d_a + d_b))
    if before:
        p = lay.placed[before[-1].id]
        return [tuple(p.c_port), tuple(np.add(p.c_port, p.exit_c * END_STUB))], 0, END_STUB, END_STUB
    if after:
        p = lay.placed[after[0].id]
        return [tuple(p.n_port), tuple(np.add(p.n_port, p.exit_n * END_STUB))], 0, END_STUB, END_STUB
    return None


def residue_marks(drawn) -> list:
    """The residue marks of a drawn figure for science.mark_findings: each linked residue's star (gid
    "linked:<chain>:<number>") and each data mark (color_by=data categories, gid "data:<chain>:<number>"), its drawn
    centre against the place of its own residue. A mark naming no residue of the layout is judged with no place."""
    lay = drawn.layout
    linked = {f"{lay.res_labels[k].chain}:{lay.res_labels[k].tag}": k for k in lay.res_linked}
    every = None  # every residue by its key, built when a data mark is met
    out = []
    for patch in drawn.figure.axes[0].patches:
        gid = patch.get_gid() or ""
        if gid.startswith("linked:"):
            k = linked.get(gid[len("linked:") :])
            if k is not None:
                out.append(
                    (gid[len("linked:") :], tuple(np.mean(patch.get_xy()[:10], axis=0)), residue_places(drawn, k))
                )
        elif gid.startswith("data:"):
            if every is None:
                every = {f"{lab.chain}:{lab.tag}": k for k, lab in enumerate(lay.res_labels)}
            k = every.get(gid[len("data:") :])
            out.append((gid, _mark_centre(patch), None if k is None else residue_places(drawn, k)))
    return out


def _mark_centre(patch) -> tuple[float, float]:
    """A mark's centre: a circle's, else its outline's points' mean (the closing point counted once)."""
    if hasattr(patch, "center"):
        return tuple(float(v) for v in patch.center)
    xy = np.asarray(patch.get_xy(), float)
    if len(xy) > 1 and np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    return tuple(float(v) for v in xy.mean(axis=0))


def science_checks(drawn) -> dict[str, list[str] | None]:
    """The science checks (science.SCIENCE) of a drawn figure, apart from the page checks: findings per check, None
    for those not yet measured (S5). Sheet block k is named "sheet <k>" (Layout.sheets, from 0). Every strand and
    ghost stands in its drawn column with its drawn span (drawn_columns, _span); S6 reads the pairings the layout
    could not draw (Layout.unshown)."""
    from .dssp import assign_dssp
    from .science import (
        NOT_YET,
        adjacent_pairs,
        bridge_counts,
        count_finding,
        mark_findings,
        neighbour_findings,
        open_finding,
        partner_findings,
    )

    lay = local_layout(drawn.layout)  # a rosette's copies each in its own frame: sheets read across their strands
    out: dict[str, list[str] | None] = {"S1": [], "S2": [], "S3": [], **dict.fromkeys(NOT_YET), "S6": []}
    out["S4"] = mark_findings(residue_marks(drawn))
    if not lay.sheets:
        return out
    count = bridge_counts(drawn.sses, assign_dssp(drawn.bb).bridges)
    unshown = {tuple(p) for p in getattr(lay, "unshown", [])}
    for k, ids in enumerate(lay.sheets):
        columns = drawn_columns(lay, k)
        ghosts = {g.sse.id: g for g in block_ghosts(lay, k)}
        spots = [(i, c, _span(ghosts[i] if g else lay.placed[i]), g) for c, col in enumerate(columns) for i, g in col]
        out["S6"] += [f"sheet {k}: {t}" for t in partner_findings(spots, lay.accessory, unshown)]
        own_ghosts = [i for i in ghosts if i in ids]
        s1, s2 = neighbour_findings(spots, count, own_ghosts)
        out["S1"] += [f"sheet {k}: {t}" for t in s1]
        out["S2"] += [f"sheet {k}: {t}" for t in s2]
        lent = [i for i in ghosts if i not in ids]
        if own_ghosts:  # a barrel: its ring is its strands (accessories aside) and the ghosts lent to it
            ring = [i for i in ids if i not in lay.accessory] + lent
            others = [i for i in ids if i in lay.accessory]
            found = count_finding(f"sheet {k}", ring, others, count, adjacent_pairs(spots))
        else:  # drawn open: no ring of _BARREL_MIN strands or more among its own strands
            found = open_finding(f"sheet {k}", ids, count)
        if found:
            out["S3"].append(found)
    return out


def science_judged(drawn) -> dict[str, int]:
    """What each science check had to judge on the drawn figure: S1 and S3 its sheet blocks, S2 its closed barrels
    (a block drawing a ghost of its own strand), S4 its residue marks, S6 its accessory strands. A check with nothing
    to judge (0) is not judged, never passed (final review B, Important 3); S5 is the campaign's."""
    lay = drawn.layout
    closed = sum(1 for k, ids in enumerate(lay.sheets) if any(g.sse.id in set(ids) for g in block_ghosts(lay, k)))
    return {
        "S1": len(lay.sheets),
        "S2": closed,
        "S3": len(lay.sheets),
        "S4": len(residue_marks(drawn)),
        "S6": len(lay.accessory),
    }


def accessories_named(layout) -> list[str]:
    """The accessory strands (Layout.accessory) the figure names rather than draws level beside their partner
    (Layout.unshown, read from the drawing), worded as the notes word them: "A1:2-6 also pairs with A1:9-23". S6
    passes such a strand because its pairing is named, so this count is the measure of accessory placement; it is
    reported (figure_checks), not a gate (plan 12 ruling)."""
    return [
        f"{a} also pairs with {b}"
        for a, b in sorted(map(tuple, getattr(layout, "unshown", [])))
        if layout.accessory.get(a) == b
    ]


# ------------------------------------------------------------------------------------------- membrane (A18)
# Plan 27 Task 2 (rulings P3, P4, P6): A18 judges each element's drawn body (Placed.rect, its full drawn extent, a
# helix's and an arrow's ends included) and every segment of every loop against the band drawn, with no tolerance but
# float snapping (BAND_EPS). What it let pass before is gone: the centre-only test (ruling R4: a soluble element was in
# the slab only when its drawn centre was), the 5 Å an element beyond a face could reach in (T5-5) and the loops'
# runs excused beside their ports (C26, the membrane slab ruling). The band's elements are membrane.band_set's, the
# one rule the layout lays the band out by.
# (membrane.BAND_EPS; the page reading of the band, page_depth and the rest, is shared with the router there)


def membrane_view(layout, sses, bb, mem, loops=None) -> dict:
    """A18 on a membrane view (layout.membrane, the membrane `mem` it was drawn for):
    - tm_helices, the helices that cross the membrane (membrane.crosses: the CA run crosses its mid-plane and spans
      half its thickness); tm_short, those that do not cross their band (crosses_band), with the share of the band
      they cover;
    - slab_helix, the helices of the band that cross no membrane, their CA median in the slab (membrane.slab_helix,
      ruling P10), whose axis (its two ports: ruling P13, the one part not judged on the body) lies off the range
      membrane.slab_range gives it (Å, the furthest): a standing one inside the slab within its CA atoms' depths, a
      flat one within them extended to its nearer face, where it may lie (Fig. 1C's NS1, NS7);
    - band_member, the elements not of the band (membrane.in_band: crossing helices, helices in the slab at any tilt,
      crossing strands) whose drawn body overlaps a band;
    - strand_in_band, those of them that are β-strands (an interfacial strand is never in the membrane);
    - outside_in_band, the elements not of the band drawn nearer the mid-plane, on their own side, than their own
      nearest CA atom (membrane.own_bound; Å nearer);
    - wrong_side, the elements not of the band whose drawn centre lies on the side opposite their own
      (membrane.own_side; Å past the mid-plane);
    - loop_in_band, the loops with a segment drawn nearer the mid-plane than their bound (loop_side_bound: their
      nearest residue, or their face when every residue lies beyond it; ruling P14) but on the stub straight out of a
      band element's port to that face, or in the lane at its own port, beside its element, from a port drawn past
      the bound to the bound, when its residues lie in the bilayer (their median in the slab; P36-P38); never through
      its own element (loops_in_band);
    - loop_crossing, the loops that join the two faces (their bound past the far face, or their end residues on
      opposite sides and their ports either side of the band) drawn across the band anywhere but in the strip at one
      of their own ports or past the band's ends (P34, P39; loops_crossing; None: `loops` not given).
    Each entry is "<id> <Å>" (findings_text: "<0.1" under 0.05 Å); band_member and strand_in_band say how far in from
    the element's own face. A barrel's ghost strand is judged as its strand, "<id> ghost <Å>".
    Loops are judged on their routed polyline. The page rounds its corners (final review M2): a rounded corner
    (render._CORNER, 0.3) and a curved one (Style loops curved, kept out of the band beside its ports' strips by
    render._beside_strips) stay within their two legs, never nearer the mid-plane than the polyline, and a corner at a
    port's stub stays inside the strip at that port, though it may cut a little inside the face there."""
    bands = layout.membrane.get("bands") or [layout.membrane["y"]]
    up = 1 if layout.membrane.get("inside", "below") == "below" else -1  # page-up is the outside's direction
    half = mem.half
    depth = mem.depth(bb.ca)
    band = band_set(sses, depth, mem)
    ghosts: dict[str, list] = {}
    for g in getattr(layout, "ghosts", []):
        ghosts.setdefault(g.sse.id, []).append(g)
    short, tm = [], 0
    member, strands, nearer, wrong, flats = [], [], [], [], []
    for s in sses:
        p = layout.placed.get(s.id)
        if p is None:
            continue
        seg = depth[s.start : s.end + 1]
        if s.kind == "H" and crosses(seg, half):
            tm += 1
            lo, hi = sorted((p.n_port[1], p.c_port[1]))
            y0, y1 = min(bands, key=lambda b: abs((b[0] + b[1]) / 2 - (lo + hi) / 2))
            if not crosses_band((lo, hi), (y0, y1)):
                short.append(f"{s.id} {(min(hi, y1) - max(lo, y0)) / max(y1 - y0, 1e-9):.2f}")
        if s.id in band:
            if slab_helix(s, depth, mem):  # its axis at its own depth range (P10, P13)
                ends = (p.n_port[1], p.c_port[1])
                b = _nearest_band(bands, sum(ends) / 2)
                lo, hi = sorted(_drawn_depth(y, b, half, up) for y in ends)
                low, high = slab_range(s, depth, mem)
                off = max(low - lo, hi - high)
                if off > _eps(b, half):
                    flats.append(f"{s.id} {_angstrom(off)}")
            continue
        side, closest = own_bound(seg)
        for q, name in [(p, s.id), *((g, f"{s.id} ghost") for g in ghosts.get(s.id, []))]:
            _, y0, _, y1 = q.rect
            into = _band_reach(y0, y1, bands, half, up, side)
            if into is not None:
                member.append(f"{name} {_angstrom(into)}")
                if s.kind == "E":
                    strands.append(f"{name} {_angstrom(into)}")
            b = _nearest_band(bands, (y0 + y1) / 2)
            eps = _eps(b, half)
            edge = min(side * _drawn_depth(y, b, half, up) for y in (y0, y1))  # its edge nearer the mid-plane
            if edge < closest - eps:
                nearer.append(f"{name} {_angstrom(closest - edge)}")
            centre = side * _drawn_depth((y0 + y1) / 2, b, half, up)
            if centre < -eps:
                wrong.append(f"{name} {_angstrom(-centre)}")
    judged = (loops, sses, depth, band, bands, mem, up, layout.placed)
    return {
        "tm_helices": tm,
        "tm_short": short,
        "slab_helix": flats,
        "band_member": member,
        "strand_in_band": strands,
        "outside_in_band": nearer,
        "wrong_side": wrong,
        "loop_in_band": None if loops is None else loops_in_band(*judged),
        "loop_crossing": None if loops is None else loops_crossing(*judged),
    }


def _angstrom(x: float) -> str:
    """A finding's Å, one decimal; one under 0.05 Å reads "<0.1", never "0.0" (it is a finding, not a pass)."""
    return "<0.1" if x < 0.05 else f"{x:.1f}"


def _loop_ends(placed, lp):
    """The drawn elements a loop joins (Layout.placed), or a ValueError naming the loop when one is not drawn (A18
    judges a loop by its own ports and boxes; route_loops routes only loops whose two ends are placed)."""
    missing = [k for k in (lp.a_id, lp.b_id) if k not in placed]
    if missing:
        raise ValueError(f"A18: loop {lp.a_id}>{lp.b_id} joins an element not drawn: {', '.join(missing)}")
    return placed[lp.a_id], placed[lp.b_id]


def loops_in_band(loops, sses, depth, band: set[str], bands, mem, up: int, placed) -> list[str]:
    """A18's loop_in_band (membrane_view; rulings P14, P36, P38): each loop (a>b) with a segment drawn nearer the
    mid-plane, on its own side, than its bound (membrane.loop_side_bound, from its residues, and its two end residues
    when they lie on opposite sides: its nearest residue, or its face when all lie beyond it), how much nearer, in Å
    (membrane.loop_excess, the measure the router's wall keeps to). Out of the port of an element of the band (`band`)
    the loop may run straight to its face along the page's normal, and cross the face line there, never through its own
    element's box. From an end the layout draws past its bound it may climb the lane at that port (membrane.loop_lanes,
    lane, in_lane: the port's x ± half its element's width and LANE, from the port to the bound; ruling P36a: a band
    element's end, or either end of a two-face loop, P36b), when its residues lie in the bilayer; nowhere else.
    `placed`: the drawn elements (Layout.placed), their ports and boxes. `up`: +1 when page-up is the outside."""
    by_id = {s.id: s for s in sses}
    out = []
    for lp in loops:
        a, b = by_id.get(lp.a_id), by_id.get(lp.b_id)
        if a is None or b is None:
            continue
        pa, pb = _loop_ends(placed, lp)
        side, bound = loop_side_bound(a, b, depth, mem)
        lanes = loop_lanes(a, b, depth, mem, ((pa.c_port, pa.rect), (pb.n_port, pb.rect)), band, bands, up)
        stubs = (lp.a_id in band, lp.b_id in band)
        worst = loop_excess(lp.points, side, bound, bands, mem.half, up, *stubs, lanes=lanes, boxes=(pa.rect, pb.rect))
        if worst is not None:
            out.append(f"{lp.a_id}>{lp.b_id} {_angstrom(worst)}")
    return out


def loops_crossing(loops, sses, depth, band: set[str], bands, mem, up: int, placed) -> list[str]:
    """A18's loop_crossing (membrane_view; rulings P34, P39): each loop that joins the two faces
    (membrane.loop_crossing: its bound past the far face, or its two end residues on opposite sides of the mid-plane
    and its ports drawn either side of the band) and runs inside the band anywhere but in the strip at one of its own
    ports across the band or past the band's ends (membrane.band_ends: the router's clearance beyond the outermost
    band element), with how far inside the band it runs there, in Å from the nearer face (membrane.crossing_excess,
    the router's own check)."""
    from .route import _INFLATE

    by_id = {s.id: s for s in sses}
    ends = band_ends([placed[k].rect for k in sorted(band) if k in placed], _INFLATE)
    out = []
    for lp in loops:
        a, b = by_id.get(lp.a_id), by_id.get(lp.b_id)
        if a is None or b is None:
            continue
        pa, pb = _loop_ends(placed, lp)
        cross, lanes = loop_crossing(a, b, depth, mem, ((pa.c_port, pa.rect), (pb.n_port, pb.rect)), band, bands, up)
        worst = crossing_excess(lp.points, lanes, ends, bands, mem.half, up) if cross else None
        if worst is not None:
            out.append(f"{lp.a_id}>{lp.b_id} {_angstrom(worst)}")
    return out


def _band_reach(y0: float, y1: float, bands, half: float, up: int, side: int) -> float | None:
    """How far (Å) a body or segment drawn from page y0 to y1 (y0 <= y1) reaches into the bands it overlaps, from the
    face on its own `side` (+1 outside), past the other face counted no further; None when it overlaps none (touching
    a face, within BAND_EPS, is not overlapping)."""
    out = None
    for b0, b1 in bands:
        if not (y1 > b0 + BAND_EPS and y0 < b1 - BAND_EPS):
            continue
        deepest = min(side * _drawn_depth(y, (b0, b1), half, up) for y in (max(y0, b0), min(y1, b1)))
        into = half - max(deepest, -half)
        out = into if out is None else max(out, into)
    return out


# ------------------------------------------------------------------------------------------- all of them


def figure_checks(drawn, path, *, assembly: str = "auto", focus: str = "auto", chains=None) -> dict:
    """Every check on one drawn figure, measured once: a count (or the largest value) under each name, the findings
    under "details", the font files the text resolved to under "font". When that font is not INK_FONT the ink
    checks (INK_CHECKS) are not run: their counts are None and "skipped" says why."""
    page = measure(drawn)
    skip = page.font_reason
    collisions = None if skip else text_collisions(drawn, page=page)
    differ: dict[tuple[str, str], str] = {}
    pairs = _copy_pairs(drawn, differ)
    details = {
        "element_overlaps": element_overlaps(drawn.layout),
        "text_collisions": collisions,
        "off_page": None if skip else off_page(drawn, page=page),
        "label_repeats": label_repeats(drawn.layout),
        "copies_alike": {
            pair: {"deviation": dev, "matched": n, "of": m, "judged": why or "yes"}
            | ({"sheets_differ": differ[pair]} if pair in differ else {})  # left out of A15 (C2b-3 R2)
            for pair, (dev, n, m, why) in pairs.items()
        },
        "copies_untied": copies_untied(drawn),  # copies with no tie when none is drawn: A15 not judged, not n/a
        "domains": domain_checks(drawn, collisions or []) if drawn.layout.domains else {},
        "stem_name_gaps": None if skip else stem_name_gaps(drawn, page=page),
        "scope_silent": scope_honest(drawn, path, assembly, focus, chains=chains),
        "order": order_checks(drawn),
        "accessories_named": accessories_named(drawn.layout),  # reported, not gated (plan 12 ruling)
    }
    domains = details["domains"]
    gaps = details["stem_name_gaps"]
    return {
        "element_overlaps": len(details["element_overlaps"]),
        "text_collisions": None if skip else len(collisions),
        "off_page": None if skip else len(details["off_page"]),
        "label_repeats": sum(len(v) for v in details["label_repeats"].values()),
        "copies_alike": max((dev for dev, *_ in pairs.values() if dev is not None), default=0.0),
        "domain_rules": sum(len(domains.get(k, [])) for k in ("cross_chain", "unalike", "small")),
        "domain_overlap": max((f for *_, f in domains.get("panel_overlaps", [])), default=0.0),
        "domain_name_collisions": None if skip else len(domains.get("name_collisions", [])),
        "stem_name_gap": None if skip else max((own for own, _ in gaps.values()), default=0.0),
        "scope_silent": len(details["scope_silent"]),
        "order": len(details["order"]),
        "accessories": len(drawn.layout.accessory),  # accessory strands drawn off a sheet's track
        "accessories_named": len(details["accessories_named"]),  # of them, named rather than drawn level: reported
        "loop_crossings": loop_crossings(drawn.loops),
        "honest_crossings": honest_crossings(drawn.layout, drawn.loops),
        "font": page.fonts,
        "skipped": dict.fromkeys(INK_CHECKS, skip) if skip else {},
        "details": details,
    }
