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

from .layout import block_ghosts, joined_pieces, lying_over

# Lines a text can be crossed by (render's gids, matched by prefix). A loop's casing (loopcase:) is drawn in the
# page colour and a chain-end stub belongs to its terminus letter: neither is read.
LINE_GIDS = ("loop:", "loop-seg:", "loop-arrow:", "na:run:", "na:tie:", "ligand-tether:", "disulfide:", "hbond:")
GLYCAN_STEM = re.compile(r"glycan:[^:]+:stem$")
# Bodies: opaque shapes a text must not sit on unless it is theirs, and that hide a line drawn under them. Element
# bodies are owned by their element id (label:<id>, resnum:<id>:<end>, ghostlabel:<id> are theirs).
ELEMENT_BODIES = ("strand:", "helix:", "helix-back:", "eta:", "ghost:")
OWN_BODIES = {"ligand:": "ligand-label:", "linked:": "linked-label:", "domain-name-plate:": "domain-label:"}
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


def element_overlaps(lay) -> list[tuple[str, str]]:
    """Pairs of element or stem boxes (their rects) that overlap, except two consecutive pieces of one bent duplex,
    which meet at their joint (layout.joined_pieces), and a helix drawn over a sheet with that sheet's strands
    (layout.lying_over). nacheck.render_checks counts these as its overlaps."""
    joined = joined_pieces(lay.stems, lay.nucleic) if lay.nucleic is not None else set()
    joined |= lying_over(lay)
    boxes = [(p.sse.id, p.rect) for p in lay.placed.values()] + [(d.id, d.rect) for d in lay.stems]
    return [
        (i, j)
        for (i, a), (j, b) in combinations(boxes, 2)
        if min(a[2], b[2]) - max(a[0], b[0]) > 1e-6
        and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6
        and frozenset((i, j)) not in joined
    ]


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
            elif gid.startswith(ELEMENT_BODIES) or gid.startswith(tuple(OWN_BODIES)) or gid.startswith("na:stem:"):
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
        its first pair's numbers and its nucleotides' letters."""
        lay = self.layout
        if lay is None or not lay.stems or self.ax is None:
            return []
        from .layout import dna_to_page

        na = lay.nucleic
        out = []
        for k, d in enumerate(lay.stems):
            ends = ((-1, -1), (1, -1), (1, 1), (-1, 1))
            corners = [dna_to_page(d, x * d.length / 2, y * d.width / 2) for x, y in ends]
            path = MPath(self.ax.transData.transform(np.array([*corners, corners[0]])), closed=True)
            own = {f"na:label:{k}", f"na:stem:{k}:first"}
            if na is not None:
                own |= {f"na:nt:{na.strands[s].labels[i].chain}:{na.strands[s].labels[i].tag}" for s, i in d.axial}
            out.append((f"na:stem:{k}", path, own))
        return out

    def own_lines(self) -> dict[str, set[str]]:
        """Text gid -> the gids of the lines it belongs to: a loop's length, a gap's label, a tie's letters, a
        ligand's tethers, a strand's letters and ends to its own drawn runs."""
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
    - text over a body (an element, a stem, a ligand or linked mark, a domain name plate) it does not belong to (its
      label, its residue numbers, its letters, its name), by more than TEXT_OVERLAP_PX2 of its ink deeper than
      TEXT_DEPTH_PX inside it, once per element (a helix
      is drawn as helix: and helix-back:); a helix drawn over a sheet and that sheet's strands are each other's
      (layout.lying_over);
    - text crossed by a line (LINE_GIDS: loops, nucleic-acid runs and ties, tethers, disulfide bars, glycan stems,
      H-bonds) whose stroke comes within `clearance_px` of its ink box where no opaque body drawn above the line
      hides it there: a line under the text still runs through its letters, a loop under the strand its letter sits
      on does not show. A text's own line (a loop's length, a gap's label, a tie's letter, a ligand's tethers, a
      strand's letters and ends) is not counted."""
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
    bodies = [(body, path, _owner(body) or body) for body, path, _, _ in page.bodies if not body.startswith("na:stem:")]
    body_rows = _extents([path.get_extents() for _, path, _ in bodies])  # each body measured once
    stem_rows = _extents([path.get_extents() for _, path, _ in page.stems])
    line_rows = _extents([extent for _, _, extent, _, _ in lines])
    pads = np.array([pad for _, _, _, pad, _ in lines], float)
    for t, name, box in page.texts:
        gid = t.get_gid() or ""
        mine = _owner(gid)
        met: set[str] = set()
        core = _core(box)
        if core is not None:  # stems: their whole body below
            for k in _meeting(core.extents, body_rows, True).tolist():
                body, path, theirs = bodies[k]
                if theirs == mine or theirs in met or frozenset((theirs, mine)) in over:
                    continue
                if _on(box, path, Bbox(body_rows[k].reshape(2, 2))):
                    met.add(theirs)
                    hits.append(f"{name} × {body.replace('helix-back:', 'helix:', 1)}")
            for k in _meeting(core.extents, stem_rows, True).tolist():
                body, path, own = page.stems[k]
                if gid not in own and _on(box, path, Bbox(stem_rows[k].reshape(2, 2))):
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
    """Chain -> the labels two or more of its drawn elements carry; "<na>" -> stem names given to two or more
    stems (ghosts aside)."""
    count: dict[str, Counter] = {}
    for p in lay.placed.values():
        count.setdefault(p.sse.chain, Counter())[p.label] += 1
    repeated = {c: sorted(k for k, n in got.items() if n > 1) for c, got in count.items()}
    names = Counter(d.name for d in lay.stems if not d.ghost and d.name)
    repeated["<na>"] = sorted(k for k, n in names.items() if n > 1)
    return {c: ks for c, ks in repeated.items() if ks}


def fit_deviation(a, b, mirrored: bool = False) -> float:
    """The largest distance between matched points a[k] and b[k] after the best rigid motion of b onto a (least
    squares: a rotation and a translation; a reflection and a translation when `mirrored`)."""
    a, b = np.asarray(a, float).reshape(-1, 2), np.asarray(b, float).reshape(-1, 2)
    if mirrored:
        b = b * [-1.0, 1.0]
    ca, cb = a.mean(axis=0), b.mean(axis=0)
    u, _, vt = np.linalg.svd((b - cb).T @ (a - ca))
    turn = u @ vt
    if np.linalg.det(turn) < 0:  # the best proper rotation
        u[:, -1] *= -1
        turn = u @ vt
    fitted = (b - cb) @ turn + ca
    return float(np.hypot(*(fitted - a).T).max())


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


def _aligned(residues: dict[str, list[str]], a: str, b: str) -> list[tuple[int, int]]:
    return aligned_positions(residues[a], residues[b])


def _copy_pairs(drawn) -> dict[tuple[str, str], tuple[float | None, int, int, str]]:
    """(copy, reference) -> (deviation, elements matched, elements of the smaller of the two, why not judged) for each
    tie of Layout.copy_ties, and for each copy no tie holds (_untied_copies). The copies' elements pair one-to-one
    by sequence (C18 ruling): the copy's residues aligned to the reference's (aligned_positions), its element ranges
    taken to the reference's residues, then model.pair_places (the same kind, sharing at least half of the shorter).
    Strands of a sheet shared between the two copies are left out (shared_sheet_elements). The deviation is that of
    fit_deviation (a reflection for a mirror tie; for an untied copy only when the figure mirrors a tie). When fewer
    than COPY_MIN_MATCHED or than half the smaller copy's elements match (two points always fit), the tie is not
    judged (deviation None) if the smaller copy has fewer than COPY_MIN_MATCHED elements to match ("n/a: k
    elements") or a shared sheet took the rest ("n/a: shared sheet"), C22 and C10 rulings; else it is inf."""
    from .model import pair_places

    lay = drawn.layout
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
    out = {}
    for copy, ref, mirrors in ties:
        shared = shared_sheet_elements(lay, copy, ref)
        a = [p for k, p in mine.get(ref, []) if k not in shared]
        b = [p for k, p in mine.get(copy, []) if k not in shared]
        to_ref = {j: i for i, j in _aligned(residues, ref, copy)} if ref in residues and copy in residues else {}

        def place(p, chain, mapping=None):
            local = range(p.sse.start - first[chain], p.sse.end - first[chain] + 1)
            got = [mapping[j] for j in local if j in mapping] if mapping is not None else list(local)
            return (p.sse.kind, min(got), max(got)) if got else (p.sse.kind, -2, -3)

        found = sorted(pair_places([place(p, ref) for p in a], [place(p, copy, to_ref) for p in b]).items())
        smaller = min(len(a), len(b))
        why = ""
        if len(found) >= COPY_MIN_MATCHED and 2 * len(found) >= smaller:
            pa = [(a[i].cx, a[i].cy) for i, _ in found]
            pb = [(b[j].cx, b[j].cy) for _, j in found]
            dev: float | None = min(fit_deviation(pa, pb, m) for m in mirrors)
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
    the gap from the name to the nearest letter."""
    lay, fig = drawn.layout, drawn.figure
    if not lay.stems or not fig.axes:
        return {}
    page = page or measure(drawn)
    to_page = fig.axes[0].transData.inverted()
    texts = {name: t for t, name, _ in page.texts}
    ink = {name: tuple(to_page.transform_bbox(box).extents) for _, name, box in page.texts}
    na = lay.nucleic
    out = {}
    for k, d in enumerate(lay.stems):
        name = ink.get(f"na:label:{k}")
        if name is None:
            continue
        tags = [f"na:nt:{na.strands[s].labels[i].chain}:{na.strands[s].labels[i].tag}" for s, i in d.axial]
        frame = Affine2D().translate(-d.cx, -d.cy).rotate(-d.angle)  # page -> (along, across) the stem
        local = {
            g: Bbox.union([p.transformed(frame).get_extents() for p in _paths_on_page(page, texts[g], to_page)])
            for g in [f"na:label:{k}", *tags]
            if g in texts
        }
        letters = [local[g] for g in tags if g in local]
        if not letters:
            own = _box_gap(name, d.rect)
        else:
            box = local[f"na:label:{k}"]
            side = 1.0 if box.y0 + box.y1 > 0 else -1.0
            row = [b for b in letters if (b.y0 + b.y1) * side > 0]
            across = (
                (box.y0 - max(b.y1 for b in row) if side > 0 else min(b.y0 for b in row) - box.y1) if row else -math.inf
            )
            beside = row and box.x1 > min(b.x0 for b in row) and box.x0 < max(b.x1 for b in row)
            own = across if beside and across > 0 else min(_box_gap(box.extents, b.extents) for b in letters)
        others = [_box_gap(name, e.rect) for j, e in enumerate(lay.stems) if j != k]
        out[f"na:label:{k}"] = (float(own), min(others) if others else math.inf)
    return out


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
    """na:label:<k> -> (gap from the name's ink box to the nearest drawn line's stroke, page units; that line's gid):
    every line figcheck draws a text against (loops, nucleic-acid runs, ties, tethers, bonds), segment by segment.
    Routes keep a stem's name layout._INK_MARGIN em clear (layout.stem_keep_box)."""
    lay, fig = drawn.layout, drawn.figure
    if not lay.stems or not fig.axes:
        return {}
    page = page or measure(drawn)
    to_page = fig.axes[0].transData.inverted()
    per_px = float(to_page.transform((1, 0))[0] - to_page.transform((0, 0))[0])  # page units per pixel
    ink = {name: tuple(to_page.transform_bbox(box).extents) for _, name, box in page.texts}
    lines = [(gid, to_page.transform(poly), pad * per_px) for gid, polys, _, pad in page.lines for poly in polys]
    out = {}
    for k in range(len(lay.stems)):
        box = ink.get(f"na:label:{k}")
        if box is None:
            continue
        best = (math.inf, "")
        for gid, poly, pad in lines:
            lo, hi = poly.min(axis=0), poly.max(axis=0)
            if _box_gap(box, (*lo, *hi)) - pad >= best[0]:
                continue
            gap = min(_segment_box_gap(p, q, box) for p, q in zip(poly, poly[1:]))
            best = min(best, (max(gap - pad, 0.0), gid))
        out[f"na:label:{k}"] = best
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
    reason says it is not drawn (random-100 C15)."""
    import gemmi

    from .scope import drawn_scope

    lay, model = drawn.layout, drawn.model
    several = getattr(model, "several", None) if assembly == "auto" else None
    read = gemmi.read_structure(str(path))[0] if several else model[0]
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
            covered |= set(expected) - set(lay.focus_chains)
    chosen = focus not in ("auto", "none")
    cropped = bool(lay.na_cropped and scope.get("crop"))
    out = []
    for chain in sorted(expected):
        if chain in drawn_chains or chain in covered or chosen or (chains and chain not in chains):
            continue
        if cropped and chain in _nucleic_chains(read):
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
      bridges, side by side in that sheet's order (either way round; a closed barrel from any strand; strands
      sharing a column, Sheet.columns, as one), each two
      neighbours pointing the same way when they pair in parallel and opposite ways when antiparallel."""
    from .dssp import assign_dssp
    from .sheets import build_sheets

    lay, sses = drawn.layout, drawn.sses
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
    sheets = {frozenset(s.id for s in sh.strands): sh for sh in build_sheets(sses, bridges)}
    at = {  # each block's strands as drawn: its own, and the ghosts standing in it for strands drawn elsewhere
        k: {**{i: lay.placed[i] for i in ids}, **{g.sse.id: g for g in block_ghosts(lay, k) if g.sse.id not in ids}}
        for k, ids in enumerate(lay.sheets)
    }
    whole = [k for k in at if frozenset(at[k]) in sheets]
    for k, ids in enumerate(lay.sheets):
        sheet = sheets.get(frozenset(at[k]))
        if sheet is None and any(set(ids) <= set(at[j]) for j in whole if j != k):
            continue  # strands lent across an unrolled ring's cut (layout._unwrap): checked where their ghosts stand
        if sheet is None:
            out.append(f"sheet {k}: {', '.join(at[k])} is no sheet of the bridges")
            continue
        want = [s.id for s in sheet.strands]
        got = _drawn_order(at[k])
        seen, wanted = got, want
        if sheet.columns:  # strands sharing a column (sheets._columns) stand one above the other: read by column
            column = {s.id: c for s, c in zip(sheet.strands, sheet.columns)}
            seen, wanted = _by_column(got, column), _by_column(want, column)
        if not _same_order(seen, wanted, sheet.closed):
            out.append(f"sheet {k}: drawn {', '.join(got)}, paired {', '.join(want)}")
        for (a, da), (b, db) in zip(zip(want, sheet.directions), zip(want[1:], sheet.directions[1:])):
            pa, pb = at[k][a], at[k][b]
            same = float(pa.direction @ pb.direction) > 0
            if same != (da == db):
                out.append(f"sheet {k}: {a} and {b} drawn {'parallel' if same else 'antiparallel'}")
    return out


def _drawn_order(at: dict) -> list[str]:
    """Strand ids (id -> its drawn strand, or the ghost standing for it) in the order they stand side by side on
    the page: across their mean direction."""
    ids, ps = list(at), list(at.values())
    ref = ps[0].direction
    mean = sum(p.direction * (1.0 if p.direction @ ref >= 0 else -1.0) for p in ps)
    across = np.array([mean[1], -mean[0]])
    return [k for _, k in sorted(((float(np.array([p.cx, p.cy]) @ across), k) for k, p in zip(ids, ps)))]


def _by_column(ids: list[str], column: dict[str, int]) -> list[str]:
    """Strand ids in order, the strands of one column (side by side in that order) as one entry, their ids
    joined in sorted order: the order of the columns, whichever way the strands sharing one stand."""
    out: list[str] = []
    for k in ids:
        if out and column[k] == column[out[-1].split("|")[0]]:
            out[-1] = "|".join(sorted([*out[-1].split("|"), k]))
        else:
            out.append(k)
    return out


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
    ref = at[0][1].direction
    mean = sum(p.direction * (1.0 if p.direction @ ref >= 0 else -1.0) for _, p, _ in at)
    across = np.array([mean[1], -mean[0]])
    cols: dict[float, list[tuple[str, bool]]] = {}
    for i, p, ghost in at:
        cols.setdefault(round(float(np.array([p.cx, p.cy]) @ across), 3), []).append((i, ghost))
    return [sorted(cols[x]) for x in sorted(cols)]


def science_checks(drawn) -> dict[str, list[str] | None]:
    """The science checks (science.SCIENCE) of a drawn figure, apart from the page checks: findings per check, None
    for those not yet measured (S4, S5). Sheet block k is named "sheet <k>" (Layout.sheets, from 0)."""
    from .dssp import assign_dssp
    from .science import NOT_YET, adjacent_pairs, bridge_counts, count_finding, neighbour_findings

    lay = drawn.layout
    out: dict[str, list[str] | None] = {"S1": [], "S2": [], "S3": [], **dict.fromkeys(NOT_YET)}
    if not lay.sheets:
        return out
    count = bridge_counts(drawn.sses, assign_dssp(drawn.bb).bridges)
    for k, ids in enumerate(lay.sheets):
        columns = drawn_columns(lay, k)
        s1, s2 = neighbour_findings(columns, count)
        out["S1"] += [f"sheet {k}: {t}" for t in s1]
        out["S2"] += [f"sheet {k}: {t}" for t in s2]
        ghosts = [(c, i) for c, col in enumerate(columns) for i, g in col if g]
        if not ghosts:
            continue
        ring = [i for i in ids if i not in lay.accessory] + [i for _, i in ghosts if i not in ids]  # lent: drawn
        # in another block, standing here for a strand of this ring
        others = [i for i in ids if i in lay.accessory]
        found = count_finding(f"sheet {k}", ring, others, count, adjacent_pairs(columns))
        if found:
            out["S3"].append(found)
    return out


# ------------------------------------------------------------------------------------------- membrane (A18)
SLAB_MARGIN = 0.3  # page units: a loop's sideways run this far inside a band's faces is in the slab
_RUN_MIN = 0.5  # page units: a shorter sideways step is no run along the slab


def crosses_band(ends: tuple[float, float], band: tuple[float, float]) -> bool:
    """Does a transmembrane element drawn from page y ends[0] to ends[1] cross the band (bottom, top) as A18 judges
    it: entering it through a face and reaching past its middle (C27 ruling)? Shared by the membrane view's own test
    (cli._off_band) and the campaign's (random100.membrane_checks)."""
    lo, hi = sorted(ends)
    y0, y1 = band
    return lo < (y0 + y1) / 2 < hi and (lo <= y0 or hi >= y1)


def membrane_view(layout, sses, bb, mem, loops=None) -> dict:
    """A18 on a membrane view (layout.membrane, the membrane `mem` it was drawn for): tm_helices, the helices whose
    CA run reaches both faces (within 2 Å); tm_short, those that do not cross their band (crosses_band), with the
    share of the band they cover; loops_in_slab, slab_runs of `loops` (None: not measured)."""
    bands = layout.membrane.get("bands") or [layout.membrane["y"]]
    depth = mem.depth(bb.ca)
    short, tm = [], 0
    for s in sses:
        p = layout.placed.get(s.id)
        if p is None or s.kind != "H":
            continue
        seg = depth[s.start : s.end + 1]
        if len(seg) < 2 or seg.min() > -mem.half + 2 or seg.max() < mem.half - 2:
            continue
        tm += 1
        lo, hi = sorted((p.n_port[1], p.c_port[1]))
        y0, y1 = min(bands, key=lambda b: abs((b[0] + b[1]) / 2 - (lo + hi) / 2))
        if not crosses_band((lo, hi), (y0, y1)):
            short.append(f"{s.id} {(min(hi, y1) - max(lo, y0)) / max(y1 - y0, 1e-9):.2f}")
    slab = None if loops is None else slab_runs(loops, bands, layout.placed)
    return {"tm_helices": tm, "tm_short": short, "loops_in_slab": slab}


def slab_runs(loops, bands, placed, margin: float = SLAB_MARGIN) -> list[str]:
    """The loops (a>b) that run sideways (more than _RUN_MIN) inside a band, more than `margin` from its faces (A18).
    Excused, at either end of a loop whose port lies inside the band: its first two segments (the C26 port rule: the
    stub out of the port and the turn after it), and its first run out of that port up to the length that clears
    the port's own element (its box grown by the router's clearance: 8W9Y's loop leaving an interfacial helix's end
    runs back above that helix to pass α8's end; membrane slab ruling); every other run in the slab counts."""
    from .route import _INFLATE

    def inside(y, m):
        return any(lo + m < y < hi - m for lo, hi in bands)

    def within(port, step, element) -> bool:  # does a run of `step` from the port stay over its element's body?
        if element is None:
            return False
        x0, _, x1, _ = element.rect
        reach = (x1 + _INFLATE - port[0]) if step > 0 else (port[0] - (x0 - _INFLATE))
        return abs(step) <= max(reach, 0.0) + 1e-9

    out = []
    for lp in loops:
        pts, n = lp.points, len(lp.points)
        segs = list(zip(pts, pts[1:]))
        runs = [k for k, ((xa, ya), (xb, yb)) in enumerate(segs) if abs(ya - yb) < 1e-9 and abs(xa - xb) > _RUN_MIN]
        start_in, end_in = inside(pts[0][1], 0.0), inside(pts[-1][1], 0.0)
        for k in runs:
            (xa, ya), (xb, _) = segs[k]
            if not inside(ya, margin):
                continue
            if (start_in and k <= 1) or (end_in and k >= n - 3):  # beside its own port (C26)
                continue
            if start_in and k == runs[0] and within(pts[0], xb - xa, placed.get(lp.a_id)):
                continue
            if end_in and k == runs[-1] and within(pts[-1], xa - xb, placed.get(lp.b_id)):
                continue
            out.append(f"{lp.a_id}>{lp.b_id}")
            break
    return out


# ------------------------------------------------------------------------------------------- all of them


def figure_checks(drawn, path, *, assembly: str = "auto", focus: str = "auto", chains=None) -> dict:
    """Every check on one drawn figure, measured once: a count (or the largest value) under each name, the findings
    under "details", the font files the text resolved to under "font". When that font is not INK_FONT the ink
    checks (INK_CHECKS) are not run: their counts are None and "skipped" says why."""
    page = measure(drawn)
    skip = page.font_reason
    collisions = None if skip else text_collisions(drawn, page=page)
    pairs = _copy_pairs(drawn)
    details = {
        "element_overlaps": element_overlaps(drawn.layout),
        "text_collisions": collisions,
        "off_page": None if skip else off_page(drawn, page=page),
        "label_repeats": label_repeats(drawn.layout),
        "copies_alike": {
            pair: {"deviation": dev, "matched": n, "of": m, "judged": why or "yes"}
            for pair, (dev, n, m, why) in pairs.items()
        },
        "domains": domain_checks(drawn, collisions or []) if drawn.layout.domains else {},
        "stem_name_gaps": None if skip else stem_name_gaps(drawn, page=page),
        "scope_silent": scope_honest(drawn, path, assembly, focus, chains=chains),
        "order": order_checks(drawn),
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
        "loop_crossings": loop_crossings(drawn.loops),
        "honest_crossings": honest_crossings(drawn.layout, drawn.loops),
        "font": page.fonts,
        "skipped": dict.fromkeys(INK_CHECKS, skip) if skip else {},
        "details": details,
    }
