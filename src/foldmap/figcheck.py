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

import numpy as np
from matplotlib.path import Path as MPath
from matplotlib.transforms import Affine2D, Bbox

from .layout import joined_pieces, lying_over

# Lines a text can be crossed by (render's gids, matched by prefix). A loop's casing (loopcase:) is drawn in the
# page colour and a chain-end stub belongs to its terminus letter: neither is read.
LINE_GIDS = ("loop:", "loop-seg:", "loop-arrow:", "na:run:", "na:tie:", "ligand-tether:", "disulfide:", "hbond:")
GLYCAN_STEM = re.compile(r"glycan:[^:]+:stem$")
# Bodies: opaque shapes a text must not sit on unless it is theirs, and that hide a line drawn under them. Element
# bodies are owned by their element id (label:<id>, resnum:<id>:<end>, ghostlabel:<id> are theirs).
ELEMENT_BODIES = ("strand:", "helix:", "helix-back:", "eta:", "ghost:")
OWN_BODIES = {"ligand:": "ligand-label:", "linked:": "linked-label:", "domain-name-plate:": "domain-label:"}
TEXT_OVERLAP_PX2 = 1.0  # text over text or a body: overlaps up to this many square pixels are anti-aliasing
OFF_PAGE_PX = 1.0  # text may reach this far past the page edge
PANEL_GRID = 0.05  # page units: the raster domain panel outlines are measured on
_STEP_PX = 0.5  # sampling step (pixels) along a line inside a text box, and across a text box over a body
_AUTO_DOMAIN = re.compile(r"D\d+( · .+)?")  # features.auto_domains names: D1, D2 ... with " · chain" when split


# ------------------------------------------------------------------------------------------- loop crossings


def loop_crossings(loops) -> int:
    """Proper crossings between the polylines of different loops (shared end points do not count)."""

    def cross(p, q, r, s) -> bool:
        def side(a, b, c):
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

        return side(r, s, p) * side(r, s, q) < 0 and side(p, q, r) * side(p, q, s) < 0

    segs = [list(zip(lp.points, lp.points[1:])) for lp in loops]
    return sum(
        cross(p, q, r, s)
        for i in range(len(segs))
        for j in range(i + 1, len(segs))
        for p, q in segs[i]
        for r, s in segs[j]
    )


def honest_crossings(lay, loops) -> int:
    """Every crossing a reader sees: loop x loop (loop_crossings), plus each loop drawn through the body of an element
    it does not join, a helix drawn over a sheet included (a loop passing under a helix not its own), except a loop
    leaving a helix drawn over a sheet across the strand its port lies on (Task 14 review)."""
    from .route import _crosses_body

    n = loop_crossings(loops)
    for lp in loops:
        segs = list(zip(lp.points, lp.points[1:]))
        a, b = lay.placed.get(lp.a_id), lay.placed.get(lp.b_id)
        ports = [p for q, p in ((a, lp.points[0]), (b, lp.points[-1])) if q is not None and q.over]
        for k, e in lay.placed.items():
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


def _renderer(fig):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas.get_renderer()


def _texts(fig) -> list:
    return [t for ax in fig.axes for t in ax.texts if t.get_visible() and t.get_text().strip()] + [
        t for t in fig.texts if t.get_visible() and t.get_text().strip()
    ]


def _name(t) -> str:
    """A text's gid, or its words for a text drawn without one (the title, legend entries)."""
    return t.get_gid() or f"text:{t.get_text()!r}"


def _ink_boxes(fig, renderer=None) -> list[tuple[object, Bbox]]:
    """(text, ink box in pixels) for every drawn text: the extents of its glyph outlines (matplotlib's TextToPath,
    placed and turned as the text is drawn), tighter than the layout box get_window_extent gives (line height,
    descent)."""
    from matplotlib.textpath import TextToPath

    renderer = renderer or _renderer(fig)
    to_path = TextToPath()
    got: list[Bbox] = []

    def record(gc, x, y, s, prop, angle, ismath=False, mtext=None):
        if not s.strip():
            return
        verts, codes = to_path.get_text_path(prop, s, ismath=ismath)
        if len(verts):
            scale = renderer.points_to_pixels(prop.get_size_in_points()) / to_path.FONT_SCALE
            if renderer.flipy():  # as RendererBase._draw_text_as_path places them: y counted from the top
                y = renderer.get_canvas_width_height()[1] - y
            place = Affine2D().scale(scale).rotate_deg(angle).translate(x, y)
            got.append(MPath(verts, codes).transformed(place).get_extents())

    renderer.draw_text = record  # this instance only: the glyph outlines are measured, nothing is drawn
    out = []
    try:
        for t in _texts(fig):
            got.clear()
            t.draw(renderer)
            if got:
                out.append((t, Bbox.union(got)))
    finally:
        del renderer.draw_text
    return out


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


class _Page:
    """The figure as drawn, for the text checks: texts with their ink boxes, bodies and lines in pixels."""

    def __init__(self, drawn):
        fig = drawn.figure
        self.fig = fig
        self.texts = _ink_boxes(fig)
        self.ax = fig.axes[0] if fig.axes else None
        self.layout = getattr(drawn, "layout", None)
        artists = [a for ax in fig.axes for a in (*ax.patches, *ax.lines)]
        self.bodies = []  # (gid, path in pixels, zorder, opaque)
        self.lines = []  # (gid, [polyline (n, 2) in pixels], zorder, half its stroke width in pixels)
        for a in artists:
            gid = a.get_gid() or ""
            if gid.startswith(LINE_GIDS) or GLYCAN_STEM.match(gid):
                self.lines.append((gid, _polylines(a), a.get_zorder(), a.get_linewidth() * fig.dpi / 72 / 2))
            elif gid.startswith(ELEMENT_BODIES) or gid.startswith(tuple(OWN_BODIES)) or gid.startswith("na:stem:"):
                path = a.get_path().transformed(a.get_transform())
                self.bodies.append((gid, path, a.get_zorder(), _opaque(a)))
        self.stems = self._stem_bodies()

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


def text_collisions(drawn) -> list[str]:
    """Every "<text gid> × <thing gid>" a reader sees on the page (a text drawn without a gid is named by its words):
    - text over text, their ink boxes overlapping by more than TEXT_OVERLAP_PX2 (a text's own parts, one gid, never);
    - text over a body (an element, a stem, a ligand or linked mark, a domain name plate) it does not belong to (its
      label, its residue numbers, its letters, its name), by more than TEXT_OVERLAP_PX2; a helix drawn over a sheet
      and that sheet's strands are each other's (layout.lying_over);
    - text crossed by a line (LINE_GIDS: loops, nucleic-acid runs and ties, tethers, disulfide bars, glycan stems,
      H-bonds) whose stroke meets its ink box where no opaque body drawn above the line hides it there: a line
      under the text still runs through its letters, a loop under the strand its letter sits on does not show. A
      text's own line (a loop's length, a gap's label, a tie's letter, a ligand's tethers, a strand's letters and
      ends) is not counted."""
    page = _Page(drawn)
    lay = page.layout
    over = {frozenset(f"el:{k}" for k in pair) for pair in lying_over(lay)} if lay is not None else set()
    own_lines = page.own_lines()
    hits: list[str] = []
    named = [(t, _name(t), box) for t, box in page.texts]
    lines = [  # each line, the extent of its drawing, half its stroke and the opaque bodies drawn above it
        (gid, polys, Bbox.union([Bbox.from_extents(*p.min(axis=0), *p.max(axis=0)) for p in polys]), pad, covers)
        for gid, polys, z, pad in page.lines
        if polys
        for covers in [[p for _, p, bz, opaque in page.bodies if opaque and bz > z]]
    ]
    for (_, ga, a), (_, gb, b) in combinations(named, 2):
        if ga != gb and _overlap(a, b) > TEXT_OVERLAP_PX2:
            hits.append(f"{ga} × {gb}")
    for _, gid, box in named:
        mine = _owner(gid)
        for body, path, _, _ in page.bodies:
            if body.startswith("na:stem:"):
                continue  # stems: their whole body below
            theirs = _owner(body)
            if theirs is not None and (theirs == mine or frozenset((theirs, mine)) in over):
                continue
            if _overlap(box, path.get_extents()) > TEXT_OVERLAP_PX2 and _area_on(box, path) > TEXT_OVERLAP_PX2:
                hits.append(f"{gid} × {body}")
        for body, path, own in page.stems:
            if gid in own:
                continue
            if _overlap(box, path.get_extents()) > TEXT_OVERLAP_PX2 and _area_on(box, path) > TEXT_OVERLAP_PX2:
                hits.append(f"{gid} × {body}")
        for line, polys, extent, pad, covers in lines:
            reach = Bbox.from_extents(box.x0 - pad, box.y0 - pad, box.x1 + pad, box.y1 + pad)  # its stroke meets ink
            if line in own_lines.get(gid, ()) or not extent.overlaps(reach):
                continue
            if _crosses(polys, reach, covers):
                hits.append(f"{gid} × {line}")
    return hits


def off_page(drawn) -> list[str]:
    """Texts (their ink) and legend marks reaching more than OFF_PAGE_PX past the page edge."""
    fig = drawn.figure
    page = fig.bbox
    renderer = _renderer(fig)
    out = [_name(t) for t, box in _ink_boxes(fig, renderer) if _beyond(box, page)]
    for ax in fig.axes:
        for a in (*ax.patches, *ax.lines):
            gid = a.get_gid() or ""
            if gid.startswith("legend") and _beyond(a.get_window_extent(renderer), page):
                out.append(gid)
    return out


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


def copies_alike(drawn) -> dict[tuple[str, str], float]:
    """(copy chain, reference chain) -> the largest element-centre deviation (page units) between the two after the
    best rigid 2D motion (fit_deviation; a reflection only for a mirror tie), for every tie the symmetry asked for
    (Layout.copy_ties). Elements are matched as one element through the copy classes (model.copy_classes over
    model.homolog_chains): elements one copy lacks are left out. A tied copy with no element matched is inf."""
    from .model import copy_classes, homolog_chains

    lay = drawn.layout
    mate = copy_classes(drawn.sses, homolog_chains(drawn.bb)).mate
    by_chain: dict[str, dict[str, tuple[float, float]]] = {}
    for k, p in lay.placed.items():
        by_chain.setdefault(p.sse.chain, {})[mate.get(k, k)] = (p.cx, p.cy)
    out = {}
    for copy, (ref, mirrored) in sorted(lay.copy_ties.items()):
        a, b = by_chain.get(ref, {}), by_chain.get(copy, {})
        keys = sorted(a.keys() & b.keys())
        out[(copy, ref)] = fit_deviation([a[k] for k in keys], [b[k] for k in keys], mirrored) if keys else math.inf
    return out


# ------------------------------------------------------------------------------------------- domains


def outline_overlap(a: MPath, b: MPath, grid: float = PANEL_GRID) -> float:
    """The area two closed outlines (page units) share, over the smaller one's area, both rasterised at `grid`."""
    area_a, area_b = _raster_area(a, a.get_extents(), grid), _raster_area(b, b.get_extents(), grid)
    if not area_a or not area_b:
        return 0.0
    ea, eb = a.get_extents(), b.get_extents()
    both = Bbox.from_extents(max(ea.x0, eb.x0), max(ea.y0, eb.y0), min(ea.x1, eb.x1), min(ea.y1, eb.y1))
    if both.x1 <= both.x0 or both.y1 <= both.y0:
        return 0.0
    pts = _grid(both, grid)
    shared = float((a.contains_points(pts) & b.contains_points(pts)).sum()) * grid * grid
    return shared / min(area_a, area_b)


def _grid(box: Bbox, step: float) -> np.ndarray:
    """Cell centres of a grid of `step`, anchored at the origin so every box shares one raster."""
    xs = (np.arange(math.floor(box.x0 / step), math.ceil(box.x1 / step)) + 0.5) * step
    ys = (np.arange(math.floor(box.y0 / step), math.ceil(box.y1 / step)) + 0.5) * step
    gx, gy = np.meshgrid(xs, ys)
    return np.c_[gx.ravel(), gy.ravel()]


def _raster_area(path: MPath, box: Bbox, grid: float) -> float:
    return float(path.contains_points(_grid(box, grid)).sum()) * grid * grid


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
    for (na, pa), (nb, pb) in combinations(outlines.items(), 2):
        f = outline_overlap(pa, pb)
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


def stem_name_gaps(drawn) -> dict[str, tuple[float, float]]:
    """na:label:<k> -> (gap from the name's ink box to the ink of stem k's own letters (na:nt: of its nucleotides),
    or to the stem's rect when its letters are not drawn; gap to the nearest other stem's rect), page units."""
    lay, fig = drawn.layout, drawn.figure
    if not lay.stems or not fig.axes:
        return {}
    to_page = fig.axes[0].transData.inverted()
    ink = {_name(t): tuple(to_page.transform_bbox(box).extents) for t, box in _ink_boxes(fig)}
    na = lay.nucleic
    out = {}
    for k, d in enumerate(lay.stems):
        name = ink.get(f"na:label:{k}")
        if name is None:
            continue
        letters = [ink.get(f"na:nt:{na.strands[s].labels[i].chain}:{na.strands[s].labels[i].tag}") for s, i in d.axial]
        letters = [b for b in letters if b is not None]
        own = min(_box_gap(name, b) for b in letters) if letters else _box_gap(name, d.rect)
        others = [_box_gap(name, e.rect) for j, e in enumerate(lay.stems) if j != k]
        out[f"na:label:{k}"] = (own, min(others) if others else math.inf)
    return out


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
    among `chains` (--chains), and anything left out by an explicit --focus. A nucleic acid cropped by na_extent auto
    is not counted until that crop is announced (plan 10 Task 8)."""
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
    covered: set[str] = set()
    if scope["reason"]:
        if several:
            covered |= set(expected) - {c.name for c in model[0]}
        if focus == "auto" and lay.focus_chains and lay.subunits:
            covered |= set(expected) - set(lay.focus_chains)
    chosen = focus not in ("auto", "none")
    out = []
    for chain in sorted(expected):
        if chain in drawn_chains or chain in covered or chosen or (chains and chain not in chains):
            continue
        if lay.na_cropped and chain in _nucleic_chains(read):
            continue
        out.append(f"chain {chain}")
    if chosen or lay.na_cropped:
        return out
    for chain in sorted(shown):
        drawn_nts = shown[chain]
        run: list[tuple[int, str]] = []
        for r in [*expected.get(chain, []), None]:
            if r is not None and r not in drawn_nts:
                run.append(r)
                continue
            if len(run) >= 4:
                out.append(f"nt {chain}:{run[0][0]}{run[0][1]}-{run[-1][0]}{run[-1][1]}")
            run = []
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
      bridges, side by side in that sheet's order (either way round; a closed barrel from any strand), each two
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
    for k, ids in enumerate(lay.sheets):
        sheet = sheets.get(frozenset(ids))
        if sheet is None:
            out.append(f"sheet {k}: {', '.join(ids)} is no sheet of the bridges")
            continue
        want = [s.id for s in sheet.strands]
        got = _drawn_order(lay, ids)
        if not _same_order(got, want, sheet.closed):
            out.append(f"sheet {k}: drawn {', '.join(got)}, paired {', '.join(want)}")
        for (a, da), (b, db) in zip(zip(want, sheet.directions), zip(want[1:], sheet.directions[1:])):
            pa, pb = lay.placed[a], lay.placed[b]
            same = float(pa.direction @ pb.direction) > 0
            if same != (da == db):
                out.append(f"sheet {k}: {a} and {b} drawn {'parallel' if same else 'antiparallel'}")
    return out


def _drawn_order(lay, ids: list[str]) -> list[str]:
    """Strand ids in the order they stand side by side on the page: across their mean direction."""
    ps = [lay.placed[k] for k in ids]
    ref = ps[0].direction
    mean = sum(p.direction * (1.0 if p.direction @ ref >= 0 else -1.0) for p in ps)
    across = np.array([mean[1], -mean[0]])
    return [k for _, k in sorted(((float(np.array([p.cx, p.cy]) @ across), k) for k, p in zip(ids, ps)))]


def _same_order(got: list[str], want: list[str], closed: bool) -> bool:
    if closed:
        n = len(want)
        turns = [want[k:] + want[:k] for k in range(n)]
        return any(got in (t, t[::-1]) for t in turns)
    return got in (want, want[::-1])


# ------------------------------------------------------------------------------------------- all of them


def figure_checks(drawn, path, *, assembly: str = "auto", focus: str = "auto", chains=None) -> dict:
    """Every check on one drawn figure: a count (or the largest value) under each name, the findings under
    "details"."""
    collisions = text_collisions(drawn)
    details = {
        "element_overlaps": element_overlaps(drawn.layout),
        "text_collisions": collisions,
        "off_page": off_page(drawn),
        "label_repeats": label_repeats(drawn.layout),
        "copies_alike": copies_alike(drawn),
        "domains": domain_checks(drawn, collisions) if drawn.layout.domains else {},
        "stem_name_gaps": stem_name_gaps(drawn),
        "scope_silent": scope_honest(drawn, path, assembly, focus, chains=chains),
        "order": order_checks(drawn),
    }
    domains = details["domains"]
    return {
        "element_overlaps": len(details["element_overlaps"]),
        "text_collisions": len(details["text_collisions"]),
        "off_page": len(details["off_page"]),
        "label_repeats": sum(len(v) for v in details["label_repeats"].values()),
        "copies_alike": max(details["copies_alike"].values(), default=0.0),
        "domain_rules": sum(len(domains.get(k, [])) for k in ("cross_chain", "unalike", "small")),
        "domain_overlap": max((f for *_, f in domains.get("panel_overlaps", [])), default=0.0),
        "domain_name_collisions": len(domains.get("name_collisions", [])),
        "stem_name_gap": max((own for own, _ in details["stem_name_gaps"].values()), default=0.0),
        "scope_silent": len(details["scope_silent"]),
        "order": len(details["order"]),
        "loop_crossings": loop_crossings(drawn.loops),
        "honest_crossings": honest_crossings(drawn.layout, drawn.loops),
        "details": details,
    }
