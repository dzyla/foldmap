# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Text and marks beside the drawn nucleic-acid strands: contact beads kept apart along runs, a letter for each
single-stranded nucleotide, the "... N nt ..." labels of cropped gaps, isolated pairs (ties) as dotted lines or
labelled marks, and the page bounds that hold them all."""

from __future__ import annotations

import math

import numpy as np

from .layout import (
    _GLYPH_W,
    _LETTER_EM,
    NT_ROW,
    Layout,
    _footprint,
    _ink_points,
    _reach,
    _stem_label_box,
    glyph,
    label_boxes,
    line_clearance,
    text_em,
)
from .na_paths import _BEAD, NAPath, NATie, _dot, _near, _seg_dist, _seg_dists, _spaced, crosses_body, in_body, midway
from .route import clip_len, clip_lengths

_TIE_SHORT = 3.0  # stem widths: a longer tie is drawn as two labelled marks, not a line
_TICK_GAP, _TICK = 0.15, 0.35  # a tie mark: a tick this far out from its nucleotide, this long
_TIE_TURNS = (0.0, 45.0, -45.0, 90.0, -90.0, 135.0, -135.0, 180.0)  # degrees from the partner a mark may point
_TICK_REACH = (1.0, 2.0, 3.0)  # times _TICK a mark's tick is drawn, the shortest that leaves its label clear


_SLIDES = tuple(f * k for k in (0.1, 0.2, 0.3, 0.4, 0.48) for f in (1, -1))  # spacings a run's nucleotides may slide


def _unbunch(layout: Layout, paths: list[NAPath]) -> None:
    """Slide each run's nucleotides along it, as little as keeps them a bead's width (_BEAD) from every nucleotide of
    another strand or stem, where its own spacing puts one on top of another's (two runs crossing or passing close);
    if no slide does, the one that keeps them furthest. Never into its own stems' bodies."""
    at = {key: d for d in layout.stems for key in d.axial}
    taken = _Points(_BEAD)
    for d in layout.stems:
        for key in d.axial:
            taken.add(d.point(key))
    runs = [p for p in paths if p.kind == "run" and p.nts]
    for p in paths:
        if p.kind != "run":
            for xy in p.nt_xy:
                taken.add(xy)

    def room(xy, own) -> float:  # how near the nearest other nucleotide comes; -1 inside its own stems
        if any(in_body(q, d) for q in xy for d in own):
            return -1.0
        return min(taken.nearest(q) for q in xy)

    for p in runs:
        own = [at[key] for key in p.between]

        best, space = p.nt_xy, room(p.nt_xy, own)
        for shift in _SLIDES if space < _BEAD else ():
            xy = _spaced(p.points, len(p.nts), shift)
            r = room(xy, own)
            if r > space:
                best, space = xy, r
            if r >= _BEAD:
                break
        p.nt_xy = best
        for xy in best:
            taken.add(xy)


_LETTER_GAP = 0.38  # page units: closest two base letters (or a letter and a nucleotide) may sit


def _place_letters(layout: Layout, paths: list[NAPath]) -> dict[tuple[int, int], tuple[float, float]]:
    """A letter for each single-stranded nucleotide beside its strand: on the side away from the stems it leaves
    and enters (outside a hairpin loop, in the row of the stem's letters along a free end), else the other, else
    further out on either side, or a little along the strand, whichever first keeps clear of other letters,
    nucleotides, end labels, elements, stem names and the strands drawn (runs, ends, gaps: a spot a line passes
    through is not clear, random-100 C14); when none is clear, the one furthest from the other letters."""
    where = {key: d.point(key) for d in layout.stems for key in d.axial}
    lines = [np.asarray(p.points, float).reshape(-1, 2) for p in paths if len(p.points) > 1]
    strand_segs = np.vstack([np.hstack([q[:-1], q[1:]]) for q in lines]) if lines else np.empty((0, 4))
    seg_lo = np.minimum(strand_segs[:, :2], strand_segs[:, 2:])
    seg_hi = np.maximum(strand_segs[:, :2], strand_segs[:, 2:])
    em = _em(layout)
    line_clear = line_clearance(layout)  # figcheck.LINE_TEXT_CLEARANCE_PX at the page's text size
    ink = {b: _reach(_ink_points(b, False, _LETTER_EM * em, 0.0)) for b in "ACGUTIN"}

    def struck(q, key, clear: float | None = None) -> bool:  # a strand through the letter's ink, or near it
        base = layout.nucleic.strands[key[0]].bases[key[1]]
        x0, y0, x1, y1 = ink.get(base, ink["N"]) + np.array([q[0], q[1], q[0], q[1]])
        c = line_clear if clear is None else clear
        box = (x0 - c, y0 - c, x1 + c, y1 + c)
        near = (seg_lo[:, 0] < box[2]) & (seg_hi[:, 0] > box[0]) & (seg_lo[:, 1] < box[3]) & (seg_hi[:, 1] > box[1])
        return bool(near.any()) and bool((clip_lengths(strand_segs[near], [box]) > 0).any())

    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    taken = _Points(1.0)  # every letter, nucleotide and end label so far: nearer than 1 is measured exactly
    for d in layout.stems:
        for key in d.axial:
            taken.add(d.point(key, NT_ROW))
    index = {key: taken.add(xy) for key, xy in where.items()}  # a letter's own nucleotide does not count
    for _, _, xy in layout.na_ends:
        taken.add(xy)
    bodies = [*layout.placed.values(), *layout.stems]  # by their true (tilted) outlines, not their boxes
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    g = 1.2 * gh  # a stem's name stands a letter and a half clear of any base letter
    boxes = [np.add(_stem_label_box(d, gw / _GLYPH_W), (-g, -g, g, g)) for d in layout.stems]
    out: dict[tuple[int, int], tuple[float, float]] = {}

    def room(q, own: int) -> float:
        if any(x0 < q[0] < x1 and y0 < q[1] < y1 for x0, y0, x1, y1 in boxes) or any(_near(q, e) for e in bodies):
            return 0.0
        return taken.nearest(q, own)

    at = {key: d for d in layout.stems for key in d.axial}
    for p in paths:
        pts = np.asarray(p.points, float)
        stems = [at[key] for key in p.between]
        centre = np.mean([[d.cx, d.cy] for d in stems], axis=0) if stems else None
        for m, (i, xy) in enumerate(zip(p.nts, p.nt_xy)):
            xy = np.asarray(xy)
            ahead = np.asarray(p.nt_xy[min(m + 1, len(p.nt_xy) - 1)]) - np.asarray(p.nt_xy[max(m - 1, 0)])
            if len(p.nt_xy) < 2:
                ahead = pts[-1] - pts[0]
            ahead = ahead / max(float(np.hypot(*ahead)), 1e-9)
            side = np.array([-ahead[1], ahead[0]])
            if centre is not None and _dot(side, xy - centre) < 0:
                side = -side
            options = [xy + side * NT_ROW * f for f in (1.0, -1.0, 1.7, -1.7, 2.4, -2.4)]
            options += [xy + side * NT_ROW * f + ahead * t for f in (1.0, -1.0) for t in (0.4, -0.4)]
            key = (p.strand, i)
            spaces = [room(q, index[key]) for q in options]
            free = [r >= _LETTER_GAP and not struck(q, key) for q, r in zip(options, spaces)]
            if not any(free):  # crowded: further out, and further along the strand
                more = [xy + side * NT_ROW * f for f in (3.1, -3.1)]
                more += [xy + side * NT_ROW * f + ahead * t for f in (1.0, -1.0, 1.7, -1.7) for t in (0.8, -0.8)]
                options += more
                spaces += [room(q, index[key]) for q in more]
                free += [r >= _LETTER_GAP and not struck(q, key) for q, r in zip(more, spaces[-len(more) :])]
            if p.kind == "hairpin":  # a hairpin loop's letters stay in its stem's footprint (layout keeps it clear):
                # a spot inside clear with the line clearance, else one no strand runs through; else outside (8SA2)
                d = at[p.between[0]]
                fx0, fy0, fx1, fy1 = _footprint(d, gw / _GLYPH_W)
                g = gh / 2
                held = [
                    fx0 <= q[0] - d.cx - g
                    and q[0] - d.cx + g <= fx1
                    and fy0 <= q[1] - d.cy - g
                    and q[1] - d.cy + g <= fy1
                    for q in options
                ]
                inked = [h and r >= _LETTER_GAP and not struck(q, key, 0.0) for q, r, h in zip(options, spaces, held)]
                if any(ok and h for ok, h in zip(free, held)):  # a clear spot inside: only those
                    free = [ok and h for ok, h in zip(free, held)]
                elif any(inked):  # inside, no strand through its ink (within the clearance only)
                    free = inked
                elif not any(free) and any(held):  # none clear anywhere: the roomiest inside
                    spaces = [r if h else -1.0 for r, h in zip(spaces, held)]
            clear = [q for q, ok in zip(options, free) if ok]  # clear of letters and of every strand drawn
            spot = clear[0] if clear else options[int(np.argmax(spaces))]  # else the furthest from other letters
            out[(p.strand, i)] = (float(spot[0]), float(spot[1]))
            taken.add(spot)
    return out


_GAP_ALONG = tuple(0.5 + s * k / 20 for k in range(9) for s in ((1,) if k == 0 else (-1, 1)))  # middle out
_GAP_OUT = (1.0, 2.0, 3.0)  # how far beside its gap a label may sit, in multiples of the least that clears it


_GAP_TEXT = 0.75  # a gap label's text size, in ems of the figure's text (render draws it at 0.75 of the font)
_GAP_CHAR = 0.5  # its mean character width, in ems of its own text ("…", digits and spaces are narrow)
_GAP_ROOM = 1.15  # room for the page to grow (its bounds widen after routing) and for a larger font


def _gap_half(label: str, em: float) -> tuple[float, float]:
    """Half the width and height of a gap label's box, one em of the figure's text being `em` page units
    (measured: 0.48 and 0.43 em per character and 1.0 and 0.94 em high on 2NZ4 and the 4YBB fixture)."""
    size = _GAP_TEXT * em * _GAP_ROOM
    return _GAP_CHAR * size * len(label) / 2 + 0.05, size / 2


def _em(layout: Layout) -> float:
    """Page units of one em of the figure's text on this layout's page (layout.text_em, 1 unit of pad each side)."""
    x0, _, x1, _ = layout.bounds
    return text_em(x1 - x0 + 2.0)


def _gap_box(p: NAPath, em: float) -> tuple[float, float, float, float]:
    """The box a gap's label takes (one em of text `em` page units): centred at label_xy, else standing on the
    gap's midpoint."""
    hw, hh = _gap_half(p.label, em)
    if p.label_xy is None:
        x, y = midway(p.points)
        return (x - hw, y, x + hw, y + 2 * hh)
    x, y = p.label_xy
    return (x - hw, y - hh, x + hw, y + hh)


def _along(pts, f: float):
    """(the point a fraction f of the way along the polyline, the unit direction there)."""
    p = np.asarray(pts, float)
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    t = f * cum[-1]
    k = min(int(np.searchsorted(cum, t, side="right")) - 1, len(seg) - 1)
    u = (p[k + 1] - p[k]) / seg[k] if seg[k] > 0 else np.array([1.0, 0.0])
    return p[k] + (t - cum[k]) * u, u


def _label_gaps(layout: Layout, paths: list[NAPath], loops=()) -> None:
    """Put each labelled gap's "... N nt ..." beside its dashed path (NAPath.label_xy), as base letters are placed:
    at a fraction of the way along it (_GAP_ALONG, the middle first), on either side, then further out, the
    first spot whose box meets nothing: no element's or stem's body, strand or protein loop (`loops`, polylines),
    base letter, nucleotide (where its bead goes), end label, name or gap label placed before it, nor its own
    dashed line; when every spot meets something, the one meeting the fewest other things, then the least of its
    own line (a label beside its own gap reads as its own)."""
    gaps = [p for p in paths if p.kind == "gap" and p.label and len(p.points) > 1]
    if not gaps:
        return
    em = _em(layout)
    bodies = [*layout.placed.values(), *layout.ghosts, *layout.stems]
    segs = [(a, b) for p in paths for a, b in zip(p.points, p.points[1:])]
    owner = [k for k, p in enumerate(paths) for _ in p.points[1:]]  # the path each segment is part of
    segs += [(a, b) for pts in loops for a, b in zip(pts, pts[1:])]
    owner += [-1] * (len(segs) - len(owner))
    seg_box = np.array([(min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1])) for a, b in segs])
    seg_box = seg_box.reshape(-1, 4)
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    g = gw / 2 + 0.05, gh / 2
    glyphs = [d.point(key, NT_ROW) for d in layout.stems for key in d.axial] + list(layout.na_letters.values())
    glyphs += [xy for _, _, xy in layout.na_ends]
    taken = [(x - g[0], y - g[1], x + g[0], y + g[1]) for x, y in glyphs] + [b for _, b in label_boxes(layout, [])]
    beads = [d.point(key) for d in layout.stems for key in d.axial] + [q for p in paths for q in p.nt_xy]
    beads = np.asarray(beads, float).reshape(-1, 2)

    def meets(box, own: int) -> tuple[int, int]:  # (other things met, segments of its own gap met)
        x0, y0, x1, y1 = box
        edges = [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]
        n = sum(
            1
            for e in bodies
            if _near(((x0 + x1) / 2, (y0 + y1) / 2), e, 0.0)
            or x0 < e.cx < x1
            and y0 < e.cy < y1
            or any(crosses_body(a, b, e) for a, b in edges)
        )
        near = np.flatnonzero((seg_box[:, 0] < x1) & (x0 < seg_box[:, 2]) & (seg_box[:, 1] < y1) & (y0 < seg_box[:, 3]))
        hit = [k for k in near.tolist() if clip_len(*segs[k], box) > 1e-9]
        mine = sum(1 for k in hit if owner[k] == own)
        n += len(hit) - mine
        n += sum(1 for b in taken if b[0] < x1 and x0 < b[2] and b[1] < y1 and y0 < b[3])
        if len(beads):
            dx = beads[:, 0] - np.clip(beads[:, 0], x0, x1)
            dy = beads[:, 1] - np.clip(beads[:, 1], y0, y1)
            n += int(np.count_nonzero(dx * dx + dy * dy < (_BEAD / 2) ** 2))
        return n, mine

    for p in gaps:
        own = next(k for k, q in enumerate(paths) if q is p)
        hw, hh = _gap_half(p.label, em)
        best = None
        for out in _GAP_OUT:
            for f in _GAP_ALONG:
                q, u = _along(p.points, f)
                normal = np.array([-u[1], u[0]])
                reach = out * (abs(normal[0]) * hw + abs(normal[1]) * hh + 0.1)  # the box just clear of the gap
                for side in (1.0, -1.0):
                    x, y = q + side * reach * normal
                    box = (x - hw, y - hh, x + hw, y + hh)
                    n = meets(box, own)
                    if best is None or n < best[0]:
                        best = (n, (float(x), float(y)), box)
                    if n == (0, 0):
                        break
                if best[0] == (0, 0):
                    break
            if best[0] == (0, 0):
                break
        p.label_xy = best[1]
        taken.append(best[2])


def _tie_name(n: int) -> str:
    """a, b, ..., z, aa, ab, ...: the label of the n-th tie drawn as marks."""
    return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n)


def _ties(layout: Layout, paths: list[NAPath], loops=()) -> list[NATie]:
    """Each isolated pair as drawn: a straight dotted line when it is short (at most _TIE_SHORT stem widths) and
    crosses no stem, element or letter, otherwise a tick at each nucleotide with the same label (a, b, ... in
    the order of the ties), each tick turned from its partner as little as keeps it and its label clear of
    bodies, letters, nucleotides, labels, strands and protein loops (`loops`, polylines); where no turn does, the
    tick is drawn longer (_TICK_REACH) and the turns tried again (random-100 C2b-4: 4P9R's and 7Y7C's labels took
    the first turn, on a strand or loop); where none does at any length, the first turn."""
    na = layout.nucleic
    where = {key: d.point(key) for d in layout.stems for key in d.axial}
    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    letters = {key: d.point(key, NT_ROW) for d in layout.stems for key in d.axial}
    letters.update({key: np.asarray(xy) for key, xy in layout.na_letters.items()})
    bodies = [*layout.placed.values(), *layout.stems]
    width = max((d.width for d in layout.stems), default=1.0)
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    boxes = [b for _, b in label_boxes(layout, [])] + [
        (x - gw / 2, y - gh / 2, x + gw / 2, y + gh / 2) for _, _, (x, y) in layout.na_ends
    ]
    segments = [(a, b) for p in paths for a, b in zip(p.points, p.points[1:])]
    segments += [(a, b) for pts in loops for a, b in zip(pts, pts[1:])]
    seg = np.asarray(segments, float).reshape(-1, 4)
    seg_lo, seg_hi = np.minimum(seg[:, :2], seg[:, 2:]), np.maximum(seg[:, :2], seg[:, 2:])
    taken = [*letters.values(), *where.values()]
    hw, hh = gw / 2 + 0.05, gh / 2

    def near_segments(box, b) -> list:  # only these can meet the box or come within 0.15 of b (protein loops: many)
        lo = np.minimum(box[:2], np.subtract(b, 0.15))
        hi = np.maximum(box[2:], np.add(b, 0.15))
        rows = np.flatnonzero(((seg_lo <= hi) & (seg_hi >= lo)).all(axis=1))
        return [segments[k] for k in rows.tolist()]

    def clear(tick, at, own) -> bool:  # own: the tie's nucleotide, which the tick starts beside
        (a, b), (x, y) = tick, at
        box = (x - hw, y - hh, x + hw, y + hh)
        if any(_near(at, e, 0.0) or crosses_body(a, b, e) for e in bodies):
            return False
        pts = np.asarray(taken, float).reshape(-1, 2)  # vectorised (batch B: 37 s of 7Y7C's 140 s was here)
        near = pts[np.hypot(*(pts - np.asarray(own, float)).T) > 1e-9]
        if (np.hypot(*(near - np.asarray(at, float)).T) < _LETTER_GAP + hh).any():
            return False
        if (_seg_dists(near, a, b) < 0.2).any():
            return False
        if any(min(box[2], o[2]) > max(box[0], o[0]) and min(box[3], o[3]) > max(box[1], o[1]) for o in boxes):
            return False
        return not any(clip_len(p, q, box) > 1e-9 or _seg_dist(b, p, q) < 0.15 for p, q in near_segments(box, b))

    out: list[NATie] = []
    for p in na.structure.ties:
        ka, kb = (p.a, p.i), (p.b, p.j)
        if ka not in where or kb not in where:
            continue  # cropped away
        pa, pb = where[ka], where[kb]
        tie = NATie(ka, kb, [tuple(map(float, pa)), tuple(map(float, pb))])
        long = math.dist(pa, pb) > _TIE_SHORT * width
        hit = any(crosses_body(pa, pb, e) for e in bodies) or any(
            _seg_dist(xy, pa, pb) < gh / 2 for key, xy in letters.items() if key not in (ka, kb)
        )
        if long or hit:
            tie.label = _tie_name(sum(1 for t in out if t.label))
            for xy, other in ((pa, pb), (pb, pa)):
                toward = np.subtract(other, xy)
                toward = toward / max(float(np.hypot(*toward)), 1e-9)
                options = []
                for reach in _TICK_REACH:
                    for deg in _TIE_TURNS:
                        c, s_ = math.cos(math.radians(deg)), math.sin(math.radians(deg))
                        u = np.array([toward[0] * c - toward[1] * s_, toward[0] * s_ + toward[1] * c])
                        a, b = xy + u * _TICK_GAP, xy + u * (_TICK_GAP + _TICK * reach)
                        options.append(((a, b), b + u * (hh + 0.1)))
                (a, b), at = next((o for o in options if clear(*o, xy)), options[0])
                tie.marks.append(((float(a[0]), float(a[1])), (float(b[0]), float(b[1])), (float(at[0]), float(at[1]))))
                boxes.append((at[0] - hw, at[1] - hh, at[0] + hw, at[1] + hh))
                taken.append(np.asarray(at))
        out.append(tie)
    return out


class _Points:
    """Points in a cell hash of cell size `reach`: the nearest one within reach is found among the 3 x 3 cells
    around a query, so placing n letters is linear in n, not quadratic."""

    def __init__(self, reach: float):
        self.reach = reach
        self.xy: list[tuple[float, float]] = []
        self.cells: dict[tuple[int, int], list[int]] = {}

    def _cell(self, x: float, y: float) -> tuple[int, int]:
        return math.floor(x / self.reach), math.floor(y / self.reach)

    def add(self, q) -> int:
        """Add a point; its index."""
        x, y = float(q[0]), float(q[1])
        self.xy.append((x, y))
        self.cells.setdefault(self._cell(x, y), []).append(len(self.xy) - 1)
        return len(self.xy) - 1

    def nearest(self, q, skip: int = -1) -> float:
        """Distance from q to the nearest point but `skip`, or `reach` when none is nearer."""
        x, y = float(q[0]), float(q[1])
        cx, cy = self._cell(x, y)
        best = self.reach
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for k in self.cells.get((cx + dx, cy + dy), ()):
                    if k != skip:
                        tx, ty = self.xy[k]
                        best = min(best, float(np.hypot(x - tx, y - ty)))
        return best


def _widen(layout: Layout, paths: list[NAPath]) -> None:
    """Grow layout.bounds to hold the strands, every base letter, the end labels, the gaps' labels and the
    ties' labels."""
    gw, gh = glyph(layout)  # a letter's room, sized from the text drawn
    pad = NT_ROW + max(gw, gh)
    boxes = [(x - pad, y - pad, x + pad, y + pad) for p in paths for x, y in p.points]
    g = max(gw, gh) / 2 + 0.05  # a letter turned any way
    letters = [d.point(key, NT_ROW) for d in layout.stems for key in d.axial] + list(layout.na_letters.values())
    letters += [xy for _, _, xy in layout.na_ends]
    boxes += [(x - g, y - g, x + g, y + g) for x, y in letters]
    boxes += [_gap_box(p, _em(layout)) for p in paths if p.label]
    boxes += [b for name, b in label_boxes(layout, []) if name.startswith("tie-label:")]
    if not boxes:
        return
    b = np.asarray(boxes, float)
    x0, y0, x1, y1 = layout.bounds
    layout.bounds = (
        min(x0, float(b[:, 0].min())),
        min(y0, float(b[:, 1].min())),
        max(x1, float(b[:, 2].max())),
        max(y1, float(b[:, 3].max())),
    )
