# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Nucleic-acid strands on the page: single strands routed like loops, hairpin arcs, free ends, cropped gaps."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

import numpy as np

from .layout import _GLYPH_H, _GLYPH_W, NT_ROW, NT_STEP, Layout, PlacedStem, _stem_label_box, label_boxes
from .route import _clear, _curve, _element_obstacles, _Grid, _Router, _simplify

log = logging.getLogger(__name__)


_BOW_MAX = 8  # internal-loop sides up to this many nucleotides are drawn as bows, not routed
SHORT_GAP = 3  # a cropped gap of fewer nucleotides is drawn as a bare break, without its "… N nt …" label
_JOINT_STRAIGHT = 1.5  # a joint between two pieces of a bent duplex shorter than this may be drawn straight


@dataclass(eq=False)
class NAPath:
    """A drawn stretch of single strand. kind: run (stem to stem), hairpin (an arc over a stem's tip), joint
    (between two pieces of a bent duplex), tail (a strand's free end), trace (a strand with no pairs) or gap
    (nucleotides cropped away, dashed; labelled "… N nt …" unless fewer than SHORT_GAP). nts: the nucleotide
    indices drawn on it, 5'->3', at nt_xy along `points`. fallback: no drawing clear of the elements was found.
    between: the (strand, nucleotide) of the stem nucleotides it leaves and enters (one for a tail)."""

    kind: str
    strand: int
    nts: list[int]
    points: list[tuple[float, float]]
    nt_xy: list[tuple[float, float]]
    label: str = ""
    fallback: bool = False
    between: tuple = ()


def path_name(path: NAPath, na) -> str:
    """How a warning names a drawn strand: "strand <chain>:<first>-<last>" by the nucleotides on it (or, with
    none, the stem nucleotides it joins; a gap by the strand piece before it)."""
    labels = na.strands[path.strand].labels
    idx = path.nts or [i for _, i in path.between] or [len(labels) - 1]
    first, last = labels[idx[0]], labels[idx[-1]]
    return f"strand {first.chain}:{first.tag}-{last.tag}"


def _use(grid: _Grid, pts) -> None:
    """Mark the cells under a drawn strand as used, so later routes avoid crossing it (at the reuse cost)."""
    cells = sorted({c for a, b in zip(pts, pts[1:]) for c in grid.cells_along(a, b)})
    for c in cells:
        grid.used[c] += 1


def _dot(a, b) -> float:
    """a . b of two 2-vectors, written out: numpy's dot goes through BLAS, whose kernels differ in the last bit."""
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1])


def _xy(pts) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in pts]


def _length(pts) -> float:
    p = np.asarray(pts, float)
    return float(np.hypot(*np.diff(p, axis=0).T).sum()) if len(p) > 1 else 0.0


def _spaced(pts, m: int) -> list[tuple[float, float]]:
    """m points evenly spaced along the polyline, ends excluded."""
    p = np.asarray(pts, float)
    if m <= 0 or len(p) < 2:
        return []
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    out = []
    for t in np.arange(1, m + 1) / (m + 1) * cum[-1]:
        k = min(int(np.searchsorted(cum, t, side="right")) - 1, len(seg) - 1)
        f = 0.0 if seg[k] == 0 else (t - cum[k]) / seg[k]
        q = p[k] + f * (p[k + 1] - p[k])
        out.append((float(q[0]), float(q[1])))
    return out


_NA_CORNER = 0.6  # radius of a routed strand's rounded corners


def _rounded(pts, radius: float = _NA_CORNER, n: int = 8) -> list[tuple[float, float]]:
    """The polyline with each corner replaced by a quadratic curve (as a protein loop's corners are drawn),
    sampled, so the nucleotides spaced along it sit on the drawn line."""
    p = [np.asarray(q, float) for q in pts]
    out = [p[0]]
    t = np.linspace(0.0, 1.0, n)[:, None]
    for a, b, c in zip(p, p[1:], p[2:]):
        la, lc = float(np.hypot(*(b - a))), float(np.hypot(*(c - b)))
        r = min(radius, la / 2, lc / 2)
        if r < 1e-9:
            out.append(b)
            continue
        q0, q2 = b - (b - a) / la * r, b + (c - b) / lc * r
        out += list((1 - t) ** 2 * q0 + 2 * t * (1 - t) * b + t**2 * q2)
    out.append(p[-1])
    return _xy(_simplify(out))


def midway(pts) -> tuple[float, float]:
    """The point half way along a polyline (where a gap's label goes)."""
    return _spaced(pts, 1)[0] if len(pts) > 1 else tuple(pts[0])


def _bow(pa, pb, need: float, n: int = 32, out=None) -> list[tuple[float, float]]:
    """A curve from pa to pb bowed out (to the left, or towards `out`) until it is about `need` long."""
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    d = pb - pa
    chord = max(float(np.hypot(*d)), 1e-9)
    left = np.array([-d[1], d[0]]) / chord
    if out is not None and _dot(left, out) < 0:
        left = -left
    h = 0.0
    while True:
        ctrl = (pa + pb) / 2 + left * h
        t = np.linspace(0, 1, n)[:, None]
        pts = (1 - t) ** 2 * pa + 2 * t * (1 - t) * ctrl + t**2 * pb
        if _length(pts) >= need or h > 4 * need:
            return _xy(pts)
        h += 0.25


def _arc(d: PlacedStem, ka, kb, m: int) -> list[tuple[float, float]]:
    """A hairpin loop: a near-circle beyond the stem's tip from nucleotide ka to kb, room for m nucleotides."""
    pa, pb = d.point(ka), d.point(kb)
    out = d.direction if d.axial[ka] >= 0 else -d.direction
    centre, r = d.hairpin_circle(ka, kb, m)
    a0 = np.arctan2(*(pa - centre)[::-1])
    a1 = np.arctan2(*(pb - centre)[::-1])
    cands = [a1 + 2 * np.pi * k for k in (-1, 0, 1) if 0 < abs(a1 + 2 * np.pi * k - a0) < 2 * np.pi]
    a1 = max(cands, key=lambda c: _dot((np.cos((a0 + c) / 2), np.sin((a0 + c) / 2)), out))
    t = np.linspace(a0, a1, max(24, 4 * (m + 2)))
    return _xy(centre + r * np.c_[np.cos(t), np.sin(t)])


_TURNS = (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4)  # a free end's bends


def _free_end(d: PlacedStem, key, m: int, five: bool, turn: float = 0.0):
    """A strand's free end of m nucleotides, straight on from the stem and away from it (after its first
    nucleotide turned by `turn` radians): (points 5'->3', nucleotide positions 5'->3', where the end label goes)."""
    _, out = d.port(key, leaving=not five)
    rail = d.point(key)
    c, s = math.cos(turn), math.sin(turn)
    bent = np.array([out[0] * c - out[1] * s, out[0] * s + out[1] * c])
    first = rail + out * (NT_STEP / 2 + NT_STEP)
    xy = [first + bent * (k - 1) * NT_STEP for k in range(1, m + 1)]  # nearest the stem first
    tip = first + bent * (m * NT_STEP + NT_ROW)
    pts = [rail] + xy
    if five:
        pts, xy = pts[::-1], xy[::-1]
    return _xy(pts), _xy(xy), (float(tip[0]), float(tip[1]))


def _tail(router: _Router, layout: Layout, d: PlacedStem, s: int, key, nts: list[int], five: bool, labels):
    """A free end drawn straight on, or bent after its first nucleotide when that is what keeps it (and its end
    label) clear of the elements, and its label clear of loops and other text: (path, end label position,
    outward direction at its tip)."""
    found = None
    for need_label in (True, False) if labels is not None else (False,):
        for turn in _TURNS:
            pts, xy, tip = _free_end(d, key, len(nts), five, turn)
            drawn = pts[::-1] if five else pts
            if _na_clear([*drawn, tip], (d,), router.elements) and (
                not need_label or _label_room(router, layout, tip, labels)
            ):
                found = pts, xy, tip
                break
        if found:
            break
    pts, xy, tip = found or _free_end(d, key, len(nts), five)
    path = NAPath("tail", s, nts, pts, xy, fallback=found is None, between=(key,))
    far = xy[0] if five else xy[-1]
    v = np.subtract(tip, far)
    return path, tip, v / max(float(np.hypot(*v)), 1e-9)


_END_SPOTS = ((1.0, 0.0), (1.0, 0.6), (2.0, 0.0), (2.0, 0.6), (1.0, -0.6), (3.0, 0.0), (0.0, 2.4), (3.0, 1.2))


def _face_label(router: _Router, layout: Layout, d: PlacedStem, key, leaving: bool, labels) -> tuple:
    """Where the end label of a strand ending at a stem goes: just beyond its rail's end, else a little out to
    the side or further on, the first spot clear of loops, elements, stems and other text."""
    port, along = d.port(key, leaving)
    side = d.normal * d.side[key]  # away from the stem's axis
    spots = [port + along * NT_ROW * a + side * _GLYPH_H * b for a, b in _END_SPOTS]
    clear = (q for q in spots if labels is not None and _label_room(router, layout, q, labels))
    spot = next(clear, spots[0])
    return float(spot[0]), float(spot[1])


def _label_room(router: _Router, layout: Layout, xy, labels: list) -> bool:
    """An end label at xy keeps clear of the loops routed so far (their grid cells), every element's and stem's
    body, the text and stub boxes, and the end labels already placed."""
    hw, hh = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2
    x0, y0, x1, y1 = xy[0] - hw, xy[1] - hh, xy[0] + hw, xy[1] + hh
    grid = router.grid
    (i0, j0), (i1, j1) = grid.cell((x0, y0)), grid.cell((x1, y1))
    i0, j0, i1, j1 = max(i0, 0), max(j0, 0), min(i1, grid.nx - 1), min(j1, grid.ny - 1)
    if i0 <= i1 and j0 <= j1 and grid.used[i0 : i1 + 1, j0 : j1 + 1].any():
        return False
    if any(bx0 < x1 and x0 < bx1 and by0 < y1 and y0 < by1 for bx0, by0, bx1, by1 in router.keep_out):
        return False
    if any(abs(x - xy[0]) < 2 * hw and abs(y - xy[1]) < 2 * hh for x, y in labels):
        return False
    return not any(_near(xy, e, hh + 0.05) for e in [*layout.placed.values(), *layout.stems])


def _face(d: PlacedStem, key, leaving: bool) -> tuple[np.ndarray, np.ndarray]:
    """Where the strand crosses the stem's end face at this (end) nucleotide, along its rail, and the outward
    direction there: the port a strand is routed from, on the stem's outline as an element's port is."""
    _, out = d.port(key, leaving)
    reach = max(d.length / 2 - d.axial[key] * _dot(out, d.direction), 0.0)
    return d.point(key) + out * reach, out


def _local(pt, d) -> tuple[float, float]:
    c, s = math.cos(d.angle), math.sin(d.angle)
    x, y = pt[0] - d.cx, pt[1] - d.cy
    return x * c + y * s, y * c - x * s


def _in_body(pt, d) -> bool:
    u, v = _local(pt, d)
    return abs(u) < d.length / 2 - 1e-9 and abs(v) < d.width / 2 - 1e-9


def _near(pt, e, margin: float = _GLYPH_W / 2) -> bool:
    """A letter at pt would touch the element's (or stem's) body."""
    u, v = _local(pt, e)
    return abs(u) < e.length / 2 + margin and abs(v) < e.width / 2 + margin


def _leaves_at(p, q, d) -> float:
    """The fraction of the way along p-q (p inside d's body) where it leaves the body."""
    (u0, v0), (u1, v1) = _local(p, d), _local(q, d)
    t = 1.0
    for a, b, h in ((u0, u1, d.length / 2), (v0, v1, d.width / 2)):
        if b > a:
            t = min(t, (h - a) / (b - a))
        elif b < a:
            t = min(t, (-h - a) / (b - a))
    return max(0.0, t)


def _outside(pts, own) -> list[tuple[float, float]]:
    """The polyline without its first and last stretches inside its own stems' bodies (a strand runs along its
    rail to the stem's end face before it is drawn among the other elements)."""

    def trim(p):
        for i in range(len(p) - 1):
            mine = [d for d in own if _in_body(p[i], d)]
            if not mine:
                return p[i:]
            if any(_in_body(p[i + 1], d) for d in own):
                continue
            t = max(_leaves_at(p[i], p[i + 1], d) for d in mine)
            q = np.add(p[i], t * np.subtract(p[i + 1], p[i]))
            return [(float(q[0]), float(q[1]))] + list(p[i + 1 :])
        return [p[-1], p[-1]]

    out = trim(trim(list(pts))[::-1])[::-1]
    return out if len(out) > 1 else [out[0], out[0]]


def _na_clear(pts, own, elements) -> bool:
    """A strand drawn clear of every element and stem (by _clear), from where it leaves its own stems."""
    return _clear(_outside(pts, own), own, *elements, own_body=True)


def route_nucleic(layout: Layout, router: _Router | None = None) -> list[NAPath]:
    """Every single-stranded stretch as drawn: runs routed around the elements like loops (as are the joints
    between pieces of a bent duplex), hairpin loops as arcs, free ends straight on, strands without pairs along
    their projection, cropped gaps dashed. Every drawing is checked by _clear; one that fails is a fallback.
    Also fills layout.na_ends (where each strand's 5' and 3' labels go) and layout.na_letters (where each
    single-stranded nucleotide's letter goes), and widens layout.bounds to hold it all."""
    na = layout.nucleic
    if router is None:
        keep_out = [b for _, b in label_boxes(layout, [])]
        router = _Router(_Grid(layout, keep_out), _element_obstacles(layout), keep_out)
    elements = router.elements
    at = {key: d for d in layout.stems for key in d.axial}
    kind_of = {(r.strand, r.lo): r.kind for r in na.structure.runs}
    ends: list[tuple[int, str, tuple[float, float]]] = []
    exits: dict[tuple[int, int], tuple] = {}  # strand end -> (port, outward direction, owner) a gap leaves from

    def checked(path: NAPath, own) -> NAPath:
        path.fallback = path.fallback or not _na_clear(path.points, own, elements)
        return path

    gapped = {(g.before, "3′") for g in na.gaps} | {(g.after, "5′") for g in na.gaps}
    labels: list[tuple[float, float]] = []  # the end labels placed so far (a cropped end into a gap has none)

    def room(s: int, e: str):  # the labels an end's own label must keep clear of; None: it has no label
        return None if (s, e) in gapped else labels

    def end(s: int, e: str, xy) -> None:
        ends.append((s, e, xy))
        if (s, e) not in gapped:
            labels.append(xy)

    fixed: dict[tuple[int, int], NAPath] = {}  # (strand, nucleotide before it) -> a path of fixed shape
    between: list[tuple] = []  # strands from stem to stem still to route, once the fixed shapes are known
    for s, strand in enumerate(na.strands):
        n = len(strand)
        if s in layout.na_traces:
            pts = _xy(layout.na_traces[s])
            fixed[(s, -1)] = checked(NAPath("trace", s, list(range(n)), pts, pts), ())
            for label, a, b in (("5′", pts[0], pts[1]), ("3′", pts[-1], pts[-2])):
                v = np.subtract(a, b)
                v = v / max(float(np.hypot(*v)), 1e-9)
                end(s, label, tuple(float(c) for c in np.add(a, v * (NT_ROW + _GLYPH_W))))
                exits[(s, 0 if label == "5′" else n - 1)] = (np.add(a, v * NT_STEP / 2), v, None)
            continue
        held = [i for i in range(n) if (s, i) in at]
        if not held:
            continue
        first, last = held[0], held[-1]
        d0 = at[(s, first)]
        if first > 0:
            path, tip, v = _tail(router, layout, d0, s, (s, first), list(range(first)), True, room(s, "5′"))
            fixed[(s, -1)] = path
            exits[(s, 0)] = (np.add(path.nt_xy[0], v * NT_STEP / 2), v, d0)
        else:
            tip = _face_label(router, layout, d0, (s, 0), False, room(s, "5′"))
            exits[(s, 0)] = (*_face(d0, (s, 0), False), d0)
        end(s, "5′", tip)
        for i, j in zip(held, held[1:]):
            da, db = at[(s, i)], at[(s, j)]
            if da is db and da.side[(s, i)] == da.side[(s, j)]:
                continue  # along one rail of a stem (a bulge is drawn on the rail)
            if da is db:  # a hairpin loop: an arc over the tip
                fixed[(s, i)] = _strand_between(router, s, i, j, da, db, None)
            else:
                between.append((s, i, j, da, db, kind_of.get((s, i + 1))))
        dn = at[(s, last)]
        if last < n - 1:
            path, tip, v = _tail(router, layout, dn, s, (s, last), list(range(last + 1, n)), False, room(s, "3′"))
            fixed[(s, last)] = path
            exits[(s, n - 1)] = (np.add(path.nt_xy[-1], v * NT_STEP / 2), v, dn)
        else:
            tip = _face_label(router, layout, dn, (s, last), True, room(s, "3′"))
            exits[(s, n - 1)] = (*_face(dn, (s, last), True), dn)
        end(s, "3′", tip)
    # routes keep off the strands' end labels (as a protein loop keeps off its chain's N and C)
    g = _GLYPH_H / 2 + 0.05
    router.grid.keep_off([(x - g, y - g, x + g, y + g) for s, e, (x, y) in ends if (s, e) not in gapped])
    # and pay to cross a hairpin loop or a free end, as a loop pays to share an earlier loop's cells
    for p in fixed.values():
        if p.kind in ("hairpin", "tail"):
            _use(router.grid, p.points)
    routed = {(s, i): _strand_between(router, s, i, j, da, db, kind) for s, i, j, da, db, kind in between}
    out = [p for _, p in sorted({**fixed, **routed}.items())]
    for g in na.gaps:  # cropped away: the 3' end of one piece to the 5' end of the next, dashed
        a, b = exits.get((g.before, len(na.strands[g.before]) - 1)), exits.get((g.after, 0))
        if a is not None and b is not None:
            out.append(_gap(router, g, a, b))
    layout.na_ends = [e for e in ends if (e[0], e[1]) not in gapped]
    layout.na_letters = _place_letters(layout, out)
    _widen(layout, out)
    if log.isEnabledFor(logging.DEBUG):
        log.debug("%d nucleic-acid paths: %d fallbacks", len(out), sum(p.fallback for p in out))
    return out


def _strand_between(router: _Router, s: int, i: int, j: int, da: PlacedStem, db: PlacedStem, kind) -> NAPath:
    """The single strand from stem nucleotide i to stem nucleotide j (the m nucleotides between drawn on it)."""
    m, nts, between = j - i - 1, list(range(i + 1, j)), ((s, i), (s, j))
    own, need = (da, db), (j - i) * NT_STEP
    ra, rb = da.point((s, i)), db.point((s, j))
    ok = lambda pts: _na_clear(pts, own, router.elements)  # noqa: E731
    if da is db:  # over the stem's tip: a hairpin loop
        pts = _arc(da, (s, i), (s, j), m)
        return NAPath("hairpin", s, nts, pts, _spaced(pts, m), fallback=not ok(pts), between=between)
    (fa, ea), (fb, eb) = _face(da, (s, i), True), _face(db, (s, j), False)
    joint = m == 0 and da.parent and da.parent == db.parent
    if joint and math.dist(fa, fb) < _JOINT_STRAIGHT and _dot(fb - fa, ea) > 0:
        pts = _xy([ra, fa, fb, rb])  # the duplex just carries on
        if ok(pts):
            for c in router.grid.cells_along(fa, fb):
                router.grid.used[c] += 1
            return NAPath("joint", s, [], pts, [], between=between)
    if kind == "internal" and m <= _BOW_MAX:
        # one side of an internal loop: bowed out on its own side, so the two sides never meet
        pts = _bow(ra, rb, need, out=da.normal * da.side[(s, i)])
        if ok(pts):
            return NAPath("run", s, nts, pts, _spaced(pts, m), between=between)
    route, _ = router.route(fa, ea, da, fb, eb, db, own_body=True)
    fallback = route is None
    if route is None:
        route = _curve(fa, ea, fb, eb)
    pts = _xy([ra, *route, rb])
    if not fallback:
        pts = next((q for q in (_rounded(pts),) if ok(q)), pts)
    if not fallback and _length(pts) < need:  # too short for its nucleotides: bowed out, if that is clear
        left = np.array([ra[1] - rb[1], rb[0] - ra[0]])
        pts = next((b for b in (_bow(ra, rb, need, out=left), _bow(ra, rb, need, out=-left)) if ok(b)), pts)
    return NAPath("joint" if joint else "run", s, nts, pts, _spaced(pts, m), fallback=fallback, between=between)


def _gap(router: _Router, g, a, b) -> NAPath:
    """A cropped gap from one strand piece's 3' end to the next one's 5' end: routed, dashed and (unless short)
    labelled with the number of nucleotides left out."""
    (pa, ea, own_a), (pb, eb, own_b) = a, b
    own = tuple(d for d in (own_a, own_b) if d is not None)
    route = None
    if own_a is not None and own_b is not None:
        route, _ = router.route(pa, ea, own_a, pb, eb, own_b, own_body=True)
        if route is None:  # an end in another element's clearance margin
            route, _ = router.route(pa, ea, own_a, pb, eb, own_b, own_body=True, loose=True)
    fallback = route is None
    pts = _xy(route if route is not None else [pa, pb])
    if not fallback:
        pts = next((q for q in (_rounded(pts),) if _na_clear(q, own, router.elements)), pts)
    label = f"… {g.count} nt …" if g.count >= SHORT_GAP else ""
    fallback = fallback or not _na_clear(pts, own, router.elements)
    return NAPath("gap", g.before, [], pts, [], label, fallback)


_LETTER_GAP = 0.38  # page units: closest two base letters (or a letter and a nucleotide) may sit


def _place_letters(layout: Layout, paths: list[NAPath]) -> dict[tuple[int, int], tuple[float, float]]:
    """A letter for each single-stranded nucleotide beside its strand: on the side away from the stems it leaves
    and enters (outside a hairpin loop, in the row of the stem's letters along a free end), else the other, else
    further out on either side, or a little along the strand, whichever first keeps clear of other letters,
    nucleotides, end labels, elements and stem names; when none is clear, the one furthest from the others."""
    where = {key: d.point(key) for d in layout.stems for key in d.axial}
    where.update({(p.strand, i): np.asarray(xy) for p in paths for i, xy in zip(p.nts, p.nt_xy)})
    taken = _Points(1.0)  # every letter, nucleotide and end label so far: nearer than 1 is measured exactly
    for d in layout.stems:
        for key in d.axial:
            taken.add(d.point(key, NT_ROW))
    index = {key: taken.add(xy) for key, xy in where.items()}  # a letter's own nucleotide does not count
    for _, _, xy in layout.na_ends:
        taken.add(xy)
    bodies = [*layout.placed.values(), *layout.stems]  # by their true (tilted) outlines, not their boxes
    g = 1.2 * _GLYPH_H  # a stem's name stands a letter and a half clear of any base letter
    boxes = [np.add(_stem_label_box(d), (-g, -g, g, g)) for d in layout.stems]
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
            spaces = [room(q, index[(p.strand, i)]) for q in options]
            clear = [q for q, r in zip(options, spaces) if r >= _LETTER_GAP]
            spot = clear[0] if clear else options[int(np.argmax(spaces))]
            out[(p.strand, i)] = (float(spot[0]), float(spot[1]))
            taken.add(spot)
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
    """Grow layout.bounds to hold the strands, every base letter, the end labels and the gaps' labels."""
    pad = NT_ROW + max(_GLYPH_W, _GLYPH_H)
    boxes = [(x - pad, y - pad, x + pad, y + pad) for p in paths for x, y in p.points]
    g = max(_GLYPH_W, _GLYPH_H) / 2 + 0.05  # a letter turned any way
    letters = [d.point(key, NT_ROW) for d in layout.stems for key in d.axial] + list(layout.na_letters.values())
    letters += [xy for _, _, xy in layout.na_ends]
    boxes += [(x - g, y - g, x + g, y + g) for x, y in letters]
    for p in paths:
        if p.label:
            x, y = midway(p.points)
            half = _GLYPH_W * len(p.label) / 2 + 0.05
            boxes.append((x - half, y, x + half, y + _GLYPH_H))
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
