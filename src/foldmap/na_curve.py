# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A bent duplex drawn as one curved ladder (Style.na_bend curved, plan 25), instead of its straight pieces joined
by routed joints. The curve runs through the pieces in sequence order: from the first piece's 5' end, through
every piece's centre, to the last piece's far end (a centripetal Catmull-Rom spline, on the layout's grid). Each
piece keeps its base pairs at its own spacing (PlacedStem.axial) around its centre on the curve, and no step is
drawn shorter than the duplex's own (pieces that would overlap along the curve slide apart, as little as that
takes: _isotonic). Where the band (rails and base letters) would come over itself and the duplex winds a full turn
or more round the middle of its pieces (the nucleosome's two gyres seen down the superhelix axis), it is drawn as
gyres one inside the other: a ring through the pieces' angles round that middle, the first gyre outside (or
inside, as most of the pieces show it), parted by the band's width (_gyres). The pieces stay in layout.stems;
render draws the curve in their place (layout.na_curves) and route_nucleic's strands that leave or enter the duplex
follow its nucleotides (follow). In projected mode build_layout lays the curve out first (lay_out): the pieces end
to end on it at the duplex's own rung spacing (the nucleosome's gyres grown round the octamer, their rungs spread
evenly), pinned there, the protein laid out round its band."""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from .layout import (
    _CONTACT_CAP,
    _GLYPH_H,
    _GLYPH_W,
    _W_CONTACT,
    DNA_RISE,
    MARGIN,
    NT_ROW,
    SCALE,
    Layout,
    Placed,
    PlacedStem,
    Walls,
    _snap,
    centred_box,
    glyph,
)

CURVE_STEP = 0.1  # page units between the curve's samples
_FAR = 2.0  # two places on the band this many band widths apart along it are different stretches (_far_gap)
_GYRE_GAP = 0.7  # page units kept between two gyres' bands (their letters) beyond the band's width: enough that the
# boxes round two of their stretches (figcheck A5: as any two stems, chords half a rail long) meet nowhere either
_GYRE_GROW, _GYRE_TRIES = 1.15, 8  # gyres parted by too little: their spacing grows by this factor, this often
_GYRE_SPREAD = math.radians(40)  # the ring's radius at an angle: its pieces' radii within about this angle
_RAMP_MIN = 0.125  # turns: the shortest stretch over which the band goes from one gyre to the other
_R_MIN = 3.0  # rail offsets: the least radius the curve bends at (the inner rail's steps keep 2/3 of their length)
_EASE_STEP, _EASE_ROUNDS = 0.5, 800  # page units between the samples eased (_eased), and at most this many rounds
_NAME_AT = (0.5, 0.42, 0.58, 0.34, 0.66, 0.26, 0.74, 0.18, 0.82)  # fractions of the duplex its name is tried at


@dataclass(eq=False)
class NACurve:
    """One bent duplex drawn as a curve. xy: its centreline, s: the arc length at each point, tangent: the unit
    direction there (all on the grid), running the way the first arm runs 5'->3'. at: each nucleotide's arc length;
    side: +1 on the rail left of the direction, -1 on the right. pieces: the layout.stems it stands for, in order
    along the first arm. gyres: drawn as gyres one inside the other; overlap: the band still comes over itself."""

    parent: str
    name: str
    pieces: list[int]
    xy: np.ndarray
    s: np.ndarray
    tangent: np.ndarray
    rail: float
    width: float
    at: dict
    side: dict
    pairs: list
    gyres: bool = False
    overlap: bool = False
    name_xy: tuple[float, float] = (0.0, 0.0)
    name_angle: float = 0.0
    band: float = 0.0  # its rails and base letters across (2 x (rail + NT_ROW + a glyph's height))
    held: list = field(default_factory=list)  # laid out by build_layout (lay_out): each piece's (cx, cy, angle) there
    named: bool = False  # its name placed (build_curves: name_xy, name_angle); text keeps off it from then on
    name_keep: tuple | None = None  # its name's ink grown by layout._INK_MARGIN em: what routes keep out of
    drawn: dict = field(default_factory=dict, repr=False, compare=False)  # its Chords, kept (chords: read often)

    def frames(self, s) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(points, unit tangents, unit left normals) at arc lengths s; past either end the curve runs straight on."""
        s = np.atleast_1d(np.asarray(s, float))
        p = np.column_stack([np.interp(s, self.s, self.xy[:, 0]), np.interp(s, self.s, self.xy[:, 1])])
        t = np.column_stack([np.interp(s, self.s, self.tangent[:, 0]), np.interp(s, self.s, self.tangent[:, 1])])
        lo, hi = s < self.s[0], s > self.s[-1]
        p[lo] = self.xy[0] + (s[lo] - self.s[0])[:, None] * self.tangent[0]
        p[hi] = self.xy[-1] + (s[hi] - self.s[-1])[:, None] * self.tangent[-1]
        t = t / np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-12)[:, None]
        return p, t, np.column_stack([-t[:, 1], t[:, 0]])

    def to_page(self, x, y) -> np.ndarray:
        """Curve-local (arc length, offset to the left) to page coordinates."""
        p, _, n = self.frames(x)
        return p + np.atleast_1d(np.asarray(y, float))[:, None] * n

    def point(self, key: tuple[int, int], off: float = 0.0) -> np.ndarray:
        """Page position of a nucleotide on its rail, or `off` further out from the axis (PlacedStem.point)."""
        return _snap(self.to_page(self.at[key], self.side[key] * (self.rail + off))[0])

    def angle(self, key: tuple[int, int]) -> float:
        """The curve's direction (radians) at a nucleotide."""
        t = self.frames(self.at[key])[1][0]
        return float(_snap(math.atan2(t[1], t[0])))

    def span(self) -> tuple[float, float]:
        """Arc lengths of the duplex's first and last nucleotides (the outermost rungs)."""
        v = list(self.at.values())
        return min(v), max(v)

    def extent(self) -> tuple[float, float]:
        """Arc lengths of the duplex's two ends: half a base-pair step beyond its outermost rungs (as a piece's length
        reaches beyond them)."""
        lo, hi = self.span()
        return lo - DNA_RISE * SCALE / 2, hi + DNA_RISE * SCALE / 2

    def offset_line(self, s0: float, s1: float, off: float) -> np.ndarray:
        """The line `off` to the left of the curve (right if negative) from arc length s0 to s1."""
        inner = self.s[(self.s > s0) & (self.s < s1)]
        return _snap(self.to_page(np.concatenate([[s0], inner, [s1]]), off))


def curvable(na, parent: str) -> bool:
    """Whether the bent duplex `parent` (a whole stem's id) may be drawn as one curve: none of its chains pairs with
    itself in a stem outside the duplex's own helix (Stem.helix: stems continued across internal loops). A duplex
    inside a folded RNA stays in pieces, which read with the stems round it; a lone duplex (DNA bent by a protein,
    one RNA helix) curves. IHF's DNA 1 curves although its chain C also pairs in DNA 2: that stem pairs C with D."""
    if na is None:
        return True
    own = next((st for st in na.structure.stems if st.id == parent), None)
    if own is None:
        return True
    chains = {na.strands[own.a].chain, na.strands[own.b].chain}
    for st in na.structure.stems:
        if st.id == parent or st.parent == parent or (own.helix >= 0 and st.helix == own.helix):
            continue
        if na.strands[st.a].chain == na.strands[st.b].chain and na.strands[st.a].chain in chains:
            return False
    return True


def _pieces_of(stems: list[PlacedStem]) -> dict[str, list[int]]:
    """Parent id -> the indices in `stems` of its pieces, for each duplex cut into two or more pieces that follow
    one another (piece numbers #1, #2, ... with none missing: a crop that left out a middle piece leaves pieces)."""
    groups: dict[str, list[int]] = {}
    for k, d in enumerate(stems):
        if d.parent:
            groups.setdefault(d.parent, []).append(k)

    def number(k: int) -> int | None:
        tail = stems[k].id.rpartition("#")[2]
        return int(tail) if tail.isdigit() else None

    out = {}
    for parent, ks in groups.items():
        nums = sorted(number(k) for k in ks) if all(number(k) is not None for k in ks) else None
        if len(ks) > 1 and (nums is None or nums == list(range(nums[0], nums[0] + len(nums)))):
            out[parent] = ks
    return out


def bent_duplexes(stems: list[PlacedStem], na) -> list[tuple[str, list[int]]]:
    """(parent id, its pieces' indices in `stems` in order along the first arm 5'->3') of each duplex cut into two
    or more pieces that is drawn as one curve (curvable)."""

    def start(k: int) -> tuple[int, int]:  # the piece's first-arm nucleotide nearest the 5' end
        return min(key for key in stems[k].axial if stems[k].side[key] == stems[k].first)

    out = [(parent, sorted(ks, key=start)) for parent, ks in _pieces_of(stems).items() if curvable(na, parent)]
    return sorted(out, key=lambda g: g[1][0])


def left_in_pieces(stems: list[PlacedStem], na) -> list[str]:
    """The names of the bent duplexes curvable keeps in pieces (Style.na_bend curved: the figure's note)."""
    names = {st.id: st.name for st in na.structure.stems} if na is not None else {}
    return [names.get(p) or p for p in _pieces_of(stems) if not curvable(na, p)]


def _arc(xy: np.ndarray) -> np.ndarray:
    steps = np.hypot(np.diff(xy[:, 0]), np.diff(xy[:, 1]))
    return _snap(np.concatenate([[0.0], np.cumsum(steps)]))


def _tangents(xy: np.ndarray) -> np.ndarray:
    t = np.gradient(xy, axis=0)
    return _snap(t / np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-12)[:, None])


def _catmull_rom(ctrl: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """A centripetal Catmull-Rom spline through the control points, sampled about every CURVE_STEP: the points and
    the index of each control point among them."""
    ext = np.vstack([2 * ctrl[0] - ctrl[1], ctrl, 2 * ctrl[-1] - ctrl[-2]])
    out, idx, count = [], [], 0
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1 : i + 3]
        chord = math.hypot(*(p2 - p1))
        t1 = max(math.hypot(*(p1 - p0)), 1e-6) ** 0.5
        t2 = t1 + max(chord, 1e-6) ** 0.5
        t3 = t2 + max(math.hypot(*(p3 - p2)), 1e-6) ** 0.5
        n = max(4, math.ceil(chord / CURVE_STEP))
        t = (t1 + (t2 - t1) * np.arange(n) / n)[:, None]
        a1 = (t1 - t) / t1 * p0 + t / t1 * p1
        a2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2
        a3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3
        b1 = (t2 - t) / t2 * a1 + t / t2 * a2
        b2 = (t3 - t) / (t3 - t1) * a2 + (t - t1) / (t3 - t1) * a3
        idx.append(count)
        out.append((t2 - t) / (t2 - t1) * b1 + (t - t1) / (t2 - t1) * b2)
        count += n
    idx.append(count)
    return np.vstack([*out, ctrl[-1:]]), idx


def _distinct(xy: np.ndarray, idx: list[int]) -> tuple[np.ndarray, list[int]]:
    """The points on the grid with repeats dropped (a zero step has no direction), and where each index went."""
    xy = _snap(xy)
    keep = np.concatenate([[True], np.any(np.diff(xy, axis=0) != 0, axis=1)])
    where = np.cumsum(keep) - 1
    return xy[keep], [int(where[i]) for i in idx]


def _far_gap(xy: np.ndarray, s: np.ndarray, band: float) -> float:
    """The least distance between two places on the curve at least _FAR band widths apart along it (inf if none):
    below `band`, the band comes over itself there."""
    if s[-1] < _FAR * band:
        return math.inf
    at = np.searchsorted(s, np.arange(0.0, float(s[-1]), 0.25))
    p, q = xy[at], s[at]
    d = np.hypot(p[:, None, 0] - p[None, :, 0], p[:, None, 1] - p[None, :, 1])
    far = np.abs(q[:, None] - q[None, :]) >= _FAR * band
    return float(d[far].min()) if far.any() else math.inf


def _isotonic(y, w=None) -> np.ndarray:
    """The nondecreasing sequence nearest y (weights w) in least squares: pool adjacent violators."""
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    blocks: list[list[float]] = []  # [mean, weight, count]
    for v, wt in zip(np.asarray(y, float), w):
        blocks.append([float(v), float(wt), 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            b = blocks.pop()
            a = blocks[-1]
            tw = a[1] + b[1]
            blocks[-1] = [(a[0] * a[1] + b[0] * b[1]) / tw, tw, a[2] + b[2]]
    return np.array([m for m, _, c in blocks for _ in range(c)])


def _smoothstep(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def _gyres(ctrl: np.ndarray, band: float) -> tuple[np.ndarray, list[float]] | None:
    """The curve as gyres one inside the other round the middle of the pieces (the box round their centres), or
    None where the pieces wind less than a full turn round it. Each control point keeps its angle there (made to
    turn one way only); the ring's radius at an angle is its pieces' radii near that angle (_GYRE_SPREAD); the first
    gyre is drawn outside the second where the pieces of the first stand further out on average (else inside),
    the two parted by the band's width (and _GYRE_GAP, grown until the bands clear) where both go round, the band
    crossing from one to the other where only one does. More than two turns: a spiral parted by that much a turn.
    Returns the curve's points and each control point's arc length."""
    inner = ctrl[1:-1]
    c = _snap([(inner[:, 0].min() + inner[:, 0].max()) / 2, (inner[:, 1].min() + inner[:, 1].max()) / 2])
    rel = ctrl - c
    phi = np.arctan2(rel[:, 1], rel[:, 0])
    d = (np.diff(phi) + np.pi) % (2 * np.pi) - np.pi
    a = np.concatenate([[0.0], np.cumsum(d)])
    sgn = 1.0 if a[-1] >= 0 else -1.0
    a = _snap(_isotonic(sgn * a))
    turns = float(a[-1] - a[0]) / (2 * np.pi)
    if turns < 1.0:
        return None
    r = np.hypot(rel[:, 0], rel[:, 1])
    t_k = (a - a[0]) / (2 * np.pi)
    first, second = r[t_k <= turns - 1], r[t_k >= 1]
    out = 1.0 if not len(first) or not len(second) or first.mean() >= second.mean() else -1.0
    angle_k = phi[0] + sgn * a

    def base(angles: np.ndarray) -> np.ndarray:  # the ring's radius at each angle
        dd = (angles[:, None] - angle_k[None, :] + np.pi) % (2 * np.pi) - np.pi
        w = np.exp(-((dd / _GYRE_SPREAD) ** 2))
        return (w * r[None, :]).sum(axis=1) / w.sum(axis=1)

    def offset(t: np.ndarray, gap: float) -> np.ndarray:  # out from the ring, the first gyre's way
        if turns > 2:
            return out * gap * (turns / 2 - t)
        half = max((2 - turns) / 2, _RAMP_MIN)
        return out * gap * (0.5 - _smoothstep((t - (turns / 2 - half)) / (2 * half)))

    n = math.ceil(float(a[-1] - a[0]) * float(r.max() + band) / CURVE_STEP) + 1
    aa = np.linspace(float(a[0]), float(a[-1]), n)
    ring = base(phi[0] + sgn * aa)
    gap = band + _GYRE_GAP
    ang = phi[0] + sgn * aa
    for _ in range(_GYRE_TRIES):
        rr = np.maximum(ring + offset((aa - a[0]) / (2 * np.pi), gap), band)
        raw = _snap(np.column_stack([c[0] + rr * np.cos(ang), c[1] + rr * np.sin(ang)]))
        keep = np.concatenate([[True], np.any(np.diff(raw, axis=0) != 0, axis=1)])
        xy = raw[keep]
        s = _arc(xy)
        if _far_gap(xy, s, band) >= band:
            break
        gap *= _GYRE_GROW
    return xy, [float(v) for v in _snap(np.interp(a, aa[keep], s))]


def _resample(xy: np.ndarray, step: float) -> tuple[np.ndarray, float]:
    """The curve sampled every `step` (about) along its length, and that step."""
    s = _arc(xy)
    n = max(2, math.ceil(float(s[-1]) / step) + 1)
    t = np.linspace(0.0, float(s[-1]), n)
    return np.column_stack([np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1])]), float(s[-1]) / (n - 1)


def _eased(xy: np.ndarray, at: list[float], r_min: float) -> tuple[np.ndarray, list[float]]:
    """The curve with every bend tighter than radius r_min eased out, so the inner rail and its letters keep room:
    sampled every _EASE_STEP, each sample on a tighter bend (and its two neighbours either way) moved half way to
    the middle of its neighbours, round after round (the ends fixed), then drawn through them again. `at`: arc
    lengths on the curve given, returned on the curve eased (by sample)."""
    p, h = _resample(xy, _EASE_STEP)
    if len(p) < 5:
        return xy, at
    for _ in range(_EASE_ROUNDS):
        a, b = p[1:-1] - p[:-2], p[2:] - p[1:-1]
        turn = np.abs(np.arctan2(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0], a[:, 0] * b[:, 0] + a[:, 1] * b[:, 1]))
        reach = (np.hypot(a[:, 0], a[:, 1]) + np.hypot(b[:, 0], b[:, 1])) / 2
        tight = np.concatenate([[False], turn * r_min * 1.1 > reach, [False]])  # 10% to spare (drawn through again)
        if not tight.any():
            break
        for _ in range(2):  # and two samples either way
            tight[1:] |= tight[:-1].copy()
            tight[:-1] |= tight[1:].copy()
        tight[0] = tight[-1] = False
        mid = (p[:-2] + p[2:]) / 2
        move = tight[1:-1]
        p[1:-1][move] += 0.5 * (mid - p[1:-1])[move]
        p = _snap(p)
    dense, idx = _distinct(*_catmull_rom(_snap(p)))
    s = _arc(dense)
    return dense, [float(v) for v in _snap(np.interp(np.asarray(at) / h, np.arange(len(idx)), s[idx]))]


def _curve_of(stems: list[PlacedStem], ks: list[int], band: float) -> tuple[np.ndarray, list[float], bool, bool]:
    """The centreline of one bent duplex (points on the grid), the arc length of each piece's centre on it, whether
    it is drawn as gyres, and whether its band still comes over itself."""
    ps = [stems[k] for k in ks]
    ctrl = _ctrl(ps)
    r_min = _R_MIN * ps[0].rail
    xy, idx = _distinct(*_catmull_rom(ctrl))
    xy, centres = _eased(xy, [float(v) for v in _arc(xy)[idx[1:-1]]], r_min)
    if _far_gap(xy, _arc(xy), band) >= band:
        return xy, centres, False, False
    wound = _gyres(ctrl, band)
    if wound is None:
        return xy, centres, False, True
    gxy, at = _eased(*wound, r_min)
    return gxy, at[1:-1], True, _far_gap(gxy, _arc(gxy), band) < band


def _name_place(c: NACurve, layout: Layout, curves: list[NACurve], gh: float) -> None:
    """Where the duplex's one name goes (c.name_xy, c.name_angle, and c.name_keep: its ink grown by _INK_MARGIN em,
    what routes keep out of, as a stem's label_keep). Tried at each of _NAME_AT along the duplex, on the side away
    from the middle of the curve first: beside its letters there as a straight stem's name stands beside its row
    (layout._name_at): its ink, turned as drawn, _INK_MARGIN em beyond the ink of the letters on that side within
    reach along the curve (figcheck.stem_name_gaps measures it so, A17); with no letters drawn, a glyph and a half
    beyond their place. The first spot where no other stretch of the band comes nearer than that and its box meets
    no element, element name, other stem's or curve's body, other stem's or curve's name; else the spot farthest
    from all of those."""
    from .layout import (
        _INK_MARGIN,
        _LETTER_EM,
        _NAME_EM,
        _box_meets_body,
        _grown,
        _ink_em,
        _ink_points,
        _reach,
        stem_keep_box,
    )

    row = getattr(layout, "na_letter_rows", (None, False))[0]
    na = layout.nucleic
    em = _ink_em(layout, layout.font_scale)
    margin = _INK_MARGIN * em
    lo, hi = c.span()
    middle = c.xy.mean(axis=0)
    hidden = {k for cc in curves for k in cc.pieces}  # every curve being drawn: its pieces have no body or name
    others = [d for k, d in enumerate(layout.stems) if k not in hidden]
    for m, cc in enumerate(curves):
        if cc is not c:
            others += chords(cc, m, layout.stems)
    rects = [p.rect for p in [*layout.placed.values(), *layout.ghosts]]
    rects += [b for name, b in _element_names(layout)]
    names = [stem_keep_box(d, glyph(layout)[0] / _GLYPH_W) for k, d in enumerate(layout.stems) if k not in hidden]
    names += [cc.name_keep for cc in curves if cc is not c and cc.named and cc.name_keep is not None]
    letters = {}
    if row is not None and na is not None:
        for key in c.at:
            ink = _ink_points(na.strands[key[0]].bases[key[1]], False, _LETTER_EM * em, c.angle(key))
            letters[key] = ink + c.point(key, row)
    best = None

    def beside(s: float, side: float):  # the name centred beside the letters at arc length s, on `side`
        p, t, nrm = (v[0] for v in c.frames(s))
        angle = float(_snap(math.atan2(t[1], t[0])))
        ink = _ink_points(c.name, True, _NAME_EM * em, angle)
        reach = float(np.hypot(*(ink.max(axis=0) - ink.min(axis=0)))) / 2 + c.band if len(ink) else c.band
        way = side * nrm
        inner = float(np.max(ink @ -way, initial=0.0))
        row_ink = [q for key, q in letters.items() if c.side[key] == side and abs(c.at[key] - s) <= reach]
        if row_ink:
            out = max(float(np.max((q - p) @ way)) for q in row_ink)
            off = math.ceil((out + margin + inner) * 1e4) / 1e4
        else:
            off = c.rail + NT_ROW + 1.5 * gh
        return _snap(p + off * way), angle, ink, off

    for f in _NAME_AT:
        s = lo + f * (hi - lo)
        p, _, nrm = (v[0] for v in c.frames(s))
        away = 1.0 if math.hypot(*(p + nrm - middle)) >= math.hypot(*(p - nrm - middle)) else -1.0
        for side in (away, -away):
            at = s
            for _ in range(4):  # measured, as figcheck does, at the curve's point nearest the name's ink
                q, angle, ink, off = beside(at, side)
                mid = q + (ink.min(axis=0) + ink.max(axis=0)) / 2 if len(ink) else q
                near = float(c.s[int(np.argmin(np.hypot(c.xy[:, 0] - mid[0], c.xy[:, 1] - mid[1])))])
                if abs(near - at) < 1e-3:
                    break
                at = near
            box = tuple(float(v) for v in _reach(ink) + np.array([q[0], q[1], q[0], q[1]]))
            room = float(np.hypot(c.xy[:, 0] - q[0], c.xy[:, 1] - q[1]).min()) - (off - 0.01)
            gaps = [_gap(box, r) for r in rects + names] + [_gap(box, d.rect) for d in others]
            hit = any(g <= 0 for g in gaps[: len(rects) + len(names)]) or any(
                _box_meets_body(box, d) for d in others if _gap(box, d.rect) <= 0
            )
            score = min([room, *gaps]) if gaps else room
            if room >= 0 and not hit:
                best = (math.inf, q, angle, ink)
                break
            if best is None or score > best[0] + 1e-9:
                best = (score, q, angle, ink)
        if best[0] == math.inf:
            break
    _, q, angle, ink = best
    c.name_xy, c.name_angle = tuple(float(v) for v in q), angle
    c.name_keep = _grown(tuple(float(v) for v in _reach(ink) + np.array([q[0], q[1], q[0], q[1]])), margin)


def _gap(a, b) -> float:
    """The distance between two boxes (x0, y0, x1, y1), 0 where they meet or overlap."""
    return float(math.hypot(max(a[0] - b[2], b[0] - a[2], 0.0), max(a[1] - b[3], b[1] - a[3], 0.0)))


def _element_names(layout: Layout) -> list:
    """The elements' name boxes (layout.label_boxes' label: entries)."""
    from .layout import label_boxes

    return [(name, b) for name, b in label_boxes(layout, []) if name.startswith("label:")]


def _ctrl(ps: list[PlacedStem]) -> np.ndarray:
    """The curve's control points: the first piece's 5' end, each piece's centre, the last piece's far end."""
    ctrl = [np.array([ps[0].cx, ps[0].cy]) - ps[0].direction * ps[0].length / 2]
    ctrl += [np.array([d.cx, d.cy]) for d in ps]
    ctrl += [np.array([ps[-1].cx, ps[-1].cy]) + ps[-1].direction * ps[-1].length / 2]
    return _snap(np.array(ctrl))


def _make(parent: str, ks: list[int], ps: list[PlacedStem], xy, centres, gyres, overlap, na, along=None) -> NACurve:
    """The NACurve of one bent duplex on the centreline xy (`centres`: each piece's centre's arc length there), its
    nucleotides placed, its first arm's side voted by the pieces. `along` None (the pieces where the layout left
    them): each piece's pairs at its own spacing round its centre's place, pieces that would overlap slid apart
    (_isotonic). Else (off, stretch), the pieces laid out on the curve (lay_out): end to end, the duplex's 5' end at
    arc length off, every step `stretch` times its own."""
    rail, width = ps[0].rail, ps[0].width
    s = _arc(xy)
    tangent = _tangents(xy)
    tmp = NACurve(parent, "", ks, xy, s, tangent, rail, width, {}, {}, [])
    lengths = [d.length for d in ps]
    U = np.cumsum([0.0, *lengths[:-1]]) + np.array(lengths) / 2  # each centre along the straightened duplex
    if along is None:
        place, stretch = U + _isotonic(np.array(centres) - U, lengths), 1.0  # each piece's place along the curve
    else:
        off, stretch = along
        place = off + stretch * U
    vote = 0.0  # the first arm's side of the curve, as the pieces show it (those seen most clearly most)
    for d, sk in zip(ps, centres):
        t = tmp.frames(sk)[1][0]
        along_t = float(t[0]) * float(d.direction[0]) + float(t[1]) * float(d.direction[1])
        vote += (d.margin + 0.01) * d.first * (1.0 if along_t >= 0 else -1.0)
    first_side = 1 if float(_snap(vote)) >= 0 else -1
    at, side, pairs = {}, {}, []
    for d, u in zip(ps, place):
        for key, x in d.axial.items():
            at[key] = float(_snap(u + stretch * x))
            side[key] = first_side if d.side[key] == d.first else -first_side
        pairs += list(d.pairs)
    whole = next((st.name for st in na.structure.stems if st.id == parent), "") if na is not None else ""
    name = whole or ps[0].name.rstrip("abcdefghijklmnopqrstuvwxyz").strip()
    return NACurve(parent, name, ks, xy, s, tangent, rail, width, at, side, pairs, gyres, overlap)


def build_curves(layout: Layout) -> list[NACurve]:
    """Every bent duplex of the layout as one curve (NACurve), its nucleotides placed on it: the curve build_layout
    laid it out on (lay_out) while its pieces stand where that put them, else one through the pieces as placed."""
    stems = layout.stems
    na = layout.nucleic
    gw, gh = glyph(layout)
    laid = {c.parent: c for c in layout.na_curves}
    out = []
    for parent, ks in bent_duplexes(stems, na):
        ps = [stems[k] for k in ks]
        band = 2 * (ps[0].rail + NT_ROW + gh)  # rails and base letters
        c = laid.get(parent)
        if c is None or c.pieces != ks or c.held != [(d.cx, d.cy, d.angle) for d in ps]:
            xy, centres, gyres, overlap = _curve_of(stems, ks, band)
            c = _make(parent, ks, ps, xy, centres, gyres, overlap, na)
        c.band = band
        c.overlap = _far_gap(c.xy, _arc(c.xy), band) < band  # the band over itself, on the curve drawn (A5)
        out.append(c)
    for c in out:  # then each name, clear of every curve's body and of the names placed before it
        _name_place(c, layout, out, gh)
        c.named = True
    return out


# ---- build_layout's part (Style.na_bend curved, projected mode): the pieces end to end, the protein round them ----

_ROOM = 1.0  # the inner gyre holds this many times the footprints (MARGIN round each) of the items inside it
_RING_TRIES = 8  # the ring is grown about its middle at most this many times until it does
_RUNG_W = 3.0  # a binding element's pull to its contact rungs at full contact, in layout._W_CONTACT


def _chain(ps: list[PlacedStem]) -> None:
    """The pieces end to end, each at its own angle (as projected) and length, the chain standing where the pieces'
    centres stood on average (weighted by their lengths)."""
    pos, j = [], np.zeros(2)
    for d in ps:
        c = j + d.direction * d.length / 2
        pos.append(c)
        j = c + d.direction * d.length / 2
    pos = np.array(pos)
    w = np.array([d.length for d in ps], float)
    old = np.array([[d.cx, d.cy] for d in ps], float)
    pos = _snap(pos + (w[:, None] * (old - pos)).sum(axis=0) / w.sum())
    for d, p in zip(ps, pos.tolist()):
        d.cx, d.cy = float(p[0]), float(p[1])


def _ring(ps: list[PlacedStem], band: float, need: float) -> tuple[np.ndarray, list[float]] | None:
    """The duplex as gyres one inside the other (_gyres) round the middle of its pieces, grown about that middle
    until the inside of the inner gyre (less half the band) holds an area `need`: the points and each control
    point's arc length (_ctrl's), or None when the pieces wind less than a full turn."""
    ctrl = _ctrl(ps)
    inner = ctrl[1:-1]
    mid = _snap([(inner[:, 0].min() + inner[:, 0].max()) / 2, (inner[:, 1].min() + inner[:, 1].max()) / 2])
    f, got = 1.0, None
    for _ in range(_RING_TRIES):
        wound = _gyres(_snap(mid + f * (ctrl - mid)), band)
        if wound is None:
            return got
        got = _eased(*wound, _R_MIN * ps[0].rail)
        r_in = float(np.hypot(got[0][:, 0] - mid[0], got[0][:, 1] - mid[1]).min()) - band / 2
        if math.pi * max(r_in, 0.0) ** 2 >= need:
            break
        f = float(_snap(f * max(1.05, math.sqrt(need / math.pi) / max(r_in, band))))
    return got


def _onto(c: NACurve, ps: list[PlacedStem], off: float, stretch: float) -> None:
    """Each piece set on the curve: centred on its place there, along the curve's direction, its rails on the
    curve's sides, its pairs `stretch` times as far apart."""
    lengths = [d.length for d in ps]
    U = np.cumsum([0.0, *lengths[:-1]]) + np.array(lengths) / 2
    for d, u in zip(ps, U):
        p, t, _ = (v[0] for v in c.frames(off + stretch * u))
        d.cx, d.cy = (float(v) for v in _snap(p))
        d.angle = float(_snap(math.atan2(t[1], t[0])))
        arm = next(key for key in sorted(d.axial) if d.side[key] == d.first)  # a nucleotide of its first arm
        d.first = c.side[arm]
        d.side = {key: c.side[key] for key in d.axial}
        if stretch != 1.0:
            d.axial = {key: float(_snap(stretch * x)) for key, x in d.axial.items()}
            d.length = float(_snap(stretch * d.length))
    c.held = [(d.cx, d.cy, d.angle) for d in ps]


def _band_boxes(c: NACurve) -> tuple[list[tuple[float, float, float, float]], list[tuple[float, ...]]]:
    """The curve's band (rails and letters) as squares two across, a quarter of the band apart along it, and for each
    the band it stands for (Walls.band: its centreline point, the unit normal there, half the band)."""
    lo, hi = c.extent()
    h = c.band / 4
    s = np.arange(lo, hi + h / 2, h)
    p, _, n = c.frames(s)
    out = []
    for q in np.vstack([p + h * n, p - h * n]).tolist():
        out.append((q[0] - h, q[1] - h, q[0] + h, q[1] + h))
    rows = [(*a, *b, c.band / 2) for a, b in zip(np.vstack([p, p]).tolist(), np.vstack([n, n]).tolist())]
    return [tuple(float(v) for v in b) for b in _snap(out)], [tuple(float(v) for v in r) for r in rows]


def lay_out(items: list, na, style, contacts: dict, walls: Walls | None, ties: dict) -> tuple[list, Walls | None, list]:
    """Style.na_bend curved in build_layout (projected mode), before _attract: each bent duplex laid out as one
    curve, the protein to be laid out round it. Its pieces (as projected) are set end to end (_chain) and a curve
    drawn through them (_curve_of); the duplex lies on it at its own rung spacing. Where it winds a full turn or
    more round the middle of its pieces (gyres, the nucleosome), the ring is grown until the items inside it (their
    homes nearer its middle than its pieces on average) fit inside the inner gyre (_ROOM), and its rungs are spread
    evenly over it. Each piece is then set on the curve (_onto) and its item pinned (_Item.pinned: held there,
    _attract and _separate move the rest round it). Returns the curves (pieces as items' members), the walls the
    other items keep off (`walls` with the curves' bands added, _band_boxes) and springs (as _stem_springs')
    pulling each element touching a curve's nucleotides to its band's edge at the longest run of rungs it touches,
    on the side most of its chain is projected on (so its loops need not cross the band; past an inner gyre where one
    stands that way: _beyond); one projected across the band from there starts beside it (the walls would push it off
    the band the way they find). An element no rung pulls whose neighbours in sequence all stand on its chain's side
    starts there too (_to_its_side). The contacts between the binding elements and the pieces are dropped from
    `contacts` (the springs pull instead). A duplex any of whose pieces a symmetry tie holds is left as pieces."""
    stems, item_of = [], []
    for k, it in enumerate(items):
        for m, (ox, oy) in zip(it.members, it.offsets):
            if isinstance(m, PlacedStem):
                m.cx, m.cy = it.cx + ox, it.cy + oy
                stems.append(m)
                item_of.append(k)
    tied = set(ties or {}) | {j for j, _ in (ties or {}).values()}
    groups = [(p, ks) for p, ks in bent_duplexes(stems, na) if not any(item_of[k] in tied for k in ks)]
    if not groups:
        return [], walls, []
    gh = _GLYPH_H * max(1.0, float(style.font_scale))
    curves = []
    for parent, ks in groups:
        ps = [stems[k] for k in ks]
        band = 2 * (ps[0].rail + NT_ROW + gh)
        _chain(ps)
        xy, centres, gyres, overlap = _curve_of(stems, ks, band)
        total = float(sum(d.length for d in ps))
        U = np.cumsum([0.0, *[d.length for d in ps][:-1]]) + np.array([d.length for d in ps]) / 2
        w = np.array([d.length for d in ps], float)
        along = (float(_snap((w * (np.array(centres) - U)).sum() / w.sum())), 1.0)
        if gyres:
            mid = _ctrl(ps)[1:-1].mean(axis=0)
            reach = float(np.mean([math.hypot(d.cx - mid[0], d.cy - mid[1]) for d in ps]))
            mine = {item_of[k] for k in ks}
            need = 0.0
            for k, it in enumerate(items):
                if k in mine or any(isinstance(m, PlacedStem) for m in it.members):
                    continue
                hx, hy = it.half
                mx, my = it.mid
                if math.hypot(it.home[0] + mx - mid[0], it.home[1] + my - mid[1]) < reach:
                    need += (2 * hx + MARGIN) * (2 * hy + MARGIN)
            ring = _ring(ps, band, _ROOM * need)
            if ring is not None:
                xy, at = ring
                centres = at[1:-1]
                overlap = _far_gap(xy, _arc(xy), band) < band
                along = (float(at[0]), float(_snap(max(1.0, (at[-1] - at[0]) / total))))
        c = _make(parent, ks, ps, xy, centres, gyres, overlap, na, along)
        c.band = band
        _onto(c, ps, *along)
        curves.append(c)
    _part(curves, stems)  # laid out each on its own and then pinned: nothing else parts two curves
    for c in curves:
        for k in c.pieces:
            it = items[item_of[k]]
            it.cx, it.cy = stems[k].cx, stems[k].cy  # (a stem is an item of its own: offset 0)
            it.home = (it.cx, it.cy)
            it.pinned = True
    pinned = {item_of[k] for c in curves for k in c.pieces}
    boxes, rows = [], []
    for c in curves:
        b, r = _band_boxes(c)
        boxes += b
        rows += r
    hard = np.ones((len(items), len(boxes)), bool)
    hard[sorted(pinned)] = False
    band = np.asarray(rows, float).reshape(-1, 5)
    if walls is not None:
        before = walls.band if walls.band is not None else np.full((len(walls.boxes), 5), np.nan)
        walls = Walls(
            np.vstack([walls.boxes, np.asarray(boxes, float).reshape(-1, 4)]),
            np.hstack([walls.hard, hard]),
            np.vstack([before, band]),
        )
    else:
        walls = Walls(np.asarray(boxes, float).reshape(-1, 4), hard, band)
    ids = {stems[k].id for c in curves for k in c.pieces}
    for key in [key for key in contacts if key & ids]:
        del contacts[key]
    elements = [  # (item, member, its offset, its projected place) of every protein element
        (k, m, off, np.add(it.home, off))
        for k, it in enumerate(items)
        if k not in pinned
        for m, off in zip(it.members, it.offsets)
        if isinstance(m, Placed) and not m.ghost
    ]
    sides = {}  # (curve, chain) -> the side of the curve most of the chain's elements are projected on
    for n, c in enumerate(curves):
        lo, hi = c.span()
        p, t, _ = c.frames(np.linspace(lo, hi, max(2, math.ceil((hi - lo) / CURVE_STEP) + 1)))
        for _, m, _, q in elements:
            near = int(np.argmin(np.hypot(p[:, 0] - q[0], p[:, 1] - q[1])))
            d = q - p[near]
            v = float(_snap(t[near, 0] * d[1] - t[near, 1] * d[0]))
            sides[(n, m.sse.chain)] = sides.get((n, m.sse.chain), 0) + (1 if v >= 0 else -1)
    springs, start = [], {}
    for k, m, off, _ in elements:
        for n, c in enumerate(curves):
            res, keys = [], []
            for r in range(m.sse.start, m.sse.end + 1):
                hit = [x for x in na.contacts.get(r, ()) if x in c.at]
                if hit:
                    res.append(r)
                    keys += hit
            if not res:
                continue
            arcs = sorted(c.at[x] for x in keys)  # the rungs it touches, in runs along the duplex: the longest run
            runs = [[arcs[0]]]
            for v in arcs[1:]:
                if v - runs[-1][-1] > c.band:
                    runs.append([v])
                else:
                    runs[-1].append(v)
            a = float(np.mean(max(runs, key=len)))
            p, _, nrm = (v[0] for v in c.frames(a))
            f = (float(np.mean(res)) - m.sse.start + 0.5) / (m.sse.end - m.sse.start + 1)
            oi = np.add(off, m.direction * (f - 0.5) * m.length)
            side = 1.0 if sides[(n, m.sse.chain)] >= 0 else -1.0  # its chain's side: its loops cross no band
            target = _beyond(c, a, p + side * nrm * c.band / 2, side * nrm)
            gaps = [math.hypot(stems[q].cx - p[0], stems[q].cy - p[1]) for q in c.pieces]
            j = item_of[c.pieces[int(np.argmin(gaps))]]  # held by the piece nearest
            oj = target - np.array([items[j].cx, items[j].cy])
            wt = _RUNG_W * _W_CONTACT * min(4 * len(res), _CONTACT_CAP) / _CONTACT_CAP
            springs.append((k, _snap(oi), j, _snap(oj), float(_snap(wt)), float(_snap(m.width / 2 + MARGIN))))
            q = np.add(items[k].home, oi)  # is it projected on its side, clear of the band?
            v = float((q - target) @ nrm) * side + c.band / 2
            reach = c.band / 2 + max(items[k].half) + MARGIN
            if v < reach:  # no: started beside the band at its rungs, on its side (the walls push it off the band
                # the way they find, its loops then crossing it)
                start.setdefault(k, []).append(target + side * nrm * (reach - c.band / 2) - oi)
    for k, at in sorted(start.items()):
        x, y = (float(v) for v in _snap(np.mean(at, axis=0)))
        items[k].cx, items[k].cy = x, y
        items[k].home = (x, y)
    _to_its_side(items, elements, curves, sides, set(start) | {s[0] for s in springs})
    return curves, walls, springs


_PART_STEP, _PART_TRIES = 0.25, 400  # page units a curve moves at a time while parted from another, at most this often


def _centreline(c: NACurve) -> np.ndarray:
    """The curve's centreline about every half page unit (for gaps between curves)."""
    return c.xy[:: max(1, int(round(0.5 / CURVE_STEP)))]


def _part(curves: list[NACurve], stems: list[PlacedStem]) -> None:
    """Each curve after the first moved whole (its pieces with it), away from the middle of each earlier curve its band
    comes within MARGIN of, by what it lacks (at least _PART_STEP) at a time, until it keeps clear of all of them: each
    curve is laid out on its own (_chain, _ring) and then pinned, so nothing else would part two (3D0U's P3 and P4
    drawn over one another)."""
    for i, c in enumerate(curves):
        for _ in range(_PART_TRIES):
            near = None
            for e in curves[:i]:
                a, b = _centreline(c), _centreline(e)
                gap = float(np.hypot(a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1]).min())
                need = (c.band + e.band) / 2 + MARGIN - gap
                if need > 1e-9:
                    near = (e, need)
                    break
            if near is None:
                break
            e, need = near
            way = c.xy.mean(axis=0) - e.xy.mean(axis=0)
            size = float(np.hypot(*way))
            way = way / size if size > 1e-9 else np.array([1.0, 0.0])
            _shift(c, [stems[k] for k in c.pieces], _snap(way * max(_PART_STEP, need)))


def _shift(c: NACurve, ps: list[PlacedStem], v) -> None:
    """The curve and its pieces moved by v, whole."""
    c.xy = _snap(c.xy + v)
    for d in ps:
        d.cx, d.cy = float(_snap(d.cx + v[0])), float(_snap(d.cy + v[1]))
    c.held = [(d.cx, d.cy, d.angle) for d in ps]
    c.drawn.clear()


def _beyond(c: NACurve, a: float, edge: np.ndarray, way: np.ndarray) -> np.ndarray:
    """`edge` (the band's edge at arc length a, `way` the unit direction out of it), or, where another stretch of
    the band lies beyond it that way (a gyre inside the one at a), the edge of the last such stretch: the nearest
    place an element touching the duplex at a can stand (the octamer stands inside the inner gyre)."""
    other = np.abs(c.s - a) > 2 * c.band
    if not other.any():
        return edge
    pts = c.xy[other]
    out = edge
    for t in np.arange(0.0, 3 * c.band, 0.25):
        q = edge + t * way
        if float(np.hypot(pts[:, 0] - q[0], pts[:, 1] - q[1]).min()) < c.band / 2:
            out = q + way * 0.25
    return _snap(out)


def _side_of(c: NACurve, q) -> tuple[float, np.ndarray, np.ndarray]:
    """How far q lies to the left of the curve (negative: right), and the curve's point and left normal nearest it."""
    lo, hi = c.span()
    p, _, n = c.frames(np.linspace(lo, hi, max(2, math.ceil((hi - lo) / CURVE_STEP) + 1)))
    near = int(np.argmin(np.hypot(p[:, 0] - q[0], p[:, 1] - q[1])))
    return float(_snap((np.asarray(q) - p[near]) @ n[near])), p[near], n[near]


def _to_its_side(items: list, elements: list, curves: list[NACurve], sides: dict, bound: set[int]) -> None:
    """An item no rung pulls (`bound`: those some do), projected across a curve from its chain, whose elements'
    neighbours in sequence outside it are all on its chain's side, is started there instead: its place mirrored
    across the curve (at least clear of the band), so its loops need not cross it (1YTB's C-terminal helix)."""
    by_chain: dict[str, list] = {}
    for k, m, _, _ in elements:
        by_chain.setdefault(m.sse.chain, []).append((m.sse.start, m.sse.end, k))
    for rows in by_chain.values():
        rows.sort(key=lambda r: (r[0], r[1]))
    mine: dict[int, list] = {}
    for k, m, _, _ in elements:
        mine.setdefault(k, []).append(m)
    for k in sorted(mine):
        if k in bound:
            continue
        it = items[k]
        for n, c in enumerate(curves):
            chains = {m.sse.chain for m in mine[k]}
            if len(chains) != 1 or (n, next(iter(chains))) not in sides:
                continue
            chain = next(iter(chains))
            want = 1.0 if sides[(n, chain)] >= 0 else -1.0
            v, p, nrm = _side_of(c, it.home)
            if v * want >= 0:
                continue
            rows = by_chain[chain]
            here = [i for i, r in enumerate(rows) if r[2] == k]
            near = [rows[i] for i in (here[0] - 1, here[-1] + 1) if 0 <= i < len(rows) and rows[i][2] != k]
            if not near or any(_side_of(c, items[r[2]].home)[0] * want < 0 for r in near):
                continue
            reach = c.band / 2 + max(it.half) + MARGIN
            x, y = (float(u) for u in _snap(np.asarray(it.home) + (max(abs(v), reach) + abs(v)) * want * nrm))
            it.cx, it.cy = x, y
            it.home = (x, y)
            break


def _fractions(points: np.ndarray, at: np.ndarray) -> np.ndarray:
    """How far along the polyline (0 to 1, by length) each of the points `at` lies (its nearest point on it)."""
    seg = np.diff(points, axis=0)
    lens = np.hypot(seg[:, 0], seg[:, 1])
    total = float(lens.sum())
    if total <= 0 or not len(at):
        return np.zeros(len(at))
    cum = np.concatenate([[0.0], np.cumsum(lens)])
    out = []
    for q in at:
        rel = q - points[:-1]
        u = np.clip((rel[:, 0] * seg[:, 0] + rel[:, 1] * seg[:, 1]) / np.maximum(lens**2, 1e-12), 0.0, 1.0)
        foot = points[:-1] + u[:, None] * seg
        k = int(np.argmin(np.hypot(foot[:, 0] - q[0], foot[:, 1] - q[1])))
        out.append((cum[k] + u[k] * lens[k]) / total)
    return np.array(out)


def _hairpin(curves: list[NACurve], p) -> None:
    """A hairpin loop over a curve's tip drawn again as route_nucleic draws one over a stem's (na_runs._arc): on a
    straight stand-in for the curve's end there, its nucleotides spaced along it."""
    from .na_paths import _spaced
    from .na_runs import _arc

    ka, kb = p.between
    c = next(c for c in curves if ka in c.at)
    tip = (c.at[ka] + c.at[kb]) / 2
    centre, t, _ = (v[0] for v in c.frames(tip))
    stand = PlacedStem(
        "",
        "",
        [ka[0], kb[0]],
        [],
        float(centre[0]),
        float(centre[1]),
        0.0,
        float(_snap(math.atan2(t[1], t[0]))),
        c.width,
        {ka: c.at[ka] - tip, kb: c.at[kb] - tip},
        {ka: c.side[ka], kb: c.side[kb]},
    )
    m = len(p.nts)
    p.points = _arc(stand, ka, kb, m)
    p.nt_xy = _spaced(p.points, m)


def follow(layout: Layout, curves: list[NACurve], paths: list, ends: list) -> list:
    """Move the strands that leave or enter a curved duplex with its nucleotides: a path between two stem
    nucleotides moves by the first's displacement at its start and the second's at its end, blended along its
    length (none for a nucleotide not on a curve); a free end moves whole with the nucleotide it leaves; each strand
    end's label with its end nucleotide. Returns the ends moved."""
    old, new = {}, {}
    for c in curves:
        for k in c.pieces:
            d = layout.stems[k]
            for key in d.axial:
                old[key] = d.point(key)
                new[key] = c.point(key)
    disp = {key: new[key] - old[key] for key in new}
    zero = np.zeros(2)
    moved: dict[tuple[int, int], np.ndarray] = dict(disp)
    for p in paths:
        if p.kind == "joint" and p.between and p.between[0] in new:
            continue  # drawn as the curve
        pts = np.asarray(p.points, float)
        if not len(pts):
            continue
        if p.kind == "hairpin" and len(p.between) == 2 and all(k in new for k in p.between):
            _hairpin(curves, p)  # drawn again over the curve's tip: its two ends may have changed sides
            continue
        if p.between:
            da, db = disp.get(p.between[0]), disp.get(p.between[-1])
            if p.kind == "tail" or len(p.between) == 1:
                db = da
        else:  # a cropped gap: by the curve nucleotides its ends were drawn from
            near = lambda q: next((disp[k] for k in old if math.hypot(*(old[k] - q)) < 1.0), None)  # noqa: E731
            da, db = near(pts[0]), near(pts[-1])
        if da is None and db is None:
            continue
        da, db = (zero if da is None else da), (zero if db is None else db)
        f = _fractions(pts, pts)[:, None]
        p.points = [tuple(float(v) for v in q) for q in _snap(pts + (1 - f) * da + f * db)]
        if p.nt_xy:
            q = np.asarray(p.nt_xy, float)
            g = _fractions(pts, q)[:, None]
            shift = (1 - g) * da + g * db
            p.nt_xy = [tuple(float(v) for v in x) for x in _snap(q + shift)]
            for i, sh in zip(p.nts, shift):
                moved[(p.strand, i)] = sh
    na = layout.nucleic
    out = []
    for s, e, xy in ends:
        key = (s, 0) if e == "5′" else (s, len(na.strands[s]) - 1)
        sh = moved.get(key)
        out.append((s, e, xy if sh is None else tuple(float(v) for v in _snap(np.add(xy, sh)))))
    return out


def widen(layout: Layout, curves: list[NACurve]) -> None:
    """Grow layout.bounds to hold each curve's letters and name."""
    if not curves:
        return
    gw, gh = glyph(layout)
    pad = max(gw, gh) / 2 + 0.05
    pts = [c.point(key, NT_ROW) for c in curves for key in c.at]
    reach = [max(gw * len(c.name) / 2, gh) + 0.05 for c in curves]
    boxes = [(x - pad, y - pad, x + pad, y + pad) for x, y in pts]
    boxes += [(c.name_xy[0] - h, c.name_xy[1] - h, c.name_xy[0] + h, c.name_xy[1] + h) for c, h in zip(curves, reach)]
    b = np.asarray(boxes, float)
    x0, y0, x1, y1 = layout.bounds
    layout.bounds = (
        min(x0, float(b[:, 0].min())),
        min(y0, float(b[:, 1].min())),
        max(x1, float(b[:, 2].max())),
        max(y1, float(b[:, 3].max())),
    )


def bend(layout: Layout, paths: list, ends: list) -> list:
    """Style.na_bend curved: layout.na_curves set, the strands leaving them moved with them (follow), the joints
    between their pieces taken out of `paths` (drawn as the curve: what comes after, and every check, reads the
    strands as drawn), the page grown to hold them. Returns the strand ends moved."""
    layout.na_curves = build_curves(layout)
    ends = follow(layout, layout.na_curves, paths, ends)
    held = {key for c in layout.na_curves for key in c.at}
    paths[:] = [p for p in paths if not (p.kind == "joint" and p.between and p.between[0] in held)]
    widen(layout, layout.na_curves)
    return ends


# ---- the curve as drawn, for what reads the stems' bodies (checks, keep-outs): plan 25 Task 3 ----


@dataclass(eq=False)
class Chord(PlacedStem):
    """A stretch of a curve's body (NACurve: its rails' width, end to end) as a straight stand-in, what checks and
    keep-outs read in place of the curve's pieces (which are not drawn). id: the piece nearest it along the curve
    (its owner, as route._oid names it); curve: the curve's index in layout.na_curves; s0, s1: its arc lengths."""

    curve: int = -1
    s0: float = 0.0
    s1: float = 0.0
    band: float = 0.0  # the curve's band (NACurve.band)


def _owner(c: NACurve, stems: list[PlacedStem], s: float) -> PlacedStem:
    """The piece of the curve whose nucleotides' arc lengths come nearest s."""
    best, out = math.inf, stems[c.pieces[0]]
    for k in c.pieces:
        at = [c.at[key] for key in stems[k].axial if key in c.at]
        if not at:
            continue
        gap = max(min(at) - s, s - max(at), 0.0)
        if gap < best:
            best, out = gap, stems[k]
    return out


def chords(c: NACurve, n: int, stems: list[PlacedStem]) -> list[Chord]:
    """The curve's body as Chords, each at most half its rail offset long along the curve (a turn of under 10 degrees
    from one to the next at the tightest bend, _R_MIN), each grown past both ends by what closes the wedge its turn
    from its neighbour opens on the bend's outside, and across by how far the curve bows off it there (its sagitta,
    each side): one band, no gap, nothing past it but those hairs. Kept on the curve, whose geometry is fixed once
    made (keepers of the result must not change them)."""
    key = (n, tuple(stems[k].id for k in c.pieces))
    if key in c.drawn:
        return c.drawn[key]
    lo, hi = c.extent()
    m = max(1, math.ceil((hi - lo) / max(c.rail / 2, 0.05)))
    edges = np.linspace(lo, hi, m + 1)
    p = c.frames(edges)[0]
    seg = np.diff(p, axis=0)
    ang = np.arctan2(seg[:, 1], seg[:, 0])
    turn = np.abs((np.diff(ang) + np.pi) % (2 * np.pi) - np.pi)
    turn = np.maximum(np.concatenate([[0.0], turn]), np.concatenate([turn, [0.0]]))  # at either end, the larger
    grow = c.width / 2 * np.tan(np.minimum(turn, np.pi / 2) / 2)
    keys = sorted(c.at, key=lambda key: (c.at[key], key))
    out = []
    for j in range(m):
        s0, s1 = float(edges[j]), float(edges[j + 1])
        mid = (p[j] + p[j + 1]) / 2
        bow = c.xy[(c.s > s0) & (c.s < s1)]  # the curve between the ends: how far it bows off the chord
        across = np.array([-seg[j, 1], seg[j, 0]]) / max(float(np.hypot(*seg[j])), 1e-12)
        sag = float(np.abs((bow - mid) @ across).max()) if len(bow) else 0.0
        mine = [key for key in keys if s0 <= c.at[key] < s1 or (j == m - 1 and c.at[key] == s1)]
        owner = _owner(c, stems, (s0 + s1) / 2)
        out.append(
            Chord(
                owner.id,
                c.name,
                list(owner.strands),
                [],
                float(mid[0]),
                float(mid[1]),
                float(np.hypot(*seg[j])) + 2 * float(grow[j]),
                float(ang[j]),
                c.width + 2 * sag,
                {key: c.at[key] - (s0 + s1) / 2 for key in mine},
                {key: c.side[key] for key in mine},
                parent=c.parent,
                curve=n,
                s0=s0,
                s1=s1,
                band=c.band,
            )
        )
    c.drawn[key] = out
    return out


def on_curves(layout) -> set[int]:
    """The indices in layout.stems of the pieces drawn as curves (layout.na_curves)."""
    return {k for c in getattr(layout, "na_curves", ()) for k in c.pieces}


def bodies(layout) -> list[PlacedStem]:
    """The stems' bodies as drawn: every stem drawn straight, then each curve's Chords (in place of its pieces)."""
    curves = getattr(layout, "na_curves", ())
    if not curves:
        return list(layout.stems)
    hidden = on_curves(layout)
    out = [d for k, d in enumerate(layout.stems) if k not in hidden]
    for n, c in enumerate(curves):
        out += chords(c, n, layout.stems)
    return out


def drawn_pieces(layout) -> list[PlacedStem]:
    """layout.stems as drawn: each piece of a curve as a straight stand-in where its nucleotides are drawn on it (its
    centre the curve's point midway between its outermost rungs, its direction from the one to the other, its rails
    on the curve's sides), every other stem itself. For the measures that read stems one by one (nacheck.shape,
    nacheck.contacts)."""
    curves = getattr(layout, "na_curves", ())
    if not curves:
        return list(layout.stems)
    curve_of = {k: c for c in curves for k in c.pieces}
    out = []
    for k, d in enumerate(layout.stems):
        c = curve_of.get(k)
        if c is None:
            out.append(d)
            continue
        at = {key: c.at[key] for key in d.axial if key in c.at}
        lo, hi = min(at.values()), max(at.values())
        mid = (lo + hi) / 2
        a, b = c.to_page([lo, hi], [0.0, 0.0])
        cx, cy = (float(v) for v in c.to_page(mid, 0.0)[0])
        angle = float(math.atan2(b[1] - a[1], b[0] - a[0])) if hi > lo else c.angle(next(iter(at)))
        arm = next(key for key in sorted(d.axial) if d.side[key] == d.first)
        out.append(
            replace(
                d,
                cx=cx,
                cy=cy,
                angle=angle,
                length=hi - lo + DNA_RISE * SCALE,
                axial={key: v - mid for key, v in at.items()},
                side={key: c.side[key] for key in at},
                first=c.side[arm],
            )
        )
    return out


def label_box(c: NACurve, k: float, keep: bool = False) -> tuple[float, float, float, float]:
    """The box round the curve's name at its place (_name_place): with `keep` what routes keep out of, its ink grown
    by _INK_MARGIN em (name_keep, as a stem's stem_keep_box); else its text box, glyphs `k` (glyph_scale) times
    _GLYPH_W x _GLYPH_H, as layout._stem_label_box sizes a stem's."""
    if keep and c.name_keep is not None:
        return c.name_keep
    along, across = _GLYPH_W * k * len(c.name) / 2 + 0.05, _GLYPH_H * k / 2
    cos, sin = abs(math.cos(c.name_angle)), abs(math.sin(c.name_angle))
    return centred_box(*c.name_xy, cos * along + sin * across, sin * along + cos * across)


def piece_bodies(layout) -> list[list[PlacedStem]]:
    """Each of layout.stems' bodies as drawn, for the measures that read a stem's body (nacheck.contacts): the stem
    itself; a piece of a curve, the curve's stretches (chords) over its nucleotides, from half a base-pair step before
    its first rung to half a step past its last."""
    curves = getattr(layout, "na_curves", ())
    out = [[d] for d in layout.stems]
    half = DNA_RISE * SCALE / 2
    for n, c in enumerate(curves):
        mine = chords(c, n, layout.stems)
        for k in c.pieces:
            at = [c.at[key] for key in layout.stems[k].axial if key in c.at]
            lo, hi = min(at) - half, max(at) + half
            out[k] = [b for b in mine if b.s1 > lo and b.s0 < hi]
    return out
