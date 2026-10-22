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
inside, as most of the pieces show it), parted by the band's width (_gyres). The pieces stay in layout.stems as
the layout placed them; render draws the curve in their place (layout.na_curves) and route_nucleic's strands that
leave or enter the duplex follow its nucleotides (follow)."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .layout import DNA_RISE, NT_ROW, SCALE, Layout, PlacedStem, _snap, glyph

CURVE_STEP = 0.1  # page units between the curve's samples
_FAR = 2.0  # two places on the band this many band widths apart along it are different stretches (_far_gap)
_GYRE_GAP = 0.2  # page units kept between two gyres' bands (their letters) beyond the band's width
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


def bent_duplexes(stems: list[PlacedStem]) -> list[tuple[str, list[int]]]:
    """(parent id, its pieces' indices in `stems` in order along the first arm 5'->3') of each duplex cut into two
    or more pieces."""
    groups: dict[str, list[int]] = {}
    for k, d in enumerate(stems):
        if d.parent:
            groups.setdefault(d.parent, []).append(k)

    def start(k: int) -> tuple[int, int]:  # the piece's first-arm nucleotide nearest the 5' end
        return min(key for key in stems[k].axial if stems[k].side[key] == stems[k].first)

    out = [(parent, sorted(ks, key=start)) for parent, ks in groups.items() if len(ks) > 1]
    return sorted(out, key=lambda g: g[1][0])


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
    ctrl = [np.array([ps[0].cx, ps[0].cy]) - ps[0].direction * ps[0].length / 2]
    ctrl += [np.array([d.cx, d.cy]) for d in ps]
    ctrl += [np.array([ps[-1].cx, ps[-1].cy]) + ps[-1].direction * ps[-1].length / 2]
    ctrl = _snap(np.array(ctrl))
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


def _name_place(c: NACurve, band: float, off: float) -> tuple[tuple[float, float], float]:
    """Where the duplex's one name goes: beside the curve, `off` from its axis, at the first of _NAME_AT along the
    duplex where no other stretch of the band comes nearer, on the side away from the middle of the curve first."""
    lo, hi = c.span()
    middle = c.xy.mean(axis=0)
    best = None
    for f in _NAME_AT:
        p, t, n = (v[0] for v in c.frames(lo + f * (hi - lo)))
        away = 1.0 if math.hypot(*(p + n - middle)) >= math.hypot(*(p - n - middle)) else -1.0
        for side in (away, -away):
            q = p + side * off * n
            room = float(np.hypot(c.xy[:, 0] - q[0], c.xy[:, 1] - q[1]).min())
            if best is None or room > best[0] + 1e-9:
                best = (room, q, t)
            if room >= off - 0.01:
                return tuple(float(v) for v in _snap(q)), float(_snap(math.atan2(t[1], t[0])))
    _, q, t = best
    return tuple(float(v) for v in _snap(q)), float(_snap(math.atan2(t[1], t[0])))


def build_curves(layout: Layout) -> list[NACurve]:
    """Every bent duplex of the layout as one curve (NACurve), its nucleotides placed on it."""
    stems = layout.stems
    na = layout.nucleic
    gw, gh = glyph(layout)
    out = []
    for parent, ks in bent_duplexes(stems):
        ps = [stems[k] for k in ks]
        rail, width = ps[0].rail, ps[0].width
        band = 2 * (rail + NT_ROW + gh)  # rails and base letters
        xy, centres, gyres, overlap = _curve_of(stems, ks, band)
        s = _arc(xy)
        tangent = _tangents(xy)
        tmp = NACurve(parent, "", ks, xy, s, tangent, rail, width, {}, {}, [])
        lengths = [d.length for d in ps]
        U = np.cumsum([0.0, *lengths[:-1]]) + np.array(lengths) / 2  # each centre along the straightened duplex
        shift = _isotonic(np.array(centres) - U, lengths)  # each piece's place along the curve less its own
        vote = 0.0  # the first arm's side of the curve, as the pieces show it (those seen most clearly most)
        for d, sk in zip(ps, centres):
            t = tmp.frames(sk)[1][0]
            along = float(t[0]) * float(d.direction[0]) + float(t[1]) * float(d.direction[1])
            vote += (d.margin + 0.01) * d.first * (1.0 if along >= 0 else -1.0)
        first_side = 1 if float(_snap(vote)) >= 0 else -1
        at, side, pairs = {}, {}, []
        for d, u, sh in zip(ps, U, shift):
            for key, x in d.axial.items():
                at[key] = float(_snap(u + sh + x))
                side[key] = first_side if d.side[key] == d.first else -first_side
            pairs += list(d.pairs)
        whole = next((st.name for st in na.structure.stems if st.id == parent), "") if na is not None else ""
        name = whole or ps[0].name.rstrip("abcdefghijklmnopqrstuvwxyz").strip()
        c = NACurve(parent, name, ks, xy, s, tangent, rail, width, at, side, pairs, gyres, overlap)
        c.name_xy, c.name_angle = _name_place(c, band, rail + NT_ROW + 1.5 * gh)
        out.append(c)
    return out


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
    """Style.na_bend curved: layout.na_curves set, the strands leaving them moved with them (follow), the page
    grown to hold them. Returns the strand ends moved."""
    layout.na_curves = build_curves(layout)
    ends = follow(layout, layout.na_curves, paths, ends)
    widen(layout, layout.na_curves)
    return ends
