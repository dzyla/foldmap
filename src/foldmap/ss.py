# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Turn per-residue DSSP codes into helix/strand SSE objects with an axis and centroid."""

from __future__ import annotations

import functools
import logging
import re

import numpy as np

from .model import SSE, Backbone, Bridge, homolog_chains, same_place
from .sheets import _MIN_PAIRS, _SPLIT_GAP, _edge_table

log = logging.getLogger(__name__)

_MIN_STRAND = 2
_MIN_HELIX_H = 4  # an alpha run needs this many H residues ...
_MIN_HELIX_ANY = 6  # ... or this many helix-class residues in total (3-10 / pi)
_MIN_310 = 3  # an isolated 3-10 run this long becomes its own small element (kind G) when asked for
_PARTNER_GAP = 4  # residues of a partner strand this close pair with the same side of it
KINK_SPLIT = 30.0  # degrees: alpha segments of one helix run meeting at a 3-10/pi kink sharper than this split it
_MIN_KINK_ARM = 7  # residues: an alpha segment shorter than this has no axis precise enough to measure a kink by


def _runs(ss: str, classes: str, prev: np.ndarray):
    """Maximal chain-contiguous runs of residues whose code is in `classes`."""
    start = None
    for i, c in enumerate(ss):
        if c in classes:
            if start is None:
                start = i
            elif not prev[i]:
                yield start, i - 1
                start = i
        elif start is not None:
            yield start, i - 1
            start = None
    if start is not None:
        yield start, len(ss) - 1


def _partners(bridges: list[Bridge]) -> dict[int, set[int]]:
    out: dict[int, set[int]] = {}
    for i, j, _ in bridges:
        out.setdefault(i, set()).add(j)
        out.setdefault(j, set()).add(i)
    return out


def _pairing(a: SSE, k: SSE, partners: dict[int, set[int]]) -> set[int]:
    """Residues of strand k bridged to strand a."""
    return {p for r in range(a.start, a.end + 1) for p in partners.get(r, ()) if k.start <= p <= k.end}


def _bulged(
    i: int, strands: list[SSE], table: dict, ss: str, bb: Backbone, partners: dict[int, set[int]], max_gap=_SPLIT_GAP
) -> bool:
    """True when strands i and i + 1 are pieces of one strand split at a beta-bulge: 1-4 (max_gap) residues apart
    in one unbroken chain, with no helix or strand residue between, not paired with each other (a hairpin), and
    paired with the same partner strand on the same side of it (disjoint partner residues, in register)."""
    x, y = strands[i], strands[i + 1]
    gap = range(x.end + 1, y.start)
    if x.chain != y.chain or not 1 <= len(gap) <= max_gap or (i, i + 1) in table:
        return False
    if not all(bb.prev[k] for k in range(x.end + 1, y.start + 1)) or any(ss[k] in "HGIE" for k in gap):
        return False
    edge = lambda a, b: (min(a, b), max(a, b)) in table  # noqa: E731
    for k in (k for k in range(len(strands)) if edge(i, k) and edge(i + 1, k)):
        px, py = _pairing(x, strands[k], partners), _pairing(y, strands[k], partners)
        if px and py and not px & py and min(abs(u - v) for u in px for v in py) <= _PARTNER_GAP:
            return True
    return False


def _span(bb: Backbone, copies: dict[str, tuple[str, int]], x: tuple[int, int], y: tuple[int, int]) -> tuple:
    """Two pieces (first and last residue index each) in the residue numbers of the first chain of their copies
    (model.homolog_chains): (that chain, x first, x last, y first, y last)."""
    root, shift = copies[bb.labels[x[0]].chain]
    n = lambda i: bb.labels[i].seq - shift  # noqa: E731
    return root, n(x[0]), n(x[1]), n(y[0]), n(y[1])


def _ends(sses: list[SSE]) -> list[tuple[int, int]]:
    return [(s.start, s.end) for s in sses]


def _alike(u: tuple, v: tuple) -> bool:
    """Are two pairs of pieces (_span) the same place in copies of one protein: x at x's place and y at y's
    (model.same_place)?"""
    return u[0] == v[0] and same_place(u[1:3], v[1:3]) and same_place(u[3:5], v[3:5])


def _join_bulged_pieces(sses: list[SSE], ss: str, bb: Backbone, bridges: list[Bridge], copies) -> list[SSE]:
    """Join strand pieces split at a beta-bulge into one strand (see _bulged), until nothing joins. Then join the
    matching pieces of every other copy of the protein (copies() as in model.homolog_chains) up to two residues
    further apart, so that identical chains are drawn alike (1EYG chain D's 4096-4102 + 4108-4111, 5 apart)."""
    partners = _partners(bridges)
    joined: list[tuple] = []  # the pairs of pieces joined so far, as _span
    for shared in (False, True):
        while True:
            strands = sorted((s for s in sses if s.kind == "E"), key=lambda s: s.start)
            table = _edge_table(strands, bridges, _MIN_PAIRS)

            pairs = range(len(strands) - 1)
            if shared:
                found = (
                    i
                    for i in pairs
                    if any(_alike(_span(bb, copies(), *_ends(strands[i : i + 2])), u) for u in joined)
                    and _bulged(i, strands, table, ss, bb, partners, max_gap=_SPLIT_GAP + 2)
                )
            else:
                found = (i for i in pairs if _bulged(i, strands, table, ss, bb, partners))
            i = next(found, None)
            if i is None:
                break
            x, y = strands[i], strands[i + 1]
            joined.append(_span(bb, copies(), *_ends([x, y])))
            sses = [_sse("E", bb, x.start, y.end) if s is x else s for s in sses if s is not y]
    return sses


def _fit_axis(ca: np.ndarray) -> np.ndarray:
    centred = ca - ca.mean(axis=0)
    axis = np.linalg.svd(centred, full_matrices=False)[2][0]
    head, tail = ca[: min(2, len(ca))].mean(axis=0), ca[-min(2, len(ca)) :].mean(axis=0)
    if np.dot(tail - head, axis) < 0:
        axis = -axis
    return axis / np.linalg.norm(axis)


def _sse(kind: str, bb: Backbone, a: int, b: int) -> SSE:
    ca = bb.ca[a : b + 1]
    return SSE(kind, bb.labels[a].chain, a, b, bb.labels[a], bb.labels[b], ca.mean(axis=0), _fit_axis(ca))


def _angle(u: np.ndarray, v: np.ndarray) -> float:
    return float(np.degrees(np.arccos(np.clip(u @ v, -1.0, 1.0))))


def _helix_axis(ca: np.ndarray) -> np.ndarray:
    """The local axis of an alpha segment, N -> C (Kahn 1989): at each inner CA the bisector of its two CA-CA
    vectors points at the axis, so the cross product of successive bisectors runs along it; their mean. A straight
    fit through the CAs of a short segment is off by up to 40 degrees."""
    b = (ca[:-2] - ca[1:-1]) + (ca[2:] - ca[1:-1])
    along = np.cross(b[:-1], b[1:])
    along /= np.linalg.norm(along, axis=1)[:, None]
    along *= np.where(along @ (ca[-1] - ca[0]) < 0, -1.0, 1.0)[:, None]
    m = along.mean(axis=0)
    return m / np.linalg.norm(m)


def _kinks(a: int, b: int, ss: str, bb: Backbone) -> list[tuple[tuple[int, int], tuple[int, int], float]]:
    """Consecutive alpha segments (H runs) of helix run a..b, both at least _MIN_KINK_ARM long, that meet through
    3-10/pi residues only, with the angle between their local axes in degrees (snapped to 1e-4, so noise in the
    fit never decides a comparison with KINK_SPLIT)."""
    hs = [(m.start() + a, m.end() - 1 + a) for m in re.finditer(f"H{{{_MIN_KINK_ARM},}}", ss[a : b + 1])]
    axis = lambda h: _helix_axis(bb.ca[h[0] : h[1] + 1])  # noqa: E731
    return [
        (x, y, round(_angle(axis(x), axis(y)), 4)) for x, y in zip(hs, hs[1:]) if set(ss[x[1] + 1 : y[0]]) <= set("GI")
    ]


def _split_kinks(a: int, b: int, cuts) -> tuple[list[tuple[int, int]], list[int]]:
    """The pieces of helix run a..b cut at each kink (x, y) in cuts, and the corner residues between x and y,
    which belong to neither piece."""
    pieces, corners, start = [], [], a
    for x, y in sorted(cuts):
        pieces.append((start, x[1]))
        corners.extend(range(x[1] + 1, y[0]))
        start = y[0]
    pieces.append((start, b))
    return pieces, corners


_BREAK_MAX = 2  # residues: a helix broken by this few in one copy of a protein is broken there in every copy


def _copy_breaks(helices: list[tuple[int, int]], bb: Backbone, copies) -> tuple[list[tuple[int, int]], int, int]:
    """Helix runs broken alike in every copy of a protein: where DSSP breaks a helix by up to _BREAK_MAX residues in
    some copies (two runs, chain-contiguous) and other copies' runs cover that place whole, the majority decides, a
    tie breaking: the whole runs are cut there, at the same residue numbers (copies() as in model.homolog_chains),
    the residues between in no piece, each piece still a helix's length (_MIN_HELIX_ANY); or the broken copies' two
    runs are joined into one (and a run broken twice so, all three). 7OKA: 225-231 + 233-238 in B and C, 225-238 in
    A: A is cut. One copy of many broken does not cut all the others (7JU4's tubulins). Returns the runs and how
    many cuts and joins were made (logged at INFO: the README tells the rule)."""
    pairs = [
        (x, y)
        for x, y in zip(helices, helices[1:])
        if bb.labels[x[0]].chain == bb.labels[y[0]].chain
        and 1 <= y[0] - x[1] - 1 <= _BREAK_MAX
        and bool(bb.prev[x[1] + 1 : y[0] + 1].all())
    ]
    if not pairs:
        return helices, 0, 0
    n_cut = n_joined = 0
    at = {(l.chain, l.seq): k for k, l in enumerate(bb.labels)}
    broken: dict[tuple, set[str]] = {}  # (root, last number before the break, first after) -> chains broken there
    for x, y in pairs:
        root, _, xb, ya, _ = _span(bb, copies(), x, y)
        broken.setdefault((root, xb, ya), set()).add(bb.labels[x[0]].chain)

    def ref(a: int, b: int) -> tuple[str, int, int]:
        root, shift = copies()[bb.labels[a].chain]
        return root, bb.labels[a].seq - shift, bb.labels[b].seq - shift

    whole: dict[tuple, set[str]] = {g: set() for g in broken}
    for a, b in helices:
        root, lo, hi = ref(a, b)
        for g in broken:
            if g[0] == root and lo <= g[1] and g[2] <= hi:
                whole[g].add(bb.labels[a].chain)
    cut = {g for g in broken if whole[g] and len(broken[g]) >= len(whole[g])}
    join = {g for g in broken if whole[g] and len(broken[g]) < len(whole[g])}
    following: dict[int, tuple[int, int]] = {}  # a broken run's start -> the run it is joined to
    for x, y in pairs:
        root, _, xb, ya, _ = _span(bb, copies(), x, y)
        if (root, xb, ya) in join:
            following[x[0]] = y
    out = []
    skip = set()
    for a, b in helices:
        if a in skip:
            continue
        if a in following:  # joined with the next run, and that with its next (x|y|z), as the other copies have it
            at = a
            while at in following:
                at, b = following[at]
                skip.add(at)
                n_joined += 1
            out.append((a, b))
            continue
        chain = bb.labels[a].chain
        root, lo, hi = ref(a, b)
        shift = copies()[chain][1]
        cuts = sorted((xb, ya) for r, xb, ya in cut if r == root and lo <= xb and ya <= hi and xb < ya)
        pieces, start = [], a
        for xb, ya in cuts:
            end, nxt = at.get((chain, xb + shift)), at.get((chain, ya + shift))
            if end is None or nxt is None or end < start or nxt > b:
                continue
            if end - start + 1 >= _MIN_HELIX_ANY and b - nxt + 1 >= _MIN_HELIX_ANY:
                pieces.append((start, end))
                start = nxt
                n_cut += 1
        pieces.append((start, b))
        out += pieces
    if n_cut or n_joined:
        log.info("helix pieces matched across copies: %d cut, %d joined", n_cut, n_joined)
    return out, n_cut, n_joined


def build_sses(
    bb: Backbone, ss: str, short_helices: bool = False, bridges: list[Bridge] | None = None, *, matched=None
) -> list[SSE]:
    """Strands (E) and helices (H); with short_helices, isolated 3-10 runs too short to count as helices
    become their own kind G (drawn as small boxes, labelled η). Pass the DSSP bridges so that strand pieces
    split at a beta-bulge are joined into one strand (without them, strands stay as DSSP runs). A helix run is
    split where two of its alpha segments meet through 3-10/pi residues at more than KINK_SPLIT degrees (the
    globin AB corner; _kinks). A kink split or a bulge join in one copy of a protein is made in its other copies
    too, whatever their angle and up to two residues further apart, so that identical chains are drawn alike; a helix
    DSSP breaks by one or two residues in some copies and not in others is drawn as most copies have it (_copy_breaks;
    `matched`, a dict, gets the numbers of cuts and joins)."""
    copies = functools.cache(lambda: homolog_chains(bb))  # computed only when some copy splits or joins
    out = [_sse("E", bb, a, b) for a, b in _runs(ss, "E", bb.prev) if b - a + 1 >= _MIN_STRAND]
    helices = [
        (a, b)
        for a, b in _runs(ss, "HGI", bb.prev)
        if ss[a : b + 1].count("H") >= _MIN_HELIX_H or b - a + 1 >= _MIN_HELIX_ANY
    ]
    helices, cut, joined = _copy_breaks(helices, bb, copies)
    if matched is not None:
        matched.update(cut=cut, joined=joined)
    kinks = {h: _kinks(*h, ss, bb) for h in helices}
    sharp = [_span(bb, copies(), x, y) for h in helices for x, y, angle in kinks[h] if angle > KINK_SPLIT]
    corners: list[int] = []  # residues of a kink between two helix pieces: no element of their own
    for a, b in helices:
        cuts = [(x, y) for x, y, _ in kinks[(a, b)] if any(_alike(_span(bb, copies(), x, y), u) for u in sharp)]
        pieces, cut = _split_kinks(a, b, cuts)
        corners += cut
        for k, (p, q) in enumerate(pieces):
            out.append(_sse("H", bb, p, q))
            out[-1].kinked = k > 0
    if bridges is not None:
        out = _join_bulged_pieces(out, ss, bb, bridges, copies)
    if short_helices:
        taken = {k for s in out for k in range(s.start, s.end + 1)} | set(corners)
        for a, b in _runs(ss, "G", bb.prev):
            if b - a + 1 >= _MIN_310 and not taken & set(range(a, b + 1)):
                out.append(_sse("G", bb, a, b))
    out.sort(key=lambda s: s.start)
    return out
