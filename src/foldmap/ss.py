# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Turn per-residue DSSP codes into helix/strand SSE objects with an axis and centroid."""

from __future__ import annotations

import re

import numpy as np

from .model import SSE, Backbone, Bridge
from .sheets import _MIN_PAIRS, _SPLIT_GAP, _edge_table

_MIN_STRAND = 2
_MIN_HELIX_H = 4  # an alpha run needs this many H residues ...
_MIN_HELIX_ANY = 6  # ... or this many helix-class residues in total (3-10 / pi)
_MIN_310 = 3  # an isolated 3-10 run this long becomes its own small element (kind G) when asked for
_PARTNER_GAP = 4  # residues of a partner strand this close pair with the same side of it
KINK_SPLIT = 30.0  # degrees: alpha segments of one helix run meeting at a 3-10/pi kink sharper than this split it


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


def _bulged(i: int, strands: list[SSE], table: dict, ss: str, bb: Backbone, partners: dict[int, set[int]]) -> bool:
    """True when strands i and i + 1 are pieces of one strand split at a beta-bulge: 1-4 residues apart in one
    unbroken chain, with no helix or strand residue between, not paired with each other (a hairpin), and paired
    with the same partner strand on the same side of it (disjoint partner residues, in register)."""
    x, y = strands[i], strands[i + 1]
    gap = range(x.end + 1, y.start)
    if x.chain != y.chain or not 1 <= len(gap) <= _SPLIT_GAP or (i, i + 1) in table:
        return False
    if not all(bb.prev[k] for k in range(x.end + 1, y.start + 1)) or any(ss[k] in "HGIE" for k in gap):
        return False
    edge = lambda a, b: (min(a, b), max(a, b)) in table  # noqa: E731
    for k in (k for k in range(len(strands)) if edge(i, k) and edge(i + 1, k)):
        px, py = _pairing(x, strands[k], partners), _pairing(y, strands[k], partners)
        if px and py and not px & py and min(abs(u - v) for u in px for v in py) <= _PARTNER_GAP:
            return True
    return False


def _join_bulged_pieces(sses: list[SSE], ss: str, bb: Backbone, bridges: list[Bridge]) -> list[SSE]:
    """Join strand pieces split at a beta-bulge into one strand (see _bulged), until nothing joins."""
    partners = _partners(bridges)
    while True:
        strands = sorted((s for s in sses if s.kind == "E"), key=lambda s: s.start)
        table = _edge_table(strands, bridges, _MIN_PAIRS)
        i = next((i for i in range(len(strands) - 1) if _bulged(i, strands, table, ss, bb, partners)), None)
        if i is None:
            return sses
        x, y = strands[i], strands[i + 1]
        ca = bb.ca[x.start : y.end + 1]
        joined = SSE("E", x.chain, x.start, y.end, x.first, y.last, ca.mean(axis=0), _fit_axis(ca))
        sses = [joined if s is x else s for s in sses if s is not y]


def _fit_axis(ca: np.ndarray) -> np.ndarray:
    centred = ca - ca.mean(axis=0)
    axis = np.linalg.svd(centred, full_matrices=False)[2][0]
    head, tail = ca[: min(2, len(ca))].mean(axis=0), ca[-min(2, len(ca)) :].mean(axis=0)
    if np.dot(tail - head, axis) < 0:
        axis = -axis
    return axis / np.linalg.norm(axis)


def _angle(u: np.ndarray, v: np.ndarray) -> float:
    return float(np.degrees(np.arccos(np.clip(u @ v, -1.0, 1.0))))


def _split_kinks(a: int, b: int, ss: str, bb: Backbone) -> tuple[list[tuple[int, int]], list[int]]:
    """The pieces of helix run a..b cut where two of its alpha segments (H runs of _MIN_HELIX_H or more) meet
    through 3-10/pi residues only at more than KINK_SPLIT degrees (the globin AB corner), and the corner residues
    between them, which belong to neither piece."""
    hs = [(m.start() + a, m.end() - 1 + a) for m in re.finditer(f"H{{{_MIN_HELIX_H},}}", ss[a : b + 1])]
    pieces, corners, start = [], [], a
    for (x0, x1), (y0, y1) in zip(hs, hs[1:]):
        if (
            set(ss[x1 + 1 : y0]) <= set("GI")
            and _angle(_fit_axis(bb.ca[x0 : x1 + 1]), _fit_axis(bb.ca[y0 : y1 + 1])) > KINK_SPLIT
        ):
            pieces.append((start, x1))
            corners.extend(range(x1 + 1, y0))
            start = y0
    pieces.append((start, b))
    return pieces, corners


def build_sses(bb: Backbone, ss: str, short_helices: bool = False, bridges: list[Bridge] | None = None) -> list[SSE]:
    """Strands (E) and helices (H); with short_helices, isolated 3-10 runs too short to count as helices
    become their own kind G (drawn as small boxes, labelled η). Pass the DSSP bridges so that strand pieces
    split at a beta-bulge are joined into one strand (without them, strands stay as DSSP runs). A helix run
    whose alpha segments turn sharply at a 3-10/pi kink is split there (_split_kinks)."""
    out: list[SSE] = []
    corners: list[int] = []  # residues of a kink between two helix pieces: no element of their own
    for kind, classes in (("E", "E"), ("H", "HGI")):
        for a, b in _runs(ss, classes, bb.prev):
            n = b - a + 1
            if kind == "E" and n < _MIN_STRAND:
                continue
            if kind == "H" and not (ss[a : b + 1].count("H") >= _MIN_HELIX_H or n >= _MIN_HELIX_ANY):
                continue
            pieces = [(a, b)]
            if kind == "H":
                pieces, cut = _split_kinks(a, b, ss, bb)
                corners += cut
            for p, q in pieces:
                ca = bb.ca[p : q + 1]
                out.append(
                    SSE(kind, bb.labels[p].chain, p, q, bb.labels[p], bb.labels[q], ca.mean(axis=0), _fit_axis(ca))
                )
    if bridges is not None:
        out = _join_bulged_pieces(out, ss, bb, bridges)
    if short_helices:
        taken = {k for s in out for k in range(s.start, s.end + 1)} | set(corners)
        for a, b in _runs(ss, "G", bb.prev):
            if b - a + 1 >= _MIN_310 and not taken & set(range(a, b + 1)):
                ca = bb.ca[a : b + 1]
                out.append(
                    SSE("G", bb.labels[a].chain, a, b, bb.labels[a], bb.labels[b], ca.mean(axis=0), _fit_axis(ca))
                )
    out.sort(key=lambda s: s.start)
    return out
