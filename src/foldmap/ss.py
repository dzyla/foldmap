# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Turn per-residue DSSP codes into helix/strand SSE objects with an axis and centroid."""

from __future__ import annotations

import numpy as np

from .model import SSE, Backbone, Bridge

_MIN_STRAND = 2
_BULGE_GAP = 1  # strand runs this many residues apart are one strand with a beta-bulge
_MIN_HELIX_H = 4  # an alpha run needs this many H residues ...
_MIN_HELIX_ANY = 6  # ... or this many helix-class residues in total (3-10 / pi)
_MIN_310 = 3  # an isolated 3-10 run this long becomes its own small element (kind G) when asked for
_PARTNER_GAP = 4  # bridge partners this close along one chain belong to the same partner strand


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


def _same_partner(a: tuple[int, int], b: tuple[int, int], partners: dict[int, set[int]], bb: Backbone) -> bool:
    """True when strand runs a and b pair with the same partner strand (so they are one strand around a bulge)."""
    pa = {p for k in range(a[0], a[1] + 1) for p in partners.get(k, ())}
    pb = {p for k in range(b[0], b[1] + 1) for p in partners.get(k, ())}
    return any(abs(u - v) <= _PARTNER_GAP and bb.labels[u].chain == bb.labels[v].chain for u in pa for v in pb)


def _merge_bulges(runs, ss: str, bb: Backbone, partners: dict[int, set[int]] | None = None):
    """Join strand runs separated by a single non-helix residue (a beta-bulge, not a hairpin turn).
    With partners (from the DSSP bridges), only runs that pair with the same partner strand are joined."""
    prev = bb.prev
    merged: list[list[int]] = []
    for a, b in runs:
        if merged:
            last = merged[-1]
            gap = a - last[1] - 1
            if (
                gap == _BULGE_GAP
                and prev[last[1] + 1]
                and prev[a]
                and ss[last[1] + 1] in "-BT"
                and (partners is None or _same_partner((last[0], last[1]), (a, b), partners, bb))
            ):
                last[1] = b
                continue
        merged.append([a, b])
    return [tuple(r) for r in merged]


def _fit_axis(ca: np.ndarray) -> np.ndarray:
    centred = ca - ca.mean(axis=0)
    axis = np.linalg.svd(centred, full_matrices=False)[2][0]
    head, tail = ca[: min(2, len(ca))].mean(axis=0), ca[-min(2, len(ca)) :].mean(axis=0)
    if np.dot(tail - head, axis) < 0:
        axis = -axis
    return axis / np.linalg.norm(axis)


def build_sses(bb: Backbone, ss: str, short_helices: bool = False, bridges: list[Bridge] | None = None) -> list[SSE]:
    """Strands (E) and helices (H); with short_helices, isolated 3-10 runs too short to count as helices
    become their own kind G (drawn as small boxes, labelled η). Pass the DSSP bridges so that strand runs one
    residue apart are joined only when they pair with the same strand."""
    partners = _partners(bridges) if bridges is not None else None
    out: list[SSE] = []
    for kind, classes in (("E", "E"), ("H", "HGI")):
        runs = list(_runs(ss, classes, bb.prev))
        if kind == "E":
            runs = _merge_bulges(runs, ss, bb, partners)
        for a, b in runs:
            n = b - a + 1
            if kind == "E" and n < _MIN_STRAND:
                continue
            if kind == "H" and not (ss[a : b + 1].count("H") >= _MIN_HELIX_H or n >= _MIN_HELIX_ANY):
                continue
            ca = bb.ca[a : b + 1]
            out.append(SSE(kind, bb.labels[a].chain, a, b, bb.labels[a], bb.labels[b], ca.mean(axis=0), _fit_axis(ca)))
    if short_helices:
        taken = {k for s in out for k in range(s.start, s.end + 1)}
        for a, b in _runs(ss, "G", bb.prev):
            if b - a + 1 >= _MIN_310 and not taken & set(range(a, b + 1)):
                ca = bb.ca[a : b + 1]
                out.append(
                    SSE("G", bb.labels[a].chain, a, b, bb.labels[a], bb.labels[b], ca.mean(axis=0), _fit_axis(ca))
                )
    out.sort(key=lambda s: s.start)
    return out
