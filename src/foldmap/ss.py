"""Turn per-residue DSSP codes into helix/strand SSE objects with an axis and centroid."""

from __future__ import annotations

import numpy as np

from .model import SSE, Backbone

_MIN_STRAND = 2
_BULGE_GAP = 1  # strand runs this many residues apart are one strand with a beta-bulge
_MIN_HELIX_H = 4  # an alpha run needs this many H residues ...
_MIN_HELIX_ANY = 6  # ... or this many helix-class residues in total (3-10 / pi)
_MIN_310 = 3  # an isolated 3-10 run this long becomes its own small element (kind G) when asked for


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


def _merge_bulges(runs, ss: str, prev: np.ndarray):
    """Join strand runs separated by a single non-helix residue (a beta-bulge, not a hairpin turn)."""
    merged: list[list[int]] = []
    for a, b in runs:
        if merged:
            last = merged[-1]
            gap = a - last[1] - 1
            if gap == _BULGE_GAP and prev[last[1] + 1] and prev[a] and ss[last[1] + 1] in "-B":
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


def build_sses(bb: Backbone, ss: str, short_helices: bool = False) -> list[SSE]:
    """Strands (E) and helices (H); with short_helices, isolated 3-10 runs too short to count as helices
    become their own kind G (drawn as small boxes, labelled η)."""
    out: list[SSE] = []
    for kind, classes in (("E", "E"), ("H", "HGI")):
        runs = list(_runs(ss, classes, bb.prev))
        if kind == "E":
            runs = _merge_bulges(runs, ss, bb.prev)
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
