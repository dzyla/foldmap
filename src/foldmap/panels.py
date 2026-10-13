# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Domain panel outlines that follow their elements (plan 10 Task 3).

A panel is the merged outline of its domain's element boxes, padded by a margin, with rounded corners, the gaps
between near elements filled (a morphological closing), and its name's box wrapped. Where elements of two domains
lie close (strands of two domains in one sheet, a chain's strand in another chain's sheet) the space between them
is split: a point belongs to the domain whose elements are nearer, a gutter kept between the two where there is
room for it, each element inside its own panel, and no panel past half way between two elements. So panels of
different domains never share area, and a panel never covers another domain's element or an element in no domain.

Each domain gets a field on one raster (GRID page units, anchored at the origin, so a figure always gets the same
raster): positive inside its panel, negative outside. render draws the zero contour; layout reads the cells inside
to place the name. Pure geometry: boxes in, fields out."""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

GRID = 0.1  # page units per raster cell
CLOSE = 2.0  # gaps between a domain's padded elements up to twice this are filled (2LZM: 1.2 left bays)
GUTTER = 0.3  # clear space kept between two domains' panels, where their elements leave room for it
NAME_PAD = 0.15  # a name's box is wrapped this far
EDGE = 0.05  # each element's own box, grown this far (never past half way to another's), is inside its panel

Box = tuple[float, float, float, float]
Turned = tuple[float, float, float, float, float]  # a turned box: centre x, y, half length, half width, angle


class Fields(NamedTuple):
    xs: np.ndarray  # cell centres along x
    ys: np.ndarray  # ... and y
    values: dict[str, np.ndarray]  # domain -> (len(ys), len(xs)) field, positive inside its panel


def _axis(lo: float, hi: float) -> np.ndarray:
    return (np.arange(math.floor(lo / GRID), math.ceil(hi / GRID)) + 0.5) * GRID


def _distance(boxes: list[Box], xs: np.ndarray, ys: np.ndarray, reach: float) -> np.ndarray:
    """Euclidean distance of every cell centre to the nearest of the boxes (inf beyond `reach` of all of them)."""
    out = np.full((len(ys), len(xs)), np.inf)
    for x0, y0, x1, y1 in boxes:
        i0, i1 = np.searchsorted(xs, x0 - reach), np.searchsorted(xs, x1 + reach, side="right")
        j0, j1 = np.searchsorted(ys, y0 - reach), np.searchsorted(ys, y1 + reach, side="right")
        if i0 >= i1 or j0 >= j1:
            continue
        gx, gy = xs[i0:i1], ys[j0:j1]
        dx = np.maximum(np.maximum(x0 - gx, gx - x1), 0.0)
        dy = np.maximum(np.maximum(y0 - gy, gy - y1), 0.0)
        np.minimum(out[j0:j1, i0:i1], np.hypot(dy[:, None], dx[None, :]), out=out[j0:j1, i0:i1])
    return out


def _turned_distance(bodies: list[Turned], xs: np.ndarray, ys: np.ndarray, reach: float) -> np.ndarray:
    """Euclidean distance of every cell centre to the nearest of the turned boxes (inf beyond `reach`)."""
    out = np.full((len(ys), len(xs)), np.inf)
    for cx, cy, hl, hw, angle in bodies:
        c, s = math.cos(angle), math.sin(angle)
        ex, ey = abs(c) * hl + abs(s) * hw + reach, abs(s) * hl + abs(c) * hw + reach
        i0, i1 = np.searchsorted(xs, cx - ex), np.searchsorted(xs, cx + ex, side="right")
        j0, j1 = np.searchsorted(ys, cy - ey), np.searchsorted(ys, cy + ey, side="right")
        if i0 >= i1 or j0 >= j1:
            continue
        gx, gy = xs[i0:i1][None, :] - cx, ys[j0:j1][:, None] - cy
        along = np.maximum(np.abs(gx * c + gy * s) - hl, 0.0)
        across = np.maximum(np.abs(-gx * s + gy * c) - hw, 0.0)
        np.minimum(out[j0:j1, i0:i1], np.hypot(along, across), out=out[j0:j1, i0:i1])
    return out


def panel_fields(
    owned: dict[str, list[Box]],
    free: list[Box],
    pad: float,
    names: dict[str, Box] | None = None,
    turned: list[Turned] = (),
) -> Fields:
    """The field of each domain's panel: `owned` holds each domain's element boxes, `free` the boxes of what lies
    in no domain (other elements), `turned` the turned boxes of what lies in none either (stems), `pad` the margin
    around elements, `names` each domain's name box (wrapped NAME_PAD, and kept by its domain against the others).
    Domains in the order given."""
    from scipy.ndimage import distance_transform_edt

    names = names or {}
    every = [b for bs in owned.values() for b in bs] + list(free) + list(names.values())
    every += [(x - hl - hw, y - hl - hw, x + hl + hw, y + hl + hw) for x, y, hl, hw, _ in turned]
    if not owned or not every:
        return Fields(np.empty(0), np.empty(0), {})
    reach = pad + CLOSE + GUTTER + 1.0
    r = np.asarray(every, float)
    xs = _axis(r[:, 0].min() - reach, r[:, 2].max() + reach)
    ys = _axis(r[:, 1].min() - reach, r[:, 3].max() + reach)
    elements = {d: _distance(bs, xs, ys, reach) for d, bs in owned.items()}
    none = np.full((len(ys), len(xs)), np.inf)
    named = {d: _distance([names[d]], xs, ys, reach) if d in names else none for d in owned}
    own = {d: np.minimum(elements[d], named[d]) for d in owned}
    loose = np.minimum(_distance(list(free), xs, ys, reach), _turned_distance(list(turned), xs, ys, reach))
    values = {}
    for d in owned:
        rival = loose.copy()
        for e in owned:
            if e != d:
                np.minimum(rival, own[e], out=rival)
        base = np.maximum(pad - elements[d], NAME_PAD - named[d])
        dilated = elements[d] <= pad + CLOSE
        closed = distance_transform_edt(dilated) * GRID - CLOSE  # >= 0: inside the dilation eroded back
        shape = np.maximum(base, closed)
        with np.errstate(invalid="ignore"):  # inf - inf far from everything: outside
            split = np.maximum((rival - own[d] - GUTTER) / 2, np.minimum(EDGE - own[d], (rival - own[d]) / 2))
        split = np.where(np.isinf(rival), np.inf, split)
        values[d] = np.clip(np.nan_to_num(np.minimum(shape, split), nan=-1.0), -1.0, 1.0)
    return Fields(xs, ys, values)
