"""Kabsch-Sander DSSP: H-bond energies, helices, bridges.

Not implemented: bulge merging and the turn/bend (T, S) classes. Output codes:
H alpha, G 3-10, I pi, E strand (ladder of >= 2 bridges), B isolated bridge, '-' other.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.spatial import cKDTree

from .model import Backbone, Bridge

_Q1Q2F = 0.084 * 332.0  # kcal/mol * Å
_HBOND_E = -0.5
_CA_CUTOFF = 9.0


class DsspResult(NamedTuple):
    ss: str
    bridges: list[Bridge]


def _hydrogen_positions(bb: Backbone) -> np.ndarray:
    n = len(bb)
    h = np.full((n, 3), np.nan)
    for i in range(1, n):
        if bb.prev[i] and bb.labels[i].name != "PRO":
            v = bb.xyz[i - 1, 2] - bb.xyz[i - 1, 3]  # C(i-1) - O(i-1)
            h[i] = bb.xyz[i, 0] + v / np.linalg.norm(v)
    return h


def _hbond_set(bb: Backbone) -> set[tuple[int, int]]:
    """Set of (acceptor, donor): CO of acceptor bonded to NH of donor."""
    xyz = bb.xyz
    h = _hydrogen_positions(bb)
    close = np.array(sorted(cKDTree(bb.ca).query_pairs(_CA_CUTOFF)))
    if close.size == 0:
        return set()
    both = np.vstack([close, close[:, ::-1]])
    don, acc = both[:, 0], both[:, 1]
    keep = ~np.isnan(h[don, 0]) & ~((acc == don - 1) & bb.prev[don])
    don, acc = don[keep], acc[keep]
    hd, nd = h[don], xyz[don, 0]
    oa, ca = xyz[acc, 3], xyz[acc, 2]
    d = lambda p, q: np.linalg.norm(p - q, axis=1)  # noqa: E731
    e = _Q1Q2F * (1 / d(oa, nd) + 1 / d(ca, hd) - 1 / d(oa, hd) - 1 / d(ca, nd))
    order = np.argsort(e)
    used_d, used_a, bonds = {}, {}, set()
    for k in order:
        if e[k] >= _HBOND_E:
            break
        di, ai = int(don[k]), int(acc[k])
        if used_d.get(di, 0) < 2 and used_a.get(ai, 0) < 2:
            used_d[di] = used_d.get(di, 0) + 1
            used_a[ai] = used_a.get(ai, 0) + 1
            bonds.add((ai, di))
    return bonds


def _resolve(n, alpha, g, pi, ladder, bridge) -> list[str]:
    """Per-residue code from overlapping evidence. DSSP priority: H > E > B > G > I."""
    out = []
    for i in range(n):
        if i in alpha:
            out.append("H")
        elif i in ladder:
            out.append("E")
        elif i in bridge:
            out.append("B")
        elif i in g:
            out.append("G")
        elif i in pi:
            out.append("I")
        else:
            out.append("-")
    return out


def assign_dssp(bb: Backbone) -> DsspResult:
    n, prev, nxt = len(bb), bb.prev, bb.nxt
    bonds = _hbond_set(bb)

    def hb(a: int | None, b: int | None) -> bool:
        return a is not None and b is not None and (a, b) in bonds

    def before(i: int) -> int | None:
        return i - 1 if prev[i] else None

    def after(i: int) -> int | None:
        return i + 1 if nxt[i] else None

    turns: dict[str, set[int]] = {"H": set(), "G": set(), "I": set()}
    for k, code in ((4, "H"), (3, "G"), (5, "I")):
        turn = [i + k < n and (i, i + k) in bonds and all(nxt[i + m] for m in range(k)) for i in range(n)]
        for i in range(1, n - k):
            if turn[i - 1] and turn[i]:
                turns[code].update(range(i, i + k))

    bridges: list[Bridge] = []
    close = cKDTree(bb.ca).query_pairs(_CA_CUTOFF)
    for i, j in sorted(close):
        if j < i + 3:
            continue
        im, ip, jm, jp = before(i), after(i), before(j), after(j)
        if (hb(im, j) and hb(j, ip)) or (hb(jm, i) and hb(i, jp)):
            bridges.append(Bridge(i, j, "P"))
        if (hb(i, j) and hb(j, i)) or (hb(im, jp) and hb(jm, ip)):
            bridges.append(Bridge(i, j, "A"))

    known = set(bridges)
    in_ladder, in_bridge = set(), set()
    for i, j, kind in bridges:
        in_bridge.update((i, j))
        if i + 1 < n and prev[i + 1]:
            nj = j + 1 if kind == "P" else j - 1
            if 0 <= nj < n and Bridge(i + 1, nj, kind) in known:
                in_ladder.update((i, j, i + 1, nj))
    ss = _resolve(n, turns["H"], turns["G"], turns["I"], in_ladder, in_bridge)
    return DsspResult("".join(ss), bridges)
