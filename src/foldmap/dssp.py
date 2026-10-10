"""Kabsch-Sander DSSP: H-bond energies, n-turns, helices, bridges and ladders.

Not implemented: beta-bulges and bends (S). Output codes: H alpha, G 3-10, I pi, E strand (ladder of >= 2
bridges), B isolated bridge, T H-bonded turn, '-' other. DSSP priority: H > E > B > G > I.
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


class Ladder(NamedTuple):
    """Bridges of one kind along two strands; i and j ascending (j runs backwards for antiparallel)."""

    kind: str
    i: list[int]
    j: list[int]


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
    used_d, bonds = {}, set()
    for k in order:
        if e[k] >= _HBOND_E:
            break
        di, ai = int(don[k]), int(acc[k])
        if used_d.get(di, 0) < 2:  # DSSP keeps each NH's two best acceptors; a CO may accept any number
            used_d[di] = used_d.get(di, 0) + 1
            bonds.add((ai, di))
    return bonds


def _turn_starts(bonds: set[tuple[int, int]], nxt: np.ndarray, n: int) -> dict[int, list[bool]]:
    """starts[k][i]: CO(i) is bonded to NH(i+k) with no chain break in between (an n-turn of length k)."""
    return {
        k: [i + k < n and (i, i + k) in bonds and all(nxt[i + m] for m in range(k)) for i in range(n)]
        for k in (3, 4, 5)
    }


def _find_bridges(bb: Backbone, bonds: set[tuple[int, int]]) -> list[Bridge]:
    """Bridges i < j (j >= i + 3). As in DSSP, i-1, i+1, j-1 and j+1 must exist and be peptide-bonded to i and j."""
    prev, nxt = bb.prev, bb.nxt
    out: list[Bridge] = []
    for i, j in sorted(cKDTree(bb.ca).query_pairs(_CA_CUTOFF)):
        if j < i + 3 or not (prev[i] and nxt[i] and prev[j] and nxt[j]):
            continue
        if ((i - 1, j) in bonds and (j, i + 1) in bonds) or ((j - 1, i) in bonds and (i, j + 1) in bonds):
            out.append(Bridge(i, j, "P"))
        if ((i, j) in bonds and (j, i) in bonds) or ((i - 1, j + 1) in bonds and (j - 1, i + 1) in bonds):
            out.append(Bridge(i, j, "A"))
    return out


def _ladders(bridges: list[Bridge], prev: np.ndarray) -> list[Ladder]:
    """Consecutive bridges of one kind form a ladder."""
    ladders: list[Ladder] = []
    for i, j, kind in sorted(bridges):
        for lad in ladders:
            if kind == "P" and lad.kind == "P" and i == lad.i[-1] + 1 and prev[i] and j == lad.j[-1] + 1:
                lad.i.append(i)
                lad.j.append(j)
                break
            if kind == "A" and lad.kind == "A" and i == lad.i[-1] + 1 and prev[i] and j == lad.j[0] - 1:
                lad.i.append(i)
                lad.j.insert(0, j)
                break
        else:
            ladders.append(Ladder(kind, [i], [j]))
    return ladders


def _assign_codes(n: int, ladders: list[Ladder], starts: dict[int, list[bool]]) -> list[str]:
    """Per-residue codes: E/B from the ladders, then H over them, then G and I on residues still '-',
    then T inside any turn that is still '-'."""
    ss = ["-"] * n
    for lad in ladders:
        code = "E" if len(lad.i) > 1 else "B"
        for lo, hi in ((lad.i[0], lad.i[-1]), (lad.j[0], lad.j[-1])):
            for r in range(lo, hi + 1):
                if ss[r] != "E":
                    ss[r] = code
    for k, code in ((4, "H"), (3, "G"), (5, "I")):
        turn = starts[k]
        for i in range(1, n - k):
            if turn[i - 1] and turn[i]:
                for r in range(i, i + k):
                    if code == "H" or ss[r] == "-":
                        ss[r] = code
    inside = {i + m for k, turn in starts.items() for i in range(n) if turn[i] for m in range(1, k)}
    return ["T" if c == "-" and i in inside else c for i, c in enumerate(ss)]


def assign_dssp(bb: Backbone) -> DsspResult:
    n = len(bb)
    bonds = _hbond_set(bb)
    bridges = _find_bridges(bb, bonds)
    ss = _assign_codes(n, _ladders(bridges, bb.prev), _turn_starts(bonds, bb.nxt, n))
    return DsspResult("".join(ss), bridges)
