"""Kabsch-Sander DSSP as mkdssp 3.x does it: H-bond energies, n-turns, helices, bridges, ladders, beta-bulges.

Not implemented: bends (S). Output codes: H alpha, G 3-10, I pi, E strand (ladder of >= 2 bridges, or ladders
joined across a beta-bulge), B isolated bridge, T H-bonded turn, '-' other. Priority as in mkdssp >= 2.1:
H replaces E/B; G and I need every residue of the segment free, except that I may replace H (prefer_pi).
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.spatial import cKDTree

from .model import Backbone, Bridge

_Q1Q2F = 0.084 * 332.0  # kcal/mol * Å
_HBOND_E = -0.5
_HBOND_E_MIN = -9.9  # kcal/mol: DSSP's floor, also the energy of any pair with atoms closer than _MIN_DIST
_MIN_DIST = 0.5  # Å
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


def _hbond_energy(d_on, d_ch, d_oh, d_cn) -> np.ndarray:
    """Electrostatic energy from the O-N, C-H, O-H and C-N distances, rounded to 0.001 kcal/mol and floored at
    -9.9 kcal/mol as in mkdssp."""
    d_on, d_ch, d_oh, d_cn = (np.asarray(x, float) for x in (d_on, d_ch, d_oh, d_cn))
    with np.errstate(divide="ignore", invalid="ignore"):
        e = _Q1Q2F * (1 / d_on + 1 / d_ch - 1 / d_oh - 1 / d_cn)
    e = np.round(e * 1000) / 1000
    clash = np.minimum.reduce([d_on, d_ch, d_oh, d_cn]) < _MIN_DIST
    return np.where(clash, _HBOND_E_MIN, np.maximum(e, _HBOND_E_MIN))


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
    e = _hbond_energy(d(oa, nd), d(ca, hd), d(oa, hd), d(ca, nd))
    order = np.lexsort((acc, e))  # equal energies: the lower acceptor wins, as in mkdssp
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
    """Bridges i < j (j >= i + 3). As in DSSP, i-1, i+1, j-1 and j+1 must exist and be peptide-bonded to i and j,
    and a pair matching both patterns is parallel."""
    prev, nxt = bb.prev, bb.nxt
    out: list[Bridge] = []
    for i, j in sorted(cKDTree(bb.ca).query_pairs(_CA_CUTOFF)):
        if j < i + 3 or not (prev[i] and nxt[i] and prev[j] and nxt[j]):
            continue
        if ((i - 1, j) in bonds and (j, i + 1) in bonds) or ((j - 1, i) in bonds and (i, j + 1) in bonds):
            out.append(Bridge(i, j, "P"))
        elif ((i, j) in bonds and (j, i) in bonds) or ((i - 1, j + 1) in bonds and (j - 1, i + 1) in bonds):
            out.append(Bridge(i, j, "A"))
    return out


def _unbroken(prev: np.ndarray, lo: int, hi: int) -> bool:
    return bool(prev[lo + 1 : hi + 1].all())


def _gap(a: int, b: int) -> int:
    """a - b as mkdssp computes it, in unsigned arithmetic: a negative gap is never small."""
    return a - b if a >= b else 1 << 30


def _ladders(bridges: list[Bridge], prev: np.ndarray) -> list[Ladder]:
    """Consecutive bridges of one kind form a ladder. Two ladders of one kind are joined across a beta-bulge
    when the gap is < 6 residues on one strand and < 3 on the other (mkdssp's rule)."""
    ladders: list[Ladder] = []
    last: dict[tuple[str, int, int], Ladder] = {}  # (kind, i, j) of each ladder's latest bridge
    for i, j, kind in sorted(bridges):
        lad = last.pop((kind, i - 1, j - 1 if kind == "P" else j + 1), None) if prev[i] else None
        if lad is None:
            lad = Ladder(kind, [i], [j])
            ladders.append(lad)
        else:
            lad.i.append(i)
            if kind == "P":
                lad.j.append(j)
            else:
                lad.j.insert(0, j)
        last[(kind, i, j)] = lad
    ladders.sort(key=lambda lad: (lad.i[0], lad.j[0]))
    a = 0
    while a < len(ladders):
        b = a + 1
        while b < len(ladders):
            x, y = ladders[a], ladders[b]
            ibi, iei, jbi, jei = x.i[0], x.i[-1], x.j[0], x.j[-1]
            ibj, iej, jbj, jej = y.i[0], y.i[-1], y.j[0], y.j[-1]
            if ibj >= iei + 6:  # sorted by first i: no later ladder can join this one either
                break
            joinable = (
                x.kind == y.kind
                and _gap(ibj, iei) < 6
                and not (iei >= ibj and ibi <= iej)
                and _unbroken(prev, min(ibi, ibj), max(iei, iej))
                and _unbroken(prev, min(jbi, jbj), max(jei, jej))
            )
            other = _gap(jbj, jei) if x.kind == "P" else _gap(jbi, jej)
            if joinable and ((other < 6 and _gap(ibj, iei) < 3) or other < 3):
                ladders[a] = Ladder(x.kind, sorted(set(x.i) | set(y.i)), sorted(set(x.j) | set(y.j)))
                del ladders[b]
            else:
                b += 1
        a += 1
    return ladders


def _assign_codes(n: int, ladders: list[Ladder], starts: dict[int, list[bool]], prefer_pi: bool = True) -> list[str]:
    """Per-residue codes in mkdssp's order: E/B from the ladders; H over them; then G and I only where every
    residue of the segment is free (I may replace H when prefer_pi); T inside any turn that is still '-'."""
    ss = ["-"] * n
    for lad in ladders:
        code = "E" if len(lad.i) > 1 else "B"
        for lo, hi in ((lad.i[0], lad.i[-1]), (lad.j[0], lad.j[-1])):
            for r in range(lo, hi + 1):
                if ss[r] != "E":
                    ss[r] = code
    for k, code, free in ((4, "H", None), (3, "G", "-G"), (5, "I", "-IH" if prefer_pi else "-I")):
        turn = starts[k]
        for i in range(1, n - k):
            if turn[i - 1] and turn[i]:
                seg = range(i, i + k)
                if free is None or all(ss[r] in free for r in seg):
                    for r in seg:
                        ss[r] = code
    inside = {i + m for k, turn in starts.items() for i in range(n) if turn[i] for m in range(1, k)}
    return ["T" if c == "-" and i in inside else c for i, c in enumerate(ss)]


def assign_dssp(bb: Backbone, prefer_pi: bool = True) -> DsspResult:
    n = len(bb)
    bonds = _hbond_set(bb)
    bridges = _find_bridges(bb, bonds)
    ss = _assign_codes(n, _ladders(bridges, bb.prev), _turn_starts(bonds, bb.nxt, n), prefer_pi)
    return DsspResult("".join(ss), bridges)
