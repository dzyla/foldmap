# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Residue alignment of one chain of a structure A (the reference) with one chain of a structure B (plan 13a).

Three routes give the residue pairs: by sequence (the global alignment A19 uses, figcheck.aligned_positions), by
structure (a TM-align-style aligner on Cα, Zhang & Skolnick 2005) or from a file (a two-row alignment, the escape
hatch for any pair the aligner gets wrong). `auto` takes the sequence route when the sequence alignment is at least
40 % identical and its pairs superpose as one fold (TM-score at least SEQUENCE_ROUTE_TM), else the structure route.
Every route then scores its pairs with the TM-score normalised by A's chain length, over the superposition of B on A
that maximises it (TM-align's search: Kabsch on fragments of the pairs, then on the pairs within a cutoff, repeated),
always a proper rotation (symmetry.kabsch, never a mirror).
element_mates then says which of B's elements are which of A's (mates, pieces of a split helix, extras) and whether
the two chains are one fold (H1).

Determinism (A12): distances are snapped to 1/GRID Å and pair scores to 1/_SCORE before any choice, so every sum
the dynamic programming and the TM-score compare is exact and independent of summation order; ties go to the
earlier candidate (seed, round, fragment, threading shift, local superposition, traceback move) in a fixed order.
The work is bounded by counts (ROUNDS rounds, a fixed set of fragments, shifts and local superpositions, the last
capped by LOCAL_CELLS), never by a clock.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .model import SSE, Backbone, ResLabel, overlap, pair_places

GRID = 1024  # distances are snapped to 1/GRID Å before scoring: exact in binary floating point, the same on any CPU
_SCORE = 2.0**20  # pair scores are snapped to 1/_SCORE: sums of up to 2^33 of them are exact in a double
GAP_OPEN = -0.6  # TM-align's gap-open penalty (Zhang & Skolnick 2005); no extension penalty
SS_GAP = -1.0  # the gap penalty of TM-align's secondary-structure seed (a match scores 1)
ROUNDS = 20  # superpose / score / re-align rounds per seed, at most; also the iterations of one fragment's search
ROUND_STEP = 40  # fragment starts every ROUND_STEP pairs in the rounds' search (TM-align's simplify_step); 1 at the end
SEQUENCE_IDENTITY = 0.4  # auto: the sequence route at this identity over the aligned residues, else the structure
# auto: ... and only when the sequence alignment's pairs superpose at this TM-score (normalised by A): a gapped
# alignment of two unrelated sequences can reach the identity (13a-10). Its own constant, the value H1's SAME_FOLD_TM
# had when it was set, so a retune of H1 does not move the routing
SEQUENCE_ROUTE_TM = 0.5
SAME_FOLD_TM = 0.5  # H1: one fold at this TM-score (normalised by A) ... (frozen on the tuning set, 13a-10)
SAME_FOLD_FRACTION = 0.5  # ... with this fraction of the chain with fewer elements matched (frozen, 13a-10)
# TM-align's local-superposition seed (get_initial5, plan 24): fragments of these lengths (at most a third, then a
# half, of the shorter chain) superposed pairwise, their starts every _jump(length) residues along each chain
LOCAL_FRAGMENTS = (20, 100)
# ... the work bounded by a count: the starts' alignments cover at most LOCAL_CELLS residue pairs in all (both jumps
# widened together, by whole multiples, until they do: TM-align's own steps up to about 520 residues a chain) ...
LOCAL_CELLS = 2**26
LOCAL_BATCH = 2**22  # ... aligned LOCAL_BATCH residue pairs at a time (memory; the result does not depend on it)
WIDER = 1.5  # get_initial5's and initial_ssplus's distance scale: d0 + WIDER Å
SS_BONUS = 0.5  # initial_ssplus: a residue pair's score gains this when the two residues' H/E/C agree
# a candidate alignment that is not one of the first three seeds' penalised refinements (those are always scored, so
# the result is never below what they alone give) is scored in full (tm_score: a fragment start at every pair) only
# when its quick score (_search at ROUND_STEP) exceeds the best alignment's quick score less SCREEN·L
SCREEN = 0.0  # 0, 0.003 and 0.01 chose the same alignments on the 55 homologue pairs (plan 24)


@dataclass(frozen=True)
class ChainAlignment:
    pairs: tuple[tuple[int, int], ...]  # (index into A's residues of chain_a, index into B's of chain_b), increasing
    tm_score: float  # normalised by A's chain length (the reference)
    rmsd: float  # over `pairs`, after their own least-squares superposition (as TM-align reports it)
    identity: float  # identical residues / len(pairs)
    method: str  # "sequence" | "structure" | "file"
    rotation: np.ndarray = field(compare=False)  # 3x3 proper rotation, A ≈ R·B + t (the TM-score's superposition)
    translation: np.ndarray = field(compare=False)  # (3,)
    # the residues `pairs` index (chain_a's of A, chain_b's of B): how element_mates finds a structure's elements,
    # whose start/end count all its residues; empty: the elements' start/end are read as these indices
    residues_a: tuple[ResLabel, ...] = field(default=(), repr=False)
    residues_b: tuple[ResLabel, ...] = field(default=(), repr=False)


def d0(length: int) -> float:
    """TM-score's distance scale d0(L) = 1.24·(L−15)^(1/3) − 1.8, at least 0.5 (Zhang & Skolnick 2004)."""
    if length <= 15:  # the cube root of a negative number is negative: below 0.5 anyway
        return 0.5
    return max(0.5, 1.24 * (length - 15) ** (1 / 3) - 1.8)


def _terms(d: np.ndarray, scale: float) -> np.ndarray:
    """Each distance's TM-score term 1/(1+(d/d0)²), snapped to 1/_SCORE."""
    return np.round(_SCORE / (1.0 + (d / scale) ** 2)) / _SCORE


def pair_distances(pa: np.ndarray, pb: np.ndarray, r: np.ndarray, t: np.ndarray) -> np.ndarray:
    """|R·b + t − a| for each pair, snapped to 1/GRID Å."""
    return np.round(np.sqrt(((pb @ r.T + t - pa) ** 2).sum(axis=1)) * GRID) / GRID


def _cutoff(scale: float) -> float:
    """TM-align's d0_search: the distance within which pairs are refitted, d0 kept within 4.5–8 Å."""
    return min(max(scale, 4.5), 8.0)


def _within(d: np.ndarray, cut: float) -> np.ndarray:
    """The pairs closer than cut, the cut widened by 0.5 Å until at least three are (or all, when fewer)."""
    sel = d < cut
    while sel.sum() < min(3, len(d)):
        cut += 0.5
        sel = d < cut
    return sel


def _fit(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """symmetry.kabsch's R, t (q ≈ p·Rᵀ + t, a proper rotation), the same numbers (a mean is the sum divided by the
    count, as numpy's), without its RMSD: the searches fit tens of thousands of selections."""
    pc, qc = p.sum(axis=0) / len(p), q.sum(axis=0) / len(q)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    r = vt.T @ u.T
    (a, b, c), (d, e, f), (g, h, i) = r.tolist()
    if a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g) < 0:  # kabsch's d: a mirror made proper
        r = (vt.T * np.array([1.0, 1.0, -1.0])) @ u.T
    return r, qc - pc @ r.T


def _search(pa: np.ndarray, pb: np.ndarray, scale: float, step: int) -> tuple[float, np.ndarray, np.ndarray]:
    """(sum of TM-score terms, R, t): the superposition of the paired points pb on pa (pa ≈ R·pb + t) with the
    highest TM-score that TM-align's search finds. Starts: all pairs, then contiguous fragments of n/2, n/4 … n/32
    pairs (at least 4) every `step` pairs; from each, Kabsch on the selection, keep the pairs closer than
    d0_search, refit, up to ROUNDS times or until the selection repeats one already fitted."""
    n = len(pa)
    starts, size = [(0, n)], n // 2
    for _ in range(5):
        if size < 4:
            break
        firsts = list(range(0, n - size + 1, step))
        if firsts[-1] != n - size:
            firsts.append(n - size)
        starts += [(k, size) for k in firsts]
        size //= 2
    cut = _cutoff(scale)
    best, best_r, best_t = -1.0, np.eye(3), np.zeros(3)
    fitted: set[bytes] = set()  # a selection already fitted leads where it led before: nothing new to find
    for k, size in starts:
        sel = np.zeros(n, bool)
        sel[k : k + size] = True
        for _ in range(ROUNDS):
            key = np.packbits(sel).tobytes()
            if key in fitted:
                break
            fitted.add(key)
            r, t = _fit(pb[sel], pa[sel])
            d = pair_distances(pa, pb, r, t)
            score = float(_terms(d, scale).sum())  # exact: snapped terms
            if score > best:
                best, best_r, best_t = score, r, t
            sel = _within(d, cut)
    return best, best_r, best_t


def tm_score(ca_a, ca_b, pairs, length_a, step: int = 1) -> tuple[float, np.ndarray, np.ndarray, float]:
    """(TM-score, R, t, RMSD) of B's Cα on A's over `pairs` ((i in A, j in B)), with A ≈ R·B + t the proper rotation
    that maximises the TM-score (_search; fragment starts every `step` pairs). Each pair scores 1/(1+(d/d0)²),
    d snapped to 1/GRID Å, summed over length_a. The RMSD is over all pairs after their own Kabsch fit."""
    from .symmetry import kabsch

    if not len(pairs):
        return 0.0, np.eye(3), np.zeros(3), 0.0
    ia, ib = (np.array(k, dtype=int) for k in zip(*pairs))
    pa, pb = np.asarray(ca_a, float)[ia], np.asarray(ca_b, float)[ib]
    score, r, t = _search(pa, pb, d0(length_a), step)
    rmsd = kabsch(pb, pa)[2]
    return score / length_a, r, t, rmsd


def _dp(score: np.ndarray, gap: float) -> list[tuple[int, int]]:
    """The pairs (i, j) of the global alignment of rows i with columns j that maximises the summed score of the
    pairs plus `gap` for each gap opened (no extension penalty; gaps at either end are free). Three states, M (i
    and j paired), X (row i left out: a gap in B) and Y (column j left out: a gap in A), filled a row at a time:
    with no extension penalty Y along a row is a running maximum. Ties: a pair, then a gap in B, then a gap in A,
    both when filling and in the traceback. Exact when the scores and gap are multiples of 1/_SCORE."""
    n, m = score.shape
    if not n or not m:
        return []
    low = -np.inf
    mm, xx, yy = (np.full((n + 1, m + 1), low) for _ in range(3))
    mm[0, :] = mm[:, 0] = 0.0  # the start: anything before the first pair is free
    pm, px, py = (np.zeros((n + 1, m + 1), np.int8) for _ in range(3))
    gaps = np.full(m, gap)
    gaps[0] = 0.0  # a gap opened from column 0 is a leading gap

    def pick(a, b, c):
        return np.where((a >= b) & (a >= c), 0, np.where(b >= c, 1, 2)).astype(np.int8)

    for i in range(1, n + 1):
        m0, x0, y0 = mm[i - 1], xx[i - 1], yy[i - 1]
        a, b, c = m0[:-1], x0[:-1], y0[:-1]
        pm[i, 1:] = pick(a, b, c)
        mm[i, 1:] = score[i - 1] + np.maximum(np.maximum(a, b), c)
        g = 0.0 if i == 1 else gap  # from row 0: a leading gap
        a, b, c = m0[1:] + g, x0[1:], y0[1:] + g
        px[i, 1:] = pick(a, b, c)
        xx[i, 1:] = np.maximum(np.maximum(a, b), c)
        a, b = mm[i, :-1] + gaps, xx[i, :-1] + gaps
        run = np.maximum.accumulate(np.maximum(a, b))  # Y[i, j] for j = 1..m
        before = np.concatenate(([low], run[:-1]))  # Y[i, j-1]
        py[i, 1:] = pick(a, b, before)
        yy[i, 1:] = run
    # the end: any cell of the last row or column (what follows it is a free end gap), the corner first
    ends = [(n, j) for j in range(m, -1, -1)] + [(i, m) for i in range(n - 1, -1, -1)]
    best = np.maximum(np.maximum(mm, xx), yy)
    i, j = ends[int(np.argmax([best[e] for e in ends]))]
    state = int(pick(mm[i, j], xx[i, j], yy[i, j]))
    out = []
    while i > 0 and j > 0:
        if state == 0:
            out.append((i - 1, j - 1))
            state = int(pm[i, j])
            i, j = i - 1, j - 1
        elif state == 1:
            state = int(px[i, j])
            i -= 1
        else:
            state = int(py[i, j])
            j -= 1
    return out[::-1]


def _matrix(ca_a: np.ndarray, ca_b: np.ndarray, r: np.ndarray, t: np.ndarray, scale: float) -> np.ndarray:
    """Every residue pair's TM-score term with B superposed by (r, t), on snapped distances: _terms of
    pair_distances' distances, the same numbers, computed in place a coordinate at a time."""
    moved = ca_b @ r.T + t
    x = np.zeros((len(ca_a), len(moved)))
    for k in range(3):  # x² + y² + z², summed in that order
        diff = np.subtract.outer(ca_a[:, k], moved[:, k])
        np.multiply(diff, diff, out=diff)
        x += diff
    np.sqrt(x, out=x)
    x *= GRID
    np.round(x, out=x)
    x /= GRID
    x /= scale
    np.multiply(x, x, out=x)
    x += 1.0
    np.divide(_SCORE, x, out=x)
    np.round(x, out=x)
    x /= _SCORE
    return x


def _quick(pa: np.ndarray, pb: np.ndarray, scale: float) -> float:
    """A gapless threading's score, as TM-align's get_score_fast: the best TM-score sum of a Kabsch fit on all
    pairs and two refits on the pairs within d0_search, then within √(d0_search² + 1)."""
    r, t = _fit(pb, pa)
    d = pair_distances(pa, pb, r, t)
    best = float(_terms(d, scale).sum())
    for cut in (_cutoff(scale), float(np.sqrt(_cutoff(scale) ** 2 + 1))):
        sel = _within(d, cut)
        if sel.all():
            continue
        r, t = _fit(pb[sel], pa[sel])
        d = pair_distances(pa, pb, r, t)
        best = max(best, float(_terms(d, scale).sum()))
    return best


def _thread(ca_a: np.ndarray, ca_b: np.ndarray, scale: float) -> list[tuple[int, int]]:
    """The best gapless threading of one chain along the other: B's residue j against A's j + k for each shift k
    that overlaps at least half the shorter chain (and at least 5 residues), scored by _quick; ties to the
    smaller k."""
    n, m = len(ca_a), len(ca_b)
    least = max(min(n, m) // 2, 5)
    best, out = -1.0, []
    for k in range(least - m, n - least + 1):
        j0, j1 = max(0, -k), min(m, n - k)
        if j1 - j0 < 3:
            continue
        score = _quick(ca_a[j0 + k : j1 + k], ca_b[j0:j1], scale)
        if score > best:
            best, out = score, [(j + k, j) for j in range(j0, j1)]
    return out


def _free_dp(score: np.ndarray) -> list[list[tuple[int, int]]]:
    """For each matrix of score (k, rows, columns), the pairs of the alignment that maximises their summed score
    with no gap penalty at all (TM-align's NWDP_TM with gap_open 0): H[i, j] = max(H[i-1, j-1] + score, H[i-1, j],
    H[i, j-1]), all k filled a row at a time (along a row, a running maximum: every score is positive). Traceback
    from the corner: a pair, then a gap in A (j-1), then a gap in B (i-1). Exact when the scores are multiples of
    1/_SCORE."""
    k, n, m = score.shape
    h = np.zeros((k, n + 1, m + 1))
    for i in range(1, n + 1):
        row = np.maximum(h[:, i - 1, :-1] + score[:, i - 1], h[:, i - 1, 1:])
        np.maximum.accumulate(row, axis=1, out=h[:, i, 1:])
    out = []
    for hb, sb in zip(h, score):
        i, j, pairs = n, m, []
        while i and j:
            here = hb.item(i, j)
            if here == hb.item(i - 1, j - 1) + sb.item(i - 1, j - 1):
                i, j = i - 1, j - 1
                pairs.append((i, j))
            elif here == hb.item(i, j - 1):
                j -= 1
            else:
                i -= 1
        out.append(pairs[::-1])
    return out


def _jump(length: int) -> int:
    """get_initial5's step between fragment starts along a chain of `length` residues, at most a third of it."""
    step = 45 if length > 250 else 35 if length > 200 else 25 if length > 150 else 15
    return max(1, min(step, length // 3))


def _local_starts(n: int, m: int) -> list[tuple[int, int, int]]:
    """The fragment pairs _local superposes, as (A's start, B's start, length), in order: LOCAL_FRAGMENTS (at most a
    third, then a half, of the shorter chain; none under 3), the starts every _jump(length) residues along each
    chain, both steps widened by the least whole multiple that keeps the work (the starts × n × m residue pairs
    aligned) within LOCAL_CELLS; none when a single start per length is already over it."""
    shorter = min(n, m)
    sizes = [s for s in (min(LOCAL_FRAGMENTS[0], shorter // 3), min(LOCAL_FRAGMENTS[1], shorter // 2)) if s >= 3]
    wider = 1
    while True:
        step_a, step_b = _jump(n) * wider, _jump(m) * wider
        out = [(i, j, s) for s in sizes for i in range(0, n - s + 1, step_a) for j in range(0, m - s + 1, step_b)]
        if len(out) * n * m <= LOCAL_CELLS:
            return out
        if step_a >= n and step_b >= m:
            return []
        wider += 1


def _local(ca_a: np.ndarray, ca_b: np.ndarray, scale: float) -> list[tuple[int, int]]:
    """TM-align's local-superposition seed (get_initial5): for each fragment pair of _local_starts, B's fragment
    superposed on A's (Kabsch), the chains aligned under that superposition by _free_dp on the TM-score terms with
    d0 + WIDER, and the alignment scored as a gapless threading is (_quick). The alignment that scores highest; ties
    to the earlier start."""
    n, m = len(ca_a), len(ca_b)
    todo = _local_starts(n, m)
    wide = scale + WIDER
    batch = max(1, LOCAL_BATCH // (n * m))
    best, out = -1.0, []
    for at in range(0, len(todo), batch):
        part = todo[at : at + batch]
        mats = np.empty((len(part), n, m))
        for k, (i, j, s) in enumerate(part):
            r, t = _fit(ca_b[j : j + s], ca_a[i : i + s])
            mats[k] = _matrix(ca_a, ca_b, r, t, wide)
        for pairs in _free_dp(mats):
            if len(pairs) < 3:
                continue
            ia, ib = (np.array(x, dtype=int) for x in zip(*pairs))
            score = _quick(ca_a[ia], ca_b[ib], scale)
            if score > best:
                best, out = score, pairs
    return out


def _ss_plus(ca_a: np.ndarray, ca_b: np.ndarray, ss_a: np.ndarray, ss_b: np.ndarray, pairs, scale: float):
    """TM-align's secondary-structure-plus-distance seed (initial_ssplus): B superposed on A by Kabsch on `pairs`
    (the best alignment found so far), each residue pair scored by its TM-score term with d0 + WIDER plus SS_BONUS
    when the two residues' H/E/C agree, aligned by _dp with SS_GAP."""
    if len(pairs) < 3:
        return []
    ia, ib = (np.array(x, dtype=int) for x in zip(*pairs))
    r, t = _fit(ca_b[ib], ca_a[ia])
    score = _matrix(ca_a, ca_b, r, t, scale + WIDER) + SS_BONUS * (ss_a[:, None] == ss_b[None, :])
    return _dp(score, SS_GAP)


def _three_state(bb: Backbone) -> np.ndarray:
    """The chain's secondary structure from foldmap's DSSP as H (H, G, I), E (E, B) or C (the rest)."""
    from .dssp import assign_dssp

    ss = assign_dssp(bb).ss
    return np.array(["H" if c in "HGI" else "E" if c in "EB" else "C" for c in ss])


def _refine(ca_a: np.ndarray, ca_b: np.ndarray, pairs, scale: float, free: bool = False):
    """From a seed alignment, up to ROUNDS rounds: superpose on the pairs (_search), score every residue pair under
    that superposition, re-align by _dp with GAP_OPEN (free: by _free_dp, no gap penalty, TM-align's second gap
    value in DP_iter); stop when the pairs repeat (unchanged, or a cycle). Returns the round with the highest
    TM-score sum, ties to the earlier: (sum, pairs, R, t)."""
    gap = round(GAP_OPEN * _SCORE) / _SCORE
    best = (-1.0, [], np.eye(3), np.zeros(3))
    seen: set[tuple] = set()
    for _ in range(ROUNDS):
        key = tuple(pairs)
        if not pairs or key in seen:
            break
        seen.add(key)
        ia, ib = (np.array(k, dtype=int) for k in zip(*pairs))
        score, r, t = _search(ca_a[ia], ca_b[ib], scale, ROUND_STEP)
        if score > best[0]:
            best = (score, list(pairs), r, t)
        terms = _matrix(ca_a, ca_b, r, t, scale)
        pairs = _free_dp(terms[None])[0] if free else _dp(terms, gap)
    return best


def _structure_pairs(a: Backbone, b: Backbone, seq_pairs) -> tuple[list[tuple[int, int]], tuple]:
    """TM-align's heuristics on Cα: five seeds, in this order: the sequence alignment, the alignment of the H/E/C
    strings, the best gapless threading (_thread), the best local superposition (_local) and the H/E/C plus
    distance alignment under the best alignment so far (_ss_plus); the last two are plan 24's. Each is refined
    (_refine with GAP_OPEN), and the refinement goes on without a gap penalty from its best round's superposition
    (_refine, free: TM-align's DP_iter tries both gap values); each of the two gives a candidate. A candidate keeps
    its pairs within TM-align's score_d8 = 1.5·L^0.3 + 3.5 Å (L: A's length) and is scored by tm_score (the new ones
    only past the SCREEN); the highest wins, ties to the earlier candidate, so the first three seeds' penalised
    refinements, scored as before, set a floor the others can only raise. d0 is A's throughout: the alignment
    maximises the TM-score it reports. Returns (pairs, tm_score's result for them)."""
    ca_a, ca_b = np.asarray(a.ca, float), np.asarray(b.ca, float)
    n = len(ca_a)
    scale = d0(n)
    ss_a, ss_b = _three_state(a), _three_state(b)
    far = 1.5 * n**0.3 + 3.5
    out, scored = [], tm_score(ca_a, ca_b, [], n)
    known: dict[tuple, tuple] = {}  # tm_score of pairs already scored: seeds often end at the same alignment
    quick: dict[tuple, float] = {}  # _search's sum at ROUND_STEP of pairs: the screen of the new candidates

    def rough(pairs) -> float:
        key = tuple(pairs)
        if key not in quick:
            ia, ib = (np.array(k, dtype=int) for k in zip(*pairs))
            quick[key] = _search(ca_a[ia], ca_b[ib], scale, ROUND_STEP)[0]
        return quick[key]

    def consider(refined, screened: bool) -> None:
        nonlocal out, scored
        _, pairs, r, t = refined
        if not pairs:
            return
        ia, ib = (np.array(k, dtype=int) for k in zip(*pairs))
        near = pair_distances(ca_a[ia], ca_b[ib], r, t) <= far
        pairs = [p for p, keep in zip(pairs, near) if keep]
        key = tuple(pairs)
        if key not in known:
            if screened and out and rough(pairs) <= rough(out) - SCREEN * n:
                return
            known[key] = tm_score(ca_a, ca_b, pairs, n)
        if known[key][0] > scored[0]:
            out, scored = pairs, known[key]

    seeds = [
        lambda: list(seq_pairs),
        lambda: _dp((ss_a[:, None] == ss_b[None, :]).astype(float), SS_GAP),
        lambda: _thread(ca_a, ca_b, scale),
        lambda: _local(ca_a, ca_b, scale),
        lambda: _ss_plus(ca_a, ca_b, ss_a, ss_b, out, scale),  # `out` as the seeds before it left it
    ]
    for k, seed in enumerate(seeds):
        refined = _refine(ca_a, ca_b, seed(), scale)
        consider(refined, screened=k >= 3)
        if refined[1]:  # TM-align's DP_iter goes on without a gap penalty from where the penalised rounds ended
            _, _, r, t = refined
            free = _refine(ca_a, ca_b, _free_dp(_matrix(ca_a, ca_b, r, t, scale)[None])[0], scale, free=True)
            consider(free, screened=True)
    return out, scored


def _chain(bb: Backbone, chain: str, side: str) -> list[int]:
    mine = [k for k, l in enumerate(bb.labels) if l.chain == chain]
    if not mine:
        found = ", ".join(dict.fromkeys(l.chain for l in bb.labels)) or "none"
        raise ValueError(f"no chain {chain!r} in {side} (chains: {found})")
    return mine


def _file_pairs(path: str, seq_a: str, seq_b: str, chain_a: str, chain_b: str) -> list[tuple[int, int]]:
    """The pairs a two-row alignment file gives: row 1 is A's chain, row 2 B's; each row without its gaps must be
    its chain's sequence (the modelled residues, one-letter codes)."""
    from .sequence import GAP, read_alignment

    if not Path(path).is_file():
        raise ValueError(f"match {path!r}: not auto, sequence, structure or an alignment file")
    aln = read_alignment(path)
    if len(aln.rows) != 2:
        extra = ", ".join(f"row {k + 1} {n!r}" for k, n in enumerate(aln.names[2:], 2))
        raise ValueError(
            f"{path}: an alignment of two chains has two rows (A's chain, then B's); got {len(aln.rows)}"
            + (f": {extra} is one too many" if extra else "")
        )
    for k, name, row, seq, chain, side in zip((1, 2), aln.names, aln.rows, (seq_a, seq_b), (chain_a, chain_b), "AB"):
        plain = row.replace(GAP, "")
        if plain != seq:
            at = next((p for p, (x, y) in enumerate(zip(plain, seq)) if x != y), min(len(plain), len(seq)))
            raise ValueError(
                f"{path}: row {k} {name!r} is not chain {chain}'s sequence in {side} "
                f"({len(plain)} residues in the row, {len(seq)} modelled; they differ from residue {at + 1})"
            )
    out, i, j = [], 0, 0
    for x, y in zip(*aln.rows):
        if x != GAP and y != GAP:
            out.append((i, j))
        i += x != GAP
        j += y != GAP
    return out


def _identity(names_a, names_b, pairs) -> float:
    return sum(names_a[i] == names_b[j] for i, j in pairs) / len(pairs) if pairs else 0.0


def align_chains(a: Backbone, chain_a: str, b: Backbone, chain_b: str, match: str = "auto") -> ChainAlignment:
    """Align chain_b of B to chain_a of A (the reference). match: auto (sequence at ≥ SEQUENCE_IDENTITY identity
    over the sequence alignment's pairs when those pairs also superpose at TM-score ≥ SEQUENCE_ROUTE_TM, else
    structure), sequence, structure, or the path of a two-row alignment (FASTA, Clustal or Stockholm; A's chain
    first). Residue indices count each chain's modelled residues from 0."""
    mine_a, mine_b = _chain(a, chain_a, "A"), _chain(b, chain_b, "B")
    names_a = [a.labels[k].name for k in mine_a]
    names_b = [b.labels[k].name for k in mine_b]
    if match in ("auto", "sequence", "structure"):
        from .figcheck import aligned_positions

        pairs, method, scored = aligned_positions(names_a, names_b), "sequence", None
        by_sequence = match == "sequence"
        if match == "auto" and _identity(names_a, names_b, pairs) >= SEQUENCE_IDENTITY:
            # a gapped alignment of two unrelated sequences can reach the identity (256B on 1MBN: 0.44 over 97
            # pairs): the sequence route is taken only when its pairs also superpose as one fold (13a-10)
            scored = tm_score(a.ca[mine_a], b.ca[mine_b], pairs, len(mine_a))
            by_sequence = scored[0] >= SEQUENCE_ROUTE_TM
        if not by_sequence:
            pairs, scored = _structure_pairs(a.only([chain_a]), b.only([chain_b]), pairs)
            method = "structure"
    else:
        from .sequence import one_letter

        pairs = _file_pairs(str(match), one_letter(names_a), one_letter(names_b), chain_a, chain_b)
        method, scored = "file", None
    pairs = tuple((int(i), int(j)) for i, j in pairs)
    tm, r, t, rmsd = scored or tm_score(a.ca[mine_a], b.ca[mine_b], pairs, len(mine_a))
    return ChainAlignment(
        pairs=pairs,
        tm_score=tm,
        rmsd=rmsd,
        identity=_identity(names_a, names_b, pairs),
        method=method,
        rotation=r,
        translation=t,
        residues_a=tuple(a.labels[k] for k in mine_a),
        residues_b=tuple(b.labels[k] for k in mine_b),
    )


@dataclass(frozen=True)
class Mates:
    mate: dict[str, str]  # B element id -> A element id (one to one)
    pieces: dict[str, str]  # B element id -> the A element it is a piece of (a split helix)
    extras: tuple[str, ...]  # B element ids with no mate, N to C
    only_a: tuple[str, ...]  # A element ids with no mate, N to C
    same_fold: bool  # tm_score >= SAME_FOLD_TM and >= SAME_FOLD_FRACTION of the smaller chain's elements matched
    smaller: int  # the elements of the chain with fewer (aligned chains only): what SAME_FOLD_FRACTION is a share of


_PAIRS_AS = {"G": "H"}  # a 3-10 helix (G) and an alpha helix (H) may mate; a strand (E) only a strand


def _chain_elements(sses: list[SSE], residues: tuple[ResLabel, ...]) -> list[tuple[int, int, SSE]]:
    """The elements of the aligned chain as (first, last residue index in that chain, element), N to C (ties: by
    last, then the given order). An element's start/end count all its structure's residues; the place of its first
    residue among `residues` is its index in the chain. Elements whose first residue is not among them (other
    chains) are left out. With no residues, start/end are read as the chain's indices."""
    if not residues:
        found = [(s.start, s.end, k, s) for k, s in enumerate(sses)]
    else:
        at: dict[ResLabel, int] = {}
        for k, label in enumerate(residues):
            at.setdefault(label, k)
        found = [(at[s.first], at[s.first] + len(s) - 1, k, s) for k, s in enumerate(sses) if s.first in at]
    return [(i, j, s) for i, j, _, s in sorted(found, key=lambda x: x[:3])]


def element_mates(sses_a: list[SSE], sses_b: list[SSE], aln: ChainAlignment) -> Mates:
    """Which elements of B's aligned chain are which of A's (the reference). Each B element's residues are taken to
    A's through aln.pairs (the first and last A residue paired with one of them; none paired: no place). Then
    model.pair_places pairs them one to one: the same kind (H and G alike), sharing at least half of the shorter
    element, the largest overlaps first, ties to the earlier element (N to C). A B helix left over whose place lies
    in a matched A helix (at least half of it, the largest overlap, then the earlier helix) is a piece of that
    helix (a helix B splits in two at a kink, A does not). The rest are B's extras and A's elements alone. One fold
    (H1): aln.tm_score >= SAME_FOLD_TM and the mates at least SAME_FOLD_FRACTION of the elements of the chain that
    has fewer. Only the aligned chains' elements count (aln.residues_a/_b; see _chain_elements)."""
    els_a, els_b = _chain_elements(sses_a, aln.residues_a), _chain_elements(sses_b, aln.residues_b)
    to_a = {j: i for i, j in aln.pairs}
    placed = []  # (index in els_b, (kind, first, last in A's indices))
    for k, (first, last, s) in enumerate(els_b):
        got = [to_a[j] for j in range(first, last + 1) if j in to_a]
        if got:
            placed.append((k, (_PAIRS_AS.get(s.kind, s.kind), min(got), max(got))))
    ys = [(_PAIRS_AS.get(s.kind, s.kind), first, last) for first, last, s in els_a]
    paired = pair_places([p for _, p in placed], ys)
    mate = {placed[x][0]: y for x, y in sorted(paired.items())}
    matched = set(mate.values())
    pieces: dict[int, int] = {}
    for k, (kind, first, last) in placed:
        if k in mate or kind != "H":
            continue
        inside = [
            (-overlap((first, last), y[1:]), n)
            for n, y in enumerate(ys)
            if n in matched and y[0] == "H" and 2 * overlap((first, last), y[1:]) >= last - first + 1
        ]
        if inside:
            pieces[k] = min(inside)[1]
    smaller = min(len(els_a), len(els_b))
    return Mates(
        mate={els_b[k][2].id: els_a[n][2].id for k, n in sorted(mate.items())},
        pieces={els_b[k][2].id: els_a[n][2].id for k, n in sorted(pieces.items())},
        extras=tuple(s.id for k, (_, _, s) in enumerate(els_b) if k not in mate and k not in pieces),
        only_a=tuple(s.id for n, (_, _, s) in enumerate(els_a) if n not in matched),
        same_fold=aln.tm_score >= SAME_FOLD_TM and len(mate) >= SAME_FOLD_FRACTION * smaller,
        smaller=smaller,
    )
