# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Residue alignment of one chain of a structure A (the reference) with one chain of a structure B (plan 13a).

Three routes give the residue pairs: by sequence (the global alignment A19 uses, figcheck.aligned_positions), by
structure (a TM-align-style aligner; not yet) or from a file (a two-row alignment, the escape hatch for any pair the
aligner gets wrong). Every route then superposes B on A over its pairs with a proper rotation (symmetry._kabsch,
never a mirror) and scores the pairs with the TM-score, normalised by A's chain length.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .model import Backbone

GRID = 1024  # distances are snapped to 1/GRID Å before scoring: exact in binary floating point, the same on any CPU


@dataclass(frozen=True)
class ChainAlignment:
    pairs: tuple[tuple[int, int], ...]  # (index into A's residues of chain_a, index into B's of chain_b), increasing
    tm_score: float  # normalised by A's chain length (the reference)
    rmsd: float  # over `pairs` after superposition
    identity: float  # identical residues / len(pairs)
    method: str  # "sequence" | "structure" | "file"
    rotation: np.ndarray = field(compare=False)  # 3x3 proper rotation, A ≈ R·B + t
    translation: np.ndarray = field(compare=False)  # (3,)


def d0(length: int) -> float:
    """TM-score's distance scale d0(L) = 1.24·(L−15)^(1/3) − 1.8, at least 0.5 (Zhang & Skolnick 2004)."""
    if length <= 15:  # the cube root of a negative number is negative: below 0.5 anyway
        return 0.5
    return max(0.5, 1.24 * (length - 15) ** (1 / 3) - 1.8)


def tm_score(ca_a, ca_b, pairs, length_a) -> tuple[float, np.ndarray, np.ndarray, float]:
    """(TM-score, R, t, RMSD) of B's Cα superposed on A's over `pairs` ((i in A, j in B)), with A ≈ R·B + t a proper
    rotation (symmetry._kabsch). Each pair scores 1/(1+(d/d0)²), d snapped to 1/GRID Å, summed over length_a."""
    from .symmetry import _kabsch

    if not len(pairs):
        return 0.0, np.eye(3), np.zeros(3), 0.0
    ia, ib = (np.array(k, dtype=int) for k in zip(*pairs))
    pa, pb = np.asarray(ca_a, float)[ia], np.asarray(ca_b, float)[ib]
    r, t, rmsd = _kabsch(pb, pa)
    d = np.round(np.linalg.norm(pb @ r.T + t - pa, axis=1) * GRID) / GRID
    scale = d0(length_a)
    tm = math.fsum((1.0 / (1.0 + (d / scale) ** 2)).tolist()) / length_a  # fsum: one rounding, any summation order
    return tm, r, t, rmsd


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


def align_chains(a: Backbone, chain_a: str, b: Backbone, chain_b: str, match: str = "auto") -> ChainAlignment:
    """Align chain_b of B to chain_a of A (the reference). match: auto, sequence, structure, or the path of a
    two-row alignment (FASTA, Clustal or Stockholm; A's chain first). Residue indices count each chain's modelled
    residues from 0."""
    mine_a, mine_b = _chain(a, chain_a, "A"), _chain(b, chain_b, "B")
    names_a = [a.labels[k].name for k in mine_a]
    names_b = [b.labels[k].name for k in mine_b]
    if match == "structure":
        raise NotImplementedError("the structure route of the aligner is not implemented yet")
    if match in ("auto", "sequence"):  # auto: the sequence route until the structure route exists
        from .figcheck import aligned_positions

        pairs, method = aligned_positions(names_a, names_b), "sequence"
    else:
        from .sequence import _one_letter

        pairs = _file_pairs(str(match), _one_letter(names_a), _one_letter(names_b), chain_a, chain_b)
        method = "file"
    pairs = tuple((int(i), int(j)) for i, j in pairs)
    tm, r, t, rmsd = tm_score(a.ca[mine_a], b.ca[mine_b], pairs, len(mine_a))
    same = sum(names_a[i] == names_b[j] for i, j in pairs)
    return ChainAlignment(
        pairs=pairs,
        tm_score=tm,
        rmsd=rmsd,
        identity=same / len(pairs) if pairs else 0.0,
        method=method,
        rotation=r,
        translation=t,
    )
