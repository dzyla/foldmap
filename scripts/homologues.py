# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The homologue-pair campaign (plan 13a, spec 4.4): structure B drawn like a reference A (`--like`), judged by
H1-H5 (foldmap.criteria.like_checks) beside B drawn on its own (needs network for the PDB entries; not run in CI).

    python scripts/homologues.py align --set docs/validation/homologues/tuning.json --out DIR
                                       [--usalign PATH|auto] [--jobs 1]
    python scripts/homologues.py usalign-report tuning=DIR1 heldout=DIR2 --out docs/validation/homologues
                                       --summary docs/validation/homologues/heldout-summary.json --commit SHA

align aligns each pair's two chains (foldmap.align.align_chains, auto route; as deposited) and element-mates
them, without drawing: per pair the method, TM-score (normalised by A), RMSD, identity, aligned pairs, mates and
same_fold, to DIR/align.json and a table on stdout. With --usalign (the binary, or auto: find_usalign) each pair
is also compared with US-align (Task 13a-3, compare_usalign): foldmap's structure route against US-align on the
same residues, both TM-scores normalised by A, and the share of US-align's residue pairs foldmap reproduces within
±1. usalign-report gathers such runs into usalign.json, usalign-report.md and the summary's 'usalign' entry
(write_usalign_report).

    python scripts/homologues.py run --set docs/validation/homologues/tuning.json --out DIR [--jobs 4]
                                     [--pairs 2VWD_A-2ZB5_A,...]

run draws every pair (or the --pairs named, by stem A_CHAIN-B_CHAIN) with `pair`, a child process each (python -I,
PYTHONHASHSEED=0, OPENBLAS_NUM_THREADS=1), then writes the report as `report` does. `pair` draws B like A as
`foldmap check B --like A --like-chain CA --chains CB --assembly asu` does (each side as deposited: the set
names deposited chains), and B alone (its same chain), judges both (criteria.check_figure: the S- and A-checks)
and the pair (criteria.like_checks, with the pair's kind for H1), and writes DIR/<stem>.json and two pages side by
side (compare's page): DIR/<stem>-like.png (A | B drawn like A) and DIR/<stem>-alone.png (A | B drawn on its own).

    python scripts/homologues.py sample --seed 20261021 --out docs/validation/homologues/heldout.json

sample draws the held-out set (spec 4.4) from CATH's domain list (release CATH_RELEASE, downloaded once into
~/.cache/foldmap/cath; CATH data CC BY 4.0, https://www.cathdb.info) with a seeded random.Random, in three strata
drawn in this order: 20 "remote" pairs (one topology; two S35 clusters, of different superfamilies, or of one
superfamily at structure-alignment identity under 30 %), 10 "close" pairs (one superfamily and S35 cluster, two S95
clusters, identity 30-90 %) and 10 negatives (one architecture, two topologies). Only whole-chain domains (CATH's
domain number 00: the chain is one domain) of 60-600 residues, of current entries, whose entries are not in the
tuning set and not already in the held-out set. Each level (topology, cluster, domain) is chosen uniformly among the
sorted keys. A drawn pair is kept when both chains are read (foldmap.fetch, --assembly asu) with 60-600 residues
and the stratum's identity holds (foldmap.align.align_chains, structure route: identical / aligned pairs), else it
is recorded under "rejected" with the reason and the next pair drawn. Positives are "remote" and "close". The set
is reproducible given the same fetches: a fetch that fails rejects that draw, so on another day the same seed can
give another set (the committed heldout.json is the one measured).

    python scripts/homologues.py report --set docs/validation/homologues/tuning.json --out DIR [--md PATH]

report reads DIR/<stem>.json of every pair and writes the table (DIR/report.md, or --md): per pair H1-H5, H2's and
H3's values drawn alike and drawn independently, the checks B drawn like A loses, and the failures by check. When
the set's pairs carry a stratum (the held-out set), the table gains the stratum, and the report a summary (passes
per check and stratum; the medians drawn independently and drawn alike), also written as JSON (DIR/summary.json).

    python scripts/homologues.py diff --out DIR [--jobs 4] [--pairs 1RM5-2HKI,...] [--md PATH]

diff (plan 13e-3) draws the two-state pairs as `foldmap diff A B` draws them (compare.diff_page: publication look,
auto route and assembly, the pair's domains), a child process each (diff-pair). The pairs are the random-100
samples' S11 pairs (seed as A, partner as B, each its first protein chain) and docs/validation/diff/hinges.json
(published domain motions; synthetic hinges on 1UBQ, written into DIR). For each pair it writes DIR/<stem>.png and
.svg, the --json output (DIR/<stem>-diff.json), and DIR/<stem>.json with: H5 and the S-/A-checks of each panel
(compare.panel_checks); the Diff's numbers at each core-fit cut of CUTS (core kept and its share per domain, domain
turns, each two domains' relative rotation fitted on their own); and the comparison with what the sample records
of the pair. Then the tables (write_diff_report; diff-report rewrites them from DIR).
"""

from __future__ import annotations

import argparse
import json
import os
import random
import statistics
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
CHECKS = ("H1", "H2", "H3", "H4", "H5")
PAIR_TIMEOUT_S = 1800  # a pair's child is killed after this many seconds
TUNING = ROOT / "docs" / "validation" / "homologues" / "tuning.json"
CATH_RELEASE = "v4_3_0"  # CATH's domain list and names of this release (fixed: the sample does not move with CATH)
CATH_URL = (
    "https://download.cathdb.info/cath/releases/all-releases/{release}/cath-classification-data/{name}-{release}.txt"
)
STRATA = (("remote", 20, "positive"), ("close", 10, "positive"), ("negative", 10, "negative"))  # in draw order
CHAIN_LENGTH = (60, 600)  # residues of a sampled chain (CATH's domain length, and as foldmap reads it)
REMOTE_IDENTITY = 0.30  # a remote pair of one superfamily is kept under this identity
CLOSE_IDENTITY = (0.30, 0.90)  # a close pair is kept at this identity (both inclusive)
MAX_DRAWS = 400  # draws per stratum at most


def read_set(path) -> list[dict]:
    """The pairs of a set file ({"pairs": [{a: {id, chain}, b: {id, chain}, kind, note}, ...]}), as written."""
    return json.loads(Path(path).read_text(encoding="utf-8"))["pairs"]


def stem(pair: dict) -> str:
    """A pair's file stem: A_CHAIN-B_CHAIN."""
    return f"{pair['a']['id']}_{pair['a']['chain']}-{pair['b']['id']}_{pair['b']['chain']}"


def _round(x) -> float | None:
    return None if x is None else float(round(float(x), 4)) + 0.0


# -------------------------------------------------------------------------------------------------------- align


def _sides(pair: dict) -> list[tuple]:
    """Each side's (elements, backbone) as `foldmap plot B --like A` reads them (make_layout's _drawn_like): each
    chain laid out alone first, for its elements (as deposited)."""
    from foldmap.cli import make_layout, resolve_structure
    from foldmap.io import read_model

    sides = []
    for side in ("a", "b"):
        path = resolve_structure(pair[side]["id"])
        model = read_model(path, "asu")
        lay, sses, bb = make_layout(path, chains=[pair[side]["chain"]], model=model, assembly="asu")
        sides.append((sses, bb))
    return sides


def _align(pair: dict, usalign: str = "") -> dict:
    """The pair aligned and mated as `foldmap plot B --like A` aligns it (auto route); with `usalign` (the US-align
    binary) also compared with US-align (compare_usalign)."""
    from foldmap.align import align_chains, element_mates

    (sses_a, bb_a), (sses_b, bb_b) = _sides(pair)
    aln = align_chains(bb_a, pair["a"]["chain"], bb_b, pair["b"]["chain"], "auto")
    mates = element_mates(sses_a, sses_b, aln)
    extra = {}
    if usalign:
        extra["usalign"] = compare_usalign(bb_a, pair["a"]["chain"], bb_b, pair["b"]["chain"], usalign, aln)
    return {
        **extra,
        "pair": stem(pair),
        "kind": pair["kind"],
        "method": aln.method,
        "tm_score": _round(aln.tm_score),
        "rmsd": _round(aln.rmsd),
        "identity": _round(aln.identity),
        "aligned": len(aln.pairs),
        "length_a": len(aln.residues_a),
        "length_b": len(aln.residues_b),
        "mates": len(mates.mate),
        "pieces": len(mates.pieces),
        "extras": len(mates.extras),
        "only_a": len(mates.only_a),
        "same_fold": mates.same_fold,
    }


def cmd_align(args) -> int:
    pairs = read_set(args.set)
    exe = ""
    if args.usalign:
        exe = find_usalign(None if args.usalign == "auto" else args.usalign)
        if not exe:
            raise SystemExit("--usalign auto: no US-align found (USALIGN, PATH, ~/foldmap-remote/usalign-env/bin)")
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        rows = list(pool.map(_align, pairs, [exe] * len(pairs)))
    for p, r in zip(pairs, rows):
        if "stratum" in p:
            r["stratum"] = p["stratum"]
        r["source_stamp"] = _stamp()  # align.json: each pair with the code that aligned it
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "align.json").write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print("pair\tkind\tmethod\tTM\tRMSD\tidentity\taligned\tmates/pieces/extras/only_a\tsame_fold")
    for r in rows:
        print(
            f"{r['pair']}\t{r['kind']}\t{r['method']}\t{r['tm_score']:.4f}\t{r['rmsd']:.2f}\t{r['identity']:.3f}\t"
            f"{r['aligned']}/{r['length_a']}\t{r['mates']}/{r['pieces']}/{r['extras']}/{r['only_a']}\t{r['same_fold']}"
        )
    if exe:
        print("\npair\tUS-align TM\tfoldmap structure TM\tpairs within ±1\tlimits")
        for r in rows:
            u = r["usalign"]
            print(f"{r['pair']}\t{u['usalign']['tm_a']:.4f}\t{u['structure']['tm_score']:.4f}\t"
                  f"{u['agreement']['within_1']:.3f}\t{u['limits']['verdict']}")  # fmt: skip
    return 0


# ------------------------------------------------------------------------------------------------------ US-align

# Task 13a-3's comparison of foldmap's structure route with US-align (Zhang et al. 2022). Both TM-scores are
# normalised by A's chain (the reference; d0 of A's length, as foldmap's). The limits apply where US-align's TM-score
# is at least USALIGN_TM: foldmap's within TM_TOLERANCE of it, and at least PAIR_AGREEMENT of US-align's residue pairs
# (i, j) reproduced, foldmap pairing A's residue i with one of B's within PAIR_SLACK residues of j.
USALIGN_TM = 0.5
TM_TOLERANCE = 0.05
PAIR_AGREEMENT = 0.85
PAIR_SLACK = 1
USALIGN_ARGS = ("-mol", "prot", "-mm", "0", "-ter", "2")  # one protein chain against one (each file holds one)


def find_usalign(path: str | None = None) -> str:
    """The US-align binary: `path`, else $USALIGN, else USalign on PATH, else zyla-office's
    ~/foldmap-remote/usalign-env/bin/USalign; '' when none is there."""
    import shutil

    for p in (path, os.environ.get("USALIGN"), shutil.which("USalign"),
              str(Path.home() / "foldmap-remote" / "usalign-env" / "bin" / "USalign")):  # fmt: skip
        if p and Path(p).is_file() and os.access(p, os.X_OK):
            return str(p)
    return ""


def write_ca_pdb(labels, ca, path: Path) -> None:
    """One chain's Cα as a PDB file US-align reads: the residues exactly as foldmap aligns them (its Backbone of the
    chain), numbered 1..n in their order (no insertion codes, nothing left out), chain A, then TER."""
    lines = [
        f"ATOM  {k + 1:5d}  CA  {(label.name or 'UNK')[:3]:>3s} A{k + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
        f"           C"
        for k, (label, (x, y, z)) in enumerate(zip(labels, ca))
    ]
    path.write_text("\n".join([*lines, "TER", "END"]) + "\n", encoding="ascii")


def parse_usalign(text: str) -> dict:
    """US-align's default output for two chains: its version, both lengths, the TM-score normalised by each
    structure, the RMSD and aligned length it reports, and the alignment: `rows` (structure 1, the marks, structure
    2) and `pairs`, the (i, j) of every column where both rows hold a residue (':' and '.' alike: the pairs it
    aligns), with `close` those marked ':' (under 5 Å). i and j count each structure's residues from 0."""
    import re

    got: dict = {"tm": {}}
    lines = text.splitlines()
    for line in lines:
        if m := re.search(r"US-align \(Version (\S+)\)", line):
            got["version"] = m.group(1)
        elif m := re.match(r"Length of Structure_([12]):\s+(\d+) residues", line):
            got[f"length_{m.group(1)}"] = int(m.group(2))
        elif m := re.match(r"Aligned length=\s*(\d+), RMSD=\s*([\d.]+)", line):
            got["aligned"], got["rmsd"] = int(m.group(1)), float(m.group(2))
        elif m := re.match(r"TM-score=\s*([\d.]+) \(normalized by length of Structure_([12])", line):
            got["tm"][int(m.group(2))] = float(m.group(1))
    at = next(k for k, line in enumerate(lines) if line.startswith('(":" denotes'))
    rows = [line.rstrip("\n") for line in lines[at + 1 : at + 4]]
    if len(rows) != 3 or len(rows[0]) != len(rows[2]):
        raise ValueError("US-align's alignment rows not found")
    rows[1] = rows[1].ljust(len(rows[0]))
    pairs, close, i, j = [], [], 0, 0
    for x, mark, y in zip(*rows):
        if x != "-" and y != "-":
            pairs.append((i, j))
            if mark == ":":
                close.append((i, j))
        i += x != "-"
        j += y != "-"
    if i != got.get("length_1") or j != got.get("length_2"):
        raise ValueError(f"US-align's rows hold {i} and {j} residues, not {got.get('length_1')}, {got.get('length_2')}")
    return {**got, "rows": rows, "pairs": pairs, "close": close}


def run_usalign(exe: str, labels_a, ca_a, labels_b, ca_b) -> dict:
    """US-align on A's chain (structure 1) and B's (structure 2), each written by write_ca_pdb; parse_usalign's
    result."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        fa, fb = Path(tmp) / "A.pdb", Path(tmp) / "B.pdb"
        write_ca_pdb(labels_a, ca_a, fa)
        write_ca_pdb(labels_b, ca_b, fb)
        done = subprocess.run([exe, str(fa), str(fb), *USALIGN_ARGS], capture_output=True, text=True, check=True)
    return parse_usalign(done.stdout)


def pair_agreement(ref, got, slack: int = PAIR_SLACK) -> float | None:
    """The share of the pairs `ref` ((i, j): A's residue i with B's j) that `got` reproduces: got pairs A's i with a
    residue of B within `slack` of j. None when ref is empty."""
    to_b = dict(got)
    if not ref:
        return None
    return sum(i in to_b and abs(to_b[i] - j) <= slack for i, j in ref) / len(ref)


def disagreements(ref, got, labels_a, slack: int = PAIR_SLACK) -> list[str]:
    """Where `got` departs from `ref`, as runs of A's residues (ref's pairs in order, A's numbers as written): 'A
    12-30: B shifted +4' (got pairs them 4 residues further along B) or 'A 12-30: not paired' (got leaves them out)."""
    to_b = dict(got)
    runs: list[list] = []  # [first i, last i, what]
    for i, j in ref:
        what = "not paired" if i not in to_b else (None if abs(to_b[i] - j) <= slack else f"B shifted {to_b[i] - j:+d}")
        if what is None:
            continue
        if runs and runs[-1][2] == what and runs[-1][1] == i - 1:
            runs[-1][1] = i
        else:
            runs.append([i, i, what])
    return [f"A {labels_a[x].tag}-{labels_a[y].tag} ({y - x + 1}): {what}" for x, y, what in runs]


def gapped_rows(pairs, seq_a: str, seq_b: str) -> list[str]:
    """Two gapped rows (A, then B) of a pairwise alignment given by its pairs, as US-align writes it: before each pair
    the residues of A, then of B, that no pair holds."""
    ra, rb, i, j = [], [], 0, 0
    for x, y in [*pairs, (len(seq_a), len(seq_b))]:
        ra.append(seq_a[i:x] + "-" * (y - j))
        rb.append("-" * (x - i) + seq_b[j:y])
        if x < len(seq_a):
            ra.append(seq_a[x])
            rb.append(seq_b[y])
        i, j = x + 1, y + 1
    return ["".join(ra), "".join(rb)]


def compare_usalign(bb_a, chain_a: str, bb_b, chain_b: str, exe: str, auto=None) -> dict:
    """foldmap's structure route (align_chains, match 'structure') on the chains against US-align on the same
    residues: both TM-scores normalised by A, foldmap's TM-score of US-align's own pairs (scoring against
    alignment), the share of US-align's pairs foldmap reproduces within PAIR_SLACK (and exactly, and of the pairs
    under 5 Å), the limits (module constants) and, where the pairs disagree, the runs of A's residues (disagreements)
    and both alignments' rows. `auto`: the auto route's alignment, whose agreement is recorded too."""
    from foldmap.align import align_chains, tm_score
    from foldmap.sequence import one_letter

    aln = align_chains(bb_a, chain_a, bb_b, chain_b, "structure")
    la, lb = list(aln.residues_a), list(aln.residues_b)
    ka = [k for k, label in enumerate(bb_a.labels) if label.chain == chain_a]
    kb = [k for k, label in enumerate(bb_b.labels) if label.chain == chain_b]
    ca_a, ca_b = bb_a.ca[ka], bb_b.ca[kb]
    us = run_usalign(exe, la, ca_a, lb, ca_b)
    ref, got = us["pairs"], list(aln.pairs)
    share = pair_agreement(ref, got)
    tm_us, tm_fm = us["tm"][1], float(aln.tm_score)
    applies = tm_us >= USALIGN_TM
    tm_ok = abs(tm_fm - tm_us) <= TM_TOLERANCE
    pairs_ok = share is not None and share >= PAIR_AGREEMENT
    verdict = "pass" if tm_ok and pairs_ok else "miss"
    seq_a, seq_b = one_letter([x.name for x in la]), one_letter([x.name for x in lb])
    return {
        "usalign": {
            "version": us.get("version"),
            "tm_a": _round(tm_us),
            "tm_b": _round(us["tm"][2]),
            "rmsd": _round(us["rmsd"]),
            "aligned": len(ref),
            "aligned_reported": us["aligned"],
            "close": len(us["close"]),
            "foldmap_tm_of_its_pairs": _round(tm_score(ca_a, ca_b, ref, len(la))[0]),
        },
        "structure": {"tm_score": _round(tm_fm), "rmsd": _round(aln.rmsd), "aligned": len(got)},
        "auto": None
        if auto is None
        else {
            "method": auto.method,
            "tm_score": _round(auto.tm_score),
            "within_1": _round(pair_agreement(ref, list(auto.pairs))),
        },  # fmt: skip
        "tm_difference": _round(tm_fm - tm_us),
        "agreement": {
            "within_1": _round(share),
            "exact": _round(pair_agreement(ref, got, 0)),
            "close_within_1": _round(pair_agreement(us["close"], got)),
            "foldmap_pairs_in_usalign": _round(pair_agreement(got, ref)),
        },
        "limits": {"applies": applies, "tm": tm_ok, "pairs": pairs_ok, "verdict": verdict if applies else "–"},
        "disagreements": disagreements(ref, got, la),
        "rows": {"usalign": [us["rows"][0], us["rows"][2]], "foldmap": gapped_rows(got, seq_a, seq_b)},
    }


def usalign_summary(rows: list[dict]) -> dict:
    """The comparison over rows (cmd_align's, each with its 'usalign' entry): the pairs; those US-align puts at TM
    USALIGN_TM or more (where the limits apply), how many of them meet the TM limit, the pair limit and both, and the
    median share of pairs reproduced; over all pairs, the median TM difference (foldmap - US-align) and how many fall
    on the same side of SAME_FOLD_TM in both aligners (H1's cut)."""
    from foldmap.align import SAME_FOLD_TM

    us = [r["usalign"] for r in rows]
    judged = [u for u in us if u["limits"]["applies"]]
    return {
        "pairs": len(us),
        "judged": len(judged),
        "tm_within": sum(u["limits"]["tm"] for u in judged),
        "pairs_agree": sum(u["limits"]["pairs"] for u in judged),
        "both": sum(u["limits"]["tm"] and u["limits"]["pairs"] for u in judged),
        "within_1_median": _median([u["agreement"]["within_1"] for u in judged]),
        "within_1_min": _round(min((u["agreement"]["within_1"] for u in judged), default=None)),
        "tm_difference_median": _median([u["tm_difference"] for u in us]),
        "tm_difference_max_abs": _round(max((abs(u["tm_difference"]) for u in us), default=0.0)),
        "same_side_of_h1": sum(
            (u["structure"]["tm_score"] >= SAME_FOLD_TM) == (u["usalign"]["tm_a"] >= SAME_FOLD_TM) for u in us
        ),  # fmt: skip
        "foldmap_lower": sum(u["tm_difference"] < 0 for u in us),
        # foldmap's TM-score of US-align's own pairs minus US-align's: the scoring, apart from the alignment search
        "scoring_difference_max_abs": _round(
            max((abs(u["usalign"]["foldmap_tm_of_its_pairs"] - u["usalign"]["tm_a"]) for u in us), default=0.0)
        ),
        # pairs where US-align's pairs score higher in foldmap's own TM-score than foldmap's pairs (by > 0.005):
        # foldmap's search stopped short of an alignment it would itself prefer
        "search_short": sum(u["usalign"]["foldmap_tm_of_its_pairs"] > u["structure"]["tm_score"] + 0.005 for u in us),
    }


def _h1_refused(rows: list[dict]) -> list[dict]:
    """The positives foldmap does not call one fold (the auto route's same_fold false: H1 refuses them)."""
    return [r for r in rows if r["kind"] == "positive" and not r["same_fold"]]


_STAMP: dict | None = None


def _stamp() -> dict:
    """The code that measured (foldmap.source_stamp, plan 28 Task 1), stamped into each result as it is written: taken
    once, when main() starts (review M8), and reused."""
    global _STAMP
    if _STAMP is None:
        from foldmap.source_stamp import source_stamp

        _STAMP = source_stamp(__file__)
    return _STAMP


def write_usalign_report(runs: dict[str, list[dict]], out: Path, summary_path: Path | None, commit: str) -> Path:
    """The comparison's files from align --usalign runs (`runs`: set name -> its align.json rows): out/usalign.json
    (every row, per set), out/usalign-report.md (per pair: both TM-scores, the pairs reproduced, the limits; the
    misses with where the alignments part; the held-out positives H1 refuses beside US-align's TM-score) and, with
    summary_path (heldout-summary.json), its 'usalign' entry (usalign_summary per set and over all, the H1 refusals)."""
    from foldmap.align import SAME_FOLD_TM
    from foldmap.source_stamp import code_mix

    every = [r for rows in runs.values() for r in rows]
    version = next((r["usalign"]["usalign"]["version"] for r in every), None)
    how = (
        f"scripts/homologues.py align --usalign at {commit}: foldmap's structure route (align_chains, match "
        f"'structure') against US-align {version} (`USalign A.pdb B.pdb {' '.join(USALIGN_ARGS)}`) on the same "
        "residues (each chain's Cα as foldmap reads it, --assembly asu). TM-scores normalised by A's chain (US-align: "
        "by Structure_1, A; d0 of A's length, as foldmap's). US-align's pairs: every column of its alignment holding "
        f"a residue of both. Limits where US-align's TM >= {USALIGN_TM}: foldmap's TM within {TM_TOLERANCE}, and >= "
        f"{PAIR_AGREEMENT:.0%} of US-align's pairs reproduced within ±{PAIR_SLACK} residue (13a-3)."
    )
    sums = {name: usalign_summary(rows) for name, rows in runs.items()}
    sums["all"] = usalign_summary(every)
    held = runs.get("heldout", [])
    refused = [
        {
            "pair": r["pair"],
            "stratum": r.get("stratum"),
            "foldmap_tm": r["tm_score"],
            "usalign_tm": r["usalign"]["usalign"]["tm_a"],
            "usalign_one_fold": r["usalign"]["usalign"]["tm_a"] >= SAME_FOLD_TM,
        }
        for r in _h1_refused(held)
    ]
    out.mkdir(parents=True, exist_ok=True)
    data = {"about": how, "measured_at": commit, "usalign_version": version, "sets": runs, "source_stamp": _stamp()}
    data |= code_mix(r.get("source_stamp") for r in every)  # the aligning code, per pair (ruling Q13)
    (out / "usalign.json").write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    lines = ["# foldmap's aligner against US-align", "", how, ""]
    lines += [
        "Columns: US-align's TM-score; foldmap's (structure route) and the difference; foldmap's TM-score of "
        "US-align's own pairs (scoring apart from search); the auto route's TM-score as the campaign used it; the "
        f"share of US-align's pairs foldmap reproduces within ±{PAIR_SLACK} (exactly; of those under 5 Å); the limits "
        "('–': US-align's TM under the limits' cut). What the numbers say (the misses, the H1 refusals, the "
        "aligner's search) is in review.md.",
        "",
    ]
    for name, rows in runs.items():
        s = sums[name]
        lines += [
            f"## {name}",
            "",
            "| pair | kind | US-align TM | foldmap TM | diff. | foldmap TM of US pairs | auto route, TM | pairs ±1 "
            "| exact | under 5 Å ±1 | limits |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        for r in rows:
            u = r["usalign"]
            kind = r["kind"] + (f" ({r['stratum']})" if r.get("stratum") else "")
            lines.append(
                f"| {r['pair']} | {kind} | {u['usalign']['tm_a']:.3f} | {u['structure']['tm_score']:.3f} | "
                f"{u['tm_difference']:+.3f} | {u['usalign']['foldmap_tm_of_its_pairs']:.3f} | {r['method']}, "
                f"{r['tm_score']:.3f} | {u['agreement']['within_1']:.3f} | {u['agreement']['exact']:.3f} | "
                f"{(u['agreement']['close_within_1'] or 0):.3f} | {u['limits']['verdict']} |"
            )
        lines += [
            "",
            f"- {s['judged']} of {s['pairs']} pairs at US-align TM >= {USALIGN_TM}: TM within {TM_TOLERANCE} "
            f"{s['tm_within']}, pairs >= {PAIR_AGREEMENT:.0%} {s['pairs_agree']}, both **{s['both']}**; pairs "
            f"reproduced (±{PAIR_SLACK}): median {s['within_1_median']}, least {s['within_1_min']}.",
            f"- All {s['pairs']}: TM difference (foldmap - US-align) median {s['tm_difference_median']}, largest "
            f"{s['tm_difference_max_abs']}; foldmap lower in {s['foldmap_lower']}; same side of H1's "
            f"{SAME_FOLD_TM} in both: {s['same_side_of_h1']}; foldmap's own score of US-align's pairs above its "
            f"own alignment's (search short): {s['search_short']}; scoring difference on US-align's pairs at most "
            f"{s['scoring_difference_max_abs']}.",
            "",
        ]
        misses = [r for r in rows if r["usalign"]["limits"]["verdict"] == "miss"]
        for r in misses:
            u = r["usalign"]
            lines += [f"### Miss: {r['pair']}", "", "Where foldmap's pairs part from US-align's (A's residues):", ""]
            lines += [f"- {d}" for d in u["disagreements"]] + ["", "```"]
            for tag in ("usalign", "foldmap"):
                ra, rb = u["rows"][tag]
                for k in range(0, len(ra), 100):
                    lines += [f"{tag:8s} A {ra[k : k + 100]}", f"{tag:8s} B {rb[k : k + 100]}", ""]
            lines += ["```", ""]
    if refused:
        under = sum(not x["usalign_one_fold"] for x in refused)
        lines += [
            "## The held-out positives H1 refuses",
            "",
            "| pair | stratum | foldmap TM (auto) | US-align TM | US-align: one fold |",
            "|---|---|---|---|---|",
            *[
                f"| {x['pair']} | {x['stratum']} | {x['foldmap_tm']:.3f} | {x['usalign_tm']:.3f} | "
                f"{'yes' if x['usalign_one_fold'] else 'no'} |"
                for x in refused
            ],  # fmt: skip
            "",
            f"US-align puts {under} of {len(refused)} under TM {SAME_FOLD_TM} too.",
            "",
        ]
    md = out / "usalign-report.md"
    md.write_text("\n".join(lines), encoding="utf-8")
    if summary_path is not None:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["usalign"] = {
            "about": how + " Added after the held-out run (plan 23 Task 3): the H1-H5 numbers of this file stand as "
            "measured at its top-level measured_at; this entry's measured_at is the comparison's. Details: "
            "usalign-report.md, usalign.json.",
            "measured_at": commit,
            "version": version,
            "limits": {
                "usalign_tm": USALIGN_TM,
                "tm_tolerance": TM_TOLERANCE,
                "pair_agreement": PAIR_AGREEMENT,
                "pair_slack": PAIR_SLACK,
            },  # fmt: skip
            **sums,
            "h1_refused": refused,
            "source_stamp": _stamp(),
        }
        summary_path.write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    return md


def cmd_usalign_report(args) -> int:
    runs = {}
    for item in args.runs:
        name, _, path = item.partition("=")
        runs[name] = json.loads((Path(path) / "align.json").read_text(encoding="utf-8"))
    summary = Path(args.summary) if args.summary else None
    write_usalign_report(runs, Path(args.out), summary, args.commit)
    return 0


# --------------------------------------------------------------------------------------------------------- pair


def _summary(found: dict) -> dict:
    """check_figure's result in short: the S-checks with problems (their counts), the criteria that fail, and the
    kind of each check named (check_figure's "kinds", science.KIND: consistency, external or readability; plan 28
    Task 2)."""
    science = {k: len(v) for k, v in sorted(found["science"].items()) if v}
    failed = sorted(k for k, v in found["criteria"].items() if v.get("pass") is False)
    kinds = found.get("kinds") or {}
    return {"science": science, "failed": failed, "kinds": {k: kinds[k] for k in [*science, *failed] if k in kinds}}


def run_pair(pair: dict, out: Path) -> dict:
    """One pair drawn like and alone, judged, its pages written (module doc)."""
    from foldmap import cli
    from foldmap.compare import _page
    from foldmap.criteria import check_figure, jsonable, like_checks
    from foldmap.render import save
    from foldmap.spec import draw_spec

    a, b = pair["a"], pair["b"]
    args = cli.build_parser().parse_args(["check", b["id"], "--like", a["id"], "--like-chain", a["chain"],
                                          "--chains", b["chain"], "--assembly", "asu"])  # fmt: skip
    args.structure = cli.resolve_structure(args.structure)
    spec = cli._spec_of(args)
    t0 = time.perf_counter()
    got: dict = {}
    drawn = draw_spec(spec, like_out=got)
    t_like = time.perf_counter() - t0
    alone_spec = replace(spec, like=None)
    alone = draw_spec(alone_spec)
    found_like = check_figure(drawn, spec, spec.structure)
    found_alone = check_figure(alone, alone_spec, spec.structure)
    ref = got.get("reference")
    note = drawn.layout.like_note
    h = like_checks(
        ref[0] if ref else None,
        drawn.layout,
        got.get("alignment"),
        got.get("mates"),
        note,
        got.get("like") is not None,
        alone.layout,
        found_like,
        found_alone,
        kind=pair["kind"],
        refused=bool(got.get("refused")),
    )
    aln, mates = got.get("alignment"), got.get("mates")
    s = stem(pair)
    if ref:
        lay_a, sses_a, bb_a, model_a, chain_a = ref
        look = spec.look
        for tag, other, title in (("like", drawn, "drawn like it"), ("alone", alone, "drawn on its own")):
            panels = (
                (lay_a, sses_a, bb_a, model_a, f"{a['id']} chain {chain_a}"),
                (other.layout, other.sses, other.bb, other.model, f"{b['id']} chain {b['chain']}, {title}"),
            )
            fig, _ = _page(panels, note if tag == "like" else "", look)
            save(fig, out / f"{s}-{tag}.png")
    result = {
        "pair": s,
        "kind": pair["kind"],
        "note": note,
        "alignment": None
        if aln is None
        else {
            "method": aln.method,
            "tm_score": _round(aln.tm_score),
            "rmsd": _round(aln.rmsd),
            "identity": _round(aln.identity),
            "aligned": len(aln.pairs),
            "length_a": len(aln.residues_a),
            "length_b": len(aln.residues_b),
        },
        "mates": None
        if mates is None
        else {
            "mates": len(mates.mate),
            "pieces": len(mates.pieces),
            "extras": len(mates.extras),
            "only_a": len(mates.only_a),
            "same_fold": mates.same_fold,
        },
        "like": jsonable(h),
        "checks_like": _summary(found_like),
        "checks_alone": _summary(found_alone),
        "seconds_like": round(t_like, 1),
    }
    result["source_stamp"] = _stamp()
    (out / f"{s}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return result


def cmd_pair(args) -> int:
    pair = {
        "a": dict(zip(("id", "chain"), args.a.split(":"))),
        "b": dict(zip(("id", "chain"), args.b.split(":"))),
        "kind": args.kind,
    }
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    r = run_pair(pair, out)
    print(f"{r['pair']}: " + ", ".join(f"{k} {_word(r['like'][k])}" for k in CHECKS))
    return 0


# ---------------------------------------------------------------------------------------------------------- run


def cmd_run(args) -> int:
    pairs = read_set(args.set)
    if args.pairs:
        wanted = set(args.pairs.split(","))
        pairs = [p for p in pairs if stem(p) in wanted]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1"}

    def one(p: dict) -> tuple[str, int]:
        s = stem(p)
        cmd = [sys.executable, "-I", str(SCRIPT), "pair", f"{p['a']['id']}:{p['a']['chain']}",
               f"{p['b']['id']}:{p['b']['chain']}", "--kind", p["kind"], "--out", str(out)]  # fmt: skip
        with open(out / f"{s}.log", "w", encoding="utf-8") as log:
            try:
                done = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, timeout=PAIR_TIMEOUT_S)
                code = done.returncode
            except subprocess.TimeoutExpired:
                code = -9
        print(f"{s}: exit {code}", flush=True)
        return s, code

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        codes = dict(pool.map(one, pairs))
    write_report(read_set(args.set), out, Path(args.md) if args.md else out / "report.md")
    return 0 if all(c == 0 for c in codes.values()) else 1


# --------------------------------------------------------------------------------------------------------- diff

DIFF_SAMPLES = tuple(ROOT / "docs" / "validation" / d / "sample.json" for d in ("random100", "random100-heldout",
                                                                              "random100-third"))  # fmt: skip
HINGES = ROOT / "docs" / "validation" / "diff" / "hinges.json"
CUTS = (2.0, 1.0, 0.5)  # the core fit's cut (diff.CORE_CUT, in d0) each pair's numbers are measured at
MOVED_A = 3.0  # Å: a residue counted as moved (the table's "moved" column)


def sample_pairs(path) -> list[dict]:
    """The two-state pairs (S11) of a random-100 sample.json: the seed as A, its partner as B, each its first protein
    chain (as random-100's pair_check reads them), with what the sample records of the pair (`known`)."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    by_id = {e["id"]: e for e in data["entries"] if e.get("pair")}
    out = []
    for e in sorted((e for e in by_id.values() if e["role"] == "seed"), key=lambda e: int(e["pair"][1:])):
        p = by_id[e["partner"]]
        out.append(
            {
                "name": f"{Path(path).parent.name} {e['pair']}",
                "a": {"id": e["id"], "chain": None},
                "b": {"id": p["id"], "chain": None},
                "type": e["pair_type"],
                "known": {
                    "ligand_touching": [e["ligand_touching"], p["ligand_touching"]],
                    "assembly1_chains": [e["assembly1_chains"], p["assembly1_chains"]],
                    "ca_rmsd": e["ca_rmsd"],
                },
            }
        )
    return out


def diff_set(samples=DIFF_SAMPLES, hinges=HINGES) -> list[dict]:
    """Every sample's pairs, then the hinge set's (published domain motions and synthetic hinges)."""
    pairs = [p for s in samples for p in sample_pairs(s)]
    return pairs + json.loads(Path(hinges).read_text(encoding="utf-8"))["pairs"]


def diff_stem(pair: dict) -> str:
    def one(side: dict) -> str:
        return side.get("stem") or side["id"] + (f"_{side['chain']}" if side.get("chain") else "")

    return one(pair["a"]) + "-" + one(pair["b"])


def _hinge_file(pair: dict, out: Path) -> Path:
    """A synthetic hinge (pair["synthetic"]): A's chain with residues `moved` turned `turn` degrees about `about`'s
    CA (axis `axis`), written as PDB."""
    import gemmi
    import numpy as np

    s = pair["synthetic"]
    st = gemmi.read_structure(str(ROOT / pair["a"]["id"]))
    chain = st[0][pair["a"]["chain"]]
    lo, hi = s["moved"]
    centre = np.array(chain[str(s["about"])][0]["CA"][0].pos.tolist())
    axis = np.asarray(s["axis"], float) / np.linalg.norm(s["axis"])
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    angle = np.radians(s["turn"])
    turn = np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k
    for res in chain:
        if lo <= res.seqid.num <= hi:
            for atom in res:
                atom.pos = gemmi.Position(*(turn @ (np.array(atom.pos.tolist()) - centre) + centre))
    path = out / f"{diff_stem(pair)}.pdb"
    st.write_pdb(str(path))
    return path


def _measure(sa, sb, domains: dict | None) -> dict:
    """The Diff's numbers at each cut of CUTS (diff.CORE_CUT set for the call): the core kept (pairs, and its share
    in each domain), the core's RMSD, the largest displacement, residues moved past MOVED_A, the domains' turns after
    the core fit and, core-free, each two domains' relative rotation (Kabsch of each on its own)."""
    import numpy as np

    from foldmap import diff as diff_module
    from foldmap.symmetry import axis_angle, kabsch

    out: dict = {}
    keep = diff_module.CORE_CUT
    try:
        for cut in CUTS:
            diff_module.CORE_CUT = cut
            d = diff_module.diff_sides(sa, sb, "auto", domains)
            aln = d.alignment
            ia, ib = (np.array(k, dtype=int) for k in zip(*aln.pairs))
            pa, pb = np.asarray(sa.bb.ca, float)[ia], np.asarray(sb.bb.ca, float)[ib]
            core = set(d.core)
            sel = np.array([p in core for p in aln.pairs], bool)
            moved = pb[sel] @ d.rotation.T + d.translation - pa[sel]
            row = {
                "core": len(d.core),
                "aligned": len(aln.pairs),
                "core_rmsd": _round(np.sqrt((moved**2).sum(axis=1).mean())),
                "max_displacement": _round(max(d.displacement.values(), default=0.0)),
                "moved": sum(v > MOVED_A for v in d.displacement.values()),
                "domain_turns": {k: _round(v) for k, v in d.domain_turns.items()},
            }
            if domains:
                mine = {k: diff_module._domain_residues(k, v, sa) for k, v in domains.items()}
                row["core_share"] = {k: sum(i in m for i, _ in d.core) for k, m in mine.items()}
                rot = {}
                for k, m in mine.items():
                    s = np.array([i in m for i in ia], bool)
                    rot[k] = kabsch(pb[s], pa[s])[0]
                names = list(domains)
                row["relative"] = {
                    f"{x}/{y}": _round(np.degrees(axis_angle(rot[x] @ rot[y].T)[0]))
                    for n, x in enumerate(names)
                    for y in names[n + 1 :]
                }
            out[f"{cut:g}"] = row
    finally:
        diff_module.CORE_CUT = keep
    out["route"] = d.alignment.method
    out["all_rmsd"] = _round(kabsch(pb, pa)[2])  # one fit on every aligned pair, as random-100's ca_rmsd
    out["counts"] = dict(sorted(diff_module.counts(d).items()))
    out["ligands"] = [list(r) for r in d.contacts["ligands"]]
    out["interfaces"] = [list(r) for r in d.contacts["interfaces"]]
    return out


def _known(pair: dict, m: dict) -> dict | None:
    """The diff against what the sample records of the pair: apo/holo, a ligand gained (or lost) on the side the
    sample says touches one; monomer/oligomer, an interface gained (or lost) toward the oligomer; open/closed, the
    RMSD of one fit on every aligned pair against the sample's ca_rmsd, and residues moved."""
    known = pair.get("known")
    if not known:
        return None
    t = pair["type"]
    if t == "apo/holo":
        want = "gained" if known["ligand_touching"][1] else "lost"
        got = [r for r in m["ligands"] if r[-1] == want]
        return {"expect": f"a ligand {want}", "agrees": bool(got), "got": m["ligands"]}
    if t == "monomer/oligomer":
        want = "gained" if known["assembly1_chains"][1] > known["assembly1_chains"][0] else "lost"
        got = [r for r in m["interfaces"] if r[-1] == want]
        return {"expect": f"an interface {want}", "agrees": bool(got), "got": m["interfaces"]}
    rmsd = m["all_rmsd"]
    return {
        "expect": f"CA RMSD {known['ca_rmsd']} Å",
        "agrees": rmsd is not None and abs(rmsd - known["ca_rmsd"]) <= 0.5,
        "got": rmsd,
    }


def diff_pair(pair: dict, out: Path) -> dict:
    """One pair: drawn as `foldmap diff A B` draws it (compare.diff_page, the publication look, auto route and
    assembly, the pair's domains), the page written (DIR/<stem>.png and .svg, the JSON as --json writes it), H5 and
    the S-/A-checks of each panel (compare.panel_checks) and the Diff's numbers at each cut of CUTS (_measure)."""
    from foldmap.compare import diff_page, panel_checks, write_json
    from foldmap.diff import CORE_CUT, read_side
    from foldmap.render import save
    from foldmap.style import resolve_style

    s = diff_stem(pair)
    a = pair["a"]["id"] if "synthetic" not in pair else str(ROOT / pair["a"]["id"])
    b = pair["b"]["id"] if "synthetic" not in pair else str(_hinge_file(pair, out))
    ca, cb = pair["a"].get("chain"), pair["b"].get("chain")
    domains = pair.get("domains") or None
    look = resolve_style("publication")
    t0 = time.perf_counter()
    made = diff_page(a, b, ca, cb, "auto", look, "auto", domains)
    seconds = time.perf_counter() - t0
    for ext in ("png", "svg"):
        save(made.figure, out / f"{s}.{ext}")
    write_json(made.data, out / f"{s}-diff.json")
    found: dict = {}
    h5 = panel_checks(made, "auto", found)
    sa = read_side(a, ca, "A", "auto", look.helices_310)
    sb = read_side(b, cb, "B", "auto", look.helices_310)
    measured = _measure(sa, sb, domains)
    result = {
        "pair": s,
        "name": pair["name"],
        "type": pair["type"],
        "chains": [sa.chain, sb.chain],
        "drawn_at_cut": CORE_CUT,
        "note": made.note,
        "h5": {k: {"pass": v["pass"], "detail": v.get("detail", "")} for k, v in h5.items()},
        "checks": {k: {"panel": _summary(x), "alone": _summary(y)} for k, (x, y) in found.items()},
        "measured": measured,
        "known": _known(pair, measured),
        "published": pair.get("published"),
        "seconds": round(seconds, 1),
    }
    result["source_stamp"] = _stamp()
    (out / f"{s}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return result


def cmd_diff_pair(args) -> int:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    r = diff_pair(json.loads(Path(args.pair).read_text(encoding="utf-8")), out)
    print(f"{r['pair']}: H5 " + ", ".join(f"{k} {_word(v)}" for k, v in r["h5"].items()))
    return 0


def cmd_diff(args) -> int:
    """Every pair of diff_set (or --pairs, by stem) in a child each, as run does; then the table (write_diff_report)."""
    pairs = diff_set()
    if args.pairs:
        wanted = set(args.pairs.split(","))
        pairs = [p for p in pairs if diff_stem(p) in wanted]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1"}

    def one(p: dict) -> tuple[str, int]:
        s = diff_stem(p)
        spec = out / f"{s}.pair.json"
        spec.write_text(json.dumps(p, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        cmd = [sys.executable, "-I", str(SCRIPT), "diff-pair", str(spec), "--out", str(out)]
        with open(out / f"{s}.log", "w", encoding="utf-8") as log:
            try:
                done = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, timeout=PAIR_TIMEOUT_S)
                code = done.returncode
            except subprocess.TimeoutExpired:
                code = -9
        print(f"{s}: exit {code}", flush=True)
        return s, code

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        codes = dict(pool.map(one, pairs))
    write_diff_report(diff_set(), out, Path(args.md) if args.md else out / "report.md")
    return 0 if all(c == 0 for c in codes.values()) else 1


def _turns(row: dict) -> str:
    return ", ".join(f"{k} {v:.1f}" for k, v in row.get("domain_turns", {}).items())


def write_diff_report(pairs: list[dict], out: Path, md: Path) -> Path:
    """The tables of a diff run: per pair the route, H5 per panel, the S-/A-checks each panel fails beside its chain
    drawn alone, the changes counted, the known difference; then each pair's numbers at both cuts."""
    from foldmap.source_stamp import code_mix, stamp_line

    measured = []  # each drawn pair's stamp (diff_pair)
    head = [
        "| pair | set | type | route | H5 A | H5 B | A fails (alone) | B fails (alone) | changes | known | agrees |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    cuts = [
        "| pair | cut | core / aligned | core RMSD | max disp. | moved > 3 Å | domain turns | relative "
        "| core per domain |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for p in pairs:
        s = diff_stem(p)
        f = out / f"{s}.json"
        if not f.is_file():
            head.append(f"| {s} | {p['name']} | {p['type']} | not drawn (see {s}.log) | | | | | | | |")
            continue
        r = json.loads(f.read_text(encoding="utf-8"))
        m, k = r["measured"], r["known"] or {}
        measured.append(r.get("source_stamp"))

        def fails(side: str, r=r) -> str:
            c = r["checks"].get(side, {})
            return f"{', '.join(c['panel']['failed']) or '-'} ({', '.join(c['alone']['failed']) or '-'})" if c else ""

        changed = ", ".join(f"{n.removeprefix('contacts.')} {v}" for n, v in m["counts"].items() if v) or "none"
        agrees = "" if not k else ("yes" if k["agrees"] else "no")
        head.append(
            f"| {s} | {r['name']} | {r['type']} | {m['route']} | {_word(r['h5'].get('a'))} | "
            f"{_word(r['h5'].get('b'))} | {fails('a')} | {fails('b')} | {changed} | {k.get('expect', '')} | {agrees} |"
        )
        for cut in CUTS:
            row = m[f"{cut:g}"]
            rel = ", ".join(f"{n} {v:.1f}" for n, v in row.get("relative", {}).items())
            share = ", ".join(f"{n} {v}" for n, v in row.get("core_share", {}).items())
            cuts.append(
                f"| {s} | {cut:g}·d0 | {row['core']} / {row['aligned']} | {row['core_rmsd']:.2f} | "
                f"{row['max_displacement']:.1f} | {row['moved']} | {_turns(row)} | {rel} | {share} |"
            )
    mix = code_mix(measured)
    code = [
        f"Report written by: {stamp_line(_stamp())}",
        f"Pairs measured by {len(mix['entry_code'])} code version(s): "
        + ", ".join(f"{c} ({n})" for c, n in mix["entry_code"].items())
        + (" -- MIXED CODE" if mix["mixed_code"] else ""),
        "",
    ]
    md.parent.mkdir(parents=True, exist_ok=True)
    md.write_text("\n".join(code + head) + "\n\n" + "\n".join(cuts) + "\n", encoding="utf-8")
    return md


def cmd_diff_report(args) -> int:
    out = Path(args.out)
    print(f"wrote {write_diff_report(diff_set(), out, Path(args.md) if args.md else out / 'report.md')}")
    return 0


# ------------------------------------------------------------------------------------------------------- sample


@dataclass(frozen=True)
class Domain:
    """One line of CATH's domain list (CathDomainList v2.0)."""

    name: str  # 1oaiA00: PDB id, chain, domain number (00: the whole chain is one domain)
    cath: tuple[int, int, int, int]  # class, architecture, topology, homologous superfamily
    s35: int
    s95: int
    length: int
    resolution: float  # 999 NMR, 1000 obsolete entry

    @property
    def pdb(self) -> str:
        return self.name[:4].upper()

    @property
    def chain(self) -> str:
        return self.name[4]

    def node(self, level: int) -> str:
        """The CATH code to `level` (2 architecture, 3 topology, 4 superfamily)."""
        return ".".join(str(x) for x in self.cath[:level])


def read_domain_list(text: str) -> list[Domain]:
    """CATH's domain list, every domain, in file order."""
    out = []
    for line in text.splitlines():
        f = line.split()
        if not f or line.startswith("#"):
            continue
        out.append(Domain(f[0], tuple(int(x) for x in f[1:5]), int(f[5]), int(f[7]), int(f[10]), float(f[11])))
    return out


def read_names(text: str) -> dict[str, str]:
    """CATH's node names (cath-names: node, representative domain, ':' name)."""
    out = {}
    for line in text.splitlines():
        if line.startswith("#") or ":" not in line:
            continue
        head, name = line.split(":", 1)
        out[head.split()[0]] = name.strip()
    return out


def eligible(domains: list[Domain], exclude: set[str]) -> list[Domain]:
    """The domains a held-out chain may be (module doc), one per entry, superfamily and S95 cluster (its first chain
    by name), sorted by name."""
    lo, hi = CHAIN_LENGTH
    kept: dict[tuple, Domain] = {}
    for d in sorted(domains, key=lambda d: d.name):
        if d.name[5:7] != "00" or d.chain == "0" or not lo <= d.length <= hi or d.resolution >= 1000:
            continue
        if d.pdb in exclude:
            continue
        kept.setdefault((d.pdb, d.cath, d.s95), d)
    return sorted(kept.values(), key=lambda d: d.name)


def _tree(domains: list[Domain], keys) -> dict:
    """Nested dicts by each key function in turn; the leaves are lists of domains in name order."""
    root: dict = {}
    for d in domains:
        node = root
        for k in keys[:-1]:
            node = node.setdefault(k(d), {})
        node.setdefault(keys[-1](d), []).append(d)
    return root


def _choose(rng: random.Random, options):
    return options[rng.randrange(len(options))]


def _down(rng: random.Random, node) -> Domain:
    """One domain under a tree node, each level chosen uniformly among its sorted keys."""
    while isinstance(node, dict):
        node = node[_choose(rng, sorted(node))]
    return _choose(rng, node)


class Draws:
    """The candidate pairs of each stratum, drawn from one seeded generator (module doc)."""

    def __init__(self, domains: list[Domain], seed: int):
        self.rng = random.Random(seed)
        topology = _tree(domains, [lambda d: d.node(3), lambda d: (d.node(4), d.s35)])
        self.remote = {t: c for t, c in topology.items() if len(c) >= 2}
        close = _tree(domains, [lambda d: d.node(4), lambda d: d.s35, lambda d: d.s95])
        self.close = {
            h: {s: c for s, c in sub.items() if len(c) >= 2}
            for h, sub in close.items()
            if any(len(c) >= 2 for c in sub.values())
        }
        arch = _tree(domains, [lambda d: d.node(2), lambda d: d.node(3), lambda d: d.node(4), lambda d: d.s35])
        self.negative = {a: t for a, t in arch.items() if len(t) >= 2}

    def draw(self, stratum: str) -> tuple[Domain, Domain]:
        rng = self.rng
        if stratum == "remote":
            clusters = self.remote[_choose(rng, sorted(self.remote))]
            one, two = rng.sample(sorted(clusters), 2)
            return _choose(rng, clusters[one]), _choose(rng, clusters[two])
        if stratum == "close":
            family = self.close[_choose(rng, sorted(self.close))]
            clusters = family[_choose(rng, sorted(family))]
            one, two = rng.sample(sorted(clusters), 2)
            return _choose(rng, clusters[one]), _choose(rng, clusters[two])
        topologies = self.negative[_choose(rng, sorted(self.negative))]
        one, two = rng.sample(sorted(topologies), 2)
        return _down(rng, topologies[one]), _down(rng, topologies[two])


def _chain_read(pdb: str, chain: str):
    """The chain's backbone as the campaign reads it (--assembly asu), and its residue count."""
    from foldmap.cli import resolve_structure
    from foldmap.io import load_backbone, read_model

    path = resolve_structure(pdb)
    bb = load_backbone(path, "asu", model=read_model(path, "asu"))
    return bb, sum(1 for label in bb.labels if label.chain == chain)


def verify(stratum: str, a: Domain, b: Domain) -> tuple[str, dict]:
    """Why a drawn pair is not kept ('' when it is), and what was measured: each chain's residues as read, and the
    structure alignment's identity, TM-score and aligned pairs."""
    from foldmap.align import align_chains

    lo, hi = CHAIN_LENGTH
    got: dict = {}
    sides = []
    for side, d in (("a", a), ("b", b)):
        try:
            bb, n = _chain_read(d.pdb, d.chain)
        except Exception as err:  # noqa: BLE001 - any failure to fetch or read rejects the pair, with its reason
            return f"{d.pdb} not read: {type(err).__name__}: {err}", got
        got[f"residues_{side}"] = n
        if not lo <= n <= hi:
            return f"{d.pdb} chain {d.chain}: {n} residues read", got
        sides.append(bb)
    aln = align_chains(sides[0], a.chain, sides[1], b.chain, "structure")
    got.update(identity=_round(aln.identity), tm_score=_round(aln.tm_score), aligned=len(aln.pairs))
    if stratum == "remote" and a.node(4) == b.node(4) and aln.identity >= REMOTE_IDENTITY:
        return f"one superfamily at identity {aln.identity:.3f} (remote: under {REMOTE_IDENTITY:g})", got
    if stratum == "close" and not CLOSE_IDENTITY[0] <= aln.identity <= CLOSE_IDENTITY[1]:
        return f"identity {aln.identity:.3f} (close: {CLOSE_IDENTITY[0]:g}-{CLOSE_IDENTITY[1]:g})", got
    return "", got


def _cath_text(name: str) -> str:
    """A CATH release file, downloaded once into ~/.cache/foldmap/cath."""
    from foldmap.fetch import AGENT, CACHE, atomic_write

    path = CACHE / "cath" / f"{name}-{CATH_RELEASE}.txt"
    if not path.is_file():
        req = urllib.request.Request(CATH_URL.format(release=CATH_RELEASE, name=name), headers={"User-Agent": AGENT})
        with urllib.request.urlopen(req, timeout=300) as r:
            atomic_write(path, r.read().decode("utf-8"))
    return path.read_text(encoding="utf-8")


def _describe(stratum: str, a: Domain, b: Domain, names: dict[str, str], got: dict) -> str:
    def named(code: str) -> str:
        return f"{code} ({names[code]})" if names.get(code) else code

    if stratum == "negative":
        head = f"negative: architecture {named(a.node(2))}; topologies {named(a.node(3))} / {named(b.node(3))}"
    elif stratum == "close":
        head = f"close: superfamily {named(a.node(4))}, S35 cluster {a.s35}, S95 clusters {a.s95} / {b.s95}"
    else:
        sf = named(a.node(4)) if a.node(4) == b.node(4) else f"{named(a.node(4))} / {named(b.node(4))}"
        head = f"remote: topology {named(a.node(3))}; superfamilies {sf}"
    return (
        f"{head}. CATH {CATH_RELEASE} domains {a.name} ({a.length} residues) and {b.name} ({b.length}); read "
        f"{got['residues_a']} and {got['residues_b']} residues; structure alignment: identity {got['identity']:.3f}, "
        f"TM-score {got['tm_score']:.3f} over {got['aligned']} pairs."
    )


def sample(domains: list[Domain], names: dict[str, str], seed: int, exclude: set[str], check=verify) -> dict:
    """The held-out set (module doc): {about, seed, cath_release, pairs (sorted as the tuning set), rejected (in
    draw order)}. `check(stratum, a, b)` -> (reason, measured) is verify; tests pass their own."""
    pool = eligible(domains, exclude)
    draws = Draws(pool, seed)
    used = set(exclude)
    pairs, rejected = [], []
    for stratum, wanted, kind in STRATA:
        kept = 0
        for _ in range(MAX_DRAWS):
            if kept == wanted:
                break
            a, b = draws.draw(stratum)
            if a.pdb == b.pdb or a.pdb in used or b.pdb in used:
                reason, got = "an entry already in the set (or both chains of one entry)", {}
            else:
                reason, got = check(stratum, a, b)
            if reason:
                rejected.append({"stratum": stratum, "a": a.name, "b": b.name, "reason": reason, **got})
                continue
            used |= {a.pdb, b.pdb}
            kept += 1
            pairs.append(
                {
                    "a": {"id": a.pdb, "chain": a.chain},
                    "b": {"id": b.pdb, "chain": b.chain},
                    "kind": kind,
                    "stratum": stratum,
                    "cath": {"a": f"{a.name} {a.node(4)}", "b": f"{b.name} {b.node(4)}"},
                    "note": _describe(stratum, a, b, names, got),
                }
            )
        if kept < wanted:
            raise RuntimeError(f"{stratum}: {kept} of {wanted} pairs in {MAX_DRAWS} draws")
    pairs.sort(key=lambda p: (p["a"]["id"], p["a"]["chain"], p["b"]["id"], p["b"]["chain"]))
    lo, hi = CHAIN_LENGTH
    return {
        "about": (
            f"Plan 13a held-out homologue set (spec 4.4), drawn by `scripts/homologues.py sample --seed {seed}` from "
            f"CATH {CATH_RELEASE} (CC BY 4.0, https://www.cathdb.info): whole-chain domains of {lo}-{hi} residues, "
            "entries of the tuning set excluded, each entry at most once. A is the reference, B is drawn like it. "
            "Strata: remote (one topology; other superfamily, or identity under 30 %), close (one superfamily, "
            "identity 30-90 %), negative (one architecture, other topology). Identity: identical / aligned pairs of "
            "foldmap's structure alignment. Reproducible given the same fetches: a draw whose entry cannot be read "
            "is rejected (listed under rejected, with the reason) and the next one drawn, so a fetch that fails on "
            "another day gives another set from the same seed. Measured once (Task 13a-11)."
        ),
        "seed": seed,
        "cath_release": CATH_RELEASE,
        "pairs": pairs,
        "rejected": rejected,
    }


def cmd_sample(args) -> int:
    exclude = {p[side]["id"].upper() for p in read_set(args.exclude) for side in ("a", "b")}
    data = sample(
        read_domain_list(_cath_text("cath-domain-list")), read_names(_cath_text("cath-names")), args.seed, exclude
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    counts = {s: sum(p["stratum"] == s for p in data["pairs"]) for s, _, _ in STRATA}
    print(f"wrote {out}: {counts}, {len(data['rejected'])} drawn pairs rejected")
    return 0


# ------------------------------------------------------------------------------------------------------- report


def _word(c: dict | None) -> str:
    if not c or c.get("pass") is None:
        return "–"
    return "pass" if c["pass"] else "**fail**"


def _place(v: dict | None) -> str:
    return "–" if not v else f"{v['median']:.2f} / {v['p90']:.2f}"


def _share(v: dict | None) -> str:
    return "–" if not v else f"{v['share']:.3f} ({v['pairs']})"


def _outcome(c: dict | None) -> str:
    if not c or c.get("pass") is None:
        return "not judged"
    return "pass" if c["pass"] else "fail"


def _median(xs: list[float]) -> float | None:
    return _round(statistics.median(xs)) if xs else None


def summarize(pairs: list[dict], results: dict[str, dict]) -> dict:
    """The held-out summary (module doc): per stratum the pairs, their results and each check's outcomes; the
    positives that pass H1-H5 and the negatives refused (H1); H2's median and H3's share over the positives with H2
    judged, drawn independently and drawn alike (the medians of the pairs' values); the alignment routes."""
    out: dict = {"strata": {}, "positives": {}, "negatives": {}}
    for stratum in dict.fromkeys(p.get("stratum", p["kind"]) for p in pairs):
        group = [p for p in pairs if p.get("stratum", p["kind"]) == stratum]
        got = [results[stem(p)] for p in group if stem(p) in results]
        out["strata"][stratum] = {
            "pairs": len(group),
            "results": len(got),
            "checks": {
                k: {o: sum(_outcome(r["like"][k]) == o for r in got) for o in ("pass", "fail", "not judged")}
                for k in CHECKS
            },
            "all_pass": sum(all(r["like"][k].get("pass") is True for k in CHECKS) for r in got),
            "routes": {
                m: sum((r["alignment"] or {}).get("method") == m for r in got) for m in ("sequence", "structure")
            },
            "tm_range": [
                _round(min(t)) if (t := [r["alignment"]["tm_score"] for r in got if r["alignment"]]) else None,
                _round(max(t)) if t else None,
            ],
        }
    pos = [p for p in pairs if p["kind"] == "positive"]
    neg = [p for p in pairs if p["kind"] == "negative"]
    rp = [results[stem(p)] for p in pos if stem(p) in results]
    rn = [results[stem(p)] for p in neg if stem(p) in results]
    judged = [r for r in rp if (r["like"]["H2"].get("value") and r["like"]["H2"].get("baseline"))]
    out["positives"] = {
        "pairs": len(pos),
        "results": len(rp),
        "all_pass": sum(all(r["like"][k].get("pass") is True for k in CHECKS) for r in rp),
        "pass": {k: sum(r["like"][k].get("pass") is True for r in rp) for k in CHECKS},
        "judged_h2": len(judged),
        "h2_median_independent": _median([r["like"]["H2"]["baseline"]["median"] for r in judged]),
        "h2_median_alike": _median([r["like"]["H2"]["value"]["median"] for r in judged]),
        "h2_p90_independent": _median([r["like"]["H2"]["baseline"]["p90"] for r in judged]),
        "h2_p90_alike": _median([r["like"]["H2"]["value"]["p90"] for r in judged]),
        "h3_share_independent": _median(
            [r["like"]["H3"]["baseline"]["share"] for r in judged if r["like"]["H3"].get("baseline")]
        ),
        "h3_share_alike": _median([r["like"]["H3"]["value"]["share"] for r in judged if r["like"]["H3"].get("value")]),
        "h2_closer_alike": sum(
            r["like"]["H2"]["value"]["median"] < r["like"]["H2"]["baseline"]["median"] for r in judged
        ),
    }
    out["negatives"] = {
        "pairs": len(neg),
        "results": len(rn),
        "refused": sum(r["like"]["H1"].get("pass") is True for r in rn),
        "h5_pass": sum(r["like"]["H5"].get("pass") is True for r in rn),
    }
    return out


def write_report(pairs: list[dict], out: Path, md: Path) -> Path:
    """The table of every pair's result (module doc)."""
    from foldmap.source_stamp import code_mix

    rows, failures, missing, results = [], {k: [] for k in CHECKS}, [], {}
    strata = any("stratum" in p for p in pairs)
    for p in pairs:
        f = out / f"{stem(p)}.json"
        if not f.is_file():
            missing.append(stem(p))
            continue
        r = json.loads(f.read_text(encoding="utf-8"))
        results[stem(p)] = r
        h, aln, m = r["like"], r["alignment"] or {}, r["mates"] or {}
        for k in CHECKS:
            if h[k].get("pass") is False:
                failures[k].append(f"{r['pair']} ({p.get('stratum', r['kind'])}): {h[k]['detail']}")
        lost = (h["H5"].get("value") or {}).get("lost") or []
        kind = f"{r['kind']} | {p['stratum']}" if strata else r["kind"]
        rows.append(
            f"| {r['pair']} | {kind} | {aln.get('method', '–')} | {aln.get('tm_score', 0):.3f} | "
            f"{m.get('mates', 0)}/{m.get('pieces', 0)}/{m.get('extras', 0)}/{m.get('only_a', 0)} | "
            + " | ".join(_word(h[k]) for k in CHECKS)
            + f" | {_place(h['H2'].get('value'))} | {_place(h['H2'].get('baseline'))} | "
            f"{_share(h['H3'].get('value'))} | {_share(h['H3'].get('baseline'))} | {', '.join(lost) or '–'} |"
        )
    lines = [
        "# Homologue pairs: drawn like the reference (H1-H5)",
        "",
        "TM: foldmap's TM-score normalised by A. Elements: mates / pieces / extras (B only) / only in A. H2: median / "
        "90th percentile of the mates' centre deviations (page units) after the best rigid 2D motion; H3: share of "
        "mate pairs whose page relation agrees (pairs judged). 'alike': B drawn like A; 'indep.': A and B drawn "
        "independently (the baseline). Lost: checks B drawn alone passes and B drawn like A does not (H5).",
        "",
        "| pair | kind | " + ("stratum | " if strata else "") + "route | TM | elements | H1 | H2 | H3 | H4 | H5 | "
        "H2 alike | H2 indep. | H3 alike | H3 indep. | lost |",
        "|---|---|" + ("---|" if strata else "") + "---|---|---|---|---|---|---|---|---|---|---|---|---|",
        *rows,
        "",
        "## Failures by check",
        "",
    ]
    for k in CHECKS:
        lines.append(f"- **{k}**: " + ("none" if not failures[k] else ""))
        lines += [f"  - {x}" for x in failures[k]]
    if missing:
        lines += ["", f"No result (the pair's child failed; see its log): {', '.join(missing)}"]
    if strata:
        summary = summarize(pairs, results)
        summary["source_stamp"] = _stamp()
        summary |= code_mix(r.get("source_stamp") for r in results.values())  # each pair's measuring code
        (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        lines += ["", "## Summary", "", "| stratum | pairs | results | all pass | " + " | ".join(CHECKS) + " |"]
        lines.append("|---|---|---|---|" + "---|" * len(CHECKS))
        for name, g in summary["strata"].items():
            cells = [f"{g['checks'][k]['pass']} / {g['checks'][k]['fail']} / {g['checks'][k]['not judged']}"
                     for k in CHECKS]  # fmt: skip
            lines.append(f"| {name} | {g['pairs']} | {g['results']} | {g['all_pass']} | " + " | ".join(cells) + " |")
        ps, ns = summary["positives"], summary["negatives"]

        def num(x):
            return "–" if x is None else f"{x:.2f}"

        lines += [
            "",
            "Cells: pass / fail / not judged.",
            "",
            f"- Positives passing H1-H5: {ps['all_pass']} of {ps['pairs']}; negatives refused (H1): {ns['refused']} "
            f"of {ns['pairs']}.",
            f"- Over the {ps['judged_h2']} positives with H2 judged (medians of the pairs' values): H2 median "
            f"drawn independently {num(ps['h2_median_independent'])}, drawn alike {num(ps['h2_median_alike'])}; "
            f"H2 90th percentile {num(ps['h2_p90_independent'])} and {num(ps['h2_p90_alike'])}; H3 share "
            f"{num(ps['h3_share_independent'])} and {num(ps['h3_share_alike'])}. Closer drawn alike (H2 median): "
            f"{ps['h2_closer_alike']} of {ps['judged_h2']}.",
        ]
    md.parent.mkdir(parents=True, exist_ok=True)
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    return md


def cmd_report(args) -> int:
    out = Path(args.out)
    write_report(read_set(args.set), out, Path(args.md) if args.md else out / "report.md")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    al = sub.add_parser("align", help="align and mate every pair, no drawing")
    al.add_argument("--set", required=True)
    al.add_argument("--out", required=True)
    al.add_argument("--usalign", default="", help="also compare with US-align: its binary, or 'auto' (find_usalign)")
    al.add_argument("--jobs", type=int, default=1)
    ur = sub.add_parser("usalign-report", help="the US-align comparison's files from align --usalign runs")
    ur.add_argument("runs", nargs="+", help="NAME=DIR of each align --usalign run (DIR/align.json)")
    ur.add_argument("--out", required=True, help="where usalign.json and usalign-report.md go")
    ur.add_argument("--summary", default="", help="a summary JSON (heldout-summary.json) whose 'usalign' is set")
    ur.add_argument("--commit", required=True, help="the commit the runs measured")
    pr = sub.add_parser("pair", help="one pair (run's child)")
    pr.add_argument("a", help="ID:CHAIN of the reference")
    pr.add_argument("b", help="ID:CHAIN drawn like it")
    pr.add_argument("--kind", choices=("positive", "negative"), required=True)
    pr.add_argument("--out", required=True)
    rn = sub.add_parser("run", help="every pair of a set, then the report")
    rn.add_argument("--set", required=True)
    rn.add_argument("--out", required=True)
    rn.add_argument("--jobs", type=int, default=4)
    rn.add_argument("--pairs", default="", help="only these stems (A_CHAIN-B_CHAIN), comma-separated")
    rn.add_argument("--md", default="")
    sm = sub.add_parser("sample", help="draw the held-out set from CATH (needs network)")
    sm.add_argument("--seed", type=int, required=True)
    sm.add_argument("--out", required=True)
    sm.add_argument("--exclude", default=str(TUNING), help="a set whose entries are left out (default: the tuning set)")
    rp = sub.add_parser("report", help="the table of a run")
    rp.add_argument("--set", required=True)
    rp.add_argument("--out", required=True)
    rp.add_argument("--md", default="")
    df = sub.add_parser("diff", help="the two-state pairs drawn as diffs (plan 13e-3), then the tables")
    df.add_argument("--out", required=True)
    df.add_argument("--jobs", type=int, default=4)
    df.add_argument("--pairs", default="", help="only these stems (A-B), comma-separated")
    df.add_argument("--md", default="")
    dp = sub.add_parser("diff-pair", help="one diff pair (diff's child)")
    dp.add_argument("pair", help="the pair as JSON (a diff_set item)")
    dp.add_argument("--out", required=True)
    dr = sub.add_parser("diff-report", help="the tables of a diff run")
    dr.add_argument("--out", required=True)
    dr.add_argument("--md", default="")
    args = ap.parse_args(argv)
    _stamp()  # the code as this process starts, written with every result it measures (review M8)
    commands = {"align": cmd_align, "usalign-report": cmd_usalign_report, "pair": cmd_pair, "run": cmd_run,
                "sample": cmd_sample, "report": cmd_report,
                "diff": cmd_diff, "diff-pair": cmd_diff_pair, "diff-report": cmd_diff_report}  # fmt: skip
    return commands[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
