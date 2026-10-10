# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Compare foldmap's DSSP with mkdssp on perturbed copies of the test structures, residue by residue.

Needs mkdssp 3.1.2 on PATH; CI does not run this. Variants are written as PDB, because mkdssp 3.1.2 fails on
mmCIF written by gemmi. Model 1 only. Modes:

    none       the structure as deposited
    delete     about one residue in 60 of each chain removed (chain breaks)
    jitter25   every atom moved by a Gaussian of 0.25 Å per axis (H-bond energies near the -0.5 cutoff)
    jitter50   the same with 0.5 Å
    deljitter  delete and jitter25 together
    icode      residues 32-45 of every chain longer than 45 renumbered 31A, 31B, ... (insertion codes)
    bulge      one residue inside a beta-bulge removed (one variant per residue a joined ladder skips)

    python scripts/fuzz_dssp.py [--seeds 3] [--modes jitter25,bulge] [1TIM 5NKT ...]

Prints every variant that disagrees, then the disagreements per mode; exits 1 if there are any.
"""

from __future__ import annotations

import argparse
import functools
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import gemmi
from make_dssp_reference import check_version, parse_dssp

from foldmap.dssp import _find_bridges, _hbond_set, _ladders, assign_dssp
from foldmap.io import load_backbone

DATA = Path(__file__).resolve().parents[1] / "tests" / "data"
MODES = ("none", "delete", "jitter25", "jitter50", "deljitter", "icode", "bulge")
SEEDLESS = {"none", "icode"}  # one variant whatever --seeds says (bulge: one per bulge residue)
DELETE_EVERY = 60  # residues per deletion
ICODE_FIRST, ICODE_COUNT = 31, 14  # residue indices renumbered with insertion codes


def mkdssp(path: Path) -> dict[tuple[str, int, str], str] | None:
    """mkdssp's code per (chain, number, insertion code), or None when it fails on the file."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.dssp"
        run = subprocess.run(["mkdssp", str(path), str(out)], capture_output=True, text=True)
        if run.returncode:
            return None
        return {(c, n, i): code for c, n, i, code in parse_dssp(out.read_text(encoding="utf-8"))}


def compare(path: Path) -> tuple[int, int, list] | None:
    """Residues, residues only one program saw, and (residue, mkdssp, foldmap) for each disagreement."""
    ref = mkdssp(path)
    if ref is None:
        return None
    bb = load_backbone(path, "asu")
    ss = assign_dssp(bb).ss
    keys = [(lab.chain, lab.seq, lab.icode) for lab in bb.labels]
    bad = [(k, ref[k], ss[i]) for i, k in enumerate(keys) if k in ref and ref[k] != ss[i]]
    return len(keys), len(set(keys) ^ set(ref)), bad


@functools.cache
def _bulge_residues(name: str) -> list[tuple[str, int, str]]:
    """Residues a joined ladder skips on one strand, in the deposited structure."""
    bb = load_backbone(DATA / f"{name}.cif", "asu")
    out = []
    for lad in _ladders(_find_bridges(bb, _hbond_set(bb)), bb.prev):
        for side in (lad.i, lad.j):
            out += [bb.labels[a + 1] for a, b in zip(side, side[1:]) if b - a > 1]
    return [(lab.chain, lab.seq, lab.icode) for lab in out]


def _jitter(model: gemmi.Model, rng: random.Random, sd: float) -> None:
    for chain in model:
        for res in chain:
            for atom in res:
                p = atom.pos
                atom.pos = gemmi.Position(p.x + rng.gauss(0, sd), p.y + rng.gauss(0, sd), p.z + rng.gauss(0, sd))


def variant(name: str, mode: str, seed: int, tmp: Path) -> Path:
    """A perturbed copy of a test structure as PDB (for bulge, seed picks the residue)."""
    rng = random.Random(seed)
    st = gemmi.read_structure(str(DATA / f"{name}.cif"))
    st.setup_entities()
    while len(st) > 1:
        del st[1]
    model = st[0]
    if mode in ("delete", "deljitter"):
        for chain in model:
            for _ in range(max(1, len(chain) // DELETE_EVERY)):
                del chain[rng.randrange(len(chain))]
    if mode in ("jitter25", "deljitter"):
        _jitter(model, rng, 0.25)
    if mode == "jitter50":
        _jitter(model, rng, 0.5)
    if mode == "icode":
        for chain in model:
            res = list(chain)
            if len(res) > ICODE_FIRST + ICODE_COUNT:
                base = res[ICODE_FIRST - 1].seqid.num
                for k, r in enumerate(res[ICODE_FIRST : ICODE_FIRST + ICODE_COUNT]):
                    r.seqid = gemmi.SeqId(base, chr(ord("A") + k))
    if mode == "bulge":
        chain_name, num, icode = _bulge_residues(name)[seed]
        chain = next(c for c in model if c.name == chain_name)
        del chain[next(k for k, r in enumerate(chain) if (r.seqid.num, r.seqid.icode.strip()) == (num, icode))]
    path = tmp / f"{name}_{mode}{seed}.pdb"
    st.write_pdb(str(path))
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("names", nargs="*", help="test structure IDs (default: all of tests/data)")
    ap.add_argument("--seeds", type=int, default=3, help="variants per structure for delete and jitter modes")
    ap.add_argument("--modes", default=",".join(MODES), help="comma-separated, from: " + " ".join(MODES))
    args = ap.parse_args()
    check_version()
    modes = args.modes.split(",")
    names = args.names or sorted(p.stem for p in DATA.glob("*.cif"))
    stats = {m: [0, 0, 0, 0, 0] for m in modes}  # variants, residues, disagreements, key mismatches, failures
    with tempfile.TemporaryDirectory() as tmp:
        for mode in modes:
            for name in names:
                count = 1 if mode in SEEDLESS else len(_bulge_residues(name)) if mode == "bulge" else args.seeds
                for seed in range(count):
                    row = stats[mode]
                    result = compare(variant(name, mode, seed, Path(tmp)))
                    if result is None:
                        row[4] += 1
                        print(name, mode, seed, "mkdssp failed")
                        continue
                    n, missing, bad = result
                    row[0] += 1
                    row[1] += n
                    row[2] += len(bad)
                    row[3] += missing
                    if bad or missing:
                        print(name, mode, seed, f"n={n} keys={missing} disagree={len(bad)}", bad[:8])
    print(f"\n{'mode':10} {'variants':>8} {'residues':>9} {'disagree':>8} {'keys':>5} {'failed':>6}")
    for mode, (v, n, bad, miss, failed) in stats.items():
        print(f"{mode:10} {v:8d} {n:9d} {bad:8d} {miss:5d} {failed:6d}")
    return int(any(row[2] or row[3] for row in stats.values()))


if __name__ == "__main__":
    sys.exit(main())
