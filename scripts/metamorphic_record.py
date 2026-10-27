# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The record of plan 28 Task 7 (review 3.6, 3.8; rulings Q34, Q35): what the metamorphic tests and the independent
kernels measure, the findings MF1-MF4 with their mechanisms, and the deep checks' external runs.

    python scripts/metamorphic_record.py [--out docs/validation/metamorphic]

Writes report.json and report.md, stamped with the code that measured (foldmap.source_stamp). It measures with the
tests' own helpers (tests/metamorphic.py, tests/test_metamorphic.py, tests/test_independent_kernels.py), on committed
fixtures only:

- exact transforms and the general rotation: the largest geometry gap of each fixture's figure, and whether its
  structure facts and `foldmap check` verdicts are equal;
- MF1: the reversed chain renaming of the assemblies of several deposited chains (1TIM, 1LMB, 1BL8);
- noise: the largest change over three seeds per fixture and sigma (DSSP codes as a fraction of residues; bridges,
  base pairs and transmembrane residues as |x △ y| / (|x| + |y|));
- the kernels: H-bonds and bridges per protein fixture; base pairs, both, only foldmap, only the kernel, and every
  disagreement by why (MF2 hbond, MF3 c1, MF4 plane, lambda, modified);
- D2 (mkdssp, when MKDSSP is set) on the protein fixtures, D3 (RNApolis, when NA_RNAPOLIS_PYTHON is set) on the
  nucleic-acid fixtures as committed and as PDB copies, and which MF2 pairs RNApolis calls. A recorded run, not a test:
  the tests that check the RNApolis reading (tests/test_deep_check.py: test_d3_with_rnapolis,
  test_mf2_pairs_are_rnapolis_pairs) are skipped where RNApolis is not set up."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))

import gemmi  # noqa: E402
import numpy as np  # noqa: E402

from foldmap import deep  # noqa: E402
from foldmap.source_stamp import source_stamp, stamp_line  # noqa: E402

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts
OUT = ROOT / "docs" / "validation" / "metamorphic"
REVERSED = ("1TIM.cif", "1LMB.cif", "1BL8.cif")
D3_ENTRIES = ("1EHZ", "2OEU", "4OO8", "1KX5", "1AOI")
FINDINGS = {
    "MF1": "Chain names choose protomer 0: symmetry._orbits walks the chains in name order (symmetry.py:234), "
    "symmetry._protomers sorts the orbits by their first chain's name (symmetry.py:255) and breaks its ties by name "
    "(symmetry.py:259); protomer 0 is drawn on the left (cli._dyad_frame) or is the copy the others mirror "
    "(cli._mirror_frame); on an unrolled frame its nucleic-acid stems are the template its copies' stems are tied to "
    "(layout._stem_ties, plan 28's NA fix F2), and when protomer 0 holds no nucleic acid no stem is tied at all, "
    "so a renaming can change that figure in kind. A renaming that reverses the chains' name order moves the figure of "
    "an assembly of several deposited chains; its structure facts and check verdicts stay equal.",
    "MF2": "Base-pair H-bond window: foldmap calls a canonical pair whose Watson-Crick-edge H-bonds reach 3.8 Å "
    "(nucleic.HBOND_MAX, nucleic.py:16, applied at :147), one stretched to 4.0 Å (HBOND_STRETCH, :17, applied at "
    ":149); the kernel holds two of them to 3.5 Å.",
    "MF3": "C1'-C1' window: foldmap accepts 9.0-12.0 Å (nucleic.C1_PAIR, nucleic.py:18, applied at :127 and :131); "
    "the kernel 9.5-11.5 Å.",
    "MF4": "Planarity: foldmap drops a pair whose base planes are more than 45 degrees apart (nucleic.PLANE_MAX, "
    "nucleic.py:19, applied at :156); the kernel has no planarity bound.",
}


def _gap(a: dict, b: dict) -> float:
    import metamorphic as M

    return M.geometry_gap(a, b)


def transforms(work: Path) -> dict:
    """{fixture: {transform: {"gap", "facts_equal", "verdicts_equal"}}} for the exact transforms, the general
    rotation and (REVERSED) the order-reversing renaming."""
    import test_metamorphic as TM

    out = {}
    names = list(dict.fromkeys([*TM.FIXTURES, *REVERSED]))
    for name in names:
        base = TM._read(work, name, "base", True)
        rows = {}
        kinds = [*TM.EXACT, "turn"] if name in TM.FIXTURES else []
        for t in kinds + (["reverse"] if name in REVERSED else []):
            got = TM._read(work, name, t, True)
            rows[t] = {
                "gap": _gap(base["geometry"], got["geometry"]),
                "facts_equal": got["facts"] == base["facts"],
                "verdicts_equal": got["verdicts"] == base["verdicts"],
            }
        out[name] = rows
    return out


def noise(work: Path) -> dict:
    import metamorphic as M
    import test_metamorphic as TM

    out = {}
    for name, kw in TM.FIXTURES.items():
        base = M.facts(M.rewrite(M.DATA / name, work / "noise-base.cif"), **kw)
        for sigma in sorted(TM.NOISE):
            worst: dict = {}
            for seed in range(3):
                got = M.facts(M.rewrite(M.DATA / name, work / "noise.cif", noise=(sigma, seed)), **kw)
                rates = {"bridges": TM._rate(base["bridges"], got["bridges"])}
                if base["ss"]:
                    rates["ss"] = float(np.mean([a != b for a, b in zip(base["ss"], got["ss"])]))
                if "base_pairs" in base:
                    rates["base_pairs"] = TM._rate(base["base_pairs"], got["base_pairs"])
                if "transmembrane" in base:
                    rates["tm_residues"] = TM._rate(base["transmembrane_residues"], got["transmembrane_residues"])
                for k, v in rates.items():
                    worst[k] = round(max(worst.get(k, 0.0), v), 4)
            out.setdefault(name, {})[str(sigma)] = {
                **worst,
                "bridges_n": len(base["bridges"]),
                "base_pairs_n": len(base.get("base_pairs", [])),
            }
    return out


def kernels() -> dict:
    import test_independent_kernels as K
    from foldmap.dssp import _hbond_set, assign_dssp

    proteins = {}
    for name in K.PROTEIN_FIXTURES:
        bb = K._protein(name)
        ks, ours = K.ks_bonds(bb), _hbond_set(bb)
        bridges = K.ks_bridges(bb, ks)
        proteins[name] = {
            "residues": len(bb),
            "kernel_bonds": len(ks),
            "foldmap_bonds": len(ours),
            "kernel_only": len(ks - ours),
            "foldmap_only": len(ours - ks),
            "bridges": len(bridges),
            "bridges_equal": bridges == {(b.i, b.j, b.kind) for b in assign_dssp(bb).bridges},
        }
    counts = {"both": 0, "only_foldmap": 0, "only_kernel": 0}
    for name in K.NA_FIXTURES:
        ours, theirs = K.foldmap_pairs(K.NA / name).keys(), K.canonical_pairs(K.NA / name).keys()
        counts["both"] += len(ours & theirs)
        counts["only_foldmap"] += len(ours - theirs)
        counts["only_kernel"] += len(theirs - ours)
    why = {w: K._any(w) for w in ("hbond", "c1", "plane", "lambda", "modified")}
    return {"dssp": proteins, "base_pairs": {**counts, "disagreements": why}}


def external(work: Path) -> dict:
    """D2 on the protein fixtures; D3 on the nucleic-acid fixtures as committed and as PDB copies, with the MF2 pairs
    RNApolis does not call."""
    import test_independent_kernels as K

    d2 = {name: deep.dssp_check(K.DATA / f"{name}.cif") for name in K.PROTEIN_FIXTURES}
    d3 = {}
    for name in D3_ENTRIES:
        src = K.NA / f"{name}.cif.gz"
        pdb = work / f"{name}.pdb"
        gemmi.read_structure(str(src)).write_pdb(str(pdb))
        as_pdb = deep.pairs_check(pdb)
        row = {"mmcif": deep.pairs_check(src), "pdb": as_pdb}
        theirs = deep.rnapolis_pairs(pdb)
        if isinstance(theirs, set):
            mf2 = K.mf2_pairs(f"{name}.cif.gz")
            row["mf2_pairs"] = len(mf2)
            row["mf2_not_in_rnapolis"] = [p for p in mf2 if frozenset(map(_key, p.split("/"))) not in theirs]
        d3[name] = row
    return {"D2": d2, "D3": d3}


def _key(text: str):
    import re

    chain, rest = text.split(":")
    number, icode = re.fullmatch(r"(-?\d+)([A-Za-z]?)", rest).groups()
    return chain, int(number), icode


def render(doc: dict) -> str:
    lines = [
        "# Metamorphic tests and independent kernels",
        "",
        f"Plan 28 Task 7 (rulings Q34, Q35). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "What tests/test_metamorphic.py, tests/test_independent_kernels.py and tests/test_deep_check.py assert, "
        "measured on committed fixtures; the findings they pin with xfail(strict). The deep-check runs below are "
        "recorded runs (they need mkdssp and RNApolis), not tests.",
        "",
        "## Findings",
        "",
        "Named MF1-MF4 (metamorphic findings), apart from plan 28's NA fixes F1-F7 (docs/validation/na-fixes).",
        "",
    ]
    for code, text in FINDINGS.items():
        lines.append(f"- **{code}.** {text}")
    rev = {n: rows["reverse"] for n, rows in doc["transforms"].items() if "reverse" in rows}
    lines += ["", "MF1, the order-reversing renaming (largest geometry gap, page units):", ""]
    lines += ["| Fixture | gap | facts equal | verdicts equal |", "|---|---|---|---|"]
    lines += [f"| {n} | {r['gap']:.4g} | {r['facts_equal']} | {r['verdicts_equal']} |" for n, r in rev.items()]
    bp = doc["kernels"]["base_pairs"]
    why = bp["disagreements"]
    lines += [
        "",
        f"MF2-MF4, base pairs over {len(doc['na_fixtures'])} fixtures: both {bp['both']}, only foldmap "
        f"{bp['only_foldmap']}, only the kernel {bp['only_kernel']}. By why: MF2 (hbond) {len(why['hbond'])}, "
        f"MF3 (c1) {len(why['c1'])}, MF4 (plane) {len(why['plane'])}, glycosidic angle {len(why['lambda'])} (never "
        f"alone), modified nucleotides {len(why['modified'])} (outside the kernel's scope; not a finding).",
        "",
    ]
    for w in ("c1", "plane", "lambda", "modified"):
        lines.append(f"- {w}: {', '.join(why[w])}")
    lines.append(f"- hbond: {', '.join(why['hbond'])}")
    lines += ["", "## Transforms (largest geometry gap; facts and verdicts equal)", ""]
    cols = ["motion", "rename", "shuffle", "pdb", "turn"]
    lines += ["| Fixture | " + " | ".join(cols) + " |", "|---" * (len(cols) + 1) + "|"]
    for n, rows in doc["transforms"].items():
        if "motion" not in rows:
            continue
        cells = [
            f"{rows[c]['gap']:.2g}" + ("" if rows[c]["facts_equal"] and rows[c]["verdicts_equal"] else " (differ)")
            for c in cols
        ]
        lines.append(f"| {n} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "1C3W's pdb gap comes from its assembly: the PDB file's BIOMT operators are written to 6 decimals "
        "(tests/test_metamorphic.py: test_pdb_assembly_operators_are_rounded).",
        "",
        "## Noise (largest change over three seeds)",
        "",
        "| Fixture | sigma Å | DSSP codes | bridges (n) | base pairs (n) | TM residues |",
        "|---|---|---|---|---|---|",
    ]
    for n, by in doc["noise"].items():
        for sigma, r in by.items():
            lines.append(
                f"| {n} | {sigma} | {r.get('ss', '-')} | {r['bridges']} ({r['bridges_n']}) | "
                f"{r.get('base_pairs', '-')} ({r['base_pairs_n']}) | {r.get('tm_residues', '-')} |"
            )
    lines += ["", "## Kernels: backbone H-bonds and bridges (Kabsch & Sander 1983)", ""]
    lines += ["| Fixture | residues | kernel bonds | foldmap bonds | kernel only | foldmap only | bridges equal |"]
    lines += ["|---|---|---|---|---|---|---|"]
    for n, r in doc["kernels"]["dssp"].items():
        lines.append(
            f"| {n} | {r['residues']} | {r['kernel_bonds']} | {r['foldmap_bonds']} | {r['kernel_only']} | "
            f"{r['foldmap_only']} | {r['bridges_equal']} ({r['bridges']}) |"
        )
    lines += [
        "",
        "Kernel only: the peptide neighbour's CO(i-1)...NH(i), which the energy formula puts below -0.5 kcal/mol and "
        "the DSSP program (and foldmap, dssp.py:69) skips.",
        "",
        "## Deep checks, recorded runs",
        "",
    ]
    for n, r in doc["external"]["D2"].items():
        lines.append(f"- D2 {n}: {r['status']} ({r['detail']})")
    for n, r in doc["external"]["D3"].items():
        extra = ""
        if "mf2_pairs" in r:
            extra = f"; MF2 pairs {r['mf2_pairs']}, not in RNApolis: {', '.join(r['mf2_not_in_rnapolis']) or 'none'}"
        lines.append(
            f"- D3 {n}: as committed {r['mmcif']['status']} ({r['mmcif']['detail']}); as PDB {r['pdb']['status']} "
            f"({r['pdb']['detail']}){extra}"
        )
    lines += [
        "",
        "Earlier RNApolis runs on whole downloaded entries, from older code: docs/validation/na/report.md.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(OUT))
    args = parser.parse_args(argv)
    import test_independent_kernels as K

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        doc = {
            "source_stamp": SOURCE_STAMP,
            "findings": FINDINGS,
            "na_fixtures": K.NA_FIXTURES,
            "transforms": transforms(work),
            "noise": noise(work),
            "kernels": kernels(),
            "external": external(work),
        }
    doc["base_pairs"] = doc["kernels"]["base_pairs"]  # where tests/test_independent_kernels.py reads the record
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8", newline="\n")
    (out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(f"wrote {out / 'report.json'} and report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
