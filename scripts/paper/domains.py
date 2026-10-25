# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Step 2 (zyla-lab): score the automatic domains against CATH on the committed fixtures, with the code at the harness
commit (build/tree/src, exported by extract.py), never the installed foldmap.

    python scripts/paper/domains.py   (step 2 of scripts/paper/run_all.sh)

Writes build/domains.json: per chain of the 33-chain reference (cath_reference.json, inputs cath_inputs.npz) and the
11-chain fresh set (cath_fresh.json, fresh_inputs.npz) its overlap, boundary recall, precision and F1
(domain_score.inputs_score) and the overlap of calling the chain one domain (the baseline), and the splits:

  tuning / held   the split the parameters were chosen on: held = random.Random("20261016:domains").sample of the
                  sorted (pdb, chain) keys, a third (11), drawn chain by chain (domain_score.held_out at 6c6c089; the
                  module's own held_out now draws by homology group, which is after the fact)
  held_nohom      the held chains whose homology group ("group") has no member among the tuning chains
  group_held      domain_score.held_out as it is now (by homology group), reported for comparison only

Pure computation on committed arrays (no structure file is read); about a minute on one core.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import BUILD, TREE, run_stamp  # noqa: E402

FIRST_SPLIT_SEED = "20261016:domains"  # domain_score.HELD_OUT_SEED; the first split drew chains with it (6c6c089)


def main() -> int:
    sys.path.insert(0, str(TREE / "src"))
    import numpy as np

    import foldmap
    from foldmap import domain_score as ds

    where = Path(foldmap.__file__).resolve()
    if TREE.resolve() not in where.parents:
        raise SystemExit(f"foldmap imported from {where}, not from the harness tree {TREE}")
    assert ds.HELD_OUT_SEED == FIRST_SPLIT_SEED, ds.HELD_OUT_SEED

    data = TREE / "tests" / "data" / "domains"
    out = {"targets": ds.TARGETS, "tolerance": ds.TOLERANCE, "foldmap_from": str(where.parent)}
    out["source_stamp"] = run_stamp(__file__)  # the code that scores, stamped as the file is written (plan 28 Task 1)
    for name, ref, fixture in (
        ("reference", "cath_reference.json", "cath_inputs.npz"),
        ("fresh", "cath_fresh.json", "fresh_inputs.npz"),
    ):
        entries = json.loads((data / ref).read_text(encoding="utf-8"))["chains"]
        inputs = ds.load_fixture(data / fixture)
        rows = []
        for e in entries:
            inp = inputs[f"{e['pdb']}_{e['chain']}"]
            got = ds.inputs_score(inp, e)
            refl = ds.reference_labels([l.seq for l in inp.labels if l.chain == e["chain"]], e["domains"])
            one = np.zeros(len(refl), dtype=int)
            got["baseline_overlap"] = ds.overlap(refl, one)
            got["group"] = e.get("group")
            got["superfamilies"] = e.get("superfamilies")
            rows.append(got)
        out[name] = rows

    keys = sorted((r["pdb"], r["chain"]) for r in out["reference"])
    held = set(random.Random(FIRST_SPLIT_SEED).sample(keys, len(keys) // 3))
    tuning_groups = {r["group"] for r in out["reference"] if (r["pdb"], r["chain"]) not in held}
    group_held = ds.held_out(ds.load_reference(data / "cath_reference.json"))
    for r in out["reference"]:
        k = (r["pdb"], r["chain"])
        r["split"] = "held" if k in held else "tuning"
        r["held_nohom"] = k in held and r["group"] not in tuning_groups
        r["group_held"] = k in group_held
    (BUILD / "domains.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    for split in ("tuning", "held"):
        part = [r for r in out["reference"] if r["split"] == split]
        print(split, len(part), {k: round(v, 3) for k, v in ds.summary(part).items() if k != "n"})
    part = [r for r in out["reference"] if r["held_nohom"]]
    print("held_nohom", len(part), {k: round(v, 3) for k, v in ds.summary(part).items() if k != "n"})
    print("fresh", len(out["fresh"]), {k: round(v, 3) for k, v in ds.summary(out["fresh"]).items() if k != "n"})
    return 0


if __name__ == "__main__":
    sys.exit(main())
