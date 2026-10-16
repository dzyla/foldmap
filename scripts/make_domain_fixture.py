# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Build the domain calibration fixtures: what automatic domains are found from (domain_score.ChainInputs: CA
coordinates, residue numbers, elements, their contacts and βα rows) for each chain of a reference set, so the ruled
targets are checked in CI without the structures (plan 12 Task 3, fix round 1).

    python scripts/make_domain_fixture.py                  # cath_reference.json -> tests/data/domains/cath_inputs.npz
    python scripts/make_domain_fixture.py --set fresh      # cath_fresh.json -> tests/data/domains/fresh_inputs.npz

The structures come from $FOLDMAP_DOMAIN_CACHE (default ~/.cache/foldmap-domains; scripts/make_domain_cache.py fills
it). Every chain is checked: its fixture must give the same domains as its structure file.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from make_domain_cache import cache_dir  # noqa: E402

from foldmap import domain_score  # noqa: E402

DATA = ROOT / "tests" / "data" / "domains"
SETS = {"reference": ("cath_reference.json", "cath_inputs.npz"), "fresh": ("cath_fresh.json", "fresh_inputs.npz")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--set", choices=sorted(SETS), default="reference")
    args = ap.parse_args(argv)
    source, target = SETS[args.set]
    entries = json.loads((DATA / source).read_text(encoding="utf-8"))["chains"]
    chains = {}
    for e in entries:
        inputs = domain_score.chain_inputs(cache_dir() / f"{e['pdb']}.cif", e["chain"])
        chains[f"{e['pdb']}_{e['chain']}"] = inputs
    domain_score.save_fixture(DATA / target, chains)
    back = domain_score.load_fixture(DATA / target)
    for key, inputs in chains.items():
        if domain_score.inputs_domains(back[key]) != domain_score.inputs_domains(inputs):
            raise SystemExit(f"{key}: the fixture's domains differ from the structure's")
    print(f"{DATA / target}: {len(chains)} chains, {(DATA / target).stat().st_size / 1024:.0f} KB, domains alike")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
