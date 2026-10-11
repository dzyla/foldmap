# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Build the protein fixtures for the published-figure comparisons in tests/data/fidelity/ (needs network).

    python scripts/make_fidelity_fixtures.py               # download from RCSB, write tests/data/fidelity/*.cif.gz
    python scripts/make_fidelity_fixtures.py --from DIR    # use DIR/<ID>.cif instead of downloading
    python scripts/make_fidelity_fixtures.py 1TEN 7TIM     # only these fixtures

Waters and hydrogens are dropped, and the water subchains leave the assembly generators (gemmi cannot build from a
subchain the file lacks); nothing else is trimmed (chains, ligands, assembly records, helix and sheet records stay).
Files are written gzipped, outside tests/data/*.cif, so scripts/bench.py never reads them.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import gemmi

OUT = Path(__file__).resolve().parent.parent / "tests" / "data" / "fidelity"
ENTRIES = {
    "1TEN": "tenascin FN3 domain: a seven-strand sandwich with a bulged strand",
    "1LDM": "lactate dehydrogenase M4 tetramer: a Rossmann fold, built from its assembly",
    "2ACY": "acylphosphatase: a ferredoxin-like fold with a helix over the sheet",
    "1MBN": "myoglobin: all helical, with a haem",
    "7TIM": "triosephosphate isomerase: a TIM barrel",
}


def _fetch(pdb: str, src: Path | None) -> str:
    if src:
        return (src / f"{pdb}.cif").read_text(encoding="utf-8")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.cif", timeout=60) as r:
        return r.read().decode("utf-8")


def build(text: str) -> bytes:
    st = gemmi.read_structure_string(text)
    st.remove_waters()
    st.remove_hydrogens()
    kept = {r.subchain for c in st[0] for r in c}
    for assembly in st.assemblies:  # gemmi refuses to build from a subchain the file no longer has (the waters)
        for g in assembly.generators:
            g.subchains = [s for s in g.subchains if s in kept]
    return gzip.compress(st.make_mmcif_document().as_string().encode("utf-8"), 9, mtime=0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("names", nargs="*", help="fixtures to build (default: all)")
    ap.add_argument("--from", dest="src", type=Path, help="read DIR/<ID>.cif instead of downloading")
    args = ap.parse_args(argv)
    unknown = set(args.names) - set(ENTRIES)
    if unknown:
        ap.error(f"unknown fixture {', '.join(sorted(unknown))}")
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    for pdb in args.names or ENTRIES:
        data = build(_fetch(pdb, args.src))
        (OUT / f"{pdb}.cif.gz").write_bytes(data)
        total += len(data)
        print(f"{pdb}.cif.gz: {len(data) / 1024:.0f} KB")
    print(f"total {total / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
