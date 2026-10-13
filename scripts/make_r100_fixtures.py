# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Build the random-100 campaign fixtures in tests/data/r100/ (needs network).

    python scripts/make_r100_fixtures.py               # download from RCSB, write tests/data/r100/*.cif.gz
    python scripts/make_r100_fixtures.py --from DIR    # use DIR/<ID>.cif instead of downloading
    python scripts/make_r100_fixtures.py 7FAB-HL       # only these fixtures

Waters and hydrogens are dropped, and the water subchains leave the assembly generators (gemmi cannot build from a
subchain the file lacks). An entry in WINDOWS keeps only the chains it names, each whole or cut to a residue window;
the helix records keep their deposited ids. Each fixture stays under 100 KB and all of them under 600 KB
(tests/test_r100_fixtures.py). Files are written gzipped, outside tests/data/*.cif, so scripts/bench.py never reads
them.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import gemmi
from make_fidelity_fixtures import _fetch, _keep_helix_ids

OUT = Path(__file__).resolve().parent.parent / "tests" / "data" / "r100"
ENTRIES = {  # fixture name -> (PDB id, what it is for)
    "7FAB-HL": ("7FAB", "Fab New, chains H and L whole: disulfide bars over strand labels (Task 6, K11)"),
}
WINDOWS: dict[str, dict[str, tuple[int, int] | None]] = {  # fixture -> {chain: (first, last residue) or None: whole}
    "7FAB-HL": {"H": None, "L": None},
}


def build(text: str, window: dict[str, tuple[int, int] | None] | None = None) -> bytes:
    st = gemmi.read_structure_string(text)
    st.remove_waters()
    st.remove_hydrogens()
    if window:
        for model in st:
            for name in [c.name for c in model if c.name not in window]:
                model.remove_chain(name)
            for c in model:
                span = window[c.name]
                if span is None:
                    continue
                for i in reversed(range(len(c))):
                    if not span[0] <= c[i].seqid.num <= span[1]:
                        del c[i]
    kept = {r.subchain for c in st[0] for r in c}
    for assembly in st.assemblies:  # gemmi refuses to build from a subchain the file no longer has (the waters)
        for g in assembly.generators:
            g.subchains = [s for s in g.subchains if s in kept]
    doc = st.make_mmcif_document()
    _keep_helix_ids(gemmi.cif.read_string(text).sole_block(), doc.sole_block(), {c.name for c in st[0]})
    return gzip.compress(doc.as_string().encode("utf-8"), 9, mtime=0)


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
    for name in args.names or ENTRIES:
        pdb, _ = ENTRIES[name]
        data = build(_fetch(pdb, args.src), WINDOWS.get(name))
        (OUT / f"{name}.cif.gz").write_bytes(data)
        total += len(data)
        print(f"{name}.cif.gz: {len(data) / 1024:.0f} KB")
    print(f"total {total / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
