# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Build the protein fixtures for the published-figure comparisons in tests/data/fidelity/ (needs network).

    python scripts/make_fidelity_fixtures.py               # download from RCSB, write tests/data/fidelity/*.cif.gz
    python scripts/make_fidelity_fixtures.py --from DIR    # use DIR/<ID>.cif instead of downloading
    python scripts/make_fidelity_fixtures.py 1TEN 7TIM     # only these fixtures

Waters and hydrogens are dropped, and the water subchains leave the assembly generators (gemmi cannot build from a
subchain the file lacks); nothing else is trimmed (chains, ligands, assembly records, helix and sheet records stay,
the helices with their deposited ids), except for the entries in WINDOWS, cut down to one stretch of one chain.
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
    "1PV7": "lactose permease, chain B 70-115 only: a pi-bulge in a straight TM helix (H 87-90, I 91-95, H 96-100)",
}
WINDOWS = {"1PV7": ("B", 70, 115)}  # chain, first and last residue kept (DSSP of the window matches the full chain)


def _fetch(pdb: str, src: Path | None) -> str:
    if src:
        return (src / f"{pdb}.cif").read_text(encoding="utf-8")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.cif", timeout=60) as r:
        return r.read().decode("utf-8")


def build(text: str, window: tuple[str, int, int] | None = None) -> bytes:
    st = gemmi.read_structure_string(text)
    st.remove_waters()
    st.remove_hydrogens()
    if window:
        chain, lo, hi = window
        for model in st:
            for name in [c.name for c in model if c.name != chain]:
                model.remove_chain(name)
            for c in model:
                for i in reversed(range(len(c))):
                    if not lo <= c[i].seqid.num <= hi:
                        del c[i]
    kept = {r.subchain for c in st[0] for r in c}
    for assembly in st.assemblies:  # gemmi refuses to build from a subchain the file no longer has (the waters)
        for g in assembly.generators:
            g.subchains = [s for s in g.subchains if s in kept]
    doc = st.make_mmcif_document()
    _keep_helix_ids(gemmi.cif.read_string(text).sole_block(), doc.sole_block(), {c.name for c in st[0]})
    return gzip.compress(doc.as_string().encode("utf-8"), 9, mtime=0)


def _keep_helix_ids(src, out, chains: set[str]) -> None:
    """The deposited _struct_conf rows (of the chains kept) in place of the ones gemmi writes: gemmi renumbers the
    helices H1, H2 ... and drops pdbx_PDB_helix_id, the deposited names (myoglobin's A-H) that --names reads."""
    table = src.find_mmcif_category("_struct_conf.")
    if not table:
        return
    tags = [tag.split(".", 1)[1] for tag in table.tags]
    rows = [list(row) for row in table if row[tags.index("beg_auth_asym_id")] in chains]
    out.set_mmcif_category("_struct_conf.", {tag: [r[k] for r in rows] for k, tag in enumerate(tags)}, raw=True)


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
        data = build(_fetch(pdb, args.src), WINDOWS.get(pdb))
        (OUT / f"{pdb}.cif.gz").write_bytes(data)
        total += len(data)
        print(f"{pdb}.cif.gz: {len(data) / 1024:.0f} KB")
    print(f"total {total / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
