# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Build the nucleic-acid test fixtures in tests/data/na/ (needs network and rnapolis; not run in CI).

    python scripts/make_na_fixtures.py                 # download from RCSB, write tests/data/na/*.cif.gz
    python scripts/make_na_fixtures.py --from DIR      # use DIR/<ID>.cif instead of downloading
    python scripts/make_na_fixtures.py 1ATP 4YBB-psi   # only these fixtures (reference.json keeps the others)

Each entry keeps its nucleic-acid chains (optionally a residue window of them), every protein chain with an atom
within 8 Å of the kept nucleic acid (whole chains, so their secondary structure is intact) and no waters. Only the
entries in ASSEMBLIES keep their assembly records (`_pdbx_struct_assembly`, `_oper_list`, `_assembly_gen`), so
`auto` builds their assembly; the others are read as deposited. Files are written gzipped; gemmi reads them directly.
reference.json holds the canonical and wobble pairs RNApolis finds in each deposited entry (pip install rnapolis),
the independent reference tests/test_na_pairs.py checks foldmap's pairing against. CONTROLS are fixtures for one
property each, with no reference pairs: chains and residue windows kept as listed, waters dropped.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from pathlib import Path

import gemmi

OUT = Path(__file__).resolve().parent.parent / "tests" / "data" / "na"
NEAR = 8.0  # Å

# id -> (nucleic-acid chains to keep or None for all, {chain: (first, last)} residue windows,
#        protein chains allowed or None for any within NEAR)
ENTRIES = {
    "30QM": (None, {}, None),  # E. coli DRT10 RT-ncRNA: 121-nt ncRNA, six stems
    "1EHZ": (None, {}, None),  # yeast tRNA-Phe: cloverleaf, modified nucleotides, G4.U69
    "4OO8": (["B", "C"], {}, ["A"]),  # Cas9-sgRNA-target DNA, first complex: RNA-DNA heteroduplex
    "1KX5": (None, {}, None),  # nucleosome core particle: 147 bp wrapped DNA
    "1IHF": (None, {}, None),  # IHF-DNA: DNA bent by more than 160 degrees
    "3CRX": (None, {}, None),  # Cre-Holliday junction: four arms
    "1RVA": (None, {}, None),  # EcoRV dimer kinking its site by about 50°: the duplex lies across the dimer axis
    "4LD0": (None, {}, None),  # RuvC dimer on a Holliday junction: each arm held by one protomer
    # E. coli 30S subunit: 16S rRNA 1-1000 (5' and central domains) and six proteins that bind them
    # (S4 AD, S8 AH, S15 AO, S16 AP, S17 AQ, S20 AT)
    "4YBB": (["AA"], {"AA": (1, 1000)}, ["AD", "AH", "AO", "AP", "AQ", "AT"]),
}
# entries whose assembly is kept: 3CRX's operator lays a copy of every DNA strand on a deposited one (the junction
# sits on the crystallographic two-fold), which io.read_model must drop
ASSEMBLIES = {"3CRX"}
# fixture name -> (id, chains to keep or None for all, {chain: (first, last)} residue windows)
CONTROLS = {
    "1ATP": ("1ATP", None, {}),  # protein kinase A with ATP: a nucleotide ligand, never a strand
    "4YBB-psi": ("4YBB", ["DA"], {"DA": (2586, 2608)}),  # 23S rRNA hairpin: pseudouridine 2604 pairs with A2590
}


def _is_nucleic(chain) -> bool:
    return chain.get_polymer().check_polymer_type() in (
        gemmi.PolymerType.Dna,
        gemmi.PolymerType.Rna,
        gemmi.PolymerType.DnaRnaHybrid,
    )


def _first_model(st: gemmi.Structure, windows) -> gemmi.Model:
    """Model 1 without waters or hydrogens, the windowed chains cut to their windows."""
    st.remove_waters()
    st.remove_hydrogens()
    while len(st) > 1:
        del st[1]
    model = st[0]
    for name, (lo, hi) in windows.items():
        chain = model[name]
        for k in reversed(range(len(chain))):
            if not lo <= chain[k].seqid.num <= hi:
                del chain[k]
    return model


def _finish(st: gemmi.Structure, keep_assemblies: bool = False) -> gemmi.Structure:
    if not keep_assemblies:
        st.assemblies.clear()
    st.remove_empty_chains()
    kept = {r.subchain for c in st[0] for r in c}
    for assembly in st.assemblies:  # gemmi refuses to build from a subchain the file no longer has (waters)
        for g in assembly.generators:
            g.subchains = [s for s in g.subchains if s in kept]
    st.setup_entities()
    return st


def control(st: gemmi.Structure, chains, windows) -> gemmi.Structure:
    model = _first_model(st, windows)
    for name in [c.name for c in model]:
        if chains is not None and name not in chains:
            model.remove_chain(name)
    _finish(st)
    kept = {r.subchain for c in model for r in c}
    names = {c.name for c in model}
    for k in reversed(range(len(st.entities))):  # entities and links of dropped chains only add bulk
        if not kept & set(st.entities[k].subchains):
            del st.entities[k]
    for k in reversed(range(len(st.connections))):
        if not {st.connections[k].partner1.chain_name, st.connections[k].partner2.chain_name} <= names:
            del st.connections[k]
    return st


def trim(st: gemmi.Structure, na_chains, windows, allowed, keep_assemblies: bool = False) -> gemmi.Structure:
    model = _first_model(st, windows)
    keep_na = [c.name for c in model if _is_nucleic(c) and (na_chains is None or c.name in na_chains)]
    ns = gemmi.NeighborSearch(model, st.cell, NEAR + 1).populate()
    near = set()
    for name in keep_na:
        for res in model[name]:
            for atom in res:
                for mark in ns.find_atoms(atom.pos, "\0", radius=NEAR):
                    near.add(mark.to_cra(model).chain.name)
    proteins = [
        c.name
        for c in model
        if c.name in near
        and (allowed is None or c.name in allowed)
        and c.get_polymer().check_polymer_type() == gemmi.PolymerType.PeptideL
    ]
    for name in [c.name for c in model]:
        if name not in keep_na and name not in proteins:
            model.remove_chain(name)
    return _finish(st, keep_assemblies)


def _reference_pairs(text: str) -> list:
    """Canonical and wobble pairs of the deposited entry by RNApolis (an independent annotator):
    [[chain, number, icode], [chain, number, icode], "WC" or "GU"], sorted."""
    import tempfile

    from rnapolis.annotator import find_pairs
    from rnapolis.parser import read_3d_structure
    from rnapolis.util import handle_input_file

    kinds = {"XIX": "WC", "XX": "WC", "XXVIII": "GU"}  # Saenger classes: G-C, A-U/T, G-U wobble
    with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False, encoding="utf-8") as f:
        f.write(text)
    found = find_pairs(read_3d_structure(handle_input_file(f.name), 1), 1)[0]
    Path(f.name).unlink()
    out = []
    for bp in found:
        kind = kinds.get(bp.saenger.value if bp.saenger is not None else None)
        if kind:
            a, b = bp.nt1.auth, bp.nt2.auth
            pair = sorted([[a.chain, a.number, (a.icode or "").strip()], [b.chain, b.number, (b.icode or "").strip()]])
            out.append([*pair, kind])
    return sorted(out)


def _fetch(pdb: str, src: Path | None) -> str:
    if src:
        return (src / f"{pdb}.cif").read_text(encoding="utf-8")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.cif.gz", timeout=60) as r:
        return gzip.decompress(r.read()).decode("utf-8")


def _write(name: str, st: gemmi.Structure) -> int:
    path = OUT / f"{name}.cif.gz"
    cif = st.make_mmcif_document().as_string()
    path.write_bytes(gzip.compress(cif.encode("utf-8"), 9, mtime=0))  # no timestamp: rerunning is reproducible
    size = path.stat().st_size
    chains = ",".join(c.name for c in st[0])
    print(f"{path.name}: chains {chains}, {size / 1024:.0f} KB")
    return size


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("names", nargs="*", help="fixtures to build (default: all)")
    ap.add_argument("--from", dest="src", type=Path, help="read DIR/<ID>.cif instead of downloading")
    args = ap.parse_args(argv)
    unknown = set(args.names) - set(ENTRIES) - set(CONTROLS)
    if unknown:
        ap.error(f"unknown fixture {', '.join(sorted(unknown))}")
    wanted = set(args.names) or set(ENTRIES) | set(CONTROLS)
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    ref_path = OUT / "reference.json"
    reference: dict[str, list] = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}
    for pdb, (na_chains, windows, allowed) in ENTRIES.items():
        if pdb not in wanted:
            continue
        text = _fetch(pdb, args.src)
        st = trim(gemmi.read_structure_string(text), na_chains, windows, allowed, pdb in ASSEMBLIES)
        kept = {(c.name, r.seqid.num, r.seqid.icode.strip()) for c in st[0] for r in c}
        reference[pdb] = [p for p in _reference_pairs(text) if tuple(p[0]) in kept and tuple(p[1]) in kept]
        total += _write(pdb, st)
    for name, (pdb, chains, windows) in CONTROLS.items():
        if name in wanted:
            total += _write(name, control(gemmi.read_structure_string(_fetch(pdb, args.src)), chains, windows))
    if wanted & set(ENTRIES):
        ref_path.write_text(json.dumps(reference, separators=(",", ":")) + "\n", encoding="utf-8")
        print(f"reference pairs in {ref_path}")
    print(f"total {total / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
