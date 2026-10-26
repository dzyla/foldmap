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
    # HIV-1 RT polypurine-tract RNA:DNA hybrid alone, around E868·F886 (one stretched hydrogen bond)
    "1HYS": (["E", "F"], {"E": (863, 876), "F": (878, 891)}, []),
    "1AOI": (["I", "J"], {}, []),  # nucleosome DNA alone: one 146-bp duplex with unpaired steps at I20-21, I24-25
}
# entries whose assembly is kept: 3CRX's operator lays a copy of every DNA strand on a deposited one (the junction
# sits on the crystallographic two-fold), which io.read_model must drop
ASSEMBLIES = {"3CRX", "3EGZ"}
# fixture name -> (id, chains to keep or None for all, {chain: (first, last)} residue windows)
CONTROLS = {
    "1ATP": ("1ATP", None, {}),  # protein kinase A with ATP: a nucleotide ligand, never a strand
    "4YBB-psi": ("4YBB", ["DA"], {"DA": (2586, 2608)}),  # 23S rRNA hairpin: pseudouridine 2604 pairs with A2590
    "2OEU": ("2OEU", None, {}),  # hammerhead ribozyme: three stems placed by the strands between them (P2 and P3)
    "2GIS": ("2GIS", ["A"], {}),  # SAM-I riboswitch, RNA only: internal loops and runs between seven stems
    "1URN": ("1URN", ["C", "R"], {}),  # U1A on an RNA hairpin whose loop is not modelled (U13-C14): a chain break
    "1QBJ": ("1QBJ", None, {}),  # Zalpha on Z-DNA: strand F pairs only with a crystal mate, drawn as a trace
    "1EYG": ("1EYG", None, {}),  # E. coli SSB tetramer with ssDNA wrapped over it: unpaired traces on their protein
}

# plan 28's nucleic-acid fixes (docs/validation/na-fixes): fixture name -> (id, chains or None, windows, slim: None or
# SLIM, the protein far from the nucleic acid cut to its backbone and CB; anisotropic B factors are dropped), written to
# tests/data/na_fixes/ (OUT_FIXES), apart from the fixtures the independent kernels compare (tests/data/na/*)
OUT_FIXES = OUT.parent / "na_fixes"
SLIM = 8.0  # Å: a protein residue farther than this from every nucleotide keeps N, CA, C, O and CB only (_slim)
FIXES = {
    # F2: a C4 assembly (kept: ASSEMBLIES) of four protein-RNA copies, each RNA held by one protomer; unrolled, each
    # RNA must stay with its protomer (one stem was unrolled into the next protomer, copy B1 torn across the seam)
    "3EGZ": ("3EGZ", None, {}, None),
    # F4: a protein (A) on its RNA (B), B:32-35 not modelled; the gap from B:31 crossed the 5' end B:36-40
    "4AS1": ("4AS1", ["A", "B"], {}, SLIM),
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


def control(st: gemmi.Structure, chains, windows, keep_assemblies: bool = False) -> gemmi.Structure:
    model = _first_model(st, windows)
    for name in [c.name for c in model]:
        if chains is not None and name not in chains:
            model.remove_chain(name)
    _finish(st, keep_assemblies)
    kept = {r.subchain for c in model for r in c}
    names = {c.name for c in model}
    for k in reversed(range(len(st.entities))):  # entities and links of dropped chains only add bulk
        if not kept & set(st.entities[k].subchains):
            del st.entities[k]
    for k in reversed(range(len(st.connections))):
        if not {st.connections[k].partner1.chain_name, st.connections[k].partner2.chain_name} <= names:
            del st.connections[k]
    return st


def _slim(st: gemmi.Structure, near: float) -> gemmi.Structure:
    """Protein residues with no atom within `near` of a nucleotide cut to N, CA, C, O and CB: secondary structure
    (DSSP reads the backbone) and every protein-nucleic-acid contact (io.CONTACT, 4 Å) as in the entry, at a fraction
    of the size."""
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, near + 1).populate()
    close = set()
    for chain in model:
        if _is_nucleic(chain):
            for res in chain:
                for atom in res:
                    for mark in ns.find_atoms(atom.pos, "\0", radius=near):
                        cra = mark.to_cra(model)
                        close.add((cra.chain.name, cra.residue.seqid.num, cra.residue.seqid.icode))
    for chain in model:
        if chain.get_polymer().check_polymer_type() != gemmi.PolymerType.PeptideL:
            continue
        for res in chain:
            if (chain.name, res.seqid.num, res.seqid.icode) not in close:
                for k in reversed(range(len(res))):
                    if res[k].name not in ("N", "CA", "C", "O", "CB"):
                        del res[k]
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


def _write(name: str, st: gemmi.Structure, out: Path = OUT, lean: bool = False) -> int:
    path = out / f"{name}.cif.gz"
    doc = st.make_mmcif_document()
    if lean:  # anisotropic B: no figure reads it
        doc.sole_block().find_mmcif_category("_atom_site_anisotrop.").erase()
    cif = doc.as_string()
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
    unknown = set(args.names) - set(ENTRIES) - set(CONTROLS) - set(FIXES)
    if unknown:
        ap.error(f"unknown fixture {', '.join(sorted(unknown))}")
    wanted = set(args.names) or set(ENTRIES) | set(CONTROLS) | set(FIXES)
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
            st = gemmi.read_structure_string(_fetch(pdb, args.src))
            total += _write(name, control(st, chains, windows, name in ASSEMBLIES))
    OUT_FIXES.mkdir(parents=True, exist_ok=True)
    for name, (pdb, chains, windows, slim) in FIXES.items():
        if name in wanted:
            st = gemmi.read_structure_string(_fetch(pdb, args.src))
            st = control(st, chains, windows, name in ASSEMBLIES)
            total += _write(name, _slim(st, slim) if slim else st, OUT_FIXES, lean=bool(slim))
    if wanted & set(ENTRIES):
        ref_path.write_text(json.dumps(reference, separators=(",", ":")) + "\n", encoding="utf-8")
        print(f"reference pairs in {ref_path}")
    print(f"total {total / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
