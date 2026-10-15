# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Build the random-100 campaign fixtures in tests/data/r100/ (needs network).

    python scripts/make_r100_fixtures.py               # download from RCSB, write tests/data/r100/*.cif.gz
    python scripts/make_r100_fixtures.py --from DIR    # use DIR/<ID>.cif instead of downloading
    python scripts/make_r100_fixtures.py 7FAB-HL       # only these fixtures

Waters and hydrogens are dropped, and the water subchains leave the assembly generators (gemmi cannot build from a
subchain the file lacks). An entry in WINDOWS keeps only the chains it names, each whole or cut to residue windows; the
helix records keep their deposited ids. TRIMS shrinks an entry further: `assembly` keeps only the subchains of that
biological assembly, `backbone` keeps only a protein's N, CA, C, O and CB atoms (`atoms`: these names instead, e.g. with
the cysteines' SG for disulfides), `first_model` keeps an NMR ensemble's model 1, `drop` erases the mmCIF categories
named, `no_protein_sequence` drops the protein entities' full sequences (a big complex's are most of a small window's
file). Each fixture stays under 100 KB and all of them under 800 KB (tests/test_r100_fixtures.py). Files are written
gzipped, outside tests/data/*.cif, so scripts/bench.py never reads them.
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
    "3OWW-AB": (
        "3OWW",
        "glycine riboswitch, chains A and B around the residues where they pair: two identical RNA entities, "
        "deposited once each, so --assembly auto counted two copies and "
        "drew assembly 1's chain A alone (random-100 C30, 3OWW ruling)",
    ),
    "7FAB-HL": ("7FAB", "Fab New, chains H and L whole: disulfide bars over strand labels (Task 6, K11)"),
    "5LJ3-V": (
        "5LJ3",
        "spliceosome: RNA V whole, Z's stem with it and the protein elements beside them on the full figure's page; "
        "its 3' end V:85-92 ran back over its own stem and buried the gap after V:92 (Task 5, K5)",
    ),
    "4V4F-ring": ("4V4F", "TRAP-RNA: assembly 1, the 11-subunit ring and its RNA, backbone only (Task 5, K24)"),
    "1IGT-AB": ("1IGT", "IgG: a light and a heavy chain, backbone and SG: disulfide bars over labels (Task 6, K11)"),
    "1GOX-A": ("1GOX", "glycolate oxidase, one subunit and its FMN: tethers on a template barrel (Task 6, K8)"),
    "9FWZ-DB": (
        "9FWZ",
        "type 1 pilus usher FimD (D, a 24-strand outer-membrane beta-barrel) and the FimA subunit in its lumen (B), "
        "backbone only: the beta-barrel membrane estimate; B spans the slab without defining it (C-MEM-BARREL)",
    ),
    "1BXW-A": ("1BXW", "OmpA, an 8-strand outer-membrane beta-barrel, backbone only: the control (C-MEM-BARREL)"),
    "2DGS-A": (
        "2DGS",
        "an NMR domain (model 1, backbone only) with a 14-residue N tail whose count text a loop ran through "
        "(random-100 C1)",
    ),
    "9ON5-A": (
        "9ON5",
        "P2X receptor chain A, residues 100-200, backbone only: a long loop's residue count written on a strand "
        "(random-100 C3)",
    ),
    "1BXI": (
        "1BXI",
        "colicin E9 DNase and Im9, backbone, CB and the side-chain atoms that hold a metal: the Ni mark was put on the "
        "loop B:22-25>B:32-33 (random-100 C5)",
    ),
    "4UM3-two": (
        "4UM3",
        "eight pentamers in the asymmetric unit, cut to two (assembly 1's P-T and K-O), residues 30-40, backbone "
        "only: assembly 1 holds two single-copy entities of the same protein (random-100 C30)",
    ),
    "6B5W": (
        "6B5W",
        "a 19-residue lasso peptide (NMR model 1, backbone only) with no helix or strand: nothing drawn, which the "
        "note must say (random-100 C15)",
    ),
    "7IMM-A": (
        "7IMM",
        "a dimer (assembly 1: chain A twice) whose two copies make one continuous A B B A sheet, backbone only: the "
        "shared sheet A15 leaves out (random-100 C10)",
    ),
    "8Q4O": (
        "8Q4O",
        "an RNA G-quadruplex (NMR model 1): three stacked G-tetrads, no Watson-Crick pair, so no stem: drawn as a "
        "trace, which the note must say (random-100 C29 ruling)",
    ),
    "3L9G-A": (
        "3L9G",
        "urate oxidase, chain A (assembly 1: four copies by symmetry), backbone only: its copies share sheets, so "
        "the domains view parted their domains for ever and exploded the page (random-100 C31)",
    ),
    "6SO9": (
        "6SO9",
        "an RRM and its bound 6-nt RNA (NMR model 1; protein backbone only): two protein loops ran across the RNA's "
        "letters, routed before the RNA was (random-100 C11)",
    ),
    "2MJW-BD": (
        "2MJW",
        "S100P dimer (chains B and D, NMR model 1, backbone only): D is numbered 95-188 where B is 1-94, so the copies "
        "shared no residue number and were never tied (random-100 C18/A15: B~D 3.34)",
    ),
    "3L2V-A": (
        "3L2V",
        "PFV intasome: the inner integrase chain A cut to three windows that touch the DNA, N, CA, C and O only, and "
        "the viral DNA C and D whole (assembly 1 doubles them): each DNA end is held by both inner chains, so "
        "unrolling the C2 tetramer drew each one turned 18-58 degrees from its 3D axis (random-100 C16)",
    ),
    "9CJI-AB": (
        "9CJI",
        "a Cas12 effector's crRNA (chain B whole) and the protein windows that touch it, backbone only: its 3' end was "
        "drawn turned onto the 5' end's line, nucleotides and beads on top of each other (random-100 C16, A17 "
        "bead_overlaps)",
    ),
    "3OWW-P5": (
        "3OWW",
        "glycine riboswitch chain B, the P5 hairpin (B:12-28), whose loop lacks B22: its 4-nt 5' piece pointed away "
        "from the 3' piece and the dashed gap went round the stem across a run (random-100 C16, A17 na_crossings)",
    ),
}
Span = tuple[int, int]
WINDOWS: dict[str, dict[str, Span | list[Span] | None]] = {  # fixture -> {chain: windows (first, last), None: whole}
    "3OWW-AB": {c: [(5, 12), (38, 45), (60, 70), (80, 88)] for c in "AB"},
    "3L9G-A": {"A": None},
    "3OWW-P5": {"B": (12, 28)},
    "9CJI-AB": {"B": None, "A": [(10, 22), (740, 765), (784, 812), (850, 875)]},
    "3L2V-A": {"A": [(60, 75), (100, 118), (205, 230)], "C": None, "D": None},
    "2MJW-BD": {"B": None, "D": None},  # where the two RNAs pair
    "7FAB-HL": {"H": None, "L": None},
    "1IGT-AB": {"A": None, "B": None},
    "9FWZ-DB": {"D": (116, 700), "B": None},
    "9ON5-A": {"A": (100, 200)},
    "4UM3-two": {c: (30, 40) for c in "KLMNOPQRST"},  # the barrel (139-659), its plug, the start of CTD1; FimA
    "5LJ3-V": {  # found by a seeded search over the page neighbours of the gap (plan 10 Task 5 report)
        "V": None,
        "Z": (1, 16),
        "A": [(208, 232), (514, 536), (1284, 1290), (1019, 1036), (906, 925), (840, 871)],
        "O": [(164, 192), (82, 93)],
        "H": (146, 159),
        "N": (121, 129),
        "T": [(434, 441), (482, 498), (406, 420)],
    },
}
TRIMS: dict[str, dict] = {
    "5LJ3-V": {"no_protein_sequence": True},
    "4V4F-ring": {"assembly": "1", "backbone": True},
    "1IGT-AB": {"atoms": ("N", "CA", "C", "O", "CB", "SG")},
    "9FWZ-DB": {"backbone": True, "no_protein_sequence": True},
    "1BXW-A": {"backbone": True},
    "7IMM-A": {"backbone": True},
    "2DGS-A": {"first_model": True, "backbone": True},
    "9ON5-A": {"backbone": True},
    "4UM3-two": {
        "backbone": True,
        "drop": ("_atom_site_anisotrop.", "_pdbx_refine_tls.", "_pdbx_refine_tls_group.", "_struct_sheet_order."),
    },
    "6B5W": {"first_model": True, "backbone": True},
    "8Q4O": {"first_model": True},
    "2MJW-BD": {"first_model": True, "backbone": True},
    "6SO9": {"first_model": True, "backbone": True},
    "3L9G-A": {"atoms": ("N", "CA", "C", "O"), "no_protein_sequence": True, "drop": ("_atom_site_anisotrop.",)},
    "3OWW-P5": {"drop": ("_atom_site_anisotrop.", "_pdbx_refine_tls.", "_pdbx_refine_tls_group.")},
    "9CJI-AB": {"backbone": True, "no_protein_sequence": True, "drop": ("_atom_site_anisotrop.",)},
    "3L2V-A": {"atoms": ("N", "CA", "C", "O"), "no_protein_sequence": True, "drop": ("_atom_site_anisotrop.",)},
    "1BXI": {
        "atoms": ("N", "CA", "C", "O", "CB", "ND1", "NE2", "SG", "OD1", "OD2", "OE1", "OE2", "NZ", "OG", "OG1", "SD")
    },
}
_BACKBONE = {"N", "CA", "C", "O", "CB"}


def _protein(r: gemmi.Residue) -> bool:
    info = gemmi.find_tabulated_residue(r.name)
    return info is not None and info.is_amino_acid()


def build(text: str, window: dict[str, Span | list[Span] | None] | None = None, trim: dict | None = None) -> bytes:
    trim = trim or {}
    st = gemmi.read_structure_string(text)
    if trim.get("first_model"):  # an NMR ensemble: model 1 only
        while len(st) > 1:
            del st[len(st) - 1]
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
                spans = span if isinstance(span, list) else [span]
                for i in reversed(range(len(c))):
                    if not any(lo <= c[i].seqid.num <= hi for lo, hi in spans):
                        del c[i]
    if "assembly" in trim:
        (assembly,) = [a for a in st.assemblies if a.name == trim["assembly"]]
        sub = {s for g in assembly.generators for s in g.subchains}
        for model in st:
            for c in model:
                for i in reversed(range(len(c))):
                    if c[i].subchain not in sub:
                        del c[i]
    keep = set(trim["atoms"]) if "atoms" in trim else _BACKBONE if trim.get("backbone") else None
    if keep is not None:
        for model in st:
            for r in (r for c in model for r in c if _protein(r)):
                for i in reversed(range(len(r))):
                    if r[i].name not in keep:
                        del r[i]
    st.remove_empty_chains()
    if trim.get("no_protein_sequence"):
        for ent in st.entities:
            if ent.polymer_type == gemmi.PolymerType.PeptideL:
                ent.full_sequence = []
    kept = {r.subchain for c in st[0] for r in c}
    for assembly in st.assemblies:  # gemmi refuses to build from a subchain the file no longer has (the waters)
        for g in assembly.generators:
            g.subchains = [s for s in g.subchains if s in kept]
    doc = st.make_mmcif_document()
    for name in trim.get("drop", ()):  # categories a fixture does not need (anisotropic B, TLS groups)
        doc.sole_block().find_mmcif_category(name).erase()
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
        data = build(_fetch(pdb, args.src), WINDOWS.get(name), TRIMS.get(name))
        (OUT / f"{name}.cif.gz").write_bytes(data)
        total += len(data)
        print(f"{name}.cif.gz: {len(data) / 1024:.0f} KB")
    print(f"total {total / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
