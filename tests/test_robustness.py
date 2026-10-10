"""Inputs and situations the pipeline must survive: missing atoms, DNA with no protein elements, assemblies
named by subchain, concurrent drawing, files that are not structures."""

from pathlib import Path

import gemmi

from foldmap.cli import main
from foldmap.frame import view_frame
from foldmap.io import load_backbone, load_nucleic
from foldmap.layout import build_layout

DATA = Path(__file__).parent / "data"


def test_dna_is_laid_out_when_the_protein_has_no_elements():
    path = DATA / "1LMB.cif"
    bb = load_backbone(path, "asu")
    na = load_nucleic(path, bb, "asu")
    for mode in ("projected", "stack"):
        lay = build_layout([], [], view_frame([]), mode, nucleic=na)
        assert len(lay.dna) == 1 and lay.bounds[2] > lay.bounds[0]


def test_peptides_on_dna_draw_the_dna(tmp_path):
    st = gemmi.read_structure(str(DATA / "1LMB.cif"))
    st.remove_ligands_and_waters()
    for c in st[0]:
        if c.name in ("3", "4"):  # the repressor chains, cut to 6 residues: no helix or strand left
            while len(c) > 6:
                del c[len(c) - 1]
    st.setup_entities()
    src = tmp_path / "peptides_dna.cif"
    st.make_mmcif_document().write_file(str(src))
    out = tmp_path / "f.svg"
    assert main(["plot", str(src), "-o", str(out), "--assembly", "asu"]) == 0
    svg = out.read_text(encoding="utf-8")
    assert 'id="nt:' in svg and "no helices or strands found" not in svg
