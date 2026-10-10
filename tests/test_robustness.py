"""Inputs and situations the pipeline must survive: missing atoms, DNA with no protein elements, assemblies
named by subchain, concurrent drawing, files that are not structures."""

from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.cli import main
from foldmap.frame import view_frame
from foldmap.io import load_backbone, load_nucleic
from foldmap.layout import build_layout
from helpers import pipeline

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


def _strip(src, tmp_path, name, keep_atom):
    """A copy of `src` (chain A) keeping only the atoms for which keep_atom(residue index, atom name) is true."""
    st = gemmi.read_structure(str(src))
    st.remove_ligands_and_waters()
    for k, r in enumerate(st[0]["A"]):
        for i in range(len(r) - 1, -1, -1):
            if not keep_atom(k, r[i].name):
                r.remove_atom(r[i].name, r[i].altloc)
    st.setup_entities()
    out = tmp_path / name
    st.make_mmcif_document().write_file(str(out))
    return out


def test_missing_carbonyl_oxygens_are_rebuilt(ubq, tmp_path):
    holed = _strip(ubq, tmp_path, "noO.cif", lambda k, atom: not (k % 5 == 2 and atom == "O"))
    full, rebuilt = load_backbone(ubq), load_backbone(holed)
    assert len(rebuilt) == 76 and rebuilt.prev.sum() == 75  # no residue lost, no false chain break
    assert np.linalg.norm(full.xyz[:, 3] - rebuilt.xyz[:, 3], axis=1).max() < 0.3  # Å, close to the real O
    assert [s.id for s in pipeline(holed)[1]] == [s.id for s in pipeline(ubq)[1]]


def test_last_residue_missing_oxygen_is_dropped(ubq, tmp_path):
    holed = _strip(ubq, tmp_path, "noO76.cif", lambda k, atom: not (k == 75 and atom == "O"))
    bb = load_backbone(holed)
    assert len(bb) == 75 and bb.labels[-1].seq == 75  # no successor to place its O from


def test_missing_oxygen_before_a_chain_break_is_not_rebuilt_across_it(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    st.remove_ligands_and_waters()
    chain = st[0]["A"]
    next(r for r in chain if r.seqid.num == 30).remove_atom("O", "*")
    del chain[next(k for k, r in enumerate(chain) if r.seqid.num == 31)]
    st.setup_entities()
    out = tmp_path / "gap.cif"
    st.make_mmcif_document().write_file(str(out))
    bb = load_backbone(out)
    assert [l.seq for l in bb.labels if 28 <= l.seq <= 33] == [28, 29, 32, 33]


def test_ca_only_model_is_explained(ubq, tmp_path):
    trace = _strip(ubq, tmp_path, "ca.cif", lambda k, atom: atom == "CA")
    with pytest.raises(ValueError, match="only Cα atoms"):
        load_backbone(trace)
