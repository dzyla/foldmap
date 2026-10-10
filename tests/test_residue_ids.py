"""Residue identity: insertion codes and negative numbers survive ids, references and layout files."""

import string

import gemmi
import pytest

from foldmap.cli import main, make_layout
from foldmap.layout import resolve
from foldmap.model import ResLabel, res_ref
from helpers import pipeline


def _renumber(src, tmp_path, name, new_seqid):
    """A copy of `src` with every residue of chain A renumbered by new_seqid(k, residue)."""
    st = gemmi.read_structure(str(src))
    st.remove_ligands_and_waters()
    for k, r in enumerate(st[0]["A"]):
        r.seqid = new_seqid(k, r)
    st.setup_entities()
    out = tmp_path / name
    st.make_mmcif_document().write_file(str(out))
    return out


def _icodes(ubq, tmp_path):
    """Ubiquitin with residues 1-17 (strands β1 and β2) all numbered 1, insertion codes A..Q (Kabat-style)."""
    return _renumber(
        ubq, tmp_path, "icodes.cif", lambda k, r: gemmi.SeqId(1, string.ascii_uppercase[k]) if k < 17 else r.seqid
    )


def _negative(ubq, tmp_path):
    return _renumber(ubq, tmp_path, "negative.cif", lambda k, r: gemmi.SeqId(r.seqid.num - 30, " "))


def test_tag_and_res_ref_carry_sign_and_insertion_code():
    a, b = ResLabel("H", 100, "A", "GLY"), ResLabel("H", 100, "K", "TYR")
    assert a.tag == "100A"
    assert res_ref("H", a, b) == "res:H:100A-100K"
    assert res_ref("A", ResLabel("A", -7, "", "ALA"), ResLabel("A", 4, "", "ALA")) == "res:A:-7-4"


def test_insertion_codes_keep_element_ids_distinct(ubq, tmp_path):
    _, sses, _ = pipeline(_icodes(ubq, tmp_path))
    ids = [s.id for s in sses]
    assert len(ids) == len(set(ids))
    assert ids[:2] == ["A:1B-1G", "A:1L-1P"]


def test_layout_with_insertion_codes_builds_and_resolves(ubq, tmp_path):
    lay, sses, _ = make_layout(_icodes(ubq, tmp_path))
    assert "A:1B-1G" in lay.placed and "A:1L-1P" in lay.placed
    assert resolve(lay, sses, "res:A:1L-1P") == ["A:1L-1P"]
    assert resolve(lay, sses, "res:A:1B-1G") == ["A:1B-1G"]


def test_ranges_without_insertion_codes_stop_before_the_insertions(ubq, tmp_path):
    lay, sses, _ = make_layout(_icodes(ubq, tmp_path))
    assert resolve(lay, sses, "res:A:1A-1Z") == ["A:1B-1G", "A:1L-1P"]
    with pytest.raises(ValueError, match="no element"):
        resolve(lay, sses, "res:A:1")  # residue 1 itself, without an insertion code, is in no element


def test_pdb_format_keeps_insertion_codes(ubq, tmp_path):
    cif = _icodes(ubq, tmp_path)
    pdb = tmp_path / "icodes.pdb"
    gemmi.read_structure(str(cif)).write_pdb(str(pdb))
    assert [s.id for s in pipeline(pdb)[1]] == [s.id for s in pipeline(cif)[1]]


def test_negative_residue_references_resolve(ubq, tmp_path):
    lay, sses, _ = make_layout(_negative(ubq, tmp_path))
    helix = next(s for s in sses if s.kind == "H")
    assert helix.id == "A:-7-4"
    assert resolve(lay, sses, "res:A:-7-4") == [helix.id]
    assert resolve(lay, sses, "res:-7--3") == [helix.id]


def test_plain_residue_references_still_work(ubq):
    lay, sses, _ = make_layout(ubq)
    assert resolve(lay, sses, "res:A:23-34") == ["A:23-34"]
    assert resolve(lay, sses, "res:30") == ["A:23-34"]


def test_layout_file_with_negative_numbers_reloads(ubq, tmp_path):
    src = _negative(ubq, tmp_path)
    saved = tmp_path / "fig.yaml"
    assert (
        main(["plot", str(src), "--rename", "α1=H1", "--save-layout", str(saved), "-o", str(tmp_path / "a.svg")]) == 0
    )
    assert "res:A:-7-4: H1" in saved.read_text(encoding="utf-8")
    assert main(["plot", str(src), "--layout-file", str(saved), "-o", str(tmp_path / "b.svg")]) == 0
    assert ">H1<" in (tmp_path / "b.svg").read_text(encoding="utf-8")


def test_layout_file_with_insertion_codes_reloads(ubq, tmp_path):
    src = _icodes(ubq, tmp_path)
    saved = tmp_path / "fig.yaml"
    assert (
        main(["plot", str(src), "--rename", "#2=B2", "--save-layout", str(saved), "-o", str(tmp_path / "a.svg")]) == 0
    )
    assert "res:A:1L-1P: B2" in saved.read_text(encoding="utf-8")
    assert main(["plot", str(src), "--layout-file", str(saved), "-o", str(tmp_path / "b.svg")]) == 0
