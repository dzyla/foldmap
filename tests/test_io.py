from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.io import load_backbone

DATA = Path(__file__).parent / "data"


def test_ubiquitin_backbone(ubq):
    bb = load_backbone(ubq)
    assert len(bb) == 76
    assert bb.xyz.shape == (76, 4, 3)
    assert bb.prev.sum() == 75 and not bb.prev[0]
    assert bb.labels[0].name == "MET" and bb.labels[0].seq == 1
    assert np.isfinite(bb.xyz).all()


def test_multichain_break_and_glycans_skipped(zs5):
    bb = load_backbone(zs5)
    assert len(bb) == 257  # glycan chains B, C contribute nothing
    assert {l.chain for l in bb.labels} == {"A", "D"}
    first_d = next(i for i, l in enumerate(bb.labels) if l.chain == "D")
    assert not bb.prev[first_d]


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_backbone(tmp_path / "nope.cif")


def test_no_protein_raises(tmp_path):
    f = tmp_path / "water.pdb"
    f.write_text("HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nEND\n")
    with pytest.raises(ValueError, match="no protein residues"):
        load_backbone(f)


def test_pdb_format_matches_mmcif(ubq, tmp_path):
    out = tmp_path / "u.pdb"
    gemmi.read_structure(str(ubq)).write_pdb(str(out))
    a, b = load_backbone(ubq), load_backbone(out)
    assert a.labels == b.labels
    assert np.allclose(a.xyz, b.xyz, atol=1e-3)


def test_only_first_model_used(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    second = st[0].clone()
    second.num = 2
    for chain in second:
        for res in chain:
            for atom in res:
                atom.pos = gemmi.Position(atom.pos.x + 100, atom.pos.y, atom.pos.z)
    st.add_model(second)
    out = tmp_path / "two_models.cif"
    st.make_mmcif_document().write_file(str(out))
    assert len(gemmi.read_structure(str(out))) == 2  # the fixture really has two models
    assert np.allclose(load_backbone(out).xyz, load_backbone(ubq).xyz, atol=1e-3)


def test_modified_residue_mid_chain_is_kept(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    res = next(r for r in st[0]["A"] if r.seqid.num == 30)
    res.name, res.het_flag = "MSE", "H"  # e.g. selenomethionine recorded as HETATM
    out = tmp_path / "mse.cif"
    st.make_mmcif_document().write_file(str(out))
    bb = load_backbone(out)
    assert len(bb) == 76
    i30 = next(i for i, l in enumerate(bb.labels) if l.seq == 30)
    assert bb.prev[i30] and bb.prev[i30 + 1]


def test_free_hetero_residue_with_backbone_atoms_is_not_a_chain_member(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    chain = st[0]["A"]
    lone = gemmi.Residue()
    lone.name, lone.het_flag, lone.seqid = "ALA", "H", gemmi.SeqId(500, " ")
    for name, x in (("N", 100.0), ("CA", 101.4), ("C", 102.8), ("O", 103.2)):
        atom = gemmi.Atom()
        atom.name, atom.element = name, gemmi.Element(name[0])
        atom.pos = gemmi.Position(x, 100, 100)
        lone.add_atom(atom)
    chain.add_residue(lone)
    out = tmp_path / "lone.cif"
    st.make_mmcif_document().write_file(str(out))
    assert len(load_backbone(out)) == 76


def test_microheterogeneity_keeps_one_residue_per_position(ubq, tmp_path):
    pdb = tmp_path / "u.pdb"
    gemmi.read_structure(str(ubq)).write_pdb(str(pdb))
    out_lines = []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and int(line[22:26]) == 30:
            out_lines.append(line[:16] + "A" + line[17:])
            out_lines.append(line[:16] + "B" + "LEU" + line[20:])
        else:
            out_lines.append(line)
    pdb.write_text("\n".join(out_lines) + "\n")
    bb = load_backbone(pdb)
    assert len(bb) == 76
    assert bb.prev.sum() == 75


def test_residue_missing_backbone_atom_is_skipped_and_breaks_chain(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    res = next(r for r in st[0]["A"] if r.seqid.num == 30)
    res.remove_atom("CA", "*")  # a missing O is rebuilt instead (test_robustness)
    out = tmp_path / "gap.cif"
    st.make_mmcif_document().write_file(str(out))
    bb = load_backbone(out)
    assert len(bb) == 75
    i31 = next(i for i, l in enumerate(bb.labels) if l.seq == 31)
    assert not bb.prev[i31]  # 29 and 31 are no longer bonded


def _stretched_ubiquitin(tmp_path, gap):
    """1UBQ with residues 30 onwards moved so that C(29)-N(30) is `gap` Å."""
    st = gemmi.read_structure(str(Path(__file__).parent / "data" / "1UBQ.cif"))
    chain = st[0][0]
    c = np.array(chain[28]["C"][0].pos.tolist())  # Ile 29
    n = np.array(chain[29]["N"][0].pos.tolist())  # Lys 30
    shift = (n - c) / np.linalg.norm(n - c) * (gap - np.linalg.norm(n - c))
    for res in list(chain)[29:]:
        for atom in res:
            atom.pos = gemmi.Position(*(np.array(atom.pos.tolist()) + shift))
    path = tmp_path / f"gap{gap}.cif"
    st.make_mmcif_document().write_file(str(path))
    return load_backbone(path, "asu")


def test_peptide_bond_up_to_2_5_angstrom_as_dssp(tmp_path):
    bb = _stretched_ubiquitin(tmp_path, 2.2)
    k = next(i for i, l in enumerate(bb.labels) if l.seq == 30)
    assert bb.prev[k]
    bb = _stretched_ubiquitin(tmp_path, 2.7)
    assert not bb.prev[next(i for i, l in enumerate(bb.labels) if l.seq == 30)]


def test_gfp_chromophore_keeps_the_chain_whole():
    from foldmap.cli import make_layout
    from foldmap.dssp import assign_dssp
    from foldmap.model import ResLabel
    from foldmap.route import route_loops
    from foldmap.ss import build_sses

    bb = load_backbone(DATA / "1EMA.cif", "asu")
    k = bb.labels.index(ResLabel("A", 66, "", "CRO"))
    assert bb.prev[k] and bb.prev[k + 1] and bb.linked == frozenset({k})
    d = assign_dssp(bb)
    assert "A:57-64" in {s.id for s in build_sses(bb, d.ss, bridges=d.bridges)}  # the central helix, one piece
    lay, sses, bb = make_layout(DATA / "1EMA.cif")
    assert not [lp for lp in route_loops(lay, sses, bb) if lp.dashed]  # no "missing residues" at 65-67
    assert lay.res_linked == [k]


def test_only_peptide_linked_nonstandard_residues_join_the_chain(tmp_path):
    st = gemmi.read_structure(str(DATA / "1EMA.cif"))
    chain = st[0]["A"]
    k = next(i for i, r in enumerate(chain) if r.name == "CRO")
    for a in chain[k]:
        a.pos = gemmi.Position(a.pos.x + 30, a.pos.y, a.pos.z)  # moved off: bonded to neither neighbour
    path = tmp_path / "loose.cif"
    st.make_mmcif_document().write_file(str(path))
    bb = load_backbone(path, "asu")
    assert not bb.linked and all(l.name != "CRO" for l in bb.labels)


def test_backbone_only_keeps_and_remaps_linked_residues():
    bb = load_backbone(DATA / "1EMA.cif", "asu")
    (k,) = bb.linked
    assert bb.only({"A"}).linked == bb.linked
    assert not bb.only({"Z"}).linked
