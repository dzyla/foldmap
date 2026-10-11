"""Drawing only some chains of a large structure."""

from pathlib import Path

import pytest

from foldmap.cli import main, make_layout

DATA = Path(__file__).parent / "data"


def test_only_the_chosen_chains_are_drawn():
    lay, sses, bb = make_layout(DATA / "2HHB.cif", assembly="asu", chains=["A", "B"])
    assert {s.chain for s in sses} == {"A", "B"} and set(lay.res_chain) == {"A", "B"}
    assert all(g.chain in ("A", "B") for g in lay.links.ligands)  # hemes of chains C and D go too


def test_unknown_chain_is_explained():
    with pytest.raises(ValueError, match="Q"):
        make_layout(DATA / "2HHB.cif", assembly="asu", chains=["Q"])


def test_cli_chains(tmp_path):
    out = tmp_path / "ab.svg"
    assert main(["plot", str(DATA / "2HHB.cif"), "--assembly", "asu", "--chains", "A", "-o", str(out)]) == 0
    assert "Chain A" in out.read_text() and "Chain B" not in out.read_text()


NA = DATA / "na"


def test_nucleic_acid_chains_can_be_chosen_too():
    """4YBB drew the neighbouring ribosome's 23S rRNA through a lattice contact: `chains` names nucleic-acid chains
    as well as protein ones. Each kind named is restricted to the chains named; a kind not named is drawn whole
    (protein chains alone keep every strand, as before)."""
    lay, sses, _ = make_layout(NA / "1EYG.cif.gz", chains=["Q"])
    assert {s.chain for s in lay.nucleic.strands} == {"Q"} and {s.chain for s in sses} == {"A", "B", "C", "D"}
    assert all(lay.nucleic.strands[k].chain == "Q" for v in lay.nucleic.contacts.values() for k, _ in v)
    lay, sses, _ = make_layout(NA / "1EYG.cif.gz", chains=["A", "R"])
    assert {s.chain for s in lay.nucleic.strands} == {"R"} and {s.chain for s in sses} == {"A"}
    assert len(lay.nucleic.strands) == 2 and len(lay.nucleic.gaps) == 1  # R's chain break is kept
    lay, sses, _ = make_layout(NA / "1EYG.cif.gz", chains=["A"])
    assert {s.chain for s in lay.nucleic.strands} == {"Q", "R"} and {s.chain for s in sses} == {"A"}


def test_an_unknown_chain_names_both_kinds_present():
    with pytest.raises(ValueError, match=r"no chain X; .*protein chains A, B, C, D.*nucleic-acid chains Q, R"):
        make_layout(NA / "1EYG.cif.gz", chains=["Q", "X"])
