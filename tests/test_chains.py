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
