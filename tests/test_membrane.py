"""Membrane bands: the lipid bilayer drawn behind transmembrane proteins, with its inside and outside."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main, make_figure, make_layout
from foldmap.io import load_backbone
from foldmap.membrane import find_membrane

DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("name", ["1BL8", "1C3W"])
def test_membrane_proteins_are_recognised(name):
    m = find_membrane(DATA / f"{name}.cif", load_backbone(DATA / f"{name}.cif", "asu"))
    assert m is not None and m.source == "estimate" and 13 <= m.half <= 17
    assert np.isclose(np.linalg.norm(m.normal), 1.0)


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "2HHB", "8UUP", "8UTF", "6ZS5", "1TIM", "AF-P04637"])
def test_soluble_proteins_are_not(name):
    assert find_membrane(DATA / f"{name}.cif", load_backbone(DATA / f"{name}.cif", "asu")) is None


def test_opm_dummy_atoms_give_the_exact_membrane(tmp_path):
    import gemmi

    st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    st.setup_entities()
    chain = gemmi.Chain("M")
    for k, (name, z) in enumerate((("N", 14.2), ("O", -14.2))):
        res = gemmi.Residue()
        res.name, res.seqid, res.het_flag = "DUM", gemmi.SeqId(900 + k, " "), "H"
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = name, gemmi.Element(name), gemmi.Position(0, 0, z)
        res.add_atom(atom)
        chain.add_residue(res)
    st[0].add_chain(chain)
    out = tmp_path / "opm.pdb"
    st.write_pdb(str(out))
    m = find_membrane(out, load_backbone(out, "asu"))
    assert m.source == "opm" and m.half == pytest.approx(14.2) and abs(m.normal[2]) == pytest.approx(1.0)


def _band(fig):
    return [a for a in fig.axes[0].patches if (a.get_gid() or "") == "membrane"]


def test_band_drawn_across_the_transmembrane_helices():
    lay, sses, _ = make_layout(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    assert lay.membrane is not None
    y0, y1 = lay.membrane["y"]
    assert y1 - y0 > 3.0
    long = [p for p in lay.placed.values() if p.sse.kind == "H" and len(p.sse) >= 18]
    crossing = [p for p in long if p.rect[1] < y0 + 0.5 and p.rect[3] > y1 - 0.5]
    assert len(crossing) >= 5  # bacteriorhodopsin's seven TM helices run through the band
    fig = make_figure(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    assert len(_band(fig)) == 1
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert {"in", "out"} <= texts


def test_positive_inside_rule_puts_bacteriorhodopsin_c_terminus_inside():
    from foldmap.layout import termini

    lay, sses, _ = make_layout(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    y0, y1 = lay.membrane["y"]
    ends = {end: port for end, _, port, _ in termini(lay, sses)}
    inside_low = lay.membrane["inside"] == "below"
    c_y, n_y = ends["C"][1], ends["N"][1]
    assert (c_y < y0) == inside_low and (n_y > y1) == inside_low  # N outside, C in the cytoplasm


def test_soluble_protein_with_membrane_auto_draws_nothing(tmp_path, capsys):
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(tmp_path / "x.svg"), "--membrane", "auto"]) == 0
    assert "membrane" in capsys.readouterr().err.lower()
    assert 'id="membrane"' not in (tmp_path / "x.svg").read_text()


def test_sides_labelled_only_when_the_positive_inside_rule_is_clear():
    fig = make_figure(DATA / "1BL8.cif", membrane="auto", assembly="asu")  # truncated KcsA: too close to call
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert len(_band(fig)) == 1 and not ({"in", "out"} & texts)


def test_outside_on_top_even_for_symmetric_assemblies():
    lay, _, _ = make_layout(DATA / "1C3W.cif", membrane="auto")  # bacteriorhodopsin trimer, unrolled about C3
    assert lay.membrane["sides"] and lay.membrane["inside"] == "below"
