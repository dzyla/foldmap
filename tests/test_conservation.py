"""Conservation colouring of the topology figure from a multiple sequence alignment (ConSurf's 9 grades)."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap.cli import main, make_figure, make_layout
from foldmap.render import CONSURF, element_colours
from foldmap.sequence import read_alignment, residue_conservation
from foldmap.style import THEMES, Style

DATA = Path(__file__).parent / "data"
UBQ, FAMILY = DATA / "1UBQ.cif", DATA / "ubq_family.fasta"


def test_scores_land_on_structure_residues():
    scores = residue_conservation(UBQ, read_alignment(FAMILY))
    assert len(scores) == 76 and all(0.0 <= v <= 1.0 for v in scores.values())
    lay, _, bb = make_layout(UBQ)
    first_gly = next(k for k, l in enumerate(bb.labels) if l.seq == 75)  # Gly75: invariant
    assert scores[first_gly] == pytest.approx(1.0)


def test_unrelated_chains_get_no_scores():
    assert residue_conservation(DATA / "2LZM.cif", read_alignment(FAMILY)) == {}


def test_elements_coloured_in_consurf_grades():
    lay, sses, _ = make_layout(UBQ, msa=FAMILY)
    assert len(lay.res_cons) == len(lay.res_chain)
    colours, residue = element_colours(lay, sses, Style(color_by="conservation"))
    grades = {to_hex(c) for c in CONSURF}
    assert {to_hex(c) for c in colours.values()} <= grades
    means = {s.id: np.nanmean(lay.res_cons[s.start : s.end + 1]) for s in sses}
    lo, hi = min(means, key=means.get), max(means, key=means.get)
    assert CONSURF.index(colours[hi]) > CONSURF.index(colours[lo])  # more conserved: further toward maroon


def test_conservation_needs_an_alignment(tmp_path, capsys):
    assert main(["plot", str(UBQ), "-o", str(tmp_path / "x.svg"), "--theme", "conservation"]) == 1
    assert "--msa" in capsys.readouterr().err
    assert main(["plot", str(UBQ), "-o", str(tmp_path / "x.svg"), "--theme", "conservation", "--msa", str(FAMILY)]) == 0


def test_theme_and_legend():
    assert THEMES["conservation"].color_by == "conservation"
    fig = make_figure(UBQ, look=THEMES["conservation"], msa=FAMILY)
    assert any("conservation" in t.get_text() for t in fig.axes[0].texts)


def test_layout_file_keeps_the_alignment(tmp_path):
    import yaml

    lf = tmp_path / "f.yaml"
    assert (
        main(
            [
                "plot",
                str(UBQ),
                "-o",
                str(tmp_path / "a.svg"),
                "--msa",
                str(FAMILY),
                "--theme",
                "conservation",
                "--save-layout",
                str(lf),
            ]
        )
        == 0
    )
    assert yaml.safe_load(lf.read_text())["layout"]["msa"] == str(FAMILY)
    assert main(["plot", str(UBQ), "-o", str(tmp_path / "b.svg"), "--layout-file", str(lf)]) == 0


def test_helix_labels_stay_readable_on_pale_grades():
    from matplotlib.colors import to_rgb

    fig = make_figure(UBQ, look=THEMES["conservation"], msa=FAMILY)
    for t in fig.axes[0].texts:
        if (t.get_gid() or "").startswith("label:") and t.get_text()[:1] in "αη":
            r, g, b = to_rgb(t.get_color())
            assert 0.2126 * r + 0.7152 * g + 0.0722 * b <= 0.5, t.get_text()
