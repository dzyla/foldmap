"""Conservation colouring of the topology figure from a multiple sequence alignment (ConSurf's 9 grades)."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap.cli import main, make_figure, make_layout
from foldmap.render import CONSURF, element_colours, grade_colour
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


def test_layout_carries_grades_and_honest_legend():
    lay, sses, _ = make_layout(UBQ, msa=FAMILY)
    assert len(lay.res_grade) == len(lay.res_chain) and set(lay.res_grade) <= set(range(10))
    assert lay.cons_label == "conservation (entropy-based, ConSurf colour scale)"
    fig = make_figure(UBQ, look=THEMES["conservation"], msa=FAMILY)
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert any("entropy-based, ConSurf colour scale" in t for t in texts)


def test_grade_colours():
    assert [grade_colour(g) for g in range(1, 10)] == list(CONSURF)
    assert to_hex(grade_colour(0)) == "#d9d9d9"


def consurf_table(tmp_path, grade_of):
    """A ConSurf grades table for ubiquitin chain A: grade_of(residue number) -> grade (synthesised in the
    documented consurf_grades.txt layout; no real file is available offline)."""
    import gemmi

    st = gemmi.read_structure(str(UBQ))
    rows = [" POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL", "    \t    \t        \t(normalized)"]
    for k, res in enumerate(st[0]["A"].get_polymer(), 1):
        one = gemmi.find_tabulated_residue(res.name).one_letter_code.upper()
        g = grade_of(res.seqid.num)
        flag = "*" if k == 1 else ""  # ConSurf's mark for a grade with too little data
        rows.append(f"{k:4d}\t   {one}\t{res.name}{res.seqid.num}:A\t{0.1 * (5 - g):.3f}\t\t  {g}{flag}\t-")
    rows.append(f"{len(rows):4d}\t   G\t      -   \t 0.125\t\t  4\t-")  # a position not in the structure
    path = tmp_path / "consurf_grades.txt"
    path.write_text("Amino Acid Conservation Scores\n=====\n\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_consurf_table_colours_by_its_own_grades(tmp_path):
    from foldmap.sequence import is_consurf, read_consurf

    table = consurf_table(tmp_path, lambda n: 9 if n <= 10 else 1)
    assert is_consurf(table) and not is_consurf(FAMILY)
    grades = read_consurf(table)
    assert grades[("A", 1, "")] == 9 and grades[("A", 70, "")] == 1 and len(grades) == 76
    lay, sses, _ = make_layout(UBQ, msa=table)
    assert lay.cons_label == "conservation (ConSurf grades)"
    colours, _ = element_colours(lay, sses, Style(color_by="conservation"))
    first = min(sses, key=lambda s: s.start)  # β1, residues 1-7: grade 9
    assert to_hex(colours[first.id]) == to_hex(CONSURF[8])


def test_consurf_table_is_not_an_alignment(tmp_path):
    with pytest.raises(ValueError, match="ConSurf"):
        read_alignment(consurf_table(tmp_path, lambda n: 5))
