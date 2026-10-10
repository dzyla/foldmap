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


def test_real_consurf_db_table_is_parsed():
    """A trimmed ConSurf-DB score file (header `ATOM`, atoms as RES:number:chain, legend and footer text)."""
    from foldmap.sequence import is_consurf, read_consurf

    table = DATA / "consurf" / "1TIM_A_grades.txt"
    assert is_consurf(table)
    grades = read_consurf(table)
    assert len(grades) == 12
    assert grades[("A", 1, "")] == 8  # '8*': below ConSurf's confidence cut-off, grade kept
    assert grades[("A", 4, "")] == 9 and grades[("A", 6, "")] == 4 and grades[("A", 13, "")] == 9
    assert ("A", 3, "") not in grades  # numbering follows the structure (POS 3 is residue 4)


def test_consurf_atoms_with_insertion_codes_and_missing_positions(tmp_path):
    from foldmap.sequence import is_consurf, read_consurf

    header = "POS\tSEQ\t     ATOM\t SCORE\tCOLOR\tCONFIDENCE INTERVAL\n"
    for atoms in (("ALA:52A:A", "ALA:-3:A"), ("ALA52A:A", "ALA-3:A")):
        path = tmp_path / "t.txt"
        path.write_text(
            "Legend: the COLOR and ATOM columns\n\n" + header + f"1\tA\t{atoms[0]}\t0.1\t 7*\t-\n"
            f"2\tA\t{atoms[1]}\t0.1\t 2\t-\n3\tG\t      -   \t0.1\t 4\t-\n",
            encoding="utf-8",
        )
        assert is_consurf(path)
        assert read_consurf(path) == {("A", 52, "A"): 7, ("A", -3, ""): 2}
    legend_only = tmp_path / "l.txt"
    legend_only.write_text("COLOR and ATOM are columns\n>a\nMKV\n", encoding="utf-8")
    assert not is_consurf(legend_only)


def test_consurf_grades_follow_identical_chains():
    """ConSurf-DB tables cover one chain of a homodimer: its identical partner takes the same grades."""
    lay, _, bb = make_layout(UBQ.parent / "1TIM.cif", msa=DATA / "consurf" / "1TIM_A_grades.txt")
    by = {(l.chain, l.seq): g for l, g in zip(bb.labels, lay.res_grade)}
    for n in (4, 6, 13):
        assert by[("A", n)] in (9, 4) and by[("B", n)] == by[("A", n)]
    assert by[("B", 100)] == 0  # beyond the trimmed fixture: no grade for either chain


def test_consurf_grades_skip_chains_with_a_different_sequence(tmp_path):
    import gemmi

    st = gemmi.read_structure(str(DATA / "1TIM.cif"))
    for res in st[0]["B"]:
        res.name = "GLY"  # chain B no longer has chain A's sequence
    mutated = tmp_path / "tim_b_gly.cif"
    st.make_mmcif_document().write_file(str(mutated))
    lay, _, bb = make_layout(mutated, msa=DATA / "consurf" / "1TIM_A_grades.txt")
    by = {(l.chain, l.seq): g for l, g in zip(bb.labels, lay.res_grade)}
    assert by[("A", 4)] == 9 and by[("B", 4)] == 0


def test_legend_does_not_widen_the_figure_beyond_the_protein():
    wide = make_figure(UBQ, look=THEMES["conservation"], msa=FAMILY)
    plain = make_figure(UBQ, look=THEMES["publication"])
    assert wide.get_size_inches()[0] <= plain.get_size_inches()[0] * 1.1


def test_long_legend_is_never_clipped():
    fig = make_figure(UBQ, look=THEMES["conservation"], msa=FAMILY)
    fig.canvas.draw()
    legend = [t for t in fig.axes[0].texts if "conservation" in t.get_text()]
    assert legend
    for t in legend:
        assert t.get_window_extent().x1 <= fig.bbox.x1
