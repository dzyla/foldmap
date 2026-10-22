"""The researcher's own data (plan 13d-4): per-residue values and crosslinks read against the structure, every row
that names no residue of it said with its row number (Review Focus 5)."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main, make_layout
from foldmap.model import Backbone, ResLabel
from foldmap.spec import figure_spec
from foldmap.userdata import ResidueData, read_crosslinks, read_residue_data, read_user_data

DATA = Path(__file__).parent / "data"
UBQ = DATA / "1UBQ.cif"


def _bb(*chains: tuple[str, list[tuple[str, str]]]) -> Backbone:
    """A backbone of the given chains, each [(residue tag, residue name)] (coordinates zero: only labels are read)."""
    labels = []
    for chain, residues in chains:
        for tag, name in residues:
            seq, icode = (int(tag[:-1]), tag[-1]) if tag[-1].isalpha() else (int(tag), "")
            labels.append(ResLabel(chain, seq, icode, name))
    n = len(labels)
    return Backbone(labels, np.zeros((n, 4, 3)), np.zeros(n, bool), {c for c, _ in chains})


CHAIN_A = [("1", "MET"), ("2", "GLN"), ("3", "ILE"), ("4", "PHE"), ("100", "GLY"), ("100A", "SER"), ("100B", "TYR")]
CHAIN_B = [("1", "ALA"), ("2", "GLY")]


def _csv(tmp_path, text: str, name: str = "data.csv") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _rows(problems) -> dict[str, str]:
    """{where: text} of the warnings."""
    assert not any(p.fatal for p in problems), problems
    return {p.where: p.text for p in problems}


def test_a_column_of_numbers_is_numbers_and_one_of_words_is_categories(tmp_path):
    bb = _bb(("A", CHAIN_A), ("B", CHAIN_B))
    path = _csv(tmp_path, "chain,residue,site,uptake\nA,1,core,0.5\nA,2,surface,NA\nB,2,core,-1.25\nA,3,,3\n")
    got, problems = read_residue_data(path, None, bb)
    assert problems == []  # site is words: the first column of numbers is uptake, the only one
    assert got == ResidueData("uptake", "number", {("A", "1"): 0.5, ("A", "3"): 3.0, ("B", "2"): -1.25})
    got, problems = read_residue_data(path, "SITE", bb)  # a column named in any case
    assert problems == []
    assert got == ResidueData("site", "category", {("A", "1"): "core", ("A", "2"): "surface", ("B", "2"): "core"})
    assert list(got.values) == [("A", "1"), ("A", "2"), ("B", "2")]  # in the structure's residue order


def test_several_columns_of_numbers_the_first_is_read_and_the_others_named(tmp_path):
    bb = _bb(("A", CHAIN_A))
    path = _csv(tmp_path, "chain\tresidue\tddG\tuptake\nA\t1\t0.5\t2\nA\t4\t1e-1\t3\n")  # tab separated
    got, problems = read_residue_data(path, None, bb)
    assert got.column == "ddG" and got.values == {("A", "1"): 0.5, ("A", "4"): 0.1}
    assert _rows(problems) == {
        f"--data {path}": "read column ddG, the first of numbers; the file also has uptake (--data-column NAME)"
    }


def test_insertion_codes_in_the_residue_cell_or_a_column_of_their_own(tmp_path):
    bb = _bb(("A", CHAIN_A))
    got, problems = read_residue_data(_csv(tmp_path, "chain,residue,v\nA,100,1\nA,100A,2\nA,100 B,3\n"), None, bb)
    assert problems == [] and got.values == {("A", "100"): 1.0, ("A", "100A"): 2.0, ("A", "100B"): 3.0}
    path = _csv(tmp_path, "chain,resnum,icode,v\nA,100,,1\nA,100,A,2\nA,100,?,4\nA,100,C,5\n", "icode.csv")
    got, problems = read_residue_data(path, None, bb)
    assert got.values == {("A", "100"): 1.0, ("A", "100A"): 2.0}
    assert _rows(problems) == {
        f"--data {path} row 4": "A 100 is given by row 2 already; the row is not used",  # ? is no insertion code
        f"--data {path} row 5": "chain A has no residue 100C (the structure file numbers its residues 1 to 100); "
        "the row is not used",
    }


def test_every_row_that_names_no_residue_is_said_with_its_row_number(tmp_path):
    bb = _bb(("A", CHAIN_A), ("B", CHAIN_B))
    text = "# HDX uptake, 30 s\nchain,residue,uptake\nA,1,0.1\nA,7,0.2\nC,1,0.3\n\nA,x1,0.4\nR,5,1\nA,1,0.9\nB,2,2\n"
    path = _csv(tmp_path, text)
    got, problems = read_residue_data(path, None, bb, na_chains=["R"])
    assert got.values == {("A", "1"): 0.1, ("B", "2"): 2.0}
    where = f"--data {path} row"
    assert _rows(problems) == {  # the file's line numbers: the comment is row 1, the header row 2
        f"{where} 4": "chain A has no residue 7 (the structure file numbers its residues 1 to 100); "
        "the row is not used",
        f"{where} 5": "no protein chain C in the structure; its protein chains are A, B; the row is not used",
        f"{where} 7": "'x1' is no residue number (45, -7 or 100A); the row is not used",
        f"{where} 8": "chain R is a nucleic acid; only protein residues are read; the row is not used",
        f"{where} 9": "A 1 is given by row 3 already; the row is not used",
    }


def test_a_residue_named_otherwise_than_in_the_structure_is_not_coloured(tmp_path):
    """A file numbered another way (here one off) is caught where its residue names differ from the structure's."""
    bb = _bb(("A", CHAIN_A))
    path = _csv(tmp_path, "chain,residue,aa,v\nA,1,M,1\nA,2,ile,2\nA,3,ILE,3\nA,4,X,4\n")
    got, problems = read_residue_data(path, None, bb)
    assert got.values == {("A", "1"): 1.0, ("A", "3"): 3.0}
    assert _rows(problems) == {
        f"--data {path} row 3": "A 2 is GLN in the structure, the row says ile; the row is not used",
        f"--data {path} row 5": "A 4 is PHE in the structure, the row says X; the row is not used",
    }


def test_an_assembly_copy_takes_its_deposited_chains_values(tmp_path):
    bb = _bb(("A1", CHAIN_A), ("B1", CHAIN_B), ("A2", CHAIN_A), ("B2", CHAIN_B))
    path = _csv(tmp_path, "chain,residue,v\nA,2,5\nA,100A,6\nB2,1,7\n")
    got, problems = read_residue_data(path, None, bb)
    assert got.values == {
        ("A1", "2"): 5.0,
        ("A1", "100A"): 6.0,
        ("A2", "2"): 5.0,
        ("A2", "100A"): 6.0,
        ("B2", "1"): 7.0,  # a copy named as drawn: that copy only
    }
    assert _rows(problems) == {
        f"--data {path}": "chain A is drawn as its assembly copies A1, A2; each copy takes chain A's rows"
    }


def test_an_assembly_built_from_the_file_maps_deposited_chains_to_every_copy(tmp_path):
    """1SLV under --assembly auto is drawn as A1, B1, A2, B2: a file of the deposited chain A colours A1 and A2."""
    from foldmap.io import load_backbone

    bb = load_backbone(DATA / "r100" / "1SLV.cif.gz", "auto")
    assert {"A1", "A2"} <= {l.chain for l in bb.labels} and "A" not in {l.chain for l in bb.labels}
    first = next(l for l in bb.labels if l.chain == "A1")
    path = _csv(tmp_path, f"chain,residue,resname,v\nA,{first.tag},{first.name},1.5\n")
    got, problems = read_residue_data(path, None, bb)
    assert got.values == {("A1", first.tag): 1.5, ("A2", first.tag): 1.5}
    assert [p.text for p in problems] == [
        "chain A is drawn as its assembly copies A1, A2; each copy takes chain A's rows"
    ]


def test_a_file_with_no_column_of_numbers_names_its_columns(tmp_path):
    bb = _bb(("A", CHAIN_A))
    path = _csv(tmp_path, "chain,residue,site\nA,1,core\n")
    got, problems = read_residue_data(path, None, bb)
    assert got.values == {} and len(problems) == 1 and problems[0].fatal
    assert str(problems[0]) == (
        f"--data {path}: no column of numbers; name the column to read with --data-column NAME; its columns are "
        "chain, residue, site"
    )
    _, problems = read_residue_data(path, "uptake", bb)
    assert [str(p) for p in problems if p.fatal] == [
        f"--data {path}: no value column 'uptake'; its columns are chain, residue, site"
    ]


@pytest.mark.parametrize(
    "text, said",
    [
        ("chain,v\nA,1\n", "no residue column (residue, resnum or resi); its columns are chain, v"),
        ("residue,v\n1,1\n", "no chain column, and the structure has protein chains A, B; its columns are residue, v"),
        (
            "chain,residue,v\nA,500,1\nC,1,2\n",
            "no row names a residue of the structure (residue numbers as the structure file writes them)",
        ),
        ("chain,residue,v\nA,1,\n", "no row gives a value in column v"),
        ("", "the file is empty; its first line is a header (chain,residue,...)"),
    ],
)
def test_a_file_that_gives_nothing_is_an_error(tmp_path, text, said):
    _, problems = read_residue_data(_csv(tmp_path, text), None, _bb(("A", CHAIN_A), ("B", CHAIN_B)))
    assert [p.text for p in problems if p.fatal] == [said]


def test_one_protein_chain_needs_no_chain_column(tmp_path):
    got, problems = read_residue_data(_csv(tmp_path, "residue;value\n2;1,5\n3;2\n"), "value", _bb(("A", CHAIN_A)))
    assert problems == [] and got.kind == "category"  # semicolons; 1,5 is no number: categories
    got, _ = read_residue_data(_csv(tmp_path, "resi;value\n2;-0.5\n"), None, _bb(("A", CHAIN_A)))
    assert got.values == {("A", "2"): -0.5}


def test_non_finite_numbers_are_said_not_read(tmp_path):
    path = _csv(tmp_path, "chain,residue,v\nA,1,inf\nA,2,1\n")
    got, problems = read_residue_data(path, None, _bb(("A", CHAIN_A)))
    assert got.kind == "number" and got.values == {("A", "2"): 1.0}
    assert _rows(problems) == {f"--data {path} row 2": "v 'inf' is not a finite number; the row is not used"}


def test_crosslinks_are_read_with_their_distances_and_unmatched_rows_named(tmp_path):
    bb = _bb(("A", CHAIN_A), ("B", CHAIN_B))
    text = "chain1,residue1,chain2,residue2,distance\nA,1,B,2,11.5\nA,100A,A,3,\nA,9,B,1,20\nA,2,B,1,far\n"
    path = _csv(tmp_path, text, "xl.csv")
    links, problems = read_crosslinks(path, bb)
    assert links == [(("A", "1"), ("B", "2"), 11.5), (("A", "100A"), ("A", "3"), None), (("A", "2"), ("B", "1"), None)]
    assert _rows(problems) == {
        f"--crosslinks {path} row 4": "chain A has no residue 9 (the structure file numbers its residues 1 to 100); "
        "the row is not used",
        f"--crosslinks {path} row 5": "distance 'far' is no number; the crosslink is read without one",
    }


def test_crosslinks_between_copies_join_the_copies_of_one_number(tmp_path):
    bb = _bb(("A1", CHAIN_A), ("B1", CHAIN_B), ("A2", CHAIN_A), ("B2", CHAIN_B), ("C", CHAIN_B))
    path = _csv(tmp_path, "chain_a,residue_a,chain_b,residue_b\nA,1,B,2\nA,3,C,1\n", "xl.csv")
    links, problems = read_crosslinks(path, bb)
    assert links == [
        (("A1", "1"), ("B1", "2"), None),
        (("A2", "1"), ("B2", "2"), None),
        (("A1", "3"), ("C", "1"), None),  # one end a chain as named: to every copy of the other
        (("A2", "3"), ("C", "1"), None),
    ]
    assert [p.text for p in problems] == [
        "chain A is drawn as its assembly copies A1, A2; each copy takes chain A's rows",
        "chain B is drawn as its assembly copies B1, B2; each copy takes chain B's rows",
    ]
    _, problems = read_crosslinks(_csv(tmp_path, "chain1,chain2,residue2\nA,B,1\n", "bad.csv"), bb)
    assert [p.text for p in problems if p.fatal] == [
        "no residue1 and residue2 columns; its columns are chain1, chain2, residue2"
    ]
    apart = _bb(("A1", CHAIN_A), ("A2", CHAIN_A), ("B3", CHAIN_B), ("C", CHAIN_B))
    path = _csv(tmp_path, "chain1,residue1,chain2,residue2\nA,1,B,2\nA,1,C,2\n", "apart.csv")
    links, problems = read_crosslinks(path, apart)
    assert links == [(("A1", "1"), ("C", "2"), None), (("A2", "1"), ("C", "2"), None)]
    assert {p.where: p.text for p in problems if p.where.endswith("row 2")} == {
        f"--crosslinks {path} row 2": "the assembly copies of its two chains share no number (A1-B1); "
        "the row is not used"
    }


def test_reading_is_deterministic(tmp_path):
    bb = _bb(("A1", CHAIN_A), ("A2", CHAIN_A))
    path = _csv(tmp_path, "chain,residue,v\nA,4,1\nA,1,2\nA,77,3\n")
    first = read_user_data(path, None, path, bb)  # the crosslinks reader refuses this file's header: said, not raised
    assert read_user_data(path, None, path, bb) == first
    assert list(first[0].values) == [("A1", "1"), ("A1", "4"), ("A2", "1"), ("A2", "4")]


def test_make_layout_keeps_what_was_read_and_says_each_row_not_used(tmp_path, capsys):
    path = _csv(tmp_path, "chain,residue,uptake\nA,1,0.5\nA,76,0.25\nA,77,1\nA,78,1\n")
    lay, _, _ = make_layout(UBQ, data=str(path))
    assert lay.residue_data == ResidueData("uptake", "number", {("A", "1"): 0.5, ("A", "76"): 0.25})
    assert lay.crosslinks == [] and [p.where for p in lay.data_problems] == [
        f"--data {path} row 4",
        f"--data {path} row 5",
    ]
    err = capsys.readouterr().err
    assert f"foldmap: warning: --data {path} row 4: chain A has no residue 77 (the structure file numbers" in err
    assert f"foldmap: --data {path}: column uptake (numbers) on 2 residues" in err
    plain, _, _ = make_layout(UBQ)  # without data: the Layout's defaults (none), nothing said
    assert (plain.residue_data, plain.crosslinks, plain.data_problems) == (None, [], [])
    assert capsys.readouterr().err == ""


def test_the_command_line_refuses_a_file_it_cannot_read_before_layout(tmp_path, capsys):
    out = tmp_path / "x.svg"
    words = _csv(tmp_path, "chain,residue,site\nA,1,core\n")
    assert main(["plot", str(UBQ), "-o", str(out), "--data", str(words)]) == 2
    err = capsys.readouterr().err
    assert err == (
        f"foldmap: error: --data {words}: no column of numbers; name the column to read with --data-column NAME; "
        "its columns are chain, residue, site\n"
    )
    assert not out.exists()
    assert main(["plot", str(UBQ), "-o", str(out), "--data", str(words), "--data-column", "site"]) == 0
    assert "column site (categories) on 1 residues" in capsys.readouterr().err
    assert main(["plot", str(UBQ), "-o", str(out), "--data-column", "site"]) == 2
    assert "--data-column needs --data FILE.csv" in capsys.readouterr().err
    assert main(["plot", str(UBQ), "-o", str(out), "--crosslinks", str(tmp_path / "none.csv")]) == 2
    assert capsys.readouterr().err == f"foldmap: error: --crosslinks {tmp_path / 'none.csv'}: no such file\n"


def test_a_layout_file_keeps_the_data_files_and_one_without_has_none(tmp_path, capsys):
    import yaml

    data = _csv(tmp_path, "chain,residue,uptake,ddG\nA,1,0.5,1\n")
    xl = _csv(tmp_path, "chain1,residue1,chain2,residue2\nA,1,A,76\n", "xl.csv")
    saved = tmp_path / "fig.yaml"
    args = ["plot", str(UBQ), "-o", str(tmp_path / "x.svg"), "--data", str(data), "--data-column", "ddG"]
    assert main([*args, "--crosslinks", str(xl), "--save-layout", str(saved)]) == 0
    layout = yaml.safe_load(saved.read_text(encoding="utf-8"))["layout"]
    assert (layout["data"], layout["data_column"], layout["crosslinks"]) == (str(data), "ddG", str(xl))
    spec = figure_spec(UBQ, layout_file=saved)
    assert (spec.data, spec.data_column, spec.crosslinks) == (str(data), "ddG", str(xl))
    assert spec.layout_kwargs()["data_column"] == "ddG"
    err = capsys.readouterr().err
    assert err.count("column ddG (numbers) on 1 residues") == 1  # read once for the figure and its layout file
    plain = tmp_path / "plain.yaml"
    assert main(["plot", str(UBQ), "-o", str(tmp_path / "y.svg"), "--save-layout", str(plain)]) == 0
    assert not {"data", "data_column", "crosslinks"} & set(yaml.safe_load(plain.read_text())["layout"])
    with pytest.raises(ValueError, match="data_column 'ddG' needs data"):
        figure_spec(UBQ, data_column="ddG")


def test_many_rows_not_used_are_said_up_to_twenty_then_counted(tmp_path, capsys):
    """13d-4 review: a large file with many residues the model lacks is said row by row up to SAID_ROWS, then 'and N
    more'; the layout keeps every one."""
    from foldmap.cli import SAID_ROWS

    rows = "".join(f"A,{n},1\n" for n in range(100, 125))  # 25 residues 1UBQ lacks
    path = _csv(tmp_path, "chain,residue,uptake\nA,1,0.5\n" + rows)
    lay, _, _ = make_layout(UBQ, data=str(path))
    assert len([p for p in lay.data_problems if " row " in p.where]) == 25
    err = capsys.readouterr().err
    assert err.count(f"foldmap: warning: --data {path} row ") == SAID_ROWS == 20
    assert f"foldmap: warning: --data {path}: and 5 more rows not used\n" in err


def test_values_on_chains_the_figure_leaves_out_are_noted(tmp_path, capsys):
    """13d-4 review: rows that match the structure but on a chain --chains (or --focus) leaves out are read and kept,
    and a note says they are not drawn."""
    hb = DATA / "2HHB.cif"
    path = _csv(tmp_path, "chain,residue,v\nA,10,1\nA,11,2\nB,10,3\n")
    xl = _csv(tmp_path, "chain1,residue1,chain2,residue2\nA,10,B,10\nA,10,A,20\n", "xl.csv")
    lay, _, _ = make_layout(hb, chains=["A"], data=str(path), crosslinks=str(xl), assembly="asu")
    notes = {p.where: p.text for p in lay.data_problems}
    assert notes[f"--data {path}"] == (
        "chain B's values (1 residue) are not drawn: the figure leaves chain B out (--chains, --focus)"
    )
    assert [p.text for p in lay.data_problems if p.where == f"--crosslinks {xl}"] == [
        "1 crosslink reaches a chain the figure leaves out (B; --chains, --focus)",
        "crosslink A:10 – B:10: Cα–Cα 32.4 Å, longer than its 30 Å span",  # measured all the same; not drawn
    ]
    assert f"foldmap: warning: --data {path}: chain B's values" in capsys.readouterr().err
    lay, _, _ = make_layout(hb, data=str(path), assembly="asu")  # every chain drawn: nothing to note
    assert lay.data_problems == []


def test_a_crafted_table_is_refused_without_filling_memory(tmp_path):
    """Section review 13b-13d, Critical 1: a few-KB table with thousands of columns and one-cell rows was padded to
    the header's width row by row (4k x 4k: 123 MB from 30 KB). Wide headers, too many rows, a broken quote and an
    unreadable path are each one fatal Problem, never a traceback or a large allocation."""
    from foldmap import userdata

    wide = tmp_path / "wide.csv"
    wide.write_text(",".join(f"c{k}" for k in range(4000)) + "\n" + "x\n" * 4000)
    header, rows, problems = userdata._rows(wide, "wide.csv")
    assert (header, rows) == ([], []) and problems and problems[0].fatal and "columns" in problems[0].text

    long = tmp_path / "long.csv"
    long.write_text("chain,residue,value\n" + "A,1,0\n" * (userdata.MAX_ROWS + 1))
    header, rows, problems = userdata._rows(long, "long.csv")
    assert rows == [] and problems[0].fatal and "rows" in problems[0].text

    quote = tmp_path / "quote.csv"
    quote.write_text('chain,residue,value\nA,"1,' + "x" * 200_000 + "\n")
    header, rows, problems = userdata._rows(quote, "quote.csv")
    assert rows == [] and problems[0].fatal

    header, rows, problems = userdata._rows(tmp_path, "a directory")
    assert rows == [] and problems[0].fatal

    short = tmp_path / "short.csv"
    short.write_text("chain,residue,value,note\nA,5,1.5\n")
    header, rows, problems = userdata._rows(short, "short.csv")
    assert problems == [] and rows == [(2, ["A", "5", "1.5", ""])]
