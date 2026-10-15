"""Science checks (plan 12 Task 1): what a figure says about the structure, judged against its coordinates apart from
the page checks (figcheck A1-A20). S1: sheet neighbours pair; S2: a barrel's closure pairs; S3: a barrel's strand
count is its pairing cycle."""

from pathlib import Path

import pytest

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import Sheet
from foldmap.science import SCIENCE, neighbour_findings, pairing_cycle, science_lines, sheet_checks
from foldmap.sheets import build_sheets
from foldmap.ss import build_sses

DATA = Path(__file__).parent / "data"
R100 = DATA / "r100"


def counts(pairs: dict) -> callable:
    table = {frozenset(p): n for p, n in pairs.items()}
    return lambda a, b: table.get(frozenset((a, b)), 0)


def col(*ids, ghost=False):
    return [(i, ghost) for i in ids]


def test_a_strand_wedged_between_two_that_pair_fails_s1():
    count = counts({("a", "b"): 4, ("b", "c"): 5, ("b", "x"): 3})  # x pairs with b only, as 9FWZ's D:662-666
    s1, s2 = neighbour_findings([col("a"), col("b"), col("x"), col("c")], count)
    assert s1 == ["x | c drawn side by side with no bridge"] and s2 == []


def test_a_strand_stacked_beside_its_parent_passes_and_one_beside_none_fails():
    cols = [col("a"), col("b"), col("c", "x")]  # x stacked with c, b's other partner, beside b
    assert neighbour_findings(cols, counts({("a", "b"): 4, ("b", "c"): 5, ("b", "x"): 3})) == ([], [])
    s1, _ = neighbour_findings(cols, counts({("a", "b"): 4, ("b", "c"): 5, ("a", "x"): 3}))  # x pairs with a only
    assert s1 == ["x drawn in a shared column beside none of its partners"]
    s1, _ = neighbour_findings([col("a"), col("b", "x")], counts({("b", "x"): 3}))  # no bridge across at all
    assert s1 == [
        "columns a | b/x share no bridge",
        "b drawn in a shared column beside none of its partners",
        "x drawn in a shared column beside none of its partners",
    ]


def test_an_unpaired_closure_fails_s2():
    count = counts({("a", "b"): 4, ("b", "c"): 5})
    s1, s2 = neighbour_findings([col("a"), col("b"), col("c"), col("a", ghost=True)], count)
    assert s1 == [] and s2 == ["c | a (ghost) drawn side by side with no bridge"]


def test_the_pairing_cycle_leaves_a_pendant_strand_out():
    ring = [f"s{k}" for k in range(24)]
    pairs = {(ring[k], ring[(k + 1) % 24]): 6 for k in range(24)}
    pairs[("s19", "acc")] = 5  # a strand pairing with one barrel strand only
    assert pairing_cycle([*ring, "acc"], counts(pairs)) == 24
    pairs[(ring[23], ring[0])] = 1  # one closing bridge counts only as the drawn closure (1TIM's)
    assert pairing_cycle([*ring, "acc"], counts(pairs)) == 0
    assert pairing_cycle([*ring, "acc"], counts(pairs), [(ring[0], ring[23])]) == 24


def wedged(sheet: Sheet) -> Sheet:
    """The sheet as the review found it drawn: each accessory strand put back into the track right after its parent
    (the old sheets._order_tree); a sheet with none is returned as it is."""
    strands, dirs = list(sheet.strands), list(sheet.directions)
    for acc in sheet.accessory:
        k = strands.index(acc.parent) + 1
        strands.insert(k, acc.strand)
        dirs.insert(k, acc.direction)
    if not sheet.accessory:
        return sheet
    return Sheet(strands, dirs, ["?"] * (len(strands) - 1), sheet.closed, True)


def _built(path):
    bb = load_backbone(path, "asu")
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss, short_helices=True, bridges=d.bridges)
    return sses, build_sheets(sses, d.bridges), d.bridges


def test_fimd_drawn_with_its_accessory_strand_in_the_ring_fails_s1_and_s3():
    """9FWZ FimD (the review's R1): D:662-666 pairs only with D:645-659, but was drawn between D:645-659 and D:139-153,
    which pair with 12 bridges, making a 25-strand barrel of the 24-strand ring."""
    sses, sheets, bridges = _built(R100 / "9FWZ-DB.cif.gz")
    (barrel,) = [sh for sh in sheets if sh.closed]
    found = sheet_checks([wedged(barrel)], sses, bridges)
    ((s1,),) = [found["S1"]]  # wedged after D:645-659, beside a strand it does not pair with (D:139-153 as drawn)
    assert s1.startswith("sheet 1: D:662-666 | D:") and s1.endswith(" drawn side by side with no bridge")
    assert found["S2"] == []  # its closure, D:197-209 | D:214-223, pairs
    assert found["S3"] == ["sheet 1: 25 barrel strands drawn, the pairing cycle has 24"]
    assert found["S4"] is None and found["S5"] is None  # hooks, not yet measured


@pytest.mark.parametrize("path", [DATA / "1EMA.cif", DATA / "1TIM.cif", R100 / "1BXW-A.cif.gz", DATA / "1UBQ.cif"])
def test_clean_sheets_pass_every_science_check(path):
    sses, sheets, bridges = _built(path)
    found = sheet_checks(sheets, sses, bridges)
    assert found["S1"] == found["S2"] == found["S3"] == []


def _drawn(path):
    from foldmap.spec import draw_spec, figure_spec

    return draw_spec(figure_spec(path))


@pytest.mark.parametrize("path", [DATA / "1EMA.cif", R100 / "1BXW-A.cif.gz", DATA / "1TIM.cif"])
def test_a_control_barrel_figure_passes_and_a_moved_strand_fails(path):
    from foldmap.figcheck import drawn_columns, science_checks

    d = _drawn(path)
    assert science_checks(d) | {"S4": None, "S5": None} == {"S1": [], "S2": [], "S3": [], "S4": None, "S5": None}
    k = next(k for k, ids in enumerate(d.layout.sheets) if len(ids) >= 8)
    cols = drawn_columns(d.layout, k)
    a, b = d.layout.placed[cols[1][0][0]], d.layout.placed[cols[4][0][0]]
    a.cx, b.cx = b.cx, a.cx  # two strands of the ring swapped: neighbours that do not pair
    found = science_checks(d)
    assert found["S1"] and all(t.startswith(f"sheet {k}: ") for t in found["S1"])
    ghost = next(g for g in d.layout.ghosts if g.sse.id in d.layout.sheets[k])
    ghost.cx = d.layout.placed[cols[3][0][0]].cx + 0.5  # the ghost moved into the ring, beside a non-partner
    assert science_checks(d)["S2"]


def test_science_lines_say_pass_count_or_not_yet_measured():
    lines = science_lines({"S1": [], "S2": ["sheet 1: c | a (ghost) drawn side by side with no bridge"], "S3": []})
    assert lines[0] == "S1 sheet neighbours pair: pass"
    assert lines[1:3] == [
        "S2 barrel closure pairs: 1 problem(s)",
        "  sheet 1: c | a (ghost) drawn side by side with no bridge",
    ]
    assert lines[-2:] == [
        "S4 residue marks at their residue: not yet measured",
        "S5 automatic domains agree with a reference: not yet measured",
    ]
    assert list(SCIENCE) == ["S1", "S2", "S3", "S4", "S5"]


def test_summary_prints_the_science_checks():
    from foldmap.cli import summarize

    text = summarize(str(R100 / "1BXW-A.cif.gz"))
    assert "S1 sheet neighbours pair: pass" in text and "S3 barrel strand count is its pairing cycle: pass" in text


def test_check_command_prints_the_figures_science_and_exits_1_on_a_problem(capsys, monkeypatch):
    from foldmap import figcheck
    from foldmap.cli import main

    assert main(["check", str(DATA / "1UBQ.cif")]) == 0
    out = capsys.readouterr().out
    assert "S1 sheet neighbours pair: pass" in out and "S4 residue marks at their residue: not yet measured" in out
    monkeypatch.setattr(
        figcheck,
        "science_checks",
        lambda d: {
            "S1": ["sheet 0: a | b drawn side by side with no bridge"],
            "S2": [],
            "S3": [],
            "S4": None,
            "S5": None,
        },
    )
    assert main(["check", str(DATA / "1UBQ.cif")]) == 1
    assert "S1 sheet neighbours pair: 1 problem(s)" in capsys.readouterr().out


@pytest.mark.parametrize(
    "path",
    [
        R100 / "9FWZ-DB.cif.gz",
        R100 / "2OMF-A.cif.gz",
        R100 / "3L9G-A.cif.gz",
        R100 / "5KQJ.cif.gz",
        DATA / "8UTF.cif",
    ],
)
def test_branched_sheets_draw_every_strand_beside_a_partner(path):
    """Plan 12 Task 2: S1-S3 hold in the drawn figure of entries with branched sheets (9FWZ's and OmpF's barrels with
    a pendant strand, 3L9G's barrels with branches of several strands, 5KQJ's stacked pair, 8UTF's partners along one
    long strand), and every accessory strand stands in a column next to its parent's."""
    from foldmap.figcheck import drawn_columns, science_checks

    d = _drawn(path)
    found = science_checks(d)
    assert found["S1"] == found["S2"] == found["S3"] == [], found
    assert d.layout.accessory
    for k in range(len(d.layout.sheets)):
        where = {i: c for c, col in enumerate(drawn_columns(d.layout, k)) for i, g in col if not g}
        for strand, parent in d.layout.accessory.items():
            if strand in where and parent in where:
                assert abs(where[strand] - where[parent]) == 1, (strand, parent)


def test_fimd_figure_has_its_24_strand_barrel_and_says_where_its_pendant_strand_is(capsys, tmp_path):
    from foldmap.cli import main

    assert main(["plot", str(R100 / "9FWZ-DB.cif.gz"), "-o", str(tmp_path / "f.svg")]) == 0
    assert "D:662-666 drawn beside its partner D:645-659" in capsys.readouterr().err
    d = _drawn(R100 / "9FWZ-DB.cif.gz")
    (k,) = [k for k, ids in enumerate(d.layout.sheets) if len(ids) > 20]
    ring = [i for i in d.layout.sheets[k] if i not in d.layout.accessory]
    assert len(ring) == 24 and d.layout.accessory == {"D:662-666": "D:645-659"}
    assert any("D:662-666 drawn beside its partner D:645-659" in n for n in d.layout.sheet_notes)


def test_every_accessory_is_named_with_its_partner_in_the_notes():
    """Fix round 1: the note names each accessory strand with its partner, all of them (3L9G-A has 16)."""
    d = _drawn(R100 / "3L9G-A.cif.gz")
    (note,) = [n for n in d.layout.sheet_notes if "off a sheet's run" in n]
    assert note.startswith("16 strand(s) off a sheet's run of paired strands: ")
    assert all(f"{a} drawn beside its partner {p}" in note for a, p in d.layout.accessory.items())


def test_summary_names_the_accessory_strands():
    from foldmap.cli import summarize

    text = summarize(str(R100 / "9FWZ-DB.cif.gz"))
    assert "sheet 1: 24 strand(s) closed ambiguous" in text
    assert "  beside D:645-659 (off the run): D:662-666" in text


def test_ssb_draws_each_subunit_its_own_ob_fold_barrel():
    """1EYG (SSB tetramer): the four OB-fold barrels pair across subunits (A:1005-1021 with B:2005-2021). They were
    drawn as two 10-strand "barrels", each with two neighbours that do not pair; a part of a component closing a ring
    of its own is a sheet of its own: four closed 5-strand barrels, S1-S3 clean."""
    from foldmap.figcheck import science_checks

    sses, sheets, bridges = _built(DATA / "na" / "1EYG.cif.gz")
    barrels = [sh for sh in sheets if sh.closed]
    assert [len(sh.strands) for sh in barrels] == [5, 5, 5, 5]
    assert [sorted({s.chain for s in sh.strands}) for sh in barrels] == [["A"], ["B"], ["C"], ["D"]]
    d = _drawn(DATA / "na" / "1EYG.cif.gz")
    found = science_checks(d)
    assert found["S1"] == found["S2"] == found["S3"] == []
    assert "barrels drawn apart, each closed, that also pair: A:1005-1021 with B:2005-2021, " in d.layout.sheet_notes[0]
