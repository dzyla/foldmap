"""Science checks (plan 12 Task 1): what a figure says about the structure, judged against its coordinates apart from
the page checks (figcheck A1-A20). S1: sheet neighbours pair; S2: a barrel's closure pairs; S3: a barrel's strand
count is its pairing cycle."""

from pathlib import Path

import numpy as np
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


def row(*ids, ghost=()):
    """Strands side by side on one row (spans -5..5), one per column from 0."""
    return [(i, k, (-5.0, 5.0), i in ghost) for k, i in enumerate(ids)]


def spot(i, c, lo=None, hi=None, ghost=False):
    return (i, c, None if lo is None else (lo, hi), ghost)


def test_a_strand_wedged_between_two_that_pair_fails_s1():
    count = counts({("a", "b"): 4, ("b", "c"): 5, ("b", "x"): 3})  # x pairs with b only, as 9FWZ's D:662-666
    s1, s2 = neighbour_findings(row("a", "b", "x", "c"), count)
    assert s1 == ["x | c drawn side by side with no bridge"] and s2 == []


def test_s1_reads_spans_strands_beside_each_other_must_pair_and_one_column_is_collinear():
    """Fix round 3 ruling: any two strands in adjacent columns whose drawn spans overlap must share a bridge; two in
    one column are collinear (end to end) and must not overlap."""
    count = counts({("a", "b"): 4, ("b", "c"): 5, ("b", "x"): 3})
    track = row("a", "b", "c")
    assert neighbour_findings([*track, spot("x", 2, 6, 9)], count) == ([], [])  # end to end with c, level with b
    s1, _ = neighbour_findings([*track, spot("x", 2, 6, 9), spot("y", 3, 6, 9)], count)
    assert s1 == ["x | y drawn side by side with no bridge"]  # beside x at its height: 0 bridges
    assert neighbour_findings([*track, spot("y", 3, 6, 9)], count)[0] == []  # beside nothing at its height
    s1, _ = neighbour_findings([*track, spot("x", 2, 4, 9)], count)
    assert s1 == ["c and x overlap in one column"]


def test_an_unpaired_closure_fails_s2():
    count = counts({("a", "b"): 4, ("b", "c"): 5})
    s1, s2 = neighbour_findings([*row("a", "b", "c"), ("a", 3, (-5.0, 5.0), True)], count)
    assert s1 == [] and s2 == ["c | a (ghost) drawn side by side with no bridge"]


def test_8utf_a_strand_beside_an_unpaired_ghost_fails_s1():
    """Fix round 3 (review): 8UTF's a:266-273, stacked above a:278-301 in its column, stood beside the ghost of
    c:440-442 (spans overlapping by 1.14, 0 bridges); the column reading passed it. S1 must flag any such pair: the
    figure is clean."""
    from foldmap.figcheck import science_checks

    found = science_checks(_drawn(DATA / "8UTF.cif"))
    assert found["S1"] == found["S2"] == [], found


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
    assert found["S4"] is None and found["S5"] is None  # S4 needs a drawn figure, S5 a reference: not measured


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
    assert science_checks(d) == {"S1": [], "S2": [], "S3": [], "S4": [], "S5": None, "S6": []}
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
    assert lines[-3:] == [
        "S4 residue marks at their residue: not yet measured",
        "S5 automatic domains agree with a reference: not yet measured",
        "S6 accessory drawn at its partner: not yet measured",
    ]
    assert list(SCIENCE) == ["S1", "S2", "S3", "S4", "S5", "S6"]


def test_summary_prints_the_science_checks():
    from foldmap.cli import summarize

    text = summarize(str(R100 / "1BXW-A.cif.gz"))
    assert "S1 sheet neighbours pair: pass" in text and "S3 barrel strand count is its pairing cycle: pass" in text
    # final review B, Important 3: a check with nothing to judge in the sheets as built is not judged, never a pass
    assert "S6 accessory drawn at its partner: not judged (no accessory strand)" in text


def test_science_lines_say_not_judged_with_nothing_to_judge():
    lines = science_lines({"S1": [], "S2": [], "S3": [], "S4": None, "S5": None, "S6": []}, {"S1": 2, "S2": 0, "S6": 0})
    assert lines[0] == "S1 sheet neighbours pair: pass"
    assert lines[1] == "S2 barrel closure pairs: not judged (no closed barrel drawn)"
    assert lines[5] == "S6 accessory drawn at its partner: not judged (no accessory strand)"


def test_check_command_prints_the_figures_science_and_exits_1_on_a_problem(capsys, monkeypatch):
    from foldmap import figcheck
    from foldmap.cli import main

    assert main(["check", str(DATA / "1UBQ.cif")]) == 0
    out = capsys.readouterr().out
    assert "S1 sheet neighbours pair: pass" in out  # 1UBQ draws no residue mark: S4 has nothing to judge
    assert "S4 residue marks at their residue: not judged (no residue mark drawn)" in out
    monkeypatch.setattr(
        figcheck,
        "science_checks",
        lambda d: {
            "S1": ["sheet 0: a | b drawn side by side with no bridge"],
            "S2": [],
            "S3": [],
            "S4": None,
            "S5": None,
            "S6": [],
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
    assert found["S1"] == found["S2"] == found["S3"] == found["S6"] == [], found
    assert d.layout.accessory
    for k in range(len(d.layout.sheets)):  # fix round 1 ruling: in the partner's own column, or beside it alone
        where = {i: c for c, col in enumerate(drawn_columns(d.layout, k)) for i, g in col if not g}
        for strand, parent in d.layout.accessory.items():
            if strand in where and parent in where:
                assert abs(where[strand] - where[parent]) <= 1, (strand, parent)


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
    away = dict(d.layout.unshown)
    assert all(
        f"{a} {'also pairs with' if away.get(a) == p else 'drawn beside its partner'} {p}" in note
        for a, p in d.layout.accessory.items()
    )


def test_summary_names_the_accessory_strands():
    from foldmap.cli import summarize

    text = summarize(str(R100 / "9FWZ-DB.cif.gz"))
    assert "sheet 1: 24 strand(s) closed ambiguous" in text
    assert "  D:662-666 drawn beside its partner D:645-659" in text  # worded as the figure's note (accessory_lines)


@pytest.mark.parametrize("name", ["2OMF-A", "3L9G-A"])
def test_summary_words_each_accessory_as_the_figure_draws_it(name):
    """Fix round 4 (the re-review's N2): summary said "drawn beside its partner" for every accessory, also those the
    figure draws away from their partner and names (OmpF's three A:2-6, 8 of 3L9G's 16). It now reads the default
    figure: a level one is "drawn beside its partner Y", the others "pairs with Y (off the run)"."""
    from foldmap.cli import summarize

    lay = _drawn(R100 / f"{name}.cif.gz").layout
    away = set(lay.unshown)
    assert away and len(away) < len(lay.accessory) or name == "2OMF-A"
    text = summarize(str(R100 / f"{name}.cif.gz"))
    for a, p in lay.accessory.items():
        if (a, p) in away:
            assert f"  {a} pairs with {p} (off the run)" in text and f"{a} drawn beside" not in text
        else:
            assert f"  {a} drawn beside its partner {p}" in text and f"{a} pairs with" not in text


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


def test_s6_an_accessory_stands_beside_its_own_partner():
    """Fix round 3 ruling: an accessory strand stands in a column next to its partner's, level with it (spans
    overlapping); one in its partner's own column is collinear with it (end to end) and does not read as pairing."""
    from foldmap.science import partner_findings

    track = [spot("p", 1, -5, 5), spot("q", 2, -5, 5)]
    assert partner_findings([*track, spot("x", 2, 6, 8)], {"x": "p"}) == [
        "x drawn beside it, not level with it, not beside its partner p"
    ]
    assert partner_findings([*track, spot("x", 2, 3, 8)], {"x": "p"}) == []  # q's column, level with p
    assert partner_findings([*track, spot("x", 0, -2, 2)], {"x": "p"}) == []  # a new column beside p
    assert partner_findings([*track, spot("x", 1, 6, 8)], {"x": "p"}) == [
        "x drawn in its partner's column, not beside its partner p"
    ]
    assert partner_findings([*track, spot("x", 3, 6, 8)], {"x": "p"}) == [
        "x drawn away from it, not beside its partner p"
    ]
    assert partner_findings([*track, spot("x", 2)], {"x": "p"}) == []  # by column only, when the height is not known


def test_fimd_pendant_strand_stands_at_its_partner():
    """Fix round 3 (the review): S31, D:662-666, lies beside its partner S30 (D:645-659) on S29's side (D:634-641),
    8.6 Å past S29's end: drawn in S29's column, end to end with it, level with S30. S6 holds for the sheets as built
    and for the figure."""
    from foldmap.figcheck import science_checks

    sses, sheets, bridges = _built(R100 / "9FWZ-DB.cif.gz")
    assert sheet_checks(sheets, sses, bridges)["S6"] == []
    assert science_checks(_drawn(R100 / "9FWZ-DB.cif.gz"))["S6"] == []


@pytest.mark.parametrize("path", [R100 / "2OMF-A.cif.gz", DATA / "8UTF.cif"])
def test_an_accessory_stands_on_the_same_side_of_its_partner_in_every_copy(path):
    """Fix round 2: copies draw an accessory strand alike: towards the same end of its partner (seen along the
    partner's own direction) in every copy (2OMF's three A:2-6, 8UTF's three 266-273 and 227-228)."""
    lay = _drawn(path).layout
    sides: dict[tuple, set] = {}
    for a, p in lay.accessory.items():
        acc, par = lay.placed[a], lay.placed[p]
        key = (acc.sse.first.seq, par.sse.first.seq)
        sides.setdefault(key, set()).add(np.sign((acc.cy - par.cy) * par.direction[1]))
    assert sides and all(len(s) == 1 for s in sides.values()), sides


def test_s6_a_pairing_not_drawn_passes_when_it_is_named():
    """Fix round 2/3 rulings: an accessory that cannot stand level with its partner (no height free beside it) is
    drawn as near as it can, and its pairing is named (Layout.unshown, in the notes)."""
    from foldmap.science import partner_findings

    spots = [spot("p", 1, -5, 5), spot("q", 2, -5, 5), spot("x", 2, 6, 8)]
    assert partner_findings(spots, {"x": "p"}) == ["x drawn beside it, not level with it, not beside its partner p"]
    assert partner_findings(spots, {"x": "p"}, {("x", "p")}) == []


def _branch_sheet():
    """Seven strands in a row (t0-t6) and a branch off t3: a pairs with t3, c1 and c2 with a only (as 7ARD's
    x:79-86 with x:97-107, x:58-59 and x:67-68): a has three partners, two ends."""
    from foldmap.model import Bridge
    from helpers import fake_sse

    def strand(i, x, y, up=True):
        return fake_sse("E", 20 * i, 20 * i + 6, (x, y, 0.0), (0.0, 1.0 if up else -1.0, 0.0))

    track = [strand(i, 4.8 * i, 0.0, i % 2 == 0) for i in range(7)]
    a, c1, c2 = strand(7, 4.8 * 3, 25.0, False), strand(8, 4.8 * 3, 50.0), strand(9, 4.8 * 4, 25.0)

    def ladder(x, y):
        return [Bridge(x.start + k, y.start + k, "A") for k in range(3)]

    bridges = [b for p, q in zip(track, track[1:]) for b in ladder(p, q)]
    bridges += ladder(track[3], a) + ladder(a, c1) + ladder(a, c2)
    return [*track, a, c1, c2], bridges


def test_a_strand_with_more_accessories_than_ends_names_the_one_drawn_away():
    """Fix round 2/3 (7ARD): a has three partners. Each accessory stands in the column beside its partner, level with
    it where a height is free there; one that cannot is drawn as near as it can and named (Layout.unshown, "also
    pairs with"). S1 (spans) and S6 pass on the drawn block."""
    from foldmap import figcheck
    from foldmap.frame import Frame
    from foldmap.layout import build_layout
    from foldmap.science import bridge_counts, neighbour_findings, partner_findings

    sses, bridges = _branch_sheet()
    sheets = build_sheets(sses, bridges)
    (sh,) = sheets
    assert [x.strand.start for x in sh.accessory] == [140, 160, 180]  # a, then c1 and c2 hanging from a
    frame = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    lay = build_layout(sses, sheets, frame, attract=False)
    (k,) = range(len(lay.sheets))
    spots = [
        (i, c, figcheck._span(lay.placed[i]), g) for c, col in enumerate(figcheck.drawn_columns(lay, k)) for i, g in col
    ]
    assert neighbour_findings(spots, bridge_counts(sses, bridges)) == ([], [])
    assert partner_findings(spots, lay.accessory, set(lay.unshown)) == []
    (note,) = [n for n in lay.sheet_notes if "off a sheet's run" in n]
    assert all(f"{x} also pairs with {p}" in note for x, p in lay.unshown)


def test_s3_flags_a_ring_of_five_drawn_open():
    """Fix round 3 (the review): S3 ran on closed sheets only; a sheet drawn open whose pairing closes a ring of 5 or
    more strands is a barrel drawn open (8C3A's 5-ring was)."""
    from foldmap.science import open_finding

    ring = [f"s{k}" for k in range(5)]
    pairs = {(ring[k], ring[(k + 1) % 5]): 4 for k in range(5)}
    assert open_finding("sheet 1", ring, counts(pairs)) == "sheet 1: drawn open, its pairing closes a ring of 5 strands"
    del pairs[(ring[4], ring[0])]
    assert open_finding("sheet 1", ring, counts(pairs)) is None


def test_check_help_names_s6():
    from foldmap.cli import build_parser

    sub = next(a for a in build_parser()._actions if a.dest == "command")
    assert "S6" in sub.choices["check"].description


def test_s4_flags_a_mark_far_along_a_long_loop_from_its_residue():
    from foldmap.science import mark_findings

    loop = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0)]  # a 20-unit route; the residue is 2 past its first element
    near = ("A:66", (0.0, 1.0), (loop, 0, 2 * 0.5, 2 * 0.5))
    far = ("A:66", (5.0, 10.0), (loop, 0, 2 * 0.5, 2 * 0.5))
    off = ("A:66", (3.0, 3.0), (loop, 0, 2 * 0.5, 2 * 0.5))
    assert mark_findings([near]) == []
    assert "along its route" in mark_findings([far])[0]
    assert "off its own place" in mark_findings([off])[0]
    assert mark_findings([("A:66", (10.0, 9.5), (loop, -1, 0.5, 0.5))]) == []  # measured from the far end


def test_s4_flags_a_mark_at_the_element_end_for_a_residue_well_into_the_loop():
    """T4 review item 2: S4 bounded a mark's arc only from above, so a mark left at the element end passed for any
    residue. The place now carries the arc its residue distance puts the mark at, and a mark short of it by more
    than MARK_SLACK is flagged: here the residue is 6 into a 20-unit route, expected 3 units along it."""
    from foldmap.science import mark_findings

    loop = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0)]
    assert mark_findings([("A:70", (0.0, 3.0), (loop, 0, 6 * 0.5, 6 * 0.5))]) == []
    got = mark_findings([("A:70", (0.0, 0.0), (loop, 0, 6 * 0.5, 6 * 0.5))])
    assert len(got) == 1 and "short of" in got[0], got


def test_s4_places_a_loop_residue_by_its_residue_distance_on_a_drawn_figure():
    """figcheck.residue_places carries the expected arc: on 1EMA's drawn loops, every loop residue's mark at
    render.residue_point passes S4, and the same mark at its nearer element's end fails once the residue is two or
    more residues in (expected arc 1.0 or more against MARK_SLACK 0.5), unless the route is too short to hold it."""
    from foldmap.figcheck import residue_places
    from foldmap.render import residue_point
    from foldmap.science import MARK_SLACK, mark_findings

    d = _drawn(DATA / "1EMA.cif")
    chains = [lab.chain for lab in d.bb.labels]
    ids = {s.id: s for s in d.sses}
    flagged = 0
    for lp in d.loops:
        a, b = ids.get(lp.a_id), ids.get(lp.b_id)
        if a is None or b is None or a.chain != b.chain:
            continue
        for k in range(a.end + 1, b.start):
            place = residue_places(d, k)
            pts, end, _, expected = place
            assert mark_findings([(str(k), residue_point(d.layout, d.loops, d.sses, chains, k), place)]) == []
            at_end = tuple(pts[0] if end == 0 else pts[-1])
            got = mark_findings([(str(k), at_end, place)])
            assert bool(got) == (expected > MARK_SLACK), (k, expected, got)
            flagged += bool(got)
    assert flagged > 10


def test_1ema_chromophore_star_passes_s4_and_the_old_fraction_would_not(monkeypatch):
    from foldmap import render
    from foldmap.figcheck import science_checks

    d = _drawn(DATA / "1EMA.cif")
    assert science_checks(d)["S4"] == []
    # the pre-R3 rule: a fraction of the routed loop's length, 8 units along a 24-unit route for CRO66
    real = render.residue_point

    def by_fraction(layout, loops, sses, bb, r):
        at = real(layout, loops, sses, bb, r)
        lp = next((lp for lp in loops if lp.a_id == "A:57-64" and lp.b_id == "A:69-71"), None)
        if lp is None or r != 63:
            return at
        pts = lp.points
        seg = [((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5 for a, b in zip(pts, pts[1:])]
        target, k = sum(seg) / 3, 0
        while target > seg[k]:
            target -= seg[k]
            k += 1
        t = target / seg[k]
        return tuple(pts[k][i] + t * (pts[k + 1][i] - pts[k][i]) for i in (0, 1))

    monkeypatch.setattr(render, "residue_point", by_fraction)
    assert science_checks(_drawn(DATA / "1EMA.cif"))["S4"]


def _unshown_case(ox, ghost_at=None):
    """An accessory x (parent p, also pairing with o), strands standing vertically (angle pi/2, length 4, width 1):
    x at column 0, p at column 5 (apart), o at column `ox`; optionally a ghost of p at column `ghost_at`."""
    from tests.test_sheets import fake

    from foldmap.layout import Placed, _unshown, accessory_lines
    from foldmap.model import Accessory

    x, p, o = fake(0, 9), fake(20, 29), fake(40, 49)
    sh = Sheet([p, o], [1, 1], ["P"], accessory=[Accessory(x, p, 0, 1, (o,))])

    def at(s, cx, ghost=False):
        return Placed(s, cx, 0.0, 4.0, np.pi / 2, 1.0, ghost=ghost)

    placed = {x.id: at(x, 0.0), p.id: at(p, 5.0), o.id: at(o, ox)}
    ghosts = [at(p, ghost_at, True)] if ghost_at is not None else []
    return sh, x, p, o, _unshown(placed, [sh], ghosts), accessory_lines


def test_a_partners_ghost_beside_the_accessory_counts_as_shown():
    """Final review 3 (8UTF c:440-442): the partner drawn as a ghost level beside the accessory is drawn beside."""
    _, x, p, o, away, _ = _unshown_case(ox=1.0)
    assert (x.id, p.id) in away and (x.id, o.id) not in away
    _, x, p, o, away, _ = _unshown_case(ox=1.0, ghost_at=1.0)
    assert away == []
    d = _drawn(DATA / "8UTF.cif")
    (note,) = [n for n in d.layout.sheet_notes if "c:440-442" in n]
    assert "c:440-442 drawn beside its partner A:35-57" in note


def test_an_accessory_pairing_with_one_of_its_others_is_named_when_not_drawn_beside():
    """Final review 4: Accessory.others pairings that are not side by side are carried into unshown and named."""
    sh, x, p, o, away, lines = _unshown_case(ox=8.0)
    assert sorted(away) == sorted([(x.id, p.id), (x.id, o.id)])
    assert f"{x.id} also pairs with {o.id}" in lines(sh, away)
    sh, x, p, o, away, lines = _unshown_case(ox=1.0)
    assert (x.id, o.id) not in away and f"{x.id} also pairs with {o.id}" not in lines(sh, away)


def test_s3_reports_an_open_sheet_whose_cycles_cannot_be_counted(monkeypatch):
    """Final review 6: pairing_cycle gives up (None) past _MAX_CYCLES; an open sheet reports it, as a barrel does."""
    from foldmap import science

    ring = [f"s{k}" for k in range(5)]
    pairs = {(ring[k], ring[(k + 1) % 5]): 4 for k in range(5)}
    monkeypatch.setattr(science, "_MAX_CYCLES", 0)
    assert science.open_finding("sheet 1", ring, counts(pairs)) == "sheet 1: too many pairing cycles to count"


def test_s4_reports_a_mark_with_no_place_found():
    """Final review 5: a drawn mark whose residue has no place to be judged against is a finding, not skipped."""
    from foldmap.science import mark_findings

    assert mark_findings([("A:66", (1.0, 1.0), None)]) == ["A:66 is drawn, but no place was found for its residue"]
