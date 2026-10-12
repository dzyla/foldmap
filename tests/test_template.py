"""The βα template: all-parallel α/β sheets drawn as one row, each crossover helix between its two strands."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_layout
from foldmap.route import route_loops
from foldmap.style import Style
from helpers import FID, loop_crossings, pipeline

DATA = Path(__file__).parent / "data"


def _row(lay, sses, chain):
    sheet = next(b for b in lay.sheets if b and b[0].startswith(chain + ":") and len(b) == 8)
    strands = [lay.placed[k] for k in sheet]
    return strands, sorted(strands, key=lambda p: p.cx)


def test_tim_barrel_is_one_row_with_each_helix_between_its_strands():
    """The (βα)8 row of the textbook TIM figures: 2 loop crossings, against 30 for the projected figure (b43695e). The
    row is wider than the projected figure (27.2 x 13.0 = 353 units², against 20.3 x 14.2 = 288): eight strands, the
    ghost and seven helices, each helix with its name beside it, side by side. The crossings, not the area, decide
    whether the figure reads (Task 12 ruling: the area limit was a target). Opening the barrel between β8 and β1
    (the textbook cut) moved the free elements: 353.10 measured, hence 354."""
    lay, sses, bb = make_layout(DATA / "1TIM.cif", chains=["A"])
    strands, by_x = _row(lay, sses, "A")
    ys = {round(p.cy, 6) for p in strands}
    assert len(ys) == 1  # one row
    order = {p.sse.id: i for i, p in enumerate(by_x)}
    seq = sorted(strands, key=lambda p: p.sse.start)
    inrow = 0
    for a, b in zip(seq, seq[1:]):
        between = [p for p in lay.placed.values() if p.sse.kind == "H" and a.sse.end < p.sse.start < b.sse.start]
        assert between
        if abs(order[a.sse.id] - order[b.sse.id]) == 1:
            lo, hi = sorted((a.cx, b.cx))
            assert any(lo < h.cx < hi and abs(h.cy - a.cy) < 1e-6 for h in between)  # in the gap, in the row
            assert any(np.cos(h.angle - a.angle) < -0.99 for h in between)  # antiparallel to the strands
            inrow += 1
    assert inrow == 7  # the barrel opens between β8 and β1, which no unit joins: every unit joins neighbours
    (ghost,) = [g for g in lay.ghosts if g.sse.chain == "A"]
    assert ghost.sse in (by_x[0].sse, by_x[-1].sse) and {by_x[0].sse.start, by_x[-1].sse.start} == {
        min(p.sse.start for p in strands),
        max(p.sse.start for p in strands),
    }
    assert [k for row in lay.ba_rows for k in row if lay.placed[k].sse.kind == "H"]  # Task 13 numbers these
    loops = route_loops(lay, sses, bb)
    assert loop_crossings(loops) <= 4  # 30 at b43695e
    assert lay.width * lay.height <= 354.0  # measured 353.10 with the β8|β1 cut (x1.22 the projected 288 units²)


def test_7tim_barrel_row_keeps_a_unit_without_helix():
    lay, sses, bb = make_layout(FID / "7TIM.cif.gz", chains=["A"])
    strands, _ = _row(lay, sses, "A")
    assert len({round(p.cy, 6) for p in strands}) == 1
    placed = list(lay.placed.values()) + lay.ghosts
    assert not [
        (a.sse.id, b.sse.id)
        for i, a in enumerate(placed)
        for b in placed[i + 1 :]
        if min(a.rect[2], b.rect[2]) - max(a.rect[0], b.rect[0]) > 1e-6
        and min(a.rect[3], b.rect[3]) - max(a.rect[1], b.rect[1]) > 1e-6
    ]  # the extra helices, placed by projection, overlap nothing
    assert loop_crossings(route_loops(lay, sses, bb)) <= 4  # 29 at b43695e


def test_rossmann_sheet_reads_321456_with_helices_between():
    lay, sses, bb = make_layout(FID / "1LDM.cif.gz", assembly="asu")
    block = next(b for b in lay.sheets if len(b) == 6)
    xs = sorted(block, key=lambda k: lay.placed[k].cx)
    rank = {k: i for i, k in enumerate(sorted(block, key=lambda k: lay.placed[k].sse.start))}
    assert [rank[k] + 1 for k in xs] in ([3, 2, 1, 4, 5, 6], [6, 5, 4, 1, 2, 3])
    assert loop_crossings(route_loops(lay, sses, bb)) <= 10  # 21 at b43695e


@pytest.mark.parametrize("name", ["1UBQ", "5NKT", "1EMA", "2LZM"])
def test_template_only_for_all_parallel_sheets(name):
    on = make_layout(DATA / f"{name}.cif")[0]
    off = make_layout(DATA / f"{name}.cif", look=Style(fold_template="off"))[0]
    assert {k: (p.cx, p.cy, p.angle) for k, p in on.placed.items()} == {
        k: (p.cx, p.cy, p.angle) for k, p in off.placed.items()
    }


def test_template_off_restores_projection():
    lay = make_layout(DATA / "1TIM.cif", chains=["A"], look=Style(fold_template="off"))[0]
    strands = [lay.placed[k] for k in next(b for b in lay.sheets if len(b) == 8)]
    helices = [p for p in lay.placed.values() if p.sse.kind == "H"]
    assert not all(abs(h.cy - strands[0].cy) < 1e-6 for h in helices)  # helices around the barrel, as before
    assert lay.ba_rows == []


@pytest.mark.parametrize(
    "path, applies",
    [
        (DATA / "1TIM.cif", [True, True]),  # both protomers' barrels
        (FID / "7TIM.cif.gz", [True, True]),
        (FID / "1LDM.cif.gz", [True]),  # the NAD-binding sheet, 4 of 5 units in the row
        (DATA / "1UBQ.cif", [False]),  # mixed
        (DATA / "5NKT.cif", [False, False]),  # antiparallel and mixed
        (FID / "2ACY.cif.gz", [False]),
    ],
)
def test_template_applies_is_a_pure_test_of_the_sheet(path, applies):
    from foldmap.layout import template_applies

    _, sses, sheets = pipeline(path)
    big = [sh for sh in sheets if len(sh.strands) >= 4]
    before = [([s.id for s in sh.strands], list(sh.directions)) for sh in big]
    assert [template_applies(sh, sses) for sh in big] == applies
    assert [([s.id for s in sh.strands], list(sh.directions)) for sh in big] == before  # no side effects


def test_build_layout_takes_the_template_choice_from_the_style_by_default():
    """build_layout(template=None) follows style.fold_template, so a style passed on its own is honoured; an explicit
    template overrides it."""
    from foldmap.frame import view_frame
    from foldmap.layout import build_layout

    _, sses, sheets = pipeline(DATA / "1TIM.cif")
    sses = [s for s in sses if s.chain == "A"]
    sheets = [sh for sh in sheets if sh.strands[0].chain == "A"]
    f = view_frame(sses)

    def rows(**kw):
        return build_layout(sses, sheets, f, "projected", attract=False, **kw).ba_rows

    assert rows() and rows(style=Style()) and rows(style=Style(fold_template="off"), template="auto")
    assert rows(style=Style(fold_template="off")) == [] == rows(template="off")
    with pytest.raises(ValueError, match="fold_template"):
        rows(template="on")


@pytest.mark.parametrize(
    "path, kw",
    [(DATA / "1TIM.cif", {"chains": ["A"]}), (DATA / "1TIM.cif", {}), (FID / "7TIM.cif.gz", {"chains": ["A"]})],
)
def test_template_barrel_opens_between_its_last_and_first_strand(path, kw):
    """The textbook cut (the Commons TIM schematic): the barrel opens between β8 and β1, so the row reads β1 α1 β2 ...
    β8 (or back) and the ghost repeats an end strand. build_sheets cut 1TIM between β7 and β8: the β7 -> β8 loop then
    ran over the whole row from the real β7 to the closing helix beside the ghost."""
    lay, sses, bb = make_layout(path, **kw)
    loops = route_loops(lay, sses, bb)
    for row in lay.ba_rows:
        strands = sorted((lay.placed[k] for k in row if lay.placed[k].sse.kind == "E"), key=lambda p: p.cx)
        rank = [sorted(strands, key=lambda p: p.sse.start).index(p) + 1 for p in strands]
        assert rank in (list(range(1, 9)), list(range(8, 0, -1))), rank
        (ghost,) = [g for g in lay.ghosts if g.sse.chain == strands[0].sse.chain]
        assert ghost.sse in (strands[0].sse, strands[-1].sse)
        x0, x1 = strands[0].cx, max(strands[-1].cx, ghost.cx)
        for lp in loops:
            if lp.a_id in row or lp.b_id in row:
                xs = [x for x, _ in lp.points]
                assert max(xs) - min(xs) < 0.5 * (x1 - x0), (lp.a_id, lp.b_id)  # no loop over the whole row


def test_a_template_row_across_the_mirror_line_takes_side_by_side_copies():
    """1TIM seen with its dyad in the page (--up 0.927,0.051,-0.371): each protomer's row is wider than its half of
    the page and straddles the mirror line, so mirroring lifted B's row 9.4 units over A's while B's free helices
    stayed mirrored (Task 12 review). Such a figure drops the mirror: the protomers sit side by side, drawn alike."""
    from foldmap.route import clip_len

    lay, sses, bb = make_layout(DATA / "1TIM.cif", up=np.array([0.927, 0.051, -0.371]))
    assert lay.frame.mirror_x is not None and len(lay.ba_rows) == 2
    placed = list(lay.placed.values()) + lay.ghosts
    assert not [
        (a.sse.id, b.sse.id)
        for i, a in enumerate(placed)
        for b in placed[i + 1 :]
        if min(a.rect[2], b.rect[2]) - max(a.rect[0], b.rect[0]) > 1e-6
        and min(a.rect[3], b.rect[3]) - max(a.rect[1], b.rect[1]) > 1e-6
    ]
    loops = route_loops(lay, sses, bb)
    for row in lay.ba_rows:
        r = np.array([lay.placed[k].rect for k in row])
        box = (r[:, 0].min(), r[:, 1].min(), r[:, 2].max(), r[:, 3].max())
        chain = row[0].split(":")[0]
        crossing = [
            (lp.a_id, lp.b_id)
            for lp in loops
            if not lp.a_id.startswith(chain + ":")
            and any(clip_len(p, q, box) > 1e-9 for p, q in zip(lp.points, lp.points[1:]))
        ]
        assert not crossing, (chain, crossing)
    a, b = ([lay.placed[k] for k in row if lay.placed[k].sse.kind == "E"] for row in lay.ba_rows)
    assert {round(p.cy, 6) for p in a} == {round(p.cy, 6) for p in b}  # side by side, at one height

    # drawn alike: the same strands in the same order, and each ghost repeats the same strand (the rows were built
    # from opposite faces, so they opened at opposite ends)
    def ranks(row):
        return [sorted(row, key=lambda p: p.sse.start).index(p) for p in sorted(row, key=lambda p: p.cx)]

    assert ranks(a) == ranks(b)
    ga, gb = sorted(lay.ghosts, key=lambda g: g.sse.chain)
    assert sorted(a, key=lambda p: p.sse.start).index(lay.placed[ga.sse.id]) == sorted(
        b, key=lambda p: p.sse.start
    ).index(lay.placed[gb.sse.id])


def test_copies_of_one_protein_share_the_template_decision():
    """Copies of one protein (model.homolog_chains) are drawn alike: when one copy's sheet takes the template, the
    matching sheets of the other copies take it too (here B lacks three of its row helices, 4 of 7 units: alone it
    would not qualify); when a matching sheet cannot be a row at all (ambiguous, mixed), none of them is one."""
    import dataclasses

    from foldmap.frame import view_frame
    from foldmap.layout import build_layout, template_choices
    from foldmap.model import homolog_chains

    bb, sses, sheets = pipeline(DATA / "1TIM.cif")
    copies = homolog_chains(bb)
    gone = {"B:18-30", "B:45-54", "B:80-86"}
    sses = [s for s in sses if s.id not in gone]
    big = [sh for sh in sheets if len(sh.strands) == 8]
    from foldmap.layout import template_applies

    assert [template_applies(sh, sses) for sh in big] == [True, False]
    assert template_choices(big, sses, copies) == [True, True]
    assert template_choices(big, sses, {}) == [True, False]  # without copies, each sheet alone
    twisted = [big[0], dataclasses.replace(big[1], ambiguous=True)]
    assert template_choices(twisted, sses, copies) == [False, False]
    lay = build_layout(sses, sheets, view_frame(sses), "projected", attract=False, copies=copies)
    assert len(lay.ba_rows) == 2


@pytest.mark.parametrize("kw", [{"chains": ["A"]}, {}, {"up": np.array([0.927, 0.051, -0.371])}])
def test_ba_rows_list_each_row_left_to_right_as_drawn(kw):
    """Layout.ba_rows reads each row as drawn, left to right, also for a copy that took its mate's offsets."""
    lay = make_layout(DATA / "1TIM.cif", **kw)[0]
    assert lay.ba_rows
    for row in lay.ba_rows:
        assert row == sorted(row, key=lambda k: lay.placed[k].cx)


def test_open_at_ends_leaves_a_barrel_of_two_chains_as_it_is():
    from foldmap.layout import _open_at_ends
    from foldmap.model import Sheet
    from helpers import fake_sse

    strands = [fake_sse("E", 10 * k, 10 * k + 4, chain="AB"[k % 2]) for k in range(6)]
    sheet = Sheet(strands, [1] * 6, ["P"] * 5, closed=True)
    assert _open_at_ends(sheet) is sheet
    one = Sheet([fake_sse("E", s, s + 4) for s in (30, 20, 10, 0, 50, 40)], [1] * 6, ["P"] * 5, closed=True)
    opened = _open_at_ends(one)  # first (0) and last (50) in sequence are neighbours in the ring: opened there
    assert [s.start for s in opened.strands] == [50, 40, 30, 20, 10, 0] and opened.closed
