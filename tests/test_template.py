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
    """The (βα)8 row of the textbook TIM figures: 3 loop crossings, against 30 for the projected figure (b43695e). The
    row is wider than the projected figure (27.2 x 13.0 = 353 units², against 20.3 x 14.2 = 288): eight strands, the
    ghost and seven helices, each helix with its name beside it, side by side. The crossings, not the area, decide
    whether the figure reads (Task 12 ruling: the area limit was a target)."""
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
    # six units join strands side by side; the seventh closes the barrel (β7 -> β8: first and last slot), its helix
    # in the row between the last strand and the ghost of the first
    assert inrow == 6
    (ghost,) = [g for g in lay.ghosts if g.sse.chain == "A"]
    last = by_x[-1]
    a, b = sorted((ghost.sse, last.sse), key=lambda s: s.start)
    closure = [p for p in lay.placed.values() if p.sse.kind == "H" and a.end < p.sse.start < b.start]
    assert any(
        last.cx < h.cx < ghost.cx and abs(h.cy - last.cy) < 1e-6 and np.cos(h.angle - last.angle) < -0.99
        for h in closure
    )
    assert [k for row in lay.ba_rows for k in row if lay.placed[k].sse.kind == "H"]  # Task 13 numbers these
    loops = route_loops(lay, sses, bb)
    assert loop_crossings(loops) <= 4  # 30 at b43695e
    assert lay.width * lay.height <= 353.0  # measured (x1.22 the projected figure's 288 units²)


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
