"""Copies drawn alike (plan 15 Task 5, random-100 held-out A15): of several copy sheets tied to one mate only one is
tied, the others are items of their own; what the tied pushes leave overlapping slides clear an orbit at a time; copies
on a duplex frame are tied side by side; copies drawn with no tie are reported as such; a letter placed again off one
loop sends the loop it lands on round it."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from foldmap import layout
from foldmap.layout import Placed, _collide, _Item
from helpers import fake_sse, r100_entry


def _item(x, y, *members, chain="A"):
    """An item at (x, y) of upright strands (chain, first, last, x offset)."""
    placed = [Placed(fake_sse("E", a, b, chain=c), 0.0, 0.0, 2.0, np.pi / 2, 0.5) for c, a, b, _ in members]
    return _Item(placed, [(ox, 0.0) for *_, ox in members], x, y, (x, y), members[0][1], chain)


def _offsets(items, ties):
    return {i: (round(items[i].cx - items[j].cx, 9), round(items[i].cy - items[j].cy, 9)) for i, (j, _) in ties.items()}


def _clear(items):
    h, m = [it.half for it in items], [it.mid for it in items]
    return not any(
        _collide(items[i], h[i], m[i], items[j], h[j], m[j])
        for i in range(len(items))
        for j in range(i + 1, len(items))
    )


# --------------------------------------------------------------------------------------- one copy item per mate


def test_of_two_copy_items_tied_to_one_mate_the_larger_keeps_its_place():
    """2BBJ E~A 1.42 (9GTP's 1K, S, 2e-2i, 3e-3i alike): DSSP split E's sheet into two blocks where A's is one, and
    _symmetry_ties tied both to A's sheet as protomer 4. Starting on one place, the two overlapped, and the pushes,
    which move an orbit as one, could only part them by halves: each stood off its mate's place by half the
    overlap, and E's other items with it were drawn unlike A's. Only the block with the most elements is tied (in
    build_layout, straight after _symmetry_ties); the other is an item of its own, which the copies' pushes leave
    where they hold."""

    def made():
        mate = _item(0.0, 0.0, ("A", 1, 6, 0.0), ("A", 10, 15, 0.6), ("A", 20, 25, 1.2))
        helix = _item(0.0, 4.0, ("A", 40, 50, 0.0))
        big = _item(20.0, 0.0, ("B", 1, 6, 0.0), ("B", 20, 25, 1.2))
        small = _item(20.0, 0.0, ("B", 10, 13, 0.0))  # the split block: on its mate's place, over the big one
        copy_helix = _item(20.0, 4.0, ("B", 40, 50, 0.0))
        return [mate, helix, big, small, copy_helix]

    ties = {2: (0, 1), 3: (0, 1), 4: (1, 1)}
    items = made()
    assert layout._one_per_mate(items, ties) == {2: (0, 1), 4: (1, 1)}
    layout._separate(items, ties, again=True)  # both tied: parted by halves, the sheet off its mate's place
    got = _offsets(items, {2: (0, 1), 4: (1, 1)})
    assert got[2] != got[4], got
    items = made()
    state = layout._separate(items, layout._one_per_mate(items, ties), again=True)
    assert state in ("alike", "wide")
    got = _offsets(items, {2: (0, 1), 4: (1, 1)})
    assert got[2] == got[4], got  # the copy's sheet and helix stand as its mate's do
    assert _clear(items)
    src = Path(layout.__file__).read_text()  # and build_layout ties them so
    assert (
        "_one_per_mate(items, ties, helices=not regroup)"
        in src[src.index("ties = _symmetry_ties(items, frame") :][:900]
    )
    # a helix cut in two: judged too, unless _regroup_copies runs next (it joins the pieces end to end)
    piece = [_item(0.0, 0.0, ("A", 1, 20, 0.0))]
    h1 = Placed(fake_sse("H", 1, 9, chain="B"), 0.0, 0.0, 2.0, 0.0, 0.5)
    h2 = Placed(fake_sse("H", 12, 20, chain="B"), 0.0, 0.0, 2.0, 0.0, 0.5)
    pieces = piece + [_Item([h], [(0.0, 0.0)], 20.0, 0.0, (20.0, 0.0), 1, "B") for h in (h1, h2)]
    assert layout._one_per_mate(pieces, {1: (0, 1), 2: (0, 1)}) == {1: (0, 1)}
    assert layout._one_per_mate(pieces, {1: (0, 1), 2: (0, 1)}, helices=False) == {1: (0, 1), 2: (0, 1)}


# --------------------------------------------------------------------------------------- orbits slide as one


@pytest.mark.parametrize("mirror", [None, 10.0])
def test_an_orbit_slides_clear_as_one(mirror):
    """13FZ C~A 17.1: the tied pushes traded a jam back and forth between items and the walls of a trace, the same
    in both copies, and gave up: every item was moved alone, and the copies were drawn unlike. What the pushes leave
    overlapping slides clear with each orbit (an item and its copies) moved by one step, a mirrored copy the other
    way across: they stay alike."""
    walls = [_item(0.0, 0.0, ("Z", 1, 6, 0.0)), _item(20.0, 0.0, ("Z", 11, 16, 0.0))]
    a = _item(0.3, 0.0, ("A", 20, 25, 0.0))
    b = _item(20.3 if mirror is None else 19.7, 0.0, ("B", 20, 25, 0.0))
    items = [*walls, a, b]
    ties = {3: (2, 1)}
    before = (b.cx - a.cx, b.cx + a.cx)
    layout._slide_clear(items, ties=ties, mirror=mirror)
    if mirror is None:
        assert b.cx - a.cx == pytest.approx(before[0]) and b.cy == a.cy
    else:  # mirror images about x = 10 stay mirror images
        assert b.cx + a.cx == pytest.approx(before[1]) and b.cy == a.cy
    assert _clear(items)


def test_the_second_chance_keeps_copies_alike_when_the_pushes_leave_a_jam(monkeypatch):
    """The second chance (_tied_again) used to give up when its pushes left an overlap between tied items: it slides
    them clear an orbit at a time and keeps the result when nothing overlaps and the room allows."""
    mate, other = _item(0.0, 0.0, ("A", 1, 6, 0.0)), _item(0.3, 0.0, ("A", 20, 25, 0.0))
    copy, copy_other = _item(20.0, 0.0, ("B", 1, 6, 0.0)), _item(20.3, 0.0, ("B", 20, 25, 0.0))
    items = [mate, other, copy, copy_other]
    ties = {2: (0, 1), 3: (1, 1)}
    monkeypatch.setattr(layout, "_push_apart", lambda *a, **k: False)  # a jam the pushes never clear
    start = np.array([(it.cx, it.cy) for it in items], float)
    halves = np.array([it.half for it in items], float)
    mids = np.array([it.mid for it in items], float)
    apart = layout._joint_matrix(len(items), set())
    assert layout._tied_again(items, ties, start, halves, mids, start.copy(), set(), apart, None, None, None)
    got = _offsets(items, ties)
    assert got[2] == got[3], got
    assert _clear(items)


# --------------------------------------------------------------------------------------- duplex frame


def test_copies_on_a_duplex_frame_are_tied_side_by_side_unless_they_share_the_duplex(monkeypatch):
    """5K7Z: two dimers, the duplex held by A and B alone; the straight duplex chose the duplex frame, whose dyad was
    not in the page, so the frame knew no protomers and the four copies were drawn unalike with no tie. The frame now
    learns them (no mirror line) unless a strand is held by two protomers (unrolled, they would pull it apart)."""
    from foldmap import cli

    sym = SimpleNamespace(protomers=[["A", "B"], ["C", "D"]])
    for shared, want in ((False, {"A": 0, "B": 0, "C": 1, "D": 1}), (True, None)):
        monkeypatch.setattr(cli, "_na_shared", lambda *a, _s=shared: _s)
        frame = SimpleNamespace(sector=None, classes=None, mirror_x=None)
        cli._side_by_side(frame, sym, None, None)
        assert frame.sector == want and frame.mirror_x is None


# --------------------------------------------------------------------------------------- A15 reports untied copies


def test_copies_drawn_with_no_tie_are_not_judged_rather_than_absent():
    """5K7Z: four copies drawn with no tie said "does not apply: no copies to compare". Copies there are: A15 says
    they were not judged and why; a figure with no copies still does not apply."""
    from foldmap.criteria import NOT_APPLICABLE, NOT_JUDGED, _copies
    from foldmap.figcheck import copies_untied
    from foldmap.model import ResLabel

    seq = ["ALA", "GLY", "SER", "LEU", "VAL", "ILE", "PRO", "THR", "ASP", "GLU"] * 3
    labels = [ResLabel(c, k + 1, "", n) for c in "ABC" for k, n in enumerate(seq if c != "C" else seq[::-1])]
    placed = {f"{c}:1-5": Placed(fake_sse("H", 1, 5, chain=c), 0.0, 0.0, 2.0, 0.0, 0.5) for c in "ABC"}
    drawn = SimpleNamespace(
        layout=SimpleNamespace(copy_ties={}, placed=placed, partial=set()), bb=SimpleNamespace(labels=labels)
    )
    assert copies_untied(drawn) == ["B~A"]
    drawn.layout.partial = {"B"}  # a focus view's lent fragment is no copy drawn
    assert copies_untied(drawn) == []

    def view(untied):
        return {"checks": {"copies_alike": 0.0, "details": {"copies_alike": {}, "copies_untied": untied}}}

    got = _copies({"default": view(["B~A"]), "focus": view([])})
    assert got["pass"] is None and got["status"] == NOT_JUDGED and "B~A" in got["detail"]
    got = _copies({"default": view([]), "focus": view([])})
    assert got["status"] == NOT_APPLICABLE and got["detail"] == "no copies to compare"


# --------------------------------------------------------------------------------------- a letter sends a loop round


def test_a_letter_placed_again_onto_another_loop_sends_that_loop_round_it(monkeypatch):
    """13FZ B:6: placed again off one loop (_letters_off_loops), its only spot not on the sheet under it lay on another
    loop, which had been routed before the letter moved there. When a letter moves, the loops are routed off the
    nucleic acid's text once more (_off_na_text), then the letters still on a loop placed again."""
    from pathlib import Path

    from foldmap import na_labels, route
    from foldmap.cli import make_layout

    calls = []
    real_text, real_letters = route._off_na_text, na_labels._letters_off_loops

    def text(*a):
        calls.append("text")
        return real_text(*a)

    def letters(*a):
        calls.append("letters")
        got = real_letters(*a)
        return got if calls.count("letters") > 1 else {(0, 0)}  # the first pass moved one

    monkeypatch.setattr(route, "_off_na_text", text)
    monkeypatch.setattr(na_labels, "_letters_off_loops", letters)
    lay, sses, bb = make_layout(Path(__file__).parent / "data" / "na" / "1EHZ.cif.gz")
    route.route_loops(lay, sses, bb)
    assert calls == ["text", "letters", "text", "letters"]


# --------------------------------------------------------------------------------------- the entries


@pytest.mark.parametrize("entry", ["2BBJ", "13FZ", "5K7Z"])
def test_heldout_copies_are_drawn_alike(entry):
    """The held-out entries: 2BBJ E~A (split sheet), 13FZ C~A (a jam against its RNA traces; and B:6's letter off
    the loops), 5K7Z (two dimers on a duplex frame, tied)."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(r100_entry(entry)))
    got = figcheck.copies_alike(d)
    assert got and max(got.values()) <= 0.5, got
    assert [h for h in figcheck.text_collisions(d) if h.startswith("na:nt:")] == []
