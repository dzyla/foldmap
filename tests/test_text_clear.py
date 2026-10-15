"""Small texts placed clear of lines, bodies and other text (random-100 batch A: C1 tail counts, C3 loop counts, C5
ligand marks, C6 a loop's own helix label, C28 glycan stems and linked-residue labels). The figures are small
fixtures cut from campaign entries (tests/data/r100); figcheck measures them as the campaign does."""

from pathlib import Path

import numpy as np
import pytest

from foldmap import figcheck
from foldmap.layout import label_boxes, overlaps_any, tail_counts
from foldmap.render import ligand_marks, linked_marks, loop_counts
from foldmap.route import clip_len
from foldmap.spec import draw_spec, figure_spec

DATA = Path(__file__).parent / "data"
R100 = DATA / "r100"
FID = Path(__file__).parent / "data" / "fidelity"


def drawn(path, **kw):
    return draw_spec(figure_spec(path, **kw))


def hits(d, kind: str) -> list[str]:
    return [h for h in figcheck.text_collisions(d) if kind in h]


@pytest.mark.parametrize("scale", [1.0, 2.0])
def test_tail_count_sits_clear_beyond_its_letter(scale):
    """C1 (random-100): a long tail's "NN aa" was written on the stub's line beyond the N/C letter with no obstacle
    test, so loops ran through it (2DGS: the loop A:156-161>A:164-173 through the N tail's "11 aa"), and on big
    pages it overlapped its own letter (8W9Y, 8QEM). It is placed after routing (a keep-out moved the
    nucleosome's loops onto its DNA) from the letter's far edge, sized from the drawn text, clear of elements,
    text and loops; the page holds it."""
    d = drawn(R100 / "2DGS-A.cif.gz", style=[f"font_scale={scale}"])
    lay = d.layout
    counts = tail_counts(lay, d.sses, d.loops)
    assert [(t["end"], t["text"]) for t in counts] == [("N", "11 aa"), ("C", "14 aa")]
    boxes = dict(label_boxes(lay, d.sses))
    texts = {t.get_gid(): t for t in d.figure.axes[0].texts}
    for t in counts:
        end = t["end"]
        assert not overlaps_any(t["box"], list(boxes.values()))
        assert not overlaps_any(t["box"], [p.rect for p in lay.placed.values()])
        assert texts[f"tail-length:{end}:A"].get_position() == pytest.approx(t["centre"])
    assert hits(d, "tail-length") == []
    xmin, xmax = d.figure.axes[0].get_xlim()
    ymin, ymax = d.figure.axes[0].get_ylim()
    assert all(xmin < t["box"][0] and t["box"][2] < xmax and ymin < t["box"][1] and t["box"][3] < ymax for t in counts)


def test_tail_count_looks_further_when_its_near_spots_are_crossed():
    """7BBH, 6JO6, 5LDW (random-100 A6 tail-length x loop, batch C2a): on a crowded page every one of the five spots
    beside the N/C letter was crossed by a loop, and the count went on the first of them, over the loop. Here a loop
    is laid through each spot the count takes, one after another; past the five near spots it looks further out
    beside the letter and along its stub, and takes a spot no loop crosses."""
    from foldmap.route import Loop

    d = drawn(R100 / "2DGS-A.cif.gz")
    lay, loops = d.layout, list(d.loops)
    for _ in range(5):
        (x0, y0, x1, y1) = next(t["box"] for t in tail_counts(lay, d.sses, loops) if t["end"] == "N")
        cy = (y0 + y1) / 2
        loops.append(Loop("x", "y", [(x0 - 0.3, cy), (x1 + 0.3, cy)], False, False))
    box = next(t["box"] for t in tail_counts(lay, d.sses, loops) if t["end"] == "N")
    lines = [(p, q) for lp in loops[len(d.loops) :] for p, q in zip(lp.points, lp.points[1:])]
    assert not any(clip_len(p, q, box) > 1e-9 for p, q in lines)


def test_loop_count_sits_clear_of_elements_lines_and_text():
    """C3 (random-100): a long loop's "NN aa" was written beside its longest straight run with no test, so it sat on
    strands (9ON5: the count of A:121-126>A:154-161 on strand 154-161), other loops and labels, and lay outside the
    page (7BRG focus). It now takes the first of a few spots (the longest runs, either side, near the middle) clear
    of elements, loops, text and earlier counts, and the page holds it."""
    d = drawn(R100 / "9ON5-A.cif.gz")
    counts = loop_counts(d.layout, d.loops, d.sses)
    assert [c["name"] for c in counts] == ["A:121-126>A:154-161"]
    assert hits(d, "loop-length") == []
    ax = d.figure.axes[0]
    texts = {t.get_gid(): t for t in ax.texts}
    for c in counts:
        assert texts[f"loop-length:{c['name']}"].get_position() == pytest.approx(c["centre"])
        x0, y0, x1, y1 = c["box"]
        assert ax.get_xlim()[0] < x0 < x1 < ax.get_xlim()[1] and ax.get_ylim()[0] < y0 < y1 < ax.get_ylim()[1]
    assert figcheck.off_page(d) == []


def _mark_box(m) -> tuple[float, float, float, float]:
    (cx, cy), (hw, hh) = m["centre"], m["half"]
    return cx - hw, cy - hh, cx + hw, cy + hh


@pytest.mark.parametrize("path, kw", [(R100 / "1BXI.cif.gz", {}), (FID / "1MBN.cif.gz", {})], ids=["1BXI", "1MBN"])
def test_ligand_marks_sit_off_loops_and_earlier_tethers(path, kw):
    """C5 (random-100): a ligand's mark was kept clear of elements, labels and other marks but not of the loops'
    lines, so 1BXI's Ni sat on the loop B:22-25>B:32-33 (and 1A0Z's hems, 3MHT's SAH, 5SB5's Ca); nor were earlier
    marks' tethers and names obstacles to later ones (6JO6: chlorophyll names crossed by other chlorophylls'
    tethers). Marks now keep off loop lines, counts, earlier tethers and names wherever a spot is free."""
    spec = figure_spec(path, **kw)
    d = draw_spec(spec)
    marks = ligand_marks(d.layout, d.loops, d.sses, spec.look)
    assert marks
    for i, m in enumerate(marks):
        box = _mark_box(m)
        assert not any(clip_len(p, q, box) > 1e-9 for lp in d.loops for p, q in zip(lp.points, lp.points[1:])), m
        for e in marks[:i]:  # no earlier tether runs through it, and its tethers through no earlier mark
            assert not any(clip_len(e["centre"], t, box) > 1e-9 for t in e["targets"])
            assert not any(clip_len(m["centre"], t, _mark_box(e)) > 1e-9 for t in m["targets"])
    assert hits(d, "ligand") == []


@pytest.mark.parametrize("scale", [1.0, 0.5])
def test_linked_label_place_is_sized_from_its_ink(scale):
    """C28 (random-100): a linked residue's label (GFP's CRO66) is drawn at no less than 7 pt, but its place was
    tested with a fixed box narrower than its ink (8K8E: "JUU15" at 7 pt on a 3.6 pt page, 2.3 page units wide
    against a 0.9 box, ran over the η1 label). The box is now its ink's, so the spot it takes is clear."""
    d = drawn(DATA / "1EMA.cif", style=[f"font_scale={scale}"])
    marks = [m for m in linked_marks(d.layout, d.loops, d.sses, []) if m["text_at"]]
    assert marks
    to_page = d.figure.axes[0].transData.inverted()
    ink = {t.get_gid(): to_page.transform([[b.x0, b.y0], [b.x1, b.y1]]) for t, b in figcheck._ink_boxes(d.figure)}
    for m in marks:
        (x0, y0), (x1, y1) = ink[f"linked-label:{m['label'].chain}:{m['label'].tag}"]
        (cx, cy), (hw, hh) = m["text_at"], m["half"]
        assert cx - hw <= x0 and x1 <= cx + hw and cy - hh <= y0 and y1 <= cy + hh, (m["text"], x0, x1, cx, hw)
    assert hits(d, "linked") == []


def _stem_and_first_sugar(fig, r):
    patches = {p.get_gid(): p for p in fig.axes[0].patches}
    at, _ = patches[f"glycan:{r}:stem"].get_path().vertices
    sugar = patches[f"glycan:{r}:0"].get_window_extent()
    centre = fig.axes[0].transData.inverted().transform([(sugar.x0 + sugar.x1) / 2, (sugar.y0 + sugar.y1) / 2])
    return at, centre


def test_glycan_stem_crosses_no_name(monkeypatch):
    """C28 (random-100): a glycan's sugars were kept off names, but not its stem (1WWL domains: the stem of the
    glycan at 254 ran across strand K's letter). A name put on a 6ZS5 stem, between the chain and the first sugar,
    now turns the stem another way."""
    import foldmap.render as render

    d = drawn(DATA / "6ZS5.cif")
    r = sorted(g.split(":")[1] for g in (p.get_gid() or "" for p in d.figure.axes[0].patches) if g.endswith(":stem"))[0]
    at, sugar = _stem_and_first_sugar(d.figure, r)
    mid = (at + sugar) / 2
    name = (mid[0] - 0.05, mid[1] - 0.05, mid[0] + 0.05, mid[1] + 0.05)
    real = render.label_boxes
    monkeypatch.setattr(render, "label_boxes", lambda *a, **k: [*real(*a, **k), ("label:test", name)])
    at, sugar = _stem_and_first_sugar(drawn(DATA / "6ZS5.cif").figure, r)
    assert clip_len(at, sugar, name) == 0


def test_a_ligand_is_never_left_without_a_mark():
    """C5 fix round 1: with every ring spot crossed by earlier tethers or loop lines, a chosen ligand still gets a
    mark where the old code put it (6JO6 lost 24 ligands on a crowded page when the tether test had no fallback)."""
    from foldmap.render import _place_mark

    anchor = np.array([0.0, 0.0])
    half = (0.3, 0.3)
    wall = [(-9.0, y, 9.0, y + 0.1) for y in np.arange(-4.0, 4.0, 0.2)]  # loop lines over every ring spot
    cross = [((-9.0, y), (9.0, y)) for y in np.arange(-4.0, 4.0, 0.2)]  # and tethers of earlier marks
    spot = _place_mark(anchor, half, [(1.0, 0.0)], [], wall, [], cross)
    assert spot is not None


def test_a_ligand_boxed_in_by_elements_and_text_still_gets_a_mark():
    """Batch B's campaign: on 7QUP's 65-chain page GTP:9E500 found no spot clear of elements and text within
    _SPOT_RADIUS and was dropped (90 ligand names -> 89). The search widens ring by ring past that radius, and a mark
    with nowhere clear at all is put at its anchor: never dropped."""
    from foldmap.render import _SPOT_RADIUS, _place_mark

    anchor, half = np.array([0.0, 0.0]), (0.3, 0.3)
    r = _SPOT_RADIUS + 0.5
    near = [(-r, -r, r, r)]  # an element box over every spot within the search radius
    spot = _place_mark(anchor, half, [(1.0, 0.0)], near, [], [], [])
    assert spot is not None and max(abs(spot[0]), abs(spot[1])) > r
    everywhere = [(-1e3, -1e3, 1e3, 1e3)]
    assert _place_mark(anchor, half, [(1.0, 0.0)], everywhere, [], [], []) is not None


def test_a_mark_with_nowhere_clear_takes_the_least_bad_spot_and_says_so():
    """Fix round 1 (minor): with nothing clear even far out, the mark went on its anchor, over whatever was there,
    silently. It takes the spot overlapping the fewest element bodies, then the fewest other boxes (here: off the
    element, on a label box), and the figure records it (Drawn.crowded_marks; the CLI warns)."""
    from foldmap.render import _SPOT_FAR, _SPOT_RADIUS, _place_mark

    anchor, half = np.array([0.0, 0.0]), (0.3, 0.3)
    far = _SPOT_FAR * _SPOT_RADIUS + 1.0
    body = (-far, -far, 0.5, far)  # an element over the left and the middle
    label = (0.5, -far, far, far)  # a label over the right
    spot = _place_mark(anchor, half, [(1.0, 0.0)], [body, label], [], [], [], [body])
    assert spot is not None and spot[0] - half[0] > 0.5  # on the label, not the element
    d = draw_spec(figure_spec(R100 / "1BXI.cif.gz"))
    assert d.crowded_marks == []


def test_every_chosen_ligand_of_a_crowded_page_has_a_mark():
    """Same, end to end: 1GOX's FMN and 1BXI's metal, with a loop line laid over every spot near each."""
    spec = figure_spec(R100 / "1BXI.cif.gz")
    d = draw_spec(spec)
    from foldmap.route import Loop

    crowd = [Loop("x", "y", [(x, y), (x + 40.0, y)], False, False) for x in (-20.0,) for y in np.arange(-20, 20, 0.2)]
    full = ligand_marks(d.layout, d.loops, d.sses, spec.look)
    jammed = ligand_marks(d.layout, [*d.loops, *crowd], d.sses, spec.look)
    assert full and {m["ligand"].name for m in jammed} == {m["ligand"].name for m in full}


def test_line_clearance_is_figchecks_pixels_at_the_drawn_text_size():
    """Fix round 1 (minor): C14's letter spots and C11's reroutes kept 0.07 and 0.1 em from NA text, which is at
    least figcheck.LINE_TEXT_CLEARANCE_PX only near 7 pt text; on a wide page (4 pt text) it fell short. They now keep
    layout.line_clearance: that many pixels (and half a line) in page units at the page's text scale, as the helix-name
    mover does."""
    from foldmap.figcheck import LINE_TEXT_CLEARANCE_PX
    from foldmap.layout import Layout, Placed, line_clearance, page_text
    from helpers import fake_sse

    for width in (10.0, 160.0):
        lay = Layout({"A:1-10": Placed(fake_sse("H", 0, 9), 0.0, 0.0, 1.0, 0.0, 0.5)}, (0.0, 0.0, width, 10.0))
        lay.placed["B:1-10"] = Placed(fake_sse("H", 20, 29, chain="B"), width - 2.0, 0.0, 1.0, 0.0, 0.5)
        per_pt = page_text(lay)[1]
        assert line_clearance(lay) == pytest.approx((LINE_TEXT_CLEARANCE_PX + 0.5) * 0.72 * per_pt)
