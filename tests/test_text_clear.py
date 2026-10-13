"""Small texts placed clear of lines, bodies and other text (random-100 batch A: C1 tail counts, C3 loop counts, C5
ligand marks, C6 a loop's own helix label, C28 glycan stems and linked-residue labels). The figures are small
fixtures cut from campaign entries (tests/data/r100); figcheck measures them as the campaign does."""

from pathlib import Path

import pytest

from foldmap import figcheck
from foldmap.layout import label_boxes, overlaps_any, tail_counts
from foldmap.render import ligand_marks, loop_counts
from foldmap.route import clip_len
from foldmap.spec import draw_spec, figure_spec

R100 = Path(__file__).parent / "data" / "r100"
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
