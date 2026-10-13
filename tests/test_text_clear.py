"""Small texts placed clear of lines, bodies and other text (random-100 batch A: C1 tail counts, C3 loop counts, C5
ligand marks, C6 a loop's own helix label, C28 glycan stems and linked-residue labels). The figures are small
fixtures cut from campaign entries (tests/data/r100); figcheck measures them as the campaign does."""

from pathlib import Path

import pytest

from foldmap import figcheck
from foldmap.layout import label_boxes, overlaps_any, tail_counts
from foldmap.spec import draw_spec, figure_spec

R100 = Path(__file__).parent / "data" / "r100"


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
