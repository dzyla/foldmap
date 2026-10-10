"""FigureSpec: one description of a figure for the command line, the app and the MCP server."""

from pathlib import Path

import pytest

from foldmap import layoutfile
from foldmap.spec import FigureSpec, draw_spec, figure_spec, layout_document
from foldmap.style import THEMES

DATA = Path(__file__).parent / "data"
UBQ = DATA / "1UBQ.cif"


def test_defaults():
    spec = figure_spec(UBQ)
    assert spec == FigureSpec(UBQ, look=THEMES["publication"])
    assert spec.mode == "projected" and spec.rotate == 0.0 and spec.edits() == {"rename": {}, "swap": [], "move": []}


def test_layout_file_then_given_values_including_falsy_ones(tmp_path):
    lf = layoutfile.save(
        tmp_path / "fig.yaml",
        {
            "foldmap": 1,
            "theme": "blueprint",
            "style": {"residue_numbers": True},
            "layout": {"rotate": 90, "focus": "auto", "chains": "A", "up": "0,0,1", "membrane": "auto"},
            "domains": {"core": ["res:A:1-40"]},
            "edits": {"rename": {"res:A:1-7": "N"}, "move": {"res:A:10-17": [1, -2]}},
        },
    )
    spec = figure_spec(
        UBQ, layout_file=lf, rotate=0.0, focus="none", rename={"res:A:23-34": "M"}, style=["legend=false"]
    )
    assert spec.theme == "blueprint" and spec.look.background == THEMES["blueprint"].background
    assert spec.look.residue_numbers is True and spec.look.legend is False
    assert spec.rotate == 0.0 and spec.focus == "none" and spec.membrane == "auto"
    assert spec.chains == ["A"] and spec.up == [0.0, 0.0, 1.0] and spec.domains == [("core", ["res:A:1-40"])]
    assert spec.rename == {"res:A:1-7": "N", "res:A:23-34": "M"} and spec.move == [("res:A:10-17", (1.0, -2.0))]


@pytest.mark.parametrize(
    "kw, message",
    [
        ({"theme": "nope"}, "unknown theme 'nope'"),
        ({"mode": "sideways"}, "unknown mode 'sideways'"),
        ({"colour": "red"}, "unknown figure option"),
        ({"up": "1,2"}, "three numbers"),
    ],
)
def test_bad_values_are_explained(kw, message):
    with pytest.raises(ValueError, match=message):
        figure_spec(UBQ, **kw)


def test_draw_and_document_use_the_same_spec():
    spec = figure_spec(UBQ, theme="rainbow", rotate=30.0, swap=[("A", "B")])
    drawn = draw_spec(spec)
    assert drawn.figure is not None and drawn.stuck == [] and len(drawn.sses) == len(drawn.layout.placed)
    doc = layout_document(spec, structure_name="1UBQ.cif")
    assert doc["structure"] == "1UBQ.cif" and doc["theme"] == "rainbow" and doc["layout"]["rotate"] == 30.0
    assert doc["edits"]["swap"] == [["res:A:2-7", "res:A:12-16"]]


def test_shorthands_membrane_and_passthrough_options():
    spec = figure_spec(
        UBQ, palette="okabe-ito", loops="curved", style=["loops=orthogonal"], membrane="flip", msa="x.grades"
    )
    assert spec.look.palette == "okabe-ito" and spec.look.loops == "orthogonal"
    assert spec.membrane == "flip" and spec.msa == "x.grades"
    with pytest.raises(ValueError, match="unknown membrane 'up'"):
        figure_spec(UBQ, membrane="up")
