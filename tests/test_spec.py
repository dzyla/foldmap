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


def test_none_and_empty_arguments_keep_the_files_values(tmp_path):
    lf = layoutfile.save(
        tmp_path / "f.yaml",
        {
            "foldmap": 1,
            "layout": {
                "focus": "none",
                "chains": "A",
                "membrane": "flip",
                "msa": "a.fa",
                "uniprot": "P0",
                "rotate": 45,
            },
        },
    )
    spec = figure_spec(
        UBQ, layout_file=lf, focus=None, chains="", membrane="", msa=None, uniprot="", rotate=None, protomers=""
    )
    assert (spec.focus, spec.chains, spec.membrane, spec.msa, spec.uniprot, spec.rotate) == (
        "none", ["A"], "flip", "a.fa", "P0", 45.0,
    )  # fmt: skip
    assert spec.protomers is None


def test_fill_and_shorthand_over_file_style(tmp_path):
    lf = layoutfile.save(tmp_path / "f.yaml", {"foldmap": 1, "style": {"fill": "bold", "palette": "greys"}})
    spec = figure_spec(UBQ, layout_file=lf, fill="outline", style=[])
    assert spec.look.fill == "outline" and spec.look.palette == "greys"
    assert figure_spec(UBQ, layout_file=lf, palette="tol-muted").look.palette == "tol-muted"


def test_flip_v_from_file_and_explicit_false(tmp_path):
    lf = layoutfile.save(tmp_path / "f.yaml", {"foldmap": 1, "layout": {"flip_v": True, "membrane": False}})
    assert figure_spec(UBQ, layout_file=lf).flip_v is True
    assert figure_spec(UBQ, layout_file=lf, flip_v=False).flip_v is False
    assert figure_spec(UBQ, layout_file=lf, flip_v=None).flip_v is True
    assert figure_spec(UBQ, layout_file=lf).membrane == "off"


def test_edit_shapes():
    spec = figure_spec(UBQ, swap=["A,B", ("C", "D")], move={"res:A:2-7": [1, 2]})
    assert spec.swap == [("A", "B"), ("C", "D")] and spec.move == [("res:A:2-7", (1.0, 2.0))]
    with pytest.raises(ValueError, match="two elements"):
        figure_spec(UBQ, swap=["A,B,C"])
    with pytest.raises(ValueError, match="two numbers"):
        figure_spec(UBQ, move=[("res:A:2-7", [1])])


def test_round_trip_through_a_layout_file(tmp_path, monkeypatch):
    from foldmap import uniprot

    ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    stub = uniprot.Entry("P0CG48", "UBC_HUMAN", "Polyubiquitin-C", "UBC", "Homo sapiens", ubiquitin)
    monkeypatch.setattr(uniprot, "entry", lambda acc, cache_dir=None: stub)  # no network
    spec = figure_spec(
        UBQ,
        theme="blueprint",
        style=["residue_numbers=true"],
        title="T",
        mode="stack",
        rotate=30.0,
        flip_v=True,
        dna=False,
        symmetry="off",
        symmetry_tol=2.5,
        protomers=[["A"]],
        assembly="asu",
        up=[0.0, 0.0, 1.0],
        view=[1.0, 0.0, 0.0],
        domains=[("core", ["res:A:2-7"])],
        msa=str(DATA / "ubq_family.fasta"),
        msa_reference="1UBQ_A",
        membrane="flip",
        uniprot="P0CG48",
        chains=["A"],
        focus="none",
        rename={"res:A:2-7": "N"},
        swap=[("res:A:2-7", "res:A:12-16")],
        move=[("res:A:2-7", (1.0, -2.0))],
    )
    doc = layout_document(spec, structure_name="1UBQ.cif")
    path = layoutfile.save(tmp_path / "rt.yaml", doc)
    assert figure_spec(UBQ, layout_file=path) == spec


def test_an_empty_title_clears_the_layout_files_title(tmp_path):
    lf = layoutfile.save(tmp_path / "f.yaml", {"foldmap": 1, "layout": {"title": "From file"}})
    assert figure_spec(UBQ, layout_file=lf).title == "From file"
    assert figure_spec(UBQ, layout_file=lf, title=None).title == "From file"
    assert figure_spec(UBQ, layout_file=lf, title="").title is None
    assert figure_spec(UBQ, layout_file=lf, title="Mine").title == "Mine"


def test_dna_text_from_a_layout_file_is_read_as_yes_or_no(tmp_path):
    lf = layoutfile.save(tmp_path / "f.yaml", {"foldmap": 1, "layout": {"dna": "no"}})
    assert figure_spec(UBQ, layout_file=lf).dna is False


def test_flip_v_accepts_booleans_and_yes_no_text_only(tmp_path):
    def spec(value):
        lf = layoutfile.save(tmp_path / "f.yaml", {"foldmap": 1, "layout": {"flip_v": value}})
        return figure_spec(UBQ, layout_file=lf)

    assert spec(True).flip_v is True and spec("Yes").flip_v is True and spec("true").flip_v is True
    assert spec("no").flip_v is False and spec("FALSE").flip_v is False
    with pytest.raises(ValueError, match="flip_v"):
        spec("maybe")
    with pytest.raises(ValueError, match="flip_v"):
        spec(1)


def test_layout_document_of_a_conservation_figure(tmp_path):
    """A conservation-coloured figure with an alignment can be saved as a layout file and reloaded."""
    msa = str(DATA / "ubq_family.fasta")
    spec = figure_spec(UBQ, theme="conservation", msa=msa)
    doc = layout_document(spec)
    assert doc["theme"] == "conservation" and doc["layout"]["msa"] == msa
    again = figure_spec(UBQ, layout_file=layoutfile.save(tmp_path / "fig.yaml", doc))
    assert again.theme == "conservation" and str(again.msa) == msa
