"""FigureSpec: one description of a figure for the command line, the app and the MCP server."""

import re
from pathlib import Path

import pytest

from foldmap import layoutfile
from foldmap.spec import FigureSpec, draw_spec, figure_spec, layout_document
from foldmap.style import THEMES, Style

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


# ------------------------------------------------------------------ references checked before drawing (Q-A19)
T4L = DATA / "2LZM.cif"  # T4 lysozyme: chain A, residues 1-164
HB = DATA / "2HHB.cif"  # haemoglobin as deposited: chains A, B, C, D


def _domain(*refs, name="N", **kw):
    from foldmap.spec import check_references

    return check_references(FigureSpec(T4L, domains=[(name, list(refs))]), **kw)


def test_each_domain_form_before_the_structure_is_read():
    from foldmap.spec import Problem

    assert _domain("A:1-60") == [Problem("--domain N=A:1-60", "read A:1-60 as the residue range res:A:1-60", False)]
    assert _domain("1-60") == [Problem("--domain N=1-60", "read 1-60 as the residue range res:1-60", False)]
    assert _domain("res:A:1-60") == [] and _domain("α1", "A:G", "#3", "res:A:-7-4", "res:H:100A-100K") == []
    (bad,) = _domain("A:x-y")
    assert bad.fatal and str(bad).startswith("--domain N=A:x-y: A:x-y is neither a residue range nor an element")
    assert "res:A:1-60" in bad.text and "α1" in bad.text  # the accepted forms
    (bad,) = _domain("res:A:x-y")
    assert bad.fatal and bad.text.startswith("res:A:x-y is no residue range; write res:150-159")
    for empty in ("", " "):
        (bad,) = _domain(empty)
        assert bad.fatal and bad.text.startswith("an empty reference")
    (bad,) = _domain()
    assert bad.fatal and bad.where == "--domain N=" and "needs elements" in bad.text
    assert [p.fatal for p in _domain("#0", "A:", ":α1")] == [True, True, True]
    assert _domain("A:1-60", elements={"A:1-60"}) == []  # an element id of the figure drawn: no warning


def test_rename_swap_move_and_mark_references_are_checked_alike():
    from foldmap.spec import check_references

    spec = FigureSpec(
        T4L,
        look=Style(mark="A:1-60=red, α1"),
        rename={"A:1-60": "N", "res:A:x": "C"},
        swap=[("A:61-70", "α2")],
        move=[("70-80", (1.0, 0.0))],
    )
    said = [(p.where, p.text, p.fatal) for p in check_references(spec)]
    assert said == [
        ("--rename A:1-60=N", "read A:1-60 as the residue range res:A:1-60", False),
        ("--rename res:A:x=C", "res:A:x is no residue range; write res:150-159, res:A:150-159, res:A:-7-4 or "
         "res:H:100A-100K", True),
        ("--swap A:61-70,α2", "read A:61-70 as the residue range res:A:61-70", False),
        ("--move 70-80", "read 70-80 as the residue range res:70-80", False),
        ("mark A:1-60=red", "read A:1-60 as the residue range res:A:1-60", False),
    ]  # fmt: skip
    # outside a domain a hyphen may be part of a name (--names, --rename): left to the layout to find
    assert check_references(FigureSpec(T4L, rename={"A:β-hairpin": "x"})) == []


def test_a_layout_files_domains_are_its_element_ids_and_not_checked():
    from foldmap.spec import SavedDomains, check_references

    assert check_references(FigureSpec(T4L, domains=SavedDomains([("N", ["A:999-1005"])]))) == []
    assert check_references(FigureSpec(T4L, domains="auto")) == []


def test_chains_and_residues_once_the_backbone_is_read():
    from foldmap.io import load_backbone
    from foldmap.spec import check_references

    bb = load_backbone(T4L)
    (bad,) = _domain("res:B:1-20", bb=bb)
    assert bad.fatal and bad.text == "no protein chain B in the structure; its protein chains are A"
    (bad,) = _domain("res:A:900-950", bb=bb)
    said = r"res:A:900-950 meets no residue of chain A: its residues are numbered 1 to 16\d"
    assert bad.fatal and re.fullmatch(said, bad.text)
    # read as a range, then checked as one: the warning, then the error
    assert [p.fatal for p in _domain("A:900-950", bb=bb)] == [False, True]
    (bad,) = _domain("900-950", bb=bb)[1:]
    assert bad.text.startswith("900-950 meets no residue of any chain drawn (numbered A 1 to 16")
    assert [p for p in _domain("res:A:1-60", "A:1-60", "res:A:150-200", "α1", bb=bb) if p.fatal] == []
    hb = load_backbone(HB, "asu")
    spec = FigureSpec(HB, domains=[("N", ["res:B:1-20"])], focus="A")
    (bad,) = check_references(spec, hb, drawn={"A"})
    assert bad.text == "chain B is not drawn: --focus A draws chain(s) A (add --focus none to draw every chain)"
    spec = FigureSpec(HB, rename={"B:α1": "x"}, chains=["A", "C"])
    (bad,) = check_references(spec, hb)
    assert bad.where == "--rename B:α1=x" and bad.text == "chain B is not drawn: --chains draws A, C only"


def test_make_layout_checks_the_references_before_laying_out():
    """Pass 2 runs in make_layout once the backbone is read and focus has chosen its chains: an error, never a deep
    ValueError or an empty domain (Review Focus 3)."""
    from foldmap.cli import BadReference, make_layout

    with pytest.raises(BadReference, match=r"^--domain N=res:A:900-950: res:A:900-950 meets no residue of chain A"):
        make_layout(T4L, domains=[("N", ["res:A:900-950"])])
    with pytest.raises(BadReference, match=r"^--domain N=res:B:1-20: chain B is not drawn: --focus A draws"):
        make_layout(HB, assembly="asu", focus="A", domains=[("N", ["res:B:1-20"])])
    with pytest.raises(BadReference, match=r"^--swap α1,res:Z:1-9: no protein chain Z .* are A, B, C, D$"):
        make_layout(HB, assembly="asu", swap=[("α1", "res:Z:1-9")])


def test_the_command_line_says_each_form(tmp_path, capsys):
    from foldmap.cli import main

    out = tmp_path / "t4l.svg"
    assert main(["plot", str(T4L), "-o", str(out), "--domain", "N=A:1-60"]) == 0
    err = capsys.readouterr().err
    assert "foldmap: warning: --domain N=A:1-60: read A:1-60 as the residue range res:A:1-60" in err
    assert "domain:N" in out.read_text()
    assert main(["plot", str(T4L), "-o", str(out), "--domain", "N=1-60"]) == 0
    assert "read 1-60 as the residue range res:1-60" in capsys.readouterr().err
    assert main(["plot", str(T4L), "-o", str(out), "--domain", "N=res:A:1-60"]) == 0
    assert "warning" not in capsys.readouterr().err
    for item, said in [
        ("N=A:x-y", "--domain N=A:x-y: A:x-y is neither a residue range nor an element"),
        ("N=res:B:1-20", "--domain N=res:B:1-20: no protein chain B in the structure; its protein chains are A"),
        ("N=res:A:900-950", "--domain N=res:A:900-950: res:A:900-950 meets no residue of chain A"),
        ("N=,", "--domain N=: a domain needs elements"),
        ("N=", "--domain needs NAME=REF[,REF...], got 'N='"),
    ]:
        assert main(["plot", str(T4L), "-o", str(tmp_path / "x.svg"), "--domain", item]) == 2, item
        assert capsys.readouterr().err.startswith(f"foldmap: error: {said}"), item
    assert not (tmp_path / "x.svg").exists()


def test_a_figure_drawn_like_a_reference_keeps_its_marks_to_itself():
    """Off means identical: the reference is laid out with this figure's look, whose marks name this figure's
    elements (chain B here); they are not checked against the reference (chain A), which never draws them."""
    from foldmap.cli import make_layout

    look = Style(color_by="bfactor", mark="res:B:1-20")  # not an identity colouring: no reference colours resolved
    like = {"structure": str(HB), "chain": "A"}
    lay, _, _ = make_layout(HB, assembly="asu", chains=["B"], look=look, like=like)
    assert lay.like_note
