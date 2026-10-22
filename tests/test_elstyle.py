"""Each element's own colour and width set by hand (plan 23 Task 1): --color / --width, the style keys
element_colors / element_widths, layout files, the app's Colour and Width editors and the MCP server."""

import asyncio
from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap import criteria, figcheck
from foldmap.cli import BadReference, build_parser, main
from foldmap.elstyle import colours, entries, joined, widths
from foldmap.spec import check_references, draw_spec, figure_spec, layout_document
from foldmap.style import LAYOUT_KEYS, THEME_KEYS, Style, resolve_style

DATA = Path(__file__).parent / "data"
UBQ, HHB, LMB = DATA / "1UBQ.cif", DATA / "2HHB.cif", DATA / "1LMB.cif"


def _gids(fig, prefix):
    ax = fig.axes[0]
    return [a for a in [*ax.patches, *ax.texts, *ax.lines] if (a.get_gid() or "").startswith(prefix)]


def _ids(drawn, label, chain="A"):
    return [k for k, p in drawn.layout.placed.items() if p.label == label and p.sse.chain == chain]


def test_the_keys_parse_and_say_what_is_wrong():
    assert Style().element_colors == Style().element_widths == ""
    assert "element_widths" in LAYOUT_KEYS and "element_colors" in THEME_KEYS
    assert colours("α1=red, loop:β2=#1F77B4,na:B=teal") == [
        ("α1", "#ff0000"),
        ("loop:β2", "#1f77b4"),
        ("na:B", "#008080"),
    ]
    assert widths("α1=2,res:A:1-10=0.5") == [("α1", 2.0), ("res:A:1-10", 0.5)]
    assert entries("a=b=c", "element_colors") == [("a=b", "c")]  # the value follows the last '='
    assert joined("α1=red", ["β1=blue"]) == "α1=red,β1=blue" and joined("", []) == ""
    for bad, words in (
        ("element_colors=α1", "REF=COLOUR"),
        ("element_colors=α1=nocolour", "not a colour"),
        ("element_widths=α1=5", "0.25 to 4"),
        ("element_widths=α1=wide", "0.25 to 4"),
    ):
        with pytest.raises(ValueError, match=words):
            resolve_style(overrides=[bad])


def test_the_command_line_writes_the_style_keys_after_the_file_s():
    args = build_parser().parse_args(
        ["plot", str(UBQ), "-o", "x.svg", "--set", "element_colors=α1=blue", "--color", "α1=red"]
        + ["--colour", "A=#00ff00", "--width", "α1=2", "--width", "loop:α1=1.5"]
    )
    from foldmap.cli import _spec_of

    look = _spec_of(args).look
    assert look.element_colors == "α1=blue,α1=red,A=#00ff00" and look.element_widths == "α1=2,loop:α1=1.5"


def test_references_are_checked_before_layout():
    spec = figure_spec(UBQ, colors=["Q:α1=red", "loop:=red", "na:=red"], widths=["res:A:900-910=2"])
    found = [str(p) for p in check_references(spec)]
    assert any(p.startswith("--color loop:") and "names nothing" in p for p in found)
    assert any(p.startswith("--color na:") for p in found)
    from foldmap.io import load_backbone

    bb = load_backbone(UBQ)
    late = [str(p) for p in check_references(spec, bb) if p.fatal]
    assert any("no protein chain Q" in p for p in late) and any("meets no residue" in p for p in late)


def test_a_reference_that_names_nothing_is_an_error_naming_it():
    with pytest.raises(ValueError, match=r"element_colors α9: no element"):
        draw_spec(figure_spec(UBQ, colors=["α9=red"]))
    with pytest.raises(ValueError, match=r"element_widths na:B: no nucleic-acid chain"):
        draw_spec(figure_spec(UBQ, widths=["na:B=2"]))
    with pytest.raises((ValueError, BadReference)):
        draw_spec(figure_spec(UBQ, colors=["loop:res:A:900-910=red"]))


def test_a_colour_set_by_hand_wins_over_the_colouring_and_the_legend_says_so():
    plain = draw_spec(figure_spec(UBQ, style=["color_by=sequence"]))
    assert not _gids(plain.figure, "legend-hand")
    drawn = draw_spec(figure_spec(UBQ, style=["color_by=sequence", "mark=α1=blue"], colors=["α1=#d1495b"]))
    (helix,) = _gids(drawn.figure, f"helix:{_ids(drawn, 'α1')[0]}")
    assert to_hex(helix.get_facecolor()) == "#d1495b"  # over color_by and the mark
    assert _gids(drawn.figure, "legend-hand")
    said = lambda d: "colours set by hand" in [t.get_text() for t in d.figure.axes[0].texts]  # noqa: E731
    assert said(drawn) and not said(plain)


def test_a_loop_s_colour_and_width():
    plain = draw_spec(figure_spec(UBQ))
    drawn = draw_spec(figure_spec(UBQ, colors=["loop:α1=#2ca02c"], widths=["loop:α1=2"]))
    a1 = _ids(drawn, "α1")[0]
    (line,) = [x for x in _gids(drawn.figure, f"loop:{a1}>") if x.get_gid().count(">") == 1]
    (before,) = [x for x in _gids(plain.figure, f"loop:{a1}>") if x.get_gid().count(">") == 1]
    assert to_hex(line.get_edgecolor()) == "#2ca02c"
    assert line.get_linewidth() == pytest.approx(2 * before.get_linewidth())
    assert drawn.layout.loop_widths == {a1: 2.0} and drawn.layout.loop_stroke == pytest.approx(2.0)
    others = [x for x in _gids(drawn.figure, "loop:") if not x.get_gid().startswith(f"loop:{a1}>")]
    assert all(x.get_linewidth() == pytest.approx(before.get_linewidth()) for x in others)


def _gap_after(label: str) -> str:
    """A residue range strictly inside the loop after element `label` of 1UBQ, as res:A:FIRST-LAST."""
    drawn = draw_spec(figure_spec(UBQ))
    placed = sorted((p.sse for p in drawn.layout.placed.values()), key=lambda s: s.start)
    k = next(i for i, s in enumerate(placed) if drawn.layout.placed[s.id].label == label)
    first, last = placed[k].last.seq + 1, placed[k + 1].first.seq - 1
    assert first <= last
    return f"res:A:{first}-{last}"


def test_a_residue_range_in_a_loop_names_that_loop():
    drawn = draw_spec(figure_spec(UBQ, colors=[f"{_gap_after('α1')}=#2ca02c"]))  # no element there
    a1 = _ids(drawn, "α1")[0]
    assert any(to_hex(x.get_edgecolor()) == "#2ca02c" for x in _gids(drawn.figure, f"loop:{a1}>"))


def test_a_wider_helix_and_strand_keep_their_room():
    plain = draw_spec(figure_spec(UBQ))
    drawn = draw_spec(figure_spec(UBQ, widths=["α1=2", "A=2"]))
    for label in ("α1", "A"):
        (k,) = _ids(drawn, label)
        assert drawn.layout.placed[k].width == pytest.approx(2 * plain.layout.placed[k].width)
    strands = sorted((p for p in drawn.layout.placed.values() if p.sse.kind == "E"), key=lambda p: p.cx)
    for a, b in zip(strands, strands[1:]):  # neighbours stay at least their half widths apart
        assert b.cx - a.cx >= (a.width + b.width) / 2 - 1e-6
    others = [k for k, p in plain.layout.placed.items() if p.label not in ("α1", "A")]
    assert all(drawn.layout.placed[k].width == plain.layout.placed[k].width for k in others)


@pytest.mark.parametrize("path, refs", [(UBQ, ["α1=2", "A=2", "loop:α1=2"]), (HHB, ["res:A:1-75=2"])])
def test_wide_elements_keep_text_and_boxes_off(path, refs):
    """A5 and A6 hold at width 2 (A6 measured only where the text is in Liberation Sans: test_criteria._font_off)."""
    spec = figure_spec(path, widths=refs)
    drawn = draw_spec(spec)
    got = criteria.check_figure(drawn, spec, path)["criteria"]
    assert got["A5"]["pass"] is True, got["A5"]
    assert got["A6"]["pass"] is (None if figcheck.ink_font_reason(drawn.figure) else True), got["A6"]


def test_a_nucleic_acid_backbone_s_colour_and_width():
    plain = draw_spec(figure_spec(LMB))
    drawn = draw_spec(figure_spec(LMB, colors=["na:1=#d1495b"], widths=["na:1=2"]))
    assert drawn.layout.na_widths == {"1": 2.0}
    lw = {
        to_hex(r.get_edgecolor()): r.get_linewidth() for r in _gids(drawn.figure, "na:stem:") if ":rails" in r.get_gid()
    }
    (base,) = {r.get_linewidth() for r in _gids(plain.figure, "na:stem:") if ":rails" in r.get_gid()}
    assert lw["#d1495b"] == pytest.approx(2 * base) and min(lw.values()) == pytest.approx(base)
    legend = [t.get_text() for t in drawn.figure.axes[0].texts]
    assert "DNA (1)" in legend and "colours set by hand" in legend


def test_a_ligand_s_colour_and_width():
    drawn = draw_spec(figure_spec(HHB, colors=["lig:HEM:A142=#2ca02c"], widths=["lig:HEM=2"]))
    plain = draw_spec(figure_spec(HHB))
    boxes = {x.get_gid(): x for x in _gids(drawn.figure, "ligand:HEM")}
    before = {x.get_gid(): x for x in _gids(plain.figure, "ligand:HEM")}
    assert to_hex(boxes["ligand:HEM:A142"].get_facecolor()) == "#2ca02c"
    assert all(to_hex(b.get_facecolor()) != "#2ca02c" for k, b in boxes.items() if k != "ligand:HEM:A142")
    assert all(b.get_linewidth() == pytest.approx(2 * before[k].get_linewidth()) for k, b in boxes.items())


def test_the_layout_file_keeps_them_as_residue_ranges_and_reads_them_back(tmp_path):
    from foldmap import layoutfile

    gap = _gap_after("α1")
    spec = figure_spec(UBQ, colors=["α1=#d1495b", "loop:α1=teal", f"{gap}=red"], widths=["α1=2", "loop:A=1.5"])
    drawn = draw_spec(spec)
    doc = layout_document(spec, drawn=drawn)
    style = doc["style"]
    assert style["element_colors"] == f"res:A:23-34=#d1495b,loop:res:A:23-34=teal,{gap}=red"
    assert style["element_widths"].startswith("res:A:23-34=2,loop:res:A:")
    path = layoutfile.save(tmp_path / "ubq.yaml", doc)
    again = draw_spec(figure_spec(UBQ, layout_file=path))
    assert {k: p.width for k, p in again.layout.placed.items()} == {k: p.width for k, p in drawn.layout.placed.items()}
    assert again.layout.loop_widths == drawn.layout.loop_widths
    colour = lambda d: {x.get_gid(): to_hex(x.get_facecolor()) for x in _gids(d.figure, "helix:")}  # noqa: E731
    assert colour(again) == colour(drawn)
    out = tmp_path / "again.yaml"
    assert main(
        ["plot", str(UBQ), "--layout-file", str(path), "-o", str(tmp_path / "a.svg"), "--save-layout", str(out)]
    ) in (0, None)
    assert layoutfile.load(out)["style"] == layoutfile.load(path)["style"]  # saved again: the same (A13)


def test_the_app_editors_write_what_the_command_line_takes():
    from foldmap.app_controls import hand_choices, reference_problems
    from foldmap.app_edits import rows_of, text_of

    rows = [{"element": "α1", "colour": "#d1495b"}, {"element": "loop:β2", "colour": "teal"}]
    text = text_of("element_colors", rows)
    cli = figure_spec(UBQ, colors=["α1=#d1495b", "loop:β2=teal"]).look.element_colors
    assert text == cli and rows_of("element_colors", text) == rows
    rows = [{"element": "α1", "width": 2.0}, {"element": "na:B", "width": 1.5}]
    text = text_of("element_widths", rows)
    assert text == "α1=2,na:B=1.5" == figure_spec(UBQ, widths=["α1=2", "na:B=1.5"]).look.element_widths
    assert rows_of("element_widths", text) == rows
    assert [str(p) for p in reference_problems({"element_widths": "α1=9"}, "element_widths")] == [
        "element_widths: element_widths: α1 has '9'; a width is a factor from 0.25 to 4"
    ]
    assert reference_problems({"element_colors": "α1=red"}, "element_colors") == []
    got = hand_choices([("res:A:1-7", "A · A 1–7"), ("res:A:23-34", "α1 · A 23–34")], ["B"], {"res:A:23-34"})
    assert ("loop:res:A:1-7", "loop after A") in got and ("na:B", "nucleic acid B") in got
    assert not any(ref == "loop:res:A:23-34" for ref, _ in got)


def test_the_mcp_server_takes_colours_and_widths(tmp_path):
    pytest.importorskip("mcp")
    from foldmap.mcp_server import server

    out = tmp_path / "ubq.svg"
    asyncio.run(
        server.call_tool(
            "draw_topology",
            {"source": str(UBQ), "output_path": str(out), "colors": {"α1": "#d1495b"}, "widths": {"α1": 2}},
        )
    )
    svg = out.read_text()
    assert "colours set by hand" in svg and "#d1495b" in svg


def test_without_them_the_layout_holds_no_widths():
    drawn = draw_spec(figure_spec(UBQ))
    lay = drawn.layout
    assert (lay.element_widths, lay.loop_widths, lay.na_widths, lay.ligand_widths) == ({}, {}, {}, [])
    assert np.isclose(lay.loop_stroke, 1.0) and not _gids(drawn.figure, "legend-hand")


def test_the_app_draws_them_and_blames_the_editor_for_a_wrong_reference():
    pytest.importorskip("streamlit")
    from app_helpers import as_text, opened, svg_of

    at = opened(UBQ)
    as_text(at, "element_colors", "α1=#d1495b")
    as_text(at, "element_widths", "α1=2,loop:α1=1.5")
    assert "figure_failed" not in at.session_state and "colours set by hand" in svg_of(at)
    assert "#d1495b" in svg_of(at)
    as_text(at, "element_widths", "α9=2")  # names nothing: the figure fails, the Width editor is named
    assert at.session_state["figure_failed"]
    assert any("Width" in e.value and "α9" in e.value for e in at.error), [e.value for e in at.error]
