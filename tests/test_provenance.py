"""Plan 12 Task 5 (review R7): one list of how a figure was made — what is inferred, ambiguous, or moved for
readability (provenance.figure_provenance) — said alike on the command line, in the app, by the MCP server and in
the figure files' own metadata."""

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from foldmap.cli import main
from foldmap.spec import draw_spec, figure_spec

DATA = Path(__file__).parent / "data"
R100 = DATA / "r100"
FID = DATA / "fidelity"
FIMD = R100 / "9FWZ-DB.cif.gz"
PENDANT = "1 strand(s) off a sheet's run of paired strands: D:662-666 drawn beside its partner D:645-659"
T4L_AUTO = "domains D1, D2 inferred from residue contacts, not a published domain assignment"


def _kinds(items) -> list[str]:
    return [i["kind"] for i in items]


def _one(items, kind) -> dict:
    (item,) = [i for i in items if i["kind"] == kind]
    return item


def test_every_item_has_its_kind_text_count_and_elements():
    items = draw_spec(figure_spec(FIMD)).provenance
    assert items and all(set(i) == {"kind", "text", "count", "elements"} for i in items)
    assert all(isinstance(i["count"], int) and i["count"] > 0 and isinstance(i["elements"], list) for i in items)


def test_a_figure_drawn_by_the_standard_method_alone_has_an_empty_list():
    assert draw_spec(figure_spec(DATA / "1UBQ.cif")).provenance == []


def test_a_pendant_strand_is_an_ambiguous_pairing_and_the_barrel_closure_a_ghost():
    """9FWZ (FimD): its pendant strand D:662-666 is drawn beside its partner, off the 24-strand ring (plan 12 T2), and
    the ring's first strand is drawn again, dashed, to close it."""
    items = draw_spec(figure_spec(FIMD)).provenance
    pairing = _one(items, "ambiguous_pairing")
    assert pairing["text"] == PENDANT and pairing["count"] == 1
    assert pairing["elements"] == ["D:662-666", "D:645-659"]
    ghost = _one(items, "ghost_strands")
    assert ghost["elements"] == ["D:214-223"] and "barrel" in ghost["text"] and "D:214-223" in ghost["text"]


def test_a_pairing_named_but_not_level_is_counted():
    """Plan 12 ruling: the count of accessories named but not drawn level with their partner (2OMF: 3) is shown."""
    items = draw_spec(figure_spec(R100 / "2OMF-A.cif.gz")).provenance
    (pairing,) = [i for i in items if i["kind"] == "ambiguous_pairing" and "off a sheet's run" in i["text"]]
    assert pairing["count"] == 3 and "3 not level with their partner, named only" in pairing["text"]


def test_gfp_closes_its_barrel_with_a_ghost_strand():
    items = draw_spec(figure_spec(DATA / "1EMA.cif")).provenance
    ghost = _one(items, "ghost_strands")
    assert ghost["elements"] == ["A:160-170"] and ghost["count"] == 1


def test_automatic_domains_are_said_to_be_inferred():
    items = draw_spec(figure_spec(DATA / "2LZM.cif", domains="auto")).provenance
    dom = _one(items, "domains")
    assert dom["text"] == T4L_AUTO and dom["elements"] == ["D1", "D2"] and dom["count"] == 2
    assert items[0]["kind"] == "domains"  # what is inferred first


def test_domains_given_by_hand_are_said_to_be_given():
    items = draw_spec(figure_spec(DATA / "2LZM.cif", domains=[("N", ["res:A:13-59"])])).provenance
    assert _one(items, "domains")["text"] == "domains N as given, not inferred"


def test_a_strand_joined_across_a_bulge_is_listed_with_the_bulge():
    """1TEN: DSSP splits A:880-890 at a beta-bulge; foldmap draws it as one strand (ss._join_bulged_pieces)."""
    items = draw_spec(figure_spec(FID / "1TEN.cif.gz")).provenance
    bulge = _one(items, "bulge_joined")
    assert bulge["elements"] == ["A:880-890"] and bulge["count"] == 1
    assert re.search(r"A:880-890 \(bulge at A:\d+(-\d+)?\)", bulge["text"]), bulge["text"]


def test_helix_pieces_matched_across_copies_are_listed():
    """7TIM: one helix DSSP breaks in one copy only is cut in the other, so the copies are drawn alike."""
    items = draw_spec(figure_spec(FID / "7TIM.cif.gz")).provenance
    matched = _one(items, "helix_matched")
    assert "1 cut" in matched["text"] and matched["elements"] and matched["count"] == len(matched["elements"])


def test_a_page_past_its_largest_is_listed(monkeypatch):
    from test_page_text import past_the_largest_page

    past_the_largest_page(monkeypatch)
    items = draw_spec(figure_spec(DATA / "1UBQ.cif")).provenance
    texts = [i["text"] for i in items if i["kind"] == "page_scaled"]
    assert len(texts) == 2 and texts[0].startswith("the page (3 x 2 in) is held at 3 in")


def test_names_moved_for_readability_are_listed():
    """2LZM with its gallery lobes: the N-lobe's name stands above its panel with a leader (test_domains)."""
    from test_domains import T4L_GALLERY

    d = draw_spec(figure_spec(DATA / "2LZM.cif", assembly="asu", domains=T4L_GALLERY))
    assert d.layout.domain_leaders
    moved = _one(d.provenance, "names_moved")
    assert set(d.layout.domain_leaders) <= set(moved["elements"]) and "leader" in moved["text"]


def test_hidden_additives_are_listed():
    items = draw_spec(figure_spec(DATA / "4CH9.cif")).provenance
    add = _one(items, "additives_hidden")
    assert add["elements"] == ["CL", "EDO", "SO4"] and "ligands=all" in add["text"]


def _desc(svg_text: str) -> str:
    root = ET.fromstring(svg_text)
    found = root.find("{http://www.w3.org/2000/svg}desc")
    return "" if found is None else found.text


def test_the_svg_desc_holds_the_list(tmp_path):
    out = tmp_path / "f.svg"
    d = draw_spec(figure_spec(FIMD))
    from foldmap.render import save, save_svg

    save(d.figure, out)
    desc = _desc(out.read_text(encoding="utf-8"))
    assert desc.startswith("How this figure was made:")
    assert all(f"- {i['text']}" in desc for i in d.provenance)
    assert _desc(save_svg(d.figure)) == desc  # the app's inline SVG says it too
    plain = tmp_path / "u.svg"
    save(draw_spec(figure_spec(DATA / "1UBQ.cif")).figure, plain)
    assert "<desc" not in plain.read_text(encoding="utf-8")  # nothing to say: the file as before


def test_pdf_and_png_carry_the_list_in_their_metadata(tmp_path):
    from foldmap.render import save

    d = draw_spec(figure_spec(FIMD))
    pdf, png = save(d.figure, tmp_path / "f.pdf"), save(d.figure, tmp_path / "f.png")
    assert b"/Subject" in pdf.read_bytes()
    from PIL import Image

    assert Image.open(png).text["Description"].startswith("How this figure was made:")


@pytest.mark.parametrize("ext", ["svg", "pdf", "png"])
def test_figure_files_stay_deterministic(tmp_path, ext):
    from foldmap.render import save

    a = save(draw_spec(figure_spec(FIMD)).figure, tmp_path / f"a.{ext}").read_bytes()
    b = save(draw_spec(figure_spec(FIMD)).figure, tmp_path / f"b.{ext}").read_bytes()
    assert a == b


def test_the_command_line_prints_one_block(tmp_path, capsys):
    assert main(["plot", str(FIMD), "-o", str(tmp_path / "f.svg")]) == 0
    err = capsys.readouterr().err
    assert err.count("How this figure was made:") == 1
    assert f"\n  - {PENDANT}\n" in err
    assert "  - warning: 1 loop(s) drawn over a label or chain end, the only route found: D:618-629>D:634-641" in err


def test_the_command_line_says_automatic_domains_are_inferred(tmp_path, capsys):
    assert main(["plot", str(DATA / "2LZM.cif"), "-o", str(tmp_path / "t.svg"), "--domains", "auto"]) == 0
    assert f"How this figure was made:\n  - {T4L_AUTO}\n" in capsys.readouterr().err


def test_the_mcp_result_carries_the_list(tmp_path):
    pytest.importorskip("mcp")
    from test_mcp import call

    r = call("draw_topology", {"source": str(FIMD), "output_path": str(tmp_path / "f.svg")})
    assert [c.type for c in r.content] == ["text", "image"]
    items = r.structured_content["provenance"]
    assert items == draw_spec(figure_spec(FIMD)).provenance
    assert f"; {PENDANT}" in r.content[0].text  # the note, as before
    r = call("draw_topology", {"source": str(DATA / "1UBQ.cif"), "output_path": str(tmp_path / "u.svg")})
    assert r.structured_content["provenance"] == []
