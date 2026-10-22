"""The MCP server: tools agents use to summarise structures and draw figures."""

import asyncio
import re
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from foldmap.mcp_server import server  # noqa: E402

DATA = Path(__file__).parent / "data"
UBQ = str(DATA / "1UBQ.cif")


@pytest.fixture(autouse=True)
def _no_output_dir(monkeypatch):
    monkeypatch.delenv("FOLDMAP_OUTPUT_DIR", raising=False)


def call(name, args):
    return asyncio.run(server.call_tool(name, args))


def test_tools_are_listed_with_descriptions():
    tools = asyncio.run(server.list_tools())
    names = {t.name for t in tools}
    assert names == {
        "summarize_structure",
        "list_styles",
        "draw_topology",
        "compare_structures",
        "diff_structures",
        "draw_sequence",
        "interactive_page",
        "uniprot_info",
        "structure_data",
    }
    assert all(t.description for t in tools)


def test_summary_lists_labels_ranges_and_ligands():
    text = call("summarize_structure", {"source": str(DATA / "1ZAA.cif")}).content[0].text
    assert "α1: C 19-30 (helix" in text and "metal ZN" in text and "held by CYS" in text


def test_draw_topology_writes_the_file_and_returns_a_preview(tmp_path):
    out = tmp_path / "ubq.svg"
    r = call(
        "draw_topology",
        {
            "source": str(DATA / "1UBQ.cif"),
            "output_path": str(out),
            "style": {"color_by": "sequence", "residue_numbers": "true"},
            "swap": [["A", "B"]],
            "domains": {"core": ["res:A:1-40"]},
        },
    )
    kinds = [c.type for c in r.content]
    assert out.is_file() and "resnum" in out.read_text() and "domain:core" in out.read_text()
    assert kinds == ["text", "image"] and str(out) in r.content[0].text


def test_errors_come_back_as_messages(tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    with pytest.raises(ToolError, match="conservation"):
        call(
            "draw_topology",
            {"source": str(DATA / "1UBQ.cif"), "output_path": str(tmp_path / "x.svg"), "theme": "conservation"},
        )


def test_draw_sequence_with_alignment(tmp_path):
    out = tmp_path / "seq.png"
    r = call(
        "draw_sequence",
        {
            "source": str(DATA / "1UBQ.cif"),
            "output_path": str(out),
            "msa_path": str(DATA / "ubq_family.fasta"),
            "columns": 40,
        },
    )
    assert out.stat().st_size > 1000 and [c.type for c in r.content] == ["text", "image"]


def test_styles_and_interactive(tmp_path):
    assert "blueprint" in call("list_styles", {}).content[0].text
    out = tmp_path / "v.html"
    call("interactive_page", {"source": str(DATA / "1UBQ.cif"), "output_path": str(out)})
    assert "<html" in out.read_text()


def test_server_speaks_mcp_over_stdio(tmp_path):
    from mcp import Client, StdioServerParameters

    async def go():
        params = StdioServerParameters(command=sys.executable, args=["-m", "foldmap.mcp_server"])
        async with Client(params) as client:
            tools = await client.list_tools()
            result = await client.call_tool("summarize_structure", {"source": str(DATA / "1UBQ.cif")})
            return {t.name for t in tools.tools}, result.content[0].text

    names, text = asyncio.run(go())
    assert "draw_topology" in names and "76 residues" in text


def test_draw_topology_can_draw_the_whole_assembly(tmp_path):
    whole, one = tmp_path / "whole.svg", tmp_path / "one.svg"
    call("draw_topology", {"source": str(DATA / "6Y7S.cif"), "output_path": str(whole), "focus": "none"})
    call("draw_topology", {"source": str(DATA / "6Y7S.cif"), "output_path": str(one)})
    assert whole.read_text().count('id="strand:') > 2 * one.read_text().count('id="strand:')


def test_bad_theme_is_a_message_not_an_exit(tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    for tool in ("draw_topology", "draw_sequence"):
        with pytest.raises(ToolError, match="unknown theme 'nope'"):
            call(tool, {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "theme": "nope"})
    with pytest.raises(ToolError, match="unknown theme 'nope'"):
        call("interactive_page", {"source": UBQ, "output_path": str(tmp_path / "x.html"), "theme": "nope"})
    with pytest.raises(ToolError, match="'weird'"):
        call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "mode": "weird"})


def test_unreadable_files_are_a_message(tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    junk = tmp_path / "notes.cif"
    junk.write_text("this is not a structure\n")
    with pytest.raises(ToolError):
        call("summarize_structure", {"source": str(junk)})


def test_server_survives_a_bad_argument_over_stdio(tmp_path):
    from mcp import Client, StdioServerParameters

    async def go():
        params = StdioServerParameters(command=sys.executable, args=["-m", "foldmap.mcp_server"])
        async with Client(params) as client:
            bad = await client.call_tool(
                "draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "theme": "nope"}
            )
            good = await client.call_tool("summarize_structure", {"source": UBQ})
            return bad, good

    bad, good = asyncio.run(asyncio.wait_for(go(), 120))
    assert bad.is_error and "unknown theme" in bad.content[0].text
    assert "76 residues" in good.content[0].text


@pytest.mark.parametrize(
    "tool, args",
    [
        ("interactive_page", {"output_path": "page.sh"}),
        ("draw_topology", {"output_path": "fig.svg", "save_layout": "layout.sh"}),
        ("draw_topology", {"output_path": "fig.html"}),
        ("draw_sequence", {"output_path": "seq.txt"}),
    ],
)
def test_outputs_must_have_the_right_extension(tool, args, tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    args = {k: str(tmp_path / v) for k, v in args.items()}
    with pytest.raises(ToolError, match="must end in"):
        call(tool, {"source": UBQ, **args})
    assert not any(tmp_path.iterdir())


def test_output_dir_confines_writes(tmp_path, monkeypatch):
    from mcp.server.mcpserver.exceptions import ToolError

    monkeypatch.setenv("FOLDMAP_OUTPUT_DIR", str(tmp_path / "allowed"))
    with pytest.raises(ToolError, match="outside FOLDMAP_OUTPUT_DIR"):
        call("draw_sequence", {"source": UBQ, "output_path": str(tmp_path / "elsewhere.svg")})
    call("draw_sequence", {"source": UBQ, "output_path": str(tmp_path / "allowed" / "s.svg")})
    assert (tmp_path / "allowed" / "s.svg").is_file()


def test_output_dir_follows_relative_paths_and_symlinks(tmp_path, monkeypatch):
    from mcp.server.mcpserver.exceptions import ToolError

    allowed = tmp_path / "allowed"
    allowed.mkdir()
    (allowed / "escape").symlink_to(tmp_path)
    monkeypatch.setenv("FOLDMAP_OUTPUT_DIR", str(allowed))
    monkeypatch.chdir(allowed)
    call("draw_sequence", {"source": UBQ, "output_path": "rel.svg"})
    assert (allowed / "rel.svg").is_file()
    with pytest.raises(ToolError, match="outside FOLDMAP_OUTPUT_DIR"):
        call("draw_sequence", {"source": UBQ, "output_path": "escape/out.svg"})
    assert not (tmp_path / "out.svg").exists()


def test_output_dir_covers_layouts_pages_and_default_outputs(tmp_path, monkeypatch):
    from mcp.server.mcpserver.exceptions import ToolError

    root = tmp_path / "root"
    monkeypatch.setenv("FOLDMAP_OUTPUT_DIR", str(root))
    outside = tmp_path / "outside"
    for tool, args in [
        ("draw_topology", {"save_layout": str(outside / "l.yaml")}),
        ("interactive_page", {"output_path": str(outside / "p.html")}),
    ]:
        with pytest.raises(ToolError, match="outside FOLDMAP_OUTPUT_DIR"):
            call(tool, {"source": UBQ, **args})
    assert not outside.exists()
    assert not root.exists()  # nothing was created by the rejected calls
    call("draw_sequence", {"source": UBQ})
    assert len(list(root.glob("foldmap-*/*.svg"))) == 1
    call("draw_sequence", {"source": UBQ, "output_path": str(root / "UPPER.SVG")})
    assert (root / "UPPER.SVG").is_file()


def test_layout_file_settings_are_not_overridden_by_tool_defaults(tmp_path):
    from foldmap.cli import main

    lf = tmp_path / "fig.yaml"
    first = ["plot", UBQ, "-o", str(tmp_path / "a.svg"), "--theme", "blueprint", "--rotate", "90"]
    assert main([*first, "--save-layout", str(lf)]) == 0
    assert main(["plot", UBQ, "-o", str(tmp_path / "cli.svg"), "--layout-file", str(lf)]) == 0
    call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "mcp.svg"), "layout_file": str(lf)})
    assert (tmp_path / "mcp.svg").read_text() == (tmp_path / "cli.svg").read_text()


def test_draw_topology_lays_out_once_and_previews_that_figure(tmp_path, monkeypatch):
    import foldmap.cli as cli

    calls = []
    real = cli.make_layout
    monkeypatch.setattr(cli, "make_layout", lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "u.svg"), "rename": {"res:A:1-7": "Zq"}})
    assert len(calls) == 1 and [c.type for c in r.content] == ["text", "image"]
    assert ">Zq<" in (tmp_path / "u.svg").read_text()


def test_draw_topology_takes_chains_and_checks_moves(tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    out = tmp_path / "hb.svg"
    call("draw_topology", {"source": str(DATA / "2HHB.cif"), "output_path": str(out), "chains": ["A"], "focus": "none"})
    assert 'id="helix:A:' in out.read_text() and 'id="helix:B:' not in out.read_text()
    with pytest.raises(ToolError, match="two numbers"):
        call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "move": {"A": [1]}})


def test_a_structure_without_elements_still_gets_a_preview(tmp_path):
    from helpers import tripeptide

    tri = tripeptide(DATA / "1UBQ.cif", tmp_path)
    r = call("draw_topology", {"source": str(tri), "output_path": str(tmp_path / "t.svg")})
    assert [c.type for c in r.content] == ["text", "image"]


def _spec_drawn_by(monkeypatch, args):
    """The FigureSpec draw_topology builds for args (the drawing itself is stopped)."""
    from mcp.server.mcpserver.exceptions import ToolError

    import foldmap.spec as spec_mod

    seen = []

    def stop(spec):
        seen.append(spec)
        raise ValueError("stopped")

    monkeypatch.setattr(spec_mod, "draw_spec", stop)
    with pytest.raises(ToolError, match="stopped"):
        call("draw_topology", {"source": UBQ, **args})
    return seen[0]


def test_explicit_domains_beat_uniprot_domains(monkeypatch):
    spec = _spec_drawn_by(
        monkeypatch, {"uniprot": "auto", "uniprot_domains": True, "domains": {"core": ["res:A:1-40"]}}
    )
    assert spec.domains == [("core", ["res:A:1-40"])]


def test_uniprot_domains_work_with_a_layout_files_uniprot(tmp_path, monkeypatch):
    lf = tmp_path / "fig.yaml"
    lf.write_text("foldmap: 1\nlayout: {uniprot: P0CG48}\n", encoding="utf-8")
    spec = _spec_drawn_by(monkeypatch, {"layout_file": str(lf), "uniprot_domains": True})
    assert spec.uniprot == "P0CG48" and spec.domains == "uniprot"
    spec = _spec_drawn_by(monkeypatch, {"layout_file": str(lf)})
    assert spec.domains is None


def test_uniprot_domains_without_any_uniprot_is_an_error(tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    with pytest.raises(ToolError, match="needs UniProt annotation"):
        call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "uniprot_domains": True})


@pytest.mark.parametrize(
    "tool, args",
    [
        ("draw_topology", {"theme": "nope"}),
        ("draw_topology", {"mode": "weird"}),
        ("draw_topology", {"move": {"A": [1]}}),
        ("draw_sequence", {"theme": "nope"}),
        ("draw_sequence", {"msa_path": "/nowhere/x.fasta"}),
        ("interactive_page", {"theme": "nope"}),
        ("interactive_page", {"viewer": "nope"}),
    ],
)
def test_rejected_calls_leave_no_temporary_folders(tool, args, tmp_path, monkeypatch):
    from mcp.server.mcpserver.exceptions import ToolError

    root = tmp_path / "root"
    monkeypatch.setenv("FOLDMAP_OUTPUT_DIR", str(root))
    with pytest.raises(ToolError):
        call(tool, {"source": UBQ, **args})
    assert not root.exists() or not any(root.iterdir())


def test_error_messages_keep_their_quotes(tmp_path):
    from mcp.server.mcpserver.exceptions import ToolError

    with pytest.raises(ToolError) as err:
        call("draw_sequence", {"source": UBQ, "output_path": str(tmp_path / "fig.txt")})
    assert "'fig.txt' must end in" in str(err.value)  # the repr's opening quote is not stripped


def test_backbone_colour_keys_reach_the_figure_through_the_tool(tmp_path):
    out = tmp_path / "lmb.svg"
    call(
        "draw_topology",
        {"source": str(DATA / "1LMB.cif"), "output_path": str(out), "style": {"na_backbone": "#123456"}},
    )
    assert "#123456" in out.read_text()
    from mcp.server.mcpserver.exceptions import ToolError

    with pytest.raises(ToolError, match="backbone_colors"):
        call(
            "draw_topology",
            {"source": str(DATA / "1LMB.cif"), "output_path": str(out), "style": {"backbone_colors": "Z=red"}},
        )


def test_draw_topology_takes_names(tmp_path, monkeypatch):
    from foldmap import spec as spec_module

    seen = []
    real = spec_module.draw_spec
    monkeypatch.setattr(spec_module, "draw_spec", lambda s: seen.append(s.names) or real(s))
    out = tmp_path / "ubq.svg"
    r = call("draw_topology", {"source": UBQ, "output_path": str(out), "names": ["beta"]})
    assert out.is_file() and str(out) in r.content[0].text and seen == [["beta"]]


def test_draw_topology_says_when_it_drew_part_of_the_file(tmp_path):
    r = call("draw_topology", {"source": str(DATA / "6Y7S.cif"), "output_path": str(tmp_path / "one.svg")})
    text = r.content[0].text
    assert '; drew one of 6 subunits of a helical filament: 1 of 6 chains; focus="none" draws all of them' in text
    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "u.svg")})
    assert "drew " not in r.content[0].text.replace("wrote", "")


def test_a_crowded_ligand_mark_is_noted(tmp_path, monkeypatch):
    """Batch B fix round 2: the MCP note says when a ligand mark had no clear spot (as it says stuck and over-label
    loops). 1MBN's marks forced onto their first target (an element)."""
    import numpy as np

    from foldmap import render

    mbn = str(Path(__file__).parent / "data" / "fidelity" / "1MBN.cif.gz")
    r = call("draw_topology", {"source": mbn, "output_path": str(tmp_path / "m.svg")})
    assert "no clear spot" not in r.content[0].text
    monkeypatch.setattr(render, "_place_mark", lambda anchor, half, targets, *a, **k: np.asarray(targets[0], float))
    r = call("draw_topology", {"source": mbn, "output_path": str(tmp_path / "n.svg")})
    assert "ligand mark(s) with no clear spot" in r.content[0].text


def test_a_page_held_at_its_largest_is_noted(tmp_path, monkeypatch):
    """random-100 C4: the MCP note says when the page was held at layout.MAX_PAGE_IN (text under 6 pt)."""
    from test_page_text import past_the_largest_page

    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "u.svg")})
    assert "held at" not in r.content[0].text
    past_the_largest_page(monkeypatch)
    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "v.svg")})
    assert "; the page (3 x 2 in) is held at 3 in on its longer side" in r.content[0].text


def test_strands_drawn_beside_their_partner_are_noted(tmp_path):
    """Plan 12 Task 2: the MCP note says which strands are drawn beside their partner, off a sheet's run of paired
    strands (9FWZ's D:662-666), as the command line and the app say it."""
    fimd = str(Path(__file__).parent / "data" / "r100" / "9FWZ-DB.cif.gz")
    r = call("draw_topology", {"source": fimd, "output_path": str(tmp_path / "f.svg")})
    assert "; 1 strand(s) off a sheet's run of paired strands: D:662-666 drawn beside its partner D:645-659" in (
        r.content[0].text
    )
    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "u.svg")})
    assert "beside its partner" not in r.content[0].text


def test_the_preview_of_a_large_page_is_capped(tmp_path, monkeypatch):
    """Batch C1 review (finding 3): a page grown for 6 pt text (up to 48 in) made a 4320 px preview at 90 dpi; the
    preview's longer side is capped at PREVIEW_PX, and a small page keeps 90 dpi."""
    import base64

    from foldmap import mcp_server

    def width(r):
        image = next(c for c in r.content if getattr(c, "type", "") == "image")
        data = base64.b64decode(image.data)
        return int.from_bytes(data[16:20], "big"), int.from_bytes(data[20:24], "big")

    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "u.svg")})
    small = width(r)
    monkeypatch.setattr(mcp_server, "PREVIEW_PX", 200)  # 1UBQ's 4.2 in page at 90 dpi is 375 px wide
    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "v.svg")})
    assert small[0] > 300 and max(width(r)) <= 200


def test_draw_topology_draws_like_a_reference(tmp_path):
    """like_source (with like_chain): the note and the structured result's provenance say what it was drawn like;
    like_chain without like_source is a message."""
    from mcp.server.mcpserver.exceptions import ToolError

    hb = str(DATA / "2HHB.cif")
    r = call(
        "draw_topology",
        {"source": hb, "output_path": str(tmp_path / "b.svg"), "chains": ["B"], "like_source": hb, "like_chain": "A"},
    )
    said = "drawn like 2HHB.cif chain A: 7 of 9 elements matched (sequence alignment, TM-score 0.89 over 139 residues)"
    like = [i for i in r.structured_content["provenance"] if i["kind"] == "like"]
    assert len(like) == 1 and like[0]["text"].startswith(said)
    assert said in r.content[0].text
    with pytest.raises(ToolError, match="like_source"):
        call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "like_chain": "A"})


def test_references_are_checked_before_drawing_as_the_command_line_checks_them(tmp_path):
    """Q-A19: a range without res: is drawn with a warning in the note; a reference that cannot be drawn is an error
    result (isError) with the command line's text, said before anything is laid out."""
    from mcp.server.mcpserver.exceptions import ToolError

    t4l, out = str(DATA / "2LZM.cif"), str(tmp_path / "t4l.svg")
    r = call("draw_topology", {"source": t4l, "output_path": out, "domains": {"N": ["A:1-60"]}})
    assert "; warning: --domain N=A:1-60: read A:1-60 as the residue range res:A:1-60" in r.content[0].text
    assert "domain:N" in Path(out).read_text()
    r = call("draw_topology", {"source": t4l, "output_path": out, "domains": {"N": ["res:A:1-60"]}})
    assert "warning" not in r.content[0].text
    for domains, said in [
        ({"N": ["A:x-y"]}, "--domain N=A:x-y: A:x-y is neither a residue range nor an element"),
        ({"N": ["res:B:1-20"]}, "--domain N=res:B:1-20: no protein chain B in the structure; its protein chains are A"),
        ({"N": ["res:A:900-950"]}, "--domain N=res:A:900-950: res:A:900-950 meets no residue of chain A"),
        ({"N": [""]}, "--domain N=: an empty reference"),
        ({"N": []}, "--domain N=: a domain needs elements"),
    ]:
        with pytest.raises(ToolError) as caught:
            call("draw_topology", {"source": t4l, "output_path": str(tmp_path / "x.svg"), "domains": domains})
        assert said in str(caught.value), (domains, str(caught.value))
    with pytest.raises(ToolError, match=r"--domain N=res:B:1-20: chain B is not drawn: --focus A draws"):
        call("draw_topology", {"source": str(DATA / "2HHB.cif"), "assembly": "asu", "focus": "A",
                               "domains": {"N": ["res:B:1-20"]}})  # fmt: skip


def test_elements_are_named_in_every_convention(tmp_path):
    """Q-A12: dep: and beta: names reach draw_topology as on the command line; a label meaning another element in
    the β-scheme is an error result listing the choices."""
    from mcp.server.mcpserver.exceptions import ToolError

    fid = DATA / "fidelity"
    out = str(tmp_path / "mb.svg")
    call("draw_topology", {"source": str(fid / "1MBN.cif.gz"), "output_path": out, "swap": [["dep:A", "dep:E"]]})
    assert Path(out).is_file()
    with pytest.raises(ToolError, match=r"α4 is the figure's α4 \(A 96–102\) and the β-scheme's α4 \(A 106–118\)"):
        call(
            "draw_topology",
            {"source": str(DATA / "1TIM.cif"), "output_path": out, "chains": ["A"], "swap": [["α4", "α1"]]},
        )


def test_draw_topology_returns_the_interaction_lines_drawn(tmp_path):
    """Plan 13c-5: the interaction style keys reach the figure; the structured result lists each line drawn as
    {kind, a, b, distance, atoms, id}, id its gid in the SVG, distance the shortest of the atoms listed; none
    without them."""
    out = tmp_path / "lmb.svg"
    style = {"interactions": "hbond,salt,cation-pi", "interactions_between": "na", "interaction_max": "50"}
    r = call("draw_topology", {"source": str(DATA / "1LMB.cif"), "output_path": str(out), "style": style})
    lines = r.structured_content["interactions"]
    gids = set(re.findall(r'<g id="((?:hbond|ia:[a-z-]+):[^"]+)">', out.read_text()))
    assert lines and {x["id"] for x in lines} == gids and {x["kind"] for x in lines} == {"hbond", "salt", "cation-pi"}
    for x in lines:
        assert set(x) == {"kind", "a", "b", "distance", "atoms", "id"}
        assert x["distance"] == min(a["distance"] for a in x["atoms"]) and 2.0 < x["distance"] <= 6.0
        assert re.fullmatch(r"[^:]+:[A-Za-z]+-?\d+\w*", x["a"]) and re.fullmatch(r"[^:]+:[A-Za-z]+-?\d+\w*", x["b"])
        assert all({"a", "b", "distance", "angle", "note"} == set(a) for a in x["atoms"])
    salt = next(x for x in lines if x["kind"] == "salt")
    assert salt["id"].startswith("ia:salt:") and salt["atoms"][0]["angle"] is None
    plain = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "u.svg")})
    assert plain.structured_content["interactions"] == []


def test_structure_data_returns_the_command_lines_dict(tmp_path, capsys):
    """structure_data (plan 13d-3): the structured result is the dict `foldmap data` writes for the same structure
    and options, the text its JSON; sections pick the sections; output_path writes the JSON, or the nucleic acid as
    dot-bracket, BPSEQ or CT, the same bytes as the command line's; a bad section or suffix is an error."""
    import json

    from mcp.server.mcpserver.exceptions import ToolError

    from foldmap.cli import main
    from foldmap.datafile import DEFAULT, dumps

    path = str(DATA / "1LMB.cif")
    assert main(["data", path]) == 0
    want = json.loads(capsys.readouterr().out)
    r = call("structure_data", {"source": path})
    assert not r.is_error and r.structured_content == want and set(want) == {"foldmap_data", "foldmap", *DEFAULT}
    assert r.content[0].text == dumps(want)
    assert any(f["measure"] == "atoms" for f in want["interfaces"])  # protein-DNA pairs, no interactions asked
    r = call("structure_data", {"source": path, "sections": ["elements", "nucleic_acid"], "chains": ["3", "1", "2"]})
    assert set(r.structured_content) == {"foldmap_data", "foldmap", "elements", "nucleic_acid"}
    assert {e["chain"] for e in r.structured_content["elements"]} == {"3"}
    for ext in ("dbn", "bpseq", "ct", "json"):
        out = tmp_path / f"x.{ext}"
        r = call("structure_data", {"source": path, "output_path": str(out)})
        assert r.content[0].text == f"wrote {out}" and r.structured_content == want
        assert main(["data", path, "-o", str(tmp_path / f"cli.{ext}")]) == 0
        assert out.read_bytes() == (tmp_path / f"cli.{ext}").read_bytes(), ext
    with pytest.raises(ToolError, match="no data section colours"):
        call("structure_data", {"source": path, "sections": ["colours"]})
    with pytest.raises(ToolError, match="must end in"):
        call("structure_data", {"source": path, "output_path": str(tmp_path / "x.csv")})
    with pytest.raises(ToolError, match="nucleic_acid"):
        call("structure_data", {"source": path, "sections": ["elements"], "output_path": str(tmp_path / "y.ct")})


def test_draw_topology_takes_your_data_and_crosslinks_and_says_their_problems(tmp_path):
    """data_path, data_column and crosslinks_path: the figure is coloured by the data and the crosslinks drawn;
    each row not used (and each crosslink longer than its span) is a warning in the note and in data_problems; a
    table none of whose rows names a residue is an error result (isError), as the command line refuses it."""
    from mcp.server.mcpserver.exceptions import ToolError

    data, xl, out = tmp_path / "d.csv", tmp_path / "xl.csv", tmp_path / "ubq.svg"
    data.write_text("chain,residue,ddG\nA,3,0.5\nA,30,-1.2\nA,999,1\n", encoding="utf-8")
    xl.write_text("chain1,residue1,chain2,residue2,distance\nA,6,A,48,\nA,1,A,72,5\n", encoding="utf-8")
    r = call(
        "draw_topology", {"source": UBQ, "output_path": str(out), "data_path": str(data), "crosslinks_path": str(xl)}
    )
    assert not r.is_error
    problems = r.structured_content["data_problems"]
    assert [p["where"] for p in problems] == [f"--data {data} row 4", f"--crosslinks {xl}"]
    assert "chain A has no residue 999" in problems[0]["text"] and "longer than its 5 Å span" in problems[1]["text"]
    assert f"; warning: --data {data} row 4: chain A has no residue 999" in r.content[0].text
    kinds = {x["id"]: x["kind"] for x in r.structured_content["interactions"]}
    assert kinds == {"ia:xl:A:6-A:48": "xl", "ia:xl:A:1-A:72": "xl-over"}
    svg = out.read_text()
    assert 'id="data-stripe:' in svg and 'id="ia:xl:A:6-A:48"' in svg
    r = call("draw_topology", {"source": UBQ, "output_path": str(out), "data_path": str(data), "data_column": "ddG"})
    assert not r.is_error
    bad = tmp_path / "bad.csv"
    bad.write_text("chain,residue,ddG\nA,900,1\n", encoding="utf-8")
    with pytest.raises(ToolError, match=re.escape(f"--data {bad}: no row names a residue of the structure")):
        call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "data_path": str(bad)})
    with pytest.raises(ToolError, match="needs data"):
        call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "x.svg"), "data_column": "ddG"})


def test_diff_structures_returns_the_preview_what_changed_and_h5(tmp_path):
    """diff_structures (plan 13e-2): writes the diff page and its JSON and returns the note, a PNG preview and what
    changed (the JSON's dict and the provenance), with check H5 of each panel; a bad domain is an error."""
    import json

    a, b = (str(DATA / "diff" / f"{p}.cif.gz") for p in ("4AKE", "1AKE"))
    out, js = tmp_path / "diff.png", tmp_path / "diff.json"
    r = call(
        "diff_structures",
        {
            "source_a": a,
            "source_b": b,
            "domains": {"LID": ["res:A:122-159"]},
            "output_path": str(out),
            "json_path": str(js),
            "check": True,
        },
    )
    assert not r.is_error and [c.type for c in r.content] == ["text", "image"]
    text = r.content[0].text
    assert out.is_file() and str(out) in text and "changes from 4AKE.cif.gz to 1AKE.cif.gz" in text
    assert "H5 left panel: pass" in text and "H5 right panel: pass" in text
    data = r.structured_content
    assert data["alignment"]["method"] == "sequence" and data["domain_turns"]["LID"] > 30
    assert data["only_a"] and data["displacement"] and set(data["h5"]) == {"a", "b"}
    assert any(i["kind"] == "diff" for i in data["provenance"])
    assert json.loads(js.read_text(encoding="utf-8")) == {
        k: v for k, v in data.items() if k not in ("provenance", "h5")
    }
    from mcp.server.mcpserver.exceptions import ToolError

    with pytest.raises(ToolError, match="domain X: 'B:1-40' is not a residue range"):
        call("diff_structures", {"source_a": a, "source_b": b, "domains": {"X": ["B:1-40"]}})


def test_diff_structures_refuses_a_slash_in_a_domain_name():
    """domain_angles names FIRST/NAME: a domain name holding "/" is refused before anything is drawn."""
    from mcp.server.mcpserver.exceptions import ToolError

    ubq = str(DATA / "1UBQ.cif")
    with pytest.raises(ToolError, match="cannot hold /"):
        call("diff_structures", {"source_a": ubq, "source_b": ubq, "domains": {"N/C": ["1-40"]}})


def test_draw_topology_says_where_a_duplex_is_bent(tmp_path):
    """Plan 25 Task 3: style {"na_bend": "curved"} draws IHF's bent duplex as one curve; the structured result says
    where a duplex is drawn bent (the choice applies), whichever way it is drawn."""
    out = tmp_path / "ihf.svg"
    r = call(
        "draw_topology",
        {"source": str(DATA / "na" / "1IHF.cif.gz"), "output_path": str(out), "style": {"na_bend": "curved"}},
    )
    assert r.structured_content["bent_duplex"] is True and 'id="na:curve:0:rails"' in out.read_text()
    r = call("draw_topology", {"source": UBQ, "output_path": str(tmp_path / "ubq.svg")})
    assert r.structured_content["bent_duplex"] is False
