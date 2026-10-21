"""The MCP server: tools agents use to summarise structures and draw figures."""

import asyncio
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
        "draw_sequence",
        "interactive_page",
        "uniprot_info",
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
