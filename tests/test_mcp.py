"""The MCP server: tools agents use to summarise structures and draw figures."""

import asyncio
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from foldmap.mcp_server import server  # noqa: E402

DATA = Path(__file__).parent / "data"
UBQ = str(DATA / "1UBQ.cif")


def call(name, args):
    return asyncio.run(server.call_tool(name, args))


def test_tools_are_listed_with_descriptions():
    tools = asyncio.run(server.list_tools())
    names = {t.name for t in tools}
    assert names == {
        "summarize_structure",
        "list_styles",
        "draw_topology",
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
