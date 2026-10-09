"""The MCP server: tools agents use to summarise structures and draw figures."""

import asyncio
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from foldmap.mcp_server import server  # noqa: E402

DATA = Path(__file__).parent / "data"


def call(name, args):
    return asyncio.run(server.call_tool(name, args))


def test_tools_are_listed_with_descriptions():
    tools = asyncio.run(server.list_tools())
    names = {t.name for t in tools}
    assert names == {"summarize_structure", "list_styles", "draw_topology", "draw_sequence", "interactive_page"}
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
