"""foldmap as an MCP server, so agents can draw protein topology diagrams and sequence figures.

Run with `foldmap mcp` (stdio). Structures are given as a file path, a PDB ID (e.g. 1LMB) or a UniProt
accession (e.g. P04637, for the AlphaFold model). Figures are written to disk; the drawing tools also return a
small PNG preview so the agent can look at the result."""

from __future__ import annotations

import functools
import tempfile
from pathlib import Path

from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import ToolError

from . import CREDIT, __version__
from .fetch import fetch

INSTRUCTIONS = (
    """foldmap draws publication-grade protein topology diagrams (2D maps of helices, strands and loops) and
sequence figures with secondary structure on top. Start with summarize_structure to see the elements and their
labels (A, B... strands; α1, α2... helices; η1... 3-10 helices), then draw_topology. Refer to elements by label,
chain:label, #k (k-th element) or residue range res:CHAIN:FIRST-LAST. Call list_styles for themes and style keys.
"""
    + CREDIT
)

server = MCPServer("foldmap", instructions=INSTRUCTIONS, version=__version__)


def tool(fn):
    """Register fn as a tool whose input problems come back to the agent as readable messages."""

    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (ValueError, FileNotFoundError, KeyError) as err:
            raise ToolError(str(err).strip("'\"")) from err

    return server.tool()(wrapped)


def _structure(source: str) -> Path:
    return fetch(source)


def _output(path: str | None, suffix: str, stem: str) -> Path:
    if path:
        out = Path(path).expanduser()
        if out.suffix.lower() not in (".svg", ".png", ".pdf", ".html"):
            raise ValueError(f"output_path must end in .svg, .png or .pdf (got {out.name!r})")
        return out
    return Path(tempfile.mkdtemp(prefix="foldmap-")) / f"{stem}{suffix}"


def _preview(fig) -> Image:
    from .render import save

    with tempfile.TemporaryDirectory() as tmp:
        data = save(fig, Path(tmp) / "preview.png", dpi=90).read_bytes()
    return Image(data=data, format="png")


@tool
def summarize_structure(source: str, assembly: str = "auto") -> str:
    """Describe a structure: chains, secondary structure elements with the labels the figures use and their
    residue ranges, sheets, symmetry, DNA/RNA, disulfides, glycans and bound ligands.

    source: file path, PDB ID (1LMB) or UniProt accession (P04637). assembly: auto, asu or an assembly id."""
    from .cli import make_layout, summarize

    path = _structure(source)
    layout, sses, bb = make_layout(path, assembly=assembly)
    lines = [summarize(path), "", "Elements (label: chain residues, kind):"]
    kinds = {"H": "helix", "E": "strand", "G": "3-10 helix"}
    for s in sses:
        p = layout.placed.get(s.id)
        label = p.label if p is not None else "?"
        lines.append(f"  {label}: {s.chain} {s.first.seq}-{s.last.seq} ({kinds.get(s.kind, s.kind)}, {len(s)} res)")
    links = layout.links
    if links is not None:
        for i, j in links.disulfides:
            a, b = bb.labels[i], bb.labels[j]
            lines.append(f"disulfide: {a.chain} Cys{a.seq} - {b.chain} Cys{b.seq}")
        for r, sugars in links.glycans:
            lines.append(f"glycan on {bb.labels[r].chain} {bb.labels[r].name}{bb.labels[r].seq}: {'-'.join(sugars)}")
        for g in links.ligands:
            held = ", ".join(f"{bb.labels[r].name}{bb.labels[r].seq}" for r in g.contacts[:8])
            extra = " (buffer additive, hidden by default)" if g.additive else ""
            lines.append(f"{'metal' if g.metal else 'ligand'} {g.name} {g.chain}{g.seq}: held by {held}{extra}")
    return "\n".join(lines)


@tool
def uniprot_info(source: str) -> str:
    """UniProt annotation for a UniProt accession (P04637), a PDB ID or a structure file: which entry each chain
    is, the protein's names, domains and regions, sites and modifications, its experimental PDB structures (widest
    coverage first, with method and resolution) and its AlphaFold model. To show UniProt sites on a figure call
    draw_topology with uniprot="auto" (or an accession); add uniprot_domains=true for UniProt domain panels."""
    from .uniprot import report

    return report(source)


@tool
def list_styles() -> str:
    """Themes (whole looks) and every style key with its allowed values."""
    from .cli import styles_help

    return styles_help()


@tool
def draw_topology(
    source: str,
    output_path: str | None = None,
    theme: str = "publication",
    style: dict[str, str] | None = None,
    title: str | None = None,
    assembly: str = "auto",
    symmetry: str = "auto",
    membrane: str = "off",
    focus: str = "auto",
    uniprot: str | None = None,
    uniprot_domains: bool = False,
    domains: dict[str, list[str]] | None = None,
    swap: list[list[str]] | None = None,
    move: dict[str, list[float]] | None = None,
    rename: dict[str, str] | None = None,
    msa_path: str | None = None,
    rotate: float = 0.0,
    mode: str = "projected",
    layout_file: str | None = None,
    save_layout: str | None = None,
) -> list:
    """Draw the topology diagram and write it to output_path (.svg, .png or .pdf; default: a temporary .svg).
    Returns where it was written plus a PNG preview.

    style: style keys to change, e.g. {"color_by": "sequence", "residue_numbers": "true", "mark": "α2=#d1495b"}.
    domains: named domain panels, e.g. {"N-lobe": ["res:A:13-59"]}. swap: pairs of elements to exchange, e.g.
    [["α1", "α3"]]. move: nudges in page units, e.g. {"α2": [1, -2]}. rename: {"res:A:167-182": "Gd"}.
    membrane: off or auto (draw the lipid bilayer). focus: auto draws filaments, cages and large assemblies as one
    subunit plus the neighbouring parts of its fold; none draws the whole assembly; or chains, e.g. "B".
    msa_path: alignment for theme 'conservation'. uniprot: "auto" or an accession, to mark UniProt active and binding
    sites (style key uniprot_sites: all, or types like "active,modified"); uniprot_domains: UniProt domain panels.
    layout_file / save_layout: load or write a YAML layout file that re-creates the figure."""
    from .cli import main

    path = _structure(source)
    out = _output(output_path, ".svg", path.stem)
    argv = [
        "plot",
        str(path),
        "-o",
        str(out),
        "--theme",
        theme,
        "--assembly",
        assembly,
        "--symmetry",
        symmetry,
        "--membrane",
        membrane,
        "--focus",
        focus,
        "--mode",
        mode,
        "--rotate",
        str(rotate),
    ]
    for key, value in (style or {}).items():
        argv += ["--set", f"{key}={value}"]
    if uniprot:
        argv += ["--uniprot", uniprot]
        if uniprot_domains:
            argv += ["--domains", "uniprot"]
    for name, refs in (domains or {}).items():
        argv += ["--domain", f"{name}={','.join(refs)}"]
    for pair in swap or []:
        argv += ["--swap", ",".join(pair)]
    for ref, (dx, dy) in (move or {}).items():
        argv += ["--move", f"{ref}={dx},{dy}"]
    for ref, name in (rename or {}).items():
        argv += ["--rename", f"{ref}={name}"]
    if title:
        argv += ["--title", title]
    if msa_path:
        argv += ["--msa", msa_path]
    if layout_file:
        argv += ["--layout-file", layout_file]
    if save_layout:
        argv += ["--save-layout", save_layout]
    _run(main, argv)
    preview = _figure_preview(path, theme, style, title, assembly, symmetry, membrane, msa_path, focus)
    note = f"wrote {out}" + (f" and layout file {save_layout}" if save_layout else "")
    return [note, preview] if preview is not None else [note]


def _run(main, argv: list[str]) -> None:
    import contextlib
    import io

    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        code = main(argv)
    if code != 0:
        raise ValueError(err.getvalue().strip() or "foldmap failed")


def _figure_preview(path, theme, style, title, assembly, symmetry, membrane, msa_path, focus="auto"):
    from .cli import make_figure
    from .style import resolve_style

    try:
        look = resolve_style(theme, None, [f"{k}={v}" for k, v in (style or {}).items()])
        fig = make_figure(
            path,
            title=title,
            look=look,
            assembly=assembly,
            symmetry=symmetry,
            membrane=membrane,
            msa=msa_path,
            focus=focus,
        )
        return _preview(fig)
    except ValueError:
        return None


@tool
def draw_sequence(
    source: str,
    output_path: str | None = None,
    msa_path: str | None = None,
    reference: str | None = None,
    chain: str | None = None,
    columns: int = 60,
    theme: str = "publication",
    ss_colour: str = "figure",
    full_sequence: bool = False,
    title: str | None = None,
) -> list:
    """Draw the sequence with its secondary structure on top (labels and colours as in the topology figure), in
    rows of `columns`. With msa_path (FASTA, Clustal or Stockholm), the whole alignment is shown, ESPript-style,
    with the structure following its reference row. Writes output_path (.svg/.png/.pdf); returns a preview."""
    from .render import save
    from .seqplot import draw_sequence as draw
    from .sequence import read_alignment
    from .style import THEMES

    path = _structure(source)
    out = _output(output_path, ".svg", f"{path.stem}-sequence")
    fig = draw(
        path,
        chain=chain,
        alignment=read_alignment(msa_path) if msa_path else None,
        reference=reference,
        columns=columns,
        look=THEMES[theme],
        full_sequence=full_sequence,
        ss_colour=ss_colour,
        title=title,
    )
    save(fig, out)
    return [f"wrote {out}", _preview(fig)]


@tool
def interactive_page(
    source: str,
    output_path: str | None = None,
    theme: str = "publication",
    assembly: str = "auto",
    viewer: str = "molstar",
) -> str:
    """Write a self-contained HTML page linking the topology, a residue contact map and the 3D model (hover
    anything: the same residues light up in all three). viewer: molstar (Mol*, default) or 3dmol. Returns the path."""
    from .interactive import write_page
    from .style import THEMES

    path = _structure(source)
    out = Path(output_path).expanduser() if output_path else _output(None, ".html", f"{path.stem}-explorer")
    return f"wrote {write_page(path, out, look=THEMES[theme], assembly=assembly, viewer=viewer)}"


def main() -> None:
    server.run()  # stdio


if __name__ == "__main__":
    main()
