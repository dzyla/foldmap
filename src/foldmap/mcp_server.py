# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""foldmap as an MCP server, so agents can draw protein topology diagrams and sequence figures.

Run with `foldmap mcp` (stdio). Structures are given as a file path, a PDB ID (e.g. 1LMB) or a UniProt
accession (e.g. P04637, for the AlphaFold model). Figures are written to disk; the drawing tools also return a
small PNG preview so the agent can look at the result."""

from __future__ import annotations

import contextlib
import functools
import os
import shutil
import tempfile
from pathlib import Path

from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import CallToolResult, TextContent

from . import CREDIT, __version__
from .fetch import fetch

INSTRUCTIONS = (
    """foldmap draws publication-grade protein topology diagrams (2D maps of helices, strands and loops) and
sequence figures with secondary structure on top. Start with summarize_structure to see the elements and their
labels (A, B... strands; α1, α2... helices; η1... 3-10 helices), then draw_topology. Refer to elements by label,
chain:label, #k (k-th element), residue range res:CHAIN:FIRST-LAST, or a name in another convention: dep:NAME (the
deposited helix name, e.g. dep:E in a globin) or beta:NAME (the β-scheme's, e.g. beta:α4); a label that means another
element in one of those is refused with the choices (fig:NAME is the label drawn). Call list_styles for themes and
style keys. structure_data gives what the figure says (elements, sheets, contacts, interfaces, a nucleic acid's base
pairs as dot-bracket) as data. compare_structures draws two structures of one fold side by side on one map;
diff_structures draws what changed between two states of one protein on that page.
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
        except KeyError as err:  # str(KeyError) is the repr of its key: quotes included
            raise ToolError(str(err.args[0]) if err.args else "not found") from err
        except (ValueError, FileNotFoundError, RuntimeError, OSError) as err:
            raise ToolError(str(err)) from err

    return server.tool()(wrapped)


def _structure(source: str) -> Path:
    return fetch(source)


FIGURE_SUFFIXES = (".svg", ".png", ".pdf")
PAGE_SUFFIXES = (".html",)
LAYOUT_SUFFIXES = (".yaml", ".yml")
JSON_SUFFIXES = (".json",)
DATA_SUFFIXES = (".json", ".dbn", ".bpseq", ".ct")


def _checked(path: str, allowed: tuple[str, ...]) -> Path:
    """The path asked for, if it has an allowed extension and lies inside FOLDMAP_OUTPUT_DIR (when that is set)."""
    out = Path(path).expanduser()
    if out.suffix.lower() not in allowed:
        names = ", ".join(allowed[:-1]) + " or " + allowed[-1] if len(allowed) > 1 else allowed[0]
        raise ValueError(f"{out.name!r} must end in {names}")
    root = os.environ.get("FOLDMAP_OUTPUT_DIR")
    if root and not out.resolve().is_relative_to(Path(root).expanduser().resolve()):
        raise ValueError(f"{out} is outside FOLDMAP_OUTPUT_DIR ({root}); write inside it")
    return out


@contextlib.contextmanager
def _output(path: str | None, suffix: str, stem: str, allowed: tuple[str, ...]):
    """Where to write: the path asked for (checked by _checked), or a file in a fresh temporary folder (inside
    FOLDMAP_OUTPUT_DIR when set). A temporary folder is removed again if the call fails."""
    if path:
        yield _checked(path, allowed)
        return
    root = os.environ.get("FOLDMAP_OUTPUT_DIR")
    folder = Path(root).expanduser() if root else None
    if folder:
        folder.mkdir(parents=True, exist_ok=True)
    made = Path(tempfile.mkdtemp(prefix="foldmap-", dir=folder))
    try:
        yield made / f"{stem}{suffix}"
    except BaseException:
        shutil.rmtree(made, ignore_errors=True)
        raise


def _theme(name: str):
    from .style import THEMES

    if name not in THEMES:
        raise ValueError(f"unknown theme {name!r}; choose from {', '.join(THEMES)}")
    return THEMES[name]


PREVIEW_DPI = 90  # the preview image's resolution ...
PREVIEW_PX = 1600  # ... with its longer side at most this many pixels (a page grown for 6 pt text: 48 in)


def _preview(fig) -> Image:
    from .render import save

    dpi = min(PREVIEW_DPI, PREVIEW_PX / max(fig.get_size_inches()))
    with tempfile.TemporaryDirectory() as tmp:
        data = save(fig, Path(tmp) / "preview.png", dpi=dpi).read_bytes()
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
        lines.append(f"  {label}: {s.chain} {s.first.tag}-{s.last.tag} ({kinds.get(s.kind, s.kind)}, {len(s)} res)")
    links = layout.links
    if links is not None:
        for i, j in links.disulfides:
            a, b = bb.labels[i], bb.labels[j]
            lines.append(f"disulfide: {a.chain} Cys{a.tag} - {b.chain} Cys{b.tag}")
        for r, sugars in links.glycans:
            lines.append(f"glycan on {bb.labels[r].chain} {bb.labels[r].name}{bb.labels[r].tag}: {'-'.join(sugars)}")
        for g in links.ligands:
            held = ", ".join(
                [f"{bb.labels[r].name}{bb.labels[r].tag}" for r in g.contacts[:8]]
                or [f"nucleotide {c}:{n}{i}" for c, n, i in g.nucleotides[:8]]  # held by the nucleic acid alone
            )
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
    theme: str | None = None,
    style: dict[str, str] | None = None,
    title: str | None = None,
    assembly: str | None = None,
    symmetry: str | None = None,
    membrane: str | None = None,
    focus: str | None = None,
    chains: list[str] | None = None,
    uniprot: str | None = None,
    uniprot_domains: bool = False,
    domains: dict[str, list[str]] | None = None,
    domains_from: str | None = None,
    swap: list[list[str]] | None = None,
    move: dict[str, list[float]] | None = None,
    rename: dict[str, str] | None = None,
    colors: dict[str, str] | None = None,
    widths: dict[str, float] | None = None,
    names: list[str] | None = None,
    msa_path: str | None = None,
    rotate: float | None = None,
    mode: str | None = None,
    layout_file: str | None = None,
    save_layout: str | None = None,
    like_source: str | None = None,
    like_chain: str | None = None,
    like_match: str | None = None,
    data_path: str | None = None,
    data_column: str | None = None,
    crosslinks_path: str | None = None,
) -> CallToolResult:
    """Draw the topology diagram and write it to output_path (.svg, .png or .pdf; default: a temporary .svg).
    Returns where it was written plus a PNG preview of the same figure; when assembly or focus auto drew part of
    the file (one copy of a complex, one subunit of a large assembly), the note says so and which arguments draw
    every chain. The note also says how the figure was made: what in it is inferred (automatic domains), ambiguous
    (strands drawn off a sheet's run) or moved for readability; the structured result's `provenance` lists the same
    as {kind, text, count, elements} (empty when nothing beyond the standard method), and its `interactions` every
    interaction line drawn as {kind, a, b, distance, atoms, id} (empty unless style asks for them: interactions
    "all" or kinds like "hbond,salt,stack,cation-pi,hydrophobic", interactions_between "auto", "na", "interchain"
    or "chains:A+B", interaction_max; the older hbonds key draws H-bond lines only).

    Options left out come from layout_file (if given), otherwise the defaults: theme publication, assembly auto,
    symmetry auto, membrane off, focus auto, mode projected, rotate 0.
    style: style keys to change, e.g. {"color_by": "sequence", "residue_numbers": "true", "mark": "α2=#d1495b"};
    {"na_bend": "curved"} draws a bent duplex (nucleosome, IHF, TBP) as one curved ladder with the protein laid out
    round it, instead of straight pieces (the result's `bent_duplex` is true where there is one).
    domains: named domain panels, e.g. {"N-lobe": ["res:A:13-59"]}. domains_from: "auto" (inferred from residue
    contacts, not a published assignment: the note says so) or "cath" (CATH's published chopping of a PDB entry,
    fetched from PDBe; the note names its domains, or says there are none). swap: pairs of elements to exchange, e.g.
    [["α1", "α3"]]. move: nudges in page units, e.g. {"α2": [1, -2]}. rename: {"res:A:167-182": "Gd"}. colors: each
    one's own colour, winning over the colouring (the legend says colours were set by hand), e.g. {"α1": "red",
    "res:A:20-25": "#1f77b4"}; widths: each one's drawn width times a factor from 0.25 to 4 (a helix's diameter, an
    arrow's width, a loop's stroke, a backbone's weight; the layout keeps the room), e.g. {"α1": 2}; there a
    reference may also be loop:REF (the loop after element REF), na:CHAIN (a nucleic-acid backbone) or lig:CODE (a
    ligand's box, lig:HEM or lig:HEM:A142). Every
    element reference is checked before drawing: one that cannot be drawn (a range that does not parse, a chain the
    structure lacks or does not draw, residues it does not have) is an error saying so; a range written without res:
    (A:1-60) is read as res:A:1-60 and the note warns that it was. An element may also be named as "names" would
    name it, without applying them: dep:NAME (deposited, e.g. dep:E) or beta:NAME (β-scheme, e.g. beta:α4); a label
    meaning another element in one of those is an error listing the choices (fig:NAME is the label drawn).
    names: conventional names, in order (later win, rename wins over all): "deposited" (the file's own helix ids,
    e.g. globin A-H), "beta" (strands β1…βn, αn after βn in a βα row) or a YAML file of element and stem names.
    membrane: off, auto (draw the lipid bilayer) or flip (the same, the two sides swapped). focus: auto draws
    filaments, cages and large assemblies as one subunit plus the neighbouring parts of its fold; none draws the
    whole assembly; or chains, e.g. "B".
    chains: draw only these chains, e.g. ["A", "B"].
    msa_path: alignment for theme 'conservation'. uniprot: "auto" or an accession, to mark UniProt active and binding
    sites (style key uniprot_sites: all, or types like "active,modified"); uniprot_domains: UniProt domain panels.
    layout_file / save_layout: load or write a YAML layout file that re-creates the figure
    (save_layout must end in .yaml or .yml).
    like_source: draw this structure like a reference of the same fold (a file path or PDB ID): each element matched
    to one of the reference's takes its place, angle, name and (coloured by identity) colour. One protein chain of
    each: `chains` names this structure's, like_chain the reference's (else each its first protein chain).
    like_match: auto (sequence at 40% identity, else structure), sequence, structure, or a two-row alignment file
    (the reference's chain first). The note, and the `provenance` item of kind "like", say how many elements
    matched, or why it was drawn on its own.
    data_path: your own per-residue values, a table (CSV, TSV) with chain, residue and value columns; the figure is
    coloured by them (color_by data) unless style asks for another colouring: numbers as a stripe per residue along
    each element (style data_style mean: one colour per element), words (categories) as a mark on each residue.
    data_column: the column read (default: the first column of numbers). crosslinks_path: a table with chain1,
    residue1, chain2, residue2 and an optional distance, the crosslinker's span in Å (default 30 Å, DSS/BS3):
    each crosslink is drawn as a solid line between its residues, or dashed when its Cα–Cα distance in the structure
    is longer than its span. A file none of whose rows names a residue of the structure (or that lacks a needed
    column) is an error; each row not used, and each crosslink longer than its span, is a warning in the note and
    in the structured result's `data_problems` [{where, text}]; crosslinks drawn are among `interactions` (kind xl,
    or xl-over when longer than their span)."""
    from . import layoutfile
    from .provenance import lines
    from .render import save
    from .scope import drawn_scope, scope_note
    from .spec import check_references, draw_spec, figure_spec, layout_document
    from .userdata import SAID_ROWS

    if domains_from not in (None, "auto", "cath"):
        raise ValueError(f"domains_from is auto or cath; got {domains_from!r}")
    if not like_source and (like_chain or like_match):
        raise ValueError("like_chain and like_match need like_source (the reference structure)")
    like = {"structure": like_source, "chain": like_chain, "match": like_match} if like_source else None
    path = _structure(source)
    keep = _checked(save_layout, LAYOUT_SUFFIXES) if save_layout else None  # checked first
    with _output(output_path, ".svg", path.stem, FIGURE_SUFFIXES) as out:
        for ref, d in (move or {}).items():
            if len(d) != 2:
                raise ValueError(f"move needs two numbers [dx, dy] for {ref!r}, got {d!r}")
        spec = figure_spec(
            path,
            theme=theme,
            style=[f"{k}={v}" for k, v in (style or {}).items()],
            layout_file=layout_file,
            title=title,
            mode=mode,
            rotate=rotate,
            assembly=assembly,
            symmetry=symmetry,
            membrane=membrane,
            focus=focus,
            chains=chains,
            uniprot=uniprot,
            domains=[(n, list(r)) for n, r in (domains or {}).items()]
            or domains_from
            or ("uniprot" if uniprot_domains else None),
            msa=msa_path,
            rename=rename,
            names=names,
            swap=[tuple(p) for p in swap or []],
            move=[(ref, (d[0], d[1])) for ref, d in (move or {}).items()],
            colors=[f"{ref}={value}" for ref, value in (colors or {}).items()],
            widths=[f"{ref}={value}" for ref, value in (widths or {}).items()],
            like=like,
            data=data_path,
            data_column=data_column,
            crosslinks=crosslinks_path,
        )
        warned = []  # element references checked before the structure is laid out (Q-A19), as the command line does
        for problem in check_references(spec):
            if problem.fatal:
                raise ValueError(str(problem))
            warned.append(problem)
        drawn = draw_spec(spec)
        save(drawn.figure, out)
        note = f"wrote {out}" + "".join(f"; warning: {problem}" for problem in warned)
        if keep is not None:
            note += f" and layout file {layoutfile.save(keep, layout_document(spec, drawn=drawn))}"
        problems = list(drawn.layout.data_problems)  # your data: each row not used, each crosslink too long
        note += "".join(f"; warning: {p}" for p in problems[:SAID_ROWS])
        if len(problems) > SAID_ROWS:
            note += f"; and {len(problems) - SAID_ROWS} more data warnings (data_problems lists each)"
        items = drawn.provenance  # how it was made: the note says it, the structured result lists it
        note += "".join(f"; {text}" for text in lines(items))
        scope = drawn_scope(drawn.model, drawn.layout, spec.assembly, spec.focus)
        said = scope_note({**scope, "additives": []}, how="mcp")  # the additives hidden are in the list
        if said:  # an automatic choice drew part of the file
            note += f"; {said}"
        return CallToolResult(
            content=[TextContent(type="text", text=note), _preview(drawn.figure).to_image_content()],
            structured_content={
                "provenance": items,
                "interactions": _interaction_lines(drawn.figure),
                "data_problems": [{"where": p.where, "text": p.text} for p in problems],
                "bent_duplex": drawn.layout.bent_duplex,  # na_bend (pieces or curved) is a choice for it
                "fallbacks": drawn.fallbacks,  # each named fallback taken -> how many times (fallback.FALLBACKS)
            },
        )


def _interaction_lines(fig) -> list[dict]:
    """The interaction lines a figure draws (style keys interactions / hbonds), one per line: {kind, a, b,
    distance, atoms, id}: a and b its residues (A:Arg19, 1:DT7), distance the shortest of the pair (Å), atoms each
    interaction behind it {a, b, distance, angle, note}, id the line's id in the SVG (hbond:..., ia:<kind>:...)."""
    overlay = getattr(fig, "foldmap_interactions", None)
    keep = ("kind", "a", "b", "distance", "atoms")
    return [{**{k: x[k] for k in keep}, "id": x["key"]} for x in (overlay.records() if overlay else [])]


@tool
def compare_structures(
    source_a: str,
    source_b: str,
    chain_a: str | None = None,
    chain_b: str | None = None,
    match: str = "auto",
    output_path: str | None = None,
    json_path: str | None = None,
    assembly: str = "auto",
) -> CallToolResult:
    """Draw two structures of one fold on one page: source_a (the reference) on the left, source_b drawn like it on
    the right (each element matched to one of A's takes its place, angle, name and colour), at one page scale, with
    one legend and the note under both. Writes output_path (.svg, .png or .pdf; default: a temporary .svg) and
    returns where, the note and a PNG preview; the structured result is what was matched: alignment {method,
    tm_score, rmsd, identity, pairs ([A's residue, B's residue])}, mates, pieces (of a helix A has whole), extras
    (B's elements A lacks), only_a (A's elements B lacks), each element by its label and residues, plus the figure's
    provenance. json_path (.json) writes the same, without provenance, to a file.

    source_a, source_b: file path, PDB ID or UniProt accession. One protein chain of each: chain_a, chain_b (else each
    its first). match: auto (sequence at 40% identity, else structure), sequence, structure, or a two-row alignment
    file (A's chain first). assembly: auto, asu (each file as deposited, under its deposited chain names) or an
    assembly id. Not the same fold: B is drawn on its own beside A, and the note says so."""
    from .compare import compare, write_json
    from .provenance import lines
    from .render import save

    keep = _checked(json_path, JSON_SUFFIXES) if json_path else None  # checked first
    with _output(output_path, ".svg", f"{Path(str(source_b)).stem}-compare", FIGURE_SUFFIXES) as out:
        made = compare(source_a, source_b, chain_a, chain_b, match or "auto", assembly=assembly or "auto")
        save(made.figure, out)
        note = f"wrote {out}"
        if keep is not None:
            note += f" and {write_json(made.data, keep)}"
        note += "".join(f"; {text}" for text in lines(made.provenance))
        return CallToolResult(
            content=[TextContent(type="text", text=note), _preview(made.figure).to_image_content()],
            structured_content={**made.data, "provenance": made.provenance},
        )


@tool
def diff_structures(
    source_a: str,
    source_b: str,
    chain_a: str | None = None,
    chain_b: str | None = None,
    match: str = "auto",
    domains: dict[str, list[str]] | None = None,
    output_path: str | None = None,
    json_path: str | None = None,
    assembly: str = "auto",
    check: bool = False,
) -> CallToolResult:
    """What changed between two structures of one protein (two states: apo and holo, open and closed), drawn on
    their compare page (compare_structures' page: source_a on the left, source_b drawn like it on the right): an
    element only one structure has as a dashed outline in the other's panel; source_b's panel coloured by each
    residue's CA displacement after superposition on the matched core (when any residue moved); and, in the note under
    both, secondary-structure changes, sheet pairing gained, lost or shifted, domain turns, interface and contact
    changes. Writes output_path (.svg, .png or .pdf; default: a temporary .svg) and returns where, the note and a PNG
    preview; the structured result is what changed: compare_structures' matching (alignment, mates, pieces, extras,
    only_a), then ss_changes [{chain, residue, a, b}] (kinds H, G, E, C none, - absent), pairing [{strands, change}],
    displacement [[chain, residue, angstroms]] (source_b's residues), core (the number of pairs the core fit kept),
    domain_turns {name: degrees} (domains given: each domain's turn after the core fit), domain_angles {"FIRST/NAME":
    degrees} (two domains or more: each after the first against it, each fitted on its own), contacts {elements,
    interfaces, interface_elements, ligands, ligand_elements: [[..., gained | lost]]} (partner chains named by
    source_a's chain in their place) and partners {source_b's partner chain: the name contacts give it}, plus the
    figure's provenance.
    json_path (.json) writes the same, without provenance.

    source_a, source_b: file path, PDB ID or UniProt accession. One protein chain of each: chain_a, chain_b (else
    each its first). match: auto (sequence when the chains are one protein, else structure), sequence, structure, or
    a two-row alignment file (A's chain first). domains: {name: residue ranges of source_a's chain, e.g.
    ["res:A:122-159"]}, each domain's rigid turn given (a name holds no "/"). assembly: auto, asu or an assembly id.
    check: also judge each panel as a figure of its own (H5: it passes every check its chain drawn alone passes), as
    h5 {a, b}."""
    from .compare import diff_page, panel_checks, write_json
    from .provenance import lines
    from .render import save

    slashed = [name for name in domains or {} if "/" in name]
    if slashed:  # domain_angles names FIRST/NAME
        raise ValueError(f"domain {slashed[0]!r}: a domain name cannot hold / (domain_angles names FIRST/NAME)")
    keep = _checked(json_path, JSON_SUFFIXES) if json_path else None  # checked first
    with _output(output_path, ".svg", f"{Path(str(source_b)).stem}-diff", FIGURE_SUFFIXES) as out:
        made = diff_page(
            source_a,
            source_b,
            chain_a,
            chain_b,
            match or "auto",
            assembly=assembly or "auto",
            domains=dict(domains) if domains else None,
        )
        save(made.figure, out)
        note = f"wrote {out}"
        if keep is not None:
            note += f" and {write_json(made.data, keep)}"
        note += "".join(f"; {text}" for text in lines(made.provenance))
        result = {**made.data, "provenance": made.provenance}
        if check:
            result["h5"] = panel_checks(made, assembly or "auto")
            note += "".join(
                f"; H5 {side} panel: {'not judged' if got['pass'] is None else 'pass' if got['pass'] else 'fail'}"
                for side, got in zip(("left", "right"), result["h5"].values())
            )
        return CallToolResult(
            content=[TextContent(type="text", text=note), _preview(made.figure).to_image_content()],
            structured_content=result,
        )


@tool
def structure_data(
    source: str,
    sections: list[str] | None = None,
    output_path: str | None = None,
    assembly: str | None = None,
    symmetry: str | None = None,
    membrane: str | None = None,
    focus: str | None = None,
    chains: list[str] | None = None,
    domains: dict[str, list[str]] | None = None,
    domains_from: str | None = None,
    rename: dict[str, str] | None = None,
    names: list[str] | None = None,
    style: dict[str, str] | None = None,
    layout_file: str | None = None,
) -> CallToolResult:
    """What the topology figure says about a structure, as data (the command line's `foldmap data`; JSON schema
    version 1): the structured result is the whole dict, the text the same as JSON (or, with output_path, where it
    was written). The structure is laid out exactly as draw_topology draws it with the same arguments, so names
    and labels are the figure's.

    sections (default all but interactions): structure; elements (kind, residues, label, page place); sheets (strand
    order as drawn, directions, each pairing's kind and residue register from the DSSP bridges, accessory strands,
    a barrel's n and shear S); disulfides; symmetry (and stoichiometry); membrane (with membrane="auto"); ligands
    (ligands, metal ions, glycans, fused residues: the residues and elements they touch); interfaces (chain pairs in
    contact: protein pairs by CA pairs, protein-nucleic acid by atom contacts, nucleic-acid pairs by base pairs;
    element pairs); nucleic_acid (strands, sequences, the base pairs the figure draws with their stem and
    pseudoknot bracket level, the dot-bracket); interactions (every interaction between two chains, only when
    named). output_path: .json writes the dict; .dbn, .bpseq or .ct the nucleic acid as dot-bracket, BPSEQ or CT.
    The other arguments are draw_topology's."""
    from .datafile import DEFAULT, NA_FORMATS, SECTIONS, dumps, spec_data
    from .spec import check_references, figure_spec

    if domains_from not in (None, "auto", "cath"):
        raise ValueError(f"domains_from is auto or cath; got {domains_from!r}")
    asked = tuple(sections) if sections else DEFAULT
    unknown = [x for x in asked if x not in SECTIONS]
    if unknown:
        raise ValueError(f"no data section {', '.join(unknown)}; the sections are {', '.join(SECTIONS)}")
    out = _checked(output_path, DATA_SUFFIXES) if output_path else None  # checked first
    form = out.suffix.lower().lstrip(".") if out is not None else "json"
    if form != "json" and "nucleic_acid" not in asked:
        raise ValueError(f"a .{form} file is the nucleic_acid section, which sections leaves out")
    path = _structure(source)
    spec = figure_spec(
        path,
        style=[f"{k}={v}" for k, v in (style or {}).items()],
        layout_file=layout_file,
        assembly=assembly,
        symmetry=symmetry,
        membrane=membrane,
        focus=focus,
        chains=chains,
        domains=[(n, list(r)) for n, r in (domains or {}).items()] or domains_from,
        rename=rename,
        names=names,
    )
    warned = []
    for problem in check_references(spec):
        if problem.fatal:
            raise ValueError(str(problem))
        warned.append(problem)
    data = spec_data(spec, asked)
    text = dumps(data)
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text if form == "json" else NA_FORMATS[form](data))
        text = f"wrote {out}"
    text += "".join(f"; warning: {problem}" for problem in warned)
    return CallToolResult(content=[TextContent(type="text", text=text)], structured_content=data)


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
    data_path: str | None = None,
    data_column: str | None = None,
) -> list:
    """Draw the sequence with its secondary structure on top (labels and colours as in the topology figure), in
    rows of `columns`. With msa_path (FASTA, Clustal or Stockholm), the whole alignment is shown, ESPript-style,
    with the structure following its reference row. data_path, data_column: your own per-residue values (as
    draw_topology takes them), drawn on the track as the topology figure draws them. Writes output_path
    (.svg/.png/.pdf); returns a preview."""
    from .render import save
    from .seqplot import draw_sequence as draw
    from .sequence import read_alignment
    from .spec import data_look

    if data_column and not data_path:
        raise ValueError(f"data_column {data_column!r} needs data_path: the file it is a column of")
    path = _structure(source)
    with _output(output_path, ".svg", f"{path.stem}-sequence", FIGURE_SUFFIXES) as out:
        fig = draw(
            path,
            chain=chain,
            alignment=read_alignment(msa_path) if msa_path else None,
            reference=reference,
            columns=columns,
            look=data_look(_theme(theme), data_path, theme),
            full_sequence=full_sequence,
            ss_colour=ss_colour,
            title=title,
            msa=msa_path,
            data=data_path,
            data_column=data_column,
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
    """Write an HTML page (output_path must end in .html) linking the topology, a residue contact map and the 3D model
    (hover anything: the same residues light up in all three). viewer: molstar (Mol*, default) or 3dmol.
    Returns the path."""
    from .interactive import write_page

    path = _structure(source)
    with _output(output_path, ".html", f"{path.stem}-explorer", PAGE_SUFFIXES) as out:
        return f"wrote {write_page(path, out, look=_theme(theme), assembly=assembly, viewer=viewer)}"


def main() -> None:
    server.run()  # stdio


if __name__ == "__main__":
    main()
