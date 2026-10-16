# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The app's sidebar sections: how each settings widget is made (its key, kept while hidden, the setting it names
when the figure fails) and the sections Style, Show, Chains & assembly, Edit layout (with the element names and the
domains) and Conservation & UniProt (with the alignment upload).

Streamlit comes in as an argument (app.main imports it); pandas only inside editor_frame (Streamlit's own
dependency, for the editors' tables)."""

from __future__ import annotations

from collections.abc import MutableMapping

from foldmap.app_edits import (
    EDITORS,
    NAME_KEYWORDS,
    SYMMETRIES,
    UPLOADS,
    apply_table,
    choose_names,
    compose_into,
    display,
    element_choices,
    hex_colour,
    shown,
    switch_domains,
    sync_composites,
    table_rows,
    upload_name,
    upload_types,
    wrong_upload,
)
from foldmap.render import SITE_GROUPS
from foldmap.style import CHOICES, THEMES


def touched(ss: MutableMapping, key: str) -> None:
    """Every settings widget's callback: the setting changed last is the one to check when the figure fails."""
    ss["last_changed"] = key


def setting(ss: MutableMapping, key: str) -> dict:
    """A settings widget's arguments: its key (the setting's own), kept while the widget is hidden, and touched."""
    return {"key": key, "on_change": touched, "args": (ss, key), "persist_state": "session"}


def part(ss: MutableMapping, key: str, canonical: str) -> dict:
    """A composite's widget arguments (`ligands__mode`, ...): it writes the canonical setting (compose_into)."""
    return {"key": key, "on_change": compose_into, "args": (ss, canonical), "persist_state": "session"}


COLOURS = ("background", "ink", "dna_color", "rna_color", "hbond_color")  # settings with a colour picker


def prepare(ss: MutableMapping) -> None:
    """Before the workspace draws: the composites' widget keys derived from their settings (sync_composites), and
    every colour setting as #rrggbb, all a colour picker takes (a layout file may say teal; the figure is the
    same, and is drawn with the value the picker shows)."""
    sync_composites(ss)
    for key in COLOURS:
        if key in ss:
            value = hex_colour(ss[key], getattr(THEMES[ss.get("theme", "publication")], key))
            if ss[key] != value:
                ss[key] = value


def colour(st, label: str, ss: MutableMapping, key: str) -> None:
    """A colour setting's picker (its value made #rrggbb by prepare)."""
    st.color_picker(label, **setting(ss, key))


def segmented(st, label: str, options, kw: dict, **more):
    """A segmented choice that always holds one option, as wide as the sidebar."""
    return st.segmented_control(label, list(options), required=True, width="stretch", **kw, **more)


def style_section(st, ss: MutableMapping) -> None:
    """Style: how the elements and loops are drawn, and the figure's two colours."""
    with st.expander("Style", key="sec_style"):
        with st.container(horizontal=True, key="style_shapes"):  # two short choices a row
            segmented(st, "Fill", CHOICES["fill"], setting(ss, "fill"))
            segmented(st, "Helix shading", CHOICES["helix_shading"], setting(ss, "helix_shading"))
        with st.container(horizontal=True, key="style_angles"):
            segmented(
                st,
                "Helix angle",
                CHOICES["helix_angle"],
                setting(ss, "helix_angle"),
                help="snap: nearest of vertical or horizontal · upright: always vertical · tilted: as in 3D",
            )
            segmented(st, "Loop shape", CHOICES["loops"], setting(ss, "loops"))
        with st.container(horizontal=True, vertical_alignment="bottom", key="style_colours"):
            st.selectbox("Loop colour", CHOICES["loop_color"], **setting(ss, "loop_color"))
            colour(st, "Background", ss, "background")
            colour(st, "Ink", ss, "ink")
        with st.container(horizontal=True, vertical_alignment="center", key="style_loops"):
            st.toggle("Loop arrows, N → C", **setting(ss, "loop_arrows"))
            st.slider("Loop weight", 0.4, 3.0, step=0.1, **setting(ss, "loop_width"))


HBONDS = {"off": "off", "na": "protein–DNA/RNA", "interchain": "interchain", "chains": "chain pairs"}


def show_section(st, ss: MutableMapping) -> None:
    """Show: optional features (bridges, sugars, panels, ligands, the membrane), hydrogen bonds, per-chain backbone
    colours (protein loops too, so always shown) and, with nucleic acid in the structure, how DNA and RNA are drawn
    and the contacts check."""
    with st.expander("Show", key="sec_show"):
        with st.container(horizontal=True, gap="xsmall", key="show_toggles"):
            st.toggle("Disulfides", **setting(ss, "disulfides"))
            st.toggle("Glycans", **setting(ss, "glycans"))
            st.toggle(
                "Sheet panels",
                help="a light panel behind each β-sheet, so a sandwich reads as two sheets",
                **setting(ss, "sheet_panels"),
            )
            st.toggle("3₁₀ helices as boxes", **setting(ss, "helices_310"))
        with st.container(horizontal=True, key="show_choices"):  # the membrane's three beside the ligands' four
            segmented(
                st,
                "Ligands and ions",
                ("auto", "all", "none", "codes"),
                part(ss, "ligands__mode", "ligands"),
                help="auto: bound ligands and ions, without buffer additives · codes: only these",
            )
            segmented(
                st,
                "Membrane",
                ("off", "auto", "flip"),
                setting(ss, "membrane"),
                help="auto: find the bilayer and draw it; flip: the same, other side on top",
            )
        if ss.get("ligands__mode") == "codes":
            st.text_input(
                "Ligand codes",
                placeholder="HEM, ZN",
                label_visibility="collapsed",
                **part(ss, "ligands__codes", "ligands"),
            )
        hbonds_block(st, ss)
        st.text_input(
            "Per-chain colours",
            placeholder="A=#1f77b4, dna=#555555",
            help="KEY=COLOUR, KEY a chain id or protein / dna / rna: its loops (or DNA/RNA backbone) in that colour",
            **setting(ss, "backbone_colors"),
        )
        if ss.get("na_chains"):
            nucleic_block(st, ss)


def hbonds_block(st, ss: MutableMapping) -> None:
    """Hydrogen bonds: which (none, protein with DNA/RNA where the structure has nucleic acid, between chains,
    chosen chain pairs), their colour and the most lines drawn."""
    mode = ss.get("hbonds__mode", "off")
    options = [k for k in HBONDS if k != "na" or ss.get("na_chains") or mode == "na"]  # a layout file's na stays
    segmented(st, "Hydrogen bonds", options, part(ss, "hbonds__mode", "hbonds"), format_func=HBONDS.get)
    if mode == "chains":
        st.text_input(
            "Chain pairs", placeholder="A+B, C+D", label_visibility="collapsed", **part(ss, "hbonds__pairs", "hbonds")
        )
    if mode != "off":
        with st.container(horizontal=True, vertical_alignment="bottom", key="hbond_look"):
            colour(st, "H-bond colour", ss, "hbond_color")
            st.number_input("Most H-bonds drawn", min_value=1, step=50, **setting(ss, "hbond_max"))


def nucleic_block(st, ss: MutableMapping) -> None:
    """DNA and RNA: base letters, the stems' style, what is shown, the backbone's colours and the contacts check."""
    st.markdown("**DNA and RNA**")
    segmented(
        st,
        "Stems",
        CHOICES["na_style"],
        setting(ss, "na_style"),
        help="ladder: rails and base-pair rungs · ribbon: double helix · block: bars · letters",
    )
    segmented(
        st,
        "Shown",
        CHOICES["na_extent"],
        setting(ss, "na_extent"),
        help="auto: everything up to 300 nt, else the stems the protein touches",
    )
    segmented(
        st,
        "Backbone",
        ("modality", "chain", "ink", "colour"),
        part(ss, "na_backbone__mode", "na_backbone"),
        help="modality: DNA and RNA colours · chain: each strand's chain colour · ink · colour: one colour",
    )
    mode = ss.get("na_backbone__mode", "modality")
    with st.container(horizontal=True, vertical_alignment="bottom", key="na_colours"):  # beside the base letters
        st.toggle("Base letters", **setting(ss, "nucleotide_labels"))
        if mode == "modality":
            colour(st, "DNA colour", ss, "dna_color")
            colour(st, "RNA colour", ss, "rna_color")
        elif mode == "colour":
            st.color_picker("Backbone colour", **part(ss, "na_backbone__colour", "na_backbone"))
    st.toggle("Check contacts against the coordinates", **setting(ss, "check_contacts"))


def capped(option, most: int = 24) -> str:
    """A picker's label for a chain or assembly id, which a file may make any length: at most `most` characters."""
    text = str(option)
    return text if len(text) <= most else text[: most - 1] + "…"


def listed(options, current) -> list:
    """A picker's options, then what the setting holds that they lack (a layout file may name anything): a widget
    never holds a value it does not offer."""
    current = [current] if isinstance(current, str) else list(current or [])
    return [*options, *(c for c in current if c not in options)]


HIGHLIGHT = ("none", "asu", "protomer", "chains")
FOCUS = {"auto": "auto", "none": "none", "chains": "these chains in full"}


def chains_section(st, ss: MutableMapping) -> None:
    """Chains & assembly: the chains drawn, the chains in colour (the rest grey), the assembly, the symmetry and
    what a large assembly shows in full. Chains are picked from the structure's own; "as text" types the
    command line's syntax instead."""
    with st.expander("Chains & assembly", key="sec_chains"):
        chains_help = "draw only these chains (protein and nucleic acid each restricted only if some are named)"
        heading(st, "Chains", "chains__as_text", chains_help)
        if ss.get("chains__as_text"):
            st.text_input("Chains", placeholder="A,B", label_visibility="collapsed", **setting(ss, "chains"))
        else:
            st.multiselect(
                "Chains",
                listed(ss.get("all_chains") or [], ss.get("chains__pick")),
                placeholder="All chains",
                format_func=capped,
                label_visibility="collapsed",
                **part(ss, "chains__pick", "chains"),
            )
        highlight_block(st, ss)
        with st.container(horizontal=True, key="assembly_row"):
            st.selectbox(
                "Assembly",
                listed(["auto", "asu", *(ss.get("assemblies") or [])], ss.get("assembly")),
                format_func=capped,
                help="auto: the file's first biological assembly when it adds copies · asu: as deposited · or an "
                "assembly of the file",
                **setting(ss, "assembly"),
            )
            st.selectbox(
                "Symmetry",
                listed(SYMMETRIES, ss.get("symmetry")),
                help="auto: found, every protomer drawn alike · off: ignored · or the symmetry to insist on",
                **setting(ss, "symmetry"),
            )
        segmented(
            st,
            "Large assemblies",
            FOCUS,
            part(ss, "focus__mode", "focus"),
            format_func=FOCUS.get,
            help="auto: filaments, cages and big assemblies as one subunit plus what its neighbours add to its "
            "fold · none: everything",
        )
        if ss.get("focus__mode") == "chains":
            st.multiselect(
                "Chains in full",
                listed(ss.get("protein_chains") or [], ss.get("focus__chains")),
                placeholder="Pick chains",
                format_func=capped,
                label_visibility="collapsed",
                **part(ss, "focus__chains", "focus"),
            )


def highlight_block(st, ss: MutableMapping) -> None:
    """Highlight: the chains in colour, the rest grey (none, the asymmetric unit, one protomer or chains picked),
    or as text."""
    highlight_help = (
        "what stays in colour, the rest in grey · none: everything · asu: the chains as deposited · protomer: the "
        "first protomer · chains: the chains picked"
    )
    heading(st, "Highlight", "highlight__as_text", highlight_help)
    if ss.get("highlight__as_text"):
        st.text_input(
            "Highlight",
            placeholder="none, asu, protomer or A,B",
            label_visibility="collapsed",
            **setting(ss, "highlight"),
        )
    else:
        segmented(st, "Highlight", HIGHLIGHT, part(ss, "highlight__mode", "highlight"), label_visibility="collapsed")
        if ss.get("highlight__mode") == "chains":
            st.multiselect(
                "Highlighted chains",
                listed(ss.get("protein_chains") or [], ss.get("highlight__chains")),
                placeholder="Pick chains",
                format_func=capped,
                label_visibility="collapsed",
                **part(ss, "highlight__chains", "highlight"),
            )


def heading(st, label: str, as_text: str, help: str) -> None:
    """A control's label with its "as text" switch on the right (the control below shows no label of its own)."""
    with st.container(horizontal=True, vertical_alignment="center", key=f"{as_text}_head"):
        st.markdown(label, help=help, width="stretch")
        st.toggle("as text", key=as_text, persist_state="session")


EDITOR_TABS = ["Swap", "Move", "Rename", "Mark", "Domains"]
FOLD_TEMPLATE = {"auto": "βα row", "off": "as projected"}
HELIX_ON_SHEET = {"off": "beside the sheet", "auto": "over the sheet"}
EXAMPLE = {  # the syntax, as text
    "swap": "α1,α3",
    "move": "α2=1,-2",
    "rename": "α1=Gd",
    "mark": "α3=#2ca02c, β1",
    "domains": "N=res:A:1-80\nC=α3,α4",
}
EDITOR_HELP = {
    "swap": "two elements trade places · as text: α1,α3 (several: separate with ;)",
    "move": "shift an element by dx, dy in page units · as text: α2=1,-2",
    "rename": "a new name for an element · as text: α1=Gd or res:A:167-182=Gd",
    "mark": "elements drawn in a colour, or the theme's mark colour when none is given · as text: "
    "G, A:G=#2ca02c, #3, res:150-159=red",
    "domains": "a domain's elements are kept together, drawn on a panel with its name · as text: one per line, "
    "NAME=REF[,REF...], e.g. N=res:A:1-80",
}


def layout_section(st, ss: MutableMapping) -> None:
    """Edit layout: the arrangement, rotation and element sizes, the two fold options, the element names (before the
    editors: a rename wins over them), then the editors (swap, move, rename, mark), each a table of drawn elements
    or, "as text", the command line's syntax, and the domains."""
    with st.expander("Edit layout", key="sec_layout"):
        with st.container(horizontal=True, key="layout_shape"):  # two short choices a row
            segmented(
                st,
                "Arrangement",
                ("projected", "stack"),
                setting(ss, "mode"),
                help="projected: as in 3D, flattened · stack: elements in sequence order",
            )
            segmented(
                st,
                "Parallel α/β sheets",
                CHOICES["fold_template"],
                setting(ss, "fold_template"),
                format_func=FOLD_TEMPLATE.get,
                help="an all-parallel α/β sheet (TIM barrel, Rossmann fold) drawn as one row of strands, its helices "
                "above and below · or as projected",
            )
        segmented(
            st,
            "Helices packed on a sheet",
            CHOICES["helix_on_sheet"],
            setting(ss, "helix_on_sheet"),
            format_func=HELIX_ON_SHEET.get,
            help="a helix packed on a sheet's face stays beside the sheet, or is drawn over it",
        )
        with st.container(horizontal=True, key="layout_rotate"):  # two sliders a row
            st.slider("Rotate (°)", -180.0, 180.0, step=0.5, **setting(ss, "rotate"))
            st.slider("Helix size", 0.5, 2.0, step=0.05, **setting(ss, "helix_scale"))
        with st.container(horizontal=True, key="layout_sizes"):
            st.slider("Strand size", 0.5, 2.0, step=0.05, **setting(ss, "strand_scale"))
            if ss.get("na_chains"):
                st.slider("DNA size", 0.5, 2.0, step=0.05, **setting(ss, "dna_scale"))
        names_block(st, ss)
        number, text = st.column_config.NumberColumn, st.column_config.TextColumn
        offsets = {"dx": number("dx", step=0.5, width=NARROW), "dy": number("dy", step=0.5, width=NARROW)}
        colour_column = text("Colour", validate=r"^(#[0-9a-fA-F]{6}|[a-z]+)?$", width=CELL)
        swap, move, rename, mark, domains = st.tabs(EDITOR_TABS, key="editors")  # one editor in view at a time
        with swap:
            editor(st, ss, "swap", "Swap", {"a": "Element", "b": "with"})
        with move:
            editor(st, ss, "move", "Move", {"element": "Element", **offsets})
        with rename:
            editor(st, ss, "rename", "Rename", {"element": "Element", "name": text("Name", width=CELL)})
        with mark:
            editor(st, ss, "mark", "Mark", {"element": "Element", "colour": colour_column})
        with domains:
            domains_block(st, ss)


NAMES = {"generic": "generic", "deposited": "deposited", "beta": "β1…βn"}  # NAMES_HELP says what each is
NAMES_HELP = (
    "conventional names · generic: foldmap's own labels (α1, A) · deposited: the helix ids the file itself gives "
    "(a globin's A–H) · β1…βn: strands numbered in sequence, each helix of a βα row named after the strand before "
    "it · a names file (YAML) names elements by residue range and nucleic-acid stems by one nucleotide they pair, "
    "and wins over the choice above; a rename wins over both"
)
NAMES_FILE_MAX = 100_000  # bytes: a names file is a few kB


def names_block(st, ss: MutableMapping) -> None:
    """Element names (--names): a keyword, and a names file uploaded (then its name and Remove instead of the
    uploader). A layout file's names file that is not on this machine is asked for under the figure."""
    segmented(
        st,
        "Element names",
        NAMES,
        {"key": "names__kw", "on_change": choose_names, "args": (ss,), "persist_state": "session"},
        format_func=NAMES.get,
        help=NAMES_HELP,
    )
    files = [n for n in ss.get("names") or [] if n not in NAME_KEYWORDS]
    if files:
        with st.container(horizontal=True, vertical_alignment="center", gap="small", key="names_file"):
            st.caption(
                f"Names file: {shown(', '.join(names_file_name(ss, f) for f in files), limit=80)}", width="content"
            )
            st.button("Remove", key="names_clear", type="tertiary", on_click=clear_names, args=(ss,), width="content")
    else:
        st.file_uploader(
            "Names file (YAML)",
            type=upload_types("names"),
            key="names_upload",
            on_change=store_names,
            args=(ss,),
            max_upload_size=1,
            help="elements: {res:A:3-18: A, ...} and stems: {A:1: acceptor, ...} (docs/names has examples)",
        )
    if ss.get("names_error"):
        st.error(shown(ss["names_error"]))


def names_file_name(ss: MutableMapping, path: str) -> str:
    """What a names file is called: the name it was uploaded with, else its file name."""
    kept, name = ss.get("names_uploaded") or ("", "")
    return name if path == kept else path.replace("\\", "/").rsplit("/", 1)[-1]


def store_names(ss: MutableMapping) -> None:
    """The names uploader's callback: an uploaded names file, untrusted, is size-capped and read with
    names.parse_names (yamlfile.read_yaml: no anchors or aliases; a mapping of references to names) before it is kept
    (fetch.keep_upload) and used: the names are the keyword, then this file. Else names_error, and nothing changes."""
    from foldmap.fetch import keep_upload
    from foldmap.names import parse_names

    ss["last_changed"], ss["names_error"] = "names", ""
    up = ss.get("names_upload")
    if up is None:
        return
    if wrong_upload(up.name, "names"):
        ss["names_error"] = wrong_upload(up.name, "names")
        return
    data = up.getvalue()
    try:
        if len(data) > NAMES_FILE_MAX:
            raise ValueError(f"{up.name}: a names file is at most {NAMES_FILE_MAX // 1000} kB")
        parse_names(data.decode("utf-8"), up.name)
        kept = str(keep_upload(up.name, data, (".yaml", ".yml")))
    except UnicodeDecodeError:
        ss["names_error"] = f"{up.name}: a names file is UTF-8 text"
        return
    except ValueError as err:
        ss["names_error"] = str(err)
        return
    except OSError:
        ss["names_error"] = "The names file could not be stored on the server; try again."
        return
    ss["names"] = [*(n for n in ss.get("names") or [] if n in NAME_KEYWORDS), kept]
    ss.update(names_uploaded=(kept, upload_name(up.name)), names_missing=[])
    ss["_synced"] = {**(ss.get("_synced") or {}), "names": ss["names"]}


def clear_names(ss: MutableMapping) -> None:
    """Remove: the names file goes, the keyword stays."""
    ss["names"] = [n for n in ss.get("names") or [] if n in NAME_KEYWORDS]
    ss.pop("names_uploaded", None)
    ss.update(names_error="", last_changed="names", _synced={**(ss.get("_synced") or {}), "names": ss["names"]})


def domains_block(st, ss: MutableMapping) -> None:
    """Domains: automatic, or a table of names and residue ranges (or, as text, the command line's lines); with
    domains, colour by domain and the domain bar."""
    st.toggle(
        "Automatic domains",
        key="domains__auto",
        on_change=switch_domains,
        args=(ss,),
        persist_state="session",
        help="inferred from residue contacts, chain by chain (copies of a protein alike): not a published "
        "domain assignment",
    )
    st.toggle(
        "CATH domains",
        help="the published CATH chopping of a PDB entry (fetched from PDBe); with none, the note says so",
        **setting(ss, "cath_domains"),
    )
    if not ss.get("domains__auto"):
        text = st.column_config.TextColumn
        ranges = text("Ranges", width=CELL + 10, help="residue ranges and elements, e.g. A:1-80, α3")
        editor(st, ss, "domains", "Domains", {"name": text("Name", width=CELL - 10), "ranges": ranges})
    if ss.get("domains") or ss.get("uniprot_domains") or ss.get("cath_domains"):
        st.toggle(
            "Colour by domain",
            help="elements in eight colour-blind-safe tones, one per domain, instead of the chain colours (a theme "
            "with its own colouring keeps it)",
            **setting(ss, "domain_colours"),
        )
        segmented(
            st,
            "Domain bar",
            CHOICES["domain_bar"],
            setting(ss, "domain_bar"),
            help="a bar under the figure with each chain's domains in residue order · auto: for 2 or more domains · "
            "on: for any · off: never",
        )


CELL, NARROW = 130, 60  # column widths (px): two cells, or an element and two numbers, fit the sidebar
ELEMENT_COLUMNS = {  # domains: names and ranges are typed
    "swap": ("a", "b"),
    "move": ("element",),
    "rename": ("element",),
    "mark": ("element",),
    "domains": (),
}


def element_options(kind: str, rows: list[dict], choices: list[tuple[str, str]]) -> list[str]:
    """What an element column offers: the drawn elements' references, then any other the rows already hold (a
    label typed as text, a residue range), so a cell never holds a value it does not offer."""
    refs = [ref for ref, _ in choices]
    held = (row.get(c) for row in rows for c in ELEMENT_COLUMNS[kind])
    return refs + [h for h in dict.fromkeys(held) if isinstance(h, str) and h and h not in refs]


def editor_frame(kind: str, rows: list[dict]):
    """An editor's rows as the table the data editor shows: its columns in order, move's offsets as numbers (a
    cell that is not one is empty)."""
    import pandas as pd

    frame = pd.DataFrame(rows, columns=list(EDITORS[kind]))
    for column in ("dx", "dy") if kind == "move" else ():
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype(float)
    return frame


def editor(st, ss: MutableMapping, kind: str, label: str, columns: dict) -> None:
    """One editor (in its tab, which names it): its "as text" switch with the editor's help, then a table (rows of
    the text, element cells picked from the drawn elements) or the text itself (labelled `label` for screen
    readers), and under it the exact text it writes. `columns`: each column's config; a str is an element column of
    that label. A table edit is applied by app_edits.apply_table."""
    st.toggle("as text", key=f"{kind}__as_text", help=EDITOR_HELP[kind], persist_state="session")
    if ss.get(f"{kind}__as_text"):
        field = st.text_area if kind == "domains" else st.text_input  # domains: one per line
        field(label, placeholder=EXAMPLE[kind], label_visibility="collapsed", **setting(ss, kind))
    else:
        rows = table_rows(ss, kind)
        choices = element_choices(ss.get("elements") or [])
        options = element_options(kind, rows, choices)
        config = {
            column: st.column_config.SelectboxColumn(
                spec, options=options, format_func=lambda ref: display(ref, choices), width=CELL
            )
            if isinstance(spec, str)
            else spec
            for column, spec in columns.items()
        }
        st.data_editor(
            editor_frame(kind, rows),
            key=f"{kind}__table_{ss.get(f'{kind}__ver', 0)}",
            num_rows="dynamic",
            hide_index=True,
            column_config=config,
            placeholder="default" if kind == "mark" else "",
            on_change=apply_table,
            args=(ss, kind),
        )
    if ss.get(kind):
        st.caption(f"Writes: {shown(ss[kind])}".replace("\n", "  \n"))  # a line each (domains)


UNIPROT = ("off", "auto", "accession")
SITES = ("auto", "all", "none", "types")
CONSERVATION_NOTE = "Then colour by 'conservation'; a ConSurf grades file colours by its own grades."
ALIGNMENT_TYPES = tuple(s.lstrip(".") for s in UPLOADS["alignment"][0])  # the uploader's files (its help lists them)
ALIGNMENT_HELP = (
    "an alignment with the structure's sequence among its rows: colour by 'conservation', and the Sequence view "
    "shows the alignment · a ConSurf grades table (consurf_grades.txt) colours by ConSurf's own grades · files: "
    + ", ".join(f".{t}" for t in ALIGNMENT_TYPES)
)


def conservation_section(st, ss: MutableMapping) -> None:
    """Conservation & UniProt: the alignment (or ConSurf grades table) upload, the alignment's row that is the
    structure's (picked from its rows), then the UniProt annotation: which entry, its domains as panels and the
    sites drawn. A layout file's alignment that is not on this machine is asked for here (and in the top bar)."""
    with st.expander("Conservation & UniProt", key="sec_conservation"):
        st.file_uploader(
            "Multiple sequence alignment or ConSurf grades table",
            type=upload_types("alignment"),
            key="msa_upload",
            on_change=store_alignment,
            args=(ss,),
            max_upload_size=ALIGNMENT_MAX // MB,
            help=ALIGNMENT_HELP,
        )
        if ss.get("msa_error"):
            st.error(shown(ss["msa_error"]))
        missing = missing_alignment(ss)
        if missing:
            st.warning(
                f"The layout file's alignment {missing} is not on this machine; upload it under Conservation & UniProt."
            )
        st.caption(CONSERVATION_NOTE)
        rows = ss.get("msa_rows") or []
        if ss.get("msa_path") and rows:  # an alignment (a ConSurf table has no rows to pick)
            st.selectbox(
                "Reference row",
                listed(["", *rows], ss.get("msa_reference")),
                format_func=lambda row: capped(row) if row else "best match",
                help="the alignment's row that is the structure's sequence · best match: the row most like it",
                **setting(ss, "msa_reference"),
            )
        segmented(
            st,
            "UniProt",
            UNIPROT,
            part(ss, "uniprot__mode", "uniprot"),
            help="off · auto: each chain's entry, found through SIFTS (or the AlphaFold model's own) · accession: "
            "this entry for every chain",
        )
        if ss.get("uniprot__mode") == "accession":
            st.text_input(
                "UniProt accession",
                placeholder="P0CG48",
                label_visibility="collapsed",
                **part(ss, "uniprot__acc", "uniprot"),
            )
        st.toggle(
            "UniProt domains as panels",
            help="the entry's domains, each a panel with its name (in place of the domains under Edit layout)",
            **setting(ss, "uniprot_domains"),
        )
        segmented(
            st,
            "UniProt sites",
            SITES,
            part(ss, "uniprot_sites__mode", "uniprot_sites"),
            help="the entry's sites marked on the figure · auto: active and binding sites · types: the kinds picked",
        )
        if ss.get("uniprot_sites__mode") == "types":
            st.multiselect(
                "Site types",
                list(SITE_GROUPS),
                format_func=SITE_GROUPS.get,
                placeholder="Pick site types",
                label_visibility="collapsed",
                **part(ss, "uniprot_sites__types", "uniprot_sites"),
            )


def missing_alignment(ss: MutableMapping) -> str:
    """A layout file's alignment that is not on this machine, as shown (hosted: never a server path); '' when there
    is none or an alignment has been uploaded since."""
    from foldmap.app import local_files

    if not ss.get("msa_missing") or ss.get("msa_path"):
        return ""
    return shown(ss["msa_missing"], paths=not local_files())


def row_names(path) -> list[str]:
    """The names of an alignment's rows, read (ValueError for a file that is not an alignment); [] for a ConSurf
    grades table, which has none (the pipeline reads it as it is)."""
    from foldmap.sequence import is_consurf, read_alignment

    return [] if is_consurf(path) else list(read_alignment(path).names)


def alignment_rows(path: str) -> list[str]:
    """The Reference row picker's rows for an alignment named by a layout file (row_names); [] for no file or a
    file that is not an alignment (its error is shown when the figure reads it)."""
    if not path:
        return []
    try:
        return row_names(path)
    except (ValueError, OSError):
        return []


MB = 1024 * 1024  # bytes in Streamlit's upload megabyte (max_upload_size)
ALIGNMENT_MAX = 10 * MB  # the alignment uploader's cap, checked again before the file is kept or read


def store_alignment(ss: MutableMapping) -> None:
    """The alignment uploader's callback: an uploaded alignment, size-capped, kept (fetch.keep_upload) and read (a
    ConSurf grades table is used as it is: the pipeline reads it), with its rows for the Reference row picker (a
    reference the new rows lack goes back to best match). Else msa_error, and no alignment."""
    from foldmap.fetch import ALIGNMENT_SUFFIXES, keep_upload

    ss["last_changed"] = "msa_path"
    up = ss.get("msa_upload")
    if up is None:
        ss.update(msa_path="", msa_rows=[])
        return
    if wrong_upload(up.name, "alignment"):
        ss.update(msa_error=wrong_upload(up.name, "alignment"), msa_path="", msa_rows=[])
        return
    target = None
    try:
        data = up.getvalue()
        if len(data) > ALIGNMENT_MAX:
            raise ValueError(f"{upload_name(up.name)}: an alignment is at most {ALIGNMENT_MAX // MB} MB")
        target = keep_upload(up.name, data, ALIGNMENT_SUFFIXES)
        rows = row_names(target)
    except (ValueError, OSError) as err:
        said = str(err) if isinstance(err, ValueError) else "the alignment could not be stored; try again"
        ss["msa_error"] = said.replace(str(target), upload_name(up.name)) if target else said  # the upload's name
        ss.update(msa_path="", msa_rows=[])
        return
    if ss.get("msa_reference") and ss["msa_reference"] not in rows:
        ss["msa_reference"] = ""
    ss.update(msa_error="", msa_path=str(target), msa_rows=rows, msa_uploaded=(str(target), upload_name(up.name)))


def contact_lines(report: dict | None, most: int = 10) -> list[str]:
    """The contacts check under the figure: the count line, then one line per problem (at most `most`, then how
    many more). Only the first line is the app's own: the others name the file's residues, so the caller shows them
    through app_edits.shown."""
    if report is None:
        return ["No contact beads in this figure to check."]
    problems = report.get("problems") or []
    lines = [f"Contacts checked against the coordinates: {report['verified']} of {report['drawn']} confirmed."]
    lines += problems[:most]
    if len(problems) > most:
        lines.append(f"… and {len(problems) - most} more")
    return lines
