# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The app's sidebar sections: how each settings widget is made (its key, kept while hidden, the setting it names
when the figure fails) and the sections Style, Show and Chains & assembly.

Streamlit comes in as an argument (app.main imports it); no pandas."""

from __future__ import annotations

from collections.abc import MutableMapping

from foldmap.app_edits import SYMMETRIES, compose_into, hex_colour, sync_composites
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
        segmented(st, "Fill", CHOICES["fill"], setting(ss, "fill"))
        segmented(st, "Helix shading", CHOICES["helix_shading"], setting(ss, "helix_shading"))
        segmented(
            st,
            "Helix angle",
            CHOICES["helix_angle"],
            setting(ss, "helix_angle"),
            help="snap: nearest of vertical or horizontal · upright: always vertical · tilted: as in 3D",
        )
        segmented(st, "Loop shape", CHOICES["loops"], setting(ss, "loops"))
        st.selectbox("Loop colour", CHOICES["loop_color"], **setting(ss, "loop_color"))
        st.toggle("Loop arrows, N → C", **setting(ss, "loop_arrows"))
        st.slider("Loop weight", 0.4, 3.0, step=0.1, **setting(ss, "loop_width"))
        with st.container(horizontal=True, key="style_colours"):
            colour(st, "Background", ss, "background")
            colour(st, "Ink", ss, "ink")


HBONDS = {"off": "off", "na": "protein–DNA/RNA", "interchain": "interchain", "chains": "chosen pairs"}


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
        segmented(
            st,
            "Ligands and ions",
            ("auto", "all", "none", "codes"),
            part(ss, "ligands__mode", "ligands"),
            help="auto: bound ligands and ions, without buffer additives · codes: only these",
        )
        if ss.get("ligands__mode") == "codes":
            st.text_input(
                "Ligand codes",
                placeholder="HEM, ZN",
                label_visibility="collapsed",
                **part(ss, "ligands__codes", "ligands"),
            )
        segmented(
            st,
            "Membrane",
            ("off", "auto", "flip"),
            setting(ss, "membrane"),
            help="auto: find the bilayer and draw it; flip: the same, other side on top",
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
    st.toggle("Base letters", **setting(ss, "nucleotide_labels"))
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
    if mode == "modality":
        with st.container(horizontal=True, key="na_colours"):
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
        if ss.get("chains__as_text"):
            st.text_input("Chains", placeholder="A,B", **setting(ss, "chains"))
        else:
            st.multiselect(
                "Chains",
                listed(ss.get("all_chains") or [], ss.get("chains__pick")),
                placeholder="All chains",
                format_func=capped,
                help="draw only these chains (protein and nucleic acid each restricted only if some are named)",
                **part(ss, "chains__pick", "chains"),
            )
        st.toggle("Chains as text", key="chains__as_text", persist_state="session")
        highlight_block(st, ss)
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
    if ss.get("highlight__as_text"):
        st.text_input("Highlight", placeholder="none, asu, protomer or A,B", **setting(ss, "highlight"))
    else:
        segmented(
            st,
            "Highlight",
            HIGHLIGHT,
            part(ss, "highlight__mode", "highlight"),
            help="what stays in colour, the rest in grey · none: everything · asu: the chains as deposited · "
            "protomer: the first protomer · chains: the chains picked",
        )
        if ss.get("highlight__mode") == "chains":
            st.multiselect(
                "Highlighted chains",
                listed(ss.get("protein_chains") or [], ss.get("highlight__chains")),
                placeholder="Pick chains",
                format_func=capped,
                label_visibility="collapsed",
                **part(ss, "highlight__chains", "highlight"),
            )
    st.toggle("Highlight as text", key="highlight__as_text", persist_state="session")


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
