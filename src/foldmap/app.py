# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Streamlit app: load a structure (drawn at once with the current theme), then fine-tune and download.

Run with `foldmap app` (or `streamlit run src/foldmap/app.py`)."""

from __future__ import annotations

import copy
import json
import logging
import os
import re
import tempfile
from collections.abc import MutableMapping
from dataclasses import fields, replace
from pathlib import Path

from foldmap import CREDIT
from foldmap.app_edits import shown
from foldmap.cli import _adjustments, _domain_spec
from foldmap.fetch import ALIGNMENT_SUFFIXES, fetch, keep_upload, source_id
from foldmap.render import save, save_svg
from foldmap.spec import FigureSpec, draw_spec, parse_flag, parse_protomers, parse_vector
from foldmap.spec import layout_document as spec_document
from foldmap.style import _KEY_ALIASES, CHOICES, THEMES, Style
from foldmap.symmetry import TOLERANCE

STYLE_FIELDS = {f.name: f.type for f in fields(Style)}
PREVIEW_THEMES = [  # shown in the style picker as static thumbnails (scripts/make_previews.py)
    "publication",
    "journal",
    "minimal",
    "print",
    "shaded",
    "trace",
    "rainbow",
    "richardson",
    "cartoon",
    "presentation",
    "flexibility",
    "hydropathy",
    "alphafold",
    "goodsell",
    "blueprint",
]
PREVIEWS = Path(__file__).parent / "assets" / "previews"


def preview_image(theme: str) -> bytes:
    """The theme's thumbnail, drawn once on a reference structure and shipped with foldmap."""
    return (PREVIEWS / f"{theme}.png").read_bytes()


LAYOUT_DEFAULTS = {
    "mode": "projected",
    "rotate": 0.0,
    "symmetry": "auto",
    "assembly": "auto",
    "rename": "",
    "swap": "",
    "move": "",
    "title": "",
    "membrane": "off",
    "focus": "auto",
    "uniprot": "off",
    "uniprot_domains": False,
    "msa_path": "",
    "msa_reference": "",
    "domains": "",
    "chains": "",
    # layout-file options the app has no control for: kept as loaded and saved again
    "flip_v": False,
    "dna": True,
    "symmetry_tol": TOLERANCE,
    "protomers": None,
    "up": None,
    "view": None,
}
HIDDEN_LAYOUT = ("flip_v", "dna", "symmetry_tol", "protomers", "up", "view")  # in layout files, no app control


def local_files() -> bool:
    """Whether the app may open files by path: yes on your own machine; set FOLDMAP_APP_LOCAL_FILES=0 when the
    app is hosted for others, so visitors can only load PDB IDs, UniProt accessions and uploads."""
    return os.environ.get("FOLDMAP_APP_LOCAL_FILES", "1").strip().lower() not in ("0", "false", "no", "off")


EXAMPLES = {"1LMB": "λ repressor–DNA", "5NKT": "FimA", "8UTF": "measles F", "1UBQ": "ubiquitin"}


# ---------------------------------------------------------------------------------------------- core
def look_of(settings: dict) -> Style:
    base = THEMES[settings.get("theme", "publication")]
    return replace(base, **{k: v for k, v in settings.items() if k in STYLE_FIELDS}).validate()


def layout_options(settings: dict) -> dict:
    class Args:  # the CLI's own parsers for rename / swap / move
        rename = [p for p in re.split(r"[;\n]", settings.get("rename", "")) if p.strip()]
        swap = [p.strip() for p in re.split(r"[;\n]", settings.get("swap", "")) if p.strip()]
        move = [p.strip() for p in re.split(r"[;\n]", settings.get("move", "")) if p.strip()]
        up = view = None

    lines = [p.strip() for p in re.split(r"[;\n]", settings.get("domains", "") or "") if p.strip()]
    domains = "auto" if lines == ["auto"] else _domain_spec(lines) or None
    if settings.get("uniprot_domains") and (settings.get("uniprot") or "off") not in ("off", ""):
        domains = "uniprot"
    return {
        "symmetry": settings.get("symmetry", "auto") or "auto",
        "assembly": settings.get("assembly", "auto") or "auto",
        "rotate": float(settings.get("rotate", 0)),
        "domains": domains,
        "membrane": settings.get("membrane", "off") or "off",
        "focus": settings.get("focus", "auto") or "auto",
        "uniprot": None if (settings.get("uniprot") or "off") in ("off", "") else settings["uniprot"],
        "msa": settings.get("msa_path") or None,
        "msa_reference": settings.get("msa_reference") or None,
        "chains": [c.strip() for c in (settings.get("chains") or "").split(",") if c.strip()] or None,
        **_adjustments(Args),
        **{k: settings[k] for k in HIDDEN_LAYOUT if k in settings},
    }


def spec_of(path, settings: dict) -> FigureSpec:
    """The figure the app's settings describe: the same FigureSpec the command line and the MCP server build."""
    return FigureSpec(
        Path(path),
        theme=settings.get("theme", "publication"),
        look=look_of(settings),
        title=settings.get("title") or None,
        mode=settings.get("mode", "projected"),
        **layout_options(settings),
    )


def structure_choices(source: str, limit: int = 12) -> list[tuple[str, str]]:
    """For a UniProt accession: its AlphaFold model and experimental structures (widest coverage first) as
    (id to load, description) pairs; [] for anything else."""
    from foldmap import uniprot

    if not uniprot.ACCESSION.fullmatch((source or "").strip().upper()):
        return []
    e = uniprot.entry(source)
    out = []
    if e.alphafold:
        out.append((f"AF-{e.alphafold}-F1", f"AlphaFold model · whole chain ({len(e.sequence)} aa)"))
    for s in e.structures[:limit]:
        res = f"{s.resolution:.2f} Å" if s.resolution is not None else ""
        out.append((s.id, f"{s.method} {res} · chains {s.chains}".strip()))
    return out


def sequence_figure(path, settings: dict, view: dict):
    """The sequence figure (with the alignment, if one is loaded). A ConSurf grades table colours the topology
    only: the sequence is drawn without an alignment. view: chain, columns, full_sequence, ss_colour."""
    from foldmap.seqplot import draw_sequence
    from foldmap.sequence import is_consurf, read_alignment

    msa = settings.get("msa_path")
    return draw_sequence(
        path,
        chain=view.get("chain") or None,
        alignment=read_alignment(msa) if msa and not is_consurf(msa) else None,
        reference=settings.get("msa_reference") or None,
        columns=int(view.get("columns", 60)),
        look=look_of(settings),
        full_sequence=bool(view.get("full_sequence", False)),
        ss_colour=view.get("ss_colour", "figure"),
        title=settings.get("title") or None,
        msa=msa or None,
    )


def sequence_files(path, settings: dict, view: dict) -> dict[str, bytes]:
    """The sequence figure (sequence_figure) as SVG, PNG and PDF."""
    fig = sequence_figure(path, settings, view)
    out = {"svg": save_svg(fig).encode()}
    with tempfile.TemporaryDirectory() as tmp:
        for ext in ("png", "pdf"):
            out[ext] = save(fig, Path(tmp) / f"sequence.{ext}").read_bytes()
    return out


def figure(path, settings: dict):
    drawn = draw_spec(spec_of(path, settings))
    return drawn.figure, drawn.loops


def figure_files(path, settings: dict) -> dict[str, bytes]:
    fig, _ = figure(path, settings)
    out = {"svg": save_svg(fig).encode()}
    with tempfile.TemporaryDirectory() as tmp:
        for ext in ("png", "pdf"):
            out[ext] = save(fig, Path(tmp) / f"figure.{ext}").read_bytes()
    return out


def layout_document(path, settings: dict) -> dict:
    """The app's settings as a layout file document (element references saved as residue ranges)."""
    return spec_document(spec_of(path, settings), structure_name=Path(path).name)


def settings_from_document(doc: dict) -> dict:
    """App settings from a layout file document: the theme, its style changes, layout options, edits as text."""
    theme = doc.get("theme") or "publication"
    style = {_KEY_ALIASES.get(k, k): v for k, v in (doc.get("style") or {}).items()}  # earlier key names too
    out = {"theme": theme, **vars(THEMES[theme]), **style}
    lay = doc.get("layout") or {}
    for key, default in LAYOUT_DEFAULTS.items():
        out[key] = default
    for key in (
        "mode",
        "rotate",
        "symmetry",
        "assembly",
        "title",
        "membrane",
        "focus",
        "msa_reference",
        *HIDDEN_LAYOUT,
    ):
        if lay.get(key) is not None:
            out[key] = lay[key]
    for key in ("symmetry", "membrane"):  # unquoted `off` in YAML is False
        if out[key] is False:
            out[key] = "off"
    out["up"], out["view"] = parse_vector(out["up"], "up"), parse_vector(out["view"], "view")
    if isinstance(out["protomers"], str):
        out["protomers"] = parse_protomers(out["protomers"])
    out["flip_v"], out["dna"] = parse_flag(out["flip_v"], "flip_v"), parse_flag(out["dna"], "dna")
    out["uniprot"] = lay.get("uniprot") or "off"
    chains = lay.get("chains")
    out["chains"] = ",".join(chains) if isinstance(chains, list) else (chains or "")
    msa = lay.get("msa")
    out["msa_path"] = msa if msa and local_files() and Path(msa).is_file() else ""  # no server paths when hosted
    out["msa_missing"] = msa if msa and not out["msa_path"] else ""  # saved on another machine
    out["title"] = out["title"] or ""
    out["rotate"] = float(out["rotate"])
    e = doc.get("edits") or {}
    out["rename"] = "; ".join(f"{k}={v}" for k, v in (e.get("rename") or {}).items())
    out["swap"] = "; ".join(f"{a},{b}" for a, b in e.get("swap") or [])
    domains = doc.get("domains") or {}
    if domains == "uniprot":  # the app's UniProt domains switch
        out["uniprot_domains"], domains = True, {}
    if isinstance(domains, str):  # auto
        out["domains"] = domains
    else:
        pairs = domains.items() if isinstance(domains, dict) else domains
        out["domains"] = "\n".join(f"{name}={','.join(refs)}" for name, refs in pairs)
    out["move"] = "; ".join(f"{k}={d[0]:g},{d[1]:g}" for k, d in (e.get("move") or {}).items())
    return out


def trouble(err: Exception) -> str:
    """What to tell the visitor: a ValueError says what is wrong with the input; anything else means the structure
    file itself is gone or unreadable (an upload is removed after a day), and server paths stay out of the message."""
    if isinstance(err, ValueError):
        return str(err)
    return "this upload is no longer available (or cannot be read); load it again"


# ---------------------------------------------------------------------------------------------- session
STRUCTURE_STYLE = ("mark", "highlight", "hbonds", "backbone_colors")  # Style fields that name this structure's parts
STRUCTURE_WIDGETS = ("seq_chain", "names__kw", "domains__auto", "names_upload")  # widget keys of one structure
_LOG = logging.getLogger("foldmap.app")


def settings_of(ss: MutableMapping) -> dict:
    """The figure settings in the session: every Style field, every LAYOUT_DEFAULTS key and the theme."""
    return {k: ss[k] for k in [*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme"] if k in ss}


def current_blob(ss: MutableMapping) -> str:
    """The settings as JSON: the key of the figure caches."""
    return json.dumps(settings_of(ss), sort_keys=True, default=str)


def apply_theme(ss: MutableMapping, name: str | None = None) -> None:
    """Make `name` (default: the session's theme) the theme, with every style value of its own; layout settings
    not yet in the session get their defaults."""
    name = name or ss["theme"]
    ss["theme"] = name
    for key, value in vars(THEMES[name]).items():
        ss[key] = copy.deepcopy(value)
    for key, value in LAYOUT_DEFAULTS.items():
        ss.setdefault(key, copy.deepcopy(value))


# interim: remove when the look widgets use persist_state="session" (Task 11)
def keep_look(ss: MutableMapping) -> None:
    """Keep the theme and its style values while their widgets are not shown (the Start screen): the next structure
    is drawn in the visitor's look. Streamlit forgets a widget's value one run after the widget stops rendering,
    unless the value is set again."""
    for key in ("theme", *STYLE_FIELDS):
        if key in ss:
            ss[key] = ss[key]


def reset_structure(ss: MutableMapping) -> None:
    """Forget everything that belonged to the last structure: layout settings, edits, marks, the alignment, notes,
    editor tables and per-structure widgets. The theme and its other style values carry over."""
    for key, value in LAYOUT_DEFAULTS.items():
        ss[key] = copy.deepcopy(value)
    look = THEMES[ss.get("theme", "publication")]
    for key in STRUCTURE_STYLE:
        ss[key] = getattr(look, key)
    for key in ("msa_path", "msa_error", "msa_reference", "msa_missing", "names_missing", "layout_error"):
        ss[key] = ""
    ss["_synced"] = {}
    from foldmap.app_edits import EDITORS

    for kind in EDITORS:  # a new table key: no edit of the last structure's table is applied to this one
        ss[f"{kind}__ver"] = ss.get(f"{kind}__ver", 0) + 1
    for key in (*STRUCTURE_WIDGETS, "drawn_blob", "slow_noted", "last_good"):
        ss.pop(key, None)
    ss["auto_redraw"] = True
    ss["pane"] = "Figure"


def picked_row(table: dict | None, choices: list[tuple[str, str]]) -> str | None:
    """The id of the structure picked in the UniProt table (its selection state), else None."""
    from foldmap.app_edits import row_source

    rows = ((table or {}).get("selection") or {}).get("rows") or []
    return row_source(list(rows), choices)


# ---------------------------------------------------------------------------------------------- UI
def main() -> None:
    import streamlit as st
    import streamlit.components.v1 as components

    from foldmap import app_jobs, uniprot

    st.set_page_config(page_title="Foldmap", page_icon="🧬", layout="wide")
    ss = st.session_state
    ss.setdefault("stage", "load")
    lim = app_jobs.limits(hosted=not local_files())

    @st.cache_data(show_spinner=False, max_entries=32)
    def render(path: str, blob: str) -> dict:
        return app_jobs.run_job(app_jobs.draw_job, path, blob, lim.slow, limits=lim)

    @st.cache_data(show_spinner=False, max_entries=64)
    def sized(path: str) -> dict:
        return app_jobs.run_job(app_jobs.size_job, path, limits=lim)

    @st.cache_data(show_spinner=False, ttl=86400, max_entries=64)
    def choices(accession: str) -> list[tuple[str, str]]:
        return structure_choices(accession)

    @st.cache_data(show_spinner=False, max_entries=16)
    def sequence(path: str, settings_json: str, view_json: str) -> dict:
        return sequence_files(path, json.loads(settings_json), json.loads(view_json))

    def store_alignment() -> None:
        up = ss.get("msa_upload")
        if up is None:
            ss["msa_path"] = ""
            return
        try:
            target = keep_upload(up.name, up.getvalue(), ALIGNMENT_SUFFIXES)
            from foldmap.sequence import is_consurf, read_alignment

            if not is_consurf(target):  # a ConSurf grades table is used as it is (make_layout reads it)
                read_alignment(target)
        except (ValueError, OSError) as err:
            ss["msa_error"] = (
                str(err) if isinstance(err, ValueError) else "the alignment could not be stored; try again"
            )
            ss["msa_path"] = ""
            return
        ss["msa_error"], ss["msa_path"] = "", str(target)

    @st.cache_data(show_spinner=False, max_entries=8)
    def interactive(path: str, settings_json: str) -> str:
        from foldmap.interactive import build_page

        settings = json.loads(settings_json)
        return build_page(path, **spec_of(path, settings).page_kwargs())

    def load(source: str | None = None, upload: bool = False) -> None:
        """Load a structure and draw its first figure; the workspace opens only when that succeeds. upload: the
        file in the Start uploader (its callback); else `source`, or the text typed on Start."""
        if upload:
            up = ss.get("upload")
            if up is None:  # the upload was removed
                return
            text, file, name = "", (up.name, up.getvalue()), re.sub(r"[^\w.-]+", "_", Path(up.name).stem)[:60]
        else:
            if source is None:  # typed: the last accession's list goes
                ss["uniprot_for"] = ""
            text, file = (source if source is not None else ss.get("source", "")).strip(), None
            if uniprot.ACCESSION.fullmatch(text.upper()):  # its structures are listed to choose from
                ss["uniprot_for"], ss["load_error"] = text.upper(), ""
                return
            name = ""
        if "theme" not in ss:  # the first structure of this visit
            apply_theme(ss, "publication")
        try:
            path = fetch(text, file, local=local_files())
        except ValueError as err:  # not an id, or the download failed
            ss["load_error"] = shown(err)
            return
        except OSError:  # an upload that could not be kept (a full disk): downloads fail as ValueError
            ss["load_error"] = "The file could not be stored on the server; try again."
            return
        name = name or source_id(text) or path.stem
        try:
            size = sized(str(path))
        except (ValueError, OSError, RuntimeError, app_jobs.JobLimit) as err:
            if isinstance(err, (OSError, RuntimeError)) and not isinstance(err, app_jobs.JobFailed):
                message = f"{name} could not be read as a structure (mmCIF or PDB)"  # gemmi cannot read it
            else:
                message = app_jobs.explain(err, name, lim)[0]
            ss["load_error"] = shown(message, path=path, name=name)
            return
        if lim.max_residues and size["residues"] > lim.max_residues:
            ss["load_error"] = shown(app_jobs.size_message(name, size, lim))
            return
        reset_structure(ss)
        settings = settings_of(ss)
        failed = app_jobs.known_failure(str(path), settings)
        try:
            if failed is not None:
                raise failed
            result = render(str(path), current_blob(ss))  # a structure that cannot be drawn never opens
        except Exception as err:  # every failure is a message here (a pipeline bug too: logged for the server)
            if isinstance(err, app_jobs.JobLimit):
                app_jobs.remember_failure(str(path), settings, err)
            elif not isinstance(err, (ValueError, OSError, RuntimeError)):
                _LOG.exception("drawing %s failed", name)
            ss["load_error"] = shown(app_jobs.explain(err, name, lim)[0], path=path, name=name)
            return
        local_path = local_files() and bool(text) and Path(text).expanduser().is_file()
        ss.update(
            path=str(path),
            name=name,
            info=result["summary"],
            elements=result["elements"],
            protein_chains=result["chains"],
            na_chains=result["na_chains"],
            all_chains=list(dict.fromkeys([*result["chains"], *result["na_chains"]])),
            assemblies=size["assemblies"],
            share_id="" if file or local_path else source_id(text) or "",
            stage="edit",
            load_error="",
            uniprot_for="",
        )
        if result["seconds"] > lim.slow:
            ss["auto_redraw"], ss["slow_noted"] = False, str(path)

    def load_example() -> None:
        if ss.get("example"):  # the pill stays selected: Start is not shown again until a new structure
            load(ss["example"])

    def pick_row() -> None:
        picked = picked_row(ss.get("uniprot_table"), choices(ss["uniprot_for"]))
        if picked:
            load(picked)

    def new_structure() -> None:
        ss["stage"], ss["share_id"] = "load", ""

    def apply_layout_file() -> None:
        from foldmap import layoutfile

        up = ss.get("layout_upload")
        if up is None:
            return
        try:
            with tempfile.TemporaryDirectory() as tmp:
                f = Path(tmp) / "layout.yaml"
                f.write_bytes(up.getvalue())
                doc = layoutfile.load(f)
            settings = settings_from_document(doc)
        except (ValueError, OSError) as err:
            ss["layout_error"] = str(err)
            return
        ss["layout_error"] = ""
        for key, value in settings.items():
            ss[key] = value

    if ss["stage"] == "load":  # a share link (?id=…&theme=…) opens its structure once
        qid = st.query_params.get("id")
        if qid and ss.get("qp_tried") != qid:
            ss["qp_tried"] = qid
            theme = st.query_params.get("theme")
            if theme in THEMES:
                apply_theme(ss, theme)
            if source_id(qid):
                load(source_id(qid))
            else:
                ss["load_error"] = "The link's id is not a PDB ID or UniProt accession."
    desired = {"id": ss["share_id"], "theme": ss["theme"]} if ss["stage"] == "edit" and ss.get("share_id") else {}
    if st.query_params.to_dict() != desired:  # uploads and paths never enter the link
        st.query_params.from_dict(desired)

    if ss["stage"] == "load":
        keep_look(ss)
        with st.container(key="start"):
            st.header("foldmap", anchor=False)
            st.caption("Protein topology diagrams from structures")
            with st.form("start", border=False, enter_to_submit=True):
                hint = "PDB ID, UniProt accession, or file path" if local_files() else "PDB ID or UniProt accession"
                st.text_input("Structure", key="source", placeholder=hint, label_visibility="collapsed")
                st.form_submit_button("Draw", key="load", type="primary", on_click=load)
            if ss.get("load_error"):
                st.error(ss["load_error"])
            st.pills(
                "Examples",
                list(EXAMPLES),
                key="example",
                format_func=lambda k: f"{k} · {EXAMPLES[k]}",
                on_change=load_example,
            )
            st.file_uploader(
                "or drop a .cif / .pdb file",
                type=["cif", "mmcif", "pdb", "ent"],
                key="upload",
                on_change=load,
                kwargs={"upload": True},
            )
            if ss.get("uniprot_for"):
                acc = ss["uniprot_for"]
                try:
                    rows = choices(acc)
                except ValueError as err:
                    st.caption(f"UniProt: {shown(err)}")
                else:
                    if rows:
                        st.caption(f"Structures of {shown(acc)} (UniProt): pick one to draw it.")
                        st.dataframe(
                            {"id": [r[0] for r in rows], "what": [r[1] for r in rows]},
                            key="uniprot_table",
                            on_select=pick_row,
                            selection_mode="single-row",
                            hide_index=True,
                        )
                    else:
                        st.caption(f"UniProt lists no structure or AlphaFold model for {shown(acc)}.")
            st.caption(CREDIT)
        return

    path, name = ss["path"], ss["name"]
    st.caption(CREDIT)  # Task 4 moves it to the bottom of the sidebar
    st.markdown("#### Foldmap · protein topology diagrams")

    # ---- editor
    with st.sidebar:
        st.selectbox("Theme", list(THEMES), key="theme", on_change=apply_theme, args=(ss,))
        with st.expander("Colours", expanded=True):
            st.selectbox("Colour by", CHOICES["color_by"], key="color_by")
            st.selectbox("Palette", CHOICES["palette"], key="palette")
            if ss["color_by"] == "sequence":
                st.selectbox("N → C ramp", CHOICES["sequence_map"], key="sequence_map")
            st.selectbox("Fill", CHOICES["fill"], key="fill")
            st.selectbox("Helix shading", CHOICES["helix_shading"], key="helix_shading")
            st.text_input("Mark elements", key="mark", help="e.g. G, A:G=#2ca02c, #3, res:150-159=red")
            st.text_input("Highlight", key="highlight", help="none, asu, protomer, or chains like A,B")
            c_bg, c_ink = st.columns(2)
            c_bg.color_picker("Background", key="background")
            c_ink.color_picker("Ink", key="ink")
        with st.expander("Loops"):
            st.selectbox("Loop shape", CHOICES["loops"], key="loops")
            st.selectbox("Loop colour", CHOICES["loop_color"], key="loop_color")
            st.checkbox("Direction arrows", key="loop_arrows")
            st.slider("Loop weight", 0.4, 3.0, step=0.1, key="loop_width")
        with st.expander("Labels"):
            st.checkbox("Element labels", key="labels")
            st.checkbox("Residue numbers", key="residue_numbers")
            st.checkbox("Legend", key="legend")
            st.checkbox("Base letters (DNA/RNA)", key="nucleotide_labels")
            st.selectbox(
                "Stems (DNA/RNA)",
                CHOICES["na_style"],
                key="na_style",
                help="ladder: rails and base-pair rungs · ribbon: double helix · block: bars · letters",
            )
            st.selectbox(
                "DNA/RNA shown",
                CHOICES["na_extent"],
                key="na_extent",
                help="auto: everything up to 300 nt, else the stems the protein touches",
            )
            st.text_input(
                "DNA/RNA backbone", key="na_backbone", help="modality, chain, ink, or one colour such as #8c564b"
            )
            st.text_input(
                "Backbone colours",
                key="backbone_colors",
                help="e.g. A=#1f77b4, dna=#555555, I=crimson (chain or protein/dna/rna)",
            )
            st.text_input(
                "Hydrogen bonds", key="hbonds", help="off, na (protein with DNA/RNA), interchain, or chains:A+B,C+D"
            )
            st.color_picker("H-bond colour", key="hbond_color")
            st.number_input("Most H-bond lines", min_value=1, step=50, key="hbond_max")
            st.slider("Text size", 0.6, 2.0, step=0.1, key="font_scale")
            st.text_input("Title", key="title")
        with st.expander("Features"):
            st.checkbox("Disulfides", key="disulfides")
            st.checkbox("Glycans", key="glycans")
            st.checkbox("Sheet panels", key="sheet_panels")
            st.text_input(
                "Ligands and ions", key="ligands", help="auto (skips buffer additives), all, none, or codes like HEM,ZN"
            )
            st.selectbox(
                "Membrane",
                ["off", "auto", "flip"],
                key="membrane",
                help="auto: find the bilayer and draw it; flip: the same, other side on top",
            )
        with st.expander("Alignment (conservation)"):
            st.file_uploader(
                "Multiple sequence alignment or ConSurf grades table",
                type=["fasta", "fa", "aln", "sto", "txt"],
                key="msa_upload",
                on_change=store_alignment,
            )
            if ss.get("msa_error"):
                st.error(shown(ss["msa_error"]))
            st.text_input(
                "Reference row", key="msa_reference", help="name of the structure's row (default: best match)"
            )
            st.caption(
                "With an alignment: colour by 'conservation', and the Sequence tab shows the alignment. "
                "A ConSurf grades table (consurf_grades.txt) colours by ConSurf's own grades."
            )
        with st.expander("Sizes and layout"):
            st.selectbox("Mode", ["projected", "stack"], key="mode")
            st.selectbox("Helix angle", CHOICES["helix_angle"], key="helix_angle")
            st.slider("Helix size", 0.5, 2.0, step=0.05, key="helix_scale")
            st.slider("Strand size", 0.5, 2.0, step=0.05, key="strand_scale")
            st.slider("DNA size", 0.5, 2.0, step=0.05, key="dna_scale")
            st.slider("Rotate (°)", -180.0, 180.0, step=0.5, key="rotate")
            st.text_input("Swap", key="swap", help="two elements, e.g. α1,α3 (several: separate with ;)")
            st.text_input("Move", key="move", help="e.g. α2=1,-2 (page units)")
            st.text_input("Rename", key="rename", help="e.g. res:167-182=Gd")
            st.text_area(
                "Domains",
                key="domains",
                height=80,
                help="one per line, NAME=REF[,REF...], e.g. ZPN=res:A:331-440; or just 'auto'",
            )
        with st.expander("UniProt annotation"):
            st.text_input(
                "UniProt", key="uniprot", help="off, auto (each chain's entry via SIFTS / AlphaFold), or an accession"
            )
            st.checkbox("UniProt domains as panels", key="uniprot_domains")
            st.text_input(
                "Sites shown",
                key="uniprot_sites",
                help="auto (active + binding), all, none, or types: active,binding,modified",
            )
        with st.expander("Symmetry and assembly"):
            st.text_input("Symmetry", key="symmetry", help="auto, off, C2, D3, helical")
            st.text_input("Assembly", key="assembly", help="auto, asu, or an assembly id")
            st.text_input(
                "Chains",
                key="chains",
                help="draw only these chains, e.g. A,B (empty: all; protein and nucleic acid each restricted if named)",
            )
            st.text_input(
                "Large assemblies",
                key="focus",
                help="auto: filaments, cages and big assemblies as one subunit plus the neighbouring parts of its "
                "fold; none: draw everything; or chains to show in full, e.g. B",
            )
        st.file_uploader("Load a layout file", type=["yaml", "yml"], key="layout_upload", on_change=apply_layout_file)
        if ss.get("layout_error"):
            st.error(shown(ss["layout_error"]))
        if ss.get("msa_missing") and not ss.get("msa_path"):
            st.warning(
                f"The layout file's alignment {shown(ss['msa_missing'], paths=False)} is not on this machine; "
                "upload it under Alignment."
            )
        st.button("← another structure", key="new_structure", on_click=new_structure)

    settings, blob = settings_of(ss), current_blob(ss)
    st.write(f"**{shown(name, limit=80)}** · {shown(ss.get('info', ''))}")
    stem = name
    tab_fig, tab_seq, tab_live = st.tabs(["Figure", "Sequence", "Interactive"])
    with tab_fig:
        try:
            failed = app_jobs.known_failure(path, settings)
            if failed is not None:
                raise failed
            result = render(path, blob)
        except (ValueError, OSError, RuntimeError, app_jobs.JobLimit) as err:
            if isinstance(err, app_jobs.JobLimit):
                app_jobs.remember_failure(path, settings, err)
            st.error(shown(app_jobs.explain(err, name, lim, ss.get("last_changed"))[0], path=path, name=name))
            return
        if result["fallback"]:
            st.warning(f"No clear route for: {shown(', '.join(result['fallback']))} (drawn as plain curves).")
        if result["over_labels"]:
            st.warning(
                f"Drawn over a label or chain end (the only route found): {shown(', '.join(result['over_labels']))}."
            )
        page_bg = settings.get("background", "#ffffff")
        st.markdown(
            f'<div style="background:{page_bg};border-radius:6px;padding:8px">{result["svg"]}</div>',
            unsafe_allow_html=True,
        )
        eager = tuple(sorted(result["files"].items()))
        c1, c2, c3 = st.columns(3)
        for col, ext, label, mime in (
            (c1, "svg", "SVG (editable)", "image/svg+xml"),
            (c2, "png", "PNG (300 dpi)", "image/png"),
            (c3, "pdf", "PDF", "application/pdf"),
        ):
            try:
                data = app_jobs.figure_file(path, blob, ext, lim, eager)
            except (ValueError, OSError, RuntimeError, app_jobs.JobLimit) as err:
                col.caption(
                    f"{ext.upper()} unavailable: {shown(app_jobs.explain(err, name, lim)[0], path=path, name=name)}"
                )
            else:
                col.download_button(label, data, f"{stem}.{ext}", mime, key=f"dl_{ext}")
        try:
            import yaml

            doc = layout_document(path, settings)
            st.download_button(
                "Layout file (YAML): re-render this figure later",
                yaml.safe_dump(doc, sort_keys=False, allow_unicode=True).encode(),
                f"{stem}-layout.yaml",
                "application/x-yaml",
                key="dl_layout",
            )
        except (ValueError, OSError, RuntimeError) as err:
            st.caption(f"Layout file unavailable: {shown(app_jobs.explain(err, name, lim)[0], path=path, name=name)}")
    with tab_seq:
        c1, c2, c3 = st.columns(3)
        c1.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns")
        c2.selectbox("Element colours", ["figure", "black"], key="seq_ss_colour")
        c3.checkbox("Include unmodelled ends", key="seq_full")
        view = {"columns": ss["seq_columns"], "ss_colour": ss["seq_ss_colour"], "full_sequence": ss["seq_full"]}
        try:
            seq = sequence(path, blob, json.dumps(view, sort_keys=True))
        except (ValueError, OSError, RuntimeError) as err:
            st.error(shown(app_jobs.explain(err, name, lim)[0], path=path, name=name))
        else:
            page_bg = settings.get("background", "#ffffff")
            st.markdown(
                f'<div style="background:{page_bg};border-radius:6px;padding:8px;overflow-x:auto">'
                f"{seq['svg'].decode()}</div>",
                unsafe_allow_html=True,
            )
            d1, d2, d3 = st.columns(3)
            d1.download_button("SVG", seq["svg"], f"{stem}-sequence.svg", "image/svg+xml", key="dl_seq_svg")
            d2.download_button("PNG (300 dpi)", seq["png"], f"{stem}-sequence.png", "image/png", key="dl_seq_png")
            d3.download_button("PDF", seq["pdf"], f"{stem}-sequence.pdf", "application/pdf", key="dl_seq_pdf")
    with tab_live:
        try:
            page = interactive(path, blob)
        except (ValueError, OSError, RuntimeError) as err:
            st.error(shown(app_jobs.explain(err, name, lim)[0], path=path, name=name))
        else:
            st.download_button(
                "Interactive page (HTML)", page.encode(), f"{stem}-explorer.html", "text/html", key="dl_html"
            )
            if hasattr(st, "iframe"):  # Streamlit >= 1.52; the page is foldmap's own output
                st.iframe(page, height=1100)
            else:
                components.html(page, height=1100, scrolling=True)


def run() -> None:
    """Entry point for `foldmap app`."""
    import subprocess
    import sys

    raise SystemExit(subprocess.call([sys.executable, "-m", "streamlit", "run", __file__]))


if __name__ == "__main__":
    main()
