# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Streamlit app: load a structure (drawn at once with the current theme), then fine-tune and download.

Run with `foldmap app` (or `streamlit run src/foldmap/app.py`)."""

from __future__ import annotations

import base64
import copy
import json
import logging
import os
import re
import tempfile
from collections.abc import MutableMapping
from dataclasses import fields, replace
from pathlib import Path

from foldmap import CREDIT, app_jobs
from foldmap.app_controls import (
    MB,
    alignment_rows,
    chains_section,
    conservation_section,
    layout_section,
    missing_alignment,
    prepare,
    setting,
    show_section,
    style_section,
    touched,
)
from foldmap.app_edits import chain_choices, shown, upload_name
from foldmap.app_views import Workspace, download_menu, explorer_view, figure_view, said, sequence_view
from foldmap.cli import _adjustments, _domain_spec
from foldmap.fetch import fetch, source_id
from foldmap.names import KEYWORDS as NAME_KEYWORDS
from foldmap.render import save, save_svg
from foldmap.spec import FigureSpec, SavedDomains, draw_spec, parse_flag, parse_protomers, parse_vector
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
ASSETS = Path(__file__).parent / "assets"  # icon.png, icon.svg (scripts/make_app_icon.py), previews/
PREVIEWS = ASSETS / "previews"


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
    "domains_from_file": "",  # the domains text a layout file gave, while unchanged: ids it names may be stale
    "chains": "",
    "names": [],  # --names: deposited, beta or a YAML file (later ones win)
    # layout-file options the app has no control for: kept as loaded and saved again
    "flip_v": False,
    "dna": True,
    "symmetry_tol": TOLERANCE,
    "protomers": None,
    "up": None,
    "view": None,
}
VIEWS = ["Figure", "Sequence", "3D"]  # the workspace view switch (key "pane")
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
    if isinstance(domains, list) and settings.get("domains_from_file") == settings.get("domains"):
        domains = SavedDomains(domains, quiet=True)  # a layout file's: ids this structure lacks are dropped; the
        # app shows the note (app_jobs.draw_job), so the pipeline prints none
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
        names=list(settings.get("names") or []),
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


def layout_document(path, settings: dict, uploads: dict[str, str] | None = None) -> dict:
    """The app's settings as a layout file document (element references saved as residue ranges). Hosted, the files
    it names (names files, the alignment) can only be the server's copies of uploads: each is named as uploaded
    (`uploads`: kept path -> uploaded name; else its file name), never by its server path."""
    doc = spec_document(spec_of(path, settings), structure_name=Path(path).name)
    if not local_files():

        def own(p: str) -> str:
            return (uploads or {}).get(p) or Path(p).name

        if doc.get("names"):
            doc["names"] = [n if n in NAME_KEYWORDS else own(n) for n in doc["names"]]
        if (doc.get("layout") or {}).get("msa"):
            doc["layout"]["msa"] = own(doc["layout"]["msa"])
    return doc


def uploaded_files(doc: dict) -> list[str]:
    """The files a hosted layout document names (they were uploads): names files, then the alignment."""
    msa = (doc.get("layout") or {}).get("msa")
    return [n for n in doc.get("names") or [] if n not in NAME_KEYWORDS] + ([msa] if msa else [])


def settings_from_document(doc: dict) -> dict:
    """App settings from a layout file document: the theme, its style changes, layout options, edits as text."""
    theme = doc.get("theme") or "publication"
    style = {_KEY_ALIASES.get(k, k): v for k, v in (doc.get("style") or {}).items()}  # earlier key names too
    out = {"theme": theme, **vars(THEMES[theme]), **style}
    lay = doc.get("layout") or {}
    for key, default in LAYOUT_DEFAULTS.items():
        out[key] = copy.deepcopy(default)  # no list shared between sessions
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
    out["assembly"] = out["assembly"] or "auto"  # `assembly: ""` draws as auto: the picker shows it so
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
        out["domains_from_file"] = out["domains"]
    out["move"] = "; ".join(f"{k}={d[0]:g},{d[1]:g}" for k, d in (e.get("move") or {}).items())
    names = [str(n) for n in doc.get("names") or []]
    out["names"] = [n for n in names if n in NAME_KEYWORDS or (local_files() and Path(n).is_file())]
    out["names_missing"] = [n for n in names if n not in out["names"]]  # saved elsewhere; never a server path hosted
    return out


# ---------------------------------------------------------------------------------------------- look
APP_THEME = {  # Okabe-Ito blue on white, and its dark twin; every colour pair is checked in tests/test_app_look.py
    "theme.light": {
        "primaryColor": "#0072B2",
        "backgroundColor": "#FFFFFF",
        "secondaryBackgroundColor": "#F5F7FA",
        "textColor": "#1F2933",
        "borderColor": "#848D97",  # widget boundaries need 3:1 (Q-A18)
        "baseRadius": "8px",
        "showWidgetBorder": True,
        "headingFontWeights": 600,
    },
    "theme.dark": {
        "primaryColor": "#56B4E9",
        "backgroundColor": "#0F1419",
        "secondaryBackgroundColor": "#1A2129",
        "textColor": "#E6EAEE",
        "borderColor": "#6B7682",
        "baseRadius": "8px",
        "showWidgetBorder": True,
        "headingFontWeights": 600,
    },
}
MUTED = {"light": "#5B6670", "dark": "#9AA5B1"}  # captions, upload limits, placeholders (the stylesheet)
CARD_BORDER = {"light": "#E3E7EC", "dark": "#2A333D"}  # the figure card's edge (decorative)
BUTTON_LABEL = {"light": "#FFFFFF", "dark": "#0F1419"}  # on the primary colour: white on #56B4E9 is only 2.31:1
APP_CONFIG = {  # what `foldmap app` passes and .streamlit/config.toml holds (Streamlit Community Cloud reads it)
    "client": {"toolbarMode": "viewer", "showErrorDetails": "none"},  # no Deploy menu; visitors see no traceback
    "server": {"maxUploadSize": 50},  # MB, every uploader (the layout file and alignment ones cap lower)
    "browser": {"gatherUsageStats": False},
    **APP_THEME,
}
SIDEBAR_WIDTH = 440  # px: the sidebar's segmented choices and toggle rows fit on one line (Streamlit's own: 336)


def _toml(value) -> str:
    return str(value).lower() if isinstance(value, bool) else f'"{value}"' if isinstance(value, str) else str(value)


def config_toml() -> str:
    """APP_CONFIG as .streamlit/config.toml: a section per entry, in order."""
    out = [
        '# Generated from foldmap.app.APP_CONFIG: python -c "from foldmap.app import config_toml; '
        "print(config_toml(), end='')\" > .streamlit/config.toml"
    ]
    for section, values in APP_CONFIG.items():
        out += ["", f"[{section}]", *(f"{key} = {_toml(value)}" for key, value in values.items())]
    return "\n".join(out) + "\n"


def run_args(app_path) -> list[str]:
    """The `streamlit run` command of `foldmap app`: APP_CONFIG as flags (the same values as the config file)."""
    import sys

    flags = [
        f"--{section}.{key}={str(value).lower() if isinstance(value, bool) else value}"
        for section, values in APP_CONFIG.items()
        for key, value in values.items()
    ]
    return [sys.executable, "-m", "streamlit", "run", str(app_path), *flags]


SIDE, PAIR = '[data-testid="stSidebarUserContent"]', '[data-testid="stElementContainer"]'  # PAIR: two a row, or
# one a row where the sidebar is too narrow for two (a phone)


def _colours(mode: str) -> str:
    """The stylesheet's colours of one theme: muted text (captions, upload limits, input hints, placeholders:
    Streamlit fades them below 4.5:1), the chosen option of a segmented control in the text colour (the primary
    colour on its tinted background is 4.22:1 in light), a switch's track when off in the border colour (Streamlit's
    grey is 1.5:1) and its knob when on in the background colour (light grey on #56B4E9 is 1.9:1), the figure
    card's edge."""
    theme = APP_THEME[f"theme.{mode}"]
    switch = '[data-testid="stCheckbox"] label:has(input[role="switch"]'
    return f"""
[data-testid="stCaptionContainer"], [data-testid="stFileUploaderDropzoneInstructions"] span,
[data-testid="InputInstructions"],
[data-testid="stApp"] input::placeholder, [data-testid="stApp"] textarea::placeholder {{opacity:1;color:{MUTED[mode]}}}
[data-testid="stButtonGroup"] button[aria-checked="true"] p {{color:{theme["textColor"]}}}
{switch}:not(:checked)) > div:not([data-testid]) {{background:{theme["borderColor"]}}}
{switch}:checked) > div:not([data-testid]) > div {{background:{theme["backgroundColor"]}}}
.st-key-canvas .fm-card {{border:1px solid {CARD_BORDER[mode]}}}
"""


STYLESHEET = f"""
.st-key-start {{max-width:640px;margin:0 auto}}
@media (min-width:640px) {{[data-testid="stMainBlockContainer"] {{padding:4rem 1.5rem 2rem}}}}
.st-key-topbar {{row-gap:0.25rem}}
.st-key-canvas .fm-card {{border-radius:8px;padding:12px;overflow-x:auto}}
.st-key-canvas svg {{display:block;width:100%;height:auto;max-height:80vh}}
[data-testid="stMainBlockContainer"]:has([data-testid="stAlert"]) .st-key-canvas svg {{max-height:calc(80vh - 5rem)}}
.st-key-seq_options .st-key-seq_full {{padding-top:2.25rem}}
.st-key-explorer iframe {{height:min(1100px,85vh)!important}}
{SIDE} {{padding-bottom:0.5rem}}
{SIDE} [data-testid="stVerticalBlock"] {{gap:0.375rem}}
{SIDE} [data-testid="stExpanderDetails"] {{padding-top:0.5rem;padding-bottom:0.75rem}}
{SIDE} [data-testid="stExpanderDetails"] [data-testid="stMarkdownContainer"] p {{font-size:14px}}
{SIDE} [data-testid="stButtonGroup"] button {{padding-left:0.5rem;padding-right:0.5rem;min-width:0}}
{SIDE} [data-testid="stFileUploaderDropzone"] {{padding:0.375rem 0.5rem;min-height:0;
  flex-direction:row;align-items:center}}
{SIDE} [data-testid="stFileUploaderDropzone"] button {{min-height:2rem}}
{SIDE} [data-testid="stFileUploaderDropzoneInstructions"] {{height:auto}}
.st-key-style_shapes {PAIR}, .st-key-style_angles {PAIR},
.st-key-layout_shape {PAIR}, .st-key-show_choices {PAIR} {{flex:1 1 9.5rem}}
.st-key-show_choices .st-key-ligands__mode {{flex-basis:11.5rem}}
.st-key-show_choices .st-key-membrane {{flex-basis:8.5rem}}
.st-key-reset_theme button {{min-height:0}}
.st-key-credit p {{font-size:12px}}
.st-key-phone_hint {{display:none}}
@media (max-width:640px) {{.st-key-phone_hint {{display:flex}}}}
{_colours("light")}"""
DARK_EXTRA = f"""{_colours("dark")}
.st-key-load button p, .st-key-redraw button p {{color:{BUTTON_LABEL["dark"]}}}
"""


MARK = (  # the icon (scripts/make_app_icon.py) as the Start screen's mark
    '<img src="data:image/png;base64,{}" width="44" height="44" alt="">'.format(
        base64.b64encode((ASSETS / "icon.png").read_bytes()).decode()
    )
)
PHONE_HINT = "Settings: open the panel with » at the top left."  # Streamlit's button shows »


def page_css(kind: str | None) -> str:
    """The stylesheet for the theme in use (st.context.theme.type): DARK_EXTRA added in the dark theme."""
    return STYLESHEET + (DARK_EXTRA if kind == "dark" else "")


def trouble(err: Exception) -> str:
    """What to tell the visitor: a ValueError says what is wrong with the input; anything else means the structure
    file itself is gone or unreadable (an upload is removed after a day), and server paths stay out of the message."""
    if isinstance(err, ValueError):
        return str(err)
    return "this upload is no longer available (or cannot be read); load it again"


# ---------------------------------------------------------------------------------------------- session
STRUCTURE_STYLE = ("mark", "highlight", "hbonds", "backbone_colors")  # Style fields that name this structure's parts
NAMES_AND_DOMAINS = (  # what the names and domains controls hold besides the settings (a layout file replaces it)
    "names__kw",
    "names_upload",
    "names_uploaded",
    "names_error",
    "domains__auto",
    "domains__manual",
)
STRUCTURE_WIDGETS = ("seq_chain", *NAMES_AND_DOMAINS)  # widget keys of one structure
_LOG = logging.getLogger("foldmap.app")


def settings_of(ss: MutableMapping) -> dict:
    """The figure settings in the session: every Style field, every LAYOUT_DEFAULTS key and the theme."""
    return {k: ss[k] for k in [*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme"] if k in ss}


def current_blob(ss: MutableMapping) -> str:
    """The settings as JSON: the key of the figure caches."""
    return json.dumps(settings_of(ss), sort_keys=True, default=str)


def apply_theme(ss: MutableMapping, name: str | None = None) -> None:
    """Make `name` (default: the session's theme) the theme, with every style value of its own except those that
    name this structure's parts (STRUCTURE_STYLE: marks, highlight, H-bond pairs, per-chain colours), which stay
    once set; layout settings not yet in the session get their defaults."""
    name = name or ss["theme"]
    ss["theme"] = name
    for key, value in vars(THEMES[name]).items():
        if key not in STRUCTURE_STYLE or key not in ss:
            ss[key] = copy.deepcopy(value)
    for key, value in LAYOUT_DEFAULTS.items():
        ss.setdefault(key, copy.deepcopy(value))


# ---------------------------------------------------------------------------------------------- sidebar
THEMES_NOTE = "Shown on ubiquitin, the reference structure (alphafold on p53). Pick one to apply it to your structure."


def theme_changed(ss: MutableMapping) -> None:
    apply_theme(ss)
    touched(ss, "theme")


def pick_theme(ss: MutableMapping, name: str) -> None:
    """A thumbnail's button in the theme dialog: that theme, with its values."""
    apply_theme(ss, name)
    ss["last_changed"] = "theme"


def reset_theme(ss: MutableMapping) -> None:
    """Back to the theme's own values (every style value changed since it was picked); what names this structure's
    parts (STRUCTURE_STYLE) and the layout settings stay."""
    apply_theme(ss, ss["theme"])
    ss["last_changed"] = "theme"


def browse_themes(st, ss: MutableMapping) -> None:
    """The theme dialog: a thumbnail per theme (static files drawn on a reference structure, so opening it draws
    nothing) above a button that applies it. A click reruns the dialog alone, which then reruns the whole app
    without it: the dialog closes on the visitor's structure in the new theme."""

    @st.dialog("Themes", width="large")
    def themes() -> None:
        st.caption(THEMES_NOTE)
        with st.container(horizontal=True, gap="xsmall", key="theme_grid"):  # two tiles across a phone
            for theme in PREVIEW_THEMES:
                with st.container(width=150, border=True, gap="xsmall", key=f"tile_{theme}"):
                    st.image(preview_image(theme), alt=f"The {theme} theme on a reference structure")
                    if st.button(theme, key=f"pick_{theme}", on_click=pick_theme, args=(ss, theme), width="stretch"):
                        st.rerun()
            with st.container(width=150, border=True, key="theme_note"):
                st.caption("conservation · needs an alignment")

    themes()


def essentials(st, ss: MutableMapping) -> None:
    """The top of the sidebar, outside any section: theme, colouring, labels, text size, title, reset."""
    from foldmap.app_edits import palette_used

    with st.container(gap="xsmall", key="essentials"):  # tight: the essentials and the sections fit one screen
        with st.container(horizontal=True, vertical_alignment="bottom", key="theme_row"):
            st.selectbox(
                "Theme", list(THEMES), key="theme", on_change=theme_changed, args=(ss,), persist_state="session"
            )
            if st.button("Browse themes", key="theme_browse", type="tertiary", icon=":material/palette:"):
                browse_themes(st, ss)
        companions = [  # beside Colour by: the palette where it colours something, the N → C ramp of `sequence`
            ("palette", "Palette", palette_used(ss["color_by"], ss["loop_color"])),
            ("sequence_map", "N → C colours", ss["color_by"] == "sequence"),
        ]
        beside = [(key, label) for key, label, used in companions if used]
        cols = st.columns(1 + len(beside))  # side by side
        cols[0].selectbox("Colour by", CHOICES["color_by"], **setting(ss, "color_by"))
        for col, (key, label) in zip(cols[1:], beside):
            col.selectbox(label, CHOICES[key], **setting(ss, key))
        with st.container(horizontal=True, gap="xsmall", key="label_toggles"):
            st.toggle("Labels", **setting(ss, "labels"))
            st.toggle("Residue numbers", **setting(ss, "residue_numbers"))
            st.toggle("Legend", **setting(ss, "legend"))
        with st.container(horizontal=True, key="text_row"):
            st.slider("Text size", 0.6, 2.0, step=0.1, **setting(ss, "font_scale"))
            st.text_input("Title", placeholder="No title", **setting(ss, "title"))
        st.button(
            "Reset to theme",
            key="reset_theme",
            type="tertiary",
            icon=":material/restart_alt:",
            on_click=reset_theme,
            args=(ss,),
            help="every style value back to the theme's own (marks, highlight, H-bond pairs and the layout stay)",
        )


def reset_structure(ss: MutableMapping) -> None:
    """Forget everything that belonged to the last structure: layout settings, edits, marks, the alignment, notes,
    editor tables and per-structure widgets. The theme and its other style values carry over."""
    for key, value in LAYOUT_DEFAULTS.items():
        ss[key] = copy.deepcopy(value)
    look = THEMES[ss.get("theme", "publication")]
    for key in STRUCTURE_STYLE:
        ss[key] = getattr(look, key)
    for key in ("msa_path", "msa_error", "msa_reference", "msa_missing", "layout_error"):
        ss[key] = ""
    ss["names_missing"], ss["msa_rows"] = [], []
    ss["_synced"] = {}
    from foldmap.app_edits import EDITORS

    for kind in EDITORS:  # a new table key: no edit of the last structure's table is applied to this one
        ss[f"{kind}__ver"] = ss.get(f"{kind}__ver", 0) + 1
        ss.pop(f"{kind}__blank", None)  # nor its rows added but not filled in
    for key in (
        *STRUCTURE_WIDGETS,
        "drawn_blob",
        "slow_noted",
        "slow_seconds",
        "last_good",
        "last_changed",
        "figure_failed",
        "seq_drawn",
        "page_drawn",
        "msa_uploaded",
    ):
        ss.pop(key, None)
    ss["auto_redraw"] = True
    ss["pane"] = "Figure"


def apply_document(ss: MutableMapping, settings: dict) -> None:
    """A layout file's settings (settings_from_document) in place of the session's: no control is blamed for what
    follows, and the rows waiting in the editors' tables and the names and domains controls' state belong to the
    settings it replaced (those controls are derived from its settings again)."""
    from foldmap.app_edits import EDITORS

    ss.pop("last_changed", None)  # many settings at once: no one control to blame
    for key, value in settings.items():
        ss[key] = value
    for kind in EDITORS:
        ss.pop(f"{kind}__blank", None)
    for key in NAMES_AND_DOMAINS:
        ss.pop(key, None)
    ss["_synced"] = {k: v for k, v in (ss.get("_synced") or {}).items() if k not in ("names", "domains")}


LAYOUT_FILE_MAX = MB  # the layout file uploader's cap (a layout file is a few kB), checked again before it is read


def apply_layout_file(ss: MutableMapping) -> None:
    """The layout file uploader's callback: the file, size-capped, read (layoutfile.load: YAML, checked) and its
    settings in place of the session's (apply_document), with the rows of the alignment it names; else
    layout_error, and nothing changes."""
    from foldmap import layoutfile

    up = ss.get("layout_upload")
    if up is None:
        return
    data = up.getvalue()
    if len(data) > LAYOUT_FILE_MAX:
        ss["layout_error"] = f"{upload_name(up.name)}: a layout file is at most {LAYOUT_FILE_MAX // MB} MB"
        return
    try:
        with tempfile.TemporaryDirectory() as tmp:
            f = Path(tmp) / "layout.yaml"
            f.write_bytes(data)
            doc = layoutfile.load(f)
        settings = settings_from_document(doc)
    except (ValueError, OSError) as err:
        ss["layout_error"] = str(err)
        return
    ss["layout_error"] = ""
    apply_document(ss, settings)
    ss["msa_rows"] = alignment_rows(ss.get("msa_path", ""))


def layout_file_section(st, ss: MutableMapping) -> None:
    """Layout file: load one (its settings replace the session's); saving one is in the Download popover."""
    with st.expander("Layout file", key="sec_layoutfile"):
        st.file_uploader(
            "Load a layout file",
            type=["yaml", "yml"],
            key="layout_upload",
            on_change=apply_layout_file,
            args=(ss,),
            max_upload_size=LAYOUT_FILE_MAX // MB,
            help="a layout file saved from the Download menu (or the command line's --save-layout): the theme, "
            "style, layout options and edits it holds replace the current ones",
        )
        if ss.get("layout_error"):
            st.error(shown(ss["layout_error"]))


def picked_row(table: dict | None, choices: list[tuple[str, str]]) -> str | None:
    """The id of the structure picked in the UniProt table (its selection state), else None."""
    from foldmap.app_edits import row_source

    rows = ((table or {}).get("selection") or {}).get("rows") or []
    return row_source(list(rows), choices)


# ---------------------------------------------------------------------------------------------- UI
def main() -> None:
    import streamlit as st

    from foldmap import uniprot

    st.set_page_config(
        page_title="foldmap",
        page_icon=str(ASSETS / "icon.png"),
        layout="wide",
        initial_sidebar_state=SIDEBAR_WIDTH,  # "auto" (collapsed on a phone), this wide
    )
    st.html(f"<style>{page_css(st.context.theme.type)}</style>")  # once per run, both screens
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
    def sequence_svg(path: str, blob: str, view_json: str) -> bytes:
        return app_jobs.run_job(app_jobs.sequence_job, path, blob, view_json, "svg", limits=lim)

    @st.cache_data(show_spinner=False, max_entries=4)
    def page(path: str, blob: str) -> str:
        return app_jobs.run_job(app_jobs.page_job, path, blob, limits=lim)

    @st.cache_data(show_spinner=False, max_entries=16)
    def contacts(path: str, key: str, _blob: str) -> dict | None:  # key: app_jobs.contacts_key; _blob not hashed
        return app_jobs.run_job(app_jobs.contacts_job, path, _blob, limits=lim)

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
            ss["load_error"] = shown(err, paths=not local_files())
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
            ss["load_error"] = said(message, path=path, name=name, paths=not local_files())
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
            ss["load_error"] = said(app_jobs.explain(err, name, lim)[0], path=path, name=name, paths=not local_files())
            return
        local_path = local_files() and bool(text) and Path(text).expanduser().is_file()
        ss.update(
            path=str(path),
            name=name,
            info=result["summary"],
            elements=result["elements"],
            protein_chains=result["chains"],
            na_chains=result["na_chains"],
            all_chains=chain_choices(size["chain_ids"], [*result["chains"], *result["na_chains"]]),
            protein_chains_total={"assembly": "auto", "count": len(result["chains"])},  # reset_structure: auto
            assemblies=size["assemblies"],
            share_id="" if file or local_path else source_id(text) or "",
            stage="edit",
            load_error="",
            uniprot_for="",
        )
        if result["seconds"] > lim.slow:
            ss.update(auto_redraw=False, slow_noted=str(path), slow_seconds=result["seconds"])

    def load_example() -> None:
        if ss.get("example"):  # the pill stays selected: Start is not shown again until a new structure
            load(ss["example"])

    def pick_row() -> None:
        picked = picked_row(ss.get("uniprot_table"), choices(ss["uniprot_for"]))
        if picked:
            load(picked)

    def new_structure() -> None:
        ss["stage"], ss["share_id"], ss["source"] = "load", "", ""

    def redraw() -> None:
        ss["drawn_blob"] = current_blob(ss)

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

    if ss["stage"] == "load":  # the look's widgets keep their values while Start shows (persist_state)
        with st.container(key="start"):
            with st.container(horizontal=True, vertical_alignment="center", gap="small", key="wordmark"):
                st.html(MARK, width="content")  # the mark beside the name: decorative (no image toolbar to tab to)
                st.header("foldmap", anchor=False, width="content")
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
            with st.container(key="credit"):
                st.caption(CREDIT)
        return

    # ---- workspace: top bar, the active view, then the sidebar (a view may switch auto-redraw off first)
    prepare(ss)
    with st.container(key="phone_hint"):  # shown under 640 px only (the stylesheet), where the sidebar starts closed
        st.caption(PHONE_HINT)
    top = st.container(key="topbar")
    with top:
        status = st.empty()
    ws = Workspace(st, ss, ss["path"], ss["name"], current_blob(ss), lim, render, sequence_svg, page, contacts, status)

    def top_bar() -> None:
        """Name and summary; on the right the view switch, Download and New structure (Redraw and the auto-redraw
        switch too, after a slow figure)."""
        with top:
            with st.container(horizontal=True, vertical_alignment="center"):
                st.markdown(f"**{shown(ws.name, limit=80)}**", width="content")
                with st.container(key="summary"):
                    st.caption(ss.get("info", ""))
                if not ss["auto_redraw"]:
                    st.button("Redraw", key="redraw", type="primary", on_click=redraw)
                    st.toggle("Auto-redraw", key="auto_redraw", persist_state="session")
                st.segmented_control(
                    "View",
                    VIEWS,
                    key="pane",
                    required=True,
                    label_visibility="collapsed",
                    persist_state="session",
                )
                with st.container(horizontal=True, width="content"):  # wraps as one: never a lone button
                    with st.popover("Download", key="downloads"):
                        download_menu(ws)
                    st.button("New structure", key="new_structure", type="tertiary", on_click=new_structure)
            if not ss["auto_redraw"]:
                st.caption(
                    f"Auto-redraw is off for this structure: drawing took {ss.get('slow_seconds', 0.0):.1f} s. "
                    "Press Redraw after your changes."
                )
            if missing_alignment(ss):
                st.caption(f"⚠ Alignment {missing_alignment(ss)} is not on this machine (Conservation & UniProt).")

    ss.setdefault("auto_redraw", True)
    if ss.get("pane") not in VIEWS:
        ss["pane"] = "Figure"
    {"Figure": figure_view, "Sequence": sequence_view, "3D": explorer_view}[ss["pane"]](ws)  # only this view
    top_bar()

    # ---- sidebar: the essentials, then six sections, then the credit
    with st.sidebar:
        essentials(st, ss)
        with st.container(gap="xxsmall", key="sections"):
            style_section(st, ss)
            show_section(st, ss)
            chains_section(st, ss)
            layout_section(st, ss)
            conservation_section(st, ss)
            layout_file_section(st, ss)
        with st.container(key="credit"):
            st.caption(CREDIT)


def run() -> None:
    """Entry point for `foldmap app`."""
    import subprocess

    raise SystemExit(subprocess.call(run_args(__file__)))


if __name__ == "__main__":
    main()
