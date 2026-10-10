"""Streamlit app: load a structure, pick a style from previews of it, then fine-tune and download.

Run with `foldmap app` (or `streamlit run src/foldmap/app.py`)."""

from __future__ import annotations

import json
import os
import re
import tempfile
from dataclasses import fields, replace
from pathlib import Path

from foldmap import CREDIT
from foldmap.cli import _adjustments, _domain_spec, make_layout
from foldmap.fetch import ALIGNMENT_SUFFIXES, fetch, keep_upload
from foldmap.render import save, save_svg
from foldmap.spec import FigureSpec, draw_spec
from foldmap.spec import layout_document as spec_document
from foldmap.style import CHOICES, THEMES, Style
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
    "rotate": 0,
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


EXAMPLES = {"1LMB": "λ repressor on DNA", "5NKT": "FimA", "8UTF": "measles F trimer", "1UBQ": "ubiquitin"}


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


def sequence_files(path, settings: dict, view: dict) -> dict[str, bytes]:
    """The sequence figure (with the alignment, if one is loaded) as SVG, PNG and PDF. A ConSurf grades table
    colours the topology only: the sequence is drawn without an alignment."""
    from foldmap.seqplot import draw_sequence
    from foldmap.sequence import is_consurf, read_alignment

    msa = settings.get("msa_path")
    fig = draw_sequence(
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
    out = {"theme": theme, **vars(THEMES[theme]), **(doc.get("style") or {})}
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
    out["uniprot"] = lay.get("uniprot") or "off"
    chains = lay.get("chains")
    out["chains"] = ",".join(chains) if isinstance(chains, list) else (chains or "")
    msa = lay.get("msa")
    out["msa_path"] = msa if msa and local_files() and Path(msa).is_file() else ""  # no server paths when hosted
    out["msa_missing"] = msa if msa and not out["msa_path"] else ""  # saved on another machine
    out["title"] = out["title"] or ""
    out["rotate"] = int(round(float(out["rotate"])))
    e = doc.get("edits") or {}
    out["rename"] = "; ".join(f"{k}={v}" for k, v in (e.get("rename") or {}).items())
    out["swap"] = "; ".join(f"{a},{b}" for a, b in e.get("swap") or [])
    out["domains"] = "\n".join(f"{name}={','.join(refs)}" for name, refs in (doc.get("domains") or {}).items())
    out["move"] = "; ".join(f"{k}={d[0]:g},{d[1]:g}" for k, d in (e.get("move") or {}).items())
    return out


def summary(path) -> str:
    from foldmap.symmetry import detect_symmetry

    layout, sses, bb = make_layout(path)
    chains = sorted({l.chain for l in bb.labels})
    helices, strands = sum(s.kind == "H" for s in sses), sum(s.kind == "E" for s in sses)
    try:
        sym = detect_symmetry(bb, sses)
    except ValueError:
        sym = None
    extra = f" · symmetry {sym.label}" if sym else ""
    return f"{len(chains)} chain(s) · {len(bb)} residues · {helices} helices · {strands} strands{extra}"


# ---------------------------------------------------------------------------------------------- UI
def main() -> None:
    import streamlit as st
    import streamlit.components.v1 as components

    st.set_page_config(page_title="Foldmap", page_icon="🧬", layout="wide")
    st.caption(CREDIT)
    ss = st.session_state
    ss.setdefault("stage", "load")

    @st.cache_data(show_spinner=False, max_entries=32)
    def render(path: str, settings_json: str) -> tuple[str, dict]:
        settings = json.loads(settings_json)
        fig, loops = figure(path, settings)
        warn = [f"{l.a_id}>{l.b_id}" for l in loops if l.fallback]
        return save_svg(fig), {"fallback": warn}

    @st.cache_data(show_spinner=False, max_entries=16)
    def files(path: str, settings_json: str) -> dict:
        return figure_files(path, json.loads(settings_json))

    @st.cache_data(show_spinner=False, max_entries=16)
    def sequence(path: str, settings_json: str, view_json: str) -> dict:
        return sequence_files(path, json.loads(settings_json), json.loads(view_json))

    def store_alignment() -> None:
        up = ss.get("msa_upload")
        if up is None:
            ss["msa_path"] = ""
            return
        target = keep_upload(up.name, up.getvalue(), ALIGNMENT_SUFFIXES)
        try:
            from foldmap.sequence import is_consurf, read_alignment

            if not is_consurf(target):  # a ConSurf grades table is used as it is (make_layout reads it)
                read_alignment(target)
        except ValueError as err:
            ss["msa_error"], ss["msa_path"] = str(err), ""
            return
        ss["msa_error"], ss["msa_path"] = "", str(target)

    @st.cache_data(show_spinner=False, max_entries=8)
    def interactive(path: str, settings_json: str) -> str:
        from foldmap.interactive import build_page

        settings = json.loads(settings_json)
        return build_page(path, **spec_of(path, settings).page_kwargs())

    def apply_theme(name: str | None = None) -> None:
        name = name or ss["theme"]
        ss["theme"] = name
        for key, value in vars(THEMES[name]).items():
            ss[key] = value
        for key, value in LAYOUT_DEFAULTS.items():
            ss.setdefault(key, value)

    def load(example: str | None = None) -> None:
        up = ss.get("upload")
        try:
            upload = (up.name, up.getvalue()) if up and not example else None
            path = fetch(example or ss.get("source", ""), upload, local=local_files())
            ss["info"] = summary(path)
        except (ValueError, FileNotFoundError, OSError, RuntimeError) as err:  # RuntimeError: gemmi cannot read it
            ss["load_error"] = str(err)
            return
        ss["path"], ss["stage"], ss["load_error"] = str(path), "style", ""
        for key in HIDDEN_LAYOUT:  # the last structure's are not this one's
            ss[key] = LAYOUT_DEFAULTS[key]
        for key in ("msa_path", "msa_error", "msa_reference", "msa_missing"):  # the last structure's alignment
            ss[key] = ""

    def pick(theme: str) -> None:
        apply_theme(theme)
        ss["stage"] = "edit"

    def restart() -> None:
        ss["stage"] = "load"

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
        except (ValueError, OSError) as err:
            ss["layout_error"] = str(err)
            return
        ss["layout_error"] = ""
        for key, value in settings_from_document(doc).items():
            ss[key] = value

    st.markdown("#### Foldmap · protein topology diagrams")

    if ss["stage"] == "load":
        paths = local_files()
        st.write(
            "Start with a structure: a PDB ID, a UniProt accession (AlphaFold model), "
            + ("a file path, or an upload." if paths else "or an upload.")
        )
        left, right = st.columns([3, 2])
        with left:
            label = "PDB ID, UniProt accession or file path" if paths else "PDB ID or UniProt accession"
            st.text_input(label, key="source", placeholder="e.g. 1LMB or P04637")
            st.button("Load structure", key="load", type="primary", on_click=load)
            st.caption("Examples")
            cols = st.columns(len(EXAMPLES))
            for col, (pdb, what) in zip(cols, EXAMPLES.items()):
                col.button(f"{pdb} · {what}", key=f"ex_{pdb}", on_click=load, args=(pdb,))
        with right:
            st.file_uploader("…or upload a file", type=["cif", "mmcif", "pdb", "ent"], key="upload")
        if ss.get("load_error"):
            st.error(ss["load_error"])
        try:
            choices = structure_choices(ss.get("source", ""))
        except ValueError as err:
            choices = []
            st.caption(f"UniProt: {err}")
        if choices:
            st.write("**Structures of this protein** (UniProt): pick one, or load the AlphaFold model above.")
            for sid, text in choices:
                c1, c2 = st.columns([1, 5])
                c1.button(sid, key=f"pdb_{sid}", on_click=load, args=(sid,))
                c2.caption(text)
        return

    path = ss["path"]
    if ss["stage"] == "style":
        st.write(f"**{Path(path).stem}** · {ss.get('info', '')}")
        st.write(
            "Pick a style to start from: the thumbnails show each theme on ubiquitin (AlphaFold confidence on "
            "p53); your structure is drawn once you choose, and every control stays adjustable afterwards."
        )
        for start in range(0, len(PREVIEW_THEMES), 5):  # rows of five, so each row lines up
            cols = st.columns(5)
            for col, theme in zip(cols, PREVIEW_THEMES[start : start + 5]):
                with col:
                    st.image(preview_image(theme), caption=theme, width="stretch")
                    st.button("Use this style", key=f"pick_{theme}", on_click=pick, args=(theme,))
        st.button("← another structure", key="back_style", on_click=restart)
        return

    # ---- editor
    with st.sidebar:
        st.selectbox("Theme", list(THEMES), key="theme", on_change=apply_theme)
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
                st.error(ss["msa_error"])
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
            st.slider("Rotate (°)", -180, 180, step=15, key="rotate")
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
            st.text_input("Chains", key="chains", help="draw only these chains, e.g. A,B (empty: all)")
            st.text_input(
                "Large assemblies",
                key="focus",
                help="auto: filaments, cages and big assemblies as one subunit plus the neighbouring parts of its "
                "fold; none: draw everything; or chains to show in full, e.g. B",
            )
        st.file_uploader("Load a layout file", type=["yaml", "yml"], key="layout_upload", on_change=apply_layout_file)
        if ss.get("layout_error"):
            st.error(ss["layout_error"])
        if ss.get("msa_missing") and not ss.get("msa_path"):
            st.warning(
                f"The layout file's alignment {ss['msa_missing']} is not on this machine; upload it under Alignment."
            )
        st.button("← another structure", key="back_edit", on_click=restart)

    settings = {k: ss[k] for k in [*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme"] if k in ss}
    blob = json.dumps(settings, sort_keys=True, default=str)
    st.write(f"**{Path(path).stem}** · {ss.get('info', '')}")
    tab_fig, tab_seq, tab_live = st.tabs(["Figure", "Sequence", "Interactive"])
    with tab_fig:
        try:
            svg, notes = render(path, blob)
        except ValueError as err:
            st.error(str(err))
            return
        if notes["fallback"]:
            st.warning("No clear route for: " + ", ".join(notes["fallback"]) + " (drawn as plain curves).")
        page_bg = settings.get("background", "#ffffff")
        st.markdown(
            f'<div style="background:{page_bg};border-radius:6px;padding:8px">{svg}</div>', unsafe_allow_html=True
        )
        out = files(path, blob)
        stem = Path(path).stem
        c1, c2, c3 = st.columns(3)
        c1.download_button("SVG (editable)", out["svg"], f"{stem}.svg", "image/svg+xml", key="dl_svg")
        c2.download_button("PNG (300 dpi)", out["png"], f"{stem}.png", "image/png", key="dl_png")
        c3.download_button("PDF", out["pdf"], f"{stem}.pdf", "application/pdf", key="dl_pdf")
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
        except ValueError as err:
            st.caption(f"Layout file unavailable: {err}")
    with tab_seq:
        c1, c2, c3 = st.columns(3)
        c1.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns")
        c2.selectbox("Element colours", ["figure", "black"], key="seq_ss_colour")
        c3.checkbox("Include unmodelled ends", key="seq_full")
        view = {"columns": ss["seq_columns"], "ss_colour": ss["seq_ss_colour"], "full_sequence": ss["seq_full"]}
        try:
            seq = sequence(path, blob, json.dumps(view, sort_keys=True))
        except ValueError as err:
            st.error(str(err))
        else:
            page_bg = settings.get("background", "#ffffff")
            st.markdown(
                f'<div style="background:{page_bg};border-radius:6px;padding:8px;overflow-x:auto">'
                f"{seq['svg'].decode()}</div>",
                unsafe_allow_html=True,
            )
            stem = Path(path).stem
            d1, d2, d3 = st.columns(3)
            d1.download_button("SVG", seq["svg"], f"{stem}-sequence.svg", "image/svg+xml", key="dl_seq_svg")
            d2.download_button("PNG (300 dpi)", seq["png"], f"{stem}-sequence.png", "image/png", key="dl_seq_png")
            d3.download_button("PDF", seq["pdf"], f"{stem}-sequence.pdf", "application/pdf", key="dl_seq_pdf")
    with tab_live:
        page = interactive(path, blob)
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
