"""Streamlit app: load a structure, pick a style from previews of it, then fine-tune and download.

Run with `foldmap app` (or `streamlit run src/foldmap/app.py`)."""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
import urllib.request
from dataclasses import fields, replace
from pathlib import Path

from foldmap.cli import _adjustments, _domain_spec, _protomers, make_figure_and_loops, make_layout
from foldmap.render import save, save_svg
from foldmap.style import CHOICES, THEMES, Style

STYLE_FIELDS = {f.name: f.type for f in fields(Style)}
PREVIEW_THEMES = ["publication", "shaded", "trace", "rainbow", "richardson", "flexibility", "hydropathy",
                  "goodsell", "minimal"]
LAYOUT_DEFAULTS = {"mode": "projected", "rotate": 0, "symmetry": "auto", "assembly": "auto", "rename": "",
                   "swap": "", "move": "", "title": "", "domains": ""}
CACHE = Path.home() / ".cache" / "foldmap"
EXAMPLES = {"1LMB": "λ repressor on DNA", "5NKT": "FimA", "8UTF": "measles F trimer", "1UBQ": "ubiquitin"}


# ---------------------------------------------------------------------------------------------- core
def fetch(source: str, upload: tuple[str, bytes] | None = None) -> Path:
    """A local structure file from an upload, a path, or a PDB ID (downloaded from RCSB once, then cached)."""
    if upload is not None:
        name, data = upload
        CACHE.mkdir(parents=True, exist_ok=True)
        path = CACHE / f"upload-{hashlib.sha1(data).hexdigest()[:12]}{Path(name).suffix or '.cif'}"
        path.write_bytes(data)
        return path
    text = (source or "").strip()
    if text and Path(text).expanduser().is_file():
        return Path(text).expanduser()
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", text):
        CACHE.mkdir(parents=True, exist_ok=True)
        path = CACHE / f"{text.upper()}.cif"
        if not path.is_file():
            try:
                with urllib.request.urlopen(f"https://files.rcsb.org/download/{text.upper()}.cif", timeout=30) as r:
                    path.write_bytes(r.read())
            except OSError as err:
                raise ValueError(f"could not download {text.upper()} from the PDB ({err})") from None
        return path
    raise ValueError("Enter a 4-character PDB ID (e.g. 1LMB), a path to a .cif or .pdb file, or upload a file.")


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
    return {"symmetry": settings.get("symmetry", "auto") or "auto", "assembly": settings.get("assembly", "auto") or "auto",
            "rotate": float(settings.get("rotate", 0)), "domains": domains, **_adjustments(Args)}


def figure(path, settings: dict):
    return make_figure_and_loops(path, settings.get("mode", "projected"), title=settings.get("title") or None,
                                 look=look_of(settings), **layout_options(settings))


def figure_files(path, settings: dict) -> dict[str, bytes]:
    fig, _ = figure(path, settings)
    out = {"svg": save_svg(fig).encode()}
    with tempfile.TemporaryDirectory() as tmp:
        for ext in ("png", "pdf"):
            out[ext] = save(fig, Path(tmp) / f"figure.{ext}").read_bytes()
    return out


def layout_document(path, settings: dict) -> dict:
    """The app's settings as a layout file document (element references saved as residue ranges)."""
    from foldmap import layoutfile

    look, opts = look_of(settings), layout_options(settings)
    edits = {k: opts.pop(k) for k in ("rename", "swap", "move")}
    layout, sses, _ = make_layout(path, settings.get("mode", "projected"), opts["rotate"], look=look,
                                  symmetry=opts["symmetry"], assembly=opts["assembly"], domains=opts["domains"])
    options = {"mode": settings.get("mode", "projected"), "rotate": opts["rotate"], "symmetry": opts["symmetry"],
               "assembly": opts["assembly"], "title": settings.get("title") or None}
    return layoutfile.document(Path(path).name, settings.get("theme", "publication"), look, options, edits, layout, sses)


def settings_from_document(doc: dict) -> dict:
    """App settings from a layout file document: the theme, its style changes, layout options, edits as text."""
    theme = doc.get("theme") or "publication"
    out = {"theme": theme, **vars(THEMES[theme]), **(doc.get("style") or {})}
    lay = doc.get("layout") or {}
    for key, default in LAYOUT_DEFAULTS.items():
        out[key] = lay.get(key, default) if key in ("mode", "rotate", "symmetry", "assembly", "title") else default
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
    ss = st.session_state
    ss.setdefault("stage", "load")

    @st.cache_data(show_spinner=False, max_entries=64)
    def preview(path: str, theme: str) -> bytes:
        fig, _ = figure(path, {"theme": theme, "legend": False})
        with tempfile.TemporaryDirectory() as tmp:
            fig.set_dpi(60)
            return save(fig, Path(tmp) / "p.png").read_bytes()

    @st.cache_data(show_spinner=False, max_entries=32)
    def render(path: str, settings_json: str) -> tuple[str, dict]:
        settings = json.loads(settings_json)
        fig, loops = figure(path, settings)
        warn = [f"{l.a_id}>{l.b_id}" for l in loops if l.fallback]
        return save_svg(fig), {"fallback": warn}

    @st.cache_data(show_spinner=False, max_entries=16)
    def files(path: str, settings_json: str) -> dict:
        return figure_files(path, json.loads(settings_json))

    @st.cache_data(show_spinner=False, max_entries=8)
    def interactive(path: str, settings_json: str) -> str:
        from foldmap.interactive import build_page

        settings = json.loads(settings_json)
        return build_page(path, settings.get("mode", "projected"), look=look_of(settings),
                          title=settings.get("title") or None, **layout_options(settings))

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
            path = fetch(example or ss.get("source", ""), (up.name, up.getvalue()) if up and not example else None)
            ss["info"] = summary(path)
        except (ValueError, FileNotFoundError, OSError) as err:
            ss["load_error"] = str(err)
            return
        ss["path"], ss["stage"], ss["load_error"] = str(path), "style", ""

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
        st.write("Start with a structure: a PDB ID, a file path, or an upload (mmCIF or PDB).")
        left, right = st.columns([3, 2])
        with left:
            st.text_input("PDB ID or file path", key="source", placeholder="e.g. 1LMB")
            st.button("Load structure", key="load", type="primary", on_click=load)
            st.caption("Examples")
            cols = st.columns(len(EXAMPLES))
            for col, (pdb, what) in zip(cols, EXAMPLES.items()):
                col.button(f"{pdb} · {what}", key=f"ex_{pdb}", on_click=load, args=(pdb,))
        with right:
            st.file_uploader("…or upload a file", type=["cif", "mmcif", "pdb", "ent"], key="upload")
        if ss.get("load_error"):
            st.error(ss["load_error"])
        return

    path = ss["path"]
    if ss["stage"] == "style":
        st.write(f"**{Path(path).stem}** · {ss.get('info', '')}")
        st.write("Pick a style to start from; every control stays adjustable afterwards.")
        cols = st.columns(3)
        for k, theme in enumerate(PREVIEW_THEMES):
            with cols[k % 3]:
                try:
                    st.image(preview(path, theme), caption=theme, width="stretch")
                except Exception as err:  # noqa: BLE001 - a preview that fails should not block the others
                    st.warning(f"{theme}: {err}")
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
            st.text_area("Domains", key="domains", height=80,
                         help="one per line, NAME=REF[,REF...], e.g. ZPN=res:A:331-440; or just 'auto'")
        with st.expander("Symmetry and assembly"):
            st.text_input("Symmetry", key="symmetry", help="auto, off, C2, D3, helical")
            st.text_input("Assembly", key="assembly", help="auto, asu, or an assembly id")
        st.file_uploader("Load a layout file", type=["yaml", "yml"], key="layout_upload", on_change=apply_layout_file)
        if ss.get("layout_error"):
            st.error(ss["layout_error"])
        st.button("← another structure", key="back_edit", on_click=restart)

    settings = {k: ss[k] for k in [*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme"] if k in ss}
    blob = json.dumps(settings, sort_keys=True, default=str)
    st.write(f"**{Path(path).stem}** · {ss.get('info', '')}")
    tab_fig, tab_live = st.tabs(["Figure", "Interactive"])
    with tab_fig:
        try:
            svg, notes = render(path, blob)
        except ValueError as err:
            st.error(str(err))
            return
        if notes["fallback"]:
            st.warning("No clear route for: " + ", ".join(notes["fallback"]) + " (drawn as plain curves).")
        st.markdown(f'<div style="background:#fff;border-radius:6px;padding:8px">{svg}</div>',
                    unsafe_allow_html=True)
        out = files(path, blob)
        stem = Path(path).stem
        c1, c2, c3 = st.columns(3)
        c1.download_button("SVG (editable)", out["svg"], f"{stem}.svg", "image/svg+xml", key="dl_svg")
        c2.download_button("PNG (300 dpi)", out["png"], f"{stem}.png", "image/png", key="dl_png")
        c3.download_button("PDF", out["pdf"], f"{stem}.pdf", "application/pdf", key="dl_pdf")
        try:
            import yaml

            doc = layout_document(path, settings)
            st.download_button("Layout file (YAML): re-render this figure later", yaml.safe_dump(
                doc, sort_keys=False, allow_unicode=True).encode(), f"{stem}-layout.yaml", "application/x-yaml",
                key="dl_layout")
        except ValueError as err:
            st.caption(f"Layout file unavailable: {err}")
    with tab_live:
        page = interactive(path, blob)
        st.download_button("Interactive page (HTML)", page.encode(), f"{stem}-explorer.html", "text/html", key="dl_html")
        components.html(page, height=1100, scrolling=True)


def run() -> None:
    """Entry point for `foldmap app`."""
    import subprocess
    import sys

    raise SystemExit(subprocess.call([sys.executable, "-m", "streamlit", "run", __file__]))


if __name__ == "__main__":
    main()
