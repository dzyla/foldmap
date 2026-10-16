# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The app's workspace views (Figure, Sequence, 3D) and its Download popover: what one run of the workspace draws
and offers, given the Workspace that app.main builds. The figure card and the message helper (said) live here too.

Streamlit comes in through Workspace.st (app.main imports it); app is imported inside the functions that need it."""

from __future__ import annotations

import html
import json
import logging
import re
from collections.abc import Callable, MutableMapping
from dataclasses import dataclass
from functools import partial
from typing import Any, NamedTuple

from foldmap import app_edits, app_jobs
from foldmap.app_controls import capped
from foldmap.app_edits import shown

_LOG = logging.getLogger("foldmap.app")  # the app's logger: one name for every app message on the server


def figure_card(svg: str, background: str, name: str, kind: str = "Topology diagram", natural: bool = False) -> str:
    """The figure as HTML for the canvas: a card on the figure's own background, the SVG named for screen readers
    ("{kind} of {name}"). A figure more than three times as wide as tall keeps its height and scrolls inside the
    card; natural: the figure keeps its own size (the sequence figure, so its letters stay readable) and scrolls
    inside the card when wider."""
    root = re.search(r"<svg\b[^>]*>", svg)
    tag = root.group() if root else ""
    attrs = {  # CSS px; matplotlib writes pt (4/3 px)
        k: float(v) * (4 / 3 if unit == "pt" else 1)
        for k, v, unit in re.findall(r'\s(width|height)="([\d.]+)(pt|px)?"', tag)
    }
    extra = f' role="img" aria-label="{html.escape(kind, quote=True)} of {html.escape(name, quote=True)}"'
    if natural and attrs.get("height") and attrs.get("width"):
        extra += f' style="width:{attrs["width"]:.0f}px;max-width:none;height:{attrs["height"]:.0f}px;max-height:none"'
    elif attrs.get("height") and attrs.get("width", 0) / attrs["height"] > 3:
        extra += f' style="width:auto;max-width:none;height:min(60vh,{attrs["height"]:.0f}px)"'
    if root:
        svg = svg[: root.start()] + "<svg" + extra + tag[4:] + svg[root.end() :]
    return f'<div class="fm-card" style="background:{html.escape(background, quote=True)}">{svg}</div>'


def said(message: str, **kw) -> str:
    """shown(message, **kw), but a bug message's closing request to report it keeps its link (outside the escaped,
    shortened text, so a long exception text never cuts it away)."""
    from foldmap.app_jobs import REPORT, REPORT_SENTENCE

    if not message.endswith(REPORT_SENTENCE):
        return shown(message, **kw)
    link = REPORT.removeprefix("https://")
    return shown(message[: -len(REPORT_SENTENCE)], **kw) + f" Please report it at [{link}]({REPORT})."


# ---------------------------------------------------------------------------------------------- views
@dataclass
class Workspace:
    """What the views of one workspace run share: Streamlit, the session, the structure, the settings JSON (the
    cache key), the job limits, the app's caches (render, sequence_svg, page, contacts) and where the spinner
    shows."""

    st: Any
    ss: MutableMapping
    path: str
    name: str
    blob: str
    lim: app_jobs.Limits
    render: Callable[[str, str], dict]
    sequence_svg: Callable[[str, str, str], bytes]
    page: Callable[[str, str], str]
    contacts: Callable[[str, str, str], dict | None]  # (path, app_jobs.contacts_key, blob)
    status: Any


def problem(ws: Workspace, err: Exception) -> str:
    """The callout for a failed job (app_jobs.explain), naming the control to check; a limit's message names its
    own."""
    from foldmap.app import local_files

    message, control = app_jobs.explain(err, ws.name, ws.lim, ws.ss.get("last_changed"))
    kept, name = ws.ss.get("names_uploaded") or ("", "")
    if kept:  # an uploaded names file's errors name it as uploaded, not by the copy kept on the server
        message = message.replace(kept, name)
    hint = f" Check {control}." if control and not isinstance(err, app_jobs.JobLimit) else ""
    if hint and not message.rstrip().endswith((".", "!", "?")):  # the hint is a sentence of its own
        message = message.rstrip() + "."
    return said(message, path=ws.path, name=ws.name, paths=not local_files()) + hint


def guarded(ws: Workspace, job: str, blob: str, make: Callable, view_json: str = ""):
    """make() for a job ('figure', 'sequence', 'page' or 'contacts'), unless the same request went over a limit
    before: then that JobLimit at once. A JobLimit is remembered for this job (and sequence view) only; a bug is
    logged."""
    settings = json.loads(blob)
    failed = app_jobs.known_failure(ws.path, settings, job, view_json)
    if failed is not None:
        raise failed
    try:
        return make()
    except app_jobs.JobLimit as err:
        app_jobs.remember_failure(ws.path, settings, err, job, view_json)
        raise
    except (ValueError, OSError, RuntimeError):
        raise
    except Exception:  # a pipeline bug: a message for the visitor, the traceback for the server
        _LOG.exception("drawing %s failed", ws.name)
        raise


def figure_view(ws: Workspace) -> None:
    """The figure in its card; on a failure the error above the last good figure (`figure_failed` names the
    settings that failed, for the Download popover)."""
    st, ss = ws.st, ws.ss
    target = ws.blob if ss["auto_redraw"] else ss.get("drawn_blob", ws.blob)
    result = None
    try:
        with ws.status, st.spinner("Drawing…"):
            result = guarded(ws, "figure", target, lambda: ws.render(ws.path, target))
    except Exception as err:  # the canvas keeps the last good figure
        ss["figure_failed"] = target
        st.error(problem(ws, err))
    if result is not None:
        ss.pop("figure_failed", None)
        ss.update(
            drawn_blob=target,
            last_good={
                "path": ws.path,
                "blob": target,
                "svg": result["svg"],
                "files": tuple(sorted(result["files"].items())) or None,
                "scope": result.get("scope"),
                "png_dpi": result.get("png_dpi"),  # the download button's label
            },
            elements=result["elements"],
            protein_chains=result["chains"],
            na_chains=result["na_chains"],
        )
        ss["info"] = drawn_summary(ss, json.loads(target), result)
        slow = result["seconds"]  # the drawing itself, as at load (not the job's start-up or a cache hit)
        if slow > ws.lim.slow and ss.get("slow_noted") != ws.path:  # the toggle is created later in this run
            ss.update(auto_redraw=False, slow_noted=ws.path, slow_seconds=slow)
    good = ss.get("last_good") or {}
    if good.get("path") != ws.path:
        return
    scope_bar(ws, good.get("scope"), good["blob"])
    with st.container(key="canvas"):
        background = json.loads(good["blob"]).get("background", "#ffffff")
        st.markdown(figure_card(good["svg"], background, ws.name), unsafe_allow_html=True)
        notes = (
            app_edits.notes_text(
                result["fallback"],
                result["over_labels"],
                result.get("crowded_marks", []),
                result.get("page_notes", []),
                result.get("sheet_notes", []),
            )
            if result
            else ""
        )
        ions = (good.get("scope") or {}).get("na_ions")
        if ions:  # never marked: said, as the command line and MCP notes say it (C33)
            from foldmap.scope import na_ions_text

            notes = "; ".join(part for part in (notes, f"Not drawn: {na_ions_text(ions)}.") if part)
        if notes:
            st.caption(shown(notes))
        if result and result.get("domain_source_note"):  # inferred domains said to be, as the command line does
            st.caption(shown(result["domain_source_note"][0].upper() + result["domain_source_note"][1:] + "."))
        if result and result.get("domains_note"):  # as the command line warns
            note = result["domains_note"]
            st.warning(shown(f"{note[0].upper()}{note[1:]}. {DOMAINS_EDIT}"))
        missing_files(ws)
        if ss.get("check_contacts") and ss.get("na_chains"):
            contacts_check(ws, good["blob"])


DOMAINS_EDIT = "Edit Domains to remove them; an edit there checks every range."  # the layout file's domains are
# read leniently only while unchanged (app.layout_options)


def missing_files(ws: Workspace) -> None:
    """Under the figure: a layout file's names files that are not on this machine (hosted: never read), to upload."""
    from foldmap.app import local_files

    for p in ws.ss.get("names_missing") or []:
        name = shown(p, paths=not local_files())
        ws.st.warning(f"The layout file's names file {name} is not on this machine; upload it under Edit layout.")


# ---------------------------------------------------------------------------------------------- scope bar
def scope_text(scope: dict) -> str:
    """The bar's line: what an automatic choice drew, how many of the chains (when not all of them), then what of
    the nucleic acid (scope.scope_what)."""
    from foldmap.scope import scope_what

    return f"Drawn: {scope_what(scope, ' — {drawn} of {total} chains')}."


def scope_flags(scope: dict) -> str:
    """The command-line flags that draw every chain."""
    return " ".join(f"--{k} {v}" for k, v in scope["all_options"].items())


def too_large_text(scope: dict, name: str) -> str:
    return (
        f"All {scope['file_chains']} chains: too large for this app — draw with foldmap plot {name} "
        f"{scope_flags(scope)}"
    )


def too_large(ws: Workspace, scope: dict, blob: str) -> bool:
    """Would drawing every chain go over this app's residue cap, or did it go over a limit before?"""
    if ws.lim.max_residues and scope.get("residues", 0) > ws.lim.max_residues:
        return True
    every = {**json.loads(blob), **scope["all_options"]}
    return app_jobs.known_failure(ws.path, every) is not None


def show_all_chains(ss: MutableMapping, options: dict) -> None:
    """The bar's "Show all chains": the settings that draw every chain, drawn at once (auto-redraw off too)."""
    from foldmap.app import current_blob

    for key, value in options.items():
        if key == "set":  # a style key, as --set takes it: key=value
            key, value = value.split("=", 1)
        ss[key] = value
        ss["last_changed"] = key  # a failure names the control that changed
    if not ss.get("auto_redraw", True):
        ss["drawn_blob"] = current_blob(ss)


def keep_scope(ss: MutableMapping) -> None:
    """The bar's "Keep this": no bar for this structure again (app.reset_structure forgets it)."""
    ss["scope_kept"] = True


def scope_bar(ws: Workspace, scope: dict | None, blob: str) -> None:
    """Above the figure, when an automatic choice drew part of the file (scope.drawn_scope): what was drawn, a
    button that draws every chain (or why this app cannot) and one that keeps the figure as it is."""
    st, ss = ws.st, ws.ss
    if not scope or not scope.get("all_options") or ss.get("scope_kept"):  # chains asked for: nothing to offer
        return
    with st.container(key="scope_bar", horizontal=True, vertical_alignment="center", gap="small"):
        st.markdown(shown(scope_text(scope)), width="stretch")
        if too_large(ws, scope, blob):
            st.caption(shown(too_large_text(scope, ws.name)), width="stretch")
        else:
            st.button("Show all chains", key="scope_all", on_click=show_all_chains, args=(ss, scope["all_options"]))
        st.button("Keep this", key="scope_keep", type="tertiary", on_click=keep_scope, args=(ss,))


def drawn_summary(ss: MutableMapping, settings: dict, result: dict) -> str:
    """The top bar's line for the figure drawn, counting protein and nucleic-acid chains as the scope bar does. A
    figure of every chain counts the chains of its assembly (`chains_total`); one of some chains says "2 of
    6 chains" when that count is of the same assembly; one that an automatic choice drew part of
    (scope.drawn_scope) says how many chains of the file it drew."""
    assembly = settings.get("assembly") or "auto"
    if not settings.get("chains"):
        ss["chains_total"] = {"assembly": assembly, "count": result["chain_count"]}
        scope = result.get("scope") or {}
        if scope.get("reason") and scope["drawn_chains"] < scope["file_chains"]:
            return app_jobs.scoped_summary(result["summary"], scope["drawn_chains"], scope["file_chains"])
        return result["summary"]
    total = ss.get("chains_total") or {}
    if total.get("assembly") != assembly:  # not counted in this assembly yet
        return result["summary"]
    return app_jobs.restricted_summary(result["summary"], total["count"])


def contacts_check(ws: Workspace, blob: str) -> None:
    """The figure drawn (`blob`) checked against the coordinates, as captions under the card (only on request: a
    second pass over the structure)."""
    from foldmap.app import local_files
    from foldmap.app_controls import contact_lines

    st = ws.st
    try:
        with ws.status, st.spinner("Checking contacts…"):
            key = app_jobs.contacts_key(json.loads(blob))  # a colour or theme change does not check again
            report = guarded(ws, "contacts", blob, lambda: ws.contacts(ws.path, key, blob))
    except Exception as err:
        message = app_jobs.explain(err, ws.name, ws.lim)[0]
        st.warning("Contacts check: " + said(message, path=ws.path, name=ws.name, paths=not local_files()))
        return
    first, *problems = contact_lines(report)
    st.caption(first)
    for line in problems:  # they name the file's residues
        st.caption(shown(line))


def sequence_options(ss: MutableMapping) -> dict:
    """The Sequence view's options: its widgets' values (their defaults before it was first shown) and the chain
    picked, else the first protein chain."""
    chains = ss.get("protein_chains") or []
    chain = ss.get("seq_chain")
    return {
        "columns": int(ss.get("seq_columns", 60)),
        "ss_colour": ss.get("seq_ss_colour") or "figure",
        "full_sequence": bool(ss.get("seq_full", False)),
        "chain": chain if chain in chains else (chains[0] if chains else None),
    }


def sequence_json(ss: MutableMapping) -> str:
    """sequence_options as JSON: the key of the sequence caches and of its failures."""
    return json.dumps(sequence_options(ss), sort_keys=True)


def sequence_view(ws: Workspace) -> None:
    """The sequence figure, its options in one row above it (a chain picker when there are several chains)."""
    st, ss = ws.st, ws.ss
    chains = ss.get("protein_chains") or []
    if ss.get("seq_chain") not in chains:  # the chain of another structure, or of other chain settings
        ss.pop("seq_chain", None)
    with st.container(horizontal=True, vertical_alignment="top", key="seq_options"):
        st.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns", width=300, persist_state="session")
        st.segmented_control(
            "Element colours",
            ["figure", "black"],
            default="figure",
            required=True,
            key="seq_ss_colour",
            persist_state="session",
        )
        st.toggle("Unmodelled ends", key="seq_full", persist_state="session")
        if len(chains) > 1:
            st.selectbox("Chain", chains, key="seq_chain", format_func=capped, width=160, persist_state="session")
    view_json = sequence_json(ss)
    try:
        with ws.status, st.spinner("Drawing…"):
            svg = guarded(ws, "sequence", ws.blob, lambda: ws.sequence_svg(ws.path, ws.blob, view_json), view_json)
    except Exception as err:
        st.error(problem(ws, err))
        return
    ss["seq_drawn"] = (ws.blob, view_json)  # what the Download popover offers
    with st.container(key="canvas"):
        background = json.loads(ws.blob).get("background", "#ffffff")
        st.markdown(
            figure_card(svg.decode(), background, ws.name, kind="Sequence", natural=True), unsafe_allow_html=True
        )


PDB_ID = re.compile(r"[0-9][A-Za-z0-9]{3}")


def explorer_view(ws: Workspace) -> None:
    """The interactive page (figure beside the 3D structure), embedded unless larger than PAGE_EMBED_MAX."""
    st, ss = ws.st, ws.ss
    try:
        with ws.status, st.spinner("Building the 3D page…"):
            page = guarded(ws, "page", ws.blob, lambda: ws.page(ws.path, ws.blob))
    except Exception as err:
        st.error(problem(ws, err))
        return
    ss["page_drawn"] = ws.blob  # what the Download popover offers (from the same cache)
    size = len(page.encode("utf-8"))
    with st.container(key="explorer"):
        if size > app_jobs.PAGE_EMBED_MAX:
            share_id = ss.get("share_id") or ""
            rcsb = (
                f" or open {shown(share_id)} at https://www.rcsb.org/3d-view/{share_id}"
                if PDB_ID.fullmatch(share_id)
                else ""
            )
            st.info(
                f"The 3D page is {size / 1e6:.0f} MB, too large to show here. Download it (Download ▸ 3D page){rcsb}."
            )
        else:
            st.iframe(page, height=1100)  # foldmap's own page


# ---------------------------------------------------------------------------------------------- downloads
class Download(NamedTuple):
    """One button of the Download popover; data is called only when the button is clicked."""

    key: str
    label: str
    data: Callable[[], bytes]
    file_name: str
    mime: str


FORMATS = (
    ("svg", "SVG · editable", "image/svg+xml"),
    ("png", "PNG · 300 dpi", "image/png"),
    ("pdf", "PDF", "application/pdf"),
)


def png_label(dpi: int | None) -> str:
    """The PNG button's label: the dpi the figure's PNG is written at (300 when not known)."""
    return f"PNG · {dpi or 300} dpi"


FORMATS_NOTE = (
    "SVG: edit in Illustrator or Inkscape · PNG: slides and web · PDF: print · Layout file: re-render this figure "
    "later with `foldmap plot --layout-file`"
)


UPLOADS_NOTE = (
    "# Files uploaded to the app, named here as uploaded: {files}.\n"
    "# To draw this layout again, put them beside it (or write their paths); in the app, upload them again.\n"
)


def uploads_of(ss: MutableMapping) -> dict[str, str]:
    """Kept path -> uploaded name of this structure's uploads that a layout file names (names file, alignment)."""
    return dict(v for v in (ss.get("names_uploaded"), ss.get("msa_uploaded")) if v)


def layout_yaml(path, blob: str, uploads: dict[str, str] | None = None) -> bytes:
    """The layout file of these settings (YAML). Hosted, uploads are named as uploaded (app.layout_document), and
    a comment above asks for them."""
    import yaml

    from foldmap.app import layout_document, local_files, uploaded_files

    doc = layout_document(path, json.loads(blob), uploads)
    text = yaml.safe_dump(doc, sort_keys=False, allow_unicode=True)
    files = [] if local_files() else uploaded_files(doc)
    from foldmap.app_edits import upload_name

    note = UPLOADS_NOTE.format(files=", ".join(upload_name(f) for f in files)) if files else ""
    return (note + text).encode()


def share_url(base: str | None, share_id: str, theme: str) -> str:
    """The link that opens this structure in this theme: the app's address (its query string dropped) with ?id and
    ?theme only."""
    return (base or "").split("?", 1)[0] + f"?id={share_id}&theme={theme}"


def figure_downloads(
    path: str,
    blob: str,
    stem: str,
    lim: app_jobs.Limits,
    eager: tuple | None = None,
    uploads: dict[str, str] | None = None,
    png_dpi: int | None = None,
) -> list[Download]:
    """The figure's files (eager: the render's files as sorted (ext, bytes) items) and its layout file (uploads:
    uploads_of the session). png_dpi: what the figure's PNG is written at (render.png_dpi: 300, less on a page past
    50 megapixels), which its button says."""
    out = [
        Download(
            f"dl_{ext}",
            png_label(png_dpi) if ext == "png" else label,
            partial(app_jobs.figure_file, path, blob, ext, lim, eager),
            f"{stem}.{ext}",
            mime,
        )
        for ext, label, mime in FORMATS
    ]
    out.append(
        Download(
            "dl_layout",
            "Layout file",
            partial(layout_yaml, path, blob, uploads),
            f"{stem}-layout.yaml",
            "application/x-yaml",
        )
    )
    return out


def sequence_downloads(path: str, blob: str, view_json: str, stem: str, lim: app_jobs.Limits) -> list[Download]:
    """The sequence figure's files, with the Sequence view's options."""
    return [
        Download(
            f"dl_seq_{ext}",
            label,
            partial(app_jobs.sequence_file, path, blob, view_json, ext, lim),
            f"{stem}-sequence.{ext}",
            mime,
        )
        for ext, label, mime in FORMATS
    ]


def page_bytes(page: Callable[[str, str], str], path: str, blob: str) -> bytes:
    """The interactive page as a file, from `page` (the app's cache, which holds the one the 3D view showed)."""
    return page(path, blob).encode("utf-8")


def page_downloads(path: str, blob: str, stem: str, page: Callable[[str, str], str]) -> list[Download]:
    """The interactive page the 3D view built (page: the app's page cache)."""
    data = partial(page_bytes, page, path, blob)
    return [Download("dl_html", "Interactive page", data, f"{stem}-explorer.html", "text/html")]


PREVIOUS_NOTE = "These are the previous figure's files: the current settings do not draw."
LAST_DRAWN_NOTE = "These are the files of the figure as last drawn; open the Figure view to draw your changes."


def figure_offered(ss: MutableMapping, path: str, blob: str) -> tuple[str, tuple | None, str] | None:
    """The figure the Download popover offers: always the last one drawn of this structure (last_good), as
    (settings JSON, eager files, a note when it is not the figure of the current settings), else None."""
    good = ss.get("last_good") or {}
    if good.get("path") != path:
        return None
    target = blob if ss.get("auto_redraw", True) else ss.get("drawn_blob", blob)
    if ss.get("figure_failed") == target:
        note = PREVIOUS_NOTE
    else:
        note = LAST_DRAWN_NOTE if good["blob"] != target else ""
    return good["blob"], good.get("files"), note


def download_menu(ws: Workspace) -> None:
    """The Download popover's body: the files of what each view drew (each made when clicked; a view not opened yet
    offers none), what the formats are for, and the share link."""
    st, ss = ws.st, ws.ss

    def buttons(files: list[Download]) -> None:
        with st.container(horizontal=True):
            for d in files:
                st.download_button(d.label, d.data, d.file_name, d.mime, key=d.key, on_click="ignore")

    st.markdown("**Figure**")
    offered = figure_offered(ss, ws.path, ws.blob)
    if offered is None:
        st.caption("Nothing to download until the figure draws.")
    else:
        blob, eager, note = offered
        if note:
            st.caption(note)
        dpi = (ss.get("last_good") or {}).get("png_dpi")
        buttons(figure_downloads(ws.path, blob, ws.name, ws.lim, eager, uploads_of(ss), png_dpi=dpi))
    st.markdown("**Sequence**")
    if ss.get("seq_drawn"):
        blob, view_json = ss["seq_drawn"]
        buttons(sequence_downloads(ws.path, blob, view_json, ws.name, ws.lim))
    else:
        st.caption("Open the Sequence view to make this file.")
    st.markdown("**3D**")
    if ss.get("page_drawn"):
        buttons(page_downloads(ws.path, ss["page_drawn"], ws.name, ws.page))
    else:
        st.caption("Open the 3D view to make this file.")
    st.caption(FORMATS_NOTE)
    if ss.get("share_id"):
        st.caption("Link to this figure")
        st.code(share_url(st.context.url, ss["share_id"], ss["theme"]), language=None)
