"""Driving the Streamlit app headlessly (streamlit.testing AppTest)."""

import shutil
from pathlib import Path

APP = Path(__file__).parents[1] / "src" / "foldmap" / "app.py"
DATA = Path(__file__).parent / "data"


def app_test(timeout: float = 120):
    from streamlit.testing.v1 import AppTest

    return AppTest.from_file(str(APP), default_timeout=timeout)


def load(at, source):
    at.run()
    at.text_input(key="source").input(str(source))
    at.button(key="load").click()
    at.run()
    return at


def opened(source):
    at = load(app_test(), source)
    assert not at.exception and at.session_state["stage"] == "edit", [e.value for e in at.error]
    return at


def svg_of(at) -> str:
    return next(m.value for m in at.markdown if "<svg" in m.value)


def as_text(at, key, value):
    """Type `value` into an editor's "as text" field (the same syntax as the command line)."""
    at.toggle(key=f"{key}__as_text").set_value(True).run()
    (at.text_area(key=key) if key == "domains" else at.text_input(key=key)).set_value(value).run()
    return at


def seed_cache(monkeypatch, tmp_path, *ids):
    """fetch() finds these PDB entries in a private cache: no network."""
    import foldmap.fetch as fetch_mod

    monkeypatch.setattr(fetch_mod, "CACHE", tmp_path)
    for i in ids:
        shutil.copy(DATA / f"{i}.cif", tmp_path / f"{i}.cif")


def widget_keys(block, skip=("expander", "popover", "dialog")):
    from streamlit.testing.v1.element_tree import Block, Widget

    out = []
    for child in getattr(block, "children", {}).values():
        if isinstance(child, Block) and child.type in skip:
            continue
        if isinstance(child, Widget):
            out.append(child.key)
        out += widget_keys(child, skip)
    return out


def fresh_caches():
    import streamlit as st

    from foldmap import app_jobs

    st.cache_data.clear()
    app_jobs.clear()


def count_calls(monkeypatch, module, name) -> list:
    calls, real = [], getattr(module, name)
    monkeypatch.setattr(module, name, lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    return calls


def visible(text) -> str:
    """What a reader sees of the app's Markdown-escaped text (app_edits.shown): escapes resolved, zero-width spaces
    gone."""
    import re

    return re.sub(r"\\(.)", r"\1", text).replace("\u200b", "")


def rendered_text(at) -> str:
    """Every text the app put on the page (markdown, captions, callouts, headings, code), one per line."""
    kinds = ("markdown", "caption", "error", "warning", "info", "success", "text", "title", "header", "subheader")
    return "\n".join(str(e.value) for kind in (*kinds, "code") for e in getattr(at, kind))


def serve(monkeypatch, files: dict[str, bytes]) -> list[str]:
    """Answer the app's downloads offline: a URL ending in a key of `files` gets its bytes, any other the PDB's 404
    (urllib raises HTTPError, as on an obsolete or mistyped id). Returns the URLs asked for."""
    import io
    import urllib.error

    import foldmap.fetch as fetch_mod

    asked = []

    class Reply(io.BytesIO):
        headers: dict = {}

    def urlopen(request, timeout=0):
        url = getattr(request, "full_url", request)
        asked.append(url)
        for end, data in files.items():
            if url.endswith(end):
                return Reply(data)
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(fetch_mod.urllib.request, "urlopen", urlopen)
    return asked


def keep_downloads(monkeypatch) -> dict:
    """The data callables of the download buttons as the app hands them to Streamlit (st.download_button with a
    callable: MediaFileManager.add_deferred), by the file id the button carries. A click runs that callable on the
    server (execute_deferred); AppTest has no server, so `press` runs it. `monkeypatch`: pytest's, or a
    pytest.MonkeyPatch.context()."""
    from streamlit.runtime.media_file_manager import MediaFileManager

    kept: dict = {}
    real = MediaFileManager.add_deferred

    def add_deferred(self, data_callable, *args, **kwargs):
        file_id = real(self, data_callable, *args, **kwargs)
        kept[file_id] = data_callable
        return file_id

    monkeypatch.setattr(MediaFileManager, "add_deferred", add_deferred)
    return kept


def press(at, key: str, kept: dict) -> bytes:
    """What a click on the download button `key` gives the visitor: its callable run, as the server runs it."""
    data = kept[at.download_button(key=key).proto.deferred_file_id]()
    return data.encode() if isinstance(data, str) else bytes(data)
