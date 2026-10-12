"""Driving the Streamlit app headlessly (streamlit.testing AppTest)."""

import shutil
from pathlib import Path

APP = Path(__file__).parents[1] / "src" / "foldmap" / "app.py"
DATA = Path(__file__).parent / "data"


def app_test():
    from streamlit.testing.v1 import AppTest

    return AppTest.from_file(str(APP), default_timeout=120)


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


def rendered_text(at) -> str:
    """Every text the app put on the page (markdown, captions, callouts, headings, code), one per line."""
    kinds = ("markdown", "caption", "error", "warning", "info", "success", "text", "title", "header", "subheader")
    return "\n".join(str(e.value) for kind in (*kinds, "code") for e in getattr(at, kind))
