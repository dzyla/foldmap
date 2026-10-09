"""The Streamlit app, driven headlessly: load a structure, pick a style, edit, download."""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = Path(__file__).parents[1] / "src" / "foldmap" / "app.py"
DATA = Path(__file__).parent / "data"


def _app():
    return AppTest.from_file(str(APP), default_timeout=120)


def _load(at, source):
    at.run()
    at.text_input(key="source").input(str(source))
    at.button(key="load").click()
    at.run()
    return at


def _svg(at) -> str:
    return next(m.value for m in at.markdown if "<svg" in m.value)


def test_first_screen_asks_for_a_structure():
    at = _app()
    at.run()
    assert not at.exception
    assert at.text_input(key="source") is not None and at.session_state["stage"] == "load"


def test_loading_leads_to_style_choice_with_previews():
    at = _load(_app(), DATA / "1UBQ.cif")
    assert not at.exception and at.session_state["stage"] == "style"
    picks = [b for b in at.button if (b.key or "").startswith("pick_")]
    assert len(picks) >= 6


def test_picking_a_style_opens_the_editor_with_that_theme():
    at = _load(_app(), DATA / "1LMB.cif")
    at.button(key="pick_trace").click()
    at.run()
    assert not at.exception and at.session_state["stage"] == "edit"
    assert at.selectbox(key="theme").value == "trace"
    assert at.session_state["loop_arrows"] is True and "<svg" in _svg(at)


def test_controls_change_the_figure():
    at = _load(_app(), DATA / "1UBQ.cif")
    at.button(key="pick_publication").click()
    at.run()
    assert "resnum:" not in _svg(at)
    at.checkbox(key="residue_numbers").check()
    at.run()
    assert "resnum:" in _svg(at)
    at.selectbox(key="mode").select("stack")
    at.run()
    assert not at.exception and "<svg" in _svg(at)


def test_bad_input_is_reported_not_raised():
    at = _load(_app(), "/nowhere/missing.cif")
    assert not at.exception and at.session_state["stage"] == "load" and at.error


def test_layout_errors_are_shown_in_the_editor():
    at = _load(_app(), DATA / "2LZM.cif")
    at.button(key="pick_publication").click()
    at.run()
    at.text_input(key="swap").input("α1")  # one element: not a swap
    at.run()
    assert not at.exception and any("swap" in e.value for e in at.error)


def test_download_buttons_are_offered():
    from foldmap.app import figure_files

    files = figure_files(DATA / "1UBQ.cif", {"theme": "publication"})
    assert set(files) == {"svg", "png", "pdf"}
    assert files["png"][:4] == b"\x89PNG" and files["pdf"][:4] == b"%PDF" and b"<svg" in files["svg"]


def test_app_passes_membrane_and_alignment_to_the_layout():
    from foldmap.app import layout_options

    opts = layout_options({"membrane": "auto", "msa_path": str(DATA / "ubq_family.fasta"), "msa_reference": "1UBQ"})
    assert opts["membrane"] == "auto" and opts["msa"] == str(DATA / "ubq_family.fasta")
    assert opts["msa_reference"] == "1UBQ"
    assert layout_options({})["msa"] is None


def test_app_sequence_files():
    from foldmap.app import sequence_files

    out = sequence_files(DATA / "1UBQ.cif", {"theme": "publication"}, {"columns": 30, "ss_colour": "figure"})
    assert out["svg"].startswith(b"<svg") and out["png"][:4] == b"\x89PNG" and out["pdf"][:4] == b"%PDF"
    msa = sequence_files(
        DATA / "1UBQ.cif", {"theme": "publication", "msa_path": str(DATA / "ubq_family.fasta")}, {"columns": 40}
    )
    assert b"strict:" in msa["svg"]


def test_app_shows_the_sequence_tab():
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.text_input(key="source").set_value(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    at.button(key="pick_publication").click().run()
    assert not at.exception
    assert any(t.label == "Sequence" for t in at.tabs)
    assert at.slider(key="seq_columns").value == 60


def test_app_passes_focus_to_the_layout():
    from foldmap.app import layout_options

    assert layout_options({"focus": "none"})["focus"] == "none"
    assert layout_options({})["focus"] == "auto"


def test_app_lists_structures_for_a_uniprot_accession(monkeypatch):
    import foldmap.app as app
    from foldmap import uniprot

    e = uniprot.parse_entry(__import__("json").loads((DATA / "uniprot" / "P07911.json").read_text()))
    monkeypatch.setattr(uniprot, "entry", lambda acc, cache_dir=None: e)
    choices = app.structure_choices("P07911")
    assert choices[0][0] == "AF-P07911-F1" and "AlphaFold" in choices[0][1]
    assert len(choices) > 2 and all(len(c[0]) == 4 for c in choices[1:])  # PDB IDs, widest coverage first
    assert app.structure_choices("1LMB") == []


def test_app_passes_uniprot_options():
    from foldmap.app import layout_options

    opts = layout_options({"uniprot": "auto", "uniprot_domains": True})
    assert opts["uniprot"] == "auto" and opts["domains"] == "uniprot"
    assert layout_options({})["uniprot"] is None
