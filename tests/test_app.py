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


def test_style_previews_are_static_files_for_every_theme():
    from foldmap.app import PREVIEW_THEMES, preview_image

    for theme in PREVIEW_THEMES:
        data = preview_image(theme)
        assert data[:8] == b"\x89PNG\r\n\x1a\n" and len(data) < 120_000, theme


def test_style_stage_draws_nothing(monkeypatch):
    """Picking a style must not render the loaded structure once per theme."""
    from streamlit.testing.v1 import AppTest

    import foldmap.spec as spec

    calls = []
    real = spec.draw_spec
    monkeypatch.setattr(spec, "draw_spec", lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    at = AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.text_input(key="source").set_value(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    assert at.session_state["stage"] == "style" and not at.exception
    assert len(at.image) == len(__import__("foldmap.app", fromlist=["x"]).PREVIEW_THEMES)
    assert calls == []


def test_file_that_is_not_a_structure_is_reported(tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text("just some notes\n", encoding="utf-8")
    at = _load(_app(), notes)
    assert not at.exception and at.session_state["stage"] == "load" and at.error


def test_app_offers_flip_for_the_membrane():
    from foldmap.app import layout_options

    opts = layout_options({"membrane": "flip"})
    assert opts["membrane"] == "flip"


def test_app_accepts_a_consurf_table():
    """A ConSurf grades table colours by conservation; the Sequence tab draws the plain sequence beside it."""
    from foldmap.app import sequence_files

    table = DATA / "consurf" / "1TIM_A_grades.txt"
    at = _load(_app(), DATA / "1TIM.cif")
    at.button(key="pick_publication").click().run()
    at.file_uploader(key="msa_upload").set_value(("1TIM_A_grades.txt", table.read_bytes(), "text/plain")).run()
    assert not at.exception and not at.error
    assert at.session_state["msa_error"] == "" and at.session_state["msa_path"].endswith(".txt")
    assert any("ConSurf" in c.value for c in at.caption)
    out = sequence_files(DATA / "1TIM.cif", {"theme": "publication", "msa_path": str(table)}, {"columns": 40})
    assert out["svg"].startswith(b"<svg") and b"strict:" not in out["svg"]


def test_app_layout_file_keeps_every_layout_option():
    from foldmap.app import layout_document, settings_from_document

    msa = str(DATA / "ubq_family.fasta")
    settings = {"theme": "publication", "focus": "none", "chains": "A", "msa_path": msa, "msa_reference": "1UBQ_A"}
    doc = layout_document(DATA / "1UBQ.cif", settings)
    assert doc["layout"]["focus"] == "none" and doc["layout"]["chains"] == ["A"] and doc["layout"]["msa"] == msa
    back = settings_from_document(doc)
    assert back["focus"] == "none" and back["chains"] == "A" and back["msa_path"] == msa
    assert back["msa_reference"] == "1UBQ_A" and back["uniprot"] == "off" and back["membrane"] == "off"
    back = settings_from_document({"layout": {"membrane": "auto", "uniprot": "P0CG48", "msa": "/nowhere/x.fasta"}})
    assert back["membrane"] == "auto" and back["uniprot"] == "P0CG48" and back["msa_path"] == ""


def test_app_layout_file_round_trips_through_the_app():
    """Download the layout file with every option set, load it into a fresh app: every option comes back."""
    import yaml

    from foldmap.app import layout_document

    msa = str(DATA / "ubq_family.fasta")
    settings = {
        "theme": "trace", "membrane": "flip", "focus": "none", "chains": "A", "msa_path": msa,
        "msa_reference": "1UBQ_A",
        "flip_v": True, "symmetry_tol": 2.5, "up": [0, 1, 0], "view": [0, 0, 1],
        "rename": "res:A:23-34=Hx", "swap": "", "move": "", "dna": False,
        "domains": "N=res:A:1-30", "rotate": 30, "title": "Round trip",
    }  # fmt: skip
    text = yaml.safe_dump(layout_document(DATA / "1UBQ.cif", settings), sort_keys=False)
    at = _load(_app(), DATA / "1UBQ.cif")
    at.button(key="pick_publication").click().run()
    at.file_uploader(key="layout_upload").set_value(("layout.yaml", text.encode(), "application/x-yaml")).run()
    assert not at.exception and not at.error
    ss = at.session_state
    assert ss["theme"] == "trace" and ss["membrane"] == "flip" and ss["focus"] == "none" and ss["chains"] == "A"
    assert ss["msa_path"] == msa and ss["msa_reference"] == "1UBQ_A" and ss["dna"] is False
    assert ss["rotate"] == 30 and ss["title"] == "Round trip" and "Hx" in ss["rename"] and "N=" in ss["domains"]
    from foldmap.app import LAYOUT_DEFAULTS, STYLE_FIELDS

    again = layout_document(
        DATA / "1UBQ.cif", {k: ss[k] for k in [*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme"] if k in ss}
    )
    first = yaml.safe_load(text)
    for key in ("membrane", "focus", "chains", "msa", "msa_reference", "dna", "flip_v", "symmetry_tol", "up", "view"):
        assert again["layout"][key] == first["layout"][key], key
    assert again["domains"] == first["domains"] and again["edits"] == first["edits"]


def test_app_sequence_tab_with_the_conservation_theme(tmp_path):
    """The conservation theme colours the sequence view's elements from the same alignment as the figure."""
    at = _load(_app(), DATA / "1UBQ.cif")
    at.button(key="pick_publication").click().run()
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    at.selectbox(key="theme").select("conservation").run()
    assert not at.exception and at.session_state["msa_path"]
    assert not at.error and any("strict:" in m.value for m in at.markdown)


def test_app_layout_file_options_without_a_control_do_not_carry_to_the_next_structure(tmp_path):
    import yaml

    doc = {"foldmap": 1, "layout": {"flip_v": True, "protomers": [["A", "B"], ["C", "D"]], "msa": "/nowhere/x.fasta"}}
    at = _load(_app(), DATA / "2HHB.cif")
    at.button(key="pick_publication").click().run()
    up = ("layout.yaml", yaml.safe_dump(doc).encode(), "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    assert at.session_state["flip_v"] is True and at.session_state["protomers"] == [["A", "B"], ["C", "D"]]
    assert any("/nowhere/x.fasta" in w.value for w in at.warning)
    at.button(key="back_edit").click().run()
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    at.button(key="pick_publication").click().run()
    assert not at.exception and not at.error
    assert at.session_state["flip_v"] is False and at.session_state["protomers"] is None
    assert not at.warning and "<svg" in _svg(at)


def test_app_edits_of_one_structure_do_not_carry_to_the_next():
    at = _load(_app(), DATA / "2HHB.cif")
    at.button(key="pick_publication").click().run()
    at.text_input(key="rename").set_value("res:B:1-10=Zz")
    at.text_input(key="chains").set_value("B,D").run()
    assert not at.exception and not at.error
    at.button(key="back_edit").click().run()
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    at.button(key="pick_publication").click().run()
    assert not at.exception and not at.error
    ss = at.session_state
    assert ss["rename"] == "" and ss["chains"] == "" and ss["domains"] == "" and "<svg" in _svg(at)
