"""The Streamlit app, driven headlessly: Start, load (drawn at once), share links, edit, download."""

import re
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from app_helpers import DATA, app_test, load, opened, seed_cache, svg_of  # noqa: E402


def _visible(text):
    """What a reader sees of the app's Markdown-escaped text: escapes resolved, zero-width spaces gone."""
    return re.sub(r"\\(.)", r"\1", text).replace("\u200b", "")


def test_first_screen_asks_for_a_structure():
    at = app_test()
    at.run()
    assert not at.exception
    assert at.text_input(key="source") is not None and at.session_state["stage"] == "load"


def test_loading_draws_at_once():
    """Was test_loading_leads_to_style_choice_with_previews: no style stage, the figure is drawn on load."""
    at = opened(DATA / "1UBQ.cif")
    assert at.session_state["theme"] == "publication" and "<svg" in svg_of(at)
    assert at.session_state["info"] == "1 chain · 76 residues · 1 helix · 5 strands"
    assert not at.image  # no thumbnails


def test_picking_a_theme_redraws_with_it():
    """Was test_picking_a_style_opens_the_editor_with_that_theme (the dialog: test_app_controls.py)."""
    at = opened(DATA / "1LMB.cif")
    at.selectbox(key="theme").select("trace").run()
    assert not at.exception and at.session_state["loop_arrows"] is True and "<svg" in svg_of(at)


def test_mode_control_redraws():
    """The mode half of test_controls_change_the_figure (its residue-numbers half is in test_app_controls.py);
    Task 9 restates it as test_layout_mode_and_sizes."""
    at = opened(DATA / "1UBQ.cif")
    at.selectbox(key="mode").select("stack").run()
    assert not at.exception and "<svg" in svg_of(at)


def test_nucleic_acid_style_control_redraws_the_stems():
    from foldmap.style import CHOICES

    at = load(app_test(), DATA / "1LMB.cif")
    style = at.selectbox(key="na_style")
    assert list(style.options) == list(CHOICES["na_style"]) and style.value == "ladder"
    assert 'id="na:stem:0:rails"' in svg_of(at) and ":block" not in svg_of(at)
    style.select("block")
    at.run()
    assert not at.exception
    assert 'id="na:stem:0:block"' in svg_of(at) and ":rails" not in svg_of(at)


def test_bad_input_is_reported_not_raised():
    at = load(app_test(), "/nowhere/missing.cif")
    assert not at.exception and at.session_state["stage"] == "load" and at.error


def test_layout_errors_are_shown_in_the_editor():
    at = load(app_test(), DATA / "2LZM.cif")
    at.text_input(key="swap").input("α1")  # one element: not a swap
    at.run()
    assert not at.exception and any("swap" in e.value for e in at.error)


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


def test_app_shows_the_sequence_view():
    at = app_test()
    at.run()
    at.text_input(key="source").set_value(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    assert not at.exception
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not at.exception and at.slider(key="seq_columns").value == 60


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


def test_file_that_is_not_a_structure_is_reported(tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text("just some notes\n", encoding="utf-8")
    at = load(app_test(), notes)
    assert not at.exception and at.session_state["stage"] == "load" and at.error
    assert at.error[0].value == "notes could not be read as a structure (mmCIF or PDB)"


def test_an_upload_that_is_not_a_structure_names_no_server_path(monkeypatch):
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("notes.pdb", b"just some notes\n", "chemical/x-pdb")).run()
    assert not at.exception and at.session_state["stage"] == "load"
    text = at.error[0].value
    assert _visible(text) == "notes: no protein residues with a complete backbone" and "uploads" not in text


def test_app_offers_flip_for_the_membrane():
    from foldmap.app import layout_options

    opts = layout_options({"membrane": "flip"})
    assert opts["membrane"] == "flip"


def test_app_accepts_a_consurf_table():
    """A ConSurf grades table colours by conservation; the Sequence tab draws the plain sequence beside it."""
    from foldmap.app import sequence_files

    table = DATA / "consurf" / "1TIM_A_grades.txt"
    at = load(app_test(), DATA / "1TIM.cif")
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
    at = load(app_test(), DATA / "1UBQ.cif")
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


def test_app_sequence_view_with_the_conservation_theme(tmp_path):
    """The conservation theme colours the sequence view's elements from the same alignment as the figure."""
    at = load(app_test(), DATA / "1UBQ.cif")
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    at.selectbox(key="theme").select("conservation").run()
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not at.exception and at.session_state["msa_path"]
    assert not at.error and any("strict:" in m.value for m in at.markdown)


def test_app_layout_file_options_without_a_control_do_not_carry_to_the_next_structure(tmp_path):
    import yaml

    doc = {"foldmap": 1, "layout": {"flip_v": True, "protomers": [["A", "B"], ["C", "D"]], "msa": "/nowhere/x.fasta"}}
    at = load(app_test(), DATA / "2HHB.cif")
    up = ("layout.yaml", yaml.safe_dump(doc).encode(), "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    assert at.session_state["flip_v"] is True and at.session_state["protomers"] == [["A", "B"], ["C", "D"]]
    assert any("/nowhere/x.fasta" in w.value for w in at.warning)
    at.button(key="new_structure").click().run()
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    assert not at.exception and not at.error
    assert at.session_state["flip_v"] is False and at.session_state["protomers"] is None
    assert not at.warning and "<svg" in svg_of(at)


def test_app_edits_of_one_structure_do_not_carry_to_the_next():
    at = load(app_test(), DATA / "2HHB.cif")
    at.text_input(key="rename").set_value("res:B:1-10=Zz")
    at.text_input(key="chains").set_value("B,D").run()
    assert not at.exception and not at.error
    at.button(key="new_structure").click().run()
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    assert not at.exception and not at.error
    ss = at.session_state
    assert ss["rename"] == "" and ss["chains"] == "" and ss["domains"] == "" and "<svg" in svg_of(at)


def test_hosted_app_refuses_file_paths(monkeypatch):
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    at = load(app_test(), DATA / "1UBQ.cif")
    assert at.session_state["stage"] == "load" and "is not a PDB ID" in at.error[0].value


def test_hosted_app_ignores_alignment_paths_in_layout_files(monkeypatch):
    from foldmap.app import settings_from_document

    doc = {"layout": {"msa": str(DATA / "ubq_family.fasta")}}
    assert settings_from_document(doc)["msa_path"] == str(DATA / "ubq_family.fasta")
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "off")
    back = settings_from_document(doc)
    assert back["msa_path"] == "" and back["msa_missing"] == str(DATA / "ubq_family.fasta")


def test_app_alignment_of_one_structure_does_not_carry_to_the_next():
    at = load(app_test(), DATA / "1UBQ.cif")
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    assert at.session_state["msa_path"]
    at.button(key="new_structure").click().run()
    at.text_input(key="source").input(str(DATA / "1LMB.cif"))
    at.button(key="load").click().run()
    assert not at.exception and not at.error
    assert at.session_state["msa_path"] == "" and at.session_state["msa_error"] == ""
    assert at.session_state["msa_reference"] == ""


def test_app_pruned_upload_is_a_message_not_a_traceback():
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "chemical/x-cif")).run()
    assert not at.exception and not at.error
    kept = Path(at.session_state["path"])
    assert kept.is_file()
    kept.unlink()  # a long hosted session: the upload was pruned
    at.selectbox(key="mode").select("stack").run()
    assert not at.exception
    assert any("load it again" in e.value and str(kept.parent) not in e.value for e in at.error)


def test_app_pruned_upload_in_the_sequence_view_is_a_message(monkeypatch):
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "chemical/x-cif")).run()
    at.segmented_control(key="pane").set_value("Sequence").run()
    Path(at.session_state["path"]).unlink()
    at.slider(key="seq_columns").set_value(80).run()  # the figure is cached; the sequence is drawn again
    assert not at.exception and any("load it again" in e.value for e in at.error)


def test_app_alignment_that_cannot_be_stored_is_a_message(monkeypatch):
    import foldmap.fetch as fetch_mod

    def fail(*a, **k):
        raise OSError("disk full: /srv/secret/path")

    at = load(app_test(), DATA / "1UBQ.cif")
    monkeypatch.setattr(fetch_mod, "atomic_write", fail)
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    assert not at.exception and at.session_state["msa_path"] == ""
    assert at.session_state["msa_error"] and "/srv/secret" not in at.session_state["msa_error"]


def test_app_reads_hand_written_layout_values_like_the_command_line():
    """Unquoted `symmetry: off` is YAML false; up, view and protomers may be written as text."""
    from foldmap.app import settings_from_document

    doc = {"layout": {"symmetry": False, "membrane": False, "up": "0,0,1", "view": [1, 0, 0], "protomers": "A,B;C,D"}}
    got = settings_from_document(doc)
    assert got["symmetry"] == "off" and got["membrane"] == "off"
    assert (
        got["up"] == [0.0, 0.0, 1.0] and got["view"] == [1.0, 0.0, 0.0] and got["protomers"] == [["A", "B"], ["C", "D"]]
    )
    assert settings_from_document({"layout": {"flip_v": "no", "dna": "no"}})["flip_v"] is False
    with pytest.raises(ValueError, match="needs three numbers"):
        settings_from_document({"layout": {"up": "1,2"}})


def test_app_loads_a_hand_written_layout_file():
    text = "foldmap: 1\nlayout:\n  symmetry: off\n  membrane: off\n  up: '0,0,1'\n  protomers: 'A'\n  flip_v: 'no'\n"
    at = load(app_test(), DATA / "1UBQ.cif")
    at.file_uploader(key="layout_upload").set_value(("layout.yaml", text.encode(), "application/x-yaml")).run()
    assert not at.exception and not at.error
    assert at.session_state["symmetry"] == "off" and at.session_state["up"] == [0.0, 0.0, 1.0]
    assert "<svg" in svg_of(at)


def test_app_reports_a_layout_file_with_a_bad_value():
    at = load(app_test(), DATA / "1UBQ.cif")
    text = b"foldmap: 1\nlayout: {up: '1,2', flip_v: maybe}\n"
    at.file_uploader(key="layout_upload").set_value(("layout.yaml", text, "application/x-yaml")).run()
    assert not at.exception and any("three numbers" in e.value for e in at.error)
    bad = b"layout: {rotate: 30\n"
    at.file_uploader(key="layout_upload").set_value(("layout.yaml", bad, "application/x-yaml")).run()
    assert not at.exception and any("not valid YAML" in e.value for e in at.error)


def test_app_layout_file_keeps_a_fractional_rotate():
    import yaml

    from foldmap.app import settings_from_document

    assert settings_from_document({"layout": {"rotate": 12.5}})["rotate"] == 12.5
    doc = {"foldmap": 1, "layout": {"rotate": 12.5}}
    at = load(app_test(), DATA / "1UBQ.cif")
    up = ("layout.yaml", yaml.safe_dump(doc).encode(), "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    assert not at.exception and not at.error
    assert at.session_state["rotate"] == 12.5
    old = {"foldmap": 1, "layout": {"rotate": 30}}
    at.file_uploader(key="layout_upload").set_value(
        ("l.yaml", yaml.safe_dump(old).encode(), "application/x-yaml")
    ).run()
    assert not at.exception and at.session_state["rotate"] == 30


def test_example_loads_in_one_action(monkeypatch, tmp_path):
    from foldmap.app import EXAMPLES

    seed_cache(monkeypatch, tmp_path, "1LMB")
    at = app_test()
    at.run()
    labels = at.pills(key="example").options  # AppTest lists the formatted labels
    assert labels[0] == "1LMB · λ repressor–DNA" and [o.split(" · ")[0] for o in labels] == list(EXAMPLES)
    at.pills(key="example").set_value("1LMB").run()
    assert not at.exception and at.session_state["stage"] == "edit" and "<svg" in svg_of(at)
    assert at.query_params["id"] == "1LMB" and at.query_params["theme"] == "publication"


def test_share_link_opens_the_workspace(monkeypatch, tmp_path):
    seed_cache(monkeypatch, tmp_path, "1LMB")
    at = app_test()
    at.query_params["id"] = "1LMB"
    at.query_params["theme"] = "journal"
    at.run()
    assert not at.exception and at.session_state["stage"] == "edit"
    assert at.session_state["theme"] == "journal" and at.session_state["share_id"] == "1LMB"


def test_uploads_and_paths_never_enter_the_link():
    at = opened(DATA / "1UBQ.cif")
    assert "id" not in at.query_params and at.session_state["share_id"] == ""
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "chemical/x-cif")).run()
    assert at.session_state["stage"] == "edit" and "id" not in at.query_params


def test_a_uniprot_accession_lists_its_structures(monkeypatch):
    import json

    from foldmap import uniprot

    e = uniprot.parse_entry(json.loads((DATA / "uniprot" / "P07911.json").read_text()))
    monkeypatch.setattr(uniprot, "entry", lambda acc, cache_dir=None: e)
    at = load(app_test(), "P07911")
    assert not at.exception and at.session_state["stage"] == "load" and not at.error
    ids = at.dataframe[0].value["id"].tolist()
    assert ids[0] == "AF-P07911-F1" and len(ids) > 2


def test_a_giant_entry_is_refused_with_a_clear_message(monkeypatch):
    monkeypatch.setenv("FOLDMAP_APP_MAX_RESIDUES", "50")
    at = load(app_test(), DATA / "1UBQ.cif")
    assert at.session_state["stage"] == "load"
    text = at.error[0].value
    assert _visible(text).startswith("1UBQ is too large for this app: 76 residues in 1 chain(s)")
    assert "foldmap plot 1UBQ -o 1UBQ.svg" in _visible(text)


def test_start_shows_the_credit_line():
    from foldmap import CREDIT

    at = app_test()
    at.run()
    assert any(c.value == CREDIT for c in at.caption)


def test_the_look_carries_to_the_next_structure_and_its_parts_do_not():
    at = opened(DATA / "2HHB.cif")
    at.selectbox(key="theme").select("trace").run()
    at.slider(key="font_scale").set_value(1.4)
    at.text_input(key="mark").set_value("res:A:1-10=red").run()
    at.segmented_control(key="pane").set_value("Sequence").run()
    at.selectbox(key="seq_chain").select("D").run()  # Sequence-view and names widgets of this structure
    at.segmented_control(key="pane").set_value("3D").run()
    assert at.session_state["seq_drawn"] and at.session_state["page_drawn"]  # what the views drew
    at.segmented_control(key="pane").set_value("Figure").run()
    at.text_input(key="swap").input("α1").run()  # a failing setting: the popover offers the last good figure
    assert at.session_state["figure_failed"]
    at.session_state["names__kw"] = ["beta"]  # (Task 10)
    at.button(key="new_structure").click().run()
    assert at.session_state["stage"] == "load" and "id" not in at.query_params
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    ss = at.session_state
    assert not at.exception and ss["stage"] == "edit" and ss["name"] == "1UBQ"
    assert ss["theme"] == "trace" and ss["font_scale"] == 1.4 and ss["mark"] == ""
    assert "seq_chain" not in ss and "names__kw" not in ss and ss["auto_redraw"] is True
    assert "figure_failed" not in ss and "seq_drawn" not in ss and "page_drawn" not in ss and ss["pane"] == "Figure"
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not at.exception and "seq_chain" not in [s.key for s in at.selectbox]


def test_a_share_link_opens_ids_only(monkeypatch):
    at = app_test()
    at.query_params["id"] = str(DATA / "1UBQ.cif")
    at.query_params["theme"] = "nonsense"
    at.run()
    assert not at.exception and at.session_state["stage"] == "load"
    assert "not a PDB ID" in at.error[0].value and "id" not in at.query_params


def test_a_picked_table_row_is_the_structure_to_load():
    from foldmap.app import picked_row

    rows = [("AF-P07911-F1", "AlphaFold model"), ("1ABC", "X-ray")]
    assert picked_row({"selection": {"rows": [1], "columns": []}}, rows) == "1ABC"
    assert picked_row({"selection": {"rows": []}}, rows) is None and picked_row(None, rows) is None


def test_a_hosted_upload_that_gemmi_refuses_names_no_server_path(monkeypatch):
    from pathlib import Path as P

    from app_helpers import rendered_text

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    for data in (b"garbage\n", b"data_a\n_cell.length_a 1\n_x\n"):
        at = app_test()
        at.run()
        at.file_uploader(key="upload").set_value(("bad.cif", data, "chemical/x-cif")).run()
        text = rendered_text(at)
        assert not at.exception and at.session_state["stage"] == "load" and at.error
        assert "uploads" not in text and str(P.home()) not in text and "pytest-of" not in text, text


def test_crafted_ids_and_file_names_carry_no_markdown(monkeypatch):
    from app_helpers import rendered_text

    at = app_test()
    at.query_params["id"] = "![x](https://example.invalid/p.png) [log in again](https://example.invalid)"
    at.run()
    assert [e.value for e in at.error] == ["The link's id is not a PDB ID or UniProt accession."]
    assert "](" not in rendered_text(at) and "![" not in rendered_text(at)
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    for data in (b"garbage\n", (DATA / "1UBQ.cif").read_bytes()):
        at = app_test()
        at.run()
        at.file_uploader(key="upload").set_value(("![x](http://e) [y](http://e).cif", data, "chemical/x-cif")).run()
        assert not at.exception and at.error or at.session_state["stage"] == "edit"
        assert "](" not in rendered_text(at) and "![" not in rendered_text(at)


def test_an_upload_that_cannot_be_stored_says_so(monkeypatch):
    import foldmap.fetch as fetch_mod

    def full(*a, **k):
        raise OSError(28, "No space left on device: /srv/secret")

    monkeypatch.setattr(fetch_mod, "atomic_write", full)
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "chemical/x-cif")).run()
    assert [e.value for e in at.error] == ["The file could not be stored on the server; try again."]


def test_a_share_link_with_an_unknown_theme_opens_in_publication(monkeypatch, tmp_path):
    seed_cache(monkeypatch, tmp_path, "1LMB")
    at = app_test()
    at.query_params["id"] = "1LMB"
    at.query_params["theme"] = "nonsense"
    at.run()
    assert not at.exception and at.session_state["stage"] == "edit" and at.session_state["theme"] == "publication"
    assert at.query_params["id"] == "1LMB" and at.query_params["theme"] == "publication"


def test_an_app_made_link_opens_the_same_figure(monkeypatch, tmp_path):
    seed_cache(monkeypatch, tmp_path, "1LMB")
    at = app_test()
    at.query_params["id"] = "1lmb"
    at.query_params["theme"] = "journal"
    at.run()
    link = {k: at.query_params[k] for k in at.query_params}
    assert link == {"id": "1LMB", "theme": "journal"}
    again = app_test()
    for k, v in link.items():
        again.query_params[k] = v
    again.run()
    assert not again.exception and again.session_state["stage"] == "edit"
    assert again.session_state["theme"] == "journal" and svg_of(again) == svg_of(at)
    assert {k: again.query_params[k] for k in again.query_params} == link


def test_a_capped_structure_is_never_drawn(monkeypatch):
    from app_helpers import count_calls, fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    calls = count_calls(monkeypatch, app_jobs, "draw_job")
    monkeypatch.setenv("FOLDMAP_APP_MAX_RESIDUES", "50")
    at = load(app_test(), DATA / "1UBQ.cif")
    assert at.session_state["stage"] == "load" and calls == []
    monkeypatch.delenv("FOLDMAP_APP_MAX_RESIDUES")
    at = load(app_test(), DATA / "1UBQ.cif")
    assert at.session_state["stage"] == "edit" and len(calls) == 1


def test_the_hosted_app_uses_the_hosted_size_cap(monkeypatch):
    """No FOLDMAP_APP_MAX_RESIDUES: the cap comes from app_jobs.HOSTED because the app is hosted."""
    from foldmap import app_jobs

    monkeypatch.delenv("FOLDMAP_APP_MAX_RESIDUES", raising=False)
    monkeypatch.setattr(app_jobs, "HOSTED", app_jobs.Limits(max_residues=50))  # in-process, small
    upload = ("ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "chemical/x-cif")
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(upload).run()
    assert at.session_state["stage"] == "edit"  # local: no cap
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(upload).run()
    assert at.session_state["stage"] == "load" and "too large for this app" in at.error[0].value


def test_a_link_to_an_accession_lists_its_structures(monkeypatch):
    import json

    from foldmap import uniprot

    e = uniprot.parse_entry(json.loads((DATA / "uniprot" / "P07911.json").read_text()))
    monkeypatch.setattr(uniprot, "entry", lambda acc, cache_dir=None: e)
    at = app_test()
    at.query_params["id"] = "p07911"
    at.run()
    assert not at.exception and at.session_state["stage"] == "load" and not at.error
    assert at.dataframe[0].value["id"].tolist()[0] == "AF-P07911-F1" and "id" not in at.query_params


# ---------------------------------------------------------------------------------------------- workspace
def test_workspace_top_bar_and_one_view():
    at = opened(DATA / "1LMB.cif")
    assert at.segmented_control(key="pane").options == ["Figure", "Sequence", "3D"]
    assert at.segmented_control(key="pane").value == "Figure" and not at.tabs
    assert "new_structure" in [b.key for b in at.button] and at.session_state["auto_redraw"] is True
    assert any(m.value == "**1LMB**" for m in at.markdown)
    assert any(c.value.startswith("2 chains · 179 residues") for c in at.caption)


def test_only_the_active_view_is_built(monkeypatch):
    import foldmap.interactive as interactive
    import foldmap.seqplot as seqplot
    from app_helpers import count_calls, fresh_caches

    fresh_caches()
    seq = count_calls(monkeypatch, seqplot, "draw_sequence")
    page = count_calls(monkeypatch, interactive, "build_page")
    at = opened(DATA / "1UBQ.cif")
    at.run()  # a rerun in the Figure view
    assert seq == [] and page == []
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not at.exception and len(seq) == 1 and page == []
    at.segmented_control(key="pane").set_value("3D").run()
    assert not at.exception and len(seq) == 1 and len(page) == 1


def test_a_figure_setting_changed_in_another_view_draws_no_figure(monkeypatch):
    from app_helpers import count_calls, fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    draws = count_calls(monkeypatch, app_jobs, "draw_job")
    at = opened(DATA / "1UBQ.cif")
    at.segmented_control(key="pane").set_value("Sequence").run()
    drawn = len(draws)
    at.selectbox(key="palette").select("tol-muted").run()
    assert not at.exception and len(draws) == drawn
    at.segmented_control(key="pane").set_value("Figure").run()
    assert len(draws) == drawn + 1


def test_a_failing_setting_keeps_the_last_figure_and_names_the_control():
    at = opened(DATA / "2LZM.cif")
    before = svg_of(at)
    at.text_input(key="swap").input("α1").run()  # one element: not a swap (Task 9: as_text)
    assert not at.exception
    assert any("swap" in e.value and "Check Swap (Edit layout)." in e.value for e in at.error)
    assert any(_visible(e.value).endswith("got 'α1'. Check Swap (Edit layout).") for e in at.error)  # two sentences
    assert svg_of(at) == before


def test_a_pipeline_crash_is_a_message_not_a_traceback(monkeypatch, caplog):
    """A domain range without res: is a ValueError naming Domains (R2); a crash inside the pipeline (an IndexError
    here, made on purpose) is a message too, logged for the server, and the last figure stays."""
    from app_helpers import fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    at = opened(DATA / "2LZM.cif")
    before = svg_of(at)
    at.text_area(key="domains").set_value("N=A:1-60").run()
    assert not at.exception and svg_of(at) == before
    assert any("Check Domains (Edit layout)." in e.value for e in at.error)
    real = app_jobs.draw_job

    def crash(path, settings_json, *a, **k):
        if "N=res:A:1-60" in settings_json:
            raise IndexError("list index out of range")
        return real(path, settings_json, *a, **k)

    monkeypatch.setattr(app_jobs, "draw_job", crash)
    at.text_area(key="domains").set_value("N=res:A:1-60").run()
    assert not at.exception and svg_of(at) == before
    text = _visible(at.error[0].value)
    assert text.startswith("foldmap could not draw this setting (IndexError: list index out of range).")
    assert "github.com/dzyla/foldmap/issues" in text and text.endswith(" Check Domains (Edit layout).")
    assert "Traceback" not in text and any(r.exc_info and r.name == "foldmap.app" for r in caplog.records)
    assert "](https://github.com/dzyla/foldmap/issues)" in at.error[0].value  # a link, outside the escaped text


def test_a_long_crash_message_keeps_the_report_link_and_the_hint(monkeypatch):
    from app_helpers import fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    real = app_jobs.draw_job

    def crash(path, settings_json, *a, **k):
        if "N=res:A:1-60" in settings_json:
            raise IndexError("x" * 600)
        return real(path, settings_json, *a, **k)

    monkeypatch.setattr(app_jobs, "draw_job", crash)
    at = opened(DATA / "2LZM.cif")
    at.text_area(key="domains").set_value("N=res:A:1-60").run()
    text = at.error[0].value
    assert not at.exception and "x" * 600 not in text
    assert text.endswith(
        " Please report it at [github.com/dzyla/foldmap/issues](https://github.com/dzyla/foldmap/issues)."
        " Check Domains (Edit layout)."
    )


def test_local_mode_shows_the_visitors_own_path_in_full(monkeypatch, tmp_path):
    """Locally the paths are the visitor's own (load, and a figure that fails); hosted they are cut to a name."""
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "1")
    at = load(app_test(), "/nowhere/missing.cif")
    assert not at.exception and "/nowhere/missing.cif" in _visible(at.error[0].value)
    junk = tmp_path / "family.fasta"
    junk.write_text("not an alignment\n", encoding="utf-8")
    at = opened(DATA / "1UBQ.cif")
    at.session_state["msa_path"] = str(junk)
    at.run()
    assert not at.exception and _visible(at.error[0].value).startswith(f"{junk}: unknown alignment format")


def test_a_slow_structure_turns_auto_redraw_off(monkeypatch):
    from app_helpers import fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    monkeypatch.setattr(app_jobs, "SLOW_SECONDS", 0.0)
    at = opened(DATA / "1UBQ.cif")
    assert at.session_state["auto_redraw"] is False and "redraw" in [b.key for b in at.button]
    assert any(c.value.startswith("Auto-redraw is off for this structure") for c in at.caption)
    before = svg_of(at)
    at.selectbox(key="palette").select("tol-bright").run()
    assert svg_of(at) == before  # not redrawn until asked
    at.button(key="redraw").click().run()
    assert svg_of(at) != before


def test_a_slow_setting_turns_auto_redraw_off_after_a_fast_first_draw(monkeypatch):
    from app_helpers import fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    real = app_jobs.draw_job

    def slow_for_tol_muted(path, settings_json, *a, **k):
        out = real(path, settings_json, *a, **k)
        return {**out, "seconds": 9.5} if "tol-muted" in settings_json else out

    monkeypatch.setattr(app_jobs, "draw_job", slow_for_tol_muted)
    at = opened(DATA / "1UBQ.cif")
    assert at.session_state["auto_redraw"] is True and "redraw" not in [b.key for b in at.button]
    first = svg_of(at)
    at.selectbox(key="palette").select("tol-muted").run()
    slow = svg_of(at)
    assert slow != first and at.session_state["auto_redraw"] is False  # the slow figure is shown, then auto stops
    assert any(c.value.startswith("Auto-redraw is off for this structure: drawing took 9.5 s") for c in at.caption)
    at.selectbox(key="palette").select("tol-bright").run()
    assert svg_of(at) == slow
    at.toggle(key="auto_redraw").set_value(True).run()  # back on: stays on for this structure
    assert at.session_state["auto_redraw"] is True and svg_of(at) != slow


def test_figure_card_is_accessible_and_scrolls_only_when_very_wide():
    from foldmap.app_views import figure_card

    card = figure_card('<svg width="400pt" height="300pt" viewBox="0 0 400 300"></svg>', "#ffffff", "1UBQ")
    assert 'role="img"' in card and 'aria-label="Topology diagram of 1UBQ"' in card and "max-width:none" not in card
    wide = figure_card('<svg width="3000pt" height="500pt" viewBox="0 0 3000 500"></svg>', "#ffffff", "x")
    assert "max-width:none" in wide and 'class="fm-card"' in wide
    assert "height:min(60vh,667px)" in wide  # 500 pt is 667 CSS px
    assert "height:min(60vh,500px)" in figure_card('<svg width="3000px" height="500px"></svg>', "#fff", "x")
    assert "background:#fff&quot;&gt;" in figure_card("<svg></svg>", '#fff">', "x")
    seq = figure_card('<svg width="300pt" height="900pt"></svg>', "#fff", "1UBQ", kind="Sequence", natural=True)
    assert 'aria-label="Sequence of 1UBQ"' in seq and "width:400px;max-width:none;height:1200px;max-height:none" in seq


def test_a_layout_file_forgets_the_control_changed_last():
    at = opened(DATA / "1UBQ.cif")
    at.text_input(key="swap").input("α1").run()
    assert at.session_state["last_changed"] == "swap"
    up = ("layout.yaml", b"foldmap: 1\nlayout: {rotate: 30}\n", "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    assert not at.exception and "last_changed" not in at.session_state


# ---------------------------------------------------------------------------------------------- views and downloads
def test_sequence_view_offers_every_protein_chain():
    at = opened(DATA / "1LMB.cif")
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert at.slider(key="seq_columns").value == 60 and at.selectbox(key="seq_chain").options == ["3", "4"]
    assert at.segmented_control(key="seq_ss_colour").value == "figure" and at.toggle(key="seq_full").value is False
    first = svg_of(at)
    assert 'aria-label="Sequence of 1LMB"' in first
    at.selectbox(key="seq_chain").select("4").run()
    assert not at.exception and svg_of(at) != first


def test_a_single_chain_has_no_chain_picker():
    at = opened(DATA / "1UBQ.cif")
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not at.exception and "seq_chain" not in [s.key for s in at.selectbox]


def test_sequence_options_do_not_redraw_the_figure(monkeypatch):
    import foldmap.spec as spec
    from app_helpers import count_calls

    at = opened(DATA / "1UBQ.cif")
    at.segmented_control(key="pane").set_value("Sequence").run()
    calls = count_calls(monkeypatch, spec, "draw_spec")
    first = svg_of(at)
    at.slider(key="seq_columns").set_value(80).run()
    assert calls == [] and not at.exception and svg_of(at) != first


def test_the_3d_page_is_embedded_unless_too_large(monkeypatch):
    from foldmap import app_jobs

    at = opened(DATA / "1UBQ.cif")
    at.segmented_control(key="pane").set_value("3D").run()
    assert not at.exception and not at.info
    monkeypatch.setattr(app_jobs, "PAGE_EMBED_MAX", 1000)
    at.run()
    assert any("too large to show here" in i.value for i in at.info)
    assert not any("rcsb.org" in i.value for i in at.info)  # a file: no PDB id to point at


def test_a_too_large_3d_page_of_a_pdb_entry_points_at_rcsb(monkeypatch, tmp_path):
    from foldmap import app_jobs

    seed_cache(monkeypatch, tmp_path, "1UBQ")
    monkeypatch.setattr(app_jobs, "PAGE_EMBED_MAX", 1000)
    at = opened("1UBQ")
    at.segmented_control(key="pane").set_value("3D").run()
    text = _visible(at.info[0].value)
    assert text.startswith("The 3D page is ") and "MB, too large to show here. Download it (Download ▸ 3D page)" in text
    assert text.endswith(" or open 1UBQ at https://www.rcsb.org/3d-view/1UBQ.")


def test_sequence_and_3d_limits_are_remembered_apart_from_the_figure(monkeypatch):
    """A sequence over a limit is remembered for its own view (R4): the figure still draws, other options still try."""
    from app_helpers import fresh_caches
    from foldmap import app_jobs

    fresh_caches()

    def too_long(*a, **k):
        raise app_jobs.JobLimit("time", "stopped after 90 s")

    at = opened(DATA / "1UBQ.cif")
    monkeypatch.setattr(app_jobs, "sequence_job", too_long)
    monkeypatch.setattr(app_jobs, "page_job", too_long)
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not at.exception and _visible(at.error[0].value).startswith("1UBQ took longer than")
    jobs = sorted(k[1] for k in app_jobs.FAILED)
    assert jobs == ["sequence"] and next(iter(app_jobs.FAILED))[3]  # keyed with the view
    at.segmented_control(key="pane").set_value("3D").run()
    assert not at.exception and _visible(at.error[0].value).startswith("1UBQ took longer than")
    assert sorted(k[1] for k in app_jobs.FAILED) == ["page", "sequence"]
    assert "dl_seq_svg" not in [b.key for b in at.download_button]  # nothing drawn: nothing to download
    assert any(c.value == "Open the Sequence view to make this file." for c in at.caption)
    at.segmented_control(key="pane").set_value("Figure").run()
    assert not at.exception and not at.error and "<svg" in svg_of(at)


def test_downloads_are_built_only_when_asked(monkeypatch):
    """Each group offers what its view drew (nothing before the view was opened), as files made on a click only."""
    from app_helpers import count_calls, fresh_caches
    from foldmap import app_jobs

    fresh_caches()
    calls = count_calls(monkeypatch, app_jobs, "files_job")
    seq = count_calls(monkeypatch, app_jobs, "sequence_job")
    page = count_calls(monkeypatch, app_jobs, "page_job")
    at = opened(DATA / "1UBQ.cif")
    figure = {"dl_svg", "dl_png", "dl_pdf", "dl_layout"}
    assert {b.key for b in at.download_button} == figure
    assert any(c.value == "Open the Sequence view to make this file." for c in at.caption)
    assert any(c.value == "Open the 3D view to make this file." for c in at.caption)
    assert calls == [] and seq == [] and page == []
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert len(seq) == 1  # the view's own figure (svg)
    assert {b.key for b in at.download_button} == figure | {"dl_seq_svg", "dl_seq_png", "dl_seq_pdf"}
    at.segmented_control(key="pane").set_value("3D").run()
    assert len(page) == 1 and {b.key for b in at.download_button} == figure | {
        "dl_seq_svg",
        "dl_seq_png",
        "dl_seq_pdf",
        "dl_html",
    }
    at.segmented_control(key="pane").set_value("Figure").run()
    assert calls == [] and len(seq) == 1 and len(page) == 1
    assert any(c.value.startswith("SVG: edit in Illustrator or Inkscape · PNG: slides and web") for c in at.caption)


def test_downloads_offer_what_was_drawn():
    """A figure setting changed in another view: the Figure files stay those of the figure drawn, and say so; the
    sequence files are those of the sequence drawn, with its options."""
    at = opened(DATA / "1UBQ.cif")
    drawn = at.session_state["last_good"]["blob"]
    at.segmented_control(key="pane").set_value("Sequence").run()
    at.slider(key="seq_columns").set_value(80).run()
    blob, view_json = at.session_state["seq_drawn"]
    assert '"columns": 80' in view_json and blob == drawn
    at.selectbox(key="palette").select("tol-muted").run()
    ss = at.session_state
    assert ss["last_good"]["blob"] == drawn and "tol-muted" in ss["seq_drawn"][0]
    assert any(c.value.startswith("These are the files of the figure as last drawn") for c in at.caption)
    at.segmented_control(key="pane").set_value("Figure").run()
    assert not any(c.value.startswith("These are the files of the figure as last") for c in at.caption)


def test_every_download_makes_its_whole_file():
    """The popover's callables (R8), called as a click calls them: each makes its whole file."""
    import json

    from foldmap import app_jobs
    from foldmap.app_views import figure_downloads, page_downloads, sequence_downloads

    path, blob, lim = str(DATA / "1UBQ.cif"), json.dumps({"theme": "publication"}), app_jobs.Limits()
    built = []  # the 3D page: the one the view built (its cache), not built again

    def page(p, b):
        built.append((p, b))
        return app_jobs.page_job(p, b)

    view = json.dumps({"chain": "A", "columns": 30, "full_sequence": False, "ss_colour": "black"}, sort_keys=True)
    files = [
        *figure_downloads(path, blob, "1UBQ", lim),
        *sequence_downloads(path, blob, view, "1UBQ", lim),
        *page_downloads(path, blob, "1UBQ", page),
    ]
    names = {d.key: d.file_name for d in files}
    assert names == {
        "dl_svg": "1UBQ.svg",
        "dl_png": "1UBQ.png",
        "dl_pdf": "1UBQ.pdf",
        "dl_layout": "1UBQ-layout.yaml",
        "dl_seq_svg": "1UBQ-sequence.svg",
        "dl_seq_png": "1UBQ-sequence.png",
        "dl_seq_pdf": "1UBQ-sequence.pdf",
        "dl_html": "1UBQ-explorer.html",
    }
    start = {"svg": b"<svg", "png": b"\x89PNG", "pdf": b"%PDF", "yaml": b"foldmap: 1\n", "html": b"<!doctype html>"}
    for d in files:
        data = d.data()
        assert data.startswith(start[d.file_name.rsplit(".", 1)[1]]), d.key
    assert built == [(path, blob)]
    assert b"columns: 30" not in files[4].data()  # the sequence file is the sequence figure, not the layout
    eager = (("pdf", b"%PDF eager"), ("png", b"\x89PNG eager"))
    kept = {d.key: d.data() for d in figure_downloads(path, blob, "x", lim, eager)}
    assert kept["dl_png"] == b"\x89PNG eager" and kept["dl_pdf"] == b"%PDF eager" and kept["dl_svg"][:4] == b"<svg"


def test_download_buttons_are_offered():
    from foldmap.app import figure_files

    files = figure_files(DATA / "1UBQ.cif", {"theme": "publication"})
    assert set(files) == {"svg", "png", "pdf"}
    assert files["png"][:4] == b"\x89PNG" and files["pdf"][:4] == b"%PDF" and b"<svg" in files["svg"]


def test_download_files_are_whole():
    import json

    from foldmap import app_jobs
    from foldmap.app_views import layout_yaml

    blob = json.dumps({"theme": "publication"})
    png = app_jobs.figure_file(str(DATA / "1UBQ.cif"), blob, "png", app_jobs.Limits())
    assert png[:4] == b"\x89PNG"
    assert layout_yaml(str(DATA / "1UBQ.cif"), blob).startswith(b"foldmap: 1\n")
    html = app_jobs.page_file(str(DATA / "1UBQ.cif"), blob, app_jobs.Limits())
    assert html.startswith(b"<!doctype html>")


def test_a_failing_setting_offers_the_last_good_figure():
    """The canvas keeps the last good figure, and so does the Download popover, which says so."""
    from foldmap.app_views import figure_offered

    at = opened(DATA / "2LZM.cif")
    good = at.session_state["last_good"]["blob"]
    at.text_input(key="swap").input("α1").run()  # one element: not a swap
    assert not at.exception and at.error
    assert {"dl_svg", "dl_png", "dl_pdf", "dl_layout"} <= {b.key for b in at.download_button}
    assert any(c.value.startswith("These are the previous figure's files") for c in at.caption)
    ss = at.session_state
    blob, eager, note = figure_offered(ss, ss["path"], ss["figure_failed"])
    assert blob == good and note.startswith("These are the previous figure's files") and eager is None
    at.text_input(key="swap").input("").run()
    assert not at.error and not any(c.value.startswith("These are the previous") for c in at.caption)


def test_figure_offered_is_the_figure_on_screen():
    from foldmap.app_views import figure_offered

    files = (("pdf", b"%PDF"), ("png", b"\x89PNG"))
    ss = {"auto_redraw": True, "last_good": {"path": "p", "blob": "B1", "svg": "<svg/>", "files": files}}
    assert figure_offered(ss, "p", "B1") == ("B1", files, "")  # drawn: its eager files come along
    blob, eager, note = figure_offered(ss, "p", "B2")  # changed in another view, not drawn yet: the drawn one
    assert (blob, eager) == ("B1", files) and note.startswith("These are the files of the figure as last drawn")
    blob, eager, note = figure_offered({**ss, "figure_failed": "B2"}, "p", "B2")
    assert (blob, eager) == ("B1", files) and note.startswith("These are the previous figure's files")
    assert figure_offered({**ss, "figure_failed": "B2"}, "q", "B2") is None  # nothing good of this structure
    assert figure_offered({**ss, "auto_redraw": False, "drawn_blob": "B1"}, "p", "B2") == ("B1", files, "")


def test_share_url_carries_id_and_theme_only():
    from foldmap.app_views import share_url

    assert share_url("https://foldmap.streamlit.app/?id=X", "1LMB", "journal") == (
        "https://foldmap.streamlit.app/?id=1LMB&theme=journal"
    )
    assert share_url(None, "1LMB", "trace") == "?id=1LMB&theme=trace"


def test_the_popover_shows_the_share_link_of_an_id_only(monkeypatch, tmp_path):
    seed_cache(monkeypatch, tmp_path, "1UBQ")
    at = opened("1UBQ")
    assert [c.value for c in at.code] == ["?id=1UBQ&theme=publication"]
    assert any(c.value == "Link to this figure" for c in at.caption)
    at = opened(DATA / "1UBQ.cif")
    assert not at.code
