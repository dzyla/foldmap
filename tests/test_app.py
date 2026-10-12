"""The Streamlit app, driven headlessly: Start, load (drawn at once), share links, edit, download."""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from app_helpers import DATA, app_test, load, opened, seed_cache, svg_of  # noqa: E402


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
    """Was test_picking_a_style_opens_the_editor_with_that_theme (the dialog arrives in Task 6)."""
    at = opened(DATA / "1LMB.cif")
    at.selectbox(key="theme").select("trace").run()
    assert not at.exception and at.session_state["loop_arrows"] is True and "<svg" in svg_of(at)


def test_controls_change_the_figure():
    at = load(app_test(), DATA / "1UBQ.cif")
    assert "resnum:" not in svg_of(at)
    at.checkbox(key="residue_numbers").check()
    at.run()
    assert "resnum:" in svg_of(at)
    at.selectbox(key="mode").select("stack")
    at.run()
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
    at = app_test()
    at.run()
    at.text_input(key="source").set_value(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
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
    assert text == "notes: no protein residues with a complete backbone" and "uploads" not in text


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


def test_app_sequence_tab_with_the_conservation_theme(tmp_path):
    """The conservation theme colours the sequence view's elements from the same alignment as the figure."""
    at = load(app_test(), DATA / "1UBQ.cif")
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    at.selectbox(key="theme").select("conservation").run()
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


def test_app_pruned_upload_in_the_sequence_tab_is_a_message(monkeypatch):
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "chemical/x-cif")).run()
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
    assert text.startswith("1UBQ is too large for this app: 76 residues in 1 chain(s)")
    assert "foldmap plot 1UBQ -o 1UBQ.svg" in text


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
    at.session_state["seq_chain"] = "B"  # Sequence-view and names widgets of this structure (Tasks 5 and 10)
    at.session_state["names__kw"] = ["beta"]
    at.button(key="new_structure").click().run()
    assert at.session_state["stage"] == "load" and "id" not in at.query_params
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    ss = at.session_state
    assert not at.exception and ss["stage"] == "edit" and ss["name"] == "1UBQ"
    assert ss["theme"] == "trace" and ss["font_scale"] == 1.4 and ss["mark"] == ""
    assert "seq_chain" not in ss and "names__kw" not in ss and ss["auto_redraw"] is True


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
