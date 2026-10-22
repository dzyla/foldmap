"""Uploads whose names Windows reads differently: Streamlit's own extension check there splits a name at its first
colon (an NTFS alternate data stream), so a crafted name such as "[n](https://evil.example).yaml" failed the check
and crashed the page. The app checks every upload's suffix itself, and keeps uploads by content hash and a
whitelisted suffix only."""

import os
from pathlib import Path

import pytest

pytest.importorskip("streamlit")

from app_helpers import DATA, app_test, opened, rendered_text, seed_cache, svg_of, visible  # noqa: E402


@pytest.fixture
def windows_names(monkeypatch):
    """Streamlit's filename check as it runs on Windows (file_uploader_utils: an alternate data stream split)."""
    from streamlit.elements.lib import file_uploader_utils as utils

    def as_on_windows(filename: str) -> tuple[str, str]:
        if ":" in filename:
            main, ads = filename.split(":", 1)
            if os.path.splitext(ads)[1]:
                return main, os.path.splitext(main)[1]
        return filename, os.path.splitext(filename)[1]

    monkeypatch.setattr(utils, "_get_main_filename_and_extension", as_on_windows)


CRAFTED = "[n](https://evil.example)"


def test_a_crafted_structure_upload_is_drawn_or_refused_with_a_message(windows_names, monkeypatch):
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    for data in (b"garbage\n", (DATA / "1UBQ.cif").read_bytes()):
        at = app_test()
        at.run()
        at.file_uploader(key="upload").set_value((f"{CRAFTED}.cif", data, "chemical/x-cif")).run()
        assert not at.exception and (at.error or at.session_state["stage"] == "edit")
        assert "](" not in rendered_text(at)


def test_crafted_names_layout_and_alignment_uploads_do_not_crash(windows_names, monkeypatch, tmp_path):
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    seed_cache(monkeypatch, tmp_path, "1UBQ")
    at = opened("1UBQ")
    at.file_uploader(key="names_upload").set_value((f"{CRAFTED}.yaml", b"elements:\n  res:A:23-34: Hx\n", "x")).run()
    assert not at.exception and len(at.session_state["names"]) == 1
    fasta = (DATA / "ubq_family.fasta").read_bytes()
    at.file_uploader(key="msa_upload").set_value((f"{CRAFTED}.fasta", fasta, "text/plain")).run()
    assert not at.exception and at.session_state["msa_path"]
    kept = (at.session_state["names"][0], at.session_state["msa_path"])
    at.file_uploader(key="layout_upload").set_value((f"{CRAFTED}.yaml", b"theme: journal\n", "x")).run()
    assert not at.exception and at.session_state["theme"] == "journal"  # its settings replace the session's
    for path in kept:
        name = os.path.basename(path)
        assert name.startswith("upload-") and "evil" not in name and "[" not in name


def test_an_upload_of_the_wrong_kind_is_a_message_not_a_crash(monkeypatch, tmp_path):
    """Streamlit no longer refuses a suffix itself (the uploaders name a MIME type too): the app does, in words."""
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("evil.sh", b"data_x\n", "text/x-sh")).run()
    assert not at.exception and at.session_state["stage"] == "load"
    said = [visible(e.value) for e in at.error]
    assert said == ["evil.sh is not a structure file: upload a .cif, .mmcif, .pdb or .ent file."]
    seed_cache(monkeypatch, tmp_path, "1UBQ")
    at = opened("1UBQ")
    at.file_uploader(key="names_upload").set_value(("names.sh", b"elements: {}\n", "x")).run()
    assert (
        not at.exception
        and at.session_state["names_error"] == "names.sh is not a names file: upload a .yaml or .yml file."
    )
    at.file_uploader(key="msa_upload").set_value(("aln.sh", b">x\nAA\n", "x")).run()
    assert not at.exception and at.session_state["msa_error"].startswith("aln.sh is not an alignment: upload a .fasta")
    at.file_uploader(key="layout_upload").set_value(("layout.sh", b"theme: journal\n", "x")).run()
    assert not at.exception and at.session_state["theme"] == "publication"
    assert at.session_state["layout_error"] == "layout.sh is not a layout file: upload a .yaml or .yml file."


def test_upload_suffixes_are_read_the_same_on_every_system():
    from foldmap.app_edits import upload_suffix

    assert upload_suffix("a.CIF") == ".cif" and upload_suffix("x:y.yaml") == ".yaml"
    assert upload_suffix("[n](https://evil.example).yaml") == ".yaml"
    assert upload_suffix("dir\\x.fasta") == ".fasta" and upload_suffix("noext") == "" and upload_suffix("a.") == ""


def test_kept_uploads_are_named_by_content_and_a_whitelisted_suffix(monkeypatch, tmp_path):
    from foldmap import fetch

    monkeypatch.setattr(fetch, "CACHE", tmp_path)
    for name in ("x:y.cif", "a\\b.pdb", "C:evil.cif", "../../e.ent", "x\n.cif", "[n](https://e.x).cif"):
        kept = fetch.keep_upload(name, b"data_x\n")
        assert kept.parent == tmp_path / "uploads" and kept.name[:7] == "upload-" and len(kept.name) <= 23
        assert kept.suffix in fetch.STRUCTURE_SUFFIXES


DDG = b"chain,residue,ddG,site\nA,3,0.5,core\nA,30,-1.2,rim\nA,999,1,core\n"
XL = b"chain1,residue1,chain2,residue2,distance\nA,6,A,48,\nA,1,A,72,5\n"


def test_a_data_upload_draws_what_the_command_line_draws(monkeypatch, tmp_path, capsys):
    """Your data: the per-residue values and the crosslinks uploaded draw the same SVG as `foldmap plot --data
    --crosslinks` on those files; Colour by offers data only once data are loaded (and turns to it, as the command
    line does); the rows not used and the crosslink longer than its span are said under the figure, each table named
    as uploaded; the share link holds the id and theme only."""
    from foldmap.app_views import figure_card
    from foldmap.cli import main

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "1")
    seed_cache(monkeypatch, tmp_path, "1UBQ")
    at = opened("1UBQ")
    assert "data" not in at.selectbox(key="color_by").options
    at.file_uploader(key="data_upload").set_value(("my ddG.csv", DDG, "text/csv")).run()
    at.file_uploader(key="crosslinks_upload").set_value(("xl.csv", XL, "text/csv")).run()
    ss = at.session_state
    assert not at.exception and ss["color_by"] == "data" and "data" in at.selectbox(key="color_by").options
    assert ss["data_read"] == "column ddG (numbers) on 2 residues; 1 row not used (named under the figure)"
    assert ss["crosslinks_read"] == "2 crosslinks" and at.selectbox(key="data_column").options == [
        "first column of numbers",
        "ddG",
        "site",
    ]
    out = tmp_path / "cli.svg"
    args = ["plot", str(tmp_path / "1UBQ.cif"), "-o", str(out), "--data", ss["data_path"]]
    assert main([*args, "--crosslinks", ss["crosslinks_path"]]) == 0
    capsys.readouterr()
    cli = out.read_text(encoding="utf-8")
    assert svg_of(at).count('id="ia:xl:') == 2 and svg_of(at) == figure_card(
        cli[cli.index("<svg") :], "#ffffff", "1UBQ"
    )
    said = visible(next(w.value for w in at.warning if "Your data" in w.value))
    assert "--data my ddG.csv row 4: chain A has no residue 999" in said and "longer than its 5 Å span" in said
    assert ss["data_path"] not in said and set(at.query_params) == {"id", "theme"}  # uploads never enter the link
    at.selectbox(key="data_column").select("site").run()
    assert not at.exception and 'id="data:A:3"' in svg_of(at)  # categories: marks
    at.file_uploader(key="data_upload").set_value(None).run()
    assert not at.exception and ss["color_by"] == "chain" and "data" not in at.selectbox(key="color_by").options


def test_a_data_upload_that_names_no_residue_is_refused_in_words(monkeypatch, tmp_path):
    """A table that names no residue is refused, in words. Section review 13b-13d, Important 7: the table and the
    structure are read in a job under the session's limits (app_jobs.table_job), as every structure read."""
    from foldmap import app_jobs

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    seed_cache(monkeypatch, tmp_path, "1UBQ")
    jobs, real = [], app_jobs.run_job
    monkeypatch.setattr(app_jobs, "run_job", lambda fn, *a, **k: (jobs.append((fn, k["limits"])), real(fn, *a, **k))[1])
    at = opened("1UBQ")
    at.file_uploader(key="data_upload").set_value(("bad.csv", b"chain,residue,x\nA,900,1\n", "text/csv")).run()
    ss = at.session_state
    assert (app_jobs.table_job, app_jobs.limits(hosted=True)) in jobs
    assert not at.exception and ss["data_path"] == "" and ss["color_by"] == "chain"
    assert ss["data_error"].startswith("bad.csv: no row names a residue of the structure")
    assert "data" not in at.selectbox(key="color_by").options
    at.file_uploader(key="crosslinks_upload").set_value(("xl.sh", XL, "text/csv")).run()
    assert ss["crosslinks_error"] == "xl.sh is not a table: upload a .csv, .tsv or .txt file."


def test_a_new_structure_after_data_draws_in_its_chain_colours(monkeypatch, tmp_path, capsys):
    """Section review 13b-13d, Important 3 and Minor 9: data loaded turn Colour by to data; a new structure forgets
    the data, and so its colouring by them (else every element of the next figure was grey), and the last
    structure's crosslinks error. The next figure is the command line's, in chain colours."""
    from foldmap import app_jobs
    from foldmap.app_views import figure_card
    from foldmap.cli import main

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "1")
    monkeypatch.setattr(app_jobs, "SLOW_SECONDS", float("inf"))
    seed_cache(monkeypatch, tmp_path, "1UBQ", "2LZM")
    at = opened("1UBQ")
    at.file_uploader(key="data_upload").set_value(("ddG.csv", DDG, "text/csv")).run()
    at.file_uploader(key="crosslinks_upload").set_value(("xl.sh", XL, "text/csv")).run()
    ss = at.session_state
    assert ss["color_by"] == "data" and ss["crosslinks_error"]
    at.button(key="new_structure").click().run()
    at.text_input(key="source").input("2LZM")
    at.button(key="load").click().run()
    assert not at.exception and ss["stage"] == "edit" and ss["name"] == "2LZM", [e.value for e in at.error]
    assert ss["color_by"] == "chain" and ss["data_path"] == "" and ss["crosslinks_error"] == ""
    assert "data" not in at.selectbox(key="color_by").options
    out = tmp_path / "2LZM.svg"
    assert main(["plot", str(tmp_path / "2LZM.cif"), "-o", str(out)]) == 0
    capsys.readouterr()
    cli = out.read_text(encoding="utf-8")
    assert svg_of(at) == figure_card(cli[cli.index("<svg") :], "#ffffff", "2LZM")


def test_colour_by_data_without_data_is_the_themes_colouring():
    """Section review 13b-13d, Important 3: the app's FigureSpec never draws color_by=data with no data loaded."""
    from foldmap.app import look_of

    assert look_of({"theme": "publication", "color_by": "data"}).color_by == "chain"
    assert look_of({"theme": "publication", "color_by": "data", "data_path": "d.csv"}).color_by == "data"
    assert look_of({"theme": "publication", "color_by": "sstype"}).color_by == "sstype"


def test_the_sequence_view_draws_the_data_too(monkeypatch, tmp_path):
    """Section review 13b-13d, Important 4: the app's Sequence view under color_by=data shows the data's stripes, as
    `foldmap sequence --data` does, not every element grey."""
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "1")
    seed_cache(monkeypatch, tmp_path, "1UBQ")
    at = opened("1UBQ")
    at.segmented_control(key="pane").set_value("Sequence").run()
    assert not any("ss-stripe:" in m.value for m in at.markdown)
    at.file_uploader(key="data_upload").set_value(("ddG.csv", DDG, "text/csv")).run()
    assert not at.exception and at.session_state["color_by"] == "data"
    assert any("<svg" in m.value and 'id="ss-stripe:' in m.value for m in at.markdown)


def test_an_upload_clears_only_its_own_tables_warning(monkeypatch, tmp_path):
    """Section review 13b-13d, Minor 10: a layout file's data and crosslinks tables not on this machine are each
    warned of until a table of that kind is uploaded."""
    from types import SimpleNamespace

    from foldmap import fetch
    from foldmap.app import settings_from_document
    from foldmap.app_controls import missing_tables, store_table

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "1")
    monkeypatch.setattr(fetch, "CACHE", tmp_path)
    given = {"data": str(tmp_path / "gone" / "d.csv"), "crosslinks": str(tmp_path / "gone" / "xl.csv")}
    settings = settings_from_document({"layout": given})
    assert settings["data_missing"] == given and settings["data_path"] == settings["crosslinks_path"] == ""
    ss = {
        **settings,
        "path": str(DATA / "1UBQ.cif"),
        "data_upload": SimpleNamespace(name="d.csv", getvalue=lambda: DDG),
    }
    store_table(ss, "data")
    assert ss["data_path"] and not ss["data_error"] and ss["data_missing"] == {"crosslinks": given["crosslinks"]}
    assert len(missing_tables(ss)) == 1 and "xl.csv" in visible(missing_tables(ss)[0])
    ss["crosslinks_upload"] = SimpleNamespace(name="xl.csv", getvalue=lambda: XL)
    store_table(ss, "crosslinks")
    assert ss["crosslinks_path"] and ss["data_missing"] == {} and missing_tables(ss) == []


def test_an_uploaded_structures_data_name_it_as_uploaded(monkeypatch, tmp_path):
    """Section review 13b-13d, Minor 21: the data downloads of an uploaded structure name it as uploaded, not by the
    copy the server keeps (upload-<hash>.cif)."""
    import json

    from foldmap import app_jobs, fetch

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    monkeypatch.setattr(fetch, "CACHE", tmp_path)
    at = app_test()
    at.run()
    at.file_uploader(key="upload").set_value(("my ubq.cif", (DATA / "1UBQ.cif").read_bytes(), "x")).run()
    ss = at.session_state
    assert not at.exception and ss["stage"] == "edit" and ss["structure_uploaded"] == "my ubq.cif"
    assert Path(ss["path"]).name.startswith("upload-")
    blob = json.dumps({"theme": "publication"})
    named = json.loads(app_jobs.data_job(ss["path"], blob, "json", ss["structure_uploaded"]))
    assert named["structure"]["name"] == "my ubq.cif"
    assert json.loads(app_jobs.data_job(ss["path"], blob, "json"))["structure"]["name"] == Path(ss["path"]).name
