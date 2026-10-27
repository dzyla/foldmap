"""Plan 29's Review focus: five inputs a visitor meets that no fix task's own tests exercise, each driven through
the app as the visitor drives it (type, press Draw, or drop a file) and judged on what the page shows. No network:
downloads are answered by app_helpers.serve (a 404, as the PDB answers an obsolete or mistyped id, for anything
else)."""

import gzip
import sys

import pytest

pytest.importorskip("streamlit")
from app_helpers import (  # noqa: E402
    DATA,
    app_test,
    count_calls,
    fresh_caches,
    keep_downloads,
    load,
    opened,
    press,
    seed_cache,
    serve,
    svg_of,
)
from app_helpers import visible as _visible  # noqa: E402

needs_fork = pytest.mark.skipif(sys.platform == "win32", reason="isolated jobs need forkserver")
HOSTED = "0"  # FOLDMAP_APP_LOCAL_FILES: the public app (foldmap.streamlit.app) is hosted


@pytest.fixture(autouse=True)
def _fresh():
    """No figure, size or failure kept from another test (the floor's test draws a page held small)."""
    fresh_caches()
    yield
    fresh_caches()


def _upload(at, name: str, data: bytes):
    at.run()
    at.file_uploader(key="upload").set_value((name, data, "chemical/x-cif")).run()
    return at


# 1. an obsolete or mistyped PDB id ---------------------------------------------------------------------------------
@pytest.mark.parametrize("hosted", ["0", "1"])
@pytest.mark.parametrize("typed", ["XXXX", "1ABC"])  # XXXX: no id at all; 1ABC: an id the PDB answers with a 404
def test_an_unknown_id_is_said_and_the_box_kept(monkeypatch, hosted, typed):
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", hosted)
    asked = serve(monkeypatch, {})
    at = load(app_test(), typed)
    assert not at.exception
    said = [_visible(e.value) for e in at.error]
    assert len(said) == 1 and "not found" in said[0].lower() and typed in said[0], said
    assert at.session_state["stage"] != "edit" and at.text_input(key="source").value == typed
    assert asked == ([] if typed == "XXXX" else ["https://files.rcsb.org/download/1ABC.cif"])


# 2. an entry over the hosted size limit ----------------------------------------------------------------------------
@needs_fork
def test_an_entry_over_the_hosted_limit_is_refused_at_once_with_its_count_and_the_limit(monkeypatch, tmp_path):
    from foldmap import app_jobs

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", HOSTED)
    monkeypatch.setenv("FOLDMAP_APP_MAX_RESIDUES", "50")  # app_jobs.limits(hosted=True), its cap held low
    assert app_jobs.limits(hosted=True).max_residues == 50 and app_jobs.limits(hosted=True).isolated
    seed_cache(monkeypatch, tmp_path, "1UBQ")  # 76 residues
    draws = count_calls(monkeypatch, app_jobs, "draw_job")
    at = load(app_test(), "1UBQ")
    assert not at.exception and at.session_state["stage"] == "load" and draws == []  # no draw started
    said = [_visible(e.value) for e in at.error]
    assert len(said) == 1 and said[0].startswith("1UBQ is too large for this app: 76 residues in 1 chain(s)")
    assert "the limit here is 50" in said[0]


# 3. a file with no protein and no nucleic acid, or a pure nucleic-acid file ----------------------------------------
@needs_fork
def test_a_file_of_waters_only_is_one_plain_sentence(monkeypatch):
    """1UBQ without its protein: its 58 waters (HETATM), its entity records kept, as a ligand-only file has them."""
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", HOSTED)
    lines = (DATA / "1UBQ.cif").read_text(encoding="utf-8").splitlines(keepends=True)
    waters = "".join(line for line in lines if not line.startswith("ATOM ")).encode()
    at = _upload(app_test(), "waters.cif", waters)
    assert not at.exception and at.session_state["stage"] == "load"
    assert [_visible(e.value) for e in at.error] == ["waters: no structure in this file (no protein, DNA or RNA chain)"]


@needs_fork
def test_a_pure_rna_file_is_drawn(monkeypatch):
    """1EHZ, tRNA-Phe: one RNA chain, no protein."""
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", HOSTED)
    at = _upload(app_test(), "1EHZ.cif", gzip.decompress((DATA / "na" / "1EHZ.cif.gz").read_bytes()))
    assert not at.exception and at.session_state["stage"] == "edit", [e.value for e in at.error]
    assert not at.error and "<svg" in svg_of(at) and at.session_state["protein_chains"] == []
    assert at.session_state["na_chains"] == ["A"] and "nucleotides" in at.session_state["info"]


# 4. an AlphaFold id typed into the same box ------------------------------------------------------------------------
@needs_fork
def test_an_alphafold_id_is_fetched_and_drawn_like_a_pdb_id(monkeypatch):
    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", HOSTED)
    asked = serve(monkeypatch, {"/AF-P04637-F1-model_v6.cif": (DATA / "AF-P04637.cif").read_bytes()})
    at = load(app_test(), "AF-P04637-F1")
    assert not at.exception and at.session_state["stage"] == "edit", [e.value for e in at.error]
    assert asked == ["https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.cif"]
    assert at.session_state["name"] == "AF-P04637-F1" and "<svg" in svg_of(at)
    assert at.query_params["id"] == "AF-P04637-F1"  # its share link, as a PDB id's


# 5. downloads of a figure below the readability floor --------------------------------------------------------------
def test_every_download_of_a_figure_below_the_floor_is_whole_and_the_note_is_shown(monkeypatch):
    """The floor's fixture (tests/test_floor.py small_page): pages held at 3 in, so 1UBQ's text falls under 6 pt.
    Local, so the draw and the downloads run in this process, where the page is held small."""
    from foldmap import layout

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "1")
    monkeypatch.setattr(layout, "MAX_PAGE_IN", 3.0)
    kept = keep_downloads(monkeypatch)
    at = opened(DATA / "1UBQ.cif")
    notes = [_visible(w.value) for w in at.warning]
    assert any(n.startswith("Below the readability floor: its text is") for n in notes), notes
    start = {"svg": b"<svg", "png": b"\x89PNG", "pdf": b"%PDF"}
    for ext, magic in start.items():
        data = press(at, f"dl_{ext}", kept)
        assert data.startswith(magic) and len(data) > 1000, ext
    for key in ("dl_layout", "dl_data_json", "dl_data_zip"):
        assert press(at, key, kept), key
    assert not at.exception and not at.error
