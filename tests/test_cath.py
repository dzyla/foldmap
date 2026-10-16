"""--domains cath: CATH's published chopping through the PDBe API (plan 12 Task 3). Offline: the responses are
fixtures in tests/data/cath, and any download attempt fails the test."""

import io
import json
import shutil
import urllib.error
from pathlib import Path

import numpy as np
import pytest

from foldmap import cath
from foldmap.cli import main, make_layout
from helpers import FID

DATA = Path(__file__).parent / "data"


@pytest.fixture
def cache(tmp_path, monkeypatch):
    """A cache holding the fixtures; any download attempt fails the test."""
    folder = tmp_path / "cath"
    folder.mkdir()
    for f in (DATA / "cath").iterdir():
        shutil.copy(f, folder / f.name)
    monkeypatch.setattr(cath, "CACHE", tmp_path)

    def offline(*a, **k):
        raise AssertionError("network used")

    monkeypatch.setattr(cath.urllib.request, "urlopen", offline)
    return folder


def test_chopping_reads_each_domain_and_its_segments(cache):
    got = cath.chopping("1LDM")
    assert [(d.domain, d.code, d.chain) for d in got] == [
        ("1ldmA01", "3.40.50.720", "A"),
        ("1ldmA02", "3.90.110.10", "A"),
    ]
    assert got[0].segments == [((1, ""), (162, ""))] and got[0].topology == "Rossmann fold"


def test_cath_domains_are_panels_named_by_their_cath_id(cache):
    lay, sses, _ = make_layout(FID / "1LDM.cif.gz", domains="cath", assembly="asu")
    names = [n for n, _ in lay.domains]
    assert names == ["1ldmA01", "1ldmA02"] and lay.domain_source == "CATH"
    by_id = {s.id: s for s in sses}
    assert all(by_id[k].last.seq <= 162 for k in dict(lay.domains)["1ldmA01"])
    assert all(by_id[k].first.seq >= 163 for k in dict(lay.domains)["1ldmA02"])
    assert "CATH" in lay.domain_note and "1ldmA01 (3.40.50.720, Rossmann fold)" in lay.domain_note


def test_copies_of_a_chain_take_its_cath_domains(cache):
    lay, _, _ = make_layout(FID / "1LDM.cif.gz", domains="cath")
    names = [n for n, _ in lay.domains]
    assert len(names) == 8 and names[0] == "1ldmA01 · A1"


def test_an_entry_without_a_cath_chopping_is_said_and_drawn_without_domains(cache, monkeypatch):
    def not_found(req, timeout=None):
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, io.BytesIO(b"{}"))

    monkeypatch.setattr(cath.urllib.request, "urlopen", not_found)
    lay, _, _ = make_layout(DATA / "1UBQ.cif", domains="cath")
    assert lay.domains == [] and "no CATH chopping for 1UBQ" in lay.domain_note
    assert json.loads((cache / "cath_1ubq.json").read_text()) == {}  # cached: not asked again


def test_an_unreachable_pdbe_is_an_error(cache, monkeypatch):
    def down(req, timeout=None):
        raise urllib.error.URLError("no route")

    monkeypatch.setattr(cath.urllib.request, "urlopen", down)
    with pytest.raises(ValueError, match="could not reach www.ebi.ac.uk"):
        cath.chopping("1ABC")


def test_cli_domains_cath_notes_the_source(cache, tmp_path, capsys):
    out = tmp_path / "t4l.svg"
    assert main(["plot", str(DATA / "2LZM.cif"), "-o", str(out), "--domains", "cath"]) == 0
    assert "domains from CATH (2LZM, via PDBe): 2lzmA00" in capsys.readouterr().err
    assert 'id="domain:2lzmA00"' in out.read_text()


def test_mcp_domains_from_cath_and_auto(monkeypatch):
    pytest.importorskip("mcp")
    from test_mcp import _spec_drawn_by

    assert _spec_drawn_by(monkeypatch, {"domains_from": "cath"}).domains == "cath"
    assert _spec_drawn_by(monkeypatch, {"domains_from": "auto"}).domains == "auto"
    assert _spec_drawn_by(monkeypatch, {"domains_from": "auto", "domains": {"N": ["α1"]}}).domains == [("N", ["α1"])]


def test_mcp_note_says_cath_domains(cache, tmp_path):
    pytest.importorskip("mcp")
    from test_mcp import call

    out = call(
        "draw_topology",
        {"source": str(DATA / "2LZM.cif"), "output_path": str(tmp_path / "t.svg"), "domains_from": "cath"},
    )
    assert "domains from CATH (2LZM, via PDBe)" in str(out)


def test_app_cath_switch_asks_for_cath_domains():
    from foldmap.app import layout_options

    assert layout_options({"cath_domains": True})["domains"] == "cath"
    assert layout_options({"cath_domains": False, "domains": "auto"})["domains"] == "auto"


def test_a_layout_file_keeps_where_its_domains_came_from(cache, tmp_path):
    """Saved, 'auto' domains are frozen as element ranges; the file says they were inferred, so the figure read
    back still says so (and CATH's still say CATH)."""
    import yaml

    for flag, source in (("auto", "inferred"), ("cath", "CATH")):
        lf = tmp_path / f"{flag}.yaml"
        out = tmp_path / f"{flag}.svg"
        assert main(["plot", str(DATA / "2LZM.cif"), "-o", str(out), "--domains", flag, "--save-layout", str(lf)]) == 0
        assert yaml.safe_load(lf.read_text())["domain_source"] == source
        from foldmap.spec import draw_spec, figure_spec

        assert draw_spec(figure_spec(DATA / "2LZM.cif", layout_file=lf)).layout.domain_source == source


def test_app_job_result_says_the_domains_are_inferred():
    from foldmap import app_jobs

    got = app_jobs.draw_job(DATA / "2LZM.cif", json.dumps({"domains": "auto"}))
    assert (
        got["provenance"][0]["text"]
        == "domains D1, D2 inferred from residue contacts, not a published domain assignment"
    )


def test_a_200_answer_that_is_not_json_is_a_clear_error(cache, monkeypatch):
    class Page:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self):
            return b"<html>maintenance</html>"

    monkeypatch.setattr(cath.urllib.request, "urlopen", lambda req, timeout=None: Page())
    with pytest.raises(ValueError, match="PDBe's CATH answer for 1abc is not JSON"):
        cath.chopping("1ABC")
    assert not (cache / "cath_1abc.json").exists()


def test_no_chopping_is_asked_again_after_30_days(cache, monkeypatch):
    import os
    import time

    (cache / "cath_9abc.json").write_text("{}")
    asked = []

    def answer(req, timeout=None):
        asked.append(req.full_url)
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, io.BytesIO(b"{}"))

    monkeypatch.setattr(cath.urllib.request, "urlopen", answer)
    assert cath.chopping("9ABC") == [] and asked == []  # fresh: the cached answer
    old = time.time() - 31 * 24 * 3600
    os.utime(cache / "cath_9abc.json", (old, old))
    assert cath.chopping("9ABC") == [] and len(asked) == 1  # a month on: CATH may have chopped it since


def test_a_chain_named_like_a_copy_takes_only_its_own_domains(cache):
    """An author chain "A1" is not a copy of chain "A" when CATH chops both: each takes its own domains, once."""
    from types import SimpleNamespace

    from foldmap.model import SSE, ResLabel

    labels = [ResLabel(c, n, "", "ALA") for c in ("A", "A1") for n in range(1, 41)]
    zero = np.zeros(3)
    sses = [
        SSE("H", c, k * 40 + 2, k * 40 + 30, labels[k * 40 + 2], labels[k * 40 + 30], zero, zero)
        for k, c in enumerate(("A", "A1"))
    ]
    (cache / "cath_1abc.json").write_text(
        json.dumps(
            {
                "1abc": {
                    "CATH": {
                        "1.10.10.10": {
                            "topology": "t",
                            "mappings": [
                                {
                                    "domain": d,
                                    "chain_id": c,
                                    "segment_id": 1,
                                    "start": {"author_residue_number": 1},
                                    "end": {"author_residue_number": 40},
                                }
                                for d, c in (("1abcA00", "A"), ("1abcA100", "A1"))
                            ],
                        }
                    }
                }
            }
        )
    )
    path = cache / "x.cif"
    path.write_text("data_1ABC\n_entry.id 1ABC\n")
    got, _ = cath.domains(path, SimpleNamespace(labels=labels), sses)
    assert got == [("1abcA00 · A", ["A:3-31"]), ("1abcA100 · A1", ["A1:3-31"])]


def test_element_id_domains_hold_each_element_once():
    from foldmap.cli import _domains

    assert _domains([("X", ["A:3-11"]), ("Y", ["A:3-11", "A:14-19"])], _t4l_sses(), {}) == [
        ("X", ["A:3-11"]),
        ("Y", ["A:14-19"]),
    ]


def _t4l_sses():
    _, sses, _ = make_layout(DATA / "2LZM.cif")
    return sses
