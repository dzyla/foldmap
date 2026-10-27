"""Plan 29's harness (scripts/app29.py) and what it adds to scripts/random100.py: the seeded tune/check split, the
app view judged as `foldmap check` judges a figure, the sampler's --no-alphafold switch (ruling Q61) and a short
sample taken and reported per stratum. No network: pools are stubbed, the figure is a fixture of tests/data."""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_t29", ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def app29():
    return _load("app29")


@pytest.fixture(scope="module")
def random100():
    return _load("random100")


def test_split_is_seeded_disjoint_and_complete(app29):
    ids = [f"{i}ABC" for i in range(100)]
    a = app29.split(ids, seed=20261027, n_check=40)
    assert a == app29.split(list(reversed(ids)), seed=20261027, n_check=40)
    assert len(a["check"]) == 40 and set(a["check"]).isdisjoint(a["tune"])
    assert sorted(a["check"] + a["tune"]) == sorted(ids)
    assert a != app29.split(ids, seed=1, n_check=40)  # the seed chooses


def test_the_default_split_is_60_40_of_the_sample(app29):
    """Ruling P29-2: check is round(0.4 x n) of the sample: 32 of the 79 drawn, 47 tune."""
    ids = [f"{i}XYZ" for i in range(79)]
    a = app29.split(ids)
    assert (len(a["tune"]), len(a["check"])) == (47, 32)
    assert a["check"] == app29.split(ids, n_check=32)["check"]  # the same ranking, cut at 32


def test_app_view_judges_the_apps_own_figure(random100, tmp_path):
    got = random100.app_view(DATA / "1UBQ.cif", "1UBQ", keep=tmp_path)
    assert got["refused"] is None and got["error"] is None
    assert set(got["criteria"]) >= {"A1", "A4", "A6"} and got["floor"]["ok"] in (True, False)
    assert set(got["criteria"].values()) <= {"pass", "fail", "not measured"}
    assert (tmp_path / "app.svg").stat().st_size > 0 and (tmp_path / "app.png").stat().st_size > 0
    assert (tmp_path / "app.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert got["seconds"] > 0 and got["peak_mb"] == "not measured"  # app_jobs hands back no peak memory
    assert "<svg" in (tmp_path / "app.svg").read_text(encoding="utf-8")


def test_the_draw_is_plan_29s_call(app29):
    """Plan 29 Task 1 Step 1, as given, with --no-alphafold (ruling Q61)."""
    d = app29.DRAW
    assert d[0] == "sample" and "--no-alphafold" in d
    for flag, value in (("--seed", "20261027"), ("--date", "2026-10-27"), ("--released-after", "2026-10-13")):
        assert d[d.index(flag) + 1] == value
    excluded = [d[k + 1] for k, x in enumerate(d) if x == "--exclude"]
    names = ("random100", "random100-heldout", "random100-third")
    assert excluded == [f"docs/validation/{n}/sample.json" for n in names]
    assert d[d.index("--out") + 1] == "docs/validation/app29/sample.json"


def _stub_network(random100, monkeypatch, asked: list):
    """Every experimental pool holds one id, 1AAA (a fixture file); the computational one holds AF_AFP00001F1."""

    def pool(q):
        content = q["request_options"]["results_content_type"][0]
        asked.append(content)
        return ["AF_AFP00001F1"] if content == "computational" else ["1AAA"]

    feats = {"protein": True, "na": False, "ca_only": False, "method": "X-RAY DIFFRACTION", "chains": 1}
    feats.update(deposited_residues=76, has=[], accessions=["P0CG48"], assembly1_chains=1, n_linked=0, rna_nt=0)
    feats.update(model_name="m")
    monkeypatch.setattr(random100, "fetch_pool", pool)
    monkeypatch.setattr(random100, "resolve", lambda i: DATA / "1UBQ.cif")
    monkeypatch.setattr(random100, "features", lambda path, contacts=False: dict(feats))


def test_no_alphafold_draws_no_alphafold_model_and_a_short_sample_is_written(random100, monkeypatch, tmp_path):
    """Ruling Q61: --no-alphafold never queries the AlphaFold DB pool (it has no release window); S12 stays short and
    says why. A pool short of its quota gives what it has: the sample is still written, its quotas say how short."""
    asked: list = []
    _stub_network(random100, monkeypatch, asked)
    out = tmp_path / "sample.json"
    argv = ["sample", "--date", "2026-10-27", "--released-after", "2026-10-13", "--no-alphafold", "--out", str(out)]
    code = random100.main(argv)
    assert code == 1  # short of its quotas: written, and said
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert "computational" not in asked and asked
    assert [e["id"] for e in doc["entries"]] == ["1AAA"]  # one id in every pool: drawn once, then "already in"
    assert doc["not_drawn"] == {"S12": random100.NO_ALPHAFOLD}
    assert doc["quotas"]["S12:S"] == {"target": 1, "reached": 0} and doc["quotas"]["S1:S"]["reached"] == 1
    s12 = next(st for st in doc["strata"] if st["name"] == "S12")
    assert s12["not_drawn"] and "query" not in s12
    asked.clear()
    with_af = random100.sample(20261027, "2026-10-27", after="2026-10-13")
    assert "computational" in asked and "not_drawn" not in with_af
    assert any(e["id"] == "AF_AFP00001F1" for e in with_af["entries"])


def test_the_report_states_the_shortfall_per_stratum(app29, monkeypatch, tmp_path):
    """run --record writes a stamped report.json and report.md; the sample's shortfall is said per stratum, and the
    entries are counted per half of the split."""
    strata = [{"name": "S1", "title": "X-ray, 1 protein chain", "n": 12}, {"name": "S6", "title": "Membrane", "n": 8}]
    strata.append({"name": "S12", "title": "AlphaFold DB models", "n": 3, "not_drawn": "not drawn (--no-alphafold)"})
    entries = [{"id": f"{k}AB{k}", "stratum": "S1"} for k in range(1, 6)]
    sample = tmp_path / "sample.json"
    doc = {"seed": 20261027, "query_date": "2026-10-27", "released_after": "2026-10-13"}
    doc.update(strata=strata, entries=entries, quotas={})
    sample.write_text(json.dumps(doc), encoding="utf-8")
    split = app29.add_split(sample, n_check=2)["split"]
    assert len(split["check"]) == 2 and len(split["tune"]) == 3

    def one(i, out):
        crit = {"A1": "pass", "A4": "fail" if i == "1AB1" else "pass"}
        app = {"refused": None, "error": None, "residues": 76, "criteria": crit, "criteria_detail": {"A4": "x"}}
        app.update(floor={"ok": i != "2AB2", "reasons": []}, seconds=1.0, peak_mb="not measured", svg_bytes=10)
        if i == "3AB3":
            app = {"refused": "time", "error": None, "residues": 9000, "criteria": {}, "floor": None}
        return {"id": i, "source_stamp": app29.SOURCE_STAMP, "app": app, "process": {"exit": 0}}

    monkeypatch.setattr(app29, "run_one", one)
    record = tmp_path / "record"
    app29.main(["run", "--out", str(tmp_path / "run"), "--sample", str(sample), "--record", str(record)])
    rep = json.loads((record / "report.json").read_text(encoding="utf-8"))
    assert rep["source_stamp"] == app29.SOURCE_STAMP and rep["mixed_code"] is False
    assert (tmp_path / "run" / "report.json").is_file()
    short = {r["stratum"]: r["short"] for r in rep["sample"]["strata"]}
    assert short == {"S1": 7, "S6": 8, "S12": 3}
    assert rep["sample"]["split"] == {"tune": 3, "check": 2}
    t = rep["totals"]["all"]
    assert (t["entries"], t["drawn"], t["refused"], t["criteria_failed"], t["floor_failed"]) == (5, 4, 1, 1, 1)
    assert rep["totals"]["tune"]["entries"] + rep["totals"]["check"]["entries"] == 5
    md = (record / "report.md").read_text(encoding="utf-8")
    assert "| S6 | Membrane | 8 | 0 | 8 |" in md and "| S12 | AlphaFold DB models | 3 | 0 | 3 (not drawn" in md
    assert "Written by:" in md and "3AB3: refused: time" in md


def test_entry_writes_its_record_and_run_ui_waits_for_task_2(app29, monkeypatch, tmp_path):
    r100 = app29.random100()
    monkeypatch.setattr(r100, "resolve", lambda i: DATA / "1UBQ.cif")
    monkeypatch.setattr(r100, "app_view", lambda path, name, keep=None: {"refused": None, "error": None, "keep": keep})
    assert app29.main(["entry", "1UBQ", "--out", str(tmp_path)]) == 0
    rec = json.loads((tmp_path / "1UBQ" / "app.json").read_text(encoding="utf-8"))
    assert rec["id"] == "1UBQ" and rec["source_stamp"] == app29.SOURCE_STAMP
    assert rec["app"]["keep"] == str(tmp_path / "1UBQ") and rec["file"]["name"] == "1UBQ.cif"
    with pytest.raises(NotImplementedError, match="Task 2"):
        app29.main(["run", "--ui", "--out", str(tmp_path), "--ids", "1UBQ"])


def test_a_failed_app_draw_fails_the_view(random100, monkeypatch, tmp_path, capsys):
    """app_view now catches a failed draw ("error"); random100's app view still exits 1 on it, so A1 fails as it did
    when the failure was a crash, and says why on stderr (A1's detail reads its tail)."""
    from foldmap import cli

    monkeypatch.setattr(cli, "draw_spec", cli.draw_spec)  # cmd_view wraps it: restored after the test
    monkeypatch.setattr(random100, "app_view", lambda path, name: {"refused": None, "error": "ValueError: no protein"})
    code = random100.main(["view", "--name", "X", "--work", str(tmp_path), "app", str(DATA / "1UBQ.cif")])
    assert code == 1 and "app: ValueError: no protein" in capsys.readouterr().err
    assert json.loads((tmp_path / "X-app.json").read_text(encoding="utf-8"))["exit"] == 1
