"""Plan 29's harness (scripts/app29.py) and what it adds to scripts/random100.py: the seeded tune/check split, the
app view judged as `foldmap check` judges a figure, the sampler's --no-alphafold switch (ruling Q61) and a short
sample taken and reported per stratum. No network: pools are stubbed, the figure is a fixture of tests/data."""

import importlib.util
import json
import os
import sys
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
    got = random100.app_view(DATA / "1UBQ.cif", "1UBQ", keep=tmp_path, judge=True)
    assert got["refused"] is None and got["error"] is None
    assert set(got["criteria"]) >= {"A1", "A4", "A6"} and got["floor"]["ok"] in (True, False)
    assert set(got["criteria"].values()) <= {"pass", "fail", "not measured"}
    assert (tmp_path / "app.svg").stat().st_size > 0 and (tmp_path / "app.png").stat().st_size > 0
    assert (tmp_path / "app.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert got["seconds"] > 0 and got["peak_mb"] == "not measured"  # app_jobs hands back no peak memory
    assert 0 <= got["draw_seconds"] <= got["seconds"]  # draw_job's own time, inside the job's wall time
    assert "<svg" in (tmp_path / "app.svg").read_text(encoding="utf-8")


def test_the_campaigns_app_view_is_not_judged(random100, monkeypatch):
    """Review finding 2 (ruling P29-3): random100's own app view draws and measures as before plan 29: no
    check_figure, no floor (judge defaults to False; only app29 asks for it)."""
    from foldmap import criteria, floor

    def never(*a, **k):
        raise AssertionError("the campaign's app view judged the figure")

    monkeypatch.setattr(criteria, "check_figure", never)
    monkeypatch.setattr(floor, "check", never)
    got = random100.app_view(DATA / "1UBQ.cif", "1UBQ")
    assert "measure_error" not in got and got["science"]["S1"] == []
    assert got["criteria"] == {} and got["floor"] is None


def test_a_judging_error_fails_closed(app29, random100, monkeypatch):
    """Review finding 5: a judge that raises leaves the figure not measured (measure_error, no criteria, no floor),
    never drawn and passed."""
    from foldmap import criteria

    def broken(*a, **k):
        raise RuntimeError("judge broke")

    monkeypatch.setattr(criteria, "check_figure", broken)
    got = random100.app_view(DATA / "1UBQ.cif", "1UBQ", judge=True)
    assert "judge broke" in got["measure_error"] and got["criteria"] == {} and got["floor"] is None
    assert app29.verdict(app29.row({"id": "1UBQ", "app": got}, "tune", "S1")) == "not measured"


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


def _planted_sample(app29, tmp_path) -> tuple[Path, list, list]:
    """A five-entry sample split 3 tune / 2 check, with a short S6 and an S12 not drawn."""
    strata = [{"name": "S1", "title": "X-ray, 1 protein chain", "n": 12}, {"name": "S6", "title": "Membrane", "n": 8}]
    strata.append({"name": "S12", "title": "AlphaFold DB models", "n": 3, "not_drawn": "not drawn (--no-alphafold)"})
    entries = [{"id": f"{k}AB{k}", "stratum": "S1"} for k in range(1, 6)]
    sample = tmp_path / "sample.json"
    doc = {"seed": 20261027, "query_date": "2026-10-27", "released_after": "2026-10-13"}
    doc.update(strata=strata, entries=entries, quotas={})
    sample.write_text(json.dumps(doc), encoding="utf-8")
    split = app29.add_split(sample, n_check=2)["split"]
    assert len(split["check"]) == 2 and len(split["tune"]) == 3
    return sample, split["tune"], split["check"]


def test_the_report_states_the_shortfall_per_stratum(app29, monkeypatch, tmp_path, capsys):
    """run --record writes a stamped report.json and report.md; the sample's shortfall is said per stratum, and the
    entries are counted per half of the split. The check half is sealed (rulings P29-4, P29-5): run's real per-entry
    writer runs (each child in-process here; only the drawing is stubbed), and no check id is named in any path or
    file under the run folder or the record but report-check.json, nor on the console. A child's stdout.txt and
    stderr.txt go to the entry's own folder: out/sealed/<key>/ for a check entry, where they may name it (a warning
    can carry the file's path), and nowhere else."""
    import contextlib

    sample, tune, check = _planted_sample(app29, tmp_path)
    r100 = app29.random100()
    cache = tmp_path / "cache"  # where resolve finds the files: outside the run folder
    cache.mkdir()
    streams = {}

    def resolve(i):
        (cache / f"{i}.cif").write_text(f"data_{i}\n", encoding="utf-8")
        return cache / f"{i}.cif"

    def app_view(path, name, keep=None, judge=False):  # the drawing, stubbed; it names its file as a warning would
        i = path.stem
        print(f"reading {path} as {name}", file=sys.stderr)
        (keep / "app.svg").write_text("<svg/>", encoding="utf-8")
        (keep / "app.png").write_bytes(b"\x89PNG\r\n\x1a\n")
        crit = {"A1": "pass", "A4": "fail" if i == check[0] else "pass"}  # a check id that fails a criterion
        got = {"refused": None, "error": None, "residues": 76, "criteria": crit}
        got.update(criteria_detail={"A4": f"{name} at {path}"}, seconds=1.0, peak_mb="not measured", svg_bytes=6)
        got.update(floor={"ok": i != tune[1], "reasons": []}, draw_seconds=0.5)
        if i == tune[0]:
            got.update(refused="time", message=f"{name} took longer than 90 s", criteria={}, floor=None)
        if i == check[1]:  # ... and one whose draw failed, saying the file's name
            got.update(error=f"ValueError: cannot read {path} ({path.name})", criteria={}, floor=None)
        return got

    def child(cmd, *, env, timeout, out, err):  # the real `entry`, in this process, its streams to run's files
        streams[cmd[cmd.index("entry") + 1]] = (out, err)
        with open(out, "w", encoding="utf-8") as fo, open(err, "w", encoding="utf-8") as fe:
            with contextlib.redirect_stdout(fo), contextlib.redirect_stderr(fe):
                code = app29.main(cmd[cmd.index("entry") :])
        return {"exit": code, "seconds": 1.0, "max_rss_mb": 50.0, "timed_out": False}

    monkeypatch.setattr(r100, "resolve", resolve)
    monkeypatch.setattr(r100, "app_view", app_view)
    monkeypatch.setattr(r100, "run_child", child)
    record, run = tmp_path / "record", tmp_path / "run"
    app29.main(["run", "--out", str(run), "--sample", str(sample), "--record", str(record)])
    rep = json.loads((record / "report.json").read_text(encoding="utf-8"))
    assert rep["source_stamp"] == app29.SOURCE_STAMP and rep["mixed_code"] is False
    short = {r["stratum"]: r["short"] for r in rep["sample"]["strata"]}
    assert short == {"S1": 7, "S6": 8, "S12": 3}
    assert rep["sample"]["split"] == {"tune": 3, "check": 2}
    t = rep["totals"]["all"]
    assert (t["entries"], t["drawn"], t["refused"], t["error"]) == (5, 3, 1, 1)
    assert (t["criteria_failed"], t["floor_failed"]) == (1, 1)
    c = rep["totals"]["check"]
    assert (c["entries"], c["error"], c["criteria_failed"]) == (2, 1, 1)  # counted ...
    md = (record / "report.md").read_text(encoding="utf-8")
    assert "| S6 | Membrane | 8 | 0 | 8 |" in md and "| S12 | AlphaFold DB models | 3 | 0 | 3 (not drawn" in md
    assert "Written by:" in md and f"{tune[0]}: refused: time" in md and f"| {tune[1]} |" in md
    said = capsys.readouterr().out
    assert f"{tune[0]}: refused" in said and all(i not in said for i in check)
    keys = {i: app29.sealed_key(i, 20261027) for i in check}
    for i in check:  # the child's streams: in the entry's sealed folder
        assert streams[i] == (run / "sealed" / keys[i] / "stdout.txt", run / "sealed" / keys[i] / "stderr.txt")
    for i in tune:
        assert streams[i] == (run / i / "stdout.txt", run / i / "stderr.txt") and (run / i / "app.png").is_file()
    for folder in (record, run):
        for f in folder.rglob("*"):
            assert all(i not in str(f.relative_to(folder)) for i in check), f  # ... never named in a path
            if not f.is_file() or f.name == "report-check.json":
                continue
            text = f.read_bytes().decode("utf-8", errors="replace")
            named = [i for i in check if i in text]
            if named:  # only a child's own streams, in its sealed folder, may name it
                assert f.name in ("stdout.txt", "stderr.txt") and f.parent == run / "sealed" / keys[named[0]], f
        sealed = json.loads((folder / "report-check.json").read_text(encoding="utf-8"))
        assert sorted(r["id"] for r in sealed["entries"]) == sorted(check)
        assert {r["id"]: r["dir"] for r in sealed["entries"]} == {i: f"sealed/{keys[i]}" for i in check}
    assert (run / "sealed" / keys[check[0]] / "app.svg").is_file()
    assert "stderr.txt" in {f.name for f in (run / "sealed" / keys[check[1]]).iterdir()}
    assert check[1] in (run / "sealed" / keys[check[1]] / "stderr.txt").read_text(encoding="utf-8")  # pinned
    assert sorted(r["id"] for r in app29.load_tune(record)["entries"]) == sorted(tune)
    both = app29.load_both(record, out=run)
    assert sorted(r["id"] for r in both["entries"]) == sorted(tune + check) and both["totals"] == rep["totals"]
    folders = {r["id"]: r["folder"] for r in both["entries"]}
    assert folders[check[0]] == str(run / "sealed" / keys[check[0]]) and folders[tune[0]] == str(run / tune[0])
    late = next(r for r in both["entries"] if r["id"] == check[1])
    assert late["error"].startswith("ValueError: cannot read <") and check[1] not in late["error"]


def test_ids_of_the_check_half_need_unseal(app29, monkeypatch, tmp_path, capsys):
    """Ruling P29-5 (N3): run --ids refuses an id of the check half, without naming it, unless --unseal (Task N+1)."""
    sample, tune, check = _planted_sample(app29, tmp_path)
    ran = []
    monkeypatch.setattr(app29, "run_one", lambda i, out, key=None: ran.append((i, key)) or {"id": i, "app": {}})
    argv = ["run", "--out", str(tmp_path / "run"), "--sample", str(sample), "--ids"]
    with pytest.raises(SystemExit, match="1 id\\(s\\) of the check half") as stop:
        app29.main(argv + [f"{tune[0]},{check[0]}"])
    assert check[0] not in str(stop.value) and ran == []
    app29.main(argv + [tune[0]])  # tune ids run as asked
    assert ran == [(tune[0], None)]
    app29.main(argv + [check[0], "--unseal"])  # Task N+1: still written sealed
    assert ran[-1] == (check[0], app29.sealed_key(check[0], 20261027))
    with pytest.raises(SystemExit):
        app29.main(["run", "--help"])
    assert "Task N+1" in " ".join(capsys.readouterr().out.split())  # --help says --unseal is Task N+1's


def test_entry_writes_its_record(app29, monkeypatch, tmp_path):
    r100 = app29.random100()
    monkeypatch.setattr(r100, "resolve", lambda i: DATA / "1UBQ.cif")
    monkeypatch.setattr(
        r100, "app_view", lambda path, name, keep=None, judge=False: {"refused": None, "keep": keep, "judge": judge}
    )
    assert app29.main(["entry", "1UBQ", "--out", str(tmp_path)]) == 0
    rec = json.loads((tmp_path / "1UBQ" / "app.json").read_text(encoding="utf-8"))
    assert rec["id"] == "1UBQ" and rec["source_stamp"] == app29.SOURCE_STAMP
    assert rec["app"]["keep"] == str(tmp_path / "1UBQ") and rec["file"]["name"] == "1UBQ.cif"
    assert rec["app"]["judge"] is True  # app29 judges; the campaign's view does not


def test_run_ui_drives_the_page_for_the_tune_half_only(app29, monkeypatch, tmp_path, capsys):
    """Plan 29 Task 2: --ui runs the page session (the `ui` child) for each tune-half id and never for a check-half
    one; the report's `ui` is its record ("not run" for the check half, in report-check.json only), its totals count
    the sessions that failed, and each row lists the sha256 of the files its folder kept."""
    sample, tune, check = _planted_sample(app29, tmp_path)
    r100 = app29.random100()
    paged = []

    def child(cmd, *, env, timeout, out, err):
        what, i = cmd[cmd.index(str(app29.SCRIPT)) + 1 :][:2]
        keep = Path(out).parent
        if what == "entry":
            rec = {"source_stamp": app29.SOURCE_STAMP, "app": {"refused": None, "error": None}}
            (keep / "app.svg").write_text("<svg/>", encoding="utf-8")
            (keep / "app.json").write_text(json.dumps(rec), encoding="utf-8")
        else:
            assert what == "ui" and env["FOLDMAP_APP_LOCAL_FILES"] == "0" and timeout == app29.UI_TIMEOUT
            paged.append(i)
            downs = {"dl_svg": {"bytes": 6, "sha256": "x", "ok": i != tune[0]}}
            ui = {"exception": None, "harness_error": None, "errors": [], "stage": "edit", "downloads": downs}
            ui["steps"] = ["load", "theme trace", "downloads"]
            (keep / "ui.svg").write_text("<svg/>", encoding="utf-8")
            (keep / "ui.json").write_text(json.dumps({"source_stamp": app29.SOURCE_STAMP, "ui": ui}), "utf-8")
        return {"exit": 0, "seconds": 1.0, "max_rss_mb": 50.0, "timed_out": False}

    monkeypatch.setattr(r100, "run_child", child)
    record, run = tmp_path / "record", tmp_path / "run"
    app29.main(["run", "--ui", "--out", str(run), "--sample", str(sample), "--record", str(record)])
    assert sorted(paged) == sorted(tune)
    rep = json.loads((record / "report.json").read_text(encoding="utf-8"))
    assert rep["mixed_code"] is False
    rows = {r["id"]: r for r in rep["entries"]}
    assert sorted(rows) == sorted(tune) and all(rows[i]["ui"]["stage"] == "edit" for i in tune)
    assert set(rows[tune[1]]["files"]) == {"app.svg", "ui.svg"}  # ui.json and the child's streams are not files kept
    assert rep["totals"]["tune"]["ui_run"] == 3 and rep["totals"]["tune"]["ui_failed"] == 1  # a download not whole
    assert rep["totals"]["check"]["ui_run"] == 0
    sealed = json.loads((record / "report-check.json").read_text(encoding="utf-8"))
    assert all(r["ui"] == "not run" for r in sealed["entries"])
    md = (record / "report.md").read_text(encoding="utf-8")
    assert "## The page (run --ui)" in md and f"| {tune[0]} | edit |" in md and all(i not in md for i in check)
    said = capsys.readouterr().out
    assert f"{tune[0]}: drawn; page failed" in said and all(i not in said for i in check)
    with pytest.raises(SystemExit, match="check half"):
        app29.main(["ui", check[0], "--out", str(run), "--sample", str(sample)])


@pytest.mark.skipif(sys.platform == "win32", reason="isolated jobs need forkserver")
def test_the_page_session_draws_changes_the_theme_and_presses_every_download(app29, monkeypatch, tmp_path):
    """ui_session on 1UBQ (from the cache, no network), hosted: drawn, the theme changed, six files, each whole."""
    from app_helpers import fresh_caches, seed_cache

    before = os.environ.get("FOLDMAP_APP_LOCAL_FILES")
    (tmp_path / "pdb").mkdir()
    seed_cache(monkeypatch, tmp_path / "pdb", "1UBQ")
    fresh_caches()
    keep = tmp_path / "keep"
    keep.mkdir()
    ui = app29.ui_session("1UBQ", keep)
    assert ui["exception"] is None and ui["harness_error"] is None and ui["errors"] == [], ui
    assert ui["stage"] == "edit" and ui["steps"][0] == "load" and ui["steps"][1].startswith("theme ")
    assert ui["steps"][-1] == "downloads" and set(ui["downloads"]) == set(app29.UI_FILES)
    assert all(d["ok"] and d["bytes"] > 0 for d in ui["downloads"].values()), ui["downloads"]
    assert sorted(f.name for f in keep.iterdir()) == sorted(app29.UI_FILES.values())
    assert os.environ.get("FOLDMAP_APP_LOCAL_FILES") == before  # hosted for the session only
    fresh_caches()


def test_a_failed_app_draw_fails_the_view(random100, monkeypatch, tmp_path, capsys):
    """app_view now catches a failed draw ("error"); random100's app view gives no answer for it, as when it crashed:
    A2b fails ("gave no answer", review finding 1), it exits 1 so A1 fails, and says why on stderr (A1's detail reads
    its tail); the error is kept under app_error."""
    from foldmap import cli

    monkeypatch.setattr(cli, "draw_spec", cli.draw_spec)  # cmd_view wraps it: restored after the test
    monkeypatch.setattr(random100, "app_view", lambda path, name: {"refused": None, "error": "ValueError: no protein"})
    code = random100.main(["view", "--name", "X", "--work", str(tmp_path), "app", str(DATA / "1UBQ.cif")])
    assert code == 1 and "app: ValueError: no protein" in capsys.readouterr().err
    rec = json.loads((tmp_path / "X-app.json").read_text(encoding="utf-8"))
    assert rec["exit"] == 1 and rec["app"] is None and rec["app_error"]["error"] == "ValueError: no protein"
    view = {"result": rec, "timed_out": False, "exit": 1, "traceback": None, "stderr_tail": ["app: ValueError"]}
    judged = random100.criteria({"views": {"app": {**view, "written": {}}}})
    assert judged["A2b"]["pass"] is False and "gave no answer" in judged["A2b"]["detail"]
    assert judged["A1"]["pass"] is False


def test_a_child_that_timed_out_keeps_the_report_quotable(app29, monkeypatch, tmp_path):
    """Review finding 4: an entry whose process was killed is recorded by run, under run's own stamp, so one timeout
    does not make the report mixed_code (never quoted)."""
    r100 = app29.random100()

    def child(cmd, *, env, timeout, out, err):
        i = cmd[cmd.index("entry") + 1]
        if i == "2BBB":
            return {"exit": -9, "seconds": timeout, "max_rss_mb": None, "timed_out": True}
        rec = {"id": i, "source_stamp": app29.SOURCE_STAMP, "app": {"refused": None, "error": None}}
        (Path(cmd[cmd.index("--out") + 1]) / i / "app.json").write_text(json.dumps(rec), encoding="utf-8")
        return {"exit": 0, "seconds": 1.0, "max_rss_mb": 100.0, "timed_out": False}

    monkeypatch.setattr(r100, "run_child", child)
    run = tmp_path / "run"
    app29.main(["run", "--out", str(run), "--ids", "1AAA,2BBB", "--sample", str(tmp_path / "none.json")])
    rep = json.loads((run / "report.json").read_text(encoding="utf-8"))
    assert rep["mixed_code"] is False
    late = next(r for r in rep["entries"] if r["id"] == "2BBB")
    assert late["error"] == "no result (timed out)" and late["process"]["timed_out"] is True


def test_the_committed_split_is_the_seeded_split():
    """Review finding 8: docs/validation/app29/sample.json holds split(its entries), 47 tune / 32 check (P29-2)."""
    path = ROOT / "docs" / "validation" / "app29" / "sample.json"
    if not path.is_file():
        pytest.skip("docs/validation is not part of the sdist")
    app29 = _load("app29")
    doc = json.loads(path.read_text(encoding="utf-8"))
    assert doc["split"] == app29.split([e["id"] for e in doc["entries"]], seed=doc["seed"])
    assert (len(doc["split"]["tune"]), len(doc["split"]["check"])) == (47, 32)


def test_random100_is_loaded_once_and_whole_from_many_threads():
    """run --jobs N calls app29.random100() from N threads at once: each gets the one module, fully run (the
    before measurement's first run died on 'random100 has no attribute run_child': the module was published before
    it ran). Its global stays unset until the module has run."""
    from concurrent.futures import ThreadPoolExecutor

    fresh = _load("app29")
    seen_while_loading = []
    real = importlib.util.module_from_spec

    def watched(spec):
        mod = real(spec)
        loader_exec = spec.loader.exec_module

        def exec_module(m):
            seen_while_loading.append(fresh._R100)
            loader_exec(m)

        spec.loader.exec_module = exec_module
        return mod

    import threading

    gate = threading.Barrier(8)

    def one(_):
        gate.wait()
        return fresh.random100()

    try:
        importlib.util.module_from_spec = watched
        with ThreadPoolExecutor(8) as pool:
            got = list(pool.map(one, range(8)))
    finally:
        importlib.util.module_from_spec = real
    assert seen_while_loading == [None]  # loaded once, and not published while it ran
    assert all(m is got[0] for m in got) and hasattr(got[0], "run_child")
