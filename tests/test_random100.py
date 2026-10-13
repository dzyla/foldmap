"""The random-100 campaign script (scripts/random100.py): the sampler's draw, strata, quotas and pairs, and the
runner's measures, contact sheet and report. No network (R10): pools are synthetic lists, the search service's
answers are recorded responses, and the runner reads the fixtures in tests/data."""

import importlib.util
import io
import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "random100.py"
DATA = Path(__file__).parent / "data"


def _load():
    spec = importlib.util.spec_from_file_location("random100", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def r100():
    return _load()


POOL = [f"{n}AB{c}" for n in range(1, 10) for c in "CDEFGH"]  # 54 synthetic ids


# ------------------------------------------------------------------------------------------------ sampler


def _draw_in_child(seed: str) -> list[str]:
    code = (
        "import importlib.util, json, sys\n"
        f"spec = importlib.util.spec_from_file_location('random100', {str(SCRIPT)!r})\n"
        "mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)\n"
        f"taken, _ = mod.draw({POOL!r}, '20261013:S1:S', 7, lambda i: (True, ''))\n"
        "print(json.dumps(taken))\n"
    )
    env = {**os.environ, "PYTHONHASHSEED": seed}
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True)
    return json.loads(done.stdout)


def test_draw_is_reproducible_and_stratified(r100):
    a, b = _draw_in_child("0"), _draw_in_child("1")
    assert a == b and len(a) == 7
    taken, skipped = r100.draw(POOL, "20261013:S1:S", 7, lambda i: (True, ""))
    assert taken == a  # and the same here, whatever this process's hash seed
    other, _ = r100.draw(POOL, "20261013:S1:M", 7, lambda i: (True, ""))
    assert other != taken  # another key, another order
    assert r100.order(POOL, "20261013:S1:S")[:7] == taken

    def accept(i):  # every id whose second character is odd is refused
        return (int(i[0]) % 2 == 0, "" if int(i[0]) % 2 == 0 else "odd")

    taken, skipped = r100.draw(POOL, "20261013:S1:S", 5, accept)
    assert len(taken) == 5 and all(int(i[0]) % 2 == 0 for i in taken)
    assert skipped and all(reason == "odd" and int(i[0]) % 2 for i, reason in skipped)
    seen = r100.order(POOL, "20261013:S1:S")
    last = seen.index(taken[-1])
    assert [i for i, _ in skipped] == [i for i in seen[:last] if int(i[0]) % 2]  # every skip before the last take


def test_pool_order_does_not_matter(r100):
    shuffled = list(POOL)
    random.Random(4).shuffle(shuffled)
    assert r100.draw(shuffled, "k", 9, lambda i: (True, ""))[0] == r100.draw(POOL, "k", 9, lambda i: (True, ""))[0]


def test_draw_looks_at_most_200_candidates(r100):
    pool = [f"{n:04d}" for n in range(1000)]
    taken, skipped = r100.draw(pool, "k", 5, lambda i: (False, "no"))
    assert taken == [] and len(skipped) == 200


def test_strata_add_up(r100):
    strata = r100.STRATA
    assert sum(s["n"] for s in strata) == 100
    for s in strata:
        if s["bins"]:
            assert sum(s["bins"].values()) == s["n"], s["name"]
    totals = {b: sum(s["bins"].get(b, 0) for s in strata) for b in ("S", "M", "L", "XL")}
    assert totals == {"S": 21, "M": 29, "L": 23, "XL": 17}
    pairs = next(s for s in strata if s["name"] == "S11")
    assert pairs["n"] == 10 and pairs["pairs"] == {"apo/holo": 2, "monomer/oligomer": 2, "open/closed": 1}
    extras = {(s["name"], x["tag"]): x["n"] for s in strata for x in s.get("extra", [])}
    assert extras == {("S4", "helical"): 1, ("S4", "big"): 2, ("S4", "cap"): 3, ("S8", "ribosome"): 1}


def test_query_has_release_cutoff(r100):
    date = {"attribute": "rcsb_accession_info.initial_release_date", "operator": "less_or_equal"}
    for s in r100.STRATA:
        for b in s["bins"] or [None]:
            q = r100.query(r100.stratum_nodes(s, b), "2026-10-13", s.get("content", "experimental"))
            assert q["query"]["logical_operator"] == "and"
            assert q["request_options"]["return_all_hits"] is True
            nodes = [n["parameters"] for n in q["query"]["nodes"] if n["type"] == "terminal"]
            assert {**date, "value": "2026-10-13"} in nodes, s["name"]
    s12 = next(s for s in r100.STRATA if s["name"] == "S12")
    assert r100.query(r100.stratum_nodes(s12, "S"), "2026-10-13", "computational")["request_options"][
        "results_content_type"
    ] == ["computational"]


class _Answer(io.BytesIO):
    def __init__(self, body: bytes, status: int = 200):
        super().__init__(body)
        self.status = status

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_fetch_pool_reads_a_recorded_answer(r100, monkeypatch):
    # recorded from search.rcsb.org (results_verbosity compact), cut to three ids
    recorded = b'{"query_id":"a1","result_type":"entry","total_count":3,"result_set":["105M","100D","101D"]}'
    sent = []

    def urlopen(req, timeout):
        sent.append(json.loads(req.data))
        return _Answer(recorded)

    monkeypatch.setattr(r100.urllib.request, "urlopen", urlopen)
    q = r100.query([], "2026-10-13")
    assert r100.fetch_pool(q) == ["100D", "101D", "105M"]
    assert sent == [q]
    monkeypatch.setattr(r100.urllib.request, "urlopen", lambda req, timeout: _Answer(b"", 204))
    assert r100.fetch_pool(q) == []  # no hits: the service answers 204 with no body


def test_pool_record_is_a_count_and_a_hash(r100):
    import hashlib

    rec = r100.pool_record(["B", "A"])
    assert rec == {"count": 2, "sha256": hashlib.sha256(b"A\nB").hexdigest()}


def _entries(features: dict, stratum_of: dict) -> list[dict]:
    return [{"id": i, "stratum": stratum_of[i], "bin": "S", "order": k} for k, i in enumerate(features)]


def test_top_up_replaces_last_lacking_entry_deterministically(r100):
    features = {"A1": set(), "A2": {"glycan"}, "A3": set(), "B1": set(), "B2": {"modified"}}
    sample = _entries(features, {"A1": "S1", "A2": "S1", "A3": "S1", "B1": "S2", "B2": "S2"})
    # the pools' seeded orders after what was drawn: S1 offers G1 (glycan), S2 offers G2 (glycan)
    rest = {("S1", "S"): ["X1", "G1", "G3"], ("S2", "S"): ["G2"]}
    has = {"X1": set(), "G1": {"glycan"}, "G2": {"glycan"}, "G3": {"glycan"}}

    def run():
        s = [dict(e) for e in sample]
        feats = {**{k: set(v) for k, v in features.items()}, **{k: set(v) for k, v in has.items()}}
        replaced = r100.top_up(
            s, feats, {k: list(v) for k, v in rest.items()}, {"glycan": 3, "modified": 1}, accept=lambda i: (True, "")
        )
        return s, replaced

    s, replaced = run()
    assert [e["id"] for e in s] == ["A1", "A2", "G1", "G2", "B2"]  # A3 (last of S1), then B1 (S2: B2 is modified)
    assert replaced == [
        {"id": "A3", "reason": "top-up: lacks glycan", "replacement": "G1"},
        {"id": "B1", "reason": "top-up: lacks glycan", "replacement": "G2"},
    ]
    assert run() == (s, replaced)


def test_top_up_never_drops_another_quota_below_its_target(r100):
    features = {"A1": {"modified"}, "A2": set()}
    sample = _entries(features, {"A1": "S1", "A2": "S1"})
    s = [dict(e) for e in sample]
    feats = {**features, "A1": {"modified"}, "A2": set(), "G1": {"glycan"}, "G2": {"glycan"}}
    replaced = r100.top_up(
        s, feats, {("S1", "S"): ["G1", "G2"]}, {"glycan": 2, "modified": 1}, accept=lambda i: (True, "")
    )
    assert [e["id"] for e in s] == ["A1", "G1"]  # A1 is the only modified entry: kept, the quota stays short
    assert [r["id"] for r in replaced] == ["A2"]


def test_pairs_fill_each_type_once(r100):
    seeds = [f"S{k}" for k in range(12)]
    partners = {s: [f"{s}p{j}" for j in range(3)] for s in seeds}
    kinds = ["apo/holo", "monomer/oligomer", "open/closed"]

    def classify(seed, partner):  # every partner of seed k is of type k % 3 (S9..S11 of no type)
        k = int(seed[1:])
        return set() if k >= 9 else {kinds[k % 3]}

    got, skipped = r100.pairs(seeds, lambda s: partners[s], classify)
    counts = {k: sum(1 for p in got if p["type"] == k) for k in kinds}
    assert counts == {"apo/holo": 2, "monomer/oligomer": 2, "open/closed": 1}
    assert all(p["partner"] in partners[p["seed"]] for p in got)
    assert got == r100.pairs(list(reversed(seeds)), lambda s: partners[s], classify)[0]  # sorted, then seeded
    order = r100.order(seeds, "20261013:S11")
    assert [p["seed"] for p in got] == [s for s in order if s in {p["seed"] for p in got}]


def test_pairs_try_at_most_eight_partners(r100):
    tried = []

    def classify(seed, partner):
        tried.append(partner)
        return set()

    r100.pairs(["S0"], lambda s: [f"P{k:02d}" for k in range(20)], classify)
    assert len(tried) == 8
    want = sorted(f"P{k:02d}" for k in range(20))
    random.Random("20261013:S11:S0").shuffle(want)
    assert tried == want[:8]


# ------------------------------------------------------------------------------------------------ runner


def test_entry_measures_a_fixture(r100, tmp_path):
    code = r100.main(["entry", str(DATA / "1UBQ.cif"), "--out", str(tmp_path), "--views", "quick"])
    saved = json.loads((tmp_path / "1UBQ.json").read_text(encoding="utf-8"))
    assert code == 0
    assert set(saved["criteria"]) == set(r100.CRITERIA) and {f"A{k}" for k in range(1, 21)} <= set(saved["criteria"])
    assert saved["pass"] is True, saved["failed"]
    assert set(saved["views"]) == set(r100.QUICK)
    c = saved["criteria"]
    assert c["A1"]["pass"] and c["A11"]["pass"] and c["A12"]["pass"] and c["A13"]["pass"]
    assert c["A16"]["pass"] is None and c["A17"]["pass"] is None and c["A19"]["pass"] is None  # does not apply
    d = saved["views"]["default"]["result"]
    assert d["checks"]["element_overlaps"] == 0 and d["signature_sha256"] and d["svg_sha256"]
    assert d["route"]["records"] >= 1 and d["route"]["loops"] > 0  # route's summary record was read
    assert (tmp_path / "work" / "1UBQ" / "1UBQ-trace.json").is_file()


def test_entry_reports_a_crash_without_dying(r100, tmp_path, monkeypatch):
    real = r100.view_command

    def crashing(view, *args, **kw):
        cmd, env = real(view, *args, **kw)
        if view == "seed1":  # this view's process dies with a traceback before drawing
            cmd = [sys.executable, "-I", "-c", "raise RuntimeError('boom')"]
        return cmd, env

    monkeypatch.setattr(r100, "view_command", crashing)
    e = r100.run_entry(str(DATA / "1UBQ.cif"), tmp_path, ("default", "seed1"))
    assert e["pass"] is False
    assert e["criteria"]["A1"]["pass"] is False and "seed1" in e["criteria"]["A1"]["detail"]
    assert "traceback" in e["criteria"]["A1"]["detail"]
    assert e["views"]["default"]["result"]["checks"]["element_overlaps"] == 0  # the other view still measured
    assert e["criteria"]["A5"]["pass"] is True and e["criteria"]["A12"]["pass"] is False


def test_criteria_fail_closed_on_unmeasured_checks(r100):
    view = {"exit": 0, "timed_out": False, "traceback": False, "python_warnings": [], "written": {}, "seconds": 1.0}
    result = {"checks": {"element_overlaps": 0, "text_collisions": None, "skipped": {"text_collisions": "font"}}}
    e = {"views": {"default": {**view, "result": result}}}
    c = r100.criteria(e)
    assert c["A5"]["pass"] is True and c["A6"]["pass"] is False and "not measured (font)" in c["A6"]["detail"]
    warned = {"views": {"default": {**view, "python_warnings": ["DeprecationWarning: x"], "result": result}}}
    assert r100.criteria(warned)["A1"]["pass"] is False


def test_a_ca_only_entry_must_be_refused_in_one_line(r100):
    view = {"timed_out": False, "python_warnings": [], "written": {}, "seconds": 1.0, "result": None}
    refused = {"exit": 1, "traceback": False, "foldmap_said": ["foldmap: x.cif: only CA atoms"]}
    e = {"refusal": True, "views": {"default": {**view, **refused}}}
    c = r100.criteria(e)
    assert c["A1c"]["pass"] is True and c["A1"]["pass"] is None
    drawn = {"exit": 0, "traceback": False, "foldmap_said": []}
    assert r100.criteria({"refusal": True, "views": {"default": {**view, **drawn}}})["A1c"]["pass"] is False
    crashed = {"exit": 1, "traceback": True, "foldmap_said": []}
    assert r100.criteria({"refusal": True, "views": {"default": {**view, **crashed}}})["A1c"]["pass"] is False


def test_dssp_agreement_against_the_reference(r100):
    rows = []
    for line in (DATA / "dssp" / "1UBQ.ss").read_text(encoding="utf-8").splitlines():
        chain, first, codes = line.split()
        rows += [(chain, int(first) + k, "", code) for k, code in enumerate(codes)]
    got = r100.dssp_agreement(DATA / "1UBQ.cif", rows)
    assert got["agree"] == 1.0 and got["residues"] == 76
    rows[10] = (*rows[10][:3], "G" if rows[10][3] != "G" else "H")
    assert r100.dssp_agreement(DATA / "1UBQ.cif", rows)["agree"] == round(1 - 1 / 76, 4)


def test_route_summary_record_is_read(r100):
    import logging

    handler = r100._RouteLog()
    log = logging.getLogger("foldmap.route")
    record = log.makeRecord(
        "foldmap.route",
        logging.DEBUG,
        "f",
        1,
        "%d loops: %d fallback curves, %d of them crossing an element, %d kink links rejected, "
        "%d drawn nearer than %.2f to an element, %d drawings rejected as crossing",
        (9, 1, 1, 5, 2, 0.15, 3),
        None,
    )
    handler.emit(record)
    handler.emit(record)
    s = handler.summary()
    assert s["loops"] == 18 and s["fallback_curves"] == 2 and s["tighter_than_inflate"] == 4
    assert s["fallback_curves_crossing_an_element"] == 2 and s["kink_links_rejected"] == 10
    assert s["rejected_as_crossing"] == 6 and s["records"] == 2 and "inflate" not in s


def test_route_summary_record_is_the_one_route_logs(r100, caplog):
    """The runner reads the record route.route_loops writes (its format and argument count)."""
    import logging

    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    caplog.set_level(logging.DEBUG, logger="foldmap.route")
    layout, sses, bb = make_layout(DATA / "1UBQ.cif")
    route_loops(layout, sses, bb)
    handler = r100._RouteLog()
    for record in caplog.records:
        handler.emit(record)
    s = handler.summary()
    assert s["records"] == 1 and s["loops"] > 0 and s["kink_links_rejected"] >= 0


def test_pair_check_matches_elements_by_residue_range(r100):
    def entry(elements):
        return {"views": {"beta": {"result": {"elements": elements}}}}

    a = [
        {"id": "A:1-8", "chain": "A", "kind": "E", "first": 1, "last": 8, "label": "β1"},
        {"id": "A:12-30", "chain": "A", "kind": "H", "first": 12, "last": 30, "label": "α1"},
    ]
    b = [
        {"id": "B:2-9", "chain": "B", "kind": "E", "first": 2, "last": 9, "label": "β1"},
        {"id": "B:14-30", "chain": "B", "kind": "H", "first": 14, "last": 30, "label": "α1"},
    ]
    assert r100.pair_check(entry(a), entry(b))["pass"] is True
    renamed = [dict(b[0]), dict(b[1], label="α2")]
    got = r100.pair_check(entry(a), entry(renamed))
    assert got["pass"] is False and "α1~α2" in got["detail"]
    assert r100.pair_check(entry(a), entry(b[:1]))["pass"] is False  # α1 has no partner


def test_report_md_lists_every_entry_and_the_limits(r100, tmp_path):
    def entry(i, ok):
        crit = {
            k: {"pass": True if ok else (False if k == "A6" else True), "detail": "" if ok else "default 2"}
            for k in r100.CRITERIA
        }
        e = {"id": i, "name": i, "stratum": "S1", "bin": "S", "views": {}, "facts": {"chains": 1}, "criteria": crit}
        return r100._verdict(e)

    entries = [entry("1ABC", True), entry("2DEF", False), entry("AF_AFP12345F1", True)]
    r100.report(entries, tmp_path, None, total=100)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "2 of 100 pass." in md
    for i in ("1ABC", "2DEF", "AF_AFP12345F1"):
        assert f"| {i} |" in md
    assert "A6: default 2" in md
    for k, v in r100.LIMITS.items():
        assert f"- {k}: {v}" in md
    doc = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert doc["limits"] == r100.LIMITS and len(doc["entries"]) == 3 and doc["foldmap_commit"]
    assert set(doc["versions"]) >= {"gemmi", "matplotlib", "numpy", "mkdssp", "rnapolis"}
    html = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert html.count("<img") == 3 and 'href="2DEF.png"' in html


def test_sheet_writes_a_png(r100, tmp_path):
    e = r100.run_entry(str(DATA / "1LMB.cif"), tmp_path, ("default",))
    drawn = r100.sheet(e, tmp_path)
    png = tmp_path / "1LMB.png"
    assert png.is_file() and png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert drawn["panels"] == 3 and drawn["protein_traces"] >= 2 and drawn["na_traces"] >= 2  # protein and DNA
    trace = json.loads((tmp_path / "work" / "1LMB" / "1LMB-trace.json").read_text(encoding="utf-8"))
    coloured = [c for ch in trace["chains"] for c in ch["colour"] if c]
    assert coloured and all(c.startswith("#") for c in coloured)  # the figure's element colours, by residue
    assert any(c for s in trace["strands"] for c in s["colour"])  # stems in their colours


def test_gitignore_keeps_images_local():
    import shutil

    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    for name in ("x.png", "x.svg", "index.html", "report-step3.md", "run-step2.log", "work/x/y.json", "1ABC.json"):
        done = subprocess.run(["git", "check-ignore", f"docs/validation/random100/{name}"], cwd=ROOT)
        assert done.returncode == 0, name
    for name in ("sample.json", "report.json", "report.md", "run.log", "review.md"):
        done = subprocess.run(["git", "check-ignore", f"docs/validation/random100/{name}"], cwd=ROOT)
        assert done.returncode == 1, name
