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

from foldmap import criteria

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "random100.py"
DATA = Path(__file__).parent / "data"
FID = DATA / "fidelity"


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


def _font_reason() -> str:
    """Why the ink checks are not run on this machine ('' on Linux CI, where text is drawn in Liberation Sans)."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    return figcheck.ink_font_reason(draw_spec(figure_spec(DATA / "1UBQ.cif")).figure)


def test_entry_measures_a_fixture(r100, tmp_path):
    """Where text is not drawn in Liberation Sans (macOS, Windows) the ink checks fail closed as not measured, and
    they are the only failures."""
    code = r100.main(["entry", str(DATA / "1UBQ.cif"), "--out", str(tmp_path), "--views", "quick"])
    saved = json.loads((tmp_path / "1UBQ.json").read_text(encoding="utf-8"))
    assert set(saved["criteria"]) == set(r100.CRITERIA) and {f"A{k}" for k in range(1, 21)} <= set(saved["criteria"])
    if _font_reason():
        assert code == 1 and saved["failed"] and all("not measured (ink checks" in f for f in saved["failed"])
        return
    assert code == 0
    assert saved["pass"] is True, saved["failed"]
    assert set(saved["views"]) == set(r100.QUICK)
    c = saved["criteria"]
    assert c["A1"]["pass"] and c["A11"]["pass"] and c["A12"]["pass"] and c["A13"]["pass"]
    assert c["A16"]["pass"] is None and c["A17"]["pass"] is None and c["A19"]["pass"] is None  # does not apply
    d = saved["views"]["default"]["result"]
    assert d["checks"]["element_overlaps"] == 0 and d["signature_sha256"] and d["svg_sha256"]
    assert d["route"]["records"] >= 1 and d["route"]["loops"] > 0  # route's summary record was read
    assert (tmp_path / "work" / "1UBQ" / "1UBQ-trace.json").is_file()


@pytest.mark.parametrize("platform", ["posix", "no wait4", "no wait4, no resource"])
def test_run_child_without_wait4(r100, tmp_path, monkeypatch, platform):
    """Windows has no os.wait4 (nor os.killpg, nor resource): the child is waited for with subprocess, its max RSS
    read from the children's rusage where there is one, else None; a hung child is still killed at the timeout."""
    if platform != "posix":
        monkeypatch.delattr(r100.os, "wait4", raising=False)
        monkeypatch.delattr(r100.os, "killpg", raising=False)
    if platform == "no wait4, no resource":
        monkeypatch.setitem(sys.modules, "resource", None)
    run = [sys.executable, "-I", "-c", "import sys; b = bytearray(64 << 20); sys.exit(3)"]
    got = r100.run_child(run, timeout=60, out=tmp_path / "o", err=tmp_path / "e")
    assert got["exit"] == 3 and not got["timed_out"] and got["seconds"] > 0
    try:
        import resource  # noqa: F401  (absent on Windows: no rusage, whatever the platform case)
    except ImportError:
        has_resource = False
    else:
        has_resource = True
    if platform == "no wait4, no resource" or not has_resource:
        assert got["max_rss_mb"] is None
    else:
        assert 60 < got["max_rss_mb"] < 4096
    hang = [sys.executable, "-I", "-c", "import time; time.sleep(60)"]
    got = r100.run_child(hang, timeout=0.5, out=tmp_path / "o", err=tmp_path / "e")
    assert got["timed_out"] and got["exit"] != 0 and got["seconds"] < 30
    e = {"views": {"default": {**got, "result": None}}}
    assert "default=" in " ".join(r100._metric_items(e)) and "max_rss_mb" in r100.reported(e)


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


def test_a17_gates_stem_names_at_the_margin_tolerance(r100):
    """A17 (ruling): a stem name may sit at most figcheck.stem_gap_tolerance(margin) from its own letters, the margin
    being layout._INK_MARGIN em in the figure's page units (measured per figure: stem_gap_limit), and nearer its
    own stem than any other."""
    view = {"exit": 0, "timed_out": False, "traceback": False, "python_warnings": [], "written": {}, "seconds": 1.0}

    def a17(gaps, limit=0.4):
        checks = {"details": {"stem_name_gaps": gaps}}
        result = {"checks": checks, "na": {"failed": []}, "stem_gap_limit": limit}
        return r100.criteria({"views": {"default": {**view, "result": result}}})["A17"]

    assert a17({"P1": [0.39, 3.0]})["pass"] is True  # 1.0 em out: beyond the old fixed 0.25, within the margin
    assert a17({"P1": [0.41, 3.0]})["pass"] is False
    assert a17({"P1": [0.3, 0.2]})["pass"] is False  # nearer another stem
    assert a17({"P1": [0.3, 3.0]}, limit=None)["pass"] is False  # names drawn but no limit measured: fail closed


def test_stem_gap_limit_is_the_margin_in_page_units(r100):
    """The limit read from a drawn stem name equals the one the layout tests read from a base letter."""
    from foldmap import figcheck
    from foldmap.layout import _INK_MARGIN
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(DATA / "na" / "1EHZ.cif.gz"))
    ax = d.figure.axes[0]
    letter = next(t for t in ax.texts if (t.get_gid() or "").startswith("na:nt:"))
    px_per_unit = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    em = letter.get_fontsize() / 0.75 * d.figure.dpi / 72 / px_per_unit
    assert criteria.stem_gap_limit(d.figure) == pytest.approx(figcheck.stem_gap_tolerance(_INK_MARGIN * em))
    assert 0.3 < criteria.stem_gap_limit(d.figure) < 1.0


def test_tether_lengths_are_measured_from_mark_to_target(r100):
    """A20 reads every ligand tether (they are patches, drawn in pieces where no body or text is; reading lines
    found none): myoglobin's heme and hydroxide, six tethers, the longest 6.06 page units."""
    from foldmap.spec import draw_spec, figure_spec
    from foldmap.style import Style

    d = draw_spec(figure_spec(DATA / "fidelity" / "1MBN.cif.gz"))
    got = criteria.tether_lengths(d, Style())
    assert len(got) == 6 and all(g.startswith("ligand-tether:") for g, *_ in got)
    assert max(n for _, n, _, _ in got) == pytest.approx(6.06, abs=0.01)


def test_a20_gates_only_tethers_to_other_chains(r100):
    """C2 (random-100): A20 bounds only a tether to another chain than the ligand's home (ruling); a long tether
    to the home chain is reported, not failed. Each tether records its home and target chains (myoglobin: one
    chain, so every tether is a home one)."""
    from foldmap.spec import draw_spec, figure_spec
    from foldmap.style import Style

    d = draw_spec(figure_spec(DATA / "fidelity" / "1MBN.cif.gz"))
    got = criteria.tether_lengths(d, Style())
    assert got and all(home == target == "A" for _, _, home, target in got)
    view = {"exit": 0, "timed_out": False, "traceback": False, "python_warnings": [], "written": {}, "seconds": 1.0}

    def a20(tethers):
        return r100.criteria({"views": {"default": {**view, "result": {"tethers": tethers}}}})["A20"]

    home = a20([["ligand-tether:URC:A300:0", 9.01, "A", "A"], ["ligand-tether:URC:A300:1", 2.0, "A", "A"]])
    assert home["pass"] is None and home["status"] == "not judged"  # no tether A20 bounds (final review B, I3)
    assert "reported: ligand-tether:URC:A300:0 9.01 (home A)" in home["detail"]
    other = a20([["ligand-tether:CLA:B609:2", 7.08, "B", "H"], ["ligand-tether:CLA:B609:0", 3.0, "B", "B"]])
    assert other["pass"] is False and other["detail"].startswith("ligand-tether:CLA:B609:2 7.08 (B to H)")
    assert a20([["ligand-tether:CLA:B609:2", 6.4, "B", "H"]])["pass"] is True
    # the detail names the first four of each kind and says how many more there are (C2b-5 fix round 1, M1)
    many = a20([[f"ligand-tether:CLA:A{k}:1", 9.0 + k, "A", "A"] for k in range(6)])
    assert many["pass"] is None and many["detail"].endswith("ligand-tether:CLA:A3:1 12.0 (home A) (+2 more)")
    crossed = a20([[f"ligand-tether:CLA:B{k}:3", 7.0 + k, "B", "H"] for k in range(5)])
    assert crossed["pass"] is False and crossed["detail"] == (
        "ligand-tether:CLA:B0:3 7.0 (B to H); ligand-tether:CLA:B1:3 8.0 (B to H); "
        "ligand-tether:CLA:B2:3 9.0 (B to H); ligand-tether:CLA:B3:3 10.0 (B to H); (+1 more)"
    )


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
    assert s["second_search_nodes"] == 0 and s["second_searches"] == 0  # ubiquitin routes within the first budget


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
    assert got["pass"] is False and "α1~α2" in got["detail"]  # nothing unpaired before it explains the name


def test_pair_check_allows_an_extra_element_and_its_numbering(r100):
    """A19 (ruling): one entry may hold an element the other lacks (8P3N's short helix B:121-125), and the generic
    names after it differ; enough elements must pair (the pairs' order is checked too; with each chain's elements
    apart in sequence, residue-range pairs cannot cross)."""

    def entry(elements):
        return {"views": {"beta": {"result": {"elements": elements}}}}

    def el(chain, kind, first, last, label):
        return {
            "id": f"{chain}:{first}-{last}",
            "chain": chain,
            "kind": kind,
            "first": first,
            "last": last,
            "label": label,
        }

    a = [el("A", "E", 1, 8, "β1"), el("A", "H", 12, 30, "α1"), el("A", "E", 34, 40, "β2"), el("A", "H", 44, 60, "α2")]
    b = [el("B", "E", 2, 9, "β1"), el("B", "H", 13, 30, "α1"), el("B", "E", 34, 41, "β2"), el("B", "H", 44, 50, "α2")]
    b += [el("B", "H", 52, 60, "α3")]  # the last helix split in two: one extra, numbering unchanged before it
    got = r100.pair_check(entry(a), entry(b))
    assert got["pass"] is True and "B:44-50 unpaired" in got["detail"]  # A:44-60 pairs the larger part
    extra = [el("B", "E", 2, 9, "β1"), el("B", "H", 10, 11, "α1")]  # an extra short helix: every later α shifts
    extra += [el("B", "H", 13, 30, "α2"), el("B", "E", 34, 41, "β2"), el("B", "H", 44, 60, "α3")]
    got = r100.pair_check(entry(a), entry(extra))
    assert got["pass"] is True and "B:10-11 unpaired" in got["detail"]
    assert r100.pair_check(entry(a), entry(b[:1]))["pass"] is False  # one of four paired: too few
    renamed = [dict(x) for x in extra]
    renamed[3]["label"] = "β3"  # a strand renamed with no extra strand before it: a real naming difference
    got = r100.pair_check(entry(a), entry(renamed))
    assert got["pass"] is False and "β2~β3 (A:34-40~B:34-41)" in got["detail"]  # final review B, Important 1


def test_pair_check_pairs_elements_by_sequence_not_numbers(r100):
    """C18 (random-100): 2HKI numbers GAPDH ahead of 1RM5 (by 1 to 3 residues), so pairing by author numbers paired
    14 of 28-30 elements. The beta view records the first chain's residues, and A19 aligns the two sequences and
    takes the second entry's element ranges to the first's numbering before pairing (ruling)."""
    import random as rnd

    names = rnd.Random("C18").choices(["ALA", "GLY", "SER", "LEU", "LYS", "GLU", "ASP", "PHE", "TRP", "VAL"], k=70)
    nums_a = list(range(1, 71))
    nums_b = list(range(1, 21)) + list(range(27, 77))  # b numbers six ahead after residue 20

    def el(chain, kind, first, last, label):
        return {"id": f"{chain}:{first}-{last}", "chain": chain, "kind": kind, "first": first, "last": last}, label

    def entry(chain, nums, elements):
        els = [dict(e, label=label) for e, label in elements]
        res = {"chain": chain, "residues": [[n, x] for n, x in zip(nums, names)]}
        return {"views": {"beta": {"result": {"elements": els, "chain_residues": res}}}}

    a = entry("A", nums_a, [el("A", "E", 5, 10, "A"), el("A", "H", 25, 35, "α1"), el("A", "E", 40, 44, "B")])
    a["views"]["beta"]["result"]["elements"].append(el("A", "H", 50, 58, "α2")[0] | {"label": "α2"})
    b = entry("B", nums_b, [el("B", "E", 5, 10, "A"), el("B", "H", 31, 41, "α1"), el("B", "E", 46, 50, "B")])
    b["views"]["beta"]["result"]["elements"].append(el("B", "H", 56, 64, "α2")[0] | {"label": "α2"})
    got = r100.pair_check(a, b)
    assert got["pass"] is True and "unpaired" not in got["detail"], got
    for e in (a, b):  # without the residues (an older report) the numbers are read as before
        del e["views"]["beta"]["result"]["chain_residues"]
    assert r100.pair_check(a, b)["pass"] is False


def test_beta_view_records_the_first_chain_residues(r100, tmp_path):
    """C18: the beta view's measure carries the residues of its first element's chain (number, name) for A19."""
    from foldmap.spec import draw_spec, figure_spec

    spec = figure_spec(DATA / "1UBQ.cif", names=["beta"])
    got = r100.measure_view("beta", draw_spec(spec), spec, DATA / "1UBQ.cif", "1UBQ", tmp_path)
    res = got["chain_residues"]
    assert res["chain"] == "A" and res["residues"][:2] == [[1, "MET"], [2, "GLN"]] and len(res["residues"]) == 76


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


def test_stem_gap_limit_reads_a_block_name_at_its_own_size(r100):
    """A block-style stem's name is drawn at 0.9 em, a ladder's at 0.8 (render): A17's em comes out the figure's em
    either way, the size an element's label is drawn at (4YBB's crop draws blocks)."""
    from foldmap import figcheck
    from foldmap.layout import _INK_MARGIN
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(DATA / "na" / "4YBB.cif.gz"))
    ax = d.figure.axes[0]
    assert any((a.get_gid() or "").endswith(":block") for a in ax.patches)
    label = next(t for t in ax.texts if (t.get_gid() or "").startswith("label:"))  # drawn at one em
    px = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    em = label.get_fontsize() * d.figure.dpi / 72 / px
    assert criteria.stem_gap_limit(d.figure) == pytest.approx(figcheck.stem_gap_tolerance(_INK_MARGIN * em))


def test_a_chain_with_nothing_to_draw_does_not_widen_the_focus_view(r100):
    """C15 fix round 1: "not drawn: chain P" is part of the scope's reason but no scope the default view left out of
    the chains, so 7ARD and 7PIU keep their single-chain focus view (area 538, not 18,984)."""
    r = r100.focus_target
    assert (
        r({"scope": {"reason": "not drawn: chain P", "chains_reason": "", "crop": ""}, "first_protein_chain": "A"})
        == "A"
    )
    assert (
        r(
            {
                "scope": {"reason": "x", "chains_reason": "assembly 1, one of 2 copies", "crop": ""},
                "first_protein_chain": "A",
            }
        )
        == "none"
    )
    assert (
        r({"scope": {"reason": "x", "chains_reason": "", "crop": "RNA cropped"}, "first_protein_chain": "A"}) == "none"
    )
    assert r({"scope": {}, "first_protein_chain": ""}) == "none"


def test_the_focus_view_draws_the_chain_with_most_elements(r100):
    """C25 (random-100): with nothing left out, the focus view drew the first protein chain, whatever it was: 7JU4's
    chain "0" is a two-helix coiled coil. It draws the protein chain with the most elements in the default view
    (the first of them on a tie), else the first protein chain."""
    els = [{"chain": "0"}, {"chain": "0"}, {"chain": "A"}, {"chain": "A"}, {"chain": "A"}, {"chain": "B"}]
    assert r100.focus_target({"scope": {}, "first_protein_chain": "0", "elements": els}) == "A"
    assert r100.focus_target({"scope": {}, "first_protein_chain": "0", "elements": els[:4]}) == "0"
    assert r100.focus_target({"scope": {}, "first_protein_chain": "0"}) == "0"


def test_a18_tm_helix_passes_when_it_reaches_both_leaflets(r100):
    """C27 (random-100, ruling): a TM helix passes A18 when it enters the band through a face and reaches past its
    middle, not by covering 75 % of the band's height: a tilted or offset helix that
    crosses one face and reaches past the middle (6PS3 A:267-298 covers 0.73, 7PIU R:48-70 0.70) passes; one
    that stays in one leaflet fails. Bacteriorhodopsin's seven helices, crossing the membrane (membrane.crosses, plan
    26), pass as drawn."""
    from dataclasses import replace

    import numpy as np

    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "1C3W.cif"
    d = draw_spec(figure_spec(path, membrane="auto", assembly="asu"))
    m = criteria.membrane_checks(d, path)
    assert m["band"] and m["tm_helices"] == 7 and m["tm_short"] == [], m
    lay = d.layout
    (y0, y1), mid = lay.membrane["y"], sum(lay.membrane["y"]) / 2
    hid = next(k for k, p in lay.placed.items() if abs(np.sin(p.angle)) > 0.9 and p.length > y1 - y0)
    p = lay.placed[hid]
    half = abs(p.c_port[1] - p.n_port[1]) / 2
    lay.placed[hid] = replace(p, cy=y0 + 0.7 * (y1 - y0) - half)  # crosses the lower face, ends at 0.70: spans
    assert criteria.membrane_checks(d, path)["tm_short"] == []
    lay.placed[hid] = replace(p, cy=mid - 0.05 - half)  # crosses the lower face, ends short of the middle
    assert [x.split()[0] for x in criteria.membrane_checks(d, path)["tm_short"]] == [hid]
    lay.placed[hid] = replace(p, cy=mid + 0.05 + half)  # the upper leaflet only
    assert [x.split()[0] for x in criteria.membrane_checks(d, path)["tm_short"]] == [hid]
    short = (y1 - y0) / 3  # a short helix centred in the band, entering it through neither face (fix round 1)
    lay.placed[hid] = replace(p, cy=mid, length=short * abs(p.length) / max(2 * half, 1e-9))
    assert [x.split()[0] for x in criteria.membrane_checks(d, path)["tm_short"]] == [hid]


def test_science_checks_are_reported_apart_from_the_criteria(r100, tmp_path):
    """Plan 12 Task 1: the default view's science checks (S1-S5) get their own section in report.md and never make an
    entry fail; the measure of a figure carries them."""
    from foldmap.spec import draw_spec, figure_spec

    spec = figure_spec(DATA / "1UBQ.cif")
    got = r100.measure_view("default", draw_spec(spec), spec, DATA / "1UBQ.cif", "1UBQ", tmp_path)
    assert got["science"] == {"S1": [], "S2": [], "S3": [], "S4": [], "S5": None, "S6": []}
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}
    found = {"S1": ["sheet 0: a | b drawn side by side with no bridge"], "S2": [], "S3": [], "S4": None, "S5": None}
    e = {"id": "1ABC", "name": "1ABC", "views": {"default": {"result": {"science": found}}}, "criteria": crit}
    e = r100._verdict(e)
    assert e["pass"] and r100.science(e)["S1"] == {
        "pass": False,
        "detail": "default: " + found["S1"][0],
        "kind": "consistency",
    }
    assert r100.science(e)["S4"]["pass"] is None
    r100.report([e], tmp_path, None, total=1)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "1 of 1 pass." in md and "## Science checks (S1-S6), apart from A1-A20" in md
    assert (
        "| S1 sheet neighbours pair | 0 | 1 | 0 | 0 |" in md and "default: S1 sheet 0: a / b drawn side by side" in md
    )
    # plan 28 Task 2 (ruling Q19): each check says what it is evidence of, in report.md and report.json
    assert "| S1 sheet neighbours pair | 0 | 1 | 0 | 0 | consistency |" in md
    assert "| S5 automatic domains agree with a reference | 0 | 0 | 0 | 1 | external |" in md
    assert r100.science(e)["S5"]["kind"] == "external"
    kinds = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))["check_kinds"]
    assert kinds["S5"] == "external" and kinds["S6"] == "consistency" and kinds["A10"] == "readability"
    assert set(kinds) == {"S1", "S2", "S3", "S4", "S5", "S6", *r100.CRITERIA}


def test_science_checks_run_in_every_drawn_view(r100, tmp_path):
    """Batch C2b-S: the science checks run on every drawn view (default, deposited, beta, focus, domains, membrane,
    and the app's figure), not the default view alone: the code that ties, regroups and shares copies' sheets is
    drawn only in some views (9P1X focus). An entry fails S_k when any view has a finding; each finding names its
    view; report.md adds a breakdown per view. S5 stays once per entry (chain-level, against CATH)."""
    from foldmap.spec import draw_spec, figure_spec

    clean = {"S1": [], "S2": [], "S3": [], "S4": [], "S5": None, "S6": []}
    for view, kw in (("focus", {"focus": "A"}), ("beta", {"names": ["beta"]})):
        spec = figure_spec(DATA / "1UBQ.cif", **kw)
        got = r100.measure_view(view, draw_spec(spec), spec, DATA / "1UBQ.cif", "1UBQ", tmp_path)
        assert got["science"] == clean, view
    assert r100.app_view(DATA / "1UBQ.cif", "1UBQ")["science"] == clean
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}
    focus = {**clean, "S1": ["sheet 27: a | b drawn side by side with no bridge"]}
    domains = {**clean, "S6": ["sheet 3: c drawn away from it, not beside its partner d"]}
    views = {
        "default": {"result": {"science": clean}},
        "focus": {"result": {"science": focus}},
        "domains": {"result": {"science": domains}},
        "app": {"result": {"app": {"refused": None}, "science": clean}},
        "seed1": {"result": {}},  # not judged for science (A12 holds it to the default view)
    }
    e = r100._verdict({"id": "1ABC", "name": "1ABC", "views": views, "criteria": crit})
    s = r100.science(e)
    assert e["pass"] and s["S1"] == {
        "pass": False,
        "detail": "focus: sheet 27: a | b drawn side by side with no bridge",
        "kind": "consistency",
    }
    assert s["S6"]["pass"] is False and s["S6"]["detail"].startswith("domains: sheet 3: c drawn away")
    assert s["S2"]["pass"] is True and s["S4"]["pass"] is True and s["S5"]["pass"] is None
    assert r100.science_views(e) == {"default": {}, "focus": {"S1": 1}, "domains": {"S6": 1}, "app": {}}
    none = r100.science({"id": "2ABC", "views": {"default": {"result": None}}})
    assert none["S1"]["pass"] is None
    plain = r100._verdict(
        {"id": "1UBQ", "name": "1UBQ", "views": {"default": {"result": {"science": clean}}}, "criteria": dict(crit)}
    )
    r100.report([e, plain], tmp_path, None, total=2)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "2 of 2 pass." in md and "| S1 sheet neighbours pair | 1 | 1 | 0 | 0 |" in md
    assert "### Science by view" in md
    assert "| focus | 1 | 1 | 0 | 0 | 0 | 0 |" in md and "| default | 2 | 0 | 0 | 0 | 0 | 0 |" in md
    assert "| domains | 1 | 0 | 0 | 0 | 0 | 1 |" in md
    assert "focus: S1 sheet 27: a / b drawn side by side" in md and "domains: S6 sheet 3: c drawn away" in md


def test_the_app_views_science_is_measured_on_the_figure_the_app_drew(r100, monkeypatch):
    """C2b-S fix round 1 (review minor 2): the app view's science is measured on a second drawing of the app's spec;
    its SVG must equal the one the app's job returned, or the view is not measured (measure_error, which fails A1)."""
    got = r100.app_view(DATA / "1UBQ.cif", "1UBQ")
    assert "measure_error" not in got and got["science"]["S1"] == []
    with pytest.raises(ValueError, match="differs from the figure the app drew"):
        r100.app_science(DATA / "1UBQ.cif", "<svg/>")
    real = r100.app_science
    monkeypatch.setattr(r100, "app_science", lambda path, svg: real(path, svg + " "))
    got = r100.app_view(DATA / "1UBQ.cif", "1UBQ")
    assert "science" not in got and "differs from the figure the app drew" in got["measure_error"]


def test_accessories_named_are_reported_per_view_and_never_gate(r100, tmp_path):
    """Plan 12 ruling (fix round 4): the count of accessories named rather than drawn level beside their partner is
    reported per entry and view (OmpF 3 of 3, 3L9G 8 of 16), apart from A1-A20: it never makes an entry fail."""
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "r100" / "2OMF-A.cif.gz"
    spec = figure_spec(path)
    got = r100.measure_view("default", draw_spec(spec), spec, path, "2OMF-A", tmp_path)
    assert (got["checks"]["accessories_named"], got["checks"]["accessories"]) == (3, 3)
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}

    def view(named, of):
        return {"result": {"checks": {"accessories_named": named, "accessories": of}}}

    e = {"id": "3L9G", "name": "3L9G", "views": {"default": view(8, 16), "focus": view(6, 16)}, "criteria": crit}
    e = r100._verdict(e)
    assert e["pass"] and r100.accessories_named(e) == {"default": (8, 16), "focus": (6, 16)}
    plain = r100._verdict({"id": "1UBQ", "name": "1UBQ", "views": {"default": view(0, 0)}, "criteria": dict(crit)})
    r100.report([e, plain], tmp_path, None, total=2)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "2 of 2 pass." in md and "## Accessories named, not drawn level beside their partner" in md
    assert "| 3L9G | default 8 of 16; focus 6 of 16 |" in md and "| 1UBQ |" not in md.split("## Accessories")[1]


def test_s5_scores_auto_domains_where_a_cath_reference_exists(r100, tmp_path):
    """Plan 12 Task 3: S5 measures an entry's automatic domains against its CATH chopping in the committed reference
    (tests/data/domains/cath_reference.json), per chain: overlap and boundary F1; a chain under the ruled targets is a
    finding. Entries without a reference are not measured. Reported, never gating."""
    from foldmap.domain_score import TARGETS

    got, note, _ = r100.domain_reference(FID / "1LDM.cif.gz", "1LDM", fetch=False)
    assert [g["chain"] for g in got] == ["A"] and got[0]["overlap"] > 0.9 and got[0]["n_ref"] == 2
    assert got[0]["source"] == "reference" and note == "a local file: CATH not fetched"
    assert r100.domain_reference(DATA / "1UBQ.cif", "1UBQ", fetch=False)[0] is None
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}
    found = {"S1": [], "S2": [], "S3": [], "S4": [], "S5": None, "S6": []}
    low = {"pdb": "1ABC", "chain": "A", "overlap": 0.6, "f1": 0.5, "recall": 1.0, "precision": 0.33, "n_ref": 2}
    e = {"id": "1ABC", "name": "1ABC", "views": {"default": {"result": {"science": found}}}, "criteria": crit}
    e = r100._verdict({**e, "domains_cath": [{**low, "n_auto": 4}]})
    assert e["pass"] and r100.science(e)["S5"]["pass"] is False
    assert "A: overlap 0.60, boundary F1 0.50, recall 1.00" in r100.science(e)["S5"]["detail"]
    ok = {**low, "overlap": TARGETS["overlap"], "f1": TARGETS["f1"], "n_auto": 2}
    assert r100.science({**e, "domains_cath": [ok]})["S5"]["pass"] is True
    assert r100.science({**e, "domains_cath": None})["S5"]["pass"] is None
    missed = {**ok, "recall": 0.5}  # a CATH boundary missed: recall under .85 is a finding too
    assert "recall 0.50" in r100.science({**e, "domains_cath": [missed]})["S5"]["detail"]
    r100.report([e], tmp_path, None, total=1)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "## S5: automatic domains against CATH (reported, not gated)" in md
    assert "the 33 chains the automatic domains were calibrated on" in md
    assert "| 1ABC | A | reference | 2 | 4 | 0.60 | 0.50 | 1.00 | 0.33 |" in md  # a report from before: the reference


# ------------------------------------------------------------------------------------- final review B (fix F1)

VIEW = {"exit": 0, "timed_out": False, "traceback": False, "python_warnings": [], "written": {}, "seconds": 1.0}


def _judged(r100, views: dict, **entry) -> dict:
    return r100.criteria({"views": {v: {**VIEW, "result": res} for v, res in views.items()}, **entry})


def test_nothing_to_judge_is_never_a_pass(r100):
    """Final review B, Important 3: a criterion with nothing to judge reads not judged (it applies, the entry gave it
    nothing) or does not apply, never pass: A15 with no tie judged, A20 with no tether to another chain, A16 on a
    domains view that drew no domain, A10 unmeasured, A2b refused over the cap."""

    def a15(ties):
        return _judged(r100, {"default": {"checks": {"copies_alike": 0.0, "details": {"copies_alike": ties}}}})["A15"]

    assert a15({}) == {"pass": None, "detail": "no copies to compare", "status": "does not apply"}
    got = a15({"B~A": {"deviation": None, "matched": 2, "of": 2, "judged": "n/a: 2 elements"}})
    assert got["pass"] is None and got["status"] == "not judged" and "n/a: 2 elements" in got["detail"]
    assert a15({"B~A": {"deviation": 0.1, "matched": 9, "of": 9, "judged": "yes"}})["pass"] is True

    def a20(tethers):
        return _judged(r100, {"default": {"tethers": tethers}})["A20"]

    assert a20([]) == {"pass": None, "detail": "no ligand tether", "status": "does not apply"}
    home = a20([["ligand-tether:HEM:A1:0", 9.0, "A", "A"]])
    assert (
        home["pass"] is None and home["status"] == "not judged" and "reported: ligand-tether:HEM:A1:0" in home["detail"]
    )
    assert a20([["ligand-tether:HEM:A1:0", 3.0, "A", "B"]])["pass"] is True

    counts = {"domain_rules": 0, "domain_overlap": 0.0, "domain_name_collisions": 0}
    none = _judged(r100, {"domains": {"checks": {**counts, "details": {"domains": {}}}}})["A16"]
    assert none["pass"] is None and none["status"] == "not judged"
    drew = {"cross_chain": [], "unalike": [], "small": [], "panel_overlaps": [], "name_collisions": []}
    assert _judged(r100, {"domains": {"checks": {**counts, "details": {"domains": drew}}}})["A16"]["pass"] is True

    a10 = _judged(r100, {}, dssp={"agree": None, "unmeasured": "no mkdssp on PATH"})["A10"]
    assert a10 == {"pass": None, "detail": "unmeasured: no mkdssp on PATH", "status": "not judged"}
    app = {"app": {"refused": "size", "residues": 31_000, "message": "too big"}}
    assert _judged(r100, {"app": app})["A2b"]["status"] == "does not apply"
    assert [criteria.outcome(x) for x in (a10, {"pass": True}, {"pass": False}, {"pass": None})] == [
        "not judged",
        "pass",
        "fail",
        "does not apply",
    ]


def test_a_count_that_is_not_a_number_fails(r100):
    """Final review B, Minor 8: a NaN count fails (NaN > limit is False); jsonable turns a numpy NaN into None."""
    import numpy as np

    assert _judged(r100, {"default": {"checks": {"element_overlaps": float("nan")}}})["A5"]["pass"] is False
    assert criteria.jsonable({"x": np.float64("nan"), "y": np.float64(1.5)}) == {"x": None, "y": 1.5}


def test_pair_check_compares_names_along_the_whole_alignment(r100):
    """Final review B, Important 1: names are compared on every pair, each number less the unpaired elements of its
    kind before it (6YPX/6YQD: an extra η before η1~η2, an extra α before α2~α1), so a pair whose first element is
    unpaired is still compared; a different naming (9YKJ's αN~α2 against 10UJ's generic α2) fails."""

    def entry(els):
        return {"views": {"beta": {"result": {"elements": els}}}}

    def el(chain, kind, first, last, label):
        return {
            "id": f"{chain}:{first}-{last}",
            "chain": chain,
            "kind": kind,
            "first": first,
            "last": last,
            "label": label,
        }

    a = [
        el("A", "H", 55, 60, "α1"),
        el("A", "E", 68, 71, "β1"),
        el("A", "E", 75, 79, "β2"),
        el("A", "G", 100, 102, "η1"),
    ]
    a += [el("A", "H", 105, 124, "α2"), el("A", "E", 131, 137, "β3")]
    b = [
        el("B", "G", 65, 67, "η1"),
        el("B", "E", 68, 71, "β1"),
        el("B", "E", 75, 79, "β2"),
        el("B", "G", 100, 102, "η2"),
    ]
    b += [el("B", "H", 105, 124, "α1"), el("B", "E", 131, 137, "β3")]
    got = r100.pair_check(entry(a), entry(b))
    assert got["pass"] is True and got["detail"].startswith("names compared on all 5 pairs"), got
    other = [dict(x) for x in b]
    other[4]["label"] = "αN"  # another naming of the same helix
    got = r100.pair_check(entry(a), entry(other))
    assert got["pass"] is False and "α2~αN (A:105-124~B:105-124)" in got["detail"], got
    assert criteria.same_name("α10a", "α11a", [], ["α3"]) and not criteria.same_name("α10a", "α10b", [], [])


def test_a17_is_judged_on_the_view_that_draws_the_nucleic_acid(r100):
    """Final review B, Important 5: 8QEM's default view draws one subunit of its C12 ring and no DNA; its focus
    view (--focus none) draws the DNA, so A17 is measured there (the focus view's process is asked to: --na). An entry
    with nucleic acid that no view measured fails when a view drew it, and is not judged when none did."""
    na = {"na": {"failed": ["pair_f1 0.5"]}, "checks": {"details": {"stem_name_gaps": {}}}, "stem_gap_limit": 0.4}
    got = _judged(r100, {"default": {}, "focus": na})["A17"]
    assert (
        got["pass"] is False and got["detail"].startswith("measured on the focus view") and "pair_f1" in got["detail"]
    )
    facts = {"na": True, "protein": True}
    drawn = _judged(r100, {"default": {}, "focus": {"na_chains_drawn": ["Y"]}}, facts=facts)["A17"]
    assert drawn["pass"] is False and "measured in none" in drawn["detail"]
    assert _judged(r100, {"default": {}}, facts=facts)["A17"]["status"] == "not judged"
    cmd, _ = r100.view_command("focus", Path("x.cif"), "x", Path("w"), "none", na=True)
    assert "--na" in cmd and cmd.index("--na") < cmd.index("focus")


def test_science_with_nothing_to_judge_is_not_judged(r100):
    """Final review B, Important 3: a science check passes only where a view had something to judge (science_judged:
    sheets, barrels, marks, accessories); with protein and nothing it is not judged, with no protein it does not
    apply."""
    clean = {"S1": [], "S2": [], "S3": [], "S4": [], "S5": None, "S6": []}
    judged = {"S1": 1, "S2": 0, "S3": 1, "S4": 0, "S6": 0}
    e = {"facts": {"protein": True}, "views": {"default": {"result": {"science": clean, "science_judged": judged}}}}
    got = r100.science(e)
    assert got["S1"]["pass"] is True and got["S3"]["pass"] is True
    assert {k: got[k].get("status") for k in ("S2", "S4", "S5", "S6")} == dict.fromkeys(
        ("S2", "S4", "S5", "S6"), "not judged"
    )
    assert got["S2"]["detail"] == "no closed barrel drawn in any view"
    rna = {
        "facts": {"protein": False},
        "views": {"default": {"result": {"science": clean, "science_judged": dict.fromkeys(judged, 0)}}},
    }
    assert {v.get("status") for v in r100.science(rna).values()} == {"does not apply"}


def test_report_md_counts_each_outcome_per_criterion(r100, tmp_path):
    """report.md gives every criterion's pass / fail / not judged / does not apply, and lists each entry a criterion did
    not judge, with why (A10's unmeasured, as the plan asks)."""
    crit = {k: {"pass": None, "detail": "", "status": "does not apply"} for k in r100.CRITERIA}
    crit["A5"] = {"pass": True, "detail": ""}
    crit["A10"] = {"pass": None, "detail": "unmeasured: no mkdssp on PATH", "status": "not judged"}
    e = r100._verdict({"id": "1ABC", "name": "1ABC", "stratum": "S1", "views": {}, "criteria": crit})
    r100.report([e], tmp_path, None, total=1)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "| criterion | pass | fail | not judged | does not apply |" in md
    assert "| A5 no element or stem boxes overlapping | 1 | 0 | 0 | 0 |" in md
    assert "| A10 secondary structure agrees with mkdssp | 0 | 0 | 1 | 0 |" in md
    assert "- **A10** (1): 1ABC (unmeasured: no mkdssp on PATH)" in md


def test_additives_are_no_ligand(r100, tmp_path):
    """Final review B, Important 4: the ligands the sampler counts (and apo/holo reads) leave out the additives
    foldmap hides (io.ADDITIVES): 6YPX was holo by glycerol alone. The same residue under another name is one."""
    import gemmi

    def with_residue(name):
        st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
        st.remove_waters()
        chain = gemmi.Chain("Z")
        res = gemmi.Residue()
        res.name, res.seqid, res.het_flag = name, gemmi.SeqId(1, " "), "H"
        ca = st[0]["A"][0]["CA"][0].pos
        for k, el in enumerate("CCCOOO"):
            atom = gemmi.Atom()
            atom.name, atom.element = f"{el}{k}", gemmi.Element(el)
            atom.pos = gemmi.Position(ca.x + 2.5 + 0.1 * k, ca.y, ca.z)
            res.add_atom(atom)
        chain.add_residue(res)
        st[0].add_chain(chain)
        st.setup_entities()
        path = tmp_path / f"{name}.cif"
        st.make_mmcif_document().write_file(str(path))
        return r100.features(path, contacts=True)

    gol = with_residue("GOL")
    assert gol["ligands"] == [] and "ligand" not in gol["has"] and gol["ligand_touching"] is False
    other = with_residue("XYZ")
    assert other["ligands"] == ["XYZ"] and "ligand" in other["has"] and other["ligand_touching"] is True


def test_ca_rmsd_aligns_residues_by_sequence(r100, monkeypatch):
    """Final review B, Minor 9: the open/closed RMSD pairs residues by sequence, not by author number (C18)."""
    import gemmi
    import numpy as np

    names = random.Random("rmsd").choices(["ALA", "GLY", "SER", "LEU", "LYS", "TRP"], k=40)
    a = [(n, gemmi.Position(1.5 * k, 2.3 * np.cos(k), 2.3 * np.sin(k))) for k, n in enumerate(names)]
    b = [("GLY", gemmi.Position(0, 50, 0))] * 3 + a  # three extra residues ahead: numbers shifted by three
    monkeypatch.setattr(r100, "ca_trace", lambda path: a if path == "a" else b)
    assert r100.ca_rmsd("a", "b") == pytest.approx(0.0, abs=1e-6)


def test_the_sampler_never_draws_an_excluded_id(r100):
    """The held-out sample (final review A, Important 2) draws none of the tuned sample's ids."""
    s = r100.Sampler(20261017, "2026-10-13", exclude=["1ABC"])
    s.feats["1ABC"] = s.feats["2DEF"] = {"protein": True, "na": False, "ca_only": False}
    assert s.check("1ABC") == (False, "in the excluded sample") and s.check("2DEF") == (True, "")


def test_sample_excludes_the_ids_of_every_sample_named(r100, tmp_path, monkeypatch):
    """Plan 15 Task 9: the third sample draws none of the tuned or the held-out ids (--exclude given twice), and a
    CA-only replacement of either counts as one of its ids."""
    tuned, held = tmp_path / "tuned.json", tmp_path / "held.json"
    tuned.write_text(json.dumps({"entries": [{"id": "1ABC"}, {"id": "2DEF"}], "replaced": []}))
    held.write_text(
        json.dumps({"entries": [{"id": "3GHI"}, {"id": "1ABC"}], "replaced": [{"id": "4JKL", "replacement": None}]})
    )
    seen = {}

    def fake(seed, date, exclude=()):
        seen["exclude"] = list(exclude)
        return {"entries": [], "quotas": {}}

    monkeypatch.setattr(r100, "sample", fake)
    out = tmp_path / "third.json"
    r100.main(["sample", "--date", "2026-10-13", "--exclude", str(tuned), "--exclude", str(held), "--out", str(out)])
    assert seen["exclude"] == ["1ABC", "2DEF", "3GHI", "4JKL"]
    r100.main(["sample", "--date", "2026-10-13", "--out", str(out)])
    assert seen["exclude"] == []


def test_an_entrys_work_folder_starts_empty(r100, tmp_path):
    """Final review B, Minor 5: a view that writes nothing is never judged on an earlier run's files."""
    work = tmp_path / "work" / "1ABC"
    work.mkdir(parents=True)
    for name in ("1ABC-default.svg", "1ABC-default.json", "entry.out", "entry.err"):
        (work / name).write_text("old")
    r100.clear_work(work)
    assert sorted(p.name for p in work.iterdir()) == ["entry.err", "entry.out"]


def test_the_commit_is_never_taken_from_the_environment(r100, monkeypatch):
    """Plan 28 ruling Q11 (replacing final review B, Minor 6's FOLDMAP_COMMIT): foldmap_commit is git's or "unknown";
    a runner's asserted commit goes into the source stamp, marked as asserted."""
    monkeypatch.setenv("FOLDMAP_COMMIT", "abc123")
    assert r100._commit() != "abc123"
    git = r100.SOURCE_STAMP["git"] or {}
    assert r100._commit() == (git["commit"] if git.get("source") == "git" else "unknown")


def test_a_gzipped_report_is_read(r100, tmp_path):
    """A report.json over 5 MB is committed gzipped; --failed-from reads it either way."""
    import gzip

    doc = {"entries": [{"id": "1ABC", "pass": False}]}
    (tmp_path / "report.json.gz").write_bytes(gzip.compress(json.dumps(doc).encode()))
    assert r100.read_report(tmp_path / "report.json") == doc == r100.read_report(tmp_path / "report.json.gz")


# ------------------------------------------------------------------------------------- F1 review, fix round 1


def _pair(els_a, els_b, seq_a=None, seq_b=None):
    def entry(chain, els, seq):
        result = {"elements": els}
        if seq is not None:
            result["chain_residues"] = {"chain": chain, "residues": [[k + 1, n] for k, n in enumerate(seq)]}
        return {"views": {"beta": {"result": result}}}

    return entry("A", els_a, seq_a), entry("B", els_b, seq_b)


def _el(chain, kind, first, last, label):
    return {"id": f"{chain}:{first}-{last}", "chain": chain, "kind": kind, "first": first, "last": last, "label": label}


def test_pair_check_equal_names_always_match(r100):
    """F1 review, Important 1 (2BBJ/2BBH): η2~η2 with an unpaired η1 drawn on one side only (the other's η1 is not
    among its drawn elements) is the same name: equal names pass before any shift."""
    a = [_el("A", "E", 10, 15, "β1"), _el("A", "E", 20, 25, "β2"), _el("A", "G", 52, 55, "η2")]
    a += [_el("A", "H", 60, 80, "α1")]
    b = [_el("B", "E", 10, 15, "β1"), _el("B", "E", 20, 25, "β2"), _el("B", "G", 47, 49, "η1")]
    b += [_el("B", "G", 51, 55, "η2"), _el("B", "H", 60, 80, "α1")]
    got = r100.pair_check(*_pair(a, b))
    assert got["pass"] is True, got
    assert criteria.same_name("η2", "η2", [], ["η1"])


def test_a_pair_of_two_proteins_is_not_judged_as_one(r100):
    """F1 review, Important 2 (held-out P3, 8U9T/7IT6: nsp5 and nsp1 of one polyprotein accession, P0DTD1): when the
    two first chains are not one protein by sequence (identity over the aligned residues under PAIR_IDENTITY), A19
    does not apply; a pair of one protein is judged."""
    rng = random.Random("nsp")
    aa = ["ALA", "GLY", "SER", "LEU", "LYS", "GLU", "ASP", "PHE", "TRP", "VAL", "ILE", "PRO", "THR", "ASN"]
    one, other = rng.choices(aa, k=120), rng.choices(aa, k=110)
    a = [_el("A", "E", 5, 10, "β1"), _el("A", "H", 20, 40, "α1"), _el("A", "E", 50, 55, "β2")]
    b = [_el("B", "H", 3, 30, "α1"), _el("B", "E", 40, 45, "β1"), _el("B", "H", 60, 90, "α2")]
    got = r100.pair_check(*_pair(a, b, one, other))
    assert got["pass"] is None and got["status"] == "does not apply" and got["detail"].startswith("not one protein")
    same = r100.pair_check(*_pair(a, [dict(x, chain="B", id="B" + x["id"][1:]) for x in a], one, one))
    assert same["pass"] is True


def test_a11_and_a8_with_nothing_drawn_are_not_passed(r100):
    """F1 review, Important 4: A11 (N->C and sheet order) on a default view with no protein element does not apply
    when the figure has no protein (8Q4O, 3OWW) and is not judged when its protein has no element drawn; A8 with no
    element and no nucleic acid drawn in any label view is not judged."""
    empty = {"checks": {"order": 0, "label_repeats": 0}, "elements": [], "na_chains_drawn": []}
    rna = {**empty, "na_chains_drawn": ["A"], "first_protein_chain": None}
    got = _judged(r100, {"default": rna})
    assert got["A11"]["status"] == "does not apply" and got["A8"]["pass"] is True  # its stem names are judged
    coil = {**empty, "first_protein_chain": "A"}
    got = _judged(r100, {"default": coil})
    assert got["A11"]["status"] == "not judged" and got["A8"]["status"] == "not judged"
    drawn = {**coil, "elements": [{"id": "A:1-5", "chain": "A", "kind": "E", "first": 1, "last": 5, "label": "β1"}]}
    assert _judged(r100, {"default": drawn})["A11"]["pass"] is True


def test_a16_names_its_domain_name_collisions(r100):
    """F1 review, Minor 6 (7JU4 "domains 2 ()"): A16's detail names what failed."""
    counts = {"domain_rules": 0, "domain_overlap": 0.0, "domain_name_collisions": 1}
    hit = "domain-label:D2 · 4 × stub:C:4"
    doms = {"cross_chain": [], "unalike": [], "small": [], "panel_overlaps": [], "name_collisions": [hit]}
    got = _judged(r100, {"domains": {"checks": {**counts, "details": {"domains": doms}}}})["A16"]
    assert got["pass"] is False and hit in got["detail"], got


# ------------------------------------------------------------------------------------------ plan 15 Task 8


def test_the_sampler_pairs_one_protein_not_one_accession(r100, monkeypatch):
    """Plan 15 Task 8 (held-out P3: two proteins of one polyprotein accession, nsp5 with nsp1, drawn as an apo/holo
    pair): Sampler.pairs keeps a pair only when its two first protein chains are one protein by sequence
    (criteria.one_protein, the rule A19 judges by), not by UniProt accession alone. A polyprotein-shaped fixture:
    three entries of one accession, the seed and one partner the same protein, the other partner another protein cut
    from the same polyprotein (ligand-bound, so it would be an apo/holo partner)."""
    import gemmi
    import numpy as np

    rng = random.Random("polyprotein")
    aa = ["ALA", "GLY", "SER", "LEU", "LYS", "GLU", "ASP", "PHE", "TRP", "VAL", "ILE", "PRO", "THR", "ASN"]
    first, second = rng.choices(aa, k=120), rng.choices(aa, k=110)  # two proteins of one polyprotein

    def trace(names):
        return [(n, gemmi.Position(1.5 * k, 2.3 * np.cos(k), 2.3 * np.sin(k))) for k, n in enumerate(names)]

    traces = {"1AAA": trace(first), "2BBB": trace(second), "3CCC": trace(first)}
    monkeypatch.setattr(r100, "resolve", lambda i: i)  # no download: the files are their traces
    monkeypatch.setattr(r100, "ca_trace", lambda path: traces[path])
    one = {"protein": True, "na": False, "ca_only": False, "accessions": ["P00001"], "assembly1_chains": 1}
    one.update(deposited_residues=130, method="X-RAY DIFFRACTION", chains=1, has=[])
    s11 = next(st for st in r100.STRATA if st["name"] == "S11")

    def drawn(partners):
        s = r100.Sampler(20261017, "2026-10-13")
        s.feats = {i: {**one, "ligand_touching": i != "1AAA"} for i in traces}
        pools = iter([["1AAA"], ["1AAA", *partners]])  # the seed pool, then the accession's entries
        monkeypatch.setattr(s, "pool", lambda nodes, content="experimental": next(pools))
        s.pairs(s11)
        return s

    other = drawn(["2BBB"])  # the only partner of the accession is another protein: no pair
    assert other.entries == [] and [x["id"] for x in other.skipped] == ["1AAA"]
    assert r100.same_protein("1AAA", "2BBB")["one"] is False and r100.same_protein("1AAA", "3CCC")["one"] is True
    both = drawn(["2BBB", "3CCC"])  # the partner that is the same protein is taken, whatever the seeded order
    assert [(e["id"], e["role"], e["pair_type"]) for e in both.entries] == [
        ("1AAA", "seed", "apo/holo"),
        ("3CCC", "partner", "apo/holo"),
    ]


@pytest.fixture
def cath_offline(tmp_path, monkeypatch):
    """foldmap.cath's cache holding the fixtures of tests/data/cath (1LDM, 2LZM); any other entry is a 404."""
    import shutil
    import urllib.error

    from foldmap import cath

    folder = tmp_path / "cath-cache" / "cath"
    folder.mkdir(parents=True)
    for f in (DATA / "cath").iterdir():
        shutil.copy(f, folder / f.name)
    monkeypatch.setattr(cath, "CACHE", tmp_path / "cath-cache")

    def not_found(req, timeout=None):
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, io.BytesIO(b"{}"))

    monkeypatch.setattr(cath.urllib.request, "urlopen", not_found)
    return cath


def test_s5_scores_every_chain_cath_chops(r100, cath_offline, monkeypatch):
    """Plan 15 Task 8: S5 is judged on every sampled entry whose chains CATH chops, not only on the 33 calibration
    chains: each chain's chopping is fetched (foldmap.cath: the PDBe API, cached) and scored as the reference chains
    are (one-domain chains too); a calibration chain keeps its committed chopping. An entry CATH chops no chain of,
    an id that is no PDB entry, or PDBe unreachable: not judged, with the reason, never failed."""
    import urllib.error

    got, note, unscored = r100.domain_reference(DATA / "2LZM.cif", "2LZM")  # not in the reference: CATH's one domain
    assert note == "" and [(g["chain"], g["source"], g["cath"], g["n_ref"]) for g in got] == [
        ("A", "CATH", ["2lzmA00"], 1)
    ]
    got, _, _ = r100.domain_reference(FID / "1LDM.cif.gz", "1LDM")  # in the reference: scored on it, once
    assert [(g["chain"], g["source"], g["n_ref"]) for g in got] == [("A", "reference", 2)]
    assert r100.cath_chains("1LDM") == [
        {"pdb": "1LDM", "chain": "A", "domains": [[[1, 162]], [[163, 329]]], "cath": ["1ldmA01", "1ldmA02"]}
    ]  # shaped as the reference's own chain
    assert r100.domain_reference(DATA / "1UBQ.cif", "1UBQ") == (None, "CATH chops none of its chains", [])
    assert r100.domain_reference(DATA / "1UBQ.cif", "AF_AFP0CG48F1") == (
        None,
        "not a PDB entry (CATH chops PDB entries)",
        [],
    )

    def down(req, timeout=None):
        raise urllib.error.URLError("no route")

    monkeypatch.setattr(cath_offline.urllib.request, "urlopen", down)
    none, why, _ = r100.domain_reference(DATA / "1UBQ.cif", "9XYZ")
    assert none is None and why.startswith("CATH chopping not fetched (ValueError: could not reach www.ebi.ac.uk")
    e = {"id": "9XYZ", "facts": {"protein": True}, "views": {"default": {"result": {"science": {}}}}}
    s5 = r100.science({**e, "domains_cath": None, "domains_cath_note": why})["S5"]
    assert s5["pass"] is None and s5["status"] == "not judged" and s5["detail"] == why


def test_report_counts_the_s5_chains_and_credits_cath(r100, tmp_path):
    """Plan 15 Task 8: report.md says how many chains S5 scored and how many of them CATH divides into two domains
    or more, with CATH's attribution (CC BY 4.0), each chain's source, and why the other protein entries were not
    judged."""
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}
    view = {"default": {"result": {"science": {}}}}
    two = {"chain": "A", "overlap": 0.9, "f1": 0.8, "recall": 1.0, "precision": 0.67, "n_ref": 2, "n_auto": 3}
    one = {**two, "chain": "B", "f1": 1.0, "precision": 1.0, "n_ref": 1, "n_auto": 1, "source": "CATH"}
    entries = [
        {
            "id": "1ABC",
            "name": "1ABC",
            "facts": {"protein": True},
            "views": view,
            "criteria": crit,
            "domains_cath": [two, one],
        },
        {
            "id": "2DEF",
            "name": "2DEF",
            "facts": {"protein": True},
            "views": view,
            "criteria": crit,
            "domains_cath": [{**two}],
        },
        {
            "id": "3GHI",
            "name": "3GHI",
            "facts": {"protein": True},
            "views": view,
            "criteria": crit,
            "domains_cath": None,
        },
    ]
    entries[1]["domains_cath"][0]["source"] = "CATH"
    entries[2]["domains_cath_note"] = "CATH chops none of its chains"
    r100.report([r100._verdict(e) for e in entries], tmp_path, None, total=3)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "3 chains scored in 2 entries; 2 of them with at least two CATH domains." in md
    assert "CC BY 4.0" in md and "https://www.cathdb.info" in md and "Sillitoe" in md
    assert "| 1ABC | B | CATH | 1 | 1 | 0.90 | 1.00 | 1.00 | 1.00 |" in md
    assert "| mean of 2 (two domains or more) |" in md and "| mean of 1 (one domain) |" in md
    assert "- CATH chops none of its chains: 1" in md


def test_a_chain_whose_s5_scoring_crashes_is_counted_in_the_report(r100, cath_offline, monkeypatch, tmp_path):
    """Task 8 fix round 1: a chain CATH chops whose scoring raises (chain_inputs or auto_domains crashing) is never
    silently left out when another chain of the entry is scored: domain_reference lists it with the exception's type,
    and report.md counts the unscored chains by reason (grouped by type, not message) and names them."""
    from foldmap import domain_score

    real = domain_score.chain_score

    def crashing(path, c):
        if c["chain"] == "B":
            raise RuntimeError(f"auto_domains broke on {path}")
        return real(path, c)

    monkeypatch.setattr(domain_score, "chain_score", crashing)
    chopped = r100.cath_chains("2LZM")[0]
    monkeypatch.setattr(
        r100, "cath_chains", lambda pdb: [chopped, {**chopped, "chain": "B"}, {**chopped, "chain": "C"}]
    )
    got, note, unscored = r100.domain_reference(DATA / "2LZM.cif", "2LZM")
    assert [g["chain"] for g in got] == ["A"]  # C has no residue in the file
    assert [(u["chain"], u["reason"]) for u in unscored] == [("B", "RuntimeError"), ("C", r100.NO_RESIDUE)]
    assert "chain B not scored (RuntimeError: auto_domains broke" in note
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}
    view = {"default": {"result": {"science": {}}}}

    def entry(i, scores, unscored):
        e = {"id": i, "name": i, "facts": {"protein": True}, "views": view, "criteria": crit, "domains_cath": scores}
        return r100._verdict({**e, "domains_cath_unscored": unscored, "domains_cath_note": "x (y)"})

    other = [{"chain": "D", "reason": "RuntimeError", "detail": "another message"}]
    r100.report([entry("2LZM", got, unscored), entry("1ABC", None, other)], tmp_path, None, total=2)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "3 chains CATH chops not scored, by reason:" in md
    assert "- RuntimeError: 2 (2LZM B, 1ABC D)" in md and f"- {r100.NO_RESIDUE}: 1 (2LZM C)" in md
    assert "- no chain CATH chops could be scored: 1" in md


# ------------------------------------------------------------------------------ plan 28 Task 1, fix round 1


def _stamped_entry(r100, i, stamp):
    crit = {k: {"pass": True, "detail": ""} for k in r100.CRITERIA}
    e = {"id": i, "name": i, "stratum": "S1", "bin": "S", "views": {}, "facts": {"chains": 1}, "criteria": crit}
    e["source_stamp"] = stamp
    return r100._verdict(e)


def test_a_report_of_entries_measured_by_different_code_says_so(r100, tmp_path):
    """Ruling Q13: entries measured at X and rerun at Y make a mixed report, never quoted as one measurement."""
    other = r100.SOURCE_STAMP | {"code_fingerprint": "0" * 16}
    r100.report(
        [_stamped_entry(r100, "1ABC", r100.SOURCE_STAMP), _stamped_entry(r100, "2DEF", other)], tmp_path, None, total=2
    )
    doc = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert doc["mixed_code"] is True and sorted(doc["entry_code"].values()) == [1, 1]
    assert "MIXED CODE" in (tmp_path / "report.md").read_text(encoding="utf-8")
    r100.report([_stamped_entry(r100, "1ABC", r100.SOURCE_STAMP)], tmp_path, None, total=1)
    assert json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))["mixed_code"] is False


def _paper_sample(tmp_path, fingerprint):
    sample = tmp_path / "sample.json"
    protocol = {"seed": 1, "source_stamp": {"code_fingerprint": fingerprint}}
    sample.write_text(json.dumps({"protocol": protocol, "entries": [], "replaced": []}), encoding="utf-8")
    return sample


def test_a_paper_sample_runs_once_whole_with_the_code_that_drew_it(r100, tmp_path):
    """Ruling Q14: a second run, --ids, --failed-from or other code is refused without --not-for-paper."""
    import argparse

    out = tmp_path / "run"
    sample = _paper_sample(tmp_path, r100.SOURCE_STAMP["code_fingerprint"])
    plain = argparse.Namespace(ids=None, failed_from=None)
    assert r100.paper_run_problem(out, sample, plain) is None
    assert "whole" in r100.paper_run_problem(out, sample, argparse.Namespace(ids="1ABC", failed_from=None))
    assert "whole" in r100.paper_run_problem(out, sample, argparse.Namespace(ids=None, failed_from="r.json"))
    other = _paper_sample(tmp_path, "0" * 16)
    assert "not the code that drew" in r100.paper_run_problem(out, other, plain)
    out.mkdir()
    (out / "run.log").write_text("1ABC: PASS\n", encoding="utf-8")
    assert "runs once" in r100.paper_run_problem(out, sample, plain)
    with pytest.raises(SystemExit, match="refused"):
        r100.main(["run", "--sample", str(sample), "--out", str(out)])
    assert not (out / "report.json").exists()
    no_protocol = tmp_path / "plain.json"
    no_protocol.write_text(json.dumps({"entries": [], "replaced": []}), encoding="utf-8")
    assert r100.paper_run_problem(out, no_protocol, plain) is None


def test_not_for_paper_marks_the_report_never_to_be_quoted(r100, tmp_path):
    out = tmp_path / "run"
    sample = _paper_sample(tmp_path, "0" * 16)  # other code: refused, unless marked
    assert r100.main(["run", "--sample", str(sample), "--out", str(out), "--not-for-paper"]) == 0
    doc = json.loads((out / "report.json").read_text(encoding="utf-8"))
    assert "not-for-paper" in doc["not_for_paper"] and doc["protocol"]["seed"] == 1
    assert "NOT FOR THE PAPER" in (out / "report.md").read_text(encoding="utf-8")
