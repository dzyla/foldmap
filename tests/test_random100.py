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
    assert r100.stem_gap_limit(d.figure) == pytest.approx(figcheck.stem_gap_tolerance(_INK_MARGIN * em))
    assert 0.3 < r100.stem_gap_limit(d.figure) < 1.0


def test_tether_lengths_are_measured_from_mark_to_target(r100):
    """A20 reads every ligand tether (they are patches, drawn in pieces where no body or text is; reading lines
    found none): myoglobin's heme and hydroxide, six tethers, the longest 6.06 page units."""
    from foldmap.spec import draw_spec, figure_spec
    from foldmap.style import Style

    d = draw_spec(figure_spec(DATA / "fidelity" / "1MBN.cif.gz"))
    got = r100.tether_lengths(d, Style())
    assert len(got) == 6 and all(g.startswith("ligand-tether:") for g, *_ in got)
    assert max(n for _, n, _, _ in got) == pytest.approx(6.06, abs=0.01)


def test_a20_gates_only_tethers_to_other_chains(r100):
    """C2 (random-100): A20 bounds only a tether to another chain than the ligand's home (ruling); a long tether
    to the home chain is reported, not failed. Each tether records its home and target chains (myoglobin: one
    chain, so every tether is a home one)."""
    from foldmap.spec import draw_spec, figure_spec
    from foldmap.style import Style

    d = draw_spec(figure_spec(DATA / "fidelity" / "1MBN.cif.gz"))
    got = r100.tether_lengths(d, Style())
    assert got and all(home == target == "A" for _, _, home, target in got)
    view = {"exit": 0, "timed_out": False, "traceback": False, "python_warnings": [], "written": {}, "seconds": 1.0}

    def a20(tethers):
        return r100.criteria({"views": {"default": {**view, "result": {"tethers": tethers}}}})["A20"]

    home = a20([["ligand-tether:URC:A300:0", 9.01, "A", "A"], ["ligand-tether:URC:A300:1", 2.0, "A", "A"]])
    assert home["pass"] is True and "reported: ligand-tether:URC:A300:0 9.01 (home A)" in home["detail"]
    other = a20([["ligand-tether:CLA:B609:2", 7.08, "B", "H"], ["ligand-tether:CLA:B609:0", 3.0, "B", "B"]])
    assert other["pass"] is False and other["detail"].startswith("ligand-tether:CLA:B609:2 7.08 (B to H)")
    assert a20([["ligand-tether:CLA:B609:2", 6.4, "B", "H"]])["pass"] is True


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
    assert got["pass"] is True  # after the extra helix B:10-11 names are not compared


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
    assert r100.stem_gap_limit(d.figure) == pytest.approx(figcheck.stem_gap_tolerance(_INK_MARGIN * em))


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
    that stays in one leaflet fails. Bacteriorhodopsin's five helices whose CA runs reach both faces pass as drawn."""
    from dataclasses import replace

    import numpy as np

    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "1C3W.cif"
    d = draw_spec(figure_spec(path, membrane="auto", assembly="asu"))
    m = r100.membrane_checks(d, path)
    assert m["band"] and m["tm_helices"] == 5 and m["tm_short"] == [], m
    lay = d.layout
    (y0, y1), mid = lay.membrane["y"], sum(lay.membrane["y"]) / 2
    hid = next(k for k, p in lay.placed.items() if abs(np.sin(p.angle)) > 0.9 and p.length > y1 - y0)
    p = lay.placed[hid]
    half = abs(p.c_port[1] - p.n_port[1]) / 2
    lay.placed[hid] = replace(p, cy=y0 + 0.7 * (y1 - y0) - half)  # crosses the lower face, ends at 0.70: spans
    assert r100.membrane_checks(d, path)["tm_short"] == []
    lay.placed[hid] = replace(p, cy=mid - 0.05 - half)  # crosses the lower face, ends short of the middle
    assert [x.split()[0] for x in r100.membrane_checks(d, path)["tm_short"]] == [hid]
    lay.placed[hid] = replace(p, cy=mid + 0.05 + half)  # the upper leaflet only
    assert [x.split()[0] for x in r100.membrane_checks(d, path)["tm_short"]] == [hid]
    short = (y1 - y0) / 3  # a short helix centred in the band, entering it through neither face (fix round 1)
    lay.placed[hid] = replace(p, cy=mid, length=short * abs(p.length) / max(2 * half, 1e-9))
    assert [x.split()[0] for x in r100.membrane_checks(d, path)["tm_short"]] == [hid]


def test_a18_a_run_into_a_port_just_inside_the_band_is_at_that_port(r100):
    """C26 (random-100): a loop runs sideways inside the band only next to one of its own ports inside the band.
    A port 0.2 inside a face (7SMR E:162-163, a helix end at 3.46 in the band 3.17..10.24) was not "inside" by the
    0.3 margin the run itself is measured with, so the run that must reach it from above counted as in the slab.
    The port is measured against the band itself; a run between two ports outside the band still counts."""
    from foldmap.route import Loop
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "1C3W.cif"
    d = draw_spec(figure_spec(path, membrane="auto", assembly="asu"))
    y0, _ = d.layout.membrane["y"]
    run = [(-30.0, y0 - 1.0), (-30.0, y0 + 0.6), (-27.0, y0 + 0.6)]
    d.loops.append(Loop("p:1-2", "p:3-4", [*run, (-27.0, y0 + 0.2)], False, False))  # into a port 0.2 inside
    d.loops.append(Loop("q:1-2", "q:3-4", [*run, (-27.0, y0 - 1.0)], False, False))  # back out below the face
    assert r100.membrane_checks(d, path)["loops_in_slab"] == ["q:1-2>q:3-4"]


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
    assert e["pass"] and r100.science(e)["S1"] == {"pass": False, "detail": found["S1"][0]}
    assert r100.science(e)["S4"]["pass"] is None
    r100.report([e], tmp_path, None, total=1)
    md = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "1 of 1 pass." in md and "## Science checks (S1-S6), apart from A1-A20" in md
    assert "| S1 sheet neighbours pair | 0 | 1 | 0 |" in md and "S1: sheet 0: a / b drawn side by side" in md


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

    got = r100.domain_reference(FID / "1LDM.cif.gz", "1LDM")
    assert [g["chain"] for g in got] == ["A"] and got[0]["overlap"] > 0.9 and got[0]["n_ref"] == 2
    assert r100.domain_reference(DATA / "1UBQ.cif", "1UBQ") is None
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
    assert "In-sample: these chains are among the 33" in md
    assert "| 1ABC | A | 2 | 4 | 0.60 | 0.50 | 1.00 | 0.33 |" in md
