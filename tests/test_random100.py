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
