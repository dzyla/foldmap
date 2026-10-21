"""The homologue tuning set (docs/validation/homologues/tuning.json) and its campaign script (scripts/homologues.py):
the set's form only (no network, no US-align), and the report's table from a synthetic result."""

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SET = ROOT / "docs" / "validation" / "homologues" / "tuning.json"
SCRIPT = ROOT / "scripts" / "homologues.py"


@pytest.fixture(scope="module")
def hom():
    spec = importlib.util.spec_from_file_location("homologues", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # its dataclass reads its own module
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def pairs(hom):
    return hom.read_set(SET)


def test_the_set_parses_and_every_pair_names_two_chains(pairs):
    assert len(pairs) >= 12
    for p in pairs:
        assert set(p) == {"a", "b", "kind", "note"}
        for side in ("a", "b"):
            assert set(p[side]) == {"id", "chain"}
            assert re.fullmatch(r"[0-9][0-9A-Z]{3}", p[side]["id"]), p
            assert re.fullmatch(r"[0-9A-Za-z]{1,4}", p[side]["chain"]), p
        assert p["note"].strip()
        assert (p["a"]["id"], p["a"]["chain"]) != (p["b"]["id"], p["b"]["chain"])


def test_positives_and_negatives_are_both_present(pairs):
    kinds = [p["kind"] for p in pairs]
    assert set(kinds) == {"positive", "negative"}
    assert kinds.count("negative") >= 4


def test_pairs_are_sorted_and_unique(pairs, hom):
    keys = [(p["a"]["id"], p["a"]["chain"], p["b"]["id"], p["b"]["chain"]) for p in pairs]
    assert keys == sorted(keys)
    assert len({hom.stem(p) for p in pairs}) == len(pairs)


def test_f7_pair_is_in_the_set(pairs):
    assert any(p["a"] == {"id": "2VWD", "chain": "A"} and p["b"] == {"id": "2ZB5", "chain": "A"} for p in pairs)


def test_report_lists_every_failure_by_check(hom, tmp_path):
    pair = {"a": {"id": "1ABC", "chain": "A"}, "b": {"id": "2DEF", "chain": "B"}, "kind": "positive", "note": "x"}
    ok = {"pass": True, "detail": "fine"}
    result = {
        "pair": hom.stem(pair),
        "kind": "positive",
        "note": "",
        "alignment": {"method": "structure", "tm_score": 0.61},
        "mates": {"mates": 9, "pieces": 1, "extras": 2, "only_a": 0, "same_fold": True},
        "like": {
            "H1": ok,
            "H2": {**ok, "value": {"median": 0.2, "p90": 0.9}, "baseline": {"median": 3.0, "p90": 5.0}},
            "H3": {"pass": False, "detail": "0.9000 of 20 pairs agree", "value": {"share": 0.9, "pairs": 20}},
            "H4": ok,
            "H5": {"pass": False, "detail": "lost A5", "value": {"alone": ["A5"], "lost": ["A5"]}},
        },
    }
    (tmp_path / f"{hom.stem(pair)}.json").write_text(json.dumps(result))
    md = hom.write_report([pair, {**pair, "b": {"id": "3GHI", "chain": "A"}}], tmp_path, tmp_path / "r.md")
    text = md.read_text()
    row = "| 1ABC_A-2DEF_B | positive | structure | 0.610 | 9/1/2/0 | pass | pass | **fail** | pass | **fail** |"
    assert row in text
    assert "0.20 / 0.90 | 3.00 / 5.00 | 0.900 (20) | – | A5 |" in text
    assert "- **H3**: \n  - 1ABC_A-2DEF_B (positive): 0.9000 of 20 pairs agree" in text
    assert "No result (the pair's child failed; see its log): 1ABC_A-3GHI_A" in text


def _domain_list() -> str:
    """A small CATH domain list: eligible whole-chain domains in two architectures, and one of each kind the sample
    leaves out (a chain of two domains, too short, too long, an obsolete entry, chain 0, a tuning-set entry)."""
    rows = []
    n = 0

    def add(cath, s35, s95, length=120, resolution=2.0, chain="A", number="00", pdb=None):
        nonlocal n
        n += 1
        name = (pdb or f"{1 + n % 9}x{n:02d}").lower() + chain + number
        rows.append(f"{name} {' '.join(map(str, cath))} {s35} 1 {s95} 1 1 {length} {resolution:.3f}")

    for h, s35 in ((720, 1), (720, 2), (300, 1), (300, 2), (300, 3)):
        for s95 in (1, 2, 3):
            add((3, 40, 50, h), s35, s95)
    for s35 in range(1, 7):
        add((3, 40, 30, 10), s35, 1)
    for s35 in (1, 2, 3, 4):
        add((2, 60, 40, 10), s35, s35)
    for s35 in range(1, 6):
        add((2, 60, 120, 10), s35, 1)
    add((3, 40, 50, 720), 1, 1, number="01")
    add((3, 40, 50, 720), 1, 1, length=59)
    add((3, 40, 50, 720), 1, 1, length=601)
    add((3, 40, 50, 720), 1, 1, resolution=1000.0)
    add((3, 40, 50, 720), 1, 1, chain="0")
    add((3, 40, 50, 720), 1, 1, pdb="1UBQ")
    return "# CathDomainList\n" + "\n".join(rows) + "\n"


def test_the_sample_is_seeded_and_keeps_each_stratum(hom, monkeypatch):
    """sample: the same seed gives the same set; remote pairs share a topology, close pairs a superfamily and an S35
    cluster in two S95 clusters, negatives an architecture in two topologies; only eligible chains, each entry once,
    the tuning set's left out; a pair the check refuses is recorded with its reason and another drawn."""
    strata = (("remote", 3, "positive"), ("close", 2, "positive"), ("negative", 2, "negative"))
    monkeypatch.setattr(hom, "STRATA", strata)
    domains = hom.read_domain_list(_domain_list())
    assert len(domains) == 36 and len(hom.eligible(domains, {"1UBQ"})) == 30
    calls = []

    def check(stratum, a, b):
        calls.append((stratum, a.name, b.name))
        got = {"residues_a": 100, "residues_b": 100, "identity": 0.5, "tm_score": 0.6, "aligned": 90}
        return ("refused" if len(calls) == 1 else ""), got

    data = hom.sample(domains, {"3.40.50": "Rossmann fold"}, 7, {"1UBQ"}, check=check)
    calls.clear()
    assert hom.sample(domains, {"3.40.50": "Rossmann fold"}, 7, {"1UBQ"}, check=check) == data
    assert data["seed"] == 7 and data["cath_release"] == hom.CATH_RELEASE
    assert len(data["rejected"]) >= 1 and data["rejected"][0]["reason"] == "refused"
    by = {}
    for p in data["pairs"]:
        by.setdefault(p["stratum"], []).append(p)
        assert p["kind"] == ("negative" if p["stratum"] == "negative" else "positive")
        assert p["note"].startswith(p["stratum"] + ":")
    assert {k: len(v) for k, v in by.items()} == {"remote": 3, "close": 2, "negative": 2}
    names = [(p["cath"]["a"].split(), p["cath"]["b"].split()) for p in data["pairs"]]
    for p, ((na, ca), (nb, cb)) in zip(data["pairs"], names):
        assert na.endswith("00") and nb.endswith("00") and na[4] != "0" and nb[4] != "0"
        top_a, top_b = ca.rsplit(".", 1)[0], cb.rsplit(".", 1)[0]
        if p["stratum"] == "remote":
            assert top_a == top_b
        elif p["stratum"] == "close":
            assert ca == cb
        else:
            assert top_a != top_b and top_a.rsplit(".", 1)[0] == top_b.rsplit(".", 1)[0]
    ids = [p[s]["id"] for p in data["pairs"] for s in ("a", "b")]
    assert len(ids) == len(set(ids)) and "1UBQ" not in ids
    keys = [(p["a"]["id"], p["a"]["chain"], p["b"]["id"], p["b"]["chain"]) for p in data["pairs"]]
    assert keys == sorted(keys)
    assert hom.sample(domains, {}, 8, {"1UBQ"}, check=check)["pairs"] != data["pairs"]  # another seed, another set


def test_the_held_out_summary_counts_by_stratum(hom, tmp_path):
    """A set whose pairs carry a stratum: the table gains the column, the report the summary and summary.json."""
    ok = {"pass": True, "detail": "fine"}

    def result(pair, h2, base, h3_pass=True):
        return {
            "pair": hom.stem(pair),
            "kind": pair["kind"],
            "note": "",
            "alignment": {"method": "structure", "tm_score": 0.6},
            "mates": {"mates": 9, "pieces": 0, "extras": 0, "only_a": 0, "same_fold": True},
            "like": {
                "H1": ok,
                "H2": {**ok, "value": {"median": h2, "p90": 2 * h2}, "baseline": {"median": base, "p90": 2 * base}},
                "H3": {
                    "pass": h3_pass,
                    "detail": "x",
                    "value": {"share": 1.0, "pairs": 10},
                    "baseline": {"share": 0.5, "pairs": 10},
                },  # fmt: skip
                "H4": ok,
                "H5": ok,
            },
        }

    def pair(a, kind, stratum):
        return {"a": {"id": a, "chain": "A"}, "b": {"id": "9ZZZ", "chain": "B"}, "kind": kind, "stratum": stratum}

    pairs = [
        pair("1AAA", "positive", "remote"),
        pair("2BBB", "positive", "close"),
        pair("3CCC", "negative", "negative"),
    ]
    (tmp_path / f"{hom.stem(pairs[0])}.json").write_text(json.dumps(result(pairs[0], 0.2, 3.0)))
    (tmp_path / f"{hom.stem(pairs[1])}.json").write_text(json.dumps(result(pairs[1], 0.4, 1.0, h3_pass=False)))
    neg = result(pairs[2], 0.0, 2.0)
    neg["like"]["H2"] = neg["like"]["H3"] = neg["like"]["H4"] = {"pass": None, "detail": "not judged"}
    (tmp_path / f"{hom.stem(pairs[2])}.json").write_text(json.dumps(neg))
    text = hom.write_report(pairs, tmp_path, tmp_path / "r.md").read_text()
    assert "| 1AAA_A-9ZZZ_B | positive | remote | structure | 0.600 |" in text
    assert "| close | 1 | 1 | 0 | 1 / 0 / 0 | 1 / 0 / 0 | 0 / 1 / 0 | 1 / 0 / 0 | 1 / 0 / 0 |" in text
    assert "Positives passing H1-H5: 1 of 2; negatives refused (H1): 1 of 1." in text
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["positives"]["h2_median_independent"] == 2.0 and summary["positives"]["h2_median_alike"] == 0.3
    assert summary["positives"]["h2_closer_alike"] == 2 and summary["strata"]["negative"]["checks"]["H2"] == {
        "pass": 0, "fail": 0, "not judged": 1}  # fmt: skip


def test_the_held_out_set_is_as_drawn(hom, pairs):
    """heldout.json (sample --seed 20261021): 20 remote, 10 close, 10 negative pairs, sorted, each entry once, none
    of the tuning set's entries."""
    data = json.loads((SET.parent / "heldout.json").read_text(encoding="utf-8"))
    held = data["pairs"]
    assert data["seed"] == 20261021 and data["cath_release"] == hom.CATH_RELEASE
    assert [sum(p["stratum"] == s for p in held) for s, _, _ in hom.STRATA] == [n for _, n, _ in hom.STRATA]
    keys = [(p["a"]["id"], p["a"]["chain"], p["b"]["id"], p["b"]["chain"]) for p in held]
    assert keys == sorted(keys)
    ids = [p[s]["id"] for p in held for s in ("a", "b")]
    assert len(ids) == len(set(ids)) == 80
    assert not set(ids) & {p[s]["id"] for p in pairs for s in ("a", "b")}
