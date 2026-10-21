"""The homologue tuning set (docs/validation/homologues/tuning.json) and its campaign script (scripts/homologues.py):
the set's form only (no network, no US-align), and the report's table from a synthetic result."""

import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SET = ROOT / "docs" / "validation" / "homologues" / "tuning.json"
SCRIPT = ROOT / "scripts" / "homologues.py"


@pytest.fixture(scope="module")
def hom():
    spec = importlib.util.spec_from_file_location("homologues", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
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
