"""Automatic domains against reference choppings (plan 12 Task 3, R2): the scoring, the reference set, and the
calibration on the 33 CATH chains (needs the structures in the domain cache: scripts/make_domain_cache.py)."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from foldmap import domain_score as ds

DATA = Path(__file__).parent / "data"
REFERENCE = DATA / "domains" / "cath_reference.json"
CACHE = Path(os.environ.get("FOLDMAP_DOMAIN_CACHE") or Path.home() / ".cache" / "foldmap-domains")


def _cached() -> bool:
    doc = json.loads(REFERENCE.read_text(encoding="utf-8"))
    return all((CACHE / f"{e['pdb']}.cif").is_file() for e in doc["chains"])


needs_cache = pytest.mark.skipif(
    not _cached(), reason=f"domain calibration structures not in {CACHE} (python scripts/make_domain_cache.py)"
)


def test_reference_set_is_the_ruled_one():
    doc = json.loads(REFERENCE.read_text(encoding="utf-8"))
    assert len(doc["chains"]) == 33
    assert all(len(e["domains"]) >= 2 for e in doc["chains"])
    assert "CC BY 4.0" in doc["license"] and "pdbe/api/mappings/cath" in doc["source"]
    lit = {e["pdb"]: e["domains"] for e in doc["literature"]}
    assert lit["2LZM"] == [[[13, 75]], [[1, 12], [76, 164]]]
    assert {e["pdb"] for e in doc["chains"]} >= {"1LDM", "3PGK", "1CZA", "7FAB"}


def test_reference_labels_hold_discontinuous_domains():
    labels = ds.reference_labels(range(1, 21), [[[5, 10]], [[1, 4], [11, 18]]])
    assert labels.tolist() == [1] * 4 + [0] * 6 + [1] * 8 + [-1, -1]


def test_element_labels_fill_loops_from_the_nearest_element():
    assert ds.element_labels(10, [(1, 2, 0), (6, 8, 1)]).tolist() == [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    assert ds.element_labels(4, []).tolist() == [0, 0, 0, 0]


def test_overlap_charges_a_split_domain_what_its_extra_piece_holds():
    ref = np.array([0] * 50 + [1] * 50)
    assert ds.overlap(ref, ref) == 1.0
    assert ds.overlap(ref, np.zeros(100, int)) == 0.5  # one domain called: half the residues are matched
    split = np.array([0] * 25 + [2] * 25 + [1] * 50)
    assert ds.overlap(ref, split) == 0.75
    unjudged = np.array([-1] * 10 + [0] * 40 + [1] * 50)
    assert ds.overlap(unjudged, np.array([1] * 50 + [0] * 50)) == 1.0  # names do not matter


def test_boundaries_within_the_tolerance():
    ref = np.array([0] * 100 + [1] * 100)
    auto = np.array([0] * 115 + [1] * 85)
    assert ds.boundary_scores(ref, auto) == (1.0, 1.0, 1.0)  # 15 residues off: found
    far = np.array([0] * 125 + [1] * 75)
    assert ds.boundary_scores(ref, far) == (0.0, 0.0, 0.0)
    extra = np.array([0] * 100 + [1] * 50 + [2] * 50)
    r, p, f1 = ds.boundary_scores(ref, extra)
    assert (r, p) == (1.0, 0.5) and f1 == pytest.approx(2 / 3)
    assert ds.boundary_scores(np.zeros(50, int), np.zeros(50, int)) == (1.0, 1.0, 1.0)
    assert ds.boundary_scores(ref, np.zeros(200, int))[1] == 0.0  # a split missed: nothing right either


def test_held_out_is_a_fixed_third():
    entries = ds.load_reference(REFERENCE)
    held = ds.held_out(entries)
    assert len(held) == 11 and held == ds.held_out(list(reversed(entries)))


def _scores(entries):
    return [ds.chain_score(CACHE / f"{e['pdb']}.cif", e) for e in entries]


@pytest.mark.xfail(strict=True, reason="auto domains over-split: plan 12 Task 3 step 2")
@needs_cache
def test_auto_domains_meet_the_ruled_targets_on_the_cath_set():
    """Plan 12 ruling R2: mean overlap >= .80, boundary F1 >= .70 and recall >= .85 on the 33 chains, and on the
    held-out third alone (never tuned on) no worse than the targets less .05."""
    entries = ds.load_reference(REFERENCE)
    held = ds.held_out(entries)
    scores = _scores(entries)
    every = ds.summary(scores)
    rest = ds.summary([s for s in scores if (s["pdb"], s["chain"]) in held])
    assert every["overlap"] >= ds.TARGETS["overlap"], every
    assert every["f1"] >= ds.TARGETS["f1"], every
    assert every["recall"] >= ds.TARGETS["recall"], every
    for key in ("overlap", "f1", "recall"):
        assert rest[key] >= ds.TARGETS[key] - 0.05, rest


@pytest.mark.xfail(strict=True, reason="auto domains over-split: plan 12 Task 3 step 2")
def test_t4_lysozyme_finds_the_literature_boundary():
    """2LZM (R2): the literature's domains are N 13-75 and C 1-12 + 76-164; the cut near 75/76 is found and none is
    made near 112/115 (the old split, 37 residues off)."""
    (entry,) = [e for e in json.loads(REFERENCE.read_text(encoding="utf-8"))["literature"] if e["pdb"] == "2LZM"]
    got = ds.chain_score(DATA / "2LZM.cif", entry)
    assert got["n_auto"] == 2 and got["recall"] >= 0.5 and got["precision"] == 1.0, got
    assert got["overlap"] >= 0.85, got
