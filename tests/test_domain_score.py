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


def test_t4_lysozyme_finds_the_literature_boundary():
    """2LZM (R2): the literature's domains are N 13-75 and C 1-12 + 76-164; the cut near 75/76 is found and none is
    made near 112/115 (the old split, 37 residues off)."""
    (entry,) = [e for e in json.loads(REFERENCE.read_text(encoding="utf-8"))["literature"] if e["pdb"] == "2LZM"]
    got = ds.chain_score(DATA / "2LZM.cif", entry)
    assert got["n_auto"] == 2 and got["recall"] >= 0.5 and got["precision"] == 1.0, got
    assert got["overlap"] >= 0.85, got


def _blob(n: int, centre) -> np.ndarray:
    """n points 3.8 Å apart on a snake through a cube about `centre`: a compact domain's CA atoms."""
    side = int(np.ceil(n ** (1 / 3)))
    pts = [
        (x if (y + z * side) % 2 == 0 else side - 1 - x, y if z % 2 == 0 else side - 1 - y, z)
        for z in range(side)
        for y in range(side)
        for x in range(side)
    ]
    return np.array(pts[:n], float) * 3.8 + np.asarray(centre, float)


def test_residue_domains_find_a_discontinuous_domain():
    """T4 lysozyme's pattern: the chain starts in one domain (30 residues), makes the other (60), and comes back
    (60): two domains, the first holding both ends."""
    from foldmap.features import residue_domains

    far = _blob(90, (0, 0, 0))
    near = _blob(60, (60, 0, 0))
    ca = np.vstack([far[:30], near, far[30:]])
    assert residue_domains(ca).tolist() == [0] * 30 + [1] * 60 + [0] * 60


def test_residue_domains_find_domains_that_alternate_twice():
    """5JDN's pattern: A, B, A, B along the chain, each run 50 residues: one cut of two runs."""
    from foldmap.features import residue_domains

    a, b = _blob(100, (0, 0, 0)), _blob(100, (60, 0, 0))
    ca = np.vstack([a[:50], b[:50], a[50:], b[50:]])
    assert residue_domains(ca).tolist() == ([0] * 50 + [1] * 50) * 2


def test_residue_domains_keep_a_group_whole_and_one_domain_alone():
    from foldmap.features import residue_domains

    ca = np.vstack([_blob(80, (0, 0, 0)), _blob(80, (60, 0, 0))])
    assert residue_domains(ca).tolist() == [0] * 80 + [1] * 80
    assert residue_domains(ca, keep=[[10, 150]]).tolist() == [0] * 160  # a βα row spanning both is never cut
    assert residue_domains(_blob(160, (0, 0, 0))).tolist() == [0] * 160


def test_a_stray_element_goes_home():
    """An element run of under DOMAIN_STRAY residues between two runs of another domain joins it (1LDM's 247-261
    left in the Rossmann domain); a run at the chain's end stays (2LZM's helix A)."""
    from foldmap.features import _no_strays

    assert _no_strays([0, 0, 1, 0, 0], [10, 10, 15, 10, 10]) == [0, 0, 0, 0, 0]
    assert _no_strays([0, 0, 1, 1, 0], [10, 10, 15, 15, 10]) == [0, 0, 1, 1, 0]  # 30 residues: a domain's piece
    assert _no_strays([1, 0, 0, 0], [12, 10, 10, 10]) == [1, 0, 0, 0]
