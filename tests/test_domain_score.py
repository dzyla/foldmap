"""Automatic domains against reference choppings (plan 12 Task 3, R2): the scoring, the reference set, and the
calibration on the 33 CATH chains. The ruled targets are checked in CI on a committed fixture of what the domains are
found from (cath_inputs.npz, scripts/make_domain_fixture.py); the same check on the structure files themselves runs
when they are in the domain cache (scripts/make_domain_cache.py).

What the numbers mean (fix round 1): the parameters were chosen on a seeded two thirds of the 33 chains, and the
other third, drawn chain by chain, held near-duplicates of training chains (6FY9/6FYF, 2YBE/3PGK, the hexokinases).
On the seven held-out chains with no training homologue: overlap .82, recall .80, precision .70, F1 .72. The
held-out third is now drawn by homology group (overlap .80, recall .81, precision .68, F1 .72), but that is after
the fact. The never-seen check is the fresh set (cath_fresh.json): eleven chains of superfamilies none of the 33
has, scored once with the final algorithm and gated only at a sanity floor."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from foldmap import domain_score as ds

DATA = Path(__file__).parent / "data"
REFERENCE = DATA / "domains" / "cath_reference.json"
FIXTURE = DATA / "domains" / "cath_inputs.npz"
FRESH = DATA / "domains" / "cath_fresh.json"
FRESH_FIXTURE = DATA / "domains" / "fresh_inputs.npz"
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


def test_held_out_is_a_fixed_third_of_whole_homology_groups():
    entries = ds.load_reference(REFERENCE)
    held = ds.held_out(entries)
    assert len(held) == 11 and held == ds.held_out(list(reversed(entries)))
    side = {}
    for e in entries:  # no group on both sides: a held-out chain has no homologue among the rest
        assert side.setdefault(e["group"], (e["pdb"], e["chain"]) in held) == ((e["pdb"], e["chain"]) in held), e
    groups = {e["group"] for e in entries if e["pdb"] in ("1CZA", "1HKB", "2YHX", "1BDG")}
    assert len(groups) == 1  # the hexokinases are one group


def _fixture_scores(fixture, entries):
    chains = ds.load_fixture(fixture)
    return [ds.inputs_score(chains[f"{e['pdb']}_{e['chain']}"], e) for e in entries]


def test_auto_domains_meet_the_ruled_targets_on_the_cath_set():
    """Plan 12 ruling R2, in CI: on the 33 chains (their committed inputs) mean overlap >= .80, boundary F1 >= .70
    and recall >= .85. The held-out third (by homology group, reported in the module docstring) only above a sanity
    floor: it was not untouched (see above)."""
    entries = ds.load_reference(REFERENCE)
    held = ds.held_out(entries)
    scores = _fixture_scores(FIXTURE, entries)
    every = ds.summary(scores)
    assert every["overlap"] >= ds.TARGETS["overlap"], every
    assert every["f1"] >= ds.TARGETS["f1"], every
    assert every["recall"] >= ds.TARGETS["recall"], every
    assert ds.summary([s for s in scores if (s["pdb"], s["chain"]) in held])["overlap"] >= 0.70


@pytest.mark.parametrize(
    "key, path",
    [
        ("1LDM_A", "fidelity/1LDM.cif.gz"),
        ("7FAB_H", "r100/7FAB-HL.cif.gz"),
        ("1SLV_B", "r100/1SLV.cif.gz"),
        ("2ILU_A", "r100/2ILU.cif.gz"),
    ],
)
def test_the_fixture_gives_the_domains_the_structure_gives(key, path):
    """The fixture cannot drift from the pipeline: for the reference chains whose structure is committed, its
    inputs give the same domains as reading the file (DSSP, elements, contacts, rows) does now."""
    inputs = ds.load_fixture(FIXTURE)[key]
    assert ds.inputs_domains(inputs) == ds.inputs_domains(ds.chain_inputs(DATA / path, key.split("_")[1]))


@needs_cache
def test_the_structures_give_the_fixtures_domains_and_meet_the_targets():
    """End to end on the structure files (the domain cache): every chain's domains as the fixture's, and the ruled
    targets on the 33."""
    entries = ds.load_reference(REFERENCE)
    chains = ds.load_fixture(FIXTURE)
    scores = []
    for e in entries:
        inputs = ds.chain_inputs(CACHE / f"{e['pdb']}.cif", e["chain"])
        assert ds.inputs_domains(inputs) == ds.inputs_domains(chains[f"{e['pdb']}_{e['chain']}"]), e["pdb"]
        scores.append(ds.inputs_score(inputs, e))
    every = ds.summary(scores)
    assert all(every[k] >= ds.TARGETS[k] for k in ("overlap", "f1", "recall")), every


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
    left in the Rossmann domain); a run at the chain's end stays."""
    from foldmap.features import _no_strays

    assert _no_strays([0, 0, 1, 0, 0], [10, 10, 15, 10, 10]) == [0, 0, 0, 0, 0]
    assert _no_strays([0, 0, 1, 1, 0], [10, 10, 15, 15, 10]) == [0, 0, 1, 1, 0]  # 30 residues: a domain's piece
    assert _no_strays([1, 0, 0, 0], [12, 10, 10, 10]) == [1, 0, 0, 0]
    assert _no_strays([0, 1, 0, 1, 0], [5, 5, 5, 5, 5]) == [0, 0, 0, 0, 0]  # the runs as they stand after each move
