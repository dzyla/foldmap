"""Base pairing against an independent annotator (RNApolis, stored in tests/data/na/reference.json)."""

import json
from pathlib import Path

import pytest

from foldmap.io import _nucleic_strands, load_backbone, load_nucleic, read_model
from foldmap.nucleic import BasePair, pair_nucleotides, parent_base

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB"]
REFERENCE = json.loads((NA / "reference.json").read_text(encoding="utf-8"))


def test_every_fixture_has_reference_pairs():
    assert sorted(REFERENCE) == sorted(FIXTURES)
    for pdb in FIXTURES:
        assert (NA / f"{pdb}.cif.gz").is_file()
        assert REFERENCE[pdb], pdb


@pytest.mark.parametrize("pdb", FIXTURES)
def test_every_reference_nucleotide_is_in_its_fixture(pdb):
    model, _ = read_model(NA / f"{pdb}.cif.gz", "asu")
    present = {(c.name, r.seqid.num, r.seqid.icode.strip()) for c in model for r in c}
    cited = {tuple(nt) for a, b, _ in REFERENCE[pdb] for nt in (a, b)}
    assert cited <= present, (pdb, sorted(cited - present))


def _found(pdb):
    model, _ = read_model(NA / f"{pdb}.cif.gz", "asu")
    return _nucleic_strands(model)


def _named(strands, pairs):
    def key(s, i):
        lab = strands[s].labels[i]
        return (lab.chain, lab.seq, lab.icode)

    return {frozenset((key(p.a, p.i), key(p.b, p.j))): p.kind for p in pairs}


@pytest.mark.parametrize("pdb", FIXTURES)
def test_pairs_agree_with_rnapolis(pdb):
    found = _found(pdb)
    strands = [s for s, _ in found]
    ours = _named(strands, pair_nucleotides([run for _, run in found]))
    theirs = {frozenset((tuple(a), tuple(b))): kind for a, b, kind in REFERENCE[pdb]}
    both = set(ours) & set(theirs)
    precision, recall = len(both) / len(ours), len(both) / len(theirs)
    assert precision >= 0.95, (pdb, sorted(map(sorted, set(ours) - set(theirs))))
    assert recall >= 0.90, (pdb, sorted(map(sorted, set(theirs) - set(ours))))
    assert 2 * precision * recall / (precision + recall) >= 0.95
    assert all(ours[k] == theirs[k] for k in both)  # a wobble is a wobble in both


def test_trna_wobble_and_modified_nucleotides_pair():
    found = _found("1EHZ")
    strands = [s for s, _ in found]
    assert len(strands) == 1 and len(strands[0]) == 76  # 5MC and YG, unknown to gemmi, keep the strand whole
    pairs = _named(strands, pair_nucleotides([run for _, run in found]))
    key = lambda a, b: frozenset((("A", a, ""), ("A", b, "")))  # noqa: E731
    assert pairs[key(4, 69)] == "GU"
    assert pairs[key(10, 25)] == "WC"  # 2MG10-C25
    assert pairs[key(30, 40)] == "WC"  # G30-5MC40
    assert pairs[key(49, 65)] == "WC"  # 5MC49-G65


def test_pair_closing_a_two_nucleotide_loop_is_rejected():
    found = _found("30QM")
    strands = [s for s, _ in found]
    pairs = _named(strands, pair_nucleotides([run for _, run in found]))
    assert frozenset((("1", 82, ""), ("1", 85, ""))) not in pairs
    assert frozenset((("1", 81, ""), ("1", 86, ""))) in pairs


def test_parent_bases_of_modified_nucleotides():
    model, _ = read_model(NA / "1EHZ.cif.gz", "asu")
    names = {r.name: parent_base(r) for r in model["A"]}
    assert names["PSU"] == "U" and names["2MG"] == "G" and names["5MC"] == "C" and names["1MA"] == "A"
    assert names["H2U"] == "U" and names["OMC"] == "C"


def test_pairs_are_mutual_and_each_nucleotide_pairs_once():
    found = _found("4YBB")
    pairs = pair_nucleotides([run for _, run in found])
    used = [(p.a, p.i) for p in pairs] + [(p.b, p.j) for p in pairs]
    assert len(used) == len(set(used))
    assert all((p.a, p.i) < (p.b, p.j) for p in pairs)
    assert all(isinstance(p, BasePair) and p.kind in ("WC", "GU") for p in pairs)


def test_no_nucleotides_no_pairs():
    assert pair_nucleotides([]) == []


def test_load_nucleic_keeps_four_tuple_pairs_and_the_kinds():
    path = NA / "30QM.cif.gz"  # protein and G·U wobbles (1EHZ, RNA only, has no Backbone to load against)
    na = load_nucleic(path, load_backbone(path, "asu"), "asu")
    assert na.pairs == [tuple(p[:4]) for p in na.base_pairs]
    assert {p.kind for p in na.base_pairs} == {"WC", "GU"}
    assert all(len(p) == 4 for p in na.pairs)
