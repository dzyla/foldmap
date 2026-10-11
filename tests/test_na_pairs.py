"""Base pairing against an independent annotator (RNApolis, stored in tests/data/na/reference.json)."""

import json
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap import nucleic
from foldmap.io import _nucleic_strands, load_backbone, load_nucleic, read_model
from foldmap.nucleic import BasePair, pair_nucleotides, parent_base

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB", "1RVA", "4LD0", "1HYS", "1AOI"]
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
    assert {"5MC", "YYG"} <= {lab.name for lab in strands[0].labels}
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


def test_pseudouridine_pairs_with_its_turned_face():
    # Ψ hangs from C5, so its Watson-Crick carbonyl is O2: A2590 N6-Ψ2604 O2 3.0 Å (23S rRNA, 4YBB)
    found = _found("4YBB-psi")
    strands = [s for s, _ in found]
    pairs = _named(strands, pair_nucleotides([run for _, run in found]))
    assert strands[0].labels[2604 - 2586].name == "PSU"
    assert pairs[frozenset((("DA", 2590, ""), ("DA", 2604, "")))] == "WC"


def test_an_untabulated_purine_pairs_as_its_parent_base():
    found = _found("4YBB")
    (run,) = [run for _, run in found]
    g7m = next(r for r in run if r.seqid.num == 527)
    info = gemmi.find_tabulated_residue(g7m.name)
    assert g7m.name == "G7M" and not (info and info.is_nucleic_acid())  # gemmi does not know it as a nucleotide
    assert parent_base(g7m) == "G"  # from its atoms: N9 and O6
    pairs = _named([s for s, _ in found], pair_nucleotides([run]))
    assert pairs[frozenset((("AA", 522, ""), ("AA", 527, "")))] == "WC"


def _tRNA_pairs():
    found = _found("1EHZ")
    run = found[0][1]
    return run, {(run[p.i].seqid.num, run[p.j].seqid.num) for p in pair_nucleotides([run])}


def test_a_trans_pair_is_not_a_helix_pair(monkeypatch):
    # The Levitt pair G15·C48 bonds like G-C but with its glycosidic bonds on opposite sides
    assert (15, 48) not in _tRNA_pairs()[1]
    monkeypatch.setattr(nucleic, "_cis", lambda *a: True)
    assert (15, 48) in _tRNA_pairs()[1]


def _tilt(res, degrees):
    """Turn the base of res (and its C1') about the N3-N4 line: the N1-N3 and O6-N4 bonds of a G-C pair hold."""
    p, q = (np.array(res.find_atom(n, "*").pos.tolist()) for n in ("N3", "N4"))
    k = (q - p) / np.linalg.norm(q - p)
    t = np.radians(degrees)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
    for a in res:
        if ("'" not in a.name and "P" not in a.name) or a.name == "C1'":
            a.pos = gemmi.Position(*(R @ (np.array(a.pos.tolist()) - p) + p))


@pytest.mark.parametrize("degrees, plane_max, paired", [(0, 45.0, True), (60, 45.0, False), (60, 70.0, True)])
def test_bases_far_from_coplanar_do_not_pair(monkeypatch, degrees, plane_max, paired):
    run = _found("1EHZ")[0][1]
    c72 = next(r for r in run if r.seqid.num == 72)
    _tilt(c72, degrees)
    monkeypatch.setattr(nucleic, "PLANE_MAX", plane_max)
    assert ((1, 72) in {(run[p.i].seqid.num, run[p.j].seqid.num) for p in pair_nucleotides([run])}) is paired


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
    path = NA / "1EHZ.cif.gz"
    na = load_nucleic(path, load_backbone(path, "asu"), "asu")
    assert na.pairs == [tuple(p[:4]) for p in na.base_pairs]
    assert {p.kind for p in na.base_pairs} == {"WC", "GU"}
    assert all(len(p) == 4 for p in na.pairs)


@pytest.mark.parametrize("pdb", FIXTURES)
def test_pair_precision_against_rnapolis_is_at_least_098_where_it_was(pdb):
    # the stretched-bond second pass (C11) adds no pair RNApolis does not call; 4YBB's window held 0.977 before it
    found = _found(pdb)
    strands = [s for s, _ in found]
    ours = _named(strands, pair_nucleotides([run for _, run in found]))
    theirs = {frozenset((tuple(a), tuple(b))) for a, b, _ in REFERENCE[pdb]}
    assert len(set(ours) & theirs) / len(ours) >= (0.977 if pdb == "4YBB" else 0.98)


def _hys():
    found = _found("1HYS")
    return [s for s, _ in found], [run for _, run in found]


def test_an_in_register_pair_with_one_stretched_bond_is_accepted():
    # E868·F886: N1-N3 3.6 Å but N6-O4 3.98 Å, stacked in register on the strict pair E867·F887
    strands, runs = _hys()
    pairs = _named(strands, pair_nucleotides(runs))
    assert pairs[frozenset((("E", 868, ""), ("F", 886, "")))] == "WC"
    assert frozenset((("E", 875, ""), ("F", 879, ""))) not in pairs  # N6-O4 4.51 Å: beyond the stretch


def test_a_stretched_pair_needs_a_strict_neighbour_in_register():
    # break E867·F887's N1-N3 bond: E868·F886 then stacks on no strict pair (E869 is a bulge) and is not taken
    strands, runs = _hys()
    (e867,) = [r for r in runs[0] if r.seqid.num == 867]
    n1 = e867.find_atom("N1", "*")
    n1.pos = gemmi.Position(n1.pos.x + 3.0, n1.pos.y, n1.pos.z)
    pairs = _named(strands, pair_nucleotides(runs))
    assert frozenset((("E", 867, ""), ("F", 887, ""))) not in pairs
    assert frozenset((("E", 868, ""), ("F", 886, ""))) not in pairs
    assert frozenset((("E", 866, ""), ("F", 888, ""))) in pairs  # the strict pairs around are untouched


def _synthetic_run(names_and_bases):
    """A gemmi chain of RNA-like residues, O3'-P linked, each with the sugar atoms and the base atoms given."""
    base_atoms = {
        "A": ("N9", "N6"),
        "C": ("N1", "N4"),
        "U": ("N1", "N3", "O4"),
    }
    chain = gemmi.Chain("A")
    for k, (name, base) in enumerate(names_and_bases):
        res = gemmi.Residue()
        res.name, res.seqid = name, gemmi.SeqId(k + 1, " ")
        x = 3.5 * k
        for a, dx, dy in (
            ("P", -1.6, 0),
            ("C4'", 0, 1),
            ("O4'", 0.5, 1.2),
            ("C1'", 1, 1),
            ("O3'", 1.6, 0),
            ("O2'", 1, 2),
        ):
            atom = gemmi.Atom()
            atom.name, atom.pos = a, gemmi.Position(x + dx, dy, 0)
            res.add_atom(atom)
        for n in base_atoms[base]:
            atom = gemmi.Atom()
            atom.name, atom.pos = n, gemmi.Position(x, 5, 0)
            res.add_atom(atom)
        chain.add_residue(res)
    model = gemmi.Model("1")
    model.add_chain(chain)
    return model


def test_strand_bases_are_the_parent_bases_of_modified_and_untabulated_nucleotides():
    # CFL and UFT are 2'-F cytidine and uridine: gemmi does not know them, and their names end in L and T
    model = _synthetic_run([("A", "A"), ("CFL", "C"), ("UFT", "U"), ("A", "A")])
    ((strand, run),) = _nucleic_strands(model)
    assert [r.name for r in run] == ["A", "CFL", "UFT", "A"]
    assert strand.bases == "ACUA"


def test_strand_bases_of_a_nucleotide_without_a_base_are_n():
    model = _synthetic_run([("A", "A"), ("ZZZ", "A"), ("A", "A")])
    for atom_name in ("N9", "N6"):
        model[0][1].remove_atom(atom_name, "*")
    ((strand, _),) = _nucleic_strands(model)
    assert strand.bases == "ANA"


def test_strand_bases_follow_the_residues_through_a_cut():
    from foldmap.model import _cut, nucleic_from_pairs

    strand = _found("1EHZ")[0][0]
    assert len(strand.bases) == len(strand) == 76 and "N" not in strand.bases  # 5MC, YYG keep their parent base
    assert strand.labels[39].name == "5MC" and strand.bases[39] == "C"
    (piece,) = _cut(nucleic_from_pairs([strand], [], {}, {}), {0: [(35, 45)]}).strands
    assert piece.bases == strand.bases[35:46] and piece.bases[4] == "C"
