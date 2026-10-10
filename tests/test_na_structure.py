"""Stems, runs, ties, bends and cropping of nucleic acids."""

from pathlib import Path

from foldmap.io import _nucleic_strands, load_backbone, load_nucleic, read_model
from foldmap.model import nucleic_from_pairs
from foldmap.nucleic import pair_nucleotides

NA = Path(__file__).parent / "data" / "na"
DATA = Path(__file__).parent / "data"


def _na(pdb):
    path = NA / f"{pdb}.cif.gz" if (NA / f"{pdb}.cif.gz").is_file() else DATA / f"{pdb}.cif"
    try:
        bb = load_backbone(path, "asu")
    except ValueError:  # nucleic acid only (1EHZ): no Backbone to load against until Task 6
        found = _nucleic_strands(read_model(path, "asu")[0])
        return nucleic_from_pairs([s for s, _ in found], pair_nucleotides([r for _, r in found]), {}, {})
    return load_nucleic(path, bb, "asu")


def _arms(na, st):
    (a, i0, i1), (b, j0, j1) = st.arm5, st.arm3
    la, lb = na.strands[a].labels, na.strands[b].labels
    return f"{la[i0].chain}:{la[i0].tag}-{la[i1].tag}/{lb[j0].chain}:{lb[j0].tag}-{lb[j1].tag}"


def test_30qm_ncrna_has_its_six_stems():
    na = _na("30QM")
    got = {st.name: (len(st), _arms(na, st)) for st in na.structure.stems}
    assert got == {
        "P1": (7, "1:1-7/1:112-118"),
        "P2": (5, "1:8-12/1:16-20"),
        "P3": (7, "1:21-27/1:33-39"),
        "P4": (5, "1:52-56/1:62-66"),
        "P5": (3, "1:67-69/1:73-75"),
        "P6": (5, "1:77-81/1:86-90"),
    }
    kinds = sorted(r.kind for r in na.structure.runs)
    assert kinds.count("hairpin") == 5 and kinds.count("junction") == 3 and kinds.count("tail") == 1


def test_trna_cloverleaf():
    na = _na("1EHZ")
    got = [(st.name, len(st), _arms(na, st)) for st in na.structure.stems]
    assert got == [
        ("P1", 7, "A:1-7/A:66-72"),  # acceptor
        ("P2", 4, "A:10-13/A:22-25"),  # D
        ("P3", 4, "A:27-30/A:40-43"),  # anticodon: A31·Ψ39 is not a Watson-Crick pair in this structure
        ("P4", 5, "A:49-53/A:61-65"),  # T
    ]
    assert sorted(r.kind for r in na.structure.runs).count("hairpin") == 3


def test_cas9_heteroduplex_is_one_stem_across_two_strands():
    na = _na("4OO8")
    first = na.structure.stems[0]
    assert (len(first), _arms(na, first)) == (20, "B:1-20/C:1-20")
    assert na.strands[first.a].rna and not na.strands[first.b].rna
    assert first.name == "P1"  # RNA in it: named like an RNA stem


def test_holliday_junction_has_four_arms_on_four_strands():
    na = _na("3CRX")
    stems = na.structure.stems
    assert len(stems) == 4
    assert {frozenset((na.strands[s.a].chain, na.strands[s.b].chain)) for s in stems} == {
        frozenset("CD"),
        frozenset("CF"),
        frozenset("DE"),
        frozenset("EF"),
    }
    assert [s.name for s in stems] == ["DNA 1", "DNA 2", "DNA 3", "DNA 4"]


def test_an_isolated_pair_is_a_tie_not_a_stem():
    na = _na("1EHZ")
    assert len(na.structure.ties) == 1
    tie = na.structure.ties[0]
    assert not any((tie.a, tie.i) in set(st.nucleotides()) for st in na.structure.stems)


def test_stem_ids_are_unique_and_name_their_arms():
    na = _na("4YBB")
    ids = [st.id for st in na.structure.stems]
    assert len(ids) == len(set(ids))
    assert ids[0].startswith("AA:") and "/" in ids[0]
