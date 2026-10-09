"""Sequence view data: the chain's sequence with its secondary structure, and alignments mapped onto it."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.sequence import (
    Alignment,
    chain_track,
    column_classes,
    conservation,
    map_alignment,
    read_alignment,
)

DATA = Path(__file__).parent / "data"
UBQ_SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def test_track_holds_sequence_numbers_and_structure():
    t = chain_track(DATA / "1UBQ.cif")
    assert t.chain == "A" and t.letters == UBQ_SEQ and t.numbers[0] == 1 and t.numbers[-1] == 76
    assert len(t.ss) == len(t.letters) and set(t.ss) <= set("HEGT- ")
    labels = [e.label for e in t.elements]
    assert labels[:2] == ["A", "B"] and "α1" in labels and "η1" in labels  # the topology figure's own labels
    a = t.elements[0]
    assert a.kind == "E" and t.ss[a.start] == "E" and a.start == 1  # strand A starts at Gln2


def test_unmodelled_residues_are_marked_and_numbered():
    t = chain_track(DATA / "2LZM.cif")
    assert all(t.observed) or not t.observed[-1]
    t = chain_track(DATA / "6ZS5.cif", chain="A")
    assert t.observed[0] and t.numbers[0] == 328  # observed span by default
    full = chain_track(DATA / "6ZS5.cif", chain="A", full_sequence=True)
    assert len(full.letters) == 640 and not full.observed[0] and full.numbers[0] == 1
    assert full.ss[0] == " " and full.numbers[327] == 328


def test_alignment_formats_round_trip(tmp_path):
    fa = read_alignment(DATA / "ubq_family.fasta")
    assert len(fa.names) == 6 and len({len(r) for r in fa.rows}) == 1 and fa.names[0] == "1UBQ_A"
    clustal = tmp_path / "x.aln"
    clustal.write_text(
        "CLUSTAL W (1.83) multiple sequence alignment\n\n"
        + "".join(f"{n:<20}{r[:40]}\n" for n, r in zip(fa.names, fa.rows))
        + "\n"
        + "".join(f"{n:<20}{r[40:]}\n" for n, r in zip(fa.names, fa.rows))
    )
    assert read_alignment(clustal).rows == fa.rows
    sto = tmp_path / "x.sto"
    sto.write_text(
        "# STOCKHOLM 1.0\n" + "".join(f"{n} {r.replace('-', '.')}\n" for n, r in zip(fa.names, fa.rows)) + "//\n"
    )
    assert read_alignment(sto).rows == fa.rows


def test_bad_alignment_is_explained(tmp_path):
    p = tmp_path / "bad.fasta"
    p.write_text(">a\nMKV\n>b\nMK\n")
    with pytest.raises(ValueError, match="same length"):
        read_alignment(p)


def test_reference_found_and_columns_mapped():
    t = chain_track(DATA / "1UBQ.cif")
    aln = read_alignment(DATA / "ubq_family.fasta")
    m = map_alignment(aln, t)
    assert m.reference == 0
    assert len(m.position) == aln.width
    assert m.position[0] == 0 and m.position[35] is None and m.position[37] == 35  # the insertion: no residue
    assert [t.letters[p] for p in m.position if p is not None] == list(t.letters)


def test_reference_by_name_or_best_identity():
    t = chain_track(DATA / "1UBQ.cif")
    aln = read_alignment(DATA / "ubq_family.fasta")
    assert map_alignment(aln, t, reference="plant").reference == 2
    shuffled = Alignment(aln.names[1:] + aln.names[:1], aln.rows[1:] + aln.rows[:1])
    assert shuffled.names[map_alignment(shuffled, t).reference] in ("1UBQ_A", "insertion_variant")
    with pytest.raises(ValueError, match="nothing"):
        map_alignment(aln, t, reference="nothing")


def test_conserved_columns():
    aln = Alignment(["a", "b", "c", "d"], ["MKLV-", "MRIV-", "MKLAG", "MKIVG"])
    strict, similar = column_classes(aln, threshold=0.7)
    assert list(strict) == [True, False, False, False, False]
    assert similar[1] and similar[2] and not similar[4]  # K/R and L/I are similar; a gapped column is not


def test_conservation_scores():
    aln = Alignment(["a", "b", "c", "d"], ["MKLV-", "MRIW-", "MKLAG", "MKIVG"])
    c = conservation(aln)
    assert c[0] == pytest.approx(1.0) and np.all((c >= 0) & (c <= 1))
    assert c[0] > c[1] > c[3] and c[4] < c[0]
