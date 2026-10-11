"""Stems, runs, ties, bends and cropping of nucleic acids."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.io import _nucleic_strands, load_backbone, load_nucleic, read_model
from foldmap.model import AUTO_NT, NAStrand, ResLabel, crop_nucleic, nucleic_from_pairs
from foldmap.nucleic import (
    BasePair,
    NARun,
    bend,
    pair_nucleotides,
    secondary_structure,
    segment_stem,
    stem_axis,
    straight_duplex,
)
from foldmap.style import THEMES, resolve_style

NA = Path(__file__).parent / "data" / "na"
DATA = Path(__file__).parent / "data"


def _na(pdb):
    path = NA / f"{pdb}.cif.gz" if (NA / f"{pdb}.cif.gz").is_file() else DATA / f"{pdb}.cif"
    try:
        bb = load_backbone(path, "asu")
    except ValueError as e:  # nucleic acid only (1EHZ): no Backbone to load against until Task 6
        if "no protein residues" not in str(e):
            raise
        found = _nucleic_strands(read_model(path, "asu")[0])
        return nucleic_from_pairs([s for s, _ in found], pair_nucleotides([r for _, r in found]), {}, {})
    return load_nucleic(path, bb, "asu")


def _arms(na, st):
    (a, i0, i1), (b, j0, j1) = st.arm5, st.arm3
    la, lb = na.strands[a].labels, na.strands[b].labels
    return f"{la[i0].chain}:{la[i0].tag}-{la[i1].tag}/{lb[j0].chain}:{lb[j0].tag}-{lb[j1].tag}"


def _runs(na):
    out = []
    for r in na.structure.runs:
        labels = na.strands[r.strand].labels
        out.append((r.kind, f"{labels[r.lo].chain}:{labels[r.lo].tag}-{labels[r.hi].tag}"))
    return out


def _strands(*lengths):
    """Synthetic RNA strands A, B ... of these lengths, numbered from 1 (nucleotide index k is residue k + 1)."""
    out = []
    for chain, n in zip("ABCDEFGH", lengths):
        out.append(
            NAStrand(
                chain, [ResLabel(chain, k + 1, "", "A") for k in range(n)], np.zeros((n, 3)), np.zeros((n, 3)), True
            )
        )
    return out


def _pairs(*ij, a=0, b=0):
    return [BasePair(a, i, b, j, "WC") for i, j in ij]


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
    assert _runs(na) == [
        ("hairpin", "1:13-15"),
        ("hairpin", "1:28-32"),
        ("junction", "1:40-51"),
        ("hairpin", "1:57-61"),
        ("hairpin", "1:70-72"),
        ("junction", "1:76-76"),
        ("hairpin", "1:82-85"),
        ("junction", "1:91-111"),
        ("tail", "1:119-121"),
    ]


def test_trna_cloverleaf():
    na = _na("1EHZ")
    got = [(st.name, len(st), _arms(na, st)) for st in na.structure.stems]
    assert got == [
        ("P1", 7, "A:1-7/A:66-72"),  # acceptor
        ("P2", 4, "A:10-13/A:22-25"),  # D
        ("P3", 4, "A:27-30/A:40-43"),  # anticodon: A31·Ψ39 is not a Watson-Crick pair in this structure
        ("P4", 5, "A:49-53/A:61-65"),  # T
    ]
    assert _runs(na) == [
        ("junction", "A:8-9"),
        ("hairpin", "A:14-21"),  # D loop: holds G19 of the G19·C56 tie
        ("junction", "A:26-26"),
        ("hairpin", "A:31-39"),  # anticodon loop
        ("junction", "A:44-48"),  # variable loop
        ("hairpin", "A:54-60"),  # T loop: holds C56
        ("tail", "A:73-76"),  # CCA end
    ]


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


def test_a_bulge_nucleotide_that_is_paired_elsewhere_breaks_the_stem():
    # (10,40) -> (11,37) would skip 38 and 39 as a 2-nt bulge, but 38 pairs with 60
    pairs = _pairs((10, 40), (11, 37), (12, 36), (38, 60))
    ss = secondary_structure(_strands(70), pairs)
    assert [st.pairs for st in ss.stems] == [_pairs((11, 37), (12, 36))]
    assert ss.ties == _pairs((10, 40), (38, 60))
    stem_nts = [n for st in ss.stems for n in st.nucleotides()]
    assert len(stem_nts) == len(set(stem_nts))
    tie_nts = {(t.a, t.i) for t in ss.ties} | {(t.b, t.j) for t in ss.ties}
    assert not tie_nts & set(stem_nts)


def test_bulges_are_sorted():
    # a 1-nt bulge on the second arm (17), then one on the first arm (3)
    ss = secondary_structure(_strands(20), _pairs((0, 19), (1, 18), (2, 16), (4, 15), (5, 14)))
    assert len(ss.stems) == 1 and len(ss.stems[0]) == 5
    assert ss.stems[0].bulges == [(0, 3), (0, 17)]


def test_a_1x1_mismatch_stays_inside_the_stem():
    ss = secondary_structure(_strands(20), _pairs((0, 19), (1, 18), (3, 16), (4, 15)))
    assert [len(st) for st in ss.stems] == [4] and not ss.ties
    assert ss.stems[0].bulges == [(0, 2), (0, 17)]
    assert ss.runs == [NARun(0, 5, 14, (0, 4), (0, 15), "hairpin")]


def test_an_internal_loop_between_two_stems():
    ss = secondary_structure(_strands(30), _pairs((0, 29), (1, 28), (2, 27), (6, 23), (7, 22), (8, 21)))
    assert [st.id for st in ss.stems] == ["A:1-3/A:28-30", "A:7-9/A:22-24"]
    assert [(r.lo, r.hi, r.kind) for r in ss.runs] == [(3, 5, "internal"), (9, 20, "hairpin"), (24, 26, "internal")]


def test_a_run_no_pair_encloses_is_a_linker():
    ss = secondary_structure(_strands(22), _pairs((0, 9), (1, 8), (12, 21), (13, 20)))
    assert [(r.lo, r.hi, r.kind) for r in ss.runs] == [(2, 7, "hairpin"), (10, 11, "linker"), (14, 19, "hairpin")]
    assert ss.runs[1].left == (0, 9) and ss.runs[1].right == (0, 12)


def test_strands_without_pairs_are_one_tail_each():
    ss = secondary_structure(_strands(5, 3), [])
    assert not ss.stems and not ss.ties
    assert ss.runs == [NARun(0, 0, 4, None, None, "tail"), NARun(1, 0, 2, None, None, "tail")]


def test_a_one_nucleotide_strand_is_a_tail():
    ss = secondary_structure(_strands(4, 1, 4), _pairs((0, 3), (1, 2), a=0, b=2))
    assert ss.stems[0].id == "A:1-2/C:3-4"
    assert ss.runs == [
        NARun(0, 2, 3, (0, 1), None, "tail"),
        NARun(1, 0, 0, None, None, "tail"),
        NARun(2, 0, 1, None, (2, 2), "tail"),
    ]


@pytest.mark.parametrize("pdb", ["1LMB", "1ZAA"])
def test_straight_duplexes_stay_whole_and_keep_the_classic_frame(pdb):
    na = _na(pdb)
    (st,) = na.structure.stems
    assert len(segment_stem(st, na.strands)) == 1
    assert straight_duplex(na.structure, na.strands) == [st.a, st.b]


def test_long_duplex_axis_agrees_with_its_backbone():
    na = _na("1LMB")
    (st,) = na.structure.stems
    _, axis = stem_axis(st, na.strands)
    pts = np.vstack([s.p for s in na.strands])
    pca = np.linalg.svd(pts - pts.mean(axis=0))[2][0]
    assert abs(axis @ pca) > np.cos(np.radians(8))


def test_nucleosome_dna_follows_the_superhelix_in_straight_pieces():
    na = _na("1KX5")
    (st,) = na.structure.stems
    pieces = segment_stem(st, na.strands)
    assert len(pieces) >= 10
    assert sum(len(p) for p in pieces) == len(st)
    axes = [stem_axis(p, na.strands)[1] for p in pieces]
    turns = [np.degrees(np.arccos(np.clip(a @ b, -1, 1))) for a, b in zip(axes, axes[1:])]
    assert 500 < sum(turns) < 700  # 1.65 turns of superhelix (594 deg), not doubled back by a flipped axis
    assert max(turns) < 60  # neighbouring pieces turn gently
    assert [p.name for p in pieces[:3]] == ["DNA 1a", "DNA 1b", "DNA 1c"]
    assert all(p.parent == st.id for p in pieces)
    assert straight_duplex(na.structure, na.strands) is None  # bent: not the classic straight picture


def test_ihf_bends_its_dna_by_more_than_150_degrees():
    na = _na("1IHF")
    assert len({st.a for st in na.structure.stems}) == 1  # both stems run 5'->3' along strand C: axes comparable
    pieces = [p for st in na.structure.stems for p in segment_stem(st, na.strands)]
    assert len(pieces) >= 2
    assert bend([pieces[0], pieces[-1]], na.strands) > 150


@pytest.mark.parametrize("max_turn", [25.0, 10.0])
@pytest.mark.parametrize("pdb", ["1IHF", "1KX5", "30QM", "4OO8", "4YBB"])
def test_pieces_keep_every_nucleotide_and_bulge_of_their_stem(pdb, max_turn):
    na = _na(pdb)
    for st in na.structure.stems:
        pieces = segment_stem(st, na.strands, max_turn)
        assert sorted(n for p in pieces for n in p.nucleotides()) == sorted(st.nucleotides())
        assert sorted(x for p in pieces for x in p.bulges) == st.bulges


@pytest.mark.parametrize("pdb", ["1EHZ", "1IHF", "1KX5", "1LMB", "30QM", "4YBB"])
def test_stem_and_piece_axes_point_along_the_first_arm(pdb):
    na = _na(pdb)
    for st in na.structure.stems:
        for p in [st] + segment_stem(st, na.strands):
            a, i0, i1 = p.arm5
            c1 = na.strands[a].c1
            assert stem_axis(p, na.strands)[1] @ (c1[i1] - c1[i0]) > 0, (p.id, len(p))


def test_a_short_stem_axis_crosses_its_pairs_and_rises_with_them():
    na = _na("1EHZ")
    (st,) = [s for s in na.structure.stems if s.name == "P3"]  # 4-bp anticodon stem: the cross-vector axis
    a = np.array([na.strands[p.a].c1[p.i] for p in st.pairs])
    b = np.array([na.strands[p.b].c1[p.j] for p in st.pairs])
    _, axis = stem_axis(st, na.strands)
    across = (b - a) / np.linalg.norm(b - a, axis=1)[:, None]
    assert np.all(np.abs(across @ axis) < 0.1)  # perpendicular to every C1'-C1' line
    rise = np.diff((a + b) / 2 @ axis)
    assert np.all((2.5 < rise) & (rise < 4.5))  # each pair one stacking step further along it (measured 2.8-3.5)


def _touched(na):
    return {(na.strands[s].labels[i].chain, na.strands[s].labels[i].seq) for v in na.contacts.values() for s, i in v}


def test_auto_extent_crops_a_large_rna_to_the_touched_stems():
    na = _na("4YBB")
    assert na.size > AUTO_NT
    cut = crop_nucleic(na, by="stems")
    assert cut.size < na.size and cut.gaps
    assert _touched(cut) == _touched(na)  # nothing the protein touches is cropped away
    assert all(g.count > 0 for g in cut.gaps)
    for st in cut.structure.stems:  # stems survive whole
        assert len(st) >= 2


def test_crop_by_contacts_is_unchanged():
    na = _na("1LMB")
    cut = crop_nucleic(na, flank=0)
    assert len(cut.strands) == len(na.strands) and not cut.gaps


def test_unknown_crop_mode_is_an_error():
    with pytest.raises(ValueError, match="contacts or stems"):
        crop_nucleic(_na("1LMB"), by="everything")


def test_focus_keeps_only_the_nucleic_acid_the_subunit_holds():
    from foldmap.model import nucleic_from_pairs, nucleic_touching

    na = _na("1KX5")
    first = {(0, i) for i in range(10)}
    lone = nucleic_from_pairs(na.strands, [], {7: first}, {7: "A"})  # unpaired, touched only on strand I by A
    kept = nucleic_touching(lone, ["A"])
    assert [s.chain for s in kept.strands] == ["I"]
    assert nucleic_touching(lone, ["B"]).strands == []
    both = nucleic_touching(nucleic_from_pairs(na.strands, na.base_pairs, {7: first}, {7: "A"}), ["A"])
    assert [s.chain for s in both.strands] == ["I", "J"]  # the paired strand comes along


def test_theme_defaults_for_stems():
    assert THEMES["publication"].na_style == "ladder"
    assert THEMES["minimal"].na_style == THEMES["print"].na_style == "block"
    assert THEMES["presentation"].na_style == THEMES["cartoon"].na_style == "ribbon"


def test_dna_extent_is_still_accepted():
    assert resolve_style(None, None, ["dna_extent=all"]).na_extent == "all"
    with pytest.raises(ValueError, match="na_extent"):
        resolve_style(None, None, ["na_extent=some"])
