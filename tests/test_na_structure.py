"""Stems, runs, ties, bends and cropping of nucleic acids."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.io import load_backbone, load_nucleic
from foldmap.model import AUTO_NT, NAStrand, ResLabel, crop_nucleic, nucleic_from_pairs
from foldmap.nucleic import (
    BasePair,
    NARun,
    bend,
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
    return load_nucleic(path, load_backbone(path, "asu"), "asu")


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
    cut_stems = [_pair_labels(cut, st.pairs) for st in cut.structure.stems]
    touched = {(s, i) for v in na.contacts.values() for s, i in v}
    for st in na.structure.stems:  # every stem the protein touches survives whole
        if touched & set(st.nucleotides()):
            assert _pair_labels(na, st.pairs) in cut_stems


def _label(na, s, i):
    return na.strands[s].labels[i][:3]


def _pair_labels(na, pairs):
    return {frozenset((_label(na, p.a, p.i), _label(na, p.b, p.j))) for p in pairs}


def test_a_stem_crop_draws_no_paired_nucleotide_single_stranded_and_no_half_stem():
    na = _na("4YBB")
    cut = crop_nucleic(na, by="stems")
    whole = [_pair_labels(na, st.pairs) for st in na.structure.stems]
    for st in cut.structure.stems:  # each stem drawn is an original stem with all its pairs
        assert _pair_labels(cut, st.pairs) in whole
    paired = {x for p in na.base_pairs for x in (_label(na, p.a, p.i), _label(na, p.b, p.j))}
    tie_ends = {x for p in cut.structure.ties for x in (_label(cut, p.a, p.i), _label(cut, p.b, p.j))}
    for r in cut.structure.runs:
        for i in range(r.lo, r.hi + 1):
            x = _label(cut, r.strand, i)
            assert x not in paired or x in tie_ends, x


def test_stem_crop_takes_neighbour_stems_on_the_same_strand_only():
    def strand(chain):
        n = 40
        labels = [ResLabel(chain, k + 1, "", "G" if k % 2 else "C") for k in range(n)]
        return NAStrand(chain, labels, np.zeros((n, 3)), np.zeros((n, 3)), True)

    def hairpin(s, i):  # two stacked pairs closing a 6-nt loop: i-(i+9), (i+1)-(i+8)
        return [BasePair(s, i, s, i + 9, "WC"), BasePair(s, i + 1, s, i + 8, "WC")]

    pairs = hairpin(0, 0) + hairpin(0, 20) + hairpin(1, 0) + hairpin(1, 20)
    na = nucleic_from_pairs([strand("R"), strand("S")], pairs, {5: {(0, 20)}}, {5: "A"})
    assert [(st.a, st.pairs[0].i) for st in na.structure.stems] == [(0, 0), (0, 20), (1, 0), (1, 20)]
    cut = crop_nucleic(na, by="stems")
    assert {s.chain for s in cut.strands} == {"R"}  # strand S's first stem is next in the list, not in sequence
    assert len(cut.structure.stems) == 2


def test_focus_on_a_subunit_keeps_only_its_contacts_and_crops_to_its_stems():
    from foldmap.cli import make_layout

    full, _, _ = make_layout(NA / "4YBB.cif.gz", focus="none")
    lay, _, _ = make_layout(NA / "4YBB.cif.gz", focus="AD")
    assert set(lay.nucleic.contact_chain.values()) <= {"AD"}
    assert lay.nucleic.size < full.nucleic.size  # auto crops to the stems the focused subunit touches


def test_focus_drops_contacts_of_chains_not_placed():
    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(NA / "1KX5.cif.gz", focus="A")
    assert lay.nucleic.contacts
    assert set(lay.nucleic.contact_chain.values()) <= {s.chain for s in sses}


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


def _duplex(n, splay=()):
    """Two antiparallel DNA strands A, B of n nucleotides, B[n-1-k] paired across from A[k] (C1' 10.5 Å apart, 3.4 Å
    rise); the nucleotides of B in `splay` swung 8 Å further out (not in the helix)."""
    a = np.array([[0.0, 0.0, 3.4 * k] for k in range(n)])
    b = np.array([[10.5 + (8.0 if n - 1 - k in splay else 0.0), 0.0, 3.4 * (n - 1 - k)] for k in range(n)])
    return [
        NAStrand(c, [ResLabel(c, k + 1, "", "A") for k in range(n)], xyz, xyz, False) for c, xyz in (("A", a), ("B", b))
    ]


@pytest.mark.parametrize("gap", [1, 2])
def test_a_duplex_stays_whole_across_unpaired_steps_in_register(gap):
    # A[k]·B[n-1-k] paired, but A[5..5+gap) unpaired on both strands: the same helix either side (1AOI I20-21)
    n = 14
    pairs = [BasePair(0, k, 1, n - 1 - k, "WC") for k in range(n) if not 5 <= k < 5 + gap]
    ss = secondary_structure(_duplex(n), pairs)
    assert [len(st) for st in ss.stems] == [n - gap] and not ss.ties


def test_a_two_by_two_loop_whose_steps_leave_the_helix_splits_the_duplex():
    n = 14
    pairs = [BasePair(0, k, 1, n - 1 - k, "WC") for k in range(n) if k not in (5, 6)]
    ss = secondary_structure(_duplex(n, splay=(5, 6)), pairs)
    assert [len(st) for st in ss.stems] == [5, 7]


def test_the_1aoi_nucleosome_dna_is_one_duplex_cut_only_where_it_bends():
    na = _na("1AOI")
    (st,) = na.structure.stems  # I20-21 and J272-273 unpaired in register: inside the stem
    assert len(st) >= 140 and st.name == "DNA 1"
    assert {(0, 20 - 1), (0, 21 - 1)} <= set(st.bulges)
    pieces = segment_stem(st, na.strands)
    assert len(pieces) > 1 and all(p.parent == st.id for p in pieces)


def test_the_sam_i_riboswitch_is_named_as_in_the_literature():
    # Roy 2019 Fig 1A: P1 1-8/86-94, P2 13-41 (kink-turn inside), P3 42-60, P4 69-82, and the pseudoknot between
    # the P2 loop and J3/4 (25-28·65-68)
    ss = _na("2GIS").structure
    assert [(st.name, st.id) for st in ss.stems] == [
        ("P1", "A:1-8/A:86-93"),
        ("P2a", "A:13-17/A:38-41"),
        ("P2b", "A:21-23/A:29-31"),
        ("PK", "A:25-28/A:65-68"),
        ("P3a", "A:42-44/A:58-60"),
        ("P3b", "A:47-49/A:54-56"),
        ("P4", "A:69-73/A:78-82"),
    ]


def test_a_helix_continued_across_an_internal_loop_keeps_its_number():
    ss = secondary_structure(_strands(30), _pairs((0, 29), (1, 28), (2, 27), (6, 23), (7, 22), (8, 21)))
    assert [st.name for st in ss.stems] == ["P1a", "P1b"]


def test_stems_branching_from_a_junction_are_numbered_apart():
    pairs = _pairs((0, 39), (1, 38), (4, 14), (5, 13), (20, 32), (21, 31))
    ss = secondary_structure(_strands(40), pairs)
    assert [st.name for st in ss.stems] == ["P1", "P2", "P3"]


@pytest.mark.parametrize("first, names", [(3, ["P1", "PK"]), (2, ["PK", "P1"])])
def test_a_pseudoknot_is_the_shorter_of_two_crossing_stems_or_the_later_one(first, names):
    # a hairpin from 0 crossed by the 3-bp stem 10-12/26-28: of equal stems the 5' one nests, a shorter one is the knot
    hairpin = _pairs(*[(k, 20 - k) for k in range(first)])
    knot = _pairs((10, 28), (11, 27), (12, 26))
    ss = secondary_structure(_strands(30), hairpin + knot)
    assert [st.name for st in ss.stems] == names
    kinds = {r.lo: r.kind for r in ss.runs}
    assert kinds[first] == ("hairpin" if first == 3 else "linker")  # the loop the nested stem closes


def _bent(n, at, degrees):
    """A synthetic n-bp duplex (strands A, B) whose axis turns by `degrees` about x after pair `at`."""
    t = np.radians(degrees)
    turn = np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
    a, b, mid = np.zeros((n, 3)), np.zeros((n, 3)), np.zeros(3)
    for k in range(n):
        frame = turn if k >= at else np.eye(3)
        u, e1, e2 = frame @ [0, 0, 1], frame @ [1, 0, 0], frame @ [0, 1, 0]
        if k:
            mid = mid + 3.4 * u
        r = 5.25 * (np.cos(np.radians(36 * k)) * e1 + np.sin(np.radians(36 * k)) * e2)
        a[k], b[n - 1 - k] = mid - r, mid + r
    return [
        NAStrand(c, [ResLabel(c, k + 1, "", "A") for k in range(n)], xyz, xyz, False) for c, xyz in (("A", a), ("B", b))
    ]


@pytest.mark.parametrize("name, marks", [("DNA 1", "abc"), ("P2a", "123")])
def test_pieces_of_a_bent_stem_are_lettered_or_numbered_after_a_letter(name, marks):
    strands = _bent(20, 10, 60)
    (st,) = secondary_structure(strands, [BasePair(0, k, 1, 19 - k, "WC") for k in range(20)]).stems
    st.name = name
    pieces = segment_stem(st, strands)
    assert len(pieces) > 1 and [p.name for p in pieces] == [name + m for m in marks[: len(pieces)]]
