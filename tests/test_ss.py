from pathlib import Path

import numpy as np
import pytest

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import Backbone, Bridge, ResLabel
from foldmap.ss import build_sses
from helpers import FID, pipeline

NA = Path(__file__).parent / "data" / "na"


def test_ubiquitin_sses(ubq):
    bb = load_backbone(ubq)
    sses = build_sses(bb, assign_dssp(bb).ss)
    assert [(s.kind, s.id) for s in sses] == [
        ("E", "A:2-7"),
        ("E", "A:12-16"),
        ("H", "A:23-34"),
        ("E", "A:41-45"),
        ("E", "A:48-49"),
        ("E", "A:66-71"),
    ]  # the 3-residue 3-10 turn at 38-40 is not drawn


def test_axis_and_centroid(ubq):
    bb = load_backbone(ubq)
    for s in build_sses(bb, assign_dssp(bb).ss):
        assert abs(np.linalg.norm(s.axis) - 1) < 1e-9
        span = bb.ca[s.end] - bb.ca[s.start]
        assert np.dot(span, s.axis) > 0  # points N -> C
        assert np.allclose(s.centroid, bb.ca[s.start : s.end + 1].mean(axis=0))


def test_sses_never_span_chain_breaks(zs5):
    bb = load_backbone(zs5)
    for ss in ("E" * len(bb), "H" * len(bb)):
        for s in build_sses(bb, ss):
            assert all(bb.labels[i].chain == s.chain for i in range(s.start, s.end + 1))
            assert all(bb.prev[i] for i in range(s.start + 1, s.end + 1))


def test_short_runs_filtered(ubq):
    bb = load_backbone(ubq)
    ss = list("-" * len(bb))
    ss[5] = "E"  # lone strand residue
    ss[10:13] = "GGG"  # 3-10 turn too short
    ss[20:23] = "HHH"  # alpha run too short
    assert build_sses(bb, "".join(ss)) == []


def test_a_classic_bulge_stays_inside_one_dssp_strand(zya):
    """6ZYA A:499-509, deposited as one strand, pairs with 563-572 except at 504 (503-568, 505-567): a classic
    β-bulge, which DSSP's own bulge rule labels E. It is one strand with no join (no bridges passed), and stays one
    with them."""
    bb = load_backbone(zya, "asu")
    d = assign_dssp(bb)
    for bridges in (None, d.bridges):
        ids = {s.id for s in build_sses(bb, d.ss, bridges=bridges)}
        assert "A:499-509" in ids
        assert "A:499-503" not in ids and "A:505-509" not in ids


def _same_side(x, y):
    """Three bridges each from strands starting at x and y to residues 51-49 and 47-45 of one partner strand."""
    return [Bridge(x + t, 51 - t, "A") for t in range(3)] + [Bridge(y + t, 47 - t, "A") for t in range(3)]


def test_strands_one_residue_apart_join_only_with_their_bridges(ubq):
    bb = load_backbone(ubq)
    one = list("-" * len(bb))
    one[10:14] = "EEEE"
    one[15:19] = "EEEE"  # one residue (index 14) between
    one[40:52] = "E" * 12  # their partner
    pieces = lambda bridges: [(s.start, s.end) for s in build_sses(bb, "".join(one), bridges=bridges)]  # noqa: E731
    assert pieces(_same_side(10, 15)) == [(10, 18), (40, 51)]
    assert pieces(None) == [(10, 13), (15, 18), (40, 51)]  # without bridges, strands stay as DSSP runs


def test_merge_never_crosses_chain_break(zs5, ubq):
    bb = load_backbone(zs5)
    i = next(k for k in range(1, len(bb)) if not bb.prev[k])  # first break
    ss = list("-" * len(bb))
    ss[i - 3 : i] = "EEE"
    ss[i + 1 : i + 4] = "EEE"  # break sits between, one residue (i) in the gap
    ss[40:52] = "E" * 12  # their partner, before the break
    for s in build_sses(bb, "".join(ss), bridges=_same_side(i - 3, i + 1)):
        assert all(bb.prev[k] for k in range(s.start + 1, s.end + 1))
    bb = load_backbone(ubq)
    bb.prev[15] = False  # a break inside one chain, at the end of the one-residue gap joined above
    one = list("-" * len(bb))
    one[10:14], one[15:19], one[40:52] = "EEEE", "EEEE", "E" * 12
    pieces = [(s.start, s.end) for s in build_sses(bb, "".join(one), bridges=_same_side(10, 15))]
    assert pieces == [(10, 13), (15, 18), (40, 51)]


def test_strands_one_residue_apart_join_only_when_they_share_a_partner():
    bb = load_backbone(Path(__file__).parent / "data" / "6ZS5.cif", "asu")
    d = assign_dssp(bb)
    ids = {s.id for s in build_sses(bb, d.ss, bridges=d.bridges)}
    assert {"A:419-428", "A:430-434"} <= ids and "A:419-434" not in ids  # each half pairs with another sheet
    r74 = load_backbone(Path(__file__).parent / "data" / "6R74.cif", "asu")
    d = assign_dssp(r74)
    ids = {s.id for s in build_sses(r74, d.ss, bridges=d.bridges)}
    assert {"A:114-123", "A:146-158"} <= ids  # 146-147 and 150-158 pair with 3-14 in register, across a 2-residue bulge
    assert "A:127-144" in ids  # pairs with 52-61 and 97-105 only (147-95 is one bridge): another sheet, as in 5NKT


def test_tenascin_g_strand_split_by_a_bulge_is_one_strand():
    """1TEN: DSSP breaks G at a 3-residue bulge (880-881 | 885-890); both pieces pair with F on the same side."""
    _, sses, _ = pipeline(FID / "1TEN.cif.gz")
    ids = [s.id for s in sses if s.kind == "E"]
    assert "A:880-890" in ids and "A:880-881" not in ids and "A:885-890" not in ids
    assert len(ids) == 7  # A B C C' E F G


def test_bulge_join_follows_the_partner_register():
    bb = load_backbone(Path(__file__).parent / "data" / "1UBQ.cif")
    ss = ["-"] * len(bb)
    ss[10:14] = "EEEE"  # x: 10-13
    ss[17:21] = "EEEE"  # y: 17-20, three residues (14-16) later
    ss[40:52] = "E" * 12  # k: 40-51
    same_side = [Bridge(10 + t, 51 - t, "A") for t in range(3)] + [Bridge(17 + t, 46 - t, "A") for t in range(3)]
    joined = build_sses(bb, "".join(ss), bridges=same_side)
    assert [(s.start, s.end) for s in joined if s.kind == "E"] == [(10, 20), (40, 51)]


def test_hairpin_and_opposite_sides_are_not_joined():
    bb = load_backbone(Path(__file__).parent / "data" / "1UBQ.cif")
    ss = ["-"] * len(bb)
    ss[10:14], ss[17:21], ss[40:52] = "EEEE", "EEEE", "E" * 12
    hairpin = [Bridge(10 + t, 20 - t, "A") for t in range(4)]
    opposite = [Bridge(10 + t, 45 - t, "A") for t in range(3)] + [Bridge(17 + t, 45 - t, "A") for t in range(3)]
    for bridges in (hairpin, opposite):  # opposite: k residues 43-45 pair with both, one on each face
        pieces = [(s.start, s.end) for s in build_sses(bb, "".join(ss), bridges=bridges) if s.kind == "E"]
        assert (10, 13) in pieces and (17, 20) in pieces


def test_globin_ab_corner_splits_the_alpha_chain_helix():
    """2HHB α: DSSP H 4-17, G 18-20, H 21-35 with the local axes 105° apart (deposited HELX_P: A 3-18, B 20-35)."""
    bb, sses, _ = pipeline(Path(__file__).parent / "data" / "2HHB.cif")
    ids = {s.id for s in sses}
    for c in "AC":
        assert {f"{c}:4-17", f"{c}:21-35"} <= ids and f"{c}:4-35" not in ids
    d = assign_dssp(bb)
    with310 = build_sses(bb, d.ss, short_helices=True, bridges=d.bridges)
    assert not any(s.kind == "G" and s.chain == "A" and s.first.seq == 18 for s in with310)  # no η at the corner


def test_pi_bulge_of_a_tm_helix_stays_one_helix():
    _, sses, _ = pipeline(Path(__file__).parent / "data" / "1C3W.cif")  # helix G: 201-213, I 214-218, 219-225, 9°
    assert "A:201-225" in {s.id for s in sses}


def _ideal_helix(n: int, bend_after: int | None = None, bend: float = 0.0) -> Backbone:
    """CAs of an ideal alpha helix (radius 2.3 Å, rise 1.5 Å, 100° a residue) along z; residues after bend_after
    turned by `bend` degrees about the x axis through that residue's point on the helix axis."""
    t = np.radians(100.0) * np.arange(n)
    ca = np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)], axis=1)
    if bend_after is not None:
        c, s = np.cos(np.radians(bend)), np.sin(np.radians(bend))
        pivot = np.array([0.0, 0.0, 1.5 * bend_after])
        ca[bend_after + 1 :] = (ca[bend_after + 1 :] - pivot) @ np.array([[1, 0, 0], [0, c, -s], [0, s, c]]).T + pivot
    labels = [ResLabel("A", i + 1, "", "ALA") for i in range(n)]
    return Backbone(labels, np.repeat(ca[:, None, :], 4, axis=1), np.arange(n) > 0)


@pytest.mark.parametrize("arm", [4, 5, 8, 12])
def test_a_pi_bulge_in_a_straight_helix_is_no_kink(arm):
    """A straight fit through 4-5 CAs is off by up to 40°, so short arms of a straight helix looked kinked."""
    ss = "H" * arm + "IIIII" + "H" * arm
    assert [s.id for s in build_sses(_ideal_helix(len(ss)), ss)] == [f"A:1-{len(ss)}"]


@pytest.mark.parametrize("bend, pieces", [(15.0, ["A:1-25"]), (45.0, ["A:1-10", "A:16-25"])])
def test_a_bent_helix_splits_above_kink_split_only(bend, pieces):
    ss = "H" * 10 + "IIIII" + "H" * 10
    assert [s.id for s in build_sses(_ideal_helix(len(ss), 12, bend), ss)] == pieces


def test_a_pi_bulge_with_short_arms_in_a_tm_helix_stays_one_helix():
    """1PV7 B (window 70-115): H 87-90, I 91-95, H 96-100. A straight fit through the CAs puts the two arms 58°
    apart; their local axes are 8° apart, and arms of 4-5 residues are too short to measure a kink by anyway."""
    from foldmap.ss import KINK_SPLIT, _angle, _helix_axis

    bb, sses, _ = pipeline(FID / "1PV7.cif.gz")
    seq = {l.seq: k for k, l in enumerate(bb.labels)}
    assert assign_dssp(bb).ss[seq[87] : seq[100] + 1] == "HHHHIIIIIHHHHH"
    assert "B:87-100" in {s.id for s in sses}  # one helix
    x, y = bb.ca[seq[87] : seq[90] + 1], bb.ca[seq[96] : seq[100] + 1]
    assert _angle(_helix_axis(x), _helix_axis(y)) < KINK_SPLIT


@pytest.mark.parametrize(
    "path, splits",
    [
        ("2HHB.cif", {"A:4-17 A:21-35", "C:4-17 C:21-35"}),  # the globin AB corner, 105° and 104°
        ("1C3W.cif", set()),  # helix G's π-bulge, 9°
        ("8UTF.cif", set()),  # a/b/c 142-181 + 187-222 across a π-bulge, 11-12°
        ("1TIM.cif", set()),
        ("2LZM.cif", set()),
        ("fidelity/1MBN.cif.gz", set()),
        ("fidelity/1PV7.cif.gz", set()),
    ],
)
def test_helix_runs_are_split_at_these_kinks_only(path, splits):
    bb, sses, _ = pipeline(Path(__file__).parent / "data" / path)
    ss = assign_dssp(bb).ss
    hs = [s for s in sses if s.kind == "H"]
    cut = {
        f"{x.id} {y.id}"
        for x, y in zip(hs, hs[1:])
        if x.chain == y.chain and y.start > x.end + 1 and set(ss[x.end + 1 : y.start]) <= set("GI")
    }
    assert cut == splits


def test_copies_of_one_protein_are_found_whatever_their_numbering():
    from foldmap.model import homolog_chains

    assert homolog_chains(load_backbone(NA / "1EYG.cif.gz")) == {
        "A": ("A", 0),
        "B": ("A", 1000),
        "C": ("A", 2000),
        "D": ("A", 3000),
    }
    assert homolog_chains(load_backbone(Path(__file__).parent / "data" / "2HHB.cif", "asu")) == {
        "A": ("A", 0),
        "B": ("B", 0),
        "C": ("A", 0),
        "D": ("B", 0),
    }  # α2β2: the α chains are copies, and the β chains


def test_a_bulge_joined_in_one_copy_is_joined_in_every_copy():
    """1EYG SSB: β5 is split at a bulge in all four chains, 4 residues apart in A-C (joined) but 5 in D."""
    _, sses, _ = pipeline(NA / "1EYG.cif.gz")
    counts = {c: sum(s.chain == c for s in sses) for c in "ABCD"}
    assert len(set(counts.values())) == 1, counts
    ids = {s.id for s in sses}
    assert "D:4096-4111" in ids and "D:4096-4102" not in ids


def test_a_kink_split_in_one_copy_is_split_in_every_copy(monkeypatch):
    """2HHB α: the AB corner turns by 105.6° in chain A and 104.0° in chain C. With the split angle between the
    two, chain C is split as chain A is."""
    import foldmap.ss

    monkeypatch.setattr(foldmap.ss, "KINK_SPLIT", 104.8)
    _, sses, _ = pipeline(Path(__file__).parent / "data" / "2HHB.cif")
    ids = {s.id for s in sses}
    assert {"A:4-17", "A:21-35", "C:4-17", "C:21-35"} <= ids
    monkeypatch.setattr(foldmap.ss, "KINK_SPLIT", 106.0)  # neither copy is split
    _, sses, _ = pipeline(Path(__file__).parent / "data" / "2HHB.cif")
    assert {"A:4-35", "C:4-35"} <= {s.id for s in sses}


def test_a_helix_broken_in_one_copy_is_broken_in_every_copy():
    """7OKA (random-100 A15 1.12): DSSP broke helix 225-238 at residue 232 in chains B and C but not in A, so B's
    and C's two pieces were each tied to A's whole helix, took its length and stood one above the other. A break of
    one or two residues in one copy of a protein is made in its other copies too, at the same residue numbers, as a
    kink split is (2HHB: chain C's helix 20-35 broken at 28 by hand; chain A follows; the beta chains do not)."""
    bb = load_backbone(Path(__file__).parent / "data" / "2HHB.cif", "asu")
    ss = list(assign_dssp(bb).ss)
    at = {(l.chain, l.seq): k for k, l in enumerate(bb.labels)}
    whole = {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))}
    helix = next(h for h in sorted(whole) if h[0] == "C" and h[1] <= 28 <= h[2])
    ss[at[("C", 28)]] = "T"
    got = {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))}
    for c in "AC":
        assert (c, helix[1], 27) in got and (c, 29, helix[2]) in got, sorted(h for h in got if h[0] == c)[:6]
    assert {h for h in got if h[0] in "BD"} == {h for h in whole if h[0] in "BD"}


def test_a_helix_broken_in_one_copy_of_three_is_joined_again():
    """The majority decides: bacteriorhodopsin's trimer (1C3W, assembly 1: three copies) with one copy's helix broken
    by hand at one residue is joined there again, as the two others have it; their helices are not cut (one noisy
    copy of many, as 7JU4's tubulins, does not split every copy)."""
    bb = load_backbone(Path(__file__).parent / "data" / "1C3W.cif", "auto")
    ss = list(assign_dssp(bb).ss)
    plain = {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))}
    chains = sorted({c for c, _, _ in plain})
    assert len(chains) == 3
    c, a, b = next(h for h in sorted(plain) if h[0] == chains[0] and h[2] - h[1] >= 20)
    k = next(i for i, l in enumerate(bb.labels) if l.chain == c and l.seq == (a + b) // 2)
    ss[k] = "T"
    assert {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))} == plain


def test_helix_pieces_matched_across_copies_are_told(caplog):
    """Fix round 1 ruling: the copy-majority cuts and joins are told: build_sses counts them (`matched`), logs
    "helix pieces matched across copies: n cut, n joined" at INFO, and the figure's scope carries the counts (7TIM:
    B's 178-203 cut like A's; 1UBQ, one chain: nothing)."""
    import logging

    from foldmap.cli import make_layout
    from foldmap.io import read_model
    from foldmap.scope import drawn_scope

    bb = load_backbone(FID / "7TIM.cif.gz")
    matched: dict = {}
    with caplog.at_level(logging.INFO, logger="foldmap.ss"):
        build_sses(bb, assign_dssp(bb).ss, matched=matched)
    assert matched == {"cut": 1, "joined": 0}
    assert "helix pieces matched across copies: 1 cut, 0 joined" in caplog.text
    lay, _, _ = make_layout(FID / "7TIM.cif.gz")
    assert drawn_scope(read_model(FID / "7TIM.cif.gz"), lay)["helix_matched"] == {"cut": 1, "joined": 0}
    plain, _, _ = make_layout(Path(__file__).parent / "data" / "1UBQ.cif")
    assert drawn_scope(read_model(Path(__file__).parent / "data" / "1UBQ.cif"), plain)["helix_matched"] == {
        "cut": 0,
        "joined": 0,
    }


def test_a_helix_broken_twice_in_one_copy_of_three_is_joined_whole():
    """Fix round 1 (minor): a run broken twice by joined breaks (x|y|z) is joined into one, not just x+y."""
    bb = load_backbone(Path(__file__).parent / "data" / "1C3W.cif", "auto")
    ss = list(assign_dssp(bb).ss)
    plain = {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))}
    chains = sorted({c for c, _, _ in plain})
    c, a, b = next(h for h in sorted(plain) if h[0] == chains[0] and h[2] - h[1] >= 20)
    for seq in (a + 7, a + 14):
        ss[next(i for i, l in enumerate(bb.labels) if l.chain == c and l.seq == seq)] = "T"
    assert {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))} == plain


def test_a_join_and_a_cut_in_one_structure():
    """Fix round 1 campaign (7JU4, 7QUP crashed: the join's loop variable shadowed the residue index): in 1C3W's
    trimer, one copy's helix broken by hand is joined again, and a later helix broken in two copies of three is cut
    in the third at the same residue."""
    bb = load_backbone(Path(__file__).parent / "data" / "1C3W.cif", "auto")
    ss = list(assign_dssp(bb).ss)
    plain = {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))}
    chains = sorted({c for c, _, _ in plain})
    first, second = [h for h in sorted(plain, key=lambda h: h[1]) if h[0] == chains[0] and h[2] - h[1] >= 20][:2]

    def broken(chain, seq):
        ss[next(i for i, l in enumerate(bb.labels) if l.chain == chain and l.seq == seq)] = "T"

    broken(chains[0], (first[1] + first[2]) // 2)  # one copy of three: joined
    cut = (second[1] + second[2]) // 2
    for c in chains[1:]:  # two copies of three: cut in the third too
        broken(c, cut)
    got = {(s.chain, s.first.seq, s.last.seq) for s in build_sses(bb, "".join(ss))}
    assert first in got
    for c in chains:
        assert (c, second[1], cut - 1) in got and (c, cut + 1, second[2]) in got, c
