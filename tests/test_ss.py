from pathlib import Path

import numpy as np
import pytest

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import Bridge
from foldmap.ss import build_sses
from helpers import FID, pipeline


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


def test_bulge_split_strand_is_one_strand(zya):
    bb = load_backbone(zya, "asu")
    ids = {s.id for s in build_sses(bb, assign_dssp(bb).ss)}
    assert "A:499-509" in ids  # deposited as one strand; DSSP splits it around a bulge at 504
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
    """2HHB α: DSSP H 4-17, G 18-20, H 21-35 with the axes 96° apart (deposited HELX_P: A 3-18, B 20-35)."""
    bb, sses, _ = pipeline(Path(__file__).parent / "data" / "2HHB.cif")
    ids = {s.id for s in sses}
    for c in "AC":
        assert {f"{c}:4-17", f"{c}:21-35"} <= ids and f"{c}:4-35" not in ids
    d = assign_dssp(bb)
    with310 = build_sses(bb, d.ss, short_helices=True, bridges=d.bridges)
    assert not any(s.kind == "G" and s.chain == "A" and s.first.seq == 18 for s in with310)  # no η at the corner


def test_pi_bulge_of_a_tm_helix_stays_one_helix():
    _, sses, _ = pipeline(Path(__file__).parent / "data" / "1C3W.cif")  # helix G: I-bulge 214-218, 22°
    assert "A:201-225" in {s.id for s in sses}


@pytest.mark.parametrize("name", ["2HHB", "1C3W", "1TIM", "2LZM", "8UTF"])
def test_no_helix_hides_a_kink_of_more_than_30_degrees(name):
    import re

    from foldmap.ss import KINK_SPLIT, _fit_axis

    bb = load_backbone(Path(__file__).parent / "data" / f"{name}.cif", "asu")
    d = assign_dssp(bb)
    for s in build_sses(bb, d.ss, bridges=d.bridges):
        seg = d.ss[s.start : s.end + 1]
        hs = [(m.start() + s.start, m.end() - 1 + s.start) for m in re.finditer(r"H{4,}", seg)]
        for (x0, x1), (y0, y1) in zip(hs, hs[1:]):
            if set(d.ss[x1 + 1 : y0]) <= set("GI"):
                cos = _fit_axis(bb.ca[x0 : x1 + 1]) @ _fit_axis(bb.ca[y0 : y1 + 1])
                assert np.degrees(np.arccos(np.clip(cos, -1, 1))) <= KINK_SPLIT, s.id
