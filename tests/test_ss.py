import numpy as np

from topoplot.dssp import assign_dssp
from topoplot.io import load_backbone
from topoplot.ss import build_sses


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


def test_strands_one_residue_apart_merge_two_apart_do_not(ubq):
    bb = load_backbone(ubq)
    one = list("-" * len(bb))
    one[10:14] = "EEEE"
    one[15:19] = "EEEE"  # one residue (index 14) between
    (merged,) = build_sses(bb, "".join(one))
    assert (merged.start, merged.end) == (10, 18)
    two = list("-" * len(bb))
    two[10:14] = "EEEE"
    two[16:20] = "EEEE"  # a 2-residue gap can be a hairpin turn
    assert [(s.start, s.end) for s in build_sses(bb, "".join(two))] == [(10, 13), (16, 19)]


def test_merge_never_crosses_chain_break(zs5):
    bb = load_backbone(zs5)
    i = next(k for k in range(1, len(bb)) if not bb.prev[k])  # first break
    ss = list("-" * len(bb))
    ss[i - 3 : i] = "EEE"
    ss[i + 1 : i + 4] = "EEE"  # break sits between, one residue (i) in the gap
    for s in build_sses(bb, "".join(ss)):
        assert all(bb.prev[k] for k in range(s.start + 1, s.end + 1))
