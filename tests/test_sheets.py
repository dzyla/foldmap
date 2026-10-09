import numpy as np

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import Bridge, ResLabel, SSE
from foldmap.sheets import build_sheets
from foldmap.ss import build_sses


def fake(start, end):
    lab = ResLabel("A", start, "", "ALA")
    return SSE("E", "A", start, end, lab, lab, np.zeros(3), np.array([1.0, 0, 0]))


def ladder(a, b, kind, n=3):
    """n bridges between strands a and b (SSE objects)."""
    return [Bridge(a.start + k, b.start + k, kind) for k in range(n)]


def test_path_order_and_directions():
    s1, s2, s3 = fake(0, 9), fake(20, 29), fake(40, 49)
    sheets = build_sheets([s1, s2, s3], ladder(s1, s2, "A") + ladder(s2, s3, "P"))
    (sh,) = sheets
    assert [s.start for s in sh.strands] == [0, 20, 40]
    assert sh.directions == [1, -1, -1]
    assert sh.pair_kinds == ["A", "P"]
    assert not sh.closed and not sh.ambiguous


def test_cycle_is_closed_and_cut():
    s = [fake(0, 9), fake(20, 29), fake(40, 49)]
    br = ladder(s[0], s[1], "A", 4) + ladder(s[1], s[2], "A", 4) + ladder(s[0], s[2], "A", 2)
    (sh,) = build_sheets(s, br)
    assert sh.closed and len(sh.strands) == 3
    assert [x.start for x in sh.strands] in ([0, 20, 40], [40, 20, 0])  # weakest edge (0-2) cut


def test_branch_is_ambiguous():
    c, a, b, d = fake(0, 9), fake(20, 29), fake(40, 49), fake(60, 69)
    br = ladder(c, a, "A") + ladder(c, b, "A") + ladder(c, d, "A")
    (sh,) = build_sheets([c, a, b, d], br)
    assert sh.ambiguous and len(sh.strands) == 4
    assert sorted(s.start for s in sh.strands) == [0, 20, 40, 60]


def test_weak_pairing_and_isolated_strand_make_separate_sheets():
    s1, s2, s3 = fake(0, 9), fake(20, 29), fake(40, 49)
    sheets = build_sheets([s1, s2, s3], ladder(s1, s2, "A", 1))  # 1 bridge < minimum
    assert [len(sh.strands) for sh in sheets] == [1, 1, 1]
    assert all(sh.pair_kinds == [] and sh.directions == [1] for sh in sheets)


def test_helices_ignored():
    h = fake(0, 9)
    h.kind = "H"
    assert build_sheets([h], []) == []


def _pipeline(path):
    bb = load_backbone(path, "asu")
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss)
    return bb, sses, build_sheets(sses, d.bridges)


def test_ubiquitin_mixed_sheet(ubq):
    _, _, (sh,) = _pipeline(ubq)
    ids = [s.id for s in sh.strands]
    known = ["A:12-16", "A:2-7", "A:66-71", "A:41-45", "A:48-49"]  # beta2-1-5-3-4
    assert ids in (known, known[::-1])
    kinds = sh.pair_kinds if ids == known else sh.pair_kinds[::-1]
    assert kinds == ["A", "P", "A", "A"]  # only the 1-5 pairing is parallel


def test_uromodulin_sheets(zs5):
    bb, sses, sheets = _pipeline(zs5)
    strands = [s for s in sses if s.kind == "E"]
    assert sum(len(sh.strands) for sh in sheets) == len(strands)  # every strand in exactly one sheet
    assert len({s.id for sh in sheets for s in sh.strands}) == len(strands)
    assert any({s.chain for s in sh.strands} == {"A", "D"} for sh in sheets)  # complementation
    first = [s.id for s in sheets[0].strands]  # matches deposited sheet AA1 in the mmCIF
    assert first == ["A:331-333", "A:339-345", "A:377-383", "A:369-371"]


def test_split_pieces_on_both_sides_of_one_partner_are_flagged():
    k, a, b = fake(0, 9), fake(20, 29), fake(32, 40)  # a, b: same chain, 2 residues apart
    (sh,) = build_sheets([k, a, b], ladder(k, a, "A") + ladder(k, b, "A"))
    assert sh.ambiguous  # two end-to-end pieces cannot flank their partner


def test_distant_strands_on_both_sides_are_not_flagged():
    k, a, b = fake(0, 9), fake(20, 29), fake(60, 69)
    (sh,) = build_sheets([k, a, b], ladder(k, a, "A") + ladder(k, b, "A"))
    assert not sh.ambiguous


def test_filament_sheet_with_bulge_strand_is_ordered_like_the_deposited_one(zya):
    _, _, sheets = _pipeline(zya)
    sh = next(s for s in sheets if any(x.id == "A:499-509" for x in s.strands))
    ids = [s.id for s in sh.strands]
    assert ids in (["A:563-572", "A:499-509", "A:520-521"], ["A:520-521", "A:499-509", "A:563-572"])
    assert not sh.ambiguous  # deposited order is 563-572 | 499-509 | 520-523
