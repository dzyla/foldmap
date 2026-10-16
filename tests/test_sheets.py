from pathlib import Path

import numpy as np
import pytest

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import SSE, Bridge, ResLabel
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


def test_a_branch_strand_is_kept_off_the_track_beside_its_parent():
    """Plan 12 (the review's R1): a strand off the pairing tree's longest path was inserted right after its parent,
    between the parent and a strand the parent pairs with. It is now an accessory strand, drawn beside its parent."""
    c, a, b, d = fake(0, 9), fake(20, 29), fake(40, 49), fake(60, 69)
    br = ladder(c, a, "A") + ladder(c, b, "A") + ladder(c, d, "A")
    (sh,) = build_sheets([c, a, b, d], br)
    assert sh.ambiguous and [s.start for s in sh.strands] == [20, 0, 40]  # a c b: each pairs with the next
    ((acc, parent, column, direction, _),) = sh.accessory
    assert acc is d and parent is c and abs(column - 1) == 1 and direction == 1  # a c d: + - +; beside c
    assert sorted(s.start for s in [*sh.strands, acc]) == [0, 20, 40, 60]


def _ring_sheet(n, extra=()):
    """n strands paired round a ring (6 bridges each, the last pair 3), plus (strand, partner index) pendants."""
    ring = [fake(20 * k, 20 * k + 9) for k in range(n)]
    for k, s in enumerate(ring):
        s.centroid = np.array([np.cos(2 * np.pi * k / n), np.sin(2 * np.pi * k / n), 0.0]) * 10
        s.axis = np.array([0.0, 0.0, 1.0])
    br = [b for k in range(n - 1) for b in ladder(ring[k], ring[k + 1], "A", 6)] + ladder(ring[-1], ring[0], "A", 3)
    pend = []
    for x, at in extra:
        x.centroid, x.axis = ring[at].centroid * 1.0 + np.array([0.0, 0.0, 12.0]), np.array([0.0, 0.0, 1.0])
        br += ladder(x, ring[at], "A", 4)
        pend.append(x)
    return ring, pend, build_sheets(ring + pend, br)


def test_a_barrel_keeps_its_ring_and_draws_a_pendant_strand_off_it():
    """9FWZ's case in small: a 6-strand ring and a short strand pairing with one ring strand only. The barrel is the
    6 ring strands (closed, its ends pair); the pendant is an accessory strand beside its partner."""
    ring, (x,), (sh,) = _ring_sheet(6, [(fake(500, 504), 2)])
    assert sh.closed and len(sh.strands) == 6 and set(sh.strands) == set(ring)
    assert {sh.strands[0].start, sh.strands[-1].start} in ({0, 100}, {100, 0})  # cut at the weakest pairing
    ((acc, parent, column, _, _),) = sh.accessory
    assert acc is x and parent is ring[2] and abs(column - sh.strands.index(ring[2])) == 1  # beside its partner


def test_a_three_strand_cycle_is_no_barrel_and_draws_open_without_a_note():
    """Three strands pairing in a triangle, then a run of three more off one of them: a cycle of 3 is no barrel
    (_BARREL_MIN 5), so the sheet is drawn open, the run in a row; no note calls it a barrel drawn open (fix round 3,
    the review: the note named a barrel here); the note names the pairing it does not draw side by side."""
    from foldmap.layout import sheet_note_items

    def sheet_notes(sheets):
        return [i["text"] for i in sheet_note_items(sheets)]

    t = [fake(20 * k, 20 * k + 9) for k in range(6)]
    br = ladder(t[0], t[1], "A", 5) + ladder(t[1], t[2], "A", 5) + ladder(t[0], t[2], "A", 3)
    br += ladder(t[2], t[3], "A", 5) + ladder(t[3], t[4], "A", 5) + ladder(t[4], t[5], "A", 5)
    (sh,) = build_sheets(t, br)
    assert not sh.closed and len(sh.strands) == 6
    assert sheet_notes([sh]) == [f"pairings not drawn side by side: {t[0].id} also pairs with {t[2].id}"]


def test_a_five_ring_with_a_tail_on_a_path_tree_draws_closed():
    """8C3A sheets 103 and 220 (fix round 3, the review): the spanning tree is a path (tail t0, then t1 .. t5), and
    t1 and t5 pair (2 bridges), closing a 5-ring: it was never searched for (only branched trees were), so the sheet
    was drawn open with a false note. Every cyclic sheet takes its longest ring: a closed 5-ring and t0 beside it."""
    t = [fake(20 * k, 20 * k + 9) for k in range(6)]
    br = [b for p, q in zip(t, t[1:]) for b in ladder(p, q, "A", 5)] + ladder(t[1], t[5], "A", 2)
    (sh,) = build_sheets(t, br)
    assert sh.closed and {s.start for s in sh.strands} == {20, 40, 60, 80, 100}
    assert [(a.strand.start, a.parent.start) for a in sh.accessory] == [(0, 20)]


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
    drawn = [s.id for sh in sheets for s in [*sh.strands, *(a.strand for a in sh.accessory)]]
    assert len(drawn) == len(set(drawn)) == len(strands)  # every strand in exactly one sheet (track or accessory)
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
    assert ids in (["A:563-572", "A:499-509", "A:520-523"], ["A:520-523", "A:499-509", "A:563-572"])
    assert not sh.ambiguous  # deposited order is 563-572 | 499-509 | 520-523


def test_tim_barrel_is_closed_by_its_single_beta7_beta8_bridge():
    """1TIM: beta8 pairs back with beta1 (two bridges) and beta7 (one bridge): a closed barrel, as published."""
    from pathlib import Path

    from foldmap.dssp import assign_dssp
    from foldmap.io import load_backbone
    from foldmap.ss import build_sses

    bb = load_backbone(Path(__file__).parent / "data" / "1TIM.cif", "asu")
    d = assign_dssp(bb)
    sheets = build_sheets(build_sses(bb, d.ss), d.bridges)
    assert len(sheets) == 2 and all(sh.closed and len(sh.strands) == 8 for sh in sheets)


def test_single_bridges_do_not_close_ordinary_sheets():
    from pathlib import Path

    from foldmap.dssp import assign_dssp
    from foldmap.io import load_backbone
    from foldmap.ss import build_sses

    for name in ("1UBQ", "5NKT", "2JTY", "6ZS5"):
        bb = load_backbone(Path(__file__).parent / "data" / f"{name}.cif", "asu")
        d = assign_dssp(bb)
        assert not any(sh.closed for sh in build_sheets(build_sses(bb, d.ss), d.bridges)), name


def test_bulged_strands_take_their_true_place_in_the_sheet():
    from helpers import FID, pipeline

    _, _, sheets = pipeline(FID / "1TEN.cif.gz")
    sh = next(s for s in sheets if any(x.id == "A:868-877" for x in s.strands))
    known = ["A:847-852", "A:832-839", "A:868-877", "A:880-890"]  # C' C F G, as Lappalainen 2008 Fig 3
    assert [s.id for s in sh.strands] in (known, known[::-1]) and not sh.ambiguous
    _, _, sheets = pipeline(Path(__file__).parent / "data" / "5NKT.cif")
    sh = next(s for s in sheets if any(x.id == "A:68-74" for x in s.strands))
    known = ["A:29-35", "A:8-14", "A:149-158", "A:68-74", "A:114-123"]  # B'/B'' one strand, beside 68-74
    assert [s.id for s in sh.strands] in (known, known[::-1]) and not sh.ambiguous


def test_a_strand_continuing_another_along_their_partner_shares_its_column():
    """5KQJ (random-100 C34): F pairs with E end to end with B, on E's same side (F at 3D y +11.8, B at -1.1, 0.8 Å
    apart along E). The sheet's one path of strands drew F as a column of its own, wedged between E and B (C D E F B
    A); C34 stacked it in B's column, end to end with B, as it lies (0.6 Å apart). Plan 12 fix round 3 ruling (one
    column means collinear): F stands in B's column, end to end with B, level with E, its partner."""
    from foldmap import figcheck
    from foldmap.cli import make_layout
    from foldmap.spec import draw_spec, figure_spec

    path = Path(__file__).parent / "data" / "r100" / "5KQJ.cif.gz"
    lay, sses, _ = make_layout(path)
    b, e, f = lay.placed["A:44-50"], lay.placed["A:83-93"], lay.placed["A:96-99"]
    assert abs(b.cx - f.cx) < 1e-9 and abs(e.cx - b.cx) > 0.5
    span = lambda p: sorted((p.n_port[1], p.c_port[1]))  # noqa: E731
    (by0, by1), (ey0, ey1), (fy0, fy1) = span(b), span(e), span(f)
    assert by1 <= fy0 or fy1 <= by0  # end to end with B
    assert min(ey1, fy1) > max(ey0, fy0)  # level with E
    block = next(ids for ids in lay.sheets if "A:96-99" in ids)
    assert len({round(lay.placed[k].cx, 6) for k in block}) == len(block) - 1
    d = draw_spec(figure_spec(path))
    assert not [v for v in figcheck.order_checks(d) if v.startswith("sheet")]


def test_the_loops_out_of_a_shared_column_never_meet_in_its_gap():
    """5KQJ (batch C2a review Q1): B's loop left its top upward and F's its bottom downward on the same line in the
    1.2 gap between them, meeting at one junction, so the column read as F running into B. The two loops leaving the
    facing ends of a stack gap (Layout.stack_gaps) start nearer those ends and turn out of the gap, the lower as low
    as it can and the upper as high, so no point of one comes within 0.05 of the other there."""
    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    def in_box(points, box):  # the polyline sampled every 0.02, the samples inside the box
        out = []
        for p, q in zip(points, points[1:]):
            for t in np.linspace(0.0, 1.0, max(2, int(np.hypot(q[0] - p[0], q[1] - p[1]) / 0.02) + 1)):
                x, y = p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t
                if box[0] - 1e-9 <= x <= box[2] + 1e-9 and box[1] - 1e-9 <= y <= box[3] + 1e-9:
                    out.append((x, y))
        return np.array(out).reshape(-1, 2)

    lay, sses, bb = make_layout(Path(__file__).parent / "data" / "r100" / "5KQJ.cif.gz")
    loops = route_loops(lay, sses, bb)
    assert lay.stack_gaps and not any(lp.fallback for lp in loops)
    for box in lay.stack_gaps:  # the loops leaving the gap's two facing ends (fix round 1: F above E, its partner;
        cx = (box[0] + box[2]) / 2  # other loops may cross the gap)
        facing = [
            lp
            for lp in loops
            if any(
                abs(p[0] - cx) < 1e-6 and min(abs(p[1] - box[1]), abs(p[1] - box[3])) < 1e-6
                for p in (lp.points[0], lp.points[-1])
            )
        ]
        assert len(facing) == 2
        a, b = (in_box(lp.points, box) for lp in facing)
        if len(a) and len(b):
            d = a[:, None, :] - b[None, :, :]
            assert np.hypot(d[..., 0], d[..., 1]).min() >= 0.05


def test_partners_far_apart_along_a_long_strand_share_its_edge_in_one_column():
    """C34's controls (batch C2a review Q3), revised by plan 12: 8UTF's a:266-273 and a:318-322 both pair with
    a:278-301 on one side, 27.6 Å apart along it. C2a kept each in a column of its own, so one stood between
    a:278-301 and the other, which do not pair (an S1 failure); the one off the track now stands in the other's
    column, end to end with it (collinear along a:278-301's edge, fix round 3 ruling). A sheet with no branch (1TIM's
    barrel) has no accessory strand and no stack gap."""
    from foldmap.cli import make_layout

    data = Path(__file__).parent / "data"
    bb = load_backbone(data / "8UTF.cif", "auto")
    dssp = assign_dssp(bb)
    sheets = build_sheets(build_sses(bb, dssp.ss, bridges=dssp.bridges), dssp.bridges)
    pair = {"a:266-273", "a:318-322"}
    (sh,) = [sh for sh in sheets if any(a.strand.id in pair for a in sh.accessory)]
    (acc,) = [a for a in sh.accessory if a.strand.id in pair]
    assert acc.parent.id == "a:278-301" and sh.strands[acc.column].id in pair  # end to end with the other
    lay, _, _ = make_layout(data / "1TIM.cif")
    assert lay.stack_gaps == [] and lay.accessory == {}
    bb = load_backbone(data / "1TIM.cif", "auto")
    dssp = assign_dssp(bb)
    assert all(sh.accessory == [] for sh in build_sheets(build_sses(bb, dssp.ss, bridges=dssp.bridges), dssp.bridges))


def test_fimd_barrel_has_24_strands_and_its_pendant_strand_beside_its_partner():
    """9FWZ FimD (the review's R1): the barrel is 24 strands, D:645-659 next to D:139-153 (12 bridges); D:662-666,
    which pairs with D:645-659 only, is an accessory strand beside it."""
    from foldmap.science import sheet_checks

    path = Path(__file__).parent / "data" / "r100" / "9FWZ-DB.cif.gz"
    bb = load_backbone(path, "asu")
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss, short_helices=True, bridges=d.bridges)
    sheets = build_sheets(sses, d.bridges)
    (barrel,) = [sh for sh in sheets if sh.closed]
    ids = [s.id for s in barrel.strands]
    assert len(ids) == 24 and abs(ids.index("D:645-659") - ids.index("D:139-153")) == 1
    assert [(a.strand.id, a.parent.id) for a in barrel.accessory] == [("D:662-666", "D:645-659")]
    assert ids[barrel.accessory[0].column] == "D:634-641"  # beside D:645-659, on D:634-641's side (the review)
    assert sheet_checks(sheets, sses, d.bridges) == dict.fromkeys(("S1", "S2", "S3", "S6"), []) | {
        "S4": None,
        "S5": None,
    }


def test_ompf_barrel_has_16_strands_its_n_terminal_strand_beside_a9_23():
    """OmpF 2OMF (plan 12 pre-scan): A:2-6 pairs with A:9-23 only and was drawn between A:9-23 and A:36-50, a
    17-strand barrel; the porin's barrel has 16 strands."""
    path = Path(__file__).parent / "data" / "r100" / "2OMF-A.cif.gz"
    bb = load_backbone(path, "asu")
    d = assign_dssp(bb)
    (barrel,) = [sh for sh in build_sheets(build_sses(bb, d.ss, bridges=d.bridges), d.bridges) if sh.closed]
    ids = [s.id for s in barrel.strands]
    assert len(ids) == 16 and abs(ids.index("A:9-23") - ids.index("A:36-50")) == 1
    assert [(a.strand.id, a.parent.id) for a in barrel.accessory] == [("A:2-6", "A:9-23")]


def test_a_rings_chords_are_named():
    """Fix round 3 (the review): 9P1X's 9-ring has chords (O:29-36/P:28-35) its drawing cannot show side by side;
    Sheet.chords holds them, and the notes and summary name them ("also pairs with")."""
    from foldmap.layout import sheet_note_items

    def sheet_notes(sheets):
        return [i["text"] for i in sheet_note_items(sheets)]

    ring, _, (sh,) = _ring_sheet(6)
    extra = ladder(ring[1], ring[4], "A", 4)
    ring2 = build_sheets(
        ring,
        [b for k in range(5) for b in ladder(ring[k], ring[k + 1], "A", 6)] + ladder(ring[-1], ring[0], "A", 3) + extra,
    )
    (sh,) = ring2
    assert [(a.start, b.start) for a, b in sh.chords] in ([(20, 80)], [(80, 20)])
    assert any("pairings not drawn side by side" in n for n in sheet_notes([sh]))


def test_notes_file_each_sheets_lines_under_their_own_heading():
    """Fix round 4 (the re-review's N1, 9P1X --focus none): sheet_notes joined every sheet's accessory lines and
    chord lines into one list and cut it at the total accessory count, so with a chord in one sheet and an accessory
    in a later one the chord was listed among the strands off the run and the accessory under "pairings not drawn
    side by side". Each line now stays under its own heading."""
    from foldmap.layout import sheet_note_items

    def sheet_notes(sheets):
        return [i["text"] for i in sheet_note_items(sheets)]

    t = [fake(20 * k, 20 * k + 9) for k in range(6)]  # sheet 1: a 3-cycle, then a run; its chord t0/t2
    br = ladder(t[0], t[1], "A", 5) + ladder(t[1], t[2], "A", 5) + ladder(t[0], t[2], "A", 3)
    br += ladder(t[2], t[3], "A", 5) + ladder(t[3], t[4], "A", 5) + ladder(t[4], t[5], "A", 5)
    c, a, b, d = fake(200, 209), fake(220, 229), fake(240, 249), fake(260, 269)  # sheet 2: d off the run, at c
    br += ladder(c, a, "A") + ladder(c, b, "A") + ladder(c, d, "A")
    one, two = build_sheets([*t, c, a, b, d], br)
    assert one.chords and not one.accessory and two.accessory and not two.chords
    assert sheet_notes([one, two]) == [
        f"1 strand(s) off a sheet's run of paired strands: {d.id} drawn beside its partner {c.id}",
        f"pairings not drawn side by side: {t[0].id} also pairs with {t[2].id}",
    ]


def test_a_sheet_drawn_open_holds_no_ring_the_builder_finds():
    """Fix round 4 (the re-review, section 3): `Sheet.unclosed` and its note ("barrel closure not supported by
    pairing: drawn open") could no longer fire. _ring is the tree path between the two ends of a pairing left out of
    the maximum spanning tree, so whenever it finds a ring of _BARREL_MIN strands or more the track's ends pair and
    the sheet is closed; and a found ring always becomes the track (it is never longer than the tree's longest path,
    and when as long it covers it). Both are gone; this pins why over random pairing networks: no sheet drawn open
    holds a ring _ring finds. (A ring that is no single tree cycle, the sum of two, is not found: science S3 flags such
    a sheet drawn open, open_finding; none in the scan set or the sample.)"""
    import random

    import networkx as nx

    import foldmap.layout as layout
    from foldmap.sheets import _edge_table, _ring

    rng = random.Random(12)
    opened = 0
    for _ in range(400):
        n = rng.randint(5, 10)
        t = [fake(20 * k, 20 * k + 9) for k in range(n)]
        for s in t:
            s.centroid, s.axis = np.array([rng.uniform(-10, 10) for _ in range(3)]), np.array([0.0, 0.0, 1.0])
        br = [b for a, c in zip(t, t[1:]) for b in ladder(a, c, "A", rng.randint(2, 6))]
        for _ in range(rng.randint(1, 4)):
            i, j = rng.sample(range(n), 2)
            br += ladder(t[i], t[j], "A", rng.randint(1, 6))
        graph = nx.Graph()
        graph.add_weighted_edges_from((a, b, sum(c.values())) for (a, b), c in _edge_table(t, br).items())
        for sh in build_sheets(t, br):
            assert not hasattr(sh, "unclosed")
            if sh.closed:
                continue
            opened += 1
            sub = graph.subgraph([t.index(x) for x in [*sh.strands, *(a.strand for a in sh.accessory)]])
            assert _ring(nx.maximum_spanning_tree(sub), sub) is None
    assert opened > 50 and not hasattr(layout, "UNCLOSED")


@pytest.mark.parametrize("weights, drawn_closed", [((4,) * 6, True), ((2, 4, 2, 2, 2, 2), False)])
def test_a_two_cycle_pairing_network_is_flagged_by_s3_when_drawn_open(weights, drawn_closed):
    """M1 (plan 12 T1+T2 re-review): pairings 0-20, 20-40, 40-60, 60-80, 0-40, 20-80 make two cycles sharing the
    pairing 20-40 (0-20-40 and 20-40-60-80), whose sum is a five-strand ring (0-20-80-60-40) that is no single tree
    cycle. `_ring` finds it only when the spanning tree happens to lie along it: with equal pairings it does (drawn
    closed, S3 passes), with 20-40 stronger it misses it (drawn open, a gap recorded for the sheet builder) and S3
    flags the ring the open sheet holds. Over all 729 weightings of three levels, no sheet is open and unflagged."""
    from foldmap.science import sheet_checks

    s = [fake(20 * k, 20 * k + 9) for k in range(5)]
    for k, x in enumerate(s):
        x.centroid, x.axis = np.array([np.cos(1.2 * k), np.sin(1.2 * k), 0.0]) * 10, np.array([0.0, 0.0, 1.0])
    pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 2), (1, 4)]
    br = [b for (i, j), w in zip(pairs, weights) for b in ladder(s[i], s[j], "A", w)]
    sheets = build_sheets(s, br)
    found = sheet_checks(sheets, s, br)
    assert any(sh.closed for sh in sheets) is drawn_closed
    assert bool(found["S3"]) is (not drawn_closed)
