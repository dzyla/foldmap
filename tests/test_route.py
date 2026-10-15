import logging
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from foldmap import route
from foldmap.frame import view_frame
from foldmap.layout import ETA_W, HELIX_W, PITCH, Layout, Placed, build_layout, label_boxes
from foldmap.route import _INFLATE, _MAX_CELLS, _Grid, route_loops
from helpers import fake_sse, pipeline

DATA = Path(__file__).parent / "data"

MODES = ["projected", "stack"]


def strand(start, cx, cy, length=3.0, angle=np.pi / 2):
    return Placed(fake_sse("E", start, start + 7), cx, cy, length, angle, PITCH, "A")


def bar(start, cx, cy, length, angle):
    return Placed(fake_sse("H", start, start + 9), cx, cy, length, angle, HELIX_W, "a")


def layout_of(*placed):
    rects = np.array([p.rect for p in placed])
    b = (rects[:, 0].min(), rects[:, 1].min(), rects[:, 2].max(), rects[:, 3].max())
    return Layout({p.sse.id: p for p in placed}, tuple(float(x) for x in b))


def stub_bb(n=200):
    return SimpleNamespace(prev=np.ones(n, bool))


def inside(pt, rect, eps=0.05):
    return rect[0] + eps < pt[0] < rect[2] - eps and rect[1] + eps < pt[1] < rect[3] - eps


def segments_clear(points, rects, step=0.05):
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        n = max(2, int(np.hypot(x1 - x0, y1 - y0) / step))
        for t in np.linspace(0, 1, n):
            pt = (x0 + t * (x1 - x0), y0 + t * (y1 - y0))
            if any(inside(pt, r) for r in rects):
                return False
    return True


def axis_aligned(points, tol=1e-9):
    return all(abs(a[0] - b[0]) < tol or abs(a[1] - b[1]) < tol for a, b in zip(points, points[1:]))


def test_two_strands_are_joined_by_an_orthogonal_path_around_the_obstacles():
    a = strand(0, 0, 0)  # points up; C port on top
    b = strand(20, 4, 0, angle=-np.pi / 2)  # points down; N port on top
    lay = layout_of(a, b)
    (loop,) = route_loops(lay, [a.sse, b.sse], stub_bb())
    assert not loop.fallback and not loop.dashed
    assert np.allclose(loop.points[0], a.c_port) and np.allclose(loop.points[-1], b.n_port)
    assert axis_aligned(loop.points)
    assert segments_clear(loop.points, [a.rect, b.rect])
    assert max(p[1] for p in loop.points) > a.c_port[1] + 0.2  # goes over the tops


def test_walled_in_port_falls_back_and_finishes():
    a = strand(0, 0, 0)
    b = strand(20, 4, 0, length=2.0, angle=-np.pi / 2)
    walls = [
        bar(40, 4, 1 + HELIX_W / 2 + 0.2, 5, 0),
        bar(50, 4, -1 - HELIX_W / 2 - 0.2, 5, 0),
        bar(60, 1.5 - HELIX_W / 2, 0, 5, np.pi / 2),
        bar(70, 6.5 + HELIX_W / 2, 0, 5, np.pi / 2),
    ]
    lay = layout_of(a, b, *walls)
    (loop,) = route_loops(lay, [a.sse, b.sse], stub_bb())
    assert loop.fallback
    assert np.allclose(loop.points[0], a.c_port) and np.allclose(loop.points[-1], b.n_port)
    assert len(loop.points) > 2


def shared_length(p, q):
    total = 0.0
    for a0, a1 in zip(p, p[1:]):
        for b0, b1 in zip(q, q[1:]):
            if abs(a0[0] - a1[0]) < 1e-9 and abs(b0[0] - b1[0]) < 1e-9 and abs(a0[0] - b0[0]) < 0.05:
                total += max(0, min(max(a0[1], a1[1]), max(b0[1], b1[1])) - max(min(a0[1], a1[1]), min(b0[1], b1[1])))
            if abs(a0[1] - a1[1]) < 1e-9 and abs(b0[1] - b1[1]) < 1e-9 and abs(a0[1] - b0[1]) < 0.05:
                total += max(0, min(max(a0[0], a1[0]), max(b0[0], b1[0])) - max(min(a0[0], a1[0]), min(b0[0], b1[0])))
    return total


def test_later_loops_do_not_run_on_top_of_earlier_ones():
    s = [
        strand(0, 0, 0),
        strand(20, 6, 0, angle=-np.pi / 2),
        strand(40, 2, 0),
        strand(60, 4, 0, angle=-np.pi / 2),
    ]
    lay = layout_of(*s)
    loops = route_loops(lay, [p.sse for p in s], stub_bb())
    assert len(loops) == 3 and not any(l.fallback for l in loops)
    for i, l1 in enumerate(loops):
        for l2 in loops[i + 1 :]:
            assert shared_length(l1.points, l2.points) < 0.3, (l1.a_id, l2.a_id)


def test_dashed_exactly_across_chain_breaks(ubq, tmp_path):
    import gemmi

    st = gemmi.read_structure(str(ubq))
    chain = st[0]["A"]
    idx = next(i for i, r in enumerate(chain) if r.seqid.num == 20)
    del chain[idx]  # residue 20 is in the loop between strand 12-16 and helix 23-34
    gapped = tmp_path / "gap.cif"
    st.make_mmcif_document().write_file(str(gapped))
    bb, sses, sheets = pipeline(gapped)
    lay = build_layout(sses, sheets, view_frame(sses), "projected")
    loops = route_loops(lay, sses, bb)
    dashed = [(l.a_id, l.b_id) for l in loops if l.dashed]
    assert dashed == [("A:12-16", "A:23-34")]
    assert len(loops) == len(sses) - 1


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("name", ["ubq", "zs5", "zya"])
def test_real_structures_loops_join_their_ports_and_stay_clear(name, mode, request):
    bb, sses, sheets = pipeline(request.getfixturevalue(name))
    lay = build_layout(sses, sheets, view_frame(sses), mode)
    loops = route_loops(lay, sses, bb)
    pairs = sum(1 for a, b in zip(sses, sses[1:]) if a.chain == b.chain)
    assert len(loops) == pairs
    for l in loops:
        a, b = lay.placed[l.a_id], lay.placed[l.b_id]
        assert np.allclose(l.points[0], a.c_port) and np.allclose(l.points[-1], b.n_port)
        if l.fallback:
            continue
        others = [p.rect for k, p in lay.placed.items() if k not in (l.a_id, l.b_id)]
        assert segments_clear(l.points, others, step=0.1), (l.a_id, l.b_id)  # never through a third element
        assert segments_clear(l.points[1:-1], [a.rect, b.rect], step=0.1), (
            l.a_id,
            l.b_id,
        )  # only the stubs touch the ends
    assert sum(l.fallback for l in loops) <= max(1, len(loops) // 4)


def test_deterministic(zs5):
    bb, sses, sheets = pipeline(zs5)
    lay = build_layout(sses, sheets, view_frame(sses), "projected")
    a = [(l.a_id, l.points) for l in route_loops(lay, sses, bb)]
    b = [(l.a_id, l.points) for l in route_loops(lay, sses, bb)]
    assert a == b


def test_loops_go_around_barrel_ghost_strands():
    bb, sses, sheets = pipeline(__import__("pathlib").Path(__file__).parent / "data" / "1EMA.cif")
    lay = build_layout(sses, sheets, view_frame(sses), "projected")
    (ghost,) = lay.ghosts
    for loop in route_loops(lay, sses, bb):
        if not loop.fallback:
            assert segments_clear(loop.points, [ghost.rect]), (loop.a_id, loop.b_id)


@pytest.mark.parametrize("pdb", ["1TIM", "2HHB"])
def test_routes_never_circle_back_over_themselves(pdb):
    from pathlib import Path

    bb, sses, sheets = pipeline(Path(__file__).parent / "data" / f"{pdb}.cif")
    lay = build_layout(sses, sheets, view_frame(sses), "projected")
    for loop in route_loops(lay, sses, bb):
        if loop.fallback:
            continue
        assert not self_crossing(loop.points), (loop.a_id, loop.b_id, loop.points)  # a curl crosses itself


def self_crossing(points, eps=1e-9):
    """True if two non-adjacent segments of the polyline touch or cross."""

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def on(p, a, b):
        return (
            min(a[0], b[0]) - eps <= p[0] <= max(a[0], b[0]) + eps
            and min(a[1], b[1]) - eps <= p[1] <= max(a[1], b[1]) + eps
        )

    def meet(a, b, c, d):
        d1, d2, d3, d4 = cross(c, d, a), cross(c, d, b), cross(a, b, c), cross(a, b, d)
        if ((d1 > eps and d2 < -eps) or (d1 < -eps and d2 > eps)) and (
            (d3 > eps and d4 < -eps) or (d3 < -eps and d4 > eps)
        ):
            return True
        return any(
            abs(v) <= eps and on(p, *seg)
            for v, p, seg in ((d1, a, (c, d)), (d2, b, (c, d)), (d3, c, (a, b)), (d4, d, (a, b)))
        )

    segs = list(zip(points, points[1:]))
    return any(meet(*segs[i], *segs[j]) for i in range(len(segs)) for j in range(i + 2, len(segs)))


def test_kinked_helix_pieces_are_joined_straight_not_by_a_fallback():
    from pathlib import Path

    from foldmap.cli import make_layout

    lay, sses, bb = make_layout(Path(__file__).parent / "data" / "8UUP.cif")
    loops = {(l.a_id, l.b_id): l for l in route_loops(lay, sses, bb)}
    for key in [("B:456-472", "B:474-478"), ("B:474-478", "B:480-485"), ("F:456-472", "F:474-477")]:
        l = loops[key]
        assert not l.fallback and len(l.points) == 2, key  # one straight link across the kink
        assert np.allclose(l.points[0], lay.placed[key[0]].c_port) and np.allclose(
            l.points[1], lay.placed[key[1]].n_port
        )


def clip_len(p, q, rect):
    """Length of segment p-q strictly inside rect (Liang-Barsky clipping); 0.0 when it misses."""
    x0, y0, x1, y1 = rect
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for pk, qk in ((-dx, p[0] - x0), (dx, x1 - p[0]), (-dy, p[1] - y0), (dy, y1 - p[1])):
        if pk == 0:
            if qk <= 0:
                return 0.0
            continue
        if pk < 0:
            t0 = max(t0, qk / pk)
        else:
            t1 = min(t1, qk / pk)
        if t0 >= t1:
            return 0.0
    return (t1 - t0) * np.hypot(dx, dy)


def segment_clear(p, q, rects, grow=_INFLATE) -> bool:
    """True when segment p-q keeps the clearance from every rect (each grown by `grow`)."""
    return all(clip_len(p, q, (r[0] - grow, r[1] - grow, r[2] + grow, r[3] + grow)) <= 1e-9 for r in rects)


@pytest.mark.parametrize(
    "third",
    [
        Placed(fake_sse("G", 0, 3, chain="B"), 0.7, 0.0, 2.0, np.pi / 2, HELIX_W * ETA_W, "η1"),
        Placed(fake_sse("E", 0, 6, chain="B"), 0.7, 0.0, 3.0, np.pi / 2, PITCH, "A"),
    ],
    ids=["3-10", "strand"],
)
def test_kink_shortcut_does_not_cross_a_third_element(third):
    a = bar(0, -2.0, 0.0, 4.0, 0.0)  # C port at (0, 0), leaving along +x
    b = bar(11, 3.4, 0.0, 4.0, 0.0)  # N port at (1.4, 0): one residue on, close enough for a kink link
    lay = layout_of(a, b, third)
    (loop,) = route_loops(lay, [a.sse, b.sse, third.sse], stub_bb())
    grown = np.add(third.rect, (-_INFLATE, -_INFLATE, _INFLATE, _INFLATE))  # with the router's clearance
    clear = all(clip_len(p, q, grown) < 1e-9 for p, q in zip(loop.points, loop.points[1:]))
    assert loop.fallback or clear  # a drawn route (not the marked fallback curve) never passes through it


def corner_chains():
    """Four lone helices far apart: they stretch the page until the router's grid cells are coarse."""
    return [
        Placed(fake_sse("H", 0, 9, chain=c), x, y, 4.0, 0.0, HELIX_W, "α")
        for c, x, y in (("W", 0, 0), ("X", 1480, 0), ("Y", 0, 1480), ("Z", 1480, 1480))
    ]


def other_rects(loop, placed):
    return [p.rect for p in placed if p.sse.id not in (loop.a_id, loop.b_id)]


@pytest.mark.parametrize("step", range(12))
def test_coarse_grid_keeps_narrow_elements(step):
    a = Placed(fake_sse("E", 0, 7), 700.0, 700.0, 3.0, np.pi / 2, PITCH, "A")  # up; C port on top
    b = Placed(fake_sse("E", 20, 27), 706.0, 700.0, 3.0, -np.pi / 2, PITCH, "B")  # down; N port on top
    eta = Placed(fake_sse("G", 0, 3, chain="C"), 702.4 + step * 0.1, 702.0, 16.0, np.pi / 2, HELIX_W * ETA_W, "η1")
    placed = [a, b, eta, *corner_chains()]
    lay = layout_of(*placed)
    assert _Grid(lay).res > 1.2  # narrower than a cell once grown by the clearance
    loops = route_loops(lay, [p.sse for p in placed], stub_bb(10_000))
    (loop,) = [l for l in loops if l.a_id == a.sse.id]
    assert not loop.fallback and all(
        segment_clear(p, q, [eta.rect], grow=0.0) for p, q in zip(loop.points, loop.points[1:])
    )


@pytest.mark.parametrize("page", ["normal", "big"])
def test_random_layouts_never_route_through_an_element(page):
    rng = np.random.default_rng(0)
    widths = {"H": HELIX_W, "E": PITCH, "G": HELIX_W * ETA_W}
    checked = fallbacks = 0
    for _ in range(150):
        off = 700.0 if page == "big" else 0.0
        placed = []
        for i in range(int(rng.integers(4, 9))):
            kind = str(rng.choice(["H", "E", "G"]))
            ang = float(rng.choice([0, np.pi / 2, np.pi, -np.pi / 2]))
            x, y, length = rng.uniform(0, 14) + off, rng.uniform(0, 10) + off, rng.uniform(2, 6)
            placed.append(Placed(fake_sse(kind, 20 * i, 20 * i + 9), x, y, length, ang, widths[kind], f"e{i}"))
        r = np.array([p.rect for p in placed])
        if any(
            min(r[i, 2], r[j, 2]) > max(r[i, 0], r[j, 0]) - 0.05
            and min(r[i, 3], r[j, 3]) > max(r[i, 1], r[j, 1]) - 0.05
            for i in range(len(placed))
            for j in range(i + 1, len(placed))
        ):
            continue  # only non-overlapping boxes; near-touching ones are allowed
        placed += corner_chains() if page == "big" else []
        lay = layout_of(*placed)
        assert (_Grid(lay).res > 1.0) == (page == "big")
        for loop in route_loops(lay, [p.sse for p in placed], stub_bb(10_000)):
            checked += 1
            fallbacks += loop.fallback
            if not loop.fallback:
                rects = other_rects(loop, placed)
                assert all(segment_clear(p, q, rects, grow=0.0) for p, q in zip(loop.points, loop.points[1:])), loop
    assert checked > 100
    # 5 and 35 of 134 now: a route can no longer be a silent crossing, so over-blocking shows only as fallbacks
    # (a port skip that ignored the coarse reach made it 100 on the big page)
    assert fallbacks <= {"normal": 8, "big": 45}[page], fallbacks


def test_5lj3_window_routes_without_fallback():
    """A long gap or free end inside a large RNP routes clear (K5). 5LJ3's RNA V with Z's stem and the protein
    elements beside them on the full figure's page (tests/data/r100/5LJ3-V): the 3' end V:85-92 was drawn back over
    its own stem, as a free end's clearance check trimmed its far end too when that lay inside the stem, and the gap
    after V:92 then had no way out of the stem ("gap after V:92", as on the full figure). Every free end now keeps
    out of its own stem past where it leaves it, and every path is drawn clear."""
    from foldmap.na_paths import in_body
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(DATA / "r100" / "5LJ3-V.cif.gz"))
    lay = d.layout
    at = {key: s for s in lay.stems for key in s.axial}
    for p in lay.na_paths:
        if p.kind == "tail":
            (key,) = p.between
            assert not any(in_body(q, at[key]) for q in p.nt_xy), p.nts
    assert any(p.kind == "gap" for p in lay.na_paths)
    assert d.stuck == []


def logged(caplog, what):
    """The count route_loops logged for `what` in its debug line."""
    (n,) = re.findall(rf"\b(\d+) {what}", caplog.text)
    return int(n)


def strands_over_the_top(*more):
    """Two strands joined over their tops (as in the first test), plus any further elements."""
    a, b = strand(0, 0, 0), strand(20, 4, 0, angle=-np.pi / 2)
    return a, b, layout_of(a, b, *more)


def test_a_route_pushed_into_an_element_after_a_star_is_redrawn_unsnapped(monkeypatch, caplog):
    """A post-processing step that moves a found route into its own elements: the route is checked again
    after it and redrawn from the A* cells without the port snap, never drawn through the elements."""
    snap = route._snap_to_port

    def snap_too_far(pts, port, exit_vec, at_start):
        snap(pts, port, exit_vec, at_start)
        for p in pts:
            p[1] -= 2.0  # down through the bodies of both strands

    monkeypatch.setattr(route, "_snap_to_port", snap_too_far)
    caplog.set_level(logging.DEBUG, logger="foldmap.route")
    a, b, lay = strands_over_the_top()
    (loop,) = route_loops(lay, [a.sse, b.sse], stub_bb())
    assert logged(caplog, "drawings rejected as crossing") == 1  # the snapped drawing
    assert not loop.fallback
    assert np.allclose(loop.points[0], a.c_port) and np.allclose(loop.points[-1], b.n_port)
    assert all(segment_clear(p, q, [a.rect, b.rect], grow=0.0) for p, q in zip(loop.points, loop.points[1:]))
    assert max(p[1] for p in loop.points) > a.c_port[1] + 0.2  # still over the tops


def test_a_route_through_an_element_that_cannot_be_redrawn_is_marked_fallback(monkeypatch, caplog):
    """A* cells straight through a third element (a broken router step): neither the snapped nor the
    unsnapped drawing is clear, so the loop is the flagged fallback curve, not a silent crossing."""
    wall = bar(40, 2, 1.0, 6.0, np.pi / 2)  # stands between the strands, up past their tops
    a, b, lay = strands_over_the_top(wall)
    monkeypatch.setattr(
        route,
        "_astar",
        lambda grid, start, sd, goal, gd, *budget: [(x, start[1]) for x in range(start[0], goal[0] + 1)],
    )
    caplog.set_level(logging.DEBUG, logger="foldmap.route")
    (loop,) = route_loops(lay, [a.sse, b.sse], stub_bb())
    assert loop.fallback
    assert logged(caplog, "fallback curves") == 1
    assert logged(caplog, "drawings rejected as crossing") == 4  # snapped and not, with and without labels


# res 0.5 and 1.3 only: below 2 * _INFLATE = 0.3 the inflated rect alone always spans a cell centre, so a
# finer grid could not tell whether the reach is there
@pytest.mark.parametrize("res", [0.5, 1.3])
@pytest.mark.parametrize("width", [0.0, 0.35])
def test_coarse_grid_blocks_the_cells_under_any_narrow_obstacle(res, width):
    side = res * np.sqrt(_MAX_CELLS)  # a page this size gets cells of this size
    lay = Layout({}, (0.0, 0.0, side - 2 * 3.0, side - 2 * 3.0))
    for k in range(8):  # slide the obstacle across one cell in both directions
        x = y = 100.0 + k * res / 8
        g = _Grid(lay, [(x, y, x + width, y + width)])
        assert g.res == pytest.approx(res)
        (i0, j0), (i1, j1) = g.cell((x, y)), g.cell((x + width, y + width))
        assert g.soft[i0 : i1 + 1, :].any() and g.soft[:, j0 : j1 + 1].any(), k  # its own column and row bands
        assert g.soft[i0 : i1 + 1, j0 : j1 + 1].any(), k  # and a cell it lies in


def kink_with(third):
    """Helix a leaves (0, 0) along +x; helix b, one residue on, is entered at (1, 1) from below: a kink."""
    a = bar(0, -2.0, 0.0, 4.0, 0.0)
    b = bar(11, 1.0, 3.0, 4.0, np.pi / 2)
    lay = layout_of(a, b, third)
    (loop,) = route_loops(lay, [a.sse, b.sse], stub_bb())  # chain B has no chain-end stubs in the way
    return a, b, lay, loop


def test_kink_blocked_by_a_small_element_is_routed_around_it():
    third = Placed(fake_sse("E", 0, 6, chain="B"), 0.5, 0.5, 0.3, np.pi / 4, 0.3, "A")  # on the diagonal
    a, b, _, loop = kink_with(third)
    assert not segment_clear(a.c_port, b.n_port, [third.rect], grow=0.0)
    assert not loop.fallback and not loop.crosses_labels and len(loop.points) > 2
    assert all(segment_clear(p, q, [third.rect]) for p, q in zip(loop.points, loop.points[1:]))


def test_debug_summary_counts_kinks(monkeypatch, caplog):
    """The routing debug line (read by the random-100 campaign, plan 08b K1) counts the kink links rejected for a
    third element in the way and the fallback curves that cross an element, beside the fallback curves."""
    third = Placed(fake_sse("E", 0, 6, chain="B"), 0.5, 0.5, 0.3, np.pi / 4, 0.3, "A")  # on the kink's diagonal
    caplog.set_level(logging.DEBUG, logger="foldmap.route")
    kink_with(third)
    assert logged(caplog, "kink links rejected") == 1
    assert logged(caplog, "fallback curves") == 0 and logged(caplog, "of them crossing an element") == 0
    caplog.clear()
    monkeypatch.setattr(route, "_astar", lambda *args: None)  # no grid route: the loop is the fallback curve
    _, _, _, loop = kink_with(third)
    assert loop.fallback and all(segment_clear(p, q, [third.rect]) for p, q in zip(loop.points, loop.points[1:]))
    assert logged(caplog, "fallback curves") == 1 and logged(caplog, "of them crossing an element") == 0
    caplog.clear()
    a, b = bar(0, -2.0, 0.0, 4.0, 0.0), bar(11, 1.0, 3.0, 4.0, np.pi / 2)  # as kink_with draws them
    mid = loop.points[len(loop.points) // 2]  # and a fourth element where that curve passes
    fourth = Placed(fake_sse("E", 20, 26, chain="B"), *mid, 0.3, 0.0, 0.3, "A")
    (loop,) = route_loops(layout_of(a, b, third, fourth), [a.sse, b.sse], stub_bb())
    assert loop.fallback and not all(segment_clear(p, q, [fourth.rect]) for p, q in zip(loop.points, loop.points[1:]))
    assert logged(caplog, "fallback curves") == 1 and logged(caplog, "of them crossing an element") == 1


def test_kink_blocked_only_by_another_elements_label_is_routed_around_it():
    third = Placed(fake_sse("H", 0, 9, chain="B"), 0.15, 0.7, 0.4, -np.pi / 2, 0.1, "α")  # label to its right
    a, b, lay, loop = kink_with(third)
    label = dict(label_boxes(lay, [a.sse, b.sse]))["label:B:0-9"]
    assert segment_clear(a.c_port, b.n_port, [third.rect])  # the element itself leaves room for the link
    assert not segment_clear(a.c_port, b.n_port, [label])  # its label does not
    assert not loop.fallback and not loop.crosses_labels and len(loop.points) > 2
    assert all(segment_clear(p, q, [label, third.rect]) for p, q in zip(loop.points, loop.points[1:]))


def test_kink_between_tilted_pieces_with_overlapping_boxes_stays_straight():
    a = bar(0, -1.3, -0.75, 3.0, np.pi / 6)  # C port near (0, 0)
    b = bar(11, a.c_port[0] + 0.6 + 1.5 * np.cos(np.pi / 3), a.c_port[1] + 1.5 * np.sin(np.pi / 3), 3.0, np.pi / 3)
    assert not segment_clear(a.c_port, b.n_port, [a.rect], grow=0.0)  # the link runs inside both boxes
    assert not segment_clear(a.c_port, b.n_port, [b.rect], grow=0.0)
    lay = layout_of(a, b)
    (loop,) = route_loops(lay, [a.sse, b.sse], stub_bb())
    assert not loop.fallback and len(loop.points) == 2


@pytest.mark.parametrize("x", [700.0, 700.3, 700.6, 700.9])
def test_coarse_grid_still_routes_an_open_loop(x):
    """On a coarse grid a port's first cells are blocked by its own element out to the cell reach; they must
    be skipped like the element's own box, not taken to wall the port in."""
    a = Placed(fake_sse("E", 0, 7), x, 700.0, 3.0, np.pi / 2, PITCH, "A")  # up; C port on top
    b = Placed(fake_sse("E", 20, 27), x + 6.0, 700.0, 3.0, -np.pi / 2, PITCH, "B")  # down; N port on top
    placed = [a, b, *corner_chains()]
    lay = layout_of(*placed)
    assert _Grid(lay).res > 1.2
    (loop,) = [l for l in route_loops(lay, [p.sse for p in placed], stub_bb(10_000)) if l.a_id == a.sse.id]
    assert not loop.fallback
    assert all(segment_clear(p, q, [a.rect, b.rect], grow=0.0) for p, q in zip(loop.points, loop.points[1:]))


def helix_piece(start, port, angle, length):
    """A helix piece whose N port is at `port`."""
    c = np.add(port, np.array([np.cos(angle), np.sin(angle)]) * length / 2)
    return Placed(fake_sse("H", start, start + 9), float(c[0]), float(c[1]), length, angle, HELIX_W, "α")


def test_kink_is_not_drawn_through_the_box_of_a_third_piece_of_the_same_helix():
    """Only the two joined pieces are exempt. A third piece, beyond the next kink, is an obstacle by its
    bounding box like any other element (the grid blocks the same box), even where its tilted body leaves
    room: the link is not drawn straight, and the loop is routed or flagged instead."""
    a = helix_piece(0, (-3 * np.cos(np.pi / 6), -3 * np.sin(np.pi / 6)), np.pi / 6, 3.0)  # C port at (0, 0)
    b = helix_piece(11, (0.5, 0.29), np.pi / 6, 0.6)  # a short piece one residue on
    c = helix_piece(22, np.add(b.c_port, (0.08, 0.1)), 3 * np.pi / 4, 4.0)  # turns back up over the a-b gap
    assert not segment_clear(a.c_port, b.n_port, [c.rect])  # within the clearance of c's box
    assert not route._crosses_body(a.c_port, b.n_port, c)  # though not through c itself
    loops = {l.a_id: l for l in route_loops(layout_of(a, b, c), [a.sse, b.sse, c.sse], stub_bb())}
    assert len(loops[a.sse.id].points) > 2  # not the straight kink link
    assert loops[b.sse.id].points == [b.c_port, c.n_port]  # the b-c kink is


@pytest.mark.parametrize("assembly", ["asu", "auto"])
def test_loops_between_tm_helices_stay_out_of_the_slab(assembly):
    """1C3W ASU: the C-D loop ran along y = -4.04 inside the band -4.40..2.96."""
    from foldmap.cli import make_layout

    lay, sses, bb = make_layout(DATA / "1C3W.cif", membrane="auto", assembly=assembly)
    y0, y1 = lay.membrane["y"]
    for lp in route_loops(lay, sses, bb):
        for (xa, ya), (xb, yb) in zip(lp.points, lp.points[1:]):
            if abs(ya - yb) < 1e-9 and abs(xa - xb) > 0.5:
                assert not y0 + 0.3 < ya < y1 - 0.3, (lp.a_id, lp.b_id, ya)
        assert not lp.fallback


@pytest.mark.parametrize("assembly", ["asu", "auto"])
def test_stacked_loops_stay_out_of_the_slab_except_at_a_port_inside_it(assembly):
    """1C3W --mode stack --membrane auto: a loop runs sideways inside a band only next to one of its own ports
    that is itself inside the band (a strand or helix ending within the bilayer): on the stub out of that port or
    the first run after it (alpha2 -> strand A passes under the hairpin, whose lower port is in the band)."""
    from foldmap.cli import make_layout

    lay, sses, bb = make_layout(DATA / "1C3W.cif", mode="stack", membrane="auto", assembly=assembly)
    bands = lay.membrane["bands"]
    assert len(bands) == len({s.chain for s in sses})  # each protomer wider than _MAX_ROW: a row each

    def inside(y):
        return any(lo + 0.3 < y < hi - 0.3 for lo, hi in bands)

    for lp in route_loops(lay, sses, bb):
        n = len(lp.points)
        for k, ((xa, ya), (xb, yb)) in enumerate(zip(lp.points, lp.points[1:])):
            if abs(ya - yb) < 1e-9 and abs(xa - xb) > 0.5 and inside(ya):
                at_port = (k <= 1 and inside(lp.points[0][1])) or (k >= n - 3 and inside(lp.points[-1][1]))
                assert at_port, (lp.a_id, lp.b_id, ya)
        assert not lp.fallback


def test_the_router_goes_over_the_slab_rather_than_along_it():
    """Two ports below a thin band, an element between them reaching up to the band's lower face: the shortest
    way runs sideways just inside the slab; with the slab cost the path climbs over the band instead."""
    lay = Layout({}, (0.0, 0.0, 20.0, 12.0))
    lay.membrane = {"y": (5.0, 6.5), "inside": "below"}
    grid = _Grid(lay, [(6.0, -5.0, 14.0, 5.0)])  # the element between the two ports, up to the lower face
    path = route._astar(grid, grid.cell((2.0, 4.5)), 0, grid.cell((18.0, 4.5)), 0)
    assert path is not None
    sideways = [a for a, b in zip(path, path[1:]) if a[1] == b[1] and grid.slab[a] and grid.slab[b]]
    assert not sideways, [grid.center(c) for c in sideways[:3]]
    assert max(grid.center(c)[1] for c in path) > 6.5  # over the band


def test_only_the_over_helix_port_walks_loose(monkeypatch):
    """A loop between a helix drawn over a sheet and another element: only the helix's own port may walk out over
    blocked cells (loose); the other port leaves as any port does (Task 14 review)."""
    from foldmap.cli import make_layout
    from foldmap.style import Style

    calls = []
    walk = _Grid.free_cell_along

    def spy(self, port, exit_vec, own, loose=False):
        calls.append((own.sse.id, own.over, loose))
        return walk(self, port, exit_vec, own, loose)

    monkeypatch.setattr(_Grid, "free_cell_along", spy)
    lay, sses, bb = make_layout(DATA / "1UBQ.cif", look=Style(helix_on_sheet="auto"))
    route_loops(lay, sses, bb)
    assert any(over for _, over, _ in calls)
    assert all(loose == over for _, over, loose in calls), [c for c in calls if c[1] != c[2]]


@pytest.mark.parametrize("path", ["1UBQ.cif", "fidelity/2ACY.cif.gz"])
def test_a_loop_leaves_a_helix_over_a_sheet_straight_across_only_the_strand_under_its_port(path):
    """From a helix drawn over a sheet, a loop leaves its port in the page direction nearest the helix's own and runs
    straight out past the sheet's edge, crossing no strand but the one its port lies on (Task 14 review)."""
    from foldmap.cli import make_layout
    from foldmap.route import _crosses_body
    from foldmap.style import Style

    lay, sses, bb = make_layout(DATA / path, look=Style(helix_on_sheet="auto"))
    loops = route_loops(lay, sses, bb)
    checked = 0
    for lp in loops:
        assert not lp.fallback, (lp.a_id, lp.b_id)
        for helix, (p, q), port in (
            (lay.placed[lp.a_id], lp.points[:2], lp.points[0]),
            (lay.placed[lp.b_id], lp.points[-2:][::-1], lp.points[-1]),
        ):
            if not helix.over:
                continue
            checked += 1
            assert abs(p[0] - q[0]) < 1e-9 or abs(p[1] - q[1]) < 1e-9, (lp.a_id, lp.b_id, p, q)  # straight out
            for k in lay.on_sheet[helix.sse.id]:
                x0, y0, x1, y1 = lay.placed[k].rect
                if _crosses_body(p, q, lay.placed[k]):
                    assert x0 <= port[0] <= x1 and y0 <= port[1] <= y1, (lp.a_id, lp.b_id, k)
    assert checked >= 2


@pytest.mark.parametrize(
    "path, kw", [("1UBQ.cif", {}), ("fidelity/2ACY.cif.gz", {}), ("2LZM.cif", {}), ("na/4OO8.cif.gz", {"dna": False})]
)
def test_a_helix_over_a_sheet_is_an_obstacle_to_every_other_loop(path, kw):
    """A helix drawn over a sheet is still in the way of loops that neither join it nor a strand of its sheet: they
    go round its box, never under it (Task 14 review)."""
    from foldmap.cli import make_layout
    from foldmap.style import Style

    lay, sses, bb = make_layout(DATA / path, look=Style(helix_on_sheet="auto"), **kw)
    assert lay.on_sheet
    for lp in route_loops(lay, sses, bb):
        if lp.fallback:
            continue
        for h, ids in lay.on_sheet.items():
            if {lp.a_id, lp.b_id} & {h, *ids}:
                continue
            assert segments_clear(lp.points, [lay.placed[h].rect]), (lp.a_id, lp.b_id, h)


def test_a_helix_over_a_sheet_keeps_its_box_as_a_soft_obstacle():
    """Its box is an obstacle owned by the helix (route._element_obstacles); in the grid the cells it covers beyond
    the strands are soft, not free: given up to a loop only when nothing else routes (Task 14 review)."""
    from foldmap.cli import make_layout
    from foldmap.route import _element_obstacles
    from foldmap.style import Style

    lay, _, _ = make_layout(DATA / "1UBQ.cif", look=Style(helix_on_sheet="auto"))
    rects, owners = _element_obstacles(lay)
    assert "A:23-34" in owners
    assert tuple(rects[:, owners.index("A:23-34")]) == lay.placed["A:23-34"].rect
    grid = _Grid(lay)
    x0, y0, x1, y1 = lay.placed["A:23-34"].rect
    inside = np.outer((grid.xs > x0) & (grid.xs < x1), (grid.ys > y0) & (grid.ys < y1))
    hard = _Grid(Layout({k: p for k, p in lay.placed.items() if not p.over}, lay.bounds)).blocked
    assert (inside & ~hard).any() and np.all(grid.soft[inside & ~hard]) and np.all(grid.blocked[inside])


def test_a_loop_out_of_a_helix_over_a_sheet_keeps_off_text_unless_flagged():
    """2ACY α2's N port lies on β1 just above β1's N-terminal letter: the walk out of the sheet stops at its edge, so
    the loop into it crosses that text only as a last resort, flagged (Loop.crosses_labels), never silently."""
    from foldmap.cli import make_layout
    from foldmap.style import Style
    from helpers import FID

    lay, sses, bb = make_layout(FID / "2ACY.cif.gz", look=Style(helix_on_sheet="auto"))
    boxes = [b for _, b in label_boxes(lay, sses)]
    for lp in route_loops(lay, sses, bb):
        if not lp.fallback and not lp.crosses_labels:
            assert segments_clear(lp.points, boxes), (lp.a_id, lp.b_id)


def _name_in_the_way(bounds=(-4.0, -3.0, 4.0, 3.0)):
    """A horizontal helix whose name sits above it, and a loop line that runs just above the name's ink: outside the
    name's keep-out box (so the router would not mind it), inside the room figcheck keeps between a line and text."""
    from foldmap.figcheck import LINE_TEXT_CLEARANCE_PX
    from foldmap.layout import page_text

    h = bar(10, 0.0, 0.0, 3.0, 0.0)
    lay = Layout({h.sse.id: h}, bounds)
    em = route._ink_em(lay, lay.font_scale)
    room = (LINE_TEXT_CLEARANCE_PX + 0.5) * 0.72 * page_text(lay)[1]
    top = route._name_ink(lay, h, em)[3]
    y = top + room / 2
    loop = route.Loop("x", "y", [(-2.0, y), (2.0, y)], False, False)
    return lay, h, loop, y, room


def test_a_helix_name_a_loop_runs_too_near_moves_to_the_clear_side():
    """random-100 C6: the route keeps out of a name's rough box (smaller than its text on a wide page) and its own
    helix's loops may leave through it, so the name is set clear of every loop's line after routing."""
    lay, h, loop, y, room = _name_in_the_way()
    assert h.label_side == 1 and label_boxes(lay, [h.sse])[0][1][3] < y  # outside the keep-out box
    route._clear_helix_names(lay, [h.sse], [loop])
    ink = route._name_ink(lay, h, route._ink_em(lay, lay.font_scale))
    assert h.label_side == -1 and ink[3] < y - room  # below the helix, with the room a line keeps from text


def test_a_helix_name_no_loop_runs_near_stays():
    lay, h, loop, y, room = _name_in_the_way()
    loop.points[:] = [(-2.0, y + 3 * room), (2.0, y + 3 * room)]
    route._clear_helix_names(lay, [h.sse], [loop])
    assert (h.label_side, h.label_along) == (1, 0.0)


def test_a_helix_name_with_no_better_spot_stays():
    lay, h, loop, y, room = _name_in_the_way(bounds=(-4.0, -0.5, 4.0, 3.0))  # no room below the helix
    other = [(-2.0, y - 0.9), (2.0, y - 0.9)]
    route._clear_helix_names(lay, [h.sse], [loop, route.Loop("x", "z", other, False, False)])
    assert h.label_side == 1
