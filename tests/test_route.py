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
from foldmap.spec import draw_spec, figure_spec
from helpers import fake_sse, pipeline, r100_entry

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


def test_a_loop_leaves_the_slab_rather_than_run_along_it_over_an_element():
    """C26 (random-100, ruling: loops inside the slab are real defects). An element reaching 0.3 into the band
    between two ports below it (7SMR B:171-173>B:190-201 ran over such an end at 0.4 inside the bottom face): the
    slab cost per cell made the run along the slab over it cheaper than the way under it; a step along the slab
    now costs _SLAB_RUN, so the route goes under. Across the slab (TM loops between faces) costs as before."""
    lay = Layout({}, (0.0, -18.0, 20.0, 14.0))
    lay.membrane = {"y": (5.0, 12.0), "inside": "below"}
    grid = _Grid(lay, [(6.0, -15.0, 14.0, 5.3)])  # the element between the ports, its top inside the band
    path = route._astar(grid, grid.cell((2.0, 4.5)), 0, grid.cell((18.0, 4.5)), 0)
    assert path is not None
    sideways = [a for a, b in zip(path, path[1:]) if a[1] == b[1] and grid.slab[a] and grid.slab[b]]
    assert not sideways, [grid.center(c) for c in sideways[:3]]
    assert min(grid.center(c)[1] for c in path) < -15.0  # under the element


def test_a_loop_between_ports_below_the_band_stays_out_of_the_slab():
    """random-100 held-out 7JI3 (A18: C:101-108>C:116-119): a loop from a strand's foot to the next strand, both ports
    below the band, ran up into the slab and over the head of a strand reaching into it, the way round the crowded
    sheet below costing more than the run along the slab. Two ports below the band and a long element between them
    whose head is 0.6 inside it: ports on one side of every band are joined with the slab closed first, so the loop
    goes the long way under the element and nothing runs in the slab."""
    from foldmap.figcheck import slab_runs

    a = bar(0, 0.0, 1.0, 4.0, np.pi / 2)  # C port at (0, 3), below the band, pointing up
    b = bar(10, 10.0, 1.0, 4.0, -np.pi / 2)  # N port at (10, 3)
    wall = bar(40, 5.0, -27.2, 65.6, np.pi / 2)  # from y = -60 up to 5.6, 0.6 inside the band's lower face
    lay = layout_of(a, b, wall)
    band = (5.0, 12.0)
    lay.membrane = {"y": band, "bands": [band], "inside": "below"}
    loops = route_loops(lay, [a.sse, b.sse, wall.sse], stub_bb())
    (ab,) = [lp for lp in loops if (lp.a_id, lp.b_id) == (a.sse.id, b.sse.id)]
    assert not ab.fallback
    assert slab_runs([ab], [band], lay.placed) == []
    assert max(y for _, y in ab.points) < band[0]


def _slab_case(case: str, dy: float):
    """The two slab-port geometries of random-100 C2b-5, shifted up by dy against a fixed grid origin (an element
    far below holds the page's bounds): `outside`, 6JO6's helix ending 0.03 below the lower face, its C port up into
    the slab, the next element's port inside the band; `inband`, 5LDW's short helix in the band, its C port up into
    a 0.81 gap under the helix stacked on it, with a tall helix 0.81 to its right whose loop enters that gap."""
    anchor = bar(90, -20.0, -10.0, 1.0, np.pi / 2)
    if case == "outside":
        els = [bar(0, 0.0, dy - 1.03, 2.0, np.pi / 2), bar(20, 3.0, dy + 3.0, 5.0, np.pi / 2)]
        band = (dy, dy + 8.0)
    else:
        els = [
            bar(0, 1.81, dy + 1.2, 6.4, -np.pi / 2),  # right of the corridor, its C port at the bottom
            bar(10, 0.0, dy + 2.71, 3.8, np.pi / 2),  # stacked above, its N port at 0.81
            bar(20, 0.0, dy - 1.035, 2.07, np.pi / 2),  # C port at (0, 0), up into the gap
            bar(30, 5.28, dy - 1.71, 1.72, np.pi / 2),  # N port below the band
        ]
        band = (dy - 2.10, dy + 4.80)
    lay = layout_of(*els, anchor)
    lay.membrane = {"y": band, "bands": [band], "inside": "below"}
    return lay, [p.sse for p in [*els, anchor]], band


@pytest.mark.parametrize("case", ["outside", "inband"])
def test_a_route_out_of_a_port_whose_first_step_lands_in_the_slab_turns_at_the_port(case):
    """Random-100 C2b-5 (6JO6 5:167-169>5:179-207; 5LDW L:448-454 in C17's second layout). Every route started 0.4
    out along its exit (_PUSH): out of a port at a face, pointing in, the turn round its element's end lay 0.42 inside
    the face; out of an in-band port facing a narrow gap, the 0.4 row was the one the loop entering the gap took,
    and the route jogged to run mid-slab (A18). Fix round 1 (review I1): whatever the grid's phase (the elements
    shifted across one cell in 0.025 steps), a port whose first step lands in the slab starts in the first free row
    past its clearance (fifth-of-a-cell steps from _INFLATE), and a port just outside a face has its first turn
    drawn at its clearance when that is clear: no loop runs in the slab, none falls back. A port whose first step
    leaves the slab, or that leaves sideways, starts as before."""
    from foldmap.figcheck import slab_runs

    for k in range(10):
        lay, sses, band = _slab_case(case, k * 0.025)
        grid = _Grid(lay, [])
        a = lay.placed[sses[0 if case == "outside" else 2].id]
        assert grid.slab[grid.cell((a.c_port[0], a.c_port[1] + 0.4))]  # the usual first step is in the slab
        start = grid.center(grid.free_cell_along(a.c_port, a.exit_c, a))[1] - a.c_port[1]
        assert _INFLATE < start <= _INFLATE + grid.res + 1e-9, (k, start)  # the first free row past the clearance
        loops = route_loops(lay, sses, stub_bb())
        assert not any(lp.fallback for lp in loops), k
        assert slab_runs(loops, [band], lay.placed) == [], (k, [lp.points for lp in loops])
    if case == "inband":
        return
    c = bar(40, 0.0, 9.0, 2.0, np.pi / 2)  # a port above the band, its first step out of it: as before
    up = _Grid(layout_of(c), [])
    port, out = c.c_port, c.exit_c
    assert up.free_cell_along(port, out, c) == up.cell((port[0] + out[0] * 0.4, port[1] + out[1] * 0.4))
    # a port leaving sideways inside the band (6JLN's c:81-83): its first step is in the slab, it starts as before
    d = bar(60, 4.0, 4.0, 1.175, 0.0)
    side = layout_of(d)
    side.membrane = {"y": (0.0, 8.0), "bands": [(0.0, 8.0)], "inside": "below"}
    grid = _Grid(side, [])
    port, out = d.c_port, d.exit_c
    assert abs((port[0] - grid.x0) % grid.res - 0.175) < 1e-9  # a nearer start would find a free cell
    assert grid.slab[grid.cell((port[0] + 0.4, port[1]))]
    assert grid.free_cell_along(port, out, d) == grid.cell((port[0] + out[0] * 0.4, port[1] + out[1] * 0.4))


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


def _name_in_the_way(bounds=(-8.0, -6.0, 8.0, 6.0)):
    """A horizontal helix whose name sits above it, and a line that runs just above the name's ink, inside the room
    figcheck keeps between a line and text: a fallback curve (no route holds it: random-100 C6, routed loops keep
    clear of names sized from the drawn text, test_routed_loops_keep_off_every_helix_name_as_drawn)."""
    from foldmap.figcheck import LINE_TEXT_CLEARANCE_PX
    from foldmap.layout import page_text

    h = bar(10, 0.0, 0.0, 3.0, 0.0)
    h.label = "α10"
    lay = Layout({h.sse.id: h}, bounds)
    lay.font_scale = 1.8  # larger text: a larger name (and keep-out box)
    em = route._ink_em(lay, lay.font_scale)
    room = (LINE_TEXT_CLEARANCE_PX + 0.5) * 0.72 * page_text(lay)[1]
    top = route._name_ink(lay, h, em)[3]
    y = top + room / 2
    loop = route.Loop("x", "y", [(-2.0, y), (2.0, y)], False, True)
    return lay, h, loop, y, room


def test_a_helix_name_a_loop_runs_too_near_moves_to_the_clear_side():
    """random-100 C6: a line no route holds (a fallback curve) that runs through a name or within figcheck's room of
    its ink moves the name clear of it, after routing (the mover is the fallback of the router's keep-outs)."""
    lay, h, loop, y, room = _name_in_the_way()
    ink = route._name_ink(lay, h, route._ink_em(lay, lay.font_scale))
    assert h.label_side == 1 and ink[3] < y < ink[3] + room  # within the room a line keeps from text
    route._clear_helix_names(lay, [h.sse], [loop])
    ink = route._name_ink(lay, h, route._ink_em(lay, lay.font_scale))
    assert h.label_side == -1 and ink[3] < y - room  # below the helix, with the room a line keeps from text


def test_a_helix_name_no_loop_runs_near_stays():
    lay, h, loop, y, room = _name_in_the_way()
    loop.points[:] = [(-2.0, y + 3 * room), (2.0, y + 3 * room)]
    route._clear_helix_names(lay, [h.sse], [loop])
    assert (h.label_side, h.label_along) == (1, 0.0)


def test_a_helix_name_with_no_better_spot_stays():
    lay, h, loop, y, room = _name_in_the_way(bounds=(-8.0, -0.5, 8.0, 6.0))  # no room below the helix
    other = [(-2.0, y - 0.9), (2.0, y - 0.9)]
    route._clear_helix_names(lay, [h.sse], [loop, route.Loop("x", "z", other, False, False)])
    assert h.label_side == 1


def test_a_helix_name_is_not_moved_onto_a_line_to_keep_further_from_others():
    """Batch C2b-1 (random-100 6YF9 focus: α2 of BN, DM, EB, EK, EN under a loop of another chain): on a page of
    small text the room a line keeps from a name is wider than the name, so every spot had a line within it; the
    mover counted those lines alike, near or through, and took the spot with the fewest: under the helix, on the
    loop. A line through the name now weighs more than any number of lines near it."""
    lay, h, loop, y, room = _name_in_the_way()
    em = route._ink_em(lay, lay.font_scale)
    top = route._name_ink(lay, h, em)[3]
    near = [
        route.Loop("x", f"n{k}", [(-6.0, top + f * room), (6.0, top + f * room)], False, True)
        for k, f in enumerate((0.75, 0.95))
    ]  # two lines near the name, above it, all along the helix
    h.label_side = -1
    below = route._name_ink(lay, h, em)
    h.label_side = 1
    through = route.Loop("x", "t", [(-6.0, (below[1] + below[3]) / 2), (6.0, (below[1] + below[3]) / 2)], False, True)
    route._clear_helix_names(lay, [h.sse], [*near, through])  # one line through every spot below it
    assert h.label_side == 1


def test_a_helix_name_keeps_off_a_nucleotide_trace():
    """Fix round 1: nucleic-acid traces count where a name may go (a name the mover put on a strand would undo
    layout._element_label_sides, 1EYG)."""
    lay, h, loop, y, room = _name_in_the_way()
    loop.points[:] = [(-2.0, 5.0), (2.0, 5.0)]  # a line far away: only the trace strikes the name
    lay.na_traces = {0: np.array([[-0.2, 1.0], [0.2, 1.0]])}  # nucleotides where the name stands above the helix
    route._clear_helix_names(lay, [h.sse], [loop])
    assert h.label_side == -1


def test_a_mirror_copys_name_follows_its_mate(caplog):
    """Fix round 1: copies drawn as mirror images are named alike (_mirror_sides): when the mate's name moves to
    the other side, its mirror copy's goes with it (and the campaign's count says two names moved)."""
    caplog.set_level(logging.DEBUG, logger="foldmap.route")
    lay, mate, loop, y, room = _name_in_the_way()
    copy = bar(30, 5.0, 0.0, 3.0, np.pi)
    copy.label = mate.label
    lay.placed[copy.sse.id] = copy
    lay.mirror_pairs = [(copy.sse.id, mate.sse.id)]
    route._clear_helix_names(lay, [mate.sse, copy.sse], [loop])
    assert mate.label_side == -1 and copy.label_side == -1
    assert "2 helix names moved off lines" in caplog.text


def test_no_helix_name_of_1eyg_sits_on_a_nucleotide_trace():
    """Fix round 1: 1EYG's single strand runs over chain A's α1 and layout._element_label_sides put the name on the
    other side; the post-route mover must not undo that."""
    from foldmap.layout import _element_label_box, _trace_hits
    from foldmap.spec import draw_spec, figure_spec

    lay = draw_spec(figure_spec(DATA / "na" / "1EYG.cif.gz")).layout
    traces = [np.asarray(p, float).reshape(-1, 2) for _, p in sorted(lay.na_traces.items())]
    assert traces
    assert [k for k, p in lay.placed.items() if p.sse.kind in "HG" and _trace_hits(_element_label_box(p), traces)] == []


def test_protein_loops_keep_off_the_nucleotide_letters():
    """C11 (random-100): protein loops are routed before the nucleic acid's single strands, whose letters and 5'/3'
    labels are placed after, so the RRM loops of 6SO9 ran across three letters and the 3' label of its bound
    6-mer. A loop crossing NA text is routed again around it, once that text is placed (none dropped, none moved)."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(DATA / "r100" / "6SO9.cif.gz"))
    assert len(d.layout.na_letters) == 6 and not any(lp.fallback for lp in d.loops)
    assert d.over_labels == []  # and no route over a name instead (the campaign's A4: a reroute flagged over)
    hits = [h for h in figcheck.text_collisions(d) if h.startswith("na:")]
    assert hits == [], hits


def test_a_reroute_off_na_text_crosses_no_more_strands():
    """Fix round 1 (minor): C11's reroute took a route that crossed less NA text and no more text in all, whatever it
    crossed of the nucleic acid's strands, routed round the old path. A reroute may cross the strands at most once
    more for each piece of text it frees (6SO9's loop: 3 letters freed, 3 more crossings of the bound RNA's trace);
    here one letter freed for two more crossings: not taken."""
    from foldmap.na_paths import _Strands

    letter = [(4.0, 1.5, 5.0, 2.5)]
    strands = _Strands(None)
    strands.add([(6.0, -5.0), (6.0, 5.0)])  # a strand across x = 6
    old = [(0.0, 0.0), (0.0, 2.0), (10.0, 2.0), (10.0, 0.0)]  # through the letter, across the strand once
    around = [(0.0, 0.0), (0.0, 3.0), (10.0, 3.0), (10.0, 0.0)]  # over the letter, across the strand once
    weave = [(0.0, 0.0), (0.0, 3.0), (7.0, 3.0), (7.0, -1.0), (5.5, -1.0), (5.5, -2.0), (10.0, -2.0), (10.0, 0.0)]
    assert route._better(around, old, letter, [], strands)
    assert route._strand_crossings(weave, strands) == 3 and not route._better(weave, old, letter, [], strands)


def test_clip_lengths_is_clip_len_for_every_pair():
    """Fix round 1 (minor): the four segment-against-box loops (text hits, letter spots, tethers over bodies, tie
    spots) share route.clip_lengths, which gives clip_len for every (segment, box) pair; boxes_meeting is the open
    box overlap."""
    rng = np.random.default_rng(7)
    segs = rng.uniform(-3, 3, (40, 4))
    segs[:5, 2] = segs[:5, 0]  # vertical
    segs[5:10, 3] = segs[5:10, 1]  # horizontal
    lo = rng.uniform(-3, 2, (30, 2))
    boxes = np.hstack([lo, lo + rng.uniform(0.1, 2.0, (30, 2))])
    got = route.clip_lengths(segs, boxes)
    want = np.array([[route.clip_len(s[:2], s[2:], b) for b in boxes] for s in segs])
    assert np.allclose(got, want) and ((got > 0) == (want > 0)).all()
    assert route.boxes_meeting(boxes, (0.0, 0.0, 1.0, 1.0)).tolist() == [
        b[0] < 1 and b[2] > 0 and b[1] < 1 and b[3] > 0 for b in boxes
    ]


def test_a_reroute_may_turn_more_to_free_text_it_ran_through():
    """Fix round 1 campaign: with the clearance in pixels, a reroute of a loop that only grazes a letter's clearance
    may add at most _REROUTE_TURNS corners (1KX5's loop wove round letters it did not touch), but one that runs
    through the letters' ink may add up to _REROUTE_TURNS_INK (7XSP's B:355-374 loop: 6 more to free two letters)."""
    ink = [(4.0, 1.8, 5.0, 2.2)]
    pad = [(3.8, 1.6, 5.2, 2.4)]
    through = [(0.0, 0.0), (0.0, 2.0), (10.0, 2.0), (10.0, 0.0)]  # through the ink
    grazing = [(0.0, 0.0), (0.0, 2.3), (10.0, 2.3), (10.0, 0.0)]  # inside the clearance pad only
    wiggle = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 4.0), (6.0, 4.0), (6.0, 3.0)]
    wiggle += [(8.0, 3.0), (8.0, 1.0), (10.0, 1.0), (10.0, 0.0)]  # 8 corners more, clear of both
    assert route._better(wiggle, through, pad, [], None, ink)
    assert not route._better(wiggle, grazing, pad, [], None, ink)


def _names_struck_by_routes(path, scale: str) -> list[tuple[str, str]]:
    """(helix, loop) for every helix name whose drawn ink (route._name_ink) a routed loop (not a fallback curve) comes
    within figcheck's line clearance of, as the loops stand when the post-route mover (_clear_helix_names) starts."""
    from foldmap.layout import _ink_em, line_clearance
    from foldmap.spec import draw_spec, figure_spec

    got = []
    mover = route._clear_helix_names

    def look(layout, sses, loops):
        em, grow = _ink_em(layout, layout.font_scale), line_clearance(layout)
        for p in layout.placed.values():
            if p.sse.kind not in ("H", "G") or not p.label:
                continue
            x0, y0, x1, y1 = route._name_ink(layout, p, em)
            room = (x0 - grow, y0 - grow, x1 + grow, y1 + grow)
            for lp in loops:
                if not lp.fallback and any(route.clip_len(a, b, room) > 1e-9 for a, b in zip(lp.points, lp.points[1:])):
                    got.append((p.sse.id, f"{lp.a_id}>{lp.b_id}"))
        mover(layout, sses, loops)

    route._clear_helix_names = look
    try:
        draw_spec(figure_spec(path, style=[f"font_scale={scale}"]))
    finally:
        route._clear_helix_names = mover
    return got


@pytest.mark.parametrize(("name", "scale"), [("8UUP", "1"), ("8UUP", "2"), ("8UTF", "2"), ("25IN", "2")])
def test_routed_loops_keep_off_every_helix_name_as_drawn(name, scale):
    """random-100 C6, router side: the boxes loops route round are sized from the text as drawn (C4: text held at
    6 pt, its em the same page units on any page; layout.glyph grows the boxes with larger text), so no routed loop
    comes within figcheck's line clearance of a helix name's ink, its own loops' ends included, and the post-route
    mover has nothing to move: it is left for lines no route holds (fallback curves, nucleic-acid paths). At f57887a
    (4 pt text on 8UUP's 84-unit page; boxes of a fixed size under larger text): 8UUP 2 names struck (one by its own
    loop), 45 at text size 2; 8UTF 10 and 25IN 2 at text size 2."""
    assert _names_struck_by_routes(DATA / f"{name}.cif", scale) == []


def test_keep_out_boxes_grow_with_larger_text():
    """random-100 C6: a name's and a chain end's keep-out box is sized from the text drawn (layout.glyph): three
    times the text, three times the glyph (the 0.05 pad aside); never smaller than at the usual size."""
    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(DATA / "2HHB.cif")
    small = dict(label_boxes(lay, sses))
    lay.font_scale = 3.0
    big = dict(label_boxes(lay, sses))
    lay.font_scale = 0.6
    less = dict(label_boxes(lay, sses))
    for key in ("N:A", next(k for k in small if k.startswith("label:"))):
        (a0, b0, a1, b1), (c0, d0, c1, d1) = small[key], big[key]
        assert d1 - d0 == pytest.approx(3 * (b1 - b0)) and c1 - c0 == pytest.approx(3 * (a1 - a0 - 0.1) + 0.1)
        assert less[key] == small[key]


def test_a_stack_gap_toll_is_never_negative(monkeypatch):
    """Fix round 4: the toll _apart_in_stack_gaps lays on a gap's cells (rising with the height, or the depth) took
    cells whose centres lie just past the gap's bottom or top, where it went below zero (8UTF: -4.9). A* then found
    loops through them cheaper than a turn: a:428-431>a:440-442 curled back over itself (test_matrix). The toll is 0
    or more, and no loop of 8UTF crosses itself."""
    from foldmap.cli import make_layout

    seen = []
    real = route._route_again

    def spy(layout, router, lp):
        seen.append(float(router.grid.toll.min()))
        return real(layout, router, lp)

    monkeypatch.setattr(route, "_route_again", spy)
    lay, sses, bb = make_layout(DATA / "8UTF.cif", "projected")
    loops = route_loops(lay, sses, bb)
    assert seen and min(seen) >= 0.0, seen
    assert not [(l.a_id, l.b_id) for l in loops if not l.fallback and self_crossing(l.points)]


def _walled(n_wall=60.0):
    """Two strands either side of a long wall (a helix), whose loop must go round the wall's end."""
    a = strand(0, 0, 0)  # C port on top
    b = strand(20, 4, 0, angle=-np.pi / 2)  # N port on top
    wall = bar(40, 2, n_wall / 2 - 3, n_wall, np.pi / 2)
    return a, b, wall


def test_a_search_out_of_budget_is_run_again_wider_when_a_route_exists(monkeypatch):
    """Random-100 C2b-2 (7JU4, 7QUP, 8QEM): on a large page a loop round a wall of long elements needs more A* nodes
    than _MAX_NODES; it was drawn as a fallback curve through them. When free cells join its two ends the search runs
    once more with _MAX_NODES_WIDE (here _MAX_NODES is cut so this small page needs the second search)."""
    a, b, wall = _walled()
    monkeypatch.setattr(route, "_MAX_NODES", 500)
    (loop,) = route_loops(layout_of(a, b, wall), [a.sse, b.sse], stub_bb())
    assert not loop.fallback
    assert segments_clear(loop.points, [wall.rect])
    assert min(p[1] for p in loop.points) < wall.rect[1] or max(p[1] for p in loop.points) > wall.rect[3]  # round
    monkeypatch.setattr(route, "_MAX_NODES_WIDE", 0)  # no second search: the first runs out
    (loop,) = route_loops(layout_of(a, b, wall), [a.sse, b.sse], stub_bb())
    assert loop.fallback


def test_no_second_search_when_no_free_cells_join_the_ends(monkeypatch):
    """A port shut in a box of elements has no route, however far the search goes: the second search is not run
    (its time is saved), and the loop stays a reported fallback. A guard: the first search does run."""
    a = strand(0, 0, 0)
    b = strand(20, 4, 0, length=2.0, angle=-np.pi / 2)  # N port on top, room above it inside the box
    walls = [bar(40, 4, 2.6, 7, 0), bar(50, 4, -2.6, 7, 0)]
    walls += [bar(60, 1.2, 0, 6, np.pi / 2), bar(70, 6.8, 0, 6, np.pi / 2)]
    searches = []
    real = route._astar

    def counted(*args):
        searches.append(args[5])  # its budget
        return real(*args)

    monkeypatch.setattr(route, "_MAX_NODES", 5)  # the first search runs out of budget, not of cells
    monkeypatch.setattr(route, "_astar", counted)
    (loop,) = route_loops(layout_of(a, b, *walls), [a.sse, b.sse], stub_bb())
    assert loop.fallback
    assert searches and route._MAX_NODES_WIDE not in searches


@pytest.mark.parametrize(
    "entry, view",
    [("8QEM", {"focus": "none"}), ("7JU4", {})],
    ids=["8QEM-focus", "7JU4-default"],
)
def test_large_pages_route_their_long_loops(entry, view):
    """Random-100 C2b-2 A3/A4: 8QEM's focus view (A:221-239>A:248-255 and its copy in B) and 7JU4's default view
    (four J and K loops round its wall of long helices, four loops over labels) ran out of the A* budget."""
    d = draw_spec(figure_spec(r100_entry(entry), **view))
    assert d.stuck == [] and d.over_labels == []


def test_a_route_found_with_the_keep_outs_open_is_over_a_label_only_if_it_crosses_one():
    """Random-100 C2b-2 A4 (8RRX's D and G loops, 9P1X Q:98-103): a label's soft cells beside a port's way out end
    its first, strict search (_Grid.free_cell_along), so the route is found with the keep-outs open; drawn, it
    passes beside the label, in no keep-out box. It is no loop over a label: the flag is the drawn route's own test,
    as _route_again's (2ACY's loop through its label, test_cli, stays flagged)."""
    a = strand(0, 0, 0)  # C port on top at (0, 1.5), leaving up
    b = strand(20, 4, 0, angle=-np.pi / 2)
    lay = layout_of(a, b)
    grid = _Grid(lay)
    cx, _ = grid.center(grid.cell((a.c_port[0], a.c_port[1] + route._PUSH[0])))  # the first cell out of a's port
    # a label beside the port's axis line, near enough that its keep-out margin covers that cell's centre
    box = (cx + 0.1, 1.7, cx + 1.0, 2.3) if cx > 0.1 else (cx - 1.0, 1.7, cx - 0.1, 2.3)
    router = route._Router(_Grid(lay, [box]), route._element_obstacles(lay), [box])
    assert router.grid.free_cell_along(a.c_port, a.exit_c, a) is None  # the strict search cannot leave the port
    pts, over = router.route(a.c_port, a.exit_c, a, b.n_port, b.exit_n, b)
    assert pts is not None and route._text_hits(pts, [box]) == 0
    assert not over


@pytest.mark.parametrize("view", [{"domains": "auto"}, {"focus": "none"}], ids=["domains", "focus"])
def test_8rrx_loops_beside_labels_are_not_over_them(view):
    """Random-100 C2b-2 A4: 8RRX's D:4985-4997>D:5000-5002 (domains) and G:4158-4166>G:4168-4174 (focus)."""
    d = draw_spec(figure_spec(r100_entry("8RRX"), **view))
    assert d.stuck == [] and d.over_labels == []


def _two_walls():
    """Two pairs of strands, each either side of its own long wall: two loops that each need the second search."""
    a, b, wall = _walled()
    c = strand(100, 20, 0)
    d = strand(120, 24, 0, angle=-np.pi / 2)
    wall2 = bar(140, 22, 27, 60, np.pi / 2)
    return [a, b, c, d], [wall, wall2]


def test_the_second_searches_share_one_budget_per_page(monkeypatch):
    """Fix round 1 (review: the wide search had no per-page cap; 7QUP's focus view took 158 s): every second
    search of one figure's routing draws on _MAX_NODES_PAGE nodes, counted as A* expands them (never time). Once
    it is spent a search that runs out falls back as before, a reported fallback curve; the order is routing
    order, so the same page always spends it alike."""
    strands, walls = _two_walls()
    pairs = [strands[0].sse, strands[1].sse], [strands[2].sse, strands[3].sse]
    spends = []
    real = route._astar

    def counted(*args):
        got = real(*args)
        if len(args) > 6:  # a second search: what it spent
            spends.append(args[6][0])
        return got

    def draw():
        lay = layout_of(*strands, *walls)
        sses = [*pairs[0], *pairs[1]]
        return [(lp.fallback, lp.points) for lp in route_loops(lay, sses, stub_bb()) if lp.a_id != "A:20-27"]

    monkeypatch.setattr(route, "_MAX_NODES", 500)
    monkeypatch.setattr(route, "_astar", counted)
    assert [fb for fb, _ in draw()] == [False, False]  # budget enough for both
    first = spends[0]
    assert len(spends) == 2 and first > 500
    monkeypatch.setattr(route, "_MAX_NODES_PAGE", first)  # enough for the first loop's second search only
    spends.clear()
    once = draw()
    assert sorted(fb for fb, _ in once) == [False, True] and spends == [first]
    assert draw() == once  # the same page spends it alike


def test_a_search_out_of_budget_spends_exactly_its_budget():
    """The page's budget is never overdrawn: a search stopped by its budget counts the nodes it expanded, its
    budget, not one more (7Y7C's domains view showed 8 000 001 of 8 000 000)."""
    a, b, wall = _walled()
    grid = _Grid(layout_of(a, b, wall))
    start, goal = grid.free_cell_along(a.c_port, a.exit_c, a), grid.free_cell_along(b.n_port, b.exit_n, b)
    spent = [0]
    assert route._astar(grid, start, 1, goal, 3, 5, spent) is None and spent == [5]


def test_no_second_search_when_the_port_stubs_cut_the_ends_apart():
    """Random-100 C2b-6 (7Y7C): a cropped gap whose two ends lie 0.57 apart on one stem's face ran four wide searches
    and spent 5.8 M of the page's 8 M nodes, though no route existed: the stubs drawn straight to its ports, which the
    search may not cross, cut its ends apart, and the free-cell test of the second search (_Router._joined) opened
    them. Driven through _route_ports (review minor 2): a port sits in the only hole of a wall, its stub fills the
    hole, the other port is beyond the wall; the test now reads the grid as searched, so no second search runs and
    the page keeps its budget."""
    lay = Layout({}, (0.0, 0.0, 10.0, 10.0))
    router = route._Router(_Grid(lay), route._element_obstacles(lay))
    grid = router.grid
    wall, row = grid.nx // 2, grid.ny // 2
    grid.blocked[wall, :] = True
    grid.blocked[wall, row] = False  # the one way through the wall, where port a sits, leaving rightwards
    pa = grid.center((wall, row))
    pb = grid.center((wall - 8, row))  # leaving leftwards, beyond the wall
    own_a = SimpleNamespace(rect=(pa[0], pa[1], pa[0], pa[1]))
    own_b = SimpleNamespace(rect=(pb[0], pb[1], pb[0], pb[1]))
    assert (wall, row) in grid.stub_cells(pa, grid.free_cell_along(pa, (1.0, 0.0), own_a))  # its stub fills the hole
    cells = route._route_ports(grid, pa, (1.0, 0.0), own_a, pb, (-1.0, 0.0), own_b, False, 10, router)
    assert cells is None
    assert router.wide_searches == 0 and router.wide_left == route._MAX_NODES_PAGE


def test_a_chain_ends_letter_a_line_runs_through_moves_round_its_stub():
    """Random-100 C2b-6 T (review zip example 08, 9FWZ --membrane auto): chain D's C letter, END_LABEL out of S31's
    stub in the stack gap under S29, stood in the D:618-629 loop that runs up into S29's N end (A6: terminus:C:D x
    loop; no route passes it). A letter a line runs through now turns round its stub's tip, or stands further out,
    to the spot lines cross least (route._clear_termini)."""
    from foldmap import figcheck

    d = draw_spec(figure_spec(Path(__file__).parent / "data" / "r100" / "9FWZ-DB.cif.gz", membrane="auto"))
    assert [h for h in figcheck.text_collisions(d) if h.startswith("terminus:")] == []
    assert d.layout.end_labels  # moved, not hidden


@pytest.mark.parametrize("first, second", [(0.1, 0.0), (0.0, 0.1)], ids=["first-higher", "second-higher"])
def test_both_turns_pulled_on_one_run_keep_both_ports_clearances(first, second):
    """Random-100 C2b-5 parked minor (a), C2b-6 review minor 1: on a route of one run between two ports at a face (up,
    across, down), both pulls move the same run (_pull_turns). It goes to the clearance that keeps both ports clear,
    the one farther out along their exits, whichever port is higher."""
    pts = [(0.0, first), (0.0, 0.5), (2.0, 0.5), (2.0, second)]  # both ports below the run, 0.1 apart in height
    out = route._pull_turns(pts, (0.0, first), (0.0, 1.0), (2.0, second), (0.0, 1.0), True, True)
    run = max(first, second) + _INFLATE
    assert out == [(0.0, first), (0.0, run), (2.0, run), (2.0, second)]


def test_a_second_ports_pull_on_five_points_is_judged_on_the_run_it_moved():
    """Random-100 C2b-6 re-review (a): on a route of five points the two pulls move different runs (the first port's
    points 1-2, the second's 3-2), so a drawing is judged on the run its own pull moved (_pull_turns), not on point 1,
    which the second port's pull leaves where it was. Here only the second port's pull applies (the first port's
    stub runs across): the drawing keeps it, its run at the second port's clearance."""
    pa, pb = (0.0, 0.0), (3.0, 0.2)
    pts = [pa, (1.0, 0.0), (1.0, 0.6), (3.0, 0.6), pb]
    out = route._pull_turns(pts, pa, (0.0, 1.0), pb, (0.0, 1.0), True, True)
    run = pb[1] + _INFLATE
    assert out == [pa, (1.0, 0.0), (1.0, run), (3.0, run), pb]


@pytest.mark.parametrize(
    "pb, run",
    [((2.0, 0.25), 0.2), ((2.0, -0.5), -0.25)],
    ids=["facing-too-close", "both-pulled"],
)
def test_no_drawing_keeps_both_clearances_falls_back_to_the_first_ports_pull(pb, run):
    """Random-100 C2b-6 re-review (b): two ports at a face whose exits point opposite ways (the first up, the second
    down) with no room for both clearances: no drawing keeps both ports clear (_pull_turns), and the run goes to the
    first port's clearance, whether only its pull applies (ports facing each other 0.25 apart, closer than two
    clearances) or both do (the second's would take the run to its own clearance, under the first port)."""
    pa = (0.0, 0.0)
    pts = [pa, (0.0, run), (2.0, run), pb]
    out = route._pull_turns(pts, pa, (0.0, 1.0), pb, (0.0, -1.0), True, True)
    assert out == [pa, (0.0, _INFLATE), (2.0, _INFLATE), pb]


def test_the_cells_a_pulled_route_uses_are_the_cells_under_it():
    """Random-100 C2b-5 parked minor (b): a route whose first turn was pulled to its port's clearance (_pull_turns)
    marked the grid cells of its search, before the pull, as used, so later routes paid to share cells it does not
    run in and not those it does. The cells marked are those under the route as drawn."""
    for k in range(10):  # every grid phase, as the slab tests sweep it
        lay, sses, _ = _slab_case("outside", k * 0.025)
        a, b = (lay.placed[s.id] for s in sses[:2])
        router = route._Router(_Grid(lay, []), route._element_obstacles(lay))
        pts = router.drawn(a.c_port, a.exit_c, a, b.n_port, b.exit_n, b)
        assert pts is not None and abs(pts[1][1] - (a.c_port[1] + _INFLATE)) < 1e-9, k  # pulled to the clearance
        grid = router.grid
        under = {c for p, q in zip(pts, pts[1:]) for c in grid.cells_along(p, q)}
        assert {tuple(c) for c in np.argwhere(grid.used > 0).tolist()} == under, k


def test_a_coarse_pages_narrow_gap_is_routed_on_a_fine_grid(monkeypatch):
    """Random-100 held-out 7BKC and 7NZ2 domains views (pages of 2269 x 116 and 430 x 358 page units): past _MAX_CELLS
    the grid is coarser than RES (0.42, 0.32), and G:550-552, standing 0.8 under G:532-544, port facing port, had no
    free cell between them: a fallback curve. A loop the page's grid cannot route is routed again on a grid at RES
    round its two ports (_fine_route), its drawing checked against every element as any route is."""
    a = bar(0, 0.0, 5.0, 4.0, -np.pi / 2)  # C port at the bottom, (0, 3)
    b = bar(20, 0.15, 1.7, 1.0, -np.pi / 2)  # N port on top, (0.15, 2.2): 0.8 under a's
    far = Placed(fake_sse("H", 100, 109, chain="B"), 4000.0, 400.0, 3.0, 0.0, HELIX_W)  # a vast page: cells of 1.03
    lay = layout_of(a, b, far)
    assert _Grid(lay).res > 1.0
    (loop,) = route_loops(lay, [a.sse, b.sse, far.sse], stub_bb())
    assert not loop.fallback
    assert np.allclose(loop.points[0], a.c_port) and np.allclose(loop.points[-1], b.n_port)
    assert route._clear(loop.points, (a, b), *route._element_obstacles(lay))
    monkeypatch.setattr(route, "_FINE_PADS", ())
    (loop,) = route_loops(layout_of(a, b, far), [a.sse, b.sse, far.sse], stub_bb())
    assert loop.fallback  # the page's grid alone finds none


def test_a_large_domain_view_routes_every_loop():
    """Random-100 held-out 7BKC's domains view (28 fallback curves, 7NZ2's 4 alike): no fallback curve; and no loop
    flagged over a label (D:183-195>D:198-201 and its copy's, through a helix name moved off them afterwards:
    _still_over)."""
    d = draw_spec(figure_spec(r100_entry("7BKC"), domains="auto"))
    assert d.stuck == [] and d.over_labels == []


def test_a_loop_is_over_a_label_only_while_the_label_is_there():
    """Random-100 held-out 7BKC's domains view: D:183-195>D:198-201 was found only with the keep-outs open, through
    D:253-260's name; _clear_helix_names then moved that name off it, and the figure's text was crossed by no line,
    yet the loop stayed flagged (A4). Loop.crosses_labels is read again once the names have moved (_still_over): a
    flagged loop through no keep-out box is over no label; one through a box stays flagged; an unflagged loop is
    never flagged."""
    from foldmap.route import Loop, _still_over

    h = bar(0, 0.0, 0.0, 4.0, 0.0)  # its name above it
    lay = layout_of(h)
    (box,) = [b for n, b in label_boxes(lay, [h.sse], route=True) if n == f"label:{h.sse.id}"]
    mid = ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2)
    through = [(mid[0], mid[1] - 2.0), (mid[0], mid[1] + 2.0)]
    beside = [(box[2] + 1.0, mid[1] - 2.0), (box[2] + 1.0, mid[1] + 2.0)]
    loops = [Loop("a", "b", beside, False, False, True), Loop("c", "d", through, False, False, True)]
    loops.append(Loop("e", "f", through, False, False, False))
    _still_over(lay, [h.sse], loops)
    assert [lp.crosses_labels for lp in loops] == [False, True, False]
