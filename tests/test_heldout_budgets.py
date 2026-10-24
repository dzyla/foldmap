"""Bounded work on very large pages (plan 15 Task 7, random-100 held-out 8D3D's focus view, killed after 900 s): the
second starts of stretched disulfides share one budget per layout; a coarse page's window routers draw on the page's
search budget, and so do the slab-closed searches of loops beside a membrane; a page grows for its smallest text to at
most a PDF reader's largest page, the text's font_scale counted; estimator, view and judge share one face reach."""

import logging

import numpy as np
import pytest

from foldmap import layout, route
from foldmap.layout import HELIX_W, Layout, Placed, _Item
from helpers import fake_sse

# ------------------------------------------------------------------------------------- disulfide second starts


def _bridged(k: int):
    """k pairs of one-helix items 30 apart, each pair bridged by a disulfide, every item on its own chain."""
    items, sses, bridges = [], [], []
    for n in range(k):
        pair = []
        for side, x in (("a", 0.0), ("b", 30.0)):
            s = fake_sse("H", 1, 10, chain=f"{n}{side}")
            p = Placed(s, 0.0, 0.0, 3.0, 0.0, HELIX_W)
            items.append(_Item([p], [(0.0, 0.0)], x, 10.0 * n, (x, 10.0 * n), 1, s.chain))
            sses.append(s)
            pair.append(s.id)
        bridges.append((pair[0], 0.5, pair[1], 0.5))
    return items, sses, bridges


def _second_starts(monkeypatch, budget):
    """Bridged pairs held near home (a strong home pull, so each bridge stays stretched): how many L-BFGS runs
    _attract made (4 for the first minimum, 3 per side of each second start), and the debug record."""
    calls = []
    real = layout._lbfgs

    def counted(fun, x0, *a, **kw):
        calls.append(1)
        return real(fun, x0, *a, **kw)

    monkeypatch.setattr(layout, "_K_HOME", 50.0)
    monkeypatch.setattr(layout, "_lbfgs", counted)
    monkeypatch.setattr(layout, "_SS_RESTART_PAIRS", budget)
    items, sses, bridges = _bridged(3)
    layout._attract(items, sses, {}, bridges=bridges)
    return len(calls), [(it.cx, it.cy) for it in items]


def test_the_second_starts_of_one_layout_share_a_budget(monkeypatch, caplog):
    """Random-100 held-out 8D3D's focus view (896 items, 848 bridges): every disulfide still stretched after the
    first minimum was started again from five places, each a solve of the whole page (about 2 s there), with no
    bound: killed after 900 s. The second starts spend at most _SS_RESTART_PAIRS pair evaluations (each energy
    evaluation counted as its items squared), counted as they run; past it a bridge keeps the first minimum. Within
    the budget nothing changes."""
    whole, placed = _second_starts(monkeypatch, 10**18)
    assert whole == 4 + 3 * 5 * 3  # three stretched bridges, five sides each
    assert _second_starts(monkeypatch, 10**18) == (whole, placed)  # the same layout every time
    caplog.set_level(logging.DEBUG, logger="foldmap.layout")
    assert _second_starts(monkeypatch, 0)[0] == 4  # none: the first minimum
    assert _second_starts(monkeypatch, 1)[0] == 4 + 3 * 5  # the first bridge's start spends it; the rest keep theirs
    assert any("not for the budget" in r.getMessage() for r in caplog.records)


def test_a_page_of_many_bridges_keeps_its_second_starts_within_the_budget():
    """At 8D3D's size (896 items) the budget allows a few bridges' second starts of its hundreds, not every one: a
    second start of the whole page costs about 750 evaluations of 896 x 896 pairs, five per bridge."""
    per_bridge = 5 * 750 * 896 * 896
    assert 1 <= layout._SS_RESTART_PAIRS / per_bridge <= 10


@pytest.mark.parametrize("seed", range(6))
def test_the_overlap_pairs_found_by_sorting_are_the_dense_ones_in_order(seed):
    """_x_overlaps finds the pairs overlapping in x among the items within each one's reach only (sorted by x, not
    all n x n): the same pairs, in the same row-major order as np.nonzero(np.triu(...)) gave, so every sum over them
    and every minimum is as before; clustered, tied and spread positions alike."""
    rng = np.random.default_rng(seed)
    n = 300
    h = np.round(rng.uniform(0.2, 3.0, (n, 2)), 4)
    want = h[:, None, :] + h[None, :, :] + layout.MARGIN
    joined = rng.random((n, n)) < 0.01
    joined |= joined.T
    want[joined] = (h[:, None, :] + h[None, :, :])[joined]
    spread = (5.0, 60.0, 400.0)[seed % 3]
    Q = np.round(rng.uniform(0.0, spread, (n, 2)), 4)
    Q[: n // 4, 0] = Q[0, 0]  # a quarter of them on one x
    dense = np.nonzero(np.triu(want[..., 0] > np.abs(Q[None, :, 0] - Q[:, None, 0]), 1))
    got = layout._x_overlaps(Q, want, want[..., 0].max(axis=1))
    assert np.array_equal(got[0], dense[0]) and np.array_equal(got[1], dense[1])


def test_outlines_are_asked_only_of_the_points_within_their_extent(monkeypatch):
    """render._covered and visible_pieces ask a body's outline only of the points within its vertices' extent (8D3D's
    focus view asked about 1,800 bodies of every point of 848 bars): what is covered, and drawn, is as before."""
    from matplotlib.path import Path as MPath

    from foldmap import render

    rng = np.random.default_rng(3)
    shapes = []
    for _ in range(40):
        c, r = rng.uniform(0, 20, 2), rng.uniform(0.3, 2.0)
        t = np.linspace(0, 2 * np.pi, 7)[:-1] + rng.uniform(0, 1)
        ring = MPath(np.c_[c[0] + r * np.cos(t), c[1] + r * np.sin(t)])
        shapes.append(MPath.circle(c, r) if rng.random() < 0.3 else ring)
    pts = rng.uniform(-1, 21, (500, 2))
    routes = [[(0.0, 0.0), (20.0, 20.0), (20.0, 0.0)], [(-1.0, 10.0), (21.0, 9.0)], [(3.0, 21.0), (3.0, -1.0)]]
    halos = [(4.0, 4.0, 6.0, 7.0)]
    fast = render._covered(pts, shapes), [render.visible_pieces(r, shapes, halos) for r in routes]
    monkeypatch.setattr(render, "_within_extents", lambda p, paths: np.ones((len(p), len(paths)), bool))
    every = render._covered(pts, shapes), [render.visible_pieces(r, shapes, halos) for r in routes]
    assert fast[0].any() and np.array_equal(fast[0], every[0])
    assert [[q.tolist() for q in r] for r in fast[1]] == [[q.tolist() for q in r] for r in every[1]]


def test_an_outline_of_no_vertices_covers_nothing():
    """Task 7 review: a Path of no vertices among the bodies (last, first or between) made the extents' reduceat raise
    IndexError; it holds no point, and the others are asked as before."""
    from matplotlib.path import Path as MPath

    from foldmap import render

    square = MPath([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)])
    empty = MPath(np.zeros((0, 2)))
    pts = np.array([[1.0, 1.0], [3.0, 3.0]])
    for shapes in ([square, empty], [empty, square], [empty, square, empty], [empty]):
        within = render._within_extents(pts, shapes)
        assert within.shape == (2, len(shapes)) and not within[:, [s is empty for s in shapes]].any()
        assert render._covered(pts, shapes).tolist() == [square in shapes, False]
        render.visible_pieces([(-1.0, 1.0), (3.0, 1.0)], shapes)
    assert render._within_extents(pts, []).shape == (2, 0)


# ----------------------------------------------------------------------------------------- one search budget


def bar(start, cx, cy, length, angle, chain="A"):
    return Placed(fake_sse("H", start, start + 9, chain=chain), cx, cy, length, angle, HELIX_W, "a")


def layout_of(*placed):
    rects = np.array([p.rect for p in placed])
    b = (rects[:, 0].min(), rects[:, 1].min(), rects[:, 2].max(), rects[:, 3].max())
    return Layout({p.sse.id: p for p in placed}, tuple(float(x) for x in b))


def _coarse_gap():
    """A loop between two helices 4 apart on a page so large its grid is coarse (test_route's vast page)."""
    a = bar(0, 0.0, 5.0, 4.0, -np.pi / 2)
    b = bar(20, 4.0, 1.7, 1.0, -np.pi / 2)
    far = Placed(fake_sse("H", 100, 109, chain="B"), 4000.0, 400.0, 3.0, 0.0, HELIX_W)
    return a, b, layout_of(a, b, far)


def test_the_window_routers_of_a_coarse_page_draw_on_its_budget(monkeypatch):
    """Plan 15 Task 7 (carried from Task 2): route._fine_route gave each window router a fresh _MAX_NODES_PAGE, so a
    page of many loops its coarse grid cannot route had no bound on its searches. A window's router starts on what
    the page's windows have left of their own pool (window_left, _MAX_NODES_WINDOWS) and gives back what it spent:
    with nothing left a window runs no second search. The pool is not the page's own budget (plan 15 final review:
    drawn from it, already spent on a huge page, 7NZ2's domains view lost three loops to fallback curves): a page
    that has spent its own still routes its windows."""
    a, b, lay = _coarse_gap()
    monkeypatch.setattr(route, "_MAX_NODES", 0)  # no first search finds the route: only a second search can
    keep = layout.label_boxes(lay, [a.sse, b.sse], route=True)
    elements = route._element_obstacles(lay)

    def fine(page_left, window_left):
        page = route._Router(route._Grid(lay, [r for _, r in keep]), elements)
        page.wide_left, page.window_left = page_left, window_left
        pts, _ = route._fine_route(lay, keep, elements, a, b, (set(), set()), [], page)
        return pts, page

    pts, page = fine(0, route._MAX_NODES_WINDOWS)  # the page's own budget spent: the windows' pool is not
    spent = route._MAX_NODES_WINDOWS - page.window_left
    assert pts is not None and spent > 0 and page.wide_searches >= 1 and page.wide_left == 0
    pts, page = fine(route._MAX_NODES_PAGE, 0)
    assert pts is None and page.window_left == 0 and page.wide_searches == 0
    assert page.wide_left == route._MAX_NODES_PAGE  # nor do the windows touch the page's own


def test_a_slab_closed_search_is_counted_and_drawn_from_the_page_budget():
    """Task 6's slab-closed first search (ports on one side of every band) ran on its own 200 000 nodes, counted
    nowhere: it is counted in the router's stats and drawn from the page's budget; with that spent it is not run."""
    a = bar(0, 0.0, 1.0, 4.0, np.pi / 2)
    b = bar(10, 10.0, 1.0, 4.0, -np.pi / 2)
    wall = bar(40, 5.0, -27.2, 65.6, np.pi / 2)
    lay = layout_of(a, b, wall)
    band = (5.0, 12.0)
    lay.membrane = {"y": band, "bands": [band], "inside": "below"}

    def routed(left):
        router = route._Router(route._Grid(lay), route._element_obstacles(lay))
        router.wide_left = left
        pts, _ = router.route(a.c_port, a.exit_c, a, b.n_port, b.exit_n, b)
        return router, pts

    router, pts = routed(route._MAX_NODES_PAGE)
    assert router.closed_searches == 1 and router.closed_nodes > 0
    assert router.wide_left == route._MAX_NODES_PAGE - router.closed_nodes
    assert max(y for _, y in pts) < band[0]  # under the wall, out of the slab
    router, _ = routed(0)
    assert router.closed_searches == 0 and router.closed_nodes == 0


# ------------------------------------------------------------------------------------------------ page size


def test_a_page_grows_for_its_smallest_text_to_a_pdf_readers_largest_page():
    """Carried from Task 2: page_scale's 1 pt floor grows a page past MAX_PAGE_IN without bound; PDF readers open no
    page past 200 in (matplotlib writes no UserUnit): held at MAX_GROWN_IN, the notes say its text is under 1 pt."""
    from foldmap.layout import _FLOOR_SCALE, MAX_GROWN_IN, MAX_PAGE_IN, page_fit, page_scale

    assert 2263.0 * page_scale(2263.0, 110.0) == pytest.approx(2263.0 * _FLOOR_SCALE)  # 7BKC's: 154 in, as before
    assert 6000.0 * page_scale(6000.0, 300.0) == pytest.approx(MAX_GROWN_IN)
    assert 6000.0 * page_scale(300.0, 6000.0) == pytest.approx(MAX_GROWN_IN)  # the longer side, either way
    note = page_fit(MAX_GROWN_IN, 10.0, 0.9)
    assert "held at 200 in" in note and "under 1 pt" in note
    assert "past 48 in, so that no text is under 1 pt" in page_fit(154.0, 9.0, 1.6)
    assert "held at 48 in" in page_fit(MAX_PAGE_IN, 9.0, 5.0)


def test_the_floor_counts_the_text_as_drawn():
    """Carried from Task 2: a Style.font_scale under 1 (the journal theme's 0.9) draws the smallest text smaller than
    the floor was set for; the floor is for the text as drawn. A font_scale of 1 or more leaves it as it was."""
    from foldmap.layout import _EM_TUNED, MPL_MIN_PT, SMALLEST_TEXT, page_scale

    smallest = SMALLEST_TEXT * _EM_TUNED * 72
    assert smallest * page_scale(2263.0, 110.0, 0.9) * 0.9 == pytest.approx(MPL_MIN_PT)
    assert page_scale(2263.0, 110.0, 1.4) == page_scale(2263.0, 110.0)
    assert page_scale(40.0, 30.0, 0.9) == page_scale(40.0, 30.0)  # a page the floor does not reach


# ------------------------------------------------------------------------------------------------ membrane


def test_estimator_view_and_judge_share_one_face_reach():
    """Task 6 review: membrane.FACE_REACH replaced the literal 2 Å in cli._transmembrane, cli._face_medians and
    figcheck.membrane_view. Plan 26 (ruling R1): the membrane view, A18 and foldmap data count a transmembrane element
    by one rule, membrane.crosses (the CA run crosses the mid-plane and spans half the thickness); FACE_REACH is the
    estimate's own (membrane._spans)."""
    import inspect

    from foldmap import cli, datafile, figcheck, membrane

    assert cli.crosses is membrane.crosses and figcheck.crosses is membrane.crosses
    for f in (cli._transmembrane, cli._face_medians, figcheck.membrane_view, datafile._membrane):
        src = inspect.getsource(f)
        assert "crosses(" in src or "crossing(" in src, f.__name__
        assert "half + 2" not in src and "half - 2" not in src
    for f in (cli._transmembrane, cli._face_medians, figcheck.membrane_view):
        assert "FACE_REACH" not in inspect.getsource(f), f.__name__
    assert "HALF + FACE_REACH" in inspect.getsource(membrane._spans)
