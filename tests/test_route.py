from types import SimpleNamespace

import numpy as np
import pytest

from foldmap.frame import view_frame
from foldmap.layout import HELIX_W, PITCH, Layout, Placed, build_layout
from foldmap.route import route_loops
from helpers import fake_sse, pipeline

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
    for key in [("B:456-478", "B:480-485"), ("F:456-472", "F:474-485")]:
        l = loops[key]
        assert not l.fallback and len(l.points) == 2, key  # one straight link across the kink
        assert np.allclose(l.points[0], lay.placed[key[0]].c_port) and np.allclose(
            l.points[1], lay.placed[key[1]].n_port
        )
