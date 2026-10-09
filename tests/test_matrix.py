"""Every structure in tests/data through the whole pipeline, in both layout modes: the invariants the
per-module tests check on one or two structures must hold on all of them (sheets, barrels, helical
bundles, membrane proteins, multi-chain assemblies, a large complex)."""

from pathlib import Path

import numpy as np
import pytest
from test_layout import overlaps
from test_route import segments_clear, self_crossing

from topoplot.cli import make_layout
from topoplot.route import route_loops

DATA = Path(__file__).parent / "data"
STRUCTURES = sorted(p.stem for p in DATA.glob("*.cif"))


@pytest.fixture(scope="module", params=[(s, m) for s in STRUCTURES for m in ("projected", "stack")],
                ids=lambda p: f"{p[0]}-{p[1]}")
def case(request):
    name, mode = request.param
    lay, sses, bb = make_layout(DATA / f"{name}.cif", mode)
    return name, mode, lay, sses, route_loops(lay, sses, bb)


def test_matrix_covers_the_test_structures():
    assert len(STRUCTURES) >= 10


def test_every_element_placed_once_and_nothing_overlaps(case):
    _, _, lay, sses, _ = case
    assert sorted(lay.placed) == sorted(s.id for s in sses)
    items = [*lay.placed.values(), *lay.ghosts]
    for i, a in enumerate(items):
        for b in items[i + 1 :]:
            assert not overlaps(a, b), (a.sse.id, b.sse.id)


def test_one_loop_per_sequence_neighbour_pair(case):
    _, _, _, sses, loops = case
    want = [(a.id, b.id) for a, b in zip(sses, sses[1:]) if a.chain == b.chain]
    assert [(l.a_id, l.b_id) for l in loops] == want


def test_routed_loops_avoid_other_elements_and_themselves(case):
    _, _, lay, _, loops = case
    for l in loops:
        if l.fallback:
            continue
        others = [p.rect for k, p in lay.placed.items() if k not in (l.a_id, l.b_id)]
        others += [g.rect for g in lay.ghosts if g.sse.id not in (l.a_id, l.b_id)]
        assert segments_clear(l.points, others), (l.a_id, l.b_id)
        assert not self_crossing(l.points), (l.a_id, l.b_id)


def test_fallback_curves_are_rare(case):
    _, _, _, _, loops = case
    if loops:
        assert sum(l.fallback for l in loops) <= max(1, 0.05 * len(loops))


def test_loops_start_and_end_at_their_ports(case):
    _, _, lay, _, loops = case
    for l in loops:
        assert np.allclose(l.points[0], lay.placed[l.a_id].c_port)
        assert np.allclose(l.points[-1], lay.placed[l.b_id].n_port)


def test_loops_keep_off_labels_and_chain_ends(case):
    from topoplot.layout import label_boxes

    _, _, lay, sses, loops = case
    boxes = label_boxes(lay, sses)
    assert len(boxes) >= sum(s.kind == "H" for s in sses)  # one per helix label, plus termini
    for l in loops:
        if not l.fallback and not l.crosses_labels:
            assert segments_clear(l.points, [b for _, b in boxes]), (l.a_id, l.b_id)
    assert sum(l.crosses_labels for l in loops) <= max(1, 0.05 * len(loops))


def test_labels_fit_inside_the_layout_bounds(case):
    from topoplot.layout import label_boxes

    _, _, lay, sses, _ = case
    x0, y0, x1, y1 = lay.bounds
    for name, (a, b, c, d) in label_boxes(lay, sses):
        assert x0 - 1e-9 <= a and c <= x1 + 1e-9 and y0 - 1e-9 <= b and d <= y1 + 1e-9, name


def test_label_keep_out_never_forces_a_fallback(case, monkeypatch):
    import topoplot.route as route

    name, mode, lay, sses, loops = case
    monkeypatch.setattr(route, "label_boxes", lambda *a: [])
    bb = make_layout(DATA / f"{name}.cif", mode)[2]
    without = {(l.a_id, l.b_id) for l in route.route_loops(lay, sses, bb) if l.fallback}
    assert {(l.a_id, l.b_id) for l in loops if l.fallback} <= without


def test_chain_end_labels_sit_clear_of_elements(case):
    from topoplot.layout import label_boxes

    _, _, lay, sses, _ = case
    rects = [p.rect for p in [*lay.placed.values(), *lay.ghosts, *lay.dna]]
    for name, box in label_boxes(lay, sses):
        if name[:2] in ("N:", "C:"):
            for r in rects:
                assert not (min(box[2], r[2]) - max(box[0], r[0]) > 1e-6 and min(box[3], r[3]) - max(box[1], r[1]) > 1e-6), name
