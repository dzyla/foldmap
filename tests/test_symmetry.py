import numpy as np
import pytest

from foldmap.symmetry import detect_symmetry
from helpers import pipeline

DATA = __import__("pathlib").Path(__file__).parent / "data"


def sym(name):
    bb, sses, _ = pipeline(DATA / f"{name}.cif")
    return bb, sses, detect_symmetry(bb, sses)


@pytest.mark.parametrize("name, n", [("8UUP", 3), ("2HHB", 2), ("1TIM", 2), ("1LMB", 2)])
def test_cyclic_assemblies_are_found(name, n):
    _, _, s = sym(name)
    assert s is not None and s.n == n and s.rmsd < 3.0
    assert abs(np.linalg.norm(s.axis) - 1) < 1e-9


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "5NKT", "1ZAA"])
def test_single_chains_have_no_symmetry(name):
    assert sym(name)[2] is None


def test_hemoglobin_protomers_pair_alpha_with_its_beta():
    _, _, s = sym("2HHB")
    assert sorted(map(sorted, s.protomers)) == [["A", "B"], ["C", "D"]]


def test_8uup_protomers_are_heterodimers_in_rotation_order():
    _, _, s = sym("8UUP")
    assert all(len(p) == 2 for p in s.protomers)
    assert sorted(c for p in s.protomers for c in p) == list("ABCDEF")


def _mates(sses, s):
    """Element of protomer 0 -> its copy in protomer k (same chain class, same residue numbers)."""
    cls = {c: (k, i) for k, p in enumerate(s.protomers) for i, c in enumerate(p)}
    by = {}
    for e in sses:
        if e.chain in cls:
            k, i = cls[e.chain]
            by[(k, i, e.first.seq, e.last.seq)] = e
    return [
        [by[(k, i, a, b)] for k in range(s.n)]
        for (k0, i, a, b) in list(by)
        if k0 == 0
        if all((k, i, a, b) in by for k in range(s.n))
    ]


@pytest.mark.parametrize("name", ["8UUP", "2HHB", "1TIM"])
def test_unrolled_projection_draws_every_protomer_the_same(name):
    from foldmap.frame import symmetric_frame

    bb, sses, s = sym(name)
    f = symmetric_frame(sses, s)
    orbits = _mates(sses, s)
    assert len(orbits) >= 3
    steps = []
    for orbit in orbits:
        xy = np.array([f.project(e.centroid)[0] for e in orbit])
        d = np.array([f.direction(e.centroid, e.axis) for e in orbit])
        steps.append(np.diff(xy, axis=0))
        assert np.allclose(d, d[0], atol=0.25)  # same orientation in every copy
    steps = np.array(steps)  # orbit, copy step, xy
    period = np.median(steps[:, :, 0])
    assert period > 5  # copies side by side (Å)
    assert np.median(np.abs(steps[:, :, 0] - period)) < 2.0 and np.median(np.abs(steps[:, :, 1])) < 2.0


def test_unrolled_frame_is_a_rotation_everywhere():
    from foldmap.frame import symmetric_frame

    _, sses, s = sym("8UUP")
    f = symmetric_frame(sses, s)
    for e in sses[::7]:
        u, v, w = f.basis(e.centroid)
        m = np.array([u, v, w])
        assert np.allclose(m @ m.T, np.eye(3), atol=1e-9) and np.linalg.det(m) == pytest.approx(1.0)


@pytest.mark.parametrize("name", ["8UUP", "2HHB", "1TIM"])
def test_symmetric_layout_repeats_the_protomer(name):
    from foldmap.cli import make_layout

    bb, sses, s = sym(name)
    lay, _, _ = make_layout(DATA / f"{name}.cif")
    orbits = [o for o in _mates(sses, s) if all(e.id in lay.placed for e in o)]
    rel = []
    for orbit in orbits:
        p = [lay.placed[e.id] for e in orbit]
        rel.append([(q.cx - p[0].cx, q.cy - p[0].cy) for q in p[1:]])
        assert len({round(q.angle % (2 * np.pi), 6) for q in p}) == 1  # same direction in each copy
    rel = np.array(rel)  # orbit, copy, xy
    typical = np.median(rel, axis=0)
    assert (np.abs(rel - typical).max(axis=2) < 1.5).mean() > 0.8  # most elements sit at the same place in each copy


@pytest.mark.parametrize("name", ["8UUP", "2HHB", "1TIM", "4CH9", "9TX7"])
def test_separation_keeps_every_copy_where_its_mate_is(name, monkeypatch):
    """Overlaps are cleared by moving each element together with its symmetry copies, so after separation
    every copy still sits at the same offset from its protomer-0 mate as the others of its protomer."""
    from foldmap import layout
    from foldmap.cli import make_layout

    seen = {}
    real = layout._symmetry_ties

    def spy(items, frame):
        seen["items"], seen["ties"] = items, real(items, frame)
        return seen["ties"]

    monkeypatch.setattr(layout, "_symmetry_ties", spy)
    make_layout(DATA / f"{name}.cif", assembly="asu" if name == "4CH9" else "auto")
    items, ties = seen["items"], seen["ties"]
    offsets: dict[int, list[tuple[float, float]]] = {}
    for i, (j, k) in ties.items():
        offsets.setdefault(k, []).append((items[i].cx - items[j].cx, items[i].cy - items[j].cy))
    assert offsets
    for o in offsets.values():
        assert np.ptp(np.array(o), axis=0).max() < 1e-9


def _box(x, y, start, chain="A", tall=2.0, angle=np.pi / 2):
    """A one-strand item at (x, y), 1 wide and `tall` long (upright; angle 0: lying)."""
    from foldmap.layout import Placed, _Item
    from helpers import fake_sse

    s = fake_sse("E", start, start + 6, chain=chain)
    return _Item([Placed(s, 0.0, 0.0, tall, angle, 1.0)], [(0.0, 0.0)], x, y, (x, y), start, chain)


def _area(items):
    x0 = min(it.cx + it.box[0] for it in items)
    y0 = min(it.cy + it.box[1] for it in items)
    return (max(it.cx + it.box[2] for it in items) - x0) * (max(it.cy + it.box[3] for it in items) - y0)


def _overlaps(items):
    from foldmap.layout import _collide

    h, m = [it.half for it in items], [it.mid for it in items]
    return sum(_collide(items[i], h[i], m[i], items[j], h[j], m[j]) for i in range(len(items)) for j in range(i))


def _pushed_apart(monkeypatch, items, ties):
    """Run _separate; how many overlaps were left for the fallback slide."""
    from foldmap import layout

    left = []
    real = layout._slide_clear
    monkeypatch.setattr(layout, "_slide_clear", lambda its, *a: left.append(_overlaps(its)) or real(its, *a))
    layout._separate(items, ties)
    assert _overlaps(items) == 0
    return left[0]


def test_a_push_in_one_protomer_moves_the_copies_too(monkeypatch):
    """Element 0 overlaps an untied neighbour (1); its copy (2) overlaps nothing but moves the same way."""
    items = [_box(0.0, 0.0, 1), _box(0.4, 0.0, 20), _box(10.0, 0.0, 1, "B")]
    assert _pushed_apart(monkeypatch, items, {2: (0, 1)}) == 0
    assert items[0].cx != 0.0
    assert (items[2].cx - items[0].cx, items[2].cy - items[0].cy) == pytest.approx((10.0, 0.0), abs=1e-12)


def test_without_ties_each_item_moves_alone(monkeypatch):
    items = [_box(0.0, 0.0, 1), _box(0.4, 0.0, 20), _box(10.0, 0.0, 1, "B")]
    assert _pushed_apart(monkeypatch, items, {}) == 0
    assert (items[2].cx, items[2].cy) == (10.0, 0.0)


def test_copies_of_one_element_that_overlap_are_pushed_apart_alone(monkeypatch):
    """Moving a whole orbit cannot part two of its own copies: one above the other, they move apart alone."""
    items = [_box(0.0, 0.0, 1, tall=4.0, angle=0.0), _box(0.3, 0.2, 1, "B", tall=4.0, angle=0.0)]
    assert _pushed_apart(monkeypatch, items, {1: (0, 1)}) == 0
    assert items[0].cx == 0.0 and items[1].cx == 0.3  # a push up and down only
    assert items[0].cy < 0.0 < 0.2 < items[1].cy


def test_overlapping_copies_side_by_side_widen_the_spread(monkeypatch):
    """Two copies of one element overlapping sideways: protomer 1 moves further along, the copy goes with it."""
    items = [_box(0.0, 0.0, 1), _box(0.3, 0.0, 1, "B")]
    assert _pushed_apart(monkeypatch, items, {1: (0, 1)}) == 0
    assert items[0].cx == 0.0 and items[1].cx > 0.3


def test_a_sideways_overlap_between_protomers_widens_the_spread(monkeypatch):
    """Protomer 1's copy of `a` touches protomer 0's `b`: protomer 1 moves over as a whole."""
    items = [_box(0.0, 0.0, 1), _box(5.0, 0.0, 20), _box(6.6, 0.0, 1, "B"), _box(11.6, 0.0, 20, "B")]
    assert _pushed_apart(monkeypatch, items, {2: (0, 1), 3: (1, 1)}) == 0
    assert items[2].cx - items[0].cx == pytest.approx(items[3].cx - items[1].cx, abs=1e-12)
    assert items[1].cx - items[0].cx == pytest.approx(5.0, abs=1e-12)  # each protomer itself unchanged
    assert items[2].cx - items[0].cx > 6.6


def test_copies_move_alone_when_moving_them_together_cannot_clear(monkeypatch):
    """A free element (u) wedged between protomer 1's copies, where protomer 0 has room: the symmetric pushes
    go round in circles. After _SEP_ORBIT_SWEEPS sweeps the items start again from where they were and move
    one by one."""
    from foldmap import layout

    sweeps = []
    real = layout._overlapping
    monkeypatch.setattr(layout, "_overlapping", lambda *a: sweeps.append(1) or real(*a))
    a, b = (0.8, 0.3), (3.4, 3.4)
    items = [_box(*a, 1), _box(*b, 20), _box(a[0] + 3, a[1], 1, "B"), _box(b[0] + 3, b[1], 20, "B"), _box(5.3, 2.8, 40)]
    alone = [_box(it.cx, it.cy, it.order, it.chain) for it in items]
    assert _pushed_apart(monkeypatch, items, {2: (0, 1), 3: (1, 1)}) == 0
    assert len(sweeps) > layout._SEP_ORBIT_SWEEPS
    layout._separate(alone)
    assert [(it.cx, it.cy) for it in items] == [(it.cx, it.cy) for it in alone]  # as if there were no copies


def test_symmetric_separation_never_costs_much_more_room_than_moving_items_alone():
    """Four copies of an 11-element protomer, one spread apart, and four free elements among them (a random
    stress case): the symmetric pushes widen the spread to four times its start. The figure keeps the
    symmetric result only if it is at most 5 % larger than separating the items one by one."""
    from foldmap import layout

    protomer = [(2.48, 9.39, 3.83, 0), (4.71, 4.16, 3.11, 0), (-0.75, 7.38, 3.56, 1), (5.13, 0.55, 1.63, 0)]
    protomer += [(5.77, 6.83, 2.73, 1), (7.57, 2.26, 2.48, 1), (1.82, 1.91, 1.51, 1), (4.78, 2.35, 2.96, 0)]
    protomer += [(7.58, 5.51, 2.41, 1), (0.41, 3.97, 3.87, 0), (1.85, 6.68, 2.81, 1)]
    free = [(23.11, 8.18, 3.02, 1), (48.59, 6.31, 1.39, 1), (25.27, 5.35, 1.53, 0), (36.83, 9.35, 2.63, 0)]
    spread, lift = 11.21, -0.64

    def build():
        items, ties = [], {}
        for k in range(4):
            for a, (x, y, tall, up) in enumerate(protomer):
                items.append(_box(x + k * spread, y + k * lift, 10 * a + 1, "ABCD"[k], tall, np.pi / 2 * up))
                if k:
                    ties[len(items) - 1] = (a, k)
        items += [_box(x, y, 500 + 10 * u, "Z", tall, np.pi / 2 * up) for u, (x, y, tall, up) in enumerate(free)]
        return items, ties

    items, ties = build()
    alone, _ = build()
    layout._separate(items, ties)
    layout._separate(alone)
    assert _overlaps(items) == 0
    assert _area(items) <= 1.05 * _area(alone)


def test_symmetry_can_be_switched_off(tmp_path):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(DATA / "2HHB.cif"), "-o", str(out), "--symmetry", "off"]) == 0
