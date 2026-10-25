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


@pytest.mark.parametrize("name", ["8UUP", "1TIM"])
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


@pytest.mark.parametrize("name", ["2HHB", "25IN"])
def test_a_two_fold_assembly_is_drawn_as_one_copy_and_its_mirror_image(name):
    """Plan 18: a C2 assembly is seen from the side of its dyad, the dyad up the page and protomer 0 on the left;
    every element of protomer 1 sits at its mate's mirror image about one vertical line (x -> 2 * mirror_x - x, the
    same y), its direction mirrored. With two_fold=mirror, whatever the protomers hold (plan 20 Task 1)."""
    from foldmap.cli import make_layout
    from foldmap.frame import Frame
    from foldmap.style import Style

    bb, sses, s = sym(name)
    assert s.kind == "C" and s.n == 2
    lay, _, _ = make_layout(DATA / f"{name}.cif", look=Style(two_fold="mirror"))
    assert type(lay.frame) is Frame and lay.frame.mirror_x is not None
    assert abs(float(lay.frame.v @ s.axis)) == pytest.approx(1.0)  # the dyad up the page
    m = lay.frame.mirror_x
    pairs = [o for o in _mates(sses, s) if all(e.id in lay.placed for e in o)]
    assert len(pairs) >= 8
    for a, b in pairs:
        p, q = lay.placed[a.id], lay.placed[b.id]
        assert q.cx == pytest.approx(2 * m - p.cx, abs=1e-6) and q.cy == pytest.approx(p.cy, abs=1e-6)
        assert np.cos(q.angle) == pytest.approx(-np.cos(p.angle), abs=1e-6)
        assert np.sin(q.angle) == pytest.approx(np.sin(p.angle), abs=1e-6)
    assert np.mean([lay.placed[a.id].cx for a, _ in pairs]) < m  # protomer 0 on the left


def test_two_fold_unrolled_draws_the_copies_side_by_side():
    """Plan 18 Task 2: TrmD (1UAM, the author's example) by default seen from the side of its dyad, one copy the
    other's mirror image; with two_fold=unrolled its copies side by side around the axis, as before plan 18 (a
    cylinder frame, no copy tied as a mirror image)."""
    from foldmap.cli import make_layout
    from foldmap.frame import CylinderFrame, Frame
    from foldmap.style import Style

    uam = DATA / "assemblies" / "1UAM.cif.gz"
    mirror, _, _ = make_layout(uam)
    assert type(mirror.frame) is Frame and mirror.frame.mirror_x is not None
    assert any(m for _, m in mirror.copy_ties.values()) and mirror.two_fold_assembly
    unrolled, _, _ = make_layout(uam, look=Style(two_fold="unrolled"))
    assert type(unrolled.frame) is CylinderFrame and unrolled.two_fold_assembly
    assert unrolled.copy_ties and not any(m for _, m in unrolled.copy_ties.values())


def test_two_fold_auto_mirrors_only_where_a_protomer_is_one_chain():
    """Plan 20 Task 1: by default (two_fold=auto) haemoglobin (2HHB, protomers α1β1 | α2β2) is unrolled side by
    side, as before plan 18, so its α1β2 interface stays together (a mirror line would put α1 and β2 in opposite
    corners); TrmD (1UAM, one chain a protomer) is drawn as one copy and its mirror image."""
    from foldmap.cli import _mirrors, make_layout
    from foldmap.frame import CylinderFrame, Frame

    bb, _, s = sym("2HHB")
    assert not _mirrors("auto", s, bb) and _mirrors("mirror", s, bb) and not _mirrors("unrolled", s, bb)
    held = __import__("types").SimpleNamespace(protomers=[["A", "X"], ["C", "Y"]])  # X, Y: nucleic-acid chains
    assert _mirrors("auto", held, bb)  # one protein chain a protomer: nucleic-acid chains do not count
    hb, _, _ = make_layout(DATA / "2HHB.cif")
    assert type(hb.frame) is CylinderFrame and hb.two_fold_assembly
    assert not any(m for _, m in hb.copy_ties.values())
    uam, _, _ = make_layout(DATA / "assemblies" / "1UAM.cif.gz")
    assert type(uam.frame) is Frame and uam.frame.mirror_x is not None
    assert any(m for _, m in uam.copy_ties.values())


def _turned(monkeypatch, path):
    """make_layout of `path`, with what _dyad_turn was given and chose."""
    import foldmap.cli as cli

    seen = {}
    real = cli._dyad_turn

    def turn(*args):
        seen["args"], seen["turn"] = args, real(*args)
        return seen["turn"]

    monkeypatch.setattr(cli, "_dyad_turn", turn)
    lay, sses, _ = cli.make_layout(path)
    return lay, sses, seen


def test_the_dyad_view_is_turned_to_keep_the_interface_together(monkeypatch):
    """Plan 19 Task 1: TrmD (1UAM) seen square-on to the line between its protomers piles the interface (α1, α5, the
    H/I hairpins of both chains) on the mirror line, and the layout scattered it (cross-chain elements closer than
    16 Å drawn at a median 4.6 times their 3D distance, rank correlation of drawn and 3D distances 0.35). Turned
    75° about the dyad, the projection keeps the distances' order best and the figure keeps the interface together;
    protomer 0's centre stays left of the mirror line."""
    from foldmap.layout import SCALE

    lay, sses, seen = _turned(monkeypatch, DATA / "assemblies" / "1UAM.cif.gz")
    assert seen["turn"] == 75 and lay.frame.mirror_x is not None
    placed = [s for s in sses if s.id in lay.placed]
    xy = {s.id: np.array([lay.placed[s.id].cx, lay.placed[s.id].cy]) / SCALE for s in placed}
    pairs = [(a, b) for i, a in enumerate(placed) for b in placed[i + 1 :]]
    d3 = np.array([np.linalg.norm(a.centroid - b.centroid) for a, b in pairs])
    d2 = np.array([np.linalg.norm(xy[a.id] - xy[b.id]) for a, b in pairs])
    near = np.array([a.chain != b.chain for a, b in pairs]) & (d3 < 16)
    assert near.sum() >= 10 and np.median(d2[near] / d3[near]) < 1.25
    ranks = [np.argsort(np.argsort(d)) for d in (d3, d2)]
    assert np.corrcoef(*ranks)[0, 1] >= 0.8
    first = [s.centroid for s in sses if s.chain == "A1"]
    assert lay.frame.project(np.mean(first, axis=0))[0, 0] * SCALE < lay.frame.mirror_x


def test_the_dyad_turn_is_the_same_on_any_cpu(monkeypatch):
    """The turn is chosen on snapped distances with tied ranks averaged: noise in the last bits of the centroids and
    axes (BLAS kernels) never changes it. Turns tied within DYAD_TURN_GAIN go to the smallest."""
    import foldmap.cli as cli

    _, _, seen = _turned(monkeypatch, DATA / "assemblies" / "1UAM.cif.gz")
    cents, up, right, face = seen["args"]
    rng = np.random.default_rng(0)
    for _ in range(4):
        nudged = [v * (1 + rng.uniform(-1e-15, 1e-15, v.shape)) for v in (cents, up, right, face)]
        assert cli._dyad_turn(*nudged) == 75
    monkeypatch.setattr(cli, "DYAD_TURN_GAIN", 1.0)
    assert cli._dyad_turn(cents, up, right, face) == 0


@pytest.mark.parametrize("name", ["8UUP", "2HHB", "1TIM", "4CH9", "9TX7"])
def test_separation_keeps_every_copy_where_its_mate_is(name, monkeypatch):
    """Overlaps are cleared by moving each element together with its symmetry copies, so after separation
    every copy still sits at the same offset from its protomer-0 mate as the others of its protomer (a mirrored copy,
    plan 18: at its mate's mirror image, the sum of their x and the difference of their y the same for all)."""
    from foldmap import layout
    from foldmap.cli import make_layout

    seen = {}
    real = layout._symmetry_ties

    def spy(items, frame, *rest):
        seen["items"], seen["ties"], seen["mirror"] = items, real(items, frame, *rest), rest[0]
        return seen["ties"]

    monkeypatch.setattr(layout, "_symmetry_ties", spy)
    make_layout(DATA / f"{name}.cif", assembly="asu" if name == "4CH9" else "auto")
    items, ties = seen["items"], seen["ties"]
    offsets: dict[int, list[tuple[float, float]]] = {}
    for i, (j, k) in ties.items():
        x = items[i].cx + items[j].cx if seen["mirror"] else items[i].cx - items[j].cx
        offsets.setdefault(k, []).append((x, items[i].cy - items[j].cy))
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
    monkeypatch.setattr(layout, "_slide_clear", lambda its, *a, **k: left.append(_overlaps(its)) or real(its, *a, **k))
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
    one by one (with no second chance, _separate's `again` False: as if there were no copies)."""
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
    stress case): the symmetric pushes widen the spread to four times its start, 3.9 times the one-by-one room. The
    figure keeps that symmetric result only if it is at most 5 % larger than separating the items one by one; else
    the second chance (_tied_again, C2b-3), the four free elements yielding to the copies, keeps them alike in
    1.15 times the room, within _SEP_TIE_ROOM."""
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
    assert layout._separate(items, ties, again=True) == "wide"
    layout._separate(alone)
    assert _overlaps(items) == 0
    assert 1.05 * _area(alone) < _area(items) <= layout._SEP_TIE_ROOM * _area(alone)
    for i, (j, _) in ties.items():
        k = ties[i][1]
        assert (items[i].cx - items[j].cx, items[i].cy - items[j].cy) == pytest.approx(
            (items[11 * k].cx - items[0].cx, items[11 * k].cy - items[0].cy), abs=1e-12
        )


def test_mirror_images_are_kept_when_they_cost_a_little_more_room():
    """9GEP (random-100 A15 E~A 1.17, F~B 0.81, batch C2a): a dimer drawn as mirror images about a dyad in the page,
    crowded by untied chains, cleared its overlaps symmetrically only after about 3,000 passes and in 1.07 times the
    one-by-one area; both over the limits, so its copies were moved one by one and drawn unalike. Mirror images
    (the figure says so) get as many passes as moving one by one and up to _SEP_MIRROR_ROOM. Three elements and
    their mirror images (a random case found to need 1.058 times the room) stay exact mirror images."""
    from foldmap import layout

    base = [(-1.46, 0.78, 2.5, 1), (-0.34, -0.17, 1.1, 1), (-0.7, -0.41, 3.8, 0)]

    def build():
        items = [_box(x, y, 10 * a + 1, "A", t, np.pi / 2 * u) for a, (x, y, t, u) in enumerate(base)]
        items += [_box(-x, y, 10 * a + 1, "B", t, np.pi / 2 * u) for a, (x, y, t, u) in enumerate(base)]
        return items, {3 + a: (a, 1) for a in range(3)}

    items, ties = build()
    alone, _ = build()
    layout._separate(items, ties, mirror=0.0)
    layout._separate(alone)
    assert _overlaps(items) == 0
    for a in range(3):
        assert (items[3 + a].cx, items[3 + a].cy) == pytest.approx((-items[a].cx, items[a].cy), abs=1e-12)
    assert 1.05 * _area(alone) < _area(items) <= layout._SEP_MIRROR_ROOM * _area(alone)


def test_symmetry_can_be_switched_off(tmp_path):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(DATA / "2HHB.cif"), "-o", str(out), "--symmetry", "off"]) == 0


def test_a_shared_sheet_between_two_copies_yields_and_the_copies_stay_alike():
    """Plan 12 fix round 2 (1SLV): a dimer's continuous sheet is one item holding both protomers' strands, with no
    copy of its own. Squeezed between elements whose copies were tied, every push of it moved an orbit (an element
    and its copy together), each against its own neighbours: at stack pad 1.45 the tied pushes never cleared, the
    items were moved one by one and the copies drawn unalike (A15 1.28; 0.0 at 1.35, 1.10 at 1.55). Now such a
    `shared` item yields to an orbit the whole way (the orbit stays put), and slides clear last (_slide_clear), so
    the copies keep their tie; 1SLV draws alike (0.0) at pads 1.35, 1.45 and 1.55."""
    import numpy as np

    from foldmap import layout

    items = [_box(0.0, 0.0, 1, tall=3.0), _box(0.6, 0.2, 30, tall=4.0), _box(5.0, 0.0, 1, "B", tall=3.0)]
    halves = np.array([it.half for it in items])
    mids = np.array([it.mid for it in items])
    xy = np.array([(it.cx, it.cy) for it in items])
    layout._push_apart(xy, halves, mids, xy.copy(), {2: (0, 1)}, 1, shared={1})
    assert xy[0].tolist() == [0.0, 0.0] and xy[2].tolist() == [5.0, 0.0]  # the orbit kept its place
    assert xy[1][0] > 1.7  # the shared item went the whole way (1.6 to clear, less the first sweeps' pull home)
    xy = np.array([(it.cx, it.cy) for it in items])
    layout._push_apart(xy, halves, mids, xy.copy(), {2: (0, 1)}, 1)  # without: the orbit moves
    assert xy[0].tolist() != [0.0, 0.0] and xy[2].tolist() != [5.0, 0.0]
    a, b = _box(0.0, 0.0, 1), _box(0.3, 0.0, 30)  # left over: the shared item slides, not the tied one
    layout._slide_clear([b, a], last=[0])
    assert (a.cx, a.cy) == (0.0, 0.0) and (b.cx, b.cy) != (0.3, 0.0)


def test_copies_that_need_more_room_are_kept_alike_when_free_items_yield():
    """Random-100 C2b-3 cause A: a free element (u, no copies) wedged between protomer 1's copies, where protomer 0
    has room, makes the tied pushes go round in circles (test_copies_move_alone_when_moving_them_together_cannot_clear).
    Run again with the free item yielding to the orbits, they clear, and the copies stay alike: the free item moves."""
    from foldmap import layout

    a, b = (0.8, 0.3), (3.4, 3.4)
    items = [_box(*a, 1), _box(*b, 20), _box(a[0] + 3, a[1], 1, "B"), _box(b[0] + 3, b[1], 20, "B"), _box(5.3, 2.8, 40)]
    assert layout._separate(items, {2: (0, 1), 3: (1, 1)}, again=True) == "wide"
    assert _overlaps(items) == 0
    for i, j in ((2, 0), (3, 1)):
        assert (items[i].cx - items[j].cx, items[i].cy - items[j].cy) == pytest.approx(
            (items[2].cx - items[0].cx, items[2].cy - items[0].cy), abs=1e-12
        )
