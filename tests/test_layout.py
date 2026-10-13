from pathlib import Path

import numpy as np
import pytest

from foldmap.frame import Frame, view_frame
from foldmap.layout import HELIX_RISE, HELIX_W, PITCH, SCALE, STRAND_RISE, build_layout
from foldmap.style import Style
from helpers import FID, fake_sse, honest_crossings, pipeline

DATA = Path(__file__).parent / "data"
ON_SHEET = Style(helix_on_sheet="auto")  # helices packed on a sheet drawn over it (opt-in, Task 14)

MODES = ["projected", "stack"]
NAMES = ["ubq", "zs5", "zya"]


def overlaps(a, b, tol=1e-6):
    ax0, ay0, ax1, ay1 = a.rect
    bx0, by0, bx1, by1 = b.rect
    return min(ax1, bx1) - max(ax0, bx0) > tol and min(ay1, by1) - max(ay0, by0) > tol


def build(path, mode, **style):
    from foldmap.style import Style

    _, sses, sheets = pipeline(path)
    f = view_frame(sses)
    return sses, sheets, f, build_layout(sses, sheets, f, mode, style=Style(**style))


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("name", NAMES)
def test_every_sse_placed_once_without_overlap(name, mode, request):
    sses, _, _, lay = build(request.getfixturevalue(name), mode)
    assert sorted(lay.placed) == sorted(s.id for s in sses)
    items = list(lay.placed.values())
    for i, a in enumerate(items):
        for b in items[i + 1 :]:
            assert not overlaps(a, b), (a.sse.id, b.sse.id)


@pytest.mark.parametrize("mode", MODES)
def test_sheet_strands_sit_at_pitch_with_matching_arrow_directions(ubq, mode):
    _, sheets, _, lay = build(ubq, mode)
    (sheet,) = sheets
    pl = [lay.placed[s.id] for s in sheet.strands]
    xs = [p.cx for p in pl]
    steps = np.diff(sorted(xs))
    assert np.allclose(steps, PITCH)  # side by side, one pitch apart
    assert len({round(p.cy, 6) for p in pl}) == 1
    for k, kind in enumerate(sheet.pair_kinds):
        d = np.dot([np.cos(pl[k].angle), np.sin(pl[k].angle)], [np.cos(pl[k + 1].angle), np.sin(pl[k + 1].angle)])
        assert (d < 0) == (kind == "A")  # antiparallel neighbours point opposite ways


@pytest.mark.parametrize("pdb", sorted(p.stem for p in (Path(__file__).parent / "data").glob("*.cif")))
@pytest.mark.parametrize("mode", MODES)
def test_sheets_are_drawn_face_on_never_mirrored(pdb, mode):
    """Snapping strands vertical is an in-plane turn of the projected sheet: the strand drawn rightmost
    must lie to the right of the drawn-up direction in the projection, else the figure shows the back face."""
    sses, sheets, f, lay = build(Path(__file__).parent / "data" / f"{pdb}.cif", mode)
    checked = 0
    for sh in sheets:
        if len(sh.strands) < 2 or sh.ambiguous:
            continue
        drawn = sorted(sh.strands, key=lambda s: lay.placed[s.id].cx)
        first, last = drawn[0], drawn[-1]
        up3d = sum(lay.placed[x.id].direction[1] * x.axis for x in sh.strands)  # 3D direction the sheet is drawn up
        up2d = f.project(f.origin + up3d)[0]
        lateral = f.project(last.centroid)[0] - f.project(first.centroid)[0]
        right = np.array([up2d[1], -up2d[0]])
        if np.linalg.norm(up2d) < 0.3 or abs(lateral @ right) < 2.0:
            continue  # strands seen end-on, or no clear sideways spread: either order is honest
        assert lateral @ right > 0, (sh.strands[0].id, round(float(lateral @ right), 2))
        checked += 1
    if pdb in ("1UBQ", "6ZS5"):
        assert checked > 0


@pytest.mark.parametrize("name", NAMES)
def test_arrow_up_down_matches_projected_axis(name, request):
    sses, sheets, f, lay = build(request.getfixturevalue(name), "projected", helix_angle="tilted")
    checked = 0
    ambiguous = {s.id for sh in sheets if sh.ambiguous for s in sh.strands}  # one "up" cannot fit a branched sheet
    for sh in sheets:  # nor neighbours twisted far apart (6ZYA A:369-371 and A:377-383, cos -0.34)
        for a, b in zip(sh.strands, sh.strands[1:]):
            if abs(float(a.axis @ b.axis)) < 0.5:
                ambiguous |= {a.id, b.id}
    for s in sses:
        up = float(s.axis @ f.v)
        p = lay.placed[s.id]
        if s.kind == "E" and abs(up) > 0.6 and s.id not in ambiguous:
            assert (np.sin(p.angle) > 0) == (up > 0), s.id
            checked += 1
        if s.kind == "H" and np.hypot(s.axis @ f.u, up) > 0.5:
            want = np.array([s.axis @ f.u, up])
            got = np.array([np.cos(p.angle), np.sin(p.angle)])
            assert np.dot(want / np.linalg.norm(want), got) > 0.999
    if name == "ubq":
        assert checked > 0  # the check must actually have run


@pytest.mark.parametrize("mode", MODES)
def test_ports_are_the_element_ends(ubq, mode):
    sses, _, _, lay = build(ubq, mode)
    for p in lay.placed.values():
        n, c = np.array(p.n_port), np.array(p.c_port)
        assert np.linalg.norm(c - n) == pytest.approx(p.length)
        assert np.allclose((c - n) / p.length, [np.cos(p.angle), np.sin(p.angle)])
        assert np.allclose(p.exit_c, [np.cos(p.angle), np.sin(p.angle)])
        assert np.allclose(p.exit_n, -np.array(p.exit_c))


@pytest.mark.parametrize("mode", MODES)
def test_deterministic(zs5, mode):
    _, _, _, a = build(zs5, mode)
    _, _, _, b = build(zs5, mode)
    assert {k: (v.cx, v.cy, v.angle) for k, v in a.placed.items()} == {
        k: (v.cx, v.cy, v.angle) for k, v in b.placed.items()
    }


def test_labels_are_per_chain_in_sequence_order(zs5):
    sses, _, _, lay = build(zs5, "projected")
    for chain in ("A", "D"):
        strands = [s for s in sses if s.chain == chain and s.kind == "E"]
        labels = [lay.placed[s.id].label for s in strands]
        assert labels[:3] == ["A", "B", "C"]
        assert len(set(labels)) == len(labels)


def test_stack_mode_is_compact_rows_in_sequence_order(zs5):
    sses, sheets, _, lay = build(zs5, "stack")
    x0, _, x1, _ = lay.bounds
    assert x1 - x0 <= 45
    blocks = []  # (anchor chain, first residue, top, bottom, x) per sheet or helix
    for sh in sheets:
        first = min(sh.strands, key=lambda s: s.start)
        ps = [lay.placed[s.id] for s in sh.strands]
        top, bottom = max(p.rect[3] for p in ps), min(p.rect[1] for p in ps)
        blocks.append((first.chain, first.start, top, bottom, float(np.mean([p.cx for p in ps]))))
    for s in sses:
        if s.kind == "H":
            p = lay.placed[s.id]
            blocks.append((s.chain, s.start, p.rect[3], p.rect[1], p.cx))
    for chain in {b[0] for b in blocks}:
        mine = sorted(b for b in blocks if b[0] == chain)  # by first residue
        for (_, _, t0, b0, x0), (_, _, t1, b1, x1) in zip(mine, mine[1:]):
            if b1 < t0 and b0 < t1:  # side by side in one row (blocks are centred, not top-aligned): left to right
                assert x1 > x0, chain
            else:  # or on the next row down
                assert t1 <= b0, chain
    assert max(b[3] for b in blocks) >= min(b[2] for b in blocks)  # it really wraps into several rows


def test_helix_pile_on_one_point_is_separated():
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    sses = [fake_sse("H", 10 * k, 10 * k + 8, centroid=(0, 0, 0), axis=(0.2, 0.1, 1)) for k in range(12)]
    for mode in MODES:
        lay = build_layout(sses, [], f, mode)
        items = list(lay.placed.values())
        assert len(items) == 12
        assert not any(overlaps(a, b) for i, a in enumerate(items) for b in items[i + 1 :])


def test_no_elements_gives_empty_layout():
    f = view_frame([])
    lay = build_layout([], [], f, "projected")
    assert lay.placed == {} and lay.bounds == (0.0, 0.0, 0.0, 0.0)


@pytest.mark.parametrize("name", NAMES + ["ubq"])
def test_element_lengths_are_true_to_scale(name, request):
    sses, _, f, lay = build(request.getfixturevalue(name), "projected")
    for s in sses:
        p = lay.placed[s.id]
        if s.kind == "E":
            assert p.length == pytest.approx(max(1.0, (len(s) - 1) * STRAND_RISE * SCALE))
        else:  # never foreshortened: a tilted helix keeps its true length, so every coil has the same pitch
            assert p.length == pytest.approx(max(HELIX_W, (len(s) - 1) * HELIX_RISE * SCALE))


def test_one_unit_is_the_same_distance_for_strands_and_helices():
    assert PITCH == pytest.approx(4.8 * SCALE, rel=0.05)  # neighbouring strands are ~4.8 A apart
    assert HELIX_W == pytest.approx(5.0 * SCALE, rel=0.15)  # a helix ribbon is ~5 A across (CA radius 2.3 A)
    from foldmap.layout import HEAD

    assert 0.9 * HEAD <= HELIX_W <= 1.4 * HEAD  # about as wide as a strand's arrowhead, as in a cartoon


def loop_length(lay, sses, sheets=None):
    """Port-to-port distance summed over loops; with `sheets`, only loops between separately placed
    items (a loop inside one sheet block is fixed by the strand order, no translation changes it)."""
    block = {s.id: k for k, sh in enumerate(sheets or []) for s in sh.strands}
    return sum(
        np.hypot(*(np.array(lay.placed[a.id].c_port) - np.array(lay.placed[b.id].n_port)))
        for a, b in zip(sses, sses[1:])
        if a.chain == b.chain and not (sheets is not None and a.id in block and block.get(a.id) == block.get(b.id))
    )


def order_kept(home, got, axis, gap=10.0, skip=None):
    """Fraction of element pairs more than `gap` Å apart along `axis` in 3D whose page order is kept.
    `skip` (bool per element) leaves out pairs involving those elements."""
    dh = home[:, None, axis] - home[None, :, axis]
    dg = got[:, None, axis] - got[None, :, axis]
    clear = dh > gap
    if skip is not None:
        clear &= ~skip[:, None] & ~skip[None, :]
    if clear.sum() < 3:
        return 1.0  # too few clearly separated pairs to say anything about order
    return float((dg[clear] > 0).mean())


@pytest.mark.parametrize("pdb", ["6ZS5", "1TIM", "1C3W", "2HHB", "1EMA", "2LZM", "8UUP"])
def test_attraction_shortens_loops_and_keeps_the_projected_fold(pdb, ubq):
    from foldmap.features import sse_contacts

    bb, sses, sheets = pipeline(ubq.parent / f"{pdb}.cif")
    f = view_frame(sses)
    base = build_layout(sses, sheets, f, "projected", attract=False)
    new = build_layout(sses, sheets, f, "projected", contacts=sse_contacts(bb, sses))
    gain = 1.0 if any(sh.closed for sh in sheets) else 0.95  # a barrel's floor is ~0.92 even unconstrained
    assert loop_length(new, sses, sheets) < gain * loop_length(base, sses, sheets)
    assert loop_length(new, sses) <= loop_length(base, sses)
    home = f.project(np.array([s.centroid for s in sses]))
    xy = lambda lay: np.array([[lay.placed[s.id].cx, lay.placed[s.id].cy] for s in sses])  # noqa: E731
    barrel = np.array([any(s in sh.strands for sh in sheets if sh.closed) for s in sses])  # unrolled: no 3D order
    for axis in (0, 1):  # the fold's left-right/up-down order survives
        was = order_kept(home, xy(base), axis, skip=barrel)
        assert order_kept(home, xy(new), axis, skip=barrel) >= min(0.8, was) - 0.15
    assert not any(
        overlaps(a, b) for i, a in enumerate(new.placed.values()) for b in list(new.placed.values())[i + 1 :]
    )


def test_attraction_stays_bounded_when_loops_pull_hard(monkeypatch):
    import foldmap.layout as L

    monkeypatch.setattr(L, "_W_SEQ", 4.0)  # strong pulls once made the relaxation diverge
    bb, sses, sheets = pipeline(Path(__file__).parent / "data" / "1TIM.cif")
    f = view_frame(sses)
    lay = build_layout(sses, sheets, f, "projected")
    home = f.project(np.array([s.centroid for s in sses])) * SCALE
    span = np.ptp(home, axis=0).max()
    x0, y0, x1, y1 = lay.bounds
    assert max(x1 - x0, y1 - y0) < 3 * span + 20


@pytest.mark.parametrize("mode", MODES)
def test_closed_barrel_repeats_its_first_strand_after_the_last(mode):
    sses, sheets, _, lay = build(Path(__file__).parent / "data" / "1EMA.cif", mode)
    (sheet,) = sheets
    assert sheet.closed
    (ghost,) = lay.ghosts
    pl = sorted((lay.placed[s.id] for s in sheet.strands), key=lambda p: p.cx)
    first, last = pl[0], pl[-1]
    assert ghost.sse is first.sse and ghost.label == first.label and ghost.angle == first.angle
    assert ghost.cx == pytest.approx(last.cx + PITCH) and ghost.cy == pytest.approx(last.cy)
    assert all(p is not ghost for p in lay.placed.values())  # placed keeps only the real strand
    others = list(lay.placed.values())
    assert not any(overlaps(ghost, o) for o in others if o is not first)
    x0, y0, x1, y1 = lay.bounds
    assert ghost.rect[2] <= x1 + 1e-9


def test_open_sheets_have_no_ghosts(ubq):
    _, _, _, lay = build(ubq, "projected")
    assert lay.ghosts == []


@pytest.mark.parametrize("pdb", ["2HHB", "1TIM", "1BL8"])
def test_stack_mode_draws_helices_full_length_and_upright(pdb):
    sses, _, f, lay = build(Path(__file__).parent / "data" / f"{pdb}.cif", "stack")
    for s in sses:
        if s.kind == "H":
            p = lay.placed[s.id]
            assert p.length == pytest.approx(max(HELIX_W, (len(s) - 1) * HELIX_RISE * SCALE)), s.id  # no 3D tilt here
            assert abs(np.cos(p.angle)) < 1e-9 and (np.sin(p.angle) > 0) == (s.axis @ f.v >= 0), s.id


@pytest.mark.parametrize("pdb", ["2HHB", "1TIM", "8UUP", "6ZS5"])
def test_stack_rows_leave_room_for_chain_end_labels(pdb):
    from foldmap.layout import label_boxes

    sses, _, _, lay = build(Path(__file__).parent / "data" / f"{pdb}.cif", "stack")
    ends = [b for name, b in label_boxes(lay, sses) if name.startswith(("N:", "C:"))]
    for i, a in enumerate(ends):
        for b in ends[i + 1 :]:
            assert min(a[2], b[2]) <= max(a[0], b[0]) or min(a[3], b[3]) <= max(a[1], b[1]), (a, b)


def test_overlap_fallback_slide_runs_and_separates(monkeypatch):
    import foldmap.layout as L

    monkeypatch.setattr(L, "_SEP_SWEEPS", 0)  # no push-apart sweeps: straight to the guaranteed slide
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    sses = [fake_sse("H", 10 * k, 10 * k + 8, centroid=(0, 0, 0), axis=(1, 0, 0)) for k in range(5)]
    lay = build_layout(sses, [], f, "projected", attract=False)
    items = list(lay.placed.values())
    assert not any(overlaps(a, b) for i, a in enumerate(items) for b in items[i + 1 :])


# 1LMB without its symmetry: under auto its dyad lies in the page and the two α5 copies are drawn as mirror images
# of each other instead, at one height (test_symmetry_scenarios::test_dimer_on_dna_is_an_exact_mirror_image). 8UUP's
# stalk of three protomers is drawn with each copy under its symmetry (test_copies::test_a_ring_stalk_...); without it
# the three chains' helices stand side by side.
@pytest.mark.parametrize(
    "pdb, member, kw", [("8UUP", "D:456-472", {"symmetry": "off"}), ("1LMB", "3:78-88", {"symmetry": "off"})]
)
def test_bundled_helices_are_drawn_side_by_side(pdb, member, kw):
    from foldmap.cli import make_layout
    from foldmap.features import helix_bundles

    lay, sses, _ = make_layout(Path(__file__).parent / "data" / f"{pdb}.cif", **kw)
    group = next(g for g in helix_bundles(sses) if member in g)
    pl = [lay.placed[k] for k in group]
    d = pl[0].direction
    assert all(abs(abs(p.direction @ d) - 1) < 1e-9 for p in pl)  # one axis for the whole bundle
    side = np.array([-d[1], d[0]])
    columns = sorted({round(float(np.array([p.cx, p.cy]) @ side), 3) for p in pl})
    assert len(columns) == len({k.split(":")[0] for k in group})  # one column per chain
    from foldmap.features import BUNDLE_DIST

    gap = columns[1] - columns[0]
    assert np.allclose(np.diff(columns), gap) and gap <= max(BUNDLE_DIST * SCALE, pl[0].width + 1.2) + 1e-6  # packed
    span: dict[str, list[float]] = {}  # per chain column: from the lowest piece start to the highest piece end
    for k, p in zip(group, pl):
        c = float(np.array([p.cx, p.cy]) @ d)
        lo, hi = span.get(k.split(":")[0], [np.inf, -np.inf])
        span[k.split(":")[0]] = [min(lo, c - p.length / 2), max(hi, c + p.length / 2)]
    assert max(a for a, _ in span.values()) < min(b for _, b in span.values())  # columns side by side: a stalk


def test_group_by_chain_keeps_chain_order_and_filters():
    from foldmap.model import group_by_chain
    from helpers import fake_sse

    b1, a1, b2 = fake_sse("E", 3, 7, chain="B"), fake_sse("H", 1, 9, chain="A"), fake_sse("E", 12, 15, chain="B")
    assert list(group_by_chain([b1, a1, b2])) == ["B", "A"]  # first appearance, not sorted
    assert group_by_chain([b1, a1, b2])["B"] == [b1, b2]
    assert group_by_chain([b1, a1, b2], keep={b2.id}) == {"B": [b2]}


def test_centred_box_and_overlaps_any():
    from foldmap.layout import centred_box, overlaps_any

    assert centred_box(1.0, 2.0, 0.5, 0.25) == (0.5, 1.75, 1.5, 2.25)
    box = (0.0, 0.0, 1.0, 1.0)
    assert overlaps_any(box, [(5, 5, 6, 6), (0.5, 0.5, 2.0, 2.0)])
    assert not overlaps_any(box, [(1.0, 0.0, 2.0, 1.0), (0.0, 1.0, 1.0, 2.0)])  # edges only touch
    assert not overlaps_any(box, [])


def test_kink_split_pieces_are_held_end_to_end():
    """2HHB α: helices A (4-17) and B (21-35) are one helix run split at its 3₁₀ corner (Task 2). They were drawn
    5 helix widths apart (α E, seen between them, wedged itself in even against a strong joint spring), so the AB
    corner did not read as one. The pieces are now one unit, B starting where the loop from A turns once into it."""
    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(Path(__file__).parent / "data" / "2HHB.cif")
    assert [s.id for s in sses if s.kinked] == ["A:21-35", "C:21-35"]  # each piece that continues a kinked helix
    loops = {(lp.a_id, lp.b_id): lp for lp in route_loops(lay, sses, bb)}
    for c in "AC":
        a, b = lay.placed[f"{c}:4-17"], lay.placed[f"{c}:21-35"]
        assert np.hypot(*np.subtract(a.c_port, b.n_port)) <= 1.5 * a.width, c  # 4.96 widths at 2d08b75
        assert not overlaps(a, b), c
        corner = loops[(a.sse.id, b.sse.id)]
        assert not corner.fallback and len(corner.points) <= 3, c  # one turn


@pytest.mark.parametrize("turn", [0.0, 90.0, 180.0])  # running on, across, turned back
def test_kinked_pieces_never_overlap_and_meet_end_to_end(turn):
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    t = np.radians(turn)
    a = fake_sse("H", 1, 14, centroid=(0, 0, 0), axis=(1, 0, 0))
    b = fake_sse("H", 18, 31, centroid=(30, 5, 0), axis=(np.cos(t), np.sin(t), 0))
    b.kinked = True
    lay = build_layout([a, b], [], f, "projected")
    pa, pb = lay.placed[a.id], lay.placed[b.id]
    assert not overlaps(pa, pb)
    assert np.hypot(*np.subtract(pa.c_port, pb.n_port)) <= 2.0 * pa.width  # held together, not where projected


def _body(p):
    """The four corners of an element's body (length x width, turned by its angle)."""
    d = np.array([np.cos(p.angle), np.sin(p.angle)])
    n = np.array([-d[1], d[0]])
    c = np.array([p.cx, p.cy])
    return [c + sx * d * p.length / 2 + sy * n * p.width / 2 for sx, sy in ((1, 1), (1, -1), (-1, -1), (-1, 1))]


def _bodies_overlap(p, q) -> bool:
    """Separating axes: do two turned bodies overlap (touching excluded)?"""
    P, Q = _body(p), _body(q)
    for corners in (P, Q):
        for i in range(4):
            e = corners[(i + 1) % 4] - corners[i]
            axis = np.array([-e[1], e[0]])
            a, b = [v @ axis for v in P], [v @ axis for v in Q]
            if max(a) <= min(b) + 1e-6 or max(b) <= min(a) + 1e-6:
                return False
    return True


@pytest.mark.parametrize("mode", ["snap", "tilted"])
@pytest.mark.parametrize("turn", [105.0, 115.0, 125.0, 150.0])
@pytest.mark.parametrize("side", [1, -1])
def test_a_kinked_helix_drawn_tilted_never_cuts_into_itself(mode, turn, side):
    """Drawn as seen (tilted), a piece turned by 114-118° started 'across' and cut into the last piece's end; from
    120° the hairpin side came from the projected centres, not the piece's own direction (Task 12 review)."""
    from foldmap.style import Style

    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    t = np.radians(turn)
    a = fake_sse("H", 1, 14, centroid=(0, 0, 0), axis=(1, 0, 0))
    b = fake_sse("H", 18, 31, centroid=(5, -side * 3, 0), axis=(np.cos(t), side * np.sin(t), 0))  # seen on the far side
    b.kinked = True
    lay = build_layout([a, b], [], f, "projected", style=Style(helix_angle=mode))
    pa, pb = lay.placed[a.id], lay.placed[b.id]
    assert not _bodies_overlap(pa, pb)
    assert np.hypot(*np.subtract(pa.c_port, pb.n_port)) <= 2.5 * pa.width


def test_kinked_pieces_join_only_across_a_corner_of_their_own_chain():
    """A piece marked kinked joins the helix piece of its own chain that ends just before it (a 3₁₀/π corner of at
    most _KINK_CORNER residues), never merely the previous element in the list: a focus subset can leave its real
    predecessor out (Task 12 review)."""
    from foldmap.layout import _KINK_CORNER, _kink_groups

    a = fake_sse("H", 1, 14)
    b = fake_sse("H", 18, 31)
    b.kinked = True
    far = fake_sse("H", 18 + 14 + _KINK_CORNER + 2, 60)
    far.kinked = True
    other = fake_sse("H", 33, 45, chain="B")
    other.kinked = True
    assert _kink_groups([a, b], set()) == [[a, b]]
    assert _kink_groups([a, far], set()) == []  # its predecessor left out: the helix before is too far
    assert _kink_groups([b, other], set()) == []  # another chain
    assert _kink_groups([a, other, b], set()) == [[a, b]]  # a chain's pieces, whatever lies between in the list


@pytest.mark.parametrize("path, hits", [("1UBQ.cif", ["A:23-34"]), ("2LZM.cif", ["A:39-50"])])
def test_helices_on_sheets_finds_the_helix_packed_on_a_sheet_face(path, hits):
    """A helix of >= 10 residues lying flat 6-14 Å over a sheet of >= 3 strands of its chain, its centroid over the
    sheet (plan 11, measured: 1UBQ α1, 2LZM A:39-50); every other helix stays beside."""
    from foldmap.features import helices_on_sheets

    bb, sses, sheets = pipeline(DATA / path, "auto")
    over = helices_on_sheets(bb, sses, sheets)
    assert sorted(over) == hits
    assert all(len(sheets[k].strands) >= 3 for k in over.values())


def test_helices_on_sheets_takes_both_ferredoxin_helices():
    from foldmap.features import helices_on_sheets

    bb, sses, sheets = pipeline(FID / "2ACY.cif.gz", "auto")
    assert sorted(helices_on_sheets(bb, sses, sheets)) == ["A:22-32", "A:55-67"]


@pytest.mark.parametrize("path, helices", [("1UBQ.cif", ["A:23-34"]), ("fidelity/2ACY.cif.gz", ["A:22-32", "A:55-67"])])
def test_helix_packing_on_a_sheet_is_drawn_over_it(path, helices):
    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(DATA / path, look=ON_SHEET)
    block = next(b for b in lay.sheets if len(b) >= 4)
    rects = np.array([lay.placed[k].rect for k in block])
    x0, y0, x1, y1 = rects[:, 0].min(), rects[:, 1].min(), rects[:, 2].max(), rects[:, 3].max()
    for h in helices:
        p = lay.placed[h]
        assert p.over and x0 <= p.cx <= x1 and y0 <= p.cy <= y1, h
    loops = route_loops(lay, sses, bb)
    assert not [lp for lp in loops if lp.fallback]


def test_ferredoxin_helices_over_the_sheet_keep_crossings_low():
    """2ACY with both helices over its sheet: at most the 2 crossings it had with them at its flanks (plan 11 Task 14
    target, measured at b43695e), counted honestly (loops under a helix not their own, stubs across a strand not
    under their port), in under 0.8 of the area."""
    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(FID / "2ACY.cif.gz", look=ON_SHEET)
    off = make_layout(FID / "2ACY.cif.gz")[0]
    assert honest_crossings(lay, route_loops(lay, sses, bb)) <= 2
    assert lay.width * lay.height < 0.8 * off.width * off.height


def test_honest_crossings_count_loops_under_a_helix_and_stubs_off_their_port():
    from foldmap.layout import Layout, Placed
    from foldmap.route import Loop

    strand = Placed(fake_sse("E", 1, 6), 0.0, 0.0, 4.0, np.pi / 2, 1.1)
    other = Placed(fake_sse("E", 20, 25), 1.1, 0.0, 4.0, -np.pi / 2, 1.1)
    helix = Placed(fake_sse("H", 8, 18), -0.5, 0.5, 2.0, 0.0, 1.0, over=True)  # its C port on `strand`
    far = Placed(fake_sse("H", 30, 40), 5.0, 0.0, 2.0, np.pi / 2, 1.0)
    lay = Layout({p.sse.id: p for p in (strand, other, helix, far)}, (-3, -3, 7, 3))
    lay.on_sheet = {helix.sse.id: [strand.sse.id, other.sse.id]}
    stub = Loop(helix.sse.id, far.sse.id, [helix.c_port, (1.6, 0.5), (1.6, 2.5), (5.0, 2.5), far.n_port], False, False)
    assert honest_crossings(lay, [stub]) == 1  # its stub crosses `other`, which its port does not lie on
    under = Loop(strand.sse.id, far.sse.id, [strand.c_port, (0.0, 3.0), (5.0, 3.0), far.n_port], False, False)
    beneath = Loop(far.sse.id, other.sse.id, [(-1.0, 0.5), (0.5, 0.5)], False, False)
    assert honest_crossings(lay, [under]) == 0 and honest_crossings(lay, [beneath]) == 2  # under the helix, `strand`


def test_symmetric_and_template_sheets_have_no_helix_over_them():
    from foldmap.cli import make_layout

    for path in ("1TIM.cif", "8UUP.cif"):
        lay, _, _ = make_layout(DATA / path, look=ON_SHEET)
        assert not any(p.over for p in lay.placed.values()), path


def test_helix_on_sheet_is_opt_in():
    """Off by default (Task 14 ruling): a helix packed on a sheet stays beside it unless helix_on_sheet=auto."""
    from foldmap.cli import make_layout
    from foldmap.style import CHOICES, Style

    assert Style().helix_on_sheet == "off" and CHOICES["helix_on_sheet"] == ("off", "auto")
    assert not make_layout(DATA / "1UBQ.cif")[0].placed["A:23-34"].over
    lay, _, _ = make_layout(DATA / "1UBQ.cif", look=Style(helix_on_sheet="off"))
    assert not lay.placed["A:23-34"].over


def test_helix_over_a_sheet_lies_as_it_packs():
    """1UBQ α1 crosses the sheet diagonally, as Burroughs 2007 draws it (Fig 1B): neither upright nor flat, its
    offsets and angle snapped; the sheet with it is narrower than the sheet with α1 at its side."""
    from foldmap.cli import make_layout
    from foldmap.layout import _snap
    from foldmap.style import Style

    lay, _, _ = make_layout(DATA / "1UBQ.cif", look=ON_SHEET)
    p = lay.placed["A:23-34"]
    tilt = abs(np.degrees(p.angle)) % 180
    assert 15 < min(tilt, 180 - tilt) and abs(tilt - 90) > 15, tilt
    assert p.angle == float(_snap(p.angle))
    off = make_layout(DATA / "1UBQ.cif", look=Style(helix_on_sheet="off"))[0]
    assert lay.width < off.width


def test_over_helices_are_dropped_where_a_layout_cannot_hold_them():
    """build_layout takes `over` (helix id -> sheet index) but never draws a helix over a βα template row, in stack
    mode, or for a piece of a kinked helix (a rigid unit, 12a, kept whole beside the sheet)."""
    from foldmap.style import Style

    bb, sses, sheets = pipeline(DATA / "1UBQ.cif", "auto")
    f = view_frame(sses)
    over = {"A:23-34": 0}
    assert build_layout(sses, sheets, f, over=over).placed["A:23-34"].over
    assert not build_layout(sses, sheets, f, "stack", over=over).placed["A:23-34"].over
    assert not build_layout(sses, sheets, f, style=Style(), over={}).placed["A:23-34"].over
    h = next(s for s in sses if s.id == "A:23-34")
    piece = fake_sse("H", h.end + 2, h.end + 4, centroid=h.centroid + 6 * h.axis, axis=h.axis)
    piece.kinked = True  # α1 split at a kink: the two pieces are one rigid unit (12a), drawn beside the sheet
    with_piece = sorted([*sses, piece], key=lambda s: (s.chain, s.start))
    lay = build_layout(with_piece, sheets, f, over=over)
    assert not lay.placed["A:23-34"].over and not lay.placed[piece.id].over


def test_no_helix_is_drawn_over_a_sheet_beside_nucleic_acid():
    """1URN's RRM α1 packs on its sheet, but the RNA's stems and free ends are placed by the protein elements they
    touch (nacheck's stems_nearer_frac and the free-end pull were measured with the helix beside): no helix over a
    sheet in a figure holding nucleic acid."""
    from foldmap.cli import make_layout
    from foldmap.features import helices_on_sheets

    bb, sses, sheets = pipeline(DATA / "na" / "1URN.cif.gz", "auto")
    assert helices_on_sheets(bb, sses, sheets)
    lay, _, _ = make_layout(DATA / "na" / "1URN.cif.gz", look=ON_SHEET)
    assert lay.stems and not lay.on_sheet and not any(p.over for p in lay.placed.values())


def test_copies_of_a_protein_draw_a_helix_over_their_sheets_alike():
    """A helix over its sheet in one copy is over it in every copy drawn (model.homolog_chains, residue numbers taken
    to the first copy's); a copy whose helix is not over its sheet keeps the others beside too."""
    from foldmap.layout import _over_alike

    h_a, h_b, h_c = (fake_sse("H", 22, 33, chain=c) for c in "ABC")
    h_c.first, h_c.last = type(h_c.first)("C", 123, "", "ALA"), type(h_c.last)("C", 134, "", "ALA")
    sses = [h_a, h_b, h_c]
    copies = {"A": ("A", 0), "B": ("A", 0), "C": ("A", 100)}
    rowed = [False, False, False]
    both = {h_a.id: 0, h_b.id: 1, h_c.id: 2}
    assert _over_alike(both, sses, rowed, set(), copies) == both
    assert _over_alike({h_a.id: 0, h_b.id: 1}, sses, rowed, set(), copies) == {}  # C's copy stays beside
    assert _over_alike({h_a.id: 0}, [h_a, h_b], rowed, set(), {}) == {h_a.id: 0}  # unrelated chains: each its own
    assert _over_alike(both, sses, [False, True, False], set(), copies) == {}  # B's sheet is a βα row
    assert _over_alike(both, sses, rowed, {h_b.id}, copies) == {}  # B's helix is held by a bundle or a kink


@pytest.mark.parametrize("path", ["1UBQ.cif", "fidelity/2ACY.cif.gz", "2LZM.cif"])
def test_the_name_of_a_helix_over_a_sheet_keeps_off_arrowheads(path):
    """Its name's plate covers no strand's arrowhead, nor the straight start of a loop out of a strand's end (Task 14
    review: 1UBQ α1's name sat on β1's head)."""
    from foldmap.cli import make_layout
    from foldmap.layout import HEAD, HEAD_LEN, _element_label_box

    lay, _, _ = make_layout(DATA / path, look=ON_SHEET)
    for h in lay.on_sheet:
        box = _element_label_box(lay.placed[h])
        for q in [*lay.placed.values(), *lay.ghosts]:
            if q.sse.kind != "E":
                continue
            tip = np.array(q.c_port)
            head = [
                tip - q.direction * HEAD_LEN + n * HEAD / 2 * q.width / PITCH
                for n in (q.direction[::-1] * (1, -1), -q.direction[::-1] * (1, -1))
            ]
            xs, ys = zip(tip, *head)
            assert not (box[0] < max(xs) and min(xs) < box[2] and box[1] < max(ys) and min(ys) < box[3]), (h, q.sse.id)


def test_layout_fields_after_bounds_are_keywords():
    """Everything after placed and bounds is passed by name: build_layout once passed a bool by position into the
    slot that is now `domains` (a field that had been removed), harmless only because domains was set again later."""
    from foldmap.layout import Layout

    with pytest.raises(TypeError):
        Layout({}, (0.0, 0.0, 0.0, 0.0), [], [], None, [], True)
    bb, sses, sheets = pipeline(DATA / "2LZM.cif")
    lay = build_layout(sses, sheets, view_frame(sses))
    assert lay.domains == [] and lay.sheets and lay.nucleic is None
