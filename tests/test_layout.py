from pathlib import Path

import numpy as np
import pytest
from helpers import fake_sse, pipeline

from topoplot.frame import Frame, view_frame
from topoplot.layout import HELIX_W, HELIX_RISE, PITCH, SCALE, STRAND_RISE, build_layout

MODES = ["projected", "stack"]
NAMES = ["ubq", "zs5", "zya"]


def overlaps(a, b, tol=1e-6):
    ax0, ay0, ax1, ay1 = a.rect
    bx0, by0, bx1, by1 = b.rect
    return min(ax1, bx1) - max(ax0, bx0) > tol and min(ay1, by1) - max(ay0, by0) > tol


def build(path, mode, **style):
    from topoplot.style import Style

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
    assert {k: (v.cx, v.cy, v.angle) for k, v in a.placed.items()} == {k: (v.cx, v.cy, v.angle) for k, v in b.placed.items()}


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
    blocks = []  # (anchor chain, first residue, row, x) per sheet or helix
    for sh in sheets:
        first = min(sh.strands, key=lambda s: s.start)
        ps = [lay.placed[s.id] for s in sh.strands]
        blocks.append((first.chain, first.start, round(-ps[0].cy, 3), float(np.mean([p.cx for p in ps]))))
    for s in sses:
        if s.kind == "H":
            p = lay.placed[s.id]
            blocks.append((s.chain, s.start, round(-p.cy, 3), p.cx))
    for chain in {b[0] for b in blocks}:
        mine = sorted(b for b in blocks if b[0] == chain)  # by first residue
        reading = [(b[2], b[3]) for b in mine]
        assert reading == sorted(reading), chain  # row by row, left to right
    assert len({b[2] for b in blocks}) > 1  # it really wraps into several rows


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
    from topoplot.layout import HEAD

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
    from topoplot.features import sse_contacts

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
    import topoplot.layout as L

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
    from topoplot.layout import label_boxes

    sses, _, _, lay = build(Path(__file__).parent / "data" / f"{pdb}.cif", "stack")
    ends = [b for name, b in label_boxes(lay, sses) if name.startswith(("N:", "C:"))]
    for i, a in enumerate(ends):
        for b in ends[i + 1 :]:
            assert min(a[2], b[2]) <= max(a[0], b[0]) or min(a[3], b[3]) <= max(a[1], b[1]), (a, b)


def test_overlap_fallback_slide_runs_and_separates(monkeypatch):
    import topoplot.layout as L

    monkeypatch.setattr(L, "_SEP_SWEEPS", 0)  # no push-apart sweeps: straight to the guaranteed slide
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    sses = [fake_sse("H", 10 * k, 10 * k + 8, centroid=(0, 0, 0), axis=(1, 0, 0)) for k in range(5)]
    lay = build_layout(sses, [], f, "projected", attract=False)
    items = list(lay.placed.values())
    assert not any(overlaps(a, b) for i, a in enumerate(items) for b in items[i + 1 :])


@pytest.mark.parametrize("pdb, member", [("8UUP", "D:456-485"), ("1LMB", "3:78-90")])
def test_bundled_helices_are_drawn_side_by_side(pdb, member):
    from topoplot.cli import make_layout
    from topoplot.features import helix_bundles

    lay, sses, _ = make_layout(Path(__file__).parent / "data" / f"{pdb}.cif")
    group = next(g for g in helix_bundles(sses) if member in g)
    pl = [lay.placed[k] for k in group]
    d = pl[0].direction
    assert all(abs(abs(p.direction @ d) - 1) < 1e-9 for p in pl)  # one axis for the whole bundle
    side = np.array([-d[1], d[0]])
    columns = sorted({round(float(np.array([p.cx, p.cy]) @ side), 3) for p in pl})
    assert len(columns) == len({k.split(":")[0] for k in group})  # one column per chain
    from topoplot.features import BUNDLE_DIST

    gap = columns[1] - columns[0]
    assert np.allclose(np.diff(columns), gap) and gap <= max(BUNDLE_DIST * SCALE, pl[0].width + 1.2) + 1e-6  # packed
    span: dict[str, list[float]] = {}  # per chain column: from the lowest piece start to the highest piece end
    for k, p in zip(group, pl):
        c = float(np.array([p.cx, p.cy]) @ d)
        lo, hi = span.get(k.split(":")[0], [np.inf, -np.inf])
        span[k.split(":")[0]] = [min(lo, c - p.length / 2), max(hi, c + p.length / 2)]
    assert max(a for a, _ in span.values()) < min(b for _, b in span.values())  # columns side by side: a stalk
