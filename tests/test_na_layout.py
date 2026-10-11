"""Nucleic acids on the page: stems placed like elements, strands routed, four drawing styles, the explorer."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_figure_and_loops, make_layout
from foldmap.layout import PlacedStem
from foldmap.route import route_loops

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB"]
SMALL = ["30QM", "1EHZ", "1IHF", "3CRX"]


def _path(pdb):
    return NA / f"{pdb}.cif.gz"


def _routed(pdb, **kw):
    lay, sses, bb = make_layout(_path(pdb), **kw)
    loops = route_loops(lay, sses, bb)
    return lay, sses, bb, loops


def _overlap(a, b):
    return min(a[2], b[2]) - max(a[0], b[0]) > 1e-6 and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6


def _gid(fig, prefix):
    ax = fig.axes[0]
    return [a for a in [*ax.patches, *ax.texts, *ax.lines, *ax.collections] if (a.get_gid() or "").startswith(prefix)]


@pytest.mark.parametrize("pdb", FIXTURES)
def test_every_stem_piece_is_placed_without_overlap(pdb):
    lay, _, _ = make_layout(_path(pdb))
    pieces = len(lay.stems)
    assert pieces >= len(lay.nucleic.structure.stems) > 0
    boxes = [p.rect for p in lay.placed.values()] + [d.rect for d in lay.stems]
    assert not any(_overlap(a, b) for k, a in enumerate(boxes) for b in boxes[k + 1 :])
    x0, y0, x1, y1 = lay.bounds
    assert all(x0 <= r[0] and r[2] <= x1 and y0 <= r[1] and r[3] <= y1 for r in boxes)


def test_stems_are_drawn_to_scale_along_their_projected_axes():
    lay, _, _ = make_layout(_path("30QM"))
    d = next(d for d in lay.stems if d.name == "P1")
    assert isinstance(d, PlacedStem)
    assert d.width == pytest.approx(20.0 * 0.23)
    firsts = sorted((key for key in d.axial if d.side[key] == d.first), key=lambda k: k[1])
    assert np.all(np.diff([d.axial[k] for k in firsts]) > 0)  # the first arm runs 5'->3' along the direction
    seconds = sorted((key for key in d.axial if d.side[key] != d.first), key=lambda k: k[1])
    assert np.all(np.diff([d.axial[k] for k in seconds]) < 0)  # the second arm back
    for a, i, b, j in d.pairs:
        assert d.axial[(a, i)] == pytest.approx(d.axial[(b, j)])  # a rung straight across


def _joints(lay):
    """(exit port, entry port, piece width) where a strand goes from one piece of a bent duplex to the next."""
    at = {key: d for d in lay.stems for key in d.axial}
    out = []
    for s, strand in enumerate(lay.nucleic.strands):
        for i in range(len(strand) - 1):
            d, e = at.get((s, i)), at.get((s, i + 1))
            if d is not None and e is not None and d is not e and d.parent and d.parent == e.parent:
                out.append((d.port((s, i), True)[0], e.port((s, i + 1), False)[0], d.width))
    return out


@pytest.mark.parametrize("pdb", ["1KX5", "1IHF"])
def test_a_bent_duplex_is_a_chain_of_pieces_end_to_end(pdb):
    """Consecutive pieces stay near each other. Box separation parts tilted pieces, so a joint can be a couple of
    duplex widths long; Task 7 routes joints as loops between stem ends (route_nucleic, _clear). This bounds the
    current lengths (1KX5 at most 2.47 widths, 1IHF 1.60) against regressions."""
    lay, _, _ = make_layout(_path(pdb))
    joints = _joints(lay)
    assert len(joints) >= (40 if pdb == "1KX5" else 4)  # both strands at every joint
    gaps = [float(np.hypot(*(q - p))) / w for p, q, w in joints]
    assert max(gaps) <= 2.6, sorted(gaps)[-5:]


@pytest.mark.parametrize("pdb", FIXTURES)
def test_both_arms_run_monotonic_along_every_stem(pdb):
    lay, _, _ = make_layout(_path(pdb))
    for d in lay.stems:
        firsts = sorted((key for key in d.axial if d.side[key] == d.first), key=lambda k: k[1])
        seconds = sorted((key for key in d.axial if d.side[key] != d.first), key=lambda k: k[1])
        assert np.all(np.diff([d.axial[k] for k in firsts]) > 0), d.name  # bulges included
        assert np.all(np.diff([d.axial[k] for k in seconds]) < 0), d.name


@pytest.mark.parametrize("pdb", FIXTURES)
def test_bulged_nucleotides_get_room_on_their_rail(pdb):
    from foldmap.layout import _MIN_RISE, SCALE

    lay, _, _ = make_layout(_path(pdb))
    for d in lay.stems:
        for sign in (1, -1):
            along = sorted(v for key, v in d.axial.items() if d.side[key] == sign)
            assert np.all(np.diff(along) >= 0.9 * _MIN_RISE * SCALE), (d.name, np.round(np.diff(along), 3))


def _upright_size(t, renderer):
    """Width and height (display units) of a text as if written level."""
    rot = t.get_rotation()
    t.set_rotation(0)
    box = t.get_window_extent(renderer)
    t.set_rotation(rot)
    return box.width, box.height


@pytest.mark.parametrize("pdb", ["1KX5", "30QM", "3CRX"])
def test_stem_names_keep_clear_of_base_letters(pdb):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig, _ = make_figure_and_loops(_path(pdb))
    renderer = FigureCanvasAgg(fig).get_renderer()
    to_px = fig.axes[0].transData.transform
    names, letters = _gid(fig, "na:label:"), _gid(fig, "na:nt:")
    assert names and letters
    for n in names:  # each letter (as a disc) against the name's own box, written along its stem
        w, h = _upright_size(n, renderer)
        c, rot = to_px(n.get_position()), np.radians(n.get_rotation())
        along, across = np.array([np.cos(rot), np.sin(rot)]), np.array([-np.sin(rot), np.cos(rot)])
        for t in letters:
            r = max(_upright_size(t, renderer)) * 1.5  # a letter and a half from its centre: read apart at a glance
            q = to_px(t.get_position()) - c
            assert abs(q @ along) >= w / 2 + r or abs(q @ across) >= h / 2 + r, (n.get_text(), t.get_gid())


def test_layout_is_the_same_whatever_the_hash_seed():
    code = (
        "from foldmap.cli import make_layout\n"
        "from foldmap.route import route_loops\n"
        f"lay, sses, bb = make_layout({str(_path('30QM'))!r})\n"
        "print([(d.id, round(d.cx, 6), round(d.cy, 6), round(d.angle, 6)) for d in lay.stems])\n"
        "route_loops(lay, sses, bb)\n"
        "print([(p.kind, p.points, p.nt_xy) for p in lay.na_paths], lay.na_ends, sorted(lay.na_letters.items()))"
    )
    runs = {
        subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONHASHSEED": seed, "OPENBLAS_NUM_THREADS": "1"},
            check=True,
        ).stdout
        for seed in ("0", "1", "2")
    }
    assert len(runs) == 1


def test_wobble_pairs_are_dotted_rungs():
    fig, _ = make_figure_and_loops(_path("30QM"))
    rungs = [a for a in _gid(fig, "na:stem:") if a.get_gid().endswith(":rung")]
    assert sum(1 for r in rungs if r.get_linestyle() != "-") == 4  # 30QM has four G·U pairs


def test_beads_mark_every_touched_nucleotide_in_its_chain_colour():
    from matplotlib.colors import to_hex

    from foldmap.palette import chain_colors

    fig, lay = make_figure_and_loops(_path("30QM"))[0], make_layout(_path("30QM"))[0]
    beads = _gid(fig, "na:contact:")
    colours = chain_colors(["A", "B"])
    for b in beads:
        chain = b.get_gid().split(":")[2]
        assert {to_hex(c) for c in b.get_facecolors()} == {colours[chain].lower()}
    touched = {x for v in lay.nucleic.contacts.values() for x in v}
    assert sum(len(b.get_offsets()) for b in beads) == len(touched)


def test_protein_free_rna_is_drawn_on_its_own():
    lay, sses, bb, loops = _routed("1EHZ")
    assert len(bb) == 0 and sses == [] and loops == []
    assert [d.name for d in lay.stems] == ["P1", "P2", "P3", "P4"]
    assert lay.na_paths


@pytest.mark.parametrize("pdb", FIXTURES)
def test_strands_route_without_fallbacks_and_every_nucleotide_is_drawn_once(pdb):
    lay, _, _, _ = _routed(pdb)
    assert not [p for p in lay.na_paths if p.fallback]
    na = lay.nucleic
    drawn = [key for d in lay.stems for key in d.axial] + [(p.strand, i) for p in lay.na_paths for i in p.nts]
    assert len(drawn) == len(set(drawn))
    assert set(drawn) == {(s, i) for s, strand in enumerate(na.strands) for i in range(len(strand))}
    for p in lay.na_paths:
        assert len(p.nt_xy) == len(p.nts)


def test_runs_go_round_the_other_elements():
    from test_route import segments_clear

    lay, _, _, _ = _routed("30QM")
    for p in lay.na_paths:
        if p.kind == "run":
            others = [q.rect for q in lay.placed.values()]
            assert segments_clear(p.points[1:-1], others), p.nts[:3]


@pytest.mark.parametrize("pdb", ["1KX5", "1IHF", "3CRX", "30QM", "4OO8"])
def test_joints_and_runs_cross_no_stem_or_element(pdb):
    """Every run and every joint between two pieces of a bent duplex, as drawn from where it leaves its own stems'
    bodies, passes the router's own final check (route._clear) against every element and stem."""
    from foldmap.route import _clear, _element_obstacles
    from foldmap.route_na import _outside

    lay, _, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    rects, owners = _element_obstacles(lay)
    kinds = {"run": 0, "joint": 0}
    for p in lay.na_paths:
        if p.kind in kinds:
            own = tuple(at[key] for key in p.between)
            assert _clear(_outside(p.points, own), own, rects, owners, own_body=True), (p.kind, p.strand, p.between)
            kinds[p.kind] += 1
    assert kinds["joint"] >= (20 if pdb == "1KX5" else 2 if pdb == "1IHF" else 0)
    assert kinds["run"] or pdb == "1KX5"


@pytest.mark.parametrize("pdb", FIXTURES)
def test_runs_keep_off_the_stem_names(pdb):
    from foldmap.layout import label_boxes
    from foldmap.route_na import _outside

    lay, sses, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    names = [b for n, b in label_boxes(lay, sses) if n.startswith("stem-label:")]
    for p in lay.na_paths:
        if p.kind in ("run", "joint"):
            pts = _outside(p.points, tuple(at[key] for key in p.between))
            assert all(segment_clear(a, b, names) for a, b in zip(pts, pts[1:])), (p.kind, p.between)


def segment_clear(p, q, rects):
    from foldmap.route import clip_len

    return all(clip_len(p, q, r) <= 1e-9 for r in rects)


def test_hairpin_loops_are_arcs_beyond_the_stem_tip():
    lay, _, _, _ = _routed("30QM")
    by_key = {key: d for d in lay.stems for key in d.axial}
    hairpins = [p for p in lay.na_paths if p.kind == "hairpin"]
    assert len(hairpins) == 5
    for p in hairpins:
        d = by_key[(p.strand, p.nts[0] - 1)]
        tip = max(d.axial[(p.strand, p.nts[0] - 1)], d.axial[(p.strand, p.nts[-1] + 1)], key=abs)
        for xy in p.nt_xy:  # every loop nucleotide beyond the stem's end
            along = (np.asarray(xy) - [d.cx, d.cy]) @ d.direction
            assert along * np.sign(tip) > abs(tip) - 1e-6


def test_every_strand_end_is_labelled_unless_cropped_into_a_gap():
    lay, _, _, _ = _routed("4YBB")
    ends = {(s, e) for s, e, _ in lay.na_ends}
    gapped = {(g.before, "3′") for g in lay.nucleic.gaps} | {(g.after, "5′") for g in lay.nucleic.gaps}
    for s in range(len(lay.nucleic.strands)):
        for e in ("5′", "3′"):
            assert ((s, e) in ends) != ((s, e) in gapped)
    assert len([p for p in lay.na_paths if p.kind == "gap"]) == len(lay.nucleic.gaps)
    assert lay.na_cropped


def test_short_gaps_are_a_bare_break_and_long_ones_say_how_many_nucleotides():
    fig, _ = make_figure_and_loops(_path("4YBB"))
    lay, _, _, _ = _routed("4YBB")
    gaps = [p for p in lay.na_paths if p.kind == "gap"]
    counts = sorted(g.count for g in lay.nucleic.gaps)
    assert counts[0] < 3 <= counts[-1]  # 4YBB cropped: 1-nt gaps from keeping stems whole, and long ones
    for p, g in zip(gaps, lay.nucleic.gaps):
        assert p.label == ("" if g.count < 3 else f"… {g.count} nt …")
    labels = sorted(t.get_text() for t in _gid(fig, "na:gap:"))
    assert labels == sorted(f"… {c} nt …" for c in counts if c >= 3)
    assert len([a for a in _gid(fig, "na:run:") if a.get_linestyle() != "-"]) == len(gaps)  # every gap is a break


@pytest.mark.parametrize("pdb", ["30QM", "4OO8", "1EHZ"])
def test_no_text_leaves_the_page_and_letters_do_not_collide(pdb):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig, _ = make_figure_and_loops(_path(pdb))
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    r = canvas.get_renderer()
    page = fig.bbox
    texts = [t for t in fig.axes[0].texts if t.get_text().strip()]
    for t in texts:
        e = t.get_window_extent(r)
        assert e.x0 >= page.x0 - 1 and e.y0 >= page.y0 - 1 and e.x1 <= page.x1 + 1 and e.y1 <= page.y1 + 1, t.get_text()
    boxes = [t.get_window_extent(r) for t in texts if (t.get_gid() or "").startswith("na:nt:")]
    clash = 0
    for k, a in enumerate(boxes):
        for b in boxes[k + 1 :]:
            w, h = min(a.x1, b.x1) - max(a.x0, b.x0), min(a.y1, b.y1) - max(a.y0, b.y0)
            clash += w > 0 and h > 0 and w * h > min(a.width * a.height, b.width * b.height) / 3
    assert clash == 0


def test_strands_without_pairs_are_traced(monkeypatch):
    import foldmap.io

    monkeypatch.setattr(foldmap.io, "pair_nucleotides", lambda runs: [])  # as for a G-quadruplex
    fig, _ = make_figure_and_loops(_path("1IHF"))
    lay, _, _, _ = _routed("1IHF")
    assert lay.stems == [] and set(lay.na_traces) == {0, 1, 2}
    assert [p.kind for p in lay.na_paths] == ["trace"] * 3
    assert len(_gid(fig, "na:nt:")) == 70 and len(_gid(fig, "na:end:")) == 6
    assert any("DNA (C, D, E)" in t.get_text() for t in fig.axes[0].texts)


def test_plot_warns_about_strands_drawn_as_fallbacks(tmp_path, capsys, monkeypatch):
    import foldmap.route_na as route_na
    from foldmap.cli import main

    real = route_na.route_nucleic

    def one_stuck(layout, *a, **k):  # every fixture routes cleanly; mark one run as a fallback
        paths = real(layout, *a, **k)
        next(p for p in paths if p.kind == "run").fallback = True
        return paths

    monkeypatch.setattr(route_na, "route_nucleic", one_stuck)
    assert main(["plot", str(_path("1EHZ")), "-o", str(tmp_path / "t.svg")]) == 0
    err = capsys.readouterr().err
    lay, _, _, _ = _routed("1EHZ")
    run = next(p for p in lay.na_paths if p.kind == "run")
    labels = lay.nucleic.strands[run.strand].labels
    name = f"strand {labels[run.nts[0]].chain}:{labels[run.nts[0]].tag}-{labels[run.nts[-1]].tag}"
    assert "no clear route" in err and name in err, err


def test_app_and_spec_report_stuck_strands(monkeypatch):
    import foldmap.route_na as route_na
    from foldmap.spec import FigureSpec, draw_spec

    real = route_na.route_nucleic

    def all_stuck(layout, *a, **k):
        paths = real(layout, *a, **k)
        for p in paths:
            p.fallback = p.kind == "hairpin"
        return paths

    monkeypatch.setattr(route_na, "route_nucleic", all_stuck)
    drawn = draw_spec(FigureSpec(_path("1EHZ")))
    assert len([s for s in drawn.stuck if s.startswith("strand A:")]) == 3  # 1EHZ's three hairpin loops


def _segments_cross(p, q, a, b) -> bool:
    """Segments p-q and a-b cross properly (not just touch at an end)."""

    def turn(u, v, w):
        return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])

    return turn(a, b, p) * turn(a, b, q) < 0 and turn(p, q, a) * turn(p, q, b) < 0


@pytest.mark.parametrize("pdb", ["1EHZ", "4OO8", "30QM"])
def test_runs_steer_round_hairpin_loops_and_free_ends(pdb):
    from foldmap.route_na import _outside

    lay, _, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    drawn = {id(p): _outside(p.points, tuple(at[key] for key in p.between)) for p in lay.na_paths}
    fixed = [p for p in lay.na_paths if p.kind in ("hairpin", "tail")]
    for p in lay.na_paths:
        if p.kind in ("run", "joint", "gap"):
            a = drawn[id(p)]
            for f in fixed:
                b = drawn[id(f)]
                assert not any(_segments_cross(*s, *t) for s in zip(a, a[1:]) for t in zip(b, b[1:])), (
                    p.kind,
                    p.nts[:2],
                    f.kind,
                    f.nts[:2],
                )


def _in_outline(pt, e) -> bool:
    c, s = np.cos(e.angle), np.sin(e.angle)
    x, y = pt[0] - e.cx, pt[1] - e.cy
    return abs(x * c + y * s) < e.length / 2 and abs(y * c - x * s) < e.width / 2


@pytest.mark.parametrize("pdb", FIXTURES)
def test_end_labels_keep_clear_of_loops_elements_and_other_labels(pdb):
    from foldmap.layout import _GLYPH_H, _GLYPH_W, centred_box, label_boxes
    from foldmap.route import clip_len

    lay, sses, _, loops = _routed(pdb)
    hw, hh = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2
    boxes = [centred_box(x, y, hw, hh) for _, _, (x, y) in lay.na_ends]
    others = [b for _, b in label_boxes(lay, sses)]
    bodies = [*lay.placed.values(), *lay.stems]
    for k, (s, end, _) in enumerate(lay.na_ends):
        box = boxes[k]
        assert not any(_overlap(box, o) for o in others + boxes[:k] + boxes[k + 1 :]), (s, end)
        for loop in loops:
            assert all(clip_len(p, q, box) <= 1e-9 for p, q in zip(loop.points, loop.points[1:])), (s, end, loop.a_id)
        samples = [(x, y) for x in np.linspace(box[0], box[2], 4) for y in np.linspace(box[1], box[3], 4)]
        assert not any(_in_outline(q, e) for q in samples for e in bodies), (s, end)


def test_letters_are_placed_quickly_on_a_whole_ribosomal_rna(monkeypatch):
    """Each letter looks only at its neighbours (a cell hash), not at every letter placed before it: on the whole
    4YBB small-subunit RNA (1000 nt, na_extent=all) placing them took 6.5 s when it was quadratic."""
    import time

    import foldmap.route_na as route_na
    from foldmap.style import Style

    lay, sses, bb = make_layout(_path("4YBB"), look=Style(na_extent="all"))
    assert lay.nucleic.size == 1000
    real, spent = route_na._place_letters, []

    def timed(*a, **k):
        t = time.perf_counter()
        out = real(*a, **k)
        spent.append(time.perf_counter() - t)
        return out

    monkeypatch.setattr(route_na, "_place_letters", timed)
    route_loops(lay, sses, bb)
    assert len(lay.na_letters) > 400 and spent[0] < 1.5, spent
