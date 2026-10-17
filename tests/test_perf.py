"""Speed-ups that must not change results: each fast path checked against the plain one it replaced."""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

import foldmap.layout as L
from foldmap import io
from foldmap.frame import Frame
from foldmap.layout import build_layout
from helpers import fake_sse

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).parent.parent

_SIGNATURE = """
import json, sys
from foldmap.cli import make_layout
out = {}
for f in sys.argv[1:]:
    lay, _, _ = make_layout(f)
    out[f] = {k: [p.cx, p.cy, p.angle] for k, p in sorted(lay.placed.items())}
print(json.dumps(out))
"""


def test_layout_does_not_depend_on_the_hash_seed():
    files = [str(DATA / f"{p}.cif") for p in ("2HHB", "1LMB", "9TX7")]
    runs = []
    for seed in ("0", "12345", "777"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run(
            [sys.executable, "-c", _SIGNATURE, *files], env=env, capture_output=True, text=True, check=True
        )
        runs.append(json.loads(done.stdout))
    assert runs[0] and all(runs[0][f] for f in runs[0])
    assert runs[0] == runs[1] == runs[2]


def _pile(n):
    """n helices projected onto one point: the push-apart must untangle them by itself."""
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    return [fake_sse("H", 10 * k, 10 * k + 8, centroid=(0, 0, 0), axis=(0.2, 0.1, 1)) for k in range(n)], f


def test_a_figure_reads_its_file_once(monkeypatch):
    from foldmap import cli

    calls = []
    real = io.read_model

    def counting(*a, **k):
        calls.append(a[0])
        return real(*a, **k)

    monkeypatch.setattr(io, "read_model", counting)
    monkeypatch.setattr(cli, "read_model", counting, raising=False)
    cli.make_layout(DATA / "1LMB.cif")  # protein, DNA and links all come from the one read
    assert len(calls) == 1
    cli.summarize(str(DATA / "1LMB.cif"))
    assert len(calls) == 2


def _tiled(pdb, copies):
    """Model 1 of a test entry, repeated `copies` times 80 Å apart along x (new chain names)."""
    import gemmi

    st = gemmi.read_structure(str(DATA / f"{pdb}.cif"))
    st.setup_entities()
    model = st[0]
    base = [c.clone() for c in model]
    for c in list(model):
        model.remove_chain(c.name)
    names = (a + b for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" for b in "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    for k in range(copies):
        for c in base:
            new = c.clone()
            new.name = next(names)
            for r in new:
                for a in r:
                    a.pos = gemmi.Position(a.pos.x + 80 * k, a.pos.y, a.pos.z)
            model.add_chain(new)
    return model


@pytest.mark.parametrize("pdb, copies", [("1C3W", 1), ("2HHB", 1), ("1BL8", 1), ("25IN", 1), ("25IN", 6)])
def test_bonded_ligand_parts_equal_the_all_against_all_search(pdb, copies):
    parts = io._ligand_parts(_tiled(pdb, copies))
    assert parts
    brute = {
        (a, b)
        for a in range(len(parts))
        for b in range(a + 1, len(parts))
        if not parts[a][3]
        and not parts[b][3]
        and any(x.pos.dist(y.pos) < io._BOND for x in parts[a][2] for y in parts[b][2])
    }
    assert io._bonded_parts(parts) == brute


def test_no_ligands_need_no_search():
    assert io._bonded_parts([]) == set()


def test_pairing_a_large_rna_is_fast():
    import time

    from foldmap.io import _nucleic_strands, read_model
    from foldmap.nucleic import pair_nucleotides

    model, _ = read_model(Path(__file__).parent / "data" / "na" / "4YBB.cif.gz", "asu")
    runs = [run for _, run in _nucleic_strands(model)]
    t = time.perf_counter()
    pairs = pair_nucleotides(runs)
    assert len(pairs) > 250 and time.perf_counter() - t < 3.0  # 1,000 nt: about 0.2 s measured


def _random_state(n, seed=0):
    rng = np.random.default_rng(seed)
    P = rng.normal(0, 3, (n, 2))
    home = rng.normal(0, 3, (n, 2))
    anchor_w = rng.uniform(0.01, 0.1, (n, 1))
    h = rng.uniform(0.3, 2, (n, 2))
    want = h[:, None] + h[None] + L.MARGIN
    W = np.triu(rng.uniform(0, 1, (n, n)) * (rng.uniform(0, 1, (n, n)) < 0.1), 1)
    W = W + W.T
    loops = (np.array([0]), np.zeros((1, 2)), np.array([n - 1]), np.zeros((1, 2)), np.ones(1), np.ones(1))
    return P, home, anchor_w, want, W, loops, rng.normal(0, 0.2, (n, 2))


def _pull(W):
    i, j = np.nonzero(np.triu(W > 0, 1))
    return i, j, W[i, j]


@pytest.mark.parametrize("n", [2, 7, 40, 300])
def test_pair_list_energy_equals_dense_energy(n):
    P, home, aw, want, W, loops, mids = _random_state(n)
    e, g = L._energy(P.ravel(), home, aw, want, _pull(W), loops, 30.0, mids)
    e0, g0 = L._energy_dense(P.ravel(), home, aw, want, W, loops, 30.0, mids)
    assert e == pytest.approx(e0, rel=1e-12)
    np.testing.assert_allclose(g, g0, rtol=1e-9, atol=1e-9)


def test_energy_gradient_matches_finite_differences():
    P, home, aw, want, W, loops, mids = _random_state(12, seed=5)
    ci, cj, _ = _pull(W)
    Q = P + mids
    gap = np.maximum(-(want[ci, cj, 0] - abs(Q[cj, 0] - Q[ci, 0])), -(want[ci, cj, 1] - abs(Q[cj, 1] - Q[ci, 1])))
    assert (gap > 0).any()  # at least one contact pair is pulling
    x = P.ravel()
    _, g = L._energy(x, home, aw, want, _pull(W), loops, 3.0, mids)
    h = 1e-6
    for k in range(len(x)):
        up, down = x.copy(), x.copy()
        up[k] += h
        down[k] -= h
        fd = (
            L._energy(up, home, aw, want, _pull(W), loops, 3.0, mids)[0]
            - L._energy(down, home, aw, want, _pull(W), loops, 3.0, mids)[0]
        ) / (2 * h)
        assert g[k] == pytest.approx(fd, rel=1e-4, abs=1e-4), k


def test_small_layouts_keep_the_dense_energy(monkeypatch):
    from foldmap.cli import make_layout

    called = []
    dense, pair = L._energy_dense, L._energy
    monkeypatch.setattr(L, "_energy_dense", lambda *a, **k: called.append("dense") or dense(*a, **k))
    monkeypatch.setattr(L, "_energy", lambda *a, **k: called.append("pair") or pair(*a, **k))
    layout = make_layout(DATA / "1BL8.cif")[0]  # 12 items, well under the threshold
    assert layout.placed and set(called) == {"dense"}
    called.clear()
    monkeypatch.setattr(L, "_PAIR_ENERGY_ABOVE", 5)
    make_layout(DATA / "1BL8.cif")
    assert set(called) == {"pair"}


def test_push_apart_clears_a_crowd_without_the_fallback_slide(monkeypatch):
    left = []
    real = L._slide_clear

    def spy(items, *a):
        h = [it.half for it in items]
        m = [it.mid for it in items]
        left.append(
            sum(L._collide(items[i], h[i], m[i], items[j], h[j], m[j]) for i in range(len(items)) for j in range(i))
        )
        real(items, *a)

    monkeypatch.setattr(L, "_slide_clear", spy)
    sses, f = _pile(60)
    build_layout(sses, [], f, "projected", attract=False)
    assert left == [0]


def test_fallback_slide_takes_the_nearest_clear_side(monkeypatch):
    monkeypatch.setattr(L, "_SEP_SWEEPS", 0)  # straight to the slide
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    a = fake_sse("H", 1, 9, centroid=(0, 0, 0), axis=(0, 1, 0))
    b = fake_sse("H", 20, 28, centroid=(-2.0, 0, 0), axis=(0, 1, 0))  # just left of a, overlapping
    lay = build_layout([a, b], [], f, "projected", attract=False)
    pa, pb = lay.placed[a.id], lay.placed[b.id]
    assert pb.cx < pa.cx  # pushed out to the left, the short way, not all the way across to the right
    assert pb.rect[2] <= pa.rect[0] + 1e-9 or pb.rect[3] <= pa.rect[1] + 1e-9 or pb.rect[1] >= pa.rect[3] - 1e-9


def test_one_element_layout_is_left_where_it_projects():
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    s = fake_sse("H", 1, 12, centroid=(4.0, -2.0, 0), axis=(0, 1, 0))
    lay = build_layout([s], [], f, "projected", domains=[("only", [s.id])])
    p = lay.placed[s.id]
    assert (p.cx, p.cy) == pytest.approx((4.0 * L.SCALE, -2.0 * L.SCALE))


@pytest.mark.slow
def test_a_32_chain_complex_lays_out_in_reasonable_time(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        from bench import big_complex
    finally:
        sys.path.pop(0)
    from foldmap.cli import make_layout

    path = big_complex(tmp_path / "big32.cif")
    t = time.perf_counter()
    lay, _, _ = make_layout(path)
    assert time.perf_counter() - t < 60  # ~7-11 s here; 96-163 s before the push-apart and energy rewrites
    rects = [p.rect for p in lay.placed.values()]
    for i, a in enumerate(rects):
        for b in rects[i + 1 :]:
            assert min(a[2], b[2]) - max(a[0], b[0]) <= 1e-6 or min(a[3], b[3]) - max(a[1], b[1]) <= 1e-6


def test_bench_leaves_no_temp_dir_behind(monkeypatch, tmp_path):
    import tempfile

    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import bench
    finally:
        sys.path.pop(0)
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    assert bench.main(["2HHB"]) == 0
    assert list(tmp_path.iterdir()) == []


def test_big_complex_has_no_duplicated_chain(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        from bench import big_complex
    finally:
        sys.path.pop(0)
    import gemmi

    names = [c.name for c in gemmi.read_structure(str(big_complex(tmp_path / "big.cif", copies=2)))[0]]
    assert len(names) == 8 and len(set(names)) == 8 and not {"B", "D"} & set(names)


# random-100 C33: the campaign's slowest views (7JU4, 7Y7C, 8VWQ at its 180 s limit; 7JU4's domains view 650 s and
# 20 GB) spent their time in work done twice or pair by pair. Each fast path is checked against the plain one.


def test_saving_the_layout_file_lays_the_figure_out_once(tmp_path, monkeypatch):
    """`plot --save-layout` laid the figure out a second time for the layout file (27 of 7JU4's 183 s): with no
    edits the file is written from the figure's own layout, and says the same."""
    from foldmap import cli, layoutfile
    from foldmap.spec import figure_spec, layout_document

    calls = []
    real = cli.make_layout
    monkeypatch.setattr(cli, "make_layout", lambda *a, **k: calls.append(1) or real(*a, **k))
    out = tmp_path / "l.yaml"
    assert cli.main(["plot", str(DATA / "2HHB.cif"), "-o", str(tmp_path / "f.svg"), "--save-layout", str(out)]) == 0
    assert len(calls) == 1
    monkeypatch.setattr(cli, "make_layout", real)
    fresh = layoutfile.save(tmp_path / "m.yaml", layout_document(figure_spec(DATA / "2HHB.cif")))
    assert out.read_text(encoding="utf-8").splitlines()[1:] == fresh.read_text(encoding="utf-8").splitlines()[1:]


def _dense_overlapping(xy, halves, mids, apart=None):
    q = xy + mids
    depth = halves[:, None, :] + halves[None, :, :] + L.MARGIN - np.abs(q[None, :, :] - q[:, None, :])
    hit = (depth[..., 0] > 1e-9) & (depth[..., 1] > 1e-9)
    if apart is not None:
        hit &= ~apart
    return np.argwhere(np.triu(hit, 1))


@pytest.mark.parametrize("seed", range(6))
def test_overlapping_pairs_of_many_footprints_equal_the_dense_test(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(L._DENSE_PAIRS + 1, 700))
    xy = L._snap(rng.uniform(-50, 50, (n, 2)))
    if seed % 2:
        xy[: n // 3, 0] = xy[0, 0]  # many footprints in one column
    halves, mids = L._snap(rng.uniform(0.2, 3.0, (n, 2))), L._snap(rng.uniform(-1, 1, (n, 2)))
    apart = rng.random((n, n)) < 0.02
    apart |= apart.T
    for a in (None, apart):
        want, got = _dense_overlapping(xy, halves, mids, a), L._overlapping(xy, halves, mids, a)
        assert len(want) > 0 and want.tolist() == got.tolist()


def _plain_panel_fields(owned, free, pad, names=None, turned=()):
    """panels.panel_fields as it was: every domain's raster held, its rivals taken one by one."""
    from scipy.ndimage import distance_transform_edt

    from foldmap import panels as P

    names = names or {}
    every = [b for bs in owned.values() for b in bs] + list(free) + list(names.values())
    every += [(x - hl - hw, y - hl - hw, x + hl + hw, y + hl + hw) for x, y, hl, hw, _ in turned]
    reach = pad + P.CLOSE + P.GUTTER + 1.0
    r = np.asarray(every, float)
    xs = P._axis(r[:, 0].min() - reach, r[:, 2].max() + reach)
    ys = P._axis(r[:, 1].min() - reach, r[:, 3].max() + reach)
    elements = {d: P._distance(bs, xs, ys, reach) for d, bs in owned.items()}
    none = np.full((len(ys), len(xs)), np.inf)
    named = {d: P._distance([names[d]], xs, ys, reach) if d in names else none for d in owned}
    own = {d: np.minimum(elements[d], named[d]) for d in owned}
    loose = np.minimum(P._distance(list(free), xs, ys, reach), P._turned_distance(list(turned), xs, ys, reach))
    values = {}
    for d in owned:
        rival = loose.copy()
        for e in owned:
            if e != d:
                np.minimum(rival, own[e], out=rival)
        base = np.maximum(pad - elements[d], P.NAME_PAD - named[d])
        closed = distance_transform_edt(elements[d] <= pad + P.CLOSE) * P.GRID - P.CLOSE
        shape = np.maximum(base, closed)
        with np.errstate(invalid="ignore"):
            split = np.maximum((rival - own[d] - P.GUTTER) / 2, np.minimum(P.EDGE - own[d], (rival - own[d]) / 2))
        split = np.where(np.isinf(rival), np.inf, split)
        values[d] = np.clip(np.nan_to_num(np.minimum(shape, split), nan=-1.0), -1.0, 1.0)
    return xs, ys, values


@pytest.mark.parametrize("seed", range(8))
def test_panel_fields_equal_the_domain_by_domain_rivals(seed):
    from foldmap.panels import panel_fields

    rng = np.random.default_rng(seed)

    def boxes(n):
        c, h = rng.uniform(0, 30, (n, 2)), rng.uniform(0.2, 2.0, (n, 2))
        return [tuple(map(float, (*(c[k] - h[k]), *(c[k] + h[k])))) for k in range(n)]

    owned = {f"D{k}": boxes(int(rng.integers(1, 6))) for k in range(7)}
    owned["D1"] = list(owned["D0"])  # two domains at the same place: ties in the nearest
    owned["far"] = [(80.0, 80.0, 81.0, 81.5)]  # a domain far from the rest: its own window of the raster
    owned["name only"], owned["nothing"] = [], []
    names = {d: boxes(1)[0] for d in list(owned)[::2]}
    names["name only"] = (60.0, 2.0, 61.5, 2.6)
    free, turned = boxes(3), [(5.0, 5.0, 2.0, 0.5, 0.3)]
    xs, ys, want = _plain_panel_fields(owned, free, 0.4, names, turned)
    got = panel_fields(owned, free, 0.4, names, turned)
    assert got.xs.tolist() == xs.tolist() and got.ys.tolist() == ys.tolist() and list(got.values) == list(want)
    for d in want:
        assert np.array_equal(got.values[d], want[d]), d


def test_panel_fields_keep_each_domain_s_window_not_the_page():
    """Final review A-I1: each domain's field was kept as a page-sized raster, every one at once (7JU4's domains
    view: 4.6 GB). Each is kept as its window and given whole on request, one at a time: 40 small domains spread
    over a page hold less than one page's raster between them, and each reads as it did (-1 outside its window)."""
    from foldmap.panels import panel_fields

    at = {f"D{k}": (10.0 * (k % 8), 10.0 * (k // 8)) for k in range(40)}
    owned = {d: [(x, y, x + 1.5, y + 1.0)] for d, (x, y) in at.items()}
    got = panel_fields(owned, [], 0.4)
    page = got.xs.size * got.ys.size * 8
    assert got.values.nbytes < page
    field = got.values["D9"]
    assert field.shape == (got.ys.size, got.xs.size) and field.max() > 0 and (field == -1.0).mean() > 0.9


def _plain_text_collisions(drawn, page):
    """figcheck.text_collisions as it was: every pair of texts, every text against every body and line."""
    from itertools import combinations

    from matplotlib.transforms import Bbox

    from foldmap import figcheck as F

    lay = page.layout
    over = {frozenset(f"el:{k}" for k in pair) for pair in F.lying_over(lay)}
    own_lines, hits = page.own_lines(), []
    lines = [
        (gid, polys, Bbox.union([Bbox.from_extents(*p.min(axis=0), *p.max(axis=0)) for p in polys]), pad, covers)
        for gid, polys, z, pad in page.lines
        if polys
        for covers in [[p for _, p, bz, opaque in page.bodies if opaque and bz > z]]
    ]
    for (ta, ga, a), (tb, gb, b) in combinations(page.texts, 2):
        if ga != gb and F._texts_meet(a, b) and F._texts_meet(*page.pair_boxes(ta, a, tb, b)):
            hits.append(f"{ga} × {gb}")
    for t, name, box in page.texts:
        gid, met = t.get_gid() or "", set()
        mine = F._owner(gid)
        for body, path, _, _ in page.bodies:
            if body.startswith("na:stem:"):
                continue
            theirs = F._owner(body) or body
            if theirs == mine or theirs in met or frozenset((theirs, mine)) in over:
                continue
            if F._on(box, path):
                met.add(theirs)
                hits.append(f"{name} × {body.replace('helix-back:', 'helix:', 1)}")
        for body, path, own in page.stems:
            if gid not in own and F._on(box, path):
                hits.append(f"{name} × {body}")
        for line, polys, extent, pad, covers in lines:
            reach = Bbox.from_extents(box.x0 - pad, box.y0 - pad, box.x1 + pad, box.y1 + pad)
            if line in own_lines.get(gid, ()) or not extent.overlaps(reach):
                continue
            if F._crosses(polys, reach, covers):
                hits.append(f"{name} × {line}")
    return hits


def _plain_crossings(lay, loops):
    from foldmap.route import _crosses_body

    def cross(p, q, r, s):
        def side(a, b, c):
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

        return side(r, s, p) * side(r, s, q) < 0 and side(p, q, r) * side(p, q, s) < 0

    segs = [list(zip(lp.points, lp.points[1:])) for lp in loops]
    n = sum(
        cross(p, q, r, s)
        for i in range(len(segs))
        for j in range(i + 1, len(segs))
        for p, q in segs[i]
        for r, s in segs[j]
    )
    honest = n
    for lp, ss in zip(loops, segs):
        a, b = lay.placed.get(lp.a_id), lay.placed.get(lp.b_id)
        ports = [p for q, p in ((a, lp.points[0]), (b, lp.points[-1])) if q is not None and q.over]
        for k, e in lay.placed.items():
            if k in (lp.a_id, lp.b_id) or not any(_crosses_body(p, q, e) for p, q in ss):
                continue
            x0, y0, x1, y1 = e.rect
            if e.sse.kind == "E" and any(x0 <= x <= x1 and y0 <= y <= y1 for x, y in ports):
                continue
            honest += 1
    return n, honest


@pytest.mark.parametrize(
    ("path", "style"),
    [("8UUP.cif", ["font_scale=4"]), ("6ZS5.cif", []), ("na/1KX5.cif.gz", []), ("2HHB.cif", ["font_scale=5"])],
)
def test_text_checks_and_crossings_equal_their_pair_by_pair_counts(path, style):
    """figcheck's text collisions and loop crossings look only at pairs whose boxes meet (7JU4: 90 s of its 183 s
    were the checks); the findings are the same, in the same order."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(DATA / path, style=style))
    page = figcheck.measure(d)
    want = _plain_text_collisions(d, page)  # text over text, bodies and lines: 15 on 8UUP, 13 on 2HHB at these sizes
    assert figcheck.text_collisions(d, page=page) == want
    assert (figcheck.loop_crossings(d.loops), figcheck.honest_crossings(d.layout, d.loops)) == _plain_crossings(
        d.layout, d.loops
    )


def test_nucleic_acid_pair_counts_equal_their_pair_by_pair_counts():
    """nacheck.render_checks' letter and bead overlaps look only at pairs that can meet (7Y7C: 30 s); the same."""
    from itertools import combinations

    from matplotlib.backends.backend_agg import FigureCanvasAgg

    from foldmap import nacheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(DATA / "na" / "1KX5.cif.gz"))
    ax = d.figure.axes[0]
    rng = np.random.default_rng(0)  # letters and beads crowded on purpose (the drawn ones never overlap)
    for k, (x, y) in enumerate(rng.uniform(0, 2, (80, 2))):
        ax.text(x, y, "G", fontsize=6).set_gid(f"na:nt:Z:{k}")
    ax.scatter(*rng.uniform(0, 3, (2, 40)), s=30).set_gid("na:contact:test")
    got = nacheck.render_checks(d.layout, d.loops, d.figure, d.sses)
    r = FigureCanvasAgg(d.figure).get_renderer()
    letters = [
        t.get_window_extent(r) for t in ax.texts if (t.get_gid() or "").startswith("na:nt:") and t.get_text().strip()
    ]
    bad = 0
    for a, b in combinations(letters, 2):
        w, h = min(a.x1, b.x1) - max(a.x0, b.x0), min(a.y1, b.y1) - max(a.y0, b.y0)
        bad += w > 0 and h > 0 and w * h > min(a.width * a.height, b.width * b.height) / 3
    centres, radii = [], []
    for c in ax.collections:
        if (c.get_gid() or "").startswith("na:contact:"):
            xy = ax.transData.transform(c.get_offsets())
            centres += list(xy)
            radii += [float(np.sqrt(c.get_sizes()[0])) / 2 * d.figure.dpi / 72] * len(xy)
    beads = sum(
        np.hypot(*(centres[i] - centres[j])) < radii[i] + radii[j] - 1e-6
        for i, j in combinations(range(len(centres)), 2)
    )
    assert len(letters) > 100 and len(centres) > 20 and bad > 5 and beads > 5
    assert (got["letter_overlaps"], got["bead_overlaps"]) == (bad, beads)


def test_parting_domains_measures_each_footprint_once(monkeypatch):
    """_part_domains moves items, never reshapes them: each item's footprint (_Item.box) is measured once, not four
    times an item for every pair of domains in every sweep (7JU4's domains view: 165 s)."""
    from foldmap.spec import draw_spec, figure_spec

    calls, inside, parts = [], [False], []
    real, part = L._Item.box, L._part_domains
    monkeypatch.setattr(L._Item, "box", property(lambda it: (inside[0] and calls.append(1)) or real.fget(it)))

    def parting(items, *a, **k):
        inside[0] = True
        try:
            return part(items, *a, **k)
        finally:
            inside[0] = False
            parts.append(len(items))

    monkeypatch.setattr(L, "_part_domains", parting)
    draw_spec(figure_spec(DATA / "r100" / "3L9G-A.cif.gz", domains="auto"))
    assert parts and len(calls) == sum(parts)
