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


def _blas_threads():
    from threadpoolctl import threadpool_info

    return [i["num_threads"] for i in threadpool_info() if i["user_api"] == "blas"]


def _pile(n):
    """n helices projected onto one point: the push-apart must untangle them by itself."""
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    return [fake_sse("H", 10 * k, 10 * k + 8, centroid=(0, 0, 0), axis=(0.2, 0.1, 1)) for k in range(n)], f


def test_layout_runs_blas_on_one_thread_and_restores_it(monkeypatch):
    from threadpoolctl import threadpool_limits

    if not _blas_threads():
        pytest.skip("no BLAS library detected")
    with threadpool_limits(2, "blas"):  # not already 1, so 'one inside, restored after' says something
        before = _blas_threads()
        seen = []
        real = L._attract

        def spy(*a, **k):
            seen.append(_blas_threads())
            return real(*a, **k)

        monkeypatch.setattr(L, "_attract", spy)
        sses, f = _pile(4)
        build_layout(sses, [], f, "projected")
        assert seen and all(t == 1 for t in seen[0])
        assert _blas_threads() == before


def test_overlapping_layouts_in_threads_keep_the_limit_and_restore_it():
    import threading

    before = _blas_threads()
    inside, first_out = threading.Barrier(2), threading.Event()
    seen = []

    def first():
        with L._one_blas_thread():
            inside.wait()
        first_out.set()

    def second():
        with L._one_blas_thread():
            inside.wait()
            first_out.wait()
            seen.append(_blas_threads())  # the other layout has finished: this one still runs on one thread

    threads = [threading.Thread(target=first), threading.Thread(target=second)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert seen and all(n == 1 for n in seen[0])
    assert _blas_threads() == before


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

    def spy(items):
        h = [it.half for it in items]
        m = [it.mid for it in items]
        left.append(
            sum(L._collide(items[i], h[i], m[i], items[j], h[j], m[j]) for i in range(len(items)) for j in range(i))
        )
        real(items)

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


def test_a_broken_blas_limit_does_not_fail_the_figure(monkeypatch):
    import threadpoolctl

    def broken(*a, **k):
        raise RuntimeError("broken BLAS")

    monkeypatch.setattr(threadpoolctl, "threadpool_limits", broken)
    with L._one_blas_thread():
        pass
    with L._one_blas_thread():  # the count stayed consistent: the next layout tries again, and also runs
        assert L._BLAS_USERS == 1
    assert L._BLAS_USERS == 0 and L._BLAS_LIMIT is None
    sses, f = _pile(3)
    build_layout(sses, [], f, "projected")


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
