"""The layout optimiser: a numpy-only L-BFGS on inputs snapped to a grid, so a layout is the same on any CPU,
BLAS kernel and numpy/scipy version."""

import inspect
from pathlib import Path

import numpy as np
import pytest

import foldmap.layout as L

DATA = Path(__file__).parent / "data"


def _quadratic(n=12, seed=3):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(n, n))
    A = q.T @ q + n * np.eye(n)  # symmetric positive definite, condition number ~ 10
    b = rng.normal(size=n)

    def f(x):
        Ax = A @ x
        return float(0.5 * x @ Ax - b @ x), Ax - b

    return f, np.linalg.solve(A, b), rng.normal(size=n)


def _rosen(x):
    """Rosenbrock function and gradient, elementwise only."""
    a, b = x[:-1], x[1:]
    f = float((100.0 * (b - a * a) ** 2 + (1.0 - a) ** 2).sum())
    g = np.zeros_like(x)
    g[:-1] = -400.0 * a * (b - a * a) - 2.0 * (1.0 - a)
    g[1:] += 200.0 * (b - a * a)
    return f, g


def _scipy(f, x0):
    from scipy.optimize import minimize

    return minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 300}).x


def test_lbfgs_minimises_a_quadratic_like_scipy():
    f, best, x0 = _quadratic()
    x = L._lbfgs(f, x0)
    np.testing.assert_allclose(x, best, atol=1e-6)
    assert f(x)[0] <= f(_scipy(f, x0))[0] + 1e-10


@pytest.mark.parametrize("x0", [[-1.2, 1.0], [-1.2, 1.0, -0.5, 0.8, 1.5, -0.3, 0.2, 2.0]])
def test_lbfgs_minimises_rosenbrock_like_scipy(x0):
    x0 = np.array(x0)
    x = L._lbfgs(_rosen, x0)
    ref = _scipy(_rosen, x0)
    assert _rosen(x)[0] < 1e-9 and _rosen(ref)[0] < 1e-9  # both reach the minimum at (1, ..., 1)
    np.testing.assert_allclose(x, np.ones_like(x0), atol=1e-4)


def test_lbfgs_stops_at_maxiter_and_leaves_its_start_alone():
    x0 = np.array([-1.2, 1.0])
    calls = []

    def counted(x):
        calls.append(x.copy())
        return _rosen(x)

    x = L._lbfgs(counted, x0, maxiter=3)
    assert x0.tolist() == [-1.2, 1.0]
    assert _rosen(x)[0] < _rosen(x0)[0]
    assert len(calls) <= 1 + 3 * 20  # three iterations, at most 20 evaluations per line search
    assert L._lbfgs(_rosen, x0, maxiter=3).tolist() == x.tolist()  # deterministic


def test_lbfgs_stops_at_once_on_a_flat_start():
    calls = []

    def flat(x):
        calls.append(1)
        return float((x * x).sum()), 2 * x

    assert L._lbfgs(flat, np.zeros(4)).tolist() == [0.0] * 4
    assert len(calls) == 1


_NO_BLAS = ("@", "np.dot", "np.matmul", "np.linalg", "np.inner", "np.vdot", "np.einsum", ".dot(", "np.hypot")


@pytest.mark.parametrize("fn", ["_lbfgs", "_wolfe", "_energy", "_energy_dense", "_pair_grad"])
def test_the_optimiser_makes_no_blas_calls(fn):
    """BLAS results depend on the CPU kernel: dot products are (a * b).sum() (numpy's own pairwise sum), and
    distances np.sqrt(x * x + y * y) (correctly rounded, unlike a platform's hypot)."""
    src = inspect.getsource(getattr(L, fn))
    assert not [s for s in _NO_BLAS if s in src]


def test_lbfgs_runs_with_numpy_blas_entry_points_removed(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("BLAS called")

    for mod, name in [(np, "dot"), (np, "matmul"), (np, "inner"), (np, "vdot"), (np.linalg, "norm")]:
        monkeypatch.setattr(mod, name, boom)
    x = L._lbfgs(_rosen, np.array([-1.2, 1.0]))
    assert np.abs(x - 1).max() < 1e-4


def test_layout_does_not_use_scipy_minimize(monkeypatch):
    import scipy.optimize

    from foldmap.cli import make_layout

    def boom(*a, **k):
        raise AssertionError("scipy.optimize.minimize called")

    monkeypatch.setattr(scipy.optimize, "minimize", boom)
    assert make_layout(DATA / "1BL8.cif")[0].placed  # dense energy
    monkeypatch.setattr(L, "_PAIR_ENERGY_ABOVE", 5)
    assert make_layout(DATA / "1BL8.cif")[0].placed  # pair-list energy


def _on_grid(v) -> bool:
    v = np.asarray(v, float) * 1e4
    return bool(np.all(np.abs(v - np.round(v)) < 1e-6))


@pytest.mark.parametrize(
    "pdb, focus, restarts",
    [("6ZS5", "auto", True), ("8UUP", "auto", True), ("6Y7S", "none", False), ("1LMB", "auto", False)],
)  # disulfide restarts; copies side by side (spread); a filament's copies rising (lift); DNA
def test_optimiser_inputs_are_snapped_to_a_grid(pdb, focus, restarts, monkeypatch):
    """Positions, sizes, spring ends, spread and lift on a 1e-4 page-unit grid: noise of ~1e-16 in the
    projection (an SVD differs between BLAS builds) cannot reach the optimiser."""
    from foldmap.cli import make_layout

    starts, firsts = [], []  # each run's start, and the energy's inputs at its first evaluation
    real_lbfgs, real_dense = L._lbfgs, L._energy_dense

    def lbfgs(f, x0, **k):
        starts.append(np.array(x0))
        return real_lbfgs(f, x0, **k)

    def dense(flat, home, anchor_w, want, W, loops, lam, mids=None):
        if len(firsts) < len(starts):
            firsts.append((lam, flat, home, want, loops[1], loops[3], mids))
        return real_dense(flat, home, anchor_w, want, W, loops, lam, mids)

    monkeypatch.setattr(L, "_lbfgs", lbfgs)
    monkeypatch.setattr(L, "_energy_dense", dense)
    make_layout(DATA / f"{pdb}.cif", focus=focus)
    assert starts and len(firsts) == len(starts)
    lams = [f[0] for f in firsts]
    # a solve raises the overlap penalty stage by stage, each stage starting where the last one stopped: the
    # first stage of every solve (the initial one and each disulfide restart) starts on the grid
    first_stages = [z for k, z in enumerate(starts) if k == 0 or lams[k] < lams[k - 1]]
    assert len(first_stages) == 1 + (len(starts) - len(L._LAMBDAS)) // (len(L._LAMBDAS) - 1)
    assert (len(first_stages) > 1) == restarts
    assert all(_on_grid(z) for z in first_stages)
    for _, _, home, want, oi, oj, mids in firsts:
        assert _on_grid(home) and _on_grid(want) and _on_grid(oi) and _on_grid(oj) and _on_grid(mids)
    assert _on_grid(firsts[0][1])  # spread and lift too: every copy starts on the grid


def test_snap_rounds_to_the_grid_and_is_idempotent():
    v = np.array([0.123456789, -3.00004999, 12.34565001, 1e-17])
    s = L._snap(v)
    assert s.tolist() == [0.1235, -3.0, 12.3457, 0.0]
    assert L._snap(s).tolist() == s.tolist()
    assert L._snap(v + 1e-15).tolist() == s.tolist()  # noise far below the grid is gone
