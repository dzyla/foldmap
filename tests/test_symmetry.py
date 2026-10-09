import numpy as np
import pytest
from helpers import pipeline

from foldmap.symmetry import detect_symmetry

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
    return [[by[(k, i, a, b)] for k in range(s.n)] for (k0, i, a, b) in list(by) if k0 == 0
            if all((k, i, a, b) in by for k in range(s.n))]


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


@pytest.mark.parametrize("name", ["8UUP", "2HHB"])
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


def test_symmetry_can_be_switched_off(tmp_path):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(DATA / "2HHB.cif"), "-o", str(out), "--symmetry", "off"]) == 0
