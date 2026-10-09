import numpy as np
import pytest

from topoplot.dssp import assign_dssp
from topoplot.frame import view_frame
from topoplot.io import load_backbone
from topoplot.ss import build_sses
from helpers import pipeline


def sses_of(path):
    bb = load_backbone(path)
    return bb, build_sses(bb, assign_dssp(bb).ss)


@pytest.mark.parametrize("name", ["ubq", "zs5", "zya"])
def test_frame_is_a_proper_rotation(name, request):
    _, sses = sses_of(request.getfixturevalue(name))
    f = view_frame(sses)
    m = np.array([f.u, f.v, f.w])
    assert np.allclose(m @ m.T, np.eye(3), atol=1e-9)
    assert np.linalg.det(m) == pytest.approx(1.0)  # never the mirror image


def test_deterministic(zs5):
    _, sses = sses_of(zs5)
    a, b = view_frame(sses), view_frame(sses)
    assert np.array_equal(a.u, b.u) and np.array_equal(a.v, b.v)


def test_project_matches_dot_products(ubq):
    _, sses = sses_of(ubq)
    f = view_frame(sses)
    c = np.array([s.centroid for s in sses])
    p = f.project(c)
    assert np.allclose(p[:, 0], (c - f.origin) @ f.u)
    assert np.allclose(p[:, 1], (c - f.origin) @ f.v)
    assert np.allclose(f.depth(c), (c - f.origin) @ f.w)


def test_rotate_90_turns_content_counter_clockwise(ubq):
    _, sses = sses_of(ubq)
    f0, f90 = view_frame(sses), view_frame(sses, rotate=90)
    pt = (f0.origin + f0.u)[None, :]
    assert np.allclose(f0.project(pt), [[1, 0]], atol=1e-9)
    assert np.allclose(f90.project(pt), [[0, 1]], atol=1e-9)
    assert np.allclose(f90.w, f0.w)


def test_flip_v_is_an_explicit_mirror(ubq):
    _, sses = sses_of(ubq)
    f0, ff = view_frame(sses), view_frame(sses, flip_v=True)
    assert np.allclose(ff.v, -f0.v)
    assert np.dot(np.cross(ff.u, ff.v), ff.w) < 0  # left-handed, only on request


def test_first_element_points_up_and_last_is_to_the_right(zs5):
    _, sses = sses_of(zs5)
    f = view_frame(sses)
    assert sses[0].axis @ f.v >= 0
    p = f.project(np.array([s.centroid for s in sses]))
    assert p[-1, 0] >= p[0, 0]


@pytest.mark.parametrize(
    "pdb,kind,floor",
    [("1C3W", "H", 0.9), ("1EMA", "E", 0.7), ("6ZS5", "E", 0.75), ("1TIM", "E", 0.6), ("2HHB", "H", 0.45)],
)
def test_page_up_is_the_dominant_element_axis(pdb, kind, floor, zs5):
    _, sses = sses_of(zs5.parent / f"{pdb}.cif")
    f = view_frame(sses)
    mine = [s for s in sses if s.kind == kind]
    got = np.mean([abs(s.axis @ f.v) for s in mine])
    scatter = sum(len(s) * np.outer(s.axis, s.axis) for s in mine)
    best_up = np.linalg.eigh(scatter)[1][:, -1]  # the best "up" for this kind alone
    best = np.mean([abs(s.axis @ best_up) for s in mine])
    assert got >= floor  # elements lie along the page (measured values: 0.98, 0.76, 0.82, 0.64, 0.52)
    assert got >= 0.95 * best  # and the shared "up" gives up almost nothing to a per-kind optimum


def test_explicit_up_vector_is_honoured(ubq):
    _, sses = sses_of(ubq)
    f = view_frame(sses, up=[0.0, 0.0, 2.0])
    assert np.allclose(f.v, [0, 0, 1])
    assert abs(f.u @ f.v) < 1e-9 and np.linalg.det(np.array([f.u, f.v, f.w])) == pytest.approx(1.0)


@pytest.mark.parametrize("count", [0, 1, 2])
def test_few_elements_still_give_a_valid_frame(ubq, count):
    _, sses = sses_of(ubq)
    f = view_frame(sses[:count])
    m = np.array([f.u, f.v, f.w])
    assert np.allclose(m @ m.T, np.eye(3), atol=1e-9)
    assert np.linalg.det(m) == pytest.approx(1.0)


def test_explicit_up_is_kept_even_against_the_first_element(ubq):
    _, sses, _ = pipeline(ubq)
    for up in ([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]):
        f = view_frame(sses, up=up)
        assert np.allclose(f.v, up)  # a membrane normal names the side that goes on top; never flip it


def test_zero_up_vector_is_rejected(ubq):
    _, sses, _ = pipeline(ubq)
    with pytest.raises(ValueError, match="up"):
        view_frame(sses, up=[0, 0, 0])
