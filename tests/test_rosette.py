"""Plan 21 Task 1: a cyclic assembly asked for as a rosette (Style.cyclic=rosette) is seen down its axis, copy k drawn
as protomer 0's drawing turned by k * 360 / n degrees about the axis and pushed out just enough to clear."""

import math
from pathlib import Path

import numpy as np
import pytest

from foldmap.frame import RosetteFrame, turn2
from foldmap.layout import _ROSETTE_GAP, _rosette_samples
from foldmap.spec import draw_spec, figure_spec
from foldmap.style import Style

DATA = Path(__file__).parent / "data"


def drawn(path, *settings, **kw):
    return draw_spec(figure_spec(path, style=list(settings), **kw))


@pytest.fixture(scope="module")
def kcsa():
    return drawn(DATA / "1BL8.cif", "cyclic=rosette")


def test_unrolled_is_the_default():
    assert Style().cyclic == "unrolled"
    d = drawn(DATA / "1BL8.cif")
    assert not d.layout.rosette and type(d.layout.frame).__name__ == "CylinderFrame"


def test_rosette_frame_turns_each_protomer_onto_protomer_0():
    axis, ref = np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0])
    f = RosetteFrame(np.zeros(3), axis, ref, 4, {"A": 0, "B": 1, "C": 2, "D": 3}, {c: 0 for c in "ABCD"})
    assert f.angles == pytest.approx([0, math.pi / 2, math.pi, 3 * math.pi / 2])
    m = np.array([f.u, f.v, f.w])
    assert np.linalg.det(m) == pytest.approx(1.0)  # a rotation, never a mirror
    p = np.array([0.3, 5.0, 2.0])  # a point of protomer 0, toward page-up
    q = np.array([-5.0, 0.3, 2.0])  # the same point a quarter turn on (anticlockwise seen from the axis' tip)
    assert f.project(p, ["A"])[0] == pytest.approx([0.3, 5.0])
    assert f.project(q, ["B"])[0] == pytest.approx(turn2(np.array([0.3, 5.0]), math.pi / 2))
    local = f.local(1)  # protomer 1 seen as protomer 0 is: its point lands where protomer 0's does
    assert local.project(q)[0] == pytest.approx([0.3, 5.0])


def test_each_copy_is_protomer_0_turned(kcsa):
    """KcsA (1BL8, C4): every copy's elements sit at protomer 0's, turned by k * 90 degrees, its angle with them."""
    lay = kcsa.layout
    assert lay.rosette == 4 and type(lay.frame).__name__ == "RosetteFrame"
    assert not lay.rosette_strays
    sector, classes = lay.frame.sector, lay.frame.classes
    zero = {(classes[p.sse.chain], p.sse.first.seq): p for p in lay.placed.values() if sector[p.sse.chain] == 0}
    checked = 0
    for p in lay.placed.values():
        k = sector[p.sse.chain]
        if not k:
            continue
        m = zero[(classes[p.sse.chain], p.sse.first.seq)]
        a = lay.frame.angles[k]
        assert [p.cx, p.cy] == pytest.approx(turn2(np.array([m.cx, m.cy]), a), abs=2e-4)
        assert math.cos(p.angle - m.angle - a) == pytest.approx(1.0, abs=1e-6)
        assert (p.length, p.label) == (m.length, m.label)
        checked += 1
    assert checked == 3 * len(zero)


def test_copies_are_pushed_just_clear(kcsa):
    """The push keeps neighbouring copies _ROSETTE_GAP apart; one step less would not."""
    lay = kcsa.layout
    own = [p for p in lay.placed.values() if lay.frame.sector[p.sse.chain] == 0]
    pts, rad = _rosette_samples(own)
    reach = rad[:, None] + rad[None, :] + _ROSETTE_GAP
    from scipy.spatial.distance import cdist

    for a in lay.frame.angles[1:]:
        assert (cdist(pts, turn2(pts, a)) - reach).min() >= -1e-6
    assert lay.rosette_push > 0


def test_rosette_is_deterministic(kcsa):
    again = drawn(DATA / "1BL8.cif", "cyclic=rosette").layout
    assert [(k, p.cx, p.cy, p.angle) for k, p in kcsa.layout.placed.items()] == [
        (k, p.cx, p.cy, p.angle) for k, p in again.placed.items()
    ]


def test_axis_symbol_and_note(kcsa):
    from foldmap.render import save_svg

    svg = save_svg(kcsa.figure)
    assert 'id="axis:4"' in svg
    items = kcsa.figure.foldmap_provenance
    assert any(i["kind"] == "cyclic" and "4-fold axis" in i["text"] for i in items)


def test_a_two_fold_assembly_ignores_the_key():
    """cyclic=rosette is for n = 3..8: a dimer is drawn as before."""
    a = drawn(DATA / "2HHB.cif", "cyclic=rosette").layout
    b = drawn(DATA / "2HHB.cif").layout
    assert not a.rosette and [(p.cx, p.cy) for p in a.placed.values()] == [(p.cx, p.cy) for p in b.placed.values()]
