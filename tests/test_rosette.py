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


# ------------------------------------------------------------------ plan 21 Task 2: routes, labels, ligands, checks


def test_each_copy_loop_is_protomer_0_loop_turned(kcsa):
    """Every copy's loop is its protomer-0 mate's route turned with the copy (routed once, among every copy's
    elements): the copies are exact images, loops included."""
    lay = kcsa.layout
    loops = {(lp.a_id, lp.b_id): lp for lp in kcsa.loops}
    checked = 0
    for (a, b), lp in loops.items():
        if a not in lay.rosette_mates:
            continue
        (ma, k), (mb, _) = lay.rosette_mates[a], lay.rosette_mates[b]
        want = turn2(np.array(loops[(ma, mb)].points), lay.frame.angles[k])
        assert np.abs(np.array(lp.points) - want).max() < 1e-3
        assert lp.fallback == loops[(ma, mb)].fallback
        checked += 1
    assert checked == 3 * sum(1 for a, _ in loops if lay.frame.sector[lay.placed[a].sse.chain] == 0)


def test_axis_ions_one_untethered_cluster_at_the_centre():
    """Ions held by every copy of a rosette (a channel's pore) are marked once each, in one cluster round the axis
    symbol at the page centre, with no tether, clear of the symbol and of each other; one held by fewer copies is
    not on the axis."""
    from types import SimpleNamespace

    from foldmap.layout import AXIS_SYMBOL_R, centred_box, overlaps_any
    from foldmap.model import Ligand
    from foldmap.render import _axis_cluster, _on_axis

    frame = SimpleNamespace(sector={"A": 0, "B": 1, "C": 2, "D": 3})
    lay = SimpleNamespace(rosette=4, frame=frame, res_chain=["A", "B", "C", "D"])
    ions = [Ligand("K", "A", 900 + n, True, "K", contacts=[0, 1, 2, 3]) for n in range(6)]
    assert all(_on_axis(lay, g) for g in ions)
    assert not _on_axis(lay, Ligand("K", "A", 990, True, "K", contacts=[0, 1, 2]))
    symbol = centred_box(0.0, 0.0, AXIS_SYMBOL_R, AXIS_SYMBOL_R)
    marks = _axis_cluster(ions, [symbol], [])
    assert len(marks) == 6
    boxes = []
    for m in marks:
        (x, y), (hw, hh) = m["centre"], m["half"]
        box = (x - hw, y - hh, x + hw, y + hh)
        assert m["targets"] == [] and not m["crowded"] and math.hypot(x, y) < 2.5
        assert not overlaps_any(box, [symbol, *boxes])
        boxes.append(box)


def test_copy_names_keep_off_the_turned_bodies(kcsa):
    """A turned copy's helix name is on a side clear of every other element's (tilted) body where one side is."""
    from foldmap.layout import _GLYPH_W, _box_on_body, _element_label_box, glyph

    lay = kcsa.layout
    k = glyph(lay)[0] / _GLYPH_W
    bodies = [*lay.placed.values(), *lay.ghosts]
    for i in lay.rosette_mates:
        p = lay.placed[i]
        if p.sse.kind != "H":
            continue
        hits = lambda p=p: any(_box_on_body(_element_label_box(p, k), q) for q in bodies if q is not p)  # noqa: E731
        if hits():
            p.label_side = -p.label_side
            assert hits(), f"{i}: its name could stand clear on its other side"
            p.label_side = -p.label_side


@pytest.mark.parametrize("name", ["1BL8", "8UUP"])
def test_rosette_passes_the_checks(name):
    """figcheck reads a rosette in each copy's own frame where it reads page axes (sheet columns and spans), and
    element bodies at their angle: no element overlaps, sheets in order, the science checks, copies alike (A15, by
    the best rigid motion, a rotation)."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec as drawn_spec

    d = drawn_spec(figure_spec(DATA / f"{name}.cif", style=["cyclic=rosette"]))
    assert d.layout.rosette
    assert figcheck.element_overlaps(d.layout) == []
    assert figcheck.order_checks(d) == []
    science = figcheck.science_checks(d)
    assert all(not science[s] for s in ("S1", "S2", "S3", "S6"))
    assert all(dev is not None and dev < 0.5 for dev in figcheck.copies_alike(d).values())


def test_cyclic_round_trips(tmp_path):
    """The cyclic choice is saved with the figure and read back (app and command line alike); the share link carries
    the theme and the structure only, as for every style key."""
    import yaml

    from foldmap.app import layout_document, settings_from_document

    doc = layout_document(DATA / "1UBQ.cif", {"theme": "publication", "cyclic": "rosette"})
    assert doc["style"]["cyclic"] == "rosette"
    assert settings_from_document(doc)["cyclic"] == "rosette"
    lf = tmp_path / "fig.yaml"
    lf.write_text(yaml.safe_dump(doc), encoding="utf-8")
    assert figure_spec(DATA / "1UBQ.cif", layout_file=lf).look.cyclic == "rosette"
    assert "cyclic" not in layout_document(DATA / "1UBQ.cif", {"theme": "publication"})["style"]  # the default


def test_styles_lists_cyclic(capsys):
    from foldmap.cli import main

    main(["styles"])
    assert "cyclic" in capsys.readouterr().out
