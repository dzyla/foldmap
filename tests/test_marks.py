"""Referring to elements (label, chain:label, #index, res:range): marking, renaming, swapping, moving."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap.cli import make_figure, make_layout
from foldmap.style import Style

DATA = Path(__file__).parent / "data"
FIMA = DATA / "2JTY.cif"  # self-complemented FimA: the donor strand is the last one, residues 167-182


def _fill(fig, sid):
    for a in fig.axes[0].patches:
        if (a.get_gid() or "") in (f"strand:{sid}", f"helix:{sid}", f"eta:{sid}"):
            return to_hex(a.get_facecolor())
    raise KeyError(sid)


def test_resolve_element_references():
    from foldmap.layout import resolve

    lay, sses, _ = make_layout(FIMA)
    donor = next(s.id for s in sses if s.first.seq == 167)
    label = lay.placed[donor].label
    assert resolve(lay, sses, label) == [donor]
    assert resolve(lay, sses, f"A:{label}") == [donor]
    assert resolve(lay, sses, f"#{[s.id for s in sses].index(donor) + 1}") == [donor]
    assert resolve(lay, sses, "res:170-175") == [donor]
    with pytest.raises(ValueError, match="Z"):
        resolve(lay, sses, "Z")


def test_mark_paints_the_named_elements_only():
    lay, sses, _ = make_layout(FIMA)
    donor = next(s.id for s in sses if s.first.seq == 167)
    plain = make_figure(FIMA)
    fig = make_figure(FIMA, look=Style(mark="res:167-182=#2ca02c"))
    assert _fill(fig, donor) == "#2ca02c"
    others = [s.id for s in sses if s.id != donor]
    assert all(_fill(fig, k) == _fill(plain, k) for k in others)


def test_mark_without_a_colour_uses_the_mark_colour():
    from foldmap.render import MARK_DEFAULT

    lay, sses, _ = make_layout(FIMA)
    donor = next(s.id for s in sses if s.first.seq == 167)
    fig = make_figure(FIMA, look=Style(mark=lay.placed[donor].label))
    assert _fill(fig, donor) == MARK_DEFAULT


def test_bad_mark_reference_is_an_error():
    with pytest.raises(ValueError, match="Q9"):
        make_figure(FIMA, look=Style(mark="Q9"))


def test_highlight_legend_is_compact():
    fig = make_figure(DATA / "assemblies" / "7A4M.cif", look=Style(highlight="asu"))
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert sum(t.startswith("Chain ") for t in texts) == 1 and any("symmetry copies" in t for t in texts)


def test_rename_sets_element_labels(tmp_path):
    from foldmap.cli import main

    lay, sses, _ = make_layout(FIMA, rename={"res:167-182": "Gd"})
    donor = next(s.id for s in sses if s.first.seq == 167)
    assert lay.placed[donor].label == "Gd"
    out = tmp_path / "f.svg"
    assert main(["plot", str(FIMA), "-o", str(out), "--rename", "res:167-182=Gd"]) == 0
    assert ">Gd<" in out.read_text()


def test_swap_exchanges_two_elements():
    lay, sses, _ = make_layout(DATA / "2LZM.cif")
    a, b = "A:3-11", "A:60-79"
    before = {k: (lay.placed[k].cx, lay.placed[k].cy) for k in (a, b)}
    lab = {k: lay.placed[k].label for k in (a, b)}
    swapped, _, _ = make_layout(DATA / "2LZM.cif", swap=[(lab[a], lab[b])])
    for k, other in ((a, b), (b, a)):  # each now sits nearer the other's old place than its own
        to_other = np.hypot(swapped.placed[k].cx - before[other][0], swapped.placed[k].cy - before[other][1])
        to_own = np.hypot(swapped.placed[k].cx - before[k][0], swapped.placed[k].cy - before[k][1])
        assert to_other < to_own


def test_move_shifts_an_element_and_nothing_overlaps():
    lay, sses, _ = make_layout(DATA / "2LZM.cif")
    k = "A:3-11"
    moved, _, _ = make_layout(DATA / "2LZM.cif", move=[(lay.placed[k].label, (0.0, 5.0))])
    assert moved.placed[k].cy - lay.placed[k].cy > 3.0
    items = list(moved.placed.values())
    for i, p in enumerate(items):
        for q in items[i + 1 :]:
            x = min(p.rect[2], q.rect[2]) - max(p.rect[0], q.rect[0])
            y = min(p.rect[3], q.rect[3]) - max(p.rect[1], q.rect[1])
            assert not (x > 1e-6 and y > 1e-6)


def test_cli_swap_move_and_view(tmp_path):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(DATA / "2LZM.cif"), "-o", str(out), "--swap", "α1,α3", "--move", "α2=1,-2"]) == 0
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(out), "--up", "0,0,1", "--view", "1,0,0"]) == 0
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(out), "--swap", "α1"]) == 1
