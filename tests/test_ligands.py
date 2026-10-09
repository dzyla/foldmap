"""Ligands and metal ions: found with the residues that hold them, drawn as markers tethered to those residues."""

from pathlib import Path

import pytest

from foldmap.cli import make_figure, make_layout
from foldmap.io import load_backbone, load_links
from foldmap.style import Style

DATA = Path(__file__).parent / "data"


def _ligands(name):
    bb = load_backbone(DATA / f"{name}.cif", "asu")
    return bb, load_links(DATA / f"{name}.cif", bb, "asu").ligands


def test_zinc_fingers_hold_zinc_by_two_cys_two_his():
    bb, ligs = _ligands("1ZAA")
    zn = [g for g in ligs if g.name == "ZN"]
    assert len(zn) == 3 and all(g.metal for g in zn)
    for g in zn:
        assert sorted(bb.labels[r].name for r in g.contacts) == ["CYS", "CYS", "HIS", "HIS"]


def test_heme_sits_on_its_histidines():
    bb, ligs = _ligands("2HHB")
    hem = [g for g in ligs if g.name == "HEM"]
    assert len(hem) == 4 and not any(g.metal for g in hem)
    assert all("HIS" in {bb.labels[r].name for r in g.contacts} for g in hem)


def test_potassium_in_the_selectivity_filter():
    _, ligs = _ligands("1BL8")
    assert [g.name for g in ligs if g.metal] and all(g.contacts for g in ligs if g.name == "K")


def test_modified_residues_are_not_ligands_and_additives_are_flagged():
    _, ligs = _ligands("1EMA")
    assert not any(g.name in ("CRO", "MSE") for g in ligs)
    _, ligs = _ligands("5NKT")
    assert ligs and all(g.additive for g in ligs)  # sulfate from the crystallisation buffer


def _markers(fig):
    return [a for a in fig.axes[0].patches if (a.get_gid() or "").startswith("ligand:")]


@pytest.mark.parametrize("name, n", [("1ZAA", 3), ("2HHB", 4), ("5NKT", 0), ("1UBQ", 0)])
def test_markers_by_default_skip_additives(name, n):
    assert len(_markers(make_figure(DATA / f"{name}.cif", assembly="asu"))) == n


def test_marker_choice():
    path = DATA / "2HHB.cif"
    assert not _markers(make_figure(path, look=Style(ligands="none"), assembly="asu"))
    everything = _markers(make_figure(path, look=Style(ligands="all"), assembly="asu"))
    assert len(everything) == 5  # 4 HEM + the one PO4 touching protein (the other floats free)
    only = _markers(make_figure(path, look=Style(ligands="PO4"), assembly="asu"))
    assert len(only) == 1


def test_markers_keep_clear_of_elements_and_each_other():
    lay, _, _ = make_layout(DATA / "2HHB.cif", assembly="asu")
    fig = make_figure(DATA / "2HHB.cif", assembly="asu")
    boxes = [m.get_extents().transformed(fig.axes[0].transData.inverted()) for m in _markers(fig)]
    for b in boxes:
        for p in lay.placed.values():
            x0, y0, x1, y1 = p.rect
            assert not (b.x0 < x1 - 1e-6 and x0 < b.x1 - 1e-6 and b.y0 < y1 - 1e-6 and y0 < b.y1 - 1e-6), p.label
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            assert not boxes[i].overlaps(boxes[j])


def test_every_marker_is_tethered():
    fig = make_figure(DATA / "1ZAA.cif", assembly="asu")
    tethers = [a for a in fig.axes[0].patches if (a.get_gid() or "").startswith("ligand-tether:")]
    for m in _markers(fig):
        key = m.get_gid().split(":", 1)[1]
        assert any(t.get_gid().startswith(f"ligand-tether:{key}:") for t in tethers)
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert "Zn" in texts


def test_ligands_style_value_checked():
    with pytest.raises(ValueError, match="ligands"):
        Style(ligands="hem,,").validate()


def test_ligands_tether_to_at_most_four_elements():
    fig = make_figure(DATA / "2HHB.cif", assembly="asu")
    tethers = [a.get_gid() for a in fig.axes[0].patches if (a.get_gid() or "").startswith("ligand-tether:HEM")]
    per = {}
    for gid in tethers:
        per.setdefault(gid.rsplit(":", 1)[0], 0)
        per[gid.rsplit(":", 1)[0]] += 1
    assert len(per) == 4 and all(1 <= n <= 4 for n in per.values())
