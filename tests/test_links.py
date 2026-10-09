"""Disulfide bridges and glycans: found from the file's geometry, drawn at their residues."""

from pathlib import Path

import numpy as np
import pytest

from topoplot.cli import make_figure, make_layout
from topoplot.io import load_backbone, load_links
from topoplot.style import Style

DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("name, n_ss", [("5NKT", 1), ("2JTY", 1), ("6ZS5", 5), ("1UBQ", 0)])
def test_disulfides_are_found(name, n_ss):
    bb = load_backbone(DATA / f"{name}.cif", "asu")
    links = load_links(DATA / f"{name}.cif", bb, "asu")
    assert len(links.disulfides) == n_ss
    for i, j in links.disulfides:
        assert bb.labels[i].name == "CYS" and bb.labels[j].name == "CYS" and i != j


def test_glycans_hang_on_asparagine_with_their_sugars_in_order():
    bb = load_backbone(DATA / "6ZS5.cif", "asu")
    links = load_links(DATA / "6ZS5.cif", bb, "asu")
    assert links.glycans
    for r, sugars in links.glycans:
        assert bb.labels[r].name in ("ASN", "SER", "THR") and sugars and sugars[0] in ("NAG", "NDG", "A2G", "NGA")


def _gids(fig, prefix):
    return [a for a in [*fig.axes[0].patches, *fig.axes[0].lines, *fig.axes[0].collections]
            if (a.get_gid() or "").startswith(prefix)]


@pytest.mark.parametrize("name", ["5NKT", "6ZS5"])
def test_disulfides_are_drawn_between_their_residues(name):
    from topoplot.render import residue_point
    from topoplot.route import route_loops

    path = DATA / f"{name}.cif"
    lay, sses, bb = make_layout(path)
    links = load_links(path, bb)
    fig = make_figure(path)
    bars = _gids(fig, "disulfide:")
    assert len(bars) == len(links.disulfides) > 0
    loops = route_loops(lay, sses, bb)
    for bar in bars:
        i, j = (int(x) for x in bar.get_gid().split(":")[1].split("-"))
        ends = bar.get_path().vertices[[0, -1]]
        want = np.array([residue_point(lay, loops, sses, bb, i), residue_point(lay, loops, sses, bb, j)])
        assert np.allclose(sorted(map(tuple, ends)), sorted(map(tuple, want)), atol=1e-6)


def test_glycans_are_drawn_with_snfg_symbols():
    fig = make_figure(DATA / "6ZS5.cif")
    glycans = _gids(fig, "glycan:")
    bb = load_backbone(DATA / "6ZS5.cif")
    assert len({g.get_gid().split(":")[1] for g in glycans}) == len(load_links(DATA / "6ZS5.cif", bb).glycans)


def test_links_can_be_switched_off():
    fig = make_figure(DATA / "6ZS5.cif", look=Style(disulfides=False, glycans=False))
    assert not _gids(fig, "disulfide:") and not _gids(fig, "glycan:")


def test_residue_point_lies_on_its_element():
    from topoplot.render import residue_point
    from topoplot.route import route_loops

    lay, sses, bb = make_layout(DATA / "1UBQ.cif")
    loops = route_loops(lay, sses, bb)
    for s in sses:
        p = lay.placed[s.id]
        for r in (s.start, s.end, (s.start + s.end) // 2):
            q = np.array(residue_point(lay, loops, sses, bb, r))
            a, b = np.array(p.n_port), np.array(p.c_port)
            t = (q - a) @ (b - a) / ((b - a) @ (b - a))
            assert -1e-9 <= t <= 1 + 1e-9 and np.allclose(a + t * (b - a), q)


@pytest.mark.parametrize("name, strict", [("2JTY", True), ("6ZS5", False)])  # 6ZS5: bridges inside one sheet
def test_layout_pulls_bridged_cysteines_together(name, strict, monkeypatch):
    import topoplot.cli as cli
    from topoplot.render import residue_point
    from topoplot.route import route_loops

    path = DATA / f"{name}.cif"

    def mean_bar():
        lay, sses, bb = make_layout(path)
        loops = route_loops(lay, sses, bb)
        ends = [(residue_point(lay, loops, sses, bb, i), residue_point(lay, loops, sses, bb, j))
                for i, j in lay.links.disulfides]
        return np.mean([np.hypot(a[0] - b[0], a[1] - b[1]) for a, b in ends])

    pulled = mean_bar()
    monkeypatch.setattr(cli, "_bridge_contacts", lambda *a: {})
    monkeypatch.setattr(cli, "_bridge_springs", lambda *a: [])
    free = mean_bar()
    assert pulled < free if strict else pulled <= free + 1e-3  # shorter where blocks can move, never longer


@pytest.mark.parametrize("name", ["2JTY", "8UTF"])
def test_cysteine_springs_bring_bridged_residues_closer_than_box_contacts(name, monkeypatch):
    import topoplot.cli as cli
    from topoplot.render import residue_point
    from topoplot.route import route_loops

    path = DATA / f"{name}.cif"

    def bars():
        lay, sses, bb = make_layout(path)
        loops = route_loops(lay, sses, bb)
        ends = [(residue_point(lay, loops, sses, bb, i), residue_point(lay, loops, sses, bb, j))
                for i, j in lay.links.disulfides]
        return np.array([np.hypot(a[0] - b[0], a[1] - b[1]) for a, b in ends if a and b])

    sprung = bars()
    monkeypatch.setattr(cli, "_bridge_springs", lambda *a: [])
    boxed = bars()
    if name == "2JTY":  # the N-terminal Cys reaches across the sandwich: relocation makes the bar short
        assert sprung.mean() < 0.5 * boxed.mean()
    else:  # many bridges sit inside one sheet or tie symmetry copies: still shorter on average and at worst
        assert sprung.mean() < boxed.mean() and sprung.max() < boxed.max()
