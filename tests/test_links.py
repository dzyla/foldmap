"""Disulfide bridges and glycans: found from the file's geometry, drawn at their residues."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_figure, make_layout
from foldmap.io import load_backbone, load_links
from foldmap.style import Style

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
    return [
        a
        for a in [*fig.axes[0].patches, *fig.axes[0].lines, *fig.axes[0].collections]
        if (a.get_gid() or "").startswith(prefix)
    ]


@pytest.mark.parametrize("name", ["5NKT", "6ZS5"])
def test_disulfides_are_drawn_between_their_residues(name):
    from foldmap.render import disulfide_points
    from foldmap.route import route_loops

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
        want = np.array(disulfide_points(lay, loops, sses, lay.res_chain, i, j))
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
    from foldmap.render import residue_point
    from foldmap.route import route_loops

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
    import foldmap.cli as cli
    from foldmap.render import residue_point
    from foldmap.route import route_loops

    path = DATA / f"{name}.cif"

    def mean_bar():
        lay, sses, bb = make_layout(path)
        loops = route_loops(lay, sses, bb)
        ends = [
            (residue_point(lay, loops, sses, bb, i), residue_point(lay, loops, sses, bb, j))
            for i, j in lay.links.disulfides
        ]
        return np.mean([np.hypot(a[0] - b[0], a[1] - b[1]) for a, b in ends])

    pulled = mean_bar()
    monkeypatch.setattr(cli, "_bridge_contacts", lambda *a: {})
    monkeypatch.setattr(cli, "_bridge_springs", lambda *a: [])
    free = mean_bar()
    assert pulled < free if strict else pulled <= free + 1e-3  # shorter where blocks can move, never longer


@pytest.mark.parametrize("name", ["2JTY", "8UTF"])
def test_cysteine_springs_bring_bridged_residues_closer_than_box_contacts(name, monkeypatch):
    import foldmap.cli as cli
    from foldmap.render import residue_point
    from foldmap.route import route_loops

    path = DATA / f"{name}.cif"

    def bars():
        lay, sses, bb = make_layout(path)
        loops = route_loops(lay, sses, bb)
        ends = [
            (residue_point(lay, loops, sses, bb, i), residue_point(lay, loops, sses, bb, j))
            for i, j in lay.links.disulfides
        ]
        return np.array([np.hypot(a[0] - b[0], a[1] - b[1]) for a, b in ends if a and b])

    sprung = bars()
    monkeypatch.setattr(cli, "_bridge_springs", lambda *a: [])
    boxed = bars()
    if name == "2JTY":  # the N-terminal Cys reaches across the sandwich: relocation makes the bar short
        assert sprung.mean() < 0.5 * boxed.mean()
    else:  # many bridges sit inside one sheet or tie symmetry copies: still shorter on average and at worst
        assert sprung.mean() < boxed.mean() and sprung.max() < boxed.max()


@pytest.mark.parametrize("name", ["6ZS5", "8UUP", "6ZYA", "8UTF"])
def test_glycans_stay_off_elements(name):
    from foldmap.cli import make_figure, make_layout

    path = Path(__file__).parent / "data" / f"{name}.cif"
    lay, _, _ = make_layout(path)
    fig = make_figure(path)
    inv = fig.axes[0].transData.inverted()
    rects = [p.rect for p in lay.placed.values()]
    bad = []
    for a in fig.axes[0].patches:
        gid = a.get_gid() or ""
        if gid.startswith("glycan:") and not gid.endswith(":stem"):
            b = a.get_window_extent().transformed(inv)
            for x0, y0, x1, y1 in rects:
                if min(b.x1, x1) - max(b.x0, x0) > 0.05 and min(b.y1, y1) - max(b.y0, y0) > 0.05:
                    bad.append(gid)
    assert not bad, bad


def _bars(name):
    from foldmap.cli import make_figure

    fig = make_figure(Path(__file__).parent / "data" / f"{name}.cif")
    out = []
    for a in fig.axes[0].patches:
        if (a.get_gid() or "").startswith("disulfide:"):
            v = a.get_path().vertices
            out.append(float(np.hypot(*(v[-1] - v[0]))))
    return out


@pytest.mark.parametrize("name", ["8UUP", "6Y7S", "6ZS5"])
def test_loop_cysteines_sit_where_their_loop_passes_closest_to_the_partner(name):
    bars = _bars(name)
    assert bars and min(bars) >= 0.49  # always a visible bar
    # 6ZS5's Cys335-Cys425 joins strand A's loop to helix H across the C-F-G sheet: its strand order is fixed by
    # backbone H-bonds, so that one bar has to span three strands; every other bar is short
    # and Cys527-Cys582 joins the G-H loop (between two sheets) to the short C-terminal stub of strand J
    allowed = 2 if name == "6ZS5" else 0
    assert sum(b > 3.0 for b in bars) <= allowed, sorted(bars)


GLY = Path(__file__).parent / "data" / "glycans"


def _links(name):
    path = GLY / f"{name}.cif"
    bb = load_backbone(path, "asu")
    return bb, load_links(path, bb, "asu")


def test_sialic_acid_joins_its_glycan_through_c2():
    """IgG1 Fc (4BYH): the N297 glycan ends Gal-a2,6-Neu5Ac; Neu5Ac's glycosidic bond leaves from C2."""
    bb, links = _links("4BYH_AC")
    ((site, names),) = links.glycans
    assert (bb.labels[site].name, bb.labels[site].seq) == ("ASN", 297)
    assert sorted(names) == sorted(["NAG", "NAG", "FUC", "BMA", "MAN", "MAN", "NAG", "NAG", "GAL", "SIA"])
    assert all(lig.name != "SIA" for lig in links.ligands)  # drawn once, in the glycan


def test_kdn_is_a_sialic_acid_of_its_glycan(tmp_path):
    """Kdn (as 4BYH's Neu5Ac renamed) joins through C2 like Neu5Ac and is drawn as SNFG's green diamond."""
    import gemmi

    from foldmap.io import ANOMERIC, SUGARS
    from foldmap.render import _SNFG

    st = gemmi.read_structure(str(GLY / "4BYH_AC.cif"))
    for res in (r for ch in st[0] for r in ch if r.name == "SIA"):
        res.name = "KDN"
    path = tmp_path / "kdn.cif"
    st.make_mmcif_document().write_file(str(path))
    bb = load_backbone(path, "asu")
    links = load_links(path, bb, "asu")
    ((_, names),) = links.glycans
    assert "KDN" in names and all("KDN" not in lig.codes for lig in links.ligands)
    assert _SNFG["KDN"] == ("diamond", "#00a651")
    assert set(ANOMERIC) <= SUGARS  # only sugars that can be in a glycan


def test_free_sugars_are_ligands():
    _, mbp = _links("1ANF")
    assert [(lig.symbol, lig.codes) for lig in mbp.ligands] == [("Glc×2", ["GLC", "GLC"])]
    assert len(mbp.ligands[0].contacts) >= 10  # the binding cleft holds it
    _, ha = _links("2WR7_site")
    (analogue,) = ha.ligands
    assert analogue.symbol == "Gal·GlcNAc·Gal·Neu5Ac" and not ha.glycans


def _seg(x0, x1, y=0.0):
    return (np.array([x0, y]), np.array([x1, y]))


def test_bars_on_top_of_each_other_take_separate_lanes():
    from foldmap.render import disulfide_lanes

    assert disulfide_lanes([_seg(0, 2), _seg(1, 3)], 0.1) == [0, 1]
    assert disulfide_lanes([_seg(0, 4), _seg(1, 3), _seg(2, 5)], 0.1) == [0, 1, 2]
    assert disulfide_lanes([_seg(3, 1), _seg(2, 0)], 0.1) == [0, 1]  # the direction of a bar does not matter


def test_bars_that_do_not_cover_each_other_stay_straight():
    from foldmap.render import disulfide_lanes

    assert disulfide_lanes([_seg(0, 1), _seg(2, 3)], 0.1) == [0, 0]  # disjoint spans, same height
    assert disulfide_lanes([_seg(0, 2), _seg(0, 2, y=0.5)], 0.1) == [0, 0]  # parallel but clearly apart
    cross = [(np.array([0.0, -1.0]), np.array([0.0, 1.0])), _seg(-1, 1)]
    assert disulfide_lanes(cross, 0.1) == [0, 0]  # crossing at 90 degrees
    slanted = (np.array([0.0, 0.0]), np.array([2.0, 0.6]))
    assert disulfide_lanes([_seg(0, 2), slanted], 0.1) == [0, 0]  # about 17 degrees apart


def test_a_chain_of_overlaps_reuses_lanes():
    from foldmap.render import disulfide_lanes

    assert disulfide_lanes([_seg(0, 2), _seg(1, 3), _seg(2.5, 4)], 0.1) == [0, 1, 0]


def test_the_bracket_leaves_on_the_side_with_less_in_the_way():
    from foldmap.render import disulfide_side

    a, b = np.array([0.0, 0.0]), np.array([3.0, 0.0])
    assert disulfide_side(a, b, 1, [], 0.1) == -1  # nothing in the way: below (the normal of a -> b points up)
    assert disulfide_side(b, a, 1, [], 0.1) == 1  # the same side whichever way the bar is read
    assert disulfide_side(a, b, 1, [(0.5, -1.0, 2.5, -0.1)], 0.1) == 1  # a box below pushes it above
    assert disulfide_side(np.array([0.0, 0.0]), np.array([0.0, 3.0]), 1, [], 0.1) == -1  # vertical: to the right


def test_overlapping_bars_are_drawn_as_brackets_that_keep_their_cysteines(monkeypatch):
    import foldmap.render as render

    path = DATA / "6ZS5.cif"
    pairs = load_links(path, load_backbone(path)).disulfides
    spots = {}
    for k, pair in enumerate(pairs):  # the first two lie on top of each other, the others far apart
        y = 0.0 if k < 2 else 50.0 * k
        spots[pair] = (np.array([1.0, y]), np.array([4.0 + k, y]))
    monkeypatch.setattr(render, "disulfide_points", lambda lay, loops, sses, chains, i, j: spots[(i, j)])
    fig = make_figure(path)
    bars = {b.get_gid(): b.get_path().vertices for b in _gids(fig, "disulfide:")}
    first, second = (bars[f"disulfide:{i}-{j}"] for i, j in pairs[:2])
    assert len(first) == 2 and len(second) == 4
    a, b = spots[pairs[1]]
    assert np.allclose(second[0], a) and np.allclose(second[-1], b)  # still ends on its cysteines
    assert abs(second[1][1] - a[1]) >= 0.1 and np.isclose(second[1][1], second[2][1])  # a lane clear of the first bar
    assert all(len(v) == 2 for g, v in bars.items() if g not in {f"disulfide:{i}-{j}" for i, j in pairs[:2]})
    dots = {d.get_gid(): d.get_offsets() for d in _gids(fig, "ss-dot:")}
    i, j = pairs[1]
    assert np.allclose(dots[f"ss-dot:{i}-{j}"], [a, b])


def test_4ch9_bond_lying_on_another_detours_in_a_bracket():
    fig = make_figure(DATA / "4CH9.cif")
    bars = {b.get_gid(): b.get_path().vertices for b in _gids(fig, "disulfide:")}
    assert [len(v) for _, v in sorted(bars.items())].count(4) == 1
    assert len(bars["disulfide:249-512"]) == 4  # nested under the longer bar, ends still on its cysteines
    assert any(d.get_gid() == "ss-dot:249-512" for d in _gids(fig, "ss-dot:"))
