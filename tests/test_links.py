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


def _on_outline(q, path, tol=0.03) -> bool:
    """Is point q within tol page units of the closed outline `path`?"""
    v = path.vertices
    best = np.inf
    for a, b in zip(v, v[1:]):
        ab = b - a
        t = float(np.clip(np.dot(q - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0))
        best = min(best, float(np.hypot(*(q - a - t * ab))))
    return best <= tol


@pytest.mark.parametrize("name", ["5NKT", "6ZS5"])
def test_disulfides_are_drawn_between_their_residues(name):
    """Each bar joins its two cysteines: one in the open (a loop) keeps its sulfur dot on its own loop; one in an
    element, whose body hides the bar's end, gets its end mark where the bar comes out from under that body, on its
    outline (plan 10 Task 6)."""
    from foldmap.render import link_bodies, loop_of, residue_point
    from foldmap.route import route_loops

    path = DATA / f"{name}.cif"
    lay, sses, bb = make_layout(path)
    links = load_links(path, bb)
    fig = make_figure(path)
    bars = _gids(fig, "disulfide:")
    assert len(bars) == len(links.disulfides) > 0
    loops = route_loops(lay, sses, bb)
    dots = {d.get_gid(): d.get_offsets() for d in _gids(fig, "ss-dot:")}
    ends = {d.get_gid(): d.get_offsets() for d in _gids(fig, "disulfide-end:")}
    outline = dict(zip([*lay.placed, *(f"ghost:{g.sse.id}" for g in lay.ghosts)], link_bodies(lay)))
    for bar in bars:
        i, j = (int(x) for x in bar.get_gid().split(":")[1].split("-"))
        drawn = bar.get_path().vertices
        marks = [*dots.get(f"ss-dot:{i}-{j}", []), *ends.get(f"disulfide-end:{i}-{j}", [])]
        assert len(marks) == 2
        for r in (i, j):
            loop = loop_of(lay, loops, sses, lay.res_chain, r)
            home = next((s.id for s in sses if s.start <= r <= s.end and s.id in lay.placed), None)
            if loop is not None or home is None:  # in a loop or past a chain end: a dot on its loop or at its place
                at = residue_point(lay, loops, sses, lay.res_chain, r)
                on = [m for m in dots[f"ss-dot:{i}-{j}"] if loop is None and np.allclose(m, at)]
                on += [m for m in dots[f"ss-dot:{i}-{j}"] if loop is not None and _on_outline(m, _Poly(loop), 1e-6)]
                assert on, (r, name)
            else:  # in an element whose body hides it: an end mark where the bar comes out (at the body's
                # edge, or just past it where text keeps the bar off), a stretch of the bar starting there; at a
                # tip the body does not cover, its dot
                at = residue_point(lay, loops, sses, lay.res_chain, r)
                dot = [m for m in dots.get(f"ss-dot:{i}-{j}", []) if np.allclose(m, at)]
                if not dot:
                    assert outline[home].contains_point(at), (r, name)
                    mark = [m for m in ends[f"disulfide-end:{i}-{j}"] if np.hypot(*(m - at)) < 1.0]
                    assert mark and any(np.allclose(mark[0], v) for v in drawn), (r, name)
                    assert not outline[home].contains_point(mark[0]) or _on_outline(mark[0], outline[home])


class _Poly:
    def __init__(self, pts):
        self.vertices = np.asarray(pts, float)


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
    """Each disulfide's length on the page, cysteine to cysteine (its route's ends: the drawn bar shows only where
    no body or text covers it)."""
    from foldmap.render import disulfide_bars
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(Path(__file__).parent / "data" / f"{name}.cif")
    loops = route_loops(lay, sses, bb)
    return [float(np.hypot(*np.subtract(b, a))) for _, _, a, b, _ in disulfide_bars(lay, loops, sses, Style(), 1, 1)]


@pytest.mark.parametrize("name", ["8UUP", "6Y7S", "6ZS5"])
def test_loop_cysteines_sit_where_their_loop_passes_closest_to_the_partner(name):
    bars = _bars(name)
    assert bars and min(bars) >= 0.49  # always a visible bar
    # 6ZS5's Cys335-Cys425 joins strand A's loop to helix H across the C-F-G sheet: its strand order is fixed by
    # backbone H-bonds, so that one bar has to span three strands; every other bar is short
    # and Cys527-Cys582 joins the G-H loop (between two sheets) to the short C-terminal stub of strand J
    allowed = 2 if name == "6ZS5" else 0
    if name == "8UUP":  # F2's Cys68 (the loop into α2) to F1's Cys195 (α5's first residue): α2 and α5 are one bundle,
        # staggered 2.7 page units as in 3D, and since each protomer is drawn alike (test_copies) the loop meets α2
        # from the side, so each copy's bar spans the bundle's stagger: three alike bars just over 3
        assert sum(b > 3.0 for b in bars) <= 3 and max(bars) <= 3.1, sorted(bars)
        return
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
    assert disulfide_side(a, b, 1, [], 0.1, 0.3) == -1  # nothing in the way: below (the normal of a -> b points up)
    assert disulfide_side(b, a, 1, [], 0.1, 0.3) == 1  # the same side whichever way the bar is read
    assert disulfide_side(a, b, 1, [(0.5, -1.0, 2.5, -0.1)], 0.1, 0.3) == 1  # a box below pushes it above
    assert disulfide_side(np.array([0.0, 0.0]), np.array([0.0, 3.0]), 1, [], 0.1, 0.3) == -1  # vertical: to the right


def test_overlapping_bars_are_drawn_as_brackets_that_keep_their_cysteines(monkeypatch):
    import foldmap.render as render

    path = DATA / "6ZS5.cif"
    pairs = load_links(path, load_backbone(path)).disulfides
    spots = {}
    for k, pair in enumerate(pairs):  # the first two lie on top of each other, the others far apart
        y = 0.0 if k < 2 else 50.0 * k
        spots[pair] = (np.array([1.0, y]), np.array([4.0 + k, y]))
    monkeypatch.setattr(render, "disulfide_points", lambda lay, loops, sses, chains, i, j, avoid=(): spots[(i, j)])
    # the lanes only: nothing hides these planted bars (what bodies and text hide is tested on real figures)
    monkeypatch.setattr(render, "visible_pieces", lambda points, bodies, halos=(): [np.asarray(points, float)])
    fig = make_figure(path)
    bars = {b.get_gid(): b.get_path().vertices for b in _gids(fig, "disulfide:")}
    first, second = (bars[f"disulfide:{i}-{j}"] for i, j in pairs[:2])
    assert len(first) == 2 and len(second) == 4
    a, b = spots[pairs[1]]
    assert np.allclose(second[0], a) and np.allclose(second[-1], b)  # still ends on its cysteines
    assert abs(second[1][1] - a[1]) >= 0.1 and np.isclose(second[1][1], second[2][1])  # a lane clear of the first bar
    assert all(len(v) == 2 for g, v in bars.items() if g not in {f"disulfide:{i}-{j}" for i, j in pairs[:2]})
    i, j = pairs[1]
    marks = [d for d in _gids(fig, "ss-dot:") + _gids(fig, "disulfide-end:") if d.get_gid().endswith(f":{i}-{j}")]
    marks = sorted(map(tuple, np.vstack([d.get_offsets() for d in marks])))  # a dot, or an end mark under a body
    assert np.allclose(marks, sorted([tuple(a), tuple(b)]))


def test_4ch9_bond_lying_on_another_detours_in_a_bracket():
    fig = make_figure(DATA / "4CH9.cif", assembly="asu")
    bars = {b.get_gid(): b.get_path().vertices for b in _gids(fig, "disulfide:")}
    assert [len(v) for _, v in sorted(bars.items())].count(4) == 1
    assert len(bars["disulfide:249-512"]) == 4  # nested under the longer bar, ends still on its cysteines
    assert any(d.get_gid() == "ss-dot:249-512" for d in _gids(fig, "ss-dot:"))
    from foldmap.render import _IN_PER_UNIT, disulfide_step

    ax = fig.axes[0]
    scale = fig.get_figwidth() / np.subtract(*ax.get_xlim()[::-1])  # inches per page unit
    lw = float(np.clip(0.9 * scale / _IN_PER_UNIT, 0.5, 1.2))
    step = disulfide_step(lw, scale * 72)
    v, straight = bars["disulfide:249-512"], bars["disulfide:226-535"]
    assert np.hypot(*(v[1] - v[0])) >= step - 1e-9  # the legs are at least a lane step long
    along = (straight[1] - straight[0]) / np.hypot(*(straight[1] - straight[0]))
    gap = abs((v[1] - straight[0]) @ np.array([-along[1], along[0]]))
    assert gap >= step - 1e-9  # the run clears the straight bar by a lane step


@pytest.mark.parametrize("inches_per_unit", [0.12, 0.3, 0.6])  # a large, a typical and a small figure
def test_the_lane_step_clears_a_sulfur_dot(inches_per_unit):
    from foldmap.render import _IN_PER_UNIT, disulfide_step

    ppu = inches_per_unit * 72
    lw = float(np.clip(0.9 * inches_per_unit / _IN_PER_UNIT, 0.5, 1.2))  # the line width draw() derives
    step = disulfide_step(lw, ppu)
    assert step > 5.5 * lw / ppu  # more than a dot's diameter
    assert step >= 3 * 2.4 * lw / ppu  # and three bar widths


def _drawn_points(fig, prefix):
    """Points along every drawn stretch of the lines whose gid starts with prefix, every 0.02 page units."""
    from foldmap.render import _dense

    out = []
    for a in fig.axes[0].patches:
        if (a.get_gid() or "").startswith(prefix):
            for poly in a.get_path().to_polygons(closed_only=False):
                out.append(_dense(poly)[0][1:-1])  # its ends meet the body edge they come out at
    return np.vstack(out) if out else np.empty((0, 2))


@pytest.mark.parametrize(
    "path",
    [DATA / "r100" / "7FAB-HL.cif.gz", DATA / "r100" / "1IGT-AB.cif.gz", DATA / "6ZS5.cif"],
    ids=["7FAB-HL", "1IGT-AB", "6ZS5"],
)
def test_disulfide_bars_cross_no_label(path):
    """A disulfide bar passes under the element bodies and the text (Q-R10): it is drawn only where no body or
    letter is, so no bar runs over a strand's letter (7FAB's C-I bar over F, E and D; 1IGT's bars over its strands'
    letters; 6ZS5's C-terminal bar through its "C" and "12 aa"). Judged by figcheck's ink at its line clearance,
    LINE_TEXT_CLEARANCE_PX (the Task 1 ruling for Task 6), and every bond still drawn."""
    from foldmap import figcheck
    from foldmap.render import link_bodies
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(path))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    hits = figcheck.text_collisions(d, clearance_px=figcheck.LINE_TEXT_CLEARANCE_PX)
    assert [h for h in hits if "disulfide" in h] == []
    bars = {a.get_gid() for a in d.figure.axes[0].patches if (a.get_gid() or "").startswith("disulfide:")}
    assert bars and len(bars) == len(d.layout.links.disulfides)
    pts = _drawn_points(d.figure, "disulfide:")
    assert not any(body.contains_points(pts).any() for body in link_bodies(d.layout))  # never over a body
