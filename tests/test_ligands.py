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


def test_covalently_joined_ligand_residues_are_one_ligand():
    """25IN: each substrate is farnesyl phosphate bonded to GlcNAc-1-P, then Rha and three Galf - one molecule."""
    bb, ligs = _ligands("25IN")
    long = [g for g in ligs if "GZL" in g.codes]
    assert len(long) == 2 and all(g.codes == ["A1E8P", "A1E9H", "RAM", "GZL", "GZL", "GZL"] for g in long)
    assert all(g.symbol == "A1E8P·GlcNAc-P·Rha·Galf×3" for g in long)
    short = [g for g in ligs if g.codes == ["A1E8P", "GN1", "RAM"]]
    assert len(short) == 2 and len(ligs) == 4  # per transporter half: one long and one short lipid-linked glycan


def test_ligand_choice_by_any_of_its_codes():
    fig = make_figure(DATA / "25IN.cif", look=Style(ligands="GZL"))
    assert len(_markers(fig)) == 2


def test_glycosylated_entries_gain_no_spurious_ligands():
    """Sugars inside glycans stay glycans: the glycosylated fixtures list the same ligands as before (none)."""
    for name in ("6ZS5", "6ZYA", "8UTF", "8UUP"):
        path = DATA / f"{name}.cif"
        bb = load_backbone(path, "auto")
        links = load_links(path, bb, "auto")
        assert links.glycans and not links.ligands, name


def test_glycans_of_chains_left_out_are_not_ligands():
    """--chains drops the other chains' glycans with their chains; they never come back as ligands."""
    for name in ("6ZS5", "6ZYA", "8UTF", "8UUP"):
        path = DATA / f"{name}.cif"
        bb = load_backbone(path, "auto")
        first = bb.labels[0].chain
        links = load_links(path, bb.only([first]), "auto")
        assert not links.ligands, (name, [lg.symbol for lg in links.ligands])
    lay, _, _ = make_layout(DATA / "6ZYA.cif", chains=["A1"])
    assert not lay.links.ligands


def _home_of(lay, g):
    """The ligand's own chain when it touches it, else the chain it touches most (ties: alphabetically first)."""
    chains = [lay.res_chain[r] for r in g.contacts]
    if g.chain in chains:
        return g.chain
    return max(sorted(set(chains)), key=chains.count)


def test_ligand_tethers_stay_in_their_own_pocket():
    """1LDM tetramer: OXM between A1 and A3 was tethered across the figure, to an A3 residue 11.0 page units away.
    The home chain is always tethered (its targets are unbounded: true contact information); a target on another
    chain is kept only within reach of the pocket."""
    import numpy as np

    from foldmap.render import _SPOT_RADIUS, _TETHER_REACH, TETHER_CAP, ligand_marks
    from foldmap.route import route_loops
    from helpers import FID

    lay, sses, bb = make_layout(FID / "1LDM.cif.gz")
    marks = ligand_marks(lay, route_loops(lay, sses, bb), sses, Style())
    assert len(marks) >= 4  # NAD and OXM of every subunit
    assert TETHER_CAP == _TETHER_REACH + _SPOT_RADIUS  # the formula, see the comment at _TETHER_REACH
    assert TETHER_CAP < 11.0  # shorter than the cross-pocket tether it removed
    shared = [m for m in marks if len({lay.res_chain[r] for r in m["ligand"].contacts}) > 1]
    assert len(shared) == 4  # the OXMs between A1/A3 and A2/A4: the case the cap is about
    for m in marks:
        home = _home_of(lay, m["ligand"])
        chains = [lay.res_chain[r] for r in m["residues"]]
        assert len(m["residues"]) == len(m["targets"])
        assert home in chains
        c = np.array(m["centre"])
        for t, chain in zip(m["targets"], chains):
            if chain != home:
                assert float(np.hypot(*(np.array(t) - c))) <= TETHER_CAP, (m["ligand"].name, chain)


def test_a_ligand_whose_home_chain_is_not_drawn_is_tethered_to_the_drawn_one():
    """Focus on one protomer of 1LDM: the OXM of chain A3 touches A1 as well; A3 is not drawn, so it is tethered to
    A1 alone (_local_targets falls back to the drawn chain rather than dropping the ligand)."""
    from foldmap.render import ligand_marks
    from foldmap.route import route_loops
    from helpers import FID

    lay, sses, bb = make_layout(FID / "1LDM.cif.gz", focus="protomer")
    marks = ligand_marks(lay, route_loops(lay, sses, bb), sses, Style())
    a3 = [m for m in marks if m["ligand"].name == "OXM" and m["ligand"].chain == "A3"]
    assert len(a3) == 1 and a3[0]["targets"]
    assert {lay.res_chain[r] for r in a3[0]["residues"]} == {"A1"}


def test_haem_tethers_keep_to_the_pocket():
    """2HHB: every haem keeps its four tethers (its own chain's targets are never dropped)."""
    from foldmap.render import ligand_marks
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(DATA / "2HHB.cif")
    marks = ligand_marks(lay, route_loops(lay, sses, bb), sses, Style())
    assert [(m["ligand"].chain, len(m["targets"])) for m in marks] == [("A", 4), ("B", 4), ("C", 4), ("D", 4)]


@pytest.mark.parametrize(
    "path, kw",
    [
        (DATA / "fidelity" / "1MBN.cif.gz", {}),
        (DATA / "r100" / "1GOX-A.cif.gz", {}),
        (DATA / "fidelity" / "1LDM.cif.gz", {"focus": "A1"}),
    ],
    ids=["1MBN", "1GOX", "1LDM-focus-A1"],
)
def test_tethers_cross_no_label(path, kw):
    """A ligand's tethers pass under the element bodies and the text (Q-R10): drawn only where no body or letter
    is, and the ligand's mark takes the nearest clear spot (within 3 units) whose tethers cross no label, lengthening
    none past the plan-11 reach: 1MBN's heme tether over α7 (the "C" of the carried item is stale), 1GOX's FMN
    tethers on the βα barrel over a helix label (K8), 1LDM --focus A1's oxamate tethers over three labels. Judged
    by figcheck's ink at LINE_TEXT_CLEARANCE_PX."""
    from foldmap import figcheck
    from foldmap.render import link_bodies
    from foldmap.spec import draw_spec, figure_spec
    from test_links import _drawn_points

    d = draw_spec(figure_spec(path, **kw))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    hits = figcheck.text_collisions(d, clearance_px=figcheck.LINE_TEXT_CLEARANCE_PX)
    assert [h for h in hits if "ligand-tether" in h] == []
    pts = _drawn_points(d.figure, "ligand-tether:")
    assert len(pts) and not any(body.contains_points(pts).any() for body in link_bodies(d.layout))
