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


@pytest.mark.parametrize("entry, ligands", [("6JLN", 3), ("8ESZ", 1)])
def test_a_tether_to_another_chain_stays_within_reach_of_the_mark(entry, ligands):
    """Random-100 C2b-5 (A20). A target on another chain is kept within _TETHER_REACH of the home chain's main contact
    (_local_targets), but the mark goes near the mean of all its targets, which home contacts drawn far away pull off
    the main one: 6JLN's CLA B609 sat 10.15 from H's helix (3.80 from the main contact; its other B contact 14.91 off),
    b617 7.58 from h, c508 7.53 from k, and 8ESZ's 3PE 5604 6.94 from B2. Such a mark now takes the first spot, out to
    _SPOT_FAR rings, that keeps every tether to another chain within _TETHER_REACH; its home tethers stay drawn."""
    import numpy as np

    from foldmap.render import _TETHER_REACH, ligand_marks
    from foldmap.route import route_loops
    from helpers import r100_entry

    lay, sses, bb = make_layout(r100_entry(entry))
    marks = ligand_marks(lay, route_loops(lay, sses, bb), sses, Style())
    shared = 0
    for m in marks:
        chains = [lay.res_chain[r] for r in m["residues"]]
        c = np.array(m["centre"])
        far = [(ch, round(float(np.hypot(*(np.array(t) - c))), 2)) for t, ch in zip(m["targets"], chains)]
        far = [(ch, d) for ch, d in far if ch != chains[0] and d > _TETHER_REACH]
        shared += len(set(chains)) > 1
        assert not far, (m["ligand"].name, m["ligand"].chain, m["ligand"].seq, far)
    assert shared >= ligands


def test_a_mark_moved_for_another_chain_keeps_its_main_contact_within_reach():
    """Random-100 C2b-5 (review of the run): held only to its tethers to other chains, the spot search took 1SLV's Cu
    B1248 (chain B1) in the beta view (its first clear spot left A1's target 6.6 away) to the one text-free spot that
    kept A1 within _TETHER_REACH, 7.35 from its anchor and 7.67 from its main contact, its tethers back across the
    sheet. A spot moved for another chain also keeps the home chain's main contact (its heaviest target) within
    reach; with none such, the first spot stays."""
    import numpy as np

    from foldmap.render import _TETHER_REACH, ligand_marks
    from foldmap.route import route_loops
    from helpers import r100_entry

    lay, sses, bb = make_layout(r100_entry("1SLV"), names=["beta"])
    marks = ligand_marks(lay, route_loops(lay, sses, bb), sses, Style())
    cus = [m for m in marks if m["ligand"].name == "CU"]
    assert [m["ligand"].chain for m in cus] == ["B1", "B2"]  # one per copy of the assembly
    for cu in cus:
        assert len({lay.res_chain[r] for r in cu["residues"]}) > 1  # held by two chains
        main = cu["targets"][0]  # the heaviest home target comes first here (two residues on B, the loop 143)
        assert float(np.hypot(*np.subtract(main, cu["centre"]))) <= _TETHER_REACH, (cu["centre"], main)


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
    from foldmap.spec import draw_spec, figure_spec
    from test_links import _drawn_points, bodies_under

    d = draw_spec(figure_spec(path, **kw))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    hits = figcheck.text_collisions(d, clearance_px=figcheck.LINE_TEXT_CLEARANCE_PX)
    assert [h for h in hits if "ligand-tether" in h] == []
    pts = _drawn_points(d.figure, "ligand-tether:")
    assert len(pts) and bodies_under(d.figure, "ligand-tether:") == {}  # never over a body, as drawn


def test_tethers_run_over_fewer_bodies():
    """C9 (random-100): a mark's spot was the first clear one on rings round its anchor, whatever its straight
    tethers ran over (5Q42's Ni tether over strand B, 9YLM's Zn tethers across the hairpins, 7LX9's YGM). Where the
    first spot's tethers cross element bodies (other than the one each target lies on), the spot a few rings out
    whose tethers cross the least is taken, lengthening no tether past the first spot's reach. 1GOX's FMN, tethered
    to four elements of its barrel: 12.0 page units of tether over bodies before; every mark still drawn."""
    from foldmap.render import _over_bodies, ligand_marks
    from foldmap.spec import draw_spec, figure_spec

    spec = figure_spec(Path(__file__).parent / "data" / "r100" / "1GOX-A.cif.gz")
    d = draw_spec(spec)
    bodies = [p.rect for p in d.layout.placed.values()]
    marks = ligand_marks(d.layout, d.loops, d.sses, spec.look)
    assert len(marks) == 4
    over = sum(_over_bodies(m["centre"], m["targets"], bodies) for m in marks)
    assert over < 11.5, over
    assert _over_bodies((0.0, 0.0), [(1.0, 0.5)], [(0.5, 0.0, 1.5, 1.0)]) == 0.0  # its own target's body
    assert _over_bodies((0.0, 0.0), [(3.0, 0.5)], [(0.5, 0.0, 1.5, 1.0)]) > 0.9


def test_ligands_held_by_the_nucleic_acid_alone_are_marked_and_its_ions_said(tmp_path, capsys):
    """8XZE (random-100 C33): a THF-II riboswitch's SPM and THG touch only RNA, so io._ligands left them out (no
    protein residue holds them) and nothing said so. A small molecule held by nucleotides alone is now a ligand with
    those nucleotides (Ligand.nucleotides), marked and tethered to the stems and strands holding it, off the RNA's
    stems, paths and text; metal ions held by nucleotides alone are counted (Links.na_ions) and the note names them
    (a ribosome holds hundreds: 7Y7C's note names 204 of its 283 Mg, the ones nucleotides alone hold)."""
    from foldmap.cli import main
    from foldmap.render import ligand_marks
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "r100" / "8XZE.cif.gz"
    spec = figure_spec(path)
    d = draw_spec(spec)
    marks = {m["ligand"].name: m for m in ligand_marks(d.layout, d.loops, d.sses, spec.look)}
    assert sorted(marks) == ["SPM", "THG"]
    assert all(m["targets"] and not m["residues"] and m["ligand"].nucleotides for m in marks.values())
    stems = [s.rect for s in d.layout.stems]
    for m in marks.values():
        (cx, cy), (hw, hh) = m["centre"], m["half"]
        assert not any(cx - hw < x1 and x0 < cx + hw and cy - hh < y1 and y0 < cy + hh for x0, y0, x1, y1 in stems)
    assert d.layout.links.na_ions == {"MG": 1}
    assert main(["plot", str(path), "-o", str(tmp_path / "f.svg")]) == 0
    assert "not drawn: metal ions held by the nucleic acid alone (MG)" in capsys.readouterr().err


def test_a_ligand_held_by_a_stem_is_tethered_to_the_rail_facing_it():
    """8XZE (batch C2a review Q5): THG's tether to P1c ran from the stem's far rail along a G-C rung, reading as an
    extra rung. A tether to a stem ends on the touched nucleotide drawn nearest the mark, on the rail facing it."""
    import numpy as np

    from foldmap.render import _na_style, ligand_marks, nucleotide_points
    from foldmap.spec import draw_spec, figure_spec

    spec = figure_spec(DATA / "r100" / "8XZE.cif.gz")
    d = draw_spec(spec)
    lay = d.layout
    where = nucleotide_points(lay, _na_style(lay, spec.look))
    key = {
        (lab.chain, lab.seq, lab.icode): (s, i)
        for s, st in enumerate(lay.nucleic.strands)
        for i, lab in enumerate(st.labels)
    }
    stem = {k: st.id for st in lay.stems for k in st.axial}
    for m in ligand_marks(lay, d.loops, d.sses, spec.look):
        held = [key[n] for n in m["ligand"].nucleotides if n in key and key[n] in where]
        c = np.asarray(m["centre"], float)
        for t in m["targets"]:
            k = min(held, key=lambda h: np.hypot(*(np.asarray(where[h], float) - t)))  # the nucleotide it ends on
            same = [h for h in held if stem.get(h, h[0]) == stem.get(k, k[0])]
            gaps = [np.hypot(*(np.asarray(where[h], float) - c)) for h in same]
            assert np.hypot(*(np.asarray(t) - c)) <= min(gaps) + 1e-9
