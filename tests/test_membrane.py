"""Membrane bands: the lipid bilayer drawn behind transmembrane proteins, with its inside and outside."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main, make_figure, make_layout
from foldmap.io import load_backbone
from foldmap.membrane import find_membrane

DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("name", ["1BL8", "1C3W"])
def test_membrane_proteins_are_recognised(name):
    m = find_membrane(DATA / f"{name}.cif", load_backbone(DATA / f"{name}.cif", "asu"))
    assert m is not None and m.source == "estimate" and 13 <= m.half <= 17
    assert np.isclose(np.linalg.norm(m.normal), 1.0)


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "2HHB", "8UUP", "8UTF", "6ZS5", "1TIM", "AF-P04637", "1EMA"])
def test_soluble_proteins_are_not(name):
    assert find_membrane(DATA / f"{name}.cif", load_backbone(DATA / f"{name}.cif", "asu")) is None


def test_a_soluble_tetramer_is_not_a_membrane():
    """1LDM (lactate dehydrogenase tetramer, the fidelity fixture): its long helices cross no slab oily enough."""
    path = DATA / "fidelity" / "1LDM.cif.gz"
    assert find_membrane(path, load_backbone(path, "auto")) is None


def test_a_slab_crossed_by_many_long_elements_may_be_more_polar():
    """6JO6 (photosystem I; random-100 C17, 6JO6 ruling): its slab is 0.174 polar (chlorophyll pockets, His ligands),
    over the 0.15 limit, so no membrane was found and the membrane view drew no band. A slab that 10 or more long
    elements cross may be up to 0.20 polar, the hydrophobic minimum (0.55) still holding. The cut keeps chains A and B
    within 18 Å of the mid-plane: 0.162 polar, 0.567 hydrophobic, 21 long elements crossing."""
    from foldmap import membrane

    path = DATA / "r100" / "6JO6-AB.cif.gz"
    bb = load_backbone(path, "auto")
    found = find_membrane(path, bb)
    assert found is not None and found.source == "estimate" and found.half == membrane.HALF
    assert membrane._MAX_POLAR < 0.174 <= membrane._MAX_POLAR_CROSSED and membrane._MANY_CROSSING == 10


def test_opm_dummy_atoms_give_the_exact_membrane(tmp_path):
    import gemmi

    st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    st.setup_entities()
    chain = gemmi.Chain("M")
    for k, (name, z) in enumerate((("N", 14.2), ("O", -14.2))):
        res = gemmi.Residue()
        res.name, res.seqid, res.het_flag = "DUM", gemmi.SeqId(900 + k, " "), "H"
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = name, gemmi.Element(name), gemmi.Position(0, 0, z)
        res.add_atom(atom)
        chain.add_residue(res)
    st[0].add_chain(chain)
    out = tmp_path / "opm.pdb"
    st.write_pdb(str(out))
    m = find_membrane(out, load_backbone(out, "asu"))
    assert m.source == "opm" and m.half == pytest.approx(14.2) and abs(m.normal[2]) == pytest.approx(1.0)


def _band(fig):
    return [a for a in fig.axes[0].patches if (a.get_gid() or "") == "membrane"]


def test_band_drawn_across_the_transmembrane_helices():
    lay, sses, _ = make_layout(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    assert lay.membrane is not None
    y0, y1 = lay.membrane["y"]
    assert y1 - y0 > 3.0
    long = [p for p in lay.placed.values() if p.sse.kind == "H" and len(p.sse) >= 18]
    crossing = [p for p in long if p.rect[1] < y0 + 1.0 and p.rect[3] > y1 - 1.0]
    assert len(crossing) >= 5  # bacteriorhodopsin's seven TM helices run through the band
    fig = make_figure(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    assert len(_band(fig)) == 1
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert {"in", "out"} <= texts


def test_positive_inside_rule_puts_bacteriorhodopsin_c_terminus_inside():
    from foldmap.layout import termini

    lay, sses, _ = make_layout(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    y0, y1 = lay.membrane["y"]
    ends = {end: port for end, _, port, _ in termini(lay, sses)}
    inside_low = lay.membrane["inside"] == "below"
    c_y, n_y = ends["C"][1], ends["N"][1]
    # C in the cytoplasm; N on the outer side (helix A starts at residue 10, so its end sits just inside the band)
    assert (c_y < y0) == inside_low and (n_y > (y0 + y1) / 2) == inside_low


def test_soluble_protein_with_membrane_auto_draws_nothing(tmp_path, capsys):
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(tmp_path / "x.svg"), "--membrane", "auto"]) == 0
    assert "membrane" in capsys.readouterr().err.lower()
    assert 'id="membrane"' not in (tmp_path / "x.svg").read_text()


def test_sides_labelled_only_when_the_positive_inside_rule_is_clear():
    fig = make_figure(DATA / "1BL8.cif", membrane="auto", assembly="asu")  # truncated KcsA: too close to call
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert len(_band(fig)) == 1 and not ({"in", "out"} & texts)


def test_outside_on_top_even_for_symmetric_assemblies():
    lay, _, _ = make_layout(DATA / "1C3W.cif", membrane="auto")  # bacteriorhodopsin trimer, unrolled about C3
    assert lay.membrane["sides"] and lay.membrane["inside"] == "below"


def test_trimer_protomers_are_drawn_side_by_side():
    lay, _, _ = make_layout(DATA / "1C3W.cif", membrane="auto")
    span = {}
    for p in lay.placed.values():
        lo, hi = span.get(p.sse.chain, (p.cx, p.cx))
        span[p.sse.chain] = (min(lo, p.cx), max(hi, p.cx))
    (a0, a1), (b0, b1), (c0, c1) = sorted(span.values())
    assert len(span) == 3 and a1 < b0 and b1 < c0  # each chain its own x-range, no interleaving


OPM = DATA / "opm" / "1BL8_opm.pdb"


def test_opm_file_gives_the_sides():
    """OPM marks the outer boundary with O dummies, the inner (cytoplasmic) one with N: KcsA's selectivity
    filter faces out, its N-terminus in."""
    bb = load_backbone(OPM, "asu")
    m = find_membrane(OPM, bb)
    assert m.source == "opm" and m.sided and m.half == pytest.approx(16.2, abs=0.05)
    filter77 = next(k for k, l in enumerate(bb.labels) if l.chain == "A" and l.seq == 77)
    n_term = next(k for k, l in enumerate(bb.labels) if l.chain == "A")
    assert m.depth(bb.ca[filter77]) == pytest.approx(12.8, abs=0.5)  # was -12.8: re-derived from Lys/Arg
    assert m.depth(bb.ca[n_term]) < -15


def test_dummy_atoms_without_sides_fall_back_to_positive_inside(tmp_path):
    """Review focus 4: dummies all of one name give the planes, not the sides."""
    import gemmi

    st = gemmi.read_structure(str(OPM))
    for chain in st[0]:
        for res in chain:
            if res.name == "DUM":
                for atom in res:
                    atom.name = "C"
    out = tmp_path / "one_name.pdb"
    st.write_pdb(str(out))
    m = find_membrane(out, load_backbone(out, "asu"))
    assert m.source == "opm" and m.half == pytest.approx(16.2, abs=0.05)
    assert abs(m.normal[2]) == pytest.approx(1.0)


def test_membrane_flip_turns_the_figure_over_without_mirroring():
    auto, _, _ = make_layout(DATA / "1C3W.cif", membrane="auto", assembly="asu")
    flip, _, _ = make_layout(DATA / "1C3W.cif", membrane="flip", assembly="asu")
    a, f = np.array(auto.page_axes), np.array(flip.page_axes)
    assert np.allclose(f[0], a[0], atol=1e-6)  # page-right kept
    assert np.allclose(f[1], -a[1], atol=1e-6) and np.allclose(f[2], -a[2], atol=1e-6)  # up and depth turned over
    assert np.linalg.det(f) > 0  # a rotation, not a mirror
    assert flip.membrane["inside"] == auto.membrane["inside"]  # the face now called "outside" is drawn on top


def test_membrane_flip_turns_an_unrolled_channel_over():
    auto, sses, _ = make_layout(DATA / "1BL8.cif", membrane="auto")  # KcsA: C4 axis along the membrane normal
    flip, _, _ = make_layout(DATA / "1BL8.cif", membrane="flip")
    top, low = "A:86-112", "A:27-51"  # the pore helix pair either side of the selectivity filter
    assert auto.placed[top].cy > auto.placed[low].cy  # the unrolled channel stands with this helix on top
    assert flip.placed[top].cy < flip.placed[low].cy  # and is turned over
    assert flip.membrane["y"] != auto.membrane["y"]


def test_membrane_flip_round_trips_through_a_layout_file(tmp_path):
    saved = tmp_path / "fig.yaml"
    assert (
        main(
            [
                "plot",
                str(DATA / "1C3W.cif"),
                "--membrane",
                "flip",
                "-o",
                str(tmp_path / "a.svg"),
                "--save-layout",
                str(saved),
            ]
        )
        == 0
    )
    assert "membrane: flip" in saved.read_text(encoding="utf-8")
    assert main(["plot", str(DATA / "1C3W.cif"), "--layout-file", str(saved), "-o", str(tmp_path / "b.svg")]) == 0
    assert (tmp_path / "a.svg").read_bytes() == (tmp_path / "b.svg").read_bytes()


def _row_band(lay, p):
    """The band of the stack row an element sits in: the one nearest its middle."""
    mid = (p.rect[1] + p.rect[3]) / 2
    return min(lay.membrane["bands"], key=lambda b: abs((b[0] + b[1]) / 2 - mid))


@pytest.mark.parametrize("assembly", ["asu", "auto"])
def test_stack_mode_keeps_every_tm_helix_in_the_bilayer(assembly):
    """1C3W --mode stack --membrane auto drew the band at one height while the rows wrapped: helix G sat below it."""
    lay, sses, bb = make_layout(DATA / "1C3W.cif", mode="stack", membrane="auto", assembly=assembly)
    tm = [s for s in sses if s.kind == "H" and len(s) >= 15]  # the seven TM helices of each chain
    assert len(tm) == 7 * len({s.chain for s in sses})
    for s in tm:
        p = lay.placed[s.id]
        y0, y1 = _row_band(lay, p)  # each protomer (36.5 wide) has a row and a band of its own
        lo, hi = sorted((p.n_port[1], p.c_port[1]))
        assert min(hi, y1) - max(lo, y0) >= 0.75 * (y1 - y0), s.id  # measured worst: α3, 0.83


@pytest.mark.parametrize(
    "path, chains, rows",
    [(DATA / "1C3W.cif", ["A1", "A2"], [["A1"], ["A2"]]), (OPM, None, [["A", "B"], ["C", "D"]])],
    ids=["1C3W-two-chains", "1BL8-opm"],
)
def test_a_wide_depth_keeping_stack_packs_whole_chains_into_rows_each_with_its_band(path, chains, rows, monkeypatch):
    """1PV7 --mode stack --membrane auto put 21 helices in one row 2160 px wide. A row wider than _MAX_ROW wraps:
    whole chains packed into rows in chain order up to _MAX_ROW (a chain wider than that alone in its row), each
    row with its own band, and each element at its depth: a whole row moves as one. KcsA: two rows of two."""
    from foldmap import layout as layout_mod

    lay, sses, _ = make_layout(path, mode="stack", membrane="auto", chains=chains)
    by_row = [sorted({lay.placed[k].sse.chain for k in ids if k in lay.placed}) for ids in lay.stack_rows]
    assert by_row == rows
    assert len(lay.membrane["bands"]) == len(rows) and lay.membrane["y"] == lay.membrane["bands"][0]
    span = [
        (min(lay.placed[k].rect[1] for k in ids), max(lay.placed[k].rect[3] for k in ids), ids)
        for ids in lay.stack_rows
    ]
    assert all(lo > hi for (lo, _, _), (_, hi, _) in zip(span, span[1:]))  # rows top to bottom, apart
    for _, _, ids in span:
        x0 = min(lay.placed[k].rect[0] for k in ids)
        x1 = max(lay.placed[k].rect[2] for k in ids)
        assert x1 - x0 <= layout_mod._MAX_ROW or len({lay.placed[k].sse.chain for k in ids}) == 1
    for p in lay.placed.values():
        if p.sse.kind == "H" and len(p.sse) >= 20:  # TM helices cross their own row's band
            y0, y1 = _row_band(lay, p)
            lo, hi = sorted((p.n_port[1], p.c_port[1]))
            assert min(hi, y1) - max(lo, y0) >= 0.75 * (y1 - y0), p.sse.id
    fig = make_figure(path, mode="stack", membrane="auto", chains=chains)
    drawn = [a for a in fig.axes[0].patches if (a.get_gid() or "").split(":")[0] == "membrane"]
    assert len(drawn) == len(rows)
    monkeypatch.setattr(layout_mod, "_MAX_ROW", 1e9)  # the same figure in one row: every element's depth
    one, _, _ = make_layout(path, mode="stack", membrane="auto", chains=chains)
    for ids in lay.stack_rows:
        shift = {round(lay.placed[k].cy - one.placed[k].cy, 6) for k in ids if k in lay.placed}
        assert len(shift) == 1, ids


def test_a_narrow_depth_keeping_stack_stays_one_row():
    lay, sses, _ = make_layout(OPM, mode="stack", membrane="auto", chains=["A", "B"])
    assert lay.width <= 30.0 and len(lay.membrane["bands"]) == 1


@pytest.mark.parametrize("kw", [{}, {"membrane": "auto"}, {"mode": "stack", "membrane": "auto"}])
def test_opm_dummy_atoms_are_never_drawn_as_ligands(kw):
    """1BL8 from OPM drew 77 'DUM' ligand boxes: the membrane boundary markers are not molecules. They still give
    the membrane (sides included); KcsA's K+ ions are still drawn."""
    lay, _, _ = make_layout(OPM, **kw)
    assert not any("DUM" in g.codes for g in lay.links.ligands)
    assert any(g.metal and g.symbol == "K" for g in lay.links.ligands)
    texts = [t.get_text() for t in make_figure(OPM, **kw).axes[0].texts]
    assert "DUM" not in texts and "K" in texts
    if kw:
        assert lay.membrane["source"] == "opm" and lay.membrane["sides"]


# ------------------------------------------------------------------ outer-membrane beta-barrels (C-MEM-BARREL)

R100 = DATA / "r100"
OMP_NORMALS = {  # the membrane normal OPM gives, in each fixture's frame (OPM tilts these barrels 5 degrees)
    "9FWZ-DB": (-0.0193, -0.2999, 0.9538),  # OPM's 3RFZ FimD superposed on 9FWZ chain D (CA rmsd 1.7 A, 678 residues)
    "1BXW-A": (-0.9302, 0.3594, -0.0749),  # OPM's 1bxw (OmpA)
}


def _closed_barrel_strands(bb):
    """The strands (>= 7 residues) of the fixture's closed beta-barrel, as foldmap's own sheets find it."""
    from foldmap.dssp import assign_dssp
    from foldmap.sheets import build_sheets
    from foldmap.ss import build_sses

    dssp = assign_dssp(bb)
    sheets = build_sheets(build_sses(bb, dssp.ss, bridges=dssp.bridges), dssp.bridges)
    (barrel,) = [sh for sh in sheets if sh.closed and len(sh.strands) >= 8]
    return [s for s in barrel.strands if len(s) >= 7]


@pytest.mark.parametrize("name", OMP_NORMALS)
def test_outer_membrane_barrels_are_recognised(name):
    """9FWZ (the FimD usher, author: "9FWZ does not draw membrane") and OmpA: the barrel's axis is the normal, the
    band is 24 A thick (an outer membrane), and every barrel strand runs through it."""
    path = R100 / f"{name}.cif.gz"
    bb = load_backbone(path, "asu")
    m = find_membrane(path, bb)
    assert m is not None and m.source == "estimate" and m.half == 12.0
    assert np.linalg.norm(m.normal) == pytest.approx(1.0, abs=1e-3)
    opm = np.array(OMP_NORMALS[name]) / np.linalg.norm(OMP_NORMALS[name])
    assert np.degrees(np.arccos(abs(m.normal @ opm) / np.linalg.norm(m.normal))) < 10.0
    strands = _closed_barrel_strands(bb)
    assert len(strands) >= 8
    for s in strands:
        d = m.depth(bb.ca[s.start : s.end + 1])
        through = min(d.max(), m.half) - max(d.min(), -m.half)
        assert through >= m.half, s.id  # at least half the slab's thickness inside it


def test_outer_membrane_barrel_normal_and_centre_are_on_the_grid():
    """Snapped like every other layout input, so the page frame built from it never carries BLAS noise."""
    path = R100 / "1BXW-A.cif.gz"
    m = find_membrane(path, load_backbone(path, "asu"))
    assert np.array_equal(m.normal, np.round(m.normal * 1e4) / 1e4)
    assert m.centre == np.round(m.centre * 1e4) / 1e4


@pytest.mark.parametrize("path", [DATA / "1EMA.cif", DATA / "1TIM.cif", DATA / "fidelity" / "1LDM.cif.gz"])
def test_soluble_barrels_are_not_membranes(path):
    """GFP is a closed 11-strand barrel of long strands, but its outside is polar: no membrane."""
    assert find_membrane(path, load_backbone(path, "asu")) is None


def test_gfp_barrel_is_rejected_by_its_outward_residues():
    from foldmap.membrane import _MIN_OUTWARD_HYDROPHOBIC, _barrel_slab

    bb = load_backbone(DATA / "1EMA.cif", "asu")
    _, _, hydrophobic, charged = _barrel_slab(bb, _closed_barrel_strands(bb))
    assert hydrophobic < _MIN_OUTWARD_HYDROPHOBIC  # measured 0.35 (outer-membrane barrels 0.76-1.00)


HELICAL = {  # find_membrane before the beta-barrel estimate: (normal, centre, half, source, sided), or None
    DATA / "1C3W.cif": (
        [0.015091081379981436, 0.11620109627440002, -0.9931110534514308],
        -1.5448989290611692,
        15.0,
        "estimate",
        True,
    ),
    DATA / "1BL8.cif": (
        [0.36292304637025796, -0.0033065376136334106, 0.9318132480397247],
        47.81279573824658,
        15.0,
        "estimate",
        False,
    ),
    OPM: ([0.0, 0.0, 1.0], 0.0, 16.200000000000003, "opm", True),
    DATA / "fidelity" / "1PV7.cif.gz": None,  # a 46-residue window of LacY: too little to call
}


@pytest.mark.parametrize("path", HELICAL, ids=lambda p: p.name)
def test_helical_membranes_are_unchanged_by_the_barrel_estimate(path):
    m = find_membrane(path, load_backbone(path, "asu"))
    if HELICAL[path] is None:
        assert m is None
        return
    normal, centre, half, source, sided = HELICAL[path]
    assert m.normal.tolist() == pytest.approx(normal, rel=1e-12, abs=1e-12)
    assert m.centre == pytest.approx(centre, rel=1e-12, abs=1e-12)
    assert (m.half, m.source, m.sided) == (half, source, sided)


@pytest.mark.parametrize("mode", ["projected", "stack"])
def test_outer_membrane_barrel_strands_cross_the_drawn_band(mode):
    """9FWZ --membrane auto: the band is drawn, and every drawn strand of the FimD barrel crosses it."""
    path = R100 / "9FWZ-DB.cif.gz"
    lay, sses, bb = make_layout(path, mode=mode, membrane="auto")
    assert lay.membrane is not None
    strands = {s.id for s in _closed_barrel_strands(bb)}
    drawn = [p for p in lay.placed.values() if p.sse.id in strands]
    assert len(drawn) >= 8
    for p in drawn:
        y0, y1 = _row_band(lay, p)
        lo, hi = sorted((p.n_port[1], p.c_port[1]))
        assert min(hi, y1) - max(lo, y0) >= 0.5 * (y1 - y0), p.sse.id
    fig = make_figure(path, mode=mode, membrane="auto")
    assert [a for a in fig.axes[0].patches if (a.get_gid() or "").split(":")[0].startswith("membrane")]


@pytest.mark.parametrize("name, sided", [("1BXW-A", True), ("9FWZ-DB", False)])
def test_outer_membrane_barrel_sides_come_from_its_loops(name, sided):
    """The positive-inside rule called OmpA's sides backwards (measured on 19 outer-membrane barrels against OPM: wrong
    on 16, 10 of them called clear). An outer-membrane barrel's extracellular loops are long and its periplasmic
    turns short: OmpA's loops (51 residues against 10) put the outside where OPM has it; FimD's (75 against 71, its
    plug capped) are too close to call, so no in/out labels, though the long-loop face still goes on top."""
    path = R100 / f"{name}.cif.gz"
    m = find_membrane(path, load_backbone(path, "asu"))
    assert m.sided == sided
    assert m.normal @ np.array(OMP_NORMALS[name]) > 0  # the outside where OPM has it


def test_band_lies_over_the_sheet_plates_and_under_the_loops():
    """9FWZ: the band hid under the barrel's sheet plate, showing only beyond its edges. It now lies over every plate
    (translucent, so the plate still shows through) and under the loops and bodies."""
    from matplotlib.colors import to_hex

    from foldmap.render import Z

    assert Z["sheet-panel"] < Z["membrane"] < Z["membrane-heads"] < Z["loop"] < Z["body"]
    ax = make_figure(R100 / "9FWZ-DB.cif.gz", membrane="auto").axes[0]
    (band,) = _band(ax.figure)
    assert band.get_zorder() == Z["membrane"] and 0 < band.get_alpha() < 1
    assert to_hex(band.get_facecolor()) == "#f7f0de"
    bx0, by0, bx1, by1 = band.get_bbox().extents
    plates = [p for p in ax.patches if (p.get_gid() or "").startswith("sheet-panel:")]
    crossed = [p for p in plates if p.get_extents().overlaps(band.get_extents())]
    assert crossed  # the unrolled barrel's plate
    for p in crossed:
        px0, _, px1, _ = p.get_bbox().extents
        assert band.get_zorder() > p.get_zorder() and bx0 <= px0 and px1 <= bx1  # over the plate, all across it
    heads = [c for c in ax.collections if (c.get_gid() or "").startswith("membrane-heads:")]
    assert len(heads) == 2 and all(c.get_zorder() == Z["membrane-heads"] for c in heads)


def test_stack_rows_get_a_band_only_from_the_elements_that_define_the_membrane():
    """9FWZ --mode stack drew two bands: one on FimD's row and one on the row of the FimA subunit whose strands lie in
    the usher's lumen and so span the slab's depth. Only the barrel's own strands give a row its band."""
    lay, sses, bb = make_layout(R100 / "9FWZ-DB.cif.gz", mode="stack", membrane="auto")
    rows = [{lay.placed[k].sse.chain for k in ids if k in lay.placed} for ids in lay.stack_rows]
    assert rows == [{"B"}, {"D"}]  # FimA and FimD each in a row of their own
    assert len(lay.membrane["bands"]) == 1 and lay.membrane["y"] == lay.membrane["bands"][0]
    y0, y1 = lay.membrane["y"]
    usher = [lay.placed[k] for k in lay.stack_rows[1] if k in lay.placed]
    assert min(p.rect[1] for p in usher) < y0 < y1 < max(p.rect[3] for p in usher)  # the band is on FimD's row
    fig = make_figure(R100 / "9FWZ-DB.cif.gz", mode="stack", membrane="auto")
    assert len([a for a in fig.axes[0].patches if (a.get_gid() or "").split(":")[0] == "membrane"]) == 1


def test_every_transmembrane_helix_crosses_the_band_in_the_membrane_view():
    """1BL8 (random-100 C17, A18 on 6JLN and 5LDW): the projected membrane view placed KcsA's inner helices (86-112)
    by their loops and contacts, with half of each above the band it drew (0.49 of the band covered, no face
    crossed past the middle). When a transmembrane element comes out off the band, the view is laid out again with
    every transmembrane element at its projected height (moved to meet a face when drawn wholly inside the slab),
    and the band is the membrane's faces as projected: each enters through a face and reaches past the middle."""
    path = DATA / "1BL8.cif"
    lay, sses, bb = make_layout(path, membrane="auto")
    mem = find_membrane(path, bb)
    depth = mem.depth(bb.ca)
    lo, hi = lay.membrane["y"]
    tm = [
        s
        for s in sses
        if depth[s.start : s.end + 1].min() <= -mem.half + 2 and depth[s.start : s.end + 1].max() >= mem.half - 2
    ]
    assert len(tm) == 8
    for s in tm:
        p = lay.placed[s.id]
        y0, y1 = sorted((p.n_port[1], p.c_port[1]))
        assert y0 < (lo + hi) / 2 < y1 and (y0 <= lo or y1 >= hi), s.id


def test_a_loop_runs_in_the_slab_only_beside_its_own_port():
    """A18's loops in the slab (figcheck.slab_runs; C26 port rule, membrane slab ruling): with a band from 0 to 10, a
    loop leaving a port inside it may run sideways on its first two segments (the stub out of the port and the turn
    after it), and on its first run as far as that stays over its own element's body (8W9Y: back above the
    interfacial helix it leaves, round α8's end). A run past the element, or any later run, is in the slab."""
    from types import SimpleNamespace

    from foldmap.figcheck import slab_runs
    from foldmap.route import Loop

    helix = SimpleNamespace(rect=(-6.0, 3.5, 0.0, 4.5))  # lying in the slab, its C end (the port) at x = 0
    placed = {"h": helix, "o": SimpleNamespace(rect=(9.0, -5.0, 10.0, 5.0))}
    band = [(0.0, 10.0)]

    def loop(*pts):
        return Loop("h", "o", list(pts), False, False)

    turn = loop((0, 4), (0, 6), (3, 6), (3, 12), (9.5, 12), (9.5, 5))  # up, then sideways on its 2nd segment
    back = loop((0, 4), (0.3, 4), (0.3, 6), (-2, 6), (-2, 12), (9.5, 12), (9.5, 5))  # 3rd segment, over the helix
    past = loop((0, 4), (0.3, 4), (0.3, 6), (-7, 6), (-7, 12), (9.5, 12), (9.5, 5))  # 3rd segment, past the helix
    later = loop((0, 4), (0, 6), (2, 6), (2, 7), (5, 7), (5, 12), (9.5, 12), (9.5, 5))  # a 2nd run in the slab
    assert slab_runs([turn, back], band, placed) == []
    assert slab_runs([past], band, placed) == ["h>o"] and slab_runs([later], band, placed) == ["h>o"]


def test_the_second_membrane_layout_is_kept_only_when_it_does_better(monkeypatch):
    """C17 fix round 1: 6JLN's second layout put every TM helix on the band but added slab loops (1 -> 3); it is now
    kept only when it has fewer helices off the band and no more slab loops or route fallbacks (cli._membrane_score,
    on a routed copy). 1BL8's does better (4 -> 0, nothing added) and is kept; scored worse, the first is kept."""
    import foldmap.cli as cli

    path = DATA / "1BL8.cif"
    kept, _, _ = make_layout(path, membrane="auto")
    scores = iter([(4, 0, 0), (0, 1, 0)])
    monkeypatch.setattr(cli, "_membrane_score", lambda *a: next(scores))
    first, _, _ = make_layout(path, membrane="auto")
    assert kept.membrane["y"] != first.membrane["y"]  # the faces as projected, against the medians of the first
    assert next(scores, None) is None


@pytest.mark.parametrize("entry", ["5LDW", "6JO6", "6JLN"])
def test_a_route_turning_out_of_a_port_at_or_in_the_slab_turns_at_the_port(entry):
    """Random-100 C2b-5 (A18). 6JO6: 5:167-169's C port lies 0.03 below the lower face and points up; its route's
    first cell, 0.4 out, lay 0.42 inside the face, and the turn round the helix's end ran there. 5LDW: in C17's second
    layout L:448-454's C port, inside the band, points up into a 0.81 gap under L:369-380; the 0.4 step met the row
    of the loop entering L:369-380, and the route jogged down to run 1.0 sideways in mid-slab, one loop in the slab
    (scored (0, 1, 0) against the first layout's (6, 0, 0), so C17 kept six TM helices off the band). A port whose
    first step lands in the slab now starts nearer (route._Grid._into_slab): 5LDW keeps the second layout, every TM
    helix across the band, and neither view has a loop in the slab. 6JLN: c:81-83's N port, inside the band, leaves
    sideways; starting it nearer moved its route's turn beside the element, 0.75 along mid-slab, so a port leaving
    sideways starts as before."""
    from foldmap.figcheck import membrane_view
    from foldmap.route import route_loops
    from helpers import r100_entry

    path = r100_entry(entry)
    lay, sses, bb = make_layout(path, membrane="auto")
    view = membrane_view(lay, sses, bb, find_membrane(path, bb), route_loops(lay, sses, bb))
    assert view["tm_helices"] >= 14 and view["tm_short"] == [] and view["loops_in_slab"] == [], view


def test_a_users_move_never_pushes_a_transmembrane_helix_off_its_height():
    """C17 fix round 1 (Q8): the transmembrane elements of a membrane view laid out at their projected height keep it
    through the user's own moves (_adjust's _separate) and the pushes off unpaired strands (_off_walls). In 1BL8, the
    pore helix moved onto the top of an inner helix pushed that helix 1.24 down."""
    path = DATA / "1BL8.cif"
    plain, _, _ = make_layout(path, membrane="auto")
    moved, _, _ = make_layout(path, membrane="auto", move=[("res:A:62-74", (2.64, 10.6))])
    assert moved.placed["A:86-112"].cy == plain.placed["A:86-112"].cy
    assert moved.placed["A:62-74"].cy > plain.placed["A:62-74"].cy + 10


def test_fimd_band_survives_its_pendant_strand_leaving_the_barrel():
    """Plan 12 Task 2: D:662-666, off FimD's 24-strand ring, is no barrel strand (membrane._barrels reads the track
    only). The band is still found, its normal within 10° of the one found before (fbcee15)."""
    path = R100 / "9FWZ-DB.cif.gz"
    found = find_membrane(path, load_backbone(path, "auto"))
    assert found is not None and found.source == "estimate"
    before = np.array([-0.0211, -0.225, 0.9741])
    assert np.degrees(np.arccos(abs(found.normal @ before) / np.linalg.norm(before))) <= 10.0


def test_a_slab_many_long_elements_cross_but_not_oily_enough_is_refused():
    """Plan 10 review (C17): the 0.20 polar allowance for a slab 10 or more long elements cross never waives the
    0.55 hydrophobic minimum. 6JO6-AB (21 crossing, 0.567 hydrophobic) with enough of its slab's oily residues made
    Ser (neither oily nor polar) to fall under 0.55 finds no membrane."""
    import dataclasses

    from foldmap import membrane

    path = DATA / "r100" / "6JO6-AB.cif.gz"
    bb = load_backbone(path, "auto")
    assert membrane._estimate(bb, membrane.build_sses(bb, membrane.assign_dssp(bb).ss)) is not None
    oily = [k for k, l in enumerate(bb.labels) if l.name in membrane._HYDROPHOBIC]
    labels = list(bb.labels)
    for k in oily[::8]:  # one in eight oily residues
        labels[k] = labels[k]._replace(name="SER")
    thin = dataclasses.replace(bb, labels=labels)
    sses = membrane.build_sses(thin, membrane.assign_dssp(thin).ss)
    long = [s for s in sses if (s.kind == "H" and len(s) >= 16) or (s.kind == "E" and len(s) >= 8)]
    assert len(long) >= membrane._MANY_CROSSING
    assert membrane._estimate(thin, sses) is None
