"""Membrane bands: the lipid bilayer drawn behind transmembrane proteins, with its inside and outside."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main, make_figure, make_layout
from foldmap.io import load_backbone
from foldmap.membrane import find_membrane
from helpers import R100_CACHE

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


class _Element:  # an element as membrane._estimate reads it: kind, length, axis, residue span and chain
    def __init__(self, kind, start, end, chain):
        self.kind, self.start, self.end, self.chain, self.axis = kind, start, end, chain, np.array([0.0, 0.0, 1.0])

    def __len__(self):
        return self.end - self.start + 1


def _synthetic(helices, others):
    """A backbone (CA and residue labels only) and its elements: `helices` are (chain, x, z0, z1, name) runs along z,
    one CA every 1.5 Å named name(z); `others` are (chain, xyz, name) residues in no element."""
    from types import SimpleNamespace

    ca, labels, elements = [], [], []
    for chain, x, z0, z1, name in helices:
        start = len(ca)
        for z in np.arange(z0, z1 + 1e-9, 1.5):
            ca.append((x, 0.0, z))
            labels.append(SimpleNamespace(chain=chain, name=name(z)))
        elements.append(_Element("H", start, len(ca) - 1, chain))
    for chain, xyz, name in others:
        ca.append(xyz)
        labels.append(SimpleNamespace(chain=chain, name=name))
    return SimpleNamespace(ca=np.array(ca, float), labels=labels), elements


def test_a_slab_no_element_spans_moves_to_one_the_elements_span():
    """random-100 held-out 6MJ2 (TRPM2): the summed score put the slab 4-5 Å into the dense cytoplasmic domain,
    where no S1-S4 helix reaches both faces, so the membrane view found none to draw its band by. Four helices along
    z from 0 to 31.5 Å and an oily blob just below them (z = -6): the best-scoring slab takes in the blob (centre 8),
    where no helix spans it (nor crosses it: refused); the best-scoring slab the helices span (centre 14) is the
    membrane."""
    from foldmap import membrane

    helices = [("A", x, 0.0, 31.5, lambda z: "LEU") for x in (0.0, 10.0, 20.0, 30.0)]
    blob = [("A", (40.0 + k % 8, float(k // 8), -6.0), "ALA") for k in range(60)]
    bb, sses = _synthetic(helices, blob)
    m = membrane._estimate(bb, sses)
    assert m is not None and m.half == membrane.HALF
    for s in sses:
        d = m.depth(bb.ca[s.start : s.end + 1])
        assert d.min() <= -m.half + membrane.FACE_REACH and d.max() >= m.half - membrane.FACE_REACH


def test_a_slab_is_judged_by_the_chains_that_cross_it():
    """random-100 held-out 7VVK (PTH1R-Gs): the hormone in the receptor's pocket and Gs's α5 in its cytoplasmic cavity
    made the slab 0.189 polar, over the 0.15 limit, and no membrane was found. A chain with no element across the slab
    is a partner reaching into it, not lipid-facing: four helices oily in the slab and charged beyond it (chain A) and
    a charged peptide lying in the mid-plane (chain P, no element) are a membrane, judged by chain A's residues."""
    from foldmap import membrane

    def name(z):
        return "LYS" if abs(z) > membrane.HALF else ("GLY" if round(z / 1.5) % 5 == 0 else "LEU")

    helices = [("A", x, -25.0, 25.0, name) for x in (0.0, 10.0, 20.0, 30.0)]
    peptide = [("P", (5.0 + k, 5.0, 0.0), "LYS" if k % 2 else "GLU") for k in range(24)]
    bb, sses = _synthetic(helices, peptide)
    polar = np.mean([l.name in membrane._POLAR for l, z in zip(bb.labels, bb.ca[:, 2]) if abs(z) < membrane.HALF])
    assert polar > membrane._MAX_POLAR  # every residue in the slab counted: refused
    m = membrane._estimate(bb, sses)
    assert m is not None and abs(m.centre) <= 1.0


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


def test_opm_dummy_atoms_after_end_give_the_membrane(tmp_path):
    """OPM files often put the DUM records after END; they give the planes, and nothing else after END is read."""
    import gemmi

    st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    out = tmp_path / "opm_after_end.pdb"
    st.write_pdb(str(out))
    body = [ln for ln in out.read_text().splitlines() if ln[:6].rstrip() not in ("END", "MASTER")]
    dum = [
        f"HETATM{9000 + k:5d} {name:<4s} DUM  {900 + k:4d}    {0:8.3f}{0:8.3f}{z:8.3f}  1.00  0.00"
        for k, (name, z) in enumerate((("N", 14.2), ("O", -14.2)))
    ]
    stray = "ATOM   9100  CA  ALA A 999       0.000   0.000   0.000  1.00  0.00           C"
    out.write_text("\n".join(body + ["END"] + dum + [stray]) + "\n")
    bb = load_backbone(out, "asu")
    assert len(bb.labels) == len(load_backbone(DATA / "1UBQ.cif", "asu").labels)  # the stray record is not in the model
    m = find_membrane(out, bb)
    assert m.source == "opm" and m.half == pytest.approx(14.2) and m.sided


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


def test_a_band_is_drawn_by_the_line_of_helices_that_stop_short_of_a_face():
    """Plan 26 (32KS: TM helices of 15-19 residues in a 30 Å slab; random-100 held-out 6MJ2, whose S1-S4 helices stop
    1-2 Å short of the inner face): an element that crosses the membrane (membrane.crosses) passes a face where its CA
    run does, or, stopping short of it, where the line through its run's first and last turn would, carried on along
    the element as drawn. One wholly on one side of the mid-plane gives no face."""
    from types import SimpleNamespace

    from foldmap.cli import _face_medians
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    z = np.arange(-13.5, 16.6, 1.5)  # 21 CAs, 1.5 Å short of the inner face, past the outer one
    near = SimpleNamespace(id="near", start=0, end=len(z) - 1)
    place = SimpleNamespace(n_port=(0.0, 0.0), c_port=(0.0, 10.0))  # N -> C up the page, 0.5 a CA
    bb = SimpleNamespace(ca=np.column_stack([np.zeros_like(z), np.zeros_like(z), z]))
    lay = SimpleNamespace(placed={"near": place})
    y_out, y_in = _face_medians(lay, [near], bb, mem)
    assert y_out == pytest.approx(10.0 * (15.0 + 13.5) / 30.0) and y_in == pytest.approx(10.0 * (-15.0 + 13.5) / 30.0)
    runs = {"short": np.arange(-10.0, 10.1, 1.0), "above": np.arange(2.0, 20.1, 1.0)}  # 20 Å across; wholly above
    elements = [near]
    for name, (k, run) in zip(runs, enumerate(runs.values())):
        start = len(bb.ca)
        bb.ca = np.vstack([bb.ca, np.column_stack([np.zeros_like(run), np.zeros_like(run), run])])
        elements.append(SimpleNamespace(id=name, start=start, end=start + len(run) - 1))
        lay.placed[name] = SimpleNamespace(n_port=(5.0 * k, 20.0), c_port=(5.0 * k, 30.0))
    y_out, y_in = _face_medians(lay, elements, bb, mem)  # short: faces at 32.5 and 17.5; above: none
    assert y_out == pytest.approx((9.5 + 32.5) / 2) and y_in == pytest.approx((-0.5 + 17.5) / 2)


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
    crossed past the middle). The view is laid out with every transmembrane element at its projected height (moved to
    meet a face when drawn wholly inside the slab; plan 26 Task 2: from the first layout on), and the band is the
    membrane's faces as projected: each enters through a face and reaches past the middle."""
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


def test_a_membrane_view_is_laid_out_once_at_its_projected_heights(monkeypatch):
    """Plan 26 Task 2 (C17's second layout made the first): a projected membrane view is laid out once, every
    transmembrane element at its projected height (`level`) against the faces as projected (`band`), and the band it
    draws is those faces. Before, a first layout drew the band at the medians of wherever its helices fell and a second
    at projected heights was kept only when A18 scored it better (7VVK R:332-345 was left short of the middle)."""
    import foldmap.cli as cli

    calls = []
    real = cli.build_layout

    def spy(*a, **k):
        calls.append(k)
        return real(*a, **k)

    monkeypatch.setattr(cli, "build_layout", spy)
    path = DATA / "1BL8.cif"
    lay, sses, bb = make_layout(path, membrane="auto")
    assert len(calls) == 1 and calls[0]["level"] and calls[0]["level"] == calls[0]["upright"]
    assert lay.membrane["y"] == calls[0]["band"]


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
    # ruling R3: 5LDW M:392-402, an 11-residue helix half the band long, reaches a face but not past the middle; it
    # is carried to Task 2 as an A18 fail, the only helix short here
    short = [x.split()[0] for x in view["tm_short"]]
    assert view["tm_helices"] >= 14 and short == (["M:392-402"] if entry == "5LDW" else []), view
    assert view["loops_in_slab"] == [], view


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


def test_a_membrane_two_fold_dimer_is_drawn_no_worse_than_unrolled(monkeypatch):
    """Plan 18 fix round 1: a C2 dimer seen from the side of its dyad is laid out again unrolled when its mirror image
    leaves a transmembrane helix off its band (A18) or draws the copies unalike (A15, judged on membrane views too);
    the unrolled layout is kept when it leaves fewer helices off the band, or as many and the copies more alike. 6JLN
    --membrane (its dyad along the normal): the drawn figure is no worse than the unrolled one on bands, then copies."""
    import foldmap.cli as cli
    from foldmap.cli import _unalike
    from foldmap.figcheck import membrane_view
    from helpers import r100_entry

    path = r100_entry("6JLN")

    def scored():
        lay, sses, bb = make_layout(path, membrane="auto")
        return len(membrane_view(lay, sses, bb, find_membrane(path, bb))["tm_short"]), _unalike(lay, bb)

    drawn = scored()
    monkeypatch.setattr(cli, "_dyad_frame", lambda *a: None)
    unrolled = scored()
    assert drawn <= unrolled, (drawn, unrolled)


# ------------------------------------------------------------------------------- plan 26: the whole structure
class _Chains:  # a backbone as membrane._orient reads it: CA atoms and residue labels (chain, number, name)
    def __init__(self, residues):
        from types import SimpleNamespace

        self.ca = np.array([xyz for _, _, xyz, _ in residues], float)
        self.labels = [SimpleNamespace(chain=c, seq=n, icode="", name=name) for c, n, _, name in residues]

    def __len__(self):
        return len(self.labels)


def _protein(chain, x, loop_kr):
    """One chain along z: an N-terminal loop of four residues 5 Å below the lower face (`loop_kr` of them Lys, the
    rest Gly), a helix across the 30 Å membrane (LEU), a 25-residue linker far above it, then a soluble domain of 20
    residues 3 Å above the upper face, 12 of them Lys. Returns its residues and its helix's (start, end) offsets."""
    out = [(chain, n + 1, (x, 0.0, -20.0), "LYS" if n < loop_kr else "GLY") for n in range(4)]
    helix = [(chain, len(out) + k + 1, (x, 0.0, z), "LEU") for k, z in enumerate(np.arange(-21.0, 21.1, 1.5))]
    span = (len(out), len(out) + len(helix) - 1)
    out += helix
    out += [(chain, len(out) + k + 1, (x + 5.0, 0.0, 60.0), "GLY") for k in range(25)]
    out += [(chain, len(out) + k + 1, (x + 5.0, float(k), 18.0), "LYS" if k < 12 else "ALA") for k in range(20)]
    return out, span


def _oriented(proteins):
    from foldmap import membrane

    residues, elements = [], []
    for chain, x, loop_kr in proteins:
        got, (a, b) = _protein(chain, x, loop_kr)
        elements.append(_Element("H", len(residues) + a, len(residues) + b, chain))
        residues += got
    bb = _Chains(residues)
    flat = membrane.Membrane(np.array([0.0, 0.0, 1.0]), 0.0, membrane.HALF, "estimate")
    return membrane._orient(flat, bb, elements), bb


def test_the_positive_inside_rule_counts_the_loops_not_a_soluble_domain_beyond():
    """Plan 26 (32KS, ComEC): the periplasmic OB and β-lactamase-like domains' own Lys/Arg within 15 Å of the face
    outvoted the loops and put them on the cytoplasmic side. The rule counts the loops of the elements that cross the
    membrane: a soluble domain on top with 12 Lys near the face, the helix's N-terminal loop under it with 4: the
    loops' side is the cytoplasm, the domain outside, and the sides are called."""
    m, bb = _oriented([("A", 0.0, 4)])
    assert m.normal @ np.array([0.0, 0.0, 1.0]) > 0 and m.sided
    domain = m.depth(bb.ca[-20:])
    assert (domain > 0).all()
    d0 = bb.ca[:, 2]  # every residue near the faces, as the rule counted before: 12 above, 4 below (the domain inside)
    basic = np.array([l.name == "LYS" for l in bb.labels])
    assert (basic & (d0 > 10) & (d0 < 30)).sum() == 12 and (basic & (d0 > -30) & (d0 < -10)).sum() == 4


def test_copies_of_one_protein_count_their_loops_once():
    """1BL8 (KcsA, four copies): its loops' Lys/Arg are 2 to 1 per copy, 8 to 4 summed, which called its sides the
    wrong way. Copies of one protein count once (their mean): a loop with 2 Lys under the membrane and none above
    calls no side (they must differ by 3), and neither do four copies of it, though summed (8 to 0) they would."""
    one, _ = _oriented([("A", 0.0, 2)])
    four, _ = _oriented([(c, 30.0 * k, 2) for k, c in enumerate("ABCD")])
    assert not one.sided and not four.sided
    assert one.normal @ np.array([0.0, 0.0, 1.0]) > 0 and np.allclose(four.normal, one.normal)


def test_the_membrane_is_found_on_every_chain_whatever_chains_draws(monkeypatch):
    """Plan 26 (32KS --chains A): the membrane was estimated on the chains --chains kept, so one protomer of a dimer
    got a plane of its own, ~22° off the dimer's. It is found on every protein chain read; the drawing narrows, the
    bilayer does not, and A18 (criteria.membrane_checks) judges the figure by the membrane it was drawn with."""
    from foldmap import criteria, membrane
    from foldmap.spec import draw_spec, figure_spec

    path = R100 / "6JO6-AB.cif.gz"
    whole = find_membrane(path, load_backbone(path, "auto"))
    d = draw_spec(figure_spec(path, membrane="auto", chains=["A"]))
    assert {l.chain for l in d.bb.labels} == {"A"}
    m = d.layout.bilayer
    assert m.normal.tolist() == whole.normal.tolist() and m.centre == whole.centre and m.sided == whole.sided

    def refound(*_):
        raise AssertionError("membrane_checks found the membrane again")

    monkeypatch.setattr(membrane, "find_membrane", refound)
    assert criteria.membrane_checks(d, path)["band"]


KS32 = R100_CACHE / "32KS.cif"  # the download cache (FOLDMAP_R100_CACHE overrides it), as the random-100 tests read it


@pytest.mark.skipif(not KS32.is_file(), reason="32KS is not in the download cache")
def test_comec_one_protomer_lies_in_the_dimers_membrane_with_its_periplasmic_domains_outside():
    """32KS (ComEC, Deselaers et al., NAR 2025, Fig. 1C) --chains A --membrane auto: the membrane is the dimer's, so
    TM1-TM3 (A:5-19, 29-43, 50-74) cross its mid-plane, inside the band (chain A alone put them below it: TM2 at 3 to
    23 Å, wholly on one side), and the OB fold (A:85-196) and the β-lactamase-like domain (A:531-800) lie outside, the
    sides called by the loops, with no --membrane flip."""
    whole = find_membrane(KS32, load_backbone(KS32, "auto"))
    lay, _, bb = make_layout(KS32, membrane="auto", chains=["A"])
    m = lay.bilayer
    assert m.normal.tolist() == whole.normal.tolist() and m.centre == whole.centre and m.sided
    d = m.depth(bb.ca)
    num = np.array([l.seq for l in bb.labels])
    for lo, hi in ((5, 19), (29, 43), (50, 74)):
        seg = d[(num >= lo) & (num <= hi)]
        assert seg.min() < 0 < seg.max() and seg.min() > -m.half - 5 and seg.max() < m.half + 5, (lo, hi)
    for lo, hi in ((85, 196), (531, 800)):
        assert d[(num >= lo) & (num <= hi)].mean() > m.half, (lo, hi)


COMEC_TMS = [(5, 19), (29, 43), (50, 74), (261, 277), (281, 298), (303, 320), (328, 339), (349, 363), (377, 399)]
COMEC_TMS += [(475, 493), (508, 526)]  # chain A's eleven upright helices, the paper's TM1-TM11 by place and tilt


@pytest.mark.skipif(not KS32.is_file(), reason="32KS is not in the download cache")
@pytest.mark.parametrize("chains", [None, ["A"]], ids=["dimer", "A"])
def test_comec_transmembrane_helices_count_as_crossing_and_stand_in_the_band(chains):
    """Plan 26 (ruling R1): the band ran through the OB barrels and A18 judged 2 of the paper's 11 TM helices per
    chain, the others stopping short of a face of the 30 Å slab (FACE_REACH). Counted by membrane.crosses, 10 of
    chain A's 11 cross (A:328-339 spans 14 Å), in the dimer and with --chains A alike, 10 helices in each chain,
    and A18 finds every one drawn across its band."""
    from foldmap import membrane
    from foldmap.figcheck import membrane_view
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(KS32, membrane="auto", chains=chains)
    m = lay.bilayer
    d = m.depth(bb.ca)
    view = membrane_view(lay, sses, bb, m, route_loops(lay, sses, bb))
    drawn = sorted({l.chain for l in bb.labels})
    for chain in drawn:  # B models no A:508-526 (TM11); its TM7, B:324-341, is 18 residues and crosses
        tm = [
            s.id
            for s in sses
            if s.chain == chain and s.kind == "H" and membrane.crosses(d[s.start : s.end + 1], m.half)
        ]
        assert len(tm) == 10, (chain, tm)
    on, num = np.array([l.chain == "A" for l in bb.labels]), np.array([l.seq for l in bb.labels])
    counted = [membrane.crosses(d[on & (num >= lo) & (num <= hi)], m.half) for lo, hi in COMEC_TMS]
    assert sum(counted) == 10 and not counted[6], counted
    assert view["tm_helices"] == 10 * len(drawn) and view["tm_short"] == [] and view["loops_in_slab"] == [], view


# ------------------------------------------------------- plan 26 fix round 1: the band test and _to_face, in CI
def test_a_tm_helix_crosses_its_band_only_past_the_middle_and_from_a_face():
    """C27 as ruling R3 reads it: an end on the middle is not past it (strict); a face is reached within
    layout.BAND_TOUCH (1e-9, float noise: 6DQN's ports ~4e-15 short), no further."""
    from foldmap.layout import BAND_TOUCH, crosses_band

    band = (0.0, 10.0)
    assert BAND_TOUCH == 1e-9
    assert not crosses_band((5.0, 10.0), band)  # from the upper face to the middle exactly: not past it
    assert crosses_band((5.0 - 1e-6, 10.0), band)
    assert crosses_band((2.0, 10.0 - 1e-12), band) and not crosses_band((2.0, 10.0 - 1e-6), band)
    assert crosses_band((1e-12, 6.0), band) and not crosses_band((1e-6, 6.0), band)
    assert crosses_band((6.0, -1.0), band)  # either way up


def _upright_item(ends):
    """An item at height 0 holding upright transmembrane helices drawn from page y lo to hi, as layout._to_face reads
    one: each member's id, angle, length and ghost flag, and its offset from the item's centre."""
    from types import SimpleNamespace

    members = [
        SimpleNamespace(sse=SimpleNamespace(id=f"t{k}"), angle=np.pi / 2, length=hi - lo, ghost=False)
        for k, (lo, hi) in enumerate(ends)
    ]
    return SimpleNamespace(members=members, offsets=[(0.0, (lo + hi) / 2) for lo, hi in ends])


def test_to_face_leaves_every_tm_helix_of_an_item_crossing_the_band():
    """6ZKT N:26-44 (plan 26): an item's first move brings the helix drawn wholly inside the band (2..8) down to the
    lower face, and carries its partner (-3..7) to end on the middle, short of it. The least move that takes ends
    layout._PAST beyond a face or the middle and leaves both crossing is taken instead: up 2.1, the inner helix past
    the upper face."""
    from foldmap.layout import _PAST, _to_face, crosses_band

    band, ends = (0.0, 10.0), [(2.0, 8.0), (-3.0, 7.0)]
    assert not crosses_band((-3.0 - 2.0, 7.0 - 2.0), band)  # the first move (down 2) leaves the partner short
    it = _upright_item(ends)
    shift = _to_face(it, {"t0", "t1"}, band, 0.0)
    assert shift == pytest.approx(10.0 + _PAST - 8.0)
    assert all(crosses_band((lo + shift, hi + shift), band) for lo, hi in ends)


def test_to_face_cannot_make_a_helix_half_the_band_long_cross_it():
    """5LDW M:392-402 (ruling R3): an 11-residue helix drawn upright is exactly half the band long, so it can reach a
    face or pass the middle, never both; _to_face leaves it where the first move put it and A18 reports it short."""
    from foldmap.layout import _to_face, crosses_band

    band = (0.0, 10.0)
    it = _upright_item([(6.0, 11.0)])
    shift = _to_face(it, {"t0"}, band, 0.0)
    assert shift == 0.0 and not crosses_band((6.0, 11.0), band)
    assert not any(crosses_band((lo, lo + 5.0), band) for lo in np.round(np.arange(-6.0, 11.0, 0.01), 2))


def test_upright_helices_stand_up_the_page_whatever_their_tilt():
    """build_layout `upright` (a membrane view's transmembrane helices, cli): a helix seen nearly across the page
    (7VVK's TM1, tilted 52°, was snapped flat across the outer face) is drawn up or down the page; without it, the
    default snapped angle lays it across."""
    from foldmap.frame import Frame
    from foldmap.layout import build_layout
    from helpers import fake_sse

    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    s = fake_sse("H", 10, 30, centroid=(0, 0, 0), axis=(1.0, 0.6, 0.0))
    across = build_layout([s], [], f, "projected").placed[s.id]
    up = build_layout([s], [], f, "projected", upright={s.id}).placed[s.id]
    assert abs(np.sin(across.angle)) < 1e-9 and abs(np.sin(up.angle)) == pytest.approx(1.0)


# --------------------------------------------- plan 26 fix round 2: soluble elements in the slab, loops stop at gaps
def test_a18_flags_a_soluble_element_drawn_in_the_slab_not_an_interfacial_or_crossing_one():
    """Ruling R4: A18's sub-check "soluble element drawn in the slab" (figcheck.membrane_view soluble_in_slab). On a
    30 Å slab drawn as the band 0..10: an element whose CA centre lies 10 Å beyond the outer face, drawn with its centre
    in the band, is flagged; the same element drawn above the band is not, nor is one lying 3 Å beyond the face
    (interfacial: within SOLUBLE_BEYOND), nor a helix crossing the membrane (drawn across the band)."""
    from types import SimpleNamespace

    from foldmap.figcheck import SOLUBLE_BEYOND, membrane_view
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    runs = {  # id: (CA depths, kind, drawn centre y)
        "far": (np.full(6, 25.0), "H", 5.0),
        "far-above": (np.full(6, 25.0), "E", 14.0),
        "face": (np.full(6, 18.0), "H", 5.0),
        "across": (np.arange(-21.0, 21.1, 1.5), "H", 5.0),
    }
    ca, sses, placed = [], [], {}
    for name, (z, kind, cy) in runs.items():
        sses.append(SimpleNamespace(id=name, kind=kind, start=len(ca), end=len(ca) + len(z) - 1))
        ca += [(0.0, 0.0, float(v)) for v in z]
        placed[name] = SimpleNamespace(cy=cy, n_port=(0.0, cy - 6.0), c_port=(0.0, cy + 6.0))
    lay = SimpleNamespace(membrane={"y": (0.0, 10.0), "bands": [(0.0, 10.0)]}, placed=placed)
    bb = SimpleNamespace(ca=np.array(ca))
    view = membrane_view(lay, sses, bb, mem)
    assert 18.0 - mem.half < SOLUBLE_BEYOND < 25.0 - mem.half
    assert view["soluble_in_slab"] == ["far +25"]
    assert view["tm_helices"] == 1 and view["tm_short"] == []


def test_the_positive_inside_loops_stop_at_a_gap_in_the_numbering():
    """Review M6 (plan 26): membrane._loops walked 20 residues by backbone index, so past an unmodelled stretch it
    reached into whatever was modelled next. A helix across the membrane, three residues after it numbered on, then a
    gap of 30 in the numbering and five Lys 5 Å beyond the face: the three are its loop, the Lys are not."""
    from foldmap import membrane

    helix = [("A", k + 1, (0.0, 0.0, z), "LEU") for k, z in enumerate(np.arange(-21.0, 21.1, 1.5))]
    after = [("A", len(helix) + 1 + k, (0.0, 1.0, 20.0), "GLY") for k in range(3)]
    beyond = [("A", len(helix) + 34 + k, (0.0, 2.0, 20.0), "LYS") for k in range(5)]
    bb = _Chains(helix + after + beyond)
    m = membrane.Membrane(np.array([0.0, 0.0, 1.0]), 0.0, membrane.HALF, "estimate")
    loops = membrane._loops(m, bb, [_Element("H", 0, len(helix) - 1, "A")])
    n = len(helix)
    assert loops[n : n + 3].all() and not loops[n + 3 :].any()


# ---------------------------------------------------- plan 26 Task 2: upright and flat helices, the band kept clear
def _page():
    """A frame whose page-up is the membrane normal (z), page x is x, and y points at the viewer."""
    from foldmap.frame import Frame

    return Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 0, 1.0]), np.array([0, -1.0, 0]))


def test_a_helix_lies_flat_by_its_tilt_to_the_normal_not_as_it_is_seen():
    """membrane.lies_flat (FLAT_TILT, 45°): a helix tilted 44° from the normal stands, one tilted 46° lies. 32KS's
    amphipathic NS1 (A:200-219, 71° from the normal) points into the page, so the snapped angle stood it upright;
    build_layout `flat` lays it along the page, as its tilt says."""
    from foldmap.layout import build_layout
    from foldmap.membrane import FLAT_TILT, lies_flat
    from helpers import fake_sse

    normal = np.array([0.0, 0.0, 1.0])

    def tilted(deg):  # an axis tilted `deg` from the normal, its sideways part into the page
        return np.array([0.0, np.sin(np.radians(deg)), np.cos(np.radians(deg))])

    assert FLAT_TILT == 45.0
    assert not lies_flat(tilted(44.0), normal) and lies_flat(tilted(46.0), normal) and lies_flat(-tilted(71.0), normal)
    s = fake_sse("H", 10, 30, centroid=(0, 0, 0), axis=tilted(71.0))
    seen = build_layout([s], [], _page(), "projected").placed[s.id]
    flat = build_layout([s], [], _page(), "projected", flat={s.id}).placed[s.id]
    assert abs(np.sin(seen.angle)) == pytest.approx(1.0) and abs(np.sin(flat.angle)) < 1e-9


def test_a_flat_helix_lies_on_its_face_and_a_slab_element_keeps_its_height():
    """build_layout `on_face` with `band` (cli._membrane_places): a flat helix at a face, projected inside the band,
    is drawn centred on that face (+1 the upper, -1 the lower); an element in the slab crossing no membrane (0) keeps
    its projected height while a free element seen at the same place is moved off it."""
    from foldmap.layout import build_layout
    from helpers import fake_sse

    band = (-3.45, 3.45)
    upper = fake_sse("H", 0, 15, centroid=(0, 0, 10.0), axis=(1, 0, 0))
    lower = fake_sse("H", 20, 35, centroid=(30, 0, -12.0), axis=(1, 0, 0))
    slab = fake_sse("H", 40, 50, centroid=(60, 0, -4.0), axis=(0, 0, 1))
    free = fake_sse("H", 60, 70, centroid=(60, 5, -4.0), axis=(0, 0, 1))
    els = [upper, lower, slab, free]
    lay = build_layout(
        els,
        [],
        _page(),
        "projected",
        flat={upper.id, lower.id},
        band=band,
        on_face={upper.id: 1, lower.id: -1, slab.id: 0},
    )
    p = lay.placed
    assert p[upper.id].cy == pytest.approx(band[1]) and p[lower.id].cy == pytest.approx(band[0])
    assert abs(np.sin(p[upper.id].angle)) < 1e-9 and abs(np.sin(p[lower.id].angle)) < 1e-9
    assert p[slab.id].cy == pytest.approx(-4.0 * 0.23)
    a, b = p[slab.id].rect, p[free.id].rect
    assert a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1]  # the free one moved off the held one


def test_a_soluble_element_is_kept_off_the_band_on_its_own_side():
    """build_layout `sides` (ruling R4, plan 26 Task 2): two soluble helices projected into the band (7PIU's Gβ was
    pulled there by its loops; 32KS chain B's OB sheet by its neighbours) are each kept off it on their own side, MARGIN
    clear of the face (layout._membrane_walls); without `sides` they stay in the band."""
    from foldmap.layout import MARGIN, build_layout
    from helpers import fake_sse

    band = (-3.45, 3.45)
    up = fake_sse("H", 0, 12, centroid=(0, 0, 2.0), axis=(1, 0, 0))
    down = fake_sse("H", 20, 32, centroid=(20, 0, -2.0), axis=(1, 0, 0))
    plain = build_layout([up, down], [], _page(), "projected", band=band).placed
    assert all(band[0] < plain[s.id].cy < band[1] for s in (up, down))
    kept = build_layout([up, down], [], _page(), "projected", band=band, sides={up.id: 1, down.id: -1}).placed
    assert kept[up.id].rect[1] >= band[1] + MARGIN - 1e-6
    assert kept[down.id].rect[3] <= band[0] - MARGIN + 1e-6


def test_an_untied_copy_item_moves_with_its_protomer():
    """layout._untied_spread (32KS chain B's OB sheet, 22 Å off its depth): in an unrolled two-fold assembly the ties
    set protomer 1's copies a spread right of their mates, wider than the frame's own period; an item of protomer 1 no
    tie places (its own sheet where DSSP split its mate's) is moved by the same amount less the frame's offset, so it
    starts beside its own protomer's elements; protomer 0's items and the tied ones stay where they are."""
    from foldmap.frame import CylinderFrame
    from foldmap.layout import SCALE, Placed, _Item, _untied_spread
    from helpers import fake_sse

    frame = CylinderFrame(np.zeros(3), np.array([0, 0, 1.0]), np.array([1.0, 0, 0]), 2, 10.0, {"A": 0, "B": 1})
    items = []
    for k, chain in enumerate("ABBA"):
        s = fake_sse("H", 10 * k, 10 * k + 8, chain=chain)
        items.append(_Item([Placed(s, 0, 0, 3.0, 0.0, 1.0)], [(0.0, 0.0)], float(k), 1.0, (float(k), 1.0), k, chain))
    spread = 50.0
    _untied_spread(items, {1: (0, 1)}, frame, spread)
    shift = spread - 10.0 * np.pi * SCALE
    assert [it.cx for it in items] == pytest.approx([0.0, 1.0, 2.0 + shift, 3.0])
    assert items[2].home == pytest.approx((2.0 + shift, 1.0)) and all(it.cy == 1.0 for it in items)
