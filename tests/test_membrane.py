"""Membrane bands: the lipid bilayer drawn behind transmembrane proteins, with its inside and outside."""

import functools
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


def test_a_membrane_view_is_laid_out_once_at_its_projected_heights(monkeypatch):
    """A guard against C17's two-layout path coming back (review-2 M2; the behaviour, every TM helix across the band at
    the faces as projected, is test_every_transmembrane_helix_crosses_the_band_in_the_membrane_view's): a projected
    membrane view is laid out once, every transmembrane element at its projected height (`level`) against the faces
    as projected (`band`), and the band it draws is those faces. Before, a first layout drew the band at the medians
    of wherever its helices fell and a second at projected heights was kept only when A18 scored it better (7VVK
    R:332-345 was left short of the middle). 1BL8 lets no copy go, so one call it is."""
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
    sideways starts as before. Plan 27 Task 2: the router is held here to its own measure (helpers.sideways_in_slab,
    A18's loop test before it was made strict); test_a18_holds_loops_beyond_a_face_out_of_the_slab holds it to A18."""
    from foldmap.figcheck import membrane_view
    from foldmap.route import route_loops
    from helpers import r100_entry, sideways_in_slab

    path = r100_entry(entry)
    lay, sses, bb = make_layout(path, membrane="auto")
    loops = route_loops(lay, sses, bb)
    view = membrane_view(lay, sses, bb, find_membrane(path, bb), loops)
    # ruling R6: 5LDW M:392-402, an 11-residue helix half the band long, is drawn a hair longer and crosses
    assert view["tm_helices"] >= 14 and view["tm_short"] == [], view
    assert sideways_in_slab(loops, lay.membrane["bands"], lay.placed) == []


# plan 27 rule P7: A18 made strict before the layout is (Tasks 4 and 5 take these marks off)
PENDING = pytest.mark.xfail(  # only a failed assertion is the expected failure, never a crash (review-2 M4)
    strict=True, raises=AssertionError, reason="plan 27 T4/T5: layout not yet by construction"
)


@pytest.mark.parametrize("entry", [pytest.param(e, marks=PENDING) for e in ("5LDW", "6JO6", "6JLN")])
def test_a18_holds_loops_beyond_a_face_out_of_the_slab(entry):
    """Plan 27 Task 2 (A18 strict, rulings P4, P14): random-100 C2b-5's membrane views (5LDW, 6JO6, 6JLN), no loop drawn
    nearer the mid-plane than its residues (its face, when all lie beyond it) but on its band element's port stub
    (loop_in_band)."""
    from foldmap.figcheck import membrane_view
    from foldmap.route import route_loops
    from helpers import r100_entry

    path = r100_entry(entry)
    lay, sses, bb = make_layout(path, membrane="auto")
    view = membrane_view(lay, sses, bb, lay.bilayer, route_loops(lay, sses, bb))
    assert view["loop_in_band"] == [], view["loop_in_band"]


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
    assert view["tm_helices"] == 10 * len(drawn) and view["tm_short"] == [], view  # A18 strict: test_a18_on_comec


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


# ------------------------------------------------ plan 27 Task 2: A18 strict, every body and loop segment, no tolerance
UPRIGHT, FLAT = (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)


def _y(depth: float) -> float:
    """Page y of a depth (Å) on _a18's page: the 30 Å slab is the band 0..10, the outside up, 3 Å a page unit."""
    return 5.0 + depth / 3.0


def _a18(runs, loops=None, defining=None):
    """figcheck.membrane_view on a 30 Å slab (normal z, mid-plane at 0) drawn as the band 0..10, the outside up (_y).
    runs: id -> (kind, CA depths, axis, (N port, C port)), in sequence: a port is a page y (at x = 0) or an (x, y)
    point; each element is drawn 1 wide about the line between its ports, its rect that body's bounding box; kind
    None: residues no element holds (a loop's). loops: (a, b, points)."""
    from types import SimpleNamespace

    from foldmap.figcheck import membrane_view
    from foldmap.membrane import Membrane
    from foldmap.route import Loop

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate", defining=defining)
    ca, sses, placed = [], [], {}
    for name, (kind, z, axis, ports) in runs.items():
        if kind is not None:
            s = SimpleNamespace(id=name, kind=kind, axis=np.array(axis), start=len(ca), end=len(ca) + len(z) - 1)
            sses.append(s)
            n, c = (np.array(q, float) if np.ndim(q) else np.array([0.0, q]) for q in ports)
            d = c - n
            d = d / np.hypot(*d) if np.hypot(*d) > 0 else np.array([0.0, 1.0])
            side = np.array([-d[1], d[0]]) / 2
            corners = np.array([n + side, n - side, c + side, c - side])
            rect = (*corners.min(axis=0), *corners.max(axis=0))
            placed[name] = SimpleNamespace(sse=s, rect=rect, n_port=tuple(n), c_port=tuple(c), cy=(n[1] + c[1]) / 2)
        ca += [(0.0, 0.0, float(v)) for v in z]
    band = (0.0, 10.0)
    lay = SimpleNamespace(membrane={"y": band, "bands": [band], "inside": "below"}, placed=placed, ghosts=[])
    drawn = None if loops is None else [Loop(a, b, list(pts), False, False) for a, b, pts in loops]
    return membrane_view(lay, sses, SimpleNamespace(ca=np.array(ca)), mem, drawn), sses, mem, np.array(ca)[:, 2]


TM = ("H", np.linspace(-18.0, 18.0, 25), UPRIGHT, (_y(-18.0), _y(18.0)))


def test_the_band_holds_crossing_helices_helices_in_the_slab_and_crossing_strands():
    """Plan 27 rulings P5, P10, P11 (membrane.in_band, the one rule A18 and the layout take the band by): the band holds
    a helix crossing the membrane; a helix or 3-10 helix whose CA median lies in the slab, at any tilt (P10: flat, or
    standing short of crossing, as ComEC's TM7); and a strand crossing the membrane, of a barrel or in a barrel's
    lumen (P11: 9FWZ's FimA, of no barrel defining the membrane). Not a helix whose CA median lies beyond a face, flat
    or standing, nor an interfacial strand, even one Membrane.defining names."""
    from foldmap.membrane import band_set, own_side, slab_helix

    runs = {
        "tm": TM,
        "flat-in": ("H", np.array([9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0]), FLAT, (0, 0)),  # median 12
        "eta-in": ("G", np.full(3, -8.0), FLAT, (0, 0)),
        "standing": ("H", np.linspace(2.0, 12.0, 8), UPRIGHT, (0, 0)),
        "flat-out": ("H", np.array([13.0, 14.0, 15.5, 16.0, 17.0]), FLAT, (0, 0)),  # median 15.5: beyond the face
        "standing-out": ("H", np.linspace(14.0, 24.0, 8), UPRIGHT, (0, 0)),
        "crossing-strand": ("E", np.linspace(-14.0, 14.0, 9), UPRIGHT, (0, 0)),
        "interfacial": ("E", np.linspace(10.4, 17.0, 3), UPRIGHT, (0, 0)),
    }
    _, sses, mem, depth = _a18(runs)
    expected = {"tm", "flat-in", "eta-in", "standing", "crossing-strand"}
    assert band_set(sses, depth, mem) == expected
    assert {s.id for s in sses if slab_helix(s, depth, mem)} == {"flat-in", "eta-in", "standing"}
    by = {s.id: s for s in sses}
    for named in ("tm", "interfacial"):  # a barrel elsewhere (the strand crossing in its lumen), or named interfacial
        barrel = frozenset(range(by[named].start, by[named].end + 1))
        _, sses, mem, depth = _a18(runs, defining=barrel)
        assert band_set(sses, depth, mem) == expected, named
    assert own_side(np.array([-1.0, 1.0])) == 1 and own_side(np.array([-2.0, 1.0])) == -1


def test_a18_holds_a_helix_in_the_slab_by_its_axis_at_its_own_depth_range():
    """Rulings P10, P13 (slab_helix): a helix of the band crossing no membrane, its CA median in the slab, is judged by
    its axis (its ports), not its body (Å off, the furthest). Standing, the axis lies inside the band and within its
    CA atoms' depths: not stretched across the band (ComEC's TM7, A:328-339, -8.8 to +5.1 Å, drawn face to face), not
    1 Å deeper than its deepest CA, not past a face though its CA atoms reach there (5LDW L:212-222, 4.0 to 18.3 Å);
    drawn to its depth span it passes, tilted on the page too (its body's corners reach further). Flat, the axis lies
    within its CA depths extended to its nearer face: on the face line (Fig. 1C's NS7, A:436-461, its CA atoms past
    it or all inside), or inside no deeper than its deepest CA; not on the face row outside the band, nor deeper. None
    is judged as an element off the band."""
    tm7 = np.linspace(-8.8, 5.1, 12)
    flat = np.array([10.5, 11.0, 12.0, 12.5, 13.0, 14.9])
    ns7 = np.array([12.5, 13.0, 13.9, 14.5, 16.0, 18.7])  # median 14.2, its outer CA atoms beyond the face

    def across(depth):  # a flat helix drawn along the page, its axis at this depth
        return ((-2.0, _y(depth)), (2.0, _y(depth)))

    runs = {
        "tm7": ("H", tm7, UPRIGHT, (_y(-8.8), _y(5.1))),
        "tilted": ("H", tm7, UPRIGHT, ((-1.0, _y(-8.8)), (1.0, _y(5.1)))),
        "stretched": ("H", tm7, UPRIGHT, (_y(-15.0), _y(15.0))),
        "deeper": ("H", tm7, UPRIGHT, (_y(-9.8), _y(4.1))),
        "past-face": ("H", np.linspace(4.0, 18.0, 8), UPRIGHT, (_y(4.0), _y(18.0))),  # at its depths, past the face
        "on-face": ("H", ns7, FLAT, across(15.0)),
        "on-face-below": ("H", np.array([11.0, 12.0, 12.5, 13.0, 13.5, 14.0]), FLAT, across(15.0)),  # CA all inside
        "inside": ("H", flat, FLAT, across(12.7)),
        "face-row": ("H", np.full(6, 12.0), FLAT, across(16.6)),
        "flat-deep": ("H", flat, FLAT, across(9.5)),
    }
    view, *_ = _a18(runs)
    assert view["slab_helix"] == ["stretched 9.9", "deeper 1.0", "past-face 3.0", "face-row 1.6", "flat-deep 1.0"]
    assert view["band_member"] == [] and view["outside_in_band"] == [] and view["tm_helices"] == 0, view


def test_a18_flags_every_element_of_no_band_whose_drawn_body_overlaps_it():
    """Rulings P3, P4 (band_member, strand_in_band): an element not of the band whose drawn body (its whole rect)
    overlaps the band, by any amount past float snapping, is flagged, how far in from its own face: the interfacial
    hairpin strand drawn at its depth (32KS B:402-404, 10.4 to 17.0 Å), a soluble helix dipping 1 Å in. Among them the
    β-strands are flagged again. Not one drawn from the face out, nor a crossing helix or a helix flat in the slab drawn
    in it."""
    runs = {
        "tm": TM,
        "hairpin": ("E", np.linspace(10.4, 17.0, 3), UPRIGHT, (_y(10.5), _y(17.1))),
        "on-face": ("E", np.linspace(10.4, 17.0, 3), UPRIGHT, (10.0, 12.2)),
        "hair": ("E", np.linspace(16.0, 20.0, 3), UPRIGHT, (10.0 - 1e-3, 11.0)),  # 0.003 Å in: no tolerance
        "snapped": ("E", np.linspace(16.0, 20.0, 3), UPRIGHT, (10.0 - 1e-7, 11.0)),  # float snapping: on the face
        "dips": ("H", np.linspace(16.0, 24.0, 6), UPRIGHT, (_y(14.0), _y(22.0))),
        "flat-in": ("H", np.full(6, 11.0), FLAT, (_y(11.0), _y(11.0))),
    }
    view, *_ = _a18(runs)
    assert view["band_member"] == ["hairpin 4.5", "hair <0.1", "dips 1.0"], view
    assert view["strand_in_band"] == ["hairpin 4.5", "hair <0.1"], view
    assert view["tm_helices"] == 1 and view["tm_short"] == []


def test_a18_flags_an_element_drawn_nearer_the_mid_plane_than_its_nearest_residue_or_on_the_wrong_side():
    """Rulings P4 (outside_in_band, wrong_side): an element not of the band is drawn no nearer the mid-plane, on its own
    side, than its nearest CA atom (Å nearer), and its drawn centre lies on its own side (membrane.own_side; Å past the
    mid-plane), outside the band too: a soluble helix 25-30 Å out drawn from 20 Å; one 31-34 Å in drawn 50-55 Å out
    (6ZKT α4, drawn below the band, the other way: both); a short strand straddling the mid-plane (-6 to +8 Å, its own
    side the outside) drawn from -5 to +1 Å, no nearer than its residues but its centre inside. Not one drawn from its
    nearest residue out, nor further out, nor one inside drawn inside."""
    runs = {
        "near": ("H", np.linspace(25.0, 30.0, 4), UPRIGHT, (_y(20.0), _y(27.0))),
        "at-bound": ("H", np.linspace(25.0, 30.0, 4), UPRIGHT, (_y(25.0), _y(32.0))),
        "far-out": ("H", np.linspace(25.0, 30.0, 4), UPRIGHT, (_y(40.0), _y(47.0))),
        "flipped": ("H", np.linspace(-34.0, -31.0, 4), UPRIGHT, (_y(50.0), _y(55.0))),
        "inner": ("H", np.linspace(-34.0, -31.0, 4), UPRIGHT, (_y(-40.0), _y(-35.0))),
        "straddle": ("E", np.linspace(-6.0, 8.0, 3), UPRIGHT, (_y(-5.0), _y(1.0))),
    }
    view, *_ = _a18(runs)
    assert view["outside_in_band"] == ["near 5.0", "flipped 86.0"], view
    assert view["wrong_side"] == ["flipped 52.5", "straddle 2.0"], view
    assert view["band_member"] == ["straddle 20.0"] and view["strand_in_band"] == ["straddle 20.0"], view


def test_a18_flags_a_loop_drawn_nearer_the_mid_plane_than_its_residues_but_on_its_ports_stub():
    """Rulings P4, P14 (loop_in_band): no loop is drawn nearer the mid-plane, on its own side, than its bound (Å
    nearer): its nearest residue, or its face where every residue lies beyond it. Between two elements outside the band
    a loop of residues beyond the face stays out of the slab (32KS A:436-461>A:467-469, dipping 3-5 Å in); out of a band
    element's port it runs only straight along the normal to that face, crossing the face line at the port (a TM
    helix ending short of the face): a turn inside the slab, or a first step away from the face, is in it. A loop one
    of whose residues lies inside the face (10 Å: a TM exit DSSP ends short) may dip to that residue, no further: drawn
    sideways along the mid-plane it is 10 Å too deep. Loops not given, the part is not measured (None)."""
    out = np.full(4, 20.0)
    runs = {
        "tm1": ("H", np.linspace(-18.0, 13.5, 22), UPRIGHT, (_y(-18.0), 9.0)),  # its C port 1.5 Å inside the face
        "l1": (None, out, None, (0, 0)),
        "s1": ("E", np.full(3, 19.0), UPRIGHT, (11.0, 12.0)),
        "l2": (None, out, None, (0, 0)),
        "s2": ("E", np.full(3, 19.0), UPRIGHT, (12.0, 11.0)),
        "l3": (None, np.array([20.0, 10.0, 20.0]), None, (0, 0)),  # one residue 5 Å inside the face
        "s3": ("E", np.full(3, 19.0), UPRIGHT, (11.0, 12.0)),
    }
    up = [(0.0, 9.0), (0.0, 10.6), (3.0, 10.6), (3.0, 11.0)]  # the stub to the face, on along it
    turn = [(0.0, 9.0), (0.0, 9.5), (2.0, 9.5), (2.0, 10.6), (3.0, 10.6), (3.0, 11.0)]  # turns in the slab
    away = [(0.0, 9.0), (0.0, 8.0), (0.0, 10.6), (3.0, 10.6), (3.0, 11.0)]  # first down, then up its stub
    dip = [(0.0, 12.0), (0.0, 9.7), (4.0, 9.7), (4.0, 12.0)]  # between two outside the band
    over = [(0.0, 12.0), (0.0, 10.0), (4.0, 10.0), (4.0, 12.0)]  # along the face line: not in the slab
    to_residue = [(0.0, 11.0), (0.0, _y(10.0)), (4.0, _y(10.0)), (4.0, 11.0)]  # down to its 10 Å residue
    mid_slab = [(0.0, 11.0), (0.0, 5.0), (4.0, 5.0), (4.0, 11.0)]  # sideways along the mid-plane
    loops = [("tm1", "s1", up), ("s1", "s2", over), ("s2", "s3", to_residue)]
    view, *_ = _a18(runs, loops)
    assert view["loop_in_band"] == [], view
    bad = [("tm1", "s1", turn), ("tm1", "s1", away), ("s1", "s2", dip), ("s2", "s3", mid_slab)]
    view, *_ = _a18(runs, bad)
    assert view["loop_in_band"] == ["tm1>s1 1.5", "tm1>s1 6.0", "s1>s2 0.9", "s2>s3 10.0"], view
    assert _a18(runs)[0]["loop_in_band"] is None


def test_a18_judges_its_parts_each_by_name_and_none_unmeasured_as_a_pass():
    """Plan 27 Task 2: A18's detail names each part's findings (criteria.MEMBRANE_WORDS), at most MEMBRANE_SHOWN per
    part and the rest counted (review-2 M3), and fails on any; a part a result does not hold (measured before plan 27)
    is not measured, never a pass. A finding under 0.05 Å reads "<0.1", never "0.0"."""
    from foldmap import criteria
    from foldmap.figcheck import _angstrom

    clean = {"band": True, "tm_helices": 7, "tm_short": [], **{k: [] for k in criteria.MEMBRANE_PARTS}}
    judged = criteria.figure_criteria({"membrane": {"result": {"membrane": clean}}})["A18"]
    assert judged["pass"] is True
    for part in criteria.MEMBRANE_PARTS:
        one = criteria.figure_criteria({"membrane": {"result": {"membrane": {**clean, part: ["A:1-9 2.5"]}}}})["A18"]
        assert one["pass"] is False and one["detail"] == f"{criteria.MEMBRANE_WORDS[part]} A:1-9 2.5 Å", part
        old = {k: v for k, v in clean.items() if k != part}
        assert criteria.figure_criteria({"membrane": {"result": {"membrane": old}}})["A18"]["pass"] is False
    many = [f"A:{k}-{k + 4} 1.0" for k in range(criteria.MEMBRANE_SHOWN + 3)]
    capped = criteria.figure_criteria({"membrane": {"result": {"membrane": {**clean, "wrong_side": many}}}})["A18"]
    shown = capped["detail"].split("; ")
    assert len(shown) == criteria.MEMBRANE_SHOWN + 1 and shown[-1] == "(+3 more wrong_side)", shown
    assert (_angstrom(0.049), _angstrom(0.05), _angstrom(2.25)) == ("<0.1", "0.1", "2.2")


# ---------------------------------------------------------------- plan 26 fix round 2: loops stop at gaps
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
    """build_layout `on_face` with `band` (cli._membrane_places, ruling R8): a flat helix at a face, projected inside
    the band, lies on that face's first row, just outside the band (+1 the upper face, -1 the lower), within its chain's
    transmembrane span; one of a chain with no transmembrane element is released to the ordinary layout and kept off
    the band like a soluble element; an element in the slab crossing no membrane (0) keeps its projected height while a
    free element seen at the same place is moved off it."""
    from foldmap.layout import MARGIN, SCALE, build_layout
    from helpers import fake_sse

    band = (-3.45, 3.45)
    tm = fake_sse("H", 0, 25, centroid=(0, 0, 0.0), axis=(0, 0, 1))
    upper = fake_sse("H", 30, 41, centroid=(5, 0, 10.0), axis=(1, 0, 0))
    lower = fake_sse("H", 45, 56, centroid=(-5, 0, -12.0), axis=(1, 0, 0))
    slab = fake_sse("H", 60, 70, centroid=(60, 0, -4.0), axis=(0, 0, 1))
    free = fake_sse("H", 80, 90, centroid=(60, 5, -4.0), axis=(0, 0, 1))
    other = fake_sse("H", 0, 11, centroid=(-40, 0, 10.0), axis=(1, 0, 0), chain="B")
    els = [tm, upper, lower, slab, free, other]
    flat = {upper.id, lower.id, other.id}
    on = {upper.id: 1, lower.id: -1, slab.id: 0, other.id: 1}
    lay = build_layout(els, [], _page(), "projected", level={tm.id}, upright={tm.id}, flat=flat, band=band, on_face=on)
    p = lay.placed
    assert band[1] + MARGIN - 1e-6 <= p[upper.id].rect[1] < band[1] + 2 * MARGIN
    assert band[0] - 2 * MARGIN < p[lower.id].rect[3] <= band[0] - MARGIN + 1e-6
    assert abs(np.sin(p[upper.id].angle)) < 1e-9 and abs(np.sin(p[lower.id].angle)) < 1e-9
    assert p[other.id].rect[1] >= band[1] + MARGIN - 1e-6  # released, off the band on its side
    assert p[slab.id].cy == pytest.approx(-4.0 * SCALE)
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


def test_membrane_places_reads_each_element_by_its_depth_and_tilt():
    """cli._membrane_places (review-2 I4) on a 30 Å slab, normal up: a lying helix 4.9 Å beyond the outer face lies on
    it, one 5.1 Å beyond is soluble; a standing helix and a strand in the slab keep their height (0), a lying helix in
    the slab lies on the face its centre is nearer; a 3-10 helix in the slab is placed as any other; an element crossing
    the membrane is none of these. The face sign follows the page: outer face up (+1) or down (-1). With no band only
    `flat` is read."""
    from types import SimpleNamespace

    from foldmap.cli import _membrane_places
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    lying, standing = np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])
    runs = {  # id: (kind, axis, CA depths)
        "face": ("H", lying, np.full(5, 19.9)),
        "beyond": ("H", lying, np.full(5, 20.1)),
        "deep": ("H", standing, np.full(5, -25.0)),
        "upright": ("H", standing, np.linspace(-6.0, 6.0, 5)),
        "strand": ("E", lying, np.full(5, 3.0)),
        "inner": ("H", lying, np.full(5, -6.0)),
        "turn": ("G", standing, np.full(3, 2.0)),
        "tm": ("H", standing, np.linspace(-16.0, 16.0, 9)),
    }
    ca, sses = [], []
    for name, (kind, axis, z) in runs.items():
        sses.append(SimpleNamespace(id=name, kind=kind, axis=axis, start=len(ca), end=len(ca) + len(z) - 1))
        ca += [(0.0, 0.0, float(v)) for v in z]
    bb = SimpleNamespace(ca=np.array(ca))
    flat, on_face, sides = _membrane_places(sses, bb, mem, {"tm"}, (3.45, -3.45))  # outer face up the page
    assert flat == {"face", "beyond", "inner"}
    assert on_face == {"face": 1, "upright": 0, "strand": 0, "inner": -1}
    assert sides == {"beyond": 1, "deep": -1}
    flat2, on_face2, sides2 = _membrane_places(sses, bb, mem, {"tm"}, (-3.45, 3.45))  # outer face down
    assert flat2 == flat and on_face2 == {"face": -1, "upright": 0, "strand": 0, "inner": 1}
    assert sides2 == {"beyond": -1, "deep": 1}
    assert _membrane_places(sses, bb, mem, {"tm"}, None) == (flat, {}, {})


def test_a_transmembrane_helix_half_the_band_long_is_drawn_a_hair_longer_and_crosses():
    """Ruling R6 (5LDW M:392-402, an 11-residue helix exactly half the band long): layout._stretch_short draws a
    transmembrane helix of at most half the band 2 * _PAST longer, so _to_face takes one end to a face and the other
    past the middle; a longer one, a lying one and one crossing no membrane are left as they are."""
    from foldmap.layout import _PAST, Placed, _Item, _stretch_short, _to_face, crosses_band
    from helpers import fake_sse

    band = (0.0, 6.9)
    shapes = {"half": (3.45, np.pi / 2), "long": (5.0, np.pi / 2), "lying": (3.0, 0.0), "other": (3.0, np.pi / 2)}
    items, members = [], {}
    for k, (name, (length, angle)) in enumerate(shapes.items()):
        s = fake_sse("H", 20 * k, 20 * k + 10)
        members[name] = Placed(s, 0, 0, length, angle, 1.0)
        items.append(_Item([members[name]], [(0.0, 0.0)], 0.0, 3.45, (0.0, 3.45), k, "A"))
    level = {members[n].sse.id for n in ("half", "long", "lying")}
    _stretch_short(items, level, band)
    assert members["half"].length == pytest.approx(3.45 + 2 * _PAST)
    box = items[0].box  # the item's box is read from its member's rect: it covers the stretched length (N3)
    assert box[3] - box[1] >= 3.45 + 2 * _PAST - 1e-9
    assert [members[n].length for n in ("long", "lying", "other")] == [5.0, 3.0, 3.0]
    shift = _to_face(items[0], level, band, 3.45)
    m = members["half"]
    assert crosses_band((3.45 + shift - m.length / 2, 3.45 + shift + m.length / 2), band)


def test_a_face_row_keeps_the_projected_order_and_stacks_outward_when_full():
    """Ruling R8 (layout._pack_row, _face_rows): a row packs its helices in their projected order, MARGIN apart, the
    least moved, round the columns standing in it; three flat helices seen on top of one another over one
    transmembrane helix cannot all lie within their chain's span on the face's first row, so one goes a row further
    out, the order along each row kept."""
    from foldmap.layout import MARGIN, _pack_row, build_layout
    from helpers import fake_sse

    at = _pack_row([0.0, 0.0], [1.0, 1.0], [(-np.inf, np.inf)])
    assert at == pytest.approx([-1.0 - MARGIN / 2, 1.0 + MARGIN / 2])
    at = _pack_row([-0.5, 0.5], [0.5, 0.5], [(-np.inf, -1.8), (1.8, np.inf)])  # a column from -1.8 to 1.8
    assert at[0] <= -2.3 + 1e-9 and at[1] >= 2.3 - 1e-9
    band = (-3.45, 3.45)
    tm = fake_sse("H", 0, 25, centroid=(0, 0, 0.0), axis=(0, 0, 1))
    flats = [
        fake_sse("H", 30 + 15 * k, 41 + 15 * k, centroid=(x, 0, 10.0 + k), axis=(1, 0, 0))
        for k, x in enumerate((-1, 0, 1))
    ]
    on = {s.id: 1 for s in flats}
    lay = build_layout(
        [tm, *flats], [], _page(), "projected", level={tm.id}, upright={tm.id}, flat=set(on), band=band, on_face=on
    )
    p = lay.placed
    rows = {}
    for s in flats:
        assert p[s.id].rect[1] >= band[1] + MARGIN - 1e-6
        rows.setdefault(round(p[s.id].rect[1], 3), []).append(s)
    assert len(rows) >= 2
    for row in rows.values():  # in each row, left to right in projected order
        assert [s.id for s in sorted(row, key=lambda s: p[s.id].cx)] == [s.id for s in row]


def test_transmembrane_columns_keep_their_projected_order():
    """Ruling R11 (review-2 I5, layout._in_order): each chain's items holding transmembrane elements are dealt back the
    places _attract found, left to right in their projected order; tied items and other chains' keep theirs."""
    from foldmap.layout import Placed, _in_order, _Item
    from helpers import fake_sse

    items = []
    for k, (chain, x) in enumerate((("A", 5.0), ("A", -5.0), ("A", 1.0), ("B", 9.0))):
        s = fake_sse("H", 10 * k, 10 * k + 20, chain=chain)
        items.append(_Item([Placed(s, 0, 0, 5.0, np.pi / 2, 1.0)], [(0.0, 0.0)], x, 0.0, (x, 0.0), k, chain))
    level = {it.members[0].sse.id for it in items}
    _in_order(items, level, [-3.0, 0.0, 2.0, 0.0], {})  # projected: item 0 leftmost of A, then 1, then 2
    assert [it.cx for it in items] == pytest.approx([-5.0, 1.0, 5.0, 9.0])
    assert [it.home[0] for it in items] == pytest.approx([-5.0, 1.0, 5.0, 9.0])


def test_a_band_side_name_moves_off_a_loop_crossing_it():
    """Ruling R9 (review-2 I1: 32KS's "out" crossed by the loop into NS1 at the face row's left end): a side name a loop
    crosses at the band's left end is drawn at its right end; with none crossing it, it stays (render._clear_of_lines,
    run on every line drawn once the page is done: _place_side_names)."""
    from types import SimpleNamespace

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    from foldmap.render import _clear_of_lines, _meets_box

    assert _meets_box((0.0, -1.0), (0.0, 1.0), (-0.5, -0.5, 0.5, 0.5))
    assert not _meets_box((1.0, -1.0), (1.0, 1.0), (-0.5, -0.5, 0.5, 0.5))
    ax = Figure().add_subplot()
    c = SimpleNamespace(pt_per_unit=20.0)
    t = ax.text(0.1, 5.35, "out", ha="left", va="bottom", fontsize=9)
    _clear_of_lines(c, t, [np.array([(0.4, 4.0), (0.4, 8.0)])], 0.1, 19.9, 1.0)
    assert t.get_position() == (19.9, 5.35) and t.get_horizontalalignment() == "right"
    t2 = ax.text(0.1, 5.35, "out", ha="left", va="bottom", fontsize=9)
    _clear_of_lines(c, t2, [np.array([(10.0, 4.0), (10.0, 8.0)])], 0.1, 19.9, 1.0)
    assert t2.get_position() == (0.1, 5.35) and t2.get_horizontalalignment() == "left"


def test_a_held_face_helix_stays_within_its_chain_span_through_the_layout():
    """Ruling R15 (re-review-2 N1: 6JLN e:42-47 ended 310 Å beyond its chain's columns after _attract and _separate):
    layout._out_of_room finds each held face helix outside its chain's padded transmembrane span as the columns now
    stand, and _into_room moves it sideways back to the nearer edge; a released one, or one inside, is left."""
    from foldmap.layout import _FACE_PAD, Placed, _into_room, _Item, _out_of_room
    from helpers import fake_sse

    tm = fake_sse("H", 0, 25)
    items = [_Item([Placed(tm, 0, 0, 6.0, np.pi / 2, 1.0)], [(0.0, 0.0)], 0.0, 0.0, (0.0, 0.0), 0, "A")]
    for k, x in enumerate((30.0, 2.0, -40.0)):
        s = fake_sse("H", 30 + 15 * k, 41 + 15 * k)
        items.append(_Item([Placed(s, 0, 0, 3.0, 0.0, 1.0)], [(0.0, 0.0)], x, 5.0, (x, 5.0), k + 1, "A"))
    faced = {1: (1.0, 4.0, 1), 2: (2.0, 4.0, 1), 3: (-1.0, 4.0, 1)}
    kept = np.array([True, True, True, False])
    level = {tm.id}
    assert _out_of_room(items, faced, kept, level) == [1]
    _into_room(items, faced, kept, level, [1])
    assert items[1].cx == pytest.approx(0.5 + _FACE_PAD) and items[1].home[0] == items[1].cx
    assert _out_of_room(items, faced, kept, level) == [] and items[3].cx == -40.0


# ------------------------------------------------------------- plan 26 Task 3: domains with a membrane
def test_auto_domains_join_the_domains_holding_the_band_and_leave_the_rest():
    """features.auto_domains `band` (ruling T3-a: contacts alone never split what the band holds; 32KS-A's competence
    domain came in five): the domains holding any of its elements are one, numbered as the first; the others stay as
    found; a chain whose every domain holds some is not split. On a committed fixture chain, 1AOZ_A (three domains)."""
    from types import SimpleNamespace

    from foldmap import domain_score as ds
    from foldmap.features import auto_domains

    inputs = ds.load_fixture(Path(__file__).parent / "data" / "domains" / "cath_inputs.npz")["1AOZ_A"]
    bb = SimpleNamespace(labels=inputs.labels, ca=inputs.ca)

    def auto(band=frozenset()):
        return auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb, band=band)

    (_, d1), (_, d2), (_, d3) = auto()
    assert auto({d2[0]}) == auto()  # a band in one domain joins nothing
    joined = auto({d1[-1], d3[0]})
    assert [name for name, _ in joined] == ["D1", "D2"]
    assert sorted(joined[0][1]) == sorted(d1 + d3) and joined[1][1] == d2
    assert auto({d1[0], d2[0], d3[0]}) == []


def test_the_band_joins_a_proteins_domains_in_every_copy_alike():
    """features._join_band: a domain holding a band element in any copy of a protein (`ref`, chain -> its reference)
    joins in every copy; the domains are numbered again from 0; a protein left with one domain is not split."""
    from foldmap.features import _join_band

    owner = {f"{c}{k}": (c, k) for c in "AB" for k in range(4)}
    owner |= {"C0": ("C", 0), "C1": ("C", 1)}
    _join_band(owner, {"A": "A", "B": "A", "C": "C"}, {"A0", "B2", "C0", "C1"})
    assert {k: d for k, (_, d) in owner.items()} == {
        "A0": 0,
        "A1": 1,
        "A2": 0,
        "A3": 2,
        "B0": 0,
        "B1": 1,
        "B2": 0,
        "B3": 2,
    }


def test_a_soluble_domain_takes_the_side_its_centroid_lies_on():
    """cli._domain_sides (ruling T3-c) on a 30 Å slab, normal up: a domain holding a transmembrane or face element
    (`held`) is none; any other goes by its residues' centroid, above or below, whatever its elements' own sides (MIX:
    one 30 Å above, one 10 Å below); the sign follows the page, outer face up (+1) or down (-1)."""
    from types import SimpleNamespace

    from foldmap.cli import _domain_sides
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    runs = {"tm": np.linspace(-16, 16, 9), "face": np.full(5, 17.0), "up": np.full(5, 25.0), "down": np.full(5, -30.0)}
    runs |= {"a": np.full(5, 30.0), "b": np.full(5, -10.0)}
    ca, sses = [], []
    for name, z in runs.items():
        sses.append(SimpleNamespace(id=name, start=len(ca), end=len(ca) + len(z) - 1))
        ca += [(0.0, 0.0, float(v)) for v in z]
    bb = SimpleNamespace(ca=np.array(ca))
    named = [("M", ["tm", "up"]), ("F", ["face"]), ("UP", ["up"]), ("DOWN", ["down"]), ("MIX", ["a", "b"])]
    assert _domain_sides(named, sses, bb, mem, {"tm", "face"}, (3.45, -3.45)) == {"UP": 1, "DOWN": -1, "MIX": 1}
    assert _domain_sides(named, sses, bb, mem, {"tm", "face"}, (-3.45, 3.45)) == {"UP": -1, "DOWN": 1, "MIX": -1}


def _domain_item(s, x, y, length, angle):
    from foldmap.layout import Placed, _Item

    return _Item([Placed(s, 0, 0, length, angle, 1.0)], [(0.0, 0.0)], x, y, (x, y), s.start, s.chain)


@pytest.mark.parametrize("case", ["band-soluble", "band-band", "soluble-soluble"])
def test_domains_part_by_the_band(case):
    """layout._part_domains in a membrane view (`kept`, `sides`; rulings T3-a to T3-c; 32KS-A's competence domain was
    parted 1-2 Å off the band, four TM helices short): a soluble domain moves off a domain in the band outward on its
    side, at its own x, and the band's domain stays; two domains in the band part sideways only, though up or down is
    cheaper; of two soluble domains on one side, the outer moves on outward and the inner stays."""
    from foldmap.layout import _DOMAIN_GAP, _part_domains
    from helpers import fake_sse

    a, b = fake_sse("H", 0, 20), fake_sse("H", 30, 40)
    up, flat = np.pi / 2, 0.0
    if case == "band-soluble":  # a TM column, and a soluble helix lying over its top
        items = [_domain_item(a, 0.0, 0.0, 8.0, up), _domain_item(b, 0.5, 1.0, 6.0, flat)]
        kept, sides = np.array([True, False]), [0, 1]
    elif case == "band-band":  # a TM column and a face helix held over its top, of two domains in the band
        items = [_domain_item(a, 0.0, 0.0, 8.0, up), _domain_item(b, 0.5, 4.5, 8.0, flat)]
        kept, sides = np.array([True, True]), [0, 0]
    else:  # two soluble helices above the band, the second a little higher
        items = [_domain_item(a, 0.0, 6.0, 6.0, flat), _domain_item(b, 0.5, 7.0, 6.0, flat)]
        kept, sides = np.array([False, False]), [1, 1]
    before = [(it.cx, it.cy) for it in items]
    assert _part_domains(items, [("P", [a.id]), ("Q", [b.id])], kept=kept, sides=sides)
    (x0, y0), (x1, y1) = before
    if case == "band-soluble":
        assert (items[0].cx, items[0].cy) == (x0, y0) and items[1].cx == x1 and items[1].cy > y1
        assert items[1].cy + items[1].box[1] >= items[0].cy + items[0].box[3] + _DOMAIN_GAP - 1e-6
    elif case == "band-band":
        assert items[0].cy == y0 and items[1].cy == y1 and abs(items[1].cx - items[0].cx) > abs(x1 - x0)
    else:
        assert (items[0].cx, items[0].cy) == (x0, y0) and items[1].cx == x1 and items[1].cy > y1


def test_a_soluble_domains_panel_stays_out_of_the_slab():
    """build_layout `domain_sides` (ruling T3-c): a soluble domain's elements keep off the band far enough that its
    panel (domain_panel: DOMAIN_PAD round its elements, its name's room above) stays clear of the slab, on either side;
    with the elements walled only MARGIN clear (no `domain_sides`), the lower domain's panel reaches into it."""
    from foldmap.layout import build_layout, domain_panel
    from helpers import fake_sse

    band = (-3.45, 3.45)
    tm = fake_sse("H", 0, 25, centroid=(0, 0, 0.0), axis=(0, 0, 1))
    hi = fake_sse("H", 30, 41, centroid=(60, 0, 2.0), axis=(1, 0, 0))
    lo = fake_sse("H", 45, 56, centroid=(-60, 0, -2.0), axis=(0, 0, 1))  # upright, its label beside it, not above
    tail = fake_sse("H", 70, 80, centroid=(120, 0, 30.0), axis=(1, 0, 0))  # the C-terminus: no end letter on `lo`
    doms = [("M", [tm.id]), ("HI", [hi.id, tail.id]), ("LO", [lo.id])]
    sides = {hi.id: 1, tail.id: 1, lo.id: -1}
    kw = dict(level={tm.id}, upright={tm.id}, flat={hi.id, tail.id}, band=band, sides=sides, domains=doms)
    els = [tm, hi, lo, tail]
    lay = build_layout(els, [], _page(), "projected", domain_sides={"HI": 1, "LO": -1}, **kw)
    assert domain_panel(lay, "HI")[1] >= band[1] and domain_panel(lay, "LO")[3] <= band[0]
    plain = build_layout(els, [], _page(), "projected", **kw)
    assert domain_panel(plain, "LO")[3] > band[0]


@pytest.mark.parametrize("step", ["parting", "widening"])
def test_the_face_rows_bound_holds_after_the_domains_are_parted(monkeypatch, step):
    """Ruling R15 through the whole layout (re-review-3's parked note): a held face helix that a step after the domains'
    parting, or a hand-widened element's separation, carries 40 units off its chain's columns is brought back within
    its padded transmembrane span (or released, at its projected place), as after _attract and _separate. The 40 are
    put down to the step itself, the room not carried by the parting (_parted_rooms, review-3 M2: tested below)."""
    import foldmap.layout as L
    from foldmap.style import Style
    from helpers import fake_sse

    band = (-3.45, 3.45)
    tm = fake_sse("H", 0, 25, centroid=(0, 0, 0.0), axis=(0, 0, 1))
    face = fake_sse("H", 30, 41, centroid=(5, 0, 10.0), axis=(1, 0, 0))

    def away(items, *a, **k):  # the step moves the face helix far off, as a big parting or widening may
        for it in items:
            if any(m.sse.id == face.id for m in it.members):
                it.cx += 40.0
                it.home = (it.cx, it.home[1])
        return True

    monkeypatch.setattr(L, "_part_domains" if step == "parting" else "_widen", away)
    monkeypatch.setattr(L, "_parted_rooms", lambda *a, **k: {})  # the rooms not carried (review-3 M2: tested below)
    lay = L.build_layout(
        [tm, face],
        [],
        _page(),
        "projected",
        level={tm.id},
        upright={tm.id},
        flat={face.id},
        band=band,
        on_face={face.id: 1},
        domains=[("M", [tm.id]), ("F", [face.id])],
        domain_sides={},
        style=Style(element_widths="α1=2") if step == "widening" else None,
    )
    p = lay.placed
    pad = max(L._FACE_PAD, p[face.id].length)
    assert p[tm.id].rect[0] - pad - 1e-6 <= p[face.id].cx <= p[tm.id].rect[2] + pad + 1e-6


@pytest.mark.skipif(not KS32.is_file(), reason="32KS is not in the download cache")
def test_comec_auto_domains_in_the_membrane_view_are_the_papers_three():
    """32KS --chains A --membrane auto --domains auto (Deselaers et al., Fig. 1C): three domains, the competence domain
    (every TM helix and the NS helices on the face, A:5-526 less the OB fold), the OB fold (A:87-195) and the
    β-lactamase-like domain (A:535-795); A18 passes, the competence panel spans the band and the other two lie above it
    (before: eight panels, four TM helices parted off the band)."""
    from foldmap import criteria
    from foldmap.layout import domain_panel
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(KS32, membrane="auto", chains=["A"], domains="auto", style=["color_by=domain"]))
    lay = d.layout
    first = {name: min(lay.placed[k].sse.first.seq for k in ids) for name, ids in lay.domains}
    last = {name: max(lay.placed[k].sse.last.seq for k in ids) for name, ids in lay.domains}
    assert (first, last) == ({"D1": 5, "D2": 87, "D3": 535}, {"D1": 526, "D2": 195, "D3": 795})
    m = criteria.membrane_checks(d, KS32)
    assert m["band"] and not m["tm_short"], m  # and the rest of A18, strict: test_a18_on_comec
    lo, hi = lay.membrane["y"]
    assert domain_panel(lay, "D1")[1] <= lo and domain_panel(lay, "D1")[3] >= hi
    assert domain_panel(lay, "D2")[1] >= hi and domain_panel(lay, "D3")[1] >= hi


def test_a_membrane_view_with_face_helices_still_counts_the_copies_it_let_go(monkeypatch):
    """Ruling R18: 76a7d7d's R15 loop rebound build_layout's `gone` (the copy items _one_per_mate let go) to its own
    last, empty, list, so Layout.let_go read 0 on every membrane view with a face helix and cli.make_layout never
    laid such a view out again without letting go (the A15 retry). Here one copy item is let go (the ties and the
    choice faked) beside a helix held on a face: the layout counts it."""
    import foldmap.layout as L
    from helpers import fake_sse

    band = (-3.45, 3.45)
    tm = fake_sse("H", 0, 25, centroid=(0, 0, 0.0), axis=(0, 0, 1))
    face = fake_sse("H", 30, 41, centroid=(5, 0, 10.0), axis=(1, 0, 0))
    copy = fake_sse("H", 0, 25, centroid=(40, 0, 0.0), axis=(0, 0, 1), chain="B")

    def tied(items, *a, **k):  # B's helix tied to A's, protomer 1
        at = {m.sse.id: n for n, it in enumerate(items) for m in it.members}
        return {at[copy.id]: (at[tm.id], 1)}

    monkeypatch.setattr(L, "_symmetry_ties", tied)
    monkeypatch.setattr(L, "_one_per_mate", lambda items, ties, **k: [])  # and let go
    lay = L.build_layout(
        [tm, face, copy],
        [],
        _page(),
        "projected",
        level={tm.id, copy.id},
        upright={tm.id, copy.id},
        flat={face.id},
        band=band,
        on_face={face.id: 1},
        alike_first=False,
    )
    assert lay.let_go == 1


# ------------------------------------------------- plan 26 Task 3 fix round 1 (review-3: R20, M1-M5)
@pytest.mark.parametrize("case", ["carried", "aside"])
def test_a_soluble_domain_keeps_to_its_chains_membrane_panel(case):
    """Ruling R20 (layout._part_way's `follow`, `own`): when two chains' band domains part sideways, a soluble domain of
    a chain moves with that chain's band domain (it stayed at its absolute x); a soluble domain overlapping another
    chain's band domain moves sideways off it where that is cheaper, its height kept (it was only lifted)."""
    from foldmap.layout import _part_domains
    from helpers import fake_sse

    up, flat = np.pi / 2, 0.0
    a, b = fake_sse("H", 0, 20), fake_sse("H", 0, 20, chain="B")
    s = fake_sse("H", 30, 40)  # chain A's soluble helix
    if case == "carried":  # A's and B's columns side by side, overlapping; A's soluble helix high above A's
        items = [
            _domain_item(a, 0.0, 0.0, 8.0, up),
            _domain_item(b, 1.0, 0.0, 8.0, up),
            _domain_item(s, -0.5, 9.0, 4.0, flat),
        ]
    else:  # A's columns far left; A's soluble helix over the edge of B's column, deeper than it is wide
        items = [
            _domain_item(a, -40.0, 0.0, 8.0, up),
            _domain_item(b, 0.0, 0.0, 8.0, up),
            _domain_item(s, 3.0, 1.0, 4.0, flat),
        ]
    before = [(it.cx, it.cy) for it in items]
    doms = [("A", [a.id]), ("B", [b.id]), ("S", [s.id])]
    assert _part_domains(items, doms, kept=np.array([True, True, False]), sides=[0, 0, 1])
    if case == "carried":
        assert items[0].cx != before[0][0] and items[2].cx - before[2][0] == pytest.approx(items[0].cx - before[0][0])
        assert items[2].cy == before[2][1]
    else:
        assert items[2].cy == before[2][1] and items[2].cx > before[2][0]
        assert items[0].cx == before[0][0] and items[1].cx == before[1][0]


def test_a_face_domain_parted_off_its_columns_is_held_where_it_was_parted(monkeypatch):
    """Review-3 M2: a user's domain of face helices parted sideways off its chain's TM panel is held within its room
    carried as far as the parting moved it (_parted_rooms), not pulled back over the panel by R15."""
    import foldmap.layout as L
    from helpers import fake_sse

    band = (-3.45, 3.45)
    tm = fake_sse("H", 0, 25, centroid=(0, 0, 0.0), axis=(0, 0, 1))
    face = fake_sse("H", 30, 41, centroid=(5, 0, 10.0), axis=(1, 0, 0))

    def apart(items, *a, **k):  # the face domain parted 40 units off the TM domain, as a wide panel would push it
        for it in items:
            if any(m.sse.id == face.id for m in it.members):
                it.cx += 40.0
                it.home = (it.cx, it.home[1])
        return True

    monkeypatch.setattr(L, "_part_domains", apart)
    lay = L.build_layout(
        [tm, face],
        [],
        _page(),
        "projected",
        level={tm.id},
        upright={tm.id},
        flat={face.id},
        band=band,
        on_face={face.id: 1},
        domains=[("TM", [tm.id]), ("NS", [face.id])],
        domain_sides={},
    )
    p = lay.placed
    assert p[face.id].cx > p[tm.id].rect[2] + max(L._FACE_PAD, p[face.id].length) + 20.0


def test_the_note_says_when_the_band_joined_every_domain():
    """Review-3 M1 (ruling R19): when residue contacts split a chain but the band joins every domain (5JDN's two TM
    halves), the note says no domains are drawn and why, not that contacts found one domain."""
    from foldmap.cli import domain_note

    note = domain_note("inferred", [], band=True, absorbed=True)
    assert note.startswith("no domains drawn") and "one domain per chain" in note and "reads as one domain" not in note
    assert domain_note("inferred", [], band=True).startswith("no domains inferred: each chain reads as one domain")
    assert (
        domain_note("inferred", [], band=False)
        == "no domains inferred: each chain reads as one domain from its residue contacts"
    )


def test_the_band_holds_a_crossing_subunit_that_does_not_define_the_membrane():
    """Review-3 M3 and M4 (cli._band_ids): 9FWZ's FimA subunit (chain B) lies in the usher's lumen, crossing the slab
    without defining the membrane (a beta-barrel's strands do): it is in no transmembrane set, but the band holds it,
    so auto domains never split it off and a domain of it is not soluble. The face and slab ids of _membrane_places
    are the same whichever way the faces lie on the page (the contract _band_ids asserts)."""
    from foldmap.cli import _band_ids, _domain_sides, _membrane_places, _transmembrane

    lay, sses, bb = make_layout(DATA / "r100" / "9FWZ-DB.cif.gz", membrane="auto")
    mem = lay.bilayer
    tm = _transmembrane(sses, bb, mem)
    held = _band_ids(sses, bb, mem, tm)
    lumen = {"B:52-61", "B:68-74", "B:127-138", "B:153-158"}
    assert lumen <= held and not lumen & tm
    up, down = (_membrane_places(sses, bb, mem, tm, f)[1] for f in ((1.0, -1.0), (-1.0, 1.0)))
    assert up.keys() == down.keys() and all(up[k] == -down[k] for k in up)
    assert _domain_sides([("FimA", sorted(lumen))], sses, bb, mem, held, (3.45, -3.45)) == {}


def test_a_comec_like_protein_draws_the_three_domains_of_the_paper():
    """Review-3 M5, on CI: a synthetic ComEC (one chain; 30 Å slab, normal up): three TM helices, a soluble block above
    them, three TM helices with a helix lying on the outer face, and another soluble block above, in sequence. Residue
    contacts find four domains (the two TM groups 40 Å apart); the band joins the TM groups and the face helix
    (features.auto_domains `band`, cli._band_ids), so three are drawn, as the paper's competence, OB and β-lactamase
    domains (cli._membrane_args, layout.build_layout): the competence panel across the band, every TM helix across its
    middle, the two soluble panels above the outer face, no two panels overlapping."""
    from types import SimpleNamespace

    from foldmap.cli import _band_ids, _faces_on_page, _membrane_args, _transmembrane
    from foldmap.features import auto_domains
    from foldmap.layout import build_layout, domain_panel
    from foldmap.membrane import Membrane
    from foldmap.model import ResLabel
    from helpers import fake_sse

    ca, sses = [], []

    def helix(start, end, n=24):
        start, end = np.asarray(start, float), np.asarray(end, float)
        first = len(ca)
        ca.extend(start + (end - start) * t for t in np.linspace(0.0, 1.0, n))
        axis = (end - start) / np.linalg.norm(end - start)
        sses.append(fake_sse("H", first, len(ca) - 1, centroid=(start + end) / 2, axis=axis))
        return sses[-1].id

    tm1 = [helix((x, 0, -17.25 * sg), (x, 0, 17.25 * sg)) for x, sg in ((0, 1), (10, -1), (20, 1))]
    ob = [helix((0, 0, z), (28.5, 0, z), 20) for z in (32, 37, 42)]
    tm2 = [helix((x, 0, -17.25 * sg), (x, 0, 17.25 * sg)) for x, sg in ((60, 1), (70, -1), (80, 1))]
    face = helix((56, 6, 17), (84.5, 6, 17), 20)
    blact = [helix((60, 0, z), (88.5, 0, z), 20) for z in (32, 37, 42)]

    bb = SimpleNamespace(ca=np.array(ca), labels=[ResLabel("A", k, "", "ALA") for k in range(len(ca))])
    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    tm = _transmembrane(sses, bb, mem)
    assert tm == set(tm1 + tm2)
    assert len(auto_domains(sses, {}, [], {}, frozenset(), bb=bb)) == 4  # contacts alone: four
    in_band = _band_ids(sses, bb, mem, tm)
    named = auto_domains(sses, {}, [], {}, frozenset(), bb=bb, band=in_band)
    assert [(n, set(ids)) for n, ids in named] == [("D1", set(tm1 + tm2 + [face])), ("D2", set(ob)), ("D3", set(blact))]
    args = _membrane_args(sses, bb, mem, tm, _faces_on_page(_page(), sses, bb, mem, tm), named, in_band)
    assert args["domain_sides"] == {"D2": 1, "D3": 1}
    lay = build_layout(sses, [], _page(), "projected", domains=named, **args)
    lo, hi = args["band"]
    panel = {n: domain_panel(lay, n) for n, _ in named}
    assert panel["D1"][1] < lo and panel["D1"][3] > hi
    assert panel["D2"][1] >= hi and panel["D3"][1] >= hi
    names = sorted(panel)
    for i, m in enumerate(names):
        for n in names[i + 1 :]:
            P, Q = panel[m], panel[n]
            assert P[2] <= Q[0] or Q[2] <= P[0] or P[3] <= Q[1] or Q[3] <= P[1], (m, n)
    for k in tm:
        r = lay.placed[k].rect
        assert r[1] < (lo + hi) / 2 < r[3]


# ------------------------------------------------- plan 26 Task 5: joins, the face's flat cut, TM1, A18, the turn
def _straight(points, chain="A"):
    """A backbone as membrane.joins and ss._copy_breaks read it, one residue per CA point (given as runs of points):
    CA atoms, labels (chain, number from 1) and every residue bonded to the one before within its chain."""
    from types import SimpleNamespace

    from foldmap.model import ResLabel

    ca = np.array(points, float)
    labels = [ResLabel(chain, k + 1, "", "ALA") for k in range(len(ca))]
    prev = np.array([k > 0 for k in range(len(ca))])
    return SimpleNamespace(ca=ca, labels=labels, prev=prev)


def _run(start, direction, n, step=1.5):
    start, direction = np.asarray(start, float), np.asarray(direction, float)
    return [start + direction / np.linalg.norm(direction) * step * k for k in range(n)]


def test_flat_helix_pieces_a_break_apart_are_joined_only_going_on_as_one():
    """Ruling T5-1, item 1 (32KS's NS7, A:436-461, broken by DSSP at 447-448, its pieces 37° apart): two helix runs
    two residues apart, both lying flat at the outer face, are one helix in a membrane view when one goes on as the
    other (membrane.JOIN_ANGLE); not when the second folds back (a hairpin, 32KS's NS1 and NS2, 151°), nor three
    residues apart, nor when one stands."""
    from foldmap.membrane import JOIN_ANGLE, Membrane, joins

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")

    def pair(turn_deg, gap=2, second=(None,)):
        first = _run((0, 0, 17), (1, 0, 0), 11)
        end = first[-1]
        t = np.radians(turn_deg)
        d = np.array([np.cos(t), np.sin(t), 0.0]) if second[0] is None else np.asarray(second, float)
        linker = _run(end + d * 1.5, d, gap)
        rest = _run(linker[-1] + d * 1.5, d, 13)
        bb = _straight(first + linker + rest)
        return bb, (0, 10), (11 + gap, 11 + gap + 12)

    assert JOIN_ANGLE == 45.0

    def joined(bb, x, y):
        return joins(bb, mem, x, y)

    assert joined(*pair(37.0))
    assert not joined(*pair(151.0))
    assert not joined(*pair(37.0, gap=3))
    assert not joined(*pair(0.0, second=(0.3, 0.0, 1.0)))


def test_join_runs_joins_each_next_run_its_predicate_names():
    """ss._join_runs: three runs, the predicate joining the first two: one run of both, reported as made; the third
    kept apart."""
    from foldmap.ss import _join_runs

    runs = [(0, 10), (13, 25), (30, 40)]
    out, made = _join_runs(runs, lambda x, y: y[0] == 13)
    assert out == [(0, 25), (30, 40)] and made == [(0, 25)]
    assert _join_runs(runs, lambda x, y: False) == (runs, [])


def test_copies_broken_where_the_joined_helix_crosses_are_joined_not_cut():
    """Ruling T5-1, item 3 (32KS: chain B's TM8 broken at B:364, chain A's whole): one copy whole and one broken by a
    residue tie, so _copy_breaks cuts the whole one; with `whole_if` (membrane.crosses_joined: the joined helix
    crosses the membrane) the broken copy is joined instead, both drawn whole."""
    from foldmap.membrane import Membrane, crosses_joined
    from foldmap.model import ResLabel
    from foldmap.ss import _copy_breaks

    helix = _run((0, 0, -18), (0, 0, 1), 25)
    bb = _straight(helix + [p + np.array([30.0, 0, 0]) for p in helix])
    bb.labels = [ResLabel("A" if k < 25 else "B", k % 25 + 1, "", "ALA") for k in range(50)]
    bb.prev = np.array([k % 25 > 0 for k in range(50)])
    copies = lambda: {"A": ("A", 0), "B": ("A", 0)}  # noqa: E731
    runs = [(0, 24), (25, 37), (39, 49)]
    cut, n_cut, _, _ = _copy_breaks(runs, bb, copies)
    assert n_cut == 1 and (0, 12) in cut and (14, 24) in cut
    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    whole, n_cut, n_joined, _ = _copy_breaks(runs, bb, copies, lambda x, y: crosses_joined(bb, mem, x, y))
    assert n_cut == 0 and n_joined == 1 and whole == [(0, 24), (25, 49)]


def test_a_helix_at_a_face_lies_flat_from_35_degrees():
    """Ruling T5-2 (32KS's NS2 and NS3, tilted 44.6° and 38.6°, centres 2.9 and 6.1 Å beyond the outer face, stood
    upright above it): in a face's interfacial zone (membrane.INTERFACE either side of it) a helix lies flat from
    FACE_FLAT_TILT; elsewhere from FLAT_TILT; a crossing helix always stands."""
    from foldmap.membrane import FACE_FLAT_TILT, FLAT_TILT, INTERFACE, Membrane, lies_flat_here

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")

    def at(deg, centre, n=5):
        axis = np.array([np.sin(np.radians(deg)), 0.0, np.cos(np.radians(deg))])
        return axis, centre + np.linspace(-1.0, 1.0, n) * np.cos(np.radians(deg)) * 3

    assert (FACE_FLAT_TILT, FLAT_TILT, INTERFACE) == (35.0, 45.0, 7.5)
    assert lies_flat_here(*at(38.6, 21.1), mem) and lies_flat_here(*at(36.0, -9.0), mem)
    assert not lies_flat_here(*at(34.0, 21.1), mem)
    assert not lies_flat_here(*at(40.0, 0.0), mem) and not lies_flat_here(*at(40.0, 30.0), mem)
    assert lies_flat_here(*at(46.0, 0.0), mem)
    axis = np.array([np.sin(np.radians(50.0)), 0.0, np.cos(np.radians(50.0))])
    assert not lies_flat_here(axis, np.linspace(-17.0, 17.0, 12), mem)  # it crosses: it stands


def test_a_tm_helix_drawn_between_the_faces_goes_to_the_face_its_run_comes_nearer():
    """Ruling T5-4 (32KS's TM1, A:5-19: its run -11.0 to +9.2 Å, drawn 21 Å long about its CA centre at +0.1 Å, went
    4.4 Å out rather than 4.6 Å in, its N end 5 Å off): a helix drawn wholly inside the band goes to the face its CA
    run comes nearer (cli._nearer_faces, layout._to_face `nearer`), not the nearer by the least move."""
    from types import SimpleNamespace

    from foldmap.cli import _nearer_faces
    from foldmap.layout import _to_face, crosses_band
    from foldmap.membrane import Membrane

    band, ends = (0.0, 10.0), (1.6, 7.6)  # the lower face 1.6 away, the upper 2.4
    it = _upright_item([ends])
    assert _to_face(it, {"t0"}, band, 0.0) == pytest.approx(-1.6)  # the least move
    shift = _to_face(it, {"t0"}, band, 0.0, {"t0": 1})
    assert shift == pytest.approx(2.4) and crosses_band((ends[0] + shift, ends[1] + shift), band)
    assert _to_face(it, {"t0"}, band, 0.0, {"t0": -1}) == pytest.approx(-1.6)
    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    z = np.array([-11.0, -7.9, -6.2, -5.9, -4.6, -1.9, -0.9, -0.6, 1.6, 3.8, 4.3, 5.0, 7.7, 9.1, 9.2])
    bb = SimpleNamespace(ca=np.stack([np.zeros(15), np.zeros(15), z], axis=1))
    tm = [SimpleNamespace(id="TM1", start=0, end=14)]
    assert _nearer_faces(tm, bb, mem, {"TM1"}, (3.45, -3.45)) == {"TM1": -1}  # outer face up: the lower face
    assert _nearer_faces(tm, bb, mem, {"TM1"}, (-3.45, 3.45)) == {"TM1": 1}


def test_an_element_wholly_beyond_a_face_is_kept_off_the_band():
    """Ruling T5-5's layout (6ZKT's M:50-51, 18-19 Å out, drawn 10 Å into the slab): cli._membrane_places puts an
    element standing wholly beyond a face, its centre within SOLUBLE_BEYOND, among the soluble ones kept off the band
    on that side (`sides`); a flat one there still lies on the face, and one reaching into the slab is placed as
    before."""
    from types import SimpleNamespace

    from foldmap.cli import _membrane_places
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    runs = {
        "strand": ("E", np.array([0.0, 0.0, 1.0]), np.linspace(18.3, 19.3, 3)),
        "flat": ("H", np.array([1.0, 0.0, 0.0]), np.full(5, 17.0)),
        "reaching": ("E", np.array([0.0, 0.0, 1.0]), np.linspace(14.0, 18.0, 3)),
    }
    ca, sses = [], []
    for name, (kind, axis, z) in runs.items():
        sses.append(SimpleNamespace(id=name, kind=kind, axis=axis, start=len(ca), end=len(ca) + len(z) - 1))
        ca += [(0.0, 0.0, float(v)) for v in z]
    _, on_face, sides = _membrane_places(sses, SimpleNamespace(ca=np.array(ca)), mem, set(), (3.45, -3.45))
    assert sides == {"strand": 1} and on_face == {"flat": 1}


def _band_helices(order_x, spread=100.0):
    """Six helices crossing a 30 Å slab (normal z) at y = 10 k in sequence and x = order_x[k], with two soluble
    helices `spread` Å apart along x above it, so view_frame's page-right is x."""
    from types import SimpleNamespace

    from helpers import fake_sse

    ca, sses = [], []
    for k, x in enumerate(order_x):
        first = len(ca)
        ca += [(x, 10.0 * k, z) for z in np.linspace(-18.0, 18.0, 25)]
        sses.append(fake_sse("H", first, len(ca) - 1, centroid=(x, 10.0 * k, 0.0), axis=(0, 0, 1)))
    for x in (-spread / 2, spread / 2):
        first = len(ca)
        ca += [(x + t, 25.0, 40.0) for t in np.linspace(-9.0, 9.0, 13)]
        sses.append(fake_sse("H", first, len(ca) - 1, centroid=(x, 25.0, 40.0), axis=(1, 0, 0)))
    return sses, SimpleNamespace(ca=np.array(ca))


def test_a_membrane_view_turns_about_the_normal_to_stand_its_helices_in_sequence():
    """Ruling T5-3 (32KS chain A: TM1 in the middle of the band, tau 0.47): of the turns about the normal in
    AZIMUTH_STEP steps, cli._in_sequence takes the one whose page x order of the band helices agrees best with their
    sequence: helices in sequence along y, seen along x in no order, are turned until all stand in order; the frame
    stays a rotation with page-up the normal. Helices already in order keep the frame as chosen."""
    from foldmap.cli import _in_sequence, _sequence_score
    from foldmap.frame import view_frame
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    sses, bb = _band_helices([3.0, -2.0, 1.0, -3.0, 2.0, 0.0])
    frame = view_frame(sses, up=mem.normal)
    cents = [[s.centroid for s in sses[:6]]]
    assert _sequence_score(frame, cents) < 15
    turned = _in_sequence(frame, sses, bb, mem)
    assert _sequence_score(turned, cents) == 15  # tau 1
    assert np.allclose(turned.v, frame.v) and np.allclose(np.cross(turned.u, turned.v), turned.w)
    assert turned is not frame
    ordered, bb2 = _band_helices([-5.0, -3.0, -1.0, 1.0, 3.0, 5.0])
    plain = view_frame(ordered, up=mem.normal)
    assert _in_sequence(plain, ordered, bb2, mem) is plain


def test_a_turn_worth_less_than_a_tenth_in_tau_is_not_taken(monkeypatch):
    """Ruling T5-3's floor (MIN_TAU_GAIN_INV, 0.1): fifteen pairs, one of them out of order as the frame was chosen;
    a turn puts all in order, tau +2/15 (0.13): taken. With the floor at 1/6 (0.17) the frame as chosen is kept."""
    import foldmap.cli as cli
    from foldmap.frame import view_frame
    from foldmap.membrane import Membrane

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    sses, bb = _band_helices([-5.0, -3.0, 1.0, -1.0, 3.0, 5.0])  # x order: one pair swapped
    frame = view_frame(sses, up=mem.normal)
    cents = [[s.centroid for s in sses[:6]]]
    best = cli._in_sequence(frame, sses, bb, mem)
    assert cli._sequence_score(frame, cents) == 13 and cli._sequence_score(best, cents) == 15
    monkeypatch.setattr(cli, "MIN_TAU_GAIN_INV", 6)
    assert cli._in_sequence(frame, sses, bb, mem) is frame


@pytest.mark.skipif(not KS32.is_file(), reason="32KS is not in the download cache")
def test_comec_reads_as_its_figure_1c():
    """Plan 26 Task 5 gates on 32KS (ComEC, Fig. 1C): with --chains A, NS7 (A:436-461) is one helix lying flat, NS2 and
    NS3 (A:222-233, A:241-250) lie flat, TM1 (A:5-19) reaches the inner face, its N terminus below the band, and the
    band helices read more nearly in sequence (tau over 0.7); in the dimer chain A's TM8 (A:349-371) is one helix
    across the band, and nothing beyond a face is drawn in the slab."""
    from foldmap.figcheck import membrane_view
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(KS32, membrane="auto", chains=["A"])
    p = lay.placed
    for k in ("A:436-461", "A:222-233", "A:241-250"):
        assert k in p and abs(np.sin(p[k].angle)) < 1e-9, k
    from foldmap.cli import _transmembrane

    lo, hi = lay.membrane["y"]
    tm1 = sorted((p["A:5-19"].n_port[1], p["A:5-19"].c_port[1]))
    assert tm1[0] == pytest.approx(lo, abs=1e-3) and p["A:5-19"].n_port[1] < p["A:5-19"].c_port[1]
    tm = _transmembrane(sses, bb, lay.bilayer)
    x = [p[s.id].cx for s in sses if s.id in tm and s.kind == "H"]
    pairs = [(i, j) for i in range(len(x)) for j in range(i + 1, len(x))]
    assert sum(np.sign(x[j] - x[i]) for i, j in pairs) / len(pairs) > 0.7
    dimer, dsses, dbb = make_layout(KS32, membrane="auto")
    assert "A:349-371" in dimer.placed and "B:349-371" in dimer.placed
    r = dimer.placed["A:349-371"].rect
    assert r[1] < sum(dimer.membrane["y"]) / 2 < r[3]
    view = membrane_view(dimer, dsses, dbb, dimer.bilayer, route_loops(dimer, dsses, dbb))
    assert view["tm_short"] == []  # and nothing beyond a face in the slab, A18 strict: test_a18_on_comec


COMEC_CASES = {  # the three ComEC figures A18 is held on (Fig. 1C): the dimer, chain A, chain A with its domains
    "dimer": {},
    "A": {"chains": ["A"]},
    "A-domains": {"chains": ["A"], "domains": "auto", "style": ["color_by=domain"]},
}
COMEC_PASSING = {"wrong_side"}  # A18's parts ComEC's figures pass today: they guard 32KS until Tasks 4 and 5


@functools.cache
def _comec_a18(case: str) -> dict:
    """criteria.membrane_checks on one ComEC figure (COMEC_CASES), drawn once per process."""
    from foldmap import criteria
    from foldmap.spec import draw_spec, figure_spec

    return criteria.membrane_checks(draw_spec(figure_spec(KS32, membrane="auto", **COMEC_CASES[case])), KS32)


@pytest.mark.skipif(not KS32.is_file(), reason="32KS is not in the download cache")
@pytest.mark.parametrize(
    "case, part",
    [
        pytest.param(case, part, marks=() if part in COMEC_PASSING else PENDING, id=f"{case}-{part}")
        for case in COMEC_CASES
        for part in ("slab_helix", "band_member", "strand_in_band", "outside_in_band", "wrong_side", "loop_in_band")
    ],
)
def test_a18_on_comec(case, part):
    """Plan 27 Task 2 (A18 strict) on 32KS, ComEC, the figure the author judged the membrane view by (Fig. 1C): in the
    dimer, chain A alone and chain A with its domains, each part of A18 passes on its own (criteria.membrane_checks
    on the drawn figure; review-2 M4: marked per part, so a Task's gate shows its own parts flip): helices in the slab
    at their depths, only the band's elements in the band, no β-strand there (the OB barrels' G/H hairpins), none
    nearer than its residues or on the wrong side, no loop nearer than its residues."""
    m = _comec_a18(case)
    assert m["band"] and not m["tm_short"], m
    assert m[part] == [], m[part]


def test_a_two_fold_dimer_about_the_normal_keeps_its_frame():
    """Ruling R23 (32KS: turned 55° for the pair's summed tau, its protomers were drawn through each other): a turn is
    judged by its worst copy (cli._worst_tau). Two copies of a protein, the second the first turned 180° about the
    normal and one band helix short (32KS chain B has no TM11): any turn that stands one in sequence reverses the
    other, so the frame as chosen is kept; summed over the chains, the turn would have won (15 - 10 pairs)."""
    from types import SimpleNamespace

    from foldmap.cli import _in_sequence, _sequence_score, _worst_tau
    from foldmap.frame import Frame, view_frame
    from foldmap.membrane import Membrane
    from helpers import fake_sse

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    ca, sses = [], []
    for chain, sign, n in (("A", 1.0, 6), ("B", -1.0, 5)):
        for k in range(n):
            first, y = len(ca), sign * (10.0 * k + 5.0)
            ca += [(0.0, y, z) for z in np.linspace(-18.0, 18.0, 25)]
            sses.append(fake_sse("H", first, len(ca) - 1, centroid=(0.0, y, 0.0), axis=(0, 0, 1), chain=chain))
    for x in (-150.0, 150.0):
        first = len(ca)
        ca += [(x + t, 0.0, 40.0) for t in np.linspace(-9.0, 9.0, 13)]
        sses.append(fake_sse("H", first, len(ca) - 1, centroid=(x, 0.0, 40.0), axis=(1, 0, 0), chain="A"))
    bb = SimpleNamespace(ca=np.array(ca))
    frame = view_frame(sses, up=mem.normal)
    a = [s.centroid for s in sses if s.chain == "A"][:6]
    b = [s.centroid for s in sses if s.chain == "B"]
    assert _worst_tau(frame, [[a, b]]) == 0 and _sequence_score(frame, [a, b]) == 0
    assert _in_sequence(frame, sses, bb, mem, {"A": ("A", 0), "B": ("A", 0)}) is frame
    turned = Frame(frame.origin, np.array([0.0, 1.0, 0.0]), frame.v, np.cross([0.0, 1.0, 0.0], frame.v))
    assert _sequence_score(turned, [a, b]) == 15 - 10 and _worst_tau(turned, [[a, b]]) == -1  # summed, it won


def test_a_turn_that_lowers_the_pooled_tau_is_not_taken():
    """Ruling R27 (6ZKT, a photosystem of many chains: turned 105° for one chain's tau, -1.0 to -0.6, the tau pooled
    over all its band helices fell): two proteins, A's six band helices in sequence along x, B's three against it.
    A quarter turn lifts the worse chain (B: -1 to 0, a gain of 1 over the 0.1 floor) but pools A's 15 pairs and B's 3
    from +12 to 0: the frame as chosen is kept."""
    from types import SimpleNamespace

    from foldmap.cli import _in_sequence, _sequence_score, _worst_tau
    from foldmap.frame import Frame, view_frame
    from foldmap.membrane import Membrane
    from helpers import fake_sse

    mem = Membrane(np.array([0.0, 0.0, 1.0]), 0.0, 15.0, "estimate")
    ca, sses = [], []
    for chain, xs in (("A", (-5.0, -3.0, -1.0, 1.0, 3.0, 5.0)), ("B", (4.0, 0.0, -4.0))):
        for x in xs:
            first = len(ca)
            ca += [(x, 0.0, z) for z in np.linspace(-18.0, 18.0, 25)]
            sses.append(fake_sse("H", first, len(ca) - 1, centroid=(x, 0.0, 0.0), axis=(0, 0, 1), chain=chain))
    for x in (-150.0, 150.0):
        first = len(ca)
        ca += [(x + t, 30.0, 40.0) for t in np.linspace(-9.0, 9.0, 13)]
        sses.append(fake_sse("H", first, len(ca) - 1, centroid=(x, 30.0, 40.0), axis=(1, 0, 0), chain="C"))
    frame = view_frame(sses, up=mem.normal)
    a = [s.centroid for s in sses if s.chain == "A"]
    b = [s.centroid for s in sses if s.chain == "B"]
    quarter = Frame(frame.origin, np.array([0.0, 1.0, 0.0]), frame.v, np.cross([0.0, 1.0, 0.0], frame.v))
    assert _worst_tau(frame, [[a], [b]]) == -1 and _worst_tau(quarter, [[a], [b]]) == 0  # the worse chain gains 1
    assert _sequence_score(frame, [a, b]) == 12 and _sequence_score(quarter, [a, b]) == 0  # pooled falls
    assert _in_sequence(frame, sses, SimpleNamespace(ca=np.array(ca)), mem) is frame


def test_the_membrane_join_note_lists_only_joins_the_copies_keep():
    """Review-5 M1: a join of one copy that _copy_breaks cuts back (the other copy broken there, the join refused for
    it: a tie, so the whole run is cut) is no join, and build_sses' `matched["membrane"]` leaves it out; joined in
    both copies, both are listed."""
    from foldmap.model import Backbone, ResLabel
    from foldmap.ss import build_sses

    t = np.radians(100.0) * np.arange(25)
    one = np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(25)], axis=1)
    ca = np.vstack([one, one + np.array([30.0, 0.0, 0.0])])
    labels = [ResLabel(c, k + 1, "", "ALA") for c in "AB" for k in range(25)]
    bb = Backbone(labels, np.repeat(ca[:, None, :], 4, axis=1), np.array([k % 25 > 0 for k in range(50)]))
    ss = ("H" * 12 + "T" + "H" * 12) * 2
    only_a: dict = {}
    sses = build_sses(bb, ss, matched=only_a, join=lambda x, y: bb.labels[x[0]].chain == bb.labels[y[0]].chain == "A")
    assert only_a["membrane"] == [] and only_a["cut"] == 1
    assert [s.id for s in sses if s.chain == "A"] == ["A:1-12", "A:14-25"]
    both: dict = {}
    build_sses(bb, ss, matched=both, join=lambda x, y: bb.labels[x[0]].chain == bb.labels[y[0]].chain)
    assert both["membrane"] == ["A:1-25", "B:1-25"] and both["cut"] == 0


def test_make_layout_wires_the_joins_the_flat_rule_and_the_turn(monkeypatch):
    """Review-5 M7, on a tracked fixture (1C3W, bacteriorhodopsin, one chain): a projected membrane view asks
    membrane.joins of its helix runs and copies' breaks (crosses_joined), lays helices flat by lies_flat_here, passes
    `nearer` to build_layout and draws the frame _in_sequence returns; with no membrane none of them is asked, and in
    stack mode neither the joins nor the turn."""
    import foldmap.cli as cli
    from foldmap.frame import Frame

    asked: dict[str, int] = {}

    def spy(name):
        real = getattr(cli, name)

        def wrapped(*a, **k):
            asked[name] = asked.get(name, 0) + 1
            return real(*a, **k)

        monkeypatch.setattr(cli, name, wrapped)

    for name in ("membrane_joins", "lies_flat_here", "_nearer_faces"):
        spy(name)
    monkeypatch.setattr(cli, "crosses_joined", lambda *a: asked.update(crosses_joined=1) or False)
    turned = {}

    def quarter(frame, *a, **k):  # the frame turned a quarter about the normal, so its use shows on the page axes
        u = np.cross(frame.v, frame.u)
        turned["u"] = u
        return Frame(frame.origin, u, frame.v, np.cross(u, frame.v))

    monkeypatch.setattr(cli, "_in_sequence", quarter)
    built = {}
    real_build = cli.build_layout

    def build(*a, **k):
        built.update(k)
        return real_build(*a, **k)

    monkeypatch.setattr(cli, "build_layout", build)
    one = Path(__file__).parent / "data" / "1C3W.cif"
    lay, _, _ = make_layout(one, membrane="auto", chains=["A1"])
    assert asked.get("membrane_joins") and asked.get("lies_flat_here") and asked.get("_nearer_faces")
    assert "nearer" in built and np.allclose(lay.page_axes[0], turned["u"])
    asked.clear(), turned.clear(), built.clear()
    make_layout(one, chains=["A1"])
    assert not asked and not turned and "nearer" not in built
    make_layout(one, "stack", membrane="auto", chains=["A1"])
    assert not asked.get("membrane_joins") and not turned


@pytest.mark.parametrize("mirror", [None, 0.0])
def test_tied_columns_are_dealt_their_projected_order_and_their_copies_follow(mirror):
    """Ruling R29 (6BO8, an unrolled tetramer: its protomer-0 columns were left to _attract, being symmetry mates,
    and came out in reverse order in every chain, tau 0.6 -> -0.867): layout._in_order deals a mate's columns their
    projected order too, and each copy tied to a mate moves as far as its mate (the other way when copies are drawn
    as mirror images), so the copies stay drawn alike. It says whether anything moved."""
    from foldmap.layout import Placed, _in_order, _Item
    from helpers import fake_sse

    items = []
    for k, (chain, x) in enumerate((("A", 5.0), ("A", -5.0), ("A", 1.0), ("B", 25.0), ("B", 15.0), ("B", 21.0))):
        s = fake_sse("H", 10 * k, 10 * k + 20, chain=chain)
        items.append(_Item([Placed(s, 0, 0, 5.0, np.pi / 2, 1.0)], [(0.0, 0.0)], x, 0.0, (x, 0.0), k, chain))
    level = {it.members[0].sse.id for it in items}
    ties = {3: (0, 1), 4: (1, 1), 5: (2, 1)}  # B's columns are A's copies, protomer 1
    across = [-3.0, 0.0, 2.0, 17.0, 20.0, 22.0]  # projected: A's item 0 leftmost, then 1, then 2
    assert _in_order(items, level, across, ties, mirror)
    assert [it.cx for it in items[:3]] == pytest.approx([-5.0, 1.0, 5.0])
    moves = np.array([-10.0, 6.0, 4.0]) * (-1 if mirror is not None else 1)
    assert [it.cx for it in items[3:]] == pytest.approx(list(np.array([25.0, 15.0, 21.0]) + moves))
    assert not _in_order(items, level, across, ties, mirror)  # in order now: nothing moves


def test_a_missing_chain_is_said_before_the_membrane_is_looked_for(monkeypatch):
    """Review-5 M3: make_layout looks for the membrane before DSSP's elements (a membrane view joins helix pieces,
    ruling T5-1), but after the cheap checks: a chain the structure lacks fails without a membrane search."""
    import foldmap.membrane as membrane

    called = []
    monkeypatch.setattr(membrane, "find_membrane", lambda *a: called.append(1))
    with pytest.raises(ValueError, match="no chain Z"):
        make_layout(Path(__file__).parent / "data" / "1C3W.cif", membrane="auto", chains=["Z"])
    assert not called
