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


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "2HHB", "8UUP", "8UTF", "6ZS5", "1TIM", "AF-P04637"])
def test_soluble_proteins_are_not(name):
    assert find_membrane(DATA / f"{name}.cif", load_backbone(DATA / f"{name}.cif", "asu")) is None


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
        y0, y1 = _row_band(lay, p)  # the trimer wraps one chain per row, each row with its band
        lo, hi = sorted((p.n_port[1], p.c_port[1]))
        assert min(hi, y1) - max(lo, y0) >= 0.75 * (y1 - y0), s.id  # measured worst: α3, 0.83


@pytest.mark.parametrize(
    "path, chains",
    [(DATA / "1C3W.cif", ["A1", "A2"]), (OPM, None)],
    ids=["1C3W-two-chains", "1BL8-opm"],
)
def test_a_wide_depth_keeping_stack_wraps_one_chain_per_row_each_with_its_band(path, chains, monkeypatch):
    """1PV7 --mode stack --membrane auto put 21 helices in one row 2160 px wide. A row wider than _MAX_ROW wraps:
    one chain per row, each row with its own band, and each element at its depth: a whole row moves as one."""
    from foldmap import layout as layout_mod

    lay, sses, _ = make_layout(path, mode="stack", membrane="auto", chains=chains)
    names = sorted({s.chain for s in sses})
    assert len(names) >= 2 and len(lay.membrane["bands"]) == len(names)
    assert lay.membrane["y"] == lay.membrane["bands"][0]
    span = {
        c: (
            min(p.rect[1] for p in lay.placed.values() if p.sse.chain == c),
            max(p.rect[3] for p in lay.placed.values() if p.sse.chain == c),
        )
        for c in names
    }
    rows = sorted(span.values(), reverse=True)
    assert all(lo > hi for (lo, _), (_, hi) in zip(rows, rows[1:]))  # one chain per row, rows apart
    for p in lay.placed.values():
        if p.sse.kind == "H" and len(p.sse) >= 20:  # TM helices cross their own row's band
            y0, y1 = _row_band(lay, p)
            lo, hi = sorted((p.n_port[1], p.c_port[1]))
            assert min(hi, y1) - max(lo, y0) >= 0.75 * (y1 - y0), p.sse.id
    fig = make_figure(path, mode="stack", membrane="auto", chains=chains)
    drawn = [a for a in fig.axes[0].patches if (a.get_gid() or "").split(":")[0] == "membrane"]
    assert len(drawn) == len(names)
    monkeypatch.setattr(layout_mod, "_MAX_ROW", 1e9)  # the same figure in one row: every element's depth
    one, _, _ = make_layout(path, mode="stack", membrane="auto", chains=chains)
    for c in names:
        shift = {round(p.cy - one.placed[k].cy, 6) for k, p in lay.placed.items() if p.sse.chain == c}
        assert len(shift) == 1, c


def test_a_narrow_depth_keeping_stack_stays_one_row():
    lay, sses, _ = make_layout(OPM, mode="stack", membrane="auto", chains=["A", "B"])
    assert lay.width <= 30.0 and len(lay.membrane["bands"]) == 1
