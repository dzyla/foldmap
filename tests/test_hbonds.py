"""Hydrogen bonds from the coordinates (foldmap.hbonds): known interactions found, selections honoured, waters
ignored, output deterministic and fast."""

import math
import os
import subprocess
import sys
import time
from functools import lru_cache
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap import hbonds
from foldmap.hbonds import HBond, Selection, find_hbonds, parse_selection
from foldmap.io import read_model
from helpers import text_centre

DATA = Path(__file__).parent / "data"


def _path(name):
    return next(p for p in (DATA / "na" / f"{name}.cif.gz", DATA / f"{name}.cif") if p.is_file())


@lru_cache(maxsize=None)
def _model(name):
    return read_model(_path(name))


@lru_cache(maxsize=None)
def _found(name, select):
    return tuple(find_hbonds(_model(name), select))


def _tag(atom):
    return f"{atom.res.chain}:{atom.res.name}{atom.res.tag} {atom.name}"


def _protein(atom) -> bool:
    info = gemmi.find_tabulated_residue(atom.res.name)
    return info is not None and info.is_amino_acid()


_BACKBONE = {"N", "CA", "C", "O", "OXT"}
_NOT_BASE = {"P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}

# (donor, acceptor, distance in Å as measured in the file)
LAMBDA_SIDE_CHAIN_TO_BASE = {
    ("1:DA4 N6", "3:GLN44 OE1", 3.12),
    ("3:GLN44 NE2", "1:DA4 N7", 3.32),
    ("2:DA24 N6", "4:GLN44 OE1", 2.94),
    ("4:GLN44 NE2", "2:DA24 N7", 2.91),
    ("3:SER45 OG", "2:DG36 O6", 2.60),
    ("3:SER45 OG", "2:DG37 O6", 2.99),
    ("1:DC5 N4", "3:SER45 OG", 3.04),
    ("4:SER45 OG", "1:DG16 O6", 2.88),
    ("4:SER45 OG", "1:DG16 N7", 3.18),
    ("3:ASN55 ND2", "2:DG34 N7", 3.06),
    ("4:ASN55 ND2", "1:DG14 N7", 2.91),
    ("4:LYS4 NZ", "1:DG13 O6", 2.80),
    ("4:LYS4 NZ", "1:DG13 N7", 3.33),
    ("4:LYS4 NZ", "1:DG14 O6", 2.90),
    ("4:LYS4 NZ", "1:DG14 N7", 3.42),
    ("4:LYS3 NZ", "1:DG11 O6", 2.57),
    ("4:LYS3 NZ", "2:DG30 O6", 2.99),
}


def test_lambda_repressor_side_chains_bond_the_major_groove_bases():
    """1LMB, lambda repressor (chains 3 and 4, residues numbered 1-92 as in the file) on its operator (DNA chains
    1 and 2). The file's primary citation is Beamer & Pabo, J. Mol. Biol. 227:177-196 (1992); its citation 1
    (Clarke et al., Science 254:267-270, 1991) is on the N-terminal arm. The interactions usually described for
    this complex, as they appear in this file (each checked against the coordinates, distances measured
    heavy atom to heavy atom):

    - Gln44 to the adenine of base pair 2: Gln44 OE1 accepts from DA N6 and Gln44 NE2 donates to DA N7, in both
      half-sites (3:GLN44 with 1:DA4; 4:GLN44 with 2:DA24).
    - Ser45 OG to guanine O6/N7 (3:SER45 with 2:DG36, 2:DG37 and, as acceptor, 1:DC5 N4; 4:SER45 with 1:DG16).
      4:SER45 OG to 1:DG16 O6 (2.88 Å) has a CB-OG-O6 angle of 86.6°: found because a hydroxyl donor, whose H
      rotates, needs only 60°.
    - Asn55 ND2 to guanine N7 (3:ASN55 with 2:DG34; 4:ASN55 with 1:DG14).
    - Lys4 NZ of the N-terminal arm to guanine O6/N7 (4:LYS4 with 1:DG13 and 1:DG14; the arm of chain 3 starts
      at residue 6 and is not ordered), and Lys3 NZ to guanine O6 (4:LYS3 with 1:DG11 and 2:DG30).
    - Gln33 reaches the phosphate backbone, not a base: its main-chain N and side-chain NE2 donate to OP2
      (3:GLN33 with 1:DT3 and 1:DA4; 4:GLN33 with 2:DT23 and 2:DA24).

    These are every side-chain N/O/S to base N/O pair of the file within 2.5-3.5 Å (each a donor with an
    acceptor), so the side-chain to base bonds found must be exactly this set."""
    found = _found("1LMB", "na")
    side_to_base = set()
    for b in found:
        prot, nt = (b.donor, b.acceptor) if _protein(b.donor) else (b.acceptor, b.donor)
        if prot.name not in _BACKBONE and nt.name not in _NOT_BASE:
            side_to_base.add((_tag(b.donor), _tag(b.acceptor), round(b.distance, 2)))
    assert side_to_base == LAMBDA_SIDE_CHAIN_TO_BASE
    tags = {(_tag(b.donor), _tag(b.acceptor)) for b in found}
    assert {("3:GLN33 N", "1:DT3 OP2"), ("3:GLN33 NE2", "1:DA4 OP2")} <= tags
    assert {("4:GLN33 N", "2:DT23 OP2"), ("4:GLN33 NE2", "2:DA24 OP2")} <= tags


def test_na_bonds_join_protein_and_nucleic_acid():
    found = _found("1LMB", "na")
    assert found and all(_protein(b.donor) != _protein(b.acceptor) for b in found)
    assert all(2.5 <= b.distance <= 3.5 for b in found)


def test_haemoglobin_alpha1_beta1_interface():
    found = _found("2HHB", "chains:A+B")
    assert len(found) > 0
    for b in found:
        assert {b.donor.res.chain, b.acceptor.res.chain} == {"A", "B"}
        assert 2.5 <= b.distance <= 3.5
    assert found == _found("2HHB", "chains:B+A")
    assert set(found) <= set(_found("2HHB", "interchain"))


def test_interchain_covers_every_chain_pair_and_no_chain_with_itself():
    found = _found("2HHB", "interchain")
    pairs = {frozenset((b.donor.res.chain, b.acceptor.res.chain)) for b in found}
    assert all(len(p) == 2 for p in pairs)
    assert {frozenset("AB"), frozenset("CD")} <= pairs


def test_waters_are_ignored():
    model, _ = _model("1LMB")
    waters = [r for c in model for r in c if r.is_water()]
    assert waters  # the file has waters, some within bonding distance of the protein
    found = _found("1LMB", "interchain")
    assert found and all("HOH" not in (b.donor.res.name, b.acceptor.res.name) for b in found)


def test_each_atom_pair_once_and_sorted():
    found = _found("2HHB", "interchain")
    assert list(found) == sorted(found)
    pairs = [frozenset((b.donor, b.acceptor)) for b in found]
    assert len(pairs) == len(set(pairs))


def test_off_finds_nothing():
    assert find_hbonds(_model("1LMB"), "off") == []


_PRINT = """
import sys
from foldmap.hbonds import find_hbonds
from foldmap.io import read_model
print(repr(find_hbonds(read_model(sys.argv[1]), sys.argv[2])))
"""


def test_output_does_not_depend_on_the_hash_seed():
    runs = set()
    for seed in ("0", "12345", "777"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run(
            [sys.executable, "-c", _PRINT, str(_path("1LMB")), "interchain"],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1 and "HBond(" in runs.pop()


@pytest.mark.parametrize("select", ["na", "interchain"])
def test_ribosome_fixture_is_fast(select):
    model = _model("4YBB")
    start = time.perf_counter()
    found = find_hbonds(model, select)
    assert time.perf_counter() - start < 2.0
    assert found


@pytest.mark.parametrize(
    "text, canonical",
    [
        ("off", "off"),
        ("na", "na"),
        (" interchain ", "interchain"),
        ("chains:A+B", "chains:A+B"),
        ("chains: A+B , C+D", "chains:A+B,C+D"),
        ("chains:B+A", "chains:A+B"),
        ("chains:A+B,B+A,C+D", "chains:A+B,C+D"),
        ("NA", "na"),
        ("Interchain", "interchain"),
        ("CHAINS:a+B", "chains:B+a"),
    ],
)
def test_selection_parses_and_prints_back(text, canonical):
    sel = parse_selection(text)
    assert isinstance(sel, Selection) and str(sel) == canonical
    assert parse_selection(str(sel)) == sel


@pytest.mark.parametrize(
    "text, words",
    [
        ("ligands", "reserved"),
        ("all", "off, na, interchain or chains:"),
        ("chains:", "chain pair"),
        ("chains:A", "chain pair"),
        ("chains:A+B,", "chain pair"),
        ("chains:A+B+C", "chain pair"),
        ("chains:A+A", "two different chains"),
    ],
)
def test_bad_selections_are_refused_clearly(text, words):
    with pytest.raises(ValueError, match="hbonds") as err:
        parse_selection(text)
    assert words in str(err.value)


def test_a_chain_the_model_lacks_is_named():
    with pytest.raises(ValueError, match="Z"):
        find_hbonds(_model("2HHB"), "chains:A+Z")


def test_hbond_carries_donor_acceptor_distance():
    b = _found("1LMB", "na")[0]
    assert isinstance(b, HBond) and b == HBond(b.donor, b.acceptor, b.distance)


def test_selection_equality_ignores_pair_order_and_kind_is_checked():
    assert parse_selection("chains:A+B") == parse_selection("chains:B+A") == Selection("chains", (("B", "A"),))
    with pytest.raises(ValueError, match="hbonds"):
        Selection("bogus")
    with pytest.raises(ValueError, match="hbonds"):
        Selection("chains")
    with pytest.raises(ValueError, match="hbonds"):
        Selection("na", (("A", "B"),))


def _roles(model, chain, seq):
    """atom name -> "D", "A" or "DA" for the donor and acceptor atoms of one residue."""
    sites = hbonds._sites(model[0] if isinstance(model, tuple) else model, None)
    out = {}
    for k, name in enumerate(sites.names):
        lab = sites.labels[sites.owner[k]]
        if (lab.chain, lab.seq) == (chain, seq):
            out[name] = ("D" if sites.donor[k] else "") + ("A" if sites.acceptor[k] else "")
    return out


def test_modified_nucleotides_keep_only_the_roles_their_chemistry_allows():
    """1EHZ, tRNA-Phe. A methyl on an amino N leaves one H (2MG N2 donates), two leave none (M2G N2); a methyl on a
    ring N takes its lone pair (7MG N7, 1MA N1 accept nothing); 2'-O-methyl leaves O2' an acceptor only (OMC);
    wybutosine (YYG) has N1, N2 and N3 all substituted: N1 (a bridgehead) and N3 (methylated) neither donate nor
    accept, and N2, closed into the extra imidazole ring as an imine, only accepts."""
    model = _model("1EHZ")
    assert _roles(model, "A", 10)["N2"] == "D"  # 2MG
    assert "N2" not in _roles(model, "A", 26)  # M2G
    assert "N7" not in _roles(model, "A", 46) and _roles(model, "A", 46)["N2"] == "D"  # 7MG
    assert "N1" not in _roles(model, "A", 58) and _roles(model, "A", 58)["N6"] == "D"  # 1MA
    assert _roles(model, "A", 32)["O2'"] == "A"  # OMC
    assert _roles(model, "A", 31)["O2'"] == "DA"  # an unmodified neighbour
    yyg = _roles(model, "A", 37)
    assert {n: r for n, r in yyg.items() if n in ("N1", "N2", "N3", "N7", "O6")} == {"N2": "A", "N7": "A", "O6": "A"}


def test_pseudouridine_n1_donates_and_uridine_n1_does_not():
    model = _model("4YBB-psi")
    sites = hbonds._sites(model[0], None)
    by_name = {}
    for k, name in enumerate(sites.names):
        by_name.setdefault(sites.labels[sites.owner[k]].name, {})[name] = bool(sites.donor[k])
    assert by_name["PSU"]["N1"] and by_name["PSU"]["N3"]
    assert "N1" not in by_name["U"] and by_name["U"]["N3"]


def test_terminal_hydroxyls_donate_and_linked_ones_only_accept():
    model = _model("1LMB")
    first, last = _roles(model, "1", 1), _roles(model, "1", 20)
    middle = _roles(model, "1", 10)
    assert first["O5'"] == "DA" and last["O3'"] == "DA"
    assert middle["O5'"] == "A" and middle["O3'"] == "A"


def _model_of(*residues):
    """A model from (chain, seq, name, {atom: (x, y, z)}), each a polymer residue, in that order."""
    model = gemmi.Model("1")
    for chain, seq, name, atoms in residues:
        if chain not in [c.name for c in model]:
            model.add_chain(gemmi.Chain(chain))
        res = gemmi.Residue()
        res.name, res.seqid, res.entity_type = name, gemmi.SeqId(seq, " "), gemmi.EntityType.Polymer
        for atom_name, xyz in atoms.items():
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = atom_name, gemmi.Element(atom_name[0]), gemmi.Position(*xyz)
            res.add_atom(atom)
        model[chain].add_residue(res)
    return model


def _at(angle, dist=3.0):
    """A point `dist` from the origin at `angle` degrees from the -x axis (where antecedents sit)."""
    t = math.radians(angle)
    return (-dist * math.cos(t), dist * math.sin(t), 0.0)


@pytest.mark.parametrize("angle, found", [(100, True), (75, True), (50, False)])
def test_a_hydroxyl_donor_needs_60_degrees(angle, found):
    model = _model_of(("A", 1, "SER", {"CB": (-1.43, 0, 0), "OG": (0, 0, 0)}), ("B", 1, "ASP", {"OD1": _at(angle)}))
    assert bool(find_hbonds(model, "interchain")) == found


@pytest.mark.parametrize("angle, found", [(100, True), (75, False)])
def test_an_sp2_donor_needs_90_degrees(angle, found):
    model = _model_of(
        ("A", 1, "TRP", {"CD1": (-0.7, 1.1, 0), "CE2": (-0.7, -1.1, 0), "NE1": (0, 0, 0)}),
        ("B", 1, "ASP", {"OD1": _at(angle)}),
    )
    assert bool(find_hbonds(model, "interchain")) == found


def test_an_unknown_antecedent_passes():
    cb = {"CB": (-1.43, 0, 0)}
    acceptor = ("B", 1, "ASP", {"OD1": _at(30)})
    assert find_hbonds(_model_of(("A", 1, "SER", {**cb, "OG": (0, 0, 0)}), acceptor), "interchain") == []
    found = find_hbonds(_model_of(("A", 1, "SER", {"OG": (0, 0, 0)}), acceptor), "interchain")
    assert [(str(b.donor), str(b.acceptor)) for b in found] == [("A:SER1 OG", "B:ASP1 OD1")]


def test_a_pair_where_both_atoms_can_donate_is_counted_once():
    model = _model_of(
        ("A", 1, "SER", {"CB": (-1.43, 0, 0), "OG": (0, 0, 0)}),
        ("B", 1, "THR", {"CB": (2.8 + 1.43, 0, 0), "OG1": (2.8, 0, 0)}),
    )
    found = find_hbonds(model, "interchain")
    assert [(str(b.donor), str(b.acceptor), round(b.distance, 2)) for b in found] == [("A:SER1 OG", "B:THR1 OG1", 2.8)]


def test_microheterogeneity_keeps_the_first_residue_of_a_number():
    ser = ("A", 1, "SER", {"CB": (-1.43, 0, 0), "OG": (0, 0, 0)})
    acceptor = ("B", 1, "ASP", {"OD1": (10.0, 0, 0)})
    thr = {"CB": (10 - 2.9 - 1.43, 0, 0), "OG1": (10 - 2.9, 0, 0)}
    assert find_hbonds(_model_of(ser, ("A", 2, "THR", thr), acceptor), "interchain")  # its own number: found
    assert find_hbonds(_model_of(ser, ("A", 1, "THR", thr), acceptor), "interchain") == []  # second at 1: dropped


def test_a_polymer_residue_type_bound_as_a_ligand_is_left_out():
    model = _model("1LMB")[0].clone()
    for res in model["4"]:
        if res.seqid.num == 44:
            res.entity_type = gemmi.EntityType.NonPolymer
    tags = {str(a) for b in find_hbonds(model, "na") for a in (b.donor, b.acceptor)}
    assert tags and not any(t.startswith("4:GLN44 ") for t in tags)
    assert any(str(a).startswith("4:GLN44 ") for b in _found("1LMB", "na") for a in (b.donor, b.acceptor))


# ---- the overlay: dashed lines between the residues' points on the page --------------------------------------------

from foldmap.cli import make_layout  # noqa: E402
from foldmap.nacheck import deposited_hydrog  # noqa: E402
from foldmap.render import (  # noqa: E402
    INTERACTION_ONE,
    draw,
    interaction_overlay,
    nucleotide_points,
    residue_name,
    residue_point,
)
from foldmap.route import route_loops  # noqa: E402
from foldmap.style import Style, resolve_style  # noqa: E402


def _hbond_overlay(layout, loops, sses, look, model):
    """The H-bond lines the look asks for: interaction_overlay's H-bond part."""
    return interaction_overlay(layout, loops, sses, look, model, kinds=("hbond",)).hbonds


@lru_cache(maxsize=None)
def _drawn(name, *settings):
    """(figure, layout, loops, sses, bb, look) of the figure drawn with these style settings."""
    look = resolve_style(None, None, list(settings))
    lay, sses, bb = make_layout(_path(name), look=look, model=_model(name))
    loops = route_loops(lay, sses, bb)
    return draw(lay, loops, sses, look=look, model=_model(name)), lay, loops, sses, bb, look


def _lines(fig):
    return [l for l in fig.axes[0].lines if (l.get_gid() or "").startswith("hbond:")]


def _legend_texts(fig):
    return [t.get_text() for t in fig.axes[0].texts if t.get_text().startswith("H-bonds")]


def test_na_overlay_draws_a_line_per_bonded_pair_between_the_residue_and_the_nucleotide():
    fig, lay, loops, sses, bb, look = _drawn("1LMB", "hbonds=na", "hbond_max=1000")
    lines = {l.get_gid(): l for l in _lines(fig)}
    index = {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}
    where = nucleotide_points(lay)
    nt = {
        (lab.chain, lab.seq, lab.icode): where[(s, i)]
        for s, strand in enumerate(lay.nucleic.strands)
        for i, lab in enumerate(strand.labels)
        if (s, i) in where
    }
    pairs = {}
    for b in find_hbonds(_model("1LMB"), "na"):
        res, base = (b.donor.res, b.acceptor.res) if _protein(b.donor) else (b.acceptor.res, b.donor.res)
        pairs.setdefault((res, base), []).append(b)
    assert len(pairs) > 20
    drawn = 0
    for res, base in pairs:
        gid = f"hbond:{res.chain}:{res.tag}-{base.chain}:{base.tag}"
        point = residue_point(lay, loops, sses, lay.res_chain, index[(res.chain, res.seq, res.icode)])
        if point is None or (base.chain, base.seq, base.icode) not in nt:
            assert gid not in lines
            continue
        xs, ys = lines[gid].get_data()
        assert (xs[0], ys[0]) == pytest.approx(point)
        assert (xs[-1], ys[-1]) == pytest.approx(tuple(nt[(base.chain, base.seq, base.icode)]))
        drawn += 1
    assert drawn == len(lines) >= 20
    assert all(l.get_linestyle() != "-" for l in lines.values())  # dashed


def test_the_legend_counts_the_lines_and_says_when_the_cap_cut_some():
    fig = _drawn("1LMB", "hbonds=na", "hbond_max=1000")[0]
    n = len(_lines(fig))
    assert _legend_texts(fig) == [f"H-bonds ({n})"]
    assert any(l.get_gid() == "legend-hbond" for l in fig.axes[0].lines)
    capped = _drawn("1LMB", "hbonds=na", "hbond_max=10")[0]
    assert len(_lines(capped)) == 10
    assert _legend_texts(capped) == [f"H-bonds (10 of {n})"]


def test_hbond_max_keeps_the_shortest_lines():
    all_ = _hbond_overlay(
        *_drawn("1LMB", "hbonds=na", "hbond_max=1000")[1:4], Style(hbonds="na", hbond_max=1000), _model("1LMB")
    )
    few = _hbond_overlay(
        *_drawn("1LMB", "hbonds=na", "hbond_max=1000")[1:4], Style(hbonds="na", hbond_max=7), _model("1LMB")
    )
    assert [l.gid for l in few.lines] == [l.gid for l in all_.lines[:7]]
    lengths = [l.length for l in all_.lines]
    assert lengths == sorted(lengths) and few.total == len(all_.lines)


def test_hbonds_without_the_structure_is_an_error_not_an_empty_figure():
    _, lay, loops, sses, _, look = _drawn("1LMB", "hbonds=na")
    with pytest.raises(ValueError, match="model="):
        draw(lay, loops, sses, look=look)
    draw(lay, loops, sses, look=Style(hbonds="off"))  # no model needed when nothing is asked
    assert not hasattr(lay, "model")


def test_the_legend_counts_pairs_that_are_not_drawn():
    from foldmap.render import HBondOverlay

    assert HBondOverlay([1] * 5, 5).legend() == "H-bonds (5)"
    assert HBondOverlay([1] * 5, 9).legend() == "H-bonds (5 of 9)"
    assert HBondOverlay([1] * 5, 5, 3).legend() == "H-bonds (5, 3 not shown)"
    assert HBondOverlay([1] * 5, 9, 2).legend() == "H-bonds (5 of 9, 2 not shown)"
    look = resolve_style(None, None, ["hbonds=interchain", "hbond_max=1000"])  # only chains A and B are drawn
    lay, sses, bb = make_layout(_path("2HHB"), look=look, model=_model("2HHB"), chains=["A", "B"])
    overlay = _hbond_overlay(lay, route_loops(lay, sses, bb), sses, look, _model("2HHB"))
    assert overlay.unplaced > 0 and overlay.legend() == f"H-bonds ({len(overlay.lines)}, {overlay.unplaced} not shown)"


def test_long_lines_fade_and_short_ones_do_not():
    from foldmap.render import hbond_alpha

    assert hbond_alpha(0.0, 100.0) == 1.0 and hbond_alpha(15.0, 100.0) == 1.0
    assert hbond_alpha(50.0, 100.0) == pytest.approx(0.4) == hbond_alpha(400.0, 100.0)
    assert 0.4 < hbond_alpha(30.0, 100.0) < 1.0
    fig = _drawn("30QM", "hbonds=na")[0]
    alphas = [l.get_alpha() for l in _lines(fig)]
    assert max(alphas) == 1.0 and min(alphas) < 0.9 and min(alphas) >= 0.4


def test_lines_sit_below_the_beads_and_the_rails():
    fig = _drawn("1LMB", "hbonds=na")[0]
    assert {l.get_zorder() for l in _lines(fig)} == {2.15}
    beads = [a for a in fig.axes[0].collections + fig.axes[0].patches if a.get_zorder() >= 2.2]
    assert beads  # the beads and rails draw above the dashes


def test_off_draws_nothing_and_the_overlay_moves_no_element():
    off, lay_off = _drawn("1LMB", "hbonds=off")[:2]
    on, lay_on = _drawn("1LMB", "hbonds=na")[:2]
    assert not _lines(off) and not _legend_texts(off) and _lines(on)
    assert lay_off.bounds == lay_on.bounds
    assert {k: (p.n_port, p.c_port) for k, p in lay_off.placed.items()} == {
        k: (p.n_port, p.c_port) for k, p in lay_on.placed.items()
    }
    assert [(d.axial, d.point(next(iter(d.axial))).tolist()) for d in lay_off.stems] == [
        (d.axial, d.point(next(iter(d.axial))).tolist()) for d in lay_on.stems
    ]


def test_chain_pairs_draw_only_lines_between_those_chains():
    fig, lay, *_ = _drawn("2HHB", "hbonds=chains:A+B")
    lines = _lines(fig)
    assert len(lines) > 3
    for l in lines:
        a, b = l.get_gid()[len("hbond:") :].split("-")
        assert {a.split(":")[0], b.split(":")[0]} == {"A", "B"}


def test_the_interactive_page_lists_the_bonds_of_each_line():
    from foldmap.interactive import _interactions

    fig = _drawn("1LMB", "hbonds=na", "hbond_max=1000")[0]
    notes = _interactions(fig)
    assert len(notes) == len(_lines(fig)) and {n["key"] for n in notes} == {l.get_gid() for l in _lines(fig)}
    by_key = {n["key"]: n for n in notes}
    gln = next(n for k, n in by_key.items() if k.startswith("hbond:3:44-1:4"))
    assert gln["text"].startswith("H-bond 3:Gln44 – 1:DA4: ") and "Å" in gln["text"] and len(gln["residues"]) == 1
    assert gln["het"] == [["1", 4]] and gln["label"] == "H-bond" and gln["kind"] == "hbond"


def test_the_interactive_page_lists_every_interaction_line_with_its_atoms():
    """Plan 13c-5: each line drawn (hbond: and ia: gids) is in the page's data under its gid, with its residues, its
    nucleotides and a note naming its atoms (in the order of its residues), distances and, for stacking and
    cation-π, the angle; a protein-protein line holds both residues."""
    from foldmap.interactive import _interactions

    fig = _drawn("1LMB", "interactions=all,hydrophobic", "interaction_max=1000", "hbond_max=1000")[0]
    notes = {n["key"]: n for n in _interactions(fig)}
    gids = {l.get_gid() for l in fig.axes[0].lines if (l.get_gid() or "").startswith(("hbond:", "ia:"))}
    assert notes and set(notes) == gids and {n["kind"] for n in notes.values()} >= {"hbond", "salt", "cation-pi"}
    for line in fig.foldmap_interactions.lines():
        note = notes[line.gid]
        first, second = (residue_name(lab) for lab in line.labels)
        assert note["text"].startswith(f"{INTERACTION_ONE[note['kind']]} {first} – {second}: ")
        assert note["residues"] == list(line.residues) and len(note["residues"]) + len(note["het"]) == 2
    salt = next(n for n in notes.values() if n["kind"] == "salt" and n["het"])
    line = next(l for l in fig.foldmap_interactions.lines() if l.gid == salt["key"])
    x = line.found[0]
    ours = x.a[2] if (x.a[0], x.a[1]) == (line.labels[0].chain, line.labels[0].tag) else x.b[2]
    assert salt["text"].split(": ", 1)[1].startswith(f"{ours} – ") and f"{x.distance:.1f} Å" in salt["text"]
    cation = next(n for n in notes.values() if n["kind"] == "cation-pi")
    assert "-ring" in cation["text"] and "°" in cation["text"]


def test_bad_hbond_keys_are_explained():
    for text, match in (
        ("hbonds=sometimes", "hbonds"),
        ("hbonds=chains:A", "chain pair"),
        ("hbond_max=0", "hbond_max"),
        ("hbond_max=lots", "hbond_max"),
        ("hbond_max=2.5", "hbond_max"),
        ("hbond_color=notacolour", "hbond_color"),
    ):
        with pytest.raises(ValueError, match=match):
            resolve_style(None, None, [text])
    assert resolve_style(None, None, ["hbonds=NA", "hbond_max=50"]).hbond_max == 50


def test_hbond_keys_travel_through_style_layout_files_and_the_app():
    from foldmap.spec import figure_spec, layout_document

    look = resolve_style(None, None, ["hbonds=chains:A+B", "hbond_color=#112233", "hbond_max=40"])
    assert (look.hbonds, look.hbond_color, look.hbond_max) == ("chains:A+B", "#112233", 40)
    doc = layout_document(figure_spec(_path("1LMB"), style=["hbonds=na", "hbond_max=40"]))
    assert doc["style"]["hbonds"] == "na" and doc["style"]["hbond_max"] == 40
    from foldmap.app import look_of

    assert look_of({"hbonds": "interchain", "hbond_max": 12}).hbond_max == 12


def test_the_overlay_agrees_with_the_deposited_hydrogen_bonds_where_the_file_has_them(capsys):
    """Reported, not a gate. Of the file's `hydrog` records: those between protein and nucleic acid, how many the
    overlay draws as a residue pair; and, as a check of the detector, those between two chains (base pairs), how many
    atom pairs `interchain` finds. Today's files annotate base pairs only, so the first count is 0 of 0."""
    rows = []
    for name in ("1LMB", "30QM", "1ZAA", "1KX5", "1IHF", "3CRX"):
        records = deposited_hydrog(_path(name))
        found = {
            frozenset((_tag(b.donor).split()[0], _tag(b.acceptor).split()[0])) for b in _found(name, Selection("na"))
        }
        protein_na = {frozenset((r.a.split()[0], r.b.split()[0])) for r in records if _is_protein_tag(r)}
        atoms = {frozenset((_tag(b.donor), _tag(b.acceptor))) for b in _found(name, Selection("interchain"))}
        across = {frozenset((r.a, r.b)) for r in records if r.a.split(":")[0] != r.b.split(":")[0]}
        rows.append(
            f"{name}: protein-NA {len(protein_na & found)} of {len(protein_na)} drawn; "
            f"base-pair atoms {len(across & atoms)} of {len(across)} found"
        )
    with capsys.disabled():
        print("\ndeposited hydrog vs overlay:\n  " + "\n  ".join(rows))


def _is_protein_tag(record) -> bool:
    names = [r.split()[0].split(":")[1].rstrip("-0123456789") for r in (record.a, record.b)]
    infos = [gemmi.find_tabulated_residue(n) for n in names]
    return any(i is not None and i.is_amino_acid() for i in infos)


def test_the_default_cap_holds_on_a_nucleosome():
    """1KX5 `interchain` has ~480 bonded residue pairs: 200 drawn, the legend says of how many."""
    fig = _drawn("1KX5", "hbonds=interchain")[0]
    assert len(_lines(fig)) == 200
    assert _legend_texts(fig)[0].startswith("H-bonds (200 of ")


@pytest.mark.parametrize("pdb", ["30QM", "1LMB"])
@pytest.mark.parametrize("style", ["ladder", "ribbon", "block", "letters"])
def test_the_nucleotide_end_of_every_line_is_a_drawn_bead_or_letter_in_every_style(pdb, style):
    fig = _drawn(pdb, "hbonds=na", "hbond_max=1000", f"na_style={style}")[0]
    ax = fig.axes[0]
    marks = [xy for c in ax.collections if (c.get_gid() or "").startswith("na:contact:") for xy in c.get_offsets()]
    marks += [text_centre(t) for t in ax.texts if (t.get_gid() or "").startswith("na:nt:")]
    marks = np.array(marks, float)
    lines = _lines(fig)
    assert lines and len(marks)
    for line in lines:
        ends = np.column_stack([line.get_xdata(), line.get_ydata()])
        near = [np.hypot(*(marks - e).T).min() for e in ends]
        assert min(near) < 1e-6, (line.get_gid(), near)


# ---- plan 13c-4: the interaction overlay (salt bridges, stacking, cation-π, hydrophobic) beside the H-bonds ----------

_TODAY = {  # sha256 of `foldmap plot X --set hbonds=na` drawn at 96aebc8, before the overlay generalised
    "matplotlib": "3.11.2",
    "1LMB": "10043b9df29b5c6a5dd6a0dd5efcc86fcea6387988590544a0cb93a1932565a4",
    "1ZAA": "15b6ae433057122da3fe67654206a4e338a6d1df916c38a00eea507c710e35c8",
}


@pytest.mark.parametrize("name", ["1LMB", "1ZAA"])
def test_hbond_figures_are_byte_identical_to_those_before_the_overlay(name, tmp_path):
    """The H-bond lines are interaction_overlay's H-bond part: an hbonds=na figure is the same file as before."""
    import hashlib

    import matplotlib

    from foldmap.cli import main

    if matplotlib.__version__ != _TODAY["matplotlib"]:
        pytest.skip(f"recorded with matplotlib {_TODAY['matplotlib']}")
    out = tmp_path / f"{name}.svg"
    assert main(["plot", str(_path(name)), "--set", "hbonds=na", "-o", str(out)]) == 0
    assert hashlib.sha256(out.read_bytes()).hexdigest() == _TODAY[name]


def _ia(fig, kind=""):
    return [l for l in fig.axes[0].lines if (l.get_gid() or "").startswith(f"ia:{kind}")]


def _points(line):
    """The drawn points of a line (page units), its breaks (NaN) left out."""
    xy = np.column_stack([line.get_xdata(), line.get_ydata()]).astype(float)
    return xy[np.isfinite(xy).all(axis=1)]


def _across(points, ends):
    """Each point's signed distance from the line through `ends`, and its fraction along it."""
    p, q = np.asarray(ends, float)
    d = (q - p) / np.hypot(*(q - p))
    w = np.asarray(points, float) - p
    return w[:, 1] * d[0] - w[:, 0] * d[1], (w @ d) / np.hypot(*(q - p))


def test_salt_bridge_lines_end_on_their_residues():
    """interactions=salt on λ repressor-DNA: ia:salt: lines between their residues, each protein end on its
    residue's place (S4's rule, figcheck.residue_places and science.mark_findings) and each nucleotide end on its
    nucleotide's point; what is drawn lies on that line, between its ends (broken only where it passes a letter)."""
    from types import SimpleNamespace

    from foldmap.figcheck import residue_places
    from foldmap.science import mark_findings

    fig, lay, loops, sses, bb, look = _drawn("1LMB", "interactions=salt")
    lines = {l.get_gid(): l for l in _ia(fig, "salt:")}
    assert len(lines) >= 4 and not _lines(fig) and not _ia(fig, "stack:")
    drawn = SimpleNamespace(layout=lay, loops=loops, sses=sses)
    index = {(l.chain, l.tag): k for k, l in enumerate(lay.res_labels)}
    where = nucleotide_points(lay)
    nts = {
        (lab.chain, lab.tag): tuple(where[(s, i)])
        for s, strand in enumerate(lay.nucleic.strands)
        for i, lab in enumerate(strand.labels)
        if (s, i) in where
    }
    overlay = fig.foldmap_interactions
    assert set(lines) == {l.gid for l in overlay.kinds[0].lines}
    marks = []
    for line in overlay.kinds[0].lines:
        for end, part in zip(line.ends, line.gid[len("ia:salt:") :].split("-")):
            chain, tag = part.split(":")
            if (chain, tag) in index:
                marks.append((part, end, residue_places(drawn, index[(chain, tag)])))
            else:
                assert end == pytest.approx(nts[(chain, tag)])
        off, along = _across(_points(lines[line.gid]), line.ends)
        assert np.abs(off).max() < 1e-9 and along.min() > -1e-9 and along.max() < 1 + 1e-9
    assert marks and mark_findings(marks) == []
    lys = next(l for l in overlay.kinds[0].lines if any(x.a[2] == "NZ" for x in l.found))
    assert all(x.kind == "salt" for x in lys.found) and [x.distance for x in lys.found] == sorted(
        x.distance for x in lys.found
    )


def test_two_kinds_between_one_pair_lie_side_by_side():
    """Two kinds between one residue pair: parallel lines a gap apart (INTERACTION_GAP line widths), one each side
    of the line between the residues' points."""
    from foldmap.render import INTERACTION_GAP

    fig = _drawn("1LMB", "interactions=hbond,salt", "hbond_max=1000")[0]
    hb = {l.get_gid()[len("hbond:") :]: l for l in _lines(fig)}
    salt = {l.get_gid()[len("ia:salt:") :]: l for l in _ia(fig, "salt:")}
    shared = sorted(hb.keys() & salt.keys())
    assert shared
    overlay = fig.foldmap_interactions
    by_gid = {l.gid: l for l in overlay.lines()}
    lw = hb[shared[0]].get_linewidth()  # points
    ax = fig.axes[0]
    pt_per_unit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0] * 72 / fig.dpi
    half = INTERACTION_GAP * lw / pt_per_unit / 2
    for pair in shared:
        centre = by_gid[f"hbond:{pair}"].ends
        assert np.allclose(by_gid[f"ia:salt:{pair}"].ends, centre)
        first, _ = _across(_points(hb[pair]), centre)
        second, _ = _across(_points(salt[pair]), centre)
        assert first == pytest.approx(np.full(len(first), -half))  # parallel, either side, a gap apart
        assert second == pytest.approx(np.full(len(second), half))
        assert overlay.slots[f"hbond:{pair}"] == (0, 2) and overlay.slots[f"ia:salt:{pair}"] == (1, 2)
    alone = [g for g in hb if g not in salt]
    assert alone and all(f"hbond:{g}" not in overlay.slots for g in alone)


def test_each_kind_has_a_legend_line_and_the_cap_says_how_many():
    from foldmap.render import INTERACTION_INK, KindOverlay

    assert KindOverlay("salt", [1] * 5, 5).legend() == "Salt bridges (5)"
    assert KindOverlay("salt", [1] * 200, 1412).legend() == "Salt bridges (200 of 1,412 shown)"
    assert KindOverlay("stack", [1] * 3, 3, 2).legend() == "π-stacking (3, 2 with an end not drawn)"
    every = _drawn("1LMB", "interactions=salt", "interaction_max=1000")[0]
    n = len(_ia(every, "salt:"))
    capped = _drawn("1LMB", "interactions=salt", "interaction_max=3")[0]
    assert len(_ia(capped, "salt:")) == 3
    texts = [t.get_text() for t in capped.axes[0].texts]
    assert f"Salt bridges (3 of {n} shown)" in texts and not any(t.startswith("H-bonds") for t in texts)
    shortest = sorted(every.foldmap_interactions.kinds[0].lines, key=lambda l: (l.length, l.gid))[:3]
    assert [l.gid for l in capped.foldmap_interactions.kinds[0].lines] == [l.gid for l in shortest]
    fig = _drawn("1LMB", "interactions=all,hydrophobic")[0]
    marks = {l.get_gid(): l for l in fig.axes[0].lines if (l.get_gid() or "").startswith("legend-")}
    drawn = {l.get_gid().split(":")[1] for l in _ia(fig)}
    assert drawn and {f"legend-ia:{k}" for k in drawn} <= set(marks)
    for kind in drawn:
        line = next(l for l in _ia(fig, f"{kind}:"))
        assert marks[f"legend-ia:{kind}"].get_color() == line.get_color() == INTERACTION_INK[kind][0]
        assert marks[f"legend-ia:{kind}"].get_linestyle() == line.get_linestyle() != "-"


def test_every_line_drawn_on_lambda_passes_the_independent_re_check():
    from foldmap.interactions import Interaction, verify_interactions

    fig = _drawn("1LMB", "interactions=all", "hbond_max=1000", "interaction_max=1000")[0]
    overlay = fig.foldmap_interactions
    assert overlay.dropped == () and {p.kind for p in overlay.kinds if p.lines} >= {"salt", "cation-pi"}
    found = [x for part in overlay.kinds for line in part.lines for x in line.found]
    found += [
        Interaction(
            "hbond",
            (b.donor.res.chain, b.donor.res.tag, b.donor.name),
            (b.acceptor.res.chain, b.acceptor.res.tag, b.acceptor.name),
            b.distance,
            None,
        )
        for line in overlay.hbonds.lines
        for b in line.bonds
    ]
    assert len(found) > 30 and verify_interactions(_model("1LMB"), found) == []
    assert not any(i["kind"] == "interactions_dropped" for i in fig.foldmap_provenance)


def test_a_line_the_re_check_flags_is_not_drawn_and_the_notes_say_so(monkeypatch):
    from foldmap import render
    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    real = render.failing_lines

    def flag_first_salt(model, drawn):
        n = next(k for k, x in enumerate(drawn) if x.kind == "salt")
        return [*real(model, drawn), (n, f"salt {drawn[n].a} - {drawn[n].b}: moved")]

    look = resolve_style(None, None, ["interactions=salt", "interaction_max=1000"])
    lay, sses, bb = make_layout(_path("1LMB"), look=look, model=_model("1LMB"))
    loops = route_loops(lay, sses, bb)
    before = render.interaction_overlay(lay, loops, sses, look, _model("1LMB"))
    monkeypatch.setattr(render, "failing_lines", flag_first_salt)
    fig = draw(lay, loops, sses, look=look, model=_model("1LMB"))
    after = fig.foldmap_interactions
    assert len(after.dropped) == 1
    assert sum(len(l.found) for l in after.kinds[0].lines) == sum(len(l.found) for l in before.kinds[0].lines) - 1
    item = next(i for i in fig.foldmap_provenance if i["kind"] == "interactions_dropped")
    assert item["count"] == 1 and "re-check" in item["text"] and "moved" in item["text"]


@pytest.mark.parametrize("pdb", ["1LMB", "1URN"])
def test_interaction_lines_keep_off_every_letter_and_a6_reads_them(pdb):
    """The lines of the new kinds are drawn after every text and broken where they would pass a letter, so A6 finds
    no text crossed by an ia: line (the H-bond lines are drawn as before plan 13c); figcheck reads ia: lines as lines
    (LINE_GIDS): a text put across one afterwards is a collision."""
    from types import SimpleNamespace

    from foldmap import figcheck

    fig, lay = _drawn(pdb, "interactions=salt,stack,cation-pi,hydrophobic", "interaction_max=199")[:2]
    drawn = SimpleNamespace(figure=fig, layout=lay)
    page = figcheck.measure(drawn)
    assert any(gid.startswith("ia:") for gid, *_ in page.lines)
    assert not [h for h in figcheck.text_collisions(drawn, page=page) if " × ia:" in h]
    line = max(_ia(fig), key=lambda l: len(_points(l)))
    at = _points(line)[:2].mean(axis=0)  # on its first drawn piece
    fig.axes[0].text(*at, "WWWW", ha="center", va="center", fontsize=9, zorder=9).set_gid("label:test")
    hits = figcheck.text_collisions(drawn)
    assert f"label:test × {line.get_gid()}" in hits


@pytest.mark.parametrize(
    "where, settings, focus",
    [
        ("6LTO.cif.gz", ["cyclic=rosette"], "none"),  # C7, each copy protomer 0 turned by k * 360 / 7 degrees
        ("1BL8.cif", ["cyclic=rosette"], None),  # C4
        ("8UTF.cif", [], None),  # a ring unrolled, a strand lent across its cut
        ("assemblies/1UAM.cif.gz", [], None),  # a C2 drawn as one copy and its mirror image (two_fold=auto)
    ],
    ids=["6LTO-rosette", "1BL8-rosette", "8UTF", "1UAM-mirror"],
)
def test_interaction_lines_follow_turned_and_mirrored_copies(where, settings, focus):
    """Review Focus 4 (plan 13c-5): on a page whose copies are turned (cyclic=rosette) or mirrored (two_fold), each
    interaction line between copies ends on its residues as drawn in their copies (S4's rule, figcheck.residue_places
    worked out apart from render.residue_point), its ink on the segment between them, and no line passes under the
    rosette's axis symbol. A6: the new kinds cross no text; with H-bonds as well (interactions=all), the only texts
    crossed are crossed by hbond: lines, drawn whole as before plan 13c (the parked lambda-hbonds problem, 13c-6)."""
    from matplotlib.path import Path as MPath

    from foldmap import figcheck
    from foldmap.science import mark_findings
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / where
    new = draw_spec(figure_spec(path, style=[*settings, "interactions=salt,stack,cation-pi,hydrophobic"], focus=focus))
    lay = new.layout
    assert lay.rosette in (4, 7) if "cyclic=rosette" in settings else lay.rosette == 0
    if "1UAM" in where:
        assert lay.frame.mirror_x is not None
    assert figcheck.text_collisions(new) == []
    every = draw_spec(figure_spec(path, style=[*settings, "interactions=all,hydrophobic"], focus=focus))
    assert all(" × hbond:" in h for h in figcheck.text_collisions(every))
    for d in (new, every):
        overlay = d.figure.foldmap_interactions
        lines = overlay.lines()
        between = [l for l in lines if l.labels[0].chain != l.labels[1].chain]
        copies = {lab.chain for line in between for lab in line.labels}  # turned copies past the first two too
        assert overlay.dropped == () and len(between) >= 4 and len(copies) >= (3 if d.layout.rosette else 2)
        index = {(l.chain, l.tag): k for k, l in enumerate(d.layout.res_labels)}
        marks = [
            (f"{lab.chain}:{lab.tag}", end, figcheck.residue_places(d, index[(lab.chain, lab.tag)]))
            for line in lines
            for end, lab in zip(line.ends, line.labels)
        ]
        assert mark_findings(marks) == []
        ink = [l for l in d.figure.axes[0].lines if (l.get_gid() or "").startswith(("ia:", "hbond:"))]
        drawn = {l.get_gid(): _points(l) for l in ink}
        for line in lines:
            if line.gid in drawn and line.gid not in overlay.slots:
                off, along = _across(drawn[line.gid], line.ends)
                assert np.abs(off).max() < 1e-9 and along.min() > -1e-9 and along.max() < 1 + 1e-9
        symbol = [p for p in d.figure.axes[0].patches if (p.get_gid() or "").startswith("axis:")]
        assert bool(symbol) == bool(d.layout.rosette)
        for p in symbol:
            shape = MPath(p.get_xy())
            for line in lines:
                (x0, y0), (x1, y1) = line.ends
                t = np.linspace(0, 1, 400)[:, None]
                assert not shape.contains_points(np.c_[x0 + t * (x1 - x0), y0 + t * (y1 - y0)]).any(), line.gid
