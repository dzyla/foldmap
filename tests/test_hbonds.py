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
import pytest

from foldmap import hbonds
from foldmap.hbonds import HBond, Selection, find_hbonds, parse_selection
from foldmap.io import read_model

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
from foldmap.render import draw, hbond_overlay, nucleotide_points, residue_point  # noqa: E402
from foldmap.route import route_loops  # noqa: E402
from foldmap.style import Style, resolve_style  # noqa: E402


@lru_cache(maxsize=None)
def _drawn(name, *settings):
    """(figure, layout, loops, sses, bb, look) of the figure drawn with these style settings."""
    look = resolve_style(None, None, list(settings))
    lay, sses, bb = make_layout(_path(name), look=look, model=_model(name))
    loops = route_loops(lay, sses, bb)
    return draw(lay, loops, sses, look=look), lay, loops, sses, bb, look


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
    all_ = hbond_overlay(*_drawn("1LMB", "hbonds=na", "hbond_max=1000")[1:4], Style(hbonds="na", hbond_max=1000))
    few = hbond_overlay(*_drawn("1LMB", "hbonds=na", "hbond_max=1000")[1:4], Style(hbonds="na", hbond_max=7))
    assert [l.gid for l in few.lines] == [l.gid for l in all_.lines[:7]]
    lengths = [l.length for l in all_.lines]
    assert lengths == sorted(lengths) and few.total == len(all_.lines)


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
    from foldmap.interactive import _hbonds

    fig = _drawn("1LMB", "hbonds=na", "hbond_max=1000")[0]
    notes = _hbonds(fig)
    assert len(notes) == len(_lines(fig))
    by_key = {n["key"]: n for n in notes}
    gln = next(n for k, n in by_key.items() if k.startswith("3:44-1:4"))
    assert "Gln44" in gln["text"] and "DA4" in gln["text"] and "Å" in gln["text"] and len(gln["residues"]) == 1


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
