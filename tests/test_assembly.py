"""Biological assemblies built from the file's operators, cubic point groups, and the highlight mode."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.io import load_backbone

DATA = Path(__file__).parent / "data"
APO = DATA / "assemblies" / "7A4M.cif"  # mouse apoferritin: one subunit deposited, 24 by the O operators


def test_auto_builds_the_24mer_and_remembers_the_deposited_subunit():
    bb = load_backbone(APO)
    assert len({l.chain for l in bb.labels}) == 24
    assert len(bb.asu_chains) == 1 and bb.asu_chains <= {l.chain for l in bb.labels}


def test_asu_keeps_the_file_as_deposited():
    bb = load_backbone(APO, assembly="asu")
    assert {l.chain for l in bb.labels} == {"A"} and bb.asu_chains == {"A"}


def test_unknown_assembly_is_an_error():
    with pytest.raises(ValueError, match="assembly"):
        load_backbone(APO, assembly="7")


def test_single_copy_files_are_unchanged(ubq):
    bb = load_backbone(ubq)
    assert {l.chain for l in bb.labels} == {"A"} and bb.asu_chains == {"A"}


def test_apoferritin_is_octahedral():
    from foldmap.symmetry import detect_symmetry

    bb = load_backbone(APO)
    s = detect_symmetry(bb)
    assert s.kind == "O" and s.label == "O" and s.n == 4
    assert len(s.protomers) == 4 and all(len(p) == 6 for p in s.protomers)


def test_tetrahedral_group_is_recognised(tmp_path):
    from itertools import product

    from foldmap.symmetry import detect_symmetry
    from test_symmetry_scenarios import build, rot

    gens = [rot((1, 1, 1), 120), rot((0, 0, 1), 180)]
    group = [np.eye(3)]
    while True:  # close the group under the generators
        new = [g @ h for g, h in product(gens, group)]
        add = [m for m in new if not any(np.allclose(m, x) for x in group)]
        if not add:
            break
        group += add
    assert len(group) == 12
    path = build(tmp_path, [(m, np.zeros(3)) for m in group], offset=(16.0, 7.0, 2.0))
    s = detect_symmetry(load_backbone(path))
    assert s.kind == "T" and s.n == 3


def test_highlight_keeps_the_asu_coloured_and_greys_the_mates(tmp_path):
    from matplotlib.colors import to_hex

    from foldmap.cli import make_figure, make_layout
    from foldmap.render import MATE_GREY
    from foldmap.style import Style

    lay, sses, bb = make_layout(APO)
    fig = make_figure(APO, look=Style(highlight="asu"))
    ax = fig.axes[0]
    fills = {
        a.get_gid().split(":", 1)[1]: to_hex(a.get_facecolor())
        for a in ax.patches
        if (a.get_gid() or "").startswith("helix:")
    }
    asu = [s.id for s in sses if s.chain in bb.asu_chains]
    mates = [s.id for s in sses if s.chain not in bb.asu_chains]
    assert asu and all(fills[k] != MATE_GREY for k in asu)
    assert mates and all(fills[k] == MATE_GREY for k in mates)


def test_highlight_by_chain_list():
    from matplotlib.colors import to_hex

    from foldmap.cli import make_figure
    from foldmap.render import MATE_GREY
    from foldmap.style import Style

    fig = make_figure(DATA / "2HHB.cif", look=Style(highlight="A,B"))
    for a in fig.axes[0].patches:
        gid = a.get_gid() or ""
        if gid.startswith("helix:"):
            chain = gid.split(":")[1]
            assert (to_hex(a.get_facecolor()) == MATE_GREY) == (chain in "CD"), gid


def test_bad_highlight_value_is_rejected():
    from foldmap.style import resolve_style

    with pytest.raises(ValueError, match="highlight"):
        resolve_style(overrides=["highlight=A;B"])


def test_cli_assembly_and_highlight(tmp_path, capsys):
    from foldmap.cli import main

    out = tmp_path / "apo.svg"
    assert main(["plot", str(APO), "-o", str(out), "--assembly", "asu"]) == 0
    assert main(["summary", str(APO)]) == 0
    text = capsys.readouterr().out
    assert "24 chain(s)" in text and "symmetry: O" in text


def _model(chains):
    """A gemmi Model of one-atom residues: {chain name: [xyz, ...]}."""
    import gemmi

    model = gemmi.Model("1")
    for name, xyz in chains.items():
        chain = gemmi.Chain(name)
        for k, p in enumerate(xyz):
            res = gemmi.Residue()
            res.name, res.seqid = "ALA", gemmi.SeqId(k + 1, " ")
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = "CA", gemmi.Element("C"), gemmi.Position(*p)
            res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    return model


def test_a_copy_lying_on_a_chain_already_present_is_dropped():
    """A copy with >= 90 % of its atoms within 1 Å of one chain already present is the same molecule seen twice
    (an operator mapping a chain on a special position onto itself or its partner); a copy that only partly
    overlaps, or lies on the union of two chains, is kept."""
    from foldmap.io import _drop_coincident

    rng = np.random.default_rng(0)
    a = rng.uniform(0, 40, (20, 3))
    b = a + 20.0
    on_a = a + 0.4  # every atom 0.69 Å from its twin
    half = np.vstack([a[:10] + 0.4, a[10:] + 6.0])  # half of it coincides: a real neighbour, kept
    split = np.vstack([a[:9] + 0.4, b[9:] + 0.4])  # 45 % on A and 55 % on B: on no single chain
    model = _model({"A1": a, "B1": b, "A2": on_a, "C2": half, "D2": split})
    assert _drop_coincident(model) == ["A2"]
    assert [c.name for c in model] == ["A1", "B1", "C2", "D2"]


def test_holliday_junction_assembly_keeps_each_dna_strand_once():
    """3CRX: the two-fold operator lays a copy of every DNA strand on a deposited one (C2 on E1, D2 on F1, ...).
    The assembly keeps the four Cre and the four strands, which pair into the junction's four arms."""
    from scipy.spatial import cKDTree

    from foldmap.io import load_nucleic, read_model

    path = DATA / "na" / "3CRX.cif.gz"
    model = read_model(path)
    bb = load_backbone(path, model=model)
    assert len({l.chain for l in bb.labels}) == 4  # four Cre
    na = load_nucleic(path, bb, model=model)
    assert len(na.strands) == 4
    for i, s in enumerate(na.strands):
        for t in na.strands[:i]:
            d, _ = cKDTree(t.p).query(s.p)
            assert (d < 1.5).mean() < 0.5, (s.chain, t.chain)
    assert len(na.structure.stems) == 4
