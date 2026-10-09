"""Biological assemblies built from the file's operators, cubic point groups, and the highlight mode."""

from pathlib import Path

import numpy as np
import pytest
from helpers import pipeline

from topoplot.io import load_backbone

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
    from topoplot.symmetry import detect_symmetry

    bb = load_backbone(APO)
    s = detect_symmetry(bb)
    assert s.kind == "O" and s.label == "O" and s.n == 4
    assert len(s.protomers) == 4 and all(len(p) == 6 for p in s.protomers)


def test_tetrahedral_group_is_recognised(tmp_path):
    from itertools import product

    from test_symmetry_scenarios import build, rot

    from topoplot.symmetry import detect_symmetry

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

    from topoplot.cli import make_figure, make_layout
    from topoplot.render import MATE_GREY
    from topoplot.style import Style

    lay, sses, bb = make_layout(APO)
    fig = make_figure(APO, look=Style(highlight="asu"))
    ax = fig.axes[0]
    fills = {a.get_gid().split(":", 1)[1]: to_hex(a.get_facecolor()) for a in ax.patches
             if (a.get_gid() or "").startswith("helix:")}
    asu = [s.id for s in sses if s.chain in bb.asu_chains]
    mates = [s.id for s in sses if s.chain not in bb.asu_chains]
    assert asu and all(fills[k] != MATE_GREY for k in asu)
    assert mates and all(fills[k] == MATE_GREY for k in mates)


def test_highlight_by_chain_list():
    from matplotlib.colors import to_hex

    from topoplot.cli import make_figure
    from topoplot.render import MATE_GREY
    from topoplot.style import Style

    fig = make_figure(DATA / "2HHB.cif", look=Style(highlight="A,B"))
    for a in fig.axes[0].patches:
        gid = a.get_gid() or ""
        if gid.startswith("helix:"):
            chain = gid.split(":")[1]
            assert (to_hex(a.get_facecolor()) == MATE_GREY) == (chain in "CD"), gid


def test_bad_highlight_value_is_rejected():
    from topoplot.style import resolve_style

    with pytest.raises(ValueError, match="highlight"):
        resolve_style(overrides=["highlight=A;B"])


def test_cli_assembly_and_highlight(tmp_path, capsys):
    from topoplot.cli import main

    out = tmp_path / "apo.svg"
    assert main(["plot", str(APO), "-o", str(out), "--assembly", "asu"]) == 0
    assert main(["summary", str(APO)]) == 0
    text = capsys.readouterr().out
    assert "24 chain(s)" in text and "symmetry: O" in text
