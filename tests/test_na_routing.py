# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Nucleic-acid routing on the cases plan 28's diagnosis traced (na-diagnosis.md, fixes F4-F6): gaps, joints and
beads, each on a committed fixture."""

from pathlib import Path

import numpy as np

from foldmap.cli import make_layout
from foldmap.io import read_model
from foldmap.nacheck import _crossings, _segments, render_checks
from foldmap.render import draw
from foldmap.route import route_loops

NA = Path(__file__).parent / "data" / "na_fixes"  # plan 28 fixtures, apart from the kernels' tests/data/na


def _drawn(name):
    path = NA / f"{name}.cif.gz"
    model = read_model(path)
    lay, sses, bb = make_layout(path, model=model)
    loops = route_loops(lay, sses, bb)
    return lay, sses, bb, loops, model, draw(lay, loops, sses, name)


def test_gap_does_not_cross_the_tail_it_joins():
    """Plan 28 F4 (M5): in 4AS1, B:32-35 are not modelled, so a gap joins B:31 to B:36, whose 5' free end B:36-40 is
    drawn first. It turned toward its touching protein, away from B:31, and the gap, arriving along the tail's outward
    heading, went up past the tail's base and hooked back into its tip, crossing it. Now the 5' end after a gap turns
    toward where the molecule comes from as well (_toward's target), and the gap may arrive along the chord: it
    crosses no free end, and the figure has no nucleic-acid crossing."""
    lay, sses, bb, loops, model, fig = _drawn("4AS1")
    na = lay.nucleic
    (g,) = [g for g in na.gaps if na.strands[g.before].labels[-1].seq == 31]
    gap = next(p for p in lay.na_paths if p.kind == "gap" and p.strand == g.before)
    tails = [p for p in lay.na_paths if p.kind == "tail"]
    assert tails and not [t.strand for t in tails if _crossings(_segments(gap.points), _segments(t.points))]
    assert render_checks(lay, loops, fig, sses, model=model, bb=bb)["na_crossings"] == 0
    assert np.isfinite(np.asarray(gap.points, float)).all()
