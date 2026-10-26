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


def test_beads_clear_where_runs_meet():
    """Plan 28 F6 (M7): _unbunch slides each run against the runs placed before it, which stay put, so where runs meet
    the later one could only keep the best spot it found (6QIK: A:G1447 and y:G85, 0.28 apart where three runs meet;
    a bead is 0.45 wide). Now the other run slides instead (na_labels._clear_junctions). Planted on 4OO8: a run redrawn
    across another's first nucleotide, its middle nucleotide 0.2 from it; afterwards no nucleotide of the two is under
    a bead's width from one of the other."""
    from foldmap.na_labels import _clear_junctions
    from foldmap.na_paths import _BEAD, _spaced

    path = Path(__file__).parent / "data" / "na" / "4OO8.cif.gz"
    lay, sses, bb = make_layout(path)
    route_loops(lay, sses, bb)
    runs = [p for p in lay.na_paths if p.kind == "run" and p.nts]
    q = runs[0]
    p = next(r for r in runs[1:] if len(r.nts) % 2 == 1 and len(r.nts) >= 3)
    x, y = q.nt_xy[0]
    half = (len(p.nts) + 1) * 0.55 / 2
    p.points = [(x - half, y + 0.2), (x + half, y + 0.2)]
    p.nt_xy = _spaced(p.points, len(p.nts))

    def nearest() -> float:
        return min(np.hypot(a[0] - b[0], a[1] - b[1]) for a in p.nt_xy for b in q.nt_xy)

    assert nearest() < _BEAD
    _clear_junctions(lay, lay.na_paths, runs, {key: d for d in lay.stems for key in d.axial})
    assert nearest() >= _BEAD - 1e-9
