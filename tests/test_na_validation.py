"""The figure against the structure, on the committed fixtures: the fast part of scripts/na_validate.py.

Limits come from the measured campaign (docs/validation/na/report.md); a change that drops a fixture below them
moved the drawing away from the molecule."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_layout
from foldmap.io import read_model
from foldmap.nacheck import THRESHOLDS, contacts, judge, render_checks, shape
from foldmap.render import draw
from foldmap.route import route_loops

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB"]
LIMITS = {  # metric: (lowest or highest allowed, "min" or "max")
    "spearman": (0.6, "min"),
    "procrustes": (0.8, "max"),
    "angle_median": (5.0, "max"),
    "mirrored": (0, "max"),
    "contact_auc": (0.55, "min"),
    "overlaps": (0, "max"),
    "text_off_page": (0, "max"),
    "letter_overlaps": (0, "max"),
    "beads_missing": (0, "max"),
    "bead_colour_wrong": (0, "max"),
    "fallbacks": (0, "max"),
}


@pytest.fixture(scope="module", params=FIXTURES)
def measured(request):
    path = NA / f"{request.param}.cif.gz"
    model = read_model(path)
    lay, sses, bb = make_layout(path, model=model)
    loops = route_loops(lay, sses, bb)
    fig = draw(lay, loops, sses, request.param)
    m = {**shape(lay), **contacts(lay, sses), **render_checks(lay, loops, fig, sses, model=model, bb=bb)}
    return request.param, m


@pytest.mark.parametrize("metric", sorted(LIMITS))
def test_fixture_is_drawn_true_to_its_structure(measured, metric):
    pdb, m = measured
    value = m.get(metric)
    if value is None:  # not measurable here: too few stems, no protein
        return
    limit, kind = LIMITS[metric]
    assert value >= limit if kind == "min" else value <= limit, f"{pdb}: {metric} {value} (limit {limit})"


def test_judge_names_each_metric_past_its_limit():
    assert judge({"spearman": 0.9, "procrustes": 0.1, "mirrored": 0, "contact_auc": None}) == []
    failed = judge({"spearman": 0.5, "procrustes": 0.9, "mirrored": 1, "letter_overlaps": 0})
    assert failed == ["spearman 0.5 (limit 0.6)", "procrustes 0.9 (limit 0.8)", "mirrored 1 (limit 0)"]
    assert set(LIMITS) <= set(THRESHOLDS)


def test_procrustes_fits_a_turned_scaled_copy_exactly_and_never_a_mirror_image():
    from foldmap.nacheck import _similarity_rmsd

    a = np.array([[0.0, 0.0], [3.0, 0.5], [1.0, 2.0], [-1.5, 1.0], [2.0, -2.0]])
    turn = np.array([[np.cos(0.7), -np.sin(0.7)], [np.sin(0.7), np.cos(0.7)]])
    assert _similarity_rmsd(a, 3 * a @ turn.T + 5) == pytest.approx(0, abs=1e-12)
    assert _similarity_rmsd(a, a * [1, -1]) > 0.5


def test_a_stem_turned_onto_its_mirror_image_is_counted():
    lay, _, _ = make_layout(NA / "1EHZ.cif.gz")
    assert shape(lay)["mirrored"] == 0
    clear = [d for d in lay.stems if d.margin >= 2.75]
    assert clear
    clear[0].first = -clear[0].first  # its first arm drawn on the side opposite the one it is seen on
    assert shape(lay)["mirrored"] == 1
