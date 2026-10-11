"""The figure against the structure, on the committed fixtures: the fast part of scripts/na_validate.py.

Limits come from the measured campaign (docs/validation/na/report.md); a change that drops a fixture below them
moved the drawing away from the molecule."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_layout
from foldmap.io import read_model
from foldmap.nacheck import SIDE_SEEN, THRESHOLDS, contacts, judge, na_crossings, render_checks, shape
from foldmap.render import draw
from foldmap.route import route_loops

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB"]
LIMITS = {  # metric: (lowest or highest allowed, "min" or "max")
    "spearman": (0.6, "min"),
    "procrustes": (0.8, "max"),
    "angle_median": (5.0, "max"),
    "angle_max": (5.0, "max"),  # fixtures 0.0, plus a margin; Task 14 tightens it
    "joint_error_max": (5.0, "max"),
    "mirrored": (0, "max"),
    "contact_auc": (0.6, "min"),
    "overlaps": (0, "max"),
    "text_off_page": (0, "max"),
    "letter_overlaps": (0, "max"),
    "beads_missing": (0, "max"),
    "bead_colour_wrong": (0, "max"),
    "beads_false": (0, "max"),
    "beads_unmatched": (0, "max"),
    "fallbacks": (0, "max"),
}
# Campaign defects whose target (nacheck.THRESHOLDS) is 0 but which some fixtures still show: each fixture is held
# at its value measured when the metric came in (absent: 0), a regression bound. Task 14 tightens them after the
# campaign (Task 13). 4YBB: run AA:26 across the free ends AA:2-8 and AA:557-566; two contact beads overlapping;
# the "... 36 nt ..." label on stem P34 and on a run.
CAMPAIGN = {
    "na_crossings": {"4YBB": 2},
    "bead_overlaps": {"4YBB": 2},
    "gap_label_collisions": {"4YBB": 2},
}
# stems_nearer_frac (target 1.0) held at its measured value (absent: 1.0): on 30QM one contacted stem of six
NEARER = {"30QM": 0.833}
# what each fixture cannot measure (None): Spearman under 5 stems, contacts without protein, bends without a
# split duplex; every other metric must be there, so a check that stops measuring fails
BENDS = {"bend_3d", "bend_drawn", "joint_error_max"}
CONTACT = {"contact_auc", "contact_precision_at_k", "stems_nearer_frac"}
UNMEASURED = {
    "30QM": BENDS,
    "1EHZ": {"spearman"} | CONTACT | BENDS,
    "4OO8": BENDS,
    "1KX5": set(),
    "1IHF": {"spearman"},
    "3CRX": {"spearman"} | BENDS,
    "4YBB": set(),
}
METRICS = sorted(set(LIMITS) | set(CAMPAIGN) | BENDS | CONTACT | {"stems_drawn", "procrustes", "na_protein_crossings"})


@pytest.fixture(scope="module", params=FIXTURES)
def measured(request):
    path = NA / f"{request.param}.cif.gz"
    model = read_model(path)
    lay, sses, bb = make_layout(path, model=model)
    loops = route_loops(lay, sses, bb)
    fig = draw(lay, loops, sses, request.param)
    m = {**shape(lay), **contacts(lay, sses), **render_checks(lay, loops, fig, sses, model=model, bb=bb)}
    return request.param, m, lay


def test_fixture_measures_every_stem_and_every_metric_it_can(measured):
    pdb, m, lay = measured
    assert m["stems_drawn"] == len(lay.stems) > 0
    assert {k for k in METRICS if m.get(k) is None} == UNMEASURED[pdb]


@pytest.mark.parametrize("metric", sorted(LIMITS))
def test_fixture_is_drawn_true_to_its_structure(measured, metric):
    pdb, m, _ = measured
    value = m.get(metric)
    if metric in UNMEASURED[pdb]:
        assert value is None
        return
    limit, kind = LIMITS[metric]
    assert value >= limit if kind == "min" else value <= limit, f"{pdb}: {metric} {value} (limit {limit})"


@pytest.mark.parametrize("metric", sorted(CAMPAIGN))
def test_fixture_has_no_more_campaign_defects_than_when_measured(measured, metric):
    pdb, m, _ = measured
    bound = CAMPAIGN[metric].get(pdb, 0)
    assert m[metric] <= bound, f"{pdb}: {metric} {m[metric]} (held at {bound} until Task 14)"


def test_contacted_stems_stay_nearer_their_touching_elements(measured):
    pdb, m, _ = measured
    if "stems_nearer_frac" in UNMEASURED[pdb]:
        return
    assert m["stems_nearer_frac"] >= NEARER.get(pdb, 1.0)


def test_crossings_count_proper_crossings_only():
    from foldmap.nacheck import _crossings, _segments

    x = _segments([(0, 0), (2, 2)])
    assert _crossings(x, _segments([(0, 2), (2, 0)])) == 1
    assert _crossings(x, _segments([(2, 2), (3, 0)])) == 0  # meeting at an end is not a crossing
    assert _crossings(x, _segments([(0, 1), (2, 1), (2, 3), (0, 3)])) == 1


def test_judge_names_each_metric_past_its_limit():
    assert judge({"spearman": 0.9, "procrustes": 0.1, "mirrored": 0, "contact_auc": None}) == []
    failed = judge({"spearman": 0.5, "procrustes": 0.9, "mirrored": 1, "letter_overlaps": 0})
    assert failed == ["spearman 0.5 (limit 0.6)", "procrustes 0.9 (limit 0.8)", "mirrored 1 (limit 0)"]
    assert set(LIMITS) <= set(THRESHOLDS)
    assert judge({"procrustes": float("nan")}) == ["procrustes nan (limit 0.8)"]


def test_procrustes_fits_a_turned_scaled_copy_exactly_and_never_a_mirror_image():
    from foldmap.nacheck import _similarity_rmsd

    a = np.array([[0.0, 0.0], [3.0, 0.5], [1.0, 2.0], [-1.5, 1.0], [2.0, -2.0]])
    turn = np.array([[np.cos(0.7), -np.sin(0.7)], [np.sin(0.7), np.cos(0.7)]])
    assert _similarity_rmsd(a, 3 * a @ turn.T + 5) == pytest.approx(0, abs=1e-12)
    assert _similarity_rmsd(a, a * [1, -1]) > 0.5


def test_a_stem_turned_onto_its_mirror_image_is_counted():
    lay, _, _ = make_layout(NA / "1EHZ.cif.gz")
    assert shape(lay)["mirrored"] == 0
    clear = [d for d in lay.stems if d.margin >= SIDE_SEEN]
    assert clear
    clear[0].first = -clear[0].first  # its first arm drawn on the side opposite the one it is seen on
    assert shape(lay)["mirrored"] == 1


def test_a_bend_drawn_the_wrong_way_is_a_joint_error():
    lay, _, _ = make_layout(NA / "1IHF.cif.gz")
    st = lay.stems
    assert shape(lay)["joint_error_max"] == 0.0
    k = max(i for i in range(1, len(st)) if st[i].parent and st[i].parent == st[i - 1].parent)
    st[k].angle = 2 * st[k - 1].angle - st[k].angle  # the last piece reflected about its predecessor
    m = shape(lay)
    assert m["joint_error_max"] == pytest.approx(39.9, abs=0.5)  # its 20-degree turn drawn as -20
    assert any(f.startswith("joint_error_max") for f in judge(m))


def test_a_stem_the_structure_does_not_have_is_an_error():
    lay, _, _ = make_layout(NA / "1EHZ.cif.gz")
    lay.stems[0].id = "nowhere"
    with pytest.raises(ValueError, match="nowhere"):
        shape(lay)


def _drawn(pdb):
    path = NA / f"{pdb}.cif.gz"
    model = read_model(path)
    lay, sses, bb = make_layout(path, model=model)
    loops = route_loops(lay, sses, bb)
    return lay, sses, bb, loops, model, draw(lay, loops, sses, pdb)


def test_planted_crossings_are_counted():
    import copy

    lay, _, _, loops, _, _ = _drawn("4OO8")
    base = na_crossings(lay, loops)["na_crossings"]
    run = next(p for p in lay.na_paths if p.kind == "run")
    tail = np.asarray(next(p for p in lay.na_paths if p.kind == "tail").points)
    a, b = tail[len(tail) // 2 - 1], tail[len(tail) // 2]
    normal = np.array([a[1] - b[1], b[0] - a[0]]) / np.hypot(*(b - a))
    d = lay.stems[0]
    for pts in (
        [tuple((a + b) / 2 - 3 * normal), tuple((a + b) / 2 + 3 * normal)],  # across the free end
        [tuple(np.add((d.cx, d.cy), -50 * d.normal)), tuple(np.add((d.cx, d.cy), 50 * d.normal))],  # a stem's body
    ):
        planted = copy.copy(run)
        planted.points, planted.between = pts, ()
        lay.na_paths.append(planted)
        assert na_crossings(lay, loops)["na_crossings"] == base + 1
        lay.na_paths.pop()


def test_planted_bead_and_gap_label_collisions_are_counted():
    lay, sses, bb, loops, model, fig = _drawn("4YBB")
    ax = fig.axes[0]

    def count(key):
        return render_checks(lay, loops, fig, sses, model=model, bb=bb)[key]

    beads, labels = count("bead_overlaps"), count("gap_label_collisions")
    sc = next(c for c in ax.collections if (c.get_gid() or "").startswith("na:contact:"))
    drawn = np.array(sc.get_offsets())
    off = drawn.copy()
    off[1] = off[0] + [0.05, 0]  # a bead on top of another
    sc.set_offsets(off)
    assert count("bead_overlaps") == beads + 1
    sc.set_offsets(drawn)
    label = next(t for t in ax.texts if (t.get_gid() or "").startswith("na:gap:") and t.get_gid() != "na:gap:105")
    helix = max(lay.placed.values(), key=lambda e: e.length)
    label.set_position((helix.cx, helix.cy - 0.1))  # on a helix
    assert count("gap_label_collisions") == labels + 1
