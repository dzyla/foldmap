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
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB", "1RVA", "4LD0"]
# Thin margins: procrustes 0.5 against 4LD0's 0.499, spearman 0.85 against 4OO8 and 1KX5 at 0.85.
# contact_auc does not respond to the contact pull (layout._W_CONTACT): test_na_layout.py's displacement test is
# the guard. The angle limits are a regression guard on plain frames, vacuous as fidelity (M2).
# metric: (lowest or highest allowed, "min" or "max"): the worst fixture value, rounded safe (Task 14)
LIMITS = {
    "spearman": (0.85, "min"),
    "procrustes": (0.5, "max"),
    "angle_median": (0.05, "max"),
    "angle_max": (0.05, "max"),  # fixtures 0.0 (reports round to 0.1)
    "joint_error_max": (0.05, "max"),
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
# Campaign defects whose target (nacheck.THRESHOLDS) is 0: a fixture could be held at a measured value here (absent:
# 0), a regression bound. The final campaign (Task 14) leaves every fixture at 0, so none is held. The 4YBB campaign
# entry (not this fixture) keeps fallbacks 2, na_crossings 1 and bead_overlaps 1 as accepted failures
# (docs/validation/na/review.md).
CAMPAIGN = {
    "na_crossings": {},  # every fixture at 0
    "bead_overlaps": {},  # every fixture at 0
    "gap_label_collisions": {},  # every fixture at 0
}
# stems_nearer_frac (target 1.0) held at its measured value (absent: 1.0). 30QM is drawn unrolled, its two protomers
# side by side (the C4 revision of Task 13b step 3: unrolling tears none of its duplexes): one stem held by both
# protomers (P3, touched by A:288-306, A:346-362 and B:37-57) sits between them, nearer neither than the median.
NEARER: dict[str, float] = {"30QM": 0.833}
# LIMITS held at a fixture's measured value until a later step: 3CRX's four-Cre assembly (C5) is seen edge-on to
# the junction plane (arms 3, 4, 1, 2 in a row), drawn as an X (C2 and C3 leave it at 0.817: it has no runs or
# joints). 1RVA's contact_auc hold (0.58) went with C3: its two pieces meet across the kink (0.609).
# 30QM, unrolled side by side (see NEARER): its stems are measured against the plain projection of the same
# orientation (nacheck.measuring_frame), so the angles and the distance ranks measure the unrolling itself, which
# turns each protomer's half to face the viewer (angle_median 35.7, angle_max 109.6, spearman 0.521, procrustes
# 0.776: the limit set from the final campaign, 0.50, excludes the unrolled figure).
HELD = {
    ("3CRX", "procrustes"): 0.817,
    ("30QM", "procrustes"): 0.776,
    ("30QM", "angle_median"): 35.7,
    ("30QM", "angle_max"): 109.6,
    ("30QM", "spearman"): 0.521,
}
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
    "1RVA": {"spearman", "procrustes"},
    "4LD0": {"spearman"} | BENDS,
}
METRICS = sorted(
    set(LIMITS)
    | set(CAMPAIGN)
    | BENDS
    | CONTACT
    | {"stems_drawn", "procrustes", "na_protein_crossings", "na_run_crossings"}
)


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
    limit = HELD.get((pdb, metric), limit)
    assert value >= limit if kind == "min" else value <= limit, f"{pdb}: {metric} {value} (limit {limit})"


@pytest.mark.parametrize("metric", sorted(CAMPAIGN))
def test_fixture_has_no_more_campaign_defects_than_when_measured(measured, metric):
    pdb, m, _ = measured
    bound = CAMPAIGN[metric].get(pdb, 0)
    assert m[metric] <= bound, f"{pdb}: {metric} {m[metric]} (held at {bound})"


def test_contacted_stems_stay_nearer_their_touching_elements(measured):
    pdb, m, _ = measured
    if "stems_nearer_frac" in UNMEASURED[pdb]:
        return
    assert m["stems_nearer_frac"] >= NEARER.get(pdb, 1.0)


def test_every_hold_still_fails_its_real_limit(measured):
    """A hold (HELD, NEARER, CAMPAIGN) excuses a fixture from nacheck.THRESHOLDS at its measured value. Once the
    drawing meets the real limit the hold must go: this fails until it is removed (strict-xfail style)."""
    pdb, m, _ = measured
    holds = [metric for (p, metric) in HELD if p == pdb]
    holds += ["stems_nearer_frac"] if pdb in NEARER else []
    holds += [metric for metric, bound in CAMPAIGN.items() if pdb in bound]
    for metric in holds:
        assert judge({metric: m[metric]}), f"{pdb}: {metric} {m[metric]} meets THRESHOLDS now: remove its hold"


def test_crossings_count_proper_crossings_only():
    from foldmap.nacheck import _crossings, _segments

    x = _segments([(0, 0), (2, 2)])
    assert _crossings(x, _segments([(0, 2), (2, 0)])) == 1
    assert _crossings(x, _segments([(2, 2), (3, 0)])) == 0  # meeting at an end is not a crossing
    assert _crossings(x, _segments([(0, 1), (2, 1), (2, 3), (0, 3)])) == 1


def test_judge_names_each_metric_past_its_limit():
    assert judge({"spearman": 0.9, "procrustes": 0.1, "mirrored": 0, "contact_auc": None}) == []
    failed = judge({"spearman": 0.5, "procrustes": 0.9, "mirrored": 1, "letter_overlaps": 0})
    assert failed == ["spearman 0.5 (limit 0.6)", "procrustes 0.9 (limit 0.5)", "mirrored 1 (limit 0)"]
    assert set(LIMITS) <= set(THRESHOLDS)
    assert judge({"procrustes": float("nan")}) == ["procrustes nan (limit 0.5)"]


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
        [tuple(np.add((d.cx, d.cy), -d.width * d.normal)), tuple(np.add((d.cx, d.cy), d.width * d.normal))],  # a body
    ):
        planted = copy.copy(run)
        planted.points, planted.between = pts, ()
        lay.na_paths.append(planted)
        assert na_crossings(lay, loops)["na_crossings"] == base + 1
        lay.na_paths.pop()


def test_runs_crossing_runs_are_counted_and_reported_only():
    """na_run_crossings: pairs of runs, joints and gaps that cross at least once (the two strands of one joint, which
    may cross where the duplex carries on, excepted). Reported, not judged (THRESHOLDS None)."""
    import copy

    lay, _, _, loops, _, _ = _drawn("4OO8")
    base = na_crossings(lay, loops)
    run = next(p for p in lay.na_paths if p.kind == "run" and len(p.points) > 3)
    a, b = np.asarray(run.points[1]), np.asarray(run.points[2])
    mid, normal = (a + b) / 2, np.array([a[1] - b[1], b[0] - a[0]]) / np.hypot(*(b - a))
    planted = copy.copy(run)
    end = tuple(mid + 0.2 * normal)
    planted.points, planted.between = [tuple(mid - 0.2 * normal), tuple(mid - 0.1 * normal), end, end], ()
    lay.na_paths.append(planted)
    after = na_crossings(lay, loops)
    lay.na_paths.pop()
    assert after["na_run_crossings"] == base["na_run_crossings"] + 1
    assert after["na_crossings"] == base["na_crossings"]
    assert THRESHOLDS["na_run_crossings"] is None and not judge({"na_run_crossings": 99})


def test_the_two_strands_of_one_joint_may_cross_without_counting():
    """na_run_crossings excepts the two strands of one joint (they cross where the duplex carries on), but not two
    joints whose pieces are not placed: a missing piece must not make two joints look like the same one."""
    import copy

    lay, _, _, loops, _, _ = _drawn("1KX5")
    base = na_crossings(lay, loops)["na_run_crossings"]
    joint = next(p for p in lay.na_paths if p.kind == "joint")
    a, b = np.asarray(joint.points[1]), np.asarray(joint.points[2])  # its first segment between the faces
    mid, normal = (a + b) / 2, np.array([a[1] - b[1], b[0] - a[0]]) / max(np.hypot(*(b - a)), 1e-9)

    def across(path, centre, n):  # a copy of path drawn as a short stroke across centre along n
        out = copy.copy(path)
        p, q = tuple(centre - 0.05 * n), tuple(centre + 0.05 * n)
        out.points = [p, p, q, q]
        return out

    lay.na_paths.append(across(joint, mid, normal))  # the same two pieces: exempt
    assert na_crossings(lay, loops)["na_run_crossings"] == base
    lay.na_paths.pop()
    far, e = np.array([1e3, 1e3]), np.array([1.0, 0.0])
    one, two = across(joint, far, e), across(joint, far, e[::-1])
    one.between = two.between = ((99, 0), (99, 1))  # pieces not placed
    lay.na_paths += [one, two]
    after = na_crossings(lay, loops)["na_run_crossings"]
    del lay.na_paths[-2:]
    assert after == base + 1


def test_a_run_through_its_own_stem_is_counted():
    """A run is trimmed only by its two rail stretches (its first and last segments), so a drawing that leaves its
    rail straight into its own stem's core and out of the far side counts as a pass through that stem. Trimmed by
    any first stretch inside its own stems, as before, the whole pass was dropped (2GIS: an X over P3 and P2)."""
    import copy

    lay, _, _, loops, _, _ = _drawn("4OO8")
    base = na_crossings(lay, loops)["na_crossings"]
    at = {key: d for d in lay.stems for key in d.axial}
    run = next(p for p in lay.na_paths if p.kind == "run")
    d = at[run.between[0]]
    centre = np.array([d.cx, d.cy])
    across = -d.normal * d.side[run.between[0]]  # from its rail over the axis to the other side
    far = centre + across * (d.width + 0.5)
    planted = copy.copy(run)
    end = tuple(far + across * 0.5)
    planted.points = [tuple(d.point(run.between[0])), tuple(centre), tuple(far), end, end]
    lay.na_paths.append(planted)
    assert na_crossings(lay, loops)["na_crossings"] == base + 1


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


def _na_validate():
    import sys

    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
    try:
        import na_validate
    finally:
        sys.path.pop(0)
    return na_validate


def _reference(pdb) -> set[frozenset]:
    import json

    ref = json.loads((NA / "reference.json").read_text(encoding="utf-8"))[pdb]
    return {frozenset((tuple(a), tuple(b))) for a, b, _ in ref}


def test_drawn_chains_map_back_to_the_deposited_chains_they_copy():
    """Each drawn chain maps to exactly one (deposited chain, copy), by its name, its residues and its shape; a
    name that fits no deposited chain, or fits two (deposited names ending in digits, identical twins), is left
    unmapped (the entry then fails)."""
    from foldmap.model import ResLabel

    nv = _na_validate()
    rng = np.random.default_rng(2)
    xyz = {"A": rng.uniform(0, 30, (4, 3)), "B": rng.uniform(0, 30, (4, 3))}
    res = {c: {(k, "", "G") for k in range(1, 5)} for c in xyz}
    c1 = {(c, k, ""): xyz[c][k - 1] for c in xyz for k in range(1, 5)}

    def drawn(name, src, shift=0.0):
        return {name: [(ResLabel(name, k, "", "G"), xyz[src][k - 1] + shift) for k in range(1, 5)]}

    assert nv._chain_map({**drawn("A", "A"), **drawn("B", "B")}, res, c1, True) == (
        {"A": ("A", ""), "B": ("B", "")},
        [],
    )
    got = nv._chain_map({**drawn("A1", "A"), **drawn("A2", "A", 40.0), **drawn("X1", "A")}, res, c1, False)
    assert got == ({"A1": ("A", "1"), "A2": ("A", "2")}, ["X1"])
    assert nv._chain_map(drawn("B1", "A"), res, c1, False) == ({}, ["B1"])  # B's residues, A's shape
    # deposited A and A1: drawn A11 is A copy 11 or A1 copy 1; by shape it is A1's copy, A's twin is ambiguous
    res1, c11 = {**res, "A1": res["A"]}, {**c1, **{("A1", k, ""): xyz["B"][k - 1] for k in range(1, 5)}}
    assert nv._chain_map(drawn("A11", "B"), res1, c11, False) == ({"A11": ("A1", "1")}, [])
    twins = {**c1, **{("A1", k, ""): xyz["A"][k - 1] + 5.0 for k in range(1, 5)}}
    assert nv._chain_map(drawn("A11", "A"), res1, twins, False) == ({}, ["A11"])


def _moved(strand, chain, shift):
    """The strand as copy `chain`, moved by `shift` Å."""
    from dataclasses import replace

    labels = [l._replace(chain=chain) for l in strand.labels]
    return replace(strand, chain=chain, labels=labels, p=strand.p + shift, c1=strand.c1 + shift)


def _two_copies(na, shift):
    """na drawn twice, as copy 1 and copy 2 (moved by `shift` Å), each with its own pairs."""
    n = len(na.strands)
    strands = [*na.strands, *(_moved(s, s.chain[:-1] + "2", shift) for s in na.strands)]
    pairs = [*na.base_pairs, *(p._replace(a=p.a + n, b=p.b + n) for p in na.base_pairs)]
    return strands, pairs


def test_a_pair_counts_only_within_one_copy_or_one_operator():
    """3CRX drawn as two copies. A pair rewired from C1·D1 to C1·D2 joins C to D, as the deposited model does,
    but to the twin of D in the other copy: wrong. Copy 2 drawn on top of copy 1 (one operator under two
    numbers, as gemmi numbers one operator's copies in two generators) pairs the same molecule: right."""
    from foldmap.model import nucleic_from_pairs

    nv = _na_validate()
    path, theirs = NA / "3CRX.cif.gz", _reference("3CRX")
    asu = read_model(path, "asu")[0]
    na = make_layout(path)[0].nucleic
    n = len(na.strands)
    strands, pairs = _two_copies(na, 60.0)
    d = next(k for k, s in enumerate(na.strands) if s.chain == "D1")
    k = next(k for k, p in enumerate(pairs) if p.b == d and strands[p.a].chain == "C1")
    pairs[k] = pairs[k]._replace(b=d + n)  # C1 paired with D2
    m = nv.drawn_pairing(nucleic_from_pairs(strands, pairs, {}, {}), theirs, asu, False)
    assert m["unmapped_chains"] == [] and m["pairs_foldmap"] == len(pairs)
    assert m["pair_precision"] == round(1 - 1 / len(pairs), 3), m
    assert len(m["extra_pairs"]) == 1 and {x[0] for x in m["extra_pairs"][0]} == {"C1", "D2"}, m
    d2 = next(k for k, s in enumerate(strands) if s.chain == "D2")
    strands[d2] = _moved(na.strands[d], "D2", 0.0)  # D2 is D1's molecule: one operator, two numbers
    m = nv.drawn_pairing(nucleic_from_pairs(strands, pairs, {}, {}), theirs, asu, False)
    assert m["pair_precision"] == 1.0, m


def test_a_copy_holding_one_nucleotide_twice_is_unmapped():
    from dataclasses import replace

    nv = _na_validate()
    path, theirs = NA / "3CRX.cif.gz", _reference("3CRX")
    na = make_layout(path)[0].nucleic
    twice = replace(na, strands=[*na.strands, na.strands[0]])
    m = nv.drawn_pairing(twice, theirs, read_model(path, "asu")[0], False)
    assert m["unmapped_chains"] == [na.strands[0].chain]


def test_pairing_is_scored_on_the_pairs_the_figure_draws(monkeypatch):
    """3CRX: the assembly figure draws the deposited pairs once per strand (C5), scored like the deposited model.
    With its coincident copies kept, the pairing alternates between twins (C1 with F2 where the deposited
    structure pairs C with D): the deposited model still scores 1.0, the drawn pairs do not."""
    import foldmap.io as fio

    nv = _na_validate()
    path, theirs = NA / "3CRX.cif.gz", _reference("3CRX")
    asu = read_model(path, "asu")[0]
    assert nv.pairing(path, theirs)["deposited_pair_f1"] == 1.0
    lay, _, _ = make_layout(path)
    m = nv.drawn_pairing(lay.nucleic, theirs, asu, False)
    assert (m["pair_f1"], m["stem_recall"], m["stems_foldmap"], m["unmapped_chains"]) == (1.0, 1.0, 4, [])
    monkeypatch.setattr(fio, "_drop_coincident", lambda model, **kw: [])
    lay, _, _ = make_layout(path)
    assert len(lay.nucleic.strands) == 8
    m = nv.drawn_pairing(lay.nucleic, theirs, asu, False)
    assert m["pair_precision"] < 0.9 and m["stem_recall"] < 0.5, m


def test_a_cropped_figure_is_scored_on_the_pairs_it_draws_both_nucleotides_of():
    """4YBB: auto crops the rRNA to the stems the proteins touch; the pairs it leaves out are not missed."""
    nv = _na_validate()
    path, theirs = NA / "4YBB.cif.gz", _reference("4YBB")
    lay, _, _ = make_layout(path)
    m = nv.drawn_pairing(lay.nucleic, theirs, read_model(path, "asu")[0], True)
    whole = nv.pairing(path, theirs)  # 0.977 precision, 0.981 recall on the whole deposited fixture
    assert m["pairs_reference"] < len(theirs) and m["pair_recall"] >= whole["deposited_pair_recall"], m
    assert m["pair_precision"] >= 0.97 and m["unmapped_chains"] == [], m


def test_an_unrolled_frame_is_measured_against_the_plain_projection(monkeypatch):
    """1RVA drawn unrolled (as before C4): its kinked duplex becomes an antiparallel hairpin, bend -179.7°. Seen
    through the unrolled frame itself the 3D bend agreed (-179.7°, joint error 0). Shape is measured through a
    plain Frame of the same orientation at the nucleic acid, which sees the kink as it is and the hairpin as an
    error (the duplex is seen foreshortened, its kink as about 104°)."""
    import foldmap.cli as cli
    from foldmap.frame import CylinderFrame, Frame
    from foldmap.nacheck import measuring_frame

    monkeypatch.setattr(cli, "_na_across_protomers", lambda *a: False)
    lay, _, _ = make_layout(NA / "1RVA.cif.gz")
    assert isinstance(lay.frame, CylinderFrame)
    plain = measuring_frame(lay)
    assert type(plain) is Frame
    centre = np.vstack([s.p for s in lay.nucleic.strands]).mean(axis=0)
    assert np.allclose([plain.u, plain.v, plain.w], lay.frame.basis(centre))
    m = shape(lay)
    assert abs(m["bend_drawn"]) > 170
    assert abs(m["bend_3d"]) < 110 and m["joint_error_max"] > 60, m  # measured: -104.0 seen, error 75.7
    assert any(f.startswith("joint_error_max") for f in judge(m))


def test_a_plain_frame_is_measured_as_drawn():
    lay, _, _ = make_layout(NA / "1EHZ.cif.gz")
    from foldmap.nacheck import measuring_frame

    assert measuring_frame(lay) is lay.frame


def _score_synthetic(ours_ij, theirs_ij, n=80):
    from foldmap.model import NAStrand, ResLabel
    from foldmap.nucleic import BasePair

    nv = _na_validate()
    strand = NAStrand("A", [ResLabel("A", k + 1, "", "G") for k in range(n)], np.zeros((n, 3)), np.zeros((n, 3)), True)
    where = {(0, i): ("", ("A", i + 1, "")) for i in range(n)}
    ours = [BasePair(0, i, 0, j, "WC") for i, j in ours_ij]
    theirs = {frozenset((("A", i + 1, ""), ("A", j + 1, ""))) for i, j in theirs_ij}
    return nv._score([strand], ours, theirs, where, drawn_only=True)


def test_stem_recall_matches_a_stem_by_most_of_its_pairs():
    """M4, 5F9R: RNApolis' A29-36/51-58 against foldmap's A29-38/49-58, which carries on across a 1x1 mismatch
    (A38·A49, the better pair by its geometry): every reference pair is drawn, yet the old rule (both arms' ends
    within one nucleotide) failed it. Stems are matched one to one, a reference stem recalled when its match holds
    at least two thirds of its pairs."""
    stem = [(28 + k, 57 - k) for k in range(8)]
    m = _score_synthetic([*stem, (37, 48)], stem)
    assert (m["stems_reference"], m["stem_recall"]) == (1, 1.0)
    m = _score_synthetic(stem[:5], stem)  # five of eight: under two thirds
    assert m["stem_recall"] == 0.0
    m = _score_synthetic(stem[:6], stem)
    assert m["stem_recall"] == 1.0


def test_one_foldmap_stem_recalls_one_reference_stem_only():
    # two reference stems split by a 3x3 loop, drawn by foldmap as one stem (as a merged helix would be)
    outer, inner = [(10 + k, 69 - k) for k in range(6)], [(19 + k, 60 - k) for k in range(6)]
    theirs = outer + inner
    ours = outer + [(16, 63), (17, 62), (18, 61)] + inner  # carried on through what the reference leaves unpaired
    m = _score_synthetic(ours, theirs)
    assert (m["stems_reference"], m["stems_foldmap"], m["stem_recall"]) == (2, 1, 0.5)


def test_stem_recall_leaves_out_reference_pairs_closing_too_short_a_loop():
    """M4, 5F9R: C73·G76 closes a 2-nt loop, which no foldmap pair may (MIN_LOOP): it does not make a reference
    stem of the isolated pair stacked on it."""
    m = _score_synthetic([], [(70, 76), (72, 75)])
    assert m["stems_reference"] == 0 and m["pairs_reference"] == 2
    m = _score_synthetic([(20, 40), (21, 39)], [(20, 40), (21, 39), (70, 76), (72, 75)])
    assert (m["stems_reference"], m["stem_recall"]) == (1, 1.0)


def test_a_polyline_is_broken_where_it_jumps_across_an_unrolled_seam():
    nv = _na_validate()
    q = np.array([[0.0, 0.0], [3.0, 0.0], [90.0, 1.0], [93.0, 1.0]])
    out = nv.seam_breaks(q, 40.0)
    assert out.shape == (5, 2) and np.isnan(out[2]).all()
    assert np.array_equal(out[[0, 1, 3, 4]], q)
    assert np.array_equal(nv.seam_breaks(q, None), q)  # a plain frame has no seam


def test_a_panel_title_is_wrapped_to_its_panel():
    nv = _na_validate()
    title = "drawn layout (stems coloured, arrows: first arm 5'->3')"
    lines = nv.wrap_title(title, 2.4, 9).split("\n")
    assert len(lines) > 1 and all(len(s) <= 2.4 / (9 * 0.6 / 72) for s in lines)
    assert nv.wrap_title(title, 8.0, 9) == title


def test_the_contact_sheet_of_an_unrolled_frame_draws_no_seam_and_fits_its_titles(monkeypatch, tmp_path):
    """C14: an unrolled (CylinderFrame) figure's 3D panel joined points on either side of the unrolling's seam by
    long lines across the panel (30QM, 1RVA, 4LD0); its titles ran into each other on narrow panels (4LVV)."""
    import matplotlib.pyplot as plt

    import foldmap.cli as cli
    from foldmap.frame import CylinderFrame
    from foldmap.layout import SCALE

    nv = _na_validate()
    monkeypatch.setattr(cli, "_na_across_protomers", lambda *a: False)
    lay, sses, bb = make_layout(NA / "1RVA.cif.gz")
    assert isinstance(lay.frame, CylinderFrame)
    fig = draw(lay, route_loops(lay, sses, bb), sses, "1RVA")
    made, figure = [], plt.figure
    monkeypatch.setattr(plt, "figure", lambda *a, **k: made.append(figure(*a, **k)) or made[-1])
    monkeypatch.setattr(plt, "close", lambda *a: None)
    assert nv.contact_sheet("1RVA", lay, sses, bb, fig, {"failed": []}, tmp_path).is_file()
    sheet = made[-1]
    a0, a1, a2 = sheet.axes
    for line in a2.get_lines():
        q = np.asarray(line.get_xydata(), float)
        step = np.linalg.norm(np.diff(q, axis=0), axis=1)
        assert not (step[np.isfinite(step)] > lay.frame.period * SCALE / 2).any()
    renderer = sheet.canvas.get_renderer()
    for a in (a0, a1, a2):
        assert a.title.get_window_extent(renderer).width <= a.get_window_extent(renderer).width * 1.05 + 1


def test_single_stranded_contacts_are_reported_not_judged():
    """ss_contacts_nearer: the fraction of touched single-stranded nucleotides (free ends, loops, runs, traces) drawn
    nearer their touching elements than the median element; reported, not judged (no stem metric saw 1URN's loop
    drawn away from U1A)."""
    lay, sses, bb = make_layout(NA / "1URN.cif.gz")
    assert contacts(lay, sses)["ss_contacts_nearer"] is None  # not routed: nothing drawn single-stranded yet
    route_loops(lay, sses, bb)
    m = contacts(lay, sses)
    assert 0.0 < m["ss_contacts_nearer"] <= 1.0
    assert "ss_contacts_nearer" not in THRESHOLDS and not judge({"ss_contacts_nearer": 0.0})


def test_contact_sheet_panels_hang_from_the_top(monkeypatch, tmp_path):
    """1LMB: the 3D panel of a straight duplex is wider than the 2.5:1 the sheet allows, so it is shrunk in
    height; it hangs from the top of its slot (anchor N), its title in line with the others, not mid-sheet."""
    import matplotlib.pyplot as plt

    nv = _na_validate()
    path = Path(__file__).parent / "data" / "1LMB.cif"
    lay, sses, bb = make_layout(path)
    fig = draw(lay, route_loops(lay, sses, bb), sses, "1LMB")
    made, figure = [], plt.figure
    monkeypatch.setattr(plt, "figure", lambda *a, **k: made.append(figure(*a, **k)) or made[-1])
    monkeypatch.setattr(plt, "close", lambda *a: None)
    nv.contact_sheet("1LMB", lay, sses, bb, fig, {"failed": []}, tmp_path)
    sheet = made[-1]
    sheet.canvas.draw()
    tops = [a.get_position().y1 for a in sheet.axes]
    assert max(tops) - min(tops) < 1e-6, tops
