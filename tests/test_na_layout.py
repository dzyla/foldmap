"""Nucleic acids on the page: stems placed like elements, strands routed, four drawing styles, the explorer."""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap.cli import make_figure_and_loops, make_layout
from foldmap.layout import PlacedStem, label_boxes
from foldmap.route import route_loops

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB"]
SMALL = ["30QM", "1EHZ", "1IHF", "3CRX"]


def _path(pdb):
    return NA / f"{pdb}.cif.gz"


def _routed(pdb, **kw):
    lay, sses, bb = make_layout(_path(pdb), **kw)
    loops = route_loops(lay, sses, bb)
    return lay, sses, bb, loops


def _overlap(a, b):
    return min(a[2], b[2]) - max(a[0], b[0]) > 1e-6 and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6


def _gid(fig, prefix):
    ax = fig.axes[0]
    return [a for a in [*ax.patches, *ax.texts, *ax.lines, *ax.collections] if (a.get_gid() or "").startswith(prefix)]


@pytest.mark.parametrize("pdb", FIXTURES)
def test_every_stem_piece_is_placed_without_overlap(pdb):
    from foldmap.layout import joined_pieces

    lay, _, _ = make_layout(_path(pdb))
    pieces = len(lay.stems)
    assert pieces >= len(lay.nucleic.structure.stems) > 0
    joined = joined_pieces(lay.stems, lay.nucleic)  # consecutive pieces of a bent duplex meet at their joint
    boxes = [(p.sse.id, p.rect) for p in lay.placed.values()] + [(d.id, d.rect) for d in lay.stems]
    assert not any(
        _overlap(a, b) and frozenset((i, j)) not in joined for k, (i, a) in enumerate(boxes) for j, b in boxes[k + 1 :]
    )
    boxes = [b for _, b in boxes]
    x0, y0, x1, y1 = lay.bounds
    assert all(x0 <= r[0] and r[2] <= x1 and y0 <= r[1] and r[3] <= y1 for r in boxes)


def test_stems_are_drawn_to_scale_along_their_projected_axes():
    lay, _, _ = make_layout(_path("30QM"))
    d = next(d for d in lay.stems if d.name == "P1")
    assert isinstance(d, PlacedStem)
    assert d.width == pytest.approx(20.0 * 0.23)
    firsts = sorted((key for key in d.axial if d.side[key] == d.first), key=lambda k: k[1])
    assert np.all(np.diff([d.axial[k] for k in firsts]) > 0)  # the first arm runs 5'->3' along the direction
    seconds = sorted((key for key in d.axial if d.side[key] != d.first), key=lambda k: k[1])
    assert np.all(np.diff([d.axial[k] for k in seconds]) < 0)  # the second arm back
    for a, i, b, j in d.pairs:
        assert d.axial[(a, i)] == pytest.approx(d.axial[(b, j)])  # a rung straight across


def _joint_steps(lay):
    """(piece, next piece, exit port, entry port) where a strand goes from one piece of a bent duplex to the next
    (layout.joined_pieces)."""
    from foldmap.layout import joined_pieces

    joined = joined_pieces(lay.stems, lay.nucleic)
    at = {key: d for d in lay.stems for key in d.axial}
    out = []
    for s, strand in enumerate(lay.nucleic.strands):
        for i in range(len(strand) - 1):
            d, e = at.get((s, i)), at.get((s, i + 1))
            if d is not None and e is not None and frozenset((d.id, e.id)) in joined:
                out.append((d, e, d.port((s, i), True)[0], e.port((s, i + 1), False)[0]))
    return out


def _joints(lay):
    """(exit port, entry port, piece width) where a strand goes from one piece of a bent duplex to the next."""
    return [(p, q, d.width) for d, _, p, q in _joint_steps(lay)]


@pytest.mark.parametrize("pdb", ["1KX5", "1IHF"])
def test_a_bent_duplex_is_a_chain_of_pieces_end_to_end(pdb):
    """Consecutive pieces stay near each other. Box separation parts tilted pieces, so a joint can be a couple of
    duplex widths long; Task 7 routes joints as loops between stem ends (route_nucleic, _clear). With no margin
    kept between consecutive pieces (C3) this bounds the current lengths (1KX5 at most 2.34 widths, its 22 pieces
    separated one by one; 1IHF 2.27; were 2.47 and 2.42) against regressions."""
    lay, _, _ = make_layout(_path(pdb))
    joints = _joints(lay)
    assert len(joints) >= (40 if pdb == "1KX5" else 4)  # both strands at every joint
    gaps = [float(np.hypot(*(q - p))) / w for p, q, w in joints]
    assert max(gaps) <= 2.4, sorted(gaps)[-5:]


def _crosses_box(p, q, box) -> bool:
    """Does the segment p-q pass through the inside of box (x0, y0, x1, y1)?"""
    from foldmap.route import clip_len

    return clip_len(np.asarray(p, float), np.asarray(q, float), box) > 1e-9


def test_few_elements_are_placed_between_consecutive_pieces_of_a_bent_duplex():
    """1IHF: the joint springs held the pieces' ends together but the overlap term pushed their boxes a margin
    apart, so six of the protein's elements came between DNA 1b and DNA 1c (on a joint, or on the line between the
    middles of its two strands' steps). Consecutive pieces of one duplex keep no margin between them and are
    separated together (chains of at most four pieces); two elements remain, a regression bound: IHF's arms wrap
    the DNA in 3D, so in projection they lie over it between its pieces (review.md; making the joints firmer
    collapses the nucleosome ring)."""
    lay, _, _ = make_layout(_path("1IHF"))
    steps: dict[tuple[str, str], list] = {}
    for d, e, p, q in _joint_steps(lay):
        steps.setdefault(tuple(sorted((d.id, e.id))), []).append((d, e, p, q))
    assert len(steps) == 2 and all(len(v) == 2 for v in steps.values())
    between = []
    for (d, e, p0, q0), (_, _, p1, q1) in steps.values():
        lines = [(p0, q0), (p1, q1), ((p0 + q1) / 2, (q0 + p1) / 2)]
        others = [(p.sse.id, p.rect) for p in lay.placed.values()] + [
            (x.name, x.rect) for x in lay.stems if x not in (d, e)
        ]
        between += [name for name, box in others if any(_crosses_box(a, b, box) for a, b in lines)]
    assert len(between) <= 2, between  # A:50-57 and A:76-83 (six before)


@pytest.mark.parametrize("pdb", ["1KX5", "1IHF", "1RVA", "4YBB"])
def test_consecutive_pieces_overlap_by_at_most_a_quarter_width(pdb):
    """Consecutive pieces of a duplex may touch, but their bodies overlap by at most a quarter of a duplex width
    (layout._JOINT_OVERLAP): with no overlap term between them at all, 1RVA's two pieces were drawn half on top of
    each other at the kink."""
    from foldmap.layout import body_overlap, joined_pieces

    lay, _, _ = make_layout(_path(pdb))
    stems = {d.id: d for d in lay.stems}
    pairs = sorted(tuple(sorted(p)) for p in joined_pieces(lay.stems, lay.nucleic))
    assert pairs
    for a, b in pairs:
        d, e = stems[a], stems[b]
        assert body_overlap(d, e) <= 0.25 * d.width, (d.name, e.name, body_overlap(d, e))


def test_the_nucleosome_dna_wraps_the_histones():
    """The 22 pieces of the nucleosome DNA moved as one group carried 12 histone elements out of the ring, whose
    loops then crossed the DNA (54 crossings). Chains of more than four pieces are separated piece by piece: every
    element lies inside the DNA's hull, and at most two protein loops cross a strand."""
    from foldmap.nacheck import na_crossings
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(_path("1KX5"))
    loops = route_loops(lay, sses, bb)
    corners = np.array([(r[x], r[y]) for d in lay.stems for r in [d.rect] for x in (0, 2) for y in (1, 3)])
    hull = _convex_hull(corners)
    outside = [
        p.sse.id
        for p in lay.placed.values()
        if not _in_hull(hull, ((p.rect[0] + p.rect[2]) / 2, (p.rect[1] + p.rect[3]) / 2))
    ]
    assert not outside, outside
    assert na_crossings(lay, loops)["na_protein_crossings"] <= 2


def _convex_hull(pts):
    pts = sorted(set(map(tuple, pts.tolist())))

    def half(seq):
        out = []
        for p in seq:
            while len(out) >= 2 and _cross2(np.subtract(out[-1], out[-2]), np.subtract(p, out[-2])) <= 0:
                out.pop()
            out.append(p)
        return out[:-1]

    return half(pts) + half(pts[::-1])


def _in_hull(hull, p) -> bool:
    return all(_cross2(np.subtract(b, a), np.subtract(p, a)) >= 0 for a, b in zip(hull, hull[1:] + hull[:1]))


def _cross2(u, v) -> float:
    return float(u[0] * v[1] - u[1] * v[0])


def _joined_both_ways(lay):
    """(stem A, stem B, (s, i, j), (t, k, l)): a strand steps from A's nucleotide i to B's j, another from B's k
    back to A's l (the next piece of a bent duplex, or the stem across an internal loop or bulge)."""
    at = {key: n for n, d in enumerate(lay.stems) for key in d.axial}
    steps = {}
    for s, strand in enumerate(lay.nucleic.strands):
        held = [i for i in range(len(strand)) if (s, i) in at]
        for i, j in zip(held, held[1:]):
            if at[(s, i)] != at[(s, j)]:
                steps.setdefault((at[(s, i)], at[(s, j)]), []).append((s, i, j))
    return [
        (lay.stems[a], lay.stems[b], there[0], steps[(b, a)][0])
        for (a, b), there in sorted(steps.items())
        if a < b and len(there) == 1 and len(steps.get((b, a), ())) == 1
    ]


def _chords_cross(a, b, there, back) -> bool:
    """Do the straight chords (rail point to rail point) of the two strands joining stems a and b cross?"""
    (s, i, j), (t, k, l) = there, back
    return _segments_cross(a.point((s, i)), b.point((s, j)), b.point((t, k)), a.point((t, l)))


@pytest.mark.parametrize("pdb", FIXTURES)
def test_strands_cross_between_joined_stems_only_where_one_is_seen_clearly(pdb):
    """A stem whose arms the projection cannot tell apart (margin < _SIDE_CLEAR) takes the side that keeps the
    strands to its neighbours uncrossed, so two stems joined both ways (a bent duplex's next piece, the stem
    across an internal loop) are drawn with crossed strands where at least one of them is seen clearly (what the
    projection shows: the strands of a helix cross every half turn), or, when a run of unclear stems lies between
    clear ones seen on opposite sides, once, at the run's least clearly seen link."""
    from foldmap.layout import _SIDE_CLEAR

    lay, _, _ = make_layout(_path(pdb))
    links = _joined_both_ways(lay)
    unclear = {d.name for d in lay.stems if d.margin < _SIDE_CLEAR}
    group = {n: {n} for n in unclear}  # unclear stems linked to one another
    for a, b, _, _ in links:
        if a.name in unclear and b.name in unclear and group[a.name] is not group[b.name]:
            merged = group[a.name] | group[b.name]
            for n in merged:
                group[n] = merged
    crossed: dict[frozenset, list] = {}
    for a, b, there, back in links:
        if _chords_cross(a, b, there, back) and max(a.margin, b.margin) < _SIDE_CLEAR:
            crossed.setdefault(frozenset(group[a.name]), []).append((a, b))
    for members, pairs in crossed.items():
        assert len(pairs) == 1, [(a.name, b.name) for a, b in pairs]
        (a, b), mine = pairs[0], [(c, d) for c, d, _, _ in links if c.name in members or d.name in members]
        assert max(a.margin, b.margin) == min(max(c.margin, d.margin) for c, d in mine), (a.name, b.name)


def test_nucleosome_strands_cross_only_where_its_pieces_are_seen_to():
    """At 15 of the 21 joints of the nucleosome DNA the strands crossed when every piece took the side it is seen
    on. 5 crossings join two pieces seen clearly on opposite sides (1e-1f, 1m-1n, 1n-1o, 1o-1p, 1s-1t) and the
    unclear pieces 1g-1l lie between 1f and 1m, which are seen on opposite sides, so 6 is the least without
    drawing a clearly seen piece mirrored; the sixth falls between the least clear pieces, 1g and 1h."""
    lay, _, _ = make_layout(_path("1KX5"))
    joints = [j for j in _joined_both_ways(lay) if j[0].parent and j[0].parent == j[1].parent]
    assert len(joints) == 21
    assert sum(_chords_cross(*j) for j in joints) <= 6
    assert [
        (a.name, b.name) for a, b, *rest in joints if _chords_cross(a, b, *rest) and a.margin < 2.75 and b.margin < 2.75
    ] == [("DNA 1g", "DNA 1h")]


@pytest.mark.parametrize("pdb", FIXTURES)
def test_stems_seen_clearly_keep_their_side_and_none_is_drawn_mirrored(pdb):
    """Only a stem whose side the projection cannot tell is turned round; the rotation-vs-reflection fit of
    every drawn stem to its projected C1' atoms never prefers the reflection (by 30 %)."""
    from foldmap.layout import _SIDE_CLEAR, SCALE, _side_signal, stem_pieces
    from foldmap.nucleic import stem_axis

    def fit(P, Q, reflect):
        P, Q = P - P.mean(0), Q - Q.mean(0)
        P = P * [1, -1] if reflect else P
        U, _, Vt = np.linalg.svd(P.T @ Q)
        if np.linalg.det(U @ Vt) < 0:
            U[:, -1] *= -1
        return np.sqrt(((P @ (U @ Vt) - Q) ** 2).sum(1).mean())

    lay, _, _ = make_layout(_path(pdb))
    na, frame = lay.nucleic, lay.frame
    pieces = {p.id: p for p in stem_pieces(na)}
    for d in lay.stems:
        signal = _side_signal(pieces[d.id], na, frame, stem_axis(pieces[d.id], na.strands)[0], d.angle)
        assert d.margin == pytest.approx(abs(signal))
        if d.margin >= _SIDE_CLEAR:
            assert d.first == (1 if signal >= 0 else -1), d.name
        keys = sorted(d.axial)
        P = np.asarray(frame.project(np.array([na.strands[s].c1[i] for s, i in keys])))[:, :2] * SCALE
        Q = np.array([d.point(k) for k in keys])
        assert fit(P, Q, True) >= 0.7 * fit(P, Q, False), d.name


def test_internal_loop_sides_cross_less_on_the_ribosome():
    """The two sides of an internal loop are bowed apart; their chords crossed at 7 of 11 such loops on 4YBB when
    every stem took the side it is seen on. A regression bound: the crossings left join stems seen clearly. (10
    loops since AA:454-462/470-479's 2x2 loop in register stays inside its stem.)"""
    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    kind = {(r.strand, r.lo): r.kind for r in lay.nucleic.structure.runs}
    sides = [p for p in lay.na_paths if p.kind == "run" and p.nts and kind.get((p.strand, p.nts[0])) == "internal"]
    loops = [
        (p, q)
        for n, p in enumerate(sides)
        for q in sides[n + 1 :]
        if {at[k].id for k in p.between} == {at[k].id for k in q.between}
    ]
    assert len(loops) == 10
    assert sum(_segments_cross(p.points[0], p.points[-1], q.points[0], q.points[-1]) for p, q in loops) <= 6


def test_short_runs_and_strand_exchanges_hold_their_stem_ends_together(monkeypatch):
    """A strand passing from one stem to another of more than two strands through at most two nucleotides (the
    four 0-nt strand exchanges of the Holliday junction) pulls the two stem ends together more firmly than a bent
    duplex's joint, rest length its nucleotides' steps: DNA 2 touches no protein and still comes to the junction.
    Not within one RNA chain (a tRNA's short junction runs)."""
    import foldmap.layout as layout

    seen = []
    real = layout._attract

    def spy(*a, **k):
        seen.append(list(k["extra"]))
        return real(*a, **k)

    monkeypatch.setattr(layout, "_attract", spy)
    make_layout(_path("3CRX"))
    make_layout(_path("1EHZ"))
    springs, trna = seen
    assert len(springs) == 4
    assert all(w == layout._W_SHORT > layout._W_JOINT and slack == layout._JOINT_SLACK for *_, w, slack in springs)
    assert trna and all(w < layout._W_SHORT for *_, w, _ in trna)


def _seen(lay):
    """Each drawn stem's centre as seen in the layout's frame, page units."""
    from foldmap.layout import SCALE, stem_pieces
    from foldmap.nucleic import stem_axis

    pieces = {p.id: p for p in stem_pieces(lay.nucleic)}
    centres = np.array([stem_axis(pieces[d.id], lay.nucleic.strands)[0] for d in lay.stems])
    return {d.name: np.asarray(c)[:2] * SCALE for d, c in zip(lay.stems, lay.frame.project(centres))}


def test_trna_variable_loop_spans_what_its_l_shape_shows():
    """1EHZ: the variable loop (5 nt from the anticodon stem P3 to the T stem P4) is a short run, and the T loop
    does not point at the anticodon loop (with every short run held like a strand exchange, P4 was pulled across
    the cloverleaf: the variable loop ran 20.7 units and the T loop pointed 32 degrees off the anticodon loop).
    The tRNA is drawn in its L shape (runs pull only beyond the distance 3D shows), so the variable loop spans what
    3D shows between P3 and P4 rather than a cloverleaf's 6 units: its ends lie no farther apart than 1.35 times
    their distance at the stems' projected places (measured 1.27), and the path drawn round the stems is at most
    2.1 times that distance long (measured 1.98)."""
    import math

    lay, _, _, _ = _routed("1EHZ")
    at = {key: d for d in lay.stems for key in d.axial}
    variable = [p for p in lay.na_paths if p.kind == "run" and [at[k].name for k in p.between] == ["P3", "P4"]]
    assert len(variable) == 1 and len(variable[0].nts) == 5
    (a, b), seen = variable[0].between, _seen(lay)
    da, db = at[a], at[b]
    drawn = math.dist(da.port(a, True)[0], db.port(b, False)[0])
    there = math.dist(
        da.port(a, True)[0] - (da.cx, da.cy) + seen["P3"], db.port(b, False)[0] - (db.cx, db.cy) + seen["P4"]
    )
    assert drawn <= 1.35 * there, (drawn, there)
    path = sum(math.dist(p, q) for p, q in zip(variable[0].points, variable[0].points[1:]))
    assert path <= 2.1 * there, (path, there)

    def loop(d):
        a, i, b, j = d.pairs[-1]
        centre, _ = d.hairpin_circle((a, i), (b, j), d.hairpin)
        return centre, d.direction if d.axial[(a, i)] >= 0 else -d.direction

    stems = {d.name: d for d in lay.stems}
    (c4, out4), (c3, _) = loop(stems["P4"]), loop(stems["P3"])
    towards = (c3 - c4) / np.linalg.norm(c3 - c4)
    assert out4 @ towards < np.cos(np.radians(60))


def test_a_trna_is_drawn_in_its_3d_arrangement():
    """A run between two stems pulls only when the drawing stretches it beyond the distance 3D shows: pulled to its
    nucleotides' length, the runs drew 1EHZ's cloverleaf away from its L shape (procrustes 0.716)."""
    from foldmap.nacheck import shape

    lay, _, _ = make_layout(_path("1EHZ"))
    assert shape(lay)["procrustes"] < 0.4


def test_stems_joined_by_runs_keep_the_order_they_are_seen_in():
    """2OEU (hammerhead ribozyme): P3 is seen left of P2; the runs pulling P2 and P3 to their nucleotides' length
    swapped them on the page."""
    lay, _, _ = make_layout(_path("2OEU"))
    seen = _seen(lay)
    drawn = {d.name: np.array([d.cx, d.cy]) for d in lay.stems}
    assert (drawn["P3"] - drawn["P2"]) @ (seen["P3"] - seen["P2"]) > 0


@pytest.mark.parametrize("pdb", ["1RVA", "4LD0"])
def test_nucleic_acid_torn_by_unrolling_is_drawn_in_the_plain_frame(pdb):
    """A symmetric dimer is unrolled one protomer beside the other. When that would tear a nucleic acid held by
    both protomers apart (1RVA: the two pieces of its kinked duplex land on either side of the seam, drawn as a
    180° hairpin; 4LD0: the junction's arms turn against one another by up to 139°), the plain frame is used, so each
    stem lands next to what touches it."""
    from foldmap.frame import Frame
    from foldmap.nacheck import contacts, shape

    lay, sses, _ = make_layout(_path(pdb))
    assert type(lay.frame) is Frame
    m = {**shape(lay), **contacts(lay, sses)}
    if pdb == "1RVA":  # a kink, drawn as a 180° hairpin when unrolled; now drawn as the frame sees it (the view
        # sees the duplex nearly end-on, so its 39° kink projects to about 102°: the frame's choice, not the drawing's)
        assert abs(m["bend_drawn"] - m["bend_3d"]) <= 10
    if pdb == "4LD0":  # each arm was drawn under the protomer that does not touch it (AUC 0.506)
        assert m["contact_auc"] >= 0.6


def test_a_nucleic_acid_the_unrolling_does_not_tear_keeps_the_side_by_side_frame():
    """30QM: one strand is held by both protomers, but unrolling tears none of its duplexes (each stem is drawn
    whole where its centre is seen) and turns no linked stems against one another by more than 57°, so the two
    protomers stay side by side, the figure the author chose (the plain frame interleaved them into one tangle)."""
    from foldmap.frame import CylinderFrame

    lay, _, _ = make_layout(_path("30QM"))
    assert type(lay.frame) is CylinderFrame


def test_a_kinked_duplex_across_a_dimer_is_no_longer_a_hairpin():
    """1RVA: unrolled, the two half-sites were drawn antiparallel side by side (bend 179.7°); projected, the drawn
    bend is the bend seen in the frame."""
    from foldmap.frame import Frame
    from foldmap.nacheck import shape

    lay, _, _ = make_layout(_path("1RVA"))
    assert type(lay.frame) is Frame
    m = shape(lay)
    assert abs(m["bend_drawn"]) < 150 and m["joint_error_max"] <= 5.0


def test_unrolled_frame_kept_when_each_strand_is_held_by_one_protomer():
    """The plain frame only for a nucleic acid with stems that some strand ties to two protomers, each touching
    it with at least NA_HOLD_MIN (3) residues: a glancing contact from the other protomer does not count."""
    from types import SimpleNamespace

    from foldmap.cli import _na_across_protomers
    from foldmap.model import ResLabel

    bb = SimpleNamespace(labels=[ResLabel(c, k, "", "ALA") for k, c in enumerate("AAABBBAB")])
    sym = SimpleNamespace(protomers=[["A"], ["B"]])
    stems = SimpleNamespace(stems=[object()])
    one_each = SimpleNamespace(structure=stems, contacts={r: {(0 if r in (0, 1, 2, 6) else 1, r)} for r in range(8)})
    shared = SimpleNamespace(structure=stems, contacts={r: {(0, r)} for r in range(6)})
    glancing = SimpleNamespace(structure=stems, contacts={0: {(0, 1)}, 1: {(0, 2)}, 2: {(0, 3)}, 3: {(0, 9)}})
    two_two = SimpleNamespace(structure=stems, contacts={0: {(0, 1)}, 1: {(0, 2)}, 3: {(0, 3)}, 4: {(0, 9)}})
    unpaired = SimpleNamespace(structure=SimpleNamespace(stems=[]), contacts=shared.contacts)
    assert not _na_across_protomers(one_each, sym, bb)
    assert _na_across_protomers(shared, sym, bb)
    assert not _na_across_protomers(glancing, sym, bb)
    assert not _na_across_protomers(two_two, sym, bb)
    assert not _na_across_protomers(unpaired, sym, bb)
    assert not _na_across_protomers(shared, None, bb)
    assert not _na_across_protomers(None, sym, bb)


def test_holliday_junction_reads_as_four_arms():
    """3CRX: the arms' inner ends gather at the junction (DNA 2, touching no protein, was 2.8 widths away), the
    arms lie in four quadrants around it (DNA 2 shared one with DNA 1), and the strand exchanges are short.
    Box separation of the tilted arms keeps the inner ends 0.9-1.5 widths from the centre. Measured on the
    deposited half (two Cre): in the four-Cre assembly two Cre helices sit between the inner ends (placement,
    plan 07 Task 13b step 2)."""
    import math

    lay, _, _, _ = _routed("3CRX", assembly="asu")
    at = {key: d for d in lay.stems for key in d.axial}
    w = lay.stems[0].width
    ends = {d.name: [] for d in lay.stems}
    for s, strand in enumerate(lay.nucleic.strands):
        held = [i for i in range(len(strand)) if (s, i) in at]
        for i, j in zip(held, held[1:]):
            if at[(s, i)] is not at[(s, j)] and j - i - 1 <= 2:
                ends[at[(s, i)].name].append(at[(s, i)].point((s, i)))
                ends[at[(s, j)].name].append(at[(s, j)].point((s, j)))
    inner = {n: np.mean(v, axis=0) for n, v in ends.items()}
    centre = np.mean(list(inner.values()), axis=0)
    assert len(inner) == 4
    assert all(math.dist(p, centre) <= 1.6 * w for p in inner.values()), inner
    assert len({(d.cx >= centre[0], d.cy >= centre[1]) for d in lay.stems}) == 4
    paths = sorted(
        sum(math.dist(a, b) for a, b in zip(p.points, p.points[1:])) / w
        for p in lay.na_paths
        if p.kind == "run" and not p.nts
    )
    assert len(paths) == 4 and paths[2] <= 1.7 and paths[3] <= 4.3, paths


@pytest.mark.parametrize("pdb", ["30QM", "1EHZ", "4OO8", "4YBB"])
def test_hairpin_loops_lie_within_their_stems_footprint(pdb):
    """The room a hairpin loop needs beyond the stem's tip, its letters included, is part of the stem's
    footprint, so no neighbour is placed where the loop is drawn."""
    from foldmap.layout import _GLYPH_H, _footprint

    lay, _, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    arcs = [p for p in lay.na_paths if p.kind == "hairpin"]
    assert arcs
    for p in arcs:
        d = at[p.between[0]]
        x0, y0, x1, y1 = _footprint(d)
        assert all(x0 - 1e-6 <= x - d.cx <= x1 + 1e-6 and y0 - 1e-6 <= y - d.cy <= y1 + 1e-6 for x, y in p.points), (
            d.name
        )
        g = _GLYPH_H / 2
        for i in p.nts:
            x, y = lay.na_letters[(p.strand, i)]
            assert x0 <= x - d.cx - g and x - d.cx + g <= x1 and y0 <= y - d.cy - g and y - d.cy + g <= y1, (d.name, i)


def test_a_long_tie_is_two_labelled_marks():
    """1EHZ's G19·C56 (D loop to T loop) would cross two ladders as a line: it is drawn as a tick at each
    nucleotide, both labelled "a"."""
    fig, _ = make_figure_and_loops(_path("1EHZ"))
    marks, labels = _gid(fig, "na:tie:19-56:mark"), _gid(fig, "na:tie:19-56:label")
    assert len(marks) == 2 and [t.get_text() for t in labels] == ["a", "a"]
    for t in labels:  # not a base: smaller than the letters, upright, in the rung grey
        base = _gid(fig, "na:nt:")[0]
        assert t.get_fontsize() == pytest.approx(0.75 * base.get_fontsize())
        assert t.get_fontstyle() == "normal" and base.get_fontstyle() == "normal"
        assert to_hex(t.get_color()) != to_hex(base.get_color())
    assert not [a for a in _gid(fig, "na:tie:") if not a.get_gid().endswith((":mark", ":label"))]


@pytest.mark.parametrize("pdb", FIXTURES)
def test_tie_lines_cross_no_stem_and_tie_labels_count_as_text(pdb):
    from foldmap.layout import label_boxes
    from foldmap.route_na import crosses_body

    lay, sses, _, _ = _routed(pdb)
    names = [n for n, _ in label_boxes(lay, sses)]
    for n, tie in enumerate(lay.na_ties):
        if tie.label:
            assert len(tie.marks) == 2 and {f"tie-label:{n}:0", f"tie-label:{n}:1"} <= set(names)
        else:
            assert not any(crosses_body(*tie.ends, e) for e in [*lay.stems, *lay.placed.values()])


@pytest.mark.parametrize("pdb", FIXTURES)
def test_both_arms_run_monotonic_along_every_stem(pdb):
    lay, _, _ = make_layout(_path(pdb))
    for d in lay.stems:
        firsts = sorted((key for key in d.axial if d.side[key] == d.first), key=lambda k: k[1])
        seconds = sorted((key for key in d.axial if d.side[key] != d.first), key=lambda k: k[1])
        assert np.all(np.diff([d.axial[k] for k in firsts]) > 0), d.name  # bulges included
        assert np.all(np.diff([d.axial[k] for k in seconds]) < 0), d.name


@pytest.mark.parametrize("pdb", FIXTURES)
def test_bulged_nucleotides_get_room_on_their_rail(pdb):
    from foldmap.layout import _MIN_RISE, SCALE

    lay, _, _ = make_layout(_path(pdb))
    for d in lay.stems:
        for sign in (1, -1):
            along = sorted(v for key, v in d.axial.items() if d.side[key] == sign)
            assert np.all(np.diff(along) >= 0.9 * _MIN_RISE * SCALE), (d.name, np.round(np.diff(along), 3))


def _upright_size(t, renderer):
    """Width and height (display units) of a text as if written level."""
    rot = t.get_rotation()
    t.set_rotation(0)
    box = t.get_window_extent(renderer)
    t.set_rotation(rot)
    return box.width, box.height


@pytest.mark.parametrize("pdb", ["1KX5", "30QM", "3CRX"])
def test_stem_names_keep_clear_of_base_letters(pdb):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig, _ = make_figure_and_loops(_path(pdb))
    renderer = FigureCanvasAgg(fig).get_renderer()
    to_px = fig.axes[0].transData.transform
    names, letters = _gid(fig, "na:label:"), _gid(fig, "na:nt:")
    assert names and letters
    for n in names:  # each letter (as a disc) against the name's own box, written along its stem
        w, h = _upright_size(n, renderer)
        c, rot = to_px(n.get_position()), np.radians(n.get_rotation())
        along, across = np.array([np.cos(rot), np.sin(rot)]), np.array([-np.sin(rot), np.cos(rot)])
        for t in letters:
            r = max(_upright_size(t, renderer)) * 1.5  # a letter and a half from its centre: read apart at a glance
            q = to_px(t.get_position()) - c
            assert abs(q @ along) >= w / 2 + r or abs(q @ across) >= h / 2 + r, (n.get_text(), t.get_gid())


def test_layout_is_the_same_whatever_the_hash_seed():
    code = (
        "from foldmap.cli import make_layout\n"
        "from foldmap.route import route_loops\n"
        f"lay, sses, bb = make_layout({str(_path('30QM'))!r})\n"
        "print([(d.id, round(d.cx, 6), round(d.cy, 6), round(d.angle, 6)) for d in lay.stems])\n"
        "route_loops(lay, sses, bb)\n"
        "print([(p.kind, p.points, p.nt_xy) for p in lay.na_paths], lay.na_ends, sorted(lay.na_letters.items()))"
    )
    runs = {
        subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONHASHSEED": seed, "OPENBLAS_NUM_THREADS": "1"},
            check=True,
        ).stdout
        for seed in ("0", "1", "2")
    }
    assert len(runs) == 1


def test_wobble_pairs_are_dotted_rungs():
    fig, _ = make_figure_and_loops(_path("30QM"))
    rungs = [a for a in _gid(fig, "na:stem:") if a.get_gid().endswith(":rung")]
    assert sum(1 for r in rungs if r.get_linestyle() != "-") == 4  # 30QM has four G·U pairs


def test_beads_mark_every_touched_nucleotide_in_its_chain_colour():
    from matplotlib.colors import to_hex

    from foldmap.palette import chain_colors

    fig, lay = make_figure_and_loops(_path("30QM"))[0], make_layout(_path("30QM"))[0]
    beads = _gid(fig, "na:contact:")
    colours = chain_colors(["A", "B"])
    for b in beads:
        chain = b.get_gid().split(":")[2]
        assert {to_hex(c) for c in b.get_facecolors()} == {colours[chain].lower()}
    touched = {x for v in lay.nucleic.contacts.values() for x in v}
    assert sum(len(b.get_offsets()) for b in beads) == len(touched)


def test_protein_free_rna_is_drawn_on_its_own():
    lay, sses, bb, loops = _routed("1EHZ")
    assert len(bb) == 0 and sses == [] and loops == []
    assert [d.name for d in lay.stems] == ["P1", "P2", "P3", "P4"]
    assert lay.na_paths


@pytest.mark.parametrize("pdb", FIXTURES)
def test_strands_route_without_fallbacks_and_every_nucleotide_is_drawn_once(pdb):
    lay, _, _, _ = _routed(pdb)
    assert not [p for p in lay.na_paths if p.fallback]
    na = lay.nucleic
    drawn = [key for d in lay.stems for key in d.axial] + [(p.strand, i) for p in lay.na_paths for i in p.nts]
    assert len(drawn) == len(set(drawn))
    assert set(drawn) == {(s, i) for s, strand in enumerate(na.strands) for i in range(len(strand))}
    for p in lay.na_paths:
        assert len(p.nt_xy) == len(p.nts)


def test_runs_go_round_the_other_elements():
    from test_route import segments_clear

    lay, _, _, _ = _routed("30QM")
    for p in lay.na_paths:
        if p.kind == "run":
            others = [q.rect for q in lay.placed.values()]
            assert segments_clear(p.points[1:-1], others), p.nts[:3]


@pytest.mark.parametrize("pdb", ["1KX5", "1IHF", "3CRX", "30QM", "4OO8"])
def test_joints_and_runs_cross_no_stem_or_element(pdb):
    """Every run and every joint between two pieces of a bent duplex, as drawn from where it leaves its own stems'
    bodies, passes the router's own final check (route._clear) against every element and stem. A joint drawn
    straight from face to face (two pieces that touch or overlap, which no overlap term keeps apart) is checked as a
    kink: clear by _INFLATE of every element and stem but its own two."""
    from foldmap.route import _INFLATE, _clear, _element_obstacles
    from foldmap.route_na import among

    lay, _, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    rects, owners = _element_obstacles(lay)
    kinds = {"run": 0, "joint": 0}
    for p in lay.na_paths:
        if p.kind in kinds:
            own = tuple(at[key] for key in p.between)
            if p.kind == "joint" and len(p.points) == 4:
                assert _clear(p.points[1:3], own, rects, owners, grow=_INFLATE, kink=True), (p.kind, p.between)
            else:
                assert _clear(among(p, at), own, rects, owners, own_body=True), (p.kind, p.between)
            kinds[p.kind] += 1
    assert kinds["joint"] >= (20 if pdb == "1KX5" else 2 if pdb == "1IHF" else 0)
    assert kinds["run"] or pdb == "1KX5"


# runs over a stem-name box, held at the measured count: on 4YBB the 1-nt run P1->P2 ((0, 23) to (0, 25)) is boxed in
# by the stems and hairpin loops drawn round it (Task 13b step 4, once AA:454-462/470-479 became one stem and the
# helices were lettered): its one route clear of the elements gives up a soft cell of P12b's name box (block style
# draws that name inside the stem, so nothing visible is crossed)
OVER_NAMES = {"4YBB": [((0, 23), (0, 25))]}


@pytest.mark.parametrize("pdb", FIXTURES)
def test_runs_keep_off_the_stem_names(pdb):
    from foldmap.layout import label_boxes
    from foldmap.route_na import among

    lay, sses, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    names = [b for n, b in label_boxes(lay, sses) if n.startswith("stem-label:")]
    over = [
        p.between
        for p in lay.na_paths
        if p.kind in ("run", "joint")
        and not all(segment_clear(a, b, names) for a, b in zip(among(p, at), among(p, at)[1:]))
    ]
    assert over == OVER_NAMES.get(pdb, [])  # a held run that keeps off the names now: remove its hold


def test_a_run_goes_round_a_stem_name_it_could_cross():
    """4YBB: with the stems where 3D has them (C2), the 1-nt run from P27 to P26 (P17b to P17a since helices are
    lettered) was bowed across P27's name. Stem names are hard obstacles for a run's first drawing; one that crosses
    a name is drawn only when nothing else is clear of the elements."""
    from foldmap.layout import label_boxes
    from foldmap.route_na import among

    lay, sses, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    names = [b for n, b in label_boxes(lay, sses) if n.startswith("stem-label:")]
    run = next(p for p in lay.na_paths if p.between == ((1, 391), (1, 393)))
    assert [at[k].name for k in run.between] == ["P17b", "P17a"] and not run.fallback
    pts = among(run, at)
    assert all(segment_clear(a, b, names) for a, b in zip(pts, pts[1:]))


def test_a_run_over_a_stem_name_is_drawn_when_nothing_else_is_clear():
    """The names are hard only in the first pass: with every drawing over a name, each run is drawn as without."""
    from foldmap.route import _element_obstacles, _Grid, _Router
    from foldmap.route_na import _strand_between

    lay, _, _, _ = _routed("1EHZ")
    at = {key: d for d in lay.stems for key in d.axial}
    runs = [p for p in lay.na_paths if p.kind == "run"]
    assert runs

    def drawn(p, names):
        (s, i), (_, j) = p.between
        router = _Router(_Grid(lay), _element_obstacles(lay))
        return _strand_between(router, s, i, j, at[(s, i)], at[(s, j)], None, names)

    for p in runs:
        plain, covered = drawn(p, ()), drawn(p, [(-1e3, -1e3, 1e3, 1e3)])
        assert covered.points == plain.points and covered.fallback == plain.fallback


def test_an_internal_loop_side_keeps_its_bow_when_every_drawing_meets_a_name():
    """One side of an internal loop is bowed out on its own side. When the bow (and every route) crosses a stem
    name, the bow is still the drawing taken (the soft choice), not a route nor a fallback (4OO8's two sides)."""
    from foldmap.route import _element_obstacles, _Grid, _Router
    from foldmap.route_na import NT_STEP, _face, _side_bow, _strand_between

    lay, _, _, _ = _routed("4OO8")
    at = {key: d for d in lay.stems for key in d.axial}
    kind = {(r.strand, r.lo): r.kind for r in lay.nucleic.structure.runs}
    sides = [p for p in lay.na_paths if p.kind == "run" and p.nts and kind.get((p.strand, p.nts[0])) == "internal"]
    assert len(sides) == 2
    for p in sides:
        (s, i), (_, j) = p.between
        da, db = at[(s, i)], at[(s, j)]
        (ra, (fa, _)), (rb, (fb, _)) = (
            (da.point((s, i)), _face(da, (s, i), True)),
            (db.point((s, j)), _face(db, (s, j), False)),
        )
        reach = (j - i) * NT_STEP - np.hypot(*(fa - ra)) - np.hypot(*(rb - fb))
        bow = [tuple(map(float, ra)), *_side_bow(da, (s, i), fa, fb, reach), tuple(map(float, rb))]
        assert p.points == bow  # drawn as its bow, from face to face
        router = _Router(_Grid(lay), _element_obstacles(lay))
        covered = _strand_between(router, s, i, j, da, db, "internal", [(-1e3, -1e3, 1e3, 1e3)])
        assert covered.points == bow and not covered.fallback and covered.kind == "run"


def test_a_straight_joint_may_cross_its_partner_strand_but_no_other():
    """A joint no route reaches is drawn straight across the meeting corner of two touching pieces only if it
    crosses no strand drawn so far, except the other strand of the same joint (3D0U's P3b-P3a: the two strands
    cross there as they do along a helix)."""
    from foldmap.route import _element_obstacles, _Grid, _Router
    from foldmap.route_na import _face, _straight_joint, _Strands

    lay, _, _, _ = _routed("1RVA")
    at = {key: d for d in lay.stems for key in d.axial}
    joint = next(p for p in lay.na_paths if p.kind == "joint" and len(p.points) == 4)
    a, b = joint.between
    da, db = at[a], at[b]
    (fa, ea), (fb, _) = _face(da, a, True), _face(db, b, False)
    mid, normal = (np.add(fa, fb)) / 2, np.array([fa[1] - fb[1], fb[0] - fa[0]])
    across = [tuple(mid + normal), tuple(mid - normal)]  # a strand crossing the joint's middle
    router = _Router(_Grid(lay), _element_obstacles(lay))
    assert _straight_joint(router, da, db, fa, ea, fb, _Strands(router.grid))
    partner, other = _Strands(router.grid), _Strands(router.grid)
    partner.add(across, frozenset((da.id, db.id)))
    other.add(across)
    assert _straight_joint(router, da, db, fa, ea, fb, partner)
    assert not _straight_joint(router, da, db, fa, ea, fb, other)


def _enters_core(p, q, d) -> bool:
    """Does segment p-q run inside the stem's core, strictly between its two rails (where its pairs are drawn)?"""
    from foldmap.route import clip_len

    c, s = np.cos(d.angle), np.sin(d.angle)
    u0, v0 = (p[0] - d.cx) * c + (p[1] - d.cy) * s, (p[1] - d.cy) * c - (p[0] - d.cx) * s
    u1, v1 = (q[0] - d.cx) * c + (q[1] - d.cy) * s, (q[1] - d.cy) * c - (q[0] - d.cx) * s
    hl, hw = d.length / 2 - 1e-6, d.rail - 1e-6
    return clip_len((u0, v0), (u1, v1), (-hl, -hw, hl, hw)) > 1e-6


@pytest.mark.parametrize("pdb", ["2GIS", "1EHZ", "30QM", "4OO8", "1IHF", "1KX5", "3CRX"])
def test_no_run_or_joint_enters_its_own_stems_core(pdb):
    """A run or joint leaves its stem along its rail to the stem's end face and is drawn from face to face; no part
    of it runs between the rails of its own stems (2GIS: the sides of the P2-P3 internal loop were bowed from rail
    point to rail point, straight through P3's core, an X over both stems)."""
    lay, _, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    bad = [
        (p.kind, p.between, d.name)
        for p in lay.na_paths
        if p.kind in ("run", "joint")
        for d in {at[k].id: at[k] for k in p.between}.values()
        if any(_enters_core(a, b, d) for a, b in zip(p.points, p.points[1:]))
    ]
    assert not bad


def segment_clear(p, q, rects):
    from foldmap.route import clip_len

    return all(clip_len(p, q, r) <= 1e-9 for r in rects)


def test_hairpin_loops_are_arcs_beyond_the_stem_tip():
    lay, _, _, _ = _routed("30QM")
    by_key = {key: d for d in lay.stems for key in d.axial}
    hairpins = [p for p in lay.na_paths if p.kind == "hairpin"]
    assert len(hairpins) == 5
    for p in hairpins:
        d = by_key[(p.strand, p.nts[0] - 1)]
        tip = max(d.axial[(p.strand, p.nts[0] - 1)], d.axial[(p.strand, p.nts[-1] + 1)], key=abs)
        for xy in p.nt_xy:  # every loop nucleotide beyond the stem's end
            along = (np.asarray(xy) - [d.cx, d.cy]) @ d.direction
            assert along * np.sign(tip) > abs(tip) - 1e-6


def test_every_strand_end_is_labelled_unless_cropped_into_a_gap():
    lay, _, _, _ = _routed("4YBB")
    ends = {(s, e) for s, e, _ in lay.na_ends}
    gapped = {(g.before, "3′") for g in lay.nucleic.gaps} | {(g.after, "5′") for g in lay.nucleic.gaps}
    for s in range(len(lay.nucleic.strands)):
        for e in ("5′", "3′"):
            assert ((s, e) in ends) != ((s, e) in gapped)
    assert len([p for p in lay.na_paths if p.kind == "gap"]) == len(lay.nucleic.gaps)
    assert lay.na_cropped


def test_short_gaps_are_a_bare_break_and_long_ones_say_how_many_nucleotides():
    fig, _ = make_figure_and_loops(_path("4YBB"))
    lay, _, _, _ = _routed("4YBB")
    gaps = [p for p in lay.na_paths if p.kind == "gap"]
    counts = sorted(g.count for g in lay.nucleic.gaps)
    assert counts[0] < 3 <= counts[-1]  # 4YBB cropped: 1-nt gaps from keeping stems whole, and long ones
    for p, g in zip(gaps, lay.nucleic.gaps):
        assert p.label == ("" if g.count < 3 else f"… {g.count} nt …")
    labels = sorted(t.get_text() for t in _gid(fig, "na:gap:"))
    assert labels == sorted(f"… {c} nt …" for c in counts if c >= 3)
    assert len([a for a in _gid(fig, "na:run:") if a.get_linestyle() != "-"]) == len(gaps)  # every gap is a break


@pytest.mark.parametrize("pdb", ["30QM", "4OO8", "1EHZ"])
def test_no_text_leaves_the_page_and_letters_do_not_collide(pdb):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig, _ = make_figure_and_loops(_path(pdb))
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    r = canvas.get_renderer()
    page = fig.bbox
    texts = [t for t in fig.axes[0].texts if t.get_text().strip()]
    for t in texts:
        e = t.get_window_extent(r)
        assert e.x0 >= page.x0 - 1 and e.y0 >= page.y0 - 1 and e.x1 <= page.x1 + 1 and e.y1 <= page.y1 + 1, t.get_text()
    boxes = [t.get_window_extent(r) for t in texts if (t.get_gid() or "").startswith("na:nt:")]
    clash = 0
    for k, a in enumerate(boxes):
        for b in boxes[k + 1 :]:
            w, h = min(a.x1, b.x1) - max(a.x0, b.x0), min(a.y1, b.y1) - max(a.y0, b.y0)
            clash += w > 0 and h > 0 and w * h > min(a.width * a.height, b.width * b.height) / 3
    assert clash == 0


def test_strands_without_pairs_are_traced(monkeypatch):
    import foldmap.io

    monkeypatch.setattr(foldmap.io, "pair_nucleotides", lambda runs: [])  # as for a G-quadruplex
    fig, _ = make_figure_and_loops(_path("1IHF"))
    lay, _, _, _ = _routed("1IHF")
    assert lay.stems == [] and set(lay.na_traces) == {0, 1, 2}
    assert [p.kind for p in lay.na_paths] == ["trace"] * 3
    assert len(_gid(fig, "na:nt:")) == 70 and len(_gid(fig, "na:end:")) == 6
    assert any("DNA (C, D, E)" in t.get_text() for t in fig.axes[0].texts)


def test_plot_warns_about_strands_drawn_as_fallbacks(tmp_path, capsys, monkeypatch):
    import foldmap.route_na as route_na
    from foldmap.cli import main

    real = route_na.route_nucleic

    def one_stuck(layout, *a, **k):  # every fixture routes cleanly; mark one run as a fallback
        paths = real(layout, *a, **k)
        next(p for p in paths if p.kind == "run").fallback = True
        return paths

    monkeypatch.setattr(route_na, "route_nucleic", one_stuck)
    assert main(["plot", str(_path("1EHZ")), "-o", str(tmp_path / "t.svg")]) == 0
    err = capsys.readouterr().err
    lay, _, _, _ = _routed("1EHZ")
    run = next(p for p in lay.na_paths if p.kind == "run")
    labels = lay.nucleic.strands[run.strand].labels
    name = f"strand {labels[run.nts[0]].chain}:{labels[run.nts[0]].tag}-{labels[run.nts[-1]].tag}"
    assert "no clear route" in err and name in err, err


def test_app_and_spec_report_stuck_strands(monkeypatch):
    import foldmap.route_na as route_na
    from foldmap.spec import FigureSpec, draw_spec

    real = route_na.route_nucleic

    def all_stuck(layout, *a, **k):
        paths = real(layout, *a, **k)
        for p in paths:
            p.fallback = p.kind == "hairpin"
        return paths

    monkeypatch.setattr(route_na, "route_nucleic", all_stuck)
    drawn = draw_spec(FigureSpec(_path("1EHZ")))
    assert len([s for s in drawn.stuck if s.startswith("strand A:")]) == 3  # 1EHZ's three hairpin loops


def _segments_cross(p, q, a, b) -> bool:
    """Segments p-q and a-b cross properly (not just touch at an end)."""

    def turn(u, v, w):
        return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])

    return turn(a, b, p) * turn(a, b, q) < 0 and turn(p, q, a) * turn(p, q, b) < 0


@pytest.mark.parametrize("pdb", ["1EHZ", "4OO8", "30QM", "4LD0", "2GIS"])
def test_runs_steer_round_hairpin_loops_and_free_ends(pdb):
    """Hairpin loops and free ends are hard obstacles for a run's first routing pass (4LD0: the run C:24-25 crossed
    the free end D:13 when they were only a cost to cross)."""
    from foldmap.route_na import among

    lay, _, _, _ = _routed(pdb)
    at = {key: d for d in lay.stems for key in d.axial}
    drawn = {id(p): among(p, at) for p in lay.na_paths}
    fixed = [p for p in lay.na_paths if p.kind in ("hairpin", "tail")]
    for p in lay.na_paths:
        if p.kind in ("run", "joint", "gap"):
            a = drawn[id(p)]
            for f in fixed:
                b = drawn[id(f)]
                assert not any(_segments_cross(*s, *t) for s in zip(a, a[1:]) for t in zip(b, b[1:])), (
                    p.kind,
                    p.nts[:2],
                    f.kind,
                    f.nts[:2],
                )


def test_a_long_free_end_is_routed_round_what_blocks_it(monkeypatch):
    """A free end of more than six nucleotides that no straight drawing keeps clear of the elements is routed out to
    free space (5W1H's 17-nt spacer was drawn straight through three helices). 4YBB's 10-nt end AA:557-566, with
    only the straight drawing allowed and an element planted across it."""
    import foldmap.route_na as route_na
    from foldmap.layout import NT_STEP
    from foldmap.route import _clear, _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    tail = next(p for p in lay.na_paths if p.kind == "tail" and len(p.nts) == 10)
    (key,) = tail.between
    d = at[key]
    five = tail.nts[-1] < key[1]
    ahead = np.asarray(tail.nt_xy[0 if five else -1])
    out = ahead - np.asarray(tail.nt_xy[-1 if five else 0])
    out /= np.hypot(*out)
    centre = ahead - out * 2.0  # a block across the end's middle, 3 wide
    block = (*np.minimum(centre - 1.5, centre + 1.5), *np.maximum(centre - 1.5, centre + 1.5))
    rects, owners = _element_obstacles(lay)
    elements = (np.hstack([rects, np.array(block, float).reshape(4, 1)]), [*owners, "block"])
    grid = _Grid(lay)
    grid.blocked |= grid._mask([block])
    router = _Router(grid, elements)
    monkeypatch.setattr(route_na, "_TURNS", (0.0,))  # the straight drawing only
    path, tip, _ = route_na._tail(router, lay, d, key[0], key, tail.nts, five, None)
    assert not path.fallback and path.nts == tail.nts and len(path.nt_xy) == 10
    assert _clear(route_na.outside_own(path.points, (d,)), (d,), *elements, own_body=True)
    steps = np.hypot(*np.diff(np.asarray(path.nt_xy), axis=0).T)
    assert np.all(steps <= NT_STEP + 1e-6) and np.all(steps > 0.5 * NT_STEP)


def test_a_free_end_keeps_out_of_its_partners_doorway(monkeypatch):
    """A free end is not turned across the way its partner nucleotide's strand leaves the same end of the stem, when
    another drawing is clear (5F9R: the non-target strand's free end D:12-20, turned across the face of DNA 1, closed
    off the run from P1b, which had to cross it). 4YBB's 10-nt end AA:557-566, offered first a turn across the
    doorway."""
    import foldmap.route_na as route_na
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    tail = next(p for p in lay.na_paths if p.kind == "tail" and len(p.nts) == 10)
    (key,) = tail.between
    d, five = at[key], tail.nts[-1] < key[1]
    door = route_na._doorway(d, key, five)
    turns = [
        t
        for t in (np.pi / 2, -np.pi / 2)
        if route_na._crosses_segment(route_na._free_end(d, key, 10, five, t)[0], *door)
    ]
    assert len(turns) == 1  # the turn towards the partner runs across its doorway
    monkeypatch.setattr(route_na, "_TURNS", (turns[0], 0.0))
    router = _Router(_Grid(lay), _element_obstacles(lay))
    path, _, _ = route_na._tail(router, lay, d, key[0], key, tail.nts, five, None)
    assert not path.fallback and not route_na._crosses_segment(path.points, *door)


@pytest.mark.parametrize("pdb, runs, protein", [("2GIS", 6, 0), ("4OO8", 4, 13)])
def test_runs_keep_off_the_strands_drawn_before_them(pdb, runs, protein):
    """A run with nucleotides on it is routed off every strand drawn before it, and every protein loop, unless that
    is a long detour (route_na._DETOUR); failing that it pays for the cells beside them too (_BESIDE_DETOUR), else
    crossing them is only a cost; neither detour more than _DETOUR_STEPS steps per nucleotide longer. Measured with
    crossing them only a cost: run-run crossings 2GIS 8, 4OO8 5 (protein loops crossed 15). (30QM, pinned here while
    drawn in the plain frame, is drawn unrolled since the C4 revision.) 4OO8 protein 12 -> 13 when its helices were
    lettered (P2a/P2b: the longer names move the placement)."""
    from foldmap.nacheck import na_crossings

    lay, _, _, loops = _routed(pdb)
    m = na_crossings(lay, loops)
    assert m["na_run_crossings"] <= runs and m["na_protein_crossings"] <= protein


@pytest.mark.parametrize("pdb", [*FIXTURES, "2GIS", "2OEU"])
def test_no_nucleotide_sits_on_another_strands(pdb):
    """No nucleotide is drawn within a contact bead's width of a nucleotide of another strand or stem: where two runs
    cross or pass close (4OO8: the two sides of the P2-P3 internal loop cross where the stems are seen to), a run's
    nucleotides slide along it (4OO8's B:27 and B:42 were 0.25 apart, their beads on top of each other)."""
    from foldmap.route_na import _BEAD

    lay, _, _, _ = _routed(pdb)
    owner = [("stem", d.id) for d in lay.stems for _ in d.axial] + [
        ("path", k) for k, p in enumerate(lay.na_paths) for _ in p.nt_xy
    ]
    xy = np.array([d.point(key) for d in lay.stems for key in d.axial] + [q for p in lay.na_paths for q in p.nt_xy])
    close = [
        (owner[a], owner[b])
        for a in range(len(xy))
        for b in np.flatnonzero(np.hypot(*(xy[a + 1 :] - xy[a]).T) < _BEAD) + a + 1
        if owner[a] != owner[b]
    ]
    assert not close


def test_a_crowded_letter_looks_further_out_and_along_its_strand():
    """A base letter with no clear spot beside its strand looks further out and further along it before settling on
    the least bad spot (5B43: the 5' end B:-19--16 squeezed between its stem and a sheet; B:-17's letter landed
    0.15 from B:-18's). 1EHZ's acceptor end with an element planted over every nearer spot of its second letter."""
    from types import SimpleNamespace

    from foldmap.layout import NT_ROW
    from foldmap.route_na import _LETTER_GAP, _near, _place_letters

    lay, _, _, _ = _routed("1EHZ")
    tail = next(p for p in lay.na_paths if p.kind == "tail")
    key = (tail.strand, tail.nts[1])
    xy = np.asarray(tail.nt_xy[1])
    ahead = np.subtract(tail.nt_xy[2], tail.nt_xy[0])
    angle = float(np.arctan2(ahead[1], ahead[0]))
    block = SimpleNamespace(cx=xy[0], cy=xy[1], length=0.5, width=2 * 2.4 * NT_ROW + 0.4, angle=angle)
    lay.placed = {**lay.placed, "block": block}
    letters = _place_letters(lay, lay.na_paths)
    spot = letters[key]
    others = [q for k, q in letters.items() if k != key]
    assert not _near(spot, block) and min(np.hypot(*(np.subtract(others, spot)).T)) >= _LETTER_GAP


def _straight(length):
    from foldmap.route_na import NAPath

    return NAPath("run", 0, [1, 2], [(0.0, 0.0), (0.0, 0.0), (length, 0.0), (length, 0.0)], [])


@pytest.mark.parametrize("stage", ["first", "beside"])
def test_a_later_stage_takes_over_only_within_its_detour_caps(stage):
    """A staged drawing of a run replaces its plain one only within both caps: `bound` times as long (_DETOUR for
    the drawing off every strand and loop, _BESIDE_DETOUR for the one paying for the cells beside them) and at most
    _DETOUR_STEPS nucleotide steps longer per nucleotide and end (4YBB's P72->P71 went 5.7 times its direct distance
    round a cluster of stems)."""
    from foldmap.route_na import _BESIDE_DETOUR, _DETOUR, _DETOUR_STEPS, NT_STEP, _length, _takes_over

    bound = _DETOUR if stage == "first" else _BESIDE_DETOUR
    m, plain = 2, _straight(20.0)
    extra = _DETOUR_STEPS * NT_STEP * (m + 2)  # 13.2: under (bound - 1) x 20 for both stages
    assert extra < (bound - 1) * _length(plain.points)
    assert _takes_over(_straight(20.0 + extra - 0.1), plain, bound, m, ())
    assert not _takes_over(_straight(20.0 + extra + 0.3), plain, bound, m, ())  # over the step cap
    short = _straight(4.0)  # here the ratio is the tighter cap
    assert not _takes_over(_straight(4.0 * bound + 0.3), short, bound, m, ())


def test_a_later_stage_never_replaces_a_clean_drawing_with_a_fallback():
    from foldmap.route_na import _DETOUR, _takes_over

    plain, staged = _straight(10.0), _straight(10.5)
    staged.fallback = True
    assert not _takes_over(staged, plain, _DETOUR, 2, ())
    plain.fallback = True
    assert _takes_over(staged, plain, _DETOUR, 2, ())  # both fallbacks: the later stage may still be taken


def test_a_staged_strand_counts_the_router_debug_counts_of_its_drawing_only(monkeypatch):
    """router.tight and router.rejected count once per strand routed: the attempts a later stage discards (here
    the plain drawing and the first stage) leave no count behind."""
    import foldmap.route_na as route_na
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("1EHZ")
    at = {key: d for d in lay.stems for key in d.axial}
    run = next(p for p in lay.na_paths if p.kind == "run" and p.nts)
    (s, i), (_, j) = run.between
    router = _Router(_Grid(lay), _element_obstacles(lay))
    plain, staged = copy.copy(run), copy.copy(run)
    plain.fallback = staged.fallback = False

    calls = []

    def drawing(router, *args):  # every attempt counts one of each; level 2 draws nothing
        router.tight += 1
        router.rejected += 1
        calls.append(args[-1])
        return None if args[-1] == 2 else plain if len(calls) == 1 else staged

    monkeypatch.setattr(route_na, "_drawing", drawing)
    strands = route_na._Strands(router.grid)
    a, b = np.asarray(plain.points[1]), np.asarray(plain.points[2])
    mid, normal = (a + b) / 2, np.array([a[1] - b[1], b[0] - a[0]])
    strands.add([tuple(mid - normal), tuple(mid + normal)])  # a strand across the plain drawing
    extra = np.zeros_like(router.grid.used)
    got = route_na._strand_between(router, s, i, j, at[(s, i)], at[(s, j)], None, (), strands, None, extra)
    assert calls == [1, 2, 1] and got is staged and (router.tight, router.rejected) == (1, 1)


def _in_outline(pt, e) -> bool:
    c, s = np.cos(e.angle), np.sin(e.angle)
    x, y = pt[0] - e.cx, pt[1] - e.cy
    return abs(x * c + y * s) < e.length / 2 and abs(y * c - x * s) < e.width / 2


@pytest.mark.parametrize("pdb", FIXTURES)
def test_end_labels_keep_clear_of_loops_elements_and_other_labels(pdb):
    from foldmap.layout import _GLYPH_H, _GLYPH_W, centred_box, label_boxes
    from foldmap.route import clip_len

    lay, sses, _, loops = _routed(pdb)
    hw, hh = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2
    boxes = [centred_box(x, y, hw, hh) for _, _, (x, y) in lay.na_ends]
    others = [b for _, b in label_boxes(lay, sses)]
    bodies = [*lay.placed.values(), *lay.stems]
    for k, (s, end, _) in enumerate(lay.na_ends):
        box = boxes[k]
        assert not any(_overlap(box, o) for o in others + boxes[:k] + boxes[k + 1 :]), (s, end)
        for loop in loops:
            assert all(clip_len(p, q, box) <= 1e-9 for p, q in zip(loop.points, loop.points[1:])), (s, end, loop.a_id)
        samples = [(x, y) for x in np.linspace(box[0], box[2], 4) for y in np.linspace(box[1], box[3], 4)]
        assert not any(_in_outline(q, e) for q in samples for e in bodies), (s, end)


def test_letters_are_placed_quickly_on_a_whole_ribosomal_rna(monkeypatch):
    """Each letter looks only at its neighbours (a cell hash), not at every letter placed before it: on the whole
    4YBB small-subunit RNA (1000 nt, na_extent=all) placing them took 6.5 s when it was quadratic."""
    import time

    import foldmap.route_na as route_na
    from foldmap.style import Style

    lay, sses, bb = make_layout(_path("4YBB"), look=Style(na_extent="all"))
    assert lay.nucleic.size == 1000
    real, spent = route_na._place_letters, []

    def timed(*a, **k):
        t = time.perf_counter()
        out = real(*a, **k)
        spent.append(time.perf_counter() - t)
        return out

    monkeypatch.setattr(route_na, "_place_letters", timed)
    route_loops(lay, sses, bb)
    assert len(lay.na_letters) > 400 and spent[0] < 1.5, spent


def test_unlabelled_ties_are_drawn_like_rungs_solid_for_watson_crick_dotted_for_gu():
    """The rung convention holds for tie lines: a WC tie is a thin solid line, a G·U tie dotted, both rung grey."""
    from foldmap import render
    from foldmap.route_na import NATie

    lay, sses, _, loops = _routed("1EHZ")
    pairs = lay.nucleic.base_pairs
    wc = next(p for p in pairs if p.kind == "WC")
    gu = next(p for p in pairs if p.kind == "GU")
    ties = []
    for p in (wc, gu):
        a, b = (p.a, p.i), (p.b, p.j)
        ties.append(NATie(a, b, [(0.0, 0.0), (1.0, 0.0)]))
    lay.na_ties = [*lay.na_ties, *ties]
    fig = render.draw(lay, loops, sses)
    drawn = {t.get_gid(): t for t in _gid(fig, "na:tie:") if not t.get_gid().endswith((":mark", ":label"))}
    tag = lambda p: f"na:tie:{lay.nucleic.strands[p.a].labels[p.i].tag}-{lay.nucleic.strands[p.b].labels[p.j].tag}"  # noqa: E731
    line_wc, line_gu = drawn[tag(wc)], drawn[tag(gu)]
    assert line_wc.get_linestyle() == "-"
    assert line_gu.get_linestyle() != "-"  # dotted, as the G·U rungs
    rung = _gid(fig, "na:stem:0:rung")[0].get_color()
    assert line_wc.get_color() == line_gu.get_color() == rung


# ---- backbone colours: modality, chain, override -------------------------------------------------------------

LMB = Path(__file__).parent / "data" / "1LMB.cif"


def _backbone(fig):
    """Colours of everything drawn as nucleic-acid backbone: stem rails and single-strand paths."""
    ax = fig.axes[0]
    rails = [
        to_hex(p.get_edgecolor())
        for p in ax.patches
        if (p.get_gid() or "").startswith("na:stem:") and "rails" in p.get_gid()
    ]
    runs = [to_hex(a.get_color()) for a in ax.lines if (a.get_gid() or "").startswith("na:run:")]
    return rails, runs


def _fig(path, *sets):
    from foldmap.style import resolve_style

    return make_figure_and_loops(path, look=resolve_style(None, None, [f"{s}" for s in sets]))[0]


def _loop_colours(fig):
    return {to_hex(p.get_edgecolor()) for p in fig.axes[0].patches if (p.get_gid() or "").startswith("loop:")}


def test_dna_and_rna_backbones_have_their_own_colours_and_protein_loops_stay_ink():
    rails, runs = _backbone(_fig(LMB))
    assert rails and set(rails + runs) == {"#8c564b"}
    assert _loop_colours(_fig(LMB)) == {"#000000"}
    rails, runs = _backbone(_fig(_path("1EHZ")))
    assert rails and set(rails + runs) == {"#6a3d9a"}


def test_every_backbone_piece_takes_the_backbone_colour_and_letters_stay_ink():
    fig = _fig(_path("1EHZ"))
    ax = fig.axes[0]
    pieces = [a for a in [*ax.patches, *ax.lines] if (a.get_gid() or "").startswith(("na:stem:", "na:run:", "na:tie:"))]
    coloured = {
        to_hex(a.get_edgecolor() if hasattr(a, "get_edgecolor") else a.get_color())
        for a in pieces
        if not a.get_gid().endswith(("rung", ":label")) and a.get_gid().count(":") >= 2
    }
    assert "#6a3d9a" in coloured
    assert {to_hex(t.get_color()) for t in _gid(fig, "na:nt:") if t.get_fontweight() == "normal"} == {"#3d4852"}
    assert {to_hex(a.get_color()) for a in _gid(fig, "na:stem:0:rung")} == {"#c3cad1"}


def test_na_backbone_ink_reproduces_the_old_look():
    rails, runs = _backbone(_fig(LMB, "na_backbone=ink"))
    assert set(rails + runs) == {"#3d4852"}


def test_na_backbone_colour_paints_all_nucleic_acid():
    rails, runs = _backbone(_fig(_path("1EHZ"), "na_backbone=#123456"))
    assert set(rails + runs) == {"#123456"}


def test_na_backbone_chain_gives_each_chain_its_own_colour_avoiding_protein_chains():
    from foldmap.palette import chain_colors, delta_e

    rails, _ = _backbone(_fig(LMB, "na_backbone=chain"))
    assert len(set(rails)) == 2  # chains 1 and 2
    protein = list(chain_colors(["3", "4"]).values())
    assert all(delta_e(c, p) >= 15 for c in rails for p in protein)


def test_backbone_colors_chain_beats_modality_and_recolours_protein_loops():
    fig = _fig(LMB, "backbone_colors=3=#1f77b4,dna=#555555,2=crimson")
    lay, sses, bb = make_layout(LMB)
    route_loops(lay, sses, bb)  # routes the strands: lay.na_paths
    chain_of = {i: s.chain for i, s in enumerate(lay.nucleic.strands)}
    runs = {chain_of[p.strand]: to_hex(a.get_color()) for p, a in zip(lay.na_paths, _gid(fig, "na:run:"))}
    assert runs == {"1": "#555555", "2": "#dc143c"}  # chain 2 beats dna; chain 1 takes the modality colour
    rails = {to_hex(p.get_edgecolor()) for p in _gid(fig, "na:stem:") if ":rails" in p.get_gid()}
    assert rails == {"#555555", "#dc143c"}
    assert _loop_colours(fig) == {"#1f77b4", "#000000"}  # chain 3 blue, chain 4 untouched


def test_backbone_colors_protein_key_recolours_every_protein_chain():
    assert _loop_colours(_fig(LMB, "backbone_colors=protein=#1f77b4")) == {"#1f77b4"}


@pytest.mark.parametrize(
    "text,named",
    [("Z=red", "Z"), ("3=notacolour", "3"), ("dna=#12", "dna"), ("justakey", "justakey")],
)
def test_bad_backbone_colors_name_the_key(text, named):
    with pytest.raises(ValueError, match="backbone_colors") as err:
        _fig(LMB, f"backbone_colors={text}")
    assert named in str(err.value)


def test_bad_na_backbone_is_explained():
    with pytest.raises(ValueError, match="na_backbone must be"):
        _fig(LMB, "na_backbone=nonsense")


def test_the_legend_entry_for_nucleic_acid_uses_the_backbone_colour():
    from foldmap import render
    from foldmap.palette import chain_colors
    from foldmap.style import Style, resolve_style

    lay, _, _ = make_layout(LMB)
    colors = chain_colors(["3", "4"])
    assert ("#8c564b", "DNA (1, 2)") in render._legend_entries(
        lay, Style(), ["3", "4"], colors, [], {"1": "#8c564b", "2": "#8c564b"}
    )
    look = resolve_style(None, None, ["backbone_colors=dna=#555555"])
    assert ("#555555", "DNA (1, 2)") in render._legend_entries(
        lay, look, ["3", "4"], colors, [], {"1": "#555555", "2": "#555555"}
    )


def test_colours_never_move_the_layout():
    a, _, _ = make_layout(LMB)
    from foldmap.style import resolve_style

    b, _, _ = make_layout(LMB, look=resolve_style(None, None, ["na_backbone=chain", "backbone_colors=3=red"]))
    assert a.bounds == b.bounds and [d.rect for d in a.stems] == [d.rect for d in b.stems]


@pytest.mark.parametrize(
    "style, parts", [("ladder", ":rails"), ("ribbon", ":front"), ("block", ":block"), ("letters", ":rung")]
)
def test_each_stem_style_draws_its_own_parts(style, parts):
    from foldmap.style import Style

    fig, _ = make_figure_and_loops(_path("30QM"), look=Style(na_style=style))
    found = _gid(fig, "na:stem:")
    assert {a.get_gid().split(":")[2] for a in found} == {str(k) for k in range(6)}
    assert any(a.get_gid().endswith(parts) for a in found)
    letters = _gid(fig, "na:nt:")
    if style == "block":
        assert not letters
    elif style == "ribbon":
        assert len(letters) == 121 - 2 * 32  # single strands only: the double helix carries no letters
    else:
        assert len(letters) == 121


def test_every_block_bar_numbers_its_first_pair_beside_it():
    from foldmap.style import Style

    look = Style(na_style="block")
    fig, _ = make_figure_and_loops(_path("30QM"), look=look)
    lay, _, _ = make_layout(_path("30QM"), look=look)
    firsts = {a.get_gid(): a for a in _gid(fig, "na:stem:") if a.get_gid().endswith(":first")}
    ladder, _ = make_figure_and_loops(_path("30QM"), look=Style(na_style="ladder"))
    letter = _gid(ladder, "na:nt:")[0].get_fontsize()  # block draws no letters: the size they would have
    lay, sses, _ = make_layout(_path("30QM"), look=Style(na_style="block"))
    boxes = dict(label_boxes(lay, sses))
    assert len(firsts) == len(lay.stems) == 6
    for k, d in enumerate(lay.stems):
        a, i, b, j = d.pairs[0]
        text = firsts[f"na:stem:{k}:first"]
        assert text.get_text() == f"{lay.nucleic.strands[a].labels[i].seq}\u00b7{lay.nucleic.strands[b].labels[j].seq}"
        assert text.get_fontsize() == pytest.approx(letter)
        x0, y0, x1, y1 = boxes[f"stem-label:{k}"]  # the slot the layout reserves beside the stem
        x, y = text.get_position()
        assert x0 <= x <= x1 and y0 <= y <= y1
        assert text.get_rotation() % 360 < 90 or text.get_rotation() % 360 > 270  # reads left to right


def test_the_layout_does_not_depend_on_na_style():
    from foldmap.style import Style

    layouts = {
        s: make_layout(_path("30QM"), look=Style(na_style=s))[0] for s in ("ladder", "ribbon", "block", "letters")
    }
    base = layouts["ladder"]
    for s, lay in layouts.items():
        assert [d.rect for d in lay.stems] == [d.rect for d in base.stems], s
        assert {k: (p.cx, p.cy) for k, p in lay.placed.items()} == {k: (p.cx, p.cy) for k, p in base.placed.items()}, s


def test_a_cropped_large_rna_switches_ladder_to_block():
    fig, _ = make_figure_and_loops(_path("4YBB"))
    assert any(a.get_gid().endswith(":block") for a in _gid(fig, "na:stem:"))
    assert _gid(fig, "na:gap:")


node = shutil.which("node")


def test_explorer_lists_stems_and_strands():
    from foldmap.interactive import build_page

    page = build_page(_path("30QM"))
    blob = page.split('<script type="application/json" id="topo-data">', 1)[1].split("</script>", 1)[0]
    data = json.loads(blob)
    keys = [x["key"] for x in data["links"]["nucleics"]]
    assert keys[:6] == [f"stem:{k}" for k in range(6)]
    p1 = data["links"]["nucleics"][0]
    assert p1["text"].startswith("P1 \u00b7 7 bp") and p1["het"] and p1["residues"]
    assert 'id="na:stem:0:rails"' in page


@pytest.mark.skipif(node is None, reason="node not installed")
def test_explorer_maps_stem_artists_to_their_stem(tmp_path):
    from foldmap.interactive import _LIB

    probe = tmp_path / "probe.js"
    probe.write_text(
        _LIB
        + """
const lib = globalThis.TopoLib;
console.log(JSON.stringify(['na:stem:3', 'na:stem:3:rails', 'na:label:3', 'na:run:5', 'na:nt:A:12', 'helix:A:1-9']
  .map(lib.nucleicKey)));
""",
        encoding="utf-8",
    )
    res = subprocess.run([node, str(probe)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert json.loads(res.stdout) == ["stem:3", "stem:3", "stem:3", "run:5", None, None]


def test_explorer_page_of_a_nucleic_acid_only_file_has_no_contact_map():
    from foldmap.interactive import build_page

    page = build_page(_path("1EHZ"))
    data = json.loads(page.split('<script type="application/json" id="topo-data">', 1)[1].split("</script>", 1)[0])
    assert data["residues"] == [] and data["ca"] == []
    assert data["links"]["nucleics"][0]["key"] == "stem:0" and data["links"]["nucleics"][0]["het"]
    app = page.split('<script id="topo-app">', 1)[1]
    assert "data.ca.length > 0" in app and "if (!hasMap) return;" in app


_HARNESS = r"""
// Runs the explorer's app script against stub DOM and viewer objects; drives viewer.on / viewer.fly.
const fs = require('fs');
const job = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const calls = [];
const generic = () => new Proxy(function () {}, {
  get(_, k) {
    if (k === 'getAttribute') return () => '0 0 10 10';
    if (k === 'textContent') return '';
    if (k === 'dataset') return {};
    if (k === 'style') return {};
    if (k === 'classList') return { add() {}, remove() {}, toggle() {} };
    if (k === 'querySelectorAll') return () => [];
    if (k === 'getContext') return () => ({ createImageData: () => ({ data: { set() {} } }), putImageData() {},
                                            fillRect() {},
                                            clearRect() {} });
    if (k === 'Symbol.toPrimitive') return undefined;
    return generic();
  },
  set() { return true; },
  apply() { return generic(); },
});
const elements = {};
globalThis.document = {
  getElementById: id => id === 'topo-data' ? { textContent: job.blob } : (elements[id] = elements[id] || generic()),
  querySelector: () => generic(), querySelectorAll: () => [], addEventListener() {},
};
globalThis.window = globalThis;
globalThis.addEventListener = () => {};
globalThis.ResizeObserver = class { observe() {} };
globalThis.URL.createObjectURL = () => 'blob:x';
globalThis.$3Dmol = { createViewer: () => new Proxy({
  setStyle: (s, st) => calls.push(['setStyle', s, st]), zoomTo: s => calls.push(['zoomTo', s]),
}, { get: (t, k) => k in t ? t[k] : () => {} }) };
globalThis.PDBeMolstarPlugin = class {
  constructor() {
    this.events = { loadComplete: { subscribe: cb => { this.done = cb; } } };
    this.visual = { select: async p => { calls.push(['select', p]); }, focus: async q => { calls.push(['focus', q]); },
                    reset() {} };
    this.plugin = { canvas3d: { handleResize() {},
                                camera: { getSnapshot: () => ({ target: [0, 0, 0], position: [1, 1, 1],
                                                                                  radius: 1 }), setState() {} } } };
  }
  render() { setTimeout(() => this.done(true), 0); }
};
eval(job.lib);
eval(job.app);
setTimeout(async () => {
  const out = [];
  for (const c of job.cases) {
    for (const op of ['on', 'fly']) {
      calls.length = 0;
      let error = null;
      const het = c.het === null ? undefined : c.het;
      try { await window.foldmapViewer[op](c.residues, het); } catch (e) { error = String(e); }
      await new Promise(r => setTimeout(r, 600));
      out.push({ case: c.name, op, error, calls: JSON.parse(JSON.stringify(calls)) });
    }
  }
  console.log(JSON.stringify(out));
  process.exit(0);
}, 100);
"""


def _explorer_cases(pdb, viewer):
    from foldmap.interactive import _APP, _LIB, build_page

    page = build_page(_path(pdb), viewer=viewer)
    blob = page.split('<script type="application/json" id="topo-data">', 1)[1].split("</script>", 1)[0]
    # the blob is HTML-safe JSON (< escapes); the script reads it back through JSON.parse unchanged
    return json.loads(blob), blob, _LIB, _APP


def _drive(tmp_path, pdb, viewer, cases):
    data, blob, lib, app = _explorer_cases(pdb, viewer)
    job = {"blob": blob, "lib": lib, "app": app, "cases": cases(data)}
    (tmp_path / "job.json").write_text(json.dumps(job), encoding="utf-8")
    (tmp_path / "harness.js").write_text(_HARNESS, encoding="utf-8")
    res = subprocess.run(
        [node, str(tmp_path / "harness.js"), str(tmp_path / "job.json")], capture_output=True, text=True
    )
    assert res.returncode == 0, res.stderr
    return data, {(r["case"], r["op"]): r for r in json.loads(res.stdout)}


def _stem_with_contacts(data):
    return next(x for x in data["links"]["nucleics"] if x["key"].startswith("stem:") and x["residues"])


@pytest.mark.skipif(node is None, reason="node is not installed; the explorer's viewer code cannot be exercised")
@pytest.mark.parametrize("viewer", ["3dmol", "molstar"])
def test_explorer_viewers_take_protein_nucleic_and_mixed_selections(tmp_path, viewer):
    def cases(data):
        helix = data["elements"][0]
        mixed = _stem_with_contacts(data)
        return [
            {"name": "protein", "residues": list(range(helix["start"], helix["end"] + 1)), "het": []},
            {"name": "mixed", "residues": mixed["residues"], "het": mixed["het"]},
            {"name": "protein-no-het-argument", "residues": [data["elements"][0]["start"]], "het": None},
        ]

    data, got = _drive(tmp_path, "1IHF", viewer, cases)
    helix = data["elements"][0]
    n = helix["end"] - helix["start"] + 1
    mixed = _stem_with_contacts(data)
    for r in got.values():
        assert r["error"] is None, r
    key = "setStyle" if viewer == "3dmol" else "select"
    for op in ("on", "fly"):
        assert got[("protein", op)]["calls"], op  # something was sent to the viewer
    if viewer == "3dmol":
        on = got[("protein", "on")]["calls"]
        assert any(c[0] == "setStyle" and len(c[1].get("or", [])) == n for c in on)
        assert not any("stick" in c[2] and c[2]["stick"].get("color") for c in on)  # no nucleotide styling
        assert len(got[("protein", "fly")]["calls"][0][1]["or"]) == n
        on = got[("mixed", "on")]["calls"]
        assert any(len(c[1].get("or", [])) == len(mixed["residues"]) for c in on if c[0] == key)
        assert sum(1 for c in on if c[0] == key and "resi" in c[1]) == len(mixed["het"])
        assert len(got[("mixed", "fly")]["calls"][0][1]["or"]) == len(mixed["residues"]) + len(mixed["het"])
        assert got[("protein-no-het-argument", "fly")]["calls"][0][1]["or"]
    else:
        parts = got[("protein", "on")]["calls"][0][1]["data"]
        assert sum(p["end_auth_seq_id"] - p["beg_auth_seq_id"] + 1 for p in parts) >= n
        assert not any(p.get("representation") for p in parts)
        parts = got[("mixed", "on")]["calls"][0][1]["data"]
        assert sum(1 for p in parts if p.get("representation") == "ball-and-stick") == len(mixed["het"])
        focus = got[("mixed", "fly")]["calls"][0][1]
        assert sum(1 for q in focus if "auth_seq_id" in q) == len(mixed["het"])
        assert got[("protein-no-het-argument", "fly")]["calls"]


@pytest.mark.skipif(node is None, reason="node is not installed; the explorer's viewer code cannot be exercised")
@pytest.mark.parametrize("viewer", ["3dmol", "molstar"])
def test_explorer_viewers_take_a_nucleic_acid_only_selection(tmp_path, viewer):
    def cases(data):
        stem = data["links"]["nucleics"][0]
        return [{"name": "rna", "residues": [], "het": stem["het"]}]

    data, got = _drive(tmp_path, "1EHZ", viewer, cases)
    het = data["links"]["nucleics"][0]["het"]
    assert got[("rna", "on")]["error"] is None and got[("rna", "fly")]["error"] is None
    if viewer == "3dmol":
        styled = [c for c in got[("rna", "on")]["calls"] if c[0] == "setStyle" and "resi" in c[1]]
        assert len(styled) == len(het)
        assert len(got[("rna", "fly")]["calls"][0][1]["or"]) == len(het)
    else:
        assert len(got[("rna", "on")]["calls"][0][1]["data"]) == len(het)
        assert len(got[("rna", "fly")]["calls"][0][1]) == len(het)


def test_explorer_texts_name_both_chains_and_single_nucleotides():
    from foldmap.interactive import build_page

    texts = {}
    pages = {}
    for pdb in ("3CRX", "1IHF"):
        pages[pdb] = build_page(_path(pdb))
        data = json.loads(
            pages[pdb].split('<script type="application/json" id="topo-data">', 1)[1].split("</script>", 1)[0]
        )
        texts[pdb] = [x["text"] for x in data["links"]["nucleics"]]
        for t in texts[pdb]:
            assert "–-" not in t, t  # never "C -50–-50"
    assert "DNA 1a · 9 bp · C -49 … E 49" in texts["1IHF"]  # a stem across two chains names both
    assert "tail C -50 (1 nt)" in texts["1IHF"]  # a single nucleotide is not a range
    assert "tail C1 1 (1 nt)" in texts["3CRX"]  # the assembly names its chains C1, D1, ...
    assert "<p>5 chains · 190 residues · 70 nt · " in pages["1IHF"]  # protein page header counts nucleotides
    rna = build_page(_path("1EHZ"))
    assert "<p>1 chain · 76 nt · hover" in rna


def test_summary_lists_the_stems(capsys):
    from foldmap.cli import main

    assert main(["summary", str(_path("30QM"))]) == 0
    out = capsys.readouterr().out
    assert "32 base pairs" in out and "P1" in out and "1:1-7 · 1:112-118" in out


@pytest.mark.parametrize("pdb", SMALL)
def test_stems_sit_nearer_the_elements_that_touch_them(pdb):
    from foldmap.nacheck import contacts

    lay, sses, _ = make_layout(_path(pdb))
    m = contacts(lay, sses)
    assert m["contact_auc"] is None or m["contact_auc"] >= 0.55  # IHF, wrapped by its arms, is the lowest


def test_a_chain_break_is_drawn_as_a_gap_with_no_ends_at_it():
    """1URN: the unmodelled U13-C14 of the hairpin loop is a dashed gap from R:12 to R:15, with no 5'/3' labels
    there; the hairpin is not drawn as the classic straight two-strand duplex."""
    from foldmap.frame import RadialFrame

    lay, _, _, _ = _routed("1URN")
    assert [s.id for s in lay.nucleic.strands] == ["R:1-12", "R:15-20"]
    assert sorted((s, e) for s, e, _ in lay.na_ends) == [(0, "5′"), (1, "3′")]
    (gap,) = [p for p in lay.na_paths if p.kind == "gap"]
    assert gap.strand == 0 and gap.label == ""  # 2 nt: a bare dashed break
    assert not isinstance(lay.frame, RadialFrame)


def test_a_trace_is_an_obstacle_the_protein_keeps_off():
    """1QBJ: hexamer F pairs only with a crystal mate, so it is drawn as a trace at its projected place; chain C's
    alpha1 was placed on it (the trace a routing fallback). The trace's box is a fixed obstacle in placement."""
    lay, _, _, _ = _routed("1QBJ")
    (pts,) = lay.na_traces.values()
    box = (*np.min(pts, axis=0), *np.max(pts, axis=0))
    assert not [p.sse.id for p in lay.placed.values() if _overlap(p.rect, box)]
    assert not [d.id for d in lay.stems if _overlap(d.rect, box)]
    assert not [p.kind for p in lay.na_paths if p.fallback]


def test_the_trace_boxes_reach_the_optimiser_on_the_grid(monkeypatch):
    import foldmap.layout as L

    seen, real = [], L._energy_dense

    def dense(*a):
        seen.append(a[8])
        return real(*a)

    monkeypatch.setattr(L, "_energy_dense", dense)
    make_layout(_path("1QBJ"))
    centres, want = seen[0]
    assert len(centres) == 1
    for v in (centres, want):
        v = np.asarray(v) * 1e4
        assert np.all(np.abs(v - np.round(v)) < 1e-6)


def _gap_label_box(p):
    from foldmap.layout import _GLYPH_H, _GLYPH_W

    x, y = p.label_xy
    hw = _GLYPH_W * len(p.label) / 2 + 0.05
    return (x - hw, y - _GLYPH_H / 2, x + hw, y + _GLYPH_H / 2)


def test_a_gap_label_is_placed_clear_of_letters_beads_bodies_and_strands():
    """4YBB: "... 21 nt ..." sat at its gap's midpoint on element AC:26-46 and two runs (C10). A gap label is placed
    along its gap, either side, where its box keeps clear of bodies, strands, letters and nucleotides."""
    from foldmap.route import clip_len
    from foldmap.route_na import crosses_body

    lay, _, _, loops = _routed("4YBB")
    labelled = [p for p in lay.na_paths if p.kind == "gap" and p.label]
    assert labelled
    letters = [*lay.na_letters.values(), *[d.point(k, 0.0) for d in lay.stems for k in d.axial]]
    segments = [(a, b) for q in lay.na_paths for a, b in zip(q.points, q.points[1:])]
    segments += [(a, b) for loop in loops for a, b in zip(loop.points, loop.points[1:])]
    for p in labelled:
        x0, y0, x1, y1 = box = _gap_label_box(p)
        edges = [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]
        assert not any(x0 < x < x1 and y0 < y < y1 for x, y in letters)
        assert not any(clip_len(a, b, box) > 1e-9 for a, b in segments)
        bodies = [*lay.placed.values(), *lay.stems]
        assert not any(crosses_body(a, b, e) for e in bodies for a, b in edges)


def test_a_gap_label_moves_off_a_letter_planted_on_its_first_spot():
    from foldmap.route_na import _label_gaps

    lay, _, _, loops = _routed("4YBB")
    p = next(p for p in lay.na_paths if p.kind == "gap" and p.label)
    first = p.label_xy
    lay.na_letters[(-1, -1)] = first  # a letter where the label went
    _label_gaps(lay, lay.na_paths, [loop.points for loop in loops])
    x0, y0, x1, y1 = _gap_label_box(p)
    assert p.label_xy != first and not (x0 < first[0] < x1 and y0 < first[1] < y1)
