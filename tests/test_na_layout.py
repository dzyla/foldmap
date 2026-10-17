"""Nucleic acids on the page: stems placed like elements, strands routed, four drawing styles, the explorer."""

import copy
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex
from matplotlib.transforms import Bbox

from foldmap.cli import make_figure_and_loops, make_layout
from foldmap.layout import PlacedStem, label_boxes
from foldmap.route import route_loops
from helpers import text_centre

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
    from foldmap.layout import joined_pieces, lying_over

    lay, _, _ = make_layout(_path(pdb))
    pieces = len(lay.stems)
    assert pieces >= len(lay.nucleic.structure.stems) > 0
    joined = joined_pieces(lay.stems, lay.nucleic)  # consecutive pieces of a bent duplex meet at their joint
    joined |= lying_over(lay)  # and a helix drawn over a sheet lies on its strands
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


def test_dna_held_by_both_protomers_is_not_unrolled():
    """3L2V (random-100 C16): each viral DNA end of the PFV intasome is gripped by both inner integrase chains.
    Unrolled side by side, each DNA sat beside one protomer, turned 18-28° (58° on the whole entry) from its 3D axis
    with the other protomer's contacts a protomer's width away. When most stem pieces are held by two protomers the
    plain frame is used (on the whole entry its dyad lies in the page: the two halves drawn as mirror images), and
    every stem is drawn along its seen axis."""
    from foldmap.frame import Frame
    from foldmap.nacheck import shape

    lay, _, _ = make_layout(NA.parent / "r100" / "3L2V-A.cif.gz")
    assert type(lay.frame) is Frame
    m = shape(lay)
    assert m["stems_drawn"] == 2 and m["angle_max"] == 0.0


def test_a_free_end_never_runs_along_another_strand():
    """9CJI (random-100 C16, A17 bead_overlaps 2): the crRNA's 7-nt 3' end was turned onto the line of its 5' end,
    which leaves the same end of the stem; nucleotides 3 and 24 were drawn 0.1 apart (their beads on top of each
    other). A free end's drawing whose nucleotides come within a bead's width of a strand drawn so far is passed
    over, as one that crosses its partner's doorway is."""
    from foldmap.na_paths import _BEAD

    lay, sses, bb = make_layout(NA.parent / "r100" / "9CJI-AB.cif.gz")
    route_loops(lay, sses, bb)
    pts = [(p.kind, k, np.asarray(xy, float)) for p in lay.na_paths for k, xy in zip(p.nts, p.nt_xy)]
    assert sum(kind == "tail" for kind, _, _ in pts) == 11
    close = [(a, b) for i, (_, a, p) in enumerate(pts) for _, b, q in pts[i + 1 :] if np.hypot(*(p - q)) < _BEAD]
    assert close == []
    assert not any(p.fallback for p in lay.na_paths)


def test_a_free_end_into_a_gap_turns_toward_where_the_molecule_goes_on():
    """3OWW (random-100 C16, A17 na_crossings 2): P5's hairpin loop lacks B22, so it is drawn as the 4-nt 3' end
    B:18-21 of one piece, a dashed gap and the 1-nt 5' end B:23 of the next. B:18-21 went straight on, away from
    B:23, and the gap came back round the stem (5.5 long), across a run on the whole figure. A 3' free end the protein
    does not touch, running into a gap, turns first toward the stem where the molecule goes on: the gap closes the
    loop beside the stem's end."""
    import math

    lay, sses, bb = make_layout(NA.parent / "r100" / "3OWW-P5.cif.gz")
    route_loops(lay, sses, bb)
    (gap,) = [p for p in lay.na_paths if p.kind == "gap"]
    length = sum(math.dist(a, b) for a, b in zip(gap.points, gap.points[1:]))
    assert length < 2.5 and not gap.fallback  # 2.1; straight on, 5.5


def test_one_shared_stem_among_several_stays_unrolled():
    """30QM: one of six stems (P3) is held by both RT protomers; the NA campaign's C4 ruling keeps it unrolled (the
    linked stems turn at most 57°): rule (c) of na_torn_by_unrolling needs most pieces shared."""
    from foldmap.frame import CylinderFrame

    lay, _, _ = make_layout(_path("30QM"))
    assert type(lay.frame) is CylinderFrame


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
    for t in labels:  # a base letter's size, upright, bold in the tick's colour; the legend says what it means
        base = _gid(fig, "na:nt:")[0]
        assert t.get_fontsize() == pytest.approx(base.get_fontsize())
        assert t.get_rotation() == 0 and t.get_fontweight() == "bold" and base.get_fontweight() == "normal"
        assert to_hex(t.get_color()) == to_hex(marks[0].get_color()) != to_hex(base.get_color())
    assert not [a for a in _gid(fig, "na:tie:") if not a.get_gid().endswith((":mark", ":label"))]


@pytest.mark.parametrize("pdb", FIXTURES)
def test_tie_lines_cross_no_stem_and_tie_labels_count_as_text(pdb):
    from foldmap.layout import label_boxes
    from foldmap.na_paths import crosses_body

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


def _page_em(d) -> float:
    """One em of a drawn figure's text in page units: a base letter is drawn at 0.75 em (render._draw_na_marks);
    points -> pixels (dpi / 72) -> page units (the axes' pixels per unit)."""
    ax = d.figure.axes[0]
    letter = next(t for t in ax.texts if (t.get_gid() or "").startswith("na:nt:"))
    px_per_unit = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    return letter.get_fontsize() / 0.75 * d.figure.dpi / 72 / px_per_unit


@pytest.mark.parametrize("pdb", ["1KX5", "30QM", "3CRX"])
def test_stem_names_keep_clear_of_base_letters(pdb):
    """No stem name sits on a base letter: its ink (figcheck, fixed hinting, Liberation Sans) keeps at least half the
    0.5 em margin from every letter's ink, its own row's (where it sits at the margin: plan 10 Task 4, which
    replaced the "a letter and a half" clearance of the one-glyph-out place) and every other stem's. A name and the
    letters turned alike with it (its own row) are compared in that turned frame (figcheck Page.pair_boxes): the
    page boxes of a long name turned 30 degrees reach a whole em past its ink."""
    from foldmap import figcheck
    from foldmap.layout import _INK_MARGIN
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(_path(pdb)))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    em, page = _page_em(d), figcheck.measure(d)
    px = d.figure.axes[0].transData.transform((1, 0))[0] - d.figure.axes[0].transData.transform((0, 0))[0]
    names = [(t, n, b) for t, n, b in page.texts if n.startswith("na:label:")]
    letters = [(t, n, b) for t, n, b in page.texts if n.startswith("na:nt:")]
    assert names and letters
    near = [
        (n, m)
        for tn, n, a in names
        for tm, m, b in letters
        if figcheck._box_gap(*(box.extents for box in page.pair_boxes(tn, a, tm, b))) < 0.5 * _INK_MARGIN * em * px
    ]
    assert not near, (pdb, near)


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
    from foldmap.na_paths import among
    from foldmap.route import _INFLATE, _clear, _element_obstacles

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


# runs over a stem-name box, held at the measured count (a held run that keeps off them must be removed): none since
# fix round 1 of Task 13b step 4 (4YBB's 1-nt run P1->P2 was boxed in and gave up a soft cell of P12b's name box)
OVER_NAMES: dict = {}


@pytest.mark.parametrize("pdb", FIXTURES)
def test_runs_keep_off_the_stem_names(pdb):
    from foldmap.layout import label_boxes
    from foldmap.na_paths import among

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
    from foldmap.na_paths import among

    lay, sses, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    names = [b for n, b in label_boxes(lay, sses) if n.startswith("stem-label:")]
    run = next(p for p in lay.na_paths if p.between == ((1, 391), (1, 393)))
    assert [at[k].name for k in run.between] == ["P17b", "P17a"] and not run.fallback
    pts = among(run, at)
    assert all(segment_clear(a, b, names) for a, b in zip(pts, pts[1:]))


def test_a_run_over_a_stem_name_is_drawn_when_nothing_else_is_clear():
    """The names are hard only in the first pass: with every drawing over a name, each run is drawn as without."""
    from foldmap.na_runs import _strand_between
    from foldmap.route import _element_obstacles, _Grid, _Router

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
    name, the bow is still the drawing taken (the soft choice), not a route nor a fallback (2GIS's P3a-P3b loop;
    4OO8's P2a-P2b sides, pinned here before, are routed since its helix is held together: their faces no longer
    see each other across the loop)."""
    from foldmap.layout import NT_STEP
    from foldmap.na_paths import _face
    from foldmap.na_runs import _side_bow, _strand_between
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("2GIS")
    at = {key: d for d in lay.stems for key in d.axial}
    kind = {(r.strand, r.lo): r.kind for r in lay.nucleic.structure.runs}
    sides = [
        p
        for p in lay.na_paths
        if p.kind == "run"
        and p.nts
        and kind.get((p.strand, p.nts[0])) == "internal"
        and {at[k].name for k in p.between} == {"P3a", "P3b"}
    ]
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
    from foldmap.na_paths import _face, _Strands
    from foldmap.na_runs import _straight_joint
    from foldmap.route import _element_obstacles, _Grid, _Router

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


def test_a_joint_no_route_reaches_is_drawn_straight_a_little_further():
    """C24 (random-100): a joint between two pieces of a bent duplex that no route reaches fell back to a plain
    curve when its faces were 1.5 or more apart and the pieces did not touch, though the straight line between them
    was clear (9C75 P2b -> P2a, faces 1.83 apart). Such a joint is drawn straight up to _JOINT_REACH, as a kink
    joins two pieces of a helix; further apart it is not."""
    from foldmap.na_paths import _face, _Strands
    from foldmap.na_runs import _JOINT_REACH, _straight_joint
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("1AOI")
    at = {key: d for d in lay.stems for key in d.axial}
    router = _Router(_Grid(lay), _element_obstacles(lay))
    reach = {}
    for p in (p for p in lay.na_paths if p.kind == "joint"):
        a, b = p.between
        (fa, ea), (fb, _) = _face(at[a], a, True), _face(at[b], b, False)
        reach[(at[a].name, at[b].name)] = (
            math.dist(fa, fb),
            _straight_joint(router, at[a], at[b], fa, ea, fb, _Strands(router.grid)),
        )
    assert any(1.5 < d < _JOINT_REACH and ok for d, ok in reach.values()), reach  # DNA 1b -> 1a, 1u -> 1t
    assert not any(ok for d, ok in reach.values() if d >= _JOINT_REACH)


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
    from foldmap.na_paths import among

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
    import foldmap.na_ends as na_ends
    from foldmap.layout import NT_STEP
    from foldmap.na_paths import outside_own
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
    monkeypatch.setattr(na_ends, "_TURNS", (0.0,))  # the straight drawing only
    path, tip, _ = na_ends._tail(router, lay, d, key[0], key, tail.nts, five, None)
    assert not path.fallback and path.nts == tail.nts and len(path.nt_xy) == 10
    assert _clear(outside_own(path.points, (d,)), (d,), *elements, own_body=True)
    steps = np.hypot(*np.diff(np.asarray(path.nt_xy), axis=0).T)
    assert np.all(steps <= NT_STEP + 1e-6) and np.all(steps > 0.5 * NT_STEP)


def test_a_routed_free_end_is_judged_along_its_line_to_its_label(monkeypatch):
    """A free end routed out to free space is checked along the line it is drawn on, out to where its end label
    goes, not along the straight chord from its last nucleotide to the label: past a corner the chord cut an element
    the line keeps clear of (5LJ3's 12-nt end E:-16--5 was drawn as a fallback). 4YBB's 10-nt end, routed (no
    straight drawing offered) along a path fixed to turn between its last nucleotide and its label, with a small
    element planted on that chord."""
    import foldmap.na_ends as na_ends
    from foldmap.na_paths import _face, _seg_dist
    from foldmap.route import _DIRS, _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    tail = next(p for p in lay.na_paths if p.kind == "tail" and len(p.nts) == 10)
    (key,) = tail.between
    d, five = at[key], tail.nts[-1] < key[1]
    face, way = _face(d, key, not five)
    grid = _Grid(lay)
    start = grid.free_cell_along(face, way, d)
    heading = min(range(4), key=lambda k: -(_DIRS[k][0] * way[0] + _DIRS[k][1] * way[1]))

    def run(more=()):  # the end drawn along 19 cells straight out, then 13 to the side that is free
        for side in ((heading + 1) % 4, (heading + 3) % 4):
            cells = [(start[0] + _DIRS[heading][0] * k, start[1] + _DIRS[heading][1] * k) for k in range(20)]
            cells += [(cells[-1][0] + _DIRS[side][0] * k, cells[-1][1] + _DIRS[side][1] * k) for k in range(1, 14)]
            if all(grid.free(c) for c in cells):
                break
        rects, owners = _element_obstacles(lay)
        elements = (np.hstack([rects, np.array(more, float).reshape(-1, 4).T]), [*owners, *["small"] * len(more)])
        monkeypatch.setattr(na_ends, "_reach_out", lambda *args, **kw: list(cells))
        return na_ends._tail(_Router(_Grid(lay), elements), lay, d, key[0], key, tail.nts, five, None)

    monkeypatch.setattr(na_ends, "_TURNS", ())  # no straight drawing: the end is routed
    first, tip, _ = run()
    last = first.nt_xy[0 if five else -1]
    mid = np.add(last, tip) / 2  # the chord's middle, off the line by its corner ...
    assert not first.fallback and min(_seg_dist(mid, a, b) for a, b in zip(first.points, first.points[1:])) > 0.2
    small = (mid[0] - 0.05, mid[1] - 0.05, mid[0] + 0.05, mid[1] + 0.05)  # ... and an element there
    path, again, _ = run([small])
    assert not path.fallback and np.allclose(path.points, first.points) and np.allclose(again, tip)


def test_a_routed_five_prime_end_is_judged_along_its_line_too(monkeypatch):
    """A routed 5' end is kept out of other hairpin loops along its line from its stem out to its label, as a 3' end
    is; it was judged by [5'->3' points, label], whose last segment is a chord from the stem back to the label
    across the whole end (fix round 1). 4YBB's 7-nt 5' end, routed (no straight drawing offered)."""
    import foldmap.na_ends as na_ends
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    tail = next(p for p in lay.na_paths if p.kind == "tail" and len(p.nts) == 7)
    (key,) = tail.between
    d, five = at[key], tail.nts[-1] < key[1]
    assert five
    seen = []
    real = na_ends._enters
    monkeypatch.setattr(
        na_ends, "_enters", lambda pts, circles, own: seen.append(np.asarray(pts)) or real(pts, circles, own)
    )
    monkeypatch.setattr(na_ends, "_TURNS", ())  # no straight drawing: the end is routed
    path, tip, _ = na_ends._tail(
        _Router(_Grid(lay), _element_obstacles(lay)), lay, d, key[0], key, tail.nts, five, None
    )
    assert not path.fallback and seen
    line = seen[-1]
    assert np.allclose(line[0], d.point(key)) and np.allclose(line[-1], tip)  # stem first, label last
    for q in path.nt_xy:  # through every nucleotide, in order from the stem
        a, b = line[:-1], line[1:]
        t = np.clip(((q - a) * (b - a)).sum(1) / np.maximum(((b - a) ** 2).sum(1), 1e-12), 0, 1)
        assert np.hypot(*(a + (b - a) * t[:, None] - q).T).min() < 1e-6


def test_a_free_end_keeps_out_of_its_partners_doorway(monkeypatch):
    """A free end is not turned across the way its partner nucleotide's strand leaves the same end of the stem, when
    another drawing is clear (5F9R: the non-target strand's free end D:12-20, turned across the face of DNA 1, closed
    off the run from P1b, which had to cross it). 4YBB's 10-nt end AA:557-566, offered first a turn across the
    doorway."""
    import foldmap.na_ends as na_ends
    from foldmap.na_paths import _crosses_segment
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    tail = next(p for p in lay.na_paths if p.kind == "tail" and len(p.nts) == 10)
    (key,) = tail.between
    d, five = at[key], tail.nts[-1] < key[1]
    door = na_ends._doorway(d, key, five)
    turns = [t for t in (np.pi / 2, -np.pi / 2) if _crosses_segment(na_ends._free_end(d, key, 10, five, t)[0], *door)]
    assert len(turns) == 1  # the turn towards the partner runs across its doorway
    monkeypatch.setattr(na_ends, "_TURNS", (turns[0], 0.0))
    router = _Router(_Grid(lay), _element_obstacles(lay))
    path, _, _ = na_ends._tail(router, lay, d, key[0], key, tail.nts, five, None)
    assert not path.fallback and not _crosses_segment(path.points, *door)


@pytest.mark.parametrize("pdb, runs, protein", [("2GIS", 6, 0), ("4OO8", 4, 13)])
def test_runs_keep_off_the_strands_drawn_before_them(pdb, runs, protein):
    """A run with nucleotides on it is routed off every strand drawn before it, and every protein loop, unless that
    is a long detour (na_runs._DETOUR); failing that it pays for the cells beside them too (_BESIDE_DETOUR), else
    crossing them is only a cost; neither detour more than _DETOUR_STEPS steps per nucleotide longer. Measured with
    crossing them only a cost: run-run crossings 2GIS 8, 4OO8 5 (protein loops crossed 15). (30QM, pinned here while
    drawn in the plain frame, is drawn unrolled since the C4 revision.) 4OO8 protein 12 -> 13 when its helices were
    lettered (P2a/P2b: the longer names move the placement)."""
    from foldmap.nacheck import na_crossings

    lay, _, _, loops = _routed(pdb)
    m = na_crossings(lay, loops)
    assert m["na_run_crossings"] <= runs and m["na_protein_crossings"] <= protein


# two strands whose nucleotides sit within a bead of each other, held at their measured pairs (a held pair that clears
# must be removed): on 4YBB the runs P6 -> P7a and P20a2 -> P20b are drawn along one another through a corridor no
# route off the strands reaches, once short internal loops hold their helices together (Task 13b step 4, fix round 1);
# their spacings differ (1.43 and 2.05), so no slide clears every nucleotide. Plan 12 Task 2: with its branch strand
# AH:86-87 stacked beside its partner the layout no longer runs them along one another: none is crowded
CROWDED: dict[str, set] = {}


@pytest.mark.parametrize("pdb", [*FIXTURES, "2GIS", "2OEU"])
def test_no_nucleotide_sits_on_another_strands(pdb):
    """No nucleotide is drawn within a contact bead's width of a nucleotide of another strand or stem: where two runs
    cross or pass close (4OO8: the two sides of the P2-P3 internal loop cross where the stems are seen to), a run's
    nucleotides slide along it (4OO8's B:27 and B:42 were 0.25 apart, their beads on top of each other)."""
    from foldmap.na_paths import _BEAD

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

    def name(o):  # a stem by its id, a path by the stem nucleotides it runs between
        return o[1] if o[0] == "stem" else lay.na_paths[o[1]].between

    assert {frozenset({name(a), name(b)}) for a, b in close} == CROWDED.get(pdb, set())


def test_a_crowded_letter_looks_further_out_and_along_its_strand():
    """A base letter with no clear spot beside its strand looks further out and further along it before settling on
    the least bad spot (5B43: the 5' end B:-19--16 squeezed between its stem and a sheet; B:-17's letter landed
    0.15 from B:-18's). 1EHZ's acceptor end with an element planted over every nearer spot of its second letter."""
    from types import SimpleNamespace

    from foldmap.layout import NT_ROW
    from foldmap.na_labels import _LETTER_GAP, _place_letters
    from foldmap.na_paths import _near

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
    from foldmap.na_paths import NAPath

    return NAPath("run", 0, [1, 2], [(0.0, 0.0), (0.0, 0.0), (length, 0.0), (length, 0.0)], [])


@pytest.mark.parametrize("stage", ["first", "beside"])
def test_a_later_stage_takes_over_only_within_its_detour_caps(stage):
    """A staged drawing of a run replaces its plain one only within both caps: `bound` times as long (_DETOUR for
    the drawing off every strand and loop, _BESIDE_DETOUR for the one paying for the cells beside them) and at most
    _DETOUR_STEPS nucleotide steps longer per nucleotide and end (4YBB's P72->P71 went 5.7 times its direct distance
    round a cluster of stems)."""
    from foldmap.layout import NT_STEP
    from foldmap.na_paths import _length
    from foldmap.na_runs import _BESIDE_DETOUR, _DETOUR, _DETOUR_STEPS, _takes_over

    bound = _DETOUR if stage == "first" else _BESIDE_DETOUR
    m, plain = 2, _straight(20.0)
    extra = _DETOUR_STEPS * NT_STEP * (m + 2)  # 13.2: under (bound - 1) x 20 for both stages
    assert extra < (bound - 1) * _length(plain.points)
    assert _takes_over(_straight(20.0 + extra - 0.1), plain, bound, m, ())
    assert not _takes_over(_straight(20.0 + extra + 0.3), plain, bound, m, ())  # over the step cap
    short = _straight(4.0)  # here the ratio is the tighter cap
    assert not _takes_over(_straight(4.0 * bound + 0.3), short, bound, m, ())


def test_a_later_stage_never_replaces_a_clean_drawing_with_a_fallback():
    from foldmap.na_runs import _DETOUR, _takes_over

    plain, staged = _straight(10.0), _straight(10.5)
    staged.fallback = True
    assert not _takes_over(staged, plain, _DETOUR, 2, ())
    plain.fallback = True
    assert _takes_over(staged, plain, _DETOUR, 2, ())  # both fallbacks: the later stage may still be taken


def test_a_staged_strand_counts_the_router_debug_counts_of_its_drawing_only(monkeypatch):
    """router.tight and router.rejected count once per strand routed: the attempts a later stage discards (here
    the plain drawing and the first stage) leave no count behind."""
    import foldmap.na_runs as na_runs
    from foldmap.na_paths import _Strands
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

    monkeypatch.setattr(na_runs, "_drawing", drawing)
    strands = _Strands(router.grid)
    a, b = np.asarray(plain.points[1]), np.asarray(plain.points[2])
    mid, normal = (a + b) / 2, np.array([a[1] - b[1], b[0] - a[0]])
    strands.add([tuple(mid - normal), tuple(mid + normal)])  # a strand across the plain drawing
    extra = np.zeros_like(router.grid.used)
    got = na_runs._strand_between(router, s, i, j, at[(s, i)], at[(s, j)], None, (), strands, None, extra)
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
    from foldmap.na_paths import NATie

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
        x, y = text_centre(text)
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
// Runs the explorer's app script against stub DOM and viewer objects; drives viewer.on and records any camera move.
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
    this.plugin = { managers: { camera: { focusLoci: () => calls.push(['focusLoci']),
                                          focusSpheres: () => calls.push(['focusSpheres']) } },
                    canvas3d: { handleResize() {},
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
    for (const op of ['on']) {
      calls.length = 0;
      let error = null;
      const het = c.het === null ? undefined : c.het;
      try { await window.foldmapViewer[op](c.residues, het); } catch (e) { error = String(e); }
      await new Promise(r => setTimeout(r, 600));
      out.push({ case: c.name, op, error, calls: JSON.parse(JSON.stringify(calls)) });
    }
  }
  calls.length = 0;  // a click in Mol*'s own canvas goes through the camera manager: it must not move the camera
  const cam = window.foldmapViewer.plugin && window.foldmapViewer.plugin.plugin.managers.camera;
  if (cam) { cam.focusLoci(); cam.focusSpheres(); }
  out.push({ case: '-', op: 'camera', error: null, fly: typeof window.foldmapViewer.fly,
             calls: JSON.parse(JSON.stringify(calls)) });
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
    assert got[("protein", "on")]["calls"]  # something was sent to the viewer
    camera = got[("-", "camera")]  # the camera moves only when the reader moves it (the author: no zoom-in)
    assert camera["fly"] == "undefined" and camera["calls"] == []
    for r in got.values():
        assert not any(c[0] in ("zoomTo", "focus") for c in r["calls"]), r
    if viewer == "3dmol":
        on = got[("protein", "on")]["calls"]
        assert any(c[0] == "setStyle" and len(c[1].get("or", [])) == n for c in on)
        assert not any("stick" in c[2] and c[2]["stick"].get("color") for c in on)  # no nucleotide styling
        on = got[("mixed", "on")]["calls"]
        assert any(len(c[1].get("or", [])) == len(mixed["residues"]) for c in on if c[0] == key)
        assert sum(1 for c in on if c[0] == key and "resi" in c[1]) == len(mixed["het"])
    else:
        parts = got[("protein", "on")]["calls"][0][1]["data"]
        assert sum(p["end_auth_seq_id"] - p["beg_auth_seq_id"] + 1 for p in parts) >= n
        assert not any(p.get("representation") for p in parts)
        parts = got[("mixed", "on")]["calls"][0][1]["data"]
        assert sum(1 for p in parts if p.get("representation") == "ball-and-stick") == len(mixed["het"])


@pytest.mark.skipif(node is None, reason="node is not installed; the explorer's viewer code cannot be exercised")
@pytest.mark.parametrize("viewer", ["3dmol", "molstar"])
def test_explorer_viewers_take_a_nucleic_acid_only_selection(tmp_path, viewer):
    def cases(data):
        stem = data["links"]["nucleics"][0]
        return [{"name": "rna", "residues": [], "het": stem["het"]}]

    data, got = _drive(tmp_path, "1EHZ", viewer, cases)
    het = data["links"]["nucleics"][0]["het"]
    assert got[("rna", "on")]["error"] is None
    if viewer == "3dmol":
        styled = [c for c in got[("rna", "on")]["calls"] if c[0] == "setStyle" and "resi" in c[1]]
        assert len(styled) == len(het)
    else:
        assert len(got[("rna", "on")]["calls"][0][1]["data"]) == len(het)


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


def test_summary_of_a_nucleic_acid_alone_counts_its_nucleotides(capsys):
    from foldmap.cli import main

    assert main(["summary", str(_path("1EHZ"))]) == 0
    first = capsys.readouterr().out.splitlines()[0]
    assert first.endswith(": 76 nucleotides, 1 chain(s)")


@pytest.mark.parametrize("pdb", SMALL)
def test_stems_sit_nearer_the_elements_that_touch_them(pdb):
    from foldmap.nacheck import contacts

    lay, sses, _ = make_layout(_path(pdb))
    m = contacts(lay, sses)
    assert m["contact_auc"] is None or m["contact_auc"] >= 0.55  # IHF, wrapped by its arms, is the lowest


def test_the_contact_pull_moves_the_stems_of_30qm(monkeypatch):
    """Guard for layout._W_CONTACT: contact_auc does not respond to the pull (30QM scores 0.823 with it and 0.865
    without), so this displacement test is the guard. Without the pull 30QM's stems move by 8.9 page units at most
    (1RVA 2.8, 1IHF 0.9, 4LD0 0.7)."""
    import foldmap.layout as layout

    with_pull, _, _ = make_layout(_path("30QM"))
    monkeypatch.setattr(layout, "_W_CONTACT", 0.0)
    without, _, _ = make_layout(_path("30QM"))
    assert [d.id for d in with_pull.stems] == [d.id for d in without.stems]
    shift = max(np.hypot(a.cx - b.cx, a.cy - b.cy) for a, b in zip(with_pull.stems, without.stems, strict=True))
    assert shift > 3.0, shift


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
    alpha1 was placed on it (the trace a routing fallback). Each trace nucleotide is a fixed obstacle in placement
    for every element that does not touch the trace."""
    from foldmap.layout import _TRACE_HALF

    lay, _, _, _ = _routed("1QBJ")
    ((s, pts),) = lay.na_traces.items()
    touching = {r for r, v in lay.nucleic.contacts.items() if any(x[0] == s for x in v)}
    for x, y in np.asarray(pts).tolist():
        box = (x - _TRACE_HALF, y - _TRACE_HALF, x + _TRACE_HALF, y + _TRACE_HALF)
        on = [p.sse.id for p in lay.placed.values() if _overlap(p.rect, box)]
        on = [i for i in on if not touching & set(range(lay.placed[i].sse.start, lay.placed[i].sse.end + 1))]
        assert not on
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
    centres, want, _ = seen[0]
    assert len(centres) == 6  # one box per nucleotide of hexamer F
    for v in (centres, want):
        v = np.asarray(v) * 1e4
        assert np.all(np.abs(v - np.round(v)) < 1e-6)


def _gap_label_box(p, lay):
    from foldmap.na_labels import _em, _gap_box

    return _gap_box(p, _em(lay))


def test_a_gap_label_is_placed_clear_of_letters_beads_bodies_and_strands():
    """4YBB: "... 21 nt ..." sat at its gap's midpoint on element AC:26-46 and two runs (C10). A gap label is placed
    along its gap, either side, where its box keeps clear of bodies, strands, letters and nucleotides."""
    from foldmap.na_paths import crosses_body
    from foldmap.route import clip_len

    lay, _, _, loops = _routed("4YBB")
    labelled = [p for p in lay.na_paths if p.kind == "gap" and p.label]
    assert labelled
    letters = [*lay.na_letters.values(), *[d.point(k, 0.0) for d in lay.stems for k in d.axial]]
    segments = [(a, b) for q in lay.na_paths for a, b in zip(q.points, q.points[1:])]
    segments += [(a, b) for loop in loops for a, b in zip(loop.points, loop.points[1:])]
    for p in labelled:
        x0, y0, x1, y1 = box = _gap_label_box(p, lay)
        edges = [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]
        assert not any(x0 < x < x1 and y0 < y < y1 for x, y in letters)
        assert not any(clip_len(a, b, box) > 1e-9 for a, b in segments)
        bodies = [*lay.placed.values(), *lay.stems]
        assert not any(crosses_body(a, b, e) for e in bodies for a, b in edges)


def test_a_gap_label_moves_off_a_letter_planted_on_its_first_spot():
    from foldmap.na_labels import _label_gaps

    lay, _, _, loops = _routed("4YBB")
    p = next(p for p in lay.na_paths if p.kind == "gap" and p.label)
    first = p.label_xy
    lay.na_letters[(-1, -1)] = first  # a letter where the label went
    _label_gaps(lay, lay.na_paths, [loop.points for loop in loops])
    x0, y0, x1, y1 = _gap_label_box(p, lay)
    assert p.label_xy != first and not (x0 < first[0] < x1 and y0 < first[1] < y1)


def test_a_gap_with_no_way_on_the_grid_is_curved_before_it_is_a_fallback():
    """4PLX: the 3-nt chain break B:55-59 found no way on the grid and was drawn straight through PK2 (a
    fallback, an na_crossing). A gap the grid cannot route is curved along its ends, or bowed either side, when
    that keeps clear; only then straight."""
    from types import SimpleNamespace

    from foldmap.na_ends import _gap
    from foldmap.route import _element_obstacles, _Grid

    lay, _, _, _ = _routed("1URN")
    (p,) = [q for q in lay.na_paths if q.kind == "gap"]
    pts = np.asarray(p.points, float)
    ea, eb = pts[1] - pts[0], pts[-2] - pts[-1]
    ea, eb = ea / np.hypot(*ea), eb / np.hypot(*eb)
    (own,) = lay.stems  # both free ends leave the hairpin's one stem
    stuck = SimpleNamespace(route=lambda *a, **k: (None, None), elements=_element_obstacles(lay), grid=_Grid(lay))
    drawn = _gap(stuck, lay.nucleic.gaps[0], (pts[0], ea, own), (pts[-1], -eb, own))
    assert not drawn.fallback and len(drawn.points) > 2


def test_a_gap_label_is_drawn_inside_the_box_placement_kept_clear():
    """The box _label_gaps keeps clear holds the drawn text, and not much more (2NZ4: a box 2.6 times too wide found
    no clear spot between P2 and its gap and put "... 6 nt ..." on the stem)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    from foldmap.na_labels import _em, _gap_box
    from foldmap.render import draw

    lay, sses, _, loops = _routed("4YBB")
    fig = draw(lay, loops, sses)
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    to_page = fig.axes[0].transData.inverted()
    texts = _gid(fig, "na:gap:")
    assert texts
    for t in texts:
        x0, y0, x1, y1 = to_page.transform_bbox(t.get_window_extent(canvas.get_renderer())).extents
        b = _gap_box(lay.na_paths[int(t.get_gid().rsplit(":", 1)[1])], _em(lay))
        assert b[0] <= x0 and x1 <= b[2] and b[1] <= y0 and y1 <= b[3]
        assert (b[2] - b[0]) < 1.6 * (x1 - x0)


def test_a_single_strand_wrapped_over_its_protein_keeps_the_protein_under_it():
    """1EYG: the ssDNA wraps over the SSB tetramer. As one bounding box the traces pushed the tetramer apart and
    below them (touching elements up to 13.4 from the nearest nucleotide, two dimers split); elements touching a
    trace are exempt from it, so every touching element stays within reach of the strand it holds."""
    lay, _, _, _ = _routed("1EYG")
    pts = np.vstack([np.asarray(v) for v in lay.na_traces.values()])
    traced = set(lay.na_traces)
    far = []
    for p in lay.placed.values():
        if not any(x[0] in traced for r in range(p.sse.start, p.sse.end + 1) for x in lay.nucleic.contacts.get(r, ())):
            continue
        x0, y0, x1, y1 = p.rect
        gap = np.hypot(
            np.clip(np.maximum(x0 - pts[:, 0], pts[:, 0] - x1), 0, None),
            np.clip(np.maximum(y0 - pts[:, 1], pts[:, 1] - y1), 0, None),
        )
        if gap.min() > 6.0:
            far.append((p.sse.id, round(float(gap.min()), 1)))
    assert not far
    # and the traces lying on the protein they touch are no routing fallbacks (na_ends._trace_clear)
    assert not [p.strand for p in lay.na_paths if p.kind == "trace" and p.fallback]


def test_a_ghost_of_a_strand_holding_a_trace_is_that_strand():
    """A ghost is the strand it copies: where that strand holds a single strand's nucleotides, the trace may lie
    on the ghost as on the strand. 1EYG: B:2051-2060 holds nucleotides 1-5 and 8 of strand 1; a ghost of it laid
    on that trace at nucleotide 3 is clear. A ghost of an element holding none of it still blocks."""
    from dataclasses import replace

    def fallbacks(ghost_of):
        lay, sses, bb = make_layout(_path("1EYG"))
        x, y = np.asarray(lay.na_traces[1])[3, :2]
        lay.ghosts.append(replace(lay.placed[ghost_of], cx=float(x), cy=float(y), ghost=True))
        route_loops(lay, sses, bb)
        return [q.strand for q in lay.na_paths if q.kind == "trace" and q.fallback]

    assert fallbacks("B:2051-2060") == []
    assert 1 in fallbacks("A:1061-1070")  # this helix, in no sheet, holds nothing of strand 1


def test_a_chain_break_with_one_end_not_drawn_is_a_short_break_mark(monkeypatch):
    """1URN with the piece after the break cut to its last nucleotide (no stem, too short for a trace): it has no
    end to route to, so the break is marked by a short dashed stub at the end that is drawn, not left unmarked."""
    import foldmap.cli as cli
    from foldmap.model import _cut

    monkeypatch.setattr(cli, "_na_extent", lambda na, extent: (_cut(na, {0: [(0, 11)], 1: [(5, 5)]}), False))
    lay, _, _, _ = _routed("1URN")
    assert lay.nucleic.gaps and 0 in lay.na_traces
    (gap,) = [p for p in lay.na_paths if p.kind == "gap"]
    assert len(gap.points) == 2 and not gap.fallback


def _piece(cx, cy, angle, length=4.0, width=2.0, ident="p"):
    from types import SimpleNamespace

    return SimpleNamespace(cx=cx, cy=cy, angle=angle, length=length, width=width, id=ident)


def test_a_straight_joint_along_its_own_piece_turns_an_elbow_round_it():
    """4YBB's P11b2 -> P11b1: the straight line from one piece's face to the next one's ran along the second piece's
    body. A joint whose straight line runs into its own piece deeper than a bead's radius is drawn out of the one
    face and into the other round a corner (na_runs._elbow), clear of both; a shallow graze of a corner is kept."""
    from types import SimpleNamespace

    from foldmap.na_paths import _BEAD, graze
    from foldmap.na_runs import _elbow, _grazes

    da, db = _piece(0.0, 0.0, 0.0, ident="a"), _piece(3.5, 0.5, np.pi / 2, ident="b")  # b stands right of a's face
    fa, ea, fb, eb = (2.0, 0.5), (1.0, 0.0), (4.0, 2.5), (0.0, 1.0)  # out of a's end, into b's top
    assert _grazes([fa, fb], (da, db))  # straight, it cuts through b
    router = SimpleNamespace(elements=(np.empty((4, 0)), []))
    pts = _elbow(router, da, db, fa, ea, fb, eb)
    assert pts is not None and pts[0] == fa and pts[-1] == fb and len(pts) > 2
    assert not _grazes(pts, (da, db))
    corner = _piece(0.0, 0.0, 0.0)  # a line clipping a corner by less than a bead's radius is no graze
    assert 0 < graze((1.7, -1.2), (2.2, -0.7), corner) < _BEAD / 2 and not _grazes(
        [(1.7, -1.2), (2.2, -0.7)], (corner,)
    )


def test_a_shallow_graze_at_a_path_end_is_not_counted_as_through_a_stem():
    from foldmap.nacheck import _through_stem

    d = _piece(0.0, 0.0, 0.0)
    shallow = np.array([[1.7, -1.2, 2.2, -0.7]])  # one segment (first and last), clipping a corner by 0.1
    assert not _through_stem(shallow, d)
    deep = np.array([[-3.0, 0.0, 3.0, 0.0]])
    assert _through_stem(deep, d)
    middle = np.array([[-3.0, -3.0, 1.7, -1.2], [1.7, -1.2, 2.2, -0.7], [2.2, -0.7, 4.0, -3.0]])
    assert _through_stem(middle, d)  # inside a path, even a graze counts


def test_a_trace_is_exempt_from_an_element_only_near_the_nucleotides_it_touches():
    """na_ends._trace_clear: an element touching a trace may lie under it where it holds it, not anywhere along
    it (an element touching the 5' end, crossed by the trace twenty nucleotides on, makes a fallback)."""
    from types import SimpleNamespace

    from foldmap.layout import Placed
    from foldmap.na_ends import _trace_clear
    from helpers import fake_sse

    def check(touched: int, crossed: str = "H") -> bool:
        """The trace's nucleotides 14-16 cross whichever of the helix (H), the strand beside it in a sheet (E) or a
        ghost of the helix (ghost) sits at y = 0; the others sit clear of it at y = 6."""
        h, e = fake_sse("H", 5, 5), fake_sse("E", 6, 6)
        y = lambda k: 0.0 if crossed == k else 6.0  # noqa: E731
        layout = SimpleNamespace(
            nucleic=SimpleNamespace(contacts={5: {(0, touched)}}),
            placed={h.id: Placed(h, 15.0, y("H"), 1.6, 0.0, 2.0), e.id: Placed(e, 15.0, y("E"), 1.6, 0.0, 2.0)},
            ghosts=[Placed(h, 15.0, y("ghost"), 1.6, 0.0, 2.0, ghost=True)],
            stems=[],
            sheets=[[h.id, e.id]] if crossed == "E" else [],
        )
        return _trace_clear(layout, 0, [(float(x), 0.0) for x in range(21)])

    assert check(15) and check(13)  # it holds the nucleotides there (or a neighbour's)
    assert not check(0)  # it holds the 5' end only: the trace crossing it far on is not clear
    assert check(15, "E") and not check(0, "E")  # its sheet likewise
    assert check(15, "ghost") and not check(0, "ghost")  # and a ghost of it


def test_trace_nucleotides_drawn_on_one_spot_are_nudged_apart():
    """1EYG: two single strands over the SSB tetramer project two nucleotides onto one spot, their contact beads on
    each other; Q:16-17 runs toward the viewer, its steps shorter than a bead (bead_overlaps 1). The later
    nucleotide moves sideways off its trace, as little as clears a bead's width."""
    from types import SimpleNamespace

    from foldmap.frame import Frame
    from foldmap.layout import SCALE, _traces
    from foldmap.model import NAStrand, ResLabel
    from foldmap.na_paths import _BEAD

    def strand(chain, xyz):
        labels = [ResLabel(chain, k + 1, "", "DT") for k in range(len(xyz))]
        p = np.asarray(xyz, float)
        return NAStrand(chain, labels, p, p, False)

    a = strand("Q", [[x * 6.0, 0.0, 0.0] for x in range(5)])  # along x
    b = strand("R", [[12.0, y * 6.0 - 12.0, 3.0] for y in range(5)])  # along y, its middle on a's (12, 0)
    na = SimpleNamespace(strands=[a, b], structure=SimpleNamespace(stems=[]))
    frame = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    tr = _traces(na, frame)
    assert np.allclose(tr[0], a.p[:, :2] * SCALE)  # the first strand stays where it is seen
    pts = np.vstack(list(tr.values()))
    d = np.hypot(*(pts[:, None] - pts[None]).transpose(2, 0, 1))
    assert d[np.triu_indices(len(pts), 1)].min() >= _BEAD - 1e-3  # on the grid
    moved = tr[1][2] - b.p[2, :2] * SCALE
    assert abs(moved[1]) < 1e-3 and _BEAD - 1e-3 <= abs(moved[0]) <= _BEAD + 0.05  # sideways, by a bead
    assert np.allclose(np.delete(tr[1], 2, axis=0), np.delete(b.p[:, :2], 2, axis=0) * SCALE)
    end_on = strand("S", [[30.0, 30.0, z * 6.0] for z in range(4)])  # seen end-on: one spot, steps 1.5 Å apart
    end_on.p[:, 0] += np.arange(4) * 1.5
    pts = _traces(SimpleNamespace(strands=[end_on], structure=SimpleNamespace(stems=[])), frame)[0]
    d = np.hypot(*(pts[:, None] - pts[None]).transpose(2, 0, 1))
    assert d[np.triu_indices(len(pts), 1)].min() >= _BEAD - 1e-3


def test_a_free_end_is_never_drawn_into_another_stems_hairpin_loop(monkeypatch):
    """4YBB (step 4 campaign): a free end was drawn inside another stem's hairpin circle, and the gap leaving it
    crossed the loop. A free end's drawing keeps out of every other stem's hairpin circle when another drawing does.
    4YBB's 10-nt end AA:557-566 with a hairpin circle planted over its straight drawing."""
    import foldmap.na_ends as na_ends
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _, _ = _routed("4YBB")
    at = {key: d for d in lay.stems for key in d.axial}
    tail = next(p for p in lay.na_paths if p.kind == "tail" and len(p.nts) == 10)
    (key,) = tail.between
    d, five = at[key], tail.nts[-1] < key[1]
    straight = np.asarray(na_ends._free_end(d, key, 10, five, 0.0)[1])
    other = next(s for s in lay.stems if s is not d)
    circle = (straight.mean(axis=0), 1.0, other)
    assert na_ends._enters(straight, [circle], d)
    monkeypatch.setattr(na_ends, "_TURNS", (0.0, np.pi / 2, -np.pi / 2))
    router = _Router(_Grid(lay), _element_obstacles(lay))
    path, _, _ = na_ends._tail(router, lay, d, key[0], key, tail.nts, five, None, circles=[circle])
    assert not path.fallback and not na_ends._enters(path.points, [circle], d)
    assert not na_ends._enters(straight, [(circle[0], 1.0, d)], d)  # its own stem's loop is not another's


def test_a_free_end_the_protein_holds_lies_toward_it():
    """1URN: U1A's RRM holds the hairpin loop R:6-12, drawn as a free end; it was drawn straight on, away from the
    sheet (R:11 5.4 from its nearest touching strand against a median element gap of 5.1). A touched free end tries
    the turns that bring it nearest its touching elements first: each loop nucleotide a strand of the RRM touches
    is drawn nearer that strand than the median element."""
    from foldmap.nacheck import _gap, ss_drawn

    lay, sses, _, _ = _routed("1URN")
    na, where = lay.nucleic, ss_drawn(lay)
    drawn = [s for s in sses if s.id in lay.placed]
    checked = 0
    for x in sorted(where):
        strands = [
            s for s in drawn if s.kind == "E" and any(x in na.contacts.get(r, ()) for r in range(s.start, s.end + 1))
        ]
        if not strands:
            continue
        pt = (*where[x], *where[x])
        gaps = [_gap(lay.placed[s.id].rect, pt) for s in drawn]
        assert min(_gap(lay.placed[s.id].rect, pt) for s in strands) < np.median(gaps), x
        checked += 1
    assert checked >= 4


def test_a_trace_seen_straight_end_on_is_nudged_apart_too():
    """Three nucleotides projected onto one spot: the middle one's neighbours coincide, so its trace has no
    direction there; it is nudged off its 5' neighbour's way, or along a fixed axis."""
    from types import SimpleNamespace

    from foldmap.frame import Frame
    from foldmap.layout import _traces
    from foldmap.model import NAStrand, ResLabel
    from foldmap.na_paths import _BEAD

    p = np.array([[0.0, 0.0, z * 6.0] for z in range(3)])
    strand = NAStrand("S", [ResLabel("S", k + 1, "", "DT") for k in range(3)], p, p, False)
    frame = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    pts = _traces(SimpleNamespace(strands=[strand], structure=SimpleNamespace(stems=[])), frame)[0]
    d = np.hypot(*(pts[:, None] - pts[None]).transpose(2, 0, 1))
    assert d[np.triu_indices(3, 1)].min() >= _BEAD - 1e-3


def test_drawn_base_letters_are_the_parent_bases():
    # 1EHZ: tRNA with 5MC, YYG, PSU, ... drawn as their parent base, not the last character of the residue name
    fig, _ = make_figure_and_loops(_path("1EHZ"))
    drawn = {t.get_gid(): t.get_text() for t in _gid(fig, "na:nt:")}
    assert len(drawn) == 76
    assert drawn["na:nt:A:40"] == "C"  # 5MC40
    assert drawn["na:nt:A:37"] == "G"  # YYG37 (wybutosine)
    assert set(drawn.values()) <= set("ACGU")


def test_a_block_bars_first_pair_names_the_partners_chain_when_it_is_another_chain():
    from foldmap.style import Style

    fig, _ = make_figure_and_loops(_path("1IHF"), look=Style(na_style="block"))
    lay, _, _ = make_layout(_path("1IHF"), look=Style(na_style="block"))
    firsts = {a.get_gid(): a.get_text() for a in _gid(fig, "na:stem:") if a.get_gid().endswith(":first")}
    strands = lay.nucleic.strands
    across = 0
    for k, d in enumerate(lay.stems):
        a, i, b, j = d.pairs[0]
        own, other = strands[a].labels[i], strands[b].labels[j]
        if other.chain == own.chain:
            assert firsts[f"na:stem:{k}:first"] == f"{own.seq}·{other.seq}"
        else:
            across += 1
            assert firsts[f"na:stem:{k}:first"] == f"{own.seq}·{other.chain}:{other.seq}"
    assert across


def test_only_traces_of_three_nucleotides_or_more_are_walls():
    from types import SimpleNamespace

    from foldmap.layout import _trace_walls

    na = SimpleNamespace(contacts={})
    short = np.array([[0.0, 0.0], [1.0, 0.0]])
    long = np.array([[5.0, 0.0], [6.0, 0.0], [7.0, 0.0]])
    assert len(_trace_walls([], {0: short}, na).boxes) == 0
    assert len(_trace_walls([], {0: short, 1: long}, na).boxes) == 3


def test_a_gap_is_named_as_a_gap_in_warnings():
    """A fallback dashed gap is warned about as "gap after <chain>:<last>", not as a one-nucleotide strand (5LJ3 read
    "strand I:10-10" for the gap skipping I:11-53)."""
    from types import SimpleNamespace

    from foldmap.na_paths import NAPath, path_name

    labels = [SimpleNamespace(chain="I", tag=str(n)) for n in range(1, 11)]
    na = SimpleNamespace(strands=[SimpleNamespace(labels=labels)])
    assert path_name(NAPath("gap", 0, [], [], []), na) == "gap after I:10"
    assert path_name(NAPath("tail", 0, [8, 9], [], []), na) == "strand I:9-10"


@pytest.mark.parametrize("pdb", ["2GIS", "30QM", "1EHZ", "4OO8", "1KX5", "3CRX"])
def test_stem_labels_sit_by_their_own_stem_clear_of_other_text(pdb):
    from foldmap import figcheck
    from foldmap.layout import stem_label_pos
    from foldmap.spec import draw_spec, figure_spec

    lay, sses, bb = make_layout(NA / f"{pdb}.cif.gz")
    route_loops(lay, sses, bb)

    def body(p, d):
        v = np.subtract(p, (d.cx, d.cy))
        a, b = abs(v @ d.direction) - d.length / 2, abs(v @ d.normal) - d.width / 2
        return float(np.hypot(max(a, 0), max(b, 0)))

    for k, d in enumerate(lay.stems):
        p = stem_label_pos(d)
        dist = [body(p, e) for e in lay.stems]
        assert dist[k] == min(dist) and dist[k] <= d.width, (pdb, d.name)
    d = draw_spec(figure_spec(NA / f"{pdb}.cif.gz"))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    page = figcheck.measure(d)  # ink boxes (fixed hinting), pixels
    for t, g, box in page.texts:
        if g.startswith("na:label:"):  # 1KX5 labels 6, 9, 12, 15 did; a name 0.5 em from its letters (plan 10
            # Task 4) meets their line boxes (get_window_extent: line height and descent) but not their ink, nor
            # their page boxes in the frame they are turned in alike (Page.pair_boxes)
            hits = [h for u, h, b in page.texts if h != g and Bbox.overlaps(*page.pair_boxes(t, box, u, b))]
            assert not hits, (pdb, g, hits)


@pytest.mark.parametrize("pdb", ["1EHZ", "2GIS"])
def test_tie_letters_are_base_letter_size_and_explained(pdb):
    from foldmap.cli import make_figure

    ax = make_figure(NA / f"{pdb}.cif.gz").axes[0]
    ties = [t for t in ax.texts if (t.get_gid() or "").endswith(":label") and (t.get_gid() or "").startswith("na:tie:")]
    letters = [t for t in ax.texts if (t.get_gid() or "").startswith("na:nt:")]
    assert ties and {t.get_fontsize() for t in ties} == {letters[0].get_fontsize()}
    assert any("tertiary base pair" in t.get_text() for t in ax.texts)
    assert [ln for ln in ax.lines if ln.get_gid() == "legend-tie"]  # drawn as the two ticks


def test_element_labels_keep_off_nucleic_acid_traces():
    """1EYG: chain A's α1 name sat under an ssDNA trace's bead and line; an element's name takes its other side
    when a trace, its beads or its letters lie on the first."""
    from foldmap.layout import _GLYPH_H, _GLYPH_W
    from foldmap.na_paths import _BEAD
    from foldmap.route import clip_len

    lay, sses, _, _ = _routed("1EYG")
    beads = [xy for p in lay.na_paths for xy in p.nt_xy]
    segs = [(a, b) for p in lay.na_paths for a, b in zip(p.points, p.points[1:])]
    r, gw, gh = _BEAD / 2, _GLYPH_W / 2, _GLYPH_H / 2
    boxes = [(n, b) for n, b in label_boxes(lay, sses) if n.startswith("label:")]
    assert boxes and segs
    for n, box in boxes:
        x0, y0, x1, y1 = box
        assert not [q for q in beads if x0 - r < q[0] < x1 + r and y0 - r < q[1] < y1 + r], n
        assert not [s for s in segs if clip_len(*s, box) > 1e-9], n
        assert not [q for q in lay.na_letters.values() if x0 - gw < q[0] < x1 + gw and y0 - gh < q[1] < y1 + gh], n


def test_stem_name_slots_do_not_depend_on_the_installed_fonts(monkeypatch):
    """Names are measured from built-in Liberation Sans Bold metrics: on a DejaVu-only machine (a bare container)
    7 of 1KX5's 22 slots and its bounds changed while the layout measured with the fonts it found."""
    from matplotlib import font_manager

    lay, _, _ = make_layout(_path("1KX5"))
    want = ([(d.label_side, d.label_out) for d in lay.stems], lay.bounds)
    dejavu = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", weight="bold"))
    monkeypatch.setattr(font_manager, "findfont", lambda *a, **k: dejavu)
    monkeypatch.setattr(font_manager.fontManager, "findfont", lambda *a, **k: dejavu)
    monkeypatch.setattr(font_manager.fontManager, "_find_fonts_by_props", lambda *a, **k: [dejavu])
    lay, _, _ = make_layout(_path("1KX5"))
    assert ([(d.label_side, d.label_out) for d in lay.stems], lay.bounds) == want


@pytest.mark.parametrize("path", [NA / "30QM.cif.gz", Path(__file__).parent / "data" / "1UBQ.cif"])
def test_no_tie_legend_without_lettered_ties(path):
    """The legend explains tie letters only when a tie is drawn as two lettered marks: 30QM has none, 1UBQ no
    nucleic acid at all."""
    from foldmap.cli import make_figure

    ax = make_figure(path).axes[0]
    assert not [t for t in ax.texts if (t.get_gid() or "").startswith("na:tie:") and t.get_gid().endswith(":label")]
    assert not [t for t in ax.texts if "tertiary base pair" in t.get_text()]
    assert not [a for a in [*ax.lines, *ax.texts] if (a.get_gid() or "").startswith("legend-tie")]


_TRNA_NAMES = Path(__file__).parents[1] / "docs" / "names" / "1EHZ-trna.yaml"


@pytest.mark.parametrize(
    "pdb, names",
    [("1EHZ", None), ("1EHZ", _TRNA_NAMES), ("2GIS", None), ("30QM", None), ("1KX5", None)],
    ids=["1EHZ", "1EHZ-trna-names", "2GIS", "30QM", "1KX5"],
)
def test_stem_name_sits_by_its_own_letters(pdb, names):
    """Each stem's name sits 1.0 em (layout._INK_MARGIN) beyond its own letter row, across the stem, centred along
    it: 1EHZ's D name sat six glyph heights above its row, T between the T arm and the D-loop run; at 0.25 em (plan
    10 Task 4) D and T read as letters of their rows, at 0.5 em T still stacked on a G of the T arm (fix round 1).
    Measured by figcheck.stem_name_gaps in the stem's frame under its fixed hinting (Liberation Sans only); the limit
    is that margin in page units plus 1e-3 (figcheck.stem_gap_tolerance, R7), and at least 0.9 of 1.0 em of the
    drawn text (the layout's em, from the page before routing and names widen it, is up to 6 % smaller than the
    drawn page's: 30QM and 1KX5 measure 0.95 em). The name is also nearer its own letters than any other stem (fix
    round 2 ruling: at 1.0 em 1KX5's DNA 1h is 0.53 from its letters and 1.05 from the next piece, past the T4
    brief's half), never on its letters, and the letters stay apart (nacheck letter_overlaps)."""
    from foldmap import figcheck
    from foldmap.layout import _INK_MARGIN
    from foldmap.nacheck import render_checks
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(_path(pdb), names=[str(names)] if names else None))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    limit = figcheck.stem_gap_tolerance(_INK_MARGIN * _page_em(d))  # the margin in page units, plus 1e-3
    gaps = figcheck.stem_name_gaps(d)
    assert gaps and len(gaps) == len([s for s in d.layout.stems if not s.ghost])
    far = {k: (round(own, 3), round(other, 3)) for k, (own, other) in gaps.items() if not 0 < own <= limit}
    assert not far, (pdb, round(limit, 3), far)
    em = _page_em(d)
    tight = {k: round(own / em, 3) for k, (own, _) in gaps.items() if own < 0.9 * 1.0 * em}  # in ems
    assert not tight, (pdb, tight)
    near_other = {k: (own, other) for k, (own, other) in gaps.items() if not own < other}
    assert not near_other, (pdb, near_other)
    assert render_checks(d.layout, d.loops, d.figure, d.sses)["letter_overlaps"] == 0
    for s in d.layout.stems:  # centred along its row: midway between the row's first and last letters
        row = [s.axial[key] for key in s.axial if s.side[key] == s.label_side]
        assert s.label_along == pytest.approx((min(row) + max(row)) / 2, abs=1e-4), (pdb, s.name)


@pytest.mark.parametrize(
    "pdb, names",
    [("1EHZ", None), ("1EHZ", _TRNA_NAMES), ("2GIS", None), ("30QM", None), ("1KX5", None), ("4OO8", None)],
    ids=["1EHZ", "1EHZ-trna-names", "2GIS", "30QM", "1KX5", "4OO8"],
)
def test_routes_keep_the_margin_from_stem_names(pdb, names):
    """No drawn line (run, joint, tail, gap, hairpin loop, protein loop, tie) passes within layout._INK_MARGIN em of
    a stem's name: routes keep out of its ink box grown by the margin (layout.stem_keep_box). 30QM's P4 had a run
    0.32 em from its ink (plan 10 Task 4 fix round 1). Measured from the name's ink box to each line's stroke,
    segment by segment (figcheck.stem_name_line_gaps), in ems of the drawn text."""
    from foldmap import figcheck
    from foldmap.layout import _INK_MARGIN
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(_path(pdb), names=[str(names)] if names else None))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    em = _page_em(d)
    gaps = figcheck.stem_name_line_gaps(d)
    assert gaps
    near = {k: (round(gap / em, 3), line) for k, (gap, line) in gaps.items() if gap < _INK_MARGIN * em}
    assert not near, (pdb, near)


def test_block_slots_hold_no_name_so_routes_keep_no_margin_from_them():
    """In block style (4YBB, cropped) a stem's name is drawn on its block and its side slot holds the first pair's
    numbers: routes keep out of that text's box only, not a name's margin (layout.stem_keep_box), which raised the
    NA campaign's 4YBB crossings from 110 to 142 (plan 10 Task 4 fix round 2)."""
    from foldmap.layout import _stem_label_box, stem_keep_box

    lay, sses, bb = make_layout(NA / "4YBB.cif.gz")
    assert lay.na_cropped and lay.stems
    assert all(stem_keep_box(d) == _stem_label_box(d) for d in lay.stems)
    lay, sses, bb = make_layout(NA / "1EHZ.cif.gz")  # letters drawn: the name's ink and its margin
    assert all(d.label_keep is not None and stem_keep_box(d) != _stem_label_box(d) for d in lay.stems)


def test_a_stem_name_stays_nearer_its_own_stem_unless_no_slot_allows_it():
    """_pick_slot takes the first slot that is clear, nearer the name's own stem than any other and has the margin
    of room; else (a fallback) the most room among the clear slots nearer its own stem, and only when no clear slot
    is nearer its own stem the most room among the clear ones (plan 10 Task 4 fix round 2 ruling)."""
    from foldmap.layout import _pick_slot

    assert _pick_slot([(True, True, 0.2), (True, True, 1.5), (True, True, 9.0)], 1.0) == (1, False)
    assert _pick_slot([(True, False, 9.0), (True, True, 0.3), (True, True, 0.6)], 1.0) == (2, True)
    assert _pick_slot([(True, False, 0.3), (True, False, 0.9), (False, True, 9.0)], 1.0) == (1, True)
    assert _pick_slot([(False, True, 2.0), (True, False, 0.1)], 1.0) == (1, True)
    assert _pick_slot([(True, True, 0.5), (True, True, 0.5)], 1.0) == (0, True)  # the first of equals


@pytest.mark.parametrize("pdb", ["1EHZ", "1KX5"])
def test_stem_names_and_letters_are_drawn_where_the_layout_measures_them(pdb):
    """The layout places a name by its ink (layout._ink_line: built-in Liberation Sans hulls and advances) and render
    starts the text where that ink lands (layout.text_anchor), so the drawn ink box (figcheck, fixed hinting) is the
    predicted one, within 5e-3 page units, and the hulls lie inside the ink (to 2e-3: matplotlib rounds kerning to
    whole pixels at its outline size, which moves the end of a kerned name such as "DNA 1a" 1.2e-3 units on 3.10).
    1KX5's names hold a space and a kerned pair ("DNA 1a"), where matplotlib 3.10's centring drifted 0.05 units."""
    from foldmap import figcheck
    from foldmap.layout import _ink_line, _ink_points, _reach, _turn
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(_path(pdb)))
    if reason := figcheck.ink_font_reason(d.figure):
        pytest.skip(reason)
    ax, em = d.figure.axes[0], _page_em(d)
    to_page = ax.transData.inverted()
    drawn = {name: (t, to_page.transform_bbox(box).extents) for t, name, box in figcheck.measure(d).texts}
    for k, s in enumerate(d.layout.stems):
        t, got = drawn[f"na:label:{k}"]
        want = _reach(_ink_points(s.name, True, 0.8 * em, s.angle)) + np.tile(text_centre(t), 2)
        assert np.all((got - want) * [-1, -1, 1, 1] > -2e-3) and np.abs(got - want).max() < 5e-3, (s.name, got - want)
    letters = [(t, box) for name, (t, box) in drawn.items() if name.startswith("na:nt:")]
    assert letters
    for t, got in letters:
        ink = _turn(_ink_line(t.get_text(), t.get_fontweight() == "bold")[0] * 0.75 * em, t.get_rotation())
        want = _reach(ink) + np.tile(text_centre(t), 2)
        assert np.all((got - want) * [-1, -1, 1, 1] > -2e-3) and np.abs(got - want).max() < 5e-3, (
            t.get_gid(),
            got - want,
        )


def test_single_stranded_letters_keep_off_the_strands():
    """C14 (random-100): a single-stranded nucleotide's letter is placed after every strand is routed, but its
    spot was tested against letters, nucleotides, names and bodies only, so a run or a cropped gap could pass
    through it (3OWW: A:87 under the dashed gap; 8SA2 had 14 such letters). Spots a strand passes through are not
    clear."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(NA.parent / "r100" / "3OWW-AB.cif.gz"))
    assert d.layout.na_letters
    hits = [h for h in figcheck.text_collisions(d) if h.startswith("na:nt:") and "× na:run:" in h]
    assert hits == [], hits


def test_a_hairpin_letter_struck_by_the_clearance_keeps_an_inside_spot_clear_of_ink(monkeypatch):
    """Batch B fix rounds 1-2 (8SA2 A:35, A:36): a hairpin loop's letter whose spots inside its stem's footprint are
    all within the line clearance of a strand takes an inside spot no strand runs through (clear of the ink), before
    the roomiest inside spot or any spot outside the footprint. 1EHZ with the clearance made huge (every spot "near"
    a strand) and a short strand laid through where one hairpin letter would go: it moves to another inside spot
    clear of the ink."""
    from foldmap import na_labels
    from foldmap.layout import _GLYPH_H, _LETTER_EM, _footprint, _ink_points, _reach
    from foldmap.na_paths import NAPath
    from foldmap.route import clip_lengths

    lay, _, _, _ = _routed("1EHZ")
    monkeypatch.setattr(na_labels, "line_clearance", lambda layout: 5.0)
    first = na_labels._place_letters(lay, lay.na_paths)
    at = {key: d for d in lay.stems for key in d.axial}
    arc = next(p for p in lay.na_paths if p.kind == "hairpin")
    key = (arc.strand, arc.nts[len(arc.nts) // 2])
    x, y = first[key]
    bar = NAPath("run", arc.strand, [], [(x - 0.3, y), (x + 0.3, y)], [])  # a strand through that letter
    letters = na_labels._place_letters(lay, [*lay.na_paths, bar])
    d = at[arc.between[0]]
    x0, y0, x1, y1 = _footprint(d)
    g = _GLYPH_H / 2
    lx, ly = letters[key]
    assert (lx, ly) != (x, y)
    assert x0 <= lx - d.cx - g and lx - d.cx + g <= x1 and y0 <= ly - d.cy - g and ly - d.cy + g <= y1
    base = lay.nucleic.strands[key[0]].bases[key[1]]
    box = _reach(_ink_points(base, False, _LETTER_EM * na_labels._em(lay), 0.0)) + np.array([lx, ly, lx, ly])
    segs = [(*a, *b) for p in [*lay.na_paths, bar] for a, b in zip(p.points, p.points[1:])]
    assert not (clip_lengths(segs, [box]) > 0).any()


def test_a_free_end_tries_every_drawing_off_the_strands_before_one_crossing_them(monkeypatch):
    """Batch C2a review Q10: _tail's first pass keeps off the strands drawn so far (9CJI); in it a routed end is tried
    only strictly (kept off their cells); one crossing strands as a cost comes only in the old order after, never
    before a straight drawing that merely runs near a strand. With every straight drawing near a strand and no route
    found, the first pass asks for strict routes only (before the fix it also asked for crossing ones)."""
    from foldmap import na_ends
    from foldmap.na_paths import _Strands

    calls = []

    def routed(*args):
        calls.append(args[-1])  # strict
        return None

    monkeypatch.setattr(na_ends, "_routed_end", routed)
    monkeypatch.setattr(_Strands, "near", lambda self, points, r: True)
    lay, sses, bb = make_layout(NA.parent / "r100" / "9CJI-AB.cif.gz")
    route_loops(lay, sses, bb)
    # the 7-nt 3' end: strict routes with and without its partner's doorway (label room or not), then the old order
    # finds a straight drawing (near no longer asked): no crossing route tried
    assert calls == [True] * 4


@pytest.mark.parametrize("path", [NA / "30QM.cif.gz", NA.parent / "r100" / "9CJI-AB.cif.gz"], ids=["30QM", "9CJI"])
def test_each_end_label_names_its_own_end(path):
    """9CJI (batch C2a review Q6): the crRNA's 3' end, drawn after the 5' label was placed, started 0.6 from that label
    and 1.1 from the label's own end, so "5′" read as naming it; 30QM's 3' label lay 0.61 from a run's nucleotide.
    Once the free ends are drawn, and again once the runs are, a free end's label that lies nearer another drawn
    nucleotide than its own end moves round its end to a spot where its own is nearest (route_na._own_ends)."""
    lay, sses, bb = make_layout(path)
    route_loops(lay, sses, bb)
    drawn = {(p.strand, k): np.asarray(xy, float) for p in lay.na_paths for k, xy in zip(p.nts, p.nt_xy)}
    drawn.update({key: np.asarray(d.point(key), float) for d in lay.stems for key in d.axial})
    tails = {(p.strand, "5′" if p.nts[0] == 0 else "3′"): p for p in lay.na_paths if p.kind == "tail"}
    checked = 0
    for s, e, xy in lay.na_ends:
        if (s, e) not in tails:
            continue
        p = tails[(s, e)]
        own = (s, p.nts[0] if e == "5′" else p.nts[-1])
        nearest = min(drawn, key=lambda k: (float(np.hypot(*(drawn[k] - np.asarray(xy)))), k))
        assert nearest == own, (s, e)
        checked += 1
    assert checked >= 1


@pytest.mark.parametrize("m", [4, 8, 13, 15, 18, 24])
def test_a_hairpin_loop_holds_its_nucleotides_a_bead_apart(m):
    """Random-100 C2b-4 (7Y7C, 8VWQ bead_overlaps): a hairpin loop's circle was sized as if its arc were 0.8 of a
    circle, but on a stem wider than the loop needs the circle is the stem's own width and the arc only about half
    of it, so 13 to 18 nucleotides stood 0.40-0.44 apart and their contact beads (one bead width, na_paths._BEAD)
    overlapped: 25 of 7Y7C's 26 overlaps. Every two neighbours on the arc, and the stem nucleotides at its feet, now
    stand at least a bead apart."""
    from foldmap.na_paths import _BEAD
    from foldmap.na_runs import _arc, _spaced

    stem = PlacedStem(
        "s",
        "P1",
        [0, 0],
        [(0, 0, 0, 40)],
        0.0,
        0.0,
        4.0,
        0.0,
        4.6,
        {(0, 0): 2.0, (0, 40): 2.0},
        {(0, 0): 1, (0, 40): -1},
    )
    pts = _arc(stem, (0, 0), (0, 40), m)
    xy = np.array([pts[0], *_spaced(pts, m), pts[-1]])
    assert np.hypot(*np.diff(xy, axis=0).T).min() >= _BEAD - 1e-6


def test_an_empty_hairpin_loop_on_a_narrow_stem_has_a_circle():
    """Fix round 1 (review M1): with no nucleotide on the loop the step between its feet is the chord whatever the
    radius, so on a stem narrower than a bead (dna_scale under about 0.36) the circle's growth never ended."""
    stem = PlacedStem(
        "s", "P1", [0, 0], [(0, 0, 0, 1)], 0.0, 0.0, 4.0, 0.0, 0.8, {(0, 0): 2.0, (0, 1): 2.0}, {(0, 0): 1, (0, 1): -1}
    )
    _, r = stem.hairpin_circle((0, 0), (0, 1), 0)
    half = float(np.hypot(*(stem.point((0, 1)) - stem.point((0, 0))))) / 2
    assert half < 0.4545 / 2 and r == pytest.approx(half * 1.02)


def test_a_nucleic_acid_without_protein_keeps_the_order_of_its_stems():
    """Random-100 C2b-4 (3OWW, A17 spearman 0.484): the domain-swapped dimer's ten stems were projected into a box
    of 0.98 times their footprints' area; _attract cleared the overlaps by pushing each stem its own way, and
    the pairwise distances lost their order (0.792 as projected, 0.493 after). The projection is now spread by one
    factor to the room the stems take (_room_for_stems) before the pull, which keeps every distance ratio."""
    from foldmap.nacheck import THRESHOLDS, shape
    from helpers import r100_entry

    lay, _, _ = make_layout(r100_entry("3OWW"))
    assert not lay.placed and len(lay.stems) == 10
    assert shape(lay)["spearman"] >= THRESHOLDS["spearman"]


def test_spreading_stems_to_their_room_keeps_every_distance_ratio():
    """_room_for_stems moves stem items by one factor about their middle, homes and places alike, and only when
    every item is a stem and their box lacks the room."""
    from foldmap.layout import _Item, _room_for_stems

    def stem_item(x, y):
        d = PlacedStem("s", "P", [0, 0], [], 0.0, 0.0, 4.0, 0.0, 4.0, {}, {})
        return _Item([d], [(0.0, 0.0)], x, y, (x, y), 0, "A")

    items = [stem_item(0.0, 0.0), stem_item(1.0, 0.0), stem_item(0.0, 2.0)]
    before = np.array([it.home for it in items])
    _room_for_stems(items)
    after = np.array([it.home for it in items])
    d0 = np.hypot(*(before[:, None] - before[None]).transpose(2, 0, 1))
    d1 = np.hypot(*(after[:, None] - after[None]).transpose(2, 0, 1))
    ratio = d1[d0 > 0] / d0[d0 > 0]
    assert ratio.min() > 1.5 and np.ptp(ratio) < 1e-3
    assert [(it.cx, it.cy) for it in items] == [tuple(h) for h in after.tolist()]
    roomy = [stem_item(0.0, 0.0), stem_item(40.0, 0.0)]
    _room_for_stems(roomy)
    assert [it.home for it in roomy] == [(0.0, 0.0), (40.0, 0.0)]


def test_a_tie_label_finds_a_clear_spot_with_a_longer_tick():
    """Random-100 C2b-4 (A6 on 7Y7C and 8VWQ: tie labels on strands and protein loops): a tertiary pair's mark took
    the first turn whose label crossed nothing, and when no turn's did, the first, on a strand (4P9R's tie 3-284 on
    run 23). The tick is now drawn longer (2 and 3 times) before giving up, and protein loops count as strands."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec
    from helpers import r100_entry

    d = draw_spec(figure_spec(r100_entry("4P9R")))
    assert [t for t in d.layout.na_ties if t.label]
    assert [h for h in figcheck.text_collisions(d) if h.startswith("na:tie:")] == []


def _two_stems(gap: float, parent: str = "", rise: float = 0.0):
    """Stems P1 (strand 0's nucleotide 0, leaving rightwards) and P2 (nucleotide 3 or, for joined pieces, 1, entered
    from the left), their facing ends `gap` apart along x and P2 `rise` higher, and a router on a page holding only
    them."""
    from foldmap.layout import Layout
    from foldmap.route import _element_obstacles, _Grid, _Router

    j = 1 if parent else 3
    da = PlacedStem("P1", "P1", [0, 0], [], 0.0, 0.0, 3.0, 0.0, 4.6, {(0, 0): 1.0}, {(0, 0): 1}, parent=parent)
    db = PlacedStem("P2", "P2", [0, 0], [], 3.0 + gap, rise, 3.0, 0.0, 4.6, {(0, j): -1.0}, {(0, j): 1}, parent=parent)
    lay = Layout({}, (-3.0, -4.0 + min(rise, 0.0), 6.0 + gap, 4.0 + max(rise, 0.0)), stems=[da, db])
    return da, db, j, _Router(_Grid(lay), _element_obstacles(lay))


def test_a_strand_whose_face_looks_into_a_margin_still_routes():
    """Random-100 C2b-4 (7Y7C a:843-844, 8VWQ 1:836): a strand's first free cell out of its stem's end face lay in
    another stem's clearance margin, so no route started and it was drawn as a plain curve. Routed again with its
    ends loose (as a cropped gap is), the drawn route's own check decides."""
    from foldmap.na_paths import _face
    from foldmap.na_runs import _strand_between
    from foldmap.route import _PUSH

    da, db, j, router = _two_stems(6.0)
    fa, ea = _face(da, (0, 0), True)
    for t in _PUSH[:3]:  # the cells in front of P1's face, taken by a margin (no body there)
        router.grid.blocked[router.grid.cell(fa + ea * t)] = True
    path = _strand_between(router, 0, 0, j, da, db, None)
    assert path.kind == "run" and not path.fallback


def test_joined_pieces_kept_end_to_end_may_join_straight_farther_apart():
    """Random-100 C2b-4: two pieces of one bent duplex 4.4 apart whose joint no route reaches were drawn as a plain
    curve: a straight joint was allowed only under 2.5. Up to _JOINT_REACH the duplex now carries on straight when the
    line leaves the face along its axis, as here, and is clear (the refusals: the next two tests)."""
    from foldmap.na_paths import _face
    from foldmap.na_runs import _straight_joint

    da, db, j, router = _two_stems(4.4, parent="P")
    fa, ea = _face(da, (0, 0), True)
    fb, _ = _face(db, (0, j), False)
    assert math.dist(fa, fb) == pytest.approx(4.4)
    assert _straight_joint(router, da, db, fa, ea, fb)


def test_a_long_straight_joint_keeps_off_names_and_protein_loops():
    """Random-100 C2b-4 fix round 1 (review I1): beyond _JOINT_NEAR a straight joint is drawn only when it keeps off
    the stems' and domains' names and the protein loops, as a routed strand does; a name box or a loop across the
    line refuses it."""
    from foldmap.na_paths import _face, _Strands
    from foldmap.na_runs import _straight_joint

    da, db, j, router = _two_stems(4.4, parent="P")
    fa, ea = _face(da, (0, 0), True)
    fb, _ = _face(db, (0, j), False)
    mid = (fa + fb) / 2
    name = (mid[0] - 0.3, mid[1] - 0.2, mid[0] + 0.3, mid[1] + 0.2)
    assert not _straight_joint(router, da, db, fa, ea, fb, names=[name])
    loops = _Strands(None)
    loops.add([(mid[0], mid[1] - 1.0), (mid[0], mid[1] + 1.0)])
    assert not _straight_joint(router, da, db, fa, ea, fb, loops=loops)


def test_a_long_straight_joint_leaves_its_face_along_its_axis():
    """Random-100 C2b-4 fix round 1 (review I1; 8VWQ P26a -> P26b, 4.41 apart and 76° off its face, had been drawn
    straight): beyond _JOINT_NEAR the line must leave its face within 60° of its axis (_JOINT_FORWARD); a partner 4.4
    away and 70° off the axis is refused."""
    from foldmap.na_paths import _face
    from foldmap.na_runs import _straight_joint

    da, db, j, router = _two_stems(1.5, parent="P", rise=4.2)
    fa, ea = _face(da, (0, 0), True)
    fb, _ = _face(db, (0, j), False)
    gap = math.dist(fa, fb)
    assert 4.0 < gap < 5.0 and float(np.dot(fb - fa, ea)) / gap < 0.5
    assert not _straight_joint(router, da, db, fa, ea, fb)


def test_strands_keep_their_drawings_off_domain_names(monkeypatch):
    """Random-100 C2b-4 (7Y7C domains view, A16: "D1 · B" under run A:1256-1257): a run's bow is not routed on the grid,
    so of the text it keeps off only what route_nucleic hands it as names, and that was the stems' names alone. The
    domains' names (their boxes as routes keep out of them, label_boxes) are now among them."""
    import foldmap.route_na as route_na
    from foldmap.layout import label_boxes

    seen = []
    real = route_na._strand_between

    def spy(router, s, i, j, da, db, kind, names=(), *rest, **kw):
        seen.append(list(names))
        return real(router, s, i, j, da, db, kind, names, *rest, **kw)

    monkeypatch.setattr(route_na, "_strand_between", spy)
    lay, sses, bb = make_layout(_path("4OO8"), domains="auto")
    route_loops(lay, sses, bb)
    boxes = [box for name, box in label_boxes(lay, []) if name.startswith("domain-label:")]
    assert boxes and seen
    assert all(b in names for names in seen if names for b in boxes)
    # and as drawn (review M5): no run or joint passes through a domain's name between its faces
    from foldmap.na_runs import _off_names

    drawn = [p for p in lay.na_paths if p.kind in ("run", "joint")]
    assert drawn and all(_off_names(p.points, boxes) for p in drawn)


def test_a_mark_held_by_nucleotides_tethers_off_the_stems_ladders():
    """Random-100 C2b-4 fix round 1 (review I2): a ligand held by nucleotides alone took the first clear spot whatever
    its tethers crossed, since its tethers were held off the stems' own boxes, which every tether to a rail meets; once
    the stems were spread (_room_for_stems), 3OWW's GLY A89 sat below P3b and its tethers crossed P3b's rungs. Its
    spot is now chosen, as a protein's ligand's is, for the tethers as drawn, inside the stems' bodies the least
    (_over_stems): each tether enters a stem by no more than its rail's inset."""
    from foldmap.layout import RAIL_IN
    from foldmap.na_paths import _local
    from foldmap.render import ligand_marks
    from foldmap.route import clip_len
    from foldmap.spec import draw_spec, figure_spec
    from helpers import r100_entry

    spec = figure_spec(r100_entry("3OWW"))
    d = draw_spec(spec)
    marks = [m for m in ligand_marks(d.layout, d.loops, d.sses, spec.look) if not m["ligand"].contacts]
    assert len(marks) == 2
    for m in marks:
        for t in m["targets"]:
            for s in d.layout.stems:
                body = (-s.length / 2, -s.width / 2, s.length / 2, s.width / 2)
                inside = clip_len(_local(m["centre"], s), _local(t, s), body)
                assert inside <= RAIL_IN + 0.1, (m["ligand"].name, s.name, inside)


def test_a_mark_held_by_nucleotides_keeps_its_tethers_off_strands():
    """Random-100 C2b-4 fix round 2 (ruling): a mark held by nucleotides judges every clear spot by the stems its
    tethers cross, then the strands (hairpin loops, runs), then their length; scored on the stems alone, it took the
    first clear spot here, its tether across the strand. A protein's mark (no `na`) keeps its first clear spot."""
    from foldmap.render import _mark_spot, _strands_crossed

    strand = np.array([[-0.3, -1.0, 0.1, -1.0]])  # a strand's segment between the anchor and the target
    target = [(0.0, -2.0)]
    assert _strands_crossed((0.0, 0.0), target, [strand]) == 1
    spot = _mark_spot(np.zeros(2), (0.1, 0.1), target, [], [], [], na=([], [strand]))
    assert _strands_crossed(spot, target, [strand]) == 0
    assert np.hypot(*np.subtract(spot, target[0])) < 0.2  # of those, the one with the shortest tether
    assert list(_mark_spot(np.zeros(2), (0.1, 0.1), target, [], [], [])) == [0.0, 0.0]


def test_a_mark_held_by_nucleotides_crosses_fewer_strands_before_it_keeps_off_bodies():
    """Random-100 C2b-4 fix round 3 (ruling): for a mark held by nucleotides the strands its tethers cross are judged
    before their length over element bodies, a continuous measure that almost never ties. Here every spot whose tether
    misses the strand runs a little over a body, and every spot clear of the body crosses the strand: the first kind
    wins."""
    from foldmap.render import _mark_spot, _over_bodies, _strands_crossed

    strand = np.array([[-5.0, -1.0, 0.1, -1.0]])
    body = (0.15, -0.8, 3.0, -0.6)
    target = [(0.0, -2.0)]
    spot = _mark_spot(np.zeros(2), (0.05, 0.05), target, [], [], [], [body], radius=1.0, na=([], [strand]))
    assert _strands_crossed(spot, target, [strand]) == 0
    assert _over_bodies(spot, target, [body]) > 0


def test_a_cropped_gap_from_a_traces_end_is_routed():
    """Random-100 C2b-6 (7Y7C's gaps after Z:37, Z:61 and V:21): a gap with an end on a strand drawn along its
    projection (a trace, which no stem holds) had no owner there, so it was never routed: a curve or a bow when one
    was clear, else drawn straight across whatever lay between (six stems and a terminus). Its trace's end is now a
    point a route leaves from, so it goes round."""
    from types import SimpleNamespace

    from foldmap.layout import Layout
    from foldmap.na_ends import _gap
    from foldmap.na_paths import _face, crosses_body
    from foldmap.route import _element_obstacles, _Grid, _Router

    wall = PlacedStem("P1", "P1", [0, 0], [], 0.0, 0.0, 3.0, 0.0, 20.0, {(0, 0): 1.0}, {(0, 0): 1})
    db = PlacedStem("P2", "P2", [0, 0], [], 6.0, 0.0, 3.0, 0.0, 4.6, {(0, 3): -1.0}, {(0, 3): 1})
    lay = Layout({}, (-9.0, -13.0, 10.0, 13.0), stems=[wall, db])
    router = _Router(_Grid(lay), _element_obstacles(lay))
    fb, eb = _face(db, (0, 3), False)
    pa, ea = np.array([-5.0, 0.0]), np.array([-1.0, 0.0])  # a trace's end left of the wall, pointing away
    path = _gap(router, SimpleNamespace(before=0, after=1, count=5), (pa, ea, None), (fb, eb, db))
    assert not path.fallback
    assert not any(crosses_body(p, q, wall) for p, q in zip(path.points, path.points[1:]))


def test_a_cropped_gap_whose_end_faces_a_stem_starts_from_the_nearest_free_cell():
    """Random-100 C2b-6 (7Y7C's domains view, gap after a:1439): the way straight out of a gap's end ran through a
    neighbouring stem's margin for the whole of _PUSH, so no route started, even with its ends loose, and with no curve
    or bow clear the gap was drawn straight across three stems. A gap may leave its end in any direction ahead of it:
    it now starts from the nearest free cell ahead whose straight stub is clear, and routes on from there."""
    from types import SimpleNamespace

    from foldmap.layout import Layout
    from foldmap.na_ends import _gap
    from foldmap.na_paths import _face, crosses_body
    from foldmap.route import _element_obstacles, _Grid, _Router

    da = PlacedStem("P1", "P1", [0, 0], [], -6.0, 0.0, 3.0, 0.0, 4.6, {(0, 0): 1.0}, {(0, 0): 1})
    wall = PlacedStem("W", "W", [1, 1], [], 0.0, 0.0, 3.0, 0.0, 30.0, {(1, 0): 1.0}, {(1, 0): 1})
    db = PlacedStem("P2", "P2", [0, 0], [], 9.0, 0.0, 3.0, 0.0, 4.6, {(0, 3): -1.0}, {(0, 3): 1})
    fb, eb = _face(db, (0, 3), False)
    front = PlacedStem("N", "N", [2, 2], [], fb[0] - 1.8, fb[1], 3.0, 0.0, 1.0, {(2, 0): 1.0}, {(2, 0): 1})
    lay = Layout({}, (-9.0, -18.0, 12.0, 18.0), stems=[da, wall, db, front])
    router = _Router(_Grid(lay), _element_obstacles(lay))
    assert router.grid.free_cell_along(fb, eb, db, True) is None  # straight out: the neighbour's margin all the way
    fa, ea = _face(da, (0, 0), True)
    path = _gap(router, SimpleNamespace(before=0, after=1, count=5), (fa, ea, da), (fb, eb, db))
    assert not path.fallback
    assert not any(crosses_body(p, q, e) for e in (wall, front) for p, q in zip(path.points, path.points[1:]))


def test_a_joint_no_route_or_line_reaches_is_drawn_as_a_clear_curve():
    """Random-100 C2b-6 (8VWQ's P26a <-> P26b 73-76° off their faces, P124a's pieces turned back from each other, A3
    in six views; 7Y7C's domains view): two pieces of one bent duplex that no route reaches and no straight line or
    elbow may join were drawn as a plain fallback curve, whatever it crossed. The joint is now drawn as a smooth curve
    leaving and entering along the faces' axes when one is clear of every element and stem (its own two by their
    bodies), of the strands but its partner, of names and protein loops; only when none is does it fall back."""
    from foldmap.na_paths import _face, in_body
    from foldmap.na_runs import _strand_between
    from foldmap.route import _PUSH

    da, db, j, router = _two_stems(1.5, parent="P", rise=4.2)  # the partner 70° off the face's axis: no straight line
    fa, ea = _face(da, (0, 0), True)
    fb, eb = _face(db, (0, j), False)
    for port, out in ((fa, ea), (fb, eb)):  # no route: every cell straight out of both faces taken
        for t in _PUSH:
            router.grid.blocked[router.grid.cell(port + out * t)] = True
    path = _strand_between(router, 0, 0, j, da, db, None)
    assert path.kind == "joint" and not path.fallback
    mid = path.points[2:-2]
    assert mid and not any(in_body(p, d) for p in mid for d in (da, db))


def test_a_strand_whose_face_looks_into_a_margin_all_along_routes_from_beside_it():
    """Random-100 C2b-6 (8VWQ's 2:834-835 between P124a's pieces turned back from each other; 7Y7C's domains view): the
    way straight out of a face ran through a neighbouring margin for the whole of _PUSH, so even loose no route
    started and the strand was a plain fallback curve. It now starts from the nearest free cell ahead of the face,
    reached by a clear straight stub (na_ends._from_beside, as a cropped gap does)."""
    from foldmap.na_paths import _face
    from foldmap.na_runs import _strand_between
    from foldmap.route import _PUSH

    da, db, j, router = _two_stems(6.0)
    fa, ea = _face(da, (0, 0), True)
    for t in _PUSH:  # every cell straight out of P1's face taken by a margin (no body there)
        router.grid.blocked[router.grid.cell(fa + ea * t)] = True
    path = _strand_between(router, 0, 0, j, da, db, None)
    assert path.kind == "run" and not path.fallback
