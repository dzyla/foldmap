"""Nucleic-acid text and lines against the protein (plan 15 Task 3, random-100 held-out A6): base letters, end labels,
tie and gap labels keep off elements, protein loops and other strands as protein labels do; a strand's route keeps
off text where it can; a protein-held ligand's mark keeps off the stems; a trace lies on its protein only where the
protein holds it."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from foldmap.cli import make_layout
from foldmap.layout import NT_ROW, PITCH, Placed, _box_meets_body, glyph
from foldmap.route import route_loops
from helpers import fake_sse

NA = Path(__file__).parent / "data" / "na"


def _routed(pdb):
    lay, sses, bb = make_layout(NA / f"{pdb}.cif.gz")
    loops = route_loops(lay, sses, bb)
    return lay, sses, bb, loops


def _block(cx, cy, w, h):
    """A body (an unnamed element) w wide and h high centred at (cx, cy), unturned."""
    return Placed(fake_sse("H", 9000, 9001, chain="Z"), cx, cy, w, 0.0, h)


def _meets(box, e) -> bool:
    from foldmap.na_labels import _outline

    return any(_box_meets_body(box, part) for part in _outline(e))


# ------------------------------------------------------------------------------------------------- traces


def test_a_trace_nucleotide_on_an_element_that_does_not_hold_it_is_nudged_off():
    """13FZ: RNA B, held by its sheet at nucleotides 5-6, was projected across A:1027-1031 of that sheet at nucleotide
    B:3, which nothing holds there (a fallback; its copy D drew clear). Once the elements are placed, a trace nucleotide
    whose segments cross an element not holding them (_trace_clear's test) is nudged off it, a bead or two, and the
    others stay where they were projected (13FZ's B:3 moved 1.13, about 2.5 beads)."""
    from foldmap.layout import _TRACE_HALF
    from foldmap.na_ends import _trace_clear, off_unheld

    h = fake_sse("H", 5, 5)
    trace = np.array([[float(x), 0.0] for x in range(13)])
    layout = SimpleNamespace(
        nucleic=SimpleNamespace(contacts={5: {(0, 12)}}),  # it holds the 3' end only
        placed={h.id: Placed(h, 3.0, 0.0, 0.6, np.pi / 2, 0.5)},  # across nucleotide 3
        ghosts=[],
        stems=[],
        sheets=[],
        na_traces={0: trace},
    )
    assert not _trace_clear(layout, 0, trace)
    moved = off_unheld(layout, 0, trace)
    assert _trace_clear(layout, 0, moved)
    shift = np.hypot(*(moved - trace).T)
    assert shift[3] <= 2.5 * _TRACE_HALF + 1e-9 and np.delete(shift, 3).max() == 0.0
    layout.nucleic.contacts = {5: {(0, 3)}}  # where it holds it, it stays under it
    assert (off_unheld(layout, 0, trace) == trace).all()


# ------------------------------------------------------------------------------------------------- letters


def test_outline_of_a_strand_is_its_shaft_and_head_not_its_pitch_wide_box():
    """2FD8: three nucleotides lying across a jelly roll had no spot clear of the strands' boxes (each the pitch wide,
    touching its neighbours); between two shafts a letter's ink is clear of both."""
    from foldmap.na_labels import _outline

    a = Placed(fake_sse("E", 1, 6), 0.0, 0.0, 4.0, np.pi / 2, PITCH)
    b = Placed(fake_sse("E", 9, 14), PITCH, 0.0, 4.0, -np.pi / 2, PITCH)
    between = (PITCH / 2 - 0.1, -0.2, PITCH / 2 + 0.1, 0.2)
    assert _meets(between, _block(0.0, 0.0, PITCH, 4.0))  # the boxes meet it
    assert not any(_box_meets_body(between, part) for e in (a, b) for part in _outline(e))  # the drawn arrows do not
    assert any(_box_meets_body((-0.1, -0.2, 0.1, 0.2), part) for part in _outline(a))  # on the shaft: met
    head = (-0.45, 1.6, -0.35, 1.7)  # beside the shaft, by the arrowhead's base: met
    assert any(_box_meets_body(head, part) for part in _outline(a))
    assert not any(_box_meets_body((-0.45, -1.0, -0.35, -0.9), part) for part in _outline(a))  # beside the shaft


def test_a_letter_keeps_off_a_protein_loop():
    """13FZ: B:2's and B:8's letters on protein loops; letters were kept off the strands drawn, not off the loops."""
    from foldmap.layout import _LETTER_EM, _ink_points, _reach
    from foldmap.na_labels import _em, _place_letters
    from foldmap.route import clip_len

    lay, _, _, _ = _routed("1EHZ")
    first = _place_letters(lay, lay.na_paths)
    key = sorted(first)[len(first) // 2]
    x, y = first[key]
    loop = [(x - 1.0, y + 0.01), (x + 1.0, y - 0.01)]  # a protein loop through that letter
    letters = _place_letters(lay, lay.na_paths, [loop])
    lx, ly = letters[key]
    assert (lx, ly) != (x, y)
    base = lay.nucleic.strands[key[0]].bases[key[1]]
    box = _reach(_ink_points(base, False, _LETTER_EM * _em(lay), 0.0)) + np.array([lx, ly, lx, ly])
    assert clip_len(*loop, tuple(box)) == 0


def test_only_the_letters_a_protein_loop_still_runs_through_move():
    """After route._off_na_text routed what it could off the letters, a letter a loop still runs through moves off
    it; every other letter stays (a letter kept off the loops from the start left a loop crossing a strand it had
    been routed off before: one more protein-loop crossing in four NA campaign entries)."""
    from foldmap.layout import _LETTER_EM, _ink_points, _reach
    from foldmap.na_labels import _em, _letters_off_loops
    from foldmap.route import clip_len

    lay, _, _, _ = _routed("1EHZ")
    first = dict(lay.na_letters)
    key = sorted(first)[len(first) // 2]
    x, y = first[key]
    loop = [(x - 0.12, y + 0.01), (x + 0.12, y - 0.01)]  # through it, short of its neighbours
    _letters_off_loops(lay, [loop])
    moved = {k for k, xy in lay.na_letters.items() if xy != first[k]}
    assert moved == {key}
    lx, ly = lay.na_letters[key]
    base = lay.nucleic.strands[key[0]].bases[key[1]]
    box = _reach(_ink_points(base, False, _LETTER_EM * _em(lay), 0.0)) + np.array([lx, ly, lx, ly])
    assert clip_len(*loop, tuple(box)) == 0


def test_a_letter_with_every_near_spot_covered_looks_round_its_nucleotide_further_out():
    """2FD8: every spot beside a trace's nucleotides lay on the sheet under it; the letter took the first, on a
    strand. With nothing clear beside it, a letter looks on rings further out (_LETTER_RINGS) for a clear spot."""
    from foldmap.na_labels import _LETTER_GAP, _place_letters
    from foldmap.na_paths import _near

    lay, _, _, _ = _routed("1EHZ")
    tail = next(p for p in lay.na_paths if p.kind == "tail")
    key = (tail.strand, tail.nts[1])
    x, y = tail.nt_xy[1]
    reach = 3.1 * NT_ROW + 0.5  # past the furthest spot beside it, either side and along
    block = _block(x, y, 2 * reach, 2 * reach)
    lay.placed = {**lay.placed, "block": block}
    letters = _place_letters(lay, lay.na_paths)
    spot = letters[key]
    others = [q for k, q in letters.items() if k != key]
    assert not _near(spot, block)
    assert min(np.hypot(*(np.subtract(others, spot)).T)) >= _LETTER_GAP


# ------------------------------------------------------------------------------------------------- tie and gap labels


def test_a_tie_label_whose_box_meets_a_body_is_not_clear():
    """6S13: a tie's label whose centre stood just off a stem had half its box on it (clear judged by the centre)."""
    from foldmap.na_labels import _ties

    lay, _, _, loops = _routed("1EHZ")
    pts = [lp.points for lp in loops]
    tie = next(t for t in lay.na_ties if t.label)
    (_, b, at) = tie.marks[0]
    gw, gh = glyph(lay)
    hw, hh = gw / 2 + 0.05, gh / 2
    u = np.subtract(at, b) / np.hypot(*np.subtract(at, b))
    n = np.array([-u[1], u[0]])
    side = 1.0
    c = np.add(at, n * (side + 0.5 * min(hw, hh)))  # its near edge half a label's half size from the label's centre
    block = _block(float(c[0]), float(c[1]), 2 * side, 2 * side)
    assert not (block.rect[0] < at[0] < block.rect[2] and block.rect[1] < at[1] < block.rect[3])  # the centre is off
    lay.placed = {**lay.placed, "block": block}
    for t in _ties(lay, lay.na_paths, pts):
        for _, _, (x, y) in t.marks:
            box = (x - hw, y - hh, x + hw, y + hh)
            assert not _meets(box, block)


def test_a_gap_label_with_no_clear_spot_beside_its_gap_looks_round_it(monkeypatch):
    """8ZTU: "... 5 nt ..." by a gap 0.54 long among runs: every spot beside it met a run, so it took the least bad,
    on one. With none clear beside it, the label looks round the gap's midpoint (_GAP_RINGS) for a clear spot."""
    from foldmap import na_labels
    from foldmap.na_labels import _along, _em, _gap_box, _gap_half, _label_gaps

    lay, _, _, loops = _routed("4YBB")
    p = next(q for q in lay.na_paths if q.kind == "gap" and q.label)
    monkeypatch.setattr(na_labels, "_GAP_OUT", (1.0,))
    monkeypatch.setattr(na_labels, "_GAP_ALONG", (0.5,))
    hw, hh = _gap_half(p.label, _em(lay))
    q, u = _along(p.points, 0.5)
    normal = np.array([-u[1], u[0]])
    reach = abs(normal[0]) * hw + abs(normal[1]) * hh + 0.1
    planted = [tuple(q + s * reach * normal) for s in (1.0, -1.0)]  # a letter on either spot beside it
    for k, xy in enumerate(planted):
        lay.na_letters[(-1, -1 - k)] = xy
    _label_gaps(lay, lay.na_paths, [lp.points for lp in loops])
    x0, y0, x1, y1 = _gap_box(p, _em(lay))
    gw, gh = glyph(lay)
    g = gw / 2 + 0.05, gh / 2  # a letter's box, as _label_gaps keeps clear of it
    assert not any(x - g[0] < x1 and x0 < x + g[0] and y - g[1] < y1 and y0 < y + g[1] for x, y in planted)


# ------------------------------------------------------------------------------------------------- end labels


def _end_spots(pdb, loops=()):
    from foldmap.route_na import route_nucleic

    lay, _, _ = make_layout(NA / f"{pdb}.cif.gz")
    route_nucleic(lay, None, loops)
    return lay, {(s, e): xy for s, e, xy in lay.na_ends}


@pytest.mark.parametrize("end", ["5′", "3′"])
def test_an_end_label_a_protein_loop_runs_through_moves_off_it(end):
    """7NZ2's domains view: K2:4-80's 5' label (beyond its stem's rail: _face_label found no clear spot) on protein
    loop A4:1270-1279>A4:1281-1292. An end label, a free end's or a stem end's, moves off a protein loop through it,
    as off another strand's line. 1EHZ: its 5' end is paired (a stem end), its 3' end a free end."""
    from foldmap.layout import line_band, lines_meeting, segments_of

    _, first = _end_spots("1EHZ")
    (s, _), (x, y) = next(((k, xy) for k, xy in first.items() if k[1] == end))
    loop = [(x - 1.5, y), (x + 1.5, y)]
    lay, moved = _end_spots("1EHZ", [loop])
    mx, my = moved[(s, end)]
    assert (mx, my) != (x, y)
    g = glyph(lay)[1] / 2 + 0.05
    assert not lines_meeting((mx - g, my - g, mx + g, my + g), segments_of([loop]), line_band(lay)).any()


def test_a_trace_end_label_on_an_element_moves_off_it():
    """2FD8: the 5' label of B:501-503, a trace lying across a sheet, stood on a strand (a trace's labels were never
    judged). It moves round its end nucleotide, further out if need be, to a spot off every element."""
    from foldmap.na_paths import _near

    lay, first = _end_spots("1QBJ")
    (key, (x, y)) = next((k, xy) for k, xy in sorted(first.items()) if k[0] in lay.na_traces)  # the trace's 5' end
    block = _block(x, y, 0.6, 0.6)
    lay2, _, _ = make_layout(NA / "1QBJ.cif.gz")
    lay2.placed = {**lay2.placed, "block": block}
    from foldmap.route_na import route_nucleic

    route_nucleic(lay2, None, ())
    spot = next(xy for s, e, xy in lay2.na_ends if (s, e) == key)
    assert not _near(spot, block, glyph(lay2)[1] / 2)


# ------------------------------------------------------------------------------------------------- routes and marks


def test_a_strand_route_through_text_pays_for_each_text_cell(monkeypatch):
    """8ZTU's domains view: a run's route found only with the keep-outs open ran straight up through α3's name, the
    shortest way. A nucleic-acid route searched with the keep-outs open pays text_toll for each text cell it enters;
    the grid's toll is as it was afterwards."""
    from foldmap import route as route_mod
    from foldmap.layout import Layout
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay = Layout({}, (0.0, 0.0, 10.0, 10.0))
    grid = _Grid(lay, [(4.0, 4.0, 6.0, 6.0)])
    router = _Router(grid, _element_obstacles(lay), [(4.0, 4.0, 6.0, 6.0)])
    seen = []

    def drawn(self, *a, **k):
        seen.append(None if self.grid.toll is None else self.grid.toll.copy())
        return None if len(seen) == 1 else [(1.0, 5.0), (9.0, 5.0)]

    monkeypatch.setattr(route_mod._Router, "drawn", drawn)
    pts, over = router.route((1.0, 5.0), (1, 0), None, (9.0, 5.0), (-1, 0), None, text_toll=25.0)
    assert seen[0] is None and seen[1] is not None
    assert (seen[1][grid.soft] == 25.0).all() and (seen[1][~grid.soft] == 0.0).all()
    assert grid.toll is None and over
    seen.clear()
    router.route((1.0, 5.0), (1, 0), None, (9.0, 5.0), (-1, 0), None)  # a protein loop: no toll
    assert seen == [None, None]


def test_nucleic_acid_routes_ask_for_the_text_toll(monkeypatch):
    from foldmap import na_runs
    from foldmap.na_paths import TEXT_TOLL

    asked = []
    router = SimpleNamespace(
        route=lambda *a, **k: (asked.append(k.get("text_toll")), ([(0, 0), (1, 1)], True))[1], grid=None
    )
    na_runs._faces_route(router, (0, 0), (1, 0), None, (1, 1), (-1, 0), None)
    assert asked == [TEXT_TOLL]


def test_a_rounded_corner_is_one_line_struck():
    """6L9Z's domains view: N of chain A, grazed by a strand's rounded corner (two short segments within the
    clearance), moved onto another strand straight through it (one segment): _clear_termini counted segments."""
    from foldmap.layout import segments_of
    from foldmap.route import _lines_struck

    corner = [(1.05 + 0.3 * np.cos(t), -0.05 + 0.3 * np.sin(t)) for t in np.linspace(np.pi / 2, np.pi, 9)]
    straight = [(0.5, -1.0), (0.5, 1.0)]
    lines = [corner, straight]
    segs = segments_of(lines)
    owner = np.repeat(np.arange(len(lines)), [len(p) - 1 for p in lines])
    box = (0.0, -0.5, 0.8, 0.3)
    assert _lines_struck(box, segs, owner, 0.07) == 2
    assert _lines_struck((0.2, -0.5, 0.79, -0.2), segs, owner, 0.07) == 1  # the straight one only


def test_a_protein_held_ligand_keeps_off_the_stems():
    """8ZTU: Mg a6201, held by a protein residue, was drawn on stem 188: a mark held by protein kept off the elements,
    labels and loops, not the nucleic acid (only a mark held by nucleotides alone did). 1IHF's Cd, held by the protein,
    with a stem moved onto its mark: the mark moves off it."""
    import dataclasses

    from foldmap.render import ligand_marks
    from foldmap.style import Style

    lay, sses, _, loops = _routed("1IHF")
    look = Style(ligands="all")
    marks = [m for m in ligand_marks(lay, loops, sses, look) if m["ligand"].contacts]
    assert marks
    (cx, cy) = marks[0]["centre"]
    d = lay.stems[0]
    lay.stems[0] = dataclasses.replace(d, cx=cx, cy=cy)
    moved = [m for m in ligand_marks(lay, loops, sses, look) if m["ligand"] is marks[0]["ligand"]][0]
    (mx, my), (hx, hy) = moved["centre"], moved["half"]
    assert not _box_meets_body((mx - hx, my - hy, mx + hx, my + hy), lay.stems[0])


def test_a_label_never_stands_across_a_body_from_what_it_names():
    """Review of held-out 2FD8: rings, outlines and far reaches put B:502's letter left of strand L naming a nucleotide
    right of it, and the 5' label past strand B. A spot found so must see what it names: the line between them crosses
    no body's drawn outline but the one the named point lies under (_in_sight)."""
    from foldmap.na_labels import _boxes_of, _in_sight

    wall = Placed(fake_sse("E", 1, 6), 0.0, 0.0, 4.0, np.pi / 2, PITCH)  # a strand standing up at x = 0
    bodies = [wall]
    rects = _boxes_of(bodies)
    assert not _in_sight((-1.0, 0.0), (1.0, 0.0), bodies, rects)  # across its shaft
    assert _in_sight((1.0, 1.0), (1.0, -1.0), bodies, rects)  # beside it
    assert _in_sight((-1.0, 0.0), (0.1, 0.0), bodies, rects)  # to a point under it: its own body
    assert _in_sight((0.6, 0.0), (0.6, 3.0), bodies, rects)  # beside the shaft and its head, past the tip


def test_a_crowded_letter_with_no_spot_in_sight_keeps_the_nearest_least_bad_one():
    """With every near spot covered and every spot further out across a frame of bodies from its nucleotide, a
    letter keeps the least bad spot beside it (the collision stays, honestly counted) rather than a far one that would
    name what lies beyond the frame."""
    from foldmap.na_labels import _boxes_of, _in_sight, _place_letters

    lay, _, _, _ = _routed("1EHZ")
    tail = next(p for p in lay.na_paths if p.kind == "tail")
    key = (tail.strand, tail.nts[1])
    x, y = tail.nt_xy[1]
    reach = 3.1 * NT_ROW + 0.5
    frame = [
        Placed(fake_sse("H", 9100 + 10 * k, 9101 + 10 * k, chain="Z"), x + dx, y + dy, 2 * reach + 1.0, a, 0.1)
        for k, (dx, dy, a) in enumerate(
            (
                (reach + 0.2, 0.0, np.pi / 2),
                (-reach - 0.2, 0.0, np.pi / 2),
                (0.0, reach + 0.2, 0.0),
                (0.0, -reach - 0.2, 0.0),
            )
        )
    ]
    block = _block(x, y, 2 * reach, 2 * reach)
    lay.placed = {**lay.placed, "block": block, **{f"frame{k}": f for k, f in enumerate(frame)}}
    spot = _place_letters(lay, lay.na_paths)[key]
    assert np.hypot(spot[0] - x, spot[1] - y) <= 3.1 * NT_ROW + 1e-6  # one of the spots beside it
    bodies = list(lay.placed.values())
    assert _in_sight(spot, (x, y), bodies, _boxes_of(bodies))


# ------------------------------------------------------------------------------------------------- na_line


LMB = Path(__file__).parent / "data" / "1LMB.cif"


def _widths(fig) -> dict[str, list[float]]:
    """Each nucleic-acid backbone part's stroke width (points, as the SVG writes it), by kind."""
    ax = fig.axes[0]
    out: dict[str, list[float]] = {"run": [], "gap": [], "rails": [], "block": [], "legend-na": []}
    for a in (*ax.lines, *ax.patches):
        gid = a.get_gid() or ""
        kind = gid.split(":")[3] if gid.startswith("na:stem:") and len(gid.split(":")) > 3 else gid
        if gid.startswith("na:run:"):
            out["gap" if a.get_linestyle() not in ("-", "solid") else "run"].append(float(a.get_linewidth()))
        elif kind in ("rails", "block", "legend-na"):
            out[kind].append(float(a.get_linewidth()))
    return out


@pytest.mark.parametrize("na_style", ["ladder", "block"])
def test_na_line_widens_runs_rails_and_block_edges_in_proportion(na_style):
    """At na_line 1.5 a run's stroke, a ladder's rails, a block's edge and the legend's backbone stroke are 1.5 times
    as wide as at 1 (a cropped gap's dashes stay); the page is the same size, so the widths compare directly."""
    from foldmap.cli import make_figure_and_loops
    from foldmap.style import resolve_style

    def widths(k):
        look = resolve_style(None, None, [f"na_style={na_style}", f"na_line={k}"])
        fig = make_figure_and_loops(LMB, look=look)[0]
        return _widths(fig), tuple(fig.get_size_inches())

    (one, size1), (wide, size2) = widths(1.0), widths(1.5)
    assert size1 == size2
    assert one["run"] and one[na_style if na_style == "block" else "rails"] and one["legend-na"]
    for kind in ("run", "rails", "block", "legend-na"):
        assert len(wide[kind]) == len(one[kind])
        assert wide[kind] == pytest.approx([1.5 * w for w in one[kind]]), kind
    assert wide["gap"] == one["gap"]


def test_na_line_widens_the_ribbon_tube():
    """The ribbon's backbone tube is na_line times as wide about the same centre line."""
    from foldmap.render import _dna_turns

    lay, _, _ = make_layout(LMB)
    d = lay.stems[0]
    for slot in (0, 1):
        front, back = _dna_turns(d, slot)
        front2, back2 = _dna_turns(d, slot, 1.5)
        for q, q2 in zip(front + back, front2 + back2):
            n = len(q) // 2  # a band: the centre line offset one way, then back the other
            centre = (q[:n] + q[n:][::-1]) / 2
            assert (q2[:n] + q2[n:][::-1]) / 2 == pytest.approx(centre)
            assert q2[:n] - centre == pytest.approx(1.5 * (q[:n] - centre))


def test_text_bands_grow_by_the_strand_stroke_only():
    """line_band and _run_band (what texts keep from a line, as A6 reads it) grow by half the extra stroke of a
    run, NA_RUN_W times na_line, and by nothing else; with a loop wider than the run, the loop's stroke rules."""
    from dataclasses import replace

    from foldmap.layout import NA_RUN_W, _run_band, line_band, line_pt, page_text
    from foldmap.style import resolve_style

    lay, sses, bb = make_layout(LMB)
    route_loops(lay, sses, bb)  # the strands' drawn paths
    assert lay.na_line == 1.0 and lay.na_paths and lay.loop_stroke < NA_RUN_W
    unit = page_text(lay)[1]
    pt = line_pt(1.0 / (72 * unit))
    for k in (1.25, 1.5, 2.0):
        wide = replace(lay, na_line=k)
        extra = pt * NA_RUN_W * (k - 1.0) / 2 * unit
        assert line_band(wide) - line_band(lay) == pytest.approx(extra)
        assert _run_band(wide) - _run_band(lay) == pytest.approx(extra)
    heavy = replace(lay, loop_stroke=3.0)
    assert line_band(replace(heavy, na_line=2.0)) == line_band(heavy)  # 2.8 < 3: the loop rules
    assert make_layout(LMB, look=resolve_style(None, None, ["na_line=2"]))[0].na_line == 2.0  # the layout carries it


def test_text_keeps_off_a_strand_drawn_twice_as_wide():
    """A6 (no text crossed by a line, bodies included) passes on 1LMB with the backbone twice as wide."""
    from foldmap import criteria, figcheck
    from foldmap.spec import draw_spec, figure_spec

    spec = figure_spec(LMB, style=["na_line=2", "na_backbone=chain"])
    drawn = draw_spec(spec)
    got = criteria.check_figure(drawn, spec, LMB)["criteria"]
    if reason := figcheck.ink_font_reason(drawn.figure):  # A6's ink is read in Liberation Sans only: not measured here
        assert got["A6"]["pass"] is None and got["A6"]["here"] == reason, got["A6"]
    else:
        assert got["A6"]["pass"], got["A6"]
