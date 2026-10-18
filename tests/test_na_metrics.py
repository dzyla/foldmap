"""Nucleic-acid metrics on the held-out sample (plan 15 Task 4, random-100 held-out A17): beads of a chain the figure
does not draw, a crop piece too short to draw, beads of runs side by side, gaps across free ends, gap labels read from
either end of their line, counts on tie lines, end labels between strand shafts."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from foldmap.layout import PITCH, SHAFT, Placed
from helpers import fake_sse

DATA = Path(__file__).parent / "data"


# ------------------------------------------------------------------------------------------------- beads


def test_a_contact_of_a_chain_the_figure_does_not_draw_is_no_missing_bead():
    """6S13: chain l (135 residues, no helix or strand) is not drawn and the figure says so (scope: not drawn, chain
    l); its 134 contacts with the RNA had no bead (no colour of its chain is drawn) and were counted as missing beads.
    A protein chain with no element (Layout.coil_chains) is "not drawn", expected, like a chain outside the focus."""
    from foldmap.cli import make_layout
    from foldmap.io import read_model
    from foldmap.nacheck import verify_contacts
    from foldmap.render import draw
    from foldmap.route import route_loops

    path = DATA / "na" / "1IHF.cif.gz"
    model = read_model(path)
    lay, sses, bb = make_layout(path, model=model)
    na = lay.nucleic
    chain = na.contact_chain[sorted(na.contacts)[0]]
    # the figure draws none of that chain's beads, as for a chain with no helix or strand
    na.contacts = {r: v for r, v in na.contacts.items() if na.contact_chain.get(r) != chain}
    fig = draw(lay, route_loops(lay, sses, bb), sses)
    report = verify_contacts(model, fig.foldmap_beads, lay, bb)
    assert {m.pair.residue.split(":")[0] for m in report.missing} == {chain}  # read as missing beads
    lay.coil_chains = [(chain, 99)]  # said not drawn: nothing of it to bead
    report = verify_contacts(model, fig.foldmap_beads, lay, bb)
    assert report.missing == []
    assert {m.reason for m in report.missed if m.pair.residue.startswith(f"{chain}:")} == {"not drawn"}


def test_beads_of_two_runs_side_by_side_are_kept_apart():
    """6S13: two runs routed a grid cell apart for several nucleotides, spaced differently (1.45 and 1.0 apart along
    them), drew three pairs of beads on each other: no slide of a whole run clears them all. Each nucleotide still on
    another's spot slides alone, within the room its neighbours leave, along its own run."""
    from foldmap.na_labels import _SLIDES, _unbunch
    from foldmap.na_paths import _BEAD, NAPath, _spaced

    first = [(0.0, 0.0), (10.0, 0.0)]
    second = [(0.0, 0.25), (10.0, 0.25)]
    a = NAPath("run", 0, list(range(1, 7)), first, _spaced(first, 6))
    b = NAPath("run", 1, list(range(1, 10)), second, _spaced(second, 9))

    def closest(p, q):
        return min(float(np.hypot(*(np.subtract(x, y)))) for x in p.nt_xy for y in q.nt_xy)

    assert closest(a, b) < _BEAD
    whole = [  # every slide of the second run as a whole leaves a pair of beads on each other
        min(float(np.hypot(*np.subtract(x, y))) for x in a.nt_xy for y in _spaced(second, 9, f)) for f in _SLIDES
    ]
    assert max(whole) < _BEAD
    _unbunch(SimpleNamespace(stems=[]), [a, b])
    assert closest(a, b) >= _BEAD - 1e-9
    xs = [x for x, _ in b.nt_xy]
    assert all(y == 0.25 for _, y in b.nt_xy) and all(q - p >= _BEAD - 1e-9 for p, q in zip(xs, xs[1:]))
    assert xs[0] >= _BEAD - 1e-9 and xs[-1] <= 10.0 - _BEAD + 1e-9  # never on its ends


# ------------------------------------------------------------------------------------------------- gaps


def test_a_gap_keeps_off_a_free_end_where_it_can():
    """6S13: the gaps after a:36 and a:462 ran across the free 5' end of a:37-72 they lead to (na_crossings 2): a gap
    paid only the reuse cost to cross a hairpin loop or free end. Like a run's first pass, it now keeps off them when a
    route does."""
    from foldmap.cli import make_layout
    from foldmap.na_ends import _gap
    from foldmap.na_paths import _Strands
    from foldmap.nucleic import NAGap
    from foldmap.route import _element_obstacles, _Grid, _Router

    lay, _, _ = make_layout(DATA / "na" / "1EHZ.cif.gz")
    probe = _Grid(lay)
    free = ~probe.blocked
    # an open stretch of the page, 8 x 4 units clear of everything
    nx, ny = int(8 / probe.res), int(4 / probe.res)
    spot = next(
        (i, j) for i in range(probe.nx - nx) for j in range(probe.ny - ny) if free[i : i + nx, j : j + ny].all()
    )
    x0, y0 = probe.center(spot)
    pa, pb = np.array([x0 + 1.0, y0 + 2.0]), np.array([x0 + 7.0, y0 + 2.0])
    end = [(x0 + 4.0, y0 + 1.0), (x0 + 4.0, y0 + 3.0)]  # a free end standing across the straight way

    def gap(strands):
        router = _Router(_Grid(lay), _element_obstacles(lay))
        held = None
        if strands:
            held = _Strands(router.grid)
            held.add(end, fixed=True)
        path = _gap(router, NAGap(0, 1, 5), (pa, np.array([1.0, 0.0]), None), (pb, np.array([-1.0, 0.0]), None), held)
        check = _Strands(router.grid)
        check.add(end, fixed=True)
        return path, check.crosses(path.points, fixed_only=True)

    plain, crossed = gap(False)
    assert crossed and not plain.fallback  # the straight way, across it
    kept, crossed = gap(True)
    assert not crossed and not kept.fallback


def test_a_gap_label_is_read_from_either_end_of_its_line():
    """8ZTU's domains view: "… 5 nt …" had one clear spot, beside the end of its dashed line but across a stem from
    the line's midpoint, and was put on a stem instead. A spot round the midpoint is taken when it is in sight of the
    midpoint or of either end."""
    from foldmap.na_labels import _boxes_of, _gap_ring

    stem = SimpleNamespace(cx=5.25, cy=2.0, length=4.5, width=2.0, angle=0.0)  # x 3..7.5, y 1..3
    line = [(0.0, 0.0), (8.0, 0.0), (8.0, 5.0)]  # bent round the stem; its midpoint (6.5, 0) below it
    hw, hh = 1.2, 0.3
    bodies = [stem]
    rects = _boxes_of(bodies)

    def only_above(box):  # the clear spots: above the stem
        return box[1] > 3.2

    found = _gap_ring(line, hw, hh, only_above, bodies, rects)
    assert found is not None
    (x, y), _ = found
    assert y - hh > 3.2 and x > 7.5  # beside the line's far end, in sight of it
    # judged from the midpoint alone it would be refused: the stem stands between them
    from foldmap.na_labels import _in_sight
    from foldmap.na_paths import midway

    assert not _in_sight((x, y), midway(line), bodies, rects)


# ------------------------------------------------------------------------------------------------- counts and ends


def test_counts_keep_off_tie_lines():
    """6S13: a loop's "32 aa" was placed on the dotted line of the tie 448-506, which loop and tail counts did not
    keep off (they kept off strands only). na_strand_lines gives them the ties' lines: a short tie's dotted line, a
    long one's two ticks."""
    from foldmap.layout import na_strand_lines
    from foldmap.na_paths import NAPath, NATie

    lay = SimpleNamespace(
        na_paths=[NAPath("run", 0, [], [(0.0, 0.0), (1.0, 0.0)], [])],
        na_ties=[
            NATie((0, 1), (0, 9), [(0.0, 1.0), (3.0, 1.0)]),
            NATie(
                (0, 2),
                (0, 8),
                [(0.0, 2.0), (9.0, 2.0)],
                "a",
                [((0, 2), (0, 2.3), (0, 2.6)), ((9, 2), (9, 2.3), (9, 2.6))],
            ),
        ],
    )
    lines = [np.asarray(q, float).tolist() for q in na_strand_lines(lay)]
    assert lines == [[[0, 0], [1, 0]], [[0, 1], [3, 1]], [[0, 2], [0, 2.3]], [[9, 2], [9, 2.3]]]


def test_an_end_label_may_stand_beside_a_strand_shaft():
    """2FD8: the 5′ label of a nucleotide lying under a sheet found no room within reach: _label_room judged the
    strands by their pitch-wide boxes. It judges them by their drawn outline (shaft and head), as letters are."""
    from foldmap.na_ends import _label_room

    strand = Placed(fake_sse("E", 1, 6), 0.0, 0.0, 4.0, np.pi / 2, PITCH)  # pointing up, its head at the top
    grid = SimpleNamespace(cell=lambda q: (0, 0), nx=1, ny=1, used=np.zeros((1, 1)))
    router = SimpleNamespace(grid=grid, keep_out=[])
    lay = SimpleNamespace(width=20.0, height=20.0, font_scale=1.0, placed={"s": strand}, stems=[])
    from foldmap.layout import glyph

    hh = glyph(lay)[1] / 2
    beside = SHAFT / 2 + hh + 0.05 + 0.01  # just clear of the shaft, inside the pitch-wide box
    assert beside < PITCH / 2 + hh + 0.05
    assert _label_room(router, lay, (beside, -1.0), [])
    assert not _label_room(router, lay, (beside - 0.05, -1.0), [])  # on the shaft's margin
    assert not _label_room(router, lay, (beside, 1.8), [])  # beside the arrowhead, which is wider
