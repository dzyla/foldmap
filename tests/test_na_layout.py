"""Nucleic acids on the page: stems placed like elements, strands routed, four drawing styles, the explorer."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_figure_and_loops, make_layout
from foldmap.layout import PlacedStem
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
    lay, _, _ = make_layout(_path(pdb))
    pieces = len(lay.stems)
    assert pieces >= len(lay.nucleic.structure.stems) > 0
    boxes = [p.rect for p in lay.placed.values()] + [d.rect for d in lay.stems]
    assert not any(_overlap(a, b) for k, a in enumerate(boxes) for b in boxes[k + 1 :])
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


def _joints(lay):
    """(exit port, entry port, piece width) where a strand goes from one piece of a bent duplex to the next."""
    at = {key: d for d in lay.stems for key in d.axial}
    out = []
    for s, strand in enumerate(lay.nucleic.strands):
        for i in range(len(strand) - 1):
            d, e = at.get((s, i)), at.get((s, i + 1))
            if d is not None and e is not None and d is not e and d.parent and d.parent == e.parent:
                out.append((d.port((s, i), True)[0], e.port((s, i + 1), False)[0], d.width))
    return out


@pytest.mark.parametrize("pdb", ["1KX5", "1IHF"])
def test_a_bent_duplex_is_a_chain_of_pieces_end_to_end(pdb):
    """Consecutive pieces stay near each other. Box separation parts tilted pieces, so a joint can be a couple of
    duplex widths long; Task 7 routes joints as loops between stem ends (route_nucleic, _clear). This bounds the
    current lengths (1KX5 at most 2.47 widths, 1IHF 1.60) against regressions."""
    lay, _, _ = make_layout(_path(pdb))
    joints = _joints(lay)
    assert len(joints) >= (40 if pdb == "1KX5" else 4)  # both strands at every joint
    gaps = [float(np.hypot(*(q - p))) / w for p, q, w in joints]
    assert max(gaps) <= 2.6, sorted(gaps)[-5:]


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
        f"lay, _, _ = make_layout({str(_path('30QM'))!r})\n"
        "print([(d.id, round(d.cx, 6), round(d.cy, 6), round(d.angle, 6)) for d in lay.stems])"
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


@pytest.mark.xfail(reason="beads on single-stranded nucleotides arrive with route_nucleic in Task 7", strict=True)
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
