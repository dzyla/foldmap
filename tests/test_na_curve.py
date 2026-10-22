"""Plan 25 Task 1: a bent duplex drawn as one curve (Style.na_bend curved, na_curve) instead of straight pieces
joined by routed joints."""

import math
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_layout
from foldmap.layout import DNA_RISE, NT_ROW, SCALE, glyph
from foldmap.na_curve import _R_MIN, _isotonic, bent_duplexes, build_curves
from foldmap.nacheck import _crossings, _segments
from foldmap.route import route_loops
from foldmap.spec import draw_spec, figure_spec
from foldmap.style import CHOICES, LAYOUT_KEYS, Style, resolve_style

NA = Path(__file__).parent / "data" / "na"


def _path(pdb):
    return NA / f"{pdb}.cif.gz"


def _routed(pdb, bend="curved"):
    lay, sses, bb = make_layout(_path(pdb), look=replace(Style(), na_bend=bend))
    loops = route_loops(lay, sses, bb)
    return lay, loops


def test_na_bend_is_a_layout_key_pieces_by_default():
    assert Style().na_bend == "pieces" and CHOICES["na_bend"] == ("pieces", "curved") and "na_bend" in LAYOUT_KEYS
    assert resolve_style(overrides=["na_bend=curved"]).na_bend == "curved"
    with pytest.raises(ValueError, match="na_bend must be one of pieces, curved"):
        resolve_style(overrides=["na_bend=smooth"])


def test_pieces_draw_no_curve():
    lay, _ = _routed("1IHF", "pieces")
    assert lay.na_bend == "pieces" and lay.na_curves == []
    assert any(p.kind == "joint" for p in lay.na_paths)


def test_isotonic_pools_adjacent_violators():
    assert np.allclose(_isotonic([1.0, 3.0, 2.0, 4.0]), [1.0, 2.5, 2.5, 4.0])
    assert np.allclose(_isotonic([3.0, 1.0], [1.0, 3.0]), [1.5, 1.5])


def test_a_bent_duplex_is_one_curve_through_its_pieces_in_order():
    """1IHF's DNA 1 (three pieces): one curve; its pieces in order along the first arm, each keeping its own rung
    spacing around its place on the curve, never closer to the next than the duplex's own step; the first arm on
    one rail, 5'->3' along the curve, the second back; each rung straight across; no bend tighter than the inner
    rail and its letters allow."""
    lay, _ = _routed("1IHF")
    assert len(lay.na_curves) == 1
    c = lay.na_curves[0]
    (parent, ks), *_ = bent_duplexes(lay.stems)
    assert c.pieces == ks and len(ks) == 3 and c.name == "DNA 1" and not c.gyres and not c.overlap
    held = {key for k in ks for key in lay.stems[k].axial}
    assert set(c.at) == held == set(c.side)
    pieces = [lay.stems[k] for k in ks]
    assert len({c.side[key] for d in pieces for key in d.axial if d.side[key] == d.first}) == 1  # one rail
    firsts = sorted((key for d in pieces for key in d.axial if d.side[key] == d.first), key=lambda k: k[1])
    seconds = sorted((key for d in pieces for key in d.axial if d.side[key] != d.first), key=lambda k: k[1])
    assert np.all(np.diff([c.at[k] for k in firsts]) > 0)  # 5'->3' along the curve
    assert np.all(np.diff([c.at[k] for k in seconds]) < 0)  # the second arm back
    for d in pieces:  # each piece at its own spacing
        keys = sorted(d.axial)
        assert np.allclose([c.at[k] - c.at[keys[0]] for k in keys], [d.axial[k] - d.axial[keys[0]] for k in keys])
    for d, e in zip(pieces, pieces[1:]):  # never closer than the duplex's own step
        assert min(c.at[k] for k in e.axial) - max(c.at[k] for k in d.axial) >= DNA_RISE * SCALE - 1e-6
    for x, i, y, j in c.pairs:
        assert c.at[(x, i)] == c.at[(y, j)] and c.side[(x, i)] == -c.side[(y, j)]
        rung = c.point((y, j)) - c.point((x, i))
        t = c.frames(c.at[(x, i)])[1][0]
        assert abs(rung @ t) < 1e-3 and math.isclose(float(np.hypot(*rung)), 2 * c.rail, abs_tol=1e-3)
    lo, hi = c.extent()
    s = np.arange(lo + 1.0, hi - 1.0, 0.25)  # the radius of the circle through the curve 1 before and after
    a, m, b = c.frames(s - 1.0)[0], c.frames(s)[0], c.frames(s + 1.0)[0]
    ab, am, mb = (np.hypot(*(q - p).T) for p, q in ((a, b), (a, m), (m, b)))
    area2 = np.abs((m[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]) - (m[:, 1] - a[:, 1]) * (b[:, 0] - a[:, 0]))
    assert float((ab * am * mb / np.maximum(2 * area2, 1e-12)).min()) >= 0.9 * _R_MIN * c.rail
    for d in pieces:  # through (near) the pieces' centres
        assert float(np.hypot(c.xy[:, 0] - d.cx, c.xy[:, 1] - d.cy).min()) < 1.5, d.name


def test_the_curve_is_drawn_in_place_of_its_pieces_and_joints():
    """The drawing: one ladder (the first piece's gids), no piece's own rails, no joint, one name; the contact beads
    and the strands' 5'/3' labels with the nucleotides on the curve."""
    drawn = draw_spec(figure_spec(_path("1IHF"), style=["na_bend=curved"]))
    lay, ax = drawn.layout, drawn.figure.axes[0]
    c = lay.na_curves[0]
    gids = {a.get_gid() for a in [*ax.patches, *ax.lines, *ax.texts, *ax.collections] if a.get_gid()}
    k0 = c.pieces[0]
    assert f"na:stem:{k0}:rails" in gids and f"na:label:{k0}" in gids
    assert not any(f"na:stem:{k}:rails" in gids or f"na:label:{k}" in gids for k in c.pieces[1:])
    joints = [n for n, p in enumerate(lay.na_paths) if p.kind == "joint"]
    assert joints and not any(f"na:run:{n}" in gids for n in joints)
    names = [t for t in ax.texts if t.get_gid() == f"na:label:{k0}"]
    assert [t.get_text() for t in names] == ["DNA 1"]
    rungs = [a for a in ax.lines if a.get_gid() == f"na:stem:{k0}:rung"]
    assert len(rungs) == len(c.pairs)
    beads = drawn.figure.foldmap_beads.marks
    on = [key for marks in beads.values() for key, _, _ in marks if key in c.at]
    assert on
    for coll in ax.collections:
        chain = (coll.get_gid() or "").removeprefix("na:contact:")
        if chain in beads:
            keys = [key for key, _, _ in beads[chain]]
            for key, xy in zip(keys, np.asarray(coll.get_offsets(), float)):
                if key in c.at:
                    assert np.allclose(xy, c.point(key))
    na = lay.nucleic
    for s, end, xy in lay.na_ends:
        key = (s, 0) if end == "5′" else (s, len(na.strands[s]) - 1)
        if key in c.at:
            assert float(np.hypot(*(np.asarray(xy) - c.point(key)))) < 3.0, (s, end)


@pytest.mark.parametrize("style", ["ribbon", "block", "letters"])
def test_every_stem_style_draws_the_curve(style):
    drawn = draw_spec(figure_spec(_path("1IHF"), style=["na_bend=curved", f"na_style={style}"]))
    c = drawn.layout.na_curves[0]
    gids = {a.get_gid() for a in [*drawn.figure.axes[0].patches, *drawn.figure.axes[0].lines] if a.get_gid()}
    k0 = c.pieces[0]
    want = {"ribbon": f"na:stem:{k0}:front", "block": f"na:stem:{k0}:block", "letters": f"na:stem:{k0}:rung"}[style]
    assert want in gids
    assert not any(g.startswith(f"na:stem:{k}:") for k in c.pieces[1:] for g in gids)


def test_the_nucleosome_gyres_are_drawn_one_inside_the_other():
    """1KX5: a curve through the pieces in order crosses itself where today's layout alternates the two gyres'
    pieces (the 13-0 mockup: 28 bp overlapping, centrelines crossing); drawn as gyres one inside the other, no two
    stretches 40 or more base pairs apart come within the band (rails and letters) of each other, and no rail
    crosses another."""
    lay, _ = _routed("1KX5")
    (c,) = lay.na_curves
    assert c.gyres and not c.overlap and len(c.pieces) == 22 and len(c.pairs) == 146
    band = 2 * (c.rail + NT_ROW + glyph(lay)[1])
    us = np.array(sorted(c.at[(a, i)] for a, i, _, _ in c.pairs))
    lo, hi = c.span()
    s = np.arange(lo, hi, 0.25)
    p = c.frames(s)[0]
    bp = np.interp(s, us, np.arange(len(us)))
    dist = np.hypot(p[:, None, 0] - p[None, :, 0], p[:, None, 1] - p[None, :, 1])
    far = np.abs(bp[:, None] - bp[None, :]) >= 40
    assert dist[far].min() >= band
    rails = []
    for side in (1, -1):
        at = [c.at[k] for k in c.at if c.side[k] == side]
        rails.append(_segments(c.offset_line(min(at), max(at), side * c.rail)))
    assert _crossings(rails[0], rails[1]) == 0 and all(_crossings(r, r) == 0 for r in rails)


def test_the_curve_is_on_the_grid_and_the_same_every_time():
    """A12: the curve's points and its nucleotides' places are on the layout's grid, and built again they are the
    same."""
    lay, _ = _routed("1IHF")
    c = lay.na_curves[0]
    for v in (c.xy.ravel(), c.s, c.tangent.ravel(), np.array(list(c.at.values()))):
        assert np.array_equal(v, np.round(v * 1e4) / 1e4)
    (again,) = build_curves(lay)
    assert np.array_equal(again.xy, c.xy) and again.at == c.at and again.name_xy == c.name_xy


def test_the_curve_does_not_depend_on_the_hash_seed():
    code = (
        "import hashlib, sys\n"
        "from dataclasses import replace\n"
        "from foldmap.cli import make_layout\n"
        "from foldmap.route import route_loops\n"
        "from foldmap.style import Style\n"
        "lay, sses, bb = make_layout(sys.argv[1], look=replace(Style(), na_bend='curved'))\n"
        "route_loops(lay, sses, bb)\n"
        "c = lay.na_curves[0]\n"
        "h = hashlib.sha256(c.xy.tobytes() + repr(sorted(c.at.items())).encode() + repr(lay.na_ends).encode())\n"
        "print(h.hexdigest())\n"
    )
    out = set()
    for seed in ("0", "1"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        run = subprocess.run([sys.executable, "-c", code, str(_path("1IHF"))], capture_output=True, text=True, env=env)
        assert run.returncode == 0, run.stderr
        out.add(run.stdout.strip())
    assert len(out) == 1


# ---- Task 2: the pieces laid out end to end on the curve, the protein round it ----


def _band_points(c, lay):
    """Points over the curve's band: along both rails and both rows of letters, every 0.1."""
    lo, hi = c.span()
    s = np.arange(lo, hi + 0.05, 0.1)
    return np.vstack([c.to_page(s, off) for off in (c.rail, -c.rail, c.rail + NT_ROW, -(c.rail + NT_ROW))])


def _under(pts, e) -> bool:
    cs, sn = math.cos(e.angle), math.sin(e.angle)
    x, y = pts[:, 0] - e.cx, pts[:, 1] - e.cy
    u, v = x * cs + y * sn, y * cs - x * sn
    return bool(((np.abs(u) < e.length / 2) & (np.abs(v) < e.width / 2)).any())


def _rails(c):
    out = []
    for side in (1, -1):
        at = [c.at[k] for k in c.at if c.side[k] == side]
        out.append(_segments(c.offset_line(min(at), max(at), side * c.rail)))
    return out


def test_the_pieces_are_laid_end_to_end_on_the_curve():
    """1IHF's DNA 1: build_layout sets its pieces end to end on one curve (na_curve.lay_out), so the curve drawn
    is the one it laid out: each piece centred on it along its direction, no bare rail between two pieces (no step
    between consecutive pairs longer than the longest a piece has itself, or the duplex's own step at a joint)."""
    lay, _ = _routed("1IHF")
    (c,) = lay.na_curves
    pieces = [lay.stems[k] for k in c.pieces]
    assert c.held == [(d.cx, d.cy, d.angle) for d in pieces]
    for d in pieces:
        keys = [(a, i) for a, i, _, _ in d.pairs]
        mid = (min(c.at[k] for k in keys) + max(c.at[k] for k in keys)) / 2
        p, t, _ = (v[0] for v in c.frames(mid))
        assert float(np.hypot(p[0] - d.cx, p[1] - d.cy)) < 0.05, d.name
        assert abs(math.atan2(t[1], t[0]) - d.angle) < 0.01, d.name
    own = max(float(np.diff(sorted(d.axial[(a, i)] for a, i, _, _ in d.pairs)).max()) for d in pieces)
    us = np.array(sorted(c.at[(a, i)] for a, i, _, _ in c.pairs))
    assert float(np.diff(us).max()) <= max(own, DNA_RISE * SCALE) + 1e-6


def test_the_protein_is_laid_out_round_the_curve():
    """1IHF: no protein element lies on the curve's band (rails and letters), no protein loop crosses its rails, and
    chain A's β-arm (A:60-63, A:70-73) stands on the inside of the bend, as in the structure, beside the band."""
    lay, loops = _routed("1IHF")
    (c,) = lay.na_curves
    band = _band_points(c, lay)
    assert not [sid for sid, e in lay.placed.items() if _under(band, e)]
    rails = _rails(c)
    assert sum(_crossings(r, _segments(lp.points)) for r in rails for lp in loops if len(lp.points) > 1) == 0
    lo, hi = c.span()
    s = np.arange(lo, hi, 0.1)
    pts, tan, nrm = c.frames(s)
    for sid in ("A:60-63", "A:70-73"):
        e = lay.placed[sid]
        k = int(np.argmin(np.hypot(pts[:, 0] - e.cx, pts[:, 1] - e.cy)))
        t0, t1 = c.frames(s[k] - 1.5)[1][0], c.frames(s[k] + 1.5)[1][0]
        bend = float(t0[0] * t1[1] - t0[1] * t1[0])
        side = float((np.array([e.cx, e.cy]) - pts[k]) @ nrm[k])
        assert side * bend > 0, sid  # on the side the curve turns toward
        assert float(np.hypot(e.cx - pts[k][0], e.cy - pts[k][1])) < c.band, sid


def test_the_nucleosome_wraps_its_octamer():
    """1KX5: the gyres grown round the octamer (every protein element inside the inner gyre, none on the band, no
    loop across a rail), the rungs spread evenly over them (each step between pairs the same multiple of its own)."""
    lay, loops = _routed("1KX5")
    (c,) = lay.na_curves
    assert c.gyres and not c.overlap
    band = _band_points(c, lay)
    assert not [sid for sid, e in lay.placed.items() if _under(band, e)]
    rails = _rails(c)
    assert sum(_crossings(r, _segments(lp.points)) for r in rails for lp in loops if len(lp.points) > 1) == 0
    mid = (c.xy.min(axis=0) + c.xy.max(axis=0)) / 2  # the ring's middle
    inner = float(np.hypot(c.xy[:, 0] - mid[0], c.xy[:, 1] - mid[1]).min())
    for sid, e in lay.placed.items():
        assert math.hypot(e.cx - mid[0], e.cy - mid[1]) < inner, sid
    pieces = [lay.stems[k] for k in c.pieces]
    for d in pieces:  # each piece's pairs as far apart on the curve as in the piece (stretched alike)
        keys = sorted((a, i) for a, i, _, _ in d.pairs)
        on, own = [c.at[k] - c.at[keys[0]] for k in keys], [d.axial[k] - d.axial[keys[0]] for k in keys]
        assert np.allclose(on, own, atol=1e-3)
    steps = [min(c.at[k] for k in e.axial) - max(c.at[k] for k in d.axial) for d, e in zip(pieces, pieces[1:])]
    assert max(steps) - min(steps) < 0.01


def test_pieces_lay_out_nothing_on_a_curve():
    """The default (pieces): build_layout lays no curve out and pins nothing."""
    lay, _, _ = make_layout(_path("1IHF"), look=Style())
    assert lay.na_curves == []
