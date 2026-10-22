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
    (parent, ks), *_ = bent_duplexes(lay.stems, lay.nucleic)
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
    """The drawing: one ladder with gids of its own (na:curve:0:..., Task 3), no piece's rails or name, no joint (the
    joints taken out of the strands drawn), one name; the contact beads and the strands' 5'/3' labels with the
    nucleotides on the curve."""
    drawn = draw_spec(figure_spec(_path("1IHF"), style=["na_bend=curved"]))
    lay, ax = drawn.layout, drawn.figure.axes[0]
    c = lay.na_curves[0]
    gids = {a.get_gid() for a in [*ax.patches, *ax.lines, *ax.texts, *ax.collections] if a.get_gid()}
    assert "na:curve:0:rails" in gids and "na:curve:0:label" in gids
    assert not any(f"na:stem:{k}:" in g or g == f"na:label:{k}" for k in c.pieces for g in gids)
    assert not [p for p in lay.na_paths if p.kind == "joint" and p.between[0] in c.at]
    names = [t for t in ax.texts if t.get_gid() == "na:curve:0:label"]
    assert [t.get_text() for t in names] == ["DNA 1"]
    rungs = [a for a in ax.lines if a.get_gid() == "na:curve:0:rung"]
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
    want = {"ribbon": "na:curve:0:front", "block": "na:curve:0:block", "letters": "na:curve:0:rung"}[style]
    assert want in gids
    assert not any(g.startswith(f"na:stem:{k}:") for k in c.pieces for g in gids)


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


@pytest.mark.parametrize("pdb", ["1IHF", "1KX5"])  # an arc; gyres round the octamer (_gyres, _ring)
def test_the_curve_does_not_depend_on_the_hash_seed(pdb):
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
        run = subprocess.run([sys.executable, "-c", code, str(_path(pdb))], capture_output=True, text=True, env=env)
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
    """1IHF: no protein element lies on the curve's band (rails and letters), no protein loop crosses its rails (where
    the text is in the font the routing keep-offs are calibrated on: figcheck.ink_font_reason, as test_criteria), and
    chain A's β-arm (A:60-63, A:70-73) stands on the inside of the bend, as in the structure, beside the band."""
    drawn = _drawn_curved("1IHF")
    lay, loops = drawn.layout, drawn.loops
    (c,) = lay.na_curves
    band = _band_points(c, lay)
    assert not [sid for sid, e in lay.placed.items() if _under(band, e)]
    rails = _rails(c)
    if not _font_off(drawn):  # routes keep off text measured in the calibrated font
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
    loop across a rail where the text is in the calibrated font), the rungs spread evenly over them (each step between
    pairs the same multiple of its own)."""
    drawn = _drawn_curved("1KX5")
    lay, loops = drawn.layout, drawn.loops
    (c,) = lay.na_curves
    assert c.gyres and not c.overlap
    band = _band_points(c, lay)
    assert not [sid for sid, e in lay.placed.items() if _under(band, e)]
    rails = _rails(c)
    if not _font_off(drawn):  # routes keep off text measured in the calibrated font
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


# ---- Task 3: the checks, keep-outs and pages read the curve as drawn, never its pieces ----


def test_the_curve_body_as_drawn_is_one_band_of_straight_stretches():
    """layout.stem_bodies: 1IHF's curve as straight stretches (Chord) in place of its pieces, which no check or
    keep-out reads; together they cover the band between the rails from end to end (no gap where a turn opens one),
    every nucleotide on a stretch of its own curve, and reach no further than half a body width off it."""
    from foldmap.layout import stem_bodies, stems_at
    from foldmap.na_curve import Chord
    from foldmap.na_paths import in_body

    lay, _ = _routed("1IHF")
    (c,) = lay.na_curves
    bodies = stem_bodies(lay)
    chords = [b for b in bodies if isinstance(b, Chord)]
    assert chords and {b.curve for b in chords} == {0}
    assert not any(b is lay.stems[k] for k in c.pieces for b in bodies)
    lo, hi = c.extent()
    s = np.linspace(lo + 0.05, hi - 0.05, 300)
    for off in (-c.width / 2 + 0.02, 0.0, c.width / 2 - 0.02):
        for q in c.to_page(s, off):
            assert any(in_body(q, b) for b in chords), (off, q)
    for b in chords:  # no stretch reaches past the band's side
        corners = [b.cx, b.cy] + np.array(
            [x * b.length / 2 * b.direction + y * b.width / 2 * b.normal for x in (-1, 1) for y in (-1, 1)]
        )
        for q in corners:
            assert float(np.hypot(c.xy[:, 0] - q[0], c.xy[:, 1] - q[1]).min()) <= c.width / 2 * 1.2 + 0.05
    at = stems_at(lay)
    for key in c.at:
        assert isinstance(at[key], tuple) and any(in_body(c.point(key), b) for b in at[key]), key


def test_the_checks_read_the_curve_not_its_pieces():
    """figcheck: element overlaps count the curve once by its name (its stretches meeting at a turn are no overlap),
    the page's stem bodies hold the curve (na:curve:0) and no piece, the stem-name gaps and line gaps are the curve's
    name's, the stem names repeated are the curve's; layout.label_boxes and na_text_boxes keep text off the curve's
    name and letters, not the pieces'; the joints are no strand drawn (provenance's stuck reads no hidden one)."""
    from foldmap.figcheck import (
        _stem_boxes,
        element_overlaps,
        label_repeats,
        measure,
        stem_name_gaps,
        stem_name_line_gaps,
    )
    from foldmap.layout import _LETTER_EM, _ink_em, _ink_points, _reach, label_boxes, na_text_boxes

    drawn = _drawn_curved("1IHF")
    lay = drawn.layout
    (c,) = lay.na_curves
    pieces = {lay.stems[k].id for k in c.pieces}
    joined: set = set()
    names = [n for n, _ in _stem_boxes(lay, joined)]
    assert not pieces & set(names) and any(n.startswith("DNA 1@") for n in names)
    assert all(not any(x.startswith("DNA 1@") for x in pair) for pair in element_overlaps(lay))
    page = measure(drawn)
    stems = [g for g, _, _ in page.stems]
    assert "na:curve:0" in stems and not any(g == f"na:stem:{k}" for k in c.pieces for g in stems)
    gaps, lines = stem_name_gaps(drawn, page=page), stem_name_line_gaps(drawn, page=page)
    assert "na:curve:0:label" in gaps and "na:curve:0:label" in lines
    assert not any(f"na:label:{k}" in gaps for k in c.pieces)
    assert label_repeats(lay) == {}
    keep = dict(label_boxes(lay, drawn.sses))
    assert "stem-label:curve:0" in keep and not any(f"stem-label:{k}" in keep for k in c.pieces)
    em = _ink_em(lay, lay.font_scale)
    na = lay.nucleic
    key = next(iter(sorted(c.at)))
    letter = tuple(
        _reach(
            _ink_points(na.strands[key[0]].bases[key[1]], False, _LETTER_EM * em, c.angle(key))
            + c.point(key, lay.na_letter_rows[0])
        )
    )
    assert letter in na_text_boxes(lay)
    assert not [p for p in lay.na_paths if p.kind == "joint" and p.between and p.between[0] in c.at]


def test_the_na_metrics_measure_the_pieces_where_the_curve_draws_them():
    """nacheck reads each piece of a curve where its nucleotides are drawn (na_curve.drawn_pieces): centred on the
    curve, from its first rung to its last, its rails on the curve's sides; its body (na_curve.piece_bodies, for the
    contact gaps) the curve's stretches that hold its nucleotides as drawn."""
    from foldmap.na_curve import drawn_pieces, piece_bodies
    from foldmap.na_paths import in_body
    from foldmap.nacheck import shape

    lay, _ = _routed("1IHF")
    (c,) = lay.na_curves
    seen, bodies = drawn_pieces(lay), piece_bodies(lay)
    assert len(seen) == len(bodies) == len(lay.stems)
    for k in c.pieces:
        d = seen[k]
        assert set(d.axial) == set(lay.stems[k].axial) and d.id == lay.stems[k].id
        assert float(np.hypot(c.xy[:, 0] - d.cx, c.xy[:, 1] - d.cy).min()) < 0.05
        ends = sorted(d.axial, key=lambda key: c.at[key])
        a, b = c.to_page([c.at[ends[0]], c.at[ends[-1]]], [0.0, 0.0])
        assert abs(math.remainder(math.atan2(b[1] - a[1], b[0] - a[0]) - d.angle, 2 * math.pi)) < 1e-6
        for key in d.axial:
            assert d.side[key] == c.side[key]
            assert any(in_body(c.point(key), e) for e in bodies[k]), key
    assert shape(lay)["stems_drawn"] == len(lay.stems)


def test_off_the_bodies_and_boxes_are_the_pieces():
    """The default (pieces): every body and box the checks and keep-outs read is the layout's own stem."""
    from foldmap.layout import stem_bodies, stem_name_keeps, stem_points, stems_at

    lay, _ = _routed("1IHF", "pieces")
    assert stem_bodies(lay) is lay.stems
    assert stems_at(lay) == {key: d for d in lay.stems for key in d.axial}
    pts = stem_points(lay)
    assert all(np.array_equal(pts[key], d.point(key)) for d in lay.stems for key in d.axial)
    assert len(stem_name_keeps(lay)) == len(lay.stems) and lay.bent_duplex


def test_the_page_hovers_and_picks_the_curve_as_one():
    """The interactive page: one nucleic-acid entry for the curve (curve:0, its nucleotides and the residues touching
    them), none for its pieces; its ink's click target na:curve:0:hit."""
    import json

    from foldmap.interactive import _hit_targets, build_page

    page = build_page(_path("1IHF"), look=replace(Style(), na_bend="curved"))
    data = json.loads(page.split('<script type="application/json" id="topo-data">', 1)[1].split("</script>", 1)[0])
    keys = [x["key"] for x in data["links"]["nucleics"]]
    lay, _ = _routed("1IHF")
    (c,) = lay.na_curves
    assert "curve:0" in keys and not any(f"stem:{k}" in keys for k in c.pieces)
    entry = next(x for x in data["links"]["nucleics"] if x["key"] == "curve:0")
    assert entry["text"].startswith("DNA 1 · ") and len(entry["het"]) == len(c.at) and entry["residues"]
    assert 'id="na:curve:0:rails"' in page and 'id="na:curve:0:hit"' in page
    svg = '<g id="na:curve:2:rails">\n <path d="M 0 0 L 1 1"/>'
    assert 'id="na:curve:2:hit"' in _hit_targets(svg + '\n<g id="loop:A>B">\n <path d="M 0 0"/>')


def _drawn_curved(pdb, *style):
    return draw_spec(figure_spec(_path(pdb), style=["na_bend=curved", *style]))


def test_na_bend_is_listed_and_a_layout_file_keeps_it(tmp_path):
    """`foldmap styles` lists the key with its choices; a layout file saved from a curved figure keeps na_bend and
    redraws the same figure (A13)."""
    from foldmap import layoutfile
    from foldmap.cli import styles_help
    from foldmap.render import save
    from foldmap.spec import layout_document

    assert "na_bend           one of pieces, curved" in styles_help()
    spec = figure_spec(_path("1IHF"), style=["na_bend=curved"])
    drawn = draw_spec(spec)
    doc = layout_document(spec, drawn=drawn)
    assert doc["style"]["na_bend"] == "curved"
    keep = layoutfile.save(tmp_path / "ihf.yaml", doc)
    again = draw_spec(figure_spec(_path("1IHF"), layout_file=keep))
    assert again.layout.na_bend == "curved" and len(again.layout.na_curves) == 1
    save(drawn.figure, tmp_path / "a.svg")
    save(again.figure, tmp_path / "b.svg")
    assert (tmp_path / "a.svg").read_bytes() == (tmp_path / "b.svg").read_bytes()


def _font_off(drawn) -> str:
    """Why the routing keep-offs (text ink) are not the calibrated ones here: figcheck.ink_font_reason ('' on Linux
    with Liberation Sans; Arial on macOS and Windows)."""
    from foldmap.figcheck import ink_font_reason

    return ink_font_reason(drawn.figure)


# ---- fix round 1: the gate, two curves apart, the band walls, the name ----


def test_a_duplex_in_a_folded_rna_stays_in_pieces_with_a_note():
    """na_curve.curvable: a bent duplex one of whose chains pairs with itself in a stem outside the duplex's own helix
    (a duplex inside a folded RNA: 4YBB's rRNA stems) stays in pieces under na_bend=curved, the figure says so, and no
    choice is offered (Layout.bent_duplex); IHF's DNA 1 still curves (its chain C pairs in DNA 2 with chain D)."""
    from foldmap.na_curve import _pieces_of, curvable
    from foldmap.provenance import _na_bend

    lay, _ = _routed("4YBB")
    assert _pieces_of(lay.stems) and lay.na_curves == [] and not lay.bent_duplex
    (note,) = _na_bend(lay)
    assert note["kind"] == "na_pieces" and note["count"] == len(_pieces_of(lay.stems))
    assert "part of a folded RNA" in note["text"]
    ihf, _ = _routed("1IHF")
    (parent,) = _pieces_of(ihf.stems)
    assert curvable(ihf.nucleic, parent) and ihf.bent_duplex and len(ihf.na_curves) == 1
    assert _na_bend(ihf) == []


def _straight_curve(x0: float, y: float, band: float = 4.0):
    """A straight NACurve from (x0, y) ten units along x (no pieces), for _part."""
    from foldmap.na_curve import NACurve

    xs = np.arange(0.0, 10.05, 0.1)
    xy = np.column_stack([x0 + xs, np.full(len(xs), y)])
    c = NACurve("p", "DNA", [], xy, xs.copy(), np.tile([1.0, 0.0], (len(xs), 1)), 1.0, 2.0, {}, {}, [])
    c.band = band
    return c


def test_two_curves_are_parted_whole():
    """na_curve._part: each curve is laid out on its own and pinned, so the later one of two that meet is moved whole
    until its band keeps MARGIN off the earlier one's; one already clear stays."""
    from foldmap.layout import MARGIN
    from foldmap.na_curve import _part

    a, b, far = _straight_curve(0.0, 0.0), _straight_curve(3.0, 1.0), _straight_curve(0.0, 40.0)
    before = far.xy.copy()
    _part([a, b, far], [])
    gap = float(np.hypot(b.xy[:, None, 0] - a.xy[None, :, 0], b.xy[:, None, 1] - a.xy[None, :, 1]).min())
    assert gap >= (a.band + b.band) / 2 + MARGIN - 1e-6
    assert np.array_equal(far.xy, before) and np.array_equal(a.xy, _straight_curve(0.0, 0.0).xy)


def test_the_band_walls_settle():
    """Walls standing for a curve's band (Walls.band, as na_curve._band_boxes builds them: squares two across) are
    left along the band's normal, the whole strip at once: a footprint on a diagonal band is pushed clear in one go
    and stays clear (no staircase of squares to ping-pong on), and one beside the band but in a square's corner (and
    its MARGIN) is not met."""
    from foldmap.layout import MARGIN, Walls, _off_walls

    n = np.array([-1.0, 1.0]) / math.sqrt(2)  # a band along the diagonal through the origin, half width 1
    boxes, band = [], []
    for t in (-1.0, 0.0, 1.0):
        for side in (1.0, -1.0):
            c = np.array([t, t]) + side * 0.5 * n
            boxes.append([c[0] - 0.5, c[1] - 0.5, c[0] + 0.5, c[1] + 0.5])
            band.append([t, t, n[0], n[1], 1.0, n[1], -n[0], np.inf, np.inf])  # a band with no end near
    walls = Walls(np.array(boxes), np.ones((1, len(boxes)), bool), np.array(band))
    group = [np.array([0])]
    xy, halves, mids = np.array([[0.2, 0.4]]), np.array([[0.3, 0.3]]), np.zeros((1, 2))
    assert _off_walls(xy, halves, mids, walls, group)
    assert float(xy[0] @ n) >= 0.3 * math.sqrt(2) + 1.0 + MARGIN - 1e-9  # off the whole strip, along its normal
    assert not _off_walls(xy, halves, mids, walls, group, move=False)
    corner = np.array([0.0, 0.0]) + 0.5 * n + 1.1 * np.array([-1.0, 1.0])  # in the square's box with MARGIN only
    assert not _off_walls(corner[None, :], np.array([[0.01, 0.01]]), mids, walls, group, move=False)


@pytest.mark.parametrize("pdb", ["1IHF", "1KX5", "1RVA", "1AOI"])
def test_the_curve_name_sits_at_the_stem_name_margin(pdb):
    """The curve's name stands as a straight stem's does (layout._name_at): its ink _INK_MARGIN em beyond its letter
    row's (A17's stem-name gap, measured by figcheck in the curve's frame), clear of every element, other name and
    other stem; routes keep off its ink grown by that margin (name_keep, label_boxes route=True)."""
    from foldmap.criteria import stem_gap_limit
    from foldmap.figcheck import ink_font_reason, measure, stem_name_gaps
    from foldmap.layout import _INK_MARGIN, _ink_em, label_boxes

    drawn = _drawn_curved(pdb)
    lay = drawn.layout
    (c,) = lay.na_curves
    keep = dict(label_boxes(lay, drawn.sses, route=True))["stem-label:curve:0"]
    assert keep == c.name_keep
    g = _INK_MARGIN * _ink_em(lay, lay.font_scale)  # the name's ink: its keep-out less the margin
    x0, y0, x1, y1 = c.name_keep[0] + g, c.name_keep[1] + g, c.name_keep[2] - g, c.name_keep[3] - g
    for p in lay.placed.values():
        r = p.rect
        assert not (r[0] < x1 and x0 < r[2] and r[1] < y1 and y0 < r[3]), p.sse.id
    if ink_font_reason(drawn.figure):  # the gap is measured on Liberation Sans's ink
        return
    own, _ = stem_name_gaps(drawn, page=measure(drawn))["na:curve:0:label"]
    assert own <= stem_gap_limit(drawn.figure)
