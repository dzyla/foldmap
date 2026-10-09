import xml.etree.ElementTree as ET

import pytest
from helpers import pipeline, tripeptide

from topoplot.cli import make_figure
from topoplot.palette import chain_colors
from topoplot.render import save


def svg_ids(path):
    root = ET.parse(path).getroot()
    return [e.get("id") for e in root.iter() if e.get("id")]


def test_palette_is_stable_and_distinct():
    colors = chain_colors(["A", "B", "C", "D"])
    assert list(colors) == ["A", "B", "C", "D"]
    assert len(set(colors.values())) == 4
    assert all(c.startswith("#") and len(c) == 7 for c in colors.values())
    assert chain_colors(["A", "B"]) == {k: colors[k] for k in "AB"}  # same chain, same colour


@pytest.mark.parametrize("mode", ["projected", "stack"])
def test_svg_has_one_named_element_per_strand_helix_and_loop(zs5, tmp_path, mode):
    _, sses, _ = pipeline(zs5)
    out = tmp_path / "f.svg"
    save(make_figure(zs5, mode=mode), out)
    ids = svg_ids(out)
    assert sum(i.startswith("strand:") for i in ids) == sum(s.kind == "E" for s in sses)
    assert sum(i.startswith("helix:") for i in ids) == sum(s.kind == "H" for s in sses)
    pairs = sum(1 for a, b in zip(sses, sses[1:]) if a.chain == b.chain)
    assert sum(i.startswith("loop:") for i in ids) == pairs


def test_helix_and_text_stay_editable_in_svg(ubq, tmp_path):
    out = tmp_path / "u.svg"
    save(make_figure(ubq), out)
    text = out.read_text()
    assert "<text" in text  # fonts are not converted to outlines
    for label in (">N<", ">C<", ">A<", "α1"):
        assert label in text
    assert "helix:A:23-34" in text


def test_png_and_pdf_are_written(ubq, tmp_path):
    fig = make_figure(ubq)
    save(fig, tmp_path / "u.png")
    save(fig, tmp_path / "u.pdf")
    assert (tmp_path / "u.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert (tmp_path / "u.pdf").read_bytes()[:5] == b"%PDF-"


def test_unknown_extension_rejected(ubq, tmp_path):
    with pytest.raises(ValueError, match="extension"):
        save(make_figure(ubq), tmp_path / "u.tiff")


def test_identical_runs_give_identical_svg(zs5, tmp_path):
    for n in ("a", "b"):
        save(make_figure(zs5), tmp_path / f"{n}.svg")
    assert (tmp_path / "a.svg").read_bytes() == (tmp_path / "b.svg").read_bytes()


def test_structure_without_secondary_structure_renders_a_note(ubq, tmp_path):
    out = tmp_path / "t.svg"
    save(make_figure(tripeptide(ubq, tmp_path)), out)
    assert "no helices or strands" in out.read_text()


def test_title_is_drawn(ubq, tmp_path):
    out = tmp_path / "u.svg"
    save(make_figure(ubq, title="Ubiquitin 1UBQ"), out)
    assert "Ubiquitin 1UBQ" in out.read_text()


def test_barrel_ghost_strand_is_drawn_and_named(tmp_path, ubq):
    out = tmp_path / "g.svg"
    save(make_figure(ubq.parent / "1EMA.cif"), out)
    assert sum(i.startswith("ghost:") for i in svg_ids(out)) == 1
    save(make_figure(ubq), out)
    assert not any(i.startswith("ghost:") for i in svg_ids(out))


def test_named_palettes_are_distinct_and_unknown_names_fail():
    from topoplot.palette import PALETTES

    assert {"okabe-ito", "tol-bright", "tol-muted", "greys"} <= set(PALETTES)
    for name in PALETTES:
        cols = chain_colors(["A", "B", "C"], name)
        assert len(set(cols.values())) == 3, name
    assert chain_colors(["A"]) == chain_colors(["A"], "okabe-ito")
    with pytest.raises(ValueError, match="palette"):
        chain_colors(["A"], "rainbow")


def _artist(fig, gid):
    for a in [*fig.axes[0].patches, *fig.axes[0].texts]:
        if a.get_gid() == gid:
            return a
    raise KeyError(gid)


def test_helix_labels_sit_outside_the_helix(ubq):
    import numpy as np

    from topoplot.cli import make_layout
    from topoplot.render import draw
    from topoplot.route import route_loops

    lay, sses, bb = make_layout(ubq.parent / "1EMA.cif")
    fig = draw(lay, route_loops(lay, sses, bb), sses)
    for p in lay.placed.values():
        if p.sse.kind == "H":
            x, y = _artist(fig, f"label:{p.sse.id}").get_position()
            d = np.array([x - p.cx, y - p.cy])
            along = d @ p.direction
            across = abs(d @ np.array([-p.direction[1], p.direction[0]]))
            outside = abs(along) > p.length / 2 or across > p.width / 2  # beyond the rounded bar
            assert outside and np.hypot(*d) > p.width / 2, p.sse.id


@pytest.mark.parametrize("style", ["bold", "pale", "outline"])
def test_styles_render_every_element(zs5, tmp_path, style):
    _, sses, _ = pipeline(zs5)
    out = tmp_path / f"{style}.svg"
    save(make_figure(zs5, style=style), out)
    ids = svg_ids(out)
    assert sum(i.startswith("strand:") for i in ids) == sum(s.kind == "E" for s in sses)


@pytest.mark.parametrize("style", ["bold", "pale", "outline"])
def test_no_style_draws_hatched_bars(ubq, tmp_path, style):
    out = tmp_path / f"{style}.svg"
    save(make_figure(ubq, style=style), out)
    assert "<pattern" not in out.read_text()


def _polys(patch):
    return [p for p in patch.get_path().to_polygons(closed_only=False) if len(p) > 2]


def test_helix_is_a_coil_built_from_one_repeated_turn(ubq):
    import numpy as np

    from topoplot.cli import make_layout

    lay, sses, _ = make_layout(ubq)
    fig = make_figure(ubq)
    p = lay.placed["A:23-34"]  # 12 residues: about 3 turns of 3.6
    front = [poly @ np.array([[np.cos(-p.angle), np.sin(-p.angle)], [-np.sin(-p.angle), np.cos(-p.angle)]])
             for poly in _polys(_artist(fig, "helix:A:23-34"))]  # page -> helix frame (rotation only)
    back = _polys(_artist(fig, "helix-back:A:23-34"))
    turns = len(back)
    assert turns == 3 and len(front) == turns + 1  # a half front band at each end
    full = front[1:-1]
    for a in full[1:]:  # every full turn is the same fragment, only shifted along the axis
        shift = a.mean(axis=0) - full[0].mean(axis=0)
        assert np.allclose(a - shift, full[0], atol=1e-6)
        assert abs(shift[1]) < 1e-6 and shift[0] > 0
    span = np.vstack(front) @ np.array([1.0, 0.0])
    assert span.max() - span.min() == pytest.approx(p.length, abs=0.6)  # the coil fills the element


def test_coil_starts_and_ends_on_the_ports(ubq):
    import numpy as np

    from topoplot.cli import make_layout

    lay, _, _ = make_layout(ubq)
    fig = make_figure(ubq)
    p = lay.placed["A:23-34"]
    pts = np.vstack(_polys(_artist(fig, "helix:A:23-34")))
    for port in (p.n_port, p.c_port):
        assert np.min(np.hypot(*(pts - port).T)) < 0.25  # a loop meets the ribbon, not empty space


def test_outline_style_leaves_strands_white(ubq):
    from matplotlib.colors import to_hex

    fig = make_figure(ubq, style="outline")
    strand = _artist(fig, "strand:A:2-7")
    assert to_hex(strand.get_facecolor()) == "#ffffff"
    assert to_hex(strand.get_edgecolor()) != "#ffffff"


def test_unknown_style_is_rejected(ubq):
    with pytest.raises(ValueError, match="style"):
        make_figure(ubq, style="neon")


def test_curved_loops_are_smooth_and_still_end_at_the_ports(ubq):
    import numpy as np

    from matplotlib.path import Path as MPath

    from topoplot.cli import make_layout

    fig = make_figure(ubq, loops="curved")
    lay, sses, _ = make_layout(ubq)
    for a, b in zip(sses, sses[1:]):
        path = _artist(fig, f"loop:{a.id}>{b.id}").get_path()
        assert np.allclose(path.vertices[0], lay.placed[a.id].c_port)
        assert np.allclose(path.vertices[-1], lay.placed[b.id].n_port)
    codes = np.concatenate([_artist(fig, f"loop:{a.id}>{b.id}").get_path().codes for a, b in zip(sses, sses[1:])])
    assert (codes == MPath.CURVE4).sum() > (codes == MPath.CURVE3).sum()  # cubic smoothing, not rounded corners
    with pytest.raises(ValueError, match="loop"):
        make_figure(ubq, loops="zigzag")


@pytest.mark.parametrize("pdb", ["1EMA", "2HHB", "1TIM"])
def test_curved_corners_stay_clear_of_elements(ubq, pdb):
    import numpy as np
    from matplotlib.path import Path as MPath

    from topoplot.cli import make_layout
    from topoplot.render import draw
    from topoplot.route import route_loops

    lay, sses, bb = make_layout(ubq.parent / f"{pdb}.cif")
    loops = route_loops(lay, sses, bb)
    fig = draw(lay, loops, sses, loop_style="curved")
    rects = [p.rect for p in [*lay.placed.values(), *lay.ghosts]]
    t = np.linspace(0, 1, 25)[:, None]
    for loop in loops:
        if loop.fallback:
            continue
        path = _artist(fig, f"loop:{loop.a_id}>{loop.b_id}").get_path()
        v, codes = path.vertices, path.codes
        for k in np.flatnonzero(codes == MPath.CURVE4)[::3]:  # each cubic: v[k-1] start, v[k], v[k+1] controls, v[k+2] end
            p0, p1, p2, p3 = v[k - 1], v[k], v[k + 1], v[k + 2]
            pts = (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3
            for x0, y0, x1, y1 in rects:
                inside = (pts[:, 0] > x0 + 1e-6) & (pts[:, 0] < x1 - 1e-6) & (pts[:, 1] > y0 + 1e-6) & (pts[:, 1] < y1 - 1e-6)
                assert not inside.any(), (loop.a_id, loop.b_id)


def test_every_corner_of_a_routed_loop_is_rounded(ubq):
    import numpy as np
    from matplotlib.path import Path as MPath

    from topoplot.cli import make_layout
    from topoplot.render import draw
    from topoplot.route import route_loops

    lay, sses, bb = make_layout(ubq.parent / "1TIM.cif")
    loops = route_loops(lay, sses, bb)
    fig = draw(lay, loops, sses)
    for loop in loops:
        if loop.fallback:
            continue
        codes = _artist(fig, f"loop:{loop.a_id}>{loop.b_id}").get_path().codes
        assert (codes == MPath.CURVE3).sum() == 2 * (len(loop.points) - 2), (loop.a_id, loop.b_id)


def test_legend_fits_a_narrow_many_chain_figure():
    import numpy as np
    from helpers import fake_sse

    from topoplot.frame import Frame
    from topoplot.layout import build_layout
    from topoplot.render import draw

    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    sses = [fake_sse("H", 1, 20, centroid=(0, 0, 0), axis=(0, 1, 0), chain=c) for c in "ABCDEF"]  # coiled-coil-like
    lay = build_layout(sses, [], f, "stack")
    fig = draw(lay, [], sses)
    fig.canvas.draw()
    page = fig.bbox
    legend = [t for t in fig.axes[0].texts if t.get_text().startswith("Chain ")]
    assert len(legend) == 6
    for t in legend:
        bb = t.get_window_extent()
        assert page.x0 <= bb.x0 and bb.x1 <= page.x1 and page.y0 <= bb.y0 and bb.y1 <= page.y1, t.get_text()


def test_each_sheet_sits_on_its_own_panel(ubq):
    from topoplot.cli import make_layout
    from topoplot.style import Style

    path = ubq.parent / "5NKT.cif"  # FimA: a two-sheet sandwich
    lay, sses, _ = make_layout(path)
    fig = make_figure(path)
    panels = _gids(fig, "sheet-panel:")
    assert len(panels) == len(lay.sheets) == 2
    for panel, ids in zip(panels, lay.sheets):
        box = panel.get_window_extent().transformed(fig.axes[0].transData.inverted())
        for k in ids:
            x0, y0, x1, y1 = lay.placed[k].rect
            assert box.x0 <= x0 and x1 <= box.x1 and box.y0 <= y0 and y1 <= box.y1, k
    assert not _gids(make_figure(path, look=Style(sheet_panels=False)), "sheet-panel:")


def test_text_uses_a_journal_font(ubq, tmp_path):
    out = tmp_path / "f.svg"
    save(make_figure(ubq), out)
    import re

    families = re.findall(r"font-family: '?([^',;]+)", out.read_text())
    assert families and all(f in ("Arial", "Helvetica", "Liberation Sans") for f in families)  # first choice


def _gids(fig, prefix):
    return [a for a in [*fig.axes[0].patches, *fig.axes[0].texts] if (a.get_gid() or "").startswith(prefix)]
