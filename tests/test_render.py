import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_figure
from foldmap.palette import chain_colors
from foldmap.render import save
from helpers import pipeline, tripeptide

DATA = Path(__file__).parent / "data"


def svg_ids(path):
    root = ET.parse(path).getroot()
    return [e.get("id") for e in root.iter() if e.get("id")]


def test_palette_is_stable_and_distinct():
    colors = chain_colors(["A", "B", "C", "D"])
    assert list(colors) == ["A", "B", "C", "D"]
    assert len(set(colors.values())) == 4
    assert all(c.startswith("#") and len(c) == 7 for c in colors.values())
    assert chain_colors(["A", "B"]) == {k: colors[k] for k in "AB"}  # same chain, same colour


def test_chains_past_the_palette_keep_distinct_colours():
    """random-100 C35: chains A and I drew one #0072B2 (5X4Z, 6JO6, 7Y7C): past the palette's first pass each chain
    takes a lightness step of the colour it repeats that stands CIEDE2000 >= 10 from every colour assigned, the
    white page and the helix back faces of its family; the first eight are the palette's, so figures of up to
    eight chains are unchanged."""
    from itertools import combinations

    from foldmap.palette import BACK_SHADE, CHAIN_DE, PALETTES, darken, delta_e

    chains = [f"c{k}" for k in range(24)]
    colours = list(chain_colors(chains).values())
    assert colours[:8] == PALETTES["okabe-ito"]
    assert min(delta_e(a, b) for a, b in combinations(colours, 2)) >= CHAIN_DE
    for k in range(8, 24):  # nor like a helix back face of its family, nor near the white page
        family = colours[k % 8 : k : 8]
        assert min(delta_e(colours[k], darken(c, BACK_SHADE)) for c in family) >= CHAIN_DE
        assert delta_e(colours[k], "#ffffff") >= CHAIN_DE
    assert chain_colors(chains) == dict(zip(chains, colours))  # deterministic
    assert len(set(chain_colors([f"c{k}" for k in range(61)]).values())) == 61  # 7JU4: every chain its own colour


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


def test_svg_has_lf_line_ends_on_every_platform(ubq, tmp_path, monkeypatch):
    """matplotlib writes an SVG given a path in text mode (CRLF on Windows), so every cross-platform SVG hash would
    differ: save() writes it from a buffer with newline='\\n', with or without provenance (<desc>)."""
    from matplotlib.figure import Figure

    targets = []
    real = Figure.savefig
    monkeypatch.setattr(Figure, "savefig", lambda self, fname, **kw: targets.append(fname) or real(self, fname, **kw))
    for n, items in (("with", None), ("without", [])):
        fig = make_figure(ubq)
        if items is not None:
            fig.foldmap_provenance = items
        save(fig, tmp_path / f"{n}.svg")
        data = (tmp_path / f"{n}.svg").read_bytes()
        assert b"\n" in data and b"\r" not in data
    assert targets and not any(isinstance(t, (str, Path)) for t in targets)  # never a path: no text-mode write


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
    from foldmap.palette import PALETTES

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

    from foldmap.cli import make_layout
    from foldmap.render import draw
    from foldmap.route import route_loops

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

    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(ubq)
    fig = make_figure(ubq)
    p = lay.placed["A:23-34"]  # 12 residues: about 3 turns of 3.6
    front = [
        poly @ np.array([[np.cos(-p.angle), np.sin(-p.angle)], [-np.sin(-p.angle), np.cos(-p.angle)]])
        for poly in _polys(_artist(fig, "helix:A:23-34"))
    ]  # page -> helix frame (rotation only)
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

    from foldmap.cli import make_layout

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

    from foldmap.cli import make_layout

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

    from foldmap.cli import make_layout
    from foldmap.render import draw
    from foldmap.route import route_loops

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
        for k in np.flatnonzero(codes == MPath.CURVE4)[
            ::3
        ]:  # each cubic: v[k-1] start, v[k], v[k+1] controls, v[k+2] end
            p0, p1, p2, p3 = v[k - 1], v[k], v[k + 1], v[k + 2]
            pts = (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3
            for x0, y0, x1, y1 in rects:
                inside = (
                    (pts[:, 0] > x0 + 1e-6)
                    & (pts[:, 0] < x1 - 1e-6)
                    & (pts[:, 1] > y0 + 1e-6)
                    & (pts[:, 1] < y1 - 1e-6)
                )
                assert not inside.any(), (loop.a_id, loop.b_id)


def test_every_corner_of_a_routed_loop_is_rounded(ubq):
    from matplotlib.path import Path as MPath

    from foldmap.cli import make_layout
    from foldmap.render import draw
    from foldmap.route import route_loops

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

    from foldmap.frame import Frame
    from foldmap.layout import build_layout
    from foldmap.render import draw
    from helpers import fake_sse

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
    from foldmap.cli import make_layout
    from foldmap.style import Style

    path = ubq.parent / "5NKT.cif"  # FimA: a two-sheet sandwich plus a small hairpin sheet
    lay, sses, _ = make_layout(path)
    fig = make_figure(path)
    panels = _gids(fig, "sheet-panel:")
    assert len(panels) == len(lay.sheets) == 3
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


@pytest.mark.parametrize("name", ["1LMB", "1ZAA", "1UBQ"])
def test_title_clears_everything_else(name):
    from pathlib import Path

    from matplotlib.backends.backend_agg import FigureCanvasAgg

    from foldmap.cli import make_figure

    title_text = "A title for the figure, as long as a real one (with its PDB id)"
    fig = make_figure(Path(__file__).parent / "data" / f"{name}.cif", title=title_text, assembly="asu")
    r = FigureCanvasAgg(fig).get_renderer()
    title = next(t for t in fig.axes[0].texts if t.get_text() == title_text)
    tb = title.get_window_extent(r)
    for t in fig.axes[0].texts:
        if t is not title and t.get_text().strip():
            assert not tb.overlaps(t.get_window_extent(r)), t.get_text()


def test_long_loops_and_tails_say_how_long_they_are():
    from pathlib import Path

    from foldmap.cli import make_figure, make_layout

    path = Path(__file__).parent / "data" / "AF-P04637.cif"  # p53: long disordered C-terminal tail
    lay, sses, _ = make_layout(path)
    fig = make_figure(path)
    texts = {(t.get_gid() or ""): t.get_text() for t in fig.axes[0].texts}
    tails = {k: v for k, v in texts.items() if k.startswith("tail-length:")}
    n_tail, c_tail = sses[0].start, len(lay.res_chain) - 1 - sses[-1].end
    assert tails["tail-length:N:A"] == f"{n_tail} aa" and tails["tail-length:C:A"] == f"{c_tail} aa"
    loops = [v for k, v in texts.items() if k.startswith("loop-length:")]
    assert all(int(v.split()[0]) >= 25 for v in loops)
    short = make_figure(Path(__file__).parent / "data" / "1UBQ.cif")
    assert not any((t.get_gid() or "").startswith(("tail-length:", "loop-length:")) for t in short.axes[0].texts)


@pytest.mark.parametrize("name", ["6Y7S", "8UTF"])
def test_legend_text_ends_before_the_next_swatch_on_a_wide_figure(name):
    """A wide figure draws legend text at its size (random-100 C4: 6 pt at least, the page growing past 7.2 in for
    it); the legend is spaced for that size, not a smaller one."""
    from pathlib import Path

    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = make_figure(Path(__file__).parent / "data" / f"{name}.cif", focus="none")
    r = FigureCanvasAgg(fig).get_renderer()
    swatches = [p.get_window_extent(r) for p in _gids(fig, "legend-chain:")]
    texts = [t.get_window_extent(r) for t in _gids(fig, "legend-chain-label:")]
    assert len(texts) > 4 and fig.get_size_inches()[0] > 7.2
    assert min(t.get_fontsize() for t in _gids(fig, "legend-chain-label:")) >= 6.0 - 1e-9
    for t in texts:
        row = [s for s in swatches if abs((s.y0 + s.y1) / 2 - (t.y0 + t.y1) / 2) < 2 and s.x0 > t.x0]
        if row:
            assert t.x1 <= min(s.x0 for s in row), (t, min(s.x0 for s in row))


def test_chromophore_is_marked_on_the_chain():
    fig = make_figure(DATA / "1EMA.cif")
    gids = {a.get_gid() for a in fig.axes[0].get_children()}
    assert "linked:A:66" in gids
    assert "CRO66" in {t.get_text() for t in fig.axes[0].texts}


def test_chromophore_star_is_easy_to_see():
    fig = make_figure(DATA / "1EMA.cif")
    ax = fig.axes[0]
    (star,) = [a for a in ax.get_children() if a.get_gid() == "linked:A:66"]
    v = star.get_xy()
    r = max(np.hypot(*(v - v.mean(axis=0)).T))
    assert r >= 0.35
    (label,) = [t for t in ax.texts if t.get_gid() == "linked-label:A:66"]
    assert label.get_fontsize() >= 7
    (swatch,) = [a for a in ax.get_children() if a.get_gid() == "legend-linked"]
    assert len(swatch.get_xy()) >= 10  # a star, not the plain square of the other entries


def test_chromophore_label_keeps_off_loops_and_ghosts():
    """1EMA: CRO66's label sat on the loop running up beside ghost strand H; its ring search now skips spots
    on a loop's path as well as on elements, ghosts and names."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    from foldmap.cli import make_figure_and_loops
    from foldmap.route import clip_len

    fig, loops = make_figure_and_loops(DATA / "1EMA.cif")
    FigureCanvasAgg(fig).draw()
    ax = fig.axes[0]
    (label,) = [t for t in ax.texts if t.get_gid() == "linked-label:A:66"]
    b = label.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
    box = (b.x0, b.y0, b.x1, b.y1)
    assert not [lp.a_id for lp in loops for p, q in zip(lp.points, lp.points[1:]) if clip_len(p, q, box) > 1e-9]
    ghosts = [a for a in ax.patches if (a.get_gid() or "").startswith("ghost:")]
    assert ghosts and not [g.get_gid() for g in ghosts if g.get_window_extent().overlaps(label.get_window_extent())]


@pytest.mark.parametrize("name", ["5NKT", "25IN", "6Y7S"])
def test_every_loop_lies_inside_the_page(name):
    """5NKT: a loop ran past the canvas edge and was clipped. Since the bulge join (Task 1) 5NKT's loops keep
    0.53 units of page; 25IN A:54-74>A:79-105 (0.13) and 6Y7S D:122-123>D:127-139 (0.19) still came this close."""
    fig = make_figure(DATA / f"{name}.cif")
    ax = fig.axes[0]
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    for p in ax.patches:
        if (p.get_gid() or "").startswith("loop:"):
            v = p.get_path().vertices
            assert v[:, 0].min() >= x0 + 0.2 and v[:, 0].max() <= x1 - 0.2, p.get_gid()
            assert v[:, 1].min() >= y0 + 0.2 and v[:, 1].max() <= y1 - 0.2, p.get_gid()


@pytest.mark.parametrize(
    "name, legend",
    [
        ("1UBQ.cif", set()),
        ("2HHB.cif", {"Chain A", "Chain B", "Chain C", "Chain D"}),
        ("fidelity/1TEN.cif.gz", set()),  # one chain each: no legend at all (COMPARE §2: wasted legends)
        ("fidelity/2ACY.cif.gz", set()),
        ("fidelity/1MBN.cif.gz", set()),
    ],
)
def test_chain_legend_only_when_there_is_more_than_one_chain(name, legend):
    ax = make_figure(DATA / name).axes[0]
    texts = {t.get_text() for t in ax.texts}
    assert {t for t in texts if t.startswith("Chain ")} == legend
    if not legend:  # nothing else to explain on these: no legend artist of any kind
        assert not [a for a in [*ax.patches, *ax.lines, *ax.texts] if (a.get_gid() or "").startswith("legend")]


def test_one_protein_chain_with_nucleic_acid_keeps_the_nucleic_acid_legend():
    """The single-chain rule drops only the "Chain X" box: the RNA backbone and H-bond swatches stay (1URN)."""
    from foldmap.style import Style

    fig = make_figure(DATA / "na" / "1URN.cif.gz", look=Style(hbonds="na"))
    gids = sorted(a.get_gid() for a in [*fig.axes[0].patches, *fig.axes[0].lines] if a.get_gid())
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert {"legend-na", "legend-hbond"} <= set(gids) and not any(g.startswith("legend-chain") for g in gids)
    assert "RNA (R)" in texts and any(t.startswith("H-bonds (") for t in texts)
    assert not any(t.startswith("Chain ") for t in texts)


def test_helix_over_a_sheet_is_drawn_on_top_see_through_with_its_loops_and_a_readable_name():
    """1UBQ α1 lies over the sheet (Task 14): drawn above the strands, its loops above the strands and below it, its
    name on a plate of the page colour, clear of every strand letter. Its colour is the strands' (one chain), so it
    takes a darker tone of it, fill at alpha 0.8, edge opaque, to stand apart from them (Task 14 review)."""
    from matplotlib.colors import to_hex

    from foldmap.style import Style

    fig = make_figure(DATA / "1UBQ.cif", look=Style(helix_on_sheet="auto"))
    ax = fig.axes[0]
    by = {a.get_gid(): a for a in [*ax.patches, *ax.texts] if a.get_gid()}
    helix, back, strand = by["helix:A:23-34"], by["helix-back:A:23-34"], by["strand:A:2-7"]
    assert to_hex(helix.get_facecolor()) != to_hex(strand.get_facecolor())
    for part in (helix, back):
        assert part.get_facecolor()[3] == pytest.approx(0.8) and part.get_edgecolor()[3] == 1.0
    assert helix.get_zorder() > back.get_zorder() > strand.get_zorder()
    for gid in ("loop:A:12-16>A:23-34", "loop:A:23-34>A:38-40"):
        assert strand.get_zorder() < by[gid].get_zorder() < helix.get_zorder(), gid
    assert by["loop:A:2-7>A:12-16"].get_zorder() < strand.get_zorder()  # every other loop stays below
    name = by["label:A:23-34"]
    assert name.get_bbox_patch() is not None and name.get_zorder() > by["label:A:2-7"].get_zorder()
    fig.canvas.draw()
    box = name.get_window_extent()
    for k, t in by.items():
        if k.startswith("label:") and k != "label:A:23-34":
            assert not box.overlaps(t.get_window_extent()), k


def test_helix_over_a_sheet_of_another_colour_is_see_through_at_0_6():
    """Coloured unlike its strands (sstype: helices and strands apart), the helix keeps its colour at alpha 0.6."""
    from matplotlib.colors import to_hex

    from foldmap.style import Style

    fig = make_figure(DATA / "1UBQ.cif", look=Style(helix_on_sheet="auto", color_by="sstype"))
    by = {a.get_gid(): a for a in fig.axes[0].patches if a.get_gid()}
    helix, strand = by["helix:A:23-34"], by["strand:A:2-7"]
    assert to_hex(helix.get_facecolor()) != to_hex(strand.get_facecolor())
    assert helix.get_alpha() == 0.6 and by["helix-back:A:23-34"].get_alpha() == 0.6


_COLOURS_240 = {  # sha256 of the 240 colours, space-joined, as chain_colors gave them before it kept them (0e51849)
    "okabe-ito": "cb016017057262915c31addcaf33fd2b6201b74327a9346d3072175cd8884f71",
    "tol-bright": "f13ab1542e15c6adce23a768d91afbb26beed5651199acc86523cd09c09bbd76",
    "tol-muted": "5e209d2210beefdbb382eb88f147986f2b92a9834b6a11055884855cc9b5ca50",
    "blueprint": "5b0e6294e1ea2ae95d129a3e9bbaeb84e3fa56f8852a17f4f8435e9e1f3ac9f1",
    "greys": "e409c154bd3c1f54c88295dd3390eead2c31a1b02146318dc0960f43f66bf846",
}


def test_chain_colours_are_worked_out_once_per_palette(monkeypatch):
    """Final review C (Minor 1): chain_colors compares each new colour with every one assigned, 3.7 s for 240 chains,
    and a figure asks for them two or three times. A palette's colours depend on the count only, each a prefix of
    the next, so they are kept: a second figure's call compares no colours again. The colours are the ones it gave
    before (6YF9's 180 chains and more)."""
    import hashlib

    from foldmap import palette

    chains = [f"c{k}" for k in range(240)]
    first = chain_colors(chains)
    compared = []
    real = palette.delta_e_lab
    monkeypatch.setattr(palette, "delta_e_lab", lambda a, b: compared.append(1) or real(a, b))
    assert chain_colors(chains[:200]) == {c: first[c] for c in chains[:200]} and chain_colors(chains) == first
    assert compared == []
    for name, digest in _COLOURS_240.items():
        got = chain_colors(chains, name)
        assert hashlib.sha256(" ".join(got.values()).encode()).hexdigest() == digest, name
