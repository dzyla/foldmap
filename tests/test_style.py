import pytest

from foldmap.style import LAYOUT_KEYS, THEME_KEYS, THEMES, Style, resolve_style

PRESETS = THEMES  # older name


def test_presets_cover_the_common_uses_and_are_valid():
    assert {"publication", "minimal", "print", "presentation", "cartoon"} <= set(PRESETS)
    for s in PRESETS.values():
        assert isinstance(s, Style)
        s.validate()


def test_default_is_publication_without_gloss_and_with_axis_snapped_helices():
    s = resolve_style()
    assert s == PRESETS["publication"]
    assert s.helix_shading == "depth" and s.helix_angle == "snap" and s.fill == "bold"


def test_overrides_apply_on_top_of_a_preset():
    s = resolve_style("minimal", overrides=["helix_scale=0.8", "loops=curved", "labels=false"])
    assert s.helix_scale == 0.8 and s.loops == "curved" and s.labels is False
    assert s.fill == PRESETS["minimal"].fill


@pytest.mark.parametrize(
    "bad, msg",
    [
        ("colour=red", "helix_scale"),
        ("helix_scale=big", "number"),
        ("loops=zigzag", "orthogonal"),
        ("helix_scale=-1", "positive"),
        ("helix_scale", "key=value"),
    ],
)
def test_bad_overrides_say_what_is_allowed(bad, msg):
    with pytest.raises(ValueError, match=msg):
        resolve_style(overrides=[bad])


def test_unknown_preset_is_rejected():
    with pytest.raises(ValueError, match="publication"):
        resolve_style("fancy")


def test_style_file_can_name_a_preset_and_override_it(tmp_path):
    f = tmp_path / "s.yaml"
    f.write_text("preset: print\nhelix_scale: 0.7\nfont_scale: 1.2\n")
    s = resolve_style(style_file=f, overrides=["font_scale=1.5"])
    assert s.palette == PRESETS["print"].palette and s.helix_scale == 0.7 and s.font_scale == 1.5


def test_style_file_with_unknown_key_fails(tmp_path):
    f = tmp_path / "s.yaml"
    f.write_text("helix_size: 2\n")
    with pytest.raises(ValueError, match="helix_size"):
        resolve_style(style_file=f)


# --- the style reaching the figure ---------------------------------------------------------------
import numpy as np  # noqa: E402

from foldmap.cli import make_figure, make_layout  # noqa: E402
from helpers import pipeline  # noqa: E402


def _artists(fig, prefix):
    ax = fig.axes[0]
    return [a for a in [*ax.patches, *ax.texts, *ax.lines] if (a.get_gid() or "").startswith(prefix)]


def test_no_glint_unless_gloss(ubq):
    from dataclasses import replace

    assert not _artists(make_figure(ubq), "helix-shine:")
    assert _artists(make_figure(ubq, look=replace(Style(), helix_shading="gloss")), "helix-shine:")


def test_shading_none_gives_the_back_face_the_front_colour(ubq):
    fig = make_figure(ubq, look=Style(helix_shading="none"))
    (front,), (back,) = _artists(fig, "helix:"), _artists(fig, "helix-back:")
    assert tuple(front.get_facecolor()) == tuple(back.get_facecolor())


@pytest.mark.parametrize("scale", [0.6, 1.5])
def test_element_scales_change_drawn_widths(zs5, scale):
    from foldmap.layout import DNA_W, HELIX_W, PITCH

    look = Style(helix_scale=scale, strand_scale=scale, dna_scale=scale)
    lay, sses, _ = make_layout(zs5.parent / "1LMB.cif", look=look)
    assert all(p.width == pytest.approx(HELIX_W * scale) for p in lay.placed.values() if p.sse.kind == "H")
    assert lay.dna[0].width == pytest.approx(DNA_W * scale)
    lay, sses, _ = make_layout(zs5, look=look)
    sheet = sorted((p for p in lay.placed.values() if p.sse.kind == "E"), key=lambda p: (round(p.cy, 3), p.cx))
    gaps = [b.cx - a.cx for a, b in zip(sheet, sheet[1:]) if abs(a.cy - b.cy) < 1e-9]
    assert any(abs(g - PITCH * scale) < 1e-6 for g in gaps)


def test_line_and_font_scales(ubq):
    base, big = make_figure(ubq), make_figure(ubq, look=Style(loop_width=2.0, font_scale=1.5))
    loop = lambda f: _artists(f, "loop:")[0].get_linewidth()  # noqa: E731
    font = lambda f: _artists(f, "label:")[0].get_fontsize()  # noqa: E731
    assert loop(big) == pytest.approx(2 * loop(base))
    assert font(big) == pytest.approx(1.5 * font(base))


def test_labels_and_legend_can_be_switched_off(ubq):
    fig = make_figure(ubq, look=Style(labels=False, legend=False))
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert not _artists(fig, "label:") and not any(t.startswith("Chain ") for t in texts)


@pytest.mark.parametrize("pdb", ["2LZM", "1TIM"])
def test_helix_angle_modes(ubq, pdb):
    path = ubq.parent / f"{pdb}.cif"
    snapped = make_layout(path)[0]
    upright = make_layout(path, look=Style(helix_angle="upright"))[0]
    tilted = make_layout(path, look=Style(helix_angle="tilted"))[0]
    helices = [k for k, p in snapped.placed.items() if p.sse.kind == "H"]
    q = lambda a: abs(np.sin(2 * a))  # noqa: E731  # 0 on the axes
    assert all(q(snapped.placed[k].angle) < 1e-9 for k in helices)
    assert any(abs(np.cos(snapped.placed[k].angle)) > 0.5 for k in helices)  # some lie horizontal
    assert all(abs(np.cos(upright.placed[k].angle)) < 1e-9 for k in helices)
    assert any(q(tilted.placed[k].angle) > 0.2 for k in helices)


def test_snapped_helix_keeps_the_nearest_axis_direction(ubq):
    from foldmap.frame import view_frame

    path = ubq.parent / "2LZM.cif"
    _, sses, _ = pipeline(path)
    f = view_frame(sses)
    lay = make_layout(path)[0]
    for s in sses:
        if s.kind == "H":
            d = np.array([s.axis @ f.u, s.axis @ f.v])
            if np.linalg.norm(d) > 0.3:
                got = lay.placed[s.id].direction
                assert got @ d / np.linalg.norm(d) >= np.sqrt(0.5) - 1e-9, s.id  # within 45° of the true tilt


def test_cli_preset_set_and_style_file(ubq, tmp_path):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(ubq), "-o", str(out), "--theme", "minimal", "--set", "helix_scale=0.8"]) == 0
    f = tmp_path / "look.yaml"
    f.write_text("theme: print\nlabels: false\n")
    assert main(["plot", str(ubq), "-o", str(out), "--style-file", str(f)]) == 0
    assert 'id="label:' not in out.read_text()  # no element labels (legend ids may contain the word)
    assert main(["plot", str(ubq), "-o", str(out), "--set", "helix_scale=big"]) == 1


def test_cli_lists_presets(capsys):
    from foldmap.cli import main

    assert main(["styles"]) == 0
    out = capsys.readouterr().out
    assert all(name in out for name in PRESETS) and "helix_scale" in out


def test_themes_only_change_how_things_look_not_where_they_go():
    default = Style()
    for name, theme in THEMES.items():
        for key in LAYOUT_KEYS:
            assert getattr(theme, key) == getattr(default, key), (name, key)
    assert "helix_angle" in LAYOUT_KEYS and "fill" in THEME_KEYS


@pytest.mark.parametrize("old, new", [("classic", "bold"), ("flat", "pale")])
def test_old_fill_names_still_work(old, new):
    assert resolve_style(overrides=[f"fill={old}"]).fill == new


def test_styles_help_separates_themes_theme_keys_and_layout_keys(capsys):
    from foldmap.cli import main

    assert main(["styles"]) == 0
    out = capsys.readouterr().out
    assert "Themes" in out and "Theme keys" in out and "Layout keys" in out


def _hex(c):
    from matplotlib.colors import to_hex

    return to_hex(c)


def test_sequence_colouring_runs_n_to_c_along_each_chain(ubq):

    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(ubq)
    fig = make_figure(ubq, look=Style(color_by="sequence"))
    cols = [
        _hex(_artists(fig, f"{ {'E': 'strand', 'G': 'eta'}.get(s.kind, 'helix') }:{s.id}")[0].get_facecolor())
        for s in sses
    ]
    assert len(set(cols)) == len(cols)  # every element its own shade
    from foldmap.render import sequence_colour

    ramp = Style().sequence_map
    assert cols[0] == sequence_colour(0.0, ramp) and cols[-1] == sequence_colour(1.0, ramp)  # N end, C end
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert any("N" in t and "C" in t and "→" in t for t in texts)  # the legend explains it


def test_sequence_colouring_restarts_for_each_chain(zs5):
    fig = make_figure(zs5, look=Style(color_by="sequence"))

    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(zs5)
    from foldmap.render import sequence_colour

    start = sequence_colour(0.0, Style().sequence_map)
    for chain in ("A", "D"):
        first = next(s for s in sses if s.chain == chain)
        art = _artists(fig, f"{ {'E': 'strand', 'G': 'eta'}.get(first.kind, 'helix') }:{first.id}")[0]
        assert _hex(art.get_facecolor()) == start, chain


def _lum(hexcol):
    from matplotlib.colors import to_rgb

    r, g, b = to_rgb(hexcol)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _element_colours(fig, sses):
    out = {}
    for s in sses:
        art = _artists(fig, f"{ {'E': 'strand', 'G': 'eta'}.get(s.kind, 'helix') }:{s.id}")
        if art:
            out[s.id] = _hex(art[0].get_facecolor())
    return out


def test_shade_colouring_darkens_each_chain_from_n_to_c(zs5):
    from foldmap.cli import make_layout

    _, sses, _ = make_layout(zs5)
    cols = _element_colours(make_figure(zs5, look=Style(color_by="shade")), sses)
    for chain in ("A", "D"):
        mine = [cols[s.id] for s in sses if s.chain == chain and s.id in cols]
        lum = [_lum(c) for c in mine]
        assert len(set(mine)) == len(mine) and all(a > b for a, b in zip(lum, lum[1:]))  # lighter at N, darker at C
    a = next(cols[s.id] for s in sses if s.chain == "A")
    d = next(cols[s.id] for s in sses if s.chain == "D")
    assert a != d  # chains keep their own hue


@pytest.mark.parametrize("mode", ["chain", "element"])
def test_loops_can_take_the_colour_of_their_chain_or_element(zs5, mode):
    from foldmap.cli import make_layout
    from foldmap.palette import chain_colors

    _, sses, _ = make_layout(zs5)
    fig = make_figure(zs5, look=Style(loop_color=mode))
    cols = _element_colours(fig, sses)
    chain = chain_colors(["A", "D"])
    by_id = {s.id: s for s in sses}
    for loop in _artists(fig, "loop:"):
        a = loop.get_gid().split(":", 1)[1].split(">")[0]
        from foldmap.render import _legible

        want = chain[by_id[a].chain] if mode == "chain" else cols[a]
        assert _hex(loop.get_edgecolor()) == _legible(want).lower()  # the colour, darkened if too pale for a line


def test_loop_arrows_point_from_n_to_c(ubq):
    from foldmap.cli import make_layout
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(ubq)
    loops = {f"{l.a_id}>{l.b_id}": l for l in route_loops(lay, sses, bb)}
    fig = make_figure(ubq, look=Style(loop_arrows=True))
    arrows = _artists(fig, "loop-arrow:")
    assert len(arrows) == len(loops)
    for art in arrows:
        loop = loops[art.get_gid().split(":", 1)[1]]
        tip, tail = art.get_path().vertices[1], (art.get_path().vertices[0] + art.get_path().vertices[2]) / 2
        pts = np.array(loop.points)
        seg = max(range(len(pts) - 1), key=lambda k: np.hypot(*(pts[k + 1] - pts[k])))  # the arrow's segment
        assert (tip - tail) @ (pts[seg + 1] - pts[seg]) > 0  # pointing on toward the C-terminal element
    assert not _artists(make_figure(ubq), "loop-arrow:")


def test_residue_numbers_mark_each_element_end(ubq):
    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(ubq)
    fig = make_figure(ubq, look=Style(residue_numbers=True))
    for s in sses:
        p = lay.placed[s.id]
        for end, port, num in (("N", p.n_port, s.first.seq), ("C", p.c_port, s.last.seq)):
            (t,) = _artists(fig, f"resnum:{s.id}:{end}")
            assert t.get_text() == str(num)
            assert np.hypot(*(np.array(t.get_position()) - port)) < p.width / 2 + 0.9
    assert not _artists(make_figure(ubq), "resnum:")


def test_coloured_loops_stay_dark_enough_to_follow(zs5):
    fig = make_figure(zs5, look=Style(color_by="shade", loop_color="element"))
    assert all(_lum(_hex(l.get_edgecolor())) <= 0.55 for l in _artists(fig, "loop:"))


@pytest.mark.parametrize("name, mode", [("1ZAA", "stack"), ("5NKT", "projected")])
def test_residue_numbers_never_sit_on_another_element(ubq, name, mode):
    from foldmap.cli import make_layout

    path = ubq.parent / f"{name}.cif"
    lay, sses, _ = make_layout(path, mode)
    fig = make_figure(path, mode, look=Style(residue_numbers=True))
    for t in _artists(fig, "resnum:"):
        own = t.get_gid().split(":")[1] + ":" + t.get_gid().split(":")[2]
        x, y = t.get_position()
        for k, p in lay.placed.items():
            if k != own:
                x0, y0, x1, y1 = p.rect
                assert not (x0 < x < x1 and y0 < y < y1), (t.get_gid(), k)


def _element_fills(fig, sses):
    return _element_colours(fig, sses)


def test_bfactor_colouring_runs_rigid_blue_to_flexible_red(ubq):
    from foldmap.cli import make_layout
    from foldmap.io import load_backbone
    from foldmap.render import property_colour

    lay, sses, _ = make_layout(ubq)
    bb = load_backbone(ubq)
    mean_b = {s.id: float(np.mean(bb.b[s.start : s.end + 1])) for s in sses}
    cols = _element_fills(make_figure(ubq, look=Style(color_by="bfactor")), sses)
    lo, hi = min(mean_b, key=mean_b.get), max(mean_b, key=mean_b.get)
    assert cols[lo] == property_colour(0.0, "bfactor") and cols[hi] == property_colour(1.0, "bfactor")
    texts = [t.get_text() for t in make_figure(ubq, look=Style(color_by="bfactor")).axes[0].texts]
    assert any("B-factor" in t for t in texts)


def test_hydropathy_colouring_follows_kyte_doolittle(ubq):
    from foldmap.cli import make_layout
    from foldmap.render import KYTE_DOOLITTLE, property_colour

    lay, sses, bb = make_layout(ubq)
    kd = {
        s.id: float(np.mean([KYTE_DOOLITTLE.get(bb.labels[r].name, 0.0) for r in range(s.start, s.end + 1)]))
        for s in sses
    }
    cols = _element_fills(make_figure(ubq, look=Style(color_by="hydropathy")), sses)
    assert cols[max(kd, key=kd.get)] == property_colour(1.0, "hydropathy")
    assert cols[min(kd, key=kd.get)] == property_colour(0.0, "hydropathy")


def test_sstype_colouring_gives_one_colour_per_kind(zs5):
    from foldmap.cli import make_layout

    _, sses, _ = make_layout(DATA_UBQ := zs5.parent / "2LZM.cif")
    cols = _element_fills(make_figure(DATA_UBQ, look=Style(color_by="sstype")), sses)
    helices = {cols[s.id] for s in sses if s.kind == "H"}
    strands = {cols[s.id] for s in sses if s.kind == "E"}
    assert len(helices) == 1 and len(strands) == 1 and helices != strands


@pytest.mark.parametrize("theme", ["flexibility", "hydropathy", "richardson", "goodsell", "journal"])
def test_new_themes_render(ubq, theme, tmp_path):
    from foldmap.render import save
    from foldmap.style import THEMES

    assert theme in THEMES
    save(make_figure(ubq, look=THEMES[theme]), tmp_path / f"{theme}.svg")


def test_pale_fill_reaches_strands_too(ubq):
    bold = make_figure(ubq)
    pale = make_figure(ubq, look=Style(fill="pale"))
    s = "strand:A:2-7"
    assert _lum(_hex(_artists(pale, s)[0].get_facecolor())) > _lum(_hex(_artists(bold, s)[0].get_facecolor())) + 0.15
