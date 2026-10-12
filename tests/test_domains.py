"""Domain panels: named groups of elements drawn on a soft panel with their name."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main, make_figure, make_layout
from helpers import FID

DATA = Path(__file__).parent / "data"
UMOD = DATA / "6ZS5.cif"
DOMAINS = [("ZPN", ["res:A:331-440"]), ("ZPC", ["res:D:447-582"])]


def _gids(fig, prefix):
    return [a for a in [*fig.axes[0].patches, *fig.axes[0].texts] if (a.get_gid() or "").startswith(prefix)]


def test_domains_resolve_to_their_elements():
    lay, sses, _ = make_layout(UMOD, domains=DOMAINS)
    got = dict(lay.domains)
    assert set(got) == {"ZPN", "ZPC"}
    assert all(k.startswith("A:") for k in got["ZPN"]) and all(k.startswith("D:") for k in got["ZPC"])
    assert len(got["ZPN"]) + len(got["ZPC"]) == len(sses)


def test_each_domain_gets_a_named_panel_around_all_its_elements():
    lay, _, _ = make_layout(UMOD, domains=DOMAINS)
    fig = make_figure(UMOD, domains=DOMAINS)
    panels = {p.get_gid().split(":", 1)[1]: p for p in _gids(fig, "domain:")}
    names = {t.get_text() for t in _gids(fig, "domain-label:")}
    assert set(panels) == {"ZPN", "ZPC"} and names == {"ZPN", "ZPC"}
    inv = fig.axes[0].transData.inverted()
    for name, ids in lay.domains:
        box = panels[name].get_window_extent().transformed(inv)
        for k in ids:
            x0, y0, x1, y1 = lay.placed[k].rect
            assert box.x0 <= x0 + 1e-6 and x1 <= box.x1 + 1e-6 and box.y0 <= y0 + 1e-6 and y1 <= box.y1 + 1e-6


T4L = DATA / "2LZM.cif"
T4L_DOMAINS = [("N", ["res:A:13-59"]), ("C", ["res:A:1-12", "res:A:60-164"])]  # the classic lobes


def test_domains_are_kept_apart():
    lay, _, _ = make_layout(T4L, domains=T4L_DOMAINS)
    boxes = []
    for _, ids in lay.domains:
        r = np.array([lay.placed[k].rect for k in ids])
        boxes.append((r[:, 0].min(), r[:, 1].min(), r[:, 2].max(), r[:, 3].max()))
    (a, b) = boxes
    ox = min(a[2], b[2]) - max(a[0], b[0])
    oy = min(a[3], b[3]) - max(a[1], b[1])
    area = lambda r: (r[2] - r[0]) * (r[3] - r[1])  # noqa: E731
    assert max(ox, 0) * max(oy, 0) < 0.15 * min(area(a), area(b))  # panels barely overlap, if at all


def test_auto_domains_split_a_two_domain_protein():
    lay, sses, _ = make_layout(T4L, domains="auto")
    assert len(lay.domains) >= 2 and all(name.startswith("D") for name, _ in lay.domains)
    seen = [k for _, ids in lay.domains for k in ids]
    assert len(seen) == len(set(seen))


def test_domains_sharing_a_sheet_still_get_panels():  # UMOD: one sheet holds strands of both chains
    fig = make_figure(UMOD, domains=DOMAINS)
    assert len(_gids(fig, "domain:")) == 2


def test_auto_domains_leave_a_single_domain_protein_alone():
    lay, _, _ = make_layout(DATA / "1UBQ.cif", domains="auto")
    assert lay.domains == []


def test_unknown_domain_reference_is_an_error():
    with pytest.raises(ValueError, match="Q"):
        make_layout(UMOD, domains=[("X", ["Q"])])


def test_a_residue_range_without_res_is_an_error_naming_the_domain():
    """A:1-60 looks like an element id but is no element: a message, not an IndexError (Q-A19)."""
    with pytest.raises(ValueError, match=r"domain 'N'.*'A:1-60'.*res:A:1-60"):
        make_layout(DATA / "2LZM.cif", domains=[("N", ["A:1-60"])])


def test_a_domain_left_without_elements_is_an_error_naming_it():
    with pytest.raises(ValueError, match=r"domain 'C'.*no elements"):
        make_layout(DATA / "2LZM.cif", domains=[("N", ["α1"]), ("C", ["α1"])])


def test_element_ids_of_the_structure_pass_the_check():
    lay, _, _ = make_layout(T4L, domains=[("N", ["A:3-11", "A:14-19"]), ("C", ["A:60-79"])])
    assert lay.domains == [("N", ["A:3-11", "A:14-19"]), ("C", ["A:60-79"])]


def test_a_layout_file_naming_elements_the_structure_lacks_drops_them_with_a_warning(tmp_path, capsys):
    """A saved file outlives its structure: stale element ids are dropped as before, now named in a warning; the
    same id given with --domain is still an error."""
    lf = tmp_path / "t4l.yaml"
    lf.write_text("foldmap: 1\ndomains:\n  N: [A:3-11, A:999-1005]\n  C: [A:60-79]\n  X: [A:2000-2010]\n")
    out = tmp_path / "t4l.svg"
    assert main(["plot", str(T4L), "-o", str(out), "--layout-file", str(lf)]) == 0
    err = capsys.readouterr().err
    assert "foldmap: warning:" in err and "A:999-1005" in err and "A:2000-2010" in err
    svg = out.read_text()
    assert "domain:N" in svg and "domain:C" in svg and "domain:X" not in svg
    assert main(["plot", str(T4L), "-o", str(out), "--domain", "N=A:3-11,A:999-1005"]) == 1
    assert "'A:999-1005' is no element" in capsys.readouterr().err


def test_cli_domains_and_layout_file(tmp_path):
    import yaml

    out, lf = tmp_path / "u.svg", tmp_path / "u.yaml"
    assert (
        main(
            [
                "plot",
                str(UMOD),
                "-o",
                str(out),
                "--domain",
                "ZPN=res:A:331-440",
                "--domain",
                "ZPC=res:D:447-582",
                "--save-layout",
                str(lf),
            ]
        )
        == 0
    )
    assert "domain:ZPN" in out.read_text()
    doc = yaml.safe_load(lf.read_text())
    assert set(doc["domains"]) == {"ZPN", "ZPC"}
    again = tmp_path / "v.svg"
    assert main(["plot", str(UMOD), "-o", str(again), "--layout-file", str(lf)]) == 0
    assert again.read_bytes() == out.read_bytes()
    assert main(["plot", str(UMOD), "-o", str(out), "--domain", "ZPN"]) == 1  # no '='


def test_app_settings_carry_domains():
    from foldmap.app import layout_options

    opts = layout_options({"domains": "ZPN=res:A:331-440\nZPC=res:D:447-582"})
    assert opts["domains"] == DOMAINS


def test_app_layout_file_keeps_domains():
    from foldmap.app import layout_document, settings_from_document

    doc = layout_document(T4L, {"theme": "publication", "domains": "N=res:A:13-59\nC=res:A:1-12,res:A:60-164"})
    assert set(doc["domains"]) == {"N", "C"}
    back = settings_from_document(doc)
    assert back["domains"].splitlines()[0].startswith("N=res:A:")


def test_loops_keep_off_domain_names():
    from foldmap.layout import label_boxes
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(T4L, domains=T4L_DOMAINS)
    names = [b for k, b in label_boxes(lay, sses) if k.startswith("domain-label:")]
    assert len(names) == 2
    for loop in route_loops(lay, sses, bb):
        pts = np.asarray(loop.points, float)
        dense = np.vstack([a + t * (b - a) for a, b in zip(pts, pts[1:]) for t in np.linspace(0, 1, 40)])
        for x0, y0, x1, y1 in names:
            inside = (dense[:, 0] > x0) & (dense[:, 0] < x1) & (dense[:, 1] > y0) & (dense[:, 1] < y1)
            assert not inside.any()


def _meet(a, b) -> bool:
    """Do two (x0, y0, x1, y1) boxes overlap with some area?"""
    return min(a[2], b[2]) > max(a[0], b[0]) and min(a[3], b[3]) > max(a[1], b[1])


def test_a_ba_row_lies_in_one_domain_panel():
    """1LDM's Rossmann sheet is drawn as one βα row (one rigid unit). Auto domains split it (D1 the first βαβαβ, D3
    the second, α2 in D2), so three panels overlapped the row and each other; a row is never split now."""
    from foldmap.layout import domain_panel
    from helpers import FID

    lay, _, _ = make_layout(FID / "1LDM.cif.gz", focus="protomer", domains="auto")
    panels = {name: domain_panel(lay, name) for name, ids in lay.domains if ids}
    assert len(panels) >= 2 and lay.ba_rows
    names = sorted(panels)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            assert not _meet(panels[a], panels[b]), (a, b)
    owner = {k: name for name, ids in lay.domains for k in ids}
    for row in lay.ba_rows:
        assert len({owner.get(k) for k in row}) == 1, row
        for k in row:
            for name, box in panels.items():
                if name != owner.get(k):
                    assert not _meet(lay.placed[k].rect, box), (k, name)


def test_domains_colour_elements_and_draw_a_bar():
    from foldmap.cli import make_figure, make_layout
    from foldmap.render import _DOMAIN_TONES, element_colours
    from foldmap.style import Style

    lay, sses, _ = make_layout(DATA / "2LZM.cif", domains="auto")
    fig = make_figure(DATA / "2LZM.cif", domains="auto")
    gids = [a.get_gid() or "" for a in fig.axes[0].get_children()]
    bars = [g for g in gids if g.startswith("domain-bar:A:")]
    assert len(bars) == len(lay.domains) >= 2
    for look in (Style(color_by="domain"), Style()):  # chain colouring turns into domain colouring
        colours, _ = element_colours(lay, sses, look)
        for k, (_, ids) in enumerate(lay.domains):
            assert {colours[i] for i in ids} == {_DOMAIN_TONES[k]}


def test_other_colourings_win_over_domains():
    from foldmap.cli import make_layout
    from foldmap.render import _DOMAIN_TONES, element_colours
    from foldmap.style import Style

    lay, sses, _ = make_layout(DATA / "2LZM.cif", domains="auto")
    for look in (Style(color_by="sstype"), Style(domain_colours=False)):
        colours, _ = element_colours(lay, sses, look)
        assert not set(colours.values()) & set(_DOMAIN_TONES), look


CAS9_FILE = Path(__file__).parents[1] / "docs" / "domains" / "4OO8-cas9.txt"


def _flag_file(path):
    """The --domain flags of a docs/domains file, parsed as the CLI parses them."""
    from foldmap.spec import parse_domains

    words = " ".join(l for l in path.read_text().splitlines() if not l.startswith("#")).split()
    return parse_domains([v for k, v in zip(words, words[1:]) if k == "--domain"])


_LOOP_SLACK = 10  # residues: measured worst 10 (RuvC-III|PI, a boundary that falls in the loop from 1098 to PI's first
# element); a boundary at a loop's middle sits off the published one by up to half that loop


def _ranges(refs):
    return [tuple(int(v) for v in r.removeprefix("res:").split(":")[-1].split("-")) for r in refs]


def test_domain_bar_segments_end_at_the_domain_boundaries():
    """Cas9 (4OO8 chain A) with Nishimasu 2014 Fig 1A's eight domains: the bar reads RuvC-I, BH, REC1, REC2, REC1,
    RuvC-II, HNH, RuvC-III, PI, and each segment starts and ends at its domain range's ends, to within the half loop
    the bar splits between two domains (_LOOP_SLACK residues), except at the chain's own ends."""
    from foldmap.render import domain_bars

    domains = _flag_file(CAS9_FILE)
    assert [n for n, _ in domains] == ["RuvC-I", "BH", "REC1", "REC2", "RuvC-II", "HNH", "RuvC-III", "PI"]
    lay, _, _ = make_layout(DATA / "na" / "4OO8.cif.gz", chains=["A"], domains=domains)
    ((chain, first, last, segments),) = domain_bars(lay)
    assert chain == "A"
    names = [s.name for s in segments]
    assert names == ["RuvC-I", "BH", "REC1", "REC2", "REC1", "RuvC-II", "HNH", "RuvC-III", "PI"]
    assert [s.gid for s in segments][4] == "domain-bar:A:REC1#2"
    ranges = {n: _ranges(refs) for n, refs in domains}
    for seg in segments:
        assert seg.lo < seg.hi
        assert any(
            (seg.lo == first or abs(seg.lo - a) <= _LOOP_SLACK) and (seg.hi == last or abs(seg.hi - b) <= _LOOP_SLACK)
            for a, b in ranges[seg.name]
        ), (seg, ranges[seg.name])
    assert all(a.hi == b.lo for a, b in zip(segments, segments[1:]))  # no gap, no overlap


def test_a_domain_in_pieces_gets_a_segment_for_each():
    """T4 lysozyme's C lobe is residues 1-12 and 60-164: its two stretches are two segments, the second gid #2."""
    from foldmap.render import domain_bars

    lay, _, _ = make_layout(T4L, domains=T4L_DOMAINS)
    ((_, _, _, segments),) = domain_bars(lay)
    assert [s.gid for s in segments] == ["domain-bar:A:C", "domain-bar:A:N", "domain-bar:A:C#2"]
    assert segments[0].tone == segments[2].tone != segments[1].tone


def test_domain_bar_choice_and_legend():
    from foldmap.style import Style

    def bar_and_legend(look, domains=T4L_DOMAINS):
        fig = make_figure(T4L, domains=domains, look=look)
        texts = {t.get_text() for t in fig.axes[0].texts}
        return bool(_gids(fig, "domain-bar:")), texts

    on, texts = bar_and_legend(Style())
    assert on and {"N", "C"} <= texts and "Chain A" not in texts  # one legend entry per domain
    assert not bar_and_legend(Style(domain_bar="off"))[0]
    one = [T4L_DOMAINS[0]]
    assert not bar_and_legend(Style(), one)[0] and bar_and_legend(Style(domain_bar="on"), one)[0]
    with pytest.raises(ValueError, match="domain_bar"):
        Style(domain_bar="maybe").validate()


def test_domain_bar_keeps_every_text_on_the_page():
    """Cas9's eight domains and 8UTF's two chains: every bar text lies inside the figure, the bar under the content
    and above the legend."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    cases = [
        (DATA / "na" / "4OO8.cif.gz", {"chains": ["A"], "domains": _flag_file(CAS9_FILE)}),
        (DATA / "8UTF.cif", {"focus": "protomer", "domains": _flag_file(CAS9_FILE.parent / "8UTF-measles-F.txt")}),
    ]
    for path, kw in cases:
        fig = make_figure(path, **kw)
        renderer = FigureCanvasAgg(fig).get_renderer()
        page = fig.bbox
        bar = [
            a
            for a in [*fig.axes[0].patches, *fig.axes[0].texts, *fig.axes[0].lines]
            if (a.get_gid() or "").startswith("domain-bar")
        ]
        assert bar
        for a in bar:
            box = a.get_window_extent(renderer)
            assert page.x0 - 0.5 <= box.x0 and box.x1 <= page.x1 + 0.5, a.get_gid()
            assert page.y0 - 0.5 <= box.y0 and box.y1 <= page.y1 + 0.5, a.get_gid()
        legend = [t for t in fig.axes[0].texts if t.get_text() in ("F2", "RuvC-I") and not (t.get_gid() or "")]
        lowest_bar = min(a.get_window_extent(renderer).y0 for a in bar)
        assert all(t.get_window_extent(renderer).y1 <= lowest_bar for t in legend)


def test_neighbours_lending_elements_get_no_bar():
    """1LDM, one subunit with its neighbours' lent elements (A4's in D1): one bar, the subunit's own chain."""
    from foldmap.render import domain_bars
    from helpers import FID

    lay, _, _ = make_layout(FID / "1LDM.cif.gz", focus="protomer", domains="auto")
    assert any(k.startswith("A4:") for _, ids in lay.domains for k in ids) and lay.partial
    assert [chain for chain, *_ in domain_bars(lay)] == ["A1"]


def test_domain_bars_share_one_residue_scale():
    """8UTF's F2 (chain A, 72 residues drawn) and F1 (chain a, 344): one residue scale for every bar, as Zyla 2024
    Fig 3A draws the precursor, so F2's bar is that much shorter than F1's; both start at the content's left edge."""
    from foldmap.render import domain_bars

    kw = {"focus": "protomer", "domains": _flag_file(CAS9_FILE.parent / "8UTF-measles-F.txt")}
    lay, _, _ = make_layout(DATA / "8UTF.cif", **kw)
    spans = {chain: last - first for chain, first, last, _ in domain_bars(lay)}
    assert set(spans) == {"A", "a"} and spans["A"] < spans["a"] / 3
    fig = make_figure(DATA / "8UTF.cif", **kw)
    bars = {}
    for p in _gids(fig, "domain-bar:"):
        chain = p.get_gid().split(":")[1]
        x0, x1 = bars.get(chain, (np.inf, -np.inf))
        bars[chain] = (min(x0, p.get_x()), max(x1, p.get_x() + p.get_width()))
    assert bars["A"][0] == pytest.approx(bars["a"][0])
    width = {c: x1 - x0 for c, (x0, x1) in bars.items()}
    assert width["A"] / width["a"] == pytest.approx(spans["A"] / spans["a"], rel=1e-6)


_CVD = {  # Machado, Oliveira and Fernandes 2009, severity 1, on linear RGB
    "deuteranopia": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
    "protanopia": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
    "tritanopia": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
}


def _seen(colour: str, matrix) -> str:
    """The colour as a viewer with that colour-vision deficiency sees it."""
    rgb = np.array([int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5)])
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    out = np.clip(np.array(matrix) @ lin, 0, 1)
    srgb = np.where(out <= 0.0031308, 12.92 * out, 1.055 * out ** (1 / 2.4) - 0.055)
    return "#" + "".join(f"{round(v * 255):02x}" for v in srgb)


def test_cas9_domains_get_eight_different_tones():
    """Eight tones, so Cas9's eight domains (Nishimasu 2014) all differ; the ninth domain takes the first again."""
    from foldmap.render import _DOMAIN_TONES, domain_tone

    domains = _flag_file(CAS9_FILE)
    lay, _, _ = make_layout(DATA / "na" / "4OO8.cif.gz", chains=["A"], domains=domains)
    tones = [domain_tone(k) for k, (_, ids) in enumerate(lay.domains) if ids]
    assert len(tones) == 8 and len(set(tones)) == 8
    assert len(_DOMAIN_TONES) == 8 and domain_tone(8) == domain_tone(0)


def test_domain_tones_stay_apart_for_colour_blind_readers_and_from_nucleic_acid():
    """Every two tones differ by CIEDE2000 >= 10 as seen with normal vision and with each dichromacy (Machado 2009),
    and every tone differs as much from the nucleic-acid backbones (Style dna_color, rna_color and the journal
    theme's RNA), the H-bond lines and the grey of other chains."""
    import itertools

    from foldmap.palette import delta_e
    from foldmap.render import _DOMAIN_TONES, MATE_GREY
    from foldmap.style import THEMES, Style

    for name, matrix in [("normal", None), *_CVD.items()]:
        seen = [c if matrix is None else _seen(c, matrix) for c in _DOMAIN_TONES]
        worst = min(delta_e(a, b) for a, b in itertools.combinations(seen, 2))
        assert worst >= 10, (name, worst)
    others = {Style().dna_color, Style().rna_color, THEMES["journal"].rna_color, Style().hbond_color, MATE_GREY}
    assert min(delta_e(a, b) for a in _DOMAIN_TONES for b in others) >= 10


def test_print_and_blueprint_keep_their_palettes_unless_domain_colours_is_set():
    """The print (greys) and blueprint (pale on dark) themes keep their own palettes with domains named; domain
    colouring applies to them only when domain_colours is set. Colour themes colour by domain as before."""
    from foldmap.render import _DOMAIN_TONES, element_colours
    from foldmap.style import resolve_style

    lay, sses, _ = make_layout(T4L, domains=T4L_DOMAINS)
    tones = set(_DOMAIN_TONES)

    def colours(theme, *overrides):
        return set(element_colours(lay, sses, resolve_style(theme, overrides=list(overrides)))[0].values())

    for theme in ("print", "blueprint"):
        assert not colours(theme) & tones, theme
        assert colours(theme, "domain_colours=true") <= tones, theme
    for theme in ("publication", "minimal", "presentation", "cartoon", "goodsell", "journal"):
        assert colours(theme) <= tones, theme


def test_domain_names_read_on_their_panels_and_on_the_page():
    """A panel's name, and a bar name drawn below the bar, keep a WCAG contrast of 4.5 or more on what they sit on
    (the panel's tint; the page), for all eight tones, on a light page and on a dark one (blueprint's)."""
    from foldmap.render import _DOMAIN_TONES, _mix, contrast_ratio, domain_name_colour
    from foldmap.style import THEMES

    for page in ("#ffffff", THEMES["blueprint"].background):
        for tone in _DOMAIN_TONES:
            panel = _mix(tone, page, 0.86)
            assert contrast_ratio(domain_name_colour(tone, panel), panel) >= 4.5, (page, tone)
            assert contrast_ratio(domain_name_colour(tone, page), page) >= 4.5, (page, tone)


def test_drawn_domain_names_have_contrast():
    """As drawn: Cas9's eight panel names on their panels, the bar names in their segments or on the page below,
    on the default page and on blueprint's dark page with domain colours asked for."""
    from matplotlib.colors import to_hex

    from foldmap.render import contrast_ratio
    from foldmap.style import resolve_style

    for look in (resolve_style("publication"), resolve_style("blueprint", overrides=["domain_colours=true"])):
        fig = make_figure(DATA / "na" / "4OO8.cif.gz", chains=["A"], domains=_flag_file(CAS9_FILE), look=look)
        ax = fig.axes[0]
        panels = {p.get_gid().split(":", 1)[1]: to_hex(p.get_facecolor()) for p in _gids(fig, "domain:")}
        segs = {p.get_gid(): (p, to_hex(p.get_facecolor())) for p in _gids(fig, "domain-bar:")}
        checked = 0
        for t in ax.texts:
            gid, ink = t.get_gid() or "", to_hex(t.get_color())
            if gid.startswith("domain-label:"):
                under = panels[gid.split(":", 1)[1]]
            elif gid.startswith("domain-bar-label:"):
                seg, fill = segs[gid.replace("domain-bar-label:", "domain-bar:", 1)]
                inside = seg.get_y() <= t.get_position()[1] <= seg.get_y() + seg.get_height()
                under = fill if inside else look.background
            else:
                continue
            assert contrast_ratio(ink, under) >= 4.5, (look.background, gid, ink, under)
            checked += 1
        assert checked >= 16


def _bar_segments(fig):
    return sorted((p for p in _gids(fig, "domain-bar:")), key=lambda p: (p.get_y(), p.get_x()))


def test_print_draws_domains_in_its_greys_with_hatches():
    """print (domain_colours off): panels and bar in the theme's greys, never a domain tone; neighbouring bar
    segments differ in grey or hatch, and their names are dark."""
    from matplotlib.colors import to_hex

    from foldmap.render import _DOMAIN_TONES, _mix, contrast_ratio, save_svg
    from foldmap.style import resolve_style

    look = resolve_style("print")
    fig = make_figure(DATA / "na" / "4OO8.cif.gz", chains=["A"], domains=_flag_file(CAS9_FILE), look=look)
    svg = save_svg(fig).lower()
    for tone in _DOMAIN_TONES:
        for shade in (tone, _mix(tone, look.background, 0.86), _mix(tone, look.background, 0.55)):
            assert shade not in svg, shade
    segs = _bar_segments(fig)
    assert len(segs) == 9
    for a, b in zip(segs, segs[1:]):
        assert to_hex(a.get_facecolor()) != to_hex(b.get_facecolor()) or a.get_hatch() != b.get_hatch()
    assert len({(to_hex(p.get_facecolor()), p.get_hatch()) for p in segs}) == 8  # one look per domain
    fills = {p.get_gid(): to_hex(p.get_facecolor()) for p in segs}
    for t in _gids(fig, "domain-bar-label:"):
        seg = fills[t.get_gid().replace("domain-bar-label:", "domain-bar:", 1)]
        assert contrast_ratio(to_hex(t.get_color()), seg) >= 4.5 and to_hex(t.get_color()) == "#000000"


def test_blueprint_draws_domains_in_its_palette():
    """blueprint (domain_colours off): panels tinted and bar segments filled from its own palette."""
    from matplotlib.colors import to_hex

    from foldmap.palette import PALETTES
    from foldmap.render import _mix
    from foldmap.style import resolve_style

    look = resolve_style("blueprint")
    fig = make_figure(T4L, domains=T4L_DOMAINS, look=look)
    pale = [c.lower() for c in PALETTES["blueprint"]]
    panels = {p.get_gid(): to_hex(p.get_facecolor()) for p in _gids(fig, "domain:")}
    assert panels == {
        "domain:N": _mix(pale[0], look.background, 0.86),
        "domain:C": _mix(pale[1], look.background, 0.86),
    }
    assert [to_hex(p.get_facecolor()) for p in _bar_segments(fig)] == [pale[1], pale[0], pale[1]]  # C, N, C


def test_tie_letters_and_a_domain_bar_share_the_page():
    """4YBB (ribosome: tertiary pairs drawn as ties with letters, and their legend entry) with a domain on each of two
    proteins: the ties' letters, the tie legend and the domain bar all lie in the page and none overlaps another."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = make_figure(DATA / "na" / "4YBB.cif.gz", domains=[("S4", ["res:AD:1-300"]), ("S8", ["res:AH:1-300"])])
    renderer = FigureCanvasAgg(fig).get_renderer()
    ax = fig.axes[0]
    drawn = [a for a in [*ax.patches, *ax.texts, *ax.lines] if a.get_gid()]
    ties = [a for a in drawn if a.get_gid().startswith(("legend-tie", "na:tie:"))]  # marks, letters, legend
    bar = [a for a in drawn if a.get_gid().startswith("domain-bar")]
    assert any(a.get_gid().startswith("na:tie:") for a in ties) and any("legend-tie" in a.get_gid() for a in ties)
    assert len(_gids(fig, "domain-bar:")) == 2
    page = fig.bbox
    boxes = {id(a): a.get_window_extent(renderer) for a in ties + bar}
    for a in ties + bar:
        b = boxes[id(a)]
        assert page.x0 - 0.5 <= b.x0 and b.x1 <= page.x1 + 0.5 and page.y0 - 0.5 <= b.y0 and b.y1 <= page.y1 + 0.5
    for a in ties:
        for b in bar:
            assert not boxes[id(a)].overlaps(boxes[id(b)]), (a.get_gid(), b.get_gid())


def _name_overlaps(fig, name):
    """What the drawn name of a domain lies on: elements, sheet plates and other domains' panels (gids)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    renderer = FigureCanvasAgg(fig).get_renderer()
    ax = fig.axes[0]
    (text,) = [t for t in ax.texts if t.get_gid() == f"domain-label:{name}"]
    box = text.get_window_extent(renderer)
    hits = []
    for a in ax.patches:
        gid = a.get_gid() or ""
        if gid.startswith(("strand:", "helix:", "sheet-panel:", "eta:")) or (
            gid.startswith("domain:") and gid != f"domain:{name}"
        ):
            if a.get_window_extent(renderer).overlaps(box):
                hits.append(gid)
    return hits


def test_a_domain_name_takes_the_top_corner_clear_of_elements_and_other_panels():
    """8UTF, one subunit with the docs domains: the DII panel's top left lies under a neighbour's grey strand (A′)
    and its plate, so DII's name goes to the panel's top right, which is clear. A name whose top left is clear stays
    there (ZPN, ZPC on 6ZS5)."""
    kw = {"focus": "protomer", "domains": _flag_file(CAS9_FILE.parent / "8UTF-measles-F.txt")}
    fig = make_figure(DATA / "8UTF.cif", **kw)
    assert _name_overlaps(fig, "DII") == []
    lay, _, _ = make_layout(UMOD, domains=DOMAINS)
    assert lay.domain_name_right == set()


# ---------------------------------------------------------------------------------------------- auto, chain by chain


def _auto_domains_alike(path, **kw):
    """lay.domains of `path` with --domains auto, after checking that each domain lies in one chain, holds at
    least DOMAIN_MIN_SSE_RESIDUES residues of elements, and that every copy of the protein (model.copy_classes)
    has as many domains as its reference and each element in the domain its mate is in."""
    from foldmap.features import DOMAIN_MIN_SSE_RESIDUES
    from foldmap.io import load_backbone
    from foldmap.model import copy_classes, homolog_chains

    lay, sses, bb = make_layout(path, domains="auto", **kw)
    by_id = {s.id: s for s in sses}
    owner = {}
    for name, ids in lay.domains:
        assert len({by_id[k].chain for k in ids}) == 1, name
        assert sum(len(by_id[k]) for k in ids) >= DOMAIN_MIN_SSE_RESIDUES, name
        owner.update({k: name.split(" · ")[0] for k in ids})
    classes = copy_classes(sses, homolog_chains(load_backbone(path, "auto")))
    count = {c: len({owner[s.id] for s in sses if s.chain == c and s.id in owner}) for c in {s.chain for s in sses}}
    for members in classes.groups.values():
        assert len({count[c] for c in members}) == 1, count
        for s in (s for s in sses if s.chain in members):
            assert owner.get(s.id) == owner.get(classes.mate[s.id]) or classes.mate[s.id] == s.id, s.id
    return lay


def test_auto_domains_of_a_tetramer_lie_each_in_one_chain_and_copies_split_alike():
    """1LDM: the domains were communities of the whole tetramer, three or four chains each. Now each chain is split
    on its own (its reference copy's split shared by the others): every copy has the NAD-binding Rossmann domain and
    the catalytic domain, named alike."""
    lay = _auto_domains_alike(FID / "1LDM.cif.gz")
    names = [name for name, _ in lay.domains]
    assert len(names) == 8 and {n.split(" · ")[1] for n in names} == {"A1", "A2", "A3", "A4"}


def test_auto_domains_fold_the_small_ones_after_the_ba_row_is_kept_whole():
    """1TIM: merging the barrel's βα row left chain A a 10-residue D2 (two short helices) and chain B none. The
    smallest domains are folded after the merge now, so neither chain is split: a TIM barrel is one domain."""
    lay = _auto_domains_alike(DATA / "1TIM.cif")
    assert lay.domains == []


def test_copies_of_a_domain_share_one_tone():
    """1LDM --domains auto: four copies of the chain, each split into D1 (the NAD-binding domain) and D2 (the
    catalytic one). Every copy of D1 is drawn in one tone and every copy of D2 in another: panels, bar segments,
    elements and one legend entry each."""
    from matplotlib.colors import to_hex

    from foldmap.render import _DOMAIN_TONES, _legend_entries, domain_bars, element_colours
    from foldmap.style import Style

    lay, sses, _ = make_layout(FID / "1LDM.cif.gz", domains="auto")
    stem = {name: name.split(" · ")[0] for name, _ in lay.domains}
    assert sorted(set(stem.values())) == ["D1", "D2"] and len(lay.domains) == 8
    bar = {seg.name: seg.tone for *_, segs in domain_bars(lay, Style()) for seg in segs}
    colours, _ = element_colours(lay, sses, Style())
    element = {name: {colours[k] for k in ids} for name, ids in lay.domains}
    fig = make_figure(FID / "1LDM.cif.gz", domains="auto")
    panel = {a.get_gid()[len("domain:") :]: to_hex(a.get_edgecolor()) for a in _gids(fig, "domain:")}
    for tones in (bar, {name: next(iter(c)) for name, c in element.items()}, panel):
        by_stem = {s: {tones[n] for n in stem if stem[n] == s} for s in ("D1", "D2")}
        assert all(len(t) == 1 for t in by_stem.values()) and by_stem["D1"] != by_stem["D2"], by_stem
    assert all(len(c) == 1 for c in element.values())
    assert {bar[n] for n in stem} == set(_DOMAIN_TONES[:2])
    chains = ["A1", "A2", "A3", "A4"]
    entries = _legend_entries(lay, Style(), chains, dict.fromkeys(chains, "#000000"), [], {})
    assert [e for e in entries if e[0] in _DOMAIN_TONES] == [(_DOMAIN_TONES[0], "D1"), (_DOMAIN_TONES[1], "D2")]
