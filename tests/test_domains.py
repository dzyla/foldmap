"""Domain panels: named groups of elements drawn on a soft panel with their name."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main, make_figure, make_layout

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
