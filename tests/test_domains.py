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
    assert main(["plot", str(UMOD), "-o", str(out), "--domain", "ZPN=res:A:331-440", "--domain",
                 "ZPC=res:D:447-582", "--save-layout", str(lf)]) == 0
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
