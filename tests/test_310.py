"""Short 3-10 helices: their own element kind (G), drawn as small boxes, labelled η1, η2..."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import make_figure, make_layout
from foldmap.style import Style

DATA = Path(__file__).parent / "data"


def _sses(name, on=True):
    from foldmap.dssp import assign_dssp
    from foldmap.io import load_backbone
    from foldmap.ss import build_sses

    bb = load_backbone(DATA / f"{name}.cif", "asu")
    d = assign_dssp(bb)
    return bb, d.ss, build_sses(bb, d.ss, short_helices=on)


@pytest.mark.parametrize("name, n", [("1UBQ", 2), ("5NKT", 1), ("1LMB", 2), ("2LZM", 0)])
def test_isolated_310_runs_become_g_elements(name, n):
    bb, ss, sses = _sses(name)
    g = [s for s in sses if s.kind == "G"]
    assert len(g) == n
    for s in g:
        assert set(ss[s.start : s.end + 1]) == {"G"} and len(s) >= 3
    spans = sorted((s.start, s.end) for s in sses)
    assert all(a[1] < b[0] for a, b in zip(spans, spans[1:]))  # never overlapping another element


def test_off_by_default_at_the_unit_level():
    _, _, sses = _sses("1UBQ", on=False)
    assert not any(s.kind == "G" for s in sses)


def test_310_layout_size_and_labels():
    from foldmap.layout import HELIX_W, SCALE, TEN_RISE

    lay, sses, _ = make_layout(DATA / "1UBQ.cif")
    g = [s for s in sses if s.kind == "G"]
    assert len(g) == 2
    for k, s in enumerate(g, 1):
        p = lay.placed[s.id]
        assert p.width < HELIX_W and p.length == pytest.approx(max(p.width, (len(s) - 1) * TEN_RISE * SCALE))
        assert p.label == f"η{k}"
    assert sum(s.kind == "H" for s in sses) == 1 and lay.placed[next(s.id for s in sses if s.kind == "H")].label == "α1"


def _gids(fig, prefix):
    return [a for a in [*fig.axes[0].patches, *fig.axes[0].texts] if (a.get_gid() or "").startswith(prefix)]


def test_310_drawn_as_boxes_and_switchable():
    from matplotlib.patches import FancyBboxPatch

    fig = make_figure(DATA / "1UBQ.cif")
    boxes = _gids(fig, "eta:")
    assert len(boxes) == 2 and all(isinstance(b, FancyBboxPatch) for b in boxes)
    assert {t.get_text() for t in _gids(fig, "label:")} >= {"η1", "η2"}
    off = make_figure(DATA / "1UBQ.cif", look=Style(helices_310=False))
    assert not _gids(off, "eta:")


def test_loops_reach_310_elements(ubq):
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(ubq)
    loops = route_loops(lay, sses, bb)
    ids = {s.id for s in sses if s.kind == "G"}
    touching = [l for l in loops if l.a_id in ids or l.b_id in ids]
    assert len(touching) == 2 * len(ids)
    for l in touching:
        if l.b_id in ids:
            assert np.allclose(l.points[-1], lay.placed[l.b_id].n_port)


def test_richardson_gives_310_its_own_colour():
    from matplotlib.colors import to_hex

    fig = make_figure(DATA / "1UBQ.cif", look=Style(color_by="sstype"))
    eta = {to_hex(b.get_facecolor()) for b in _gids(fig, "eta:")}
    helix = {to_hex(b.get_facecolor()) for b in _gids(fig, "helix:")}
    assert len(eta) == 1 and eta != helix


def test_summary_counts_310(capsys):
    from foldmap.cli import main

    assert main(["summary", str(DATA / "1UBQ.cif")]) == 0
    assert "3₁₀ helices: 2" in capsys.readouterr().out


def test_interactive_page_lists_310_elements():
    import json
    import re

    from foldmap.interactive import build_page

    page = build_page(DATA / "1UBQ.cif")
    data = json.loads(re.search(r'id="topo-data">(.*?)</script>', page, re.S).group(1))
    g = [e for e in data["elements"] if e["kind"] == "G"]
    assert len(g) == 2 and all(f'id="eta:{e["id"]}"' in page for e in g)
