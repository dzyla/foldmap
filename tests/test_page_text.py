"""random-100 C4: a figure's text is never drawn under layout.MIN_TEXT_PT (6 pt): a page too wide for MAX_WIDTH_IN at
that size grows instead, up to MAX_PAGE_IN, past which the whole figure is scaled to fit (and the notes say so).
The text's em is then the same page units on every page, so on a wide page no name's ink reaches past the box
loops keep out of (on 8UUP, an 84-unit page, 4 pt text took 0.65 page units an em and ran 0.08 past its box), and the
legend's rows are sized from the text drawn in them (layout.glyph_scale)."""

from pathlib import Path

import pytest

from foldmap import figcheck, layout
from foldmap.cli import make_figure, make_layout
from foldmap.layout import label_boxes
from foldmap.route import route_loops
from foldmap.style import Style

DATA = Path(__file__).parent / "data"
SIX_PT = 6.0  # the ruling's floor (layout.MIN_TEXT_PT)


def _texts(fig, prefix):
    return [t for t in fig.axes[0].texts if (t.get_gid() or "").startswith(prefix)]


def test_a_wide_page_keeps_its_text_at_six_points_and_grows_instead():
    fig = make_figure(DATA / "8UUP.cif")  # 84 page units wide: 7.2 in would need 4 pt text
    names = _texts(fig, "label:")  # helix names at the figure's text size; strand letters at 0.9 of it
    assert names and max(t.get_fontsize() for t in names) >= SIX_PT - 1e-9
    assert fig.get_size_inches()[0] > layout.MAX_WIDTH_IN
    assert fig.foldmap_page_fit == ""  # under MAX_PAGE_IN: nothing to say


def test_a_small_page_is_drawn_as_before():
    fig = make_figure(DATA / "1UBQ.cif")  # 12 units: IN_PER_UNIT, 7.4 pt text
    assert fig.get_size_inches()[0] < layout.MAX_WIDTH_IN
    assert max(t.get_fontsize() for t in _texts(fig, "label:")) == pytest.approx(0.55 * layout.SHAFT * 0.3 * 72)


def test_every_name_and_chain_end_ink_lies_inside_its_route_keep_out_on_a_wide_page():
    path = DATA / "8UUP.cif"
    fig = make_figure(path)
    inv = fig.axes[0].transData.inverted()
    lay, sses, bb = make_layout(path)
    route_loops(lay, sses, bb)
    boxes = dict(label_boxes(lay, sses, route=True))
    seen = 0
    for t, ink in figcheck._ink_boxes(fig):
        gid = t.get_gid() or ""
        key = gid if gid.startswith("label:") else ":".join(gid.split(":")[1:]) if gid.startswith("terminus:") else ""
        if key not in boxes:
            continue
        (a0, b0), (a1, b1) = inv.transform(ink.get_points())
        x0, y0, x1, y1 = boxes[key]
        assert x0 <= a0 and y0 <= b0 and a1 <= x1 and b1 <= y1, (gid, (a0, b0, a1, b1), boxes[key])
        seen += 1
    assert seen > 40


def test_legend_rows_are_sized_from_their_text(monkeypatch):
    """Text drawn four times its usual size (Style.font_scale): the legend's row step (layout.glyph_scale) grows with
    it, so its rows do not run into each other (twelve chain entries on 1UBQ's small page wrap to several rows; at
    f57887a rows stood 1.0 page unit apart under 1.4-unit text)."""
    from foldmap import render

    def twelve(*a, **k):  # twelve chain entries, so the legend wraps
        return [("#0072b2", f"Chain {c}") for c in "ABCDEFGHIJKL"]

    monkeypatch.setattr(render, "_legend_entries", twelve)
    fig = make_figure(DATA / "1UBQ.cif", look=Style(font_scale=4.0))
    inv = fig.axes[0].transData.inverted()
    inks = [inv.transform(b.get_points())[:, 1] for t, b in figcheck._ink_boxes(fig) if _legend(t)]
    rows: dict[float, tuple[float, float]] = {}  # a row (its entries share a baseline): its ink's low and high y
    for lo, hi in inks:
        key = next((k for k in rows if abs(k - (lo + hi) / 2) < 0.3), (lo + hi) / 2)
        a, b = rows.get(key, (lo, hi))
        rows[key] = (min(a, lo), max(b, hi))
    spans = sorted(rows.values(), reverse=True)
    assert len(spans) >= 2
    for (low, _), (_, high) in zip(spans, spans[1:]):  # each row's ink ends above the next row's
        assert low > high, spans


def _legend(t) -> bool:
    return (t.get_gid() or "").startswith("legend-chain-label:")


def test_past_the_largest_page_the_whole_figure_is_scaled_and_said(monkeypatch):
    from foldmap.layout import page_scale, text_em

    assert layout.MIN_TEXT_PT == SIX_PT
    assert text_em(400.0, 60.0) == pytest.approx(text_em(10.0))  # the text keeps its page units
    assert 400.0 * page_scale(400.0, 60.0) == pytest.approx(layout.MAX_PAGE_IN)
    assert 60.0 * page_scale(10.0, 60.0) <= layout.MAX_PAGE_IN + 1e-9
    monkeypatch.setattr(layout, "MAX_PAGE_IN", 3.0)  # 1UBQ is 4 in wide: past it
    fig = make_figure(DATA / "1UBQ.cif")
    assert max(fig.get_size_inches()) == pytest.approx(3.0)
    assert "held at 3 in" in fig.foldmap_page_fit and "under 6 pt" in fig.foldmap_page_fit
    assert "its text is 5.9 pt, under 6 pt" in layout.page_fit(48.0, 44.0, 5.97)  # not "6.0 pt, under 6" (7Y7C)
    assert layout.page_fit(48.0, 44.0, 6.0) == ""


def past_the_largest_page(monkeypatch):
    from foldmap import render

    monkeypatch.setattr(layout, "MAX_PAGE_IN", 3.0)  # 1UBQ's page (4 x 2.6 in) is past it
    monkeypatch.setattr(render, "MAX_PNG_PIXELS", 200_000)  # and its PNG past this at 300 dpi


def test_a_page_held_at_its_largest_and_a_lower_png_dpi_are_told_on_the_command_line(tmp_path, monkeypatch, capsys):
    from foldmap.cli import main

    out = tmp_path / "u.png"
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(out)]) == 0
    assert "held at" not in capsys.readouterr().err
    past_the_largest_page(monkeypatch)
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(out)]) == 0
    err = capsys.readouterr().err
    assert "How this figure was made:\n  - the page (3 x 2 in) is held at 3 in on its longer side" in err
    assert "its PNG is " in err and "dpi, not 300" in err
    w, h = int.from_bytes(out.read_bytes()[16:20], "big"), int.from_bytes(out.read_bytes()[20:24], "big")
    assert w * h <= 200_000


def test_glyph_room_without_a_layout_follows_the_text_size(monkeypatch):
    """Batch C1 review (finding 5): an item's footprint and a stem's name box hold no layout; while build_layout runs
    they take the glyph scale of its style's text (layout._GLYPH_K), as layout.glyph does for a layout, and 1
    outside it (so at text size 1 nothing changes: bench and snapshots identical)."""
    seen = []
    real = layout._footprint
    monkeypatch.setattr(layout, "_footprint", lambda m, k=None: seen.append(layout._GLYPH_K.get()) or real(m, k))
    make_layout(DATA / "1UBQ.cif", look=Style(font_scale=2.0))
    assert seen and set(seen) == {2.0} and layout._GLYPH_K.get() == 1.0
    seen.clear()
    make_layout(DATA / "1UBQ.cif")
    assert set(seen) == {1.0}


def test_stem_name_boxes_grow_with_the_text():
    """Batch C1 review (finding 5): a stem's name box (its keep-out where its ink is not measured) is sized from the
    text drawn, through layout.glyph."""
    lay, sses, _ = make_layout(DATA / "1LMB.cif")
    small = {n: b for n, b in label_boxes(lay, sses) if n.startswith("stem-label:")}
    lay.font_scale = 2.0
    big = {n: b for n, b in label_boxes(lay, sses) if n.startswith("stem-label:")}
    assert small and all((big[n][2] - big[n][0]) > (small[n][2] - small[n][0]) for n in small)
