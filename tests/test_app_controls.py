"""The sidebar: a short essentials list, then sections for everything else."""

import pytest

pytest.importorskip("streamlit")

from app_helpers import DATA, count_calls, opened, svg_of, widget_keys  # noqa: E402


def test_first_view_shows_at_most_twelve_controls():
    at = opened(DATA / "1UBQ.cif")
    side, main = widget_keys(at.sidebar), widget_keys(at.main)
    assert side == [
        "theme", "theme_browse", "color_by", "palette", "labels", "residue_numbers", "legend",
        "font_scale", "title", "reset_theme",
    ]  # fmt: skip
    assert main == ["pane", "new_structure"] and len(side) + len(main) <= 12  # the view switch is `pane` (R1)


def test_palette_gives_way_where_it_colours_nothing():
    at = opened(DATA / "1UBQ.cif")
    at.selectbox(key="color_by").select("sequence").run()
    keys = widget_keys(at.sidebar)
    assert "sequence_map" in keys and "palette" not in keys


def test_theme_dialog_applies_a_theme_and_draws_nothing_new(monkeypatch):
    import foldmap.spec as spec
    from foldmap.app import PREVIEW_THEMES

    at = opened(DATA / "1LMB.cif")
    calls = count_calls(monkeypatch, spec, "draw_spec")
    at.button(key="theme_browse").click().run()
    assert len(at.image) == len(PREVIEW_THEMES) and calls == []  # thumbnails are static files
    assert any(c.value.startswith("Shown on ubiquitin, the reference structure") for c in at.caption)
    at.button(key="pick_trace").click().run()
    assert not at.exception and at.selectbox(key="theme").value == "trace"
    assert at.session_state["loop_arrows"] is True and not at.image and "<svg" in svg_of(at)


def test_reset_restores_the_theme():
    from foldmap.style import THEMES

    at = opened(DATA / "1UBQ.cif")
    at.selectbox(key="palette").select("tol-bright").run()
    at.toggle(key="legend").set_value(False).run()
    at.button(key="reset_theme").click().run()
    assert at.session_state["palette"] == THEMES["publication"].palette and at.session_state["legend"] is True


def test_controls_change_the_figure():
    """Was in test_app.py (its residue-numbers half; the mode half is in Task 9)."""
    at = opened(DATA / "1UBQ.cif")
    assert "resnum:" not in svg_of(at)
    at.toggle(key="residue_numbers").set_value(True).run()
    assert "resnum:" in svg_of(at)
