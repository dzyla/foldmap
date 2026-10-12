"""How the app looks: one theme source, contrast, a small stylesheet, the icon, the hosting files."""

import re
import tomllib
from pathlib import Path

import pytest

from foldmap.app import APP_CONFIG, APP_THEME, MUTED, STYLESHEET, config_toml, run_args
from foldmap.render import contrast_ratio

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / ".streamlit" / "config.toml"
in_repo = pytest.mark.skipif(not CONFIG.exists(), reason="the Streamlit config is not part of the sdist")


@in_repo
def test_config_file_is_generated_from_app_config():
    assert CONFIG.read_text(encoding="utf-8") == config_toml()
    parsed = tomllib.loads(CONFIG.read_text(encoding="utf-8"))
    assert parsed["client"]["toolbarMode"] == "viewer" and parsed["theme"]["light"]["primaryColor"] == "#0072B2"
    assert parsed["theme"]["dark"]["backgroundColor"] == "#0F1419" and parsed["server"]["maxUploadSize"] == 50


def test_foldmap_app_passes_the_same_values():
    args = run_args("app.py")
    assert args[1:4] == ["-m", "streamlit", "run"] and "--client.toolbarMode=viewer" in args
    assert "--theme.light.primaryColor=#0072B2" in args and "--theme.dark.primaryColor=#56B4E9" in args
    assert len([a for a in args if a.startswith("--")]) == sum(len(v) for v in APP_CONFIG.values())


def test_hosted_visitors_see_no_traceback_and_uploads_are_capped():
    assert APP_CONFIG["client"]["showErrorDetails"] == "none"  # defence in depth: the app catches its errors
    assert APP_CONFIG["server"]["maxUploadSize"] == 50 and APP_CONFIG["browser"]["gatherUsageStats"] is False
    assert "--client.showErrorDetails=none" in run_args("app.py")
    assert "--browser.gatherUsageStats=false" in run_args("app.py")


@pytest.mark.parametrize("mode", ["light", "dark"])
def test_theme_colours_meet_contrast(mode):
    t = APP_THEME[f"theme.{mode}"]
    for bg in (t["backgroundColor"], t["secondaryBackgroundColor"]):
        assert contrast_ratio(t["textColor"], bg) >= 7
        assert contrast_ratio(t["primaryColor"], bg) >= 4.5  # tightest: #0072B2 on #F5F7FA, 4.83
        assert contrast_ratio(MUTED[mode], bg) >= 4.5
        assert contrast_ratio(t["borderColor"], bg) >= 3  # Q-A18


def test_primary_button_labels_meet_contrast_in_both_themes():
    from foldmap.app import BUTTON_LABEL, page_css

    for mode in ("light", "dark"):
        assert contrast_ratio(BUTTON_LABEL[mode], APP_THEME[f"theme.{mode}"]["primaryColor"]) >= 4.5
    assert contrast_ratio("#FFFFFF", APP_THEME["theme.dark"]["primaryColor"]) < 4.5  # why dark needs its own
    assert BUTTON_LABEL["dark"] in page_css("dark") and BUTTON_LABEL["dark"] not in page_css("light")
    assert MUTED["dark"] in page_css("dark") and MUTED["dark"] not in page_css("light")
    assert page_css(None) == page_css("light") and page_css("light").startswith(STYLESHEET)


def test_stylesheet_is_small_and_keeps_focus():
    lines = [l for l in STYLESHEET.splitlines() if l.strip()]
    assert len(lines) <= 60
    assert "outline: none" not in STYLESHEET and "outline:0" not in STYLESHEET.replace(" ", "")
    selectors = re.findall(r"([^{}]+)\{", STYLESHEET)
    for sel in selectors:
        for part in sel.split(","):
            part = part.strip()
            if part.startswith("@media"):
                continue
            assert re.match(r'(\.st-key-[\w-]+|\.fm-[\w-]+|\[data-testid="\w+"\])', part), part


def test_icon_assets_ship_with_the_package():
    from foldmap.app import ASSETS

    png = (ASSETS / "icon.png").read_bytes()
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) < 8_000
    assert (ASSETS / "icon.svg").read_text(encoding="utf-8").lstrip().startswith(("<?xml", "<svg"))


@in_repo
def test_hosting_installs_the_figure_font():
    assert "fonts-liberation" in (ROOT / "packages.txt").read_text(encoding="utf-8").split()
