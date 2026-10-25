# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Liberation Sans forced, for the cross-platform determinism run (determinism.py, same directory).

render._FONTS lists Arial first, so macOS and Windows draw with Arial even where Liberation Sans 2.1.5 is installed
(CI installs it on all three). foldmap writes SVG text as text (svg.fonttype "none"), not as glyph paths, but the
layout measures every label in the resolved font and the SVG names it: another font gives other positions and other
SVG hashes. The determinism run must therefore force Liberation Sans before any figure is drawn, on every platform,
and assert what the text resolved to; font_report() records the font file's revision too.

render binds the list at import: _MEASURE = FontProperties(family=_FONTS), _RC["font.sans-serif"] is _FONTS, and
set_fontfamily(_FONTS) is called with it. So the list is mutated IN PLACE (rebinding the name would leave those
holding the old order), and _MEASURE's family is set again in case FontProperties copied the list.

Use, in the process that draws (before the first figure):

    import fonts
    fonts.force_liberation_sans()          # raises if Liberation Sans is not what text resolves to
    ... draw with foldmap ...
    fonts.assert_figure_font(fig)          # every text of the drawn figure resolved to Liberation Sans

and record fonts.font_report() with each platform's hashes.
"""

from __future__ import annotations

FONT = "Liberation Sans"  # figcheck.INK_FONT: the font CI, the hosted app and the ink gates use


def force_liberation_sans() -> str:
    """Put Liberation Sans first in render._FONTS (in place) and check it resolves; returns the font file path."""
    from matplotlib.font_manager import findfont
    from matplotlib.ft2font import FT2Font

    from foldmap import figcheck, render

    assert figcheck.INK_FONT == FONT, figcheck.INK_FONT
    order = [FONT] + [f for f in render._FONTS if f != FONT]
    render._FONTS[:] = order  # in place: _RC and every caller holding the list see the new order
    render._MEASURE.set_family(render._FONTS)
    assert render._RC["font.sans-serif"] is render._FONTS or render._RC["font.sans-serif"] == order
    render._RC["font.sans-serif"] = render._FONTS
    path = findfont(render._MEASURE, fallback_to_default=False)
    family = FT2Font(path).family_name
    if family != FONT:
        raise RuntimeError(f"text resolves to {family} ({path}), not {FONT}: install fonts-liberation 2.1.5")
    return path


def assert_figure_font(fig) -> None:
    """Every text of a drawn matplotlib figure resolved to Liberation Sans (figcheck.ink_font_reason is empty)."""
    from foldmap import figcheck

    why = figcheck.ink_font_reason(fig)
    if why:
        raise RuntimeError(why)


def font_report() -> dict:
    """What to store beside each platform's layout signatures and SVG hashes."""
    import platform

    import matplotlib
    from matplotlib.font_manager import findfont
    from matplotlib.ft2font import FT2Font

    from foldmap import render

    path = findfont(render._MEASURE, fallback_to_default=False)
    font = FT2Font(path)
    head = font.get_sfnt_table("head") or {}
    major, minor = head.get("fontRevision", (0, 0))
    version = font.get_sfnt().get((3, 1, 0x409, 5)) or font.get_sfnt().get((1, 0, 0, 5)) or b""  # name ID 5
    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "font_order": list(render._FONTS),
        "font_file": path,
        "font_family": font.family_name,
        "font_style": font.style_name,
        "font_revision": f"{major + minor / 65536:.3f}",  # head.fontRevision (16.16)
        "font_version": version.decode("utf-16-be" if b"\x00" in version else "latin-1").strip(),
    }
