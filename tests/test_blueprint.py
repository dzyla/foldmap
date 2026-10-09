"""Dark backgrounds: the blueprint theme, and background/ink as style keys, with everything legible."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex, to_rgb

from foldmap.cli import make_figure
from foldmap.render import save
from foldmap.seqplot import draw_sequence
from foldmap.style import THEMES, Style

DATA = Path(__file__).parent / "data"
BLUE = THEMES["blueprint"]


def _lum(c):
    r, g, b = (x / 12.92 if x <= 0.03928 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(c))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a, b):
    la, lb = sorted((_lum(a), _lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def test_theme_is_dark():
    assert _lum(BLUE.background) < 0.05 and contrast(BLUE.ink, BLUE.background) > 7


@pytest.mark.parametrize("name", ["1UBQ", "1LMB", "2HHB", "5NKT"])
def test_everything_legible_on_the_dark_page(name, tmp_path):
    fig = make_figure(DATA / f"{name}.cif", look=BLUE, title="t", assembly="asu")
    assert to_hex(fig.get_facecolor()) == to_hex(BLUE.background)
    ax = fig.axes[0]
    for t in ax.texts:
        if t.get_text().strip() and not (t.get_gid() or "").startswith(("ligand-label", "label:")):
            assert contrast(t.get_color(), BLUE.background) >= 3, (t.get_text(), t.get_color())
    for p in ax.patches:
        gid = p.get_gid() or ""
        if gid.startswith(("loop:", "stub:")):
            assert contrast(to_hex(p.get_edgecolor()), BLUE.background) >= 3, gid
        if gid.startswith("loopcase:"):
            assert to_hex(p.get_edgecolor()) == to_hex(BLUE.background)
    png = save(fig, tmp_path / "x.png", dpi=50)
    from matplotlib.image import imread

    assert np.allclose(imread(png)[1, 1, :3], to_rgb(BLUE.background), atol=0.02)


def test_background_and_ink_are_style_keys():
    look = Style(background="#202020", ink="#f0f0f0").validate()
    fig = make_figure(DATA / "1UBQ.cif", look=look)
    assert to_hex(fig.get_facecolor()) == "#202020"
    with pytest.raises(ValueError, match="background"):
        Style(background="nope").validate()


def test_sequence_view_follows_the_theme():
    fig = draw_sequence(DATA / "1UBQ.cif", look=BLUE)
    assert to_hex(fig.get_facecolor()) == to_hex(BLUE.background)
    letters = [t for t in fig.axes[0].texts if (t.get_gid() or "").startswith("res:")]
    assert all(contrast(t.get_color(), BLUE.background) >= 3 for t in letters)
