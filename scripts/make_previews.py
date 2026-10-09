"""Render the theme thumbnails shown in the app's style picker (src/foldmap/assets/previews/*.png).

They are drawn once, on ubiquitin, and shipped with the package, so picking a style is instant whatever
structure is loaded. Re-run after changing a theme:

    python scripts/make_previews.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from foldmap.app import PREVIEW_THEMES  # noqa: E402
from foldmap.cli import make_figure  # noqa: E402
from foldmap.render import save  # noqa: E402
from foldmap.style import THEMES, resolve_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src" / "foldmap" / "assets" / "previews"
SOURCE = {"alphafold": ROOT / "tests" / "data" / "AF-P04637.cif"}  # pLDDT needs a predicted model


CANVAS = (420, 280)  # every thumbnail the same size, so the picker's rows line up


def _letterbox(path: Path, background: str) -> None:
    """Fit the image inside CANVAS, centred on the theme's page colour."""
    from matplotlib.colors import to_rgb
    from PIL import Image

    img = Image.open(path).convert("RGB")
    img.thumbnail((CANVAS[0] - 16, CANVAS[1] - 16), Image.LANCZOS)
    page = Image.new("RGB", CANVAS, tuple(int(255 * c) for c in to_rgb(background)))
    page.paste(img, ((CANVAS[0] - img.width) // 2, (CANVAS[1] - img.height) // 2))
    page.save(path, optimize=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for theme in PREVIEW_THEMES:
        assert theme in THEMES, theme
        path = SOURCE.get(theme, ROOT / "tests" / "data" / "1UBQ.cif")
        fig = make_figure(path, look=resolve_style(theme, None, ["legend=false"]))
        out = save(fig, OUT / f"{theme}.png", dpi=110)
        _letterbox(out, THEMES[theme].background)
        print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
