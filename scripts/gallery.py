"""Render the README gallery into docs/images from the structures in tests/data.

python scripts/gallery.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from foldmap.cli import make_figure  # noqa: E402
from foldmap.render import save  # noqa: E402
from foldmap.style import resolve_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / "tests" / "data", ROOT / "docs" / "images"

# (output name, structure, title, theme, style overrides, layout options)
FIGURES = [
    ("ubiquitin", "1UBQ.cif", "Ubiquitin (1UBQ)", "publication", ["residue_numbers=true"], {"assembly": "asu"}),
    ("lambda-repressor-dna", "1LMB.cif", "λ repressor on DNA (1LMB)", "publication", [], {"assembly": "asu"}),
    ("fima-richardson", "5NKT.cif", "FimA (5NKT) · richardson", "richardson", [], {"assembly": "asu"}),
    (
        "t4-lysozyme-domains",
        "2LZM.cif",
        "T4 lysozyme (2LZM) · domains",
        "journal",
        [],
        {"assembly": "asu", "domains": [("N-lobe", ["res:A:13-59"]), ("C-lobe", ["res:A:1-12", "res:A:60-164"])]},
    ),
    (
        "measles-f-trimer",
        "8UTF.cif",
        "Measles F trimer (8UTF) · C3, one protomer highlighted",
        "publication",
        ["highlight=asu", "loop_color=chain"],
        {},
    ),
    ("tim-barrel-rainbow", "1TIM.cif", "TIM barrel (1TIM) · N→C", "rainbow", [], {"assembly": "asu"}),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, structure, title, theme, overrides, options in FIGURES:
        fig = make_figure(DATA / structure, title=title, look=resolve_style(theme, None, overrides), **options)
        path = OUT / f"{name}.png"
        save(fig, path, dpi=110)
        print(f"wrote {path.relative_to(ROOT)} ({path.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
