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

EXAMPLES = ROOT / "docs" / "examples"

# (output name, structure, title, theme, style overrides, layout options)
FIGURES = [
    ("lambda-repressor-dna", "1LMB.cif", "λ repressor on DNA (1LMB)", "publication", [], {"assembly": "asu"}),
    (
        "t4-lysozyme-domains",
        "2LZM.cif",
        "T4 lysozyme (2LZM) · domains",
        "journal",
        [],
        {"assembly": "asu", "domains": [("N-lobe", ["res:A:13-59"]), ("C-lobe", ["res:A:1-12", "res:A:60-164"])]},
    ),
    ("fima-richardson", "5NKT.cif", "FimA (5NKT) · richardson", "richardson", [], {"assembly": "asu"}),
    ("p53-alphafold", "AF-P04637.cif", "p53 AlphaFold model (AF-P04637) · pLDDT", "alphafold", [], {}),
    (
        "measles-f-trimer",
        "8UTF.cif",
        "Measles F trimer (8UTF) · C3, one protomer highlighted",
        "publication",
        ["highlight=asu", "loop_color=chain"],
        {},
    ),
    ("tim-barrel-rainbow", "1TIM.cif", "TIM barrel (1TIM) · N→C", "rainbow", [], {"assembly": "asu"}),
    (
        "bacteriorhodopsin-membrane",
        "1C3W.cif",
        "Bacteriorhodopsin (1C3W) · membrane, retinal",
        "publication",
        [],
        {"membrane": "auto"},
    ),
    ("zinc-fingers-blueprint", "1ZAA.cif", "Zif268 zinc fingers on DNA (1ZAA) · blueprint", "blueprint", [], {}),
    ("haemoglobin-hemes", "2HHB.cif", "Haemoglobin (2HHB) · hemes", "journal", [], {"assembly": "asu"}),
    (
        "ubiquitin-conservation",
        "1UBQ.cif",
        "Ubiquitin (1UBQ) · conservation in the ubiquitin-like family",
        "conservation",
        [],
        {"msa": EXAMPLES / "ubiquitin_family.fasta"},
    ),
]

# (output name, structure, title, alignment or None, columns, theme)
SEQUENCES = [
    ("ubiquitin-sequence", "1UBQ.cif", "Ubiquitin (1UBQ)", None, 40, "publication"),
    (
        "ubiquitin-family-alignment",
        "1UBQ.cif",
        "Ubiquitin-like family · structure of 1UBQ on top",
        EXAMPLES / "ubiquitin_family.fasta",
        41,
        "publication",
    ),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, structure, title, theme, overrides, options in FIGURES:
        fig = make_figure(DATA / structure, title=title, look=resolve_style(theme, None, overrides), **options)
        path = OUT / f"{name}.png"
        save(fig, path, dpi=110)
        print(f"wrote {path.relative_to(ROOT)} ({path.stat().st_size // 1024} KB)")
    from foldmap.seqplot import draw_sequence
    from foldmap.sequence import read_alignment
    from foldmap.style import THEMES

    for name, structure, title, msa, columns, theme in SEQUENCES:
        fig = draw_sequence(
            DATA / structure,
            alignment=read_alignment(msa) if msa else None,
            columns=columns,
            look=THEMES[theme],
            title=title,
        )
        path = OUT / f"{name}.png"
        save(fig, path, dpi=130)
        print(f"wrote {path.relative_to(ROOT)} ({path.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
