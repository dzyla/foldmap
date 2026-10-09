"""Build the foldmap website (GitHub Pages) into a folder: the page in site/ plus every figure, rendered by
foldmap itself so the site always shows what the current code draws.

    python scripts/build_site.py --out _site          # full site
    python scripts/build_site.py --out _site --quick  # a small subset, for tests
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from foldmap import __version__  # noqa: E402
from foldmap.cli import make_figure  # noqa: E402
from foldmap.interactive import write_page  # noqa: E402
from foldmap.render import save  # noqa: E402
from foldmap.seqplot import draw_sequence  # noqa: E402
from foldmap.sequence import read_alignment  # noqa: E402
from foldmap.style import THEMES, resolve_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA, SITE, EXAMPLES = ROOT / "tests" / "data", ROOT / "site", ROOT / "docs" / "examples"
REPO = "https://github.com/dzyla/foldmap"

THEME_NOTES = {
    "publication": "The default: bold fills, coiled-ribbon helices with depth shading, orthogonal loops.",
    "journal": "Muted palette, thinner lines and smaller text, sized for a journal column.",
    "print": "Greyscale with pale fills: reads well in black-and-white print.",
    "minimal": "Outlines only, no shading: the cleanest possible map.",
    "presentation": "Glossy helices, curved loops and large text for slides.",
    "cartoon": "Bright palette, glossy helices and curved loops.",
    "richardson": "Coloured by element type, after Jane Richardson's ribbon drawings.",
    "rainbow": "Each element coloured by its place in the chain, N (blue) to C (red).",
    "shaded": "Each chain in its own hue, elements darkening from N to C.",
    "trace": "N→C shading with coloured loops and direction arrows: easy to follow the chain.",
    "flexibility": "Mean B-factor per element: rigid (blue) to flexible (red).",
    "hydropathy": "Kyte–Doolittle hydropathy: hydrophilic to hydrophobic.",
    "goodsell": "Pale, illustrative colours in the spirit of David Goodsell.",
    "alphafold": "AlphaFold pLDDT confidence bands, loops coloured residue by residue.",
    "conservation": "ConSurf conservation grades from a multiple sequence alignment.",
    "blueprint": "A dark navy page with pale outlined elements.",
}

STRUCTURES = [  # theme playground
    ("1UBQ", "Ubiquitin", "β-grasp fold, one chain"),
    ("5NKT", "FimA pilin", "Ig-like donor-strand fold"),
    ("1LMB", "λ repressor on DNA", "helix-turn-helix dimer bound to its operator"),
    ("1TIM", "Triosephosphate isomerase", "(βα)₈ TIM barrel dimer"),
]
MSA = EXAMPLES / "ubiquitin_family.fasta"

# modes: (key, title, text, structure file, [(caption, theme, overrides, options, command)])
MODES = [
    (
        "symmetry",
        "Symmetric assemblies",
        "Cyclic, dihedral, cubic and helical symmetry is detected automatically; every protomer is drawn alike, side "
        "by side. Highlight the asymmetric unit or one protomer and the copies turn grey.",
        "8UTF.cif",
        [
            (
                "Measles F trimer, one protomer highlighted",
                "publication",
                ["highlight=asu", "loop_color=chain"],
                {},
                "foldmap plot 8UTF --set highlight=asu --set loop_color=chain",
            )
        ],
    ),
    (
        "filaments",
        "Filaments and large assemblies",
        "Filaments, cages and big complexes are drawn one subunit at a time, together with the neighbours' strands "
        "and helices that complete its fold (grey, primed). --focus none draws everything.",
        "6Y7S.cif",
        [
            ("Type 1 pilus rod: one FimA subunit and its neighbours", "publication", [], {}, "foldmap plot 6Y7S"),
            (
                "The same rod, every subunit drawn",
                "publication",
                [],
                {"focus": "none"},
                "foldmap plot 6Y7S --focus none",
            ),
        ],
    ),
    (
        "membrane",
        "Membrane proteins",
        "The lipid bilayer is found from the structure (or from OPM dummy atoms) and drawn as a band, the page turned "
        "so the membrane lies flat; in and out follow the positive-inside rule.",
        "1C3W.cif",
        [
            (
                "Bacteriorhodopsin trimer with its retinal",
                "publication",
                [],
                {"membrane": "auto"},
                "foldmap plot 1C3W --membrane auto",
            )
        ],
    ),
    (
        "dna",
        "Protein–DNA complexes",
        "Duplexes are drawn as double helices with base letters and contact beads, unrolled around the DNA axis so "
        "binding helices sit under the bases they read.",
        "1ZAA.cif",
        [("Zif268 zinc fingers on DNA (blueprint theme)", "blueprint", [], {}, "foldmap plot 1ZAA --theme blueprint")],
    ),
    (
        "alphafold",
        "AlphaFold models",
        "Give a UniProt accession and foldmap downloads the AlphaFold model; the alphafold theme colours elements, "
        "loops and tails by pLDDT confidence.",
        "AF-P04637.cif",
        [
            (
                "p53 (AF-P04637): ordered core, disordered tails",
                "alphafold",
                [],
                {},
                "foldmap plot P04637 --theme alphafold",
            )
        ],
    ),
    (
        "conservation",
        "Conservation from an alignment",
        "Colour the structure by how conserved each position is in a multiple sequence alignment (FASTA, Clustal or "
        "Stockholm), in ConSurf's nine grades.",
        "1UBQ.cif",
        [
            (
                "Ubiquitin coloured by the ubiquitin-like family",
                "conservation",
                [],
                {"msa": MSA},
                "foldmap plot 1UBQ --theme conservation --msa family.fasta",
            )
        ],
    ),
    (
        "domains",
        "Domains and edits",
        "Name domains (or find them automatically) and they are drawn on their own panels; swap, move or rename any "
        "element, and save it all to a layout file that re-creates the figure.",
        "2LZM.cif",
        [
            (
                "T4 lysozyme with its two lobes",
                "journal",
                [],
                {
                    "assembly": "asu",
                    "domains": [("N-lobe", ["res:A:13-59"]), ("C-lobe", ["res:A:1-12", "res:A:60-164"])],
                },
                "foldmap plot 2LZM --domain N-lobe=res:A:13-59 --domain C-lobe=res:A:1-12,res:A:60-164",
            )
        ],
    ),
    (
        "ligands",
        "Ligands, metals, glycans, disulfides",
        "Bound ligands and metal ions are tethered to the residues that hold them; glycans are drawn with SNFG "
        "symbols; disulfide partners are pulled close in the layout.",
        "2HHB.cif",
        [("Haemoglobin and its four hemes", "publication", [], {"assembly": "asu"}, "foldmap plot 2HHB")],
    ),
]
MODES.append(
    (
        "uniprot",
        "UniProt annotation",
        "Each chain is matched to its UniProt entry (SIFTS, or the AlphaFold accession) and its features mapped onto "
        "the structure: domains become panels, active and binding sites are marked. foldmap uniprot lists every "
        "structure of the protein and its AlphaFold model.",
        "2HHB.cif",
        [
            (
                "Haemoglobin: UniProt Globin domains and heme-binding histidines",
                "publication",
                [],
                {"assembly": "asu", "uniprot": "auto", "domains": "uniprot"},
                "foldmap plot 2HHB --uniprot auto --domains uniprot",
            )
        ],
    )
)

SEQUENCES = [  # (key, caption, alignment or None, columns, command)
    ("sequence", "Ubiquitin: sequence with secondary structure", None, 40, "foldmap sequence 1UBQ --columns 40"),
    (
        "alignment",
        "The ubiquitin-like family with 1UBQ's structure on top",
        MSA,
        41,
        "foldmap sequence 1UBQ --msa family.fasta",
    ),
]
EXPLORERS = [("1LMB", "λ repressor on DNA"), ("8UTF", "Measles F trimer"), ("5NKT", "FimA pilin")]


def _svg(fig, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    save(fig, path)
    return path.as_posix()


def _offline_uniprot(folder: Path) -> None:
    """Point UniProt lookups at the saved test fixtures, so the site builds without reaching UniProt or PDBe."""
    import foldmap.uniprot as uniprot

    cache = folder / "uniprot"
    cache.mkdir(parents=True, exist_ok=True)
    for f in (DATA / "uniprot").iterdir():
        shutil.copy(f, cache / f.name)
    uniprot.CACHE = folder


def build(out: Path, quick: bool = False) -> dict:
    out = Path(out)
    import tempfile

    _offline_uniprot(Path(tempfile.mkdtemp(prefix="foldmap-site-")))
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(SITE, out)
    assets = out / "assets"
    themes = list(THEMES)[:2] if quick else list(THEMES)
    structures = STRUCTURES[:1] if quick else STRUCTURES
    modes = MODES[:1] if quick else MODES
    explorers = EXPLORERS[:1] if quick else EXPLORERS

    data: dict = {
        "version": __version__,
        "repo": REPO,
        "themes": [],
        "structures": [],
        "modes": [],
        "sequences": [],
        "explorers": [],
    }
    for name in themes:
        data["themes"].append({"name": name, "note": THEME_NOTES.get(name, ""), "dark": name == "blueprint"})
    for pdb, title, note in structures:
        figs = {}
        for name in themes:
            extra = {"msa": MSA} if name == "conservation" else {}
            if name == "conservation" and pdb != "1UBQ":
                continue  # the example alignment is of ubiquitin-like proteins
            fig = make_figure(DATA / f"{pdb}.cif", title=None, look=THEMES[name], assembly="asu", **extra)
            _svg(fig, assets / "themes" / pdb / f"{name}.svg")
            figs[name] = f"assets/themes/{pdb}/{name}.svg"
        data["structures"].append({"id": pdb, "title": title, "note": note, "figures": figs})
        print(f"themes: {pdb} ({len(figs)})", flush=True)
    for key, title, text, structure, panels in modes:
        rendered = []
        for k, (caption, theme, overrides, options, command) in enumerate(panels):
            fig = make_figure(DATA / structure, look=resolve_style(theme, None, overrides), **options)
            _svg(fig, assets / "modes" / f"{key}-{k}.svg")
            rendered.append(
                {
                    "caption": caption,
                    "src": f"assets/modes/{key}-{k}.svg",
                    "command": command,
                    "dark": theme == "blueprint",
                }
            )
        data["modes"].append({"key": key, "title": title, "text": text, "panels": rendered})
        print(f"mode: {key}", flush=True)
    if not quick:
        for key, caption, msa, columns, command in SEQUENCES:
            fig = draw_sequence(DATA / "1UBQ.cif", alignment=read_alignment(msa) if msa else None, columns=columns)
            _svg(fig, assets / "sequence" / f"{key}.svg")
            data["sequences"].append({"caption": caption, "src": f"assets/sequence/{key}.svg", "command": command})
    for pdb, title in explorers:
        page = assets / "explorers" / f"{pdb}.html"
        page.parent.mkdir(parents=True, exist_ok=True)
        write_page(DATA / f"{pdb}.cif", page, title=f"{title} ({pdb})")
        data["explorers"].append({"id": pdb, "title": title, "src": f"assets/explorers/{pdb}.html"})
        print(f"explorer: {pdb}", flush=True)

    index = out / "index.html"
    blob = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    index.write_text(index.read_text().replace("{{SITE_DATA}}", blob).replace("{{VERSION}}", __version__))
    (out / ".nojekyll").write_text("")
    return data


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="_site")
    ap.add_argument("--quick", action="store_true", help="a small subset (tests)")
    args = ap.parse_args()
    data = build(Path(args.out), args.quick)
    n = sum(len(s["figures"]) for s in data["structures"]) + sum(len(m["panels"]) for m in data["modes"])
    print(f"wrote {args.out}: {n} figures, {len(data['explorers'])} explorers")


if __name__ == "__main__":
    main()
