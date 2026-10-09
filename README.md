# foldmap

[![tests](https://github.com/dzyla/foldmap/actions/workflows/tests.yml/badge.svg)](https://github.com/dzyla/foldmap/actions/workflows/tests.yml)
[![license: AGPL-3.0](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
![python](https://img.shields.io/badge/python-3.11%2B-blue.svg)

**Publication-grade protein topology diagrams, straight from a PDB/mmCIF file, a PDB ID or an AlphaFold model.**

**Website, live theme gallery and interactive examples: <https://dzyla.github.io/foldmap/>**

```bash
pip install git+https://github.com/dzyla/foldmap.git
foldmap plot 1LMB -o lambda.svg
```

foldmap reads the structure, assigns secondary structure, finds β-sheets, helix bundles and symmetry, and lays the
elements out as a clean 2D map: coiled-ribbon helices, strand arrows and 3₁₀ helices on one consistent scale,
with loops routed around everything else. It also draws the sequence with the secondary structure on top,
alone or over a multiple sequence alignment.

<p align="center">
<img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/lambda-repressor-dna.png" width="560"
alt="Topology of the lambda repressor dimer bound to its operator DNA">
</p>

## Gallery

<table>
<tr>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/t4-lysozyme-domains.png" alt="T4 lysozyme with N- and C-lobe domain panels"><br>
<sub><b>Domains.</b> T4 lysozyme (2LZM), journal theme: named domain panels, each kept together.</sub></td>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/measles-f-trimer.png" alt="Measles fusion protein trimer with one protomer highlighted"><br>
<sub><b>Symmetry.</b> Measles F trimer (8UTF): C3 found, protomers drawn alike, one highlighted, disulfides and glycans marked.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/bacteriorhodopsin-membrane.png" alt="Bacteriorhodopsin trimer in a membrane band"><br>
<sub><b>Membranes.</b> Bacteriorhodopsin (1C3W): the bilayer found from the structure, sides by the positive-inside rule, retinal bound.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/zinc-fingers-blueprint.png" alt="Zif268 zinc fingers on DNA in the dark blueprint theme"><br>
<sub><b>Metals, DNA, dark theme.</b> Zif268 zinc fingers (1ZAA), blueprint theme: each Zn tethered to its two Cys and two His.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/p53-alphafold.png" alt="p53 AlphaFold model coloured by pLDDT"><br>
<sub><b>AlphaFold models.</b> p53 (AF-P04637), alphafold theme: elements, loops and tails in pLDDT confidence bands.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ubiquitin-conservation.png" alt="Ubiquitin coloured by conservation"><br>
<sub><b>Conservation.</b> Ubiquitin (1UBQ) coloured from an alignment of the ubiquitin-like family, ConSurf grades.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/tim-barrel-rainbow.png" alt="TIM barrel dimer coloured N to C"><br>
<sub><b>Barrels.</b> Triosephosphate isomerase (1TIM), rainbow theme: the closed (βα)₈ barrel, first strand repeated as a ghost.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/haemoglobin-hemes.png" alt="Haemoglobin tetramer with hemes"><br>
<sub><b>Ligands.</b> Haemoglobin (2HHB): each heme placed beside, and tethered to, the helices that hold it.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/fima-richardson.png" alt="FimA pilin coloured by secondary structure type"><br>
<sub><b>Disulfides and 3₁₀ helices.</b> FimA (5NKT), Richardson colouring: the disulfide as a bar, the 3₁₀ helix as η1.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ubiquitin-sequence.png" alt="Ubiquitin sequence with secondary structure"><br>
<sub><b>Sequence view.</b> The same labels and colours as the topology, turns (TT) and residue numbers.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/fima-pilus-rod.png" alt="One FimA subunit of the type 1 pilus rod with its neighbours' strands"><br>
<sub><b>Filaments.</b> Type 1 pilus rod (6Y7S): one FimA subunit in full; the next subunit's donor strand A″ completes its Ig fold, its own strand A pairs with the previous subunit.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ferritin-cage.png" alt="One apoferritin subunit with the neighbouring helices of the four-fold channel"><br>
<sub><b>Cages.</b> Apoferritin (7A4M), 24 subunits: one chain in full, with the neighbours' α5 helices that line the four-fold channel.</sub></td>
</tr>
</table>

### Large assemblies

Filaments, cages and other large symmetric assemblies are drawn one subunit at a time (`--focus auto`): the
smallest repeating unit (one copy of each distinct chain) in full, plus the parts of its neighbours that belong to
its fold: strands paired with its own by backbone hydrogen bonds (donor strands, domain swaps, β-augmentation,
cross-β stacking) and helices bundled with its own. Neighbours' elements are grey and primed (A′, A″), carry no
termini or loops of their own, and a line under the legend says what is shown, e.g. *1 of 6 subunits shown ·
helical, 115.0° twist, 7.9 Å rise per subunit*. `--focus none` draws everything; `--focus B` picks the subunit.

### Sequence and alignment

<p align="center">
<img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ubiquitin-family-alignment.png" width="640"
alt="Ubiquitin-like family alignment with the structure of 1UBQ on top">
</p>

The structure of one member (here 1UBQ) drawn over an alignment of the family, ESPript-style: identical columns
white on red, similar columns red and framed, conservation bars underneath
([alignment](docs/examples/ubiquitin_family.fasta): human, yeast and barley ubiquitin, NEDD8, ISG15 and SUMO1
from UniProt, aligned with FAMSA).

### Interactive explorer

<p align="center">
<img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/interactive.png" width="760"
alt="Interactive page: topology, 3D model and contact map side by side">
</p>

`foldmap interactive` writes one HTML page with the topology, the 3D model and a residue contact map side by
side.

| Action | Effect |
| --- | --- |
| Hover an element, loop, map cell or residue in 3D | preview: the same residues light up in all three views |
| Click | select it (stays lit) and fly the 3D view to it |
| Click it again, click empty paper, or press Esc | clear the selection |
| Shift-click | add to (or remove from) the selection |
| Click a chain in the legend, or its N/C terminus | select the whole chain |
| Mouse wheel / drag / `fit` on the topology | zoom / pan / fit |

The info bar lists what is selected and how many residues. The 3D model starts in the same orientation as the figure. The 3D panel uses
[Mol*](https://molstar.org) (the PDBe build); `--viewer 3dmol` switches to the lighter
[3Dmol.js](https://3dmol.csb.pitt.edu/). Either loads from a CDN, so the page needs an internet connection.

## Install

```bash
pip install "foldmap @ git+https://github.com/dzyla/foldmap.git"              # command-line tool and Python package
pip install "foldmap[app] @ git+https://github.com/dzyla/foldmap.git"         # + the Streamlit app
pip install "foldmap[mcp] @ git+https://github.com/dzyla/foldmap.git"         # + the MCP server for AI agents
```

Python 3.11 or newer. From a clone: `pip install -e ".[dev,app,mcp]"`. (PyPI release coming: then simply
`pip install foldmap`.)

## Quick start

Every command takes a file path, a PDB ID (downloaded from RCSB) or a UniProt accession (the AlphaFold model).
Downloads are cached in `~/.cache/foldmap`.

```bash
foldmap summary 1UBQ                              # elements, sheets, symmetry, ligands
foldmap plot 1UBQ -o ubiquitin.svg                # .svg (editable text), .png (300 dpi) or .pdf; -o repeats
foldmap plot 1UBQ -o fig.png --theme richardson --set residue_numbers=true
foldmap plot P04637 -o p53.svg --theme alphafold  # AlphaFold model, pLDDT colours
foldmap sequence 1UBQ -o seq.svg --columns 50     # sequence with secondary structure
foldmap sequence 1UBQ -o aln.svg --msa family.fasta
foldmap interactive 8UTF -o explorer.html         # linked topology / 3D / contact map
foldmap app                                       # the Streamlit app
foldmap styles                                    # every theme and style key
```

## What it handles

- **Secondary structure** by a DSSP implementation: α-helices, β-strands with bulges, 3₁₀ helices, turns; β-sheets
  with strand order and direction, closed barrels shown with a ghost of the first strand.
- **Layout** from the 3D structure: elements projected onto the page, then arranged so sheets stay in order,
  contacting elements and disulfide partners stay close, and nothing overlaps; loops routed orthogonally around
  elements and labels.
- **Assemblies and symmetry**: biological assemblies are built from the file; cyclic, dihedral, cubic and helical
  symmetry is detected and every protomer drawn alike; highlight one protomer or the asymmetric unit; filaments,
  cages and large assemblies drawn as one subunit with the neighbouring parts that complete its fold.
- **Protein–DNA/RNA**: duplexes drawn as double helices with base letters and contact beads, cropped to the
  stretch the protein touches.
- **Annotations**: disulfides, glycans (SNFG symbols), ligands and metal ions, membrane bands, named domain panels,
  residue numbers, marked elements.
- **Colouring**: by chain, N→C, per-element shade, secondary-structure type, B-factor, hydropathy, AlphaFold pLDDT
  or conservation from an alignment.

## Themes

Themes change only how things are drawn, never where. `foldmap styles` lists them with every key.

| Theme | Look |
| --- | --- |
| `publication` | the default: bold fills, ribbon helices with depth shading |
| `journal` | muted palette, thinner lines, smaller text |
| `print` / `minimal` | greyscale with pale fills / outlines only |
| `presentation`, `cartoon` | glossy helices, curved loops, larger text |
| `richardson` | coloured by element type |
| `rainbow`, `shaded`, `trace` | N→C colouring; `trace` also colours loops and adds direction arrows |
| `flexibility`, `hydropathy` | B-factor or Kyte–Doolittle hydropathy |
| `alphafold` | pLDDT confidence bands, loops and tails residue by residue |
| `conservation` | ConSurf grades from `--msa` |
| `blueprint` | dark navy page, pale outlined elements |
| `goodsell` | pale, illustrative |

Change any key with `--set KEY=VALUE` (repeatable) or a YAML file for `--style-file`. Commonly used keys:
`color_by`, `palette`, `fill`, `loops`, `loop_color` (`residue` colours loops residue by residue),
`residue_numbers`, `highlight` (`asu`, `protomer` or chains), `mark` (elements in a colour of their own),
`ligands` (`auto`, `all`, `none` or codes like `HEM,ZN`), `disulfides`, `glycans`, `helices_310`, `dna_extent`
(`contacts` or `all`), `background`, `ink`, `helix_scale`, `strand_scale`, `font_scale`.

## Shaping a figure

| Option | Example |
| --- | --- |
| Swap two elements | `--swap α1,α3` |
| Nudge an element | `--move α2=1,-2` |
| Rename an element | `--rename res:A:167-182=Gd` |
| Highlight elements | `--set mark=α4=#2ca02c` |
| Named domain panels | `--domain N-lobe=res:A:13-59` (repeatable) or `--domains auto` |
| Only some chains | `--chains A,B` |
| Large assemblies | `--focus auto\|none\|protomer\|B` (one subunit plus what its neighbours add to its fold) |
| Membrane | `--membrane auto` (OPM dummy atoms if present, otherwise estimated) |
| Conservation | `--theme conservation --msa family.fasta` |
| Orientation | `--up X,Y,Z`, `--view X,Y,Z`, `--rotate DEG` |
| Symmetry and assembly | `--symmetry auto\|off\|C3\|D2\|helical`, `--protomers "A,B;C,D"`, `--assembly auto\|asu\|ID` |

Elements are named by label (`A`, `α2`, `η1`), `chain:label`, position (`#3`) or residue range (`res:A:10-25`).

**Layout files.** `--save-layout fig.yaml` records everything that defines a figure: theme, style changes,
layout options, domains and edits, with element references stored as residue ranges so they survive a new model.
`foldmap plot model.cif --layout-file fig.yaml -o fig.svg` re-renders it; the command line overrides the file.

## Sequence view

```bash
foldmap sequence 1UBQ -o seq.svg --columns 60                         # one chain
foldmap sequence 1UBQ -o aln.svg --msa family.fasta --reference UBC   # alignment (FASTA, Clustal, Stockholm)
```

Options: `--chain`, `--columns` (residues per row), `--full-sequence` (unmodelled ends too, shown in grey),
`--ss-colour figure|black`, `--similarity` (fraction that must agree for a similar column, default 0.7),
`--no-conservation-bar`, `--theme`/`--set` (element colours follow the topology figure).

## How the topology is made

Every figure (command line, app, MCP or Python) goes through the same steps:

1. **Read the structure.** gemmi reads the mmCIF/PDB file (or the file downloaded for a PDB ID or UniProt
   accession). With `--assembly auto` the biological assembly is built from the file's own operators.
2. **Secondary structure.** A DSSP implementation assigns hydrogen bonds from the backbone and calls α-helices,
   3₁₀ helices, β-strands (with bulges) and turns.
3. **Sheets and bundles.** Strands joined by backbone H-bond ladders form sheets; their order and directions come
   from the ladder network, and a sheet whose ends pair is a closed barrel. Helices packed along their length
   across chains form bundles.
4. **Context.** Symmetry (cyclic, dihedral, cubic, helical) is found by superposing chains of the same sequence;
   disulfides, glycans and ligands are read from the atoms; a membrane is estimated on request.
5. **Projection.** A frame turns 3D into page coordinates: the dominant element direction points up and the
   widest spread runs left–right (single chains); a symmetric assembly is unrolled around its axis, protein–DNA
   complexes around the DNA axis, membrane proteins with the bilayer horizontal. Never mirrored.
6. **Layout.** Sheets (as rigid blocks), helices and DNA are placed by minimising one energy: springs to their
   projected positions, short loops between sequence neighbours, attraction between elements in contact, a pull
   between disulfide-bonded cysteines, and an overlap penalty; then any remaining overlaps are cleared.
   Symmetry copies are tied to look identical.
7. **Loops.** Each loop is routed on a grid with A* search: orthogonal lines that go around elements and labels.
8. **Drawing.** matplotlib draws the elements to scale (one helix turn tiled along each helix, arrows for
   strands), then loops, annotations, legend and title, and writes SVG (text kept as text), PNG or PDF.

Steps 1–7 depend only on the structure and the layout options; themes and style keys only change step 8, so a
figure keeps its shape while you restyle it.

## App

`foldmap app` opens a Streamlit app in the browser. It runs the same pipeline as the command line:

1. **Load**: type a PDB ID, a UniProt accession (the AlphaFold model) or a file path, upload a file, or pick an
   example. A summary shows chains, residues, element counts and the symmetry found.
2. **Pick a style**: previews of *your* structure in each theme; choose one to start from.
3. **Edit**: every option sits in the sidebar (colours, loops, labels, features such as disulfides, glycans,
   ligands and membranes, sizes and layout, symmetry and assembly, large-assembly focus, alignment upload). The
   figure redraws as you change anything; results are cached, so switching back is instant.
4. **Export**: the Figure tab downloads SVG, PNG (300 dpi), PDF and a layout file that recreates the figure; the
   Sequence tab draws the sequence (or the uploaded alignment) with secondary structure on top; the Interactive
   tab embeds the explorer page and downloads it as HTML.

## MCP server for AI agents

`foldmap mcp` runs a [Model Context Protocol](https://modelcontextprotocol.io) server over stdio, so agents can
summarise structures and draw figures. Tools: `summarize_structure`, `draw_topology`, `draw_sequence`,
`interactive_page`, `list_styles`; the drawing tools return a PNG preview along with the file they wrote.

```bash
claude mcp add foldmap -- foldmap mcp        # Claude Code
```

```json
{ "mcpServers": { "foldmap": { "command": "foldmap", "args": ["mcp"] } } }
```

(the second form for clients configured by JSON, such as Claude Desktop).

## Python

```python
from foldmap.cli import make_figure
from foldmap.fetch import fetch
from foldmap.render import save
from foldmap.style import resolve_style

path = fetch("1UBQ")  # or a local file
fig = make_figure(path, title="Ubiquitin", look=resolve_style("journal", None, ["residue_numbers=true"]))
save(fig, "ubiquitin.pdf")
```

## Development

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev,app,mcp]"
.venv/bin/python -m pytest              # about 850 tests, ~1.5 min
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/python scripts/gallery.py     # re-render the images above
.venv/bin/python scripts/sweep.py --sample 60 --out sweep/   # render random PDB entries and check them
```

## Citing

If you use foldmap in published work, please cite it. GitHub's "Cite this repository" button gives the reference
in APA and BibTeX form (from [`CITATION.cff`](CITATION.cff)).

## License

foldmap is free software under the [GNU Affero General Public License v3.0 or later](LICENSE), with one additional
attribution term (see [`NOTICE`](NOTICE)). You may use, study, change and share it; any modified version you
distribute, or run as a service for others, must be released under the same license with its source code and keep
the foldmap attribution. Copyright © 2026 Dawid Zyla.
