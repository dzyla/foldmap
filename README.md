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
<sub><b>Conservation.</b> Ubiquitin (1UBQ) coloured from an alignment of the ubiquitin-like family: weighted entropy, in ConSurf's colour scale.</sub></td>
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
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ncrna-30qm.png" alt="DRT10 RT-ncRNA complex with six RNA stems drawn as ladders between the two protomers"><br>
<sub><b>RNA folds.</b> DRT10 RT–ncRNA (30QM): six stems as ladders between the two RT protomers, contacts as beads.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/nucleosome-1kx5.png" alt="Nucleosome core particle with the DNA as a chain of straight duplex pieces round the histones"><br>
<sub><b>Bent DNA.</b> Nucleosome (1KX5): 146 base pairs in 22 straight pieces round the histones.</sub></td>
</tr>
</table>

### Large assemblies

Filaments, cages and other large symmetric assemblies are drawn one subunit at a time (`--focus auto`): the
smallest repeating unit (one copy of each distinct chain) in full, plus the parts of its neighbours that belong to
its fold: strands paired with its own by backbone hydrogen bonds (donor strands, domain swaps, β-augmentation,
cross-β stacking) and helices bundled with its own. Neighbours' elements are grey and primed (A′, A″), carry no
termini or loops of their own, and a line under the legend says what is shown, e.g. *1 of 6 subunits shown ·
helical, 115.0° twist, 7.9 Å rise per subunit*. A left-handed filament has a negative twist and says so:
*−115.0° twist (left-handed)*. `--focus none` draws everything; `--focus B` picks the subunit.

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

- **Secondary structure** by a DSSP implementation that matches mkdssp 3.1.2 residue for residue on the test set:
  α-, 3₁₀- and π-helices (π before α, as in DSSP ≥ 2.1), β-strands with bulges, turns; β-sheets
  with strand order and direction, closed barrels shown with a ghost of the first strand.
- **Layout** from the 3D structure: elements projected onto the page, then arranged so sheets stay in order,
  contacting elements and disulfide partners stay close, and nothing overlaps; loops routed orthogonally around
  elements and labels.
- **Assemblies and symmetry**: biological assemblies are built from the file; cyclic, dihedral, cubic and helical
  symmetry is detected and every protomer drawn alike; highlight one protomer or the asymmetric unit; filaments,
  cages and large assemblies drawn as one subunit with the neighbouring parts that complete its fold.
- **DNA and RNA as their fold**: base pairs (Watson–Crick and G·U) from the coordinates, stems, hairpins, junctions
  and single strands placed beside the protein that holds them; bent duplexes as chains of straight pieces;
  `na_style` `ladder` (default), `ribbon`, `block` or `letters`; `na_extent` `auto` crops RNA over 300 nt to the
  stems the protein touches.
- **Annotations**: disulfides, glycans (SNFG symbols, sialic acids included), ligands and metal ions, membrane bands,
  named domain panels, residue numbers, marked elements.
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
| `conservation` | conservation from `--msa`: Henikoff-weighted Shannon entropy, z-score grades 1–9 in ConSurf's colours (not Rate4Site), or ConSurf's own grades table |
| `blueprint` | dark navy page, pale outlined elements |
| `goodsell` | pale, illustrative |

Change any key with `--set KEY=VALUE` (repeatable) or a YAML file for `--style-file`. Commonly used keys:
`color_by`, `palette`, `fill`, `loops`, `loop_color` (`residue` colours loops residue by residue),
`residue_numbers`, `highlight` (`asu`, `protomer` or chains), `mark` (elements in a colour of their own),
`ligands` (`auto`, `all`, `none` or codes like `HEM,ZN`), `na_backbone`, `backbone_colors`, `dna_color`, `rna_color`,
`hbonds`, `hbond_color`, `hbond_max`, `disulfides`, `glycans`, `helices_310`, `na_extent`, `na_style`,
`background`, `ink`, `helix_scale`, `strand_scale`, `font_scale` (`dna_extent` still works on the command line and in style and layout files).

## Shaping a figure

| Option | Example |
| --- | --- |
| Swap two elements | `--swap α1,α3` |
| Nudge an element | `--move α2=1,-2` |
| Rename an element | `--rename res:A:167-182=Gd` |
| Conventional names | `--names deposited` (the file's helix ids, e.g. globin A–H), `--names beta` (β1…βn) or `--names docs/names/2HHB-globin.yaml` (repeatable) |
| Highlight elements | `--set mark=α4=#2ca02c` |
| Named domain panels | `--domain N-lobe=res:A:13-59` (repeatable) or `--domains auto` |
| Only some chains | `--chains A,B` |
| Large assemblies | `--focus auto\|none\|protomer\|B` (one subunit plus what its neighbours add to its fold) |
| Membrane | `--membrane auto` (OPM dummy atoms if present, otherwise estimated); `--membrane flip` puts the other side on top |
| Conservation | `--theme conservation --msa family.fasta`, or `--msa consurf_grades.txt` for ConSurf's own grades |
| Orientation | `--up X,Y,Z`, `--view X,Y,Z`, `--rotate DEG` |
| Symmetry and assembly | `--symmetry auto\|off\|C3\|D2\|helical`, `--protomers "A,B;C,D"`, `--assembly auto\|asu\|ID` |

Elements are named by label (`A`, `α2`, `η1`), `chain:label`, position (`#3`) or residue range (`res:A:10-25`). Negative numbers and insertion codes are written as in the file: `res:A:-7-4`, `res:H:100A-100K`.

A label reference reads the labels the figure shows: the label copies of one protein share, the primed labels of
elements a neighbour lends in focus mode (`α2′`) and any `--names`; `--rename` is applied last and wins over all of them.

**Layout files.** `--save-layout fig.yaml` records everything that defines a figure: theme, style changes,
layout options, domains and edits, with element references stored as residue ranges so they survive a new model.
`foldmap plot model.cif --layout-file fig.yaml -o fig.svg` re-renders it; the command line overrides the file.

**What the fold drawing shows.** Base pairs are found from the coordinates (Watson–Crick and G·U, with modified
nucleotides paired by their parent base) and grouped into stems, named `P1`, `P2`, ... for RNA and `DNA 1`, ... for
DNA. Each stem is placed as a rigid item beside the protein that holds it; a duplex that bends (nucleosome, IHF, TBP)
is cut where its axis turns by 25° into straight pieces placed end to end. Single strands, hairpin loops (as arcs)
and the links between stems are routed like protein loops, with both ends marked 5′ and 3′, base letters, and contact
beads in the colour of the protein chain that touches them. A nucleic acid with no pairs at all (a G-quadruplex, one
strand on a protein) is drawn as a projected trace with letters. `foldmap summary` lists each stem with its base
pairs and ranges, e.g. `P1  7 bp  1:1-7 · 1:112-118`.

`na_style` chooses how stems are drawn (layout is unaffected): `ladder` (rails and base-pair rungs; the default),
`ribbon` (a double helix; the `presentation` and `cartoon` themes), `block` (bars; `minimal` and `print`) or
`letters` (only the base letters). A large RNA cropped by `na_extent auto` switches `ladder` to `block`. `na_extent`
chooses what is drawn and placed (a layout key): `auto` (everything up to 300 nt, otherwise the stems the protein
touches, with gaps marked), `contacts` (the stretch the protein touches, with 6 nt beyond) or `all`. `dna_extent` is
an alias of `na_extent`.

**Nucleic-acid colours.** DNA backbones are brown (`dna_color`, `#8c564b`) and RNA backbones purple (`rna_color`,
`#6a3d9a`), so they read apart from the black protein loops; contact beads keep the colour of the protein chain that
touches them. `na_backbone` chooses the scheme: `modality` (the default: `dna_color` and `rna_color`), `chain` (one
colour per nucleic-acid chain), `ink` (as the loops) or a single colour for all nucleic acid (`#rrggbb` or a name).
`backbone_colors` overrides single backbones: a comma-separated list of `KEY=COLOUR`, where KEY is a chain id, or
`protein`, `dna` or `rna`. A chain key beats a modality key, and `protein` or a protein chain id recolours the loops.
Chain ids must exist in the figure (an unknown one is an error that lists the chains). Example, on a figure with
protein chain A and DNA chains I and J (1KX5):

```bash
foldmap plot 1KX5 --set 'backbone_colors=A=#1f77b4,dna=#555555,I=crimson' -o 1kx5.svg
```

**Hydrogen bonds.** `--set hbonds=na` draws dashed lines for the hydrogen bonds between protein and DNA/RNA,
`interchain` between all different chains, `chains:A+B,C+D` between the listed chain pairs (default `off`). Bonds are
found from the heavy atoms (donor and acceptor 2.5-3.5 Å apart, with the donor's hydrogen able to point at the
acceptor), then drawn as one thin line per residue pair, from the residue's point on its element or loop to the
nucleotide's (or partner residue's) point, under the contact beads, elements and letters. Long lines fade (down to
40% opacity beyond half the page diagonal). At most `hbond_max` lines are drawn (default 200); when more are found
the shortest *on the page* are kept (not the shortest in Å) and the legend reads "H-bonds (n of m)", with
", k not shown" for pairs whose residue is not drawn (cropped nucleic acid, other chains). `hbond_color` sets the
colour (the default is chosen per theme to stand clear of the chain colours, the backbones and the ink). Only drawing
changes, never the layout. Example: `foldmap plot 1LMB.cif --set hbonds=na -o 1lmb.svg`.

On large complexes the lines can pile up into a hairball (30QM with `hbonds=na` draws 79). Restrict them with
`chains:` (`--set hbonds=chains:A+B`) or lower `hbond_max` (`--set hbond_max=30`).

**Checking the contact beads.** `foldmap contacts 1LMB.cif` draws the figure `plot` would (it takes the same options),
checks every contact bead against the deposited atoms with an independent search, and prints how many beads were
drawn, verified, false, missing, unmatched, in the wrong colour or not on the first toucher, plus a line for each
problem bead. It exits 0 when clean and 1 on any false, unmatched or wrong-colour bead; `-o fig.svg` also saves the
figure.

## UniProt annotation

`--uniprot auto` finds each chain's UniProt entry (PDBe SIFTS mapping for PDB entries, the accession for AlphaFold
models and files that name one; `.pdb` files map the same way, with or without SEQRES; or give one,
`--uniprot P07911`) and maps its features onto the structure, correcting for numbering offsets. Active and binding sites are marked on the figure (`--set uniprot_sites=all` for every site and
modification), `--domains uniprot` turns UniProt domains into panels, and `foldmap sequence ... --uniprot auto`
adds a feature track.

```bash
foldmap uniprot P07911                             # domains, sites, every PDB structure, the AlphaFold model
foldmap plot 2HHB --uniprot auto --domains uniprot -o hb.svg
foldmap sequence 6ZS5 --chain A --uniprot auto -o umod.svg
```

In the app, typing an accession lists its AlphaFold model and experimental structures (widest coverage first,
with method and resolution) to load with one click.

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
2. **Secondary structure.** A DSSP implementation assigns hydrogen bonds from the backbone and calls α-, 3₁₀- and
   π-helices, β-strands (with bulges) and turns, exactly as mkdssp 3.1.2 does on every test structure.
3. **Sheets and bundles.** Strands joined by backbone H-bond ladders form sheets; their order and directions come
   from the ladder network, and a sheet whose ends pair is a closed barrel. Helices packed along their length
   across chains form bundles.
4. **Context.** Symmetry (cyclic, dihedral, cubic, helical) is found by superposing chains of the same sequence;
   disulfides, glycans and ligands are read from the atoms; a membrane is estimated on request. Nucleic acids:
   base pairs grouped into stems; long bent duplexes cut where the axis turns by 25°.
5. **Projection.** A frame turns 3D into page coordinates: the dominant element direction points up and the
   widest spread runs left–right (single chains); a symmetric assembly is unrolled around its axis, protein–DNA
   complexes around the DNA axis, membrane proteins with the bilayer horizontal. Never mirrored.
6. **Layout.** Sheets (as rigid blocks), helices and DNA are placed by minimising one energy: springs to their
   projected positions, short loops between sequence neighbours, attraction between elements in contact, a pull
   between disulfide-bonded cysteines, and an overlap penalty; then any remaining overlaps are cleared.
   Symmetry copies are tied to look identical; stems are placed as rigid items.
7. **Loops.** Each loop is routed on a grid with A* search: orthogonal lines that go around elements and labels;
   single strands are routed the same way, hairpin loops as arcs.
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

Hosting the app for other people: set `FOLDMAP_APP_LOCAL_FILES=0` so the app accepts PDB IDs, UniProt accessions
and uploads but never opens a file path on the server (including alignment paths in uploaded layout files).
Uploads are kept in `~/.cache/foldmap/uploads` for a day.

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

Figures must end in `.svg`, `.png` or `.pdf`, pages in `.html` and layout files in `.yaml` or `.yml`. Set
`FOLDMAP_OUTPUT_DIR` in the server's environment to confine the figure, page and layout files the server writes (default outputs included) to one folder;
the structure download cache (`~/.cache/foldmap`) is not covered.

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
.venv/bin/python scripts/bench.py                                # time layouts (add --record/--check to compare runs)
```

## Citing

If you use foldmap in published work, please cite it:

> Zyla, D. *foldmap: publication-grade protein topology diagrams* (2026). https://github.com/dzyla/foldmap

```bibtex
@software{zyla_foldmap_2026,
  author  = {Zyla, Dawid},
  title   = {foldmap: publication-grade protein topology diagrams},
  year    = {2026},
  url     = {https://github.com/dzyla/foldmap},
  license = {AGPL-3.0-or-later}
}
```

GitHub's "Cite this repository" button gives the same reference (from [`CITATION.cff`](CITATION.cff)). Each
release is archived on Zenodo with its own DOI; please cite the version you used.

## License

foldmap is free software under the [GNU Affero General Public License v3.0 or later](LICENSE), with one additional
attribution term (see [`NOTICE`](NOTICE)). You may use, study, change and share it; any modified version you
distribute, or run as a service for others, must be released under the same license with its source code and keep
the foldmap attribution. Copyright © 2026 Dawid Zyla.

Figures and interactive pages that foldmap generates are yours to license as you like (see [`NOTICE`](NOTICE), output
exception); the credit line on interactive pages stays.
