# foldmap

[![tests](https://github.com/dzyla/foldmap/actions/workflows/tests.yml/badge.svg)](https://github.com/dzyla/foldmap/actions/workflows/tests.yml)
[![license: AGPL-3.0](https://img.shields.io/badge/license-AGPL--3.0-blue.svg)](LICENSE)
![python](https://img.shields.io/badge/python-3.11%2B-blue.svg)

**Publication-grade topology diagrams of proteins, DNA and RNA, straight from a PDB/mmCIF file, a PDB ID or an
AlphaFold model.**

**Try it in the browser: <https://foldmap.streamlit.app>** · **Website and theme gallery:
<https://dzyla.github.io/foldmap/>**

```bash
pip install git+https://github.com/dzyla/foldmap.git
foldmap plot 1LMB -o lambda.svg
```

foldmap reads the structure, assigns secondary structure, finds β-sheets, helix bundles, base pairs and symmetry,
and lays the elements out as a clean 2D map: coiled-ribbon helices, strand arrows and 3₁₀ helices on one consistent
scale, DNA and RNA as stems and ladders beside the protein that holds them, with loops and strands routed around
everything else. The same input always gives the same figure. It also draws the sequence with the secondary
structure on top, alone or over a multiple sequence alignment.

New in this version:

- **DNA and RNA drawn as their fold**: stems, hairpins, pseudoknots and junctions, nucleosome DNA round its histones,
  contact beads checked against the file, and protein–nucleic acid hydrogen bonds on request.
- **Figures that follow the published conventions**: TIM barrels and Rossmann folds as a βα row, kinked helices
  split (globin A/B), bulged strands kept whole, chromophores kept in the chain, symmetric dimers drawn as mirror
  images (a two-fold assembly of one-chain copies seen from the side of its dyad, one copy the other's mirror image), conventional names (`--names`) and domain colours with a domain bar.
- **A redesigned app** at <https://foldmap.streamlit.app>: one box to start, share links, light and dark themes,
  every option in tidy sections, downloads built from the figure on screen.

<table>
<tr>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/lambda-repressor-dna.png" alt="Topology of the lambda repressor dimer bound to its operator DNA"></td>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ncrna-30qm.png" alt="DRT10 reverse transcriptase with its non-coding RNA drawn as six stems"></td>
</tr>
</table>

## Gallery

<table>
<tr>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/t4-lysozyme-domains.png" alt="T4 lysozyme with N- and C-lobe domain panels"><br>
<sub><b>Domains.</b> T4 lysozyme (2LZM), journal theme: named domain panels, each kept together, elements in domain colours and a domain bar in residue order.</sub></td>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/kcsa-c4.png" alt="KcsA potassium channel tetramer with one protomer highlighted"><br>
<sub><b>Symmetry.</b> KcsA (1BL8): C4 found, the four protomers drawn alike, one highlighted, the K⁺ ions of the filter tethered to the chains that hold them.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/bacteriorhodopsin-membrane.png" alt="Bacteriorhodopsin trimer in a membrane band"><br>
<sub><b>Membranes.</b> Bacteriorhodopsin (1C3W): the bilayer found from the structure, sides by the positive-inside rule, retinal bound.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/zinc-fingers-blueprint.png" alt="Zif268 zinc fingers on DNA in the dark blueprint theme"><br>
<sub><b>Metals, DNA, dark theme.</b> Zif268 zinc fingers (1ZAA), blueprint theme: each Zn tethered to its two Cys and two His.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/p53-alphafold.png" alt="p53 AlphaFold model coloured by pLDDT"><br>
<sub><b>AlphaFold models.</b> p53 (AF-P04637), alphafold theme: elements, loops and tails in pLDDT confidence bands; α1 moved clear of the core with <code>--move</code>.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ubiquitin-conservation.png" alt="Ubiquitin coloured by conservation"><br>
<sub><b>Conservation.</b> Ubiquitin (1UBQ) coloured from an alignment of the ubiquitin-like family: weighted entropy, in ConSurf's colour scale.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/tim-barrel-rainbow.png" alt="TIM barrel chain A drawn as one beta-alpha row with strands named beta1 to beta8"><br>
<sub><b>βα rows.</b> Triosephosphate isomerase (1TIM, chain A), rainbow theme, <code>--names beta</code>: the closed (βα)₈ barrel as one row, each helix between the strands it joins, opened between β8 and β1, β1 repeated as a ghost after β8.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/myoglobin-names.png" alt="Myoglobin with its helices named A to H"><br>
<sub><b>Conventional names.</b> Myoglobin (1MBN), <code>--names deposited</code>: helices A–H as the file names them (a short helix it leaves unnamed keeps its generic α4), the heme tethered to the helices that hold it.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/haemoglobin-hemes.png" alt="Haemoglobin tetramer with hemes"><br>
<sub><b>Ligands.</b> Haemoglobin (2HHB): each heme placed beside, and tethered to, the helices that hold it.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/gfp-chromophore.png" alt="Green fluorescent protein barrel with the chromophore starred"><br>
<sub><b>Chromophores.</b> GFP (1EMA): the chromophore CRO66 kept in the chain and marked with a star; the 11-strand barrel closed by a ghost of strand H.</sub></td>
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
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/trna-1ehz.png" alt="Yeast tRNA-Phe with its four arms as ladders, named acceptor, D, anticodon and T"><br>
<sub><b>RNA folds.</b> tRNA-Phe (1EHZ), arms named from <code>docs/names/1EHZ-trna.yaml</code>: four stems as ladders with their loops and base letters, the tertiary pair marked by tie letters.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ncrna-30qm.png" alt="DRT10 RT-ncRNA complex with six RNA stems drawn as ladders beside and between the two protomers"><br>
<sub><b>RNA–protein complexes.</b> DRT10 RT–ncRNA (30QM): the 121-nt ncRNA's six stems (P1–P6) as ladders beside and between the two RT protomers, contact beads in the colour of the chain that touches them.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/lambda-hbonds.png" alt="Lambda repressor dimer on DNA with protein-DNA hydrogen bonds as dashed lines"><br>
<sub><b>Hydrogen bonds.</b> λ repressor (1LMB), <code>--set hbonds=na --set hbond_max=16</code>: the 16 shortest of 30 protein–DNA hydrogen bonds as dashed lines.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/nucleosome-1kx5.png" alt="Nucleosome core particle with the DNA as a chain of straight duplex pieces round the histones"><br>
<sub><b>Bent DNA.</b> Nucleosome (1KX5): 146 base pairs in 22 straight pieces round the histones.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/measles-f-trimer.png" alt="Measles fusion protein trimer with one protomer highlighted"><br>
<sub><b>Symmetry.</b> Measles F trimer (8UTF): C3 found, protomers drawn alike, one highlighted, disulfides and glycans marked.</sub></td>
<td></td>
</tr>
</table>

### Large assemblies

Filaments, cages and other large symmetric assemblies are drawn one subunit at a time (`--focus auto`): the
smallest repeating unit (one copy of each distinct chain) in full, plus the parts of its neighbours that belong to
its fold: strands paired with its own by backbone hydrogen bonds (donor strands, domain swaps, β-augmentation,
cross-β stacking) and helices bundled with its own. Neighbours' elements are grey and primed (A′, A″), carry no
termini or loops of their own, and a line under the legend says what is shown, e.g. *1 of 6 subunits shown ·
helical, 115.0° twist, 7.9 Å rise per subunit*. A left-handed filament has a negative twist and says so:
*−115.0° twist (left-handed)*. `--focus none` draws everything; `--focus B` picks the subunit. A two-subunit
assembly (a C2 histone octamer) is never cut in half, nor a ring whose subunits share one DNA or RNA (nucleosome DNA):
those are drawn whole.

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
| Click | select it (stays lit); a click in the figure also brings the 3D view to it, a click in 3D keeps the camera |
| Click within 4 px of a loop or a DNA/RNA strand or stem | select it, at any zoom (where a loop runs beside a strand, the loop) |
| Hover or click an interaction line (H-bond, salt bridge, π-stacking, cation-π, hydrophobic) | its residues, atoms, distances (and angle, for stacking and cation-π) in the info bar; a click selects the line (within 4 px) and lights its residue and nucleotide in 3D |
| Click a map cell | select the whole helix or strand of each of its two residues (for a loop residue, the nearest one along its chain) and bring the 3D view to them |
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
foldmap summary 1UBQ                              # elements, sheets, science checks, symmetry
foldmap plot 1UBQ -o ubiquitin.svg                # .svg (editable text), .png (300 dpi, less past 50 MP) or .pdf
foldmap plot 1UBQ -o fig.png --theme richardson --set residue_numbers=true
foldmap plot P04637 -o p53.svg --theme alphafold  # AlphaFold model, pLDDT colours
foldmap sequence 1UBQ -o seq.svg --columns 50     # sequence with secondary structure
foldmap sequence 1UBQ -o aln.svg --msa family.fasta
foldmap interactive 8UTF -o explorer.html         # linked topology / 3D / contact map
foldmap contacts 1LMB                             # check the contact beads against the atoms
foldmap check 9FWZ                                # check the figure's science and page (S1-S6, A1-A20)
foldmap data 1TIM -o tim.json                     # the figure as data: elements, labels, sheets, barrel n and S
foldmap app                                       # the Streamlit app
foldmap styles                                    # every theme and style key
```

## What it handles

- **Secondary structure** by a DSSP implementation checked against mkdssp 3.1.2. It assigns no bends (S); with
  mkdssp's bends read as none, 268,834 of the 268,835 residues of the random-100 campaign's 96 protein entries got
  mkdssp's code (one residue differs, in 8K8E; the campaign requires 99% per entry):
  α-, 3₁₀- and π-helices (π before α, as in DSSP ≥ 2.1), β-strands with bulges, turns; β-sheets
  with strand order and direction, closed barrels shown with a ghost of the first strand. A helix kinked by more
  than 30° is split into pieces drawn end to end; a strand that a bulge splits stays one strand when both pieces
  pair with the same partner. Copies of one protein are drawn alike: a kink split or bulge join in one copy is made
  in all, and a helix that DSSP breaks by one or two residues in some copies and not in others is drawn as most
  copies have it (a tie breaks it); the figure's scope facts count these cuts and joins.
- **Layout** from the 3D structure: elements projected onto the page, then arranged so sheets stay in order,
  contacting elements and disulfide partners stay close, and nothing overlaps; loops routed orthogonally around
  elements and labels. The layout is deterministic: the same figure on every machine where the same font is
  installed (the legend is sized from the font the figure's text resolves to, which `foldmap check` records; CI and
  the hosted app use Liberation Sans). A figure is at most 7.2 in wide (double column) while its text stays at 6 pt
  or larger; a larger figure's page grows instead (up to 48 in on its longer side, past which the whole figure is
  scaled to fit and the notes say its text is under 6 pt), and every keep-out round text and every legend row is
  sized from the text as drawn. Its PNG is 300 dpi up to 50
  megapixels and takes a lower dpi past that (said in the notes); SVG and PDF keep the page whole. The 6 pt floor is
  at text size 1: a smaller text size (`--set font_scale=0.8`, the app's Text size, the journal theme's 0.9) draws
  the text that much smaller.
- **How each figure was made** comes with it: one list of what it infers, leaves ambiguous or moves for
  readability — domains inferred from contacts (or CATH's, UniProt's, or given), strands drawn off a sheet's run and
  pairings only named, strands joined across a β-bulge, helix pieces matched across copies, ghost strands, names
  moved off their usual spot, additives not marked, ligand marks with no clear spot, loops with no clear route or
  over a label, a page scaled down or a PNG under 300 dpi. The command line prints it ("How this figure was made:"),
  the app shows it under the figure, the MCP result carries it as structured `provenance`, and the files keep it
  (the SVG's `<desc>`, the PDF's Subject, the PNG's Description); in Python it is `Drawn.provenance`. An empty list
  means the standard method alone.
- **Fold conventions** of the published figures: all-parallel α/β sheets (TIM barrels, Rossmann folds) as one βα
  row, conventional names (`--names`), helices drawn over the sheet they pack on (opt-in); see
  [Conventions of the published figures](#conventions-of-the-published-figures).
- **Assemblies and symmetry**: `--assembly auto` builds the file's first biological assembly when it adds copies
  (up to 24 chains) or when the file holds several copies of one complex, and never one that leaves out the DNA or
  RNA; cyclic, dihedral, cubic and helical symmetry is detected and the protomers are drawn alike (a two-fold
  assembly of one-chain copies seen from the side of its dyad, the two copies drawn as mirror images; interleaved protomers such as the
  measles F trimer's too, a strand lent across the cut of the unrolled ring drawn past its end, with a dashed ghost
  where it pairs, or without with `--set copy_ghosts=off`); highlight one protomer or the asymmetric unit;
  filaments, cages and large assemblies drawn as one subunit with the neighbouring parts that complete its fold.
  When either automatic choice draws part of the file, foldmap says so (a note on stderr, a bar in the app, a line
  in the MCP result) and names the setting that draws everything (`--assembly asu`, `--focus none`).
- **DNA and RNA as their fold**: base pairs from the coordinates, stems, hairpins, pseudoknots, junctions and single
  strands placed beside the protein that holds them; bent duplexes as chains of straight pieces; protein–DNA/RNA
  hydrogen bonds on request; see [RNA and DNA](#rna-and-dna).
- **Annotations**: disulfides, glycans (SNFG symbols, sialic acids included), ligands and metal ions, chromophores
  and other residues fused into the chain (starred), membrane bands, named domain panels with a domain bar, residue
  numbers, marked elements. A small molecule held by the nucleic acid alone (a riboswitch's ligand) is marked and
  tethered to the stems that hold it; metal ions held by the nucleic acid alone are not marked (a ribosome holds
  hundreds), and the note names them.
- **Colouring**: by chain, N→C, per-element shade, secondary-structure type, domain, B-factor, hydropathy,
  AlphaFold pLDDT or conservation from an alignment. Chains take the palette's colours in order; past its first
  pass (the ninth chain on the default palette) each chain takes a lighter or darker step of the colour it repeats,
  chosen to stand clear (CIEDE2000 ≥ 10, as far as the palette's hues allow: 24 chains on the default palette) of
  every colour already used, the white page and the helices' back faces, so chains do not share a colour.

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
`ligands` (`auto`, `all`, `none` or codes like `HEM,ZN`), `fold_template`, `helix_on_sheet`, `two_fold`, `cyclic`,
`copy_ghosts`, `domain_bar`, `domain_colours`, `na_backbone`, `backbone_colors`, `dna_color`, `rna_color`, `hbonds`,
`hbond_color`, `hbond_max`, `interactions`, `interactions_between`, `interaction_max`, `disulfides`, `glycans`, `helices_310`, `na_extent`, `na_style`, `background`, `ink`,
`helix_scale`, `strand_scale`, `font_scale`, `na_line` (`dna_extent` still works on the command line and in style and layout
files).

## Shaping a figure

| Option | Example |
| --- | --- |
| Swap two elements | `--swap α1,α3` |
| Nudge an element | `--move α2=1,-2` |
| Rename an element | `--rename res:A:167-182=Gd` |
| Conventional names | `--names deposited` (the file's helix ids, e.g. globin A–H), `--names beta` (β1…βn) or `--names docs/names/2HHB-globin.yaml` (repeatable) |
| Highlight elements | `--set mark=α4=#2ca02c` |
| Named domain panels | `--domain N-lobe=res:A:13-59` (repeatable), `--domains auto` (inferred), `--domains cath` (CATH's published domains, needs network) or `--domains uniprot` |
| Domain colours and bar | `--set color_by=domain`, `--set domain_bar=on\|off`, `--set domain_colours=false` |
| βα row | `--set fold_template=off` draws an all-parallel α/β sheet as projected |
| Helices over a sheet | `--set helix_on_sheet=auto` (opt-in) |
| Two-fold assemblies | `--set two_fold=mirror` or `unrolled`: a C2 assembly's copies always as mirror images, or side by side (default `auto`: mirrored where each copy is one protein chain) |
| Cyclic assemblies | `--set cyclic=rosette` (opt-in): a four- or seven-fold ring (C4, C7) drawn whole, seen down its axis, each copy the first one's drawing turned about it (default `unrolled`: side by side) |
| Strands lent across a ring's cut | `--set copy_ghosts=off` draws such a strand once, past the row's end, with no dashed copies of it or its partners |
| Hydrogen bonds | `--set hbonds=na`, `interchain` or `chains:A+B`; `--set hbond_max=30` |
| Interactions | `--set interactions=hbond,salt` (or `all`), `--set interactions_between=chains:A+B` |
| Only some chains | `--chains A,B` |
| Draw like a reference | `--like 2HHB --like-chain B` (one chain of each; see below) |
| Large assemblies | `--focus auto\|none\|protomer\|B` (one subunit plus what its neighbours add to its fold) |
| Membrane | `--membrane auto` (OPM dummy atoms if present, otherwise estimated); `--membrane flip` puts the other side on top. When the layout leaves a transmembrane element off the band, the view is laid out again with every one at its projected height and the band drawn at the membrane's faces |
| Conservation | `--theme conservation --msa family.fasta`, or `--msa consurf_grades.txt` for ConSurf's own grades |
| Orientation | `--up X,Y,Z`, `--view X,Y,Z`, `--rotate DEG` |
| Symmetry and assembly | `--symmetry auto\|off\|C3\|D2\|helical`, `--protomers "A,B;C,D"`, `--assembly auto\|asu\|ID` |

Elements are named by label (`A`, `α2`, `η1`), `chain:label`, position (`#3`) or residue range (`res:A:10-25`). Negative numbers and insertion codes are written as in the file: `res:A:-7-4`, `res:H:100A-100K`.

References are checked before anything is drawn, the same way on the command line, in the app (under the field) and
in the MCP server: a range written without `res:` (`--domain N=A:1-60`) is read as `res:A:1-60` with a warning; a
range that does not parse (`A:x-y`), a chain the structure lacks (the message lists its chains), a chain left out by
`--chains` or `--focus` (the message says which), or a range that meets none of the chain's residues
(`res:A:900-950` on a 164-residue chain) is an error (`foldmap: error: --domain N=A:x-y: …`, exit 2).

A label reference reads the labels the figure shows: the label copies of one protein share, the primed labels of
elements a neighbour lends in focus mode (`α2′`) and any `--names`; `--rename` is applied last and wins over all of them.

An element can also be named as `--names` would name it, without applying those names: `dep:E` (the deposited helix
name: myoglobin's helix E) or `beta:α4` (the β-scheme's), each also with a chain (`beta:A:α4`), in `--rename`,
`--swap`, `--move`, `--domain` and `mark`, and in the MCP server. A name without a prefix is the figure's label first
(what is on screen), then the deposited name, then the β-scheme's; one that means different elements in two of them is
refused with the choices: on 1TIM, `--swap α4,α1` gives `foldmap: error: α4 is the figure's α4 (A 96–102) and the
β-scheme's α4 (A 106–118); write fig:α4 or beta:α4` (exit 2; `fig:` is the label drawn). What a layout file stores
does not change: residue ranges. In `--domain`, a βα-row helix has no β-scheme name yet (the rows are laid out
after the domains are read); use its deposited name, its label or its range. In the app, an element picker's options
show every name an element has (1TIM: `α5 · helix αD2 · β-scheme α4 · A 106–118`), so typing any of them finds it;
the cell keeps the residue range.

**Automatic domains are inferred, not published.** `--domains auto` cuts each chain where it has few residue
contacts across (normalized cuts over CA atoms within 12 Å; a domain may be two pieces of the chain, as in 5JDN,
whose two domains alternate along it). It is a guess from the coordinates, and the figure says so: the domain bar reads
"domains inferred from contacts", and the command line, app and MCP notes say "inferred from residue contacts, not a
published domain assignment". On 33 multi-domain chains chopped by CATH it matches 85 % of the residues to the right
domain (calling each chain one domain matches 54 %) and finds 88 % of CATH's domain boundaries within 20 residues
(boundary F1 0.79). Its parameters were tuned on 22 of those chains, the other 11 held out (drawn chain by chain, so
some had a homologue among the 22). Eleven more chains, from CATH superfamilies none of the 33 has, were picked
before foldmap ran on them and never tuned on: on those it matches 86 % of the residues and finds 85 % of the
boundaries (F1 0.78). The calibration is `tests/test_domain_score.py`. For a
published assignment, `--domains cath` draws CATH's chopping of a PDB entry (fetched once from the PDBe API, panels named by CATH id such as `1ldmA01`,
captioned "domains: CATH"); when CATH has none the note says so and the figure is drawn without domains.

**Domains and βα rows.** `--domains auto` never splits a sheet drawn as a βα row (a TIM barrel, a Rossmann fold):
the row and its helices fall in one domain. A `--domain` set of your own that splits a row is drawn as given, but the
row is drawn as one unit, so its domains' panels overlap it. To let the parts move apart, draw the sheet as
projected with `--set fold_template=off`.

### Drawing a structure like another

`--like REFERENCE` draws a structure the way a reference of the same fold is drawn, so that two related proteins
can be compared side by side. The reference is laid out first; its chain is aligned to this structure's; each
element matched to one of the reference's takes that element's place, angle and name, and also its colour when
elements are coloured by identity (`chain`, `sequence`, `shade`, `sstype`, `domain`). Colourings by data
(`bfactor`, `plddt`, `conservation`, `hydropathy`) keep this structure's own values. Each sheet keeps its own strand
order, and overlaps are cleared as usual.

```bash
foldmap plot 1MBN --like 2HHB --like-chain B -o myoglobin-like-beta.svg     # myoglobin drawn like haemoglobin β
foldmap plot 2HHB --chains A --like 2HHB --like-chain B -o alpha-like-beta.svg
```

- **One chain of each.** The chain drawn is the one `--chains` names, else the structure's first protein chain (the
  notes say it was chosen, not asked for). The reference's chain is `--like-chain`, else its first protein chain
  (the note says "chain A (its first protein chain)"). Naming two protein chains with `--chains` or `--like-chain` is
  a usage error (exit 2).
- **Alignment.** `--like-match auto` (the default) aligns by sequence at 40 % identity or more, and otherwise by
  structure (TM-align's method). You can also force `sequence` or `structure`, or give a two-row alignment file with
  the reference's chain first.
- **The reference's own layout.** `--like-layout A.yaml` lays the reference out with a saved layout's options and
  edits.
- **Names.** An element only this structure has is named after the matched element of its type (strand, α- or
  3-10 helix) before it, with a letter (α3a, α3b). After a strand named by a capital letter it takes a prime instead
  (C′, C″), as Ig strands do. Where this structure splits a helix the reference has whole, the second part is primed
  (α6, α6′). A helix matched to one of the other type keeps its own type's sign with its mate's number (an α-helix
  matched to η1 is α1, never η1); where the reference already uses that name, it is named as an extra is, and the
  note says which mate it stands for ("named by their own helix type: α2a for η1").
- **The notes** say what was done, for example "drawn like 2HHB.cif chain A: 7 of 9 elements matched (sequence
  alignment, TM-score 0.89 over 139 residues); only in this structure: α3a, α3b".
- **When it is refused.** A structure that is not the reference's fold (TM-score under 0.5, or fewer than half the
  elements matched) is drawn on its own, and the note says why: "not drawn like 1UBQ.cif: not the same fold
  (TM-score 0.14; 4 of 13 elements matched)". When fewer than 70 % of the shorter chain's residues align, the note
  says how much did ("only 62% of the shorter chain aligned": a hinge between domains, or a circular permutation).
- **Frames that cannot be carried.** A figure is drawn like a reference through the reference's page frame, so only
  a plain one carries (a membrane view, `--up`, `--view`, `--rotate` and `--flip-v` do). With `--mode stack`, or a
  chain drawn on a straight duplex, on either side, the structure is drawn on its own, and the note says why and
  what would draw it: "not drawn like 1UBQ.cif: this structure's frame (--mode stack) cannot be carried; draw it
  with --mode projected".
- **Saved and shared.** `--save-layout` records the reference as given, its chain and the match (`like:` under
  `layout:`), so `--layout-file` draws the same figure again. The MCP tool `draw_topology` takes `like_source`,
  `like_chain` and `like_match`.

#### Two structures on one page

`foldmap compare A B` draws one protein chain of A on the left and one of B drawn like it on the right, both at one
page scale (the larger panel's), with one legend and the note under both:

```bash
foldmap compare 2HHB 1MBN --chain-a B -o globins.svg --json mates.json   # haemoglobin β, myoglobin drawn like it
```

`--chain-a` and `--chain-b` name the chains (else each its first protein chain), `--match` aligns them as
`--like-match` does, `--assembly` reads both files as `plot --assembly` does (`asu` keeps the deposited chain
names: under `auto` a chain of a built assembly is A1, A2, ...), and `--theme`, `--style-file` and `--set` style both
panels. A pair that is not one fold is
drawn as two figures of their own, side by side, with the note saying why. Each panel keeps its own ids, prefixed
`a/` and `b/`. `--json` writes what was matched: `alignment` (`method`, `tm_score`, `rmsd`, `identity`, and `pairs`
of residue numbers, A's then B's), `mates`, `pieces` (of a helix A has whole), `extras` (B's elements A lacks) and
`only_a` (A's elements B lacks), each element by its label as drawn and its residues; numbers to 0.001. The MCP
tool `compare_structures` (with `assembly`) returns the same page's preview and the same record.

### The figure as data

`foldmap data` lays the structure out as `plot` would (without drawing it) and writes what the figure says about it
as JSON, under the figure's own names:

```bash
foldmap data 1TIM -o tim.json                        # every section
foldmap data 1C3W --membrane auto --sections membrane,elements
foldmap data --schema                                # the format, as a JSON Schema
```

- `structure`: the file, the protein and nucleic-acid chains read, the subunit drawn under `--focus`;
- `elements`: each helix, 3-10 helix and strand with its first and last residue (chain, number, insertion code), its
  label as the figure writes it, its page place (`x`, `y`, `angle` N→C in degrees, `length`, `width`, page units),
  its sheet and domain;
- `sheets`: each sheet's strands left to right as drawn with their directions, each neighbouring pair's kind (`P`,
  `A`) and residue register (the DSSP bridges, a residue of each strand), accessory strands beside their partner,
  and for a closed sheet the barrel's strand number `n` and shear number `S` (a residue followed across each
  pairing's register once round the barrel comes back `S` residues along: 1TIM's (βα)8 barrels give n 8, S 8);
- `disulfides`, each end's residue and element;
- `symmetry`: the symmetry found (group, protomers, CA RMSD, axis; twist and rise of a helical one), how the figure
  draws the copies, and the stoichiometry (protein chains as copies of one protein, e.g. `A2B2`);
- `membrane` (with `--membrane`): the transmembrane elements (whose CA run reaches both faces, as A18 counts them)
  with the side each enters from and leaves to, every other element's side (`out`, `in`, `slab`), each chain's
  termini, and the band on the page; `out` is the side the membrane's normal points to (`sided`: whether the
  positive-inside rule settled it).

It takes `plot`'s `--chains`, `--assembly`, `--symmetry`, `--focus`, `--names`, `--rename`, `--domain`,
`--domains`, `--uniprot`, `--membrane`, `--layout-file` and `--set` (layout keys), so its labels and places are
those of the figure drawn with the same options. `--sections` keeps some of the sections; `-o` writes a file (else
standard output). The format is versioned (`"foldmap_data": 1`); keys are sorted, lists come in a fixed order and
numbers are rounded to 0.001, so the same structure and options give the same bytes.

### Conventions of the published figures

foldmap follows the conventions of the hand-drawn topology figures readers know (Richardson 1981, the TIM and
Rossmann schematics, Burroughs 2007, Nishimasu 2014). Each is a layout or style key; the defaults are below.

**Conventional names** (`--names`, opt-in). `--names deposited` takes the helix ids the file itself gives, so a
globin's helices read A–H (an element the file gives no id keeps its generic label). `--names beta` numbers the
strands β1…βn in sequence and names each helix of a βα row after the strand before it (β1 α1 β2 α2 …); the
chain's other helices take the name of the strand before them primed (α4′, α4″, α4‴), or lettered (α10a, α10b …)
when a strand is followed by more helices than three primes can tell apart.
`--names FILE.yaml` names elements by residue range and, under `stems:`, nucleic-acid stems by one nucleotide they
pair. The option repeats, and a later spec wins. `docs/names/` ships four files: globin helices (2HHB), FN3
strands (1TEN), histone-fold helices (1KX5) and tRNA arms (1EHZ). Each file says where its ranges come from.

```bash
foldmap plot 1MBN --names deposited -o myoglobin.svg
foldmap plot 1TIM --names beta -o tim.svg
foldmap plot 1KX5 --focus A,B --names docs/names/1KX5-histones.yaml -o nucleosome.svg
```

**The βα row** (`fold_template`, default `auto`). An all-parallel α/β sheet of 4 to 10 strands, such as a TIM
barrel or a Rossmann fold, is drawn as one row, with each crossover helix between the two strands it joins. A barrel
opens between its last and first strand and shows the first again as a dashed ghost. Copies of one protein share
the choice. `--set fold_template=off` draws the sheet as projected.

**Helices over a sheet** (`helix_on_sheet`, default `off`). With `--set helix_on_sheet=auto`, a helix that packs on
a sheet's face is drawn over the sheet plate, semi-transparent, as it lies (ferredoxin-like folds, ubiquitin's α1),
instead of beside it. The option is off for symmetric figures, for sheets drawn as a βα row and whenever nucleic acid
is drawn. It suits small α/β folds; on others it can add crossings, so it is opt-in.

```bash
foldmap plot 2ACY --set helix_on_sheet=auto -o acylphosphatase.svg
```

**Domain colours and the domain bar.** When domains are named (`--domain`, `--domains auto` or `--domains uniprot`),
elements are coloured by domain in eight colour-blind-safe tones. This replaces the default chain colouring
(`domain_colours`, on by default); a theme with its own colouring (rainbow, conservation, ...) keeps it, and
`--set color_by=domain` asks for the domain tones anyway. A linear bar under the figure shows each chain's domains in
residue order, on one residue scale for every chain (`domain_bar`: `auto` draws it for 2 or more domains, `on` for
any, `off` never). The print and blueprint themes keep their own palettes: greys with a hatch per domain, or the
blueprint tints. `--set domain_colours=true` asks for the tones there. A domain's name sits at the top left of its
panel, or at the top right when the left is taken. Each domain is kept together, but when domains interleave in 3D
(Cas9's RuvC parts, measles F's DI–DIII) their panels can overlap. `docs/domains/` holds the `--domain` flags for
Cas9 (4OO8, Nishimasu 2014 Fig 1A) and measles F (8UTF, Zyla 2024 Fig 3A):

```bash
foldmap plot 4OO8 --assembly 1 $(grep -v '^#' docs/domains/4OO8-cas9.txt) -o cas9.svg
foldmap plot 8UTF --focus protomer $(grep -v '^#' docs/domains/8UTF-measles-F.txt) -o measles-F.svg
```

**Chromophores and other residues fused into the chain.** A peptide-linked non-standard residue without the usual
backbone atom names (GFP's CRO66) keeps the chain whole: no dashed "missing" loop, and the helix through it stays one
element. It is marked with a star on the chain and the legend entry "chromophore / linked residue".

**Mirrored dimers.** Two copies related by a two-fold axis lying in the page (λ repressor on its operator, 1LMB) are
drawn as exact mirror images about that axis, names included, even where that takes up to a tenth more room than placing
the copies one by one (other copies: a twentieth). By default (no `--up`, `--view`, `--rotate` or `--flip-v` asked for,
and no straight duplex to draw along) a two-fold (C2) assembly is seen from the side of its dyad for this: the dyad
points up the page and the first protomer sits on the left, so the copies meet at the line where their interface lies,
when each copy is one protein chain (nucleic-acid chains do not count). A copy of two or more chains is unrolled side by
side instead: interfaces between different chains cross the dyad, and no mirror line keeps them together (haemoglobin,
2HHB, whose α1β2 contact would put α1 and β2 in opposite corners).
Of the views around the dyad (every 15°, up to 75° either way from square-on to the line between the protomers'
centres), foldmap takes the one whose projection best keeps the order of the distances between elements (their rank
correlation with the 3D distances), the smallest turn of those within 0.01 of the best, so that elements close in 3D
across the two copies stay close on the page (TrmD, 1UAM, is seen turned by 75°). Where no mirror image can be drawn (a
βα row such as a TIM barrel's would cross the line), or where the mirror image draws the copies unalike or leaves a
membrane dimer's transmembrane helices off their band and unrolling does better, the copies are unrolled side by side
instead. `--set two_fold=unrolled` (the app's "Two-fold assemblies" choice, offered when one is drawn: automatic,
mirror image, side by side) always draws a two-fold assembly that way, its copies side by side around the axis as any
other assembly's; `two_fold=mirror` seeks the mirror image for any two-fold assembly, whatever its copies hold (with the
same fallbacks); `two_fold=auto`, the default, chooses by the chains of a copy as above. Other symmetric copies,
and the unrolled cylinders of rings and filaments, are drawn alike.

A ring of copies (a Cn assembly) is unrolled side by side by default. For a four- or seven-fold ring (C4, C7) drawn
whole, `--set cyclic=rosette` (the app's "Cyclic assemblies" choice, offered only then: side by side, down the axis)
draws it seen down its n-fold axis instead, marked by a filled n-gon at the centre: the first copy is laid out in its
own frame, and every other copy is that drawing, loops included, turned by a whole multiple of 360/n degrees about
the axis, all pushed out from it just enough to clear their neighbours. Names stay upright, each on the side of its
helix that is clear. Each copy's sheets are its own: a pairing between two copies' strands is not drawn. Ligands held
by every copy (the ions of a channel's pore) are drawn once, as one untethered cluster at the centre. A ring with DNA
or RNA, a chain outside the copies, a stack row, a membrane view or domains (`--domains`, whose panels do not turn
with a copy) stays unrolled, and the notes say why. Where the figure draws one subunit of a large assembly (`--focus
auto` or `protomer`) or only the chains named (`--focus A,B`), no rosette is drawn and a note says it needs the whole
assembly (`--focus none`). Edits (move, swap, rename) name the first copy's elements, which every copy follows; an edit
naming another copy's element is not applied, and a note says so. Other orders are not offered: in a look test on 56
assemblies, three-, five- and six-fold rosettes, and D2 assemblies drawn as two mirrored rings, read worse than their
unrolled figures on most entries, with more check failures (route fallbacks, text collisions, sheet order); no
eight-fold assembly was at hand. With `cyclic=rosette` those are drawn as the default figure, with a note saying so.

**Layout files.** `--save-layout fig.yaml` records everything that defines a figure: theme, style changes,
layout options, domains and edits, with element references stored as residue ranges so they survive a new model.
`foldmap plot model.cif --layout-file fig.yaml -o fig.svg` re-renders it; the command line overrides the file.

### RNA and DNA

DNA and RNA are drawn as their fold. Base pairs are found from the coordinates (Watson–Crick and G·U, with
modified nucleotides paired by their parent base) and grouped into stems, named `P1`, `P2`, ... for RNA and `DNA 1`,
... for DNA. Each stem is placed as a rigid item beside the protein that holds it; a duplex that bends (nucleosome, IHF, TBP)
is cut where its axis turns by 25° into straight pieces placed end to end. Single strands, hairpin loops (as arcs)
and the links between stems are routed like protein loops, with both ends marked 5′ and 3′, base letters, and contact
beads in the colour of the protein chain that touches them. A nucleic acid with no pairs at all (a G-quadruplex, one
strand on a protein) is drawn as a projected trace with letters. `foldmap summary` lists each stem with its base
pairs and ranges, e.g. `P1  7 bp  1:1-7 · 1:112-118`.

RNA stems are numbered by their nesting and lettered across internal loops (`P2a`, `P2b`); a stem that crosses a
longer one is a pseudoknot, named `PK` (`PK1`, `PK2`, ... when there are several). A tertiary base pair outside the
stems is marked at both ends by a tie letter (`a`…`a`). A chain break inside one molecule is drawn as a dashed gap
labelled with the number of nucleotides left out, not as new 5′ and 3′ ends.

```bash
foldmap plot 1EHZ --names docs/names/1EHZ-trna.yaml -o trna.svg      # tRNA arms named acceptor, D, anticodon, T
foldmap plot 30QM --set na_style=block -o 30qm.svg                   # stems as bars
foldmap plot 1LMB --set hbonds=na --set hbond_max=16 -o 1lmb.svg     # protein–DNA hydrogen bonds
```

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
`na_line` sets the backbone's weight (1 = the default; `1.5` draws runs, ladder rails, ribbon tubes and block edges
half as thick again). It is a layout key: base letters and labels keep off a thicker strand, so they may move.
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

**Interactions.** `--set interactions=hbond,salt,stack,cation-pi,hydrophobic` (any of these; `all` is every kind but
hydrophobic, `all,hydrophobic` adds it; default `off`) chooses the kinds of interaction drawn, and
`--set interactions_between=auto|na|interchain|chains:A+B,C+D` their partners (`auto`, the default: `na` when the
figure holds nucleic acid, else `interchain`). `interaction_max` caps the lines of each kind (default 200; H-bonds keep
`hbond_max`). `hbonds=na` keeps working and means `interactions=hbond` with `interactions_between=na`, so old layout
files and app sessions draw the same figure; when both keys are given, `interactions` wins and a warning says so.
In the app's Show section one set of controls covers both: a checkbox per kind, the partners, and the caps (an old
file with `hbonds: na` opens as H-bonds ticked between protein and DNA/RNA, and is saved with the interaction keys
once those controls change). The MCP tool `draw_topology` takes the keys in `style` and lists every line drawn in its
structured result (`interactions`: kind, the two residues, the shortest distance, the atoms and the line's SVG id).

Each kind is a thin line between the two residues, as the H-bond lines are, one line per kind and residue pair (a
pair's several salt-bridge atom pairs, or a Trp's two stacked rings, make one line), with its own colour and dash:
H-bonds magenta short dashes, salt bridges crimson long dashes, π-stacking dark green dense dots, cation-π blue
dash-dot, hydrophobic contacts dark olive sparse dots. The colours stand clear of every chain palette, the H-bond
magenta and the DNA/RNA backbones in normal and colour-blind sight; the dashes alone tell the kinds apart in print
(the `print` and `blueprint` themes draw every kind in `hbond_color`, by dash only). Two kinds between one pair lie
side by side. The legend names each kind drawn with its count ("Salt bridges (200 of 1,412 shown)" when capped).
Criteria, from the heavy atoms: salt bridge, an anionic O (Asp, Glu, C-terminus, phosphate) within 4.0 Å of a
cationic N (Lys, Arg, N-terminus; His "if protonated"); π-stacking, ring centroids within 5.5 Å, parallel (≤ 30°)
or T-shaped (60–90°), offset ≤ 2.0 Å, never two nucleobases (the ladder shows those); cation-π, Lys NZ or the Arg
guanidinium centre within 6.0 Å of a ring centroid, offset ≤ 2.0 Å; hydrophobic, carbons bonded only to carbon
within 4.0 Å, never two nucleotides (across a duplex those are pyrimidine carbons of neighbouring base pairs, which
the ladder shows). Before drawing, every interaction found is re-measured by an independent check
(`interactions.verify_interactions`); one it flags is not drawn and the figure's notes say so. Example:
`foldmap plot 1LMB.cif --set interactions=all -o 1lmb.svg`.

**Checking the science.** A figure can look clean and still say something false about the structure, so foldmap
checks what it draws against the coordinates, apart from how the page looks. `foldmap check 9FWZ` draws the figure
`plot` would (same options, `-o` to save it) and prints one line per check and one per problem, exiting 1 on any:
S1, two strands drawn side by side in a sheet (adjacent columns, level with each other) are joined by at least one
backbone bridge, and two strands in one column (collinear, end to end) never overlap; S2, a closed barrel's closure
(the ghost beside the last strand) is joined; S3, a barrel's strand count is the cycle its pairing makes, and no sheet
drawn open holds a ring of 5 strands or more; S6, a strand off a sheet's run of paired strands stands beside its
partner, level with it, or the notes name the pairing ("x:67-68 also pairs with x:79-86"). S4, a linked residue's star (a chromophore, say) stands at its residue on the
drawn figure; it judges only linked-residue stars, not ligand tethers, glycans, disulfide ends or UniProt site marks,
and `foldmap summary` (which draws no marks) reports it not measured. S5 (automatic domains against a reference) is
measured only in the campaign. `foldmap summary` prints the same checks of the sheets as built.

Under them `foldmap check` prints the figure criteria A1-A20, the page checks the random-100 campaign judges every
figure by, with the same limits (`src/foldmap/criteria.py`): no route fallbacks (A3) or loops over labels (A4), no
overlapping elements (A5), no text collisions (A6) or text off the page (A7), labels unique per chain (A8), nothing
left out without a note (A9), N→C and sheet order (A11), copies drawn alike (A15), domains with `--domains` (A16),
the membrane band with `--membrane` (A18) and local ligand tethers (A20); each line reads pass or fail with what
failed, "not judged" with why when the figure gives it nothing to judge (no copies, no tether to another chain), or
"does not apply". The science checks read "not judged" likewise when nothing is drawn for them (no sheet, no barrel,
no residue mark, no accessory strand). A `--membrane` or `--domains` figure is judged as the campaign judges its membrane or
domains view: the criteria it keeps for the default figure (A8, A9 on a membrane figure, A11, A15, A20) read "not
measured here", to be checked on the figure drawn without that option. The criteria that need more than this one drawing or a tool foldmap does not use (the
process's exit and time, the app's hosted job, mkdssp, a second hash seed, a saved-layout round trip, RNApolis base
pairs, a two-state partner) read "not measured here" with the reason, as do the text checks when the text is not
set in the font they are calibrated on (Liberation Sans; A16 then still judges its domain rules and panel overlaps,
only the domain names go unmeasured); they are never counted as passing. A failed criterion
also exits 1.

A figure drawn like a reference (`foldmap check 1MBN --like 2HHB --like-chain B`) is also judged on H1-H5, after
A1-A20: H1, it is drawn alike exactly when the two are one fold (TM-score at least 0.5 and half the smaller chain's
elements matched), else refused with the note; H2, matched elements sit in the same places (the median and 90th
percentile of the mates' centre deviations after the best rigid 2D motion, at most 1.0 and 2.5 page units); H3, the
map reads the same way (of every two mates more than 1 unit apart, at least 95 % keep the signs of their page dx and
dy, and each mated sheet keeps its strands' order unless the note names a pairing difference); H4, mates carry one
name (A19's rule); H5, the figure still passes every S- and A-check the structure drawn alone passes. H2 and H3 also
print the same measures for the two drawn independently, today's behaviour, for comparison. A check with nothing to
judge (fewer than three mates, a pair refused) reads "not judged". The H1-H3 limits and the pull that holds a
matched element at its reference place were set on a tuning set of 15 homologue pairs
(`docs/validation/homologues/tuning.json`, `scripts/homologues.py`) and frozen before any held-out pair was drawn.
Measured once on 40 held-out pairs drawn from CATH by a seeded script (`docs/validation/homologues/heldout.json`:
20 of one topology at under 30 % identity or across superfamilies, 10 within a superfamily at 30-90 %, 10 negatives
of one architecture), 16 of the 30 positives pass H1-H5 and 9 of the 10 negatives are refused (the tenth, whose
reference is drawn as a DNA duplex, is not judged). Seven cross-superfamily pairs score under TM 0.5 and are refused;
one close pair has a sheet carried mirrored, and others differ in strand pairing (`heldout-report.md`, `review.md`).
Over the positives drawn alike, the mates' median deviation from the reference is 0.11 page units, against 1.89 for
the two drawn independently.

**Checking the contact beads.** `foldmap contacts 1LMB.cif` draws the figure `plot` would (it takes the same options),
checks every contact bead against the deposited atoms with an independent search, and prints how many beads were
drawn, verified, false, missing, unmatched, in the wrong colour or not on the first toucher, plus a line for each
problem bead. It exits 0 when clean and 1 on any false, unmatched or wrong-colour bead; `-o fig.svg` also saves the
figure.

**Validation.** Thirty-nine DNA and RNA entries (tRNAs, riboswitches, ribozymes, CRISPR complexes, RNPs,
polymerases, a ribosome, bent protein–DNA complexes, junctions and non-canonical structures) are drawn and scored
by `scripts/na_validate.py`: base pairs and stems against RNApolis, the drawn layout against the 3D structure,
contact beads against the atoms, and overlaps, crossings and text on the page. 34 of 39 pass every limit; each of
the five that fail (30QM, 4PJO, 3KTQ, 1YTB, 3CRX) is explained with its cause in
[docs/validation/na/review.md](docs/validation/na/review.md). The full table is in
[docs/validation/na/report.md](docs/validation/na/report.md).

**Random-100.** One hundred PDB entries and AlphaFold models drawn at random in twelve strata (X-ray, cryo-EM and
NMR; oligomers and hetero-complexes; membrane proteins; protein–DNA, protein–RNA and RNA alone; glycoproteins; one
protein in two states; AlphaFold models; up to ribosomes near the app's 30,000-residue cap) are drawn by
`scripts/random100.py` in every view and judged on the twenty criteria A1–A20 that `foldmap check` prints (A1c, a
CA-only entry refused in one line, met no such entry), from crashes, time, route fallbacks and text collisions to
secondary structure against mkdssp (268,834 of 268,835 residues agree), determinism and copies drawn alike. Each
criterion reads pass, fail, not judged (nothing to judge in that entry: no copies, no tether to another chain, no
domain drawn) or does not apply, and an entry passes when none fails; nothing to judge is never counted as a pass
(A15, copies drawn alike, judged 37 entries and passed all 37). The layout's constants were tuned on this first
sample, and a second sample (new seed, same strata, none of the first's entries), held out at first, became the input
of the plan that fixed its failures, so a third sample was drawn after those fixes, with none of the other two's
entries, and run once with no fix. **On structures foldmap was not tuned on (the third sample), 83 of 100 pass**,
in the same run as 93 of 100 on the tuned sample and 95 of 100 on the second sample (80 of 100 before it was
fixed against). The third sample's failures are listed by check in
[docs/validation/random100-third/review.md](docs/validation/random100-third/review.md); the tuned sample's seven
(the ribosomes 7Y7C and 8VWQ, 9YLM and 8BW5 as ruled, the two-state pair 9YKJ/10UJ named two ways, and 8QEM's DNA
shape in its whole-ring view) are explained in
[docs/validation/random100/review.md](docs/validation/random100/review.md), with the full table in
[report.md](docs/validation/random100/report.md); the second sample's in
[docs/validation/random100-heldout/review.md](docs/validation/random100-heldout/review.md).

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

Helices, strands and 3₁₀ helices are drawn over their residues as the figure draws them: the helix ribbon (one
turn per 3.6 residues), the arrow with its head on the last residue and the rounded 3₁₀ box, with the figure's
`fill` and `helix_shading`.

## How the topology is made

Every figure (command line, app, MCP or Python) goes through the same steps:

1. **Read the structure.** gemmi reads the mmCIF/PDB file (or the file downloaded for a PDB ID or UniProt
   accession). With `--assembly auto` the biological assembly is built from the file's own operators.
2. **Secondary structure.** A DSSP implementation assigns hydrogen bonds from the backbone and calls α-, 3₁₀- and
   π-helices, β-strands (with bulges) and turns, following mkdssp 3.1.2. Bends (S) are not assigned; with mkdssp's
   bends read as none, 268,834 of 268,835 residues agree over the random-100 campaign's 96 protein entries (it
   requires 99% per entry).
3. **Sheets and bundles.** Strands joined by backbone H-bond ladders form sheets; their order and directions come
   from the ladder network, and a sheet whose ends pair is a closed barrel. Each strand drawn side by side with
   another pairs with it: a branched sheet's run of paired strands is its ladder network's longest path (a closed
   barrel's: its longest pairing cycle, so FimD 9FWZ has 24 barrel strands and OmpF 16), and a strand off that run
   (a short strand pairing with one barrel strand, or a branch of several) is drawn beside its partner on the side it
   lies on in 3D, level with it: in the next strand's column, end to end with that strand (one column means
   collinear), or in a new column past the sheet's end; never between two strands that pair. When no height level with
   its partner is free there, it stands as near as it can and the note names the pairing ("also pairs with"); so do
   pairings a ring cannot draw side by side (its chords). The figure's note names each such strand with its partner
   ("D:662-666 drawn beside its partner D:645-659"); `foldmap summary` words each as the default figure draws it,
   level ("drawn beside its partner") or not ("A1:2-6 pairs with A1:9-23 (off the run)"). Barrels that pair with one another (the
   four OB folds of the SSB tetramer 1EYG) are each drawn closed, as a sheet of its own; the note and the summary name
   the pairings joining them. Helices packed along their length across chains
   form bundles.
4. **Context.** Symmetry (cyclic, dihedral, cubic, helical) is found by superposing chains that are copies of one another (at least 90% of the shorter sequence matches the
   other by difflib's matching blocks, so point mutants and different constructs count);
   disulfides, glycans and ligands are read from the atoms; a membrane is estimated on request (the normal is the dominant axis of the long helices and strands, or a barrel's
   axis; a fixed-thickness slab is slid along it to the most hydrophobic, least charged position; no tilt search). Nucleic acids:
   base pairs grouped into stems; long bent duplexes cut where the axis turns by 25°.
5. **Projection.** A frame turns 3D into page coordinates: the dominant element direction points up and the
   widest spread runs left–right (single chains); a two-fold assembly of one-chain copies is seen from the side of its
   dyad, the dyad up the page (its copies drawn as mirror images, step 6); any other symmetric assembly is unrolled around its axis,
   protein–DNA complexes around the DNA axis, membrane proteins with the bilayer horizontal. The frame itself is
   never a mirror. An assembly whose nucleic acid unrolling would tear apart stays in the plain projected view: a duplex kinked across protomers,
   linked stems that would turn against one another, or most stems each held by two protomers (the PFV intasome
   3L2V, whose viral DNA ends are gripped by both inner integrase chains).
6. **Layout.** Sheets (as rigid blocks), helices and DNA are placed by minimising one energy: springs to their
   projected positions, short loops between sequence neighbours, attraction between elements in contact, a pull
   between disulfide-bonded cysteines, and an overlap penalty; then any remaining overlaps are cleared.
   Symmetry copies are tied to look identical (two copies across a dyad in the page as mirror images); an
   all-parallel α/β sheet is placed as one βα row; stems are placed as rigid items. The optimiser is
   deterministic, so a figure comes out the same on every machine where the same font is installed.
7. **Loops.** Each loop is routed on a grid with A* search: orthogonal lines that go around elements and labels;
   single strands are routed the same way, hairpin loops as arcs.
8. **Drawing.** matplotlib draws the elements to scale (one helix turn tiled along each helix, arrows for
   strands), then loops, annotations, legend and title, and writes SVG (text kept as text), PNG or PDF.

Steps 1–7 depend only on the structure and the layout options; themes and style keys only change step 8, so a
figure keeps its shape while you restyle it.

## App

`foldmap app` opens a Streamlit app in the browser. It runs the same pipeline as the command line, so a figure
from the app is the figure `foldmap plot` draws with the same settings. It is also online at
<https://foldmap.streamlit.app>.

**Start.** One box takes a PDB ID or a UniProt accession (and, on your own machine, a file path); press Enter or
Draw. An accession lists its AlphaFold model and the experimental structures (widest coverage first); pick one to
draw it. The examples are one click, and a `.cif` or `.pdb` file can be dropped in.

**Workspace.** The figure is drawn at once in the current theme. The top bar names the structure, summarises it,
and switches between three views: Figure, Sequence (the sequence, or an uploaded alignment, with the secondary
structure on top) and 3D (the interactive page: click an element of the figure, or a cell of the contact map, and the 3D view glides to it). Only the chosen view is drawn. Download makes each file when you
click it, from what was drawn: SVG, PNG (300 dpi, less past 50 megapixels; the button says which), PDF and the
layout file for the figure, the sequence figure, and the 3D page as HTML. A setting that fails keeps the last good figure on screen and names the control to check.
**New structure** returns to Start; the theme and style stay for the next structure, while the last structure's
chains, edits, marks and alignment go.

The sidebar starts with the theme (Browse themes shows them all), the colouring, labels, legend, text size and
title. The rest is in sections:

- **Style**: fill, helix shading, loops, background and ink.
- **Show**: ligands, disulfides, glycans, membrane, DNA and RNA, interactions (H-bonds, salt bridges, π-stacking, cation-π, hydrophobic; their partners and caps), and the contacts check.
- **Chains & assembly**: chains to draw, highlight, assembly, symmetry, and how large assemblies are drawn.
- **Edit layout**: arrangement, rotation, sizes, the βα template, helices packed on a sheet, swap, move, rename
  and mark (as tables or as text), element names, and domains with their colours and bar.
- **Conservation & UniProt**: an alignment for conservation colours, and UniProt sites and domains.
- **Layout file**: load one saved from Download (or by `foldmap plot --save-layout`) to restore every setting.

**Part of the file.** When foldmap draws part of the file on its own (assembly 1 when the file holds several
copies of a complex, or one subunit of a cage, a filament or a large ring), a bar above the figure says what was
drawn and how many of the file's chains, as in "Drawn: one of 24 subunits of an O cage — 1 of 24 chains."
**Show all chains** draws every chain; **Keep this** hides the bar for that structure. When every chain is too large for the online app, the bar gives the command that draws it instead.

**Share links.** A structure loaded by PDB ID or accession puts its id and theme in the address, as in
`https://foldmap.streamlit.app/?id=1LMB&theme=journal`, so the link opens the same structure in the same theme.
Uploads and file paths never enter the link.

**Slow structures.** A figure that takes longer than 4 s to draw turns automatic redrawing off for that
structure: change what you want, then press Redraw (or switch Auto-redraw back on).

### Hosting the app

The live app runs on Streamlit Community Cloud from this repository, with `src/foldmap/app.py` as the main file
and Python 3.11 or newer. Cloud runs it from the repository root and reads these files there:

- `requirements.txt` installs foldmap from the repository, editable (`-e .[app]`), with its `app` extra:
  Streamlit 1.65 or newer, which the app needs (it sets the sidebar's width in pixels). Editable, because Cloud
  pulls each push into its checkout without reinstalling: the app then runs the code it was pushed with. Cloud
  reads this file before `pyproject.toml`, which it would give to Poetry. Poetry would not install the `app` extra,
  and it refuses a version read from the code (`foldmap/__init__.py`).
- `.streamlit/config.toml` holds the app's theme, the viewer toolbar (no Deploy menu), 50 MB uploads and no error
  details for visitors. It is generated from `foldmap.app.APP_CONFIG` (the command is its first line), and
  `streamlit run` reads it when it starts from the repository root. `foldmap app` passes the same values as flags.
- `packages.txt` installs `fonts-liberation`: Liberation Sans, the font the figures and the layout's text
  metrics are made for.

On Community Cloud (the app runs from `/mount/src`) the app is hosted unless told otherwise; elsewhere set
`FOLDMAP_APP_LOCAL_FILES = "0"` (on Cloud, as a root-level secret in the app's settings: Cloud passes root-level
secrets to the app as environment variables). Only `FOLDMAP_APP_LOCAL_FILES = "1"` makes a Cloud app open server
paths. The server's log says which mode it runs in. Hosted:

- the app never opens a file path on the server, including alignment and names paths inside uploaded layout
  files, and files it offers for download name uploads, never server paths;
- uploads are kept in `~/.cache/foldmap/uploads` for a day;
- every drawing runs in its own process, within limits.

| Variable | Hosted default | |
| --- | --- | --- |
| `FOLDMAP_APP_TIMEOUT` | 90 | seconds a drawing may take |
| `FOLDMAP_APP_MEMORY_MB` | 1024 | memory of one drawing, MB (two drawings and the app fit Community Cloud's memory) |
| `FOLDMAP_APP_MAX_RESIDUES` | 30,000 | larger structures are refused before any wait, with the command to draw them locally |
| `FOLDMAP_APP_JOBS` | 2 | drawings at a time |
| `FOLDMAP_APP_SLOW` | 4 | seconds; a figure slower than this turns automatic redrawing off |

A structure that went over a limit is remembered, so the next visitor who asks for it gets the message at once.
On your own machine there are no limits except `FOLDMAP_APP_SLOW`; each variable sets its limit there too.

## MCP server for AI agents

`foldmap mcp` runs a [Model Context Protocol](https://modelcontextprotocol.io) server over stdio, so agents can
summarise structures and draw figures. Tools: `summarize_structure`, `draw_topology`, `compare_structures`,
`draw_sequence`, `interactive_page`, `list_styles`; the drawing tools return a PNG preview along with the file they
wrote.

```bash
claude mcp add foldmap -- foldmap mcp        # Claude Code
```

```json
{ "mcpServers": { "foldmap": { "command": "foldmap", "args": ["mcp"] } } }
```

(the second form for clients configured by JSON, such as Claude Desktop).

Figures must end in `.svg`, `.png` or `.pdf`, pages in `.html`, layout files in `.yaml` or `.yml` and
`compare_structures`' record in `.json`. Set
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
.venv/bin/python -m pytest              # about 2200 tests, ~5 min
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/python scripts/gallery.py     # re-render the images above
.venv/bin/python scripts/sweep.py --sample 60 --out sweep/   # render random PDB entries and check them
.venv/bin/python scripts/bench.py                                # time layouts (add --record/--check to compare runs)
```

**The website** (<https://dzyla.github.io/foldmap/>) is built by `scripts/build_site.py`, which renders every figure
with foldmap itself, so the site shows what the current code draws. The `website` workflow
(`.github/workflows/pages.yml`) rebuilds and publishes it on every push to `main`. To check it before a push, build it
locally and open `_site/index.html`:

```bash
.venv/bin/python scripts/build_site.py --out _site            # the full site
.venv/bin/python scripts/build_site.py --out _site --quick    # a small subset
```

The README images come from `scripts/gallery.py` (into `docs/images/`).

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
