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
  images, conventional names (`--names`) and domain colours with a domain bar.
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
<sub><b>βα rows.</b> Triosephosphate isomerase (1TIM, chain A), rainbow theme, <code>--names beta</code>: the closed (βα)₈ barrel as one row, each helix between the strands it joins, opened between β8 and β1, β8 repeated as a ghost.</sub></td>
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
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/ncrna-30qm.png" alt="DRT10 RT-ncRNA complex with six RNA stems drawn as ladders under the two protomers"><br>
<sub><b>RNA–protein complexes.</b> DRT10 RT–ncRNA (30QM): the 121-nt ncRNA's six stems (P1–P6) as ladders under the two RT protomers, contact beads in the colour of the chain that touches them.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/lambda-hbonds.png" alt="Lambda repressor dimer on DNA with protein-DNA hydrogen bonds as dashed lines"><br>
<sub><b>Hydrogen bonds.</b> λ repressor (1LMB), <code>--set hbonds=na --set hbond_max=16</code>: the 16 shortest of 30 protein–DNA hydrogen bonds as dashed lines.</sub></td>
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
foldmap contacts 1LMB                             # check the contact beads against the atoms
foldmap app                                       # the Streamlit app
foldmap styles                                    # every theme and style key
```

## What it handles

- **Secondary structure** by a DSSP implementation that matches mkdssp 3.1.2 residue for residue on the test set:
  α-, 3₁₀- and π-helices (π before α, as in DSSP ≥ 2.1), β-strands with bulges, turns; β-sheets
  with strand order and direction, closed barrels shown with a ghost of the first strand. A helix kinked by more
  than 30° is split into pieces drawn end to end; a strand that a bulge splits stays one strand when both pieces
  pair with the same partner.
- **Layout** from the 3D structure: elements projected onto the page, then arranged so sheets stay in order,
  contacting elements and disulfide partners stay close, and nothing overlaps; loops routed orthogonally around
  elements and labels. The layout is deterministic: the same figure on every machine.
- **Fold conventions** of the published figures: all-parallel α/β sheets (TIM barrels, Rossmann folds) as one βα
  row, conventional names (`--names`), helices drawn over the sheet they pack on (opt-in); see
  [Conventions of the published figures](#conventions-of-the-published-figures).
- **Assemblies and symmetry**: `--assembly auto` builds the file's first biological assembly when it adds copies
  (up to 24 chains) or when the file holds several copies of one complex, and never one that leaves out the DNA or
  RNA; cyclic, dihedral, cubic and helical symmetry is detected and the protomers are placed alike (two copies across a
  dyad in the page as mirror images; large, interleaved protomers such as the measles F trimer can still differ); highlight one protomer or the asymmetric unit; filaments, cages and large
  assemblies drawn as one subunit with the neighbouring parts that complete its fold. When either automatic choice
  draws part of the file, foldmap says so (a note on stderr, a bar in the app, a line in the MCP result) and names
  the setting that draws everything (`--assembly asu`, `--focus none`).
- **DNA and RNA as their fold**: base pairs from the coordinates, stems, hairpins, pseudoknots, junctions and single
  strands placed beside the protein that holds them; bent duplexes as chains of straight pieces; protein–DNA/RNA
  hydrogen bonds on request; see [RNA and DNA](#rna-and-dna).
- **Annotations**: disulfides, glycans (SNFG symbols, sialic acids included), ligands and metal ions, chromophores
  and other residues fused into the chain (starred), membrane bands, named domain panels with a domain bar, residue
  numbers, marked elements.
- **Colouring**: by chain, N→C, per-element shade, secondary-structure type, domain, B-factor, hydropathy,
  AlphaFold pLDDT or conservation from an alignment.

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
`ligands` (`auto`, `all`, `none` or codes like `HEM,ZN`), `fold_template`, `helix_on_sheet`, `domain_bar`,
`domain_colours`, `na_backbone`, `backbone_colors`, `dna_color`, `rna_color`, `hbonds`, `hbond_color`, `hbond_max`,
`disulfides`, `glycans`, `helices_310`, `na_extent`, `na_style`, `background`, `ink`, `helix_scale`,
`strand_scale`, `font_scale` (`dna_extent` still works on the command line and in style and layout files).

## Shaping a figure

| Option | Example |
| --- | --- |
| Swap two elements | `--swap α1,α3` |
| Nudge an element | `--move α2=1,-2` |
| Rename an element | `--rename res:A:167-182=Gd` |
| Conventional names | `--names deposited` (the file's helix ids, e.g. globin A–H), `--names beta` (β1…βn) or `--names docs/names/2HHB-globin.yaml` (repeatable) |
| Highlight elements | `--set mark=α4=#2ca02c` |
| Named domain panels | `--domain N-lobe=res:A:13-59` (repeatable), `--domains auto` or `--domains uniprot` |
| Domain colours and bar | `--set color_by=domain`, `--set domain_bar=on\|off`, `--set domain_colours=false` |
| βα row | `--set fold_template=off` draws an all-parallel α/β sheet as projected |
| Helices over a sheet | `--set helix_on_sheet=auto` (opt-in) |
| Hydrogen bonds | `--set hbonds=na`, `interchain` or `chains:A+B`; `--set hbond_max=30` |
| Only some chains | `--chains A,B` |
| Large assemblies | `--focus auto\|none\|protomer\|B` (one subunit plus what its neighbours add to its fold) |
| Membrane | `--membrane auto` (OPM dummy atoms if present, otherwise estimated); `--membrane flip` puts the other side on top |
| Conservation | `--theme conservation --msa family.fasta`, or `--msa consurf_grades.txt` for ConSurf's own grades |
| Orientation | `--up X,Y,Z`, `--view X,Y,Z`, `--rotate DEG` |
| Symmetry and assembly | `--symmetry auto\|off\|C3\|D2\|helical`, `--protomers "A,B;C,D"`, `--assembly auto\|asu\|ID` |

Elements are named by label (`A`, `α2`, `η1`), `chain:label`, position (`#3`) or residue range (`res:A:10-25`). Negative numbers and insertion codes are written as in the file: `res:A:-7-4`, `res:H:100A-100K`.

A label reference reads the labels the figure shows: the label copies of one protein share, the primed labels of
elements a neighbour lends in focus mode (`α2′`) and any `--names`; `--rename` is applied last and wins over all of them.

**Domains and βα rows.** `--domains auto` never splits a sheet drawn as a βα row (a TIM barrel, a Rossmann fold):
the row and its helices fall in one domain. A `--domain` set of your own that splits a row is drawn as given, but the
row is drawn as one unit, so its domains' panels overlap it. To let the parts move apart, draw the sheet as
projected with `--set fold_template=off`.

### Conventions of the published figures

foldmap follows the conventions of the hand-drawn topology figures readers know (Richardson 1981, the TIM and
Rossmann schematics, Burroughs 2007, Nishimasu 2014). Each is a layout or style key; the defaults are below.

**Conventional names** (`--names`, opt-in). `--names deposited` takes the helix ids the file itself gives, so a
globin's helices read A–H (an element the file gives no id keeps its generic label). `--names beta` numbers the
strands β1…βn in sequence and names each helix of a βα row after the strand before it (β1 α1 β2 α2 …).
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
drawn as exact mirror images about that axis, names included. Other symmetric copies, and the unrolled cylinders of
filaments, are drawn alike.

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

**Validation.** Thirty-nine DNA and RNA entries (tRNAs, riboswitches, ribozymes, CRISPR complexes, RNPs,
polymerases, a ribosome, bent protein–DNA complexes, junctions and non-canonical structures) are drawn and scored
by `scripts/na_validate.py`: base pairs and stems against RNApolis, the drawn layout against the 3D structure,
contact beads against the atoms, and overlaps, crossings and text on the page. 33 of 39 pass every limit; each of
the six that fail (30QM, 4PJO, 3KTQ, 4YBB, 1YTB, 3CRX) is explained with its cause in
[docs/validation/na/review.md](docs/validation/na/review.md). The full table is in
[docs/validation/na/report.md](docs/validation/na/report.md).

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
   Symmetry copies are tied to look identical (two copies across a dyad in the page as mirror images); an
   all-parallel α/β sheet is placed as one βα row; stems are placed as rigid items. The optimiser is
   deterministic, so a figure comes out the same on every machine.
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
structure on top) and 3D (the interactive page). Only the chosen view is drawn. Download makes each file when you
click it, from what was drawn: SVG, PNG (300 dpi), PDF and the layout file for the figure, the sequence figure,
and the 3D page as HTML. A setting that fails keeps the last good figure on screen and names the control to check.

The sidebar starts with the theme (Browse themes shows them all), the colouring, labels, legend, text size and
title. The rest is in sections:

- **Style**: fill, helix shading, loops, background and ink.
- **Show**: ligands, disulfides, glycans, membrane, DNA and RNA, hydrogen bonds, and the contacts check.
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
