# foldmap

**Publication-grade protein topology diagrams from PDB/mmCIF structures.**

foldmap reads a structure, assigns secondary structure, finds β-sheets, helix bundles and symmetry, and lays the
elements out as a clean 2D topology diagram. Helices are drawn as coiled ribbons, strands as arrows and short 3₁₀
helices as small boxes, all on one consistent scale, with loops routed around everything else. It handles
protein–DNA/RNA complexes, symmetric assemblies from cryo-EM and crystallography, AlphaFold models, disulfides and
glycans.

<table>
<tr>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/lambda-repressor-dna.png" alt="λ repressor dimer bound to DNA"><br>
<sub><b>Protein–DNA.</b> λ repressor (1LMB): the duplex is unrolled above its binding helices, with base letters and contact beads.</sub></td>
<td width="50%"><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/t4-lysozyme-domains.png" alt="T4 lysozyme with N- and C-lobe domain panels"><br>
<sub><b>Domain panels.</b> T4 lysozyme (2LZM), journal theme: named domains kept together on their own panels.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/fima-richardson.png" alt="FimA pilin coloured by secondary structure type"><br>
<sub><b>Disulfides and 3₁₀ helices.</b> FimA (5NKT), Richardson colouring: the disulfide is drawn as a bar, the 3₁₀ helix as η1.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/p53-alphafold.png" alt="p53 AlphaFold model coloured by pLDDT confidence"><br>
<sub><b>AlphaFold models.</b> p53 (AF-P04637), alphafold theme: elements, loops and tails in pLDDT confidence bands.</sub></td>
</tr>
<tr>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/measles-f-trimer.png" alt="Measles fusion protein trimer with one protomer highlighted"><br>
<sub><b>Symmetry.</b> Measles F trimer (8UTF): C3 detected, protomers drawn alike, one highlighted and its copies greyed.</sub></td>
<td><img src="https://raw.githubusercontent.com/dzyla/foldmap/main/docs/images/tim-barrel-rainbow.png" alt="Triosephosphate isomerase dimer coloured N to C"><br>
<sub><b>N→C colouring.</b> Triosephosphate isomerase (1TIM), rainbow theme: each element shaded by its place in the chain.</sub></td>
</tr>
</table>

## Install

```bash
pip install .            # command-line tool
pip install ".[app]"     # plus the Streamlit app
```

Requires Python 3.11 or newer.

## Quick start

```bash
foldmap summary model.cif                     # elements, sheets, symmetry
foldmap plot model.cif -o figure.svg          # .svg, .png (300 dpi) or .pdf; -o is repeatable
foldmap plot model.cif -o fig.png --theme richardson --set residue_numbers=true
foldmap interactive model.cif -o view.html    # linked topology / contact map / 3D view
foldmap app                                   # Streamlit app: PDB IDs, UniProt accessions (AlphaFold), files
```

SVG output keeps text as text, so labels stay editable in Illustrator or Inkscape.

## Themes and styles

`foldmap styles` lists every theme and style key. Themes change only how things are drawn, never where:

| Theme | Look |
| --- | --- |
| `publication` | the default: bold fills, ribbon helices with depth shading |
| `journal` | muted palette, thinner lines, smaller text |
| `print` | greyscale, pale fills |
| `minimal` | outlines only |
| `presentation`, `cartoon` | glossy helices, curved loops, larger text |
| `richardson` | coloured by element type |
| `rainbow`, `shaded`, `trace` | N→C colouring; `trace` also colours loops and adds direction arrows |
| `flexibility`, `hydropathy` | coloured by B-factor or Kyte–Doolittle hydropathy |
| `alphafold` | AlphaFold pLDDT confidence bands, with loops and tails coloured residue by residue |
| `goodsell` | pale, illustrative |

Change any single key with `--set KEY=VALUE` (repeatable), or collect them in a YAML file for `--style-file`.
Useful keys: `color_by`, `palette`, `fill`, `loops`, `loop_color` (`residue` colours loops residue by
residue), `residue_numbers`, `highlight`
(`asu`, `protomer` or a chain list), `mark` (highlight chosen elements), `disulfides`, `glycans`,
`helices_310`, `helix_scale`, `strand_scale`, `font_scale`.

## Shaping a figure

The layout starts from the 3D structure: elements are projected onto the page, then arranged so that sheets stay
in order, contacting elements stay close and nothing overlaps. You can steer it:

| Option | Example |
| --- | --- |
| Swap two elements | `--swap α1,α3` |
| Nudge an element | `--move α2=1,-2` |
| Rename an element | `--rename res:A:167-182=Gd` |
| Highlight elements | `--set mark=α4=#2ca02c` |
| Named domain panels | `--domain N-lobe=res:A:13-59` or `--domains auto` |
| Orientation | `--up X,Y,Z`, `--view X,Y,Z`, `--rotate DEG` |
| Symmetry | `--symmetry auto\|off\|C3\|D2\|helical`, `--protomers "A,B;C,D"` |
| Assembly | `--assembly auto\|asu\|ID` |

Elements are referenced by label (`A`, `α2`, `η1`), `chain:label`, position (`#3`) or residue range
(`res:A:10-25`).

### Layout files

`--save-layout fig.yaml` records everything that defines a figure: theme, style changes, layout options, domains
and edits. Element references are stored as residue ranges, so they keep pointing at the same part of the chain
when a new model gains or loses an element. Re-render with:

```bash
foldmap plot model.cif --layout-file fig.yaml -o fig.svg
```

Anything given on the command line overrides the file. The Streamlit app saves and loads the same files.

## Interactive view

`foldmap interactive` writes one self-contained HTML page with the topology, a residue contact map and the 3D model
side by side. Hover any element, loop, map cell or atom and the same residues light up in all three. The 3D viewer
([3Dmol.js](https://3dmol.csb.pitt.edu/)) loads from a CDN, so the page needs an internet connection.

## Development

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev,app]"
.venv/bin/python -m pytest          # full suite, about a minute
.venv/bin/ruff check . && .venv/bin/ruff format --check .
.venv/bin/python scripts/gallery.py # re-render the images above
```

Test structures live in `tests/data`.

## Citing

If you use foldmap in published work, please cite it. GitHub's "Cite this repository" button gives the
reference in APA and BibTeX form (from [`CITATION.cff`](CITATION.cff)).

## License

foldmap is free software under the [GNU Affero General Public License v3.0 or later](LICENSE), with one
additional attribution term (see [`NOTICE`](NOTICE)). In short: you may use, study, change and share it; any
modified version you distribute, or run as a service for others, must be released under the same license with
its source code, and must keep the foldmap attribution. Copyright © 2026 Dawid Zyla.
