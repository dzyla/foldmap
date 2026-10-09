# foldmap

Publication-grade protein topology diagrams from PDB/mmCIF files.

foldmap reads a structure, assigns secondary structure, finds β-sheets, helix bundles and symmetry, and lays
the elements out as a clean 2D topology diagram: coiled-ribbon helices, strand arrows and 3₁₀ helices on one
consistent scale, with orthogonally routed loops. Protein–DNA/RNA complexes, symmetric assemblies (cryo-EM and
crystallographic), disulfides and glycans are supported.

## Install

```bash
pip install .            # command-line tool
pip install ".[app]"     # plus the Streamlit app
```

Requires Python ≥ 3.11.

## Use

```bash
foldmap summary model.cif                     # elements, sheets, symmetry
foldmap plot model.cif -o figure.svg          # SVG, PNG or PDF
foldmap plot model.cif -o fig.png --theme richardson --set residue_numbers=true
foldmap styles                                # list themes and style keys
foldmap interactive model.cif -o view.html    # linked topology / contact map / 3D view
foldmap app                                   # Streamlit app in the browser
```

### Shaping a figure

| Option | Example |
| --- | --- |
| Swap two elements | `--swap α1,α3` |
| Nudge an element | `--move α2=1,-2` |
| Rename an element | `--rename res:A:167-182=Gd` |
| Highlight elements | `--set mark=α4=#2ca02c` |
| Named domain panels | `--domain N-lobe=res:A:13-59` or `--domains auto` |
| Orientation | `--up X,Y,Z`, `--view X,Y,Z`, `--rotate DEG` |
| Symmetry | `--symmetry auto\|off\|C3\|D2\|helical`, `--assembly auto\|asu\|ID` |

Elements are referenced by label (`A`, `α2`), `chain:label`, index (`#3`) or residue range (`res:A:10-25`).

### Layout files

`--save-layout fig.yaml` records the theme, style changes, layout options, domains and edits; element
references are stored as residue ranges so they survive model changes. Re-render with
`foldmap plot model.cif --layout-file fig.yaml`; anything on the command line overrides the file.

## Development

```bash
python -m venv .venv && .venv/bin/pip install -e ".[test,app]"
.venv/bin/python -m pytest
```
