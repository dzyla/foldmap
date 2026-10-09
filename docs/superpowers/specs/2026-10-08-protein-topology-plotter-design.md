# topoplot: publication-grade protein topology diagrams

Date: 2026-10-08. Status: draft for review.

## 1. Goal

A local Python tool (CLI + API) that turns a PDB/mmCIF file into a vector 2D topology figure good enough to put in a paper without redrawing. Target quality: the hand-drawn Fig 2B of Stanisich, Zyla et al., eLife 2020 (uromodulin filament), or better. Existing tools (TOPS, Pro-origami, PDBsum, 2DProts, SSDraw) are single-chain, dated, or not topology plots.

Must show:
- Global orientation of helices and strands (true up/down direction).
- Strands in real sheet order, including sheet complementation across chains.
- Domains, chains, intra-chain connectivity (loops).
- Inter-chain interfaces.
- Disulfide bonds.
- N-glycans (SNFG symbols) at their Asn residues.
- Membrane proteins: orientation is meaningful (vertical axis = membrane normal, membrane drawn as a band).

Scope: general proteins (soluble, membrane, glycoproteins, multi-chain assemblies).

### Reference conventions (from the eLife Fig 2B)
Strands = arrows with letter labels; helices = hatched cylinders; one color per chain; loops = black connector lines; N/C labels; gray dashed boxes around shared sheets; brackets for interfaces; domain labels at the left. The one thing to improve is loop routing (hand-drawn loops cross and tangle).

### Test data
`data/raw_dl/6ZS5.cif` (3.5 Å; chains A, D protein; B, C glycans; 11 disulfides) and `6ZYA.cif` (extended filament; 14 disulfides). The reference figure is the user's Fig 2B (image copied to `reference/` by the user).

## 2. Non-goals (v1)
Interactive GUI/web app; predicting membrane placement (user supplies it); 3D rendering; conservation/B-factor coloring (planned later); NMR ensembles (first model only).

## 3. Architecture

Pipeline of independent stages with plain dataclass interfaces:

```
parse -> assign -> features -> layout -> route -> render
```

| Module | Responsibility | Depends on |
|---|---|---|
| `io` | Read PDB/mmCIF with gemmi; residues, coordinates, struct_conn, deposited SS records, OPM dummy atoms | gemmi |
| `ss` | Secondary structure: always the built-in numpy DSSP (Kabsch-Sander), because sheet pairing needs backbone H-bonds and deposited sheet records can reuse strands across sheets (6ZS5); the deposited annotation is only a test cross-check; merge to SSE objects (kind, chain, residue range, axis vector, centroid, direction) | numpy |
| `sheets` | Strand pairing graph from backbone H-bonds (parallel/antiparallel); order strands within each sheet; handle branched sheets and barrels by falling back to a cut-and-unroll ordering | networkx |
| `features` | Disulfides (SG-SG < 2.5 Å plus struct_conn), interfaces (SSE-SSE contacts < 4.5 Å across chains), glycans (attach HETATM/branch chains to Asn via struct_conn, map CCD codes to SNFG), domains (user ranges or file) | scipy |
| `layout` | Place sheet blocks and helices (section 4) | scipy |
| `route` | Loop and connector routing (section 5) | numpy |
| `render` | Draw to matplotlib; SVG with editable text and named layers (`gid`), PDF with embedded Type 42 fonts | matplotlib |
| `cli`, `config` | `topoplot in.cif -o fig.svg [--pdf] [--config layout.yaml]`; YAML for domains, colors, view axis, membrane normal, manual nudges | |

Each stage is testable on its own. `layout` and `route` never touch gemmi objects; they take SSE/sheet/feature dataclasses.

## 4. Layout (approach C: hybrid)

1. **Global frame.** Default: PCA of SSE centroids gives page x/y, with the third axis as view direction. `--view` accepts a vector. **Membrane mode**: `--membrane-normal x,y,z` or OPM dummy atoms; page y = normal, membrane band drawn, element tilt preserved.
2. **Blocks.** Each sheet becomes a rigid block: strands side by side in pairing order, arrows pointing along the real direction. The block's long axis is aligned to its projected mean strand direction and snapped to vertical in default mode. Helices are single elements drawn as hatched cylinders (tilt kept in membrane mode).
3. **Initial positions** from projected 3D centroids of blocks/helices.
4. **Constraint relaxation** (iterative, scipy): remove rectangle overlaps; preserve the rank order of projected positions along x and y where possible; pull sequence neighbors together; keep interface partners near each other; minimize total loop length. Deterministic seed.
5. **Layout modes:** `projected` (above) and `stack` (domains/chains as horizontal rows ordered along the sequence, as in the eLife figure). Default: `stack` for soluble multi-domain proteins, `projected` in membrane mode. User-overridable.
6. **Manual overrides.** YAML `nudge:` entries (element id -> dx, dy, flip) applied after the solver, so a figure is reproducible.

## 5. Routing and annotation

- **Loops:** orthogonal routing with rounded corners on a coarse grid (A*, penalties for crossings and bends), attached to strand/helix end ports. Crossings are drawn with a small gap in the lower line. Fall back to a Bézier when no route fits.
- **Disulfides:** short zig-zag or dashed arcs in a fixed color between the two Cys-bearing elements; label with residue numbers when requested.
- **Interfaces:** brackets on the side of the interacting elements, labeled I_A, I_B (names from config or auto); optional faint edges between contacting elements.
- **Shared sheets:** dashed gray box around a sheet with strands from more than one chain, labeled CS_x.
- **Glycans:** SNFG symbols (square/circle/triangle/diamond, standard colors) drawn as a small tree at the Asn site; legend in the figure margin.
- **Domains:** shaded rounded rectangles behind elements, label at the left.
- **Chains:** one color each from a colorblind-safe palette; overridable.
- **Labels:** strand letters/numbers auto-assigned per sheet (A, B, C... in sequence order), overridable in config; N/C termini.
- **Style:** vector only, sans-serif font, line widths and sizes in points, canvas width selectable (single/double column), so text is legible after journal scaling.

## 6. Build order (each step ships something usable)

1. `io`, `ss`, `sheets`; CLI prints SSEs and sheet ordering for a structure. Tested on 6ZS5/6ZYA.
2. `layout`, `route`, basic `render` (strands, helices, loops, chain colors). Compare to Fig 2B.
3. Feature layers: disulfides, interfaces, shared-sheet boxes, domains, glycans/SNFG.
4. Membrane mode and a small set of helical membrane protein tests.

## 7. Testing

- **Unit:** DSSP output vs the deposited SS annotation (agreement threshold set on first run, not guessed now); sheet ordering on a known antiparallel sheet; disulfide count equals 11 for 6ZS5; glycan-to-Asn assignment on 6ZS5; interface detection across chains.
- **Layout invariants (property tests):** no overlapping blocks; every SSE placed exactly once; routed loops end at the right ports; crossing count not worse than a naive straight-line baseline.
- **Rendering:** SVG parses; expected layer ids present; text stays as text; PDF embeds fonts.
- **Visual:** render 6ZS5/6ZYA to PNG and compare against Fig 2B by eye with the user. This is the acceptance test for step 2.

## 8. Decisions made (defaults, changeable)
- Name/package: `topoplot`. Python 3.13; deps: gemmi, numpy, scipy, networkx, matplotlib, PyYAML. No `mkdssp` dependency (not installed here).
- Drawing through matplotlib, not a hand-written SVG writer.
- Uses first model only; altlocs: first altloc.
- Project is not a git repository; I will not run `git init` unless asked.

## 9. Risks
- Loop routing quality is the hardest part; the Bézier fallback and manual nudges bound the worst case.
- Barrels and highly branched sheets have no unique 2D order; cut-and-unroll is a heuristic and needs a user override.
- Built-in DSSP may differ slightly from `mkdssp`; deposited annotation is used as a test cross-check (6ZS5: every strand residue called is in a deposited strand).
