# topoplot Step 2: Layout, Loop Routing, First Rendering

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (chosen by the user for this project). Checkbox steps; every task is test-first (write test, watch it fail, implement, watch it pass, run the whole suite).

**Goal:** `topoplot plot FILE -o fig.svg -o fig.png` draws strands (arrows), helices (hatched bars), sequence loops, N/C termini and chain colours as a 2D topology figure, in two layout modes the user can compare: `projected` (3D-faithful) and `stack` (rows in sequence order, like the eLife Fig 2B).

**Architecture:** `frame` (global page axes from 3D) -> `layout` (place sheet blocks and helices, resolve overlaps) -> `route` (orthogonal A* loop routing around obstacles) -> `render` (matplotlib, named SVG ids, editable text). All consume the step-1 dataclasses (`SSE`, `Sheet`, `Backbone`).

**Tech Stack:** numpy, matplotlib (Agg for PNG; svg.fonttype none; pdf.fonttype 42), pillow (matplotlib runtime dependency, installed into the venv).

**Spec:** `docs/superpowers/specs/2026-10-08-protein-topology-plotter-design.md` (sections 4, 5 for chains/loops/style, 6 step 2). Reference style: `reference/uromodulin_fig2_manual.png`.

## Global Constraints

- Page units: 1 unit = one strand pitch (1.1 units between neighbouring strand centres, arrow shaft 0.62, head 1.0 wide, 0.7 long). 3D to page scale 0.25 units/Å in `projected` mode.
- A strand arrow points to its C-terminal end; a helix is drawn as a hatched bar whose N->C direction is its projected axis.
- Sheet orientation comes from 3D: strands drawn "up" when the sheet's mean strand direction projects onto page-up >= 0; left-to-right order follows the lateral 3D direction projected on page-x.
- Colours: one per chain from the colour-blind-safe Okabe-Ito set; loops black; text stays text in SVG; fixed sans-serif font.
- Deterministic output (no random seeds); same input gives byte-identical SVG.
- Loops only join consecutive SSEs of the same chain; dashed when a chain break (missing residues) lies between them. No inter-chain, disulfide, domain or glycan drawing yet (step 3).
- Strand letters A, B, C... follow sequence order **per chain** (not per sheet; the eLife figure letters across both sheets of a domain). Helices are labelled `αN` per chain. Rulings noted in the ledger.
- No git repo: tasks end with the whole suite green and a ledger line.

## Review Focus

1. A one-strand sheet, a sheet with only helices elsewhere, and a protein with zero SSEs (tripeptide) must still render (empty figure with a note, not a crash).
2. Boxed-in ports (a loop whose exit is blocked) must fall back to a visible curve and be reported, never hang or silently disappear.
3. Mirror-image safety: the page frame must be a proper rotation (right-handed), otherwise the figure shows the enantiomer.
4. Overlap resolution must terminate and leave no overlapping elements even when the 3D projection stacks many elements on one point (e.g. a helical bundle viewed down its axis).
5. Determinism: repeated runs are identical, so figures are reproducible.

---

### Task 1: View frame

**Files:** Create `src/topoplot/frame.py`, `tests/test_frame.py`.
**Interfaces:** Produces `frame.Frame(origin, u, v, w)` (3-vectors; `project(points) -> (n,2)`, `depth(points)`), `frame.view_frame(sses, rotate=0.0, flip_v=False) -> Frame`.
Behaviour: u,v = top two principal axes of SSE centroids, `w = u x v` (right-handed); sign of u so first->last SSE centroid has u>=0; sign of v so the first SSE sits at v >= mean (N-terminus up); `rotate` (degrees, counter-clockwise in the page) rotates u,v about w; `flip_v` mirrors on request only.
Tests: orthonormal and `det([u,v,w]) = +1` on 1UBQ/6ZS5; deterministic; rotate 90 maps old u to old v; `project` of a centroid equals `[(c-origin).u, (c-origin).v]`; fewer than 3 SSEs (1 or 2) still gives a valid frame.

### Task 2: Layout

**Files:** Create `src/topoplot/layout.py`, `tests/test_layout.py`.
**Interfaces:** Consumes `Frame`, `list[SSE]`, `list[Sheet]`, `Backbone`. Produces `layout.Placed(sse, cx, cy, length, angle, width, label)` with `n_port`, `c_port` (page xy), `exit_n`, `exit_c` (unit vectors pointing away from the element), `rect` (axis-aligned bbox); `layout.Layout(placed: dict[str, Placed], width, height)` ; `layout.build_layout(bb, sses, sheets, frame, mode="projected"|"stack") -> Layout`. Keys are `SSE.id`.
Behaviour: strands in a sheet at 1.1 pitch, drawn up/down from 3D (see constraints); helix angle = atan2(axis.v, axis.u), length `max(0.9, n*1.5*0.25*|proj|)`; `projected` places blocks/helices at projected centroids x0.25 then pushes overlapping rects apart (margin 0.8, bounded iterations, spring back toward start), `stack` packs elements in sequence order into rows of width <= 30 units, new row per chain; letters/labels assigned per chain.
Tests: every SSE placed exactly once; no two rects overlap (both modes; ubiquitin, 6ZS5, 6ZYA, and a synthetic pile of 12 helices at one point); antiparallel neighbours in a sheet have opposite arrow angles, parallel the same; left-to-right order matches projected lateral direction; ports: `|c_port - n_port| = length` and direction equals `angle`; deterministic; empty SSE list gives an empty layout.

### Task 3: Loop routing

**Files:** Create `src/topoplot/route.py`, `tests/test_route.py`.
**Interfaces:** Consumes `Layout`. Produces `route.Loop(a_id, b_id, points, dashed, fallback)` and `route.route_loops(layout, sses, bb) -> list[Loop]` (one per consecutive same-chain SSE pair; dashed when `bb.prev` has a gap between them).
Behaviour: grid 0.25 units over the layout bbox + 3 margin; strand/helix rects are obstacles; A* with 4 directions, turn penalty 2.5, penalty 6 for cells already used by earlier loops, shortest loops routed first; each path starts at the source port exiting along `exit_c` and ends entering the target port against `exit_n`; collinear points simplified; if no path (node cap 200k) a smooth curve (`fallback=True`).
Tests: two side-by-side strands joined by an orthogonal path that touches no obstacle rect and starts/ends at the ports; a walled-in port produces `fallback=True` and finishes; a second loop avoids the first's cells when an equal-length alternative exists; dashed flag set across a chain break (6ZS5 chain A has gaps); loop count on ubiquitin equals SSEs-1.

### Task 4: Rendering and `plot` command

**Files:** Create `src/topoplot/render.py`, `src/topoplot/palette.py`; modify `src/topoplot/cli.py`; create `tests/test_render.py`; extend `tests/test_cli.py`.
**Interfaces:** `palette.chain_colors(chains: list[str]) -> dict[str,str]`; `render.draw(layout, loops, sses, bb, title=None) -> matplotlib Figure`; `render.save(fig, path)` by extension (svg/pdf/png); CLI `topoplot plot FILE -o OUT [-o OUT...] [--mode projected|stack] [--rotate DEG] [--flip-v]`.
Behaviour: strand = arrow polygon with letter in white; helix = hatched rounded bar with `αN` label beside; loops = black round-capped lines with a white casing so later loops cut gaps in earlier ones (dashed for breaks, fallback as curve); N and C labels with short stub at chain termini; chain label and colour legend; empty layout draws a note. SVG: `svg.fonttype=none`, artists carry `gid`s `strand:<id>`, `helix:<id>`, `loop:<a>><b>`.
Tests: SVG parses as XML; contains one `strand:` gid per strand, one `helix:` per helix, one `loop:` per loop; contains `<text` elements (text not outlined); PNG and PDF files have correct magic bytes; two runs give byte-identical SVG; tripeptide renders without error; CLI writes both files, exits 1 on a missing file.

### Task 5: Acceptance renders

Render 1UBQ, 6ZS5, 6ZYA in both modes to `output/` as PNG + SVG, view them, compare with `reference/uromodulin_fig2_manual.png`, fix visible defects found (each fix gets a test where it is testable), then hand the PNGs to the user.

## Self-review
Spec coverage (step 2): layout (4.1-4.6 except YAML nudge, which is step 3+), routing, chain colours, labels, termini, vector SVG/PDF. Not in this step: interfaces, shared-sheet boxes, disulfides, domains, glycans, membrane. Known limit: ambiguous (branched) sheets keep the step-1 heuristic order.
