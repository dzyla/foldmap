# Random-100, run 1: triage and visual review (plan 10 Task 11)

Run 1 drew the 100 sampled entries at the frozen commit `44a504f`, in all views, with RNApolis. **24 of 100
pass.** The 76 failing entries fail 14 checks: A6 56, A20 30, A15 16, A16 16, A4 11, A3 9, A17 8, A18 7, A7 4,
A9 4, A2 2, A19 2, A2b 1, A5 1.

This file sorts every failure into causes (C1…C34). A cause is one code-level reason, not one check. The causes
are ordered by the number of entries they touch, and each has a class:
- **(a)** a real figure defect;
- **(b)** a measurement or check defect: the check is wrong or too strict, with the evidence given;
- **(c)** an inherent limit that needs a ruling.

The visual review (34 contact sheets: the 24 passing entries and the 10 worst by check count) adds causes the
checks miss. They are marked "visual" and placed in the same order.

Notes on method:
- **Reproduction.** Every mechanism was traced on the downloaded file. Scratch probes ran `foldmap` from the
  frozen tree (`git archive 44a504f`). They wrap `cli.draw_spec` as `random100.py view` does, then recompute one
  measure with more detail: tether home and target chains, per-element copy residuals, page width and em. The
  probes are in `/home/dzyla/.claude/jobs/23e71d62/tmp/t11/probe/` and are not committed.
- **Line numbers** are at `2dcc0b2`/`f12bfdd`. They are the same as at `44a504f` for every file cited except
  `render.py` (bars and tethers); each cite also names the function.
- **"Page width"** is in layout page units. The text em in page units is `layout.text_em(width)`.

## What the commits after run 1 change

The commits `8ac832b..2dcc0b2` fixed:
- bars and tethers over bodies (closed strand outlines);
- the scope wording;
- the 5′ chord;
- A17's em for block-style stem names, and A19's comment.

`f12bfdd` (C-MEM-BARREL) added outer-membrane β-barrels.

Twenty failing entries were re-measured at HEAD `f12bfdd`, all views, with RNApolis (output in
`/home/dzyla/.claude/jobs/23e71d62/tmp/t11/head1/`):
- **Entries:** 6B5W, 7PIU, 7MW8, 8BW5, 9YLM, 9C75, 3OWW, 29RY, 7SFH, 8K8E, 9YKJ, 13MV, 6IDJ, 7SMR, 6JO6, 9CJI,
  7ARD, 2DGS, 8PQV, 1WWL.
- **Result:** every one fails the same checks with the same counts. Only the order of hit strings changed, plus
  one A6 hit fewer in 6JO6 domains (611 against 612; a ligand label × tether hit that a body now hides).
- **Coverage:** these entries cover:
  - A9: 6B5W, 7PIU, 7MW8, 7ARD. The scope-wording commit does not touch a chain that draws nothing, so C15
    stands.
  - A17: 8BW5, 9YLM, 9C75, 3OWW, 9CJI. The block-em fix changes no gate here.
  - the A6 kinds tail-length, loop-length, NA, own-loop and glycan;
  - A20.

**No run-1 failure is fixed by those commits** among the 20 re-measured. For the 56 not re-measured:
- **Possibly fixed:** only ligand-label × disulfide hits where a body now hides the bar (4UM3, 7AOD). 7SMR's
  domains view, re-measured, keeps all 12 of its hits. Neither is either entry's last failure.
- **Not checked at HEAD:** 9FWZ is not in the sample. `f12bfdd`'s β-barrel membrane changes A18 only where a
  closed barrel of ≥ 8 strands is in the membrane view (none of the seven A18 failures is one).

## Cause table

Ordered by entries. "Entries" counts entries with at least one failure of that cause. "Last" counts entries for
which it is the only cause left, so fixing it alone makes them pass. The "visual" rows count the entries a reviewer
named.

| Cause | Class | Entries | Last | Checks | Effort | Bench impact (likely) |
|---|---|---|---|---|---|---|
| C1 tail-length text placed blind, not a keep-out | a | 34 | 6 | A6 | S | none: text only (1UBQ-type long tails: AF-P04637, 2JTY may move a route) |
| C2 A20 gates same-chain tethers, against the ruling | b | 30 (28 clear) | 6 | A20 | S | none (runner only) |
| C3 loop-length and gap-count text placed blind, not a keep-out or page bound | a | 27 (+ visual 1) | 1 | A6, A7 | S | none: text only (a route may move where it joins the keep-outs) |
| C4 big pages: text held at 4 pt outgrows the fixed-size keep-outs | a | 26 (+ visual 7) | 1 | A6, A7, A16 names | L | 6Y7S, 8UTF, 8UUP, 25IN, 4OO8-size figures (every page wider than 44 units) |
| C5 ligand mark spot ignores loops, other ligands' tethers and labels, bars | a | 19 (+ visual 6) | 0 | A6 | S–M | ligand entries (1MBN-type, 2HHB, 4CH9): mark positions |
| C6 a loop's end segments may cross its own element's label | a | 23 | 0 | A6 | S–M | most entries with helices: loop or label positions shift |
| C7 visual: parallel loops 10–20 px apart cannot be told apart | a | visual 15 | – | (none) | M–L | most entries |
| C8 crowded assemblies: route fallbacks, loops over labels, element overlap | a/c | 12 | 0 | A3, A4, A5 | L | large entries only |
| C9 visual: straight tethers cross bodies and run along loops | a | visual 12 | – | (none) | M | ligand entries |
| C10 copies alike counts intermolecular-sheet strands | b (ruling) | 11 | 4 | A15 | S | none (check only) |
| C11 protein loops route before the NA strands and cross their letters | a | 8 | 3 | A6, A17 | M | NA–protein entries (1LMB, 4OO8 snapshots); NA campaign rerun |
| C12 tied copies drift apart in large assemblies and filaments | a | 6 (+ visual 2) | 0 | A15 | L | 6Y7S, 8UTF, 8UUP, 25IN |
| C13 visual: arrangement off the 3D contacts or the familiar published form (GFP central helix, globin α chain, Fab VL/VH, 5KQJ, 3MHT, 8P3N, 8OH7) | a | visual 7 | – | (none) | L | 1EMA, 2HHB |
| C14 NA runs cross NA letters | a | 6 | 1 | A6, A17 | M | NA entries; NA campaign rerun |
| C15 a protein chain with no helix or strand draws nothing, silently | a | 4 | 2 | A9 | S (note) / M (draw) | none |
| C16 NA layout metrics on small protein–NA entries | a/b/c | 4 | 2 | A17 | M | NA entries |
| C17 large membrane complexes: TM helices off the band (or no band) | a | 3 (+ visual 1) | 0 | A18 | L | 1C3W-type unchanged; large membrane entries |
| C18 pairs and copies matched by author numbers; untied copies go unchecked | b + a | 2 (+ visual 2) | 1 | A19 | S (A19) / M (ties) | none (A19); 2MJW-type ties |
| C19 visual: two-state pairs drawn and named independently | c (ruling) | visual 4 | – | (none) | M | none |
| C20 visual: additives listed in the facts, hidden from the figure with no note | c (ruling) | visual 3 | – | (none) | S | none |
| C21 domain panels of one chain overlap | a | 3 | 0 | A16 | M | 1LDM / 8UTF domain views |
| C22 copies with fewer than 3 elements judged inf | b | 3 | 0 | A15 | S | none |
| C23 automatic domains span chains | a | 3 | 0 | A16 | S–M | none |
| C24 NA route fallbacks | a | 3 | 0 | A3, A17 | L | NA entries; NA campaign rerun |
| C25 visual: neighbour ghosts not deduplicated; focus picks an uninformative chain | a | visual 3 | – | (none) | M | focus views |
| C26 loops run inside the membrane slab | a/c | 2 | 0 | A18 | M | membrane entries |
| C27 A18 requires 75 % band cover from tilted or offset TM helices | b/c (ruling) | 2 | 1 | A18 | S | none |
| C28 glycan stem and linked-residue label over labels | a | 2 | 0 | A6 | S | glycan entries (6ZS5 snapshot) |
| C29 visual: G-quadruplex tetrads not drawn; NA letters detached from the trace | a/c | visual 2 | – | (none) | M–L | NA entries; NA campaign rerun |
| C30 `--assembly auto` keeps a 40-chain ASU when assembly 1's pentamer holds two single-copy entities | a | 1 | 0 | A2, A2b, A4, A6, A16 | S | none |
| C31 domains view explodes the page when domains interlock across copies | a | 1 | 0 | A3, A4, A6, A7, A16 | M | none |
| C32 61-chain default view over 180 s | c/a | 1 | 0 | A2 | M | none |
| C33 visual: ligands that touch only nucleic acid are never drawn, silently | a | visual 1 | – | (none) | M | none (NA-only bench entries have no ligands) |
| C34 visual: a bulged strand drawn as its own column between paired strands | a | visual 1 | – | (none) | M | sheet entries |

Bench note: the bench is the 22 entries of `.superpowers/sdd/2026-10-10-08-deterministic-layout/baseline.json`.
- **Text-only causes:** C1, C3 and C20, unless the text also becomes a routing keep-out. Each keep-out added to
  `label_boxes(route=True)` moves routes, and so the signature, wherever a loop passed the text's spot.
- **Check-only causes:** C2, C10, C22, C18 (A19 part) and C27. They change no figure.

### Entries recovered per effort

Effort weights: S = 1, M = 2, L = 4. Each step takes the cause with the most new entries per effort, given the
causes already taken. "New passes" counts the automated checks only; the visual causes (C7, C9, C13, C19, C20, C25,
C29, C33, C34) are not in it.

| Step | Fix | New passes | Effort | Total passing |
|---|---|---|---|---|
| 1 | C1 tail-length text | 6 | S | 30 |
| 2 | C2 A20 per ruling | 7 | S | 37 |
| 3 | C3 loop-length text | 7 | S | 44 |
| 4 | C10 intermolecular sheets in A15 (ruling) | 6 | S | 50 |
| 5 | C18 A19 by sequence | 2 | S | 52 |
| 6 | C28 glycan and linked labels | 2 | S | 54 |
| 7 | C11 NA letters as obstacles for protein loops | 4 | M | 58 |
| 8 | C16 NA metrics (two to fix, two to rule) | 3 | M | 61 |
| 9 | C27 A18 tolerance (ruling) | 1 | S | 62 |
| 10 | C15 coil-only chains | 3 | M | 65 |
| 11 | C14 NA runs over letters | 2 | M | 67 |

After that, every remaining entry needs C4 (big pages), C8 (crowded routing), C12 (copy drift) or C17
(membrane complexes). Those are L-effort causes that touch the large assemblies, which make up most of S4, S6
and S8.

## Causes

### C1 (a): tail-length text placed blind and not a keep-out (34 entries, 6 last)

**Entries:**

| Group | Entries |
|---|---|
| Small pages | 29RY, 2DGS, 4XEI, 4XG1, 6FY9, 6FYF, 8R5X, 1SLV, 1WWL, 2B9V, 5X4Z, 6AD0, 7BRG, 7SMR, 7XSP, 8W9Y, 9GEP, AF_AFA0A4E9ES35F1, 7BBH, 8OZJ, 9P1X |
| Large pages | 5JB1, 5LDW, 6JLN, 6JO6, 6MHS, 7AOD, 7ARD, 7JU4, 7Y7C, 8ESZ, 8QEM, 8RVL, 8VWQ |

**Checks:** A6, as `tail-length:<end>:<chain> × loop:…`. Also `terminus × tail-length`, where a chain's own letter
and its own count touch: 8W9Y C:A–F and 8QEM N:M–X.

**Mechanism:**
- `render._draw_termini` (render.py:2976–3031) writes "NN aa" at `port + ex·(END_LABEL + 0.62)`. The `ex` is the
  stub direction, so the text lies on the stub's line beyond the N/C letter. It runs no obstacle test.
- `layout.label_boxes` (layout.py:3345–3373) keeps the stub and the N/C letter out of routes, but not the count.
  Loops therefore route straight through it. 29RY shows this: η1→A passes the "13 aa" under the N.
- Along a horizontal stub, the letter box (`_GLYPH_W/2 + 0.05` = 0.21) plus half of "123 aa" at 0.72 em
  (≈ 0.45) is wider than the 0.62 gap. The count then overlaps its own letter (8W9Y, 8QEM).

**Fix:** treat the count as part of the terminus. Add its ink box (at 0.72 em of `text_em(width)`) to
`label_boxes` under `{end}-len:{chain}`, and to `render` page bounds. Measure the gap from the letter's far edge,
not from the port. Pick the side (along the stub, else perpendicular) clear of loops, as `ligand_marks` picks a
spot.

**Fixture:** 29RY or 2DGS cut to the N-terminal 60 residues (< 100 KB). Test: no `tail-length × loop` in
`figcheck.text_collisions`. Add 8W9Y chain A for the letter overlap.

**Effort:** S. **Bench:** a keep-out moves routes where a loop crossed a long tail's count. AF-P04637 (p53, long
tails) and 2JTY are the likely members; the others have no tail ≥ 10.

**Batch C2a result (6JO6, 7BBH, 5LDW: tail-length × loop left after C1):** when all five near spots meet something
(crowded pages), the count now looks further out beside its letter and along the stub (`layout._TAIL_COUNT_FAR`,
up to 1.2 units), and when nothing is clear takes the spot meeting the fewest things instead of the first. 6JO6: 2
hits -> 0 in every view. 7BBH (focus none) and 5LDW keep one hit each per chain: the N letter is boxed in by a glycan,
a disulfide, strands and loops; the nearest clear spot is 2.7 units away (8 em), too far to read as the letter's
count, so it is not taken (page crowding, C4/C8). Bench identical; no snapshot changed.

### C2 (b): A20 gates same-chain tethers, against the ruling (30 entries, 6 last)

**Entries (all 30 A20 failures):** 13MV, 2YBE, 6IDJ, 6KTL, 7GLY, 5JDN (last cause), 2C4M, 2ILU, 3L9G, 4UM3,
4XG1, 5DWZ, 5LDW, 5SB5, 6DPV, 6JLN, 6JO6, 7ARD, 7JU4, 7OKA, 7QUP, 7SFH, 7SMR, 8ESZ, 8K8E, 8OZJ, 8PQV, 8RRX,
9GEP, 9P1X.

**Checks:** A20.

**Mechanism:**
- The ruling after Tasks 5–8 says A20 bounds only tethers to *other* chains; same-chain lengths are reported.
- `scripts/random100.py:criteria` (line 1696–1698) gates every `tether_lengths` entry over 6.5.
- `tether_lengths` (random100.py:1156–1166) does not record the target's chain, so the ruling cannot be applied
  as written.

**Evidence:** every A20 entry was recomputed with `render.ligand_marks` and `_home_chain` at `44a504f`. Of the
long tethers, 28 entries have only home-chain targets. Two have a real other-chain tether over 6.5:
- 6JLN `CLA:B609` → chain H, 7.08;
- 8ESZ `3PE:5604` → chain B2, 6.94.

Both stay failing on other causes.

**Fix:** `tether_lengths` returns `(gid, length, home chain, target chain)` from `_local_targets`' residues. A20
gates only `target ≠ home`; the rest go into the report as "reported". A test uses 13MV (home only, passes) and
a synthetic other-chain target.

**Effort:** S. **Bench:** none. The same-chain lengths stay a visual matter (C9).

### C3 (a): loop-length and gap-count text placed blind (27 entries + visual 1, 1 last)

**Entries:** 9ON5 (last), 2B9V, 2HKI, 5JB1, 5LDW, 5X4Z, 6AD0, 6JLN, 6JO6, 6MHS, 7AOD, 7ARD, 7BBH, 7BRG, 7JU4,
7QUP, 7SFH, 7XSP, 7Y7C, 8ESZ, 8PQV, 8RRX, 8VWQ, 9GEP, 9P1X, 9YKJ, AF_AFA0A4E9ES35F1. Visual: 3MHT (its grey "12 aa" gap
count sits between the D→η5 and F→η5 runs, so the loop it labels is unclear: the same blind placement).

**Checks:** A6, as `loop-length:… × loop/strand/helix/label/eta/ligand`. A7: 7BRG focus, where the count lies off
the page.

**Mechanism:**
- `render` (the loop drawing, render.py:2473–2489) writes "NN aa" beside the loop's longest straight run. It
  places it 0.15 to the right of a vertical run or 0.22 above a horizontal one.
- It tests nothing: not elements (9ON5: the count of `A:107-126>A:154-161` sits on strand 154-161), not other
  loops, and not labels.
- It is not in `label_boxes`, so the page bounds do not include it (7BRG, A7).

**Fix:**
- Choose among a few candidates: the longest runs, either side, near the run's middle. Take the first whose ink
  box clears element rects, other loops (as polylines) and `label_boxes`.
- Add the chosen box to the page bounds.
- Placement is after routing, so the box need not be a routing keep-out (bench-neutral). Optionally add it as a
  soft keep-out to later strands only.

**Fixture:** 9ON5 chain A (P2X2, about 330 residues; cut to the 100–170 window plus a neighbour). Test: no
`loop-length × *` collisions, and 7BRG's count inside the page.

**Effort:** S. **Bench:** none if render-only. The signature is layout and loops, and this text is in neither.

### C4 (a): big pages, where text held at 4 pt outgrows fixed-size keep-outs (26 entries, 1 last)

**Entries:**
- 2C4M (domains), 2QNA (focus, the last cause), 3L9G (domains), 4UM3, 5JB1, 5LDW, 6DPV, 6JLN, 6JO6, 6YF9 (focus),
  7AOD, 7ARD, 7BBH (focus), 7JU4, 7QUP (focus), 7SMR, 7Y7C, 8ESZ, 8OZJ, 8QEM (focus), 8RRX, 8RVL, 8VWQ, 8W9Y,
  9GEP, 9P1X.
- Visual (legend rows overlap; strand letters fuse into words): 4UM3, 7ARD, 6JLN, 7QUP, 8QEM, 7Y7C, 8VWQ.

**Checks:**
- A6: label × loop, label × label, terminus × loop, legend × legend, glycan on big pages;
- A16 name collisions: 4UM3, 5LDW, 7AOD, 7ARD, 7JU4, 7Y7C, 8ESZ, 8OZJ, 8RVL, 8VWQ;
- A7: 6YF9 and 7QUP focus.

**Mechanism:**
- `layout.font_pt` and `text_em` (layout.py:3376–3387) scale a page wider than `MAX_WIDTH_IN / IN_PER_UNIT` (24
  units) down to 7.2 in. The text is held at **at least 4 pt**.
- Past about 44 page units, the em in page units grows with the width. Measured with the probe:

  | Entry | Width (units) | Font (pt) | em (units) |
  |---|---|---|---|
  | 29RY | 25 | 7.0 | 0.34 |
  | 5X4Z | 68 | 4 | 0.52 |
  | 6JO6 | 161 | 4 | 1.24 |

- Every keep-out and spacing is fixed in page units: `_GLYPH_W, _GLYPH_H = 0.32, 0.42` (layout.py:44),
  `_element_label_box` (layout.py:3340), and the terminus boxes in `label_boxes`. So is the legend's
  `_LEGEND_ROW`: `render._draw_legend` (render.py:3041) puts rows one fixed step apart, so rows of 4 pt text
  taller than the step overlap.
- On a 160-unit page a 2-letter label is about 3.6× wider than its keep-out, so loops route through labels.
  Neighbouring labels and strand letters overlap too: 4UM3 shows "CFGEBA" fused across narrow strands.
- The A6 count rises from 0–10 below 44 units to hundreds above about 80 (report.json, area against count).

**Fix (two parts):**
1. **Keep-outs and the legend in ems** of the final page. The router runs after `build_layout`, so
   `layout.width` is known. `label_boxes(route=True)`, `ligand_marks`, glycans and the legend take
   `em = text_em(layout.width)` and size glyph boxes as `k·em`.
2. **Spacing in ems.** Element gaps (the `_separate` clearance) would also need the em to stop labels overlapping
   labels. That changes every large layout. The cheaper alternative is a page that grows in inches past
   `MAX_WIDTH_IN` (poster mode) so the font stays at 7 pt, which needs a ruling on maximum figure size.

**Fixture:** 2QNA focus (2 chains, importin-β; cut to chain A 600–760) for label × label at 48 units, plus a
synthetic wide page. Test: `text_collisions` empty on the fixture; keep-out width ≥ the label's ink width at
any page width.

**Effort:** L (part 1 alone is M). **Bench:** every bench entry wider than 44 units moves (6Y7S, 8UTF, 8UUP,
25IN, 6ZYA among them). They need re-recording with this reason.

**Batch C1 result (C4):** ruled (progress.md): text has a floor of 6 pt and the page grows instead.
- `layout.page_scale` sets a page's inches per page unit: `IN_PER_UNIT`, less on a page wider than `MAX_WIDTH_IN`
  (7.2 in), never so little that the text falls under `MIN_TEXT_PT` (6 pt), and never so much that a side passes
  `MAX_PAGE_IN` (48 in, about A0's long side). The text's em is then 0.341 page units on every page (it was up to
  1.24 on 6JO6), so the keep-outs sized in page units hold it again. Past 48 in the whole figure, text with it, is
  scaled to fit (no keep-out changes; only its print size falls under 6 pt), and the CLI, app and MCP notes say so.
- A PNG is 300 dpi up to 50 megapixels, and takes a lower dpi past that (said in the notes); SVG and PDF keep the
  page whole.
- The legend's row step is sized from its text (`layout.glyph_scale`, the one helper the keep-outs go through).
- "6 pt" is the figure's text size (helix names, terminus letters, legend); texts drawn at a set fraction of it
  (strand letters 0.9, counts 0.72) keep their fractions, as on a small page.
- Bench: both baselines `identical` (a protein layout does not depend on the em; the NA entries are under 29.5 page
  units, where nothing changes). Pages under about 29.5 units are drawn exactly as before.
- Snapshots re-recorded, each because its page grew (C4): `2HHB-journal`, `2HHB-shaded`, `2HHB-trace`,
  `2HHB-highlight` (40 units: 7.2 -> 10.2 in, nothing else moves once scaled back), `1TIM-symmetry`, `1TIM-shaded`
  (7.2 -> 14.2 in; text 4 -> 6 pt, so the legend entries are narrower), `1C3W-membrane` (7.2 -> 13.8 in; legend
  and the RET marks, placed from the text's size), `1KX5-nucleosome` and `30QM-ncrna` (7.2 -> 17.7 and 17.4 in;
  base letters, stem names and the strands routed round them are sized from the text, smaller in page units now).

### C5 (a): ligand marks ignore loops, later tethers, other ligand labels and bars (19 entries + visual 6)

**Entries:**
- From the checks: 2C4M, 3L9G, 4UM3, 5LDW, 5SB5, 6DPV, 6JLN, 6JO6, 7AOD, 7ARD, 7JU4, 7QUP, 7SMR, 7Y7C, 8ESZ, 8RRX,
  8VWQ, 9GEP, 9P1X.
- Visual: 1A0Z (the HEM boxes of B and D sit on the loop leaving α2), 1BXI (Ni on the α1→B loop), 9C0M (three
  Ca beads on loop lines), 3MHT (SAH on the η3 loop), 8P3N (ADZ on loops in both copies), 8OYQ (W6T over the
  α4→B link).

**Checks:** A6 `ligand-label × loop / ligand-tether / disulfide / ligand-label`. 6JO6 has 2,600: chlorophyll
labels crossed by other chlorophylls' tethers.

**Mechanism:** `render.ligand_marks` (render.py:1217–1275) runs after routing and is given `loops`. But its
obstacles are only `taken`:
- element rects;
- ghosts;
- `label_boxes`;
- the boxes of marks already placed.

It never tests loop polylines. Its "no tether crosses text" test covers labels and strand letters, but not the
labels of ligands placed earlier or the tethers that later ligands draw. Disulfide bars are not obstacles either.

**Fix:**
- Add loop segments (thin boxes, as the glycan placement already does at render.py:1712), disulfide bars and
  earlier marks' tethers to the spot test.
- When a later ligand's tether would cross an earlier label, try its next spot.

**Fixture:** 5SB5 chain C (one CA mark on loop C:10-28>C:48-51 in the domains view) or 1BXI (Ni). **Effort:** S–M.
**Bench:** entries with ligands: 1MBN-type, 2HHB, 4CH9, 6ZS5. Mark positions only, so the signature changes only
if the marks are in it (they are not: `bench.signature` is layout + loops). Snapshots of ligand figures change.

### C6 (a): a loop's end segments may cross its own element's label (23 entries)

**Entries:** 2C4M, 3L9G, 4UM3, 5JB1, 5LDW, 6DPV, 6JLN, 6JO6, 6MHS, 6YF9, 7AOD, 7ARD, 7JU4, 7QUP, 7SMR, 7Y7C, 8ESZ,
8OZJ, 8QEM, 8RRX, 8RVL, 8VWQ, 9P1X.

**Checks:** A6 `label:X:id × loop:…>X:id` (or `X:id>…`). Examples:
- 6MHS `label:B:275-283 × loop:B:265-271>B:275-283` in all four copies;
- 7SMR domains `label:A:237-239 × loop:A:213-234>A:237-239`.

**Mechanism:**
- `route._clear` (route.py:217–221): `mine = owners[k] in skip`. Rects owned by the loop's own elements are
  skipped for the first and last segment.
- `route_loops` (route.py:474) gives each `label:<id>` keep-out its element as owner. So a loop entering or
  leaving a helix may run straight through that helix's own label.
- `helix_label_pos` (layout.py:2859) puts the label on the bar's upper side whatever side the ports use.
- On small pages a short port stub rarely reaches the label. On big pages (C4) the label is larger than its box,
  so most hits are there.

**Fix:** exempt only the own *body* for end segments, not the own label. In `route_loops`, give `label:` rects
the owner `None`, or a separate `own-label` owner that `_clear` does not skip. Add `label_side` flipping in
layout when the label's side holds a port.

**Fixture:** 6MHS chain A window 228–306 (4 elements). **Effort:** S–M. **Bench:** any bench helix whose own
loop now has to go round its label. Measure first.
**Batch A result (interim):** `route._clear_helix_names` moves a helix name off the lines that strike it, after
routing (bench unchanged). Snapshots re-recorded: 1KX5-nucleosome and 30QM-ncrna (batch A, helix names only; 1KX5
re-recorded again in fix round 1 when NA traces and paths became obstacles). The router-side keep-out sized from the
drawn text moves to batch C with C4. No fixture under 150 KB reproduces the defect on a real page: it needs a page
wide enough that 2 px is a large part of a loop spacing (6MHS windows, 8UUP and 1GOX cuts all draw clean), so the tests
are synthetic.

**Batch C1 result (C6, router side):** every keep-out a loop routes round is sized from the text as drawn:
`layout.glyph` (the `_GLYPH_W` x `_GLYPH_H` glyph grown by `glyph_scale` of the drawn em) sizes the helix-name,
chain-end and tie-label boxes in `label_boxes`, the strands' 5′/3′ keep-outs (`route_na`, `na_ends`) and the tail
count's letter room; with C4's 6 pt floor the em is the same page units on every page. Routed loops that come within
figcheck's line clearance of a helix name's ink, before the mover runs: 8UUP 2 (one its own loop) and 45 at text
size 2, 8UTF 10 and 25IN 2 at text size 2 at `f57887a`; none now. The post-route mover (`route._clear_helix_names`)
stays as the fallback for lines no route holds; it now logs how many names it moved (the campaign reports
`names_moved`). Across 97 sample entries (default view; 7JU4, 7Y7C, 8VWQ left out for time) it moved 177 names at
`f57887a` (7AOD 39, 8RRX 27, 6DPV 25, 5JB1 19, 9GEP 16, 6MHS 15, ...) and 1 now: 6JLN's η4 of chain B, struck by
the fallback curve `E:17-39>E:42-47` (no route exists there, C8). The figure is unchanged at text size 1 (bench and
snapshots identical).

### C7 (a, visual): parallel loops 10–20 px apart cannot be told apart (15 entries)

**Entries:** 5Q42, 5HI7, 1F09, 9C0M, 7FDN (chain B), 1BXI, 2I5Y, 3L9G, 7JU4, 7ARD, 8OH7, 3MHT, 10UJ, AF_AFO42816F1,
8OYQ. Examples:
- 5Q42 three horizontals under the K/O/Q tips at y≈670/688/708;
- 9C0M two verticals for about 600 px down the left of the sheet;
- 1F09 five loops under the barrel with jogs and swaps.

**Checks:** none. No figcheck measure counts run-together length; `loop_crossings` counts crossings only.

**Mechanism:**
- The grid router (`route._Grid`/`_Router`) spaces routes by cell (`RES`) with `grid.used` as a soft cost.
- Two loops in neighbouring cells run one cell (about 0.2 page units, 15–20 px at 300 dpi) apart along their
  whole shared channel.
- Nothing penalises the length of a shared run, only crossings.

**Fix:**
- Raise the cost of a cell next to a used cell, scaled by run length, so later loops take a parallel channel
  two cells away or a different side.
- Add a figcheck metric: the length over which two loops stay within 0.25 units. Report it first, then gate.

**Effort:** M–L. **Bench:** most entries (routes move). A ruling is needed before it gates.

### C8 (a/c): crowded assemblies, with route fallbacks, loops over labels and element overlap (12 entries)

**Entries:** 4UM3, 6JLN (A5), 7ARD, 7JU4, 7QUP (focus), 7Y7C, 8OZJ (domains), 8QEM (focus), 8RVL (domains), 8VWQ,
9GEP, 9YKJ (beta).

**Checks:** A3 protein loops (`drawn as plain curves`), A4, A5 (6JLN `b:398-402 × b:429-433`).

**Mechanism:**
- `route_loops` (route.py:493–517) falls back to `_curve` when the grid router finds no path. It accepts a route
  over a label (`over`) when the keep-out cells are the only way (`_Grid.soft`, route.py:85–87).
- In 40–60-chain assemblies the channels between packed elements are fewer than the loops needing them.
- On big pages C4 makes this worse (text that has no keep-out room).
- 9GEP (12 chains) and 9YKJ (2 chains, beta view) are mid-size cases: one dense spot each.

**Fix:** after C4, measure again. Options:
- a wider `_separate` clearance in dense rows;
- a second routing pass that widens only the gaps a fallback needed;
- for the largest, a ruling that the default for > ~30 chains is a focus view (plan 07's ruling covers 6JLN-size
  figures only partly).

**Effort:** L. **Bench:** large entries.

### C9 (a, visual): straight tethers cross bodies and run along loops (12 entries)

**Entries:** 5Q42, 5HI7, 7LX9, 9C0M, 1BXI, 2I5Y, 3L9G, 7QUP, 6JLN, 9YLM, 3MHT (a SAH tether about 500 px to the α7
loop), 8P3N (an ADZ tether about 300 px). Examples:
- 5Q42 Ni tether about 400 px over strand B and A's arrowhead;
- 7LX9 YGM→α6 tether on top of the dashed α5→α6 gap loop for about 230 px;
- 9YLM Zn tethers across the A/B and C/D hairpins;
- 2I5Y NAG stalk across four gp120 strands.

**Checks:** none. A20 measures only length. A6 counts tethers only where they cross text.

**Mechanism:**
- `render` draws each tether as one straight segment from mark to target.
- `ligand_marks` picks the spot so tethers cross no *text* (render.py:1255–1257), but not so they avoid
  elements or run clear of loops.
- 8ac832b now hides the part under an opaque body, but the line still visibly enters and leaves the body.

**Fix:**
- Score candidate spots by the tether length crossing element bodies and by the length running within about
  0.15 of a loop. Pick the lowest.
- Drop a target whose tether would cross a body when another target of the same element is reachable.
- Add a reported metric (tether × body length) to figcheck.

**Effort:** M. **Bench:** none (marks are not in the signature); ligand snapshots.

**Batch B result:** the spot of a mark is chosen among the clear spots up to 4 rings past the first by the length of
its tethers over element bodies (`render._over_bodies`; the body a target lies on not counted), no tether lengthened
past the first spot's reach. Tether over bodies: 7LX9 1.15 -> 0, 3MHT 0.56 -> 0, 9YLM 3.15 -> 2.5, 1GOX-A 12.0 ->
11.3; 5Q42, 5HI7 and 8P3N unchanged (no spot in reach does better). Running along loops is not scored yet.

**Batch C2a result (C9, running along loops):** a spot's score adds the length its tethers run along a loop's line
(`render._along_loops`: beside a segment within 20 degrees of its direction, both its ends within 0.15 of the tether's
line) to the length over bodies. Along loops, as drawn: 5HI7 1.5 -> 0, 8P3N 0 -> 0 (its ADZ tether over loops is a
crossing, scored as bodies are not), 5Q42 0.72 -> 0.72 (its Ni reaches a residue on the horizontal α1 -> C loop
from the only clear spot level with it; no spot within the 4 rings past it does better). Snapshot re-recorded:
1KX5-nucleosome (the MN:E3132 tether ran 0.54 along a loop; its mark moved 0.8 units). Bench identical.

### C10 (b, ruling): copies alike counts strands of intermolecular sheets (11 entries, 4 last)

**Entries:** 7IMM, 9PXS, 6YPX, 6YQD (last), 1SLV, 1WWL, 2ILU, 3L9G, 5DWZ, 8QEM (focus, part), 9P1X (focus, part).

**Checks:** A15.

**Evidence:** per-element residuals of each failing tie at `44a504f`:

| Entry | Tie | Deviation | Elements | In a cross-chain sheet | Without those |
|---|---|---|---|---|---|
| 7IMM | A2~A1 | 5.63 | 10 | 2 | 0.04 |
| 9PXS | A2~A1 | 14.21 | 24 | 1 (strand C) | 0.00 |
| 2ILU | A2~A1 | 27.4 | 38 | 10 | 0.00 |
| 2ILU | A4~A3 | 29.1 | 38 | 10 | 0.00 |
| 3L9G | A2~A1 | 24.2 | 20 | 12 | 0.00 |
| 3L9G | A3~A4 | 22.2 | 20 | 12 | 0.00 |
| 1WWL | B~A | 12.3 | 18 | 12 | 0.00 |
| 5DWZ | E~B | 11.4 | 22 | 2 | 0.12 |
| 5DWZ | F~C | 13.5 | 30 | 3 | 0.01 |
| 6YPX | BBB~AAA | 4.57 | 8 | 5 | 0.00 |

- 6YQD and 1SLV A2~A1: every matched element is in a cross-chain sheet.
- 1SLV B2~B1 is 0.00.

**Mechanism:**
- A sheet whose strands come from two tied copies (a dimer's continuous sheet: 7IMM's A B B A, 9PXS's swapped
  strand C) is one layout item. Task 2 rules that it is never tied: it is drawn once, in the middle, faithful to
  3D.
- Its strands therefore cannot sit at their copy's tied place.
- `figcheck._copy_pairs` (figcheck.py:544–564) still matches and fits them. On the page every *other* element of
  each pair is drawn exactly alike.

**Ruling needed:**
- **Proposed:** A15 leaves out elements whose sheet (`lay.sheets`) holds strands of more than one tied copy, and
  still needs ≥ 3 matched elements (else, with C22's rule, report "n/a: shared sheet").
- **Visual counterpart:** a shared sheet drawn once is what published figures do (e.g. the continuous
  β-sandwich of dimeric enzymes).
- **Alternative:** draw such a 2-fold as a stated mirror (`mirror` in `_symmetry_ties`). The `dev_other_hand`
  column shows that would not give 0 either.

**Effort:** S (check). **Bench:** none.

### C11 (a): protein loops route before the NA single strands and cross their letters (8 entries, 3 last)

**Entries:** 6SO9, 9CJI, 9K9V (domains) are last. Also 3L2V (domains), 7XSP, 7Y7C, 8VWQ, 9GEP (domains, a stem
name).

**Checks:** A6:
- `na:nt:<chain>:<n> × loop:…`. 6SO9: RRM loops `A:522-525>A:534-544` and `A:581-583>A:586-592` cross U/C/A of
  the bound 6-mer.
- `na:end:… × loop` (9CJI 3′).

A17 `bead_overlaps` in 9CJI.

**Mechanism:**
- `route_loops` (route.py:543–546) routes every protein loop first.
- Then `route_na.route_nucleic` (route_na.py:64–100) routes the single strands and only then places their
  letters (`_place_letters`, route_na.py:187) and ends.
- So protein loops never see NA letters or 5′/3′ labels. `route_nucleic` weighs crossing a protein loop only
  softly against crossing strands (route_na.py:92–100).

**Fix:**
- Route NA single strands, letters and ends first. Their letter boxes become keep-outs (`grid.keep_off`) for
  protein loops.
- Or, after NA placement, re-route the protein loops that cross a letter box.

**Fixture:** 6SO9 (RRM + 6-nt RNA, 2 chains, small). **Effort:** M. **Bench:** 1LMB, 1ZAA, 4OO8 snapshots; NA
campaign rerun (constraint).

**Batch B result:** after `route_nucleic` places the NA text, a protein loop crossing it (letters, stem letters on
their row, 5'/3' labels: `layout.na_text_boxes`) is routed again with that text tolled on the grid, and takes the new
route only when it crosses less NA text and no more text in all. 6SO9 4 hits -> 0; 9CJI and 7XSP default views
clear. No snapshot changed.

### C12 (a): tied copies drift apart in large assemblies and filaments (6 entries + visual 7FDN)

**Entries:** 7OKA (B~A 1.12, C~A 1.12, one element B:225-231), 9GEP (E~A 1.17, F~B 0.81, mirror ties), 6JLN, 9P1X
(focus none), 8QEM (focus none), 7QUP (focus none). Visual: 7FDN (the copies' loops routed differently), 8QEM
(one chain-N copy lacks strand G).

**Checks:** A15. The deviations left after removing C10's cross-chain strands:

| Entry | Tie | Deviation left |
|---|---|---|
| 6JLN | e~E | 15.1 |
| 6JLN | d~D | 6.1 |
| 6JLN | b~B | 5.3 |
| 9P1X | C~V | 1.67 |
| 9P1X | D~X | 2.26 |
| 8QEM | per copy | 0.4–10.7 |
| 7QUP | all 60 ties | up to 44 |

**Mechanism:**
- `layout._copy_ties` (layout.py:2135–2149) records the ties the symmetry asked for, "whether or not _separate
  keeps them tied".
- `_separate` (layout.py:756) and the push-apart passes move elements individually in crowded rows.
- In a filament or ring drawn whole (`--focus none`), copies at the seam are placed by their own projection.
  7QUP shows two forms of one tubulin: one spread out, one compact.
- 7OKA's single element and 9GEP's small deviations are the same push acting on one element.

**Fix:**
- Make `_separate` move tied orbits as units (as `_part_domains` does for domains, layout.py:2535).
- For whole filaments, tie each copy to its neighbour's drawn shape, not to protomer 0's projection.

**Effort:** L. **Bench:** 6Y7S, 8UTF, 8UUP, 25IN (ties) re-recorded.

**Batch C2a result (9GEP):** not a push on one element: the nucleosome and its bound N/M chains are drawn as mirror
images about a dyad in the page, and `_separate`'s symmetric pushes cleared the overlaps only after about 3,000
passes (the limit is `_SEP_ORBIT_SWEEPS`, 1,000) and in 1.07 times the one-by-one area (limit `_SEP_SYM_ROOM`, 1.05),
so the copies were moved one by one and drawn unalike. With a mirror (the figure says the copies are mirror images)
the symmetric pushes now run as many passes as the one-by-one ones (`_SEP_SWEEPS`) and may take up to
`_SEP_MIRROR_ROOM` (1.10); other ties keep both limits. 9GEP copies 1.17/0.81 -> 0 (every tie below 1e-12), entry
passes. Bench identical (1LMB, the bench's mirror, unchanged); NA campaign identical; no snapshot changed. The test
is synthetic (tests/test_symmetry.py, like the other separation tests): a real cut reproducing it needs the
nucleosome, its DNA and the 469-residue chain N (over 180 KB).

### C13 (a, visual): arrangement off the 3D contacts or the familiar published form (7 entries)

- **1F09 (GFP):** the central helix (57–64) and its chromophore sit above the sheet, and the 64–69 loop detours
  round the right edge. Published GFP topologies draw the central helix through the barrel (COMPARE §2.4).
- **1A0Z (haemoglobin):** in the α chains, G and H (a packed pair in 3D) are about 450 px apart with an empty box
  between them, while the β chains keep G/H stacked. The two globin chains are drawn differently (COMPARE §2.8).
- **2I5Y:** the Fab's VL sits beside CH1 instead of on VH; L's V–C linker crosses the whole page.
- **5KQJ:** α1 is drawn right of the sheet beside A. In 3D it is left and below, under C/D, touching α3 (6.1 Å),
  which the page puts at the far left.
- **3MHT:** α8 (308–323) is drawn upper right above α4. In 3D it is the lowest element, level with strands A/B.
- **8P3N:** the two CheA copies are stacked front and back along the 2-fold (depth ±16 Å). Chain A is drawn as a
  translated copy of B, so its left/right order is reversed against 3D (x rank correlation −0.52), with no note.
- **8OH7:** the repeats drift down about 360 px from left to right, where in 3D they sit level.

**Mechanism:**
- The projected layout keeps each chain's own projection.
- When two elements overlap in projection, `_separate` moves one aside. It keeps neither the moved element's
  contacts nor its vertical order (5KQJ, 3MHT).
- There is no fold template, and nothing holds a helix inside a barrel.
- A copy behind its mate along the view is translated, not mirrored (8P3N: the same choice as C10's dyads in the
  page).

**Effort:** L. These are V3/V8 items, not gates. A reported metric (rank correlation of page against projected
3D position, per chain, as the batch-2 reviewer computed) would make them measurable.

### C14 (a): NA runs cross NA letters (6 entries, 1 last)

**Entries:** 8SA2 (last), 3OWW, 9C75, 8QEM, 7Y7C, 8VWQ.

**Checks:**
- A6 `na:nt:… × na:run:…`. Examples: 8SA2 A:26-28 under run 18, 3OWW B:18/23/33, 9C75 R:38-40.
- A17 `na_crossings`: 3OWW 2.

**Mechanism:**
- `route_na.route_nucleic` routes runs around stem bodies (`own_body`, route._clear).
- It places letters afterwards (`_place_letters`, route_na.py:187), so a later run is not kept off earlier
  letters, and letters are not pulled off runs drawn later.
- nacheck's `letter_overlaps` counts letter on letter only, so A17 passes 8SA2 while figcheck counts the line
  through the letter.

**Fix:**
- Place each run's letters as soon as it is routed, and add their boxes as soft keep-outs (`_Router.keep_off`)
  before the next run.
- Or re-route runs that cross letters.

**Fixture:** 3OWW chain B (glycine riboswitch, 2 chains). **Effort:** M. **Bench:** NA snapshots; NA campaign
rerun.

**Batch B result:** `na_labels._place_letters` counts a spot a drawn strand (run, end, gap) passes through, within
the line clearance of the letter's ink, as not clear, so the letter takes its next clear spot (with none clear, the
roomiest, as before: making a struck spot roomless crowded 1HQ1's B:148 onto another letter) (8SA2 14 hits -> 0;
3OWW-AB fixture A:87 -> 0; 30QM 1:59/1:60). Snapshots re-recorded: 1EHZ-trna (A:48's letter moves off the corner
of the strand to the P4 loop) and 30QM-ncrna (1:59 and 1:60 off run 2). Fix round 1: the letter and reroute clearances are figcheck's 2 px
(and half a line) at the page's text size (`layout.line_clearance`), not ems; hairpin-loop letters stay in their
stem's footprint. Snapshots re-recorded again: 1EHZ-trna (A:14, A:26) and 30QM-ncrna (1:41, 1:61, 1:104, 1:111), letters
a few pixels further from strands.

### C15 (a): a protein chain with no helix or strand draws nothing, silently (4 entries, 2 last)

**Entries:** 6B5W (the whole entry is one 19-residue lasso peptide: an empty page, area 0), 7MW8 chain K, 7PIU
chain P (a peptide), 7ARD chain q.

**Checks:** A9 (`figcheck.scope_honest`, figcheck.py:785: chain read, not drawn, not said).

**Mechanism:**
- The layout draws only elements (`layout.placed`). Termini and loops come from elements (`termini`,
  layout.py:2868).
- A chain with no SSE after `features` has no item, and neither the CLI note nor `scope.drawn_scope`
  (scope.py:30) lists it.

**Fix:**
1. **Note (S):** `drawn_scope` lists "chain P (19 residues, no helix or strand) not drawn" and the CLI says it.
2. **Draw (M):** a coil-only chain gets a short wavy "peptide" glyph with N/C and its chain name, placed by its
   contacts like a ligand mark. This is the published convention for bound peptides.

An empty figure for 6B5W should at least say why.

**Effort:** S / M. **Bench:** none.

### C16 (a/b/c): NA layout metrics on small protein–NA entries (4 entries, 2 last)

| Entry | Failure | Note |
|---|---|---|
| 9YLM | stems_nearer_frac 0.0 | (a or c) Three Cys2His2 zinc fingers lie in the major groove along one 12-bp duplex. The figure packs the protein below the duplex as one compact group. The touching helices (α1, α2, α3) average farther than the median element, because the small sheets sit between them and the DNA. Published ZF figures spread the fingers along the DNA. Fix: contact-driven placement for proteins whose elements touch a stem at many points along it (M); else a ruling like M1 (1ZBI). |
| 8BW5 | pair_f1 0.8, stem_recall 0.5 | (b/c) A duplex/quadruplex aptamer. The G-quadruplex tetrads are G·G Hoogsteen pairs, which foldmap's pair reading does not form as stems (RNApolis lists them). Ruling: count quadruplex tetrads out of pair_f1, as M4 does for MIN_LOOP. |
| 3L2V | angle_median 40.6°, angle_max 58.3°, contact_auc 0.48, stems_nearer 0.5 | (a) The PFV intasome holds two viral DNA ends and target DNA with a sharp bend. Stems are drawn at angles 40–58° away from 3D. Under investigation: the frame (`na_torn_by_unrolling`, layout.py:3394) may unroll the tetramer and turn the duplexes. |
| 3OWW | spearman 0.484, na_crossings 2 | (a) Glycine riboswitch: the stem order along the page does not follow 3D, plus C14's crossings. |

**Effort:** M overall. The NA campaign rerun applies.

**Batch B classification (per entry, probed at HEAD):**
- **8BW5 (c, inherent):** RNApolis' two extra canonical pairs, F5·F37 and F10·F33, are isolated A·T pairs with a
  stretched bond (3.44/4.33 Å and 3.87/3.33 Å, C1'-C1' 10.5 and 12.1 Å) beside the quadruplex. foldmap takes a
  stretched pair only stacked in register on a strict one (`nucleic.pair_nucleotides`, tuned on the 39 NA entries),
  so it draws 4 of 6 pairs (pair_f1 0.8, stem_recall 0.5). The G-tetrads are not the cause (RNApolis' reference counts
  canonical pairs only); C29's note now names them. Needs a ruling: accept, or count isolated stretched pairs.
- **9YLM (a, real):** the three zinc fingers' recognition helices touch the duplex, but the layout puts their small
  β-hairpins between them and the DNA (stems_nearer_frac 0). Needs contact-driven placement along the duplex; not
  done in batch B.
- **3L2V (a, real):** the C2 intasome is unrolled protomer by protomer, which turns each viral DNA 40-58° from its 3D
  axis (angle metrics measure the unrolling; contact_auc 0.48, stems_nearer 0.5). Not done in batch B.
- **3OWW (a, real):** the domain-swapped dimer (both chains, 3OWW ruling) has its stems in an order that follows 3D
  poorly (spearman 0.484) and two run-run crossings (na_crossings 2, not letters: C14 cleared its letter hit). Not done.

**Batch C2a result (3L2V):** `layout.na_torn_by_unrolling` gains rule (c): when most stem pieces (more than half)
are each held by two protomers (≥ 3 residues each, `_held_by_two`), unrolling would set each piece beside one protomer
with the other's contacts a protomer's width away, so the plain frame is used. 3L2V is drawn in the plain frame, which
sees its dyad in the page: the two halves are mirror images and the two viral DNA ends make a V above the inner
chains (angle_median 40.6 -> 0, angle_max 58.3 -> 0, contact_auc 0.481 -> 0.833, stems_nearer 0.5 -> 1.0; entry
passes in every view). 30QM (one of six stems shared) stays unrolled, as the NA campaign's C4 ruling keeps it; over
the 39 NA entries and the 18 random-100 NA entries only 3L2V changes frame. Fixture 3L2V-A (24 KB). Bench identical.

**Batch C2a result (9CJI):** bead_overlaps 2 was not a placement of the beads: the crRNA's 7-nt 3' end (B:22-28)
was turned onto the line of its 5' end (B:4-7), which leaves the same end of the stem, so nucleotides B4/B25 and
B5/B26 were drawn 0.1 apart with their beads on top of each other (a free end's straight or turned drawing was
checked against elements, hairpin loops, its label and its partner's doorway, never against the strands drawn
before it). `na_ends._tail` now takes, before any other, a drawing whose nucleotides keep a bead's width (`_BEAD`) off
every strand drawn so far (`_Strands.near`); only when none does are the old drawings tried. 9CJI passes in every
view. NA campaign: every metric of the 39 entries identical; bench identical. Fixture 9CJI-AB (21 KB).

**Batch C2a result (3OWW, na_crossings):** both crossings came from P5's hairpin loop, which lacks B22: it is drawn as
the 4-nt 3' end B:18-21 of one piece, a dashed gap and the 1-nt 5' end B:23 of the next. B:18-21 went straight on,
away from B:23, and the gap came back round the stem across the run B:75-81 and its own end. A 3' free end the
protein does not touch, running into a gap, now tries first the turns that bring its tip nearest the stem where the
molecule goes on (`na_ends._toward`, `target`): na_crossings 2 -> 0 (na_run_crossings, reported, 2 -> 3). Only 3'
ends: turning the 5' end after the gap too pulled 1URN's C15 onto C12 (tried, dropped). NA campaign: every metric
of the 39 entries identical; bench identical. Fixture 3OWW-P5 (9 KB).

**Batch C2a, 3OWW spearman 0.484 (for a ruling):** the drawn stem centres follow the projection poorly
(projected vs 3D 0.792, drawn vs projected 0.592): the stems are drawn about twice their projected size (letter rows,
rungs at the schematic step), so `_separate` pushes them apart by whichever way frees an overlap, and P4 (B:1-7 with
B:82-88) lands 11 units below its projected place under P5. Holding the projected order through the expansion is a
change to the shared placement (`_attract`/`_separate`) that all 39 NA entries go through; not done here.

### C17 (a): large membrane complexes, with TM helices off the band or no band (3 entries, + visual 1)

**Entries:**
- 5LDW (complex I: about 25 TM helices wrapped into a second row above the band, fractions −1.34 to 0.71);
- 6JLN (photosystem II dimer: the membrane view is the default layout with a band laid over it, the band about
  half the TM-helix length, single-TM subunits above it);
- 6JO6 (no band drawn);
- visual: 7ARD (complex I with lipids is not detected as a membrane protein, so it has no membrane view).

**Checks:** A18.

**Mechanism:**
- `cli._membrane_band` (cli.py:442–474) builds one band per stack row from `_face_medians`.
- The membrane view keeps the projected layout. Only the stack (`layout._stack`, layout.py:1360) aligns depth.
  When the in-band row overflows its width, `_stack` wraps it into rows, and TM helices land in a row with no
  band (5LDW).
- 6JO6: `membrane.find_membrane` returned no TM element crossing both faces. This is under change by
  C-MEM-BARREL; re-measure after it lands.
- 7ARD: no membrane is estimated at all (A18 "–").

**Fix:**
- For the membrane view, place every TM element in one band row and wrap soluble items outside it.
- Size the band from all TM elements, not medians per row.
- Re-measure 6JO6/7ARD after C-MEM-BARREL.

**Effort:** L. **Bench:** 1C3W, 1BL8, 1PV7 membrane fixtures must stay as they are.

**Batch C2a result (C17, projected membrane view):** after the layout `cli._off_band` checks the transmembrane
elements against the band they make, as A18 judges them; when one is off it (entering through no face, or not
reaching past the middle), the view is laid out again (`build_layout(level=, band=)`): the items holding
transmembrane elements keep their projected height (pushed and slid sideways only, `_push_apart`/`_slide_clear`
`level`), each moved the least way to meet a face when drawn wholly inside the slab (`_to_face`), and the band is
drawn at the membrane's faces as projected (`_faces_on_page`). Views whose elements already cross the band are not
touched (1C3W, 6PS3, 7PIU, 5JDN and the 1C3W-membrane snapshot unchanged). 6JLN: TM helices off the band 20 -> 0
(3 slab loops instead of 1, 6 fallbacks instead of 4: its other causes are C8/C12); 5LDW 6 -> 0 (one slab loop,
L:448-454>L:457-462); 1BL8 (fixture) 4 -> 0. Bench identical. 6JO6 still draws no band: its slab is 0.174 polar
(limit 0.15, `membrane._MAX_POLAR`), with 46 long elements crossing it; the soluble entries of the sample with 10 or
more crossing elements are at most 0.39 hydrophobic (limit 0.55) -- a rule "polar up to 0.20 when 10 or more long
elements cross" would find it, and needs a ruling (a test needs a 70 KB cut).

**Batch C2a, loops in the slab (for a ruling):** the three left are jogs beside an element end, not runs along the
slab: 7SMR B:171-173>B:190-201 (1.75 sideways 0.38 inside the face, over an element end at the face), 8W9Y
C/F:225-249>271-280 (k = 2 after a 0.31 port stub: 1.25 sideways round α8's end, from an interfacial helix's port
inside the slab), 6JO6 (with a band) 5:167-169>5:179-207 (0.72 sideways 0.42 inside, from a port 0.03 outside the
face). A18 excuses a sideways run only as one of a loop's first two segments from a port inside the band: a ruling
could excuse the first sideways run of more than 0.5 (a stub under 0.5 not counted) before the loop first leaves
the band.

### C18 (b + a): pairs and copies matched by author residue numbers; untied copies go unchecked (2 + visual 2)

**Entries:** 1RM5 and 2HKI (A19, one pair, both fail). Visual: 2MJW (S100P B and D not tied), 8QEM (one copy
lacks strand G).

**Evidence:**
- **1RM5 and 2HKI** (GAPDH S188A against A4): the elements are the same, but the numbering runs ahead in 2HKI by
  1→3 residues: 1-6/2-7, 25-31/28-34, 63-66/66-69, 115-118/118-121, 315-330/318-333. `model.pair_places` pairs by
  overlapping author numbers, so only 14 of 28–30 pair. `random100.pair_check` (random100.py:1702) fails it.
- **2MJW:** `homolog_chains` finds `D = B + 94` (chain D is numbered 95–188). DSSP on the NMR model gives D two
  2-residue strands and a split helix, so `_symmetry_ties` ties A~C but not B~D. A15 judges only requested ties,
  so the unlike drawing of B and D goes unmeasured.

**Fix:**
- A19: map residue numbers by a sequence alignment of the two chains (gemmi `align_sequence_to_polymer` or a
  plain Needleman–Wunsch on one-letter codes) before `pair_places` (S).
- Ties: let `_symmetry_ties` tie copies whose elements pair one-to-one for ≥ 80 % after dropping elements under
  4 residues (M).
- figcheck: report homolog copies of a symmetric assembly that are not tied.

**Bench:** none for A19; 2MJW-type ties.

### C19 (c, visual, ruling): two-state pairs drawn and named independently (4 entries)

**Entries:** 10UJ/9YKJ and 8P3N/8OYQ (S11 pairs). One small DSSP difference shifts every later generic name:
- residues 37–50 are α2 in 10UJ and α1 in 9YKJ;
- 128–136 are α6 in 8P3N and α7 in 8OYQ.

Each entry also picks its own page axes, so the same CheA domain is drawn rotated: α3 horizontal above the sheet
in 8P3N, vertical beside it in 8OYQ.

**Checks:** A19 passes both pairs, by the A19 ruling: names are compared only up to the first unpaired element.
V9 ("names and order match, drawings differ only where 3D differs") fails.

**Ruling needed:** V9 against the A19 ruling. If V9 stands, the fix is a pair mode:
- name the second entry's elements by the first's (residue-range pairing, as A19 does);
- draw it with the first's page axes after superposition.

**Effort:** M.

### C20 (c, visual, ruling): additives listed in the facts, hidden from the figure with no note (3 entries)

**Entries:** 5Q42 (MLI), 9B7C (GOL), 7LX9 (TRS).

**Checks:** none.

**Mechanism:** `io.ADDITIVES` (io.py:582) is hidden by `ligand_marks`' `not g.additive`. The contact sheet's
facts list every ligand. Neither the figure nor the CLI says an additive was left out. The reviewers fail V1/V7
on it.

**Ruling needed:** either additives are not "ligands the facts list" (the facts line filters them, the reviewer
checklist reads V1 that way), or the CLI note says "hidden as additives: GOL (--ligands all shows them)".

**Effort:** S. **Bench:** none. 9B7C passes the visual review otherwise.

### C21 (a): domain panels of one chain overlap (3 entries)

**Entries:** 2ILU (D1·A3/D4·A3 0.0675, D1·A4/D4·A4 0.0688), 5X4Z (D4·A/D9·B 0.026), 7JU4 (D1·7/D3·7 0.060).

**Checks:** A16 `panel_overlaps` > 0.02.

**Mechanism:** `panels.py` splits interleaved space by the nearest element. Where one domain's element sits
between two of another's, the outlines share a corner. 2ILU's D1/D4 is a Rossmann + oligomerisation-domain
interleave.

**Fix:** after outlining, shrink a panel locally where two overlap (clip by the bisector between the two
domains' nearest elements). **Effort:** M. **Bench:** domain views (1LDM, 8UTF, 4OO8 snapshots).

**Batch B result:** check defect. A panel with another domain's element interleaved in it has a hole (an inner ring
turning the other way); Agg fills by the nonzero rule, so the hole is empty on the page, but `figcheck._cells` used
`Path.contains_points`, which tests each ring alone and counted the hole as inside. Cells are now counted by nonzero
winding (the panel fields already share no cell, asserted on 3L9G). 5X4Z, 7JU4 and 3L9G overlaps -> none; 2ILU's went
with C31.

### C22 (b): copies with fewer than 3 elements judged inf (3 entries)

**Entries:** 7AOD (its only A15 failure: `D~P` inf, 2 of 2 matched, a two-helix chain), 6JLN (`f~F`, `i~I`,
`j~J`, `k~K` inf with 1 of 1: single-TM subunits), 9P1X.

**Checks:** A15.

**Mechanism:** `figcheck._copy_pairs` (figcheck.py:561) requires `COPY_MIN_MATCHED = 3` and half the smaller
copy, else inf (ruling Imp3). A chain with one or two elements can never meet it.

**Fix:** when the *smaller copy itself* has < 3 elements and all of them match, compare those (2 points fit any
two points by rotation, so the deviation is the difference of their separations). Or report "n/a (k
elements)".

**Effort:** S. **Bench:** none. Recovers 7AOD's A15. 7AOD still fails C1/C3/C4/C6/C5.

### C23 (a): automatic domains span chains (3 entries)

**Entries:** 8RRX (D19), 8QEM (D1·A, D3·N, D4·N), 7BBH (D4, D6).

**Checks:** A16 `cross_chain`.

**Mechanism:**
- `features.auto_domains` (features.py:166–206) assigns a copy's elements the domain of their mate in the
  reference copy. `owner[s.id] = near if chain in lent else (chain, near[1])` (features.py:195): a lent chain's
  elements join the mate's (chain, d) and so another chain's domain.
- `_fold_small_copies` (features.py:197) then merges small copies into a neighbour's domain without checking the
  chain.
- The name has no " · chain" when only one chain is split (8RRX D19, 7BBH D4), so the cross-chain domain is easy
  to miss on the page.

**Fix:** `_fold_small_copies` and the lent branch keep the chain: a lent element starts or joins a domain of its
own chain. **Effort:** S–M. **Bench:** none (domains are a view).

**Batch B result:** lent elements (focus mode, drawn grey and primed) are no longer members of a domain
(`Layout.domains`); `Layout.domain_lent` keeps them per domain so the panel, its name field and the sheet tint are
drawn round them as before (`layout.panel_members`; 1LDM protomer figure identical). 8RRX, 7BBH and 8QEM cross_chain
-> none. Making them a domain of their own chain was tried: one-strand domains, panels torn round the grey strands.

### C24 (a): NA route fallbacks (3 entries)

**Entries:** 9C75 (`strand R:28-29`), 7Y7C (`gap after a:1439`, Z:27/37/61), 8VWQ (strand 1:662-663, 934-935,
943-945, 1506-1507).

**Checks:** A3; A17 `fallbacks`.

**Mechanism:** `route_na.route_nucleic` gives up on a 2-nt run between stems placed tight together (9C75), and
on gaps across the 4,000-nt crops. **Effort:** L for the large ones, M for 9C75. **Bench:** NA entries; NA
campaign rerun.

**Batch B result:** a joint between two pieces of a bent duplex that no route reaches is drawn straight up to
`na_runs._JOINT_REACH` (2.5) under the same tests (forward, clear by _INFLATE, crossing no strand but its partner):
9C75 P2b -> P2a (1.83 apart) no longer falls back. 7Y7C and 8VWQ (gaps across 4,000-nt crops) are not this case and
stay with C4/C8 (batch C).

### C25 (a, visual): neighbour ghosts not deduplicated; focus picks an uninformative chain (3 entries)

- **8QEM default:** each of the 11 symmetry mates' α11 is drawn as its own labelled ghost "α11‴" in a row. The
  domains view even gives each a D4·N hull.
- **7JU4 focus:** the view picks chain "0", a two-helix coiled coil, apparently by sort order.
- **7ARD focus:** a small 4-helix subunit with 9 neighbour helices.

**Mechanism:**
- Ghosts come from `_symmetry_ties`/`block_ghosts` per mate, with no merge of mates that lend the same element.
- Focus chooses the first subunit when no symmetry ranks them.

**Fix:** merge ghosts of one element class lent by several mates into one "×11" ghost. Focus picks the largest or
most-contacted chain. **Effort:** M.

**Batch B result:** a focused subunit borrows the helices packed directly along its own
(`features.helix_partners`), not every helix of the bundle: 8QEM's ring of twelve C-terminal helices lends two, not
eleven (15 lent elements -> 6). The campaign's focus view takes the protein chain with the most elements drawn in
the default view (`random100.focus_target`), not the first (7JU4's chain "0", a two-helix coiled coil).

### C26 (a/c): loops run inside the membrane slab (2 entries)

**Entries:** 7SMR (B:117-129>B:141-150, B:171-173>B:190-201, E:156-160>E:162-163, E:166-167>E:184-195), 8W9Y
(C/F:225-249>…:271-280).

**Checks:** A18 `loops in the slab`.

**Mechanism:** the router's slab cells (route.py:88–91) are a soft cost. A loop between two elements whose ports
are both inside the slab (a short re-entrant or amphipathic helix, 8W9Y C:225-249) has no outside path shorter
than the slab one.

**Ruling or fix:** decide whether a loop joining two in-slab ports counts (it is in the slab in 3D too), or
route such loops out to the nearer face. **Effort:** M.

**Batch B result:** a step along the slab costs `route._SLAB_RUN` (12) where a step across it costs `_SLAB_COST`
(3), so a loop leaves the slab at a face instead of running along it over an element that reaches into the band
(7SMR B:117-129). A18's port rule measures the port against the band itself, not its 0.3 margin (7SMR E:162-163 and
E:166-167 end 0.29 inside the face). 7SMR: 4 loops in the slab -> 1 (B:171-173>B:190-201 still dips 0.4 into the
band over an element end at the face). At 24 that loop went out too, but another of chain B's loops then crossed the
whole slab twice to run along the far face: worse, so 12 was kept. 8W9Y keeps C/F:225-249>271-280: the loop leaves an
interfacial helix's end inside the slab with α8 right above it; the way out through the near face (tried at cost 30)
is a long detour that strikes the α7 and α10 names. Snapshot re-recorded: 1C3W-membrane (one A1 loop above the band
takes the other side of a strand end; crossings 6 -> 3).

### C27 (b/c, ruling): A18 requires 75 % band cover from TM helices that are tilted or offset (2 entries, 1 last)

**Entries:** 6PS3 (β2AR, TM A:267-298 covers 0.73 of the band), 7PIU (R:48-70 0.70).

**Mechanism:**
- The band is the median face crossing over all TM elements (`cli._face_medians`).
- A TM helix crossing at a 30° tilt projects about 0.87 of its height, and its ends sit off the medians.
- `random100.membrane_checks` asks for ≥ 0.75 of the band's height.
- In both entries the helix visibly spans the band with one end a little short.

**Ruling:** 0.70, or "reaches within 0.5 units of both faces". The alternative fix is to snap TM-helix ends to
the faces in the membrane view. **Effort:** S.

**Batch B result (ruling C27):** A18 passes a TM helix that spans the band through both leaflets: its drawn extent
reaches past the band's middle from a face (6PS3 A:267-298, 0.73; 7PIU R:48-70, 0.70). A helix in one leaflet or off
the band still fails. Read of "spans the band (enters both faces)": the strict reading, both face lines reached,
fails nearly every TM helix of the passing membrane entries (5JDN, 7SMR, 8W9Y: 0.3-1.6 units short), so it was not
taken; the controller should confirm.

### C28 (a): glycan stems and linked-residue labels over labels (2 entries)

**Entries:** 1WWL domains (`label:A:259-261 × glycan:254:stem`), 8K8E (`label:H:10-12 × linked-label:H:15`).
Glycan hits on pages over 44 units (7SMR, 8OZJ, 4UM3, 7BBH) are counted under C4.

**Mechanism:**
- Glycan stems: `render` (render.py:1700–1760) avoids `label_boxes` and loops, but scores only overlap count and
  takes the least-bad lead. A short helix label in the way is accepted.
- `linked_marks` (render.py:1340) puts the linked residue's label at `text_at` without testing helix labels.

**Effort:** S. **Bench:** glycan snapshots (6ZS5).

**Batch C2a result (7SMR, label x glycan stem):** the glycan at D:1358 hangs on strand G mid-sheet; left and right
are strands, up are loops, so its stem ran down the strand over the letter G (the "fewest crossings" fallback finds
no way without one). A glycan stem is now drawn as a ligand's tether is, only where it shows: off the halo of every
text (`render.visible_pieces` with `text_halos`). 7SMR's A6 hit -> 0. Bench identical; no snapshot changed.

### C29 (a/c, visual): G-quadruplex tetrads not drawn; NA letters detached from the trace (2 entries)

- **8Q4O (RNA G-quadruplex, NMR):** no canonical pairs, so no stems. The figure is a bare projected P trace that
  crosses itself about 4 times, with letters 30–60 px off it and two "C"s struck by the backbone. Published
  figures draw stacked G-tetrad planes.
- **8XZE:** at the P1c/P1b junction six nucleotides fan round a 50-px backbone jog out of reading order, and two
  "A"s float inside the ladder.

**Mechanism:** an unpaired trace places letters by projected offset (`route_na._place_letters`, route_na.py:187)
with no backbone obstacle test. Hoogsteen tetrads are not read as pairs (as C16's 8BW5).

**Ruling or fix:** a tetrad glyph (M–L), or a note "G-quadruplex: tetrads not drawn". Letter placement on the
trace with a backbone keep-out (M).

**Batch B result (ruling C29):** `nucleic.g_tetrads` reads G-tetrads (four guanines, each Hoogsteen-bonded to the
next: N1-O6 and N2-N7 within 3.5 Å); `scope_note` says "not drawn: 3 G-tetrads of a G-quadruplex (their Gs drawn
unpaired)" (8Q4O: 3; 8BW5: 2). 8XZE's junction letters detached from the trace (the second half of C29) are not
done.

### C30 (a): `--assembly auto` keeps a 40-chain ASU (1 entry)

**Entry:** 4UM3. **Checks:** A2 (301 s), A2b (app over 90 s), A4, A6, A16. 4UM3 also has C4/C6/C5 hits.

**Mechanism:**
- 4UM3 is eight pentamers in the ASU. Assembly 1 is one pentamer (P–T): identity operators.
- `io._one_complex_of_several` (io.py:228–248) requires every entity of the assembly to be deposited at least
  twice as often as it holds it. Entities 2 and 3 are each one chain in the whole ASU: 38 chains are entity 1,
  and R and S are their own entities, probably a modified or alternate copy. So the test fails.
- `read_model` (io.py:290–299) then draws the ASU (copies 5 ≤ 40 deposited).
- The 40 chains are untied (copy_ties 0), so copies are drawn in two mirror forms (visual).

**Fix:**
- Compare by polymer *sequence* (or count the assembly's chains: the ASU holds ≥ 2× the assembly's polymer
  chains and every assembly chain's sequence occurs ≥ 2× in the ASU), not by entity id.
- Fixture: a 2-pentamer cut of 4UM3 (< 100 KB, CA+backbone).

**Effort:** S. **Bench:** none (no bench entry has this shape). Expected to turn A2/A2b/A16 and most of 4UM3's
A6 hits.

**Batch B (NA campaign):** batch A's fix round 2 (`io.PAIRED_MIN` 2, base pairs only) was not followed by an NA
run; 4PLX's three copies touch in the crystal by one swapped pair of bases (2 residue pairs), were made one object and
failed (na_crossings 1). PAIRED_MIN is 3 (3OWW pairs in 4).

### C31 (a): the domains view explodes the page when domains interlock across copies (1 entry)

**Entry:** 3L9G (urate oxidase tetramer). The domains view is 2160×292 px, area 174,652 against 1,599 default.
A1+A4 sit at the far left and A2+A3 at the far right, joined by about 2000-px loops. The D1/D2/D3 legend is off
the page.

**Checks:** A3, A4, A6 (1,393), A7, A16 (120 name collisions).

**Mechanism:**
- `layout._part_domains` (layout.py:2535–2560) parts "two domains of one unit, copies of each other … by
  widening the spread between protomers (protomer k moves k times the step)".
- 3L9G's copies share sheets (C10: 12 of 20 elements of each tie are in cross-chain sheets). Domains of different
  copies are one item and can never part.
- So every one of the 60 sweeps widens the spread again.

**Fix:**
- Stop widening when an overlap does not shrink between sweeps.
- Treat domains that share an item as one (do not part them).

**Fixture:** 3L9G cut to 2 chains. **Effort:** M. **Bench:** domain snapshots unchanged (no interlocked copies
there), to be measured.

**Batch B result:** when `_part_domains`' sweeps leave the domains unparted and the items over twice as wide and high,
it starts over with domains that share an item (a sheet shared by copies) never parted, stopping when a sweep no
longer shrinks the overlap (joining them always changed 1WWL, 2ILU, 7SMR, 6JLN and 5LDW, which parted fine). 3L9G's domains view: 2160 x 292 px -> about the
default's size, no fallback, the legend on the page.

### C32 (c/a): a 61-chain default view over 180 s (1 entry)

**Entry:** 7JU4 (192.7 s, limit 180 s for > 20 chains). **Check:** A2.

**Mechanism:** all 61 chains are drawn (assembly 1 = all chains, no symmetry); routing about 1,000 loops
dominates. **Ruling or fix:** a profile first. 7JU4's figure is unusable anyway (C4/C8/C7), so a focus default
above about 40 chains of no symmetry would fix both. **Effort:** M.

**Batch C1 result (C32; the time cause the Task 12 rulings call C33):** profiled, the slow views spent their time in
work done twice or pair by pair, not in drawing. Removed, with the output the same (SVG, layout file, figcheck
findings and NA metrics byte-identical on 13 views: 29RY, 6JO6, 6DPV, 3OWW, 8BW5, 1KX5 and the domains views of 2C4M,
3L9G-A, 5X4Z, 9GEP, 6JLN, 8RVL, 2ILU; both bench baselines identical):
- `plot --save-layout` laid the figure out a second time for the layout file (27 of 7JU4's 183 s); with no edits it
  is written from the figure's own layout (`spec.layout_document(..., drawn=)`).
- `figcheck.text_collisions` measured every text against every body and line (`path.get_extents` 2 M times):
  each body is measured once, and only pairs whose boxes meet are tested (90 s -> 0.5 s on 7JU4).
  `loop_crossings`, `honest_crossings` and `nacheck.render_checks`' letter and bead overlaps likewise.
- `layout._overlapping` measured all n x n footprint pairs each sweep; past 256 footprints only pairs near enough
  across the page (same arithmetic, same pairs and order).
- `layout._part_domains` measured every item's footprint four times for every pair of domains in every sweep
  (90,032 footprints on 3L9G-A for 34 items): once now.
- `panels.panel_fields` held every domain's whole raster and took each domain's rivals one by one: the nearest and
  second-nearest domain per cell are kept once, and each domain is worked on its own raster window.

Times of the campaign's view process (`random100.py view`, eight run side by side, same load before and after):
7JU4 default 202 -> 52 s, 7Y7C default 193 -> 82 s, 8VWQ default 194 -> 76 s, 7JU4 domains 632 -> 225 s (its
figure alone: 20.3 GB -> 6.2 GB). What is left is foldmap's own routing (A* on 7JU4's 61-chain domains page: 106 s
profiled) and NA letter placement (7Y7C).

### C33 (a, visual): ligands that touch only nucleic acid are never drawn, silently (1 entry)

**Entry:** 8XZE (THF-II riboswitch: the facts list SPM and THG; neither is drawn, and nothing is said).

**Checks:** none. No criterion compares the facts' ligands (less additives) with the marks drawn.

**Mechanism:** a ligand's `contacts` are protein residues of the backbone. `render._tether_targets`
(render.py:1170–1190) finds no target for a ligand that touches only nucleotides, and `ligand_marks` skips it
(`if not targets: continue`, render.py:1241). The same holds for a nucleotide-only pocket in a protein–RNA entry.

**Fix:**
- Tether ligands to the nucleotides they touch (`residue_point` for NA letters, or the stem's rung).
- Add a check: every non-additive ligand in the facts has a mark, else the note says why.

**Fixture:** 8XZE (1 RNA chain, small). **Effort:** M. **Bench:** none; NA snapshots with ligands, if any.

**Batch C2a result (C33):** `io._ligands` keeps a small molecule touching nucleotides but no protein residue, with
those nucleotides (`Ligand.nucleotides`); it is marked (`render._na_targets`: one tether per stem or strand it sits
on, up to four) and kept off the nucleic acid's stems, paths, letters and names (`_na_spots`). Metal ions held by
nucleotides alone are counted (`Links.na_ions`), not marked (7Y7C holds 283 Mg, 8VWQ 586), and the CLI, MCP and app
notes name them ("not drawn: metal ions held by the nucleic acid alone (MG)"). Drawn tethers now keep off stem names'
margin and the nucleic acid's letters, which are drawn after them (2GIS's SAM tether by P1, 8XZE's THG by A21).
Iridium and cobalt hexammine (IRI, NCO: RNA crystallisation and phasing agents) are additives. 8XZE passes in every
view. Figures that gain marks in the NA campaign (metrics identical): 1Y26 ADE, 2GDI TPP, 2GIS SAM, 2NZ4 GLP, 3D0U
LYS, 4LVV FFO, 4YBB PUT and SPD; their notes name the ions. Bench identical; no snapshot changed.

### C34 (a, visual): a bulged strand drawn as its own column between paired strands (1 entry)

**Entry:** 5KQJ. The sheet is drawn C D E F B A. B–E are paired (4.2 Å), but F is wedged between them. F pairs with
E (4.3 Å) and sits 7.1 Å from B: F and B share E's free side end to end (F at 3D y = +11.8, B at −1.1).

**Mechanism:** `sheets.build_sheets` (sheets.py:70) orders each sheet as one path of strands. Two strands that
both pair with the same side of a third (a β-bulge or end-to-end stack, COMPARE §3.1) get consecutive columns.

**Fix:** detect two strands on one side of a strand with no residue overlap along it, and stack them in one
column. **Effort:** M. **Bench:** sheet entries where it occurs (to be measured).

**Batch B result (C34):** not done. A sheet block is one row of columns (`layout._sheet_item`, offsets per slot,
all centred vertically); putting two end-to-end strands in one column needs per-strand vertical offsets, A11's
order check and the sheet plates to follow. Larger than an M change; left for a later batch.

**Batch B copy defects:**
- **2MJW B~D 3.34 -> 0:** `_symmetry_ties` matched copies by shared residue numbers; chain D is numbered 95-188 where
  B is 1-94, so they shared none and were never tied. Numbers are taken to the first copy's by homolog_chains'
  offset.
- **7OKA 1.12 -> 0:** DSSP broke helix 225-238 at 232 in B and C, not in A; each piece was tied to A's whole helix and
  took its length. Where DSSP breaks a helix by one or two residues in some copies, the majority decides
  (`ss._copy_breaks`; a tie cuts): the whole copies are cut there, or the broken ones joined. 7TIM's B:178-203 is now
  broken like A's (both draw α8 and α9); 5JB1, 7JU4, 8QEM and 9P1X join a few broken pieces (fewer elements).

### C35 (a, batch C1 review): chain colours repeat past the palette (3 entries seen)

**Entries:** 5X4Z (A = I = `#0072b2`, B = J, C = K, D = L), 6JO6 (three passes), 7Y7C. **Checks:** none (visual).

**Mechanism:** `palette.chain_colors` cycled `cols[i % len(cols)]`.

**Batch C1 fix round 1 result:** from the second pass each chain takes a lightness step of the colour it repeats
(tint 0.45, darken 0.75, tint 0.7, ... as the pass's own step first), the first that stands CIEDE2000 >= 10 from every
colour already assigned, the white page and its family's helix back faces (`palette.BACK_SHADE`); else the step of
its L* that stands furthest. 24 chains on okabe-ito: pairwise >= 10.2; 61 chains: all distinct (min 4.1, the
palette's eight hues allow no more). The first eight are the palette's: figures of up to eight chains are unchanged
(no snapshot has more; none re-recorded; bench identical).

### C36 (a, batch C2a): `--names beta` writes a run of primes past ‴ (1 entry)

**Entry:** 9YKJ (A4 in the beta view: A:194-198>A:204-210 and A:664-674>A:683-696 drawn over labels). Chain A has
nine helices after β10 (and chain B nine after β8). `names._primed` gives αn, αn′, αn″, αn‴ and then "′" × k, so the
beta view wrote α10′′′′′′′′′ (nine primes): names up to three times as wide as α10, whose label boxes left the loops
no route.

**Batch C2a result:** when a strand is followed by more helices than αn … αn‴ can name, `names._row_helices` letters
them in sequence instead, αna, αnb … (the row helix keeps αn); four or fewer keep their primes, so 1TIM's α4′ and
α5′ stay. 9YKJ passes in every view. README: `--names beta` paragraph. Bench identical (names are opt-in).

## Visual review

**Who reviewed:**
- Three reviewers read the contact sheet and the full-size default figure of 34 entries: the 24 passing and the
  10 worst by check count.
- Each answered V1–V9 (constraints, "Visual review: what 'pass' means").
- Batch 3's three NA and membrane sheets were cross-checked by the triager (legend crops of 6JLN and 7QUP focus,
  6SO9, 9YLM, 7IMM, 9PXS, 29RY).
- The calibration re-review (10 seeded entries) was not run in this task.

**Strictness:**
- The reviewers were strict. Only 4 of the 24 check-passing entries pass every item: 3RVC, 1MG8, 8B4R and
  AF_AFB1Q234F1.
- The frequent "no" answers are V2:
  - parallel loops (C7);
  - ligand marks and tethers on lines and bodies (C5, C9);
  - loops crossing at an arrow tip.
- V8 failed on more than two edits.
- These are real and are now causes C9 and C7. A ruling should say whether V2's "two lines run together" is
  gated on its own or only through V8.

### Passing entries (24)

Batch 1:
- **5Q42:** V1 n (MLI hidden, C20); V2 n (three horizontals under the K/O/Q tips 18 px apart; M/L tips; Ni tether
  about 400 px over B and A's head, C9/C7); V3 y; V4 n/a; V5 y; V6 n/a; V7 n (C20); V8 n (37 crossings; α3/α8/η3/η4
  loops cut both lower sheets); V9 n/a. **FAIL**
- **5HI7:** V1 y; V2 n (three verticals 17 px apart by α9/α10; two loops together under α13/α14; 62X tether over
  strand A; Zn tether over α8: C9, C7); V3 y; V4 n/a; V5 y (the 62X tethers reach η2 about 330 px away); V6 n/a;
  V7 y; V8 n (15 crossings, tangled α13–α15 routing). **FAIL**
- **3RVC:** V1–V8 y (edits: move α3 so the G→α3 loop clears α2's tip; move the C→D loop off E's tip); V9 n/a.
  **PASS**
- **1F09:** V1 y (CR2-66 star); V2 n (I→J, G→H, B→η2, η2→C, N→A run 15–20 px apart under the sheet, C7; the star
  touches H's dashed closure); V3 y (cycle 8-9-4-5-6-1-2-3-11-10-7 correct); V4 n/a; V5 y; V6 n/a; V7 y; V8 n
  (central helix above the sheet, C13). **FAIL**
- **9C0M:** V1 y; V2 n (three Ca beads on loop lines, C5; a Ca tether along the G loop for about 280 px, C9; two
  loops together for about 600 px, C7; α7→α6 and α5 loops ambiguous at α5's tip); V3 y; V4 n/a; V5 y; V6 n/a; V7 y;
  V8 n (22 crossings). **FAIL**
- **9B7C:** V1 n (GOL hidden, C20); V2 y; V3 y; V4 n/a; V5 n/a; V6 n/a; V7 n (C20); V8 y (rename α1–α4 to A–D,
  optional). **FAIL on C20 only**
- **7LX9:** V1 n (TRS hidden, C20); V2 n (YGM→α6 tether on the dashed α5→α6 gap loop for about 230 px, C9); V3 y;
  V4 n/a; V5 y; V6 n/a; V7 n (C20); V8 y (move YGM; rename A–J). **FAIL**
- **7FDN:** V1 y; V2 n (chain B: two loops together from x≈1475–1640; α3→C crosses the dashed gap at B's tip); V3 y
  (unsure); V4 n (elements alike, loops routed differently in A and B, C12); V5 n/a; V6 n/a; V7 y; V8 n. **FAIL**
- **2I5Y:** V1 y; V2 n (a NAG stalk about 290 px across J/O/M/L, C9; dashed gaps and loops together under the
  gp120 sheet; "10 aa" unexplained); V3 y; V4 n (VL beside CH1, C13); V5 n (NAG marks unlabelled, two floating
  at the page bottom); V7 n (no page note for "one of 2 copies in the file"; BIF/DPR/VLM unmarked); V8 n (60
  crossings). **FAIL**
- **1A0Z:** V1 y; V2 n (HEM of B and D on the α2 loop, C5; B's HEM→α7 tether over α5, C9); V3 n (α chains: G/H
  450 px apart, C13); V4 y; V5 y; V6 n/a; V7 y; V8 n (move HEM; bring α8 to α7; rename A–H). **FAIL**
- **1BXI:** V1 y; V2 n (Ni on the α1→B loop, C5; Ni→α6 tether over D, C9; loops together by A and η3, C7); V3 y;
  V4 y; V5 y; V6 n/a; V7 y; V8 y (move Ni; move α6). **FAIL (V2)**
- **2MJW:** V1 y; V2 y; V3 y; V4 n (S100P B and D drawn completely differently: not tied, C18); V5–V7 y/n/a; V8 n
  (redraw D as B). **FAIL**

Batch 2 (positions in pixels of the full-size default PNG; arrangement checked by projecting 3D centroids onto
the page axes; pairing by minimum CA–CA distances):
- **8Q4O:** V1 y; V2 n (two "C"s struck by the backbone; letters 30–60 px off the trace, which crosses itself about
  4×, C29); V3 y; V6 n (no G-tetrads, C29); V7 n (non-canonical pairing left out, unsaid); V8 n. **FAIL**
- **1MG8:** V1–V8 y (edits: rename A–D to β1/β2/β3/β5, which also removes the clash of the C-terminus letter with
  strand C; move α1 onto the sheet face); V9 n/a. **PASS**
- **5KQJ:** V1 y; V2 y; V3 n (F wedged between paired E and B, C34; α1 right of the sheet but left/under C/D in 3D,
  C13); V7 y; V8 n (also the α2→B loop detours about 1,600 px round the page bottom). **FAIL**
- **8OH7:** V1 y; V2 n (α2→α3→α4 loops parallel for about 190 px into a hooked corner; α15→α17 loops together for
  about 100 px; α12's return loop hugs its edge, C7); V3 y; V8 n (repeats drift about 360 px down, C13). **FAIL**
- **8B4R:** V1–V8 y (edits: merge the two sheet plates of the one kinked hairpin; straighten the C→D jog). **PASS**
- **3MHT:** V1 y; V2 n (SAH box on the η3 loop, C5; "12 aa" gap count between two runs, C3; α1/α2 loops parallel for
  about 400 px, C7; knot round η1/η2); V3 n (α8 above α4, lowest in 3D, C13); V5 y (one SAH tether about 500 px,
  C9); V6 y; V7 y; V8 n. **FAIL**
- **8XZE:** V1 n (SPM, THG not drawn, C33); V2 n (junction nucleotides fanned out of order; "A"s and "G A" floating
  inside the ladder, C29; all letters turned about 90°); V3 y; V5 n (C33); V6 y; V7 n (C33); V8 n. **FAIL**
- **10UJ:** V1 y; V2 n (3–4 loops 10–15 px apart under sheet A J E B D C; bundles round α15; 83 crossings, C7); V3 y
  (rank correlation x 0.84, y 0.92); V8 n; V9 n (names and layout not aligned with 9YKJ, C19). **FAIL**
- **8P3N:** V1 y; V2 n (ADZ on loops in both copies, C5; ADZ at 3.5 pt); V3 n (copy A translated from B though
  stacked behind it: left/right reversed, C13); V4 n (ADZ placed differently in the two copies; α5→α6 differs by a
  real DSSP difference); V8 n; V9 n (C19). **FAIL**
- **8OYQ:** V1 y; V2 n (α2→A and α3→α4 links parallel about 160 px then share one vertical next to W6T, C7; W6T
  over the α4→B link, C5); V3 y; V5 y; V7 y; V8 y (two moves); V9 n (C19). **FAIL**
- **AF_AFO42816F1:** V1 y; V2 n (two loops together for about 450 px down the sheet's right edge, past "N"; more
  parallels round α15/η1, C7); V3 y (x 0.96, y 0.96); V8 n (27 crossings; wants the βα-unit layout, COMPARE §3.3).
  **FAIL**
- **AF_AFB1Q234F1:** V1–V8 y (edit: straighten the α1→α2 jog). **PASS**

**Tally of the 24 check-passing entries:** 4 pass every item: 3RVC, 1MG8, 8B4R and AF_AFB1Q234F1. 9B7C passes
apart from C20. The other 19 fail on the visual causes above, mostly C7 (loops together), C5 (marks on loops)
and C9 (tethers).

### Worst 10 by check count

- **3L9G:** V1 y; V2 n (4–6 connectors side by side along both margins and between the sheets, C7; MRD/URC
  tethers across both sheets, C9); V3 y; V4 n (copies' sheets pool foreign strands: the C10 shared sheets, with
  an unexplained ghost "I" and "D"); V5 n (URC/MRD 9–16 units from their contacts, C2 reported lengths); V6 n/a;
  V7 y; V8 n. Domains view: the C31 explosion. **FAIL**
- **4UM3:** V1 y; V2 n (strand letters fuse into words over 8-px strands, C4; two legend rows overlap, C4); V3 n (a
  40-chain ASU mosaic, C30); V4 n (two mirror forms with no stated mirror: copies untied, C30); V5 n; V6 n/a; V7 y;
  V8 n. Focus view (one jelly-roll chain) close to publishable. **FAIL**
- **7JU4:** V1 n (N/C labels about 5 px, termini inside the bundle); V2 n (about 10 parallel connectors down the
  centre, C7; fused letters, C4); V3 n; V4 n (KJLIHGFAD sheets in both hands); V5 n (GDP/GTP tethers 10–23 units);
  V7 y; V8 n (unusable; domains view a hairball). **FAIL**
- **7ARD:** V1 n; V2 n (fused strand words, C4; ≥ 10 concentric loops round a group, C8/C7; legend rows overlap,
  C4); V3 y (coarse); V4 n (hetero chains interleaved); V5 n (not detected as membrane, C17); V7 y; V8 n. **FAIL**
- **6JLN:** V1 y; V2 n (about 100 dashed tethers over the TM region, C9; CLA/BCR tags about 6 px, C4; legend rows
  overlap, C4); V3 y; V4 n (the two monomers differ, C12); V5 n (band half the helix length, single-TM subunits
  above it, C17); V7 y; V8 n. **FAIL**
- **7QUP:** default (α/β-tubulin) clean: V1 y; V2 n (GTP tether about 400 px across α5/α9, C9; three loops
  together by 7B α10/α11, C7); V3 y; V4 y; V5 y; V7 n (reviewer: the note "1 of 33 subunits" does not say why 65
  chains make 33; the CLI note does say "--focus none draws all of them"; see the ruling list); V8 y (two moves).
  Focus none: one tubulin in two forms, C12. **FAIL**
- **8QEM:** V1 y; V2 n ("35 aa" under N's terminus, C1); V3 y; V4 y; V6 n (DNA Y/Z not in the default; the note
  says "2 of 26 chains … --focus none draws all", not that DNA is left out: ruling); V7 n (same); V8 n (11 ghost
  α11‴, C25). **FAIL**
- **7Y7C:** V1 y; V2 n (stem blocks overlap, rotated range labels on blocks: C4/C14; fused letters; 4–5 legend rows
  overlap, C4; purple connectors in perimeter bundles, C7); V3 n (no recognisable LSU domain layout); V6 n
  (contact beads about 3 px; blocks overlap); V7 n (crop said only in the facts/CLI: ruling); V8 n. **FAIL**
- **8VWQ:** as 7Y7C, and stem names `na:label:189`/`221` 0.000 from their letters; protein groups stranded
  outside the RNA ("EDABC", "DCBA"). **FAIL**
- **5LDW:** V1 y; V2 n (SF4 tags over strand letters "CEF"; "44 aa" over strand B, C1); V3 y; V4 y; V5 n (TM
  helices wrapped into a row above the band, C17); V7 y; V8 n. **FAIL**

### Published-figure comparison (COMPARE.md folds)

- **1F09 (GFP, §2.4):** the barrel order matches Bignon 2023 up to the cut point (8-9-4-5-6-1-2-3-11-10-7 against
  1-2-3-11-10-7-8-9-4-5-6). The §2.4 chromophore bug is fixed: CR2 kept, chain unbroken, starred glyph. The
  central helix is still not drawn through the barrel, and the tangle under the sheet is worse than the
  reference (C13, C7).
- **7LX9 (T4 lysozyme, §2.7):** much closer to Patzer 2012 than 2LZM was: N-lobe top-left, α3 horizontal across
  the middle. The C-lobe splits into a lower group and a right column, not one row; names are α1–α9, not A–J.
- **1A0Z (haemoglobin, §2.8):** both §2.8 problems are fixed: A/B split at the kink, HEM with contacts to E and F.
  Names α1–α8, not A–H. The α chains sprawl (G/H apart) where the β chains are compact (C13).
- **1MG8 (ubiquitin-like β-grasp, §2.2):** close to Burroughs and Pro-origami. The order β2-β1-β5-β3 and the
  parallel β1/β5 pair are right. As in §2.2, α1 sits beside the sheet, not across its face, and strands are
  lettered A–D, not β1–β5.
- **8XZE (riboswitch, §2.13):** the hairpin-rod drawing matches the planar convention better than the 2GIS
  projection did: no backbone crossings, stem names within one stem width. Published riboswitch figures always
  show the ligand in its binding loop; here THF and spermine are missing (C33), and the junction nucleotides are
  not in a readable order (C29).
- **8Q4O (G-quadruplex):** published figures draw stacked tetrad planes linked by loops. foldmap gives an
  unpaired projected trace with floating letters. Below the convention, and §2.13/2.14 offer no fallback (C29).
- **Rossmann (2ILU, 2HKI/1RM5, §2.5) and RRM (2DGS, 6SO9, §2.6):**
  - 2ILU's default and 2HKI's figure keep the Rossmann sheet in 3D order. Their failures are C10/C3/C18, not the
    fold.
  - 6SO9 (RRM + RNA) draws the β1-α1-β2-β3-α2-β4 sheet (order D E A C B, 4-1-3-2 with β-hairpin D) with the RNA
    over the sheet face, as published RRM–RNA figures do. Its failure is C11 (loops over the RNA letters).

### Rulings needed

1. **C10:** leave strands of shared (intermolecular) sheets out of A15, or require a stated mirror.
2. **C2:** confirm A20's implementation of the ruling (other-chain tethers only), and keep same-chain lengths
   as reported.
3. **C20:** additives in the facts against the figure (V1/V7).
4. **C27:** A18 tolerance for tilted TM helices.
5. **C26:** loops between two in-slab ports.
6. **C16:** 8BW5 quadruplex pairs in pair_f1; 9YLM like M1 or a layout fix.
7. **C22:** copies with < 3 elements.
8. **Scope on the page:**
   - V7 for 7QUP, 8QEM, 7Y7C and 8VWQ: the CLI note says what was left out and how to draw it, but the figure
     carries no note. Is the CLI note enough (it was for plan 11)?
   - 8QEM's dropped DNA is covered only by the chain count.
9. **C7:** whether V2's "two lines run together" fails an entry on its own.
10. **C4:** maximum figure size (a poster mode past 7.2 in) against em-scaled spacing.
11. **C8/C32:** default to a focus view above about 30–40 untied chains.
12. **C19:** V9 "names match" against the A19 ruling (generic numbering may differ after an extra element).
13. **C29:** draw G-tetrads, or say they are not drawn.
14. **C33:** a criterion that every non-additive ligand in the facts is drawn or explained.

## Per-entry triage (76 failing entries)

Format: `<ID>: <criterion> <value> → <cause>`. For A6 the value is the hit count over the figure views, and the
causes are given with their hit counts. Visual-only findings for an entry are in the visual review above.

- **13MV** (S1, 1 ch): A20 ligand-tether:A1CUK:A402:3 7.644 → C2
- **29RY** (S1, 1 ch): A6 3 hits (beta 1, default 1, deposited 1) → C1 3
- **2YBE** (S1, 1 ch): A20 ligand-tether:ALF:A1420:1 9.283; ligand-tether:LA8:A1419:3 12.492 → C2
- **6IDJ** (S1, 1 ch): A20 ligand-tether:9AU:A406:3 11.544; ligand-tether:FMN:A404:0 6.983; ligand-tether:FMN:A404:1 … → C2
- **7IMM** (S1, 1 ch): A15 default 5.626768267951088 → C10
- **2B9V** (S2, 16 ch): A6 11 hits (beta 1, default 1, deposited 1, focus 8) → C3 7, C1 4
- **2C4M** (S2, 4 ch): A6 13 hits (domains 13) → C4 9, C6 3, C5 1; A20 ligand-tether:PLP:A1634:1 12.241; ligand-tether:PLP:A1634:2 6.755; ligand-tether:PLP:B1634… → C2
- **2ILU** (S2, 1 ch): A15 default 29.076857965920958 → C10; A16 domains 0.0688 () → C21; A20 ligand-tether:NDP:A3482:3 7.916; ligand-tether:NDP:A4482:3 7.916 → C2
- **3L9G** (S2, 1 ch): A3 domains: A2:27-40>A2:43-47, A2:76-90>A2:94-104, A2:178-188>A2:194-199, A2:227-241>A2:245-2… → C31; A4 domains: A1:27-40>A1:43-47, A1:43-47>A1:51-53, A1:51-53>A1:57-70, A1:57-70>A1:76-90 → C31; A6 1393 hits (domains 1393) → C4 1225, C6 146, C5 22; A7 domains 3 (text#147:'D1'; text#148:'D2'; text#149:'D3') → C31; A15 default 24.227503589905304 → C10; A16 domains 120 () → C31; A20 ligand-tether:URC:A2300:0 9.011; ligand-tether:URC:A2300:1 9.367; ligand-tether:URC:A2300:… → C2
- **4XG1** (S2, 4 ch): A6 8 hits (beta 2, default 2, deposited 2, focus 2) → C1 8; A20 ligand-tether:LLP:A508:3 12.42; ligand-tether:LLP:B508:3 12.432 → C2
- **6KTL** (S2, 4 ch): A20 ligand-tether:INS:D402:3 7.436 → C2
- **6YF9** (S2, 180 ch): A6 10319 hits (focus 10319) → C4 8584, C6 1735; A7 focus 8 (label:FQ:100-115; label:GF:100-115; label:CZ:100-115) → C4
- **7GLY** (S2, 2 ch): A20 ligand-tether:RDX:A407:3 7.41; ligand-tether:RDX:B405:3 8.594 → C2
- **7OKA** (S2, 30 ch): A15 default 1.122071319273961; focus 1.122071319273961 → C12; A20 ligand-tether:VFT:B302:2 8.908 → C2
- **7SFH** (S2, 1 ch): A6 9 hits (beta 3, default 3, deposited 3) → C3 9; A20 ligand-tether:91I:A1401:3 6.983; ligand-tether:91I:A2401:3 6.936 → C2
- **9PXS** (S2, 1 ch): A15 default 14.208779293414773 → C10
- **1SLV** (S3, 2 ch): A6 8 hits (beta 2, default 2, deposited 2, domains 1, focus 1) → C1 8; A15 default 6.474416199657219 → C10
- **2QNA** (S3, 2 ch): A6 1 hits (focus 1) → C4 1
- **4UM3** (S3, 40 ch): A2 301.2 s (limit 180 s) → C30; A4 domains: M:3-13>M:27-38, T:150-153>T:171-184, W:27-38>W:47-59, W:86-88>W:96-97 → C30, C8; A6 1701 hits (beta 517, default 354, deposited 354, domains 476) → C6 941, C4 569, C5 191; A16 domains 36 () → C4, C30; A20 ligand-tether:09R:E301:2 8.706; ligand-tether:09R:B301:2 7.379; ligand-tether:09R:C301:1 6… → C2; A2b time: 4UM3 took longer than 90 s to draw, the limit of this app. Draw part of it (Chains &… → C30
- **4XEI** (S3, 7 ch): A6 4 hits (beta 1, default 1, deposited 1, domains 1) → C1 4
- **5DWZ** (S3, 8 ch): A15 default 13.473711234569112; focus 13.473711234569112 → C10; A20 ligand-tether:MPO:C406:0 6.767; ligand-tether:MPO:C406:1 6.922 → C2
- **5SB5** (S3, 6 ch): A6 1 hits (domains 1) → C5 1; A20 ligand-tether:GTP:A501:1 8.362; ligand-tether:GTP:A501:3 11.455; ligand-tether:GTP:C501:3 … → C2
- **5X4Z** (S3, 24 ch): A6 57 hits (beta 11, default 11, deposited 11, domains 13, focus 11) → C3 34, C1 23; A16 domains 0.0261 () → C21
- **5JB1** (S4, 6 ch): A6 145 hits (beta 61, default 28, deposited 28, domains 28) → C4 66, C3 40, C1 28, C6 11
- **6AD0** (S4, 6 ch): A6 15 hits (beta 3, default 3, deposited 3, domains 6) → C3 9, C1 6
- **6DPV** (S4, 12 ch): A6 81 hits (beta 6, default 7, deposited 7, domains 61) → C5 33, C6 27, C4 21; A20 ligand-tether:GTP:A501:3 10.243; ligand-tether:GDP:B501:2 6.621; ligand-tether:GDP:B501:3 … → C2
- **6MHS** (S4, 4 ch): A6 38 hits (beta 8, default 8, deposited 8, domains 14) → C3 22, C6 12, C1 4
- **7AOD** (S4, 24 ch): A6 802 hits (beta 300, default 121, deposited 121, domains 257, focus 3) → C4 357, C6 272, C5 78, C1 49, C3 46; A15 default inf → C22; A16 domains 4 () → C4
- **7ARD** (S4, 51 ch): A3 default: H:96-105>H:107-118, H:107-118>H:121-151, H:121-151>H:163-168; deposited: H:96-105… → C8; A4 domains: D:11-16>D:28-34, D:37-44>D:52-56, G:346-348>G:349-353, G:359-372>G:376-379 → C8; A6 2795 hits (beta 949, default 596, deposited 596, domains 654) → C6 1282, C4 1070, C1 288, C5 131, C3 24; A9 default 1 (chain q); domains 1 (chain q) → C15; A16 domains 14 () → C4; A20 ligand-tether:PC7:N403:2 8.217; ligand-tether:PTY:H301:1 7.076; ligand-tether:PTY:H301:2 8… → C2
- **7JU4** (S4, 61 ch): A2 192.74 s (limit 180 s) → C32; A3 default: J:132-139>J:145-160, J:165-171>J:183-195, J:206-215>J:224-243, J:224-243>J:252-25… → C8; A4 domains: 6:44-49>6:51-53, 6:51-53>6:59-61, D:10-27>D:65-69, H:252-259>H:267-273 → C8; A6 2351 hits (beta 459, default 359, deposited 359, domains 1174) → C4 1365, C6 613, C5 335, C1 23, C3 15; A16 domains 0.0598 (); domains 77 () → C4, C21; A20 ligand-tether:GDP:6501:0 12.513; ligand-tether:GDP:6501:2 22.862; ligand-tether:GDP:6501:3… → C2
- **7QUP** (S4, 65 ch): A3 focus: 9C:312-320>9C:325-336 → C8; A6 5630 hits (focus 5630) → C4 3391, C6 1976, C5 249, C3 14; A7 focus 1 (label:1B:405-424) → C4; A15 focus 44.23032373350822 → C12; A20 ligand-tether:GTP:7A500:3 9.859; ligand-tether:GDP:7B600:3 8.4 → C2
- **8ESZ** (S4, 43 ch): A6 919 hits (beta 264, default 212, deposited 212, domains 231) → C4 308, C6 259, C5 152, C1 149, C3 51; A16 domains 1 () → C4; A20 ligand-tether:3PE:S7302:1 7.1; ligand-tether:WSF:6201:1 6.814; ligand-tether:3PE:4506:1 12… → C2
- **8PQV** (S4, 1 ch): A6 3 hits (beta 1, default 1, deposited 1) → C3 3; A20 ligand-tether:ADP:A4701:3 7.349; ligand-tether:ATP:A4702:3 6.978; ligand-tether:ANP:A4704:… → C2
- **8R5X** (S4, 4 ch): A6 14 hits (beta 4, default 4, deposited 4, domains 2) → C1 14
- **8RRX** (S4, 12 ch): A6 2701 hits (beta 9, default 14, deposited 14, domains 77, focus 2587) → C4 1868, C6 694, C3 103, C5 36; A16 domains 1 () → C23; A20 ligand-tether:CFF:D5104:0 8.658 → C2
- **8RVL** (S4, 34 ch): A4 domains: R:106-119>R:131-138, T:78-99>T:105-118, 7:178-181>7:193-194, 8:44-46>8:67-68 → C8; A6 589 hits (beta 196, default 104, deposited 104, domains 185) → C4 301, C6 177, C1 111; A16 domains 6 () → C4
- **9ON5** (S4, 3 ch): A6 6 hits (beta 1, default 1, deposited 1, domains 3) → C3 6
- **9P1X** (S4, 30 ch): A6 173 hits (beta 3, default 3, deposited 3, domains 5, focus 159) → C5 42, C3 37, C1 35, C4 31, C6 28; A15 focus 12.765786219803081 → C10, C22, C12; A20 ligand-tether:ICS:A602:1 6.928; ligand-tether:HCA:A603:1 6.623; ligand-tether:HCA:A603:2 7… → C2
- **2DGS** (S5, 1 ch): A6 3 hits (beta 1, default 1, deposited 1) → C1 3
- **6B5W** (S5, 1 ch): A9 default 1 (chain A) → C15
- **5JDN** (S6, 1 ch): A20 ligand-tether:MYS:A407:1 7.095; ligand-tether:SR:A403:2 6.659 → C2
- **5LDW** (S6, 45 ch): A6 801 hits (beta 178, default 153, deposited 190, domains 177, membrane 103) → C6 218, C1 211, C4 157, C3 144, C5 71; A16 domains 2 () → C4; A18 TM helix A:84-103 0.71; TM helix L:85-106 0.40; TM helix M:24-41 -0.39; TM helix M:89-111 … → C17; A20 ligand-tether:FMN:F501:1 6.791; ligand-tether:SF4:F502:3 16.485; ligand-tether:NAP:P501:2 … → C2
- **6JLN** (S6, 39 ch): A3 default: E:17-39>E:42-47, E:42-47>E:72-81, b:392-395>b:398-402, b:398-402>b:410-411; depos… → C8; A5 default 2 (['b:398-402', 'b:429-433']; ['b:429-433', 'b:336-347']); deposited 2 (['b:398-4… → C8; A6 2307 hits (beta 462, default 407, deposited 407, domains 621, membrane 410) → C5 1119, C6 618, C4 489, C1 49, C3 32; A15 default inf → C22, C12; A18 TM helix A:110-136 0.24; TM helix A:143-165 -0.03; TM helix B:93-116 0.60; TM helix B:234-… → C17; A20 ligand-tether:BCR:A410:1 8.118; ligand-tether:BCR:A410:2 7.186; ligand-tether:CLA:C505:0 8… → C2
- **6JO6** (S6, 19 ch): A6 3154 hits (beta 646, default 632, deposited 632, domains 612, membrane 632) → C5 2623, C6 208, C4 199, C1 114, C3 10; A18 no band drawn → C17; A20 ligand-tether:CLA:A805:0 7.599; ligand-tether:CLA:A805:1 7.879; ligand-tether:CLA:A842:0 6… → C2
- **6PS3** (S6, 1 ch): A18 TM helix A:267-298 0.73 → C27
- **7PIU** (S6, 6 ch): A9 default 1 (chain P); domains 1 (chain P) → C15; A18 TM helix R:48-70 0.70 → C27
- **7SMR** (S6, 5 ch): A6 34 hits (beta 6, default 6, deposited 6, domains 12, membrane 4) → C1 20, C6 5, C5 5, C4 4; A18 loop in the slab B:117-129>B:141-150; loop in the slab B:171-173>B:190-201; loop in the sl… → C26; A20 ligand-tether:CCE:A601:1 6.814; ligand-tether:CLR:C503:2 11.13 → C2
- **8W9Y** (S6, 6 ch): A6 34 hits (beta 5, default 5, deposited 5, domains 13, membrane 6) → C1 29, C4 5; A18 loop in the slab C:225-249>C:271-280; loop in the slab F:225-249>F:271-280 → C26
- **3L2V** (S7, 4 ch): A6 2 hits (domains 2) → C11 2; A17 angle_median 40.6 (limit 0.05); angle_max 58.3 (limit 0.05); contact_auc 0.481 (limit 0.6)… → C16
- **8BW5** (S7, 3 ch): A17 pair_f1 0.8 (limit 0.95); stem_recall 0.5 (limit 0.85) → C16
- **8QEM** (S7, 26 ch): A3 focus: A:221-239>A:248-255, B:221-239>B:248-255, strand Z:28-29 → C8; A4 focus: C:221-239>C:248-255 → C8; A6 382 hits (beta 1, default 1, deposited 1, domains 1, focus 378) → C6 194, C4 136, C1 46, C14 6; A15 focus 64.76575120350533 → C10, C12; A16 domains 3 () → C23
- **9GEP** (S7, 12 ch): A4 default: D:53-80>D:85-86, H:53-80>H:85-86, N:564-566>N:573-575; deposited: D:53-80>D:85-86… → C8; A6 30 hits (beta 7, default 7, deposited 7, domains 9) → C1 13, C5 7, C3 6, C4 3, C11 1; A15 default 1.168226716003388 → C12; A20 ligand-tether:HEM:N603:3 7.3; ligand-tether:NAG:N602:1 8.861 → C2
- **9K9V** (S7, 5 ch): A6 1 hits (domains 1) → C11 1
- **9YLM** (S7, 3 ch): A17 stems_nearer_frac 0.0 (limit 1.0) → C16
- **6SO9** (S8, 2 ch): A6 16 hits (beta 4, default 4, deposited 4, focus 4) → C11 16
- **7MW8** (S8, 12 ch): A9 default 1 (chain K) → C15
- **7XSP** (S8, 3 ch): A6 20 hits (beta 5, default 5, deposited 5, focus 5) → C3 8, C11 8, C1 4
- **7Y7C** (S8, 55 ch): A3 default: gap after a:1439, gap after Z:27, gap after Z:37, gap after Z:61; deposited: gap … → C24; A4 default: c:208-213>c:222-224; deposited: c:208-213>c:222-224; beta: c:208-213>c:222-224; f… → C8; A6 7254 hits (beta 1698, default 1363, deposited 1363, domains 1467, focus 1363) → C14 3429, C4 1747, C11 932, C6 824, C1 222, C5 59, C3 41; A16 domains 26 () → C4; A17 fallbacks 5 (limit 0); na_crossings 16 (limit 0); bead_overlaps 31 (limit 0); gap_label_co… → C4, C11, C14, C24
- **8VWQ** (S8, 56 ch): A3 default: strand 1:662-663, strand 1:934-935, strand 1:943-945, strand 1:1506-1507; deposit… → C24; A4 default: M:110-121>M:128-132, O:3-21>O:25-30, P:37-41>P:44-51, P:58-65>P:68-75; deposited:… → C8; A6 7600 hits (beta 1753, default 1380, deposited 1380, domains 1707, focus 1380) → C14 3549, C4 1927, C6 826, C11 758, C1 309, C3 124, C5 107; A16 domains 11 () → C4; A17 stems_nearer_frac 0.993 (limit 1.0); fallbacks 5 (limit 0); na_crossings 9 (limit 0); bead… → C4, C11, C14, C24
- **9C75** (S8, 3 ch): A3 default: strand R:28-29; deposited: strand R:28-29; beta: strand R:28-29; domains: strand … → C24; A6 6 hits (beta 1, default 1, deposited 1, domains 1, focus 2) → C14 6; A17 fallbacks 1 (limit 0) → C24
- **9CJI** (S8, 4 ch): A6 4 hits (beta 1, default 1, deposited 1, focus 1) → C11 4; A17 bead_overlaps 2 (limit 0) → C11
- **3OWW** (S9, 2 ch): A6 6 hits (default 3, focus 3) → C14 6; A17 spearman 0.484 (limit 0.6); na_crossings 2 (limit 0) → C14, C16
- **8SA2** (S9, 2 ch): A6 28 hits (default 14, focus 14) → C14 28
- **1WWL** (S10, 2 ch): A6 6 hits (beta 1, default 1, deposited 1, domains 2, focus 1) → C1 5, C28 1; A15 default 12.261689266655294 → C10
- **7BBH** (S10, 3 ch): A6 34 hits (beta 1, default 1, deposited 1, domains 2, focus 29) → C1 17, C4 11, C3 6; A16 domains 2 () → C23
- **7BRG** (S10, 3 ch): A6 17 hits (beta 5, default 5, deposited 5, domains 1, focus 1) → C3 13, C1 4; A7 focus 1 (loop-length:A:237-241>A:267-276) → C3
- **8K8E** (S10, 5 ch): A6 4 hits (beta 1, default 1, deposited 1, domains 1) → C28 4; A20 ligand-tether:CLR:C301:1 8.543 → C2
- **8OZJ** (S10, 11 ch): A4 domains: D:747-750>D:755-765, D:769-777>D:790-795, H:747-750>H:755-765, H:769-777>H:790-79… → C8; A6 25 hits (beta 5, default 3, deposited 3, domains 14) → C4 13, C1 11, C6 1; A16 domains 6 () → C4; A20 ligand-tether:PTY:H1101:1 6.62; ligand-tether:PTY:D1101:1 7.642 → C2
- **1RM5** (S11, 3 ch): A19 only 14 of 30 and 28 elements paired; A1:63-66 unpaired; A1:79-81 unpaired; A1:84-87 unpai… → C18
- **2HKI** (S11, 1 ch): A6 7 hits (beta 1, default 1, deposited 1, domains 4) → C3 7; A19 only 14 of 28 and 30 elements paired; A1:66-69 unpaired; A1:87-90 unpaired; A1:94-97 unpai… → C18
- **6FY9** (S11, 1 ch): A6 6 hits (beta 2, default 2, deposited 2) → C1 6
- **6FYF** (S11, 4 ch): A6 4 hits (beta 1, default 1, deposited 1, focus 1) → C1 4
- **6YPX** (S11, 4 ch): A15 default 4.567999924904073; focus 4.567999924904073 → C10
- **6YQD** (S11, 2 ch): A15 default 4.4242737012349895 → C10
- **9YKJ** (S11, 2 ch): A4 beta: A:194-198>A:204-210, A:664-674>A:683-696 → C8; A6 3 hits (beta 1, default 1, deposited 1) → C3 3
- **AF_AFA0A4E9ES35F1** (S12, 1 ch): A6 9 hits (beta 3, default 3, deposited 3) → C1 6, C3 3

