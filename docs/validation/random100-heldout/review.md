# Random-100 held-out sample: review

**This sample is no longer held out.** Its first run (80 of 100, below) was the input of plan 15
(`docs/superpowers/plans/2026-10-17-15-heldout-failures.md`), whose Tasks 1-8 fixed the mechanisms behind its
failures and checked each fix on its entries. Its pass rate now measures fit, as the tuned sample's does. The count
for structures foldmap was not tuned on is the third sample's, drawn after those fixes and run once:
[docs/validation/random100-third/](../random100-third/review.md).

## The plan 15 run: 95 of 100

Run 2026-10-18 at 7696947 on zyla-lab (32 cores), every view, `--jobs 10`, 862 s wall:

    python scripts/random100.py run --views all --jobs 10 \
        --sample docs/validation/random100-heldout/sample.json --out DIR

`report.md`, `report.json.gz` and `run.log` are that run, judged by the same criteria and limits as before (no limit in
`criteria.py` or `random100.py` was raised in plan 15). A1, A2, A2b, A4, A5, A7-A14 and A18-A20 pass on every entry that
has something to judge; A3, A6, A15, A16 and A17 fail on the five entries below. The run itself recorded its commit as
"unknown" (zyla-lab's synced copy of the tree has no .git); 7696947 was stamped into the reports afterwards from that
synced tree, and the reports were judged again from the measurements, with the same counts.

### The 20 failures of the first run

| id | first run fails | now | what changed it (plan 15) |
|---|---|---|---|
| 2FR2 | A3 | pass | Task 1: an accessory whose stack gap stays shut stands on its partner's other side (`_column_heights`) |
| 2NAD | A3, A5, A6 | pass | Task 1: a copy's element its mate lacks is slid off the shared ones after the take (`_slid_off`) |
| 8BCU | A3 | pass | Task 1: closed pockets between accessory stacks (`way_out` follows runs to the sheet's edge, `unpocketed`); Task 2: the domains-view A4 that Task 1's layout showed (a domain leader on a loop's way out) |
| 12BM | A3, A6, A16 | pass | Task 1 (closed pockets); the domain-name collision went with the rerouted loop |
| 8D3D | A1 (focus killed after 900 s), A5-A8, A15 | pass | Task 7: the disulfide second starts share one work budget per layout; focus view 63.7 s in this run (10 entries at once; 47.6 s alone in Task 7) |
| 7BKC | A3, A4, A6, A16 (domains) | pass | Task 2: domain leaders keep off loops' ways out, names keep room |
| 6MJ2 | A18 no band | pass | Task 6: a band from helices that only reach a face (`cli._face_medians`) |
| 7JI3 | A18 loop in the slab | pass | Task 6: loops on one side stay out of the slab |
| 7VVK | A18 no band | pass | Task 6: a slab no element spans moves to one they do, judged by the chains that cross it (polar fraction 0.149 against 0.15) |
| 6ZKT | A4, A6 (domains) | pass | Task 2: a domain name's leader kept off a loop's way out |
| 6L9Z | A6 (domains) | pass | Task 3: a terminus label grazed by a joint moved off the run |
| 5K7Z | A3, A15, A17 | pass | Task 1: one conformer letter per model (`io._one_conformer`) makes E/F a duplex (pair_f1 0.182 → 1.0); Task 5: copies on a duplex frame are tied (`cli._side_by_side`) |
| 2FD8 | A6 | pass | Task 3: letters judged against element outlines, end labels moved when crossed; Task 4: end-label room by outline |
| 13FZ | A3, A6, A15, A17 | pass | Task 3: nucleotides nudged off elements that do not hold them, letters off loops; Task 5: tied orbits slid clear as one (`_tied_again`) |
| 2BBJ | A15 | pass | Task 5: of two copy items tied to one mate in one protomer, the smaller is let go (`_one_per_mate`) |
| 9GTP | A15 (focus 52.11) | **fail** (focus 41.39) | Task 5, counted: M~Q 25.77 inherent (Q's 46-group strands share P's fused strands; refit 0 without them); c~a, e~a, f~a 41.39 a layout cause, parked (`_host` picks a shared sheet's host protomer by residue count and flips between copies) |
| 9K4F | A17 stems_nearer_frac 0.875 | **fail** (0.875) | Task 4, counted: layout displacement, deferred (needs `_attract` work with its own bench) |
| 6S13 | A3, A6, A17 | **fail**: A3, A6 (domains), A17 pair_f1 0.944 | Task 4: beads_missing 134 was measured wrong (fixed in `verify_contacts`), runs side by side and gaps across a free end fixed, the tie label off stem 166; counted: domains-view split hairpin stem overlap (A3, and its 3 A6 first-pair numbers), pair_f1 0.944 inherent (3.58 Å model; the depositors' own H-bond records agree with foldmap on 86 % of the disputed pairs) |
| 8ZTU | A3, A6, A16, A17 | **fail**: A3 (domains) | Tasks 3-4 fixed its A6, A16 and A17; counted: the domains-view split hairpin stem overlap (strands a:2701-2711) |
| 7NZ2 | A3, A4, A6, A16 (domains), A17 | **fail**: A3, A6, A16 (domains), A17 spearman 0.446 | counted: spearman 0.446 (Task 4: anisotropic stretch of a 44-chain complex, deferred). **Regression**, see below: its domains view passed A3/A4/A6/A16 from Task 2 to Task 6 and fails A3, A6 and A16 again since Task 7 |

Fifteen of the twenty pass; the five that fail are counted failures (rulings in the plan 15 ledger), and no entry
that passed in the first run fails now.

### 7NZ2's domains view: a regression from Task 7

At 7696947 the domains view of 7NZ2 has three loops with no route (B3:211-221>B3:227-408, B3:1051-1123>B3:1136-1140,
F3:182-186>F3:197-205), and the plain curves cross three domain names (D7 · A1, D7 · B3, D1 · A3: A6 and A16). Tasks
2-6 drew this view with none. Drawn alone at each commit (`random100.py view domains`, zyla-lab), the change comes
with 1bf5f36 (Task 7):

| commit | layout area | fallback curves | second searches refused | route nodes |
|---|---|---|---|---|
| b2f8ebb (Task 5), a963bb6 (Task 6) | 149,316 | 0 | 6 | 7,849,567 |
| 1bf5f36 (Task 7), 7696947 | 149,316 | 3 | 18 | 7,849,567 |

The page is the same size (3456 × 3143 pt, below the 200 in cap Task 7 added), and the refused second searches
tripled, which points at Task 7's shared page search budget (window routers now draw on the page's budget instead of
a fresh one each), exercised in full only by a full run. Not isolated further and not fixed here: 7NZ2 already
counts as a failure (A17), so the count does not change; the domains-view failure is new since Task 6.

**After the run (2026-10-18, plan 15's final fix round).** The cause was Task 7's window budget: a coarse page's fine
windows (`route._fine_route`) drew on the page's own second-search budget, which a huge page has mostly spent by the
time its windows run, so every window search was refused. The windows now share a pool of their own
(`_MAX_NODES_WINDOWS`, 8,000,000 nodes per page). Drawn again (zyla-lab), 7NZ2's domains view has its Task 6 drawing
back: signature 4987d5a72ecd, 0 fallback curves, 6 second searches refused, no text collision; 7NZ2 still fails
A17 (spearman 0.446). The counts above are as measured at 7696947 and are not re-quoted.

### Changes against the first run (field by field, `report.json`)

- **A10** (DSSP agreement, after Task 1's one-conformer read): identical detail on all 100 entries.
- **A15 "does not apply" → "not judged"** on 15 entries (Task 5: copies that exist but are drawn untied are not
  judged, with the reason): 5S5Z, 6XER, 3L4G, 6ET5, 8BCU, 6V6S, 7ZDP, 8WMV, 6ZKT, 8WH0, 6L9Z, 2VBO, 7NZ2, 9YER, 9QK1.
- **Domains view** (A4 loops over labels, A6 text collisions, A16 name collisions; leaders counted in the SVG):
  7BKC A4 2 → 0, A6 70 → 0, A16 6 → 0; 7NZ2 A4 5 → 0, A6 9 → 3, A16 6 → 3, leaders 36 → 31; 8ZTU A6 7 → 0,
  A16 1 → 0; 12BM A6 1 → 0, A16 1 → 0, leaders 11 → 8; 6ZKT A4 1 → 0, A6 1 → 0; 2NAD A6 2 → 0; 6L9Z A6 1 → 0.
  Leaders only: 2BBJ 4 → 10, 7KH0 4 → 3, 8BCU 6 → 5, 7ZDP 3 → 2, 8WMV 4 → 3, 13FZ 2 → 1, 8FLU 1 → 0, 7VVK 1 → 0.
  Domain-name counts are unchanged on every entry.
- **A20 tethers** (default view): lengths moved on 9 entries (1PTY, 2NAD, 12BM, 5K7Z, 2FD8, 7NZ2, 13FZ, 8ZTU, 3U3H),
  as their layouts changed; A20's verdict is the same on all of them (pass on 5K7Z, 7NZ2, 8ZTU; not judged on the
  others). The longest tether of 7NZ2 (11.575) and of 8ZTU (36.762) is unchanged.
- **Nucleic-acid metrics** (A17's view): 5K7Z pair_f1 0.182 → 1.0, stem_recall 0 → 1, fallbacks 1 → 0; 13FZ fallbacks
  1 → 0, na_protein_crossings 9 → 7; 6S13 beads_missing 134 → 0, na_crossings 2 → 0, bead_overlaps 7 → 0,
  na_protein_crossings 46 → 48; 8ZTU beads_missing 1 → 0, bead_overlaps 2 → 0, gap_label_collisions 1 → 0. No other
  entry's NA metrics moved, including the NA entries where Task 5 lays out twice (let-go copy items).
- **sheets_differ** (ties left out of A15 because the copies' sheets differ): 24 over all views, as before, entry by
  entry.
- **Route fallbacks**: 2FR2 14 → 0, 2NAD 16 → 0, 8BCU 32 → 0, 12BM 48 → 0, 7BKC 28 → 0, 7NZ2 4 → 3 (summed over
  views). No other entry has one.
- **Time** (view processes, 10 entries at once): 8D3D focus 900 s (killed) → 63.7 s; 9GTP focus 29.2 → 69.5 s
  (Task 5 lays it out twice; 5.6 M of its 8 M route budget); 7BKC domains 85.8 → 67.6 s with 1,256,079 nodes, as in
  Task 7; the slowest view is 7NZ2 domains, 192.8 s. A2 (the default view, limit 180 s for XL) passes on every
  entry; the slowest default is 8ZTU, 60.7 s.

S5 (automatic domains against CATH, reported, not gated and not in the count): 164 chains scored in 37 entries.
One-domain chains: 98 pass, 20 fail (the automatic domains split them); chains with two or more CATH domains: 18
pass, 28 fail.

## The first run: 80 of 100

The layout's constants were tuned on the first random-100 sample (`docs/validation/random100/`), so this sample was
drawn with a new seed and none of its ids and run once at 5650dbd, every view, no fix; after the F1 review the eight
entries whose result that check fix could change were drawn again at 44a1c1b. 80 of 100 passed, against 92 of 100 on
the tuned sample; the twenty failures are the table above. That run's `report.md` and `report.json.gz` are in git
history (as committed at 1531604).

## The two-state pairs

- **P3, 8U9T/7IT6, is no two-state pair.** Both entries carry one UniProt accession, P0DTD1, the SARS-CoV-2 pp1ab
  polyprotein, but 8U9T's chain is nsp5 (Mpro, 306 residues) and 7IT6's is nsp1 (113 residues): 50% identity over
  the 113 aligned residues. A19 and the pair's type read "does not apply: not one protein" (criteria.one_protein).
  Plan 15 Task 8 made the sampler pair on one protein by sequence, not on the accession alone; this sample is not
  redrawn.
- **2BBJ/2BBH (P5):** equal names always match; A19 passes on all pairs.
- The pairs' types as their files read now are in `report.md` ("Two-state pairs").
