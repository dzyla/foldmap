# Plan 29 triage: the tune half's failures and faults, by mechanism

Plan 29 Task 3, written by hand at f538981 (read-only: no code changed). Inputs:

- the before record's tune half, read through `scripts/app29.py` (load_tune), and nothing of the check half;
  `docs/validation/app29/before/report.json`: 47 drawn, 0 refused, 0 errors, UI 47/47 clean, `A20` failed on
  `9X4L`, `R1` failed on `31FD`;
- a contact sheet of the 47 tune figures (the app's PNG download), every one looked at, the doubtful ones at full
  size;
- `foldmap check --deep --json` of the 47 tune entries and of the plan 28 leads, drawn on zyla-office, and the
  figure notes (the SVG's "How this figure was made");
- three read-only probes on zyla-office (the figure drawn the app's way, then measured): ligand marks crossed by a
  loop line, disulfide bars over other elements' bodies, empty page share, the focus decision of the floor
  entries, and a profile of each entry over the hosted slow limit (`app_jobs.SLOW_SECONDS` (4.0 s)).

Visibility: 1 a crash or refusal, 2 a failed check or page criterion, 3 a look-only fault. Score: entries × weight,
the weight 3, 2 or 1 for visibility 1, 2 or 3 (so a crash outranks a look); the entries are the tune half's and the
leads' together, each counted in its own column.

## The faults, one row each

| id | what the user sees | check that sees it | mechanism |
|---|---|---|---|
| `31FD` | a 48 x 9 in page, text 5.0 pt; the floor's note offers only "chain A in full" | `R1` | `cli._focus_chains`: 'auto' never takes one subunit of two (`cli._subunits` 2), so this C2 photosystem dimer (30 chains, 233 elements) is drawn whole; with focus protomer it is above the floor. The app has no protomer control: `floor._phrase` drops that option for the app |
| `6JLN`, `9FZL` (leads) | the same: 48 x 9 and 48 x 10 in, 4.5 and 2.0 pt | `R1` | the same; both pass the floor with focus protomer (C2, 2 sectors) |
| `6ET5`, `7JU4`, `7NZ2`, `8WMV` (leads) | text 2.4-5.5 pt on a page held at 48 in | `R1` | no symmetry found (no sector), so no subunit to draw; `layout.page_fit` holds the longer side at `layout.MAX_PAGE_IN` (48.0) and these projections are 8:1 to 1.3:1 |
| `9X4L` | ADP 32001's tether to chain 5 is 7.89 units long | `A20` | `render._place_mark`: its kept step looks for a spot only under the keep-outs the first spot met, or stricter (out to `render._SPOT_FAR` x `render._SPOT_RADIUS`) for a spot that holds chain 5's target within `render._TETHER_REACH` (6.5), finds none on this crowded page, and the spot stays (2.42 from its main contact) |
| `9X4L` | nothing on the page | `S3` (foldmap check only; the app does not show it) | `figcheck.science_checks` counts every ghost lent to a barrel block as a ring strand. Block 7 (chain 3's 5-strand barrel) holds the ghost of 7:7-8, which `layout._unwrap` lent beside the accessory 3:208-213; it pairs with no ring strand. Ring 6, `science.pairing_cycle` 5: a false finding |
| `9X4L` | chain 7's first strand drawn alone at the page's right edge beside a dashed copy, its loop crossing the page | none | `layout._unwrap` (by design: a strand lent across an unrolled ring's cut is drawn one past the row's end) |
| `9SXF`, `9T57`, `29TB`, `31FD`, `35VC`, `9WNE`, `9YOA` | a ligand or ion box sitting on a loop line: the line stops at the box, as if the loop ran into the ligand | none (a line under an opaque mark is hidden, and `A4` judges element labels) | `render._place_mark`: when no spot within `render._SPOT_RADIUS` is clear of the loop lines and the tethers drawn before, the next keep-out layers leave the loop lines out; none of these marks is crowded (`render._least_bad_spot` not used) |
| `9Z0S`, `9XIC`, `9TMF`, `9YOA`, `9ZSJ`, `12KF`, `34QD`, `23FR` | a disulfide bar drawn straight across strands or helices it does not join (`9Z0S`: 3.3 units over four strands, both Ig domains of both chains) | none (layout.md: the length of a disulfide bar is seen by no check) | `render.disulfide_points` gives the two cysteines' points; the bar is drawn straight between them over whatever lies between |
| `23FR` | 17 s for one subunit (11 x 10 units): automatic redrawing turns off | none (time is said, not judged; `A2b` passes under 90 s) | `symmetry._helical`: every pair of copies is tried as the screw, and each candidate is tested on every copy against every member of its cluster, about n^4/4 pair fits (1.26 million `symmetry._Assembly.pair` calls for 42 copies; 11.6 of 12.4 s) |
| `9YPK`, `9YPN` | 21 and 18 s | none | `route.route_loops` 33 and 27 s of 40 and 34 s under the profiler; mechanism below it not found |
| `31FD`, `9X4L`, `9XIC` | 12, 14 and 9 s | none | `layout._separate` (5.6 and 8.6 s) and `layout._attract` (9.6 s for `9XIC`); mechanism below them not found |
| `7F5S`, `7JU4`, `7NZ2`, `9FZL` (leads) | 40-56 s, half of the hosted 90 s | none | `route.route_loops` 28-57 s; mechanism below it not found |
| `9NMN` (lead) | the loop into helix α153 runs under its name | `A4`, `A6` | `route._Router.route` takes `loop_over_text` (the figure's notes): the name sits in the only gap the loop can pass, between α153 and the next helix. The gap is narrower than a name and a loop lane; why the layout leaves it so is not traced. Plan 28 found that `layout._STACK_NUDGE` off removes it, with no mechanism joining the two |
| `9FZL` (lead, domains view) | domain names D1, D2 crossed by loops | `A4`, `A6`, `A16`, `A3` | the page's route budget is spent (14 wider searches not run): `route._Router.route` takes `loop_over_text`. `5NF8` passes; `7F5S` fails on nucleic-acid tie labels, not names. The `layout:_NAME_OPEN` prune is not the cause (the prune record: keys identical) |
| `4OO8` (lead, the Cas9 supplementary figure) | strand letters collide in small sheets; the P1 stem's base letters overlap | `A6` | only at font_scale 3 or more: 0 text collisions at 1 and 2, 8 at 3, 157 at 6.25 (101 base letter pairs, 22 strand letter pairs). Text grows with font_scale; the strand pitch and the stem's rung spacing do not. The app's Text size stops at 2.0: the app never draws it. A manuscript figure fault, not an app one |
| `26EL` | 66% of the page empty: three chains on a diagonal | none | not found (the probe's empty share; the layout's projected places were not traced further) |
| `22MZ` | a pseudoknotted riboswitch whose backbone crosses itself four times, ligand boxes on its lines | none (`D3` passes: 25 of 25 pairs) | not found |
| `44LS` | the two copies of a homodimer drawn with different sheets | none (`A15` passes) | not found (copies' DSSP not compared) |
| `9WXY` | a page with one sentence, "no helices or strands found" | none | by design: a 16-residue peptide with no helix or strand; nothing is drawn |
| every tune id | a ResourceWarning in the page session's stderr | none | not found: the app's own temporary folders (`app_jobs.figure_bytes`, the layout upload) are closed by their with blocks; the warning comes from the session harness or Streamlit's AppTest |

Leads closed without a fault on the app's path: the latent `layout._open_at_ends` / `like.Like.order` bug (no trigger in
the tune half; the freeze rule waits for an unseen case), and the `layout:_NAME_OPEN` prune above.

## The groups, ranked

| rank | mechanism | tune entries | leads | visibility | score | proposed |
|---|---|---|---|---|---|---|
| 1 | slow layout or routing (`route.route_loops`, `layout._separate`, `layout._attract`) | 5 | 4 | 3 | 9 | no: no mechanism pinned below these |
| 2 | disulfide bar drawn straight over bodies (`render.disulfide_points`) | 8 | 0 | 3 | 8 | 6 |
| 3 | the floor of large asymmetric complexes (`layout.page_fit`) | 0 | 4 | 2 | 8 | no: no fix short of a new arrangement |
| 4 | ligand mark on a loop line (`render._place_mark`, the layers without loop lines) | 7 | 0 | 3 | 7 | 1 |
| 5 | a two-copy complex drawn whole under the floor (`cli._focus_chains`; no protomer control in the app) | 1 | 2 | 2 | 6 | 2 |
| 6 | a lent ghost counted as a barrel strand (`figcheck.science_checks`) | 1 | 0 | 2 | 2 | 4 |
| 7 | a bound tether left long on a crowded page (`render._place_mark`, its kept step) | 1 | 0 | 2 | 2 | 5 |
| 8 | the loop under α153 (`route._Router.route`; the gap not traced) | 0 | 1 | 2 | 2 | no: half traced |
| 9 | the screw search's n^4 fits (`symmetry._helical`) | 1 | 0 | 3 | 1 | 3 |

## Proposals (the controller rules)

Ordered by what a fix buys on the tune half and how sure its mechanism is: the disulfide bars score highest of the
proposed but come last, as their fault status is in doubt (see 6).

1. **Ligand marks off the loop lines** (`render._place_mark`). Failing test: `9SXF` (295 residues), the A1JRG box
   crossed by a loop. Should change: `9SXF`, `9T57`, `29TB`, `31FD`, `35VC`, `9WNE`, `9YOA`. Risk: marks move
   further out and their tethers grow on every entry with a crowded pocket; `A20` bounds tethers to other chains, so
   it must hold on both halves.
2. **The two-copy complex under the floor** (`cli._focus_chains`, `floor._phrase`). Failing test: `31FD` (the floor
   fails; one subunit passes). Should change: `31FD`; leads `6JLN`, `9FZL`. Either the app offers one subunit
   (`floor._phrase` then names it: no figure changes), or 'auto' takes one of two when the whole falls under the
   floor (a figure changes for every C2 complex under the floor, a second drawing's time). Risk: a C2 complex above
   the floor must keep its whole figure.
3. **The screw search** (`symmetry._helical`). Failing test: `23FR`, its time (the answer must be the same).
   Should change: `23FR`'s time only; every figure byte-identical on both halves, the third sample and the bench.
   Risk: a candidate tested only against its nearest copies (`symmetry._Assembly.nearest`) can pick another step
   than the first copy in cluster order; any figure that moves rejects the fix. A filament of twice the copies would
   take about 16 times as long, past the hosted 90 s: a refusal.
4. **A lent ghost is no barrel strand** (`figcheck.science_checks`). Failing test: `9X4L`, `S3` (a check fix: no
   figure changes). Should change: `9X4L`'s `S3` only. Risk: a strand lent across a ring's cut that is a ring strand
   (the docstring's measles F sheet) must still count.
5. **A bound tether kept within reach** (`render._place_mark`, its kept step). Failing test: `9X4L`, `A20`. Should change:
   `9X4L`. Risk: looser keep-outs put the mark on loop lines or text, against proposal 1; the two meet in one function
   and should be ruled together.
6. **Disulfide bars around bodies** (`render.disulfide_points`). Failing test: `9Z0S`. Should change: `9Z0S`,
   `9XIC`, `9TMF`, `9YOA`, `9ZSJ`, `12KF`, `34QD`, `23FR`. Risk: high. Hand-drawn Ig figures draw the B-F bridge
   across the flattened sandwich, so this may be no fault. A routed bar is a new drawing rule (the freeze rule), and
   no check measures it yet.

Not proposed: the slow layout and routing with no mechanism found below `route.route_loops`, `layout._separate` and
`layout._attract`; the floor of asymmetric complexes; `9NMN`; the Cas9 figure (font_scale 6.25 is the manuscript's,
past the app's 2.0); `26EL`, `22MZ`, `44LS` (not found); `9X4L`'s lent strand and `9WXY` (by design).
