# Plan 29 Task 7: a lent ghost is no barrel strand (triage proposal 4, ruling P29-7)

Written by scripts/app29_light.py report at e8f1e47 (clean); ruling P29-13's light loop: the tune half and the bench, the check half, the third sample and the NA controls left to Task N+1.

**Verdict: kept.** Expected to change: 9X4L's S3 only (foldmap check; the app does not show S3); every figure byte-identical.

## Mechanism

figcheck.science_checks counted every ghost lent to a barrel block (layout._unwrap: a neighbouring copy's strand drawn across an unrolled ring's cut leaves its ghost where it pairs) as a ring strand. 9X4L's sheet 7 (chain 3's 5-strand barrel) holds the ghost of 7:7-8, lent beside the accessory 3:208-213 and paired with it alone: ring 6 against science.pairing_cycle's 5. science.barrel_finding now counts a lent ghost in the ring only when it is joined to a ring strand in the pairing graph (pairing_cycle's own edge rule, science._joined), directly or through another such ghost (science.lent_ring); the others are counted as accessories are. A check fix: no figure may change.

## The tune half (scripts/app29.py run, the app's way)

- before: `docs/validation/app29/fix-floor/run` (git 975f336); after: `docs/validation/app29/fix-ghost/run` (git e8f1e47)
- entries 47 -> 47; missing: none
- the app's SVG changed on: none; its PNG on: none
- criteria, floor, error or refusal changed on: none
- draw seconds, all tune entries: 118.76 -> 109.74 (8 entries drawn at once: indicative)
- 9X4L: draw 14.61 -> 14.42 s; the entry's process (draw, PNG, checks) 16.85 -> 16.59 s

## bench --check --exact

- 22 of 22 identical, exit 0

## The fix's own measure (science.json)

- 139 entries judged of 139 (the tune half and tests/data, drawn as foldmap check draws), 0 not drawn
- science findings changed on: 9X4L; new findings: none
- 9X4L S3: ['sheet 7: 6 barrel strands drawn, the pairing cycle has 5'] -> none
- S3 findings left after: none

### The control

No barrel block in tests/data or the tune half holds a ring strand lent across a cut: the measles F sheets (8UTF) and the other blocks with lent ghosts (8UUP, 1SLV, 4V4F-ring, 14CU) are drawn open, judged by science.open_finding, which this fix leaves as it was. The guard is a synthetic ring (tests/test_science.py::test_a_ring_strand_lent_across_the_cut_still_counts); the barrels whose S3 passes (1EMA, 1TIM, 1BXW-A: tests/test_science.py::test_a_control_barrel_figure_passes_and_a_moved_strand_fails) and 8UTF's S3 are pinned.
