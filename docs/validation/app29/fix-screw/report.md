# Plan 29 Task 6: the screw search without n^4 pair fits (triage proposal 3, ruling P29-7)

Written by scripts/app29_light.py report at 4ccb38f (clean); ruling P29-13's light loop: the tune half and the bench, the check half, the third sample and the NA controls left to Task N+1.

**Verdict: kept.** Expected to change: 23FR's time only; every figure byte-identical.

## Mechanism

symmetry._helical tried every pair of copies as the screw and tested each on every copy against every member of its cluster: n^4/2 pair fits (23FR's 42 copies: 1.26 million _Assembly.pair calls). A fit's RMSD is at least the distance between the moved and the target centroids of its matched atoms, so a member whose centroid lands farther than the tolerance (plus 1e-6 Å for rounding) is passed over unfitted; the first member within the tolerance, in cluster order, is the one the full search took.

## The tune half (scripts/app29.py run, the app's way)

- before: `docs/validation/app29/fix-floor/run` (git 975f336); after: `docs/validation/app29/fix-screw/run` (git 9fd9c0e)
- entries 47 -> 47; missing: none
- the app's SVG changed on: none; its PNG on: none
- criteria, floor, error or refusal changed on: none
- draw seconds, all tune entries: 118.76 -> 109.75 (8 entries drawn at once: indicative)
- 23FR: draw 8.82 -> 0.69 s; the entry's process (draw, PNG, checks) 17.2 -> 0.98 s

## bench --check --exact

- 22 of 22 identical, exit 0

## The fix's own measure (helical.json)

- 664 structures measured of 664 (tests/data, the tune half, the download cache less the check half), 0 skipped: the return of _helical before and after differs on 0 (MIN_HELIX_COPIES and MIN_HELIX_COPIES_ASKED; every float compared bit for bit)
- a screw found on 20 (4 copies or more) and 20 (3 or more)
- _Assembly.pair calls, all structures: 7,690,410 before, 47,790 after; _helical's seconds: 51.38 before, 3.6 after (run in parallel: indicative)
- more pair calls after on: 6TIK (a screw passes, so every copy's images are gathered at once)

| entry | where | largest cluster | screw (copies) | pair calls | _helical s |
|---|---|---|---|---|---|
| 23FR | tune | 42 | 7 | 1,261,426 -> 2,625 | 8.517 -> 0.243 |
| 26EL | tune | 20 | 10 | 29,302 -> 590 | 0.563 -> 0.105 |
| 32IL | tune | 10 | 5 | 1,903 -> 145 | 0.023 -> 0.009 |
| 6Y7S | tests/data | 6 | 6 | 470 -> 51 | 0.005 -> 0.004 |
| 7A4M | tests/data | 24 | - | 276 -> 276 | 0.041 -> 0.041 |
| 3CRX | tests/data | 4 | - | 67 -> 22 | 0.005 -> 0.006 |
