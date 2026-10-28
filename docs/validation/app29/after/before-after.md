# Plan 29: before and after, both halves

Written by scripts/app29_after.py report at git b84ebc3 (clean). The check half is opened here, once (ruling P29-5): it was measured, never tuned on; nothing was changed because of it.

Runs (scripts/app29.py run --ui, the app's code path; the page session on the tune half only):

- before: `docs/validation/app29/before`: the committed before record; git 549bef9 (check rows 549bef9), Linux-7.0.0-38-generic-x86_64-with-glibc2.43
- rerun: `docs/validation/app29/after/before-rerun`: the before commit rerun on the after run's host, same settings; git f538981 (check rows f538981), Linux-6.8.0-138-generic-x86_64-with-glibc2.35
- after: `docs/validation/app29/after`: the after record; git 162e067 (check rows 162e067), Linux-6.8.0-138-generic-x86_64-with-glibc2.35

Pass: drawn, judged, no criterion failed (A1-A20) and above the readability floor (R1). Crash: no figure (an error, a job killed). Floor exit and criteria failed: of the drawn. Page failed: the tune half's page session raised, stopped or a download gave no whole file. Counts carry a Wilson 95% interval of their share; seconds (the draw job's wall time, what a visitor waits; draw s: draw_job's own time) are of the drawn entries, each with a distribution-free 95% interval (the order statistics at n x the Wilson interval of the quantile's share). Seconds compare only within one host: zyla-office's cores are about twice as fast as gpu-desktop's, and the runs on gpu-desktop drew 96 entries at once.

## The table

### tune half

| run | entries | pass | crash | refused | floor exit | criteria failed | page failed | seconds median | seconds p95 | draw s median | draw s p95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| before (committed, zyla-office) | 47 | 45/47 [86-99%] | 0/47 [0-8%] | 0/47 [0-8%] | 1/47 [0-11%] | 1/47 [0-11%] | 0/47 [0-8%] | 0.82 [0.36-1.19] | 17.131 [2.7-23.66] | 0.79 [0.34-1.13] | 13.694 [2.63-21.29] |
| before (rerun, gpu-desktop) | 47 | 45/47 [86-99%] | 0/47 [0-8%] | 0/47 [0-8%] | 1/47 [0-11%] | 1/47 [0-11%] | 0/47 [0-8%] | 2.86 [1.5-4.12] | 45.742 [11.58-66.17] | 2.78 [1.44-3.96] | 38.973 [9.05-61.66] |
| after (gpu-desktop) | 47 | 47/47 [92-100%] | 0/47 [0-8%] | 0/47 [0-8%] | 0/47 [0-8%] | 0/47 [0-8%] | 0/47 [0-8%] | 2.8 [1.65-4.49] | 43.943 [7.16-66.39] | 2.41 [1.52-3.98] | 40.072 [5.44-61.77] |

### check half

| run | entries | pass | crash | refused | floor exit | criteria failed | page failed | seconds median | seconds p95 | draw s median | draw s p95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| before (committed, zyla-office) | 32 | 27/32 [68-93%] | 0/32 [0-11%] | 1/32 [1-16%] | 1/32 [1-16%] | 3/32 [3-24%] | 0/32 [0-11%] | 0.7 [0.28-1.17] | 14.315 [2.68-42.19] | 0.62 [0.25-1.14] | 12.24 [2.61-39.22] |
| before (rerun, gpu-desktop) | 32 | 27/32 [68-93%] | 0/32 [0-11%] | 2/32 [2-20%] | 1/32 [1-16%] | 2/32 [2-20%] | 0/32 [0-11%] | 2.615 [0.83-4.38] | 21.347 [9.88-56.42] | 2.335 [0.7-3.8] | 18.39 [6.16-52.05] |
| after (gpu-desktop) | 32 | 28/32 [72-95%] | 0/32 [0-11%] | 2/32 [2-20%] | 0/32 [0-11%] | 2/32 [2-20%] | 0/32 [0-11%] | 2.415 [0.86-5.34] | 22.116 [9.81-65.26] | 2.19 [0.69-4.49] | 19.065 [6.04-61.79] |

### both halves

| run | entries | pass | crash | refused | floor exit | criteria failed | page failed | seconds median | seconds p95 | draw s median | draw s p95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| before (committed, zyla-office) | 79 | 72/79 [83-96%] | 0/79 [0-5%] | 1/79 [0-7%] | 2/79 [1-9%] | 4/79 [2-12%] | 0/79 [0-5%] | 0.785 [0.47-1.01] | 17.747 [7.52-23.66] | 0.735 [0.42-0.96] | 15.141 [5.61-21.29] |
| before (rerun, gpu-desktop) | 79 | 72/79 [83-96%] | 0/79 [0-5%] | 2/79 [1-9%] | 2/79 [1-9%] | 3/79 [1-11%] | 0/79 [0-5%] | 2.86 [1.69-3.61] | 43.762 [20.71-60.16] | 2.54 [1.58-3.34] | 33.378 [16.19-55.69] |
| after (gpu-desktop) | 79 | 75/79 [88-98%] | 0/79 [0-5%] | 2/79 [1-9%] | 0/79 [0-5%] | 2/79 [1-9%] | 0/79 [0-5%] | 2.68 [1.73-3.88] | 41.298 [14.25-65.26] | 2.36 [1.6-3.48] | 38.692 [13.02-61.77] |

## What changed (rerun -> after, same host)

| id | half | outcome | moved | seconds | page failed | cause |
|---|---|---|---|---|---|---|
| 26EL | tune | pass -> pass | app.svg, app.png | 6.92 -> 4.49 | False -> False | Task 8b Fix A (fix-copies): focus auto draws A alone, not A, F and M on a diagonal; page 1091 x 1414 -> 525 x 668 pt |
| 29TB | tune | pass -> pass | app.svg, app.png | 1.4 -> 1.25 | False -> False | Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's) |
| 31FD | tune | fail floor -> pass | app.svg, app.png | 34.22 -> 40.24 | False -> False | Task 5 (fix-floor): the two-copy complex was below the floor whole, so one subunit is drawn; page 3456 x 619 -> 1477 x 883 pt; floor exit -> pass |
| 9SVL | check | pass -> pass | app.svg, app.png | 1.02 -> 1.1 | False -> False | Task 4 (fix-ligand): 11 ligand groups moved and the page widened 361.14 -> 374.38 pt with them, so the rest shifts |
| 9SXF | tune | pass -> pass | app.svg, app.png | 0.77 -> 0.86 | False -> False | Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's) |
| 9WNE | tune | pass -> pass | app.svg, app.png | 6.67 -> 6.47 | False -> False | Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (4 SVG groups, every one a ligand's) |
| 9X4L | tune | fail A20 -> pass | app.svg, app.png | 46.93 -> 45.53 | False -> False | Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's); its A20 (ligand tethers local) now passes |
| 9X4M | check | fail floor -> pass | app.svg, app.png | 56.42 -> 65.26 | False -> False | Task 5 (fix-floor): the two-copy complex was below the floor whole, so one subunit is drawn; page 3456 x 1203 -> 1746 x 1466 pt; floor exit -> pass |
| 9YN3 | check | pass -> pass | app.svg, app.png | 20.72 -> 20.48 | False -> False | Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's) |

## The committed before record against its rerun (another host)

| id | half | outcome | moved |
|---|---|---|---|
| 9SY6 | check | fail A3,A6 -> refused | app.svg, app.png |

## Over the hosted job limit on gpu-desktop: the same entries alone on zyla-office

The hosted draw job is killed at 90 s (refused: time). An entry near it on zyla-office passes it on gpu-desktop's slower cores, so each such entry (and each lead that timed out) was also run at both commits on zyla-office, a few at once (before git f538981, after git 279289f). Not in the table above.

| id | half | committed before | before | after | seconds | moved |
|---|---|---|---|---|---|---|
| 9FZL | lead | - | fail floor | pass | 55.56 -> 74.14 | app.svg, app.png |
| 9SY6 | check | fail A3,A6 | fail A3,A6 | fail A3,A6 | 45.31 -> 47.04 | - |
| 9SYG | check | refused | refused | refused | 90.12 -> 90.06 | - |

## The parked measures (P29-13: the third sample, the NA controls, the science measure)

### the third sample (random100-third, default view): 64b03ad (main before plan 29's fixes) against the final commit

- before git 64b03ad, after git 162e067; entries 100 -> 100; pass 82 -> 82; keys: moved
- worse: none
- better: 9FZL: fallback counts_hidden 2 -> 1; 9FZL: fallback helix_name_moved 22 -> 0; 9FZL: fallback page_text_small 1 -> 0
- moved (no judge orders it): 7SPI: area 426.01 -> 338.03; 9FZL: area 65458.72 -> 14561.27
- the SVG moved on: 1NID, 3ABR, 3G6E, 3WAU, 4GWY, 5S4U, 6IZP, 6PO4, 6VBG, 7F5S, 7SPI, 7XBK, 8EM7, 9FZL
- **1NID**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (9 SVG groups, every one a ligand's); no key moved
- **3ABR**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's); no key moved
- **3G6E**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (15 SVG groups, every one a ligand's); no key moved
- **3WAU**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's); no key moved
- **4GWY**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (10 SVG groups, every one a ligand's); no key moved
- **5S4U**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (6 SVG groups, every one a ligand's); no key moved
- **6IZP**: Task 4 (fix-ligand): the two B2R marks moved inward, so the page is 206.74 -> 203.96 pt wide and every element shifts with it (10 ligand groups moved, the rest translated); no key moved
- **6PO4**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (12 SVG groups, every one a ligand's); no key moved
- **6VBG**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (12 SVG groups, every one a ligand's); no key moved
- **7F5S**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (3 SVG groups, every one a ligand's); no key moved
- **7SPI**: Task 8b Fix A (fix-copies: copies modelled in part are one protein to focus auto): TraV's A1 drawn once, B1 lent where it completes A1's sheet; area 426.01 -> 338.03; pass, criteria and fallbacks identical
- **7XBK**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (21 SVG groups, every one a ligand's); no key moved
- **8EM7**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (7 SVG groups, every one a ligand's); no key moved
- **9FZL**: Task 5 (fix-floor: focus auto draws one subunit of a two-copy complex whose whole figure is below the readability floor): one subunit, area 65458.72 -> 14561.27, the fallbacks counts_hidden 2 -> 1, helix_name_moved 22 -> 0 and page_text_small 1 -> 0 (better)

### the third sample: its committed report (7696947, views all; read-only) against 64b03ad (code before plan 29's fixes)

- before git (the committed report), after git 64b03ad; entries 100 -> 100; pass 83 -> 82; keys: worse
- worse: 2YEJ: PASS -> FAIL; 2YEJ: A19 now fails; 3FT8: PASS -> FAIL; 3FT8: A19 now fails; 3G6E: fallback glyph_width_guessed 0 -> 52; 3P2W: PASS -> FAIL; 3P2W: A19 now fails; 4AS1: fallback glyph_width_guessed 0 -> 4; 4B2T: fallback counts_hidden 0 -> 1; 5U3R: PASS -> FAIL; 5U3R: A19 now fails; 5XMX: PASS -> FAIL; 5XMX: A19 now fails; 6LTO: PASS -> FAIL; 6LTO: A19 now fails; 6M67: PASS -> FAIL; 6M67: A19 now fails; 6QIK: fallback counts_hidden 0 -> 1; 6QIK: fallback glyph_width_guessed 0 -> 52; 6QIK: fallback na_run_unclear 0 -> 1; 6ZTQ: fallback counts_hidden 0 -> 1; 7F5S: fallback counts_hidden 0 -> 5; 7F5S: fallback glyph_width_guessed 0 -> 68; 8CRC: PASS -> FAIL; 8CRC: A19 now fails; 8XJV: fallback counts_hidden 0 -> 5; 8XJV: fallback glyph_width_guessed 0 -> 76; 8XJV: fallback helix_name_moved 0 -> 19; 8XJV: fallback na_gap_unclear 0 -> 3; 8XJV: fallback na_tail_straight 0 -> 4; 8XJV: fallback na_trace_unclear 0 -> 3; 8XJV: fallback page_text_small 0 -> 1; 9FZL: fallback counts_hidden 0 -> 2; 9FZL: fallback helix_name_moved 0 -> 22; 9FZL: fallback page_text_small 0 -> 1; 9NMN: fallback loop_over_text 0 -> 1
- better: 4AS1: FAIL -> PASS; 4AS1: A17 now passes; 4YJ5: FAIL -> PASS; 4YJ5: A15 now passes; 6DQN: FAIL -> PASS; 6DQN: A18 now passes; 6VTX: FAIL -> PASS; 6VTX: A17 now passes; 7F5S: FAIL -> PASS; 7F5S: A3 now passes; 7F5S: A6 now passes; 8PIN: FAIL -> PASS; 8PIN: A15 now passes; 9FZL: FAIL -> PASS; 9FZL: A16 now passes; 9FZL: A3 now passes; 9FZL: A4 now passes; 9FZL: A6 now passes
- moved (no judge orders it): 1RKA: area 1121.25 -> 679.88; 1VSK: area 352.68 -> 361.36; 3EGZ: area 1478.07 -> 1564.02; 3FT8: area 508.73 -> 433.34; 4GWY: area 1942.44 -> 1633.59; 6PO4: area 524.5 -> 643.66; 6QIK: area 27610.17 -> 27350.66; 7Y63: area 1783.14 -> 1764.17; 8EM7: area 8203.86 -> 8239.72; 8U1N: area 1723.55 -> 1117.64; 9DU3: area 552.51 -> 654.76
- the SVG moved on: 1RKA, 1VSK, 3EGZ, 3FT8, 3G6E, 4AS1, 4GWY, 6PO4, 6QIK, 7F5S, 7Y63, 8EM7, 8U1N, 8XJV, 9DU3
- for reference, not judged: 30 entries differ (older code or another view)

### the third sample: its committed report (7696947, views all; read-only) against the final commit

- before git (the committed report), after git 162e067; entries 100 -> 100; pass 83 -> 82; keys: worse
- worse: 2YEJ: PASS -> FAIL; 2YEJ: A19 now fails; 3FT8: PASS -> FAIL; 3FT8: A19 now fails; 3G6E: fallback glyph_width_guessed 0 -> 52; 3P2W: PASS -> FAIL; 3P2W: A19 now fails; 4AS1: fallback glyph_width_guessed 0 -> 4; 4B2T: fallback counts_hidden 0 -> 1; 5U3R: PASS -> FAIL; 5U3R: A19 now fails; 5XMX: PASS -> FAIL; 5XMX: A19 now fails; 6LTO: PASS -> FAIL; 6LTO: A19 now fails; 6M67: PASS -> FAIL; 6M67: A19 now fails; 6QIK: fallback counts_hidden 0 -> 1; 6QIK: fallback glyph_width_guessed 0 -> 52; 6QIK: fallback na_run_unclear 0 -> 1; 6ZTQ: fallback counts_hidden 0 -> 1; 7F5S: fallback counts_hidden 0 -> 5; 7F5S: fallback glyph_width_guessed 0 -> 68; 8CRC: PASS -> FAIL; 8CRC: A19 now fails; 8XJV: fallback counts_hidden 0 -> 5; 8XJV: fallback glyph_width_guessed 0 -> 76; 8XJV: fallback helix_name_moved 0 -> 19; 8XJV: fallback na_gap_unclear 0 -> 3; 8XJV: fallback na_tail_straight 0 -> 4; 8XJV: fallback na_trace_unclear 0 -> 3; 8XJV: fallback page_text_small 0 -> 1; 9FZL: fallback counts_hidden 0 -> 1; 9NMN: fallback loop_over_text 0 -> 1
- better: 4AS1: FAIL -> PASS; 4AS1: A17 now passes; 4YJ5: FAIL -> PASS; 4YJ5: A15 now passes; 6DQN: FAIL -> PASS; 6DQN: A18 now passes; 6VTX: FAIL -> PASS; 6VTX: A17 now passes; 7F5S: FAIL -> PASS; 7F5S: A3 now passes; 7F5S: A6 now passes; 8PIN: FAIL -> PASS; 8PIN: A15 now passes; 9FZL: FAIL -> PASS; 9FZL: A16 now passes; 9FZL: A3 now passes; 9FZL: A4 now passes; 9FZL: A6 now passes
- moved (no judge orders it): 1RKA: area 1121.25 -> 679.88; 1VSK: area 352.68 -> 361.36; 3EGZ: area 1478.07 -> 1564.02; 3FT8: area 508.73 -> 433.34; 4GWY: area 1942.44 -> 1633.59; 6PO4: area 524.5 -> 643.66; 6QIK: area 27610.17 -> 27350.66; 7SPI: area 426.01 -> 338.03; 7Y63: area 1783.14 -> 1764.17; 8EM7: area 8203.86 -> 8239.72; 8U1N: area 1723.55 -> 1117.64; 9DU3: area 552.51 -> 654.76; 9FZL: area 65458.72 -> 14561.27
- the SVG moved on: 1NID, 1RKA, 1VSK, 3ABR, 3EGZ, 3FT8, 3G6E, 3WAU, 4AS1, 4GWY, 5S4U, 6IZP, 6PO4, 6QIK, 6VBG, 7F5S, 7SPI, 7XBK, 7Y63, 8EM7, 8U1N, 8XJV, 9DU3, 9FZL
- for reference, not judged: 38 entries differ (older code or another view)

### the NA set (na_fixes_record's cases and controls and the NA campaign's entries, 69, default view): 64b03ad against the final commit

- before git 64b03ad, after git 162e067; entries 69 -> 69; pass 54 -> 54; keys: identical
- worse: none
- better: none
- moved (no judge orders it): none
- the SVG moved on: 1KX5, 4YBB, 6IZP
- A17 measures (69 entries): changed on none; worse: none; better: none
- **1KX5**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (3 SVG groups, every one a ligand's); no key moved
- **4YBB**: Task 4 (fix-ligand: render._place_mark keeps a mark clear of loop lines; the tether-reach step): only its ligand marks, tethers and labels moved (28 SVG groups, every one a ligand's); no key moved
- **6IZP**: Task 4 (fix-ligand): the two B2R marks moved inward, so the page is 206.74 -> 203.96 pt wide and every element shifts with it (10 ligand groups moved, the rest translated); no key moved

### the NA campaign (scripts/na_validate.py, docs/validation/na)

- before git 64b03ad, after git 162e067; entries 39 -> 39; pass 34 -> 34
- measures changed on: none; worse: none; better: none

### The science measure (`scripts/app29_light.py science --cache`, figcheck at 3d6c798 against now)

- 630 structures judged of 630 (cache 491, tests/data 92, tune 47), 0 not drawn
- science findings changed (Task 7's barrel_finding against the old count) on: 9X4L; new findings: none
  - 9X4L S3: ['sheet 7: 6 barrel strands drawn, the pairing cycle has 5'] -> none
- barrel blocks holding a ghost lent by another copy: 1, in 1 structures: 9X4L
- a ring strand lent across the cut (science.lent_ring takes a lent ghost into the ring): none found
- the probe failed on: none
