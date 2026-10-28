# Plan 29 fix record: fix-floor

Written by: foldmap 0.1.0, git 9869e37568f4, code fingerprint 5ef865ec11fd82a1, scripts/app29_fix.py sha256 f554f34e1941; Python 3.12.9 on Linux-6.18.33.1-microsoft-standard-WSL2-x86_64-with-glibc2.39; gemmi 0.7.5, numpy 2.5.3, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

**Verdict: better** (identical or better on every measure keeps the fix).

- worse, tune: none
- worse, check: none
- worse, third: none
- worse, na: none
- worse, bench: none
- worse, floor: none
- better, tune: 31FD: now above the floor
- better, check: floor_failed
- better, third: 9FZL: fallback counts_hidden 2 -> 1, 9FZL: fallback helix_name_moved 22 -> 0, 9FZL: fallback page_text_small 1 -> 0
- better, na: none
- better, floor: 31FD: now above the floor, 6JLN: now above the floor, 9FZL: now above the floor

## Tune half (per entry; the check half in its totals only)

| id | criteria failed | floor ok | svg changed | seconds |
|---|---|---|---|---|
| `9SXY` | - -> - | True -> True | no | 0.27 -> 0.26 |
| `9SXU` | - -> - | True -> True | no | 0.28 -> 0.27 |
| `9SXF` | - -> - | True -> True | no | 0.25 -> 0.24 |
| `9SXG` | - -> - | True -> True | no | 0.26 -> 0.25 |
| `29LB` | - -> - | True -> True | no | 0.24 -> 0.24 |
| `29TB` | - -> - | True -> True | no | 0.36 -> 0.34 |
| `9T57` | - -> - | True -> True | no | 0.36 -> 0.34 |
| `9ZFX` | - -> - | True -> True | no | 0.57 -> 0.59 |
| `14CU` | - -> - | True -> True | no | 1.24 -> 1.22 |
| `30CY` | - -> - | True -> True | no | 0.93 -> 0.92 |
| `44LS` | - -> - | True -> True | no | 1.17 -> 1.21 |
| `9YL9` | - -> - | True -> True | no | 0.63 -> 0.62 |
| `20XT` | - -> - | True -> True | no | 0.74 -> 0.73 |
| `9WVR` | - -> - | True -> True | no | 0.29 -> 0.27 |
| `9Z0S` | - -> - | True -> True | no | 1.25 -> 1.28 |
| `9TCF` | - -> - | True -> True | no | 0.91 -> 0.89 |
| `9YOA` | - -> - | True -> True | no | 3.61 -> 3.63 |
| `38NR` | - -> - | True -> True | no | 0.28 -> 0.27 |
| `23FR` | - -> - | True -> True | no | 17.24 -> 16.99 |
| `31FD` | - -> - | False -> True | yes | 12.97 -> 15.09 |
| `22VF` | - -> - | True -> True | no | 1.52 -> 1.49 |
| `12XK` | - -> - | True -> True | no | 1.18 -> 1.13 |
| `35VC` | - -> - | True -> True | no | 1.16 -> 1.11 |
| `23CX` | - -> - | True -> True | no | 1.28 -> 1.27 |
| `28HY` | - -> - | True -> True | no | 0.46 -> 0.45 |
| `26EL` | - -> - | True -> True | no | 1.57 -> 1.54 |
| `9X4L` | - -> - | True -> True | no | 16.46 -> 16.46 |
| `9NYA` | - -> - | True -> True | no | 0.14 -> 0.13 |
| `45JF` | - -> - | True -> True | no | 0.1 -> 0.1 |
| `34XH` | - -> - | True -> True | no | 0.13 -> 0.13 |
| `9WXY` | - -> - | True -> True | no | 0.05 -> 0.05 |
| `9ZSJ` | - -> - | True -> True | no | 2.66 -> 2.64 |
| `9YPK` | - -> - | True -> True | no | 23.51 -> 23.33 |
| `9YPN` | - -> - | True -> True | no | 20.49 -> 20.42 |
| `9WNE` | - -> - | True -> True | no | 1.44 -> 1.45 |
| `22MZ` | - -> - | True -> True | no | 0.42 -> 0.43 |
| `9XIC` | - -> - | True -> True | no | 9.53 -> 9.51 |
| `9TMF` | - -> - | True -> True | no | 3.41 -> 3.33 |
| `34QD` | - -> - | True -> True | no | 0.78 -> 0.79 |
| `12KF` | - -> - | True -> True | no | 1.48 -> 1.46 |
| `30CU` | - -> - | True -> True | no | 0.96 -> 0.95 |
| `30CV` | - -> - | True -> True | no | 1.0 -> 1.0 |
| `9YLA` | - -> - | True -> True | no | 0.68 -> 0.68 |
| `9YLG` | - -> - | True -> True | no | 0.66 -> 0.65 |
| `32IL` | - -> - | True -> True | no | 0.17 -> 0.17 |
| `9SXD` | - -> - | True -> True | no | 0.27 -> 0.26 |
| `9YLD` | - -> - | True -> True | no | 0.26 -> 0.26 |

## Check half (totals)

| total | before | after |
|---|---|---|
| entries | 32 | 32 |
| drawn | 31 | 31 |
| refused | 1 | 1 |
| error | 0 | 0 |
| not_measured | 0 | 0 |
| png_error | 0 | 0 |
| criteria_failed | 3 | 3 |
| floor_failed | 1 | 0 |
| ui_run | 0 | 0 |
| ui_failed | 0 | 0 |

## Third sample, default view (keys: pass, failing criteria, named fallbacks, area)

- 100 entries; pass 82 -> 82; verdict moved
- worse: none
- better: 9FZL: fallback counts_hidden 2 -> 1; 9FZL: fallback helix_name_moved 22 -> 0; 9FZL: fallback page_text_small 1 -> 0
- moved: 9FZL: area 65458.72 -> 14561.27
- SVG moved (no key reads it): 9FZL
- the committed report (views all, an older code) against the after run, default view: for reference: pass 83 -> 82; 34 worse, 18 better, 11 moved

## NA cases and controls (plan 28's NA measure)

- 46 entries; verdict identical
- worse: none
- better: none
- moved: none
- SVG moved (no key reads it): none

## bench --check --exact

- 22 entries; exit 0 -> 0
- changed: none

## The readability floor and focus auto (plan 29 Task 5)

- Tune half under the floor: 1 -> 0.
- Tune entries focus auto lays out as an assembly of two, drawn whole unless that figure is under the floor: 12. Under the floor whole, so drawn again as one subunit: 31FD; every other one is drawn whole as before (figure unchanged: 12XK, 20XT, 30CU, 30CV, 30CY, 44LS, 9X4L, 9YL9, 9YLA, 9YLG, 9ZFX; changed: none).

Leads and controls, each drawn the app's way (random100.app_view, the hosted limits; seconds: the job's, both drawings when it draws twice) and here (draw_spec, focus auto: seconds, text pt, chains drawn in full, of all; 0 is the whole). Before -> after.

| entry | role | app floor | app seconds | auto floor | auto seconds | text pt | in full |
|---|---|---|---|---|---|---|---|
| `31FD` | lead | under -> above | 13.17 -> 14.62 | under -> above | 10.71 -> 13.75 | 5.026 -> 6.0 | 0 -> 15 of 30 |
| `6JLN` | lead | under -> above | 19.4 -> 21.4 | under -> above | 17.89 -> 22.41 | 4.528 -> 6.0 | 0 -> 19 of 39 |
| `9FZL` | lead | under -> above | 48.93 -> 62.67 | under -> above | 45.39 -> 58.88 | 2.049 -> 6.0 | 0 -> 49 of 98 |
| `3UTB-protein` | control | - -> - | - -> - | above -> above | 0.33 -> 0.32 | 6.0 -> 6.0 | 0 -> 0 of 8 |
| `3UTB` | control | above -> above | 1.16 -> 1.02 | above -> above | 0.88 -> 0.81 | 6.0 -> 6.0 | 0 -> 0 of 8 |

The two drawings apart, after (focus none: the whole; focus protomer: the one subunit), seconds and text pt:

- `31FD`: whole 11.73 s, 5.026 pt (under); one subunit 2.62 s, 6.0 pt (above)
- `6JLN`: whole 17.82 s, 4.528 pt (under); one subunit 4.38 s, 6.0 pt (above)
- `9FZL`: whole 47.33 s, 2.049 pt (under); one subunit 13.44 s, 6.0 pt (above)
- `3UTB-protein`: whole 0.32 s, 6.0 pt (above); one subunit 0.21 s, 6.0 pt (above)
- `3UTB`: whole 0.89 s, 6.0 pt (above); one subunit 0.69 s, 6.0 pt (above)

Predicting the floor from the layout alone (layout.page_scale of its extent plus render's pad, before routing and drawing): the drawn floor's answer on 20 of 20 drawings measured; the largest gap between the predicted and the drawn text, 0.015 pt. The route budget, the floor's other reason, is known only once the loops are routed. spec.draw_spec draws the whole, judges it and draws again instead.

## Pictures (gitignored, for looking)

- `31FD`: before docs/superpowers/work/app29/fix-floor/png/31FD-before.png, after docs/superpowers/work/app29/fix-floor/png/31FD-after.png
