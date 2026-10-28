# The app on a fresh random sample

Plan 29. Written by: foldmap 0.1.0, git f538981b8cfd (asserted by runner, dirty False), code fingerprint f624dcc2403b2888, scripts/app29.py sha256 cdf60fa1e052; Python 3.13.2 on Linux-6.8.0-138-generic-x86_64-with-glibc2.35; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Every entry drawn the app's way (scripts/app29.py run: random100.app_view, the size check and app_jobs.draw_job under the hosted limits, the app's default settings, the PNG download) and the figure judged as `foldmap check` judges it: A1-A20 and R1, the readability floor.

Sample: `docs/validation/app29/sample.json` (sha256 b3ed519636c1), seed 20261027, released after 2026-10-13 and on or before 2026-10-27; split: tune 47, check 32.

## The sample per stratum

| stratum | title | target | drawn | short |
|---|---|---|---|---|
| S1 | X-ray, 1 protein chain | 12 | 12 | 0 |
| S2 | X-ray homo-oligomer (symmetric) | 12 | 10 | 2 |
| S3 | X-ray hetero-complex | 10 | 6 | 4 |
| S4 | cryo-EM | 15 | 15 | 0 |
| S5 | NMR, multi-model | 8 | 7 | 1 |
| S6 | Membrane | 8 | 0 | 8 |
| S7 | Protein-DNA | 7 | 6 | 1 |
| S8 | Protein-RNA | 7 | 6 | 1 |
| S9 | RNA only | 3 | 3 | 0 |
| S10 | Glycoproteins | 5 | 4 | 1 |
| S11 | One protein, two states (5 pairs) | 10 | 10 | 0 |
| S12 | AlphaFold DB models | 3 | 0 | 3 (not drawn (--no-alphafold): an AlphaFold DB pool has no release window (ruling Q61)) |

## Totals

| half | entries | drawn | refused | error | not measured | PNG error | criteria failed | floor failed | page run | page failed |
|---|---|---|---|---|---|---|---|---|---|---|
| tune | 47 | 47 | 0 | 0 | 0 | 0 | 1 | 1 | 47 | 0 |
| check | 32 | 30 | 2 | 0 | 0 | 0 | 2 | 1 | 0 | 0 |
| all | 79 | 77 | 2 | 0 | 0 | 0 | 3 | 2 | 47 | 0 |

Peak memory of the draw job: not measured unless app_jobs hands it back; `process` in report.json is the whole entry process's (its wall time and max RSS include the PNG, the redraw and the checks).

## Entries

The check half is reported in the totals only: its ids are never shown beside an outcome here, nor on the console, so no fix is tuned on them.

| id | stratum | half | residues | outcome | criteria failed | floor | seconds | draw s |
|---|---|---|---|---|---|---|---|---|
| 9SXY | S1 | tune | 339 | drawn | - | pass | 1.02 | 0.94 |
| 9SXU | S1 | tune | 342 | drawn | - | pass | 0.85 | 0.78 |
| 9SXF | S1 | tune | 295 | drawn | - | pass | 0.77 | 0.7 |
| 9SXG | S1 | tune | 292 | drawn | - | pass | 0.95 | 0.88 |
| 29LB | S1 | tune | 265 | drawn | - | pass | 0.8 | 0.72 |
| 29TB | S1 | tune | 257 | drawn | - | pass | 1.4 | 1.33 |
| 9T57 | S1 | tune | 258 | drawn | - | pass | 1.2 | 1.12 |
| 9ZFX | S2 | tune | 1335 | drawn | - | pass | 1.92 | 1.8 |
| 14CU | S2 | tune | 3404 | drawn | - | pass | 6.95 | 4.44 |
| 30CY | S2 | tune | 704 | drawn | - | pass | 3.13 | 3.03 |
| 44LS | S2 | tune | 1048 | drawn | - | pass | 5.46 | 4.14 |
| 9YL9 | S2 | tune | 343 | drawn | - | pass | 2.12 | 2.03 |
| 20XT | S2 | tune | 806 | drawn | - | pass | 2.39 | 2.29 |
| 9WVR | S2 | tune | 170 | drawn | - | pass | 1.02 | 0.94 |
| 9Z0S | S3 | tune | 2682 | drawn | - | pass | 5.38 | 4.57 |
| 9TCF | S3 | tune | 1828 | drawn | - | pass | 2.96 | 2.81 |
| 9YOA | S3 | tune | 1048 | drawn | - | pass | 13.64 | 12.5 |
| 38NR | S3 | tune | 987 | drawn | - | pass | 0.89 | 0.79 |
| 23FR | S4 | tune | 5754 | drawn | - | pass | 42.97 | 23.87 |
| 31FD | S4 | tune | 4174 | drawn | - | fail | 34.22 | 31.14 |
| 22VF | S4 | tune | 4194 | drawn | - | pass | 7.25 | 4.98 |
| 12XK | S4 | tune | 1470 | drawn | - | pass | 4.12 | 3.96 |
| 35VC | S4 | tune | 1560 | drawn | - | pass | 3.81 | 3.68 |
| 23CX | S4 | tune | 1085 | drawn | - | pass | 5.7 | 4.37 |
| 28HY | S4 | tune | 972 | drawn | - | pass | 1.69 | 1.58 |
| 26EL | S4 | tune | 12238 | drawn | - | pass | 6.92 | 4.37 |
| 9X4L | S4 | tune | 8692 | drawn | A20 | pass | 46.93 | 42.33 |
| 9NYA | S5 | tune | 86 | drawn | - | pass | 0.63 | 0.55 |
| 45JF | S5 | tune | 32 | drawn | - | pass | 0.34 | 0.26 |
| 34XH | S5 | tune | 23 | drawn | - | pass | 0.57 | 0.49 |
| 9WXY | S5 | tune | 16 | drawn | - | pass | 0.16 | 0.09 |
| 9ZSJ | S7 | tune | 1729 | drawn | - | pass | 11.58 | 9.05 |
| 9YPK | S8 | tune | 6687 | drawn | - | pass | 66.17 | 61.66 |
| 9YPN | S8 | tune | 6481 | drawn | - | pass | 60.16 | 55.69 |
| 9WNE | S8 | tune | 1310 | drawn | - | pass | 6.67 | 4.97 |
| 22MZ | S9 | tune | 76 | drawn | - | pass | 1.5 | 1.44 |
| 9XIC | S10 | tune | 5496 | drawn | - | pass | 27.98 | 26.09 |
| 9TMF | S10 | tune | 2793 | drawn | - | pass | 13.46 | 12.06 |
| 34QD | S10 | tune | 502 | drawn | - | pass | 2.86 | 2.78 |
| 12KF | S10 | tune | 1476 | drawn | - | pass | 6.6 | 5.04 |
| 30CU | S11 | tune | 1401 | drawn | - | pass | 3.24 | 3.12 |
| 30CV | S11 | tune | 1403 | drawn | - | pass | 3.47 | 3.34 |
| 9YLA | S11 | tune | 343 | drawn | - | pass | 2.41 | 2.33 |
| 9YLG | S11 | tune | 341 | drawn | - | pass | 2.61 | 2.53 |
| 32IL | S11 | tune | 730 | drawn | - | pass | 0.55 | 0.45 |
| 9SXD | S11 | tune | 627 | drawn | - | pass | 0.88 | 0.78 |
| 9YLD | S11 | tune | 341 | drawn | - | pass | 0.87 | 0.8 |

## The page (run --ui)

Each tune-half id typed into the hosted app's page, Draw pressed, the theme changed once and every download button pressed (AppTest; the check half is never driven through the page). Page failed: it raised, the session stopped, or a download gave no whole file.

| id | stage | steps | exception | errors shown | downloads whole | failed |
|---|---|---|---|---|---|---|
| 9SXY | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9SXU | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9SXF | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9SXG | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 29LB | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 29TB | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9T57 | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9ZFX | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 14CU | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 30CY | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 44LS | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YL9 | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 20XT | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9WVR | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9Z0S | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9TCF | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YOA | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 38NR | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 23FR | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 31FD | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 22VF | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 12XK | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 35VC | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 23CX | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 28HY | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 26EL | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9X4L | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9NYA | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 45JF | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 34XH | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9WXY | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9ZSJ | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YPK | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YPN | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9WNE | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 22MZ | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9XIC | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9TMF | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 34QD | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 12KF | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 30CU | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 30CV | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YLA | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YLG | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 32IL | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9SXD | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
| 9YLD | edit | load, theme minimal, downloads | - | 0 | 6/6 | - |
