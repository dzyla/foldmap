# The app on a fresh random sample

Plan 29. Written by: foldmap 0.1.0, git 6ebc5c847704 (asserted by runner, dirty False), code fingerprint 43e1d11adf24d91b, scripts/app29.py sha256 cdf60fa1e052; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

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
| tune | 47 | 47 | 0 | 0 | 0 | 0 | 0 | 1 | 47 | 0 |
| check | 32 | 31 | 1 | 0 | 0 | 0 | 3 | 1 | 0 | 0 |
| all | 79 | 78 | 1 | 0 | 0 | 0 | 3 | 2 | 47 | 0 |

Peak memory of the draw job: not measured unless app_jobs hands it back; `process` in report.json is the whole entry process's (its wall time and max RSS include the PNG, the redraw and the checks).

## Entries

The check half is reported in the totals only: its ids are never shown beside an outcome here, nor on the console, so no fix is tuned on them.

| id | stratum | half | residues | outcome | criteria failed | floor | seconds | draw s |
|---|---|---|---|---|---|---|---|---|
| 9SXY | S1 | tune | 339 | drawn | - | pass | 0.25 | 0.22 |
| 9SXU | S1 | tune | 342 | drawn | - | pass | 0.25 | 0.23 |
| 9SXF | S1 | tune | 295 | drawn | - | pass | 0.24 | 0.21 |
| 9SXG | S1 | tune | 292 | drawn | - | pass | 0.25 | 0.22 |
| 29LB | S1 | tune | 265 | drawn | - | pass | 0.23 | 0.2 |
| 29TB | S1 | tune | 257 | drawn | - | pass | 0.35 | 0.32 |
| 9T57 | S1 | tune | 258 | drawn | - | pass | 0.34 | 0.31 |
| 9ZFX | S2 | tune | 1335 | drawn | - | pass | 0.61 | 0.57 |
| 14CU | S2 | tune | 3404 | drawn | - | pass | 1.28 | 1.15 |
| 30CY | S2 | tune | 704 | drawn | - | pass | 0.96 | 0.92 |
| 44LS | S2 | tune | 1048 | drawn | - | pass | 1.24 | 1.2 |
| 9YL9 | S2 | tune | 343 | drawn | - | pass | 0.65 | 0.62 |
| 20XT | S2 | tune | 806 | drawn | - | pass | 0.73 | 0.7 |
| 9WVR | S2 | tune | 170 | drawn | - | pass | 0.29 | 0.25 |
| 9Z0S | S3 | tune | 2682 | drawn | - | pass | 1.31 | 1.27 |
| 9TCF | S3 | tune | 1828 | drawn | - | pass | 0.94 | 0.89 |
| 9YOA | S3 | tune | 1048 | drawn | - | pass | 3.63 | 3.58 |
| 38NR | S3 | tune | 987 | drawn | - | pass | 0.27 | 0.24 |
| 23FR | S4 | tune | 5754 | drawn | - | pass | 17.24 | 8.83 |
| 31FD | S4 | tune | 4174 | drawn | - | fail | 13.05 | 11.55 |
| 22VF | S4 | tune | 4194 | drawn | - | pass | 1.51 | 1.4 |
| 12XK | S4 | tune | 1470 | drawn | - | pass | 1.2 | 1.14 |
| 35VC | S4 | tune | 1560 | drawn | - | pass | 1.14 | 1.09 |
| 23CX | S4 | tune | 1085 | drawn | - | pass | 1.28 | 1.24 |
| 28HY | S4 | tune | 972 | drawn | - | pass | 0.47 | 0.42 |
| 26EL | S4 | tune | 12238 | drawn | - | pass | 1.63 | 1.43 |
| 9X4L | S4 | tune | 8692 | drawn | - | pass | 16.52 | 14.31 |
| 9NYA | S5 | tune | 86 | drawn | - | pass | 0.14 | 0.11 |
| 45JF | S5 | tune | 32 | drawn | - | pass | 0.11 | 0.08 |
| 34XH | S5 | tune | 23 | drawn | - | pass | 0.13 | 0.11 |
| 9WXY | S5 | tune | 16 | drawn | - | pass | 0.05 | 0.03 |
| 9ZSJ | S7 | tune | 1729 | drawn | - | pass | 2.64 | 2.58 |
| 9YPK | S8 | tune | 6687 | drawn | - | pass | 23.77 | 21.39 |
| 9YPN | S8 | tune | 6481 | drawn | - | pass | 20.87 | 18.57 |
| 9WNE | S8 | tune | 1310 | drawn | - | pass | 1.48 | 1.44 |
| 22MZ | S9 | tune | 76 | drawn | - | pass | 0.43 | 0.41 |
| 9XIC | S10 | tune | 5496 | drawn | - | pass | 9.65 | 8.72 |
| 9TMF | S10 | tune | 2793 | drawn | - | pass | 3.4 | 3.34 |
| 34QD | S10 | tune | 502 | drawn | - | pass | 0.78 | 0.75 |
| 12KF | S10 | tune | 1476 | drawn | - | pass | 1.48 | 1.44 |
| 30CU | S11 | tune | 1401 | drawn | - | pass | 0.94 | 0.9 |
| 30CV | S11 | tune | 1403 | drawn | - | pass | 1.03 | 0.98 |
| 9YLA | S11 | tune | 343 | drawn | - | pass | 0.69 | 0.66 |
| 9YLG | S11 | tune | 341 | drawn | - | pass | 0.69 | 0.66 |
| 32IL | S11 | tune | 730 | drawn | - | pass | 0.17 | 0.14 |
| 9SXD | S11 | tune | 627 | drawn | - | pass | 0.27 | 0.23 |
| 9YLD | S11 | tune | 341 | drawn | - | pass | 0.26 | 0.23 |

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
