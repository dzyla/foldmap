# Plan 29 fix record: fix-ligand

Written by: foldmap 0.1.0, git 552687f251ea, code fingerprint 43e1d11adf24d91b, scripts/app29_fix.py sha256 f9ba7e369dcc; Python 3.12.9 on Linux-6.18.33.1-microsoft-standard-WSL2-x86_64-with-glibc2.39; gemmi 0.7.5, numpy 2.5.3, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

**Verdict: better** (identical or better on every measure keeps the fix).

- worse, tune: none
- worse, check: none
- worse, third: none
- worse, na: none
- worse, bench: none
- worse, marks_on_loop: none
- worse, marks_crowded: none
- worse, marks_c9: none
- worse, marks_elsewhere: none
- better, tune: 9X4L: now passes A20
- better, check: none
- better, third: none
- better, na: none
- better, marks_on_loop: 9SXF, 29TB, 31FD, 9WNE
- better, marks_c9: 9SXF, 29TB, 31FD, 9X4L

## Tune half (per entry; the check half in its totals only)

| id | criteria failed | floor ok | svg changed | seconds |
|---|---|---|---|---|
| `9SXY` | - -> - | True -> True | no | 0.27 -> 0.25 |
| `9SXU` | - -> - | True -> True | no | 0.28 -> 0.25 |
| `9SXF` | - -> - | True -> True | yes | 0.25 -> 0.24 |
| `9SXG` | - -> - | True -> True | no | 0.27 -> 0.25 |
| `29LB` | - -> - | True -> True | no | 0.24 -> 0.23 |
| `29TB` | - -> - | True -> True | yes | 0.36 -> 0.35 |
| `9T57` | - -> - | True -> True | no | 0.34 -> 0.34 |
| `9ZFX` | - -> - | True -> True | no | 0.61 -> 0.61 |
| `14CU` | - -> - | True -> True | no | 1.26 -> 1.28 |
| `30CY` | - -> - | True -> True | no | 0.94 -> 0.96 |
| `44LS` | - -> - | True -> True | no | 1.22 -> 1.24 |
| `9YL9` | - -> - | True -> True | no | 0.64 -> 0.65 |
| `20XT` | - -> - | True -> True | no | 0.75 -> 0.73 |
| `9WVR` | - -> - | True -> True | no | 0.28 -> 0.29 |
| `9Z0S` | - -> - | True -> True | no | 1.32 -> 1.31 |
| `9TCF` | - -> - | True -> True | no | 0.92 -> 0.94 |
| `9YOA` | - -> - | True -> True | no | 3.68 -> 3.63 |
| `38NR` | - -> - | True -> True | no | 0.28 -> 0.27 |
| `23FR` | - -> - | True -> True | no | 17.26 -> 17.24 |
| `31FD` | - -> - | False -> False | yes | 13.14 -> 13.05 |
| `22VF` | - -> - | True -> True | no | 1.52 -> 1.51 |
| `12XK` | - -> - | True -> True | no | 1.19 -> 1.2 |
| `35VC` | - -> - | True -> True | no | 1.01 -> 1.14 |
| `23CX` | - -> - | True -> True | no | 1.29 -> 1.28 |
| `28HY` | - -> - | True -> True | no | 0.47 -> 0.47 |
| `26EL` | - -> - | True -> True | no | 1.59 -> 1.63 |
| `9X4L` | A20 -> - | True -> True | yes | 16.83 -> 16.52 |
| `9NYA` | - -> - | True -> True | no | 0.14 -> 0.14 |
| `45JF` | - -> - | True -> True | no | 0.11 -> 0.11 |
| `34XH` | - -> - | True -> True | no | 0.14 -> 0.13 |
| `9WXY` | - -> - | True -> True | no | 0.05 -> 0.05 |
| `9ZSJ` | - -> - | True -> True | no | 2.7 -> 2.64 |
| `9YPK` | - -> - | True -> True | no | 23.66 -> 23.77 |
| `9YPN` | - -> - | True -> True | no | 20.87 -> 20.87 |
| `9WNE` | - -> - | True -> True | yes | 1.47 -> 1.48 |
| `22MZ` | - -> - | True -> True | no | 0.36 -> 0.43 |
| `9XIC` | - -> - | True -> True | no | 9.77 -> 9.65 |
| `9TMF` | - -> - | True -> True | no | 3.37 -> 3.4 |
| `34QD` | - -> - | True -> True | no | 0.82 -> 0.78 |
| `12KF` | - -> - | True -> True | no | 1.52 -> 1.48 |
| `30CU` | - -> - | True -> True | no | 0.95 -> 0.94 |
| `30CV` | - -> - | True -> True | no | 1.01 -> 1.03 |
| `9YLA` | - -> - | True -> True | no | 0.68 -> 0.69 |
| `9YLG` | - -> - | True -> True | no | 0.68 -> 0.69 |
| `32IL` | - -> - | True -> True | no | 0.18 -> 0.17 |
| `9SXD` | - -> - | True -> True | no | 0.27 -> 0.27 |
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
| floor_failed | 1 | 1 |
| ui_run | 0 | 0 |
| ui_failed | 0 | 0 |

## Third sample, default view (keys: pass, failing criteria, named fallbacks, area)

- 100 entries; pass 82 -> 82; verdict identical
- worse: none
- better: none
- moved: none
- SVG moved (no key reads it): 3ABR, 3G6E, 3WAU, 5S4U, 6IZP, 6QIK, 6VBG, 7F5S, 7XBK, 8EM7
- the committed report (views all, an older code) against the after run, default view: for reference: pass 83 -> 82; 36 worse, 18 better, 10 moved

## NA cases and controls (plan 28's NA measure)

- 46 entries; verdict identical
- worse: none
- better: none
- moved: none
- SVG moved (no key reads it): 1KX5, 6IZP, 6QIK

## bench --check --exact

- 22 entries; exit 0 -> 0
- changed: none

## Ligand marks

Over the 47 tune entries measured both ways: marks 289 -> 289, on a loop line 13 -> 6, crowded 0 -> 0; C9 (the tethers' length over bodies and along loop lines, page units: render._mark_spot's ranking) 625.147 -> 618.28, of it over bodies 594.605 -> 587.131; marks moved in a figure that did not change: none.

The third sample's and the NA figures that moved (11, drawn the app's way): marks 184 -> 184, on a loop line 19 -> 5, crowded 0 -> 0; C9 248.402 -> 220.908, of it over bodies 209.243 -> 200.269.

- `1KX5`: on a loop line 1 -> 0; C9 0.0 -> 0.0 (over bodies 0.0 -> 0.0, along loop lines 0.0 -> 0.0); 1 moved
- `3ABR`: on a loop line 1 -> 0; C9 7.466 -> 5.989 (over bodies 6.686 -> 5.989, along loop lines 0.78 -> 0.0); 1 moved
- `3G6E`: on a loop line 4 -> 0; C9 19.022 -> 9.601 (over bodies 9.527 -> 8.802, along loop lines 9.495 -> 0.799); 4 moved
- `3WAU`: on a loop line 1 -> 0; C9 5.049 -> 3.014 (over bodies 5.049 -> 3.014, along loop lines 0.0 -> 0.0); 1 moved
- `5S4U`: on a loop line 1 -> 0; C9 21.882 -> 20.989 (over bodies 20.973 -> 20.08, along loop lines 0.909 -> 0.909); 1 moved
- `6IZP`: on a loop line 0 -> 0; C9 0.0 -> 0.0 (over bodies 0.0 -> 0.0, along loop lines 0.0 -> 0.0); 1 moved
- `6QIK`: on a loop line 1 -> 0; C9 12.817 -> 7.596 (over bodies 6.529 -> 7.354, along loop lines 6.288 -> 0.242); 1 moved
- `6VBG`: on a loop line 1 -> 0; C9 20.134 -> 19.823 (over bodies 20.134 -> 19.823, along loop lines 0.0 -> 0.0); 2 moved
- `7F5S`: on a loop line 2 -> 1; C9 13.993 -> 13.993 (over bodies 12.491 -> 12.491, along loop lines 1.502 -> 1.502); 1 moved
- `7XBK`: on a loop line 6 -> 4; C9 51.738 -> 44.071 (over bodies 48.148 -> 43.086, along loop lines 3.59 -> 0.985); 3 moved
- `8EM7`: on a loop line 1 -> 0; C9 96.301 -> 95.832 (over bodies 79.706 -> 79.63, along loop lines 16.595 -> 16.202); 1 moved

The tune entries whose figure changed:

| id | marks | on a loop line | C9 | longest tether to another chain | moved |
|---|---|---|---|---|---|
| `9SXF` | 1 -> 1 | 1 -> 0 | 4.6 -> 4.435 | None -> None | 1 |
| `29TB` | 3 -> 3 | 1 -> 0 | 14.378 -> 13.368 | None -> None | 1 |
| `31FD` | 108 -> 108 | 4 -> 0 | 318.723 -> 313.709 | 6.227 -> 6.227 | 12 |
| `9X4L` | 30 -> 30 | 0 -> 0 | 41.962 -> 41.284 | 7.885 -> 3.115 | 1 |
| `9WNE` | 3 -> 3 | 1 -> 0 | 2.23 -> 2.23 | None -> None | 1 |

Each moved mark: its place, the line it holds, its tethers over bodies and along loop lines (page units), and its tethers' lengths, before -> after.

- `9SXF` A1JRG:A1001: [-4.868, -1.36] -> [-3.749, -1.395]; on a loop line 0.54 -> 0.0; over bodies 4.6 -> 4.435; along loop lines 0.0 -> 0.0; tethers [5.171, 5.644, 5.182, 4.773] -> [6.277, 6.129, 6.295, 3.729]
- `29TB` A1J4K:A306: [2.269, 3.76] -> [2.666, 3.814]; on a loop line 0.54 -> 0.0; over bodies 11.759 -> 10.749; along loop lines 0.0 -> 0.0; tethers [0.97, 7.252, 4.307, 7.184] -> [1.219, 7.568, 4.689, 7.312]
- `31FD` CLA:C509: [-10.642, 7.717] -> [-12.882, 13.121]; on a loop line 0.54 -> 0.0; over bodies 6.504 -> 3.515; along loop lines 0.0 -> 0.0; tethers [1.904, 4.219, 8.081, 9.276] -> [6.761, 8.22, 5.897, 13.341]
- `31FD` CLA:C511: [-9.732, 10.684] -> [-9.732, 7.287]; on a loop line 0.0 -> 0.0; over bodies 6.549 -> 7.683; along loop lines 0.771 -> 0.0; tethers [3.693, 8.5, 7.089, 5.919] -> [2.354, 8.336, 7.814, 8.704]
- `31FD` CLA:C516: [-9.778, 6.974] -> [-9.928, 3.867]; on a loop line 0.0 -> 0.0; over bodies 0.717 -> 0.727; along loop lines 0.0 -> 0.0; tethers [2.761, 2.108, 4.27, 4.186] -> [2.59, 5.022, 1.509, 1.087]
- `31FD` CLA:c516: [104.462, 5.945] -> [106.889, 9.717]; on a loop line 0.54 -> 0.0; over bodies 5.422 -> 4.156; along loop lines 0.0 -> 0.0; tethers [6.211, 6.781, 8.643] -> [8.431, 4.46, 4.5]
- `31FD` PHO:a409: [112.324, 7.286] -> [112.238, 6.739]; on a loop line 0.0 -> 0.0; over bodies 9.319 -> 10.327; along loop lines 0.909 -> 1.975; tethers [7.224, 3.512, 24.63, 16.61] -> [7.38, 3.23, 24.64, 16.481]
- `31FD` CLA:c511: [102.467, 7.743] -> [105.467, 14.208]; on a loop line 0.54 -> 0.0; over bodies 7.118 -> 3.19; along loop lines 0.0 -> 0.0; tethers [2.314, 8.532, 7.574, 8.432] -> [7.374, 7.862, 10.663, 1.536]
- `31FD` CLA:c507: [102.62, 12.895] -> [102.991, 13.961]; on a loop line 0.0 -> 0.0; over bodies 0.0 -> 0.0; along loop lines 0.0 -> 0.0; tethers [1.326, 4.928, 2.765] -> [1.716, 6.049, 1.769]
- `31FD` CLA:c505: [101.939, 5.32] -> [102.731, 7.366]; on a loop line 0.54 -> 0.0; over bodies 0.385 -> 0.713; along loop lines 0.0 -> 0.0; tethers [1.517, 3.902, 3.013, 3.808] -> [3.045, 1.765, 4.642, 4.907]
- `31FD` CLA:c509: [107.146, 5.888] -> [107.297, 3.901]; on a loop line 0.0 -> 0.0; over bodies 6.901 -> 4.568; along loop lines 0.0 -> 0.0; tethers [7.081, 3.968, 8.747, 4.698] -> [7.989, 3.766, 8.675, 3.857]
- `31FD` CLA:c508: [107.942, 4.641] -> [112.371, 4.029]; on a loop line 0.0 -> 0.0; over bodies 3.066 -> 4.337; along loop lines 0.0 -> 0.312; tethers [3.064, 6.639, 5.315, 3.783] -> [1.415, 2.913, 8.07, 7.833]
- `31FD` CLA:c503: [107.953, 10.057] -> [109.652, 12.054]; on a loop line 0.0 -> 0.0; over bodies 1.0 -> 1.192; along loop lines 0.0 -> 0.0; tethers [1.731, 6.046, 5.636] -> [1.36, 5.976, 4.417]
- `31FD` CLA:c515: [109.521, 11.601] -> [112.121, 8.693]; on a loop line 0.0 -> 0.0; over bodies 0.0 -> 0.952; along loop lines 0.0 -> 0.0; tethers [4.302, 1.458, 1.635] -> [1.634, 5.331, 2.417]
- `9X4L` ADP:32001: [0.487, 21.602] -> [5.825, 19.393]; on a loop line 0.0 -> 0.0; over bodies 6.568 -> 5.89; along loop lines 0.0 -> 0.0; tethers [2.418, 8.493, 8.218, 7.885] -> [4.467, 13.881, 2.919, 3.115]
- `9WNE` MG:A1302: [-11.932, 2.923] -> [-7.576, 3.606]; on a loop line 0.6 -> 0.0; over bodies 0.0 -> 0.0; along loop lines 0.0 -> 0.0; tethers [0.788, 0.982] -> [4.75, 3.65]

## Pictures (gitignored, for looking)

- `9SXF`: before docs/superpowers/work/app29/fix-ligand/png/9SXF-before.png, after docs/superpowers/work/app29/fix-ligand/png/9SXF-after.png, pair docs/superpowers/work/app29/fix-ligand/png/9SXF-pair.png
- `29TB`: before docs/superpowers/work/app29/fix-ligand/png/29TB-before.png, after docs/superpowers/work/app29/fix-ligand/png/29TB-after.png, pair docs/superpowers/work/app29/fix-ligand/png/29TB-pair.png
- `31FD`: before docs/superpowers/work/app29/fix-ligand/png/31FD-before.png, after docs/superpowers/work/app29/fix-ligand/png/31FD-after.png, pair docs/superpowers/work/app29/fix-ligand/png/31FD-pair.png
- `9X4L`: before docs/superpowers/work/app29/fix-ligand/png/9X4L-before.png, after docs/superpowers/work/app29/fix-ligand/png/9X4L-after.png, pair docs/superpowers/work/app29/fix-ligand/png/9X4L-pair.png
- `9WNE`: before docs/superpowers/work/app29/fix-ligand/png/9WNE-before.png, after docs/superpowers/work/app29/fix-ligand/png/9WNE-after.png, pair docs/superpowers/work/app29/fix-ligand/png/9WNE-pair.png
