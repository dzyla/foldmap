# foldmap random-100 campaign

85 of 101 pass.

Report written by: foldmap 0.2.0, git e1392f79daf7, code fingerprint 7879eea05bd9301c, scripts/random100.py sha256 b552d3c3d4b3; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Entries measured by 1 code version(s): 7879eea05bd9301c+b552d3c3d4b358ca (101)

An entry passes when no criterion fails. A criterion with nothing to judge in an entry is not judged (or does not apply), never counted as a pass: see the outcomes per criterion below.

Sample: seed 20261025, release-date cutoff 2026-10-13.

| stratum | entries | pass |
|---|---|---|
| S1 | 12 | 12 |
| S2 | 12 | 11 |
| S3 | 10 | 10 |
| S4 | 15 | 13 |
| S5 | 8 | 7 |
| S6 | 8 | 6 |
| S7 | 7 | 5 |
| S8 | 7 | 2 |
| S9 | 3 | 2 |
| S10 | 5 | 5 |
| S11 | 10 | 8 |
| S12 | 3 | 3 |
| A1c | 1 | 1 |

## Outcomes per criterion

✓ pass, ✗ fail, ○ not judged (the criterion applies but the entry gave it nothing to judge: no tie, no tether to another chain, no domain drawn, mkdssp unmeasured ...), - does not apply.

| criterion | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| A1 no crash: exit 0, files written, no traceback, no Python warning | 100 | 0 | 0 | 1 |
| A1c a CA-only or P-only entry exits non-zero with foldmap's one-line message | 1 | 0 | 0 | 100 |
| A2 default view time | 100 | 0 | 0 | 1 |
| A2b the app draws it within the hosted limits | 99 | 1 | 0 | 1 |
| A3 no route fallbacks | 93 | 7 | 0 | 1 |
| A4 no loop over a label | 99 | 1 | 0 | 1 |
| A5 no element or stem boxes overlapping | 100 | 0 | 0 | 1 |
| A6 no text collisions | 93 | 7 | 0 | 1 |
| A7 text and legend inside the page | 100 | 0 | 0 | 1 |
| A8 labels unique per chain, stem names per figure | 100 | 0 | 0 | 1 |
| A9 nothing left out without a note | 100 | 0 | 0 | 1 |
| A10 secondary structure agrees with mkdssp | 97 | 0 | 0 | 4 |
| A11 N->C and sheet order | 97 | 0 | 0 | 4 |
| A12 same figure for any PYTHONHASHSEED | 100 | 0 | 0 | 1 |
| A13 a saved layout redraws the same figure | 100 | 0 | 0 | 1 |
| A14 interactive page and sequence plot written | 100 | 0 | 0 | 1 |
| A15 copies drawn alike | 40 | 1 | 14 | 46 |
| A16 domains: one chain, alike, not tiny, panels apart, names clear | 50 | 1 | 8 | 42 |
| A17 nucleic-acid metrics and stem names by their letters | 9 | 9 | 0 | 83 |
| A18 membrane: only the band set in the band (TM helices across it, slab helices at their depths, crossing strands), no element nearer than its residues or on the wrong side, no loop nearer than its residues but out of a band element's port to its face or, for a loop in the bilayer, beside its own port to its residues' depth; a loop joining the faces crosses the band only beside its own port or past the band's ends | 6 | 2 | 0 | 93 |
| A19 two states: same names in the same order | 8 | 2 | 0 | 91 |
| A20 ligand tethers local | 8 | 0 | 54 | 39 |

### Not judged, per criterion

- **A15** (14): 5S58 (copies drawn with no symmetry tie: default C~A, D~B); 6BYU (copies drawn with no symmetry tie: default B~A; focus B~A); 9Y9L (copies drawn with no symmetry tie: default B~A, C~A, D~A, E~A, b~a, c~a, d~a); 8C5C (copies drawn with no symmetry tie: default b~a, c~a, d~a, e~a, f~a, g~a, h~a, i~a, j~a, k~a, l~a, m~a, O~N, P~N, Q~N, R~N, S~N, T~N, U~N, V~N, W~N, X~N, Y~N, Z~); 7YF0 (copies drawn with no symmetry tie: default 2~1, 3~1, 4~1, 5~1, B~A, C~A, D~A, E~A, I~H, J~H, K~H, L~H, a~A, b~A, c~A, d~A, e~A); 2MF1 (copies drawn with no symmetry tie: default B~A, C~A, D~A, E~A, F~A); 6R0Z (copies drawn with no symmetry tie: default B~A, C~A, E~D, F~D, K~I, L~J, P~O, Q~O, R~O, S~O, T~O, U~O, V~O, W~O, X~O, Y~O, Z~O); 7ZDM (copies drawn with no symmetry tie: default X~j); 8W2Q (copies drawn with no symmetry tie: default B~A); 9QNC (copies drawn with no symmetry tie: default C~A, E~A, G~A, I~A, K~A); 9QMS (copies drawn with no symmetry tie: default G~C, H~D, E~A, F~B, c~C, d~D, a~A, b~B, o~O, N~L, g~C, h~D, e~A, f~B, M~K, k~K, n~l, m~K; focus G~C, H~D, E~A, F~B, c); 9ZP3 (copies drawn with no symmetry tie: default G~A); 7EB3 (copies drawn with no symmetry tie: default B~A, C~A); 9W9N (copies drawn with no symmetry tie: default B~A, C~A)
- **A16** (8): 4OXB (the domains view drew no domain); 5ZJU (the domains view drew no domain); 6T4C (the domains view drew no domain); 8KE9 (the domains view drew no domain); 8FTF (the domains view drew no domain); 2MF1 (the domains view drew no domain); 9W9N (the domains view drew no domain); 2H4E (the domains view drew no domain)
- **A20** (54): 4DEH (no tether to another chain; reported: ligand-tether:0JK:A1401:2 7.333 (home A)); 7KUV (no tether to another chain); 8EYE (no tether to another chain); 1JMH (no tether to another chain; reported: ligand-tether:UMP:A2317:3 11.328 (home A2); ligand-tether:UMP:A1317:3 11.328 (home A1)); 7VOD (no tether to another chain; reported: ligand-tether:7RU:A1208:2 8.242 (home A)); 2ATZ (no tether to another chain; reported: ligand-tether:DGT:A201:2 6.752 (home A)); 9AZ7 (no tether to another chain; reported: ligand-tether:HC4:A201:2 7.025 (home A)); 6NJW (no tether to another chain); 6HJZ (no tether to another chain); 3SXN (no tether to another chain); 6QQQ (no tether to another chain); 5QQQ (no tether to another chain); 4MIN (no tether to another chain); 1IVX (no tether to another chain); 6BWQ (no tether to another chain; reported: ligand-tether:CTP:B1504:2 8.002 (home B1); ligand-tether:CTP:B1504:3 9.541 (home B1); ligand-tether:CTP:A1504:0 7.827 (hom); 2FYN (no tether to another chain; reported: ligand-tether:LOP:A504:2 8.915 (home A); ligand-tether:LOP:A504:3 7.904 (home A); ligand-tether:HEM:A501:2 8.145 (home A);); 5S58 (no tether to another chain; reported: ligand-tether:GTP:A501:3 10.525 (home A); ligand-tether:GTP:C501:3 8.615 (home C); ligand-tether:NUY:D503:1 7.785 (home D)); 6BYU (no tether to another chain; reported: ligand-tether:ECJ:C1401:0 7.829 (home C); ligand-tether:ECJ:C1401:1 7.829 (home C); ligand-tether:ZN:D2003:2 7.111 (home D); 5CX1 (no tether to another chain; reported: ligand-tether:HCA:A501:2 8.151 (home A); ligand-tether:HCA:A501:3 6.72 (home A); ligand-tether:ICS:A502:1 7.771 (home A); ); 7TCP (no tether to another chain); 9Y9L (no tether to another chain; reported: ligand-tether:ADP:A400:0 7.712 (home A); ligand-tether:ADP:A400:1 7.783 (home A); ligand-tether:ADP:B400:1 8.503 (home B);); 6O20 (no tether to another chain); 9VK1 (no tether to another chain; reported: ligand-tether:3VV:A601:2 11.54 (home A); ligand-tether:3VV:B601:2 11.54 (home B)); 8C5C (no tether to another chain; reported: ligand-tether:GDP:a602:1 6.627 (home a); ligand-tether:GDP:a602:3 6.669 (home a); ligand-tether:TA1:a601:0 14.873 (home a)); 8UXG (no tether to another chain); 7YF0 (no tether to another chain); 2KZ2 (no tether to another chain); 2EGQ (no tether to another chain); 2JZI (no tether to another chain); 6R0Z (no tether to another chain; reported: ligand-tether:ADP:A603:0 7.448 (home A); ligand-tether:ADP:A603:1 6.863 (home A); ligand-tether:ADP:A603:2 6.835 (home A);); 8TPL (no tether to another chain); 5EQI (no tether to another chain; reported: ligand-tether:5RH:A501:3 8.21 (home A)); 5V5Z (no tether to another chain); 8W2Q (no tether to another chain; reported: ligand-tether:SAM:B1001:0 7.182 (home B); ligand-tether:SAH:A1001:1 6.841 (home A)); 9QNC (no tether to another chain; reported: ligand-tether:ATP:E402:2 6.537 (home E); ligand-tether:ATP:E402:3 6.747 (home E); ligand-tether:ATP:G402:2 9.681 (home G);); 6PL8 (no tether to another chain; reported: ligand-tether:1FZ:A501:0 8.628 (home A)); 9CJ9 (no tether to another chain); 6TBZ (no tether to another chain); 9ZP3 (no tether to another chain; reported: ligand-tether:A1C49:J1503:0 9.746 (home J); ligand-tether:A1C49:J1503:1 14.792 (home J); ligand-tether:A1C49:J1503:2 19.2 ); 8U8V (no tether to another chain; reported: ligand-tether:APC:E1301:3 13.553 (home E)); 9NL4 (no tether to another chain; reported: ligand-tether:ZN:A2002:3 7.393 (home A)); 8SSW (no tether to another chain); 7SA4 (no tether to another chain; reported: ligand-tether:ZN:a101:0 6.556 (home a); ligand-tether:ZN:a101:1 6.652 (home a); ligand-tether:ZN:a101:2 17.773 (home a)); 31FZ (no tether to another chain); 6D8D (no tether to another chain; reported: ligand-tether:GAL+NAG+SIA:G2:0 7.514 (home A); ligand-tether:GAL+NAG+SIA:G2:1 6.922 (home A); ligand-tether:GAL+NAG+SIA:G2); 9N3U (no tether to another chain; reported: ligand-tether:CU:A1002:2 8.61 (home A)); 6ZQ1 (no tether to another chain; reported: ligand-tether:EDG:AAA701:0 10.333 (home AAA); ligand-tether:EDG:AAA701:3 12.794 (home AAA)); 8JX8 (no tether to another chain); 9JTM (no tether to another chain); 9J7U (no tether to another chain); 8V9G (no tether to another chain; reported: ligand-tether:MER:B301:3 7.645 (home B)); 1W1P (no tether to another chain); 9S44 (no tether to another chain); 5DY5 (no tether to another chain; reported: ligand-tether:5GR:A1402:0 10.253 (home A1); ligand-tether:5GR:A2402:0 11.469 (home A2); ligand-tether:5GR:A2402:1 8.327 (h)

## Entries

process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 gates; drawing s: the drawing inside it.

| id | stratum | bin | chains | process s | drawing s | A1 | A1c | A2 | A2b | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 | A11 | A12 | A13 | A14 | A15 | A16 | A17 | A18 | A19 | A20 | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4DEH | S1 | M | 1 | 1.91 | 0.5 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 7KUV | S1 | M | 1 | 2.36 | 0.89 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 8EYE | S1 | M | 1 | 2.66 | 1.1 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 1JMH | S1 | M | 1 | 2.96 | 1.38 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7VOD | S1 | M | 1 | 1.96 | 0.54 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 9BNW | S1 | M | 1 | 2.91 | 1.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 2ATZ | S1 | S | 1 | 1.71 | 0.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2CXH | S1 | S | 1 | 1.61 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9AZ7 | S1 | S | 1 | 1.72 | 0.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3MPL | S1 | S | 1 | 1.81 | 0.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6NJW | S1 | S | 1 | 1.76 | 0.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2EY5 | S1 | S | 1 | 1.67 | 0.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6HJZ | S2 | L | 10 | 4.22 | 2.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6NSY | S2 | L | 4 | 2.46 | 0.94 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ |  |
| 3ZC1 | S2 | L | 8 | 1.76 | 0.43 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - |  |
| 3SXN | S2 | L | 6 | 2.16 | 0.74 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ○ |  |
| 6QQQ | S2 | M | 2 | 2.91 | 1.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 5QQQ | S2 | M | 2 | 3.37 | 1.77 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 4MIN | S2 | M | 4 | 4.72 | 2.97 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 1IVX | S2 | M | 2 | 3.72 | 2.04 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6BWQ | S2 | S | 2 | 4.72 | 2.99 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 4OXB | S2 | S | 2 | 4.52 | 3.07 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5ZJU | S2 | XL | 60 | 4.22 | 2.23 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - | A3: focus: C:157-161>C:164-165, F:173-174>F:179-181, G:157-161>G:164-165, I:173-174>I:185-186 |
| 3ANZ | S2 | XL | 28 | 2.21 | 0.93 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 2FYN | S3 | L | 18 | 3.71 | 1.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 4G4S | S3 | L | 16 | 6.57 | 4.22 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 5S58 | S3 | L | 6 | 5.22 | 3.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 3BLV | S3 | L | 8 | 2.76 | 1.23 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 6T4C | S3 | M | 4 | 2.12 | 0.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ✓ |  |
| 5TDO | S3 | M | 4 | 2.71 | 1.6 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 1TZY | S3 | M | 8 | 2.06 | 0.9 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 4LZF | S3 | S | 2 | 1.16 | 0.18 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6BYU | S3 | XL | 12 | 6.93 | 4.77 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 5CX1 | S3 | XL | 16 | 4.6 | 2.83 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8VU3 | S4 | L | 36 | 9.69 | 5.77 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ |  |
| 7TCP | S4 | L | 8 | 2.21 | 0.9 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9Y9L | S4 | L | 11 | 5.56 | 3.7 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 6O20 | S4 | L | 6 | 4.43 | 2.92 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9CV6 | S4 | L | 2 | 3.67 | 2.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 9VK1 | S4 | M | 2 | 2.26 | 1.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8KE9 | S4 | M | 3 | 1.66 | 0.56 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 8C5C | S4 | XL | 52 | 71.54 | 55.93 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✗ | - | - | - | ○ | A3: domains: c:198-202>c:205-214, c:250-255>c:265-271, c:265-271>c:286-294, c:305-307>c:310-318; A6: domains 5 (domain-label:D2 · c × loop:c:198-202>c:205-214; domain-label:D1 · i × loop:S:405-409>S:416-435; domain-label:D2 · S × loop:c:396-399>c:405-426); A16: domains 3 (domain-label:D2 · c × loop:c:198-202>c:205-214; domain-label:D1 · i × loop:S:405-409>S:416-435; domain-label:D2 · S × loop:c:396-399>c:405-426) |
| 8FTF | S4 | XL | 33 | 2.06 | 0.74 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 8UXG | S4 | XL | 8 | 11.47 | 8.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8WLE | S4 | XL | 52 | 4.37 | 2.24 | ✓ | - | ✓ | ✓ | ✗ | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | - | - | A3: focus: A:186-188>A:194-196, A:197-199>A:204-210, B:186-188>B:194-196, C:186-188>C:194-196; A4: focus: A:194-196>A:197-199, B:124-155>B:161-171, B:194-196>B:197-199, C:124-155>C:161-171; A6: focus 29 (terminus:N:B × loop:D:124-155>D:161-171; terminus:N:D × loop:C:194-196>C:197-199; terminus:N:G × loop:I:124-155>I:161-171); A15: focus 83.05590925422942 |
| 7YF0 | S4 | XL | 22 | 81.37 | 73.83 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 8OPJ | S4 | XL | 48 | 2.86 | 1.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - |  |
| 8VJJ | S4 | XL | 8 | 14.82 | 11.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ |  |
| 8FM9 | S4 | XL | 24 | 3.12 | 1.53 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 2MF1 | S5 | M | 7 | 9.67 | 1.08 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | ✗ | - | - | - | A17: procrustes 0.535 (limit 0.5) |
| 2MSV | S5 | S | 1 | 1.21 | 0.24 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2KZ2 | S5 | S | 1 | 1.25 | 0.23 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2EGQ | S5 | S | 1 | 1.21 | 0.22 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2LLX | S5 | S | 1 | 1.31 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2JZI | S5 | S | 2 | 1.51 | 0.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 2FJ3 | S5 | S | 1 | 1.41 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 1FVL | S5 | S | 1 | 1.06 | 0.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6R0Z | S6 | L | 26 | 9.75 | 7.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✓ | - | ○ |  |
| 8TPL | S6 | L | 7 | 2.06 | 0.83 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ○ |  |
| 8JD0 | S6 | L | 2 | 7.02 | 5.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ✓ |  |
| 8IAK | S6 | L | 8 | 4.55 | 3.03 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - |  |
| 5EQI | S6 | M | 1 | 1.71 | 0.65 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 6C1E | S6 | M | 2 | 2.51 | 1.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ |  |
| 5V5Z | S6 | M | 1 | 1.92 | 0.77 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✗ | - | ○ | A18: no band drawn |
| 7ZDM | S6 | XL | 45 | 37.55 | 32.91 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✗ | - | ✓ | A3: membrane: d:192-194>d:195-197; A6: membrane 2 (label:1:231-237 × loop:d:192-194>d:195-197; label:d:207-219 × loop:d:192-194>d:195-197); A18: TM helix A:56-72 0.80; loop drawn nearer the mid-plane than its residues d:192-194>d:195-197 1.7 Å |
| 8W2Q | S7 | L | 5 | 8.5 | 4.06 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✗ | - | - | ○ | A3: domains: strand E:21-22; A6: default 1 (na:nt:D:12 × stub:C:A); deposited 1 (na:nt:D:12 × stub:C:A); beta 1 (na:nt:D:12 × stub:C:A); domains 1 (na:nt:D:32 × na:run:5); A17: pair_f1 0.741 (limit 0.95); stem_recall 0.667 (limit 0.85) |
| 9QNC | S7 | L | 13 | 10.4 | 6.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | - | - | ○ |  |
| 7PSX | S7 | M | 12 | 3.22 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 6PL8 | S7 | M | 3 | 3.78 | 0.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 9CJ9 | S7 | M | 3 | 3.71 | 0.7 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 6TBZ | S7 | S | 3 | 3.24 | 0.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 9QMS | S7 | XL | 30 | 53.63 | 40.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✗ | - | - | - | A6: default 2 (na:tie:26-127:label × stub:N:C; terminus:C:F × na:tie:58-96:mark); deposited 2 (na:tie:26-127:label × stub:N:C; terminus:C:F × na:tie:58-96:mark); beta 2 (na:tie:26-127:label × stub:N:C; terminus:C:F × na:tie:58-96:mark); focus 2 (na:tie:26-127:label × stub:N:C; terminus:C:F × na:tie:58-96:mark); domains 1 (na:end:I:9-119:5′ × na:tie:9-148); A17: pair_f1 0.892 (limit 0.95); stem_recall 0.714 (limit 0.85) |
| 9ZP3 | S8 | L | 8 | 9.6 | 4.48 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | - | - | ○ |  |
| 8U8V | S8 | L | 6 | 5.57 | 2.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | - | - | ○ | A17: procrustes 0.518 (limit 0.5) |
| 9NL4 | S8 | M | 5 | 5.14 | 1.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | ○ | A17: pair_f1 0.865 (limit 0.95); stem_recall 0.667 (limit 0.85) |
| 8SSW | S8 | M | 4 | 4.67 | 1.6 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | - | - | ○ | A17: pair_f1 0.875 (limit 0.95) |
| 6JVX | S8 | S | 2 | 2.96 | 0.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 7SA4 | S8 | XL | 58 | 106.5 | 88.34 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ○ | A3: default: strand 1:836-836; deposited: strand 1:836-836; beta: strand 1:836-836; focus: strand 1:836-836; A6: default 3 (na:gap:513 × terminus:C:x; na:gap:513 × stub:C:x; terminus:N:T × na:run:93); deposited 3 (na:gap:513 × terminus:C:x; na:gap:513 × stub:C:x; terminus:N:T × na:run:93); beta 3 (na:gap:513 × terminus:C:x; na:gap:513 × stub:C:x; terminus:N:T × na:run:93); focus 3 (na:gap:513 × terminus:C:x; na:gap:513 × stub:C:x; terminus:N:T × na:run:93); A17: stems_nearer_frac 0.994 (limit 1.0); fallbacks 1 (limit 0); na_crossings 3 (limit 0) |
| 31FZ | S8 | XL | 78 | 156.18 | 125.87 | ✓ | - | ✓ | ✗ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ○ | A2b: time: 31FZ took longer than 90 s to draw, the limit of this app. Draw part of it (Chains & assembly), or on your own computer: pip install foldmap, then foldmap plot 31FZ -o 31FZ.svg.; A3: default: strand BC:867-876, strand CI:1430-1431, strand CI:3720-3730, gap after CI:3599; deposited: strand BC:867-876, strand CI:1430-1431, strand CI:3720-3730, gap after CI:3599; beta: strand BC:867-876, strand CI:1430-1431, strand CI:3720-3730, gap after CI:3599; focus: strand BC:867-876, strand CI:1430-1431, strand CI:3720-3730, gap after CI:3599; A6: default 3 (na:tie:928-1015:label × terminus:N:Bf; |
| 9ZHB | S9 | M | 1 | 5.64 | 1.68 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 6WLQ | S9 | S | 1 | 5.64 | 0.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 7EDN | S9 | S | 1 | 2.77 | 0.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✗ | - | - | - | A17: pair_f1 0.0 (limit 0.95) |
| 6D8D | S10 | L | 6 | 4.27 | 2.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7EB3 | S10 | L | 3 | 25.63 | 23.99 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | - |  |
| 9N3U | S10 | M | 1 | 3.22 | 1.99 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6ZQ1 | S10 | M | 1 | 1.76 | 0.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 8JX8 | S10 | XL | 10 | 35.87 | 33.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9JTM | S11 | M | 1 | 1.56 | 0.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 9J7U | S11 | M | 1 | 1.81 | 0.75 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 8V9G | S11 | M | 2 | 2.16 | 0.92 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ○ |  |
| 6TS9 | S11 | M | 2 | 1.71 | 0.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | - |  |
| 9W9N | S11 | M | 3 | 1.31 | 0.34 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | - | - | ✓ | - |  |
| 2H4E | S11 | S | 2 | 1.81 | 0.81 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | ✓ | - |  |
| 1GPF | S11 | M | 2 | 2.81 | 1.51 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - |  |
| 1W1P | S11 | M | 2 | 1.91 | 0.76 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 9S44 | S11 | M | 1 | 1.4 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✗ | ○ | A19: names compared on all 21 pairs, 8 differ; αN~α1 (A:64-72~A1:64-72); α1′~α3 (A:121-126~A1:121-126); α1″~α4 (A:129-137~A1:129-138); α1~α5 (A:147-157~A1:147-157); α2~α6 (A:172-175~A1:172-175); α5~α7 (A:206-214~A1:206-214); α6~α8 (A:241-250~A1:241-253); α9′~α10 (A:338-354~A1:338-354) (+6 more) |
| 5DY5 | S11 | M | 1 | 2.51 | 1.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | ○ | A19: names compared on all 21 pairs, 8 differ; α1~αN (A1:64-72~A:64-72); α3~α1′ (A1:121-126~A:121-126); α4~α1″ (A1:129-138~A:129-137); α5~α1 (A1:147-157~A:147-157); α6~α2 (A1:172-175~A:172-175); α7~α5 (A1:206-214~A:206-214); α8~α6 (A1:241-253~A:241-250); α10~α9′ (A1:338-354~A:338-354) (+6 more) |
| AF_AFQ8I574F1 | S12 | M | 1 | 2.51 | 1.24 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFQ6PHI8F1 | S12 | M | 1 | 2.01 | 0.85 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFQ9UZ20F1 | S12 | S | 1 | 1.36 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2CSE | A1c |  | 27 | 0.95 | 0.11 | - | ✓ | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |  |

## Two-state pairs (S11)

| pair | seed | partner | sampled as | types now | ligand touching | assembly 1 chains | CA RMSD Å | A19 |
|---|---|---|---|---|---|---|---|---|
| P1 | 9JTM | 9J7U | apo/holo | apo/holo | [True, False] | [1, 1] | 0.972 | ✓ names compared on all 16 pairs; A:40-44 unpaired; A:280-282 unpaired; A:93-95 unpaired; A:236-238 unpaired; A:278-284 unpaired |
| P2 | 8V9G | 6TS9 | apo/holo | apo/holo, monomer/oligomer | [True, False] | [2, 1] | 0.28 | ✓ names compared on all 21 pairs; A:216-219 unpaired; A:216-219 unpaired |
| P3 | 9W9N | 2H4E | open/closed | open/closed | [False, False] | [3, 4] | 21.887 | ✓ names compared on all 6 pairs; B:19-21 unpaired; B:59-62 unpaired; B:65-68 unpaired; B:76-80 unpaired; A1:23-24 unpaired; A1:41-48 unpaired; A1:54-55 unpaired; A1:75-81 unpaired |
| P4 | 1GPF | 1W1P | monomer/oligomer | apo/holo, monomer/oligomer | [False, True] | [2, 1] | 0.273 | ✓ names compared on all 36 pairs; A:445-447 unpaired |
| P5 | 9S44 | 5DY5 | monomer/oligomer | apo/holo, monomer/oligomer | [False, True] | [1, 2] | 1.737 | ✗ names compared on all 21 pairs, 8 differ; αN~α1 (A:64-72~A1:64-72); α1′~α3 (A:121-126~A1:121-126); α1″~α4 (A:129-137~A1:129-138); α1~α5 (A:147-157~A1:147-157); α2~α6 (A:172-175~A1:172-175); α5~α7 (A:206-214~A1:206-214); α6~α8 (A:241-250~A1:241-253); α9′~α10 (A:338-354~A1:338-354) (+6 more) |

## Reported, not gated

| id | honest_crossings | loop_crossings | na_run_crossings | fallback_curves | fallback_curves_crossing_an_element | kink_links_rejected | tighter_than_inflate | rejected_as_crossing | names_moved | area | max_rss_mb |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4DEH | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 335.01 | 169.9 |
| 7KUV | 47 | 47 | None | 0 | 0 | 0 | 0 | 0 | 0 | 814.09 | 173.4 |
| 8EYE | 51 | 51 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1407.44 | 182.1 |
| 1JMH | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 650.0 | 212.5 |
| 7VOD | 5 | 5 | None | 0 | 0 | 0 | 0 | 0 | 0 | 443.14 | 170.1 |
| 9BNW | 36 | 36 | None | 0 | 0 | 0 | 0 | 0 | 0 | 674.29 | 204.1 |
| 2ATZ | 7 | 7 | None | 0 | 0 | 0 | 0 | 0 | 0 | 179.28 | 168.9 |
| 2CXH | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 178.38 | 167.3 |
| 9AZ7 | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 120.59 | 167.4 |
| 3MPL | 13 | 13 | None | 0 | 0 | 0 | 0 | 0 | 0 | 401.85 | 168.2 |
| 6NJW | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 295.38 | 168.9 |
| 2EY5 | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 142.53 | 167.7 |
| 6HJZ | 110 | 110 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2024.06 | 244.7 |
| 6NSY | 30 | 30 | None | 0 | 0 | 0 | 0 | 0 | 0 | 721.71 | 229.1 |
| 3ZC1 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 343.51 | 198.4 |
| 3SXN | 35 | 35 | None | 0 | 0 | 0 | 0 | 0 | 0 | 639.73 | 229.5 |
| 6QQQ | 29 | 29 | None | 0 | 0 | 0 | 0 | 0 | 0 | 520.37 | 196.0 |
| 5QQQ | 55 | 55 | None | 0 | 0 | 0 | 0 | 0 | 0 | 942.54 | 251.4 |
| 4MIN | 104 | 104 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2624.35 | 324.6 |
| 1IVX | 73 | 73 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1639.23 | 284.3 |
| 6BWQ | 81 | 81 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4716.54 | 368.1 |
| 4OXB | 30 | 30 | None | 0 | 0 | 0 | 0 | 0 | 0 | 816.45 | 173.5 |
| 5ZJU | 24 | 24 | None | 0 | 0 | 0 | 0 | 0 | 0 | 235.63 | 444.6 |
| 3ANZ | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 440.35 | 320.1 |
| 2FYN | 49 | 49 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1489.93 | 256.2 |
| 4G4S | 234 | 234 | None | 0 | 0 | 0 | 0 | 0 | 0 | 5493.25 | 341.7 |
| 5S58 | 149 | 149 | None | 0 | 0 | 0 | 1 | 0 | 0 | 4427.37 | 373.5 |
| 3BLV | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 968.34 | 255.5 |
| 6T4C | 80 | 80 | None | 0 | 0 | 0 | 1 | 0 | 0 | 913.74 | 175.6 |
| 5TDO | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 499.74 | 172.3 |
| 1TZY | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1133.42 | 193.7 |
| 4LZF | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 93.44 | 167.2 |
| 6BYU | 315 | 315 | None | 0 | 0 | 0 | 3 | 0 | 0 | 4497.83 | 455.9 |
| 5CX1 | 155 | 155 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2321.88 | 378.5 |
| 8VU3 | 77 | 77 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2880.33 | 354.7 |
| 7TCP | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 882.99 | 218.8 |
| 9Y9L | 193 | 193 | None | 0 | 0 | 0 | 1 | 0 | 0 | 4293.47 | 411.4 |
| 6O20 | 47 | 47 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2905.71 | 372.4 |
| 9CV6 | 204 | 204 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2189.15 | 277.2 |
| 9VK1 | 20 | 20 | None | 0 | 0 | 0 | 0 | 0 | 0 | 920.14 | 213.8 |
| 8KE9 | 42 | 42 | None | 0 | 0 | 0 | 0 | 0 | 0 | 954.14 | 172.6 |
| 8C5C | 2884 | 2884 | None | 0 | 0 | 0 | 5 | 0 | 0 | 28429.83 | 762.8 |
| 8FTF | 5 | 5 | None | 0 | 0 | 0 | 0 | 0 | 0 | 344.62 | 335.9 |
| 8UXG | 311 | 311 | None | 0 | 0 | 0 | 1 | 0 | 0 | 8316.78 | 497.5 |
| 8WLE | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1539.1 | 486.5 |
| 7YF0 | 1465 | 1465 | None | 0 | 0 | 0 | 5 | 0 | 0 | 21761.5 | 626.6 |
| 8OPJ | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 159.4 | 420.6 |
| 8VJJ | 267 | 267 | None | 0 | 0 | 0 | 0 | 0 | 0 | 9731.95 | 771.0 |
| 8FM9 | 19 | 19 | None | 0 | 0 | 0 | 1 | 0 | 0 | 862.8 | 352.4 |
| 2MF1 | 18 | 18 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 1028.52 | 379.6 |
| 2MSV | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 121.93 | 246.7 |
| 2KZ2 | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 74.34 | 196.1 |
| 2EGQ | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 36.24 | 173.9 |
| 2LLX | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 125.53 | 209.6 |
| 2JZI | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 196.94 | 250.7 |
| 2FJ3 | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 78.01 | 193.1 |
| 1FVL | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 14.52 | 168.1 |
| 6R0Z | 336 | 336 | None | 0 | 0 | 0 | 2 | 0 | 0 | 7955.85 | 465.5 |
| 8TPL | 25 | 25 | None | 0 | 0 | 0 | 0 | 0 | 0 | 951.35 | 251.6 |
| 8JD0 | 55 | 55 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1764.18 | 215.4 |
| 8IAK | 156 | 156 | None | 0 | 0 | 0 | 1 | 0 | 0 | 2968.06 | 306.9 |
| 5EQI | 13 | 13 | None | 0 | 0 | 0 | 0 | 0 | 0 | 677.06 | 177.1 |
| 6C1E | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1164.38 | 217.9 |
| 5V5Z | 30 | 30 | None | 0 | 0 | 0 | 0 | 0 | 0 | 582.46 | 172.2 |
| 7ZDM | 264 | 264 | None | 0 | 0 | 0 | 1 | 0 | 0 | 14294.85 | 532.3 |
| 8W2Q | 79 | 79 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2795.9 | 365.0 |
| 9QNC | 257 | 257 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3100.66 | 343.8 |
| 7PSX | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 154.44 | 179.3 |
| 6PL8 | 41 | 41 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 726.3 | 182.2 |
| 9CJ9 | 21 | 21 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 541.19 | 183.0 |
| 6TBZ | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 233.22 | 175.7 |
| 9QMS | 1190 | 1190 | 12 | 0 | 0 | 0 | 3 | 0 | 0 | 23752.21 | 557.3 |
| 9ZP3 | 282 | 282 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 4370.36 | 468.9 |
| 8U8V | 69 | 69 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 2063.6 | 261.9 |
| 9NL4 | 50 | 50 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2038.7 | 201.1 |
| 8SSW | 104 | 104 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1363.8 | 225.1 |
| 6JVX | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 101.18 | 172.1 |
| 7SA4 | 507 | 507 | 833 | 0 | 0 | 0 | 1 | 0 | 1 | 38765.61 | 687.0 |
| 31FZ | 1004 | 1004 | 603 | 0 | 0 | 0 | 3 | 0 | 1 | 43670.16 | 926.6 |
| 9ZHB | 0 | 0 | 17 | 0 | 0 | 0 | 0 | 0 | 0 | 1981.72 | 199.4 |
| 6WLQ | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 806.49 | 289.7 |
| 7EDN | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 117.39 | 170.5 |
| 6D8D | 131 | 131 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1291.05 | 203.8 |
| 7EB3 | 213 | 213 | None | 0 | 0 | 0 | 2 | 0 | 0 | 2849.84 | 348.9 |
| 9N3U | 41 | 41 | None | 0 | 0 | 0 | 0 | 0 | 0 | 707.64 | 177.0 |
| 6ZQ1 | 47 | 47 | None | 0 | 0 | 0 | 1 | 0 | 0 | 664.62 | 176.6 |
| 8JX8 | 97 | 97 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3316.41 | 525.3 |
| 9JTM | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 454.53 | 170.2 |
| 9J7U | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 582.19 | 169.7 |
| 8V9G | 24 | 24 | None | 0 | 0 | 0 | 0 | 0 | 0 | 603.27 | 209.1 |
| 6TS9 | 17 | 17 | None | 0 | 0 | 0 | 0 | 0 | 0 | 344.41 | 169.3 |
| 9W9N | 24 | 24 | None | 0 | 0 | 0 | 0 | 0 | 0 | 254.36 | 169.0 |
| 2H4E | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 398.13 | 171.0 |
| 1GPF | 93 | 93 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1187.18 | 271.2 |
| 1W1P | 38 | 38 | None | 0 | 0 | 0 | 2 | 0 | 0 | 616.65 | 174.0 |
| 9S44 | 15 | 15 | None | 0 | 0 | 0 | 1 | 0 | 0 | 370.28 | 170.6 |
| 5DY5 | 70 | 70 | None | 0 | 0 | 0 | 0 | 0 | 0 | 814.59 | 223.1 |
| AF_AFQ8I574F1 | 68 | 68 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1925.51 | 192.8 |
| AF_AFQ6PHI8F1 | 61 | 61 | None | 0 | 0 | 0 | 0 | 0 | 0 | 804.27 | 173.0 |
| AF_AFQ9UZ20F1 | 13 | 13 | None | 0 | 0 | 0 | 2 | 0 | 0 | 165.76 | 168.2 |
| 2CSE | None | None | None | 0 | 0 | 0 | 0 | 0 | 0 | None | 158.9 |

## Fallbacks taken

Each named fallback (foldmap.fallback.FALLBACKS, docs/validation/inventory/inventory.md) counted as the figures were made: reported, not gated.

| fallback | entries (default view) | times (default view) | times (every view) |
|---|---|---|---|
| loop_curve | 0 | 0 | 110 |
| loop_over_text | 0 | 0 | 17 |
| loop_wall_detour | 0 | 0 | 1 |
| loop_wall_unkept | 0 | 0 | 1 |
| route_unsnapped | 0 | 0 | 0 |
| route_second_search_refused | 0 | 0 | 200 |
| helix_name_moved | 2 | 2 | 99 |
| na_tail_straight | 1 | 2 | 12 |
| na_gap_unclear | 1 | 1 | 6 |
| na_hairpin_unclear | 0 | 0 | 0 |
| na_run_unclear | 2 | 2 | 13 |
| na_trace_unclear | 0 | 0 | 0 |
| stem_name_least_bad | 0 | 0 | 0 |
| domain_name_least_bad | 0 | 0 | 0 |
| glyph_width_guessed | 4 | 248 | 1740 |
| mark_least_bad_spot | 0 | 0 | 0 |
| counts_hidden | 6 | 9 | 83 |
| interactions_dropped | 0 | 0 | 0 |
| interactions_hidden | 0 | 0 | 0 |
| crosslinks_unplaced | 0 | 0 | 0 |
| page_text_small | 3 | 3 | 31 |
| s3_cycles_capped | 0 | 0 | 0 |

| id | default view |
|---|---|
| 8VU3 | counts_hidden 1 |
| 7YF0 | counts_hidden 1 |
| 7ZDM | counts_hidden 1 |
| 9QMS | counts_hidden 3; glyph_width_guessed 28; page_text_small 1 |
| 8U8V | glyph_width_guessed 4 |
| 7SA4 | counts_hidden 1; glyph_width_guessed 88; helix_name_moved 1; na_run_unclear 1; page_text_small 1 |
| 31FZ | counts_hidden 2; glyph_width_guessed 128; helix_name_moved 1; na_gap_unclear 1; na_run_unclear 1; na_tail_straight 2; page_text_small 1 |

## Science checks (S1-S6), apart from A1-A20

Judged in every drawn view (default, deposited, beta, focus, domains, membrane, app); S5 once per entry.

Check kinds: S1, S2, S3, S4, S6 structural consistency (the figure against its own structure, internal); S5 external validation (against CATH); A1-A20 readability (the page).

A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet drawn, no barrel, no residue mark, no accessory strand, no CATH chopping) is not judged, one with no protein does not apply.

| check | pass | fail | not judged | does not apply | kind |
|---|---|---|---|---|---|
| S1 sheet neighbours pair | 87 | 0 | 11 | 3 | consistency |
| S2 barrel closure pairs | 13 | 0 | 85 | 3 | consistency |
| S3 barrel strand count is its pairing cycle | 86 | 1 | 11 | 3 | consistency |
| S4 residue marks at their residue | 0 | 0 | 98 | 3 | consistency |
| S5 automatic domains agree with a reference | 16 | 16 | 66 | 3 | external |
| S6 accessory drawn at its partner | 28 | 0 | 70 | 3 | consistency |

### Science by view

Entries measured in each view, and how many of them fail each check there.

| view | measured | S1 | S2 | S3 | S4 | S6 |
|---|---|---|---|---|---|---|
| default | 100 | 0 | 0 | 0 | 0 | 0 |
| deposited | 97 | 0 | 0 | 0 | 0 | 0 |
| beta | 97 | 0 | 0 | 0 | 0 | 0 |
| focus | 68 | 0 | 0 | 1 | 0 | 0 |
| domains | 59 | 0 | 0 | 0 | 0 | 0 |
| membrane | 8 | 0 | 0 | 0 | 0 | 0 |
| app | 99 | 0 | 0 | 0 | 0 | 0 |

| id | S1 | S2 | S3 | S4 | S5 | S6 | findings |
|---|---|---|---|---|---|---|---|
| 4DEH | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7KUV | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8EYE | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 1JMH | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.56, boundary F1 0.00, recall 1.00 (1 reference domains, 3 automatic) |
| 7VOD | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9BNW | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 2ATZ | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2CXH | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 9AZ7 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3MPL | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6NJW | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2EY5 | ✓ | ✓ | ✓ | ○ | ✓ | ✓ |  |
| 6HJZ | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6NSY | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.86, boundary F1 0.67, recall 0.50 (3 reference domains, 2 automatic) |
| 3ZC1 | ○ | ○ | ○ | ○ | ✗ | ○ | S5 A: overlap 0.51, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); B: overlap 0.52, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); C: overlap 0.52, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); D: overlap 0.52, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); E: overlap 0.51, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.52, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); G: overlap 0.52, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); H: overlap 0.52, |
| 3SXN | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6QQQ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5QQQ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4MIN | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.73, boundary F1 0.67, recall 0.50 (2 reference domains, 2 automatic); B: overlap 0.73, boundary F1 0.67, recall 0.50 (2 reference domains, 2 automatic); C: overlap 0.73, boundary F1 0.67, recall 0.50 (2 reference domains, 2 automatic); D: overlap 0.73, boundary F1 0.67, recall 0.50 (2 reference domains, 2 automatic) |
| 1IVX | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 6BWQ | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 4OXB | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5ZJU | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 3ANZ | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); B: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); C: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); D: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); E: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); G: overlap 0.85, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); H: overlap 0.85, |
| 2FYN | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.70, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); C: overlap 0.82, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); D: overlap 0.70, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.82, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); G: overlap 0.70, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); H: overlap 0.86, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); I: overlap 0.82, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); J: overlap 0.70, |
| 4G4S | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.64, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); E: overlap 0.67, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.67, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); G: overlap 0.65, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); L: overlap 0.56, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 5S58 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3BLV | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.54, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); B: overlap 0.55, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); C: overlap 0.56, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); D: overlap 0.56, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); E: overlap 0.57, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.55, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); G: overlap 0.58, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); H: overlap 0.55, |
| 6T4C | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 5TDO | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 B: overlap 0.98, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); D: overlap 0.98, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 1TZY | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 B: overlap 0.60, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.59, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 4LZF | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.54, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 6BYU | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 5CX1 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 B: overlap 0.91, boundary F1 0.67, recall 0.75 (4 reference domains, 4 automatic); D: overlap 0.91, boundary F1 0.67, recall 0.75 (4 reference domains, 4 automatic); F: overlap 0.91, boundary F1 0.67, recall 0.75 (4 reference domains, 4 automatic); J: overlap 0.91, boundary F1 0.67, recall 0.75 (4 reference domains, 4 automatic); L: overlap 0.91, boundary F1 0.67, recall 0.75 (4 reference domains, 4 automatic) |
| 8VU3 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7TCP | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9Y9L | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6O20 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9CV6 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9VK1 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8KE9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8C5C | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8FTF | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 8UXG | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8WLE | ✓ | ✓ | ✗ | ○ | ○ | ✓ | focus: S3 sheet 86: 12 barrel strands drawn, the pairing cycle has 17 |
| 7YF0 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8OPJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8VJJ | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8FM9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2MF1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2MSV | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 2KZ2 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 2EGQ | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2LLX | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2JZI | ○ | ○ | ○ | ○ | ✗ | ○ | S5 A: overlap 0.52, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 2FJ3 | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 1FVL | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6R0Z | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8TPL | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8JD0 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8IAK | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5EQI | ○ | ○ | ○ | ○ | ✗ | ○ | S5 A: overlap 0.55, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 6C1E | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 5V5Z | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7ZDM | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8W2Q | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9QNC | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7PSX | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 6PL8 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9CJ9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6TBZ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9QMS | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 9ZP3 | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8U8V | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9NL4 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8SSW | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6JVX | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7SA4 | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 31FZ | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 9ZHB | - | - | - | - | - | - |  |
| 6WLQ | - | - | - | - | - | - |  |
| 7EDN | - | - | - | - | - | - |  |
| 6D8D | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7EB3 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9N3U | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6ZQ1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8JX8 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9JTM | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 9J7U | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8V9G | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6TS9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9W9N | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2H4E | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 1GPF | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.65, boundary F1 0.67, recall 0.67 (3 reference domains, 3 automatic); B: overlap 0.65, boundary F1 0.67, recall 0.67 (3 reference domains, 3 automatic) |
| 1W1P | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.65, boundary F1 0.67, recall 0.67 (3 reference domains, 3 automatic); B: overlap 0.65, boundary F1 0.67, recall 0.67 (3 reference domains, 3 automatic) |
| 9S44 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5DY5 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ8I574F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ6PHI8F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ9UZ20F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2CSE | ○ | ○ | ○ | ○ | ○ | ○ |  |

## S5: automatic domains against CATH (reported, not gated)

Source: CATH domain boundaries (https://www.cathdb.info) through the PDBe API (EMBL-EBI), licensed under CC BY 4.0; Sillitoe I. et al. (2021) CATH: increased structural coverage of functional space. Nucleic Acids Res 49:D266-D273.

215 chains scored in 32 entries; 59 of them with at least two CATH domains. Every chain CATH chops is scored, one-domain chains too (a one-domain chain keeps boundary F1 1 only when the automatic domains leave it whole). 0 of the chains are in the calibration set (source "reference": in-sample, tests/data/domains/cath_reference.json, the 33 chains the automatic domains were calibrated on); the never-seen check is the fresh set (cath_fresh.json).

| id | chain | source | CATH domains | automatic | overlap | boundary F1 | recall | precision |
|---|---|---|---|---|---|---|---|---|
| 4DEH | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 1JMH | A | CATH | 1 | 3 | 0.56 | 0.00 | 1.00 | 0.00 |
| 2ATZ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2CXH | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3MPL | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 2EY5 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6NSY | A | CATH | 3 | 2 | 0.86 | 0.67 | 0.50 | 1.00 |
| 6NSY | B | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 6NSY | C | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 6NSY | D | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 3ZC1 | A | CATH | 1 | 2 | 0.51 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | B | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | C | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | D | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | E | CATH | 1 | 2 | 0.51 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | F | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | G | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 3ZC1 | H | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 3SXN | A | CATH | 3 | 3 | 0.96 | 1.00 | 1.00 | 1.00 |
| 3SXN | B | CATH | 3 | 3 | 0.96 | 1.00 | 1.00 | 1.00 |
| 3SXN | C | CATH | 3 | 3 | 0.93 | 0.80 | 1.00 | 0.67 |
| 3SXN | D | CATH | 3 | 3 | 0.93 | 0.80 | 1.00 | 0.67 |
| 3SXN | E | CATH | 3 | 3 | 0.93 | 0.80 | 1.00 | 0.67 |
| 3SXN | F | CATH | 3 | 3 | 0.93 | 0.80 | 1.00 | 0.67 |
| 4MIN | A | CATH | 2 | 2 | 0.73 | 0.67 | 0.50 | 1.00 |
| 4MIN | B | CATH | 2 | 2 | 0.73 | 0.67 | 0.50 | 1.00 |
| 4MIN | C | CATH | 2 | 2 | 0.73 | 0.67 | 0.50 | 1.00 |
| 4MIN | D | CATH | 2 | 2 | 0.73 | 0.67 | 0.50 | 1.00 |
| 1IVX | A | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 1IVX | B | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 4OXB | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4OXB | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 1 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 2 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 3 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 4 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 5 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 6 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 7 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 8 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | 9 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | G | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | I | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | J | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | L | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | M | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | N | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | O | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | P | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | Q | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | R | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | S | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | T | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | U | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | V | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | W | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | X | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | Y | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | Z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | a | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | b | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | c | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | d | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | e | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | f | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | g | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | h | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | i | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | j | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | k | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | l | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | m | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | n | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | o | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | p | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | q | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | r | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | s | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | t | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | u | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | v | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | w | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | x | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5ZJU | y | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3ANZ | A | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | B | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | C | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | D | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | E | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | F | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | G | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | H | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | I | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | J | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | K | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | L | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | M | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | N | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | O | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | P | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | Q | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | R | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | S | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | T | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | U | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | V | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | W | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | X | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | Y | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | Z | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | a | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 3ANZ | b | CATH | 1 | 2 | 0.85 | 0.00 | 1.00 | 0.00 |
| 2FYN | A | CATH | 1 | 2 | 0.70 | 0.00 | 1.00 | 0.00 |
| 2FYN | B | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2FYN | C | CATH | 2 | 1 | 0.82 | 0.00 | 0.00 | 0.00 |
| 2FYN | D | CATH | 1 | 2 | 0.70 | 0.00 | 1.00 | 0.00 |
| 2FYN | E | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2FYN | F | CATH | 2 | 1 | 0.82 | 0.00 | 0.00 | 0.00 |
| 2FYN | G | CATH | 1 | 2 | 0.70 | 0.00 | 1.00 | 0.00 |
| 2FYN | H | CATH | 2 | 1 | 0.86 | 0.00 | 0.00 | 0.00 |
| 2FYN | I | CATH | 2 | 1 | 0.82 | 0.00 | 0.00 | 0.00 |
| 2FYN | J | CATH | 1 | 2 | 0.70 | 0.00 | 1.00 | 0.00 |
| 2FYN | K | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2FYN | L | CATH | 2 | 1 | 0.82 | 0.00 | 0.00 | 0.00 |
| 2FYN | M | CATH | 1 | 2 | 0.70 | 0.00 | 1.00 | 0.00 |
| 2FYN | N | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2FYN | O | CATH | 2 | 1 | 0.82 | 0.00 | 0.00 | 0.00 |
| 2FYN | P | CATH | 1 | 2 | 0.70 | 0.00 | 1.00 | 0.00 |
| 2FYN | Q | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2FYN | R | CATH | 2 | 1 | 0.82 | 0.00 | 0.00 | 0.00 |
| 4G4S | A | CATH | 1 | 2 | 0.64 | 0.00 | 1.00 | 0.00 |
| 4G4S | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | E | CATH | 1 | 2 | 0.67 | 0.00 | 1.00 | 0.00 |
| 4G4S | F | CATH | 1 | 2 | 0.67 | 0.00 | 1.00 | 0.00 |
| 4G4S | G | CATH | 1 | 2 | 0.65 | 0.00 | 1.00 | 0.00 |
| 4G4S | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | I | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | J | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | L | CATH | 1 | 2 | 0.56 | 0.00 | 1.00 | 0.00 |
| 4G4S | M | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | N | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | O | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4G4S | P | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3BLV | A | CATH | 1 | 2 | 0.54 | 0.00 | 1.00 | 0.00 |
| 3BLV | B | CATH | 1 | 2 | 0.55 | 0.00 | 1.00 | 0.00 |
| 3BLV | C | CATH | 1 | 2 | 0.56 | 0.00 | 1.00 | 0.00 |
| 3BLV | D | CATH | 1 | 2 | 0.56 | 0.00 | 1.00 | 0.00 |
| 3BLV | E | CATH | 1 | 2 | 0.57 | 0.00 | 1.00 | 0.00 |
| 3BLV | F | CATH | 1 | 2 | 0.55 | 0.00 | 1.00 | 0.00 |
| 3BLV | G | CATH | 1 | 2 | 0.58 | 0.00 | 1.00 | 0.00 |
| 3BLV | H | CATH | 1 | 2 | 0.55 | 0.00 | 1.00 | 0.00 |
| 5TDO | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5TDO | B | CATH | 1 | 2 | 0.98 | 0.00 | 1.00 | 0.00 |
| 5TDO | C | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5TDO | D | CATH | 1 | 2 | 0.98 | 0.00 | 1.00 | 0.00 |
| 1TZY | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TZY | B | CATH | 1 | 2 | 0.60 | 0.00 | 1.00 | 0.00 |
| 1TZY | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TZY | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TZY | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TZY | F | CATH | 1 | 2 | 0.59 | 0.00 | 1.00 | 0.00 |
| 1TZY | G | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TZY | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4LZF | A | CATH | 1 | 2 | 0.54 | 0.00 | 1.00 | 0.00 |
| 6BYU | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6BYU | B | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6BYU | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6BYU | G | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6BYU | H | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6BYU | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5CX1 | A | CATH | 3 | 3 | 0.95 | 1.00 | 1.00 | 1.00 |
| 5CX1 | B | CATH | 4 | 4 | 0.91 | 0.67 | 0.75 | 0.60 |
| 5CX1 | C | CATH | 3 | 3 | 0.95 | 1.00 | 1.00 | 1.00 |
| 5CX1 | D | CATH | 4 | 4 | 0.91 | 0.67 | 0.75 | 0.60 |
| 5CX1 | E | CATH | 3 | 3 | 0.95 | 1.00 | 1.00 | 1.00 |
| 5CX1 | F | CATH | 4 | 4 | 0.91 | 0.67 | 0.75 | 0.60 |
| 5CX1 | G | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 5CX1 | H | CATH | 4 | 4 | 0.98 | 0.89 | 1.00 | 0.80 |
| 5CX1 | I | CATH | 3 | 3 | 0.97 | 1.00 | 1.00 | 1.00 |
| 5CX1 | J | CATH | 4 | 4 | 0.91 | 0.67 | 0.75 | 0.60 |
| 5CX1 | K | CATH | 3 | 3 | 0.95 | 1.00 | 1.00 | 1.00 |
| 5CX1 | L | CATH | 4 | 4 | 0.91 | 0.67 | 0.75 | 0.60 |
| 5CX1 | M | CATH | 3 | 3 | 0.95 | 1.00 | 1.00 | 1.00 |
| 5CX1 | N | CATH | 4 | 4 | 0.98 | 0.89 | 1.00 | 0.80 |
| 5CX1 | O | CATH | 3 | 3 | 0.95 | 1.00 | 1.00 | 1.00 |
| 5CX1 | P | CATH | 4 | 4 | 0.98 | 0.89 | 1.00 | 0.80 |
| 2EGQ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2LLX | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2JZI | A | CATH | 1 | 2 | 0.52 | 0.00 | 1.00 | 0.00 |
| 2FJ3 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1FVL | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5EQI | A | CATH | 1 | 2 | 0.55 | 0.00 | 1.00 | 0.00 |
| 6D8D | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6D8D | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6D8D | C | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6D8D | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6D8D | E | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6D8D | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2H4E | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2H4E | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1GPF | A | CATH | 3 | 3 | 0.65 | 0.67 | 0.67 | 0.67 |
| 1GPF | B | CATH | 3 | 3 | 0.65 | 0.67 | 0.67 | 0.67 |
| 1W1P | A | CATH | 3 | 3 | 0.65 | 0.67 | 0.67 | 0.67 |
| 1W1P | B | CATH | 3 | 3 | 0.65 | 0.67 | 0.67 | 0.67 |
| mean of 215 (all) | | | | | 0.89 | 0.65 | 0.94 | 0.65 |
| mean of 59 (two domains or more) | | | | | 0.91 | 0.78 | 0.80 | 0.79 |
| mean of 156 (one domain) | | | | | 0.89 | 0.60 | 1.00 | 0.60 |

Protein entries with no chain scored (S5 not judged), by reason:

- CATH chops none of its chains: 62
- not a PDB entry (CATH chops PDB entries): 3
- no CATH chopping: 1

## Accessories named, not drawn level beside their partner (reported, not gated)

| id | named of accessory strands, per view |
|---|---|
| 9BNW | default 0 of 2; deposited 0 of 2; beta 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2; interactive 0 of 2 |
| 2EY5 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 6HJZ | default 5 of 5; deposited 5 of 5; beta 5 of 5; focus 5 of 5; domains 5 of 5; seed1 5 of 5; roundtrip 5 of 5 |
| 1IVX | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 1; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2; interactive 0 of 2 |
| 6BWQ | default 2 of 6; deposited 2 of 6; beta 2 of 6; focus 0 of 0; domains 2 of 6; seed1 2 of 6; roundtrip 2 of 6 |
| 5ZJU | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 57 of 115; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2 |
| 3ANZ | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 7 of 7; domains 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 4G4S | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 5S58 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 3BLV | default 0 of 0; deposited 0 of 0; beta 0 of 0; focus 2 of 2; domains 0 of 0; seed1 0 of 0; roundtrip 0 of 0 |
| 6T4C | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 0 of 0; domains 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
| 6BYU | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 2 of 3; domains 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
| 9Y9L | default 4 of 11; deposited 4 of 11; beta 4 of 11; focus 0 of 1; domains 4 of 11; seed1 4 of 11; roundtrip 4 of 11 |
| 8UXG | default 4 of 7; deposited 4 of 7; beta 4 of 7; focus 16 of 28; domains 4 of 7; seed1 4 of 7; roundtrip 4 of 7 |
| 8WLE | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 27 of 37; domains 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
| 7YF0 | default 3 of 6; deposited 3 of 6; beta 3 of 6; focus 0 of 0; domains 3 of 6; seed1 3 of 6; roundtrip 3 of 6 |
| 8VJJ | default 2 of 5; deposited 2 of 5; beta 2 of 5; focus 8 of 20; domains 2 of 5; seed1 2 of 5; roundtrip 2 of 5 |
| 6R0Z | default 2 of 8; deposited 2 of 8; beta 2 of 8; focus 0 of 2; domains 2 of 8; membrane 2 of 8; seed1 2 of 8; roundtrip 2 of 8 |
| 6PL8 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 9QMS | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 9ZP3 | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 1 of 1; domains 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
| 9NL4 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 7SA4 | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 31FZ | default 0 of 6; deposited 0 of 6; beta 0 of 6; focus 0 of 6; domains 0 of 6; seed1 0 of 6; roundtrip 0 of 6 |
| 7EB3 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 2H4E | default 0 of 4; deposited 0 of 4; beta 0 of 4; focus 0 of 1; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |
| 1GPF | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 1 of 1; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2; interactive 2 of 2 |
| 1W1P | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |

## Copies' sheets that differ, left out of A15 (ruling C2b-3/R2: reported, not judged)

| id | view | tie | sheets left out | judged |
|---|---|---|---|---|
| 5ZJU | focus | 2~Z | copies' sheets differ: 2 sheets 81+82 vs Z sheet 75; 14 strands not judged | yes, 8 matched |
| 5ZJU | focus | H~e | copies' sheets differ: H sheets 23+24 vs e sheet 117; 14 strands not judged | yes, 7 matched |
| 5ZJU | focus | U~u | copies' sheets differ: U sheets 59+60 vs u sheet 162; 14 strands not judged | yes, 9 matched |
| 5ZJU | focus | b~L | copies' sheets differ: b sheets 110+111 vs L sheet 36; 14 strands not judged | yes, 8 matched |
| 5ZJU | focus | g~Z | copies' sheets differ: g sheets 123+124 vs Z sheet 75; 14 strands not judged | yes, 8 matched |

## Criteria

- **A1** no crash: exit 0, files written, no traceback, no Python warning
- **A1c** a CA-only or P-only entry exits non-zero with foldmap's one-line message
- **A2** default view time
- **A2b** the app draws it within the hosted limits
- **A3** no route fallbacks
- **A4** no loop over a label
- **A5** no element or stem boxes overlapping
- **A6** no text collisions
- **A7** text and legend inside the page
- **A8** labels unique per chain, stem names per figure
- **A9** nothing left out without a note
- **A10** secondary structure agrees with mkdssp
- **A11** N->C and sheet order
- **A12** same figure for any PYTHONHASHSEED
- **A13** a saved layout redraws the same figure
- **A14** interactive page and sequence plot written
- **A15** copies drawn alike
- **A16** domains: one chain, alike, not tiny, panels apart, names clear
- **A17** nucleic-acid metrics and stem names by their letters
- **A18** membrane: only the band set in the band (TM helices across it, slab helices at their depths, crossing strands), no element nearer than its residues or on the wrong side, no loop nearer than its residues but out of a band element's port to its face or, for a loop in the bilayer, beside its own port to its residues' depth; a loop joining the faces crosses the band only beside its own port or past the band's ends
- **A19** two states: same names in the same order
- **A20** ligand tethers local

## Limits

- view_timeout_s: 900
- default_s: 60
- default_s_many_chains: 180
- many_chains: 20
- app_s: 90
- app_memory_mb: 1024
- app_max_residues: 30000
- dssp_agree: 0.99
- copies_alike: 0.5
- domain_overlap: 0.02
- stem_name_gap: figcheck.stem_gap_tolerance(layout._INK_MARGIN em)
- tether_length: 6.5
- like_place_median: 1.0
- like_place_p90: 2.5
- like_relations: 0.95
- documented warnings: none
