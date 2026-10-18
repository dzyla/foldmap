# foldmap random-100 campaign

95 of 100 pass.

An entry passes when no criterion fails. A criterion with nothing to judge in an entry is not judged (or does not apply), never counted as a pass: see the outcomes per criterion below.

Sample: seed 20261017, release-date cutoff 2026-10-13.

| stratum | entries | pass |
|---|---|---|
| S1 | 12 | 12 |
| S2 | 12 | 12 |
| S3 | 10 | 10 |
| S4 | 15 | 14 |
| S5 | 8 | 8 |
| S6 | 8 | 8 |
| S7 | 7 | 6 |
| S8 | 7 | 4 |
| S9 | 3 | 3 |
| S10 | 5 | 5 |
| S11 | 10 | 10 |
| S12 | 3 | 3 |

## Outcomes per criterion

✓ pass, ✗ fail, ○ not judged (the criterion applies but the entry gave it nothing to judge: no tie, no tether to another chain, no domain drawn, mkdssp unmeasured ...), - does not apply.

| criterion | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| A1 no crash: exit 0, files written, no traceback, no Python warning | 100 | 0 | 0 | 0 |
| A1c a CA-only or P-only entry exits non-zero with foldmap's one-line message | 0 | 0 | 0 | 100 |
| A2 default view time | 100 | 0 | 0 | 0 |
| A2b the app draws it within the hosted limits | 100 | 0 | 0 | 0 |
| A3 no route fallbacks | 97 | 3 | 0 | 0 |
| A4 no loop over a label | 100 | 0 | 0 | 0 |
| A5 no element or stem boxes overlapping | 100 | 0 | 0 | 0 |
| A6 no text collisions | 98 | 2 | 0 | 0 |
| A7 text and legend inside the page | 100 | 0 | 0 | 0 |
| A8 labels unique per chain, stem names per figure | 99 | 0 | 1 | 0 |
| A9 nothing left out without a note | 100 | 0 | 0 | 0 |
| A10 secondary structure agrees with mkdssp | 97 | 0 | 0 | 3 |
| A11 N->C and sheet order | 96 | 0 | 1 | 3 |
| A12 same figure for any PYTHONHASHSEED | 100 | 0 | 0 | 0 |
| A13 a saved layout redraws the same figure | 100 | 0 | 0 | 0 |
| A14 interactive page and sequence plot written | 100 | 0 | 0 | 0 |
| A15 copies drawn alike | 34 | 1 | 15 | 50 |
| A16 domains: one chain, alike, not tiny, panels apart, names clear | 50 | 1 | 12 | 37 |
| A17 nucleic-acid metrics and stem names by their letters | 14 | 3 | 0 | 83 |
| A18 membrane: TM helices in the band, loops out of the slab | 8 | 0 | 0 | 92 |
| A19 two states: same names in the same order | 8 | 0 | 0 | 92 |
| A20 ligand tethers local | 12 | 0 | 53 | 35 |

### Not judged, per criterion

- **A8** (1): 9BVV (no element or stem drawn)
- **A11** (1): 9BVV (no protein element drawn)
- **A15** (15): 5S5Z (copies drawn with no symmetry tie: default C~A, D~B); 6XER (copies drawn with no symmetry tie: default C~A, D~B); 3L4G (copies drawn with no symmetry tie: default D~B; focus D~B); 6ET5 (copies drawn with no symmetry tie: default F~z, K~z, P~z, S~z, V~z, Y~z, b~z, e~z, h~z, k~z, n~z, q~z, t~z, w~z, 3~z, 6~z, G~1, N~1, Q~1, T~1, W~1, Z~1, c~1, f~); 8BCU (copies drawn with no symmetry tie: default D~C, A~C, B~C, G~H, E~C, F~C, I~H); 6V6S (copies drawn with no symmetry tie: default C~A, D~B, E~A, F~B, G~A, H~B, K~I, M~A, O~N, Q~N, R~N, S~N, T~B, Y~V, b~a, c~a, d~a, e~a, f~a, g~a, h~a, i~a, j~a, k~); 7ZDP (copies drawn with no symmetry tie: default j~X); 8WMV (copies drawn with no symmetry tie: default f~m, j~m); 6ZKT (copies drawn with no symmetry tie: default j~X); 8WH0 (copies drawn with no symmetry tie: default B~A, C~A, D~A, E~A, F~A); 6L9Z (copies drawn with no symmetry tie: default E~A, F~B, G~C, H~D, K~A, L~B, M~C, N~D, O~A, P~B, Q~C, R~D; focus E~A, F~B, G~C, H~D, K~A, L~B, M~C, N~D, O~A, P~B, Q); 2VBO (copies drawn with no symmetry tie: default B~A); 7NZ2 (copies drawn with no symmetry tie: default A2~A1, A3~A1, A4~A1, B1~A1, B2~A1, B3~A1, B4~A1, C2~C1, D1~C1, D2~C1, E2~E1, E3~E1, E4~E1, F1~E1, F2~E1, F3~E1, F4~E1); 9YER (copies drawn with no symmetry tie: default C~A); 9QK1 (copies drawn with no symmetry tie: default C~B, D~B, E~B, H~G, I~G)
- **A16** (12): 3KSI (the domains view drew no domain); 5YUH (the domains view drew no domain); 6B7J (the domains view drew no domain); 2O4P (the domains view drew no domain); 3B6N (the domains view drew no domain); 3KTK (the domains view drew no domain); 2YJA (the domains view drew no domain); 3JBE (the domains view drew no domain); 1TKW (the domains view drew no domain); 2VBO (the domains view drew no domain); 5K7Z (the domains view drew no domain); 9WB3 (the domains view drew no domain)
- **A20** (53): 1PTY (no tether to another chain; reported: ligand-tether:PTR:A323:2 7.807 (home A)); 1HMY (no tether to another chain); 6EJ3 (no tether to another chain; reported: ligand-tether:B7T:A501:0 7.443 (home A); ligand-tether:B7T:A501:1 6.648 (home A)); 5RFU (no tether to another chain); 8F03 (no tether to another chain); 5YUH (no tether to another chain); 2W6W (no tether to another chain); 3LIL (no tether to another chain); 7P46 (no tether to another chain; reported: ligand-tether:HEM+CYN:B301:1 7.452 (home B); ligand-tether:HEM+CYN:A301:2 6.984 (home A); ligand-tether:KYN:D303:0 7.308 (); 8HRP (no tether to another chain; reported: ligand-tether:NAD:A401:2 7.085 (home A); ligand-tether:G3H:C402:1 8.142 (home C); ligand-tether:NAD:D401:2 8.264 (home D)); 6SNU (no tether to another chain; reported: ligand-tether:PLP:A502:1 7.162 (home A); ligand-tether:PLP:A502:3 6.731 (home A); ligand-tether:PLP:B504:1 7.148 (home B);); 3AMZ (no tether to another chain; reported: ligand-tether:FAD:A3005:2 10.299 (home A); ligand-tether:FES:A3001:2 19.25 (home A); ligand-tether:ZXF:A3009:2 7.767 (home); 8ZT9 (no tether to another chain); 7X8H (no tether to another chain); 2NAD (no tether to another chain); 2O4P (no tether to another chain; reported: ligand-tether:TPV:A300:0 7.055 (home A)); 3B6N (no tether to another chain); 5EUJ (no tether to another chain; reported: ligand-tether:TPP:A602:2 6.904 (home A); ligand-tether:TPP:A602:3 6.548 (home A)); 9SLB (no tether to another chain); 6HW8 (no tether to another chain; reported: ligand-tether:MG:G301:2 9.455 (home G); ligand-tether:GT8:H301:0 7.81 (home H); ligand-tether:GT8:V301:0 6.703 (home V); l); 6XER (no tether to another chain; reported: ligand-tether:GTP:A501:1 9.82 (home A); ligand-tether:GTP:A501:2 6.856 (home A); ligand-tether:GDP:B501:3 7.63 (home B); l); 3ULQ (no tether to another chain); 2UZU (no tether to another chain); 6N9A (no tether to another chain; reported: ligand-tether:ZN:D1402:2 8.962 (home D1); ligand-tether:KG4:D1401:3 8.415 (home D1); ligand-tether:ZN:D2402:2 8.983 (home ); 2YJA (no tether to another chain; reported: ligand-tether:EST:B11550:0 7.343 (home B1); ligand-tether:EST:B11550:3 9.438 (home B1); ligand-tether:EST:B21550:0 7.343 (); 4N9F (no tether to another chain); 3L4G (no tether to another chain); 8D3D (no tether to another chain; reported: ligand-tether:CA:A2001:4 8.03 (home A)); 6V6S (no tether to another chain; reported: ligand-tether:ADP:U800:0 15.067 (home U); ligand-tether:ADP:U800:1 25.094 (home U); ligand-tether:ADP:U800:2 12.144 (home ); 12BM (no tether to another chain; reported: ligand-tether:ADP:A501:2 9.217 (home A); ligand-tether:ADP:A501:3 7.765 (home A); ligand-tether:AGS:C501:1 18.867 (home C)); 7BKC (no tether to another chain; reported: ligand-tether:FAD:A702:3 13.647 (home A); ligand-tether:FAD:a702:0 7.359 (home a); ligand-tether:FAD:a702:3 11.157 (home a); 9Y05 (no tether to another chain; reported: ligand-tether:JDP:A901:2 6.657 (home A)); 1TKW (no tether to another chain); 6MJ2 (no tether to another chain; reported: ligand-tether:CA:A1601:2 6.912 (home A)); 6QFA (no tether to another chain); 8JWF (no tether to another chain; reported: ligand-tether:YMZ:A1500:1 6.779 (home A); ligand-tether:YMZ:A1500:2 8.191 (home A)); 8WH0 (no tether to another chain; reported: ligand-tether:ANP:A802:3 8.21 (home A); ligand-tether:ANP:C802:1 7.546 (home C); ligand-tether:ANP:F802:1 6.752 (home F); ); 4NLG (no tether to another chain); 2FD8 (no tether to another chain); 13FZ (no tether to another chain; reported: ligand-tether:AGS:A2000:1 6.788 (home A); ligand-tether:AGS:A2000:2 13.886 (home A); ligand-tether:AGS:A2000:3 8.4 (home A); 4E8Q (no tether to another chain); 9IWF (no tether to another chain); 9YER (no tether to another chain); 9QK1 (no tether to another chain); 3IB0 (no tether to another chain; reported: ligand-tether:CO3:A687:3 6.71 (home A)); 6A4Y (no tether to another chain; reported: ligand-tether:HEM:A601:0 8.952 (home A)); 8JXD (no tether to another chain); 1C6A (no tether to another chain); 244L (no tether to another chain); 8U9T (no tether to another chain; reported: ligand-tether:W1L:A1401:3 7.845 (home A1); ligand-tether:W1L:A2401:3 7.884 (home A2)); 3U3H (no tether to another chain; reported: ligand-tether:03W+FMT:A1400:1 10.047 (home A1); ligand-tether:MRD:A12000:2 7.095 (home A1); ligand-tether:MRD:A32000:2 7.0); 6QND (no tether to another chain; reported: ligand-tether:MG:A402:0 8.701 (home A); ligand-tether:MG:A402:2 10.774 (home A); ligand-tether:MG:A402:3 7.584 (home A); l); 2BBH (no tether to another chain; reported: ligand-tether:DMU:A402:2 6.502 (home A); ligand-tether:DMU:A402:3 9.734 (home A))

## Entries

process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 gates; drawing s: the drawing inside it.

| id | stratum | bin | chains | process s | drawing s | A1 | A1c | A2 | A2b | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 | A11 | A12 | A13 | A14 | A15 | A16 | A17 | A18 | A19 | A20 | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1PTY | S1 | M | 1 | 0.9 | 0.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6TX2 | S1 | M | 1 | 0.85 | 0.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 7N71 | S1 | M | 1 | 1.25 | 0.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 1HMY | S1 | M | 1 | 0.9 | 0.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6EJ3 | S1 | M | 1 | 1.25 | 0.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5RFU | S1 | M | 1 | 1.05 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2FR2 | S1 | S | 1 | 0.7 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8F03 | S1 | S | 1 | 1.1 | 0.48 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3KSI | S1 | S | 1 | 0.95 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5YUH | S1 | S | 1 | 0.85 | 0.23 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 2W6W | S1 | S | 1 | 0.75 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3LIL | S1 | S | 1 | 0.8 | 0.18 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 7P46 | S2 | L | 8 | 1.7 | 0.9 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8HRP | S2 | L | 8 | 2.05 | 1.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6SNU | S2 | L | 4 | 1.7 | 0.96 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 3AMZ | S2 | L | 2 | 3.46 | 2.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8ZT9 | S2 | M | 2 | 1.35 | 0.71 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7X8H | S2 | M | 1 | 1.5 | 0.79 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2NAD | S2 | M | 2 | 1.31 | 0.6 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6B7J | S2 | M | 2 | 0.85 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 2O4P | S2 | S | 2 | 0.8 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 3B6N | S2 | S | 1 | 0.9 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 5EUJ | S2 | XL | 24 | 1.65 | 0.79 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9SLB | S2 | XL | 10 | 2.1 | 1.08 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ○ |  |
| 3KTK | S3 | L | 28 | 1.2 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5S5Z | S3 | L | 6 | 3.21 | 2.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 6HW8 | S3 | L | 28 | 7.11 | 5.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6XER | S3 | L | 5 | 2.61 | 1.59 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 3ULQ | S3 | M | 2 | 1.15 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 2UZU | S3 | M | 2 | 0.9 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 6N9A | S3 | M | 3 | 1.9 | 1.04 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2YJA | S3 | S | 2 | 1.05 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 4N9F | S3 | XL | 60 | 1.6 | 0.76 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 3L4G | S3 | XL | 16 | 2.4 | 1.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 8FLU | S4 | L | 6 | 2.2 | 1.54 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 6ET5 | S4 | L | 54 | 5.56 | 3.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 7KH0 | S4 | L | 7 | 3.21 | 2.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 8BCU | S4 | L | 9 | 3.51 | 2.56 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | - |  |
| 7ZZ1 | S4 | L | 2 | 1.1 | 0.43 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 3JBE | S4 | M | 5 | 1.65 | 0.96 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | - |  |
| 6WUT | S4 | M | 3 | 1.3 | 0.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 8D3D | S4 | XL | 16 | 17.28 | 15.48 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6V6S | S4 | XL | 38 | 24.84 | 21.0 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 12BM | S4 | XL | 19 | 11.57 | 9.7 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 9GTP | S4 | XL | 60 | 4.16 | 2.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | - | - | A15: focus 41.39338891387734 |
| 7BKC | S4 | XL | 26 | 13.02 | 9.25 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7SPB | S4 | XL | 78 | 1.6 | 0.75 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 7ZDP | S4 | XL | 45 | 12.52 | 9.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 9Y05 | S4 | XL | 12 | 1.7 | 0.85 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 1TKW | S5 | M | 2 | 1.0 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 2MRJ | S5 | S | 1 | 0.7 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 5GIW | S5 | S | 1 | 0.6 | 0.07 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2LZJ | S5 | S | 1 | 0.7 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2N5J | S5 | S | 1 | 0.65 | 0.1 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8SKX | S5 | S | 1 | 0.6 | 0.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 1RCH | S5 | S | 1 | 0.7 | 0.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9BVV | S5 | S | 1 | 0.55 | 0.03 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | ○ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9CZH | S6 | L | 8 | 2.91 | 1.92 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ |  |
| 6MJ2 | S6 | L | 4 | 1.95 | 1.05 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ○ |  |
| 6QFA | S6 | L | 10 | 1.65 | 0.93 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ○ |  |
| 8WMV | S6 | L | 28 | 22.24 | 15.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✓ | - | ✓ |  |
| 8JWF | S6 | M | 1 | 1.65 | 0.9 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 7JI3 | S6 | M | 3 | 1.2 | 0.58 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - |  |
| 7VVK | S6 | M | 6 | 2.25 | 1.57 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | - |  |
| 6ZKT | S6 | XL | 45 | 16.38 | 13.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✓ | - | ✓ |  |
| 8WH0 | S7 | L | 7 | 5.66 | 2.34 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | - | - | ○ |  |
| 6L9Z | S7 | L | 19 | 6.81 | 3.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | - | - | ✓ |  |
| 4NLG | S7 | M | 3 | 2.3 | 0.58 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 2VBO | S7 | M | 4 | 2.15 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | ✓ | - | - | ✓ |  |
| 5K7Z | S7 | M | 6 | 3.06 | 0.84 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | ✓ |  |
| 2FD8 | S7 | S | 2 | 1.85 | 0.25 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 7NZ2 | S7 | XL | 44 | 57.55 | 43.7 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✗ | ✗ | - | - | ✓ | A3: domains: B3:211-221>B3:227-408, B3:1051-1123>B3:1136-1140, F3:182-186>F3:197-205; A6: domains 3 (domain-label:D7 · A1 × loop:B3:211-221>B3:227-408; domain-label:D1 · A3 × loop:F3:182-186>F3:197-205; domain-label:D7 · B3 × loop:B3:211-221>B3:227-408); A16: domains 3 (domain-label:D7 · A1 × loop:B3:211-221>B3:227-408; domain-label:D1 · A3 × loop:F3:182-186>F3:197-205; domain-label:D7 · B3 × loop:B3:211-221>B3:227-408); A17: spearman 0.446 (limit 0.6) |
| 9K4F | S8 | L | 4 | 3.46 | 1.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | - | A17: stems_nearer_frac 0.875 (limit 1.0) |
| 13FZ | S8 | L | 4 | 8.06 | 5.55 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ |  |
| 9WB3 | S8 | M | 4 | 2.4 | 0.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | ✓ | - | - | - |  |
| 7F3J | S8 | M | 8 | 1.7 | 0.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 2XLJ | S8 | S | 2 | 1.85 | 0.22 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 6S13 | S8 | XL | 52 | 47.43 | 32.1 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | - | A3: domains: strand A:2729-2729, strand A:2732-2733, strand A:2737-2738; A6: domains 3 (na:stem:80:first × na:run:239; na:stem:129:first × na:run:239; na:stem:129:first × na:run:240); A17: pair_f1 0.944 (limit 0.95) |
| 8ZTU | S8 | XL | 56 | 60.74 | 43.48 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | ✓ | A3: domains: strand a:2701-2702, strand a:2705-2705, strand a:2710-2711 |
| 4E8Q | S9 | M | 1 | 3.81 | 1.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 9IWF | S9 | S | 2 | 2.0 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 1Z43 | S9 | S | 1 | 1.9 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 9YER | S10 | L | 5 | 5.06 | 4.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 9QK1 | S10 | L | 9 | 7.96 | 7.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 3IB0 | S10 | M | 1 | 2.86 | 2.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6A4Y | S10 | M | 1 | 2.1 | 1.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 8JXD | S10 | XL | 2 | 4.16 | 3.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 4P39 | S11 | M | 4 | 0.65 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | - |  |
| 3HQB | S11 | S | 2 | 0.8 | 0.23 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | - |  |
| 1C6A | S11 | S | 1 | 0.75 | 0.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 244L | S11 | S | 1 | 0.7 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 8U9T | S11 | M | 1 | 1.1 | 0.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7IT6 | S11 | S | 1 | 0.75 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 3U3H | S11 | M | 1 | 2.05 | 1.18 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ○ |  |
| 6QND | S11 | M | 1 | 0.95 | 0.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 2BBJ | S11 | L | 5 | 2.2 | 1.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - |  |
| 2BBH | S11 | S | 1 | 0.85 | 0.25 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| AF_AFQ387E4F1 | S12 | M | 1 | 0.8 | 0.25 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFQ8W0Z7F1 | S12 | M | 1 | 1.05 | 0.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFC7J156F1 | S12 | S | 1 | 0.65 | 0.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |

## Two-state pairs (S11)

| pair | seed | partner | sampled as | types now | ligand touching | assembly 1 chains | CA RMSD Å | A19 |
|---|---|---|---|---|---|---|---|---|
| P1 | 3HQB | 4P39 | open/closed | open/closed | [False, False] | [1, 1] | 4.297 | ✓ names compared on all 3 pairs; A:3-11 unpaired |
| P2 | 1C6A | 244L | apo/holo | apo/holo | [False, True] | [1, 1] | 0.219 | ✓ names compared on all 13 pairs |
| P3 | 8U9T | 7IT6 | apo/holo | not one protein (50% identity over 113 aligned residues) | None | None | None | - not one protein: 50% identity over 113 aligned residues (306 and 113 residues) |
| P4 | 3U3H | 6QND | monomer/oligomer | apo/holo, monomer/oligomer | [True, False] | [4, 1] | 0.439 | ✓ names compared on all 30 pairs; A1:146-148 unpaired |
| P5 | 2BBJ | 2BBH | monomer/oligomer | monomer/oligomer, open/closed | [False, False] | [5, 1] | 3.962 | ✓ names compared on all 13 pairs; A:193-201 unpaired; A:244-273 unpaired; A:275-299 unpaired; A:302-310 unpaired; A:330-344 unpaired; A:47-49 unpaired; A:161-163 unpaired |

## Reported, not gated

| id | honest_crossings | loop_crossings | na_run_crossings | fallback_curves | fallback_curves_crossing_an_element | kink_links_rejected | tighter_than_inflate | rejected_as_crossing | names_moved | area | max_rss_mb |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1PTY | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 282.83 | 162.1 |
| 6TX2 | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 358.13 | 160.6 |
| 7N71 | 29 | 29 | None | 0 | 0 | 0 | 0 | 0 | 0 | 952.48 | 193.8 |
| 1HMY | 14 | 14 | None | 0 | 0 | 0 | 0 | 0 | 0 | 333.48 | 162.8 |
| 6EJ3 | 26 | 26 | None | 0 | 0 | 0 | 0 | 0 | 0 | 660.62 | 164.6 |
| 5RFU | 16 | 16 | None | 0 | 0 | 0 | 0 | 0 | 0 | 654.74 | 175.0 |
| 2FR2 | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 158.5 | 159.7 |
| 8F03 | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 194.5 | 161.5 |
| 3KSI | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 322.04 | 162.6 |
| 5YUH | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 152.61 | 161.4 |
| 2W6W | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 178.43 | 160.3 |
| 3LIL | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 106.9 | 160.7 |
| 7P46 | 19 | 19 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1187.83 | 209.6 |
| 8HRP | 143 | 143 | None | 0 | 0 | 0 | 2 | 0 | 0 | 2037.26 | 231.2 |
| 6SNU | 100 | 100 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1551.74 | 320.5 |
| 3AMZ | 257 | 257 | None | 0 | 0 | 0 | 1 | 0 | 0 | 3464.26 | 344.4 |
| 8ZT9 | 19 | 19 | None | 0 | 0 | 0 | 0 | 0 | 0 | 725.24 | 175.7 |
| 7X8H | 52 | 52 | None | 0 | 0 | 0 | 2 | 0 | 0 | 875.38 | 194.8 |
| 2NAD | 51 | 51 | None | 0 | 0 | 0 | 0 | 0 | 0 | 832.99 | 181.1 |
| 6B7J | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 524.75 | 161.7 |
| 2O4P | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 381.13 | 161.5 |
| 3B6N | 12 | 12 | None | 0 | 0 | 0 | 1 | 0 | 0 | 616.9 | 163.8 |
| 5EUJ | 20 | 20 | None | 0 | 0 | 0 | 1 | 0 | 0 | 729.18 | 336.6 |
| 9SLB | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1270.52 | 411.7 |
| 3KTK | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 285.14 | 227.4 |
| 5S5Z | 159 | 159 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4315.14 | 359.2 |
| 6HW8 | 358 | 358 | None | 0 | 0 | 0 | 2 | 0 | 0 | 8425.19 | 475.2 |
| 6XER | 128 | 128 | None | 0 | 0 | 0 | 1 | 0 | 0 | 3598.61 | 346.2 |
| 3ULQ | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 646.25 | 164.5 |
| 2UZU | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 346.75 | 162.2 |
| 6N9A | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1652.02 | 222.6 |
| 2YJA | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 664.18 | 165.5 |
| 4N9F | 15 | 15 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1238.42 | 308.6 |
| 3L4G | 144 | 144 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2497.48 | 311.6 |
| 8FLU | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1666.39 | 218.4 |
| 6ET5 | 37 | 37 | None | 0 | 0 | 0 | 1 | 0 | 1 | 25390.44 | 331.2 |
| 7KH0 | 48 | 48 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1786.7 | 227.8 |
| 8BCU | 307 | 307 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3672.39 | 391.1 |
| 7ZZ1 | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 685.99 | 170.7 |
| 3JBE | 64 | 64 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1076.66 | 173.4 |
| 6WUT | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1058.6 | 212.5 |
| 8D3D | 45 | 45 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1628.23 | 554.1 |
| 6V6S | 1237 | 1237 | None | 0 | 0 | 0 | 5 | 0 | 0 | 21320.9 | 619.7 |
| 12BM | 405 | 405 | None | 0 | 0 | 0 | 2 | 0 | 0 | 10598.39 | 573.3 |
| 9GTP | 138 | 138 | None | 0 | 0 | 0 | 2 | 0 | 0 | 4545.1 | 481.9 |
| 7BKC | 620 | 620 | None | 0 | 0 | 0 | 7 | 0 | 0 | 13871.11 | 1456.9 |
| 7SPB | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 420.27 | 321.1 |
| 7ZDP | 293 | 293 | None | 0 | 0 | 0 | 3 | 0 | 0 | 15620.27 | 514.5 |
| 9Y05 | 26 | 26 | None | 0 | 0 | 0 | 1 | 0 | 0 | 947.57 | 358.5 |
| 1TKW | 21 | 21 | None | 0 | 0 | 0 | 0 | 0 | 0 | 336.08 | 248.7 |
| 2MRJ | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 67.93 | 189.8 |
| 5GIW | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 8.58 | 159.2 |
| 2LZJ | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 65.53 | 194.1 |
| 2N5J | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 45.29 | 161.4 |
| 8SKX | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 54.42 | 165.8 |
| 1RCH | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 136.41 | 166.6 |
| 9BVV | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 157.4 |
| 9CZH | 42 | 42 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1787.96 | 304.4 |
| 6MJ2 | 60 | 60 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1561.16 | 267.0 |
| 6QFA | 18 | 18 | None | 0 | 0 | 0 | 0 | 0 | 0 | 643.24 | 223.6 |
| 8WMV | 228 | 228 | None | 0 | 0 | 0 | 0 | 0 | 0 | 11009.35 | 444.4 |
| 8JWF | 58 | 58 | None | 0 | 0 | 0 | 3 | 0 | 0 | 1625.4 | 184.5 |
| 7JI3 | 61 | 61 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1219.04 | 192.0 |
| 7VVK | 32 | 32 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1729.5 | 222.8 |
| 6ZKT | 279 | 279 | None | 0 | 0 | 0 | 1 | 0 | 0 | 16192.6 | 433.0 |
| 8WH0 | 249 | 249 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 3254.26 | 372.8 |
| 6L9Z | 138 | 138 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 9369.61 | 373.7 |
| 4NLG | 24 | 24 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 709.03 | 173.3 |
| 2VBO | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 672.1 | 176.8 |
| 5K7Z | 11 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1226.5 | 197.9 |
| 2FD8 | 9 | 9 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 180.42 | 166.6 |
| 7NZ2 | 2240 | 2240 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 38235.07 | 1004.3 |
| 9K4F | 44 | 44 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 2245.59 | 199.2 |
| 13FZ | 235 | 235 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2425.15 | 345.5 |
| 9WB3 | 25 | 25 | 4 | 0 | 0 | 0 | 1 | 0 | 0 | 997.81 | 176.9 |
| 7F3J | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 152.4 | 169.1 |
| 2XLJ | 4 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 272.1 | 166.8 |
| 6S13 | 111 | 111 | 326 | 0 | 0 | 0 | 0 | 0 | 0 | 23431.42 | 736.3 |
| 8ZTU | 293 | 293 | 714 | 0 | 0 | 0 | 2 | 0 | 0 | 35176.31 | 632.8 |
| 4E8Q | 0 | 0 | 10 | 0 | 0 | 0 | 0 | 0 | 0 | 2288.24 | 199.0 |
| 9IWF | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 232.25 | 170.6 |
| 1Z43 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 341.76 | 167.7 |
| 9YER | 99 | 99 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2856.02 | 341.3 |
| 9QK1 | 92 | 92 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3730.69 | 311.4 |
| 3IB0 | 27 | 27 | None | 0 | 0 | 0 | 0 | 0 | 0 | 429.17 | 163.8 |
| 6A4Y | 41 | 41 | None | 0 | 0 | 0 | 1 | 0 | 0 | 558.23 | 165.9 |
| 8JXD | 59 | 59 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1079.91 | 214.5 |
| 4P39 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 63.27 | 158.9 |
| 3HQB | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 66.99 | 158.1 |
| 1C6A | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 190.57 | 160.3 |
| 244L | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 190.57 | 159.9 |
| 8U9T | 20 | 20 | None | 0 | 0 | 0 | 0 | 0 | 0 | 652.69 | 170.1 |
| 7IT6 | 10 | 10 | None | 0 | 0 | 0 | 1 | 0 | 0 | 101.47 | 159.7 |
| 3U3H | 71 | 71 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2235.22 | 255.6 |
| 6QND | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 627.23 | 164.0 |
| 2BBJ | 32 | 32 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2577.72 | 269.4 |
| 2BBH | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 331.67 | 161.6 |
| AF_AFQ387E4F1 | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 373.18 | 162.7 |
| AF_AFQ8W0Z7F1 | 32 | 32 | None | 0 | 0 | 0 | 0 | 0 | 0 | 719.99 | 164.2 |
| AF_AFC7J156F1 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 113.42 | 158.9 |

## Science checks (S1-S6), apart from A1-A20

Judged in every drawn view (default, deposited, beta, focus, domains, membrane, app); S5 once per entry.

A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet drawn, no barrel, no residue mark, no accessory strand, no CATH chopping) is not judged, one with no protein does not apply.

| check | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| S1 sheet neighbours pair | 86 | 0 | 11 | 3 |
| S2 barrel closure pairs | 10 | 0 | 87 | 3 |
| S3 barrel strand count is its pairing cycle | 86 | 0 | 11 | 3 |
| S4 residue marks at their residue | 0 | 0 | 97 | 3 |
| S5 automatic domains agree with a reference | 19 | 18 | 60 | 3 |
| S6 accessory drawn at its partner | 28 | 0 | 69 | 3 |

### Science by view

Entries measured in each view, and how many of them fail each check there.

| view | measured | S1 | S2 | S3 | S4 | S6 |
|---|---|---|---|---|---|---|
| default | 100 | 0 | 0 | 0 | 0 | 0 |
| deposited | 97 | 0 | 0 | 0 | 0 | 0 |
| beta | 97 | 0 | 0 | 0 | 0 | 0 |
| focus | 64 | 0 | 0 | 0 | 0 | 0 |
| domains | 63 | 0 | 0 | 0 | 0 | 0 |
| membrane | 8 | 0 | 0 | 0 | 0 | 0 |
| app | 100 | 0 | 0 | 0 | 0 | 0 |

| id | S1 | S2 | S3 | S4 | S5 | S6 | findings |
|---|---|---|---|---|---|---|---|
| 1PTY | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6TX2 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7N71 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 1HMY | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.70, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 6EJ3 | ✓ | ✓ | ✓ | ○ | ✓ | ✓ |  |
| 5RFU | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 2FR2 | ✓ | ✓ | ✓ | ○ | ✓ | ✓ |  |
| 8F03 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3KSI | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5YUH | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2W6W | ○ | ○ | ○ | ○ | ✓ | ○ |  |
| 3LIL | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7P46 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8HRP | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6SNU | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3AMZ | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.79, boundary F1 0.92, recall 0.92 (10 reference domains, 8 automatic); B: overlap 0.79, boundary F1 0.92, recall 0.92 (10 reference domains, 8 automatic) |
| 8ZT9 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7X8H | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2NAD | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6B7J | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2O4P | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 3B6N | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5EUJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9SLB | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3KTK | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 5S5Z | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6HW8 | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 B: overlap 0.66, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); C: overlap 0.65, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); F: overlap 0.64, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); G: overlap 0.65, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); K: overlap 0.56, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); P: overlap 0.66, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); Q: overlap 0.65, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); T: overlap 0.64, |
| 6XER | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3ULQ | ○ | ○ | ○ | ○ | ✗ | ○ | S5 A: overlap 0.35, boundary F1 0.00, recall 1.00 (1 reference domains, 4 automatic) |
| 2UZU | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 E: overlap 0.85, boundary F1 0.50, recall 0.33 (2 reference domains, 2 automatic) |
| 6N9A | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 D: overlap 0.73, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic) |
| 2YJA | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4N9F | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 3: overlap 0.85, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); 9: overlap 0.85, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); C: overlap 0.81, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); I: overlap 0.85, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); O: overlap 0.85, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); U: overlap 0.85, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); V: overlap 0.85, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic); f: overlap 0.85, |
| 3L4G | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.66, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); C: overlap 0.91, boundary F1 0.89, recall 0.80 (4 reference domains, 3 automatic); E: overlap 0.82, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); I: overlap 0.82, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); K: overlap 0.82, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); O: overlap 0.82, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 8FLU | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6ET5 | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 H: overlap 0.86, boundary F1 0.00, recall 0.00 (2 reference domains, 2 automatic); L: overlap 0.75, boundary F1 0.40, recall 1.00 (2 reference domains, 3 automatic); M: overlap 0.84, boundary F1 0.00, recall 0.00 (2 reference domains, 2 automatic) |
| 7KH0 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8BCU | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 7ZZ1 | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 3JBE | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6WUT | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8D3D | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6V6S | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 12BM | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 9GTP | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7BKC | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7SPB | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7ZDP | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9Y05 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 1TKW | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 B: overlap 0.76, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 2MRJ | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5GIW | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 2LZJ | ✓ | ✓ | ✓ | ○ | ✓ | ○ |  |
| 2N5J | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8SKX | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1RCH | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 9BVV | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 9CZH | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 6MJ2 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6QFA | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8WMV | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8JWF | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7JI3 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7VVK | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6ZKT | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8WH0 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6L9Z | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4NLG | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2VBO | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5K7Z | ○ | ○ | ○ | ○ | ✗ | ○ | S5 B: overlap 0.78, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); C: overlap 0.77, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); D: overlap 0.78, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 2FD8 | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7NZ2 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9K4F | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 13FZ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9WB3 | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 7F3J | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 2XLJ | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.88, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 6S13 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8ZTU | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 4E8Q | - | - | - | - | - | - |  |
| 9IWF | - | - | - | - | - | - |  |
| 1Z43 | - | - | - | - | - | - |  |
| 9YER | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9QK1 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3IB0 | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.91, boundary F1 0.80, recall 0.67 (2 reference domains, 2 automatic) |
| 6A4Y | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 8JXD | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4P39 | ○ | ○ | ○ | ○ | ✓ | ○ |  |
| 3HQB | ○ | ○ | ○ | ○ | ✓ | ○ |  |
| 1C6A | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.51, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 244L | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.51, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 8U9T | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7IT6 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3U3H | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.83, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 6QND | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2BBJ | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.89, boundary F1 0.50, recall 0.50 (3 reference domains, 3 automatic) |
| 2BBH | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.99, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| AF_AFQ387E4F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ8W0Z7F1 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| AF_AFC7J156F1 | ○ | ○ | ○ | ○ | ○ | ○ |  |

## S5: automatic domains against CATH (reported, not gated)

Source: CATH domain boundaries (https://www.cathdb.info) through the PDBe API (EMBL-EBI), licensed under CC BY 4.0; Sillitoe I. et al. (2021) CATH: increased structural coverage of functional space. Nucleic Acids Res 49:D266-D273.

164 chains scored in 37 entries; 46 of them with at least two CATH domains. Every chain CATH chops is scored, one-domain chains too (a one-domain chain keeps boundary F1 1 only when the automatic domains leave it whole). 0 of the chains are in the calibration set (source "reference": in-sample, tests/data/domains/cath_reference.json, the 33 chains the automatic domains were calibrated on); the never-seen check is the fresh set (cath_fresh.json).

| id | chain | source | CATH domains | automatic | overlap | boundary F1 | recall | precision |
|---|---|---|---|---|---|---|---|---|
| 1PTY | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1HMY | A | CATH | 2 | 1 | 0.70 | 0.00 | 0.00 | 0.00 |
| 6EJ3 | A | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2FR2 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3KSI | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2W6W | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3LIL | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3AMZ | A | CATH | 10 | 8 | 0.79 | 0.92 | 0.92 | 0.92 |
| 3AMZ | B | CATH | 10 | 8 | 0.79 | 0.92 | 0.92 | 0.92 |
| 2NAD | A | CATH | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 2NAD | B | CATH | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 6B7J | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6B7J | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2O4P | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2O4P | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3B6N | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | B | CATH | 1 | 2 | 0.66 | 0.00 | 1.00 | 0.00 |
| 6HW8 | C | CATH | 1 | 2 | 0.65 | 0.00 | 1.00 | 0.00 |
| 6HW8 | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | F | CATH | 1 | 2 | 0.64 | 0.00 | 1.00 | 0.00 |
| 6HW8 | G | CATH | 1 | 2 | 0.65 | 0.00 | 1.00 | 0.00 |
| 6HW8 | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | I | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | J | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | K | CATH | 1 | 2 | 0.56 | 0.00 | 1.00 | 0.00 |
| 6HW8 | L | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | M | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | N | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | O | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | P | CATH | 1 | 2 | 0.66 | 0.00 | 1.00 | 0.00 |
| 6HW8 | Q | CATH | 1 | 2 | 0.65 | 0.00 | 1.00 | 0.00 |
| 6HW8 | R | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | S | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | T | CATH | 1 | 2 | 0.64 | 0.00 | 1.00 | 0.00 |
| 6HW8 | U | CATH | 1 | 2 | 0.65 | 0.00 | 1.00 | 0.00 |
| 6HW8 | V | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | W | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | X | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | Y | CATH | 1 | 2 | 0.56 | 0.00 | 1.00 | 0.00 |
| 6HW8 | Z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | a | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6HW8 | b | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3ULQ | A | CATH | 1 | 4 | 0.35 | 0.00 | 1.00 | 0.00 |
| 3ULQ | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2UZU | E | CATH | 2 | 2 | 0.85 | 0.50 | 0.33 | 1.00 |
| 6N9A | D | CATH | 2 | 3 | 0.73 | 0.67 | 1.00 | 0.50 |
| 4N9F | 0 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | 3 | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | 4 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | 5 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | 6 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | 9 | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | C | CATH | 2 | 3 | 0.81 | 0.67 | 1.00 | 0.50 |
| 4N9F | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | I | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | J | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | L | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | N | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | O | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | P | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | Q | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | R | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | T | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | U | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | V | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | W | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | X | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | Y | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | Z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | a | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | c | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | e | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | f | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | g | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | h | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | i | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | k | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | l | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | m | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | n | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | o | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | r | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | s | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | t | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | u | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | w | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | x | CATH | 2 | 3 | 0.85 | 0.67 | 1.00 | 0.50 |
| 4N9F | y | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4N9F | z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | A | CATH | 1 | 2 | 0.66 | 0.00 | 1.00 | 0.00 |
| 3L4G | B | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | C | CATH | 4 | 3 | 0.91 | 0.89 | 0.80 | 1.00 |
| 3L4G | D | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | E | CATH | 1 | 2 | 0.82 | 0.00 | 1.00 | 0.00 |
| 3L4G | F | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | G | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | H | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | I | CATH | 1 | 2 | 0.82 | 0.00 | 1.00 | 0.00 |
| 3L4G | J | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | K | CATH | 1 | 2 | 0.82 | 0.00 | 1.00 | 0.00 |
| 3L4G | L | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | M | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | N | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L4G | O | CATH | 1 | 2 | 0.82 | 0.00 | 1.00 | 0.00 |
| 3L4G | P | CATH | 4 | 4 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | 3 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | 6 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | C | CATH | 2 | 2 | 0.96 | 0.80 | 1.00 | 0.67 |
| 6ET5 | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | H | CATH | 2 | 2 | 0.86 | 0.00 | 0.00 | 0.00 |
| 6ET5 | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | L | CATH | 2 | 3 | 0.75 | 0.40 | 1.00 | 0.25 |
| 6ET5 | M | CATH | 2 | 2 | 0.84 | 0.00 | 0.00 | 0.00 |
| 6ET5 | P | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | S | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | V | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | Y | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | b | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | e | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | h | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | k | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | n | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | q | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | t | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | w | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6ET5 | z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TKW | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1TKW | B | CATH | 2 | 1 | 0.76 | 0.00 | 0.00 | 0.00 |
| 2MRJ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2LZJ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1RCH | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4NLG | A | CATH | 4 | 4 | 0.96 | 1.00 | 1.00 | 1.00 |
| 2VBO | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2VBO | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5K7Z | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5K7Z | B | CATH | 2 | 1 | 0.78 | 0.00 | 0.00 | 0.00 |
| 5K7Z | C | CATH | 2 | 1 | 0.77 | 0.00 | 0.00 | 0.00 |
| 5K7Z | D | CATH | 2 | 1 | 0.78 | 0.00 | 0.00 | 0.00 |
| 2FD8 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2XLJ | A | CATH | 2 | 1 | 0.88 | 0.00 | 0.00 | 0.00 |
| 3IB0 | A | CATH | 2 | 2 | 0.91 | 0.80 | 0.67 | 1.00 |
| 6A4Y | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4P39 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4P39 | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4P39 | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4P39 | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3HQB | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3HQB | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1C6A | A | CATH | 1 | 2 | 0.51 | 0.00 | 1.00 | 0.00 |
| 244L | A | CATH | 1 | 2 | 0.51 | 0.00 | 1.00 | 0.00 |
| 3U3H | A | CATH | 1 | 2 | 0.83 | 0.00 | 1.00 | 0.00 |
| 2BBJ | A | CATH | 3 | 3 | 0.89 | 0.50 | 0.50 | 0.50 |
| 2BBJ | B | CATH | 3 | 3 | 0.98 | 1.00 | 1.00 | 1.00 |
| 2BBJ | D | CATH | 3 | 3 | 0.99 | 1.00 | 1.00 | 1.00 |
| 2BBJ | E | CATH | 3 | 3 | 0.98 | 1.00 | 1.00 | 1.00 |
| 2BBJ | F | CATH | 3 | 3 | 0.98 | 1.00 | 1.00 | 1.00 |
| 2BBH | A | CATH | 1 | 2 | 0.99 | 0.00 | 1.00 | 0.00 |
| mean of 164 (all) | | | | | 0.93 | 0.79 | 0.94 | 0.78 |
| mean of 46 (two domains or more) | | | | | 0.89 | 0.68 | 0.79 | 0.65 |
| mean of 118 (one domain) | | | | | 0.94 | 0.83 | 1.00 | 0.83 |

Protein entries with no chain scored (S5 not judged), by reason:

- CATH chops none of its chains: 57
- not a PDB entry (CATH chops PDB entries): 3

## Accessories named, not drawn level beside their partner (reported, not gated)

| id | named of accessory strands, per view |
|---|---|
| 6EJ3 | default 2 of 7; deposited 2 of 7; beta 2 of 7; seed1 2 of 7; roundtrip 2 of 7 |
| 5RFU | default 0 of 2; deposited 0 of 2; beta 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 2FR2 | default 0 of 1; deposited 0 of 1; beta 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 8F03 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 8HRP | default 0 of 4; deposited 0 of 4; beta 0 of 4; focus 0 of 4; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |
| 8ZT9 | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 1; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 2O4P | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 1; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 9SLB | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 5 of 10; seed1 0 of 1; roundtrip 0 of 1 |
| 3KTK | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 14; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 5S5Z | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 6HW8 | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 1; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 4N9F | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 1 of 1; domains 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 6ET5 | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 1 of 1; domains 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 8BCU | default 3 of 18; deposited 3 of 18; beta 3 of 18; focus 0 of 3; domains 3 of 18; seed1 3 of 18; roundtrip 3 of 18 |
| 3JBE | default 1 of 3; deposited 1 of 3; beta 1 of 3; focus 1 of 2; domains 1 of 3; seed1 1 of 3; roundtrip 1 of 3 |
| 6WUT | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 6V6S | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 12BM | default 4 of 25; deposited 4 of 25; beta 4 of 25; focus 0 of 0; domains 4 of 25; seed1 4 of 25; roundtrip 4 of 25 |
| 9GTP | default 0 of 0; deposited 0 of 0; beta 0 of 0; focus 0 of 3; domains 0 of 0; seed1 0 of 0; roundtrip 0 of 0 |
| 6QFA | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 10; domains 0 of 2; membrane 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 7JI3 | default 1 of 4; deposited 1 of 4; beta 1 of 4; focus 1 of 2; domains 1 of 4; membrane 1 of 4; seed1 1 of 4; roundtrip 1 of 4 |
| 9K4F | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 8ZTU | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 9QK1 | default 0 of 10; deposited 0 of 10; beta 0 of 10; focus 0 of 2; domains 0 of 10; seed1 0 of 10; roundtrip 0 of 10 |
| 3IB0 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 8U9T | default 0 of 2; deposited 0 of 2; beta 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 7IT6 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| AF_AFQ8W0Z7F1 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |

## Copies' sheets that differ, left out of A15 (ruling C2b-3/R2: reported, not judged)

| id | view | tie | sheets left out | judged |
|---|---|---|---|---|
| 9GTP | focus | 1G~1l | copies' sheets differ: 1G sheets 92+93 vs 1l sheet 171; 8 strands not judged | yes, 37 matched |
| 9GTP | focus | 1H~1l | copies' sheets differ: 1H sheet 126 vs 1l sheets 164+165; 10 strands not judged | yes, 37 matched |
| 9GTP | focus | 1J~1l | copies' sheets differ: 1J sheets 207+208 vs 1l sheet 171; 9 strands not judged | yes, 35 matched |
| 9GTP | focus | 1K~1l | copies' sheets differ: 1K sheet 7 vs 1l sheets 164+165; 9 strands not judged | yes, 35 matched |
| 9GTP | focus | 1L~1l | copies' sheets differ: 1L sheets 54+55 vs 1l sheet 171; 8 strands not judged | yes, 37 matched |
| 9GTP | focus | 2H~2l | copies' sheets differ: 2H sheet 137 vs 2l sheets 175+177; 12 strands not judged | yes, 20 matched |
| 9GTP | focus | 2J~2l | copies' sheets differ: 2J sheet 211 vs 2l sheets 175+177; 2J sheet 213 vs 2l sheet 178; 18 strands not judged | yes, 17 matched |
| 9GTP | focus | 2K~2l | copies' sheets differ: 2K sheet 21 vs 2l sheet 178; 6 strands not judged | yes, 23 matched |
| 9GTP | focus | 2L~2l | copies' sheets differ: 2L sheet 62 vs 2l sheet 178; 6 strands not judged | yes, 22 matched |
| 9GTP | focus | 2e~2j | copies' sheets differ: 2e sheets 111+113 vs 2j sheet 72; 16 strands not judged | yes, 16 matched |
| 9GTP | focus | 2f~2j | copies' sheets differ: 2f sheets 149+151 vs 2j sheet 72; 16 strands not judged | yes, 15 matched |
| 9GTP | focus | 2g~2j | copies' sheets differ: 2g sheets 187+189 vs 2j sheet 72; 17 strands not judged | yes, 13 matched |
| 9GTP | focus | 2h~2j | copies' sheets differ: 2h sheets 221+223 vs 2j sheet 72; 16 strands not judged | yes, 16 matched |
| 9GTP | focus | 2i~2j | copies' sheets differ: 2i sheets 31+33 vs 2j sheet 72; 16 strands not judged | yes, 15 matched |
| 9GTP | focus | 3e~3j | copies' sheets differ: 3e sheets 115+116 vs 3j sheet 75; 10 strands not judged | yes, 18 matched |
| 9GTP | focus | 3h~3j | copies' sheets differ: 3h sheets 225+226 vs 3j sheet 75; 11 strands not judged | yes, 18 matched |
| 9GTP | focus | 3i~3j | copies' sheets differ: 3i sheets 35+36 vs 3j sheet 75; 11 strands not judged | yes, 16 matched |
| 9GTP | focus | M~Q | copies' sheets differ: M sheet 135 vs Q sheet 58; 6 strands not judged | yes, 15 matched |
| 9GTP | focus | N~Q | copies' sheets differ: N sheets 135+173 vs Q sheets 17+58; 13 strands not judged | yes, 12 matched |
| 9GTP | focus | b~a | copies' sheets differ: b sheets 30+70 vs a sheet 229; 9 strands not judged | yes, 4 matched |
| 9GTP | focus | l~k | copies' sheets differ: l sheets 119+120 vs k sheet 78; 12 strands not judged | n/a: sheets differ (2 matched) |
| 9GTP | focus | m~k | copies' sheets differ: m sheets 156+158 vs k sheet 78; m sheet 120 vs k sheet 41; 13 strands not judged | n/a: sheets differ (2 matched) |
| 9GTP | focus | n~k | copies' sheets differ: n sheets 1+195 vs k sheet 78; 11 strands not judged | n/a: sheets differ (2 matched) |
| 9GTP | focus | o~k | copies' sheets differ: o sheets 0+2 vs k sheet 78; o sheet 1 vs k sheet 41; 13 strands not judged | n/a: sheets differ (1 matched) |

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
- **A18** membrane: TM helices in the band, loops out of the slab
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
- documented warnings: none
