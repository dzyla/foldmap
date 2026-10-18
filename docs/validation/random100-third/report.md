# foldmap random-100 campaign

83 of 100 pass.

An entry passes when no criterion fails. A criterion with nothing to judge in an entry is not judged (or does not apply), never counted as a pass: see the outcomes per criterion below.

Sample: seed 20261018, release-date cutoff 2026-10-13.

| stratum | entries | pass |
|---|---|---|
| S1 | 12 | 12 |
| S2 | 12 | 10 |
| S3 | 10 | 9 |
| S4 | 15 | 13 |
| S5 | 8 | 8 |
| S6 | 8 | 7 |
| S7 | 7 | 3 |
| S8 | 7 | 2 |
| S9 | 3 | 3 |
| S10 | 5 | 5 |
| S11 | 10 | 8 |
| S12 | 3 | 3 |

## Outcomes per criterion

✓ pass, ✗ fail, ○ not judged (the criterion applies but the entry gave it nothing to judge: no tie, no tether to another chain, no domain drawn, mkdssp unmeasured ...), - does not apply.

| criterion | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| A1 no crash: exit 0, files written, no traceback, no Python warning | 100 | 0 | 0 | 0 |
| A1c a CA-only or P-only entry exits non-zero with foldmap's one-line message | 0 | 0 | 0 | 100 |
| A2 default view time | 100 | 0 | 0 | 0 |
| A2b the app draws it within the hosted limits | 99 | 1 | 0 | 0 |
| A3 no route fallbacks | 96 | 4 | 0 | 0 |
| A4 no loop over a label | 98 | 2 | 0 | 0 |
| A5 no element or stem boxes overlapping | 100 | 0 | 0 | 0 |
| A6 no text collisions | 94 | 6 | 0 | 0 |
| A7 text and legend inside the page | 100 | 0 | 0 | 0 |
| A8 labels unique per chain, stem names per figure | 100 | 0 | 0 | 0 |
| A9 nothing left out without a note | 100 | 0 | 0 | 0 |
| A10 secondary structure agrees with mkdssp | 97 | 0 | 0 | 3 |
| A11 N->C and sheet order | 97 | 0 | 0 | 3 |
| A12 same figure for any PYTHONHASHSEED | 100 | 0 | 0 | 0 |
| A13 a saved layout redraws the same figure | 100 | 0 | 0 | 0 |
| A14 interactive page and sequence plot written | 100 | 0 | 0 | 0 |
| A15 copies drawn alike | 40 | 3 | 13 | 44 |
| A16 domains: one chain, alike, not tiny, panels apart, names clear | 48 | 1 | 13 | 38 |
| A17 nucleic-acid metrics and stem names by their letters | 9 | 7 | 1 | 83 |
| A18 membrane: TM helices in the band, loops out of the slab | 7 | 1 | 0 | 92 |
| A19 two states: same names in the same order | 8 | 2 | 0 | 90 |
| A20 ligand tethers local | 7 | 0 | 59 | 34 |

### Not judged, per criterion

- **A15** (13): 5S4U (copies drawn with no symmetry tie: default C~A, D~B); 4B2T (copies drawn with no symmetry tie: default z~Z); 7XBK (copies drawn with no symmetry tie: default B~A, C~A, D~A, E~A, F~A, G~A, H~A, I~A); 6TIK (copies drawn with no symmetry tie: default D~C, B~C, A~C); 9BZD (copies drawn with no symmetry tie: default B~A, D~C); 6ZTQ (copies drawn with no symmetry tie: default U~T); 5NF8 (copies drawn with no symmetry tie: default B~A); 6LAB (copies drawn with no symmetry tie: default E~A, F~B, G~C, H~D, K~A, L~B, M~C, N~D, O~A, P~B, Q~C, R~D, V~U; focus E~A, F~B, G~C, H~D, K~A, L~B, M~C, N~D, O~A, P); 6ZY7 (copies drawn with no symmetry tie: default A~B); 8XJV (copies drawn with no symmetry tie: default B~A, C~A, D~A, E~A, F~A, G~A, H~A, I~A, K~A, L~A, M~J, N~J, O~J, P~J, Q~J, R~J, S~J, T~J, U~J, V~J, W~J, Y~X, Z~X, a~); 6JCY (copies drawn with no symmetry tie: default B~A); 7KMZ (copies drawn with no symmetry tie: default A~C, B~C, D~E); 5XMX (copies drawn with no symmetry tie: default B~A)
- **A16** (13): 9PDH (the domains view drew no domain); 9UY1 (the domains view drew no domain); 4YP6 (the domains view drew no domain); 6PO4 (the domains view drew no domain); 3KGS (the domains view drew no domain); 1VSK (the domains view drew no domain); 6A7T (the domains view drew no domain); 6TIK (the domains view drew no domain); 8XJV (the domains view drew no domain); 8S8X (the domains view drew no domain); 3EGZ (the domains view drew no domain); 3FT8 (the domains view drew no domain); 5XMX (the domains view drew no domain)
- **A17** (1): 8S8X (no view draws the nucleic acid)
- **A20** (59): 9DU3 (no tether to another chain; reported: ligand-tether:A1BCX:A1401:3 7.214 (home A1)); 9SP9 (no tether to another chain); 10LI (no tether to another chain; reported: ligand-tether:AG2:A1509:0 6.504 (home A1); ligand-tether:AG2:A1509:1 7.822 (home A1); ligand-tether:AG2:A1509:2 10.767 (ho); 5ES7 (no tether to another chain; reported: ligand-tether:VAL:A702:0 7.373 (home A)); 4GWY (no tether to another chain; reported: ligand-tether:AMP:A1406:3 7.816 (home A1); ligand-tether:MG:A1404:2 9.803 (home A1); ligand-tether:F6P:A1401:2 9.048 (home); 2AIQ (no tether to another chain); 9GET (no tether to another chain); 1N6K (no tether to another chain); 2LYO (no tether to another chain); 5JNN (no tether to another chain; reported: ligand-tether:6LM:A800:2 8.175 (home A)); 3CYR (no tether to another chain); 9PDH (no tether to another chain; reported: ligand-tether:ADN:A502:1 11.389 (home A); ligand-tether:ADN:A502:2 7.952 (home A); ligand-tether:ADN:A502:3 6.776 (home A)); 9UY1 (no tether to another chain); 4YJ5 (no tether to another chain; reported: ligand-tether:SER:A603:0 7.442 (home A); ligand-tether:SER:A603:1 9.117 (home A)); 3WAU (no tether to another chain; reported: ligand-tether:M1P:A1404:2 6.582 (home A1)); 4YP6 (no tether to another chain); 6PO4 (no tether to another chain; reported: ligand-tether:BO3:A303:2 13.562 (home A); ligand-tether:BO3:B303:2 13.565 (home B)); 1NID (no tether to another chain; reported: ligand-tether:CU:A1501:2 7.209 (home A1); ligand-tether:CU:A2501:2 7.794 (home A2); ligand-tether:CU:A3501:2 7.794 (home A); 1GQ2 (no tether to another chain; reported: ligand-tether:OXL:A1583:1 7.004 (home A); ligand-tether:NAP:A1581:3 8.816 (home A)); 5S4U (no tether to another chain; reported: ligand-tether:GTP:A501:3 9.933 (home A); ligand-tether:GDP:B501:3 8.126 (home B); ligand-tether:GDP:D501:3 6.77 (home D); ); 7F1V (no tether to another chain; reported: ligand-tether:HCS:D401:0 6.521 (home D); ligand-tether:HCS:D401:2 6.607 (home D); ligand-tether:HCS:D401:3 6.529 (home D)); 3ABR (no tether to another chain); 6A7T (no tether to another chain); 1XZ4 (no tether to another chain); 4ODQ (no tether to another chain); 4ZH4 (no tether to another chain; reported: ligand-tether:4OE:C2001:1 7.004 (home C)); 9Y70 (no tether to another chain; reported: ligand-tether:AQV:A601:1 10.343 (home A); ligand-tether:AQV:D601:1 9.345 (home D); ligand-tether:AQV:B601:1 9.345 (home B)); 9BZD (no tether to another chain; reported: ligand-tether:TTP:A803:1 8.167 (home A); ligand-tether:TTP:A803:2 6.654 (home A); ligand-tether:TTP:B801:3 11.809 (home B)); 8U1N (no tether to another chain; reported: ligand-tether:ATP:A801:3 10.187 (home A); ligand-tether:ATP:B801:3 7.863 (home B)); 7U9Z (no tether to another chain; reported: ligand-tether:CA:A5003:1 6.564 (home A); ligand-tether:CA:A5003:3 8.338 (home A); ligand-tether:XAN:A5004:0 11.578 (home A); 9HEO (no tether to another chain; reported: ligand-tether:CFF:A5103:0 6.814 (home A); ligand-tether:CFF:A5103:3 6.824 (home A); ligand-tether:ATP:A5102:3 7.943 (home ); 7Y63 (no tether to another chain); 6QUE (no tether to another chain; reported: ligand-tether:NAD:A1102:1 7.189 (home A); ligand-tether:NAD:A1102:2 7.101 (home A); ligand-tether:NAD:B1101:2 7.791 (home ); 6UW4 (no tether to another chain); 6VBG (no tether to another chain; reported: ligand-tether:BNG:A702:3 11.791 (home A); ligand-tether:YIO+GAL:E1:1 10.932 (home A); ligand-tether:YIO+GAL:E1:3 12.813 (h); 4D5B (no tether to another chain); 8GG4 (no tether to another chain; reported: ligand-tether:G1I:R501:3 7.037 (home R)); 6DQN (no tether to another chain; reported: ligand-tether:I3P:A3002:2 7.381 (home A); ligand-tether:I3P:A3002:3 6.831 (home A)); 6ZY7 (no tether to another chain; reported: ligand-tether:ANP:A1601:1 9.928 (home A); ligand-tether:ANP:A1601:3 8.56 (home A); ligand-tether:ANP:B1601:1 7.373 (home B); 8VBG (no tether to another chain); 3W96 (no tether to another chain); 6VTX (no tether to another chain); 6JCY (no tether to another chain); 3G6E (no tether to another chain; reported: ligand-tether:CD:U8701:2 9.386 (home U)); 8S8X (no tether to another chain); 3EGZ (no tether to another chain; reported: ligand-tether:CTC:B1601:1 26.812 (home B1)); 7F5S (no tether to another chain; reported: ligand-tether:ZN:Lg202:0 9.772 (home Lg); ligand-tether:ZN:Lg202:1 8.999 (home Lg); ligand-tether:ZN:Lg202:2 7.845 (home L); 6QIK (no tether to another chain; reported: ligand-tether:ZN:g201:2 6.751 (home g); ligand-tether:ZN:g201:3 7.433 (home g)); 2FCY (no tether to another chain); 6IZP (no tether to another chain); 7UKO (no tether to another chain; reported: ligand-tether:MN:B2003:4 6.523 (home B)); 7G33 (no tether to another chain; reported: ligand-tether:XHQ:A902:1 6.728 (home A)); 8EM7 (no tether to another chain); 8QX1 (no tether to another chain; reported: ligand-tether:LMR:A11001:0 19.576 (home A1); ligand-tether:LMR:A11001:1 15.953 (home A1); ligand-tether:LMR:A11001:2 13.38); 8CRC (no tether to another chain); 2YEJ (no tether to another chain; reported: ligand-tether:ZZ3:A1226:1 8.076 (home A)); 3FT8 (no tether to another chain; reported: ligand-tether:MOJ:A1237:0 11.229 (home A1); ligand-tether:MOJ:A2237:0 11.229 (home A2)); 5U3R (no tether to another chain; reported: ligand-tether:7V1:A502:3 8.193 (home A); ligand-tether:B7G:A503:2 7.297 (home A)); 5XMX (no tether to another chain; reported: ligand-tether:89L:B900:2 6.775 (home B))

## Entries

process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 gates; drawing s: the drawing inside it.

| id | stratum | bin | chains | process s | drawing s | A1 | A1c | A2 | A2b | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 | A11 | A12 | A13 | A14 | A15 | A16 | A17 | A18 | A19 | A20 | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 9DU3 | S1 | M | 1 | 1.05 | 0.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9SP9 | S1 | M | 1 | 0.85 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 10LI | S1 | M | 1 | 1.55 | 0.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 5ES7 | S1 | M | 1 | 1.2 | 0.54 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 4GWY | S1 | M | 1 | 1.9 | 1.06 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 1RKA | S1 | M | 1 | 1.35 | 0.76 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 2AIQ | S1 | S | 1 | 0.9 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 9GET | S1 | S | 1 | 0.7 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 1N6K | S1 | S | 1 | 0.75 | 0.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2LYO | S1 | S | 1 | 1.05 | 0.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5JNN | S1 | S | 1 | 0.75 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3CYR | S1 | S | 1 | 0.7 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 9PDH | S2 | L | 4 | 2.51 | 1.53 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ○ | - | - | - | ○ | A15: default 9.67976255751613 |
| 9UY1 | S2 | L | 18 | 1.9 | 1.09 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 4YJ5 | S2 | L | 4 | 1.25 | 0.51 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | - | ○ | A15: focus 4.185786228121103 |
| 3EOM | S2 | L | 4 | 2.05 | 1.29 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 3WAU | S2 | M | 2 | 1.1 | 0.48 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ○ |  |
| 4YP6 | S2 | M | 3 | 1.6 | 0.74 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 6PO4 | S2 | M | 6 | 1.0 | 0.38 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 1NID | S2 | M | 1 | 1.5 | 0.76 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 3KGS | S2 | S | 2 | 0.95 | 0.34 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 1VSK | S2 | S | 1 | 0.85 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5C1A | S2 | XL | 12 | 2.1 | 1.06 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ |  |
| 1GQ2 | S2 | XL | 16 | 1.5 | 0.71 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ○ |  |
| 5S4U | S3 | L | 6 | 2.91 | 1.94 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 7F1V | S3 | L | 4 | 1.85 | 1.09 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8PIN | S3 | L | 7 | 1.1 | 0.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | - | ✓ | A15: focus 18.877391431648665 |
| 3ABR | S3 | L | 4 | 1.65 | 0.82 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6A7T | S3 | M | 2 | 1.35 | 0.79 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 1XZ4 | S3 | M | 4 | 1.05 | 0.43 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2YNR | S3 | M | 3 | 1.0 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 4ODQ | S3 | S | 2 | 0.7 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 4B2T | S3 | XL | 16 | 7.66 | 5.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | - |  |
| 4ZH4 | S3 | XL | 12 | 4.36 | 3.08 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 7XBK | S4 | L | 10 | 3.66 | 2.53 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 9Y70 | S4 | L | 4 | 2.3 | 1.5 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6TIK | S4 | L | 4 | 1.75 | 1.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | - | - | - | - |  |
| 12DV | S4 | L | 9 | 1.15 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 9BZD | S4 | L | 4 | 2.81 | 1.84 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 8U1N | S4 | M | 2 | 1.4 | 0.7 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9N2E | S4 | M | 2 | 1.0 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 6ZTQ | S4 | XL | 45 | 26.09 | 23.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 9NMN | S4 | XL | 8 | 9.52 | 7.33 | ✓ | - | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | A4: default: A:4908-4954>A:4963-4967; deposited: A:4908-4954>A:4963-4967; beta: A:4908-4954>A:4963-4967; A6: default 1 (label:A:4963-4967 × loop:A:4908-4954>A:4963-4967); deposited 1 (label:A:4963-4967 × loop:A:4908-4954>A:4963-4967); beta 1 (label:A:4963-4967 × loop:A:4908-4954>A:4963-4967) |
| 8FZN | S4 | XL | 60 | 8.41 | 6.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 7SPI | S4 | XL | 78 | 1.7 | 0.79 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 7U9Z | S4 | XL | 8 | 7.11 | 5.03 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9HEO | S4 | XL | 12 | 8.46 | 6.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7ADZ | S4 | XL | 30 | 2.35 | 1.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 9FZL | S4 | XL | 98 | 59.45 | 52.2 | ✓ | - | ✓ | ✓ | ✗ | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | - | - | - | ✓ | A3: domains: 4E:88-89>4E:119-131, 4K:15-17>4K:33-35, 4e:88-89>4e:119-131, 4k:15-17>4k:33-35; A4: domains: 2M:113-132>2M:147-154, 3A:274-282>3A:291-301, 3B:520-525>3B:533-541, 3C:234-267>3C:270-272; A6: domains 5 (domain-label:D1 · 3A × loop:3A:274-282>3A:291-301; domain-label:D2 · 3B × loop:3B:520-525>3B:533-541; domain-label:D1 · 3a × loop:3a:274-282>3a:291-301); A16: domains 4 (domain-label:D1 · 3A × loop:3A:274-282>3A:291-301; domain-label:D2 · 3B × loop:3B:520-525>3B:533-541; domain-label:D1 · 3a × loop:3a:274-282>3a:291-301) |
| 5NF8 | S5 | M | 2 | 1.0 | 0.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | - |  |
| 1UG2 | S5 | S | 1 | 0.7 | 0.13 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2LLF | S5 | S | 1 | 0.65 | 0.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2KQV | S5 | S | 1 | 0.8 | 0.18 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2MAH | S5 | S | 1 | 1.2 | 0.59 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2KVK | S5 | S | 1 | 0.7 | 0.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2MXC | S5 | S | 1 | 0.8 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6Y1H | S5 | S | 1 | 0.65 | 0.1 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 7Y63 | S6 | L | 2 | 3.01 | 2.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ○ |  |
| 6BO8 | S6 | L | 4 | 2.3 | 1.46 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - |  |
| 6QUE | S6 | L | 2 | 2.5 | 1.67 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 6UW4 | S6 | L | 4 | 2.55 | 1.57 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ○ |  |
| 6VBG | S6 | M | 4 | 1.45 | 0.77 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 4D5B | S6 | M | 2 | 0.81 | 0.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 8GG4 | S6 | M | 4 | 1.25 | 0.58 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 6DQN | S6 | XL | 4 | 3.56 | 2.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | - | ○ | A18: no band drawn |
| 6LAB | S7 | L | 22 | 6.41 | 2.99 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | - | - | ✓ |  |
| 6ZY7 | S7 | L | 6 | 4.71 | 1.93 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✗ | - | - | ○ | A17: stem_recall 0.667 (limit 0.85) |
| 9XSF | S7 | M | 10 | 3.56 | 1.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - |  |
| 8VBG | S7 | M | 3 | 2.66 | 0.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ○ | A17: contact_auc 0.505 (limit 0.6); stems_nearer_frac 0.5 (limit 1.0) |
| 3W96 | S7 | M | 10 | 3.91 | 1.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ |  |
| 6VTX | S7 | S | 3 | 1.9 | 0.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | ○ | A17: stems_nearer_frac 0.0 (limit 1.0) |
| 8XJV | S7 | XL | 110 | 129.02 | 114.5 | ✓ | - | ✓ | ✗ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | ✗ | - | - | - | A2b: time: 8XJV took longer than 90 s to draw, the limit of this app. Draw part of it (Chains & assembly), or on your own computer: pip install foldmap, then foldmap plot 8XJV -o 8XJV.svg.; A3: default: strand Au:348-503, strand Au:536-615, strand Au:886-1050, strand Au:1271-1381; deposited: strand Au:348-503, strand Au:536-615, strand Au:886-1050, strand Au:1271-1381; beta: strand Au:348-503, strand Au:536-615, strand Au:886-1050, strand Au:1271-1381; focus: strand Au:348-503, strand Au:536-615, strand Au:886-1050, strand Au:1271-1381; domains: strand Au:348-503, strand Au:536-615, strand Au: |
| 6JCY | S8 | L | 9 | 5.61 | 2.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | - | - | ○ |  |
| 3G6E | S8 | L | 31 | 35.46 | 23.89 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | ○ | A6: default 1 (na:tie:168-219:label × na:tie:164-171); deposited 1 (na:tie:168-219:label × na:tie:164-171); beta 1 (na:tie:168-219:label × na:tie:164-171); focus 1 (na:tie:168-219:label × na:tie:164-171) |
| 8S8X | S8 | M | 3 | 1.1 | 0.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | ○ | - | - | ○ |  |
| 4AS1 | S8 | M | 2 | 3.16 | 0.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | - | A17: na_crossings 1 (limit 0) |
| 3EGZ | S8 | S | 2 | 3.16 | 1.08 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✗ | - | - | ○ | A17: pair_f1 0.941 (limit 0.95); spearman 0.347 (limit 0.6); procrustes 0.982 (limit 0.5); angle_median 22.0 (limit 0.05); angle_max 125.7 (limit 0.05); mirrored 1 (limit 0) |
| 7F5S | S8 | XL | 45 | 58.1 | 41.44 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | ○ | A3: domains: strand L5:694-695; A6: domains 2 (na:tie:2575-2757:label × na:tie:2738-2742; terminus:N:Lr × na:tie:494-659:mark) |
| 6QIK | S8 | XL | 47 | 45.58 | 31.61 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ○ | A3: default: strand y:8-9; deposited: strand y:8-9; beta: strand y:8-9; focus: strand y:8-9; domains: strand y:8-9; A6: default 3 (na:tie:42-92:label × na:tie:38-45; na:tie:190-221:label × na:tie:189-205; terminus:N:i × na:tie:22-35:mark); deposited 3 (na:tie:42-92:label × na:tie:38-45; na:tie:190-221:label × na:tie:189-205; terminus:N:i × na:tie:22-35:mark); beta 3 (na:tie:42-92:label × na:tie:38-45; na:tie:190-221:label × na:tie:189-205; terminus:N:i × na:tie:22-35:mark); focus 3 (na:tie:42-92:label × na:tie:38-45; na:tie:190-221:label × na:tie:189-205; terminus:N:i × na:tie:22-35:mark); A17 |
| 9OWM | S9 | M | 2 | 3.91 | 1.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 2FCY | S9 | S | 2 | 1.85 | 0.18 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 6IZP | S9 | S | 1 | 2.3 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 7UKO | S10 | L | 8 | 6.81 | 5.95 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 7KMZ | S10 | L | 5 | 35.12 | 33.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | - |  |
| 7G33 | S10 | M | 1 | 4.01 | 3.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3O6N | S10 | M | 1 | 1.25 | 0.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8EM7 | S10 | XL | 2 | 36.86 | 30.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8OJZ | S11 | L | 2 | 1.9 | 1.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✗ | - | A19: names compared on all 48 pairs, 14 differ; αNa~αNb (A1:32-49~A1:30-49); αNb~αNc (A1:51-70~A1:51-70); αNc~αNd (A1:74-84~A1:73-84); αNd~αNe (A1:88-118~A1:88-115); αNe~αNf (A1:139-148~A1:139-148); αNf~αNg (A1:154-162~A1:154-163); α1~α1″ (A1:203-221~A1:203-221); α1″~α1‴ (A1:233-241~A1:233-240) (+11 more) |
| 8QX1 | S11 | L | 3 | 1.65 | 0.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✗ | ○ | A19: names compared on all 48 pairs, 14 differ; αNb~αNa (A1:30-49~A1:32-49); αNc~αNb (A1:51-70~A1:51-70); αNd~αNc (A1:73-84~A1:74-84); αNe~αNd (A1:88-115~A1:88-118); αNf~αNe (A1:139-148~A1:139-148); αNg~αNf (A1:154-163~A1:154-162); α1″~α1 (A1:203-221~A1:203-221); α1‴~α1″ (A1:233-240~A1:233-241) (+11 more) |
| 3P2W | S11 | S | 1 | 0.8 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | - |  |
| 8CRC | S11 | S | 1 | 0.8 | 0.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 2YEJ | S11 | S | 1 | 0.8 | 0.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 3FT8 | S11 | S | 1 | 1.05 | 0.34 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | ✓ | ○ |  |
| 5U3R | S11 | M | 2 | 1.0 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 5XMX | S11 | M | 2 | 1.15 | 0.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | - | - | ✓ | ○ |  |
| 6M67 | S11 | L | 7 | 1.2 | 0.56 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - |  |
| 6LTO | S11 | L | 7 | 1.0 | 0.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - |  |
| AF_AFQ8ZPJ0F1 | S12 | M | 1 | 1.2 | 0.43 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFQ5AGV5F1 | S12 | M | 1 | 1.35 | 0.66 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFI1K7U1F1 | S12 | S | 1 | 0.75 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |

## Two-state pairs (S11)

| pair | seed | partner | sampled as | types now | ligand touching | assembly 1 chains | CA RMSD Å | A19 |
|---|---|---|---|---|---|---|---|---|
| P1 | 8OJZ | 8QX1 | apo/holo | apo/holo | [False, True] | [4, 4] | 3.415 | ✗ names compared on all 48 pairs, 14 differ; αNa~αNb (A1:32-49~A1:30-49); αNb~αNc (A1:51-70~A1:51-70); αNc~αNd (A1:74-84~A1:73-84); αNd~αNe (A1:88-118~A1:88-115); αNe~αNf (A1:139-148~A1:139-148); αNf~αNg (A1:154-162~A1:154-163); α1~α1″ (A1:203-221~A1:203-221); α1″~α1‴ (A1:233-241~A1:233-240) (+11 more |
| P2 | 3P2W | 8CRC | apo/holo | apo/holo | [False, True] | [1, 1] | 3.839 | ✓ names compared on all 18 pairs |
| P3 | 2YEJ | 3FT8 | monomer/oligomer | monomer/oligomer | [True, True] | [1, 2] | 0.285 | ✓ names compared on all 16 pairs |
| P4 | 5U3R | 5XMX | monomer/oligomer | monomer/oligomer | [True, True] | [1, 2] | 0.426 | ✓ names compared on all 15 pairs; A:268-271 unpaired; A:268-272 unpaired |
| P5 | 6M67 | 6LTO | open/closed | open/closed | [False, False] | [7, 7] | 8.376 | ✓ names compared on all 11 pairs; A:104-106 unpaired; A:249-255 unpaired; A:279-285 unpaired; A:288-295 unpaired; A:297-299 unpaired; A:312-315 unpaired; A:317-319 unpaired; A:122-130 unpaired (+3 more) |

## Reported, not gated

| id | honest_crossings | loop_crossings | na_run_crossings | fallback_curves | fallback_curves_crossing_an_element | kink_links_rejected | tighter_than_inflate | rejected_as_crossing | names_moved | area | max_rss_mb |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 9DU3 | 20 | 20 | None | 0 | 0 | 0 | 0 | 0 | 0 | 552.51 | 167.6 |
| 9SP9 | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 374.1 | 162.5 |
| 10LI | 51 | 51 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1279.58 | 205.5 |
| 5ES7 | 47 | 47 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1004.49 | 172.6 |
| 4GWY | 74 | 74 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1942.44 | 211.3 |
| 1RKA | 108 | 108 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1121.25 | 192.1 |
| 2AIQ | 5 | 5 | None | 0 | 0 | 0 | 0 | 0 | 0 | 225.48 | 160.6 |
| 9GET | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 104.86 | 159.9 |
| 1N6K | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 205.53 | 161.2 |
| 2LYO | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 105.39 | 160.1 |
| 5JNN | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 178.45 | 160.9 |
| 3CYR | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 64.64 | 159.9 |
| 9PDH | 101 | 101 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2809.4 | 227.3 |
| 9UY1 | 156 | 156 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1575.42 | 233.0 |
| 4YJ5 | 15 | 15 | None | 0 | 0 | 0 | 0 | 0 | 0 | 689.31 | 260.7 |
| 3EOM | 67 | 67 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2719.57 | 296.0 |
| 3WAU | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 664.0 | 208.3 |
| 4YP6 | 19 | 19 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1338.08 | 173.9 |
| 6PO4 | 18 | 18 | None | 0 | 0 | 0 | 0 | 0 | 0 | 524.5 | 164.9 |
| 1NID | 103 | 103 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1197.48 | 205.8 |
| 3KGS | 26 | 26 | None | 0 | 0 | 0 | 0 | 0 | 0 | 462.27 | 163.5 |
| 1VSK | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 352.68 | 162.7 |
| 5C1A | 47 | 47 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1526.82 | 425.5 |
| 1GQ2 | 52 | 52 | None | 0 | 0 | 0 | 1 | 0 | 0 | 676.79 | 261.5 |
| 5S4U | 140 | 140 | None | 0 | 0 | 0 | 2 | 0 | 0 | 4109.34 | 341.9 |
| 7F1V | 130 | 130 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1817.19 | 244.7 |
| 8PIN | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 454.23 | 198.9 |
| 3ABR | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1132.73 | 308.6 |
| 6A7T | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 344.61 | 162.3 |
| 1XZ4 | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 693.86 | 167.0 |
| 2YNR | 8 | 8 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1076.03 | 165.5 |
| 4ODQ | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 109.5 | 160.9 |
| 4B2T | 681 | 681 | None | 0 | 0 | 0 | 1 | 0 | 0 | 8940.54 | 557.6 |
| 4ZH4 | 316 | 316 | None | 0 | 0 | 0 | 1 | 0 | 0 | 5011.79 | 464.6 |
| 7XBK | 247 | 247 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3654.12 | 393.8 |
| 9Y70 | 78 | 78 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3779.14 | 285.5 |
| 6TIK | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 672.19 | 164.7 |
| 12DV | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 695.16 | 245.4 |
| 9BZD | 214 | 214 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2387.88 | 272.8 |
| 8U1N | 44 | 44 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1723.55 | 231.2 |
| 9N2E | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 720.24 | 179.5 |
| 6ZTQ | 251 | 251 | None | 0 | 0 | 0 | 1 | 0 | 0 | 16353.65 | 522.3 |
| 9NMN | 300 | 300 | None | 0 | 0 | 0 | 1 | 0 | 0 | 10388.35 | 496.9 |
| 8FZN | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 507.36 | 734.6 |
| 7SPI | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 426.01 | 320.4 |
| 7U9Z | 301 | 301 | None | 0 | 0 | 0 | 2 | 0 | 0 | 6539.89 | 491.2 |
| 9HEO | 243 | 243 | None | 0 | 0 | 0 | 0 | 0 | 0 | 7962.12 | 504.8 |
| 7ADZ | 67 | 67 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2187.08 | 317.2 |
| 9FZL | 707 | 707 | None | 0 | 0 | 0 | 11 | 0 | 22 | 65458.72 | 2297.0 |
| 5NF8 | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 371.87 | 249.9 |
| 1UG2 | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 75.21 | 188.8 |
| 2LLF | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 129.18 | 163.5 |
| 2KQV | 7 | 7 | None | 0 | 0 | 0 | 0 | 0 | 0 | 110.91 | 200.4 |
| 2MAH | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 113.08 | 201.0 |
| 2KVK | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 140.51 | 169.4 |
| 2MXC | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 108.62 | 241.9 |
| 6Y1H | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 43.1 | 160.6 |
| 7Y63 | 64 | 64 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1783.14 | 246.9 |
| 6BO8 | 42 | 42 | None | 0 | 0 | 0 | 4 | 0 | 0 | 2350.51 | 315.9 |
| 6QUE | 158 | 158 | None | 0 | 0 | 0 | 1 | 0 | 0 | 3428.0 | 324.1 |
| 6UW4 | 16 | 16 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3260.28 | 270.7 |
| 6VBG | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 834.19 | 177.4 |
| 4D5B | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 321.83 | 164.3 |
| 8GG4 | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1223.52 | 191.8 |
| 6DQN | 73 | 73 | None | 0 | 0 | 0 | 1 | 0 | 0 | 4941.13 | 462.2 |
| 6LAB | 125 | 125 | 31 | 0 | 0 | 0 | 1 | 0 | 0 | 7149.65 | 376.5 |
| 6ZY7 | 154 | 154 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 3857.54 | 354.5 |
| 9XSF | 45 | 45 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 2382.08 | 212.3 |
| 8VBG | 28 | 28 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 1285.96 | 224.5 |
| 3W96 | 37 | 37 | 5 | 0 | 0 | 0 | 2 | 0 | 0 | 3042.22 | 236.2 |
| 6VTX | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 151.56 | 165.9 |
| 8XJV | 1271 | 1271 | 191 | 0 | 0 | 0 | 1 | 0 | 19 | 43551.94 | 591.4 |
| 6JCY | 245 | 245 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4916.37 | 449.7 |
| 3G6E | 382 | 382 | 395 | 0 | 0 | 0 | 0 | 0 | 0 | 19741.12 | 603.1 |
| 8S8X | 38 | 38 | None | 0 | 0 | 0 | 0 | 0 | 0 | 549.55 | 165.8 |
| 4AS1 | 29 | 29 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1609.99 | 193.6 |
| 3EGZ | 8 | 8 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 1478.07 | 182.4 |
| 7F5S | 458 | 458 | 503 | 0 | 0 | 0 | 2 | 0 | 0 | 29136.7 | 843.8 |
| 6QIK | 425 | 425 | 406 | 0 | 0 | 0 | 2 | 0 | 0 | 27610.17 | 716.6 |
| 9OWM | 0 | 0 | 23 | 0 | 0 | 0 | 0 | 0 | 0 | 2626.22 | 209.3 |
| 2FCY | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 167.21 | 164.5 |
| 6IZP | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 85.29 | 179.2 |
| 7UKO | 89 | 89 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1698.78 | 233.5 |
| 7KMZ | 262 | 262 | None | 0 | 0 | 0 | 1 | 0 | 0 | 4784.68 | 450.6 |
| 7G33 | 48 | 48 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1088.86 | 169.6 |
| 3O6N | 106 | 106 | None | 0 | 0 | 0 | 0 | 0 | 0 | 255.94 | 163.9 |
| 8EM7 | 364 | 364 | None | 0 | 0 | 0 | 2 | 0 | 0 | 8203.86 | 517.0 |
| 8OJZ | 55 | 55 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1611.47 | 308.9 |
| 8QX1 | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1513.7 | 274.9 |
| 3P2W | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 180.15 | 161.2 |
| 8CRC | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 182.87 | 162.2 |
| 2YEJ | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 211.7 | 160.6 |
| 3FT8 | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 508.73 | 163.7 |
| 5U3R | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 364.98 | 164.9 |
| 5XMX | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 940.22 | 165.3 |
| 6M67 | 5 | 5 | None | 0 | 0 | 1 | 0 | 0 | 0 | 305.89 | 188.8 |
| 6LTO | 3 | 3 | None | 0 | 0 | 0 | 1 | 0 | 0 | 432.47 | 195.9 |
| AF_AFQ8ZPJ0F1 | 39 | 39 | None | 0 | 0 | 0 | 0 | 0 | 0 | 421.83 | 162.9 |
| AF_AFQ5AGV5F1 | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1956.97 | 186.3 |
| AF_AFI1K7U1F1 | 5 | 5 | None | 0 | 0 | 0 | 0 | 0 | 0 | 150.43 | 159.9 |

## Science checks (S1-S6), apart from A1-A20

Judged in every drawn view (default, deposited, beta, focus, domains, membrane, app); S5 once per entry.

A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet drawn, no barrel, no residue mark, no accessory strand, no CATH chopping) is not judged, one with no protein does not apply.

| check | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| S1 sheet neighbours pair | 88 | 0 | 9 | 3 |
| S2 barrel closure pairs | 16 | 0 | 81 | 3 |
| S3 barrel strand count is its pairing cycle | 88 | 0 | 9 | 3 |
| S4 residue marks at their residue | 0 | 0 | 97 | 3 |
| S5 automatic domains agree with a reference | 23 | 12 | 62 | 3 |
| S6 accessory drawn at its partner | 25 | 0 | 72 | 3 |

### Science by view

Entries measured in each view, and how many of them fail each check there.

| view | measured | S1 | S2 | S3 | S4 | S6 |
|---|---|---|---|---|---|---|
| default | 100 | 0 | 0 | 0 | 0 | 0 |
| deposited | 97 | 0 | 0 | 0 | 0 | 0 |
| beta | 97 | 0 | 0 | 0 | 0 | 0 |
| focus | 69 | 0 | 0 | 0 | 0 | 0 |
| domains | 62 | 0 | 0 | 0 | 0 | 0 |
| membrane | 8 | 0 | 0 | 0 | 0 | 0 |
| app | 99 | 0 | 0 | 0 | 0 | 0 |

| id | S1 | S2 | S3 | S4 | S5 | S6 | findings |
|---|---|---|---|---|---|---|---|
| 9DU3 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9SP9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 10LI | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 5ES7 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 4GWY | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 1RKA | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.75, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 2AIQ | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.52, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 9GET | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1N6K | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2LYO | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5JNN | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 3CYR | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 9PDH | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9UY1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4YJ5 | ✓ | ✓ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.74, boundary F1 0.40, recall 0.25 (3 reference domains, 2 automatic); B: overlap 0.74, boundary F1 0.40, recall 0.25 (3 reference domains, 2 automatic); C: overlap 0.74, boundary F1 0.40, recall 0.25 (3 reference domains, 2 automatic); D: overlap 0.74, boundary F1 0.40, recall 0.25 (3 reference domains, 2 automatic) |
| 3EOM | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.69, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic); B: overlap 0.70, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic); C: overlap 0.51, boundary F1 0.40, recall 0.50 (3 reference domains, 3 automatic); D: overlap 0.70, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic) |
| 3WAU | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.43, boundary F1 0.00, recall 1.00 (1 reference domains, 3 automatic); B: overlap 0.44, boundary F1 0.00, recall 1.00 (1 reference domains, 3 automatic) |
| 4YP6 | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6PO4 | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 1NID | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.83, boundary F1 0.00, recall 0.00 (2 reference domains, 2 automatic) |
| 3KGS | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 1VSK | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5C1A | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 1GQ2 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.88, boundary F1 0.67, recall 1.00 (2 reference domains, 2 automatic); B: overlap 0.54, boundary F1 0.29, recall 1.00 (2 reference domains, 5 automatic); C: overlap 0.88, boundary F1 0.67, recall 1.00 (2 reference domains, 2 automatic); D: overlap 0.50, boundary F1 0.25, recall 1.00 (2 reference domains, 5 automatic); E: overlap 0.70, boundary F1 0.40, recall 1.00 (2 reference domains, 4 automatic); F: overlap 0.50, boundary F1 0.25, recall 1.00 (2 reference domains, 5 automatic); G: overlap 0.72, boundary F1 0.40, recall 1.00 (2 reference domains, 4 automatic); H: overlap 0.71, |
| 5S4U | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7F1V | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8PIN | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3ABR | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.85, boundary F1 0.67, recall 0.50 (3 reference domains, 2 automatic); C: overlap 0.53, boundary F1 0.40, recall 0.50 (3 reference domains, 2 automatic) |
| 6A7T | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 1XZ4 | ○ | ○ | ○ | ○ | ✗ | ○ | S5 B: overlap 0.59, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); D: overlap 0.59, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 2YNR | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 4ODQ | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 4B2T | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4ZH4 | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 F: overlap 0.69, boundary F1 0.57, recall 1.00 (3 reference domains, 5 automatic); L: overlap 0.70, boundary F1 0.57, recall 1.00 (3 reference domains, 5 automatic) |
| 7XBK | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9Y70 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6TIK | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 12DV | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 9BZD | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8U1N | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9N2E | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6ZTQ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9NMN | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8FZN | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7SPI | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7U9Z | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 9HEO | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 7ADZ | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9FZL | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 5NF8 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1UG2 | ○ | ○ | ○ | ○ | ✓ | ○ |  |
| 2LLF | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2KQV | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2MAH | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2KVK | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2MXC | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6Y1H | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 7Y63 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6BO8 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6QUE | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6UW4 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6VBG | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4D5B | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 8GG4 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6DQN | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 6LAB | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6ZY7 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9XSF | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8VBG | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3W96 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 D: overlap 0.61, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); H: overlap 0.60, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 6VTX | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8XJV | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6JCY | ✓ | ✓ | ✓ | ○ | ✓ | ✓ |  |
| 3G6E | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.50, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); B: overlap 0.73, boundary F1 0.89, recall 0.80 (3 reference domains, 2 automatic); L: overlap 0.65, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic); P: overlap 0.77, boundary F1 0.67, recall 0.50 (3 reference domains, 2 automatic); W: overlap 0.64, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 8S8X | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4AS1 | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 3EGZ | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7F5S | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 6QIK | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 9OWM | - | - | - | - | - | - |  |
| 2FCY | - | - | - | - | - | - |  |
| 6IZP | - | - | - | - | - | - |  |
| 7UKO | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7KMZ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7G33 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3O6N | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8EM7 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8OJZ | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 8QX1 | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 3P2W | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 8CRC | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2YEJ | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 3FT8 | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5U3R | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5XMX | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 6M67 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6LTO | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ8ZPJ0F1 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| AF_AFQ5AGV5F1 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| AF_AFI1K7U1F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |

## S5: automatic domains against CATH (reported, not gated)

Source: CATH domain boundaries (https://www.cathdb.info) through the PDBe API (EMBL-EBI), licensed under CC BY 4.0; Sillitoe I. et al. (2021) CATH: increased structural coverage of functional space. Nucleic Acids Res 49:D266-D273.

109 chains scored in 35 entries; 47 of them with at least two CATH domains. Every chain CATH chops is scored, one-domain chains too (a one-domain chain keeps boundary F1 1 only when the automatic domains leave it whole). 0 of the chains are in the calibration set (source "reference": in-sample, tests/data/domains/cath_reference.json, the 33 chains the automatic domains were calibrated on); the never-seen check is the fresh set (cath_fresh.json).

| id | chain | source | CATH domains | automatic | overlap | boundary F1 | recall | precision |
|---|---|---|---|---|---|---|---|---|
| 4GWY | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 1RKA | A | CATH | 1 | 2 | 0.75 | 0.00 | 1.00 | 0.00 |
| 2AIQ | A | CATH | 2 | 1 | 0.52 | 0.00 | 0.00 | 0.00 |
| 1N6K | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2LYO | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5JNN | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3CYR | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4YJ5 | A | CATH | 3 | 2 | 0.74 | 0.40 | 0.25 | 1.00 |
| 4YJ5 | B | CATH | 3 | 2 | 0.74 | 0.40 | 0.25 | 1.00 |
| 4YJ5 | C | CATH | 3 | 2 | 0.74 | 0.40 | 0.25 | 1.00 |
| 4YJ5 | D | CATH | 3 | 2 | 0.74 | 0.40 | 0.25 | 1.00 |
| 3EOM | A | CATH | 3 | 2 | 0.69 | 0.50 | 0.50 | 0.50 |
| 3EOM | B | CATH | 3 | 2 | 0.70 | 0.50 | 0.50 | 0.50 |
| 3EOM | C | CATH | 3 | 3 | 0.51 | 0.40 | 0.50 | 0.33 |
| 3EOM | D | CATH | 3 | 2 | 0.70 | 0.50 | 0.50 | 0.50 |
| 3WAU | A | CATH | 1 | 3 | 0.43 | 0.00 | 1.00 | 0.00 |
| 3WAU | B | CATH | 1 | 3 | 0.44 | 0.00 | 1.00 | 0.00 |
| 4YP6 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4YP6 | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4YP6 | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1NID | A | CATH | 2 | 2 | 0.83 | 0.00 | 0.00 | 0.00 |
| 3KGS | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3KGS | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1VSK | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1GQ2 | A | CATH | 2 | 2 | 0.88 | 0.67 | 1.00 | 0.50 |
| 1GQ2 | B | CATH | 2 | 5 | 0.54 | 0.29 | 1.00 | 0.17 |
| 1GQ2 | C | CATH | 2 | 2 | 0.88 | 0.67 | 1.00 | 0.50 |
| 1GQ2 | D | CATH | 2 | 5 | 0.50 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | E | CATH | 2 | 4 | 0.70 | 0.40 | 1.00 | 0.25 |
| 1GQ2 | F | CATH | 2 | 5 | 0.50 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | G | CATH | 2 | 4 | 0.72 | 0.40 | 1.00 | 0.25 |
| 1GQ2 | H | CATH | 2 | 4 | 0.71 | 0.40 | 1.00 | 0.25 |
| 1GQ2 | I | CATH | 2 | 5 | 0.49 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | J | CATH | 2 | 5 | 0.50 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | K | CATH | 2 | 5 | 0.50 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | L | CATH | 2 | 2 | 0.88 | 0.67 | 1.00 | 0.50 |
| 1GQ2 | M | CATH | 2 | 5 | 0.50 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | N | CATH | 2 | 5 | 0.50 | 0.25 | 1.00 | 0.14 |
| 1GQ2 | O | CATH | 2 | 5 | 0.56 | 0.33 | 1.00 | 0.20 |
| 1GQ2 | P | CATH | 2 | 4 | 0.71 | 0.40 | 1.00 | 0.25 |
| 3ABR | A | CATH | 3 | 2 | 0.85 | 0.67 | 0.50 | 1.00 |
| 3ABR | B | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3ABR | C | CATH | 3 | 2 | 0.53 | 0.40 | 0.50 | 0.33 |
| 3ABR | D | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1XZ4 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1XZ4 | B | CATH | 1 | 2 | 0.59 | 0.00 | 1.00 | 0.00 |
| 1XZ4 | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1XZ4 | D | CATH | 1 | 2 | 0.59 | 0.00 | 1.00 | 0.00 |
| 4ODQ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | A | CATH | 3 | 3 | 0.99 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | B | CATH | 2 | 2 | 0.96 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | F | CATH | 3 | 5 | 0.69 | 0.57 | 1.00 | 0.40 |
| 4ZH4 | G | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | H | CATH | 2 | 2 | 0.96 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4ZH4 | L | CATH | 3 | 5 | 0.70 | 0.57 | 1.00 | 0.40 |
| 1UG2 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2LLF | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2KQV | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2MAH | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2KVK | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2MXC | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | D | CATH | 1 | 2 | 0.61 | 0.00 | 1.00 | 0.00 |
| 3W96 | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | G | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3W96 | H | CATH | 1 | 2 | 0.60 | 0.00 | 1.00 | 0.00 |
| 6JCY | A | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 6JCY | B | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 3G6E | 1 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | 2 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | 3 | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | A | CATH | 2 | 1 | 0.50 | 0.00 | 0.00 | 0.00 |
| 3G6E | B | CATH | 3 | 2 | 0.73 | 0.89 | 0.80 | 1.00 |
| 3G6E | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | E | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 3G6E | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | I | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | J | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | L | CATH | 2 | 1 | 0.65 | 0.00 | 0.00 | 0.00 |
| 3G6E | M | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | N | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | O | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | P | CATH | 3 | 2 | 0.77 | 0.67 | 0.50 | 1.00 |
| 3G6E | Q | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | R | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | S | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | T | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | U | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | V | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | W | CATH | 2 | 1 | 0.64 | 0.00 | 0.00 | 0.00 |
| 3G6E | X | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3G6E | Z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4AS1 | A | CATH | 5 | 6 | 0.81 | 0.86 | 1.00 | 0.75 |
| 3EGZ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3P2W | A | CATH | 2 | 2 | 0.96 | 1.00 | 1.00 | 1.00 |
| 2YEJ | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3FT8 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5U3R | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5U3R | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5XMX | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5XMX | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| mean of 109 (all) | | | | | 0.86 | 0.73 | 0.89 | 0.74 |
| mean of 47 (two domains or more) | | | | | 0.74 | 0.53 | 0.75 | 0.54 |
| mean of 62 (one domain) | | | | | 0.95 | 0.89 | 1.00 | 0.89 |

Protein entries with no chain scored (S5 not judged), by reason:

- CATH chops none of its chains: 59
- not a PDB entry (CATH chops PDB entries): 3

## Accessories named, not drawn level beside their partner (reported, not gated)

| id | named of accessory strands, per view |
|---|---|
| 9DU3 | default 0 of 2; deposited 0 of 2; beta 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 10LI | default 0 of 4; deposited 0 of 4; beta 0 of 4; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |
| 5ES7 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 2AIQ | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 3EOM | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 0 of 0; domains 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 6PO4 | default 2 of 10; deposited 2 of 10; beta 2 of 10; focus 2 of 10; domains 2 of 10; seed1 2 of 10; roundtrip 2 of 10 |
| 1NID | default 0 of 3; deposited 0 of 3; beta 0 of 3; domains 0 of 3; seed1 0 of 3; roundtrip 0 of 3 |
| 3KGS | default 0 of 4; deposited 0 of 4; beta 0 of 4; focus 0 of 1; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |
| 5S4U | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 4ZH4 | default 3 of 3; deposited 3 of 3; beta 3 of 3; focus 3 of 3; domains 3 of 3; seed1 3 of 3; roundtrip 3 of 3 |
| 9Y70 | default 0 of 4; deposited 0 of 4; beta 0 of 4; focus 0 of 1; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |
| 9BZD | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 0 of 0; domains 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 9NMN | default 4 of 9; deposited 4 of 9; beta 4 of 9; focus 16 of 36; domains 4 of 9; seed1 4 of 9; roundtrip 4 of 9 |
| 8FZN | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 120; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 7U9Z | default 2 of 6; deposited 2 of 6; beta 2 of 6; focus 8 of 24; domains 2 of 6; seed1 2 of 6; roundtrip 2 of 6 |
| 9HEO | default 1 of 2; deposited 1 of 2; beta 1 of 2; focus 4 of 8; domains 1 of 2; seed1 1 of 2; roundtrip 1 of 2 |
| 7ADZ | default 1 of 2; deposited 1 of 2; beta 1 of 2; focus 6 of 30; domains 1 of 2; seed1 1 of 2; roundtrip 1 of 2 |
| 9FZL | default 4 of 6; deposited 4 of 6; beta 4 of 6; focus 1 of 2; domains 4 of 6; seed1 4 of 6; roundtrip 4 of 6 |
| 6DQN | default 3 of 6; deposited 3 of 6; beta 3 of 6; focus 12 of 24; domains 3 of 6; membrane 3 of 6; seed1 3 of 6; roundtrip 3 of 6 |
| 6JCY | default 3 of 3; deposited 3 of 3; beta 3 of 3; focus 2 of 2; domains 3 of 3; seed1 3 of 3; roundtrip 3 of 3 |
| 3G6E | default 1 of 2; deposited 1 of 2; beta 1 of 2; focus 1 of 2; domains 1 of 2; seed1 1 of 2; roundtrip 1 of 2 |
| 4AS1 | default 1 of 2; deposited 1 of 2; beta 1 of 2; focus 1 of 2; seed1 1 of 2; roundtrip 1 of 2 |
| 7F5S | default 0 of 3; deposited 0 of 3; beta 0 of 3; focus 0 of 3; domains 0 of 3; seed1 0 of 3; roundtrip 0 of 3 |
| 6QIK | default 0 of 3; deposited 0 of 3; beta 0 of 3; focus 0 of 3; domains 0 of 3; seed1 0 of 3; roundtrip 0 of 3 |
| AF_AFQ8ZPJ0F1 | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |

## Copies' sheets that differ, left out of A15 (ruling C2b-3/R2: reported, not judged)

| id | view | tie | sheets left out | judged |
|---|---|---|---|---|
| 3EOM | default | D~A | copies' sheets differ: D sheet 9 vs A sheets 0+1; 14 strands not judged | yes, 18 matched |
| 5C1A | focus | F~A | copies' sheets differ: F sheet 29 vs A sheet 0; 12 strands not judged | yes, 50 matched |
| 8U1N | default | B~A | copies' sheets differ: B sheet 9 vs A sheet 6; 6 strands not judged | yes, 28 matched |

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
