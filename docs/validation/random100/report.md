# foldmap random-100 campaign

91 of 100 pass.

Report written by: foldmap 0.2.0, git e1392f79daf7, code fingerprint 7879eea05bd9301c, scripts/random100.py sha256 b552d3c3d4b3; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Entries measured by 1 code version(s): 7879eea05bd9301c+b552d3c3d4b358ca (100)

An entry passes when no criterion fails. A criterion with nothing to judge in an entry is not judged (or does not apply), never counted as a pass: see the outcomes per criterion below.

Sample: seed 20261013, release-date cutoff 2026-10-13.

| stratum | entries | pass |
|---|---|---|
| S1 | 12 | 12 |
| S2 | 12 | 12 |
| S3 | 10 | 10 |
| S4 | 15 | 14 |
| S5 | 8 | 8 |
| S6 | 8 | 7 |
| S7 | 7 | 4 |
| S8 | 7 | 5 |
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
| A2b the app draws it within the hosted limits | 100 | 0 | 0 | 0 |
| A3 no route fallbacks | 98 | 2 | 0 | 0 |
| A4 no loop over a label | 100 | 0 | 0 | 0 |
| A5 no element or stem boxes overlapping | 100 | 0 | 0 | 0 |
| A6 no text collisions | 97 | 3 | 0 | 0 |
| A7 text and legend inside the page | 100 | 0 | 0 | 0 |
| A8 labels unique per chain, stem names per figure | 99 | 0 | 1 | 0 |
| A9 nothing left out without a note | 100 | 0 | 0 | 0 |
| A10 secondary structure agrees with mkdssp | 96 | 0 | 0 | 4 |
| A11 N->C and sheet order | 95 | 0 | 1 | 4 |
| A12 same figure for any PYTHONHASHSEED | 100 | 0 | 0 | 0 |
| A13 a saved layout redraws the same figure | 100 | 0 | 0 | 0 |
| A14 interactive page and sequence plot written | 100 | 0 | 0 | 0 |
| A15 copies drawn alike | 37 | 0 | 11 | 52 |
| A16 domains: one chain, alike, not tiny, panels apart, names clear | 48 | 0 | 14 | 38 |
| A17 nucleic-acid metrics and stem names by their letters | 12 | 5 | 1 | 82 |
| A18 membrane: only the band set in the band (TM helices across it, slab helices at their depths, crossing strands), no element nearer than its residues or on the wrong side, no loop nearer than its residues but out of a band element's port to its face or, for a loop in the bilayer, beside its own port to its residues' depth; a loop joining the faces crosses the band only beside its own port or past the band's ends | 8 | 0 | 0 | 92 |
| A19 two states: same names in the same order | 8 | 2 | 0 | 90 |
| A20 ligand tethers local | 10 | 0 | 51 | 39 |

### Not judged, per criterion

- **A8** (1): 6B5W (no element or stem drawn)
- **A11** (1): 6B5W (no protein element drawn)
- **A15** (11): 5SB5 (copies drawn with no symmetry tie: default C~A, D~B); 6DPV (copies drawn with no symmetry tie: default C~A, D~B, E~A, F~B, G~B, H~B, I~B, J~A, K~A, L~A); 5JB1 (copies drawn with no symmetry tie: default B~A, C~A, D~A, E~A, F~A); 7JU4 (copies drawn with no symmetry tie: default 4~0, 8~6, 9~7, C~6, D~7, G~6, H~7, I~6, J~7, K~6, L~7, N~M, O~6, P~7, Q~6, R~7, S~6, T~7, W~6, X~7, Y~6, Z~7, b~a, c~); 8ESZ (copies drawn with no symmetry tie: default AB~AC); 7ARD (copies drawn with no symmetry tie: default u~b); 7SMR (copies drawn with no symmetry tie: default D~A); 6JO6 (copies drawn with no symmetry tie: default Z~1); 5LDW (copies drawn with no symmetry tie: default U~T, l~Q, r~i); 8OZJ (copies drawn with no symmetry tie: default H~D, L~D, E~A, I~A, G~C, K~C); 6YQD (default BBB~AAA n/a: 2 matched)
- **A16** (14): 7IMM (the domains view drew no domain); 7FDN (the domains view drew no domain); 6YF9 (the domains view drew no domain); 1SLV (the domains view drew no domain); 1BXI (the domains view drew no domain); 4UM3 (the domains view drew no domain); 8R5X (the domains view drew no domain); 6AD0 (the domains view drew no domain); 2MJW (the domains view drew no domain); 8W9Y (the domains view drew no domain); 8BW5 (the domains view drew no domain); 6YPX (the domains view drew no domain); 6YQD (the domains view drew no domain); 8P3N (the domains view drew no domain)
- **A17** (1): 7MW8 (no view draws the nucleic acid)
- **A20** (51): 13MV (no tether to another chain); 5Q42 (no tether to another chain); 5HI7 (no tether to another chain); 2YBE (no tether to another chain; reported: ligand-tether:3PG:A1421:1 7.644 (home A); ligand-tether:ALF:A1420:1 9.371 (home A); ligand-tether:LA8:A1419:3 12.74 (home ); 6IDJ (no tether to another chain; reported: ligand-tether:9AU:A406:3 11.81 (home A); ligand-tether:FMN:A404:0 9.631 (home A); ligand-tether:FMN:A404:1 8.066 (home A);); 9C0M (no tether to another chain); 7LX9 (no tether to another chain); 4XG1 (no tether to another chain; reported: ligand-tether:LLP:A508:1 8.779 (home A); ligand-tether:LLP:A508:3 16.432 (home A); ligand-tether:LLP:B508:3 12.858 (home B); 2C4M (no tether to another chain; reported: ligand-tether:PLP:A1634:1 11.42 (home A); ligand-tether:PLP:B1634:1 8.015 (home B); ligand-tether:PLP:B1634:2 7.715 (home ); 6KTL (no tether to another chain; reported: ligand-tether:INS:D402:0 7.349 (home D)); 7OKA (no tether to another chain); 2ILU (no tether to another chain; reported: ligand-tether:NDP:A4482:0 6.603 (home A4)); 7GLY (no tether to another chain; reported: ligand-tether:RDX:A407:3 8.009 (home A); ligand-tether:RDX:B405:3 7.947 (home B)); 7SFH (no tether to another chain; reported: ligand-tether:91I:A1401:3 7.595 (home A1); ligand-tether:91I:A2401:3 7.59 (home A2)); 3L9G (no tether to another chain; reported: ligand-tether:MRD:A1850:2 8.072 (home A1); ligand-tether:URC:A2300:0 14.581 (home A2); ligand-tether:URC:A2300:1 9.621 (ho); 5SB5 (no tether to another chain; reported: ligand-tether:GTP:A501:1 8.179 (home A); ligand-tether:GTP:A501:3 11.837 (home A); ligand-tether:GDP:B501:1 6.905 (home B)); 5DWZ (no tether to another chain; reported: ligand-tether:MPO:C406:0 7.113 (home C); ligand-tether:MPO:C406:1 7.113 (home C)); 1A0Z (no tether to another chain); 1BXI (no tether to another chain); 5X4Z (no tether to another chain); 4UM3 (no tether to another chain; reported: ligand-tether:09R:T301:2 7.762 (home T); ligand-tether:09R:Q301:2 7.762 (home Q); ligand-tether:09R:S301:2 7.762 (home S)); 6DPV (no tether to another chain; reported: ligand-tether:GTP:A501:3 10.062 (home A); ligand-tether:GDP:B501:3 16.068 (home B); ligand-tether:GTP:K501:3 9.767 (home K); 8PQV (no tether to another chain; reported: ligand-tether:ADP:A4701:3 7.463 (home A); ligand-tether:ATP:A4702:0 8.758 (home A); ligand-tether:ATP:A4702:1 7.845 (home ); 6AD0 (no tether to another chain); 9ON5 (no tether to another chain); 7AOD (no tether to another chain); 7QUP (no tether to another chain; reported: ligand-tether:GTP:7A500:3 10.232 (home 7A)); 8RRX (no tether to another chain; reported: ligand-tether:ATP:D5103:1 6.841 (home D); ligand-tether:CFF:D5104:0 7.01 (home D); ligand-tether:CFF:D5104:2 8.007 (home D); 7SMR (no tether to another chain; reported: ligand-tether:CLR:A602:0 11.57 (home A); ligand-tether:CLR:A602:1 8.692 (home A); ligand-tether:CLR:C503:2 9.12 (home C); ); 5JDN (no tether to another chain); 6PS3 (no tether to another chain); 7PIU (no tether to another chain); 5LDW (no tether to another chain; reported: ligand-tether:FMN:F501:1 8.116 (home F); ligand-tether:SF4:F502:3 17.97 (home F)); 9GEP (no tether to another chain; reported: ligand-tether:HEM:N603:3 7.367 (home N); ligand-tether:NAG:N602:1 9.928 (home N)); 3L2V (no tether to another chain); 8BW5 (no tether to another chain); 3MHT (no tether to another chain); 9YLM (no tether to another chain); 7XSP (no tether to another chain); 7MW8 (no tether to another chain); 8VWQ (no tether to another chain); 8SA2 (no tether to another chain; reported: ligand-tether:B1Z:A301:0 6.645 (home A); ligand-tether:B1Z:A301:1 13.548 (home A); ligand-tether:B1Z:A301:2 9.75 (home A);); 3OWW (no tether to another chain); 8XZE (no tether to another chain); 8K8E (no tether to another chain; reported: ligand-tether:CLR:C301:1 8.241 (home C)); 1RM5 (no tether to another chain); 6YPX (no tether to another chain); 8P3N (no tether to another chain; reported: ligand-tether:ADZ:B201:0 7.167 (home B); ligand-tether:ADZ:B201:3 8.085 (home B); ligand-tether:ADZ:A201:1 6.661 (home A)); 8OYQ (no tether to another chain); 6FYF (no tether to another chain); 6FY9 (no tether to another chain)

## Entries

process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 gates; drawing s: the drawing inside it.

| id | stratum | bin | chains | process s | drawing s | A1 | A1c | A2 | A2b | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 | A11 | A12 | A13 | A14 | A15 | A16 | A17 | A18 | A19 | A20 | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 13MV | S1 | M | 1 | 1.35 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5Q42 | S1 | M | 1 | 1.1 | 0.37 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5HI7 | S1 | M | 1 | 1.3 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2YBE | S1 | M | 1 | 1.45 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 29RY | S1 | M | 1 | 1.2 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6IDJ | S1 | M | 1 | 1.45 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3RVC | S1 | S | 1 | 0.9 | 0.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 1F09 | S1 | S | 1 | 1.05 | 0.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9C0M | S1 | S | 1 | 1.35 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 7IMM | S1 | S | 1 | 1.55 | 0.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ✓ |  |
| 9B7C | S1 | S | 1 | 1.15 | 0.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 7LX9 | S1 | S | 1 | 1.25 | 0.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 4XG1 | S2 | L | 4 | 2.71 | 1.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2C4M | S2 | L | 4 | 3.81 | 2.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6KTL | S2 | L | 4 | 2.71 | 1.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7OKA | S2 | L | 30 | 2.31 | 1.18 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2ILU | S2 | M | 1 | 4.17 | 2.94 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7GLY | S2 | M | 2 | 1.75 | 0.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7FDN | S2 | M | 2 | 1.4 | 0.53 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 7SFH | S2 | M | 1 | 1.85 | 0.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9PXS | S2 | S | 1 | 1.81 | 0.72 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 3L9G | S2 | S | 1 | 3.68 | 2.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2B9V | S2 | XL | 16 | 2.26 | 1.08 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - |  |
| 6YF9 | S2 | XL | 180 | 7.68 | 5.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5SB5 | S3 | L | 6 | 3.97 | 2.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 4XEI | S3 | L | 7 | 2.62 | 1.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 5DWZ | S3 | L | 8 | 3.01 | 1.83 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2I5Y | S3 | L | 8 | 5.48 | 4.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 1SLV | S3 | M | 2 | 3.37 | 2.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ✓ |  |
| 2QNA | S3 | M | 2 | 1.76 | 0.79 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 1A0Z | S3 | M | 4 | 1.66 | 0.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 1BXI | S3 | S | 2 | 1.11 | 0.29 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 5X4Z | S3 | XL | 24 | 5.91 | 3.91 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 4UM3 | S3 | XL | 40 | 1.95 | 1.01 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 6DPV | S4 | L | 12 | 6.86 | 4.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ○ |  |
| 8R5X | S4 | L | 4 | 1.5 | 0.68 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | - |  |
| 6MHS | S4 | L | 4 | 3.21 | 2.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 5JB1 | S4 | L | 6 | 5.31 | 4.0 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | - |  |
| 8PQV | S4 | L | 1 | 4.51 | 2.82 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6AD0 | S4 | M | 6 | 1.6 | 0.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 9ON5 | S4 | M | 3 | 3.31 | 2.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7JU4 | S4 | XL | 61 | 58.52 | 51.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 7AOD | S4 | XL | 24 | 12.63 | 9.02 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7QUP | S4 | XL | 65 | 7.37 | 4.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9P1X | S4 | XL | 30 | 3.66 | 1.95 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ |  |
| 8RRX | S4 | XL | 12 | 12.15 | 8.77 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ | A6: default 1 (label:G:4806-4819 × helix:G:4806-4819); deposited 1 (label:G:4806-4819 × helix:G:4806-4819); beta 1 (label:G:4806-4819 × helix:G:4806-4819); domains 1 (label:G:4806-4819 × helix:G:4806-4819) |
| 8RVL | S4 | XL | 34 | 8.83 | 5.96 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 8ESZ | S4 | XL | 43 | 29.5 | 24.89 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 7ARD | S4 | XL | 51 | 23.28 | 18.65 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 2MJW | S5 | M | 4 | 1.25 | 0.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 2DGS | S5 | S | 1 | 1.15 | 0.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8Q4O | S5 | S | 1 | 2.56 | 0.13 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 1MG8 | S5 | S | 1 | 1.13 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 5KQJ | S5 | S | 1 | 1.26 | 0.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8OH7 | S5 | S | 1 | 1.31 | 0.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6B5W | S5 | S | 1 | 1.02 | 0.06 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | ✓ | ○ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8B4R | S5 | S | 1 | 1.02 | 0.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6JLN | S6 | L | 39 | 31.98 | 29.02 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ |  |
| 7SMR | S6 | L | 5 | 6.13 | 4.89 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✓ | - | ○ |  |
| 6JO6 | S6 | L | 19 | 21.11 | 13.65 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✓ | - | ✓ |  |
| 8W9Y | S6 | L | 6 | 4.22 | 3.05 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | ✓ | - | - |  |
| 5JDN | S6 | M | 1 | 1.25 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 6PS3 | S6 | M | 1 | 1.25 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 7PIU | S6 | M | 6 | 2.96 | 1.95 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 5LDW | S6 | XL | 45 | 29.57 | 26.73 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | ✓ | - | ○ | A6: deposited 6 (label:G:265-268 × helix:G:265-268; label:G:296-308 × helix:G:296-308; label:G:389-401 × helix:G:389-401) |
| 9K9V | S7 | L | 5 | 5.01 | 1.51 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | - |  |
| 9GEP | S7 | L | 12 | 9.02 | 5.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ |  |
| 3L2V | S7 | M | 4 | 3.86 | 1.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ |  |
| 8BW5 | S7 | M | 3 | 3.61 | 1.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | ✗ | - | - | ○ | A17: pair_f1 0.8 (limit 0.95); stem_recall 0.5 (limit 0.85) |
| 3MHT | S7 | M | 3 | 3.41 | 0.59 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 9YLM | S7 | S | 3 | 2.51 | 0.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | ○ | A17: stems_nearer_frac 0.0 (limit 1.0) |
| 8QEM | S7 | XL | 26 | 2.56 | 1.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | - | - | - | A17: measured on the focus view, as the default draws no nucleic acid: angle_median 3.2 (limit 0.05); angle_max 4.6 (limit 0.05); joint_error_max 3.4 (limit 0.05) |
| 7XSP | S8 | L | 3 | 4.99 | 1.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 7MW8 | S8 | L | 12 | 1.77 | 0.6 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ | - | - | ○ |  |
| 9C75 | S8 | M | 3 | 5.35 | 2.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | - |  |
| 9CJI | S8 | M | 4 | 5.39 | 1.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 6SO9 | S8 | S | 2 | 3.96 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 7Y7C | S8 | XL | 55 | 72.86 | 53.34 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ✓ | A3: domains: strand a:2705-2705, gap after a:1439; A6: domains 1 (na:stem:228:first × na:run:443); A17: na_crossings 2 (limit 0) |
| 8VWQ | S8 | XL | 56 | 76.85 | 58.42 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ○ | A3: domains: strand 1:836-836; A17: na_crossings 1 (limit 0) |
| 8SA2 | S9 | M | 2 | 5.41 | 1.71 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 3OWW | S9 | S | 2 | 4.16 | 1.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 8XZE | S9 | S | 1 | 2.81 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 8K8E | S10 | L | 5 | 5.61 | 4.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 7BBH | S10 | L | 3 | 4.66 | 3.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 1WWL | S10 | M | 2 | 2.91 | 1.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 7BRG | S10 | M | 3 | 2.81 | 1.71 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 8OZJ | S10 | XL | 11 | 19.79 | 18.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ✓ | - | - | - | ✓ |  |
| 1RM5 | S11 | M | 3 | 3.06 | 1.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ○ |  |
| 2HKI | S11 | M | 1 | 3.06 | 1.78 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - |  |
| 9YKJ | S11 | L | 2 | 2.3 | 1.09 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | ✗ | - | A19: names compared on all 43 pairs, 27 differ; αN~α2 (A:37-47~A:37-50); αN′~α3 (A:72-85~A:72-85); α1~α4 (A:100-108~A:98-111); α2~α7 (A:164-175~A:164-174); α3~α8 (A:184-192~A:184-198); α3″~α9 (A:204-210~A:203-210); α4~α10 (A:221-226~A:221-226); α4′~α11 (A:228-236~A:228-237) (+49 more) |
| 10UJ | S11 | M | 1 | 2.11 | 0.99 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✗ | - | A19: names compared on all 43 pairs, 27 differ; α2~αN (A:37-50~A:37-47); α3~αN′ (A:72-85~A:72-85); α4~α1 (A:98-111~A:100-108); α7~α2 (A:164-174~A:164-175); α8~α3 (A:184-198~A:184-192); α9~α3″ (A:203-210~A:204-210); α10~α4 (A:221-226~A:221-226); α11~α4′ (A:228-237~A:228-236) (+49 more) |
| 6YPX | S11 | M | 4 | 1.4 | 0.47 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | ✓ | ○ |  |
| 6YQD | S11 | M | 2 | 1.5 | 0.6 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | - | - | ✓ | - |  |
| 8P3N | S11 | M | 2 | 1.75 | 0.67 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | ✓ | ○ |  |
| 8OYQ | S11 | M | 2 | 1.2 | 0.29 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 6FYF | S11 | L | 4 | 1.3 | 0.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 6FY9 | S11 | M | 1 | 1.4 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| AF_AFA0A4E9ES35F1 | S12 | M | 1 | 1.35 | 0.47 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFO42816F1 | S12 | M | 1 | 1.6 | 0.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFB1Q234F1 | S12 | S | 1 | 0.85 | 0.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |

## Two-state pairs (S11)

| pair | seed | partner | sampled as | types now | ligand touching | assembly 1 chains | CA RMSD Å | A19 |
|---|---|---|---|---|---|---|---|---|
| P1 | 1RM5 | 2HKI | apo/holo | apo/holo | [True, False] | [4, 4] | 0.785 | ✓ names compared on all 27 pairs; A1:79-81 unpaired; A1:102-105 unpaired; A1:107-110 unpaired; A1:108-111 unpaired |
| P2 | 9YKJ | 10UJ | monomer/oligomer | monomer/oligomer | [False, False] | [2, 1] | 2.233 | ✗ names compared on all 43 pairs, 27 differ; αN~α2 (A:37-47~A:37-50); αN′~α3 (A:72-85~A:72-85); α1~α4 (A:100-108~A:98-111); α2~α7 (A:164-175~A:164-174); α3~α8 (A:184-192~A:184-198); α3″~α9 (A:204-210~A:203-210); α4~α10 (A:221-226~A:221-226); α4′~α11 (A:228-236~A:228-237) (+49 more) |
| P3 | 6YPX | 6YQD | apo/holo | none | [False, False] | [2, 2] | 2.539 | ✓ names compared on all 7 pairs; AAA:55-60 unpaired; AAA:65-67 unpaired |
| P4 | 8P3N | 8OYQ | monomer/oligomer | monomer/oligomer | [True, True] | [2, 1] | 0.58 | ✓ names compared on all 12 pairs; B:121-125 unpaired; B:120-123 unpaired |
| P5 | 6FYF | 6FY9 | open/closed | none | [True, True] | [1, 1] | 0.627 | ✓ names compared on all 34 pairs; A:261-263 unpaired; A:363-365 unpaired; A:412-413 unpaired; A:103-105 unpaired; A:195-197 unpaired |

## Reported, not gated

| id | honest_crossings | loop_crossings | na_run_crossings | fallback_curves | fallback_curves_crossing_an_element | kink_links_rejected | tighter_than_inflate | rejected_as_crossing | names_moved | area | max_rss_mb |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 13MV | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 435.02 | 169.4 |
| 5Q42 | 37 | 37 | None | 0 | 0 | 0 | 0 | 0 | 0 | 424.94 | 168.4 |
| 5HI7 | 15 | 15 | None | 0 | 0 | 0 | 0 | 0 | 0 | 592.21 | 173.1 |
| 2YBE | 35 | 35 | None | 0 | 0 | 0 | 1 | 0 | 0 | 521.66 | 170.1 |
| 29RY | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 464.93 | 170.4 |
| 6IDJ | 17 | 17 | None | 0 | 0 | 0 | 0 | 0 | 0 | 511.73 | 169.5 |
| 3RVC | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 167.63 | 165.7 |
| 1F09 | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 252.21 | 167.1 |
| 9C0M | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 331.75 | 167.9 |
| 7IMM | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 334.14 | 182.8 |
| 9B7C | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 106.97 | 166.6 |
| 7LX9 | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 159.92 | 166.7 |
| 4XG1 | 138 | 138 | None | 0 | 0 | 0 | 2 | 0 | 0 | 1421.95 | 228.3 |
| 2C4M | 132 | 132 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2453.41 | 300.3 |
| 6KTL | 96 | 96 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1785.27 | 222.0 |
| 7OKA | 223 | 223 | None | 0 | 0 | 0 | 0 | 0 | 0 | 799.41 | 307.5 |
| 2ILU | 311 | 311 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3320.71 | 323.8 |
| 7GLY | 24 | 24 | None | 0 | 0 | 0 | 0 | 0 | 0 | 730.43 | 198.6 |
| 7FDN | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 357.62 | 179.3 |
| 7SFH | 28 | 28 | None | 0 | 0 | 0 | 0 | 0 | 0 | 745.45 | 204.8 |
| 9PXS | 34 | 34 | None | 0 | 0 | 0 | 0 | 0 | 0 | 572.75 | 200.4 |
| 3L9G | 213 | 213 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1629.75 | 217.4 |
| 2B9V | 61 | 61 | None | 0 | 0 | 0 | 0 | 0 | 0 | 770.37 | 276.0 |
| 6YF9 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 138.13 | 586.5 |
| 5SB5 | 129 | 129 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3663.2 | 362.8 |
| 4XEI | 44 | 44 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2585.62 | 270.4 |
| 5DWZ | 137 | 137 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1541.88 | 218.3 |
| 2I5Y | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1336.58 | 208.5 |
| 1SLV | 69 | 69 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1200.7 | 183.1 |
| 2QNA | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1580.29 | 210.0 |
| 1A0Z | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 735.75 | 183.2 |
| 1BXI | 8 | 8 | None | 0 | 0 | 0 | 1 | 0 | 0 | 279.91 | 167.9 |
| 5X4Z | 310 | 310 | None | 0 | 0 | 0 | 1 | 0 | 0 | 4803.34 | 435.6 |
| 4UM3 | 35 | 35 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1190.0 | 307.3 |
| 6DPV | 427 | 427 | None | 0 | 0 | 0 | 5 | 0 | 0 | 7264.62 | 517.7 |
| 8R5X | 53 | 53 | None | 0 | 0 | 0 | 0 | 0 | 0 | 946.98 | 173.6 |
| 6MHS | 37 | 37 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3381.91 | 332.1 |
| 5JB1 | 486 | 486 | None | 0 | 0 | 0 | 2 | 0 | 0 | 3335.07 | 475.2 |
| 8PQV | 201 | 201 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4188.14 | 256.2 |
| 6AD0 | 66 | 66 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1133.97 | 176.1 |
| 9ON5 | 127 | 127 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1605.1 | 219.6 |
| 7JU4 | 2571 | 2571 | None | 0 | 0 | 0 | 6 | 0 | 0 | 29826.53 | 618.4 |
| 7AOD | 499 | 499 | None | 0 | 0 | 0 | 0 | 0 | 0 | 8034.48 | 520.2 |
| 7QUP | 25 | 25 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1077.17 | 713.9 |
| 9P1X | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1693.39 | 420.5 |
| 8RRX | 300 | 300 | None | 0 | 0 | 0 | 0 | 0 | 0 | 9410.73 | 506.2 |
| 8RVL | 483 | 483 | None | 0 | 0 | 0 | 2 | 0 | 0 | 11585.47 | 521.9 |
| 8ESZ | 279 | 279 | None | 0 | 0 | 0 | 4 | 0 | 0 | 14940.61 | 481.0 |
| 7ARD | 598 | 598 | None | 0 | 0 | 0 | 1 | 0 | 0 | 19197.26 | 512.7 |
| 2MJW | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 402.21 | 289.5 |
| 2DGS | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 49.82 | 195.2 |
| 8Q4O | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 57.85 | 180.9 |
| 1MG8 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 36.96 | 166.8 |
| 5KQJ | 5 | 5 | None | 0 | 0 | 0 | 0 | 0 | 0 | 384.61 | 213.7 |
| 8OH7 | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 429.12 | 184.5 |
| 6B5W | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 163.0 |
| 8B4R | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 12.71 | 165.1 |
| 6JLN | 52 | 52 | None | 0 | 0 | 0 | 6 | 0 | 0 | 4253.35 | 481.6 |
| 7SMR | 85 | 85 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3148.74 | 343.6 |
| 6JO6 | 180 | 180 | None | 0 | 0 | 0 | 2 | 0 | 0 | 6463.28 | 421.2 |
| 8W9Y | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1327.36 | 213.6 |
| 5JDN | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 317.04 | 168.8 |
| 6PS3 | 9 | 9 | None | 0 | 0 | 0 | 1 | 0 | 0 | 612.53 | 171.5 |
| 7PIU | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1324.58 | 201.3 |
| 5LDW | 145 | 145 | None | 0 | 0 | 0 | 0 | 0 | 0 | 10373.41 | 461.3 |
| 9K9V | 127 | 127 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2388.63 | 280.5 |
| 9GEP | 35 | 35 | 3 | 0 | 0 | 0 | 2 | 0 | 0 | 4094.72 | 277.9 |
| 3L2V | 51 | 51 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2143.77 | 222.5 |
| 8BW5 | 8 | 8 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 720.85 | 181.3 |
| 3MHT | 18 | 18 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 396.7 | 179.0 |
| 9YLM | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 155.92 | 172.5 |
| 8QEM | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1473.25 | 387.9 |
| 7XSP | 101 | 101 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1757.19 | 200.6 |
| 7MW8 | 34 | 34 | None | 0 | 0 | 0 | 2 | 0 | 0 | 490.56 | 172.5 |
| 9C75 | 25 | 25 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 1141.8 | 185.7 |
| 9CJI | 56 | 56 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2034.39 | 203.2 |
| 6SO9 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 88.68 | 249.1 |
| 7Y7C | 335 | 335 | 722 | 0 | 0 | 0 | 1 | 0 | 0 | 31653.95 | 637.4 |
| 8VWQ | 443 | 443 | 638 | 0 | 0 | 0 | 1 | 0 | 0 | 32645.44 | 675.1 |
| 8SA2 | 0 | 0 | 13 | 0 | 0 | 0 | 0 | 0 | 0 | 1973.64 | 200.4 |
| 3OWW | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 695.11 | 186.7 |
| 8XZE | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 122.15 | 172.2 |
| 8K8E | 99 | 99 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2311.75 | 244.9 |
| 7BBH | 52 | 52 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1415.51 | 260.9 |
| 1WWL | 174 | 174 | None | 0 | 0 | 0 | 0 | 0 | 0 | 514.28 | 217.6 |
| 7BRG | 52 | 52 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1078.11 | 219.5 |
| 8OZJ | 242 | 242 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4678.64 | 440.9 |
| 1RM5 | 97 | 97 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1541.86 | 226.6 |
| 2HKI | 87 | 87 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1526.49 | 215.0 |
| 9YKJ | 85 | 85 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1707.13 | 245.6 |
| 10UJ | 83 | 83 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1232.29 | 174.9 |
| 6YPX | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 209.8 | 178.6 |
| 6YQD | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 229.54 | 178.2 |
| 8P3N | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 364.39 | 189.3 |
| 8OYQ | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 154.1 | 168.8 |
| 6FYF | 36 | 36 | None | 0 | 0 | 0 | 0 | 0 | 0 | 466.19 | 172.1 |
| 6FY9 | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 593.76 | 170.2 |
| AF_AFA0A4E9ES35F1 | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 555.38 | 172.2 |
| AF_AFO42816F1 | 27 | 27 | None | 0 | 0 | 0 | 0 | 0 | 0 | 820.0 | 171.5 |
| AF_AFB1Q234F1 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 176.37 | 164.6 |

## Fallbacks taken

Each named fallback (foldmap.fallback.FALLBACKS, docs/validation/inventory/inventory.md) counted as the figures were made: reported, not gated.

| fallback | entries (default view) | times (default view) | times (every view) |
|---|---|---|---|
| loop_curve | 0 | 0 | 0 |
| loop_over_text | 0 | 0 | 0 |
| loop_wall_detour | 0 | 0 | 0 |
| loop_wall_unkept | 0 | 0 | 0 |
| route_unsnapped | 0 | 0 | 0 |
| route_second_search_refused | 0 | 0 | 0 |
| helix_name_moved | 0 | 0 | 49 |
| na_tail_straight | 0 | 0 | 0 |
| na_gap_unclear | 0 | 0 | 1 |
| na_hairpin_unclear | 0 | 0 | 0 |
| na_run_unclear | 0 | 0 | 2 |
| na_trace_unclear | 0 | 0 | 0 |
| stem_name_least_bad | 0 | 0 | 0 |
| domain_name_least_bad | 0 | 0 | 0 |
| glyph_width_guessed | 2 | 176 | 1232 |
| mark_least_bad_spot | 0 | 0 | 0 |
| counts_hidden | 8 | 13 | 83 |
| interactions_dropped | 0 | 0 | 0 |
| interactions_hidden | 0 | 0 | 0 |
| crosslinks_unplaced | 0 | 0 | 0 |
| page_text_small | 1 | 1 | 13 |
| s3_cycles_capped | 0 | 0 | 0 |

| id | default view |
|---|---|
| 2B9V | counts_hidden 1 |
| 7JU4 | page_text_small 1 |
| 8RRX | counts_hidden 1 |
| 8RVL | counts_hidden 1 |
| 8ESZ | counts_hidden 1 |
| 7SMR | counts_hidden 1 |
| 6JO6 | counts_hidden 1 |
| 7Y7C | counts_hidden 4; glyph_width_guessed 100 |
| 8VWQ | counts_hidden 3; glyph_width_guessed 76 |

## Science checks (S1-S6), apart from A1-A20

Judged in every drawn view (default, deposited, beta, focus, domains, membrane, app); S5 once per entry.

Check kinds: S1, S2, S3, S4, S6 structural consistency (the figure against its own structure, internal); S5 external validation (against CATH); A1-A20 readability (the page).

A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet drawn, no barrel, no residue mark, no accessory strand, no CATH chopping) is not judged, one with no protein does not apply.

| check | pass | fail | not judged | does not apply | kind |
|---|---|---|---|---|---|
| S1 sheet neighbours pair | 90 | 0 | 6 | 4 | consistency |
| S2 barrel closure pairs | 13 | 0 | 83 | 4 | consistency |
| S3 barrel strand count is its pairing cycle | 90 | 0 | 6 | 4 | consistency |
| S4 residue marks at their residue | 2 | 0 | 94 | 4 | consistency |
| S5 automatic domains agree with a reference | 12 | 18 | 66 | 4 | external |
| S6 accessory drawn at its partner | 37 | 0 | 59 | 4 | consistency |

### Science by view

Entries measured in each view, and how many of them fail each check there.

| view | measured | S1 | S2 | S3 | S4 | S6 |
|---|---|---|---|---|---|---|
| default | 100 | 0 | 0 | 0 | 0 | 0 |
| deposited | 96 | 0 | 0 | 0 | 0 | 0 |
| beta | 96 | 0 | 0 | 0 | 0 | 0 |
| focus | 67 | 0 | 0 | 0 | 0 | 0 |
| domains | 62 | 0 | 0 | 0 | 0 | 0 |
| membrane | 8 | 0 | 0 | 0 | 0 | 0 |
| app | 100 | 0 | 0 | 0 | 0 | 0 |

| id | S1 | S2 | S3 | S4 | S5 | S6 | findings |
|---|---|---|---|---|---|---|---|
| 13MV | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 5Q42 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.77, boundary F1 0.80, recall 1.00 (2 reference domains, 3 automatic) |
| 5HI7 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.78, boundary F1 0.89, recall 0.80 (5 reference domains, 4 automatic) |
| 2YBE | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.95, boundary F1 0.67, recall 1.00 (2 reference domains, 2 automatic) |
| 29RY | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6IDJ | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 3RVC | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 1F09 | ✓ | ✓ | ✓ | ✓ | ✓ | ○ |  |
| 9C0M | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7IMM | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9B7C | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7LX9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4XG1 | ✓ | ✓ | ✓ | ○ | ✓ | ✓ |  |
| 2C4M | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.62, boundary F1 0.57, recall 1.00 (2 reference domains, 5 automatic); B: overlap 0.62, boundary F1 0.57, recall 1.00 (2 reference domains, 5 automatic); C: overlap 0.62, boundary F1 0.57, recall 1.00 (2 reference domains, 5 automatic); D: overlap 0.63, boundary F1 0.57, recall 1.00 (2 reference domains, 5 automatic) |
| 6KTL | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7OKA | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2ILU | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7GLY | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7FDN | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7SFH | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9PXS | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3L9G | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.54, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 2B9V | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.76, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); B: overlap 0.76, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); C: overlap 0.76, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); D: overlap 0.74, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); E: overlap 0.76, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); F: overlap 0.76, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); G: overlap 0.74, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); H: overlap 0.76, |
| 6YF9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5SB5 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 4XEI | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.92, boundary F1 0.83, recall 0.83 (4 reference domains, 4 automatic); C: overlap 0.47, boundary F1 0.00, recall 1.00 (1 reference domains, 3 automatic); D: overlap 0.81, boundary F1 0.67, recall 1.00 (2 reference domains, 3 automatic) |
| 5DWZ | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 2I5Y | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 1SLV | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 B: overlap 0.54, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 2QNA | ○ | ○ | ○ | ○ | ✗ | ○ | S5 A: overlap 0.24, boundary F1 0.00, recall 1.00 (1 reference domains, 6 automatic) |
| 1A0Z | ○ | ○ | ○ | ○ | ✗ | ○ | S5 B: overlap 0.58, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); D: overlap 0.58, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 1BXI | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5X4Z | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 4UM3 | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 6DPV | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.72, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic); B: overlap 0.80, boundary F1 0.40, recall 0.50 (3 reference domains, 3 automatic); C: overlap 0.72, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic); D: overlap 0.80, boundary F1 0.40, recall 0.50 (3 reference domains, 3 automatic); E: overlap 0.72, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic); F: overlap 0.80, boundary F1 0.40, recall 0.50 (3 reference domains, 3 automatic); G: overlap 0.80, boundary F1 0.40, recall 0.50 (3 reference domains, 3 automatic); H: overlap 0.80, |
| 8R5X | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6MHS | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5JB1 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8PQV | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6AD0 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9ON5 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7JU4 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7AOD | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7QUP | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9P1X | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8RRX | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8RVL | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 8ESZ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7ARD | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 2MJW | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 2DGS | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 8Q4O | - | - | - | - | - | - |  |
| 1MG8 | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 5KQJ | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.66, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 8OH7 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 6B5W | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8B4R | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6JLN | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 7SMR | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6JO6 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8W9Y | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 5JDN | ○ | ○ | ○ | ○ | ✗ | ○ | S5 A: overlap 0.90, boundary F1 0.67, recall 0.67 (2 reference domains, 2 automatic) |
| 6PS3 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7PIU | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5LDW | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9K9V | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 9GEP | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3L2V | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.84, boundary F1 0.86, recall 0.75 (4 reference domains, 3 automatic) |
| 8BW5 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3MHT | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.92, boundary F1 0.50, recall 0.50 (2 reference domains, 2 automatic) |
| 9YLM | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8QEM | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 7XSP | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7MW8 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9C75 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9CJI | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6SO9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7Y7C | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8VWQ | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8SA2 | - | - | - | - | - | - |  |
| 3OWW | - | - | - | - | - | - |  |
| 8XZE | - | - | - | - | - | - |  |
| 8K8E | ✓ | ○ | ✓ | ✓ | ○ | ○ |  |
| 7BBH | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 1WWL | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.53, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic); B: overlap 0.54, boundary F1 0.00, recall 1.00 (1 reference domains, 2 automatic) |
| 7BRG | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8OZJ | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 1RM5 | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 2HKI | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 9YKJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 10UJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6YPX | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6YQD | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8P3N | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8OYQ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6FYF | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.75, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); C: overlap 0.75, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic); D: overlap 0.75, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic) |
| 6FY9 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.60, boundary F1 0.83, recall 1.00 (3 reference domains, 5 automatic) |
| AF_AFA0A4E9ES35F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFO42816F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFB1Q234F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |

## S5: automatic domains against CATH (reported, not gated)

Source: CATH domain boundaries (https://www.cathdb.info) through the PDBe API (EMBL-EBI), licensed under CC BY 4.0; Sillitoe I. et al. (2021) CATH: increased structural coverage of functional space. Nucleic Acids Res 49:D266-D273.

124 chains scored in 30 entries; 59 of them with at least two CATH domains. Every chain CATH chops is scored, one-domain chains too (a one-domain chain keeps boundary F1 1 only when the automatic domains leave it whole). 18 of the chains are in the calibration set (source "reference": in-sample, tests/data/domains/cath_reference.json, the 33 chains the automatic domains were calibrated on); the never-seen check is the fresh set (cath_fresh.json).

| id | chain | source | CATH domains | automatic | overlap | boundary F1 | recall | precision |
|---|---|---|---|---|---|---|---|---|
| 5Q42 | A | reference | 2 | 3 | 0.77 | 0.80 | 1.00 | 0.67 |
| 5HI7 | A | reference | 5 | 4 | 0.78 | 0.89 | 0.80 | 1.00 |
| 2YBE | A | reference | 2 | 2 | 0.95 | 0.67 | 1.00 | 0.50 |
| 3RVC | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1F09 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4XG1 | A | reference | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 4XG1 | B | CATH | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 4XG1 | C | CATH | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 4XG1 | D | CATH | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 2C4M | A | reference | 2 | 5 | 0.62 | 0.57 | 1.00 | 0.40 |
| 2C4M | B | CATH | 2 | 5 | 0.62 | 0.57 | 1.00 | 0.40 |
| 2C4M | C | CATH | 2 | 5 | 0.62 | 0.57 | 1.00 | 0.40 |
| 2C4M | D | CATH | 2 | 5 | 0.63 | 0.57 | 1.00 | 0.40 |
| 2ILU | A | reference | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3L9G | A | CATH | 1 | 2 | 0.54 | 0.00 | 1.00 | 0.00 |
| 2B9V | A | reference | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | B | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | C | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | D | CATH | 3 | 4 | 0.74 | 0.86 | 1.00 | 0.75 |
| 2B9V | E | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | F | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | G | CATH | 3 | 4 | 0.74 | 0.86 | 1.00 | 0.75 |
| 2B9V | H | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | I | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | J | CATH | 3 | 4 | 0.74 | 0.86 | 1.00 | 0.75 |
| 2B9V | K | CATH | 3 | 4 | 0.74 | 0.86 | 1.00 | 0.75 |
| 2B9V | L | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | M | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | N | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 2B9V | O | CATH | 3 | 4 | 0.74 | 0.86 | 1.00 | 0.75 |
| 2B9V | P | CATH | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 4XEI | A | reference | 4 | 4 | 0.92 | 0.83 | 0.83 | 0.83 |
| 4XEI | C | CATH | 1 | 3 | 0.47 | 0.00 | 1.00 | 0.00 |
| 4XEI | D | CATH | 2 | 3 | 0.81 | 0.67 | 1.00 | 0.50 |
| 4XEI | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4XEI | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4XEI | G | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5DWZ | C | reference | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5DWZ | D | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5DWZ | F | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5DWZ | H | CATH | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 1SLV | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1SLV | B | reference | 2 | 1 | 0.54 | 0.00 | 0.00 | 0.00 |
| 2QNA | A | CATH | 1 | 6 | 0.24 | 0.00 | 1.00 | 0.00 |
| 1A0Z | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1A0Z | B | CATH | 1 | 2 | 0.58 | 0.00 | 1.00 | 0.00 |
| 1A0Z | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1A0Z | D | CATH | 1 | 2 | 0.58 | 0.00 | 1.00 | 0.00 |
| 1BXI | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1BXI | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | E | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | F | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | G | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | H | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | I | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | J | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | K | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | L | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | M | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | N | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | O | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | P | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | Q | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | R | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | S | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | T | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | U | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | V | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | W | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | X | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | Y | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | Z | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | a | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | b | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | c | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | d | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | e | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | f | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | g | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | h | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | i | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | j | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | k | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | l | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | m | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 4UM3 | n | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6DPV | A | reference | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 6DPV | B | CATH | 3 | 3 | 0.80 | 0.40 | 0.50 | 0.33 |
| 6DPV | C | CATH | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 6DPV | D | CATH | 3 | 3 | 0.80 | 0.40 | 0.50 | 0.33 |
| 6DPV | E | CATH | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 6DPV | F | CATH | 3 | 3 | 0.80 | 0.40 | 0.50 | 0.33 |
| 6DPV | G | CATH | 3 | 3 | 0.80 | 0.40 | 0.50 | 0.33 |
| 6DPV | H | CATH | 3 | 3 | 0.80 | 0.40 | 0.50 | 0.33 |
| 6DPV | I | CATH | 3 | 3 | 0.80 | 0.40 | 0.50 | 0.33 |
| 6DPV | J | CATH | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 6DPV | K | CATH | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 6DPV | L | CATH | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 2MJW | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2MJW | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2MJW | C | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2MJW | D | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2DGS | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1MG8 | A | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5KQJ | A | CATH | 1 | 2 | 0.66 | 0.00 | 1.00 | 0.00 |
| 5JDN | A | reference | 2 | 2 | 0.90 | 0.67 | 0.67 | 0.67 |
| 3L2V | A | reference | 4 | 3 | 0.84 | 0.86 | 0.75 | 1.00 |
| 3L2V | B | CATH | 1 | 1 | 1.00 | 1.00 | 1.00 | 1.00 |
| 3MHT | A | reference | 2 | 2 | 0.92 | 0.50 | 0.50 | 0.50 |
| 1WWL | A | CATH | 1 | 2 | 0.53 | 0.00 | 1.00 | 0.00 |
| 1WWL | B | CATH | 1 | 2 | 0.54 | 0.00 | 1.00 | 0.00 |
| 1RM5 | A | reference | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1RM5 | B | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1RM5 | O | CATH | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2HKI | A | reference | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6FYF | A | reference | 3 | 4 | 0.75 | 0.86 | 1.00 | 0.75 |
| 6FYF | B | CATH | 3 | 3 | 0.98 | 1.00 | 1.00 | 1.00 |
| 6FYF | C | CATH | 3 | 4 | 0.75 | 0.86 | 1.00 | 0.75 |
| 6FYF | D | CATH | 3 | 4 | 0.75 | 0.86 | 1.00 | 0.75 |
| 6FY9 | A | reference | 3 | 5 | 0.60 | 0.83 | 1.00 | 0.71 |
| mean of 124 (all) | | | | | 0.88 | 0.82 | 0.93 | 0.79 |
| mean of 59 (two domains or more) | | | | | 0.81 | 0.76 | 0.86 | 0.70 |
| mean of 65 (one domain) | | | | | 0.94 | 0.88 | 1.00 | 0.88 |

Protein entries with no chain scored (S5 not judged), by reason:

- CATH chops none of its chains: 63
- not a PDB entry (CATH chops PDB entries): 3

## Accessories named, not drawn level beside their partner (reported, not gated)

| id | named of accessory strands, per view |
|---|---|
| 13MV | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 3RVC | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 4XG1 | default 2 of 6; deposited 2 of 6; beta 2 of 6; focus 2 of 6; domains 2 of 6; seed1 2 of 6; roundtrip 2 of 6; interactive 2 of 6 |
| 2C4M | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2; interactive 2 of 2 |
| 7GLY | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 1; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2; interactive 0 of 2 |
| 7SFH | default 0 of 2; deposited 0 of 2; beta 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2; interactive 0 of 2 |
| 9PXS | default 2 of 2; deposited 2 of 2; beta 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2; interactive 2 of 2 |
| 3L9G | default 8 of 16; deposited 8 of 16; beta 8 of 16; domains 8 of 16; seed1 8 of 16; roundtrip 8 of 16 |
| 2B9V | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 4; seed1 0 of 1; roundtrip 0 of 1 |
| 5SB5 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 4XEI | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 5DWZ | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2; interactive 2 of 2 |
| 2I5Y | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2 |
| 1SLV | default 0 of 6; deposited 0 of 6; beta 0 of 6; focus 0 of 1; domains 0 of 6; seed1 0 of 6; roundtrip 0 of 6; interactive 0 of 6 |
| 5X4Z | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 4UM3 | default 0 of 5; deposited 0 of 5; beta 0 of 5; focus 0 of 5; domains 0 of 5; seed1 0 of 5; roundtrip 0 of 5 |
| 8R5X | default 3 of 4; deposited 3 of 4; beta 3 of 4; focus 0 of 0; domains 3 of 4; seed1 3 of 4; roundtrip 3 of 4 |
| 5JB1 | default 2 of 12; deposited 2 of 12; beta 2 of 12; focus 0 of 2; domains 2 of 12; seed1 2 of 12; roundtrip 2 of 12 |
| 6AD0 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 9ON5 | default 6 of 15; deposited 6 of 15; beta 6 of 15; focus 2 of 5; domains 6 of 15; seed1 6 of 15; roundtrip 6 of 15 |
| 7JU4 | default 3 of 7; deposited 3 of 7; beta 3 of 7; focus 0 of 1; domains 3 of 7; seed1 3 of 7; roundtrip 3 of 7 |
| 7AOD | default 4 of 6; deposited 4 of 6; beta 4 of 6; focus 0 of 0; domains 4 of 6; seed1 4 of 6; roundtrip 4 of 6; interactive 4 of 6 |
| 9P1X | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 8 of 26; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 8RRX | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 8 of 12; domains 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
| 8RVL | default 0 of 6; deposited 0 of 6; beta 0 of 6; focus 0 of 0; domains 0 of 6; seed1 0 of 6; roundtrip 0 of 6; interactive 0 of 6 |
| 7ARD | default 2 of 8; deposited 2 of 8; beta 2 of 8; focus 0 of 0; domains 2 of 8; seed1 2 of 8; roundtrip 2 of 8 |
| 5KQJ | default 0 of 1; deposited 0 of 1; beta 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 6JLN | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 2 of 3; domains 2 of 3; membrane 1 of 1; seed1 1 of 1; roundtrip 1 of 1; interactive 1 of 1 |
| 7SMR | default 1 of 1; deposited 1 of 1; beta 1 of 1; focus 0 of 0; domains 1 of 1; membrane 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 9K9V | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 7XSP | default 1 of 2; deposited 1 of 2; beta 1 of 2; focus 1 of 2; seed1 1 of 2; roundtrip 1 of 2 |
| 7Y7C | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 8VWQ | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 7BBH | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 3; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 8OZJ | default 0 of 6; deposited 0 of 6; beta 0 of 6; focus 0 of 2; domains 0 of 6; seed1 0 of 6; roundtrip 0 of 6 |
| 1RM5 | default 0 of 4; deposited 0 of 4; beta 0 of 4; focus 0 of 4; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |
| 2HKI | default 0 of 4; deposited 0 of 4; beta 0 of 4; domains 0 of 4; seed1 0 of 4; roundtrip 0 of 4 |

## Copies' sheets that differ, left out of A15 (ruling C2b-3/R2: reported, not judged)

| id | view | tie | sheets left out | judged |
|---|---|---|---|---|
| 7QUP | focus | 10B~1B | copies' sheets differ: 10B sheets 3+5 vs 1B sheet 62; 21 strands not judged | yes, 21 matched |
| 7QUP | focus | 10D~1D | copies' sheets differ: 10D sheets 9+11 vs 1D sheet 67; 20 strands not judged | yes, 21 matched |
| 7QUP | focus | 11D~1D | copies' sheets differ: 11D sheets 23+25 vs 1D sheet 67; 20 strands not judged | yes, 22 matched |
| 7QUP | focus | 12D~1D | copies' sheets differ: 12D sheets 38+40 vs 1D sheet 67; 21 strands not judged | yes, 21 matched |
| 7QUP | focus | 13B~1B | copies' sheets differ: 13B sheets 47+49 vs 1B sheet 62; 20 strands not judged | yes, 21 matched |
| 7QUP | focus | 13D~1D | copies' sheets differ: 13D sheets 53+55 vs 1D sheet 67; 20 strands not judged | yes, 21 matched |
| 7QUP | focus | 2B~1B | copies' sheets differ: 2B sheets 75+77 vs 1B sheet 62; 21 strands not judged | yes, 22 matched |
| 7QUP | focus | 2D~1D | copies' sheets differ: 2D sheets 81+83 vs 1D sheet 67; 21 strands not judged | yes, 22 matched |
| 7QUP | focus | 3B~1B | copies' sheets differ: 3B sheets 90+92 vs 1B sheet 62; 21 strands not judged | yes, 22 matched |
| 7QUP | focus | 3D~1D | copies' sheets differ: 3D sheets 96+98 vs 1D sheet 67; 21 strands not judged | yes, 22 matched |
| 7QUP | focus | 4D~1D | copies' sheets differ: 4D sheets 110+112 vs 1D sheet 67; 20 strands not judged | yes, 22 matched |
| 7QUP | focus | 5B~1B | copies' sheets differ: 5B sheets 119+121 vs 1B sheet 62; 21 strands not judged | yes, 22 matched |
| 7QUP | focus | 5D~1D | copies' sheets differ: 5D sheets 125+127 vs 1D sheet 67; 21 strands not judged | yes, 21 matched |
| 7QUP | focus | 7B~1B | copies' sheets differ: 7B sheets 147+149 vs 1B sheet 62; 20 strands not judged | yes, 21 matched |
| 7QUP | focus | 7D~1D | copies' sheets differ: 7D sheets 153+155 vs 1D sheet 67; 21 strands not judged | yes, 21 matched |
| 7QUP | focus | 8B~1B | copies' sheets differ: 8B sheets 162+164 vs 1B sheet 62; 20 strands not judged | yes, 20 matched |
| 7QUP | focus | 8D~1D | copies' sheets differ: 8D sheets 168+170 vs 1D sheet 67; 20 strands not judged | yes, 20 matched |
| 7QUP | focus | 9B~1B | copies' sheets differ: 9B sheets 177+179 vs 1B sheet 62; 21 strands not judged | yes, 22 matched |
| 7QUP | focus | 9D~1D | copies' sheets differ: 9D sheets 183+185 vs 1D sheet 67; 21 strands not judged | yes, 22 matched |
| 9P1X | focus | C~V | copies' sheets differ: C sheet 7 vs V sheet 22; 5 strands not judged | yes, 35 matched |
| 9P1X | focus | G~a | copies' sheets differ: G sheet 33 vs a sheet 39; 11 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | H~b | copies' sheets differ: H sheet 33 vs b sheet 39; 9 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | I~c | copies' sheets differ: I sheet 33 vs c sheet 39; 12 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | L~V | copies' sheets differ: L sheets 41+44 vs V sheets 22+25; 12 strands not judged | yes, 31 matched |
| 9P1X | focus | O~E | copies' sheets differ: O sheet 36 vs E sheets 34+35; 10 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | P~F | copies' sheets differ: P sheet 36 vs F sheet 35; 8 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | Q~a | copies' sheets differ: Q sheet 38 vs a sheet 39; 10 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | R~J | copies' sheets differ: R sheet 36 vs J sheets 34+35; 12 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | S~c | copies' sheets differ: S sheets 37+38 vs c sheet 39; 10 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | T~b | copies' sheets differ: T sheets 37+38 vs b sheet 39; 10 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | Y~E | copies' sheets differ: Y sheet 40 vs E sheets 34+35; 10 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | Z~F | copies' sheets differ: Z sheet 40 vs F sheet 35; 8 strands not judged | n/a: sheets differ (2 matched) |
| 9P1X | focus | d~J | copies' sheets differ: d sheet 40 vs J sheets 34+35; 12 strands not judged | n/a: sheets differ (2 matched) |
| 6JLN | focus | b~B | copies' sheets differ: b sheet 18 vs B sheets 4+5; 15 strands not judged | yes, 25 matched |

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
