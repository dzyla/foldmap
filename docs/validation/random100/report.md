# foldmap random-100 campaign

92 of 100 pass.

An entry passes when no criterion fails. A criterion with nothing to judge in an entry is not judged (or does not apply), never counted as a pass: see the outcomes per criterion below.

Sample: seed 20261013, release-date cutoff 2026-10-13.

| stratum | entries | pass |
|---|---|---|
| S1 | 12 | 12 |
| S2 | 12 | 12 |
| S3 | 10 | 10 |
| S4 | 15 | 14 |
| S5 | 8 | 8 |
| S6 | 8 | 8 |
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
| A8 labels unique per chain, stem names per figure | 100 | 0 | 0 | 0 |
| A9 nothing left out without a note | 100 | 0 | 0 | 0 |
| A10 secondary structure agrees with mkdssp | 96 | 0 | 0 | 4 |
| A11 N->C and sheet order | 100 | 0 | 0 | 0 |
| A12 same figure for any PYTHONHASHSEED | 100 | 0 | 0 | 0 |
| A13 a saved layout redraws the same figure | 100 | 0 | 0 | 0 |
| A14 interactive page and sequence plot written | 100 | 0 | 0 | 0 |
| A15 copies drawn alike | 37 | 0 | 1 | 62 |
| A16 domains: one chain, alike, not tiny, panels apart, names clear | 47 | 1 | 14 | 38 |
| A17 nucleic-acid metrics and stem names by their letters | 12 | 5 | 1 | 82 |
| A18 membrane: TM helices in the band, loops out of the slab | 8 | 0 | 0 | 92 |
| A19 two states: same names in the same order | 8 | 2 | 0 | 90 |
| A20 ligand tethers local | 9 | 0 | 52 | 39 |

### Not judged, per criterion

- **A15** (1): 6YQD (default BBB~AAA n/a: 2 matched)
- **A16** (14): 7IMM (the domains view drew no domain); 7FDN (the domains view drew no domain); 6YF9 (the domains view drew no domain); 1SLV (the domains view drew no domain); 1BXI (the domains view drew no domain); 4UM3 (the domains view drew no domain); 8R5X (the domains view drew no domain); 6AD0 (the domains view drew no domain); 2MJW (the domains view drew no domain); 8W9Y (the domains view drew no domain); 8BW5 (the domains view drew no domain); 6YPX (the domains view drew no domain); 6YQD (the domains view drew no domain); 8P3N (the domains view drew no domain)
- **A17** (1): 7MW8 (no view draws the nucleic acid)
- **A20** (52): 13MV (no tether to another chain); 5Q42 (no tether to another chain); 5HI7 (no tether to another chain); 2YBE (no tether to another chain; reported: ligand-tether:3PG:A1421:1 7.644 (home A); ligand-tether:ALF:A1420:1 9.371 (home A); ligand-tether:LA8:A1419:3 12.74 (home ); 6IDJ (no tether to another chain; reported: ligand-tether:9AU:A406:3 11.81 (home A); ligand-tether:FMN:A404:0 9.631 (home A); ligand-tether:FMN:A404:1 8.066 (home A);); 9C0M (no tether to another chain); 7IMM (no tether to another chain; reported: ligand-tether:A1CNB:A2401:1 9.18 (home A2)); 7LX9 (no tether to another chain); 4XG1 (no tether to another chain; reported: ligand-tether:LLP:A508:1 6.642 (home A); ligand-tether:LLP:A508:3 8.534 (home A); ligand-tether:TME:A507:2 7.343 (home A);); 2C4M (no tether to another chain; reported: ligand-tether:PLP:A1634:1 10.956 (home A); ligand-tether:PLP:B1634:1 8.328 (home B); ligand-tether:PLP:B1634:2 8.366 (home); 6KTL (no tether to another chain; reported: ligand-tether:INS:D402:0 7.349 (home D)); 7OKA (no tether to another chain); 2ILU (no tether to another chain; reported: ligand-tether:NDP:A4482:0 6.603 (home A4)); 7GLY (no tether to another chain; reported: ligand-tether:RDX:A407:3 7.197 (home A); ligand-tether:RDX:B405:3 7.289 (home B)); 7SFH (no tether to another chain); 3L9G (no tether to another chain; reported: ligand-tether:MRD:A1850:2 8.072 (home A1); ligand-tether:URC:A2300:0 14.581 (home A2); ligand-tether:URC:A2300:1 9.621 (ho); 5SB5 (no tether to another chain; reported: ligand-tether:GTP:A501:1 8.179 (home A); ligand-tether:GTP:A501:3 11.837 (home A); ligand-tether:GDP:B501:1 6.905 (home B)); 5DWZ (no tether to another chain; reported: ligand-tether:MPO:C406:0 7.113 (home C); ligand-tether:MPO:C406:1 7.113 (home C)); 1A0Z (no tether to another chain); 1BXI (no tether to another chain); 5X4Z (no tether to another chain); 4UM3 (no tether to another chain; reported: ligand-tether:09R:T301:2 7.762 (home T); ligand-tether:09R:Q301:2 7.762 (home Q); ligand-tether:09R:S301:2 7.762 (home S)); 6DPV (no tether to another chain; reported: ligand-tether:GTP:A501:3 10.062 (home A); ligand-tether:GDP:B501:3 16.068 (home B); ligand-tether:GTP:K501:3 9.767 (home K); 8PQV (no tether to another chain; reported: ligand-tether:ADP:A4701:3 7.463 (home A); ligand-tether:ATP:A4702:2 7.417 (home A); ligand-tether:ATP:A4702:3 6.918 (home ); 6AD0 (no tether to another chain); 9ON5 (no tether to another chain); 7AOD (no tether to another chain); 7QUP (no tether to another chain; reported: ligand-tether:GTP:7A500:3 10.232 (home 7A)); 8RRX (no tether to another chain; reported: ligand-tether:ATP:D5103:1 6.841 (home D); ligand-tether:CFF:D5104:0 7.01 (home D); ligand-tether:CFF:D5104:2 8.007 (home D); 7SMR (no tether to another chain; reported: ligand-tether:CLR:A602:0 11.57 (home A); ligand-tether:CLR:A602:1 8.692 (home A); ligand-tether:CLR:C503:2 9.12 (home C); ); 5JDN (no tether to another chain); 6PS3 (no tether to another chain); 7PIU (no tether to another chain); 5LDW (no tether to another chain; reported: ligand-tether:FMN:F501:1 8.116 (home F); ligand-tether:SF4:F502:3 17.97 (home F)); 9GEP (no tether to another chain; reported: ligand-tether:HEM:N603:3 7.367 (home N); ligand-tether:NAG:N602:1 9.928 (home N)); 3L2V (no tether to another chain); 8BW5 (no tether to another chain); 3MHT (no tether to another chain); 9YLM (no tether to another chain); 7XSP (no tether to another chain); 7MW8 (no tether to another chain); 8VWQ (no tether to another chain); 8SA2 (no tether to another chain; reported: ligand-tether:B1Z:A301:0 6.645 (home A); ligand-tether:B1Z:A301:1 13.548 (home A); ligand-tether:B1Z:A301:2 9.75 (home A);); 3OWW (no tether to another chain); 8XZE (no tether to another chain); 8K8E (no tether to another chain; reported: ligand-tether:CLR:C301:1 8.241 (home C)); 1RM5 (no tether to another chain); 6YPX (no tether to another chain; reported: ligand-tether:CAC:AAA201:3 6.596 (home AAA)); 8P3N (no tether to another chain; reported: ligand-tether:ADZ:B201:3 9.125 (home B)); 8OYQ (no tether to another chain); 6FYF (no tether to another chain); 6FY9 (no tether to another chain)

## Entries

process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 gates; drawing s: the drawing inside it.

| id | stratum | bin | chains | process s | drawing s | A1 | A1c | A2 | A2b | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 | A11 | A12 | A13 | A14 | A15 | A16 | A17 | A18 | A19 | A20 | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 13MV | S1 | M | 1 | 1.3 | 0.46 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5Q42 | S1 | M | 1 | 1.0 | 0.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5HI7 | S1 | M | 1 | 1.2 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 2YBE | S1 | M | 1 | 1.1 | 0.46 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 29RY | S1 | M | 1 | 0.9 | 0.29 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6IDJ | S1 | M | 1 | 1.0 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3RVC | S1 | S | 1 | 1.0 | 0.22 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 1F09 | S1 | S | 1 | 1.05 | 0.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9C0M | S1 | S | 1 | 0.9 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 7IMM | S1 | S | 1 | 1.0 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 9B7C | S1 | S | 1 | 0.95 | 0.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 7LX9 | S1 | S | 1 | 1.0 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 4XG1 | S2 | L | 4 | 2.25 | 1.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2C4M | S2 | L | 4 | 2.81 | 1.67 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6KTL | S2 | L | 4 | 2.3 | 1.3 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7OKA | S2 | L | 30 | 2.46 | 1.57 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2ILU | S2 | M | 1 | 4.16 | 2.78 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7GLY | S2 | M | 2 | 1.25 | 0.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7FDN | S2 | M | 2 | 1.0 | 0.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 7SFH | S2 | M | 1 | 1.75 | 0.96 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9PXS | S2 | S | 1 | 1.35 | 0.54 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 3L9G | S2 | S | 1 | 3.66 | 2.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2B9V | S2 | XL | 16 | 2.0 | 1.09 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - |  |
| 6YF9 | S2 | XL | 180 | 6.26 | 4.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5SB5 | S3 | L | 6 | 3.56 | 2.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 4XEI | S3 | L | 7 | 2.46 | 1.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 5DWZ | S3 | L | 8 | 2.4 | 1.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2I5Y | S3 | L | 8 | 4.76 | 3.91 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 1SLV | S3 | M | 2 | 2.91 | 2.0 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ✓ |  |
| 2QNA | S3 | M | 2 | 1.4 | 0.59 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 1A0Z | S3 | M | 4 | 1.3 | 0.58 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 1BXI | S3 | S | 2 | 1.05 | 0.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 5X4Z | S3 | XL | 24 | 5.61 | 4.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 4UM3 | S3 | XL | 40 | 2.56 | 1.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 6DPV | S4 | L | 12 | 7.17 | 5.07 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 8R5X | S4 | L | 4 | 1.3 | 0.62 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | - |  |
| 6MHS | S4 | L | 4 | 2.91 | 1.91 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 5JB1 | S4 | L | 6 | 4.81 | 3.82 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 8PQV | S4 | L | 1 | 3.81 | 2.48 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6AD0 | S4 | M | 6 | 1.95 | 0.98 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 9ON5 | S4 | M | 3 | 4.36 | 3.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7JU4 | S4 | XL | 61 | 55.06 | 48.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | - | - | - | ✓ | A6: domains 2 (domain-label:D2 · 4 × stub:C:4; domain-label:D3 · A × stub:C:A); A16: domains 2 () |
| 7AOD | S4 | XL | 24 | 12.27 | 9.1 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7QUP | S4 | XL | 65 | 29.2 | 26.8 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9P1X | S4 | XL | 30 | 3.01 | 1.66 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ |  |
| 8RRX | S4 | XL | 12 | 9.37 | 7.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8RVL | S4 | XL | 34 | 8.47 | 6.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 8ESZ | S4 | XL | 43 | 24.69 | 21.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 7ARD | S4 | XL | 51 | 21.59 | 17.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 2MJW | S5 | M | 4 | 1.25 | 0.46 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 2DGS | S5 | S | 1 | 0.8 | 0.13 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8Q4O | S5 | S | 1 | 2.2 | 0.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 1MG8 | S5 | S | 1 | 0.85 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 5KQJ | S5 | S | 1 | 0.9 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8OH7 | S5 | S | 1 | 1.25 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6B5W | S5 | S | 1 | 0.9 | 0.06 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8B4R | S5 | S | 1 | 0.7 | 0.08 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 6JLN | S6 | L | 39 | 31.51 | 27.33 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ |  |
| 7SMR | S6 | L | 5 | 5.51 | 4.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 6JO6 | S6 | L | 19 | 19.49 | 15.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ✓ |  |
| 8W9Y | S6 | L | 6 | 4.06 | 3.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | ✓ | - | - |  |
| 5JDN | S6 | M | 1 | 1.15 | 0.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 6PS3 | S6 | M | 1 | 1.1 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 7PIU | S6 | M | 6 | 2.86 | 1.99 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 5LDW | S6 | XL | 45 | 22.49 | 20.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ○ |  |
| 9K9V | S7 | L | 5 | 4.86 | 1.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | - |  |
| 9GEP | S7 | L | 12 | 8.37 | 4.83 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ |  |
| 3L2V | S7 | M | 4 | 3.71 | 1.51 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ |  |
| 8BW5 | S7 | M | 3 | 3.81 | 1.43 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | ✗ | - | - | ○ | A17: pair_f1 0.8 (limit 0.95); stem_recall 0.5 (limit 0.85) |
| 3MHT | S7 | M | 3 | 2.91 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 9YLM | S7 | S | 3 | 1.9 | 0.23 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | ○ | A17: stems_nearer_frac 0.0 (limit 1.0) |
| 8QEM | S7 | XL | 26 | 2.36 | 1.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | - | - | - | A17: measured on the focus view, as the default draws no nucleic acid: angle_median 3.2 (limit 0.05); angle_max 4.6 (limit 0.05); joint_error_max 3.4 (limit 0.05) |
| 7XSP | S8 | L | 3 | 3.66 | 1.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 7MW8 | S8 | L | 12 | 1.3 | 0.47 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ○ | - | - | ○ |  |
| 9C75 | S8 | M | 3 | 4.38 | 1.74 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | - |  |
| 9CJI | S8 | M | 4 | 4.48 | 1.5 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 6SO9 | S8 | S | 2 | 3.73 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 7Y7C | S8 | XL | 55 | 72.71 | 53.67 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ✓ | A3: domains: strand a:2705-2705, gap after a:1439; A6: domains 1 (na:stem:228:first × na:run:443); A17: na_crossings 4 (limit 0); bead_overlaps 1 (limit 0) |
| 8VWQ | S8 | XL | 56 | 72.02 | 52.41 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | ○ | A3: domains: strand 1:836-836; A6: default 1 (ligand-label:LYS:5101 × na:stem:114); deposited 1 (ligand-label:LYS:5101 × na:stem:114); beta 1 (ligand-label:LYS:5101 × na:stem:114); focus 1 (ligand-label:LYS:5101 × na:stem:114); domains 1 (na:tie:1401-1501:label × na:stem:209); A17: na_crossings 1 (limit 0); bead_overlaps 1 (limit 0) |
| 8SA2 | S9 | M | 2 | 5.43 | 1.68 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 3OWW | S9 | S | 2 | 4.32 | 1.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 8XZE | S9 | S | 1 | 2.1 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 8K8E | S10 | L | 5 | 5.26 | 4.32 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 7BBH | S10 | L | 3 | 4.36 | 3.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 1WWL | S10 | M | 2 | 2.46 | 1.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 7BRG | S10 | M | 3 | 2.25 | 1.5 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 8OZJ | S10 | XL | 11 | 18.24 | 17.0 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 1RM5 | S11 | M | 3 | 2.76 | 1.73 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ○ |  |
| 2HKI | S11 | M | 1 | 3.31 | 2.16 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - |  |
| 9YKJ | S11 | L | 2 | 2.71 | 1.62 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | ✗ | - | A19: names compared on all 43 pairs, 27 differ; αN~α2 (A:37-47~A:37-50); αN′~α3 (A:72-85~A:72-85); α1~α4 (A:100-108~A:98-111); α2~α7 (A:164-175~A:164-174); α3~α8 (A:184-192~A:184-198); α3″~α9 (A:204-210~A:203-210); α4~α10 (A:221-226~A:221-226); α4′~α11 (A:228-236~A:228-237) (+49 more) |
| 10UJ | S11 | M | 1 | 1.9 | 0.95 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✗ | - | A19: names compared on all 43 pairs, 27 differ; α2~αN (A:37-50~A:37-47); α3~αN′ (A:72-85~A:72-85); α4~α1 (A:98-111~A:100-108); α7~α2 (A:164-174~A:164-175); α8~α3 (A:184-198~A:184-192); α9~α3″ (A:203-210~A:204-210); α10~α4 (A:221-226~A:221-226); α11~α4′ (A:228-237~A:228-236) (+49 more) |
| 6YPX | S11 | M | 4 | 1.2 | 0.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | ✓ | ○ |  |
| 6YQD | S11 | M | 2 | 0.75 | 0.2 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | ○ | - | - | ✓ | - |  |
| 8P3N | S11 | M | 2 | 1.25 | 0.45 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | ✓ | ○ |  |
| 8OYQ | S11 | M | 2 | 0.85 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 6FYF | S11 | L | 4 | 1.35 | 0.55 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 6FY9 | S11 | M | 1 | 1.3 | 0.46 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| AF_AFA0A4E9ES35F1 | S12 | M | 1 | 1.05 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFO42816F1 | S12 | M | 1 | 1.35 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFB1Q234F1 | S12 | S | 1 | 0.95 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |

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
| 13MV | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 435.02 | 164.0 |
| 5Q42 | 37 | 37 | None | 0 | 0 | 0 | 0 | 0 | 0 | 424.94 | 163.3 |
| 5HI7 | 15 | 15 | None | 0 | 0 | 0 | 0 | 0 | 0 | 592.21 | 168.1 |
| 2YBE | 35 | 35 | None | 0 | 0 | 0 | 1 | 0 | 0 | 521.66 | 165.8 |
| 29RY | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 464.93 | 164.3 |
| 6IDJ | 17 | 17 | None | 0 | 0 | 0 | 0 | 0 | 0 | 511.73 | 163.1 |
| 3RVC | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 167.63 | 160.0 |
| 1F09 | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 252.21 | 161.5 |
| 9C0M | 22 | 22 | None | 0 | 0 | 0 | 0 | 0 | 0 | 331.75 | 162.8 |
| 7IMM | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 343.58 | 162.5 |
| 9B7C | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 106.97 | 161.7 |
| 7LX9 | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 159.92 | 160.3 |
| 4XG1 | 131 | 131 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1294.48 | 197.9 |
| 2C4M | 147 | 147 | None | 0 | 0 | 0 | 3 | 0 | 0 | 2326.98 | 306.5 |
| 6KTL | 96 | 96 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1785.27 | 218.0 |
| 7OKA | 223 | 223 | None | 0 | 0 | 0 | 0 | 0 | 0 | 799.41 | 299.6 |
| 2ILU | 311 | 311 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3320.71 | 333.3 |
| 7GLY | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 565.53 | 171.7 |
| 7FDN | 15 | 15 | None | 0 | 0 | 0 | 0 | 0 | 0 | 418.44 | 163.0 |
| 7SFH | 20 | 20 | None | 0 | 0 | 0 | 1 | 0 | 0 | 705.48 | 177.2 |
| 9PXS | 19 | 19 | None | 0 | 0 | 0 | 0 | 0 | 0 | 785.08 | 192.3 |
| 3L9G | 213 | 213 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1629.75 | 211.5 |
| 2B9V | 61 | 61 | None | 0 | 0 | 0 | 0 | 0 | 0 | 770.37 | 270.9 |
| 6YF9 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 138.13 | 582.4 |
| 5SB5 | 129 | 129 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3663.2 | 356.8 |
| 4XEI | 44 | 44 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2585.62 | 264.8 |
| 5DWZ | 137 | 137 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1541.88 | 207.6 |
| 2I5Y | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1336.58 | 199.2 |
| 1SLV | 69 | 69 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1200.7 | 169.7 |
| 2QNA | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1580.29 | 195.7 |
| 1A0Z | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 735.75 | 171.4 |
| 1BXI | 8 | 8 | None | 0 | 0 | 0 | 1 | 0 | 0 | 279.91 | 161.8 |
| 5X4Z | 310 | 310 | None | 0 | 0 | 0 | 1 | 0 | 0 | 4803.34 | 455.8 |
| 4UM3 | 35 | 35 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1190.0 | 304.4 |
| 6DPV | 427 | 427 | None | 0 | 0 | 0 | 5 | 0 | 0 | 7264.62 | 514.5 |
| 8R5X | 53 | 53 | None | 0 | 0 | 0 | 0 | 0 | 0 | 946.98 | 167.5 |
| 6MHS | 37 | 37 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3381.91 | 345.1 |
| 5JB1 | 486 | 486 | None | 0 | 0 | 0 | 2 | 0 | 0 | 3335.07 | 431.3 |
| 8PQV | 201 | 201 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4188.14 | 250.8 |
| 6AD0 | 66 | 66 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1133.97 | 171.4 |
| 9ON5 | 127 | 127 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1605.1 | 215.1 |
| 7JU4 | 2571 | 2571 | None | 0 | 0 | 0 | 6 | 0 | 0 | 29826.53 | 594.8 |
| 7AOD | 499 | 499 | None | 0 | 0 | 0 | 0 | 0 | 0 | 8034.48 | 515.6 |
| 7QUP | 25 | 25 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1077.17 | 710.3 |
| 9P1X | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1693.39 | 421.0 |
| 8RRX | 300 | 300 | None | 0 | 0 | 0 | 0 | 0 | 0 | 9410.73 | 498.7 |
| 8RVL | 483 | 483 | None | 0 | 0 | 0 | 2 | 0 | 0 | 11585.47 | 481.3 |
| 8ESZ | 279 | 279 | None | 0 | 0 | 0 | 4 | 0 | 0 | 14940.61 | 464.6 |
| 7ARD | 598 | 598 | None | 0 | 0 | 0 | 1 | 0 | 0 | 19197.26 | 499.4 |
| 2MJW | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 402.21 | 247.9 |
| 2DGS | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 49.82 | 189.9 |
| 8Q4O | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 57.85 | 174.9 |
| 1MG8 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 36.96 | 161.5 |
| 5KQJ | 5 | 5 | None | 0 | 0 | 0 | 0 | 0 | 0 | 384.61 | 208.9 |
| 8OH7 | 3 | 3 | None | 0 | 0 | 0 | 0 | 0 | 0 | 429.12 | 164.5 |
| 6B5W | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 157.5 |
| 8B4R | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 12.71 | 159.6 |
| 6JLN | 137 | 137 | None | 0 | 0 | 0 | 4 | 0 | 0 | 12050.57 | 502.8 |
| 7SMR | 85 | 85 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3148.74 | 333.9 |
| 6JO6 | 180 | 180 | None | 0 | 0 | 0 | 2 | 0 | 0 | 6463.28 | 415.2 |
| 8W9Y | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1327.36 | 196.5 |
| 5JDN | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 317.04 | 162.5 |
| 6PS3 | 9 | 9 | None | 0 | 0 | 0 | 1 | 0 | 0 | 612.53 | 165.3 |
| 7PIU | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1324.58 | 196.2 |
| 5LDW | 145 | 145 | None | 0 | 0 | 0 | 0 | 0 | 0 | 10373.41 | 418.6 |
| 9K9V | 127 | 127 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2388.63 | 275.6 |
| 9GEP | 35 | 35 | 3 | 0 | 0 | 0 | 2 | 0 | 0 | 4094.72 | 272.3 |
| 3L2V | 51 | 51 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2143.77 | 217.7 |
| 8BW5 | 8 | 8 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 720.85 | 174.5 |
| 3MHT | 18 | 18 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 396.7 | 173.0 |
| 9YLM | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 155.92 | 166.2 |
| 8QEM | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1473.25 | 382.6 |
| 7XSP | 101 | 101 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1757.19 | 194.6 |
| 7MW8 | 34 | 34 | None | 0 | 0 | 0 | 2 | 0 | 0 | 490.56 | 168.8 |
| 9C75 | 25 | 25 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 1141.8 | 179.7 |
| 9CJI | 53 | 53 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2034.39 | 199.0 |
| 6SO9 | 4 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 86.66 | 240.7 |
| 7Y7C | 335 | 335 | 724 | 0 | 0 | 0 | 1 | 0 | 0 | 31653.95 | 629.2 |
| 8VWQ | 443 | 443 | 637 | 0 | 0 | 0 | 1 | 0 | 0 | 32645.44 | 665.3 |
| 8SA2 | 0 | 0 | 13 | 0 | 0 | 0 | 0 | 0 | 0 | 1973.64 | 193.5 |
| 3OWW | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 695.11 | 179.5 |
| 8XZE | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 122.15 | 164.9 |
| 8K8E | 99 | 99 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2311.75 | 240.4 |
| 7BBH | 52 | 52 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1415.51 | 255.6 |
| 1WWL | 184 | 184 | None | 0 | 0 | 0 | 0 | 0 | 0 | 539.7 | 175.3 |
| 7BRG | 71 | 71 | None | 0 | 0 | 0 | 0 | 0 | 0 | 869.36 | 176.7 |
| 8OZJ | 242 | 242 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4678.64 | 440.1 |
| 1RM5 | 97 | 97 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1541.86 | 220.7 |
| 2HKI | 87 | 87 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1526.49 | 209.4 |
| 9YKJ | 85 | 85 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1707.13 | 230.8 |
| 10UJ | 83 | 83 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1232.29 | 168.6 |
| 6YPX | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 186.98 | 161.6 |
| 6YQD | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 205.22 | 162.7 |
| 8P3N | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 375.61 | 167.2 |
| 8OYQ | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 154.1 | 162.7 |
| 6FYF | 36 | 36 | None | 0 | 0 | 0 | 0 | 0 | 0 | 466.19 | 166.8 |
| 6FY9 | 33 | 33 | None | 0 | 0 | 0 | 0 | 0 | 0 | 593.76 | 164.5 |
| AF_AFA0A4E9ES35F1 | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 555.38 | 166.8 |
| AF_AFO42816F1 | 27 | 27 | None | 0 | 0 | 0 | 0 | 0 | 0 | 820.0 | 164.9 |
| AF_AFB1Q234F1 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 176.37 | 159.5 |

## Science checks (S1-S6), apart from A1-A20

Judged in every drawn view (default, deposited, beta, focus, domains, membrane, app); S5 once per entry.

A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet drawn, no barrel, no residue mark, no accessory strand, no CATH reference) is not judged, one with no protein does not apply.

| check | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| S1 sheet neighbours pair | 90 | 0 | 6 | 4 |
| S2 barrel closure pairs | 13 | 0 | 83 | 4 |
| S3 barrel strand count is its pairing cycle | 90 | 0 | 6 | 4 |
| S4 residue marks at their residue | 2 | 0 | 94 | 4 |
| S5 automatic domains agree with a reference | 5 | 13 | 78 | 4 |
| S6 accessory drawn at its partner | 37 | 0 | 59 | 4 |

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
| 3RVC | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 1F09 | ✓ | ✓ | ✓ | ✓ | ○ | ○ |  |
| 9C0M | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7IMM | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9B7C | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7LX9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4XG1 | ✓ | ✓ | ✓ | ○ | ✓ | ✓ |  |
| 2C4M | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.62, boundary F1 0.57, recall 1.00 (2 reference domains, 5 automatic) |
| 6KTL | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7OKA | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2ILU | ✓ | ○ | ✓ | ○ | ✓ | ○ |  |
| 7GLY | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7FDN | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7SFH | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 9PXS | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3L9G | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 2B9V | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.76, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic) |
| 6YF9 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5SB5 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 4XEI | ✓ | ○ | ✓ | ○ | ✗ | ✓ | S5 A: overlap 0.92, boundary F1 0.83, recall 0.83 (4 reference domains, 4 automatic) |
| 5DWZ | ✓ | ○ | ✓ | ○ | ✓ | ✓ |  |
| 2I5Y | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 1SLV | ✓ | ✓ | ✓ | ○ | ✗ | ✓ | S5 B: overlap 0.54, boundary F1 0.00, recall 0.00 (2 reference domains, 1 automatic) |
| 2QNA | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1A0Z | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1BXI | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5X4Z | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 4UM3 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6DPV | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.72, boundary F1 0.50, recall 0.50 (3 reference domains, 2 automatic) |
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
| 2MJW | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2DGS | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8Q4O | - | - | - | - | - | - |  |
| 1MG8 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5KQJ | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
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
| 1WWL | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
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
| 6FYF | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.75, boundary F1 0.86, recall 1.00 (3 reference domains, 4 automatic) |
| 6FY9 | ✓ | ○ | ✓ | ○ | ✗ | ○ | S5 A: overlap 0.60, boundary F1 0.83, recall 1.00 (3 reference domains, 5 automatic) |
| AF_AFA0A4E9ES35F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFO42816F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFB1Q234F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |

## S5: automatic domains against CATH (reported, not gated)

In-sample: these chains are among the 33 the automatic domains were calibrated on (tests/data/domains/cath_reference.json); the never-seen check is the fresh set (cath_fresh.json).

| id | chain | CATH domains | automatic | overlap | boundary F1 | recall | precision |
|---|---|---|---|---|---|---|---|
| 5Q42 | A | 2 | 3 | 0.77 | 0.80 | 1.00 | 0.67 |
| 5HI7 | A | 5 | 4 | 0.78 | 0.89 | 0.80 | 1.00 |
| 2YBE | A | 2 | 2 | 0.95 | 0.67 | 1.00 | 0.50 |
| 4XG1 | A | 2 | 2 | 0.98 | 1.00 | 1.00 | 1.00 |
| 2C4M | A | 2 | 5 | 0.62 | 0.57 | 1.00 | 0.40 |
| 2ILU | A | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2B9V | A | 3 | 4 | 0.76 | 0.86 | 1.00 | 0.75 |
| 4XEI | A | 4 | 4 | 0.92 | 0.83 | 0.83 | 0.83 |
| 5DWZ | C | 2 | 2 | 0.99 | 1.00 | 1.00 | 1.00 |
| 1SLV | B | 2 | 1 | 0.54 | 0.00 | 0.00 | 0.00 |
| 6DPV | A | 3 | 2 | 0.72 | 0.50 | 0.50 | 0.50 |
| 5JDN | A | 2 | 2 | 0.90 | 0.67 | 0.67 | 0.67 |
| 3L2V | A | 4 | 3 | 0.84 | 0.86 | 0.75 | 1.00 |
| 3MHT | A | 2 | 2 | 0.92 | 0.50 | 0.50 | 0.50 |
| 1RM5 | A | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2HKI | A | 2 | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| 6FYF | A | 3 | 4 | 0.75 | 0.86 | 1.00 | 0.75 |
| 6FY9 | A | 3 | 5 | 0.60 | 0.83 | 1.00 | 0.71 |
| mean of 18 | | | | 0.84 | 0.77 | 0.84 | 0.74 |

## Accessories named, not drawn level beside their partner (reported, not gated)

| id | named of accessory strands, per view |
|---|---|
| 13MV | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 3RVC | default 1 of 1; deposited 1 of 1; beta 1 of 1; seed1 1 of 1; roundtrip 1 of 1 |
| 4XG1 | default 2 of 6; deposited 2 of 6; beta 2 of 6; focus 2 of 6; domains 2 of 6; seed1 2 of 6; roundtrip 2 of 6 |
| 2C4M | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2 |
| 7GLY | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 1; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 7SFH | default 0 of 2; deposited 0 of 2; beta 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 9PXS | default 2 of 2; deposited 2 of 2; beta 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2 |
| 3L9G | default 8 of 16; deposited 8 of 16; beta 8 of 16; domains 8 of 16; seed1 8 of 16; roundtrip 8 of 16 |
| 2B9V | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 4; seed1 0 of 1; roundtrip 0 of 1 |
| 5SB5 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 4XEI | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 1; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 5DWZ | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2 |
| 2I5Y | default 2 of 2; deposited 2 of 2; beta 2 of 2; focus 2 of 2; domains 2 of 2; seed1 2 of 2; roundtrip 2 of 2 |
| 1SLV | default 0 of 6; deposited 0 of 6; beta 0 of 6; focus 0 of 1; domains 0 of 6; seed1 0 of 6; roundtrip 0 of 6 |
| 5X4Z | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 0 of 2; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 4UM3 | default 0 of 5; deposited 0 of 5; beta 0 of 5; focus 0 of 5; domains 0 of 5; seed1 0 of 5; roundtrip 0 of 5 |
| 8R5X | default 3 of 4; deposited 3 of 4; beta 3 of 4; focus 0 of 0; domains 3 of 4; seed1 3 of 4; roundtrip 3 of 4 |
| 5JB1 | default 2 of 12; deposited 2 of 12; beta 2 of 12; focus 0 of 2; domains 2 of 12; seed1 2 of 12; roundtrip 2 of 12 |
| 6AD0 | default 0 of 1; deposited 0 of 1; beta 0 of 1; focus 0 of 0; domains 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 9ON5 | default 6 of 15; deposited 6 of 15; beta 6 of 15; focus 2 of 5; domains 6 of 15; seed1 6 of 15; roundtrip 6 of 15 |
| 7JU4 | default 3 of 7; deposited 3 of 7; beta 3 of 7; focus 0 of 1; domains 3 of 7; seed1 3 of 7; roundtrip 3 of 7 |
| 7AOD | default 4 of 6; deposited 4 of 6; beta 4 of 6; focus 0 of 0; domains 4 of 6; seed1 4 of 6; roundtrip 4 of 6 |
| 9P1X | default 0 of 2; deposited 0 of 2; beta 0 of 2; focus 8 of 26; domains 0 of 2; seed1 0 of 2; roundtrip 0 of 2 |
| 8RRX | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 8 of 12; domains 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
| 8RVL | default 0 of 6; deposited 0 of 6; beta 0 of 6; focus 0 of 0; domains 0 of 6; seed1 0 of 6; roundtrip 0 of 6 |
| 7ARD | default 2 of 8; deposited 2 of 8; beta 2 of 8; focus 0 of 0; domains 2 of 8; seed1 2 of 8; roundtrip 2 of 8 |
| 5KQJ | default 0 of 1; deposited 0 of 1; beta 0 of 1; seed1 0 of 1; roundtrip 0 of 1 |
| 6JLN | default 2 of 3; deposited 2 of 3; beta 2 of 3; focus 0 of 0; domains 2 of 3; membrane 2 of 3; seed1 2 of 3; roundtrip 2 of 3 |
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
| 6JLN | default | b~B | copies' sheets differ: b sheet 18 vs B sheets 4+5; 15 strands not judged | yes, 25 matched |

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
