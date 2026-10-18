# foldmap random-100 campaign

77 of 100 pass.

An entry passes when no criterion fails. A criterion with nothing to judge in an entry is not judged (or does not apply), never counted as a pass: see the outcomes per criterion below.

Sample: seed 20261017, release-date cutoff 2026-10-13.

| stratum | entries | pass |
|---|---|---|
| S1 | 12 | 11 |
| S2 | 12 | 11 |
| S3 | 10 | 10 |
| S4 | 15 | 10 |
| S5 | 8 | 8 |
| S6 | 8 | 4 |
| S7 | 7 | 3 |
| S8 | 7 | 3 |
| S9 | 3 | 3 |
| S10 | 5 | 5 |
| S11 | 10 | 6 |
| S12 | 3 | 3 |

## Outcomes per criterion

✓ pass, ✗ fail, ○ not judged (the criterion applies but the entry gave it nothing to judge: no tie, no tether to another chain, no domain drawn, mkdssp unmeasured ...), - does not apply.

| criterion | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| A1 no crash: exit 0, files written, no traceback, no Python warning | 99 | 1 | 0 | 0 |
| A1c a CA-only or P-only entry exits non-zero with foldmap's one-line message | 0 | 0 | 0 | 100 |
| A2 default view time | 100 | 0 | 0 | 0 |
| A2b the app draws it within the hosted limits | 100 | 0 | 0 | 0 |
| A3 no route fallbacks | 90 | 10 | 0 | 0 |
| A4 no loop over a label | 97 | 3 | 0 | 0 |
| A5 no element or stem boxes overlapping | 98 | 2 | 0 | 0 |
| A6 no text collisions | 89 | 11 | 0 | 0 |
| A7 text and legend inside the page | 99 | 1 | 0 | 0 |
| A8 labels unique per chain, stem names per figure | 99 | 1 | 0 | 0 |
| A9 nothing left out without a note | 100 | 0 | 0 | 0 |
| A10 secondary structure agrees with mkdssp | 97 | 0 | 0 | 3 |
| A11 N->C and sheet order | 100 | 0 | 0 | 0 |
| A12 same figure for any PYTHONHASHSEED | 100 | 0 | 0 | 0 |
| A13 a saved layout redraws the same figure | 100 | 0 | 0 | 0 |
| A14 interactive page and sequence plot written | 100 | 0 | 0 | 0 |
| A15 copies drawn alike | 30 | 5 | 0 | 65 |
| A16 domains: one chain, alike, not tiny, panels apart, names clear | 47 | 4 | 12 | 37 |
| A17 nucleic-acid metrics and stem names by their letters | 11 | 6 | 0 | 83 |
| A18 membrane: TM helices in the band, loops out of the slab | 5 | 3 | 0 | 92 |
| A19 two states: same names in the same order | 6 | 4 | 0 | 90 |
| A20 ligand tethers local | 12 | 0 | 53 | 35 |

### Not judged, per criterion

- **A16** (12): 3KSI (the domains view drew no domain); 5YUH (the domains view drew no domain); 6B7J (the domains view drew no domain); 2O4P (the domains view drew no domain); 3B6N (the domains view drew no domain); 3KTK (the domains view drew no domain); 2YJA (the domains view drew no domain); 3JBE (the domains view drew no domain); 1TKW (the domains view drew no domain); 2VBO (the domains view drew no domain); 5K7Z (the domains view drew no domain); 9WB3 (the domains view drew no domain)
- **A20** (53): 1PTY (no tether to another chain; reported: ligand-tether:PTR:A323:1 7.246 (home A)); 1HMY (no tether to another chain); 6EJ3 (no tether to another chain; reported: ligand-tether:B7T:A501:0 7.443 (home A); ligand-tether:B7T:A501:1 6.648 (home A)); 5RFU (no tether to another chain); 8F03 (no tether to another chain); 5YUH (no tether to another chain); 2W6W (no tether to another chain); 3LIL (no tether to another chain); 7P46 (no tether to another chain; reported: ligand-tether:HEM+CYN:B301:1 7.452 (home B); ligand-tether:HEM+CYN:A301:2 6.984 (home A); ligand-tether:KYN:D303:0 7.308 (); 8HRP (no tether to another chain; reported: ligand-tether:NAD:A401:2 7.085 (home A); ligand-tether:G3H:C402:1 8.142 (home C); ligand-tether:NAD:D401:2 8.264 (home D)); 6SNU (no tether to another chain; reported: ligand-tether:PLP:A502:1 7.162 (home A); ligand-tether:PLP:A502:3 6.731 (home A); ligand-tether:PLP:B504:1 7.148 (home B);); 3AMZ (no tether to another chain; reported: ligand-tether:FAD:A3005:2 10.299 (home A); ligand-tether:FES:A3001:2 19.25 (home A); ligand-tether:ZXF:A3009:2 7.767 (home); 8ZT9 (no tether to another chain); 7X8H (no tether to another chain); 2NAD (no tether to another chain); 2O4P (no tether to another chain; reported: ligand-tether:TPV:A300:0 7.055 (home A)); 3B6N (no tether to another chain); 5EUJ (no tether to another chain; reported: ligand-tether:TPP:A602:2 6.904 (home A); ligand-tether:TPP:A602:3 6.548 (home A)); 9SLB (no tether to another chain); 6HW8 (no tether to another chain; reported: ligand-tether:MG:G301:2 9.455 (home G); ligand-tether:GT8:H301:0 7.81 (home H); ligand-tether:GT8:V301:0 6.703 (home V); l); 6XER (no tether to another chain; reported: ligand-tether:GTP:A501:1 9.82 (home A); ligand-tether:GTP:A501:2 6.856 (home A); ligand-tether:GDP:B501:3 7.63 (home B); l); 3ULQ (no tether to another chain); 2UZU (no tether to another chain); 6N9A (no tether to another chain; reported: ligand-tether:ZN:D1402:2 8.962 (home D1); ligand-tether:KG4:D1401:3 8.415 (home D1); ligand-tether:ZN:D2402:2 8.983 (home ); 2YJA (no tether to another chain; reported: ligand-tether:EST:B11550:0 7.343 (home B1); ligand-tether:EST:B11550:3 9.438 (home B1); ligand-tether:EST:B21550:0 7.343 (); 4N9F (no tether to another chain); 3L4G (no tether to another chain); 8D3D (no tether to another chain; reported: ligand-tether:CA:A2001:4 8.03 (home A)); 6V6S (no tether to another chain; reported: ligand-tether:ADP:U800:0 15.067 (home U); ligand-tether:ADP:U800:1 25.094 (home U); ligand-tether:ADP:U800:2 12.144 (home ); 12BM (no tether to another chain; reported: ligand-tether:ADP:A501:2 10.784 (home A); ligand-tether:ADP:A501:3 7.069 (home A); ligand-tether:AGS:C501:1 17.008 (home C); 7BKC (no tether to another chain; reported: ligand-tether:FAD:A702:3 13.647 (home A); ligand-tether:FAD:a702:0 7.359 (home a); ligand-tether:FAD:a702:3 11.157 (home a); 9Y05 (no tether to another chain; reported: ligand-tether:JDP:A901:2 6.657 (home A)); 1TKW (no tether to another chain); 6MJ2 (no tether to another chain; reported: ligand-tether:CA:A1601:2 6.912 (home A)); 6QFA (no tether to another chain); 8JWF (no tether to another chain; reported: ligand-tether:YMZ:A1500:1 6.779 (home A); ligand-tether:YMZ:A1500:2 8.191 (home A)); 8WH0 (no tether to another chain; reported: ligand-tether:ANP:A802:3 8.21 (home A); ligand-tether:ANP:C802:1 7.546 (home C); ligand-tether:ANP:F802:1 6.752 (home F); ); 4NLG (no tether to another chain); 2FD8 (no tether to another chain); 13FZ (no tether to another chain; reported: ligand-tether:AGS:A2000:2 11.751 (home A); ligand-tether:AGS:C2000:2 15.164 (home C); ligand-tether:AGS:C2000:3 10.606 (ho); 4E8Q (no tether to another chain); 9IWF (no tether to another chain); 9YER (no tether to another chain); 9QK1 (no tether to another chain); 3IB0 (no tether to another chain; reported: ligand-tether:CO3:A687:3 6.71 (home A)); 6A4Y (no tether to another chain; reported: ligand-tether:HEM:A601:0 8.952 (home A)); 8JXD (no tether to another chain); 1C6A (no tether to another chain); 244L (no tether to another chain); 8U9T (no tether to another chain; reported: ligand-tether:W1L:A1401:3 7.845 (home A1); ligand-tether:W1L:A2401:3 7.884 (home A2)); 3U3H (no tether to another chain; reported: ligand-tether:03W+FMT:A1400:1 10.047 (home A1); ligand-tether:MRD:A12000:2 7.095 (home A1); ligand-tether:MRD:A32000:2 7.0); 6QND (no tether to another chain; reported: ligand-tether:MG:A402:0 8.701 (home A); ligand-tether:MG:A402:2 10.774 (home A); ligand-tether:MG:A402:3 7.584 (home A); l); 2BBH (no tether to another chain; reported: ligand-tether:DMU:A402:2 6.502 (home A); ligand-tether:DMU:A402:3 9.734 (home A))

## Entries

process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 gates; drawing s: the drawing inside it.

| id | stratum | bin | chains | process s | drawing s | A1 | A1c | A2 | A2b | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 | A11 | A12 | A13 | A14 | A15 | A16 | A17 | A18 | A19 | A20 | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1PTY | S1 | M | 1 | 0.85 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6TX2 | S1 | M | 1 | 1.25 | 0.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 7N71 | S1 | M | 1 | 1.4 | 0.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 1HMY | S1 | M | 1 | 1.15 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6EJ3 | S1 | M | 1 | 1.55 | 0.74 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 5RFU | S1 | M | 1 | 1.45 | 0.62 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2FR2 | S1 | S | 1 | 0.8 | 0.22 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - | A3: default: A:89-100>A:103-109, A:115-118>A:125-135; deposited: A:89-100>A:103-109, A:115-118>A:125-135; beta: A:89-100>A:103-109, A:115-118>A:125-135 |
| 8F03 | S1 | S | 1 | 1.45 | 0.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3KSI | S1 | S | 1 | 1.35 | 0.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5YUH | S1 | S | 1 | 0.9 | 0.24 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 2W6W | S1 | S | 1 | 0.9 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 3LIL | S1 | S | 1 | 0.85 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 7P46 | S2 | L | 8 | 2.35 | 1.51 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8HRP | S2 | L | 8 | 2.25 | 1.34 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6SNU | S2 | L | 4 | 1.95 | 1.12 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 3AMZ | S2 | L | 2 | 4.41 | 2.98 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 8ZT9 | S2 | M | 2 | 1.7 | 0.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 7X8H | S2 | M | 1 | 1.65 | 0.88 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2NAD | S2 | M | 2 | 1.65 | 0.87 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ | A3: default: B:250-253>B:267-270, B:267-270>B:277-281; deposited: B:250-253>B:267-270, B:267-270>B:277-281; beta: B:250-253>B:267-270, B:267-270>B:277-281; domains: B:250-253>B:267-270, B:267-270>B:277-281; A5: default 2 (['B:250-253', 'B:267-270']; ['B:267-270', 'B:277-281']); deposited 2 (['B:250-253', 'B:267-270']; ['B:267-270', 'B:277-281']); beta 2 (['B:250-253', 'B:267-270']; ['B:267-270', 'B:277-281']); domains 2 (['B:250-253', 'B:267-270']; ['B:267-270', 'B:277-281']); A6: default 2 (label:B:267-270 × label:B:277-281; label:B:267-270 × strand:B:277-281); deposited 2 (label:B:267-270 ×  |
| 6B7J | S2 | M | 2 | 0.95 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 2O4P | S2 | S | 2 | 0.85 | 0.27 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 3B6N | S2 | S | 1 | 1.2 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 5EUJ | S2 | XL | 24 | 2.35 | 1.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 9SLB | S2 | XL | 10 | 2.36 | 1.23 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ○ |  |
| 3KTK | S3 | L | 28 | 1.25 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | - |  |
| 5S5Z | S3 | L | 6 | 4.11 | 2.75 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 6HW8 | S3 | L | 28 | 8.62 | 5.91 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 6XER | S3 | L | 5 | 3.66 | 2.42 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 3ULQ | S3 | M | 2 | 1.15 | 0.43 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 2UZU | S3 | M | 2 | 1.05 | 0.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 6N9A | S3 | M | 3 | 2.51 | 1.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 2YJA | S3 | S | 2 | 1.1 | 0.41 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ○ | - | - | - | ○ |  |
| 4N9F | S3 | XL | 60 | 1.85 | 0.82 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 3L4G | S3 | XL | 16 | 2.96 | 1.69 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 8FLU | S4 | L | 6 | 2.71 | 1.81 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 6ET5 | S4 | L | 54 | 6.91 | 4.49 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 7KH0 | S4 | L | 7 | 3.56 | 2.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 8BCU | S4 | L | 9 | 7.37 | 6.17 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - | A3: default: G:92-99>G:111-114, G:119-122>G:135-140, G:150-158>G:166-173, G:193-200>G:209-215; deposited: G:92-99>G:111-114, G:119-122>G:135-140, G:150-158>G:166-173, G:193-200>G:209-215; beta: G:92-99>G:111-114, G:119-122>G:135-140, G:150-158>G:166-173, G:193-200>G:209-215; domains: G:92-99>G:111-114, G:119-122>G:135-140, G:150-158>G:166-173, G:193-200>G:209-215 |
| 7ZZ1 | S4 | L | 2 | 1.15 | 0.44 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 3JBE | S4 | M | 5 | 2.0 | 1.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | - |  |
| 6WUT | S4 | M | 3 | 1.4 | 0.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | - |  |
| 8D3D | S4 | XL | 16 | 22.09 | 19.43 | ✗ | - | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | - | ○ | A1: focus: killed after 900 s; A5: focus not measured (not measured); A6: focus not measured (not measured); A7: focus not measured (not measured); A8: focus not measured (not measured); A15: focus not measured (not measured) |
| 6V6S | S4 | XL | 38 | 31.06 | 26.34 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 12BM | S4 | XL | 19 | 12.32 | 9.79 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | - | - | - | ○ | A3: default: A:93-101>A:111-115, B:98-101>B:107-116, C:21-66>C:71-79, C:84-89>C:94-98; deposited: A:93-101>A:111-115, B:98-101>B:107-116, C:21-66>C:71-79, C:84-89>C:94-98; beta: A:93-101>A:111-115, B:98-101>B:107-116, C:21-66>C:71-79, C:84-89>C:94-98; domains: A:93-101>A:111-115, B:98-101>B:107-116, C:21-66>C:71-79, C:84-89>C:94-98; A6: domains 1 (domain-label:D1 · B × loop:B:98-101>B:107-116); A16: domains 1 () |
| 9GTP | S4 | XL | 60 | 5.91 | 3.96 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | - | - | A15: focus 52.108405597395176 |
| 7BKC | S4 | XL | 26 | 18.39 | 13.76 | ✓ | - | ✓ | ✓ | ✗ | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | - | - | - | ○ | A3: domains: E:15-20>E:26-37, E:188-192>E:203-206, C:27-37>C:41-43, C:49-53>C:57-59; A4: domains: D:183-195>D:198-201, d:183-195>d:198-201; A6: domains 70 (domain-label:D2 · A × loop:A:351-353>A:365-367; domain-label:D2 · A × loop:A:368-375>A:396-403; domain-label:D2 · D × loop:D:170-174>D:183-195); A16: domains 6 () |
| 7SPB | S4 | XL | 78 | 2.15 | 1.02 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - |  |
| 7ZDP | S4 | XL | 45 | 14.02 | 11.25 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ✓ |  |
| 9Y05 | S4 | XL | 12 | 2.31 | 1.13 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 1TKW | S5 | M | 2 | 1.1 | 0.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | - | - | - | ○ |  |
| 2MRJ | S5 | S | 1 | 1.05 | 0.21 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 5GIW | S5 | S | 1 | 0.7 | 0.07 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2LZJ | S5 | S | 1 | 0.75 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 2N5J | S5 | S | 1 | 0.85 | 0.14 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 8SKX | S5 | S | 1 | 0.85 | 0.11 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 1RCH | S5 | S | 1 | 1.05 | 0.24 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9BVV | S5 | S | 1 | 0.65 | 0.03 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| 9CZH | S6 | L | 8 | 3.01 | 2.09 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ |  |
| 6MJ2 | S6 | L | 4 | 2.15 | 1.04 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | - | ○ | A18: no band drawn |
| 6QFA | S6 | L | 10 | 1.85 | 0.95 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ○ |  |
| 8WMV | S6 | L | 28 | 32.16 | 25.02 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ✓ |  |
| 8JWF | S6 | M | 1 | 2.0 | 1.09 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ✓ | - | ○ |  |
| 7JI3 | S6 | M | 3 | 1.3 | 0.61 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | - | - | A18: loop in the slab C:101-108>C:116-119 |
| 7VVK | S6 | M | 6 | 2.51 | 1.59 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✗ | - | - | A18: no band drawn |
| 6ZKT | S6 | XL | 45 | 17.18 | 14.39 | ✓ | - | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | ✓ | - | ✓ | A4: domains: i:3-16>i:19-30; A6: domains 1 (label:i:3-16 × loop:i:3-16>i:19-30) |
| 8WH0 | S7 | L | 7 | 6.21 | 2.92 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | ○ |  |
| 6L9Z | S7 | L | 19 | 7.67 | 3.79 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | - | - | ✓ | A6: domains 1 (terminus:N:A × na:run:79) |
| 4NLG | S7 | M | 3 | 3.01 | 0.73 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 2VBO | S7 | M | 4 | 2.66 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | ✓ | - | - | ✓ |  |
| 5K7Z | S7 | M | 6 | 3.31 | 1.08 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ○ | ✗ | - | - | ✓ | A3: default: strand F:1-16; deposited: strand F:1-16; beta: strand F:1-16; domains: strand F:1-16; A15: default 1.106238016241444; A17: pair_f1 0.182 (limit 0.95); stem_recall 0.0 (limit 0.85); fallbacks 1 (limit 0) |
| 2FD8 | S7 | S | 2 | 2.05 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ | A6: default 5 (na:nt:B:501 × strand:A:116-122; na:nt:B:502 × strand:A:204-210; na:nt:B:503 × ligand:AKG:A400); deposited 5 (na:nt:B:501 × strand:A:116-122; na:nt:B:502 × strand:A:204-210; na:nt:B:503 × ligand:AKG:A400); beta 3 (na:nt:B:501 × strand:A:116-122; na:nt:B:502 × strand:A:204-210; na:end:B:501-503:5′ × strand:A:60-65); focus 5 (na:nt:B:501 × strand:A:116-122; na:nt:B:502 × strand:A:204-210; na:nt:B:503 × ligand:AKG:A400) |
| 7NZ2 | S7 | XL | 44 | 74.99 | 55.73 | ✓ | - | ✓ | ✓ | ✗ | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | ✗ | - | - | ✓ | A3: domains: B1:211-221>B1:227-408, B3:211-221>B3:227-408, B3:1051-1123>B3:1136-1140, F3:182-186>F3:197-205; A4: domains: A1:1180-1196>A1:1208-1262, A4:217-221>A4:226-400, A4:986-1045>A4:1051-1124, B1:934-971>B1:983-1044; A6: domains 9 (domain-label:D7 · A1 × loop:B3:211-221>B3:227-408; domain-label:D12 · A1 × loop:A1:1180-1196>A1:1208-1262; domain-label:D1 · A3 × loop:F3:182-186>F3:197-205); A16: domains 6 (); A17: spearman 0.446 (limit 0.6) |
| 9K4F | S8 | L | 4 | 5.11 | 1.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | - | - | - | A17: stems_nearer_frac 0.875 (limit 1.0) |
| 13FZ | S8 | L | 4 | 10.84 | 7.52 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | - | - | ○ | A3: default: strand B:1-13; deposited: strand B:1-13; beta: strand B:1-13; A6: default 3 (na:nt:B:2 × loop:A:1013-1019>A:1027-1031; na:nt:B:3 × strand:A:614-618; na:nt:B:8 × loop:A:1003-1006>A:1013-1019); deposited 3 (na:nt:B:2 × loop:A:1013-1019>A:1027-1031; na:nt:B:3 × strand:A:614-618; na:nt:B:8 × loop:A:1003-1006>A:1013-1019); beta 4 (na:nt:B:2 × loop:A:1013-1019>A:1027-1031; na:nt:B:3 × strand:A:614-618; na:nt:B:8 × loop:A:1003-1006>A:1013-1019); A15: default 17.0586549756951; A17: fallbacks 1 (limit 0) |
| 9WB3 | S8 | M | 4 | 3.01 | 0.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ○ | ✓ | - | - | - |  |
| 7F3J | S8 | M | 8 | 2.81 | 0.19 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 2XLJ | S8 | S | 2 | 3.06 | 0.36 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 6S13 | S8 | XL | 52 | 54.59 | 34.71 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✗ | - | - | - | A3: domains: strand A:2729-2729, strand A:2732-2733, strand A:2737-2738; A6: default 1 (na:tie:448-506:label × na:stem:166); deposited 1 (na:tie:448-506:label × na:stem:166); beta 1 (na:tie:448-506:label × na:stem:166); focus 1 (na:tie:448-506:label × na:stem:166); domains 3 (na:stem:80:first × na:run:239; na:stem:129:first × na:run:239; na:stem:129:first × na:run:240); A17: pair_f1 0.944 (limit 0.95); beads_missing 134 (limit 0); na_crossings 2 (limit 0); bead_overlaps 7 (limit 0) |
| 8ZTU | S8 | XL | 56 | 75.91 | 53.28 | ✓ | - | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✗ | ✗ | - | - | ✓ | A3: domains: strand a:2701-2702, strand a:2705-2705, strand a:2710-2711; A6: default 2 (na:gap:499 × na:run:477; ligand-label:MG:a6201 × na:stem:188); deposited 2 (na:gap:499 × na:run:477; ligand-label:MG:a6201 × na:stem:188); beta 2 (na:gap:499 × na:run:477; ligand-label:MG:a6201 × na:stem:188); focus 2 (na:gap:499 × na:run:477; ligand-label:MG:a6201 × na:stem:188); domains 7 (domain-label:D2 · H × na:run:485; na:gap:499 × na:stem:244; na:gap:499 × na:run:371); A16: domains 1 (); A17: beads_missing 1 (limit 0); bead_overlaps 2 (limit 0); gap_label_collisions 1 (limit 0) |
| 4E8Q | S9 | M | 1 | 4.01 | 1.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 9IWF | S9 | S | 2 | 2.41 | 0.38 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | ○ |  |
| 1Z43 | S9 | S | 1 | 2.31 | 0.38 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | - | - | - |  |
| 9YER | S10 | L | 5 | 6.06 | 4.86 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 9QK1 | S10 | L | 9 | 10.42 | 9.1 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | ✓ | - | - | - | ○ |  |
| 3IB0 | S10 | M | 1 | 3.01 | 2.31 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 6A4Y | S10 | M | 1 | 3.51 | 2.4 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | ○ |  |
| 8JXD | S10 | XL | 2 | 5.01 | 3.87 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | ○ |  |
| 4P39 | S11 | M | 4 | 0.9 | 0.17 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | - |  |
| 3HQB | S11 | S | 2 | 1.26 | 0.35 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | - |  |
| 1C6A | S11 | S | 1 | 1.0 | 0.25 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 244L | S11 | S | 1 | 1.0 | 0.28 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 8U9T | S11 | M | 1 | 1.3 | 0.52 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✗ | ○ | A19: names compared on all 3 pairs; only 3 of 24 and 11 elements paired; A1:11-14 unpaired; A1:17-22 unpaired; A1:35-39 unpaired; A1:48-50 unpaired; A1:54-59 unpaired; A1:63-65 unpaired; A1:67-70 unpaired (+22 more) |
| 7IT6 | S11 | S | 1 | 0.85 | 0.15 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✗ | - | A19: names compared on all 4 pairs; only 4 of 11 and 24 elements paired; A:35-48 unpaired; A:61-63 unpaired; A:68-73 unpaired; A:80-82 unpaired; A:95-97 unpaired; A:103-109 unpaired; A:117-124 unpaired (+20 more) |
| 3U3H | S11 | M | 1 | 3.66 | 2.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | ✓ | ○ |  |
| 6QND | S11 | M | 1 | 1.4 | 0.54 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✓ | ○ |  |
| 2BBJ | S11 | L | 5 | 2.76 | 1.64 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | - | - | ✗ | - | A15: default 1.4176624736566856; A19: names compared on all 13 pairs, 1 differ; η2~η2 (A:52-55~A:51-55); A:193-201 unpaired; A:244-273 unpaired; A:275-299 unpaired; A:302-310 unpaired; A:330-344 unpaired; A:47-49 unpaired; A:161-163 unpaired |
| 2BBH | S11 | S | 1 | 1.25 | 0.39 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | ✗ | ○ | A19: names compared on all 13 pairs, 1 differ; η2~η2 (A:51-55~A:52-55); A:47-49 unpaired; A:161-163 unpaired; A:193-201 unpaired; A:244-273 unpaired; A:275-299 unpaired; A:302-310 unpaired; A:330-344 unpaired |
| AF_AFQ387E4F1 | S12 | M | 1 | 0.9 | 0.26 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFQ8W0Z7F1 | S12 | M | 1 | 1.55 | 0.63 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |
| AF_AFC7J156F1 | S12 | S | 1 | 0.75 | 0.13 | ✓ | - | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | - | - | - | - | - | - |  |

## Two-state pairs (S11)

| pair | seed | partner | sampled as | types now | ligand touching | assembly 1 chains | CA RMSD Å | A19 |
|---|---|---|---|---|---|---|---|---|
| P1 | 3HQB | 4P39 | open/closed | open/closed | [False, False] | [1, 1] | 4.297 | ✓ names compared on all 3 pairs; A:3-11 unpaired |
| P2 | 1C6A | 244L | apo/holo | apo/holo | [False, True] | [1, 1] | 0.219 | ✓ names compared on all 13 pairs |
| P3 | 8U9T | 7IT6 | apo/holo | apo/holo, monomer/oligomer | [True, False] | [2, 1] | 20.979 | ✗ names compared on all 3 pairs; only 3 of 24 and 11 elements paired; A1:11-14 unpaired; A1:17-22 unpaired; A1:35-39 unpaired; A1:48-50 unpaired; A1:54-59 unpaired; A1:63-65 unpaired; A1:67-70 unpaired (+22 more) |
| P4 | 3U3H | 6QND | monomer/oligomer | apo/holo, monomer/oligomer | [True, False] | [4, 1] | 0.439 | ✓ names compared on all 30 pairs; A1:146-148 unpaired |
| P5 | 2BBJ | 2BBH | monomer/oligomer | monomer/oligomer, open/closed | [False, False] | [5, 1] | 3.962 | ✗ names compared on all 13 pairs, 1 differ; η2~η2 (A:52-55~A:51-55); A:193-201 unpaired; A:244-273 unpaired; A:275-299 unpaired; A:302-310 unpaired; A:330-344 unpaired; A:47-49 unpaired; A:161-163 unpaired |

## Reported, not gated

| id | honest_crossings | loop_crossings | na_run_crossings | fallback_curves | fallback_curves_crossing_an_element | kink_links_rejected | tighter_than_inflate | rejected_as_crossing | names_moved | area | max_rss_mb |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1PTY | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 282.83 | 161.9 |
| 6TX2 | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 358.13 | 162.3 |
| 7N71 | 29 | 29 | None | 0 | 0 | 0 | 0 | 0 | 0 | 952.48 | 193.6 |
| 1HMY | 14 | 14 | None | 0 | 0 | 0 | 0 | 0 | 0 | 333.48 | 162.4 |
| 6EJ3 | 26 | 26 | None | 0 | 0 | 0 | 0 | 0 | 0 | 660.62 | 164.4 |
| 5RFU | 16 | 16 | None | 0 | 0 | 0 | 0 | 0 | 0 | 654.74 | 175.9 |
| 2FR2 | 4 | 2 | None | 2 | 2 | 0 | 0 | 0 | 0 | 156.55 | 160.7 |
| 8F03 | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 194.5 | 162.2 |
| 3KSI | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 322.04 | 162.9 |
| 5YUH | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 152.61 | 162.2 |
| 2W6W | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 178.43 | 161.2 |
| 3LIL | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 106.9 | 160.6 |
| 7P46 | 19 | 19 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1187.83 | 209.2 |
| 8HRP | 143 | 143 | None | 0 | 0 | 0 | 2 | 0 | 0 | 2037.26 | 231.5 |
| 6SNU | 100 | 100 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1551.74 | 319.9 |
| 3AMZ | 257 | 257 | None | 0 | 0 | 0 | 1 | 0 | 0 | 3464.26 | 345.4 |
| 8ZT9 | 19 | 19 | None | 0 | 0 | 0 | 0 | 0 | 0 | 725.24 | 177.0 |
| 7X8H | 52 | 52 | None | 0 | 0 | 0 | 2 | 0 | 0 | 875.38 | 196.4 |
| 2NAD | 76 | 76 | None | 2 | 0 | 0 | 0 | 0 | 0 | 730.3 | 176.8 |
| 6B7J | 4 | 4 | None | 0 | 0 | 0 | 0 | 0 | 0 | 524.75 | 163.6 |
| 2O4P | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 381.13 | 161.7 |
| 3B6N | 12 | 12 | None | 0 | 0 | 0 | 1 | 0 | 0 | 616.9 | 163.7 |
| 5EUJ | 20 | 20 | None | 0 | 0 | 0 | 1 | 0 | 0 | 729.18 | 337.8 |
| 9SLB | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1270.52 | 413.9 |
| 3KTK | 10 | 10 | None | 0 | 0 | 0 | 0 | 0 | 0 | 285.14 | 228.1 |
| 5S5Z | 159 | 159 | None | 0 | 0 | 0 | 0 | 0 | 0 | 4315.14 | 360.5 |
| 6HW8 | 358 | 358 | None | 0 | 0 | 0 | 2 | 0 | 0 | 8425.19 | 462.6 |
| 6XER | 128 | 128 | None | 0 | 0 | 0 | 1 | 0 | 0 | 3598.61 | 346.9 |
| 3ULQ | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 646.25 | 165.2 |
| 2UZU | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 346.75 | 163.3 |
| 6N9A | 58 | 58 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1652.02 | 224.7 |
| 2YJA | 8 | 8 | None | 0 | 0 | 0 | 0 | 0 | 0 | 664.18 | 165.3 |
| 4N9F | 15 | 15 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1238.42 | 309.7 |
| 3L4G | 144 | 144 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2497.48 | 311.2 |
| 8FLU | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1666.39 | 218.9 |
| 6ET5 | 37 | 37 | None | 0 | 0 | 0 | 1 | 0 | 1 | 25390.44 | 328.7 |
| 7KH0 | 48 | 48 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1786.7 | 228.7 |
| 8BCU | 297 | 289 | None | 4 | 4 | 0 | 0 | 0 | 0 | 3906.15 | 388.9 |
| 7ZZ1 | 9 | 9 | None | 0 | 0 | 0 | 0 | 0 | 0 | 685.99 | 171.6 |
| 3JBE | 64 | 64 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1076.66 | 174.6 |
| 6WUT | 31 | 31 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1058.6 | 213.3 |
| 8D3D | 45 | 45 | None | 0 | 0 | 0 | 1 | 0 | 0 | 1628.23 | 518.8 |
| 6V6S | 1237 | 1237 | None | 0 | 0 | 0 | 5 | 0 | 0 | 21320.9 | 591.7 |
| 12BM | 597 | 591 | None | 6 | 6 | 0 | 0 | 0 | 0 | 12593.17 | 554.6 |
| 9GTP | 138 | 138 | None | 0 | 0 | 0 | 2 | 0 | 0 | 4545.1 | 483.2 |
| 7BKC | 620 | 620 | None | 0 | 0 | 0 | 7 | 0 | 0 | 13871.11 | 1450.7 |
| 7SPB | 11 | 11 | None | 0 | 0 | 0 | 0 | 0 | 0 | 420.27 | 321.0 |
| 7ZDP | 293 | 293 | None | 0 | 0 | 0 | 3 | 0 | 0 | 15620.27 | 515.0 |
| 9Y05 | 26 | 26 | None | 0 | 0 | 0 | 1 | 0 | 0 | 947.57 | 359.9 |
| 1TKW | 21 | 21 | None | 0 | 0 | 0 | 0 | 0 | 0 | 336.08 | 249.0 |
| 2MRJ | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 67.93 | 190.2 |
| 5GIW | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 8.58 | 159.5 |
| 2LZJ | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 65.53 | 193.8 |
| 2N5J | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 45.29 | 161.6 |
| 8SKX | 1 | 1 | None | 0 | 0 | 0 | 0 | 0 | 0 | 54.42 | 165.6 |
| 1RCH | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 136.41 | 166.5 |
| 9BVV | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 157.3 |
| 9CZH | 42 | 42 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1787.96 | 305.3 |
| 6MJ2 | 60 | 60 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1561.16 | 268.3 |
| 6QFA | 18 | 18 | None | 0 | 0 | 0 | 0 | 0 | 0 | 643.24 | 224.2 |
| 8WMV | 228 | 228 | None | 0 | 0 | 0 | 0 | 0 | 0 | 11009.35 | 444.3 |
| 8JWF | 58 | 58 | None | 0 | 0 | 0 | 3 | 0 | 0 | 1625.4 | 186.4 |
| 7JI3 | 61 | 61 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1219.04 | 193.5 |
| 7VVK | 32 | 32 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1729.5 | 223.8 |
| 6ZKT | 279 | 279 | None | 0 | 0 | 0 | 1 | 0 | 0 | 16192.6 | 486.5 |
| 8WH0 | 245 | 245 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 3254.26 | 372.7 |
| 6L9Z | 138 | 138 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 9369.61 | 396.2 |
| 4NLG | 24 | 24 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 709.03 | 173.2 |
| 2VBO | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 672.1 | 177.4 |
| 5K7Z | 10 | 10 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1080.67 | 198.9 |
| 2FD8 | 11 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 180.42 | 167.5 |
| 7NZ2 | 2240 | 2240 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 38235.07 | 999.9 |
| 9K4F | 44 | 44 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 2245.59 | 199.8 |
| 13FZ | 241 | 241 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 2346.01 | 346.4 |
| 9WB3 | 25 | 25 | 4 | 0 | 0 | 0 | 1 | 0 | 0 | 997.81 | 177.1 |
| 7F3J | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 152.4 | 169.3 |
| 2XLJ | 4 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 272.1 | 167.0 |
| 6S13 | 111 | 111 | 328 | 0 | 0 | 0 | 0 | 0 | 0 | 23431.42 | 736.8 |
| 8ZTU | 293 | 293 | 712 | 0 | 0 | 0 | 2 | 0 | 0 | 35176.31 | 633.4 |
| 4E8Q | 0 | 0 | 10 | 0 | 0 | 0 | 0 | 0 | 0 | 2288.24 | 198.7 |
| 9IWF | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 232.25 | 170.8 |
| 1Z43 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 341.76 | 168.8 |
| 9YER | 99 | 99 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2856.02 | 343.2 |
| 9QK1 | 92 | 92 | None | 0 | 0 | 0 | 0 | 0 | 0 | 3730.69 | 311.8 |
| 3IB0 | 27 | 27 | None | 0 | 0 | 0 | 0 | 0 | 0 | 429.17 | 164.5 |
| 6A4Y | 41 | 41 | None | 0 | 0 | 0 | 1 | 0 | 0 | 558.23 | 166.5 |
| 8JXD | 59 | 59 | None | 0 | 0 | 0 | 0 | 0 | 0 | 1079.91 | 216.4 |
| 4P39 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 63.27 | 159.3 |
| 3HQB | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 66.99 | 158.9 |
| 1C6A | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 190.57 | 161.3 |
| 244L | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 190.57 | 160.8 |
| 8U9T | 20 | 20 | None | 0 | 0 | 0 | 0 | 0 | 0 | 652.69 | 169.7 |
| 7IT6 | 10 | 10 | None | 0 | 0 | 0 | 1 | 0 | 0 | 101.47 | 159.6 |
| 3U3H | 71 | 71 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2235.22 | 256.0 |
| 6QND | 12 | 12 | None | 0 | 0 | 0 | 0 | 0 | 0 | 627.23 | 164.2 |
| 2BBJ | 45 | 45 | None | 0 | 0 | 0 | 0 | 0 | 0 | 2544.46 | 311.8 |
| 2BBH | 2 | 2 | None | 0 | 0 | 0 | 0 | 0 | 0 | 331.67 | 162.2 |
| AF_AFQ387E4F1 | 6 | 6 | None | 0 | 0 | 0 | 0 | 0 | 0 | 373.18 | 162.9 |
| AF_AFQ8W0Z7F1 | 32 | 32 | None | 0 | 0 | 0 | 0 | 0 | 0 | 719.99 | 165.0 |
| AF_AFC7J156F1 | 0 | 0 | None | 0 | 0 | 0 | 0 | 0 | 0 | 113.42 | 159.2 |

## Science checks (S1-S6), apart from A1-A20

Judged in every drawn view (default, deposited, beta, focus, domains, membrane, app); S5 once per entry.

A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet drawn, no barrel, no residue mark, no accessory strand, no CATH reference) is not judged, one with no protein does not apply.

| check | pass | fail | not judged | does not apply |
|---|---|---|---|---|
| S1 sheet neighbours pair | 86 | 0 | 11 | 3 |
| S2 barrel closure pairs | 10 | 0 | 87 | 3 |
| S3 barrel strand count is its pairing cycle | 86 | 0 | 11 | 3 |
| S4 residue marks at their residue | 0 | 0 | 97 | 3 |
| S5 automatic domains agree with a reference | 0 | 0 | 97 | 3 |
| S6 accessory drawn at its partner | 28 | 0 | 69 | 3 |

### Science by view

Entries measured in each view, and how many of them fail each check there.

| view | measured | S1 | S2 | S3 | S4 | S6 |
|---|---|---|---|---|---|---|
| default | 100 | 0 | 0 | 0 | 0 | 0 |
| deposited | 97 | 0 | 0 | 0 | 0 | 0 |
| beta | 97 | 0 | 0 | 0 | 0 | 0 |
| focus | 63 | 0 | 0 | 0 | 0 | 0 |
| domains | 63 | 0 | 0 | 0 | 0 | 0 |
| membrane | 8 | 0 | 0 | 0 | 0 | 0 |
| app | 100 | 0 | 0 | 0 | 0 | 0 |

| id | S1 | S2 | S3 | S4 | S5 | S6 | findings |
|---|---|---|---|---|---|---|---|
| 1PTY | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6TX2 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7N71 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 1HMY | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6EJ3 | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 5RFU | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 2FR2 | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 8F03 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3KSI | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5YUH | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2W6W | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 3LIL | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7P46 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8HRP | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6SNU | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3AMZ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8ZT9 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7X8H | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2NAD | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6B7J | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2O4P | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3B6N | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5EUJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9SLB | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3KTK | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 5S5Z | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6HW8 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6XER | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 3ULQ | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 2UZU | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6N9A | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2YJA | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4N9F | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3L4G | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8FLU | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6ET5 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
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
| 1TKW | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2MRJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5GIW | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 2LZJ | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 2N5J | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 8SKX | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1RCH | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
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
| 4NLG | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2VBO | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 5K7Z | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 2FD8 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 7NZ2 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9K4F | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 13FZ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9WB3 | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 7F3J | ✓ | ✓ | ✓ | ○ | ○ | ○ |  |
| 2XLJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6S13 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8ZTU | ✓ | ✓ | ✓ | ○ | ○ | ✓ |  |
| 4E8Q | - | - | - | - | - | - |  |
| 9IWF | - | - | - | - | - | - |  |
| 1Z43 | - | - | - | - | - | - |  |
| 9YER | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 9QK1 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3IB0 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 6A4Y | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8JXD | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 4P39 | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 3HQB | ○ | ○ | ○ | ○ | ○ | ○ |  |
| 1C6A | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 244L | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 8U9T | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 7IT6 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| 3U3H | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 6QND | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2BBJ | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| 2BBH | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ387E4F1 | ✓ | ○ | ✓ | ○ | ○ | ○ |  |
| AF_AFQ8W0Z7F1 | ✓ | ○ | ✓ | ○ | ○ | ✓ |  |
| AF_AFC7J156F1 | ○ | ○ | ○ | ○ | ○ | ○ |  |

## S5: automatic domains against CATH (reported, not gated)

In-sample: these chains are among the 33 the automatic domains were calibrated on (tests/data/domains/cath_reference.json); the never-seen check is the fresh set (cath_fresh.json).

No entry is in the CATH reference.

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
