# Sensitivity of the tuned constants

Plan 28 Task 9.2 (ruling Q29). Written by: foldmap 0.1.0, git 7f04c895e355 (asserted by runner, dirty False), code fingerprint 4b9c836b01cb8d00, scripts/sensitivity.py sha256 a4f0c262dad3; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Every module constant the inventory gives the origin `measured` or `named-entry` was set to x0.75 and x1.25 of its value, one at a time, and the random-100 campaign keys were measured on the third random sample (none of these constants was tuned on it), default view only, against an unperturbed baseline run the same way. The perturbation reaches every reader (scripts/perturb.py: the module is rewritten as it is imported, so `from ... import` binders, defaults and constants derived from it at import follow; `perturb.py reach` checked each run). Keys: each entry's pass or fail, the failing count per A- and S-criterion, the named fallbacks and the figure area of the default view; time is recorded, not judged, and so is A2 (the default view within 60 s of wall time: it moved between two identical baselines with the load of the machine), so an entry passes here when no other A-criterion fails. **stable**: no key moved; **sensitive**: some key moved; **fragile**: the pass count moved by 2 or more, or a criterion failing in no baseline entry failed in some entry. A constant takes the worse of its two perturbations, from stage B (the whole sample) when it ran there, else from stage A (its subset only).

The baseline draws the default view only, so its pass count is not the third sample's published one (every view drawn, `docs/validation/random100-third/report.md`): criteria needing other views are not judged here.

## Counts

| class | constants |
|---|---|
| stable | 44 |
| sensitive | 11 |
| fragile | 1 |

56 constants perturbed, 7 not perturbable, 1 with both perturbations skipped; zl time 1.0 h (stage A 0.4 h, stage B 0.6 h).

Of the stable constants, 25 were read by no line that ran in the default view (the reads run: a line inside a function naming the constant, perturb.install_trace), so their stability here is no evidence either way: `align:SAME_FOLD_TM`, `align:SAME_FOLD_FRACTION`, `cli:MIN_TAU_GAIN_INV`, `io:PAIRED_MIN`, `io:_PAIR_DIST`, `layout:_SEP_TIE_ROOM`, `layout:_POCKET_REACH`, `layout:_PART_BLOWUP`, `layout:_NAME_OPEN`, `membrane:_MAX_POLAR_CROSSED`, `membrane:_MANY_CROSSING`, `membrane:_FLANK`, `membrane:_FLANK_REACH`, `membrane:_MIN_OUTWARD_HYDROPHOBIC`, `membrane:_MAX_OUTWARD_CHARGED`, `membrane:_LOOP_CAP`, `membrane:_SIDED_LOOPS`, `membrane:MIN_CROSSING_HELICES`, `membrane:FACE_FLAT_TILT`, `membrane:JOIN_GAP`, `na_runs:_JOINT_REACH`, `na_runs:_BESIDE_TOLL`, `panels:CLOSE`, `route:_NA_TEXT_TOLL`, `route:_REROUTE_TURNS_INK`.

## Limitations: the fragile constants

- `layout:_SIDE_CLEAR` = 2.75 (measured), stage B: x0.75 (2.0625): sensitive; pass -1; A17 7->9; counts_hidden 15->20; helix_name_moved 41->46; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->3; area 8 entries; 8 entries moved; x1.25 (3.4375): **fragile**; pass -5; A17 7->13; A3 2->3; A6 4->3; counts_hidden 15->21; helix_name_moved 41->39; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->2; area 8 entries; 8 entries moved

## Stage A: the subset

41 entries of the third sample (seed 20261025, target 40; per stratum in proportion to its eligible entries, largest remainder, a pair one unit): 1RKA, 5ES7, 4GWY, 3CYR, 5JNN, 3O6N, 3P2W, 8CRC, 2YEJ, 3FT8, 5U3R, 5XMX, AF_AFQ5AGV5F1, 1VSK, 1GQ2, 3WAU, 4YP6, 4YJ5, 3ABR, 8PIN, 2YNR, 4ODQ, 9N2E, 7U9Z, 9BZD, 7SPI, 7XBK, 2MXC, 2KQV, 5NF8, 1UG2, 6VBG, 4D5B, 6BO8, 6QUE, 3W96, 6VTX, 8VBG, 3EGZ, 6JCY, 9OWM.

Left out, recorded slower than 120 s in the third run (every view): 8XJV (1092 s), 9FZL (706 s), 7F5S (483 s), 8EM7 (388 s), 6QIK (380 s), 7KMZ (357 s), 3G6E (291 s), 6ZTQ (271 s), 8FZN (146 s), 9NMN (121 s).

Stage A: 41 entries, baseline 32 pass; failing per criterion A17 3, A19 6, S5 6. The two baselines are identical (keys and SVGs). A2 (time, not judged) failed in none and in none of the two baselines. 111 perturbation runs; 0.42 h of zl time. The reads run (unperturbed, traced) changed no key.

Stage B: 100 entries, baseline 80 pass; failing per criterion A15 1, A17 7, A19 10, A3 2, A4 1, A6 4, S5 12. The two baselines are identical (keys and SVGs). A2 (time, not judged) failed in 8EM7 and in 7KMZ, 8EM7 of the two baselines. 24 perturbation runs; 0.62 h of zl time. The reads run (unperturbed, traced) changed no key.

## Per constant

read in: the entries of the deciding stage a line naming the constant ran in (default view), + when the judging read it.

| constant | origin | value | read in | x0.75 | x1.25 | class |
|---|---|---|---|---|---|---|
| `align:SAME_FOLD_TM` | measured | 0.5 | 0/41 | 0.375: A stable | 0.625: A stable | stable (A) |
| `align:SAME_FOLD_FRACTION` | measured | 0.5 | 0/41 | 0.375: A stable | 0.625: A stable | stable (A) |
| `cli:MIN_TAU_GAIN_INV` | measured | 10 | 0/41 | 8: A stable | 13: A stable | stable (A) |
| `features:DOMAIN_CONTACT` | measured | 12.0 | 35/100 | 9.0: A sensitive; S5 6->9; 3 entries moved; B sensitive; S5 12->15; 5 entries moved | 15.0: A sensitive; S5 6->4; 2 entries moved; B sensitive; S5 12->9; 3 entries moved | sensitive (B) |
| `features:DOMAIN_SEPARATION` | measured | 6 | 35/100 | 5: A sensitive; S5 6->5; 1 entry moved; B sensitive; S5 12->11; 1 entry moved | 8: A stable; B sensitive; S5 12->11; 1 entry moved | sensitive (B) |
| `features:DOMAIN_NCUT` | measured | 0.37 | 35/100 | 0.2775: A sensitive; 2 entries moved; B sensitive; S5 12->11; 3 entries moved | 0.4625: A sensitive; S5 6->12; 6 entries moved; B sensitive; S5 12->19; 7 entries moved | sensitive (B) |
| `features:DOMAIN_MIN_RESIDUES` | measured | 40 | 35/100 | 30: A sensitive; S5 6->5; 1 entry moved; B sensitive; S5 12->11; 1 entry moved | 50: A sensitive; S5 6->5; 1 entry moved; B sensitive; S5 12->11; 1 entry moved | sensitive (B) |
| `features:DOMAIN_MIN_SEGMENT` | measured | 20 | 22/41 | 15: A stable | 25: A stable | stable (A) |
| `figcheck:STEM_GAP_FLOOR` | measured | 0.08 | 6/41 | 0.06: A stable | 0.1: A stable | stable (A) |
| `io:PAIRED_MIN` | named-entry:4PLX,3OWW | 3 | 0/41 | 2: A stable | 4: A stable | stable (A) |
| `io:_PAIR_DIST` | named-entry:3OWW | 3.4 | 0/41 | 2.55: A stable | 4.25: A stable | stable (A) |
| `layout:_TAIL_COUNT_NEAR` | named-entry:7SMR | 1.5 | 13/41 | 1.125: A stable | 1.875: A stable | stable (A) |
| `layout:_SS_RESTART_PAIRS` | measured | 5000000000 | 2/41 | 3750000000: A stable | 6250000000: A stable | stable (A) |
| `layout:_K_LIKE` | measured | 2.56 | 41/41 | 1.92: A stable | 3.2: A stable | stable (A) |
| `layout:_SEP_MIRROR_ROOM` | named-entry:9GEP | 1.1 | 10/100 | 0.825: A sensitive; area 1 entry; 1 entry moved; B sensitive; area 1 entry; 1 entry moved | 1.375: A stable; B stable | sensitive (B) |
| `layout:_SEP_TIE_ROOM` | named-entry:8QEM,6JLN,7QUP,9P1X | 2.0 | 0/41 | 1.5: A stable | 2.5: A stable | stable (A) |
| `layout:_STACK_PAD` | measured | 1.3 | 37/100 | 0.975: A sensitive; area 5 entries; 5 entries moved; B sensitive; A17 7->9; A4 1->0; counts_hidden 15->14; helix_name_moved 41->39; loop_over_text 1->0; na_run_unclear 1->2; area 27 entries; 27 entries moved | 1.625: A sensitive; area 5 entries; 5 entries moved; B sensitive; A17 7->9; A3 2->3; A4 1->0; A6 4->3; counts_hidden 15->18; helix_name_moved 41->34; loop_over_text 1->0; na_run_unclear 1->2; area 26 entries; 26 entries moved | sensitive (B) |
| `layout:_STACK_SHUT` | named-entry:9P1X,8OZJ | 0.6 | 4/41 | 0.45: A stable | 0.75: A stable | stable (A) |
| `layout:_STACK_NUDGE` | named-entry:9P1X,8OZJ | 1.0 | 4/41 | 0.75: A stable | 1.25: A stable | stable (A) |
| `layout:_POCKET_REACH` | named-entry:8BCU,12BM | 2.0 | 0/41 | 1.5: A stable | 2.5: A stable | stable (A) |
| `layout:_END_ROOM` | named-entry:9FWZ | 0.9999999999999999 | 25/100 | 0.75: A sensitive; area 1 entry; 1 entry moved; B sensitive; pass -1; A15 1->2; helix_name_moved 41->19; area 5 entries; 5 entries moved | 1.25: A sensitive; area 1 entry; 1 entry moved; B sensitive; counts_hidden 15->16; helix_name_moved 41->32; area 5 entries; 5 entries moved | sensitive (B) |
| `layout:_BA_MAX` | named-entry:2BNH | 10 | 36/41 | 8: A stable | 13: A stable | stable (A) |
| `layout:_STEM_ROOM` | measured | 1.5 | 3/100 | 1.125: A **fragile**; pass -1; new failing A6; area 1 entry; 1 entry moved; B sensitive; pass -1; A6 4->5; area 2 entries; 2 entries moved | 1.875: A sensitive; area 1 entry; 1 entry moved; B sensitive; area 2 entries; 2 entries moved | sensitive (B) |
| `layout:_SIDE_CLEAR` | measured | 2.75 | 16/100 | 2.0625: A sensitive; area 2 entries; 2 entries moved; B sensitive; pass -1; A17 7->9; counts_hidden 15->20; helix_name_moved 41->46; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->3; area 8 entries; 8 entries moved | 3.4375: A **fragile**; pass -2; A17 3->5; area 2 entries; 2 entries moved; B **fragile**; pass -5; A17 7->13; A3 2->3; A6 4->3; counts_hidden 15->21; helix_name_moved 41->39; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->2; area 8 entries; 8 entries moved | fragile (B) |
| `layout:_HELIX_LOOP` | measured | 3 | 6/100 | 2: A sensitive; area 1 entry; 1 entry moved; B sensitive; pass -1; A17 7->9; A6 4->3; counts_hidden 15->17; helix_name_moved 41->37; na_gap_unclear 3->4; na_tail_straight 4->5; area 6 entries; 6 entries moved | 4: A stable; B sensitive; pass -1; A17 7->9; A3 2->4; counts_hidden 15->23; helix_name_moved 41->49; na_gap_unclear 3->7; na_run_unclear 1->3; na_tail_straight 4->5; route_second_search_refused 0->5; area 4 entries; 4 entries moved | sensitive (B) |
| `layout:_PART_BLOWUP` | named-entry:1WWL,2ILU,7SMR,6JLN,5LDW | 2.0 | 0/41 | 1.5: A stable | 2.5: A stable | stable (A) |
| `layout:_NAME_OPEN` | named-entry:7NZ2 | 0.25 | 0/41 | 0.1875: A stable | 0.3125: A stable | stable (A) |
| `membrane:_MAX_POLAR_CROSSED` | named-entry:6JO6 | 0.2 | 0/41 | 0.15: A stable | 0.25: A stable | stable (A) |
| `membrane:_MANY_CROSSING` | named-entry:6JO6 | 10 | 0/41 | 8: A stable | 13: A stable | stable (A) |
| `membrane:_FLANK` | measured | 20 | 0/41 | 15: A stable | 25: A stable | stable (A) |
| `membrane:_FLANK_REACH` | measured | 15.0 | 0/41 | 11.25: A stable | 18.75: A stable | stable (A) |
| `membrane:_MIN_OUTWARD_HYDROPHOBIC` | measured | 0.6 | 0/41 | 0.45: A stable | 0.75: A stable | stable (A) |
| `membrane:_MAX_OUTWARD_CHARGED` | measured | 0.15 | 0/41 | 0.1125: A stable | 0.1875: A stable | stable (A) |
| `membrane:_LOOP_CAP` | measured | 20 | 0/41 | 15: A stable | 25: A stable | stable (A) |
| `membrane:_SIDED_LOOPS` | measured | 10 | 0/41 | 8: A stable | 13: A stable | stable (A) |
| `membrane:MIN_CROSSING_HELICES` | named-entry:7PIU,7VVK | 3 | 0/41 | 2: A stable | 4: A stable | stable (A) |
| `membrane:FACE_FLAT_TILT` | named-entry:32KS | 35.0 | 0/41 | 26.25: A stable | 43.75: A stable | stable (A) |
| `membrane:JOIN_GAP` | named-entry:32KS,6UW4 | 2 | 0/41 | skipped | 3: A stable | stable (A) |
| `na_runs:_JOINT_NEAR` | named-entry:9C75 | 2.5 | 2/41 | 1.875: A stable | 3.125: A stable | stable (A) |
| `na_runs:_JOINT_REACH` | named-entry:8VWQ | 5.0 | 0/41 | 3.75: A stable | 6.25: A stable | stable (A) |
| `na_runs:_BESIDE_TOLL` | named-entry:7Y7C | 50.0 | 0/41 | 37.5: A stable | 62.5: A stable | stable (A) |
| `na_runs:_DETOUR` | measured | 3.0 | 4/41 | 2.25: A stable | 3.75: A stable | stable (A) |
| `na_runs:_DETOUR_STEPS` | measured | 6.0 | 4/41 | 4.5: A stable | 7.5: A stable | stable (A) |
| `na_runs:_FIRST_NODES` | measured | 5000 | 3/41 | 3750: A stable | 6250: A stable | stable (A) |
| `na_runs:_BESIDE_DETOUR` | measured | 2.0 | 12/100 | 1.5: A stable; B sensitive; counts_hidden 15->14; 1 entry moved | 2.5: A sensitive; area 1 entry; 1 entry moved; B sensitive; area 2 entries; 2 entries moved | sensitive (B) |
| `nucleic:HBOND_STRETCH` | named-entry:1HYS,1AOI | 4.0 | 16/100 | 3.0: A sensitive; area 1 entry; 1 entry moved; B sensitive; A17 7->8; A6 4->3; counts_hidden 15->17; helix_name_moved 41->36; na_gap_unclear 3->5; na_tail_straight 4->5; page_text_small 2->3; area 7 entries; 7 entries moved | 5.0: A sensitive; area 2 entries; 2 entries moved; B sensitive; pass -1; A17 7->9; A6 4->5; counts_hidden 15->22; na_gap_unclear 3->4; page_text_small 2->3; area 8 entries; 8 entries moved | sensitive (B) |
| `panels:CLOSE` | named-entry:2LZM | 2.0 | 0/41 | 1.5: A stable | 2.5: A stable | stable (A) |
| `render:_SPOT_FAR` | named-entry:7QUP | 3 | 1/41 | 2: A stable | 4: A stable | stable (A) |
| `render:_SPOT_MORE_RINGS` | named-entry:7LX9,9YLM,3MHT,1GOX | 4 | 21/41 | 3: A stable | 5: A stable | stable (A) |
| `route:_SLAB_RUN` | named-entry:7SMR | 12.0 | 41/41 | 9.0: A stable | 15.0: A stable | stable (A) |
| `route:_MAX_NODES_WIDE` | measured | 1600000 | 41/41 | 1200000: A stable | 2000000: A stable | stable (A) |
| `route:_MAX_NODES_PAGE` | named-entry:7QUP | 8000000 | 41/41 | 6000000: A stable | 10000000: A stable | stable (A) |
| `route:_MAX_NODES_WINDOWS` | named-entry:7NZ2,9FZL | 8000000 | 41/41 | 6000000: A stable | 10000000: A stable | stable (A) |
| `route:_NA_TEXT_TOLL` | named-entry:6SO9,9CJI,7XSP | 25.0 | 0/41 | 18.75: A stable | 31.25: A stable | stable (A) |
| `route:_REROUTE_TURNS_INK` | named-entry:7XSP | 12 | 0/41 | 9: A stable | 15: A stable | stable (A) |
| `symmetry:MIN_HELIX_COPIES` | named-entry:1URN | 4 | 25/41 | 3: A stable | 5: A stable | stable (A) |

## Skipped perturbations

- `align:SCREEN` = 0.0 x0.75: x0.75 gives the same value (0.0)
- `align:SCREEN` = 0.0 x1.25: x1.25 gives the same value (0.0)
- `membrane:JOIN_GAP` = 2 x0.75: x0.75 gives the same value (2)

## Not perturbable

- `cli:ROSETTE_ORDERS` = (4, 7): a tuple: no single factor
- `na_ends:_OFF_TURNS` = (np.float64(0.0), np.float64(0.08726646259971647), np.float64(0.174532925199432…: a tuple: no single factor
- `na_labels:_GAP_OUT` = (1.0, 2.0, 3.0, 4.0, 6.0): a tuple: no single factor
- `na_labels:_GAP_TURNS` = (0.0, 0.08726646259971647, 0.17453292519943295, 0.2617993877991494, 0.349065850…: a tuple: no single factor
- `na_runs:_JOINT_BENDS` = (1.0, 1.5, 2.0, 0.5, 3.0): a tuple: no single factor
- `na_runs:_JOINT_SWINGS` = (0.0, 1.0, 2.0, -1.0, -2.0): a tuple: no single factor
- `route:_NAME_ALONG` = (0.0, -0.25, 0.25, -0.4, 0.4): a tuple: no single factor
