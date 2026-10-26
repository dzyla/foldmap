# Sensitivity of the tuned constants

Plan 28 Task 9.2 (rulings Q29, Q30). Written by: foldmap 0.1.0, git aba1dd66c1a0 (asserted by runner, dirty False), code fingerprint 4b9c836b01cb8d00, scripts/sensitivity.py sha256 2d5a2976de8d; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Every module constant the inventory gives the origin `measured` or `named-entry` was set to x0.75 and x1.25 of its value, one at a time, and the random-100 campaign keys were measured on the third random sample, in the default view, against an unperturbed baseline run the same way: stage A on a subset, stage B on the whole sample for every constant. The third sample is not the tuning set (the first sample), but it has been run and looked at since plan 15, and some constants were set on its entries since (listed below: `named in the sample`); it is no unseen sample. The perturbation reaches every reader (scripts/perturb.py: the module is rewritten as it is imported, so `from ... import` binders, defaults and constants derived from it at import follow; `perturb.py reach` checked each run). Keys: each entry's pass or fail, the failing count per A- and S-criterion, the named fallbacks and the figure area of the view. **stable**: no key moved; **sensitive**: some key moved; **fragile**: the pass count moved by 2 or more, or a criterion failing in no baseline entry failed in some entry. A constant takes the worse of its two perturbations.

A reads run (unperturbed, traced: perturb.install_trace) records which entries ran a line, inside a function, that names each constant. A constant stable in the default view that no line of it read was measured again in the view that reads it (stage V: its entries, the default view drawn too, keyed on that view), found by a traced run of every view on the whole sample (stage X); **not exercised** when no measurable view reads it, or its stage V reads run saw no read: its stability would be no evidence.

**Time is recorded, not judged, and so is A2** (the default view drawn within 60 s of wall time, a criterion of random100's pass): it moved between two runs of one code on one machine with nothing else changed (see stage B below), so here an entry passes when no other A-criterion fails. A pass count random100 reports includes this load-dependent criterion.

The baselines draw the default view (and stage V's view) only, so their pass counts are not the third sample's published ones (every view drawn, `docs/validation/random100-third/report.md`).

**Limits of the method.** One factor at a time: no interaction between constants is measured. Two points per constant (x0.75, x1.25): a constant can move at another value and not at these. One sample (the third, 100 entries): a class is this sample's. Stage A's subset leaves out the slowest entries, where the route search budgets (`route:_MAX_NODES_*`) bind; stage B, on the whole sample, does not. A read is a line naming the constant that ran, not a decision it made. A run that crashed or whose value did not reach its readers is incomplete and gets no class.

**Named in the sample**: the constants whose inventory origin names an entry of the third sample: `membrane:JOIN_GAP` (6UW4), `route:_MAX_NODES_WINDOWS` (9FZL).

## Counts

| class | constants |
|---|---|
| fragile | 1 |
| sensitive | 21 |
| stable | 29 |
| not exercised | 5 |

56 constants perturbed, 7 not perturbable, 1 with both perturbations skipped; zl time 3.5 h (stage A 0.42 h, stage B 2.64 h, stage V (membrane) 0.04 h, stage V (domains) 0.45 h).

## Limitations: the fragile constants

- `layout:_SIDE_CLEAR` = 2.75 (measured), stage B: x0.75 (2.0625): sensitive; pass -1; A17 7->9; counts_hidden 15->20; helix_name_moved 41->46; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->3; area 8 entries; 8 entries moved; x1.25 (3.4375): **fragile**; pass -5; A17 7->13; A3 2->3; A6 4->3; counts_hidden 15->21; helix_name_moved 41->39; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->2; area 8 entries; 8 entries moved

## Not exercised

- `align:SAME_FOLD_TM` = 0.5 (measured): read by no view of random100 measured here (survey: no process read it)
- `align:SAME_FOLD_FRACTION` = 0.5 (measured): read by no view of random100 measured here (survey: no process read it)
- `io:PAIRED_MIN` = 3 (named-entry:4PLX,3OWW): read by no view of random100 measured here (survey: no process read it)
- `io:_PAIR_DIST` = 3.4 (named-entry:3OWW): read by no view of random100 measured here (survey: no process read it)
- `layout:_POCKET_REACH` = 2.0 (named-entry:8BCU,12BM): read by no view of random100 measured here (survey: no process read it)

## Stages

Stage A: 41 entries of the third sample (seed 20261025, target 40 entries; the quota per stratum counts entries, in proportion to its eligible entries, largest remainder; units are drawn whole, a pair being one unit): 1RKA, 5ES7, 4GWY, 3CYR, 5JNN, 3O6N, 3P2W, 8CRC, 2YEJ, 3FT8, 5U3R, 5XMX, AF_AFQ5AGV5F1, 1VSK, 1GQ2, 3WAU, 4YP6, 4YJ5, 3ABR, 8PIN, 2YNR, 4ODQ, 9N2E, 7U9Z, 9BZD, 7SPI, 7XBK, 2MXC, 2KQV, 5NF8, 1UG2, 6VBG, 4D5B, 6BO8, 6QUE, 3W96, 6VTX, 8VBG, 3EGZ, 6JCY, 9OWM. Left out, recorded slower than 120 s in the third run (every view), or the partner of such an entry: 8XJV (1092 s), 9FZL (706 s), 7F5S (483 s), 8EM7 (388 s), 6QIK (380 s), 7KMZ (357 s), 3G6E (291 s), 6ZTQ (271 s), 8FZN (146 s), 9NMN (121 s).

Stage A (default, keyed on the default view): 41 entries, baseline 32 pass; failing per criterion A17 3, A19 6, S5 6. The two baselines are identical (keys and SVGs). 111 perturbation runs; 0.42 h of zl time. The reads run (unperturbed, traced) changed no key. With A2 judged, as random100 judges it, the unperturbed runs of this one code pass baseline-1 32, baseline-2 32, reads 32 of 41 (A2 failed in none), their keys and SVGs otherwise identical.

Stage B (default, keyed on the default view): 100 entries, baseline 80 pass; failing per criterion A15 1, A17 7, A19 10, A3 2, A4 1, A6 4, S5 12. The two baselines are identical (keys and SVGs). 111 perturbation runs; 2.64 h of zl time. The reads run (unperturbed, traced) changed no key. With A2 judged, as random100 judges it, the unperturbed runs of this one code pass baseline-1 79, baseline-2 78, reads 80 of 100 (baseline-1: 8EM7 68.53 s (limit 60 s); baseline-2: 7KMZ 61.58 s (limit 60 s), 8EM7 70.39 s (limit 60 s)), their keys and SVGs otherwise identical.

Stage V (membrane) (default,membrane, keyed on the membrane view): 8 entries, baseline 7 pass; failing per criterion A18 1. The two baselines are identical (keys and SVGs). 23 perturbation runs; 0.04 h of zl time. The reads run (unperturbed, traced) changed no key. With A2 judged, as random100 judges it, the unperturbed runs of this one code pass baseline-1 7, baseline-2 7, reads 7 of 8 (A2 failed in none), their keys and SVGs otherwise identical.

Stage V (domains) (default,domains, keyed on the domains view): 49 entries, baseline 42 pass; failing per criterion A16 1, A17 3, A3 3, A4 2, A6 5, S5 9. The two baselines are identical (keys and SVGs). 6 perturbation runs; 0.45 h of zl time. The reads run (unperturbed, traced) changed no key. With A2 judged, as random100 judges it, the unperturbed runs of this one code pass baseline-1 40, baseline-2 40, reads 40 of 49 (baseline-1: 7KMZ 66.7 s (limit 60 s), 8EM7 86.06 s (limit 60 s); baseline-2: 7KMZ 64.7 s (limit 60 s), 8EM7 83.26 s (limit 60 s); reads: 7KMZ 74.49 s (limit 60 s), 8EM7 69.1 s (limit 60 s)), their keys and SVGs otherwise identical.

Stage V (membrane) measured `cli:MIN_TAU_GAIN_INV`, `membrane:FACE_FLAT_TILT`, `membrane:JOIN_GAP`, `membrane:MIN_CROSSING_HELICES`, `membrane:_FLANK`, `membrane:_FLANK_REACH`, `membrane:_LOOP_CAP`, `membrane:_MANY_CROSSING`, `membrane:_MAX_OUTWARD_CHARGED`, `membrane:_MAX_POLAR_CROSSED`, `membrane:_MIN_OUTWARD_HYDROPHOBIC`, `membrane:_SIDED_LOOPS` on 4D5B, 6BO8, 6DQN, 6QUE, 6UW4, 6VBG, 7Y63, 8GG4.

Stage V (domains) measured `layout:_NAME_OPEN`, `layout:_PART_BLOWUP`, `panels:CLOSE` on 10LI, 12DV, 1NID, 1RKA, 1XZ4, 3ABR, 3EOM, 3G6E, 3W96, 4B2T, 4GWY, 4YJ5, 4ZH4, 5C1A, 5NF8, 5S4U, 6BO8, 6DQN, 6JCY, 6LAB, 6QIK, 6QUE, 6UW4, 6VBG, 6ZTQ, 6ZY7, 7ADZ, 7F1V, 7F5S, 7KMZ, 7SPI, 7U9Z, 7UKO, 7XBK, 7Y63, 8EM7, 8FZN, 8GG4, 8PIN, 8U1N, 8VBG, 9BZD, 9DU3, 9FZL, 9HEO, 9N2E, 9NMN, 9XSF, 9Y70.

## Per constant

read in: the entries of the deciding stage (its view) a line naming the constant ran in, + when the judging read it. Each cell: the perturbed value, then each stage's class and what moved.

| constant | origin | value | read in | x0.75 | x1.25 | class (stage) |
|---|---|---|---|---|---|---|
| `align:SAME_FOLD_TM` | measured | 0.5 | 0/100 | 0.375: A stable; B stable | 0.625: A stable; B stable | not exercised (B) |
| `align:SAME_FOLD_FRACTION` | measured | 0.5 | 0/100 | 0.375: A stable; B stable | 0.625: A stable; B stable | not exercised (B) |
| `cli:MIN_TAU_GAIN_INV` | measured | 10 | 4/8 | 8: A stable; B stable; membrane sensitive; area 1 entry; 1 entry moved | 13: A stable; B stable; membrane stable | sensitive (membrane) |
| `features:DOMAIN_CONTACT` | measured | 12.0 | 35/100 | 9.0: A sensitive; S5 6->9; 3 entries moved; B sensitive; S5 12->15; 5 entries moved | 15.0: A sensitive; S5 6->4; 2 entries moved; B sensitive; S5 12->9; 3 entries moved | sensitive (B) |
| `features:DOMAIN_SEPARATION` | measured | 6 | 35/100 | 5: A sensitive; S5 6->5; 1 entry moved; B sensitive; S5 12->11; 1 entry moved | 8: A stable; B sensitive; S5 12->11; 1 entry moved | sensitive (B) |
| `features:DOMAIN_NCUT` | measured | 0.37 | 35/100 | 0.2775: A sensitive; 2 entries moved; B sensitive; S5 12->11; 3 entries moved | 0.4625: A sensitive; S5 6->12; 6 entries moved; B sensitive; S5 12->19; 7 entries moved | sensitive (B) |
| `features:DOMAIN_MIN_RESIDUES` | measured | 40 | 35/100 | 30: A sensitive; S5 6->5; 1 entry moved; B sensitive; S5 12->11; 1 entry moved | 50: A sensitive; S5 6->5; 1 entry moved; B sensitive; S5 12->11; 1 entry moved | sensitive (B) |
| `features:DOMAIN_MIN_SEGMENT` | measured | 20 | 35/100 | 15: A stable; B stable | 25: A stable; B sensitive; S5 12->11; 1 entry moved | sensitive (B) |
| `figcheck:STEM_GAP_FLOOR` | measured | 0.08 | 16/100 | 0.06: A stable; B stable | 0.1: A stable; B stable | stable (B) |
| `io:PAIRED_MIN` | named-entry:4PLX,3OWW | 3 | 0/100 | 2: A stable; B stable | 4: A stable; B stable | not exercised (B) |
| `io:_PAIR_DIST` | named-entry:3OWW | 3.4 | 0/100 | 2.55: A stable; B stable | 4.25: A stable; B stable | not exercised (B) |
| `layout:_TAIL_COUNT_NEAR` | named-entry:7SMR | 1.5 | 38/100 | 1.125: A stable; B stable | 1.875: A stable; B stable | stable (B) |
| `layout:_SS_RESTART_PAIRS` | measured | 5000000000 | 17/100 | 3750000000: A stable; B stable | 6250000000: A stable; B sensitive; area 1 entry; 1 entry moved | sensitive (B) |
| `layout:_K_LIKE` | measured | 2.56 | 100/100 | 1.92: A stable; B stable | 3.2: A stable; B stable | stable (B) |
| `layout:_SEP_MIRROR_ROOM` | named-entry:9GEP | 1.1 | 10/100 | 0.825: A sensitive; area 1 entry; 1 entry moved; B sensitive; area 1 entry; 1 entry moved | 1.375: A stable; B stable | sensitive (B) |
| `layout:_SEP_TIE_ROOM` | named-entry:8QEM,6JLN,7QUP,9P1X | 2.0 | 2/100 | 1.5: A stable; B sensitive; pass -1; A15 1->2; counts_hidden 15->16; helix_name_moved 41->19; area 1 entry; 1 entry moved | 2.5: A stable; B stable | sensitive (B) |
| `layout:_STACK_PAD` | measured | 1.3 | 37/100 | 0.975: A sensitive; area 5 entries; 5 entries moved; B sensitive; A17 7->9; A4 1->0; counts_hidden 15->14; helix_name_moved 41->39; loop_over_text 1->0; na_run_unclear 1->2; area 27 entries; 27 entries moved | 1.625: A sensitive; area 5 entries; 5 entries moved; B sensitive; A17 7->9; A3 2->3; A4 1->0; A6 4->3; counts_hidden 15->18; helix_name_moved 41->34; loop_over_text 1->0; na_run_unclear 1->2; area 26 entries; 26 entries moved | sensitive (B) |
| `layout:_STACK_SHUT` | named-entry:9P1X,8OZJ | 0.6 | 25/100 | 0.45: A stable; B sensitive; area 1 entry; 1 entry moved | 0.75: A stable; B stable | sensitive (B) |
| `layout:_STACK_NUDGE` | named-entry:9P1X,8OZJ | 1.0 | 25/100 | 0.75: A stable; B stable | 1.25: A stable; B stable | stable (B) |
| `layout:_POCKET_REACH` | named-entry:8BCU,12BM | 2.0 | 0/100 | 1.5: A stable; B stable | 2.5: A stable; B stable | not exercised (B) |
| `layout:_END_ROOM` | named-entry:9FWZ | 0.9999999999999999 | 25/100 | 0.75: A sensitive; area 1 entry; 1 entry moved; B sensitive; pass -1; A15 1->2; helix_name_moved 41->19; area 5 entries; 5 entries moved | 1.25: A sensitive; area 1 entry; 1 entry moved; B sensitive; counts_hidden 15->16; helix_name_moved 41->32; area 5 entries; 5 entries moved | sensitive (B) |
| `layout:_BA_MAX` | named-entry:2BNH | 10 | 88/100 | 8: A stable; B stable | 13: A stable; B stable | stable (B) |
| `layout:_STEM_ROOM` | measured | 1.5 | 3/100 | 1.125: A **fragile**; pass -1; new failing A6; area 1 entry; 1 entry moved; B sensitive; pass -1; A6 4->5; area 2 entries; 2 entries moved | 1.875: A sensitive; area 1 entry; 1 entry moved; B sensitive; area 2 entries; 2 entries moved | sensitive (B) |
| `layout:_SIDE_CLEAR` | measured | 2.75 | 16/100 | 2.0625: A sensitive; area 2 entries; 2 entries moved; B sensitive; pass -1; A17 7->9; counts_hidden 15->20; helix_name_moved 41->46; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->3; area 8 entries; 8 entries moved | 3.4375: A **fragile**; pass -2; A17 3->5; area 2 entries; 2 entries moved; B **fragile**; pass -5; A17 7->13; A3 2->3; A6 4->3; counts_hidden 15->21; helix_name_moved 41->39; na_gap_unclear 3->5; na_run_unclear 1->2; na_tail_straight 4->2; area 8 entries; 8 entries moved | fragile (B) |
| `layout:_HELIX_LOOP` | measured | 3 | 6/100 | 2: A sensitive; area 1 entry; 1 entry moved; B sensitive; pass -1; A17 7->9; A6 4->3; counts_hidden 15->17; helix_name_moved 41->37; na_gap_unclear 3->4; na_tail_straight 4->5; area 6 entries; 6 entries moved | 4: A stable; B sensitive; pass -1; A17 7->9; A3 2->4; counts_hidden 15->23; helix_name_moved 41->49; na_gap_unclear 3->7; na_run_unclear 1->3; na_tail_straight 4->5; route_second_search_refused 0->5; area 4 entries; 4 entries moved | sensitive (B) |
| `layout:_PART_BLOWUP` | named-entry:1WWL,2ILU,7SMR,6JLN,5LDW | 2.0 | 49/49 | 1.5: A stable; B stable; domains sensitive; area 3 entries; 3 entries moved | 2.5: A stable; B stable; domains stable | sensitive (domains) |
| `layout:_NAME_OPEN` | named-entry:7NZ2 | 0.25 | 49/49 | 0.1875: A stable; B stable; domains stable | 0.3125: A stable; B stable; domains sensitive; domain_name_least_bad 0->1; 1 entry moved | sensitive (domains) |
| `membrane:_MAX_POLAR_CROSSED` | named-entry:6JO6 | 0.2 | 7/8 | 0.15: A stable; B stable; membrane stable | 0.25: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_MANY_CROSSING` | named-entry:6JO6 | 10 | 7/8 | 8: A stable; B stable; membrane stable | 13: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_FLANK` | measured | 20 | 7/8 | 15: A stable; B stable; membrane stable | 25: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_FLANK_REACH` | measured | 15.0 | 7/8 | 11.25: A stable; B stable; membrane stable | 18.75: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_MIN_OUTWARD_HYDROPHOBIC` | measured | 0.6 | 1/8 | 0.45: A stable; B stable; membrane stable | 0.75: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_MAX_OUTWARD_CHARGED` | measured | 0.15 | 1/8 | 0.1125: A stable; B stable; membrane stable | 0.1875: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_LOOP_CAP` | measured | 20 | 1/8 | 15: A stable; B stable; membrane stable | 25: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:_SIDED_LOOPS` | measured | 10 | 1/8 | 8: A stable; B stable; membrane stable | 13: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:MIN_CROSSING_HELICES` | named-entry:7PIU,7VVK | 3 | 7/8 | 2: A stable; B stable; membrane stable | 4: A stable; B stable; membrane stable | stable (membrane) |
| `membrane:FACE_FLAT_TILT` | named-entry:32KS | 35.0 | 8/8 | 26.25: A stable; B stable; membrane sensitive; area 4 entries; 4 entries moved | 43.75: A stable; B stable; membrane sensitive; area 4 entries; 4 entries moved | sensitive (membrane) |
| `membrane:JOIN_GAP` | named-entry:32KS,6UW4 | 2 | 7/8 | skipped | 3: A stable; B stable; membrane stable | stable (membrane) |
| `na_runs:_JOINT_NEAR` | named-entry:9C75 | 2.5 | 6/100 | 1.875: A stable; B stable | 3.125: A stable; B stable | stable (B) |
| `na_runs:_JOINT_REACH` | named-entry:8VWQ | 5.0 | 2/100 | 3.75: A stable; B stable | 6.25: A stable; B stable | stable (B) |
| `na_runs:_BESIDE_TOLL` | named-entry:7Y7C | 50.0 | 1/100 | 37.5: A stable; B stable | 62.5: A stable; B stable | stable (B) |
| `na_runs:_DETOUR` | measured | 3.0 | 12/100 | 2.25: A stable; B stable | 3.75: A stable; B stable | stable (B) |
| `na_runs:_DETOUR_STEPS` | measured | 6.0 | 12/100 | 4.5: A stable; B sensitive; counts_hidden 15->14; 1 entry moved | 7.5: A stable; B stable | sensitive (B) |
| `na_runs:_FIRST_NODES` | measured | 5000 | 8/100 | 3750: A stable; B stable | 6250: A stable; B stable | stable (B) |
| `na_runs:_BESIDE_DETOUR` | measured | 2.0 | 12/100 | 1.5: A stable; B sensitive; counts_hidden 15->14; 1 entry moved | 2.5: A sensitive; area 1 entry; 1 entry moved; B sensitive; area 2 entries; 2 entries moved | sensitive (B) |
| `nucleic:HBOND_STRETCH` | named-entry:1HYS,1AOI | 4.0 | 16/100 | 3.0: A sensitive; area 1 entry; 1 entry moved; B sensitive; A17 7->8; A6 4->3; counts_hidden 15->17; helix_name_moved 41->36; na_gap_unclear 3->5; na_tail_straight 4->5; page_text_small 2->3; area 7 entries; 7 entries moved | 5.0: A sensitive; area 2 entries; 2 entries moved; B sensitive; pass -1; A17 7->9; A6 4->5; counts_hidden 15->22; na_gap_unclear 3->4; page_text_small 2->3; area 8 entries; 8 entries moved | sensitive (B) |
| `panels:CLOSE` | named-entry:2LZM | 2.0 | 49/49 | 1.5: A stable; B stable; domains sensitive; A16 1->2; counts_hidden 19->20; helix_name_moved 45->43; area 1 entry; 3 entries moved | 2.5: A stable; B stable; domains sensitive; counts_hidden 19->20; helix_name_moved 45->43; loop_over_text 11->10; route_second_search_refused 14->15; area 1 entry; 3 entries moved | sensitive (domains) |
| `render:_SPOT_FAR` | named-entry:7QUP | 3 | 4/100 | 2: A stable; B stable | 4: A stable; B stable | stable (B) |
| `render:_SPOT_MORE_RINGS` | named-entry:7LX9,9YLM,3MHT,1GOX | 4 | 55/100 | 3: A stable; B stable | 5: A stable; B stable | stable (B) |
| `route:_SLAB_RUN` | named-entry:7SMR | 12.0 | 100/100 | 9.0: A stable; B stable | 15.0: A stable; B stable | stable (B) |
| `route:_MAX_NODES_WIDE` | measured | 1600000 | 100/100 | 1200000: A stable; B stable | 2000000: A stable; B stable | stable (B) |
| `route:_MAX_NODES_PAGE` | named-entry:7QUP | 8000000 | 100/100 | 6000000: A stable; B stable | 10000000: A stable; B stable | stable (B) |
| `route:_MAX_NODES_WINDOWS` | named-entry:7NZ2,9FZL | 8000000 | 100/100 | 6000000: A stable; B stable | 10000000: A stable; B stable | stable (B) |
| `route:_NA_TEXT_TOLL` | named-entry:6SO9,9CJI,7XSP | 25.0 | 2/100 | 18.75: A stable; B stable | 31.25: A stable; B stable | stable (B) |
| `route:_REROUTE_TURNS_INK` | named-entry:7XSP | 12 | 2/100 | 9: A stable; B stable | 15: A stable; B stable | stable (B) |
| `symmetry:MIN_HELIX_COPIES` | named-entry:1URN | 4 | 54/100 | 3: A stable; B stable | 5: A stable; B stable | stable (B) |

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
