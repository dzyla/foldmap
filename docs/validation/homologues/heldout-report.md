# Homologue held-out set: measured once (plan 13a, Task 13a-11)

The set is `heldout.json`, drawn by `python scripts/homologues.py sample --seed 20261021 --out heldout.json` from
CATH v4_3_0 (CC BY 4.0, https://www.cathdb.info) before any held-out pair was drawn by foldmap (commit c099787):

- 20 **remote** pairs: one CATH topology, two S35 clusters, of different superfamilies or of one superfamily at
  structure-alignment identity under 30 %;
- 10 **close** pairs: one superfamily and S35 cluster, two S95 clusters, identity 30-90 %;
- 10 **negative** pairs: one architecture, two topologies.

Only whole-chain CATH domains (domain number 00) of 60-600 residues; no entry of the tuning set; each entry once.
16 drawn pairs were rejected (identity outside the stratum, or an entry already used); `heldout.json` lists them
with their reasons.

The run: `python scripts/homologues.py run --set docs/validation/homologues/heldout.json --out DIR --jobs 20` on
zyla-office at c099787, once, with the limits frozen at 76fd936 (SAME_FOLD_TM 0.5, SAME_FOLD_FRACTION 0.5, H2 1.0 /
2.5, H3 0.95, `_K_LIKE` 128). All 40 pairs ran to the end (exit 0). Nothing was changed after the run; every
failure below stands as measured. The summary is also in `heldout-summary.json`. Images stay local.

**Positives:** 16 of 30 pass H1-H5 (close 8 of 10, remote 8 of 20). **Negatives:** 9 of 10 refused with the note;
1 not judged (1RPU-1ZWT: the reference's frame, a straight duplex, cannot be carried, so B was not aligned).
**Drawn independently: H2 median 1.89, H3 0.71; drawn alike: H2 median 0.11, H3 1.00** (medians over the 22
positives with H2 judged); 22 of 22 sit closer to A drawn alike. **The aligner against US-align was not measured
in this run** (US-align was not installed); it was measured afterwards, with nothing here changed: `usalign-report.md`
(plan 23 Task 3; US-align puts all seven H1 refusals under TM 0.5 too).

TM: foldmap's TM-score normalised by A. Elements: mates / pieces / extras (B only) / only in A. H2: median / 90th percentile of the mates' centre deviations (page units) after the best rigid 2D motion; H3: share of mate pairs whose page relation agrees (pairs judged). 'alike': B drawn like A; 'indep.': A and B drawn independently (the baseline). Lost: checks B drawn alone passes and B drawn like A does not (H5).

| pair | kind | stratum | route | TM | elements | H1 | H2 | H3 | H4 | H5 | H2 alike | H2 indep. | H3 alike | H3 indep. | lost |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1CO6_A-1JDL_A | positive | close | sequence | 0.822 | 5/0/1/0 | pass | pass | pass | pass | pass | 0.12 / 0.82 | 1.87 / 3.13 | 1.000 (10) | 0.700 (10) | – |
| 1G9R_A-1JYL_A | positive | remote | structure | 0.502 | 14/2/5/11 | pass | pass | pass | pass | pass | 0.31 / 0.86 | 3.24 / 5.19 | 1.000 (91) | 0.110 (91) | – |
| 1JBI_A-3U7Z_A | positive | remote | structure | 0.246 | 2/0/9/7 | **fail** | – | – | – | pass | – | – | – | – | – |
| 1JK4_A-2HNW_A | positive | close | sequence | 0.894 | 8/1/0/2 | pass | **fail** | **fail** | pass | pass | 1.12 / 1.74 | 2.38 / 4.83 | 0.964 (28) | 0.250 (28) | – |
| 1LSH_B-5MLC_3 | negative | negative | structure | 0.190 | 3/0/1/10 | pass | – | – | – | pass | – | 0.88 / 1.17 | – | 0.000 (3) | – |
| 1NKP_A-1AM9_A | positive | remote | structure | 0.742 | 2/0/0/1 | pass | – | – | pass | pass | – | – | – | – | – |
| 1PZW_A-1QHL_A | negative | negative | structure | 0.284 | 2/0/14/2 | pass | – | – | – | pass | – | – | – | – | – |
| 1RPU_A-1ZWT_A | negative | negative | – | 0.000 | 0/0/0/0 | – | – | – | – | pass | – | – | – | – | – |
| 1T5B_A-2Z9C_A | positive | close | sequence | 0.977 | 15/0/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.00 | 0.01 / 0.03 | 1.000 (105) | 1.000 (105) | – |
| 1TYY_A-3IQ0_A | positive | remote | structure | 0.887 | 27/0/1/1 | pass | pass | pass | pass | pass | 0.09 / 1.04 | 3.33 / 7.48 | 0.997 (351) | 0.724 (351) | – |
| 1XU6_A-1PS2_A | positive | remote | structure | 0.223 | 1/0/2/3 | **fail** | – | – | – | pass | – | – | – | – | – |
| 2C3E_A-1HZF_A | negative | negative | structure | 0.283 | 8/1/6/4 | pass | – | – | – | pass | – | 6.59 / 9.78 | – | 0.214 (28) | – |
| 2D46_A-3F42_B | positive | remote | structure | 0.352 | 2/0/6/3 | **fail** | – | – | – | pass | – | – | – | – | – |
| 2EWR_A-2FCL_A | positive | close | sequence | 0.996 | 14/0/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.06 | 0.01 / 0.03 | 1.000 (91) | 1.000 (91) | – |
| 2EX5_A-4NTQ_A | negative | negative | structure | 0.148 | 3/0/4/10 | pass | – | – | – | pass | – | 0.30 / 0.49 | – | 1.000 (3) | – |
| 2F20_A-1ZN6_A | positive | remote | structure | 0.773 | 18/0/2/8 | pass | **fail** | **fail** | pass | pass | 2.44 / 6.18 | 4.41 / 11.31 | 0.673 (153) | 0.379 (153) | – |
| 2H5F_A-4Z42_B | negative | negative | structure | 0.291 | 3/0/6/3 | pass | – | – | – | pass | – | 1.10 / 1.48 | – | 0.667 (3) | – |
| 2IAY_A-4XZD_A | negative | negative | structure | 0.316 | 3/0/12/7 | pass | – | – | – | pass | – | 0.00 / 0.00 | – | 0.000 (3) | – |
| 2IJE_S-3T6A_A | positive | remote | structure | 0.731 | 13/2/5/1 | pass | pass | pass | pass | pass | 0.83 / 2.13 | 6.21 / 10.17 | 0.987 (78) | 0.167 (78) | – |
| 2JCA_A-3LQ9_A | positive | remote | structure | 0.369 | 6/0/5/5 | **fail** | – | – | – | pass | – | 2.57 / 2.72 | – | 0.067 (15) | – |
| 2LSR_A-4YL6_A | positive | remote | structure | 0.806 | 5/1/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.00 | 2.82 / 4.70 | 1.000 (10) | 0.700 (10) | – |
| 2MJ7_A-2EHB_D | positive | remote | structure | 0.484 | 7/0/4/2 | **fail** | – | – | – | pass | – | 1.07 / 2.11 | – | 0.000 (21) | – |
| 2OKF_A-2INB_A | positive | close | sequence | 0.976 | 13/0/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.00 | 0.01 / 0.02 | 1.000 (78) | 1.000 (78) | – |
| 2Q82_A-2JTR_A | positive | remote | structure | 0.392 | 6/1/3/5 | **fail** | – | – | – | pass | – | 2.63 / 4.78 | – | 0.600 (15) | – |
| 3BUU_A-3VDI_A | negative | negative | structure | 0.267 | 6/0/21/13 | pass | – | – | – | pass | – | 1.61 / 3.07 | – | 0.000 (15) | – |
| 3IE4_A-3AQX_A | positive | close | sequence | 0.957 | 7/0/1/0 | pass | pass | pass | pass | pass | 0.47 / 0.62 | 0.51 / 2.14 | 1.000 (21) | 0.143 (21) | – |
| 3N4H_A-4JCU_A | positive | remote | structure | 0.738 | 8/0/2/4 | pass | pass | pass | pass | pass | 0.09 / 0.16 | 0.37 / 1.43 | 1.000 (28) | 1.000 (28) | – |
| 3O7I_A-2O8I_A | positive | remote | structure | 0.722 | 7/0/2/1 | pass | **fail** | pass | pass | pass | 1.16 / 1.51 | 3.13 / 5.04 | 0.952 (21) | 0.191 (21) | – |
| 3OQL_A-4WWJ_A | positive | close | sequence | 0.988 | 13/0/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.00 | 2.45 / 5.18 | 1.000 (78) | 0.821 (78) | – |
| 3VDJ_A-2JXW_A | negative | negative | structure | 0.188 | 1/0/5/5 | pass | – | – | – | pass | – | – | – | – | – |
| 3WZ0_A-2AV5_A | positive | close | sequence | 0.970 | 7/0/1/1 | pass | **fail** | **fail** | pass | pass | 0.93 / 3.26 | 1.90 / 3.46 | 0.714 (21) | 0.476 (21) | – |
| 4E2U_A-2LWY_A | positive | remote | structure | 0.695 | 14/0/0/3 | pass | **fail** | **fail** | pass | pass | 0.07 / 7.58 | 5.92 / 11.19 | 0.758 (91) | 0.330 (91) | – |
| 4IPV_A-4N2P_A | negative | negative | structure | 0.260 | 5/1/3/6 | pass | – | – | – | pass | – | 2.38 / 3.07 | – | 0.200 (10) | – |
| 4U12_A-4NPX_A | positive | remote | structure | 0.763 | 4/0/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.00 | 0.21 / 0.47 | 1.000 (6) | 1.000 (6) | – |
| 5E4E_A-3TGX_B | positive | remote | structure | 0.609 | 6/0/0/2 | pass | pass | **fail** | pass | pass | 0.29 / 1.04 | 2.35 / 2.92 | 0.933 (15) | 0.533 (15) | – |
| 5GL3_A-1WPN_A | positive | remote | structure | 0.431 | 11/0/3/10 | **fail** | – | – | – | pass | – | 6.40 / 9.67 | – | 0.473 (55) | – |
| 5GY2_C-2D4X_A | positive | remote | structure | 0.785 | 6/0/1/2 | pass | pass | pass | pass | pass | 0.18 / 0.46 | 1.78 / 3.62 | 1.000 (15) | 0.800 (15) | – |
| 5KDI_A-4NT2_A | positive | remote | structure | 0.793 | 11/0/1/3 | pass | pass | pass | pass | pass | 0.80 / 1.28 | 1.51 / 4.92 | 1.000 (55) | 0.891 (55) | – |
| 5LNK_A-6G2J_A | positive | close | sequence | 0.698 | 4/0/0/0 | pass | pass | pass | pass | pass | 0.00 / 0.00 | 0.51 / 0.54 | 1.000 (6) | 1.000 (6) | – |
| 5U5S_A-2WP1_A | positive | close | sequence | 0.937 | 7/0/0/0 | pass | pass | pass | pass | pass | 0.06 / 0.21 | 0.36 / 0.68 | 1.000 (21) | 1.000 (21) | – |

## Failures by check

- **H1**: 
  - 1JBI_A-3U7Z_A (remote): a positive pair refused with the note: TM-score 0.25 (one fold at 0.5), 2 mates of the 9 elements of the chain with fewer (at 0.5)
  - 1XU6_A-1PS2_A (remote): a positive pair refused with the note: TM-score 0.22 (one fold at 0.5), 1 mates of the 3 elements of the chain with fewer (at 0.5)
  - 2D46_A-3F42_B (remote): a positive pair refused with the note: TM-score 0.35 (one fold at 0.5), 2 mates of the 5 elements of the chain with fewer (at 0.5)
  - 2JCA_A-3LQ9_A (remote): a positive pair refused with the note: TM-score 0.37 (one fold at 0.5), 6 mates of the 11 elements of the chain with fewer (at 0.5)
  - 2MJ7_A-2EHB_D (remote): a positive pair refused with the note: TM-score 0.48 (one fold at 0.5), 7 mates of the 9 elements of the chain with fewer (at 0.5)
  - 2Q82_A-2JTR_A (remote): a positive pair refused with the note: TM-score 0.39 (one fold at 0.5), 6 mates of the 10 elements of the chain with fewer (at 0.5)
  - 5GL3_A-1WPN_A (remote): a positive pair refused with the note: TM-score 0.43 (one fold at 0.5), 11 mates of the 14 elements of the chain with fewer (at 0.5)
- **H2**: 
  - 1JK4_A-2HNW_A (close): median 1.117, 90th percentile 1.741 page units over 8 mates (limits 1.0, 2.5); farthest H~H 3.031; drawn independently: median 2.380, 90th percentile 4.825
  - 2F20_A-1ZN6_A (remote): median 2.445, 90th percentile 6.179 page units over 18 mates (limits 1.0, 2.5); farthest E~E 11.758; drawn independently: median 4.406, 90th percentile 11.306
  - 3O7I_A-2O8I_A (remote): median 1.159, 90th percentile 1.508 page units over 7 mates (limits 1.0, 2.5); farthest α6~α6 1.572; drawn independently: median 3.132, 90th percentile 5.042
  - 3WZ0_A-2AV5_A (close): median 0.929, 90th percentile 3.264 page units over 7 mates (limits 1.0, 2.5); farthest B~B 3.496; drawn independently: median 1.900, 90th percentile 3.464
  - 4E2U_A-2LWY_A (remote): median 0.073, 90th percentile 7.580 page units over 14 mates (limits 1.0, 2.5); farthest I~I 13.905; drawn independently: median 5.918, 90th percentile 11.188
- **H3**: 
  - 1JK4_A-2HNW_A (close): 0.9643 of 28 pairs agree (limit 0.95); sheet order G, H for H, G; drawn independently: 0.2500 of 28 pairs
  - 2F20_A-1ZN6_A (remote): 0.6732 of 153 pairs agree (limit 0.95); 50 turned: A-α1, A-α4, A-E, A-H, A-K, A-α7; sheet order A, H, B, M for M, H, A, B; sheet order C, D, K, J, E, I, F for F, I, E, J, K, D, C; drawn independently: 0.3791 of 153 pairs
  - 3WZ0_A-2AV5_A (close): 0.7143 of 21 pairs agree (limit 0.95); 6 turned: A-B, A-C, A-D, B-C, B-D, C-D; sheet order D, A, C, B for B, C, A, D; drawn independently: 0.4762 of 21 pairs
  - 4E2U_A-2LWY_A (remote): 0.7582 of 91 pairs agree (limit 0.95); 22 turned: A-F, A-H, A-I, A-M, B-I, C-I; sheet order A, M, L, H, I for H, I, M, A, L; drawn independently: 0.3297 of 91 pairs
  - 5E4E_A-3TGX_B (remote): 0.9333 of 15 pairs agree (limit 0.95); 1 turned: α2-α5; drawn independently: 0.5333 of 15 pairs
- **H4**: none
- **H5**: none

## Summary

| stratum | pairs | results | all pass | H1 | H2 | H3 | H4 | H5 |
|---|---|---|---|---|---|---|---|---|
| close | 10 | 10 | 8 | 10 / 0 / 0 | 8 / 2 / 0 | 8 / 2 / 0 | 10 / 0 / 0 | 10 / 0 / 0 |
| remote | 20 | 20 | 8 | 13 / 7 / 0 | 9 / 3 / 8 | 9 / 3 / 8 | 13 / 0 / 7 | 20 / 0 / 0 |
| negative | 10 | 10 | 0 | 9 / 0 / 1 | 0 / 0 / 10 | 0 / 0 / 10 | 0 / 0 / 10 | 10 / 0 / 0 |

Cells: pass / fail / not judged.

- Positives passing H1-H5: 16 of 30; negatives refused (H1): 9 of 10.
- Over the 22 positives with H2 judged (medians of the pairs' values): H2 median drawn independently 1.89, drawn alike 0.11; H2 90th percentile 3.54 and 0.72; H3 share 0.71 and 1.00. Closer drawn alike (H2 median): 22 of 22.
