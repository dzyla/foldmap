# The topology carried by --like

Plan 28 Task 6 (the review's 3.5, ruling Q44). Written by: foldmap 0.2.0, git e1392f79daf7 (asserted by runner, dirty False), code fingerprint 7879eea05bd9301c, scripts/topology_transfer.py sha256 731bc36e96e0; Python 3.13.13 on Linux-6.8.0-138-generic-x86_64-with-glibc2.35; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Whether a sheet drawn like a reference keeps the reference's strand order and direction. Each pair is drawn like A twice: **carried**, as foldmap draws it (each of B's sheets reads left to right as A reads its mates, Like.order, where like._order says), and **withheld**, the same with like._order answering None, so each sheet's order is its own (layout._sheet_order's default). Both are judged by H1-H5 as the homologue campaign judges them (scripts/homologues.py; limits criteria.LIMITS, frozen at 76fd936). The record is made by scripts/topology_transfer.py; report.json holds every pair and sheet.

## The sets

- **hom_tuning** (tuning): plan 13a's homologue tuning pairs (docs/validation/homologues/tuning.json: 11 positive, 4 negative); the H1-H5 limits (criteria.LIMITS like_*) were frozen at 76fd936 while looking at this set's results (docs/validation/homologues/tuning-report.md).
- **hom_heldout** (unseen at c099787, historical): plan 13a's held-out homologue pairs (docs/validation/homologues/heldout.json: 30 positive, 10 negative); one run at c099787 with the limits frozen at 76fd936 on the tuning set; nothing changed after the run; its count describes c099787, not the release (reruns since give 17 of 30 positives against the frozen 16; docs/validation/homologues/review.md).
- Neither set is unseen for this record: hom_tuning is the set the limits were chosen on, and the fix measured here (Like.order, 1a3f9d6, plan 13a-12) was aimed at 3WZ0-2AV5, a hom_heldout pair, after its one run. The carried column is a regression measure on both sets, not evidence on unseen structures.

- Withheld is today's code with like._order answering None, judged by today's criteria, not the code before 1a3f9d6: that commit also removed H3's pairing excuse and made H1's refusal a flag; withheld reproduces the c099787 counts (16 of 30, H2 17, H3 17) all the same.

## The mechanism (3WZ0-2AV5)

TM-score 0.97 by sequence alignment. A's four-strand sheet is seen nearly edge-on through the carried frame, and B's sheet has a fifth strand at one end that A lacks (C′).

| way | Like.order | mates read / point as A's | H1 H2 H3 H4 H5 | H2 median / p90 | H3 share | sheet in another order |
|---|---|---|---|---|---|---|
| withheld | None (withheld) | no / yes | pass **fail** **fail** pass pass | 0.93 / 3.26 | 0.714 | D, A, C, B for B, C, A, D |
| carried | True | yes / yes | pass pass pass pass pass | 0.23 / 0.72 | 1.000 | none |

- **Where the order is lost** (withheld): layout._sheet_order draws a sheet's track left to right across its drawn-up direction, and reverses it when its last track strand's projected centroid lies left of its first. For a sheet seen nearly edge-on that step is a fraction of an ångström, and in B the last track strand is C′, the strand A lacks, so B's own end strands decide an order A's mates never had: B C A D comes out D A C B, each strand pointing as its mate, the mirror image of A's sheet. Nothing carried from A entered the order: like.like_from carried places, angles (the sheet's way up, from its first mated strand) and names, not the reading.
- **The fix** (1a3f9d6): like._order reads A's drawn order of the mates along B's track (Like.order: True left to right, False right to left), and layout._sheet_order takes it in place of the end strands' step, the way up still from A's angle. It is set only where B's strands are mapped and paired as in A: two or more track strands mated in one of A's sheets, pointing as B's pairing has them, in one order there; else None and the sheet's own order (the reasons below).

## Both sets, carried and withheld

| set | way | positives passing H1-H5 | H1 | H2 | H3 | H4 | H5 | negatives refused | sheets: order set | ... read and point as A's | pairs with a sheet in another order |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hom_tuning | carried | 8 of 11 | 11 | 8 | 8 | 11 | 11 | 4 of 4 | 25 of 31 | 25 of 25 | 3 |
| hom_tuning | withheld | 8 of 11 | 11 | 8 | 8 | 11 | 11 | 4 of 4 | 0 of 31 | 0 of 0 | 3 |
| hom_heldout | carried | 17 of 30 | 23 | 18 | 19 | 23 | 30 | 9 of 10 | 21 of 25 | 21 of 21 | 2 |
| hom_heldout | withheld | 16 of 30 | 23 | 17 | 17 | 23 | 30 | 9 of 10 | 0 of 25 | 0 of 0 | 4 |

- **hom_tuning**: 2 pair(s) whose verdicts or sheet orders differ carried and withheld: 2ZB5_A-1E8T_A (sheet orders only); 4CH9_A-2VWD_A (sheet orders only). Pairs losing an H check carried: none.
- **hom_heldout**: 2 pair(s) whose verdicts or sheet orders differ carried and withheld: 1JK4_A-2HNW_A (H3); 3WZ0_A-2AV5_A (H2, H3). Pairs losing an H check carried: none.

## Sheets still in another order (carried)

Each pair with a sheet of B that H3 finds in another order than its mates in A (criteria.sheet_orders), and B's sheets with two mated strands or more there: Like.order's entry or why it is None, and whether the mated track strands read (r) and point (p) as their mates, and the accessory strands read with them (a).

| set | pair | TM | H2 / H3 | sheets in another order (B's for A's) | B's sheets: Like.order, r p a |
|---|---|---|---|---|---|
| hom_tuning | 2VWD_A-2ZB5_A | 0.72 | **fail** / **fail** | L, M, O, N, K for K, M, L, N, O | True, y y y; False, y y y; None: mates in two of A's sheets, n y n; True, y y y; False, y y y |
| hom_tuning | 2ZB5_A-1E8T_A | 0.70 | **fail** / **fail** | K, M, N, L, O for L, M, O, K, N | False, y y y; None: mates in two of A's sheets, y y y; False, y y y; False, y y y; None: mates not in one order across A's sheet, n y n; True, y y y |
| hom_tuning | 4CH9_A-2VWD_A | 0.69 | **fail** / **fail** | M, P, L, Q for L, M, P, Q | True, y y y; False, y y y; False, y y y; True, y y n; True, y y y; True, y y y |
| hom_heldout | 2F20_A-1ZN6_A | 0.78 | **fail** / **fail** | A, H, B, M for M, H, A, B; C, D, K, J, E, I, F for F, I, E, J, K, D, C | None: mates in two of A's sheets, n y n |
| hom_heldout | 4E2U_A-2LWY_A | 0.70 | **fail** / **fail** | A, M, L, H, I for H, I, M, A, L | None: mates in two of A's sheets, n n n; False, y y y; False, y y y; False, y y y |

## Every pair

| set | pair | kind | TM | carried: H1 H2 H3 H4 H5 | H2 | H3 | withheld: H1 H2 H3 H4 H5 | H2 | H3 |
|---|---|---|---|---|---|---|---|---|---|
| hom_tuning | 1EMA_A-1G7K_A | positive | 0.89 | pass pass pass pass pass | 0.28 / 0.36 | 1.000 | pass pass pass pass pass | 0.28 / 0.36 | 1.000 |
| hom_tuning | 1HDC_A-1CYD_A | positive | 0.86 | pass pass pass pass pass | 0.18 / 0.33 | 1.000 | pass pass pass pass pass | 0.18 / 0.33 | 1.000 |
| hom_tuning | 1MBN_A-1ECA_A | positive | 0.80 | pass pass pass pass pass | 0.09 / 0.17 | 1.000 | pass pass pass pass pass | 0.09 / 0.17 | 1.000 |
| hom_tuning | 1REI_A-1TNF_A | negative | 0.42 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_tuning | 1SHG_A-1CKA_A | positive | 0.84 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_tuning | 1TIM_A-1HDC_A | negative | 0.40 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_tuning | 1TIM_A-1IGS_A | positive | 0.62 | pass pass pass pass pass | 0.00 / 0.03 | 1.000 | pass pass pass pass pass | 0.00 / 0.03 | 1.000 |
| hom_tuning | 1UBQ_A-1WM3_A | positive | 0.80 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_tuning | 1UBQ_A-2ACY_A | negative | 0.35 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_tuning | 256B_A-1MBN_A | negative | 0.44 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_tuning | 2ACY_A-1URN_A | positive | 0.52 | pass pass pass pass pass | 0.08 / 0.19 | 1.000 | pass pass pass pass pass | 0.08 / 0.19 | 1.000 |
| hom_tuning | 2VWD_A-2X9M_A | positive | 0.98 | pass pass pass pass pass | 0.00 / 0.06 | 1.000 | pass pass pass pass pass | 0.00 / 0.06 | 1.000 |
| hom_tuning | 2VWD_A-2ZB5_A | positive | 0.72 | pass **fail** **fail** pass pass | 0.27 / 4.37 | 0.889 | pass **fail** **fail** pass pass | 0.27 / 4.37 | 0.889 |
| hom_tuning | 2ZB5_A-1E8T_A | positive | 0.70 | pass **fail** **fail** pass pass | 2.69 / 6.35 | 0.844 | pass **fail** **fail** pass pass | 3.10 / 5.80 | 0.876 |
| hom_tuning | 4CH9_A-2VWD_A | positive | 0.69 | pass **fail** **fail** pass pass | 2.10 / 5.22 | 0.948 | pass **fail** **fail** pass pass | 1.98 / 5.22 | 0.918 |
| hom_heldout | 1CO6_A-1JDL_A | close | 0.82 | pass pass pass pass pass | 0.12 / 0.82 | 1.000 | pass pass pass pass pass | 0.12 / 0.82 | 1.000 |
| hom_heldout | 1G9R_A-1JYL_A | remote | 0.51 | pass pass pass pass pass | 0.31 / 0.85 | 1.000 | pass pass pass pass pass | 0.31 / 0.85 | 1.000 |
| hom_heldout | 1JBI_A-3U7Z_A | remote | 0.27 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 1JK4_A-2HNW_A | close | 0.89 | pass **fail** pass pass pass | 1.11 / 2.00 | 1.000 | pass **fail** **fail** pass pass | 1.12 / 1.74 | 0.964 |
| hom_heldout | 1LSH_B-5MLC_3 | negative | 0.19 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 1NKP_A-1AM9_A | remote | 0.75 | pass – – pass pass | – | – | pass – – pass pass | – | – |
| hom_heldout | 1PZW_A-1QHL_A | negative | 0.32 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 1RPU_A-1ZWT_A | negative | – | – – – – pass | – | – | – – – – pass | – | – |
| hom_heldout | 1T5B_A-2Z9C_A | close | 0.98 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_heldout | 1TYY_A-3IQ0_A | remote | 0.89 | pass pass pass pass pass | 0.09 / 1.04 | 0.997 | pass pass pass pass pass | 0.09 / 1.04 | 0.997 |
| hom_heldout | 1XU6_A-1PS2_A | remote | 0.26 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 2C3E_A-1HZF_A | negative | 0.29 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 2D46_A-3F42_B | remote | 0.38 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 2EWR_A-2FCL_A | close | 1.00 | pass pass pass pass pass | 0.00 / 0.06 | 1.000 | pass pass pass pass pass | 0.00 / 0.06 | 1.000 |
| hom_heldout | 2EX5_A-4NTQ_A | negative | 0.16 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 2F20_A-1ZN6_A | remote | 0.78 | pass **fail** **fail** pass pass | 2.44 / 6.18 | 0.673 | pass **fail** **fail** pass pass | 2.44 / 6.18 | 0.673 |
| hom_heldout | 2H5F_A-4Z42_B | negative | 0.30 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 2IAY_A-4XZD_A | negative | 0.33 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 2IJE_S-3T6A_A | remote | 0.74 | pass pass pass pass pass | 0.83 / 2.13 | 0.987 | pass pass pass pass pass | 0.83 / 2.13 | 0.987 |
| hom_heldout | 2JCA_A-3LQ9_A | remote | 0.40 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 2LSR_A-4YL6_A | remote | 0.81 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_heldout | 2MJ7_A-2EHB_D | remote | 0.49 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 2OKF_A-2INB_A | close | 0.98 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_heldout | 2Q82_A-2JTR_A | remote | 0.40 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 3BUU_A-3VDI_A | negative | 0.33 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 3IE4_A-3AQX_A | close | 0.96 | pass pass pass pass pass | 0.47 / 0.62 | 1.000 | pass pass pass pass pass | 0.47 / 0.62 | 1.000 |
| hom_heldout | 3N4H_A-4JCU_A | remote | 0.74 | pass pass pass pass pass | 0.04 / 0.12 | 1.000 | pass pass pass pass pass | 0.04 / 0.12 | 1.000 |
| hom_heldout | 3O7I_A-2O8I_A | remote | 0.72 | pass **fail** pass pass pass | 1.16 / 1.51 | 0.952 | pass **fail** pass pass pass | 1.16 / 1.51 | 0.952 |
| hom_heldout | 3OQL_A-4WWJ_A | close | 0.99 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_heldout | 3VDJ_A-2JXW_A | negative | 0.26 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 3WZ0_A-2AV5_A | close | 0.97 | pass pass pass pass pass | 0.23 / 0.72 | 1.000 | pass **fail** **fail** pass pass | 0.93 / 3.26 | 0.714 |
| hom_heldout | 4E2U_A-2LWY_A | remote | 0.70 | pass **fail** **fail** pass pass | 0.07 / 7.58 | 0.758 | pass **fail** **fail** pass pass | 0.07 / 7.58 | 0.758 |
| hom_heldout | 4IPV_A-4N2P_A | negative | 0.30 | pass – – – pass | – | – | pass – – – pass | – | – |
| hom_heldout | 4U12_A-4NPX_A | remote | 0.77 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_heldout | 5E4E_A-3TGX_B | remote | 0.62 | pass pass **fail** pass pass | 0.29 / 1.04 | 0.933 | pass pass **fail** pass pass | 0.29 / 1.04 | 0.933 |
| hom_heldout | 5GL3_A-1WPN_A | remote | 0.45 | **fail** – – – pass | – | – | **fail** – – – pass | – | – |
| hom_heldout | 5GY2_C-2D4X_A | remote | 0.80 | pass pass pass pass pass | 0.30 / 0.44 | 1.000 | pass pass pass pass pass | 0.30 / 0.44 | 1.000 |
| hom_heldout | 5KDI_A-4NT2_A | remote | 0.80 | pass pass pass pass pass | 0.80 / 1.28 | 1.000 | pass pass pass pass pass | 0.80 / 1.28 | 1.000 |
| hom_heldout | 5LNK_A-6G2J_A | close | 0.70 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 | pass pass pass pass pass | 0.00 / 0.00 | 1.000 |
| hom_heldout | 5U5S_A-2WP1_A | close | 0.94 | pass pass pass pass pass | 0.06 / 0.21 | 1.000 | pass pass pass pass pass | 0.06 / 0.21 | 1.000 |
