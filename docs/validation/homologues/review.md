# Homologue held-out set: review (plan 13a, Task 13a-11)

One run at c099787 with the limits frozen at 76fd936; the table is `heldout-report.md`, the numbers
`heldout-summary.json`. Nothing was changed after the run. Every failure below is reported as measured, with what
the run's own output and pages show; none is fixed here.

## The numbers

| | pairs | H1 | H2 | H3 | H4 | H5 | all five |
|---|---|---|---|---|---|---|---|
| close (positive) | 10 | 10 pass | 8 pass, 2 fail | 8 pass, 2 fail | 10 pass | 10 pass | 8 |
| remote (positive) | 20 | 13 pass, 7 fail | 9 pass, 3 fail, 8 not judged | 9 pass, 3 fail, 8 not judged | 13 pass, 7 not judged | 20 pass | 8 |
| negative | 10 | 9 refused, 1 not judged | not judged | not judged | not judged | 10 pass | – |

- Positives passing H1-H5: **16 of 30** (tuning set: 8 of 11). Per check over the 30 positives: H1 23, H2 17,
  H3 17, H4 23, H5 30 pass. A pair refused (H1 fail) leaves H2-H4 not judged; 1NKP-1AM9 (two helices, 2 mates)
  passes H1 and has too few mates for H2/H3.
- Negatives refused: **9 of 10**; 1RPU-1ZWT not judged (below). No negative was forced alike.
- **Drawn independently: H2 median 1.89 page units, H3 0.71; drawn alike: H2 median 0.11, H3 1.00** (medians of
  the pairs' values over the 22 positives with H2 judged; H2 90th percentile 3.54 and 0.72). Every one of the 22 is
  closer to A drawn alike than drawn independently, including the five that fail H2 or H3.
- H5: no pair lost a check by being drawn alike (40 of 40).
- Routes: all 10 close pairs took the sequence route (TM 0.70-1.00); all remote and negative pairs the structure
  route.
- **The aligner against US-align: not measured.** US-align is not installed on zyla-office (ledger: the install
  waits for the author). The H1 failures below are where that comparison matters most.

## Failures by check

### H1: 7 remote positives refused (TM under 0.5)

| pair | CATH | foldmap TM | mates |
|---|---|---|---|
| 1JBI_A-3U7Z_A | 2.170.130, superfamilies .20 / .30 | 0.25 | 2 of 9 |
| 1XU6_A-1PS2_A | 4.10.110, superfamilies .20 / .10 | 0.22 | 1 of 3 |
| 2D46_A-3F42_B | 3.30.1310, superfamilies .30 / .10 | 0.35 | 2 of 5 |
| 2JCA_A-3LQ9_A | 3.90.470, superfamilies .20 / .40 | 0.37 | 6 of 11 |
| 2MJ7_A-2EHB_D | 3.30.310, superfamilies .10 / .80 | 0.48 | 7 of 9 |
| 2Q82_A-2JTR_A | 3.40.250, superfamilies .20 / .10 | 0.39 | 6 of 10 |
| 5GL3_A-1WPN_A | 3.90.1640, superfamilies .20 / .10 | 0.43 | 11 of 14 |

All seven are cross-superfamily pairs of one topology, so the "positive" label is CATH's topology, not a measured
TM-score. Two readings are possible and this run cannot tell them apart without US-align: the pairs really are
under TM 0.5 (CATH topologies hold pairs that far apart, and then H1's refusal is right and the label is too
generous), or foldmap's aligner under-scores them. 2MJ7-2EHB (0.48, 7 of 9 matched) and 5GL3-1WPN (0.43, 11 of 14)
are the ones nearest the limit. The US-align comparison of Task 13a-3 is what decides it.

### H2 and H3: 5 positives fail both or one, H2 1 more

| pair | stratum | TM | H2 alike (indep.) | H3 alike (indep.) | what the run shows |
|---|---|---|---|---|---|
| 3WZ0_A-2AV5_A | close | 0.97 | 0.93 / 3.26 (1.90 / 3.46) | 0.714 (0.476) | the sheet B C A D is drawn D A C B in B: the whole sheet mirrored |
| 4E2U_A-2LWY_A | remote | 0.70 | 0.07 / 7.58 (5.92 / 11.19) | 0.758 (0.330) | B's main sheet pairs its strands in another order (pairing difference) |
| 2F20_A-1ZN6_A | remote | 0.77 | 2.44 / 6.18 (4.41 / 11.31) | 0.673 (0.379) | two sheets in another strand order (pairing difference) |
| 1JK4_A-2HNW_A | close | 0.89 | 1.12 / 1.74 (2.38 / 4.83) | 0.964 (0.250) | H2 median over the limit by 0.12; H3 fails on one swapped strand pair G, H |
| 3O7I_A-2O8I_A | remote | 0.72 | 1.16 / 1.51 (3.13 / 5.04) | 0.952 (0.191) | H2 median over the limit by 0.16; H3 passes |
| 5E4E_A-3TGX_B | remote | 0.61 | 0.29 / 1.04 (2.35 / 2.92) | 0.933 (0.533) | H3 only: one turned helix pair α2-α5 (14 of 15 agree) |

- **3WZ0-2AV5 is a mechanism failure**, not a difference in the structures. At TM 0.97 by sequence, B's
  four-stranded sheet is drawn as the mirror of A's (D A C B for B C A D), its connections crossing the sheet, and A's
  α2 has no mate (B shows an extra strand C′ instead). A sheet carried in reverse order is what H3 is meant to catch
  ("a sheet that is only drawn flipped must still fail H3", 13a-10). Not fixed: held-out rule.
- **4E2U-2LWY and 2F20-1ZN6** are pairing differences of the kind the tuning set's propellers showed: B's sheets
  pair their strands in another order (H3's detail: sheet order A, M, L, H, I for H, I, M, A, L in 4E2U-2LWY; two
  sheets reordered in 2F20-1ZN6), so each sheet keeps its own order and those strands sit away from their mates;
  4E2U-2LWY's H2 median is 0.07, its 90th percentile 7.58. H3's "pairing" excuse is still unused (13a-10 concern 2),
  so these count as failures; the spec counts them as failures of H3's share in any case.
- **1JK4-2HNW and 3O7I-2O8I** fail by small margins (H2 median 1.12 and 1.16 against 1.0). **5E4E-3TGX** is one
  turned pair of 15, the same borderline kind as 2HHB A like B in 13a-9 (14 of 15 = 0.933 under 0.95).

### H4, H5

None. Every judged positive names its mates alike; no pair lost an S- or A-check by being drawn like A.

### Not judged: the negative 1RPU-1ZWT

1RPU is drawn with its DNA as a straight duplex, a frame `--like` does not carry (ruling 7), so B was drawn on its
own with the note "its frame (a straight duplex) cannot be carried; draw it with --no-dna". No alignment, so H1 has
nothing to judge. The campaign draws each pair with foldmap's defaults; this is the documented behaviour.

## What this says about 13a

- Where the pair is one fold by foldmap's measure (23 positives drawn alike), 16 pass everything and all 22 judged
  pairs are drawn much closer to their reference than independently: the H2 median goes from 1.89 to 0.11.
- The weak points, in order: the remote pairs that H1 refuses (7 of 20; needs US-align to know whether the aligner
  or the labels are at fault), one sheet carried mirrored (3WZ0-2AV5), pairing differences that H3 does not yet
  excuse, and three near-limit pairs.
- The held-out pass rate (16 of 30 positives, 9 of 10 negatives refused) is lower than the tuning set's (8 of 11,
  4 of 4), mostly from the remote stratum's H1.
