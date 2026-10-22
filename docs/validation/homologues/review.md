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
- **The aligner against US-align: not measured in this run.** US-align was not installed on zyla-office then. It
  was measured afterwards, without changing anything above (plan 23 Task 3, the last section below and
  `usalign-report.md`).

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
are the ones nearest the limit. The US-align comparison of Task 13a-3 is what decides it. (Decided afterwards: US-align
puts all seven under 0.5 too; see the last section.)

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
- The weak points, in order: the remote pairs that H1 refuses (7 of 20; US-align, measured afterwards, puts all
  seven under TM 0.5 as well, so the labels and not the aligner are at fault), one sheet carried mirrored
  (3WZ0-2AV5), pairing differences that H3 does not yet excuse, and three near-limit pairs.
- The held-out pass rate (16 of 30 positives, 9 of 10 negatives refused) is lower than the tuning set's (8 of 11,
  4 of 4), mostly from the remote stratum's H1.

## The aligner against US-align (plan 23 Task 3, added afterwards)

Measured at c22ba5f with US-align 20241108 (bioconda usalign 20241201) on zyla-office, on the 15 tuning and 40
held-out pairs. Nothing above changed: H1-H5, the limits and the routes stand as measured at c099787. The tables
are in `usalign-report.md`, every pair's numbers and both alignments in `usalign.json`, and the counts in
`heldout-summary.json` under `usalign`.

**How.** `scripts/homologues.py align --usalign auto` aligns each pair with foldmap's structure route
(`align_chains(..., "structure")`). It writes each chain's Cα, as foldmap reads it under `--assembly asu`, to a
one-chain PDB file, and runs `USalign A.pdb B.pdb -mol prot -mm 0 -ter 2`. That is a single-chain, sequence-independent
alignment, with A as Structure_1.
- **Normalisation.** Both TM-scores are normalised by A, the reference: US-align's "normalized by length of
  Structure_1", with d0 of A's length, the same as foldmap's d0.
- **US-align's pairs.** Every column of its alignment that holds a residue of both chains. On all 55 pairs, that
  count equals the "Aligned length" US-align prints.
- **Agreement.** A US-align pair (i, j) is reproduced when foldmap pairs A's residue i with a residue of B within
  ±1 of j. The 13a-3 limits apply where US-align's TM ≥ 0.5: foldmap's TM is within 0.05 of US-align's, and ≥ 85 % of
  US-align's pairs are reproduced.

**Numbers.**
- **34 pairs with US-align TM ≥ 0.5** (11 tuning, 23 held-out). foldmap's TM is within 0.05 on all 34; the largest
  gap is 0.016, on 1TIM-1IGS. Pair agreement is ≥ 85 % on 32. **32 of 34 meet both limits**, with a median of
  98.9 % of pairs reproduced.
  - The held-out set alone: **23 of 23**. The least agreement is 86.7 %, on 4U12-4NPX.
  - The tuning set: 9 of 11.
- **All 55 pairs** fall on the same side of H1's TM 0.5 in both aligners.
  - foldmap is lower on 43 of 55: by 0.005 at the median, and by at most 0.08.
  - Every gap over 0.03 (five pairs) is a negative that both put far under 0.5. The largest: 3BUU-3VDI 0.267
    against 0.348, 3VDJ-2JXW 0.188 against 0.244, 4IPV-4N2P 0.260 against 0.309.
- **Scoring against search.**
  - foldmap's TM-score of US-align's own pairs equals US-align's TM-score within 0.0015 on every pair. The scoring,
    d0 and superposition agree.
  - So every difference between the two lies in the alignment search. On 28 of 55 pairs, foldmap's search stops
    short: US-align's alignment scores higher under foldmap's own TM-score than foldmap's alignment does (by more
    than 0.005). On the rest foldmap matches US-align, or beats it slightly by the score normalised by A, which
    foldmap optimises; TM-align's search, which US-align keeps, uses the shorter chain's d0.

**The two misses** both fail the pair limit, with TM within the limit. Both are internally repeated folds (a
propeller, a βα barrel), whose alignments can slip by a repeat or part of one at almost the same TM-score.
- **4CH9-2VWD** (Kelch against Nipah G, both six-bladed propellers): 0 % agreement, TM 0.681 (foldmap) against 0.678
  (US-align).
  - foldmap pairs each Kelch blade with the Nipah G blade one along from US-align's choice. B is shifted by 50-107
    residues, about one blade.
  - Under foldmap's own score, US-align's pairs give 0.678, slightly lower. The two registers are equally good
    superpositions of a pseudo-symmetric fold, so this is an inherent ambiguity, not an aligner error.
  - For `--like` it decides which blade mates which. The tuning run's H2/H3 failure of this pair (13a-10: blade
    blocks too large for A's places) may partly come from that register as well. It was not re-examined here.
- **1TIM-1IGS** (two TIM barrels): 79.6 % agreement, TM 0.605 against 0.620.
  - foldmap's pairs are shifted 2-4 residues (once 7) through A 16-54, the first βα units, and at A 241-246, the
    C-terminal end.
  - US-align's pairs score 0.620 under foldmap's TM-score as well. foldmap's search stopped 0.016 short here.

**The seven H1 refusals: US-align agrees with foldmap.**

| pair | foldmap TM | US-align TM (by A) | US-align TM by B |
|---|---|---|---|
| 1JBI_A-3U7Z_A | 0.246 | 0.243 | 0.247 |
| 1XU6_A-1PS2_A | 0.223 | 0.241 | 0.288 |
| 2D46_A-3F42_B | 0.352 | 0.380 | 0.309 |
| 2JCA_A-3LQ9_A | 0.369 | 0.394 | 0.384 |
| 2MJ7_A-2EHB_D | 0.484 | 0.492 | 0.540 |
| 2Q82_A-2JTR_A | 0.392 | 0.392 | 0.469 |
| 5GL3_A-1WPN_A | 0.431 | 0.446 | 0.598 |

- Normalised by the reference, as H1 is, US-align puts all seven under 0.5 (0.24-0.49). foldmap is lower by at most
  0.028, so the aligner does not under-score them.
- By H1's measure, CATH's topology label is generous for these cross-superfamily pairs, and the refusals stand.
- **Two qualifications.**
  - 2MJ7-2EHB is at 0.492, within 0.01 of the cut.
  - Normalised by B, the smaller chain, 2MJ7-2EHB (0.54) and 5GL3-1WPN (0.60) would reach 0.5. B's fold is largely
    contained in A's, and H1 asks whether B covers the reference.

**Proposed follow-up (its own task; not done here: align.py unchanged).**
- foldmap's structure search has three seeds: sequence, secondary structure, and gapless threading. TM-align has
  more:
  - a local-superposition seed (`get_initial5`: fragment superpositions scored by the DP);
  - a secondary-structure-plus-distance seed (`get_initial_ssplus`);
  - fragment gapless threading (`get_initial_fgt`).
- Adding the first two would likely close most of the 28 search-short cases. No H1 verdict on these sets would
  change: the nearest is 2MJ7-2EHB at 0.484, and US-align's best is 0.492.
- The change would alter `--like`, `compare` and `diff` on low-TM pairs, but not the bench, since `plot` does not
  align.
- To measure it: `homologues.py align --usalign` on both sets, with the search-short count, the median TM gap and the
  limits as the before/after numbers.
