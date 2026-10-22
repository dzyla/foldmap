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

## TM-align's other seeds (plan 24 Task 1, post hoc)

Done in plan 24 (code 28ab8e2): `align.py` gained TM-align's local-superposition seed (`get_initial5`), its
secondary-structure-plus-distance seed (`initial_ssplus`), and the second gap value of TM-align's `DP_iter` (each
seed's refinement goes on without a gap penalty). Measured afterwards, on zyla-office with US-align 20241108.
**Nothing frozen changes**: the H limits, `heldout-report.md`, `heldout-summary.json`, `tuning-report.md`,
`usalign.json` and `usalign-report.md` stand as measured. The post-hoc tables, the same as plan 23's with the new
aligner, are in `aligner-seeds/usalign-report.md` and `aligner-seeds/usalign.json`
(`homologues.py usalign-report ... --out docs/validation/homologues/aligner-seeds --commit 28ab8e2`, no `--summary`).

**Against US-align** (before: the committed `usalign.json`, 7ac7c0e; a re-run reproduced it):

| | before | after |
|---|---|---|
| search short (US-align's pairs score > 0.005 higher under foldmap's own TM-score) | 28 of 55 | 5 of 55 |
| ... of those 28: foldmap now within 0.002 of US-align's TM, or above it | – | 23 |
| pairs whose TM-score dropped | – | 0 |
| foldmap's TM lower than US-align's | 43 of 55 | 14 of 55 |
| TM difference (foldmap - US-align), median | -0.005 | 0.000 |
| US-align TM ≥ 0.5: TM within 0.05 | 34 of 34 | 34 of 34 |
| ... both 13a-3 limits | 32 of 34 (tuning 9/11, held-out 23/23) | 33 of 34 (tuning 10/11, held-out 23/23) |
| ... pairs reproduced (±1), median: tuning; held-out | 0.978; 0.990 | 1.000; 1.000 |
| ... least pairs reproduced, held-out | 0.867 | 0.968 |
| the same side of H1's 0.5 as US-align | 55 of 55 | 55 of 55 |

- The gap-free rounds do most of it: with the two seeds alone (measured during the work) 2 of the 28 closed; with
  the gap-free rounds as well, 23.
- **1TIM-1IGS**, one of plan 23's two misses, now passes: 97.4 % of US-align's pairs (was 79.6 %), TM 0.620 against
  0.620.
- **4CH9-2VWD keeps its blade register**: 0 % agreement, B one blade along from US-align's choice, TM 0.694 against
  US-align's 0.678 (was 0.681). foldmap's register still scores higher under its own TM-score. Its H2/H3 drawn alike
  are unchanged: 2.10 / 5.22 and 0.948 (231 mate pairs), still failing.
- **Still short (5)**, all under TM 0.42 and four of them negatives (foldmap / US-align): 1TIM-1HDC 0.399 / 0.411,
  1PZW-1QHL 0.316 / 0.323, 2H5F-4Z42 0.296 / 0.323, 3BUU-3VDI 0.330 / 0.348, 4IPV-4N2P 0.303 / 0.309. TM-align's
  fragment gapless threading (`get_initial_fgt`) is the seed still missing.
- **Agreement lower on 8 pairs**, each now scoring at or above US-align's TM normalised by A, the score foldmap
  maximises (US-align searches with the shorter chain's d0): 1REI-1TNF 0.747 to 0.735, 1JK4-2HNW 1.000 to 0.987,
  2ACY-1URN 0.951 to 0.926, 1XU6-1PS2 0.839 to 0.742, 2C3E-1HZF 0.423 to 0.408, 2Q82-2JTR 0.986 to 0.958, 3N4H-4JCU
  1.000 to 0.982, 3VDJ-2JXW 0.387 to 0 (TM 0.260 against 0.244). The three judged by the limits (1JK4, 2ACY, 3N4H)
  still pass both.

**Time** (the structure route alone, one pair at a time, zyla-office, one thread): median 0.14 s before, 0.19 s
after; mean 0.20 s and 0.28 s; slowest 0.90 s and 1.55 s (2ZB5-1E8T, 424 x 446 residues). `test_budget`'s
500-residue bundle (3 s limit): 0.46 s before, 0.77 s after. The local seed's work is capped by a count
(`LOCAL_CELLS`: TM-align's own fragment starts up to about 520 residues a chain, fewer beyond), and a new candidate
is scored in full only when its quick score can win (`SCREEN`); the first three seeds' penalised refinements are
scored as before, so no TM-score drops.

**H1-H5, post hoc.** `homologues.py run` on both sets, with 7ac7c0e's `align.py` and with plan 24's. **No verdict
changes**: every pair passes and fails the same checks on both sets.

| | frozen | re-run, 7ac7c0e | re-run, plan 24 |
|---|---|---|---|
| held-out positives passing H1-H5 | 16 of 30 (c099787) | 17 of 30 | 17 of 30 |
| held-out negatives refused (H1) | 9 of 10 | 9 of 10 | 9 of 10 |
| held-out positives refused by H1 | 7 | 7 | 7 |
| held-out H2 median, drawn alike / independently | 0.11 / 1.89 | 0.11 / 1.89 | 0.11 / 2.05 |
| held-out H3 share, drawn alike / independently | 1.00 / 0.71 | 1.00 / 0.71 | 1.00 / 0.70 |
| tuning positives passing H1-H5 | 8 of 11 | 8 of 11 | 8 of 11 |
| tuning negatives refused (H1) | 4 of 4 | 4 of 4 | 4 of 4 |

- Frozen to 7ac7c0e is not the aligner: two close pairs on the sequence route (1JK4-2HNW's H3, 3WZ0-2AV5's H2 and
  H3) pass at 7ac7c0e after the layout work since c099787.
- The new aligner moves values, not verdicts. Drawn alike: 3N4H-4JCU's H2 0.09 / 0.16 to 0.04 / 0.12, 5GY2-2D4X's
  0.18 / 0.46 to 0.30 / 0.44 (5 mates, was 6). The independent baseline moves where the mates change (5GY2-2D4X,
  2Q82-2JTR, the refused and negative pairs): the held-out median 1.89 to 2.05.
- The seven H1 refusals stay under 0.5: 1JBI-3U7Z 0.271, 1XU6-1PS2 0.264, 2D46-3F42 0.383, 2JCA-3LQ9 0.395,
  2MJ7-2EHB 0.492 (now US-align's own 0.492), 2Q82-2JTR 0.397, 5GL3-1WPN 0.446.
- `foldmap compare`/`--like` on 2LZM against 1UBQ (tests): TM 0.33 with 2 of 13 elements matched (was 0.32, 3),
  still not one fold.
