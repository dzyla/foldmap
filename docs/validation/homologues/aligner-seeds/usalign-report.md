# foldmap's aligner against US-align

scripts/homologues.py align --usalign at 28ab8e2: foldmap's structure route (align_chains, match 'structure') against US-align 20241108 (`USalign A.pdb B.pdb -mol prot -mm 0 -ter 2`) on the same residues (each chain's Cα as foldmap reads it, --assembly asu). TM-scores normalised by A's chain (US-align: by Structure_1, A; d0 of A's length, as foldmap's). US-align's pairs: every column of its alignment holding a residue of both. Limits where US-align's TM >= 0.5: foldmap's TM within 0.05, and >= 85% of US-align's pairs reproduced within ±1 residue (13a-3).

Columns: US-align's TM-score; foldmap's (structure route) and the difference; foldmap's TM-score of US-align's own pairs (scoring apart from search); the auto route's TM-score as the campaign used it; the share of US-align's pairs foldmap reproduces within ±1 (exactly; of those under 5 Å); the limits ('–': US-align's TM under the limits' cut). What the numbers say (the misses, the H1 refusals, the aligner's search) is in review.md.

## tuning

| pair | kind | US-align TM | foldmap TM | diff. | foldmap TM of US pairs | auto route, TM | pairs ±1 | exact | under 5 Å ±1 | limits |
|---|---|---|---|---|---|---|---|---|---|---|
| 1EMA_A-1G7K_A | positive | 0.887 | 0.887 | +0.000 | 0.887 | structure, 0.887 | 1.000 | 1.000 | 1.000 | pass |
| 1HDC_A-1CYD_A | positive | 0.861 | 0.861 | +0.000 | 0.861 | structure, 0.861 | 1.000 | 1.000 | 1.000 | pass |
| 1MBN_A-1ECA_A | positive | 0.798 | 0.798 | +0.000 | 0.798 | structure, 0.798 | 1.000 | 1.000 | 1.000 | pass |
| 1REI_A-1TNF_A | negative | 0.414 | 0.419 | +0.005 | 0.414 | structure, 0.419 | 0.735 | 0.289 | 0.804 | – |
| 1SHG_A-1CKA_A | positive | 0.843 | 0.843 | -0.000 | 0.843 | structure, 0.843 | 1.000 | 1.000 | 1.000 | pass |
| 1TIM_A-1HDC_A | negative | 0.411 | 0.399 | -0.012 | 0.411 | structure, 0.399 | 0.577 | 0.182 | 0.549 | – |
| 1TIM_A-1IGS_A | positive | 0.620 | 0.620 | -0.001 | 0.620 | structure, 0.620 | 0.974 | 0.974 | 0.988 | pass |
| 1UBQ_A-1WM3_A | positive | 0.804 | 0.804 | +0.000 | 0.804 | structure, 0.804 | 1.000 | 1.000 | 1.000 | pass |
| 1UBQ_A-2ACY_A | negative | 0.322 | 0.353 | +0.031 | 0.322 | structure, 0.353 | 0.000 | 0.000 | 0.000 | – |
| 256B_A-1MBN_A | negative | 0.435 | 0.435 | +0.000 | 0.434 | structure, 0.435 | 0.878 | 0.851 | 0.944 | – |
| 2ACY_A-1URN_A | positive | 0.518 | 0.525 | +0.006 | 0.518 | structure, 0.525 | 0.926 | 0.778 | 0.972 | pass |
| 2VWD_A-2X9M_A | positive | 0.979 | 0.979 | +0.000 | 0.979 | sequence, 0.978 | 1.000 | 1.000 | 1.000 | pass |
| 2VWD_A-2ZB5_A | positive | 0.723 | 0.723 | +0.001 | 0.723 | structure, 0.723 | 0.987 | 0.933 | 0.993 | pass |
| 2ZB5_A-1E8T_A | positive | 0.698 | 0.698 | +0.000 | 0.698 | structure, 0.698 | 0.997 | 0.979 | 1.000 | pass |
| 4CH9_A-2VWD_A | positive | 0.678 | 0.694 | +0.015 | 0.678 | structure, 0.694 | 0.000 | 0.000 | 0.000 | miss |

- 11 of 15 pairs at US-align TM >= 0.5: TM within 0.05 11, pairs >= 85% 10, both **10**; pairs reproduced (±1): median 1.0, least 0.0.
- All 15: TM difference (foldmap - US-align) median 0.0, largest 0.0313; foldmap lower in 3; same side of H1's 0.5 in both: 15; foldmap's own score of US-align's pairs above its own alignment's (search short): 1; scoring difference on US-align's pairs at most 0.001.

### Miss: 4CH9_A-2VWD_A

Where foldmap's pairs part from US-align's (A's residues):

- A 300-301 (2): B shifted -38
- A 302-322 (21): not paired
- A 323-324 (2): B shifted -59
- A 325-335 (11): B shifted -63
- A 336-336 (1): B shifted -64
- A 337-356 (20): B shifted -60
- A 357-358 (2): B shifted -58
- A 359-359 (1): B shifted -57
- A 360-360 (1): B shifted -56
- A 361-371 (11): B shifted -65
- A 372-378 (7): B shifted -67
- A 379-380 (2): B shifted -68
- A 381-381 (1): B shifted -69
- A 382-382 (1): B shifted -73
- A 383-383 (1): B shifted -72
- A 384-394 (11): B shifted -70
- A 395-403 (9): B shifted -71
- A 404-404 (1): B shifted -101
- A 405-407 (3): B shifted -106
- A 408-416 (9): B shifted -97
- A 417-417 (1): B shifted -96
- A 418-418 (1): B shifted -90
- A 419-425 (7): B shifted -93
- A 426-427 (2): B shifted -92
- A 428-428 (1): B shifted -91
- A 429-429 (1): B shifted -89
- A 430-439 (10): B shifted -90
- A 440-441 (2): B shifted -91
- A 442-450 (9): B shifted -90
- A 453-453 (1): B shifted -54
- A 455-467 (13): B shifted -53
- A 468-468 (1): B shifted -50
- A 469-477 (9): B shifted -51
- A 478-478 (1): B shifted -49
- A 479-480 (2): B shifted -48
- A 481-488 (8): B shifted -67
- A 489-490 (2): B shifted -66
- A 491-497 (7): B shifted -69
- A 498-499 (2): B shifted -70
- A 500-501 (2): not paired
- A 502-509 (8): B shifted -72
- A 510-511 (2): B shifted -71
- A 515-515 (1): B shifted -68
- A 516-522 (7): B shifted -67
- A 523-525 (3): B shifted -66
- A 526-526 (1): B shifted -70
- A 527-535 (9): B shifted -52
- A 536-536 (1): B shifted -53
- A 537-537 (1): B shifted -56
- A 538-544 (7): B shifted -53
- A 545-547 (3): B shifted -52
- A 548-558 (11): B shifted -55
- A 559-561 (3): not paired
- A 562-563 (2): B shifted -58

```
usalign  A --------------------------------------PKVMIVVGGQA---P-KA-IRSVECYD-F-EE----DRWDQIAELPS-R--RCRAGVVFMAG
usalign  B NICLQKTSNQILKPKLISTLPVVGQSGTCITDPLLAMDEGYFAYSHLERIGSCGVSKQRIIGVGEVLDRGDEVPSLFMTNVWTPPNPNTVYHCSAVYNNE

usalign  A HVYAVGGFN-GSLR---------VRTVDVYDGVK-------DQWTSIA-SM-Q----ERRSTLGAAVLND-LLYAVGGFD--------------------
usalign  B FYYVLCAVSTVGDPILNSTYWSGSLMMTRLAVKPKSNGGGYNQHQLALRSIEKGRYDKVMPYGPSGIKQGDTLYFPAVGFLVRTEFKYNDSNCPITKCQY

usalign  A -----------G-----STGLASVEAYSYKT---NEWFFVAPMN--T-RRSSVGVGVV-EGKLYAVGGYDGASRQCLSTVEQYNPATN-EWIYVADMSTR
usalign  B SKPENCRLSMGIRPNSHYILRSGLLKYNLSDGENPKVVFIEISDQRLSIGSPSKIYDSLGQPVFYQASFS--W-DTMIKFGDVLTVNPLVVNWRNNTVIS

usalign  A R-------------------SGAGVGVLSG---QLYATGG-HDGPLVRKSVEVYDPGTNTWKQVADMNM----CRRNAGVCAV-N---GLLYVVGGDDG-
usalign  B RPGQSQCPRFNTCPEICWEGVYNDAFLIDRINWISAGVFLDSNQTAENPVFTVFK---DNEILYRAQLASEDTNAQKTITNCFLLKNKIWCISLVEIYDT

usalign  A --SCNLASVEYYNPVTDKWTLLPTNMSTGRSYAGVAVIHK-
usalign  B GDNVIRPKLFAVKIPEQC----------------------T

foldmap  A PKVMIVVGGQAPKAIRSVECYDF-------EEDRWDQIAELPSR------RCRAGVVFMAGHVYAVGGFN---GS-L-RVRTVDVYDGVK-----DQWTS
foldmap  B NI---------------------CLQKTSNQILKPKLISTLPVVGQSGTCITDPLLAMDEGYFAYSHLERIGSCGVSKQRIIGVGEVLDRGDEVPSLFMT

foldmap  A IASMQE-R--RSTLGAAVLNDLLYAVGGFD-GSTG---------LASVEAYSY-K------TNEWFFVA-PM-N----TRRSSVGVGVVEG-KLYAVGGY
foldmap  B NVWTPPNPNTVYHCSAVYNNEFYYVLCAVSTVGDPILNSTYWSGSLMMTRLAVKPKSNGGGYNQHQLALRSIEKGRYDKVMPYGPSGIKQGDTLYFPAVG

foldmap  A DGA----------------------------------SRQCLSTVEQYNPAT---NEWIYVADMS--T-RRSGAGVGVL-SGQLYATGGHDGPLVRKSVE
foldmap  B FLVRTEFKYNDSNCPITKCQYSKPENCRLSMGIRPNSHYILRSGLLKYNLSDGENPKVVFIEISDQRLSIGSPSKIYDSLGQPVFYQASFS--WDTMIKF

foldmap  A V-YDPGTN-TWKQVAD-MNMC------------------RRNAGVCAVNG---LLYVVGG-DDGSCNLASVEYYNPVTDKWTLLPTNMS---TGRSYAGV
foldmap  B GDVLTVNPLVVNWRNNTVISRPGQSQCPRFNTCPEICWEGVYNDAFLIDRINWISAGVFLDSNQTAENPVFTVFK---DNEILYRAQLASEDTNAQKTIT

foldmap  A AVI-H----K-----------------------------
foldmap  B NCFLLKNKIWCISLVEIYDTGDNVIRPKLFAVKIPEQCT

```

## heldout

| pair | kind | US-align TM | foldmap TM | diff. | foldmap TM of US pairs | auto route, TM | pairs ±1 | exact | under 5 Å ±1 | limits |
|---|---|---|---|---|---|---|---|---|---|---|
| 1CO6_A-1JDL_A | positive (close) | 0.883 | 0.883 | -0.000 | 0.883 | sequence, 0.822 | 1.000 | 1.000 | 1.000 | pass |
| 1G9R_A-1JYL_A | positive (remote) | 0.510 | 0.511 | +0.001 | 0.510 | structure, 0.511 | 0.979 | 0.926 | 0.987 | pass |
| 1JBI_A-3U7Z_A | positive (remote) | 0.243 | 0.271 | +0.028 | 0.243 | structure, 0.271 | 0.327 | 0.058 | 0.281 | – |
| 1JK4_A-2HNW_A | positive (close) | 0.894 | 0.895 | +0.001 | 0.894 | sequence, 0.894 | 0.987 | 0.962 | 0.987 | pass |
| 1LSH_B-5MLC_3 | negative (negative) | 0.192 | 0.193 | +0.001 | 0.192 | structure, 0.193 | 1.000 | 0.625 | 1.000 | – |
| 1NKP_A-1AM9_A | positive (remote) | 0.749 | 0.749 | +0.000 | 0.749 | structure, 0.749 | 1.000 | 1.000 | 1.000 | pass |
| 1PZW_A-1QHL_A | negative (negative) | 0.323 | 0.316 | -0.006 | 0.323 | structure, 0.316 | 0.583 | 0.062 | 0.676 | – |
| 1RPU_A-1ZWT_A | negative (negative) | 0.331 | 0.339 | +0.008 | 0.331 | structure, 0.339 | 0.639 | 0.454 | 0.714 | – |
| 1T5B_A-2Z9C_A | positive (close) | 0.977 | 0.977 | +0.000 | 0.977 | sequence, 0.977 | 1.000 | 1.000 | 1.000 | pass |
| 1TYY_A-3IQ0_A | positive (remote) | 0.888 | 0.888 | +0.000 | 0.888 | structure, 0.888 | 1.000 | 1.000 | 1.000 | pass |
| 1XU6_A-1PS2_A | positive (remote) | 0.241 | 0.264 | +0.024 | 0.239 | structure, 0.264 | 0.742 | 0.258 | 0.783 | – |
| 2C3E_A-1HZF_A | negative (negative) | 0.271 | 0.291 | +0.020 | 0.271 | structure, 0.291 | 0.408 | 0.269 | 0.448 | – |
| 2D46_A-3F42_B | positive (remote) | 0.380 | 0.383 | +0.003 | 0.380 | structure, 0.383 | 0.750 | 0.650 | 0.867 | – |
| 2EWR_A-2FCL_A | positive (close) | 0.996 | 0.996 | +0.000 | 0.996 | sequence, 0.996 | 1.000 | 1.000 | 1.000 | pass |
| 2EX5_A-4NTQ_A | negative (negative) | 0.164 | 0.165 | +0.001 | 0.164 | structure, 0.165 | 0.750 | 0.523 | 0.829 | – |
| 2F20_A-1ZN6_A | positive (remote) | 0.782 | 0.782 | +0.000 | 0.782 | structure, 0.782 | 1.000 | 1.000 | 1.000 | pass |
| 2H5F_A-4Z42_B | negative (negative) | 0.323 | 0.296 | -0.028 | 0.322 | structure, 0.296 | 0.469 | 0.225 | 0.455 | – |
| 2IAY_A-4XZD_A | negative (negative) | 0.330 | 0.329 | -0.001 | 0.329 | structure, 0.329 | 0.833 | 0.439 | 0.900 | – |
| 2IJE_S-3T6A_A | positive (remote) | 0.739 | 0.739 | -0.001 | 0.739 | structure, 0.739 | 0.991 | 0.986 | 1.000 | pass |
| 2JCA_A-3LQ9_A | positive (remote) | 0.394 | 0.395 | +0.001 | 0.394 | structure, 0.395 | 0.943 | 0.874 | 0.966 | – |
| 2LSR_A-4YL6_A | positive (remote) | 0.811 | 0.811 | +0.000 | 0.811 | structure, 0.811 | 1.000 | 1.000 | 1.000 | pass |
| 2MJ7_A-2EHB_D | positive (remote) | 0.492 | 0.492 | +0.000 | 0.492 | structure, 0.492 | 1.000 | 1.000 | 1.000 | – |
| 2OKF_A-2INB_A | positive (close) | 0.976 | 0.976 | +0.000 | 0.976 | sequence, 0.976 | 1.000 | 1.000 | 1.000 | pass |
| 2Q82_A-2JTR_A | positive (remote) | 0.392 | 0.397 | +0.005 | 0.393 | structure, 0.397 | 0.958 | 0.732 | 0.983 | – |
| 3BUU_A-3VDI_A | negative (negative) | 0.348 | 0.330 | -0.018 | 0.347 | structure, 0.330 | 0.354 | 0.277 | 0.500 | – |
| 3IE4_A-3AQX_A | positive (close) | 0.958 | 0.958 | +0.000 | 0.958 | sequence, 0.957 | 1.000 | 1.000 | 1.000 | pass |
| 3N4H_A-4JCU_A | positive (remote) | 0.742 | 0.742 | +0.000 | 0.742 | structure, 0.742 | 0.982 | 0.973 | 1.000 | pass |
| 3O7I_A-2O8I_A | positive (remote) | 0.724 | 0.724 | -0.001 | 0.725 | structure, 0.724 | 0.986 | 0.957 | 1.000 | pass |
| 3OQL_A-4WWJ_A | positive (close) | 0.992 | 0.992 | +0.000 | 0.992 | sequence, 0.988 | 1.000 | 1.000 | 1.000 | pass |
| 3VDJ_A-2JXW_A | negative (negative) | 0.244 | 0.260 | +0.016 | 0.244 | structure, 0.260 | 0.000 | 0.000 | 0.000 | – |
| 3WZ0_A-2AV5_A | positive (close) | 0.970 | 0.970 | +0.000 | 0.970 | sequence, 0.970 | 1.000 | 1.000 | 1.000 | pass |
| 4E2U_A-2LWY_A | positive (remote) | 0.701 | 0.701 | -0.000 | 0.701 | structure, 0.701 | 1.000 | 1.000 | 1.000 | pass |
| 4IPV_A-4N2P_A | negative (negative) | 0.309 | 0.303 | -0.006 | 0.309 | structure, 0.303 | 0.099 | 0.070 | 0.051 | – |
| 4U12_A-4NPX_A | positive (remote) | 0.766 | 0.768 | +0.002 | 0.766 | structure, 0.768 | 0.987 | 0.960 | 0.986 | pass |
| 5E4E_A-3TGX_B | positive (remote) | 0.618 | 0.618 | -0.000 | 0.618 | structure, 0.618 | 0.968 | 0.947 | 0.976 | pass |
| 5GL3_A-1WPN_A | positive (remote) | 0.446 | 0.446 | +0.000 | 0.446 | structure, 0.446 | 1.000 | 1.000 | 1.000 | – |
| 5GY2_C-2D4X_A | positive (remote) | 0.798 | 0.799 | +0.001 | 0.798 | structure, 0.799 | 0.988 | 0.969 | 0.987 | pass |
| 5KDI_A-4NT2_A | positive (remote) | 0.797 | 0.797 | +0.000 | 0.797 | structure, 0.797 | 1.000 | 1.000 | 1.000 | pass |
| 5LNK_A-6G2J_A | positive (close) | 0.768 | 0.767 | -0.001 | 0.768 | sequence, 0.698 | 0.991 | 0.991 | 1.000 | pass |
| 5U5S_A-2WP1_A | positive (close) | 0.939 | 0.939 | +0.000 | 0.939 | sequence, 0.937 | 1.000 | 1.000 | 1.000 | pass |

- 23 of 40 pairs at US-align TM >= 0.5: TM within 0.05 23, pairs >= 85% 23, both **23**; pairs reproduced (±1): median 1.0, least 0.9684.
- All 40: TM difference (foldmap - US-align) median 0.0, largest 0.0282; foldmap lower in 11; same side of H1's 0.5 in both: 40; foldmap's own score of US-align's pairs above its own alignment's (search short): 4; scoring difference on US-align's pairs at most 0.0015.

## The held-out positives H1 refuses

| pair | stratum | foldmap TM (auto) | US-align TM | US-align: one fold |
|---|---|---|---|---|
| 1JBI_A-3U7Z_A | remote | 0.271 | 0.243 | no |
| 1XU6_A-1PS2_A | remote | 0.264 | 0.241 | no |
| 2D46_A-3F42_B | remote | 0.383 | 0.380 | no |
| 2JCA_A-3LQ9_A | remote | 0.395 | 0.394 | no |
| 2MJ7_A-2EHB_D | remote | 0.492 | 0.492 | no |
| 2Q82_A-2JTR_A | remote | 0.397 | 0.392 | no |
| 5GL3_A-1WPN_A | remote | 0.446 | 0.446 | no |

US-align puts 7 of 7 under TM 0.5 too.
