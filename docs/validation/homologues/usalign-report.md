# foldmap's aligner against US-align

scripts/homologues.py align --usalign at c22ba5f: foldmap's structure route (align_chains, match 'structure') against US-align 20241108 (`USalign A.pdb B.pdb -mol prot -mm 0 -ter 2`) on the same residues (each chain's Cα as foldmap reads it, --assembly asu). TM-scores normalised by A's chain (US-align: by Structure_1, A; d0 of A's length, as foldmap's). US-align's pairs: every column of its alignment holding a residue of both. Limits where US-align's TM >= 0.5: foldmap's TM within 0.05, and >= 85% of US-align's pairs reproduced within ±1 residue (13a-3).

Columns: US-align's TM-score; foldmap's (structure route) and the difference; foldmap's TM-score of US-align's own pairs (scoring apart from search); the auto route's TM-score as the campaign used it; the share of US-align's pairs foldmap reproduces within ±1 (exactly; of those under 5 Å); the limits ('–': US-align's TM under the limits' cut). What the numbers say (the misses, the H1 refusals, the aligner's search) is in review.md.

## tuning

| pair | kind | US-align TM | foldmap TM | diff. | foldmap TM of US pairs | auto route, TM | pairs ±1 | exact | under 5 Å ±1 | limits |
|---|---|---|---|---|---|---|---|---|---|---|
| 1EMA_A-1G7K_A | positive | 0.887 | 0.887 | -0.000 | 0.887 | structure, 0.887 | 0.995 | 0.995 | 0.995 | pass |
| 1HDC_A-1CYD_A | positive | 0.861 | 0.857 | -0.004 | 0.861 | structure, 0.857 | 0.992 | 0.957 | 0.996 | pass |
| 1MBN_A-1ECA_A | positive | 0.798 | 0.795 | -0.003 | 0.798 | structure, 0.795 | 0.978 | 0.963 | 0.978 | pass |
| 1REI_A-1TNF_A | negative | 0.414 | 0.405 | -0.010 | 0.414 | structure, 0.405 | 0.747 | 0.373 | 0.804 | – |
| 1SHG_A-1CKA_A | positive | 0.843 | 0.843 | -0.000 | 0.843 | structure, 0.843 | 1.000 | 1.000 | 1.000 | pass |
| 1TIM_A-1HDC_A | negative | 0.411 | 0.387 | -0.024 | 0.411 | structure, 0.387 | 0.565 | 0.171 | 0.549 | – |
| 1TIM_A-1IGS_A | positive | 0.620 | 0.605 | -0.016 | 0.620 | structure, 0.605 | 0.796 | 0.628 | 0.820 | miss |
| 1UBQ_A-1WM3_A | positive | 0.804 | 0.804 | +0.000 | 0.804 | structure, 0.804 | 1.000 | 1.000 | 1.000 | pass |
| 1UBQ_A-2ACY_A | negative | 0.322 | 0.337 | +0.015 | 0.322 | structure, 0.337 | 0.000 | 0.000 | 0.000 | – |
| 256B_A-1MBN_A | negative | 0.435 | 0.417 | -0.018 | 0.434 | structure, 0.417 | 0.878 | 0.595 | 0.944 | – |
| 2ACY_A-1URN_A | positive | 0.518 | 0.512 | -0.007 | 0.518 | structure, 0.512 | 0.951 | 0.741 | 1.000 | pass |
| 2VWD_A-2X9M_A | positive | 0.979 | 0.979 | +0.000 | 0.979 | sequence, 0.978 | 1.000 | 1.000 | 1.000 | pass |
| 2VWD_A-2ZB5_A | positive | 0.723 | 0.711 | -0.012 | 0.723 | structure, 0.711 | 0.952 | 0.822 | 0.963 | pass |
| 2ZB5_A-1E8T_A | positive | 0.698 | 0.691 | -0.007 | 0.698 | structure, 0.691 | 0.960 | 0.773 | 1.000 | pass |
| 4CH9_A-2VWD_A | positive | 0.678 | 0.681 | +0.003 | 0.678 | structure, 0.681 | 0.000 | 0.000 | 0.000 | miss |

- 11 of 15 pairs at US-align TM >= 0.5: TM within 0.05 11, pairs >= 85% 9, both **9**; pairs reproduced (±1): median 0.9779, least 0.0.
- All 15: TM difference (foldmap - US-align) median -0.0037, largest 0.0242; foldmap lower in 11; same side of H1's 0.5 in both: 15; foldmap's own score of US-align's pairs above its own alignment's (search short): 7; scoring difference on US-align's pairs at most 0.001.

### Miss: 1TIM_A-1IGS_A

Where foldmap's pairs part from US-align's (A's residues):

- A 16-16 (1): B shifted -4
- A 17-19 (3): not paired
- A 20-29 (10): B shifted -4
- A 32-34 (3): B shifted -2
- A 45-45 (1): B shifted -7
- A 46-54 (9): B shifted -3
- A 99-99 (1): not paired
- A 140-140 (1): not paired
- A 154-154 (1): not paired
- A 201-202 (2): not paired
- A 240-240 (1): not paired
- A 241-246 (6): B shifted -2

```
usalign  A ------------------------------------APRK----FFVGGNWKMN-----GKRKSLGELIHTLDGAKLSADTEVVCGAP--------SIYL
usalign  B PRYLKGWLKDVVQLSLRRPSFRASRQRPIISLNERILEFNKRNITAIIAEYKRKSPSGLDVERDPIEYSKFME--RYA---VGLSILTEEKYFNGSYETL

usalign  A DFARQKLDAKIGVAAQN-CYKVPKGAFTGEISPAMIKDIGAAWVILGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFQET
usalign  B RKIASS-VS-IPILMKDFIVKE--------SQIDDAYNLGADTVLLIV---K---ILTERELESLLEYARSYGMEPLIEIND----------E-NDLDIA

usalign  A KAIADNVKDWSKVVLAYEPVWA-IGTGKTATPQQAQEVHEKLRGWLKTHVSDAVAVQSRIIYGGSVTGGNCKELASQHDVDGFLVGGASL-KPEFVDIIN
usalign  B LRIG--A-----RFIGINSRDLETLEIN-KE--NQRKLISMIP-----------SN-VVKVAESGISERNEIEELRKLGVNAFLIGSSLMRNPEKIKEFI

usalign  A AKH
usalign  B L--

foldmap  A ------------------------------------APRK-----FFVGGNWKMNGKRK---SLGELIHTLDGAKLSADTEVVCGAPS----IYLDFARQ
foldmap  B PRYLKGWLKDVVQLSLRRPSFRASRQRPIISLNERILEFNKRNITAIIAEYKRKSP---SGLDVERDPIEYSKFMERYAVGLSILTEEKYFNGSYETLRK

foldmap  A K--LDAKIGVAAQNCYKVPKGAFTGEISPAMIKDIGAAWVILGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFQETKAIA
foldmap  B IASSVSIPILMKDFIVKE--------SQIDDAYNLGADTVLLIVK------ILTERELESLLEYARSYGMEPLIEIND-----------ENDLDIALRIG

foldmap  A DNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKTHVSDAVAVQSRIIYGGSVTGGNCKELASQHDVDGFLVGGASLKPEFVDIINAKH
foldmap  B A-------RFIGINSRDLETLEINKE--NQRKLISMIPS------------NVVKVAESGISERNEIEELRKLGVNAFLIGSSLMRNP-EKIKEFIL

```

### Miss: 4CH9_A-2VWD_A

Where foldmap's pairs part from US-align's (A's residues):

- A 300-322 (23): not paired
- A 323-324 (2): B shifted -59
- A 325-335 (11): B shifted -63
- A 336-336 (1): B shifted -64
- A 337-356 (20): B shifted -60
- A 357-360 (4): B shifted -57
- A 361-371 (11): B shifted -66
- A 372-372 (1): B shifted -73
- A 373-378 (6): B shifted -66
- A 379-380 (2): B shifted -67
- A 381-381 (1): B shifted -68
- A 382-383 (2): B shifted -72
- A 384-394 (11): B shifted -70
- A 395-403 (9): B shifted -71
- A 404-404 (1): B shifted -102
- A 405-407 (3): B shifted -107
- A 408-417 (10): B shifted -98
- A 418-418 (1): B shifted -90
- A 419-428 (10): B shifted -93
- A 429-429 (1): B shifted -89
- A 430-439 (10): B shifted -90
- A 440-441 (2): B shifted -91
- A 442-450 (9): B shifted -90
- A 453-453 (1): B shifted -54
- A 455-467 (13): B shifted -53
- A 468-468 (1): B shifted -50
- A 469-477 (9): B shifted -51
- A 478-480 (3): B shifted -48
- A 481-490 (10): B shifted -67
- A 491-497 (7): B shifted -70
- A 498-511 (14): B shifted -71
- A 515-515 (1): B shifted -68
- A 516-525 (10): B shifted -67
- A 526-527 (2): B shifted -71
- A 528-535 (8): B shifted -52
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

foldmap  A PKVMIVVGGQAPKAIRSVECYDF---------EEDRWDQIAELPSR------RCRAGVVFMAGHVYAVGGFN----GSLRVRTVDVYDGVKD-------Q
foldmap  B -----------------------NICLQKTSNQILKPKLISTLPVVGQSGTCITDPLLAMDEGYFAYSHLERIGSCGVSKQRIIGVGEVLDRGDEVPSLF

foldmap  A WTSIASMQER--RSTLGAAVLNDLLYAVGGFDGSTG---------LASVEAYSYK--------TNEWFFVAPMN------TRRSSVGVGVVEG-KLYAVG
foldmap  B MTNVWTPPNPNTVYHCSAVYNNEFYYVLCAVSTVGDPILNSTYWSGSLMMTRLAVKPKSNGGGYNQHQLALRSIEKGRYDKVMPYGPSGIKQGDTLYFPA

foldmap  A GYDGA----------------------------------SRQCLSTVEQYNPAT---NEWIYVADMS---TRRSGAGVGVLSGQLYATGGHDGPLVRKSV
foldmap  B VGFLVRTEFKYNDSNCPITKCQYSKPENCRLSMGIRPNSHYILRSGLLKYNLSDGENPKVVFIEISDQRLSIGSPSKIYDSLGQPVFYQASFSWDTMIKF

foldmap  A EVYDPGTN-TWKQVADMNMCR-------------------RNAGVCAVNG---LLYVVGG-DDGSCNLASVEYYNPVTDKWTLLPTNMS---TGRSYAGV
foldmap  B GDVLTVNPLVVNWRNNTVISRPGQSQCPRFNTCPEICWEGVYNDAFLIDRINWISAGVFLDSNQTAENPVFTVFK---DNEILYRAQLASEDTNAQKTIT

foldmap  A AVIHK----------------------------------
foldmap  B NCFLLKNKIWCISLVEIYDTGDNVIRPKLFAVKIPEQCT

```

## heldout

| pair | kind | US-align TM | foldmap TM | diff. | foldmap TM of US pairs | auto route, TM | pairs ±1 | exact | under 5 Å ±1 | limits |
|---|---|---|---|---|---|---|---|---|---|---|
| 1CO6_A-1JDL_A | positive (close) | 0.883 | 0.882 | -0.001 | 0.883 | sequence, 0.822 | 0.980 | 0.980 | 0.990 | pass |
| 1G9R_A-1JYL_A | positive (remote) | 0.510 | 0.502 | -0.009 | 0.510 | structure, 0.502 | 0.926 | 0.741 | 0.967 | pass |
| 1JBI_A-3U7Z_A | positive (remote) | 0.243 | 0.246 | +0.003 | 0.243 | structure, 0.246 | 0.058 | 0.058 | 0.094 | – |
| 1JK4_A-2HNW_A | positive (close) | 0.894 | 0.894 | -0.000 | 0.894 | sequence, 0.894 | 1.000 | 1.000 | 1.000 | pass |
| 1LSH_B-5MLC_3 | negative (negative) | 0.192 | 0.190 | -0.003 | 0.192 | structure, 0.190 | 0.896 | 0.625 | 0.919 | – |
| 1NKP_A-1AM9_A | positive (remote) | 0.749 | 0.742 | -0.006 | 0.749 | structure, 0.742 | 0.987 | 0.949 | 0.987 | pass |
| 1PZW_A-1QHL_A | negative (negative) | 0.323 | 0.284 | -0.038 | 0.323 | structure, 0.284 | 0.000 | 0.000 | 0.000 | – |
| 1RPU_A-1ZWT_A | negative (negative) | 0.331 | 0.324 | -0.007 | 0.331 | structure, 0.324 | 0.105 | 0.035 | 0.107 | – |
| 1T5B_A-2Z9C_A | positive (close) | 0.977 | 0.977 | +0.000 | 0.977 | sequence, 0.977 | 1.000 | 1.000 | 1.000 | pass |
| 1TYY_A-3IQ0_A | positive (remote) | 0.888 | 0.887 | -0.002 | 0.888 | structure, 0.887 | 0.997 | 0.979 | 0.996 | pass |
| 1XU6_A-1PS2_A | positive (remote) | 0.241 | 0.223 | -0.017 | 0.239 | structure, 0.223 | 0.839 | 0.452 | 0.913 | – |
| 2C3E_A-1HZF_A | negative (negative) | 0.271 | 0.283 | +0.012 | 0.271 | structure, 0.283 | 0.423 | 0.238 | 0.433 | – |
| 2D46_A-3F42_B | positive (remote) | 0.380 | 0.352 | -0.028 | 0.380 | structure, 0.352 | 0.275 | 0.000 | 0.300 | – |
| 2EWR_A-2FCL_A | positive (close) | 0.996 | 0.996 | +0.000 | 0.996 | sequence, 0.996 | 1.000 | 1.000 | 1.000 | pass |
| 2EX5_A-4NTQ_A | negative (negative) | 0.164 | 0.148 | -0.016 | 0.164 | structure, 0.148 | 0.636 | 0.341 | 0.714 | – |
| 2F20_A-1ZN6_A | positive (remote) | 0.782 | 0.773 | -0.009 | 0.782 | structure, 0.773 | 0.975 | 0.905 | 0.975 | pass |
| 2H5F_A-4Z42_B | negative (negative) | 0.323 | 0.291 | -0.032 | 0.322 | structure, 0.291 | 0.408 | 0.204 | 0.394 | – |
| 2IAY_A-4XZD_A | negative (negative) | 0.330 | 0.316 | -0.013 | 0.329 | structure, 0.316 | 0.818 | 0.379 | 0.860 | – |
| 2IJE_S-3T6A_A | positive (remote) | 0.739 | 0.731 | -0.008 | 0.739 | structure, 0.731 | 0.967 | 0.929 | 0.989 | pass |
| 2JCA_A-3LQ9_A | positive (remote) | 0.394 | 0.369 | -0.025 | 0.394 | structure, 0.369 | 0.736 | 0.391 | 0.759 | – |
| 2LSR_A-4YL6_A | positive (remote) | 0.811 | 0.806 | -0.005 | 0.811 | structure, 0.806 | 1.000 | 0.961 | 1.000 | pass |
| 2MJ7_A-2EHB_D | positive (remote) | 0.492 | 0.484 | -0.008 | 0.492 | structure, 0.484 | 0.954 | 0.874 | 0.963 | – |
| 2OKF_A-2INB_A | positive (close) | 0.976 | 0.976 | +0.000 | 0.976 | sequence, 0.976 | 1.000 | 1.000 | 1.000 | pass |
| 2Q82_A-2JTR_A | positive (remote) | 0.392 | 0.392 | -0.001 | 0.393 | structure, 0.392 | 0.986 | 0.775 | 0.983 | – |
| 3BUU_A-3VDI_A | negative (negative) | 0.348 | 0.267 | -0.081 | 0.347 | structure, 0.267 | 0.000 | 0.000 | 0.000 | – |
| 3IE4_A-3AQX_A | positive (close) | 0.958 | 0.958 | +0.000 | 0.958 | sequence, 0.957 | 1.000 | 1.000 | 1.000 | pass |
| 3N4H_A-4JCU_A | positive (remote) | 0.742 | 0.738 | -0.004 | 0.742 | structure, 0.738 | 1.000 | 0.929 | 1.000 | pass |
| 3O7I_A-2O8I_A | positive (remote) | 0.724 | 0.722 | -0.003 | 0.725 | structure, 0.722 | 0.986 | 0.971 | 1.000 | pass |
| 3OQL_A-4WWJ_A | positive (close) | 0.992 | 0.992 | +0.000 | 0.992 | sequence, 0.988 | 1.000 | 1.000 | 1.000 | pass |
| 3VDJ_A-2JXW_A | negative (negative) | 0.244 | 0.188 | -0.056 | 0.244 | structure, 0.188 | 0.387 | 0.097 | 0.348 | – |
| 3WZ0_A-2AV5_A | positive (close) | 0.970 | 0.970 | +0.000 | 0.970 | sequence, 0.970 | 1.000 | 1.000 | 1.000 | pass |
| 4E2U_A-2LWY_A | positive (remote) | 0.701 | 0.695 | -0.005 | 0.701 | structure, 0.695 | 0.978 | 0.948 | 0.992 | pass |
| 4IPV_A-4N2P_A | negative (negative) | 0.309 | 0.260 | -0.049 | 0.309 | structure, 0.260 | 0.000 | 0.000 | 0.000 | – |
| 4U12_A-4NPX_A | positive (remote) | 0.766 | 0.763 | -0.004 | 0.766 | structure, 0.763 | 0.867 | 0.867 | 0.903 | pass |
| 5E4E_A-3TGX_B | positive (remote) | 0.618 | 0.609 | -0.009 | 0.618 | structure, 0.609 | 0.926 | 0.874 | 0.964 | pass |
| 5GL3_A-1WPN_A | positive (remote) | 0.446 | 0.431 | -0.015 | 0.446 | structure, 0.431 | 0.641 | 0.571 | 0.722 | – |
| 5GY2_C-2D4X_A | positive (remote) | 0.798 | 0.785 | -0.013 | 0.798 | structure, 0.785 | 0.877 | 0.749 | 0.907 | pass |
| 5KDI_A-4NT2_A | positive (remote) | 0.797 | 0.793 | -0.004 | 0.797 | structure, 0.793 | 0.990 | 0.969 | 0.989 | pass |
| 5LNK_A-6G2J_A | positive (close) | 0.768 | 0.759 | -0.009 | 0.768 | sequence, 0.698 | 0.927 | 0.899 | 0.969 | pass |
| 5U5S_A-2WP1_A | positive (close) | 0.939 | 0.936 | -0.003 | 0.939 | sequence, 0.937 | 0.991 | 0.973 | 1.000 | pass |

- 23 of 40 pairs at US-align TM >= 0.5: TM within 0.05 23, pairs >= 85% 23, both **23**; pairs reproduced (±1): median 0.9897, least 0.8667.
- All 40: TM difference (foldmap - US-align) median -0.0057, largest 0.0805; foldmap lower in 32; same side of H1's 0.5 in both: 40; foldmap's own score of US-align's pairs above its own alignment's (search short): 21; scoring difference on US-align's pairs at most 0.0015.

## The held-out positives H1 refuses

| pair | stratum | foldmap TM (auto) | US-align TM | US-align: one fold |
|---|---|---|---|---|
| 1JBI_A-3U7Z_A | remote | 0.246 | 0.243 | no |
| 1XU6_A-1PS2_A | remote | 0.223 | 0.241 | no |
| 2D46_A-3F42_B | remote | 0.352 | 0.380 | no |
| 2JCA_A-3LQ9_A | remote | 0.369 | 0.394 | no |
| 2MJ7_A-2EHB_D | remote | 0.484 | 0.492 | no |
| 2Q82_A-2JTR_A | remote | 0.392 | 0.392 | no |
| 5GL3_A-1WPN_A | remote | 0.431 | 0.446 | no |

US-align puts 7 of 7 under TM 0.5 too.
