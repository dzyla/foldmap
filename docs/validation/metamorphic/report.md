# Metamorphic tests and independent kernels

Plan 28 Task 7 (rulings Q34, Q35). Written by: foldmap 0.1.0, git ddc408cc31fe (asserted by runner, dirty False), code fingerprint 72fb15138b17e48a, scripts/metamorphic_record.py sha256 666c5be1bb47; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

What tests/test_metamorphic.py, tests/test_independent_kernels.py and tests/test_deep_check.py assert, measured on committed fixtures; the findings they pin with xfail(strict). The deep-check runs below are recorded runs (they need mkdssp and RNApolis), not tests.

## Findings

- **F1.** Chain names choose protomer 0: symmetry._orbits walks the chains in name order (symmetry.py:234), symmetry._protomers sorts the orbits by their first chain's name (symmetry.py:255) and breaks its ties by name (symmetry.py:259); protomer 0 is drawn on the left (cli._dyad_frame, cli.py:359) or is the copy the others mirror (cli._mirror_frame, cli.py:257). A renaming that reverses the chains' name order moves the figure of an assembly of several deposited chains; its structure facts and check verdicts stay equal.
- **F2.** Base-pair H-bond window: foldmap calls a canonical pair whose Watson-Crick-edge H-bonds reach 3.8 Å (nucleic.HBOND_MAX, nucleic.py:16, applied at :147), one stretched to 4.0 Å (HBOND_STRETCH, :17, applied at :149); the kernel holds two of them to 3.5 Å.
- **F3.** C1'-C1' window: foldmap accepts 9.0-12.0 Å (nucleic.C1_PAIR, nucleic.py:18, applied at :127 and :131); the kernel 9.5-11.5 Å.
- **F4.** Planarity: foldmap drops a pair whose base planes are more than 45 degrees apart (nucleic.PLANE_MAX, nucleic.py:19, applied at :156); the kernel has no planarity bound.

F1, the order-reversing renaming (largest geometry gap, page units):

| Fixture | gap | facts equal | verdicts equal |
|---|---|---|---|
| 1TIM.cif | 29.05 | True | True |
| 1LMB.cif | 0.75 | True | True |
| 1BL8.cif | 15.8 | True | True |

F2-F4, base pairs over 19 fixtures: both 820, only foldmap 54, only the kernel 5. By why: F2 (hbond) 44, F3 (c1) 6, F4 (plane) 5, glycosidic angle 1 (never alone), modified nucleotides 7 (outside the kernel's scope; not a finding).

- c1: 1AOI.cif.gz: I:11/J:282, 1AOI.cif.gz: I:28/J:265, 1HYS.cif.gz: E:867/F:887, 4LD0.cif.gz: D:12/F:22, 4YBB.cif.gz: AA:112/AA:330, 4YBB.cif.gz: AA:582/AA:760
- plane: 1IHF.cif.gz: C:-37/E:37, 4YBB.cif.gz: AA:456/AA:476, 4YBB.cif.gz: AA:457/AA:475, 4YBB.cif.gz: AA:80/AA:89, 4YBB.cif.gz: AA:82/AA:87
- lambda: 1AOI.cif.gz: I:25/J:268
- modified: 1EHZ.cif.gz: A:10/A:25, 1EHZ.cif.gz: A:30/A:40, 1EHZ.cif.gz: A:49/A:65, 2OEU.cif.gz: A:1/B:20, 2OEU.cif.gz: A:10/B:11, 4YBB-psi.cif.gz: DA:2590/DA:2604, 4YBB.cif.gz: AA:522/AA:527
- hbond: 1AOI.cif.gz: I:10/J:283, 1AOI.cif.gz: I:105/J:188, 1AOI.cif.gz: I:124/J:169, 1AOI.cif.gz: I:126/J:167, 1AOI.cif.gz: I:19/J:274, 1AOI.cif.gz: I:24/J:269, 1AOI.cif.gz: I:25/J:268, 1AOI.cif.gz: I:26/J:267, 1AOI.cif.gz: I:28/J:265, 1AOI.cif.gz: I:32/J:261, 1AOI.cif.gz: I:34/J:259, 1AOI.cif.gz: I:39/J:254, 1AOI.cif.gz: I:55/J:238, 1AOI.cif.gz: I:56/J:237, 1AOI.cif.gz: I:69/J:224, 1AOI.cif.gz: I:91/J:202, 1HYS.cif.gz: E:867/F:887, 1HYS.cif.gz: E:868/F:886, 1KX5.cif.gz: I:-15/J:15, 1KX5.cif.gz: I:-19/J:19, 1KX5.cif.gz: I:-31/J:31, 1KX5.cif.gz: I:-42/J:42, 1KX5.cif.gz: I:32/J:-32, 1KX5.cif.gz: I:39/J:-39, 1KX5.cif.gz: I:40/J:-40, 1KX5.cif.gz: I:42/J:-42, 1URN.cif.gz: R:19/R:2, 3CRX.cif.gz: C:18/D:19, 3CRX.cif.gz: C:20/F:17, 3CRX.cif.gz: C:23/F:14, 3CRX.cif.gz: C:24/F:13, 3CRX.cif.gz: C:32/F:5, 3CRX.cif.gz: D:13/E:24, 3CRX.cif.gz: D:14/E:23, 3CRX.cif.gz: D:17/E:20, 3CRX.cif.gz: D:5/E:32, 4OO8.cif.gz: B:84/B:94, 4YBB.cif.gz: AA:131/AA:231, 4YBB.cif.gz: AA:449/AA:486, 4YBB.cif.gz: AA:454/AA:479, 4YBB.cif.gz: AA:459/AA:473, 4YBB.cif.gz: AA:461/AA:471, 4YBB.cif.gz: AA:582/AA:760, 4YBB.cif.gz: AA:75/AA:95

## Transforms (largest geometry gap; facts and verdicts equal)

| Fixture | motion | rename | shuffle | pdb | turn |
|---|---|---|---|---|---|
| 1UBQ.cif | 0 | 0 | 0 | 0 | 0 |
| 1TIM.cif | 0 | 0 | 0 | 0 | 0 |
| 1C3W.cif | 0 | 0 | 0 | 0.0024 | 0 |
| 1LMB.cif | 0 | 0 | 0 | 0 | 0 |
| na/1EHZ.cif.gz | 0 | 0 | 0 | 0 | 0.00057 |

1C3W's pdb gap comes from its assembly: the PDB file's BIOMT operators are written to 6 decimals (tests/test_metamorphic.py: test_pdb_assembly_operators_are_rounded).

## Noise (largest change over three seeds)

| Fixture | sigma Å | DSSP codes | bridges (n) | base pairs (n) | TM residues |
|---|---|---|---|---|---|
| 1UBQ.cif | 0.05 | 0.0132 | 0.0 (16) | - (0) | - |
| 1UBQ.cif | 0.1 | 0.0263 | 0.0303 (16) | - (0) | - |
| 1TIM.cif | 0.05 | 0.0182 | 0.0161 (61) | - (0) | - |
| 1TIM.cif | 0.1 | 0.0445 | 0.0339 (61) | - (0) | - |
| 1C3W.cif | 0.05 | 0.009 | 0.0909 (6) | - (0) | 0.0 |
| 1C3W.cif | 0.1 | 0.018 | 0.0909 (6) | - (0) | 0.006 |
| 1LMB.cif | 0.05 | 0.0335 | 0.0 (0) | 0.0 (19) | - |
| 1LMB.cif | 0.1 | 0.0447 | 0.0 (0) | 0.0 (19) | - |
| na/1EHZ.cif.gz | 0.05 | - | 0.0 (0) | 0.0 (21) | - |
| na/1EHZ.cif.gz | 0.1 | - | 0.0 (0) | 0.0233 (21) | - |

## Kernels: backbone H-bonds and bridges (Kabsch & Sander 1983)

| Fixture | residues | kernel bonds | foldmap bonds | kernel only | foldmap only | bridges equal |
|---|---|---|---|---|---|---|
| 1UBQ | 76 | 126 | 54 | 72 | 0 | True (16) |
| 1TIM | 494 | 835 | 357 | 478 | 0 | True (61) |
| 1EMA | 226 | 384 | 169 | 215 | 0 | True (94) |
| 2LZM | 164 | 282 | 122 | 160 | 0 | True (7) |
| 1C3W | 222 | 402 | 191 | 211 | 0 | True (6) |
| 6ZS5 | 257 | 392 | 147 | 245 | 0 | True (73) |

Kernel only: the peptide neighbour's CO(i-1)...NH(i), which the energy formula puts below -0.5 kcal/mol and the DSSP program (and foldmap, dssp.py:69) skips.

## Deep checks, recorded runs

- D2 1UBQ: pass (1.0000 of 76 residues agree (limit 0.99))
- D2 1TIM: pass (1.0000 of 494 residues agree (limit 0.99))
- D2 1EMA: pass (1.0000 of 217 residues agree (limit 0.99))
- D2 2LZM: pass (1.0000 of 164 residues agree (limit 0.99))
- D2 1C3W: pass (1.0000 of 222 residues agree (limit 0.99))
- D2 6ZS5: pass (1.0000 of 257 residues agree (limit 0.99))
- D3 1EHZ: as committed unavailable (RNApolis could not read it: IndexError: string index out of range); as PDB pass (F1 1.000 (21 of 21 pairs, 21 in RNApolis; limit 0.95)); F2 pairs 0, not in RNApolis: none
- D3 2OEU: as committed unavailable (RNApolis could not read it: IndexError: string index out of range); as PDB pass (F1 1.000 (20 of 20 pairs, 20 in RNApolis; limit 0.95)); F2 pairs 0, not in RNApolis: none
- D3 4OO8: as committed unavailable (RNApolis could not read it: IndexError: string index out of range); as PDB pass (F1 0.988 (42 of 42 pairs, 43 in RNApolis; limit 0.95)); F2 pairs 1, not in RNApolis: none
- D3 1KX5: as committed unavailable (RNApolis could not read it: IndexError: string index out of range); as PDB pass (F1 1.000 (146 of 146 pairs, 146 in RNApolis; limit 0.95)); F2 pairs 8, not in RNApolis: none
- D3 1AOI: as committed unavailable (RNApolis could not read it: IndexError: string index out of range); as PDB pass (F1 0.986 (139 of 140 pairs, 142 in RNApolis; limit 0.95)); F2 pairs 16, not in RNApolis: I:28/J:265

Earlier RNApolis runs on whole downloaded entries, from older code: docs/validation/na/report.md.
