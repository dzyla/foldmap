# Secondary structure (DSSP) and elements

## Scope

- `dssp.assign_dssp`: Kabsch-Sander DSSP as mkdssp 3.x does it. Backbone H-bond energies (`dssp._hbond_energy`),
  n-turns, helices, bridges, ladders and beta-bulges, to one code per residue.
- `ss.build_sses`: the codes turned into elements (`model.SSE`): strands (E), helices (H) and, as an option, short
  3-10 helices (G). Each element gets an axis and a centroid.

## Assumes

- The backbone of [reading.md](reading.md): N, CA, C, O per residue; chain breaks where the peptide bond is longer than
  2.5 Å.
- An H-bond is an energy below `dssp._HBOND_E` (-0.5 kcal/mol), rounded and floored as mkdssp does.
- The code priorities of mkdssp 2.1 and later: H replaces E and B; G and I need every residue free, and I may replace
  H.
- Bends (S) are not assigned.
- A helix run is split where two alpha segments meet at more than `ss.KINK_SPLIT` (30 degrees) through 3-10 or pi
  residues. Strand pieces split by a bulge are joined when their bridges agree (`ss._join_bulged_pieces`).
- Copies of one protein are drawn alike: a kink split, a bulge join or a helix break made in one copy is made in the
  others (`ss._copy_breaks`).

## Guarantees

- Residue by residue equal to mkdssp 3.1.2 on every structure in tests/data with a reference (21), except the
  residues next to a chromophore mkdssp skips: `tests/test_dssp_reference.py::test_agrees_with_mkdssp`.
- `A10` on the first random sample: mkdssp agrees on 268,834 of 268,835 residues
  (`docs/validation/random100/review.md`).
  The one residue is in 8K8E (0.9992 of 1318, `docs/validation/random100/report.json`).
- `D2` recorded runs: 1.0000 agreement on 1UBQ, 1TIM, 1EMA, 2LZM, 1C3W and 6ZS5
  (`docs/validation/metamorphic/report.md`).
- Bridges, ladders and sheets equal those of an independent kernel written from the 1983 paper, on six fixtures:
  `tests/test_independent_kernels.py::test_bridges_ladders_and_sheets_are_the_papers`. The kernel finds more H-bonds
  (the CO(i-1) to NH(i) neighbour, which DSSP skips); the bridges are equal.
- Exact under rigid motion, renaming, atom shuffling and the PDB writer:
  `tests/test_metamorphic.py::test_structure_facts_are_exact`.
- Copies stay alike: `tests/test_ss.py::test_a_kink_split_in_one_copy_is_split_in_every_copy`,
  `tests/test_ss.py::test_a_helix_broken_in_one_copy_is_broken_in_every_copy`.

## Fails

- Noise moves codes. At 0.1 Å Gaussian noise on every coordinate, the largest change over three seeds is 0.0445 of
  1TIM's codes and 0.0447 of 1LMB's (`docs/validation/metamorphic/report.md`). Whether mkdssp
  moves as much under the same noise is not measured.
- 1EMA (GFP): foldmap keeps the chromophore CRO66 in the chain (`io._substitute_backbone`); mkdssp breaks the chain
  there. The residues within four of it are not compared with mkdssp.
- 8K8E: one residue of 1318 disagrees with mkdssp. The cause is not investigated.
- The kink rule and the copy rules are foldmap's, not DSSP's: a helix may be drawn in pieces DSSP keeps whole, or whole
  where one copy's DSSP breaks it. They are said in the figure's notes. Whether they match an author's choice is not
  measured.

## Seen by

- `A10`: residue agreement with mkdssp 3.1.2 at `criteria.LIMITS` `dssp_agree` 0.99 (campaign only; it needs mkdssp).
- `D2`: the same measure in `foldmap check --deep`, when the MKDSSP environment variable names mkdssp.
- `S1`, `S2`, `S3`: a missed or extra bridge shows when a drawn sheet disagrees with DSSP's pairing (see
  [sheets.md](sheets.md)).
- No check sees an element split or joined by the kink or copy rules against what a reader would choose. The notes
  name them; nothing judges them.
- No check sees DSSP's noise sensitivity on a single entry. It is measured on five fixtures only.

## Evidence

- [metamorphic/report.md](../validation/metamorphic/report.md): kernels, noise and the `D2` runs.
- [random100/review.md](../validation/random100/review.md): `A10` over the first sample.
- [random100-third/review.md](../validation/random100-third/review.md): `A10` passes on every entry with something to
  judge.
- [inventory/inventory.md](../validation/inventory/inventory.md): the DSSP constants (origin: biology, DSSP's
  energies).
