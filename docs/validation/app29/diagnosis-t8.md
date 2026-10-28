# Plan 29 Task 8: 26EL and 22MZ, a time-boxed diagnosis (ruling P29-7)

Written by hand at cf16ce9 (clean src). This task was read-only: no code was changed. These are the two triage rows
marked "not found" (`docs/validation/app29/triage.md`). Both entries are from the tune half.

The probes are in `docs/superpowers/work/app29/t8/` (gitignored), and their outputs are in `out/`. Each one drew the
figure the app's way (`app_jobs.draw_settings` with `random100.app_settings()`) on zyla-office through
`.superpowers/remote/zl`, then measured the Layout. Where a probe tried a change, it patched the code in memory only.

## 26EL: 66% of the page empty, three chains on a diagonal

**Pinned.** `cli._one_of_each` takes three copies of the same protein as three distinct sequences. They do not touch
each other, and each copy is drawn at its own projected place, 132 to 301 Å apart.

What was tried and measured (`p26el.py`, `out/26el.json`):

- The structure: 20 chains, all one protein (residues 1-682). Fifteen chains are modelled whole (675 residues). Five
  are modelled in part: F (1-291), M (292-682), N (1-289), E (115-682) and S.
- `symmetry.detect_symmetry` finds C2 with two protomers of 10 chains. Protomer 0 is A B C D M F L H O P.
- `cli._focus_plan('auto')` takes `cli._one_of_each(protomer 0)` = {A, F, M}. `_subunits` gives 20 chains / 3 =
  6.67, which rounds to 7. That is the note "1 of 7 subunits shown · C2 symmetry".
- `_one_of_each` groups chains by `difflib.SequenceMatcher` over the residues each chain has modelled, and requires a
  ratio above 0.9. The ratios are A-F 0.592, A-M 0.734 and F-M 0.098. F holds only the first half of the protein and
  M only the second half, so each one forms a group of its own. Each group gives the one chain of protomer 0 it has.
- These three chains do not touch. There are 0 CA pairs within 8 Å between A, F and M. F's contacts are with C and
  E; M's are with N and K.
- Centroid distances: A-F 132.3 Å, A-M 192.0 Å, F-M 301.2 Å. Since 132 + 192 = 324 ≈ 301, A lies almost on the line
  between F and M.
- The layout keeps each chain at its projected place (projected mode). That gives three islands on one line through
  the page: F at the bottom left, A in the middle, M at the top right (page units F x -30..-7, A -17..8, M 9..28).
- The page is 60.0 x 75.3 units. The chains' boxes (each element's centre ± 0.5) cover 31.6% of it, which matches
  the triage's 66% empty.
- The figure therefore draws one protein three times: chain A whole (58 elements) and two half copies, F (23) and
  M (36). `model.homolog_chains` already calls all three copies of one protein: the figure gives F and M the same
  labels as A (A-S, S27-S41).

Mechanism: `_one_of_each` compares the whole modelled sequence. A copy modelled in part therefore falls under 0.9 and
counts as a different protein. The rest of foldmap judges copies by `model.homolog_chains`: at least 20 residues at
the same numbers, 90% of them alike, with the numbering offset allowed. `_one_of_each` disagrees with it.

**Fix proposed, by cause:** in `_one_of_each`, put two chains in one group when `model.homolog_chains` gives them
the same root, as well as when the sequence ratio is above 0.9. This adds no constant and no entry branch; it reuses
the copy rule the labels already use.

Measured with the change patched in memory:

- Reach, from `p_reach.py` (`out/reach.json`). On the 25 tune entries with a symmetry, only 26EL's focus changes,
  from {A, F, M} to {A}. The other 24 keep their chains and their subunit count. The other 22 tune entries have no
  symmetry, and `_one_of_each` is not reached for them.
- 26EL after, from `p26el_after.py` (`out/26el-after.json`):
  - The note reads "1 of 20 subunits shown · C2 symmetry", which is correct.
  - The page goes from 60.0 x 75.3 units to 27.9 x 32.9.
  - The chain box covers 93.4% of the page instead of 31.6%.
  - Chain A keeps its 58 elements.

Risk:

- Any symmetric entry whose protomer holds a copy modelled in part draws one chain fewer, and its note says more
  subunits. That is intended, but the third sample, the random-100 samples and the bench were not measured.
- `homolog_chains` could merge two different proteins that share 20 or more residue numbers with 90% the same
  names, such as near-identical paralogs or constructs in one numbering. A ratio above 0.9 merges most of these
  already.
- The `_focus_subset` guests are unchanged in kind. For 26EL after the change, no partial chain was lent.

Next step if taken: the light loop's tune half, the bench and the third sample (ruling P29-13). After that, a test
pinning 26EL-like input: a protomer with a chain modelled in part gives one chain per protein.

## 22MZ: unpaired base letters off the backbone, a backbone crossing itself four times

**Pinned, in two parts. A partial fix is proposed for the second part only.** No letter is misplaced against its own
nucleotide: each one sits where `_place_letters` is designed to put it.

What was tried and measured (`p22mz.py`, `out/22mz.json`, `out/22mz.xy.json`; the letters mapped onto the before
PNG through the SVG's run paths, with a fit residual under 1e-6):

- Five stems: P1a with 7 pairs, P1c 9, PK2 5, and P1b and PK1 with 2 each. There are nine runs, and 26 unpaired
  letters.
- Every nucleotide point lies exactly on its own run (distance 0.00). 23 of 26 letters sit 1.00 `NT_ROW` from their
  nucleotide, which is `_place_letters`' first option. 0:31 sits at 1.24, and 0:50 and 0:51 at 1.7 (crowded). The
  median distance from a letter to its own strand is 1.00 `NT_ROW` (0.55 units, about 50 px in the app's PNG), the
  same offset the stem letters have from their rails.

What makes them read as floating off:

1. **Runs much longer than their nucleotides.**
   - `na_runs._drawing` builds a run's `nt_xy` with `_spaced(points, m)`: m points spread evenly over the whole
     routed length (na_runs.py:320 and 378). The length is set by where the stems sit. `need = (j - i) * NT_STEP` is
     only its least.
   - Here the stems of a pseudoknot sit at their projected places far apart. The runs take 1.43 to 9.1 units per
     nucleotide, against `NT_STEP` 0.55: run 2 has 6 nucleotides over 22.2 units, and run 4 has one over 9.1 units.
   - Each letter is then one isolated glyph beside a long bare line, 3 to 16 steps from the next. Nothing on the line
     marks the nucleotide, so the letter does not read as attached to it.
   - The same projected stem places make the four self-crossings. Those were not traced further.
2. **7 of 26 letters sit nearer another strand than their own:** 0:8, 0:31, 0:50, 0:51, 0:56, 0:57 and 0:67. For
   example, 0:57 (U, run 6) is 0.87 `NT_ROW` from its own strand and 0.35 from run 2's line, on run 2's far side.
   - Runs 2 and 6 are routed side by side about 0.65 `NT_ROW` (0.36 units) apart, so no letter fits between them.
     The side away from the stems' centre is taken by PK2's stem letters.
   - `_place_letters` takes the first option that is clear: `room` (letters, bodies, names) and `struck` (no line
     through the letter's ink plus the line clearance). Nothing in it asks whether the letter is nearer its own
     strand than another.
   - The side choice is also badly conditioned for run 6. Its centre test is `_dot(side, xy - centre)` with the
     stems' mean centre about 4 px off the run's line at 0:57, so the sign is decided by noise.

**Fix proposed for part 2, by cause:** in `_place_letters`, a spot counts as free only when no other strand's line is
nearer it than its own strand. This adds no constant and no entry branch. Patched in memory (`p22mz_after.py`,
`out/22mz-after.json`):

- The letters nearer another strand go from 7 to 3. The 3 left (0:50, 0:56, 0:67) had no option that passed, and
  fell back to the old choice.
- Moved: 0:8 to 1.24 `NT_ROW`, 0:57 to 1.7 and 0:51 to 1.24.

Risk:

- The fix is partial, measured on 22MZ alone.
- More letters may move out to 1.7-3.1 `NT_ROW` or to the rings, and `A6` and `A17` on the nucleic-acid entries are
  unmeasured. `_letters_off_loops`' second pass places letters through the same code.

**Part 1 has no fix short of a new arrangement:**

- Bunching a run's nucleotides at `NT_STEP` instead of spreading them is a design choice about what a run's length
  means.
- Moving the stems of a pseudoknot off their projected places goes against the projection rule.

Both are for the author.

Not re-checked here: the ligand boxes on 22MZ's lines belong to the triage's `render._place_mark` group (rank 4,
`docs/validation/app29/fix-ligand`).
