# The checks

## Scope

- `science`: the science checks `S1` to `S6` (`science.SCIENCE`), each with its kind (`science.KIND`): structural
  consistency (`S1`-`S4`, `S6`) or external validation (`S5`). `science.sheet_checks` reads the sheets as built.
- `figcheck`: the page checks of a drawn figure (`figcheck.figure_checks`) and the science checks read from the drawn
  figure (`figcheck.science_checks`).
- `criteria`: the page criteria `A1` to `A20` (`criteria.CRITERIA`), their limits (`criteria.LIMITS`), and `H1` to `H5`
  for a figure drawn like a reference (`criteria.LIKE_CHECKS`). `criteria.check_figure` judges one figure as
  `foldmap check` does.
- `floor`: `R1`, the readability floor (`floor.readability_floor`, `floor.CODE`).
- `deep`: `D1` to `D3` (`deep.DEEP`), run by `foldmap check --deep`, and the `--json` report (`deep.report`).
- `nacheck`: a nucleic-acid figure against the coordinates (`nacheck.THRESHOLDS`, `nacheck.judge`), judged by `A17`
  in the campaign.
- `fallback`: the named fallbacks counted per figure (`fallback.FALLBACKS`).
- Tests only: the metamorphic tests (`tests/test_metamorphic.py`) and independent kernels
  (`tests/test_independent_kernels.py`). They are not part of `foldmap check`.

## Assumes

- A check reads what `spec.draw_spec` returns (`spec.Drawn`: the figure, loops, layout, elements, backbone, model). Page
  measures are in page units; pixel measures use the figure as drawn.
- The science checks take foldmap's own DSSP bridges and pairing as the truth. They are internal: a figure true to a
  wrong DSSP call passes them.
- Some criteria need more than one process or an outside program: `A1`, `A1c`, `A2`, `A2b`, `A10`, `A12`, `A13`, `A14`,
  `A17`, `A19` (`criteria.NOT_HERE`). `foldmap check` names them not measured; only the random-100 campaign judges
  them. `D2` in `foldmap check --deep` takes `A10`'s measure.
- An unavailable check is never counted as passed. The `--json` verdict is "incomplete" when any is unavailable.

## Guarantees

- Each check names what it is evidence of: `tests/test_science.py::test_each_check_names_what_it_is_evidence_of`.
- `foldmap check` exits 1 on a failing figure criterion:
  `tests/test_criteria.py::test_check_exits_1_when_a_figure_criterion_fails`; an unmeasured criterion never prints as a
  pass: `tests/test_criteria.py::test_an_unmeasured_criterion_never_prints_as_a_pass`.
- The package and the campaign script judge one figure alike:
  `tests/test_criteria.py::test_the_script_and_the_package_judge_one_figure_alike`.
- An unavailable deep check is never counted as passed:
  `tests/test_deep_check.py::test_unavailable_is_never_counted_as_passed`.
- `R1` fails a page held small and passes an ordinary figure:
  `tests/test_floor.py::test_a_page_held_small_is_below_the_floor`,
  `tests/test_floor.py::test_an_ordinary_figure_is_above_the_floor`; time is a note, never a reason:
  `tests/test_floor.py::test_time_is_a_note_never_a_reason`.
- Every `fallback.note` call in the code is listed in the inventory with the check that reports it; the inventory
  test fails on a call the table lacks (`tests/test_inventory.py`, `docs/validation/inventory/inventory.md`).
- Planted defects are counted: `tests/test_na_validation.py::test_planted_crossings_are_counted`,
  `tests/test_science.py::test_a_strand_wedged_between_two_that_pair_fails_s1`.

## Fails

- `A2` is load dependent: two runs of one code on one machine gave 79 and 78 passes of 100 with `A2` judged, the same
  keys otherwise (`docs/validation/sensitivity/report.md`, stage B). The sensitivity study records time and does not
  judge it.
- `R1` changes exit codes: 6 of the 300 sample default figures (9FZL, 7JU4, 6JLN, 6ET5, 8WMV, 7NZ2) now exit 1 on `R1`
  alone; 8XJV and 8ZTU failed other checks before (`docs/validation/readability-floor/report.md`). That every sample
  default under 6 pt is among those measured is inferred from the samples' reports, not run at this code on all 300.
- `R1`'s route-budget reason fired on none of the figures the floor record measured and none of the samples' default
  views, where it is tested on synthetic measures only (`docs/validation/readability-floor/report.md`). It does fire
  past the default views: 7NZ2's domains view spends the page's budget, 6 searches refused
  (`docs/validation/prune/report.md`, held), so `foldmap check 7NZ2 --domains auto` fails `R1` on it. `R1`'s
  exit-code change reaches past the default views.
- `S5` is measured only where the committed CATH reference holds the entry. On most entries it is not measured.
- `D3` cannot read the committed .cif.gz fixtures (RNApolis fails on them); it ran on PDB copies
  (`docs/validation/metamorphic/report.md`).
- The H-check limits were frozen on the homologue tuning set; the three cross-family six-bladed propeller pairs fail
  them and stay listed failures (`criteria.LIMITS`, `docs/validation/homologues/tuning-report.md`).

## Seen by

A check's own failure is seen by its tests and by the records that measured it, not by another check.

The failures no check sees, collected from the other pages:

- Reading ([reading.md](reading.md)): a plausible but wrong assembly choice; which altloc is drawn.
- DSSP ([dssp.md](dssp.md)): an element split or joined by the kink or copy rules against what a reader would choose.
- Sheets ([sheets.md](sheets.md)): the register of a pairing on the page; a one-bridge pairing a reader would call a
  sheet (`sheets._MIN_PAIRS` decides); `S6` passing an accessory whose pairing is only named in the notes.
- Symmetry ([symmetry.md](symmetry.md)): chain names choose protomer 0 (metamorphic finding MF1; only an xfail test
  sees it); a wrong symmetry call that still draws its copies alike.
- Layout ([layout.md](layout.md)): the latent `layout._open_at_ends` bug under `--like`; a drawing that moves while
  every key stays equal (the prune record shows such moves); the length of a disulfide bar; text drawn under 6 pt by
  design (residue counts, base letters, ligand names).
- Routing ([routing.md](routing.md)): loops crossing loops (`figcheck.loop_crossings`), reported, not judged.
- Fallbacks counted and gated by no check (the inventory's check column): `route_unsnapped`, `glyph_width_guessed`,
  `helix_name_moved`; and those said in the notes only: `counts_hidden`, `mark_least_bad_spot`,
  `interactions_dropped`, `interactions_hidden`, `crosslinks_unplaced`. `counts_hidden` fires often: 15 in the
  sensitivity record's stage B baseline (`docs/validation/sensitivity/report.md`).
- Nucleic acids ([nucleic-acids.md](nucleic-acids.md)): na_run_crossings and na_protein_crossings, reported, not
  judged; angle_median, angle_max and joint_error_max on a plain frame, zero by construction
  (`nacheck.ZERO_BY_CONSTRUCTION`); 3EGZ's spearman per copy (not measured, too few stems); whether a pair taken under
  foldmap's wider windows is real, where RNApolis was not run.
- Membrane ([membrane.md](membrane.md)): the membrane placed against OPM's (not measured); a wrong cytoplasmic side.
- Domains ([domains.md](domains.md)): a wrong split on an entry with no CATH reference; the residue cuts the element
  rules undo.
- Every science check is internal: a figure true to a wrong DSSP call passes `S1` to `S4` and `S6`. Only `A10` and
  `D2` hold DSSP to mkdssp, and only `S5` holds anything to an outside chopping.

## Evidence

- [inventory/inventory.md](../validation/inventory/inventory.md): every fallback and the check that reports it.
- [sensitivity/report.md](../validation/sensitivity/report.md): `A2`'s load dependence; the method's limits.
- [prune/report.md](../validation/prune/report.md): identical keys are not identical drawings.
- [metamorphic/report.md](../validation/metamorphic/report.md): metamorphic findings MF1-MF4, kernels, `D2` and `D3`
  runs.
- [readability-floor/report.md](../validation/readability-floor/report.md): `R1` and its census.
- [homologues/tuning-report.md](../validation/homologues/tuning-report.md): the `H2` and `H3` limits.
- [random100-third/review.md](../validation/random100-third/review.md): every check on the third sample.
