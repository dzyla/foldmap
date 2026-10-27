# Membrane

## Scope

- `membrane.find_membrane`: where the lipid bilayer is: its normal, its slab and its sides.
- `membrane._refit`: the normal refitted to the helices crossing the slab, until they stop changing.
- `membrane.crosses`, `membrane.slab_range`, `membrane.joins`: which elements cross the band, where a slab helix
  lies, which helix pieces are joined.
- `cli._face_medians`: a band from helices that only reach a face.
- In the membrane view, `layout.build_layout` places the band first: transmembrane elements at one height, the rest
  around them, each loop held by a wall (`route._walled_fallback`).

## Assumes

- When the file has OPM's dummy atoms, they give the membrane exactly: planes, thickness and sides.
- Otherwise it is estimated. An outer-membrane β-barrel (closed, oily outside) gives the normal as its axis and a 24 Å
  slab (`membrane.BARREL_HALF`). Else the long helices and strands give a first normal and a 30 Å slab
  (`membrane.HALF`) is slid to the most hydrophobic, least charged place; the normal is refitted to the crossing
  helices when there are at least `membrane.MIN_CROSSING_HELICES` (3).
- The normal comes from the elements' axes. No search over tilt angles is made.
- The cytoplasmic side comes from the positive-inside rule: more Lys and Arg on the loops of the crossing elements.
- A helix lies flat from 45 degrees to the normal, or from `membrane.FACE_FLAT_TILT` (35 degrees) within 7.5 Å of a
  face.
- The band holds only transmembrane helices, slab helices at their depths and crossing strands (the rule `A18`
  judges).

## Guarantees

- Membrane proteins are recognised and soluble ones are not:
  `tests/test_membrane.py::test_membrane_proteins_are_recognised`, `tests/test_membrane.py::test_soluble_proteins_are_not`,
  `tests/test_membrane.py::test_gfp_barrel_is_rejected_by_its_outward_residues`.
- OPM dummy atoms give the exact membrane: `tests/test_membrane.py::test_opm_dummy_atoms_give_the_exact_membrane`.
- The positive-inside rule puts bacteriorhodopsin's C terminus inside:
  `tests/test_membrane.py::test_positive_inside_rule_puts_bacteriorhodopsin_c_terminus_inside`; sides are labelled
  only when the rule is clear: `tests/test_membrane.py::test_sides_labelled_only_when_the_positive_inside_rule_is_clear`.
- Every transmembrane helix crosses the band in the membrane view:
  `tests/test_membrane.py::test_every_transmembrane_helix_crosses_the_band_in_the_membrane_view`.
- Under 0.1 Å noise, 1C3W's transmembrane residues change by at most 0.006 (three seeds,
  `docs/validation/metamorphic/report.md`).
- The membrane constants measured were stable at x0.75 and x1.25 on the 8 membrane entries of the third sample,
  except `membrane.FACE_FLAT_TILT` and `cli.MIN_TAU_GAIN_INV` (the view's turn about the normal), both sensitive
  (`docs/validation/sensitivity/report.md`, stage V).

## Fails

- 6DQN (third sample, at 7696947): no band drawn; `A18` fails (`docs/validation/random100-third/review.md`). In the
  sensitivity and prune runs, 7 of the 8 membrane entries pass and one fails `A18`.
- The membrane stage is no unseen set: the membrane constants were tuned in plan 26 with 6DQN, 6BO8, 6UW4 and 6QUE,
  4 of its 8 entries, among those looked at (`docs/validation/sensitivity/report.md`).
- Four constants (`membrane:_MIN_OUTWARD_HYDROPHOBIC`, `membrane:_MAX_OUTWARD_CHARGED`, `membrane:_LOOP_CAP`,
  `membrane:_SIDED_LOOPS`) were read in only one entry of the membrane stage, so their "stable" rests on that entry.
- Several constants were set for named entries: `membrane.MIN_CROSSING_HELICES` (7PIU, 7VVK),
  `membrane.FACE_FLAT_TILT` (32KS), `membrane.JOIN_GAP` (32KS, 6UW4), `membrane._MAX_POLAR_CROSSED` and
  `membrane._MANY_CROSSING` (6JO6). All are kept: turning each off moved membrane drawings or failed tests
  (`docs/validation/prune/report.md`).
- The estimated membrane against OPM's: not measured. No record compares the estimated normal, centre or thickness
  with OPM's over a set.
- A side called by the positive-inside rule on a protein where the rule fails: how often is not measured.

## Seen by

- `A18`: only the band set in the band (transmembrane helices across it, slab helices at their depths, crossing
  strands), no element nearer the mid-plane than its residues or on the wrong side, and loops held to their ports and
  faces. It judges the figure against foldmap's own membrane, not against OPM's.
- `A3`: a loop's fallback curve, including `loop_wall_detour` and `loop_wall_unkept`.
- Not seen by any check: a membrane placed wrongly but drawn consistently. `A18` holds the figure to the membrane
  foldmap found; nothing holds that membrane to OPM.
- Not seen by any check: a wrong side (cytoplasm on top when it is below).

## Evidence

- [sensitivity/report.md](../validation/sensitivity/report.md): stage V (membrane), its 8 entries and its limits.
- [prune/report.md](../validation/prune/report.md): the named-entry membrane constants, all kept.
- [random100-third/review.md](../validation/random100-third/review.md): 6DQN's `A18`.
- [metamorphic/report.md](../validation/metamorphic/report.md): the transmembrane residues under noise.
- [inventory/inventory.md](../validation/inventory/inventory.md): the membrane constants and their origins.
