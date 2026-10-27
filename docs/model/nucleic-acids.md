# Nucleic acids

## Scope

- `io.load_nucleic`: the DNA and RNA strands, their base pairs and the protein residues touching them.
- `nucleic.secondary_structure`: base pairs, stems, single-stranded runs and bends from the coordinates.
- `layout.build_layout`: stems placed along their projected axes; `na_curve` draws a bent duplex as one curve.
- `route_na.route_nucleic`: single strands, joints, free ends, gaps and their text (modules `na_runs`, `na_ends`,
  `na_labels`, `na_paths`).
- `cli._na_holders` and `layout._stem_ties`: on an unrolled symmetric frame, each strand drawn with the protomer
  holding it, its stems tied to protomer 0's copy (fix F2).
- `nacheck.shape`, `nacheck.contacts`, `nacheck.verify_contacts`, `nacheck.judge`: the drawing against the
  coordinates, judged against `nacheck.THRESHOLDS`.

## Assumes

- A canonical or wobble pair has Watson-Crick-edge H-bonds within `nucleic.HBOND_MAX` (3.8 Å), one of them stretched to
  `nucleic.HBOND_STRETCH` (4.0 Å) when the pair is in register on a strict pair; C1'-C1' within `nucleic.C1_PAIR`
  (9.0-12.0 Å); base planes within `nucleic.PLANE_MAX` (45 degrees).
- A stem is drawn along its projected 3D axis, on the side it is seen when its arms are at least `layout._SIDE_CLEAR`
  (2.75 Å) apart.
- A nucleotide touched by protein is marked with a bead in that chain's colour, backed by a heavy-atom pair within
  `io.CONTACT` (4.0 Å).
- A large nucleic acid may be cropped to the stems the protein touches (`cli._na_extent`); the notes say so.

## Guarantees

- Base pairs against an independent kernel on 19 fixtures: 820 found by both, 54 by foldmap only, 5 by the kernel
  only; every disagreement is explained by a named threshold or by modified nucleotides:
  `tests/test_independent_kernels.py::test_every_base_pair_disagreement_is_explained`,
  `tests/test_independent_kernels.py::test_the_findings_are_exactly_these`.
- `D3` recorded runs against RNApolis, as PDB: pair F1 1.000 (1EHZ, 2OEU, 1KX5), 0.988 (4OO8), 0.986 (1AOI)
  (`docs/validation/metamorphic/report.md`).
- Pairs agree with RNApolis on the committed fixtures: `tests/test_na_pairs.py::test_pairs_agree_with_rnapolis`.
- A fixture is drawn true to its structure: `tests/test_na_validation.py::test_fixture_is_drawn_true_to_its_structure`.
- Contacted stems stay nearer their touching elements:
  `tests/test_na_validation.py::test_contacted_stems_stay_nearer_their_touching_elements`.
- The NA fixes measured (`docs/validation/na-fixes/report.md`): F7, F2, F4 and F6 committed; no control of the 39 got
  worse from the base to the fix round's head. These are regression sets, not unseen ones.

## Fails

- F1 (pairs missed at low resolution) deferred past the freeze: 8XJV's pair recall, 3EGZ's A2.U64 (pair_f1 0.941) and
  6ZY7's stem_recall 0.667 stay failures.
- F3 (a stem inside the protein's projection pushed out) reverted: its target 8VBG barely moved and ten controls got
  worse. 8VBG keeps contact_auc 0.505 and stems_nearer_frac 0.5.
- F5 (a short duplex kept whole where projection makes the bend) reverted: it made control 1RVA worse. 6QIK keeps its
  fallback, crossing and `A3` strand y:8-9.
- 3EGZ, drawn unrolled per copy (F2): procrustes 0.572 per copy (limit 0.5) still fails. spearman is not measured per
  copy (3 stems, under the 5 a rank correlation needs); through one view of the whole assembly it was 0.347, a failure.
  It no longer fails only because it is not measured. Its angles (10.3, 10.4 degrees) pass only under the unrolled
  limit of 45; against the plain-frame limit 0.05 they fail.
- 8XJV: bead_overlaps went 125 to 129 with F4 and to 127 with F6, net +2 from the base; accepted because 8XJV is a
  case, not a control. Its fallbacks 10, na_crossings 10 and 125 bead overlaps are F1's.
- Pairs foldmap takes that a stricter kernel does not (metamorphic findings F2-F4): 44 by the H-bond window (to 3.8
  and 4.0 Å, the kernel 3.5), 6 by the C1'-C1' window (9.0-12.0 Å, the kernel 9.5-11.5), 5 by planarity (the kernel
  has no bound). On the `D3` runs RNApolis has 24 of the 25 H-bond-window pairs those entries hold; 1AOI I:28/J:265
  is foldmap's only. Whether the rest are real pairs is not measured.
- `nucleic.HBOND_STRETCH` is sensitive, and without it 7F5S fails `A17` and `A6` (`docs/validation/prune/report.md`).
- `layout._SIDE_CLEAR` is the one fragile constant: at x1.25 the third sample's pass count falls by 5 and `A17`
  failures go 7 to 13 (`docs/validation/sensitivity/report.md`). It was set on test structures (4YBB P5 at 2.81 Å is
  the closest stem that does not fit turned round).
- Tie labels on other ties, a collision the first two samples do not show: 3G6E, 6QIK, 7F5S (`A6`, third sample).
- Other `A17` failures on the third sample at 7696947: 6VTX (stems_nearer_frac, now not applicable after F7), 4AS1
  (na_crossings 1, 0 after F4) (`docs/validation/random100-third/review.md`).
- `D3` on the committed .cif.gz fixtures is unavailable (RNApolis cannot read them); it ran on PDB copies.

## Seen by

- `A17`: the nacheck measures past their `nacheck.THRESHOLDS`, and stem names by their letters (campaign only; it
  needs RNApolis). `foldmap check` names it not measured (`criteria.NOT_HERE`).
- `D3`: base pairs against RNApolis, when NA_RNAPOLIS_PYTHON names a Python with rnapolis.
- `A3`: strand fallbacks (`na_tail_straight`, `na_gap_unclear`, `na_hairpin_unclear`, `na_run_unclear`,
  `na_trace_unclear`).
- `A6`: tie labels and other text colliding.
- Not judged: na_run_crossings and na_protein_crossings are reported and have no limit (`nacheck.THRESHOLDS` holds
  None for them). 8XJV's 725 NA x protein crossings fail nothing.
- Not evidence of fidelity: angle_median, angle_max and joint_error_max are zero by construction on a plain frame
  (`nacheck.ZERO_BY_CONSTRUCTION`); they guard against a regression only.
- Not seen by any check: 3EGZ's spearman per copy (not measured, too few stems).
- Not seen by any check: whether a pair taken under foldmap's wider windows is real, on an entry RNApolis was not run
  on.

## Evidence

- [na-fixes/report.md](../validation/na-fixes/report.md): F2-F7 measured, F3 and F5 reverted, F1 deferred, the
  limitations.
- [metamorphic/report.md](../validation/metamorphic/report.md): findings F2-F4 and the `D3` runs.
- [sensitivity/report.md](../validation/sensitivity/report.md): `layout:_SIDE_CLEAR` fragile;
  `nucleic:HBOND_STRETCH`, `na_runs:_DETOUR_STEPS` and `na_runs:_BESIDE_DETOUR` sensitive.
- [prune/report.md](../validation/prune/report.md): `na_runs:_JOINT_NEAR` pruned; the others kept.
- [na/report.md](../validation/na/report.md): the earlier nucleic-acid campaign, older code.
- [random100-third/review.md](../validation/random100-third/review.md): the `A17` failures as first measured.
