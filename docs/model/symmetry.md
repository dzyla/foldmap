# Symmetry and copies

## Scope

- `symmetry.detect_symmetry`: the assembly's global symmetry (cyclic, dihedral, cubic, helical) or none, and its
  protomers.
- `symmetry._clusters`, `symmetry._orbits`, `symmetry._protomers`: copies grouped by sequence, orbits under the main
  axis, and protomer 0 built from one subunit of each orbit.
- `frame.symmetric_frame`, `frame.CylinderFrame`: the page frame of a symmetric assembly, side by side or unrolled.
- `cli._dyad_frame`, `cli._mirror_frame`: a two-fold assembly seen from its dyad, protomer 0 on the left or the copy
  the others mirror.
- `layout._separate`: copies kept alike while items are pushed apart.

## Assumes

- Two chains are copies when their aligned identity is at least `symmetry.SAME` (0.9). Chains shorter than
  `symmetry.MIN_RES` are left out.
- A rotation is a symmetry when it maps every clustered subunit onto a copy within `symmetry.TOLERANCE` (3.0 Å CA
  RMSD). Exact in cryo-EM models refined with symmetry; approximate in X-ray models.
- A screw is called helical under `auto` from `symmetry.MIN_HELIX_COPIES` (4) copies. Three copies on a screw are
  crystal packing (1URN) unless the user asks for helical.
- The user can ask for a symmetry, give the protomers, or switch symmetry off.

## Guarantees

- Exact cyclic and dihedral assemblies are found, with any axis:
  `tests/test_symmetry_scenarios.py::test_exact_cyclic_assemblies`,
  `tests/test_symmetry_scenarios.py::test_axis_need_not_be_z_or_through_the_origin`,
  `tests/test_symmetry_scenarios.py::test_dihedral_assemblies`.
- X-ray-like noise is tolerated: `tests/test_symmetry_scenarios.py::test_xray_like_noise_is_tolerated`.
- Translated copies are not symmetry: `tests/test_symmetry_scenarios.py::test_translated_copies_are_not_symmetry`.
- Three copies on a screw are helical only when asked:
  `tests/test_symmetry_scenarios.py::test_three_copies_on_a_screw_are_helical_only_when_asked`.
- The symmetry found does not depend on chain names:
  `tests/test_symmetry_scenarios.py::test_chain_names_do_not_matter`.
- A chain renaming that reverses the name order keeps the structure facts and the check verdicts:
  `tests/test_metamorphic.py::test_reversed_chain_names_keep_facts_and_verdicts`.
- Copies are drawn alike: `tests/test_copies.py::test_symmetric_controls_drawn_alike`.
- `symmetry.MIN_HELIX_COPIES` was stable at x0.75 and x1.25 (`docs/validation/sensitivity/report.md`) and is kept by
  the test suite (`docs/validation/prune/report.md`).

## Fails

- Chain names choose protomer 0 (metamorphic finding F1). `symmetry._orbits` walks the chains in name order and
  `symmetry._protomers` sorts the orbits by their first chain's name. Protomer 0 is drawn on the left (`cli._dyad_frame`)
  or is the copy the others mirror (`cli._mirror_frame`). On an unrolled frame its nucleic-acid stems are the template
  the copies' stems are tied to (`layout._stem_ties`). So a renaming that reverses the name order moves the figure of
  an assembly of several deposited chains: 1TIM by 29.05 page units, 1BL8 by 15.8, 1LMB by 0.75. Pinned as
  xfail(strict): `tests/test_metamorphic.py::test_reversed_chain_names_draw_the_same_figure`.
- Copies not drawn alike (`A15`) on the third random sample, at 7696947: 8PIN (focus view, B~C 18.88), 9PDH (default,
  B~C 9.68), 4YJ5 (focus, C~B 4.19) (`docs/validation/random100-third/review.md`). The causes are not investigated.
- With `layout._SEP_TIE_ROOM` at 0 (the second chance for copies drawn alike), 9FZL fails `A15`
  (`docs/validation/prune/report.md`). The constant is sensitive (`docs/validation/sensitivity/report.md`).
- A wrong symmetry call (a pseudo-symmetric complex called symmetric, or a real one missed at 3 Å) on an entry not
  among the tests: how often this happens is not measured.

## Seen by

- `A15`: copies drawn alike (deviation at most `criteria.LIMITS` `copies_alike`, 0.5 page units).
- `A12`: the same figure under any PYTHONHASHSEED, which catches an order that depends on set iteration.
- `D1`: the same figure drawn twice in one process.
- No check sees the dependence on chain names (finding F1). Only the xfail test above does; `foldmap check` and the
  campaigns draw each entry under its own names.
- No check sees a symmetry call that is wrong but still draws its copies alike.

## Evidence

- [metamorphic/report.md](../validation/metamorphic/report.md): finding F1 and the renaming gaps.
- [random100-third/review.md](../validation/random100-third/review.md): the `A15` failures.
- [sensitivity/report.md](../validation/sensitivity/report.md) and [prune/report.md](../validation/prune/report.md):
  `symmetry:MIN_HELIX_COPIES`, `layout:_SEP_TIE_ROOM`, `layout:_SEP_MIRROR_ROOM`.
