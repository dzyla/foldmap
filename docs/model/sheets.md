# Sheets

## Scope

- `sheets.build_sheets`: strands grouped into sheets by their DSSP bridges, each sheet's strands put in order side by
  side (its track), closed barrels found, and strands off the track kept as accessory strands.
- `model.Sheet`: the result: strands in track order, their directions, the pairing kind of each neighbour pair,
  `closed`, `ambiguous` and `accessory`.
- `science.sheet_checks`: the science checks on the sheets as built (`foldmap summary`).

## Assumes

- Two strands are neighbours when at least `sheets._MIN_PAIRS` (2) DSSP bridges join them. A single bridge joins
  nothing, except the one closing a barrel.
- A sheet's track is the longest path of the maximum spanning tree of its pairing graph (`sheets._order_tree`).
- A closed sheet is a cycle of at least `sheets._BARREL_MIN` (5) strands (`sheets._ring`). A shorter cycle is no
  barrel; it is drawn open.
- A strand off the track is an accessory strand. It is drawn beside the strand it hangs from, never between two strands
  that pair.
- Bridges come from [dssp.md](dssp.md). A bridge DSSP misses, the sheets miss.

## Guarantees

- A TIM barrel is closed by its single β7-β8 bridge, and ordinary sheets are not closed by single bridges:
  `tests/test_sheets.py::test_tim_barrel_is_closed_by_its_single_beta7_beta8_bridge`,
  `tests/test_sheets.py::test_single_bridges_do_not_close_ordinary_sheets`.
- FimD (9FWZ) has a 24-strand barrel and its pendant strand beside its partner:
  `tests/test_sheets.py::test_fimd_barrel_has_24_strands_and_its_pendant_strand_beside_its_partner`.
- OmpF has a 16-strand barrel: `tests/test_sheets.py::test_ompf_barrel_has_16_strands_its_n_terminal_strand_beside_a9_23`.
- A sheet drawn open holds no ring the builder finds:
  `tests/test_sheets.py::test_a_sheet_drawn_open_holds_no_ring_the_builder_finds`.
- A strand wedged between two strands that pair fails `S1`:
  `tests/test_science.py::test_a_strand_wedged_between_two_that_pair_fails_s1`.
- On the third random sample, at 7696947, `S1`, `S2`, `S3`, `S4` and `S6` found nothing in any view
  (`docs/validation/random100-third/review.md`).
- The bridges equal an independent kernel's on six fixtures
  (`tests/test_independent_kernels.py::test_bridges_ladders_and_sheets_are_the_papers`).

## Fails

- 9FWZ (FimD), the external review's case: a strand was wedged between two strands that pair. Fixed in plan 12; `S1`
  and `S6` were written for it.
- 8UTF: a strand stood beside the ghost of an unpaired strand. Fixed; held by
  `tests/test_science.py::test_8utf_a_strand_beside_an_unpaired_ghost_fails_s1`.
- Noise changes bridges: at 0.1 Å, 0.0909 of 1C3W's 6 bridges and 0.0339 of 1TIM's 61 move
  (`docs/validation/metamorphic/report.md`). A sheet with few bridges may then change shape. How often this changes a
  drawn sheet is not measured.
- `S3` stops counting pairing cycles at `science._MAX_CYCLES` and says so (the fallback `s3_cycles_capped`). On such a
  sheet the barrel count is not checked.
- `sheets._ring` takes the longest cycle that one pairing left out of the tree closes. A ring made of two such cycles
  may be missed (from the code: only the longest single such cycle is taken); `S3` flags it when drawn open. How
  often it happens is not measured.

## Seen by

- `S1`: strands drawn side by side must share a bridge; two strands in one column must not overlap in height.
- `S2`: a closed barrel's closure must pair.
- `S3`: a barrel's strand count is its pairing cycle; a sheet drawn open holds no cycle of 5 or more.
- `S6`: an accessory strand is drawn at its partner, or its pairing is named in the notes.
- `A11`: each sheet block holds one sheet's strands in that sheet's order (either way round; a closed barrel from any
  strand), neighbours pointing as their pairing has them.
- No check sees the register of a pairing (which residues pair) on the page.
- No check sees a pairing of one bridge that a reader would call a sheet; `sheets._MIN_PAIRS` decides it.
- `S6` passes an accessory whose pairing is named in the notes but not drawn beside it. The reader then has only the
  note.

## Evidence

- [random100-third/review.md](../validation/random100-third/review.md): the science checks over the third sample.
- [metamorphic/report.md](../validation/metamorphic/report.md): bridge kernels and noise.
- [inventory/inventory.md](../validation/inventory/inventory.md): `sheets._MIN_PAIRS` and `sheets._BARREL_MIN` and
  their origins.
