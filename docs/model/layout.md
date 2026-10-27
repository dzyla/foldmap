# Layout

## Scope

- `layout.build_layout`: every element placed on the page, in page units (one unit is one strand pitch).
- `frame.view_frame`: the view; which 3D directions become page right, page up and the view axis.
- `layout._sheet_item` and `layout._sheet_order`: a sheet as one block, its strands left to right as seen.
- `layout._ba_template`, `layout._open_at_ends`, `layout.template_choices`: an all-parallel α/β sheet (TIM barrel,
  Rossmann fold) drawn as one βα row.
- `layout._separate` and `layout._attract`: items pushed apart until nothing overlaps, then pulled along their loops.
- `layout.page_fit`: the page size, and the scale when the page is at its largest.
- `like.like_from` and `like.Like`: a structure drawn like a reference (`--like`), its mates at the reference's places.

## Assumes

- The layout is a projection of the 3D structure through one view frame. No textbook template is imposed, except the
  βα row, which keeps each strand's projected order.
- Elements are drawn at their true length. A sheet is drawn face-on, never as its mirror image.
- A βα row is taken by a sheet that is all parallel, has 4 to `layout._BA_MAX` (10) strands, has an order that is not a
  guess (`layout.can_be_row`), and has at least 3/5 of its βα units in the row (`layout.template_applies`). Copies of a
  protein share the decision.
- A closed barrel of one chain drawn as a row is opened between its last and first strand in sequence, when those two
  pair (`layout._open_at_ends`).
- The page grows to keep the figure's text at `layout.MIN_TEXT_PT` (6 pt) or more, up to `layout.MAX_PAGE_IN` (48 in).
  Past that the whole figure is scaled down.

## Guarantees

- Every element is placed once, with no overlap: `tests/test_layout.py::test_every_sse_placed_once_without_overlap`.
- Elements are true to scale: `tests/test_layout.py::test_element_lengths_are_true_to_scale`.
- Sheets are face-on, never mirrored: `tests/test_layout.py::test_sheets_are_drawn_face_on_never_mirrored`.
- The layout is deterministic: `tests/test_layout.py::test_deterministic`; `A12` passes on every entry of the third
  sample with something to judge (`docs/validation/random100-third/review.md`).
- The bench's 22 entries draw the same layouts at the prune record's code (`docs/validation/prune/report.md`, bench
  --exact).
- A TIM barrel is one row, a helix between each pair of strands, opened between its last and first strand:
  `tests/test_template.py::test_tim_barrel_is_one_row_with_each_helix_between_its_strands`,
  `tests/test_template.py::test_template_barrel_opens_between_its_last_and_first_strand`.
- A structure drawn like itself is its own layout: `tests/test_like.py::test_like_itself_reproduces_the_layout`. A
  mated sheet reads as the reference's: `tests/test_like.py::test_a_mated_sheet_reads_as_the_reference`.
- Drawn like a reference, 16 of 30 held-out positive pairs pass `H1` to `H5` at c099787; 22 of 22 judged positives
  sit closer to the reference than drawn alone (`docs/validation/homologues/review.md`). With the order carried, 17 of
  30 (`docs/validation/homologues/topology-transfer/report.md`). The held-out pairs are no longer unseen.
- On the third sample at 7696947, `A5` (no overlaps) and `A7` (text on the page) pass on every entry
  (`docs/validation/random100-third/review.md`).

## Fails

- **Latent, not seen by any check: `layout._open_at_ends` and the reference's reading.**
  - What: in `layout.build_layout`, each sheet's `like.Like.order` entry (`read_as`) is computed by `like._order`
    along the sheet's own track. True means the reference draws the mates left to right in track order.
  - But `layout._ba_template` is handed `layout._open_at_ends(sh)`, not the sheet. For a closed one-chain barrel whose
    first and last strands in sequence are ring neighbours (and not already the track's two ends), that is the same
    ring cyclically shifted to start between them.
  - `layout._sheet_order` then applies the reference's reading to the shifted list. The row is opened at foldmap's
    seam in the reference's direction, but the mates no longer stand in the reference's left-to-right order: the
    reading wraps at the seam, or reads reversed there.
  - Trigger: a closed one-chain sheet drawn as a βα row, its first and last strands in sequence ring neighbours, drawn
    `--like` a reference that sets its order (two or more mates in one of the reference's sheets, in one order there).
  - It occurred in none of the 46 sheets with an order set on the homologue sets (25 tuning, 21 held-out;
    `docs/validation/homologues/topology-transfer/report.md`). No test reaches the trigger: 1TIM drawn like itself
    (`tests/test_like.py::test_like_itself_reproduces_the_layout`) has its ring shifted but no order set for it
    (looked at once for this page: its one sheet's `like.Like.order` entry is None).
  - Not fixed: a fix waits for an unseen case.
- 9NMN (third sample): helix α153's name lies on the loop into it, failing `A4` and `A6` in three views.
- 8OJZ/8QX1 (third sample): one state draws a helix the other does not, so its helix names shift by one letter; `A19`
  fails. 9YKJ/10UJ on the first sample is the same kind (`docs/validation/random100/review.md`).
- Large figures: the page reaches 48 in and the text falls under 6 pt. On the three samples' default views: 9FZL
  (2.0 pt), 6ET5 (2.4), 6JLN (4.5), 7NZ2 (5.2), 8XJV (5.4), 7JU4 (5.5), 8WMV (5.5) and 8ZTU (5.7)
  (`docs/validation/readability-floor/report.md`).
- A sheet seen nearly edge-on has its left-right order decided by a fraction of an ångström between its end strands
  (`layout._sheet_order`). Under `--like` the reference's order is carried (3WZ0-2AV5). On a default figure, how often
  a sheet sits on that margin is not measured.
- Constants: `layout._SIDE_CLEAR` (a nucleic-acid stem's side, see [nucleic-acids.md](nucleic-acids.md)) is the one
  fragile constant; at x1.25 the third sample's pass count falls by 5. Sensitive: `layout._STACK_PAD`,
  `layout._HELIX_LOOP`, `layout._STEM_ROOM`, `layout._SS_RESTART_PAIRS`, `layout._SEP_MIRROR_ROOM`,
  `layout._SEP_TIE_ROOM`, `layout._STACK_SHUT`, `layout._END_ROOM`, `layout._PART_BLOWUP`
  (`docs/validation/sensitivity/report.md`). One factor at a time, two values each; interactions are not measured.
- Pruned on evidence: `layout:_NAME_OPEN` (7NZ2). Its removal left every key equal but moved three domains-view
  drawings (5NF8, 7F5S, 9FZL) that no key orders (`docs/validation/prune/report.md`).

## Seen by

- `A5`: element or stem boxes overlapping.
- `A6`, `A7`: text collisions and text off the page.
- `A8`: labels unique per chain.
- `A11`: N to C order and sheet order as drawn.
- `A12`, `A13`, `D1`: the same figure under another hash seed, from a saved layout, and drawn twice.
- `A15`: copies drawn alike (see [symmetry.md](symmetry.md)).
- `A19`: two states named alike.
- `R1`: text under 6 pt on its page (`floor.readability_floor`).
- `H1` to `H5`: a figure drawn like a reference (`foldmap check --like` and the homologue campaign only).
- Not seen by any check: the `layout._open_at_ends` bug above. `A11` passes it (a closed barrel may read from any
  strand, either way round) and `S1` to `S3` pass it (the pairing is unchanged). If it fired, `H2` and `H3` are where it
  would show; whether they catch it is not measured.
- Not seen by any check: a drawing that moves while every key stays equal. The prune record shows such moves.
- Not seen by any check: the length of a disulfide bar.
- Not seen by any check: residue counts, base letters and ligand names under 6 pt. They are drawn at 0.6 to 0.9 of the
  figure's text by design; `R1` records their size and does not judge it.

## Evidence

- [homologues/topology-transfer/report.md](../validation/homologues/topology-transfer/report.md): `like.Like.order`,
  the 46 sheets, 3WZ0-2AV5.
- [homologues/review.md](../validation/homologues/review.md): `H1` to `H5` on the held-out pairs.
- [sensitivity/report.md](../validation/sensitivity/report.md): the layout constants by class.
- [prune/report.md](../validation/prune/report.md): `layout:_NAME_OPEN` pruned; what kept the others.
- [readability-floor/report.md](../validation/readability-floor/report.md): figures under the floor.
- [random100-third/review.md](../validation/random100-third/review.md): the third sample's failures by check.
