# Pruning on evidence

Plan 28 Task 9.3 (ruling Q46). Written by: foldmap 0.1.0, git b730a05051f4, code fingerprint 1a96dacb9a78f068, scripts/prune.py sha256 59997af73e91; Python 3.12.9 on Linux-6.18.33.1-microsoft-standard-WSL2-x86_64-with-glibc2.39; gemmi 0.7.5, numpy 2.5.3, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

The review asked whether foldmap's special cases "may have emerged from fixing individual examples". Each special case the inventory (`docs/validation/inventory/inventory.md`) says was set for a named entry, and each constant the sensitivity study found not exercised, was turned off by one value, one at a time (`scripts/perturb.py`'s hook rewrites the module as it is imported, so every reader takes the value and no code is edited), and measured against an unperturbed baseline: the third random sample's campaign keys in the default view (pass, failing criteria, named fallbacks, figure area), the view that reads it where the default view does not (sensitivity's stage V entries), the nucleic-acid cases and controls of plan 28's NA fixes, the domain sets (dom_cal, dom_held, dom_fresh: each chain's overlap, F1 and recall against CATH), `scripts/bench.py --exact`, and the full test suite with the value installed in every process. **A special case is removed only when every measure is identical or better and no test fails** (`scripts/prune.py`'s docstring defines each verdict); otherwise it is kept, and the table says which key or test kept it. The third sample has been run and looked at since plan 15: it is no unseen sample. The fourth sample was not touched.

**Limits.** One value per special case, one at a time: no interaction is measured. Turning off is not always the code as it was before the special case (each value's reason says what it turns off); three values are no faithful revert: `layout:_SEP_MIRROR_ROOM` (the extra room only; the extra push passes stay), `layout:_STACK_SHUT` (the nudge, the move to the partner's other side and the pocket go together) and `render:_SPOT_FAR` (it also narrows `_least_bad_spot`'s rings). A view other than the default is measured only for a constant the default view never reads, on the entries sensitivity's stage V drew it in. A test that runs foldmap in a child process (the CLI, random100) runs it unperturbed. **Identical means identical keys and checks** (each entry's pass, failing criteria, named fallbacks and figure area; the NA measures; the domain scores; the bench's layouts; the tests), not identical drawings: a drawing can move where no key sees it. Each campaign cell names the entries whose SVG moved (no verdict reads it), and each prune below says whether its drawings stayed the same. Identical keys say nothing of an entry no measure holds, such as the special case's named entry.

## Counts

| | candidates |
|---|---|
| all | 45 |
| tried (one value each, 1 measured with another's value) | 32 |
| kept | 28 |
| pruned | 2 |
| met the condition, held | 2 |
| incomplete | 0 |
| not tried: needs code removal | 13 |

Kept by measure (a candidate kept by several counts in each): domains 2, membrane 3, na 4, r100 7, suite 27.

## Baselines

- r100: 100 entries, 82 pass (A2, wall time, not judged); the two baselines identical (keys and SVGs).
- membrane: 8 entries, 7 pass.
- domains: 49 entries, 42 pass.
- na: 46 entries; failing A17: 1YTB, 30QM, 3CRX, 3EGZ, 3KTQ, 6QIK, 6ZY7, 8BW5, 8VBG, 8XJV, 9K4F, 9YLM.
- domsets: dom_cal 22 chains, dom_held 11 chains, dom_fresh 11 chains.
- bench --exact against `.superpowers/sdd/2026-10-17-15-heldout-failures/bench-before.json`: exit 0, 22 entries, changed: none.
- suite: 4619 tests, 0 failed, 56 skipped, exit 0.

zl time: r100 0.69 h, membrane 0.01 h, domains 0.18 h, na 0.39 h, domsets 0.00 h, bench 0.01 h, suite 0.88 h.

## Pruned

- `layout:_NAME_OPEN` (removed from the code): keys and checks identical on every measure; live logic: its removal moved drawings no key orders, domains 3 SVGs (5NF8, 7F5S, 9FZL); SVGs identical on r100, na.
- `na_runs:_JOINT_NEAR` (removed from the code): keys and checks identical on every measure; SVGs identical too on r100, na: no measured drawing took it.

### The pruned code draws what was measured

The runs measured foldmap code 4eb4cad74365b1ed at 46f5ebf, before the prune. The same measures' baselines run again at c64ea7a (code 8f878e5bd79166bf), with the pruned special cases removed from the code, against each pruned candidate's "-off" run:

- r100: `layout:_NAME_OPEN`-off 100/100 SVGs byte-identical; `na_runs:_JOINT_NEAR`-off 100/100 SVGs byte-identical; moved from the measured baseline: none.
- domains: `layout:_NAME_OPEN`-off 49/49 SVGs byte-identical; moved from the measured baseline: 5NF8, 7F5S, 9FZL.
- na: `layout:_NAME_OPEN`-off 46/46 SVGs byte-identical; `na_runs:_JOINT_NEAR`-off 46/46 SVGs byte-identical; moved from the measured baseline: none.

## Met the condition, held

- `route:_MAX_NODES_PAGE`: it bounds the time a page's second route searches take (7QUP's focus view took 158 s of its 180), and no measure judges time; identical says these sets never spend the budget. Its named entry 7QUP's focus view drew identically without it, but 7NZ2's domains view (the held-out sample, which spends it: 6 searches refused) did not: A4 and A6 newly failed, loop_over_text 0 -> 1, helix_name_moved 0 -> 1 (one random100 run each way at 46f5ebf, outside Q46's measures). Removing it also removes route_second_search_refused and the readability floor's route_budget reason

## Named entries looked at (informational, not a measure of Q46)

- layout:_NAME_OPEN, 7NZ2 (held-out sample) domains view: keys identical (pass, failing A17 only, fallbacks, area); its SVG differs (the name's place)
- na_runs:_JOINT_NEAR, 9C75 default view: identical, SVG too (9C75 is also an NA control, measured above)
- route:_MAX_NODES_PAGE, 7QUP focus view: identical, SVG too, 161.7 s both; 7NZ2 domains view: worse (see held)

## Tried

Each cell: the verdict, how many moves, the first.

| candidate | origin | value -> off | decision | kept by | r100 | view | na | domsets | bench | suite |
|---|---|---|---|---|---|---|---|---|---|---|
| `layout:_NAME_OPEN` | named-entry:7NZ2 | _NAME_ROUTE + 0.1 -> 0.15 | **prune** | - | identical | domains: identical; 3 SVGs moved (5NF8, 7F5S, 9FZL) | identical | identical, not imported | identical | identical |
| `na_runs:_JOINT_NEAR` | named-entry:9C75 | 2.5 -> 0.0 | **prune** | - | identical | - | identical | identical, not imported | identical | identical |
| `align:SAME_FOLD_TM` | measured | 0.5 -> 0.0 | **keep** | suite | identical | - | identical | identical, not imported | identical, not imported | worse (3: tests.test_align::test_auto_takes_the_sequence_only_when_its_pairs_superpose; …) |
| `align:SAME_FOLD_FRACTION` | measured | 0.5 -> 0.0 | **keep** | suite | identical | - | identical | identical, not imported | identical, not imported | worse (1: tests.test_align::test_same_fold_needs_the_tm_score_and_the_elements) |
| `io:PAIRED_MIN` | named-entry:4PLX,3OWW | 3 -> 1000000000 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (4: tests.test_assembly::test_a_lattice_contact_of_two_base_pairs_is_no_pairing; …) |
| `io:_PAIR_DIST` | named-entry:3OWW | 3.4 -> 0.0 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (4: tests.test_assembly::test_a_lattice_contact_of_two_base_pairs_is_no_pairing; …) |
| `layout:_TAIL_COUNT_NEAR` | named-entry:7SMR | 1.5 -> 1000000000.0 | **keep** | suite | better (3: 6QIK: fallback counts_hidden 1 -> 0; …); 3 SVGs moved (6QIK, 7F5S, 8XJV) | - | better (2: 6QIK: named_fallbacks 54 -> 53; …); 2 SVGs moved (6QIK, 8XJV) | identical, not imported | identical | worse (1: tests.test_text_clear::test_tail_count_sits_clear_beyond_its_letter[2.0]) |
| `layout:_SEP_MIRROR_ROOM` | named-entry:9GEP | 1.10 -> 0.0 | **keep** | na, r100, suite | moved (1: 3W96: area 3042.22 -> 3115.51); 1 SVG moved (3W96) | - | worse (4: 3W96: spearman 0.892 -> 0.891; …); 2 SVGs moved (30QM, 3W96) | identical, not imported | identical | worse (2: tests.test_snapshots::test_figure_matches_snapshot[30QM-ncrna]; …) |
| `layout:_SEP_TIE_ROOM` | named-entry:8QEM,6JLN,7QUP,9P1X | 2.0 -> 0.0 | **keep** | r100, suite | worse (5: 9FZL: PASS -> FAIL; …); 2 SVGs moved (9FZL, 9Y70) | - | identical | identical, not imported | identical | worse (12: tests.test_app_controls::test_two_fold_choice_is_offered_for_a_two_fold_assembly; …) |
| `layout:_STACK_SHUT` | named-entry:9P1X,8OZJ | 0.6 -> 0.0 | **keep** | r100, suite | worse (4: 7U9Z: PASS -> FAIL; …); 2 SVGs moved (7U9Z, 9NMN) | - | identical | identical, not imported | identical | worse (3: tests.test_domains::test_domain_names_and_leaders_leave_every_loop_its_way[8BCU]; …) |
| `layout:_STACK_NUDGE` | named-entry:9P1X,8OZJ | 1.0 -> 0.0 | **keep** | r100 | moved (1: 9NMN: area 10388.35 -> 9556.46); 2 SVGs moved (7U9Z, 9NMN) | - | identical | identical, not imported | identical | identical |
| `layout:_POCKET_REACH` | named-entry:8BCU,12BM | 2.0 -> 0.0 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (2: tests.test_domains::test_domain_names_and_leaders_leave_every_loop_its_way[8BCU]; …) |
| `layout:_END_ROOM` | named-entry:9FWZ | END_LABEL - 0.15 -> 0.0 | **keep** | r100, suite | worse (7: 6DQN: PASS -> FAIL; …); 5 SVGs moved (2AIQ, 5ES7, 6DQN, …) | - | identical | identical, not imported | identical | worse (5: tests.test_figcheck::test_accessories_named_not_drawn_level_are_reported[2OMF-A-3-3]; …) |
| `layout:_BA_MAX` | named-entry:2BNH | 10 -> 1000000000 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (2: tests.test_template::test_a_row_holds_at_most_ten_strands[10-True]; …) |
| `layout:_PART_BLOWUP` | named-entry:1WWL,2ILU,7SMR,6JLN,5LDW | 2.0 -> 1000000000.0 | **keep** | domains, suite | identical | domains: moved (1: 8U1N: area 1331.46 -> 8059.7); 1 SVG moved (8U1N) | identical | identical, not imported | identical | worse (5: tests.test_copy_labels::test_no_chain_draws_one_label_twice[8UTF-auto-none-auto]; …) |
| `membrane:_MAX_POLAR_CROSSED` | named-entry:6JO6 | 0.20 -> 0.15 | **keep** | suite | identical | membrane: identical | identical | identical, not imported | identical | worse (3: tests.test_membrane::test_a_slab_crossed_by_many_long_elements_may_be_more_polar; …) |
| `membrane:_MANY_CROSSING` | named-entry:6JO6 | 10 -> 1000000000 | **keep** | suite | identical | membrane: identical | identical | identical, not imported | identical | worse (3: tests.test_membrane::test_a_slab_crossed_by_many_long_elements_may_be_more_polar; …) |
| `membrane:MIN_CROSSING_HELICES` | named-entry:7PIU,7VVK | 3 -> 1000000000 | **keep** | membrane, suite | identical | membrane: moved (7: 6BO8: area 2905.01 -> 3003.33; …); 7 SVGs moved (6BO8, 6DQN, 6QUE, …) | identical | identical, not imported | identical | worse (7: tests.test_membrane::test_a_soluble_helix_bundle_beside_the_membrane_does_not_tilt_its_normal; …) |
| `membrane:FACE_FLAT_TILT` | named-entry:32KS | 35.0 -> 45.0 | **keep** | membrane, suite | identical | membrane: moved (4: 6BO8: area 2905.01 -> 2891.84; …); 4 SVGs moved (6BO8, 6DQN, 6QUE, …) | identical | identical, not imported | identical | worse (2: tests.test_membrane::test_a_helix_at_a_face_lies_flat_from_35_degrees; …) |
| `membrane:JOIN_GAP` | named-entry:32KS,6UW4 | 2 -> 0 | **keep** | membrane, suite | identical | membrane: moved (2: 6DQN: area 3776.94 -> 3538.28; …); 2 SVGs moved (6DQN, 6QUE) | identical | identical, not imported | identical | worse (3: tests.test_membrane::test_comec_reads_as_its_figure_1c; …) |
| `na_runs:_JOINT_REACH` | named-entry:8VWQ | 5.0 -> 0.0 | **keep** | suite | identical; 2 SVGs moved (6QIK, 7F5S) | - | identical; 1 SVG moved (6QIK) | identical, not imported | identical | worse (2: tests.test_na_layout::test_a_joint_no_route_reaches_is_drawn_straight_a_little_further; …) |
| `na_runs:_BESIDE_TOLL` | named-entry:7Y7C | 50.0 -> 0.0 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (1: tests.test_na_layout::test_a_route_from_beside_a_face_keeps_off_the_strands_drawn_before_it) |
| `nucleic:HBOND_STRETCH` | named-entry:1HYS,1AOI | 4.0 -> 3.8 | **keep** | na, r100, suite | worse (7: 7F5S: PASS -> FAIL; …); 7 SVGs moved (3G6E, 6LAB, 6QIK, …) | - | worse (33: 6QIK: stem_recall 1.0 -> 0.994; …); 10 SVGs moved (1HYS, 1Z43, 5K7Z, …) | identical | identical | worse (6: tests.test_deep_check::test_f2_pairs_are_rnapolis_pairs[1AOI-outside0-agreement0]; …) |
| `panels:CLOSE` | named-entry:2LZM | 2.0 -> 0.0 | **keep** | domains, suite | identical | domains: worse (2: 6QIK: fallback counts_hidden 3 -> 4; …); 49 SVGs moved (10LI, 12DV, 1NID, …) | identical | identical, not imported | identical | worse (9: tests.test_perf::test_panel_fields_equal_the_domain_by_domain_rivals[0]; …) |
| `render:_SPOT_FAR` | named-entry:7QUP | 3 -> 1 | **keep** | r100, suite | worse (7: 3G6E: fallback mark_least_bad_spot 0 -> 1; …); 4 SVGs moved (3G6E, 4D5B, 6IZP, …) | - | identical; 1 SVG moved (6IZP) | identical, not imported | identical | worse (1: tests.test_text_clear::test_a_ligand_boxed_in_by_elements_and_text_still_gets_a_mark) |
| `render:_SPOT_MORE_RINGS` | named-entry:7LX9,9YLM,3MHT,1GOX | 4 -> 0 | **keep** | suite | identical; 41 SVGs moved (1GQ2, 1N6K, 1NID, …) | - | identical; 12 SVGs moved (1IHF, 1ZAA, 3KTQ, …) | identical, not imported | identical | worse (2: tests.test_hbonds::test_hbond_figures_are_byte_identical_to_those_before_the_overlay[1ZAA]; …) |
| `route:_SLAB_RUN` | named-entry:7SMR | 12.0 -> 3.0 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (1: tests.test_route::test_a_loop_leaves_the_slab_rather_than_run_along_it_over_an_element) |
| `route:_MAX_NODES_PAGE` | named-entry:7QUP | 8_000_000 -> 1000000000000 | **held** | - | identical | - | identical | identical, not imported | identical | identical |
| `route:_NA_TEXT_TOLL` | named-entry:6SO9,9CJI,7XSP | 25.0 -> 0.0 | **keep** | na, suite | identical; 1 SVG moved (6ZY7) | - | worse (2: 1IHF: na_protein_crossings 4 -> 5; …); 7 SVGs moved (1IHF, 30QM, 3CRX, …) | identical, not imported | identical | worse (1: tests.test_snapshots::test_figure_matches_snapshot[30QM-ncrna]) |
| `route:_REROUTE_TURNS_INK` | named-entry:7XSP | 12 -> 4 | **keep** | na, suite | identical | - | worse (1: 1IHF: na_protein_crossings 4 -> 5); 1 SVG moved (1IHF) | identical, not imported | identical | worse (1: tests.test_route::test_a_reroute_may_turn_more_to_free_text_it_ran_through) |
| `symmetry:MIN_HELIX_COPIES` | named-entry:1URN | 4 -> 3 | **keep** | suite | identical | - | identical | identical, not imported | identical | worse (2: tests.test_assembly::test_crystal_packing_of_copies_in_the_asymmetric_unit_is_not_reported_as_symmetry; …) |

### How each was turned off

- `layout:_NAME_OPEN` = 0.15: = _NAME_ROUTE: a domain name keeps no more than its routing keep-out from a loop's way out (the 7NZ2 0.1 off)
- `na_runs:_JOINT_NEAR` = 0.0: a near joint no route reaches is held to the far joints' test (na_runs._straight_joint)
- `align:SAME_FOLD_TM` = 0.0: any aligned pair is one fold by its TM-score (align.element_mates; H1)
- `align:SAME_FOLD_FRACTION` = 0.0: any aligned pair is one fold by its matched fraction (align.element_mates; H1)
- `io:PAIRED_MIN` = 1000000000: no two chains are ever held together by base pairs (io._paired never true)
- `io:_PAIR_DIST` = 0.0: no base atoms pair within 0 Å: io._paired never true
- `layout:_TAIL_COUNT_NEAR` = 1000000000.0: a tail count need not stand beside its letter (layout.tail_counts)
- `layout:_SEP_MIRROR_ROOM` = 0.0: mirror-image copies get no more room than other copies drawn alike (layout._separate: max(room, 0)); their extra push passes stay, as they are code
- `layout:_SEP_TIE_ROOM` = 0.0: the second chance for copies drawn alike is never kept (layout._tied_again)
- `layout:_STACK_SHUT` = 0.0: a stack gap needs no width to have a way out, so none is shut: no nudge, no move to the partner's other side, no pocket (layout._column_heights: they all ask way_out by _STACK_SHUT)
- `layout:_STACK_NUDGE` = 0.0: a strand in a shut stack gap is never nudged (layout._column_heights.opened)
- `layout:_POCKET_REACH` = 0.0: an accessory closing a pocket never moves (layout._column_heights.unpocketed)
- `layout:_END_ROOM` = 0.0: a stack gap a chain end faces gets no more room (layout._column_heights.room)
- `layout:_BA_MAX` = 1000000000: a βα row may hold any number of parallel strands (layout.can_be_row)
- `layout:_PART_BLOWUP` = 1000000000.0: parting domains never starts over (layout._part_domains)
- `membrane:_MAX_POLAR_CROSSED` = 0.15: = _MAX_POLAR: a slab many elements cross is held to the ordinary polarity
- `membrane:_MANY_CROSSING` = 1000000000: no slab is crossed by enough elements to relax its polarity limit
- `membrane:MIN_CROSSING_HELICES` = 1000000000: the normal is never refitted to the crossing helices (membrane._refit)
- `membrane:FACE_FLAT_TILT` = 45.0: = FLAT_TILT: a helix at a face lies flat from the ordinary tilt
- `membrane:JOIN_GAP` = 0: no two helix runs are joined in a membrane view (membrane._apart: 1 <= gap <= 0)
- `na_runs:_JOINT_REACH` = 0.0: no far joint is drawn straight (na_runs._straight_joint)
- `na_runs:_BESIDE_TOLL` = 0.0: a route from beside a face pays nothing for cells under strands drawn before
- `nucleic:HBOND_STRETCH` = 3.8: = HBOND_MAX: no stretched pair is taken in the second pass (nucleic pairs)
- `panels:CLOSE` = 0.0: no gap between a domain's padded elements is filled (panels.panel_fields)
- `render:_SPOT_FAR` = 1: a ligand mark is never looked for beyond _SPOT_RADIUS (render._place_mark), and _least_bad_spot's rings are narrowed to _SPOT_RADIUS too (not a faithful revert: that sweep reads the same constant)
- `render:_SPOT_MORE_RINGS` = 0: no ring past the first clear spot is searched for fewer crossed bodies
- `route:_SLAB_RUN` = 3.0: = _SLAB_COST: a step along the slab costs what a step across it does
- `route:_MAX_NODES_PAGE` = 1000000000000: a page's second route searches have no shared budget (never refused)
- `route:_NA_TEXT_TOLL` = 0.0: a loop routed again off nucleic-acid text pays nothing for that text's cells
- `route:_REROUTE_TURNS_INK` = 4: = _REROUTE_TURNS: a reroute off text's ink may add no more corners than off its clearance
- `symmetry:MIN_HELIX_COPIES` = 3: = MIN_HELIX_COPIES_ASKED: auto calls a screw of three copies helical

### What kept each

- `align:SAME_FOLD_TM`: suite worse (3: tests.test_align::test_auto_takes_the_sequence_only_when_its_pairs_superpose; tests.test_align::test_same_fold_needs_the_tm_score_and_the_elements; tests.test_homologues::test_the_comparison_counts_the_limits)
- `align:SAME_FOLD_FRACTION`: suite worse (1: tests.test_align::test_same_fold_needs_the_tm_score_and_the_elements)
- `io:PAIRED_MIN`: suite worse (4: tests.test_assembly::test_a_lattice_contact_of_two_base_pairs_is_no_pairing; tests.test_assembly::test_chains_that_pair_are_one_object_and_auto_keeps_both; tests.test_assembly::test_paired_chains_must_be_assembled_together)
- `io:_PAIR_DIST`: suite worse (4: tests.test_assembly::test_a_lattice_contact_of_two_base_pairs_is_no_pairing; tests.test_assembly::test_chains_that_pair_are_one_object_and_auto_keeps_both; tests.test_assembly::test_paired_chains_must_be_assembled_together)
- `layout:_TAIL_COUNT_NEAR`: suite worse (1: tests.test_text_clear::test_tail_count_sits_clear_beyond_its_letter[2.0])
- `layout:_SEP_MIRROR_ROOM`: r100 moved (1: 3W96: area 3042.22 -> 3115.51) · na worse (4: 3W96: spearman 0.892 -> 0.891; 3W96: procrustes 0.308 -> 0.309; 3W96: na_run_crossings 5 -> 6) · suite worse (2: tests.test_snapshots::test_figure_matches_snapshot[30QM-ncrna]; tests.test_symmetry::test_mirror_images_are_kept_when_they_cost_a_little_more_room)
- `layout:_SEP_TIE_ROOM`: r100 worse (5: 9FZL: PASS -> FAIL; 9FZL: A15 now fails; 9FZL: fallback counts_hidden 2 -> 3) · suite worse (12: tests.test_app_controls::test_two_fold_choice_is_offered_for_a_two_fold_assembly; tests.test_hbonds::test_interaction_lines_follow_turned_and_mirrored_copies[1UAM-mirror]; tests.test_heldout_copies::test_heldout_copies_are_drawn_alike[13FZ])
- `layout:_STACK_SHUT`: r100 worse (4: 7U9Z: PASS -> FAIL; 7U9Z: A3 now fails; 7U9Z: fallback counts_hidden 0 -> 1) · suite worse (3: tests.test_domains::test_domain_names_and_leaders_leave_every_loop_its_way[8BCU]; tests.test_layout::test_stacked_accessories_leave_every_gap_a_way_out[2FR2]; tests.test_layout::test_stacked_accessories_leave_every_gap_a_way_out[8BCU])
- `layout:_STACK_NUDGE`: r100 moved (1: 9NMN: area 10388.35 -> 9556.46)
- `layout:_POCKET_REACH`: suite worse (2: tests.test_domains::test_domain_names_and_leaders_leave_every_loop_its_way[8BCU]; tests.test_layout::test_stacked_accessories_leave_every_gap_a_way_out[8BCU])
- `layout:_END_ROOM`: r100 worse (7: 6DQN: PASS -> FAIL; 6DQN: A4 now fails; 6DQN: fallback loop_over_text 0 -> 1) · suite worse (5: tests.test_figcheck::test_accessories_named_not_drawn_level_are_reported[2OMF-A-3-3]; tests.test_provenance::test_a_pairing_named_but_not_level_is_counted; tests.test_provenance::test_the_command_line_prints_one_block)
- `layout:_BA_MAX`: suite worse (2: tests.test_template::test_a_row_holds_at_most_ten_strands[10-True]; tests.test_template::test_a_row_holds_at_most_ten_strands[11-False])
- `layout:_PART_BLOWUP`: domains moved (1: 8U1N: area 1331.46 -> 8059.7) · suite worse (5: tests.test_copy_labels::test_no_chain_draws_one_label_twice[8UTF-auto-none-auto]; tests.test_domains::test_domains_sharing_a_sheet_of_copies_are_not_parted; tests.test_domains::test_parting_domains_sums_no_overlap_unless_it_may_stop)
- `membrane:_MAX_POLAR_CROSSED`: suite worse (3: tests.test_membrane::test_a_slab_crossed_by_many_long_elements_may_be_more_polar; tests.test_membrane::test_a_slab_many_long_elements_cross_but_not_oily_enough_is_refused; tests.test_membrane::test_the_membrane_is_found_on_every_chain_whatever_chains_draws)
- `membrane:_MANY_CROSSING`: suite worse (3: tests.test_membrane::test_a_slab_crossed_by_many_long_elements_may_be_more_polar; tests.test_membrane::test_a_slab_many_long_elements_cross_but_not_oily_enough_is_refused; tests.test_membrane::test_the_membrane_is_found_on_every_chain_whatever_chains_draws)
- `membrane:MIN_CROSSING_HELICES`: membrane moved (7: 6BO8: area 2905.01 -> 3003.33; 6DQN: area 3776.94 -> 3393.39; 6QUE: area 3673.36 -> 4215.02) · suite worse (7: tests.test_membrane::test_a_soluble_helix_bundle_beside_the_membrane_does_not_tilt_its_normal; tests.test_membrane::test_fewer_than_three_crossing_helices_keep_the_first_estimate; tests.test_membrane::test_helical_membranes_are_unchanged_by_the_barrel_estimate[1BL8.cif])
- `membrane:FACE_FLAT_TILT`: membrane moved (4: 6BO8: area 2905.01 -> 2891.84; 6DQN: area 3776.94 -> 3444.04; 6QUE: area 3673.36 -> 3819.28) · suite worse (2: tests.test_membrane::test_a_helix_at_a_face_lies_flat_from_35_degrees; tests.test_membrane::test_comec_reads_as_its_figure_1c)
- `membrane:JOIN_GAP`: membrane moved (2: 6DQN: area 3776.94 -> 3538.28; 6QUE: area 3673.36 -> 4118.51) · suite worse (3: tests.test_membrane::test_comec_reads_as_its_figure_1c; tests.test_membrane::test_copies_broken_where_the_joined_helix_crosses_are_joined_not_cut; tests.test_membrane::test_flat_helix_pieces_a_break_apart_are_joined_only_going_on_as_one)
- `na_runs:_JOINT_REACH`: suite worse (2: tests.test_na_layout::test_a_joint_no_route_reaches_is_drawn_straight_a_little_further; tests.test_na_layout::test_joined_pieces_kept_end_to_end_may_join_straight_farther_apart)
- `na_runs:_BESIDE_TOLL`: suite worse (1: tests.test_na_layout::test_a_route_from_beside_a_face_keeps_off_the_strands_drawn_before_it)
- `nucleic:HBOND_STRETCH`: r100 worse (7: 7F5S: PASS -> FAIL; 7F5S: A17 now fails; 7F5S: A6 now fails) · na worse (33: 6QIK: stem_recall 1.0 -> 0.994; 6QIK: spearman 0.752 -> 0.735; 6QIK: procrustes 0.217 -> 0.286) · suite worse (6: tests.test_deep_check::test_f2_pairs_are_rnapolis_pairs[1AOI-outside0-agreement0]; tests.test_independent_kernels::test_the_findings_are_exactly_these; tests.test_independent_kernels::test_the_record_names_these_findings)
- `panels:CLOSE`: domains worse (2: 6QIK: fallback counts_hidden 3 -> 4; 8EM7: fallback counts_hidden 0 -> 1) · suite worse (9: tests.test_perf::test_panel_fields_equal_the_domain_by_domain_rivals[0]; tests.test_perf::test_panel_fields_equal_the_domain_by_domain_rivals[1]; tests.test_perf::test_panel_fields_equal_the_domain_by_domain_rivals[2])
- `render:_SPOT_FAR`: r100 worse (7: 3G6E: fallback mark_least_bad_spot 0 -> 1; 4D5B: PASS -> FAIL; 4D5B: A6 now fails) · suite worse (1: tests.test_text_clear::test_a_ligand_boxed_in_by_elements_and_text_still_gets_a_mark)
- `render:_SPOT_MORE_RINGS`: suite worse (2: tests.test_hbonds::test_hbond_figures_are_byte_identical_to_those_before_the_overlay[1ZAA]; tests.test_snapshots::test_figure_matches_snapshot[1C3W-membrane])
- `route:_SLAB_RUN`: suite worse (1: tests.test_route::test_a_loop_leaves_the_slab_rather_than_run_along_it_over_an_element)
- `route:_NA_TEXT_TOLL`: na worse (2: 1IHF: na_protein_crossings 4 -> 5; 30QM: na_protein_crossings 15 -> 17) · suite worse (1: tests.test_snapshots::test_figure_matches_snapshot[30QM-ncrna])
- `route:_REROUTE_TURNS_INK`: na worse (1: 1IHF: na_protein_crossings 4 -> 5) · suite worse (1: tests.test_route::test_a_reroute_may_turn_more_to_free_text_it_ran_through)
- `symmetry:MIN_HELIX_COPIES`: suite worse (2: tests.test_assembly::test_crystal_packing_of_copies_in_the_asymmetric_unit_is_not_reported_as_symmetry; tests.test_symmetry_scenarios::test_three_copies_on_a_screw_are_helical_only_when_asked)

## Measured with another's value

- `route:_Router.second:fallback:route_second_search_refused`: with `route:_MAX_NODES_PAGE`, **held**

## Not tried: needs code removal

- `layout:_place_domain_names:fallback:domain_name_least_bad` (fallback, named-entry:7NZ2): layout._place_domain_names: the name with no clear place must go somewhere or be left out; either is code
- `layout:_place_domain_names:fallback-marker` (fallback-marker, named-entry:7NZ2): layout._place_domain_names: the marker of domain_name_least_bad (the same code)
- `na_ends:_OFF_TURNS` (constant, named-entry:13FZ): na_ends.off_unheld: a tuple; the hook sets an int or a float only (scripts/perturb.py, ruling Q29): turning the nudge off removes its search
- `na_ends:off_unheld:fallback-marker` (fallback-marker, named-entry:13FZ): na_ends.off_unheld, called by layout.build_layout: turning the 13FZ nudge off removes that call (its search constants are tuples)
- `na_labels:_GAP_OUT` (constant, named-entry:7Y7C): na_labels._label_gaps: a tuple; the hook sets an int or a float only (scripts/perturb.py, ruling Q29): dropping 4 and 6 (7Y7C) edits the tuple
- `na_labels:_GAP_TURNS` (constant, named-entry:8ZTU): na_labels._gap_ring: a tuple; the hook sets an int or a float only (scripts/perturb.py, ruling Q29): turning the ring off removes _gap_ring's call
- `na_runs:_JOINT_BENDS` (constant, named-entry:8VWQ): na_runs._joint_bend: a tuple; the hook sets an int or a float only (scripts/perturb.py, ruling Q29): turning the bends off removes _joint_bend
- `na_runs:_JOINT_SWINGS` (constant, named-entry:8VWQ,7Y7C): na_runs._joint_bend: a tuple; the hook sets an int or a float only (scripts/perturb.py, ruling Q29): turning the swings off edits the tuple
- `render:_draw:fallback:mark_least_bad_spot` (fallback, named-entry:7QUP): render._place_mark and _least_bad_spot: the mark with no clear spot went on its anchor before; that is code
- `route:_MAX_NODES_WINDOWS` (constant, named-entry:7NZ2,9FZL): route._Router.__init__ and route_loops: the special case is the windows' own pool (drawing from the page's budget gave 7NZ2 3 fallback curves); no value of it draws them from the page's budget again, which is code
- `route:_NAME_ALONG` (constant, named-entry:6MHS,7SMR,8UUP,1KX5): route._clear_helix_names: a tuple; the hook sets an int or a float only (scripts/perturb.py, ruling Q29): no value turns the moves along a helix off
- `route:_clear_helix_names:fallback:helix_name_moved` (fallback, named-entry:6MHS,7SMR,8UUP,1KX5): route._clear_helix_names, called by route.route_loops: turning the moves off removes that call
- `route:<module>:fallback-marker` (fallback-marker, named-entry:7QUP): route: route_unsnapped (_TOUCH), route_second_search_refused (measured by route:_MAX_NODES_PAGE) and the windows' pool: code paths
