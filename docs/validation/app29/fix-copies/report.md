# Plan 29 Task 8b: copies modelled in part are one protein to focus auto (26EL, ruling P29-14)

Written by scripts/app29_light.py report at 1a15466 (clean); ruling P29-13's light loop: the tune half, the bench and the measures ruling P29-14 names; the check half left to Task N+1.

**Verdict: kept.** Expected to change: 26EL's figure (A alone, '1 of 20 subunits shown'); elsewhere only symmetric entries whose protomer holds a copy modelled in part.

## Mechanism

cli._one_of_each, which picks focus auto's chains (one of each distinct protein in a protomer), grouped chains by the whole modelled sequence (difflib ratio > 0.9). A copy modelled in part falls under it: 26EL's protomer holds A (1-682), F (1-291) and M (292-682), ratios 0.59 and 0.73, so A, F and M were drawn, 132-301 Å apart on a diagonal, said as '1 of 7 subunits'. Two chains are now also one protein when model.homolog_chains (the copy rule the labels already share by: at least 20 residues at the same numbers, 90 % alike, the numbering offset allowed) gives them one root. No new constant.

## The tune half (scripts/app29.py run, the app's way)

- before: `docs/validation/app29/fix-ghost/run` (git e8f1e47); after: `docs/validation/app29/fix-copies/run` (git 66c1a3d)
- entries 47 -> 47; missing: none
- the app's SVG changed on: 26EL; its PNG on: 26EL
- criteria, floor, error or refusal changed on: none; worse on: none
- draw seconds, all tune entries: 109.74 -> 108.5 (8 entries drawn at once: indicative)
- 26EL: draw 1.41 -> 1.06 s; the entry's process (draw, PNG, checks) 1.6 -> 1.24 s
- **26EL**: the named entry: focus auto draws A alone (58 elements) instead of A, F and M on a diagonal; the note reads '1 of 20 subunits shown · C2 symmetry' (was '1 of 7'); the page 60.0 x 75.3 -> 27.9 x 32.9 units (diagnosis-t8.md); no criterion moved

## bench --check --exact

- 22 of 22 identical, exit 0

## The fix's own measure (A-third.json)

### the third sample (random100-third, default view)

- before git cdd315f, after git 66c1a3d; entries 100 -> 100; pass 82 -> 82; keys: moved
- worse: none
- better: none
- moved (no judge orders it): 7SPI: area 426.01 -> 338.03
- the SVG moved on: 7SPI
- A17 measures (100 entries): changed on none; worse: none; better: none
- failing A17 before: 3EGZ, 6QIK, 6ZY7, 8VBG, 8XJV; after: 3EGZ, 6QIK, 6ZY7, 8VBG, 8XJV
- **7SPI**: a C13 type IV secretion core whose protomer holds two copies of TraV, A1 (150-204) and B1 (91-203), each modelled in part (ratio 0.90, not above it): both were drawn in full, with the same labels B-E twice. They are one protein by their numbers (54 residues at the same numbers, all alike), so A1 is drawn once and B1 only where it completes A1's sheet (D' and E', grey); the note says 1 of 26 subunits (78 chains over A1, C1, E1; C1/D1 and E1/F1 were already one protein each, by the ratio), was 1 of 20. Area 426.01 -> 338.03 (no judge orders it); pass, criteria and fallbacks identical


## The fix's own measure (A-na.json)

### the NA cases, controls and campaign entries (default view)

- before git cdd315f, after git 66c1a3d; entries 69 -> 69; pass 54 -> 54; keys: identical
- worse: none
- better: none
- moved (no judge orders it): none
- the SVG moved on: none
- A17 measures (69 entries): changed on none; worse: none; better: none
- failing A17 before: 1YTB, 30QM, 3CRX, 3EGZ, 3KTQ, 4PJO, 4YBB, 6QIK, 6ZY7, 8BW5, 8VBG, 8XJV, 9K4F, 9YLM; after: 1YTB, 30QM, 3CRX, 3EGZ, 3KTQ, 4PJO, 4YBB, 6QIK, 6ZY7, 8BW5, 8VBG, 8XJV, 9K4F, 9YLM


## The fix's own measure (A-navalidate.json)

### the NA campaign (scripts/na_validate.py, docs/validation/na)

- before git cdd315f, after git 66c1a3d; entries 39 -> 39; pass 34 -> 34
- measures changed on: none; worse: none; better: none


## The fix's own measure (A-third-committed.json)

### the third sample: its committed report (7696947, views all; read-only) against the fix

For reference, not judged: the before is older code, so its differences are not this fix's.

- before git ?, after git 66c1a3d; entries 100 -> 100; pass 83 -> 82; keys: worse
- worse: 2YEJ: PASS -> FAIL; 2YEJ: A19 now fails; 3FT8: PASS -> FAIL; 3FT8: A19 now fails; 3G6E: fallback glyph_width_guessed 0 -> 52; 3P2W: PASS -> FAIL; 3P2W: A19 now fails; 4AS1: fallback glyph_width_guessed 0 -> 4; 4B2T: fallback counts_hidden 0 -> 1; 5U3R: PASS -> FAIL; 5U3R: A19 now fails; 5XMX: PASS -> FAIL; 5XMX: A19 now fails; 6LTO: PASS -> FAIL; 6LTO: A19 now fails; 6M67: PASS -> FAIL; 6M67: A19 now fails; 6QIK: fallback counts_hidden 0 -> 1; 6QIK: fallback glyph_width_guessed 0 -> 52; 6QIK: fallback na_run_unclear 0 -> 1; 6ZTQ: fallback counts_hidden 0 -> 1; 7F5S: fallback counts_hidden 0 -> 5; 7F5S: fallback glyph_width_guessed 0 -> 68; 8CRC: PASS -> FAIL; 8CRC: A19 now fails; 8XJV: fallback counts_hidden 0 -> 5; 8XJV: fallback glyph_width_guessed 0 -> 76; 8XJV: fallback helix_name_moved 0 -> 19; 8XJV: fallback na_gap_unclear 0 -> 3; 8XJV: fallback na_tail_straight 0 -> 4; 8XJV: fallback na_trace_unclear 0 -> 3; 8XJV: fallback page_text_small 0 -> 1; 9FZL: fallback counts_hidden 0 -> 1; 9NMN: fallback loop_over_text 0 -> 1
- better: 4AS1: FAIL -> PASS; 4AS1: A17 now passes; 4YJ5: FAIL -> PASS; 4YJ5: A15 now passes; 6DQN: FAIL -> PASS; 6DQN: A18 now passes; 6VTX: FAIL -> PASS; 6VTX: A17 now passes; 7F5S: FAIL -> PASS; 7F5S: A3 now passes; 7F5S: A6 now passes; 8PIN: FAIL -> PASS; 8PIN: A15 now passes; 9FZL: FAIL -> PASS; 9FZL: A16 now passes; 9FZL: A3 now passes; 9FZL: A4 now passes; 9FZL: A6 now passes
- moved (no judge orders it): 1RKA: area 1121.25 -> 679.88; 1VSK: area 352.68 -> 361.36; 3EGZ: area 1478.07 -> 1564.02; 3FT8: area 508.73 -> 433.34; 6PO4: area 524.5 -> 643.66; 6QIK: area 27610.17 -> 27350.66; 7SPI: area 426.01 -> 338.03; 7Y63: area 1783.14 -> 1764.17; 8EM7: area 8203.86 -> 8239.72; 8U1N: area 1723.55 -> 1117.64; 9DU3: area 552.51 -> 654.76; 9FZL: area 65458.72 -> 14561.27
- the SVG moved on: 1NID, 1RKA, 1VSK, 3ABR, 3EGZ, 3FT8, 3G6E, 3WAU, 4AS1, 4GWY, 5S4U, 6IZP, 6PO4, 6QIK, 6VBG, 7F5S, 7SPI, 7XBK, 7Y63, 8EM7, 8U1N, 8XJV, 9DU3, 9FZL


## The control

No near-identical paralogs sharing residue numbering exist among the 710 structures at hand (tests/data, the tune half, the zl download cache less the check half): every pair of chains of two polymer entities alike at 70 % or more of 20 or more shared numbers is 96.7 % alike or more, and homolog_chains joins each (so no real entry shows it keeping such paralogs apart). Three of those pairs the fix newly joins (whole-sequence ratio 0.9 or less): HIV-1 reverse transcriptase's p66 and p51 in 8VBG and 1HYS (one gene, p51 its first 440 residues, folded otherwise) and 6V6S's poly-UNK chains of different entities (unassigned poly-alanine models, their residues all UNK). None is symmetric, so focus auto never reaches _one_of_each on them, and their figures do not change (8VBG and 1HYS are NA controls: identical). Pinned as it is: tests/test_focus.py::test_one_protein_by_its_numbers_is_today_one_protein_even_as_two_entities (8VBG captured); the rule's edge on a synthetic pair (85 % alike: two proteins): test_near_identical_paralogs_below_the_copy_rule_stay_two_proteins. The risk, named: a symmetric assembly of p66/p51 heterodimers, or of poly-UNK models of different proteins, would draw one chain where it drew two. Reach: of the 13 structures at hand where the fix joins a new pair, 4 are symmetric and reach it: 26EL and 7SPI (both measured above), 7SPB (7SPI's sibling: the same TraV pair) and 9GTP (C6; its '1l' and '2l' are copies of one baseplate protein: 8 chains -> 7, '1 of 8' -> '1 of 9'); 9TX7 (tests/data, C3) changes _one_of_each but its focus auto stays whole (6 chains).

## Notes

The third sample's committed report (7696947, views all) against this fix's run: 83 -> 82 pass, 38 entries differ; 37 of them differ identically against the BASE run (cdd315f), so they come from the code since 7696947, not from this fix; the 38th is 7SPI (area only). Read-only: measures/A-third-committed.json.
