# Plan 29 Task 8b: a base letter sits nearest its own strand (22MZ, ruling P29-14)

Written by scripts/app29_light.py report at 8919aa3 (clean); ruling P29-13's light loop: the tune half, the bench and the measures ruling P29-14 names; the check half left to Task N+1.

**Verdict: rejected.** Expected to change: 22MZ's letters (7 nearer another strand -> 3); elsewhere only letters whose first clear spot was nearer another strand.

## Mechanism

na_labels._place_letters took a single-stranded nucleotide's first spot (beside its strand, further out, or a little along it) clear of text and lines, whichever strand it was nearest. 22MZ's pseudoknot routes runs 2 and 6 side by side 0.65 NT_ROW apart: U57's spot beside run 6 lay 0.35 NT_ROW from run 2's line, on run 2's far side, and 7 of its 26 unpaired letters were nearer another strand than their own. A spot beside or along the strand is now clear only where no other strand's drawn line is nearer it than its own strand's (na_labels._own_nearest; a tie is its own). The later passes (by the drawn outlines, the rings) are unchanged. No new constant.

## The tune half (scripts/app29.py run, the app's way)

- before: `docs/validation/app29/fix-copies/run` (git 66c1a3d); after: `docs/validation/app29/fix-letters/run` (git 87952c1)
- entries 47 -> 47; missing: none
- the app's SVG changed on: 22MZ, 9YPK, 9YPN; its PNG on: 22MZ, 9YPK, 9YPN
- criteria, floor, error or refusal changed on: none; worse on: none
- draw seconds, all tune entries: 108.5 -> 120.12 (8 entries drawn at once: indicative)
- 22MZ: draw 0.39 -> 0.4 s; the entry's process (draw, PNG, checks) 0.41 -> 0.42 s
- **22MZ**: the named entry: letters nearer another strand 7 -> 3 (0:50 G, 0:56 C, 0:67 G stay: no spot beside or along the strand both clear and nearest their own); no criterion moved
- **9YPK**: letters nearer another strand 408 -> 130 of 1196 (262 moved; 2 newly nearer another strand); moved letters' median distance to their nucleotide 1.0 -> 1.24 NT_ROW; no criterion moved
- **9YPN**: letters nearer another strand 347 -> 103 of 1144 (219 moved; 2 newly nearer another strand); no criterion moved

## bench --check --exact

- 22 of 22 identical, exit 0

## The fix's own measure (B-na.json)

### the NA cases, controls and campaign entries (default view)

- before git 66c1a3d, after git 87952c1; entries 69 -> 69; pass 54 -> 52; keys: worse
- worse: 2GDI: PASS -> FAIL; 2GDI: A6 now fails; 9WB3: PASS -> FAIL; 9WB3: A6 now fails
- better: none
- moved (no judge orders it): 6IZP: area 85.29 -> 77.33; 9IWF: area 232.25 -> 229.37
- the SVG moved on: 1EHZ, 1GID, 1HQ1, 1MSW, 1QTQ, 1RLG, 1URN, 1Y26, 1Z43, 2GDI, 2GIS, 2NZ4, 2OEU, 30QM, 3D0U, 3EGZ, 3IGI, 4AS1, 4LVV, 4OO8, 4P9R, 4PJO, 4YBB, 5F9R, 6IZP, 6QIK, 7XSP, 8XJV, 9IWF, 9K4F, 9OWM, 9WB3
- A17 measures (69 entries): changed on 1URN, 5F9R; worse: 1URN ['na_protein_crossings']; better: 5F9R ['na_protein_crossings']
- failing A17 before: 1YTB, 30QM, 3CRX, 3EGZ, 3KTQ, 4PJO, 4YBB, 6QIK, 6ZY7, 8BW5, 8VBG, 8XJV, 9K4F, 9YLM; after: 1YTB, 30QM, 3CRX, 3EGZ, 3KTQ, 4PJO, 4YBB, 6QIK, 6ZY7, 8BW5, 8VBG, 8XJV, 9K4F, 9YLM
  - 1URN: na_protein_crossings 1 -> 2
  - 5F9R: na_protein_crossings 7 -> 6
- **1EHZ**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1GID**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1HQ1**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1MSW**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1QTQ**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1RLG**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1URN**: WORSE: na_protein_crossings 1 -> 2 (a protein loop routed across a moved letter's strand)
- **1Y26**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **1Z43**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **2GDI**: WORSE: A6 now fails, a letter on a tie (na:nt:X:45 x na:tie:20-45): the letter moved to a spot nearest its own strand, where the tie's label is placed after it
- **2GIS**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **2NZ4**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **2OEU**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **30QM**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **3D0U**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **3EGZ**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **3IGI**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **4AS1**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **4LVV**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **4OO8**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **4P9R**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **4PJO**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **4YBB**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **5F9R**: better: na_protein_crossings 7 -> 6
- **6IZP**: area 85.29 -> 77.33 (moved); letters moved
- **6QIK**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **7XSP**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **8XJV**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **9IWF**: area 232.25 -> 229.37 (moved); letters moved
- **9K4F**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **9OWM**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **9WB3**: WORSE: A6 now fails, a letter on a tie mark (na:nt:R:75 x na:tie:43-75:mark), as 2GDI


## The fix's own measure (B-navalidate.json)

### the NA campaign (scripts/na_validate.py, docs/validation/na)

- before git 66c1a3d, after git 87952c1; entries 39 -> 39; pass 34 -> 34
- measures changed on: 1URN, 5F9R; worse: 1URN ['na_protein_crossings']; better: 5F9R ['na_protein_crossings']
  - 1URN: na_protein_crossings 1 -> 2
    WORSE: na_protein_crossings 1 -> 2 (a protein loop routed across a moved letter's strand)
  - 5F9R: na_protein_crossings 7 -> 6
    better: na_protein_crossings 7 -> 6


## The fix's own measure (B-third.json)

### the third sample (random100-third, default view)

- before git 66c1a3d, after git 87952c1; entries 100 -> 100; pass 82 -> 83; keys: moved
- worse: none
- better: 3G6E: FAIL -> PASS; 3G6E: A6 now passes
- moved (no judge orders it): 6IZP: area 85.29 -> 77.33
- the SVG moved on: 3EGZ, 3G6E, 4AS1, 6IZP, 6QIK, 7F5S, 8XJV, 9OWM
- A17 measures (100 entries): changed on none; worse: none; better: none
- failing A17 before: 3EGZ, 6QIK, 6ZY7, 8VBG, 8XJV; after: 3EGZ, 6QIK, 6ZY7, 8VBG, 8XJV
- **3EGZ**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **3G6E**: better: A6 now passes (FAIL -> PASS)
- **4AS1**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **6IZP**: area 85.29 -> 77.33 (moved); letters moved
- **6QIK**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **7F5S**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **8XJV**: base letters moved off spots nearer another strand (the SVG only); no key moved
- **9OWM**: base letters moved off spots nearer another strand (the SVG only); no key moved


## Notes

Rejected and reverted (43e7637); the code is kept as reverted.patch beside this record. Measured at 87952c1 against Fix A's runs (749648c's code), on zyla-office: the NA set (na_fixes_record's cases and controls plus the NA campaign's entries, 69, default view) 54 -> 52 pass, 2GDI and 9WB3 failing A6, 1URN one more NA-protein crossing; so the keys are not identical or better (ruling P29-14). Also against it: tests/test_snapshots.py's 1EHZ-trna (na:nt:A:15 moves about 40 px) and 30QM-ncrna (7 nucleotides move) fail, and the moves do not read clearly better; _place_letters on the whole 4YBB SSU RNA takes 1.30 -> 2.40 s on one host (gpu-desktop, 6 runs each), over test_na_layout's 1.5 s limit there (_own_nearest measures every strand segment for every spot; a cell hash would bound it, but not the A6 failures). The outline-pass variant (3 -> 0 on 22MZ, letters 3.8 rows out) was not taken. 22MZ's letters stay as they were (7 of 26 nearer another strand), for the author.
