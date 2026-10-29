# The readability floor

Plan 28 Task 4 (the review's 3.3, rulings Q39-Q41). Written by: foldmap 0.2.0, git e1392f79daf7 (asserted by runner, dirty False), code fingerprint 7879eea05bd9301c, scripts/readability_floor.py sha256 7d8f76b7e7f5; Python 3.13.13 on Linux-6.8.0-138-generic-x86_64-with-glibc2.35; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

A figure below the readability floor says so (floor.readability_floor): in `foldmap plot`'s notes (and with `--strict-floor` it writes no figure and exits 1), in `foldmap check` as R1 (readability; a fail exits 1) and its --json report, in the MCP result and under the app's figure, naming the existing options that draw less. The figure itself is unchanged. This record is what the reasons and their limits rest on.

## The reasons and their limits

- **text_small**: the figure's own text (element names, chain ends) under 6 pt on its page (layout.MIN_TEXT_PT, the random-100 C4 floor): the page grows to keep it at 6 pt until layout.MAX_PAGE_IN (48 in), and past that the whole figure is scaled down. Measured here at the default page size, as `foldmap plot` writes it. The smallest text drawn (residue counts, base letters, ligand names, 0.6-0.9 of the figure's text by design) is recorded but not judged: it is under 6 pt on 17 of the 18 default figures measured, those far above the floor included.
- **route_budget**: the page's budget for second, wider route searches (route._MAX_NODES_PAGE) spent, counted as the named fallbacks route_second_search_refused. It fired on none of the figures measured here and on none of the samples' default views (census below: their reports count none refused), so it is tested on synthetic measures only here. It does fire past the default views: 7NZ2's domains view (the held-out sample) spends the budget, 6 searches refused (docs/validation/prune/report.md, held), so `foldmap check 7NZ2 --domains auto` fails R1 on it and R1's exit-code change reaches past the default views. A stuck loop of another kind is A3's, as it happens on small pages too.
- **time**: reported, never judged (ruling Q40). A figure that took longer than the hosted app's 90 s (criteria.LIMITS app_s) to make, measured as spec.draw_spec's wall time, is said so in the notes, the check's text, the MCP note and the app; it never fails R1, never changes `foldmap check`'s exit code or --strict-floor, and is not in the check's --json report (which holds no wall time). It depends on the machine and its load: the figures here were drawn 64 at a time. Time notes here: 9FZL (116 s), 8XJV (208 s), 8VWQ (91 s).

## Census of the three samples' default views

From each sample's committed campaign report, at the commit it names (not this code): the default views whose notes say the page held its text under 6 pt, those over 60 s (A2) and over 90 s (A2b), and those with a second route search refused for the page's budget.

| sample | commit | default views | text under 6 pt | over 60 s | over 90 s | route budget spent |
|---|---|---|---|---|---|---|
| random100 | 76969474ca5b | 100 | 7JU4 (5.5), 6JLN (4.5) | 7Y7C (66.92), 8VWQ (64.51) | none | none |
| random100-heldout | 76969474ca5b | 100 | 6ET5 (2.4), 8WMV (5.5), 7NZ2 (5.2), 8ZTU (5.7) | 8ZTU (60.74) | none | none |
| random100-third | 76969474ca5b | 100 | 9FZL (2.0), 8XJV (5.4) | 8XJV (129.02) | 8XJV (129.02) | none |

## What each option achieves on the figures below the floor

Below the floor at this code: 6 of the 18 figures measured, 8XJV (text_small), 7JU4 (text_small), 6ET5 (text_small), 8WMV (text_small), 7NZ2 (text_small), 8ZTU (text_small).

| option | tried on | brought above | median seconds | median text pt | entries |
|---|---|---|---|---|---|
| `--focus CHAIN` | 6 | 6 | 4.4 | 6.00 | 8XJV, 7JU4, 6ET5, 8WMV, 7NZ2, 8ZTU |
| `--chains CHAIN` | 6 | 6 | 2.8 | 6.00 | 8XJV, 7JU4, 6ET5, 8WMV, 7NZ2, 8ZTU |
| `the preset measured (not an option)` | 6 | 0 | 59.3 | 5.48 | 8XJV, 7JU4, 6ET5, 8WMV, 7NZ2, 8ZTU |

The preset measured is `--set nucleotide_labels=false --set ligands=none --set glycans=false --set helices_310=false` (existing style keys only), not an option: ruling Q40 dropped the --reduced-detail option this measured, as it brought none of the figures below the floor above it. It leaves out low-priority annotations, and so most of the text drawn under 6 pt, but the page's size is set by its elements, so it does not raise the figure's text; on 7Y7C (6.00 pt to 5.37 pt) it took a figure above the floor to below it.

## What `foldmap check` now says

R1 fails, and `foldmap check` exits 1, on the figures below the floor. Before R1 these default figures' other checks (the S-checks and A1-A20 judged here, criteria.check_figure) passed on 4 of them, which now exit 1 on R1 alone: 7JU4, 6ET5, 8WMV, 7NZ2; and already failed on 2, which exited 1 before: 8XJV (A3, A6), 8ZTU (A6). The census above is of the samples' reports at their own commit, so that every sample default whose text is under 6 pt is among those measured here is inferred from it, not run at this code on all 300.

## Every figure

Page in inches; the figure's text and its smallest text in points at that page size; texts drawn under 6 pt of all texts; seconds to make (draw_spec); named fallbacks taken; A5 element overlaps / A6 text collisions; the floor's reasons. Offered: the options the notes name for that figure.

### 9FZL

Notes: it took 116 s to make, over the hosted app's 90 s (said, not judged)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 45.5 x 20.9 | 6.00 | 3.60 | 179 / 750 | 116.0 | counts_hidden 1 | 0 / 0 | above |
| the preset measured (not an option) | 44.7 x 21.5 | 6.00 | 4.32 | 52 / 620 | 108.2 | counts_hidden 1 | 0 / 0 | above |

### 8XJV

Notes: below the readability floor: its text is 5.4 pt, under 6 pt, on its 47 x 48 in page (its smallest text 3.9 pt); drawing less may bring it above: --focus A (chain A in full, with what its neighbours add to its fold); --chains A (chain A only); it took 208 s to make, over the hosted app's 90 s (said, not judged)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 46.9 x 48.0 | 5.48 | 3.95 | 1288 / 1288 | 207.6 | counts_hidden 5, glyph_width_guessed 76, helix_name_moved 19, na_gap_unclear 3, na_tail_straight 4, na_trace_unclear 3, page_text_small 1 | 0 / 3 | text_small |
| --focus A (offered) | 9.0 x 10.2 | 6.00 | 4.32 | 22 / 42 | 5.8 | glyph_width_guessed 20 | 0 / 0 | above |
| --chains A (offered) | 8.1 x 10.5 | 6.00 | 4.32 | 21 / 33 | 3.7 | glyph_width_guessed 20 | 0 / 0 | above |
| the preset measured (not an option) | 45.6 x 48.0 | 5.19 | 3.74 | 1257 / 1257 | 210.1 | counts_hidden 13, glyph_width_guessed 57, helix_name_moved 15, na_gap_unclear 4, na_tail_straight 3, na_trace_unclear 3, page_text_small 1 | 0 / 5 | text_small |

### 7QUP (29,575 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 7.8 x 10.0 | 6.00 | 3.60 | 12 / 73 | 7.1 | none | 0 / 0 | above |
| the preset measured (not an option) | 7.7 x 10.7 | 6.00 | 5.10 | 1 / 62 | 7.4 | none | 0 / 0 | above |

### 6YF9 (25,200 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 5.2 x 4.3 | 7.37 | 6.26 | 0 / 15 | 8.6 | none | 0 / 0 | above |
| the preset measured (not an option) | 5.2 x 4.3 | 7.37 | 6.26 | 0 / 15 | 8.9 | none | 0 / 0 | above |

### 7JU4 (23,990 residues deposited)

Notes: below the readability floor: its text is 5.5 pt, under 6 pt, on its 33 x 48 in page (its smallest text 3.3 pt); drawing less may bring it above: --focus 0 (chain 0 in full, with what its neighbours add to its fold); --chains 0 (chain 0 only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 32.8 x 48.0 | 5.54 | 3.32 | 1298 / 1298 | 65.6 | page_text_small 1 | 0 / 0 | text_small |
| --focus 0 (offered) | 2.5 x 12.2 | 7.37 | 5.30 | 1 / 9 | 3.7 | none | 0 / 0 | above |
| --chains 0 (offered) | 1.8 x 11.0 | 7.37 | 5.30 | 1 / 5 | 1.2 | none | 0 / 0 | above |
| the preset measured (not an option) | 34.3 x 48.0 | 5.87 | 4.22 | 1162 / 1162 | 50.8 | page_text_small 1 | 0 / 0 | text_small |

### 6JLN

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 26.4 x 11.3 | 6.00 | 3.60 | 116 / 300 | 47.9 | none | 0 / 0 | above |
| the preset measured (not an option) | 25.0 x 12.9 | 6.00 | 4.32 | 13 / 197 | 19.8 | none | 0 / 0 | above |

### 7Y7C

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 47.1 x 41.8 | 6.00 | 3.72 | 654 / 1201 | 88.6 | counts_hidden 4, glyph_width_guessed 100 | 0 / 0 | above |
| the preset measured (not an option) | 48.0 x 38.9 | 5.37 | 3.87 | 1155 / 1155 | 96.3 | counts_hidden 4, glyph_width_guessed 75, na_run_unclear 2, page_text_small 1 | 0 / 1 | text_small |

### 8VWQ

Notes: it took 91 s to make, over the hosted app's 90 s (said, not judged)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 42.3 | 6.00 | 3.60 | 657 / 1243 | 90.7 | counts_hidden 3, glyph_width_guessed 76 | 0 / 0 | above |
| the preset measured (not an option) | 46.4 x 41.8 | 6.00 | 4.32 | 603 / 1189 | 89.9 | counts_hidden 2, glyph_width_guessed 57, helix_name_moved 1 | 0 / 2 | above |

### 8D3D (23,504 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 15.2 x 8.1 | 6.00 | 3.72 | 20 / 128 | 38.3 | none | 0 / 0 | above |
| the preset measured (not an option) | 15.7 x 8.9 | 6.00 | 4.32 | 5 / 113 | 34.4 | none | 0 / 0 | above |

### 7SPB (23,478 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 6.8 x 7.2 | 7.37 | 5.30 | 2 / 39 | 1.4 | none | 0 / 0 | above |
| the preset measured (not an option) | 6.1 x 7.1 | 7.37 | 5.30 | 2 / 38 | 1.5 | none | 0 / 0 | above |

### 9GTP (22,650 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 15.1 x 14.2 | 6.00 | 4.32 | 20 / 148 | 3.8 | none | 0 / 0 | above |
| the preset measured (not an option) | 15.1 x 13.4 | 6.00 | 4.32 | 11 / 139 | 3.7 | none | 0 / 0 | above |

### 6ET5

Notes: below the readability floor: its text is 2.4 pt, under 6 pt, on its 48 x 6 in page (its smallest text 1.4 pt); drawing less may bring it above: --focus C (chain C in full, with what its neighbours add to its fold); --chains C (chain C only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 5.8 | 2.46 | 1.48 | 422 / 422 | 5.5 | helix_name_moved 1, page_text_small 1 | 0 / 0 | text_small |
| --focus C (offered) | 7.0 x 5.7 | 7.37 | 4.42 | 6 / 30 | 1.1 | none | 0 / 0 | above |
| --chains C (offered) | 5.4 x 6.1 | 7.37 | 4.42 | 6 / 27 | 0.7 | none | 0 / 0 | above |
| the preset measured (not an option) | 48.0 x 8.1 | 3.30 | 2.38 | 344 / 344 | 3.5 | counts_hidden 1, page_text_small 1 | 0 / 0 | text_small |

### 8WMV

Notes: below the readability floor: its text is 5.5 pt, under 6 pt, on its 48 x 13 in page (its smallest text 3.3 pt); drawing less may bring it above: --focus A (chain A in full, with what its neighbours add to its fold); --chains A (chain A only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 12.8 | 5.53 | 3.31 | 737 / 737 | 35.4 | page_text_small 1 | 0 / 0 | text_small |
| --focus A (offered) | 9.1 x 6.8 | 6.00 | 3.60 | 79 / 120 | 2.8 | none | 0 / 0 | above |
| --chains A (offered) | 9.1 x 6.1 | 6.00 | 3.60 | 78 / 118 | 2.2 | none | 0 / 0 | above |
| the preset measured (not an option) | 48.0 x 13.0 | 5.77 | 4.16 | 322 / 322 | 15.4 | page_text_small 1 | 0 / 0 | text_small |

### 7NZ2

Notes: below the readability floor: its text is 5.2 pt, under 6 pt, on its 48 x 37 in page (its smallest text 3.1 pt); drawing less may bring it above: --focus A1 (chain A1 in full, with what its neighbours add to its fold); --chains A1 (chain A1 only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 37.3 | 5.21 | 3.13 | 1029 / 1029 | 81.2 | page_text_small 1 | 0 / 0 | text_small |
| --focus A1 (offered) | 20.1 x 22.2 | 6.00 | 3.60 | 21 / 96 | 5.2 | none | 0 / 0 | above |
| --chains A1 (offered) | 20.1 x 22.3 | 6.00 | 3.60 | 20 / 90 | 3.4 | none | 0 / 0 | above |
| the preset measured (not an option) | 48.0 x 36.2 | 4.95 | 3.57 | 929 / 929 | 67.8 | page_text_small 1 | 0 / 0 | text_small |

### 8ZTU

Notes: below the readability floor: its text is 5.7 pt, under 6 pt, on its 48 x 42 in page (its smallest text 3.4 pt); drawing less may bring it above: --focus B (chain B in full, with what its neighbours add to its fold); --chains B (chain B only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 41.7 | 5.75 | 3.45 | 1209 / 1209 | 88.4 | counts_hidden 1, glyph_width_guessed 76, page_text_small 1 | 0 / 1 | text_small |
| --focus B (offered) | 11.0 x 9.7 | 6.00 | 4.50 | 31 / 53 | 5.1 | glyph_width_guessed 40 | 0 / 0 | above |
| --chains B (offered) | 11.0 x 9.4 | 6.00 | 4.50 | 30 / 51 | 4.3 | glyph_width_guessed 40 | 0 / 0 | above |
| the preset measured (not an option) | 48.0 x 41.6 | 5.89 | 4.24 | 1159 / 1159 | 79.7 | glyph_width_guessed 57, helix_name_moved 1, page_text_small 1 | 0 / 1 | text_small |

### 8FZN (30,000 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 6.8 x 9.0 | 7.37 | 5.30 | 2 / 52 | 12.2 | none | 0 / 0 | above |
| the preset measured (not an option) | 5.2 x 9.0 | 7.37 | 5.30 | 2 / 48 | 12.2 | none | 0 / 0 | above |

### 7SPI (23,478 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 7.2 x 5.3 | 6.79 | 4.89 | 3 / 39 | 1.4 | none | 0 / 0 | above |
| the preset measured (not an option) | 7.2 x 5.3 | 6.79 | 4.89 | 3 / 38 | 1.5 | none | 0 / 0 | above |

### 9HEO (21,040 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 20.3 x 25.3 | 6.00 | 3.60 | 28 / 307 | 10.7 | none | 0 / 0 | above |
| the preset measured (not an option) | 21.3 x 25.1 | 6.00 | 4.32 | 6 / 285 | 14.6 | none | 0 / 0 | above |
