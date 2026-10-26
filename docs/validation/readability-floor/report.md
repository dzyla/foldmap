# The readability floor

Plan 28 Task 4 (the review's 3.3, ruling Q39). Written by: foldmap 0.1.0, git eac6756ee494 (asserted by runner, dirty False), code fingerprint 0fe34f04e08547af, scripts/readability_floor.py sha256 4e6ec669cc96; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

A figure below the readability floor says so (floor.readability_floor): in `foldmap plot`'s notes (and with `--strict-floor` it writes no figure and exits 1), in `foldmap check` as R1 (readability) and its --json report, in the MCP result and under the app's figure, naming the existing options that draw less. The figure itself is unchanged. This record is what the reasons and their limits rest on.

## The reasons and their limits

- **text_small**: the figure's own text (element names, chain ends) under 6 pt on its page (layout.MIN_TEXT_PT, the random-100 C4 floor): the page grows to keep it at 6 pt until layout.MAX_PAGE_IN (48 in), and past that the whole figure is scaled down. Measured here at the default page size, as `foldmap plot` writes it. The smallest text drawn (residue counts, base letters, ligand names, 0.6-0.9 of the figure's text by design) is recorded but not judged: it is under 6 pt on 17 of the 18 default figures measured, those far above the floor included.
- **route_budget**: the page's budget for second, wider route searches (route._MAX_NODES_PAGE) spent, counted as the named fallbacks route_second_search_refused. It fired on none of the figures measured here and on none of the samples' default views (census below: their reports count none refused), so it is tested on synthetic measures only; a stuck loop of another kind is A3's, as it happens on small pages too.
- **time**: longer than the hosted app's 90 s (criteria.LIMITS app_s) to make, measured as spec.draw_spec's wall time. It depends on the machine and its load: the figures here were drawn 12 at a time.

## Census of the three samples' default views

From each sample's committed campaign report, at the commit it names (not this code): the default views whose notes say the page held its text under 6 pt, those over 60 s (A2) and over 90 s (A2b), and those with a second route search refused for the page's budget.

| sample | commit | default views | text under 6 pt | over 60 s | over 90 s | route budget spent |
|---|---|---|---|---|---|---|
| random100 | 76969474ca5b | 100 | 7JU4 (5.5), 6JLN (4.5) | 7Y7C (66.92), 8VWQ (64.51) | none | none |
| random100-heldout | 76969474ca5b | 100 | 6ET5 (2.4), 8WMV (5.5), 7NZ2 (5.2), 8ZTU (5.7) | 8ZTU (60.74) | none | none |
| random100-third | 76969474ca5b | 100 | 9FZL (2.0), 8XJV (5.4) | 8XJV (129.02) | 8XJV (129.02) | none |

## What each option achieves on the figures below the floor

Below the floor at this code: 8 of the 18 figures measured, 9FZL (text_small), 8XJV (text_small, time), 7JU4 (text_small), 6JLN (text_small), 6ET5 (text_small), 8WMV (text_small), 7NZ2 (text_small), 8ZTU (text_small).

| option | tried on | brought above | median seconds | median text pt | entries |
|---|---|---|---|---|---|
| `--focus protomer` | 2 | 2 | 9.0 | 6.00 | 9FZL, 6JLN |
| `--focus CHAIN` | 8 | 8 | 2.3 | 6.00 | 9FZL, 8XJV, 7JU4, 6JLN, 6ET5, 8WMV, 7NZ2, 8ZTU |
| `--chains CHAIN` | 8 | 8 | 0.9 | 6.00 | 9FZL, 8XJV, 7JU4, 6JLN, 6ET5, 8WMV, 7NZ2, 8ZTU |
| `--reduced-detail` | 8 | 0 | 33.8 | 5.24 | 9FZL, 8XJV, 7JU4, 6JLN, 6ET5, 8WMV, 7NZ2, 8ZTU |

`--reduced-detail` is the preset `--set nucleotide_labels=false --set ligands=none --set glycans=false --set helices_310=false` (existing style keys only). It leaves out low-priority annotations, and so most of the text drawn under 6 pt, but the page's size is set by its elements, so it does not raise the figure's text: the notes never offer it for the floor.

## Every figure

Page in inches; the figure's text and its smallest text in points at that page size; texts drawn under 6 pt of all texts; seconds to make (draw_spec); named fallbacks taken; A5 element overlaps / A6 text collisions; the floor's reasons. Offered: the options the notes name for that figure.

### 9FZL

Notes: below the readability floor: its text is 2.0 pt, under 6 pt, on its 48 x 10 in page (its smallest text 1.2 pt); drawing less may bring it above: --focus protomer (one subunit in full); --focus 2M (chain 2M in full, with what its neighbours add to its fold); --chains 2M (chain 2M only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 9.8 | 2.05 | 1.23 | 1442 / 1442 | 53.6 | counts_hidden 2, helix_name_moved 22, page_text_small 1 | 0 / 0 | text_small |
| --focus protomer (offered) | 45.5 x 20.9 | 6.00 | 3.60 | 179 / 750 | 14.1 | counts_hidden 1 | 0 / 0 | above |
| --focus 2M (offered) | 8.8 x 7.2 | 6.00 | 3.60 | 10 / 52 | 2.4 | none | 0 / 0 | above |
| --chains 2M (offered) | 8.8 x 6.5 | 6.00 | 3.60 | 9 / 50 | 0.9 | none | 0 / 0 | above |
| --reduced-detail | 48.0 x 9.3 | 2.06 | 1.48 | 1195 / 1195 | 50.9 | counts_hidden 7, helix_name_moved 11, page_text_small 1 | 0 / 0 | text_small |

### 8XJV

Notes: below the readability floor: its text is 5.4 pt, under 6 pt, on its 47 x 48 in page (its smallest text 3.9 pt); it took 110 s to make, over the hosted app's 90 s; drawing less may bring it above: --focus A (chain A in full, with what its neighbours add to its fold); --chains A (chain A only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 46.9 x 48.0 | 5.48 | 3.95 | 1288 / 1288 | 109.9 | counts_hidden 5, glyph_width_guessed 76, helix_name_moved 19, na_gap_unclear 3, na_tail_straight 4, na_trace_unclear 3, page_text_small 1 | 0 / 3 | text_small, time |
| --focus A (offered) | 9.0 x 10.2 | 6.00 | 4.32 | 22 / 42 | 13.4 | glyph_width_guessed 20 | 0 / 0 | above |
| --chains A (offered) | 8.1 x 10.5 | 6.00 | 4.32 | 21 / 33 | 2.0 | glyph_width_guessed 20 | 0 / 0 | above |
| --reduced-detail | 45.6 x 48.0 | 5.19 | 3.74 | 1257 / 1257 | 113.7 | counts_hidden 13, glyph_width_guessed 57, helix_name_moved 15, na_gap_unclear 4, na_tail_straight 3, na_trace_unclear 3, page_text_small 1 | 0 / 5 | text_small, time |

### 7QUP (29,575 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 7.8 x 10.0 | 6.00 | 3.60 | 12 / 73 | 23.9 | none | 0 / 0 | above |
| --reduced-detail | 7.7 x 10.7 | 6.00 | 5.10 | 1 / 62 | 24.9 | none | 0 / 0 | above |

### 6YF9 (25,200 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 5.2 x 4.3 | 7.37 | 6.26 | 0 / 15 | 3.9 | none | 0 / 0 | above |
| --reduced-detail | 5.2 x 4.3 | 7.37 | 6.26 | 0 / 15 | 4.3 | none | 0 / 0 | above |

### 7JU4 (23,990 residues deposited)

Notes: below the readability floor: its text is 5.5 pt, under 6 pt, on its 33 x 48 in page (its smallest text 3.3 pt); drawing less may bring it above: --focus 0 (chain 0 in full, with what its neighbours add to its fold); --chains 0 (chain 0 only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 32.8 x 48.0 | 5.54 | 3.32 | 1298 / 1298 | 38.1 | page_text_small 1 | 0 / 0 | text_small |
| --focus 0 (offered) | 2.5 x 12.2 | 7.37 | 5.30 | 1 / 9 | 2.9 | none | 0 / 0 | above |
| --chains 0 (offered) | 1.8 x 11.0 | 7.37 | 5.30 | 1 / 5 | 0.9 | none | 0 / 0 | above |
| --reduced-detail | 34.3 x 48.0 | 5.87 | 4.22 | 1162 / 1162 | 30.7 | page_text_small 1 | 0 / 0 | text_small |

### 6JLN

Notes: below the readability floor: its text is 4.5 pt, under 6 pt, on its 48 x 9 in page (its smallest text 2.7 pt); drawing less may bring it above: --focus protomer (one subunit in full); --focus A (chain A in full, with what its neighbours add to its fold); --chains A (chain A only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 9.3 | 4.53 | 2.72 | 584 / 584 | 17.9 | page_text_small 1 | 0 / 0 | text_small |
| --focus protomer (offered) | 26.4 x 11.3 | 6.00 | 3.60 | 116 / 300 | 4.0 | none | 0 / 0 | above |
| --focus A (offered) | 7.2 x 5.1 | 6.00 | 3.60 | 29 / 49 | 0.7 | none | 0 / 0 | above |
| --chains A (offered) | 6.3 x 5.4 | 7.37 | 4.42 | 23 / 42 | 0.5 | none | 0 / 0 | above |
| --reduced-detail | 48.0 x 14.4 | 5.29 | 3.81 | 389 / 389 | 6.5 | page_text_small 1 | 0 / 0 | text_small |

### 7Y7C

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 47.1 x 41.8 | 6.00 | 3.72 | 654 / 1201 | 45.0 | counts_hidden 4, glyph_width_guessed 100 | 0 / 0 | above |
| --reduced-detail | 48.0 x 38.9 | 5.37 | 3.87 | 1155 / 1155 | 50.6 | counts_hidden 4, glyph_width_guessed 75, na_run_unclear 2, page_text_small 1 | 0 / 1 | text_small |

### 8VWQ

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 42.3 | 6.00 | 3.60 | 657 / 1243 | 45.5 | counts_hidden 3, glyph_width_guessed 76 | 0 / 0 | above |
| --reduced-detail | 46.4 x 41.8 | 6.00 | 4.32 | 603 / 1189 | 45.8 | counts_hidden 2, glyph_width_guessed 57, helix_name_moved 1 | 0 / 2 | above |

### 8D3D (23,504 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 15.2 x 8.1 | 6.00 | 3.72 | 20 / 128 | 16.1 | none | 0 / 0 | above |
| --reduced-detail | 15.7 x 8.9 | 6.00 | 4.32 | 5 / 113 | 14.4 | none | 0 / 0 | above |

### 7SPB (23,478 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 6.7 x 8.1 | 7.37 | 5.30 | 3 / 48 | 0.7 | none | 0 / 0 | above |
| --reduced-detail | 6.8 x 8.0 | 7.37 | 5.30 | 3 / 47 | 0.7 | none | 0 / 0 | above |

### 9GTP (22,650 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 19.4 x 15.6 | 6.00 | 4.32 | 27 / 194 | 2.2 | none | 0 / 0 | above |
| --reduced-detail | 15.8 x 15.0 | 6.00 | 4.32 | 15 / 182 | 1.9 | none | 0 / 0 | above |

### 6ET5

Notes: below the readability floor: its text is 2.4 pt, under 6 pt, on its 48 x 6 in page (its smallest text 1.4 pt); drawing less may bring it above: --focus C (chain C in full, with what its neighbours add to its fold); --chains C (chain C only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 5.8 | 2.46 | 1.48 | 422 / 422 | 2.9 | helix_name_moved 1, page_text_small 1 | 0 / 0 | text_small |
| --focus C (offered) | 7.0 x 5.7 | 7.37 | 4.42 | 6 / 30 | 0.8 | none | 0 / 0 | above |
| --chains C (offered) | 5.4 x 6.1 | 7.37 | 4.42 | 6 / 27 | 0.3 | none | 0 / 0 | above |
| --reduced-detail | 48.0 x 8.1 | 3.30 | 2.38 | 344 / 344 | 1.8 | counts_hidden 1, page_text_small 1 | 0 / 0 | text_small |

### 8WMV

Notes: below the readability floor: its text is 5.5 pt, under 6 pt, on its 48 x 13 in page (its smallest text 3.3 pt); drawing less may bring it above: --focus A (chain A in full, with what its neighbours add to its fold); --chains A (chain A only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 12.8 | 5.53 | 3.31 | 737 / 737 | 14.1 | page_text_small 1 | 0 / 0 | text_small |
| --focus A (offered) | 9.1 x 6.8 | 6.00 | 3.60 | 79 / 120 | 1.2 | none | 0 / 0 | above |
| --chains A (offered) | 9.1 x 6.1 | 6.00 | 3.60 | 78 / 118 | 0.9 | none | 0 / 0 | above |
| --reduced-detail | 48.0 x 13.0 | 5.77 | 4.16 | 322 / 322 | 6.9 | page_text_small 1 | 0 / 0 | text_small |

### 7NZ2

Notes: below the readability floor: its text is 5.2 pt, under 6 pt, on its 48 x 37 in page (its smallest text 3.1 pt); drawing less may bring it above: --focus A1 (chain A1 in full, with what its neighbours add to its fold); --chains A1 (chain A1 only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 37.3 | 5.21 | 3.13 | 1029 / 1029 | 43.8 | page_text_small 1 | 0 / 0 | text_small |
| --focus A1 (offered) | 20.1 x 22.2 | 6.00 | 3.60 | 21 / 96 | 2.6 | none | 0 / 0 | above |
| --chains A1 (offered) | 20.1 x 22.3 | 6.00 | 3.60 | 20 / 90 | 1.4 | none | 0 / 0 | above |
| --reduced-detail | 48.0 x 36.2 | 4.95 | 3.57 | 929 / 929 | 36.9 | page_text_small 1 | 0 / 0 | text_small |

### 8ZTU

Notes: below the readability floor: its text is 5.7 pt, under 6 pt, on its 48 x 42 in page (its smallest text 3.4 pt); drawing less may bring it above: --focus B (chain B in full, with what its neighbours add to its fold); --chains B (chain B only)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 48.0 x 41.7 | 5.75 | 3.45 | 1209 / 1209 | 44.3 | counts_hidden 1, glyph_width_guessed 76, page_text_small 1 | 0 / 1 | text_small |
| --focus B (offered) | 11.0 x 9.7 | 6.00 | 4.50 | 31 / 53 | 2.3 | glyph_width_guessed 40 | 0 / 0 | above |
| --chains B (offered) | 11.0 x 9.4 | 6.00 | 4.50 | 30 / 51 | 1.8 | glyph_width_guessed 40 | 0 / 0 | above |
| --reduced-detail | 48.0 x 41.6 | 5.89 | 4.24 | 1159 / 1159 | 39.7 | glyph_width_guessed 57, helix_name_moved 1, page_text_small 1 | 0 / 1 | text_small |

### 8FZN (30,000 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 6.8 x 9.0 | 7.37 | 5.30 | 2 / 52 | 6.3 | none | 0 / 0 | above |
| --reduced-detail | 5.2 x 9.0 | 7.37 | 5.30 | 2 / 48 | 6.6 | none | 0 / 0 | above |

### 7SPI (23,478 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 6.7 x 8.2 | 7.37 | 5.30 | 3 / 48 | 0.6 | none | 0 / 0 | above |
| --reduced-detail | 7.2 x 6.3 | 6.33 | 4.56 | 4 / 47 | 0.6 | none | 0 / 0 | above |

### 9HEO (21,040 residues deposited)

| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions | floor |
|---|---|---|---|---|---|---|---|---|
| (default) | 20.3 x 25.3 | 6.00 | 3.60 | 28 / 307 | 5.6 | none | 0 / 0 | above |
| --reduced-detail | 21.3 x 25.1 | 6.00 | 4.32 | 6 / 285 | 7.1 | none | 0 / 0 | above |
