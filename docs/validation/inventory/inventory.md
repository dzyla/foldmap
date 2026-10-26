# foldmap inventory: thresholds, special cases and fallbacks

Plan 28 Task 9.1 (ruling Q23), the review's "stabilize and characterize": every number the code decides by, and every
fallback it takes, listed with where it comes from. `scripts/inventory.py` scans `src/foldmap` (an AST scan) and writes
`inventory.json`; this table is that scan, completed by hand. `tests/test_inventory.py` fails when the scan finds an
item this table lacks, or this table names one the scan no longer finds: run `python scripts/inventory.py` after a
change, then complete the new rows. The merge keeps rows by id with their hand columns; a row whose line's code was
edited (its id changed) is carried to the new id when the two alone share `module:function:kind:value` ("carried");
other new rows take the origin their number's role fixes (`role_origin`), else `unknown`; gone rows are dropped. Each
is printed.

`inventory.json` holds every item with its line, comments and tests as of the commit its source stamp names
(`code_fingerprint` identifies the code exactly); only ids, kinds and values are checked to be fresh, so its line
numbers drift as the code moves until the script is run again.

## What is listed

- **constant**: a module-level assignment holding a number (or a class's UPPER_CASE one), id `module:NAME`.
- **default**: a function parameter's non-trivial numeric default, id `module:function:default:param`.
- **local**: a name bound in a function to numbers alone, id `module:function:local:name`.
- **comparison**, **argument**: a non-trivial number inside a comparison or a call's argument (with items and except
  clauses included), id `module:function:kind:value@anchor`. The anchor is 8 hex digits of the SHA-1 of the item's
  line of code with its comment cut off and all whitespace and a trailing comma dropped, so the id survives the line
  moving, its comment being reworded and reformatting within the line; `#k` for the k-th repeat in one function.
- **fallback**: a `fallback.note(NAME)` call, id `module:function:fallback:NAME`: a named fallback, counted per figure.
  A note whose name is no string literal fails the test (the scan could not read it).
- **fallback-marker**: a function whose text speaks of falling back (fallback, falls back, fell back, last resort,
  gives up, degraded), id `module:function:fallback-marker`. Its evidence says which named fallback it holds, or why
  none.

Trivial numbers (0, 1, 2, -1, 0.5) are not listed inline. Not listed at all (the scan's limits): numbers in arithmetic
outside a comparison or a call (a bare `x * 0.35`, a return value), indices and slices, numbers in f-strings, a
class's keywords, match patterns, dataclass field defaults (the style keys: `style.py` documents each), and
string-valued switches (a special case written as `if kind == "tail"` holds no number). A special case is therefore
listed by the numbers it compares, not by its condition. Fallbacks are counted only in figures made by
`spec.draw_spec` (the command line's plot, the app's figure, the MCP server, random-100): compare and diff panels,
interactive pages and `cli.make_figure_and_loops` count none.

## Columns

- **value**: as the source writes it (a constant's whole value, shortened).
- **at**: the item's source line, shortened.
- **origin**, by the evidence at hand, never guessed:
  - `biology`: a physical, chemical or structural-biology fact, or a parameter of a published method (DSSP's energies,
    TM-align's penalties, bond lengths, hydrogen-bond geometry), cited at the item;
  - `geometry`: follows from geometry or from other constants (a half cell, a right angle, a copy of another module's
    value), as the comment derives it;
  - `measured`: the comment says it was measured, calibrated or frozen on a set (the set's entries in the evidence);
  - `named-entry:<ids>`: the comment names the structure(s) the value was set or fixed on, or (evidence starting
    `git:`) the message of the commit that set this value ties it to them; a commit that only names entries in its
    subject, without tying this value to them, leaves the origin unknown and stays in the evidence;
  - `convention`: fixed by its role, not chosen to make a figure come out: units and unit conversions, numeric
    precision and float-noise tolerances, array shapes, rounding digits, software limits (sizes, timeouts, caches,
    work bounds), schema versions, colours and drawing order;
  - `unknown`: none of these could be established from the comment, the ledgers or git. A comment that cites only a
    campaign finding, a ruling or a review item (random-100 C26, ruling P35) leaves the origin unknown; the evidence
    names what it cites.
- **evidence**: the words of the comment the origin rests on (quoted, shortened), and for a constant the earliest
  commit that has its name in this file (`git log -S`: the commit that introduced it, or a rename or split that moved
  it here).
- **tests naming it**: test files that name both the item's module and the constant or the function it is in, or (a
  fallback) its quoted name; blank for a name under 4 characters or a common word. Grep evidence that a test exercises
  it, not proof that a test pins the value.
- **check**: for a fallback, the check or note that reports when it fires; `-` where none does or it does not apply.

## How the origins were filled

By rule from the comment at each item, then corrected by hand for every constant and fallback: `measured` when the
comment says measured on/over/at/in/against, calibrated, frozen or tuning set; else `named-entry` when it names PDB
entries; `convention` or `geometry` by role (`role_origin` in `scripts/inventory.py`, applied to every new row too: a
float tolerance of 1e-5 or less, rounding digits, an array's shape, `axis=`, `zorder=`, a JSON indent, 72 points to
the inch, byte sizes; 90, 180 and 360 degrees turned to radians); `unknown` otherwise. The constants' `biology`,
`geometry` and `convention` origins, the git-derived `named-entry` rows (each read against its commit's message) and
every correction of a rule are by hand.

<!-- inventory:generated (everything below is rewritten by scripts/inventory.py) -->

## Counts

| kind | items | biology | geometry | measured | named-entry | convention | unknown |
|---|---|---|---|---|---|---|---|
| constant | 590 | 62 | 25 | 28 | 36 | 100 | 339 |
| default | 33 | 0 | 0 | 0 | 0 | 3 | 30 |
| local | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| comparison | 297 | 0 | 0 | 0 | 0 | 154 | 143 |
| argument | 663 | 0 | 0 | 0 | 0 | 182 | 481 |
| fallback | 24 | 0 | 0 | 0 | 4 | 1 | 19 |
| fallback-marker | 29 | 0 | 0 | 0 | 3 | 0 | 26 |
| all | 1637 | 62 | 25 | 28 | 43 | 440 | 1039 |

Of the 43 named-entry rows, 32 rest on the comment at the item and 11 on git (evidence starting `git:`: a commit whose message ties the value to the entries).

## Named fallbacks (foldmap.fallback.FALLBACKS)

| name | sites | what it is |
|---|---|---|
| `loop_curve` | 1 | a loop no orthogonal route was found for, drawn as a smooth curve (route.route_loops; counted over the loops as drawn, a rosette's turned copies included) |
| `loop_over_text` | 1 | a loop whose only route runs over a label or a chain end (route._Router.route, its keep-outs opened; counted over the loops as drawn, once the names have moved off lines) |
| `loop_wall_detour` | 1 | a membrane view's fallback curve that would cross its loop's wall, drawn instead as a detour out of its ports and round the page (route._walled_fallback, ruling P35) |
| `loop_wall_unkept` | 1 | ... and no detour keeps to the wall either: the last one tried is drawn, logged (route._walled_fallback) |
| `route_unsnapped` | 1 | a route whose drawing snapped to its ports crossed an element, drawn unsnapped (route._Router.drawn; every drawing counted, those of searches set aside included) |
| `route_second_search_refused` | 1 | a wider second route search not run because the page's node budget was spent (route._Router.second, _MAX_NODES_PAGE) |
| `helix_name_moved` | 1 | a helix name a line runs through or near, moved to the spot lines run nearest least (route._clear_helix_names, the fallback for lines no route holds) |
| `na_tail_straight` | 1 | a free nucleic-acid end with no clear drawing, drawn straight on (na_ends._tail) |
| `na_gap_unclear` | 1 | a cropped gap's dashed line with no drawing clear of the elements (na_ends._gap) |
| `na_hairpin_unclear` | 1 | a hairpin loop's arc not clear of the elements (na_runs._strand_between) |
| `na_run_unclear` | 1 | a single strand or joint between stems with no clear route, drawn as a curve (na_runs._drawing) |
| `na_trace_unclear` | 1 | an unpaired strand's trace not clear of the elements, a nucleotide no nudge clears kept (route_na.route_nucleic, na_ends.off_unheld) |
| `stem_name_least_bad` | 1 | a stem's name with no clear slot nearer its own stem with room, put in the slot with the most room (layout._choose_label_sides, _pick_slot; logged) |
| `domain_name_least_bad` | 1 | a domain's name with no clear place, put where it leaves the loops' ways out the most room (layout._place_domain_names) |
| `glyph_width_guessed` | 3 | a text measured with a character the built-in metrics lack, its width guessed (layout._text_width and _ink_line: _FALLBACK_WIDTH; render._clear_of_lines: half an em); every measuring counted |
| `mark_least_bad_spot` | 1 | a ligand mark with no clear spot, drawn over other ink at the least bad one (render._least_bad_spot; counted over the marks drawn) |
| `counts_hidden` | 1 | a residue count of a long loop or chain end with no clear spot, left out (render.tail_counts, loop_counts) |
| `interactions_dropped` | 1 | an interaction line the independent re-check (interactions.verify_interactions) found breaking its rule, not drawn |
| `interactions_hidden` | 1 | an interaction line that would run wholly under text, not drawn (render._draw_interactions) |
| `crosslinks_unplaced` | 1 | a crosslink with an end not drawn or both ends at one point, not drawn (render) |
| `page_text_small` | 1 | a page held at its largest size, its text under 6 pt (layout.page_fit, render._canvas) |
| `s3_cycles_capped` | 1 | S3 stopped counting a sheet's pairing cycles at _MAX_CYCLES: the barrel check says there are too many to count (science.pairing_cycle) |

## align

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `align:GRID` | constant | 1024 | GRID = 1024 # distances are snapped to 1/GRID Å before scoring: exact in binary floating point, the same on any CPU | convention | distances snapped to 1/GRID Å: exact in binary floating point; earliest in this file: b183f60 2026-10-21 Chain aligner, sequence and file routes (plan 13a Task 1): align.py with ChainAl… | align | - |
| `align:_SCORE` | constant | 2.0**20 | _SCORE = 2.0**20 # pair scores are snapped to 1/_SCORE: sums of up to 2^33 of them are exact in a double | convention | scores snapped to 1/_SCORE: exact sums in a double; earliest in this file: a4c21a0 2026-10-21 Chain aligner, structure route (plan 13a Task 2): TM-align's heuristics on Cα in… | align | - |
| `align:GAP_OPEN` | constant | -0.6 | GAP_OPEN = -0.6 # TM-align's gap-open penalty (Zhang & Skolnick 2005); no extension penalty | biology | "TM-align's gap-open penalty (Zhang & Skolnick 2005)"; earliest in this file: a4c21a0 2026-10-21 Chain aligner, structure route (plan 13a Task 2): TM-align's heuristics on Cα in… | - | - |
| `align:SS_GAP` | constant | -1.0 | SS_GAP = -1.0 # the gap penalty of TM-align's secondary-structure seed (a match scores 1) | biology | "the gap penalty of TM-align's secondary-structure seed"; earliest in this file: a4c21a0 2026-10-21 Chain aligner, structure route (plan 13a Task 2): TM-align's heuristics on Cα in… | - | - |
| `align:ROUNDS` | constant | 20 | ROUNDS = 20 # superpose / score / re-align rounds per seed, at most; also the iterations of one fragment's search | unknown | earliest in this file: a4c21a0 2026-10-21 Chain aligner, structure route (plan 13a Task 2): TM-align's heuristics on Cα in… | - | - |
| `align:ROUND_STEP` | constant | 40 | ROUND_STEP = 40 # fragment starts every ROUND_STEP pairs in the rounds' search (TM-align's simplify_step); 1 at the end | biology | "TM-align's simplify_step"; earliest in this file: a4c21a0 2026-10-21 Chain aligner, structure route (plan 13a Task 2): TM-align's heuristics on Cα in… | - | - |
| `align:SEQUENCE_IDENTITY` | constant | 0.4 | SEQUENCE_IDENTITY = 0.4 # auto: the sequence route at this identity over the aligned residues, else the structure | unknown | earliest in this file: a4c21a0 2026-10-21 Chain aligner, structure route (plan 13a Task 2): TM-align's heuristics on Cα in… | - | - |
| `align:SEQUENCE_ROUTE_TM` | constant | 0.5 | SEQUENCE_ROUTE_TM = 0.5 | unknown | earliest in this file: 1a3f9d6 2026-10-21 Like.order: a mated sheet reads left to right as the reference reads its mates (… | - | - |
| `align:SAME_FOLD_TM` | constant | 0.5 | SAME_FOLD_TM = 0.5 # H1: one fold at this TM-score (normalised by A) ... (frozen on the tuning set, 13a-10) | measured | "H1: one fold at this TM-score (normalised by A) ... (frozen on the tuning set, 13a-10)"; earliest in this file: ff1d90f 2026-10-21 Element mates (plan 13a Task 4): Mates and element_mates in align.py | align, sensitivity | - |
| `align:SAME_FOLD_FRACTION` | constant | 0.5 | SAME_FOLD_FRACTION = 0.5 # ... with this fraction of the chain with fewer elements matched (frozen, 13a-10) | measured | "…ith this fraction of the chain with fewer elements matched (frozen, 13a-10)"; earliest in this file: ff1d90f 2026-10-21 Element mates (plan 13a Task 4): Mates and element_mates in align.py | - | - |
| `align:LOCAL_FRAGMENTS` | constant | (20, 100) | LOCAL_FRAGMENTS = (20, 100) | biology | "TM-align's local-superposition seed (get_initial5)"; earliest in this file: 28ab8e2 2026-10-22 TM-align's other seeds in the aligner (plan 24 Task 1): local superposition, SS … | - | - |
| `align:LOCAL_CELLS` | constant | 2**26 | LOCAL_CELLS = 2**26 | convention | a work bound; earliest in this file: 28ab8e2 2026-10-22 TM-align's other seeds in the aligner (plan 24 Task 1): local superposition, SS … | align | - |
| `align:LOCAL_BATCH` | constant | 2**22 | LOCAL_BATCH = 2**22 # ... aligned LOCAL_BATCH residue pairs at a time (memory; the result does not depend on it) | convention | "memory; the result does not depend on it"; earliest in this file: 28ab8e2 2026-10-22 TM-align's other seeds in the aligner (plan 24 Task 1): local superposition, SS … | - | - |
| `align:WIDER` | constant | 1.5 | WIDER = 1.5 # get_initial5's and initial_ssplus's distance scale: d0 + WIDER Å | biology | "get_initial5's and initial_ssplus's distance scale" (TM-align); earliest in this file: 28ab8e2 2026-10-22 TM-align's other seeds in the aligner (plan 24 Task 1): local superposition, SS … | - | - |
| `align:SS_BONUS` | constant | 0.5 | SS_BONUS = 0.5 # initial_ssplus: a residue pair's score gains this when the two residues' H/E/C agree | biology | "initial_ssplus" (TM-align); earliest in this file: 28ab8e2 2026-10-22 TM-align's other seeds in the aligner (plan 24 Task 1): local superposition, SS … | - | - |
| `align:SCREEN` | constant | 0.0 | SCREEN = 0.0 # 0, 0.003 and 0.01 chose the same alignments on the 55 homologue pairs (plan 24) | measured | "0, 0.003 and 0.01 chose the same alignments on the 55 homologue pairs (plan 24)"; earliest in this file: 28ab8e2 2026-10-22 TM-align's other seeds in the aligner (plan 24 Task 1): local superposition, SS … | sensitivity | - |
| `align:d0:comparison:15@8d42232c` | comparison | 15 | if length <= 15: # the cube root of a negative number is negative: below 0.5 anyway | unknown | - | - | - |
| `align:d0:argument:1.24@abb321f0` | argument | 1.24 | return max(0.5, 1.24 * (length - 15) ** (1 / 3) - 1.8) | unknown | - | - | - |
| `align:d0:argument:15@abb321f0` | argument | 15 | return max(0.5, 1.24 * (length - 15) ** (1 / 3) - 1.8) | unknown | - | - | - |
| `align:d0:argument:3@abb321f0` | argument | 3 | return max(0.5, 1.24 * (length - 15) ** (1 / 3) - 1.8) | unknown | - | - | - |
| `align:d0:argument:1.8@abb321f0` | argument | 1.8 | return max(0.5, 1.24 * (length - 15) ** (1 / 3) - 1.8) | unknown | - | - | - |
| `align:_cutoff:argument:4.5@680a5cef` | argument | 4.5 | return min(max(scale, 4.5), 8.0) | unknown | - | - | - |
| `align:_cutoff:argument:8.0@680a5cef` | argument | 8.0 | return min(max(scale, 4.5), 8.0) | unknown | - | - | - |
| `align:_within:argument:3@44dc4035` | argument | 3 | while sel.sum() < min(3, len(d)): | unknown | - | - | - |
| `align:_search:argument:5@bb64cdca` | argument | 5 | for _ in range(5): | unknown | - | - | - |
| `align:_search:comparison:4@c7b04a98` | comparison | 4 | if size < 4: | unknown | - | - | - |
| `align:_search:argument:3@ee7c0ed2` | argument | 3 | best, best_r, best_t = -1.0, np.eye(3), np.zeros(3) | convention | an array's shape | - | - |
| `align:_search:argument:3@ee7c0ed2#2` | argument | 3 | best, best_r, best_t = -1.0, np.eye(3), np.zeros(3) | convention | an array's shape | - | - |
| `align:tm_score:argument:3@e5b4163b` | argument | 3 | return 0.0, np.eye(3), np.zeros(3), 0.0 | convention | an array's shape | align, criteria, homologues, like | - |
| `align:tm_score:argument:3@e5b4163b#2` | argument | 3 | return 0.0, np.eye(3), np.zeros(3), 0.0 | convention | an array's shape | align, criteria, homologues, like | - |
| `align:_dp:argument:3@b31491d7` | argument | 3 | mm, xx, yy = (np.full((n + 1, m + 1), low) for _ in range(3)) | unknown | - | - | - |
| `align:_dp:argument:3@133c6233` | argument | 3 | pm, px, py = (np.zeros((n + 1, m + 1), np.int8) for _ in range(3)) | unknown | - | - | - |
| `align:_matrix:argument:3@7b80932d` | argument | 3 | for k in range(3): # x² + y² + z², summed in that order | unknown | - | align | - |
| `align:_thread:argument:5@4fcb2cd8` | argument | 5 | least = max(min(n, m) // 2, 5) | unknown | - | - | - |
| `align:_thread:comparison:3@48488f3f` | comparison | 3 | if j1 - j0 < 3: | unknown | - | - | - |
| `align:_jump:comparison:250@e5f4c14e` | comparison | 250 | step = 45 if length > 250 else 35 if length > 200 else 25 if length > 150 else 15 | unknown | - | - | - |
| `align:_jump:comparison:200@e5f4c14e` | comparison | 200 | step = 45 if length > 250 else 35 if length > 200 else 25 if length > 150 else 15 | unknown | - | - | - |
| `align:_jump:comparison:150@e5f4c14e` | comparison | 150 | step = 45 if length > 250 else 35 if length > 200 else 25 if length > 150 else 15 | unknown | - | - | - |
| `align:_jump:argument:3@b86a6971` | argument | 3 | return max(1, min(step, length // 3)) | unknown | - | - | - |
| `align:_local_starts:argument:3@4ac04bad` | argument | 3 | sizes = [s for s in (min(LOCAL_FRAGMENTS[0], shorter // 3), min(LOCAL_FRAGMENTS[1], shorter // 2)) if s >= 3] | unknown | - | align | - |
| `align:_local_starts:comparison:3@4ac04bad` | comparison | 3 | sizes = [s for s in (min(LOCAL_FRAGMENTS[0], shorter // 3), min(LOCAL_FRAGMENTS[1], shorter // 2)) if s >= 3] | unknown | - | align | - |
| `align:_local:comparison:3@8fa7bcc5` | comparison | 3 | if len(pairs) < 3: | unknown | - | align | - |
| `align:_ss_plus:comparison:3@8fa7bcc5` | comparison | 3 | if len(pairs) < 3: | unknown | - | align | - |
| `align:_refine:argument:3@850c98f4` | argument | 3 | best = (-1.0, [], np.eye(3), np.zeros(3)) | convention | an array's shape | - | - |
| `align:_refine:argument:3@850c98f4#2` | argument | 3 | best = (-1.0, [], np.eye(3), np.zeros(3)) | convention | an array's shape | - | - |
| `align:_structure_pairs:comparison:3@0238545a` | comparison | 3 | consider(refined, screened=k >= 3) | unknown | - | - | - |

## app

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `app:LAYOUT_DEFAULTS` | constant | { "mode": "projected", "rotate": 0.0, "symmetry": "auto", "assembly": "auto", "rename": "", "swap": "", "move": "", "ti… | LAYOUT_DEFAULTS = { | convention | the app's default settings; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | app, app_controls, app_edits | - |
| `app:TITLE_MAX` | constant | 200 | TITLE_MAX = 200 # characters of a layout file's title | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `app:TEXT_MAX` | constant | 2000 | TEXT_MAX = 2000 # characters of any other text a layout file gives (a path, a mark list, a reference) | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `app:IDS_MAX` | constant | 20_000 | IDS_MAX = 20_000 # characters of a list of chain ids, joined (a capsid has hundreds of chains) | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `app:ROTATE_MAX` | constant | 360.0 | ROTATE_MAX = 360.0 # degrees, each way | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `app:VECTOR_MAX` | constant | 1000.0 | VECTOR_MAX = 1000.0 # each number of up and view | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `app:STYLE_NUMBERS` | constant | {"float": (0.0, 10.0), "int": (1, 100_000)} | STYLE_NUMBERS = {"float": (0.0, 10.0), "int": (1, 100_000)} # scales and line widths; hbond_max | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `app:APP_THEME` | constant | { # Okabe-Ito blue on white, and its dark twin; every colour pair is checked in tests/test_app_look.py "theme.light": {… | APP_THEME = { # Okabe-Ito blue on white, and its dark twin; every colour pair is checked in tests/test_app_look.py | convention | colours; earliest in this file: af6271c 2026-10-12 App look: Okabe-Ito blue light and dark themes, no Deploy menu, small stylesheet… | app_look, app_scope | - |
| `app:APP_CONFIG` | constant | { # what `foldmap app` passes and .streamlit/config.toml holds (Streamlit Community Cloud reads it) "client": {"toolbar… | APP_CONFIG = { # what `foldmap app` passes and .streamlit/config.toml holds (Streamlit Community Cloud reads it) | convention | app settings; earliest in this file: af6271c 2026-10-12 App look: Okabe-Ito blue light and dark themes, no Deploy menu, small stylesheet… | app_look | - |
| `app:SIDEBAR_WIDTH` | constant | 440 | SIDEBAR_WIDTH = 440 # px: the sidebar's segmented choices and toggle rows fit on one line (Streamlit's own: 336) | convention | a UI width; earliest in this file: af6271c 2026-10-12 App look: Okabe-Ito blue light and dark themes, no Deploy menu, small stylesheet… | - | - |
| `app:CACHES` | constant | { # the app process's figure caches: (entries, MB one entry holds at most for the largest hosted figure) "render": (12,… | CACHES = { # the app process's figure caches: (entries, MB one entry holds at most for the largest hosted figure) | convention | cache sizes (sized on 6XII, the largest hosted figure); earliest in this file: 8dbc8db 2026-10-12 App caches sized to about 300 MB at worst: render 12, 3D page 2, sequence 8, kep… | app_jobs | - |
| `app:structure_choices:default:limit` | default | 12 | def structure_choices(source: str, limit: int = 12) -> list[tuple[str, str]]: | unknown | - | app | - |
| `app:sequence_figure:argument:60@f72106f6` | argument | 60 | columns=int(view.get("columns", 60)), | unknown | - | app_jobs | - |
| `app:_vector:comparison:3@04be6430` | comparison | 3 | if not isinstance(parts, list) or len(parts) != 3: | unknown | - | - | - |
| `app:settings_from_document:argument:100.0@46fcb822` | argument | 100.0 | out["symmetry_tol"] = _number("symmetry_tol", given.get("symmetry_tol", TOLERANCE), 0.0, 100.0) | unknown | - | app, app_uploads, domains, layoutfile +2 | - |
| `app:browse_themes.themes:argument:150@ec36494e` | argument | 150 | with st.container(width=150, border=True, gap="xsmall", key=f"tile_{theme}"): | unknown | - | app_scope, domains, style | - |
| `app:browse_themes.themes:argument:150@ff71e954` | argument | 150 | with st.container(width=150, border=True, key="theme_note"): | unknown | - | app_scope, domains, style | - |
| `app:essentials:argument:0.6@4cb10ff9` | argument | 0.6 | st.slider("Text size", 0.6, 2.0, step=0.1, **setting(ss, "font_scale")) | unknown | - | app_controls | - |
| `app:essentials:argument:0.1@4cb10ff9` | argument | 0.1 | st.slider("Text size", 0.6, 2.0, step=0.1, **setting(ss, "font_scale")) | unknown | - | app_controls | - |
| `app:main.sized_of:argument:64@de8fac91` | argument | 64 | @st.cache_data(show_spinner=False, max_entries=64) | unknown | - | - | - |
| `app:main.choices_of:argument:86400@b18715ed` | argument | 86400 | @st.cache_data(show_spinner=False, ttl=86400, max_entries=64) | unknown | - | - | - |
| `app:main.choices_of:argument:64@b18715ed` | argument | 64 | @st.cache_data(show_spinner=False, ttl=86400, max_entries=64) | unknown | - | - | - |
| `app:main.contacts_of:argument:16@e5fbc74d` | argument | 16 | @st.cache_data(show_spinner=False, max_entries=16) | unknown | - | - | - |

## app_controls

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `app_controls:NAMES_FILE_MAX` | constant | 100_000 | NAMES_FILE_MAX = 100_000 # bytes: a names file is a few kB | convention | an upload limit; earliest in this file: c955bea 2026-10-12 App: conventional names (keywords or an uploaded YAML) and domains with colours … | - | - |
| `app_controls:CELL` | constant | 130 | CELL, NARROW = 130, 60 # column widths (px): two cells, or an element and two numbers, fit the sidebar | convention | a UI width; earliest in this file: bc2e50a 2026-10-12 App section Edit layout: arrangement options and table editors for swap, move, r… | - | - |
| `app_controls:NARROW` | constant | 60 | CELL, NARROW = 130, 60 # column widths (px): two cells, or an element and two numbers, fit the sidebar | convention | a UI width; earliest in this file: bc2e50a 2026-10-12 App section Edit layout: arrangement options and table editors for swap, move, r… | - | - |
| `app_controls:MB` | constant | 1024 * 1024 | MB = 1024 * 1024 # bytes in Streamlit's upload megabyte (max_upload_size) | convention | bytes; earliest in this file: e6caacf 2026-10-12 App uploads: layout files capped at 1 MB and alignments at 10 MB, checked before… | - | - |
| `app_controls:ALIGNMENT_MAX` | constant | 10 * MB | ALIGNMENT_MAX = 10 * MB # the alignment uploader's cap, checked again before the file is kept or read | convention | an upload limit; earliest in this file: e6caacf 2026-10-12 App uploads: layout files capped at 1 MB and alignments at 10 MB, checked before… | app | - |
| `app_controls:TABLE_MAX` | constant | 5 * MB | TABLE_MAX = 5 * MB # the data and crosslinks uploaders' cap, checked again before the file is kept or read | convention | an upload limit; earliest in this file: 243697e 2026-10-21 Crosslinks drawn, the app's data upload and MCP data parameters (plan 13d-6): ea… | - | - |
| `app_controls:LIKE_MAX` | constant | 50 * MB | LIKE_MAX = 50 * MB # the reference uploader's cap: the Start uploader's (APP_CONFIG maxUploadSize), checked again | convention | an upload limit; earliest in this file: 6ef541c 2026-10-21 The app's Draw like… and diff (plan 13e-4): reference by PDB id or upload, Compa… | - | - |
| `app_controls:style_section:argument:0.4@db596796` | argument | 0.4 | st.slider("Loop weight", 0.4, 3.0, step=0.1, **setting(ss, "loop_width")) | unknown | - | - | - |
| `app_controls:style_section:argument:3.0@db596796` | argument | 3.0 | st.slider("Loop weight", 0.4, 3.0, step=0.1, **setting(ss, "loop_width")) | unknown | - | - | - |
| `app_controls:style_section:argument:0.1@db596796` | argument | 0.1 | st.slider("Loop weight", 0.4, 3.0, step=0.1, **setting(ss, "loop_width")) | unknown | - | - | - |
| `app_controls:interactions_block:argument:50@929242e9` | argument | 50 | st.number_input("Most H-bonds drawn", min_value=1, step=50, **setting(ss, "hbond_max")) | unknown | - | - | - |
| `app_controls:interactions_block:argument:50@52d2aae6` | argument | 50 | st.number_input("Most lines of each interaction", min_value=1, step=50, **setting(ss, "interaction_max")) | unknown | - | - | - |
| `app_controls:nucleic_block:argument:3.0@bdc1408a` | argument | 3.0 | 3.0, | unknown | - | - | - |
| `app_controls:nucleic_block:argument:0.05@3875bbe2` | argument | 0.05 | step=0.05, | unknown | - | - | - |
| `app_controls:capped:default:most` | default | 24 | def capped(option, most: int = 24) -> str: | unknown | - | app_controls | - |
| `app_controls:layout_section:argument:-180.0@c81b2a87` | argument | -180.0 | st.slider("Rotate (°)", -180.0, 180.0, step=0.5, **setting(ss, "rotate")) | unknown | - | - | - |
| `app_controls:layout_section:argument:180.0@c81b2a87` | argument | 180.0 | st.slider("Rotate (°)", -180.0, 180.0, step=0.5, **setting(ss, "rotate")) | unknown | - | - | - |
| `app_controls:layout_section:argument:0.05@25cc6aac` | argument | 0.05 | st.slider("Helix size", 0.5, 2.0, step=0.05, **setting(ss, "helix_scale")) | unknown | - | - | - |
| `app_controls:layout_section:argument:0.05@86ef1feb` | argument | 0.05 | st.slider("Strand size", 0.5, 2.0, step=0.05, **setting(ss, "strand_scale")) | unknown | - | - | - |
| `app_controls:layout_section:argument:0.05@ae69d91b` | argument | 0.05 | st.slider("DNA size", 0.5, 2.0, step=0.05, **setting(ss, "dna_scale")) | unknown | - | - | - |
| `app_controls:layout_section:argument:0.25@fe4bf90f` | argument | 0.25 | width_column = number("Width", min_value=0.25, max_value=4.0, step=0.25, width=NARROW) | unknown | - | - | - |
| `app_controls:layout_section:argument:4.0@fe4bf90f` | argument | 4.0 | width_column = number("Width", min_value=0.25, max_value=4.0, step=0.25, width=NARROW) | unknown | - | - | - |
| `app_controls:layout_section:argument:0.25@fe4bf90f#2` | argument | 0.25 | width_column = number("Width", min_value=0.25, max_value=4.0, step=0.25, width=NARROW) | unknown | - | - | - |
| `app_controls:domains_block:argument:10@f7c6e7b8` | argument | 10 | ranges = text("Ranges", width=CELL + 10, help="residue ranges and elements, e.g. A:1-80, α3") | unknown | - | - | - |
| `app_controls:domains_block:argument:10@92bfdd58` | argument | 10 | editor(st, ss, "domains", "Domains", {"name": text("Name", width=CELL - 10), "ranges": ranges}) | unknown | - | - | - |
| `app_controls:contact_lines:default:most` | default | 10 | def contact_lines(report: dict \| None, most: int = 10) -> list[str]: | unknown | - | app_controls | - |

## app_edits

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `app_edits:SYMMETRIES` | constant | ( "auto", "off", "helical", *(f"C{n}" for n in range(2, 13)), *(f"D{n}" for n in range(2, 9)), ) | SYMMETRIES = ( | convention | the choices offered; earliest in this file: cb3c663 2026-10-12 app_edits: structured app controls read and write foldmap's reference syntax | app_controls, app_edits | - |
| `app_edits:upload_name:default:limit` | default | 120 | def upload_name(name, limit: int = 120) -> str: | unknown | - | app_edits | - |
| `app_edits:shown:default:limit` | default | 400 | def shown(text, limit: int = 400, path=None, name: str = "", paths: bool = True) -> str: | unknown | - | app, app_controls, app_edits, app_scope | - |

## app_jobs

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `app_jobs:SLOW_SECONDS` | constant | 4.0 | SLOW_SECONDS = 4.0 # a figure slower than this turns automatic redrawing off and brings its files with it (Q-A6) | convention | an app behaviour limit (ruling Q-A6); earliest in this file: e68ea41 2026-10-12 app_jobs: render jobs with time, memory and concurrency limits; failures remembe… | app, app_controls, app_edits, app_jobs +1 | - |
| `app_jobs:PAGE_EMBED_MAX` | constant | 10_000_000 | PAGE_EMBED_MAX = 10_000_000 # a 3D page larger than this is offered as a download only, never embedded | convention | a size limit; earliest in this file: e68ea41 2026-10-12 app_jobs: render jobs with time, memory and concurrency limits; failures remembe… | app | - |
| `app_jobs:BUSY_WAIT` | constant | 30.0 | BUSY_WAIT = 30.0 # seconds a job waits for a free slot before the app says it is busy | convention | a timeout; earliest in this file: e68ea41 2026-10-12 app_jobs: render jobs with time, memory and concurrency limits; failures remembe… | app_jobs | - |
| `app_jobs:ISOLATED_TIMEOUT` | constant | 90.0 | ISOLATED_TIMEOUT = 90.0 # the timeout of a job given a memory cap and no timeout: an isolated job is always timed | convention | a timeout; earliest in this file: 183b4ff 2026-10-12 app_jobs: isolated jobs always timed, deaths told apart, results bounded, missin… | app_jobs | - |
| `app_jobs:RESULT_MAX` | constant | 100_000_000 | RESULT_MAX = 100_000_000 # bytes a job process may hand back (a giant's 3D page is 20-48 MB) | convention | a size limit; earliest in this file: 183b4ff 2026-10-12 app_jobs: isolated jobs always timed, deaths told apart, results bounded, missin… | app_jobs | - |
| `app_jobs:HOSTED` | constant | Limits(timeout=90, memory_mb=1024, max_residues=30_000, jobs=2) | HOSTED = Limits(timeout=90, memory_mb=1024, max_residues=30_000, jobs=2) # Q-A11; memory: 2 jobs x 1 GB and the | convention | the hosted limits (ruling Q-A11); earliest in this file: e68ea41 2026-10-12 app_jobs: render jobs with time, memory and concurrency limits; failures remembe… | app, app_jobs | - |
| `app_jobs:_DRAWN_MAX` | constant | 4 | _DRAWN_MAX = 4 | convention | a cache size; earliest in this file: e68ea41 2026-10-12 app_jobs: render jobs with time, memory and concurrency limits; failures remembe… | - | - |
| `app_jobs:_FAILED_MAX` | constant | 256 | _FAILED_MAX = 256 | convention | a cache size; earliest in this file: e68ea41 2026-10-12 app_jobs: render jobs with time, memory and concurrency limits; failures remembe… | - | - |
| `app_jobs:_isolated:argument:5@23e7a8e9` | argument | 5 | proc.join(5) | unknown | - | - | - |
| `app_jobs:_isolated:argument:5@23e7a8e9#2` | argument | 5 | proc.join(5) # it has answered and is closing | unknown | - | - | - |
| `app_jobs:_isolated:argument:5@23e7a8e9#3` | argument | 5 | proc.join(5) | unknown | - | - | - |
| `app_jobs:figure_file:argument:8@69a3429d` | argument | 8 | @functools.lru_cache(maxsize=8) # app.CACHES counts these too | convention | a cache size | app, app_jobs | - |
| `app_jobs:data_file:argument:8@69a3429d` | argument | 8 | @functools.lru_cache(maxsize=8) | convention | a cache size | - | - |
| `app_jobs:sequence_file:argument:8@69a3429d` | argument | 8 | @functools.lru_cache(maxsize=8) | convention | a cache size | app, app_jobs | - |

## app_views

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `app_views:figure_card:comparison:3@46537647` | comparison | 3 | elif attrs.get("height") and attrs.get("width", 0) / attrs["height"] > 3: | unknown | - | app, app_uploads | - |
| `app_views:figure_view:argument:2000@7bd8c347` | argument | 2000 | st.caption(shown(result["like_note"], limit=2000, paths=not local_files())) | unknown | - | - | - |
| `app_views:sequence_options:argument:60@8cac53ab` | argument | 60 | "columns": int(ss.get("seq_columns", 60)), | unknown | - | - | - |
| `app_views:sequence_view:argument:20@bd227f19` | argument | 20 | st.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns", width=300, persist_state="session") | unknown | - | - | - |
| `app_views:sequence_view:argument:150@bd227f19` | argument | 150 | st.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns", width=300, persist_state="session") | unknown | - | - | - |
| `app_views:sequence_view:argument:60@bd227f19` | argument | 60 | st.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns", width=300, persist_state="session") | unknown | - | - | - |
| `app_views:sequence_view:argument:5@bd227f19` | argument | 5 | st.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns", width=300, persist_state="session") | unknown | - | - | - |
| `app_views:sequence_view:argument:300@bd227f19` | argument | 300 | st.slider("Columns per row", 20, 150, 60, step=5, key="seq_columns", width=300, persist_state="session") | unknown | - | - | - |
| `app_views:sequence_view:argument:160@06c3de4c` | argument | 160 | st.selectbox("Chain", chains, key="seq_chain", format_func=capped, width=160, persist_state="session") | unknown | - | - | - |
| `app_views:explorer_view:argument:1100@82aecbfc` | argument | 1100 | st.iframe(page, height=1100) # foldmap's own page | unknown | - | - | - |

## cath

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `cath:NO_CHOPPING_MAX_AGE` | constant | 30 * 24 * 3600 | NO_CHOPPING_MAX_AGE = 30 * 24 * 3600 # seconds a cached "no CATH chopping" (404) is trusted: CATH grows | convention | a cache age; earliest in this file: eb9f304 2026-10-16 Domains review minors: strays re-read, CATH fetch edges and chain names, S5 reca… | - | - |
| `cath:chopping:argument:30@f5c94768` | argument | 30 | with urllib.request.urlopen(req, timeout=30) as r: | unknown | - | cath, random100 | - |
| `cath:chopping:comparison:404@a49298a0` | comparison | 404 | if err.code != 404: | unknown | - | cath, random100 | - |

## cli

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `cli:NA_HOLD_MIN` | constant | 3 | NA_HOLD_MIN = 3 # residues of one protomer touching a strand for that protomer to hold it | unknown | earliest in this file: 8007c6a 2026-10-11 C4 fix: two protomers hold a strand only with three residues each | focus, na_layout | - |
| `cli:DYAD_IN_PAGE` | constant | np.radians(10) | DYAD_IN_PAGE = np.radians(10) # a two-fold axis within this of the page plane is drawn as a mirror line | unknown | earliest in this file: 2946ad1 2026-10-11 A dimer whose dyad lies in the page is laid out as one protomer and its mirror i… | - | - |
| `cli:ROSETTE_ORDERS` | constant | (4, 7) | ROSETTE_ORDERS = (4, 7) # cyclic=rosette: Cn assemblies of these orders are drawn seen down their axis. Plan 21's | measured | "Plan 21's look test kept these two: C3, C5, C6 and C8 rosettes ... read worse than their unrolled figures on most entries tried, with more check failures"; earliest in this file: 80a381c 2026-10-20 Rosette view of a cyclic assembly, opt-in (style cyclic=rosette, default unrolle… | sensitivity | - |
| `cli:DYAD_TURNS` | constant | (0, 15, -15, 30, -30, 45, -45, 60, -60, 75, -75) | DYAD_TURNS = (0, 15, -15, 30, -30, 45, -45, 60, -60, 75, -75) # degrees about the dyad, in order of preference | unknown | earliest in this file: 10ea3fa 2026-10-20 A two-fold assembly's dyad view turned to keep its interface together (plan 19 T… | - | - |
| `cli:DYAD_TURN_GAIN` | constant | 0.01 | DYAD_TURN_GAIN = 0.01 # turns this close to the best rank correlation count as tied: the earliest is taken | unknown | earliest in this file: 10ea3fa 2026-10-20 A two-fold assembly's dyad view turned to keep its interface together (plan 19 T… | symmetry | - |
| `cli:AZIMUTH_STEP` | constant | 5 | AZIMUTH_STEP = 5 # degrees: the turns about the membrane normal a membrane view's page-right is chosen from (T5-3) | unknown | cites T5-3; earliest in this file: 102a6f2 2026-10-24 Membrane view refinements from the ComEC look test (plan 26 Task 5, rulings T5-1… | membrane | - |
| `cli:MIN_TAU_GAIN_INV` | constant | 10 | MIN_TAU_GAIN_INV = 10 | measured | "every turn worth less on the membrane campaign grew its page past the campaign's +15% gate" (entries 6JO6, 6ZKT, 30BU); earliest in this file: 102a6f2 2026-10-24 Membrane view refinements from the ComEC look test (plan 26 Task 5, rulings T5-1… | membrane | - |
| `cli:FOCUS_CHAINS` | constant | 6 | FOCUS_CHAINS, FOCUS_ELEMENTS = 6, 160 # symmetric assemblies larger than this are drawn one subunit at a time | unknown | earliest in this file: 11243f3 2026-10-09 Large assemblies: draw one subunit with what its neighbours add to its fold | focus, scope | - |
| `cli:FOCUS_ELEMENTS` | constant | 160 | FOCUS_CHAINS, FOCUS_ELEMENTS = 6, 160 # symmetric assemblies larger than this are drawn one subunit at a time | unknown | earliest in this file: 11243f3 2026-10-09 Large assemblies: draw one subunit with what its neighbours add to its fold | - | - |
| `cli:_DOMAIN_PULL` | constant | 6 | _DOMAIN_PULL = 6 # CA-pair-equivalents of attraction between any two elements of one domain | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `cli:_any_centre:argument:3@f684d89c` | argument | 3 | return np.zeros(3) # nothing to draw (a nucleic acid with --no-dna) | convention | an array's shape | - | - |
| `cli:_within_protomer:argument:-2@a5990100` | argument | -2 | return [b for b in bridges if sector.get(bb.labels[b.i].chain, -1) == sector.get(bb.labels[b.j].chain, -2)] | unknown | - | - | - |
| `cli:_dyad_turn:comparison:3@7abf1d4b` | comparison | 3 | if len(cents) < 3: | unknown | - | symmetry | - |
| `cli:_dyad_turn:argument:9@e1d90e4a` | argument | 9 | score[turn] = round(float((page * in3d).sum() / norm), 9) if norm > 0 else 0.0 | convention | rounding digits | symmetry | - |
| `cli:make_layout:comparison:0.9@a523a9c0` | comparison | 0.9 | along = sym is not None and abs(float(np.dot(sym.axis, mem.normal))) > 0.9 # already stands upright) | unknown | - | 310, assembly, cath, chains +51 | - |
| `cli:make_layout.matches:argument:20@52bac9c2` | argument | 20 | enough = not borrowed or shared >= min(20, 0.5 * min(len(rows[c]), len(names[chain]))) | unknown | - | na_contacts, na_layout, sheets | - |
| `cli:make_layout.matches:comparison:0.9@771a2cb3` | comparison | 0.9 | return shared > 0 and same >= 0.9 * shared and enough | unknown | - | na_contacts, na_layout, sheets | - |
| `cli:_in_sequence:argument:-180@6d2655ae` | argument | -180 | turns = sorted(range(-180 + AZIMUTH_STEP, 181, AZIMUTH_STEP), key=lambda a: (abs(a), a < 0)) | unknown | - | membrane | - |
| `cli:_in_sequence:argument:181@6d2655ae` | argument | 181 | turns = sorted(range(-180 + AZIMUTH_STEP, 181, AZIMUTH_STEP), key=lambda a: (abs(a), a < 0)) | unknown | - | membrane | - |
| `cli:_face_medians.page_y:comparison:1e-9@0f0573ec` | comparison | 1e-9 | elif abs(slope) > 1e-9: | convention | a floating-point tolerance | - | - |
| `cli:_one_of_each:comparison:0.9@92fb00cc` | comparison | 0.9 | if SequenceMatcher(None, seq[g[0]], seq[c], autojunk=False).ratio() > 0.9: | unknown | - | - | - |
| `cli:build_parser:argument:60@29d007de` | argument | 60 | sq.add_argument("--columns", type=int, default=60, help="residues (alignment columns) per row (default 60)") | unknown | - | cli, elstyle, names, science | - |
| `cli:build_parser:argument:0.7@afd84f49` | argument | 0.7 | default=0.7, | unknown | - | cli, elstyle, names, science | - |

## compare

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `compare:PANEL_GAP` | constant | 0.5 | PANEL_GAP = 0.5 # page units between the two panels' pages (each keeps its own pad) | unknown | earliest in this file: 15f0f69 2026-10-21 foldmap compare: A and B drawn like A on one page, one scale, one legend, the no… | - | - |
| `compare:NOTE_SIZE` | constant | 0.85 | NOTE_SIZE = 0.85 # the note's text, as a fraction of the figure's (the context caption's) | unknown | earliest in this file: 15f0f69 2026-10-21 foldmap compare: A and B drawn like A on one page, one scale, one legend, the no… | - | - |
| `compare:NOTE_LINE` | constant | 1.25 | NOTE_LINE = 1.25 # its line step, in its own ems | unknown | earliest in this file: 15f0f69 2026-10-21 foldmap compare: A and B drawn like A on one page, one scale, one legend, the no… | - | - |
| `compare:NOTE_PAD` | constant | 0.3 | NOTE_PAD = 0.3 # page units above and below it | unknown | earliest in this file: 15f0f69 2026-10-21 foldmap compare: A and B drawn like A on one page, one scale, one legend, the no… | - | - |
| `compare:MOVED` | constant | 0.1 | MOVED = 0.1 # Å: B's panel is coloured by displacement when a residue moved this far (coordinates' precision, 1e-3) | unknown | earliest in this file: 6380df1 2026-10-21 The diff drawn (plan 13e-2): foldmap diff A B -o diff.svg [--json] [--domain NAM… | - | - |
| `compare:_page:argument:30@a10d0bb6` | argument | 30 | per_line = max(30, int((width - 2) * scale * 72 / (font * 0.52))) # characters per line, as the context caption | unknown | - | - | - |
| `compare:_page:argument:72@a10d0bb6` | argument | 72 | per_line = max(30, int((width - 2) * scale * 72 / (font * 0.52))) # characters per line, as the context caption | convention | 72 points to the inch | - | - |
| `compare:_page:argument:0.52@a10d0bb6` | argument | 0.52 | per_line = max(30, int((width - 2) * scale * 72 / (font * 0.52))) # characters per line, as the context caption | unknown | - | - | - |
| `compare:_page.at:comparison:1e-9@d4cc3f31` | comparison | 1e-9 | if abs(w_u - w) > 1e-9 or abs(h_u - h) > 1e-9: # drawn as measured (render.page_size) | convention | a floating-point tolerance | - | - |
| `compare:_page.at:comparison:1e-9@d4cc3f31#2` | comparison | 1e-9 | if abs(w_u - w) > 1e-9 or abs(h_u - h) > 1e-9: # drawn as measured (render.page_size) | convention | a floating-point tolerance | - | - |
| `compare:_round:argument:3@2143302c` | argument | 3 | return float(round(float(x), 3)) + 0.0 # to 1e-3, and never -0.0 | convention | rounding digits | - | - |

## criteria

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `criteria:LIMITS` | constant | { "view_timeout_s": 900, # a view's process is killed after this many seconds (A1 then fails) "default_s": 60, # A2: wa… | LIMITS = { | unknown | cites Q-R3; Q-R7; Q-R12; earliest in this file: 55d0e3b 2026-10-17 The figure criteria A1-A20 live in the package (foldmap.criteria), and `foldmap … | random100 | - |
| `criteria:MEMBRANE_SHOWN` | constant | 4 | MEMBRANE_SHOWN = 4 # findings A18's detail shows per part, the rest counted (as the tethers' detail does) | convention | how many findings a detail shows; earliest in this file: b3657c8 2026-10-25 A18 fix round 1: every loop held to its residues, rules shared in membrane.py (p… | membrane | - |
| `criteria:_BLOCK_NAME_EM` | constant | 0.9 | _BLOCK_NAME_EM = 0.9 # render._draw_block draws a stem's name on its block at this many ems of the figure's text | geometry | "render._draw_block draws a stem's name on its block at this many ems" (a copy of render's); earliest in this file: 55d0e3b 2026-10-17 The figure criteria A1-A20 live in the package (foldmap.criteria), and `foldmap … | - | - |
| `criteria:PAIR_IDENTITY` | constant | 0.9 | PAIR_IDENTITY = 0.9 # A19: two entries' first chains are one protein at this identity over their aligned residues | unknown | earliest in this file: 44a1c1b 2026-10-17 F1 fix round 1: equal names always match (A19); a pair of two proteins is not ju… | random100 | - |
| `criteria:PAIR_COVERAGE` | constant | 0.5 | PAIR_COVERAGE = 0.5 # ... with this fraction of the shorter chain's residues identical | unknown | earliest in this file: 44a1c1b 2026-10-17 F1 fix round 1: equal names always match (A19); a pair of two proteins is not ju… | - | - |
| `criteria:LIKE_MIN_MATES` | constant | 3 | LIKE_MIN_MATES = 3 # H2, H3: judged on this many mates drawn at least (as figcheck.COPY_MIN_MATCHED: 2 always fit) | unknown | earliest in this file: a1c3bdf 2026-10-21 Checks H1-H5: a structure drawn like a reference, judged (plan 13a Task 9) | diff | - |
| `criteria:LIKE_APART` | constant | 1.0 | LIKE_APART = 1.0 # H3: two mates more than this many page units apart, in both figures, are compared ... | unknown | earliest in this file: a1c3bdf 2026-10-21 Checks H1-H5: a structure drawn like a reference, judged (plan 13a Task 9) | - | - |
| `criteria:LIKE_AXIS` | constant | 0.5 | LIKE_AXIS = 0.5 # ... on each page axis that parts them by more than this in both figures | unknown | earliest in this file: a1c3bdf 2026-10-21 Checks H1-H5: a structure drawn like a reference, judged (plan 13a Task 9) | - | - |
| `criteria:tether_lengths:argument:3@d48b8088` | argument | 3 | gid, length = f"ligand-tether:{g.name}:{g.chain}{g.seq}:{k}", round(float(np.hypot(tx - cx, ty - cy)), 3) | convention | rounding digits | random100 | - |
| `criteria:dssp_agreement:argument:4@6cb5d594` | argument | 4 | near = {keys[j] for k in bb.linked for j in range(k - 4, k + 5) if 0 <= j < len(keys)} - stood_in | unknown | - | random100 | - |
| `criteria:dssp_agreement:argument:5@6cb5d594` | argument | 5 | near = {keys[j] for k in bb.linked for j in range(k - 4, k + 5) if 0 <= j < len(keys)} - stood_in | unknown | - | random100 | - |
| `criteria:dssp_agreement:argument:4@b45e2c6f` | argument | 4 | "agree": round(1 - len(differ) / len(compared), 4), | convention | rounding digits | random100 | - |
| `criteria:figure_criteria:comparison:4@14003088` | comparison | 4 | more = lambda items: [f"(+{len(items) - 4} more)"] if len(items) > 4 else [] # noqa: E731 | unknown | - | criteria, membrane | - |
| `criteria:one_protein_names:argument:3@2bb1304e` | argument | 3 | return {"one": one, "identity": round(identity, 3), "aligned": len(pairs), "lengths": [len(xs), len(ys)]} | convention | rounding digits | - | - |
| `criteria:pair_check:argument:3@b7342df0` | argument | 3 | if len(pairs) < max(min(3, len(ea), len(eb)), max(len(ea), len(eb)) / 2): | unknown | - | random100 | - |
| `criteria:pair_check:comparison:8@d2340d28` | comparison | 8 | detail += "; " + "; ".join(shown[:8]) + (f" (+{len(shown) - 8} more)" if len(shown) > 8 else "") if shown else "" | unknown | - | random100 | - |
| `criteria:place_deviation:argument:3@4811ff2a` | argument | 3 | "median": round(float(np.median(dev)), 3), | convention | rounding digits | criteria | - |
| `criteria:place_deviation:argument:90@9c0a948c` | argument | 90 | "p90": round(float(np.percentile(dev, 90)), 3), | unknown | - | criteria | - |
| `criteria:place_deviation:argument:3@9c0a948c` | argument | 3 | "p90": round(float(np.percentile(dev, 90)), 3), | convention | rounding digits | criteria | - |
| `criteria:relations:argument:4@5d0d7bbe` | argument | 4 | "share": round(1 - len(turned) / judged, 4), | convention | rounding digits | criteria | - |
| `criteria:sheet_orders:comparison:-1e-9@ec9c5f34` | comparison | -1e-9 | if across[0] < -1e-9 or (abs(across[0]) <= 1e-9 and across[1] < 0): # read left to right (or upward) | convention | a floating-point tolerance | like | - |
| `criteria:sheet_orders:comparison:1e-9@ec9c5f34` | comparison | 1e-9 | if across[0] < -1e-9 or (abs(across[0]) <= 1e-9 and across[1] < 0): # read left to right (or upward) | convention | a floating-point tolerance | like | - |
| `criteria:h1_matched:argument:4@d333af19` | argument | 4 | value = {"tm_score": round(float(aln.tm_score), 4), "same_fold": one, "alike": alike, "refused": refused} | convention | rounding digits | criteria | - |

## datafile

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `datafile:SCHEMA_VERSION` | constant | 1 | SCHEMA_VERSION = 1 | convention | a schema version; earliest in this file: f3d88c4 2026-10-21 foldmap data: elements and sheets as versioned JSON (schema v1, --schema); regis… | datafile | - |
| `datafile:BRACKETS` | constant | (("(", ")"), ("[", "]"), ("{", "}"), ("<", ">"), *((chr(c), chr(c + 32)) for c in range(65, 91))) | BRACKETS = (("(", ")"), ("[", "]"), ("{", "}"), ("<", ">"), *((chr(c), chr(c + 32)) for c in range(65, 91))) | convention | dot-bracket notation; earliest in this file: b0806a6 2026-10-21 foldmap data: the nucleic acid (nucleic_acid section: strands, sequences, the ba… | datafile | - |
| `datafile:_num:argument:3@d1c6dd36` | argument | 3 | return round(v, 3) + 0.0 | convention | rounding digits | - | - |
| `datafile:_angle:comparison:-180.0@79839bed` | comparison | -180.0 | while v <= -180.0: | unknown | - | - | - |
| `datafile:_angle:argument:360.0@95af8aa9` | argument | 360.0 | v = _num(v + 360.0) | unknown | - | - | - |
| `datafile:_angle:comparison:180.0@5e415b76` | comparison | 180.0 | while v > 180.0: | unknown | - | - | - |
| `datafile:_angle:argument:360.0@cdb49665` | argument | 360.0 | v = _num(v - 360.0) | unknown | - | - | - |
| `datafile:zip_bytes:argument:1980@88a63680` | argument | 1980 | info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0)) | unknown | - | datafile | - |
| `datafile:_bridges:argument:-2@64f2b762` | argument | -2 | bridges = [b for b in bridges if sector.get(labels[b.i].chain, -1) == sector.get(labels[b.j].chain, -2)] | unknown | - | - | - |
| `datafile:_drawn_order:argument:9@e1942166` | argument | 9 | if round(float(across[0]), 9) < 0 or (round(float(across[0]), 9) == 0 and across[1] < 0): | convention | rounding digits | - | - |
| `datafile:_drawn_order:argument:9@e1942166#2` | argument | 9 | if round(float(across[0]), 9) < 0 or (round(float(across[0]), 9) == 0 and across[1] < 0): | convention | rounding digits | - | - |
| `datafile:_drawn_order:argument:6@a9ddd69f` | argument | 6 | return [i for _, i in sorted((round(float(np.array([p.cx, p.cy]) @ across), 6), i) for i, p in zip(ids, at))] | convention | rounding digits | - | - |
| `datafile:_letters:argument:26@ab84d5dc` | argument | 26 | k, r = divmod(k - 1, 26) | unknown | - | - | - |
| `datafile:_letters:argument:65@8ce5c395` | argument | 65 | out = chr(65 + r) + out | unknown | - | - | - |
| `datafile:_ligands:argument:5@c5dac063` | argument | 5 | search = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate() | unknown | - | - | - |
| `datafile:_na_faces:argument:5@09d602da` | argument | 5 | ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate() | unknown | - | - | - |

## deep

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `deep:VERSION` | constant | 1 | VERSION = 1 | convention | the check report schema version (deep.py docstring) | deep_check | - |
| `deep:TIMEOUT` | constant | 600 | TIMEOUT = 600 # s for an external program | convention | a time limit for mkdssp or RNApolis on one file, as random100.mkdssp_rows and na_validate use | - | - |
| `deep:parse_dssp:comparison:17@1840e38e` | comparison | 17 | if not started or len(line) < 17 or line[13] == "!": | convention | the DSSP file format: the secondary-structure code is column 17 | deep_check | - |
| `deep:mkdssp_rows:argument:60@6360e5a9` | argument | 60 | done = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60) | convention | a time limit for `mkdssp --version` | deep_check | - |
| `deep:pairs_check:argument:4@f10d1b07` | argument | 4 | "precision": round(precision, 4), | convention | rounding digits | deep_check | - |
| `deep:pairs_check:argument:4@7adbda5e` | argument | 4 | "recall": round(recall, 4), | convention | rounding digits | deep_check | - |
| `deep:pairs_check:argument:4@d0160d6a` | argument | 4 | "f1": round(f1, 4), | convention | rounding digits | deep_check | - |

## diff

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `diff:CORE_CUT` | constant | 1.0 | CORE_CUT = 1.0 # core pairs farther apart than CORE_CUT·d0 after a fit are dropped from the next (_core_fit) | unknown | earliest in this file: b565388 2026-10-21 The diff engine (plan 13e-1): diff.diff_structures(a, chain_a, b, chain_b, match… | - | - |
| `diff:PARTNER_SHIFT` | constant | 10.0 | PARTNER_SHIFT = 10.0 # Å: B's chain is in the place of A's when their centres lie this close after the core fit | unknown | earliest in this file: bb620d1 2026-10-22 foldmap diff open minors (plan 13e): partners paired by place, H5 names its pane… | - | - |
| `diff:PARTNER_SIZE` | constant | 0.5 | PARTNER_SIZE = 0.5 # ... and the shorter of the two has at least this fraction of the longer's residues | unknown | earliest in this file: bb620d1 2026-10-22 foldmap diff open minors (plan 13e): partners paired by place, H5 names its pane… | - | - |
| `diff:Diff:argument:3@ef3af752` | argument | 3 | rotation: np.ndarray = field(default_factory=lambda: np.eye(3), compare=False) # A ≈ R·B + t, the core fit | convention | an array's shape | app, app_jobs, diff | - |
| `diff:Diff:argument:3@9d6bd4af` | argument | 3 | translation: np.ndarray = field(default_factory=lambda: np.zeros(3), compare=False) | convention | an array's shape | app, app_jobs, diff | - |
| `diff:_core_fit:argument:3@07277a6c` | argument | 3 | return (), np.eye(3), np.zeros(3), np.zeros(0) | convention | an array's shape | - | - |
| `diff:_core_fit:argument:3@07277a6c#2` | argument | 3 | return (), np.eye(3), np.zeros(3), np.zeros(0) | convention | an array's shape | - | - |
| `diff:_core_fit:comparison:3@2106b6cb` | comparison | 3 | if core.sum() < 3: | unknown | - | - | - |
| `diff:_core_fit:comparison:3@2106b6cb#2` | comparison | 3 | if core.sum() < 3: # one or two pairs: no rotation to fit | unknown | - | - | - |
| `diff:_core_fit:argument:3@a8c69eaf` | argument | 3 | r, t = np.eye(3), (pa - pb).mean(axis=0) | convention | an array's shape | - | - |
| `diff:_core_fit:comparison:3@43eaecb4` | comparison | 3 | if keep.sum() < 3: | unknown | - | - | - |
| `diff:_core_fit:comparison:3@7a92f438` | comparison | 3 | if near.sum() < 3 or np.array_equal(near, keep): | unknown | - | - | - |
| `diff:_degrees:argument:3@0176be06` | argument | 3 | return round(float(np.degrees(axis_angle(rot)[0])), 3) + 0.0 | convention | rounding digits | - | - |
| `diff:_domain_turns:comparison:3@5889ce77` | comparison | 3 | if sel.sum() < 3: | unknown | - | - | - |

## domain_score

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `domain_score:TOLERANCE` | constant | 20 | TOLERANCE = 20 # residues: a boundary this near the other side's counts as found | unknown | earliest in this file: 6c6c089 2026-10-16 Domain calibration: the 33-chain CATH reference set, overlap and boundary scores… | - | - |
| `domain_score:TARGETS` | constant | {"overlap": 0.80, "f1": 0.70, "recall": 0.85} | TARGETS = {"overlap": 0.80, "f1": 0.70, "recall": 0.85} # plan 12 ruling R2, means over the CATH set | unknown | cites ruling R2; earliest in this file: 6c6c089 2026-10-16 Domain calibration: the 33-chain CATH reference set, overlap and boundary scores… | domain_score, random100 | - |
| `domain_score:fixture_arrays:argument:1000@315d04bb` | argument | 1000 | "ca": np.round(np.asarray(inputs.ca) * 1000).astype(np.int32), | convention | rounding digits | - | - |
| `domain_score:fixture_arrays:argument:3@6eff76dc` | argument | 3 | "contacts": np.array(contacts, dtype=np.int32).reshape(-1, 3), | convention | an array's shape | - | - |
| `domain_score:load_fixture:argument:3@b62bc6e6` | argument | 3 | zero = np.zeros(3) | convention | an array's shape | domain_score, membrane | - |
| `domain_score:load_fixture:argument:1000@774c6d43` | argument | 1000 | out[key] = ChainInputs(labels, a["ca"].astype(float) / 1000, sses, contacts, rows) | unknown | - | domain_score, membrane | - |
| `domain_score:held_out:comparison:3@7fef68c2` | comparison | 3 | if len(held) >= len(entries) // 3: | unknown | - | domain_score | - |

## dssp

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `dssp:_Q1Q2F` | constant | 0.084 * 332.0 | _Q1Q2F = 0.084 * 332.0 # kcal/mol * Å | biology | DSSP's electrostatic constant (Kabsch & Sander 1983); earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | dssp | - |
| `dssp:_HBOND_E` | constant | -0.5 | _HBOND_E = -0.5 | biology | DSSP's H-bond energy cutoff (Kabsch & Sander 1983); earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | dssp | - |
| `dssp:_HBOND_E_MIN` | constant | -9.9 | _HBOND_E_MIN = -9.9 # kcal/mol: DSSP's floor, also the energy of any pair with atoms closer than _MIN_DIST | biology | "DSSP's floor" (as mkdssp); earliest in this file: 4aeff2a 2026-10-09 DSSP: 2.5 Å peptide-bond cutoff and the -9.9 kcal/mol energy floor, as mkdssp | - | - |
| `dssp:_MIN_DIST` | constant | 0.5 | _MIN_DIST = 0.5 # Å | biology | DSSP's closest-atom distance (as mkdssp); earliest in this file: 4aeff2a 2026-10-09 DSSP: 2.5 Å peptide-bond cutoff and the -9.9 kcal/mol energy floor, as mkdssp | - | - |
| `dssp:_CA_CUTOFF` | constant | 9.0 | _CA_CUTOFF = 9.0 | biology | DSSP's CA search distance (as mkdssp); earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `dssp:_hydrogen_positions:argument:3@76229b13` | argument | 3 | h = np.full((n, 3), np.nan) | convention | an array's shape | - | - |
| `dssp:_hbond_energy:argument:1000@3d71898f` | argument | 1000 | e = np.round(e * 1000) / 1000 | convention | rounding digits | dssp | - |
| `dssp:_find_bridges:comparison:3@ec839397` | comparison | 3 | if j < i + 3 or not (prev[i] and nxt[i] and prev[j] and nxt[j]): | unknown | - | dssp | - |
| `dssp:_ladders:comparison:6@123fa1eb` | comparison | 6 | if ibj >= iei + 6: # sorted by first i: no later ladder can join this one either | unknown | - | dssp | - |
| `dssp:_ladders:comparison:6@30f16b0e` | comparison | 6 | and _gap(ibj, iei) < 6 | unknown | - | dssp | - |
| `dssp:_ladders:comparison:6@3d60ab27` | comparison | 6 | if joinable and ((other < 6 and _gap(ibj, iei) < 3) or other < 3): | unknown | - | dssp | - |
| `dssp:_ladders:comparison:3@3d60ab27` | comparison | 3 | if joinable and ((other < 6 and _gap(ibj, iei) < 3) or other < 3): | unknown | - | dssp | - |
| `dssp:_ladders:comparison:3@3d60ab27#2` | comparison | 3 | if joinable and ((other < 6 and _gap(ibj, iei) < 3) or other < 3): | unknown | - | dssp | - |

## elstyle

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `elstyle:WIDTH_MIN` | constant | 0.25 | WIDTH_MIN, WIDTH_MAX = 0.25, 4.0 # a width factor's range: thinner than a quarter vanishes, wider than 4 swamps | unknown | earliest in this file: 2fd5f74 2026-10-22 Each element's own colour and width (plan 23 Task 1): --color/--width, element_c… | - | - |
| `elstyle:WIDTH_MAX` | constant | 4.0 | WIDTH_MIN, WIDTH_MAX = 0.25, 4.0 # a width factor's range: thinner than a quarter vanishes, wider than 4 swamps | unknown | earliest in this file: 2fd5f74 2026-10-22 Each element's own colour and width (plan 23 Task 1): --color/--width, element_c… | - | - |

## features

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `features:CONTACT_CUTOFF` | constant | 8.0 | CONTACT_CUTOFF = 8.0 # Å between CA atoms | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `features:MIN_CONTACT_PAIRS` | constant | 3 | MIN_CONTACT_PAIRS = 3 # CA pairs needed to call two elements neighbours | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | features | - |
| `features:BUNDLE_MIN` | constant | 10 | BUNDLE_MIN = 10 # residues: shorter helices are turns at an interface, not a coiled-coil stalk | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `features:BUNDLE_COS` | constant | 0.75 | BUNDLE_COS = 0.75 # \|cos\| between axes: parallel or antiparallel | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `features:BUNDLE_DIST` | constant | 12.0 | BUNDLE_DIST = 12.0 # Å between axis segments | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout | - |
| `features:BUNDLE_COVER` | constant | 0.8 | BUNDLE_COVER = 0.8 # fraction of the shorter helix that must run that close | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `features:KINK_GAP` | constant | 2 | KINK_GAP = 2 # residues: helix pieces this close in one chain, roughly collinear, count as one helix | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `features:KINK_COS` | constant | 0.5 | KINK_COS = 0.5 # a kink bends the axis by up to 60 degrees | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout | - |
| `features:KINK_SHORT` | constant | 7 | KINK_SHORT = 7 # residues: a piece this short has no reliable axis of its own, so it continues its neighbour | unknown | earliest in this file: a74e88e 2026-10-09 Bundles: keep short kinked helix pieces with their helix; compare whole kinked h… | layout | - |
| `features:ON_SHEET_MIN` | constant | 10 | ON_SHEET_MIN = 10 # residues: a shorter helix is a turn at the sheet's edge, not one packed on its face | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:ON_SHEET_STRANDS` | constant | 3 | ON_SHEET_STRANDS = 3 # strands: a smaller sheet has no face to pack on | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:ON_SHEET_REACH` | constant | 10.0 | ON_SHEET_REACH = 10.0 # Å between a helix CA and the nearest sheet CA, counting it as touching the sheet | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:ON_SHEET_TOUCH` | constant | 0.3 | ON_SHEET_TOUCH = 0.3 # fraction of the helix's CAs that must touch the sheet | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:ON_SHEET_FLAT` | constant | 0.5 | ON_SHEET_FLAT = 0.5 # \|cos\| between the helix axis and the sheet normal: the helix lies along the face | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:ON_SHEET_HEIGHT` | constant | (6.0, 14.0) | ON_SHEET_HEIGHT = (6.0, 14.0) # Å from the sheet plane to the helix centroid: one helix packed on the face | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:ON_SHEET_EDGE` | constant | 2.0 | ON_SHEET_EDGE = 2.0 # Å: the centroid, projected on the sheet plane, lies within the sheet's extent plus this | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `features:DOMAIN_MIN_SSE_RESIDUES` | constant | 30 | DOMAIN_MIN_SSE_RESIDUES = 30 # residues in helices/strands an automatic domain needs | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | domains | - |
| `features:DOMAIN_CONTACT` | constant | 12.0 | DOMAIN_CONTACT = 12.0 # Å between CA atoms: two residues in contact, counted for the split | measured | "The residue-level split (residue_domains), calibrated on 22 of the 33 CATH chains" (plan 12 Task 3); earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | - | - |
| `features:DOMAIN_SEPARATION` | constant | 6 | DOMAIN_SEPARATION = 6 # residues apart in sequence at least: nearer pairs are the chain, not the fold | measured | in the block "calibrated on 22 of the 33 CATH chains" (plan 12 Task 3); earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | - | - |
| `features:DOMAIN_NCUT` | constant | 0.37 | DOMAIN_NCUT = 0.37 # a part splits when its best cut's normalized cut is below this | measured | in the block "calibrated on 22 of the 33 CATH chains" (plan 12 Task 3); earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | - | - |
| `features:DOMAIN_MIN_RESIDUES` | constant | 40 | DOMAIN_MIN_RESIDUES = 40 # residues a domain needs | measured | in the block "calibrated on 22 of the 33 CATH chains" (plan 12 Task 3); earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | - | - |
| `features:DOMAIN_MIN_SEGMENT` | constant | 20 | DOMAIN_MIN_SEGMENT = 20 # residues each piece of a cut needs (a domain may be two pieces: 5JDN's alternate) | measured | in the block "calibrated on 22 of the 33 CATH chains" (plan 12 Task 3); "5JDN's alternate"; earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | - | - |
| `features:DOMAIN_MAX_POINTS` | constant | 800 | DOMAIN_MAX_POINTS = 800 # longer chains are split on blocks of consecutive residues, this many blocks at most | convention | a work bound: "longer chains are split on blocks of consecutive residues, this many blocks at most"; earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | - | - |
| `features:DOMAIN_STRAY` | constant | 20 | DOMAIN_STRAY = 20 # residues of elements: a run of one domain's elements this small, between two runs of another | unknown | earliest in this file: c7e1b70 2026-10-16 Auto domains from residue contacts: two-run normalized cuts, strays taken home (… | domain_score | - |
| `features:_segment:default:n` | default | 20 | def _segment(s: SSE, n: int = 20) -> np.ndarray: | unknown | - | - | - |
| `features:_best_cut.ncut:argument:9@fdc98c6f` | argument | 9 | return np.where(ok & (volume > 0) & (total - volume > 0), np.round(value, 9), np.inf) | convention | rounding digits | - | - |

## fetch

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `fetch:AFDB_VERSION` | constant | 6 | AFDB_VERSION = 6 # AlphaFold DB model version to fetch | convention | the AlphaFold DB's model version; earliest in this file: 6a76046 2026-10-09 Add MCP server; fetch structures by ID everywhere | fetch | - |
| `fetch:UPLOAD_MAX_AGE` | constant | 24 * 3600 | UPLOAD_MAX_AGE = 24 * 3600 # seconds an uploaded file is kept in the cache | convention | a cache age; earliest in this file: 95c80bb 2026-10-10 App: FOLDMAP_APP_LOCAL_FILES=0 refuses file paths; uploads kept by hash for a day | fetch | - |
| `fetch:MB` | constant | 1024 * 1024 | MB = 1024 * 1024 | convention | bytes; earliest in this file: 6a76046 2026-10-09 Add MCP server; fetch structures by ID everywhere | - | - |
| `fetch:DOWNLOAD_MAX` | constant | 100 * MB | DOWNLOAD_MAX = 100 * MB # bytes the app downloads at most (fetch max_bytes; 3J3Q, the largest entry, is 242 MB) | convention | a download limit (3J3Q, the largest entry, is 242 MB); earliest in this file: 0e44066 2026-10-12 Downloads stream to disk with a 100 MB cap; the download cache is pruned (a week… | fetch | - |
| `fetch:DOWNLOAD_MAX_AGE` | constant | 7 * 24 * 3600 | DOWNLOAD_MAX_AGE = 7 * 24 * 3600 # seconds a downloaded structure is kept unused in the cache | convention | a cache age; earliest in this file: 0e44066 2026-10-12 Downloads stream to disk with a 100 MB cap; the download cache is pruned (a week… | - | - |
| `fetch:DOWNLOAD_CACHE_MAX` | constant | 2 * 1024 * MB | DOWNLOAD_CACHE_MAX = 2 * 1024 * MB # bytes of downloaded structures kept; the oldest go first | convention | a cache size; earliest in this file: 0e44066 2026-10-12 Downloads stream to disk with a 100 MB cap; the download cache is pruned (a week… | fetch | - |
| `fetch:PRUNE_EVERY` | constant | 3600 | PRUNE_EVERY = 3600 # seconds: the cache is pruned at most this often, from the download path | convention | a cache interval; earliest in this file: 0e44066 2026-10-12 Downloads stream to disk with a 100 MB cap; the download cache is pruned (a week… | - | - |
| `fetch:_PRUNED` | constant | [0.0] | _PRUNED = [0.0] # when the download cache was last pruned (time.monotonic), process-wide | convention | a timestamp holder; earliest in this file: 0e44066 2026-10-12 Downloads stream to disk with a 100 MB cap; the download cache is pruned (a week… | fetch | - |
| `fetch:_download:argument:64@17f6f24e` | argument | 64 | f.seek(max(0, path.stat().st_size - 64)) | unknown | - | - | - |
| `fetch:_stream:argument:30@e9cb25f3` | argument | 30 | with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": AGENT}), timeout=30) as r: | unknown | - | - | - |

## figcheck

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `figcheck:TEXT_OVERLAP_PX2` | constant | 1.0 | TEXT_OVERLAP_PX2 = 1.0 # text over text or a body: overlaps up to this many square pixels are anti-aliasing | unknown | earliest in this file: e690730 2026-10-12 figcheck: page checks for any figure (overlaps, text, scope, copies, domains, st… | - | - |
| `figcheck:TEXT_DEPTH_PX` | constant | 2.0 | TEXT_DEPTH_PX = 2.0 # ... and a text or body reaching no deeper than this into a text's ink is a near-touch (R13): | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | figcheck | - |
| `figcheck:OFF_PAGE_PX` | constant | 1.0 | OFF_PAGE_PX = 1.0 # text may reach this far past the page edge | unknown | earliest in this file: e690730 2026-10-12 figcheck: page checks for any figure (overlaps, text, scope, copies, domains, st… | - | - |
| `figcheck:PANEL_GRID` | constant | 0.05 | PANEL_GRID = 0.05 # page units: the raster domain panel outlines are measured on | unknown | a raster's cell size ("measured on" names the raster, not a measurement); earliest in this file: e690730 2026-10-12 figcheck: page checks for any figure (overlaps, text, scope, copies, domains, st… | - | - |
| `figcheck:_STEP_PX` | constant | 0.5 | _STEP_PX = 0.5 # sampling step (pixels) along a line inside a text box, and across a text box over a body | unknown | earliest in this file: e690730 2026-10-12 figcheck: page checks for any figure (overlaps, text, scope, copies, domains, st… | - | - |
| `figcheck:COLUMN_TOL` | constant | 1e-3 | COLUMN_TOL = 1e-3 # page units (in strand pitches across a sheet, _across): two strands this close share a column | convention | a float tolerance (1e-3 page units); earliest in this file: 7acb510 2026-10-21 D2 rosette, opt-in (style cyclic=rosette): a D2 drawn as two rings of two seen d… | rosette | - |
| `figcheck:ROSETTE_TURN_TOL` | constant | 1e-3 | ROSETTE_TURN_TOL = 1e-3 # page units: a turned copy back within this of its mate stands exactly where its mate does | convention | a float tolerance (1e-3 page units); earliest in this file: 3eb60d0 2026-10-20 Rosette routes, labels, ligands and checks (cyclic=rosette, opt-in): protomer 0'… | - | - |
| `figcheck:STEM_GAP_FLOOR` | constant | 0.08 | STEM_GAP_FLOOR = 0.08 # page units | measured | "the same figures measured with matplotlib 3.10.9 ... and 3.11.2 ... stem-name gaps by up to 0.047 page units (1KX5; 0.015 elsewhere)"; earliest in this file: 4c56160 2026-10-12 figcheck: ink measured under fixed hinting, one draw; every check pinned by a te… | figcheck | - |
| `figcheck:LINE_TEXT_CLEARANCE_PX` | constant | 2.0 | LINE_TEXT_CLEARANCE_PX = 2.0 | unknown | earliest in this file: 4c56160 2026-10-12 figcheck: ink measured under fixed hinting, one draw; every check pinned by a te… | figcheck, ligands, links, route +1 | - |
| `figcheck:COPY_MIN_MATCHED` | constant | 3 | COPY_MIN_MATCHED = 3 # a tied copy is compared on at least this many elements, and on half the smaller copy's | unknown | earliest in this file: 4c56160 2026-10-12 figcheck: ink measured under fixed hinting, one draw; every check pinned by a te… | figcheck | - |
| `figcheck:COPY_IDENTITY` | constant | 0.9 | COPY_IDENTITY = 0.9 # two chains are copies of one protein when this fraction of their aligned residues match ... | unknown | introduced by f2d7043 (C18: 1RM5, 2HKI, 2MJW), whose message does not tie this value to those entries; earliest in this file: f2d7043 2026-10-13 C18 (random-100): A19 and copy matching by sequence alignment; untied copies che… | - | - |
| `figcheck:COPY_LENGTH` | constant | 0.8 | COPY_LENGTH = 0.8 # ... and the shorter is at least this fraction of the longer | unknown | introduced by 2a8ae8e (C30: 3OWW), whose message ties 0.8 to "a bound peptide", no entry; earliest in this file: 2a8ae8e 2026-10-14 C30 (random-100 fix round 1): a peptide is not a copy of its receptor; chains th… | - | - |
| `figcheck:COPY_ALIGNED` | constant | 20 | COPY_ALIGNED = 20 # ... and at least this many do (as model.homolog_chains asks of residues at the same numbers) | geometry | "as model.homolog_chains asks of residues at the same numbers" (model.py: at least 20); earliest in this file: f2d7043 2026-10-13 C18 (random-100): A19 and copy matching by sequence alignment; untied copies che… | - | - |
| `figcheck:honest_crossings:argument:4@ac61fc22` | argument | 4 | rects = np.array([e.rect for _, e in placed], float).reshape(-1, 4) | convention | an array's shape | figcheck, layout, perf | - |
| `figcheck:local_layout.back:comparison:0.99@b954f5b1` | comparison | 0.99 | and math.cos(q.angle - m.angle) > 0.99 | unknown | - | datafile, diff, domains, figcheck +16 | - |
| `figcheck:_bodies_meet:default:touch` | default | 1e-6 | def _bodies_meet(p, q, touch: float = 1e-6) -> bool: | convention | a floating-point tolerance | - | - |
| `figcheck:_rects_meet:comparison:1e-6@8863f609` | comparison | 1e-6 | return min(a[2], b[2]) - max(a[0], b[0]) > 1e-6 and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6 | convention | a floating-point tolerance | - | - |
| `figcheck:_rects_meet:comparison:1e-6@8863f609#2` | comparison | 1e-6 | return min(a[2], b[2]) - max(a[0], b[0]) > 1e-6 and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6 | convention | a floating-point tolerance | - | - |
| `figcheck:<module>:comparison:3@3f5ea4b7` | comparison | 3 | if _MPL < (3, 11): # 3.11 deprecates the factor and hints at 1 always | unknown | - | - | - |
| `figcheck:<module>:comparison:11@3f5ea4b7` | comparison | 11 | if _MPL < (3, 11): # 3.11 deprecates the factor and hints at 1 always | unknown | - | - | - |
| `figcheck:stem_gap_tolerance:argument:1e-3@826e3658` | argument | 1e-3 | return max(margin + 1e-3, STEM_GAP_FLOOR) | unknown | - | figcheck, na_layout, random100 | - |
| `figcheck:Page.__init__:argument:72@e57090b2` | argument | 72 | self.lines.append((gid, _polylines(a), a.get_zorder(), a.get_linewidth() * fig.dpi / 72 / 2)) | convention | 72 points to the inch | links, membrane, random100 | - |
| `figcheck:Page.pair_boxes:comparison:1e-6@ae29b663` | comparison | 1e-6 | if abs(turn - tb.get_rotation()) > 1e-6 or turn % 90 == 0: | convention | a floating-point tolerance | figcheck, na_layout, perf | - |
| `figcheck:Page.pair_boxes:comparison:90@ae29b663` | comparison | 90 | if abs(turn - tb.get_rotation()) > 1e-6 or turn % 90 == 0: | unknown | - | figcheck, na_layout, perf | - |
| `figcheck:Page._stem_bodies:argument:0.1@3a5ceb56` | argument | 0.1 | edges = np.linspace(lo, hi, max(1, math.ceil((hi - lo) / max(c.band, 0.1))) + 1) | unknown | - | - | - |
| `figcheck:_opaque:comparison:0.99@caeed06d` | comparison | 0.99 | return a.get_facecolor()[3] >= 0.99 and (alpha is None or alpha >= 0.99) | unknown | - | - | - |
| `figcheck:_opaque:comparison:0.99@caeed06d#2` | comparison | 0.99 | return a.get_facecolor()[3] >= 0.99 and (alpha is None or alpha >= 0.99) | unknown | - | - | - |
| `figcheck:_clip:comparison:1e-12@8170f50e` | comparison | 1e-12 | if abs(dv) < 1e-12: | convention | a floating-point tolerance | - | - |
| `figcheck:_extents:argument:4@21b58bb3` | argument | 4 | e = np.array([b.extents for b in boxes], float).reshape(-1, 4) | convention | an array's shape | - | - |
| `figcheck:_depth:argument:3@ea2168fa` | argument | 3 | return min(TEXT_DEPTH_PX, *(min(b.width, b.height) / 3 for b in boxes)) | unknown | - | - | - |
| `figcheck:_cells:argument:32@f0a6dfb7` | argument | 32 | return np.sort((gi[hit].astype(np.int64) << 32) + (gj[hit].astype(np.int64) + (1 << 31))) | unknown | - | - | - |
| `figcheck:_cells:argument:31@f0a6dfb7` | argument | 31 | return np.sort((gi[hit].astype(np.int64) << 32) + (gj[hit].astype(np.int64) + (1 << 31))) | unknown | - | - | - |
| `figcheck:domain_checks:argument:4@49dedcf5` | argument | 4 | overlaps.append((na, nb, round(f, 4))) | convention | rounding digits | criteria, domains, figcheck | - |
| `figcheck:_segment_box_gap:comparison:1e-24@384de5ac` | comparison | 1e-24 | if n < 1e-24: | convention | a floating-point tolerance | figcheck | - |
| `figcheck:scope_honest:comparison:4@e572b229` | comparison | 4 | if len(run) >= 4: | unknown | - | figcheck, scope | - |
| `figcheck:order_checks:argument:-2@75c81070` | argument | -2 | bridges = [b for b in bridges if sector.get(chain[b.i].chain, -1) == sector.get(chain[b.j].chain, -2)] | unknown | - | copies, figcheck, rosette, sheets | - |
| `figcheck:drawn_columns:argument:3@849de405` | argument | 3 | cols.setdefault(round(float(np.array([p.cx, p.cy]) @ across), 3), []).append((i, ghost)) | convention | rounding digits | science | - |
| `figcheck:_angstrom:comparison:0.05@e7f1471b` | comparison | 0.05 | return "<0.1" if x < 0.05 else f"{x:.1f}" | unknown | - | membrane | - |

## frame

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `frame:CylinderFrame.basis:comparison:1e-6@f5556a78` | comparison | 1e-6 | out = rho / np.linalg.norm(rho) if np.linalg.norm(rho) > 1e-6 else self.ref | convention | a floating-point tolerance | na_validation, nucleic, symmetry | - |
| `frame:rosette_frame:argument:3@862d3b9c` | argument | 3 | first = np.mean(own, axis=0) - sym.centre if own else np.zeros(3) | convention | an array's shape | - | - |
| `frame:rosette_frame:comparison:1e-6@c443ec78` | comparison | 1e-6 | ref = _snap(_unit(ref)) if np.linalg.norm(ref) > 1e-6 else _any_perpendicular(axis) | convention | a floating-point tolerance | - | - |
| `frame:symmetric_frame:argument:3@862d3b9c` | argument | 3 | first = np.mean(own, axis=0) - sym.centre if own else np.zeros(3) | convention | an array's shape | robustness, symmetry | - |
| `frame:symmetric_frame:comparison:1e-6@b3aafe59` | comparison | 1e-6 | ref = ref / np.linalg.norm(ref) if np.linalg.norm(ref) > 1e-6 else _any_perpendicular(sym.axis) | convention | a floating-point tolerance | robustness, symmetry | - |
| `frame:symmetric_frame:argument:5.0@c64fa255` | argument | 5.0 | max(radius, 5.0), | unknown | - | robustness, symmetry | - |
| `frame:_any_perpendicular:argument:3@c02aeb00` | argument | 3 | e = np.eye(3)[int(np.argmin(np.abs(v)))] | convention | an array's shape | - | - |
| `frame:view_frame:argument:3@b7e1ab47` | argument | 3 | return Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0])) | convention | an array's shape | cli, frame, layout, like +7 | - |
| `frame:view_frame:comparison:3@9d28ae63` | comparison | 3 | if up.shape != (3,) or not np.isfinite(up).all() or np.linalg.norm(up) < 1e-9: | unknown | - | cli, frame, layout, like +7 | - |
| `frame:view_frame:comparison:1e-9@9d28ae63` | comparison | 1e-9 | if up.shape != (3,) or not np.isfinite(up).all() or np.linalg.norm(up) < 1e-9: | convention | a floating-point tolerance | cli, frame, layout, like +7 | - |
| `frame:view_frame:comparison:3@5d30a4f0` | comparison | 3 | if w.shape != (3,) or not np.isfinite(w).all() or np.linalg.norm(w) < 1e-9: | unknown | - | cli, frame, layout, like +7 | - |
| `frame:view_frame:comparison:1e-9@5d30a4f0` | comparison | 1e-9 | if w.shape != (3,) or not np.isfinite(w).all() or np.linalg.norm(w) < 1e-9: | convention | a floating-point tolerance | cli, frame, layout, like +7 | - |
| `frame:view_frame:comparison:1e-6@ca2bed81` | comparison | 1e-6 | v = _unit(v) if np.linalg.norm(v) > 1e-6 else _any_perpendicular(w) | convention | a floating-point tolerance | cli, frame, layout, like +7 | - |
| `frame:view_frame:comparison:1e-9@ff535bef` | comparison | 1e-9 | if np.linalg.norm(flat) < 1e-9: | convention | a floating-point tolerance | cli, frame, layout, like +7 | - |
| `frame:dna_view_frame:comparison:1e-6@b5167068` | comparison | 1e-6 | v = _unit(v) if np.linalg.norm(v) > 1e-6 else _any_perpendicular(u) | convention | a floating-point tolerance | nucleic | - |
| `frame:RadialFrame._outward:comparison:1e-6@950142ad` | comparison | 1e-6 | return rho / np.linalg.norm(rho) if np.linalg.norm(rho) > 1e-6 else self.down | convention | a floating-point tolerance | - | - |
| `frame:dna_radial_frame:comparison:1e-6@3c622536` | comparison | 1e-6 | down = down / np.linalg.norm(down) if np.linalg.norm(down) > 1e-6 else _any_perpendicular(axis) | convention | a floating-point tolerance | nucleic | - |
| `frame:<module>:fallback-marker` | fallback-marker | 1 | down: np.ndarray # fallback outward direction on the axis itself (toward the protein) | unknown | not a fallback: `down`, the outward direction for a point on the axis itself (a degenerate case), no worse result | - | - |

## hbonds

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `hbonds:MIN_DIST` | constant | 2.5 | MIN_DIST = 2.5 # Å, donor-acceptor | biology | donor-acceptor distance, heavy-atom H-bond (module docstring: 2.5-3.5 Å); earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:MAX_DIST` | constant | 3.5 | MAX_DIST = 3.5 # sulfur H-bonds (Cys SG, Met SD) usually run 3.3-3.9 Å; capped here at 3.5 Å by design | biology | "sulfur H-bonds (Cys SG, Met SD) usually run 3.3-3.9 Å; capped here at 3.5 Å by design"; earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:SP2_ANGLE` | constant | 90.0 | SP2_ANGLE = 90.0 # degrees, antecedent-donor-acceptor, for a donor whose H lies in a plane | biology | antecedent-donor-acceptor angle of an sp2 donor (module docstring); earliest in this file: d540a8f 2026-10-10 fix: hbonds roles follow each residue's own bonds (modified nucleotides, termina… | - | - |
| `hbonds:ROTATING_ANGLE` | constant | 60.0 | ROTATING_ANGLE = 60.0 # ... for a donor whose H rotates about its bond to the antecedent | biology | antecedent-donor-acceptor angle of a rotating donor (module docstring); earliest in this file: d540a8f 2026-10-10 fix: hbonds roles follow each residue's own bonds (modified nucleotides, termina… | - | - |
| `hbonds:_BOND` | constant | 1.7 | _BOND = 1.7 # Å, heavy atoms bonded | biology | "heavy atoms bonded" (a covalent bond length); earliest in this file: d540a8f 2026-10-10 fix: hbonds roles follow each residue's own bonds (modified nucleotides, termina… | - | - |
| `hbonds:_PEPTIDE_BOND` | constant | 2.5 | _PEPTIDE_BOND = 2.5 # Å, C(i-1)-N(i), as io | biology | "C(i-1)-N(i)" peptide bond, as io; earliest in this file: d540a8f 2026-10-10 fix: hbonds roles follow each residue's own bonds (modified nucleotides, termina… | - | - |
| `hbonds:_N_TERMINAL` | constant | _Site(D, ("CA",), 3, True) | _N_TERMINAL = _Site(D, ("CA",), 3, True) | biology | donor table (chemistry); earliest in this file: d540a8f 2026-10-10 fix: hbonds roles follow each residue's own bonds (modified nucleotides, termina… | - | - |
| `hbonds:_SIDE_CHAIN` | constant | { "ARG": {"NE": _Site(D, ("CD", "CZ")), "NH1": _Site(D, ("CZ",), 2), "NH2": _Site(D, ("CZ",), 2)}, "LYS": {"NZ": _Site(… | _SIDE_CHAIN = { | biology | donor and acceptor table (chemistry); earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:_HYDROXYL` | constant | {"O5'": (_Site(DA, ("C5'",), 1, True), "P"), "O3'": (_Site(DA, ("C3'",), 1, True), "P")} | _HYDROXYL = {"O5'": (_Site(DA, ("C5'",), 1, True), "P"), "O3'": (_Site(DA, ("C3'",), 1, True), "P")} | biology | donor and acceptor table (chemistry); earliest in this file: d540a8f 2026-10-10 fix: hbonds roles follow each residue's own bonds (modified nucleotides, termina… | - | - |
| `hbonds:_GUANINE` | constant | {"N1": _Site(D, ("C2", "C6")), "N2": _Site(D, ("C2",), 2), "O6": _ACC, "N3": _ACC, "N7": _ACC} | _GUANINE = {"N1": _Site(D, ("C2", "C6")), "N2": _Site(D, ("C2",), 2), "O6": _ACC, "N3": _ACC, "N7": _ACC} | biology | donor and acceptor table (chemistry); earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:_ADENINE` | constant | {"N6": _Site(D, ("C6",), 2), "N1": _ACC, "N3": _ACC, "N7": _ACC} | _ADENINE = {"N6": _Site(D, ("C6",), 2), "N1": _ACC, "N3": _ACC, "N7": _ACC} | biology | donor and acceptor table (chemistry); earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:_CYTOSINE` | constant | {"N4": _Site(D, ("C4",), 2), "N3": _ACC, "O2": _ACC} | _CYTOSINE = {"N4": _Site(D, ("C4",), 2), "N3": _ACC, "O2": _ACC} | biology | donor and acceptor table (chemistry); earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:PROTEIN` | constant | 0 | PROTEIN, NUCLEIC = 0, 1 | convention | an enumeration; earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:NUCLEIC` | constant | 1 | PROTEIN, NUCLEIC = 0, 1 | convention | an enumeration; earliest in this file: 0158d07 2026-10-10 feat: hydrogen-bond detection (foldmap.hbonds): donor/acceptor tables, 2.5-3.5 Å… | - | - |
| `hbonds:_roles:comparison:3@969e2017` | comparison | 3 | if acceptor and name.startswith("N") and len(neighbours) >= 3: | unknown | - | hbonds | - |
| `hbonds:_sites:argument:3@08750f1a` | argument | 3 | np.array([x[0] for x in rows], float).reshape(-1, 3), | convention | an array's shape | hbonds | - |
| `hbonds:_sites:argument:3@6cf3230c` | argument | 3 | np.array([x[1] for x in rows], float).reshape(-1, 3), | convention | an array's shape | hbonds | - |
| `hbonds:_donates:comparison:1e-12@f147b953` | comparison | 1e-12 | angle_ok = np.isnan(cos) \| (cos <= sites.max_cos[d] + 1e-12) | convention | a floating-point tolerance | - | - |

## interactions

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `interactions:SALT_MAX` | constant | 4.0 | SALT_MAX = 4.0 # Å, anionic O to cationic N | biology | module docstring: Barlow & Thornton 1983; earliest in this file: d68bd78 2026-10-21 interactions.py: salt bridges and cation-pi from the heavy atoms (plan 13c-1) | interactions | - |
| `interactions:CATION_PI_MAX` | constant | 6.0 | CATION_PI_MAX = 6.0 # Å, cation to ring centroid | biology | module docstring: PLIP's PICATION_DIST_MAX; Gallivan & Dougherty 1999; earliest in this file: d68bd78 2026-10-21 interactions.py: salt bridges and cation-pi from the heavy atoms (plan 13c-1) | interactions | - |
| `interactions:CATION_PI_OFFSET` | constant | 2.0 | CATION_PI_OFFSET = 2.0 # Å, ring centroid to the cation's projection on the ring plane | biology | module docstring: PLIP's PISTACK_OFFSET_MAX; earliest in this file: d68bd78 2026-10-21 interactions.py: salt bridges and cation-pi from the heavy atoms (plan 13c-1) | interactions | - |
| `interactions:STACK_MAX` | constant | 5.5 | STACK_MAX = 5.5 # Å, ring centroid to ring centroid | biology | module docstring: PLIP's PISTACK_DIST_MAX; earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | interactions | - |
| `interactions:STACK_PARALLEL` | constant | 30.0 | STACK_PARALLEL = 30.0 # degrees between the ring normals, at most, for a parallel stack | biology | module docstring: PLIP's PISTACK_ANG_DEV; earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | interactions | - |
| `interactions:STACK_T_MIN` | constant | 60.0 | STACK_T_MIN = 60.0 # ... at least, for a T-shaped one (up to 90) | biology | module docstring: PLIP's PISTACK_ANG_DEV (T-shaped 60-90°); earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | interactions | - |
| `interactions:STACK_OFFSET` | constant | 2.0 | STACK_OFFSET = 2.0 # Å, centroid to the other centroid's projection on its ring plane (the smaller of the two) | biology | module docstring: PLIP's PISTACK_OFFSET_MAX; earliest in this file: d68bd78 2026-10-21 interactions.py: salt bridges and cation-pi from the heavy atoms (plan 13c-1) | interactions | - |
| `interactions:HYDROPHOBIC_MAX` | constant | 4.0 | HYDROPHOBIC_MAX = 4.0 # Å, carbon to carbon | biology | module docstring: PLIP's HYDROPH_DIST_MAX; earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | interactions | - |
| `interactions:_COVALENT` | constant | 2.0 | _COVALENT = 2.0 # Å, two heavy atoms of one residue this close are bonded | biology | module docstring: "C-S is 1.8 Å, C-Se 1.95 Å; atoms two bonds apart sit 2.2 Å or more"; earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | - | - |
| `interactions:_GRID` | constant | 1e3 | _GRID = 1e3 # measured values are snapped to 1/_GRID (Å, degrees) before any limit sees them | convention | "measured values are snapped to 1/_GRID (Å, degrees) before any limit sees them" (numeric precision); earliest in this file: 15b7a3c 2026-10-21 interactions: measured values snapped to a 1e-3 grid before any limit; numeric s… | - | - |
| `interactions:_REACH` | constant | 0.5 / _GRID | _REACH = 0.5 / _GRID # a tree search reaches this far past a limit, so a value that snaps onto it is not lost | geometry | half a grid step (0.5 / _GRID); earliest in this file: 15b7a3c 2026-10-21 interactions: measured values snapped to a 1e-3 grid before any limit; numeric s… | - | - |
| `interactions:_V_SLACK` | constant | 0.5 / _GRID + 1e-5 | _V_SLACK = 0.5 / _GRID + 1e-5 # a measured value the finders' grid could place on a limit passes it | geometry | half a grid step and a float margin; earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | - | - |
| `interactions:_V_BONDED` | constant | _COVALENT + 0.5 / _GRID - 1e-5 | _V_BONDED = _COVALENT + 0.5 / _GRID - 1e-5 # surely bonded, as the finders judge it | geometry | _COVALENT and half a grid step; earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | - | - |
| `interactions:_V_RINGS` | constant | { ("PHE", 6): {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}, ("TYR", 6): {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}, ("TRP", 5… | _V_RINGS = { | biology | ring atoms (chemistry); earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | - | - |
| `interactions:_V_BASE_RINGS` | constant | {6: {"N1", "C2", "N3", "C4", "C5", "C6"}, 5: {"C4", "C5", "N7", "C8", "N9"}} | _V_BASE_RINGS = {6: {"N1", "C2", "N3", "C4", "C5", "C6"}, 5: {"C4", "C5", "N7", "C8", "N9"}} | biology | ring atoms (chemistry); earliest in this file: 130c13b 2026-10-21 interactions.py: pi-stacking, hydrophobic contacts and verify_interactions, an i… | - | - |
| `interactions:_hydrophobic_carbons:argument:3@79318fae` | argument | 3 | return np.zeros((0, 3)), np.zeros(0, int), [] | convention | an array's shape | - | - |

## interactive

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `interactive:HIT_PX` | constant | 8 | HIT_PX = 8 | convention | a click target's width; earliest in this file: e2ddff4 2026-10-21 Interactive page: a click within 4 px of a DNA/RNA strand or stem selects it, at… | browser, interactive | - |
| `interactive:_links:comparison:6@ffef8bd1` | comparison | 6 | held = ", ".join(names[:6]) + (" …" if len(names) > 6 else "") | unknown | - | - | - |

## io

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `io:_PEPTIDE_BOND_MAX` | constant | 2.5 | _PEPTIDE_BOND_MAX = 2.5 # Å, C(i-1)-N(i), as in DSSP | biology | "C(i-1)-N(i), as in DSSP"; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `io:_CO_BOND` | constant | 1.231 | _CO_BOND = 1.231 # Å, carbonyl C=O | biology | "carbonyl C=O" bond length; earliest in this file: 824fbad 2026-10-09 Rebuild missing carbonyl oxygens; explain CA-only models | - | - |
| `io:_CO_MAX` | constant | 1.5 | _CO_MAX = 1.5 # Å, the farthest an oxygen still counts as the carbonyl one | biology | the carbonyl oxygen's reach; earliest in this file: 4fda3ae 2026-10-11 Peptide-linked residues without N/CA/C names (GFP's CRO66) stay in the chain and… | - | - |
| `io:AUTO_MAX` | constant | 24 | AUTO_MAX = 24 # `auto` builds assemblies up to this many chains (a capsid is better asked for by name) | convention | an assembly size limit; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `io:COINCIDENT_FRAC` | constant | 0.90 | COINCIDENT_FRAC = 0.90 # of a copy's atoms ... | unknown | earliest in this file: 13139b4 2026-10-11 C5: drop assembly copies laid exactly on a chain already present | - | - |
| `io:COINCIDENT_DIST` | constant | 1.0 | COINCIDENT_DIST = 1.0 # ... within this many Å of one chain already present: the same molecule seen twice | unknown | earliest in this file: 13139b4 2026-10-11 C5: drop assembly copies laid exactly on a chain already present | - | - |
| `io:SAME_PROTEIN` | constant | 0.9 | SAME_PROTEIN = 0.9 # two polymer entities are one protein when their full sequences align this identical (of the | unknown | introduced by 391d211 (4UM3), whose message states the rule's 90 % and then 4UM3's outcome, not a value set on 4UM3 (as figcheck:COPY_IDENTITY); its message: "entities whose full sequences align at least 90 % identical ... as one protein ... 4UM3 draws assembly 1's pentamer" | - | - |
| `io:SAME_LENGTH` | constant | 0.8 | SAME_LENGTH = 0.8 | unknown | introduced by 2a8ae8e (C30: 3OWW), whose message ties 0.8 to "a bound peptide", no entry; earliest in this file: 2a8ae8e 2026-10-14 C30 (random-100 fix round 1): a peptide is not a copy of its receptor; chains th… | - | - |
| `io:PAIRED_MIN` | constant | 3 | PAIRED_MIN = 3 # residue pairs of two chains' nucleotides with polar base atoms this near: held together (_paired) | named-entry:4PLX,3OWW | git: ddcbcbf's message (the value's last change): "PAIRED_MIN is 3: 3OWW's dimer pairs in 4 and stays A+B ... 4PLX draws assembly 1 again" | assembly | - |
| `io:_PAIR_DIST` | constant | 3.4 | _PAIR_DIST = 3.4 # angstrom | named-entry:3OWW | git: 2a8ae8e's message: "3OWW ruling: chains held together by base pairs ... within 3.4 A ... are one object" | - | - |
| `io:_NA_LINK_MAX` | constant | 2.0 | _NA_LINK_MAX = 2.0 # Å, O3'(i-1)-P(i) | biology | "O3'(i-1)-P(i)" bond; earliest in this file: ee95752 2026-10-09 Name the two link cutoffs apart; fix routed type hint; empty first protomer | - | - |
| `io:CONTACT` | constant | 4.0 | CONTACT = 4.0 # Å, heavy atoms of a protein residue and a nucleotide | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | datafile, na_contacts, na_validation | - |
| `io:_SS_MAX` | constant | 2.5 | _SS_MAX = 2.5 # Å between cysteine SG atoms | biology | "between cysteine SG atoms" (a disulfide); earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `io:_GLYCAN_LINK_MAX` | constant | 2.0 | _GLYCAN_LINK_MAX = 2.0 # Å between a sugar's anomeric carbon and the atom it is bonded to | biology | "a sugar's anomeric carbon and the atom it is bonded to" (a covalent bond); earliest in this file: ee95752 2026-10-09 Name the two link cutoffs apart; fix routed type hint; empty first protomer | - | - |
| `io:_METAL_REACH` | constant | 3.1 | _METAL_REACH = 3.1 # Å, metal to a coordinating atom (K+ to carbonyl O is ~2.8-3.0) | biology | "K+ to carbonyl O is ~2.8-3.0"; earliest in this file: df8de19 2026-10-09 Add ligand and metal ion markers | - | - |
| `io:_LIGAND_REACH` | constant | 4.0 | _LIGAND_REACH = 4.0 # Å, ligand heavy atom to a residue atom | unknown | earliest in this file: df8de19 2026-10-09 Add ligand and metal ion markers | - | - |
| `io:_BOND` | constant | 1.9 | _BOND = 1.9 # Å between heavy atoms of two residues that are covalently joined | biology | "heavy atoms of two residues that are covalently joined"; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | perf | - |
| `io:_drop_coincident:argument:3@b9797183` | argument | 3 | xyz = np.array([a.pos.tolist() for r in chain for a in r], float).reshape(-1, 3) | convention | an array's shape | assembly, na_validation | - |
| `io:_operator_tags:comparison:1e-3@561161d9` | comparison | 1e-3 | if not len(c) or not len(c[0]) or c[0][0].pos.dist(p) > 1e-3: | unknown | - | - | - |
| `io:_same_proteins:argument:20@7ee7e2fd` | argument | 20 | or min(len(seq), len(ref)) < max(20, SAME_LENGTH * max(len(seq), len(ref))) | unknown | line last changed by 2a8ae8e C30 (random-100 fix round 1): a peptide is not a copy of its receptor; chains that pair are one object (3OWW) (git blame; evidence only) | assembly | - |
| `io:read_model:argument:1e-3@bc170007` | argument | 1e-3 | if name not in claimed and len(o) == len(xyz) and np.allclose(o, xyz, atol=1e-3): | unknown | - | app_scope, assembly, datafile, hbonds +11 | - |
| `io:load_backbone:argument:4@3fcf5457` | argument | 4 | return Backbone([], np.zeros((0, 4, 3)), np.zeros(0, bool), set(), np.zeros(0)) # drawn on its own | convention | an array's shape | 310, align, assembly, datafile +26 | - |
| `io:load_backbone:argument:3@3fcf5457` | argument | 3 | return Backbone([], np.zeros((0, 4, 3)), np.zeros(0, bool), set(), np.zeros(0)) # drawn on its own | convention | an array's shape | 310, align, assembly, datafile +26 | - |
| `io:load_nucleic:argument:5@0a63058f` | argument | 5 | ns = gemmi.NeighborSearch(model, st.cell, 5).populate() | unknown | - | assembly, datafile, independent_kernels, na_pairs +4 | - |
| `io:_ligands:argument:5@c5dac063` | argument | 5 | search = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate() | unknown | - | ligands | - |

## layout

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `layout:PITCH` | constant | 1.1 | PITCH = 1.1 # distance between neighbouring strand centres | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | copies, domains, layout, na_metrics +3 | - |
| `layout:SHAFT` | constant | 0.62 | SHAFT = 0.62 # arrow shaft width | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_metrics, page_text | - |
| `layout:SCALE` | constant | 0.23 | SCALE = 0.23 # page units per Å (PITCH / 4.8 Å strand spacing) | geometry | "PITCH / 4.8 Å strand spacing"; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | 310, layout, like, membrane +6 | - |
| `layout:HEAD` | constant | 1.0 | HEAD = 1.0 # strand arrowhead width | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout | - |
| `layout:HEAD_LEN` | constant | 0.7 | HEAD_LEN = 0.7 # strand arrowhead length (at most 0.45 of a short strand) | unknown | earliest in this file: 2473be0 2026-10-11 A helix over a sheet of its own colour takes a darker tone; its name keeps off a… | layout | - |
| `layout:HELIX_W` | constant | 1.0 | HELIX_W = 1.0 # helix ribbon diameter, ~4.4 Å: one arrowhead wide, as in a cartoon (Style.helix_scale) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | 310, copies, domains, heldout_budgets +5 | - |
| `layout:STRAND_RISE` | constant | 3.3 | STRAND_RISE = 3.3 # Å per residue along a strand | biology | Å per residue along a strand; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout | - |
| `layout:HELIX_RISE` | constant | 1.5 | HELIX_RISE = 1.5 # Å per residue along a helix | biology | Å per residue along an α-helix; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout | - |
| `layout:TEN_RISE` | constant | 2.0 | TEN_RISE = 2.0 # Å per residue along a 3-10 helix | biology | Å per residue along a 3-10 helix; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | 310 | - |
| `layout:ETA_W` | constant | 0.7 | ETA_W = 0.7 # 3-10 box width as a fraction of the helix width | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | route | - |
| `layout:MARGIN` | constant | 0.8 | MARGIN = 0.8 # clear space kept around each element so loops can pass | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | copies, heldout_budgets, membrane, na_curve +2 | - |
| `layout:_ROW_GAP` | constant | 2.5 | _ROW_GAP = 2.5 # between neighbours in a stack row | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_ROW_SPACING` | constant | 3.2 | _ROW_SPACING = 3.2 # between stack rows: room for a chain-end label below one row and above the next | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:LABEL_GAP` | constant | 0.3 | LABEL_GAP = 0.3 # a helix label's centre sits this far beyond the side of the element box (times glyph_scale) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | label_scale | - |
| `layout:_NAME_PAD` | constant | 0.05 | _NAME_PAD = 0.05 # a helix name's keep-out box reaches this far past its glyphs (times glyph_scale) | unknown | earliest in this file: ae43b7f 2026-10-21 Helix names keep their gap to their bars at any font_scale (plan 22 Task 1) | label_scale | - |
| `layout:END_STUB` | constant | 0.7 | END_STUB = 0.7 # chain-end stub length | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:END_LABEL` | constant | 1.15 | END_LABEL = 1.15 # distance from a chain-end port to its N/C letter | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | route, text_clear | - |
| `layout:LONG_TAIL` | constant | 10 | LONG_TAIL = 10 # a chain end with this many residues or more past its last element gets a count ("14 aa") | unknown | introduced by bc5578c (C1: 29RY, 2DGS, 8W9Y), whose message does not tie this value to those entries; earliest in this file: bc5578c 2026-10-13 C1 (random-100): tail counts placed clear of loops, elements and text, from the … | - | - |
| `layout:TAIL_COUNT_SIZE` | constant | 0.72 | TAIL_COUNT_SIZE = 0.72 # that count's text size, as a fraction of the figure's | unknown | introduced by bc5578c (C1: 29RY, 2DGS, 8W9Y), whose message gives "0.72 em" but ties it to no entry; earliest in this file: bc5578c 2026-10-13 C1 (random-100): tail counts placed clear of loops, elements and text, from the … | text_clear | - |
| `layout:_TAIL_COUNT_GAP` | constant | 0.08 | _TAIL_COUNT_GAP = 0.08 # page units between the N/C letter's ink and the count's | unknown | introduced by bc5578c (C1: 29RY, 2DGS, 8W9Y), whose message does not tie this value to those entries; earliest in this file: bc5578c 2026-10-13 C1 (random-100): tail counts placed clear of loops, elements and text, from the … | - | - |
| `layout:_TAIL_COUNT_NEAR` | constant | 1.5 | _TAIL_COUNT_NEAR = 1.5 # em: the most a tail count's ink stands from its letter's (it reads as the letter's) | named-entry:7SMR | git: 0503fe6's message: "within _TAIL_COUNT_NEAR (1.5 em) of the letter's ink ... 7SMR's C:E "18 aa", 2.8 page units below its letter" | text_clear | - |
| `layout:_TAIL_COUNT_FAR` | constant | (0.3, 0.6, 0.9, 1.2) | _TAIL_COUNT_FAR = (0.3, 0.6, 0.9, 1.2) # page units further out a count looks when its near spots are all taken | unknown | earliest in this file: 3124a38 2026-10-15 C1 (random-100, batch C2a): a tail count with every near spot taken looks furthe… | - | - |
| `layout:_GLYPH_W` | constant | 0.32 | _GLYPH_W, _GLYPH_H = 0.32, 0.42 # rough label glyph size in page units, for keep-out boxes (at _EM_TUNED: glyph_scale) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | label_scale, na_layout, rosette | - |
| `layout:_GLYPH_H` | constant | 0.42 | _GLYPH_W, _GLYPH_H = 0.32, 0.42 # rough label glyph size in page units, for keep-out boxes (at _EM_TUNED: glyph_scale) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_layout | - |
| `layout:_GLYPH_K` | constant | ContextVar("foldmap_glyph_scale", default=1.0) | _GLYPH_K: ContextVar[float] = ContextVar("foldmap_glyph_scale", default=1.0) # glyph_scale of the text the layout | convention | a context default (glyph scale 1); earliest in this file: f22c038 2026-10-15 C6 fix round 1 (findings 4, 5, 8): every glyph-sized room goes through layout.gl… | label_scale, page_text | - |
| `layout:IN_PER_UNIT` | constant | 0.30 | IN_PER_UNIT = 0.30 # render: inches per page unit, at most | unknown | earliest in this file: 7104faa 2026-10-11 C6/C10 follow-ups from the rerun: a gap the grid cannot route is curved before i… | page_text | - |
| `layout:MAX_WIDTH_IN` | constant | 7.2 | MAX_WIDTH_IN = 7.2 # render: a figure is drawn at most this wide (double column), inches, while its text stays ... | convention | "double column" page width; earliest in this file: 7104faa 2026-10-11 C6/C10 follow-ups from the rerun: a gap the grid cannot route is curved before i… | page_text | - |
| `layout:MIN_TEXT_PT` | constant | 6.0 | MIN_TEXT_PT = 6.0 # ... this size or larger: a wider page grows past MAX_WIDTH_IN instead (random-100 C4 ruling) | convention | the smallest text size printed (random-100 C4 ruling); earliest in this file: a3cf24d 2026-10-15 C4 (random-100, batch C1): text has a 6 pt floor and the page grows instead; leg… | page_text | - |
| `layout:MAX_PAGE_IN` | constant | 48.0 | MAX_PAGE_IN = 48.0 # render: the longest side a page grows to (A0's, about); past it the whole figure, its text | convention | "A0's, about"; earliest in this file: a3cf24d 2026-10-15 C4 (random-100, batch C1): text has a 6 pt floor and the page grows instead; leg… | app_jobs, heldout_budgets, mcp, page_text | - |
| `layout:_MAX_ROW` | constant | 30.0 | _MAX_ROW = 30.0 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | membrane, route | - |
| `layout:_MIN_STRAND` | constant | 1.0 | _MIN_STRAND = 1.0 # shortest drawn strand (room for the arrow head) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_PAIR_ENERGY_ABOVE` | constant | 200 | _PAIR_ENERGY_ABOVE = 200 | unknown | earliest in this file: d1fa11e 2026-10-10 Tests compare the pair-list energy with the production dense energy | optimiser, perf | - |
| `layout:_W_SEQ` | constant | 1.0 | _W_SEQ = 1.0 # pull between the two ends of a loop | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout | - |
| `layout:_LOOP_SLACK` | constant | 1.5 | _LOOP_SLACK = 1.5 # loop ends closer than this are left alone | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | domains | - |
| `layout:_W_SS` | constant | 2.0 | _W_SS = 2.0 # pull between two bridged cysteines (twice a loop: disulfides should read as short bars) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_SS_SLACK` | constant | 0.6 | _SS_SLACK = 0.6 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_SS_RESTART_PAIRS` | constant | 5_000_000_000 | _SS_RESTART_PAIRS = 5_000_000_000 | measured | "Measured at plan 15 Task 7 on the default and focus views of both random-100 samples" (8D3D); earliest in this file: 1bf5f36 2026-10-18 Bounded work on very large pages (plan 15 task 7, A1): disulfide second starts s… | heldout_budgets | - |
| `layout:_W_CONTACT` | constant | 0.3 | _W_CONTACT = 0.3 # pull between elements that touch in 3D (at _CONTACT_CAP CA pairs) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_layout, na_validation | - |
| `layout:_CONTACT_CAP` | constant | 20 | _CONTACT_CAP = 20 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_K_HOME` | constant | 0.02 | _K_HOME = 0.02 # pull back toward the projected position, keeps the fold recognisable | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | heldout_budgets, sensitivity | - |
| `layout:_K_LIKE` | constant | 128 * _K_HOME | _K_LIKE = 128 * _K_HOME | measured | "Frozen on the homologue tuning set (13a-10): 2, 8, 32, 64, 128 tried" (2ZB5, 1E8T); earliest in this file: c2470a8 2026-10-21 A layout drawn like a reference (plan 13a Task 5): frame_like, Like, homes | sensitivity | - |
| `layout:_STEP` | constant | 0.15 | _STEP = 0.15 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_SWEEPS` | constant | 400 | _SWEEPS = 400 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_SEP_SWEEPS` | constant | 5000 | _SEP_SWEEPS = 5000 # push-apart passes before the guaranteed slide | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | layout, perf | - |
| `layout:_SEP_ORBIT_SWEEPS` | constant | 1000 | _SEP_ORBIT_SWEEPS = 1000 # push-apart passes that keep symmetry copies alike, before they may move alone | unknown | earliest in this file: 9064749 2026-10-10 fix: symmetry copies are separated together; pilus check uses each chain's media… | symmetry | - |
| `layout:_SEP_SYM_ROOM` | constant | 1.05 | _SEP_SYM_ROOM = 1.05 # copies drawn alike may take at most this much more area than items moved one by one | unknown | earliest in this file: 441ae48 2026-10-10 fix: symmetric separation kept only if it costs at most 5 % more area than movin… | copies | - |
| `layout:_SEP_MIRROR_ROOM` | constant | 1.10 | _SEP_MIRROR_ROOM = 1.10 | named-entry:9GEP | "…s moved one by one: the figure says they are mirror images (9GEP's nucleosome needed 1.07 and about 3,000 passes)"; earliest in this file: 9eee90e 2026-10-15 C12 (random-100, batch C2a): mirror-image copies are kept for a little more room | symmetry | - |
| `layout:_SEP_TIE_ROOM` | constant | 2.0 | _SEP_TIE_ROOM = 2.0 | named-entry:8QEM,6JLN,7QUP,9P1X | "… second chance) may take up to this much (random-100 C2b-3: 8QEM's twelve ring copies drawn whole need 1.61 times the room of moving items one by on…"; earliest in this file: 0bf7a26 2026-10-16 Copies that need more room are still drawn alike: a second chance with the copy-… | copies, symmetry | - |
| `layout:_SEP_PUSH` | constant | 0.005 | _SEP_PUSH = 0.005 # each push clears an overlap by this much more, so touching boxes cannot trade it back and forth | unknown | earliest in this file: 5fce03d 2026-10-10 Push-apart that converges on crowded layouts; fallback slide to the nearest clea… | - | - |
| `layout:_SLIDE_STEP` | constant | 0.5 | _SLIDE_STEP = 0.5 # the guaranteed slide moves an item by this much at a time | unknown | earliest in this file: 5fce03d 2026-10-10 Push-apart that converges on crowded layouts; fallback slide to the nearest clea… | - | - |
| `layout:_LAMBDAS` | constant | (0.3, 3.0, 30.0, 300.0) | _LAMBDAS = (0.3, 3.0, 30.0, 300.0) # overlap penalty, stage by stage | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | optimiser | - |
| `layout:_MAX_MOVE` | constant | 0.5 | _MAX_MOVE = 0.5 # largest step an item takes in one sweep | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:_GRID` | constant | 1e4 | _GRID = 1e4 # optimiser inputs and item geometry are rounded to 1/_GRID page units (radians for angles) | convention | "optimiser inputs and item geometry are rounded to 1/_GRID page units" (determinism); earliest in this file: d067d4e 2026-10-10 fix: deterministic layout optimiser: numpy-only L-BFGS on inputs snapped to a 1e… | - | - |
| `layout:DNA_RISE` | constant | 3.4 | DNA_RISE = 3.4 # Å per base pair (B-DNA) | biology | "Å per base pair (B-DNA)"; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_curve | - |
| `layout:DNA_W` | constant | 20.0 * SCALE | DNA_W = 20.0 * SCALE # B-DNA diameter | biology | "B-DNA diameter" (20 Å); earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | style | - |
| `layout:DNA_PITCH` | constant | 10.5 * DNA_RISE * SCALE | DNA_PITCH = 10.5 * DNA_RISE * SCALE # one turn, page units | biology | one B-DNA turn: 10.5 base pairs; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:DNA_MINOR` | constant | 0.37 | DNA_MINOR = 0.37 # the second backbone trails the first by this fraction of a turn (minor groove) | biology | "the second backbone trails the first by this fraction of a turn (minor groove)"; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layout:RAIL_IN` | constant | 0.6 | RAIL_IN = 0.6 # ladder rails sit this far inside the stem's long sides, base letters outside them | unknown | earliest in this file: 72ceeed 2026-10-10 Stems placed as layout elements from the 3D structure and drawn as ladders | na_layout | - |
| `layout:NT_STEP` | constant | 0.55 | NT_STEP = 0.55 # spacing of unpaired nucleotides along a drawn single strand (schematic, not to scale) | unknown | earliest in this file: 72ceeed 2026-10-10 Stems placed as layout elements from the 3D structure and drawn as ladders | na_layout | - |
| `layout:_HAIRPIN_STEP` | constant | 0.4545 | _HAIRPIN_STEP = 0.4545 # page units: the least step between neighbours on a hairpin loop, a bead's width (_TRACE_HALF, | geometry | "a bead's width (_TRACE_HALF" twice; earliest in this file: ed8d7f2 2026-10-17 A hairpin loop's circle grows until its nucleotides stand at least a contact bea… | - | - |
| `layout:PIECE_GAP` | constant | 0.3 | PIECE_GAP = 0.3 # page units between two pieces of one helix drawn end to end (_bundle_item, _regroup_copies, | unknown | introduced by 4691767 (C2b-3: 6JLN, 9P1X, 7QUP), whose message calls it "the gap _bundle_item already used"; earliest in this file: 4691767 2026-10-16 A copy's sheets are never regrouped to mimic its mate; lone pieces stay with the… | - | - |
| `layout:_STACK_PAD` | constant | 1.3 | _STACK_PAD = 1.3 # page units between two strands drawn end to end in one column: room for the loops leaving | measured | "swept 1.2-1.6 on the 37 sample entries with accessory strands (5 views) and the 30QM and 4YBB fixtures"; earliest in this file: a2fadd7 2026-10-15 C34 (random-100, batch C2a): a bulged continuation shares its strand's column | layout, sensitivity | - |
| `layout:_STACK_SHUT` | constant | 0.6 | _STACK_SHUT = 0.6 | named-entry:9P1X,8OZJ | "… either side standing past the gap's ends but for a sliver: 9P1X 0.15, 8OZJ 0.54), the loops leaving its facing ends have no way out, or squeeze thr…"; earliest in this file: 0df3a22 2026-10-16 A stack gap shut on both sides opens: the strand stacked last stands further off… | layout | - |
| `layout:_STACK_OPEN` | constant | 0.85 | _STACK_OPEN = 0.85 | geometry | "two loops' grid cells, route.RES, each _INFLATE clear of the strands"; earliest in this file: 0df3a22 2026-10-16 A stack gap shut on both sides opens: the strand stacked last stands further off… | - | - |
| `layout:_STACK_NUDGE` | constant | 1.0 | _STACK_NUDGE = 1.0 | named-entry:9P1X,8OZJ | git: 0df3a22's message: "by at most _STACK_NUDGE 1.0 in 0.05 steps" (9P1X, 8OZJ: slivers of 0.15 and 0.54) | - | - |
| `layout:_NUDGE_STEP` | constant | 0.05 | _NUDGE_STEP = 0.05 | unknown | introduced by 0df3a22 (9P1X, 8OZJ), whose message names it only as the search step, not a value set on those entries; its message: "by at most _STACK_NUDGE 1.0 in 0.05 steps" (9P1X, 8OZJ: slivers of 0.15 and 0.54) | - | - |
| `layout:_POCKET_REACH` | constant | 2.0 | _POCKET_REACH = 2.0 | named-entry:8BCU,12BM | "…s partner or not, until it closes none (random-100 held-out 8BCU: G:111-114, level with G:119-122, closed the run four loops needed; 12BM: C:71-79 c…"; earliest in this file: 99d40e4 2026-10-17 Held-out A3/A5 (plan 15 task 1): route fallbacks and overlapping boxes, by mecha… | - | - |
| `layout:_TOUCH` | constant | 1e-6 | _TOUCH = 1e-6 # page units two drawn spans must share to stand level (science._TOUCH) | convention | a float tolerance (1e-6); earliest in this file: 7b0292a 2026-10-17 An accessory strand never stands level with a strand it does not pair with in th… | - | - |
| `layout:_CLEAR` | constant | 1e-4 | _CLEAR = 1e-4 # page units past a block's end a cleared height stands, one snap step (_column_heights) | geometry | "one snap step" (1/_GRID); earliest in this file: a287c05 2026-10-10 fix: stems the projection cannot side keep strands uncrossed; short runs hold st… | - | - |
| `layout:_END_ROOM` | constant | END_LABEL - 0.15 | _END_ROOM = END_LABEL - 0.15 | named-entry:9FWZ | "…ces: its stub and N/C letter. Fix round 3: at most this, so 9FWZ's D:662-666, end to end with D:634-641, still stands level with its partner D:645-6…"; earliest in this file: 97e3f70 2026-10-15 Plan 12 T2: off-path strands out of the sheet track; barrels are their pairing r… | layout | - |
| `layout:_BA_GAP` | constant | 0.3 | _BA_GAP = 0.3 # page units between the footprints of neighbours in a βα template row | unknown | earliest in this file: 1a9918c 2026-10-11 βα template: all-parallel α/β sheets drawn as one row with each crossover helix … | - | - |
| `layout:_BA_MAX` | constant | 10 | _BA_MAX = 10 # strands: a longer parallel sheet (2BNH's 17-strand leucine-rich repeat) is no α/β barrel or Rossmann | named-entry:2BNH | "strands: a longer parallel sheet (2BNH's 17-strand leucine-rich repeat) is no α/β barrel or Rossmann"; earliest in this file: 85200d0 2026-10-11 A βα row holds at most ten strands | template | - |
| `layout:_BA_SHARE` | constant | (3, 5) | _BA_SHARE = (3, 5) # the template applies when at least this share of a sheet's βα units are in its row | unknown | earliest in this file: 1a9918c 2026-10-11 βα template: all-parallel α/β sheets drawn as one row with each crossover helix … | - | - |
| `layout:_KINK_CORNER` | constant | 5 | _KINK_CORNER = 5 # residues: the longest corner between two pieces of a kinked helix (a π turn) | unknown | earliest in this file: a975518 2026-10-11 Kinked helix pieces join only across a short corner of their own chain | layout | - |
| `layout:_KINK_ROOM` | constant | 0.25 | _KINK_ROOM = 0.25 # one routing cell (route.RES): the corner loop's room on each leg between two kinked pieces | geometry | "one routing cell (route.RES)"; earliest in this file: 46b2e90 2026-10-11 Kinked helix pieces are drawn end to end as one unit, so the globin AB corner re… | - | - |
| `layout:_DENSE_PAIRS` | constant | 256 | _DENSE_PAIRS = 256 # _overlapping tests every pair of up to this many footprints at once; more, only near ones | convention | a work bound: every pair tested at once up to it, near ones only past it; earliest in this file: ef3f6f0 2026-10-15 C33 (random-100, batch C1): the slow views no longer lay out twice nor test ever… | perf | - |
| `layout:_JOINT_GROUP` | constant | 4 | _JOINT_GROUP = 4 # at most this many consecutive pieces of a duplex are separated as one group | unknown | earliest in this file: 5377142 2026-10-11 C3 fix: move joined pieces together only in chains of at most four; gate straigh… | - | - |
| `layout:_JOINT_OVERLAP` | constant | 0.25 | _JOINT_OVERLAP = 0.25 # widths two joined pieces' bodies may overlap and still be moved together | unknown | earliest in this file: 5377142 2026-10-11 C3 fix: move joined pieces together only in chains of at most four; gate straigh… | na_layout | - |
| `layout:_STEM_ROOM` | constant | 1.5 | _STEM_ROOM = 1.5 # the room a nucleic acid's stems take on the page: their footprints' area (with half a MARGIN round | measured | "Measured on the NA-only entries of the NA campaign and random-100" (3OWW, 1GID, 3IGI, 4LVV, 4P9R); earliest in this file: ea0765b 2026-10-17 A nucleic acid drawn without protein is spread by one factor to the room its ste… | - | - |
| `layout:_MIN_RISE` | constant | 2.6 | _MIN_RISE = 2.6 # Å per base pair, the least a stem is drawn with (A-form rises 2.8, B-form 3.4) | unknown | earliest in this file: 72ceeed 2026-10-10 Stems placed as layout elements from the 3D structure and drawn as ladders | na_layout | - |
| `layout:_SIDE_CLEAR` | constant | 2.75 | _SIDE_CLEAR = 2.75 | measured | "On the test structures every stem below it fits its projection about as well turned round; the closest one that does not is 4YBB P5 at 2.81 Å"; earliest in this file: a287c05 2026-10-10 fix: stems the projection cannot side keep strands uncrossed; short runs hold st… | na_layout | - |
| `layout:_CYCLE_MAX` | constant | 12 | _CYCLE_MAX = 12 # unclear stems linked round a cycle: up to this many are settled by trying every turning | unknown | earliest in this file: 03d6d23 2026-10-10 fix: strand-exchange springs only; hairpin room holds its letters; exact side co… | - | - |
| `layout:_TRACE_HALF` | constant | 0.45 | _TRACE_HALF = ( | unknown | earliest in this file: 334ee0b 2026-10-11 C9 revised: each trace nucleotide is a fixed obstacle, and an element touching t… | na_layout, na_text_keepoff | - |
| `layout:_WALL_MIN` | constant | 3 | _WALL_MIN = 3 # a trace of fewer nucleotides (a dinucleotide ligand-like oligo) is drawn but is no wall | unknown | earliest in this file: 996b6f6 2026-10-11 Trace walls only for traces of at least 3 nucleotides | - | - |
| `layout:_NUDGES` | constant | tuple(f * k for k in (1.05, 1.25, 1.5, 2.0) for f in (1, -1)) | _NUDGES = tuple(f * k for k in (1.05, 1.25, 1.5, 2.0) for f in (1, -1)) # sideways steps, in beads (_TRACE_HALF) | unknown | earliest in this file: 32878dd 2026-10-11 Trace nucleotides projected within a bead of one drawn before them are nudged si… | - | - |
| `layout:_W_JOINT` | constant | 4.0 | _W_JOINT = 4.0 # pull between consecutive pieces of one bent duplex: the chain follows the duplex | unknown | earliest in this file: 72ceeed 2026-10-10 Stems placed as layout elements from the 3D structure and drawn as ladders | na_layout | - |
| `layout:_JOINT_SLACK` | constant | 0.3 | _JOINT_SLACK = 0.3 | unknown | earliest in this file: 72ceeed 2026-10-10 Stems placed as layout elements from the 3D structure and drawn as ladders | na_layout | - |
| `layout:_SHORT_RUN` | constant | 2 | _SHORT_RUN = 2 | unknown | earliest in this file: a287c05 2026-10-10 fix: stems the projection cannot side keep strands uncrossed; short runs hold st… | - | - |
| `layout:_W_SHORT` | constant | 10.0 | _W_SHORT = 10.0 | unknown | earliest in this file: a287c05 2026-10-10 fix: stems the projection cannot side keep strands uncrossed; short runs hold st… | na_layout | - |
| `layout:_HELIX_LOOP` | constant | 3 | _HELIX_LOOP = 3 | measured | "Measured over the 39 campaign entries" (4YBB, 5F9R); earliest in this file: 0e226b2 2026-10-11 Stems of one helix across an internal loop of at most 3 nt a side are held close… | - | - |
| `layout:_JOIN_TRIES` | constant | 40 | _JOIN_TRIES = 40 # _SLIDE_STEP moves an unpaired member joining its regrouped copy may take to clear its members | unknown | introduced by 4691767 (C2b-3: 6JLN, 9P1X, 7QUP), whose message does not tie this value to those entries; earliest in this file: 4691767 2026-10-16 A copy's sheets are never regrouped to mimic its mate; lone pieces stay with the… | layout | - |
| `layout:_WIDEN_GAP` | constant | 0.1 | _WIDEN_GAP = 0.1 # page units kept between two members of one item that a width set by hand pushed apart | unknown | earliest in this file: 2fd5f74 2026-10-22 Each element's own colour and width (plan 23 Task 1): --color/--width, element_c… | - | - |
| `layout:BAND_TOUCH` | constant | 1e-9 | BAND_TOUCH = 1e-9 # page units: a drawn end this near a face reaches it (float noise only: 6DQN's ports, ~4e-15 short) | convention | "float noise only: 6DQN's ports, ~4e-15 short"; earliest in this file: 24ad3d2 2026-10-24 One crossing rule for the membrane view, A18 and data; sides from the loops (pla… | membrane | - |
| `layout:_PAST` | constant | 0.1 | _PAST = 0.1 # page units: _to_face's own moves take an end this far past a face or the middle (0.4 Å at SCALE) | unknown | earliest in this file: ecd734d 2026-10-24 Task 1 fix round 1: strict middle, soluble elements in the slab, CI tests for th… | membrane | - |
| `layout:_BOUND_SLACK` | constant | 1e-3 | _BOUND_SLACK = 1e-3 # page units (0.004 Å) kept inside a bound, so a place snapped to the grid (_GRID) stays within it | unknown | cites rulings P23; earliest in this file: ffde4b5 2026-10-25 Membrane band laid out first, every other element held by a one-sided bound (pla… | - | - |
| `layout:_FACE_ROWS` | constant | 3 | _FACE_ROWS = 3 # rows a face holds its flat helices in, the first on the face, the others stacked outward | unknown | introduced by 0a35b86 (plan 26: 5LDW), whose message ties 5LDW to _PAST, not this value; earliest in this file: 0a35b86 2026-10-24 Task 2 fix round 1: face rows, columns in order, 5LDW crossing, side names off l… | - | - |
| `layout:_ROW_TOUCH` | constant | 1e-3 | _ROW_TOUCH = ( | convention | a float tolerance (1e-3 page units); earliest in this file: 0a35b86 2026-10-24 Task 2 fix round 1: face rows, columns in order, 5LDW crossing, side names off l… | - | - |
| `layout:_FACE_PAD` | constant | 5.0 | _FACE_PAD = 5.0 # page units: a chain's transmembrane span, padded by this or a face helix's length, holds its helices | unknown | introduced by 0a35b86 (plan 26: 5LDW), whose message ties 5LDW to _PAST, not this value; earliest in this file: 0a35b86 2026-10-24 Task 2 fix round 1: face rows, columns in order, 5LDW crossing, side names off l… | membrane | - |
| `layout:_ORDER_ROUNDS` | constant | 3 | _ORDER_ROUNDS = 3 # build_layout: deal the band columns back their projected order and separate, at most this often | unknown | earliest in this file: 26b9ad4 2026-10-24 Task 5 fix round 1: band columns kept in projected order, review-5 minors (plan … | - | - |
| `layout:AXIS_SYMBOL_R` | constant | 0.45 | AXIS_SYMBOL_R = 0.45 # page units: a rosette's n-fold symbol at the page centre (render), corner to centre | unknown | earliest in this file: 80a381c 2026-10-20 Rosette view of a cyclic assembly, opt-in (style cyclic=rosette, default unrolle… | rosette | - |
| `layout:_ROSETTE_STEP` | constant | 0.25 | _ROSETTE_STEP = 0.25 # page units: the rosette's radial push grows by this much until neighbouring copies clear | unknown | earliest in this file: 80a381c 2026-10-20 Rosette view of a cyclic assembly, opt-in (style cyclic=rosette, default unrolle… | - | - |
| `layout:_ROSETTE_GAP` | constant | 2 * MARGIN | _ROSETTE_GAP = 2 * MARGIN # clear space kept between two copies' elements (with their names): room for loops | geometry | 2 * MARGIN; earliest in this file: 80a381c 2026-10-20 Rosette view of a cyclic assembly, opt-in (style cyclic=rosette, default unrolle… | rosette | - |
| `layout:_ROSETTE_SAMPLE` | constant | 0.4 | _ROSETTE_SAMPLE = 0.4 # page units between the points an element is sampled at along its axis, for that clearance | unknown | earliest in this file: 80a381c 2026-10-20 Rosette view of a cyclic assembly, opt-in (style cyclic=rosette, default unrolle… | - | - |
| `layout:_ROSETTE_MAX_PUSH` | constant | 200.0 | _ROSETTE_MAX_PUSH = 200.0 # page units: the push never grows past this (a copy farther out clears by then anyway) | unknown | earliest in this file: 80a381c 2026-10-20 Rosette view of a cyclic assembly, opt-in (style cyclic=rosette, default unrolle… | - | - |
| `layout:_UNIT_AGREES` | constant | 1e-9 | _UNIT_AGREES = 1e-9 # page units: a Like.unit is read only when it gives its target's place to this | convention | a float tolerance (1e-9); earliest in this file: 790f263 2026-10-21 Drawn like a reference, fix round 1 (plan 13a Task 5 review) | - | - |
| `layout:DOMAIN_PAD` | constant | 0.7 | DOMAIN_PAD = 0.7 # panel margin around a domain's elements | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | membrane | - |
| `layout:_DOMAIN_NAME_ROOM` | constant | 0.6 | _DOMAIN_NAME_ROOM = 0.6 # a domain panel's room above its elements for its name (domain_panel) | unknown | earliest in this file: 6785a8e 2026-10-24 Domains with a membrane: one membrane domain per chain, soluble domains on their… | - | - |
| `layout:_PANEL_CLEAR` | constant | 0.1 | _PANEL_CLEAR = 0.1 # page units a soluble domain's panel keeps clear of a membrane face (_limits) | unknown | earliest in this file: 6785a8e 2026-10-24 Domains with a membrane: one membrane domain per chain, soluble domains on their… | - | - |
| `layout:_DOMAIN_GAP` | constant | 2 * DOMAIN_PAD + 0.6 | _DOMAIN_GAP = 2 * DOMAIN_PAD + 0.6 # clear space between two domains' elements: both panels plus a gutter | geometry | "both panels plus a gutter" (2 * DOMAIN_PAD + 0.6); earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | membrane | - |
| `layout:_PART_BLOWUP` | constant | 2.0 | _PART_BLOWUP = 2.0 # parting that left domains unparted and the items' width + height this many times: undone | named-entry:1WWL,2ILU,7SMR,6JLN,5LDW | git: 34a54df's message: "only when they leave the domains unparted and the items over _PART_BLOWUP (2) times as wide and high" (1WWL, 2ILU, 7SMR, 6JLN, 5LDW; 3L9G unchanged) | - | - |
| `layout:_PLATE_PAD` | constant | 0.3 | _PLATE_PAD = 0.3 # a sheet's plate reaches this far beyond its strands (render._draw_sheet_panels) | unknown | earliest in this file: 7fc8801 2026-10-12 A domain's name takes the top right of its panel when the top left is taken | - | - |
| `layout:DOMAIN_NAME_SIZE` | constant | 1.05 | DOMAIN_NAME_SIZE = 1.05 # a domain's name, in ems of the figure's text (render draws it bold at this size) | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | domains | - |
| `layout:_NAME_HEIGHT` | constant | 1.0 | _NAME_HEIGHT = 1.0 # ems of its own size a name's box takes up and down (capitals, a descender, the plate's room) | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `layout:_NAME_INSET` | constant | 0.35 | _NAME_INSET = 0.35 # page units a name keeps from its panel's side | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | domains | - |
| `layout:_NAME_SINK` | constant | 0.1 | _NAME_SINK = 0.1 # page units a name's middle sits inside its panel's top edge: it lies on the edge, on its plate | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | domains | - |
| `layout:_NAME_STEP` | constant | 0.25 | _NAME_STEP = 0.25 # page units between the places tried along the top edge | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `layout:_NAME_CLEAR` | constant | 0.08 | _NAME_CLEAR = 0.08 # page units a name keeps from anything else | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `layout:_NAME_KEEP` | constant | 0.3 | _NAME_KEEP = 0.3 # page units loops and nucleic-acid strands keep from a name: room kept round it where there is | unknown | earliest in this file: 94ddff9 2026-10-12 Domain names keep room from loops and the nucleic acid (4OO8 RuvC-II) | domains | - |
| `layout:_NAME_ROUTE` | constant | 0.15 | _NAME_ROUTE = 0.15 # (a wider keep-out forces loops over it where room is short) | unknown | earliest in this file: 94ddff9 2026-10-12 Domain names keep room from loops and the nucleic acid (4OO8 RuvC-II) | domains | - |
| `layout:_NAME_OPEN` | constant | _NAME_ROUTE + 0.1 | _NAME_OPEN = _NAME_ROUTE + 0.1 # page units a name keeps from a loop's way out of an element's end (_ways_out, 0.2 | named-entry:7NZ2 | "random-100 held-out 7NZ2's domains view: "D1 · B1" 0.08 from A1:1330-1337's way out" (comment carried under it); earliest in this file: bc0e375 2026-10-17 Held-out domain names and loops (plan 15 task 2): names and leaders off the ways… | domains | - |
| `layout:_NAME_NEAR_NT` | constant | 0.3 | _NAME_NEAR_NT = 0.3 # page units a name keeps from a nucleotide of an unpaired strand, beyond its letter | unknown | "a leader's keep-out keeps the same 0.1 page units a name keeps from a nucleotide"; earliest in this file: 94ddff9 2026-10-12 Domain names keep room from loops and the nucleic acid (4OO8 RuvC-II) | - | - |
| `layout:_NAME_STUB` | constant | 0.8 | _NAME_STUB = 0.8 # page units: the straight run of a loop out of an element's end, kept clear of a name | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | domains | - |
| `layout:_SHORT_LOOP` | constant | 6.0 | _SHORT_LOOP = 6.0 # page units: a loop whose ends are this near takes the short way between them, kept clear too | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `layout:_LEADER` | constant | 0.8 | _LEADER = 0.8 # page units a name moved outside its panel stands above it, joined to it by a short leader | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | domains | - |
| `layout:NT_ROW` | constant | 0.55 | NT_ROW = 0.55 # a base letter sits this far out from its nucleotide | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_curve, na_layout, na_text_keepoff | - |
| `layout:_SLOTS` | constant | (0.25, 0.5, 0.75) | _SLOTS = (0.25, 0.5, 0.75) # page units a stem's name may move further out than beside its letters, when neither | unknown | earliest in this file: 341318d 2026-10-11 Tie letters at base-letter size with a legend entry; stem labels take the side c… | - | - |
| `layout:_INK_MARGIN` | constant | 1.0 | _INK_MARGIN = 1.0 # ems kept clear between a stem's name and every base letter (its own letter row, or its rail | unknown | earliest in this file: e9a3b53 2026-10-12 Stem names keep a 0.25 em margin clear of base letters | figcheck, na_curve, na_layout, random100 | - |
| `layout:_NAME_EM` | constant | 0.8 | _NAME_EM, _LETTER_EM = 0.8, 0.75 # ems: render draws a stem's name, and a base letter, at these sizes | geometry | "render draws a stem's name, and a base letter, at these sizes" (a copy of render's); earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | - | - |
| `layout:_LETTER_EM` | constant | 0.75 | _NAME_EM, _LETTER_EM = 0.8, 0.75 # ems: render draws a stem's name, and a base letter, at these sizes | geometry | "render draws a stem's name, and a base letter, at these sizes" (a copy of render's); earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | na_curve, na_layout, na_text_keepoff | - |
| `layout:DNA_TUBE` | constant | 0.42 | DNA_TUBE = 0.42 # a drawn double helix's backbone tube width, page units (render's ribbon) | unknown | earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | - | - |
| `layout:_FALLBACK_WIDTH` | constant | 0.6 | _FALLBACK_WIDTH = 0.6 # ems: a character the metric table lacks (about a bold lower-case letter) | unknown | earliest in this file: af1d469 2026-10-11 Stem names are measured from built-in Liberation Sans Bold metrics, not installe… | - | - |
| `layout:RUN_EDGE_W` | constant | 1.8 | RUN_EDGE_W = 1.8 # render: the dark edge under residue-coloured loops, in loop widths | unknown | earliest in this file: c45f456 2026-10-16 Counts and helix names test a line as A6 reads it: one exact segment test at hal… | - | - |
| `layout:NA_RUN_W` | constant | 1.4 | NA_RUN_W = 1.4 # render: a nucleic-acid strand's run, in line_pt | unknown | earliest in this file: c45f456 2026-10-16 Counts and helix names test a line as A6 reads it: one exact segment test at hal… | na_text_keepoff | - |
| `layout:_OVER_STUB` | constant | 0.6 | _OVER_STUB = 0.6 # page units: the straight run of a loop out of a strand's end, kept clear of a name over the sheet | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `layout:_OVER_ALONG` | constant | (0.0, -1 / 3, 1 / 3) | _OVER_ALONG = (0.0, -1 / 3, 1 / 3) # where a name over a sheet may sit along its helix: helix lengths from its middle | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `layout:COUNT_ROOM` | constant | 1.0 | COUNT_ROOM = 1.0 # page units a count may stand past the layout's bounds and lines (render's page margin, `pad`) | unknown | earliest in this file: c45f456 2026-10-16 Counts and helix names test a line as A6 reads it: one exact segment test at hal… | - | - |
| `layout:_EM_TUNED` | constant | 0.55 * SHAFT | _EM_TUNED = 0.55 * SHAFT # page units per em of the figure's text (before Style.font_scale), on every page | unknown | earliest in this file: a3cf24d 2026-10-15 C4 (random-100, batch C1): text has a 6 pt floor and the page grows instead; leg… | heldout_budgets, page_text | - |
| `layout:_MIN_SCALE` | constant | MIN_TEXT_PT / (_EM_TUNED * 72) | _MIN_SCALE = MIN_TEXT_PT / (_EM_TUNED * 72) # inches per page unit at which the text is MIN_TEXT_PT | geometry | MIN_TEXT_PT / (_EM_TUNED * 72); earliest in this file: a3cf24d 2026-10-15 C4 (random-100, batch C1): text has a 6 pt floor and the page grows instead; leg… | - | - |
| `layout:MPL_MIN_PT` | constant | 1.0 | MPL_MIN_PT = 1.0 # matplotlib draws no text smaller than this (FontProperties.set_size floors a size under 1 pt) | convention | "matplotlib draws no text smaller than this"; earliest in this file: bc0e375 2026-10-17 Held-out domain names and loops (plan 15 task 2): names and leaders off the ways… | heldout_budgets, page_text | - |
| `layout:SMALLEST_TEXT` | constant | 0.6 | SMALLEST_TEXT = 0.6 # ems of the figure's text: the smallest render draws (a ligand's name, a residue number) | unknown | earliest in this file: bc0e375 2026-10-17 Held-out domain names and loops (plan 15 task 2): names and leaders off the ways… | heldout_budgets, page_text | - |
| `layout:_FLOOR_SCALE` | constant | MPL_MIN_PT / (SMALLEST_TEXT * _EM_TUNED * 72) | _FLOOR_SCALE = MPL_MIN_PT / (SMALLEST_TEXT * _EM_TUNED * 72) # inches per page unit at which that is MPL_MIN_PT | geometry | MPL_MIN_PT / (SMALLEST_TEXT * _EM_TUNED * 72); earliest in this file: bc0e375 2026-10-17 Held-out domain names and loops (plan 15 task 2): names and leaders off the ways… | heldout_budgets | - |
| `layout:MAX_GROWN_IN` | constant | 200.0 | MAX_GROWN_IN = 200.0 # the longest side that floor grows a page to: PDF readers open no page past 200 in (14 400 pt; | convention | "PDF readers open no page past 200 in"; earliest in this file: 1bf5f36 2026-10-18 Bounded work on very large pages (plan 15 task 7, A1): disulfide second starts s… | heldout_budgets | - |
| `layout:NA_TORN_TURN` | constant | 90.0 | NA_TORN_TURN = 90.0 # degrees: linked stems turned against one another by more than this when unrolled are torn | unknown | earliest in this file: 7baebdb 2026-10-11 Move _na_torn_by_unrolling to layout.py, on _components and _SHORT_RUN | - | - |
| `layout:NA_HELD` | constant | 3 | NA_HELD = 3 # residues of one protomer touching a stem piece for that protomer to hold it (as cli.NA_HOLD_MIN) | geometry | "as cli.NA_HOLD_MIN" (a copy); earliest in this file: 3be482a 2026-10-15 C16 (random-100, batch C2a): a nucleic acid most of whose stems are held by two … | - | - |
| `layout:PlacedStem.rail:argument:4@be533680` | argument | 4 | return max(self.width / 2 - RAIL_IN, self.width / 4) | unknown | - | ligands, na_contacts, na_curve, na_layout +2 | - |
| `layout:PlacedStem.hairpin_circle:argument:1.02@a34d9b16` | argument | 1.02 | r = max(half * 1.02, (m + 1) * NT_STEP / (1.6 * np.pi)) | unknown | - | na_layout | - |
| `layout:PlacedStem.hairpin_circle:argument:1.6@a34d9b16` | argument | 1.6 | r = max(half * 1.02, (m + 1) * NT_STEP / (1.6 * np.pi)) | unknown | - | na_layout | - |
| `layout:PlacedStem.hairpin_circle:argument:0.01@cfc56112` | argument | 0.01 | r = float(_snap(r + 0.01)) | unknown | line last changed by ed8d7f2 A hairpin loop's circle grows until its nucleotides stand at least a contact bead apart (C2b-4 cause B: 7Y7C … (git blame; evidence only) | na_layout | - |
| `layout:_sheet_order:comparison:1e-6@d2a4ba49` | comparison | 1e-6 | right = np.array([up2d[1], -up2d[0]]) if np.linalg.norm(up2d) > 1e-6 else np.array([1.0, 0.0]) | convention | a floating-point tolerance | copies | - |
| `layout:_stack_gaps:comparison:1e-6@bab5b4aa` | comparison | 1e-6 | above = [q for q in strands if q is not p and abs(q.cx - p.cx) < 1e-6 and q.cy > p.cy] | convention | a floating-point tolerance | - | - |
| `layout:_column_heights.way_out:comparison:1e-9@b6c36b98` | comparison | 1e-9 | if min(hi, b) - max(lo, a) < need - 1e-9: | convention | a floating-point tolerance | - | - |
| `layout:_column_heights.shut:comparison:1e-9@f92b269f` | comparison | 1e-9 | if 0 < (y - half) - top <= room(q, k) + _STACK_NUDGE + 1e-9: | convention | a floating-point tolerance | route | - |
| `layout:_column_heights.shut:comparison:1e-9@9513ea38` | comparison | 1e-9 | elif 0 < bottom - (y + half) <= room(k, q) + _STACK_NUDGE + 1e-9: | convention | a floating-point tolerance | route | - |
| `layout:_column_heights.meets:comparison:1e-9@14fab0f5` | comparison | 1e-9 | height[q] - length(q) / 2 - room(k, q) - half + 1e-9 | convention | a floating-point tolerance | app, domains, elstyle, links +7 | - |
| `layout:_column_heights.meets:comparison:1e-9@31e056c7` | comparison | 1e-9 | < height[q] + length(q) / 2 + room(q, k) + half - 1e-9 | convention | a floating-point tolerance | app, domains, elstyle, links +7 | - |
| `layout:_open_at_ends:comparison:3@8f5bc337` | comparison | 3 | if not sheet.closed or n < 3 or len({s.chain for s in sheet.strands}) != 1: | unknown | - | template | - |
| `layout:can_be_row:comparison:4@66f511f3` | comparison | 4 | return 4 <= len(sheet.strands) <= _BA_MAX and not sheet.ambiguous and all(k == "P" for k in sheet.pair_kinds) | unknown | - | template | - |
| `layout:_kink_item:comparison:-0.5@5e4070af` | comparison | -0.5 | if turn < -0.5: # turned back: a hairpin beside the last piece, on the side the piece turns to | unknown | - | - | - |
| `layout:_kink_item:comparison:1e-9@326ce42a` | comparison | 1e-9 | if back > 1e-9: | convention | a floating-point tolerance | - | - |
| `layout:_bundle_item:argument:1.2@2c76c149` | argument | 1.2 | spacing = max(width + 1.2, float(np.mean(gaps)) * SCALE if gaps else 0.0) | unknown | - | layout | - |
| `layout:_overlapping:comparison:1e-9@76c6fe72` | comparison | 1e-9 | hit = (depth[..., 0] > 1e-9) & (depth[..., 1] > 1e-9) | convention | a floating-point tolerance | membrane, perf, symmetry | - |
| `layout:_overlapping:comparison:1e-9@76c6fe72#2` | comparison | 1e-9 | hit = (depth[..., 0] > 1e-9) & (depth[..., 1] > 1e-9) | convention | a floating-point tolerance | membrane, perf, symmetry | - |
| `layout:_overlapping:comparison:1e-9@b88670a4` | comparison | 1e-9 | hit = (dx > 1e-9) & (dy > 1e-9) | convention | a floating-point tolerance | membrane, perf, symmetry | - |
| `layout:_overlapping:comparison:1e-9@b88670a4#2` | comparison | 1e-9 | hit = (dx > 1e-9) & (dy > 1e-9) | convention | a floating-point tolerance | membrane, perf, symmetry | - |
| `layout:_band_out:comparison:1e-9@d606e0be` | comparison | 1e-9 | if abs(s) >= across - 1e-9 or u - along >= ahead - 1e-9 or u + along <= -back + 1e-9: | convention | a floating-point tolerance | na_curve | - |
| `layout:_band_out:comparison:1e-9@d606e0be#2` | comparison | 1e-9 | if abs(s) >= across - 1e-9 or u - along >= ahead - 1e-9 or u + along <= -back + 1e-9: | convention | a floating-point tolerance | na_curve | - |
| `layout:_band_out:comparison:1e-9@d606e0be#3` | comparison | 1e-9 | if abs(s) >= across - 1e-9 or u - along >= ahead - 1e-9 or u + along <= -back + 1e-9: | convention | a floating-point tolerance | na_curve | - |
| `layout:_band_out:comparison:1e-6@79dcdf0a` | comparison | 1e-6 | ways = [(g / abs(v[0]), np.array([np.sign(v[0]), 0.0])) for g, v in ways if abs(v[0]) > 1e-6] | convention | a floating-point tolerance | na_curve | - |
| `layout:_push_apart:comparison:100@bd70c568` | comparison | 100 | if not len(pairs) and not _off_walls(xy, halves, mids, walls, group, move=False) and sweep >= 100: | unknown | - | heldout_copies, symmetry, symmetry_scenarios | - |
| `layout:_push_apart:comparison:1e-9@15c144c6` | comparison | 1e-9 | if px <= 1e-9 or py <= 1e-9: # an earlier push this sweep already cleared it | convention | a floating-point tolerance | heldout_copies, symmetry, symmetry_scenarios | - |
| `layout:_push_apart:comparison:1e-9@15c144c6#2` | comparison | 1e-9 | if px <= 1e-9 or py <= 1e-9: # an earlier push this sweep already cleared it | convention | a floating-point tolerance | heldout_copies, symmetry, symmetry_scenarios | - |
| `layout:_push_apart:comparison:100@205197c3` | comparison | 100 | if sweep < 100: | unknown | - | heldout_copies, symmetry, symmetry_scenarios | - |
| `layout:_off_walls:comparison:1e-9@5d0a7078` | comparison | 1e-9 | hits = np.argwhere((pen[..., 0] > 1e-9) & (pen[..., 1] > 1e-9) & walls.hard) | convention | a floating-point tolerance | membrane, na_curve | - |
| `layout:_off_walls:comparison:1e-9@5d0a7078#2` | comparison | 1e-9 | hits = np.argwhere((pen[..., 0] > 1e-9) & (pen[..., 1] > 1e-9) & walls.hard) | convention | a floating-point tolerance | membrane, na_curve | - |
| `layout:_off_walls:comparison:1e-9@15c144c6` | comparison | 1e-9 | if px <= 1e-9 or py <= 1e-9: | convention | a floating-point tolerance | membrane, na_curve | - |
| `layout:_off_walls:comparison:1e-9@15c144c6#2` | comparison | 1e-9 | if px <= 1e-9 or py <= 1e-9: | convention | a floating-point tolerance | membrane, na_curve | - |
| `layout:_slide_clear.meets:comparison:4@05f1e9ad` | comparison | 4 | if len(o) < 4: | unknown | - | app, domains, elstyle, links +7 | - |
| `layout:_mass:argument:4@3d2e2b08` | argument | 4 | return sum(len(m.axial) / 4 if isinstance(m, PlacedStem) else 1 for m in it.members) | unknown | - | - | - |
| `layout:_energy:argument:1e-12@cbf7475b` | argument | 1e-12 | dist = np.maximum(np.sqrt(g[:, 0] * g[:, 0] + g[:, 1] * g[:, 1]), 1e-12) | convention | a floating-point tolerance | optimiser, perf | - |
| `layout:_energy_dense:argument:1e-12@cbf7475b` | argument | 1e-12 | dist = np.maximum(np.sqrt(g[:, 0] * g[:, 0] + g[:, 1] * g[:, 1]), 1e-12) | convention | a floating-point tolerance | na_layout, optimiser, perf | - |
| `layout:_lbfgs:default:maxiter` | default | 300 | def _lbfgs(fun, x0, maxiter=300, m=10, pgtol=1e-5, ftol=2.2e-9): | unknown | - | heldout_budgets, optimiser, symmetry_scenarios | - |
| `layout:_lbfgs:default:m` | default | 10 | def _lbfgs(fun, x0, maxiter=300, m=10, pgtol=1e-5, ftol=2.2e-9): | unknown | - | heldout_budgets, optimiser, symmetry_scenarios | - |
| `layout:_lbfgs:default:pgtol` | default | 1e-5 | def _lbfgs(fun, x0, maxiter=300, m=10, pgtol=1e-5, ftol=2.2e-9): | convention | a floating-point tolerance | heldout_budgets, optimiser, symmetry_scenarios | - |
| `layout:_lbfgs:default:ftol` | default | 2.2e-9 | def _lbfgs(fun, x0, maxiter=300, m=10, pgtol=1e-5, ftol=2.2e-9): | convention | a floating-point tolerance | heldout_budgets, optimiser, symmetry_scenarios | - |
| `layout:_lbfgs:argument:1e10@19805495` | argument | 1e10 | step = 1.0 if S else min(1.0 / np.sqrt((d * d).sum()), 1e10) | unknown | - | heldout_budgets, optimiser, symmetry_scenarios | - |
| `layout:_wolfe:default:c1` | default | 1e-3 | def _wolfe(fun, x, f0, g0, d, step, c1=1e-3, c2=0.9, evals=20): | unknown | - | optimiser | - |
| `layout:_wolfe:default:c2` | default | 0.9 | def _wolfe(fun, x, f0, g0, d, step, c1=1e-3, c2=0.9, evals=20): | unknown | - | optimiser | - |
| `layout:_wolfe:default:evals` | default | 20 | def _wolfe(fun, x, f0, g0, d, step, c1=1e-3, c2=0.9, evals=20): | unknown | - | optimiser | - |
| `layout:_wolfe:comparison:0.1@13289bf0` | comparison | 0.1 | if min(a, b) + 0.1 * abs(b - a) <= c <= max(a, b) - 0.1 * abs(b - a): | unknown | - | optimiser | - |
| `layout:_wolfe:comparison:0.1@13289bf0#2` | comparison | 0.1 | if min(a, b) + 0.1 * abs(b - a) <= c <= max(a, b) - 0.1 * abs(b - a): | unknown | - | optimiser | - |
| `layout:_attract:argument:4@a26058a4` | argument | 4 | boxes = np.asarray(fixed.boxes if fixed is not None else [], float).reshape(-1, 4) | convention | an array's shape | copies, heldout_budgets, like, membrane +2 | - |
| `layout:_room_for_stems:argument:4.0@c8dc89f9` | argument | 4.0 | need = float(sum(4.0 * float(x) * float(y) for x, y in h.tolist())) | unknown | line last changed by ea0765b A nucleic acid drawn without protein is spread by one factor to the room its stems take before the pull, keep… (git blame; evidence only) | na_layout | - |
| `layout:_collide:comparison:1e-9@511e5cc6` | comparison | 1e-9 | return ha[0] + hb[0] + MARGIN - abs(dx) > 1e-9 and ha[1] + hb[1] + MARGIN - abs(dy) > 1e-9 | convention | a floating-point tolerance | heldout_copies, layout, perf, symmetry | - |
| `layout:_collide:comparison:1e-9@511e5cc6#2` | comparison | 1e-9 | return ha[0] + hb[0] + MARGIN - abs(dx) > 1e-9 and ha[1] + hb[1] + MARGIN - abs(dy) > 1e-9 | convention | a floating-point tolerance | heldout_copies, layout, perf, symmetry | - |
| `layout:_stem_item:comparison:1e-6@2caeff4e` | comparison | 1e-6 | angle = 0.0 if flat or np.linalg.norm(a2) < 1e-6 else _half_turn(_snap(np.arctan2(a2[1], a2[0]))) | convention | a floating-point tolerance | - | - |
| `layout:_trace_walls:argument:4@cf5b02a8` | argument | 4 | return Walls(np.asarray(boxes, float).reshape(-1, 4), hard) | convention | an array's shape | na_layout | - |
| `layout:_traces:comparison:1e-9@e6d4e716` | comparison | 1e-9 | if float(np.hypot(*along)) < 1e-9: # its neighbours on one spot: its 5' neighbour's way, or across | convention | a floating-point tolerance | na_layout | - |
| `layout:_traces:comparison:1e-9@e6d4e716#2` | comparison | 1e-9 | if float(np.hypot(*along)) < 1e-9: | convention | a floating-point tolerance | na_layout | - |
| `layout:_drawn_wrong:argument:3@32223231` | argument | 3 | spots.append((m.sse.id, round(float(ox), 3), (oy - half, oy + half), m.ghost, m.width)) | convention | rounding digits | - | - |
| `layout:_partners:comparison:1e-6@72b46712` | comparison | 1e-6 | and 1e-6 < abs(it.offsets[q][0] - x) <= max(mp.width, it.members[q].width) + 1e-6 | convention | a floating-point tolerance | - | - |
| `layout:_partners:comparison:1e-6@72b46712#2` | comparison | 1e-6 | and 1e-6 < abs(it.offsets[q][0] - x) <= max(mp.width, it.members[q].width) + 1e-6 | convention | a floating-point tolerance | - | - |
| `layout:_widen:comparison:0.99@673631b5` | comparison | 0.99 | upright = all(abs(np.sin(it.members[k].angle)) > 0.99 for k in side) | unknown | - | membrane | - |
| `layout:_widen:comparison:0.99@63ed3cc4` | comparison | 0.99 | lying = all(abs(np.cos(it.members[k].angle)) > 0.99 for k in side) | unknown | - | membrane | - |
| `layout:_widen:argument:6@af80b3a1` | argument | 6 | at = sorted({round(it.offsets[k][axis], 6) for k in side}) | convention | rounding digits | membrane | - |
| `layout:_widen:argument:6@71780f66` | argument | 6 | grown = {c: max(extra[k] for k in side if round(it.offsets[k][axis], 6) == c) for c in at} | convention | rounding digits | membrane | - |
| `layout:_widen:comparison:1e-6@943dd31f` | comparison | 1e-6 | if over_x <= 1e-6 or over_y <= 1e-6 or met: # touching is not meeting (a sheet's pitch) | convention | a floating-point tolerance | membrane | - |
| `layout:_widen:comparison:1e-6@943dd31f#2` | comparison | 1e-6 | if over_x <= 1e-6 or over_y <= 1e-6 or met: # touching is not meeting (a sheet's pitch) | convention | a floating-point tolerance | membrane | - |
| `layout:_into_slab:comparison:1e-9@d11bd627` | comparison | 1e-9 | if rise > 1e-9: | convention | a floating-point tolerance | membrane | - |
| `layout:_pack_run:comparison:1e-9@48949366` | comparison | 1e-9 | if least > most + 1e-9: | convention | a floating-point tolerance | - | - |
| `layout:_in_order:comparison:1e-9@f64385d2` | comparison | 1e-9 | return any(abs(d) > 1e-9 for d in moved.values()) | convention | a floating-point tolerance | membrane | - |
| `layout:_parted_rooms:comparison:1e-9@e39a2fe0` | comparison | 1e-9 | if abs(d) > 1e-9: | convention | a floating-point tolerance | membrane | - |
| `layout:_out_of_room:comparison:1e-9@401054c9` | comparison | 1e-9 | if not lo - 1e-9 <= items[k].cx <= hi + 1e-9: | convention | a floating-point tolerance | membrane | - |
| `layout:_out_of_room:comparison:1e-9@401054c9#2` | comparison | 1e-9 | if not lo - 1e-9 <= items[k].cx <= hi + 1e-9: | convention | a floating-point tolerance | membrane | - |
| `layout:_stretch_short:comparison:1e-9@71a6a0cd` | comparison | 1e-9 | and abs(abs(float(np.sin(m.angle))) - 1) < 1e-9 | convention | a floating-point tolerance | membrane | - |
| `layout:_stretch_short:comparison:1e-6@8aa0fc36` | comparison | 1e-6 | and m.length <= half + 1e-6 | convention | a floating-point tolerance | membrane | - |
| `layout:_rosette_push:comparison:-1e-9@64589916` | comparison | -1e-9 | if all((cdist(q, _snap(turn2(q, a))) - reach).min() >= -1e-9 for a in angles[1:]): | convention | a floating-point tolerance | - | - |
| `layout:_box_on_body:comparison:1e-6@194b13b8` | comparison | 1e-6 | if min(a.max(), b.max()) - max(a.min(), b.min()) <= 1e-6: | convention | a floating-point tolerance | rosette | - |
| `layout:_wrap_angle:comparison:1e-9@17005f7c` | comparison | 1e-9 | return float(_HALF_TURN) if a <= -_HALF_TURN + 1e-9 else a + 0.0 | convention | a floating-point tolerance | - | - |
| `layout:_rosette:comparison:1e-6@42142d6c` | comparison | 1e-6 | out = _snap(mid / np.linalg.norm(mid)) if np.linalg.norm(mid) > 1e-6 else np.array([0.0, 1.0]) | convention | a floating-point tolerance | - | - |
| `layout:_rosette:comparison:1e-6@a1fb05c6` | comparison | 1e-6 | radial = radial / np.linalg.norm(radial) if np.linalg.norm(radial) > 1e-6 else turn2(out, a) | convention | a floating-point tolerance | - | - |
| `layout:_rosette:argument:4@4100ee4e` | argument | 4 | rects = np.array([p.rect for p in [*lay.placed.values(), *lay.ghosts]]).reshape(-1, 4) | convention | an array's shape | - | - |
| `layout:_rosette:argument:4@2b642db6` | argument | 4 | rects = np.vstack([rects, [b for _, b in label_boxes(lay, sses)] or np.empty((0, 4))]) | convention | an array's shape | - | - |
| `layout:build_layout:argument:4@572130e5` | argument | 4 | rects = np.array([p.rect for p in placed + ghosts + dna]).reshape(-1, 4) | convention | an array's shape | cli, heldout_copies, layout, like +10 | - |
| `layout:build_layout:argument:4@2b642db6` | argument | 4 | rects = np.vstack([rects, [b for _, b in label_boxes(lay, sses)] or np.empty((0, 4))]) | convention | an array's shape | cli, heldout_copies, layout, like +10 | - |
| `layout:_part_domains:default:sweeps` | default | 60 | sweeps: int = 60, | unknown | - | copies, domains, membrane, perf | - |
| `layout:_part_sweeps:comparison:1e-6@d7181a85` | comparison | 1e-6 | if now >= last - 1e-6: # the moves no longer help | convention | a floating-point tolerance | domains | - |
| `layout:_level_beside:comparison:1e-6@040360a0` | comparison | 1e-6 | beside = 1e-6 < abs(x.cx - p.cx) <= max(x.width, p.width) + 1e-6 | convention | a floating-point tolerance | - | - |
| `layout:_level_beside:comparison:1e-6@040360a0#2` | comparison | 1e-6 | beside = 1e-6 < abs(x.cx - p.cx) <= max(x.width, p.width) + 1e-6 | convention | a floating-point tolerance | - | - |
| `layout:_level_beside:comparison:1e-6@7e316208` | comparison | 1e-6 | return beside and min(x1, p1) - max(x0, p0) > 1e-6 | convention | a floating-point tolerance | - | - |
| `layout:_ways_out:argument:0.2@f24e3306` | argument | 0.2 | ways.append((min(ax, bx) - 0.2, min(ay, by) - 0.2, max(ax, bx) + 0.2, max(ay, by) + 0.2)) | unknown | - | domains | - |
| `layout:_ways_out:argument:0.2@f24e3306#2` | argument | 0.2 | ways.append((min(ax, bx) - 0.2, min(ay, by) - 0.2, max(ax, bx) + 0.2, max(ay, by) + 0.2)) | unknown | - | domains | - |
| `layout:_ways_out:argument:0.2@f24e3306#3` | argument | 0.2 | ways.append((min(ax, bx) - 0.2, min(ay, by) - 0.2, max(ax, bx) + 0.2, max(ay, by) + 0.2)) | unknown | - | domains | - |
| `layout:_ways_out:argument:0.2@f24e3306#4` | argument | 0.2 | ways.append((min(ax, bx) - 0.2, min(ay, by) - 0.2, max(ax, bx) + 0.2, max(ay, by) + 0.2)) | unknown | - | domains | - |
| `layout:_stem_ways:argument:0.4@4a988c10` | argument | 0.4 | out.append(_Way(float(c[0]), float(c[1]), d.direction, d.normal, _NAME_STUB, d.width + 0.4)) | unknown | - | domains | - |
| `layout:_name_obstacles:argument:0.2@f8c827df` | argument | 0.2 | ways.append((*(ends.min(axis=0) - 0.2), *(ends.max(axis=0) + 0.2))) | unknown | - | domains | - |
| `layout:_name_obstacles:argument:0.2@f8c827df#2` | argument | 0.2 | ways.append((*(ends.min(axis=0) - 0.2), *(ends.max(axis=0) + 0.2))) | unknown | - | domains | - |
| `layout:_place_domain_names.clear:argument:0.05@3961c6c1` | argument | 0.05 | if overlaps_any((x - 0.05, min(ya, yb), x + 0.05, max(ya, yb)), solid): | unknown | - | - | - |
| `layout:_place_domain_names.clear:argument:0.05@3961c6c1#2` | argument | 0.05 | if overlaps_any((x - 0.05, min(ya, yb), x + 0.05, max(ya, yb)), solid): | unknown | - | - | - |
| `layout:_place_domain_names.clear:argument:0.1@a0ebc66f` | argument | 0.1 | keep = _grown((x - 0.1, min(ya, yb), x + 0.1, max(ya, yb)), max(margin, _NAME_OPEN - _NAME_ROUTE)) | unknown | - | - | - |
| `layout:_place_domain_names.clear:argument:0.1@a0ebc66f#2` | argument | 0.1 | keep = _grown((x - 0.1, min(ya, yb), x + 0.1, max(ya, yb)), max(margin, _NAME_OPEN - _NAME_ROUTE)) | unknown | - | - | - |
| `layout:_place_domain_names.room:argument:0.1@c732ebcc` | argument | 0.1 | keep.append((x - 0.1, min(ya, yb), x + 0.1, max(ya, yb))) | unknown | - | - | - |
| `layout:_place_domain_names.room:argument:0.1@c732ebcc#2` | argument | 0.1 | keep.append((x - 0.1, min(ya, yb), x + 0.1, max(ya, yb))) | unknown | - | - | - |
| `layout:_place_domain_names:argument:0.1@14901bcb` | argument | 0.1 | spots = list(_name_spots(value > 0, fields.xs, fields.ys, size * _text_width(name) + 0.1, size * _NAME_HEIGHT)) | unknown | - | domains | - |
| `layout:_place_domain_names:fallback:domain_name_least_bad` | fallback | domain_name_least_bad | fallback.note("domain_name_least_bad") | named-entry:7NZ2 | docstring: "the one that leaves the ways out the most (random-100 held-out 7NZ2: a name falling back to the top left of a panel round one helix stood on its way out)" | - | A16 judges domain names clear (gated) |
| `layout:_name_spots:argument:0.05@1e78927a` | argument | 0.05 | ends = (max(sx.min() + _NAME_INSET, box[0] + 0.05), min(sx.max() - _NAME_INSET, box[2] - 0.05)) | unknown | - | domains | - |
| `layout:_name_spots:argument:0.05@1e78927a#2` | argument | 0.05 | ends = (max(sx.min() + _NAME_INSET, box[0] + 0.05), min(sx.max() - _NAME_INSET, box[2] - 0.05)) | unknown | - | domains | - |
| `layout:_name_spots:comparison:1e-6@421ef3f6` | comparison | 1e-6 | if abs(x - mid) > 1e-6 and ends[0] <= ends[1]: | convention | a floating-point tolerance | domains | - |
| `layout:_ink_line:fallback:glyph_width_guessed` | fallback | glyph_width_guessed | fallback.note("glyph_width_guessed") | unknown | _FALLBACK_WIDTH: "a character the metric table lacks (about a bold lower-case letter)" | - | none |
| `layout:_reach:argument:4@4ae348b6` | argument | 4 | return np.zeros(4) | convention | an array's shape | na_curve, na_layout, na_text_keepoff, text_clear | - |
| `layout:upright_deg:comparison:90@c1899080` | comparison | 90 | return rot - 180 if rot > 90 else rot + 180 if rot < -90 else rot | unknown | - | - | - |
| `layout:upright_deg:comparison:-90@c1899080` | comparison | -90 | return rot - 180 if rot > 90 else rot + 180 if rot < -90 else rot | unknown | - | - | - |
| `layout:_text_width:fallback:glyph_width_guessed` | fallback | glyph_width_guessed | fallback.note("glyph_width_guessed") | unknown | _FALLBACK_WIDTH: "a character the metric table lacks (about a bold lower-case letter)" | - | none |
| `layout:line_pt:argument:0.9@694aa4f3` | argument | 0.9 | return float(np.clip(0.9 * scale / IN_PER_UNIT, 0.5, 1.2)) | unknown | - | na_text_keepoff, text_clear | - |
| `layout:line_pt:argument:1.2@694aa4f3` | argument | 1.2 | return float(np.clip(0.9 * scale / IN_PER_UNIT, 0.5, 1.2)) | unknown | - | na_text_keepoff, text_clear | - |
| `layout:line_band:argument:72@7b85cc51` | argument | 72 | return (line_pt(1.0 / (72 * unit)) * stroke / 2 + 0.72) * unit | convention | 72 points to the inch | na_text_keepoff, text_clear | - |
| `layout:lines_meeting:argument:4@dcbfb270` | argument | 4 | box = np.asarray(boxes, float).reshape(-1, 4) + np.array([-band, -band, band, band]) | convention | an array's shape | na_text_keepoff, text_clear | - |
| `layout:lines_meeting:argument:4@a9e132e9` | argument | 4 | seg = np.asarray(segments, float).reshape(-1, 4) | convention | an array's shape | na_text_keepoff, text_clear | - |
| `layout:lines_meeting:comparison:1e-9@a8168bc1` | comparison | 1e-9 | out[near] = clip_lengths(seg[near], box) > 1e-9 | convention | a floating-point tolerance | na_text_keepoff, text_clear | - |
| `layout:spans_crossed:argument:4@f13a4ba9` | argument | 4 | span = np.asarray(spans, float).reshape(-1, 4) | convention | an array's shape | text_clear | - |
| `layout:spans_crossed:argument:4@a9e132e9` | argument | 4 | seg = np.asarray(segments, float).reshape(-1, 4) | convention | an array's shape | text_clear | - |
| `layout:spans_crossed:comparison:1e-9@0e109098` | comparison | 1e-9 | if length < 1e-9: | convention | a floating-point tolerance | text_clear | - |
| `layout:spans_crossed:argument:4@5de955d1` | argument | 4 | out[k] = lines_meeting((0.0, 0.0, length, 0.0), flat.reshape(-1, 4), 1e-6).any() | convention | an array's shape | text_clear | - |
| `layout:spans_crossed:argument:1e-6@5de955d1` | argument | 1e-6 | out[k] = lines_meeting((0.0, 0.0, length, 0.0), flat.reshape(-1, 4), 1e-6).any() | convention | a floating-point tolerance | text_clear | - |
| `layout:spans_crossed:argument:4@e15f9be9` | argument | 4 | rect = np.asarray(bodies, float).reshape(-1, 4) | convention | an array's shape | text_clear | - |
| `layout:spans_crossed:comparison:1e-9@11646f24` | comparison | 1e-9 | out \|= (clip_lengths(span, rect) > 1e-9).any(axis=1) | convention | a floating-point tolerance | text_clear | - |
| `layout:segments_of:argument:4@69329d44` | argument | 4 | return np.asarray(segs, float).reshape(-1, 4) | convention | an array's shape | na_text_keepoff, text_clear | - |
| `layout:_over_label_sides:argument:0.05@cf299152` | argument | 0.05 | letters = [centred_box(q.cx, q.cy, gw * len(q.label) / 2 + 0.05, gh / 2) for q in strands] | unknown | - | - | - |
| `layout:_over_label_sides:argument:0.15@e9c0e7e0` | argument | 0.15 | ends.append((min(ax, bx) - 0.15, min(ay, by) - 0.15, max(ax, bx) + 0.15, max(ay, by) + 0.15)) | unknown | - | - | - |
| `layout:_over_label_sides:argument:0.15@e9c0e7e0#2` | argument | 0.15 | ends.append((min(ax, bx) - 0.15, min(ay, by) - 0.15, max(ax, bx) + 0.15, max(ay, by) + 0.15)) | unknown | - | - | - |
| `layout:_over_label_sides:argument:0.15@e9c0e7e0#3` | argument | 0.15 | ends.append((min(ax, bx) - 0.15, min(ay, by) - 0.15, max(ax, bx) + 0.15, max(ay, by) + 0.15)) | unknown | - | - | - |
| `layout:_over_label_sides:argument:0.15@e9c0e7e0#4` | argument | 0.15 | ends.append((min(ax, bx) - 0.15, min(ay, by) - 0.15, max(ax, bx) + 0.15, max(ay, by) + 0.15)) | unknown | - | - | - |
| `layout:_over_label_sides:argument:0.45@d4a53a5f` | argument | 0.45 | base = tip - q.direction * min(HEAD_LEN, q.length * 0.45) | unknown | - | - | - |
| `layout:_choose_label_sides:argument:4@3e779994` | argument | 4 | letters = {id(d): np.array([_reach(p) for _, p in letter_ink[id(d)]]).reshape(-1, 4) for d in layout.stems} | convention | an array's shape | names | - |
| `layout:_choose_label_sides:argument:4@f0ad813c` | argument | 4 | id(d): np.vstack([letters[id(e)] for e in layout.stems if e is not d] or [np.empty((0, 4))]) | convention | an array's shape | names | - |
| `layout:_choose_label_sides:argument:4@b16b056f` | argument | 4 | leaving = _leaving_stretches(layout) if block else np.empty((0, 4)) | convention | an array's shape | names | - |
| `layout:_choose_label_sides:fallback:stem_name_least_bad` | fallback | stem_name_least_bad | fallback.note("stem_name_least_bad") | unknown | _pick_slot: "else (a fallback, logged) the one with the most room" | - | an INFO log; A17 judges stem names by their letters (gated) |
| `layout:_leaving_stretches:argument:4@68aeb51c` | argument | 4 | return np.empty((0, 4)) | convention | an array's shape | na_layout | - |
| `layout:_leaving_stretches:argument:4@75150913` | argument | 4 | return np.asarray(out, float).reshape(-1, 4) | convention | an array's shape | na_layout | - |
| `layout:_run_band:argument:72@7b85cc51` | argument | 72 | return (line_pt(1.0 / (72 * unit)) * stroke / 2 + 0.72) * unit | convention | 72 points to the inch | na_text_keepoff | - |
| `layout:thin_boxes:default:pad` | default | 0.06 | def thin_boxes(points, pad: float = 0.06) -> list[tuple[float, float, float, float]]: | unknown | - | text_clear | - |
| `layout:count_half:default:pad` | default | 0.04 | def count_half(text: str, size: float, pad: float = 0.04) -> tuple[float, float]: | unknown | - | text_clear | - |
| `layout:count_page:argument:4@a9e132e9` | argument | 4 | seg = np.asarray(segments, float).reshape(-1, 4) | convention | an array's shape | - | - |
| `layout:clear_inks:argument:4@45e95a8f` | argument | 4 | ink = np.asarray(inks, float).reshape(-1, 4) | convention | an array's shape | text_clear | - |
| `layout:clear_inks:argument:4@c14dbcda` | argument | 4 | rect = np.asarray(boxes, float).reshape(-1, 4) | convention | an array's shape | text_clear | - |
| `layout:tail_counts:argument:0.05@8e83b09d` | argument | 0.05 | own = (max(gw / 2 + 0.05, em * _text_width(end) / 2), max(gh / 2, em * 0.4)) # its keep-out (label_boxes) | unknown | - | text_clear | - |
| `layout:tail_counts:argument:0.4@8e83b09d` | argument | 0.4 | own = (max(gw / 2 + 0.05, em * _text_width(end) / 2), max(gh / 2, em * 0.4)) # its keep-out (label_boxes) | unknown | - | text_clear | - |
| `layout:label_boxes:argument:0.15@47fade12` | argument | 0.15 | a, b = np.add(port, ex * 0.15), np.add(port, ex * END_STUB) | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.05@614c6f0d` | argument | 0.05 | (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05), | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.05@614c6f0d#2` | argument | 0.05 | (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05), | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.05@614c6f0d#3` | argument | 0.05 | (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05), | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.05@614c6f0d#4` | argument | 0.05 | (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05), | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.05@d0276357` | argument | 0.05 | out.append((f"{end}:{chain}", centred_box(float(t[0]), float(t[1]), gw / 2 + 0.05, gh / 2))) | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.05@14ad571b` | argument | 0.05 | out.append((f"tie-label:{n}:{end}", centred_box(x, y, gw / 2 + 0.05, gh / 2))) | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.1@877d99ef` | argument | 0.1 | out.append((f"domain-leader:{name}", (x - 0.1, min(y0, y1), x + 0.1, max(y0, y1)))) | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:label_boxes:argument:0.1@877d99ef#2` | argument | 0.1 | out.append((f"domain-leader:{name}", (x - 0.1, min(y0, y1), x + 0.1, max(y0, y1)))) | unknown | - | domains, heldout_budgets, layout, matrix +6 | - |
| `layout:font_pt:argument:72@02145b1a` | argument | 72 | return float(min(_EM_TUNED * scale * 72, 9.0)) | convention | 72 points to the inch | - | - |
| `layout:font_pt:argument:9.0@02145b1a` | argument | 9.0 | return float(min(_EM_TUNED * scale * 72, 9.0)) | unknown | - | - | - |
| `layout:page_fit:comparison:1e-6@76d04fa7` | comparison | 1e-6 | if font >= MIN_TEXT_PT - 1e-6: | convention | a floating-point tolerance | heldout_budgets, page_text | - |
| `layout:page_fit:comparison:1e-6@fbdc5ce9` | comparison | 1e-6 | if longest <= MAX_PAGE_IN + 1e-6: | convention | a floating-point tolerance | heldout_budgets, page_text | - |
| `layout:page_fit:comparison:1e-6@64e8f504` | comparison | 1e-6 | elif longest < MAX_GROWN_IN - 1e-6: | convention | a floating-point tolerance | heldout_budgets, page_text | - |
| `layout:_part_domains:fallback-marker` | fallback-marker | 1 | for it in items: # BLACT helix B:575-585, walled 6 units over its home, fell back to it alone, off its panel) | unknown | not a fallback: the comment tells of a past bug (32KS chain B's BLACT helix) | copies, domains, membrane, perf | - |
| `layout:_place_domain_names:fallback-marker` | fallback-marker | 1 | (else of all), the one that leaves the ways out the most (random-100 held-out 7NZ2: a name falling back to the | named-entry:7NZ2 | holds domain_name_least_bad | domains | - |
| `layout:_ink_line:fallback-marker` | fallback-marker | 1 | pen += advance.get(c, _FALLBACK_WIDTH * tm.UNITS) | unknown | holds glyph_width_guessed | na_layout, text_clear | - |
| `layout:<module>:fallback-marker` | fallback-marker | 1 | _FALLBACK_WIDTH = 0.6 # ems: a character the metric table lacks (about a bold lower-case letter) | unknown | _FALLBACK_WIDTH: the width glyph_width_guessed takes | - | - |
| `layout:_text_width:fallback-marker` | fallback-marker | 1 | guess = _FALLBACK_WIDTH * text_metrics.UNITS | unknown | holds glyph_width_guessed | domains, inventory, names | - |
| `layout:_pick_slot:fallback-marker` | fallback-marker | 2 | room; else (a fallback, logged) the one with the most room among the clear and nearer ones, else among the clear | unknown | the choice behind stem_name_least_bad (counted by its caller, _choose_label_sides) | na_layout | - |

## layoutfile

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `layoutfile:VERSION` | constant | 1 | VERSION = 1 | convention | a schema version; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `layoutfile:NAME_MAX` | constant | 80 | NAME_MAX = 80 # characters of an element's name (edits: rename) | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `layoutfile:REF_MAX` | constant | 200 | REF_MAX = 200 # characters of an element reference (res:A:3-18, α3, A:G) | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |
| `layoutfile:MOVE_MAX` | constant | 100.0 | MOVE_MAX = 100.0 # page units an edit may move an element by, each way (one strand spacing is about 1.1) | convention | an input limit; earliest in this file: 81aec46 2026-10-12 User YAML: no anchors or aliases in layout, names and style files; every layout … | - | - |

## like

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `like:EXTRA_TINT` | constant | 0.35 | EXTRA_TINT = 0.35 # an extra element's colour: its N-terminal neighbour's, mixed this much with white | unknown | earliest in this file: c2470a8 2026-10-21 A layout drawn like a reference (plan 13a Task 5): frame_like, Like, homes | - | - |
| `like:LOW_COVERAGE` | constant | 0.7 | LOW_COVERAGE = 0.7 # aligned residues under this fraction of the shorter chain: the note says how few aligned | unknown | earliest in this file: ea982cc 2026-10-21 --like on the command line: names, colours and the note (plan 13a Task 6) | - | - |

## mcp_server

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `mcp_server:PREVIEW_DPI` | constant | 90 | PREVIEW_DPI = 90 # the preview image's resolution ... | convention | a preview's resolution; earliest in this file: b613d0e 2026-10-15 C4 fix round 1 (finding 3): the MCP preview of a grown page is capped at 1600 px | - | - |
| `mcp_server:PREVIEW_PX` | constant | 1600 | PREVIEW_PX = 1600 # ... with its longer side at most this many pixels (a page grown for 6 pt text: 48 in) | convention | a preview's size; earliest in this file: b613d0e 2026-10-15 C4 fix round 1 (finding 3): the MCP preview of a grown page is capped at 1600 px | mcp | - |
| `mcp_server:draw_sequence:default:columns` | default | 60 | columns: int = 60, | unknown | - | mcp, uniprot | - |

## membrane

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `membrane:HALF` | constant | 15.0 | HALF = 15.0 # Å, half the hydrophobic thickness of a typical bilayer | biology | "half the hydrophobic thickness of a typical bilayer"; earliest in this file: 3bfd82c 2026-10-09 Add membrane bands | heldout_budgets, membrane | - |
| `membrane:_MIN_HYDROPHOBIC` | constant | 0.55 | _MIN_HYDROPHOBIC, _MAX_POLAR = 0.55, 0.15 # a slab this oily and this uncharged is membrane | unknown | earliest in this file: 3bfd82c 2026-10-09 Add membrane bands | - | - |
| `membrane:_MAX_POLAR` | constant | 0.15 | _MIN_HYDROPHOBIC, _MAX_POLAR = 0.55, 0.15 # a slab this oily and this uncharged is membrane | unknown | earliest in this file: 3bfd82c 2026-10-09 Add membrane bands | membrane | - |
| `membrane:_MAX_POLAR_CROSSED` | constant | 0.20 | _MAX_POLAR_CROSSED, _MANY_CROSSING = 0.20, 10 | named-entry:6JO6 | "... or this polar when this many long elements cross it (6JO6, photosystem I: 0.174 polar, 46 crossing; cofactor pockets in the slab; the sample's so…"; earliest in this file: 9d3f3aa 2026-10-15 C17 (random-100, batch C2a fix round 1): a slab many long elements cross may be … | membrane | - |
| `membrane:_MANY_CROSSING` | constant | 10 | _MAX_POLAR_CROSSED, _MANY_CROSSING = 0.20, 10 | named-entry:6JO6 | "... or this polar when this many long elements cross it (6JO6, photosystem I: 0.174 polar, 46 crossing; cofactor pockets in the slab; the sample's so…"; earliest in this file: 9d3f3aa 2026-10-15 C17 (random-100, batch C2a fix round 1): a slab many long elements cross may be … | membrane | - |
| `membrane:_SIDED_RATIO` | constant | 2.0 | _SIDED_RATIO, _SIDED_DIFF = 2.0, 3 # Lys/Arg counts must differ this clearly to call the sides | unknown | earliest in this file: 3bfd82c 2026-10-09 Add membrane bands | - | - |
| `membrane:_SIDED_DIFF` | constant | 3 | _SIDED_RATIO, _SIDED_DIFF = 2.0, 3 # Lys/Arg counts must differ this clearly to call the sides | unknown | earliest in this file: 3bfd82c 2026-10-09 Add membrane bands | - | - |
| `membrane:_MIN_CROSSING` | constant | 10.0 | _MIN_CROSSING = 10.0 # Å an element must reach beyond the centre on both sides to cross the membrane | unknown | earliest in this file: 3bfd82c 2026-10-09 Add membrane bands | - | - |
| `membrane:_FLANK` | constant | 20 | _FLANK, _FLANK_REACH = 20, 15.0 | measured | "Against OPM on the manuscript's membrane campaign (17 helical entries): 15 sided, 14 right" (32KS, 6QFA, 6ZKT); earliest in this file: 24ad3d2 2026-10-24 One crossing rule for the membrane view, A18 and data; sides from the loops (pla… | - | - |
| `membrane:_FLANK_REACH` | constant | 15.0 | _FLANK, _FLANK_REACH = 20, 15.0 | measured | "Against OPM on the manuscript's membrane campaign (17 helical entries): 15 sided, 14 right" (32KS, 6QFA, 6ZKT); earliest in this file: 24ad3d2 2026-10-24 One crossing rule for the membrane view, A18 and data; sides from the loops (pla… | - | - |
| `membrane:BARREL_HALF` | constant | 12.0 | BARREL_HALF = 12.0 # Å | biology | "OPM gives their hydrophobic half thickness as 11.7-12.6 Å (FimD 3RFZ 11.7, OmpX 1QJ8 11.8, OmpA 1BXW 12.6, OmpF 2OMF 12.6)"; earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:_BARREL_STRANDS` | constant | 8 | _BARREL_STRANDS, _BARREL_STRAND_LEN = 8, 7 # a closed sheet with this many strands this long is a barrel | unknown | earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:_BARREL_STRAND_LEN` | constant | 7 | _BARREL_STRANDS, _BARREL_STRAND_LEN = 8, 7 # a closed sheet with this many strands this long is a barrel | unknown | earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:_GIRDLE` | constant | 4.0 | _GIRDLE = 4.0 # Å either side of a membrane face where an outward Tyr/Trp counts as its aromatic girdle | unknown | earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:_MIN_OUTWARD_HYDROPHOBIC` | constant | 0.6 | _MIN_OUTWARD_HYDROPHOBIC, _MAX_OUTWARD_CHARGED = 0.6, 0.15 | measured | "Outward residues in the slab, measured: outer-membrane barrels 0.76-1.00 oily and 0.00-0.07 charged (9FWZ, 1BXW, 1QJ8, 2OMF, 2POR, 1QD6, 1THQ, 2F1C, …" (entries 9FWZ, 1BXW, 1QJ8, 2OMF, 2POR, 1QD6, 1THQ, 2F1C, 1I78, 1FEP, 2GSK, 3EMN, 7AHL, 1EK9, 1K24, 1TLY, 2MPR, 1QD5, 1P4T, 3PRN, 1EMA, 1GFL, 1RBP); earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | membrane | - |
| `membrane:_MAX_OUTWARD_CHARGED` | constant | 0.15 | _MIN_OUTWARD_HYDROPHOBIC, _MAX_OUTWARD_CHARGED = 0.6, 0.15 | measured | "Outward residues in the slab, measured: outer-membrane barrels 0.76-1.00 oily and 0.00-0.07 charged (9FWZ, 1BXW, 1QJ8, 2OMF, 2POR, 1QD6, 1THQ, 2F1C, …" (entries 9FWZ, 1BXW, 1QJ8, 2OMF, 2POR, 1QD6, 1THQ, 2F1C, 1I78, 1FEP, 2GSK, 3EMN, 7AHL, 1EK9, 1K24, 1TLY, 2MPR, 1QD5, 1P4T, 3PRN, 1EMA, 1GFL, 1RBP); earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:_LOOP_CAP` | constant | 20 | _LOOP_CAP, _SIDED_LOOPS = 20, 10 | measured | "Against OPM's sides (19 of the barrels above ...) ... call 17 right and 13 clearly"; earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:_SIDED_LOOPS` | constant | 10 | _LOOP_CAP, _SIDED_LOOPS = 20, 10 | measured | "Against OPM's sides (19 of the barrels above ...) ... call 17 right and 13 clearly"; earliest in this file: f12bfdd 2026-10-13 C-MEM-BARREL: outer-membrane beta-barrels get a membrane | - | - |
| `membrane:MIN_CROSSING_HELICES` | constant | 3 | MIN_CROSSING_HELICES = 3 | named-entry:7PIU,7VVK | "… soluble ones included, and against OPM it was 29.1° off on 7PIU and 32.2° on 7VVK (PTH1R-Gs: the G protein's helices; its normal 28.8° from the rec…"; earliest in this file: 43fd959 2026-10-25 Membrane normal refitted to the crossing helices, until they stop changing (plan… | membrane | - |
| `membrane:REFIT_ROUNDS` | constant | 8 | REFIT_ROUNDS = 8 | unknown | earliest in this file: 43fd959 2026-10-25 Membrane normal refitted to the crossing helices, until they stop changing (plan… | - | - |
| `membrane:FACE_REACH` | constant | 2.0 | FACE_REACH = 2.0 # Å: a long element whose CA run comes this near both faces spans the slab (_spans, the estimate) | unknown | earliest in this file: a963bb6 2026-10-18 Membrane (plan 15 task 6, A18): a slab no element spans moves to one they do, ju… | datafile, heldout_budgets, membrane | - |
| `membrane:SOLUBLE_BEYOND` | constant | 5.0 | SOLUBLE_BEYOND = 5.0 # Å: an element whose CA centre lies this far beyond a face is soluble, not interfacial (R4) | unknown | earliest in this file: 2e535a1 2026-10-24 Helices upright or flat by their tilt, soluble elements off the band (plan 26 Ta… | - | - |
| `membrane:FLAT_TILT` | constant | 45.0 | FLAT_TILT = ( | unknown | earliest in this file: 2e535a1 2026-10-24 Helices upright or flat by their tilt, soluble elements off the band (plan 26 Ta… | membrane | - |
| `membrane:FACE_FLAT_TILT` | constant | 35.0 | FACE_FLAT_TILT = 35.0 # degrees | named-entry:32KS | "…hite & Wimley 1999, Annu Rev Biophys Biomol Struct 28:319). 32KS's NS2 and NS3 (A:222-233, A:241-250), tilted 44.6° and 38.6° with their centres 2.9…"; earliest in this file: 102a6f2 2026-10-24 Membrane view refinements from the ComEC look test (plan 26 Task 5, rulings T5-1… | membrane | - |
| `membrane:INTERFACE` | constant | 7.5 | INTERFACE = 7.5 # Å either side of a face | biology | the bilayer's interfacial region (White & Wimley 1999, cited at FACE_FLAT_TILT); earliest in this file: 102a6f2 2026-10-24 Membrane view refinements from the ComEC look test (plan 26 Task 5, rulings T5-1… | membrane | - |
| `membrane:BAND_EPS` | constant | 1e-6 | BAND_EPS = 1e-6 # page units: float snapping, and no more | convention | a float tolerance (1e-6); earliest in this file: bab696f 2026-10-25 The membrane slab a wall for each loop, not a cost (plan 27 T5, P33-P35) | - | - |
| `membrane:LANE` | constant | 0.65 | LANE = 0.65 # page units beside an element's box a loop may climb in (rulings P34, P36a): route._INFLATE + 2 * RES | geometry | "route._INFLATE + 2 * RES"; earliest in this file: 2af59b6 2026-10-25 T5 fix round 1: a loop climbs the lane beside its own element, never through it … | membrane, route | - |
| `membrane:JOIN_GAP` | constant | 2 | JOIN_GAP = 2 # residues | named-entry:32KS,6UW4 | "…the membrane campaign's 21 views would have been: hairpins (32KS's NS1 and NS2, A:200-219 and A:222-233, their axes 151° apart; 6UW4's A:119-129 and…"; earliest in this file: 102a6f2 2026-10-24 Membrane view refinements from the ComEC look test (plan 26 Task 5, rulings T5-1… | sensitivity | - |
| `membrane:JOIN_ANGLE` | constant | 45.0 | JOIN_ANGLE = 45.0 # degrees between the two pieces' axes, N to C | unknown | earliest in this file: 102a6f2 2026-10-24 Membrane view refinements from the ComEC look test (plan 26 Task 5, rulings T5-1… | membrane | - |
| `membrane:_cylinder:comparison:0.9@483550b7` | comparison | 0.9 | e1 = np.cross(guess, [1.0, 0.0, 0.0] if abs(guess[0]) < 0.9 else [0.0, 1.0, 0.0]) | unknown | - | - | - |
| `membrane:_barrel_slab:comparison:10@fcb07d04` | comparison | 10 | if (out & inside).sum() < 10: | unknown | - | membrane | - |
| `membrane:_long:comparison:16@67ccdf8c` | comparison | 16 | return [s for s in sses if (s.kind == "H" and len(s) >= 16) or (s.kind == "E" and len(s) >= 8)] | unknown | - | - | - |
| `membrane:_long:comparison:8@67ccdf8c` | comparison | 8 | return [s for s in sses if (s.kind == "H" and len(s) >= 16) or (s.kind == "E" and len(s) >= 8)] | unknown | - | - | - |
| `membrane:_slab:comparison:20@9d669d15` | comparison | 20 | if inside.sum() < 20: | unknown | - | - | - |
| `membrane:_slab:argument:1.5@9b35601a` | argument | 1.5 | slabs.append((hyd[inside].sum() - 1.5 * pol[inside].sum(), c, inside)) | unknown | - | - | - |
| `membrane:page_depth:argument:1e-9@97115d88` | argument | 1e-9 | return up * (y - (b0 + b1) / 2) * 2 * half / max(b1 - b0, 1e-9) | convention | a floating-point tolerance | - | - |
| `membrane:band_eps:argument:1e-9@ea514d88` | argument | 1e-9 | return BAND_EPS * 2 * half / max(band[1] - band[0], 1e-9) | convention | a floating-point tolerance | - | - |
| `membrane:meets_box:comparison:1e-12@6503eb66` | comparison | 1e-12 | if abs(pk) < 1e-12: | convention | a floating-point tolerance | - | - |
| `membrane:loop_crossing:comparison:1e-6@2cfbacd1` | comparison | 1e-6 | cross = bound <= -m.half or (two and at[0] * at[1] < 0 and min(abs(at[0]), abs(at[1])) >= m.half - 1e-6) | convention | a floating-point tolerance | membrane, route | - |
| `membrane:_pieces:comparison:1e-12@4799bc14` | comparison | 1e-12 | if abs(dk) > 1e-12: | convention | a floating-point tolerance | - | - |
| `membrane:_orient:argument:5@4ef08ac3` | argument | 5 | up_all = int(band(m.half - 5, m.half + 15, every).sum()) | unknown | - | membrane | - |
| `membrane:_orient:argument:15@4ef08ac3` | argument | 15 | up_all = int(band(m.half - 5, m.half + 15, every).sum()) | unknown | - | membrane | - |
| `membrane:_orient:argument:15@de48b723` | argument | 15 | down_all = int(band(-m.half - 15, -m.half + 5, every).sum()) | unknown | - | membrane | - |
| `membrane:_orient:argument:5@de48b723` | argument | 5 | down_all = int(band(-m.half - 15, -m.half + 5, every).sum()) | unknown | - | membrane | - |

## model

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `model:AUTO_NT` | constant | 300 | AUTO_NT = 300 # na_extent auto: larger nucleic acids are cropped to the stems the protein touches | unknown | earliest in this file: 22ace5f 2026-10-10 na_extent (auto crops large RNA to touched stems), na_style; focus keeps the sub… | na_structure | - |
| `model:_MERGE_GAP` | constant | 12 | _MERGE_GAP = 12 # crop by stems: kept stretches of one strand closer than this many nucleotides are joined | unknown | earliest in this file: 22ace5f 2026-10-10 na_extent (auto crops large RNA to touched stems), na_style; focus keeps the sub… | - | - |
| `model:homolog_chains.alike:comparison:20@b282cc0a` | comparison | 20 | return len(common) >= 20 and same >= 0.9 * len(common) | unknown | - | diff, domains, heldout_copies, membrane +6 | - |
| `model:homolog_chains.alike:comparison:0.9@b282cc0a` | comparison | 0.9 | return len(common) >= 20 and same >= 0.9 * len(common) | unknown | - | diff, domains, heldout_copies, membrane +6 | - |
| `model:generic_label:comparison:26@0840011b` | comparison | 26 | return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[n] if n < 26 else f"S{n + 1}" | unknown | - | - | - |
| `model:crop_nucleic:default:flank` | default | 6 | def crop_nucleic(na: Nucleic, flank: int = 6, *, by: str = "contacts") -> Nucleic: | unknown | - | na_structure, nucleic | - |
| `model:crop_nucleic:argument:3@f8579048` | argument | 3 | for _ in range(3): # partners of every kept nucleotide are kept too | unknown | - | na_structure, nucleic | - |

## na_curve

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `na_curve:CURVE_STEP` | constant | 0.1 | CURVE_STEP = 0.1 # page units between the curve's samples | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_FAR` | constant | 2.0 | _FAR = 2.0 # two places on the band this many band widths apart along it are different stretches (_far_gap) | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_GYRE_GAP` | constant | 0.7 | _GYRE_GAP = 0.7 # page units kept between two gyres' bands (their letters) beyond the band's width: enough that the | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_GYRE_GROW` | constant | 1.15 | _GYRE_GROW, _GYRE_TRIES = 1.15, 8 # gyres parted by too little: their spacing grows by this factor, this often | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_GYRE_TRIES` | constant | 8 | _GYRE_GROW, _GYRE_TRIES = 1.15, 8 # gyres parted by too little: their spacing grows by this factor, this often | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_GYRE_SPREAD` | constant | math.radians(40) | _GYRE_SPREAD = math.radians(40) # the ring's radius at an angle: its pieces' radii within about this angle | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_RAMP_MIN` | constant | 0.125 | _RAMP_MIN = 0.125 # turns: the shortest stretch over which the band goes from one gyre to the other | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_R_MIN` | constant | 3.0 | _R_MIN = 3.0 # rail offsets: the least radius the curve bends at (the inner rail's steps keep 2/3 of their length) | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | na_curve | - |
| `na_curve:_EASE_STEP` | constant | 0.5 | _EASE_STEP, _EASE_ROUNDS = 0.5, 800 # page units between the samples eased (_eased), and at most this many rounds | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_EASE_ROUNDS` | constant | 800 | _EASE_STEP, _EASE_ROUNDS = 0.5, 800 # page units between the samples eased (_eased), and at most this many rounds | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_NAME_AT` | constant | (0.5, 0.42, 0.58, 0.34, 0.66, 0.26, 0.74, 0.18, 0.82) | _NAME_AT = (0.5, 0.42, 0.58, 0.34, 0.66, 0.26, 0.74, 0.18, 0.82) # fractions of the duplex its name is tried at | unknown | earliest in this file: c2531bb 2026-10-22 A bent duplex as one curve (plan 25 Task 1): na_bend style key, curve through th… | - | - |
| `na_curve:_ROOM` | constant | 1.0 | _ROOM = 1.0 # the inner gyre holds this many times the footprints (MARGIN round each) of the items inside it | unknown | earliest in this file: 29ea07a 2026-10-22 The protein laid out round a curved duplex (plan 25 Task 2): pieces end to end o… | - | - |
| `na_curve:_RING_TRIES` | constant | 8 | _RING_TRIES = 8 # the ring is grown about its middle at most this many times until it does | unknown | earliest in this file: 29ea07a 2026-10-22 The protein laid out round a curved duplex (plan 25 Task 2): pieces end to end o… | - | - |
| `na_curve:_RUNG_W` | constant | 3.0 | _RUNG_W = 3.0 # a binding element's pull to its contact rungs at full contact, in layout._W_CONTACT | unknown | earliest in this file: 29ea07a 2026-10-22 The protein laid out round a curved duplex (plan 25 Task 2): pieces end to end o… | - | - |
| `na_curve:_PART_STEP` | constant | 0.25 | _PART_STEP, _PART_TRIES = 0.25, 400 # page units a curve moves at a time while parted from another, at most this often | unknown | earliest in this file: b0e1602 2026-10-22 Curved duplexes fix round 1 (plan 25): lone duplexes only, checks no looser than… | - | - |
| `na_curve:_PART_TRIES` | constant | 400 | _PART_STEP, _PART_TRIES = 0.25, 400 # page units a curve moves at a time while parted from another, at most this often | unknown | earliest in this file: b0e1602 2026-10-22 Curved duplexes fix round 1 (plan 25): lone duplexes only, checks no looser than… | - | - |
| `na_curve:NACurve.frames:argument:1e-12@d785952e` | argument | 1e-12 | t = t / np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-12)[:, None] | convention | a floating-point tolerance | na_curve | - |
| `na_curve:_tangents:argument:1e-12@81f3494e` | argument | 1e-12 | return _snap(t / np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-12)[:, None]) | convention | a floating-point tolerance | - | - |
| `na_curve:_catmull_rom:argument:1e-6@6ccd6d76` | argument | 1e-6 | t1 = max(math.hypot(*(p1 - p0)), 1e-6) ** 0.5 | convention | a floating-point tolerance | - | - |
| `na_curve:_catmull_rom:argument:1e-6@16c67625` | argument | 1e-6 | t2 = t1 + max(chord, 1e-6) ** 0.5 | convention | a floating-point tolerance | - | - |
| `na_curve:_catmull_rom:argument:1e-6@4b773a27` | argument | 1e-6 | t3 = t2 + max(math.hypot(*(p3 - p2)), 1e-6) ** 0.5 | convention | a floating-point tolerance | - | - |
| `na_curve:_catmull_rom:argument:4@d82562d4` | argument | 4 | n = max(4, math.ceil(chord / CURVE_STEP)) | unknown | - | - | - |
| `na_curve:_far_gap:argument:0.25@641c7888` | argument | 0.25 | at = np.searchsorted(s, np.arange(0.0, float(s[-1]), 0.25)) | unknown | - | - | - |
| `na_curve:_eased:comparison:5@6c9af2c0` | comparison | 5 | if len(p) < 5: | unknown | - | - | - |
| `na_curve:_eased:comparison:1.1@6513bb92` | comparison | 1.1 | tight = np.concatenate([[False], turn * r_min * 1.1 > reach, [False]]) # 10% to spare (drawn through again) | unknown | - | - | - |
| `na_curve:_name_place.beside:argument:1e4@17ff7647` | argument | 1e4 | off = math.ceil((out + margin + inner) * 1e4) / 1e4 | unknown | - | na_curve | - |
| `na_curve:_name_place:argument:4@2cc924ae` | argument | 4 | for _ in range(4): # measured, as figcheck does, at the curve's point nearest the name's ink | unknown | - | - | - |
| `na_curve:_name_place:comparison:1e-3@97774bb3` | comparison | 1e-3 | if abs(near - at) < 1e-3: | unknown | - | - | - |
| `na_curve:_name_place:comparison:1e-9@9b42c260` | comparison | 1e-9 | if best is None or score > best[0] + 1e-9: | convention | a floating-point tolerance | - | - |
| `na_curve:_ring:argument:1.05@9ef82ccf` | argument | 1.05 | f = float(_snap(f * max(1.05, math.sqrt(need / math.pi) / max(r_in, band)))) | unknown | - | na_curve | - |
| `na_curve:lay_out:argument:9@b6bea5c3` | argument | 9 | band = np.asarray(rows, float).reshape(-1, 9) | unknown | - | na_curve | - |
| `na_curve:lay_out:argument:9@d814be2f` | argument | 9 | before = walls.band if walls.band is not None else np.full((len(walls.boxes), 9), np.nan) | unknown | - | na_curve | - |
| `na_curve:lay_out:argument:4@c3d38b87` | argument | 4 | np.vstack([walls.boxes, np.asarray(boxes, float).reshape(-1, 4)]), | convention | an array's shape | na_curve | - |
| `na_curve:lay_out:argument:4@52dc8f3a` | argument | 4 | walls = Walls(np.asarray(boxes, float).reshape(-1, 4), hard, band) | convention | an array's shape | na_curve | - |
| `na_curve:lay_out:argument:4@3d11c03a` | argument | 4 | wt = _RUNG_W * _W_CONTACT * min(4 * len(res), _CONTACT_CAP) / _CONTACT_CAP | unknown | - | na_curve | - |
| `na_curve:_boxes_meet:argument:4@e5d5fc96` | argument | 4 | a = np.array([ch.rect for ch in chords(c, n, stems)], float).reshape(-1, 4) | convention | an array's shape | na_curve | - |
| `na_curve:_boxes_meet:argument:4@307a871b` | argument | 4 | b = np.array([ch.rect for ch in chords(e, m, stems)], float).reshape(-1, 4) | convention | an array's shape | na_curve | - |
| `na_curve:_boxes_meet:comparison:1e-6@3cdd2367` | comparison | 1e-6 | return bool(((w > 1e-6) & (h > 1e-6)).any()) | convention | a floating-point tolerance | na_curve | - |
| `na_curve:_boxes_meet:comparison:1e-6@3cdd2367#2` | comparison | 1e-6 | return bool(((w > 1e-6) & (h > 1e-6)).any()) | convention | a floating-point tolerance | na_curve | - |
| `na_curve:_part:comparison:1e-9@5b0a4082` | comparison | 1e-9 | if need > 1e-9 or _boxes_meet(c, i, e, j, stems): | convention | a floating-point tolerance | na_curve | - |
| `na_curve:_part:comparison:1e-9@e0f939ed` | comparison | 1e-9 | way = way / size if size > 1e-9 else np.array([1.0, 0.0]) | convention | a floating-point tolerance | na_curve | - |
| `na_curve:_beyond:argument:3@61cd7f2a` | argument | 3 | for t in np.arange(0.0, 3 * c.band, 0.25): | unknown | - | - | - |
| `na_curve:_beyond:argument:0.25@61cd7f2a` | argument | 0.25 | for t in np.arange(0.0, 3 * c.band, 0.25): | unknown | - | - | - |
| `na_curve:_fractions:argument:1e-12@17882bc2` | argument | 1e-12 | u = np.clip((rel[:, 0] * seg[:, 0] + rel[:, 1] * seg[:, 1]) / np.maximum(lens**2, 1e-12), 0.0, 1.0) | convention | a floating-point tolerance | - | - |
| `na_curve:chords:argument:0.05@5b48374a` | argument | 0.05 | m = max(1, math.ceil((hi - lo) / max(c.rail / 2, 0.05))) | unknown | - | na_curve | - |
| `na_curve:chords:argument:1e-12@ba39f391` | argument | 1e-12 | across = np.array([-seg[j, 1], seg[j, 0]]) / max(float(np.hypot(*seg[j])), 1e-12) | convention | a floating-point tolerance | na_curve | - |

## na_ends

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `na_ends:SHORT_GAP` | constant | 3 | SHORT_GAP = 3 # a cropped gap of fewer nucleotides is drawn as a bare break, without its "… N nt …" label | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_ends:_TAIL_ROUTED` | constant | 6 | _TAIL_ROUTED = 6 # a free end of more nucleotides that no straight drawing keeps clear is routed out to free space | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_ends:_TURNS` | constant | (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4) | _TURNS = (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4) # a free end's bends | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_layout | - |
| `na_ends:_END_SPOTS` | constant | ((1.0, 0.0), (1.0, 0.6), (2.0, 0.0), (2.0, 0.6), (1.0, -0.6), (3.0, 0.0), (0.0, 2.4), (3.0, 1.2)) | _END_SPOTS = ((1.0, 0.0), (1.0, 0.6), (2.0, 0.0), (2.0, 0.6), (1.0, -0.6), (3.0, 0.0), (0.0, 2.4), (3.0, 1.2)) | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_ends:_TRACE_NEAR` | constant | 2 | _TRACE_NEAR = 2 # nucleotides: an element touching a trace lies under it this far either side of what it holds | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_ends:_OFF_STEPS` | constant | tuple(1.0 + k / 4 for k in range(7)) | _OFF_STEPS = tuple(1.0 + k / 4 for k in range(7)) # beads (_TRACE_HALF) a trace nucleotide may be nudged off an | unknown | earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | - | - |
| `na_ends:_OFF_TURNS` | constant | tuple(np.radians(a) for a in range(0, 360, 5)) | _OFF_TURNS = tuple(np.radians(a) for a in range(0, 360, 5)) # in any direction, the first that clears | named-entry:13FZ | "…t holding it, nearest first: up to 2.5 (random-100 held-out 13FZ: B:3 lay 1.07 inside A:1027-1031's end) in any direction, the first that clears"; earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | sensitivity | - |
| `na_ends:_GAP_BOWS` | constant | (1.2, 1.6, 2.2) | _GAP_BOWS = (1.2, 1.6, 2.2) # a gap with no way on the grid may bow out to this many times its chord | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_ends:_enters:argument:1e-12@2180f5c2` | argument | 1e-12 | t = np.clip(((centre - a) * ab).sum(axis=1) / np.maximum((ab * ab).sum(axis=1), 1e-12), 0.0, 1.0) | convention | a floating-point tolerance | na_layout | - |
| `na_ends:_tail:argument:1e-9@84bc8c5e` | argument | 1e-9 | return path, tip, v / max(float(np.hypot(*v)), 1e-9) | convention | a floating-point tolerance | na_layout | - |
| `na_ends:_toward:comparison:1e-9@38b8f27c` | comparison | 1e-9 | ahead = sorted((t for t in _TURNS if abs(t) <= np.pi / 2 + 1e-9), key=lambda t: _snap_gap(gap(t))) | convention | a floating-point tolerance | - | - |
| `na_ends:_toward:comparison:1e-9@3ac952a4` | comparison | 1e-9 | return (*ahead, *(t for t in _TURNS if abs(t) > np.pi / 2 + 1e-9)) # folding back along the stem only after | convention | a floating-point tolerance | - | - |
| `na_ends:_snap_gap:argument:6@82cb558f` | argument | 6 | return round(v, 6) # ties (and their order, _TURNS') kept the same on every machine | convention | rounding digits | - | - |
| `na_ends:_reach_out:default:max_nodes` | default | 50_000 | def _reach_out(grid: _Grid, start, heading: int, steps: int, max_nodes: int = 50_000): | unknown | - | na_layout | - |
| `na_ends:_reach_out:argument:4@354e2d77` | argument | 4 | for nd in range(4): | unknown | - | na_layout | - |
| `na_ends:_reach_out:comparison:4@dbeae7e2` | comparison | 4 | if nd == (d + 2) % 4: | unknown | - | na_layout | - |
| `na_ends:_routed_end:argument:4@985edad3` | argument | 4 | heading = min(range(4), key=lambda k: -_dot(_DIRS[k], out)) | unknown | - | na_layout | - |
| `na_ends:_routed_end:argument:4@526b62ee` | argument | 4 | cells = _reach_out(grid, start, heading, int(math.ceil(sum(gaps) / grid.res)) + 4) | unknown | - | na_layout | - |
| `na_ends:_routed_end:comparison:0.9@e59df786` | comparison | 0.9 | if min(math.dist(p, q) for k, p in enumerate(marks) for q in marks[k + 1 :]) < 0.9 * NT_STEP: | unknown | - | na_layout | - |
| `na_ends:_label_room:argument:0.05@e6b2a12c` | argument | 0.05 | return not any(_near(xy, e, hh + 0.05) for e in parts) | unknown | - | na_metrics | - |
| `na_ends:_side_start:argument:6@4f83f474` | argument | 6 | found.append((round(d, 6), c)) | convention | rounding digits | - | - |
| `na_ends:_tail:fallback-marker` | fallback-marker | 2 | drawn so far (`strands`); else it is drawn straight, a fallback. A drawing that keeps out of its partner's | unknown | a free end drawn straight: counted as na_tail_straight (route_na.route_nucleic, over the paths as drawn) | na_layout | - |
| `na_ends:off_unheld:fallback-marker` | fallback-marker | 1 | fallback). Run once the elements are placed (layout.layout): random-100 held-out 13FZ's RNA B, held by its sheet at | named-entry:13FZ | a trace nucleotide no nudge clears stays: counted as na_trace_unclear | na_text_keepoff | - |
| `na_ends:_gap:fallback-marker` | fallback-marker | 4 | fallback = route is None | unknown | counted as na_gap_unclear | na_layout, na_metrics | - |

## na_labels

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `na_labels:_TIE_SHORT` | constant | 3.0 | _TIE_SHORT = 3.0 # stem widths: a longer tie is drawn as two labelled marks, not a line | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_TICK_GAP` | constant | 0.15 | _TICK_GAP, _TICK = 0.15, 0.35 # a tie mark: a tick this far out from its nucleotide, this long | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_TICK` | constant | 0.35 | _TICK_GAP, _TICK = 0.15, 0.35 # a tie mark: a tick this far out from its nucleotide, this long | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_TIE_TURNS` | constant | (0.0, 45.0, -45.0, 90.0, -90.0, 135.0, -135.0, 180.0) | _TIE_TURNS = (0.0, 45.0, -45.0, 90.0, -90.0, 135.0, -135.0, 180.0) # degrees from the partner a mark may point | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_TICK_REACH` | constant | (1.0, 2.0, 3.0) | _TICK_REACH = (1.0, 2.0, 3.0) # times _TICK a mark's tick is drawn, the shortest that leaves its label clear | unknown | earliest in this file: 116618c 2026-10-17 A tie mark whose label finds no clear turn is drawn with a longer tick before ta… | - | - |
| `na_labels:_SLIDES` | constant | tuple(f * k for k in (0.1, 0.2, 0.3, 0.4, 0.48) for f in (1, -1)) | _SLIDES = tuple(f * k for k in (0.1, 0.2, 0.3, 0.4, 0.48) for f in (1, -1)) # spacings a run's nucleotides may slide | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_metrics | - |
| `na_labels:_EACH_STEP` | constant | 0.05 | _EACH_STEP = 0.05 # page units: the steps a nucleotide still on another's spot slides alone along its run (_slide_each) | unknown | earliest in this file: 460d6b6 2026-10-18 Held-out NA metrics (plan 15 task 4): beads, crop, gaps, counts on ties, end lab… | - | - |
| `na_labels:_LETTER_GAP` | constant | 0.38 | _LETTER_GAP = 0.38 # page units: closest two base letters (or a letter and a nucleotide) may sit | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_layout, na_text_keepoff | - |
| `na_labels:_LETTER_RINGS` | constant | tuple(NT_ROW * f for f in (3.8, 4.6, 5.5, 6.5)) | _LETTER_RINGS = tuple(NT_ROW * f for f in (3.8, 4.6, 5.5, 6.5)) # how far out a crowded letter looks last, nearest | unknown | earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | na_text_keepoff | - |
| `na_labels:_RING_TURNS` | constant | tuple(np.radians(a) for a in range(0, 360, 15)) | _RING_TURNS = tuple(np.radians(a) for a in range(0, 360, 15)) # first, every 15° round its nucleotide | unknown | earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | - | - |
| `na_labels:_GAP_ALONG` | constant | tuple(0.5 + s * k / 20 for k in range(9) for s in ((1,) if k == 0 else (-1, 1))) | _GAP_ALONG = tuple(0.5 + s * k / 20 for k in range(9) for s in ((1,) if k == 0 else (-1, 1))) # middle out | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_text_keepoff | - |
| `na_labels:_GAP_OUT` | constant | (1.0, 2.0, 3.0, 4.0, 6.0) | _GAP_OUT = (1.0, 2.0, 3.0, 4.0, 6.0) # how far beside its gap a label may sit, in multiples of the least that clears | named-entry:7Y7C | "…les of the least that clears it (4 and 6, random-100 C2b-6: 7Y7C's "... 62 nt ..." beside a gap 0.57 across, between stems)"; earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_layout, na_text_keepoff, sensitivity | - |
| `na_labels:_GAP_RINGS` | constant | tuple(1.0 + k / 4 for k in range(9)) | _GAP_RINGS = tuple(1.0 + k / 4 for k in range(9)) # then round its midpoint, in multiples of its half width and half | unknown | earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | na_text_keepoff | - |
| `na_labels:_GAP_TURNS` | constant | tuple(math.radians(a) for a in range(0, 360, 5)) | _GAP_TURNS = tuple(math.radians(a) for a in range(0, 360, 5)) # height, every 5°, nearest first (random-100 held-out | named-entry:8ZTU | "height, every 5°, nearest first (random-100 held-out 8ZTU: "... 5 nt ..." by a gap 0.54 long among runs, where every spot beside it met one; the one …"; earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | sensitivity | - |
| `na_labels:_GAP_TEXT` | constant | 0.75 | _GAP_TEXT = 0.75 # a gap label's text size, in ems of the figure's text (render draws it at 0.75 of the font) | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_GAP_CHAR` | constant | 0.5 | _GAP_CHAR = 0.5 # its mean character width, in ems of its own text ("…", digits and spaces are narrow) | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_GAP_ROOM` | constant | 1.15 | _GAP_ROOM = 1.15 # room for the page to grow (its bounds widen after routing) and for a larger font | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_labels:_outline:argument:0.45@1a594b62` | argument | 0.45 | head = min(HEAD_LEN, e.length * 0.45) | unknown | - | na_text_keepoff | - |
| `na_labels:_boxes_of:argument:4@68aeb51c` | argument | 4 | return np.empty((0, 4)) | convention | an array's shape | na_metrics, na_text_keepoff | - |
| `na_labels:_place_letters:argument:4@6e8955a0` | argument | 4 | strand_segs = np.vstack([np.hstack([q[:-1], q[1:]]) for q in lines]) if lines else np.empty((0, 4)) | convention | an array's shape | na_layout, na_text_keepoff | - |
| `na_labels:_place_letters.outline_room:argument:0.05@fd913a9d` | argument | 0.05 | box = ink.get(base, ink["N"]) + np.array([q[0] - 0.05, q[1] - 0.05, q[0] + 0.05, q[1] + 0.05]) | unknown | - | - | - |
| `na_labels:_place_letters.outline_room:argument:0.05@fd913a9d#2` | argument | 0.05 | box = ink.get(base, ink["N"]) + np.array([q[0] - 0.05, q[1] - 0.05, q[0] + 0.05, q[1] + 0.05]) | unknown | - | - | - |
| `na_labels:_place_letters.outline_room:argument:0.05@fd913a9d#3` | argument | 0.05 | box = ink.get(base, ink["N"]) + np.array([q[0] - 0.05, q[1] - 0.05, q[0] + 0.05, q[1] + 0.05]) | unknown | - | - | - |
| `na_labels:_place_letters.outline_room:argument:0.05@fd913a9d#4` | argument | 0.05 | box = ink.get(base, ink["N"]) + np.array([q[0] - 0.05, q[1] - 0.05, q[0] + 0.05, q[1] + 0.05]) | unknown | - | - | - |
| `na_labels:_place_letters:argument:1e-9@d757559a` | argument | 1e-9 | ahead = ahead / max(float(np.hypot(*ahead)), 1e-9) | convention | a floating-point tolerance | na_layout, na_text_keepoff | - |
| `na_labels:_label_gaps:argument:4@19ae71cd` | argument | 4 | seg_box = seg_box.reshape(-1, 4) | convention | an array's shape | na_layout, na_text_keepoff | - |
| `na_labels:_label_gaps.meets:comparison:1e-9@a5335b2e` | comparison | 1e-9 | hit = [k for k in near.tolist() if clip_len(*segs[k], box) > 1e-9] | convention | a floating-point tolerance | na_layout | - |
| `na_labels:_tie_name:comparison:26@d46be215` | comparison | 26 | return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n) | unknown | - | - | - |
| `na_labels:_tie_name:argument:26@d46be215` | argument | 26 | return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n) | unknown | - | - | - |
| `na_labels:_tie_name:argument:97@d46be215` | argument | 97 | return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n) | unknown | - | - | - |
| `na_labels:_tie_name:argument:26@d46be215#2` | argument | 26 | return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n) | unknown | - | - | - |
| `na_labels:_tie_name:argument:97@d46be215#2` | argument | 97 | return _tie_name(n // 26 - 1) + chr(97 + n % 26) if n >= 26 else chr(97 + n) | unknown | - | - | - |
| `na_labels:_ties:argument:4@a9e132e9` | argument | 4 | seg = np.asarray(segments, float).reshape(-1, 4) | convention | an array's shape | na_layout, na_text_keepoff | - |
| `na_labels:_ties.near_segments:argument:0.15@9ebe60a8` | argument | 0.15 | lo = np.minimum(box[:2], np.subtract(b, 0.15)) | unknown | - | - | - |
| `na_labels:_ties.near_segments:argument:0.15@f355cefd` | argument | 0.15 | hi = np.maximum(box[2:], np.add(b, 0.15)) | unknown | - | - | - |
| `na_labels:_ties.clear:comparison:1e-9@7dd31f5a` | comparison | 1e-9 | near = pts[np.hypot(*(pts - np.asarray(own, float)).T) > 1e-9] | convention | a floating-point tolerance | - | - |
| `na_labels:_ties.clear:comparison:0.2@244477d5` | comparison | 0.2 | if (_seg_dists(near, a, b) < 0.2).any(): | unknown | - | - | - |
| `na_labels:_ties.clear:comparison:1e-9@bb96f383` | comparison | 1e-9 | return not any(clip_len(p, q, box) > 1e-9 or _seg_dist(b, p, q) < 0.15 for p, q in near_segments(box, b)) | convention | a floating-point tolerance | - | - |
| `na_labels:_ties.clear:comparison:0.15@bb96f383` | comparison | 0.15 | return not any(clip_len(p, q, box) > 1e-9 or _seg_dist(b, p, q) < 0.15 for p, q in near_segments(box, b)) | unknown | - | - | - |
| `na_labels:_ties.misses:comparison:1e-9@618208a9` | comparison | 1e-9 | through = sum(1 for p, q in near_segments(box, b) if clip_len(p, q, box) > 1e-9) # lines across its label | convention | a floating-point tolerance | na_layout | - |
| `na_labels:_ties.misses:comparison:1e-9@7dd31f5a` | comparison | 1e-9 | near = pts[np.hypot(*(pts - np.asarray(own, float)).T) > 1e-9] | convention | a floating-point tolerance | na_layout | - |
| `na_labels:_ties:argument:1e-9@fdac1dac` | argument | 1e-9 | toward = toward / max(float(np.hypot(*toward)), 1e-9) | convention | a floating-point tolerance | na_layout, na_text_keepoff | - |
| `na_labels:_ties:argument:0.1@7ea7b703` | argument | 0.1 | options.append(((a, b), b + u * (hh + 0.1))) | unknown | - | na_layout, na_text_keepoff | - |

## na_paths

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `na_paths:_NA_CORNER` | constant | 0.6 | _NA_CORNER = 0.6 # radius of a routed strand's rounded corners | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_paths:TEXT_TOLL` | constant | 25.0 | TEXT_TOLL = 25.0 # the cost of a text cell a strand's route enters once no route keeps off text (_Router.route's | unknown | earliest in this file: 1e40581 2026-10-18 Held-out NA text against the protein (plan 15 task 3): letters, end, tie, gap an… | na_text_keepoff | - |
| `na_paths:_rounded:default:n` | default | 8 | def _rounded(pts, radius: float = _NA_CORNER, n: int = 8) -> list[tuple[float, float]]: | unknown | - | - | - |
| `na_paths:_rounded:comparison:1e-9@816a462d` | comparison | 1e-9 | if r < 1e-9: | convention | a floating-point tolerance | - | - |
| `na_paths:_bow:default:n` | default | 32 | def _bow(pa, pb, need: float, n: int = 32, out=None) -> list[tuple[float, float]]: | unknown | - | - | - |
| `na_paths:_bow:argument:1e-9@2dcf3c24` | argument | 1e-9 | chord = max(float(np.hypot(*d)), 1e-9) | convention | a floating-point tolerance | - | - |
| `na_paths:_bow:comparison:4@567ec757` | comparison | 4 | if _length(pts) >= need or h > 4 * need: | unknown | - | - | - |
| `na_paths:_crosses_segment:comparison:-1e-12@598e4285` | comparison | -1e-12 | if turn(p, q, a) * turn(p, q, b) < -1e-12 and turn(a, b, p) * turn(a, b, q) < -1e-12: | convention | a floating-point tolerance | na_layout | - |
| `na_paths:_crosses_segment:comparison:-1e-12@598e4285#2` | comparison | -1e-12 | if turn(p, q, a) * turn(p, q, b) < -1e-12 and turn(a, b, p) * turn(a, b, q) < -1e-12: | convention | a floating-point tolerance | na_layout | - |
| `na_paths:in_body:comparison:1e-9@61eda2fc` | comparison | 1e-9 | return abs(u) < d.length / 2 - 1e-9 and abs(v) < d.width / 2 - 1e-9 | convention | a floating-point tolerance | na_curve, na_layout, route | - |
| `na_paths:in_body:comparison:1e-9@61eda2fc#2` | comparison | 1e-9 | return abs(u) < d.length / 2 - 1e-9 and abs(v) < d.width / 2 - 1e-9 | convention | a floating-point tolerance | na_curve, na_layout, route | - |
| `na_paths:_Strands.near:argument:1e-12@225a1965` | argument | 1e-12 | t = ((q[None] - a) * ab).sum(axis=2) / np.maximum((ab * ab).sum(axis=2), 1e-12) | convention | a floating-point tolerance | - | - |
| `na_paths:_Strands.crossings:comparison:-1e-12@b421cb47` | comparison | -1e-12 | crossed = (turn(p0, p1, u) * turn(p0, p1, v) < -1e-12) & (turn(u, v, p0) * turn(u, v, p1) < -1e-12) | convention | a floating-point tolerance | layout, na_layout, route | - |
| `na_paths:_Strands.crossings:comparison:-1e-12@b421cb47#2` | comparison | -1e-12 | crossed = (turn(p0, p1, u) * turn(p0, p1, v) < -1e-12) & (turn(u, v, p0) * turn(u, v, p1) < -1e-12) | convention | a floating-point tolerance | layout, na_layout, route | - |
| `na_paths:graze:argument:33@ee3e54d9` | argument | 33 | t = np.linspace(0.0, 1.0, 33) | unknown | - | na_layout | - |
| `na_paths:crosses_body:comparison:1e-9@70fe6bf3` | comparison | 1e-9 | return clip_len((u0, v0), (u1, v1), (-e.length / 2, -e.width / 2, e.length / 2, e.width / 2)) > 1e-9 | convention | a floating-point tolerance | na_layout | - |
| `na_paths:<module>:fallback-marker` | fallback-marker | 2 | indices drawn on it, 5'->3', at nt_xy along `points`. fallback: no drawing clear of the elements was found. | unknown | NAPath.fallback, the flag the na_* counts read | - | - |

## na_runs

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `na_runs:_BOW_MAX` | constant | 8 | _BOW_MAX = 8 # internal-loop sides up to this many nucleotides are drawn as bows, not routed | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_runs:_JOINT_STRAIGHT` | constant | 1.5 | _JOINT_STRAIGHT = 1.5 # a joint between two pieces of a bent duplex shorter than this may be drawn straight | unknown | earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_runs:_JOINT_NEAR` | constant | 2.5 | _JOINT_NEAR = 2.5 # ... and, when no route reaches it, one shorter than this (9C75 P2b -> P2a, 1.83: random-100 C24), | named-entry:9C75 | "... and, when no route reaches it, one shorter than this (9C75 P2b -> P2a, 1.83: random-100 C24),"; earliest in this file: 0bb8b58 2026-10-17 A straight joint beyond 2.5 page units must leave its face within 60 degrees of … | na_layout | - |
| `na_runs:_JOINT_REACH` | constant | 5.0 | _JOINT_REACH = 5.0 | named-entry:8VWQ | "…r piece: the duplex reads as carrying on (random-100 C2b-4: 8VWQ's P26a -> P26b, 4.41 apart but 76° off its face, stays a curve)"; earliest in this file: 324f501 2026-10-14 C24 (random-100): a joint no route reaches is drawn straight up to 2.5 units, no… | na_layout | - |
| `na_runs:_JOINT_FORWARD` | constant | 0.5 | _JOINT_FORWARD = 0.5 | unknown | earliest in this file: 0bb8b58 2026-10-17 A straight joint beyond 2.5 page units must leave its face within 60 degrees of … | na_layout | - |
| `na_runs:_JOINT_BENDS` | constant | (1.0, 1.5, 2.0, 0.5, 3.0) | _JOINT_BENDS = (1.0, 1.5, 2.0, 0.5, 3.0) | named-entry:8VWQ | "…irst clear taken: the duplex bends round; random-100 C2b-6 (8VWQ's P26a <-> P26b, 73-76° off their faces; P124a's pieces turned back)"; earliest in this file: 1ed1b1b 2026-10-17 A joint between two pieces of a bent duplex that no route or straight line reach… | sensitivity | - |
| `na_runs:_JOINT_SWINGS` | constant | (0.0, 1.0, 2.0, -1.0, -2.0) | _JOINT_SWINGS = (0.0, 1.0, 2.0, -1.0, -2.0) | named-entry:8VWQ,7Y7C | "… its rail's side first (pieces turned back from each other, 8VWQ's P124a and 7Y7C's P103, meet past the corners of their bodies, which reach beyond …"; earliest in this file: 1ed1b1b 2026-10-17 A joint between two pieces of a bent duplex that no route or straight line reach… | na_layout, sensitivity | - |
| `na_runs:_BESIDE_TOLL` | constant | 50.0 | _BESIDE_TOLL = 50.0 | named-entry:7Y7C | "… stack gap's cells (route._GAP_TOLL); C2b-6 review minor 3: 7Y7C's a:520-521 ran 1.2 along a strand, through its nucleotide"; earliest in this file: 0521d64 2026-10-17 A route from the free cells beside a face, in the last pass, pays a high toll fo… | na_layout | - |
| `na_runs:_DETOUR` | constant | 3.0 | _DETOUR = 3.0 | measured | "…wing, which keeps off the hairpin loops and free ends only. Measured on 24 campaign entries, crossing a strand only a cost: run-run crossings 250 ->…" (entries 4YBB, 2GIS, 4P9R); earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_layout | - |
| `na_runs:_DETOUR_STEPS` | constant | 6.0 | _DETOUR_STEPS = 6.0 | measured | "… drawn 5.7 times its direct distance, 1GID's P15->P16 3.6). Measured on all 39 entries: 6 keeps every entry's protein-loop crossings at or below ste…" (entries 4YBB, 1GID, 5B43); earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_layout | - |
| `na_runs:_FIRST_NODES` | constant | 5_000 | _FIRST_NODES = 5_000 | measured | "5000 (was 20000) keeps 4YBB's routing within 1.5x of step 2's time and changes no entry's crossings but 1Y26's"; earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | - | - |
| `na_runs:_BESIDE_DETOUR` | constant | 2.0 | _BESIDE_DETOUR = 2.0 | measured | "… when at most this many times as long as its plain drawing. Measured on all 39 campaign entries with the first stage: run-run crossings 255 -> 237, …" (entries 4YBB, 5B43, 3CRX); earliest in this file: b43695e 2026-10-11 Split route_na.py (1480 lines) into na_paths, na_runs, na_ends and na_labels | na_layout | - |
| `na_runs:_arc:argument:24@316c694c` | argument | 24 | t = np.linspace(a0, a1, max(24, 4 * (m + 2))) | unknown | - | na_layout | - |
| `na_runs:_arc:argument:4@316c694c` | argument | 4 | t = np.linspace(a0, a1, max(24, 4 * (m + 2))) | unknown | - | na_layout | - |
| `na_runs:_elbow:comparison:1e-9@d5b4aeb4` | comparison | 1e-9 | if abs(det) > 1e-9: # a + ea t = b + eb u, both ahead | convention | a floating-point tolerance | na_layout | - |
| `na_runs:_bend:default:n` | default | 24 | def _bend(fa, ea, fb, eb, h: float, n: int = 24, sa=(0.0, 0.0), sb=(0.0, 0.0)) -> list[tuple[float, float]]: | unknown | line last changed by 1ed1b1b A joint between two pieces of a bent duplex that no route or straight line reaches is drawn as a smooth bend … (git blame; evidence only) | - | - |
| `na_runs:_straight_joint:comparison:-1e-9@db752c4c` | comparison | -1e-9 | touching = min(a[2], b[2]) - max(a[0], b[0]) >= -1e-9 and min(a[3], b[3]) - max(a[1], b[1]) >= -1e-9 | convention | a floating-point tolerance | na_layout | - |
| `na_runs:_straight_joint:comparison:-1e-9@db752c4c#2` | comparison | -1e-9 | touching = min(a[2], b[2]) - max(a[0], b[0]) >= -1e-9 and min(a[3], b[3]) - max(a[1], b[1]) >= -1e-9 | convention | a floating-point tolerance | na_layout | - |
| `na_runs:_takes_over:comparison:1e-9@c44d0f8c` | comparison | 1e-9 | _length(path.points) <= longest + 1e-9 | convention | a floating-point tolerance | na_layout | - |
| `na_runs:_off_names:comparison:1e-9@8f01be24` | comparison | 1e-9 | if b[0] < hi_x and lo_x < b[2] and b[1] < hi_y and lo_y < b[3] and clip_len(p, q, b) > 1e-9: | convention | a floating-point tolerance | na_layout | - |
| `na_runs:_strand_between:fallback-marker` | fallback-marker | 1 | "hairpin", s, nts, pts, _spaced(pts, m), fallback=not _na_clear(pts, own, router.elements), between=between | unknown | a hairpin arc not clear: counted as na_hairpin_unclear | na_layout | - |
| `na_runs:_takes_over:fallback-marker` | fallback-marker | 2 | """May a later stage's drawing of a strand with m nucleotides replace its plain drawing? Only if it is no fallback | unknown | not a fallback: weighs two drawings' fallback flags | na_layout | - |
| `na_runs:_drawing:fallback-marker` | fallback-marker | 7 | fallback = route is None | unknown | counted as na_run_unclear | inventory, na_layout | - |

## nacheck

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `nacheck:THRESHOLDS` | constant | { "pair_f1": 0.95, # canonical + wobble pairs against RNApolis (final campaign: 0.96-1.00; started at 0.90) "stem_recal… | THRESHOLDS = { | unknown | earliest in this file: 473f099 2026-10-11 nacheck: shape, contact and rendering fidelity of nucleic-acid figures; fixture … | na_layout, na_validation | - |
| `nacheck:SIDE_SEEN` | constant | 2.75 | SIDE_SEEN = 2.75 | geometry | "Equal to the layout's own cut, layout._SIDE_CLEAR"; earliest in this file: 473f099 2026-10-11 nacheck: shape, contact and rendering fidelity of nucleic-acid figures; fixture … | na_validation | - |
| `nacheck:_heavy:argument:3@773a76ee` | argument | 3 | return np.array(xyz).reshape(-1, 3), np.array(owner, int), names | convention | an array's shape | - | - |
| `nacheck:_similarity_rmsd:argument:1e-12@83c6c6d2` | argument | 1e-12 | residual = max(aa - (dot * dot + cross * cross) / max(bb, 1e-12), 0.0) | convention | a floating-point tolerance | na_validation | - |
| `nacheck:_similarity_rmsd:argument:1e-9@2d68066b` | argument | 1e-9 | return float(np.sqrt(residual / len(a)) / max(np.sqrt(aa / len(a)), 1e-9)) | convention | a floating-point tolerance | na_validation | - |
| `nacheck:_angle:argument:1e-12@11ff906f` | argument | 1e-12 | cos = float((u * v).sum()) / max(float(np.sqrt((u * u).sum() * (v * v).sum())), 1e-12) | convention | a floating-point tolerance | - | - |
| `nacheck:shape:comparison:5@d3512330` | comparison | 5 | if len(stems) >= 5: # fewer stems give too few distances for a rank correlation | unknown | - | hbonds, na_curve, na_layout, na_validation +1 | - |
| `nacheck:shape:argument:3@3cdc9535` | argument | 3 | out["spearman"] = round(float(spearmanr(d3, d2).statistic), 3) | convention | rounding digits | hbonds, na_curve, na_layout, na_validation +1 | - |
| `nacheck:shape:comparison:3@9704058b` | comparison | 3 | if len(stems) >= 3: | unknown | - | hbonds, na_curve, na_layout, na_validation +1 | - |
| `nacheck:shape:argument:3@ceff3d3c` | argument | 3 | out["procrustes"] = round(_similarity_rmsd(frame.project(c3) * SCALE, drawn), 3) | convention | rounding digits | hbonds, na_curve, na_layout, na_validation +1 | - |
| `nacheck:contacts:argument:3@92bf1457` | argument | 3 | "contact_auc": round(float(np.mean(aucs)), 3) if aucs else None, | convention | rounding digits | layout, na_contacts, na_layout, na_metrics +1 | - |
| `nacheck:contacts:argument:3@7b296ef8` | argument | 3 | "contact_precision_at_k": round(float(np.mean(precs)), 3) if precs else None, | convention | rounding digits | layout, na_contacts, na_layout, na_metrics +1 | - |
| `nacheck:contacts:argument:3@d75cbf4b` | argument | 3 | "stems_nearer_frac": round(sum(nearer) / len(nearer), 3) if nearer else None, | convention | rounding digits | layout, na_contacts, na_layout, na_metrics +1 | - |
| `nacheck:_ss_contacts_nearer:argument:3@48960175` | argument | 3 | return round(sum(nearer) / len(nearer), 3) | convention | rounding digits | - | - |
| `nacheck:_box_hits:comparison:1e-9@9895dca8` | comparison | 1e-9 | return any(clip_len(s[:2], s[2:], box) > 1e-9 for s in seg) | convention | a floating-point tolerance | - | - |
| `nacheck:_whole:argument:1e-9@1125114d` | argument | 1e-9 | return ("curve", n, int(((d.s0 + d.s1) / 2) // max(_FAR * d.band, 1e-9))) | convention | a floating-point tolerance | - | - |
| `nacheck:render_checks:argument:4@4ced30c9` | argument | 4 | ink = np.array([e.extents for e in letters], float).reshape(-1, 4) | convention | an array's shape | na_layout, na_validation, perf | - |
| `nacheck:render_checks:comparison:3@e64f5bc1` | comparison | 3 | if w > 0 and h > 0 and w * h > min(a.width * a.height, b.width * b.height) / 3: | unknown | - | na_layout, na_validation, perf | - |
| `nacheck:render_checks:comparison:1e-6@1561dea6` | comparison | 1e-6 | bead_overlaps += int(np.count_nonzero(np.hypot(d[:, 0], d[:, 1]) < size[i] + size[i + 1 :] - 1e-6)) | convention | a floating-point tolerance | na_layout, na_validation, perf | - |
| `nacheck:<module>:fallback-marker` | fallback-marker | 1 | "fallbacks": 0, # strand or loop routes that fell back to a curve | unknown | not a fallback: the check counting the loop and nucleic-acid fallbacks (nacheck fallbacks) | - | - |
| `nacheck:render_checks:fallback-marker` | fallback-marker | 1 | "fallbacks": sum(p.fallback for p in lay.na_paths) + sum(l.fallback for l in loops), | unknown | not a fallback: the check counting the loop and nucleic-acid fallbacks (nacheck fallbacks) | na_layout, na_validation, perf | - |

## names

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `names:_PIECE_GAP` | constant | 2 | _PIECE_GAP = 2 # residues: helices of one chain this close are pieces of one helix (ss: a kink, a copy's break) | geometry | "(ss: a kink, a copy's break)": the two-residue gap of features.KINK_GAP and ss._BREAK_MAX; earliest in this file: 9a6ff77 2026-10-15 7OKA (random-100, fix round 1): the copy-majority helix rule is told (README, sc… | - | - |
| `names:_helix_rows:argument:4@decc1a40` | argument | 4 | if row.has(4) and not row[4].startswith("HELX"): | unknown | - | - | - |

## nucleic

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `nucleic:HBOND_MAX` | constant | 3.8 | HBOND_MAX = 3.8 # Å, donor-acceptor distance of a base-pair hydrogen bond (stretched pairs at 2-3 Å resolution) | biology | "donor-acceptor distance of a base-pair hydrogen bond (stretched pairs at 2-3 Å resolution)"; earliest in this file: 9dc5dfc 2026-10-10 Base pairs from hydrogen bonds and base geometry, checked against RNApolis | independent_kernels | T7 finding F2: the base-pair kernel holds the bonds to 3.5 Å (tests/test_independent_kernels.py, xfail strict) |
| `nucleic:HBOND_STRETCH` | constant | 4.0 | HBOND_STRETCH = 4.0 # Å, the one longer bond a pair stacked in register on a strict pair may have (1HYS, 1AOI) | named-entry:1HYS,1AOI | "… bond a pair stacked in register on a strict pair may have (1HYS, 1AOI)"; earliest in this file: 5b20853 2026-10-11 C11: a pair stacked in register on a strict pair is accepted with its second hyd… | independent_kernels | T7 finding F2: the base-pair kernel holds the bonds to 3.5 Å (tests/test_independent_kernels.py, xfail strict) |
| `nucleic:C1_PAIR` | constant | (9.0, 12.0) | C1_PAIR = (9.0, 12.0) # Å, C1'-C1' across a Watson-Crick or wobble pair (10.5 +- 1.5) | biology | "C1'-C1' across a Watson-Crick or wobble pair (10.5 +- 1.5)"; earliest in this file: 9dc5dfc 2026-10-10 Base pairs from hydrogen bonds and base geometry, checked against RNApolis | independent_kernels | T7 finding F3: the kernel takes 9.5-11.5 Å (tests/test_independent_kernels.py, xfail strict) |
| `nucleic:PLANE_MAX` | constant | 45.0 | PLANE_MAX = 45.0 # degrees between the two base planes (TBP and nucleosome DNA reach 40 of propeller twist) | biology | "TBP and nucleosome DNA reach 40 of propeller twist"; earliest in this file: 9dc5dfc 2026-10-10 Base pairs from hydrogen bonds and base geometry, checked against RNApolis | independent_kernels, na_pairs | T7 finding F4: the kernel has no planarity bound (tests/test_independent_kernels.py, xfail strict) |
| `nucleic:MIN_LOOP` | constant | 3 | MIN_LOOP = 3 # a hairpin closes at least this many unpaired nucleotides | biology | a hairpin loop closes at least three unpaired nucleotides; earliest in this file: 9dc5dfc 2026-10-10 Base pairs from hydrogen bonds and base geometry, checked against RNApolis | na_validation | - |
| `nucleic:SEGMENT_TURN` | constant | 25.0 | SEGMENT_TURN = 25.0 # degrees the duplex axis may turn before a new straight piece starts | unknown | earliest in this file: 8f3f2b7 2026-10-10 Stem axes; bent duplexes cut into straight pieces | - | - |
| `nucleic:PIECE_MIN` | constant | 5 | PIECE_MIN = 5 # base pairs in the shortest straight piece | unknown | earliest in this file: 8f3f2b7 2026-10-10 Stem axes; bent duplexes cut into straight pieces | - | - |
| `nucleic:SEGMENT_MIN` | constant | 2 * PIECE_MIN | SEGMENT_MIN = 2 * PIECE_MIN # base pairs in the shortest stem that can split (into two shortest pieces) | geometry | 2 * PIECE_MIN; earliest in this file: 8f3f2b7 2026-10-10 Stem axes; bent duplexes cut into straight pieces | - | - |
| `nucleic:AXIS_WINDOW` | constant | 11 | AXIS_WINDOW = 11 # base pairs (about one helical turn) whose C1' atoms give the local duplex axis | biology | "about one helical turn"; earliest in this file: 8f3f2b7 2026-10-10 Stem axes; bent duplexes cut into straight pieces | - | - |
| `nucleic:_PAIRS` | constant | { ("G", "C"): ((("N1", "N3"), ("O6", "N4"), ("N2", "O2")), 2, "WC"), ("I", "C"): ((("N1", "N3"), ("O6", "N4")), 2, "WC"… | _PAIRS = { | biology | base-pair hydrogen-bond tables (chemistry); earliest in this file: 9dc5dfc 2026-10-10 Base pairs from hydrogen bonds and base geometry, checked against RNApolis | - | - |
| `nucleic:TETRAD_HBOND` | constant | 3.5 | TETRAD_HBOND = 3.5 # Å, N1···O6 and N2···N7 of a Hoogsteen G·G pair in a G-tetrad (2.8-3.1 in solved quadruplexes) | biology | "Hoogsteen G·G pair in a G-tetrad (2.8-3.1 in solved quadruplexes)"; earliest in this file: 8bbf623 2026-10-14 C29 (random-100): G-tetrads are read and the note says they are not drawn (8Q4O,… | - | - |
| `nucleic:_base_atom:comparison:1.7@ebfe9ef3` | comparison | 1.7 | if c1 is not None and a is not None and c1.pos.dist(a.pos) < 1.7: | unknown | - | - | - |
| `nucleic:_cis:comparison:-0.5@7a93c8c1` | comparison | -0.5 | return float(n @ m) < -0.5 * float(np.linalg.norm(n) * np.linalg.norm(m)) | unknown | - | na_pairs | - |
| `nucleic:_normal:comparison:4@a907f548` | comparison | 4 | if len(pts) < 4: | unknown | - | - | - |
| `nucleic:g_tetrads:comparison:4@1230e2ae` | comparison | 4 | if a in nxt.get(d, ()) and len({a, b, c, d}) == 4: | unknown | line last changed by 8bbf623 C29 (random-100): G-tetrads are read and the note says they are not drawn (8Q4O, 8BW5) (git blame; evidence only) | - | - |
| `nucleic:_local_axes:argument:3@9f927513` | argument | 3 | out = np.zeros((n, 3)) | convention | an array's shape | - | - |
| `nucleic:_letter:argument:26@adc1d80c` | argument | 26 | n, r = divmod(n - 1, 26) | unknown | - | na_structure | - |
| `nucleic:straight_duplex:comparison:0.8@bd2ec60d` | comparison | 0.8 | or len(biggest) < 0.8 * total | unknown | - | na_structure | - |

## palette

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `palette:BACK_SHADE` | constant | 0.72 | BACK_SHADE = 0.72 # a helix's back face is its colour darkened this much (render) | unknown | earliest in this file: 0f7fac0 2026-10-15 C35 (random-100, batch C1 review): chains past the palette take a distinct light… | render | - |
| `palette:CHAIN_DE` | constant | 10.0 | CHAIN_DE = 10.0 # CIEDE2000: a chain's colour past the palette's first pass stands this far from every other | unknown | earliest in this file: 0f7fac0 2026-10-15 C35 (random-100, batch C1 review): chains past the palette take a distinct light… | render | - |
| `palette:_PASS_STEPS` | constant | (("tint", 0.45), ("darken", 0.75), ("tint", 0.7), ("darken", 0.55), ("tint", 0.25), ("darken", 0.85)) | _PASS_STEPS = (("tint", 0.45), ("darken", 0.75), ("tint", 0.7), ("darken", 0.55), ("tint", 0.25), ("darken", 0.85)) | unknown | earliest in this file: 0f7fac0 2026-10-15 C35 (random-100, batch C1 review): chains past the palette take a distinct light… | - | - |
| `palette:_LAB_STEPS` | constant | tuple(range(-60, 61, 3)) | _LAB_STEPS = tuple(range(-60, 61, 3)) # L* steps tried when none of _PASS_STEPS stands CHAIN_DE apart | unknown | earliest in this file: 0f7fac0 2026-10-15 C35 (random-100, batch C1 review): chains past the palette take a distinct light… | - | - |
| `palette:_more_colours:comparison:8@e4169485` | comparison | 8 | ladder = [_from_lab(at[0] + d, at[1], at[2]) for d in _LAB_STEPS if 8 <= at[0] + d <= 92] | unknown | - | - | - |
| `palette:_more_colours:comparison:92@e4169485` | comparison | 92 | ladder = [_from_lab(at[0] + d, at[1], at[2]) for d in _LAB_STEPS if 8 <= at[0] + d <= 92] | unknown | - | - | - |
| `palette:_rgb:argument:16@9059d0b0` | argument | 16 | return tuple(int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)) # type: ignore[return-value] | unknown | - | - | - |
| `palette:_rgb:argument:255@9059d0b0` | argument | 255 | return tuple(int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)) # type: ignore[return-value] | unknown | - | - | - |
| `palette:_rgb:argument:3@9059d0b0` | argument | 3 | return tuple(int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)) # type: ignore[return-value] | unknown | - | - | - |
| `palette:_rgb:argument:5@9059d0b0` | argument | 5 | return tuple(int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)) # type: ignore[return-value] | unknown | - | - | - |
| `palette:darken:default:factor` | default | 0.6 | def darken(hex_color: str, factor: float = 0.6) -> str: | unknown | - | render | - |
| `palette:tint:default:amount` | default | 0.6 | def tint(hex_color: str, amount: float = 0.6) -> str: | unknown | - | domains | - |
| `palette:text_color_on:comparison:0.2126@587f2001` | comparison | 0.2126 | return "black" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.6 else "white" | unknown | - | - | - |
| `palette:text_color_on:comparison:0.7152@587f2001` | comparison | 0.7152 | return "black" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.6 else "white" | unknown | - | - | - |
| `palette:text_color_on:comparison:0.0722@587f2001` | comparison | 0.0722 | return "black" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.6 else "white" | unknown | - | - | - |
| `palette:text_color_on:comparison:0.6@587f2001` | comparison | 0.6 | return "black" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.6 else "white" | unknown | - | - | - |
| `palette:_lab:comparison:0.04045@98aefd01` | comparison | 0.04045 | lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(colour)] | unknown | - | - | - |
| `palette:_lab.f:comparison:216@0c97366a` | comparison | 216 | return t ** (1 / 3) if t > 216 / 24389 else (24389 / 27 * t + 16) / 116 | unknown | - | - | - |
| `palette:_lab.f:comparison:24389@0c97366a` | comparison | 24389 | return t ** (1 / 3) if t > 216 / 24389 else (24389 / 27 * t + 16) / 116 | unknown | - | - | - |
| `palette:_from_lab.finv:comparison:3@777311c0` | comparison | 3 | return t**3 if t**3 > 216 / 24389 else (116 * t - 16) / (24389 / 27) | unknown | - | - | - |
| `palette:_from_lab.finv:comparison:216@777311c0` | comparison | 216 | return t**3 if t**3 > 216 / 24389 else (116 * t - 16) / (24389 / 27) | unknown | - | - | - |
| `palette:_from_lab.finv:comparison:24389@777311c0` | comparison | 24389 | return t**3 if t**3 > 216 / 24389 else (116 * t - 16) / (24389 / 27) | unknown | - | - | - |
| `palette:_from_lab.gamma:comparison:0.0031308@1f23bbdc` | comparison | 0.0031308 | return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055 | unknown | - | - | - |
| `palette:delta_e_lab:argument:7@14bd4860` | argument | 7 | g = 0.5 * (1 - math.sqrt(((c1 + c2) / 2) ** 7 / (((c1 + c2) / 2) ** 7 + 25**7))) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:7@14bd4860#2` | argument | 7 | g = 0.5 * (1 - math.sqrt(((c1 + c2) / 2) ** 7 / (((c1 + c2) / 2) ** 7 + 25**7))) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:25@14bd4860` | argument | 25 | g = 0.5 * (1 - math.sqrt(((c1 + c2) / 2) ** 7 / (((c1 + c2) / 2) ** 7 + 25**7))) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:7@14bd4860#3` | argument | 7 | g = 0.5 * (1 - math.sqrt(((c1 + c2) / 2) ** 7 / (((c1 + c2) / 2) ** 7 + 25**7))) | unknown | - | render, style | - |
| `palette:delta_e_lab:comparison:180@2f66c817` | comparison | 180 | elif abs(h2p - h1p) <= 180: | unknown | - | render, style | - |
| `palette:delta_e_lab:comparison:180@499448a8` | comparison | 180 | elif abs(h1p - h2p) <= 180: | unknown | - | render, style | - |
| `palette:delta_e_lab:comparison:360@de23a158` | comparison | 360 | hm = (h1p + h2p + (360 if h1p + h2p < 360 else -360)) / 2 | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:30@ddd93084` | argument | 30 | - 0.17 * math.cos(math.radians(hm - 30)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:3@eecd9779` | argument | 3 | + 0.32 * math.cos(math.radians(3 * hm + 6)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:6@eecd9779` | argument | 6 | + 0.32 * math.cos(math.radians(3 * hm + 6)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:4@971406e0` | argument | 4 | - 0.20 * math.cos(math.radians(4 * hm - 63)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:63@971406e0` | argument | 63 | - 0.20 * math.cos(math.radians(4 * hm - 63)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:20@cd06e7d2` | argument | 20 | sl = 1 + 0.015 * (lm - 50) ** 2 / math.sqrt(20 + (lm - 50) ** 2) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:50@cd06e7d2` | argument | 50 | sl = 1 + 0.015 * (lm - 50) ** 2 / math.sqrt(20 + (lm - 50) ** 2) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:7@6bdca369` | argument | 7 | rc = 2 * math.sqrt(cm**7 / (cm**7 + 25**7)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:7@6bdca369#2` | argument | 7 | rc = 2 * math.sqrt(cm**7 / (cm**7 + 25**7)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:25@6bdca369` | argument | 25 | rc = 2 * math.sqrt(cm**7 / (cm**7 + 25**7)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:7@6bdca369#3` | argument | 7 | rc = 2 * math.sqrt(cm**7 / (cm**7 + 25**7)) | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:60@9777187f` | argument | 60 | rt = -math.sin(math.radians(60 * math.exp(-(((hm - 275) / 25) ** 2)))) * rc | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:275@9777187f` | argument | 275 | rt = -math.sin(math.radians(60 * math.exp(-(((hm - 275) / 25) ** 2)))) * rc | unknown | - | render, style | - |
| `palette:delta_e_lab:argument:25@9777187f` | argument | 25 | rt = -math.sin(math.radians(60 * math.exp(-(((hm - 275) / 25) ** 2)))) * rc | unknown | - | render, style | - |

## panels

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `panels:GRID` | constant | 0.1 | GRID = 0.1 # page units per raster cell | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | domains, perf | - |
| `panels:CLOSE` | constant | 2.0 | CLOSE = 2.0 # gaps between a domain's padded elements up to twice this are filled (2LZM: 1.2 left bays) | named-entry:2LZM | "…een a domain's padded elements up to twice this are filled (2LZM: 1.2 left bays)"; earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | perf | - |
| `panels:GUTTER` | constant | 0.3 | GUTTER = 0.3 # clear space kept between two domains' panels, where their elements leave room for it | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | perf | - |
| `panels:NAME_PAD` | constant | 0.15 | NAME_PAD = 0.15 # a name's box is wrapped this far | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | perf | - |
| `panels:EDGE` | constant | 0.05 | EDGE = 0.05 # each element's own box, grown this far (never past half way to another's), is inside its panel | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `panels:panel_fields.window:argument:4@9776f952` | argument | 4 | b = np.asarray(list(owned[d]) + ([names[d]] if d in names else []), float).reshape(-1, 4) | convention | an array's shape | membrane, perf | - |

## provenance

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `provenance:stuck:fallback-marker` | fallback-marker | 2 | out = [f"{l.a_id}>{l.b_id}" for l in loops if l.fallback] | unknown | not a fallback: the notes' list of loop_curve and na_* fallbacks | copies, mcp, na_curve | - |
| `provenance:over_labels:fallback-marker` | fallback-marker | 1 | return [f"{l.a_id}>{l.b_id}" for l in loops if l.crosses_labels and not l.fallback] | unknown | not a fallback: the notes' list of loop_over_text | app_jobs | - |

## render

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `render:_HEAD_HALF` | constant | HEAD / 2 | _HEAD_HALF = HEAD / 2 | geometry | HEAD / 2; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_CORNER` | constant | 0.3 | _CORNER = 0.3 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_MEASURE` | constant | FontProperties(family=_FONTS, size=10) | _MEASURE = FontProperties(family=_FONTS, size=10) # text widths are measured at this size and scaled | convention | "text widths are measured at this size and scaled" (the result does not depend on it); earliest in this file: d306041 2026-10-10 fix: text uses foldmap's fonts however the figure is saved; legend measured by a… | - | - |
| `render:_LIPID_ALPHA` | constant | 0.8 | _LIPID_ALPHA = 0.8 # the band over a sheet's plate lets the plate show through | unknown | earliest in this file: 42c3147 2026-10-13 C-MEM-BARREL fix 1: band over the sheet plates; only defining elements band a row | - | - |
| `render:_SHORT_HELIX` | constant | 1.6 | _SHORT_HELIX = 1.6 # helices/3-10 boxes shorter than this hold their residue numbers past the ends | unknown | earliest in this file: 1e7b1f3 2026-10-09 Fix overlapping residue numbers on short elements and the sstype legend | - | - |
| `render:LONG_LOOP` | constant | 25 | LONG_LOOP = 25 # loops at least this many residues get a residue-count label (chain ends: layout.LONG_TAIL) | unknown | earliest in this file: 546bcfb 2026-10-09 Label long loops and chain ends with their residue count | - | - |
| `render:_DOMAIN_BAR_H` | constant | 1.1 | _DOMAIN_BAR_H = 1.1 # page units each chain's domain bar takes below the figure (at least; more for larger text) | unknown | earliest in this file: bee6cf8 2026-10-11 Domain colouring and a linear domain bar when domains are given | - | - |
| `render:Z` | constant | { # the drawing order, back to front (plan 10 R8: panel < plate < membrane band < loops < bodies < bars/tethers # < lab… | Z = { # the drawing order, back to front (plan 10 R8: panel < plate < membrane band < loops < bodies < bars/tethers | convention | drawing order; earliest in this file: df8de19 2026-10-09 Add ligand and metal ion markers | - | - |
| `render:_BAR_TICK` | constant | 0.3 | _BAR_TICK = 0.3 # page units from a narrow segment to its name below it | unknown | earliest in this file: bee6cf8 2026-10-11 Domain colouring and a linear domain bar when domains are given | - | - |
| `render:_LEGEND_ROW` | constant | 1.0 | _LEGEND_ROW = 1.0 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_LOOP_ROOM` | constant | 0.25 | _LOOP_ROOM = 0.25 # page units kept between a loop point and the page edge: the casing at the largest line width | unknown | earliest in this file: 23c51a3 2026-10-11 The page holds every loop; a single-chain figure has no chain legend | - | - |
| `render:_RES_PER_TURN` | constant | 3.6 | _RES_PER_TURN = 3.6 | biology | residues per α-helix turn; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_MIN_PITCH` | constant | 0.5 | _MIN_PITCH = 0.5 # page units; an end-on helix keeps legible turns instead of a dark smear | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_BAND` | constant | 0.36 | _BAND = 0.36 # ribbon width along the axis, as a fraction of one turn | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | route | - |
| `render:_SHINE` | constant | 0.4 | _SHINE = 0.4 # highlight width, as a fraction of the ribbon | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_DNA_SEG` | constant | 16 | _DNA_SEG = 16 # points along each half-turn | convention | sampling points; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_TIE_LABEL` | constant | 0.75 | _TIE_LABEL = 0.75 # a tie label's size: a base letter's (0.75 of the text size) | geometry | "a base letter's (0.75 of the text size)"; earliest in this file: 57d8694 2026-10-10 feat: nucleic-acid backbones in their own colour by modality, chain or override;… | - | - |
| `render:_CURVE_MAX` | constant | 3.0 | _CURVE_MAX = 3.0 # longest stretch of a leg a curved corner may take | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_T` | constant | np.linspace(0.0, 1.0, 25)[:, None] | _T = np.linspace(0.0, 1.0, 25)[:, None] | convention | sampling points; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_RAMP` | constant | {"turbo": (0.1, 0.92), "viridis": (0.0, 0.92), "plasma": (0.0, 0.9), "cividis": (0.0, 1.0)} | _RAMP = {"turbo": (0.1, 0.92), "viridis": (0.0, 0.92), "plasma": (0.0, 0.9), "cividis": (0.0, 1.0)} # skip muddy ends | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:_NAME_CONTRAST` | constant | 4.5 | _NAME_CONTRAST = 4.5 # WCAG AA for normal text: a domain's name against what it is drawn on | convention | "WCAG AA for normal text"; earliest in this file: 96a55f3 2026-10-12 Domain names keep a WCAG contrast of 4.5 on what they sit on | - | - |
| `render:_RESIDUE_STEP` | constant | 0.5 | _RESIDUE_STEP = 0.5 # page units along a loop's route per residue from the nearer element end (R3) | unknown | earliest in this file: 2d7344d 2026-10-16 Plan 12 T4: a loop residue's mark sits by residue distance from the nearer eleme… | - | - |
| `render:_SS_MIN` | constant | 0.5 | _SS_MIN = 0.5 # shortest disulfide bar drawn (page units) | unknown | earliest in this file: 20f0723 2026-10-09 Interactive selection; shorter disulfide bars; README on how it works | - | - |
| `render:_SS_PARALLEL` | constant | np.cos(np.radians(10)) | _SS_PARALLEL = np.cos(np.radians(10)) # bars closer in direction than this can lie on top of each other | unknown | earliest in this file: d29b7e2 2026-10-10 Disulfides: bars lying on top of each other detour into separate lanes as bracke… | - | - |
| `render:_SS_CLEAR` | constant | 2.0 | _SS_CLEAR = 2.0 # bars nearer than this many bar widths (side by side) cannot be told apart | unknown | earliest in this file: d29b7e2 2026-10-10 Disulfides: bars lying on top of each other detour into separate lanes as bracke… | - | - |
| `render:_SS_STEP` | constant | 3.0 | _SS_STEP = 3.0 # lane spacing of a detouring bar: at least this many bar widths... | unknown | earliest in this file: d29b7e2 2026-10-10 Disulfides: bars lying on top of each other detour into separate lanes as bracke… | - | - |
| `render:_SS_DOT_STEP` | constant | 2.5 | _SS_DOT_STEP = 2.5 # ...and this many sulfur-dot radii, so the legs always stand out from a dot | unknown | earliest in this file: d0ac28d 2026-10-10 Disulfides: lane step is at least three bar widths and 2.5 dot radii, so a brack… | - | - |
| `render:_LIGAND_TETHERS` | constant | 4 | _LIGAND_TETHERS = 4 # a ligand is tethered to the elements it touches most, at most this many | unknown | earliest in this file: df8de19 2026-10-09 Add ligand and metal ion markers | - | - |
| `render:_TETHER_REACH` | constant | 6.5 | _TETHER_REACH = 6.5 # page units: a tether target farther than this from the heaviest one is another pocket | unknown | earliest in this file: 8b5b8bb 2026-10-11 Ligand tethers go to the ligand's own chain and stay within reach of its pocket … | ligands | - |
| `render:_SPOT_RADIUS` | constant | 3.0 | _SPOT_RADIUS = 3.0 # the spot search stops at this ring radius | unknown | earliest in this file: 8b5b8bb 2026-10-11 Ligand tethers go to the ligand's own chain and stay within reach of its pocket … | ligands, text_clear | - |
| `render:_SPOT_FAR` | constant | 3 | _SPOT_FAR = 3 # ... and, when nothing is clear within it, at this many times it (then the mark goes on its anchor) | named-entry:7QUP | git: b98171e's message: "the search widens to _SPOT_FAR (3) times the radius ... 7QUP: GTP:9E500 marked again" | ligands, text_clear | - |
| `render:_SPOT_MORE_RINGS` | constant | 4 | _SPOT_MORE_RINGS = 4 # rings past the first clear spot's where a spot whose tethers cross fewer bodies is looked for | named-entry:7LX9,9YLM,3MHT,1GOX | git: ac828a3's message: "_mark_spot looks up to 4 rings further ... 7LX9 1.15 -> 0, 3MHT 0.56 -> 0, 9YLM 3.15 -> 2.5, 1GOX-A 12.0 -> 11.3" | - | - |
| `render:_ALONG_NEAR` | constant | 0.15 | _ALONG_NEAR, _ALONG_TURN = ( | unknown | earliest in this file: 2317e7e 2026-10-15 C9 (random-100, batch C2a): a ligand's spot counts the length its tethers run al… | - | - |
| `render:_ALONG_TURN` | constant | np.radians(20) | _ALONG_NEAR, _ALONG_TURN = ( | unknown | earliest in this file: 2317e7e 2026-10-15 C9 (random-100, batch C2a): a ligand's spot counts the length its tethers run al… | - | - |
| `render:_RAIL_GRAZE` | constant | 0.3 | _RAIL_GRAZE = 0.3 # page units a tether may run inside a stem past the entry its target needs and not count as crossing | unknown | earliest in this file: 7e0f1e8 2026-10-17 A ligand held by nucleotides alone judges every clear spot by the stems its teth… | - | - |
| `render:_AXIS_MARK_GAP` | constant | 0.12 | _AXIS_MARK_GAP = 0.12 # page units between two marks of a rosette's axis cluster | unknown | earliest in this file: 3eb60d0 2026-10-20 Rosette routes, labels, ligands and checks (cyclic=rosette, opt-in): protomer 0'… | - | - |
| `render:_LINKED_STAR` | constant | 0.38 | _LINKED_STAR = 0.38 # radius of the mark of a residue fused into the chain | unknown | earliest in this file: 4fda3ae 2026-10-11 Peptide-linked residues without N/CA/C names (GFP's CRO66) stay in the chain and… | - | - |
| `render:_LINKED_LABEL` | constant | 0.6 | _LINKED_LABEL = 0.6 # its label, as a fraction of the font | unknown | earliest in this file: 4fda3ae 2026-10-11 Peptide-linked residues without N/CA/C names (GFP's CRO66) stay in the chain and… | - | - |
| `render:_LINKED_LABEL_MIN` | constant | 7.0 | _LINKED_LABEL_MIN = 7.0 # but never smaller than this (pt) | unknown | earliest in this file: 4505c9f 2026-10-11 Task 4 review follow-ups: guard test, 1EMA helix pinned, linked-residue star lar… | - | - |
| `render:_HIDE_STEP` | constant | 0.02 | _HIDE_STEP = 0.02 # page units: a bar or a tether is sampled this finely where a body or a text hides it | unknown | earliest in this file: 36859e8 2026-10-13 Disulfide bars and ligand tethers pass under labels; bond ends marked on the str… | - | - |
| `render:_SS_PIECE_MIN` | constant | 0.15 | _SS_PIECE_MIN = 0.15 # page units: a stretch of a disulfide bar shorter than this between bodies or text is not drawn | unknown | earliest in this file: 8ac832b 2026-10-13 Bars and tethers hide under the bodies as drawn: closed outlines, one set of bar… | - | - |
| `render:KYTE_DOOLITTLE` | constant | { "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": … | KYTE_DOOLITTLE = { | biology | the Kyte-Doolittle hydropathy scale; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | style | - |
| `render:_PROPERTY` | constant | { "bfactor": ("RdYlBu_r", 0.06, 0.94, "B-factor: rigid → flexible"), "hydropathy": ("BrBG_r", 0.08, 0.92, "hydropathy: … | _PROPERTY = { | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `render:PLDDT_BANDS` | constant | ( # AlphaFold's own confidence colours: (lower bound, colour, legend text) (90.0, "#0053d6", "very high (> 90)"), (70.0… | PLDDT_BANDS = ( # AlphaFold's own confidence colours: (lower bound, colour, legend text) | convention | "AlphaFold's own confidence colours"; earliest in this file: 81f0592 2026-10-09 Add AlphaFold pLDDT colouring and residue-coloured loops | - | - |
| `render:_DIVERGING` | constant | (0.04, 0.96) | _DIVERGING = (0.04, 0.96) # the part of the diverging map used: its darkest ends left out | unknown | earliest in this file: e1c36ea 2026-10-21 Drawing the researcher's data (color_by=data, plan 13d-5): numbers as a stripe p… | - | - |
| `render:_DATA_MARK` | constant | 0.15 | _DATA_MARK = 0.15 # page units: a category mark's radius ... | unknown | earliest in this file: e1c36ea 2026-10-21 Drawing the researcher's data (color_by=data, plan 13d-5): numbers as a stripe p… | - | - |
| `render:_DATA_RING` | constant | 0.25 | _DATA_RING = 0.25 # ... and its ring's | unknown | earliest in this file: e1c36ea 2026-10-21 Drawing the researcher's data (color_by=data, plan 13d-5): numbers as a stripe p… | - | - |
| `render:INTERACTION_INK` | constant | { # kind -> (colour, dash) of its lines (plan 13c-4). H-bonds keep look.hbond_color (None here) # and their short dashe… | INTERACTION_INK = { # kind -> (colour, dash) of its lines (plan 13c-4). H-bonds keep look.hbond_color (None here) | convention | colours and dashes; earliest in this file: 96aebc8 2026-10-21 Style keys interactions, interactions_between, interaction_max; hbonds=na unchan… | hbonds, render | - |
| `render:XL_WIDTH` | constant | 1.4 | XL_WIDTH = 1.4 # a crosslink line is this many interaction lines wide: the researcher's own data, read first | unknown | earliest in this file: 243697e 2026-10-21 Crosslinks drawn, the app's data upload and MCP data parameters (plan 13d-6): ea… | - | - |
| `render:INTERACTION_GAP` | constant | 2.0 | INTERACTION_GAP = 2.0 # line widths between the parallel lines of two kinds drawn between one residue pair | unknown | earliest in this file: 79c95f1 2026-10-21 interaction_overlay: salt bridges, stacking, cation-pi and hydrophobic lines (ia… | hbonds | - |
| `render:_OVER_ALPHA` | constant | 0.6 | _OVER_ALPHA = 0.6 # its fill and edge: the strands beneath show through | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `render:_OVER_TONE` | constant | 0.7 | _OVER_TONE = 0.7 # coloured as its strands: darkened this much, to stand apart from them ... | unknown | earliest in this file: 2473be0 2026-10-11 A helix over a sheet of its own colour takes a darker tone; its name keeps off a… | - | - |
| `render:_OVER_SAME_ALPHA` | constant | 0.8 | _OVER_SAME_ALPHA = 0.8 # ... its fill at this alpha, its edge (darker still) opaque | unknown | earliest in this file: 2473be0 2026-10-11 A helix over a sheet of its own colour takes a darker tone; its name keeps off a… | - | - |
| `render:_COUNT_GAP` | constant | 0.15 | _COUNT_GAP = 0.15 # page units from a loop's line to its residue count's ink | unknown | introduced by 0f268ac (C3: 9ON5, 7BRG, 2B9V), whose message does not tie this value to those entries; earliest in this file: 0f268ac 2026-10-13 C3 (random-100): loop counts placed clear of elements, lines and text, inside th… | - | - |
| `render:_COUNT_RUNS` | constant | 4 | _COUNT_RUNS = 4 # the loop's longest straight runs tried first for its count | unknown | introduced by 0f268ac (C3: 9ON5, 7BRG, 2B9V), whose message does not tie this value to those entries; earliest in this file: 0f268ac 2026-10-13 C3 (random-100): loop counts placed clear of elements, lines and text, inside th… | text_clear | - |
| `render:_COUNT_AT` | constant | (0.5, 1 / 3, 2 / 3) | _COUNT_AT = (0.5, 1 / 3, 2 / 3) # where along a run its count is tried: its middle and thirds ... | unknown | earliest in this file: 340f18f 2026-10-16 Counts never drawn over ink: loop counts look along every run, a count with no c… | text_clear | - |
| `render:_COUNT_AT_ENDS` | constant | (1 / 6, 5 / 6) | _COUNT_AT_ENDS = (1 / 6, 5 / 6) # ... then nearer its ends | unknown | earliest in this file: 340f18f 2026-10-16 Counts never drawn over ink: loop counts look along every run, a count with no c… | - | - |
| `render:_SIDE_TRIES` | constant | 3 | _SIDE_TRIES = 3 # side names tried this many label heights out from the face, at the band's left end, then its right | unknown | introduced by 0a35b86 (plan 26: 5LDW), whose message ties 5LDW to _PAST, not this value; earliest in this file: 0a35b86 2026-10-24 Task 2 fix round 1: face rows, columns in order, 5LDW crossing, side names off l… | - | - |
| `render:_NAME_PLATE_ALPHA` | constant | 0.85 | _NAME_PLATE_ALPHA = 0.85 # a domain name's plate: the page colour at this alpha over the panel's edge | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `render:_NAME_PLATE_PAD` | constant | 0.06 | _NAME_PLATE_PAD = 0.06 # page units the plate reaches beyond the name's box | unknown | earliest in this file: 498a910 2026-10-12 Domain panels follow their elements and keep apart; names on clear plates (4OO8,… | - | - |
| `render:_COPY_DASHES` | constant | ("-", (0, (5, 2)), (0, (1.5, 1.5)), (0, (5, 2, 1.5, 2))) | _COPY_DASHES = ("-", (0, (5, 2)), (0, (1.5, 1.5)), (0, (5, 2, 1.5, 2))) # copies of one domain: one edge each | convention | dash patterns; earliest in this file: 7ec7e32 2026-10-12 Domains stay visible: tinted or edged sheet plates, a copy's own edge, leader na… | - | - |
| `render:ABSENT_DASH` | constant | (0, (3, 2)) | ABSENT_DASH = (0, (3, 2)) # ... dash ... | convention | a dash pattern; earliest in this file: 6380df1 2026-10-21 The diff drawn (plan 13e-2): foldmap diff A B -o diff.svg [--json] [--domain NAM… | - | - |
| `render:ABSENT_LW` | constant | 0.9 | ABSENT_LW = 0.9 # ... this many line widths wide | unknown | earliest in this file: 6380df1 2026-10-21 The diff drawn (plan 13e-2): foldmap diff A B -o diff.svg [--json] [--domain NAM… | - | - |
| `render:PNG_DPI` | constant | 300 | PNG_DPI = 300 # a PNG's resolution ... | convention | a PNG's resolution; earliest in this file: a3cf24d 2026-10-15 C4 (random-100, batch C1): text has a 6 pt floor and the page grows instead; leg… | - | - |
| `render:MAX_PNG_PIXELS` | constant | 50_000_000 | MAX_PNG_PIXELS = 50_000_000 # ... up to this many pixels (200 MB to draw): a larger page's PNG takes a lower dpi | convention | "200 MB to draw"; earliest in this file: a3cf24d 2026-10-15 C4 (random-100, batch C1): text has a 6 pt floor and the page grows instead; leg… | app_jobs, page_text | - |
| `render:_band:default:n` | default | 12 | def _band(t0: float, t1: float, n: int = 12, taper: int = 0, width: float = _BAND) -> np.ndarray: | unknown | - | membrane, plddt | - |
| `render:_turn:argument:-0.25@5bebbe12` | argument | -0.25 | _band(-0.25, 0.25, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:0.25@5bebbe12` | argument | 0.25 | _band(-0.25, 0.25, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:0.25@a3e0d802` | argument | 0.25 | _band(0.25, 0.75, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:0.75@a3e0d802` | argument | 0.75 | _band(0.25, 0.75, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:0.25@a4d70d51` | argument | 0.25 | _band(0.0, 0.25, 9, taper=+1, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:9@a4d70d51` | argument | 9 | _band(0.0, 0.25, 9, taper=+1, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:-0.25@ad5390b5` | argument | -0.25 | _band(-0.25, 0.0, 9, taper=-1, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:9@ad5390b5` | argument | 9 | _band(-0.25, 0.0, 9, taper=-1, width=band), | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:-0.16@0cbc2b1a` | argument | -0.16 | _band(-0.16, 0.16, 14, taper=2, width=band * _SHINE) + (band * 0.1, 0), # glint on the near face | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:0.16@0cbc2b1a` | argument | 0.16 | _band(-0.16, 0.16, 14, taper=2, width=band * _SHINE) + (band * 0.1, 0), # glint on the near face | unknown | - | diff, like, na_layout | - |
| `render:_turn:argument:14@0cbc2b1a` | argument | 14 | _band(-0.16, 0.16, 14, taper=2, width=band * _SHINE) + (band * 0.1, 0), # glint on the near face | unknown | - | diff, like, na_layout | - |
| `render:_dna_turns:comparison:1e-6@8483b1e8` | comparison | 1e-6 | if b - a > 1e-6: | convention | a floating-point tolerance | na_text_keepoff | - |
| `render:_rung_lines:argument:1e-9@b4b83675` | argument | 1e-9 | v = (q - p) / max(float(np.hypot(*(q - p))), 1e-9) | convention | a floating-point tolerance | - | - |
| `render:_rung_lines:argument:1.4@61cc0b8b` | argument | 1.4 | lw=ink.lw * 1.4, | unknown | - | - | - |
| `render:_rung_lines:argument:1.5@24372526` | argument | 1.5 | ls=(0, (1, 1.5)) if pair in ink.gu else "-", | unknown | - | - | - |
| `render:_ladder:argument:0.3@96e0981f` | argument | 0.3 | np.array([d.point(key, 0.3 if key in d.bulges else 0.0) for key in mine]) | unknown | - | - | - |
| `render:_ladder:argument:2.2@b066c232` | argument | 2.2 | lw=ink.lw * 2.2 * weight, | unknown | - | - | - |
| `render:_draw_stem_label:argument:0.8@689a8d37` | argument | 0.8 | *text_anchor(d.name, True, stem_label_pos(d), ink.em * 0.8, _upright(d.angle)), | unknown | - | - | - |
| `render:_draw_stem_label:argument:0.8@d613b6a1` | argument | 0.8 | fontsize=ink.font * 0.8, | unknown | - | - | - |
| `render:_ribbon:argument:0.8@c609e7ea` | argument | 0.8 | lw=ink.lw * 0.8, | unknown | - | - | - |
| `render:_ribbon:argument:2.2@156c968f` | argument | 2.2 | _compound(rungs, closed=False), fc="none", ec=ink.rung, lw=ink.lw * 2.2, capstyle="round", zorder=Z["na-rung"] | unknown | - | - | - |
| `render:_block:argument:1.2@7beb4240` | argument | 1.2 | lw=ink.lw * 1.2 * max(ink.weight(key[0]) for key in d.axial), | unknown | - | na_text_keepoff | - |
| `render:_block:argument:0.9@3c190eb3` | argument | 0.9 | fontsize=ink.font * 0.9, | unknown | - | na_text_keepoff | - |
| `render:_block:argument:0.75@1fe018cb` | argument | 0.75 | *text_anchor(first, False, stem_label_pos(d), ink.em * 0.75, _upright(d.angle)), | unknown | - | na_text_keepoff | - |
| `render:_block:argument:0.75@248ad675` | argument | 0.75 | fontsize=ink.font * 0.75, # a base letter's size: the stem's name slot, empty in this style | unknown | - | na_text_keepoff | - |
| `render:_letters:argument:0.3@35f266b2` | argument | 0.3 | _rung_lines(ink, k, d, inset=0.3) | unknown | - | - | - |
| `render:_curve_rails:argument:2.2@b066c232` | argument | 2.2 | lw=ink.lw * 2.2 * weight, | unknown | - | - | - |
| `render:_curve_rungs:argument:1e-9@b4b83675` | argument | 1e-9 | v = (q - p) / max(float(np.hypot(*(q - p))), 1e-9) | convention | a floating-point tolerance | - | - |
| `render:_curve_rungs:argument:1.4@61cc0b8b` | argument | 1.4 | lw=ink.lw * 1.4, | unknown | - | - | - |
| `render:_curve_rungs:argument:1.5@24372526` | argument | 1.5 | ls=(0, (1, 1.5)) if pair in ink.gu else "-", | unknown | - | - | - |
| `render:_curve_letters:argument:0.3@f2367497` | argument | 0.3 | _curve_rungs(ink, g, c, inset=0.3) | unknown | - | - | - |
| `render:_curve_turns:comparison:1e-6@8483b1e8` | comparison | 1e-6 | if b - a > 1e-6: | convention | a floating-point tolerance | - | - |
| `render:_curve_ribbon:argument:0.8@c609e7ea` | argument | 0.8 | lw=ink.lw * 0.8, | unknown | - | - | - |
| `render:_curve_ribbon:argument:2.2@156c968f` | argument | 2.2 | _compound(rungs, closed=False), fc="none", ec=ink.rung, lw=ink.lw * 2.2, capstyle="round", zorder=Z["na-rung"] | unknown | - | - | - |
| `render:_curve_block:argument:1.2@d83afa85` | argument | 1.2 | lw=ink.lw * 1.2 * max(ink.weight(key[0]) for key in c.at), | unknown | - | - | - |
| `render:_curve_block:argument:0.9@3c190eb3` | argument | 0.9 | fontsize=ink.font * 0.9, | unknown | - | - | - |
| `render:_curve_block:argument:0.75@c9128220` | argument | 0.75 | *text_anchor(first, False, c.name_xy, ink.em * 0.75, rot), | unknown | - | - | - |
| `render:_curve_block:argument:0.75@248ad675` | argument | 0.75 | fontsize=ink.font * 0.75, | unknown | - | - | - |
| `render:_draw_curve_label:argument:0.8@689140e9` | argument | 0.8 | *text_anchor(c.name, True, c.name_xy, ink.em * 0.8, rot), | unknown | - | - | - |
| `render:_draw_curve_label:argument:0.8@d613b6a1` | argument | 0.8 | fontsize=ink.font * 0.8, | unknown | - | - | - |
| `render:_draw_strands:argument:1.5@e42592a8` | argument | 1.5 | ls=((0, (2, 2)) if path.label else (0, (1, 1.5))) if gap else "-", | unknown | - | - | - |
| `render:_draw_strands:argument:0.75@ea89db02` | argument | 0.75 | t = ink.ax.text(*xy, path.label, ha="center", va=va, fontsize=ink.font * 0.75, color=ink.front) | unknown | - | - | - |
| `render:_draw_na_marks:argument:1.5@52e8e8e8` | argument | 1.5 | ls=(0, (1, 1.5)) if wobble else "-", | unknown | - | na_layout | - |
| `render:_draw_na_marks:argument:1.2@df43c0ee` | argument | 1.2 | lw=ink.lw * 1.2, | unknown | - | na_layout | - |
| `render:_draw_na_marks:argument:0.75@2382bf18` | argument | 0.75 | *text_anchor(base, bool(hit), xy, ink.em * 0.75, rot), # centred on xy, as the layout measures it | unknown | - | na_layout | - |
| `render:_draw_na_marks:argument:0.75@248ad675` | argument | 0.75 | fontsize=ink.font * 0.75, | unknown | - | na_layout | - |
| `render:_draw_na_marks:argument:0.6@96c9594b` | argument | 0.6 | linewidths=ink.lw * 0.6, | unknown | - | na_layout | - |
| `render:strand_arrow:argument:0.45@9c7d05de` | argument | 0.45 | tip, hh = min(tip, length * 0.45), _HEAD_HALF * k | unknown | - | - | - |
| `render:_smooth_path:comparison:3@edf9e943` | comparison | 3 | if len(pts) < 3: | unknown | - | - | - |
| `render:_smooth_path:comparison:0.02@2d8e413b` | comparison | 0.02 | if r < 0.02 or _corner_clear(ctrl, rects): | unknown | - | - | - |
| `render:_shade:argument:0.55@84307f85` | argument | 0.55 | return tint(colour, 0.55 * (1 - 2 * t)) if t < 0.5 else darken(colour, 1 - 0.4 * (2 * t - 1)) | unknown | - | - | - |
| `render:_shade:argument:0.4@84307f85` | argument | 0.4 | return tint(colour, 0.55 * (1 - 2 * t)) if t < 0.5 else darken(colour, 1 - 0.4 * (2 * t - 1)) | unknown | - | - | - |
| `render:_relative_luminance:comparison:0.04045@98aefd01` | comparison | 0.04045 | lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(colour)] | unknown | - | - | - |
| `render:domain_name_colour:argument:0.75@611492b0` | argument | 0.75 | start = darken(tone, 0.75) if ink_on(on) == "#000000" else tone | unknown | - | domains | - |
| `render:domain_name_colour:argument:21@592ae84b` | argument | 21 | for k in range(21): | unknown | - | domains | - |
| `render:domain_name_colour:argument:20@e87c05dc` | argument | 20 | colour = _mix(start, ink_on(on), k / 20) | unknown | - | domains | - |
| `render:_label_ink:argument:0.8@8cd27593` | argument | 0.8 | return legible(darken(colour, 0.8)) | unknown | - | - | - |
| `render:legible:argument:12@42886ae7` | argument | 12 | for _ in range(12): | unknown | - | blueprint, style | - |
| `render:legible:comparison:0.55@ef7eef24` | comparison | 0.55 | if _luminance(colour) >= 0.55: | unknown | - | blueprint, style | - |
| `render:legible:argument:0.2@57f75c5c` | argument | 0.2 | colour = _mix(to_hex(to_rgb(colour)), "#ffffff", 0.2) | unknown | - | blueprint, style | - |
| `render:legible:argument:12@42886ae7#2` | argument | 12 | for _ in range(12): | unknown | - | blueprint, style | - |
| `render:legible:comparison:0.2126@12c966fc` | comparison | 0.2126 | if 0.2126 * r + 0.7152 * g + 0.0722 * b <= ceiling: | unknown | - | blueprint, style | - |
| `render:legible:comparison:0.7152@12c966fc` | comparison | 0.7152 | if 0.2126 * r + 0.7152 * g + 0.0722 * b <= ceiling: | unknown | - | blueprint, style | - |
| `render:legible:comparison:0.0722@12c966fc` | comparison | 0.0722 | if 0.2126 * r + 0.7152 * g + 0.0722 * b <= ceiling: | unknown | - | blueprint, style | - |
| `render:legible:argument:0.85@fe40df72` | argument | 0.85 | colour = darken(to_hex((r, g, b)), 0.85) | unknown | - | blueprint, style | - |
| `render:residue_point:argument:1e-9@28829227` | argument | 1e-9 | into = (target - (run[k] - seg[k])) / max(seg[k], 1e-9) | convention | a floating-point tolerance | hbonds, links, render, science | - |
| `render:_nearest_on:argument:1e-12@2423a05e` | argument | 1e-12 | t = float(np.clip(np.dot(q - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0)) | convention | a floating-point tolerance | - | - |
| `render:disulfide_points:argument:3@f8579048` | argument | 3 | for _ in range(3): # both in loops: settle by alternating | unknown | - | links | - |
| `render:disulfide_side:comparison:1e-9@03be3852` | comparison | 1e-9 | if abs(n[1]) > 1e-9: | convention | a floating-point tolerance | links | - |
| `render:disulfide_side:comparison:1e-9@06bc9cab` | comparison | 1e-9 | return below if abs(less) < 1e-9 else (1 if less < 0 else -1) | convention | a floating-point tolerance | links | - |
| `render:_nearest_beyond:argument:30@4158c5d9` | argument | 30 | dense = np.vstack([a + t * (b - a) for a, b in zip(pts, pts[1:]) for t in np.linspace(0, 1, 30)]) | unknown | - | - | - |
| `render:_symbol:argument:0.6@9d283a1f` | argument | 0.6 | patch.set(facecolor=colour, edgecolor=legible("#2b2b2b"), linewidth=lw * 0.6, zorder=Z["mark"]) | unknown | - | - | - |
| `render:_tether_targets:comparison:0.3@5b50eb64` | comparison | 0.3 | if at is not None and all(np.hypot(at[0] - t[0], at[1] - t[1]) > 0.3 for t in targets): | unknown | - | - | - |
| `render:_local_targets:comparison:0.3@dbff7687` | comparison | 0.3 | np.hypot(t[0] - u[0], t[1] - u[1]) > 0.3 for u in targets | unknown | - | ligands | - |
| `render:_over_bodies:argument:4@e25b3f3d` | argument | 4 | boxes = np.asarray(bodies, float).reshape(-1, 4) | convention | an array's shape | ligands, na_layout | - |
| `render:_along_loops:argument:4@f02fb400` | argument | 4 | segs = np.asarray(runs, float).reshape(-1, 4) | convention | an array's shape | - | - |
| `render:_along_loops:comparison:1e-9@0e109098` | comparison | 1e-9 | if length < 1e-9: | convention | a floating-point tolerance | - | - |
| `render:_along_loops:argument:1e-12@89bd794a` | argument | 1e-12 | sin = np.abs(u[0] * e[:, 1] - u[1] * e[:, 0]) / np.maximum(el, 1e-12) | convention | a floating-point tolerance | - | - |
| `render:_along_loops:comparison:1e-9@51746fed` | comparison | 1e-9 | beside = (el > 1e-9) & (sin < np.sin(_ALONG_TURN)) & (off_a <= _ALONG_NEAR) & (off_b <= _ALONG_NEAR) | convention | a floating-point tolerance | - | - |
| `render:_strands_crossed:comparison:-1e-12@ed07140b` | comparison | -1e-12 | if ((turn(p, q, a) * turn(p, q, b) < -1e-12) & (turn(a, b, p) * turn(a, b, q) < -1e-12)).any(): | convention | a floating-point tolerance | na_layout | - |
| `render:_strands_crossed:comparison:-1e-12@ed07140b#2` | comparison | -1e-12 | if ((turn(p, q, a) * turn(p, q, b) < -1e-12) & (turn(a, b, p) * turn(a, b, q) < -1e-12)).any(): | convention | a floating-point tolerance | na_layout | - |
| `render:_mark_spot:argument:0.35@762e4173` | argument | 0.35 | for n, ring in enumerate(np.arange(0.0, radius, 0.35)): | unknown | line last changed by b98171e C9 (random-100, follow-up): a ligand with nothing clear near it is marked further out, else on its anchor (7Q… (git blame; evidence only) | na_layout | - |
| `render:_mark_spot:argument:10@277fa249` | argument | 10 | for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False): | unknown | line last changed by cae6ae1 C5 (random-100): ligand marks keep off loop lines, counts, earlier tethers and names (1BXI, 1A0Z, 5SB5, 2HHB … (git blame; evidence only) | na_layout | - |
| `render:_mark_spot:comparison:1e-9@4097f1e5` | comparison | 1e-9 | if overlaps_any(box, taken) or any(clip_len(a, b, box) > 1e-9 for a, b in tethers): | convention | a floating-point tolerance | na_layout | - |
| `render:_mark_spot:comparison:1e-9@c8844eae` | comparison | 1e-9 | clip_len(c, t, b) > 1e-9 for t in ts for b in texts | convention | a floating-point tolerance | na_layout | - |
| `render:_mark_spot:argument:4@9679854b` | argument | 4 | round(over, 4), | convention | rounding digits | na_layout | - |
| `render:_mark_spot:argument:4@0843070d` | argument | 4 | round(length, 4), | convention | rounding digits | na_layout | - |
| `render:_mark_spot:comparison:1e-9@59447e8d` | comparison | 1e-9 | if best is None or over < best[0] - 1e-9: | convention | a floating-point tolerance | na_layout | - |
| `render:_mark_spot:comparison:1e-9@e03c15f2` | comparison | 1e-9 | if over <= 1e-9: | convention | a floating-point tolerance | na_layout | - |
| `render:_mark_spot:comparison:1e-9@9091e225` | comparison | 1e-9 | if best is not None and na is None and best[0] <= 1e-9: | convention | a floating-point tolerance | na_layout | - |
| `render:_least_bad_spot:argument:4@0dc5edc0` | argument | 4 | body = np.asarray(bodies, float).reshape(-1, 4) | convention | an array's shape | - | - |
| `render:_least_bad_spot:argument:4@53a718a6` | argument | 4 | other = np.asarray(taken, float).reshape(-1, 4) | convention | an array's shape | - | - |
| `render:_least_bad_spot:argument:0.35@812d8b4f` | argument | 0.35 | for ring in np.arange(0.0, _SPOT_FAR * _SPOT_RADIUS, 0.35): | unknown | - | - | - |
| `render:_least_bad_spot:argument:10@277fa249` | argument | 10 | for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False): | unknown | - | - | - |
| `render:_na_targets:comparison:0.3@7c311ec1` | comparison | 0.3 | if all(np.hypot(x - t[0], y - t[1]) > 0.3 for t in targets): | unknown | - | - | - |
| `render:_axis_cluster:argument:6@3d5a75b0` | argument | 6 | key=lambda s: (round(float(np.hypot(*s)), 6), round(float(np.arctan2(s[1], s[0]) % (2 * np.pi)), 6)), | convention | rounding digits | rosette | - |
| `render:_axis_cluster:argument:6@3d5a75b0#2` | argument | 6 | key=lambda s: (round(float(np.hypot(*s)), 6), round(float(np.arctan2(s[1], s[0]) % (2 * np.pi)), 6)), | convention | rounding digits | rosette | - |
| `render:ligand_marks:argument:0.12@3e6c57d3` | argument | 0.12 | centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25) | unknown | line last changed by cae6ae1 C5 (random-100): ligand marks keep off loop lines, counts, earlier tethers and names (1BXI, 1A0Z, 5SB5, 2HHB … (git blame; evidence only) | ligands, na_layout, na_text_keepoff, text_clear | - |
| `render:ligand_marks:argument:0.1@3e6c57d3` | argument | 0.1 | centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25) | unknown | line last changed by cae6ae1 C5 (random-100): ligand marks keep off loop lines, counts, earlier tethers and names (1BXI, 1A0Z, 5SB5, 2HHB … (git blame; evidence only) | ligands, na_layout, na_text_keepoff, text_clear | - |
| `render:ligand_marks:argument:0.25@3e6c57d3` | argument | 0.25 | centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25) | unknown | line last changed by cae6ae1 C5 (random-100): ligand marks keep off loop lines, counts, earlier tethers and names (1BXI, 1A0Z, 5SB5, 2HHB … (git blame; evidence only) | ligands, na_layout, na_text_keepoff, text_clear | - |
| `render:_draw_ligands:argument:0.7@bc69ecf2` | argument | 0.7 | keep = link_keep(lw * 0.7, pt_per_unit, ax.figure.dpi) | unknown | - | - | - |
| `render:_draw_ligands:argument:0.7@1a10d27d` | argument | 0.7 | lw=lw * 0.7, | unknown | - | - | - |
| `render:_draw_ligands:argument:1.2@b61fa496` | argument | 1.2 | ls=(0, (1.2, 1.4)), | unknown | - | - | - |
| `render:_draw_ligands:argument:1.4@b61fa496` | argument | 1.4 | ls=(0, (1.2, 1.4)), | unknown | - | - | - |
| `render:_draw_ligands:argument:0.7@8da1d427` | argument | 0.7 | marker = Circle((cx, cy), hw, fc=fill, ec=darken(fill, 0.7), lw=edge, zorder=Z["mark"]) | unknown | - | - | - |
| `render:_draw_ligands:argument:0.7@d831aa66` | argument | 0.7 | ec=darken(own, 0.7) if own else LIGAND_EDGE, | unknown | - | - | - |
| `render:_draw_ligands:argument:0.62@52ee1c7e` | argument | 0.62 | fontsize=font * (0.62 if g.metal else 0.6), | unknown | - | - | - |
| `render:_draw_ligands:argument:0.6@52ee1c7e` | argument | 0.6 | fontsize=font * (0.62 if g.metal else 0.6), | unknown | - | - | - |
| `render:linked_marks:argument:0.1@098b55c5` | argument | 0.1 | for ring in np.arange(_LINKED_STAR + 0.1, 12.0, 0.35): | unknown | - | text_clear | - |
| `render:linked_marks:argument:12.0@098b55c5` | argument | 12.0 | for ring in np.arange(_LINKED_STAR + 0.1, 12.0, 0.35): | unknown | - | text_clear | - |
| `render:linked_marks:argument:0.35@098b55c5` | argument | 0.35 | for ring in np.arange(_LINKED_STAR + 0.1, 12.0, 0.35): | unknown | - | text_clear | - |
| `render:linked_marks:argument:10@277fa249` | argument | 10 | for angle in np.linspace(0, 2 * np.pi, max(1, int(ring * 10)), endpoint=False): | unknown | - | text_clear | - |
| `render:linked_marks:comparison:1e-9@a8a612db` | comparison | 1e-9 | ) and not any(clip_len(p, q, box) > 1e-9 for p, q in segments): | convention | a floating-point tolerance | text_clear | - |
| `render:_star:argument:10@5ee26aa7` | argument | 10 | angles = np.pi / 2 + np.arange(10) * np.pi / 5 | unknown | - | - | - |
| `render:_star:argument:10@4dff427a` | argument | 10 | radii = np.where(np.arange(10) % 2 == 0, r, 0.4 * r) | unknown | - | - | - |
| `render:_star:argument:0.4@4dff427a` | argument | 0.4 | radii = np.where(np.arange(10) % 2 == 0, r, 0.4 * r) | unknown | - | - | - |
| `render:_draw_linked:argument:0.7@1a10d27d` | argument | 0.7 | lw=lw * 0.7, | unknown | - | - | - |
| `render:_draw_uniprot:argument:0.27@4d39e91e` | argument | 0.27 | halo = Circle((x, y), 0.27, fc=_BACKGROUND.get(), ec=legible("#1d1d1f"), lw=lw * 0.9, zorder=Z["feature-ring"]) | unknown | - | - | - |
| `render:_draw_uniprot:argument:0.9@4d39e91e` | argument | 0.9 | halo = Circle((x, y), 0.27, fc=_BACKGROUND.get(), ec=legible("#1d1d1f"), lw=lw * 0.9, zorder=Z["feature-ring"]) | unknown | - | - | - |
| `render:_draw_uniprot:argument:0.17@e797ea8a` | argument | 0.17 | dot = Circle((x, y), 0.17, fc=colour, ec="none", zorder=Z["feature"]) | unknown | - | - | - |
| `render:_in_boxes:argument:4@c905b135` | argument | 4 | b = np.asarray(boxes, float).reshape(-1, 4) | convention | an array's shape | - | - |
| `render:strand_label_spot:argument:0.45@13d2860f` | argument | 0.45 | room = p.length / 2 - min(HEAD_LEN, p.length * 0.45) - max(abs(y0), abs(y1)) # the shaft, either side | unknown | - | - | - |
| `render:strand_label_spot:argument:0.05@c2431173` | argument | 0.05 | for k in range(int(max(room, 0.0) / 0.05) + 1): | unknown | - | - | - |
| `render:strand_label_spot:comparison:1e-9@a06457c9` | comparison | 1e-9 | if not any(clip_len(u, v, box) > 1e-9 for q in lines for u, v in zip(q, q[1:])): | convention | a floating-point tolerance | - | - |
| `render:_draw_links:argument:2.4@152ba74b` | argument | 2.4 | keep = link_keep(lw * 2.4, pt_per_unit, ax.figure.dpi) | unknown | - | - | - |
| `render:_draw_links:argument:2.4@23103a85` | argument | 2.4 | lw=lw * 2.4, | unknown | - | - | - |
| `render:_draw_links:argument:5.5@23c1b757` | argument | 5.5 | s=(lw * 5.5) ** 2, | unknown | - | - | - |
| `render:_draw_links:argument:0.8@03368fc1` | argument | 0.8 | stem_halos = text_halos(ax, pt_per_unit, link_keep(lw * 0.8, pt_per_unit, ax.figure.dpi)) | unknown | - | - | - |
| `render:_draw_links:argument:0.12@3e6c57d3` | argument | 0.12 | centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25) | unknown | cites C28; line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:0.1@3e6c57d3` | argument | 0.1 | centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25) | unknown | cites C28; line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:0.25@3e6c57d3` | argument | 0.25 | centred_box(p.cx, p.cy, 0.12 * len(p.label) + 0.1, 0.25) | unknown | cites C28; line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:5@1d777e03` | argument | 5 | count = min(len(sugars), 5) | unknown | - | - | - |
| `render:_draw_links:argument:0.35@b6ba9042` | argument | 0.35 | for lead in np.arange(0.35, 6.0, 0.1): # the shortest stem that clears everything | unknown | line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:6.0@b6ba9042` | argument | 6.0 | for lead in np.arange(0.35, 6.0, 0.1): # the shortest stem that clears everything | unknown | line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:0.1@b6ba9042` | argument | 0.1 | for lead in np.arange(0.35, 6.0, 0.1): # the shortest stem that clears everything | unknown | line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:0.03@70519a75` | argument | 0.03 | centred_box(*(at + out * (lead + step * k)), size / 2 + 0.03, size / 2 + 0.03) | unknown | line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:argument:0.03@70519a75#2` | argument | 0.03 | centred_box(*(at + out * (lead + step * k)), size / 2 + 0.03, size / 2 + 0.03) | unknown | line last changed by 1167733 C28 (random-100): glycan stems cross no name; a linked residue's label place sized from its ink (1WWL, 8K8E) (git blame; evidence only) | - | - |
| `render:_draw_links:comparison:1e-9@945eb9fb` | comparison | 1e-9 | if clear_stem and any(clip_len(at, stem_at, b) > 1e-9 for b in texts): | convention | a floating-point tolerance | - | - |
| `render:_draw_links:argument:0.35@fd85d0ca` | argument | 0.35 | for lead in np.arange(0.35, 6.0, 0.25): | unknown | - | - | - |
| `render:_draw_links:argument:6.0@fd85d0ca` | argument | 6.0 | for lead in np.arange(0.35, 6.0, 0.25): | unknown | - | - | - |
| `render:_draw_links:argument:0.25@fd85d0ca` | argument | 0.25 | for lead in np.arange(0.35, 6.0, 0.25): | unknown | - | - | - |
| `render:_draw_links:argument:0.8@40c0e6ea` | argument | 0.8 | lw=lw * 0.8, | unknown | - | - | - |
| `render:_chevron:comparison:3@93deae6c` | comparison | 3 | if np.hypot(*(b - a)) < 3 * size: | unknown | - | - | - |
| `render:_chevron:argument:0.8@ae6cfe83` | argument | 0.8 | np.array([m - u * size + n * size * 0.8, m + u * size * 0.6, m - u * size - n * size * 0.8]), | unknown | - | - | - |
| `render:_chevron:argument:0.6@ae6cfe83` | argument | 0.6 | np.array([m - u * size + n * size * 0.8, m + u * size * 0.6, m - u * size - n * size * 0.8]), | unknown | - | - | - |
| `render:_chevron:argument:0.8@ae6cfe83#2` | argument | 0.8 | np.array([m - u * size + n * size * 0.8, m + u * size * 0.6, m - u * size - n * size * 0.8]), | unknown | - | - | - |
| `render:_flatten:default:per_curve` | default | 12 | def _flatten(path: Path, per_curve: int = 12) -> np.ndarray: | unknown | - | - | - |
| `render:grade_colour:argument:9@b686040c` | argument | 9 | return CONSURF[int(np.clip(grade, 1, 9)) - 1] if grade else NO_DATA | unknown | - | conservation | - |
| `render:_by_data.element:argument:9@a4b4b17e` | argument | 9 | return scale.colour(round(math.fsum(mine) / len(mine), 9)) if mine else DATA_NONE | convention | rounding digits | 310, app, app_controls, app_jobs +28 | - |
| `render:_draw_data_marks:argument:0.6@82c41a8d` | argument | 0.6 | ring = Circle(m["at"], _DATA_RING, fc=background, ec=legible("#1d1d1f"), lw=lw * 0.6, zorder=Z["data-ring"]) | unknown | - | - | - |
| `render:_draw_lines:argument:1.9@83c79263` | argument | 1.9 | c.ax.plot(*xy, color=_BACKGROUND.get(), lw=lw * 1.9, solid_capstyle="butt", zorder=Z["xl-over"] - 0.01) | unknown | - | - | - |
| `render:_draw_lines:argument:0.01@83c79263` | argument | 0.01 | c.ax.plot(*xy, color=_BACKGROUND.get(), lw=lw * 1.9, solid_capstyle="butt", zorder=Z["xl-over"] - 0.01) | convention | drawing order | - | - |
| `render:_draw_interactions:argument:0.85@521d3394` | argument | 0.85 | w: text_halos(c.ax, c.pt_per_unit, link_keep(c.loop_lw * 0.85 * w, c.pt_per_unit, c.fig.dpi), turned=True) | unknown | - | - | - |
| `render:_empty_figure:argument:4.0@811e688a` | argument | 4.0 | fig = Figure(figsize=(4.0, 1.5)) | unknown | - | - | - |
| `render:_empty_figure:argument:1.5@811e688a` | argument | 1.5 | fig = Figure(figsize=(4.0, 1.5)) | unknown | - | - | - |
| `render:_empty_figure:argument:9@d4c93ec6` | argument | 9 | ax.text(0.5, 0.5, "no helices or strands found", ha="center", va="center", fontsize=9, transform=ax.transAxes) | unknown | - | - | - |
| `render:_canvas:comparison:0.9@e4c0f102` | comparison | 0.9 | if x0 + entry - 0.9 > x1: # never narrower than one legend entry; a wider page may shrink, so measure again | unknown | - | - | - |
| `render:_canvas:argument:0.9@b6a10e93` | argument | 0.9 | entry = legend_entry(x0 + entry - 0.9) | unknown | - | - | - |
| `render:_canvas:argument:0.9@b03d82a5` | argument | 0.9 | x1 = max(x1, x0 + entry - 0.9) | unknown | - | - | - |
| `render:_canvas:argument:0.9@e2de3741` | argument | 0.9 | per_row = max(1, int((x1 - x0 + 0.9) // entry)) | unknown | - | - | - |
| `render:_canvas:argument:2.4@536d099c` | argument | 2.4 | bar_h = max(_DOMAIN_BAR_H, 2.4 * em + _BAR_TICK) if bars else 0.0 | unknown | - | - | - |
| `render:_canvas:fallback:page_text_small` | fallback | page_text_small | fallback.note("page_text_small") | unknown | cites random-100 C4 (layout.page_fit) | - | notes: page_scaled |
| `render:page_legend:argument:0.9@4eaaade1` | argument | 0.9 | per_row = max(1, int((width - 2 + 0.9) // entry)) | unknown | - | - | - |
| `render:_foot:argument:1e-12@a107771a` | argument | 1e-12 | t = float(np.clip(np.dot(c - p, d) / max(float(np.dot(d, d)), 1e-12), 0.0, 1.0)) | convention | a floating-point tolerance | - | - |
| `render:_loop_boxes:argument:4@f9fa9685` | argument | 4 | boxes = np.array([p.rect for p in [*layout.placed.values(), *layout.ghosts]]).reshape(-1, 4) # none: DNA only | convention | an array's shape; line last changed by 4d11ab3 Counts keep off the loops as drawn, corners rounded or smoothed, not their routes (C2b-1 fix round 1, 6FY9) (git blame; evidence only) | route | - |
| `render:_beside_strips:local:far` | local | 1e9 | far = 1e9 | unknown | - | - | - |
| `render:_beside_strips:argument:4@287e48ca` | argument | 4 | return np.array(boxes, float).reshape(-1, 4) | convention | an array's shape | - | - |
| `render:_draw_loops:argument:3.0@8a0a1a07` | argument | 3.0 | lw=loop_lw * 3.0, | unknown | - | text_clear | - |
| `render:_draw_loops:argument:0.22@404498e2` | argument | 0.22 | arrow = _chevron(loop.points, 0.22 * c.look.loop_width**0.5) | unknown | - | text_clear | - |
| `render:_draw_loops:argument:7e-4@2422e1ca` | argument | 7e-4 | zorder=z + 7e-4, | convention | drawing order | text_clear | - |
| `render:_draw_loops:argument:3@a4db0f3d` | argument | 3 | ls=(0, (3, 2)) if loop.dashed else "-", | unknown | - | text_clear | - |
| `render:_draw_loops:argument:5e-4@aa85d2e2` | argument | 5e-4 | zorder=z + 5e-4, | convention | drawing order | text_clear | - |
| `render:_draw_loops:argument:3@a4db0f3d#2` | argument | 3 | ls=(0, (3, 2)) if loop.dashed else "-", | unknown | - | text_clear | - |
| `render:_draw_loops:argument:6e-4@883d9352` | argument | 6e-4 | zorder=z + 6e-4, | convention | drawing order | text_clear | - |
| `render:_draw_membrane:argument:0.2@ba74cdd6` | argument | 0.2 | xs = np.arange(left + 0.2, right, 0.42) | unknown | - | - | - |
| `render:_draw_membrane:argument:0.42@ba74cdd6` | argument | 0.42 | xs = np.arange(left + 0.2, right, 0.42) | unknown | - | - | - |
| `render:_draw_membrane:argument:4.2@eec1c44c` | argument | 4.2 | s=(c.lw * 4.2) ** 2, | unknown | - | - | - |
| `render:_draw_membrane:argument:0.4@aba2ff49` | argument | 0.4 | linewidths=c.lw * 0.4, | unknown | - | - | - |
| `render:_draw_membrane:argument:0.1@c14a0ae7` | argument | 0.1 | t = c.ax.text(left + 0.1, y, text, ha="left", va=va, fontsize=c.font * 0.9, style="italic", color="#7a6a45") | unknown | - | - | - |
| `render:_draw_membrane:argument:0.9@c14a0ae7` | argument | 0.9 | t = c.ax.text(left + 0.1, y, text, ha="left", va=va, fontsize=c.font * 0.9, style="italic", color="#7a6a45") | unknown | - | - | - |
| `render:_clear_of_lines:fallback:glyph_width_guessed` | fallback | glyph_width_guessed | fallback.note("glyph_width_guessed") | unknown | half an em for a character the metrics lack | - | none |
| `render:_meets_box:default:pad` | default | 0.05 | def _meets_box(p, q, box, pad: float = 0.05) -> bool: | unknown | line last changed by 0a35b86 Task 2 fix round 1: face rows, columns in order, 5LDW crossing, side names off lines (plan 26) (git blame; evidence only) | membrane | - |
| `render:_meets_box:comparison:1e-12@65884cf9` | comparison | 1e-12 | if abs(d) < 1e-12: | convention | a floating-point tolerance | membrane | - |
| `render:_draw_domains:argument:0.86@182d89e9` | argument | 0.86 | fill = _mix(tone, c.look.background, 0.86) | unknown | - | - | - |
| `render:_draw_domains:argument:0.15@e7974bc0` | argument | 0.15 | ec=_mix(tone, c.look.background, (0.15 + 0.15 * (copy % 3)) if many else 0.55), | unknown | - | - | - |
| `render:_draw_domains:argument:0.15@e7974bc0#2` | argument | 0.15 | ec=_mix(tone, c.look.background, (0.15 + 0.15 * (copy % 3)) if many else 0.55), | unknown | - | - | - |
| `render:_draw_domains:argument:3@e7974bc0` | argument | 3 | ec=_mix(tone, c.look.background, (0.15 + 0.15 * (copy % 3)) if many else 0.55), | unknown | - | - | - |
| `render:_draw_domains:argument:0.55@e7974bc0` | argument | 0.55 | ec=_mix(tone, c.look.background, (0.15 + 0.15 * (copy % 3)) if many else 0.55), | unknown | - | - | - |
| `render:_draw_domains:argument:0.8@84dffb85` | argument | 0.8 | lw=c.lw * (0.8 if many else 0.6), | unknown | - | - | - |
| `render:_draw_domains:argument:0.6@84dffb85` | argument | 0.6 | lw=c.lw * (0.8 if many else 0.6), | unknown | - | - | - |
| `render:_draw_domains:argument:0.3@fcec8a1f` | argument | 0.3 | [lx, lx], [ly0, ly1], color=_mix(tone, c.look.background, 0.3), lw=c.lw * 0.8, zorder=Z["domain-leader"] | unknown | - | - | - |
| `render:_draw_domains:argument:0.8@fcec8a1f` | argument | 0.8 | [lx, lx], [ly0, ly1], color=_mix(tone, c.look.background, 0.3), lw=c.lw * 0.8, zorder=Z["domain-leader"] | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.08@60dee2b0` | argument | 0.08 | t = c.ax.text(x, y + em + 0.08, f"{text:g}", ha=ha, color=c.look.ink, **small) | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.08@e4acc971` | argument | 0.08 | y + em + 0.08, | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.6@09ddcfdc` | argument | 0.6 | lw=c.lw * 0.6, | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.8@2013b430` | argument | 0.8 | rect.set_hatchcolor(darken(seg.tone, 0.8)) | unknown | - | - | - |
| `render:_draw_domain_bar:comparison:0.2@bd427ec5` | comparison | 0.2 | if width + 0.2 <= xb - xa: | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.7@994e7c5a` | argument | 0.7 | fontsize=c.font * 0.7, | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.1@7fceaff5` | argument | 0.1 | xs = _spread([b[0] for b in below], [b[1] for b in below], c.page[0] + 0.1, c.page[1] - 0.1, 0.2) | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.1@7fceaff5#2` | argument | 0.1 | xs = _spread([b[0] for b in below], [b[1] for b in below], c.page[0] + 0.1, c.page[1] - 0.1, 0.2) | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.2@7fceaff5` | argument | 0.2 | xs = _spread([b[0] for b in below], [b[1] for b in below], c.page[0] + 0.1, c.page[1] - 0.1, 0.2) | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.6@eb74cece` | argument | 0.6 | (tick,) = c.ax.plot([xc, xt], [y, y - _BAR_TICK], color=darken(seg.tone), lw=c.lw * 0.6) | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.04@e09f0fc8` | argument | 0.04 | y - _BAR_TICK - 0.04, | unknown | - | - | - |
| `render:_draw_domain_bar:argument:0.7@994e7c5a#2` | argument | 0.7 | fontsize=c.font * 0.7, | unknown | - | - | - |
| `render:_draw_sheet_panels:argument:0.07@dffdd8a1` | argument | 0.07 | fill = _mix(c.look.background, c.look.ink, 0.07) if c.look.background != "#ffffff" else _PANEL | unknown | - | - | - |
| `render:_draw_sheet_panels:argument:0.78@f229afa9` | argument | 0.78 | fill = _mix(domain_paint(tone_of(layout, mine[0]), c.look)[0], c.look.background, 0.78) | unknown | - | - | - |
| `render:_draw_sheet_panels:argument:0.3@80959fcb` | argument | 0.3 | ec=_mix(tone, c.look.background, 0.3), | unknown | - | - | - |
| `render:_draw_sheet_panels:argument:0.8@1efc0a3f` | argument | 0.8 | lw=c.lw * 0.8, | unknown | - | - | - |
| `render:body_paint:argument:0.45@c88ad471` | argument | 0.45 | fill = look.background if outline else tint(color, 0.45) if look.fill == "pale" else color | unknown | - | - | - |
| `render:body_paint:argument:0.35@fe610406` | argument | 0.35 | fill = {"bold": color, "pale": tint(color, 0.35), "outline": look.background}[look.fill] | unknown | - | - | - |
| `render:body_paint:argument:0.15@88b441e8` | argument | 0.15 | "pale": (tint(color, 0.15), tint(color, 0.55)), | unknown | - | - | - |
| `render:body_paint:argument:0.55@88b441e8` | argument | 0.55 | "pale": (tint(color, 0.15), tint(color, 0.55)), | unknown | - | - | - |
| `render:body_paint:argument:0.7@78f052a2` | argument | 0.7 | "outline": (look.background, _mix(color, look.background, 0.7)), | unknown | - | - | - |
| `render:_draw_elements:argument:2.4@1f2ef171` | argument | 2.4 | keep = link_keep(c.lw * 2.4, c.pt_per_unit, c.fig.dpi) | unknown | - | - | - |
| `render:_draw_elements:argument:1e-9@9d8ab27c` | argument | 1e-9 | t = float(np.dot(np.subtract(spot, p.n_port), p.direction)) / max(p.length, 1e-9) | convention | a floating-point tolerance | - | - |
| `render:_draw_elements:argument:0.9@e0020d97` | argument | 0.9 | fontsize=c.font * 0.9, | unknown | - | - | - |
| `render:_draw_elements:argument:0.45@0796ac08` | argument | 0.45 | edge = to_rgba(darken(color, 0.45), 1.0) | unknown | - | - | - |
| `render:_draw_elements:argument:0.8@66e6ecf0` | argument | 0.8 | lw=edge_w * 0.8, | unknown | - | - | - |
| `render:_draw_elements:argument:0.4@adff9fa3` | argument | 0.4 | glint, see = ("#ffffff", 0.4) if striped else (tint(fill, 0.4), 0.85) # over stripes: a white glint | unknown | - | - | - |
| `render:_draw_elements:argument:0.85@517acde2` | argument | 0.85 | t.set_bbox({"boxstyle": "round,pad=0.12", "fc": c.look.background, "ec": "none", "alpha": 0.85}) | unknown | - | - | - |
| `render:_draw_residue_numbers:comparison:3.0@f5e807d0` | comparison | 3.0 | elif p.length >= 3.0: # room for both numbers and the strand letter between them | unknown | - | - | - |
| `render:_draw_residue_numbers:argument:0.45@b93ceb8f` | argument | 0.45 | reach = 0.35 if end == "N" else min(_HEAD_LEN, p.length * 0.45) + 0.3 # clear of the head | unknown | - | - | - |
| `render:_draw_residue_numbers:argument:0.62@fdcc23c4` | argument | 0.62 | *xy, str(num), ha="center", va="center", fontsize=c.font * 0.62, color=ink, zorder=Z["resnum"] | unknown | - | - | - |
| `render:absent_outline:argument:5@6512d4cc` | argument | 5 | t = np.radians(start + np.arange(5) * 22.5) | unknown | - | diff | - |
| `render:absent_outline:argument:22.5@6512d4cc` | argument | 22.5 | t = np.radians(start + np.arange(5) * 22.5) | unknown | - | diff | - |
| `render:_draw_ghosts:argument:0.8@ac917343` | argument | 0.8 | fc=tint(color, 0.8), | unknown | - | - | - |
| `render:_draw_ghosts:argument:0.8@1efc0a3f` | argument | 0.8 | lw=c.lw * 0.8, | unknown | - | - | - |
| `render:_draw_ghosts:argument:1.5@3b9fda03` | argument | 1.5 | ls=(0, (2, 1.5)), | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@4e7a6792` | argument | 0.6 | c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:9@4e7a6792` | argument | 9 | c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@4e7a6792` | argument | 0.3 | c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@4e7a6792#2` | argument | 0.6 | c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:9@4e7a6792#2` | argument | 9 | c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@4e7a6792#3` | argument | 0.6 | c.ax.add_patch(Rectangle((x + q * 0.6 / 9, y - 0.3), 0.6 / 9, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@410ef14b` | argument | 0.3 | c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5)) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@410ef14b` | argument | 0.6 | c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5)) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@410ef14b#2` | argument | 0.6 | c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5)) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@348c4584` | argument | 0.6 | step = Rectangle((x + q * 0.6 / steps, y - 0.3), 0.6 / steps, 0.6, fc=scale.colour(v), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@348c4584` | argument | 0.3 | step = Rectangle((x + q * 0.6 / steps, y - 0.3), 0.6 / steps, 0.6, fc=scale.colour(v), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@348c4584#2` | argument | 0.6 | step = Rectangle((x + q * 0.6 / steps, y - 0.3), 0.6 / steps, 0.6, fc=scale.colour(v), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@348c4584#3` | argument | 0.6 | step = Rectangle((x + q * 0.6 / steps, y - 0.3), 0.6 / steps, 0.6, fc=scale.colour(v), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@410ef14b#2` | argument | 0.3 | c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5)) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@410ef14b#3` | argument | 0.6 | c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5)) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@410ef14b#4` | argument | 0.6 | c.ax.add_patch(Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5)) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@622ff27c` | argument | 0.3 | ring = Circle((x + 0.3, y), _DATA_RING * 1.2, fc=c.look.background, ec=legible("#1d1d1f"), lw=c.lw * 0.6) | unknown | - | - | - |
| `render:_draw_legend:argument:1.2@622ff27c` | argument | 1.2 | ring = Circle((x + 0.3, y), _DATA_RING * 1.2, fc=c.look.background, ec=legible("#1d1d1f"), lw=c.lw * 0.6) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@622ff27c` | argument | 0.6 | ring = Circle((x + 0.3, y), _DATA_RING * 1.2, fc=c.look.background, ec=legible("#1d1d1f"), lw=c.lw * 0.6) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@cd800417` | argument | 0.3 | mark = data_symbol(colour[2], (x + 0.3, y), _DATA_MARK * 1.2) | unknown | - | - | - |
| `render:_draw_legend:argument:1.2@cd800417` | argument | 1.2 | mark = data_symbol(colour[2], (x + 0.3, y), _DATA_MARK * 1.2) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@d0cbe34d` | argument | 0.3 | (x, y - 0.3), 0.6, 0.6, fc="none", ec=legible(ABSENT_INK), lw=c.lw * ABSENT_LW, ls=ABSENT_DASH | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@d0cbe34d` | argument | 0.6 | (x, y - 0.3), 0.6, 0.6, fc="none", ec=legible(ABSENT_INK), lw=c.lw * ABSENT_LW, ls=ABSENT_DASH | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@d0cbe34d#2` | argument | 0.6 | (x, y - 0.3), 0.6, 0.6, fc="none", ec=legible(ABSENT_INK), lw=c.lw * ABSENT_LW, ls=ABSENT_DASH | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@7044cee0` | argument | 0.3 | c.ax.add_patch(Rectangle((x + q * w, y - 0.3), w, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@7044cee0` | argument | 0.6 | c.ax.add_patch(Rectangle((x + q * w, y - 0.3), w, 0.6, fc=shade, ec="none")) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@1859f282` | argument | 0.3 | box = Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@1859f282` | argument | 0.6 | box = Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@1859f282#2` | argument | 0.6 | box = Rectangle((x, y - 0.3), 0.6, 0.6, fc="none", ec="#9a9a9a", lw=c.lw * 0.5) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@b9327179` | argument | 0.6 | (mark,) = c.ax.plot([x, x + 0.6], [y, y], color=ink, lw=c.loop_lw * 0.85 * wide, ls=dash) | unknown | - | - | - |
| `render:_draw_legend:argument:0.85@b9327179` | argument | 0.85 | (mark,) = c.ax.plot([x, x + 0.6], [y, y], color=ink, lw=c.loop_lw * 0.85 * wide, ls=dash) | unknown | - | - | - |
| `render:_draw_legend:argument:6@a2c48e24` | argument | 6 | for q in range(6): | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@35371108` | argument | 0.1 | Rectangle((x + q * 0.1, y - 0.3), 0.1, 0.6, fc=property_colour(q / 5, colour[1]), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@35371108` | argument | 0.3 | Rectangle((x + q * 0.1, y - 0.3), 0.1, 0.6, fc=property_colour(q / 5, colour[1]), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@35371108#2` | argument | 0.1 | Rectangle((x + q * 0.1, y - 0.3), 0.1, 0.6, fc=property_colour(q / 5, colour[1]), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@35371108` | argument | 0.6 | Rectangle((x + q * 0.1, y - 0.3), 0.1, 0.6, fc=property_colour(q / 5, colour[1]), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:5@35371108` | argument | 5 | Rectangle((x + q * 0.1, y - 0.3), 0.1, 0.6, fc=property_colour(q / 5, colour[1]), ec="none") | unknown | - | - | - |
| `render:_draw_legend:argument:6@a2c48e24#2` | argument | 6 | for q in range(6): | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@95205bc7` | argument | 0.1 | (x + q * 0.1, y - 0.3), 0.1, 0.6, fc=sequence_colour(q / 5, c.look.sequence_map), ec="none" | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@95205bc7` | argument | 0.3 | (x + q * 0.1, y - 0.3), 0.1, 0.6, fc=sequence_colour(q / 5, c.look.sequence_map), ec="none" | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@95205bc7#2` | argument | 0.1 | (x + q * 0.1, y - 0.3), 0.1, 0.6, fc=sequence_colour(q / 5, c.look.sequence_map), ec="none" | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@95205bc7` | argument | 0.6 | (x + q * 0.1, y - 0.3), 0.1, 0.6, fc=sequence_colour(q / 5, c.look.sequence_map), ec="none" | unknown | - | - | - |
| `render:_draw_legend:argument:5@95205bc7` | argument | 5 | (x + q * 0.1, y - 0.3), 0.1, 0.6, fc=sequence_colour(q / 5, c.look.sequence_map), ec="none" | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@05edabce` | argument | 0.6 | (dash,) = c.ax.plot([x, x + 0.6], [y, y], color=colour, lw=c.loop_lw * 0.85, ls=(0, (3, 2))) | unknown | - | - | - |
| `render:_draw_legend:argument:0.85@05edabce` | argument | 0.85 | (dash,) = c.ax.plot([x, x + 0.6], [y, y], color=colour, lw=c.loop_lw * 0.85, ls=(0, (3, 2))) | unknown | - | - | - |
| `render:_draw_legend:argument:3@05edabce` | argument | 3 | (dash,) = c.ax.plot([x, x + 0.6], [y, y], color=colour, lw=c.loop_lw * 0.85, ls=(0, (3, 2))) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@2d7d2561` | argument | 0.6 | [x, x + 0.6], [y, y], color=colour, lw=c.lw * 2.4 * c.look.na_line, solid_capstyle="butt" | unknown | - | - | - |
| `render:_draw_legend:argument:2.4@2d7d2561` | argument | 2.4 | [x, x + 0.6], [y, y], color=colour, lw=c.lw * 2.4 * c.look.na_line, solid_capstyle="butt" | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@40a4eac9` | argument | 0.1 | [x + 0.1, x + 0.1, np.nan, x + 0.5, x + 0.5], | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@40a4eac9#2` | argument | 0.1 | [x + 0.1, x + 0.1, np.nan, x + 0.5, x + 0.5], | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@4c9d1ea7` | argument | 0.3 | [y - 0.3, y, np.nan, y - 0.3, y], | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@4c9d1ea7#2` | argument | 0.3 | [y - 0.3, y, np.nan, y - 0.3, y], | unknown | - | - | - |
| `render:_draw_legend:argument:1.2@8bedbfb8` | argument | 1.2 | lw=c.lw * 1.2, | unknown | - | - | - |
| `render:_draw_legend:argument:0.1@a08ba808` | argument | 0.1 | for n, tx in enumerate((x + 0.1, x + 0.5)): | unknown | - | - | - |
| `render:_draw_legend:argument:0.22@c03baab1` | argument | 0.22 | y + 0.22, | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@2f3c21a0` | argument | 0.3 | star = Polygon(_star(x + 0.3, y, 0.42), closed=True, fc=colour, ec=LIGAND_EDGE, lw=c.lw * 0.7) | unknown | - | - | - |
| `render:_draw_legend:argument:0.42@2f3c21a0` | argument | 0.42 | star = Polygon(_star(x + 0.3, y, 0.42), closed=True, fc=colour, ec=LIGAND_EDGE, lw=c.lw * 0.7) | unknown | - | - | - |
| `render:_draw_legend:argument:0.7@2f3c21a0` | argument | 0.7 | star = Polygon(_star(x + 0.3, y, 0.42), closed=True, fc=colour, ec=LIGAND_EDGE, lw=c.lw * 0.7) | unknown | - | - | - |
| `render:_draw_legend:argument:0.3@306ed5ee` | argument | 0.3 | swatch = Rectangle((x, y - 0.3), 0.6, 0.6, fc=colour, ec=darken(colour), lw=c.lw * 0.8) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@306ed5ee` | argument | 0.6 | swatch = Rectangle((x, y - 0.3), 0.6, 0.6, fc=colour, ec=darken(colour), lw=c.lw * 0.8) | unknown | - | - | - |
| `render:_draw_legend:argument:0.6@306ed5ee#2` | argument | 0.6 | swatch = Rectangle((x, y - 0.3), 0.6, 0.6, fc=colour, ec=darken(colour), lw=c.lw * 0.8) | unknown | - | - | - |
| `render:_draw_legend:argument:0.8@306ed5ee` | argument | 0.8 | swatch = Rectangle((x, y - 0.3), 0.6, 0.6, fc=colour, ec=darken(colour), lw=c.lw * 0.8) | unknown | - | - | - |
| `render:_draw_legend:argument:0.85@f6e03f7b` | argument | 0.85 | t = c.ax.text(x + 0.85, y, text, ha="left", va="center", fontsize=c.font) | unknown | - | - | - |
| `render:_draw_legend:argument:30@90ec6acb` | argument | 30 | per_line = max(30, int((xmax - xmin - 2 * c.pad) * c.pt_per_unit / (c.font * 0.85 * 0.52))) | unknown | - | - | - |
| `render:_draw_legend:argument:0.85@90ec6acb` | argument | 0.85 | per_line = max(30, int((xmax - xmin - 2 * c.pad) * c.pt_per_unit / (c.font * 0.85 * 0.52))) | unknown | - | - | - |
| `render:_draw_legend:argument:0.52@90ec6acb` | argument | 0.52 | per_line = max(30, int((xmax - xmin - 2 * c.pad) * c.pt_per_unit / (c.font * 0.85 * 0.52))) | unknown | - | - | - |
| `render:_draw_legend:argument:0.55@5f772696` | argument | 0.55 | ymin + 0.55, | unknown | - | - | - |
| `render:_draw_legend:argument:0.85@c515e8c2` | argument | 0.85 | fontsize=c.font * 0.85, | unknown | - | - | - |
| `render:_draw_title:argument:0.8@10fa0669` | argument | 0.8 | c.ax.text((xmin + xmax) / 2, ymax - 0.8, title, ha="center", va="center", fontsize=c.font * 1.3, fontweight="bold") | unknown | - | - | - |
| `render:_draw_title:argument:1.3@10fa0669` | argument | 1.3 | c.ax.text((xmin + xmax) / 2, ymax - 0.8, title, ha="center", va="center", fontsize=c.font * 1.3, fontweight="bold") | unknown | - | - | - |
| `render:_draw:argument:2.4@1f2ef171` | argument | 2.4 | keep = link_keep(c.lw * 2.4, c.pt_per_unit, c.fig.dpi) | unknown | - | robustness, text_clear | - |
| `render:_draw:fallback:counts_hidden` | fallback | counts_hidden | fallback.note("counts_hidden", len(hidden)) # what the notes list (provenance), counted per figure | unknown | "counts with no clear spot, left out (provenance: counts_hidden)" | text_clear | notes: counts_hidden |
| `render:_draw:fallback:mark_least_bad_spot` | fallback | mark_least_bad_spot | fallback.note("mark_least_bad_spot", len(c.fig.foldmap_crowded_marks)) | named-entry:7QUP | render._place_mark: "nothing clear within _SPOT_RADIUS (7QUP's crowded page dropped GTP:9E500): further out, else the least bad" | - | notes: crowded_marks (a warning); random-100 reports crowded_marks (not gated) |
| `render:_draw:fallback:interactions_dropped` | fallback | interactions_dropped | fallback.note("interactions_dropped", len(hb.dropped)) | unknown | interactions.verify_interactions, the independent re-check | hbonds | notes: interactions_dropped (a warning) |
| `render:_draw:fallback:interactions_hidden` | fallback | interactions_hidden | fallback.note("interactions_hidden", len(c.fig.foldmap_interactions_hidden)) | unknown | "lines wholly under text" | - | notes: interactions_hidden |
| `render:_draw:fallback:crosslinks_unplaced` | fallback | crosslinks_unplaced | fallback.note("crosslinks_unplaced", len(hb.xl_unplaced)) | unknown | an end not drawn, or both ends at one point | render | notes: crosslinks_unplaced |
| `render:save:argument:42@fa5aa9ca` | argument | 42 | with rc_context({"svg.fonttype": "none", "pdf.fonttype": 42, "svg.hashsalt": "foldmap", **_RC}): | unknown | - | - | - |
| `render:_loop_path:fallback-marker` | fallback-marker | 2 | """A loop's line as drawn: a fallback curve as it is, else the route with smooth (Style.loops curved) or rounded | unknown | not a fallback: draws a loop_curve as it is | route | - |

## route

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `route:RES` | constant | 0.25 | RES = 0.25 # grid cell size, page units | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `route:_INFLATE` | constant | 0.15 | _INFLATE = 0.15 # keep loops this far from arrows | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_layout, route | - |
| `route:_TOUCH` | constant | 0.02 | _TOUCH = 0.02 | geometry | "_INFLATE - RES/2 = 0.025: 0.02 leaves a margin"; earliest in this file: 727a379 2026-10-10 fix: every drawn loop is validated against the elements; a crossing route is red… | - | - |
| `route:_TURN` | constant | 2.5 | _TURN = 2.5 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `route:_REUSE` | constant | 6.0 | _REUSE = 6.0 # extra cost for a cell an earlier loop already uses | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `route:_SLAB_COST` | constant | 3.0 | _SLAB_COST = 3.0 # extra cost for a cell inside the membrane slab: loops between TM helices go round it | unknown | earliest in this file: bd3e716 2026-10-11 Loops pay to cross the membrane slab: connections between TM helices stay on the… | - | - |
| `route:_SLAB_RUN` | constant | 12.0 | _SLAB_RUN = 12.0 # ... a step along the slab, not across it: a loop leaves the slab at a face (random-100 C26) | named-entry:7SMR | git: 2a9da2a's message (the value's last change): "_SLAB_RUN 12, not 24 (7SMR) ... At 12 that loop runs along the face line" | route | - |
| `route:_PAD` | constant | 3.0 | _PAD = 3.0 | unknown | cites plan 27 Task 5; final review M6; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `route:_MAX_NODES` | constant | 200_000 | _MAX_NODES = 200_000 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | heldout_budgets, route | - |
| `route:_MAX_NODES_WIDE` | constant | 1_600_000 | _MAX_NODES_WIDE = 1_600_000 | measured | "Measured on the batch's pages: 7JU4's, 7QUP's and 8QEM's stuck loops ... need at most 800 000; 8VWQ's and 7Y7C's longest nucleic-acid runs 1 600 000"; earliest in this file: 79adc49 2026-10-16 A search out of A* budget runs again wider when free cells join its ends (C2b-2,… | route | - |
| `route:_MAX_NODES_PAGE` | constant | 8_000_000 | _MAX_NODES_PAGE = 8_000_000 | named-entry:7QUP | "…t runs out of _MAX_NODES falls back as before (fix round 1: 7QUP's focus view took 158 s of its 180). About 25 s at most on a large page"; earliest in this file: 1253bbb 2026-10-16 The second A* searches of one figure share a node budget (C2b-2 fix round 1) | heldout_budgets, route | - |
| `route:_MAX_NODES_WINDOWS` | constant | 8_000_000 | _MAX_NODES_WINDOWS = 8_000_000 | named-entry:7NZ2,9FZL | "…ent its budget (plan 15 final review: drawing them from it, 7NZ2's domains view went from 0 to 3 fallback curves, 9FZL's had 5), and a fresh budget …"; earliest in this file: e5b4bae 2026-10-18 Plan 15 final fix round: a coarse page's fine windows share a pool of their own … | heldout_budgets, sensitivity | - |
| `route:_MAX_CELLS` | constant | 1_500_000 | _MAX_CELLS = 1_500_000 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | route | - |
| `route:_DIRS` | constant | [(1, 0), (0, 1), (-1, 0), (0, -1)] | _DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)] | convention | the grid's four directions; earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_layout | - |
| `route:_PUSH` | constant | [0.4 + 0.25 * k for k in range(12)] | _PUSH = [0.4 + 0.25 * k for k in range(12)] # how far out along the exit the path may start | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | na_layout, route | - |
| `route:_PUSH_NEAR` | constant | [0.15 + 0.25 * k for k in range(13)] | _PUSH_NEAR = [0.15 + 0.25 * k for k in range(13)] # ... out of a port facing a stack gap (Grid.near_ports): nearer | unknown | earliest in this file: 2f2779d 2026-10-15 C34 (random-100, batch C2a fix round 1): the loops out of a shared column never … | - | - |
| `route:_PUSH_OVER` | constant | [0.4 + 0.25 * k for k in range(60)] | _PUSH_OVER = [0.4 + 0.25 * k for k in range(60)] # from a helix drawn over a sheet: out past the sheet's edge | unknown | earliest in this file: 85ed1e6 2026-10-11 Helices packing on a sheet face are drawn over the sheet plate, semi-transparent… | - | - |
| `route:_NAME_INK_UP` | constant | 0.46 | _NAME_INK_UP, _NAME_INK_DOWN = 0.46, 0.30 # em: a helix name's ink above and below its centre (Liberation Sans) | convention | font ink extents (Liberation Sans); earliest in this file: 803038c 2026-10-14 C6 (random-100): a helix name a loop's line runs through or too near moves to a … | - | - |
| `route:_NAME_INK_DOWN` | constant | 0.30 | _NAME_INK_UP, _NAME_INK_DOWN = 0.46, 0.30 # em: a helix name's ink above and below its centre (Liberation Sans) | convention | font ink extents (Liberation Sans); earliest in this file: 803038c 2026-10-14 C6 (random-100): a helix name a loop's line runs through or too near moves to a … | - | - |
| `route:_NAME_ALONG` | constant | (0.0, -0.25, 0.25, -0.4, 0.4) | _NAME_ALONG = (0.0, -0.25, 0.25, -0.4, 0.4) # where along its helix a name that must move may sit: helix lengths | named-entry:6MHS,7SMR,8UUP,1KX5 | git: 803038c's message: "puts a name that a line strikes ... on the spot where fewest lines strike it: the other side, or along the helix (_NAME_ALONG)" (6MHS, 7SMR, 8UUP, 1KX5) | sensitivity | - |
| `route:_END_TURNS` | constant | (0.0, 45.0, -45.0, 90.0, -90.0) | _END_TURNS = (0.0, 45.0, -45.0, 90.0, -90.0) # degrees a chain end's letter may turn round its stub's tip ... | unknown | introduced by ce76192 (C2b-6: 9FWZ), whose message gives no value; earliest in this file: ce76192 2026-10-17 A chain end's N or C letter that a line runs through turns round its stub's tip,… | - | - |
| `route:_END_OUT` | constant | (0.0, 0.3, 0.6) | _END_OUT = (0.0, 0.3, 0.6) # ... and page units further out than END_LABEL, nearest first | unknown | introduced by ce76192 (C2b-6: 9FWZ), whose message gives no value; earliest in this file: ce76192 2026-10-17 A chain end's N or C letter that a line runs through turns round its stub's tip,… | - | - |
| `route:_FALLBACK_OUT` | constant | 1.0 | _FALLBACK_OUT = 1.0 # page units past the page's edge a walled loop's fallback detours by | unknown | earliest in this file: bab696f 2026-10-25 The membrane slab a wall for each loop, not a cost (plan 27 T5, P33-P35) | - | - |
| `route:_FINE_PADS` | constant | (3.0, 8.0) | _FINE_PADS = (3.0, 8.0) # page units round a loop's two ports: the windows _fine_route searches at RES, nearest first | unknown | earliest in this file: 99d40e4 2026-10-17 Held-out A3/A5 (plan 15 task 1): route fallbacks and overlapping boxes, by mecha… | route | - |
| `route:_NA_TEXT_TOLL` | constant | 25.0 | _NA_TEXT_TOLL = 25.0 # the extra cost of a cell of that text, while such a loop is routed again | named-entry:6SO9,9CJI,7XSP | git: 80a1e65's message: "routed again with the text tolled on the grid (_NA_TEXT_TOLL per cell) ... 6SO9 4 hits -> 0; 9CJI and 7XSP default views clear" | - | - |
| `route:_GAP_TOLL` | constant | 50.0 | _GAP_TOLL = 50.0 # the extra cost of a cell a loop must not take in a stack gap (_apart_in_stack_gaps) | unknown | earliest in this file: 2f2779d 2026-10-15 C34 (random-100, batch C2a fix round 1): the loops out of a shared column never … | - | - |
| `route:_REROUTE_TURNS` | constant | 4 | _REROUTE_TURNS = 4 # corners a route off NA text may add to the one it replaces when it only grazed the clearance | unknown | earliest in this file: d614150 2026-10-15 C14/C11 (random-100, fix round 1): text clearances are figcheck's pixels at the … | route | - |
| `route:_REROUTE_TURNS_INK` | constant | 12 | _REROUTE_TURNS_INK = 12 # ... and when it ran through the text's ink (7XSP's B:355-374 loop: 6 more to free 2 letters) | named-entry:7XSP | "... and when it ran through the text's ink (7XSP's B:355-374 loop: 6 more to free 2 letters)"; earliest in this file: 4f821ec 2026-10-15 C11/C14 (random-100, fix round 1): with pixel clearances, a reroute through lett… | route | - |
| `route:_element_obstacles:argument:4@80acc229` | argument | 4 | return np.array([r for _, r in own], float).reshape(-1, 4).T.copy(), [k for k, _ in own] | convention | an array's shape | heldout_budgets, na_layout, na_metrics, na_text_keepoff +1 | - |
| `route:_Grid.__init__:comparison:0.1@bc4eff2f` | comparison | 0.1 | self.slab[:, (self.ys > lo + 0.1) & (self.ys < hi - 0.1)] = True | unknown | - | fidelity_fixtures, links, membrane, random100 +1 | - |
| `route:_Grid.__init__:comparison:0.1@bc4eff2f#2` | comparison | 0.1 | self.slab[:, (self.ys > lo + 0.1) & (self.ys < hi - 0.1)] = True | unknown | - | fidelity_fixtures, links, membrane, random100 +1 | - |
| `route:_Grid.cells_along:argument:4@7129d827` | argument | 4 | n = max(2, int(math.dist(p, q) / (self.res / 4))) | unknown | - | route | - |
| `route:_Grid.one_side:comparison:0.1@8f6b0163` | comparison | 0.1 | (pa[1] <= lo + 0.1 and pb[1] <= lo + 0.1) or (pa[1] >= hi - 0.1 and pb[1] >= hi - 0.1) | unknown | - | - | - |
| `route:_Grid.one_side:comparison:0.1@8f6b0163#2` | comparison | 0.1 | (pa[1] <= lo + 0.1 and pb[1] <= lo + 0.1) or (pa[1] >= hi - 0.1 and pb[1] >= hi - 0.1) | unknown | - | - | - |
| `route:_Grid.one_side:comparison:0.1@8f6b0163#3` | comparison | 0.1 | (pa[1] <= lo + 0.1 and pb[1] <= lo + 0.1) or (pa[1] >= hi - 0.1 and pb[1] >= hi - 0.1) | unknown | - | - | - |
| `route:_Grid.one_side:comparison:0.1@8f6b0163#4` | comparison | 0.1 | (pa[1] <= lo + 0.1 and pb[1] <= lo + 0.1) or (pa[1] >= hi - 0.1 and pb[1] >= hi - 0.1) | unknown | - | - | - |
| `route:_Grid.free_cell_along:argument:6@7103c3a5` | argument | 6 | near = (round(float(port[0]), 6), round(float(port[1]), 6)) in self.near_ports | convention | rounding digits | na_layout, route | - |
| `route:_Grid.free_cell_along:argument:6@7103c3a5#2` | argument | 6 | near = (round(float(port[0]), 6), round(float(port[1]), 6)) in self.near_ports | convention | rounding digits | na_layout, route | - |
| `route:_Grid.free_cell_along:argument:5@7b92a556` | argument | 5 | steps = _INFLATE + self.res / 5 * np.arange(int((_PUSH[-1] - _INFLATE) / (self.res / 5)) + 1) | unknown | - | na_layout, route | - |
| `route:_Wall.depth_rows:argument:1e-9@72bc1277` | argument | 1e-9 | scale = 2 * self.half / np.maximum(b1 - b0, 1e-9) | convention | a floating-point tolerance | - | - |
| `route:_Wall.room:comparison:1e-9@7736f6ec` | comparison | 1e-9 | if abs(out[0]) > 1e-9 or out[1] * self.side * self.up >= 0 or not self.reaches(port): | convention | a floating-point tolerance | - | - |
| `route:_walls:comparison:1e-6@baef1a44` | comparison | 1e-6 | stub = s.id in band and abs(out_vec[0]) < 1e-6 and out_vec[1] * side * up > 0 | convention | a floating-point tolerance | membrane, route | - |
| `route:clip_lengths:argument:4@a9e132e9` | argument | 4 | seg = np.asarray(segments, float).reshape(-1, 4) | convention | an array's shape | na_layout, route | - |
| `route:clip_lengths:argument:4@69a6a215` | argument | 4 | box = np.asarray(boxes, float).reshape(-1, 4) | convention | an array's shape | na_layout, route | - |
| `route:boxes_meeting:argument:4@c905b135` | argument | 4 | b = np.asarray(boxes, float).reshape(-1, 4) | convention | an array's shape | route | - |
| `route:_crosses_body:comparison:1e-9@e7148cf1` | comparison | 1e-9 | return clip_len((u0, v0), (u1, v1), (-hl, -hw, hl, hw)) > 1e-9 | convention | a floating-point tolerance | perf, route | - |
| `route:_clear:comparison:1e-9@76cae61c` | comparison | 1e-9 | and clip_len(p, q, box) > 1e-9 | convention | a floating-point tolerance | heldout_copies, na_layout, route | - |
| `route:_astar:argument:4@63718480` | argument | 4 | for d in range(4): # leave in any heading but back into the element; a sideways start pays its turn | unknown | - | route | - |
| `route:_astar:comparison:4@d6772617` | comparison | 4 | if d != (start_dir + 2) % 4: | unknown | - | route | - |
| `route:_astar:comparison:4@c542471c` | comparison | 4 | if (x, y) == goal and d != (goal_dir + 2) % 4: | unknown | - | route | - |
| `route:_astar:argument:4@354e2d77` | argument | 4 | for nd in range(4): | unknown | - | route | - |
| `route:_astar:comparison:4@dbeae7e2` | comparison | 4 | if nd == (d + 2) % 4: | unknown | - | route | - |
| `route:_snap_to_port:comparison:1e-6@6e0f3ce6` | comparison | 1e-6 | if min(abs(exit_vec[0]), abs(exit_vec[1])) > 1e-6: | convention | a floating-point tolerance | route | - |
| `route:_snap_to_port:comparison:1e-6@a8d24538` | comparison | 1e-6 | perp = 0 if abs(exit_vec[0]) < 1e-6 else 1 | convention | a floating-point tolerance | route | - |
| `route:_snap_to_port:comparison:1e-9@988f8900` | comparison | 1e-9 | if abs(pts[i][perp] - ref) > 1e-9: | convention | a floating-point tolerance | route | - |
| `route:_simplify:comparison:1e-9@36c9c86e` | comparison | 1e-9 | if abs(cross) > 1e-9: | convention | a floating-point tolerance | - | - |
| `route:_curve:default:n` | default | 24 | def _curve(pa, ea, pb, eb, n=24): | unknown | - | route | - |
| `route:_route_ports:comparison:1e-6@12f5d19b` | comparison | 1e-6 | out = lambda e: abs(e[0]) < 1e-6 and e[1] * wall.side * wall.up > 0 # noqa: E731 - a port facing its side | convention | a floating-point tolerance | route | - |
| `route:_pull_turns.pull:comparison:4@00d78cff` | comparison | 4 | if len(p) < 4 or abs(out[0]) > 1e-9: | unknown | - | heldout_budgets, heldout_copies, layout, ligands +3 | - |
| `route:_pull_turns.pull:comparison:1e-9@00d78cff` | comparison | 1e-9 | if len(p) < 4 or abs(out[0]) > 1e-9: | convention | a floating-point tolerance | heldout_budgets, heldout_copies, layout, ligands +3 | - |
| `route:_pull_turns.pull:comparison:1e-9@607e22d0` | comparison | 1e-9 | if abs(x1 - x0) > 1e-9 or abs(y2 - y1) > 1e-9 or abs(x3 - x2) > 1e-9 or abs(y1 - y0) <= r + 1e-9: | convention | a floating-point tolerance | heldout_budgets, heldout_copies, layout, ligands +3 | - |
| `route:_pull_turns.pull:comparison:1e-9@607e22d0#2` | comparison | 1e-9 | if abs(x1 - x0) > 1e-9 or abs(y2 - y1) > 1e-9 or abs(x3 - x2) > 1e-9 or abs(y1 - y0) <= r + 1e-9: | convention | a floating-point tolerance | heldout_budgets, heldout_copies, layout, ligands +3 | - |
| `route:_pull_turns.pull:comparison:1e-9@607e22d0#3` | comparison | 1e-9 | if abs(x1 - x0) > 1e-9 or abs(y2 - y1) > 1e-9 or abs(x3 - x2) > 1e-9 or abs(y1 - y0) <= r + 1e-9: | convention | a floating-point tolerance | heldout_budgets, heldout_copies, layout, ligands +3 | - |
| `route:_pull_turns.pull:comparison:1e-9@607e22d0#4` | comparison | 1e-9 | if abs(x1 - x0) > 1e-9 or abs(y2 - y1) > 1e-9 or abs(x3 - x2) > 1e-9 or abs(y1 - y0) <= r + 1e-9: | convention | a floating-point tolerance | heldout_budgets, heldout_copies, layout, ligands +3 | - |
| `route:_pull_turns:comparison:6@949e987b` | comparison | 6 | if at_a and at_b and len(out) < 6: # both pulls move the same run: the clearance that keeps both ports clear | unknown | - | route | - |
| `route:_pull_turns.keeps:comparison:1e-9@bf6d1401` | comparison | 1e-9 | return all((y - port[1]) * e[1] >= r - 1e-9 for port, e, r in ((pa, ea, reach[0]), (pb, eb, reach[1]))) | convention | a floating-point tolerance | app, browser, copies, datafile +22 | - |
| `route:_Router._joined:argument:16@94b3a97e` | argument | 16 | key = blake2b(np.packbits(free).tobytes(), digest_size=16).digest() | unknown | - | route | - |
| `route:_Router.second:fallback:route_second_search_refused` | fallback | route_second_search_refused | fallback.note("route_second_search_refused", int(self.wide_nodes > max_nodes)) | named-entry:7QUP | _MAX_NODES_PAGE: "once spent, a search that runs out of _MAX_NODES falls back as before (fix round 1: 7QUP's focus view took 158 s of its 180)" | - | random-100 reports second_searches_refused (not gated) |
| `route:_Router.drawn:comparison:1e-9@e06060ef` | comparison | 1e-9 | into = lambda p, e: abs(e[0]) < 1e-9 and e[1] * wall.side * wall.up < 0 # noqa: E731 - into the wall | convention | a floating-point tolerance | 310, app, browser, cli +36 | - |
| `route:_Router.drawn:argument:1e-6@9f109e39` | argument | 1e-6 | and _clear(pulled, (own_a, own_b), *self.elements, grow=min(reach) - 1e-6, **kw) | convention | a floating-point tolerance | 310, app, browser, cli +36 | - |
| `route:_Router.drawn:argument:1e-6@a3e6346a` | argument | 1e-6 | pts, (own_a, own_b), *self.elements, grow=_INFLATE - 1e-6, under=under, exempt=exempt | convention | a floating-point tolerance | 310, app, browser, cli +36 | - |
| `route:_Router.drawn:fallback:route_unsnapped` | fallback | route_unsnapped | fallback.note("route_unsnapped") | unknown | _TOUCH: "A drawing that still breaches it is redrawn unsnapped, then becomes the flagged fallback" | - | none |
| `route:_clear_helix_names:fallback:helix_name_moved` | fallback | helix_name_moved | fallback.note("helix_name_moved", moved) | named-entry:6MHS,7SMR,8UUP,1KX5 | git: 803038c's message: "route._clear_helix_names, after routing, puts a name that a line strikes ... (6MHS B:275-283)" | inventory | - |
| `route:_clear_termini:argument:0.05@8c1a698f` | argument | 0.05 | half = (max(gw / 2 + 0.05, em * _text_width(end) / 2), max(gh / 2, em * 0.4)) | unknown | line last changed by ce76192 A chain end's N or C letter that a line runs through turns round its stub's tip, or stands further out, to th… (git blame; evidence only) | na_text_keepoff, route | - |
| `route:_clear_termini:argument:0.4@8c1a698f` | argument | 0.4 | half = (max(gw / 2 + 0.05, em * _text_width(end) / 2), max(gh / 2, em * 0.4)) | unknown | line last changed by ce76192 A chain end's N or C letter that a line runs through turns round its stub's tip, or stands further out, to th… (git blame; evidence only) | na_text_keepoff, route | - |
| `route:_walled_fallback:fallback:loop_wall_detour` | fallback | loop_wall_detour | fallback.note("loop_wall_detour") | unknown | cites ruling P35 (plan 27 Task 5) | - | A3 (still a fallback curve, gated); A18 judges the loop against the band (gated) |
| `route:_walled_fallback.leg:comparison:1e-6@fd25fb47` | comparison | 1e-6 | if abs(out[0]) < 1e-6 and out[1] * (to_y - port[1]) < 0: | convention | a floating-point tolerance | - | - |
| `route:_walled_fallback:fallback:loop_wall_unkept` | fallback | loop_wall_unkept | fallback.note("loop_wall_unkept") | unknown | cites ruling P35 (plan 27 Task 5) | - | a WARNING log; A3 and A18 (gated) |
| `route:_runs:comparison:1e-9@baed2711` | comparison | 1e-9 | if abs(cross) > 1e-9 or on < 0 or b == a: | convention | a floating-point tolerance | - | - |
| `route:route_loops:argument:4@27d70797` | argument | 4 | rects = np.hstack([elements[0], np.array([r for _, r in keep_out], float).reshape(-1, 4).T]) | convention | an array's shape | 310, cli, domains, focus +24 | - |
| `route:route_loops:comparison:1.5@b669155f` | comparison | 1.5 | if kinked and np.hypot(*gap) < 1.5 and gap @ pa.exit_c > 0: | unknown | line last changed by 63274b1 A copy's elements are grouped as its mate's where DSSP grouped them otherwise (C2b-3, A15: 9P1X, 7QUP, 6JLN) (git blame; evidence only) | 310, cli, domains, focus +24 | - |
| `route:route_loops:comparison:1e-9@c1d86acf` | comparison | 1e-9 | if pts is None and grid.res > RES + 1e-9: # a page so large its grid is coarse: a fine grid round the ports | convention | a floating-point tolerance | 310, cli, domains, focus +24 | - |
| `route:route_loops:fallback:loop_curve` | fallback | loop_curve | fallback.note("loop_curve", sum(lp.fallback for lp in out)) # the loops as drawn (A3, the notes) | unknown | - | inventory, sensitivity | A3 "no route fallbacks" (gated); notes: loops_unrouted |
| `route:route_loops:fallback:loop_over_text` | fallback | loop_over_text | fallback.note("loop_over_text", sum(lp.crosses_labels for lp in out)) # ... (A4, the notes) | unknown | - | - | A4 "no loop over a label" (gated); notes: loops_over_text |
| `route:_in_box:argument:0.05@a4d098f3` | argument | 0.05 | n = max(1, int(math.dist(p, q) / 0.05)) | unknown | - | - | - |
| `route:_in_box:comparison:1e-9@48e36be7` | comparison | 1e-9 | if box[0] - 1e-9 <= x <= box[2] + 1e-9 and box[1] - 1e-9 <= y <= box[3] + 1e-9: | convention | a floating-point tolerance | - | - |
| `route:_in_box:comparison:1e-9@48e36be7#2` | comparison | 1e-9 | if box[0] - 1e-9 <= x <= box[2] + 1e-9 and box[1] - 1e-9 <= y <= box[3] + 1e-9: | convention | a floating-point tolerance | - | - |
| `route:_in_box:comparison:1e-9@48e36be7#3` | comparison | 1e-9 | if box[0] - 1e-9 <= x <= box[2] + 1e-9 and box[1] - 1e-9 <= y <= box[3] + 1e-9: | convention | a floating-point tolerance | - | - |
| `route:_in_box:comparison:1e-9@48e36be7#4` | comparison | 1e-9 | if box[0] - 1e-9 <= x <= box[2] + 1e-9 and box[1] - 1e-9 <= y <= box[3] + 1e-9: | convention | a floating-point tolerance | - | - |
| `route:_in_box:argument:3@09dec32e` | argument | 3 | out.append((round(float(x), 3), round(float(y), 3))) | convention | rounding digits | - | - |
| `route:_in_box:argument:3@09dec32e#2` | argument | 3 | out.append((round(float(x), 3), round(float(y), 3))) | convention | rounding digits | - | - |
| `route:_meet_in:comparison:0.05@3b9e1f74` | comparison | 0.05 | return bool((np.hypot(d[..., 0], d[..., 1]) < 0.05).any()) | unknown | - | - | - |
| `route:_apart_in_stack_gaps:comparison:1e-6@323783a0` | comparison | 1e-6 | if abs(p[0] - cx) < 1e-6 and (abs(p[1] - y0) < 1e-6 or abs(p[1] - y1) < 1e-6): | convention | a floating-point tolerance | route | - |
| `route:_apart_in_stack_gaps:comparison:1e-6@323783a0#2` | comparison | 1e-6 | if abs(p[0] - cx) < 1e-6 and (abs(p[1] - y0) < 1e-6 or abs(p[1] - y1) < 1e-6): | convention | a floating-point tolerance | route | - |
| `route:_apart_in_stack_gaps:comparison:1e-6@323783a0#3` | comparison | 1e-6 | if abs(p[0] - cx) < 1e-6 and (abs(p[1] - y0) < 1e-6 or abs(p[1] - y1) < 1e-6): | convention | a floating-point tolerance | route | - |
| `route:_apart_in_stack_gaps:comparison:1e-6@99b62185` | comparison | 1e-6 | ends["lower" if abs(p[1] - y0) < 1e-6 else "upper"] = k | convention | a floating-point tolerance | route | - |
| `route:_apart_in_stack_gaps:argument:6@78072191` | argument | 6 | grid.near_ports = {(round(cx, 6), round(y0, 6)), (round(cx, 6), round(y1, 6))} | convention | rounding digits | route | - |
| `route:_apart_in_stack_gaps:argument:6@78072191#2` | argument | 6 | grid.near_ports = {(round(cx, 6), round(y0, 6)), (round(cx, 6), round(y1, 6))} | convention | rounding digits | route | - |
| `route:_apart_in_stack_gaps:argument:6@78072191#3` | argument | 6 | grid.near_ports = {(round(cx, 6), round(y0, 6)), (round(cx, 6), round(y1, 6))} | convention | rounding digits | route | - |
| `route:_apart_in_stack_gaps:argument:6@78072191#4` | argument | 6 | grid.near_ports = {(round(cx, 6), round(y0, 6)), (round(cx, 6), round(y1, 6))} | convention | rounding digits | route | - |
| `route:_apart_in_stack_gaps:argument:1e-9@ed1421d2` | argument | 1e-9 | grid.toll = np.where(near, _GAP_TOLL * np.clip((rows - y0) / max(y1 - y0, 1e-9), 0.0, None), 0.0) | convention | a floating-point tolerance | route | - |
| `route:_apart_in_stack_gaps:argument:1e-9@38bdcc81` | argument | 1e-9 | grid.toll = np.where(near, _GAP_TOLL * np.clip((y1 - rows) / max(y1 - y0, 1e-9), 0.0, None), 0.0) | convention | a floating-point tolerance | route | - |
| `route:<module>:fallback-marker` | fallback-marker | 5 | # 0.02 leaves a margin. A drawing that still breaches it is redrawn unsnapped, then becomes the flagged fallback. | named-entry:7QUP | the comments of route_unsnapped (_TOUCH), route_second_search_refused (_MAX_NODES_PAGE) and _MAX_NODES_WINDOWS | - | - |
| `route:_clear_helix_names:fallback-marker` | fallback-marker | 2 | from the drawn text, label_boxes and layout.glyph, at the 6 pt floor); this is the fallback for lines no route | unknown | holds helix_name_moved | route | - |
| `route:_turned_over:fallback-marker` | fallback-marker | 1 | lp.crosses_labels = not lp.fallback and bool(len(segs)) and bool(lines_meeting(boxes, segs, band).any()) | unknown | not a fallback: a fallback loop is never flagged over labels | - | - |
| `route:_walled_fallback:fallback-marker` | fallback-marker | 6 | def _walled_fallback(layout: Layout, wall: _Wall, pa: Placed, pb: Placed, curve) -> list: | unknown | holds loop_wall_detour and loop_wall_unkept (ruling P35) | route | - |
| `route:_walled_fallback.edge:fallback-marker` | fallback-marker | 1 | return y1 + _FALLBACK_OUT if y >= sum(nearest_band(wall.bands, y)) / 2 else y0 - _FALLBACK_OUT | unknown | part of loop_wall_detour | - | - |
| `route:route_loops:fallback-marker` | fallback-marker | 5 | routed: dict[tuple[str, str], tuple[list, bool, bool]] = {} # points, fallback curve, over labels | unknown | holds loop_curve and loop_over_text | 310, cli, domains, focus +24 | - |
| `route:_apart_in_stack_gaps:fallback-marker` | fallback-marker | 1 | if len(ends) < 2 or loops[ends["lower"]].fallback or loops[ends["upper"]].fallback: | unknown | not a fallback: leaves fallback loops as they are | route | - |
| `route:_off_na_text:fallback-marker` | fallback-marker | 2 | when it crosses less of it (random-100 C11: 6SO9's RRM loops ran over the bound RNA's letters). Fallback | unknown | not a fallback: reroutes loops that are no fallback | heldout_copies, na_text_keepoff | - |

## route_na

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `route_na:_END_TURNS` | constant | (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4, np.pi) | _END_TURNS = (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4, np.pi) | unknown | earliest in this file: fbcee15 2026-10-15 C16 (random-100, batch C2a fix round 1): an end label names its own end | - | - |
| `route_na:_END_FINE` | constant | tuple( # every 5° between those, the least turn first, when none of them fits x for k in range(1, 36) if k % 9 for x in… | _END_FINE = tuple( # every 5° between those, the least turn first, when none of them fits | unknown | earliest in this file: 97e3f70 2026-10-15 Plan 12 T2: off-path strands out of the sheet track; barrels are their pairing r… | - | - |
| `route_na:_END_REACH` | constant | (1.0, 1.4, 1.8, 2.4, 3.0) | _END_REACH = (1.0, 1.4, 1.8, 2.4, 3.0) # times the label's distance from its end: the same first, then further out | unknown | earliest in this file: 97e3f70 2026-10-15 Plan 12 T2: off-path strands out of the sheet track; barrels are their pairing r… | - | - |
| `route_na:route_nucleic:argument:1e-9@aef251e6` | argument | 1e-9 | v = v / max(float(np.hypot(*v)), 1e-9) | convention | a floating-point tolerance | inventory, na_layout, na_text_keepoff | - |
| `route_na:route_nucleic:fallback:na_tail_straight` | fallback | na_tail_straight | fallback.note("na_tail_straight", fell.count("tail")) | unknown | na_ends._tail: "else it is drawn straight, a fallback" | - | A3 (gated); nacheck fallbacks; notes: loops_unrouted |
| `route_na:route_nucleic:fallback:na_gap_unclear` | fallback | na_gap_unclear | fallback.note("na_gap_unclear", fell.count("gap")) | unknown | na_ends._gap | - | A3 (gated); nacheck fallbacks; notes: loops_unrouted |
| `route_na:route_nucleic:fallback:na_hairpin_unclear` | fallback | na_hairpin_unclear | fallback.note("na_hairpin_unclear", fell.count("hairpin")) | unknown | na_runs._strand_between | - | A3 (gated); nacheck fallbacks; notes: loops_unrouted |
| `route_na:route_nucleic:fallback:na_run_unclear` | fallback | na_run_unclear | fallback.note("na_run_unclear", fell.count("run") + fell.count("joint")) | unknown | na_runs._drawing | - | A3 (gated); nacheck fallbacks; notes: loops_unrouted |
| `route_na:route_nucleic:fallback:na_trace_unclear` | fallback | na_trace_unclear | fallback.note("na_trace_unclear", fell.count("trace")) | unknown | route_na: "Every drawing is checked by _clear; one that fails is a fallback" | - | A3 (gated); nacheck fallbacks; notes: loops_unrouted |
| `route_na:route_nucleic:fallback-marker` | fallback-marker | 4 | their projection, cropped gaps dashed. Every drawing is checked by _clear; one that fails is a fallback. | unknown | holds the na_* counts | inventory, na_layout, na_text_keepoff | - |

## science

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `science:RESIDUE_STEP` | constant | 0.5 | RESIDUE_STEP = 0.5 # page units along a loop's route per residue from the nearer element end (render.residue_point) | geometry | "(render.residue_point)" (a copy of render's _RESIDUE_STEP); earliest in this file: 2d7344d 2026-10-16 Plan 12 T4: a loop residue's mark sits by residue distance from the nearer eleme… | - | - |
| `science:MARK_SLACK` | constant | 0.5 | MARK_SLACK = 0.5 # page units a mark may stand off its place and its bound (a star is drawn 0.38 across) | unknown | earliest in this file: 2d7344d 2026-10-16 Plan 12 T4: a loop residue's mark sits by residue distance from the nearer eleme… | science | - |
| `science:_MAX_CYCLES` | constant | 20_000 | _MAX_CYCLES = 20_000 # simple cycles looked at before S3 gives up on a sheet (sheet graphs are sparse) | convention | a work bound (S3 says when it is reached); earliest in this file: 2b117d8 2026-10-15 Plan 12 T1: science checks S1-S3, apart from the page checks | science | - |
| `science:_TOUCH` | constant | 1e-6 | _TOUCH = 1e-6 # page units two spans must share to overlap | convention | a float tolerance (1e-6); earliest in this file: e63f831 2026-10-15 Plan 12 T2 fix round 3 (item 2): S1 reads spans, S6 an accessory level with its … | - | - |
| `science:pairing_cycle:fallback:s3_cycles_capped` | fallback | s3_cycles_capped | fallback.note("s3_cycles_capped") | convention | _MAX_CYCLES, a work bound | - | S3 says "too many pairing cycles to count" |
| `science:along:argument:1e-12@2423a05e` | argument | 1e-12 | t = float(np.clip(np.dot(q - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0)) | convention | a floating-point tolerance | copies, criteria, hbonds, science +1 | - |
| `science:<module>:fallback-marker` | fallback-marker | 1 | _MAX_CYCLES = 20_000 # simple cycles looked at before S3 gives up on a sheet (sheet graphs are sparse) | unknown | _MAX_CYCLES: s3_cycles_capped (science.pairing_cycle) | - | - |

## scope

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `scope:na_outside_text:comparison:6@4917c73a` | comparison | 6 | shown = chains[:6] + ([f"+{len(chains) - 6} more"] if len(chains) > 6 else []) | unknown | - | - | - |

## seqplot

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `seqplot:_LABEL_H` | constant | 1.1 | _LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5 # in column widths | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `seqplot:_GLYPH_H` | constant | 1.5 | _LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5 # in column widths | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `seqplot:_NUM_H` | constant | 1.1 | _LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5 # in column widths | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `seqplot:_ROW_H` | constant | 1.45 | _LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5 # in column widths | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `seqplot:_CONS_H` | constant | 1.7 | _LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5 # in column widths | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `seqplot:_GAP_H` | constant | 1.5 | _LABEL_H, _GLYPH_H, _NUM_H, _ROW_H, _CONS_H, _GAP_H = 1.1, 1.5, 1.1, 1.45, 1.7, 1.5 # in column widths | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `seqplot:_FEAT_H` | constant | 1.25 | _FEAT_H, _SITE_H, _MAX_LANES = 1.25, 1.0, 4 | unknown | earliest in this file: 4ebce52 2026-10-09 UniProt annotation; citation and Zenodo metadata | - | - |
| `seqplot:_SITE_H` | constant | 1.0 | _FEAT_H, _SITE_H, _MAX_LANES = 1.25, 1.0, 4 | unknown | earliest in this file: 4ebce52 2026-10-09 UniProt annotation; citation and Zenodo metadata | - | - |
| `seqplot:_MAX_LANES` | constant | 4 | _FEAT_H, _SITE_H, _MAX_LANES = 1.25, 1.0, 4 | unknown | earliest in this file: 4ebce52 2026-10-09 UniProt annotation; citation and Zenodo metadata | - | - |
| `seqplot:_BODY` | constant | 1.25 | _BODY = 1.25 # a helix ribbon's diameter and an arrowhead's width, column widths (the figure: HELIX_W = HEAD) | unknown | earliest in this file: feaaac2 2026-10-20 Sequence view draws helices, strands and 3-10 helices as the figure does: render… | - | - |
| `seqplot:_RIBBON` | constant | 0.28 | _RIBBON = 0.28 # the ribbon's width along its axis, a fraction of a turn: a turn here is 3.6 columns long, longer | unknown | earliest in this file: feaaac2 2026-10-20 Sequence view draws helices, strands and 3-10 helices as the figure does: render… | - | - |
| `seqplot:_MARK` | constant | 0.24 | _MARK, _MARK_RING = 0.24, 0.36 # a data mark's symbol and ring, column widths | unknown | earliest in this file: e1c36ea 2026-10-21 Drawing the researcher's data (color_by=data, plan 13d-5): numbers as a stripe p… | - | - |
| `seqplot:_MARK_RING` | constant | 0.36 | _MARK, _MARK_RING = 0.24, 0.36 # a data mark's symbol and ring, column widths | unknown | earliest in this file: e1c36ea 2026-10-21 Drawing the researcher's data (color_by=data, plan 13d-5): numbers as a stripe p… | - | - |
| `seqplot:draw_sequence:default:columns` | default | 60 | columns: int = 60, | unknown | - | app, blueprint, seqplot, uniprot | - |
| `seqplot:draw_sequence:default:threshold` | default | 0.7 | threshold: float = 0.7, | unknown | - | app, blueprint, seqplot, uniprot | - |
| `seqplot:draw_sequence:comparison:10@4fbcc5cb` | comparison | 10 | if columns < 10: | unknown | - | app, blueprint, seqplot, uniprot | - |
| `seqplot:_draw:argument:0.4@b2d521da` | argument | 0.4 | ax.text(-margin + 0.4, -1.1, title, ha="left", va="center", fontsize=font * 1.25, fontweight="bold", color=ink) | unknown | - | - | - |
| `seqplot:_draw:argument:-1.1@b2d521da` | argument | -1.1 | ax.text(-margin + 0.4, -1.1, title, ha="left", va="center", fontsize=font * 1.25, fontweight="bold", color=ink) | unknown | - | - | - |
| `seqplot:_draw:argument:1.25@b2d521da` | argument | 1.25 | ax.text(-margin + 0.4, -1.1, title, ha="left", va="center", fontsize=font * 1.25, fontweight="bold", color=ink) | unknown | - | - | - |
| `seqplot:_draw:argument:0.3@e27a38c8` | argument | 0.3 | -margin + 0.3, | unknown | - | - | - |
| `seqplot:_draw:argument:0.85@9d6d3b31` | argument | 0.85 | fontsize=font * 0.85, | unknown | - | - | - |
| `seqplot:_draw:argument:0.8@37059541` | argument | 0.8 | (run - lo - 0.5, y_bot), c - run, y_top - y_bot, fill=False, ec=frame_c, lw=0.8, zorder=3 | unknown | - | - | - |
| `seqplot:_draw:argument:3@37059541` | argument | 3 | (run - lo - 0.5, y_bot), c - run, y_top - y_bot, fill=False, ec=frame_c, lw=0.8, zorder=3 | convention | drawing order | - | - |
| `seqplot:_draw:comparison:10@613f085a` | comparison | 10 | if p is not None and track.numbers[p] % 10 == 0: | unknown | - | - | - |
| `seqplot:_draw:argument:0.72@e2b2d1e9` | argument | 0.72 | fontsize=font * 0.72, | unknown | - | - | - |
| `seqplot:_draw:argument:0.42@5fc46e39` | argument | 0.42 | bar = Rectangle((x - 0.42, base), 0.84, max(cons[c], 0.02) * (_CONS_H - 0.4), fc=BAR, ec="none") | unknown | - | - | - |
| `seqplot:_draw:argument:0.84@5fc46e39` | argument | 0.84 | bar = Rectangle((x - 0.42, base), 0.84, max(cons[c], 0.02) * (_CONS_H - 0.4), fc=BAR, ec="none") | unknown | - | - | - |
| `seqplot:_draw:argument:0.02@5fc46e39` | argument | 0.02 | bar = Rectangle((x - 0.42, base), 0.84, max(cons[c], 0.02) * (_CONS_H - 0.4), fc=BAR, ec="none") | unknown | - | - | - |
| `seqplot:_draw:argument:0.4@5fc46e39` | argument | 0.4 | bar = Rectangle((x - 0.42, base), 0.84, max(cons[c], 0.02) * (_CONS_H - 0.4), fc=BAR, ec="none") | unknown | - | - | - |
| `seqplot:_draw:argument:0.45@d176c7cf` | argument | 0.45 | (lo - b * columns - 0.45, y), hi - lo + 0.9, 0.7, fc=colour, ec="none", alpha=0.85, zorder=2 | unknown | - | - | - |
| `seqplot:_draw:argument:0.9@d176c7cf` | argument | 0.9 | (lo - b * columns - 0.45, y), hi - lo + 0.9, 0.7, fc=colour, ec="none", alpha=0.85, zorder=2 | unknown | - | - | - |
| `seqplot:_draw:argument:0.7@d176c7cf` | argument | 0.7 | (lo - b * columns - 0.45, y), hi - lo + 0.9, 0.7, fc=colour, ec="none", alpha=0.85, zorder=2 | unknown | - | - | - |
| `seqplot:_draw:argument:0.85@d176c7cf` | argument | 0.85 | (lo - b * columns - 0.45, y), hi - lo + 0.9, 0.7, fc=colour, ec="none", alpha=0.85, zorder=2 | unknown | - | - | - |
| `seqplot:_draw:comparison:0.55@33439df6` | comparison | 0.55 | text = f.label if len(f.label) * 0.55 < room else f.label[: max(0, int(room / 0.55) - 1)] + "…" | unknown | - | - | - |
| `seqplot:_draw:argument:0.55@33439df6` | argument | 0.55 | text = f.label if len(f.label) * 0.55 < room else f.label[: max(0, int(room / 0.55) - 1)] + "…" | unknown | - | - | - |
| `seqplot:_draw:comparison:2.5@e67cdae8` | comparison | 2.5 | if room > 2.5: | unknown | - | - | - |
| `seqplot:_draw:argument:0.35@9ae3bcb4` | argument | 0.35 | y + 0.35, | unknown | - | - | - |
| `seqplot:_draw:argument:0.68@958af6bc` | argument | 0.68 | fontsize=font * 0.68, | unknown | - | - | - |
| `seqplot:_draw:argument:3@f168af61` | argument | 3 | zorder=3, | convention | drawing order | - | - |
| `seqplot:_draw:argument:0.3@504b83a8` | argument | 0.3 | dot = Circle((c - b * columns, y), 0.3, fc=SITE_COLOURS.get(f.type, "#555555"), ec="none", zorder=3) | unknown | - | - | - |
| `seqplot:_draw:argument:3@504b83a8` | argument | 3 | dot = Circle((c - b * columns, y), 0.3, fc=SITE_COLOURS.get(f.type, "#555555"), ec="none", zorder=3) | convention | drawing order | - | - |
| `seqplot:_draw:argument:0.45@0d7afae6` | argument | 0.45 | top(b) - _LABEL_H * 0.45, | unknown | - | - | - |
| `seqplot:_draw:argument:0.9@4ac75e98` | argument | 0.9 | fontsize=font * 0.9, | unknown | - | - | - |
| `seqplot:_draw:argument:0.85@608cdb95` | argument | 0.85 | color=legible(darken(colour, 0.85)) if colour != "#000000" else ink, | unknown | - | - | - |
| `seqplot:_draw:argument:0.6@aa3718ba` | argument | 0.6 | ring = Circle(at, _MARK_RING, fc=look.background, ec=legible(INK), lw=0.6 * font / 8 * 0.6, zorder=3) | unknown | - | - | - |
| `seqplot:_draw:argument:8@aa3718ba` | argument | 8 | ring = Circle(at, _MARK_RING, fc=look.background, ec=legible(INK), lw=0.6 * font / 8 * 0.6, zorder=3) | unknown | - | - | - |
| `seqplot:_draw:argument:0.6@aa3718ba#2` | argument | 0.6 | ring = Circle(at, _MARK_RING, fc=look.background, ec=legible(INK), lw=0.6 * font / 8 * 0.6, zorder=3) | unknown | - | - | - |
| `seqplot:_draw:argument:3@aa3718ba` | argument | 3 | ring = Circle(at, _MARK_RING, fc=look.background, ec=legible(INK), lw=0.6 * font / 8 * 0.6, zorder=3) | convention | drawing order | - | - |
| `seqplot:_draw:argument:3.1@2f820255` | argument | 3.1 | mark.set(facecolor=colour, edgecolor="none", zorder=3.1) | convention | drawing order | - | - |
| `seqplot:_draw:argument:0.75@b11abb18` | argument | 0.75 | c - b * columns, glyph_y(b), "TT", ha="center", va="center", fontsize=font * 0.75, color=ink | unknown | - | - | - |
| `seqplot:_body:argument:0.55@b7366036` | argument | 0.55 | shade = _mix(colour, look.background, 0.55) # the plain black style: a grey back face keeps the twist | unknown | - | - | - |
| `seqplot:_body:argument:3.6@e395b595` | argument | 3.6 | turns = max(1, round(length / 3.6)) | convention | rounding digits | - | - |
| `seqplot:_body:argument:0.8@7d439fde` | argument | 0.8 | behind = PathPatch(back, fc=shade, ec=edge, lw=edge_w * 0.8, joinstyle="round", zorder=1.9) | unknown | - | - | - |
| `seqplot:_body:argument:1.9@7d439fde` | argument | 1.9 | behind = PathPatch(back, fc=shade, ec=edge, lw=edge_w * 0.8, joinstyle="round", zorder=1.9) | convention | drawing order | - | - |
| `seqplot:_body:argument:0.4@9b01d2b3` | argument | 0.4 | glint, see = ("#ffffff", 0.4) if painted else (tint(fill, 0.4), 0.85) # over stripes: a white glint | unknown | - | - | - |
| `seqplot:_body:argument:2.05@4bb181f1` | argument | 2.05 | gloss = PathPatch(shine, fc=glint, ec="none", alpha=see, zorder=2.05) | convention | drawing order | - | - |
| `seqplot:_stripes:argument:0.01@aaf19063` | argument | 0.01 | edge = PathPatch(body.get_path(), fc="none", ec=body.get_edgecolor(), lw=body.get_linewidth(), zorder=z + 0.01) | convention | drawing order | - | - |

## sequence

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `sequence:GRADE_STEP` | constant | 0.5 | GRADE_STEP = 0.5 # standard deviations per conservation grade | unknown | earliest in this file: bf7a59d 2026-10-09 Conservation: Henikoff-weighted entropy, z-score grades, labelled as what it is | - | - |
| `sequence:MIN_IDENTITY` | constant | 0.5 | MIN_IDENTITY = 0.5 # a chain takes an alignment's scores only if it is at least this identical to the reference | unknown | earliest in this file: 7a98708 2026-10-09 Add sequence view with alignments, DSSP turns, conservation colouring | - | - |
| `sequence:read_consurf:comparison:5@037ac2fa` | comparison | 5 | if len(parts) < 5 or not parts[0].isdigit(): | unknown | - | conservation | - |
| `sequence:read_consurf:comparison:9@c67ad2f5` | comparison | 9 | if atom and grade.isdigit() and 1 <= int(grade) <= 9: | unknown | - | conservation | - |
| `sequence:read_consurf:argument:4@ee7c92b2` | argument | 4 | out[(atom.group(4), int(atom.group(2)), atom.group(3))] = (atom.group(1), int(grade)) | unknown | - | conservation | - |
| `sequence:read_consurf:argument:3@ee7c92b2` | argument | 3 | out[(atom.group(4), int(atom.group(2)), atom.group(3))] = (atom.group(1), int(grade)) | unknown | - | conservation | - |
| `sequence:column_classes:default:threshold` | default | 0.7 | def column_classes(aln: Alignment, threshold: float = 0.7) -> tuple[np.ndarray, np.ndarray]: | unknown | - | seqplot, sequence | - |
| `sequence:conservation:argument:20@faf57ff1` | argument | 20 | out[c] = (1.0 - h / np.log(20)) * present | unknown | - | app, app_jobs, conservation, mcp +1 | - |
| `sequence:grades:comparison:1e-9@b24bd655` | comparison | 1e-9 | if sd < 1e-9: | convention | a floating-point tolerance | app, app_controls, conservation, sequence | - |
| `sequence:grades:argument:5@3c23e7a8` | argument | 5 | out[ok] = np.clip(np.floor(z / GRADE_STEP + 0.5).astype(int) + 5, 1, 9) | unknown | - | app, app_controls, conservation, sequence | - |
| `sequence:grades:argument:9@3c23e7a8` | argument | 9 | out[ok] = np.clip(np.floor(z / GRADE_STEP + 0.5).astype(int) + 5, 1, 9) | unknown | - | app, app_controls, conservation, sequence | - |

## sheets

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `sheets:_SPLIT_GAP` | constant | 4 | _SPLIT_GAP = 4 # strands of one chain this close end to end are pieces of one long strand | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `sheets:_BARREL_MIN` | constant | 5 | _BARREL_MIN = 5 # strands a sheet needs before a lone end-to-end bridge may close it into a barrel | unknown | earliest in this file: 24fa8fe 2026-10-09 Fix bugs found by a 60-entry PDB sweep; add scripts/sweep.py | sheets | - |
| `sheets:_MIN_PAIRS` | constant | 2 | _MIN_PAIRS = 2 # bridges needed between two strands to count them as neighbours | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |

## source_stamp

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `source_stamp:DIRTY_LISTED` | constant | 50 | DIRTY_LISTED = 50 | convention | a listing limit; earliest in this file: 7fcf7a5 2026-10-25 Evidence integrity: source stamps, the paper pipeline in scripts/paper, a pinned… | - | - |
| `source_stamp:_git.run:argument:60@b23462db` | argument | 60 | done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, timeout=60) | unknown | - | - | - |

## spec

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `spec:parse_vector:comparison:3@7275d69b` | comparison | 3 | if len(v) != 3: | unknown | - | - | - |

## ss

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `ss:_MIN_STRAND` | constant | 2 | _MIN_STRAND = 2 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `ss:_MIN_HELIX_H` | constant | 4 | _MIN_HELIX_H = 4 # an alpha run needs this many H residues ... | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `ss:_MIN_HELIX_ANY` | constant | 6 | _MIN_HELIX_ANY = 6 # ... or this many helix-class residues in total (3-10 / pi) | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `ss:_MIN_310` | constant | 3 | _MIN_310 = 3 # an isolated 3-10 run this long becomes its own small element (kind G) when asked for | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `ss:_PARTNER_GAP` | constant | 4 | _PARTNER_GAP = 4 # residues of a partner strand this close pair with the same side of it | unknown | earliest in this file: 4f83002 2026-10-09 Strands one residue apart are one strand only when they pair with the same strand | - | - |
| `ss:KINK_SPLIT` | constant | 30.0 | KINK_SPLIT = 30.0 # degrees: alpha segments of one helix run meeting at a 3-10/pi kink sharper than this split it | unknown | earliest in this file: bdf507c 2026-10-11 Helix runs split where α segments meet at a 3₁₀/π kink of more than 30° (2HHB α … | ss | - |
| `ss:_MIN_KINK_ARM` | constant | 7 | _MIN_KINK_ARM = 7 # residues: an alpha segment shorter than this has no axis precise enough to measure a kink by | unknown | earliest in this file: d8c0fbb 2026-10-11 Kinks are measured on local helix axes between arms of 7+ residues (π-bulges sta… | - | - |
| `ss:_BREAK_MAX` | constant | 2 | _BREAK_MAX = 2 # residues: a helix broken by this few in one copy of a protein is broken there in every copy | unknown | introduced by c18a18d (7OKA), whose message ties 2 to the kink rule (one or two residues; 7OKA's break is one), not to 7OKA; its message: "breaks a helix run where DSSP breaks it by one or two residues in another copy ... 7OKA: B~A, C~A 0.000" | - | - |
| `ss:_kinks:argument:4@b5b03a7a` | argument | 4 | (x, y, round(_angle(axis(x), axis(y)), 4)) for x, y in zip(hs, hs[1:]) if set(ss[x[1] + 1 : y[0]]) <= set("GI") | convention | rounding digits | - | - |

## style

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `style:THEMES` | constant | { "publication": Style(), "minimal": Style(fill="outline", helix_shading="none", loop_width=0.8, na_style="block"), "pr… | THEMES = { | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | app_controls, blueprint, conservation, domains +6 | - |

## symmetry

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `symmetry:TOLERANCE` | constant | 3.0 | TOLERANCE = 3.0 # Å, CA RMSD between a copy and its image | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `symmetry:SAME` | constant | 0.9 | SAME = 0.9 # aligned identity for two chains to be copies | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `symmetry:MIN_RES` | constant | 10 | MIN_RES = 10 # chains shorter than this are left out | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `symmetry:MIN_TURN` | constant | np.radians(10) | MIN_TURN = np.radians(10) # smaller rotations are lattice-like translations, not symmetry | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `symmetry:MIN_RISE` | constant | 2.0 | MIN_RISE = 2.0 # Å along the axis per copy for a helical filament | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `symmetry:MIN_TWIST` | constant | np.radians(3.0) | MIN_TWIST = np.radians(3.0) # a helical filament turns at least this much per copy (less: a translation) | unknown | earliest in this file: 24e5d40 2026-10-09 Symmetry: a helical filament needs a real twist (CI found pure translations read… | - | - |
| `symmetry:MIN_HELIX_COPIES` | constant | 4 | MIN_HELIX_COPIES = 4 # copies on a screw before auto calls it helical (three: crystal packing, as in 1URN) | named-entry:1URN | "…before auto calls it helical (three: crystal packing, as in 1URN)"; earliest in this file: 19b1916 2026-10-11 Protein-only ASUs with several copies draw assembly 1; a 3-copy screw is no long… | - | - |
| `symmetry:MIN_HELIX_COPIES_ASKED` | constant | 3 | MIN_HELIX_COPIES_ASKED = 3 # ... and when the user asks for "helical" | unknown | earliest in this file: b1c2790 2026-10-11 Tests at the helical threshold and the kind gate; the explicit-helical copy coun… | - | - |
| `symmetry:_PSEUDO_21` | constant | 170.0 | _PSEUDO_21 = 170.0 # degrees: a helix turning more than this per copy is a pseudo-2₁ screw, its hand not named | unknown | earliest in this file: d184032 2026-10-09 Symmetry: a pseudo-2₁ screw reports its signed twist without naming a hand | - | - |
| `symmetry:MAX_ORDER` | constant | 24 | MAX_ORDER = 24 | unknown | earliest in this file: 0438931 2026-10-09 Rename to foldmap; add README | - | - |
| `symmetry:_NEAREST` | constant | 3 | _NEAREST = 3 # copies tried as the image of each subunit: those whose centroids land closest | unknown | earliest in this file: 25288c8 2026-10-09 Symmetry search: cached superpositions, nearest images only (icosahedral 35 s to… | - | - |
| `symmetry:_Assembly.pair:argument:3@09793763` | argument | 3 | np.array([self.atoms[a][k] for k in keys]).reshape(-1, 3), | convention | an array's shape | - | - |
| `symmetry:_Assembly.pair:argument:3@b9f1aad9` | argument | 3 | np.array([self.atoms[b][k] for k in keys]).reshape(-1, 3), | convention | an array's shape | - | - |
| `symmetry:axis_angle:comparison:1e-6@dedecfeb` | comparison | 1e-6 | if np.linalg.norm(axis) < 1e-6: # 0 or 180 degrees: the eigenvector of eigenvalue +1 | convention | a floating-point tolerance | copies | - |
| `symmetry:axis_angle:comparison:1e-6@87e463de` | comparison | 1e-6 | if axis[np.flatnonzero(np.abs(axis) > 1e-6)[0]] < 0: | convention | a floating-point tolerance | copies | - |
| `symmetry:_turn:argument:3@76664890` | argument | 3 | return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k | convention | an array's shape | align, like | - |
| `symmetry:_screw:argument:3@afce1fe1` | argument | 3 | centre = np.linalg.lstsq(np.eye(3) - r, perp, rcond=None)[0] | convention | an array's shape | - | - |
| `symmetry:_candidate_axes:comparison:0.99@77c23950` | comparison | 0.99 | if rmsd <= tol and angle >= MIN_TURN and not any(abs(ax @ x) > 0.99 for x in axes): | unknown | - | - | - |
| `symmetry:_protomers:argument:8.0@897c1375` | argument | 8.0 | first.append(max(orbit, key=lambda c: (sum(tree[c].count_neighbors(tree[x], 8.0) for x in first), c))) | unknown | - | metamorphic | - |
| `symmetry:_cubic:comparison:3@35218ade` | comparison | 3 | if n == 3 and not any(abs(ax @ t) > 0.99 for t in threes): | unknown | - | - | - |
| `symmetry:_cubic:comparison:0.99@35218ade` | comparison | 0.99 | if n == 3 and not any(abs(ax @ t) > 0.99 for t in threes): | unknown | - | - | - |
| `symmetry:_cubic:comparison:5@20ef9c0f` | comparison | 5 | return "I" if 5 in found else "O" if 4 in found else "T" | unknown | - | - | - |
| `symmetry:_cubic:comparison:4@20ef9c0f` | comparison | 4 | return "I" if 5 in found else "O" if 4 in found else "T" | unknown | - | - | - |
| `symmetry:_dihedral:comparison:0.15@50676257` | comparison | 0.15 | if abs(ax @ axis) < 0.15 and _cyclic(asm, ax, centre, 2, tol) is not None: | unknown | - | - | - |
| `symmetry:_around:comparison:1e-6@cd227575` | comparison | 1e-6 | while out and a < last - 1e-6: | convention | a floating-point tolerance | - | - |
| `symmetry:_from_user:argument:15@a7596102` | argument | 15 | if abs(rise) >= MIN_RISE and abs(angle - 2 * np.pi / n) > np.radians(15): | unknown | - | - | - |
| `symmetry:detect_symmetry:comparison:3@9c5bcb2e` | comparison | 3 | cubic = _cubic(asm, centre, tolerance) if dihedral or order >= 3 else None | unknown | - | assembly, copies, focus, symmetry +1 | - |

## text_metrics

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `text_metrics:UNITS` | constant | 2048 | UNITS = 2048 # font units per em | convention | "font units per em" (Liberation Sans); earliest in this file: af1d469 2026-10-11 Stem names are measured from built-in Liberation Sans Bold metrics, not installe… | text_clear | - |
| `text_metrics:ASCENT` | constant | 1491 | ASCENT = 1491 # OS/2 typographic ascender and descender, both weights: the line a text is centred on | convention | "OS/2 typographic ascender" (Liberation Sans); earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | text_clear | - |
| `text_metrics:DESCENT` | constant | 431 | DESCENT = 431 | convention | typographic descender (Liberation Sans); earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | text_clear | - |
| `text_metrics:ADVANCE` | constant | { # character: advance width, font units (Bold) " ": 569, "!": 682, '"': 971, "#": 1139, "$": 1139, "%": 1821, "&": 147… | ADVANCE = { # character: advance width, font units (Bold) | convention | font metrics (Liberation Sans Bold); earliest in this file: af1d469 2026-10-11 Stem names are measured from built-in Liberation Sans Bold metrics, not installe… | names | - |
| `text_metrics:KERN` | constant | { # character pair: kerning added between them, font units (Bold) "AT": -152, "AV": -152, "AW": -113, "AY": -188, "A’":… | KERN = { # character pair: kerning added between them, font units (Bold) | convention | font metrics (Liberation Sans Bold); earliest in this file: af1d469 2026-10-11 Stem names are measured from built-in Liberation Sans Bold metrics, not installe… | - | - |
| `text_metrics:ADVANCE_NORMAL` | constant | { # Regular " ": 569, "!": 569, '"': 727, "#": 1139, "$": 1139, "%": 1821, "&": 1366, "'": 391, "(": 682, ")": 682, "*"… | ADVANCE_NORMAL = { # Regular | convention | font metrics (Liberation Sans Regular); earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | - | - |
| `text_metrics:KERN_NORMAL` | constant | { # Regular "AT": -152, "AV": -152, "AW": -76, "AY": -152, "A’": -152, "A ": -113, "Av": -37, "Aw": -37, "Ay": -37, "ΑΟ… | KERN_NORMAL = { # Regular | convention | font metrics (Liberation Sans Regular); earliest in this file: d644c3b 2026-10-13 Stem names sit 0.25 em beyond their own letter row (1EHZ D and T); NA campaign r… | - | - |

## uniprot

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `uniprot:MIN_IDENTITY` | constant | 0.5 | MIN_IDENTITY = 0.5 # a chain matched by alignment must be at least this identical to the UniProt sequence | unknown | earliest in this file: 4ebce52 2026-10-09 UniProt annotation; citation and Zenodo metadata | - | - |
| `uniprot:Feature.label:comparison:40@0d67f60f` | comparison | 40 | return text if len(text) <= 40 else text[:38] + "…" | unknown | - | - | - |
| `uniprot:_get:argument:30@f5c94768` | argument | 30 | with urllib.request.urlopen(req, timeout=30) as r: | unknown | - | - | - |
| `uniprot:parse_entry:argument:99.0@a4ed7e13` | argument | 99.0 | structures.sort(key=lambda s: (-s.coverage, s.resolution if s.resolution is not None else 99.0, s.id)) | unknown | - | app, app_controls | - |
| `uniprot:describe:default:limit` | default | 12 | def describe(e: Entry, limit: int = 12) -> str: | unknown | - | uniprot | - |

## userdata

| id | kind | value | at | origin | evidence | tests naming it | check |
|---|---|---|---|---|---|---|---|
| `userdata:CROSSLINK_SPAN` | constant | 30.0 | CROSSLINK_SPAN = 30.0 | biology | "30 Å, the usual limit for the lysine-reactive NHS esters (DSS, BS3"; earliest in this file: 243697e 2026-10-21 Crosslinks drawn, the app's data upload and MCP data parameters (plan 13d-6): ea… | render | - |
| `userdata:SAID_ROWS` | constant | 20 | SAID_ROWS = 20 # rows of a file said one by one (each a Problem), the rest counted in one; also the most lines of | convention | a listing limit; earliest in this file: 112b30d 2026-10-21 Your data, fixed after the 13b-13d section review (Important 3-7, Minor 8-21) | userdata | - |
| `userdata:ECHO_MAX` | constant | 40 | ECHO_MAX = 40 # characters of a cell said back in a problem (a crafted cell is not repeated whole) | convention | a listing limit; earliest in this file: 112b30d 2026-10-21 Your data, fixed after the 13b-13d section review (Important 3-7, Minor 8-21) | userdata | - |
| `userdata:MAX_CATEGORIES` | constant | 28 | MAX_CATEGORIES = 28 # categories a column may hold: render.DATA_CATEGORY_COLOURS x DATA_SHAPES marks told apart | convention | "render.DATA_CATEGORY_COLOURS x DATA_SHAPES marks told apart"; earliest in this file: 112b30d 2026-10-21 Your data, fixed after the 13b-13d section review (Important 3-7, Minor 8-21) | userdata | - |
| `userdata:CHAINS_SAID` | constant | 12 | CHAINS_SAID = 12 # the structure's protein chains named when a row's chain is none of them (a capsid has hundreds) | convention | a listing limit; earliest in this file: 112b30d 2026-10-21 Your data, fixed after the 13b-13d section review (Important 3-7, Minor 8-21) | userdata | - |
| `userdata:MAX_COLUMNS` | constant | 64 | MAX_COLUMNS = 64 # a data or crosslink table has a few columns; more is refused before any row is kept | convention | an input limit; earliest in this file: 9b120a2 2026-10-21 Uploaded tables bounded: at most 64 columns and 100,000 rows, short rows padded … | - | - |
| `userdata:MAX_ROWS` | constant | 100_000 | MAX_ROWS = 100_000 # rows kept per table (a protein has far fewer residues; a crafted file must not fill memory) | convention | an input limit; earliest in this file: 9b120a2 2026-10-21 Uploaded tables bounded: at most 64 columns and 100,000 rows, short rows padded … | userdata | - |
| `userdata:crosslink_lengths:argument:3@10d4d233` | argument | 3 | out.append(round(d, 3)) | convention | rounding digits | - | - |
