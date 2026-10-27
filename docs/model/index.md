# foldmap: a model of the system

Plan 28 Task 9.4 (the review's "stabilize and characterize", ruling Q51). One page per subsystem. Each page says what
the subsystem assumes, what it guarantees, how it fails and which checks see the failure. The architecture section of
the paper and the developer docs are drawn from these pages.

These are developer docs. They are not in the sdist (MANIFEST.in grafts only docs/examples, docs/names and
docs/domains).

## The pages

The order is the pipeline's: a file is read, secondary structure assigned, strands grouped into sheets, symmetry found,
elements laid out, loops routed. Nucleic acids, the membrane and domains are layers on that path. The checks judge the
result.

| page | subsystem | main entry |
|---|---|---|
| [reading.md](reading.md) | reading a structure file | `io.read_model`, `io.load_backbone` |
| [dssp.md](dssp.md) | secondary structure (DSSP) and elements | `dssp.assign_dssp`, `ss.build_sses` |
| [sheets.md](sheets.md) | strands grouped into sheets | `sheets.build_sheets` |
| [symmetry.md](symmetry.md) | symmetry and copies | `symmetry.detect_symmetry` |
| [layout.md](layout.md) | elements placed on the page | `layout.build_layout` |
| [routing.md](routing.md) | loops routed between elements | `route.route_loops` |
| [nucleic-acids.md](nucleic-acids.md) | DNA and RNA | `nucleic.secondary_structure`, `route_na.route_nucleic` |
| [membrane.md](membrane.md) | the membrane band | `membrane.find_membrane` |
| [domains.md](domains.md) | automatic domains | `features.auto_domains` |
| [checks.md](checks.md) | the checks themselves | `criteria.check_figure`, `figcheck.science_checks` |

## Sections of each page

- **Scope**: the modules and the main entry functions.
- **Assumes**: what the subsystem takes as given about its input.
- **Guarantees**: only what a named test or a committed record shows. Anything else is marked "not measured".
- **Fails**: known failure modes, with named PDB entries.
- **Seen by**: the check ids that see each failure, and the failures no check sees.
- **Evidence**: links to the committed records.

## How names are written

The anti-drift test `tests/test_model_docs.py` reads every backticked name on these pages:

- **module.symbol** is a name in live code, in the package `foldmap` (`layout._open_at_ends`, `like.Like.order`). The
  test imports the module and looks the name up. Line numbers are never used: they drift.
- **module:NAME** is a record id, as the inventory, sensitivity and prune records write it. It is used where a page
  points into those records, and for a name the code no longer has (a pruned constant such as `layout:_NAME_OPEN`).
  The test looks it up in those records, not in the code.
- Check ids: `S1` to `S6` (science checks, `science.SCIENCE`), `A1` to `A20` with `A1c` and `A2b` (page criteria,
  `criteria.CRITERIA`), `R1` (the readability floor, `floor.CODE`), `D1` to `D3` (deep checks, `deep.DEEP`), `H1` to
  `H5` (a figure drawn like a reference, `criteria.LIKE_CHECKS`). nacheck's measures are named as keys of
  `nacheck.THRESHOLDS`.
- Paths (`docs/validation/prune/report.md`) and tests (`tests/test_metamorphic.py::test_structure_facts_are_exact`)
  must exist.

## The records

- `docs/validation/inventory/inventory.md`: every constant, threshold and fallback, with its origin.
- `docs/validation/sensitivity/report.md`: each tuned constant at x0.75 and x1.25.
- `docs/validation/prune/report.md`: special cases turned off one at a time; two pruned.
- `docs/validation/metamorphic/report.md`: transforms that should change nothing, and independent kernels.
- `docs/validation/na-fixes/report.md`: the nucleic-acid fixes F2-F7, measured before and after.
- `docs/validation/domain-uncertainty/report.md`: why automatic domains fail, and what the split confidence is worth.
- `docs/validation/readability-floor/report.md`: figures too large to read, and the options that draw less.
- `docs/validation/homologues/review.md` and `docs/validation/homologues/topology-transfer/report.md`: figures drawn
  like a reference.
- `docs/validation/random100/review.md`, `docs/validation/random100-heldout/review.md`,
  `docs/validation/random100-third/review.md`: the three random samples. Only the third was unseen when it ran (at
  7696947); none is unseen now. `docs/validation/fourth-sample-protocol.md` fixes the unseen fourth sample.
- `docs/validation/na/report.md`: the nucleic-acid campaign (older code).

## Freeze rule

Plan 28 Task 9.5: no new feature and no new special case without a named mechanism and a measurement on an unseen
set. A page changes when the code or a record changes; `tests/test_model_docs.py` fails when a quoted name, check,
record or test is gone, and when a constant quoted with its value (`layout._SIDE_CLEAR` (2.75 Å)) no longer has
that value in the code. Numbers quoted from records are not pinned: a regenerated record means re-reading the pages
that cite it.
