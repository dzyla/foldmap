# Contributing to foldmap

foldmap is frozen for stabilisation (plan 28, the answer to an external review that asked us to stabilise and
characterise the code, not to rewrite it). This page says what a change needs. The model of the system is in
[docs/model/index.md](docs/model/index.md); read the page of the subsystem you touch first.

## The freeze rule

No new feature and no new special case without a named mechanism and a measurement on an unseen set. This applies to
fixes as much as to features.

**A named mechanism** is the code path and the reason the structure triggers it: which function decides, on which
property of the structure, and why that property makes the drawing wrong. "Entry 1ABC looked bad" is not a mechanism.
A new constant is a special case unless its value follows from biology, geometry or a convention you can cite.

**An unseen set** is a sample you did not use to tune the change. The samples that exist are in `docs/validation/`,
and none of them is unseen now:

- `random100/`: the first random sample; the layout's constants were tuned on it, so it is a tuning set.
- `random100-heldout/`: the second; held out at first, then the input of plan 15, so it is a tuning set now.
- `random100-third/`: the third; drawn and run once at 7696947 with no fix aimed at it, so its count is a historical
  prospective evaluation at that commit. It has been looked at since (its entries 9FZL and 6UW4 are named by tuned
  constants, and plan 28 diagnosed every multi-domain failure on it), so it is no unseen set for a change now
  (`docs/model/index.md`, `docs/validation/sensitivity/report.md`).
- the fourth sample, fixed in `docs/validation/fourth-sample-protocol.md`: drawn and run once at the frozen release,
  the paper's headline count. Do not draw it, run it or tune on it.
- the other campaigns (`na/`, `domain-uncertainty/`, `homologues/` and so on) cover their own subsystems.

**Until the fourth sample has been drawn and run**, no new random sample is drawn from entries released on or before
its cutoff (2026-10-13): such a draw could hold the fourth sample's entries, and looking at foldmap's results on any of
them voids it. A change made before then is measured on the regression sets and says so, as plan 28's NA fixes did
(ruling Q37). A development sample (plan 29 onward) is drawn from entries released **after** the cutoff, so it is
disjoint from the fourth sample by construction (ruling Q59), and excludes the earlier samples' ids too:

    python scripts/random100.py sample --seed N --released-after 2026-10-13 --out docs/validation/<name>/sample.json \
        --exclude docs/validation/random100/sample.json --exclude docs/validation/random100-heldout/sample.json \
        --exclude docs/validation/random100-third/sample.json

Such a sample skews to recent entries, and its pools are small for some time: on 2026-10-27 the window held no
membrane (S6) or AlphaFold-model (S12) entry, so a draw then comes up short of those quotas; report the shortfall, do
not widen the window back past the cutoff. Once the fourth sample has run, add `--exclude
docs/validation/random100-fourth/sample.json` as well.

**The measurement** must show the campaign keys and checks identical or better than before the change: the same
entries pass, no check newly fails, the fallback counts and areas do not get worse. Follow `docs/validation/prune/`,
which measures a change on the same entries before and after and lists every difference.

**Where it is recorded:** a committed record under `docs/validation/<name>/` (a `report.md` and its data), stamped
with the source it was made at (`foldmap.source_stamp`). Commit the code first, then make the record at a clean stamp
(no dirty paths), so the record names the commit it measured.

## The inventory

`docs/validation/inventory/inventory.md` lists every constant, threshold and fallback in `src/foldmap` with its
origin, evidence and check (`scripts/inventory.py` scans the code; `tests/test_inventory.py` fails when the table and
the code disagree). Run `python scripts/inventory.py` after a change, then fill the new rows by hand: the origin
(one of biology, geometry, measured, named-entry, convention; `unknown` is not an answer), the evidence (a quoted
comment, a record path, a commit) and the check that would see it go wrong.

A row whose origin is `named-entry:...`, `measured` or `unknown` is a special case or a tuned number.
`docs/validation/inventory/freeze-baseline.json` lists those that existed when the freeze began (by the id less its
line anchor, so reformatting a line does not matter, with how many share each), and `tests/test_freeze.py` fails on
any other such row unless its evidence names an existing `docs/validation/...` record. That record is your
measurement. The record must be a file, and it must name the row (its `module:NAME`, or its constant's or function's
name). Do not edit the baseline to make a failure pass; an id enters it only with the same record. When a special
case is removed, shrink the baseline after the table: run `python scripts/freeze_baseline.py`, which only lowers
counts and drops entries and never adds one. A regeneration that adds an id is the edit the rule forbids;
the test asks for the shrink when the baseline names ids the table no longer has. Within a key the count is applied
first, by id order, and then the record: give a new row its own record.

The scan does not list trivial literals (0, 1, 2, -1, 0.5 and the like; `scripts/inventory.py`), so a threshold
written with one, such as `2 * grid.res` as a reach (plan 28's F4, `na_ends.py`), gets no row and no test sees it.
Review checks every threshold a change writes, trivial literal or not, against this rule.

The test gates only named-entry, measured and unknown. The origin you write for a new row (biology, geometry,
convention) is checked in review, so give a quoted comment or a citation as evidence; the test fails a new row with
such an origin and no evidence.

## Changes admitted before an unseen measurement

Plan 28's nucleic-acid fixes (F2, F4, F6 and F7, `docs/validation/na-fixes/report.md`) each name their mechanism, but
were admitted under ruling Q37 with measurements on the regression sets only (the seven NA cases and 39 controls).
Their unseen measurement is the fourth sample at the frozen release, and one pass count there cannot isolate any one of
them. They are in the freeze baseline because they predate it, not because they met this rule; the na-fixes record
says so. "Every special case has an unseen measurement" is not true of them yet.

## Tests and style

`python -m pytest` and `ruff check . && ruff format --check .` (see the README's Development section). The model
pages are pinned by `tests/test_model_docs.py`.
