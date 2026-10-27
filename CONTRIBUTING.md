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

**An unseen set** is a sample you did not use to tune the change. The samples that exist are in `docs/validation/`:

- `random100/`: the first random sample; the layout's constants were tuned on it, so it is a tuning set.
- `random100-heldout/`: the second; held out at first, then the input of plan 15, so it is a tuning set now.
- `random100-third/`: the third; drawn and run once at 7696947 with no fix aimed at it. It was unseen then; if you
  look at its failures to shape a change, it is a tuning set for that change.
- the fourth sample, fixed in `docs/validation/fourth-sample-protocol.md`: no one has seen it, and it is the paper's
  headline count. Do not tune on it.
- the other campaigns (`na/`, `domain-uncertainty/`, `homologues/` and so on) cover their own subsystems.

For a change, draw a fresh sample the same way (`scripts/random100.py sample`, a new seed, none of the earlier ids),
or use one of the above that your change was not tuned on, and say which. Add `--exclude` for each earlier
sample, for example `python scripts/random100.py sample --seed N --exclude docs/validation/random100/sample.json
--exclude docs/validation/random100-heldout/sample.json --exclude docs/validation/random100-third/sample.json`.

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
measurement. Do not edit the baseline to make a failure pass; an id enters it only with the same record. Removing a
special case needs no baseline edit beyond regenerating the table; the test will say when the baseline names ids the
table no longer has.

The test gates only named-entry, measured and unknown. The origin you write for a new row (biology, geometry,
convention) is checked in review, so give a quoted comment or a citation as evidence; the test fails a new row with
such an origin and no evidence.

## Tests and style

`python -m pytest` and `ruff check . && ruff format --check .` (see the README's Development section). The model
pages are pinned by `tests/test_model_docs.py`.
