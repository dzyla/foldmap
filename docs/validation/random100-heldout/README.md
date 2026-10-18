# Random-100: the second sample (held out until plan 15)

The layout's constants were tuned on the first random-100 sample (`docs/validation/random100/`, plan 10 Task 12), so
its pass rate measures fit, not how well the drawing generalises (final review A, Important 2). This sample was drawn
with the same strata, quotas, size bins and search pools (release-date cutoff 2026-10-13), with a new seed and none of
the first sample's 100 ids:

    python scripts/random100.py sample --seed 20261017 --date 2026-10-13 \
        --exclude docs/validation/random100/sample.json --out docs/validation/random100-heldout/sample.json

and run once with every view and no fix: 80 of 100 passed. Plan 15 then fixed the mechanisms behind those twenty
failures, checking each fix on these entries, so **this sample is no longer held out**: its count now measures fit
too. The count for structures foldmap was not tuned on is the third sample's (`docs/validation/random100-third/`),
drawn after plan 15's fixes with none of this sample's or the first's ids.

    python scripts/random100.py run --views all --sample docs/validation/random100-heldout/sample.json --out DIR

`report.md` and `report.json.gz` (gzipped: `random100.py run --failed-from` reads it as it is) are plan 15's full run
(95 of 100), judged by the same criteria and limits as the other samples; `review.md` says which of the first run's
twenty failures pass now, why the other five fail, and what else changed.
