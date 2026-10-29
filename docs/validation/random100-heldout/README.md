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

`report.md`, `report.json.gz` (gzipped: `random100.py run --failed-from` reads it as it is) and `run.log` are the run at
the paper's release, v0.2.0 (e1392f7, code fingerprint 7879eea05bd9301c), in a clean checkout of the tag on
zyla-office, judged by the same criteria and limits as the other samples: 94 of 100 pass. Against plan 15's full run
(95), one entry fails that passed: 13FZ, A3 in the domains view (strands B:1-13 and D:1-13 routed by the fallback),
listed for the release after v0.2.0 and not fixed in it (the freeze rule). `review.md` is plan 15's review of plan
15's run: which of the first run's twenty failures passed then, why the other five failed, and what else changed.
