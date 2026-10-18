# Random-100: the held-out sample

The layout's constants were tuned on the first random-100 sample (`docs/validation/random100/`, plan 10 Task 12), so
its pass rate measures fit, not how well the drawing generalises (final review A, Important 2). This sample was drawn
with the same strata, quotas, size bins and search pools (release-date cutoff 2026-10-13), with a new seed and none of
the first sample's 100 ids:

    python scripts/random100.py sample --seed 20261017 --date 2026-10-13 \
        --exclude docs/validation/random100/sample.json --out docs/validation/random100-heldout/sample.json

and run once, with every view and no fix:

    python scripts/random100.py run --views all --sample docs/validation/random100-heldout/sample.json --out DIR

`report.md` and `report.json.gz` (gzipped: `random100.py run --failed-from` reads it as it is) are that run, judged
by the same criteria and limits as the tuned sample; `review.md` lists each failure with the checks it fails. The
failures are not investigated or fixed here: they are the next plan's input.
