# Random-100: the third sample

The first random-100 sample (`docs/validation/random100/`) is the one the layout's constants were tuned on, and the
second (`docs/validation/random100-heldout/`), held out at first, became the input of plan 15, which fixed the
mechanisms behind its failures. Neither count measures how foldmap does on structures it was not tuned on any more.
This sample does: it was drawn after plan 15's fixes, with the same strata, quotas, size bins and search pools
(release-date cutoff 2026-10-13), a new seed and none of the other two samples' 200 ids:

    python scripts/random100.py sample --seed 20261018 --date 2026-10-13 \
        --exclude docs/validation/random100/sample.json \
        --exclude docs/validation/random100-heldout/sample.json \
        --out docs/validation/random100-third/sample.json

(drawn on zyla-lab, where `.git` is not synced; `foldmap_commit` in `sample.json` was set by hand to 7696947, the
commit whose code drew it), and run once with every view and no fix:

    python scripts/random100.py run --views all --jobs 10 \
        --sample docs/validation/random100-third/sample.json --out DIR

`report.md` and `report.json.gz` (gzipped: `random100.py run --failed-from` reads it as it is) are that run, judged
by the same criteria and limits as the other two samples; `review.md` lists each failure with the checks it fails.
The failures are not investigated or fixed here.
