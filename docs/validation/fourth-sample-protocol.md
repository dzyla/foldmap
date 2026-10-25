# Random-100: the fourth sample (protocol)

Written 2026-10-25, before any draw, for plan 28 Task 1 (the external review's point 3.9). It fixes what the fourth
random sample is, how it is drawn and run, and which number the paper quotes from it. Nothing here may change after
the draw.

## Why a fourth sample

The paper reports how foldmap does on structures it was not tuned on. The three random samples so far cannot do that
any more:

- the **first** (`random100/`, seed 20261013) is the one the layout's constants were tuned on;
- the **second** (`random100-heldout/`, seed 20261017) was held out at first, then became the input of plan 15, which
  fixed the mechanisms behind its failures;
- the **third** (`random100-third/`, seed 20261018) was drawn and run once at 7696947, with no fix aimed at it. Code
  changed after that run (e5b4bae, then plans 16-28), so its count describes 7696947, not the release.

The third sample is kept and reported as what it is: a **historical prospective evaluation at 7696947**. The paper's
headline count comes from a fourth sample that no one has seen, run once on the frozen release.

## The sample

| | |
|---|---|
| Sampler | `scripts/random100.py`'s stratified sampler, unchanged at the release: the same strata, quotas, size bins, RCSB search pools and look-ahead (`LOOK`) as samples 1-3 |
| Seed | **20261025**, fixed now in `scripts/paper/fourth_sample.py` (`SEED`) |
| Release-date cutoff | **2026-10-13**, the cutoff of samples 1-3, so the population is the same, less their ids |
| Excluded | every id of samples 1-3 (entries and their CA-only / P-only replacements), read from their committed `sample.json` files, which must still hash to the SHA-256 values pinned in the script (`SAMPLE_SHA256`) |
| Output | `docs/validation/random100-fourth/sample.json`, with a `protocol` record (seed, cutoff, rule, each excluded file with its hash and id count, the release tag, and `foldmap.source_stamp` of the code that drew) |

Homologues of the excluded entries are not excluded: the exclusion is by PDB id, as for the third sample. A drawn
entry may therefore share a protein with an entry of samples 1-3; the paper states this limit beside the count.

## When and how it runs

Once, at the frozen release (the tagged commit the paper describes), after the release gate has passed
(`.github/workflows/release.yml`: tests and the cross-platform determinism workflow on that commit):

1. Check the choices without drawing (no network, writes nothing):

       python scripts/paper/fourth_sample.py plan

2. Draw it (network). The script refuses if the sample file exists, so it cannot be re-drawn:

       python scripts/paper/fourth_sample.py draw --release vX.Y.Z

3. Run it with every view, on the same release, in the paper environment (`scripts/paper/requirements.lock`):

       python scripts/random100.py run --views all --jobs N \
           --sample docs/validation/random100-fourth/sample.json --out docs/validation/random100-fourth

4. Commit `sample.json`, `report.json` (gzipped when over 5 MB, as for the other samples), `report.md` and the run log
   as they are. Every report and entry carries its own `source_stamp` (written when it is written).

## Rules

- **No fix after the run.** No code, limit or criterion changes between the draw and the quoted count. A failure is
  reported with its cause in `review.md`; it is not fixed and rerun for the paper.
- **Reruns** are allowed only for an entry whose run did not finish for a reason outside foldmap (a download that
  failed, the machine stopping). Each such rerun is listed in `review.md` with the reason; the entry is rerun on the
  same release.
- **Criteria and limits** are those of `foldmap.criteria` at the release (`CRITERIA`, `LIMITS`), the same as the
  report records.
- **The quoted number** is this sample's pass count out of its entries, with its 95% Wilson interval (plan 28 Task 2),
  and the per-criterion outcomes beside it. Samples 1-3 are reported apart, each labelled: tuning, regression input,
  historical prospective evaluation at 7696947.
- If the draw comes up short of a quota or of 100 entries, the sample is still the one run and reported, with the
  shortfall stated; it is not re-drawn with another seed.
