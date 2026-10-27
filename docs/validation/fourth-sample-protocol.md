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
| Release-date cutoff | **2026-10-13**, the cutoff of samples 1-3: the same rule. The pools are what the RCSB search service returns on the draw date (they drift: samples 1 and 3 already differ by 1-2 ids in several pools), recorded with their counts and hashes in `sample.json` |
| Excluded | every id of samples 1-3 (entries and their CA-only / P-only replacements), read from their committed `sample.json` files, which must still hash to the SHA-256 values pinned in the script (`SAMPLE_SHA256`) |
| Output | `docs/validation/random100-fourth/sample.json`, with a `protocol` record (seed, cutoff, rule, each excluded file with its hash and id count, the release tag, and `foldmap.source_stamp` of the code that drew) |

Homologues of the excluded entries are not excluded: the exclusion is by PDB id, as for the third sample. A drawn
entry may therefore share a protein with an entry of samples 1-3; the paper states this limit beside the count.

## When and how it runs

Once, at the frozen release (the tagged commit the paper describes), after the release gate has passed
(`.github/workflows/release.yml`: tests and the cross-platform determinism workflow on that commit):

1. Check the choices without drawing (no network, writes nothing):

       python scripts/paper/fourth_sample.py plan

2. Draw it (network), in a clean git checkout whose HEAD is the release tag's commit, with that checkout's foldmap.
   The script refuses anywhere else, writes only to its one fixed path, and refuses if the sample file exists, so it
   cannot be re-drawn:

       python scripts/paper/fourth_sample.py draw --release vX.Y.Z

3. Run it with every view, on the same release, in the paper environment (`scripts/paper/requirements.lock`):

       python scripts/random100.py run --views all --jobs N \
           --sample docs/validation/random100-fourth/sample.json --out docs/validation/random100-fourth

   For a sample carrying a protocol record, `random100.py` refuses, unless `--not-for-paper` is given: the sample read
   from any path but its fixed one (`random100.PAPER_SAMPLES`, so a copy elsewhere is refused too); a run written
   anywhere but beside it, in `docs/validation/random100-fourth`; a second run (a report or run log already in that
   folder, whatever `--out` says); a partial run (`--ids`); a rerun (`--failed-from`); fewer views than all; and code
   other than the draw's (by code fingerprint). `--not-for-paper` marks the report never to be quoted, and
   `scripts/paper/extract.py` refuses it, as it refuses a report whose entries were measured by different code
   (`mixed_code`). The guard holds in one checkout; a run in another clone is outside it, and the rule that the one
   run is made in the release checkout and committed as it is (step 4) is the guard there.

4. Commit `sample.json`, `report.json` (gzipped when over 5 MB, as for the other samples), `report.md` and the run log
   as they are. Every report and entry carries its own `source_stamp` (taken as its process starts).

## Rules

- **No fix after the run.** No code, limit or criterion changes between the draw and the quoted count. A failure is
  reported with its cause in `review.md`; it is not fixed and rerun for the paper.
- **No reruns for the count.** An entry whose run did not finish (a failed download, the machine stopping)
  keeps its recorded failure in the quoted count, with the reason in `review.md`. A `--not-for-paper` rerun may
  diagnose it, but its report is never quoted.
- **Criteria and limits** are those of `foldmap.criteria` at the release (`CRITERIA`, `LIMITS`), the same as the
  report records.
- **The quoted number** is this sample's pass count out of its entries, with its 95% Wilson interval (plan 28 Task 2),
  and the per-criterion outcomes beside it. Samples 1-3 are reported apart, each labelled: tuning, regression input,
  historical prospective evaluation at 7696947.
- If the draw comes up short of a quota or of 100 entries, the sample is still the one run and reported, with the
  shortfall stated; it is not re-drawn with another seed.
- **No other draw from its pools before it runs.** Until this sample is drawn and run, no other random sample is drawn
  from releases on or before its cutoff (2026-10-13): such a draw could hold its entries, and looking at foldmap's
  results on them voids it (next rule). A development sample (plan 29 onward) is drawn from entries released after
  the cutoff (`random100.py sample --released-after 2026-10-13`), so it is disjoint from this one by construction
  (ruling Q59; `CONTRIBUTING.md`).
- **The seed is public.** It is fixed here, in the open, so anyone can compute the sample before the release with
  `random100.py sample --seed 20261025 --exclude ...`. Drawing it, or looking at foldmap's results on any entry of it,
  before the release by any route voids the sample: the protocol is then rewritten with a new seed, and the reason is
  stated. (The draw and run are recorded with their source stamps, and the draw refuses anything but the release
  checkout.)
