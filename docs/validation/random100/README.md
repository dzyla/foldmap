# Random-100 campaign

100 PDB entries and AlphaFold models of every kind foldmap claims to handle, drawn at random in strata
(plan 10, `docs/superpowers/plans/2026-10-13-10-random-100.md`), each drawn in every view and checked.

## The sample

`sample.json` is written by

    python scripts/random100.py sample --date 2026-10-13

and holds the seed (20261013), the release-date cutoff, the foldmap commit it was drawn at, each stratum's
search query with its pools (count and SHA-256 of the sorted ids, per size bin and per extra quota), the 100
entries (id, stratum, bin, position in the pool's seeded order, extra quota or pair), the replacements and the
skipped candidates with their reasons, and every quota with the value reached.

The draw is reproducible: each pool is sorted and drawn in the order `random.Random("20261013:<stratum>:<bin>")`
gives it (a string seed is hashed with SHA-512, so `PYTHONHASHSEED` does not change it). Drawing again with the
same date gives the same draw as long as the pools and the downloaded files are unchanged (a candidate's file
decides whether it is CA-only, holds an N-linked glycan or counts for a quota); the file itself also records the
foldmap commit it was drawn at, so it is not the same byte for byte. Samples drawn since the final review also
record the gemmi version and the SHA-256 of every candidate file read (this one predates that).
`python scripts/random100.py sample --check` lists the sampled ids that are obsolete now, without changing
the file.

## The held-out and third samples

The layout's constants were tuned on this sample (plan 10, Task 12), so its pass rate measures fit. A second sample
that no constant was tuned on is in `docs/validation/random100-heldout/`: drawn with

    python scripts/random100.py sample --seed 20261017 --date 2026-10-13 \
        --exclude docs/validation/random100/sample.json --out docs/validation/random100-heldout/sample.json

(the same strata, quotas and pools, none of this sample's ids), and run once with every view and no fix (80 of 100).
Plan 15 then fixed the mechanisms behind its failures, so it is held out no longer. A third sample, drawn after
those fixes with none of either sample's ids (`--exclude` given twice, seed 20261018), is the one foldmap was not
tuned on: `docs/validation/random100-third/`. In plan 15's full run (7696947) the three pass 93 (this sample), 95
(the second) and 83 (the third) of 100.

`report.md`, `report.json` and `run.log` are the run at the paper's release, v0.2.0 (e1392f7, code fingerprint
7879eea05bd9301c), in a clean checkout of the tag on zyla-office: 91 of 100 pass. Against plan 15's run (93), 10UJ
and 9YKJ pass now, and 8RRX and 5LDW fail A6: a label collides with the helix it names (8RRX G:4806-4819 in every
view; 5LDW G:265-268, G:296-308 and G:389-401 in the deposited-names view). Both are listed for the release after
v0.2.0, not fixed in it (the freeze rule). `review.md` is plan 15's review of plan 15's run.

## The run

    NA_RNAPOLIS_PYTHON=/path/to/venv-with-rnapolis/bin/python python scripts/random100.py run [--jobs 1]
    python scripts/random100.py run --failed-from docs/validation/random100/report.json   # again, the failures
    python scripts/random100.py run --ids 1ABC,2DEF --views quick --out /tmp/r100           # a few, quickly

Each entry is drawn in every view that applies (default, `--names deposited`, `--names beta`, focus,
`--domains auto`, `--membrane auto` for membrane proteins, `PYTHONHASHSEED=1`, a saved layout redrawn, the
interactive page, the sequence plot and the app's hosted job), each view in its own process, and judged on the
criteria A1-A20 of the plan with the limits in `src/foldmap/criteria.py: LIMITS` (`foldmap check` judges one figure
by the same criteria). Each criterion, and each science check S1-S6, reads pass, fail, not judged (it applies, but
the entry gave it nothing to judge: no tie between copies, no ligand tether to another chain, a domains view that drew
no domain, mkdssp unmeasured, no sheet drawn) or does not apply (no nucleic acid for A17, no membrane for A18). An
entry passes when no criterion fails; nothing to judge is never counted as a pass, and `report.md` gives each
criterion's four counts and lists every case not judged.

What the checks do not measure, or measure partly by construction:

- A1c (a CA-only or P-only entry refused in one line) met no such entry in the sample: 0 entries.
- A6 reads text against loops, nucleic-acid runs and ties, tethers, disulfide bars, glycan stems, H-bonds, chain-end
  stubs and domain leaders; a domain name over another domain's tinted panel is not measured (A16 measures panels
  apart and names against text; ruling 39 leaves the rest to the visual review).
- A18 is also the test the `--membrane` view lays itself out by (it keeps whichever layout scores better), so its
  passes are partly by construction.
- A1's Python warnings are those of each view's own process; the app view draws in the app's job process, whose
  warnings are not known to reach the view's stderr.
- A2 gates the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis); `report.md` gives
  that time and the drawing time inside it apart.

Written here:

- `report.md` and `report.json` (committed): the pass count, a table per stratum and per entry, what failed, the
  statistics reported but not gated, the limits, versions, the foldmap commit and the SHA-256 of `sample.json`;
- `run.log` (committed): one line per entry and run;
- `<id>.png`, `index.html`, `<id>.json` and `work/` (kept local, gitignored): the contact sheets (the figure, the
  3D traces seen through the figure's page axes in its element colours, the same turned 90° about page-up, the
  partner's figure for a pair, then the facts and every criterion, failures in red), a grid linking them, each
  entry's measures and each view's own files.
