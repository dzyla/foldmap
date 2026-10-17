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
same date gives the same file byte for byte, as long as the pools are unchanged.
`python scripts/random100.py sample --check` lists the sampled ids that are obsolete now, without changing
the file.

## The run

    NA_RNAPOLIS_PYTHON=/path/to/venv-with-rnapolis/bin/python python scripts/random100.py run [--jobs 1]
    python scripts/random100.py run --failed-from docs/validation/random100/report.json   # again, the failures
    python scripts/random100.py run --ids 1ABC,2DEF --views quick --out /tmp/r100           # a few, quickly

Each entry is drawn in every view that applies (default, `--names deposited`, `--names beta`, focus,
`--domains auto`, `--membrane auto` for membrane proteins, `PYTHONHASHSEED=1`, a saved layout redrawn, the
interactive page, the sequence plot and the app's hosted job), each view in its own process, and judged on the
criteria A1-A20 of the plan with the limits in `src/foldmap/criteria.py: LIMITS` (`foldmap check` judges one figure
by the same criteria). Written here:

- `report.md` and `report.json` (committed): the pass count, a table per stratum and per entry, what failed, the
  statistics reported but not gated, the limits, versions, the foldmap commit and the SHA-256 of `sample.json`;
- `run.log` (committed): one line per entry and run;
- `<id>.png`, `index.html`, `<id>.json` and `work/` (kept local, gitignored): the contact sheets (the figure, the
  3D traces seen through the figure's page axes in its element colours, the same turned 90° about page-up, the
  partner's figure for a pair, then the facts and every criterion, failures in red), a grid linking them, each
  entry's measures and each view's own files.
