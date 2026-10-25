# The paper pipeline (scripts/paper)

Regenerates every number key, table and figure of the foldmap paper from one commit (a release passes its tag).
Moved here from the local manuscript folder by plan 28 Task 1 (plan 14 M1-M4 built it). The manuscript text, its
inputs and every output stay outside version control, in the manuscript folder: `docs/superpowers/work/manuscript/`
of this checkout (gitignored), or the folder `$FOLDMAP_MANUSCRIPT` names.

| Where | What |
|---|---|
| `MANUSCRIPT/main.tex`, `supplement.tex`, `refs.bib` | the text (read by `fill.py`, `crosscheck.py`, `render.sh`) |
| `MANUSCRIPT/data/` | raw measurement outputs (`m2/`: DSSP, membrane, NA rerun, determinism; `opm/`: OPM's files) |
| `MANUSCRIPT/figures/data/` | figure inputs (the CI determinism artifacts Fig. 3d reads) |
| `MANUSCRIPT/build/` | every output: `numbers.json`, `tables/`, `figures/`, the filled `.tex`, `pdf/` |

## Commands

| Output | Command |
|---|---|
| everything except the figures, in order | `scripts/paper/run_all.sh [COMMIT]` (`M2_RERUN=1` reruns the M2 measurements first) |
| `build/sources/`, `build/tree/`, `build/provenance.json` | `python scripts/paper/extract.py --commit COMMIT` (needs `.git`) |
| `build/domains.json` (CATH domain scores, code from `build/tree`) | `python scripts/paper/domains.py` |
| `build/runtime.*`, `build/tables/runtime.csv` | `python scripts/paper/runtime.py` |
| `build/numbers.json`, `build/tables/{r100,na,domains,dssp,membrane}.tex` | `python scripts/paper/make_numbers.py` |
| `build/tables/params.tex`, `params.json` (Table S1 from the code's constants) | `python scripts/paper/table_s1.py` |
| `build/main-filled.tex`, `supplement-filled.tex`, `keys.json` | `python scripts/paper/fill.py` |
| `build/crosscheck.json` (numbers in the prose against the sources) | `python scripts/paper/crosscheck.py` |
| `build/figures/` (Figs. 1-5, the supplementary gallery, `checks.json`) | `scripts/paper/figures.sh [TARGET ...]` |
| `build/pdf/` | `scripts/paper/render.sh` |
| M2: `data/opm/` | `python scripts/paper/opm_fetch.py` (network) |
| M2: `data/m2/membrane/results.json` | `python scripts/paper/membrane_opm.py` |
| M2: `data/m2/dssp/` | `python scripts/paper/dssp_mkdssp.py` (mkdssp 3.1.2 as `$MKDSSP`) |
| cross-platform determinism | `scripts/paper/determinism.py` (`.github/workflows/determinism.yml`; see below) |
| the fourth random sample, once, at the release | `python scripts/paper/fourth_sample.py plan`, then `draw --release TAG` |
| whether a report's stamps name the code of the commit they claim | `python scripts/paper/verify_stamp.py REPORT` (needs `.git`) |

The steps that compute or draw (domains, runtime, membrane, DSSP, the NA rerun, the determinism Linux leg, the
figures, the typesetting) run here by default. With `PAPER_ZL=/path/to/wrapper` they run through a wrapper with the
interface `zl [--push PATH]... [--pull PATH]... -- COMMAND` (sync this checkout to the same path on another machine,
run COMMAND there in the checkout, copy each `--pull` path back); exit 75 means busy and is retried.

`extract.py` reads git objects at the commit, never the working tree. Every other step reads what it wrote. The steps
that run foldmap put `build/tree/src` (the code at that commit) first on the path and assert it is the foldmap they
imported. `make_numbers.py` and `table_s1.py` read code constants with `ast` and never import foldmap.

## Provenance: every report says which code produced it

`foldmap.source_stamp.source_stamp(__file__)` (plan 28 Task 1) is taken when a campaign process starts and written
with each report, entry and record it measures, under the key `source_stamp`. It holds:

- the code fingerprint of the foldmap that ran (`app_jobs.code_fingerprint`, a hash of the package's sources, with or
  without git), and whether it equals the fingerprint of the script's checkout's `src/foldmap`
  (`installed_matches_checkout`): a non-editable install that differs from the checkout is never taken for its commit;
- `git`, found from the writing script's own path (never from where foldmap is imported): HEAD's commit and a dirty
  flag over `src/`, `scripts/`, `pyproject.toml` and `tests/data` ("source": "git"); on a runner's synced copy with
  no `.git`, the commit and dirty state the runner exports as `FOLDMAP_GIT_COMMIT` / `FOLDMAP_GIT_DIRTY`
  ("source": "asserted by runner");
- the script's path and SHA-256, Python, the platform and the runtime dependencies' versions.

Where it is written:

- `scripts/random100.py`: `sample.json`, each entry (also one synthesised for an entry process that died) and
  `report.json`; the report says which codes measured its entries (`entry_code`) and `mixed_code` when more than one
  did; `report.md` names both;
- `scripts/na_validate.py`: `report.json` is `{"source_stamp", "entry_code", "mixed_code", "entries"}` with a stamp
  per entry (the bare list before; readers take both, `common.na_rows`);
- `scripts/bench.py`: `_meta` of a recorded baseline;
- `scripts/homologues.py`: each row of `align.json`, each pair's result, `usalign.json` and `summary.json` (with
  `entry_code` / `mixed_code`);
- `scripts/sweep.py`: `report.json`;
- here: `domains.py`, `dssp_mkdssp.py` (`run.json` and each entry's file), `membrane_opm.py` (via `common.run_stamp`:
  the harness commit of `build/tree` and the tree's stamp, computed by this checkout's code for a tree older than the
  stamp), `determinism.py` (each run record), `make_figures.py` (each panel in `panels.json`, and `figures.json` with
  the local HEAD and dirty state `figures.sh` records and the external tools' versions).

`tests/test_source_stamp.py` finds every JSON or Markdown writer in `scripts/` by its call shape (`json.dump`, and
`write_text` / `write_bytes` of JSON or of a `.json` / `.md` path) and fails unless the writing function takes or
uses the stamp, or is on its short list of non-report writers, each with its reason.

**Checking a stamp.** `python scripts/paper/verify_stamp.py REPORT [--commit SHA]` exports `src/foldmap` at the commit
each stamp names (`git archive`) and computes its fingerprint with the same function; a stamp is verified when the
two agree (and the script's text at that commit hashes the same). `extract.py` reads the code a source was measured
with from its stamps only, verifies each, and stops on a stamp that does not verify, on `mixed_code` and on a report
marked `not_for_paper`; `make_numbers.py` refuses the last two as well. Files from before the stamp say so: a commit
written into them by hand after the run is kept as legacy and never used.

## The environment

`requirements.lock` pins, with hashes, every Python package the paper's results were measured and drawn with
(Python 3.13.13). Install it, then foldmap from the release checkout:

    pip install --require-hashes -r scripts/paper/requirements.lock && pip install --no-deps .

Rebuild it after changing `requirements.in`:

    uv pip compile scripts/paper/requirements.in --universal --python-version 3.13 --generate-hashes -o scripts/paper/requirements.lock

CI's determinism workflow installs the same lock. The tools outside Python, not in the lock:

| Tool | Version | Used by | Recorded in |
|---|---|---|---|
| mkdssp | 3.1.2 (`$MKDSSP`) | random100 A10, `dssp_mkdssp.py` | `report.json` versions, `run.json` (asserted 3.1.2) |
| RNApolis | its venv's (`$NA_RNAPOLIS_PYTHON`) | `na_validate.py`, random100 A17 | random100 `report.json` versions |
| US-align | as found (`$USALIGN`) | `homologues.py align --usalign` | `usalign.json` |
| Tectonic | 0.15.0 | `render.sh`, `make_figures.py` | `figures.json` tools |
| poppler (pdftocairo, pdftoppm) | the runner's | `make_figures.py` | `figures.json` tools |
| setuptools | `>=77` (pyproject's build backend, not pinned) | `pip install .` | the package version |

## Determinism and the release gate

`determinism.py` draws every committed test structure's default figure with Liberation Sans forced (`fonts.py`) and
compares layout signatures (as text) and SVG hashes across runs. `.github/workflows/determinism.yml` runs it on
Linux, macOS and Windows on every push to main, and on the tagged commit of every release: `release.yml` (started by
pushing a `v*` tag; see `RELEASING.md`) builds, publishes and creates the GitHub release (which Zenodo archives) only
after it and the tests pass. The Linux leg's three configurations (hash seeds, CPU features) run with
`determinism.py --configurations` (`M2_RERUN=1 run_all.sh`). Read a CI run into the numbers with
`gh run download RUN -p 'determinism-*' -D MANUSCRIPT/data/m2/determinism`, then `make_numbers.py`.

## Table S1

`table_s1.py` writes the parameter table from the constants of the code at the commit (`build/tree`): every value is
read by name (a module constant, a function's default argument, or a number inline in a named function, matched
exactly once), and a value that cannot be found stops the run. `--check-against MANUSCRIPT/supplement.tex` lists
where a typed table differs from the code. `\input{tables/params.tex}` replaces the typed table.

## The random samples

The fourth sample's protocol is `docs/validation/fourth-sample-protocol.md`: seed 20261025, cutoff 2026-10-13, no id
of samples 1-3, drawn once in a clean checkout of the release tag (the sampler refuses anywhere else) and run once,
whole, with that code (`random100.py run` refuses a second run, `--ids`, `--failed-from` or other code on a sample
carrying a protocol record, unless `--not-for-paper`, which marks the report never to be quoted); its pass count is
the one the paper quotes. The third sample stays as a historical prospective evaluation at 7696947.
