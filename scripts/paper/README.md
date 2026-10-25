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

The steps that compute or draw (domains, runtime, membrane, DSSP, the NA rerun, the determinism Linux leg, the
figures, the typesetting) run here by default. With `PAPER_ZL=/path/to/wrapper` they run through a wrapper with the
interface `zl [--push PATH]... [--pull PATH]... -- COMMAND` (sync this checkout to the same path on another machine,
run COMMAND there in the checkout, copy each `--pull` path back); exit 75 means busy and is retried.

`extract.py` reads git objects at the commit, never the working tree. Every other step reads what it wrote. The steps
that run foldmap put `build/tree/src` (the code at that commit) first on the path and assert it is the foldmap they
imported. `make_numbers.py` and `table_s1.py` read code constants with `ast` and never import foldmap.

## Provenance: every report says which code produced it

`foldmap.source_stamp.source_stamp()` (plan 28 Task 1) returns the git commit and a dirty flag (over `src/`,
`scripts/` and `pyproject.toml`) when the package sits in a git checkout, the code fingerprint
(`app_jobs.code_fingerprint`, a hash of the package's sources, taken with or without git), the writing script's
SHA-256, Python, the platform and the runtime dependencies' versions. Each campaign writes it into its report when it
writes the report, under the key `source_stamp`:

- `scripts/random100.py`: `sample.json`, each entry and `report.json` (and a `Code:` line in `report.md`);
- `scripts/na_validate.py`: `report.json` is `{"source_stamp", "entries"}` (the bare list before; readers take
  both, `common.na_rows`), and `report.md` names the code;
- `scripts/bench.py`: `_meta` of a recorded baseline;
- `scripts/homologues.py`: each pair's result, `usalign.json` and `summary.json`;
- `scripts/sweep.py`: `report.json`;
- here: `domains.py` (`build/domains.json`), `dssp_mkdssp.py` (`run.json`), `membrane_opm.py` (`results.json`), via
  `common.run_stamp` (the harness commit of `build/tree`, plus the tree's own stamp), and `make_figures.py`
  (`build/figures/figures.json`, beside the local HEAD and dirty state `figures.sh` records).

`tests/test_source_stamp.py` checks that each of these writers calls the stamp.

## The environment

`requirements.lock` pins, with hashes, every package the paper's results were measured and drawn with (Python
3.13.13). Install it, then foldmap from the release checkout:

    pip install --require-hashes -r scripts/paper/requirements.lock && pip install --no-deps .

Rebuild it after changing `requirements.in`:

    uv pip compile scripts/paper/requirements.in --universal --python-version 3.13 --generate-hashes -o scripts/paper/requirements.lock

`determinism-constraints.txt` pins the same figure-relevant versions for CI's determinism workflow.

## Determinism and the release gate

`determinism.py` draws every committed test structure's default figure with Liberation Sans forced (`fonts.py`) and
compares layout signatures (as text) and SVG hashes across runs. `.github/workflows/determinism.yml` runs it on
Linux, macOS and Windows on every push to main, and `release.yml` calls it on the commit being released: nothing is
built or published unless it passes there, as the tests must. The Linux leg's three configurations (hash seeds, CPU
features) run with `determinism.py --configurations` (`M2_RERUN=1 run_all.sh`). Read a CI run into the numbers with
`gh run download RUN -p 'determinism-*' -D MANUSCRIPT/data/m2/determinism`, then `make_numbers.py`.

## Table S1

`table_s1.py` writes the parameter table from the constants of the code at the commit (`build/tree`): every value is
read by name (a module constant, a function's default argument, or a number inline in a named function, matched
exactly once), and a value that cannot be found stops the run. `--check-against MANUSCRIPT/supplement.tex` lists
where a typed table differs from the code. `\input{tables/params.tex}` replaces the typed table.

## The random samples

The fourth sample's protocol is `docs/validation/fourth-sample-protocol.md`: seed 20261025, cutoff 2026-10-13, no id
of samples 1-3, drawn and run once at the frozen release; its pass count is the one the paper quotes. The third
sample stays as a historical prospective evaluation at 7696947.
