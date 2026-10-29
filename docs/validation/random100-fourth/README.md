# Random-100: the fourth sample

The paper's headline count. The first three random samples cannot measure foldmap on structures it was not tuned on
any more. The first tuned the layout's constants. The second became the input of plan 15. The third was run once at
7696947, and the code changed after it. This sample was drawn and run once at the frozen release, under
`docs/validation/fourth-sample-protocol.md`. It uses the same strata, quotas, size bins and search pools as the other
three (release-date cutoff 2026-10-13), the seed 20261025 fixed in the protocol, and none of the 300 ids of samples
1-3.

It was drawn on zyla-office on 2026-10-28, in a clean shallow clone of tag v0.2.0 (commit e1392f7), after the release
workflow (tests and the three-OS determinism check) had passed on the tag:

    python scripts/paper/fourth_sample.py draw --release v0.2.0

`sample.json` records the protocol: the seed, the cutoff, each excluded file with its hash, the pools' counts and
hashes, and the source stamp of the code that drew it. The sample was then run once with every view, in the same
checkout:

    python scripts/random100.py run --views all --jobs 28 \
        --sample docs/validation/random100-fourth/sample.json --out docs/validation/random100-fourth

The environment was Python 3.13.13 with `scripts/paper/requirements.lock` installed (`--require-hashes`), and foldmap
installed with `--no-deps` from that checkout. The report's source stamp says the installed code matches the
checkout: fingerprint 7879eea05bd9301c, git e1392f7, not dirty. The run used mkdssp 3.1.2 and RNApolis 0.16.17 from
its own venv (`NA_RNAPOLIS_PYTHON`). The run did not record its wall time. File times put it at about 1,190 s.

These files were committed as the run left them (3990050):

- `report.md`;
- `report.json.gz`, the 20 MB `report.json` gzipped with `gzip -9n`;
- `run.log`, one line per entry as it finished;
- `sample.json`.

85 of the 101 entries pass. The 101 are the 100 drawn entries and 2CSE, a CA-only entry run once for A1c.
`review.md` lists each failure by check, with its cause where one was established. Nothing was fixed or rerun after
the run. The run's figures and per-entry files are on zyla-office, in the run folder, and are not committed (1.6 GB).
