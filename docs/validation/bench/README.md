# Bench baseline at v0.2.0

The gate's layout bench (`scripts/bench.py --check BASELINE --exact`: every test structure and the synthetic
32-chain complex laid out, each layout's values compared with a kept one) compared against a baseline recorded in
plan 15 (`.superpowers/sdd/2026-10-17-15-heldout-failures/bench-before.json`, local, foldmap 0.1.0). Plan 29 let a
fix differ from it only on bench entries the fix's record named and explained (ruling P29-1), and regenerates it once,
at the freeze, here.

- `check-plan15.txt`: v0.2.0 against the plan 15 baseline, `--exact`: 22 of 22 identical (exit 0). Plan 29's fix
  records explained no bench entry, since none differed.
- `baseline.json`: the baseline recorded at v0.2.0 (`bench.py --record`, code fingerprint 7879eea05bd9301c, stamped in
  `_meta`), checked against a second run of the same code: 22 of 22 identical. It is the gate's baseline from now on
  (copied to the local path above).

The two baselines are equal value for value. As text they differ in 49 signed zeros (−0.0 in plan 15's, 0.0 now),
which `--exact` counts equal (`bench.py --help`); `scripts/paper/determinism.py` compares signatures as text and
would tell them apart.

Run on zyla-office (gpu-desktop unreachable), 2026-10-29.
