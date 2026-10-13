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
