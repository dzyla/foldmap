# Freeze checklist

What must be done at the frozen release, before the paper quotes a number, and what plan 28's final review accepted
and wrote down instead of fixing (plan 28 final review, rulings Q59 and Q60). The release itself is `RELEASING.md`;
the fourth sample is `docs/validation/fourth-sample-protocol.md`; the paper pipeline is `scripts/paper/README.md`.

## 1. Records regenerated at the release

Each record below was made at a commit before the release; code changed after several of them. Each is made again
at the release commit, the code committed first and the record at a clean stamp in its own commit, as for every
record. A record whose number the paper quotes is regenerated before the paper quotes it.

| record | made at | why it is regenerated | how |
|---|---|---|---|
| `random100-fourth/` | drawn and run once at v0.2.0 (85 of 101) | the paper's headline count | `fourth_sample.py draw --release TAG`, then one `random100.py run --views all`, in the release checkout (the protocol) |
| `sensitivity/` | runs at 5f09684's code | the NA fixes, the prune and the rest changed what it measures; the default baseline moved 80 to 82 (its "Code since this run" section, final review F1) | rerun stages B, V and X at the release commit (about 3.9 h on zl), then `sensitivity.py report` |
| `readability-floor/` | e6fcf36 | its census and measures predate later code; its generator now says R1's route-budget reason fires on 7NZ2's domains view (F9) | `scripts/readability_floor.py` |
| `prune/` | runs at 46f5ebf, record b730a05 | its next render adds one sentence: the record's stamp is b730a05 while the equivalence runs were at c64ea7a, and `git diff c64ea7a b730a05 -- src` is comments and docstrings only (parked T9.3) | `scripts/prune.py` report at the release |
| `metamorphic/` | this round's re-render | stamped at a pre-release commit | `scripts/metamorphic_record.py` |
| `na-fixes/` | 6ca884c | stamped at a pre-release commit; `na_runs:_JOINT_NEAR` was pruned after it (the prune measured the NA set identical) | `scripts/na_fixes_record.py` |
| `domain-uncertainty/` | 42e3fe8 | stamped at a pre-release commit | `scripts/domain_uncertainty.py` |
| `homologues/topology-transfer/` | cdd29dd | H1-H5 not remeasured since; low risk (domain-name placement and NA joints only changed) | `scripts/topology_transfer.py`, with the other homologue campaigns |
| `random100/`, `random100-heldout/` | 7696947; rerun at v0.2.0 (91 and 94 of 100) | "on this release" in the paper is false until they are rerun (claims R8, R9) | `random100.py run` at the tag, not for the headline |
| `random100-third/` | 7696947 | kept as it is: a historical prospective evaluation at 7696947 | none |
| `na/` (and the M2 NA rerun) | older code | the NA campaign the paper quotes | `M2_RERUN=1 scripts/paper/run_all.sh TAG HARNESS` |
| `inventory/` and its freeze baseline | b730a05 | the table must match the release's code (`tests/test_inventory.py` fails otherwise) | `scripts/inventory.py`; `scripts/freeze_baseline.py` only shrinks |
| the paper's numbers, tables and figures | 838c579 and later | everything quoted | `scripts/paper/run_all.sh TAG HARNESS` (HARNESS: the commit adding the fourth sample), `figures.sh`; inputs outside git rebuilt first (`scripts/paper/README.md`) |

### Plan 29's records (`docs/validation/app29/`)

Plan 29 (the app on a fresh random sample, 2026-10-27/28) added these records. The paper quotes none of them; each
is a before/after measurement of one plan, kept as it is at its commit (as `random100-third/` is), and rerun at the
release only if the paper comes to quote it (`scripts/app29.py run --ui --record`, then the record's builder).

| record | what | made at (record commit, stamp) |
|---|---|---|
| `sample.json` | the sample (79 entries, released 2026-10-13 to 10-27) and its tune/check split | 851144d |
| `before/` | the app's code path on both halves before the fixes (the check half sealed in `report-check.json`) | f538981, 549bef9 |
| `triage.md` | the tune half's faults by mechanism, ranked; ruling P29-7's fix list | 2319029 |
| `fix-ligand/` | Task 4: ligand marks clear of loop lines; the tether-reach step | 5880df6, fa95efc |
| `fix-floor/` | Task 5: focus auto takes one subunit of a two-copy complex below the floor | 3d6c798, f62c899 |
| `fix-screw/` | Task 6: the screw search without n^4 pair fits (figures byte-identical) | 6b83673, 4ccb38f |
| `fix-ghost/` | Task 7: a lent ghost is no barrel strand (S3 only) | cf16ce9, e8f1e47 |
| `diagnosis-t8.md` | Task 8: 26EL's empty page and 22MZ's letters, diagnosed | cdd315f |
| `fix-copies/` | Task 8b: copies modelled in part are one protein to focus auto (26EL) | 749648c, 1a15466 |
| `fix-letters/` | Task 8b: a base letter nearest its own strand (22MZ): rejected and reverted, the patch kept | f8388fc, 8919aa3 |
| `after/` | Task N+1: the after run, the before commit rerun on the same host, the before/after table per half (the check half opened once) and the parked measures (third sample, NA controls, science over the download cache) | 657e127; the run at 162e067, its table at b84ebc3 |
| `look.md` | Task N+1: the look test's before/after pairs, their sha256, the controller's verdicts | 657e127, b84ebc3 (verdicts by hand) |

**The bench baseline.** `bench-before.json` (`.superpowers/sdd/2026-10-17-15-heldout-failures/`, the gate's
`bench.py --check --exact`) is regenerated once, at the freeze, in its own commit with a record listing every entry
the fix records explained (ruling P29-1); plan 29's fix records explained none, since the bench stayed 22 of 22
identical through every fix and at the after measurement.

`tests/test_domain_confidence.py` pins the README's and the MCP docstring's "S5 failures from 12 to 19" to the
sensitivity record: a rerun that moves it fails there until every copy is changed.

## 2. Accepted and written down

Rows the final review's parked table ruled "accept" or "accept and document", and two findings accepted as limits.
Each is a known gap of the tooling, not of foldmap's figures.

- **Inventory ids repeat by order (T9.1).** Two equal literals in one function take `#k` by their order; swapping
  them swaps their ids and so their hand columns. This is inherent to line-anchored ids: the merge's `carried` print
  shows every carried row, and review catches a swap.
- **Delete and add of one key (T9.1).** Removing a literal and adding another with the same key (module, function,
  kind and value) in one change carries the old row's hand columns to the new one. Inherent to the key (ruling Q25);
  review reads the carried rows.
- **The freeze gate orders twins by id (T9.5).** Within a key the baseline's count is applied first, by id order, and
  then the record, so a new twin of a baseline row needs its own record. Written in `CONTRIBUTING.md`.
- **Trivial literals are not listed (F13).** 0, 1, 2, -1 and 0.5 get no inventory row, so a threshold written with
  one (such as F4's `2 * grid.res` reach in `na_ends.py`) is seen by no test. Review checks such thresholds
  (`CONTRIBUTING.md`).
- **Changes admitted before an unseen measurement (F14).** The NA fixes F2, F4, F6 and F7 were admitted on
  regression-set measurements (ruling Q37); their unseen measurement is the fourth sample, which cannot isolate them.
  Written in `CONTRIBUTING.md` and the na-fixes record.
- **MF2 is pinned per fixture, not per pair (T7).** `tests/test_independent_kernels.py` pins MF2's count per fixture;
  the record's test pins the pairs themselves where the docs ship. One more test per pair adds little.
- **The prune's equivalence runs are not committed (T9.3).** Their JSON is in the record; the run folders are bulky.
  The record says JSON only.
- **The prune test pins fingerprints, not commits (T9.3).** A code fingerprint is the stronger identity, and
  `scripts/paper/verify_stamp.py` ties the record's to b730a05 (verified).
- **The prune record's stamp and its equivalence commit (T9.3).** The record is stamped b730a05; the equivalence ran
  at c64ea7a; between them src changed only in comments and docstrings. One sentence at its next render (section 1).
- **make_figures nits (T8).** The panel-letter offsets repeat, two FS_MAX comments lag the caps, and no test calls
  the changed figure functions. Figures are regenerated and looked at at the freeze; plan 29 owns the look test.

Closed: T2's M9 (Q5 applied per entry) is applied per stem since NA fix F7 (24676d3).
