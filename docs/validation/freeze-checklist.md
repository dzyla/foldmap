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
| `random100-fourth/` | not drawn | the paper's headline count | `fourth_sample.py draw --release TAG`, then one `random100.py run --views all`, in the release checkout (the protocol) |
| `sensitivity/` | runs at 5f09684's code | the NA fixes, the prune and the rest changed what it measures; the default baseline moved 80 to 82 (its "Code since this run" section, final review F1) | rerun stages B, V and X at the release commit (about 3.9 h on zl), then `sensitivity.py report` |
| `readability-floor/` | e6fcf36 | its census and measures predate later code; its generator now says R1's route-budget reason fires on 7NZ2's domains view (F9) | `scripts/readability_floor.py` |
| `prune/` | runs at 46f5ebf, record b730a05 | its next render adds one sentence: the record's stamp is b730a05 while the equivalence runs were at c64ea7a, and `git diff c64ea7a b730a05 -- src` is comments and docstrings only (parked T9.3) | `scripts/prune.py` report at the release |
| `metamorphic/` | this round's re-render | stamped at a pre-release commit | `scripts/metamorphic_record.py` |
| `na-fixes/` | 6ca884c | stamped at a pre-release commit; `na_runs:_JOINT_NEAR` was pruned after it (the prune measured the NA set identical) | `scripts/na_fixes_record.py` |
| `domain-uncertainty/` | 42e3fe8 | stamped at a pre-release commit | `scripts/domain_uncertainty.py` |
| `homologues/topology-transfer/` | cdd29dd | H1-H5 not remeasured since; low risk (domain-name placement and NA joints only changed) | `scripts/topology_transfer.py`, with the other homologue campaigns |
| `random100/`, `random100-heldout/` | 7696947 | "on this release" in the paper is false until they are rerun (claims R8, R9) | `random100.py run` at the tag, not for the headline |
| `random100-third/` | 7696947 | kept as it is: a historical prospective evaluation at 7696947 | none |
| `na/` (and the M2 NA rerun) | older code | the NA campaign the paper quotes | `M2_RERUN=1 scripts/paper/run_all.sh TAG HARNESS` |
| `inventory/` and its freeze baseline | b730a05 | the table must match the release's code (`tests/test_inventory.py` fails otherwise) | `scripts/inventory.py`; `scripts/freeze_baseline.py` only shrinks |
| the paper's numbers, tables and figures | 838c579 and later | everything quoted | `scripts/paper/run_all.sh TAG HARNESS` (HARNESS: the commit adding the fourth sample), `figures.sh`; inputs outside git rebuilt first (`scripts/paper/README.md`) |

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
