# Automatic domains

## Scope

- `features.auto_domains`: each chain's domains under `--domains auto`, as groups of elements.
- `features.residue_domains`: the residue-level split: a part is cut at its lowest normalized cut (`features._best_cut`,
  `features._lowest_cut`) while that cut is below `features.DOMAIN_NCUT`.
- `features._chain_domains`: the residue domains turned into element domains: each element in the domain holding most
  of its residues, stray elements taken home (`features._no_strays`), small domains folded into their best-connected
  neighbour, a βα row kept whole.
- `features.split_confidence`, `features.domain_splits`, `domain_confidence`: each kept split said with its
  confidence.
- `domain_score`: automatic domains scored against a reference chopping (`S5`). `cath.chopping`: CATH's chopping from
  the PDBe API (`--domains cath`).
- `layout._part_domains`, `layout._place_domain_names`: the domains view, one panel per domain.

## Assumes

- Two residues are in contact when their CA atoms are within `features.DOMAIN_CONTACT` (12 Å) and at least
  `features.DOMAIN_SEPARATION` (6) residues apart in sequence.
- A part splits when its best cut's normalized cut is below `features.DOMAIN_NCUT` (0.37). The value was chosen on the
  dom_cal chains.
- A domain needs `features.DOMAIN_MIN_RESIDUES` (40) residues, each piece of a cut `features.DOMAIN_MIN_SEGMENT` (20),
  and `features.DOMAIN_MIN_SSE_RESIDUES` (30) residues in helices and strands. One cut takes a domain of at most two
  runs.
- A domain lies in one chain. Copies of a protein are split alike.
- A βα row (TIM barrel, Rossmann sheet) is never split.

## Guarantees

- On the 33-chain CATH calibration set the means meet the ruled targets (overlap 0.80, F1 0.70, recall 0.85;
  `domain_score.TARGETS`): `tests/test_domain_score.py::test_auto_domains_meet_the_ruled_targets_on_the_cath_set`.
  This is in-sample: the parameters were chosen on these chains. Per chain, 67 of the third sample's 109 scored chains
  meet the targets (`docs/validation/domain-uncertainty/report.md`, t = 0).
- Exact under rigid motion, renaming, shuffling and the PDB writer:
  `tests/test_metamorphic.py::test_automatic_domains_are_exact`.
- Copies of a tetramer split alike, each domain in one chain:
  `tests/test_domains.py::test_auto_domains_of_a_tetramer_lie_each_in_one_chain_and_copies_split_alike`.
- A βα row lies in one domain panel: `tests/test_domains.py::test_a_ba_row_lies_in_one_domain_panel`.
- The split confidence is 1 - (normalized cut / `features.DOMAIN_NCUT`), the margin below the limit:
  `tests/test_domain_confidence.py::test_confidence_is_the_margin_below_the_limit`. Saying it never changes the
  figure: `tests/test_domain_confidence.py::test_the_figure_is_the_same_with_or_without_the_record`.
- What the confidence is worth against CATH (`docs/validation/domain-uncertainty/report.md`): AUC 0.75 on the two
  regression sets (44 splits), 0.61 on dom_cal where the limit was chosen, 0.87 on the third sample's 99 splits (0.74
  counting each entry once). n is small; the Wilson intervals are in the record. It is a margin, not a probability.

## Fails

- On the third sample, 47 chains have two or more CATH domains: 12 meet the targets and 35 fail, in 8 entries. Seven
  single-domain chains are split (1RKA A, 3WAU A and B, 1XZ4 B and D, 3W96 D and H)
  (`docs/validation/domain-uncertainty/report.md`).
- Boundary missed (15 chains): a CATH domain packed too tightly to cut (its own normalized cut at 0.37 or above: 2AIQ A,
  3ABR A and C, 3EOM A, B and D, 3G6E W), or a residue cut made and then undone. The undoing step: the `fold` step of
  `features._chain_domains`, under `features.DOMAIN_MIN_SSE_RESIDUES`, in 7 chains (3G6E A, L and P; 4YJ5 A, B, C
  and D); `features._no_strays` in 1 (3G6E B).
- Boundary moved (5 chains): the lowest cut lies elsewhere than CATH's (1GQ2 A, C and L; 1NID A; 3EOM C).
- Extra split (15 chains): a cut below the limit inside one CATH domain, a lobe CATH keeps whole (1GQ2's other 13
  chains; 4ZH4 F and L).
- Keeping weakly supported splits undivided: every threshold above 0 lowers some dom_cal chain (record, section 3).
  Only the third sample would gain: at t = 0.30, 67 to 83 chains meet the targets (17 newly, 12 of them 1GQ2's
  copies; 4GWY A no longer does; `docs/validation/domain-uncertainty/report-chains.json`). It cannot choose t. No
  threshold is drawn from the confidence.
- Every domain constant is sensitive (`features.DOMAIN_NCUT`, `features.DOMAIN_CONTACT`,
  `features.DOMAIN_SEPARATION`, `features.DOMAIN_MIN_RESIDUES`, `features.DOMAIN_MIN_SEGMENT`;
  `docs/validation/sensitivity/report.md`). `features.DOMAIN_NCUT` at x1.25 takes the third sample's `S5` failures from
  12 to 19 entries. `panels.CLOSE` (the gaps a domain's panel fills, set on 2LZM) is fragile: at x0.75 6QIK's
  domains view fails `A16`, which no entry of the domains stage failed unperturbed (run at v0.2.0).
- Domain names on loops in the domains view: 9FZL (`A16`, `A6`, third sample); 7NZ2 (second sample). `layout:_NAME_OPEN`,
  set for 7NZ2, was pruned with every key equal (`docs/validation/prune/report.md`).

## Seen by

- `S5`: automatic domains against CATH's chopping, where the committed reference holds the entry (campaign only).
  Elsewhere `foldmap check` reports it not measured, with each split's confidence as its measures.
- `A16`: domains in one chain, alike across copies, not tiny, panels apart, names clear.
- Not seen by any check: a wrong split on an entry with no CATH reference. The confidence is said in the notes; nothing
  judges it.
- Not seen by any check: the residue cuts the element rules undo. The record's trace shows them; `foldmap check` does
  not.

## Evidence

- [domain-uncertainty/report.md](../validation/domain-uncertainty/report.md): the failure classes, the mechanisms per
  chain, the confidence against CATH, and the undivided-threshold sweep.
- [sensitivity/report.md](../validation/sensitivity/report.md): the domain constants, all sensitive; `panels:CLOSE` fragile.
- [prune/report.md](../validation/prune/report.md): `layout:_NAME_OPEN` pruned; `layout:_PART_BLOWUP` and
  `panels:CLOSE` kept.
- [random100-third/review.md](../validation/random100-third/review.md): `S5` and `A16` on the third sample.
