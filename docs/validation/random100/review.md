# Random-100 review

Causes (C1…Cn) are written here by Task 12, each with the entries it explains and the test that pins its fix.
Measure notes (M1…) record where a metric, not the drawing, is what falls short.

## Measure notes

**M1. 1ZBI `contact_auc` 0.65 (limit 0.60, higher is better): the metric, not the drawing** (plan 10 Task 7(h),
carried item K18).
- 1ZBI is RNase H bound to one 12-bp RNA/DNA duplex (one stem, `C:1-12/D:1-12`). Ten protein elements touch
  it, five of each of the two protein copies, which hold the duplex from both sides in 3D.
- The figure draws the two copies alike, side by side, with the duplex between them. Each copy's elements that
  touch the far face of the duplex in 3D are drawn on its own side of the page, away from the duplex: `A:105-122`
  at 7.75 page units, `B:178-180` at 8.48. Seven non-touching elements of the compact domains lie nearer
  (`A:79-87` 4.75, `A:93-103` 5.85, `B:93-103` 3.17, ...).
- A page cannot put both faces of one duplex next to both copies. Every stem the protein touches still has its
  touching elements nearer than the median (`stems_nearer_frac` 1.0). No drawing defect; the limit is unchanged.

## The final run (plan 10 Task 13): 96 of 100, and the four that fail

`report.md`, `report.json` and `run.log` here are the final run: every view of all 100 entries at fce7d7c, on zyla-lab,
from an empty download cache. Every downloaded file's SHA-256 matches the previous run's, and every criterion,
science check and view result is the same as that run's, apart from A10, which is now measured: mkdssp 3.1.2 agrees
on 268,834 of 268,835 residues over the 96 entries with protein (one residue of 8K8E, 0.9992 of 1318; the four
RNA-only entries have no protein to compare). The causes of the earlier failures, C1…C34, and how each was fixed, are
in `run1/review.md`. Four entries still fail, each with a ruling:

- **7Y7C and 8VWQ** (ribosomes, about 10,700 residues in 55 chains). Fixed down to residuals that are carried, not
  ruled inherent: in the `--domains` view a strand and a gap fall back to a plain curve and a stem's first-pair number
  touches a run (7Y7C); a ligand label sits on stem 114 in four views and a tie label on stem 209 (8VWQ); a free end
  crossing a gap path, and one bead overlap each (A3, A6, A17). Each has a mechanism (routing through the margins,
  hard tails for gaps), left for a later plan.
- **9YLM** (three zinc fingers along one 12-bp duplex), A17 `stems_nearer_frac` 0.0. Ruled inherent to the measure:
  the recognition helices miss by 0.09 page units, because the protein's own layout puts the fingers' small
  β-hairpins between them and the DNA, and every variant that passes it breaks 1IHF, 1ZBI, 3KTQ or 4OO8. A
  placement inside a domain that keeps the order is carried as a follow-up.
- **8BW5** (a DNA aptamer beside a G-quadruplex), A17 `pair_f1` 0.8 and `stem_recall` 0.5. Ruled inherent:
  RNApolis' two extra canonical pairs (F5·F37, F10·F33) are isolated A·T pairs with stretched bonds (C1'-C1' 10.5 and
  12.1 Å), which foldmap takes only when stacked in register on a strict pair, so it draws 4 of 6 pairs
  (`run1/review.md`, batch B).
