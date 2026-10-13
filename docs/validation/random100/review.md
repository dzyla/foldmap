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
