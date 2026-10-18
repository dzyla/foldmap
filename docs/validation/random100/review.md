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

## The final run (final review fix F1): 92 of 100, and the eight that fail

`report.md`, `report.json` and `run.log` here are the final run: every view of all 100 entries at 5650dbd, on
zyla-lab, from an empty download cache. Every downloaded file's SHA-256 matches the previous run's (Task 13, 96 of
100), and every view's SVG and layout signature are the same as that run's but one, 1RM5's focus view
(below). The count fell from 96 because the checks now measure what they claim (final review B):

- **Nothing to judge is never a pass.** Each criterion reads pass, fail, not judged or does not apply; `report.md`
  counts all four and lists every case not judged. The headline counts that were vacuous before: A15 37 judged and
  passed (1 not judged, 62 with no copies), not 100; A20 9 (52 with no tether to another chain, 39 with no ligand),
  not 100; A16 47 passed and 1 failed of the 62 domains views (14 drew no domain); A10 96 measured (mkdssp 3.1.2
  agrees on 268,834 of 268,835 residues), the 4 RNA-only entries do not apply. Science: S1 90 pass, 6 not judged
  (no sheet drawn), 4 do not apply (no protein); S2 13 pass, 83 not judged (no closed barrel); S4 2 pass, 94 not
  judged (no residue mark); S6 37 pass, 59 not judged (no accessory strand).
- **A19 compares every pair's names**, each generic number less the unpaired elements of its kind before it (it
  compared none in 9YKJ/10UJ and 6YPX/6YQD before).
- **A9 counts the file's chains with gemmi**, apart from the assembly built. 1RM5's assembly 1 is built from chains A
  and B; chain O is assembly 2's. Drawing assembly 1 is right; the note was wrong (it said nothing was left out). It
  now reads "assembly 1, which holds no copy of chain O of the file: 4 of 5 chains; --assembly asu --focus none draws
  all of them", and A9 passes on it. Because the default now draws part of the file, the focus view draws the whole
  file (`--focus none`): the one figure that changed.
- **A17 is measured on the view that draws the nucleic acid** when the default view draws none (8QEM).
- **A6 reads chain-end stubs and domain leaders** as lines (their own letter, tail length or name exempt).

The eight that fail:

- **7Y7C and 8VWQ** (ribosomes, about 10,700 residues in 55 chains). Unchanged: residuals carried, not ruled inherent.
  In the `--domains` view a strand and a gap fall back to a plain curve and a stem's first-pair number touches a run
  (7Y7C); a ligand label sits on stem 114 in four views and a tie label on stem 209 (8VWQ); a free end crossing a gap
  path, and one bead overlap each (A3, A6, A17). Each has a mechanism (routing through the margins, hard tails for
  gaps), left for a later plan.
- **9YLM** (three zinc fingers along one 12-bp duplex), A17 `stems_nearer_frac` 0.0. Unchanged; ruled inherent to the
  measure (R1): the recognition helices miss by 0.09 page units, because the protein's own layout puts the fingers'
  small β-hairpins between them and the DNA, and every variant that passes it breaks 1IHF, 1ZBI, 3KTQ or 4OO8.
- **8BW5** (a DNA aptamer beside a G-quadruplex), A17 `pair_f1` 0.8 and `stem_recall` 0.5. Unchanged; ruled
  inherent: RNApolis' two extra canonical pairs (F5·F37, F10·F33) are isolated A·T pairs with stretched bonds
  (C1'-C1' 10.5 and 12.1 Å), which foldmap takes only when stacked in register on a strict pair, so it draws 4 of 6
  pairs (`run1/review.md`, batch B).
- **9YKJ and 10UJ** (pair P2, one protein as monomer and oligomer), A19: 27 of 43 paired elements are named
  differently. 9YKJ's chain A is drawn with a βα row, so `--names beta` names its helices after the strand before
  them (αN, αN′, α1, α3′, α10a …); 10UJ's chain A has no βα row, so its helices keep generic numbers (α2 … α34). The
  same protein is named two ways in its two states: a real inconsistency, newly measured, not fixed here.
- **7JU4** (`--domains auto` view), A6 and A16: the domain names "D2 · 4" and "D3 · A" sit on the chain-end stubs
  C:4 and C:A. Newly measured (A6 did not read stubs before); a real collision, not fixed here. A16's detail reads
  "domains 2 ()": its two collisions are the A6 ones.
- **8QEM** (a C12 ring of 24 protein chains around a 63-bp DNA duplex), A17, first measured. Its default view draws
  one protomer (A, N) of the ring, as `--focus auto` does for an assembly this large: the DNA touches seven of the
  twelve protomers with only 1 or 2 residues each, so no protomer holds it (`NA_HOLD_MIN` 3) and the protomer drawn
  touches none of it. That choice is right by the rule, and the note says 6 of 26 chains were drawn, so A17 is
  measured on the focus view (`--focus none`), which draws the DNA. There the duplex's stems stand off their
  projected 3D axes: `angle_median` 3.2°, `angle_max` 4.6°, `joint_error_max` 3.4° (limits 0.05°, 0 on a plain
  frame). Not investigated here.

### The two-state pairs, re-typed

sample.json keeps the type each pair was drawn as; `report.md` ("Two-state pairs") gives the type its files read
now, with the additives foldmap hides left out of the ligands and the CA RMSD over residues aligned by sequence
(final review B, Important 4 and Minor 9). The sample is not redrawn. Two pairs change:

- **P3, 6YPX/6YQD, drawn as apo/holo, is neither.** 6YPX's only ligand is glycerol (GOL), an additive foldmap hides;
  without it neither entry has a ligand touching the protein. Both are 2-chain assemblies, CA RMSD 2.54 Å: one state,
  one of them soaked in glycerol.
- **P5, 6FYF/6FY9, drawn as open/closed, is neither.** Its 3.78 Å RMSD paired residues by author number; 6FY9
  numbers the protein one behind 6FYF (A:8-17 against A:7-16). Aligned by sequence the RMSD is 0.63 Å: one state,
  both with FAD.

P1 (1RM5/2HKI, apo/holo), P2 (9YKJ/10UJ) and P4 (8P3N/8OYQ, monomer/oligomer) keep their types. A19 judges every
pair as two entries of one protein, whatever their states.
