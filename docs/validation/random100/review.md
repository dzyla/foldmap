# Random-100 review

Causes (C1…Cn) are written here by Task 12, each with the entries it explains and the test that pins its fix.
Measure notes (M1…) record where a metric, not the drawing, is what falls short.

## The plan 15 run: 93 of 100

Run 2026-10-18 at 7696947 on zyla-lab, every view, `--jobs 10`, 864 s wall, after plan 15 fixed the second sample's
failures (`docs/superpowers/plans/2026-10-17-15-heldout-failures.md`); `report.md`, `report.json` and `run.log` are
that run. The layout's constants were tuned on this sample, so 93 measures fit. The count for structures foldmap was
not tuned on is the third sample's, 83 of 100, in the same run of the same code
([docs/validation/random100-third/review.md](../random100-third/review.md)); the second sample, which plan 15 fixed
against, passes 95 of 100.

Against the final run below (92 of 100 at 44a1c1b), field by field in `report.json`:

- **One entry changed verdict: 7JU4 now passes** (its two domain names on chain-end stubs are gone: A6 and A16
  2 → 0; leaders 14 → 16). No entry that passed fails. Not attributed to one plan 15 task.
- **The seven that fail**, all with the same checks except 8VWQ: 7Y7C (A3, A6, A17: na_crossings 4 → 1,
  bead_overlaps 1 → 0), 8VWQ (A3, A17; its A6 is gone: the LYS:5101 ligand label off stem 114 in four views and the
  tie label off stem 209 in the domains view; bead_overlaps 1 → 0), 9YLM and 8BW5 (A17, as ruled), 9YKJ/10UJ (A19,
  27 of 43 names differ, as before) and 8QEM (A17 on the focus view, the same values). The causes below stand.
- **A10** (DSSP agreement, after plan 15 Task 1's one-conformer read): identical detail on all 100 entries; the
  268,834 of 268,835 residues stand.
- **A15 "does not apply" → "not judged"** on 10 entries (Task 5: copies drawn untied are not judged, with the reason):
  5SB5, 6DPV, 5JB1, 7JU4, 8ESZ, 7ARD, 7SMR, 6JO6, 5LDW, 8OZJ. A15 still judges 37 entries and passes all 37.
- **Domains view**: A4, A6 and A16 change only on 7JU4 and 8VWQ (above); domain-name counts are unchanged on every
  entry; leader counts moved on 1A0Z (4 → 2), 7ARD (8 → 7), 9GEP (5 → 4), 8VWQ (10 → 9) and 7JU4 (14 → 16).
- **A20 tethers** (default view): lengths moved on 7OKA, 6JLN (489 → 493 tethers), 7Y7C and 8VWQ (on the two
  ribosomes the magnesium tethers grew from 0.7-1.05 to 1.05-2.8, consistent with Task 3 moving protein-held ligand
  marks off the nucleic acid; the other two are not attributed); A20's verdict is unchanged on all four (pass on
  6JLN and 7Y7C, not judged on 7OKA and 8VWQ).
- **Nucleic-acid metrics**: only 7Y7C (na_crossings 4 → 1, na_run_crossings 724 → 723, bead_overlaps 1 → 0) and 8VWQ
  (na_run_crossings 637 → 639, bead_overlaps 1 → 0) moved.
- **sheets_differ**: 40 over all views, as before, entry by entry.
- **Route fallbacks**: none on any entry, before or after. Second-search nodes moved on 7JU4 domains (+681),
  9P1X focus (694,501 → 816,673), 5LDW membrane (0 → 415,299), 7Y7C and 8VWQ (Task 7 now counts window and
  slab-closed searches in these totals).
- **Time** (10 entries at once): A2 passes on every entry, the slowest default view being 7Y7C's, 66.9 s (limit 180 s
  for XL); the slowest view is 7QUP focus, 148.5 s.

S5 (automatic domains against CATH, reported, not gated, not in the count): 124 chains scored in 30 entries.
One-domain chains: 57 pass, 8 fail (split by the automatic domains); chains with two or more CATH domains: 14 pass,
45 fail.

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
zyla-lab, from an empty download cache, then (F1 fix round 1, at 44a1c1b) every view of the eight entries whose
result the fixes could change (8QEM, 1RM5, 7JU4, and the entries with no protein element 8Q4O, 6B5W, 8SA2, 3OWW,
8XZE) drawn again, and every entry judged again from what was measured (`random100.py run --ids`, whose report
judges all entries). Their figures are the same as in the full run; only notes and judgements changed. Every downloaded file's SHA-256 matches the previous run's (Task 13, 96 of
100), and every view's SVG and layout signature are the same as that run's but one, 1RM5's focus view
(below). The count fell from 96 because the checks now measure what they claim (final review B):

- **Nothing to judge is never a pass.** Each criterion reads pass, fail, not judged or does not apply; `report.md`
  counts all four and lists every case not judged. The headline counts that were vacuous before: A15 37 judged and
  passed (1 not judged, 62 with no copies), not 100; A20 9 (52 with no tether to another chain, 39 with no ligand),
  not 100; A11 (N->C and sheet order) 95, not 100 (the 4 entries with no protein do not apply, and 6B5W, a
  protein with no element drawn, is not judged; A8 likewise not judged for 6B5W); A16 47 passed and 1 failed of the 62 domains views (14 drew no domain); A10 96 measured (mkdssp 3.1.2
  agrees on 268,834 of 268,835 residues), the 4 RNA-only entries do not apply. Science: S1 90 pass, 6 not judged
  (no sheet drawn), 4 do not apply (no protein); S2 13 pass, 83 not judged (no closed barrel); S4 2 pass, 94 not
  judged (no residue mark); S6 37 pass, 59 not judged (no accessory strand).
- **A19 compares every pair's names**, each generic number less the unpaired elements of its kind before it (it
  compared none in 9YKJ/10UJ and 6YPX/6YQD before).
- **A9 counts the file's chains with gemmi**, apart from the assembly built. 1RM5's assembly 1 is built from chains A
  and B; chain O is assembly 2's. Drawing assembly 1 is right; the note was wrong (it said nothing was left out). It
  now reads "drew assembly 1, which holds no copy of chain O of the file: 2 of 3 chains; --assembly asu --focus none
  draws all of them" (the count is of the file's chains, as `--assembly asu` draws them), and A9 passes on it.
  Because the default now draws part of the file, the campaign's focus view is `--focus none`: the one figure that
  changed. That view still builds assembly 1 (A1, B1, A2, B2), so it repeats the default figure; chain O is drawn in
  no view of the campaign.
- **A subunit focus names the nucleic acid it leaves out**, and A9 counts a nucleic-acid chain outside the focus
  unless the note names it. 8QEM's default note read "drew one of 12 subunits (C12 symmetry): 2 of 26 chains", with
  its 63-bp duplex hidden in the 26; it now reads "drew one of 12 subunits (C12 symmetry), without the DNA Z, Y: 2 of
  26 chains; --focus none draws all of them".
- **A11 and A8 with nothing drawn** are not judged (a protein with no element: 6B5W) or do not apply (A11 with no
  protein), never passed.
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
  C:4 and C:A. Newly measured (A6 did not read stubs before); a real collision, not fixed here. It is one defect
  counted by two criteria: A16's "names clear" reads the same text collisions.
- **8QEM** (a C12 ring of 24 protein chains around a 63-bp DNA duplex), A17, first measured. Its default view draws
  one protomer (A, N) of the ring, as `--focus auto` does for an assembly this large: the DNA touches seven of the
  twelve protomers with only 1 or 2 residues each, so no protomer holds it (`NA_HOLD_MIN` 3) and the protomer drawn
  touches none of it. That choice is right by the rule, and the note now names the duplex it leaves out (above),
  so A17 is measured on the focus view (`--focus none`), the only view that draws the DNA. There the duplex's stems
  stand off their projected 3D axes: `angle_median` 3.2°, `angle_max` 4.6°, `joint_error_max` 3.4°. The 0.05°
  limits (`nacheck.THRESHOLDS`) are a regression guard for plain frames, where these angles are 0 by construction
  (one rounding step); a nonzero value means the whole ring is laid out on an unrolled frame, where they measure the
  unrolling against the plain projection of the same orientation, as ruled for 30QM in the NA campaign
  (`docs/validation/na/review.md`, M2 and the C4 ruling). 3.2°, 4.6° and 3.4° would pass the earlier 5° limits. It
  fails by the plan's limits; it is not a sign that the DNA is drawn wrong.

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

### Counting failures per criterion

An entry's count is exact; a per-criterion count of failures is not a count of defects. One defect can fail
several criteria: 7JU4's two domain names on stubs fail A6 and A16; a view that times out or is not measured fails
every criterion read from it (fail closed, as 8D3D's focus view in the held-out sample fails A1, A5-A8 and A15).
