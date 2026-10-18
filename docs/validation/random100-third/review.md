# Random-100 third sample: review

**83 of 100 pass.** This is the count for structures foldmap was not tuned on: the sample was drawn after plan 15's
fixes, with none of the tuned or the second sample's ids, and run once with every view and no fix (2026-10-18, at
7696947, zyla-lab, `--jobs 10`, 1329 s wall), judged by the same criteria and limits as the other two samples. In the
same run of the same code the tuned sample (`docs/validation/random100/`) passes 93 of 100 and the second sample
(`docs/validation/random100-heldout/`, which plan 15 fixed against) 95 of 100. The run itself recorded its commit as
"unknown" (zyla-lab's synced copy of the tree has no .git); 7696947 was stamped into the reports afterwards from that
synced tree, and the reports were judged again from the measurements, with the same counts.

| stratum | entries | pass | | stratum | entries | pass |
|---|---|---|---|---|---|---|
| S1 X-ray, one chain | 12 | 12 | | S7 protein-DNA | 7 | 3 |
| S2 X-ray homo-oligomer | 12 | 10 | | S8 protein-RNA | 7 | 2 |
| S3 X-ray hetero-complex | 10 | 9 | | S9 RNA only | 3 | 3 |
| S4 cryo-EM | 15 | 13 | | S10 glycoproteins | 5 | 5 |
| S5 NMR | 8 | 8 | | S11 one protein, two states | 10 | 8 |
| S6 membrane | 8 | 7 | | S12 AlphaFold models | 3 | 3 |

Entries failing each criterion: A2b 1, A3 4, A4 2, A6 6, A15 3, A16 1, A17 7, A18 1, A19 2. A1, A2, A5, A7-A14 and
A20 pass on every entry that has something to judge (A1: no crash, no traceback, no Python warning in any view of any
entry). One defect can fail several criteria (a domain name on a loop fails A6 and A16).

The 17 failures are listed below by check, as measured. **None is investigated or fixed here**; for each, the cause
is not investigated unless a line says otherwise (a resemblance to a known cause is noted, not established).

## By check

### A17, nucleic-acid metrics (7 entries)

| id | stratum, bin | what failed |
|---|---|---|
| 8VBG | S7 M (cryo-EM, 3 chains) | contact_auc 0.505 (limit 0.6); stems_nearer_frac 0.5 (limit 1.0) |
| 6VTX | S7 S (X-ray, 108 residues) | stems_nearer_frac 0.0 (limit 1.0) |
| 6ZY7 | S7 L (cryo-EM, 6 chains) | stem_recall 0.667 (limit 0.85) |
| 8XJV | S7 XL (cryo-EM, 110 chains, 18,408 residues) | pair_f1 0.859; stem_recall 0.644; stems_nearer_frac 0.969; fallbacks 10; na_crossings 10; bead_overlaps 125 |
| 3EGZ | S8 S (X-ray, 2 chains) | pair_f1 0.941; spearman 0.347; procrustes 0.982; angle_median 22.0; angle_max 125.7; mirrored 1 |
| 4AS1 | S8 M (X-ray, 2 chains) | na_crossings 1 |
| 6QIK | S8 XL (cryo-EM, 47 chains) | fallbacks 1; na_crossings 1; bead_overlaps 1 |

stems_nearer_frac below 1 on a small protein-DNA complex is the measure 9YLM fails on the tuned sample (ruled inherent
to the measure there); 3EGZ's angle and mirror values are of the kind 30QM shows in the NA campaign (an unrolled
frame measured against the plain projection). Neither is established for these entries.

### A6, text collisions (6 entries)

| id | views | collisions |
|---|---|---|
| 9NMN | default, deposited, beta | label:A:4963-4967 × loop:A:4908-4954>A:4963-4967 (helix α153's name on the loop into it; also A4) |
| 9FZL | domains | 5: domain names D1 · 3A, D2 · 3B, D1 · 3a ... on loops (also A3, A4, A16) |
| 8XJV | default, deposited, beta, focus, domains | 4: na:stem:6 first-pair number × na:run:128; na:gap:140 × stub:C:Ao; terminus:N:P × na:run:27 |
| 3G6E | default, deposited, beta, focus | na:tie:168-219 label × na:tie:164-171 (a tie's label on another tie) |
| 6QIK | default, deposited, beta, focus | 3: tie labels 42-92 and 190-221 on other ties; terminus:N:i × a tie mark |
| 7F5S | domains | 2: tie label 2575-2757 × tie 2738-2742; terminus:N:Lr × a tie mark |

Tie labels on other ties (3G6E, 6QIK, 7F5S) are a collision the tuned and second samples do not show.

### A3, route fallbacks (4 entries)

| id | views | fallbacks |
|---|---|---|
| 8XJV | default, deposited, beta, focus, domains | strands Au:348-503, Au:536-615, Au:886-1050, Au:1271-1381 (nucleic-acid strands drawn as plain curves) |
| 6QIK | default, deposited, beta, focus, domains | strand y:8-9 |
| 7F5S | domains | strand L5:694-695 |
| 9FZL | domains | loops 4E:88-89>4E:119-131, 4K:15-17>4K:33-35 and their copies 4e, 4k |

9FZL's domains view, as measured (98 chains): the page is held at the 200 in cap (200 x 11 in; its notes: "some text
is under 1 pt and drawn at 1 pt, larger than the room kept for it", the figure's text 1.6 pt); its route record is
878 loops, 5 fallback curves, each crossing an element (the four loops above and 4x:94-101>4x:104-105), 9 second
searches run and 30 refused for the budget, 7,976,402 of the page's 8,000,000 second-search nodes spent.

**After the run (2026-10-18, plan 15's final fix round): 9FZL has now been looked at.** Its route record is the
profile of the regression found on the second sample's 7NZ2 (`docs/validation/random100-heldout/review.md`): a coarse
page's fine windows drew on the page's own second-search budget, already spent. With the windows given a pool of
their own (`route._MAX_NODES_WINDOWS`, 8,000,000 nodes), its domains view drawn again (zyla-lab) has 2 fallback curves
(4K:15-17>4K:33-35 and its copy 4k; 4E, 4e and 4x now route), 14 searches refused, the window pool spent; A4 is
unchanged (the same 11 loops over a label), A6 has 6 collisions (one new: 3A:274-282>3A:291-301 across the N end of
4G). 9FZL still fails; the 83 of 100 above stays as measured at 7696947.

### A15, copies drawn alike (3 entries)

| id | view | worst tie |
|---|---|---|
| 8PIN | focus | B~C 18.88 (D~A 16.83), all 30 elements matched |
| 9PDH | default | B~C 9.68 (33 of 34 matched) |
| 4YJ5 | focus | C~B 4.19, D~A 4.17, all 37 matched |

No sheets differ on these ties (they are judged).

### A4, a loop over a label (2 entries)

- 9NMN (default, deposited, beta): A:4908-4954>A:4963-4967 (α152 to α153) over α153's name (the same defect as its A6).
- 9FZL (domains): 2M:113-132>2M:147-154, 3A:274-282>3A:291-301, 3B:520-525>3B:533-541, 3C:234-267>3C:270-272.

### A16, domains (1 entry)

- 9FZL: 4 domain-name collisions with loops (D1 · 3A, D2 · 3B, D1 · 3a ...), the same text as its A6. A domains-view
  name on a loop is what Task 2 of plan 15 fixed for 7BKC and 6ZKT and what 7NZ2 shows again since Task 7
  (`docs/validation/random100-heldout/review.md`); whether 9FZL shares that cause is not investigated.

### A19, two states named alike (1 pair, 2 entries)

- 8OJZ/8QX1 (P1, apo/holo, CA RMSD 3.4 Å): names compared on all 48 pairs, 14 differ. 8QX1 draws a helix,
  A1:9-19 (αNa), that 8OJZ does not (its first is A1:32-49), so 8QX1's N-helices are named one letter later (8OJZ's
  αNa is 8QX1's αNb, and so on), and later names shift too (α1~α1″, α1″~α1‴). The same kind of inconsistency as
  9YKJ/10UJ on the tuned sample: one protein named two ways in its two states.

### A18, membrane (1 entry)

- 6DQN (S6 XL, cryo-EM, 4 chains, 10,684 residues): no band drawn.

### A2b, the app's hosted limits (1 entry)

- 8XJV: the app's job took longer than 90 s (the hosted limit) and was stopped, with foldmap's message (draw part of it,
  or locally). Its default view took 129 s in this run (A2's limit for XL entries is 180 s, so A2 passes), so the app
  limit is not met by load alone; the run drew 10 entries at once.

## Science checks and S5 (reported, not gated)

S1-S4 and S6 find nothing in any view. S5 (the automatic domains against CATH, not counted in the 83): 109 chains
scored in 35 entries. One-domain chains: 55 pass, 7 fail (the automatic domains split a chain CATH leaves whole);
chains with two or more CATH domains: 12 pass, 35 fail (mean boundary F1 0.53, recall 0.75). Entry level: S5 passes
on 23 entries and fails on 12.

## Time

Wall 1329 s for the sample at `--jobs 10`. The slowest entry is 8XJV (1092 s over its ten views; interactive page
138 s, default 129 s), then 9FZL (706 s; domains view 234 s) and 7F5S (483 s). A2 passes on every entry; the
slowest default view is 8XJV's, 129 s against 180 s.
