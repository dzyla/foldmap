# Nucleic-acid fixes, measured

Plan 28 (ruling Q37). Written by: foldmap 0.1.0, git 3fb5df1ecaa9 (asserted by runner, dirty False), code fingerprint be5353d8a5d01df6, scripts/na_fixes_record.py sha256 8757ecc4696c; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Each fix measured before and after on the 7 failing cases of the third random sample (8VBG, 6VTX, 6ZY7, 3EGZ, 4AS1, 6QIK, 8XJV) and 39 nucleic-acid controls that passed, in the default view, with scripts/random100.py's entry runs (A17 as random100 takes it). **These are regression sets, not unseen ones**: the failures were diagnosed on these entries and the fixes made for them. The unseen measurement is the fourth sample at the frozen release. A fix is committed only when its target moves and no control gets worse on any judged A17 measure, fallback count, crossing count or bead overlap (ruling Q37).

Runs (the code each step measured; git is the commit checked out, the step's fix uncommitted on it, so the code fingerprint names the code: a committed fix's equals its commit's src):

- F2: code 9b4ffd133cae3036, git 85481f388d6a
- F3: code 18715a42781107ab, git 218f4bb26642
- F4: code 71ec92ce87407de7, git 218f4bb26642
- F5: code 403da61da25c40ed, git 3380435a1b82
- F6: code be5353d8a5d01df6, git 3380435a1b82
- F7: code 41c64420bd651819, git 5726aeca9fe8
- base: code 72fb15138b17e48a, git 5726aeca9fe8

## F7: committed in 24676d3

- Mechanism: M3: stems_nearer_frac scored where every drawn element touches the stem.
- Change: nacheck.contacts scores a stem in stems_nearer_frac only when an element does not touch it; such stems are counted apart (stems_nearer_not_applicable) and reported not applicable (ruling Q5).
- Bench: every entry identical (the bench's structures hold no case these fixes act on).
- Steps: base -> F7.

Target:

| entry | measure | before | after |
|---|---|---|---|
| 6VTX | stems_nearer_frac | 0 | None |
| 6VTX | stems_nearer_not_applicable | None | 1 |

Cases that moved:

- 6VTX: stems_nearer_frac 0 -> None; stems_nearer_not_applicable None -> 1 (worse: stems_nearer_frac)

Controls that moved:

- none

Controls worse: none

## F2: committed in 36adc3b

- Mechanism: M4: unrolled RNA placed by the protomer it faces, measured against one plain view.
- Change: a CylinderFrame unrolls each nucleic-acid chain one protomer holds (cli._na_holders: at least NA_HOLD_MIN residues touching it and more than every other protomer together) with that protomer, its stems tied to protomer 0's copy (layout._stem_ties); the shape measures take each copy in a plain frame at its own centre, Spearman and Procrustes per copy (the worst kept), the angles judged against half the turn between neighbouring copies.
- Bench: every entry identical (the bench's structures hold no case these fixes act on).
- Steps: F7 -> F2.

Target:

| entry | measure | before | after |
|---|---|---|---|
| 3EGZ | spearman | 0.347 | None |
| 3EGZ | procrustes | 0.982 | 0.572 |
| 3EGZ | angle_median | 22 | 10.3 |
| 3EGZ | angle_max | 126 | 10.4 |
| 3EGZ | mirrored | 1 | 0 |
| 3EGZ | na_run_crossings | 11 | 8 |

Cases that moved:

- 3EGZ: spearman 0.347 -> None; procrustes 0.982 -> 0.572; angle_median 22 -> 10.3; angle_max 126 -> 10.4; mirrored 1 -> 0; na_run_crossings 11 -> 8; ss_contacts_nearer 0.8 -> 0.95 (worse: spearman)

Controls that moved:

- none

Controls worse: none

## F4: committed in 3380435

- Mechanism: M5: a gap meets a free end's tip along the end's own outward heading.
- Change: the 5' free end after a gap turns toward the gap's start as well as its touching protein; the gap tries a free arrival along the chord, then the clear straight chord, before an unblocked route; a 3' end the protein touches keeps lying along it.
- Bench: every entry identical (the bench's structures hold no case these fixes act on).
- Steps: F2 -> F4.

Target:

| entry | measure | before | after |
|---|---|---|---|
| 4AS1 | na_crossings | 1 | 0 |
| 8XJV | na_crossings | 10 | 10 |

Cases that moved:

- 4AS1: na_crossings 1 -> 0
- 6QIK: na_run_crossings 406 -> 404; na_protein_crossings 362 -> 360
- 8XJV: bead_overlaps 125 -> 129; na_protein_crossings 731 -> 725 (worse: bead_overlaps)

Controls that moved:

- none

Controls worse: none

## F5: **not committed** (reverted)

- Mechanism: M6: a short duplex cut where projection, not 3D, makes the bend.
- Change: the layout (layout._kept_whole) keeps a duplex under 2 * SEGMENT_MIN pairs whole when its 3D bend is under 60 degrees but its pieces would be drawn turned more than 90 degrees (Nucleic.whole, read by stem_pieces).
- Not committed: it moved its target (6QIK fallbacks 1 -> 0, na_crossings 1 -> 0, bead_overlaps 1 -> 0) but made control 1RVA worse: its 10-bp EcoRV site (kinked 50 degrees in 3D, under WHOLE_BEND, 104 on the page) was drawn whole, its kink gone (bend_3d -102 -> None; angle_median, angle_max, joint_error_max 0 -> not measured); 6QIK itself fell to stems_nearer_frac 0.992 and gained 100 run x run and 96 NA x protein crossings. No 3D bend separates 6QIK (52.4) from 1RVA (50): its case stays a limitation. The change as measured: reverted/F5.patch.
- Steps: F4 -> F5.

Target:

| entry | measure | before | after |
|---|---|---|---|
| 6QIK | fallbacks | 1 | 0 |
| 6QIK | na_crossings | 1 | 0 |

Cases that moved:

- 6ZY7: procrustes 0.227 -> 0.217; contact_auc 0.823 -> 0.784; na_protein_crossings 1 -> 0; stems_drawn 4 -> 3; contact_precision_at_k 0.071 -> 0.083; bend_3d 150 -> 52.1; bend_drawn 150 -> 52.1 (worse: contact_auc)
- 6QIK: spearman 0.752 -> 0.737; procrustes 0.217 -> 0.269; contact_auc 0.948 -> 0.928; stems_nearer_frac 1 -> 0.992; fallbacks 1 -> 0; na_crossings 1 -> 0; bead_overlaps 1 -> 0; na_run_crossings 404 -> 504; na_protein_crossings 360 -> 456; stems_drawn 176 -> 175; contact_precision_at_k 0.102 -> 0.093; ss_contacts_nearer 0.998 -> 1; bend_3d -53.3 -> 57.8; bend_drawn -53.3 -> 57.7 (worse: spearman, procrustes, contact_auc, stems_nearer_frac, na_run_crossings, na_protein_crossings)

Controls that moved:

- 1RVA: angle_median 0 -> None; angle_max 0 -> None; joint_error_max 0 -> None; contact_auc 0.609 -> 0.749; stems_drawn 2 -> 1; contact_precision_at_k 0.25 -> 0.538; ss_contacts_nearer 1 -> 0.5; bend_3d -102 -> None; bend_drawn -102 -> None (worse: angle_median, angle_max, joint_error_max)

Controls worse: 1RVA (angle_median, angle_max, joint_error_max)

## F6: committed in 218f4bb

- Mechanism: M7: bead clearance best effort.
- Change: a run nucleotide _unbunch leaves short of a bead's width from another run's has the other run slide (na_labels._clear_junctions, _slide_from); routing runs for (m + 1) * max(NT_STEP, bead) changes nothing, as NT_STEP (0.55) is wider than a bead (0.45), so it is not made.
- Bench: every entry identical (the bench's structures hold no case these fixes act on).
- Steps: F4 -> F6.

Target:

| entry | measure | before | after |
|---|---|---|---|
| 6QIK | bead_overlaps | 1 | 0 |
| 8XJV | bead_overlaps | 129 | 127 |

Cases that moved:

- 6QIK: bead_overlaps 1 -> 0
- 8XJV: bead_overlaps 129 -> 127

Controls that moved:

- none

Controls worse: none

## F3: **not committed** (reverted)

- Mechanism: M2: a stem inside the protein's projection pushed out by _attract.
- Change: in a plain Frame, a stem whose projected centre lies inside the hull of its touching elements gets room before _attract (its overlapping items moved straight away from it) and a 10 times stronger contact pull.
- Not committed: its target barely moved (8VBG contact_auc 0.505 -> 0.521, stems_nearer_frac 0.5 unchanged: of its two hairpin pieces only #1 lies inside the hull of its touching elements) while ten controls got worse: the hull test also takes 1IHF's bent DNA (procrustes 0.274 -> 0.389, contact_auc 0.786 -> 0.647), 3KTQ (contact_auc 0.582 -> 0.296, stems_nearer_frac 1 -> 0), 9K4F, 4NLG, 2VBO, 5W1H (contact_auc), and the stronger pull added NA x protein crossings (6JCY, 7XSP, 9CJI, 4OO8). Its case stays a limitation. The change as measured: reverted/F3.patch.
- Steps: F6 -> F3.

Target:

| entry | measure | before | after |
|---|---|---|---|
| 8VBG | contact_auc | 0.505 | 0.521 |
| 8VBG | stems_nearer_frac | 0.5 | 0.5 |
| 6QIK | na_crossings | 1 | 3 |

Cases that moved:

- 8VBG: contact_auc 0.505 -> 0.521; na_protein_crossings 0 -> 1; contact_precision_at_k 0 -> 0.125; ss_contacts_nearer 1 -> 0.667; named_fallbacks 0 -> 1 (worse: na_protein_crossings, named_fallbacks)
- 6ZY7: procrustes 0.227 -> 0.179; contact_auc 0.823 -> 0.835; na_protein_crossings 1 -> 9; contact_precision_at_k 0.071 -> 0; ss_contacts_nearer 0 -> 1 (worse: na_protein_crossings)
- 6QIK: spearman 0.752 -> 0.722; procrustes 0.217 -> 0.301; contact_auc 0.948 -> 0.914; stems_nearer_frac 1 -> 0.984; na_crossings 1 -> 3; na_run_crossings 404 -> 425; na_protein_crossings 360 -> 478; contact_precision_at_k 0.102 -> 0.1; ss_contacts_nearer 0.998 -> 0.995; named_fallbacks 54 -> 56 (worse: spearman, procrustes, contact_auc, stems_nearer_frac, na_crossings, na_run_crossings, na_protein_crossings, named_fallbacks)
- 8XJV: spearman 0.715 -> 0.732; procrustes 0.409 -> 0.388; contact_auc 0.955 -> 0.956; stems_nearer_frac 0.969 -> 1; fallbacks 10 -> 13; na_crossings 10 -> 11; bead_overlaps 127 -> 175; na_run_crossings 191 -> 212; na_protein_crossings 725 -> 884; contact_precision_at_k 0.115 -> 0.103; ss_contacts_nearer 0.965 -> 0.979; named_fallbacks 111 -> 112 (worse: fallbacks, na_crossings, bead_overlaps, na_run_crossings, na_protein_crossings, named_fallbacks)

Controls that moved:

- 6JCY: contact_auc 0.778 -> 0.843; na_protein_crossings 8 -> 10; contact_precision_at_k 0.091 -> 0.22; ss_contacts_nearer 0.917 -> 1 (worse: na_protein_crossings)
- 4NLG: contact_auc 0.808 -> 0.695; na_protein_crossings 4 -> 0; contact_precision_at_k 0.318 -> 0.227 (worse: contact_auc)
- 2VBO: contact_auc 0.778 -> 0.738; contact_precision_at_k 0.723 -> 0.585 (worse: contact_auc)
- 9K4F: spearman 0.869 -> 0.794; procrustes 0.358 -> 0.384; contact_auc 0.845 -> 0.824; stems_nearer_frac 0.875 -> 1; na_run_crossings 4 -> 5; contact_precision_at_k 0.271 -> 0.229; ss_contacts_nearer 0.895 -> 0.947 (worse: spearman, procrustes, contact_auc, na_run_crossings)
- 7XSP: contact_auc 0.801 -> 0.849; na_protein_crossings 0 -> 2; ss_contacts_nearer 0.611 -> 0.889 (worse: na_protein_crossings)
- 9CJI: procrustes 0.271 -> 0.236; contact_auc 0.826 -> 0.862; na_protein_crossings 0 -> 7; contact_precision_at_k 0.25 -> 0.188 (worse: na_protein_crossings)
- 1KX5: spearman 0.881 -> 0.891; procrustes 0.332 -> 0.314; contact_auc 0.92 -> 0.924; na_run_crossings 6 -> 5; contact_precision_at_k 0.318 -> 0.326
- 1IHF: procrustes 0.274 -> 0.389; contact_auc 0.786 -> 0.647; na_protein_crossings 4 -> 0; ss_contacts_nearer 0 -> 1 (worse: procrustes, contact_auc)
- 1RVA: contact_auc 0.609 -> 0.612; contact_precision_at_k 0.25 -> 0.125
- 1HYS: contact_auc 0.836 -> 0.902; contact_precision_at_k 0.222 -> 0.556
- 3KTQ: contact_auc 0.582 -> 0.296; stems_nearer_frac 1 -> 0; na_protein_crossings 4 -> 2; contact_precision_at_k 0.125 -> 0 (worse: contact_auc, stems_nearer_frac)
- 4OO8: spearman 0.85 -> 0.854; procrustes 0.332 -> 0.203; contact_auc 0.845 -> 0.888; na_run_crossings 4 -> 3; na_protein_crossings 11 -> 12; contact_precision_at_k 0.167 -> 0.217; ss_contacts_nearer 1 -> 0.895 (worse: na_protein_crossings)
- 3CRX: procrustes 0.817 -> 0.765; contact_auc 0.687 -> 0.752; na_run_crossings 5 -> 4; na_protein_crossings 6 -> 5; contact_precision_at_k 0.207 -> 0.239
- 5W1H: contact_auc 0.765 -> 0.681; ss_contacts_nearer 0.828 -> 0.897 (worse: contact_auc)

Controls worse: 6JCY (na_protein_crossings); 4NLG (contact_auc); 2VBO (contact_auc); 9K4F (spearman, procrustes, contact_auc, na_run_crossings); 7XSP (na_protein_crossings); 9CJI (na_protein_crossings); 1IHF (procrustes, contact_auc); 3KTQ (contact_auc, stems_nearer_frac); 4OO8 (na_protein_crossings); 5W1H (contact_auc)

## Limitations

- F1 (M1, pairs missed at low resolution) deferred past the freeze: 8XJV pair recall, 3EGZ A2.U64 (pair_f1 0.941) and 6ZY7 stem_recall 0.667 (ruling Q4) stay failures.
- 3EGZ procrustes 0.572 per copy (limit 0.5): the stems of one copy span 49 degrees round the axis, unrolled flat.
- 6QIK (F5 not committed): fallback, crossing and the A3 strand y:8-9 remain.
- 8VBG (F3 not committed): contact_auc 0.505, stems_nearer_frac 0.5 remain; 9YLM, 9K4F, 30QM, 3KTQ, 1YTB as before.
- 8XJV: everything M1 explains (fallbacks 10, na_crossings 10, bead_overlaps 127, time).
