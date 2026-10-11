# Nucleic-acid campaign review (Task 13, phase A)

Run 2026-10-11, `--require-rnapolis`, RNApolis venv `/home/dzyla/.claude/jobs/23e71d62/tmp/na-venv`: **26/39 pass**.
Every contact sheet was opened; figures were also re-rendered at 120 dpi (scratch, not committed) to read letters and strands.
Cause tags (C1–C14, M1–M5) are explained at the end. One line per finding: `<ID>: <what is wrong> → <suspected cause>`.

## Findings

- 1EHZ: stems drawn rotated against 3D (procrustes 0.716; P3 left of P2 instead of above it) → C2
- 1EHZ: tie labels "a" (G18·Ψ55, G19·C56) about half a base letter's size, unreadable → C13
- 1QTQ: tRNA drawn as a block above the protein; P1 and the acceptor end touched by the catalytic domain are drawn on the far side of the tRNA (contact precision@k 0.25) → C2 (weak home pull), accepted otherwise
- 1QTQ: tie labels unreadably small → C13
- 4LVV: all six stems stacked in one column; in 3D P3/P4 sit beside P2a/P2b/P1 → C2
- 4LVV: P3 drawn between P2a and P2b, the two pieces of one duplex 8.6 apart instead of end to end → C2, C3
- 4LVV: strand A:78 passes through its own stem body (P2a/P4); 7 run–run crossings (e.g. below P3) → C1, C12
- 4LVV: contact-sheet panel titles overlap (narrow panels) → C14
- 2GIS: strands from P2 to P3 (A:18-20, A:32-37) drawn as straight diagonals through P3's own body, an X across P3 and P2 → C1
- 2GIS: P6 drawn under P1 and P2 under P3 (3D: P6 bottom-left, P2 left); Spearman 0.601 just passes → C2
- 2GIS: 6 run–run crossings (e.g. the P7 run across the P1–P5 run), two strands drawn 1 mm apart over 25 nt below P5 → C12
- 3D0U: P6a drawn three stems below P6b (3D: adjacent, beside P2); the joint is two long vertical strands → C2, C3
- 3D0U: 10 run–run crossings → C12
- 2GDI: run X:39-47 bows through P2's body with its letters "A A A" inside the stem; X:19-20 through P2/P3; P8→P9 diagonal through P9's letters (9 own-body passes) → C1
- 2GDI: both crystal copies (X, Y) drawn end to end as if one molecule → C13 (assembly auto keeps the ASU)
- 1Y26: junction strand under P2/P4 runs along P1's top face with "U U U C" squeezed between it and P1; tie labels tiny → C12, C13
- 1GID: X of strands through P3 (A:122-126, A:196-199) and P13 (B copy); letters A:125 and A:197 collide there (the letter_overlaps 1) → C1
- 1GID: na_crossings 2: run B:139-140 (P14→P15) crosses tail B:103-106; run B:215-219 (P11→P17) crosses tail B:256-260 → C8
- 1GID: 19 run–run crossings, unmeasured → C12, M3
- 1GID: both crystal copies drawn → C13
- 3IGI: 6 own-body passes (P8a, P7/P6, P15, P14, P1) and 19 run–run crossings; narrow spike loop at P16 → C1, C12
- 4P9R: chain A drawn as three molecules (5′/3′ at A:22/A:216 and A:347/A:542, 1–2 nt unmodelled) → C6
- 4P9R: 11 run–run crossings; tie labels tiny → C12, C13
- 2OEU: P2 and P3 swapped left/right against 3D (procrustes 0.747 passes) → C2
- 2OEU: run A:4-8 bows through P1/P2; letters "A G U C" crowded in the P1–P2 junction → C1
- 2NZ4: beads_missing 2: R:U110 is left by the crop as a one-nucleotide strand piece that is never drawn, so C:LYS23/ASP24 have no bead → C7
- 2NZ4: gap_label_collisions 3: two "… 6 nt …" labels on run Q:32-34, on element C:40-44 and on a protein loop → C10
- 2NZ4: block first-pair labels ("12·22", "29·45") float between stems, ambiguous; four crystal copies drawn → C13
- 4OO8: na_crossings 1: run B:81 (P5→P6) crosses the 3′ tail B:97 → C8
- 4OO8: the two Cas9 copies are segmented differently (heteroduplex P1 straight, its copy split P7a/P7b) and named P1..P12 across copies → C13
- 5B43: heteroduplex split into P2a/P2b for a 3D turn whose in-page bend is 4.1°, pieces drawn 10 apart with protein between → C3
- 5W1H: fallbacks 1: the 17-nt crRNA spacer (tail B:-3..13) is drawn as a straight line through α28/α29/α52 → C8
- 5F9R: stem_recall 0.714: A71-73/A76-79 (C73·G76 closes a 2-nt loop, excluded by MIN_LOOP) and A29-36/51-58 vs foldmap A29-38/49-58 (A38·A49 vs RNApolis A39·A49) → M4, C11
- 5F9R: na_crossings 1: run C:10-11 (P1b→DNA 1) crosses the non-target tail D:12-20 → C8
- 5F9R: P1a and P1b (one heteroduplex) drawn 10 apart with α34 between, the joint routed through protein → C3
- 30QM: stems_nearer_frac 0.833: P3 (touched by A:288-306, A:346-362, B:37-57) drawn 18–24 away from all three; the two protomers are pushed to the page edges → C4
- 30QM: contact sheet 3D panel shows seam lines across the page (unrolled frame) → C14
- 1URN: missing U13–C14 of chain R splits the hairpin into two strands; it is drawn as the classic straight two-strand duplex with two 5′ and two 3′ ends → C6
- 4PJO: fallbacks 1: the 16-nt tail 1:123-138 is straight through protein; na_crossings 1: run 3:11 (P5→P6) crosses tail 3:123-137 → C8
- 4PJO: four snRNP copies drawn, figure unreadable at page size → C13
- 1HQ1: the symmetric internal loop between P2 and P3 is drawn as an X of its two sides (B:147-150, B:161-164 through P2/P3), beads and letters piled on the crossing → C1
- 1RLG: k-turn bulge C:5-7/D:5-7 bows through P1/P2 and P3/P4; P2 drawn above P1 where 3D has it to the left (both copies, procrustes 0.711) → C1, C2
- 1ZDI: hairpins drawn beside the coat-protein sheets they lie on (contact precision@k 0.0) → C2 (weak), accepted
- 1HYS: stem_recall 0.0: the 22-bp PPT hybrid is split into E855-867 and E870-876 because E868·F886 (N6–O4 3.98 Å) and E875·F879 (4.51 Å) fail HBOND_MAX 3.8 → C11
- 1MSW: pair F1 0.917 (N23·T123, N29·T129 missed); drawing fine → C11 (not judged)
- 3KTQ: contact_auc 0.582: primer/template split into DNA 1a/1b (in-page bend 3.6°) and drawn 7 apart with α19 between; one joint strand detours under α19 → C3
- 4YBB: fallbacks 3, na_crossings 6: the crop leaves AA:1010-1015 as an unpaired piece, drawn as a trace at its projected place far left; its two gaps go straight through P59, P60, P61 → C7
- 4YBB: gap label "… 21 nt …" on element AC:26-46 and runs AA:991-996, AA:1211-1214; bead_overlaps 4 (AC/AL, AC, AD, AS) → C10, C13
- 4YBB: 23S rRNA of the neighbouring ribosome (DA) drawn, touched through a lattice contact by the selected 30S proteins → C13 (campaign option: `chains` cannot name RNA chains)
- 4YBB: 31 own-body passes, 116 run–run crossings → C1, C12
- 1KX5: superhelix drawn as two concentric rings of pieces with long joint loops; faithful, 11 run–run crossings → accepted, C12
- 1IHF: U-turn around the protein, good; joint DNA 1b→1c routed under the β-arms → accepted
- 1YTB: each TBP β-sheet drawn between the two halves of its DNA (DNA 1a/1b, 2a/2b 13 apart), joints run over the sheet panels (12 protein crossings) → C3
- 1LMB, 1ZAA, 1BNA, 1ZBI: classic picture with the first arm at the bottom: the top strand reads 3′→5′ left to right → C13 (check intended)
- 1RVA: a ~50° kinked duplex drawn as two antiparallel pieces side by side (bend 180°) joined by long loops below → C4; the metrics agree with it by construction → M2
- 3CRX: the assembly lays a copy of every DNA strand on top of a deposited one (C1≡E2, D1≡F2 …, all nucleotides within 1.5 Å); pairing then alternates between twins and the four arms become 22 stems of 2–5 bp scattered through the protein; 18 run–run crossings; bead_overlaps 1 → C5; metrics pass → M1
- 4LD0: procrustes 0.996, contact_auc 0.506, stems_nearer_frac 0.5: each arm drawn under the protomer that does not touch it (DNA 1 touched by B under A; DNA 3 touched by A under B), DNA clustered below both → C4
- 1ZBI: "3′" end label touched by a protein loop → C13 (minor)
- 1KF1: trace self-crosses (faithful); a few letters sit on the trace → accepted
- 4PLX: chains B and C break at 55/59 and 54/60: four extra 5′/3′ labels piled between P7 and P4 → C6
- 4PLX: run B:47-49 through P4, C:47-50 through P7 → C1
- 1QBJ: fallbacks 1: unpaired hexamer F (its partner is a crystal mate) drawn as a projected trace over chain C α1, beads on the helix → C9
- all: angle_median, angle_max and joint_error_max are 0.0 on all 39 entries → M2

## Causes

- **C1** A strand drawn through its own stem's core. route_na._strand_between bows an internal loop from rail point to rail point (and a route may end at a far face); _na_clear checks outside_own(pts), which trims the whole first and last stretch inside an own stem, so a bow that enters a stem at its near side and runs to a rail point at the far end is never seen. nacheck._trimmed trims the same way.
- **C2** Stem positions follow the strand springs, not 3D. layout._stem_springs gives every run weight _W_SEQ = 1 and slack 1.5; home pull is _K_HOME 0.02 × mass. Run springs are already stretched at the projected positions (4LVV P3→P4 7.5 vs slack 1.5), so the optimiser contracts the RNA into a sequence-order chain.
- **C3** Consecutive pieces of one bent duplex cannot sit end to end: _attract's overlap term uses axis-aligned footprints (+ MARGIN, label and hairpin room) between joint-linked pieces. Measured centre gaps: home −2.9…−0.5, after _attract +6.5…+12.1 (3KTQ, 5B43, 1YTB, 4LVV). _W_JOINT 4 → 40 does not close them. segment_stem also splits on a 3D turn whose in-page bend is a few degrees (5B43, 3KTQ).
- **C4** The unrolled cylinder frame (symmetric dimers) cuts a nucleic acid lying across the axis: each half is unrolled with its protomer (1RVA, 4LD0, 30QM: the three CylinderFrame entries).
- **C5** read_model(auto) builds an assembly whose operator copies DNA onto itself (DNA on a crystallographic axis). Nothing removes coincident copies.
- **C6** An internal chain break (unmodelled nucleotides) starts a new strand with its own 5′/3′ labels; nothing marks it as a gap.
- **C7** The auto crop (_na_extent) keeps strand pieces holding no stem: a 1-nt piece (2NZ4) is not drawn at all; a longer one (4YBB) becomes a trace at its projected place, unrelated to the layout.
- **C8** Free ends are drawn straight with at most one bend, whatever their length, and placed before the runs, which only pay a cost to cross them.
- **C9** Unpaired strands (traces) stay at their projected place; the protein is not kept off them.
- **C10** Gap labels sit at the gap's midpoint with no collision search.
- **C11** pair_nucleotides needs every hydrogen bond ≤ HBOND_MAX 3.8 Å; stretched A·T/A·U pairs with one bond at 3.9–4.0 Å are lost and split stems (1HYS).
- **C12** Run–run crossings are neither avoided (grid reuse is only a cost) nor measured.
- **C13** Cosmetic or campaign-option issues, listed per entry.
- **C14** Contact-sheet tool (scripts/na_validate.contact_sheet): titles overlap on narrow panels; unrolled frames draw seams in the 3D panel.

## Metrics that are wrong rather than the drawing

- **M1** pair_f1 and stem_recall are measured on the ASU (`read_model(path, "asu")`), not on the model that is drawn. 3CRX passes with stems_foldmap 4 while the drawn assembly has 22 shattered stems.
- **M2** shape() measures through lay.frame, the same projection the layout uses, so angle and joint errors are 0.0 by construction on all 39 entries. On an unrolled frame it agrees with the distortion: 1RVA bend_3d = bend_drawn = −179.7° for a ~50° kink.
- **M3** na_crossings counts only (run/joint/gap × hairpin/tail) and passes through other stems. It misses run–run crossings (1GID 19, 3IGI 19, 3CRX 18, 1KX5 11, 4P9R 11, 3D0U 10, 4LVV 7, 2GIS 6, 4OO8 6, 4YBB 116). Because of _trimmed it also misses strands through their own stems (C1).
- **M4** stem_recall: reference pairs that break foldmap's own MIN_LOOP rule (5F9R C73·G76, 2-nt loop) are grouped into reference stems, and the ±1 end tolerance fails a stem that differs only in a 1×1-mismatch extension (5F9R A29-36 vs A29-38, where foldmap's A38·A49 has the better geometry: N1–N3 3.43, N6–O4 3.21 Å against 3.18/4.59 Å).
- **M5** (for the record) contact_auc and stems_nearer_frac on 30QM, 4LD0 and 3KTQ are right: they see C4 and C3.

## Task 13b step 1 (C5, C4, M1, M2): rerun of all 39 entries, 26/39 pass (run-step1.log)

- 3CRX: C5 fixed: the assembly keeps 4 Cre and 4 strands (copies C2–F2 dropped), drawn as four arms of 16–17 bp (was 22 stems of 2–5 bp); procrustes 0.817 now fails (limit 0.8): the projected frame sees the junction edge-on (arms 3, 4, 1, 2 in a row) while the figure draws an X → frame orientation, needs a ruling; fixture held at 0.817 in tests/test_na_validation.py HELD
- 30QM: C4 fixed: plain frame, stems_nearer_frac 1.0 (was 0.833), Spearman 0.893, procrustes 0.214, contact_auc 0.935
- 4LD0: C4 fixed: each arm beside the protomer that touches it, procrustes 0.496, contact_auc 0.839, stems_nearer_frac 1.0
- 1RVA: C4 removes the 180° hairpin, but the projected frame sees the duplex nearly end-on (projected axis lengths 0.50 and 0.34): its 39° kink (angle between piece axes) is drawn as its 102° projection, and the DNA lands beside the dimer instead of in its cleft, contact_auc 0.58 (was 0.686) → frame orientation, needs a ruling; fixture held at 0.58 in HELD, the ruled "< 90°" pin is a strict xfail
- 1URN: pair_f1 0.857 (drawn pairs, M1; deposited 0.96): focus draws one hairpin with 4 reference pairs and foldmap misses A2·U19 (N1–N3 3.57, N6–O4 3.86 Å) → C11 (surfaced by M1, not a new defect)
- M1: drawn pairs equal the deposited ones on every other entry; cropped 4YBB 0.982 (deposited 0.985), 2NZ4 and 3KTQ now 1.0 (their unmatched pairs are outside the drawn model)
- M2: no entry other than 1RVA (before C4) was unrolled with nucleic acid; the RadialFrame entries (1LMB, 1ZAA, 1BNA, 1ZBI …) measure the same through the plain frame at the duplex centre

## Task 13b step 2 (C2, C3): rerun of the 13 affected entries (run-step2.log, report-step2.md; report.md merges them into the 39-entry table, the other 26 rows are from step 1)

- C2 (runs pull only beyond the distance 3D shows): procrustes 1EHZ 0.716 → 0.188 (the tRNA in its L shape, not a cloverleaf), 2OEU 0.747 → 0.193 (P2/P3 in their 3D order), 1RLG 0.711 → 0.177, 4LVV 0.477 → 0.200 (Spearman 0.811 → 0.961), 3D0U 0.403 → 0.352
- C2 side effects, routing (step 3): 1RLG na_crossings 0 → 4 (runs across hairpin arcs of the stems now side by side as in 3D); fixtures 4OO8 bead_overlaps 1 and 4LD0 na_crossings 1 held in tests/test_na_validation.py CAMPAIGN ("routing, step 3"). Stem names are now hard obstacles for a run's first drawing (4YBB's P27–P26 run went across P27's name)
- C3 (consecutive pieces of one duplex keep no margin between them; chains of at most four are separated together; joint weight unchanged): elements between consecutive pieces 1IHF 6 → 2, 3KTQ 2 → 1, 5B43 3 → 4, 1YTB 4 → 0, 4LVV and 3D0U 1 → 0 (P3/P5a no longer between P2a/P2b and P6a/P6b); 1RVA's two pieces meet across the kink, contact_auc 0.58 → 0.609 (passes, hold removed)
- C3 accepted as is: protein that wraps DNA in 3D projects between its segments (1IHF arms, 3KTQ α19, 5B43). Holding the joints firmer (joint weight raised with the overlap penalty) removes every such element but collapses the nucleosome ring (1KX5 contact_auc 0.922 → 0.561, histones pushed outside) — the ring with the histones inside is the truer picture. 1IHF pinned at ≤ 2 as a regression bound
- C3 side effects: 5B43 letter_overlaps 0 → 1; 1YTB contact_auc 0.872 → 0.767 (still passes); 1YTB's and 5B43's halves still sit apart, joined by long routed joints. (Fix round 1: the first C3 version moved all 22 nucleosome pieces as one group and pushed 12 histone elements out of the ring, NA×protein crossings 0 → 54; that was the grouping, not histone tails. Groups are now chains of at most four pieces: 1KX5 has every element inside the DNA, na_crossings 0, NA×protein crossings 2)
- Metric: nacheck `overlaps` no longer counts two consecutive pieces of one duplex (layout.joined_pieces): they meet at their joint
- Follow-up: 3KTQ contact_auc 0.582 (limit 0.6) is unchanged by C3: with α19 moved out from between DNA 1a and 1b the protein that grips the DNA sits farther from it (0.371–0.555 in the firmer-joint variants), so the AUC fails either way; needs a placement idea for a protein that encircles DNA
- 3CRX procrustes 0.817 unchanged (no runs or joints: strand exchanges only), still held; frame orientation as ruled after step 1
- Fix round 1 rerun (run-step2.log): 8/13 pass (3D0U and 1KX5 pass again); a straight joint across touching pieces is allowed only near or touching, forward, clear by _INFLATE and crossing no strand but its partner at the same joint

## Task 13b step 3 (C1, C8, C12, M3): routing; rerun of the 14 affected entries (run-step3.log, report-step3.md: 10/14 pass) and of all 39 (run-step3-all.log; report.md is now the full 39-entry table): **30/39 pass** (was 24/39 at the start of step 3, every row re-measured)

- C1 (a run or joint runs along its rail to its stem's end face and is drawn from face to face; internal-loop sides are bowed face to face; between the faces its own stems are obstacles; route_na.among trims exactly the two rail stretches, as nacheck._trimmed now does): the X shapes through own stems are gone (2GIS P3/P2, 1GID P3/P13, 2GDI P2, 1HQ1 P2/P3, 1RLG k-turn); 1GID letter_overlaps 1 -> 0 (A:125/A:197 sat in the X), 1RLG na_crossings 4 -> 0. Pinned on 2GIS (new RNA-only fixture, 42 KB) and six fixtures: no run or joint enters an own stem's core
- C8 (hairpin loops and free ends are hard obstacles in a run's first routing pass; a free end of more than 6 nt that no straight or turned drawing keeps clear is routed out to free space on the grid; a free end keeps out of its partner's doorway): 5W1H fallback 1 -> 0 (the 17-nt spacer now runs along free space), 4PJO fallback 1 -> 0, 4LD0 na_crossings 1 -> 0 (fixture hold removed), 5F9R and 3IGI na_crossings 1 -> 0, 2NZ4 1 -> 0
- 4PJO na_crossings 1 stays: P6's 15-nt 3' end 3:123-137 is boxed in by protein loops; the drawing that keeps out of the doorway of the run 3:11 crosses 5 protein loops (NA x protein 19 -> 24), so the end is drawn straight across the doorway and the run crosses it (accepted trade: NA x protein crossings must not grow)
- C12 (a run with nucleotides is drawn off every strand drawn before it and every protein loop when that is at most 3x as long as its plain drawing; else one that also pays for the cells beside strands and protein loops, at most 2x; strand exchanges and joints keep their shortest way, paying only beside protein loops): over all 39 entries run-run crossings 255 -> 237 and NA x protein crossings 192 -> 158; no entry's NA x protein crossings grew (30QM 17 -> 6, 4OO8 23 -> 16, 5F9R 23 -> 18, 4YBB 51 -> 48, 1KX5 2 -> 1). Run-run crossings grew on 3IGI 16 -> 18, 4OO8 4 -> 5, 5F9R 3 -> 4 (reported, not judged)
- M3: na_run_crossings (runs, joints and gaps crossing one another, the two strands of one joint excepted) is reported, not judged; strands through their own stems now count in na_crossings (4YBB: the straight joint AA:253-254 clips its own piece, the fixture's one na_crossing)
- Carried from step 2: 4OO8 bead overlap (the two sides of the P2-P3 internal loop cross where the stems are seen to; a run's nucleotides now slide along it off another strand's, hold removed); 5B43 letter_overlaps 1 -> 0 (the 5' end B:-19--16 between its stem and a sheet: a crowded letter looks further out and along its strand); 1RLG bead overlap that C1 exposed (two runs one grid cell apart) -> 0 (beside-cell cost); 4YBB bead_overlaps 6 -> 1, fixture hold removed
- Still failing, other steps: 2NZ4 beads_missing/gap label (C7, C10), 5F9R stem_recall (M4), 1URN pair_f1 and 1HYS stem_recall (C11, C6), 3KTQ contact_auc (C3 follow-up), 4YBB fallbacks 5 / na_crossings 4 / gap label (C7, C10), 3CRX procrustes (held, frame), 1QBJ fallback (C9)
- Contact sheets: the text area now grows with the metric lines (the extra metric made narrow sheets' metrics run into the "failed:" line)
- Fix round 1 (run-step3-fix1.log for the 10 affected entries; full rerun run-step3-all.log): **29/39 pass**
  - Detours capped absolutely: a staged drawing is at most 6 nucleotide steps longer per nucleotide (and end) than the plain one, and never a fallback over a clean drawing. Before, 4YBB's P72->P71 was drawn at 5.7x its direct distance and 1GID's P15->P16 at 3.6x; now no run on any entry exceeds 3.2x. NA x protein crossings are 192 -> 156 over the step, and no entry is above its step-2 value. Run-run crossings: 244.
  - C4 revised (controller ruling): an assembly stays unrolled side by side unless that tears its nucleic acid apart. That means consecutive pieces of one duplex facing different protomers (1RVA), or linked stems turned against one another by more than 90 degrees between the two frames (4LD0 139; 30QM at most 57).
  - 30QM is drawn side by side again (NA x protein 6 -> 1). It fails stems_nearer_frac 0.833, accepted with cause: one stem held by both protomers (P3) sits between them. It also fails spearman 0.521, angle_median 35.7 and angle_max 109.6, accepted with cause: these measure the unrolling itself against the plain projection of the same orientation (M2). All four are held for the fixture in tests/test_na_validation.py (NEARER and HELD).

## Task 13b step 4 (C6, C7, C9, C10, C11, M4, COMPARE #10 and #6): rerun of the 11 affected entries (run-step4.log, report-step4.md: 8/11 pass, before the two follow-ups) and of all 39 (run-step4-all.log; report.md is the full table): **33/39 pass** (29/39 after step 3)

- C11 (second pairing pass: one bond within 3.8 Å, the other within 4.0 Å, only stacked in register on a strict pair): over the 39 entries and 1AOI it adds 28 pairs, 27 of them RNApolis pairs (the one extra, BA:1006·BA:1023, is in 4YBB's second ribosome, outside the drawing); no entry's precision drops (deposited 4YBB 0.992 -> 0.992, 3IGI 0.983 -> 0.983). 1HYS E868·F886 joins the PPT stem (stem_recall 0.0 -> 1.0), 1URN A2·U19 (pair_f1 0.857 -> 1.0), 1MSW pair_f1 0.917 -> 0.96. Already below 0.98 before this step and unchanged by it: 4LVV 0.968 (A26·A53), 1HQ1 0.941 (B149·B163), the 4YBB fixture window 0.977
- COMPARE #10: a 2x2 loop in register whose unpaired steps still sit across the helix (C1'-C1' within 9-12 Å) stays inside its stem: 1AOI is one 146-bp duplex "DNA 1" (was DNA 1/2/3: I20-21 unpaired, I24·J269 and I25·J268 recovered by C11), cut only where it bends (ad hoc run: PASS, pair precision 0.993, stem_recall 1.0); also merges one stem pair on 2NZ4 and eight on 4YBB
- COMPARE #6: a pseudoknotted stem (crossing a longer one within one chain; of equal ones the later) is named PK (PK1, PK2 ... when several) and nests nothing; RNA helices are numbered P1, P2 ... in 5' order and lettered across internal loops and bulges (P2a, P2b); pieces of a lettered stem are numbered (P2a1). 2GIS reads P1, P2a/P2b, PK, P3a/P3b, P4, as in Roy 2019. 4PLX gets PK1-PK3 (the A-rich 3' tract pairs across the hairpins), 2NZ4 PK1-PK3
- C6: an internal chain break is a gap of the molecule (counted by label_seq: 4P9R A:22 -> 216 is 1 nt), dashed, with no 5'/3' at it; the classic straight picture never takes two strands of one chain. 1URN: one hairpin with a bare break over its loop (was a classic two-strand duplex with two 5' and two 3' ends); 4P9R one molecule; 4PLX no stray end labels. The 4PLX break B:55-59 found no way on the grid and went straight through PK2: a gap the grid cannot route is now curved along its ends, or bowed, first (follow-up)
- C7: the crop to stems drops pieces holding no stem and merges their gaps: 2NZ4 beads_missing 2 -> 0 (R:U110's touchers cropped, not beadless), 4YBB fallbacks 5 -> 0 and na_crossings 4 -> 1
- C9: traces are fixed obstacles in placement: 1QBJ fallback 1 -> 0 (chain C's alpha1 off hexamer F); contact_auc 0.893 -> 0.786 (passes): the protein moved off the trace it touches
- C10: gap labels placed along their gap, either side, clear of bodies, strands, letters, beads, names: 4YBB gap_label_collisions 1 -> 0; 2NZ4 1 -> 0 once the label box was sized from the page's text (follow-up: a letter glyph per character made it 2.6 times the drawn text)
- M4 (metric, not drawing; fix round 1: one to one): stem_recall leaves out reference pairs closing a loop under MIN_LOOP and matches reference and foldmap stems one to one, most shared pairs first, a reference stem recalled when its match holds at least two thirds of its pairs. The old rule (both arms' ends within one nucleotide) was wrong for 5F9R: RNApolis' A29-36/51-58 is wholly inside foldmap's A29-38/49-58, which carries on across a 1x1 mismatch with A38·A49 (N1-N3 3.43, N6-O4 3.21 Å), a better pair than RNApolis' A39·A49 (N6-O4 4.59 Å); and C73·G76 closes a 2-nt loop no foldmap pair may close. 5F9R 0.714 -> 1.0, 4LVV 0.8 -> 1.0, 1GID 0.938 (15/16), 3IGI 0.818 -> 0.909 (20/22). (The first version, any foldmap stem holding half the pairs, let one foldmap stem recall several reference stems: 1GID and 3IGI read 1.0.)
- 4OO8: na_crossings 0 -> 1, new: the straight joint P6a -> P6b (pieces that touch, step 2's rule) cuts the corner of its own piece P6a on the way to its face; _straight_joint checks every body but its own two. Requiring it to keep off its own two (tried) turns this joint and 4YBB's into fallbacks, so it is left -> routing follow-up. 4YBB's one na_crossing is the same kind (joint P11b1 -> P11b2)
- NA x protein crossings 156 -> 174 over the 39 (reported, not judged): 4YBB 48 -> 76 (placement of 82 cropped stems moved with the pairs C11 adds, 48 -> 60, and the longer helix names, 59 -> 83; C7 brings it to 76), 5F9R 18 -> 8. Run-run crossings 244 -> 242
- Still failing, with cause: 30QM (held, unrolled), 3CRX procrustes (held, frame), 3KTQ contact_auc (C3 follow-up), 4PJO na_crossings (accepted trade, step 3), 4OO8 and 4YBB na_crossings 1 (straight joint through its own piece, above)

### Task 13b step 4, fix round 1 (run-step4-fix1-all.log; report.md re-measured): **34/39 pass**
- C9 revised: each trace nucleotide (a bead either side) is the obstacle, and an element touching the trace is exempt (as is the sheet holding it). 1EYG (new ad hoc entry and 82 KB fixture): the SSB tetramer stays under its ssDNA (the box version split it and pushed it below); its traces are no fallbacks. 1QBJ contact_auc 0.786 -> 0.857 (0.893 before C9). 1EYG still fails bead_overlaps 1: two traces' nucleotides coincide where the projected strands meet, which a trace cannot move to avoid
- M4 revised (one to one, two thirds): see the M4 line above; 1GID 0.938, 3IGI 0.909
- Helices held together across loops of at most 3 nt a side (strand-exchange spring): NA x protein crossings 174 -> 152 (4YBB 76 -> 54), run-run 242 -> 232; 2GIS P3a/P3b drawn end to end (loop sides 12.8/13.6 -> 4.4/4.7). Run-run 3IGI 16 -> 18 (reported)
- Straight joints: a shallow graze of their own piece's corner is kept and not counted (4OO8 P6b -> P6a: na_crossings 1 -> 0, PASS); a deeper one turns an elbow where one clears
- 4YBB na_crossings 1 remains, now a different case: a gap leaves a free end drawn inside another stem's hairpin loop (strand AA cropped piece 4 into P-hairpin 631-636). 4PJO, 3KTQ, 3CRX, 30QM as before

## Task 13b step 5 (C13, C14, COMPARE #4; carried from step 4: free ends and hairpin loops, coincident trace nucleotides, per-segment trace exemption)

- COMPARE #4: when the asymmetric unit holds several copies of a nucleic-acid complex and assembly 1 is one of them (identity operator only, a strict subset, every polymer entity deposited at least twice as often), `auto` draws assembly 1 under its deposited chain names; symmetry is then looked for among its copies only. Changed entries: 1GID, 2GDI, 2NZ4, 4OO8 (one Cas9 complex, A, B, C), 4PJO, 1RLG, 1URN (one U1A-RNA complex, P + A: no "helical" assembly; was a 120.3 degree screw of three crystal copies), 1YTB, 4PLX (and ad hoc 1EYG). Protein-only structures keep the asymmetric unit (4CH9, two copies of a dimer, would change; the protein bench is identical). na_validate scores these against the deposited model as deposited chains (as_deposited: every drawn chain a deposited chain under its own name)
- C13: `chains` names nucleic-acid chains as well as protein chains; each kind is restricted only when some of it is named. 4YBB's campaign option now names AA with its proteins: the neighbouring ribosome's 23S rRNA (DA, a lattice contact) is no longer drawn. Tie labels stay at 0.75x (ruling); no per-entry single-copy options (COMPARE #4 does it for every entry)
- C14: contact-sheet titles are wrapped to their panel (4LVV, 3D0U, 1GID, 2NZ4 now wrap); an unrolled frame's 3D panel breaks each polyline where a step is longer than half the protomer period (the seam) and projects the protein by its protomer, as the figure does (30QM, 1ZDI: no lines across the panel)
- Free ends: a free end's drawing never comes into another stem's hairpin circle when another drawing keeps out of it (4YBB step 4: the gap leaving a free end inside a hairpin loop crossed it)
- Traces: a trace nucleotide projected within a bead's width of one drawn before it (another strand's, or its own neighbour seen end-on) is nudged sideways off its trace, as little as clears it, on the grid. The 1EYG fixture's bead_overlaps 1 -> 0 (Q:16-17 runs toward the viewer, 0.31 page units apart)
- route_na._trace_clear: an element touching a trace is exempt only for the segments within two nucleotides of those it touches, not for the whole trace (no fixture or campaign trace became a fallback)

## Final run (2026-10-11, after step 5b; run.log, report.md, report.json): **33/39 pass**

Every contact sheet was viewed again. Two runs of all 39 gave identical reports. Over the 39: NA x protein crossings 152 -> 108 and run-run crossings 232 -> 206; most of the drop comes from drawing one complex instead of every crystal copy. New review lines:
- 1URN: one complex (RNA P with U1A chain A). The chain break's dashed gap runs round the protein's loops (NA x protein 0 -> 2)
- 1GID: one P4-P6 domain; spearman 0.601 and procrustes 0.405 (both pass), so the single copy sits close to the Spearman limit
- 1LMB, 1ZAA, 1BNA (radial or straight frames): the 3D panel of a duplex wider than 2.5:1 is shrunk in height, so its title sits mid-sheet. This affects the contact sheet only
- Task 14 limits: with the new limits (procrustes 0.50), 30QM also exceeds procrustes (0.776); it is held as part of its accepted unrolled figure (like spearman and the angles)
- 1EYG (ad hoc, not in the 39): drawn as assembly 1, it fails fallbacks 1: the chain-break gap after R:103-116 has no clear way. The fixture (asymmetric unit) passes
- Step 5b: ss_contacts_nearer (new, reported, not judged) is the fraction of touched single-stranded nucleotides drawn nearer their touching elements than the median element. 21 entries have touched single strands. It is 1.0 on 13 of them and below 1.0 on 1IHF 0.0, 2NZ4 0.2, 4PJO 0.273, 1URN 0.6, 1QTQ 0.667, 1MSW 0.667, 5W1H 0.828 and 30QM 0.929. These have not been examined one by one; they are follow-ups for Task 14, not failures
- Step 5b, 1URN: the loop R:6-12 held by U1A's RRM was drawn straight on, away from the sheet. A touched free end now tries first the turns (up to a right angle) that bring it nearest its touchers, so the loop lies toward the RRM (ss_contacts_nearer 0.2 -> 0.6; each loop nucleotide a β-strand touches is nearer that strand than the median element). The routed (reach-out) drawing of a long touched end still aims at free space
- Step 5b: contact-sheet panels hang from the top of their slots, so 1LMB, 1ZAA and 1BNA titles are in line (the line above about mid-sheet titles is resolved)

Remaining failures:

| Entry | Metric | Cause | Why it is accepted |
|---|---|---|---|
| 30QM | spearman 0.521, angle_median 35.7, angle_max 109.6, stems_nearer_frac 0.833 | The two protomers are drawn side by side (unrolled). P3 is held by both protomers and sits between them. | Ruled after step 3. The unrolled figure is the readable one, and the shape metrics measure the unrolling against the plain projection (M2). Held in the fixture test. |
| 3CRX | procrustes 0.817 | The four-Cre assembly is seen edge-on to the junction plane: arms 3, 4, 1 and 2 lie in a row in projection but are drawn as an X. | Ruled after step 1: the frame's choice. Users can use --rotate or --view. Held in the fixture test. |
| 3KTQ | contact_auc 0.582 | Klentaq grips the primer/template from both sides. With alpha19 out from between DNA 1a and 1b (C3), the gripping elements are drawn farther from the DNA. | C3 follow-up: no placement idea yet for a protein that encircles DNA. Firmer joints gave 0.371-0.555. |
| 4PJO | stems_nearer_frac 0.5 (new) | One snRNP is now drawn (assembly 1, was all four crystal copies). P2 has a single toucher, the helix C:4-8, which makes one residue–nucleotide contact. That helix is flipped across its own Sm sheet, away from P2, because C:76-79 took the space below; it ends 13.2 from P2 against a median of 9.5. This is placement of a one-contact element, not ring geometry. P1 passes: L:11-53 is at 0.8-1.3. | Accepted (controller review, step 5b). With two stems the fraction is either 0.5 or 1. A single contact weighs too little against the sheet's own packing. With four copies drawn the metric had been inflated by far non-touchers (asu: AUC 0.872, nearer 1.0). |
| 1YTB | contact_auc 0.544, stems_nearer_frac 0.5 (new) | One TBP-DNA complex is now drawn (assembly 1, was two). The frame looks down onto the saddle: TBP's sheet item projects at (0.0, 0.0), between DNA 1a (home y -2.0) and 1b (y 2.4), whose projected centres are 4.4 apart. In _attract the sheet's overlap term (189 on each piece at the final stage) beats the joint springs (2 x 4.0), so the halves end 16.5 apart (2.2 stem widths) with the sheet between them. DNA 1b's touchers A:102-126 end up 7.9-10.2 away. Step 5b tried closing the joint after _attract, sweeping the later piece back along the joint and stopping at the first new collision. The sheet lies in the corridor, so it stops at once (16.5 -> 14.7) and was not kept. 5B43 P2a/P2b behave the same (11.5 apart, held by A:888-931 between them) but pass. | Accepted with cause: the frame looks down onto the saddle. Follow-up idea: a view choice that keeps a protein–DNA interface side-on (the saddle seen edge-on), so the protein projects beside the DNA rather than onto it. |
| 4YBB | fallbacks 2, na_crossings 1, bead_overlaps 1 (new) | Dropping DA (C13) moves the placement. P46c's two pieces are drawn with their strands on swapped sides (the duplex turns over between them), so both joints AA:1037-1533 (2,402)->(2,403) and (2,424)->(2,425) must cross to the other rail. No elbow or bow clears the pieces, and one bow runs through P46c1. Separately, two runs (AA:959 and AA:976, held by AN and AS) pass within 0.42 of each other, and no slide gives a bead's width. | The same kind as step 4's P11b and the CROWDED fixture hold: placement on 84 cropped stems, with no routing rule left that clears it. NA x protein 54 -> 44. **Needs a ruling** (the earlier step-4 crossing is gone). |
