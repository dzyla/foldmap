# Random-100 fourth sample: review

**85 of 101 pass** (95% Wilson interval 75.8-90.0%). This is the count for the frozen release, v0.2.0 (e1392f7), on
a sample no one had seen: drawn after the release workflow passed on the tag, with the protocol's seed 20261025,
release-date cutoff 2026-10-13 and none of the 300 ids of samples 1-3 (`docs/validation/fourth-sample-protocol.md`),
and run once with every view and no fix (2026-10-28, zyla-office, `--jobs 28`, about 1,190 s wall). The 101 are the
100 drawn entries and 2CSE, a CA-only entry the sampler skipped in the S4 XL pool and runs once for A1c (it passes).
Every entry was measured by one code version (fingerprint 7879eea05bd9301c, git e1392f7, clean). No entry's run
stopped: no failed download, no view over its time limit. Each failure below is a measured outcome, not a run fault.

| stratum | entries | pass | | stratum | entries | pass |
|---|---|---|---|---|---|---|
| S1 X-ray, one chain | 12 | 12 | | S7 protein-DNA | 7 | 5 |
| S2 X-ray homo-oligomer | 12 | 11 | | S8 protein-RNA | 7 | 2 |
| S3 X-ray hetero-complex | 10 | 10 | | S9 RNA only | 3 | 2 |
| S4 cryo-EM | 15 | 13 | | S10 glycoproteins | 5 | 5 |
| S5 NMR | 8 | 7 | | S11 one protein, two states | 10 | 8 |
| S6 membrane | 8 | 6 | | S12 AlphaFold models | 3 | 3 |
| A1c CA-only (2CSE) | 1 | 1 | | | | |

Entries failing each criterion: A2b 1, A3 7, A4 1, A6 7, A15 1, A16 1, A17 9, A18 2, A19 2. A1, A2, A5, A7-A14 and
A20 pass on every entry that has something to judge (A1: no crash, no traceback, no Python warning in any view of any
entry). One defect can fail several criteria (one loop with no route fails A3, A6 and A18 on 7ZDM).

The 16 failures are listed below by check, as measured. Nothing was fixed or rerun (the protocol's rule: no fix after
the run). Causes were looked for by reading the run's measures and figures (on zyla-office) and, where a line says
so, by read-only probes that call the release's own functions on the entry's file (no figure drawn again, no code
changed; the probes are in `~/foldmap-paper/diag/r4probe/` on zyla-office, not committed). Each failure says how far
its cause is known:

- **established**: the failing value is traced to a rule or limit and the property of the input that meets it;
- **in part**: traced to a step of the drawing (a fallback taken, a naming scheme applied, a sheet drawn in two), but
  not why that step was reached;
- **not established**. A resemblance to a known case is named as a resemblance, not a cause.

## By check

### A17, nucleic-acid metrics (9 entries)

| id | stratum, bin | what failed |
|---|---|---|
| 7EDN | S9 S (X-ray 2.6 Å, 1 chain, 22 nt) | pair_f1 0.0 (16 pairs drawn, 0 reference pairs) |
| 8W2Q | S7 L (cryo-EM 3.06 Å, 5 chains) | pair_f1 0.741; stem_recall 0.667 |
| 9QMS | S7 XL (cryo-EM 5.0 Å, 30 chains) | pair_f1 0.892; stem_recall 0.714 |
| 9NL4 | S8 M (cryo-EM 4.6 Å, 5 chains) | pair_f1 0.865; stem_recall 0.667 |
| 8SSW | S8 M (X-ray 2.4 Å, 4 chains) | pair_f1 0.875 |
| 31FZ | S8 XL (cryo-EM 2.89 Å, 78 chains, 5,525 nt) | pair_f1 0.943; fallbacks 4; na_crossings 5; bead_overlaps 8 |
| 7SA4 | S8 XL (cryo-EM 2.55 Å, 58 chains, 4,724 nt) | stems_nearer_frac 0.994; fallbacks 1; na_crossings 3 |
| 2MF1 | S5 M (NMR, 7 chains, 72 nt) | procrustes 0.535 (limit 0.5) |
| 8U8V | S8 L (cryo-EM 2.74 Å, 6 chains) | procrustes 0.518 (limit 0.5) |

The pair causes come from a probe that ran foldmap's pair test (`nucleic.pair_nucleotides` and its rules) on each
pair the report lists as missed (it keeps 20 per entry) or extra.

- **7EDN: established. A limit of the measure, not of the figure.** The deposited model holds one 22-nt strand (the
  decoding-site RNA with Ag(I)). Assembly 1 adds its crystal-symmetry mate, and foldmap draws the duplex of the two
  copies (A1, A2: 16 pairs in 3 stems). A17's reference pairs come from RNApolis run on the deposited model
  (`scripts/na_validate.py`, `reference_pairs`: model 1 of the file). That model has no pair: a pair between a strand
  and its symmetry mate is not in the file. With no reference pair, every drawn pair counts as extra (precision 0),
  so pair_f1 is 0. The deposited-model measures (`deposited_pair_f1` 1.0) pass.
- **8W2Q: established** (the pairs foldmap's test leaves out: M1 of `docs/validation/na-fixes/report.md`, whose fix,
  F1, was deferred past the freeze). Of its 10 missed pairs, 5 have a hydrogen bond longer than 3.8 Å
  (`nucleic.HBOND_MAX`; N6-O4 4.06-5.31 Å), 3 have base planes 55-75° apart (`PLANE_MAX` 45°), and 2 (D9-E44, D24-E29)
  pass the test but each nucleotide keeps a better partner one step along (D9-E43, D24-E30). Those two and D10-E42
  are 3 of the 4 extra pairs. 1 of the 3 reference stems is not recalled.
- **9QMS: established** (M1). 20 missed pairs: 19 have a bond longer than 3.8 Å (4.0-7.1 Å; some also have planes
  over 45° or C1'-C1' under 9 Å), and 1 (i99-j56) loses to a better partner (i99-j55). 6 of 21 reference stems are not
  recalled.
- **9NL4: established.** Of its 4 missed pairs, 3 (B46-T25, B47-T24, B48-T23) pass foldmap's test, but foldmap pairs
  the same nucleotides one step along (B45-T26, B47-T25, B48-T24). Both registers pass on this 4.6 Å model, and each
  nucleotide keeps its best partner (`_mutual_best`). The 4th (B55-T17) has N6-O4 4.36 Å. All 6 extra pairs pass the
  test; on two runs (B45-B48 against T23-T26, B54-B56 against T15-T17) they lie one step off RNApolis' pairs.
- **8SSW: established.** Its 4 missed pairs (C6-D21, C9-D18, C18-D9, C21-D6) are G·U pairs that RNApolis calls wobble
  (Saenger XXVIII). Their wobble hydrogen-bond atoms (G N1-U O2, G O6-U N3) are 4.09-5.52 Å apart, past 3.8 Å.
- **31FZ: in part.** pair_f1 is established (M1): 110 of 1,242 reference pairs are missed. Of the 20 listed, 17 have a
  bond longer than 3.8 Å that the stretched second pass does not take either, 1 (BC:1115-1120) has C1'-C1' 12.01 Å
  (`C1_PAIR` 9-12 Å), and 2 (BC:374-392, BC:375-391) lose to a better partner (BC:375-392). The fallbacks are the
  nucleic-acid strands of A3 below. The crossings and bead overlaps are not established.
- 7SA4: not established. Its pairs pass (pair_f1 0.993). One stem is not drawn nearer its protein than the rest
  (stems_nearer_frac 0.994). 1 strand is a fallback (A3). 3 nucleic-acid crossings. 31FZ and 7SA4 resemble 8XJV on
  the third sample (a large cryo-EM nucleic-acid complex with fallbacks and crossings). That is a resemblance.
- 2MF1, 8U8V: not established. Both draw every reference pair and stem (pair_f1 1.0, stem_recall 1.0). Their
  Procrustes values are 0.035 and 0.018 over the limit, on 4 and 3 stems.

### A3, route fallbacks (7 entries)

| id | stratum, bin | views | fallbacks (the report names up to 4 per view) |
|---|---|---|---|
| 5ZJU | S2 XL (X-ray, 60 chains, a virus-like particle) | focus | 17 loops: C:157-161>C:164-165, F:173-174>F:179-181, G:157-161>G:164-165, I:173-174>I:185-186 ... |
| 8C5C | S4 XL (cryo-EM, 52 chains, 22,620 residues) | domains | 25 loops: c:198-202>c:205-214, c:250-255>c:265-271, c:265-271>c:286-294, c:305-307>c:310-318 ... |
| 8WLE | S4 XL (cryo-EM, 52 chains, the flagellar LP ring, C13) | focus | 67 loops: A:186-188>A:194-196, A:197-199>A:204-210, B:186-188>B:194-196, C:186-188>C:194-196 ... |
| 7ZDM | S6 XL (cryo-EM, 45 chains, complex I) | membrane | loop d:192-194>d:195-197 |
| 8W2Q | S7 L | domains | nucleic-acid strand E:21-22 |
| 7SA4 | S8 XL | default, deposited, beta, focus | nucleic-acid strand 1:836-836 |
| 31FZ | S8 XL | default, deposited, beta, focus | strands BC:867-876, CI:1430-1431, CI:3720-3730; gap after CI:3599 |

- 8WLE focus, 8C5C domains: **in part.** Their route record is that of a spent page budget. This is the profile of
  9FZL on the third sample. 8WLE: 7,849,719 of the page's 8,000,000 second-search nodes spent, 151 searches refused,
  67 loops drawn as curves. 8C5C: all 8,000,000 spent, 49 refused, 25 curves. The record does not say which loops
  were refused, so which curves a wider search would have routed is not established.
- 7ZDM membrane: **in part.** The loop found no route. It was drawn by the membrane view's walled fallback. The notes
  count loop_wall_detour 1 and loop_wall_unkept 1: no detour kept to the loop's wall, so the last one tried is drawn,
  a curve that runs round the page. This curve causes the 2 A6 collisions and 1 of the 2 A18 findings below. Why no
  route was found is not established.
- 5ZJU focus: not established. No second search was run or refused (0 nodes spent), so the budget is not the cause.
  The stuck loops join short strands (157-161 to 164-165; 173-174 to 179-186) of the whole cage, drawn at once.
- 8W2Q, 7SA4, 31FZ: nucleic-acid strands and a gap drawn as plain curves (the notes: na_run_unclear, na_gap_unclear,
  na_tail_straight). Not established.

5ZJU, 8C5C, 8WLE and 7ZDM fail only in a view that draws every chain on a page held at 48 in. foldmap's own notes put
each of those views below the readability floor (text 3.0-4.0 pt) and name the options that draw less. The default
view of each passes. 5ZJU and 8WLE draw one subunit there; 8C5C and 7ZDM draw all chains.

### A6, text collisions (7 entries)

| id | views | collisions |
|---|---|---|
| 8WLE | focus | 29: 9 with loops drawn as curves (A3), 17 with loops routed over labels (A4), 3 with routed loops (label:D:111-113 × loop:F:124-155>F:161-171; label:D:118-120 × loop:D:111-113>D:118-120; label:T:221-229 × loop:T:186-188>T:194-196) |
| 8C5C | domains | 5: domain names D2 · c, D1 · i, D2 · S and chain k's N end, each on a loop drawn as a curve (A3); label:Y:75-77 × loop:Y:382-398>Y:405-409 |
| 7ZDM | membrane | 2: labels 1:231-237 and d:207-219 under the curve of loop d:192-194>d:195-197 (A3) |
| 8W2Q | default, deposited, beta; domains | D:12's base letter on chain A's C-end stub; domains: D:32's base letter on run 5 |
| 9QMS | default, deposited, beta, focus; domains | tie 26-127's label × chain C's N-end stub; chain F's C end × tie 58-96's mark; domains: I:9-119's 5′ end × tie 9-148 |
| 7SA4 | default, deposited, beta, focus | gap 513's line × chain x's C end and its stub; chain T's N end × run 93 |
| 31FZ | default, deposited, beta, focus; domains | tie 928-1015's label × chain Bf's N end; stem 199's first-pair number × run 396; chain AS's N end × tie 415-19's mark; domains 5 (also gap 561 × stem 105, two tie labels on a tie and a run, CJ's 3′ end × chain AD's stub) |

- 8WLE, 8C5C, 7ZDM: **in part.** The collisions with curves and with loops routed over labels are their A3 and A4
  loops. 26 of 8WLE's 29 are these, 4 of 8C5C's 5, and both of 7ZDM's. The other 4 (8WLE 3, 8C5C 1) are not
  established.
- 8W2Q, 9QMS, 7SA4, 31FZ: not established. On 8W2Q the letter T of D:12 (x 608.6 pt) sits on chain A's C-end stub
  (x 610.5 pt, y 684.9-697.3 pt). Why the letter's spot was not kept clear of the stub is not traced. Tie labels and
  marks on other marks (9QMS, 31FZ) are the kind of collision the third sample showed on 3G6E, 6QIK and 7F5S. That is
  a resemblance.

### A18, membrane (2 entries)

- **5V5Z** (S6 M, X-ray, 1 chain, CYP51 of *Candida albicans*): no band drawn. The membrane view found no membrane
  ("no membrane found ...; drawing it without one"). **Established** (a probe of `membrane.find_membrane`'s steps on
  the file). The file has no OPM dummy atoms. The model starts at residue 25, and its N-terminal anchor helix is
  A:28-40, 13 residues. The estimate reads only helices of 16 residues or more (`membrane._long`). The six it reads
  (17-32 residues, from A:152 on) belong to the soluble P450 domain. No slab across them is membrane, so the estimate
  returns none. A single-pass protein whose membrane helix is modelled only in part gets no band.
- **7ZDM** (S6 XL, cryo-EM, 45 chains, ovine complex I), membrane view, two findings. **In part.** (1) The TM helix
  A:56-72 (17 residues) is drawn wholly inside the band and covers 0.80 of it, reaching neither face (seen on the
  figure). The view moves a TM element drawn inside the band to a face (`layout._to_face`). Why this one stayed
  inside is not established. (2) Loop d:192-194>d:195-197 is drawn 1.7 Å nearer the mid-plane than its residues.
  This is the walled-fallback curve of A3.

### A19, two states named alike (1 pair, 2 entries)

- 9S44/5DY5 (P5, monomer/oligomer, CA RMSD 1.74 Å): names compared on all 21 pairs, 8 differ. **In part.** A19 reads
  the beta views. 9S44's chain has a βα row there, so its helices are named after the strand before them (αN, α1′,
  α1″, α1, α2, α5, α6, α7, α9, α9′; `names.beta_names`). 5DY5 is drawn as its 2-chain assembly (chain A1) and has no
  row, so its helices keep their numbers in sequence (α1-α10). The same helices carry names from two schemes. Why the
  row is found in one state and not in the other is not established.

### A2b, the app's hosted limits (1 entry)

- 31FZ (S8 XL, 78 chains, 16,665 residues and 5,525 nt): the app's job was stopped at its 90 s limit (90.11 s), with
  foldmap's message (draw part of it, or locally). Its default view took 156 s (A2's limit for 20 chains or more is
  180 s, so A2 passes). **In part.** Its app view ran last. The last other entry, 7SA4, had finished at 22:42:54;
  31FZ finished at 22:48:36, 91 s after its app view began (file times on zyla-office). So no other job slowed it: the
  drawing alone takes longer than the hosted 90 s on zyla-office. Where the time goes is not traced. The third
  sample's 8XJV failed A2b too, under load.

### A4, a loop over a label (1 entry)

- 8WLE (focus): A:194-196>A:197-199, B:124-155>B:161-171, B:194-196>B:197-199, C:124-155>C:161-171 (17 in the view;
  the notes count loop_over_text 17). **In part.** These loops took the fallback `route._Router.route` takes when the
  only route found runs over a label. They are on the page whose budget was spent (A3).

### A15, copies drawn alike (1 entry)

- 8WLE (focus: the whole C13 ring, 52 chains): worst tie Y~A 83.06 (limit 0.5). **In part.** The 26 upper-case
  chains' strands are drawn in two sheets. One holds 128 strands (chains A-V, the long sheet at the foot of the page).
  The other holds 28 (chains V-Z, the upright sheet at the top). Chain V has strands in both (read from the SVG's
  sheet panels). A15 leaves out the strands of a sheet the two copies share. So the 12 ties whose copies lie in one
  drawn sheet (C, E ... U against A; V and X against Z) compare 6 elements each and pass (deviation 0). The 12 ties
  whose copies lie in different drawn sheets (B, D, F, H, J, L, N, P, R, T against Z; W and Y against A) compare all
  12 elements, strands included. They deviate 54-83 page units. The lower-case chains (a-z, 23 elements each) are
  drawn alike. Why the strands are drawn in two sheets cut at chain V, and why ties cross that cut, is not
  established.

### A16, domains (1 entry)

- 8C5C (domains): 3 domain names (D2 · c, D1 · i, D2 · S) on loops c:198-202>c:205-214, S:405-409>S:416-435 and
  c:396-399>c:405-426. Each is a curve of its A3 (**in part**, as A3). A domains-view name crossed by a loop on a page
  whose second-search budget is spent is 9FZL's profile on the third sample. That is a resemblance.

## How far each cause is known

Counted per failed criterion of an entry (31 in all):

| criterion | failed | established | in part | not established |
|---|---|---|---|---|
| A17 | 9 | 5 (7EDN, 8W2Q, 9QMS, 9NL4, 8SSW) | 1 (31FZ) | 3 (7SA4, 2MF1, 8U8V) |
| A3 | 7 | 0 | 3 (8WLE, 8C5C, 7ZDM) | 4 (5ZJU, 8W2Q, 7SA4, 31FZ) |
| A6 | 7 | 0 | 3 (8WLE, 8C5C, 7ZDM) | 4 (8W2Q, 9QMS, 7SA4, 31FZ) |
| A18 | 2 | 1 (5V5Z) | 1 (7ZDM) | 0 |
| A19 | 2 | 0 | 2 (9S44, 5DY5) | 0 |
| A2b, A4, A15, A16 | 4 | 0 | 4 (31FZ; 8WLE; 8WLE; 8C5C) | 0 |
| all | 31 | 6 | 14 | 11 |

By entry: 4 of the 16 failing entries have every failure established (7EDN, 9NL4, 8SSW, 5V5Z). 4 have none traced
(5ZJU, 2MF1, 8U8V, 7SA4). The other 8 are traced in part.

Not fixed, for the next release to weigh:

- A17's reference reads only the deposited model, so a duplex that a symmetry mate completes has no reference pairs
  (7EDN: a figure that draws the duplex fails with pair_f1 0.0). This is a fault of the measure
  (`scripts/na_validate.py`), not of the figure.
- 8WLE's whole ring ties copies across the cut between its two drawn sheets (A15 54-83), and S3 finds a fault in that
  view. It is a candidate layout fault, cause not traced.
- 7ZDM's TM helix A:56-72 is left inside the band, which `layout._to_face` exists to prevent. It is a candidate
  fault, cause not traced.
- The pair test's 3.8 Å bond limit (M1, F1 deferred) causes 4 of the 9 A17 failures (8W2Q, 9QMS, 31FZ, 8SSW). One of
  them is 8SSW's G·U wobbles in a 2.4 Å crystal structure.
- The membrane estimate's 16-residue minimum leaves a single-pass protein with a part-modelled anchor helix (5V5Z)
  with no band. This is a limit of the estimate.

## Science checks and S5 (reported, not gated)

S1, S2, S4 and S6 find nothing in any view. S3 finds one fault, in 8WLE's focus view: sheet 86 has 12 barrel strands
drawn, but its pairing cycle has 17. It is not investigated here. S5 compares the automatic domains with CATH and is
not counted in the 85. 215 chains were scored in 32 entries. One-domain chains: 93 pass, 63 fail (the automatic
domains split a chain CATH leaves whole). Chains with two or more CATH domains: 38 pass, 21 fail (mean boundary F1
0.78, recall 0.80). At entry level, S5 passes on 16 entries and fails on 16.

## Time

The run did not record its wall time. The file times on zyla-office put it at about 1,190 s (first entry folder
22:29:32, report written 22:49:24, -0600) at `--jobs 28` on 32 cores. The slowest entry is 31FZ: 1,066 s over its
ten views (default view 156 s; its app view stopped at 90 s). Next come 7SA4 (726 s; default 107 s), 8C5C (596 s;
domains view 208 s) and 7YF0 (551 s). A2 passes on every entry. The slowest default view is 31FZ's, 156 s against
180 s.
