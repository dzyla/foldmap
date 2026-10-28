# Plan 29: the look test

Written by scripts/app29_after.py sheets at git b84ebc3 (clean): every entry whose app SVG changed between the before commit rerun and the after run (both halves, gpu-desktop), the leads, and the named crops, each drawn as the app's PNG before beside after at one scale (`px/pt`: pixels per point of the page). Before: before (f538981); after: after (162e067). The sheets are under docs/superpowers/work/app29/after/sheets/ (not committed; the sha256 names each). The verdict (better, same or worse) and its one-line reason are the controller's, written by hand after looking (plan 29 Task N+1 Step 3); any "worse" is reverted or explained before main.

| # | pair | half | why listed | outcome before -> after | px/pt | sheet sha256 | before SVG | after SVG | verdict | reason |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 26EL | tune | its SVG changed; lead: Task 8b named entry: focus auto draws one copy (A alone) | pass -> pass | 1.1315 | da009fff03d2c699b30e7992cfa54c1e95037a9340f02bab5f308b3441b3d89a | 06964d4e7ba9 | ff51eedcb967 |  |  |
| 2 | 29TB | tune | its SVG changed | pass -> pass | 3.6202 | 91c39eb34bb62b0ee052414fc591d9d5b216ca14624b4670e991972db60e9022 | 68540b113d37 | 9ec3c71fe08d |  |  |
| 3 | 31FD | tune | its SVG changed; lead: Task 5 named entry: a two-copy complex under the floor, now one subunit | fail floor -> pass | 0.463 | ff7472d4e6be8f2a0a3c1be0af252e94affbe674d61d6cd12418ab780eb3d133 | e9f0b5fd6cfd | c1d847a64133 |  |  |
| 4 | 9SVL | check | its SVG changed | pass -> pass | 3.6195 | c894c5d16381589b918bb7d04fb836f6613afa2f434af022951982bcafe31636 | cbc5a7235170 | dcdfab4be9bd |  |  |
| 5 | 9SXF | tune | its SVG changed | pass -> pass | 3.0864 | 1a17e127871415e7c19bb6db0a6d7b5f2cd9b9695a27b174540b9825606991b3 | 8896426b730f | 7ae2b1157205 |  |  |
| 6 | 9WNE | tune | its SVG changed; lead: Task 4: the Mg moved off the loop crossing | pass -> pass | 1.654 | 4fb29ec6cda4df4a1c89b9b72120338d9b7da828d48ef06e7fd69f072a4ddbaa | c40bb7fc7732 | 138613e329e5 |  |  |
| 7 | 9X4L | tune | its SVG changed; lead: Task 4 named entry (A20) and Task 7 (S3); a giant assembly | fail A20 -> pass | 0.4805 | 72816b6b00ab376dd9af46133f872d2ab626857c3542fa5a1c686b267be53f74 | 40184e89b350 | 606eecd5cd80 |  |  |
| 8 | 9X4M | check | its SVG changed | fail floor -> pass | 0.463 | dc81f783d16786a780ac985df7a7db1339759821175e5879f04424128c04023d | b28990077787 | b775282bc7b7 |  |  |
| 9 | 9YN3 | check | its SVG changed | pass -> pass | 0.8445 | 8b02d86c5a8a39b3bf71d1eb34e86886f144c0b30d6c596c2f0a95196576468b | 731c565a8191 | d4462c94c7e6 |  |  |
| 10 | 6JLN | lead | lead: Task 5 lead: a two-copy complex under the floor (plan 28's six sample defaults) | - | 0.463 | 379798ade9b7c38f0f7616022790213337eeeccc5432a27679855c8dd8d1c926 | d76c0330d6f9 | 35bfa63d8975 |  |  |
| 11 | 9FZL | lead | lead: Task 5 lead: under the floor (plan 28's six sample defaults; third sample) (zyla-office: over the hosted job limit on the run's host) | - | 0.463 | cd5b08b9852eaec176a139c2cdcfe55e7b7c7c33242f2f763b3ebfd6ec09198a | 2fbb74c7de72 | 3da210018e74 |  |  |
| 12 | 7SPI | lead | lead: Task 8b: TraV copies modelled in part, one protein (third sample) | - | 2.7019 | dc6e2ef8ee722a9d8964521e1e571ce522d57a9973cc964eaa08d20128587c9f | 713c1ced9e66 | cb3a48c69eff |  |  |
| 13 | 1KX5 | lead | lead: Task 4: the nucleosome's Mn mark moved (NA control) | - | 1.2583 | 73a12e689f36b317d8dcbff137212c1de1905d605d147f9a52619c4980ec606f | 4e72b11c3e59 | d1957b3409aa |  |  |
| 14 | 9WNE-crop | crop | P29-10: the A1301/A1302 Mg tethers lying along each other (Task 4), cropped to the two marks | pass -> pass | 4.1661 | 2f95f318e17547ef57197fd107a7a9394399a65e2b418807cba45dc124998246 | c40bb7fc7732 | 138613e329e5 |  |  |
