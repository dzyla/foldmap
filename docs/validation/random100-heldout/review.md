# Random-100 held-out sample: review

77 of 100 pass, against 92 of 100 on the tuned sample (`docs/validation/random100/`), judged by the same criteria and limits, at the same commit (5650dbd), with every view and no fix. The failures below are listed with the checks they fail, as measured; none is investigated or fixed here: they are the next plan's input.

| id | stratum | bin | fails | what failed (first view's detail, cut) |
|---|---|---|---|---|
| 2FR2 | S1 | S | A3 | A3: default: A:89-100>A:103-109, A:115-118>A:125-135; deposited: A:89-100>A:103-109, A:115-118>A:125-135; beta: A:89-100>A:103-109, A:115-118>A:125-135 |
| 2NAD | S2 | M | A3, A5, A6 | A3: default: B:250-253>B:267-270, B:267-270>B:277-281; deposited: B:250-253>B:267-270, B:267-270>B:277-281; beta: B:250-253>B:267-270, B:267-270>B:277-281; domains: B:250; A5: default 2 (['B:250-253', 'B:267-270']; ['B:267-270', 'B:277-281']); deposited 2 (['B:250-253', 'B:267-270']; ['B:267-270', 'B:277-281']); beta 2 (['B:250-253', 'B:267; A6: default 2 (label:B:267-270 × label:B:277-281; label:B:267-270 × strand:B:277-281); deposited 2 (label:B:267-270 × label:B:277-281; label:B:267-270 × strand:B:277-281) |
| 8BCU | S4 | L | A3 | A3: default: G:92-99>G:111-114, G:119-122>G:135-140, G:150-158>G:166-173, G:193-200>G:209-215; deposited: G:92-99>G:111-114, G:119-122>G:135-140, G:150-158>G:166-173, G:1 |
| 8D3D | S4 | XL | A1, A5, A6, A7, A8, A15 | A1: focus: killed after 900 s; A5: focus not measured (not measured); A6: focus not measured (not measured); A7: focus not measured (not measured); A8: focus not measured (not measured); A15: focus not measured (not measured) |
| 12BM | S4 | XL | A3, A6, A16 | A3: default: A:93-101>A:111-115, B:98-101>B:107-116, C:21-66>C:71-79, C:84-89>C:94-98; deposited: A:93-101>A:111-115, B:98-101>B:107-116, C:21-66>C:71-79, C:84-89>C:94-98; A6: domains 1 (domain-label:D1 · B × loop:B:98-101>B:107-116); A16: domains 1 () |
| 9GTP | S4 | XL | A15 | A15: focus 52.108405597395176 |
| 7BKC | S4 | XL | A3, A4, A6, A16 | A3: domains: E:15-20>E:26-37, E:188-192>E:203-206, C:27-37>C:41-43, C:49-53>C:57-59; A4: domains: D:183-195>D:198-201, d:183-195>d:198-201; A6: domains 70 (domain-label:D2 · A × loop:A:351-353>A:365-367; domain-label:D2 · A × loop:A:368-375>A:396-403; domain-label:D2 · D × loop:D:170-174>D:183-195); A16: domains 6 () |
| 6MJ2 | S6 | L | A18 | A18: no band drawn |
| 7JI3 | S6 | M | A18 | A18: loop in the slab C:101-108>C:116-119 |
| 7VVK | S6 | M | A18 | A18: no band drawn |
| 6ZKT | S6 | XL | A4, A6 | A4: domains: i:3-16>i:19-30; A6: domains 1 (label:i:3-16 × loop:i:3-16>i:19-30) |
| 6L9Z | S7 | L | A6 | A6: domains 1 (terminus:N:A × na:run:79) |
| 5K7Z | S7 | M | A3, A15, A17 | A3: default: strand F:1-16; deposited: strand F:1-16; beta: strand F:1-16; domains: strand F:1-16; A15: default 1.106238016241444; A17: pair_f1 0.182 (limit 0.95); stem_recall 0.0 (limit 0.85); fallbacks 1 (limit 0) |
| 2FD8 | S7 | S | A6 | A6: default 5 (na:nt:B:501 × strand:A:116-122; na:nt:B:502 × strand:A:204-210; na:nt:B:503 × ligand:AKG:A400); deposited 5 (na:nt:B:501 × strand:A:116-122; na:nt:B:502 ×  |
| 7NZ2 | S7 | XL | A3, A4, A6, A16, A17 | A3: domains: B1:211-221>B1:227-408, B3:211-221>B3:227-408, B3:1051-1123>B3:1136-1140, F3:182-186>F3:197-205; A4: domains: A1:1180-1196>A1:1208-1262, A4:217-221>A4:226-400, A4:986-1045>A4:1051-1124, B1:934-971>B1:983-1044; A6: domains 9 (domain-label:D7 · A1 × loop:B3:211-221>B3:227-408; domain-label:D12 · A1 × loop:A1:1180-1196>A1:1208-1262; domain-label:D1 · A3 × loop:F3:182-186>F3:197-20; A16: domains 6 (); A17: spearman 0.446 (limit 0.6) |
| 9K4F | S8 | L | A17 | A17: stems_nearer_frac 0.875 (limit 1.0) |
| 13FZ | S8 | L | A3, A6, A15, A17 | A3: default: strand B:1-13; deposited: strand B:1-13; beta: strand B:1-13; A6: default 3 (na:nt:B:2 × loop:A:1013-1019>A:1027-1031; na:nt:B:3 × strand:A:614-618; na:nt:B:8 × loop:A:1003-1006>A:1013-1019); deposited 3 (na:nt:B:2 × loop:A:1013-101; A15: default 17.0586549756951; A17: fallbacks 1 (limit 0) |
| 6S13 | S8 | XL | A3, A6, A17 | A3: domains: strand A:2729-2729, strand A:2732-2733, strand A:2737-2738; A6: default 1 (na:tie:448-506:label × na:stem:166); deposited 1 (na:tie:448-506:label × na:stem:166); beta 1 (na:tie:448-506:label × na:stem:166); focus 1 (na:tie:448-506; A17: pair_f1 0.944 (limit 0.95); beads_missing 134 (limit 0); na_crossings 2 (limit 0); bead_overlaps 7 (limit 0) |
| 8ZTU | S8 | XL | A3, A6, A16, A17 | A3: domains: strand a:2701-2702, strand a:2705-2705, strand a:2710-2711; A6: default 2 (na:gap:499 × na:run:477; ligand-label:MG:a6201 × na:stem:188); deposited 2 (na:gap:499 × na:run:477; ligand-label:MG:a6201 × na:stem:188); beta 2 (na:gap:4; A16: domains 1 (); A17: beads_missing 1 (limit 0); bead_overlaps 2 (limit 0); gap_label_collisions 1 (limit 0) |
| 8U9T | S11 | M | A19 | A19: names compared on all 3 pairs; only 3 of 24 and 11 elements paired; A1:11-14 unpaired; A1:17-22 unpaired; A1:35-39 unpaired; A1:48-50 unpaired; A1:54-59 unpaired; A1 |
| 7IT6 | S11 | S | A19 | A19: names compared on all 4 pairs; only 4 of 11 and 24 elements paired; A:35-48 unpaired; A:61-63 unpaired; A:68-73 unpaired; A:80-82 unpaired; A:95-97 unpaired; A:103-1 |
| 2BBJ | S11 | L | A15, A19 | A15: default 1.4176624736566856; A19: names compared on all 13 pairs, 1 differ; η2~η2 (A:52-55~A:51-55); A:193-201 unpaired; A:244-273 unpaired; A:275-299 unpaired; A:302-310 unpaired; A:330-344 unpaired |
| 2BBH | S11 | S | A19 | A19: names compared on all 13 pairs, 1 differ; η2~η2 (A:51-55~A:52-55); A:47-49 unpaired; A:161-163 unpaired; A:193-201 unpaired; A:244-273 unpaired; A:275-299 unpaired;  |

Entries failing each criterion: A1 1, A3 10, A4 3, A5 2, A6 11, A7 1, A8 1, A15 5, A16 4, A17 6, A18 3, A19 4.

Notes on what the numbers measure:

- 8D3D's focus view was killed after 900 s (A1), so its A5-A8 and A15 read not measured there, which fails them closed.
- The two-state pair P3 (8U9T/7IT6) pairs only 3 and 4 of its elements: A19 reads each entry's first protein chain, and the pair's first chains may not be the shared protein (not checked here).
- The pair types are re-derived in `report.md` ("Two-state pairs"), as for the tuned sample.
