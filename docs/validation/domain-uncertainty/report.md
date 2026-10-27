# The automatic domains with their uncertainty

Plan 28 Task 5 (the review's 3.4, rulings Q42, Q43). Written by: foldmap 0.1.0, git 42e3fe80e84b (asserted by runner, dirty False), code fingerprint 4eb4cad74365b1ed, scripts/domain_uncertainty.py sha256 7bf6bc844b98; Python 3.13.13 on Linux-7.0.0-38-generic-x86_64-with-glibc2.43; gemmi 0.7.5, numpy 2.5.4, scipy 1.18.1, networkx 3.7, matplotlib 3.11.2, PyYAML 6.0.3

Sets, tagged as scripts/paper/sets.json tags them: dom_cal, 22 chains (tuning); dom_held, 11 chains (regression); dom_fresh, 11 chains (regression); r100_third, 109 chains (unseen at 7696947, historical; looked at again here). dom_cal and dom_held are the reference chains split chain by chain (seed 20261016:domains), as the parameters were chosen. A chain meets the targets (S5) when its overlap, boundary F1 and recall are each at least 0.8, 0.7 and 0.85 against CATH (domain_score, boundaries within 20 residues). DOMAIN_NCUT is 0.37. report-chains.json holds each chain's cuts, trace and scores at every t.

Domains against features.py at b72f65d (sha256 1596359ddc81): 153 chains compared (every chain of every set), 0 differ.

## 1. Why the multi-domain chains of the third sample fail (read-only)

97 PDB entries in the sample; 35 with a chain CATH chops scored; 109 chains scored, 47 of them with two CATH domains or more: 12 meet the targets, 35 fail, in 8 entries (copies of one protein in an entry count once each: 1GQ2 has 16 chains). Single-domain chains failing (an extra split): 7 (1RKA A, 3WAU A, 3WAU B, 1XZ4 B, 1XZ4 D, 3W96 D, 3W96 H).

What went wrong (exclusive, by the automatic and CATH domain counts):

| class | chains (entries) | what it means | the mechanism in features.py |
|---|---|---|---|
| boundary missed | 15 (5) | fewer automatic domains than CATH's | the part's lowest cut at DOMAIN_NCUT or above (residue_domains, _best_cut), or no cut keeps DOMAIN_MIN_RESIDUES / DOMAIN_MIN_SEGMENT, or a cut made and undone by a step of _chain_domains (traced below) |
| boundary moved | 5 (3) | as many, with a boundary elsewhere | the lowest normalized cut is not CATH's boundary (the criterion prefers another), a discontinuous CATH domain one cut cannot take (two runs at most, _lowest_cut), or the element rules moved a right residue cut |
| extra split | 15 (2) | more automatic domains than CATH's | a cut below DOMAIN_NCUT inside one CATH domain (a lobe or subdomain CATH keeps whole) |

Missed boundaries whose residue cut was made: 8 chains in 2 entries, 7 of them meeting the targets at the residue level; the first step of _chain_domains that brings the count below CATH's: fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) 7; _no_strays (src/foldmap/features.py:650) 1.

What the failing chains are like (tags, not exclusive): discontinuous domain 7, tightly packed 11, large 14, none of them 17. Discontinuous: a CATH domain in two runs or more along the chain; tightly packed: a CATH domain whose normalized cut against the rest of the chain (the contacts residue_domains weighs) is at DOMAIN_NCUT or above; large: three CATH domains or more, or more residues than DOMAIN_MAX_POINTS. Their CATH boundaries: found 35, moved 5, missed 23, extra 75 (moved: none within the tolerance, but an automatic boundary no CATH one claims lies between its two neighbours; extra: automatic boundaries no CATH one claims, left over). 7 of the failures meet the targets at the residue level, before the element rules.

| id | chain | residues | CATH | auto | overlap | F1 | recall | class | tags | boundaries f/m/x/e | CATH domains' ncut | mechanism |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2AIQ | A | 231 | 2 | 1 | 0.52 | 0.00 | 0.00 | boundary missed | discontinuous domain, tightly packed | 0/0/3/0 | 0.41, 0.41 | not cut: a CATH domain's own normalized cut is 0.41, at or above DOMAIN_NCUT 0.37 |
| 3ABR | A | 453 | 3 | 2 | 0.85 | 0.67 | 0.50 | boundary missed | tightly packed, large | 1/0/1/0 | 0.41, 0.49, 0.35 | not cut: a CATH domain's own normalized cut is 0.49, at or above DOMAIN_NCUT 0.37 |
| 3ABR | C | 453 | 3 | 2 | 0.53 | 0.40 | 0.50 | boundary missed | tightly packed, large | 1/0/1/2 | 0.41, 0.49, 0.35 | not cut: a CATH domain's own normalized cut is 0.49, at or above DOMAIN_NCUT 0.37 |
| 3EOM | A | 387 | 3 | 2 | 0.69 | 0.50 | 0.50 | boundary missed | tightly packed, large | 1/1/0/0 | 0.50, 0.34, 0.36 | not cut: a CATH domain's own normalized cut is 0.50, at or above DOMAIN_NCUT 0.37 |
| 3EOM | B | 378 | 3 | 2 | 0.70 | 0.50 | 0.50 | boundary missed | tightly packed, large | 1/1/0/0 | 0.51, 0.34, 0.37 | not cut: a CATH domain's own normalized cut is 0.51, at or above DOMAIN_NCUT 0.37 |
| 3EOM | D | 381 | 3 | 2 | 0.70 | 0.50 | 0.50 | boundary missed | tightly packed, large | 1/1/0/0 | 0.51, 0.33, 0.35 | not cut: a CATH domain's own normalized cut is 0.51, at or above DOMAIN_NCUT 0.37 |
| 3G6E | A | 237 | 2 | 1 | 0.50 | 0.00 | 0.00 | boundary missed | tightly packed | 0/0/1/0 | 0.45, 0.45 | the residue cut was made (2 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 2 |
| 3G6E | B | 337 | 3 | 2 | 0.73 | 0.89 | 0.80 | boundary missed | discontinuous domain, tightly packed, large | 4/0/1/0 | 0.46, 0.31, 0.28 | the residue cut was made (3 residue domains); _no_strays (src/foldmap/features.py:650) brought the count below CATH's 3 |
| 3G6E | L | 145 | 2 | 1 | 0.65 | 0.00 | 0.00 | boundary missed | - | 0/0/1/0 | 0.00, 0.00 | the residue cut was made (2 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 2 |
| 3G6E | P | 143 | 3 | 2 | 0.77 | 0.67 | 0.50 | boundary missed | large | 1/0/1/0 | 0.00, 0.12, 0.04 | the residue cut was made (3 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 3 |
| 3G6E | W | 154 | 2 | 1 | 0.64 | 0.00 | 0.00 | boundary missed | discontinuous domain, tightly packed | 0/0/2/0 | 0.38, 0.38 | not cut: a CATH domain's own normalized cut is 0.38, at or above DOMAIN_NCUT 0.37 |
| 4YJ5 | A | 518 | 3 | 2 | 0.74 | 0.40 | 0.25 | boundary missed | discontinuous domain, large | 1/0/3/0 | 0.20, 0.23, 0.09 | the residue cut was made (4 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 3 |
| 4YJ5 | B | 518 | 3 | 2 | 0.74 | 0.40 | 0.25 | boundary missed | discontinuous domain, large | 1/0/3/0 | 0.21, 0.23, 0.09 | the residue cut was made (4 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 3 |
| 4YJ5 | C | 518 | 3 | 2 | 0.74 | 0.40 | 0.25 | boundary missed | discontinuous domain, large | 1/0/3/0 | 0.21, 0.23, 0.09 | the residue cut was made (4 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 3 |
| 4YJ5 | D | 518 | 3 | 2 | 0.74 | 0.40 | 0.25 | boundary missed | discontinuous domain, large | 1/0/3/0 | 0.20, 0.25, 0.11 | the residue cut was made (4 residue domains); fold (DOMAIN_MIN_SSE_RESIDUES) (src/foldmap/features.py:652) brought the count below CATH's 3 |
| 1GQ2 | A | 555 | 2 | 2 | 0.88 | 0.67 | 1.00 | boundary moved | - | 1/0/0/1 | 0.26, 0.26 | the lowest cut lies elsewhere: foldmap's cut at normalized cut 0.18, CATH's domains at 0.26, 0.26 (the criterion prefers another boundary) |
| 1GQ2 | C | 555 | 2 | 2 | 0.88 | 0.67 | 1.00 | boundary moved | - | 1/0/0/1 | 0.26, 0.26 | the lowest cut lies elsewhere: foldmap's cut at normalized cut 0.18, CATH's domains at 0.26, 0.26 (the criterion prefers another boundary) |
| 1GQ2 | L | 555 | 2 | 2 | 0.88 | 0.67 | 1.00 | boundary moved | - | 1/0/0/1 | 0.26, 0.26 | the lowest cut lies elsewhere: foldmap's cut at normalized cut 0.18, CATH's domains at 0.26, 0.26 (the criterion prefers another boundary) |
| 1NID | A | 333 | 2 | 2 | 0.83 | 0.00 | 0.00 | boundary moved | tightly packed | 0/1/0/0 | 0.41, 0.41 | the lowest cut lies elsewhere: foldmap's cut at normalized cut 0.33, CATH's domains at 0.41, 0.41 (the criterion prefers another boundary) |
| 3EOM | C | 381 | 3 | 3 | 0.51 | 0.40 | 0.50 | boundary moved | tightly packed, large | 1/1/0/1 | 0.52, 0.35, 0.36 | the lowest cut lies elsewhere: foldmap's cut at normalized cut 0.28, CATH's domains at 0.52, 0.35, 0.36 (the criterion prefers another boundary) |
| 1GQ2 | B | 555 | 2 | 5 | 0.54 | 0.29 | 1.00 | extra split | - | 1/0/0/5 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.15, 0.11, 0.00) |
| 1GQ2 | D | 555 | 2 | 5 | 0.50 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.15, 0.10, 0.00) |
| 1GQ2 | E | 555 | 2 | 4 | 0.70 | 0.40 | 1.00 | extra split | - | 1/0/0/3 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.15, 0.10) |
| 1GQ2 | F | 555 | 2 | 5 | 0.50 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.16, 0.07, 0.01) |
| 1GQ2 | G | 555 | 2 | 4 | 0.72 | 0.40 | 1.00 | extra split | - | 1/0/0/3 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.16, 0.07) |
| 1GQ2 | H | 555 | 2 | 4 | 0.71 | 0.40 | 1.00 | extra split | - | 1/0/0/3 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.16, 0.10) |
| 1GQ2 | I | 555 | 2 | 5 | 0.49 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.16, 0.09, 0.01) |
| 1GQ2 | J | 555 | 2 | 5 | 0.50 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.16, 0.06, 0.00) |
| 1GQ2 | K | 555 | 2 | 5 | 0.50 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.16, 0.07, 0.01) |
| 1GQ2 | M | 555 | 2 | 5 | 0.50 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.14, 0.09, 0.01) |
| 1GQ2 | N | 555 | 2 | 5 | 0.50 | 0.25 | 1.00 | extra split | - | 1/0/0/6 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.15, 0.10, 0.01) |
| 1GQ2 | O | 555 | 2 | 5 | 0.56 | 0.33 | 1.00 | extra split | - | 1/0/0/4 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.15, 0.08, 0.01) |
| 1GQ2 | P | 555 | 2 | 4 | 0.71 | 0.40 | 1.00 | extra split | - | 1/0/0/3 | 0.26, 0.26 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.42, 0.15, 0.12) |
| 4ZH4 | F | 542 | 3 | 5 | 0.69 | 0.57 | 1.00 | extra split | large | 2/0/0/3 | 0.05, 0.09, 0.00 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.95, 0.88) |
| 4ZH4 | L | 539 | 3 | 5 | 0.70 | 0.57 | 1.00 | extra split | large | 2/0/0/3 | 0.05, 0.09, 0.00 | a cut below DOMAIN_NCUT 0.37 inside a CATH domain (its confidence: 0.98, 0.87) |

## 2. The confidence of each split

A split's confidence is 1 - (its normalized cut / DOMAIN_NCUT): 0 at the limit, where the part would have stayed whole, 1 for a part with no contact across (features.split_confidence, said by domain_confidence in the notes, `foldmap check`, its --json report, the MCP result and the app; never drawn). It is the margin below a tuned limit, not a probability, and moves with that limit. A split is right here when each place its two sides meet lies within 20 residues of a CATH boundary. Per group: right of all in each bin, with the Wilson 95% interval of the fraction right (scripts/paper/common.wilson); the AUC is the chance a right split's confidence is above a wrong one's.

| confidence | dom_cal (tuning: DOMAIN_NCUT was chosen on it) | dom_held + dom_fresh (regression) | r100_third, every chain | r100_third, its first chain per entry |
|---|---|---|---|---|
| 0.00-0.25 | 5/12 [0.19, 0.68] | 3/11 [0.10, 0.57] | 3/49 [0.02, 0.17] | 3/8 [0.14, 0.69] |
| 0.25-0.50 | 9/14 [0.39, 0.84] | 9/19 [0.27, 0.68] | 6/24 [0.12, 0.45] | 3/5 [0.23, 0.88] |
| 0.50-0.75 | 3/5 [0.23, 0.88] | 7/9 [0.45, 0.94] | 11/15 [0.48, 0.89] | 3/5 [0.23, 0.88] |
| 0.75-1.00 | 2/2 [0.34, 1.00] | 5/5 [0.57, 1.00] | 7/11 [0.35, 0.85] | 3/3 [0.44, 1.00] |
| AUC (right vs wrong) | 0.61 | 0.75 | 0.87 | 0.74 |
| splits (right); median confidence right / wrong | 33 (19); 0.33 / 0.29 | 44 (24); 0.49 / 0.27 | 99 (27); 0.63 / 0.16 | 21 (12); 0.55 / 0.23 |

On the third sample, 35 of the 46 wrong splits under 0.25 are 1GQ2's (copies of one protein); counted once per entry (its first chain) the bins hold a few splits each. dom_cal is where DOMAIN_NCUT was chosen, so its column is in-sample; the regression column and the per-entry column are the independent views.

## 3. Keeping a weakly supported region undivided

The reading measured: a region whose cut has confidence t or less is not cut, and never weighed again, so its stronger sub-cuts go with it (the literal "kept undivided"). That is a cut made only below DOMAIN_NCUT x (1 - t), as the confidence is the margin below that limit (up to its rounding to 0.01), measured by setting features.DOMAIN_NCUT so. Merging only across a weak drawn split, keeping the cuts inside it, is another reading, not measured. For each t on every set (chains meeting the targets of all chains; mean overlap, F1, recall):

| t | DOMAIN_NCUT x (1 - t) | dom_cal meet | dom_cal overlap / F1 / recall | dom_held meet | dom_held overlap / F1 / recall | dom_fresh meet | dom_fresh overlap / F1 / recall | r100_third meet | r100_third overlap / F1 / recall | third: entries failing S5 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 0.3700 | 9/22 | 0.874 / 0.821 / 0.889 | 2/11 | 0.806 / 0.731 / 0.871 | 7/11 | 0.862 / 0.779 / 0.849 | 67/109 | 0.859 / 0.735 / 0.893 | 12 |
| 0.05 | 0.3515 | 8/22 | 0.850 / 0.743 / 0.798 | 2/11 | 0.839 / 0.745 / 0.826 | 7/11 | 0.872 / 0.742 / 0.758 | 69/109 | 0.885 / 0.765 / 0.893 | 11 |
| 0.10 | 0.3330 | 8/22 | 0.850 / 0.743 / 0.798 | 3/11 | 0.859 / 0.763 / 0.826 | 7/11 | 0.872 / 0.742 / 0.758 | 69/109 | 0.887 / 0.761 / 0.883 | 11 |
| 0.15 | 0.3145 | 8/22 | 0.830 / 0.673 / 0.707 | 3/11 | 0.859 / 0.763 / 0.826 | 8/11 | 0.892 / 0.773 / 0.758 | 72/109 | 0.896 / 0.779 / 0.874 | 12 |
| 0.20 | 0.2960 | 8/22 | 0.830 / 0.673 / 0.707 | 3/11 | 0.859 / 0.763 / 0.826 | 7/11 | 0.874 / 0.742 / 0.712 | 80/109 | 0.914 / 0.816 / 0.874 | 12 |
| 0.25 | 0.2775 | 8/22 | 0.807 / 0.577 / 0.605 | 2/11 | 0.831 / 0.712 / 0.741 | 7/11 | 0.874 / 0.742 / 0.712 | 82/109 | 0.911 / 0.808 / 0.849 | 11 |
| 0.30 | 0.2590 | 7/22 | 0.763 / 0.472 / 0.507 | 4/11 | 0.766 / 0.530 / 0.541 | 7/11 | 0.874 / 0.742 / 0.712 | 83/109 | 0.915 / 0.827 / 0.858 | 10 |
| 0.35 | 0.2405 | 7/22 | 0.751 / 0.403 / 0.407 | 4/11 | 0.766 / 0.530 / 0.541 | 7/11 | 0.874 / 0.742 / 0.712 | 83/109 | 0.915 / 0.827 / 0.858 | 10 |
| 0.40 | 0.2220 | 6/22 | 0.742 / 0.357 / 0.361 | 3/11 | 0.714 / 0.479 / 0.465 | 7/11 | 0.842 / 0.742 / 0.712 | 83/109 | 0.915 / 0.827 / 0.858 | 10 |
| 0.45 | 0.2035 | 4/22 | 0.690 / 0.251 / 0.261 | 3/11 | 0.713 / 0.485 / 0.465 | 7/11 | 0.821 / 0.697 / 0.682 | 69/109 | 0.845 / 0.682 / 0.700 | 12 |
| 0.50 | 0.1850 | 4/22 | 0.695 / 0.251 / 0.261 | 2/11 | 0.669 / 0.394 / 0.374 | 3/11 | 0.694 / 0.424 / 0.409 | 70/109 | 0.837 / 0.667 / 0.682 | 12 |

Chosen on dom_cal (the largest t at which no dom_cal chain's overlap, F1 or recall drops): none; every t above 0 lowers some dom_cal chain.
Read most favourably (the t with most dom_cal chains meeting the targets): t = 0, as drawn; no t above 0 has as many.
What each set loses is said below at the smallest step, t = 0.05, and at any t a reading chose.

At t = 0.05: some set is worse. What fell, per set ("worse" is any of these):

- dom_cal: 9 -> 8 meet; mean overlap 0.874 -> 0.850; mean F1 0.821 -> 0.743; mean recall 0.889 -> 0.798; lowered: 5DWZ C, 5JDN A, 6FYF A
- dom_held: mean recall 0.871 -> 0.826; lowered: 4XEI A, 6FY9 A
- dom_fresh: mean F1 0.779 -> 0.742; mean recall 0.849 -> 0.758; lowered: 4APE A
- r100_third: nothing fell

Chains whose scores change:

- dom_cal 5DWZ C: domains 2 -> 1, overlap 0.99 -> 0.56, F1 1.00 -> 0.00, recall 1.00 -> 0.00, meets True -> False
- dom_cal 5JDN A: domains 2 -> 1, overlap 0.90 -> 0.75, F1 0.67 -> 0.00, recall 0.67 -> 0.00, meets False -> False
- dom_cal 6FYF A: domains 4 -> 2, overlap 0.75 -> 0.80, F1 0.86 -> 0.80, recall 1.00 -> 0.67, meets False -> False
- dom_held 2C4M A: domains 5 -> 3, overlap 0.62 -> 0.80, F1 0.57 -> 0.80, recall 1.00 -> 1.00, meets False -> False
- dom_held 4XEI A: domains 4 -> 3, overlap 0.92 -> 0.91, F1 0.83 -> 0.80, recall 0.83 -> 0.67, meets False -> False
- dom_held 6FY9 A: domains 5 -> 2, overlap 0.60 -> 0.79, F1 0.83 -> 0.80, recall 1.00 -> 0.67, meets False -> False
- dom_fresh 4APE A: domains 3 -> 2, overlap 0.76 -> 0.87, F1 0.40 -> 0.00, recall 1.00 -> 0.00, meets False -> False
- r100_third 3EOM C: domains 3 -> 2, overlap 0.51 -> 0.69, F1 0.40 -> 0.50, recall 0.50 -> 0.50, meets False -> False
- r100_third 1GQ2 B: domains 5 -> 4, overlap 0.54 -> 0.69, F1 0.29 -> 0.33, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 D: domains 5 -> 4, overlap 0.50 -> 0.72, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 F: domains 5 -> 4, overlap 0.50 -> 0.72, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 I: domains 5 -> 4, overlap 0.49 -> 0.71, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 J: domains 5 -> 4, overlap 0.50 -> 0.72, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 K: domains 5 -> 4, overlap 0.50 -> 0.72, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 M: domains 5 -> 4, overlap 0.50 -> 0.72, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 N: domains 5 -> 4, overlap 0.50 -> 0.72, F1 0.25 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1GQ2 O: domains 5 -> 4, overlap 0.56 -> 0.71, F1 0.33 -> 0.40, recall 1.00 -> 1.00, meets False -> False
- r100_third 1XZ4 B: domains 2 -> 1, overlap 0.59 -> 1.00, F1 0.00 -> 1.00, recall 1.00 -> 1.00, meets False -> True
- r100_third 1XZ4 D: domains 2 -> 1, overlap 0.59 -> 1.00, F1 0.00 -> 1.00, recall 1.00 -> 1.00, meets False -> True

The third sample alone would gain: more of its chains meet the targets at t = 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, most at t = 0.30 (67 -> 83), where 12 of its 17 new passes are 1GQ2's copies. It cannot choose t: it is no tuning set (sets.json), and one entry's copies decide its count.

Ruling Q42 commits the change only when the t chosen on dom_cal leaves no set worse: it is not committed; the confidence ships alone.

## Limits

- The third sample has been run and looked at since plan 15; its failures are described here, never tuned on. dom_cal is a tuning set, dom_held and dom_fresh regression sets (sets.json).
- CATH is one chopping; another (SCOP, ECOD) may draw some boundaries elsewhere. A wrong split here is one CATH does not draw.
- The classes are read from domain counts and boundaries; a chain's mechanism is the first that its cut record explains, not every one that may contribute.
- The confidence measures how far a cut is from not being made; it says nothing of where else a boundary could lie. It moves with DOMAIN_NCUT, a sensitive constant (docs/validation/sensitivity: x1.25 moves the third sample's S5 failures from 12 to 19).
