# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""What each numbers.json key is evidence of (plan 28 Task 2, ruling Q19), for make_numbers.py.

Every key gets a validation layer (LAYERS: structural consistency, external validation, readability, regression guard
or fidelity; null for what is no measurement or mixes layers, with a layer_note saying which) and the benchmark set it
was measured on, with that set's tag from sets.json (tuning, regression or unseen), by the pattern table FAMILIES. A
key no family names, or a set sets.json lacks, stops the run (label).
"""

from __future__ import annotations

import json
import re

from common import HARNESS, ZERO_BY_CONSTRUCTION
from m2 import THICKNESS_FIXED

LAYERS = {
    "consistency": "structural consistency: the figure against the structure's own coordinates, DSSP and pairing "
    "(internal: S1-S4, S6, the homologue checks H1-H5)",
    "external": "external validation: against a reference foldmap does not compute (CATH, RNApolis, mkdssp, OPM, "
    "US-align, published figures)",
    "readability": "readability: the page criteria A1-A20 and the gates built on them",
    "regression-guard": "a regression guard: zero or identical by construction where the method works as designed "
    "(the NA angle and joint errors on an ordinary projection, cross-platform determinism); no evidence of fidelity, "
    "so it carries no interval",
    "fidelity": "fidelity: how far the page keeps the 3D arrangement (Spearman, Procrustes, mirrored stems, contact "
    "proximity)",
}
GATE = (
    "the entry gate (no criterion A1-A20 failed) is readability (science.KIND); it includes A10 (secondary structure "
    "against mkdssp) and A17 (the NA metrics, pairing against RNApolis), external comparisons quoted on their own "
    "(dssp_*, bp_*, stem_agreement)"
)
# (key pattern, layer, set, note): the first pattern matching the whole key labels it; "{1}" in a set is the match's
# first group. A layer of None: no measurement (seeds, entry counts, tables, timings) or one mixing layers (its note
# says which); a set of None: no benchmark set.
FAMILIES = [
    (r"r100_(tuned|held|third|fourth)_(seed|cutoff)", None, "r100_{1}", None),
    (r"r100_held_pass_first", "readability", "r100_held_first", GATE),
    (r"r100_(tuned|held|third|fourth)_pass", "readability", "r100_{1}", GATE),
    (r"r100_third_(fail_by_check|fail_by_stratum|not_judged|a12_pass|a13_pass)", "readability", "r100_third", None),
    (
        r"r100_third_science(_fail)?",
        None,
        "r100_third",
        "mixes layers: S1-S4 and S6 are structural consistency checks, S5 external validation (against CATH), quoted "
        "on its own as dom_r100_*",
    ),
    (r"dom_r100_(n|pass|per_sample)", "external", "r100_pooled", "S5: automatic domains against CATH's chopping"),
    (r"limit_dssp_agree", "external", None, "A10's limit: secondary structure against mkdssp"),
    (r"dssp_r100_(agree|residues|entries|tuned)", "external", "r100_tuned", "A10: against mkdssp"),
    (r"dssp_r100_(held|third)", "external", "r100_{1}", "A10: against mkdssp"),
    (r"supp_table_r100|runtime_all|runtime_fit|supp_fig_runtime", None, "r100_pooled", None),
    (r"runtime_(median_s|p95_s)", None, "r100_third", None),
    (r"na_n_entries|supp_table_na", None, "na_benchmark", None),
    (r"bp_precision|bp_recall|stem_agreement|na_macro_means", "external", "na_benchmark", "against RNApolis"),
    (
        r"na_campaign_pass",
        "readability",
        "na_benchmark",
        "A17's gate: every nacheck.THRESHOLDS limit, which mixes external (pairing against RNApolis), fidelity, "
        "regression guards, consistency (beads) and readability (nacheck.ROLE; extra.na_fidelity)",
    ),
    (r"na_regression_guards", "regression-guard", "na_benchmark", ZERO_BY_CONSTRUCTION),
    (r"na_(spearman|procrustes|mirrored|contact_auc|stems_nearer)_\w+|na_fidelity", "fidelity", "na_benchmark", None),
    (r"pk_\w+", "external", "na_benchmark", "against RNApolis"),
    (r"dom_cal_\w+", "external", "dom_cal", "against CATH"),
    (r"dom_held_nohom_\w+", "external", "dom_held_nohom", "against CATH"),
    (r"dom_held_\w+|dom_group_held", "external", "dom_held", "against CATH"),
    (r"dom_fresh_\w+", "external", "dom_fresh", "against CATH"),
    (r"dom_all33", "external", "dom_reference", "against CATH"),
    (r"supp_table_domains", None, "dom_reference", None),
    (r"hom_usalign_h1_under", "external", "hom_heldout", "foldmap's aligner against US-align"),
    (r"hom_usalign\w*", "external", "hom_usalign", "foldmap's aligner against US-align"),
    (r"hom_(seed|n_pairs|n_pos|n_neg)", None, "hom_heldout", None),
    (r"hom_\w+", "consistency", "hom_heldout", "H1-H5: B drawn like A, judged against A and B themselves (internal)"),
    (r"dssp_\w+", "external", "dssp_m2", "against mkdssp 3.1.2"),
    (r"supp_table_dssp", None, "dssp_m2", None),
    (
        r"membrane_thickness_mad",
        None,
        "membrane_opm",
        f"OPM's thickness against a constant: the thickness is {THICKNESS_FIXED}, so this measures OPM's spread "
        "around it, not an external validation of an estimate",
    ),
    (r"membrane_\w+", "external", "membrane_opm", "against OPM"),
    (r"determinism_(n|platforms)", None, "determinism", "counts of the determinism run"),
    (r"determinism_\w+", "regression-guard", "determinism", "identical by design on every platform: no interval"),
    (r"pub_\w+", "external", "pub", "against published hand-drawn figures"),
    (r"release_tag|zenodo_doi", None, None, None),
]
RATE = re.compile(r"(.+)_(n|ci_lo|ci_hi)")  # a rate's denominator and interval keys (make_numbers.Numbers.rate)


def load_sets() -> dict:
    """scripts/paper/sets.json: each benchmark set's tag (tuning, regression or unseen) and why."""
    doc = json.loads((HARNESS / "sets.json").read_text(encoding="utf-8"))
    bad = {k: s["tag"] for k, s in doc["sets"].items() if s["tag"].split()[0] not in doc["_legend"]}
    if bad:
        raise SystemExit(f"sets.json: tags must start with tuning, regression or unseen: {bad}")
    return doc


def family(key: str, known=frozenset()) -> tuple[str | None, str | None, str | None]:
    """(layer, set, note) of a key by FAMILIES; a rate's _n, _ci_lo and _ci_hi keys (`known` holds their base) take
    their base key's. KeyError for a key no family names."""
    m = RATE.fullmatch(key)
    if m and m.group(1) in known:
        key = m.group(1)
    for pattern, layer, where, note in FAMILIES:
        got = re.fullmatch(pattern, key)
        if got:
            return layer, where.replace("{1}", got.group(1)) if where and "{1}" in where else where, note
    raise KeyError(key)


def label(nb, sets: dict) -> None:
    """Every key and extra of `nb` (make_numbers.Numbers: .keys, .extra) gets its layer, set and the set's tag
    (FAMILIES, sets.json); a key no family names, or a set sets.json lacks, stops the run."""
    known = set(nb.keys) | set(nb.extra)
    unknown, missing = [], set()
    for where in (nb.keys, nb.extra):
        for key, rec in where.items():
            try:
                layer, s, note = family(key, known)
            except KeyError:
                unknown.append(key)
                continue
            if s is not None and s not in sets["sets"]:
                missing.add(s)
                continue
            rec["layer"] = layer
            rec["set"] = s
            rec["set_tag"] = sets["sets"][s]["tag"] if s else None
            if note:
                rec["layer_note"] = note
    if unknown or missing:
        raise SystemExit(f"keys no family names: {sorted(unknown)}; sets not in sets.json: {sorted(missing)}")
