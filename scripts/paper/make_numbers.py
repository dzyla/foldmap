# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Step 3 (light; local or zyla-lab): every \\NUM key the existing campaign outputs can supply, into build/numbers.json,
with the supplementary tables under build/tables/.

    python scripts/paper/make_numbers.py

Reads only what extract.py wrote (build/sources/, build/tree/, build/provenance.json) and domains.py's
build/domains.json. Every value names its source files, the commit each was read at and the commit the measurement
in it was made at. Keys that need runs not made yet are null with "pending": "M2" / "M3" / "M8" / "freeze".
"""

from __future__ import annotations

import ast
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import date

from common import (
    BUILD,
    TABLES,
    TREE,
    load_source,
    median,
    na_rows,
    percentile,
    quotable,
    stamp,
    tex_escape,
    tex_frac,
    tex_int,
)

SAMPLES = {"tuned": "r100_tuned", "held": "r100_held", "third": "r100_third"}
OUTCOMES = ("pass", "fail", "not judged", "does not apply")

PENDING = {
    "release_tag": ("freeze", "the frozen git tag of the release; no tag points at the harness commit"),
    "hom_usalign_agreement": (
        "freeze",
        "foldmap's aligner against US-align (plan 13a-3, measured in plan 23 Task 3): "
        "heldout-summary.json's usalign entry is read once the harness commit includes it",
    ),
    "zenodo_doi": ("M8", "Zenodo release of the frozen tag (author)"),
    "dssp_set_description": ("M2", "DSSP vs mkdssp 3.1.2 on a stated set"),
    "dssp_n_chains": ("M2", "DSSP vs mkdssp run"),
    "dssp_n_residues": ("M2", "DSSP vs mkdssp run"),
    "dssp_q3_mean": ("M2", "DSSP vs mkdssp run"),
    "dssp_q8_mean": ("M2", "DSSP vs mkdssp run"),
    "dssp_bridge_agreement": ("M2", "DSSP vs mkdssp run"),
    "dssp_n_diff_residues": ("M2", "DSSP vs mkdssp run"),
    "supp_table_dssp": ("M2", "DSSP vs mkdssp run: differing residues with causes"),
    "pk_found": ("M2", "na_validate.py's report.json has no pseudoknot count: the NA rerun must add one"),
    "pk_reference": ("M2", "as pk_found"),
    "membrane_n_entries": ("M2", "membrane vs OPM, at or after a963bb6"),
    "membrane_thickness_mad": ("M2", "membrane vs OPM"),
    "membrane_tilt_median": ("M2", "membrane vs OPM"),
    "membrane_tm_recall": ("M2", "membrane vs OPM"),
    "membrane_side_correct": ("M2", "membrane vs OPM"),
    "membrane_side_called": ("M2", "membrane vs OPM"),
    "determinism_n": (
        "CI",
        "needs CI run: scripts/paper/determinism.py on Linux, macOS and Windows with Liberation "
        "Sans forced (.github/workflows/determinism.yml); the Linux leg is in extra.determinism_runs",
    ),
    "determinism_platforms": ("CI", "as determinism_n"),
    "determinism_identical": ("CI", "as determinism_n"),
    "pub_n_structures": ("M3", "redo of work/published-compare at the frozen version"),
    "pub_agree": ("M3", "as pub_n_structures"),
    "pub_total": ("M3", "as pub_n_structures"),
}


# ---- helpers --------------------------------------------------------------------------------------------------------


def constant(relpath: str, name: str):
    """A module-level literal of the code at the harness commit, read with ast (never imported)."""
    tree = ast.parse((TREE / relpath).read_text(encoding="utf-8"))
    for node in tree.body:
        targets = (
            node.targets if isinstance(node, ast.Assign) else [node.target] if isinstance(node, ast.AnnAssign) else []
        )
        if any(isinstance(t, ast.Name) and t.id == name for t in targets):
            return ast.literal_eval(node.value)
    raise KeyError(f"{name} not in {relpath}")


def outcome(c: dict) -> str:
    if c.get("pass") is True:
        return "pass"
    if c.get("pass") is False:
        return "fail"
    return c.get("status") or "not judged"


def top(counts: dict[str, int], ranks: int) -> list[tuple[str, int]]:
    """The items holding the `ranks` highest distinct counts (ties kept), highest first, then by name order."""
    nonzero = {k: v for k, v in counts.items() if v}
    levels = sorted(set(nonzero.values()), reverse=True)[:ranks]
    return [(k, v) for k, v in sorted(nonzero.items(), key=lambda kv: (-kv[1], order_key(kv[0]))) if v in levels]


def order_key(code: str):
    """A1 < A1c < A2 < A2b < ... < A20; S1 < ... < S12."""
    import re

    m = re.match(r"([A-Z]+)(\d+)(\w*)", code)
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (code, 0, "")


def join_and(parts: list[str]) -> str:
    return parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]


# ---- the sources ----------------------------------------------------------------------------------------------------


class Numbers:
    def __init__(self):
        self.prov = stamp()
        self.keys: dict[str, dict] = {}
        self.extra: dict[str, dict] = {}

    def src(self, *names: str) -> list[dict]:
        out = []
        for n in names:
            rec = self.prov["sources"].get(n)
            if rec is None:  # code or build file, not a campaign output
                out.append({"file": n, "read_at": self.prov["harness_commit_short"]})
                continue
            m = rec.get("measured_at", {})
            out.append(
                {
                    "file": rec["path"],
                    "read_at": rec["read_at"],
                    "measured_at": (m.get("report") or m.get("sample") or "")[:7] or None,
                    "measured_at_from": rec.get("measured_at_from"),
                }
            )
        return out

    def put(self, key: str, value, text: str, sources: list[str], how: str, extra: bool = False):
        (self.extra if extra else self.keys)[key] = {
            "value": value,
            "text": text,
            "sources": self.src(*sources),
            "how": how,
        }


def sample_summary(report: dict) -> dict:
    entries = report["entries"]
    crit = list(report["criteria"])
    by_crit = {c: Counter(outcome(e["criteria"][c]) for e in entries) for c in crit}
    strata = defaultdict(lambda: [0, 0])
    for e in entries:
        strata[e["stratum"]][0] += 1
        strata[e["stratum"]][1] += not e["pass"]
    science_keys = sorted({k for e in entries for k in e["science"]}, key=order_key)
    by_sci = {k: Counter(outcome(e["science"][k]) for e in entries) for k in science_keys}
    sci_fail_entries = [e["id"] for e in entries if any(outcome(s) == "fail" for s in e["science"].values())]
    recount = sum(1 for e in entries if e["pass"])
    if recount != report["passed"]:
        raise SystemExit(f"report says {report['passed']} passed, its entries {recount}")
    return {
        "passed": report["passed"],
        "n": len(entries),
        "criteria": crit,
        "by_criterion": by_crit,
        "strata": dict(strata),
        "science": by_sci,
        "science_fail_entries": sci_fail_entries,
    }


def build_r100(nb: Numbers) -> dict:
    reports = {s: load_source(f"{p}_report") for s, p in SAMPLES.items()}
    samples = {s: load_source(f"{p}_sample") for s, p in SAMPLES.items()}
    for s, p in SAMPLES.items():
        quotable(f"{p}_report", reports[s])  # never a report of mixed code or one not for the paper
    summ = {s: sample_summary(r) for s, r in reports.items()}

    for s, p in SAMPLES.items():
        nb.put(f"r100_{s}_seed", samples[s]["seed"], str(samples[s]["seed"]), [f"{p}_sample"], "sample.json seed")
        nb.put(
            f"r100_{s}_cutoff",
            samples[s]["query_date"],
            samples[s]["query_date"],
            [f"{p}_sample"],
            "sample.json query_date (release-date cutoff)",
        )
    for s, p in SAMPLES.items():
        nb.put(
            f"r100_{s}_pass",
            summ[s]["passed"],
            str(summ[s]["passed"]),
            [f"{p}_report"],
            "entries with no criterion failed (report.json passed, recounted from the entries)",
        )

    first = load_source("r100_held_first_report")
    quotable("r100_held_first_report", first)
    fs = sample_summary(first)
    nb.put(
        "r100_held_pass_first",
        fs["passed"],
        str(fs["passed"]),
        ["r100_held_first_report"],
        "the held-out sample's first full run as judged after F1 fix round 1 (report.json.gz at 1ef0706; controller "
        "ruling 2026-10-17). Its entries were measured at 5650dbd (92) and 44a1c1b (8 redrawn): see _meta",
    )

    t = summ["third"]
    fails = {c: t["by_criterion"][c]["fail"] for c in t["criteria"]}
    picked = top(fails, 3)
    nb.put(
        "r100_third_fail_by_check",
        dict(sorted(fails.items(), key=lambda kv: order_key(kv[0]))),
        join_and([f"{c} ({n} entries)" if i == 0 else f"{c} ({n})" for i, (c, n) in enumerate(picked)]),
        ["r100_third_report"],
        "entries failing each criterion; text: the three highest counts, ties kept",
    )
    by_stratum = {k: v[1] for k, v in t["strata"].items()}
    picked = top(by_stratum, 2)
    nb.put(
        "r100_third_fail_by_stratum",
        {k: {"entries": v[0], "failed": v[1]} for k, v in t["strata"].items()},
        join_and(
            [
                f"{k} ({n} of {t['strata'][k][0]} entries)" if i == 0 else f"{k} ({n} of {t['strata'][k][0]})"
                for i, (k, n) in enumerate(picked)
            ]
        ),
        ["r100_third_report"],
        "failed entries per stratum; text: the two highest counts, ties kept",
    )
    nj = {c: t["by_criterion"][c]["not judged"] for c in t["criteria"]}
    picked = top(nj, 2)
    nb.put(
        "r100_third_not_judged",
        dict(sorted(nj.items(), key=lambda kv: order_key(kv[0]))),
        join_and([f"{c} ({n} entries)" if i == 0 else f"{c} ({n})" for i, (c, n) in enumerate(picked)]),
        ["r100_third_report"],
        "entries not judged per criterion; text: the two highest counts, ties kept",
    )
    for c in ("A12", "A13"):
        n = t["by_criterion"][c]["pass"]
        nb.put(f"r100_third_{c.lower()}_pass", n, str(n), ["r100_third_report"], f"entries passing {c}")
    nb.put(
        "r100_third_science_fail",
        len(t["science_fail_entries"]),
        str(len(t["science_fail_entries"])),
        ["r100_third_report"],
        "entries failing at least one science check S1-S6 (entry science outcomes)",
    )
    failing = {k: v["fail"] for k, v in t["science"].items() if v["fail"]}
    clean = [k for k, v in t["science"].items() if not v["fail"]]
    parts = [f"{k} in {n} entries" for k, n in sorted(failing.items(), key=lambda kv: order_key(kv[0]))]
    text = join_and(parts) if parts else "none"
    if clean:
        text += f"; {join_and(clean)} in none"
    nb.put(
        "r100_third_science",
        {k: dict(v) for k, v in t["science"].items()},
        text,
        ["r100_third_report"],
        "per science check, entries by outcome; text: checks with failures and the checks with none",
    )

    # S5 over the three samples: every chain CATH chops, scored at the run (domains_cath), within the targets?
    targets = constant("src/foldmap/domain_score.py", "TARGETS")
    chains = [(s, e["id"], c) for s, r in reports.items() for e in r["entries"] for c in e.get("domains_cath") or []]
    within = [x for x in chains if all(x[2][k] >= targets[k] for k in ("overlap", "f1", "recall"))]
    nb.put(
        "dom_r100_n",
        len(chains),
        str(len(chains)),
        [f"{p}_report" for p in SAMPLES.values()],
        "chains scored by S5 (domains_cath) over the three samples: every chain CATH chops in a PDB entry",
    )
    nb.put(
        "dom_r100_pass",
        len(within),
        str(len(within)),
        [f"{p}_report" for p in SAMPLES.values()] + ["src/foldmap/domain_score.py"],
        f"of those, chains with overlap, boundary F1 and recall all at or above domain_score.TARGETS {targets}",
    )
    per = {}
    for s in SAMPLES:
        cs = [x for x in chains if x[0] == s]
        ok = [x for x in within if x[0] == s]
        multi = [x for x in cs if x[2]["n_ref"] >= 2]
        per[s] = {
            "chains": len(cs),
            "within_targets": len(ok),
            "multi_domain_chains": len(multi),
            "entries": len({x[1] for x in cs}),
            "from_reference_set": sum(x[2].get("source") == "reference" for x in cs),
            "S5_entries": dict(summ[s]["science"].get("S5", {})),
        }
    nb.put("dom_r100_per_sample", per, "", [f"{p}_report" for p in SAMPLES.values()], "S5 per sample", extra=True)

    # A10 over the tuned sample (README quotes it)
    lim = constant("src/foldmap/criteria.py", "LIMITS")
    if lim["dssp_agree"] != reports["third"]["limits"]["dssp_agree"]:
        raise SystemExit("criteria.LIMITS dssp_agree differs from the limit the runs used")
    nb.put(
        "limit_dssp_agree",
        lim["dssp_agree"],
        f"{lim['dssp_agree']:.2f}",
        ["src/foldmap/criteria.py", "r100_third_report"],
        "criteria.LIMITS['dssp_agree'] at the harness commit (equal to the runs' limits)",
    )
    for s in ("tuned", "held", "third"):
        measured = [e for e in reports[s]["entries"] if (e.get("dssp") or {}).get("agree") is not None]
        res = sum(e["dssp"]["residues"] for e in measured)
        agree = sum(e["dssp"]["residues"] - len(e["dssp"]["differ"]) for e in measured)
        check = sum(round(e["dssp"]["agree"] * e["dssp"]["residues"]) for e in measured)
        if check != agree:
            raise SystemExit(f"{s}: dssp differ lists ({agree}) disagree with agree x residues ({check})")
        if s == "tuned":
            src = ["r100_tuned_report"]
            nb.put("dssp_r100_agree", agree, tex_int(agree), src, "sum over measured entries of residues - len(differ)")
            nb.put("dssp_r100_residues", res, tex_int(res), src, "sum of dssp.residues over measured entries")
            nb.put("dssp_r100_entries", len(measured), str(len(measured)), src, "entries with dssp.agree measured")
        nb.put(
            f"dssp_r100_{s}",
            {"agree": agree, "residues": res, "entries": len(measured)},
            f"{tex_int(agree)} of {tex_int(res)}",
            [f"{SAMPLES[s]}_report"],
            "A10 residues, per sample",
            extra=True,
        )

    write_r100_tables(summ, reports, fs)
    nb.put(
        "supp_table_r100",
        "tables/r100.tex",
        r"\input{tables/r100.tex}",
        [f"{p}_report" for p in SAMPLES.values()],
        "outcome counts per criterion and science check for the three samples; per entry and criterion in "
        "tables/r100-entries.csv",
    )
    return reports


def write_r100_tables(summ: dict, reports: dict, first: dict):
    TABLES.mkdir(parents=True, exist_ok=True)
    head = " & ".join(rf"\multicolumn{{4}}{{c}}{{{s} ({summ[s]['passed']}/{summ[s]['n']} pass)}}" for s in SAMPLES)
    lines = [
        r"\scriptsize",
        r"\begin{tabular}{@{}l" + "rrrr" * len(SAMPLES) + "@{}}",
        r"\toprule",
        r" & " + head + r" \\",
        r"Check & " + " & ".join(["P", "F", "NJ", "NA"] * len(SAMPLES)) + r" \\",
        r"\midrule",
    ]
    for c in summ["third"]["criteria"]:
        cells = [str(summ[s]["by_criterion"][c][o]) for s in SAMPLES for o in OUTCOMES]
        lines.append(f"{c} & " + " & ".join(cells) + r" \\")
    lines.append(r"\midrule")
    for k in summ["third"]["science"]:
        cells = [str(summ[s]["science"][k][o]) for s in SAMPLES for o in OUTCOMES]
        lines.append(f"{k} & " + " & ".join(cells) + r" \\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        "",
        r"\noindent P pass, F fail, NJ not judged, NA does not apply. Science checks (S1--S6) are not part of "
        r"an entry's pass. Per entry and criterion: \texttt{r100-entries.csv}.",
    ]
    (TABLES / "r100.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    crit = summ["third"]["criteria"]
    sci = list(summ["third"]["science"])
    with open(TABLES / "r100-entries.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample", "id", "stratum", "bin", "pass", *crit, *sci, "failed"])
        for s, r in reports.items():
            for e in r["entries"]:
                w.writerow(
                    [
                        s,
                        e["id"],
                        e["stratum"],
                        e["bin"],
                        e["pass"],
                        *(outcome(e["criteria"][c]) for c in crit),
                        *(outcome(e["science"][k]) for k in sci),
                        "; ".join(e["failed"])[:300],
                    ]
                )


def build_na(nb: Numbers):
    raw = load_source("na_report")
    quotable("na_report", raw)
    rep = na_rows(raw)
    src = ["na_report"]
    n = len(rep)
    nb.put("na_n_entries", n, str(n), src, "entries in report.json")
    tp = sum(e["pairs_foldmap"] - len(e["extra_pairs"]) for e in rep)
    tp2 = sum(e["pairs_reference"] - len(e["missed_pairs"]) for e in rep)
    if tp != tp2:
        raise SystemExit(f"NA: true pairs {tp} from extra_pairs, {tp2} from missed_pairs")
    drawn = sum(e["pairs_foldmap"] for e in rep)
    ref = sum(e["pairs_reference"] for e in rep)
    nb.put(
        "bp_precision",
        tp / drawn,
        tex_frac(tp / drawn),
        src,
        f"pooled over entries, pairs as drawn: {tp} true of {drawn} drawn (canonical and wobble)",
    )
    nb.put("bp_recall", tp / ref, tex_frac(tp / ref), src, f"pooled: {tp} found of {ref} RNApolis pairs")
    macro = {k: sum(e[k] for e in rep) / n for k in ("pair_precision", "pair_recall", "pair_f1", "stem_recall")}
    nb.put("na_macro_means", macro, "", src, "means over entries (the pooled numbers are the keys)", extra=True)
    stems_ref = sum(e["stems_reference"] for e in rep)
    stems_found = sum(round(e["stem_recall"] * e["stems_reference"]) for e in rep)
    nb.put(
        "stem_agreement",
        stems_found / stems_ref,
        tex_frac(stems_found / stems_ref),
        src,
        f"pooled: {stems_found} of {stems_ref} reference stems found (stem_recall x stems_reference per entry)",
    )
    ok = [e["entry"] for e in rep if not e["failed"]]
    nb.put(
        "na_campaign_pass", len(ok), str(len(ok)), src, "entries whose failed list is empty (every THRESHOLDS limit)"
    )
    write_na_table(rep)
    nb.put("supp_table_na", "tables/na.tex", r"\input{tables/na.tex}", src, "per entry: pairs, stems, pass or failures")


def write_na_table(rep: list[dict]):
    lines = [
        r"\scriptsize",
        r"\begin{tabular}{@{}lrrrrrrp{6.2cm}@{}}",
        r"\toprule",
        r"Entry & Pairs drawn & Reference & Precision & Recall & Stems (ref.) & Stem recall & Outcome \\",
        r"\midrule",
    ]
    for e in rep:
        out = "pass" if not e["failed"] else tex_escape("; ".join(e["failed"]))
        lines.append(
            f"{e['entry']} & {e['pairs_foldmap']} & {e['pairs_reference']} & {e['pair_precision']:.3f} & "
            f"{e['pair_recall']:.3f} & {e['stems_reference']} & {e['stem_recall']:.3f} & {out} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    (TABLES / "na.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_domains(nb: Numbers):
    path = BUILD / "domains.json"
    if not path.exists():
        for k in ("cal", "held", "held_nohom", "fresh"):
            for m in ("overlap", "recall", "f1"):
                nb.keys[f"dom_{k}_{m}"] = {"value": None, "pending": "harness", "why": "run domains.py (zyla-lab)"}
        return
    d = json.loads(path.read_text(encoding="utf-8"))
    src = [
        "tests/data/domains/cath_reference.json",
        "tests/data/domains/cath_inputs.npz",
        "src/foldmap/domain_score.py",
        "build/domains.json (domains.py on zyla-lab, code from build/tree)",
    ]
    fresh_src = [
        "tests/data/domains/cath_fresh.json",
        "tests/data/domains/fresh_inputs.npz",
        "src/foldmap/domain_score.py",
        "build/domains.json (domains.py on zyla-lab, code from build/tree)",
    ]

    def mean(rows, k):
        return sum(r[k] for r in rows) / len(rows)

    sets = {
        "cal": (
            [r for r in d["reference"] if r["split"] == "tuning"],
            src,
            "the 22 tuning chains: cath_reference.json minus the chain-by-chain held third (seed 20261016:domains)",
        ),
        "held": (
            [r for r in d["reference"] if r["split"] == "held"],
            src,
            "the 11 chains held out chain by chain (random.Random('20261016:domains').sample of the sorted keys)",
        ),
        "held_nohom": (
            [r for r in d["reference"] if r["held_nohom"]],
            src,
            "held chains whose CATH homology group (cath_reference.json 'group': shared superfamily, "
            "directly or through another chain) has no tuning chain",
        ),
        "fresh": (d["fresh"], fresh_src, "cath_fresh.json, the 11 chains from superfamilies none of the 33 has"),
    }
    for k, (rows, s, how) in sets.items():
        for m in ("overlap", "recall", "f1", "precision"):
            v = mean(rows, m)
            nb.put(
                f"dom_{k}_{m}",
                v,
                tex_frac(v, 2),
                s,
                f"mean {m} over {len(rows)} chains: {how}",
                extra=(m == "precision"),
            )
    nohom = sets["held_nohom"][0]
    nb.put(
        "dom_held_nohom_n",
        len(nohom),
        str(len(nohom)),
        src,
        "held chains with no homologue among the tuning chains: "
        + ", ".join(f"{r['pdb']} {r['chain']}" for r in nohom),
    )
    base = mean(sets["cal"][0], "baseline_overlap")
    nb.put(
        "dom_cal_baseline_overlap",
        base,
        tex_frac(base, 2),
        src,
        "mean overlap over the 22 tuning chains when each chain is called one domain",
    )
    # what the committed docstrings call "no training homologue" (seven chains: held minus the named near-duplicates)
    nodup = [r for r in d["reference"] if r["split"] == "held" and r["pdb"] not in ("6FY9", "2YBE", "1CZA", "1HKB")]
    nb.put(
        "dom_held_nodup",
        {m: mean(nodup, m) for m in ("overlap", "recall", "precision", "f1")} | {"n": len(nodup)},
        "",
        src,
        "held chains minus the near-duplicates named in test_domain_score.py (6FY9, 2YBE, 1CZA, 1HKB): "
        "the committed 'seven with no training homologue'",
        extra=True,
    )
    every = d["reference"]
    nb.put(
        "dom_all33",
        {m: mean(every, m) for m in ("overlap", "recall", "precision", "f1")},
        "",
        src,
        "all 33 reference chains",
        extra=True,
    )
    grp = [r for r in d["reference"] if r["group_held"]]
    nb.put(
        "dom_group_held",
        {m: mean(grp, m) for m in ("overlap", "recall", "precision", "f1")} | {"n": len(grp)},
        "",
        src,
        "domain_score.held_out as now (by homology group, after the fact)",
        extra=True,
    )
    write_domain_table(d)
    nb.put(
        "supp_table_domains",
        "tables/domains.tex",
        r"\input{tables/domains.tex}",
        src + fresh_src[:2],
        "per chain: set, CATH domains, automatic domains, overlap, boundary recall, precision, F1",
    )


def write_domain_table(d: dict):
    lines = [
        r"\scriptsize",
        r"\begin{tabular}{@{}llrrrrrr@{}}",
        r"\toprule",
        r"Chain & Set & CATH & Auto & Overlap & Recall & Precision & F1 \\",
        r"\midrule",
    ]
    rows = sorted(d["reference"], key=lambda r: (r["split"] != "tuning", r["pdb"])) + d["fresh"]
    for r in rows:
        if "split" in r:
            what = "tuning" if r["split"] == "tuning" else ("held-out" if r["held_nohom"] else "held-out$^\\dagger$")
        else:
            what = "fresh"
        lines.append(
            f"{r['pdb']} {r['chain']} & {what} & {r['n_ref']} & {r['n_auto']} & {r['overlap']:.2f} & "
            f"{r['recall']:.2f} & {r['precision']:.2f} & {r['f1']:.2f} \\\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        "",
        r"\noindent $^\dagger$ shares a CATH homology group with a tuning chain.",
    ]
    (TABLES / "domains.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_homologues(nb: Numbers):
    """Plan 13a's held-out homologue pairs (docs/validation/homologues/heldout-summary.json, Task 13a-11)."""
    try:
        h = load_source("hom_summary")
    except (OSError, KeyError, ValueError):
        return  # the harness commit predates the held-out run: the keys stay unset
    src = ["hom_summary"]
    pos, neg = h["positives"], h["negatives"]
    nb.put("hom_seed", h["seed"], str(h["seed"]), src, "homologues.py sample --seed")
    nb.put("hom_n_pairs", pos["pairs"] + neg["pairs"], str(pos["pairs"] + neg["pairs"]), src, "pairs in heldout.json")
    nb.put("hom_n_pos", pos["pairs"], str(pos["pairs"]), src, "positive pairs (remote + close strata)")
    nb.put("hom_n_neg", neg["pairs"], str(neg["pairs"]), src, "negative pairs (one architecture, two topologies)")
    nb.put("hom_pos_pass", pos["all_pass"], str(pos["all_pass"]), src, "positives passing H1-H5")
    nb.put("hom_neg_refused", neg["refused"], str(neg["refused"]), src, "negatives refused with the note (H1 pass)")
    for k in ("H1", "H2", "H3", "H4", "H5"):
        nb.put(f"hom_pos_{k.lower()}", pos["pass"][k], str(pos["pass"][k]), src, f"positives passing {k}")
    for key, field in (
        ("hom_h2_indep", "h2_median_independent"),
        ("hom_h2_alike", "h2_median_alike"),
        ("hom_h3_indep", "h3_share_independent"),
        ("hom_h3_alike", "h3_share_alike"),
    ):
        v = pos[field]
        nb.put(
            key,
            v,
            "--" if v is None else f"{v:.2f}",
            src,
            f"median over the {pos['judged_h2']} positives with H2 judged of the pairs' {field}",
        )
    nb.put("hom_strata", h["strata"], "", src, "per stratum: pairs, results, each check's outcomes", extra=True)
    us = h.get("usalign")  # plan 23 Task 3: the aligner against US-align, added to the summary after the run
    if not us:
        return  # the harness commit predates the comparison: hom_usalign_agreement stays pending
    a, lim = us["all"], us["limits"]
    scope = (
        f"tuning + held-out pairs ({a['pairs']}) at US-align TM-score >= {lim['usalign_tm']}, both normalised by "
        f"the reference A; foldmap's structure route against US-align {us['version']} (measured at "
        f"{us['measured_at']})"
    )
    nb.put(
        "hom_usalign_agreement",
        a["both"],
        f"{a['both']} of {a['judged']}",
        src,
        f"pairs meeting both 13a-3 limits (TM-score within {lim['tm_tolerance']} of US-align's and >= "
        f"{lim['pair_agreement']:.0%} of US-align's residue pairs reproduced within ±{lim['pair_slack']}), of the "
        f"{scope}",
    )
    nb.put("hom_usalign_judged", a["judged"], str(a["judged"]), src, f"the {scope}")
    nb.put(
        "hom_usalign_tm_within",
        a["tm_within"],
        f"{a['tm_within']} of {a['judged']}",
        src,
        f"of the {scope}: foldmap's TM-score within {lim['tm_tolerance']} of US-align's",
    )
    nb.put(
        "hom_usalign_pairs_median",
        a["within_1_median"],
        f"{100 * a['within_1_median']:.0f}\\%",
        src,
        f"of the {scope}: median share of US-align's residue pairs foldmap reproduces within ±{lim['pair_slack']}",
    )
    nb.put(
        "hom_usalign_tm_diff_median",
        a["tm_difference_median"],
        f"{a['tm_difference_median']:.3f}",
        src,
        f"over all {a['pairs']} pairs: median of foldmap's TM-score minus US-align's",
    )
    refused = us["h1_refused"]
    under = sum(not x["usalign_one_fold"] for x in refused)
    nb.put(
        "hom_usalign_h1_under",
        under,
        f"{under} of {len(refused)}",
        src,
        "held-out positives H1 refuses (foldmap TM-score under SAME_FOLD_TM) that US-align also puts under it",
    )
    nb.put("hom_usalign", us, "", src, "the comparison per set and over all (heldout-summary.json usalign)", extra=True)


def build_runtime(nb: Numbers, reports: dict):
    from runtime import rows_from_reports

    rows = rows_from_reports(reports)
    third = [r["process_s"] for r in rows if r["sample"] == "third"]
    srcs = ["r100_third_report"]
    nb.put(
        "runtime_median_s",
        median(third),
        f"{median(third):.1f}",
        srcs,
        "third sample: median over 100 entries of the default view's process wall time (views.default.seconds, "
        "what A2 gates: drawing, measuring and for nucleic acid RNApolis), run with --jobs 10 on zyla-lab",
    )
    nb.put(
        "runtime_p95_s",
        percentile(third, 95),
        f"{percentile(third, 95):.1f}",
        srcs,
        "third sample: 95th percentile (linear interpolation) of the same",
    )
    allp = [r["process_s"] for r in rows]
    alld = [r["drawing_s"] for r in rows]
    thd = [r["drawing_s"] for r in rows if r["sample"] == "third"]
    nb.put(
        "runtime_all",
        {
            "n": len(allp),
            "process_median_s": median(allp),
            "process_p95_s": percentile(allp, 95),
            "drawing_median_s": median(alld),
            "drawing_p95_s": percentile(alld, 95),
            "third_drawing_median_s": median(thd),
            "third_drawing_p95_s": percentile(thd, 95),
            "max_process_s": max(allp),
        },
        "",
        [f"{p}_report" for p in SAMPLES.values()],
        "pooled over the three samples (300 entries)",
        extra=True,
    )
    fig = BUILD / "runtime.pdf"
    if fig.exists():
        nb.put(
            "supp_fig_runtime",
            "runtime.pdf",
            r"\includegraphics[height=2.8cm]{runtime.pdf}",
            [f"{p}_report" for p in SAMPLES.values()],
            "runtime.py (zyla-lab): default view time vs residues",
        )
    else:
        nb.keys["supp_fig_runtime"] = {"value": None, "pending": "harness", "why": "run runtime.py (zyla-lab)"}
    rt = BUILD / "runtime.json"
    if rt.exists():
        nb.put(
            "runtime_fit",
            json.loads(rt.read_text(encoding="utf-8")),
            "",
            [f"{p}_report" for p in SAMPLES.values()],
            "runtime.py: power-law fit and per-bin medians",
            extra=True,
        )


def main() -> int:
    import m2

    nb = Numbers()
    reports = build_r100(nb)
    build_na(nb)
    build_domains(nb)
    build_homologues(nb)
    build_runtime(nb, reports)
    m2.build_dssp(nb, reports)
    m2.build_membrane(nb)
    m2.build_pk(nb, na_rows(load_source("na_report")))
    m2.build_determinism(nb)
    tags = nb.prov.get("tags_at_commit") or []
    for key, (stage, why) in PENDING.items():
        if key == "release_tag" and tags:
            nb.put(key, tags[0], tags[0], [], "git tag pointing at the harness commit")
            continue
        nb.keys.setdefault(key, {"value": None, "pending": stage, "why": why})

    prov = nb.prov
    notes = [
        "The three random-100 reports record foldmap_commit 7696947 for the report and every entry. Per "
        "docs/validation/random100/review.md (838c579), the run itself recorded 'unknown' (zyla-lab's synced tree has "
        "no .git); 7696947 was stamped into the reports afterwards from that synced tree and the reports were judged "
        "again from the measurements, with the same counts.",
        "One code commit lies between 7696947 and the harness commit 838c579: e5b4bae (route._MAX_NODES_WINDOWS, a "
        "coarse page's fine windows share a pool of their own). Per its message and the reviews' after-the-run notes "
        "it changes 7NZ2's domains view (held-out sample; 7NZ2 still fails) and 9FZL's (third sample; still fails, 5 "
        "-> 2 fallback curves). The counts here are as measured at 7696947; they were not rerun at 838c579.",
        "The NA report.json carries no commit; its measurement commit (7696947) is read from review.md's last rerun "
        "heading. e5b4bae (routing of huge pages' windows) is not expected to touch the 39 NA entries, but they were "
        "not rerun at 838c579.",
        "r100_held_pass_first is read at 1ef0706: the held-out first run (drawings at 5650dbd for 92 entries, 44a1c1b "
        "for the 8 redrawn after the check corrections of F1 fix round 1). It is a historical count by design.",
        "Domain scores (dom_cal/held/fresh) are computed by domains.py from the committed fixtures with the code at "
        "the harness commit; the fixtures were built at bbe1bd8 and CI checks four of them against their structures.",
    ]
    doc = {
        "_meta": {
            "harness_commit": prov["harness_commit"],
            "harness_commit_short": prov["harness_commit_short"],
            "harness_commit_subject": prov["harness_commit_subject"],
            "tags_at_commit": tags,
            "generated": date.today().isoformat(),
            "generator": "scripts/paper/make_numbers.py",
            "sources": prov["sources"],
            "source_problems": prov["problems"],
            "reviews_after_the_run": prov["reviews_after_the_run"],
            "notes": notes,
            "pending_legend": {
                "M2": "validation still missing (plan 14 M2)",
                "CI": "needs a CI run on other platforms (job in scripts/paper/README.md)",
                "M3": "comparisons (plan 14 M3)",
                "M8": "author hand-off",
                "freeze": "the frozen release tag",
                "harness": "a harness step on zyla-lab not run yet",
            },
        },
        "keys": dict(sorted(nb.keys.items())),
        "extra": dict(sorted(nb.extra.items())),
    }
    (BUILD / "numbers.json").write_text(json.dumps(doc, indent=1, default=dict) + "\n", encoding="utf-8")
    filled = [k for k, v in nb.keys.items() if v.get("value") is not None]
    pending = [k for k, v in nb.keys.items() if v.get("value") is None]
    print(f"numbers.json: {len(filled)} keys filled, {len(pending)} pending, {len(nb.extra)} extra")
    for k in sorted(filled):
        print(f"  {k:28s} {nb.keys[k]['text'][:90]}")
    for k in sorted(pending):
        print(f"  {k:28s} pending: {nb.keys[k]['pending']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
