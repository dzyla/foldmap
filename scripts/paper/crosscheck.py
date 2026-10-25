# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Step 6 (local, light): numbers written in the prose of main.tex / supplement.tex (not \\NUM) against the sources.

    python scripts/paper/crosscheck.py

Each check finds its claim in the text by a pattern (a claim no longer found is reported, so a check never passes by
silence) and compares it with the value read from the sources at the harness commit. Writes build/crosscheck.json and
exits 1 when any check disagrees or a claim is not found. The README's "Disagreements" section lists what it found.
"""

from __future__ import annotations

import ast
import json
import re
import sys

from common import BUILD, MANUSCRIPT, TREE, load_source

TEXT = {n: (MANUSCRIPT / f"{n}.tex").read_text(encoding="utf-8") for n in ("main", "supplement")}
FLAT = {n: " ".join(t.split()) for n, t in TEXT.items()}


def constant(relpath: str, name: str):
    tree = ast.parse((TREE / relpath).read_text(encoding="utf-8"))
    for node in tree.body:
        targets = (
            node.targets if isinstance(node, ast.Assign) else [node.target] if isinstance(node, ast.AnnAssign) else []
        )
        if any(isinstance(t, ast.Name) and t.id == name for t in targets):
            return ast.literal_eval(node.value)
    raise KeyError(name)


results = []


def check(name: str, where: str, pattern: str, source_value, text_value=None, source: str = "", note: str = ""):
    """pattern: regex in the flattened text; its first group (if any) is the text's value unless text_value given."""
    m = re.search(pattern, FLAT[where])
    if not m:
        results.append(
            {
                "check": name,
                "where": where,
                "status": "claim not found",
                "pattern": pattern,
                "source": source,
                "source_value": source_value,
            }
        )
        return
    tv = text_value if text_value is not None else (m.group(1) if m.groups() else m.group(0))
    ok = str(tv) == str(source_value)
    results.append(
        {
            "check": name,
            "where": f"{where}.tex",
            "text": m.group(0)[:200],
            "text_value": tv,
            "source_value": source_value,
            "source": source,
            "status": "agrees" if ok else "DISAGREES",
            "note": note,
        }
    )


def main() -> int:
    third = load_source("r100_third_report")
    samples = {s: load_source(f"r100_{s}_sample") for s in ("tuned", "held", "third")}
    reports = {s: load_source(f"r100_{s}_report") for s in ("tuned", "held", "third")}
    ref = json.loads((TREE / "tests/data/domains/cath_reference.json").read_text(encoding="utf-8"))["chains"]
    fresh = json.loads((TREE / "tests/data/domains/cath_fresh.json").read_text(encoding="utf-8"))["chains"]
    targets = constant("src/foldmap/domain_score.py", "TARGETS")
    tol = constant("src/foldmap/domain_score.py", "TOLERANCE")
    seed = constant("src/foldmap/domain_score.py", "HELD_OUT_SEED")
    limits = constant("src/foldmap/criteria.py", "LIMITS")
    thr = constant("src/foldmap/nacheck.py", "THRESHOLDS")
    views = constant("scripts/random100.py", "VIEWS")
    quotas = constant("scripts/random100.py", "QUOTAS")
    ribo = constant("scripts/random100.py", "RIBOSOME_NT")
    science = constant("src/foldmap/science.py", "SCIENCE")

    # domains
    check(
        "reference chains",
        "supplement",
        r"holds (\d+) multi-domain chains with CATH choppings",
        len(ref),
        source="tests/data/domains/cath_reference.json",
    )
    check(
        "tuning chains",
        "main",
        r"On the (\d+) tuning chains",
        len(ref) - len(ref) // 3,
        source="cath_reference.json: 33 minus a third",
    )
    check("held-out chains", "main", r"The other (\d+) chains were held out", len(ref) // 3, source="as above")
    check(
        "fresh chains",
        "main",
        r"a further (\d+) chains from CATH superfamilies",
        len(fresh),
        source="tests/data/domains/cath_fresh.json",
    )
    check("boundary tolerance", "main", r"boundaries within (\d+) residues", tol, source="domain_score.TOLERANCE")
    check("held-out seed", "supplement", r"seed \\texttt\{([^}]*)\}", seed, source="domain_score.HELD_OUT_SEED")
    check(
        "domain targets",
        "supplement",
        r"targets \(overlap ([0-9.]+), boundary F1 ([0-9.]+), boundary recall ([0-9.]+)\)",
        f"{targets['overlap']:.2f}/{targets['f1']:.2f}/{targets['recall']:.2f}",
        text_value="/".join(
            re.search(
                r"targets \(overlap ([0-9.]+), boundary F1 ([0-9.]+), boundary recall ([0-9.]+)\)", FLAT["supplement"]
            ).groups()
        )
        if re.search(r"targets \(overlap", FLAT["supplement"])
        else None,
        source="domain_score.TARGETS",
    )
    s5_sources = sorted(
        {c.get("source") for r in reports.values() for e in r["entries"] for c in e.get("domains_cath") or []}
    )
    check(
        "S5 scope in the random samples",
        "main",
        r"Check S5 applies these scores in the campaign only to entries of the 33-chain reference",
        "S5 scores every chain CATH chops (sources: " + ", ".join(s5_sources) + ")",
        text_value="only entries of the 33-chain reference",
        source="report.json domains_cath[].source; random100.domain_reference (Task 8: CATH fetched per entry)",
        note="S5 now fetches CATH's chopping for every PDB entry; the sentence (and NOTES disagreement 1) predate it",
    )
    held_nohom_text = re.search(
        r"seven held-out chains with no training homologue give overlap (\.\d+), recall (\.\d+), "
        r"precision (\.\d+), F1 (\.\d+)",
        " ".join((TREE / "tests/data/domains/SOURCE.txt").read_text(encoding="utf-8").split()),
    )
    groups_tuning = None
    if held_nohom_text:
        nums = json.loads((BUILD / "numbers.json").read_text(encoding="utf-8"))["keys"]
        got = (
            f"{nums['dom_held_nohom_n']['value']} chains: overlap {nums['dom_held_nohom_overlap']['value']:.2f}, "
            f"recall {nums['dom_held_nohom_recall']['value']:.2f}, F1 {nums['dom_held_nohom_f1']['value']:.2f}"
        )
        results.append(
            {
                "check": "held-out chains with no homologue (committed numbers)",
                "where": "tests/data/domains/SOURCE.txt, tests/test_domain_score.py docstring",
                "text_value": "7 chains: overlap "
                + held_nohom_text.group(1)
                + ", recall "
                + held_nohom_text.group(2)
                + ", F1 "
                + held_nohom_text.group(4),
                "source_value": got,
                "source": "domains.py: held chains whose CATH homology group has no tuning chain",
                "status": "DISAGREES",
                "note": "the committed seven drop only the named near-duplicates (6FY9, 2YBE, 1CZA, 1HKB); 4XEI, "
                "2C4M and 1RM5 share a CATH superfamily (their 'group') with tuning chains (1BDG/2YHX, "
                "1GPB, 1LDM/2HKI). The text's wording 'no homologue among the tuning chains' fits the "
                "four; the committed numbers reproduce exactly for the seven (extra.dom_held_nodup).",
            }
        )
    del groups_tuning

    # random samples
    check(
        "strata",
        "main",
        r"in (twelve) strata",
        "twelve"
        if len(third["entries"]) and len(samples["third"]["strata"]) == 12
        else str(len(samples["third"]["strata"])),
        source="sample.json strata",
    )
    pairs = {s: len({e.get("pair") for e in samples[s]["entries"] if e.get("pair")}) for s in samples}
    check(
        "two-state pairs",
        "main",
        r"(five) pairs of one protein in two states",
        "five" if set(pairs.values()) == {5} else str(pairs),
        source="sample.json entries[].pair",
    )
    for s in samples:
        n = len(reports[s]["entries"])
        if n != 100:
            results.append({"check": f"{s} sample size", "status": "DISAGREES", "text_value": 100, "source_value": n})
    rows = re.findall(r"(S\d+) & [^&]+ & (\d+) & ([^&]+) &", FLAT["supplement"])
    strata = {st["name"]: st for st in samples["third"]["strata"]}
    for name, _n, _bins in rows:
        st = strata.get(name)
        if not st:
            continue
        src_bins = ", ".join(f"{b} {v['n']}" for b, v in st["bins"].items()) if st.get("bins") else "--"
        check(
            f"Table S5 {name} n", "supplement", rf"{name} & [^&]+ & (\d+) &", st["n"], source="sample.json strata[].n"
        )
        check(
            f"Table S5 {name} bins",
            "supplement",
            rf"{name} & [^&]+ & \d+ & ([^&]+?) &",
            src_bins,
            source="sample.json strata[].bins",
        )
    check(
        "ligand/glycan/modified quotas",
        "supplement",
        r"at least (\d+) entries with a ligand, (\d+) with a glycan and (\d+) with a modified residue",
        f"{quotas['ligand']}/{quotas['glycan']}/{quotas['modified']}",
        text_value="/".join(
            re.search(
                r"at least (\d+) entries with a ligand, (\d+) with a glycan and (\d+) with a modified",
                FLAT["supplement"],
            ).groups()
        ),
        source="random100.QUOTAS",
    )
    bins = samples["third"]["bins"]
    check("size bins", "supplement", r"S below (\d+)", bins["S"][1] + 1, source="sample.json bins")
    check("XL upper bound", "supplement", r"XL 8,000--(30,000)", f"{bins['XL'][1]:,}", source="sample.json bins")
    check(
        "ribosome-scale RNA",
        "supplement",
        r"ribosome-scale RNA \(\$\\ge (2\{,\}500)\$ nt\)",
        "2{,}500" if ribo == 2500 else ribo,
        source="random100.RIBOSOME_NT",
    )
    check("view timeout", "supplement", r"killed after (\d+)\\,s", limits["view_timeout_s"], source="criteria.LIMITS")
    check(
        "mkdssp version", "main", r"mkdssp (3\.1\.2)", third["versions"]["mkdssp"].split()[-1], source="report versions"
    )
    check(
        "views per entry",
        "main",
        r"Each entry is drawn in every view a user may ask for",
        len(views),
        text_value=11,
        source="random100.VIEWS",
        note="the sentence lists default, two naming schemes, focus, domains, membrane, seed, layout file, "
        "interactive, sequence, app = 11",
    )
    check(
        "science checks",
        "main",
        r"(six) science checks",
        "six" if len(science) == 6 else len(science),
        source="science.SCIENCE",
    )

    # NA limits (Table S4)
    na = {
        "pair F1": (r"pair F1 & [^&]+ & \$\\ge (0\.\d+)\$", f"{thr['pair_f1']:.2f}"),
        "stem recall": (r"stem recall & [^&]+ & \$\\ge (0\.\d+)\$", f"{thr['stem_recall']:.2f}"),
        "Spearman": (r"Spearman & [^&]+ & \$\\ge (0\.\d+)\$", f"{thr['spearman']:.2f}"),
        "Procrustes": (r"Procrustes & [^&]+ & \$\\le (0\.\d+)\$", f"{thr['procrustes']:.2f}"),
        "angle": (
            r"angle median, angle max, joint error & [^&]+ & \$\\le (0\.\d+)\^\\circ\$",
            f"{thr['angle_median']:.2f}"
            if thr["angle_median"] == thr["angle_max"] == thr["joint_error_max"]
            else "differ",
        ),
        "contact AUC": (r"contact AUC & [^&]+ & \$\\ge (0\.\d+)\$", f"{thr['contact_auc']:.2f}"),
        "stems nearer": (r"stems nearer & [^&]+ & (\d\.\d+)", f"{thr['stems_nearer_frac']:.1f}"),
        "seconds": (r"seconds & [^&]+ & \$\\le (\d+)\$", int(thr["seconds"])),
    }
    for name, (pat, val) in na.items():
        check(f"Table S4 {name}", "supplement", pat, val, source="nacheck.THRESHOLDS")

    # claims about when things ran
    check(
        "third sample: no change after its run",
        "main",
        r"was run once on the release described here, with no fix after it",
        "e5b4bae (route._MAX_NODES_WINDOWS) changed src/ after the run at 7696947; it changes 9FZL's domains view "
        "(third sample, still fails)",
        text_value="no fix after it",
        source="build/provenance.json",
        note="true of fixes aimed at the third sample; not true of the code: the release must be 7696947 or the "
        "sample rerun at the frozen tag",
    )
    check(
        "membrane CHECK (stale)",
        "main",
        r"the measurement waits for plan 15 Task 6",
        "plan 15 Task 6 landed at a963bb6 (NOTES.md)",
        text_value="waits for Task 6",
        source="git log a963bb6",
        note="CHECK text only; the measurement itself is M2",
    )
    check(
        "font CHECK (stale)",
        "supplement",
        r"the CI's macOS and Windows jobs do not install the font today",
        "CI installs Liberation Sans 2.1.5 on all three since plan 15 Task 8; render._FONTS still puts Arial first",
        text_value="macOS/Windows do not install it",
        source="NOTES.md (plan 15 Task 8 ruling)",
        note="revise S4.6's CHECK: the harness must force the font (fonts.py)",
    )

    bad = [r for r in results if r["status"] != "agrees"]
    (BUILD / "crosscheck.json").write_text(json.dumps({"checks": results}, indent=1) + "\n", encoding="utf-8")
    print(f"{len(results)} checks, {len(results) - len(bad)} agree, {len(bad)} disagree or not found")
    for r in bad:
        print(f"  {r['status']:16s} {r['check']}: text {r.get('text_value')!r} | source {r.get('source_value')!r}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
