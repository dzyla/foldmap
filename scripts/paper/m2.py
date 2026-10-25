# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The M2 keys (plan 14 M2), computed from the raw outputs saved under ../data/m2/ (and ../data/opm/), for
make_numbers.py. Nothing here is typed in: every value is counted from those files, and each key's `how` says how.

  DSSP vs mkdssp   data/m2/dssp/          (dssp_mkdssp.py on zyla-lab)
  membrane vs OPM  data/m2/membrane/      (membrane_opm.py on zyla-lab; OPM files from opm_fetch.py in data/opm/)
  pseudoknots      data/m2/na/report.json (scripts/na_validate.py; its code from its own verified stamp, else
                                          labelled "no source stamp", as the current, pre-stamp file is)
  determinism      data/m2/determinism/   (scripts/paper/determinism.py: the Linux leg on zyla-lab; the CI artifacts)
"""

from __future__ import annotations

import gzip
import json
from collections import Counter, defaultdict

from common import MANUSCRIPT, TABLES, median, na_rows, quotable, stamped_commit, tex_escape, tex_frac, tex_int

DATA = MANUSCRIPT / "data"
M2 = DATA / "m2"
Q3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}  # everything else (T, S read as '-', '-') is coil


def pct(x: float) -> str:
    """A fraction as a percentage, rounded to the fewest decimal places (one at least) at which it does not round
    to 100 (99.99950 -> 99.9995, not 100.0 and not a truncated 99.9)."""
    if x >= 1:
        return r"100\%"
    for places in range(1, 9):
        v = round(x * 100, places)
        if v < 100:
            return f"{v:.{places}f}" + r"\%"
    return f"{x * 100:.8f}" + r"\%"


def q3(c: str) -> str:
    return Q3.get(c, "C")


# ---- DSSP vs mkdssp --------------------------------------------------------------------------------------------------


def build_dssp(nb, reports: dict) -> None:
    d = M2 / "dssp"
    run_path = d / "run.json"
    if not run_path.exists():
        return
    run = json.loads(run_path.read_text(encoding="utf-8"))
    if run["failures"]:
        raise SystemExit(f"dssp run has failures: {run['failures']}")
    entries = {}
    for f in sorted(d.glob("*.json.gz")):
        doc = json.loads(gzip.decompress(f.read_bytes()))
        entries[doc["id"]] = doc
    src = [
        "data/m2/dssp/*.json.gz, raw/*.dssp.gz, run.json (scripts/paper/dssp_mkdssp.py on zyla-lab: mkdssp 3.1.2, "
        "foldmap from build/tree)"
    ]

    def symmetric(rows):
        """mkdssp's bridges read both ways: a residue whose BP column names a partner is that partner's partner too
        (mkdssp keeps two partner slots per residue, so a third bridge shows on the partner's line only)."""
        sym = defaultdict(set)
        for c, n, ic, _mk, _fm, mkp, _fmp in rows:
            k = f"{c}:{n}:{ic}"
            sym[k].update(mkp)
            for q in mkp:
                sym[q].add(k)
        return sym

    def tally(ids):
        t = Counter()
        chains = set()
        diffs = []
        for i in ids:
            sym = symmetric(entries[i]["rows"])
            for c, n, ic, mk, fm, mkp, fmp in entries[i]["rows"]:
                ms = sym.get(f"{c}:{n}:{ic}", set())
                if ms or fmp:
                    t["bridged_symmetric"] += 1
                    t["bridge_same_symmetric"] += ms == set(fmp)
                t["residues"] += 1
                chains.add((i, c))
                t["q8"] += mk == fm
                t["q3"] += q3(mk) == q3(fm)
                if mkp or fmp:
                    t["bridged"] += 1
                    t["bridge_same"] += mkp == fmp
                if mk != fm or mkp != fmp:
                    diffs.append(
                        {
                            "id": i,
                            "chain": c,
                            "number": n,
                            "icode": ic,
                            "mkdssp": mk,
                            "foldmap": fm,
                            "mkdssp_partners": mkp,
                            "foldmap_partners": fmp,
                            "one_sided": mk == fm and ms == set(fmp),
                        }
                    )
        t["chains"] = len(chains)
        return t, diffs

    # the A10 counts of each sample's report, reproduced
    for s, rep in reports.items():
        ids = run["samples"][s]
        t, _ = tally(ids)
        measured = [e for e in rep["entries"] if (e.get("dssp") or {}).get("agree") is not None]
        res = sum(e["dssp"]["residues"] for e in measured)
        agree = sum(e["dssp"]["residues"] - len(e["dssp"]["differ"]) for e in measured)
        if (t["residues"], t["q8"], len(ids)) != (res, agree, len(measured)):
            raise SystemExit(
                f"DSSP {s}: this run {t['residues']}/{t['q8']}/{len(ids)} vs A10 {res}/{agree}/{len(measured)}"
            )
    ids = sorted(entries)
    t, diffs = tally(ids)
    n = len(ids)
    desc = f"the {n} protein entries of the three random samples (all whose file mkdssp reads)"
    nb.put(
        "dssp_set_description",
        {"entries": n, "ids": ids},
        desc,
        src,
        "every entry of the three random samples whose report measured A10 (dssp.agree not None); the per-sample "
        "residue and agreement counts equal the A10 counts of each report (asserted)",
    )
    nb.put(
        "dssp_n_chains",
        t["chains"],
        tex_int(t["chains"]),
        src,
        "distinct (entry, author chain) among compared residues",
    )
    nb.put(
        "dssp_n_residues",
        t["residues"],
        tex_int(t["residues"]),
        src,
        "residues compared: read by both programs, model as deposited, chromophore stand-ins and four either side "
        "left out (criteria.dssp_agreement)",
    )
    nb.put(
        "dssp_q8_mean",
        t["q8"] / t["residues"],
        pct(t["q8"] / t["residues"]),
        src,
        f"pooled over residues: {t['q8']} of {t['residues']} with the same DSSP code (mkdssp S and blank read as "
        "'-'); text rounded to the fewest places that do not round to 100",
    )
    nb.put(
        "dssp_q3_mean",
        t["q3"] / t["residues"],
        pct(t["q3"] / t["residues"]),
        src,
        f"pooled: {t['q3']} of {t['residues']} in the same class (helix H, G, I; strand E, B; other the rest)",
    )
    nb.put(
        "dssp_bridge_agreement",
        t["bridge_same"] / t["bridged"],
        pct(t["bridge_same"] / t["bridged"]),
        src,
        f"pooled: {t['bridge_same']} of {t['bridged']} residues with a bridge partner in either program have the "
        "same set of partners in both (mkdssp BP1/BP2; foldmap DsspResult.bridges)",
    )
    nb.put(
        "dssp_n_diff_residues",
        len(diffs),
        str(len(diffs)),
        src,
        "residues whose code or bridge-partner set differs (Table S6 lists each)",
    )
    # per-chain means and per-sample tallies, for the record
    per_chain = defaultdict(lambda: [0, 0, 0])
    for i in ids:
        for c, _n, _ic, mk, fm, _a, _b in entries[i]["rows"]:
            pc = per_chain[(i, c)]
            pc[0] += 1
            pc[1] += mk == fm
            pc[2] += q3(mk) == q3(fm)
    nb.put(
        "dssp_chain_means",
        {
            "chains": len(per_chain),
            "q8": sum(v[1] / v[0] for v in per_chain.values()) / len(per_chain),
            "q3": sum(v[2] / v[0] for v in per_chain.values()) / len(per_chain),
            "chains_all_agree_q8": sum(v[1] == v[0] for v in per_chain.values()),
        },
        "",
        src,
        "means over chains (the keys are pooled over residues)",
        extra=True,
    )
    nb.put(
        "dssp_bridge_agreement_symmetric",
        t["bridge_same_symmetric"] / t["bridged_symmetric"],
        pct(t["bridge_same_symmetric"] / t["bridged_symmetric"]),
        src,
        f"chosen after the first count was seen: mkdssp's partners read both ways (symmetric): "
        f"{t['bridge_same_symmetric']} of {t['bridged_symmetric']} residues",
        extra=True,
    )
    nb.put("dssp_totals", dict(t) | {"entries": n, "mkdssp": run["mkdssp"]}, "", src, "all tallies", extra=True)
    causes = {}
    cpath = M2 / "dssp" / "causes.json"
    if cpath.exists():
        causes = json.loads(cpath.read_text(encoding="utf-8"))
    one_sided = (
        "same code; mkdssp lists the extra bridge on the partner residue's line only (two BP columns per "
        "residue): read both ways, the partners agree"
    )
    for x in diffs:
        x["cause"] = (
            one_sided
            if x["one_sided"]
            else causes.get(f"{x['id']} {x['chain']}:{x['number']}:{x['icode']}", "not diagnosed")
        )
    nb.put(
        "dssp_differing",
        diffs,
        "",
        src,
        "every differing residue, with the cause from data/m2/dssp/causes.json "
        "(diagnosed by reading both outputs; text, not a number)",
        extra=True,
    )
    write_dssp_table(diffs)
    nb.put(
        "supp_table_dssp",
        "tables/dssp.tex",
        r"\input{tables/dssp.tex}",
        src,
        "every residue whose code or bridge partners differ, with its cause",
    )


def write_dssp_table(diffs: list[dict]) -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    lines = [
        r"\scriptsize",
        r"\begin{tabular}{@{}llllp{7.5cm}@{}}",
        r"\toprule",
        r"Entry & Residue & mkdssp & \foldmap & Cause \\",
        r"\midrule",
    ]
    for x in diffs:
        what_m = x["mkdssp"] + (
            "" if x["mkdssp_partners"] == x["foldmap_partners"] else " (" + ", ".join(x["mkdssp_partners"]) + ")"
        )
        what_f = x["foldmap"] + (
            "" if x["mkdssp_partners"] == x["foldmap_partners"] else " (" + ", ".join(x["foldmap_partners"]) + ")"
        )
        lines.append(
            f"{x['id']} & {x['chain']} {x['number']}{x['icode']} & {tex_escape(what_m)} & {tex_escape(what_f)} & "
            f"{tex_escape(x['cause'])} \\\\"
        )
    if not diffs:
        lines.append(r"\multicolumn{5}{c}{none} \\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        "",
        r"\noindent Codes as DSSP's, mkdssp's bend and blank read as `-'; bridge partners (chain:number:"
        r"insertion code) in brackets where they differ.",
    ]
    (TABLES / "dssp.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---- membrane vs OPM -------------------------------------------------------------------------------------------------


def build_membrane(nb) -> None:
    path = M2 / "membrane" / "results.json"
    if not path.exists():
        return
    doc = json.loads(path.read_text(encoding="utf-8"))
    rows = doc["entries"]
    bad = [e["id"] for e in rows if e.get("error")]
    if bad:
        raise SystemExit(f"membrane: entries with errors {bad}")
    src = [
        "data/m2/membrane/results.json (scripts/paper/membrane_opm.py on zyla-lab, foldmap from build/tree)",
        "data/opm/*.json, *.pdb.gz, selection.json (scripts/paper/opm_fetch.py: OPM's API and "
        "coordinate files, opm-back.cc.lehigh.edu and opm-assets.storage.googleapis.com)",
    ]
    sel = doc["selection"]
    n = len(rows)
    found = [e for e in rows if e["estimate"] is not None]
    refused = [e["id"] for e in rows if e["estimate"] is None]
    nb.put(
        "membrane_n_entries",
        n,
        str(n),
        src,
        f"rule: {sel['rule']}; {len(sel['candidates'])} candidates, {n} taken, left out (not in OPM): "
        + ", ".join(x["id"] for x in sel["left_out"]),
    )
    td = [abs(e["estimate"]["thickness_diff"]) for e in found]
    nb.put(
        "membrane_thickness_mad",
        sum(td) / len(td),
        f"{sum(td) / len(td):.1f}",
        src,
        f"mean |2 x half - OPM's DUM-plane distance| over the {len(found)} entries foldmap finds a membrane in "
        f"(refused: {', '.join(refused) or 'none'})",
    )
    ang = [e["estimate"]["normal_angle_deg"] for e in found]
    nb.put(
        "membrane_tilt_median",
        median(ang),
        f"{median(ang):.0f}",
        src,
        f"median angle between foldmap's normal and OPM's (as axes, 0-90 deg) over the {len(found)} found; OPM's "
        "normal taken into the drawn frame by a Kabsch fit of its CA atoms",
    )
    segs = [(e["id"], s) for e in rows for s in e["segments"] if "crosses_opm" in s]
    hit = [x for x in segs if x[1].get("crosses_estimate")]
    nb.put(
        "membrane_tm_recall",
        len(hit) / len(segs),
        tex_frac(len(hit) / len(segs), 2),
        src,
        f"pooled over OPM's transmembrane segments with >= 3 CA atoms in OPM's file, every entry: {len(hit)} of "
        f"{len(segs)} cross the estimated slab (CA run reaching membrane._MIN_CROSSING = 10 A beyond its mid-plane "
        "on both sides); a refused entry's segments count as not crossing",
    )
    called = [e for e in found if e["estimate"]["sided"]]
    right = [e for e in called if e["estimate"]["side_correct"]]
    nb.put(
        "membrane_side_called",
        len(called),
        str(len(called)),
        src,
        "found entries whose sides foldmap calls (positive-inside rule, or a barrel's loops)",
    )
    nb.put(
        "membrane_side_correct",
        len(right),
        str(len(right)),
        src,
        "of those, outward normal on OPM's outer side (dot with OPM's outward normal > 0; OPM 'O' DUM atoms outer, "
        "'N' inner)",
    )
    # extras: refusals, found-only recall, OPM's own ceiling, offsets, by class, attribution
    segs_found = [x for x in segs if any(e["id"] == x[0] for e in found)]
    ceiling = [x for x in segs if x[1]["crosses_opm"]]
    off = [e["estimate"]["midplane_offset"] for e in found]
    by_class = {}
    for cls in sorted({e["classtype"] for e in rows}):
        es = [e for e in rows if e["classtype"] == cls]
        fs = [e for e in es if e["estimate"] is not None]
        ss = [(e["id"], s) for e in es for s in e["segments"] if "crosses_opm" in s]
        by_class[cls] = {
            "entries": [e["id"] for e in es],
            "found": len(fs),
            "thickness_mad": (sum(abs(e["estimate"]["thickness_diff"]) for e in fs) / len(fs)) if fs else None,
            "angle_median": median([e["estimate"]["normal_angle_deg"] for e in fs]) if fs else None,
            "recall": (sum(bool(s.get("crosses_estimate")) for _, s in ss) / len(ss)) if ss else None,
            "segments": len(ss),
        }
    topo = [e.get("opm_topology_check", {}).get("agrees") for e in rows]
    nb.put(
        "membrane_detail",
        {
            "found": len(found),
            "refused": refused,
            "recall_found_only": (sum(bool(s.get("crosses_estimate")) for _, s in segs_found) / len(segs_found))
            if segs_found
            else None,
            "segments_found_only": len(segs_found),
            "opm_own_ceiling": len(ceiling) / len(segs),
            "segments": len(segs),
            "recall_of_opm_crossing": (sum(bool(s.get("crosses_estimate")) for _, s in ceiling) / len(ceiling)),
            "segments_opm_crossing": len(ceiling),
            "span_overlap_post_hoc": overlap(rows, found),
            "segments_too_few_ca": sum(1 for e in rows for s in e["segments"] if "crosses_opm" not in s),
            "midplane_offset_median": median(off) if off else None,
            "angle_max": max(ang) if ang else None,
            "angles": {e["id"]: e["estimate"]["normal_angle_deg"] for e in found},
            "fit_rmsd_max": max(e["fit"]["rmsd"] for e in rows),
            "opm_topology_check_agrees": Counter(str(x) for x in topo),
            "by_class": by_class,
        },
        "",
        src,
        "membrane vs OPM: refusals, found-only recall, OPM's own planes as the ceiling, mid-plane offsets, "
        "per OPM class, frame-fit RMSD, OPM's DUM side convention checked against its topology field",
        extra=True,
    )
    attr = {e["id"]: e["attribution"] for e in rows}
    for x in doc.get("attribution_extra", []):
        attr[x["id"]] = x.get("attribution") or x.get("error")
    nb.put(
        "membrane_rule_attribution",
        attr,
        "",
        src,
        "the estimate with each a963bb6 rule on or off (fallback: best spanned slab; chain rule: oily/polar over "
        "crossing chains only): None = refused, else the slab centre",
        extra=True,
    )
    write_membrane_table(rows)


def overlap(rows: list[dict], found: list[dict]) -> float | None:
    """Chosen AFTER the results were seen (the pre-stated crossing rule punishes a slab a few A off): the mean over
    OPM segments of the fraction of each segment's depth span (min to max CA depth under the estimate) that lies
    inside the estimated slab. Reported as an extra, never as the key."""
    half = {e["id"]: e["estimate"]["half"] for e in found}
    fr = []
    for e in rows:
        for s in e["segments"]:
            if "depth_estimate" not in s:
                continue
            lo, hi = s["depth_estimate"]
            h = half[e["id"]]
            fr.append(max(0.0, min(hi, h) - max(lo, -h)) / (hi - lo) if hi > lo else 0.0)
    return sum(fr) / len(fr) if fr else None


def write_membrane_table(rows: list[dict]) -> None:
    lines = [
        r"\scriptsize",
        r"\begin{tabular}{@{}lllrrrrll@{}}",
        r"\toprule",
        r"Entry & OPM class & Found & Thickness diff. (\AA) & Normal angle ($^\circ$) & Mid-plane offset (\AA) & "
        r"TM recall & Sides & Fit RMSD (\AA) \\",
        r"\midrule",
    ]
    for e in rows:
        est = e["estimate"]
        segs = [s for s in e["segments"] if "crosses_opm" in s]
        rec = f"{sum(bool(s.get('crosses_estimate')) for s in segs)}/{len(segs)}"
        cls = "barrel" if "barrel" in e["classtype"].lower() else "helical"
        if est is None:
            lines.append(f"{e['id']} & {cls} & refused & & & & {rec} & & {e['fit']['rmsd']:.2f} \\\\")
            continue
        side = "not called" if not est["sided"] else ("right" if est["side_correct"] else "wrong")
        lines.append(
            f"{e['id']} & {cls} & yes & {est['thickness_diff']:+.1f} & {est['normal_angle_deg']:.1f} & "
            f"{est['midplane_offset']:.1f} & {rec} & {side} & {e['fit']['rmsd']:.2f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    (TABLES / "membrane.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---- pseudoknots -----------------------------------------------------------------------------------------------------


def _na_code(raw) -> str:
    """The NA rerun's code, as a source label: its own stamp, checked by verify_stamp (ruling Q11; Q18). A stamp that
    does not verify (or a file from before the stamps) is labelled "no source stamp ...", never quoted at a commit."""
    if not (isinstance(raw, dict) and raw.get("source_stamp")):
        return stamped_commit(raw)
    from verify_stamp import verify_report

    checked = verify_report(raw)
    if checked and all(r["verified"] for r in checked):
        return f"{stamped_commit(raw)}, verified"
    return "no source stamp verified (" + "; ".join(r["why"] for r in checked if not r["verified"]) + ")"


def build_pk(nb, committed_na: list[dict]) -> None:
    path = M2 / "na" / "report.json"
    if not path.exists():
        return
    raw = json.loads(path.read_text(encoding="utf-8"))
    quotable("data/m2/na/report.json", raw)
    rep = na_rows(raw)
    passed = sum(not e["failed"] for e in rep)
    old = {e["entry"]: e for e in committed_na}
    skip = ("seconds", "source_stamp")
    same = all(
        {k: v for k, v in e.items() if k not in skip and not k.startswith(("pk_", "deposited_pk_"))}
        == {k: v for k, v in old[e["entry"]].items() if k not in skip}
        for e in rep
    )
    if not same or passed != sum(not e["failed"] for e in committed_na):
        raise SystemExit("NA rerun differs from the committed report beyond pk_* and seconds")
    commit = _na_code(raw)  # from the report's own stamp, verified (rulings Q11, Q18; no COMMIT file)
    src = [
        f"data/m2/na/report.json (scripts/na_validate.py --require-rnapolis, {commit}; every other "
        "field equal to docs/validation/na/report.json, asserted)"
    ]
    ref = sum(e["pk_reference"] for e in rep)
    got = sum(e["pk_found"] for e in rep)
    nb.put(
        "pk_reference",
        ref,
        str(ref),
        src,
        "pooled over entries: reference stems (RNApolis pairs grouped as foldmap groups its own) that foldmap's "
        "knot rule (nucleic._knots) calls pseudoknotted, on the pairs the figure draws",
    )
    nb.put("pk_found", got, str(got), src, "of those, recalled by a drawn stem (the two-thirds match of stem_recall)")
    nb.put(
        "pk_deposited",
        {
            "found": sum(e["deposited_pk_found"] for e in rep),
            "reference": sum(e["deposited_pk_reference"] for e in rep),
        },
        "",
        src,
        "the same on the model as deposited (deposited_pk_*: every reference pair, figure crops aside)",
        extra=True,
    )
    nb.put(
        "pk_per_entry",
        {
            e["entry"]: [e["pk_found"], e["pk_reference"], e["deposited_pk_found"], e["deposited_pk_reference"]]
            for e in rep
            if e["pk_reference"] or e["deposited_pk_reference"]
        },
        "",
        src,
        "[found, reference] as drawn, then on the model as deposited",
        extra=True,
    )


# ---- determinism -----------------------------------------------------------------------------------------------------


def build_determinism(nb) -> None:
    d = M2 / "determinism"
    # the Linux leg (DIR/linux-*.json) and CI's artifacts, as `gh run download -p 'determinism-*'` lays them out
    # (DIR/determinism-<os>/<os>.json, beside <os>-svg/)
    runs = [json.loads(f.read_text(encoding="utf-8")) for f in sorted(d.rglob("*.json"))] if d.exists() else []
    if not runs:
        return
    src = ["data/m2/determinism/**/*.json (scripts/paper/determinism.py: Liberation Sans forced, fonts.py)"]
    figs = sorted({f["structure"] for r in runs for f in r["figures"]})

    def same(name):
        got = {
            (f.get("signature_sha256"), f.get("svg_sha256"))
            for r in runs
            for f in r["figures"]
            if f["structure"] == name
        }
        return len(got) == 1 and None not in next(iter(got))

    identical = [f for f in figs if same(f)]

    def one_hash(name, key):
        got = {f.get(key) for r in runs for f in r["figures"] if f["structure"] == name}
        return len(got) == 1 and None not in got

    svg_same = [f for f in figs if one_hash(f, "svg_sha256")]
    sig_same = [f for f in figs if one_hash(f, "signature_sha256")]
    systems = sorted({r["system"] for r in runs})
    fonts_ok = all(r["font"]["font_family"] == "Liberation Sans" for r in runs)
    detail = {
        "runs": [
            {k: r[k] for k in ("label", "system", "machine", "python", "numpy", "matplotlib", "env")}
            | {"font_file": r["font"]["font_file"]}
            for r in runs
        ],
        "figures": len(figs),
        "identical": len(identical),
        "differ": [f for f in figs if f not in identical],
        "svg_identical": len(svg_same),
        "signature_identical": len(sig_same),
        "signature_differs_in": {
            f: sorted(
                r["label"]
                for r in runs
                for x in r["figures"]
                if x["structure"] == f
                and x["signature_sha256"]
                != runs[0]["figures"][[y["structure"] for y in runs[0]["figures"]].index(f)]["signature_sha256"]
            )
            for f in figs
            if f not in sig_same
        },
        "systems": systems,
        "all_liberation_sans": fonts_ok,
    }
    nb.put("determinism_runs", detail, "", src, "every run of determinism.py found", extra=True)
    if len(systems) < 2:
        return  # one platform: the keys stay pending (needs the CI run)
    nb.put(
        "determinism_n",
        len(figs),
        str(len(figs)),
        src,
        "figures drawn on every platform (default view of each committed test structure)",
    )
    nb.put(
        "determinism_platforms",
        systems,
        f"{len(systems)} platforms ({', '.join(systems)})",
        src,
        "distinct platform.system() among the runs",
    )
    nb.put(
        "determinism_identical",
        len(identical),
        str(len(identical)),
        src,
        "figures whose layout signature and SVG hash are the same in every run",
    )
