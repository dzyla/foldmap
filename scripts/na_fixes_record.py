# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The record of plan 28's nucleic-acid fixes (ruling Q37): each fix measured before and after on the seven failing
cases of the third random sample and the 39 nucleic-acid controls that passed (na-diagnosis.md, "Fixes").

    python scripts/na_fixes_record.py measure LABEL [--jobs 24] [--work DIR]   # one step: runs/LABEL.json
    python scripts/na_fixes_record.py report                                 # report.json and report.md

`measure` draws every entry in the default view with scripts/random100.py's own entry runs (`run --views default
--not-for-paper`) and keeps, per entry, the A17 measures exactly as random100 takes them (nacheck through
random100.na_metrics), the figure's named fallbacks and the view's time; runs/LABEL.json is stamped with the code that
measured (foldmap.source_stamp). `report` reads the runs named in FIXES and writes, per fix, its target before and
after, every control that moved and every measure that got worse; committed and reverted fixes alike. These entries
are regression sets (the failures were diagnosed on them, the controls chosen from them), never unseen ones: the
unseen measurement is the fourth sample at the frozen release (ruling Q37)."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from foldmap.nacheck import _HIGH_IS_GOOD, THRESHOLDS
from foldmap.source_stamp import source_stamp, stamp_line

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "validation" / "na-fixes"
RANDOM100 = ROOT / "scripts" / "random100.py"
# the seven failures of the third random sample (na-diagnosis.md, per-case table)
CASES = ("8VBG", "6VTX", "6ZY7", "3EGZ", "4AS1", "6QIK", "8XJV")
# the 39 controls that passed (na-diagnosis.md, "Fixes"): 7 of the third sample, 8 held-out, 9 tuned, and the 15
# protein-nucleic-acid entries of the NA campaign plus 1ZBI. 1IHF, 1ZAA, 1ZBI, 1ZDI and 9C75 are F3's risk set.
CONTROLS = (
    "6LAB", "9XSF", "3W96", "6JCY", "9OWM", "2FCY", "6IZP", "4NLG", "2VBO", "5K7Z", "9WB3", "2XLJ", "9K4F", "1Z43",
    "9IWF", "9K9V", "3L2V", "3MHT", "9YLM", "7XSP", "9C75", "9CJI", "8BW5", "1LMB", "1ZAA", "1KX5", "1IHF", "1YTB",
    "1RVA", "1HYS", "3KTQ", "4OO8", "30QM", "4LD0", "3CRX", "1URN", "1ZDI", "5W1H", "1ZBI",
)  # fmt: skip
F3_RISK = ("1IHF", "1ZAA", "1ZBI", "1ZDI", "9C75")
# scripts/bench.py --check (heldout-failures bench-before.json) --exact, run at each committed fix (the gate)
BENCH = "every entry identical (the bench's structures hold no case these fixes act on)"
# what stays a limitation (ruling Q37; Q4: 6ZY7)
LIMITS = (
    "F1 (M1, pairs missed at low resolution) deferred past the freeze: 8XJV pair recall, 3EGZ A2.U64 (pair_f1 0.941) "
    "and 6ZY7 stem_recall 0.667 (ruling Q4) stay failures",
    "3EGZ procrustes 0.572 per copy (limit 0.5): the stems of one copy span 49 degrees round the axis, unrolled flat",
    "6QIK (F5 not committed): fallback, crossing and the A3 strand y:8-9 remain",
    "8VBG (F3 not committed): contact_auc 0.505, stems_nearer_frac 0.5 remain; 9YLM, 9K4F, 30QM, 3KTQ, 1YTB as before",
    "8XJV: everything M1 explains (fallbacks 10, na_crossings 10, bead_overlaps 127, time)",
)
# what a control must not get worse on: every judged A17 measure (time apart: it is the machine's), and the crossings
# reported but not judged
JUDGED = tuple(dict.fromkeys([*(k for k in THRESHOLDS if k != "seconds"), "na_run_crossings", "na_protein_crossings"]))
KEPT = JUDGED + (
    "pair_precision",
    "pair_recall",
    "pairs_foldmap",
    "stems_drawn",
    "contact_precision_at_k",
    "stems_nearer_not_applicable",
    "ss_contacts_nearer",
    "bend_3d",
    "bend_drawn",
)
# The fixes in the order made (ruling Q37). before/after: run labels; committed: the commit, or None when reverted
# (its case stays a limitation); target: the entries and measures it was made to move.
FIXES = [
    {
        "fix": "F7",
        "mechanism": "M3: stems_nearer_frac scored where every drawn element touches the stem",
        "change": "nacheck.contacts scores a stem in stems_nearer_frac only when an element does not touch it; such "
        "stems are counted apart (stems_nearer_not_applicable) and reported not applicable (ruling Q5)",
        "before": "base",
        "after": "F7",
        "committed": "24676d3",
        "target": {"6VTX": ("stems_nearer_frac", "stems_nearer_not_applicable")},
    },
    {
        "fix": "F2",
        "mechanism": "M4: unrolled RNA placed by the protomer it faces, measured against one plain view",
        "change": "a CylinderFrame unrolls each nucleic-acid chain one protomer holds (cli._na_holders: at least "
        "NA_HOLD_MIN residues touching it and more than every other protomer together) with that protomer, its stems "
        "tied to protomer 0's copy (layout._stem_ties); the shape "
        "measures take each copy in a plain frame at its own centre, Spearman and Procrustes per copy (the worst "
        "kept), the angles judged against half the turn between neighbouring copies",
        "before": "F7",
        "after": "F2",
        "committed": "36adc3b",
        "target": {"3EGZ": ("spearman", "procrustes", "angle_median", "angle_max", "mirrored", "na_run_crossings")},
    },
    {
        "fix": "F4",
        "mechanism": "M5: a gap meets a free end's tip along the end's own outward heading",
        "change": "the 5' free end after a gap turns toward the gap's start as well as its touching protein; the gap "
        "tries a free arrival along the chord, then the clear straight chord, before an unblocked route; a 3' end the "
        "protein touches keeps lying along it",
        "before": "F2",
        "after": "F4",
        "committed": "3380435",
        "target": {"4AS1": ("na_crossings",), "8XJV": ("na_crossings",)},
    },
    {
        "fix": "F5",
        "mechanism": "M6: a short duplex cut where projection, not 3D, makes the bend",
        "change": "the layout (layout._kept_whole) keeps a duplex under 2 * SEGMENT_MIN pairs whole when its 3D bend "
        "is under 60 degrees but its pieces would be drawn turned more than 90 degrees (Nucleic.whole, read by "
        "stem_pieces)",
        "before": "F4",
        "after": "F5",
        "committed": None,
        "why_not": "it moved its target (6QIK fallbacks 1 -> 0, na_crossings 1 -> 0, bead_overlaps 1 -> 0) but made "
        "control 1RVA worse: its 10-bp EcoRV site (kinked 50 degrees in 3D, under WHOLE_BEND, 104 on the page) was "
        "drawn whole, its kink gone (bend_3d -102 -> None; angle_median, angle_max, joint_error_max 0 -> not "
        "measured); "
        "6QIK "
        "itself fell to stems_nearer_frac 0.992 and gained 100 run x run and 96 NA x protein crossings. No 3D bend "
        "separates 6QIK (52.4) from 1RVA (50): its case stays a limitation",
        "target": {"6QIK": ("fallbacks", "na_crossings")},
    },
    {
        "fix": "F6",
        "mechanism": "M7: bead clearance best effort",
        "change": "a run nucleotide _unbunch leaves short of a bead's width from another run's has the other run slide "
        "(na_labels._clear_junctions, _slide_from); routing runs for (m + 1) * max(NT_STEP, bead) changes nothing, as "
        "NT_STEP (0.55) is wider than a bead (0.45), so it is not made",
        "before": "F4",
        "after": "F6",
        "committed": "218f4bb",
        "target": {"6QIK": ("bead_overlaps",), "8XJV": ("bead_overlaps",)},
    },
    {
        "fix": "F3",
        "mechanism": "M2: a stem inside the protein's projection pushed out by _attract",
        "change": "in a plain Frame, a stem whose projected centre lies inside the hull of its touching elements gets "
        "room before _attract (its overlapping items moved straight away from it) and a 10 times stronger contact pull",
        "before": "F6",
        "after": "F3",
        "committed": None,
        "why_not": "its target barely moved (8VBG contact_auc 0.505 -> 0.521, stems_nearer_frac 0.5 unchanged: of its "
        "two "
        "hairpin pieces only #1 lies inside the hull of its touching elements) while ten controls got worse: the hull "
        "test also takes 1IHF's bent DNA (procrustes 0.274 -> 0.389, contact_auc 0.786 -> 0.647), 3KTQ (contact_auc "
        "0.582 -> 0.296, stems_nearer_frac 1 -> 0), 9K4F, 4NLG, 2VBO, 5W1H (contact_auc), and the stronger pull added "
        "NA x "
        "protein crossings (6JCY, 7XSP, 9CJI, 4OO8). Its case stays a limitation",
        "target": {"8VBG": ("contact_auc", "stems_nearer_frac"), "6QIK": ("na_crossings",)},
    },
]


def _row(e: dict) -> dict:
    """One entry's default view as random100 measured it: the A17 measures (KEPT), what failed, the named fallbacks
    and the view's exit and time."""
    view = e.get("views", {}).get("default") or {}
    res = view.get("result") or {}
    na = res.get("na") or {}
    row = {k: na.get(k) for k in KEPT}
    row["na_failed"] = na.get("failed")
    row["named_fallbacks"] = res.get("fallbacks") or {}
    row["exit"] = view.get("exit")
    row["seconds"] = view.get("seconds")
    row["error"] = e.get("error") or res.get("measure_error")
    return row


def measure(label: str, jobs: int, work: Path) -> dict:
    ids = [*CASES, *CONTROLS]
    work.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-I", str(RANDOM100), "run", "--ids", ",".join(ids), "--views", "default"]
    cmd += ["--jobs", str(jobs), "--out", str(work), "--sample", str(work / "no-sample.json"), "--not-for-paper"]
    subprocess.run(cmd, check=False)
    rows = {}
    for i in ids:
        f = work / f"{i}.json"
        rows[i] = _row(json.loads(f.read_text(encoding="utf-8"))) if f.is_file() else {"error": "no entry file"}
    doc = {"label": label, "source_stamp": SOURCE_STAMP, "entries": rows}
    (OUT / "runs").mkdir(parents=True, exist_ok=True)
    (OUT / "runs" / f"{label}.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return doc


def _num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def worse(key: str, a, b) -> bool:
    """Whether a measure went from a to b the wrong way (a measured value lost counts as worse: the check is gone)."""
    if _num(a) and _num(b):
        return b < a if key in _HIGH_IS_GOOD else b > a
    return _num(a) and b is None


def compare(before: dict, after: dict) -> dict:
    """Per entry: every kept measure that changed (before, after), and the judged ones that got worse, with the named
    fallbacks' total among them."""
    out = {}
    for i in [*CASES, *CONTROLS]:
        a, b = before["entries"].get(i, {}), after["entries"].get(i, {})
        # stems_nearer_not_applicable is new with F7: None before it is 0, not a change
        changed = {k: [a.get(k), b.get(k)] for k in KEPT if a.get(k) != b.get(k) and (a.get(k), b.get(k)) != (None, 0)}
        bad = [k for k in JUDGED if worse(k, a.get(k), b.get(k))]
        fa, fb = sum((a.get("named_fallbacks") or {}).values()), sum((b.get("named_fallbacks") or {}).values())
        if fa != fb:
            changed["named_fallbacks"] = [fa, fb]
            if fb > fa:
                bad.append("named_fallbacks")
        if a.get("error") != b.get("error"):
            changed["error"] = [a.get("error"), b.get("error")]
            if b.get("error"):
                bad.append("error")
        if changed:
            out[i] = {"changed": changed, "worse": bad}
    return out


def report() -> dict:
    runs = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted((OUT / "runs").glob("*.json"))}
    fixes = []
    for f in FIXES:
        if f["before"] not in runs or f["after"] not in runs:
            fixes.append({**f, "measured": False})
            continue
        before, after = runs[f["before"]], runs[f["after"]]
        diff = compare(before, after)
        target = {
            i: {k: [before["entries"][i].get(k), after["entries"][i].get(k)] for k in keys}
            for i, keys in f["target"].items()
        }
        controls = {i: d for i, d in diff.items() if i in CONTROLS}
        fixes.append(
            {
                **f,
                "measured": True,
                "target": target,
                "cases_changed": {i: d for i, d in diff.items() if i in CASES},
                "controls_changed": controls,
                "controls_worse": {i: d["worse"] for i, d in controls.items() if d["worse"]},
                "code": {k: after["source_stamp"].get(k) for k in ("code_fingerprint", "git")},
            }
        )
    doc = {
        "source_stamp": SOURCE_STAMP,
        "sets": {"cases": list(CASES), "controls": list(CONTROLS), "f3_risk": list(F3_RISK)},
        "set_kind": "regression (the failures were diagnosed on these entries); never unseen",
        "judged": list(JUDGED),
        "bench": BENCH,
        "limitations": list(LIMITS),
        "runs": {
            k: {"code": {c: v["source_stamp"].get(c) for c in ("code_fingerprint", "git")}} for k, v in runs.items()
        },
        "fixes": fixes,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8")
    (OUT / "report.md").write_text("\n".join(lines(doc, runs)) + "\n", encoding="utf-8")
    return doc


def _v(x) -> str:
    return "None" if x is None else (f"{x:.3g}" if isinstance(x, float) else str(x))


def lines(doc: dict, runs: dict) -> list[str]:
    out = [
        "# Nucleic-acid fixes, measured",
        "",
        f"Plan 28 (ruling Q37). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        f"Each fix measured before and after on the {len(CASES)} failing cases of the third random sample "
        f"({', '.join(CASES)}) and {len(CONTROLS)} nucleic-acid controls that passed, in the default view, with "
        "scripts/random100.py's entry runs (A17 as random100 takes it). **These are regression sets, not unseen "
        "ones**: the failures were diagnosed on these entries and the fixes made for them. The unseen measurement is "
        "the fourth sample at the frozen release. A fix is committed only when its target moves and no control gets "
        "worse on any judged A17 measure, fallback count, crossing count or bead overlap (ruling Q37).",
        "",
        "Runs (the code each step measured; git is the commit checked out, the step's fix uncommitted on it, so the "
        "code fingerprint names the code: a committed fix's equals its commit's src):",
        "",
    ]
    for k, v in doc["runs"].items():
        g = v["code"].get("git") or {}
        out.append(f"- {k}: code {v['code'].get('code_fingerprint')}, git {str(g.get('commit'))[:12]}")
    for f in doc["fixes"]:
        state = f"committed in {f['committed']}" if f["committed"] else "**not committed** (reverted)"
        out += ["", f"## {f['fix']}: {state}", "", f"- Mechanism: {f['mechanism']}.", f"- Change: {f['change']}."]
        if f.get("why_not"):
            out.append(f"- Not committed: {f['why_not']}. The change as measured: reverted/{f['fix']}.patch.")
        else:
            out.append(f"- Bench: {BENCH}.")
        if not f.get("measured"):
            out.append("- Not measured.")
            continue
        out += [f"- Steps: {f['before']} -> {f['after']}.", "", "Target:", "", "| entry | measure | before | after |"]
        out.append("|---|---|---|---|")
        for i, ks in f["target"].items():
            out += [f"| {i} | {k} | {_v(a)} | {_v(b)} |" for k, (a, b) in ks.items()]
        out += ["", "Cases that moved:", ""]
        out += _changes(f["cases_changed"])
        out += ["", "Controls that moved:", ""]
        out += _changes(f["controls_changed"])
        worse = f["controls_worse"]
        out += ["", "Controls worse: " + ("; ".join(f"{i} ({', '.join(k)})" for i, k in worse.items()) or "none")]
    out += ["", "## Limitations", ""] + [f"- {x}." for x in LIMITS]
    return out


def _changes(d: dict) -> list[str]:
    if not d:
        return ["- none"]
    return [
        f"- {i}: "
        + "; ".join(f"{k} {_v(a)} -> {_v(b)}" for k, (a, b) in x["changed"].items())
        + (f" (worse: {', '.join(x['worse'])})" if x["worse"] else "")
        for i, x in d.items()
    ]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("measure", help="draw and measure every case and control at this code: runs/LABEL.json")
    m.add_argument("label")
    m.add_argument("--jobs", type=int, default=24)
    m.add_argument("--work", type=Path, default=Path("na-fixes-work"))
    sub.add_parser("report", help="report.json and report.md from the runs FIXES names")
    c = sub.add_parser("compare", help="print what moved from run BEFORE to run AFTER, and what got worse")
    c.add_argument("before")
    c.add_argument("after")
    args = ap.parse_args(argv)
    if args.cmd == "measure":
        measure(args.label, args.jobs, args.work)
    elif args.cmd == "compare":
        runs = [json.loads((OUT / "runs" / f"{k}.json").read_text(encoding="utf-8")) for k in (args.before, args.after)]
        diff = compare(*runs)
        for i, d in diff.items():
            print(("CASE " if i in CASES else "CONTROL ") + _changes({i: d})[0][2:])
        print("controls worse:", {i: d["worse"] for i, d in diff.items() if i in CONTROLS and d["worse"]} or "none")
    else:
        report()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
