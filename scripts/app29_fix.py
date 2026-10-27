# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A plan 29 fix's record (plan 29 Task 4 on): the fix measured against the code before it, every measure the plan
names, written as report.json and report.md, stamped with the code that wrote them (foldmap.source_stamp).

    python scripts/app29_fix.py marks --ids A,B --out FILE [--src DIR]
    python scripts/app29_fix.py report --name NAME --before RECORD --after RECORD --runs-before DIR --runs-after DIR
                                       --marks-before FILE --marks-after FILE --pngs DIR [--patch FILE] --out DIR

marks draws each id the app's way (random100.app_settings, app_jobs.draw_settings: no job process, no limits) and
lists every ligand mark: its place, how much of the loops' lines its box holds (`on_loop`, page units: a line under
an opaque mark seems to end at it), whether it is crowded, and each tether's length with its home and target chain.
--src puts that folder's foldmap first (the code before the fix, from a worktree), so the same script measures both.

report reads, read-only:
  - the before and after records of scripts/app29.py run (load_tune: report.json alone; the check half in its totals,
    never report-check.json): each tune entry's failing criteria, its floor, and whether its figure (app.svg's sha256)
    or its page session changed; the check half's totals;
  - the third sample's keys (pass, failing criteria, named fallbacks, area; scripts/prune.py's _campaign_keys and
    campaign_verdict) and the NA cases and controls (prune's na measure, na_verdict), each run by scripts/prune.py run
    --labels baseline, before (RUNS_BEFORE) and after (RUNS_AFTER); the entries whose SVG moved; the third sample's
    committed report (views all, its default view) against the after run, for reference (it is older code);
  - bench --check --exact (bench.txt in each runs folder): the entries whose layout moved, each explained by
    BENCH_WHY or said unexplained;
  - the marks of the tune entries whose figure changed, before and after (marks's output);
  - the PNGs kept for those entries (PNGS: <id>-before.png, <id>-after.png), named by path.
The verdict is the plan's: kept only when every measure is identical or better and no test fails (the suite's result
is the gate's, said in the task's report, not here)."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
THIRD = ROOT / "docs" / "validation" / "random100-third"

# the bench entries a fix moves, each with why (filled once the fix's bench run is read; empty: none moved)
BENCH_WHY: dict[str, str] = {}

# figures whose move no key judges, looked at (plan 29 Task 4 review, rulings P29-9 and P29-10): id -> what is seen
LOOK_CONCERNS: dict[str, str] = {
    "9WNE": "the Mg's two tethers, to R's and W's tips, which lie level, run along one line, one over the other: the "
    "far step's first zero-C9 spot is at angle 0 on its ring, and no check sees tethers drawn along each other "
    "(C9 counts loop lines only; parked by ruling P29-10 for the Task N+1 look test)",
    "1KX5": "the Mn left a loop line for a spot further out; its tether is several times longer (the move below), "
    "won on a C9 tie (0 against 0)",
    "9X4L": "ADP 32001's home tether to the target 15 units off is 13.9 units long, drawn in pieces under alpha17 "
    "and alpha18: ruled policy, home tethers are unpriced (plan 11; A20 bounds tethers to other chains only)",
}

LINE_EPS = 1e-9  # page units: a box holding less of a line than this holds none
C9_EPS = 2e-3  # page units: C9 as the rows keep it (3 decimals per mark, summed): differences within it are rounding


def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"app29_fix_{name}", SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _stamp() -> dict:
    if str(ROOT / "src") not in sys.path and (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    from foldmap.source_stamp import source_stamp

    return source_stamp(__file__)


# ------------------------------------------------------------------------------------------------- marks


def marks_of(i: str) -> dict:
    """Every ligand mark of `i` drawn the app's way (its download in the random-100 cache)."""
    import numpy as np

    from foldmap import app, app_jobs, render
    from foldmap.route import clip_lengths

    random100 = _script("random100")
    settings = json.loads(random100.app_settings())
    path = Path.home() / ".cache" / "foldmap" / f"{i}.cif"  # where app29's run downloaded it (foldmap.fetch)
    drawn = app_jobs.draw_settings(path, settings)[0]
    look = app.spec_of(path, settings).look
    lay, loops = drawn.layout, drawn.loops
    place, off_c9, crossings = render._place_mark, {}, {}

    def placed(
        anchor, half, targets, taken, lines, texts, tethers, bodies=(), runs=(), na=None, retarget=None, bound=()
    ):
        """The real placement; for a mark left on a loop line, also the C9 of the spot clear of the lines further
        out that _place_mark weighed it against (None: there is none), keyed by the mark's place."""
        spot = place(anchor, half, targets, taken, lines, texts, tethers, bodies, runs, na, retarget, bound)
        ts = targets if retarget is None else retarget(spot)
        crossings[(round(float(spot[0]), 3), round(float(spot[1]), 3))] = int(
            sum(render.clip_len(spot, t, b) > 1e-9 for t in ts for b in texts)
        )
        box = (spot[0] - half[0], spot[1] - half[1], spot[0] + half[0], spot[1] + half[1])
        if any(box[0] < t[2] and t[0] < box[2] and box[1] < t[3] and t[1] < box[3] for t in lines):
            far = render._SPOT_FAR * render._SPOT_RADIUS
            o = render._mark_spot(anchor, half, targets, taken + lines, texts, tethers, bodies, far, runs, na, retarget)
            off_c9[(round(float(spot[0]), 3), round(float(spot[1]), 3))] = (
                None
                if o is None
                else round(float(render._over_bodies(o, targets, bodies) + render._along_loops(o, targets, runs)), 3)
            )
        return spot

    render._place_mark = placed
    try:
        marks = render.ligand_marks(lay, loops, drawn.sses, look)
    finally:
        render._place_mark = place
    segs = np.array([(*p, *q) for lp in loops for p, q in zip(lp.points, lp.points[1:])], float).reshape(-1, 4)
    bodies = [p.rect for p in lay.placed.values()] + [g.rect for g in lay.ghosts]  # as ligand_marks keeps them
    rows = []
    for m in marks:
        g, (cx, cy), (hx, hy) = m["ligand"], m["centre"], m["half"]
        box = np.array([[cx - hx, cy - hy, cx + hx, cy + hy]])
        on = float(clip_lengths(segs, box).max()) if len(segs) else 0.0
        chains = [lay.res_chain[r] for r in m["residues"]] if m["residues"] else [g.chain] * len(m["targets"])
        rows.append(
            {
                "mark": f"{g.name}:{g.chain}{g.seq}",
                "centre": [round(float(cx), 3), round(float(cy), 3)],
                "on_loop": round(on, 3),
                "crowded": bool(m["crowded"]),
                "over_bodies": round(float(render._over_bodies(m["centre"], m["targets"], bodies)), 3),  # C9
                "along_loops": round(float(render._along_loops(m["centre"], m["targets"], segs)), 3),  # C9 (5HI7)
                "off_line_c9": off_c9.get((round(float(cx), 3), round(float(cy), 3)), "not on a line"),
                "text_crossings": crossings.get((round(float(cx), 3), round(float(cy), 3))),  # (tether, text) pairs
                "tethers": [
                    [round(float(np.hypot(tx - cx, ty - cy)), 3), chains[0], c]
                    for (tx, ty), c in zip(m["targets"], chains)
                ],
            }
        )
    return {"id": i, "marks": rows}


def cmd_marks(args) -> int:
    if args.src:
        sys.path.insert(0, str(Path(args.src).resolve()))
    elif (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    import foldmap

    out = {"source_stamp": _stamp(), "foldmap": str(Path(foldmap.__file__).resolve().parent), "entries": {}}
    for i in [x.strip() for x in args.ids.split(",") if x.strip()]:
        out["entries"][i] = marks_of(i)["marks"]
        print(i, len(out["entries"][i]), flush=True)
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    return 0


# ------------------------------------------------------------------------------------------------- the record


def _on(rows: list[dict]) -> list[str]:
    return [r["mark"] for r in rows if r["on_loop"] > LINE_EPS]


def _over(rows: list[dict]) -> float:
    return round(sum(r["over_bodies"] for r in rows), 3)


def _along(rows: list[dict]) -> float:
    return round(sum(r["along_loops"] for r in rows), 3)


def _text(rows: list[dict]) -> int:
    return sum(r.get("text_crossings") or 0 for r in rows)


def marks_order(m: dict) -> int:
    """One entry's marks, after against before, in _place_mark's precedence (`chosen`): crowded marks, then tethers
    crossing text (the gate a ranked spot meets), then C9 (beyond the 3 decimals the rows keep), then marks on a loop
    line: 1 worse, -1 better, 0 alike."""
    for a, b, eps in (
        (m["crowded"][0], m["crowded"][1], 0),
        (m["text_crossings"][0], m["text_crossings"][1], 0),
        (m["c9"][0], m["c9"][1], C9_EPS),
        (len(m["on_loop"][0]), len(m["on_loop"][1]), 0),
    ):
        if abs(b - a) > eps:
            return 1 if b > a else -1
    return 0


def _c9(rows: list[dict]) -> float:
    """C9's measure, as render._mark_spot ranks spots: the tethers' length over bodies and along loop lines."""
    return round(_over(rows) + _along(rows), 3)


def _cross(rows: list[dict]) -> list[float]:
    return [t[0] for r in rows for t in r["tethers"] if t[1] != t[2]]


def mark_changes(before: list[dict], after: list[dict]) -> dict:
    """The marks of one entry, before and after: those on a loop line, the longest tether to another chain, and each
    mark that moved (paired in order: ligand_marks's order is the ligands', not their places)."""
    moved = []
    for a, b in zip(before, after):
        if a["centre"] != b["centre"]:
            moved.append(
                {
                    "mark": b["mark"],
                    "from": a["centre"],
                    "to": b["centre"],
                    "on_loop": [a["on_loop"], b["on_loop"]],
                    "over_bodies": [a["over_bodies"], b["over_bodies"]],
                    "along_loops": [a["along_loops"], b["along_loops"]],
                    "longest_tether": [
                        max((t[0] for t in a["tethers"]), default=None),
                        max((t[0] for t in b["tethers"]), default=None),
                    ],  # fmt: skip
                    "tethers": [[t[0] for t in a["tethers"]], [t[0] for t in b["tethers"]]],
                }
            )
    xa, xb = _cross(before), _cross(after)
    return {
        "marks": [len(before), len(after)],
        "over_bodies": [_over(before), _over(after)],
        "along_loops": [_along(before), _along(after)],
        "c9": [_c9(before), _c9(after)],
        "text_crossings": [_text(before), _text(after)],
        "on_loop": [_on(before), _on(after)],
        "crowded": [sum(r["crowded"] for r in before), sum(r["crowded"] for r in after)],
        "longest_cross_chain_tether": [max(xa, default=None), max(xb, default=None)],
        "moved": moved,
    }


def tune_changes(before: dict, after: dict) -> list[dict]:
    """Each tune entry of two app29 records: failing criteria, floor, figure (app.svg) and page session."""
    b = {e["id"]: e for e in before["entries"]}
    out = []
    for e in after["entries"]:
        a = b.get(e["id"])
        if a is None:
            out.append({"id": e["id"], "note": "not in the before record"})
            continue
        sa, sb = (a.get("files") or {}).get("app.svg"), (e.get("files") or {}).get("app.svg")
        ua = (a.get("ui") or {}) if isinstance(a.get("ui"), dict) else {}
        ub = (e.get("ui") or {}) if isinstance(e.get("ui"), dict) else {}
        out.append(
            {
                "id": e["id"],
                "criteria": [a["criteria_failed"], e["criteria_failed"]],
                "floor_ok": [a["floor_ok"], e["floor_ok"]],
                "svg_changed": sa != sb,
                "refused": [a.get("refused"), e.get("refused")],
                "error": [a.get("error"), e.get("error")],
                "ui_errors": [len(ua.get("errors") or []), len(ub.get("errors") or [])],
                "ui_exception": [ua.get("exception"), ub.get("exception")],
                "seconds": [a.get("seconds"), e.get("seconds")],
            }
        )
    return out


def _worse_tune(rows: list[dict]) -> list[str]:
    bad = []
    for r in rows:
        if "criteria" not in r:
            bad.append(f"{r['id']}: {r['note']}")
            continue
        new = sorted(set(r["criteria"][1]) - set(r["criteria"][0]))
        if new:
            bad.append(f"{r['id']}: now fails {', '.join(new)}")
        if r["floor_ok"][0] and not r["floor_ok"][1]:
            bad.append(f"{r['id']}: now under the floor")
        if not r["refused"][0] and r["refused"][1]:
            bad.append(f"{r['id']}: now refused")
        if not r["error"][0] and r["error"][1]:
            bad.append(f"{r['id']}: now an error")
        if r["ui_errors"][1] > r["ui_errors"][0] or (r["ui_exception"][1] and not r["ui_exception"][0]):
            bad.append(f"{r['id']}: its page session now errs")
    return bad


def _better_tune(rows: list[dict]) -> list[str]:
    good = []
    for r in rows:
        if "criteria" not in r:
            continue
        gone = sorted(set(r["criteria"][0]) - set(r["criteria"][1]))
        if gone:
            good.append(f"{r['id']}: now passes {', '.join(gone)}")
        if not r["floor_ok"][0] and r["floor_ok"][1]:
            good.append(f"{r['id']}: now above the floor")
    return good


def check_totals(before: dict, after: dict) -> dict:
    a, b = before["totals"]["check"], after["totals"]["check"]
    worse = [
        k for k in ("refused", "error", "not_measured", "png_error", "criteria_failed", "floor_failed") if b[k] > a[k]
    ]
    worse += ["drawn"] if b["drawn"] < a["drawn"] else []
    better = [
        k for k in ("refused", "error", "not_measured", "png_error", "criteria_failed", "floor_failed") if b[k] < a[k]
    ]
    return {"before": a, "after": b, "worse": worse, "better": better}


def bench_notes(path: Path) -> dict:
    """bench.txt (bench.py --check --exact's output, then "bench exit N"): each entry's note and the exit."""
    notes, code = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("bench exit "):
            code = int(line.split()[-1])
            continue
        parts = line.split()
        if parts and ("identical" in line or "changed, area" in line or "MISSING" in line):  # an entry's line
            notes[parts[0]] = "identical" if line.rstrip().endswith("identical") else line.split("  ")[-1].strip()
    return {"exit": code, "notes": notes}


def bench_diff(before: dict, after: dict) -> dict:
    changed = sorted(n for n, v in after["notes"].items() if v != before["notes"].get(n))
    return {
        "exit": [before["exit"], after["exit"]],
        "entries": len(after["notes"]),
        "changed": [{"entry": n, "note": after["notes"][n], "why": BENCH_WHY.get(n, "unexplained")} for n in changed],
        "worse": [n for n in changed if n not in BENCH_WHY],
    }


def build(args) -> dict:
    app29 = _script("app29")
    prune = _script("prune")
    rd = _script("random100")
    sens = _script("sensitivity")
    before, after = app29.load_tune(args.before), app29.load_tune(args.after)
    tune = tune_changes(before, after)
    rb, ra = Path(args.runs_before), Path(args.runs_after)
    third_b = prune._campaign_keys(rb / "r100", "baseline", "r100")
    third_a = prune._campaign_keys(ra / "r100", "baseline", "r100")
    third = prune.campaign_verdict(third_b["entries"], third_a["entries"])
    third["svg_moved"] = prune.svg_moved(third_b["svg"], third_a["svg"])
    third["entries"] = len(third_a["entries"])
    third["passed"] = [sens.passed(third_b["entries"]), sens.passed(third_a["entries"])]
    committed = sens.keys_of(rd.read_report(THIRD / "report.json"), "default")
    ref = prune.campaign_verdict(committed, third_a["entries"])
    third["committed_reference"] = {
        "note": "the committed report (views all, an older code) against the after run, default view: for reference",
        "passed": [sens.passed(committed), sens.passed(third_a["entries"])],
        "worse": len(ref["worse"]),
        "better": len(ref["better"]),
        "moved": len(ref["moved"]),
    }
    na_b = prune._campaign_keys(rb / "na", "baseline", "na")
    na_a = prune._campaign_keys(ra / "na", "baseline", "na")
    na = prune.na_verdict({"entries": na_b["entries"]}, {"entries": na_a["entries"]})
    na["svg_moved"] = prune.svg_moved(na_b["svg"], na_a["svg"])
    na["entries"] = len(na_a["entries"])
    bench = bench_diff(bench_notes(rb / "bench.txt"), bench_notes(ra / "bench.txt"))
    mb = json.loads(Path(args.marks_before).read_text(encoding="utf-8"))
    ma = json.loads(Path(args.marks_after).read_text(encoding="utf-8"))
    changed = [r["id"] for r in tune if r.get("svg_changed")]
    marks = {i: mark_changes(mb["entries"].get(i, []), ma["entries"].get(i, [])) for i in changed}
    both = set(mb["entries"]) & set(ma["entries"])  # every entry measured both ways
    every = [r["id"] for r in tune if r["id"] in both]  # the tune half's

    def total(ids: list[str]) -> dict:
        return {
            "entries": len(ids),
            "marks": [sum(len(mb["entries"][i]) for i in ids), sum(len(ma["entries"][i]) for i in ids)],
            "on_loop": [sum(len(_on(mb["entries"][i])) for i in ids), sum(len(_on(ma["entries"][i])) for i in ids)],
            "crowded": [
                sum(r["crowded"] for i in ids for r in mb["entries"][i]),
                sum(r["crowded"] for i in ids for r in ma["entries"][i]),
            ],  # fmt: skip
            "over_bodies": [
                round(sum(_over(mb["entries"][i]) for i in ids), 3),
                round(sum(_over(ma["entries"][i]) for i in ids), 3),
            ],  # fmt: skip
            "text_crossings": [sum(_text(mb["entries"][i]) for i in ids), sum(_text(ma["entries"][i]) for i in ids)],
            "c9": [
                round(sum(_c9(mb["entries"][i]) for i in ids), 3),
                round(sum(_c9(ma["entries"][i]) for i in ids), 3),
            ],
        }

    marks_total = total(every)
    marks_total["moved_in_unchanged_figures"] = [
        i for i in every if i not in changed and mb["entries"][i] != ma["entries"][i]
    ]
    elsewhere = sorted({*third["svg_moved"], *na["svg_moved"]} & both)  # third sample and NA figures that moved
    marks_elsewhere = {i: mark_changes(mb["entries"][i], ma["entries"][i]) for i in elsewhere}
    marks_elsewhere_total = total(elsewhere)
    pngs = {
        i: {
            k: f"{args.pngs}/{i}-{k}.png"
            for k in ("before", "after", "pair")
            if (Path(args.pngs) / f"{i}-{k}.png").is_file()
        }
        for i in changed
    }
    worse = {
        "tune": _worse_tune(tune),
        "check": check_totals(before, after)["worse"],
        "third": third["worse"],
        "na": na["worse"],
        "bench": bench["worse"],
        "marks": [i for i, m in {**marks, **marks_elsewhere}.items() if marks_order(m) > 0],
    }
    better = {
        "tune": _better_tune(tune),
        "check": check_totals(before, after)["better"],
        "third": third["better"],
        "na": na["better"],
        "marks": [i for i, m in {**marks, **marks_elsewhere}.items() if marks_order(m) < 0],
    }
    patch = None
    if args.patch:
        p = Path(args.patch)
        patch = {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
    return {
        "name": args.name,
        "inputs": {
            "before_record": str(args.before),
            "before_stamp": before.get("source_stamp"),
            "after_record": str(args.after),
            "after_stamp": after.get("source_stamp"),
            "runs_before_stamp": third_b.get("report_stamp"),
            "runs_after_stamp": third_a.get("report_stamp"),
            "marks_before_foldmap": mb.get("foldmap"),
            "marks_before_stamp": mb.get("source_stamp"),
            "marks_after_foldmap": ma.get("foldmap"),
            "marks_after_stamp": ma.get("source_stamp"),
        },
        "tune": {"entries": len(tune), "svg_changed": changed, "rows": tune},
        "check": check_totals(before, after),
        "third": third,
        "na": na,
        "bench": bench,
        "marks": marks,
        "marks_total": marks_total,
        "marks_elsewhere": marks_elsewhere,
        "marks_elsewhere_total": marks_elsewhere_total,
        "not_changed": [
            {
                "id": i,
                "mark": r["mark"],
                "on_line_c9": round(r["over_bodies"] + r["along_loops"], 3),
                "off_line_c9": r.get("off_line_c9"),
            }
            for i in every
            for r in ma["entries"][i]
            if r["on_loop"] > LINE_EPS
        ],
        "look_concerns": [
            {
                "id": i,
                "seen": LOOK_CONCERNS[i],
                "moves": [
                    {
                        "mark": x["mark"],
                        "moved": round(float(sum((p - q) ** 2 for p, q in zip(x["from"], x["to"])) ** 0.5), 3),
                        "tethers": x["tethers"],
                    }
                    for x in (marks.get(i) or marks_elsewhere.get(i) or {"moved": []})["moved"]
                ],
            }
            for i in LOOK_CONCERNS
            if i in marks or i in marks_elsewhere
        ],
        "since": since_last(args, after, third_a, na_a, ma, prune),
        "pngs": pngs,
        "worse": worse,
        "better": better,
        "verdict": "worse" if any(worse.values()) else "better" if any(better.values()) else "identical",
        "patch": patch,
    }


def since_last(args, after: dict, third_a: dict, na_a: dict, ma: dict, prune) -> dict | None:
    """The figures that changed since the fix's last record (`--since-*`: that record's tune run, its third-sample and
    NA runs, and the marks its code drew for the figures that moved), each with its cause, read from its marks: which
    mark moved, on or off a loop line, its tethers crossing text and its C9, then -> now."""
    if not args.since_runs:
        return None
    prev = {e["id"]: e for e in json.loads(Path(args.since_record).read_text(encoding="utf-8"))["entries"]}
    ids = {
        e["id"]
        for e in after["entries"]
        if (e.get("files") or {}).get("app.svg") != (prev.get(e["id"], {}).get("files") or {}).get("app.svg")
    }
    for m, now in (("r100", third_a), ("na", na_a)):
        ids |= set(prune.svg_moved(prune._campaign_keys(Path(args.since_runs) / m, "baseline", m)["svg"], now["svg"]))
    mp = json.loads(Path(args.since_marks).read_text(encoding="utf-8"))["entries"]
    rows = []
    for i in sorted(ids):
        was, now = mp.get(i), ma["entries"].get(i)
        if was is None or now is None:
            rows.append({"id": i, "cause": ["not measured"], "text_crossings": None, "c9": None})
            continue
        c = mark_changes(was, now)
        why = []
        for a, b in zip(was, now):
            if a["centre"] == b["centre"]:
                continue
            off, on = a["on_loop"] > LINE_EPS and b["on_loop"] <= LINE_EPS, a["on_loop"] <= LINE_EPS < b["on_loop"]
            fewer = (b.get("text_crossings") or 0) < (a.get("text_crossings") or 0)
            if off and fewer:
                why.append(f"{b['mark']}: its spot on the line was _mark_spot's fallback (a tether crossed text); the "
                           "ranked spot off the line wins (fix round 1, I1)")  # fmt: skip
            elif on and fewer:
                why.append(
                    f"{b['mark']}: the spot off the line it took at {args.since_label} was _mark_spot's fallback (a "
                    "tether crossed text); the ranked spot on the line wins, as before the fix (fix round 1, I1)"
                )
            else:
                why.append(f"{b['mark']}: moved with the marks above (their boxes and tethers keep it off)")
        rows.append({"id": i, "cause": why, "text_crossings": c["text_crossings"], "c9": c["c9"]})
    return {"label": args.since_label, "figures": rows}


def _pair(x) -> str:
    return f"{x[0]} -> {x[1]}"


def render(doc: dict) -> str:
    from foldmap.source_stamp import stamp_line

    lines = [f"# Plan 29 fix record: {doc['name']}", "", f"Written by: {stamp_line(doc['source_stamp'])}", ""]
    lines += [f"**Verdict: {doc['verdict']}** (identical or better on every measure keeps the fix).", ""]
    for k, v in doc["worse"].items():
        lines.append(f"- worse, {k}: {', '.join(map(str, v)) if v else 'none'}")
    for k, v in doc["better"].items():
        lines.append(f"- better, {k}: {', '.join(map(str, v)) if v else 'none'}")
    lines += ["", "## Tune half (per entry; the check half in its totals only)", ""]
    lines += ["| id | criteria failed | floor ok | svg changed | seconds |", "|---|---|---|---|---|"]
    for r in doc["tune"]["rows"]:
        if "criteria" not in r:
            lines.append(f"| `{r['id']}` | {r['note']} | | | |")
            continue
        crit = _pair([",".join(c) or "-" for c in r["criteria"]])
        svg = "yes" if r["svg_changed"] else "no"
        lines.append(f"| `{r['id']}` | {crit} | {_pair(r['floor_ok'])} | {svg} | {_pair(r['seconds'])} |")
    c = doc["check"]
    lines += ["", "## Check half (totals)", "", "| total | before | after |", "|---|---|---|"]
    for k in c["before"]:
        lines.append(f"| {k} | {c['before'][k]} | {c['after'][k]} |")
    t = doc["third"]
    lines += ["", "## Third sample, default view (keys: pass, failing criteria, named fallbacks, area)", ""]
    lines.append(f"- {t['entries']} entries; pass {_pair(t['passed'])}; verdict {t['verdict']}")
    for k in ("worse", "better", "moved"):
        lines.append(f"- {k}: {'; '.join(t[k]) if t[k] else 'none'}")
    lines.append(f"- SVG moved (no key reads it): {', '.join(t['svg_moved']) if t['svg_moved'] else 'none'}")
    r = t["committed_reference"]
    lines.append(
        f"- {r['note']}: pass {_pair(r['passed'])}; {r['worse']} worse, {r['better']} better, {r['moved']} moved"
    )
    n = doc["na"]
    lines += [
        "",
        "## NA cases and controls (plan 28's NA measure)",
        "",
        f"- {n['entries']} entries; verdict {n['verdict']}",
    ]
    for k in ("worse", "better", "moved"):
        lines.append(f"- {k}: {'; '.join(n[k]) if n[k] else 'none'}")
    lines.append(f"- SVG moved (no key reads it): {', '.join(n['svg_moved']) if n['svg_moved'] else 'none'}")
    b = doc["bench"]
    lines += ["", "## bench --check --exact", "", f"- {b['entries']} entries; exit {_pair(b['exit'])}"]
    if not b["changed"]:
        lines.append("- changed: none")
    for x in b["changed"]:
        lines.append(f"- `{x['entry']}`: {x['note']}; why: {x['why']}")
    t = doc["marks_total"]
    lines += ["", "## Ligand marks", ""]
    lines.append(
        f"Over the {t['entries']} tune entries measured both ways: marks {_pair(t['marks'])}, on a loop line "
        f"{_pair(t['on_loop'])}, crowded {_pair(t['crowded'])}, tethers crossing text "
        f"{_pair(t['text_crossings'])}; C9 (the tethers' length over bodies and along loop "
        f"lines, page units: render._mark_spot's ranking) {_pair(t['c9'])}, of it over bodies "
        f"{_pair(t['over_bodies'])}; marks moved in a figure that did not change: "
        f"{', '.join(t['moved_in_unchanged_figures']) or 'none'}."
    )
    e = doc["marks_elsewhere_total"]
    lines += [
        "",
        f"The third sample's and the NA figures that moved ({e['entries']}, drawn the app's way): marks "
        f"{_pair(e['marks'])}, on a loop line {_pair(e['on_loop'])}, crowded {_pair(e['crowded'])}, tethers "
        "crossing text "
        f"{_pair(e['text_crossings'])}; C9 "
        f"{_pair(e['c9'])}, of it over bodies {_pair(e['over_bodies'])}.",
        "",
    ]
    for i, m in doc["marks_elsewhere"].items():
        lines.append(
            f"- `{i}`: on a loop line {_pair([len(x) for x in m['on_loop']])}; tethers crossing text "
            f"{_pair(m['text_crossings'])}; C9 {_pair(m['c9'])} (over bodies "
            f"{_pair(m['over_bodies'])}, along loop lines {_pair(m['along_loops'])}); {len(m['moved'])} moved"
        )
    lines += ["", "The tune entries whose figure changed:", ""]
    lines += [
        "| id | marks | on a loop line | tethers crossing text | C9 | longest tether to another chain | moved |",
        "|---|---|---|---|---|---|---|",
    ]
    for i, m in doc["marks"].items():
        on, c9, cross = _pair([len(x) for x in m["on_loop"]]), _pair(m["c9"]), _pair(m["longest_cross_chain_tether"])
        tx = _pair(m["text_crossings"])
        lines.append(f"| `{i}` | {_pair(m['marks'])} | {on} | {tx} | {c9} | {cross} | {len(m['moved'])} |")
    lines += [
        "",
        "Each moved mark: its place, the line it holds, its tethers over bodies and along loop lines (page units), and "
        "its tethers' lengths, before -> after.",
        "",
    ]
    for i, m in doc["marks"].items():
        for x in m["moved"]:
            lines.append(
                f"- `{i}` {x['mark']}: {x['from']} -> {x['to']}; on a loop line {_pair(x['on_loop'])}; over bodies "
                f"{_pair(x['over_bodies'])}; along loop lines {_pair(x['along_loops'])}; tethers {x['tethers'][0]} -> "
                f"{x['tethers'][1]}"
            )
    lines += ["", "## Not changed, by design", ""]
    lines.append(
        "Tune marks still on a loop line: the spot clear of the lines further out (to _SPOT_FAR x _SPOT_RADIUS) would "
        "run its tethers over more bodies and along more loop lines (C9) than the spot on the line, or there is none "
        "('none'). C9 in page units, on the line -> the best off the line:"
    )
    lines.append("")
    for x in doc["not_changed"]:
        off = "none" if x["off_line_c9"] is None else x["off_line_c9"]
        lines.append(f"- `{x['id']}` {x['mark']}: C9 {x['on_line_c9']} on the line, {off} off it")
    lines += [
        "",
        "## Precedence, step by step",
        "",
        "- `chosen` (a spot found after the first layer, against the spot clear of the lines further out): a spot "
        "that met the text and reach gates beats one returned for want of any; then C9 (within 1e-9, a tie to the "
        "spot off the lines); then the line. A bound target left beyond _TETHER_REACH by the spot off the lines, "
        "and kept within it by the other, keeps the other (A20 before a line).",
        "- `kept` (the bound's reach): the layers' order, the line first; then C9 within a layer.",
        "",
        "## Look concerns (no key judges them)",
        "",
    ]
    for x in doc["look_concerns"]:
        moves = "; ".join(
            f"{m['mark']} moved {m['moved']}, tethers {m['tethers'][0]} -> {m['tethers'][1]}" for m in x["moves"]
        )
        lines.append(f"- `{x['id']}`: {x['seen']}. {moves}")
    if doc.get("since"):
        lines += ["", f"## Changed since {doc['since']['label']} (the last record), and why", ""]
        lines.append(
            "Each figure that moved since that record, with its tethers crossing text and its C9 then -> now. C9 can "
            "rise while the verdict stays better: _mark_spot's own order puts the text gate (a spot whose tethers "
            "cross no text is ranked) before C9, and these marks lost a tether crossing text."
        )
        lines.append("")
        for x in doc["since"]["figures"]:
            nums = (
                ""
                if x["c9"] is None
                else f" Tethers crossing text {_pair(x['text_crossings'])}; C9 {_pair(x['c9'])}"
                + (" (C9 rises: the text gate comes first)." if x["c9"][1] > x["c9"][0] + C9_EPS else ".")
            )
            lines.append(f"- `{x['id']}`: {'; '.join(x['cause'])}.{nums}")
    lines += ["", "## Pictures (gitignored, for looking)", ""]
    for i, p in doc["pngs"].items():
        lines.append(f"- `{i}`: " + ", ".join(f"{k} {v}" for k, v in p.items()))
    if doc["patch"]:
        lines += ["", f"Patch: {doc['patch']['path']} (sha256 {doc['patch']['sha256'][:12]})"]
    return "\n".join(lines) + "\n"


def cmd_report(args) -> int:
    doc = {"source_stamp": _stamp(), **build(args)}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    (out / "report.md").write_text(render(doc), encoding="utf-8")
    print(f"{doc['name']}: {doc['verdict']}; worse {sum(map(len, doc['worse'].values()))}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("marks", help="each id's ligand marks, drawn the app's way")
    m.add_argument("--ids", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--src", help="a folder holding the foldmap package to measure (default: this checkout's)")
    r = sub.add_parser("report", help="the fix's record")
    r.add_argument("--name", required=True)
    for k in ("before", "after", "runs-before", "runs-after", "marks-before", "marks-after", "pngs", "out"):
        r.add_argument(f"--{k}", required=True)
    r.add_argument("--patch", help="the fix as a patch, kept when the fix is reverted")
    r.add_argument("--since-label", help="the fix's last record, by its commit (a fix round's record)")
    r.add_argument("--since-record", help="that record's tune run report.json")
    r.add_argument("--since-runs", help="that record's runs folder (r100, na)")
    r.add_argument("--since-marks", help="the marks that record's code drew (marks --src) for the figures that moved")
    args = parser.parse_args(argv)
    return {"marks": cmd_marks, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
