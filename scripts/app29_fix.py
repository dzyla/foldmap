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

LINE_EPS = 1e-9  # page units: a box holding less of a line than this holds none


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
    segs = np.array([(*p, *q) for lp in loops for p, q in zip(lp.points, lp.points[1:])], float).reshape(-1, 4)
    bodies = [p.rect for p in lay.placed.values()] + [g.rect for g in lay.ghosts]  # as ligand_marks keeps them
    rows = []
    for m in render.ligand_marks(lay, loops, drawn.sses, look):
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

    out = {"foldmap": str(Path(foldmap.__file__).resolve().parent), "entries": {}}
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
        "marks_on_loop": [i for i, m in marks.items() if len(m["on_loop"][1]) > len(m["on_loop"][0])],
        "marks_crowded": [i for i, m in marks.items() if m["crowded"][1] > m["crowded"][0]],
        "marks_over_bodies": [i for i, m in marks.items() if m["over_bodies"][1] > m["over_bodies"][0]],
        "marks_elsewhere": [
            i
            for i, m in marks_elsewhere.items()
            if len(m["on_loop"][1]) > len(m["on_loop"][0])
            or m["crowded"][1] > m["crowded"][0]
            or m["over_bodies"][1] > m["over_bodies"][0]
        ],  # fmt: skip
    }
    better = {
        "tune": _better_tune(tune),
        "check": check_totals(before, after)["better"],
        "third": third["better"],
        "na": na["better"],
        "marks_on_loop": [i for i, m in marks.items() if len(m["on_loop"][1]) < len(m["on_loop"][0])],
        "marks_over_bodies": [i for i, m in marks.items() if m["over_bodies"][1] < m["over_bodies"][0]],
    }
    patch = None
    if args.patch:
        p = Path(args.patch)
        patch = {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
    return {
        "name": args.name,
        "source_stamp": _stamp(),
        "inputs": {
            "before_record": str(args.before),
            "before_stamp": before.get("source_stamp"),
            "after_record": str(args.after),
            "after_stamp": after.get("source_stamp"),
            "runs_before_stamp": third_b.get("report_stamp"),
            "runs_after_stamp": third_a.get("report_stamp"),
            "marks_before_foldmap": mb.get("foldmap"),
            "marks_after_foldmap": ma.get("foldmap"),
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
        "pngs": pngs,
        "worse": worse,
        "better": better,
        "verdict": "worse" if any(worse.values()) else "better" if any(better.values()) else "identical",
        "patch": patch,
    }


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
        f"{_pair(t['on_loop'])}, crowded {_pair(t['crowded'])}, tethers over bodies (C9, page units) "
        f"{_pair(t['over_bodies'])}; marks moved in a figure that did not change: "
        f"{', '.join(t['moved_in_unchanged_figures']) or 'none'}."
    )
    e = doc["marks_elsewhere_total"]
    lines += ["", f"The third sample's and the NA figures that moved ({e['entries']}, drawn the app's way): marks "
              f"{_pair(e['marks'])}, on a loop line {_pair(e['on_loop'])}, crowded {_pair(e['crowded'])}, tethers "
              f"over bodies {_pair(e['over_bodies'])}."]  # fmt: skip
    for i, m in doc["marks_elsewhere"].items():
        lines.append(
            f"- `{i}`: on a loop line {_pair([len(x) for x in m['on_loop']])}, over bodies {_pair(m['over_bodies'])}, "
            f"{len(m['moved'])} moved"
        )
    lines += ["", "The tune entries whose figure changed:", ""]
    lines += ["| id | marks | on a loop line | longest tether to another chain | moved |", "|---|---|---|---|---|"]
    for i, m in doc["marks"].items():
        on = _pair([len(x) for x in m["on_loop"]])
        lines.append(
            f"| `{i}` | {_pair(m['marks'])} | {on} | {_pair(m['longest_cross_chain_tether'])} | {len(m['moved'])} |"
        )
    lines += [
        "",
        "Each moved mark: its place, the line it holds (page units) and its tethers' lengths, before -> after.",
        "",
    ]
    for i, m in doc["marks"].items():
        for x in m["moved"]:
            lines.append(
                f"- `{i}` {x['mark']}: {x['from']} -> {x['to']}; on a loop line {_pair(x['on_loop'])}; "
                f"over bodies {_pair(x['over_bodies'])}; tethers {x['tethers'][0]} -> {x['tethers'][1]}"
            )
    lines += ["", "## Pictures (gitignored, for looking)", ""]
    for i, p in doc["pngs"].items():
        lines.append(f"- `{i}`: " + ", ".join(f"{k} {v}" for k, v in p.items()))
    if doc["patch"]:
        lines += ["", f"Patch: {doc['patch']['path']} (sha256 {doc['patch']['sha256'][:12]})"]
    return "\n".join(lines) + "\n"


def cmd_report(args) -> int:
    doc = build(args)
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
    args = parser.parse_args(argv)
    return {"marks": cmd_marks, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
