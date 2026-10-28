# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A plan 29 fix's short record under ruling P29-13 (the lighter loop, Tasks 6 on): a fix that must change no figure,
measured on the tune half and the bench only, plus its own measure, written as report.json and report.md, stamped
with the code that wrote them (foldmap.source_stamp).

    python scripts/app29_light.py helical --old FILE --out FILE [--jobs N] [--max-copies 64]
    python scripts/app29_light.py report --name NAME --before RUN --after RUN --bench FILE --measure FILE --out DIR

helical (Task 6) runs symmetry._helical as it was (--old: that symmetry.py, e.g. `git show <commit>:src/foldmap/
symmetry.py`) and as it is on the same _Assembly of every structure: every .cif/.cif.gz under tests/data, plan 29's
tune half and every cached download (~/.cache/foldmap/*.cif) less plan 29's check half, which it never names. It
compares the two returns exactly (floats bit for bit), counts _Assembly.pair calls and times each, with both
MIN_HELIX_COPIES and MIN_HELIX_COPIES_ASKED. A structure whose largest cluster holds more than --max-copies copies is
skipped for time (the old search is n^4). Exit 1 when any return differs.

report reads two scripts/app29.py run records of the tune half (load_tune: report.json alone; RUN, a folder with it),
the bench's --check --exact output (bench.txt, then "bench exit N") and the fix's own measure (--measure, a JSON
file this script's other commands wrote; it is copied into the record), and says: every tune entry whose figure (the
app's SVG and PNG, sha256) or criteria changed, the draw seconds of the entries the fix names (NAMED), the bench's
verdict, and the measure's summary. The verdict is P29-13's: kept only when every figure is byte-identical and the
bench is identical on every entry."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "docs" / "validation" / "app29" / "sample.json"

# each fix: what it is, the entries whose draw seconds the record names, and what it was expected to change
FIXES: dict[str, dict] = {
    "fix-screw": {
        "title": "Plan 29 Task 6: the screw search without n^4 pair fits (triage proposal 3, ruling P29-7)",
        "mechanism": "symmetry._helical tried every pair of copies as the screw and tested each on every copy against "
        "every member of its cluster: n^4/2 pair fits (23FR's 42 copies: 1.26 million _Assembly.pair calls). A fit's "
        "RMSD is at least the distance between the moved and the target centroids of its matched atoms, so a member "
        "whose centroid lands farther than the tolerance (plus 1e-6 Å for rounding) is passed over unfitted; the "
        "first member within the tolerance, in cluster order, is the one the full search took.",
        "named": ["23FR"],
        "expected": "23FR's time only; every figure byte-identical",
    },
}


def _stamp() -> dict:
    if str(ROOT / "src") not in sys.path and (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    from foldmap.source_stamp import source_stamp

    return source_stamp(__file__)


def _halves() -> tuple[list[str], set[str]]:
    split = json.loads(SAMPLE.read_text(encoding="utf-8"))["split"]
    return list(split["tune"]), set(split["check"])


def _module(path: Path, name: str):
    """A foldmap module as it was (its file), imported inside the foldmap package so its relative imports resolve."""
    import foldmap  # noqa: F401

    spec = importlib.util.spec_from_file_location(f"foldmap.{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def structures() -> dict[str, str]:
    """path -> where it comes from: tests/data, the tune half (fetched), the download cache (less the check half)."""
    from foldmap.fetch import fetch

    tune, check = _halves()
    paths: dict[str, str] = {}
    for p in sorted((ROOT / "tests" / "data").rglob("*")):
        if p.is_file() and p.name.endswith((".cif", ".cif.gz")):
            paths[str(p)] = "tests/data"
    for i in tune:
        paths.setdefault(str(fetch(i)), "tune")
    for p in sorted((Path.home() / ".cache" / "foldmap").glob("*.cif")):
        if p.stem.upper() not in check:
            paths.setdefault(str(p), "cache")
    return paths


# ------------------------------------------------------------------------------------------------ helical


def _same_helical(a, b) -> bool:
    import numpy as np

    if a is None or b is None:
        return a is None and b is None
    (sa, ta, axa, ca, ra, cha, rma), (sb, tb, axb, cb, rb, chb, rmb) = a, b
    return (sa, ta, ra, cha, rma) == (sb, tb, rb, chb, rmb) and np.array_equal(axa, axb) and np.array_equal(ca, cb)


def _counting(asm) -> list[int]:
    calls, real = [0], asm.pair

    def pair(a, b):
        calls[0] += 1
        return real(a, b)

    asm.pair = pair
    return calls


def helical_one(path: str, old_file: str, max_copies: int) -> dict:
    from foldmap import symmetry as new
    from foldmap.io import load_backbone, read_model

    old = _module(Path(old_file), "_symmetry_before")
    name = Path(path).name.split(".")[0]
    try:
        bb = load_backbone(path, "auto", model=read_model(path, "auto"))
    except Exception as e:  # not a structure foldmap reads: said, not measured
        return {"name": name, "skipped": f"{type(e).__name__}: {e}"[:200]}
    asm_new, asm_old = new._Assembly(bb), new._Assembly(bb)
    largest = max((len(c) for c in asm_new.clusters), default=0)
    out: dict = {"name": name, "chains": len(asm_new.ca), "largest_cluster": largest}
    if largest > max_copies:
        out["skipped"] = f"largest cluster {largest} > {max_copies}"
        return out
    for mc in (new.MIN_HELIX_COPIES, new.MIN_HELIX_COPIES_ASKED):
        n_calls, o_calls = _counting(asm_new), _counting(asm_old)
        t0 = time.perf_counter()
        b_new = new._helical(asm_new, new.TOLERANCE, mc)
        t1 = time.perf_counter()
        b_old = old._helical(asm_old, new.TOLERANCE, mc)
        t2 = time.perf_counter()
        out[f"min{mc}"] = {
            "same": _same_helical(b_new, b_old),
            "copies": len(b_new[5][0]) if b_new is not None else None,
            "pair_calls_before": o_calls[0],
            "pair_calls_after": n_calls[0],
            "seconds_before": round(t2 - t1, 3),
            "seconds_after": round(t1 - t0, 3),
        }
        del asm_new.pair, asm_old.pair  # the class's own pair again for the next count
    return out


def cmd_helical(args) -> int:
    paths = structures()
    rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = [(pool.submit(helical_one, p, args.old, args.max_copies), src) for p, src in paths.items()]
        for f, src in futs:
            rows.append({**f.result(), "source": src})
    runs = [(r, r[k]) for r in rows for k in ("min4", "min3") if k in r]
    summary = {
        "structures": len(rows),
        "measured": sum(1 for r in rows if "min4" in r),
        "skipped": sorted(f"{r['name']}: {r['skipped']}" for r in rows if "skipped" in r),
        "differ": sorted({r["name"] for r, m in runs if not m["same"]}),
        "helical_found": {k: sum(1 for r in rows if r.get(k, {}).get("copies")) for k in ("min4", "min3")},
        "pair_calls": [sum(m["pair_calls_before"] for _, m in runs), sum(m["pair_calls_after"] for _, m in runs)],
        "seconds": [
            round(sum(m["seconds_before"] for _, m in runs), 2),
            round(sum(m["seconds_after"] for _, m in runs), 2),
        ],
        "more_pair_calls": sorted({r["name"] for r, m in runs if m["pair_calls_after"] > m["pair_calls_before"]}),
        "slowest_by": max(
            ((round(m["seconds_after"] - m["seconds_before"], 3), r["name"]) for r, m in runs), default=None
        ),
    }
    doc = {"kind": "helical", "source_stamp": _stamp(), "summary": summary, "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 1 if summary["differ"] else 0


def _helical_lines(m: dict, named: list[str]) -> list[str]:
    s = m["summary"]
    lines = [
        f"- {s['measured']} structures measured of {s['structures']} (tests/data, the tune half, the download cache "
        f"less the check half), {len(s['skipped'])} skipped: the return of _helical before and after differs on "
        f"{len(s['differ'])}{': ' + ', '.join(s['differ']) if s['differ'] else ''} (MIN_HELIX_COPIES and "
        "MIN_HELIX_COPIES_ASKED; every float compared bit for bit)",
        f"- a screw found on {s['helical_found']['min4']} (4 copies or more) and {s['helical_found']['min3']} (3 or "
        "more)",
        f"- _Assembly.pair calls, all structures: {s['pair_calls'][0]:,} before, {s['pair_calls'][1]:,} after; "
        f"_helical's seconds: {s['seconds'][0]} before, {s['seconds'][1]} after (run in parallel: indicative)",
        f"- more pair calls after on: {', '.join(s['more_pair_calls']) or 'none'} (a screw passes, so every copy's "
        "images are gathered at once)",
    ]
    rows = {r["name"]: r for r in m["rows"] if r["source"] in ("tune", "tests/data")}
    top = sorted((r for r in rows.values() if "min4" in r), key=lambda r: r["min4"]["pair_calls_before"], reverse=True)[
        :6
    ]
    lines += [
        "",
        "| entry | where | largest cluster | screw (copies) | pair calls | _helical s |",
        "|---|---|---|---|---|---|",
    ]
    for r in sorted(top, key=lambda r: r["name"] not in named):
        x = r["min4"]
        lines.append(
            f"| {r['name']} | {r['source']} | {r['largest_cluster']} | {x['copies'] or '-'} | "
            f"{x['pair_calls_before']:,} -> {x['pair_calls_after']:,} | {x['seconds_before']} -> {x['seconds_after']} |"
        )
    return lines


# ------------------------------------------------------------------------------------------------ report


def bench(path: Path) -> dict:
    """bench.txt: each entry's note ("identical" or what changed) and the exit."""
    notes, code = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("bench exit "):
            code = int(line.split()[-1])
            continue
        parts = line.split()
        if parts and ("identical" in line or "changed" in line or "MISSING" in line):
            notes[parts[0]] = "identical" if line.rstrip().endswith("identical") else line.split("  ")[-1].strip()
    return {
        "exit": code,
        "entries": len(notes),
        "identical": sum(v == "identical" for v in notes.values()),
        "notes": notes,
    }


def tune(before: dict, after: dict, named: list[str]) -> dict:
    a = {e["id"]: e for e in before["entries"] if e.get("half") == "tune"}
    b = {e["id"]: e for e in after["entries"] if e.get("half") == "tune"}
    common = sorted(a.keys() & b.keys())

    def moved(i, f):
        return a[i].get("files", {}).get(f) != b[i].get("files", {}).get(f)

    return {
        "entries": [len(a), len(b)],
        "missing": sorted(a.keys() ^ b.keys()),
        "svg_changed": [i for i in common if moved(i, "app.svg")],
        "png_changed": [i for i in common if moved(i, "app.png")],
        "criteria_changed": [
            i
            for i in common
            if (a[i]["criteria_failed"], a[i]["floor_ok"], a[i]["error"], a[i]["refused"])
            != (b[i]["criteria_failed"], b[i]["floor_ok"], b[i]["error"], b[i]["refused"])
        ],
        "named_seconds": {
            i: {"draw": [a[i]["draw_seconds"], b[i]["draw_seconds"]], "entry": [a[i]["seconds"], b[i]["seconds"]]}
            for i in named
            if i in a and i in b
        },
        "draw_seconds_total": [
            round(sum(a[i]["draw_seconds"] or 0 for i in common), 2),
            round(sum(b[i]["draw_seconds"] or 0 for i in common), 2),
        ],
        "stamps": [before.get("source_stamp", {}).get("git"), after.get("source_stamp", {}).get("git")],
    }


def cmd_report(args) -> int:
    fix = FIXES[args.name]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    measure = json.loads(Path(args.measure).read_text(encoding="utf-8"))
    kept = out / f"{measure['kind']}.json"
    if Path(args.measure).resolve() != kept.resolve():
        shutil.copyfile(args.measure, kept)
    load = lambda d: json.loads((Path(d) / "report.json").read_text(encoding="utf-8"))  # noqa: E731
    t = tune(load(args.before), load(args.after), fix["named"])
    b = bench(Path(args.bench))
    identical = not (t["svg_changed"] or t["png_changed"] or t["missing"])
    bench_ok = b["exit"] == 0 and b["identical"] == b["entries"]
    own_ok = not measure["summary"].get("differ")
    verdict = "kept" if identical and bench_ok and own_ok else "rejected"
    rep = {
        "name": args.name,
        "source_stamp": _stamp(),
        "fix": fix,
        "before": str(args.before),
        "after": str(args.after),
        "tune": t,
        "bench": b,
        "measure": {"file": kept.name, "kind": measure["kind"], "summary": measure["summary"]},
        "verdict": verdict,
    }
    (out / "report.json").write_text(json.dumps(rep, indent=1) + "\n", encoding="utf-8")
    (out / "report.md").write_text(render(rep, measure), encoding="utf-8")
    print(f"{args.name}: {verdict}")
    return 0


def render(rep: dict, measure: dict) -> str:
    t, b, fix = rep["tune"], rep["bench"], rep["fix"]
    st = rep["source_stamp"].get("git", {})
    lines = [
        f"# {fix['title']}",
        "",
        f"Written by scripts/app29_light.py report at {st.get('commit', '?')[:7]} "
        f"({'dirty' if st.get('dirty') else 'clean'}); ruling P29-13's light loop: the tune half and the bench, the "
        "check half, the third sample and the NA controls left to Task N+1.",
        "",
        f"**Verdict: {rep['verdict']}.** Expected to change: {fix['expected']}.",
        "",
        "## Mechanism",
        "",
        fix["mechanism"],
        "",
        "## The tune half (scripts/app29.py run, the app's way)",
        "",
        f"- before: `{rep['before']}` (git {(t['stamps'][0] or {}).get('commit', '?')[:7]}); after: `{rep['after']}` "
        f"(git {(t['stamps'][1] or {}).get('commit', '?')[:7]})",
        f"- entries {t['entries'][0]} -> {t['entries'][1]}; missing: {', '.join(t['missing']) or 'none'}",
        f"- the app's SVG changed on: {', '.join(t['svg_changed']) or 'none'}; its PNG on: "
        f"{', '.join(t['png_changed']) or 'none'}",
        f"- criteria, floor, error or refusal changed on: {', '.join(t['criteria_changed']) or 'none'}",
        f"- draw seconds, all tune entries: {t['draw_seconds_total'][0]} -> {t['draw_seconds_total'][1]} (8 entries "
        "drawn at once: indicative)",
    ]
    for i, s in t["named_seconds"].items():
        lines.append(
            f"- {i}: draw {s['draw'][0]} -> {s['draw'][1]} s; the entry's process (draw, PNG, checks) "
            f"{s['entry'][0]} -> {s['entry'][1]} s"
        )
    lines += [
        "",
        "## bench --check --exact",
        "",
        f"- {b['identical']} of {b['entries']} identical, exit {b['exit']}",
    ]
    lines += [f"- {n}: {v}" for n, v in sorted(b["notes"].items()) if v != "identical"]
    lines += ["", f"## The fix's own measure ({rep['measure']['file']})", ""]
    if measure["kind"] == "helical":
        lines += _helical_lines(measure, fix["named"])
    lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("helical", help="symmetry._helical before and after on every structure at hand")
    h.add_argument("--old", required=True, help="symmetry.py as it was")
    h.add_argument("--out", required=True)
    h.add_argument("--jobs", type=int, default=24)
    h.add_argument("--max-copies", type=int, default=64)
    r = sub.add_parser("report", help="the short record")
    r.add_argument("--name", required=True, choices=sorted(FIXES))
    r.add_argument("--before", required=True)
    r.add_argument("--after", required=True)
    r.add_argument("--bench", required=True)
    r.add_argument("--measure", required=True)
    r.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    return {"helical": cmd_helical, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
