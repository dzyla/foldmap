# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""How fragile foldmap's tuned constants are (plan 28 Task 9.2, ruling Q29): each module constant the inventory
(docs/validation/inventory/inventory.md) gives the origin "measured" or "named-entry" is perturbed to x0.75 and x1.25,
one at a time, and the random-100 campaign keys are measured on the third random sample, which none of them was
tuned on (docs/validation/random100-third; never the fourth sample, whose protocol any look voids).

    python scripts/sensitivity.py plan                              # the constants, perturbations and stage A subset
    python scripts/sensitivity.py run --stage A --out DIR [--labels baselines|perturbations|L1,L2] [--jobs 8]
                                      [--parallel 4]
    python scripts/sensitivity.py run --stage B --out DIR                # the whole sample, every perturbation (Q30)
    python scripts/sensitivity.py report --a DIR [--b DIR]          # docs/validation/sensitivity/report.{md,json}

A perturbation is the value times the factor, an int rounded half away from zero; it is skipped and listed when that
gives the same value or one that cannot be meant (an int count of 1 or more falling under 1). A tuple, a bool or any
other non-number is listed as not perturbable. Each run is `random100.py run --not-for-paper --views default` in a
child process with FOLDMAP_PERTURB=module:NAME=VALUE: scripts/perturb.py, installed by random100.py before it imports
foldmap, rewrites that module as it is imported, so the run's every process (the run, each entry, each view) draws
and judges with the value, and so does every reader of it (from-import binders, defaults and constants derived from
it at import); `perturb.py reach` checks it after each run. Runs write to DIR (outside docs/validation/random100-
third); a run whose keys file exists is not run again, so a stopped stage resumes.

Stage A: a fixed stratified subset of the third sample (SEED; units are entries, a pair is one unit; an entry whose
recorded time in the third run exceeds SLOW_S is left out, listed), the default view, two baselines (they must be
identical: the measurement is deterministic), a reads run and every perturbation. Stage B: the whole third sample,
default view, the same unperturbed runs and every perturbation (ruling Q30; Q29 ran only those that moved in A). The
reads run is unperturbed with FOLDMAP_TRACE_READS set (perturb.install_trace): which entries ran a line naming each
constant, so a stable constant no line read is told apart (its stability is no evidence); it must move no key.

Keys, per run: each entry's pass or fail; the failing count per A- and S-criterion; each entry's named fallbacks
(default view); each entry's figure area; time is recorded, not judged, and so are the time criteria (TIMED: A2, the
default view's wall time within 60 s, moved with the load of the machine between two identical baselines; an entry
passes here when no other A-criterion fails). A perturbation is stable when no key moves
from the baseline, sensitive when any does, fragile when the pass count moves by FRAGILE_PASS or more or a criterion
that fails in no baseline entry fails in some entry; a constant takes the worse of its two perturbations, from stage B
where it ran there. Fragile constants are the limitations list for the paper."""

from __future__ import annotations

import argparse
import datetime
import functools
import hashlib
import importlib
import importlib.util
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    sys.path.insert(0, str(ROOT / "src"))  # the checkout's foldmap, as random100.py

from foldmap.source_stamp import source_stamp, stamp_line  # noqa: E402

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts (plan 28 Task 1)
SCRIPTS = ROOT / "scripts"
THIRD = ROOT / "docs" / "validation" / "random100-third"
INVENTORY = ROOT / "docs" / "validation" / "inventory" / "inventory.md"
OUT = ROOT / "docs" / "validation" / "sensitivity"
ENV = "FOLDMAP_PERTURB"  # scripts/perturb.py's
FACTORS = (0.75, 1.25)  # ruling Q29: +-25 %
ORIGINS = ("measured", "named-entry")  # ruling Q29: the inventory origins perturbed
SEED = 20261025  # stage A's subset draw (ruling Q29: seed recorded)
SUBSET = 40  # stage A's entries, fixed before the baselines ran
SLOW_S = 120  # ruling Q29: an entry recorded slower than this is left out of stage A
FRAGILE_PASS = 2  # ruling Q29: a pass count moving by this many is fragile
TIMED = ("A2", "A2b")  # time criteria (the default view's and the app's wall time): recorded, not judged (Q29)
BASELINES = ("baseline-1", "baseline-2")
READS = "reads"  # an unperturbed run recording which constants each entry reads (perturb.install_trace)
TRACE_ENV = "FOLDMAP_TRACE_READS"  # scripts/perturb.py's
CLASSES = ("stable", "sensitive", "fragile")


def _squash(text: str, most: int = 80) -> str:
    return text if len(text) <= most else text[: most - 1] + "…"


@functools.cache  # loaded once (the runs' threads share it)
def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"foldmap_{name}", SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------- what is perturbed


def constants(md: Path = INVENTORY) -> list[dict]:
    """The inventory's module constants whose origin is measured or named-entry: [{"id", "module", "name",
    "origin", "source"}] (source: the value as the code writes it)."""
    inv = _script("inventory")
    rows = inv.read_table(md.read_text(encoding="utf-8"))
    out = []
    for i, r in rows.items():
        if r["kind"] != "constant" or inv.origin_kind(r["origin"]) not in ORIGINS:
            continue
        module, name = i.split(":", 1)
        out.append({"id": i, "module": module, "name": name, "origin": r["origin"], "source": r["value"]})
    return out


def perturbed(value, factor: float) -> tuple[int | float | None, str | None]:
    """(the perturbed value, None) or (None, why it is not perturbed)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None, f"not perturbable: a {type(value).__name__}"
    raw = value * factor
    if isinstance(value, int):
        new = int(math.copysign(math.floor(abs(raw) + 0.5), raw))
        if value >= 1 and new < 1:
            return None, f"an int count {value} would fall to {new}"
    else:
        new = round(raw, 12)  # no float noise in the literal (1.1 * 1.25)
    if new == value:
        return None, f"x{factor} gives the same value ({new})"
    return new, None


def label(cid: str, factor: float) -> str:
    return f"{cid.replace(':', '.')}x{factor}"


def plan(items: list[dict] | None = None) -> dict:
    """Each constant's runtime value (the foldmap imported here, never perturbed) and its perturbations."""
    items = constants() if items is None else items
    runs, skipped, fixed = [], [], []
    for c in items:
        value = getattr(importlib.import_module(f"foldmap.{c['module']}"), c["name"])
        c = {**c, "value": value if isinstance(value, (int, float)) else _squash(repr(value))}
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            fixed.append({**c, "why": f"a {type(value).__name__}: no single factor"})
            continue
        for f in FACTORS:
            new, why = perturbed(value, f)
            if why:
                skipped.append({**c, "factor": f, "why": why})
            else:
                text = f"{c['module']}:{c['name']}={new!r}"
                runs.append({**c, "factor": f, "new": new, "label": label(c["id"], f), "perturb": text})
    return {"constants": len(items), "perturbations": runs, "skipped": skipped, "not_perturbable": fixed}


# ------------------------------------------------------------------------------------------------- the subset


def recorded_seconds(report: dict) -> dict[str, float]:
    """Each entry's recorded time in a random-100 report (its entry process's wall seconds)."""
    return {e["id"]: float((e.get("process") or {}).get("seconds") or 0) for e in report["entries"]}


def units(sample: dict) -> list[list[dict]]:
    """The sample's entries as units: a pair (seed and partner) is one unit, any other entry one alone."""
    entries = list(sample["entries"]) + [
        {"id": r["id"], "stratum": "A1c", "refusal": True}
        for r in sample.get("replaced", [])
        if r["replacement"] is None
    ]
    by_id = {e["id"]: e for e in entries}
    out, seen = [], set()
    for e in entries:
        if e["id"] in seen:
            continue
        unit = [e] + ([by_id[e["partner"]]] if e.get("partner") in by_id else [])
        seen.update(x["id"] for x in unit)
        out.append(unit)
    return out


def subset(sample: dict, seconds: dict[str, float], size: int = SUBSET, seed: int = SEED) -> dict:
    """Stage A's entries: per stratum a share of `size` in proportion to its eligible entries (largest remainder,
    at least one unit per stratum), drawn in the order random.Random(f"{seed}:{stratum}") shuffles its units (sorted
    by first id) into; a unit any entry of which was recorded slower than SLOW_S is left out."""
    slow = {i: s for i, s in seconds.items() if s > SLOW_S}
    eligible: dict[str, list[list[dict]]] = {}
    for u in units(sample):
        if not any(e["id"] in slow for e in u):
            eligible.setdefault(u[0]["stratum"], []).append(u)
    total = sum(len(u) for us in eligible.values() for u in us)
    share = {s: size * sum(len(u) for u in us) / total for s, us in eligible.items()}
    quota = {s: math.floor(v) for s, v in share.items()}
    for s in sorted(share, key=lambda s: (-(share[s] - quota[s]), s))[: size - sum(quota.values())]:
        quota[s] += 1
    picked = []
    for s in sorted(eligible):
        us = sorted(eligible[s], key=lambda u: u[0]["id"])
        random.Random(f"{seed}:{s}").shuffle(us)
        got = 0
        for u in us:
            if got >= max(1, quota[s]):
                break
            picked += [e["id"] for e in u]
            got += len(u)
    left_out = {e["id"]: seconds[e["id"]] for u in units(sample) for e in u if any(x["id"] in slow for x in u)}
    return {"seed": seed, "size": size, "ids": picked, "quota": quota, "left_out": left_out}


# ------------------------------------------------------------------------------------------------- keys


def keys_of(report: dict) -> dict[str, dict]:
    """The campaign keys of every entry of a random-100 report: pass, the failing A- and S-criteria, the default
    view's named fallbacks and figure area, and (recorded, not judged) its seconds, its failing time criteria
    (TIMED: a wall-clock limit moves with the machine's load, not the code) and SVG hash. pass is no A-criterion but
    a timed one failing."""
    out = {}
    for e in report["entries"]:
        crit, sci = e.get("criteria") or {}, e.get("science") or {}
        d = ((e.get("views") or {}).get("default") or {}).get("result") or {}
        judged = [k for k, c in crit.items() if c.get("pass") is False and k not in TIMED]
        out[e["id"]] = {
            "pass": not judged,
            "failing": sorted(judged + [k for k, c in sci.items() if c.get("pass") is False]),
            "fallbacks": dict(sorted((d.get("fallbacks") or {}).items())),
            "area": d.get("area"),
            "seconds": (e.get("process") or {}).get("seconds"),
            "timed_failing": sorted(k for k, c in crit.items() if c.get("pass") is False and k in TIMED),
            "svg_sha256": d.get("svg_sha256"),
        }
    return out


def failing_counts(entries: dict[str, dict]) -> Counter:
    return Counter(k for e in entries.values() for k in e["failing"])


def passed(entries: dict[str, dict]) -> int:
    return sum(e["pass"] for e in entries.values())


def moved(base: dict[str, dict], other: dict[str, dict]) -> list[str]:
    """Every key that differs between two runs' entries, one line each (time and SVG hashes not judged)."""
    out = []
    for i in sorted(set(base) | set(other)):
        a, b = base.get(i), other.get(i)
        if a is None or b is None:
            out.append(f"{i}: in one run only")
            continue
        if a["pass"] != b["pass"]:
            out.append(f"{i}: {'PASS' if a['pass'] else 'FAIL'} -> {'PASS' if b['pass'] else 'FAIL'}")
        for k in sorted(set(a["failing"]) ^ set(b["failing"])):
            out.append(f"{i}: {k} {'now fails' if k in b['failing'] else 'now passes'}")
        for n in sorted(set(a["fallbacks"]) | set(b["fallbacks"])):
            x, y = a["fallbacks"].get(n, 0), b["fallbacks"].get(n, 0)
            if x != y:
                out.append(f"{i}: fallback {n} {x} -> {y}")
        if a["area"] != b["area"]:
            out.append(f"{i}: area {a['area']} -> {b['area']}")
    return out


def judge(base: dict[str, dict], run: dict[str, dict]) -> dict:
    """A perturbation's class against its stage's baseline (ruling Q29), with what moved."""
    lines = moved(base, run)
    was, now = failing_counts(base), failing_counts(run)
    new_failing = sorted(k for k in now if now[k] and not was.get(k))
    delta = passed(run) - passed(base)
    cls = "fragile" if abs(delta) >= FRAGILE_PASS or new_failing else "sensitive" if lines else "stable"
    fb_was = Counter(n for e in base.values() for n, c in e["fallbacks"].items() for _ in range(c))
    fb_now = Counter(n for e in run.values() for n, c in e["fallbacks"].items() for _ in range(c))
    return {
        "class": cls,
        "passed": passed(run),
        "pass_delta": delta,
        "new_failing": new_failing,
        "criteria": {k: [was.get(k, 0), now.get(k, 0)] for k in sorted(set(was) | set(now)) if was[k] != now[k]},
        "fallbacks": {
            n: [fb_was.get(n, 0), fb_now.get(n, 0)] for n in sorted(set(fb_was) | set(fb_now)) if fb_was[n] != fb_now[n]
        },
        "entries_moved": sorted({line.split(":", 1)[0] for line in lines}),
        "area_moved": sorted(line.split(":", 1)[0] for line in lines if ": area " in line),
        "moved": lines,
    }


def worst(classes) -> str:
    return max(classes, key=CLASSES.index, default="stable")


# ------------------------------------------------------------------------------------------------- running


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_traces(folder: Path, report: dict) -> dict[str, list[str]]:
    """Which entries read each constant, from the processes' records (perturb.install_trace): a view's or an entry's
    process counts for its entry, the run's own (the judging) as "(judging)"."""
    ids = {e.get("name") or e["id"]: e["id"] for e in report["entries"]}
    out: dict[str, set[str]] = {}
    for p in sorted(folder.glob("*.json")):
        if p.name == "names.json":
            continue
        doc = json.loads(p.read_text(encoding="utf-8"))
        argv = doc["argv"]
        if len(argv) > 3 and argv[1] == "view" and argv[2] == "--name":
            who = ids.get(argv[3], argv[3])
        elif len(argv) > 2 and argv[1] == "entry":
            who = argv[2]
        else:
            who = "(judging)"
        for c in doc["reads"]:
            out.setdefault(c, set()).add(who)
    return {c: sorted(v) for c, v in sorted(out.items())}


def run_one(
    name: str,
    perturb: str | None,
    ids: list[str] | None,
    out: Path,
    jobs: int,
    sample: Path,
    trace: list[str] | None = None,
) -> dict:
    """One random-100 run (default view, --not-for-paper) with `perturb` set, then its keys and reach; with `trace`
    (constant ids), unperturbed, which entries read each (perturb.install_trace). The run's figures and sheets are
    deleted, its report kept."""
    where = out / "runs" / name
    if where.exists():
        shutil.rmtree(where)  # a run stopped half way: run it again whole
    where.mkdir(parents=True)
    cmd = [sys.executable, "-I", str(SCRIPTS / "random100.py"), "run", "--not-for-paper", "--views", "default"]
    cmd += ["--jobs", str(jobs), "--sample", str(sample), "--out", str(where)]
    cmd += ["--ids", ",".join(ids)] if ids else []
    env = {k: v for k, v in os.environ.items() if k not in (ENV, TRACE_ENV)}
    if perturb:
        env[ENV] = perturb
    if trace:
        (where / "reads").mkdir()
        (where / "reads" / "names.json").write_text(json.dumps(trace) + "\n", encoding="utf-8")
        env[TRACE_ENV] = str(where / "reads")
    t0 = time.perf_counter()
    with open(where / "sensitivity.log", "w", encoding="utf-8") as log:
        code = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    seconds = round(time.perf_counter() - t0, 1)
    rep = _script("random100").read_report(where / "report.json") if (where / "report.json").is_file() else None
    reach = None
    if perturb:
        got = subprocess.run(
            [sys.executable, "-I", str(SCRIPTS / "perturb.py"), "reach"], env=env, capture_output=True, text=True
        )
        reach = json.loads(got.stdout) if got.returncode in (0, 1) and got.stdout.strip() else {"error": got.stderr}
    shutil.rmtree(where / "work", ignore_errors=True)
    for p in [*where.glob("*.png"), where / "index.html"]:
        p.unlink(missing_ok=True)
    doc = {
        "source_stamp": SOURCE_STAMP,  # the code that ran and keyed this run (sensitivity.py)
        "run": name,
        "perturb": perturb,
        "ids": ids,
        "exit": code,
        "seconds": seconds,
        "reach": reach,
        "report_stamp": (rep or {}).get("source_stamp"),
        "report_perturbed": (rep or {}).get("perturbed"),
        "mixed_code": (rep or {}).get("mixed_code"),
        "entries": keys_of(rep) if rep else None,
        "reads": read_traces(where / "reads", rep) if trace and rep else None,
    }
    (out / "keys").mkdir(parents=True, exist_ok=True)
    (out / "keys" / f"{name}.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    return doc


def _third() -> tuple[dict, dict]:
    sample = json.loads((THIRD / "sample.json").read_text(encoding="utf-8"))
    return sample, _script("random100").read_report(THIRD / "report.json")


def stage_runs(stage: str, the_plan: dict) -> list[tuple[str, str | None]]:
    """(name, FOLDMAP_PERTURB) of every run of stage A or B: the baselines and the reads run first, then every
    perturbation (stage B: all of them, ruling Q30, not only those that moved in stage A)."""
    runs = [(b, None) for b in (*BASELINES, READS)]
    return runs + [(p["label"], p["perturb"]) for p in the_plan["perturbations"]]


def read_keys(where: Path) -> dict[str, dict]:
    """A stage's keys files, each run's entries keyed again (keys_of) from its kept report: the keys are this code's,
    whichever code wrote the keys file."""
    out = {}
    for p in sorted((where / "keys").glob("*.json")):
        doc = json.loads(p.read_text(encoding="utf-8"))
        rep = where / "runs" / p.stem / "report.json"
        if rep.is_file():
            doc["entries"] = keys_of(_script("random100").read_report(rep))
        out[p.stem] = doc
    return out


def moved_constants(where: Path, the_plan: dict) -> set[str]:
    """The constants either perturbation of which moved any key in a stage."""
    keys = read_keys(where)
    base = keys[BASELINES[0]]["entries"]
    out = set()
    for p in the_plan["perturbations"]:
        k = keys.get(p["label"])
        if k is None:
            raise SystemExit(f"{where}: stage A has no run {p['label']}; finish it first")
        if k["entries"] is None or moved(base, k["entries"]):
            out.add(p["id"])
    return out


def cmd_plan(args) -> int:
    the_plan = plan()
    sample, rep = _third()
    sub = subset(sample, recorded_seconds(rep), args.size)
    default_s = {e["id"]: ((e.get("views") or {}).get("default") or {}).get("seconds") or 0 for e in rep["entries"]}
    print(f"{the_plan['constants']} constants, {len(the_plan['perturbations'])} perturbations")
    for p in the_plan["perturbations"]:
        print(f"  {p['label']}: {p['value']} -> {p['new']}")
    for s in the_plan["skipped"]:
        print(f"  skipped {s['id']} x{s['factor']}: {s['why']}")
    for s in the_plan["not_perturbable"]:
        print(f"  not perturbable {s['id']} ({s['value']}): {s['why']}")
    print(f"stage A: {len(sub['ids'])} entries {sub['ids']}; quota {sub['quota']}; left out {sub['left_out']}")
    picked = sum(default_s[i] for i in sub["ids"])
    print(
        f"recorded default-view seconds (third run, every view drawn): subset {picked:.0f},"
        f" whole sample {sum(default_s.values()):.0f}"
    )
    return 0


def cmd_run(args) -> int:
    out = Path(args.out)
    if out.resolve().is_relative_to(THIRD.resolve()):
        raise SystemExit(f"refused: {out} is inside {THIRD}; sensitivity runs write elsewhere")
    the_plan = plan()
    sample, rep = _third()
    out.mkdir(parents=True, exist_ok=True)
    ids = None
    if args.stage == "A":  # the subset, fixed in the stage's folder by its first run
        sub = subset(sample, recorded_seconds(rep), args.size)
        held = out / "subset.json"
        if held.is_file() and json.loads(held.read_text(encoding="utf-8"))["ids"] != sub["ids"]:
            raise SystemExit(f"{held} holds another subset (another --size?)")
        held.write_text(json.dumps({**sub, "source_stamp": SOURCE_STAMP}, indent=1) + "\n", encoding="utf-8")
        ids = sub["ids"]
    runs = stage_runs(args.stage, the_plan)
    if args.labels == "baselines":
        runs = [r for r in runs if r[0] in (*BASELINES, READS)]
    elif args.labels == "perturbations":
        runs = [r for r in runs if r[0] not in (*BASELINES, READS)]
    elif args.labels:
        wanted = set(args.labels.split(","))
        runs = [r for r in runs if r[0] in wanted]
    todo = [r for r in runs if not (out / "keys" / f"{r[0]}.json").is_file()]
    print(f"stage {args.stage}: {len(runs)} runs, {len(todo)} to run", flush=True)
    t0 = time.perf_counter()

    def one(r):
        trace = sorted({p["id"] for p in the_plan["perturbations"]}) if r[0] == READS else None
        doc = run_one(r[0], r[1], ids, out, args.jobs, THIRD / "sample.json", trace)
        n = doc["entries"]
        print(f"{r[0]}: {passed(n) if n else '-'} pass, {doc['seconds']} s, exit {doc['exit']}", flush=True)

    with ThreadPoolExecutor(max_workers=max(1, args.parallel)) as pool:
        list(pool.map(one, todo))
    with open(out / "stages.jsonl", "a", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {"stage": args.stage, "runs": [r[0] for r in todo], "seconds": round(time.perf_counter() - t0, 1)}
            )
            + "\n"
        )
    return 0


# ------------------------------------------------------------------------------------------------- the report


def _stage(where: Path, the_plan: dict) -> dict:
    """A stage's runs judged: determinism (the two baselines) and every perturbation against baseline-1."""
    keys = read_keys(where)
    base = keys[BASELINES[0]]["entries"]
    other = keys[BASELINES[1]]["entries"]
    svg = sorted(i for i in base if base[i]["svg_sha256"] != (other.get(i) or {}).get("svg_sha256"))
    timing = [json.loads(line) for line in (where / "stages.jsonl").read_text(encoding="utf-8").splitlines() if line]
    out = {
        "entries": sorted(base),
        "baseline_passed": passed(base),
        "baseline_failing": dict(sorted(failing_counts(base).items())),
        "baselines_timed_failing": [  # recorded, not judged (TIMED)
            sorted(i for i, e in run.items() if e.get("timed_failing")) for run in (base, other)
        ],
        "determinism": {"identical": not moved(base, other) and not svg, "moved": moved(base, other), "svg": svg},
        "seconds": round(sum(t["seconds"] for t in timing), 1),
        "run_seconds": round(sum(k["seconds"] for k in keys.values()), 1),
        "runs": {},
    }
    traced = keys.get(READS)
    if traced is not None:  # the reads run draws unperturbed: it must change no key either
        out["reads"] = traced["reads"] or {}
        out["reads_run_moved"] = moved(base, traced["entries"]) if traced["entries"] else ["no report"]
    for p in the_plan["perturbations"]:
        k = keys.get(p["label"])
        if k is None:
            continue
        if k["entries"] is None:
            got = {"class": "fragile", "error": f"no report (exit {k['exit']})", "moved": []}
        else:
            got = judge(base, k["entries"])
        got["reached"] = bool((k.get("reach") or {}).get("reached"))
        got["seconds"] = k["seconds"]
        out["runs"][p["label"]] = got
    return out


def build(a: Path, b: Path | None, stamp: dict) -> dict:
    """The report of stage A (folder `a`) and B (`b`, or None before it runs), stamped with `stamp`."""
    the_plan = plan()
    stage_a = _stage(a, the_plan)
    stage_b = _stage(b, the_plan) if b else None
    rows = {}
    for p in the_plan["perturbations"]:
        row = rows.setdefault(p["id"], {"id": p["id"], "origin": p["origin"], "value": p["value"], "runs": {}})
        got = {"new": p["new"], "A": stage_a["runs"].get(p["label"])}
        if stage_b and p["label"] in stage_b["runs"]:
            got["B"] = stage_b["runs"][p["label"]]
        row["runs"][f"x{p['factor']}"] = got
    for row in rows.values():
        in_b = all("B" in r for r in row["runs"].values())
        stage = "B" if in_b else "A"
        judged = [r[stage]["class"] for r in row["runs"].values() if r.get(stage)]
        row["stage"] = stage
        row["class"] = worst(judged) if len(judged) == len(row["runs"]) else None  # None: a run missing
        for name, st in (("A", stage_a), ("B", stage_b)):
            if st and "reads" in st:  # the entries a line naming it ran in (default view), the judging apart
                who = st["reads"].get(row["id"], [])
                row.setdefault("read_in", {})[name] = len([w for w in who if w != "(judging)"])
                row.setdefault("read_judging", {})[name] = "(judging)" in who
    counts = Counter(r["class"] for r in rows.values())
    unread = [r["id"] for r in rows.values() if r["class"] == "stable" and not (r.get("read_in") or {}).get(r["stage"])]
    unread = [i for i in unread if not (rows[i].get("read_judging") or {}).get(rows[i]["stage"])]
    return {
        "what": "foldmap's measured and named-entry constants perturbed x0.75 and x1.25 on the third random sample "
        "(scripts/sensitivity.py, plan 28 Task 9.2, ruling Q29)",
        "date": datetime.date.today().isoformat(),
        "source_stamp": stamp,
        "sample": "docs/validation/random100-third/sample.json",
        "sample_sha256": _sha256(THIRD / "sample.json"),
        "inventory_sha256": _sha256(INVENTORY),
        "factors": list(FACTORS),
        "rules": {"slow_s": SLOW_S, "fragile_pass": FRAGILE_PASS, "subset": SUBSET, "seed": SEED},
        "subset": {
            k: v for k, v in json.loads((a / "subset.json").read_text(encoding="utf-8")).items() if k != "source_stamp"
        },
        "stage_a": stage_a,
        "stage_b": stage_b,
        "counts": {c: counts.get(c, 0) for c in CLASSES} | {"incomplete": counts.get(None, 0)},
        "constants": list(rows.values()),
        "limitations": [r["id"] for r in rows.values() if r["class"] == "fragile"],
        "stable_unread": unread,  # stable, and no line naming it ran: no evidence either way
        "no_perturbation": sorted({s["id"] for s in the_plan["skipped"]} - set(rows)),
        "skipped": [{k: s[k] for k in ("id", "value", "factor", "why")} for s in the_plan["skipped"]],
        "not_perturbable": [{k: s[k] for k in ("id", "value", "why")} for s in the_plan["not_perturbable"]],
    }


def _said(run: dict | None) -> str:
    if not run:
        return "-"
    if run.get("error"):
        return f"**{run['class']}**: {run['error']}"
    parts = [run["class"] if run["class"] != "fragile" else "**fragile**"]
    if run["pass_delta"]:
        parts.append(f"pass {run['pass_delta']:+d}")
    if run["new_failing"]:
        parts.append("new failing " + ", ".join(run["new_failing"]))
    crit = [f"{k} {a}->{b}" for k, (a, b) in run["criteria"].items() if k not in run["new_failing"]]
    parts += crit
    parts += [f"{n} {a}->{b}" for n, (a, b) in run["fallbacks"].items()]
    if run["area_moved"]:
        parts.append(f"area {len(run['area_moved'])} entr{'y' if len(run['area_moved']) == 1 else 'ies'}")
    if run["entries_moved"]:  # some keys move without a count changing (one entry's S5 fails, another's passes)
        parts.append(f"{len(run['entries_moved'])} entr{'y' if len(run['entries_moved']) == 1 else 'ies'} moved")
    if not run.get("reached", True):
        parts.append("NOT REACHED")
    return "; ".join(parts)


def render(doc: dict) -> str:
    a, b = doc["stage_a"], doc["stage_b"]
    hours = (a["seconds"] + (b["seconds"] if b else 0)) / 3600
    lines = [
        "# Sensitivity of the tuned constants",
        "",
        f"Plan 28 Task 9.2 (ruling Q29). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "Every module constant the inventory gives the origin `measured` or `named-entry` was set to x0.75 and x1.25 "
        "of its value, one at a time, and the random-100 campaign keys were measured on the third random sample (none "
        "of these constants was tuned on it), default view only, against an unperturbed baseline run the same way. "
        "The perturbation reaches every reader (scripts/perturb.py: the module is rewritten as it is imported, so "
        "`from ... import` binders, defaults and constants derived from it at import follow; `perturb.py reach` "
        "checked each run). Keys: each entry's pass or fail, the failing count per A- and S-criterion, the named "
        "fallbacks and the figure area of the default view; time is recorded, not judged, and so is A2 (the default "
        "view within 60 s of wall time: it moved between two identical baselines with the load of the machine), so an "
        "entry passes here when no other A-criterion fails. **stable**: no key moved; "
        f"**sensitive**: some key moved; **fragile**: the pass count moved by {doc['rules']['fragile_pass']} or more, "
        "or a criterion failing in no baseline entry failed in some entry. A constant takes the worse of its two "
        "perturbations, from stage B (the whole sample) when it ran there, else from stage A (its subset only).",
        "",
        "The baseline draws the default view only, so its pass count is not the third sample's published one (every "
        "view drawn, `docs/validation/random100-third/report.md`): criteria needing other views are not judged here.",
        "",
        "## Counts",
        "",
        "| class | constants |",
        "|---|---|",
        *[f"| {c} | {n} |" for c, n in doc["counts"].items() if c != "incomplete" or n],
        "",
        f"{len(doc['constants'])} constants perturbed, {len(doc['not_perturbable'])} not perturbable, "
        f"{len(doc['no_perturbation'])} with both perturbations skipped; zl time "
        f"{hours:.1f} h (stage A {a['seconds'] / 3600:.1f} h"
        + (f", stage B {b['seconds'] / 3600:.1f} h" if b else ", stage B not run")
        + ").",
        "",
        f"Of the stable constants, {len(doc['stable_unread'])} were read by no line that ran in the default view "
        "(the reads run: a line inside a function naming the constant, perturb.install_trace), so their stability "
        "here is no evidence either way: " + (", ".join(f"`{i}`" for i in doc["stable_unread"]) or "none") + ".",
        "",
        "## Limitations: the fragile constants",
        "",
    ]
    fragile = [r for r in doc["constants"] if r["class"] == "fragile"]
    if not fragile:
        lines.append("None.")
    for r in fragile:
        said = "; ".join(f"{f} ({r['runs'][f]['new']}): {_said(r['runs'][f].get(r['stage']))}" for f in r["runs"])
        lines.append(f"- `{r['id']}` = {r['value']} ({r['origin']}), stage {r['stage']}: {said}")
    sub = doc["subset"]
    lines += [
        "",
        "## Stage A: the subset",
        "",
        f"{len(sub['ids'])} entries of the third sample (seed {sub['seed']}, target {sub['size']}; per stratum in "
        "proportion to its eligible entries, largest remainder, a pair one unit): " + ", ".join(sub["ids"]) + ".",
        "",
        f"Left out, recorded slower than {doc['rules']['slow_s']} s in the third run (every view): "
        + ", ".join(f"{i} ({s:.0f} s)" for i, s in sorted(sub["left_out"].items(), key=lambda x: -x[1]))
        + ".",
        "",
    ]
    for name, st in (("A", a), ("B", b)):
        if not st:
            continue
        det = st["determinism"]
        lines += [
            f"Stage {name}: {len(st['entries'])} entries, baseline {st['baseline_passed']} pass; failing per criterion "
            + (", ".join(f"{k} {n}" for k, n in st["baseline_failing"].items()) or "none")
            + f". The two baselines are {'identical (keys and SVGs)' if det['identical'] else 'NOT identical'}"
            + ("" if det["identical"] else f": {det['moved'][:10]} SVG {det['svg'][:10]}")
            + f". A2 (time, not judged) failed in {', '.join(st['baselines_timed_failing'][0]) or 'none'} and in "
            + f"{', '.join(st['baselines_timed_failing'][1]) or 'none'} of the two baselines"
            + f". {len(st['runs'])} perturbation runs; {st['seconds'] / 3600:.2f} h of zl time."
            + (
                ""
                if "reads_run_moved" not in st
                else " The reads run (unperturbed, traced) "
                + ("changed no key." if not st["reads_run_moved"] else f"MOVED: {st['reads_run_moved'][:10]}.")
            ),
            "",
        ]
    lines += [
        "## Per constant",
        "",
        "read in: the entries of the deciding stage a line naming the constant ran in (default view), + when the "
        "judging read it.",
        "",
        "| constant | origin | value | read in | x0.75 | x1.25 | class |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in doc["constants"]:
        cells = []
        for f in ("x0.75", "x1.25"):
            run = r["runs"].get(f)
            if run is None:
                cells.append("skipped")
                continue
            said = f"{run['new']}: A {_said(run.get('A'))}"
            if run.get("B"):
                said += f"; B {_said(run['B'])}"
            cells.append(said.replace("|", "\\|"))
        st = a if r["stage"] == "A" else b
        read = (r.get("read_in") or {}).get(r["stage"])
        read = "-" if read is None else f"{read}/{len(st['entries'])}" + ("+" if r["read_judging"][r["stage"]] else "")
        lines.append(
            f"| `{r['id']}` | {r['origin']} | {r['value']} | {read} | {cells[0]} | {cells[1]} | "
            f"{r['class'] or 'incomplete'} ({r['stage']}) |"
        )
    lines += ["", "## Skipped perturbations", ""]
    lines += [f"- `{s['id']}` = {s['value']} x{s['factor']}: {s['why']}" for s in doc["skipped"]] or ["None."]
    lines += ["", "## Not perturbable", ""]
    lines += [f"- `{s['id']}` = {s['value']}: {s['why']}" for s in doc["not_perturbable"]] or ["None."]
    return "\n".join(lines) + "\n"


def cmd_report(args) -> int:
    doc = build(Path(args.a), Path(args.b) if args.b else None, SOURCE_STAMP)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    (out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(f"{doc['counts']}; fragile: {doc['limitations']}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan", help="the constants, their perturbations and stage A's subset")
    p.add_argument("--size", type=int, default=SUBSET)
    p.set_defaults(func=cmd_plan)
    r = sub.add_parser("run", help="run a stage's runs not run yet")
    r.add_argument("--stage", choices=("A", "B"), required=True)
    r.add_argument("--out", required=True, help="the stage's folder (never inside docs/validation/random100-third)")
    r.add_argument("--labels", help="baselines, perturbations or a comma-separated list of run names")
    r.add_argument("--size", type=int, default=SUBSET)
    r.add_argument("--jobs", type=int, default=8, help="entries drawn at once in each run")
    r.add_argument("--parallel", type=int, default=1, help="runs at once")
    r.set_defaults(func=cmd_run)
    w = sub.add_parser("report", help="judge the stages and write report.md and report.json")
    w.add_argument("--a", required=True)
    w.add_argument("--b")
    w.add_argument("--out", default=str(OUT))
    w.set_defaults(func=cmd_report)
    args = ap.parse_args(argv)
    if os.environ.get(ENV, "").strip():
        ap.error(f"{ENV} is set: this process must run unperturbed (it sets {ENV} for each run's child)")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
