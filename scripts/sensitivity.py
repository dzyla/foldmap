# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""How fragile foldmap's tuned constants are (plan 28 Task 9.2, rulings Q29, Q30): each module constant the inventory
(docs/validation/inventory/inventory.md) gives the origin "measured" or "named-entry" is perturbed to x0.75 and x1.25,
one at a time, and the random-100 campaign keys are measured on the third random sample. That sample is not the
tuning set (the first), but it has been run and looked at since plan 15: two of its entries are named by constants
measured here (9FZL for route._MAX_NODES_WINDOWS, 6UW4 for membrane.JOIN_GAP; the
report computes this from the origins that name entries, a lower bound), and the membrane constants were tuned in
plan 26 with 6DQN, 6BO8, 6UW4 and 6QUE among the entries looked at. Never the fourth sample, whose protocol any look
voids.

    python scripts/sensitivity.py plan                              # the constants, perturbations and stage A subset
    python scripts/sensitivity.py run --stage A --out DIR [--labels baselines|perturbations|L1,L2] [--jobs 8]
                                      [--parallel 4]
    python scripts/sensitivity.py run --stage B --out DIR           # the whole sample, every perturbation (Q30)
    python scripts/sensitivity.py run --stage X --out DIR           # which view reads each constant (one traced run)
    python scripts/sensitivity.py run --stage V --view membrane --survey X --default-stage B --out DIR
    python scripts/sensitivity.py report --a A --b B [--survey X] [--view V ...] [--plan-from REPORT]
                                                                    # docs/validation/sensitivity/

A perturbation is the value times the factor, an int rounded half away from zero; it is skipped and listed when that
gives the same value or one that cannot be meant (an int count of 1 or more falling under 1). A tuple, a bool or any
other non-number is listed as not perturbable. Each run is `random100.py run --not-for-paper --views V` (the default
view; stage X every view but the app's, stage V the default and its own) in a
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
constant, and in which view; it must move no key. Stage X: one such traced run of every view but the app's (its job
process is not traced) on the whole sample.
Stage V (one per view): the constants stable in stage B that no line of the default view read and that the view
(VIEW_STAGES) reads in the most entries in stage X, on the entries it reads them in, the default view drawn first
(focus and domains follow from it) and the keys taken from the view; two baselines, a reads run and their
perturbations. A constant no stage V reads is "not exercised" (ruling Q30), never stable.

Keys, per run: each entry's pass or fail; the failing count per A- and S-criterion; each entry's named fallbacks
(the keyed view); each entry's figure area; time is recorded, not judged, and so are the time criteria (TIMED: A2,
the default view's wall time within 60 s, moved with the load of the machine between two identical baselines; an
entry passes here when no other A-criterion fails). A perturbation is stable when no key moves from the baseline,
sensitive when any does, fragile when the pass count moves by FRAGILE_PASS or more or a criterion that fails in no
baseline entry fails in some entry; a constant takes the worse of its two perturbations, from stage B where it ran
there (its stage V for one the default view does not read). Fragile constants are the limitations list for the
paper.

The report says which code changed after its runs (CODE_SINCE, with the prune record's pruned constants and default
baseline: plan 28 final review F1); its classes describe the runs' code until stages B, V and X are rerun at the
release (docs/validation/freeze-checklist.md). `report --plan-from` re-renders it from the kept runs with the plan as
the record holds it (constants the code no longer has, pruned since, included), so nothing is run again."""

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
PRUNE = ROOT / "docs" / "validation" / "prune" / "report.json"
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
CLASSES = ("stable", "sensitive", "fragile")  # by worsening
NOT_EXERCISED = "not exercised"  # ruling Q30: stable, and no view measured here read it: no evidence
# the entries plan 26 tuned the membrane constants on (docs/superpowers/plans/2026-10-24-26-membrane-panel.md, Task 2)
PLAN26_MEMBRANE = ("6DQN", "6BO8", "6UW4", "6QUE")
# Plan 28 final review F1: the src commits after the runs' code (the src of MEASURED_AT, its fingerprint in the report)
# up to the record's re-render, each with what it changed; the runs are not made again until the freeze
MEASURED_AT = "5f09684"
CODE_SINCE = (
    ("ee5f34d", "deep checks (foldmap check --deep, --json): said, never drawn"),
    ("ddc408c", "deep checks' fix round"),
    ("24676d3", "NA fix F7 (nacheck: stems_nearer_frac only where an element can be nearer)"),
    ("36adc3b", "NA fix F2 (layout, frame, cli, nacheck: each nucleic acid unrolled with the protomer holding it)"),
    ("3380435", "NA fix F4 (na_ends, route_na: a gap meets a free end from the side it comes from)"),
    ("218f4bb", "NA fix F6 (na_labels: where runs meet, the other run slides)"),
    ("4aec81a", "NA fix round 1 (the unrolled judgement)"),
    ("e060370", "NA fix round 1 (F4's straight chord limited to its reach)"),
    ("17ddbbc", "NA fix round 1 (F6's slide loop)"),
    ("eac6756", "the readability floor: said, never drawn"),
    ("e6fcf36", "the readability floor's fix round"),
    ("ba724ad", "the domain split confidence: said, never drawn"),
    ("42e3fe8", "the domain split confidence's fix round"),
    ("c64ea7a", "the prune: two special cases measured here removed"),
    ("b730a05", "the prune's fix round (comments and docstrings)"),
    ("702b4de", "the final review: --strict-floor fails closed, check's help text"),
)


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
    """The run name, file and folder name of a perturbation: "module:NAME" with its ":" (a name Windows forbids in a
    path) made ".", then x and the factor ("layout._STACK_PADx0.75"). Every path named after a constant's run is
    this."""
    return safe_name(f"{cid.replace(':', '.')}x{factor}")


def safe_name(name: str) -> str:
    """`name`, refused (ValueError) unless it is a file name every OS takes (no : \\ / * ? \" < > |)."""
    if not name or any(c in name for c in ':\\/*?"<>|'):
        raise ValueError(f'{name!r} is no file name on every OS (Windows forbids : \\ / * ? " < > |)')
    return name


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


def units_entries(sample: dict) -> list[dict]:
    """Every entry of a sample (its refusals included)."""
    return [e for u in units(sample) for e in u]


def named_entries(origin: str) -> set[str]:
    """The entries an inventory origin names ("named-entry:9FZL,7NZ2" -> {"9FZL", "7NZ2"}); none for any other."""
    kind, _, ids = origin.partition(":")
    return {i.strip() for i in ids.split(",") if i.strip()} if kind.strip() == "named-entry" else set()


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


def keys_of(report: dict, view: str = "default") -> dict[str, dict]:
    """The campaign keys of every entry of a random-100 report: pass, the failing A- and S-criteria, `view`'s named
    fallbacks and figure area, and (recorded, not judged) its seconds, its failing time criteria (TIMED: a wall-clock
    limit moves with the machine's load, not the code, with their details) and SVG hash. pass is no A-criterion but a
    timed one failing."""
    out = {}
    for e in report["entries"]:
        crit, sci = e.get("criteria") or {}, e.get("science") or {}
        d = ((e.get("views") or {}).get(view) or {}).get("result") or {}
        judged = [k for k, c in crit.items() if c.get("pass") is False and k not in TIMED]
        timed = {k: c.get("detail") for k, c in crit.items() if c.get("pass") is False and k in TIMED}
        out[e["id"]] = {
            "pass": not judged,
            "failing": sorted(judged + [k for k, c in sci.items() if c.get("pass") is False]),
            "fallbacks": dict(sorted((d.get("fallbacks") or {}).items())),
            "area": d.get("area"),
            "seconds": (e.get("process") or {}).get("seconds"),
            "timed_failing": sorted(timed),
            "timed_detail": dict(sorted(timed.items())),
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


def worst(classes) -> str | None:
    """The worst class of a constant's perturbations; None (incomplete) when one has no class (crashed, unreached)."""
    classes = list(classes)
    return None if None in classes else max(classes, key=CLASSES.index, default="stable")


# ------------------------------------------------------------------------------------------------- running


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _view_of(argv: list[str]) -> str:
    """The view a random100.py view process draws (`view --name N --work W [--na] VIEW PATH -- ...`)."""
    k = argv.index("--work") + 2
    return argv[k + 1] if argv[k] == "--na" else argv[k]


def read_traces(folder: Path, report: dict, by_view: bool = False) -> dict:
    """Which entries read each constant, from the processes' records (perturb.install_trace): a view's or an entry's
    process counts for its entry, the run's own (the judging) as "(judging)". -> {constant: [entries]}; with
    `by_view`, {constant: {view: [entries]}}, an entry's own process (its judging, mkdssp, CATH) as view "(entry)"
    and the run's as "(judging)"."""
    ids = {e.get("name") or e["id"]: e["id"] for e in report["entries"]}
    out: dict[str, dict[str, set[str]]] = {}
    for p in sorted(folder.glob("*.json")):
        if p.name == "names.json":
            continue
        doc = json.loads(p.read_text(encoding="utf-8"))
        argv = doc["argv"]
        if len(argv) > 3 and argv[1] == "view" and argv[2] == "--name":
            who, view = ids.get(argv[3], argv[3]), _view_of(argv)
        elif len(argv) > 2 and argv[1] == "entry":
            who, view = argv[2], "(entry)"
        else:
            who = view = "(judging)"
        for c in doc["reads"]:
            out.setdefault(c, {}).setdefault(view, set()).add(who)
    if by_view:
        return {c: {v: sorted(w) for v, w in sorted(vs.items())} for c, vs in sorted(out.items())}
    return {c: sorted(set().union(*vs.values())) for c, vs in sorted(out.items())}


def trace_records(folder: Path, report: dict) -> dict:
    """How many view-process records each entry of a traced run left (perturb.install_trace writes one at a process's
    exit; a view killed at its timeout writes none), and the entries with none: their reads are unknown, not none."""
    ids = {e.get("name") or e["id"]: e["id"] for e in report["entries"]}
    count: Counter = Counter()
    for p in sorted(folder.glob("*.json")):
        if p.name == "names.json":
            continue
        argv = json.loads(p.read_text(encoding="utf-8"))["argv"]
        if len(argv) > 3 and argv[1] == "view" and argv[2] == "--name":
            count[ids.get(argv[3], argv[3])] += 1
    every = [e["id"] for e in report["entries"]]
    return {"per_entry": {i: count.get(i, 0) for i in every}, "no_view_record": [i for i in every if not count[i]]}


def run_one(
    name: str,
    perturb: str | None,
    ids: list[str] | None,
    out: Path,
    jobs: int,
    sample: Path,
    trace: list[str] | None = None,
    views: str = "default",
    keyed: str = "default",
) -> dict:
    """One random-100 run (`views`, default the default view; --not-for-paper) with `perturb` set, then its keys
    (view `keyed`) and reach; with `trace` (constant ids), unperturbed, which entries read each, and in which view
    (perturb.install_trace). The run's figures and sheets are deleted, its report kept."""
    where = out / "runs" / safe_name(name)  # a run name is a path: never a "module:NAME" (label)
    if where.exists():
        shutil.rmtree(where)  # a run stopped half way: run it again whole
    where.mkdir(parents=True)
    cmd = [sys.executable, "-I", str(SCRIPTS / "random100.py"), "run", "--not-for-paper", "--views", views]
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
        "views": views,
        "entries": keys_of(rep, keyed) if rep else None,
        "reads": read_traces(where / "reads", rep) if trace and rep else None,
        "reads_by_view": read_traces(where / "reads", rep, by_view=True) if trace and rep else None,
        "trace_records": trace_records(where / "reads", rep) if trace and rep else None,
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
    keyed = stage_info(where).get("keyed", "default")
    out = {}
    for p in sorted((where / "keys").glob("*.json")):
        doc = json.loads(p.read_text(encoding="utf-8"))
        rep = where / "runs" / p.stem / "report.json"
        if rep.is_file():
            report = _script("random100").read_report(rep)
            doc["entries"] = keys_of(report, keyed)
            traced = where / "runs" / p.stem / "reads"
            if traced.is_dir():  # the process records, read again by this code
                doc["reads"] = read_traces(traced, report)
                doc["reads_by_view"] = read_traces(traced, report, by_view=True)
                doc["trace_records"] = trace_records(traced, report)
        out[p.stem] = doc
    return out


def stage_info(where: Path) -> dict:
    """A stage folder's stage.json (its stage, views, keyed view, entries and constants); {} for stages A and B run
    before it was written (the default view, keyed there)."""
    p = where / "stage.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {}


# the views a view stage may measure a constant in: a drawn figure with its own fallbacks and area (criteria's
# FIGURE_VIEWS less the default); seed1 and roundtrip redraw the default figure, and the app's job process is not traced
VIEW_STAGES = ("deposited", "beta", "focus", "domains", "membrane")
SURVEY = "reads-all"  # stage X's one run: every view, unperturbed, traced


def unread_in_default(default_reads: dict[str, list[str]], the_plan: dict) -> list[str]:
    """The perturbable constants no process of the whole-sample default-view reads run read (stage B's)."""
    return sorted({p["id"] for p in the_plan["perturbations"]} - {c for c, who in default_reads.items() if who})


def view_targets(the_plan: dict, default_reads: dict, survey: dict) -> dict[str, list[str]]:
    """Each constant unread in the default view, given to the view of VIEW_STAGES that read it in the most entries
    in the survey (stage X; ties in VIEW_STAGES order) -> {view: [constants]}; a constant no such view read is in
    none (not exercised)."""
    out: dict[str, list[str]] = {}
    for c in unread_in_default(default_reads, the_plan):
        seen = {v: len(survey.get(c, {}).get(v, [])) for v in VIEW_STAGES}
        best = max(VIEW_STAGES, key=lambda v: (seen[v], -VIEW_STAGES.index(v)))
        if seen[best]:
            out.setdefault(best, []).append(c)
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
    ids, views, keyed = None, "default", "default"
    every = sorted({p["id"] for p in the_plan["perturbations"]})
    traced = every  # the constants the reads run traces
    if args.stage == "A":  # the subset, fixed in the stage's folder by its first run
        sub = subset(sample, recorded_seconds(rep), args.size)
        held = out / "subset.json"
        if held.is_file() and json.loads(held.read_text(encoding="utf-8"))["ids"] != sub["ids"]:
            raise SystemExit(f"{held} holds another subset (another --size?)")
        held.write_text(json.dumps({**sub, "source_stamp": SOURCE_STAMP}, indent=1) + "\n", encoding="utf-8")
        ids = sub["ids"]
    if args.stage in ("A", "B"):
        runs = stage_runs(args.stage, the_plan)
    elif args.stage == "X":  # the survey: which view reads each constant
        # every view but the app's, whose job process the trace never reaches (random100 refuses it traced)
        runs, views = [(SURVEY, None)], ",".join(v for v in _script("random100").VIEWS if v != "app")
    else:  # V: the constants unread in the default view, in the view that reads them, on the entries it reads them in
        if not (args.view and args.survey and args.default_stage):
            raise SystemExit("stage V needs --view, --survey (stage X's folder) and --default-stage (stage B's)")
        survey = read_keys(Path(args.survey))[SURVEY]["reads_by_view"]
        default_reads = read_keys(Path(args.default_stage))[READS]["reads"]
        traced = view_targets(the_plan, default_reads, survey).get(args.view, [])
        if not traced:
            raise SystemExit(f"no constant unread in the default view is read in view {args.view}")
        ids = sorted({w for c in traced for w in survey[c].get(args.view, [])})
        views, keyed = f"default,{args.view}", args.view  # the default view first: focus and domains follow from it
        runs = [(b, None) for b in (*BASELINES, READS)]
        runs += [(p["label"], p["perturb"]) for p in the_plan["perturbations"] if p["id"] in traced]
    info = {"stage": args.stage, "views": views, "keyed": keyed, "ids": ids, "constants": traced}
    held = out / "stage.json"
    if held.is_file() and {k: v for k, v in json.loads(held.read_text(encoding="utf-8")).items() if k in info} != info:
        raise SystemExit(f"{held} holds another stage's settings")
    held.write_text(json.dumps({**info, "source_stamp": SOURCE_STAMP}, indent=1) + "\n", encoding="utf-8")
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
        trace = traced if r[0] in (READS, SURVEY) else None
        doc = run_one(r[0], r[1], ids, out, args.jobs, THIRD / "sample.json", trace, views, keyed)
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
    info = stage_info(where)
    base = keys[BASELINES[0]]["entries"]
    other = keys[BASELINES[1]]["entries"]
    svg = sorted(i for i in base if base[i]["svg_sha256"] != (other.get(i) or {}).get("svg_sha256"))
    timing = [json.loads(line) for line in (where / "stages.jsonl").read_text(encoding="utf-8").splitlines() if line]
    unperturbed = [r for r in (*BASELINES, READS) if keys.get(r) and keys[r]["entries"]]
    out = {
        "views": info.get("views", "default"),
        "keyed": info.get("keyed", "default"),
        "constants": info.get("constants"),
        "entries": sorted(base),
        "baseline_passed": passed(base),
        "baseline_failing": dict(sorted(failing_counts(base).items())),
        "baselines_timed_failing": [  # recorded, not judged (TIMED)
            sorted(i for i, e in run.items() if e.get("timed_failing")) for run in (base, other)
        ],
        # random100's own pass, a TIMED criterion judged too, for each unperturbed run: what A2 does to a count
        "passed_with_time": {
            r: sum(e["pass"] and not e.get("timed_failing") for e in keys[r]["entries"].values()) for r in unperturbed
        },
        "timed_detail": {
            r: {i: e["timed_detail"] for i, e in keys[r]["entries"].items() if e.get("timed_detail")}
            for r in unperturbed
        },
        "determinism": {"identical": not moved(base, other) and not svg, "moved": moved(base, other), "svg": svg},
        "seconds": round(sum(t["seconds"] for t in timing), 1),
        "run_seconds": round(sum(k["seconds"] for k in keys.values()), 1),
        "runs": {},
        # which code made each run's random-100 report (its stamp): build() refuses runs of more than one foldmap
        "codes": {name: _code(k) for name, k in keys.items()},
    }
    traced = keys.get(READS)
    if traced is not None:  # the reads run draws unperturbed: it must change no key either
        out["reads"] = traced["reads"] or {}
        out["reads_by_view"] = traced.get("reads_by_view") or {}
        out["reads_run_moved"] = moved(base, traced["entries"]) if traced["entries"] else ["no report"]
    for p in the_plan["perturbations"]:
        k = keys.get(p["label"])
        if k is None:
            continue
        if k["entries"] is None:  # the run crashed: nothing measured, so no class (incomplete)
            got = {"class": None, "error": f"no report (exit {k['exit']})", "moved": []}
        else:
            got = judge(base, k["entries"])
        got["reached"] = bool((k.get("reach") or {}).get("reached"))
        if not got["reached"]:  # the value did not reach its readers: what moved is not its doing (incomplete)
            got["judged"], got["class"] = got["class"], None
        got["seconds"] = k["seconds"]
        out["runs"][p["label"]] = got
    if traced is not None:
        out["trace_records"] = traced.get("trace_records")
    return out


def _code(keys_doc: dict) -> dict:
    """The code a run's random-100 report was made by: its foldmap code fingerprint and the commit it names."""
    stamp = keys_doc.get("report_stamp") or {}
    return {"code_fingerprint": stamp.get("code_fingerprint"), "commit": (stamp.get("git") or {}).get("commit")}


def _read(stage: dict, cid: str, view: str | None = None) -> tuple[int, bool]:
    """(the entries of `stage` a line naming the constant ran in, whether the judging read it); with `view`, in that
    view's processes only."""
    if view is None:
        who = stage.get("reads", {}).get(cid, [])
    else:
        who = stage.get("reads_by_view", {}).get(cid, {}).get(view, [])
    return len([w for w in who if w != "(judging)"]), "(judging)" in who


def settle(row: dict, stage_v: dict[str, dict]) -> None:
    """A constant's class once the view stages are in (ruling Q30): one that moved, or that a line of its deciding
    stage read, keeps its class; one stable and unread there takes the class of the view stage measuring it, and is
    NOT_EXERCISED when none does or that stage's reads run saw no read of it in its view either."""
    if row["class"] != "stable" or row["read_in"] or row["read_judging"]:
        return
    view = next((v for v, s in stage_v.items() if row["id"] in (s["constants"] or [])), None)
    if view is None:
        row["class"] = NOT_EXERCISED
        return
    judged = [r[view]["class"] for r in row["runs"].values() if r.get(view)]
    row["stage"], (row["read_in"], row["read_judging"]) = view, _read(stage_v[view], row["id"], view)
    row["class"] = worst(judged) if len(judged) == len(row["runs"]) else None
    if row["class"] == "stable" and not row["read_in"]:
        row["class"] = NOT_EXERCISED


def build(
    a: Path,
    b: Path | None,
    stamp: dict,
    views: list[Path] = (),
    survey: Path | None = None,
    the_plan: dict | None = None,
    sample: dict | None = None,
    prune: Path | None = PRUNE,
) -> dict:
    """The report of stage A (folder `a`), B (`b`, or None before it runs) and the view stages (`views`), with stage
    X's survey (`survey`) of the views that read each constant, stamped with `stamp`. A constant's class is from
    stage B (else A); a constant stable there that no line of the default view read takes its view stage's class,
    and is "not exercised" when no view stage covers it or its view stage's reads run saw no read either. the_plan and
    sample (sample.json) default to this checkout's (tests pass their own)."""
    the_plan = plan() if the_plan is None else the_plan
    if sample is None:
        sample = json.loads((THIRD / "sample.json").read_text(encoding="utf-8"))
    in_sample = {e["id"] for e in units_entries(sample)}
    stage_a = _stage(a, the_plan)
    stage_b = _stage(b, the_plan) if b else None
    stage_v = {}
    for where in views:
        st = _stage(where, the_plan)
        stage_v[st["keyed"]] = st
    survey_keys = read_keys(survey) if survey else {}
    surveyed = survey_keys[SURVEY]["reads_by_view"] if survey else None
    survey_seconds = (
        round(
            sum(
                json.loads(t)["seconds"]
                for t in (survey / "stages.jsonl").read_text(encoding="utf-8").splitlines()
                if t
            ),
            1,
        )
        if survey and (survey / "stages.jsonl").is_file()
        else None
    )
    # one foldmap made every run: the runs are re-keyed by this code, so their reports must measure the same foldmap
    codes = {f"{n}/{r}": c for n, st in (("A", stage_a), ("B", stage_b), *stage_v.items()) if st for r, c in
             st["codes"].items()}  # fmt: skip
    codes |= {f"X/{r}": _code(k) for r, k in survey_keys.items()}
    fingerprints = {c["code_fingerprint"] for c in codes.values()}
    if len(fingerprints) != 1 or None in fingerprints:
        raise SystemExit(f"refused: the runs' reports name {len(fingerprints)} foldmap code fingerprints "
                         f"({sorted(map(str, fingerprints))}); a report must rest on one foldmap")  # fmt: skip
    rows = {}
    for p in the_plan["perturbations"]:
        row = rows.setdefault(p["id"], {"id": p["id"], "origin": p["origin"], "value": p["value"], "runs": {}})
        row["named_in_sample"] = sorted(named_entries(p["origin"]) & in_sample)
        got = {"new": p["new"], "A": stage_a["runs"].get(p["label"])}
        if stage_b and p["label"] in stage_b["runs"]:
            got["B"] = stage_b["runs"][p["label"]]
        for view, st in stage_v.items():
            if p["label"] in st["runs"]:
                got[view] = st["runs"][p["label"]]
        row["runs"][f"x{p['factor']}"] = got
    for row in rows.values():
        stage = "B" if all("B" in r for r in row["runs"].values()) else "A"
        judged = [r[stage]["class"] for r in row["runs"].values() if r.get(stage)]
        row["stage"] = stage
        row["class"] = worst(judged) if len(judged) == len(row["runs"]) else None  # None: a run missing
        st = stage_b if stage == "B" else stage_a
        row["read_in"], row["read_judging"] = _read(st, row["id"]) if "reads" in st else (None, False)
        if surveyed is not None:
            row["survey"] = {v: len(w) for v, w in surveyed.get(row["id"], {}).items()}
        settle(row, stage_v)
    names = (*reversed(CLASSES), NOT_EXERCISED)
    counts = Counter(r["class"] for r in rows.values())
    return {
        "what": "foldmap's measured and named-entry constants perturbed x0.75 and x1.25 on the third random sample "
        "(scripts/sensitivity.py, plan 28 Task 9.2, rulings Q29, Q30)",
        "date": datetime.date.today().isoformat(),
        "source_stamp": stamp,
        "sample": "docs/validation/random100-third/sample.json",
        "sample_sha256": _sha256(THIRD / "sample.json") if (THIRD / "sample.json").is_file() else None,
        "inventory_sha256": _sha256(INVENTORY) if INVENTORY.is_file() else None,
        "factors": list(FACTORS),
        "rules": {"slow_s": SLOW_S, "fragile_pass": FRAGILE_PASS, "subset": SUBSET, "seed": SEED, "timed": TIMED},
        "subset": {
            k: v for k, v in json.loads((a / "subset.json").read_text(encoding="utf-8")).items() if k != "source_stamp"
        },
        # the sample's entries a constant was set on (its inventory origin): the third sample is not the tuning set,
        # but it has been run and looked at since plan 15, and a few constants were set on its entries since
        "named_in_sample": {
            c["id"]: hits
            for c in [*the_plan["perturbations"], *the_plan["skipped"], *the_plan["not_perturbable"]]
            if (hits := sorted(named_entries(c["origin"]) & in_sample))
        },
        "stage_a": stage_a,
        "stage_b": stage_b,
        "survey_seconds": survey_seconds,
        "code_fingerprint": fingerprints.pop(),  # every run's report, checked above
        "run_commits": sorted({c["commit"] for c in codes.values() if c["commit"]}),
        "view_stages": stage_v,
        "counts": {c: counts.get(c, 0) for c in names} | {"incomplete": counts.get(None, 0)},
        "constants": list(rows.values()),
        "limitations": [r["id"] for r in rows.values() if r["class"] == "fragile"],
        "not_exercised": [r["id"] for r in rows.values() if r["class"] == NOT_EXERCISED],
        "no_perturbation": sorted({s["id"] for s in the_plan["skipped"]} - set(rows)),
        "skipped": [{k: s[k] for k in ("id", "value", "factor", "why")} for s in the_plan["skipped"]],
        "not_perturbable": [{k: s[k] for k in ("id", "value", "why")} for s in the_plan["not_perturbable"]],
        "code_since": code_since(rows, stage_b, prune),
    }


def code_since(rows: dict, stage_b: dict | None, prune: Path | None = PRUNE) -> dict:
    """What changed in the code after the runs (plan 28 final review F1): CODE_SINCE, the measured constants the prune
    record removed since, and its default baseline beside stage B's (the same sample, view and convention: A2 not
    judged)."""
    out: dict = {"measured_at": MEASURED_AT, "commits": [{"commit": c, "what": w} for c, w in CODE_SINCE]}
    if prune is not None and prune.is_file():
        rec = json.loads(prune.read_text(encoding="utf-8"))
        out["pruned"] = [c for c in rec.get("pruned", []) if c in rows]
        out["prune_commit"] = ((rec.get("source_stamp") or {}).get("git") or {}).get("commit")
        out["prune_baseline"] = (rec.get("baselines") or {}).get("r100_passed")
        out["prune_entries"] = (rec.get("baselines") or {}).get("r100_entries")
    if stage_b:
        out["baseline"], out["entries"] = stage_b["baseline_passed"], len(stage_b["entries"])
    return out


def plan_of_record(old: dict, items: list[dict] | None = None) -> dict:
    """The plan as a committed report holds it (its constants' perturbations, the skipped and the not perturbable),
    to re-render it from the kept runs when the code no longer has every constant it measured (a pruned one): the
    values and perturbed values are the record's, never the code's now. Origins are the record's rows', else the
    inventory's (`items`, constants() by default) for a skipped or not perturbable one; one in neither stops it."""
    origins = {c["id"]: c["origin"] for c in (constants() if items is None else items)}
    origins |= {r["id"]: r["origin"] for r in old["constants"]}
    runs = []
    for r in old["constants"]:
        module, name = r["id"].split(":", 1)
        for key, run in r["runs"].items():
            f = float(key[1:])
            runs.append(
                {
                    "id": r["id"],
                    "module": module,
                    "name": name,
                    "origin": r["origin"],
                    "value": r["value"],
                    "factor": f,
                    "new": run["new"],
                    "label": label(r["id"], f),
                    "perturb": f"{module}:{name}={run['new']!r}",
                }  # fmt: skip
            )
    missing = sorted({s["id"] for s in [*old["skipped"], *old["not_perturbable"]]} - set(origins))
    if missing:
        raise SystemExit(f"no origin for {missing}: neither the record nor the inventory has them")
    return {
        "constants": len({*(r["id"] for r in old["constants"]), *old["no_perturbation"]}) + len(old["not_perturbable"]),
        "perturbations": runs,
        "skipped": [{**s, "origin": origins[s["id"]]} for s in old["skipped"]],
        "not_perturbable": [{**s, "origin": origins[s["id"]]} for s in old["not_perturbable"]],
    }


def _said(run: dict | None) -> str:
    if not run:
        return "-"
    if run.get("error"):
        return f"**incomplete**: {run['error']}"
    shown = run["class"] or f"**incomplete** (would be {run.get('judged')})"
    parts = [shown if shown != "fragile" else "**fragile**"]
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


def _a2(st: dict) -> str:
    """What the time criterion A2 did to the unperturbed runs' pass counts (random100 judges it; this report not)."""
    counts = ", ".join(f"{r} {n}" for r, n in st["passed_with_time"].items())
    details = "; ".join(
        f"{r}: " + ", ".join(f"{i} {d.get('A2', d)}" for i, d in sorted(got.items()))
        for r, got in st["timed_detail"].items()
        if got
    )
    return (
        f"With A2 judged, as random100 judges it, the unperturbed runs of this one code pass {counts} of "
        f"{len(st['entries'])} ({details or 'A2 failed in none'}), their keys and SVGs otherwise identical."
    )


def _stage_line(name: str, st: dict) -> str:
    det = st["determinism"]
    return (
        f"Stage {name} ({st['views']}, keyed on the {st['keyed']} view): {len(st['entries'])} entries, baseline "
        f"{st['baseline_passed']} pass; failing per criterion "
        + (", ".join(f"{k} {n}" for k, n in st["baseline_failing"].items()) or "none")
        + f". The two baselines are {'identical (keys and SVGs)' if det['identical'] else 'NOT identical'}"
        + ("" if det["identical"] else f": {det['moved'][:10]} SVG {det['svg'][:10]}")
        + f". {len(st['runs'])} perturbation runs; {st['seconds'] / 3600:.2f} h of zl time."
        + (
            ""
            if "reads_run_moved" not in st
            else " The reads run (unperturbed, traced) "
            + ("changed no key." if not st["reads_run_moved"] else f"MOVED: {st['reads_run_moved'][:10]}.")
        )
        + " "
        + _a2(st)
    )


def _since_lines(doc: dict) -> list[str]:
    """The "code since this run" section (plan 28 final review F1): [] for a report without code_since."""
    since = doc.get("code_since")
    if not since:
        return []
    by = ", ".join(f"{c['commit']} ({c['what']})" for c in since["commits"])
    pruned = since.get("pruned") or []
    lines = [
        "## Code since this run",
        "",
        f"These runs measured foldmap code {doc['code_fingerprint']}, the src of {since['measured_at']}. The record "
        "was re-rendered from the kept runs after the code changed, and nothing was run again. Commits to src since: "
        f"{by}.",
        "",
    ]
    if pruned:
        lines += [
            f"Measured here, pruned since: {', '.join(f'`{c}`' for c in pruned)} (removed from the code by the prune "
            f"record, docs/validation/prune/report.md, at {str(since.get('prune_commit'))[:7]}); their rows below say "
            "so.",
            "",
        ]
    fragile = [r["id"] for r in doc["constants"] if r["class"] == "fragile"]
    if "layout:_SIDE_CLEAR" in fragile:
        lines += [
            "The fragile constant `layout:_SIDE_CLEAR` places nucleic-acid stems, and what moved under it is mostly "
            "nucleic-acid fallbacks (na_gap_unclear, na_tail_straight, na_run_unclear): the NA fixes F4 and F6 above "
            "changed exactly those paths, so its class is this code's, not the release's.",
            "",
        ]
    if since.get("prune_baseline") is not None and since.get("baseline") is not None:
        lines += [
            f"The default baseline moved: stage B here passes {since['baseline']} of {since['entries']} (A2 not "
            f"judged); the prune record's baseline, the same sample, view and convention on later code, passes "
            f"{since['prune_baseline']} of {since['prune_entries']}.",
            "",
        ]
    lines += [
        "The classes below describe the runs' code. Stages B, V and X are rerun at the release commit before the paper "
        "quotes any class (docs/validation/freeze-checklist.md).",
        "",
    ]
    return lines


def render(doc: dict) -> str:
    a, b, vs = doc["stage_a"], doc["stage_b"], doc["view_stages"]
    stages = [("A", a), ("B", b), *((f"V ({v})", s) for v, s in vs.items())]
    hours = (sum(s["seconds"] for _, s in stages if s) + (doc["survey_seconds"] or 0)) / 3600
    thin = {  # view-stage constants read in one entry only: their class rests on that entry
        v: [r["id"] for r in doc["constants"] if r["stage"] == v and r["read_in"] == 1] for v in vs
    }
    membrane = vs.get("membrane", {}).get("entries", [])
    plan26 = [i for i in PLAN26_MEMBRANE if i in membrane]
    lines = [
        "# Sensitivity of the tuned constants",
        "",
        f"Plan 28 Task 9.2 (rulings Q29, Q30). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "Every module constant the inventory gives the origin `measured` or `named-entry` was set to x0.75 and x1.25 "
        "of its value, one at a time, and the random-100 campaign keys were measured on the third random sample, in "
        "the default view, against an unperturbed baseline run the same way: stage A on a subset, stage B on the "
        "whole sample for every constant. The third sample is not the tuning set (the first), but it has been run "
        "and looked at since plan 15, and some constants were set with its entries among those looked at (below: "
        "`named in the sample`, and the membrane constants); it is no unseen sample. The perturbation reaches every "
        "reader "
        "(scripts/perturb.py: the module is rewritten as it is imported, so `from ... import` binders, defaults and "
        "constants derived from it at import follow; `perturb.py reach` checked each run). Keys: each entry's pass "
        "or fail, the failing count per A- and S-criterion, the named fallbacks and the figure area of the view. "
        f"**stable**: no key moved; **sensitive**: some key moved; **fragile**: the pass count moved by "
        f"{doc['rules']['fragile_pass']} or more, or a criterion failing in no baseline entry failed in some entry. "
        "A constant takes the worse of its two perturbations.",
        "",
        "A reads run (unperturbed, traced: perturb.install_trace) records which entries ran a line, inside a "
        "function, that names each constant. A constant stable in the default view that no line of it read was "
        "measured again in the view that reads it (stage V: its entries, the default view drawn too, keyed on that "
        "view), found by a traced run of every view but the app's (its job process is not traced) on the whole sample "
        "(stage X); **not exercised** when no "
        "measurable view reads it, or its stage V reads run saw no read: its stability would be no evidence.",
        "",
        "**Time is recorded, not judged, and so is A2** (the default view drawn within 60 s of wall time, a criterion "
        "of random100's pass): it moved between two runs of one code on one machine with nothing else changed (see "
        "stage B below), so here an entry passes when no other A-criterion fails. A pass count random100 reports "
        "includes this load-dependent criterion.",
        "",
        "The baselines draw the default view (and stage V's view) only, so their pass counts are not the third "
        "sample's published ones (every view drawn, `docs/validation/random100-third/report.md`).",
        "",
        "**Limits of the method.** One factor at a time: no interaction between constants is measured. Two points "
        "per constant (x0.75, x1.25): a constant can move at another value and not at these. One sample (the third, "
        "100 entries): a class is this sample's. Stage A's subset leaves out the slowest entries, where the route "
        "search budgets (`route:_MAX_NODES_*`) may bind; stage B, on the whole sample, does not. A read is a line "
        "naming "
        "the constant that ran, not a decision it made. A run that crashed or whose value did not reach its readers "
        "is incomplete and gets no class. "
        + "".join(
            f"Stage V ({v}) has {len(st['entries'])} entries"
            + (
                f"; {len(thin[v])} of its constants were read in 1 of them ({', '.join(f'`{c}`' for c in thin[v])}), "
                "so their class rests on that one entry"
                if thin[v]
                else ""
            )
            + ". "
            for v, st in vs.items()
        )
        + f"The runs were made by {len(doc['run_commits'])} commits of the tooling ("
        + ", ".join(c[:7] for c in doc["run_commits"])
        + "), "
        "re-keyed by this report's code; every run's report names the one foldmap code fingerprint "
        f"{doc['code_fingerprint']} (checked).",
        "",
        "**Named in the sample**: the constants whose inventory origin *names* an entry of the third sample, a lower "
        "bound on the overlap (a `measured` origin names no entries): "
        + (", ".join(f"`{c}` ({', '.join(ids)})" for c, ids in doc["named_in_sample"].items()) or "none")
        + "."
        + (
            f" Besides, the membrane constants were tuned in plan 26 with {', '.join(plan26)} among the entries looked "
            f"at: {len(plan26)} of the membrane stage's {len(membrane)} entries."
            if plan26
            else ""
        ),
        "",
        *_since_lines(doc),
        "## Counts",
        "",
        "| class | constants |",
        "|---|---|",
        *[f"| {c} | {n} |" for c, n in doc["counts"].items() if c != "incomplete" or n],
        "",
        f"{len(doc['constants'])} constants perturbed, {len(doc['not_perturbable'])} not perturbable, "
        f"{len(doc['no_perturbation'])} with both perturbations skipped; zl time {hours:.1f} h ("
        + ", ".join(f"stage {n} {s['seconds'] / 3600:.2f} h" for n, s in stages if s)
        + (f", stage X {doc['survey_seconds'] / 3600:.2f} h" if doc["survey_seconds"] is not None else "")
        + ").",
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
    lines += ["", "## Not exercised", ""]
    lines += [
        f"- `{r['id']}` = {r['value']} ({r['origin']}): read by no view of random100 measured here"
        + (
            " (survey: " + ", ".join(f"{v} {n}" for v, n in sorted(r["survey"].items())) + ")"
            if r.get("survey")
            else " (survey: no process read it)"
            if "survey" in r
            else ""
        )
        for r in doc["constants"]
        if r["class"] == NOT_EXERCISED
    ] or ["None."]
    sub = doc["subset"]
    lines += [
        "",
        "## Stages",
        "",
        f"Stage A: {len(sub['ids'])} entries of the third sample (seed {sub['seed']}, target {sub['size']} entries; "
        "the quota per stratum counts entries, in proportion to its eligible entries, largest remainder; units are "
        "drawn whole, a pair being one unit): "
        + ", ".join(sub["ids"])
        + f". Left out, recorded slower than {doc['rules']['slow_s']} s in the third run (every view), or the partner "
        "of such an entry: "
        + ", ".join(
            f"{i} ({s:.0f} s{'' if s > doc['rules']['slow_s'] else ', partner'})"
            for i, s in sorted(sub["left_out"].items(), key=lambda x: -x[1])
        )
        + ".",
        "",
    ]
    for name, st in stages:
        if st:
            lines += [_stage_line(name, st), ""]
    if doc["survey_seconds"] is not None:
        lines += [
            f"Stage X: one unperturbed run of every view but the app's (its job process is not traced), traced, on the "
            f"whole sample; {doc['survey_seconds'] / 3600:.2f} h of zl time. It gave each constant the default view "
            "never reads to the view that reads it in the most entries (stage V), or to none (not exercised).",
            "",
        ]
    for v, st in vs.items():
        lines += [
            f"Stage V ({v}) measured {', '.join(f'`{c}`' for c in st['constants'])} on {', '.join(st['entries'])}.",
            "",
        ]
    for name, st in stages:
        rec = (st or {}).get("trace_records")
        if rec and rec["no_view_record"]:  # a view killed at its timeout leaves no record: its reads are unknown
            lines += [f"Stage {name}'s reads run left no view record for {', '.join(rec['no_view_record'])}.", ""]
    lines += [
        "## Per constant",
        "",
        "read in: the entries of the deciding stage (its view) a line naming the constant ran in, + when the judging "
        "read it. Each cell: the perturbed value, then each stage's class and what moved.",
        "",
        "| constant | origin | value | read in | x0.75 | x1.25 | class (stage) |",
        "|---|---|---|---|---|---|---|",
    ]
    pruned = set((doc.get("code_since") or {}).get("pruned", []))
    for r in doc["constants"]:
        cells = []
        for f in ("x0.75", "x1.25"):
            run = r["runs"].get(f)
            if run is None:
                cells.append("skipped")
                continue
            said = f"{run['new']}: " + "; ".join(f"{s} {_said(run[s])}" for s in run if s != "new" and run[s])
            cells.append(said.replace("|", "\\|"))
        st = {"A": a, "B": b}.get(r["stage"]) or vs.get(r["stage"])
        read = "-" if r["read_in"] is None else f"{r['read_in']}/{len(st['entries'])}" + "+" * r["read_judging"]
        gone = " (pruned after this study: prune record)" if r["id"] in pruned else ""
        lines.append(
            f"| `{r['id']}`{gone} | {r['origin']} | {r['value']} | {read} | {cells[0]} | {cells[1]} | "
            f"{r['class'] or 'incomplete'} ({r['stage']}) |"
        )
    lines += ["", "## Skipped perturbations", ""]
    lines += [f"- `{s['id']}` = {s['value']} x{s['factor']}: {s['why']}" for s in doc["skipped"]] or ["None."]
    lines += ["", "## Not perturbable", ""]
    lines += [f"- `{s['id']}` = {s['value']}: {s['why']}" for s in doc["not_perturbable"]] or ["None."]
    return "\n".join(lines) + "\n"


def cmd_report(args) -> int:
    views = [Path(v) for v in args.view or []]
    old = json.loads(Path(args.plan_from).read_text(encoding="utf-8")) if args.plan_from else None
    doc = build(
        Path(args.a),
        Path(args.b) if args.b else None,
        SOURCE_STAMP,
        views,
        Path(args.survey) if args.survey else None,
        the_plan=plan_of_record(old) if old else None,
    )
    if old:  # the plan, and so the inventory it was made from, is the record's
        doc["plan_from"] = {"report": str(args.plan_from), "inventory_sha256": old.get("inventory_sha256")}
        doc["inventory_sha256"] = old.get("inventory_sha256")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    (out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(f"{doc['counts']}; fragile: {doc['limitations']}; not exercised: {doc['not_exercised']}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan", help="the constants, their perturbations and stage A's subset")
    p.add_argument("--size", type=int, default=SUBSET)
    p.set_defaults(func=cmd_plan)
    r = sub.add_parser("run", help="run a stage's runs not run yet")
    r.add_argument("--stage", choices=("A", "B", "X", "V"), required=True)
    r.add_argument("--view", choices=VIEW_STAGES, help="stage V: the view")
    r.add_argument("--survey", help="stage V: stage X's folder")
    r.add_argument("--default-stage", help="stage V: stage B's folder (what the default view reads)")
    r.add_argument("--out", required=True, help="the stage's folder (never inside docs/validation/random100-third)")
    r.add_argument("--labels", help="baselines, perturbations or a comma-separated list of run names")
    r.add_argument("--size", type=int, default=SUBSET)
    r.add_argument("--jobs", type=int, default=8, help="entries drawn at once in each run")
    r.add_argument("--parallel", type=int, default=1, help="runs at once")
    r.set_defaults(func=cmd_run)
    w = sub.add_parser("report", help="judge the stages and write report.md and report.json")
    w.add_argument("--a", required=True)
    w.add_argument("--b")
    w.add_argument("--view", action="append", help="a stage V folder (may be given more than once)")
    w.add_argument("--survey", help="stage X's folder")
    w.add_argument("--out", default=str(OUT))
    w.add_argument(
        "--plan-from",
        metavar="REPORT",
        help="a committed report.json whose plan (constants, values, perturbations) the runs were made with: "
        "re-render from the kept runs when the code no longer has every constant measured (a pruned one)",
    )
    w.set_defaults(func=cmd_report)
    args = ap.parse_args(argv)
    for var in (ENV, TRACE_ENV):  # this process runs unperturbed and untraced; it sets them for each run's child
        if os.environ.get(var, "").strip():
            ap.error(f"{var} is set: this process must run unperturbed and untraced (it sets {var} per child)")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
