# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Prune only on evidence (plan 28 Task 9.3, ruling Q46): each special case set for a named entry is turned off, one
at a time, and measured; it is removed by a code commit only when, turned off, every measure is identical or better
and no test fails. Otherwise it is kept, and the record says which test or key kept it.

    python scripts/prune.py plan                                  # the candidates and how each is turned off
    python scripts/prune.py run --measure M --out DIR [--labels baseline|candidates|L1,L2] [--jobs 8] [--parallel 4]
    python scripts/prune.py report --runs DIR [--out docs/validation/prune]

The candidates (Q46): every inventory row (docs/validation/inventory/inventory.md) whose origin is named-entry (its
constants and its fallback rows), and the constants the sensitivity study (docs/validation/sensitivity/report.json)
found not exercised. Each is turned off by one value (DISABLE: the value that turns the special case off, set by
scripts/perturb.py's hook as the module is imported, so no code is edited), shares another's (SHARED: a fallback that
only that value's special case takes), or is not tried (NOT_TRIED: the hook sets an int or a float, ruling Q29, so a
tuple, a fallback or a code path is turned off only by removing code; the site is named).

Measures (M), each against an unperturbed baseline run the same way, all with the same disabling change:
  r100      the third random sample, default view, every entry: scripts/sensitivity.py's run_one (random100.py run
            --not-for-paper, FOLDMAP_PERTURB set; the keys: pass, the failing criteria, the named fallbacks, area), two
            baselines (they must be identical);
  membrane, domains   the view sensitivity's stage V measured a constant in (one the default view never reads), on the
            entries of that stage, the default view drawn first;
  na        the nucleic-acid controls and cases of plan 28's NA fixes (scripts/na_fixes_record.py: CASES, CONTROLS),
            default view, their A17 measures and named fallbacks (na_fixes_record._row and compare);
  domsets   the domain sets dom_cal, dom_held and dom_fresh (scripts/domain_uncertainty.py: their fixtures), each
            chain's overlap, F1 and recall against CATH as drawn;
  bench     scripts/bench.py --check BENCH_BEFORE --exact (every layout's values);
  suite     the full test suite (pytest -n), the value installed in every process by a plugin before foldmap is
            imported (FOLDMAP_PRUNE_PERTURB, never FOLDMAP_PERTURB: tests that run random100 run it as they do).
A measure's verdict: identical; better (something moved, and every move is one its own judge calls better: an entry
passing, a criterion no longer failing, fewer fallbacks, a judged NA measure the good way, a higher domain score, a
test passing); worse (any move the other way); moved (a change no judge orders: a figure's area, a layout's values,
an NA measure not judged); incomplete (a run with no report, or the value did not reach its readers). Pruned only
when every measure is identical or better. The fourth sample is never touched; the third is looked at, not unseen."""

from __future__ import annotations

import argparse
import contextlib
import datetime
import functools
import importlib.util
import io
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
INVENTORY = ROOT / "docs" / "validation" / "inventory" / "inventory.md"
SENSITIVITY = ROOT / "docs" / "validation" / "sensitivity" / "report.json"
THIRD = ROOT / "docs" / "validation" / "random100-third"
FOURTH = ROOT / "docs" / "validation" / "random100-fourth"
OUT = ROOT / "docs" / "validation" / "prune"
BENCH_BEFORE = ROOT / ".superpowers" / "sdd" / "2026-10-17-15-heldout-failures" / "bench-before.json"
PRUNE_ENV = "FOLDMAP_PRUNE_PERTURB"  # the suite's plugin reads it; FOLDMAP_PERTURB stays unset for the tests
MEASURES = ("r100", "membrane", "domains", "na", "domsets", "bench", "suite")
VIEWS = ("membrane", "domains")  # sensitivity's stage V views measured here
BASELINES = {"r100": ("baseline", "baseline-2")}  # every other measure: one baseline
OK = ("identical", "better")
VERDICTS = ("identical", "better", "moved", "worse", "incomplete")

# The value that turns each special case off, and why (Q46). "= X" names the rule the special case relaxes or tightens.
DISABLE = {
    "io:PAIRED_MIN": (10**9, "no two chains are ever held together by base pairs (io._paired never true)"),
    "io:_PAIR_DIST": (0.0, "no base atoms pair within 0 Å: io._paired never true"),
    "layout:_TAIL_COUNT_NEAR": (1e9, "a tail count need not stand beside its letter (layout.tail_counts)"),
    "layout:_SEP_MIRROR_ROOM": (
        0.0,
        "mirror-image copies get no more room than other copies drawn alike (layout._separate: max(room, 0)); their"
        " extra push passes stay, as they are code",
    ),
    "layout:_SEP_TIE_ROOM": (0.0, "the second chance for copies drawn alike is never kept (layout._tied_again)"),
    "layout:_STACK_SHUT": (
        0.0,
        "a stack gap needs no width to have a way out, so none is shut: no nudge, no move to the partner's other side,"
        " no pocket (layout._column_heights: they all ask way_out by _STACK_SHUT)",
    ),
    "layout:_STACK_NUDGE": (0.0, "a strand in a shut stack gap is never nudged (layout._column_heights.opened)"),
    "layout:_POCKET_REACH": (0.0, "an accessory closing a pocket never moves (layout._column_heights.unpocketed)"),
    "layout:_END_ROOM": (0.0, "a stack gap a chain end faces gets no more room (layout._column_heights.room)"),
    "layout:_BA_MAX": (10**9, "a βα row may hold any number of parallel strands (layout.can_be_row)"),
    "layout:_PART_BLOWUP": (1e9, "parting domains never starts over (layout._part_domains)"),
    "layout:_NAME_OPEN": (
        0.15,
        "= _NAME_ROUTE: a domain name keeps no more than its routing keep-out from a loop's way out (the 7NZ2 0.1 off)",
    ),
    "membrane:_MAX_POLAR_CROSSED": (0.15, "= _MAX_POLAR: a slab many elements cross is held to the ordinary polarity"),
    "membrane:_MANY_CROSSING": (10**9, "no slab is crossed by enough elements to relax its polarity limit"),
    "membrane:MIN_CROSSING_HELICES": (10**9, "the normal is never refitted to the crossing helices (membrane._refit)"),
    "membrane:FACE_FLAT_TILT": (45.0, "= FLAT_TILT: a helix at a face lies flat from the ordinary tilt"),
    "membrane:JOIN_GAP": (0, "no two helix runs are joined in a membrane view (membrane._apart: 1 <= gap <= 0)"),
    "na_runs:_JOINT_NEAR": (
        0.0,
        "a near joint no route reaches is held to the far joints' test (na_runs._straight_joint)",
    ),
    "na_runs:_JOINT_REACH": (0.0, "no far joint is drawn straight (na_runs._straight_joint)"),
    "na_runs:_BESIDE_TOLL": (0.0, "a route from beside a face pays nothing for cells under strands drawn before"),
    "nucleic:HBOND_STRETCH": (3.8, "= HBOND_MAX: no stretched pair is taken in the second pass (nucleic pairs)"),
    "panels:CLOSE": (0.0, "no gap between a domain's padded elements is filled (panels.panel_fields)"),
    "render:_SPOT_FAR": (1, "a ligand mark is never looked for beyond _SPOT_RADIUS (render._place_mark)"),
    "render:_SPOT_MORE_RINGS": (0, "no ring past the first clear spot is searched for fewer crossed bodies"),
    "route:_SLAB_RUN": (3.0, "= _SLAB_COST: a step along the slab costs what a step across it does"),
    "route:_MAX_NODES_PAGE": (10**12, "a page's second route searches have no shared budget (never refused)"),
    "route:_NA_TEXT_TOLL": (0.0, "a loop routed again off nucleic-acid text pays nothing for that text's cells"),
    "route:_REROUTE_TURNS_INK": (
        4,
        "= _REROUTE_TURNS: a reroute off text's ink may add no more corners than off its clearance",
    ),
    "symmetry:MIN_HELIX_COPIES": (3, "= MIN_HELIX_COPIES_ASKED: auto calls a screw of three copies helical"),
    "align:SAME_FOLD_TM": (0.0, "any aligned pair is one fold by its TM-score (align.element_mates; H1)"),
    "align:SAME_FOLD_FRACTION": (0.0, "any aligned pair is one fold by its matched fraction (align.element_mates; H1)"),
}
# a fallback only one DISABLE value's special case takes: measured by that value's runs
SHARED = {"route:_Router.second:fallback:route_second_search_refused": "route:_MAX_NODES_PAGE"}
_HOOK = "the hook sets an int or a float only (scripts/perturb.py, ruling Q29)"
NOT_TRIED = {
    "na_ends:_OFF_TURNS": f"na_ends.off_unheld: a tuple; {_HOOK}: turning the nudge off removes its search",
    "na_labels:_GAP_OUT": f"na_labels._label_gaps: a tuple; {_HOOK}: dropping 4 and 6 (7Y7C) edits the tuple",
    "na_labels:_GAP_TURNS": f"na_labels._gap_ring: a tuple; {_HOOK}: turning the ring off removes _gap_ring's call",
    "na_runs:_JOINT_BENDS": f"na_runs._joint_bend: a tuple; {_HOOK}: turning the bends off removes _joint_bend",
    "na_runs:_JOINT_SWINGS": f"na_runs._joint_bend: a tuple; {_HOOK}: turning the swings off edits the tuple",
    "route:_NAME_ALONG": f"route._clear_helix_names: a tuple; {_HOOK}: no value turns the moves along a helix off",
    "route:_MAX_NODES_WINDOWS": "route._Router.__init__ and route_loops: the special case is the windows' own pool "
    "(drawing from the page's budget gave 7NZ2 3 fallback curves); no value of it draws them from the page's budget "
    "again, which is code",
    "layout:_place_domain_names:fallback:domain_name_least_bad": "layout._place_domain_names: the name with no clear "
    "place must go somewhere or be left out; either is code",
    "layout:_place_domain_names:fallback-marker": "layout._place_domain_names: the marker of domain_name_least_bad "
    "(the same code)",
    "na_ends:off_unheld:fallback-marker": "na_ends.off_unheld, called by layout.build_layout: turning the 13FZ nudge "
    "off removes that call (its search constants are tuples)",
    "render:_draw:fallback:mark_least_bad_spot": "render._place_mark and _least_bad_spot: the mark with no clear spot "
    "went on its anchor before; that is code",
    "route:_clear_helix_names:fallback:helix_name_moved": "route._clear_helix_names, called by route.route_loops: "
    "turning the moves off removes that call",
    "route:<module>:fallback-marker": "route: route_unsnapped (_TOUCH), route_second_search_refused (measured by "
    "route:_MAX_NODES_PAGE) and the windows' pool: code paths",
}


# Candidates this task removed from the code (their inventory rows as measured, kept so the record still lists them)
PRUNED = {
    "layout:_NAME_OPEN": {"kind": "constant", "value": "_NAME_ROUTE + 0.1", "origin": "named-entry:7NZ2",
                          "tests": "domains"},
    "na_runs:_JOINT_NEAR": {"kind": "constant", "value": "2.5", "origin": "named-entry:9C75", "tests": "na_layout"},
}  # fmt: skip
# Candidates that met the prune condition and were not removed, and why (Q46's condition is necessary, not sufficient)
HELD = {
    "route:_MAX_NODES_PAGE": "it bounds the time a page's second route searches take (7QUP's focus view took 158 s of "
    "its 180), and no measure judges time; identical says these sets never spend the budget. Its named entry 7QUP's "
    "focus view drew identically without it, but 7NZ2's domains view (the held-out sample, which spends it: 6 "
    "searches refused) did not: A4 and A6 newly failed, loop_over_text 0 -> 1, helix_name_moved 0 -> 1 (one "
    "random100 run each way at 46f5ebf, outside Q46's measures). Removing it also removes route_second_search_refused "
    "and the readability floor's route_budget reason",
}
# The named entries of the candidates that met the condition, looked at after the measures (informational, not Q46's)
NAMED_LOOKS = (
    "layout:_NAME_OPEN, 7NZ2 (held-out sample) domains view: keys identical (pass, failing A17 only, fallbacks, "
    "area); its SVG differs (the name's place)",
    "na_runs:_JOINT_NEAR, 9C75 default view: identical, SVG too (9C75 is also an NA control, measured above)",
    "route:_MAX_NODES_PAGE, 7QUP focus view: identical, SVG too, 161.7 s both; 7NZ2 domains view: worse (see held)",
)


@functools.cache  # loaded once (the runs' threads share it)
def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"foldmap_prune_{name}", SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _stamp() -> dict:
    """This process's source stamp, taken when asked (a child installs its perturbation first)."""
    if (ROOT / "src" / "foldmap").is_dir() and str(ROOT / "src") not in sys.path:
        sys.path.insert(0, str(ROOT / "src"))  # the checkout's foldmap, as random100.py
    from foldmap.source_stamp import source_stamp

    return source_stamp(__file__)


# ------------------------------------------------------------------------------------------------- the candidates


def candidates(md: Path = INVENTORY, sensitivity: Path = SENSITIVITY) -> list[dict]:
    """Q46's candidates: the inventory's named-entry rows and the sensitivity study's not-exercised constants, with
    those this task pruned from the code (PRUNED, no longer in the inventory) first:
    [{"id", "kind", "value", "origin", "tests", "not_exercised"}]."""
    inv = _script("inventory")
    rows = inv.read_table(md.read_text(encoding="utf-8"))
    idle = json.loads(sensitivity.read_text(encoding="utf-8"))["not_exercised"]
    out = [
        {"id": i, **r, "not_exercised": i in idle, "pruned_in_code": True} for i, r in PRUNED.items() if i not in rows
    ]
    for i, r in rows.items():
        if inv.origin_kind(r["origin"]) == "named-entry" or i in idle:
            out.append(
                {
                    "id": i,
                    "kind": r["kind"],
                    "value": r["value"],
                    "origin": r["origin"],
                    "tests": r.get("tests naming it", ""),
                    "not_exercised": i in idle,
                }
            )
    return out


def label(cid: str) -> str:
    """A candidate's run name (a path: sensitivity.safe_name), "module.NAME-off"."""
    return _script("sensitivity").safe_name(f"{cid.replace(':', '.')}-off")


def perturb_text(cid: str) -> str:
    """FOLDMAP_PERTURB's text for a DISABLE candidate, checked by perturb.parse."""
    module, name = cid.split(":", 1)
    text = f"{module}:{name}={DISABLE[cid][0]!r}"
    _script("perturb").parse(text)
    return text


def view_of(cid: str, sensitivity: Path = SENSITIVITY) -> str | None:
    """The view sensitivity's stage V measured a constant in (one the default view never reads), or None."""
    stages = json.loads(sensitivity.read_text(encoding="utf-8"))["view_stages"]
    return next((v for v in VIEWS if cid in (stages.get(v, {}).get("constants") or [])), None)


def view_entries(view: str, sensitivity: Path = SENSITIVITY) -> list[str]:
    return list(json.loads(sensitivity.read_text(encoding="utf-8"))["view_stages"][view]["entries"])


def plan(items: list[dict] | None = None) -> dict:
    """Each candidate with how it is turned off: {"tried": [...], "shared": [...], "not_tried": [...], "unplanned"}."""
    items = candidates() if items is None else items
    tried, shared, not_tried, unplanned = [], [], [], []
    for c in items:
        if c["id"] in DISABLE:
            value, why = DISABLE[c["id"]]
            tried.append({**c, "off": value, "why": why, "label": label(c["id"]), "perturb": perturb_text(c["id"]),
                          "view": view_of(c["id"])})  # fmt: skip
        elif c["id"] in SHARED:
            shared.append({**c, "with": SHARED[c["id"]]})
        elif c["id"] in NOT_TRIED:
            not_tried.append({**c, "site": NOT_TRIED[c["id"]]})
        else:
            unplanned.append(c["id"])
    return {"tried": tried, "shared": shared, "not_tried": not_tried, "unplanned": unplanned}


def cmd_plan(args) -> int:
    p = plan()
    for t in p["tried"]:
        print(f"  {t['id']} = {t['value']} -> {t['off']!r}{' [' + t['view'] + ']' if t['view'] else ''}: {t['why']}")
    for s in p["shared"]:
        print(f"  {s['id']}: measured with {s['with']}")
    for n in p["not_tried"]:
        print(f"  not tried {n['id']}: {n['site']}")
    print(f"{len(p['tried'])} tried, {len(p['shared'])} shared, {len(p['not_tried'])} not tried")
    if p["unplanned"]:
        print("UNPLANNED:", p["unplanned"])
        return 1
    return 0


# ------------------------------------------------------------------------------------------------- children


PLUGIN = """# written by scripts/prune.py: the suite's processes (and xdist's workers) install the perturbation first
import importlib.util
import os

_spec = importlib.util.spec_from_file_location("foldmap_prune_perturb", {perturb!r})
_perturb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_perturb)
_perturb.install(os.environ.get({env!r}, ""))
"""


def _install(text: str) -> None:
    """Install `text`'s perturbation in this process before foldmap's module is imported ("" installs nothing)."""
    if text:
        _script("perturb").install(text)


def child_domsets(text: str) -> dict:
    """Each chain of the domain sets scored as drawn (domain_uncertainty's fixtures, _domains, _score), with `text`
    installed first."""
    _install(text)
    if str(ROOT / "src") not in sys.path and (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    du = _script("domain_uncertainty")
    ds = du.ds
    out = {}
    for name, pairs in du.fixture_sets().items():
        got = {}
        for inputs, e in pairs:
            chain = e["chain"]
            ref = ds.reference_labels([lab.seq for lab in inputs.labels if lab.chain == chain], e["domains"])
            auto = ds.layout_labels(du._domains(inputs), inputs.sses, inputs.labels, chain)
            s = du._score(ref, auto)
            got[f"{e['pdb']}_{chain}"] = {k: s[k] for k in ("overlap", "f1", "recall", "meets")}
        out[name] = got
    return {"sets": out}


def child_bench(text: str) -> dict:
    """scripts/bench.py --check BENCH_BEFORE --exact with `text` installed first: its exit and each entry's note."""
    _install(text)
    bench = _script("bench")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
        code = bench.main(["--check", str(BENCH_BEFORE), "--exact"])
    notes = {}
    for line in buf.getvalue().splitlines():
        parts = line.split()
        if parts:
            notes[parts[0]] = "identical" if line.rstrip().endswith("identical") else line.split("  ")[-1].strip()
    return {"exit": code, "changed": sorted(n for n, v in notes.items() if v != "identical"), "entries": len(notes)}


def cmd_child(args) -> int:
    t0 = time.perf_counter()
    got = {"domsets": child_domsets, "bench": child_bench}[args.what](args.perturb)
    module = args.perturb.split(":", 1)[0] if args.perturb else None
    got["reached"] = None if module is None else f"foldmap.{module}" in sys.modules
    got["perturb"] = args.perturb or None
    got["seconds"] = round(time.perf_counter() - t0, 1)
    got["source_stamp"] = _stamp()
    Path(args.to).write_text(json.dumps(got, indent=1) + "\n", encoding="utf-8")
    return 0


# ------------------------------------------------------------------------------------------------- running


def _run_child(what: str, name: str, text: str | None, out: Path) -> dict:
    """One domsets or bench run in its own process, unperturbed but for `text` (FOLDMAP_PERTURB never set)."""
    out.mkdir(parents=True, exist_ok=True)
    to = out / f"{name}.json"
    env = {k: v for k, v in os.environ.items() if k not in ("FOLDMAP_PERTURB", PRUNE_ENV)}
    cmd = [sys.executable, "-I", str(Path(__file__).resolve()), "child", what, "--to", str(to)]
    cmd += ["--perturb", text or ""]
    t0 = time.perf_counter()
    with open(out / f"{name}.log", "w", encoding="utf-8") as log:
        code = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT).returncode
    doc = json.loads(to.read_text(encoding="utf-8")) if code == 0 and to.is_file() else {"error": f"exit {code}"}
    doc["wall"] = round(time.perf_counter() - t0, 1)
    if "source_stamp" not in doc:  # the child stamps what it measured; a child that failed did not
        doc["source_stamp"] = _stamp()
    to.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    return doc


def _run_suite(name: str, text: str | None, out: Path, jobs: int) -> dict:
    """The full test suite with `text` installed in every process by a plugin (PLUGIN), xdist over `jobs` workers."""
    out.mkdir(parents=True, exist_ok=True)
    plug = out / "plugin"
    plug.mkdir(exist_ok=True)
    (plug / "foldmap_prune_plugin.py").write_text(
        PLUGIN.format(perturb=str(SCRIPTS / "perturb.py"), env=PRUNE_ENV), encoding="utf-8"
    )
    env = {k: v for k, v in os.environ.items() if k not in ("FOLDMAP_PERTURB", PRUNE_ENV, "PYTHONPATH")}
    env["PYTHONPATH"] = str(plug)
    if text:
        env[PRUNE_ENV] = text
    xml = out / f"{name}.xml"
    cmd = [sys.executable, "-m", "pytest", "-q", "-n", str(jobs), "-p", "no:cacheprovider", "-p"]
    cmd += ["foldmap_prune_plugin", f"--junitxml={xml}", "-o", "junit_family=xunit2"]
    t0 = time.perf_counter()
    with open(out / f"{name}.log", "w", encoding="utf-8") as log:
        code = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT).returncode
    doc = {"exit": code, "perturb": text, "seconds": round(time.perf_counter() - t0, 1), **junit(xml)}
    doc["source_stamp"] = _stamp()  # the tests ran the checkout's code as this process found it
    (out / f"{name}.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    xml.unlink(missing_ok=True)
    return doc


def junit(path: Path) -> dict:
    """A JUnit XML file's counts and the ids of the tests that failed or erred ({"error"} when there is none)."""
    if not path.is_file():
        return {"error": "no junit xml"}
    root = ET.parse(path).getroot()
    failed, counts = [], {"tests": 0, "failed": 0, "skipped": 0}
    for case in root.iter("testcase"):
        counts["tests"] += 1
        tid = f"{case.get('classname')}::{case.get('name')}"
        if case.find("failure") is not None or case.find("error") is not None:
            failed.append(tid)
        elif case.find("skipped") is not None:
            counts["skipped"] += 1
    counts["failed"] = len(failed)
    return {**counts, "failures": sorted(failed)}


def runs_of(measure: str, the_plan: dict) -> list[tuple[str, str | None]]:
    """(name, perturbation) of every run of a measure: its baselines, then each candidate it measures."""
    runs = [(b, None) for b in BASELINES.get(measure, ("baseline",))]
    for t in the_plan["tried"]:
        if measure in VIEWS and t["view"] != measure:
            continue
        runs.append((t["label"], t["perturb"]))
    return runs


def cmd_run(args) -> int:
    out = Path(args.out).resolve()
    for refused in (THIRD, FOURTH):
        if out.is_relative_to(refused.resolve()):
            raise SystemExit(f"refused: {out} is inside {refused}")
    the_plan = plan()
    if the_plan["unplanned"]:
        raise SystemExit(f"refused: candidates with no plan: {the_plan['unplanned']}")
    runs = runs_of(args.measure, the_plan)
    if args.labels == "baseline":
        runs = [r for r in runs if r[1] is None]
    elif args.labels == "candidates":
        runs = [r for r in runs if r[1] is not None]
    elif args.labels:
        wanted = set(args.labels.split(","))
        runs = [r for r in runs if r[0] in wanted]
    where = out / args.measure
    keyed = args.measure in ("r100", "na", *VIEWS)  # sensitivity.run_one writes keys/NAME.json
    todo = [r for r in runs if not (where / "keys" / f"{r[0]}.json" if keyed else where / f"{r[0]}.json").is_file()]
    print(f"{args.measure}: {len(runs)} runs, {len(todo)} to run", flush=True)
    sens = _script("sensitivity")

    def one(r):
        name, text = r
        if args.measure == "r100":
            doc = sens.run_one(name, text, None, where, args.jobs, THIRD / "sample.json")
            said = f"{sens.passed(doc['entries']) if doc['entries'] else '-'} pass, exit {doc['exit']}"
        elif args.measure in VIEWS:
            ids = view_entries(args.measure)
            doc = sens.run_one(name, text, ids, where, args.jobs, THIRD / "sample.json", None,
                               f"default,{args.measure}", args.measure)  # fmt: skip
            said = f"{sens.passed(doc['entries']) if doc['entries'] else '-'} pass, exit {doc['exit']}"
        elif args.measure == "na":
            na = _script("na_fixes_record")
            ids = [*na.CASES, *na.CONTROLS]
            doc = sens.run_one(name, text, ids, where, args.jobs, where / "no-sample.json")
            said = f"exit {doc['exit']}"
        elif args.measure == "suite":
            doc = _run_suite(name, text, where, args.jobs)
            said = f"exit {doc['exit']}, {doc.get('failed')} failed of {doc.get('tests')}"
        else:
            doc = _run_child(args.measure, name, text, where)
            said = doc.get("error") or f"reached {doc.get('reached')}"
        print(f"{args.measure} {name}: {said}", flush=True)

    t0 = time.perf_counter()
    parallel = 1 if args.measure == "suite" else max(1, args.parallel)
    with ThreadPoolExecutor(max_workers=parallel) as pool:
        list(pool.map(one, todo))
    where.mkdir(parents=True, exist_ok=True)
    with open(where / "stages.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps({"runs": [r[0] for r in todo], "seconds": round(time.perf_counter() - t0, 1)}) + "\n")
    return 0


# ------------------------------------------------------------------------------------------------- verdicts


def _verdict(worse: list, better: list, moved: list) -> str:
    return "worse" if worse else "moved" if moved else "better" if better else "identical"


def campaign_verdict(base: dict, run: dict) -> dict:
    """A random-100 run's keys (sensitivity.keys_of) against its baseline's: an entry passing or failing, a criterion
    failing or no longer, a named fallback's count, the figure's area (moved: no judge orders it)."""
    worse, better, moved = [], [], []
    for i in sorted(set(base) | set(run)):
        a, b = base.get(i), run.get(i)
        if a is None or b is None:
            worse.append(f"{i}: in one run only")
            continue
        if a["pass"] != b["pass"]:
            (better if b["pass"] else worse).append(f"{i}: {'FAIL -> PASS' if b['pass'] else 'PASS -> FAIL'}")
        for k in sorted(set(a["failing"]) ^ set(b["failing"])):
            now = k in b["failing"]
            (worse if now else better).append(f"{i}: {k} {'now fails' if now else 'now passes'}")
        for n in sorted(set(a["fallbacks"]) | set(b["fallbacks"])):
            x, y = a["fallbacks"].get(n, 0), b["fallbacks"].get(n, 0)
            if x != y:
                (worse if y > x else better).append(f"{i}: fallback {n} {x} -> {y}")
        if a["area"] != b["area"]:
            moved.append(f"{i}: area {a['area']} -> {b['area']}")
    return {"verdict": _verdict(worse, better, moved), "worse": worse, "better": better, "moved": moved}


def na_verdict(before: dict, after: dict) -> dict:
    """The NA entries (na_fixes_record._row of each) against the baseline's: na_fixes_record.compare's worse measures
    are worse; a judged measure that moved the good way, fewer fallbacks or an error gone are better; any other
    change (a measure not judged, or one appearing or going) moved."""
    na = _script("na_fixes_record")
    worse, better, moved = [], [], []
    for i, d in na.compare(before, after).items():
        for k, (a, b) in d["changed"].items():
            line = f"{i}: {k} {a} -> {b}"
            if k in d["worse"]:
                worse.append(line)
            elif k == "named_fallbacks" or (k == "error" and a and not b):
                better.append(line)
            elif k in na.JUDGED and na._num(a) and na._num(b) and na.worse(k, b, a):
                better.append(line)
            else:
                moved.append(line)
    return {"verdict": _verdict(worse, better, moved), "worse": worse, "better": better, "moved": moved}


def domsets_verdict(base: dict, run: dict) -> dict:
    """Each chain's overlap, F1 and recall against the baseline's: lower is worse, higher better."""
    worse, better = [], []
    for name, chains in base["sets"].items():
        for c, a in chains.items():
            b = run["sets"].get(name, {}).get(c)
            if b is None:
                worse.append(f"{name} {c}: not scored")
                continue
            for k in ("overlap", "f1", "recall"):
                if b[k] != a[k]:
                    (worse if b[k] < a[k] else better).append(f"{name} {c}: {k} {a[k]} -> {b[k]}")
    return {"verdict": _verdict(worse, better, []), "worse": worse, "better": better, "moved": []}


def bench_verdict(base: dict, run: dict) -> dict:
    """bench --exact: every layout as recorded (identical), else the entries whose layout moved (no judge orders a
    layout's values)."""
    moved = [f"{n}: layout changed" for n in run["changed"] if n not in base["changed"]]
    worse = [] if moved or run.get("exit") == base.get("exit") else [f"bench exit {run.get('exit')}"]
    return {"verdict": _verdict(worse, [], moved), "worse": worse, "better": [], "moved": moved}


def suite_verdict(base: dict, run: dict) -> dict:
    """Tests failing that pass in the baseline are worse; failing in the baseline and passing here, better."""
    was, now = set(base["failures"]), set(run["failures"])
    worse = sorted(now - was)
    better = sorted(was - now)
    if run.get("tests", 0) < base.get("tests", 0) and not worse:
        worse = [f"{base['tests'] - run['tests']} tests not run (collection)"]
    return {"verdict": _verdict(worse, better, []), "worse": worse, "better": better, "moved": []}


def _load(path: Path) -> dict | None:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def _campaign_keys(where: Path, name: str, measure: str) -> dict | None:
    """A random-100 run's keys file, its entries keyed again by sensitivity.keys_of from its kept report (the NA
    measure: na_fixes_record._row of each entry)."""
    doc = _load(where / "keys" / f"{name}.json")
    if doc is None:
        return None
    rep = where / "runs" / name / "report.json"
    rd = _script("random100").read_report(rep) if rep.is_file() or Path(f"{rep}.gz").is_file() else None
    if rd is None:
        doc["entries"] = None
    elif measure == "na":
        na = _script("na_fixes_record")
        doc["entries"] = {e["id"]: na._row(e) for e in rd["entries"]}
    else:
        doc["entries"] = _script("sensitivity").keys_of(rd, "default" if measure == "r100" else measure)
    return doc


def judge_measure(measure: str, where: Path, name: str, base_name: str = "baseline") -> dict | None:
    """A run's verdict against the measure's baseline; None when the run is missing."""
    if measure in ("r100", "na", *VIEWS):
        base, run = _campaign_keys(where, base_name, measure), _campaign_keys(where, name, measure)
        if run is None:
            return None
        if base is None or base["entries"] is None:
            return {"verdict": "incomplete", "why": "no baseline"}
        if run["entries"] is None:
            return {"verdict": "incomplete", "why": f"no report (exit {run['exit']})"}
        if run.get("perturb") and not (run.get("reach") or {}).get("reached"):
            return {"verdict": "incomplete", "why": "the value did not reach its readers"}
        if measure == "na":
            got = na_verdict({"entries": base["entries"]}, {"entries": run["entries"]})
        else:
            got = campaign_verdict(base["entries"], run["entries"])
        got["code"] = (run.get("report_stamp") or {}).get("code_fingerprint")
        return got
    base, run = _load(where / f"{base_name}.json"), _load(where / f"{name}.json")
    if run is None:
        return None
    if base is None or base.get("error"):
        return {"verdict": "incomplete", "why": "no baseline"}
    if run.get("error"):
        return {"verdict": "incomplete", "why": run["error"]}
    got = {"domsets": domsets_verdict, "bench": bench_verdict, "suite": suite_verdict}[measure](base, run)
    if measure in ("domsets", "bench"):
        got["reached"] = run.get("reached")
        got["code"] = (run.get("source_stamp") or {}).get("code_fingerprint")
    return got


# ------------------------------------------------------------------------------------------------- the record


def decide(measures: dict[str, dict | None], wanted: list[str]) -> dict:
    """Q46's condition over one candidate's measures (`wanted`: the measures it must have): prune when every one
    was run and is identical or better; else keep, with what kept it (or incomplete when a measure is missing)."""
    missing = [m for m in wanted if measures.get(m) is None]
    bad = {m: v for m, v in measures.items() if v is not None and v["verdict"] not in OK}
    if bad:
        kept_by = []
        for m, v in bad.items():
            what = v.get("worse") or v.get("moved") or [v.get("why", v["verdict"])]
            kept_by.append({"measure": m, "verdict": v["verdict"], "first": what[:3], "count": len(what)})
        return {"decision": "keep", "kept_by": kept_by, "missing": missing}
    if missing:
        return {"decision": "incomplete", "missing": missing}
    return {"decision": "prune", "kept_by": [], "missing": []}


def wanted_measures(t: dict) -> list[str]:
    return ["r100", *([t["view"]] if t["view"] else []), "na", "domsets", "bench", "suite"]


def _baseline_checks(runs: Path) -> dict:
    """The baselines: r100's two identical (the measure is deterministic), the bench as recorded in BENCH_BEFORE, the
    suite passing, each measure's baseline present."""
    out = {}
    r1, r2 = _campaign_keys(runs / "r100", "baseline", "r100"), _campaign_keys(runs / "r100", "baseline-2", "r100")
    if r1 and r2 and r1["entries"] and r2["entries"]:
        lines = _script("sensitivity").moved(r1["entries"], r2["entries"])
        e1, e2 = r1["entries"], r2["entries"]
        svg = sorted(i for i in e1 if e1[i]["svg_sha256"] != (e2.get(i) or {}).get("svg_sha256"))
        out["r100_deterministic"] = {"identical": not lines and not svg, "moved": lines, "svg": svg}
        out["r100_passed"] = _script("sensitivity").passed(r1["entries"])
        out["r100_entries"] = len(r1["entries"])
    for v in VIEWS:
        b = _campaign_keys(runs / v, "baseline", v)
        if b and b["entries"]:
            out[f"{v}_passed"], out[f"{v}_entries"] = _script("sensitivity").passed(b["entries"]), len(b["entries"])
    na = _campaign_keys(runs / "na", "baseline", "na")
    if na and na["entries"]:
        out["na_entries"] = len(na["entries"])
        out["na_failed"] = sorted(i for i, r in na["entries"].items() if r.get("na_failed"))
    bench = _load(runs / "bench" / "baseline.json")
    if bench:
        out["bench"] = {"exit": bench.get("exit"), "changed": bench.get("changed"), "entries": bench.get("entries")}
    dom = _load(runs / "domsets" / "baseline.json")
    if dom and "sets" in dom:
        out["domsets_chains"] = {k: len(v) for k, v in dom["sets"].items()}
    suite = _load(runs / "suite" / "baseline.json")
    if suite:
        out["suite"] = {k: suite.get(k) for k in ("exit", "tests", "failed", "skipped", "failures", "seconds")}
    return out


def _seconds(runs: Path) -> dict:
    out = {}
    for m in MEASURES:
        f = runs / m / "stages.jsonl"
        if f.is_file():
            out[m] = round(sum(json.loads(t)["seconds"] for t in f.read_text(encoding="utf-8").splitlines() if t), 1)
    return out


def build(runs: Path, stamp: dict, the_plan: dict | None = None) -> dict:
    """The record: each candidate's measures, verdicts and decision; the shared and the not tried; the baselines."""
    the_plan = plan() if the_plan is None else the_plan
    rows, codes = [], set()
    for t in the_plan["tried"]:
        measures = {}
        for m in wanted_measures(t):
            got = judge_measure(m, runs / m, t["label"])
            measures[m] = got
            if got and got.get("code"):
                codes.add(got["code"])
        row = {k: t[k] for k in ("id", "kind", "value", "origin", "tests", "not_exercised", "off", "why", "perturb")}
        row |= {"view": t["view"], "measures": measures, **decide(measures, wanted_measures(t))}
        if row["decision"] == "prune" and t["id"] in HELD:
            row |= {"decision": "held", "held": HELD[t["id"]]}
        row["removed"] = t["id"] in PRUNED
        rows.append(row)
    by_id = {r["id"]: r for r in rows}
    shared = [
        {**{k: s[k] for k in ("id", "kind", "value", "origin", "tests", "not_exercised")}, "with": s["with"],
         "decision": by_id[s["with"]]["decision"], "kept_by": by_id[s["with"]].get("kept_by", [])}
        for s in the_plan["shared"]
    ]  # fmt: skip
    not_tried = [{k: n[k] for k in ("id", "kind", "value", "origin", "tests", "not_exercised", "site")}
                 for n in the_plan["not_tried"]]  # fmt: skip
    if len(codes) > 1:
        raise SystemExit(f"refused: the runs name {len(codes)} foldmap code fingerprints {sorted(codes)}")
    kept_by = {}
    for r in [*rows, *shared]:
        if r["decision"] == "keep":
            for m in sorted({k["measure"] for k in r["kept_by"]}):
                kept_by[m] = kept_by.get(m, 0) + 1
    decisions = [r["decision"] for r in [*rows, *shared]]
    return {
        "what": "special cases set for a named entry, turned off one at a time and measured; pruned only when every "
        "measure is identical or better and no test fails (scripts/prune.py, plan 28 Task 9.3, ruling Q46)",
        "date": datetime.date.today().isoformat(),
        "source_stamp": stamp,
        "code_fingerprint": codes.pop() if codes else None,
        "inventory_sha256": _sha256(INVENTORY),
        "sensitivity_sha256": _sha256(SENSITIVITY),
        "bench_before": str(BENCH_BEFORE.relative_to(ROOT)),
        "counts": {
            "candidates": len(rows) + len(shared) + len(not_tried),
            "tried": len(rows) + len(shared),
            "tried_values": len(rows),
            "shared": len(shared),
            "kept": decisions.count("keep"),
            "pruned": decisions.count("prune"),
            "held": decisions.count("held"),
            "incomplete": decisions.count("incomplete"),
            "not_tried": len(not_tried),
            "kept_by_measure": dict(sorted(kept_by.items())),
        },
        "baselines": _baseline_checks(runs),
        "seconds": _seconds(runs),
        "candidates": rows,
        "shared": shared,
        "not_tried": not_tried,
        "pruned": [r["id"] for r in rows if r["decision"] == "prune"],
        "removed": sorted(PRUNED),
        "held": {i: HELD[i] for i in HELD if by_id.get(i, {}).get("decision") == "held"},
        "named_looks": list(NAMED_LOOKS),
    }


def _sha256(path: Path) -> str | None:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _cell(v: dict | None) -> str:
    if v is None:
        return "-"
    said = v["verdict"]
    what = v.get("worse") or v.get("moved") or v.get("better") or []
    if v.get("why"):
        said += f" ({v['why']})"
    elif what:
        said += f" ({len(what)}: {what[0]}{'; …' if len(what) > 1 else ''})"
    if v.get("reached") is False:
        said += ", not imported"
    return said.replace("|", "\\|")


def render(doc: dict) -> str:
    from foldmap.source_stamp import stamp_line

    c = doc["counts"]
    lines = [
        "# Pruning on evidence",
        "",
        f"Plan 28 Task 9.3 (ruling Q46). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        'The review asked whether foldmap\'s special cases "may have emerged from fixing individual examples". Each '
        "special case the inventory (`docs/validation/inventory/inventory.md`) says was set for a named entry, and "
        "each constant the sensitivity study found not exercised, was turned off by one value, one at a time "
        "(`scripts/perturb.py`'s hook rewrites the module as it is imported, so every reader takes the value and no "
        "code is edited), and measured against an unperturbed baseline: the third random sample's campaign keys in "
        "the default view (pass, failing criteria, named fallbacks, figure area), the view that reads it where the "
        "default view does not (sensitivity's stage V entries), the nucleic-acid cases and controls of plan 28's NA "
        "fixes, the domain sets (dom_cal, dom_held, dom_fresh: each chain's overlap, F1 and recall against CATH), "
        "`scripts/bench.py --exact`, and the full test suite with the value installed in every process. **A special "
        "case is removed only when every measure is identical or better and no test fails** (`scripts/prune.py`'s "
        "docstring defines each verdict); otherwise it is kept, and the table says which key or test kept it. The "
        "third sample has been run and looked at since plan 15: it is no unseen sample. The fourth sample was not "
        "touched.",
        "",
        "**Limits.** One value per special case, one at a time: no interaction is measured. Turning off is not "
        "always the code as it was before the special case (each value's reason says what it turns off). A view "
        "other than the default is measured only for a constant the default view never reads, on the entries "
        "sensitivity's stage V drew it in. A test that runs foldmap in a child process (the CLI, random100) runs it "
        "unperturbed. Identical on every measure means the measured sets never take the special case, not that its "
        "named entry needs it no longer.",
        "",
        "## Counts",
        "",
        "| | candidates |",
        "|---|---|",
        f"| all | {c['candidates']} |",
        f"| tried (one value each, {c['shared']} measured with another's value) | {c['tried']} |",
        f"| kept | {c['kept']} |",
        f"| pruned | {c['pruned']} |",
        f"| met the condition, held | {c['held']} |",
        f"| incomplete | {c['incomplete']} |",
        f"| not tried: needs code removal | {c['not_tried']} |",
        "",
        "Kept by measure (a candidate kept by several counts in each): "
        + (", ".join(f"{m} {n}" for m, n in c["kept_by_measure"].items()) or "none")
        + ".",
        "",
    ]
    b = doc["baselines"]
    lines += ["## Baselines", ""]
    if "r100_deterministic" in b:
        d = b["r100_deterministic"]
        lines.append(
            f"- r100: {b['r100_entries']} entries, {b['r100_passed']} pass (A2, wall time, not judged); the two "
            f"baselines {'identical (keys and SVGs)' if d['identical'] else 'DIFFER: ' + '; '.join(d['moved'][:5])}."
        )
    for v in VIEWS:
        if f"{v}_passed" in b:
            lines.append(f"- {v}: {b[f'{v}_entries']} entries, {b[f'{v}_passed']} pass.")
    if "na_entries" in b:
        lines.append(f"- na: {b['na_entries']} entries; failing A17: {', '.join(b['na_failed']) or 'none'}.")
    if "domsets_chains" in b:
        lines.append("- domsets: " + ", ".join(f"{k} {n} chains" for k, n in b["domsets_chains"].items()) + ".")
    if "bench" in b:
        lines.append(
            f"- bench --exact against `{doc['bench_before']}`: exit {b['bench']['exit']}, {b['bench']['entries']} "
            f"entries, changed: {', '.join(b['bench']['changed']) or 'none'}."
        )
    if "suite" in b:
        s = b["suite"]
        lines.append(
            f"- suite: {s['tests']} tests, {s['failed']} failed, {s['skipped']} skipped, exit {s['exit']}"
            + (f" (failing: {', '.join(s['failures'])})" if s["failures"] else "")
            + "."
        )
    secs = doc["seconds"]
    lines += ["", "zl time: " + ", ".join(f"{m} {s / 3600:.2f} h" for m, s in secs.items()) + ".", ""]
    lines += ["## Pruned", ""]
    lines += [f"- `{i}`" + (": removed from the code" if i in doc["removed"] else "") for i in doc["pruned"]] or [
        "None."
    ]
    lines += ["", "## Met the condition, held", ""]
    lines += [f"- `{i}`: {why}" for i, why in doc["held"].items()] or ["None."]
    lines += ["", "## Named entries looked at (informational, not a measure of Q46)", ""]
    lines += [f"- {t}" for t in doc["named_looks"]]
    lines += ["", "## Tried", "", "Each cell: the verdict, how many moves, the first.", ""]
    lines += ["| candidate | origin | value -> off | decision | kept by | r100 | view | na | domsets | bench | suite |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]  # fmt: skip
    for r in doc["candidates"]:
        m = r["measures"]
        view = f"{r['view']}: {_cell(m.get(r['view']))}" if r["view"] else "-"
        kept = ", ".join(sorted({k["measure"] for k in r.get("kept_by", [])})) or "-"
        cells = [_cell(m.get(k)) for k in ("r100", "na", "domsets", "bench", "suite")]
        lines.append(
            f"| `{r['id']}` | {r['origin']} | {r['value']} -> {r['off']!r} | **{r['decision']}** | {kept} | "
            f"{cells[0]} | {view} | {' | '.join(cells[1:])} |"
        )
    lines += ["", "### How each was turned off", ""]
    lines += [f"- `{r['id']}` = {r['off']!r}: {r['why']}" for r in doc["candidates"]]
    lines += ["", "### What kept each", ""]
    for r in doc["candidates"]:
        if r["decision"] != "keep":
            continue
        parts = [f"{k['measure']} {k['verdict']} ({k['count']}: {'; '.join(k['first'])})" for k in r["kept_by"]]
        lines.append(f"- `{r['id']}`: " + " · ".join(parts))
    lines += ["", "## Measured with another's value", ""]
    lines += [f"- `{s['id']}`: with `{s['with']}`, **{s['decision']}**" for s in doc["shared"]] or ["None."]
    lines += ["", "## Not tried: needs code removal", ""]
    lines += [f"- `{n['id']}` ({n['kind']}, {n['origin']}): {n['site']}" for n in doc["not_tried"]] or ["None."]
    return "\n".join(lines) + "\n"


def cmd_report(args) -> int:
    stamp = _stamp()
    doc = build(Path(args.runs), stamp)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8", newline="\n")
    (out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(json.dumps(doc["counts"]), "pruned:", doc["pruned"])
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("plan", help="the candidates and how each is turned off").set_defaults(func=cmd_plan)
    r = sub.add_parser("run", help="run a measure's runs not run yet")
    r.add_argument("--measure", choices=MEASURES, required=True)
    r.add_argument("--out", required=True, help="the runs' folder (a folder per measure inside it)")
    r.add_argument("--labels", help="baseline, candidates or a comma-separated list of run names")
    r.add_argument("--jobs", type=int, default=8, help="entries drawn at once in a run (suite: xdist workers)")
    r.add_argument("--parallel", type=int, default=1, help="runs at once (the suite: always one)")
    r.set_defaults(func=cmd_run)
    w = sub.add_parser("report", help="judge the runs and write report.md and report.json")
    w.add_argument("--runs", required=True)
    w.add_argument("--out", default=str(OUT))
    w.set_defaults(func=cmd_report)
    c = sub.add_parser("child", help="(internal) one domsets or bench run, the perturbation installed first")
    c.add_argument("what", choices=("domsets", "bench"))
    c.add_argument("--perturb", default="")
    c.add_argument("--to", required=True)
    c.set_defaults(func=cmd_child)
    args = ap.parse_args(argv)
    if args.command != "child":
        for var in ("FOLDMAP_PERTURB", "FOLDMAP_TRACE_READS", PRUNE_ENV):
            if os.environ.get(var, "").strip():
                ap.error(f"{var} is set: this process runs unperturbed (it sets the perturbation per child)")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
