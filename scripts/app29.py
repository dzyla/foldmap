# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Plan 29's harness: a fresh random sample drawn through the app's own code path (plan 29 Task 1; needs network for
the draw and for downloads; not run in CI). A thin driver over scripts/random100.py: its sampler and its app view.

    python scripts/app29.py draw
    python scripts/app29.py split [--sample docs/validation/app29/sample.json]
    python scripts/app29.py run --out DIR [--ids A,B] [--sample ...] [--jobs 1] [--record RECORD_DIR] [--ui]
    python scripts/app29.py ui ID --out DIR   (run --ui's child: one tune-half id through the app's page)

draw runs random100's sampler with plan 29's call (DRAW: seed 20261027, entries released after 2026-10-13 and on or
before 2026-10-27, the three earlier samples' ids excluded, no AlphaFold DB model: ruling Q61) into
docs/validation/app29/sample.json, then splits it. A pool short of its quota gives what it has (the window is young:
no membrane entry yet, and S5, S8 and S9 are small); the sampler records how short (its quotas), and report.md says it
per stratum.

split adds `"split": {"tune": [...], "check": [...]}` to the sample: split() sorts the ids by
sha256(f"{seed}:{id}") and takes the first round(CHECK_SHARE x n) as check (32 of the 79 drawn: ruling P29-2, the
plan's 60/40), the rest as tune. Fixes are tuned on `tune` and measured
on `check`, which the fix was not tuned on.

run draws each id (default: every entry of the sample) the app's way, one child process per entry (`entry`, run with
python -I, killed after ENTRY_TIMEOUT seconds; its wall time, exit and max RSS from random100.run_child):
random100.app_view, which takes the size check, app_jobs.draw_job under the hosted limits with the app's default
settings, the PNG download (app_jobs.figure_file, files_job) and judges the figure the app drew as `foldmap check`
does (A1-A20 and R1, the readability floor). Each entry writes DIR/<id>/app.json, app.svg, app.png and its child's
stdout.txt and stderr.txt; run writes DIR/report.json and DIR/report-check.json. With --record RECORD_DIR it also
writes report.json, report-check.json and report.md there, stamped with the code that ran (foldmap.source_stamp).
--ui adds, for each tune-half id, the visitor's session (plan 29 Task 2, `ui`, a child process per entry, killed
after UI_TIMEOUT seconds): the app's own page driven headless (streamlit AppTest, through tests/app_helpers.py, as the
app's tests drive it), hosted (FOLDMAP_APP_LOCAL_FILES=0): the id typed and Draw pressed, the theme changed once, then
every download button pressed (the callable a click runs on the server: app_helpers.press). It writes DIR/<id>/ui.json
and the files the buttons gave (ui.svg, ui.png, ui.pdf, ui-layout.yaml, ui-data.json, ui-data.zip); the report's
`ui` per entry is that record, "not run" without --ui and always for the check half (never driven through the page).
Each entry's row lists the sha256 of every figure and file kept in its folder ("files").

The check half is sealed (rulings P29-4, P29-5): report.json, report.md and the console name no check-half id (it is
in the totals only); its rows are in report-check.json alone, read only by load_both (Task N+1's before/after;
load_tune reads report.json alone); a check-half entry's own files go to DIR/sealed/<sha256(f"{seed}:{id}")[:12]>/,
whose app.json names no id (report-check.json is the only map back); and --ids refuses check-half ids unless
--unseal (Task N+1 only).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    # the foldmap of this checkout, not whichever one is installed (children run with python -I: no PYTHONPATH)
    sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from foldmap.source_stamp import code_mix, source_stamp, stamp_line  # noqa: E402

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts
SCRIPT = Path(__file__).resolve()
RANDOM100 = SCRIPT.parent / "random100.py"
OUT = ROOT / "docs" / "validation" / "app29"
SAMPLE = OUT / "sample.json"
SEED = 20261027
CHECK_SHARE = 0.4  # of the sample measured on, never tuned on: the plan's 60/40 (ruling P29-2)
ENTRY_TIMEOUT = 900.0  # s: one entry's process (the draw job's 90 s, the PNG, the redraw and the checks)
UI_TIMEOUT = 1800.0  # s: one entry's page session (two draws and six downloads, each job under the hosted 90 s)
UI_STEP = 600.0  # s: AppTest's wait for one run of the page (a draw job is killed at 90 s)
UI_FILES = {  # download button key -> the file kept of what it gave
    "dl_svg": "ui.svg",
    "dl_png": "ui.png",
    "dl_pdf": "ui.pdf",
    "dl_layout": "ui-layout.yaml",
    "dl_data_json": "ui-data.json",
    "dl_data_zip": "ui-data.zip",
}
UI_MAGIC = {"dl_svg": b"<svg", "dl_png": b"\x89PNG", "dl_pdf": b"%PDF", "dl_data_zip": b"PK"}  # how each file starts
DRAW = [  # plan 29 Task 1 Step 1, as given; --no-alphafold: ruling Q61
    "sample",
    "--seed",
    str(SEED),
    "--date",
    "2026-10-27",
    "--released-after",
    "2026-10-13",
    "--exclude",
    "docs/validation/random100/sample.json",
    "--exclude",
    "docs/validation/random100-heldout/sample.json",
    "--exclude",
    "docs/validation/random100-third/sample.json",
    "--no-alphafold",
    "--out",
    "docs/validation/app29/sample.json",
]
_R100 = None
_R100_LOCK = threading.Lock()


def random100():
    """scripts/random100.py as a module (loaded once, without putting scripts/ on sys.path). Thread-safe: run --jobs
    calls it from several threads at once, and the module is published only once it has run (a thread never sees it
    half loaded, without run_child)."""
    global _R100
    with _R100_LOCK:
        if _R100 is None:
            spec = importlib.util.spec_from_file_location("random100", RANDOM100)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            _R100 = mod
    return _R100


# ------------------------------------------------------------------------------------------------ sample


def split(ids, seed: int = SEED, n_check: int | None = None) -> dict[str, list[str]]:
    """The ids sorted by sha256(f"{seed}:{id}"): the first n_check (default round(CHECK_SHARE x n)) are "check", the
    rest "tune". Whatever order the ids come in, the same split."""
    ranked = sorted(set(ids), key=lambda i: hashlib.sha256(f"{seed}:{i}".encode()).hexdigest())
    n_check = round(CHECK_SHARE * len(ranked)) if n_check is None else n_check
    return {"tune": ranked[n_check:], "check": ranked[:n_check]}


def add_split(path: Path, seed: int = SEED, n_check: int | None = None) -> dict:
    """Write the split of the sample's entries into the sample file (its "split"), with how it was made."""
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    doc["split"] = split([e["id"] for e in doc["entries"]], seed, n_check)
    n = len(doc["split"]["check"])
    doc["split_rule"] = f"sorted by sha256(f'{{seed}}:{{id}}'), seed {seed}: the first {n} check, the rest tune"
    Path(path).write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    return doc


def strata_shortfall(doc: dict) -> list[dict]:
    """Per stratum of a sample: its target, the entries drawn, how many short, and why when it was not drawn."""
    drawn: dict[str, int] = {}
    for e in doc["entries"]:
        drawn[e["stratum"]] = drawn.get(e["stratum"], 0) + 1
    rows = []
    for st in doc["strata"]:
        got = drawn.get(st["name"], 0)
        row = {"stratum": st["name"], "title": st["title"], "target": st["n"], "drawn": got}
        row["short"] = max(st["n"] - got, 0)
        if st.get("not_drawn"):
            row["not_drawn"] = st["not_drawn"]
        rows.append(row)
    return rows


def cmd_draw(args) -> int:
    """random100's sampler with DRAW (a child, as the command line it is), then the split."""
    out = ROOT / DRAW[DRAW.index("--out") + 1]
    t0 = time.time()
    cmd = [sys.executable, str(RANDOM100), *DRAW]
    print(" ".join(cmd[1:]), flush=True)
    code = subprocess.run(cmd, cwd=ROOT).returncode  # 1: short of a quota, which a young window is (written)
    if code not in (0, 1) or not out.is_file() or out.stat().st_mtime < t0:
        print(f"the sampler failed (exit {code}); no sample written", file=sys.stderr)
        return code or 1
    doc = add_split(out)
    _print_draw(doc)
    return 0


def cmd_split(args) -> int:
    _print_draw(add_split(Path(args.sample)))
    return 0


def _print_draw(doc: dict) -> None:
    for r in strata_shortfall(doc):
        print(f"{r['stratum']:<4} {r['target']:>3} {r['drawn']:>3}" + (f"  short {r['short']}" if r["short"] else ""))
    print(f"{len(doc['entries'])} entries; split: tune {len(doc['split']['tune'])}, check {len(doc['split']['check'])}")


# --------------------------------------------------------------------------------------------------- run


SEALED = "sealed"  # out/sealed/<key>/: a check-half entry's own files, named by no id (ruling P29-5)


def sealed_key(i: str, seed: int = SEED) -> str:
    """The name of a check-half entry's folder: sha256(f"{seed}:{id}")[:12]. Only report-check.json maps it back."""
    return hashlib.sha256(f"{seed}:{i}".encode()).hexdigest()[:12]


def _scrub(x, words: list[str], to: str):
    """x with every string's `words` replaced by `to` (a sealed entry's id, its file's path and name)."""
    if isinstance(x, dict):
        return {k: _scrub(v, words, to) for k, v in x.items()}
    if isinstance(x, list):
        return [_scrub(v, words, to) for v in x]
    if isinstance(x, str):
        for w in words:
            x = x.replace(w, to)
    return x


def cmd_entry(args) -> int:
    """One id drawn the app's way (run's child): OUT/<id>/app.json, app.svg, app.png; with --sealed KEY (a check-half
    entry) OUT/sealed/KEY/, where nothing it writes names the id: not the folder, not app.json (the id, the file's
    name and path in a message are replaced), and the app's messages name the entry by KEY."""
    r100 = random100()
    keep = Path(args.out) / SEALED / args.sealed if args.sealed else Path(args.out) / args.id
    keep.mkdir(parents=True, exist_ok=True)
    name = args.sealed or args.id
    rec: dict = {"sealed": args.sealed} if args.sealed else {"id": args.id}
    rec["source_stamp"] = SOURCE_STAMP
    words = [args.id]
    try:
        path = r100.resolve(args.id)
    except Exception as err:  # obsolete or unreachable: no figure, said
        rec["app"] = {"refused": None, "error": f"download failed: {type(err).__name__}: {err}"[:2000]}
    else:
        words = [str(path), path.name, args.id]
        rec["file"] = {"sha256": r100._sha256(path.read_bytes())} | ({} if args.sealed else {"name": path.name})
        rec["app"] = r100.app_view(path, name, keep=keep, judge=True)
    rec = r100._jsonable(rec)
    if args.sealed:
        rec = _scrub(json.loads(json.dumps(rec, default=str)), words, f"<{args.sealed}>")
    text = json.dumps(rec, indent=1, default=str) + "\n"
    (keep / "app.json").write_text(text, encoding="utf-8", newline="\n")
    return 0


def run_one(i: str, out: Path, key: str | None = None) -> dict:
    """The entry's child process, then what it wrote (or why it wrote nothing). key: a check-half entry's sealed
    name; its folder, app.json and the child's stdout.txt and stderr.txt are under out/sealed/<key>/."""
    r100 = random100()
    keep = out / SEALED / key if key else out / i
    keep.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-I", "-W", "default", str(SCRIPT), "entry", i, "--out", str(out)]
    cmd += ["--sealed", key] if key else []
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "PYTHONHASHSEED": "0"}
    (keep / "app.json").unlink(missing_ok=True)
    proc = r100.run_child(cmd, env=env, timeout=ENTRY_TIMEOUT, out=keep / "stdout.txt", err=keep / "stderr.txt")
    got = keep / "app.json"
    if got.is_file():
        rec = json.loads(got.read_text(encoding="utf-8"))
    else:  # killed or crashed: this process says so, under its own stamp (it wrote this verdict)
        err = keep / "stderr.txt"
        text = err.read_text(encoding="utf-8", errors="replace") if err.is_file() else ""
        tail = text.strip().splitlines()[-1:] or [""]
        why = "timed out" if proc["timed_out"] else f"exit {proc['exit']}: {tail[0][:300]}"
        rec = {"source_stamp": SOURCE_STAMP, "app": {"refused": None, "error": f"no result ({why})"}}
    rec["id"] = i  # in this process only: a sealed row reaches report-check.json alone
    rec["dir"] = f"{SEALED}/{key}" if key else i  # the entry's folder under out
    rec["process"] = proc
    return rec


def _app_helpers():
    """tests/app_helpers.py as a module: the page driven as the app's tests drive it."""
    spec = importlib.util.spec_from_file_location("app_helpers_app29", ROOT / "tests" / "app_helpers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ui_session(i: str, keep: Path) -> dict:
    """One visitor's session on the app's page (hosted): type `i`, press Draw; change the theme once; press every
    download button. {"exception": the page's uncaught exception (str) or None, "errors": [{"step", "text"}] (the
    page's st.error), "warnings": the same of st.warning, "steps": what was done, "stage": the page's last,
    "downloads": {key: {"bytes", "sha256", "ok"} or {"error"}}, "harness_error": what stopped the session outside
    the page (an AppTest timeout), else None}. Each download's file is kept in `keep` (UI_FILES)."""
    import pytest

    helpers = _app_helpers()
    rec: dict = {"exception": None, "errors": [], "warnings": [], "steps": [], "stage": None, "downloads": {}}
    rec["harness_error"] = None

    def seen(at, step: str) -> bool:
        """The page after `step`: its errors and warnings kept; False when it raised (nothing more is done)."""
        rec["steps"].append(step)
        rec["stage"] = at.session_state["stage"] if "stage" in at.session_state else None
        rec["errors"] += [{"step": step, "text": str(e.value)} for e in at.error]
        rec["warnings"] += [{"step": step, "text": str(w.value)} for w in at.warning]
        if at.exception:
            rec["exception"] = "\n".join(str(e.value) for e in at.exception)[:4000]
            return False
        return True

    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("FOLDMAP_APP_LOCAL_FILES", "0")  # the public app's mode: hosted limits, jobs in their own process
        kept = helpers.keep_downloads(mp)
        try:
            at = helpers.load(helpers.app_test(UI_STEP), i)
            if not seen(at, "load") or rec["stage"] != "edit":
                return rec  # not drawn: what the page said is the record
            box = at.selectbox(key="theme")
            other = next(t for t in box.options if t != box.value)
            box.select(other).run()
            if not seen(at, f"theme {other}"):
                return rec
            for button in at.download_button:
                key = button.key
                try:
                    data = helpers.press(at, key, kept)
                except Exception as err:  # the button gave no file: said, the others still pressed
                    rec["downloads"][key] = {"error": f"{type(err).__name__}: {err}"[:2000]}
                    continue
                name = UI_FILES.get(key, f"ui-{key}")
                (keep / name).write_bytes(data)
                ok = bool(data) and data.startswith(UI_MAGIC.get(key, b""))
                rec["downloads"][key] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "ok": ok}
            rec["steps"].append("downloads")
        except Exception as err:  # outside the page: AppTest's own timeout, a widget not found
            rec["harness_error"] = f"{type(err).__name__}: {err}"[:2000]
    return rec


def cmd_ui(args) -> int:
    """One tune-half id through the app's page (run --ui's child): OUT/<id>/ui.json and the files it downloaded. Never
    a check-half id (its entries are never driven through the page: ruling P29-5)."""
    sample = Path(args.sample)
    doc = json.loads(sample.read_text(encoding="utf-8")) if sample.is_file() else {}
    if args.id in ((doc.get("split") or {}).get("check") or []):
        raise SystemExit("ui: that id is in the check half, which is never driven through the page")
    keep = Path(args.out) / args.id
    keep.mkdir(parents=True, exist_ok=True)
    rec = {"id": args.id, "source_stamp": SOURCE_STAMP, "ui": ui_session(args.id, keep)}
    text = json.dumps(rec, indent=1, default=str) + "\n"
    (keep / "ui.json").write_text(text, encoding="utf-8", newline="\n")
    return 0


def run_ui(i: str, out: Path, sample: Path) -> tuple[dict, dict | None]:
    """The entry's page session in its own process (cmd_ui): (its ui record, the stamp it ran under, or None when
    run wrote the verdict itself)."""
    r100 = random100()
    keep = out / i
    keep.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-I", "-W", "default", str(SCRIPT), "ui", i, "--out", str(out), "--sample", str(sample)]
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "PYTHONHASHSEED": "0", "FOLDMAP_APP_LOCAL_FILES": "0"}
    got = keep / "ui.json"
    got.unlink(missing_ok=True)
    proc = r100.run_child(cmd, env=env, timeout=UI_TIMEOUT, out=keep / "ui-stdout.txt", err=keep / "ui-stderr.txt")
    if got.is_file():
        rec = json.loads(got.read_text(encoding="utf-8"))
        return {**rec["ui"], "process": proc}, rec.get("source_stamp")
    err = keep / "ui-stderr.txt"
    text = err.read_text(encoding="utf-8", errors="replace") if err.is_file() else ""
    tail = text.strip().splitlines()[-1:] or [""]
    why = "timed out" if proc["timed_out"] else f"exit {proc['exit']}: {tail[0][:300]}"
    return {"exception": None, "harness_error": f"no result ({why})", "process": proc}, None


def kept_files(folder: Path) -> dict[str, str]:
    """name -> sha256 of the figures and files an entry kept (app.svg, app.png, the page session's ui files)."""
    if not folder.is_dir():
        return {}
    names = sorted(f.name for f in folder.iterdir() if f.is_file() and f.name.startswith(("app.", "ui")))
    skip = {"app.json", "ui.json", "ui-stdout.txt", "ui-stderr.txt"}
    return {n: hashlib.sha256((folder / n).read_bytes()).hexdigest() for n in names if n not in skip}


def ui_failed(ui) -> bool:
    """A page session that raised, stopped, or whose download gave no whole file (an error shown is not a failure
    by itself: a refusal is said with one)."""
    if not isinstance(ui, dict):
        return False
    bad = any(not d.get("ok") for d in (ui.get("downloads") or {}).values())
    return bool(ui.get("exception") or ui.get("harness_error") or bad)


def _last_line(text: str | None) -> str | None:
    """A traceback's last line (what was raised), None without one."""
    lines = (text or "").strip().splitlines()
    return lines[-1] if lines else None


def row(rec: dict, half: str, stratum: str) -> dict:
    """An entry as the report lists it."""
    a = rec.get("app") or {}
    crit = a.get("criteria") or {}
    floor = a.get("floor")
    return {
        "id": rec["id"],
        "dir": rec.get("dir"),
        "stratum": stratum,
        "half": half,
        "residues": a.get("residues"),
        "refused": a.get("refused"),
        "error": a.get("error"),
        "measure_error": _last_line(a.get("measure_error")),
        "png_error": a.get("png_error"),
        "criteria_failed": sorted(k for k, v in crit.items() if v == "fail"),
        "criteria_detail": {k: v for k, v in (a.get("criteria_detail") or {}).items() if crit.get(k) == "fail"},
        "floor_ok": None if floor is None else floor["ok"],
        "floor_reasons": [] if floor is None else floor["reasons"],
        "seconds": a.get("seconds"),
        "draw_seconds": a.get("draw_seconds"),
        "peak_mb": a.get("peak_mb", "not measured"),
        "svg_bytes": a.get("svg_bytes"),
        "png_bytes": a.get("png_bytes"),
        "process": rec.get("process"),
        "files": rec.get("files", {}),
        "ui": rec.get("ui", "not run"),
    }


def verdict(r: dict) -> str:
    """drawn (a figure judged), refused (a limit), error, not measured (drawn, but the judging failed)."""
    if r["error"]:
        return "error"
    if r["refused"]:
        return f"refused ({r['refused']})"
    if r["measure_error"]:
        return "not measured"
    return "drawn"


def build_report(recs: list[dict], doc: dict | None, sample_path: Path | None) -> dict:
    halves = {i: h for h, ids in ((doc or {}).get("split") or {}).items() for i in ids}
    strata = {e["id"]: e["stratum"] for e in (doc or {}).get("entries", [])}
    rows = [row(r, halves.get(r["id"], "-"), strata.get(r["id"], "-")) for r in recs]
    totals = {}
    for half in ("tune", "check", "all"):
        mine = [r for r in rows if half == "all" or r["half"] == half]
        drawn = [r for r in mine if verdict(r) == "drawn"]
        totals[half] = {
            "entries": len(mine),
            "drawn": len(drawn),
            "refused": sum(1 for r in mine if r["refused"]),
            "error": sum(1 for r in mine if r["error"]),
            "not_measured": sum(1 for r in mine if verdict(r) == "not measured"),
            "png_error": sum(1 for r in mine if r["png_error"]),
            "criteria_failed": sum(1 for r in drawn if r["criteria_failed"]),
            "floor_failed": sum(1 for r in drawn if r["floor_ok"] is False),
            "ui_run": sum(1 for r in mine if isinstance(r["ui"], dict)),
            "ui_failed": sum(1 for r in mine if ui_failed(r["ui"])),
        }
    out = {
        "source_stamp": SOURCE_STAMP,
        **code_mix([r.get("source_stamp") for r in recs] + [r["ui_stamp"] for r in recs if r.get("ui_stamp")]),
        "sample": None,
        "totals": totals,
        "entries": rows,
    }
    if doc is not None:
        out["sample"] = {
            "path": sample_path.resolve().relative_to(ROOT).as_posix()
            if sample_path.resolve().is_relative_to(ROOT)
            else str(sample_path),
            "sha256": hashlib.sha256(sample_path.read_bytes()).hexdigest(),
            "seed": doc.get("seed"),
            "query_date": doc.get("query_date"),
            "released_after": doc.get("released_after"),
            "split": {h: len(ids) for h, ids in (doc.get("split") or {}).items()},
            "strata": strata_shortfall(doc),
        }
    return out


def render(rep: dict) -> str:
    t = rep["totals"]
    lines = [
        "# The app on a fresh random sample",
        "",
        f"Plan 29. Written by: {stamp_line(rep['source_stamp'])}",
        "",
        "Every entry drawn the app's way (scripts/app29.py run: random100.app_view, the size check and "
        "app_jobs.draw_job under the hosted limits, the app's default settings, the PNG download) and the figure "
        "judged as `foldmap check` judges it: A1-A20 and R1, the readability floor.",
        "",
    ]
    if rep.get("mixed_code"):
        lines += ["**Mixed code: the entries were measured by different code; never quote this report.**", ""]
    s = rep.get("sample")
    if s:
        lines += [
            f"Sample: `{s['path']}` (sha256 {s['sha256'][:12]}), seed {s['seed']}, released after "
            f"{s['released_after']} and on or before {s['query_date']}; split: "
            + ", ".join(f"{h} {n}" for h, n in s["split"].items())
            + ".",
            "",
            "## The sample per stratum",
            "",
            "| stratum | title | target | drawn | short |",
            "|---|---|---|---|---|",
        ]
        for r in s["strata"]:
            why = f" ({r['not_drawn']})" if r.get("not_drawn") else ""
            lines.append(f"| {r['stratum']} | {r['title']} | {r['target']} | {r['drawn']} | {r['short']}{why} |")
        lines.append("")
    lines += [
        "## Totals",
        "",
        "| half | entries | drawn | refused | error | not measured | PNG error | criteria failed | floor failed "
        "| page run | page failed |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for half, n in t.items():
        lines.append(
            f"| {half} | {n['entries']} | {n['drawn']} | {n['refused']} | {n['error']} | {n['not_measured']} | "
            f"{n['png_error']} | {n['criteria_failed']} | {n['floor_failed']} | {n.get('ui_run', 0)} | "
            f"{n.get('ui_failed', 0)} |"
        )
    lines += [
        "",
        "Peak memory of the draw job: not measured unless app_jobs hands it back; `process` in report.json is the "
        "whole entry process's (its wall time and max RSS include the PNG, the redraw and the checks).",
        "",
        "## Entries",
        "",
        "The check half is reported in the totals only: its ids are never shown beside an outcome here, nor on the "
        "console, so no fix is tuned on them.",
        "",
        "| id | stratum | half | residues | outcome | criteria failed | floor | seconds | draw s |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    shown = [r for r in rep["entries"] if r["half"] != "check"]
    for r in shown:
        floor = "-" if r["floor_ok"] is None else "pass" if r["floor_ok"] else "fail"
        lines.append(
            f"| {r['id']} | {r['stratum']} | {r['half']} | {r['residues']} | {verdict(r)} | "
            f"{', '.join(r['criteria_failed']) or '-'} | {floor} | {r['seconds']} | {r.get('draw_seconds')} |"
        )
    said = [r for r in shown if r["error"] or r["refused"] or r["measure_error"] or r["png_error"]]
    if said:
        lines += ["", "## Refusals and errors", ""]
        for r in said:
            what = r["error"] or r["measure_error"] or r["png_error"] or f"refused: {r['refused']}"
            lines.append(f"- {r['id']}: {str(what).splitlines()[0][:300]}")
    ran = [r for r in shown if isinstance(r["ui"], dict)]
    if ran:
        lines += [
            "",
            "## The page (run --ui)",
            "",
            "Each tune-half id typed into the hosted app's page, Draw pressed, the theme changed once and every "
            "download button pressed (AppTest; the check half is never driven through the page). Page failed: it "
            "raised, the session stopped, or a download gave no whole file.",
            "",
            "| id | stage | steps | exception | errors shown | downloads whole | failed |",
            "|---|---|---|---|---|---|---|",
        ]
        for r in ran:
            u = r["ui"]
            downs = u.get("downloads") or {}
            whole = f"{sum(1 for d in downs.values() if d.get('ok'))}/{len(downs)}"
            exc = _last_line(u.get("exception")) or _last_line(u.get("harness_error")) or "-"
            lines.append(
                f"| {r['id']} | {u.get('stage')} | {', '.join(u.get('steps') or []) or '-'} | {exc[:200]} | "
                f"{len(u.get('errors') or [])} | {whole} | {'yes' if ui_failed(u) else '-'} |"
            )
    return "\n".join(lines) + "\n"


CHECK_FILE = "report-check.json"  # the check half's rows, sealed apart from report.json (ruling P29-4)


def write_report(rep: dict, folder: Path, md: bool) -> None:
    """report.json (and with `md` report.md): the tune half per id, the check half in the totals only; the check
    half's rows go to report-check.json alone (ruling P29-4), read back only by load_both."""
    folder.mkdir(parents=True, exist_ok=True)
    shown = {**rep, "source_stamp": SOURCE_STAMP, "entries": [r for r in rep["entries"] if r["half"] != "check"]}
    sealed = {
        "source_stamp": SOURCE_STAMP,
        "what": "the check half's rows of report.json: for plan 29's before/after measure only (load_both)",
        "entries": [r for r in rep["entries"] if r["half"] == "check"],
    }
    (folder / "report.json").write_text(json.dumps(shown, indent=1) + "\n", encoding="utf-8", newline="\n")
    (folder / CHECK_FILE).write_text(json.dumps(sealed, indent=1) + "\n", encoding="utf-8", newline="\n")
    if md:
        (folder / "report.md").write_text(render(shown), encoding="utf-8", newline="\n")


def load_tune(record_dir) -> dict:
    """A run's report as it may be read while fixing: report.json only (the tune half per id, the check half in the
    totals)."""
    return json.loads((Path(record_dir) / "report.json").read_text(encoding="utf-8"))


def load_both(record_dir, out=None) -> dict:
    """For Task N+1 only (plan 29's before/after measure, once the fixes are done): report.json joined with the check
    half's rows from report-check.json. Never read while a fix is being tuned. Each row's "dir" is its folder under
    the run's out folder (a check row's is sealed/<key>: report-check.json is the only map from key to id); with
    `out`, "folder" is that path."""
    rep = load_tune(record_dir)
    sealed = json.loads((Path(record_dir) / CHECK_FILE).read_text(encoding="utf-8"))
    entries = rep["entries"] + sealed["entries"]
    if out is not None:
        entries = [{**r, "folder": str(Path(out) / r["dir"]) if r.get("dir") else None} for r in entries]
    return {**rep, "entries": entries, "check_source_stamp": sealed["source_stamp"]}


def cmd_run(args) -> int:
    out = Path(args.out)
    sample_path = Path(args.sample)
    doc = json.loads(sample_path.read_text(encoding="utf-8")) if sample_path.is_file() else None
    if args.ids:
        ids = [i.strip() for i in args.ids.split(",") if i.strip()]
    elif doc is not None:
        ids = [e["id"] for e in doc["entries"]]
    else:
        raise SystemExit(f"no sample at {sample_path}; give --ids")
    out.mkdir(parents=True, exist_ok=True)

    halves = {i: h for h, part in ((doc or {}).get("split") or {}).items() for i in part}
    seed = (doc or {}).get("seed", SEED)
    asked = sum(1 for i in ids if halves.get(i) == "check")
    if args.ids and asked and not args.unseal:  # the check half is run whole, sealed; never picked by id (P29-5)
        raise SystemExit(
            f"--ids names {asked} id(s) of the check half; they are measured only in a whole run, sealed, or with "
            "--unseal (Task N+1 only)"
        )

    def one(i: str) -> dict:
        check = halves.get(i) == "check"
        rec = run_one(i, out, sealed_key(i, seed) if check else None)
        if args.ui and not check:  # the page session: the tune half only (P29-5: the check half never)
            rec["ui"], rec["ui_stamp"] = run_ui(i, out, sample_path)
        rec["files"] = kept_files(out / rec["dir"]) if rec.get("dir") else {}
        if not check:  # the check half is counted, never named beside an outcome
            ui = f"; page {'failed' if ui_failed(rec['ui']) else 'ok'}" if isinstance(rec.get("ui"), dict) else ""
            print(f"{i}: {verdict(row(rec, '-', '-'))}{ui}", flush=True)
        return rec

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        recs = list(pool.map(one, ids))
    rep = build_report(recs, doc, sample_path if doc is not None else None)
    write_report(rep, out, md=False)
    if args.record:
        write_report(rep, Path(args.record), md=True)
    t = rep["totals"]["all"]
    print(
        f"{t['entries']} entries: {t['drawn']} drawn, {t['refused']} refused, {t['error']} errors, "
        f"{t['not_measured']} not measured; {t['criteria_failed']} with a criterion failed, "
        f"{t['floor_failed']} under the floor"
    )
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("draw", help="draw plan 29's sample (network) and split it")
    s = sub.add_parser("split", help="write the tune/check split into a sample")
    s.add_argument("--sample", default=str(SAMPLE))
    r = sub.add_parser("run", help="draw each entry the app's way and report")
    r.add_argument("--out", required=True, help="where each entry's app.json, app.svg, app.png and report.json go")
    r.add_argument("--ids", help="comma-separated ids (default: every entry of the sample)")
    r.add_argument("--sample", default=str(SAMPLE))
    r.add_argument("--jobs", type=int, default=1, help="entries drawn at once")
    r.add_argument("--record", help="also write a stamped report.json and report.md into this folder")
    r.add_argument("--ui", action="store_true", help="also drive the app headless (AppTest; plan 29 Task 2)")
    r.add_argument(
        "--unseal",
        action="store_true",
        help="let --ids name check-half ids: for Task N+1 (the before/after measure) only, never while tuning a fix",
    )
    u = sub.add_parser("ui", help="one tune-half id through the app's page, headless (run --ui's child)")
    u.add_argument("id")
    u.add_argument("--out", required=True)
    u.add_argument("--sample", default=str(SAMPLE))
    e = sub.add_parser("entry", help="one id drawn the app's way (run's child)")
    e.add_argument("id")
    e.add_argument("--out", required=True)
    e.add_argument("--sealed", metavar="KEY", help="a check-half entry: write under OUT/sealed/KEY, naming no id")
    args = parser.parse_args(argv)
    return {"draw": cmd_draw, "split": cmd_split, "run": cmd_run, "entry": cmd_entry, "ui": cmd_ui}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
