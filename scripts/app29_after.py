# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Plan 29 Task N+1: the after measurement's record and the look test's sheets, written by script and stamped with
the code that wrote them (foldmap.source_stamp).

    python scripts/app29_after.py report --before REC --rerun REC --after REC [--measure FILE]... --out FILE.md
    python scripts/app29_after.py sheets --rerun REC --rerun-out DIR --after REC --after-out DIR
        --leads-before DIR --leads-after DIR --out DIR --look FILE.md

report reads three scripts/app29.py run records whole (app29.load_both: the check half is opened here, once; ruling
P29-5): the committed before record (zyla-office), the before commit rerun on the after run's host (the seconds are
comparable only on one host), and the after record. Per half (tune, check, all) it writes the pass count, crashes,
refusals, floor exits and failed criteria, each with a Wilson 95% interval (scripts/paper/common.wilson), and the
median and p95 of the draw job's wall time ("seconds": what a visitor waits) and of draw_job's own time
("draw_seconds"), each with a distribution-free 95% interval: the order statistics at the ranks n x the Wilson
interval of the quantile's share. It names every entry whose outcome or figure changed between the rerun and the
after run, and every entry where the committed before record and the rerun differ (the host). Each --measure (a JSON
file of scripts/app29_light.py: campaign, navalidate, science) is summarised, every changed entry with its cause
(CAUSES; "NOT EXPLAINED" when none is written). It writes FILE.md and FILE.json beside it.

sheets draws, for every entry whose app.svg changed between the rerun and the after run (tune and check), for the
leads (LEADS) and for the named crops (CROPS), the app's PNG before beside after at one scale (pixels per point of
the page, the smaller of the two PNGs' and of what fits SHEET_PX), into --out, and writes --look: one row per pair
with the sheet's sha256 and both SVGs', and an empty verdict column that the controller fills by hand after looking
(plan 29 Task N+1 Step 3).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HALVES = ("tune", "check", "all")
LEADS = {  # the look test's leads (plan 29 Task N+1, the controller's additions), with why each is looked at
    "31FD": "Task 5 named entry: a two-copy complex under the floor, now one subunit",
    "6JLN": "Task 5 lead: a two-copy complex under the floor (plan 28's six sample defaults)",
    "9FZL": "Task 5 lead: under the floor (plan 28's six sample defaults; third sample)",
    "9X4L": "Task 4 named entry (A20) and Task 7 (S3); a giant assembly",
    "26EL": "Task 8b named entry: focus auto draws one copy (A alone)",
    "9WNE": "Task 4: the Mg moved off the loop crossing",
    "7SPI": "Task 8b: TraV copies modelled in part, one protein (third sample)",
    "1KX5": "Task 4: the nucleosome's Mn mark moved (NA control)",
}
CROPS = {  # entry -> (the SVG ids cropped to, why): the look test's named rows
    "9WNE": (
        r"ligand(-tether|-label)?:MG:A13",
        "P29-10: the A1301/A1302 Mg tethers lying along each other (Task 4), cropped to the two marks",
    ),
}
SHEET_PX = 1600  # a panel's largest side, px: the sheet's scale is the smaller of this and the PNGs' own
CAUSES: dict[str, str] = {}  # entry -> the cause of its change in the third-sample and NA measures (after the run)


def _module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _app29():
    return _module(ROOT / "scripts" / "app29.py", "app29_after_app29")


def _common():
    return _module(ROOT / "scripts" / "paper" / "common.py", "app29_after_common")


def _stamp() -> dict:
    if str(ROOT / "src") not in sys.path and (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    from foldmap.source_stamp import source_stamp

    return source_stamp(__file__)


# ------------------------------------------------------------------------------------------------ the table


def passed(r: dict) -> bool:
    """Drawn, judged, no criterion failed and above the readability floor."""
    return not (r["error"] or r["refused"] or r["measure_error"]) and not r["criteria_failed"] and r["floor_ok"] is True


def drawn(r: dict) -> bool:
    return not (r["error"] or r["refused"] or r["measure_error"])


def quantile_interval(xs: list[float], q: float, wilson) -> dict:
    """The q-quantile of xs (common.percentile's interpolation) and a distribution-free 95% interval: the order
    statistics at ranks floor(n x lo) and ceil(n x hi), (lo, hi) the Wilson interval of the share q of n, clamped to
    1..n."""
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return {"n": 0}
    k = (n - 1) * q
    lo_i, hi_i = math.floor(k), math.ceil(k)
    value = xs[lo_i] + (xs[hi_i] - xs[lo_i]) * (k - lo_i)
    pl, ph = wilson(q * n, n)
    a, b = max(1, math.floor(n * pl)), min(n, math.ceil(n * ph))
    return {"n": n, "value": round(value, 3), "ci": [xs[a - 1], xs[b - 1]]}


def half_stats(rows: list[dict], wilson) -> dict:
    n = len(rows)

    def count(k: int) -> dict:
        ci = wilson(k, n)
        return {"k": k, "n": n, "ci": None if ci is None else [round(ci[0], 4), round(ci[1], 4)]}

    ok = [r for r in rows if drawn(r)]
    out = {
        "entries": n,
        "pass": count(sum(passed(r) for r in rows)),
        "crash": count(sum(1 for r in rows if r["error"])),
        "refused": count(sum(1 for r in rows if r["refused"])),
        "floor": count(sum(1 for r in ok if r["floor_ok"] is False)),
        "criteria": count(sum(1 for r in ok if r["criteria_failed"])),
        "page_failed": count(sum(1 for r in rows if _ui_failed(r.get("ui")))),
    }
    for key in ("seconds", "draw_seconds"):
        xs = [float(r[key]) for r in ok if isinstance(r.get(key), (int, float))]
        out[key] = {"median": quantile_interval(xs, 0.5, wilson), "p95": quantile_interval(xs, 0.95, wilson)}
    return out


def _ui_failed(ui) -> bool:
    if not isinstance(ui, dict):
        return False
    bad = any(not d.get("ok") for d in (ui.get("downloads") or {}).values())
    return bool(ui.get("exception") or ui.get("harness_error") or bad)


def outcome(r: dict) -> str:
    if r["error"]:
        return "error"
    if r["refused"]:
        return "refused"
    if r["measure_error"]:
        return "not measured"
    bits = list(r["criteria_failed"]) + ([] if r["floor_ok"] else ["floor"])
    return "pass" if not bits else "fail " + ",".join(bits)


def changes(a: dict[str, dict], b: dict[str, dict]) -> list[dict]:
    """Entries of both runs whose outcome, SVG or PNG differ (a -> b)."""
    out = []
    for i in sorted(a.keys() & b.keys()):
        x, y = a[i], b[i]
        moved = [f for f in ("app.svg", "app.png") if (x.get("files") or {}).get(f) != (y.get("files") or {}).get(f)]
        if outcome(x) != outcome(y) or moved:
            out.append(
                {
                    "id": i,
                    "half": y.get("half"),
                    "outcome": [outcome(x), outcome(y)],
                    "moved": moved,
                    "seconds": [x.get("seconds"), y.get("seconds")],
                    "ui_failed": [_ui_failed(x.get("ui")), _ui_failed(y.get("ui"))],
                }
            )
    return out


def _rows(rec: dict) -> dict[str, dict]:
    return {r["id"]: r for r in rec["entries"]}


def table(before: dict, rerun: dict, after: dict) -> dict:
    common = _common()
    runs = {"before": _rows(before), "rerun": _rows(rerun), "after": _rows(after)}
    stats = {
        name: {
            h: half_stats([r for r in rows.values() if h == "all" or r.get("half") == h], common.wilson) for h in HALVES
        }
        for name, rows in runs.items()
    }
    return {
        "stats": stats,
        "changed": changes(runs["rerun"], runs["after"]),
        "host_differs": changes(runs["before"], runs["rerun"]),
        "missing": {k: sorted(set(runs["before"]) ^ set(v)) for k, v in runs.items() if k != "before"},
        "stamps": {
            name: {
                "git": (rec.get("source_stamp") or {}).get("git"),
                "check": ((rec.get("check_source_stamp") or {}).get("git") or {}).get("commit"),
                "platform": ((rec.get("source_stamp") or {}).get("platform") or {}).get("platform"),
                "mixed_code": rec.get("mixed_code"),
            }
            for name, rec in (("before", before), ("rerun", rerun), ("after", after))
        },
    }


def _c(x: dict) -> str:
    if x["ci"] is None:
        return f"{x['k']}/{x['n']}"
    return f"{x['k']}/{x['n']} [{x['ci'][0] * 100:.0f}-{x['ci'][1] * 100:.0f}%]"


def _q(x: dict) -> str:
    if not x.get("n"):
        return "-"
    return f"{x['value']:g} [{x['ci'][0]:g}-{x['ci'][1]:g}]"


def table_lines(t: dict) -> list[str]:
    s = t["stats"]
    names = {
        "before": "before (committed, zyla-office)",
        "rerun": "before (rerun, gpu-desktop)",
        "after": "after (gpu-desktop)",
    }
    lines = []
    for h in HALVES:
        lines += [
            f"### {h} half" if h != "all" else "### both halves",
            "",
            "| run | entries | pass | crash | refused | floor exit | criteria failed | page failed | seconds median "
            "| seconds p95 | draw s median | draw s p95 |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        for run in ("before", "rerun", "after"):
            x = s[run][h]
            lines.append(
                f"| {names[run]} | {x['entries']} | {_c(x['pass'])} | {_c(x['crash'])} | {_c(x['refused'])} | "
                f"{_c(x['floor'])} | {_c(x['criteria'])} | {_c(x['page_failed'])} | "
                f"{_q(x['seconds']['median'])} | {_q(x['seconds']['p95'])} | {_q(x['draw_seconds']['median'])} | "
                f"{_q(x['draw_seconds']['p95'])} |"
            )
        lines.append("")
    return lines


# ------------------------------------------------------------------------------------------------ measures


def _cause(i: str) -> str:
    return CAUSES.get(i, "NOT EXPLAINED")


def measure_lines(m: dict) -> list[str]:
    s, kind = m["summary"], m["kind"]
    if kind == "science":
        src = ", ".join(f"{k} {v}" for k, v in s.get("sources", {}).items())
        lines = [
            "### The science measure (`scripts/app29_light.py science --cache`, figcheck at 3d6c798 against now)",
            "",
            f"- {s['judged']} structures judged of {s['entries']} ({src}), {len(s['errors'])} not drawn",
            f"- science findings changed (Task 7's barrel_finding against the old count) on: "
            f"{', '.join(s['differ']) or 'none'}; new findings: {', '.join(s['new_findings']) or 'none'}",
        ]
        for name, ch in s["changes"].items():
            for k, (a, b) in ch.items():
                lines.append(f"  - {name} {k}: {a or 'none'} -> {b or 'none'}")
        lines += [
            f"- barrel blocks holding a ghost lent by another copy: {s.get('lent_blocks', 0)}, in "
            f"{len(s.get('lent_entries', []))} structures: {', '.join(s.get('lent_entries', [])) or 'none'}",
            "- a ring strand lent across the cut (science.lent_ring takes a lent ghost into the ring): "
            + ("; ".join(f"{n}: {', '.join(v)}" for n, v in s.get("ring_lent", {}).items()) or "none found"),
            f"- the probe failed on: {', '.join(s.get('lent_errors', [])) or 'none'}",
        ]
        if s["errors"]:
            lines.append(f"- not drawn: {'; '.join(s['errors'])}")
        return lines + [""]
    st = [((x or {}).get("commit") or "committed")[:7] for x in s["stamps"]]
    lines = [f"### {s['label']}", "", f"- before git {st[0]}, after git {st[1]}; entries {s['entries'][0]} -> "]
    lines[-1] += f"{s['entries'][1]}; pass {s['passed'][0]} -> {s['passed'][1]}"
    if kind == "campaign":
        lines[-1] += f"; keys: {s['verdict']}"
        lines += [
            f"- worse: {'; '.join(s['worse']) or 'none'}",
            f"- better: {'; '.join(s['better']) or 'none'}",
            f"- moved (no judge orders it): {'; '.join(s['moved']) or 'none'}",
            f"- the SVG moved on: {', '.join(s['svg_moved']) or 'none'}",
        ]
        if "na" in s:
            n = s["na"]
            lines.append(
                f"- A17 measures ({n['entries'][1]} entries): changed on {', '.join(n['changed']) or 'none'}; worse: "
                f"{'; '.join(f'{i} {v}' for i, v in n['worse'].items()) or 'none'}; better: "
                f"{'; '.join(f'{i} {v}' for i, v in n['better'].items()) or 'none'}"
            )
    else:
        lines.append(
            f"- measures changed on: {', '.join(s['differ']) or 'none'}; worse: "
            f"{'; '.join(f'{i} {v}' for i, v in s['worse'].items()) or 'none'}; better: "
            f"{'; '.join(f'{i} {v}' for i, v in s['better'].items()) or 'none'}"
        )
    if s.get("reference"):
        lines.append(f"- for reference, not judged: {len(s['differ'])} entries differ (older code or another view)")
        return lines + [""]
    lines += [f"- **{i}**: {_cause(i)}" for i in s["differ"]]
    return lines + [""]


# ------------------------------------------------------------------------------------------------ report


def cmd_report(args) -> int:
    app29 = _app29()
    before, rerun, after = (app29.load_both(d) for d in (args.before, args.rerun, args.after))
    t = table(before, rerun, after)
    measures = [json.loads(Path(f).read_text(encoding="utf-8")) for f in args.measure or []]
    doc = {
        "source_stamp": _stamp(),
        "records": {"before": args.before, "rerun": args.rerun, "after": args.after},
        "table": t,
        "measures": [{"kind": m["kind"], "summary": m["summary"]} for m in measures],
        "causes": CAUSES,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    out.write_text(render(doc, measures), encoding="utf-8", newline="\n")
    print(f"wrote {out} and {out.with_suffix('.json')}")
    return 0


def _stamp_line(st: dict) -> str:
    g = st.get("git") or {}
    return f"git {(g.get('commit') or '?')[:7]} ({'dirty' if g.get('dirty') else 'clean'})"


def render(doc: dict, measures: list[dict]) -> str:
    t = doc["table"]
    st = t["stamps"]
    lines = [
        "# Plan 29: before and after, both halves",
        "",
        f"Written by scripts/app29_after.py report at {_stamp_line(doc['source_stamp'])}. The check half is opened "
        "here, once (ruling P29-5): it was measured, never tuned on; nothing was changed because of it.",
        "",
        "Runs (scripts/app29.py run --ui, the app's code path; the page session on the tune half only):",
        "",
    ]
    for name, what in (
        ("before", f"`{doc['records']['before']}`: the committed before record"),
        ("rerun", f"`{doc['records']['rerun']}`: the before commit rerun on the after run's host, same settings"),
        ("after", f"`{doc['records']['after']}`: the after record"),
    ):
        s = st[name]
        lines.append(
            f"- {name}: {what}; git {((s['git'] or {}).get('commit') or '?')[:7]} (check rows "
            f"{(s['check'] or '?')[:7]}), {s['platform']}{'; MIXED CODE' if s['mixed_code'] else ''}"
        )
    lines += [
        "",
        "Pass: drawn, judged, no criterion failed (A1-A20) and above the readability floor (R1). Crash: no figure "
        "(an error, a job killed). Floor exit and criteria failed: of the drawn. Page failed: the tune half's page "
        "session raised, stopped or a download gave no whole file. Counts carry a Wilson 95% interval of their share; "
        "seconds (the draw job's wall time, what a visitor waits; draw s: draw_job's own time) are of the drawn "
        "entries, each with a distribution-free 95% interval (the order statistics at n x the Wilson interval of the "
        "quantile's share). Seconds compare only within one host: zyla-office's cores are about twice as fast as "
        "gpu-desktop's, and the runs on gpu-desktop drew 96 entries at once.",
        "",
        "## The table",
        "",
        *table_lines(t),
        "## What changed (rerun -> after, same host)",
        "",
    ]
    if t["changed"]:
        lines += ["| id | half | outcome | moved | seconds | page failed |", "|---|---|---|---|---|---|"]
        for c in t["changed"]:
            lines.append(
                f"| {c['id']} | {c['half']} | {c['outcome'][0]} -> {c['outcome'][1]} | {', '.join(c['moved']) or '-'} "
                f"| {c['seconds'][0]} -> {c['seconds'][1]} | {c['ui_failed'][0]} -> {c['ui_failed'][1]} |"
            )
    else:
        lines.append("Nothing.")
    lines += ["", "## The committed before record against its rerun (another host)", ""]
    if t["host_differs"]:
        lines += ["| id | half | outcome | moved |", "|---|---|---|---|"]
        for c in t["host_differs"]:
            lines.append(
                f"| {c['id']} | {c['half']} | {c['outcome'][0]} -> {c['outcome'][1]} | {', '.join(c['moved']) or '-'} |"
            )
    else:
        lines.append("Identical: every outcome, SVG and PNG.")
    if any(t["missing"].values()):
        lines += ["", f"Entries missing from a run: {t['missing']}"]
    if measures:
        lines += ["", "## The parked measures (P29-13: the third sample, the NA controls, the science measure)", ""]
        for m in measures:
            lines += measure_lines(m)
    return "\n".join(lines).rstrip() + "\n"


# ------------------------------------------------------------------------------------------------ sheets


def svg_box(path: Path) -> tuple[float, float, float]:
    """(width, height) of an SVG page in points and its viewBox units per point."""
    head = path.read_text(encoding="utf-8", errors="replace")[:2000]
    tag = re.search(r"<svg\b[^>]*>", head).group(0)

    def pt(name: str) -> float:
        v = re.search(rf'\b{name}="([0-9.]+)(pt|px|in)?"', tag)
        x, unit = float(v.group(1)), v.group(2) or "px"
        return x * {"pt": 1.0, "px": 0.75, "in": 72.0}[unit]

    w, h = pt("width"), pt("height")
    vb = re.search(r'viewBox="([-0-9. ]+)"', tag)
    per = float(vb.group(1).split()[2]) / w if vb else 1.0
    return w, h, per


def svg_bbox(path: Path, pattern: str) -> tuple[float, float, float, float] | None:
    """The box (viewBox units) of every element whose id matches `pattern`: its paths' points and texts' anchors."""
    text = path.read_text(encoding="utf-8", errors="replace")
    xs, ys = [], []
    for m in re.finditer(r'<g id="([^"]+)">(.*?)</g>', text, re.S):
        if not re.search(pattern, m.group(1)):
            continue
        body = m.group(2)
        for d in re.findall(r'\bd="([^"]+)"', body):
            nums = [float(v) for v in re.findall(r"-?[0-9]+(?:\.[0-9]+)?", d)]
            xs += nums[0::2]
            ys += nums[1::2]
        for x, y in re.findall(r'\bx="(-?[0-9.]+)" y="(-?[0-9.]+)"', body):
            xs.append(float(x))
            ys.append(float(y))
    if not xs:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def _panel(folder: Path | None, crop: tuple | None):
    """(image, its page box in points, px per point) of a folder's app.png, or a grey note when it has none."""
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    png, svg = (folder / "app.png", folder / "app.svg") if folder else (None, None)
    if not png or not png.is_file() or not svg.is_file():
        return None
    img = Image.open(png).convert("RGB")
    w, h, _ = svg_box(svg)
    per_pt = img.width / w
    if crop is not None:
        x0, y0, x1, y1 = crop[0], crop[1], min(crop[2], w), min(crop[3], h)
        img = img.crop((int(x0 * per_pt), int(y0 * per_pt), int(math.ceil(x1 * per_pt)), int(math.ceil(y1 * per_pt))))
        w, h = x1 - x0, y1 - y0
    return img, (w, h), per_pt


def pair_sheet(before: Path | None, after: Path | None, out: Path, title: str, labels: tuple[str, str], crop_re=None):
    """before | after at one scale (px per point of the page) into `out`; the scale used."""
    from PIL import Image, ImageDraw, ImageFont

    crop = None
    if crop_re:
        boxes = [svg_bbox(f / "app.svg", crop_re) for f in (before, after) if f and (f / "app.svg").is_file()]
        boxes = [b for b in boxes if b]
        if boxes:
            per = svg_box((before or after) / "app.svg")[2]
            x0, y0 = min(b[0] for b in boxes) / per - 30, min(b[1] for b in boxes) / per - 30
            x1, y1 = max(b[2] for b in boxes) / per + 30, max(b[3] for b in boxes) / per + 30
            crop = (max(0.0, x0), max(0.0, y0), x1, y1)
    panels = [_panel(f, crop) for f in (before, after)]
    have = [p for p in panels if p]
    big = max((max(p[1]) for p in have), default=1.0)
    scale = min([p[2] for p in have] + [SHEET_PX / big])
    head, gap = 44, 24
    sized = []
    for p in panels:
        if p is None:
            sized.append(None)
            continue
        img, (w, h), _ = p
        sized.append(img.resize((max(1, round(w * scale)), max(1, round(h * scale))), Image.LANCZOS))
    widths = [s.width if s else 400 for s in sized]
    height = max([s.height for s in sized if s] + [200]) + head
    sheet = Image.new("RGB", (sum(widths) + gap * 3, height + gap), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.load_default(size=20)
    except TypeError:  # Pillow before 10.1
        font = ImageFont.load_default()
    x = gap
    for s, w, lab in zip(sized, widths, labels):
        draw.text((x, 8), f"{title}  {lab}", fill="black", font=font)
        if s is None:
            draw.text((x, head + 20), "no figure", fill="grey", font=font)
        else:
            sheet.paste(s, (x, head))
        x += w + gap
    out.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(out, optimize=True)
    return scale


def _sha(path: Path | None) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path and path.is_file() else None


def cmd_sheets(args) -> int:
    app29 = _app29()
    rerun, after = app29.load_both(args.rerun, args.rerun_out), app29.load_both(args.after, args.after_out)
    a, b = _rows(rerun), _rows(after)
    commit = lambda rec: (((rec.get("source_stamp") or {}).get("git") or {}).get("commit") or "?")[:7]  # noqa: E731
    labels = (f"before ({commit(rerun)})", f"after ({commit(after)})")
    pairs = []  # (name, id, half, why, before folder, after folder, crop)
    for i in sorted(a.keys() & b.keys()):
        if (a[i].get("files") or {}).get("app.svg") != (b[i].get("files") or {}).get("app.svg"):
            why = "its SVG changed" + (f"; lead: {LEADS[i]}" if i in LEADS else "")
            pairs.append((i, i, b[i].get("half"), why, Path(a[i]["folder"]), Path(b[i]["folder"]), None))
    listed = {p[0] for p in pairs}
    for i, why in LEADS.items():
        if i in listed:
            continue
        if i in a and i in b:
            pairs.append((i, i, b[i].get("half"), f"lead: {why}", Path(a[i]["folder"]), Path(b[i]["folder"]), None))
        else:
            fa, fb = Path(args.leads_before) / i, Path(args.leads_after) / i
            pairs.append((i, i, "lead", f"lead: {why}", fa, fb, None))
    for i, (pattern, why) in CROPS.items():
        if i in a and i in b:
            fa, fb = Path(a[i]["folder"]), Path(b[i]["folder"])
        else:
            fa, fb = Path(args.leads_before) / i, Path(args.leads_after) / i
        pairs.append((f"{i}-crop", i, "crop", why, fa, fb, pattern))
    out = Path(args.out)
    rows = []
    for name, i, half, why, fa, fb, crop in pairs:
        target = out / f"{name}.png"
        scale = pair_sheet(fa, fb, target, name, labels, crop)
        rows.append(
            {
                "pair": name,
                "id": i,
                "half": half,
                "why": why,
                "sheet": target.relative_to(ROOT).as_posix() if target.resolve().is_relative_to(ROOT) else str(target),
                "sheet_sha256": _sha(target),
                "px_per_pt": round(scale, 4),
                "before_svg": _sha(fa / "app.svg"),
                "after_svg": _sha(fb / "app.svg"),
                "outcome": [outcome(a[i]), outcome(b[i])] if i in a and i in b else None,
            }
        )
        print(f"{name}: {target}", flush=True)
    look = Path(args.look)
    doc = {"source_stamp": _stamp(), "labels": labels, "pairs": rows}
    look.with_suffix(".json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    look.write_text(render_look(doc), encoding="utf-8", newline="\n")
    return 0


def render_look(doc: dict) -> str:
    lines = [
        "# Plan 29: the look test",
        "",
        f"Written by scripts/app29_after.py sheets at {_stamp_line(doc['source_stamp'])}: every entry whose app SVG "
        "changed between the before commit rerun and the after run (both halves, gpu-desktop), the leads, and the "
        "named crops, each drawn as the app's PNG before beside after at one scale (`px/pt`: pixels per point of "
        f"the page). Before: {doc['labels'][0]}; after: {doc['labels'][1]}. The sheets are under "
        "docs/superpowers/work/app29/after/sheets/ (not committed; the sha256 names each). The verdict (better, same "
        "or worse) and its one-line reason are the controller's, written by hand after looking (plan 29 Task N+1 "
        'Step 3); any "worse" is reverted or explained before main.',
        "",
        "| # | pair | half | why listed | outcome before -> after | px/pt | sheet sha256 | before SVG | after SVG "
        "| verdict | reason |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for k, r in enumerate(doc["pairs"], 1):
        oc = " -> ".join(r["outcome"]) if r["outcome"] else "-"
        lines.append(
            f"| {k} | {r['pair']} | {r['half']} | {r['why']} | {oc} | {r['px_per_pt']} | {r['sheet_sha256']} | "
            f"{(r['before_svg'] or 'none')[:12]} | {(r['after_svg'] or 'none')[:12]} |  |  |"
        )
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("report", help="the before/after table per half and the parked measures")
    r.add_argument("--before", required=True, help="the committed before record (app29 run --record)")
    r.add_argument("--rerun", required=True, help="the before commit rerun on the after run's host")
    r.add_argument("--after", required=True, help="the after record")
    r.add_argument("--measure", nargs="*", help="app29_light.py measures (campaign, navalidate, science)")
    r.add_argument("--out", required=True, help="FILE.md (and FILE.json beside it)")
    s = sub.add_parser("sheets", help="the look test's before/after sheets and look.md")
    s.add_argument("--rerun", required=True)
    s.add_argument("--rerun-out", required=True, help="the rerun's run folder (its figures)")
    s.add_argument("--after", required=True)
    s.add_argument("--after-out", required=True, help="the after run's folder (its figures)")
    s.add_argument("--leads-before", required=True, help="app29 run --ids of the leads outside the sample, before")
    s.add_argument("--leads-after", required=True, help="the same, after")
    s.add_argument("--out", required=True, help="the sheets' folder")
    s.add_argument("--look", required=True, help="look.md (and look.json beside it)")
    args = ap.parse_args(argv)
    return {"report": cmd_report, "sheets": cmd_sheets}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
