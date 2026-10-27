# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The record of the readability floor (plan 28 Task 4, the review's 3.3, rulings Q39-Q41): the largest figures
measured at the default page size, which floor reasons fire on them, and what each existing option that draws less
achieves.

    python scripts/readability_floor.py [--out docs/validation/readability-floor] [--jobs 12] [--ids 9FZL,8XJV]

Writes report.json and report.md, stamped with the code that measured (foldmap.source_stamp). Needs network the first
time (the entries are fetched into ~/.cache/foldmap); not run in CI.

- census: from the three samples' committed campaign reports (docs/validation/random100*/report.json[.gz], at the
  commit each names): the default views whose notes say the page held its text under 6 pt (layout.page_fit), those
  over 60 s (A2) and 90 s (the hosted limit, A2b), and those whose route debug record counts a second route search
  refused for the page's budget; the largest entries of each sample by its residues (sample.json).
- entries: 9FZL and 8XJV (the review's), the LARGEST largest of each sample, and every default view the census names,
  each drawn as `foldmap plot ID` draws it (spec.draw_spec, the default options) in a process of its own, JOBS at a
  time (so times include that load), measured by floor.measures (page, the figure's text and its smallest text at
  the page size, named fallbacks, seconds), with A5's element overlaps and A6's text collisions (figcheck) and the
  count of texts drawn under 6 pt; its floor reasons (floor.reasons_of), its time note (floor.time_notes, never
  judged), the options the notes offer (floor.suggestions) and its other checks' verdicts (criteria.check_figure:
  the S-checks with findings and the A-criteria that fail), which say whether R1 alone makes `foldmap check` exit 1.
- tried: each option offered, drawn and measured the same way, and the preset measured (PRESET_KEYS: four existing
  style keys that leave out low-priority annotations; not an option: ruling Q40 dropped --reduced-detail after this
  measurement) on every entry."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    sys.path.insert(0, str(ROOT / "src"))  # the foldmap of this checkout (child processes inherit sys.path)

from foldmap import floor  # noqa: E402
from foldmap.criteria import LIMITS  # noqa: E402
from foldmap.source_stamp import source_stamp, stamp_line  # noqa: E402

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts
OUT = ROOT / "docs" / "validation" / "readability-floor"
SAMPLES = ("random100", "random100-heldout", "random100-third")
NAMED = ("9FZL", "8XJV")  # the review's 3.3
LARGEST = 3  # the largest entries of each sample measured (by sample.json residues)
SMALL_NOTE = re.compile(r"so its text is ([0-9.]+) pt, under 6 pt")
PRESET = "preset"  # the preset measured, not an option (ruling Q40): these style keys, given as --set
PRESET_KEYS = ("nucleotide_labels=false", "ligands=none", "glycans=false", "helices_310=false")


def _report(name: str) -> dict:
    plain = ROOT / "docs" / "validation" / name / "report.json"
    if plain.is_file():
        return json.loads(plain.read_text(encoding="utf-8"))
    return json.loads(gzip.decompress((plain.parent / "report.json.gz").read_bytes()))


def census() -> dict:
    """Per sample, from its committed report and sample.json: the default views below or near the floor's limits."""
    out = {}
    for name in SAMPLES:
        doc = _report(name)
        sample = json.loads((ROOT / "docs" / "validation" / name / "sample.json").read_text(encoding="utf-8"))
        size = {e["id"]: e.get("residues") or 0 for e in sample["entries"]}
        small, slow, over, refused = {}, {}, {}, {}
        for e in doc["entries"]:
            view = (e.get("views") or {}).get("default") or {}
            said = SMALL_NOTE.search(" ".join(view.get("stderr_tail") or []))
            if said:
                small[e["id"]] = float(said.group(1))
            seconds = view.get("seconds") or 0
            if seconds > floor.TIME_LIMIT_S:
                over[e["id"]] = seconds
            if seconds > LIMITS["default_s"]:  # A2
                slow[e["id"]] = seconds
            route = (view.get("result") or {}).get("route") or {}
            if route.get("second_searches_refused"):
                refused[e["id"]] = route["second_searches_refused"]
        largest = sorted(size, key=lambda i: (-size[i], i))[:LARGEST]
        out[name] = {
            "commit": doc.get("foldmap_commit"),
            "default_views": len(doc["entries"]),
            "text_under_6pt": small,
            "over_60_s": slow,
            "over_90_s": over,
            "route_budget_refused": refused,
            "largest": {i: size[i] for i in largest},
        }
    return out


def picks(found: dict) -> list[str]:
    ids = list(NAMED)
    for name in SAMPLES:
        got = found[name]
        ids += [*got["largest"], *got["text_under_6pt"], *got["over_60_s"], *got["over_90_s"]]
    return list(dict.fromkeys(ids))


def _spec_kw(options: dict) -> dict:
    """figure_spec's keywords for an option set: {focus: X}, {chains: X} or {PRESET: None}."""
    kw: dict = {}
    for k, v in options.items():
        if k == PRESET:
            kw["style"] = list(PRESET_KEYS)
        elif k == "chains":
            kw["chains"] = [v]
        else:
            kw[k] = v
    return kw


def phrase(options: dict) -> str:
    if PRESET in options:
        return "the preset measured (not an option)"
    return " ".join(f"--{k} {v}" for k, v in options.items()) or "(default)"


def figure(task: tuple[str, dict, bool]) -> dict:
    """One figure drawn and measured (module docstring); a failure is recorded, never raised."""
    entry, options, offered = task
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from foldmap import figcheck
    from foldmap.cli import resolve_structure
    from foldmap.spec import draw_spec, figure_spec

    row: dict = {"options": options, "said": phrase(options), "offered": offered}
    try:
        path = resolve_structure(entry)
        spec = figure_spec(path, **_spec_kw(options))
        drawn = draw_spec(spec)
    except Exception as err:  # noqa: BLE001 (recorded: the record says which figure could not be drawn)
        return {**row, "error": f"{type(err).__name__}: {err}"[:400]}
    m = floor.measures(drawn)
    reasons = floor.reasons_of(m)
    from matplotlib.text import Text

    sizes = [t.get_fontsize() for t in drawn.figure.findobj(Text) if t.get_visible() and t.get_text().strip()]
    t0 = time.perf_counter()
    checks = figcheck.figure_checks(drawn, path, assembly=spec.assembly, focus=spec.focus, chains=spec.chains)
    lay = drawn.layout
    strands = lay.nucleic.strands if lay.nucleic is not None else []
    row.update(
        measures=m,
        reasons=[r["reason"] for r in reasons],
        texts=len(sizes),
        texts_under_6pt=sum(s < floor.MIN_TEXT_PT for s in sizes),
        element_overlaps=checks.get("element_overlaps"),
        text_collisions=checks.get("text_collisions"),
        stuck=len(drawn.stuck),
        elements=len(lay.placed),
        chains_drawn=len({p.sse.chain for p in lay.placed.values()} | {s.chain for s in strands}),
        check_seconds=round(time.perf_counter() - t0, 1),
    )
    if not options:
        row["suggestions"] = floor.suggestions(drawn) if reasons else []  # the notes offer them below it only
        row["notes"] = [t["reason"] for t in floor.time_notes(m)]  # said, never judged (Q40)
        row["note"] = floor.note(reasons, row["suggestions"], timed=floor.time_notes(m))
        from foldmap.criteria import check_figure

        found = check_figure(drawn, spec, path)  # what `foldmap check` judges besides R1 (its exit code)
        row["science_findings"] = sorted(k for k, v in found["science"].items() if v)
        row["criteria_failed"] = sorted(k for k, v in found["criteria"].items() if v["pass"] is False)
    plt.close(drawn.figure)
    return row


def measure(ids: list[str], jobs: int) -> list[dict]:
    with ProcessPoolExecutor(jobs) as ex:
        defaults = list(ex.map(figure, [(i, {}, False) for i in ids]))
        tasks = []
        for i, d in zip(ids, defaults):
            if "error" in d:
                continue
            tasks += [(i, s["options"], True) for s in d["suggestions"]]
            tasks.append((i, {PRESET: None}, False))
        tried = list(ex.map(figure, tasks))
    rows = []
    for i, d in zip(ids, defaults):
        rows.append({"id": i, "default": d, "tried": [t for (e, _, _), t in zip(tasks, tried) if e == i]})
    return rows


def effects(rows: list[dict]) -> dict:
    """Per option kind (focus protomer, focus CHAIN, chains CHAIN, the preset): on how many figures below the floor
    it was tried, on how many it brought the figure above, and the median seconds and text size it gave."""
    out: dict = {}
    for r in rows:
        if "error" in r["default"] or not r["default"]["reasons"]:
            continue
        for t in r["tried"]:
            if "error" in t:
                continue
            k = next(iter(t["options"]))
            kind = (
                phrase(t["options"])
                if k == PRESET
                else f"--{k} " + ("protomer" if t["options"][k] == "protomer" else "CHAIN")
            )
            got = out.setdefault(kind, {"tried_below": 0, "brought_above": 0, "entries": [], "seconds": [], "pt": []})
            got["tried_below"] += 1
            got["brought_above"] += not t["reasons"]
            got["entries"].append(r["id"])
            got["seconds"].append(t["measures"]["seconds"])
            got["pt"].append(t["measures"]["text_pt"])
    for got in out.values():
        got["median_seconds"] = statistics.median(got.pop("seconds"))
        got["median_text_pt"] = statistics.median(got.pop("pt"))
    return out


def _fallbacks(m: dict) -> str:
    fired = m["fallbacks"]
    return ", ".join(f"{k} {n}" for k, n in fired.items()) or "none"


def _line(fig: dict) -> str:
    if "error" in fig:
        return f"| {fig['said']} | error: {fig['error']} |||||||||"
    m = fig["measures"]
    w, h = m["page_in"]
    return (
        f"| {fig['said']} | {w:.1f} x {h:.1f} | {m['text_pt']:.2f} | {m['smallest_pt']:.2f} | "
        f"{fig['texts_under_6pt']} / {fig['texts']} | {m['seconds']:.1f} | {_fallbacks(m)} | "
        f"{fig['element_overlaps']} / {fig['text_collisions']} | {', '.join(fig['reasons']) or 'above'} |"
    )


def _preset_worse(rows: list[dict]) -> str:
    """The figures the preset took from above the floor to below it, said (', and on X it took ...'), or ''."""
    worse = []
    for r in rows:
        d = r["default"]
        if "error" in d or d["reasons"]:
            continue
        for t in r["tried"]:
            if PRESET in t["options"] and "error" not in t and t["reasons"]:
                was, now = d["measures"]["text_pt"], t["measures"]["text_pt"]
                worse.append(f"{r['id']} ({was:.2f} pt to {now:.2f} pt)")
    if not worse:
        return ""
    return "; on " + ", ".join(worse) + " it took a figure above the floor to below it"


def _flips(rows: list[dict]) -> str:
    """Which default figures below the floor passed every other check (R1 alone now makes check exit 1) and which
    already failed one, said as a sentence."""
    alone, already = [], []
    for r in rows:
        d = r["default"]
        if "error" in d or not d["reasons"]:
            continue
        other = d.get("science_findings", []) + d.get("criteria_failed", [])
        (already if other else alone).append(f"{r['id']} ({', '.join(other)})" if other else r["id"])
    return (
        f"passed on {len(alone)} of them, which now exit 1 on R1 alone: {', '.join(alone) or 'none'}; and already "
        f"failed on {len(already)}, which exited 1 before: {', '.join(already) or 'none'}."
    )


def render(doc: dict) -> str:
    rows, found = doc["entries"], doc["census"]
    below = [r for r in rows if r["default"].get("reasons")]
    lines = [
        "# The readability floor",
        "",
        f"Plan 28 Task 4 (the review's 3.3, rulings Q39-Q41). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "A figure below the readability floor says so (floor.readability_floor): in `foldmap plot`'s notes (and with "
        "`--strict-floor` it writes no figure and exits 1), in `foldmap check` as R1 (readability; a fail exits 1) "
        "and its --json report, in the MCP result and under the app's figure, naming the existing options that draw "
        "less. The figure itself is unchanged. This record is what the reasons and their limits rest on.",
        "",
        "## The reasons and their limits",
        "",
        f"- **text_small**: the figure's own text (element names, chain ends) under {floor.MIN_TEXT_PT:g} pt on its "
        "page (layout.MIN_TEXT_PT, the random-100 C4 floor): the page grows to keep it at 6 pt until "
        "layout.MAX_PAGE_IN (48 in), and past that the whole figure is scaled down. Measured here at the default "
        "page size, as `foldmap plot` writes it. The smallest text drawn (residue counts, base letters, ligand "
        "names, 0.6-0.9 of the figure's text by design) is recorded but not judged: it is under 6 pt on "
        f"{sum(1 for r in rows if 'error' not in r['default'] and r['default']['measures']['smallest_pt'] < 6)} of "
        f"the {sum(1 for r in rows if 'error' not in r['default'])} default figures measured, those far above the "
        "floor included.",
        "- **route_budget**: the page's budget for second, wider route searches (route._MAX_NODES_PAGE) spent, "
        "counted as the named fallbacks route_second_search_refused. It fired on none of the figures measured here "
        "and on none of the samples' default views (census below: their reports count none refused), so it is "
        "tested on synthetic measures only here. It does fire past the default views: 7NZ2's domains view (the "
        "held-out sample) spends the budget, 6 searches refused (docs/validation/prune/report.md, held), so "
        "`foldmap check 7NZ2 --domains auto` fails R1 on it and R1's exit-code change reaches past the default views. "
        "A stuck loop of another kind is A3's, as it happens on small pages too.",
        f"- **time**: reported, never judged (ruling Q40). A figure that took longer than the hosted app's "
        f"{floor.TIME_LIMIT_S:g} s (criteria.LIMITS app_s) to make, measured as spec.draw_spec's wall time, is said "
        "so in the notes, the check's text, the MCP note and the app; it never fails R1, never changes `foldmap "
        "check`'s exit code or --strict-floor, and is not in the check's --json report (which holds no wall time). "
        "It depends on the machine and its load: the figures here were drawn "
        f"{doc['jobs']} at a time. Time notes here: "
        + (
            ", ".join(
                f"{r['id']} ({r['default']['measures']['seconds']:.0f} s)" for r in rows if r["default"].get("notes")
            )
            or "none"
        )
        + ".",
        "",
        "## Census of the three samples' default views",
        "",
        "From each sample's committed campaign report, at the commit it names (not this code): the default views "
        "whose notes say the page held its text under 6 pt, those over 60 s (A2) and over 90 s (A2b), and those "
        "with a second route search refused for the page's budget.",
        "",
        "| sample | commit | default views | text under 6 pt | over 60 s | over 90 s | route budget spent |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, got in found.items():

        def ids(d: dict) -> str:
            return ", ".join(f"{i} ({v})" for i, v in d.items()) or "none"

        lines.append(
            f"| {name} | {str(got['commit'])[:12]} | {got['default_views']} | {ids(got['text_under_6pt'])} | "
            f"{ids(got['over_60_s'])} | {ids(got['over_90_s'])} | {ids(got['route_budget_refused'])} |"
        )
    lines += [
        "",
        "## What each option achieves on the figures below the floor",
        "",
        f"Below the floor at this code: {len(below)} of the {len(rows)} figures measured, "
        + ", ".join(f"{r['id']} ({', '.join(r['default']['reasons'])})" for r in below)
        + ".",
        "",
        "| option | tried on | brought above | median seconds | median text pt | entries |",
        "|---|---|---|---|---|---|",
    ]
    for kind, got in doc["effects"].items():
        lines.append(
            f"| `{kind}` | {got['tried_below']} | {got['brought_above']} | {got['median_seconds']:.1f} | "
            f"{got['median_text_pt']:.2f} | {', '.join(got['entries'])} |"
        )
    lines += [
        "",
        f"The preset measured is `{' '.join('--set ' + s for s in PRESET_KEYS)}` (existing style keys "
        "only), not an option: ruling Q40 dropped the --reduced-detail option this measured, as it brought none of "
        "the figures below the floor above it. It leaves out low-priority annotations, and so most of the text drawn "
        "under 6 pt, but the page's size is set by its elements, so it does not raise the figure's text"
        + _preset_worse(rows)
        + ".",
        "",
        "## What `foldmap check` now says",
        "",
        "R1 fails, and `foldmap check` exits 1, on the figures below the floor. Before R1 these default figures' "
        "other checks (the S-checks and A1-A20 judged here, criteria.check_figure) "
        + _flips(rows)
        + " The census above is of the samples' reports at their own commit, so that every sample default whose "
        "text is under 6 pt is among those measured here is inferred from it, not run at this code on all 300.",
        "",
        "## Every figure",
        "",
        "Page in inches; the figure's text and its smallest text in points at that page size; texts drawn under 6 pt "
        "of all texts; seconds to make (draw_spec); named fallbacks taken; A5 element overlaps / A6 text "
        "collisions; the floor's reasons. Offered: the options the notes name for that figure.",
    ]
    for r in rows:
        d = r["default"]
        size = next((got["largest"][r["id"]] for got in found.values() if r["id"] in got["largest"]), None)
        lines += ["", f"### {r['id']}" + (f" ({size:,} residues deposited)" if size else ""), ""]
        if d.get("note"):
            lines += [f"Notes: {d['note']}", ""]
        lines += [
            "| options | page in | text pt | smallest pt | under 6 pt | s | fallbacks | overlaps / collisions "
            "| floor |",
            "|---|---|---|---|---|---|---|---|---|",
            _line(d),
        ]
        lines += [_line({**t, "said": t["said"] + (" (offered)" if t["offered"] else "")}) for t in r["tried"]]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(OUT))
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--ids", help="only these entries (comma-separated): a trial, never the record")
    args = parser.parse_args(argv)
    found = census()
    ids = args.ids.split(",") if args.ids else picks(found)
    rows = measure(ids, args.jobs)
    doc = {
        "source_stamp": SOURCE_STAMP,
        "jobs": args.jobs,
        "limits": {"text_pt": floor.MIN_TEXT_PT, "route_budget_refused": 0},
        "time_note_s": floor.TIME_LIMIT_S,  # said, never judged (Q40)
        "preset_measured": list(PRESET_KEYS),  # not an option (Q40)
        "census": found,
        "entries": rows,
        "effects": effects(rows),
    }
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8", newline="\n")
    (out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(f"wrote {out / 'report.json'} and report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
