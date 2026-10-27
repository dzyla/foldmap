# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The record of the topology carried by `--like` (plan 28 Task 6, the review's 3.5, ruling Q44): whether a sheet
drawn like a reference keeps the reference's strand order and direction, measured on the homologue tuning and
held-out sets, with the reference's order carried (Like.order, as foldmap draws) and withheld (the mechanism behind
3WZ0-2AV5's mirrored sheet in the held-out run at c099787).

    python scripts/topology_transfer.py [--out docs/validation/homologues/topology-transfer] [--jobs 16]
                                        [--pairs 3WZ0_A-2AV5_A,...]

Writes report.json and report.md, stamped with the code that measured (foldmap.source_stamp). Needs network the first
time (the entries are fetched into ~/.cache/foldmap); not run in CI.

The sets, with their tags read from scripts/paper/sets.json: hom_tuning (docs/validation/homologues/tuning.json, the
pairs the H1-H5 limits were frozen on) and hom_heldout (heldout.json, run once at c099787; 3WZ0-2AV5 is one of its
pairs, and the fix measured here, 1a3f9d6, was aimed at it after that run).

Each pair is drawn as scripts/homologues.py's `pair` draws it (B like A as `foldmap check B --like A --like-chain CA
--chains CB --assembly asu`, and B alone), in a child process of its own (python -I, PYTHONHASHSEED=0), twice like A:
- carried: as foldmap draws it, each of B's sheets reading left to right as A reads its mates where like._order says
  (Like.order);
- withheld: the same with like._order answering None for every sheet, so each sheet's order is its own
  (layout._sheet_order's default: decided by its end strands' projected centroids).
Each drawing is judged by H1-H5 (criteria.like_checks, as the campaign judges it). For every sheet of B (as
like.like_from receives it, read by wrapping like_from; nothing else is changed): its strands and mated strands,
Like.order's entry, why it is None (_why_none: the conditions of like._order, in its order), and, where it is set,
whether B's drawn mated strands read left to right and point up or down as their mates do in A.

Withheld is today's code with like._order answering None, judged by today's criteria, not the code before 1a3f9d6:
that commit also removed H3's pairing excuse and made H1's refusal a flag; withheld reproduces the c099787 counts
(16 of 30, H2 17, H3 17) all the same."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    sys.path.insert(0, str(ROOT / "src"))  # the foldmap of this checkout

from foldmap.source_stamp import source_stamp, stamp_line  # noqa: E402

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts
HOM = ROOT / "docs" / "validation" / "homologues"
OUT = HOM / "topology-transfer"
SETS = {"hom_tuning": HOM / "tuning.json", "hom_heldout": HOM / "heldout.json"}  # sets.json key -> set file
TAGS = ROOT / "scripts" / "paper" / "sets.json"
CHECKS = ("H1", "H2", "H3", "H4", "H5")
REGRESSION = "3WZ0_A-2AV5_A"  # the review's 3.5: B C A D drawn D A C B (heldout-report.md, c099787)
WAYS = ("carried", "withheld")
PAIR_TIMEOUT_S = 1800  # a pair's child is killed after this many seconds (scripts/homologues.py's)


def stem(pair: dict) -> str:
    """A pair's file stem: A_CHAIN-B_CHAIN (scripts/homologues.py's)."""
    return f"{pair['a']['id']}_{pair['a']['chain']}-{pair['b']['id']}_{pair['b']['chain']}"


def read_set(path) -> list[dict]:
    return json.loads(Path(path).read_text(encoding="utf-8"))["pairs"]


def _why_none(lay_a, sheet, mates) -> str:
    """Why like._order says nothing of this B sheet, its conditions in its order: fewer than two track strands mated
    in one of A's sheet blocks, mates in two blocks (B pairs strands A does not: a pairing difference), mates
    pointing otherwise than B's pairing has them, or mates not in one order across A's sheet."""
    import numpy as np

    placed = lay_a.placed
    block = {k: n for n, ids in enumerate(lay_a.sheets) for k in ids}
    mated = [
        (k, mates.mate[s.id])
        for k, s in enumerate(sheet.strands)
        if s.id in mates.mate and mates.mate[s.id] in placed and mates.mate[s.id] in block
    ]
    if len(mated) < 2:
        return "fewer than two mated"
    if len({block[a] for _, a in mated}) != 1:
        return "mates in two of A's sheets"
    ups = [(1 if np.sin(placed[a].angle) > 0 else -1) * sheet.directions[k] for k, a in mated]
    if len(set(ups)) != 1:
        return "mates point otherwise than B's pairing"
    return "mates not in one order across A's sheet"


def sheet_rows(lay_a, lay_b, sheets, like, mates, withheld: bool = False) -> list[dict]:
    """Each B sheet: strands, mated track strands, closed, Like.order's entry and why it is None ("withheld" where
    like._order was made to answer None: its own conditions were not asked); and, with two mated
    or more, whether B's mated track strands read left to right (`reads`) and point up or down (`points`) as their
    mates do in A, and whether its mated accessory strands (off the track, layout._drawn_strands) read so too, with
    the track's (`reads_all`)."""
    import numpy as np

    def reads(ids: list[str]) -> bool:
        xb = [lay_b.placed[x].cx for x in ids]
        xa = [lay_a.placed[mates.mate[x]].cx for x in ids]
        return [int(i) for i in np.argsort(xb, kind="stable")] == [int(i) for i in np.argsort(xa, kind="stable")]

    def mated(strands) -> list[str]:
        return [s.id for s in strands if mates.mate.get(s.id) in lay_a.placed and s.id in lay_b.placed]

    rows = []
    for k, sheet in enumerate(sheets):
        track = mated(sheet.strands)
        every = mated([*sheet.strands, *(x.strand for x in sheet.accessory)])
        row = {"strands": len(sheet.strands), "mated": len(track), "closed": bool(sheet.closed)}
        row["order"] = like.order[k] if like.order is not None else None
        if row["order"] is None:
            row["why"] = "withheld" if withheld else _why_none(lay_a, sheet, mates)
        if len(track) >= 2:
            row["reads"] = reads(track)
            row["points"] = all(
                (np.sin(lay_b.placed[x].angle) > 0) == (np.sin(lay_a.placed[mates.mate[x]].angle) > 0) for x in track
            )
            row["reads_all"] = reads(every)
        rows.append(row)
    return rows


def _verdicts(h: dict) -> dict:
    """H1-H5's verdicts (True, False, None: not judged), H2's median and p90, H3's share and its sheet orders."""
    out = {k: (h.get(k) or {}).get("pass") for k in CHECKS}
    h2 = (h.get("H2") or {}).get("value")
    h3 = (h.get("H3") or {}).get("value")
    out["h2"] = [h2["median"], h2["p90"]] if isinstance(h2, dict) and "median" in h2 else None
    out["h3"] = h3["share"] if isinstance(h3, dict) and "share" in h3 else None
    out["sheet_orders"] = list(h3.get("sheets") or []) if isinstance(h3, dict) else []
    return out


def _capturing(like_from, seen: dict):
    """like.like_from, unchanged, recording B's sheets as it receives them (`sheets_b`) into seen['sheets']."""

    def wrapped(*a, **k):
        seen["sheets"] = k["sheets_b"] if "sheets_b" in k else (a[5] if len(a) > 5 else None)
        return like_from(*a, **k)

    return wrapped


def measure_pair(pair: dict) -> dict:
    """One pair drawn like A carried and withheld and alone, judged (module doc)."""
    import foldmap.like as like_mod
    from foldmap import cli
    from foldmap.criteria import check_figure, jsonable, like_checks
    from foldmap.spec import draw_spec

    a, b = pair["a"], pair["b"]
    args = cli.build_parser().parse_args(["check", b["id"], "--like", a["id"], "--like-chain", a["chain"],
                                          "--chains", b["chain"], "--assembly", "asu"])  # fmt: skip
    args.structure = cli.resolve_structure(args.structure)
    spec = cli._spec_of(args)
    alone_spec = replace(spec, like=None)
    alone = draw_spec(alone_spec)
    found_alone = check_figure(alone, alone_spec, spec.structure)
    out = {"pair": stem(pair), "kind": pair["kind"], "stratum": pair.get("stratum")}
    like_from, order = like_mod.like_from, like_mod._order
    for way in WAYS:
        seen: dict = {}
        like_mod.like_from = _capturing(like_from, seen)
        if way == "withheld":
            like_mod._order = lambda *a, **k: None
        try:
            got: dict = {}
            drawn = draw_spec(spec, like_out=got)
        finally:
            like_mod.like_from, like_mod._order = like_from, order
        ref, made, mates, aln = got.get("reference"), got.get("like"), got.get("mates"), got.get("alignment")
        h = like_checks(
            ref[0] if ref else None,
            drawn.layout,
            aln,
            mates,
            drawn.layout.like_note,
            made is not None,
            alone.layout,
            check_figure(drawn, spec, spec.structure),
            found_alone,
            kind=pair["kind"],
            refused=bool(got.get("refused")),
        )
        result = {"alike": made is not None, **_verdicts(jsonable(h))}
        if made is not None and seen.get("sheets") is not None:
            result["sheets"] = sheet_rows(ref[0], drawn.layout, seen["sheets"], made, mates, way == "withheld")
        if aln is not None:
            out["alignment"] = {"method": aln.method, "tm_score": round(float(aln.tm_score), 4)}
        out[way] = result
    return out


def tags() -> dict[str, dict]:
    """sets.json's entry ({tag, what, reason}) of each set measured."""
    sets = json.loads(TAGS.read_text(encoding="utf-8"))["sets"]
    return {name: sets[name] for name in SETS}


def summary(rows: list[dict]) -> dict:
    """Per way: positives passing H1-H5 and each check, negatives refused (H1 passes), the sheets Like.order sets and
    how many of those read and point as A's, the None reasons, and the pairs with a sheet in another order; then the
    pairs whose verdicts differ between the two ways and the checks each loses carried."""
    out: dict = {}
    pos = [r for r in rows if r["kind"] == "positive"]
    neg = [r for r in rows if r["kind"] == "negative"]
    for way in WAYS:
        sheets = [s for r in rows for s in r[way].get("sheets", [])]
        hinted = [s for s in sheets if s["order"] is not None]
        why: dict[str, int] = {}
        for s in sheets:
            if s["order"] is None:
                why[s["why"]] = why.get(s["why"], 0) + 1
        out[way] = {
            "positives": len(pos),
            "all_pass": sum(all(r[way][k] is True for k in CHECKS) for r in pos),
            "pass": {k: sum(r[way][k] is True for r in pos) for k in CHECKS},
            "negatives": len(neg),
            "refused": sum(r[way]["H1"] is True for r in neg),
            "sheets": len(sheets),
            "order_set": len(hinted),
            "reads_and_points": sum(bool(s.get("reads") and s.get("points")) for s in hinted),
            "none_why": dict(sorted(why.items())),
            "pairs_with_sheet_orders": sorted(r["pair"] for r in rows if r[way]["sheet_orders"]),
        }
    changed = []
    for r in rows:
        diff = [k for k in CHECKS if r["carried"][k] != r["withheld"][k]]
        if diff or r["carried"]["sheet_orders"] != r["withheld"]["sheet_orders"]:
            lost = [k for k in diff if r["withheld"][k] is True and r["carried"][k] is not True]
            changed.append({"pair": r["pair"], "checks": diff, "lost": lost})
    out["changed"] = changed
    out["lost_any"] = sorted(c["pair"] for c in changed if c["lost"])
    return out


def _run_child(pair: dict, tmp: Path) -> dict:
    s = stem(pair)
    path = tmp / f"{s}.json"
    env = {**os.environ, "PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1"}
    cmd = [sys.executable, "-I", str(SCRIPT), "pair", json.dumps(pair), "--to", str(path)]
    try:
        done = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=PAIR_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return {"pair": s, "kind": pair["kind"], "error": f"timed out after {PAIR_TIMEOUT_S} s"}
    if done.returncode or not path.is_file():
        return {"pair": s, "kind": pair["kind"], "error": (done.stderr or done.stdout).strip()[-400:]}
    return json.loads(path.read_text(encoding="utf-8"))


def measure(jobs: int, only: set[str] | None) -> dict[str, list[dict]]:
    out = {}
    with tempfile.TemporaryDirectory() as tmp, ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for name, path in SETS.items():
            pairs = [p for p in read_set(path) if only is None or stem(p) in only]
            out[name] = list(pool.map(lambda p: _run_child(p, Path(tmp)), pairs))
            for r in out[name]:
                print(f"{name} {r['pair']}: {'error' if 'error' in r else 'ok'}", flush=True)
    return out


def _cell(v) -> str:
    return "–" if v is None else ("pass" if v else "**fail**")


def _orders(v: dict) -> str:
    return ", ".join(v["sheet_orders"]) or "none"


def _rp(s: dict) -> str:
    """A sheet row's reads / points, or – with fewer than two mated."""
    if "reads" not in s:
        return "–"
    return f"{'yes' if s['reads'] else 'no'} / {'yes' if s['points'] else 'no'}"


def _h(v: dict) -> str:
    h2 = "–" if v["h2"] is None else f"{v['h2'][0]:.2f} / {v['h2'][1]:.2f}"
    h3 = "–" if v["h3"] is None else f"{v['h3']:.3f}"
    return f"{' '.join(_cell(v[k]) for k in CHECKS)} | {h2} | {h3}"


def render(doc: dict) -> str:
    sets, tagged = doc["sets"], doc["tags"]
    reg = next((r for rows in sets.values() for r in rows if r["pair"] == REGRESSION and "error" not in r), None)
    lines = [
        "# The topology carried by --like",
        "",
        f"Plan 28 Task 6 (the review's 3.5, ruling Q44). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "Whether a sheet drawn like a reference keeps the reference's strand order and direction. Each pair is drawn "
        "like A twice: **carried**, as foldmap draws it (each of B's sheets reads left to right as A reads its mates, "
        "Like.order, where like._order says), and **withheld**, the same with like._order answering None, so each "
        "sheet's order is its own (layout._sheet_order's default). Both are judged by H1-H5 as the homologue "
        "campaign judges them (scripts/homologues.py; limits criteria.LIMITS, frozen at 76fd936). The record is "
        "made by scripts/topology_transfer.py; report.json holds every pair and sheet.",
        "",
        "## The sets",
        "",
    ]
    for name, t in tagged.items():
        lines.append(f"- **{name}** ({t['tag']}): {t['what']}; {t['reason']}.")
    lines += [
        "- Neither set is unseen for this record: hom_tuning is the set the limits were chosen on, and the fix "
        "measured here (Like.order, 1a3f9d6, plan 13a-12) was aimed at 3WZ0-2AV5, a hom_heldout pair, after its "
        "one run. The carried column is a regression measure on both sets, not evidence on unseen structures.",
        "",
        "- Withheld is today's code with like._order answering None, judged by today's criteria, not the code "
        "before 1a3f9d6: that commit also removed H3's pairing excuse and made H1's refusal a flag; withheld "
        "reproduces the c099787 counts (16 of 30, H2 17, H3 17) all the same.",
        "",
        "## The mechanism (3WZ0-2AV5)",
        "",
    ]
    if reg is not None:
        c, w = reg["carried"], reg["withheld"]
        ws = next((s for s in w.get("sheets", [])), {})
        cs = next((s for s in c.get("sheets", [])), {})
        lines += [
            f"TM-score {reg['alignment']['tm_score']:.2f} by {reg['alignment']['method']} alignment. A's four-strand "
            "sheet is seen nearly edge-on through the carried frame, and B's sheet has a fifth strand at one end "
            "that A lacks (C′).",
            "",
            "| way | Like.order | mates read / point as A's | H1 H2 H3 H4 H5 | H2 median / p90 | H3 share | "
            "sheet in another order |",
            "|---|---|---|---|---|---|---|",
            f"| withheld | None ({ws.get('why', '–')}) | {_rp(ws)} | {_h(w)} | {_orders(w)} |",
            f"| carried | {cs.get('order')} | {_rp(cs)} | {_h(c)} | {_orders(c)} |",
            "",
            "- **Where the order is lost** (withheld): layout._sheet_order draws a sheet's track left to right across "
            "its drawn-up direction, and reverses it when its last track strand's projected centroid lies left of its "
            "first. For a sheet seen nearly edge-on that step is a fraction of an ångström, and in B the last track "
            "strand is C′, the strand A lacks, so B's own end strands decide an order A's mates never had: B C A D "
            "comes out D A C B, each strand pointing as its mate, the mirror image of A's sheet. Nothing carried "
            "from A entered the order: like.like_from carried places, angles (the sheet's way up, from its first "
            "mated strand) and names, not the reading.",
            "- **The fix** (1a3f9d6): like._order reads A's drawn order of the mates along B's track (Like.order: "
            "True left to right, False right to left), and layout._sheet_order takes it in place of the end strands' "
            "step, the way up still from A's angle. It is set only where B's strands are mapped and paired as in A: "
            "two or more track strands mated in one of A's sheets, pointing as B's pairing has them, in one order "
            "there; else None and the sheet's own order (the reasons below).",
        ]
    else:
        lines.append(f"{REGRESSION} was not measured in this run.")
    lines += ["", "## Both sets, carried and withheld", ""]
    lines += [
        "| set | way | positives passing H1-H5 | H1 | H2 | H3 | H4 | H5 | negatives refused | sheets: order set | "
        "... read and point as A's | pairs with a sheet in another order |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name, s in doc["summary"].items():
        for way in WAYS:
            v = s[way]
            lines.append(
                f"| {name} | {way} | {v['all_pass']} of {v['positives']} | "
                + " | ".join(str(v["pass"][k]) for k in CHECKS)
                + f" | {v['refused']} of {v['negatives']} | {v['order_set']} of {v['sheets']} | "
                f"{v['reads_and_points']} of {v['order_set']} | {len(v['pairs_with_sheet_orders'])} |"
            )
    lines += [""]
    for name, s in doc["summary"].items():
        changed = s["changed"]
        lines.append(
            f"- **{name}**: {len(changed)} pair(s) whose verdicts or sheet orders differ carried and withheld"
            + (
                ": " + "; ".join(f"{c['pair']} ({', '.join(c['checks']) or 'sheet orders only'})" for c in changed)
                if changed
                else ""
            )
            + f". Pairs losing an H check carried: {', '.join(s['lost_any']) or 'none'}."
        )
    lines += ["", "## Sheets still in another order (carried)", ""]
    lines += [
        "Each pair with a sheet of B that H3 finds in another order than its mates in A (criteria.sheet_orders), "
        "and B's sheets with two mated strands or more there: Like.order's entry or why it is None, and whether the "
        "mated track strands read (r) and point (p) as their mates, and the accessory strands read with them (a).",
        "",
        "| set | pair | TM | H2 / H3 | sheets in another order (B's for A's) | B's sheets: Like.order, r p a |",
        "|---|---|---|---|---|---|",
    ]
    for name, rows in sets.items():
        for r in rows:
            if "error" in r or not r["carried"]["sheet_orders"]:
                continue
            c = r["carried"]
            each = [
                f"{'None: ' + s['why'] if s['order'] is None else s['order']}, "
                + " ".join("y" if s.get(k) else "n" for k in ("reads", "points", "reads_all"))
                for s in c.get("sheets", [])
                if s["mated"] >= 2
            ]
            lines.append(
                f"| {name} | {r['pair']} | {r['alignment']['tm_score']:.2f} | {_cell(c['H2'])} / {_cell(c['H3'])} | "
                f"{'; '.join(c['sheet_orders'])} | {'; '.join(each)} |"
            )
    lines += ["", "## Every pair", ""]
    lines += [
        "| set | pair | kind | TM | carried: H1 H2 H3 H4 H5 | H2 | H3 | withheld: H1 H2 H3 H4 H5 | H2 | H3 |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name, rows in sets.items():
        for r in rows:
            if "error" in r:
                lines.append(f"| {name} | {r['pair']} | {r['kind']} | error: {r['error'][:80]} | | | | | | |")
                continue
            tm = r.get("alignment", {}).get("tm_score")
            lines.append(
                f"| {name} | {r['pair']} | {r.get('stratum') or r['kind']} | "
                f"{'–' if tm is None else f'{tm:.2f}'} | {_h(r['carried'])} | {_h(r['withheld'])} |"
            )
    errors = [r["pair"] for rows in sets.values() for r in rows if "error" in r]
    if errors:
        lines += ["", f"Errors: {', '.join(errors)}."]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if argv[:1] == ["pair"]:  # a child: one pair, its result written to --to
        child = argparse.ArgumentParser()
        child.add_argument("pair")
        child.add_argument("--to", required=True)
        a = child.parse_args(argv[1:])
        Path(a.to).write_text(json.dumps(measure_pair(json.loads(a.pair)), sort_keys=True), encoding="utf-8")
        return 0
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(OUT))
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--pairs", help="only these pairs (A_CHAIN-B_CHAIN, comma-separated): a trial, not the record")
    args = parser.parse_args(argv)
    sets = measure(args.jobs, set(args.pairs.split(",")) if args.pairs else None)
    doc = {
        "source_stamp": SOURCE_STAMP,
        "tags": tags(),
        "summary": {name: summary([r for r in rows if "error" not in r]) for name, rows in sets.items()},
        "sets": sets,
    }
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8", newline="\n")
    (out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(f"wrote {out / 'report.json'} and report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
