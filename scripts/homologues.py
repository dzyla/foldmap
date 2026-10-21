# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The homologue-pair campaign (plan 13a, spec 4.4): structure B drawn like a reference A (`--like`), judged by
H1-H5 (foldmap.criteria.like_checks) beside B drawn on its own (needs network for the PDB entries; not run in CI).

    python scripts/homologues.py align --set docs/validation/homologues/tuning.json --out DIR

align aligns each pair's two chains (foldmap.align.align_chains, auto route; as deposited) and element-mates
them, without drawing: per pair the method, TM-score (normalised by A), RMSD, identity, aligned pairs, mates and
same_fold, to DIR/align.json and a table on stdout. (The comparison with US-align, Task 13a-3, waits for the
author's approval to install it on the remote machine; the task report says what it needs.)

    python scripts/homologues.py run --set docs/validation/homologues/tuning.json --out DIR [--jobs 4]
                                     [--pairs 2VWD_A-2ZB5_A,...]

run draws every pair (or the --pairs named, by stem A_CHAIN-B_CHAIN) with `pair`, a child process each (python -I,
PYTHONHASHSEED=0, OPENBLAS_NUM_THREADS=1), then writes the report as `report` does. `pair` draws B like A as
`foldmap check B --like A --like-chain CA --chains CB --assembly asu` does (each side as deposited: the set
names deposited chains), and B alone (its same chain), judges both (criteria.check_figure: the S- and A-checks)
and the pair (criteria.like_checks, with the pair's kind for H1), and writes DIR/<stem>.json and two pages side by
side (compare's page): DIR/<stem>-like.png (A | B drawn like A) and DIR/<stem>-alone.png (A | B drawn on its own).

    python scripts/homologues.py report --set docs/validation/homologues/tuning.json --out DIR [--md PATH]

report reads DIR/<stem>.json of every pair and writes the table (DIR/report.md, or --md): per pair H1-H5, H2's and
H3's values drawn alike and drawn independently, the checks B drawn like A loses, and the failures by check.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
CHECKS = ("H1", "H2", "H3", "H4", "H5")
PAIR_TIMEOUT_S = 1800  # a pair's child is killed after this many seconds


def read_set(path) -> list[dict]:
    """The pairs of a set file ({"pairs": [{a: {id, chain}, b: {id, chain}, kind, note}, ...]}), as written."""
    return json.loads(Path(path).read_text(encoding="utf-8"))["pairs"]


def stem(pair: dict) -> str:
    """A pair's file stem: A_CHAIN-B_CHAIN."""
    return f"{pair['a']['id']}_{pair['a']['chain']}-{pair['b']['id']}_{pair['b']['chain']}"


def _round(x) -> float | None:
    return None if x is None else float(round(float(x), 4)) + 0.0


# -------------------------------------------------------------------------------------------------------- align


def _align(pair: dict) -> dict:
    """The pair aligned and mated as `foldmap plot B --like A` aligns it (make_layout's _drawn_like): each chain laid
    out alone first, for its elements."""
    from foldmap.align import align_chains, element_mates
    from foldmap.cli import _resolve_structure, make_layout
    from foldmap.io import read_model

    sides = []
    for side in ("a", "b"):
        path = _resolve_structure(pair[side]["id"])
        model = read_model(path, "asu")
        lay, sses, bb = make_layout(path, chains=[pair[side]["chain"]], model=model, assembly="asu")
        sides.append((sses, bb))
    (sses_a, bb_a), (sses_b, bb_b) = sides
    aln = align_chains(bb_a, pair["a"]["chain"], bb_b, pair["b"]["chain"], "auto")
    mates = element_mates(sses_a, sses_b, aln)
    return {
        "pair": stem(pair),
        "kind": pair["kind"],
        "method": aln.method,
        "tm_score": _round(aln.tm_score),
        "rmsd": _round(aln.rmsd),
        "identity": _round(aln.identity),
        "aligned": len(aln.pairs),
        "length_a": len(aln.residues_a),
        "length_b": len(aln.residues_b),
        "mates": len(mates.mate),
        "pieces": len(mates.pieces),
        "extras": len(mates.extras),
        "only_a": len(mates.only_a),
        "same_fold": mates.same_fold,
    }


def cmd_align(args) -> int:
    rows = [_align(p) for p in read_set(args.set)]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "align.json").write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print("pair\tkind\tmethod\tTM\tRMSD\tidentity\taligned\tmates/pieces/extras/only_a\tsame_fold")
    for r in rows:
        print(
            f"{r['pair']}\t{r['kind']}\t{r['method']}\t{r['tm_score']:.4f}\t{r['rmsd']:.2f}\t{r['identity']:.3f}\t"
            f"{r['aligned']}/{r['length_a']}\t{r['mates']}/{r['pieces']}/{r['extras']}/{r['only_a']}\t{r['same_fold']}"
        )
    return 0


# --------------------------------------------------------------------------------------------------------- pair


def _summary(found: dict) -> dict:
    """check_figure's result in short: the S-checks with problems (their counts) and the criteria that fail."""
    return {
        "science": {k: len(v) for k, v in sorted(found["science"].items()) if v},
        "failed": sorted(k for k, v in found["criteria"].items() if v.get("pass") is False),
    }


def run_pair(pair: dict, out: Path) -> dict:
    """One pair drawn like and alone, judged, its pages written (module doc)."""
    from foldmap import cli
    from foldmap.compare import _page
    from foldmap.criteria import check_figure, jsonable, like_checks
    from foldmap.render import save
    from foldmap.spec import draw_spec

    a, b = pair["a"], pair["b"]
    args = cli.build_parser().parse_args(["check", b["id"], "--like", a["id"], "--like-chain", a["chain"],
                                          "--chains", b["chain"], "--assembly", "asu"])  # fmt: skip
    args.structure = cli._resolve_structure(args.structure)
    spec = cli._spec_of(args)
    t0 = time.perf_counter()
    got: dict = {}
    drawn = draw_spec(spec, like_out=got)
    t_like = time.perf_counter() - t0
    alone_spec = replace(spec, like=None)
    alone = draw_spec(alone_spec)
    found_like = check_figure(drawn, spec, spec.structure)
    found_alone = check_figure(alone, alone_spec, spec.structure)
    ref = got.get("reference")
    note = drawn.layout.like_note
    h = like_checks(
        ref[0] if ref else None,
        drawn.layout,
        got.get("alignment"),
        got.get("mates"),
        note,
        got.get("like") is not None,
        alone.layout,
        found_like,
        found_alone,
        kind=pair["kind"],
    )
    aln, mates = got.get("alignment"), got.get("mates")
    s = stem(pair)
    if ref:
        lay_a, sses_a, bb_a, model_a, chain_a = ref
        look = spec.look
        for tag, other, title in (("like", drawn, "drawn like it"), ("alone", alone, "drawn on its own")):
            panels = (
                (lay_a, sses_a, bb_a, model_a, f"{a['id']} chain {chain_a}"),
                (other.layout, other.sses, other.bb, other.model, f"{b['id']} chain {b['chain']}, {title}"),
            )
            fig, _ = _page(panels, note if tag == "like" else "", look)
            save(fig, out / f"{s}-{tag}.png")
    result = {
        "pair": s,
        "kind": pair["kind"],
        "note": note,
        "alignment": None
        if aln is None
        else {
            "method": aln.method,
            "tm_score": _round(aln.tm_score),
            "rmsd": _round(aln.rmsd),
            "identity": _round(aln.identity),
            "aligned": len(aln.pairs),
            "length_a": len(aln.residues_a),
            "length_b": len(aln.residues_b),
        },
        "mates": None
        if mates is None
        else {
            "mates": len(mates.mate),
            "pieces": len(mates.pieces),
            "extras": len(mates.extras),
            "only_a": len(mates.only_a),
            "same_fold": mates.same_fold,
        },
        "like": jsonable(h),
        "checks_like": _summary(found_like),
        "checks_alone": _summary(found_alone),
        "seconds_like": round(t_like, 1),
    }
    (out / f"{s}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return result


def cmd_pair(args) -> int:
    pair = {
        "a": dict(zip(("id", "chain"), args.a.split(":"))),
        "b": dict(zip(("id", "chain"), args.b.split(":"))),
        "kind": args.kind,
    }
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    r = run_pair(pair, out)
    print(f"{r['pair']}: " + ", ".join(f"{k} {_word(r['like'][k])}" for k in CHECKS))
    return 0


# ---------------------------------------------------------------------------------------------------------- run


def cmd_run(args) -> int:
    pairs = read_set(args.set)
    if args.pairs:
        wanted = set(args.pairs.split(","))
        pairs = [p for p in pairs if stem(p) in wanted]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1"}

    def one(p: dict) -> tuple[str, int]:
        s = stem(p)
        cmd = [sys.executable, "-I", str(SCRIPT), "pair", f"{p['a']['id']}:{p['a']['chain']}",
               f"{p['b']['id']}:{p['b']['chain']}", "--kind", p["kind"], "--out", str(out)]  # fmt: skip
        with open(out / f"{s}.log", "w", encoding="utf-8") as log:
            try:
                done = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, timeout=PAIR_TIMEOUT_S)
                code = done.returncode
            except subprocess.TimeoutExpired:
                code = -9
        print(f"{s}: exit {code}", flush=True)
        return s, code

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        codes = dict(pool.map(one, pairs))
    write_report(read_set(args.set), out, Path(args.md) if args.md else out / "report.md")
    return 0 if all(c == 0 for c in codes.values()) else 1


# ------------------------------------------------------------------------------------------------------- report


def _word(c: dict | None) -> str:
    if not c or c.get("pass") is None:
        return "–"
    return "pass" if c["pass"] else "**fail**"


def _place(v: dict | None) -> str:
    return "–" if not v else f"{v['median']:.2f} / {v['p90']:.2f}"


def _share(v: dict | None) -> str:
    return "–" if not v else f"{v['share']:.3f} ({v['pairs']})"


def write_report(pairs: list[dict], out: Path, md: Path) -> Path:
    """The table of every pair's result (module doc)."""
    rows, failures, missing = [], {k: [] for k in CHECKS}, []
    for p in pairs:
        f = out / f"{stem(p)}.json"
        if not f.is_file():
            missing.append(stem(p))
            continue
        r = json.loads(f.read_text(encoding="utf-8"))
        h, aln, m = r["like"], r["alignment"] or {}, r["mates"] or {}
        for k in CHECKS:
            if h[k].get("pass") is False:
                failures[k].append(f"{r['pair']} ({r['kind']}): {h[k]['detail']}")
        lost = (h["H5"].get("value") or {}).get("lost") or []
        rows.append(
            f"| {r['pair']} | {r['kind']} | {aln.get('method', '–')} | {aln.get('tm_score', 0):.3f} | "
            f"{m.get('mates', 0)}/{m.get('pieces', 0)}/{m.get('extras', 0)}/{m.get('only_a', 0)} | "
            + " | ".join(_word(h[k]) for k in CHECKS)
            + f" | {_place(h['H2'].get('value'))} | {_place(h['H2'].get('baseline'))} | "
            f"{_share(h['H3'].get('value'))} | {_share(h['H3'].get('baseline'))} | {', '.join(lost) or '–'} |"
        )
    lines = [
        "# Homologue pairs: drawn like the reference (H1-H5)",
        "",
        "TM: foldmap's TM-score normalised by A. Elements: mates / pieces / extras (B only) / only in A. H2: median / "
        "90th percentile of the mates' centre deviations (page units) after the best rigid 2D motion; H3: share of "
        "mate pairs whose page relation agrees (pairs judged). 'alike': B drawn like A; 'indep.': A and B drawn "
        "independently (the baseline). Lost: checks B drawn alone passes and B drawn like A does not (H5).",
        "",
        "| pair | kind | route | TM | elements | H1 | H2 | H3 | H4 | H5 | H2 alike | H2 indep. | H3 alike | "
        "H3 indep. | lost |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        *rows,
        "",
        "## Failures by check",
        "",
    ]
    for k in CHECKS:
        lines.append(f"- **{k}**: " + ("none" if not failures[k] else ""))
        lines += [f"  - {x}" for x in failures[k]]
    if missing:
        lines += ["", f"No result (the pair's child failed; see its log): {', '.join(missing)}"]
    md.parent.mkdir(parents=True, exist_ok=True)
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    return md


def cmd_report(args) -> int:
    out = Path(args.out)
    write_report(read_set(args.set), out, Path(args.md) if args.md else out / "report.md")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    al = sub.add_parser("align", help="align and mate every pair, no drawing")
    al.add_argument("--set", required=True)
    al.add_argument("--out", required=True)
    pr = sub.add_parser("pair", help="one pair (run's child)")
    pr.add_argument("a", help="ID:CHAIN of the reference")
    pr.add_argument("b", help="ID:CHAIN drawn like it")
    pr.add_argument("--kind", choices=("positive", "negative"), required=True)
    pr.add_argument("--out", required=True)
    rn = sub.add_parser("run", help="every pair of a set, then the report")
    rn.add_argument("--set", required=True)
    rn.add_argument("--out", required=True)
    rn.add_argument("--jobs", type=int, default=4)
    rn.add_argument("--pairs", default="", help="only these stems (A_CHAIN-B_CHAIN), comma-separated")
    rn.add_argument("--md", default="")
    rp = sub.add_parser("report", help="the table of a run")
    rp.add_argument("--set", required=True)
    rp.add_argument("--out", required=True)
    rp.add_argument("--md", default="")
    args = ap.parse_args(argv)
    return {"align": cmd_align, "pair": cmd_pair, "run": cmd_run, "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
