"""Render a random sample of PDB entries with the default style and check every figure.

    python scripts/sweep.py --sample 60 --seed 20261009 --out sweep/

For each entry: time, errors, element overlaps (must be none), loops that fell back to plain curves, loops
crossing labels; large entries (> --max-residues) are drawn one or two chains at a time. Writes report.json,
report.md and contact sheets (sheet-*.png) of the figures for a visual check.
"""

from __future__ import annotations

import argparse
import json
import random
import time
import traceback
import urllib.request
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

QUERY = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id exptl{method}
  rcsb_entry_info{ deposited_polymer_monomer_count polymer_entity_count_protein } struct{title} } }"""


def sample(n: int, seed: int) -> list[dict]:
    ids = json.load(urllib.request.urlopen("https://data.rcsb.org/rest/v1/holdings/current/entry_ids", timeout=60))
    random.seed(seed)
    pick = random.sample(ids, n * 2)
    body = json.dumps({"query": QUERY, "variables": {"ids": pick}}).encode()
    req = urllib.request.Request(
        "https://data.rcsb.org/graphql", data=body, headers={"Content-Type": "application/json"}
    )
    entries = [e for e in json.load(urllib.request.urlopen(req, timeout=60))["data"]["entries"] if e]
    out = []
    for e in entries:
        info = e.get("rcsb_entry_info") or {}
        if (info.get("polymer_entity_count_protein") or 0) > 0:
            out.append(
                {
                    "id": e["rcsb_id"],
                    "method": (e.get("exptl") or [{}])[0].get("method"),
                    "residues": info.get("deposited_polymer_monomer_count") or 0,
                    "title": (e.get("struct") or {}).get("title", "")[:80],
                }
            )
    return out[:n]


def _overlaps(layout) -> int:
    rects = [p.rect for p in layout.placed.values()]
    bad = 0
    for i in range(len(rects)):
        for j in range(i + 1, len(rects)):
            a, b = rects[i], rects[j]
            if min(a[2], b[2]) - max(a[0], b[0]) > 0.05 and min(a[3], b[3]) - max(a[1], b[1]) > 0.05:
                bad += 1
    return bad


def check(entry: dict, out_dir: str, max_residues: int) -> dict:
    from foldmap.cli import make_layout
    from foldmap.fetch import fetch
    from foldmap.io import load_backbone
    from foldmap.render import draw, save
    from foldmap.route import route_loops

    row = dict(entry)
    start = time.time()
    try:
        path = fetch(entry["id"])
        opts: dict = {}
        if entry["residues"] > max_residues:  # a large entry: one or two chains of the deposited model
            bb = load_backbone(path, "asu")
            chains, total = [], 0
            for c in dict.fromkeys(l.chain for l in bb.labels):
                chains.append(c)
                total += sum(l.chain == c for l in bb.labels)
                if total >= 300 or len(chains) == 2:
                    break
            opts = {"assembly": "asu", "chains": chains}
            row["drawn"] = f"chains {','.join(chains)}"
        layout, sses, bb = make_layout(path, **opts)
        loops = route_loops(layout, sses, bb)
        fig = draw(layout, loops, sses, entry["id"])
        save(fig, Path(out_dir) / f"{entry['id']}.png", dpi=60)
        row.update(
            ok=True,
            residues_drawn=len(bb),
            chains=len(set(l.chain for l in bb.labels)),
            elements=len(layout.placed),
            overlaps=_overlaps(layout),
            fallback=sum(l.fallback for l in loops),
            loops=len(loops),
            crossing_labels=sum(bool(getattr(l, "crosses_labels", False)) for l in loops),
        )
    except Exception as err:  # noqa: BLE001 - the report records every failure
        row.update(ok=False, error=f"{type(err).__name__}: {err}", trace=traceback.format_exc()[-1500:])
    row["seconds"] = round(time.time() - start, 1)
    return row


def sheets(rows: list[dict], out: Path, per: int = 12) -> None:
    from matplotlib.figure import Figure
    from matplotlib.image import imread

    good = [r for r in rows if r.get("ok")]
    for k in range(0, len(good), per):
        chunk = good[k : k + per]
        fig = Figure(figsize=(16, 12))
        for n, r in enumerate(chunk):
            ax = fig.add_subplot(3, 4, n + 1)
            ax.imshow(imread(out / f"{r['id']}.png"))
            ax.set_title(f"{r['id']} · {r['elements']} el · fb {r['fallback']} · ov {r['overlaps']}", fontsize=9)
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(out / f"sheet-{k // per + 1}.png", dpi=110)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20261009)
    ap.add_argument("--out", default="sweep")
    ap.add_argument("--max-residues", type=int, default=2500)
    ap.add_argument("--timeout", type=int, default=240)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    entries = sample(args.sample, args.seed)
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [(e, pool.submit(check, e, str(out), args.max_residues)) for e in entries]
        for e, f in futures:
            try:
                rows.append(f.result(timeout=args.timeout))
            except TimeoutError:
                rows.append({**e, "ok": False, "error": f"timeout after {args.timeout} s"})
            r = rows[-1]
            print(
                f"{r['id']}: {'ok' if r.get('ok') else 'FAIL'} {r.get('seconds', '')}s {r.get('error', '')[:100]}",
                flush=True,
            )
    (out / "report.json").write_text(json.dumps(rows, indent=1))
    ok = [r for r in rows if r.get("ok")]
    lines = [
        f"# Sweep: {len(ok)}/{len(rows)} rendered",
        "",
        "| id | method | residues | drawn | elements | overlaps | fallback loops | s | note |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['id']} | {r.get('method')} | {r.get('residues')} | {r.get('drawn', 'all')} | "
            f"{r.get('elements', '')} | {r.get('overlaps', '')} | {r.get('fallback', '')} | "
            f"{r.get('seconds', '')} | {r.get('error', '')[:80]} |"
        )
    (out / "report.md").write_text("\n".join(lines) + "\n")
    sheets(rows, out)


if __name__ == "__main__":
    main()
