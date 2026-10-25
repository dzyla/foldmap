# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Step 4 (zyla-lab: it draws a figure): default-view time against residue count over the three random-100 samples.

    python scripts/paper/runtime.py   (step 4 of scripts/paper/run_all.sh)

No new run: the times are the campaign's own (report.json: views.default.seconds, the default view's process wall time
that A2 gates; views.default.result.seconds_inside, the drawing inside it), the residue counts the entry facts
(facts.residues: polymer residues modelled in the file's first model, as random100.features counts them). The runs
drew ten entries at once on zyla-lab (--jobs 10), so the times include contention.

Writes build/tables/runtime.csv (one row per entry), build/runtime.json (power-law fit of time on residues, medians
per size bin) and build/runtime.png / runtime.pdf.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BUILD, TABLES, load_source, percentile  # noqa: E402

SAMPLES = {"tuned": "r100_tuned", "held": "r100_held", "third": "r100_third"}
COLOURS = {"tuned": "#0072B2", "held": "#E69F00", "third": "#009E73"}  # Okabe-Ito
LABELS = {"tuned": "tuned", "held": "second (held-out)", "third": "third"}
BINS = {"S": (0, 299), "M": (300, 1499), "L": (1500, 7999), "XL": (8000, 10**9)}


def rows_from_reports(reports: dict) -> list[dict]:
    rows = []
    for s, rep in reports.items():
        for e in rep["entries"]:
            v = e["views"]["default"]
            rows.append(
                {
                    "sample": s,
                    "id": e["id"],
                    "stratum": e["stratum"],
                    "bin": e["bin"],
                    "residues": e["facts"]["residues"],
                    "chains": e["facts"]["chains"],
                    "process_s": v["seconds"],
                    "drawing_s": (v.get("result") or {}).get("seconds_inside"),
                    "pass": e["pass"],
                }
            )
    return rows


def fit(rows: list[dict], key: str) -> dict:
    """Least squares of log10(time) on log10(residues): time ~ c * residues^k."""
    pts = [(math.log10(r["residues"]), math.log10(r[key])) for r in rows if r["residues"] > 0 and (r[key] or 0) > 0]
    n = len(pts)
    mx = sum(x for x, _ in pts) / n
    my = sum(y for _, y in pts) / n
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    k = sxy / sxx
    a = my - k * mx
    ss_res = sum((y - (a + k * x)) ** 2 for x, y in pts)
    ss_tot = sum((y - my) ** 2 for _, y in pts)
    return {"n": n, "exponent": k, "c_seconds": 10**a, "r2": 1 - ss_res / ss_tot}


def main() -> int:
    reports = {s: load_source(f"{p}_report") for s, p in SAMPLES.items()}
    rows = rows_from_reports(reports)
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "runtime.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    per_bin = {}
    for b, (lo, hi) in BINS.items():
        sel = [r for r in rows if lo <= r["residues"] <= hi]
        if sel:
            per_bin[b] = {
                "residues": [lo, hi],
                "n": len(sel),
                "process_median_s": percentile([r["process_s"] for r in sel], 50),
                "drawing_median_s": percentile([r["drawing_s"] for r in sel], 50),
                "process_max_s": max(r["process_s"] for r in sel),
            }
    out = {
        "process": fit(rows, "process_s"),
        "drawing": fit(rows, "drawing_s"),
        "per_residue_bin": per_bin,
        "residue_range": [min(r["residues"] for r in rows), max(r["residues"] for r in rows)],
    }
    (BUILD / "runtime.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharex=True, sharey=True, constrained_layout=True)
    for ax, key, title in (
        (axes[0], "process_s", "a  default view, whole process (A2)"),
        (axes[1], "drawing_s", "b  drawing only"),
    ):
        for s in SAMPLES:
            sel = [r for r in rows if r["sample"] == s]
            ax.scatter(
                [r["residues"] for r in sel],
                [r[key] for r in sel],
                s=9,
                lw=0.3,
                edgecolor="white",
                color=COLOURS[s],
                label=f"{LABELS[s]} (n={len(sel)})",
                alpha=0.9,
            )
        f_ = out["process" if key == "process_s" else "drawing"]
        xs = [10 ** (1 + i * 0.05) for i in range(int((math.log10(out["residue_range"][1] * 1.2) - 1) / 0.05) + 1)]
        ax.plot(
            xs,
            [f_["c_seconds"] * x ** f_["exponent"] for x in xs],
            color="#444444",
            lw=0.8,
            ls="-",
            label=f"fit: t ∝ n^{f_['exponent']:.2f} (R² {f_['r2']:.2f})",
        )
        if key == "process_s":
            ax.axhline(60, color="#999999", lw=0.6, ls="--")
            ax.text(12, 66, "A2 limit 60 s (180 s above 20 chains)", fontsize=7, color="#666666")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("polymer residues in the file (model 1)")
        ax.set_title(title, loc="left", fontsize=8)
        ax.grid(True, which="major", lw=0.3, color="#dddddd")
        ax.legend(fontsize=7, frameon=False, loc="lower right")
    axes[0].set_ylabel("wall time (s)")
    for ext in ("png", "pdf"):
        fig.savefig(BUILD / f"runtime.{ext}", dpi=300)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
