# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The same figure on any platform, CPU and hash seed: draw a fixed set of figures with Liberation Sans forced, record
each figure's layout signature and SVG hash, and compare runs.

    python scripts/paper/determinism.py --out RUN.json [--label NAME] [--svg-dir DIR] [--any-font] [STRUCTURE ...]
    python scripts/paper/determinism.py --configurations DIR [--any-font] [STRUCTURE ...]
    python scripts/paper/determinism.py --compare [--systems N] RUN.json|DIR|GLOB ...

Draw (--out): the default view (`foldmap plot STRUCTURE -o x.svg`, in this process, cli.draw_spec wrapped to keep the
drawn layout) of every committed test structure (tests/data/*.cif, sorted) unless structures are given. Each figure
gets
  signature          bench.signature(layout) as JSON text (sorted keys), as random100's A12/A13 compare it
  signature_sha256   sha256 of that text
  svg_sha256         sha256 of the SVG foldmap wrote (kept in --svg-dir if given)
and the run records fonts.font_report() (the font resolved after fonts.force_liberation_sans(), asserted per figure
with fonts.assert_figure_font; with the font file's revision), the platform, the Python, numpy, scipy, networkx,
matplotlib and gemmi versions, PYTHONHASHSEED and the CPU-feature environment. Exits 1 when a figure was not drawn.
--any-font draws with whatever font resolves (for comparing runs on one machine without Liberation Sans; never for
the cross-platform run).

Configurations (--configurations): one draw per CONFIGURATIONS entry, each in its own process, written to
DIR/<system>-<name>.json, then compared. On x86-64 Linux this is the paper's Linux leg: PYTHONHASHSEED 0 and 1, and 0
with NumPy's AVX dispatch off and OpenBLAS's oldest x86-64 kernels.

Compare (--compare): the runs are the JSON files given, found under a directory given, or matched by a pattern given
(none found: "0 run(s)"). Exits 1 unless every run drew the same figures with the same signature TEXT and SVG hash
(text, not parsed floats: 0.0 == -0.0 would hide a difference), and (--systems N) on at least N platform.system()
values. Runs whose package versions or font (family, revision) differ are "not comparable": that is the verdict,
and a figure difference is listed as "(versions differ)", not taken for the platform's. CI's determinism workflow
(.github/workflows/determinism.yml, pinned by determinism-constraints.txt) draws on Linux, macOS and Windows and
compares; the compare needs nothing but the standard library.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

NO_AVX = {"NPY_DISABLE_CPU_FEATURES": "AVX512F AVX512_SKX AVX2 FMA3", "OPENBLAS_CORETYPE": "Prescott"}
CONFIGURATIONS = {
    "seed0": {"PYTHONHASHSEED": "0"},
    "seed1": {"PYTHONHASHSEED": "1"},
    "seed0-noavx": {"PYTHONHASHSEED": "0", **NO_AVX},
}
VERSIONS = ("python", "numpy", "scipy", "networkx", "matplotlib", "gemmi")  # differing: not comparable (or the font)
TIMEOUT = 1800  # s for one configuration's draw (~10 s on zyla-lab)


def _fonts():
    spec = importlib.util.spec_from_file_location("paper_fonts", HERE / "fonts.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def liberation_sans_installed() -> bool:
    """Liberation Sans is what a font search for it resolves to here."""
    from matplotlib.font_manager import findfont
    from matplotlib.ft2font import FT2Font

    font = _fonts().FONT
    try:
        path = findfont(font, fallback_to_default=False)
    except ValueError:
        return False
    return FT2Font(path).family_name == font


def _jsonable(o):
    return o.tolist() if hasattr(o, "tolist") else str(o)


def signature_text(signature: dict) -> str:
    return json.dumps(signature, sort_keys=True, default=_jsonable)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def figure_row(structure: str, code: int, signature: str | None, svg: bytes | None) -> dict:
    """One figure of a run: its signature text and both hashes, or only the exit code when it was not drawn."""
    row = {"structure": structure, "exit": code}
    if signature is not None and svg is not None:
        row |= {"signature": signature, "signature_sha256": _sha256(signature.encode()), "svg_sha256": _sha256(svg)}
    return row


def _drawn(row: dict | None) -> bool:
    return row is not None and row["exit"] == 0 and "svg_sha256" in row


def draw(paths: list[Path], label: str = "", svg_dir: Path | None = None, force_font: bool = True) -> dict:
    """Draw each structure's default figure in this process; returns the run's record."""
    fonts = _fonts()
    font_file = fonts.force_liberation_sans() if force_font else None  # before anything is drawn

    import gemmi
    import matplotlib
    import matplotlib.pyplot as plt
    import networkx
    import numpy as np
    import scipy

    sys.path.insert(0, str(REPO / "scripts"))
    import bench

    from foldmap import cli

    caught: dict = {}
    real = cli.draw_spec

    def keep(spec):
        caught["drawn"] = real(spec)
        return caught["drawn"]

    cli.draw_spec = keep
    rows = []
    try:
        with tempfile.TemporaryDirectory() as tmp:
            for p in paths:
                caught.clear()
                svg = (svg_dir or Path(tmp)) / f"{p.stem}.svg"
                svg.parent.mkdir(parents=True, exist_ok=True)
                code = cli.main(["plot", str(p), "-o", str(svg)])
                drawn = caught.get("drawn")
                if drawn is not None and svg.is_file():
                    if force_font:
                        fonts.assert_figure_font(drawn.figure)
                    row = figure_row(p.name, code, signature_text(bench.signature(drawn.layout)), svg.read_bytes())
                else:
                    row = figure_row(p.name, code, None, None)
                rows.append(row)
                print(p.name, row.get("signature_sha256", "")[:12], row.get("svg_sha256", "")[:12], flush=True)
                plt.close("all")
    finally:
        cli.draw_spec = real
    return {
        "label": label,
        "system": platform.system(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "networkx": networkx.__version__,
        "matplotlib": matplotlib.__version__,
        "gemmi": gemmi.__version__,
        "font": fonts.font_report() | {"forced_file": font_file},
        "env": {k: os.environ.get(k) for k in ("PYTHONHASHSEED", *NO_AVX, "OPENBLAS_NUM_THREADS")},
        "figures": rows,
    }


def configuration_env(extra: dict) -> dict:
    """This process's environment with one configuration's variables (and no stray CPU-feature switches)."""
    env = {k: v for k, v in os.environ.items() if k not in NO_AVX}
    return env | {"OPENBLAS_NUM_THREADS": "1", **extra}


def run_configurations(out_dir: Path, structures: list[Path], force_font: bool = True) -> list[dict]:
    """Draw once per CONFIGURATIONS entry, each in a fresh process; returns the runs (also written to out_dir)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    runs = []
    for name, extra in CONFIGURATIONS.items():
        label = f"{platform.system().lower()}-{name}"
        out = out_dir / f"{label}.json"
        cmd = [sys.executable, str(Path(__file__).resolve()), "--out", str(out), "--label", label]
        cmd += [] if force_font else ["--any-font"]
        subprocess.run([*cmd, *map(str, structures)], env=configuration_env(extra), check=True, timeout=TIMEOUT)
        runs.append(json.loads(out.read_text(encoding="utf-8")))
    return runs


def _first_difference(a: str, b: str) -> str:
    i = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
    return f"first difference at character {i}: {a[max(0, i - 30) : i + 30]!r} vs {b[max(0, i - 30) : i + 30]!r}"


def _font(run: dict) -> str:
    return f"{run['font']['font_family']} {run['font'].get('font_revision', '?')}"


def not_comparable(runs: list[dict]) -> list[str]:
    """What differs between the runs other than the figures: package versions, the font (family and revision)."""
    labels = [r.get("label") or f"run {i}" for i, r in enumerate(runs)]
    out = []
    for key in (*VERSIONS, "font"):
        got = [_font(r) if key == "font" else r.get(key) for r in runs]
        if len(set(got)) > 1:
            out.append(f"{key} differs: " + ", ".join(f"{lab} {v}" for lab, v in zip(labels, got)))
    return out


def compare(runs: list[dict], systems: int = 0) -> list[str]:
    """Every difference between the runs ([] when they drew the same figures identically): first what makes them
    not comparable ("not comparable: ..."), then the figures' (marked "(versions differ)" when there is any)."""
    labels = [r.get("label") or f"run {i}" for i, r in enumerate(runs)]
    problems = []
    if len(runs) < 2:
        problems.append(f"{len(runs)} run(s): nothing to compare")
    found = sorted({r["system"] for r in runs})
    if len(found) < systems:
        problems.append(f"{len(found)} platform(s) ({', '.join(found)}), expected {systems}")
    apart = not_comparable(runs)
    problems += [f"not comparable: {why}" for why in apart]
    mark = " (versions differ)" if apart else ""
    names = [[f["structure"] for f in r["figures"]] for r in runs]
    for lab, n in zip(labels[1:], names[1:]):
        if n != names[0]:
            problems.append(f"{lab} drew {n}, {labels[0]} drew {names[0]}")
    for name in sorted(set().union(*names)) if runs else []:
        rows = {lab: next((f for f in r["figures"] if f["structure"] == name), None) for lab, r in zip(labels, runs)}
        missing = [lab for lab, f in rows.items() if not _drawn(f)]
        if missing:
            problems.append(f"{name}: not drawn in {', '.join(missing)}")
            continue
        for key in ("signature", "signature_sha256", "svg_sha256"):
            groups: dict = {}
            for lab, f in rows.items():
                groups.setdefault(f.get(key), []).append(lab)
            if len(groups) > 1:
                why = " vs ".join("[" + ", ".join(labs) + "]" for labs in groups.values())
                if key == "signature" and None not in groups:
                    a, b = list(groups)[:2]
                    why += "; " + _first_difference(a, b)
                problems.append(f"{name}: {key} differs{mark}: {why}")
    return problems


def find_runs(args: list[str]) -> list[Path]:
    """The run files: each argument a file, a directory (every *.json under it) or a pattern; nothing found, none."""
    out: list[Path] = []
    for a in args:
        p = Path(a)
        out += sorted(p.rglob("*.json")) if p.is_dir() else sorted(Path(m) for m in glob.glob(a))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--out", type=Path, help="draw, and write the run here")
    mode.add_argument("--configurations", type=Path, metavar="DIR", help="draw once per configuration, then compare")
    mode.add_argument("--compare", action="store_true", help="compare the runs given as paths")
    ap.add_argument("--label", default="")
    ap.add_argument("--svg-dir", type=Path, help="keep the SVGs here")
    ap.add_argument("--any-font", action="store_true", help="do not force Liberation Sans")
    ap.add_argument("--systems", type=int, default=0, help="--compare: fail with fewer platforms than this")
    ap.add_argument("paths", nargs="*", type=Path, help="structures to draw (default tests/data/*.cif), or runs")
    args = ap.parse_args(argv)
    if args.compare:
        runs = [json.loads(p.read_text(encoding="utf-8")) for p in find_runs(list(map(str, args.paths)))]
    else:
        structures = args.paths or sorted((REPO / "tests" / "data").glob("*.cif"))
        if args.out:
            doc = draw(structures, args.label, args.svg_dir, force_font=not args.any_font)
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
            return 0 if all(_drawn(f) for f in doc["figures"]) else 1
        runs = run_configurations(args.configurations, structures, force_font=not args.any_font)
    problems = compare(runs, args.systems)
    for line in problems:
        print(line)
    systems = sorted({r["system"] for r in runs})
    figures = len(runs[0]["figures"]) if runs else 0
    apart = [line.removeprefix("not comparable: ") for line in problems if line.startswith("not comparable: ")]
    verdict = f"not comparable: {'; '.join(apart)}" if apart else "differ" if problems else "identical"
    print(f"{len(runs)} run(s) on {', '.join(systems) or 'no platform'}: {figures} figures, {verdict}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
