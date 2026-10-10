# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Time the layout of the test structures (and a synthetic 32-chain complex) and compare layouts between runs.

python scripts/bench.py                          # every tests/data structure plus the 32-chain complex
python scripts/bench.py 6ZYA 8UTF                # some of them
python scripts/bench.py --record before.json     # keep each layout (positions, bounds)
python scripts/bench.py --check before.json --exact      # fail unless every layout is bit-identical
python scripts/bench.py --check before.json --area 1.01  # fail if a figure's area grew by more than 1 %

A signature holds every placed element, ghost and duplex with all its geometry, plus sheets, domains, unjoined
pairs and nt_rows. Loop points are derived from the layout at render time, so identical layouts give identical loops.
"""

from __future__ import annotations

import argparse
import json
import string
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import gemmi

import foldmap
from foldmap import cli, io, layout
from foldmap.cli import make_layout

DATA = Path(__file__).resolve().parent.parent / "tests" / "data"
PHASES = {  # timed wherever they are looked up
    "read_model": (io, cli),
    "pair_nucleotides": (io,),
    "_ligands": (io,),
    "_attract": (layout,),
    "_separate": (layout,),
}
SPENT: dict[str, float] = {}
OPTIONS = {  # every make_layout option the figures depend on, passed explicitly so a changed default cannot go unseen
    "mode": "projected",
    "rotate": 0.0,
    "flip_v": False,
    "dna": True,
    "symmetry": "auto",
    "assembly": "auto",
    "membrane": "off",
    "focus": "auto",
}


def big_complex(path: Path, copies: int = 8) -> Path:
    """2HHB tiled `copies` times, 70 Å apart: 4 * copies chains, no assembly records."""
    st = gemmi.read_structure(str(DATA / "2HHB.cif"))
    st.remove_ligands_and_waters()
    model = st[0]
    base = [c.clone() for c in model]
    names = iter(a + b for a in string.ascii_uppercase for b in string.ascii_uppercase)
    for c in list(model):
        model.remove_chain(c.name)
    for k in range(copies):
        for c in base:
            new = c.clone()
            new.name = next(names)
            for r in new:
                for a in r:
                    a.pos = gemmi.Position(a.pos.x + 70 * (k % 4), a.pos.y + 70 * (k // 4), a.pos.z)
            model.add_chain(new)
    st.assemblies.clear()
    st.setup_entities()
    st.make_mmcif_document().write_file(str(path))
    return path


def _timed(name: str, f):
    def run(*a, **k):
        t = time.perf_counter()
        try:
            return f(*a, **k)
        finally:
            SPENT[name] = SPENT.get(name, 0.0) + time.perf_counter() - t

    return run


def _meta() -> dict:
    try:
        done = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=DATA, check=True)
        commit = done.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    return {"options": OPTIONS, "commit": commit, "version": foldmap.__version__}


def signature(lay) -> dict:
    sig = {k: [p.cx, p.cy, p.angle, p.length, p.width, p.ends, p.ghost, p.label] for k, p in sorted(lay.placed.items())}
    sig["_ghosts"] = [[g.cx, g.cy, g.angle, g.length, g.width, g.ends, g.ghost, g.label] for g in lay.ghosts]
    sig["_dna"] = [
        [
            d.id,
            d.strands,
            d.cx,
            d.cy,
            d.angle,
            d.length,
            d.width,
            d.ghost,
            sorted([list(k), v] for k, v in d.axial.items()),
        ]
        for d in lay.dna
    ]
    sig["_bounds"] = list(lay.bounds)
    sig["_sheets"] = [list(x) for x in lay.sheets]
    sig["_domains"] = [[n, list(ids)] for n, ids in lay.domains]
    sig["_unjoined"] = sorted(sorted(pair) for pair in lay.unjoined)
    sig["_nt_rows"] = lay.nt_rows
    return sig


def _area(bounds) -> float:
    return (bounds[2] - bounds[0]) * (bounds[3] - bounds[1])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*", help="tests/data entries (e.g. 6ZYA) or 'big32'; default: all")
    ap.add_argument("--record", help="write every layout to this JSON file")
    ap.add_argument("--check", help="compare with layouts recorded earlier")
    ap.add_argument("--exact", action="store_true", help="with --check: fail unless layouts are bit-identical")
    ap.add_argument("--area", type=float, help="with --check: fail if a figure's area grew by more than this factor")
    args = ap.parse_args(argv)
    if args.check and not (args.exact or args.area):
        ap.error("--check needs --exact or --area, otherwise it cannot fail")
    with tempfile.TemporaryDirectory() as tmp:
        return _run(args, Path(tmp))


def _run(args, tmp: Path) -> int:
    names = args.names or [p.stem for p in sorted(DATA.glob("*.cif"))] + ["big32"]
    files = [big_complex(tmp / "big32.cif") if n == "big32" else DATA / f"{n}.cif" for n in names]
    for name, mods in PHASES.items():
        wrapped = _timed(name, getattr(mods[0], name))
        for mod in mods:
            setattr(mod, name, wrapped)
    ref = json.loads(Path(args.check).read_text(encoding="utf-8")) if args.check else {}
    meta = _meta()
    out, failed = {"_meta": meta}, []
    if ref and ref.get("_meta", {}).get("options") != OPTIONS:
        print("failed: make_layout options differ from the recorded ones", file=sys.stderr)
        return 1
    for f in files:
        SPENT.clear()
        t = time.perf_counter()
        lay, _, _ = make_layout(f, **OPTIONS)
        total = time.perf_counter() - t
        sig = out[f.stem] = signature(lay)
        b = sig["_bounds"]
        note = ""
        if args.check and f.stem not in ref:
            note = "MISSING from baseline"
            failed.append(f.stem)
        elif f.stem in ref:
            r = ref[f.stem]
            same = r == sig
            grew = _area(b) / max(_area(r["_bounds"]), 1e-9)
            note = "identical" if same else f"changed, area x{grew:.3f}"
            if (args.exact and not same) or (args.area and grew > args.area):
                failed.append(f.stem)
        phases = " ".join(f"{k}={v:.2f}" for k, v in SPENT.items() if v >= 0.005)
        size = f"{b[2] - b[0]:6.1f} x {b[3] - b[1]:5.1f}"
        print(f"{f.stem:10s} {len(lay.placed):4d} el {total:7.2f} s  {size}  {phases}  {note}")
    if ref:
        unmeasured = sorted(set(ref) - set(out))
        if unmeasured:
            print("not measured:", " ".join(unmeasured), file=sys.stderr)
            if not args.names:  # a full run must cover the whole baseline
                failed.extend(unmeasured)
    if args.record:
        Path(args.record).write_text(json.dumps(out), encoding="utf-8")
    if failed:
        print("failed:", " ".join(failed), file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
