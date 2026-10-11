# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Check foldmap's nucleic-acid figures against the structures they come from (needs network; not run in CI).

    python scripts/na_validate.py                       # the whole benchmark: downloads into ~/.cache/foldmap
    python scripts/na_validate.py --ids 30QM,1EHZ       # some entries
    python scripts/na_validate.py --files tests/data/na/*.cif.gz --out na-report-fixtures

RNApolis is a development tool, never a dependency: install it in a throwaway venv and point at its Python,
    python -m venv /tmp/na-venv && /tmp/na-venv/bin/pip install rnapolis gemmi
    python scripts/na_validate.py --rnapolis-python /tmp/na-venv/bin/python   (or set NA_RNAPOLIS_PYTHON)
Without it the pairing columns are empty and the other checks still run.

For every entry it measures
  pairing   foldmap's canonical and wobble pairs against RNApolis (run from its own venv): precision, recall, F1;
            stems against stems grouped the same way from the RNApolis pairs
  shape     drawn stem centres against the 3D centres: Spearman rho of pairwise distances, similarity-Procrustes
            RMSD (no reflection) normalised by the spread; stem angle against the projected 3D axis; bends between
            the pieces of a bent duplex; whether any stem's arms are drawn mirrored
  contacts  whether the protein elements touching a stem are drawn nearest it (AUC, precision at k); beads drawn
            for every touched nucleotide, in the touching chain's colour
  render    overlapping boxes, text off the page, overlapping base letters, routing fallbacks, run time
and writes a contact sheet per entry (the figure, the drawn layout with stems coloured, the 3D structure projected
the same way with the same colours) plus report.md and report.json with pass/fail against THRESHOLDS.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from foldmap.cli import make_layout  # noqa: E402
from foldmap.fetch import fetch  # noqa: E402
from foldmap.io import _nucleic_strands, read_model  # noqa: E402
from foldmap.layout import SCALE, stem_pieces  # noqa: E402
from foldmap.nacheck import THRESHOLDS, contacts, judge, render_checks, shape  # noqa: E402
from foldmap.nucleic import BasePair, pair_nucleotides, secondary_structure, stem_axis  # noqa: E402
from foldmap.palette import chain_colors  # noqa: E402
from foldmap.render import draw, save  # noqa: E402
from foldmap.route import route_loops  # noqa: E402

# id -> (category, why it is in the set)
BENCHMARK = {
    "1EHZ": ("tRNA", "tRNA-Phe alone: cloverleaf, ten kinds of modified nucleotide, G4·U69"),
    "1QTQ": ("tRNA", "GlnRS-tRNA: a tRNA held by its synthetase, acceptor end unwound into the active site"),
    "4LVV": ("riboswitch", "THF riboswitch: three-way junction with a pseudoknot"),
    "2GIS": ("riboswitch", "SAM-I riboswitch: four-way junction, kink-turn, pseudoknot"),
    "3D0U": ("riboswitch", "lysine riboswitch: five-way junction"),
    "2GDI": ("riboswitch", "TPP riboswitch: two parallel stacked arms"),
    "1Y26": ("riboswitch", "adenine riboswitch: small three-way junction, loop-loop kissing"),
    "1GID": ("ribozyme", "group I intron P4-P6 domain: long coaxial stacks, A-rich bulge, tetraloop receptor"),
    "3IGI": ("ribozyme", "group II intron: about 390 nt, many junctions, 3.1 Å"),
    "4P9R": ("ribozyme", "lariat-capping ribozyme: group I-like core, multiple pseudoknots"),
    "2OEU": ("ribozyme", "full-length hammerhead ribozyme: three stems, tertiary loop contact"),
    "2NZ4": ("ribozyme", "glmS ribozyme with U1A crystallisation module: double pseudoknot"),
    "4OO8": ("CRISPR", "Cas9-sgRNA-target DNA: RNA-DNA heteroduplex, repeat:antirepeat, stem loops"),
    "5B43": ("CRISPR", "Cas12a-crRNA-target DNA: pseudoknot handle, heteroduplex, PAM duplex"),
    "5W1H": ("CRISPR", "Cas13a-crRNA: repeat stem-loop held in a protein groove"),
    "5F9R": ("CRISPR", "Cas9 R-loop: target strand in a heteroduplex, non-target strand displaced"),
    "30QM": ("RNP", "DRT10 RT-ncRNA: 121-nt ncRNA, six stems between two protomers (the reported case)"),
    "1URN": ("RNP", "U1A-RNA hairpin: protein on a single-stranded loop"),
    "4PJO": ("RNP", "minimal U1 snRNP at 3.3 Å: U1 snRNA with Sm ring, U1-70K and U1-C"),
    "1HQ1": ("RNP", "SRP M domain-4.5S RNA: symmetric and asymmetric internal loops in one stem"),
    "1RLG": ("RNP", "L7Ae on a box C/D kink-turn"),
    "1ZDI": ("RNP", "MS2 coat protein dimer on an RNA hairpin"),
    "1HYS": ("polymerase", "HIV-1 RT on a polypurine-tract RNA:DNA hybrid"),
    "1MSW": ("polymerase", "T7 RNA polymerase elongation complex: DNA duplex, RNA-DNA hybrid, single strands"),
    "3KTQ": ("polymerase", "Klentaq ternary complex: primer/template with a single-stranded template overhang"),
    "4YBB": ("ribosome", "E. coli ribosome: 16S and 23S rRNA, about 4,500 nt, cropped to what the proteins touch"),
    "1KX5": ("protein-DNA", "nucleosome: 147 bp wrapped 1.65 turns round the histone octamer"),
    "1IHF": ("protein-DNA", "IHF: DNA bent by more than 160° at two kinks"),
    "1YTB": ("protein-DNA", "TBP: TATA box bent by about 80° and unwound"),
    "1LMB": ("protein-DNA", "lambda repressor on a straight operator: the classic picture must not change"),
    "1ZAA": ("protein-DNA", "Zif268 zinc fingers in the major groove of a straight duplex"),
    "1RVA": ("protein-DNA", "EcoRV: restriction enzyme kinking its site by about 50°"),
    "3CRX": ("junction", "Cre-Holliday junction: four arms on four strands"),
    "4LD0": ("junction", "RuvC on a Holliday junction at 3.75 Å: four arms held by a resolvase dimer"),
    "1ZBI": ("hybrid", "RNase H on a 12-bp RNA:DNA hybrid"),
    "1KF1": ("non-canonical", "human telomeric G-quadruplex: no Watson-Crick pairs at all"),
    "4PLX": ("non-canonical", "MALAT1 triple helix: Hoogsteen third strand on a duplex"),
    "1QBJ": ("non-canonical", "Zalpha on left-handed Z-DNA"),
    "1BNA": ("DNA only", "Dickerson dodecamer: one straight B-DNA duplex, no protein"),
}


def _key(label) -> tuple[str, int, str]:
    return (label.chain, label.seq, label.icode)


# run by the throwaway venv's Python (python -I -c CODE FILE): the canonical and wobble pairs RNApolis finds in the
# first model, as JSON [[chain, number, icode], [chain, number, icode]] pairs
_REFERENCE = """
import gzip, json, sys, tempfile
from pathlib import Path
from rnapolis.annotator import find_pairs
from rnapolis.parser import read_3d_structure
from rnapolis.util import handle_input_file
src = Path(sys.argv[1])
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "model.cif"
    path.write_bytes(gzip.decompress(src.read_bytes()) if src.suffix == ".gz" else src.read_bytes())
    found = find_pairs(read_3d_structure(handle_input_file(str(path)), 1), 1)[0]
out = []
for bp in found:
    if bp.saenger is not None and bp.saenger.value in ("XIX", "XX", "XXVIII"):  # G-C, A-U/T, G-U wobble
        a, b = bp.nt1.auth, bp.nt2.auth
        out.append(sorted([[a.chain, a.number, (a.icode or "").strip()], [b.chain, b.number, (b.icode or "").strip()]]))
print(json.dumps(out))
"""
RNAPOLIS_PYTHON = os.environ.get("NA_RNAPOLIS_PYTHON", "")


def reference_pairs(path: Path) -> set[frozenset]:
    """The pairs RNApolis finds, from the Python of the venv that has it. Raises RuntimeError with the reason."""
    if not RNAPOLIS_PYTHON:
        raise RuntimeError("no --rnapolis-python")
    done = subprocess.run(
        [RNAPOLIS_PYTHON, "-I", "-c", _REFERENCE, str(path)], capture_output=True, text=True, timeout=600
    )
    if done.returncode:
        tail = (done.stderr.strip().splitlines() or ["failed"])[-1]
        raise RuntimeError(f"rnapolis could not read {path.name}: {tail}")
    return {frozenset(tuple(x) for x in pair) for pair in json.loads(done.stdout)}


def pairing(path) -> dict:
    """foldmap's pairs and stems against RNApolis' on the deposited model."""
    path = Path(path)
    try:
        theirs = reference_pairs(path)
    except (RuntimeError, subprocess.TimeoutExpired) as err:
        return {"pair_f1": None, "stem_recall": None, "note": str(err)}
    model, _ = read_model(path, "asu")
    found = _nucleic_strands(model)
    strands = [s for s, _ in found]
    ours = pair_nucleotides([run for _, run in found])
    mine = {frozenset((_key(strands[p.a].labels[p.i]), _key(strands[p.b].labels[p.j]))) for p in ours}
    tp = len(mine & theirs)
    prec = tp / len(mine) if mine else 1.0
    rec = tp / len(theirs) if theirs else 1.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    # stems: group the reference pairs exactly as foldmap groups its own, then match arms within one nucleotide
    index = {_key(l): (s, i) for s, strand in enumerate(strands) for i, l in enumerate(strand.labels)}
    ref_pairs = []
    for pair in sorted(theirs, key=sorted):
        if all(k in index for k in pair):
            x, y = sorted(index[k] for k in pair)
            ref_pairs.append(BasePair(x[0], x[1], y[0], y[1], "WC"))
    ref = secondary_structure(strands, ref_pairs).stems
    got = secondary_structure(strands, ours).stems

    def span(st):
        return st.arm5, st.arm3

    def close(a, b):
        return all(x[0] == y[0] and abs(x[1] - y[1]) <= 1 and abs(x[2] - y[2]) <= 1 for x, y in zip(span(a), span(b)))

    matched = sum(1 for r in ref if any(close(r, g) for g in got))
    return {
        "pairs_foldmap": len(mine),
        "pairs_reference": len(theirs),
        "pair_precision": round(prec, 3),
        "pair_recall": round(rec, 3),
        "pair_f1": round(f1, 3),
        "stems_foldmap": len(got),
        "stems_reference": len(ref),
        "stem_recall": round(matched / len(ref), 3) if ref else 1.0,
        "missed_pairs": sorted(sorted(p) for p in theirs - mine)[:20],
        "extra_pairs": sorted(sorted(p) for p in mine - theirs)[:20],
    }


def contact_sheet(name, lay, sses, bb, fig, metrics, out: Path) -> Path:
    """Figure | drawn layout with each stem in its own colour | 3D projected through the same page axes."""
    import io

    na = lay.nucleic
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    buf.seek(0)
    img = plt.imread(buf)
    cmap = plt.get_cmap("tab20")
    colour = {d.id: cmap(k % 20) for k, d in enumerate(lay.stems)}
    sheet = plt.figure(figsize=(18, 7.5))
    a0 = sheet.add_axes([0.0, 0.18, 0.36, 0.8])
    a0.imshow(img)
    a0.set_axis_off()
    a0.set_title(f"{name}: foldmap", loc="left")
    a1 = sheet.add_axes([0.37, 0.18, 0.3, 0.8])
    for p in lay.placed.values():
        x0, y0, x1, y1 = p.rect
        a1.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="#dddddd", ec="#999999", lw=0.5))
    for d in lay.stems:
        corners = [
            np.array([d.cx, d.cy]) + sx * d.direction * d.length / 2 + sy * d.normal * d.width / 2
            for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))
        ]
        a1.add_patch(plt.Polygon(corners, fc=colour[d.id], ec="k", lw=0.5, alpha=0.85))
        a1.annotate(d.name, (d.cx, d.cy), ha="center", va="center", fontsize=7)
        tip = np.array([d.cx, d.cy]) + d.direction * d.length / 2
        a1.annotate("", tip, (d.cx, d.cy), arrowprops={"arrowstyle": "->", "lw": 0.8})
    for p in lay.na_paths:
        q = np.asarray(p.points)
        a1.plot(q[:, 0], q[:, 1], lw=0.6, color="#555555")
    a1.set_aspect("equal")
    a1.autoscale()
    a1.set_title("drawn layout (stems coloured, arrows: first arm 5'→3')", loc="left", fontsize=9)
    a1.set_axis_off()
    a2 = sheet.add_axes([0.69, 0.18, 0.3, 0.8])
    frame = lay.frame  # the figure's own frame, so unrolled symmetric and radial views are projected as drawn

    def proj(xyz):
        return frame.project(np.atleast_2d(xyz)) * SCALE

    if len(bb):
        q = proj(bb.ca)
        chains = np.array([l.chain for l in bb.labels])
        cc = chain_colors(list(dict.fromkeys(chains)), "okabe-ito")
        for c in dict.fromkeys(chains):
            m = chains == c
            a2.plot(q[m, 0], q[m, 1], lw=0.5, color=cc[c], alpha=0.5)
    for strand in na.strands:
        q = proj(strand.p)
        a2.plot(q[:, 0], q[:, 1], lw=0.8, color="#333333")
    pieces = {p.id: p for p in stem_pieces(na)}
    for d in lay.stems:
        piece = pieces.get(d.id)
        if piece is None:
            continue
        pts = np.array([na.strands[k].c1[i] for k, i in piece.nucleotides()])
        q = proj(pts)
        a2.scatter(q[:, 0], q[:, 1], s=8, color=colour[d.id])
        c, ax = stem_axis(piece, na.strands)
        q0 = proj(c)[0]
        q1 = q0 + frame.direction(c, ax) * 8 * SCALE
        a2.annotate(d.name, q0, fontsize=7)
        a2.annotate("", q1, q0, arrowprops={"arrowstyle": "->", "lw": 0.8, "color": colour[d.id]})
    touched = {x for v in na.contacts.values() for x in v}
    if touched:
        q = proj(np.array([na.strands[k].c1[i] for k, i in touched]))
        a2.scatter(q[:, 0], q[:, 1], s=3, color="red", zorder=5)
    a2.set_aspect("equal")
    a2.set_title("3D, same page axes (red: touched nucleotides)", loc="left", fontsize=9)
    a2.set_axis_off()

    text = "  ".join(f"{k}={v}" for k, v in metrics.items() if k not in ("missed_pairs", "extra_pairs", "failed"))
    fails = ", ".join(metrics.get("failed", [])) or "none"
    sheet.text(0.01, 0.1, f"failed: {fails}", fontsize=10, color="red" if metrics.get("failed") else "green")
    sheet.text(
        0.01, 0.02, "\n".join(text[i : i + 230] for i in range(0, len(text), 230)), fontsize=7, family="monospace"
    )
    path = out / f"{name}.png"
    sheet.savefig(path, dpi=90)
    plt.close(sheet)
    return path


# per-entry make_layout options, for downloaded entries only (--files are cropped already): the ribosome is
# drawn as its small subunit (16S rRNA and its proteins)
OPTIONS = {"4YBB": {"assembly": "asu", "chains": [f"A{c}" for c in "BCDEFGHIJKLMNOPQRSTU"]}}
TIMEOUT = 900  # seconds an entry may take before it is reported as a failure


def run(name: str, path: Path, out: Path, options: dict | None = None) -> dict:
    m: dict = {"entry": name}
    m.update(pairing(path))
    t = time.time()
    options = options or {}
    model = read_model(path, options.get("assembly", "auto"))  # read once; the bead checks compare with it
    lay, sses, bb = make_layout(path, model=model, **options)
    loops = route_loops(lay, sses, bb)
    fig = draw(lay, loops, sses, name, model=model)
    m["seconds"] = round(time.time() - t, 1)
    if lay.nucleic is not None:
        m.update(shape(lay))
        m.update(contacts(lay, sses))
        m.update(render_checks(lay, loops, fig, sses, model=model, bb=bb))
    m["failed"] = judge(m)
    save(fig, out / f"{name}-figure.svg")
    if lay.nucleic is not None:
        contact_sheet(name, lay, sses, bb, fig, m, out)
    plt.close("all")
    return m


def _one(name, path, out, options, queue) -> None:
    try:
        queue.put(run(name, path, out, options))
    except Exception as err:  # one broken entry must not stop the campaign
        import traceback

        traceback.print_exc()
        queue.put({"entry": name, "failed": [f"crashed: {type(err).__name__}: {err}"]})


def _isolated(name: str, path: Path, out: Path, options: dict | None = None) -> dict:
    """run() in a child process, so an entry that hangs or exhausts memory is reported, not fatal."""
    import multiprocessing as mp

    queue = mp.get_context("fork").Queue()
    child = mp.get_context("fork").Process(target=_one, args=(name, path, out, options, queue))
    child.start()
    child.join(TIMEOUT)
    if child.is_alive():
        child.terminate()
        child.join()
        return {"entry": name, "failed": [f"timed out after {TIMEOUT} s"]}
    return queue.get() if not queue.empty() else {"entry": name, "failed": [f"died (exit code {child.exitcode})"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ids", help="comma-separated entries (default: the whole benchmark)")
    ap.add_argument("--files", nargs="*", type=Path, help="local files instead of downloaded entries")
    ap.add_argument("--out", type=Path, default=Path("na-report"))
    ap.add_argument("--rnapolis-python", help="Python of a venv that has rnapolis (default: $NA_RNAPOLIS_PYTHON)")
    args = ap.parse_args(argv)
    global RNAPOLIS_PYTHON
    RNAPOLIS_PYTHON = args.rnapolis_python or RNAPOLIS_PYTHON
    args.out.mkdir(parents=True, exist_ok=True)
    jobs = (
        [(f.name.split(".")[0], f) for f in args.files]
        if args.files
        else [(i, Path(fetch(i))) for i in (args.ids.split(",") if args.ids else BENCHMARK)]
    )
    rows = []
    for name, path in jobs:
        m = _isolated(name, path, args.out, None if args.files else OPTIONS.get(name))
        rows.append(m)
        print(f"{name}: {'PASS' if not m['failed'] else 'FAIL ' + '; '.join(m['failed'])}", flush=True)
    (args.out / "report.json").write_text(json.dumps(rows, indent=1, default=str), encoding="utf-8")
    cols = [
        "entry",
        "pairs_foldmap",
        "pairs_reference",
        "pair_precision",
        "pair_recall",
        "pair_f1",
        "stems_foldmap",
        "stems_reference",
        "stem_recall",
        "stems_drawn",
        "spearman",
        "procrustes",
        "angle_median",
        "angle_max",
        "mirrored",
        "bend_3d",
        "bend_drawn",
        "joint_error_max",
        "contact_auc",
        "contact_precision_at_k",
        "stems_nearer_frac",
        "beads_missing",
        "bead_colour_wrong",
        "beads_false",
        "beads_unmatched",
        "overlaps",
        "text_off_page",
        "letter_overlaps",
        "fallbacks",
        "na_crossings",
        "na_protein_crossings",
        "bead_overlaps",
        "gap_label_collisions",
        "seconds",
    ]
    lines = ["| " + " | ".join(cols + ["result"]) + " |", "|" + "---|" * (len(cols) + 1)]
    for m in rows:
        cells = ["" if m.get(c) is None else str(m[c]) for c in cols]
        lines.append(
            "| " + " | ".join(cells + ["pass" if not m["failed"] else "FAIL: " + "; ".join(m["failed"])]) + " |"
        )
    why = [f"- **{k}** ({c}): {w}" for k, (c, w) in BENCHMARK.items()]
    passed = sum(not m["failed"] for m in rows)
    notes = [f"- {m['entry']}: {m['note']}" for m in rows if m.get("note")]
    (args.out / "report.md").write_text(
        f"# foldmap nucleic-acid validation\n\n{passed} of {len(rows)} entries pass.\n\n"
        + "\n".join(lines)
        + "\n\nThresholds: "
        + ", ".join(f"{k} {v}" for k, v in THRESHOLDS.items())
        + ("\n\nNotes:\n" + "\n".join(notes) if notes else "")
        + "\n\n## Entries\n\n"
        + "\n".join(why)
        + "\n",
        encoding="utf-8",
    )
    print(f"{passed}/{len(rows)} pass; report in {args.out}/report.md")
    return 0 if passed == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
