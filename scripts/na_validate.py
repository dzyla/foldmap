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
Pairing is measured on the pairs the figure draws (assembly copies mapped back to their deposited chains by
stripping the copy number; a cropped figure against the reference pairs it draws both nucleotides of), and is
also reported on the model as deposited (deposited_pair_f1, deposited_stem_recall, not judged).
With --require-rnapolis an entry whose pairing cannot be measured fails; without it a warning is printed and the
report says how many entries were not validated. The exit code is 0 only if every entry passes. The run time
(`seconds`) is gated but is a timing, not a deterministic quantity.
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
    name = src.name[:-3] if src.suffix == ".gz" else src.name
    path = Path(tmp) / ("model.pdb" if name.lower().endswith((".pdb", ".ent")) else "model.cif")
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
REQUIRE_RNAPOLIS = False


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


def _same_shape(p: np.ndarray, q: np.ndarray, tol: float = 0.05) -> bool:
    """Whether two point sets (same order) differ only by a rigid motion: every pairwise distance within tol Å."""
    dp = np.sqrt(((p[:, None, :] - p[None, :, :]) ** 2).sum(-1))
    dq = np.sqrt(((q[:, None, :] - q[None, :, :]) ** 2).sum(-1))
    return bool(np.abs(dp - dq).max() <= tol) if len(p) else False


def _rigid(src: np.ndarray, dst: np.ndarray):
    """(R, t) of the proper rotation and translation taking src onto dst (Kabsch); None under 3 points."""
    if len(src) < 3:
        return None
    a, b = src.mean(axis=0), dst.mean(axis=0)
    u, _, vt = np.linalg.svd((src - a).T @ (dst - b))
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, b - r @ a


def _chain_map(drawn: dict, residues: dict, c1: dict, as_deposited: bool):
    """(drawn chain -> (deposited chain, copy), [drawn chains left unmapped]). drawn: chain -> [(ResLabel, C1'
    xyz)]; residues: deposited chain -> {(number, icode, name)}; c1: deposited (chain, number, icode) -> C1' xyz.
    A model as deposited keeps its names (copy ""). An assembly names each copy by its deposited chain plus the
    operator's number (gemmi AddNumber): every split of the name into a deposited name and digits is a candidate,
    kept only if that deposited chain holds every drawn residue (number, icode and name) in the same shape (a
    rigid copy). Exactly one candidate maps the chain; none, or more than one (deposited names ending in digits,
    identical twins), leaves it unmapped: the entry fails rather than guessing."""
    out, unmapped = {}, []
    for c, items in sorted(drawn.items()):
        if as_deposited:
            cands = [(c, "")] if c in residues else []
        else:
            cands = [(d, c[len(d) :]) for d in sorted(residues) if c[len(d) :].isdigit() and c.startswith(d)]
        fits = []
        for d, n in cands:
            keys = [(d, l.seq, l.icode) for l, _ in items]
            if all((l.seq, l.icode, l.name) in residues[d] for l, _ in items) and all(k in c1 for k in keys):
                if _same_shape(np.array([x for _, x in items]), np.array([c1[k] for k in keys])):
                    fits.append((d, n))
        if len(fits) == 1:
            out[c] = fits[0]
        else:
            unmapped.append(c)
    return out, unmapped


def _score(strands, ours, theirs: set[frozenset], where: dict, drawn_only: bool, same_operator=None) -> dict:
    """Pairs and stems of `ours` (BasePairs on `strands`) against the reference pairs `theirs` (frozensets of
    deposited (chain, number, icode) keys). where: (strand, nucleotide) -> (copy, deposited key), one nucleotide
    per (copy, key). A reference pair is expected wherever both its nucleotides are drawn in one copy, or in two
    copies that same_operator(x, y) says one operator made (gemmi numbers each generator's operators apart, so
    one operator's copies may carry two numbers); a drawn pair counts only as such an expected pair, so a pair
    drawn between a nucleotide and its partner's twin in another copy is wrong. drawn_only: recall counts only
    the expected pairs (a cropped figure is not short of the pairs it leaves out); otherwise also every
    reference pair not drawn at all."""
    by_key: dict = {}
    for x, (c, k) in sorted(where.items()):
        by_key.setdefault(k, []).append((c, x))
    expected: dict[frozenset, frozenset] = {}  # drawn nucleotide pair -> reference pair
    for pair in sorted(theirs, key=sorted):
        k1, k2 = sorted(pair)
        for c1_, x in by_key.get(k1, []):
            for c2_, y in by_key.get(k2, []):
                if c1_ == c2_ or (same_operator is not None and same_operator(x, y)):
                    expected[frozenset((x, y))] = pair
    mine = {frozenset(((p.a, p.i), (p.b, p.j))) for p in ours}
    tp = len(mine & expected.keys())
    shown = set(expected.values())
    absent = 0 if drawn_only else len(theirs - shown)
    prec = tp / len(mine) if mine else 1.0
    rec = tp / (len(expected) + absent) if expected or absent else 1.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    # stems: group the expected reference pairs exactly as foldmap groups its own, then match arms within one
    # nucleotide
    ref_pairs = []
    for xy in sorted(expected, key=sorted):
        x, y = sorted(xy)
        ref_pairs.append(BasePair(x[0], x[1], y[0], y[1], "WC"))
    ref = secondary_structure(strands, ref_pairs).stems
    got = secondary_structure(strands, ours).stems

    def span(st):
        return st.arm5, st.arm3

    def close(a, b):
        return all(x[0] == y[0] and abs(x[1] - y[1]) <= 1 and abs(x[2] - y[2]) <= 1 for x, y in zip(span(a), span(b)))

    def named(xy):  # a drawn pair as [[chain, number, icode], ...], chains as drawn
        return sorted([strands[s].chain, strands[s].labels[i].seq, strands[s].labels[i].icode] for s, i in xy)

    matched = sum(1 for r in ref if any(close(r, g) for g in got))
    missed = sorted(named(xy) for xy in expected.keys() - mine)
    if not drawn_only:  # reference pairs with a nucleotide foldmap has no strand for
        missed += sorted(sorted(p) for p in theirs - shown)
    return {
        "pairs_foldmap": len(mine),
        "pairs_reference": len(expected) + absent,
        "pair_precision": round(prec, 3),
        "pair_recall": round(rec, 3),
        "pair_f1": round(f1, 3),
        "stems_foldmap": len(got),
        "stems_reference": len(ref),
        "stem_recall": round(matched / len(ref), 3) if ref else 1.0,
        "missed_pairs": missed[:20],
        "extra_pairs": sorted(named(xy) for xy in mine - expected.keys())[:20],
    }


def pairing(path, theirs: set[frozenset]) -> dict:
    """foldmap's pairs and stems against RNApolis' on the deposited model (the asymmetric unit), every key
    prefixed `deposited_`; reported, not judged (the figure may be drawn from an assembly: drawn_pairing)."""
    model, _ = read_model(Path(path), "asu")
    found = _nucleic_strands(model)
    strands = [s for s, _ in found]
    ours = pair_nucleotides([run for _, run in found])
    where = {(s, i): ("", _key(l)) for s, strand in enumerate(strands) for i, l in enumerate(strand.labels)}
    m = _score(strands, ours, theirs, where, drawn_only=False)
    return {f"deposited_{k}": v for k, v in m.items() if k not in ("missed_pairs", "extra_pairs")}


def drawn_pairing(na, theirs: set[frozenset], asu_model, as_deposited: bool) -> dict:
    """The pairs and stems the figure draws (`na`, the layout's nucleic acid: assembly copies, cropped) against
    RNApolis' on the deposited model `asu_model`, each drawn chain mapped back to the deposited one it copies
    (_chain_map; as_deposited: the figure is drawn from the deposited model itself). A true pair joins two
    nucleotides of one copy, or of two copies one operator made (the same rigid motion of their deposited
    chains); a pair drawn between a nucleotide and its partner's twin in another copy is wrong. Drawn chains
    that cannot be mapped (unmapped_chains), or a copy holding one deposited nucleotide twice, fail the entry."""
    residues = {c.name: {(r.seqid.num, r.seqid.icode.strip(), r.name) for r in c} for c in asu_model}
    c1 = {_key(l): s.c1[i] for s, _ in _nucleic_strands(asu_model) for i, l in enumerate(s.labels)}
    drawn: dict = {}
    for strand in na.strands:
        drawn.setdefault(strand.chain, []).extend(zip(strand.labels, strand.c1))
    chains, unmapped = _chain_map(drawn, residues, c1, as_deposited)
    where, seen, twice = {}, {}, set()
    for s, strand in enumerate(na.strands):
        if strand.chain not in chains:
            continue
        d, n = chains[strand.chain]
        for i, l in enumerate(strand.labels):
            at = (n, (d, l.seq, l.icode))
            if at in seen:
                twice |= {strand.chain, na.strands[seen[at][0]].chain}
            seen[at] = (s, i)
            where[(s, i)] = at
    if twice:  # one copy holding a deposited nucleotide twice: no single answer
        where = {x: v for x, v in where.items() if na.strands[x[0]].chain not in twice}
        unmapped = sorted(set(unmapped) | twice)
    motion = {}  # drawn chain -> rigid motion of its deposited chain onto it
    for c, items in drawn.items():
        if c in chains and c not in twice:
            d = chains[c][0]
            motion[c] = _rigid(np.array([c1[(d, l.seq, l.icode)] for l, _ in items]), np.array([x for _, x in items]))

    def same_operator(x, y) -> bool:
        a, b = motion.get(na.strands[x[0]].chain), motion.get(na.strands[y[0]].chain)
        return a is not None and b is not None and np.abs(a[0] - b[0]).max() < 1e-3 and _norm3(a[1] - b[1]) < 0.05

    ours = [p for p in na.base_pairs if (p.a, p.i) in where and (p.b, p.j) in where]
    m = _score(na.strands, ours, theirs, where, drawn_only=True, same_operator=same_operator)
    m["unmapped_chains"] = unmapped
    return m


def _norm3(v) -> float:
    return float(np.sqrt((np.asarray(v, float) ** 2).sum()))


def _bounds(points, margin: float = 0.06):
    """(x0, x1, y0, y1) of the points, widened by a fraction of the larger side."""
    q = np.vstack([np.atleast_2d(p) for p in points])
    (x0, y0), (x1, y1) = q.min(axis=0), q.max(axis=0)
    pad = margin * max(x1 - x0, y1 - y0, 1e-6)
    return x0 - pad, x1 + pad, y0 - pad, y1 + pad


def wrap_metrics(items: list[str], width: int) -> list[str]:
    """Lines of at most `width` characters (more only for a single longer item), broken between items."""
    lines, line = [], ""
    for item in items:
        if line and len(line) + 2 + len(item) > width:
            lines.append(line)
            line = ""
        line = f"{line}  {item}" if line else item
    return lines + ([line] if line else [])


def contact_sheet(name, lay, sses, bb, fig, metrics, out: Path) -> Path:
    """Figure | drawn layout with each stem in its own colour | 3D seen through the figure's own frame, cropped to
    the nucleic acid. The sheet's shape follows the panels' content."""
    import io

    na = lay.nucleic
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    buf.seek(0)
    img = plt.imread(buf)
    cmap = plt.get_cmap("tab20")
    colour = {d.id: cmap(k % 20) for k, d in enumerate(lay.stems)}
    frame = lay.frame  # unrolled symmetric and radial views are projected as drawn

    def proj(xyz):
        return frame.project(np.atleast_2d(xyz)) * SCALE

    corners = {
        d.id: [
            np.array([d.cx, d.cy]) + sx * d.direction * d.length / 2 + sy * d.normal * d.width / 2
            for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))
        ]
        for d in lay.stems
    }
    drawn_pts = [np.array([[p.rect[0], p.rect[1]], [p.rect[2], p.rect[3]]]) for p in lay.placed.values()]
    drawn_pts += [np.array(c) for c in corners.values()] + [np.asarray(p.points)[:, :2] for p in lay.na_paths]
    drawn_box = _bounds(drawn_pts) if drawn_pts else (0, 1, 0, 1)
    seen_pts = [proj(s.p) for s in na.strands] + [proj(s.c1) for s in na.strands]
    seen_box = _bounds(seen_pts)
    boxes = [(0, img.shape[1], 0, img.shape[0]), drawn_box, seen_box]
    aspects = [min(max((b[1] - b[0]) / max(b[3] - b[2], 1e-9), 0.4), 2.5) for b in boxes]
    panel_h, gap, text_h = 6.0, 0.25, 1.5
    widths = [panel_h * a for a in aspects]
    total_w = sum(widths) + gap * (len(widths) + 1)
    sheet = plt.figure(figsize=(total_w, panel_h + text_h + 0.5))
    axes, x = [], gap
    for w in widths:
        axes.append(
            sheet.add_axes(
                [x / total_w, text_h / (panel_h + text_h + 0.5), w / total_w, panel_h / (panel_h + text_h + 0.5)]
            )
        )
        x += w + gap
    a0, a1, a2 = axes
    a0.imshow(img)
    a0.set_title(f"{name}: foldmap", loc="left", fontsize=10)
    for p in lay.placed.values():
        x0, y0, x1, y1 = p.rect
        a1.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="#dddddd", ec="#999999", lw=0.5))
    for d in lay.stems:
        a1.add_patch(plt.Polygon(corners[d.id], fc=colour[d.id], ec="k", lw=0.5, alpha=0.85))
        a1.annotate(d.name, (d.cx, d.cy), ha="center", va="center", fontsize=7)
        tip = np.array([d.cx, d.cy]) + d.direction * d.length / 2
        a1.annotate("", tip, (d.cx, d.cy), arrowprops={"arrowstyle": "->", "lw": 0.8})
    for p in lay.na_paths:
        q = np.asarray(p.points)
        a1.plot(q[:, 0], q[:, 1], lw=0.6, color="#555555")
    a1.set_xlim(*drawn_box[:2])
    a1.set_ylim(*drawn_box[2:])
    a1.set_title("drawn layout (stems coloured, arrows: first arm 5'->3')", loc="left", fontsize=9)
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
    touched = sorted({x for v in na.contacts.values() for x in v})
    if touched:
        q = proj(np.array([na.strands[k].c1[i] for k, i in touched]))
        a2.scatter(q[:, 0], q[:, 1], s=3, color="red", zorder=5)
    a2.set_xlim(*seen_box[:2])
    a2.set_ylim(*seen_box[2:])
    a2.set_title("3D, figure's frame (red: touched nucleotides)", loc="left", fontsize=9)
    for a in (a1, a2):
        a.set_aspect("equal", adjustable="box")
    for a in axes:
        a.set_axis_off()
    items = [f"{k}={v}" for k, v in metrics.items() if k not in ("missed_pairs", "extra_pairs", "failed", "note")]
    fails = "; ".join(metrics.get("failed", [])) or "none"
    size = 9
    chars = int((total_w - 0.4) / (size * 0.62 / 72))  # monospace glyph width in inches
    lines = wrap_metrics(items, chars)
    sheet.text(
        0.15 / total_w,
        1.25 / (panel_h + text_h + 0.5),
        f"failed: {fails}"[: chars * 2],
        fontsize=10,
        color="red" if metrics.get("failed") else "green",
    )
    sheet.text(
        0.15 / total_w, 0.1 / (panel_h + text_h + 0.5), "\n".join(lines), fontsize=size, family="monospace", va="bottom"
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
    m: dict = {"entry": name, "pair_f1": None, "stem_recall": None}
    try:
        theirs = reference_pairs(path)
    except (RuntimeError, subprocess.TimeoutExpired) as err:
        theirs, m["note"] = None, str(err)
    else:
        m.update(pairing(path, theirs))
    t = time.time()
    options = options or {}
    model = read_model(path, options.get("assembly", "auto"))  # read once; the bead checks compare with it
    lay, sses, bb = make_layout(path, model=model, **options)
    loops = route_loops(lay, sses, bb)
    fig = draw(lay, loops, sses, name, model=model)
    m["seconds"] = round(time.time() - t, 1)
    if lay.nucleic is not None and theirs is not None:
        asu_model, deposited = read_model(path, "asu")
        as_deposited = model[1] == deposited and {c.name for c in model[0]} == deposited
        m.update(drawn_pairing(lay.nucleic, theirs, asu_model, as_deposited))
    if lay.nucleic is not None:
        m.update(shape(lay))
        m.update(contacts(lay, sses))
        m.update(render_checks(lay, loops, fig, sses, model=model, bb=bb))
    m["failed"] = judge(m)
    if lay.nucleic is None:  # nothing to validate is not a pass
        m["failed"].append("no nucleic acid drawn")
    if m.get("unmapped_chains"):  # drawn pairs that cannot be scored are not a pass
        m["failed"].append(
            f"pairing: drawn chains not mapped to the deposited model: {', '.join(m['unmapped_chains'])}"
        )
    if REQUIRE_RNAPOLIS and m.get("pair_f1") is None:
        m["failed"].append(f"pairing not measured ({m.get('note', 'no reason given')})")
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
    try:  # read before join: a child blocked on a full pipe cannot exit
        result = queue.get(timeout=TIMEOUT)
    except Exception:
        result = None
    if result is None:
        child.terminate()
    child.join()
    if result is not None:
        return result
    if child.exitcode in (None, -15):
        return {"entry": name, "failed": [f"timed out after {TIMEOUT} s"]}
    return {"entry": name, "failed": [f"died (exit code {child.exitcode})"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ids", help="comma-separated entries (default: the whole benchmark)")
    ap.add_argument("--files", nargs="*", type=Path, help="local files instead of downloaded entries")
    ap.add_argument("--out", type=Path, default=Path("na-report"))
    ap.add_argument("--rnapolis-python", help="Python of a venv that has rnapolis (default: $NA_RNAPOLIS_PYTHON)")
    ap.add_argument("--require-rnapolis", action="store_true", help="fail entries whose pairing is not measured")
    args = ap.parse_args(argv)
    global RNAPOLIS_PYTHON, REQUIRE_RNAPOLIS
    RNAPOLIS_PYTHON = args.rnapolis_python or RNAPOLIS_PYTHON
    REQUIRE_RNAPOLIS = args.require_rnapolis
    args.out.mkdir(parents=True, exist_ok=True)
    if args.files is not None:
        jobs = [(f.name.split(".")[0], f) for f in args.files]
    else:
        jobs = [(i, None) for i in (args.ids.split(",") if args.ids else BENCHMARK)]
    if not RNAPOLIS_PYTHON:
        print("warning: no --rnapolis-python; pairing is not validated", file=sys.stderr)
    rows = []
    for name, path in jobs:
        if path is None:
            try:
                path = Path(fetch(name))
            except Exception as err:  # a failed download fails that entry, not the campaign
                m = {"entry": name, "failed": [f"download failed: {type(err).__name__}: {err}"]}
            else:
                m = _isolated(name, path, args.out, OPTIONS.get(name))
        else:
            m = _isolated(name, path, args.out, None)
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
        "deposited_pair_f1",
        "deposited_stem_recall",
        "unmapped_chains",
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
    unvalidated = sum(1 for m in rows if m.get("pair_f1") is None)
    if unvalidated:
        print(f"warning: pairing not validated for {unvalidated} entries", file=sys.stderr)
    (args.out / "report.md").write_text(
        f"# foldmap nucleic-acid validation\n\n{passed} of {len(rows)} entries pass.\n\n"
        + (f"pairing not validated for {unvalidated} entries\n\n" if unvalidated else "")
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
