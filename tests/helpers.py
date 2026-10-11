"""Shared test helpers: run the step-1 pipeline, build synthetic SSEs."""

from pathlib import Path

import numpy as np

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import SSE, ResLabel
from foldmap.sheets import build_sheets
from foldmap.ss import build_sses

FID = Path(__file__).parent / "data" / "fidelity"  # plan-11 protein fixtures (.cif.gz, outside the bench)


def pipeline(path, assembly="asu"):
    """Step-1 pipeline on the file as deposited (figures build the biological assembly; unit tests don't)."""
    bb = load_backbone(path, assembly)
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss, bridges=d.bridges)
    return bb, sses, build_sheets(sses, d.bridges)


def fake_sse(kind, start, end, centroid=(0, 0, 0), axis=(1, 0, 0), chain="A"):
    lab = ResLabel(chain, start, "", "ALA")
    last = ResLabel(chain, end, "", "ALA")
    a = np.asarray(axis, float)
    return SSE(kind, chain, start, end, lab, last, np.asarray(centroid, float), a / np.linalg.norm(a))


def tripeptide(src, tmp_path):
    """First three residues of a structure, as a PDB file."""
    import gemmi

    full = tmp_path / "full.pdb"
    gemmi.read_structure(str(src)).write_pdb(str(full))
    keep = [l for l in full.read_text().splitlines() if l.startswith("ATOM") and int(l[22:26]) <= 3]
    out = tmp_path / "tri.pdb"
    out.write_text("\n".join(keep) + "\nEND\n")
    return out


def loop_crossings(loops) -> int:
    """Proper crossings between the polylines of different loops (shared end points do not count)."""

    def cross(p, q, r, s) -> bool:
        def side(a, b, c):
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

        return side(r, s, p) * side(r, s, q) < 0 and side(p, q, r) * side(p, q, s) < 0

    segs = [list(zip(lp.points, lp.points[1:])) for lp in loops]
    return sum(
        cross(p, q, r, s)
        for i in range(len(segs))
        for j in range(i + 1, len(segs))
        for p, q in segs[i]
        for r, s in segs[j]
    )
