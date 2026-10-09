"""Shared test helpers: run the step-1 pipeline, build synthetic SSEs."""

import numpy as np

from topoplot.dssp import assign_dssp
from topoplot.io import load_backbone
from topoplot.model import ResLabel, SSE
from topoplot.sheets import build_sheets
from topoplot.ss import build_sses


def pipeline(path, assembly="asu"):
    """Step-1 pipeline on the file as deposited (figures build the biological assembly; unit tests don't)."""
    bb = load_backbone(path, assembly)
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss)
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
