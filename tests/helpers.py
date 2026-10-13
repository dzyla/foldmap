"""Shared test helpers: run the step-1 pipeline, build synthetic SSEs."""

from pathlib import Path

import numpy as np

from foldmap.dssp import assign_dssp
from foldmap.figcheck import honest_crossings, loop_crossings  # noqa: F401 - moved to figcheck (plan 10)
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


def text_centre(t) -> np.ndarray:
    """The page point a stem name or base letter is centred on: render starts it at layout.text_anchor (ha "left",
    va "baseline", rotation_mode "anchor", plan 10 Task 4), so its position is where its pen starts; other text is
    drawn centred on its position."""
    from foldmap.layout import text_anchor

    xy = np.asarray(t.get_position(), float)
    if t.get_horizontalalignment() != "left" or t.get_rotation_mode() != "anchor":
        return xy
    ax = t.axes
    px_per_unit = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    size = t.get_fontsize() * ax.figure.dpi / 72 / px_per_unit  # page units to the em
    bold = t.get_fontweight() in ("bold", 700)
    return xy - np.asarray(text_anchor(t.get_text(), bold, (0.0, 0.0), size, t.get_rotation()))


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
