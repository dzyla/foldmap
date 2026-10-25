"""Shared test helpers: run the step-1 pipeline, build synthetic SSEs."""

import os
from pathlib import Path

import numpy as np
import pytest

from foldmap.dssp import assign_dssp
from foldmap.figcheck import honest_crossings, loop_crossings  # noqa: F401 - moved to figcheck (plan 10)
from foldmap.io import load_backbone
from foldmap.model import SSE, ResLabel
from foldmap.sheets import build_sheets
from foldmap.ss import build_sses

FID = Path(__file__).parent / "data" / "fidelity"  # plan-11 protein fixtures (.cif.gz, outside the bench)
# the random-100 campaign's downloads (scripts/random100.py run, through foldmap.fetch): whole entries too large for a
# fixture, used where they are; a test on one is skipped where they are not
R100_CACHE = Path(os.environ.get("FOLDMAP_R100_CACHE") or Path.home() / ".cache" / "foldmap")


def r100_entry(entry: str) -> Path:
    """The random-100 campaign's download of `entry` (R100_CACHE), or skip the test."""
    path = R100_CACHE / f"{entry}.cif"
    if not path.is_file():
        pytest.skip(f"{entry} not downloaded to {R100_CACHE} (python scripts/random100.py run --ids {entry})")
    return path


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


def sideways_in_slab(loops, bands, placed, margin: float = 0.3) -> list[str]:
    """The router's own measure of a loop run along the membrane slab (A18's loop test until plan 27 made A18 strict,
    figcheck.loops_in_band): the loops (a>b) that run sideways (more than 0.5) inside a band, more than `margin` from
    its faces, but beside an in-band port (its first two segments, and its first run as far as it stays over its own
    element's body grown by the router's clearance). The router's tests (C2b-5, 7JI3) hold it to this; A18 holds
    every loop to more."""
    from foldmap.route import _INFLATE

    def inside(y, m):
        return any(lo + m < y < hi - m for lo, hi in bands)

    def within(port, step, element) -> bool:  # does a run of `step` from the port stay over its element's body?
        if element is None:
            return False
        x0, _, x1, _ = element.rect
        reach = (x1 + _INFLATE - port[0]) if step > 0 else (port[0] - (x0 - _INFLATE))
        return abs(step) <= max(reach, 0.0) + 1e-9

    out = []
    for lp in loops:
        pts, n = lp.points, len(lp.points)
        segs = list(zip(pts, pts[1:]))
        runs = [k for k, ((xa, ya), (xb, yb)) in enumerate(segs) if abs(ya - yb) < 1e-9 and abs(xa - xb) > 0.5]
        start_in, end_in = inside(pts[0][1], 0.0), inside(pts[-1][1], 0.0)
        for k in runs:
            (xa, ya), (xb, _) = segs[k]
            if not inside(ya, margin):
                continue
            if (start_in and k <= 1) or (end_in and k >= n - 3):  # beside its own port (C26)
                continue
            if start_in and k == runs[0] and within(pts[0], xb - xa, placed.get(lp.a_id)):
                continue
            if end_in and k == runs[-1] and within(pts[-1], xa - xb, placed.get(lp.b_id)):
                continue
            out.append(f"{lp.a_id}>{lp.b_id}")
            break
    return out
