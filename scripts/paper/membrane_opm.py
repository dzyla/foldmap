# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""M2 membrane vs OPM, step 2 (zyla-lab): foldmap's membrane estimate against OPM's published orientation.

    python scripts/paper/membrane_opm.py   (M2_RERUN=1 scripts/paper/run_all.sh)

Entries: ../data/opm/selection.json (opm_fetch.py; the rule is stated there and in scripts/paper/README.md). For each:

  estimate  foldmap.membrane.find_membrane on the PDB entry's mmCIF (the file `foldmap plot --membrane auto` reads,
            from ~/.cache/foldmap, no DUM atoms in it) with the backbone the membrane view draws (assembly auto),
            with the code at the harness commit (build/tree/src; asserted)
  OPM       OPM's oriented file (../data/opm/<id>.pdb.gz): planes from its DUM atoms ('N' inner, 'O' outer: checked
            per entry against the record's topology where it can be), and the record's transmembrane segments
            (the primary structure's record for a secondary representation; a segment counts only when at least
            three of its CA atoms are in the file)
  frame     OPM's coordinates are the deposited ones turned and moved: a Kabsch fit of OPM's CA atoms onto the drawn
            backbone's (each OPM chain paired with the drawn chain whose residues it fits best, pairs within 3 A
            kept, refitted) takes OPM's normal and planes into the frame foldmap estimated in

Measures per entry: thickness difference (2 x half against OPM's plane distance), angle between the normals (as
axes, 0-90 deg), offset of the mid-planes at the centroid of OPM's transmembrane CA atoms, sides (when foldmap calls
them: does its outward normal point to OPM's outer side), and for each OPM segment whether its CA run crosses the
estimated slab by the rule foldmap uses for its own elements (reaching membrane._MIN_CROSSING = 10 A beyond the
mid-plane on both sides), with the same rule against OPM's own planes as the ceiling.

Rule attribution (plan 15 Task 6, a963bb6): for every entry, and for 7KH0 (the Task 6 sweep's other refusal that
turned into a find; not in the selection), the estimate is recomputed with each of the two a963bb6 rules on or off:
  fallback    a best slab no long element spans moves to the best spanned membrane-like slab
  chain rule  the oily/polar fractions count only the chains with a long element across the slab
(both on reproduces membrane._estimate exactly; asserted; both off is the code before a963bb6).

Writes ../data/m2/membrane/results.json (everything, per entry and per segment).
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import MANUSCRIPT, TREE, run_stamp  # noqa: E402

OPM = MANUSCRIPT / "data" / "opm"
OUT = MANUSCRIPT / "data" / "m2" / "membrane"
EXTRA_ATTRIBUTION = ["7KH0"]  # Task 6 sweep: refused at a963bb6^, found at a963bb6 (NOTES.md 2026-10-18)
MIN_SEG_CA = 3
PAIR_CUT = 3.0  # A: a fitted CA pair this far apart is not the same atom (another copy)


def kabsch(src, dst):
    import numpy as np

    a, b = src.mean(axis=0), dst.mean(axis=0)
    u, _, vt = np.linalg.svd((src - a).T @ (dst - b))
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = b - r @ a
    rmsd = float(np.sqrt((((src @ r.T + t) - dst) ** 2).sum(1).mean()))
    return r, t, rmsd


def read_opm(path: Path):
    """CA atoms {(chain, num, icode): (resname, xyz)} and the DUM planes of OPM's oriented file."""
    import gemmi
    import numpy as np

    with tempfile.TemporaryDirectory() as tmp:
        plain = Path(tmp) / "opm.pdb"
        plain.write_bytes(gzip.decompress(path.read_bytes()))
        st = gemmi.read_structure(str(plain))
    # DUM atoms are read from the text: in some OPM files they follow the END record, where gemmi stops reading
    # (7SMR, 6PS3, ...); foldmap's membrane._opm reads them with gemmi and would miss those (scripts/paper/README.md)
    text = gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
    dum = [
        (line[12:16].strip(), float(line[46:54]))
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[17:20] == "DUM"
    ]
    ca = {}
    for ch in st[0]:
        for r in ch:
            if r.name == "DUM":
                continue
            a = r.find_atom("CA", "*")
            if (
                a is not None
                and gemmi.find_tabulated_residue(r.name)
                and gemmi.find_tabulated_residue(r.name).is_amino_acid()
            ):
                ca[(ch.name, r.seqid.num, r.seqid.icode.strip())] = (r.name, np.array([a.pos.x, a.pos.y, a.pos.z]))
    zn = [z for n, z in dum if n == "N"]
    zo = [z for n, z in dum if n == "O"]
    return ca, {
        "inner_z": float(np.mean(zn)) if zn else None,
        "outer_z": float(np.mean(zo)) if zo else None,
        "n_dum": len(dum),
    }


def fit_frame(ca_opm: dict, bb):
    """(R, t, rmsd, n_pairs, chain_map) with OPM ~ R @ drawn + t."""
    import numpy as np

    drawn = {}
    for k, l in enumerate(bb.labels):
        drawn.setdefault(l.chain, {})[(l.seq, l.icode)] = (l.name, bb.ca[k])
    by_chain = {}
    for (c, n, i), v in ca_opm.items():
        by_chain.setdefault(c, {})[(n, i)] = v

    def pairs(x, y):
        common = [k for k in by_chain[x] if k in drawn[y] and drawn[y][k][0] == by_chain[x][k][0]]
        return common

    # per OPM chain, the drawn chain it fits best on its own
    best = {}
    for x in by_chain:
        for y in drawn:
            common = pairs(x, y)
            if len(common) < 10:
                continue
            src = np.array([drawn[y][k][1] for k in common])
            dst = np.array([by_chain[x][k][1] for k in common])
            _, _, rmsd = kabsch(src, dst)
            cand = (rmsd > 1.5, -len(common), rmsd, y != x, y)
            if x not in best or cand < best[x][0]:
                best[x] = (cand, y)
    if not best:
        return None
    # the transform from the OPM chain with most well-fitting pairs; then each OPM chain's drawn copy under it
    anchor = min(best, key=lambda x: best[x][0])
    y0 = best[anchor][1]
    common = pairs(anchor, y0)
    r, t, _ = kabsch(np.array([drawn[y0][k][1] for k in common]), np.array([by_chain[anchor][k][1] for k in common]))
    src, dst, cmap = [], [], {}
    for x in by_chain:
        pick = None
        for y in drawn:
            common = pairs(x, y)
            if len(common) < 10:
                continue
            s = np.array([drawn[y][k][1] for k in common])
            d = np.array([by_chain[x][k][1] for k in common])
            dev = np.sqrt((((s @ r.T + t) - d) ** 2).sum(1))
            ok = dev < PAIR_CUT
            if pick is None or ok.sum() > pick[0]:
                pick = (int(ok.sum()), y, s[ok], d[ok])
        if pick and pick[0] >= 10:
            cmap[x] = pick[1]
            src.append(pick[2])
            dst.append(pick[3])
    src, dst = np.concatenate(src), np.concatenate(dst)
    r, t, rmsd = kabsch(src, dst)
    return r, t, rmsd, len(src), cmap


def segments(rec: dict) -> list[tuple[str, int, int]]:
    out = []
    for sub in rec.get("subunits") or []:
        for a, b in re.findall(r"\(\s*(-?\d+)\s*-\s*(-?\d+)\s*\)", sub.get("segment") or ""):
            out.append((sub["protein_letter"], int(a), int(b)))
    return out


def variant(bb, sses, fallback: bool, chain_rule: bool):
    """membrane._estimate with each a963bb6 rule on or off: (centre, normal) or None."""
    import numpy as np

    from foldmap import membrane as m

    long = [s for s in sses if (s.kind == "H" and len(s) >= 16) or (s.kind == "E" and len(s) >= 8)]
    if not long:
        return None
    normal = np.linalg.eigh(sum(len(s) * np.outer(s.axis, s.axis) for s in long))[1][:, -1]
    z = bb.ca @ normal
    hyd = np.array([l.name in m._HYDROPHOBIC for l in bb.labels])
    pol = np.array([l.name in m._POLAR for l in bb.labels])
    chain = np.array([l.chain for l in bb.labels]) if chain_rule else np.array(["*"] * len(bb.labels))
    slabs = []
    for c in np.arange(z.min(), z.max() + 1.0, 1.0):
        inside = np.abs(z - c) < m.HALF
        if inside.sum() < 20:
            continue
        slabs.append((hyd[inside].sum() - 1.5 * pol[inside].sum(), c, inside))
    if not slabs:
        return None
    best = max(slabs, key=lambda t: t[0])

    def like(slab):
        if chain_rule:
            return m._membrane_like(long, z, slab, hyd, pol, chain)
        # before a963bb6: every residue in the slab
        _, centre, inside = slab
        crossing = [
            s
            for s in long
            if z[s.start : s.end + 1].min() < centre - m._MIN_CROSSING
            and z[s.start : s.end + 1].max() > centre + m._MIN_CROSSING
        ]
        polar = m._MAX_POLAR_CROSSED if len(crossing) >= m._MANY_CROSSING else m._MAX_POLAR
        return bool(crossing) and hyd[inside].mean() >= m._MIN_HYDROPHOBIC and pol[inside].mean() <= polar

    if fallback and not m._spans(long, z, best[1]):
        spanned = [t for t in slabs if m._spans(long, z, t[1]) and like(t)]
        best = max(spanned, key=lambda t: t[0]) if spanned else best
    if not like(best):
        return None
    return float(best[1]), normal


def attribution(path, bb) -> dict:
    """The estimate with each a963bb6 rule on or off (non-barrel path only; a barrel never reaches _estimate)."""
    from foldmap import membrane as m
    from foldmap.dssp import assign_dssp
    from foldmap.ss import build_sses

    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, bridges=dssp.bridges)
    if m._barrel_estimate(bb, sses, dssp.bridges) is not None:
        return {"barrel": True}
    out = {"barrel": False}
    for name, fb, cr in (
        ("neither", False, False),
        ("fallback_only", True, False),
        ("chain_rule_only", False, True),
        ("both", True, True),
    ):
        v = variant(bb, sses, fb, cr)
        out[name] = None if v is None else {"centre": round(v[0], 4)}
        if name == "both":
            both = v
    real = m._estimate(bb, sses)
    assert (real is None) == (both is None) and (real is None or abs(real.centre - both[0]) < 1e-9), (real, both)
    return out


def depth_range(xyz, normal, centre):
    import numpy as np

    d = np.asarray(xyz) @ normal - centre
    return float(d.min()), float(d.max())


def one(sel: dict) -> dict:
    import numpy as np

    from foldmap.fetch import fetch
    from foldmap.io import load_backbone
    from foldmap.membrane import _MIN_CROSSING, find_membrane

    pid = sel["id"]
    t0 = time.time()
    path = Path(fetch(pid))
    bb = load_backbone(path, "auto")
    mem = find_membrane(path, bb)
    rec = json.loads((OPM / f"{pid}.json").read_text(encoding="utf-8"))
    ca_opm, planes = read_opm(OPM / f"{pid}.pdb.gz")
    e = {
        "id": pid,
        "sample": sel["sample"],
        "classtype": sel["classtype"],
        "opm_primary": sel["opm_primary"],
        "secondary": sel["secondary"],
        "file": path.name,
        "residues_drawn": len(bb),
        "chains_drawn": sorted({l.chain for l in bb.labels}),
        "opm_record_thickness": rec["thickness"],
        "opm_record_tilt": rec["tilt"],
        "opm_membrane": rec["membrane"]["name"],
        "opm_topology_in": rec["membrane"].get("topology_in"),
        "opm_topology_out": rec["membrane"].get("topology_out"),
    }
    zi, zo = planes["inner_z"], planes["outer_z"]
    e["opm_planes"] = planes
    if zi is None or zo is None:
        e["error"] = "OPM file has no inner/outer DUM planes"
        return e
    s = 1.0 if zo > zi else -1.0  # OPM's outward direction along its z
    zc = (zi + zo) / 2
    e["opm_thickness"] = round(abs(zo - zi), 3)
    fit = fit_frame(ca_opm, bb)
    if fit is None:
        e["error"] = "no OPM chain fits a drawn chain"
        return e
    r, t, rmsd, npairs, cmap = fit
    e["fit"] = {"rmsd": round(rmsd, 3), "pairs": npairs, "opm_ca": len(ca_opm), "chain_map": cmap}
    n_opm = s * (r.T @ np.array([0.0, 0.0, 1.0]))
    c_opm = s * (zc - t[2])
    # OPM's own topology check: the N-terminal CA of the topology subunit is on the 'in' side when topology_show_in
    sub = rec.get("topology_subunit")
    first = sorted(k for k in ca_opm if k[0] == sub)
    if first:
        zfirst = ca_opm[first[0]][1][2]
        side_in = s * (zfirst - zc) < 0
        e["opm_topology_check"] = {
            "subunit": sub,
            "show_in": rec.get("topology_show_in"),
            "n_terminus_inner_side": bool(side_in),
            "agrees": bool(side_in) == bool(rec.get("topology_show_in")),
        }
    segs = []
    for chain, a, b in segments(rec):
        xyz = [ca_opm[k][1] for k in ca_opm if k[0] == chain and a <= k[1] <= b]
        row = {"chain": chain, "from": a, "to": b, "ca": len(xyz)}
        if len(xyz) >= MIN_SEG_CA:
            xyz = np.array(xyz)
            lo, hi = depth_range(xyz, np.array([0.0, 0.0, s]), s * zc)  # against OPM's own planes, OPM frame
            row["crosses_opm"] = lo < -_MIN_CROSSING and hi > _MIN_CROSSING
            drawn_xyz = (xyz - t) @ r  # into the drawn frame: R^T (x - t)
            row["xyz_drawn_centroid"] = drawn_xyz.mean(0).round(3).tolist()
            if mem is not None:
                lo, hi = depth_range(drawn_xyz, mem.normal, mem.centre)
                row["crosses_estimate"] = lo < -_MIN_CROSSING and hi > _MIN_CROSSING
                row["depth_estimate"] = [round(lo, 2), round(hi, 2)]
        segs.append(row)
    e["segments"] = segs
    e["estimate"] = None
    if mem is not None:
        cosang = float(np.dot(mem.normal, n_opm))
        tm = np.array([x["xyz_drawn_centroid"] for x in segs if "xyz_drawn_centroid" in x])
        g = tm.mean(0) if len(tm) else bb.ca.mean(0)
        d_est = float(g @ mem.normal - mem.centre)
        d_opm = float(g @ n_opm - c_opm)
        e["estimate"] = {
            "source": mem.source,
            "normal": [round(float(x), 6) for x in mem.normal],
            "centre": round(float(mem.centre), 4),
            "half": mem.half,
            "thickness": 2 * mem.half,
            "sided": bool(mem.sided),
            "barrel": mem.defining is not None,
            "normal_angle_deg": round(float(np.degrees(np.arccos(min(1.0, abs(cosang))))), 3),
            "normal_dot": round(cosang, 6),
            "thickness_diff": round(2 * mem.half - e["opm_thickness"], 3),
            "midplane_offset": round(abs(d_est - np.sign(cosang) * d_opm), 3),
            "side_correct": (cosang > 0) if mem.sided else None,
        }
    e["opm_normal_drawn"] = [round(float(x), 6) for x in n_opm]
    e["opm_centre_drawn"] = round(float(c_opm), 4)
    e["attribution"] = attribution(path, bb)
    e["seconds"] = round(time.time() - t0, 1)
    return e


def main() -> int:
    sys.path.insert(0, str(TREE / "src"))
    import foldmap

    where = Path(foldmap.__file__).resolve()
    assert str(where).startswith(str(TREE.resolve())), f"foldmap imported from {where}, not build/tree"
    stamp = run_stamp(__file__)  # the code as the step starts (plan 28 Task 1)
    sel = json.loads((OPM / "selection.json").read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(min(12, len(sel["taken"]))) as ex:
        rows = list(ex.map(one, sel["taken"]))
    for e in rows:
        est = e.get("estimate")
        print(
            e["id"],
            e.get("error", ""),
            e.get("fit", {}).get("rmsd"),
            "refused" if est is None else (est["normal_angle_deg"], est["thickness_diff"], est["side_correct"]),
            flush=True,
        )
    extra = []
    from foldmap.fetch import fetch
    from foldmap.io import load_backbone

    for pid in EXTRA_ATTRIBUTION:
        try:
            path = Path(fetch(pid))
            extra.append({"id": pid, "attribution": attribution(path, load_backbone(path, "auto"))})
        except Exception as err:  # noqa: BLE001
            extra.append({"id": pid, "error": f"{type(err).__name__}: {err}"[:200]})
        print(extra[-1], flush=True)
    doc = {
        "generator": "scripts/paper/membrane_opm.py",
        "source_stamp": stamp,
        "foldmap_from": str(where),
        "selection": sel,
        "entries": rows,
        "attribution_extra": extra,
    }
    (OUT / "results.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8")
    print(f"wrote {OUT / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
