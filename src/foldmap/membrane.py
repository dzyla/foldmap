# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Where the lipid bilayer is: exactly from OPM's dummy atoms when the file has them (planes, thickness and
sides), otherwise estimated from the structure. An outer-membrane beta-barrel (a closed barrel whose outside is
oily) gives the normal as its axis, a 24 Å slab slid along it to its most hydrophobic outward-facing residues and
Tyr/Trp girdles, and the sides from its loops (long outside, short turns in the periplasm). Otherwise the long
helices and strands give the normal, a 30 Å slab is slid to the most hydrophobic, least charged position, and the
cytoplasmic side comes from the positive-inside rule (more Lys/Arg on the loops of the elements crossing it)."""

from __future__ import annotations

import dataclasses
import gzip
from dataclasses import dataclass

import gemmi
import numpy as np

from .dssp import assign_dssp
from .layout import _snap, _unit
from .model import SSE, Backbone, Bridge, homolog_chains
from .sheets import build_sheets
from .ss import _fit_axis, build_sses

HALF = 15.0  # Å, half the hydrophobic thickness of a typical bilayer
_HYDROPHOBIC = {"ALA", "ILE", "LEU", "MET", "PHE", "VAL", "TRP", "CYS"}
_POLAR = {"LYS", "ARG", "ASP", "GLU", "ASN", "GLN", "HIS"}
_MIN_HYDROPHOBIC, _MAX_POLAR = 0.55, 0.15  # a slab this oily and this uncharged is membrane
# ... or this polar when this many long elements cross it (6JO6, photosystem I: 0.174 polar, 46 crossing; cofactor
# pockets in the slab; the sample's soluble entries with 10 or more crossing are at most 0.39 hydrophobic)
_MAX_POLAR_CROSSED, _MANY_CROSSING = 0.20, 10
_SIDED_RATIO, _SIDED_DIFF = 2.0, 3  # Lys/Arg counts must differ this clearly to call the sides
_MIN_CROSSING = 10.0  # Å an element must reach beyond the centre on both sides to cross the membrane
# The positive-inside rule's loops: the chain this many residues past either end of an element that crosses the
# membrane (crosses its mid-plane, spanning half its thickness), their Lys/Arg counted up to this far beyond its
# face. Soluble domains beyond are left out (plan 26, 32KS ComEC: the periplasmic OB and β-lactamase-like domains'
# own Lys/Arg near the face called their side the cytoplasm, 49 to 23; its loops say 4.5 to 9.5 per protomer).
# Against OPM on the manuscript's membrane campaign (17 helical entries): 15 sided, 14 right (6QFA wrong, as
# before), where every residue near the face sided 11, 10 right; no entry's side worse (plan 26 Task 1 report).
_FLANK, _FLANK_REACH = 20, 15.0

# Outer-membrane beta-barrels. OPM gives their hydrophobic half thickness as 11.7-12.6 Å (FimD 3RFZ 11.7, OmpX
# 1QJ8 11.8, OmpA 1BXW 12.6, OmpF 2OMF 12.6): an outer membrane is thinner than a cytoplasmic one.
BARREL_HALF = 12.0  # Å
_BARREL_STRANDS, _BARREL_STRAND_LEN = 8, 7  # a closed sheet with this many strands this long is a barrel
# Lipid-facing residues of outer-membrane barrels: the aliphatics and aromatics, and Gly, common there.
_OUTWARD_HYDROPHOBIC = _HYDROPHOBIC | {"TYR", "GLY"}
_CHARGED = {"LYS", "ARG", "ASP", "GLU"}
_GIRDLE = 4.0  # Å either side of a membrane face where an outward Tyr/Trp counts as its aromatic girdle
# Outward residues in the slab, measured: outer-membrane barrels 0.76-1.00 oily and 0.00-0.07 charged (9FWZ, 1BXW,
# 1QJ8, 2OMF, 2POR, 1QD6, 1THQ, 2F1C, 1I78, 1FEP, 2GSK, 3EMN, 7AHL, 1EK9, 1K24, 1TLY, 2MPR, 1QD5, 1P4T, 3PRN);
# soluble barrels 0.35-0.44 and 0.26-0.30 (GFP 1EMA and 1GFL, retinol-binding protein 1RBP).
_MIN_OUTWARD_HYDROPHOBIC, _MAX_OUTWARD_CHARGED = 0.6, 0.15
# Sides of a one-chain barrel: its extracellular loops are long, its periplasmic turns short. Against OPM's sides
# (19 of the barrels above; 2F1C's chains do not match OPM's), the positive-inside rule was wrong on 16, 10 of them
# called clear; the loop residues on each face, each loop capped (FimD's 88-residue plug), call 17 right and 13
# clearly (_SIDED_RATIO and this many residues apart), none of those wrongly.
_LOOP_CAP, _SIDED_LOOPS = 20, 10


@dataclass
class Membrane:
    normal: np.ndarray  # unit vector pointing to the outside (non-cytoplasmic side)
    centre: float  # position of the mid-plane along the normal (Å)
    half: float  # half thickness (Å)
    source: str  # 'opm' (dummy atoms in the file) or 'estimate'
    sided: bool = False  # True when the positive-inside rule clearly says which side is the cytoplasm
    # residues (indices) of the elements that define it, a beta-barrel's strands: only those give a page row its band
    # (a subunit in the barrel's lumen spans the slab without being in the membrane); None: any element crossing it
    defining: frozenset[int] | None = None

    def depth(self, xyz: np.ndarray) -> np.ndarray:
        """Signed distance from the mid-plane, positive outside."""
        return np.asarray(xyz) @ self.normal - self.centre


def _dummies_after_end(path) -> list[tuple[str, float]]:
    """DUM atoms (name, z) of a PDB file's records after its first END record. Many OPM files put their dummy atoms
    there, where gemmi stops reading; nothing else after END is taken, so the model is unchanged."""
    p = str(path)
    try:
        opener = gzip.open if p.endswith(".gz") else open
        with opener(p, "rt", errors="replace") as fh:
            lines = fh.read().splitlines()
    except (OSError, EOFError):
        return []
    out, ended = [], False
    for line in lines:
        if not ended:
            ended = line[:6].rstrip() == "END"
        elif line[:6] in ("ATOM  ", "HETATM") and line[17:20] == "DUM":
            try:
                out.append((line[12:16].strip(), float(line[46:54])))
            except ValueError:
                continue
    return out


def _opm(path) -> Membrane | None:
    """The bilayer from OPM/PPM dummy atoms. Those files are in the membrane frame (normal along z) and mark the
    hydrophobic boundaries with DUM atoms: 'O' on the outer (extracellular, periplasmic, lumenal) boundary and
    'N' on the inner (cytoplasmic) one, so they give the sides as well as the planes."""
    st = gemmi.read_structure(str(path))
    dummies = [(a.name, a.pos.z) for ch in st[0] for r in ch if r.name == "DUM" for a in r]
    dummies += _dummies_after_end(path)
    if not dummies:
        return None
    inner = [z for name, z in dummies if name == "N"]
    outer = [z for name, z in dummies if name == "O"]
    if inner and outer and abs(np.mean(outer) - np.mean(inner)) > 1.0:
        zi, zo = float(np.mean(inner)), float(np.mean(outer))
        normal = np.array([0.0, 0.0, 1.0 if zo > zi else -1.0])
        return Membrane(normal, (zi + zo) / 2 * normal[2], abs(zo - zi) / 2, "opm", sided=True)
    z = np.array([z for _, z in dummies])
    return Membrane(np.array([0.0, 0.0, 1.0]), 0.0, float(np.abs(z).max()), "opm")


def find_membrane(path, bb: Backbone) -> Membrane | None:
    """The membrane, or None for a soluble protein."""
    found, sses = _opm(path), None
    if found is None:
        dssp = assign_dssp(bb)
        sses = build_sses(bb, dssp.ss, bridges=dssp.bridges)
        barrel = _barrel_estimate(bb, sses, dssp.bridges)
        if barrel is not None:
            return barrel  # sided by its loops, never by the positive-inside rule
        found = _estimate(bb, sses)
    if found is None:
        return None
    if found.sided:  # OPM's own sides, when its dummy atoms give them
        return found
    if sses is None:  # OPM's planes without its sides
        dssp = assign_dssp(bb)
        sses = build_sses(bb, dssp.ss, bridges=dssp.bridges)
    return _orient(found, bb, sses)


def _barrels(sses: list[SSE], bridges: list[Bridge]) -> list[list[SSE]]:
    """The long strands of each closed beta-barrel (sheets.build_sheets closes them), largest barrel first."""
    out = []
    for sheet in build_sheets(sses, bridges):
        strands = [s for s in sheet.strands if len(s) >= _BARREL_STRAND_LEN]
        if sheet.closed and len(strands) >= _BARREL_STRANDS:
            out.append(sorted(strands, key=lambda s: s.start))
    return sorted(out, key=lambda strands: (-len(strands), strands[0].start))


def _pseudo_cb(xyz: np.ndarray) -> np.ndarray:
    """Where each residue's CB is (or would be, for Gly) from its N, CA and C: ideal tetrahedral geometry."""
    n, ca, c = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    b, cc = ca - n, c - ca
    return -0.58273431 * np.cross(b, cc) + 0.56802827 * b - 0.54067466 * cc + ca


def _cylinder(ca: np.ndarray, guess: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The axis (unit vector) and a point on it of the cylinder through these CA atoms (least squares, from
    `guess`)."""
    from scipy.optimize import least_squares

    e1 = np.cross(guess, [1.0, 0.0, 0.0] if abs(guess[0]) < 0.9 else [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(guess, e1)
    mid = ca.mean(axis=0)

    def unpack(x):
        axis = guess + x[0] * e1 + x[1] * e2
        return axis / np.linalg.norm(axis), mid + x[2] * e1 + x[3] * e2

    def residuals(x):
        axis, point = unpack(x)
        v = ca - point
        return np.linalg.norm(v - np.outer(v @ axis, axis), axis=1) - x[4]

    v = ca - mid
    radius = np.linalg.norm(v - np.outer(v @ guess, guess), axis=1).mean()
    return unpack(least_squares(residuals, [0.0, 0.0, 0.0, 0.0, radius]).x)


def _barrel_slab(bb: Backbone, strands: list[SSE]) -> tuple[np.ndarray, float, float, float]:
    """The barrel's axis (on the layout grid), the centre along it of the BARREL_HALF slab holding the most oily,
    least charged outward-facing residues with Tyr/Trp girdles at its faces, and the oily and charged fractions of
    the outward residues in that slab. A residue faces out when its (pseudo-)CB points away from the axis."""
    idx = np.concatenate([np.arange(s.start, s.end + 1) for s in strands])
    ca = bb.ca[idx]
    top = np.linalg.eigh(sum(len(s) * np.outer(s.axis, s.axis) for s in strands))[1][:, -1]
    guess = sum(len(s) * np.sign(s.axis @ top) * s.axis for s in strands)  # up and down strands made parallel
    axis, point = _cylinder(ca, guess / np.linalg.norm(guess))
    normal = _unit(axis)
    v = ca - point
    radial = v - np.outer(v @ normal, normal)
    out = ((_pseudo_cb(bb.xyz[idx]) - ca) * radial).sum(axis=1) > 0
    names = [bb.labels[k].name for k in idx]
    oily = np.array([n in _OUTWARD_HYDROPHOBIC for n in names]) & out
    charged = np.array([n in _CHARGED for n in names]) & out
    aromatic = np.array([n in ("TYR", "TRP") for n in names]) & out
    z = ca @ normal
    best = None
    for c in np.arange(np.floor(z.min()), z.max() + 1.0, 1.0):
        inside = np.abs(z - c) < BARREL_HALF
        if (out & inside).sum() < 10:
            continue
        girdle = aromatic & (np.abs(np.abs(z - c) - BARREL_HALF) < _GIRDLE)
        score = oily[inside].sum() - 1.5 * charged[inside].sum() + girdle.sum()
        if best is None or score > best[0]:
            best = (score, c, out & inside)
    if best is None:
        return normal, 0.0, 0.0, 1.0
    _, centre, facing = best
    n = facing.sum()
    return normal, float(_snap(centre)), float(oily[facing].sum() / n), float(charged[facing].sum() / n)


def _barrel_estimate(bb: Backbone, sses: list[SSE], bridges: list[Bridge]) -> Membrane | None:
    """The outer membrane around the first closed beta-barrel whose outside is oily and uncharged and whose strands
    cross the slab."""
    for strands in _barrels(sses, bridges):
        normal, centre, oily, charged = _barrel_slab(bb, strands)
        crossing = [
            s
            for s in strands
            if (bb.ca[s.start : s.end + 1] @ normal).min() < centre - _MIN_CROSSING
            and (bb.ca[s.start : s.end + 1] @ normal).max() > centre + _MIN_CROSSING
        ]
        if oily >= _MIN_OUTWARD_HYDROPHOBIC and charged <= _MAX_OUTWARD_CHARGED and crossing:
            barrel = frozenset(k for s in strands for k in range(s.start, s.end + 1))
            return _barrel_sides(Membrane(normal, centre, BARREL_HALF, "estimate", defining=barrel), bb, strands)
    return None


def _barrel_sides(m: Membrane, bb: Backbone, strands: list[SSE]) -> Membrane:
    """Point the normal to the face with the longer loops between consecutive barrel strands (the outside); sided
    only for a one-chain barrel whose two faces' loops differ clearly."""
    up = down = 0
    for a, b in zip(strands, strands[1:]):
        if a.chain == b.chain:
            loop = min(b.start - a.end - 1, _LOOP_CAP)
            if m.depth(bb.ca[[a.end, b.start]]).mean() > 0:
                up += loop
            else:
                down += loop
    one_chain = len({s.chain for s in strands}) == 1
    sided = one_chain and max(up, down) >= _SIDED_RATIO * min(up, down) and abs(up - down) >= _SIDED_LOOPS
    if down > up:
        return dataclasses.replace(m, normal=-m.normal, centre=-m.centre, sided=sided)
    return dataclasses.replace(m, sided=sided)


def _estimate(bb: Backbone, sses: list[SSE]) -> Membrane | None:
    """The best-scoring HALF slab across the long elements' common axis, if it is membrane (_membrane_like). When no
    long element spans that slab face to face (_spans), so the view would find no transmembrane element to draw
    its band by, the best-scoring slab that one spans and that is membrane takes its place (random-100 held-out
    6MJ2, TRPM2: the dense cytoplasmic domain pulls the summed score 4-5 Å below the S1-S4 helices)."""
    long = [s for s in sses if (s.kind == "H" and len(s) >= 16) or (s.kind == "E" and len(s) >= 8)]
    if not long:
        return None
    normal = np.linalg.eigh(sum(len(s) * np.outer(s.axis, s.axis) for s in long))[1][:, -1]
    z = bb.ca @ normal
    hyd = np.array([l.name in _HYDROPHOBIC for l in bb.labels])
    pol = np.array([l.name in _POLAR for l in bb.labels])
    chain = np.array([l.chain for l in bb.labels])
    slabs = []  # (score, centre, inside), bottom up
    for c in np.arange(z.min(), z.max() + 1.0, 1.0):
        inside = np.abs(z - c) < HALF
        if inside.sum() < 20:
            continue
        slabs.append((hyd[inside].sum() - 1.5 * pol[inside].sum(), c, inside))
    if not slabs:
        return None
    best = max(slabs, key=lambda t: t[0])  # the first of equal scores, as before
    if not _spans(long, z, best[1]):
        spanned = [t for t in slabs if _spans(long, z, t[1]) and _membrane_like(long, z, t, hyd, pol, chain)]
        best = max(spanned, key=lambda t: t[0]) if spanned else best
    if not _membrane_like(long, z, best, hyd, pol, chain):
        return None
    return Membrane(normal, float(best[1]), HALF, "estimate")


FACE_REACH = 2.0  # Å: a long element whose CA run comes this near both faces spans the slab (_spans, the estimate)


def crosses(depth: np.ndarray, half: float) -> bool:
    """Does an element whose CA atoms lie at these depths cross the membrane: its run crosses the mid-plane and spans
    at least half the thickness (`half`)? The one rule the membrane view lays its transmembrane elements out by
    (cli._transmembrane, _face_medians), A18 judges them by (figcheck.membrane_view), foldmap data reports
    (crosses_midplane) and the positive-inside rule takes their loops from (_loops). Plan 26: the rule before (the run
    within FACE_REACH of both faces) counted 2 of 32KS's 11 transmembrane helices per chain, their DSSP runs spanning
    14-25 Å of a 30 Å slab; OPM's membranes are as thick around helices as short (6DQN, 6BO8, 6MJ2). The predicate is
    shared, the sets it is put to are not: A18 judges helices only; the layout keeps to a barrel's strands
    (Membrane.defining) besides."""
    seg = np.asarray(depth)
    return len(seg) >= 2 and float(seg.min()) < 0.0 < float(seg.max()) and float(seg.max() - seg.min()) >= half


SOLUBLE_BEYOND = 5.0  # Å: an element whose CA centre lies this far beyond a face is soluble, not interfacial (R4)
FLAT_TILT = (
    45.0  # degrees: a helix crossing no membrane whose axis is tilted this far or more from the normal lies flat
)


# Plan 26 Task 5 ruling T5-2: at a face a helix lies flat from FACE_FLAT_TILT. The interfacial zone of a face is
# INTERFACE either side of it: the bilayer's interfacial region, its headgroups and water, is about 15 Å thick on
# each side of the hydrocarbon core (White & Wimley 1999, Annu Rev Biophys Biomol Struct 28:319). 32KS's NS2 and NS3
# (A:222-233, A:241-250), tilted 44.6° and 38.6° with their centres 2.9 and 6.1 Å beyond the outer face, stood
# upright 12-36 Å above it; ComEC's Fig. 1C draws them flat on the face.
FACE_FLAT_TILT = 35.0  # degrees
INTERFACE = 7.5  # Å either side of a face


def lies_flat(axis: np.ndarray, normal: np.ndarray, tilt: float = FLAT_TILT) -> bool:
    """Is a helix with this axis drawn lying along the membrane, not standing across it: tilted `tilt` (FLAT_TILT)
    or more from the normal (plan 26 Task 2 ruling: the nearer of upright and flat, by the angle in 3D, not as
    projected; 32KS's amphipathic NS1, A:200-219, tilted 71°, pointed into the page and was snapped upright)? The test
    is on the axis' component along the normal, on the grid, so noise never decides it."""
    along = abs(float(_snap(np.dot(_unit(np.asarray(axis, float)), normal))))
    return along <= float(_snap(np.cos(np.radians(tilt))))


def at_face(depth: np.ndarray, half: float) -> bool:
    """Does an element whose CA atoms lie at these depths lie in a face's interfacial zone: its CA centre within
    INTERFACE of either face (ruling T5-2)?"""
    return abs(abs(float(np.mean(depth))) - half) <= INTERFACE


def lies_flat_here(axis: np.ndarray, depth: np.ndarray, m: Membrane) -> bool:
    """Is a helix with this axis, its CA atoms at these depths, drawn lying flat in a membrane view: crossing no
    membrane and tilted FLAT_TILT or more from the normal, or FACE_FLAT_TILT or more when it lies in a face's
    interfacial zone (at_face; ruling T5-2); a crossing helix always stands (ruling R7)."""
    if crosses(depth, m.half):
        return False
    return lies_flat(axis, m.normal, FACE_FLAT_TILT if at_face(depth, m.half) else FLAT_TILT)


# Ruling T5-1: in a membrane view two helix runs of one chain at most JOIN_GAP non-helix residues apart are joined
# only when both lie flat or the joined helix crosses the membrane. Joined as those words allow, 199 pairs of the
# membrane campaign's 21 views would have been: hairpins (32KS's NS1 and NS2, A:200-219 and A:222-233, their axes
# 151° apart; 6UW4's A:119-129 and A:132-145, 163°), and the cytoplasmic CH1 into TM11 (A:501-505, A:508-526), TM6
# into TM7 of 32KS chain B. So the flat pieces join only when one goes on as the other (their axes within
# JOIN_ANGLE: NS7's 37°), and a crossing helix only where copies of the protein disagree (ss._copy_breaks: joined, not
# cut, when the joined helix crosses).
JOIN_GAP = 2  # residues
JOIN_ANGLE = 45.0  # degrees between the two pieces' axes, N to C


def _apart(bb: Backbone, x: tuple[int, int], y: tuple[int, int]) -> bool:
    """Are two helix runs (first and last residue index each) of one chain, 1 to JOIN_GAP residues apart in an
    unbroken chain?"""
    gap = y[0] - x[1] - 1
    if bb.labels[x[0]].chain != bb.labels[y[0]].chain or not 1 <= gap <= JOIN_GAP:
        return False
    return bool(np.asarray(bb.prev[x[1] + 1 : y[0] + 1]).all())


def joins(bb: Backbone, m: Membrane, x: tuple[int, int], y: tuple[int, int]) -> bool:
    """Ruling T5-1, item 1: in a membrane view two helix runs of one chain (first and last residue index each), at
    most JOIN_GAP non-helix residues apart in an unbroken chain, are drawn as one helix when both lie flat
    (lies_flat_here) and one goes on as the other (their axes within JOIN_ANGLE): 32KS's NS7, A:436-461, which DSSP
    breaks at 447-448 in chain A only, its pieces 37° apart."""
    if not _apart(bb, x, y):
        return False
    d = m.depth(bb.ca)
    ax, ay = (_fit_axis(bb.ca[a : b + 1]) for a, b in (x, y))
    if float(_snap(np.dot(ax, ay))) < float(_snap(np.cos(np.radians(JOIN_ANGLE)))):
        return False
    return all(lies_flat_here(axis, d[a : b + 1], m) for axis, (a, b) in ((ax, x), (ay, y)))


def crosses_joined(bb: Backbone, m: Membrane, x: tuple[int, int], y: tuple[int, int]) -> bool:
    """Ruling T5-1, item 3: does a helix DSSP breaks in two (x, y: first and last residue index each, at most JOIN_GAP
    apart) cross the membrane joined? Copies of the protein that disagree are then joined, not cut (ss._copy_breaks
    `whole`): 32KS chain B's TM8, B:349-371, broken at 364, cut chain A's whole TM8 in two, its C-terminal piece lying
    flat under the inner face."""
    return _apart(bb, x, y) and crosses(m.depth(bb.ca)[x[0] : y[1] + 1], m.half)


def _spans(long: list[SSE], z: np.ndarray, centre: float) -> bool:
    """Does a long element's CA run reach both faces of the HALF slab at `centre` (within FACE_REACH)? The estimate's
    own test (random-100 held-out 6MJ2), kept as it was measured when the view and A18 took `crosses` (plan 26): the
    planes of the manuscript's membrane campaign do not move."""
    return any(
        z[s.start : s.end + 1].min() <= centre - HALF + FACE_REACH
        and z[s.start : s.end + 1].max() >= centre + HALF - FACE_REACH
        for s in long
    )


def _membrane_like(long, z, slab, hyd, pol, chain) -> bool:
    """Is this slab membrane: crossed by a long element (_MIN_CROSSING beyond its centre both ways), and oily and
    uncharged enough over the residues of the chains that cross it? A chain with no element across the slab is a
    soluble partner reaching into it, not facing lipid (random-100 held-out 7VVK, PTH1R-Gs: the hormone in the
    receptor's pocket and Gs's α5 in its cytoplasmic cavity make the slab 0.189 polar; the receptor's own residues
    are 0.149)."""
    _, centre, inside = slab
    crossing = [
        s
        for s in long
        if z[s.start : s.end + 1].min() < centre - _MIN_CROSSING
        and z[s.start : s.end + 1].max() > centre + _MIN_CROSSING
    ]
    if not crossing:
        return False
    held = inside & np.isin(chain, sorted({s.chain for s in crossing}))
    polar = _MAX_POLAR_CROSSED if len(crossing) >= _MANY_CROSSING else _MAX_POLAR
    return hyd[held].mean() >= _MIN_HYDROPHOBIC and pol[held].mean() <= polar


def _loops(m: Membrane, bb: Backbone, sses: list[SSE]) -> np.ndarray:
    """The residues (a mask) within _FLANK of either end of an element that crosses the membrane (membrane.crosses),
    along the element's chain and up to the first gap in its numbering (an unmodelled stretch: the residues past it
    are no part of this loop)."""
    d = m.depth(bb.ca)
    chain = [l.chain for l in bb.labels]
    seq = [l.seq for l in bb.labels]
    out = np.zeros(len(bb), dtype=bool)
    for s in sses:
        if not crosses(d[s.start : s.end + 1], m.half):
            continue
        for step, j in ((-1, s.start - 1), (1, s.end + 1)):
            for _ in range(_FLANK):
                if not (0 <= j < len(bb)) or chain[j] != s.chain or abs(seq[j] - seq[j - step]) > 1:
                    break
                out[j] = True
                j += step
    return out


def _orient(m: Membrane, bb: Backbone, sses: list[SSE]) -> Membrane:
    """Point the normal outward: the side with fewer Lys/Arg on the loops of the membrane-crossing elements, within
    _FLANK_REACH of its face (positive-inside rule, _loops), copies of one protein counted once (their mean: KcsA's
    four copies of 2 to 1 are 2 to 1, not 8 to 4); on a tie, fewer Lys/Arg of any residue just beyond the membrane.
    Sided when the loops' counts differ clearly, or when every residue's do and lean the same way (6BO8, TRPV6: its
    loops 8 to 5 per copy, its ankyrin domains 72 to 28)."""
    d = m.depth(bb.ca)
    basic = np.array([l.name in ("LYS", "ARG") for l in bb.labels])
    band = lambda lo, hi, on: basic & on & (d > lo) & (d < hi)  # noqa: E731
    clear = lambda a, b: max(a, b) >= _SIDED_RATIO * min(a, b) and abs(a - b) >= _SIDED_DIFF  # noqa: E731
    chain = np.array([l.chain for l in bb.labels])
    copies: dict[str, list[str]] = {}
    for c, (first, _) in homolog_chains(bb).items():
        copies.setdefault(first, []).append(c)

    def per_protein(on: np.ndarray) -> float:  # summed over the proteins, each its copies' mean
        return float(sum(np.mean([on[chain == c].sum() for c in cs]) for cs in copies.values()))

    loops = _loops(m, bb, sses)
    up = per_protein(band(m.half, m.half + _FLANK_REACH, loops))
    down = per_protein(band(-m.half - _FLANK_REACH, -m.half, loops))
    every = np.ones(len(bb), dtype=bool)
    up_all = int(band(m.half - 5, m.half + 15, every).sum())
    down_all = int(band(-m.half - 15, -m.half + 5, every).sum())
    flip = up > down if up != down else up_all > down_all
    sided = clear(up, down) or (clear(up_all, down_all) and (up_all > down_all) == flip)
    if flip:
        return dataclasses.replace(m, normal=-m.normal, centre=-m.centre, sided=sided)
    return dataclasses.replace(m, sided=sided)
