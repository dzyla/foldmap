# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Checks of a nucleic-acid figure against the deposited coordinates.

verify_contacts: every contact bead a figure draws (a nucleotide touched by protein, in that chain's colour)
must be backed by real heavy atoms.
shape, contacts, render_checks: fidelity of the drawing to the structure (stem order and shape, protein contacts,
rendering), judged against THRESHOLDS by judge(). Used by tests/test_na_validation.py on the fixtures and by
scripts/na_validate.py on the benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

import gemmi
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from .frame import Frame
from .io import CONTACT
from .layout import SCALE, stem_bodies, stem_pieces, stems_at
from .na_curve import drawn_pieces, piece_bodies
from .na_paths import _BEAD, among, crosses_body, graze, in_body
from .nucleic import stem_axis
from .route import clip_len


@dataclass(frozen=True)
class Pair:
    """The closest heavy-atom pair between a nucleotide and a protein residue."""

    residue: str  # "A:ARG42"
    atom: str  # protein atom, "NH1"
    nucleotide: str  # "C:DG7"
    nt_atom: str  # "N7"
    dist: float

    def __str__(self) -> str:
        return f"{self.residue} {self.atom} – {self.nucleotide} {self.nt_atom} {self.dist:.1f} Å"


@dataclass(frozen=True)
class Bead:
    """One bead as drawn and what the coordinates say of it. problem: "" (verified); "unmatched" (the nucleotide or
    the residue it names is not in the model); "false" (no heavy-atom pair within CONTACT, or the residue is of
    another chain than the bead's); "wrong colour" (not the figure's colour for the residue's chain);
    "not first" (the bead does not name the first toucher the drawing rule picks)."""

    key: tuple[int, int]  # (strand, nucleotide) in the layout
    nucleotide: str  # "C:DG7"
    chain: str  # protein chain whose colour the bead carries
    colour: str
    residue: str  # the protein residue named as the toucher, "A:ARG42"
    pair: Pair | None  # closest atom pair of that residue and the nucleotide, when within CONTACT
    problem: str = ""

    @property
    def verified(self) -> bool:
        return not self.problem


@dataclass(frozen=True)
class Missed:
    """A protein residue touching a nucleotide in the coordinates that has no bead naming it."""

    nucleotide: str
    pair: Pair
    reason: str  # "missing" (a defect), or why none is expected: "focus", "no backbone", "not drawn", "first toucher",
    # "cropped"

    @property
    def counted(self) -> bool:
        return self.reason == "missing"


@dataclass(frozen=True)
class Hydrog:
    """A `hydrog` record of the file's _struct_conn: the deposited hydrogen bonds."""

    a: str  # "A:DG7 N1"
    b: str
    dist: float | None


@dataclass
class ContactReport:
    beads: list[Bead] = field(default_factory=list)
    missed: list[Missed] = field(default_factory=list)  # every unbeaded contact of a drawn nucleotide
    cropped: list[Missed] = field(default_factory=list)  # touched nucleotides the figure leaves out (crop, focus)
    hydrog: list[Hydrog] = field(default_factory=list)

    @property
    def verified(self) -> list[Bead]:
        return [b for b in self.beads if b.verified]

    def with_problem(self, problem: str) -> list[Bead]:
        return [b for b in self.beads if b.problem == problem]

    @property
    def false(self) -> list[Bead]:
        return self.with_problem("false")

    @property
    def unmatched(self) -> list[Bead]:
        return self.with_problem("unmatched")

    @property
    def wrong_colour(self) -> list[Bead]:
        return self.with_problem("wrong colour")

    @property
    def not_first(self) -> list[Bead]:
        return self.with_problem("not first")

    @property
    def missing(self) -> list[Missed]:
        return [m for m in self.missed if m.counted]

    @property
    def excluded(self) -> list[Missed]:
        return [m for m in self.missed if not m.counted] + self.cropped


def _tag(chain: str, name: str, seq: int, icode: str) -> str:
    return f"{chain}:{name}{seq}{icode.strip()}"


def _is_protein(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return info is not None and info.is_amino_acid()


def _is_nucleic(res) -> bool:
    """A nucleotide not drawn: tabulated as one, or a polymer residue with a sugar (C1', C4', O4')."""
    info = gemmi.find_tabulated_residue(res.name)
    if info is not None and (info.is_nucleic_acid() or info.is_amino_acid()):
        return info.is_nucleic_acid()
    return res.entity_type == gemmi.EntityType.Polymer and all(
        res.find_atom(n, "*") is not None for n in ("C1'", "C4'", "O4'")
    )


def deposited_hydrog(path) -> list[Hydrog]:
    """The hydrogen-bond records (_struct_conn, conn_type_id hydrog) of the file, when it has any."""
    out = []
    for con in gemmi.read_structure(str(path)).connections:
        if con.type == gemmi.ConnectionType.Hydrog:
            a, b = (
                f"{p.chain_name}:{p.res_id.name}{p.res_id.seqid.num}{p.res_id.seqid.icode.strip()} {p.atom_name}"
                for p in (con.partner1, con.partner2)
            )
            out.append(Hydrog(a, b, con.reported_distance or None))
    return out


def _heavy(residues):
    """(xyz (n, 3), owner (n,), names) of the heavy atoms (not H or D) of these gemmi residues."""
    xyz, owner, names = [], [], []
    for k, res in enumerate(residues):
        for atom in res:
            if atom.element.name not in ("H", "D"):
                xyz.append((atom.pos.x, atom.pos.y, atom.pos.z))
                owner.append(k)
                names.append(atom.name)
    return np.array(xyz).reshape(-1, 3), np.array(owner, int), names


def verify_contacts(model, drawn, layout, bb, *, path=None) -> ContactReport:
    """Check the contact beads of a drawn figure against `model`, the gemmi model (or read_model result) it was
    built from. `drawn` is what draw() leaves on the figure (`fig.foldmap_beads`, a render.ContactBeads: the beads
    per chain and the colour function the figure used); `layout` and `bb` are the layout and Backbone it was
    drawn from.
    Shared with io.load_nucleic: only the cutoff io.CONTACT and the input model. Independent: the nucleotides are
    the figure's drawn ones, each fetched from the model by its chain, number and name (a bead that cannot be is
    "unmatched", never dropped); the contacts come from a scipy cKDTree over the heavy atoms (hydrogens and
    deuteriums left out by element), distances taken in numpy, not from gemmi's NeighborSearch.
    A bead is verified when the residue it names lies in the bead's chain with a heavy-atom pair within CONTACT of
    its nucleotide, its colour is drawn.colour(residue index, that chain), and it names the first (lowest
    Backbone index) of the residues touching the nucleotide that the figure colours. Every other protein residue
    touching a drawn nucleotide is in `missed`: reason "missing" is a defect; "focus" (a chain outside
    layout.focus_chains), "no backbone" (a residue the figure does not hold), "not drawn" (a protein chain with no
    helix or strand, Layout.coil_chains: the figure draws nothing of it and says so, so no colour of its can mark a
    bead) and "first toucher" (the nucleotide's bead names an earlier residue) are expected. Touched nucleotides the
    figure does not draw are in `cropped`.
    path: also read the file's deposited hydrogen bonds into the report."""
    if isinstance(model, tuple):
        model = model[0]
    report = ContactReport(hydrog=deposited_hydrog(path) if path is not None else [])
    na = layout.nucleic
    if na is None or not na.strands or drawn is None:
        return report
    index = {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}
    by_id: dict[tuple, list] = {}  # (chain, seq, icode) -> residues of the model with that number
    for chain in model:
        for res in chain:
            by_id.setdefault((chain.name, res.seqid.num, res.seqid.icode.strip()), []).append(res)

    def fetch(ident, name):
        return next((r for r in by_id.get(ident, []) if r.name == name), None)

    # the nucleotides: those drawn (by label), and the model's others
    drawn_res: dict[tuple, tuple] = {}  # layout key -> (ident, residue or None)
    for s, strand in enumerate(na.strands):
        for i, lab in enumerate(strand.labels):
            ident = (lab.chain, lab.seq, lab.icode)
            drawn_res[(s, i)] = (ident, fetch(ident, lab.name), lab)
    drawn_ids = {ident for ident, _, _ in drawn_res.values()}
    others = [
        (ident, res) for ident, group in by_id.items() if ident not in drawn_ids for res in group if _is_nucleic(res)
    ]
    nts = [(key, ident, res) for key, (ident, res, _) in drawn_res.items() if res is not None]
    nts += [(None, ident, res) for ident, res in others]
    proteins = [
        (ident, res) for ident, group in by_id.items() for res in group if ident in index or _is_protein(res.name)
    ]
    nt_xyz, nt_owner, nt_names = _heavy([r for _, _, r in nts])
    pr_xyz, pr_owner, pr_names = _heavy([r for _, r in proteins])
    near: list[dict[tuple, Pair]] = [{} for _ in nts]  # per nucleotide: protein residue id -> closest Pair
    if len(nt_xyz) and len(pr_xyz):
        hits = cKDTree(nt_xyz).query_ball_tree(cKDTree(pr_xyz), CONTACT)
        for a, found in enumerate(hits):
            if not found:
                continue
            b = np.array(found)
            d = np.linalg.norm(pr_xyz[b] - nt_xyz[a], axis=1)
            _, ident, res = nts[nt_owner[a]]
            nt = _tag(ident[0], res.name, ident[1], ident[2])
            for j, dist in zip(b, d):
                if dist > CONTACT:
                    continue
                pid, pres = proteins[pr_owner[j]]
                if pid == ident:
                    continue
                best = near[nt_owner[a]].get(pid)
                if best is None or dist < best.dist:
                    near[nt_owner[a]][pid] = Pair(
                        _tag(pid[0], pres.name, pid[1], pid[2]), pr_names[j], nt, nt_names[a], float(dist)
                    )
    around = {key: near[n] for n, (key, _, _) in enumerate(nts) if key is not None}
    focus = layout.focus_chains
    undrawn = {c for c, _ in getattr(layout, "coil_chains", None) or []}
    colour_of = drawn.colour
    for chain, items in sorted(drawn.marks.items()):
        for key, r, colour in items:
            ident, res, lab = drawn_res[key]
            nt = _tag(ident[0], lab.name, ident[1], ident[2])
            who = bb.labels[r]
            named = (who.chain, who.seq, who.icode)
            tag = _tag(who.chain, who.name, who.seq, who.icode)
            if res is None or fetch(named, who.name) is None:
                report.beads.append(Bead(key, nt, chain, colour, tag, None, "unmatched"))
                continue
            pair = around[key].get(named)
            if pair is None or who.chain != chain:
                problem = "false"
            elif (colour_of(r, who.chain) or "").lower() != colour.lower():
                problem = "wrong colour"
            else:
                first = min(
                    (
                        index[p]
                        for p in around[key]
                        if p in index and (not focus or p[0] in focus) and colour_of(index[p], p[0])
                    ),
                    default=r,
                )
                problem = "" if first == r else "not first"
            report.beads.append(Bead(key, nt, chain, colour, tag, pair, problem))
    beaded = {key: bb.labels[r] for items in drawn.marks.values() for key, r, _ in items}
    for key, found in sorted(around.items()):
        ident, res, lab = drawn_res[key]
        nt = _tag(ident[0], lab.name, ident[1], ident[2])
        who = beaded.get(key)
        named = (who.chain, who.seq, who.icode) if who else None
        for pid, pair in sorted(found.items()):
            if pid == named:
                continue
            if pid not in index:
                reason = "no backbone"
            elif pid[0] in undrawn:
                reason = "not drawn"
            elif focus and pid[0] not in focus:
                reason = "focus"
            else:
                reason = "first toucher" if who else "missing"
            report.missed.append(Missed(nt, pair, reason))
    for n, (key, _, _) in enumerate(nts):
        if key is None:
            for pid, pair in sorted(near[n].items()):
                report.cropped.append(
                    Missed(pair.nucleotide, pair, "focus" if focus and pid[0] not in focus else "cropped")
                )
    return report


# pass/fail limits, set from the measured campaign (see the plan's Task 13); `None` disables a check
THRESHOLDS = {
    "pair_f1": 0.95,  # canonical + wobble pairs against RNApolis (final campaign: 0.96-1.00; started at 0.90)
    "stem_recall": 0.85,  # RNApolis stems matched one to one by a stem holding 2/3 of their pairs (0.857-1)
    # pairwise stem-centre distances, drawn against 3D, entries with >= 5 stems (0.601-0.96; thin: 1GID 0.601)
    "spearman": 0.60,
    # normalised RMSD of drawn stem centres after a similarity fit, no reflection (0.09-0.499; thin: 4LD0 0.499)
    "procrustes": 0.50,
    # angle and joint errors: 0 by construction on a plain frame: a regression guard on plain frames, vacuous as
    # fidelity (M2). The final campaign has 0.0 on every entry but 30QM, so the limit is one rounding step
    "angle_median": 0.05,  # degrees, drawn stem direction against its projected 3D axis
    "angle_max": 0.05,  # degrees, the worst stem; campaign 0.0 (stems drawn along their seen axes)
    "joint_error_max": 0.05,  # degrees, signed drawn bend against signed seen bend at the worst joint (campaign 0.0)
    "mirrored": 0,  # stems seen clearly (SIDE_SEEN) drawn with their arms swapped
    # elements touching a stem drawn nearer it than the rest (0.609-1.00). It does not respond to the contact pull:
    # tests/test_na_layout.py's displacement test is the guard
    "contact_auc": 0.60,
    "stems_nearer_frac": 1.0,  # contacted stems drawn nearer their touching elements than the median element
    "beads_missing": 0,  # protein residues touching a drawn nucleotide with no bead (verify_contacts' missing)
    "bead_colour_wrong": 0,  # beads not in their toucher's chain colour (verify_contacts' wrong colour)
    "beads_false": 0,  # beads no heavy-atom pair of their chain backs (verify_contacts' false)
    "beads_unmatched": 0,  # beads naming a nucleotide or residue the model lacks (verify_contacts' unmatched)
    "overlaps": 0,  # element or stem boxes overlapping
    "text_off_page": 0,
    "letter_overlaps": 0,  # base letters overlapping each other by more than a third of a glyph box
    "fallbacks": 0,  # strand or loop routes that fell back to a curve
    # targets; the fixture gate holds each fixture at its measured value until the campaign (Task 14 tightens)
    "na_crossings": 0,  # runs, joints and gaps crossing hairpin arcs, free ends or stem bodies
    "bead_overlaps": 0,  # contact beads overlapping one another
    "gap_label_collisions": 0,  # "... N nt ..." labels on each other, a stem, element, letter, bead or path
    "na_run_crossings": None,  # runs, joints and gaps crossing one another: reported, not judged
    "na_protein_crossings": None,  # nucleic-acid paths crossing protein loops: reported, not judged
    "seconds": 120.0,
}
_HIGH_IS_GOOD = ("pair_f1", "stem_recall", "spearman", "contact_auc", "stems_nearer_frac")
# What each measure is evidence of (plan 28 Task 2, ruling Q19), for the reports: a label only, no limit or verdict
# depends on it. "external": against RNApolis; "fidelity": how far the page keeps the 3D arrangement; "regression
# guard": zero by construction on an ordinary (plain) projection, since stems are drawn along their axes as seen
# (shape()), so they guard against a regression there and are no evidence of fidelity (on an unrolled frame they
# measure the unrolling); "consistency": beads against the model's own atoms; "readability": the page.
ZERO_BY_CONSTRUCTION = "zero by construction (regression guard)"
ROLE = {
    "pair_f1": "external",
    "stem_recall": "external",
    "spearman": "fidelity",
    "procrustes": "fidelity",
    "mirrored": "fidelity",
    "contact_auc": "fidelity",
    "contact_precision_at_k": "fidelity",
    "stems_nearer_frac": "fidelity",
    "ss_contacts_nearer": "fidelity",
    "angle_median": "regression-guard",
    "angle_max": "regression-guard",
    "joint_error_max": "regression-guard",
    "beads_missing": "consistency",
    "bead_colour_wrong": "consistency",
    "beads_false": "consistency",
    "beads_unmatched": "consistency",
    "overlaps": "readability",
    "text_off_page": "readability",
    "letter_overlaps": "readability",
    "fallbacks": "readability",
    "na_crossings": "readability",
    "bead_overlaps": "readability",
    "gap_label_collisions": "readability",
    "na_run_crossings": "readability",
    "na_protein_crossings": "readability",
    "seconds": "readability",  # the run time, judged with the page criteria (A2's kind)
}
# Å: a stem whose first arm is seen at least this far beside its second, across the drawn axis, has a side the
# drawing must keep. Below it the arms are seen nearly on top of each other and both sides fit the projection
# about equally. Equal to the layout's own cut, layout._SIDE_CLEAR (τ, set by the plan-07 Task 7p ruling), but
# kept separate so that raising that one shows here as mirrored stems.
SIDE_SEEN = 2.75


def _similarity_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD of 2D points b fitted onto a by translation, rotation (never reflection) and uniform scale, over a's
    spread. In closed form (the best turn is atan2 of the summed cross and dot products), elementwise sums only,
    so no BLAS kernel can change the result."""
    a0, b0 = a - a.mean(axis=0), b - b.mean(axis=0)
    dot = float((a0 * b0).sum())
    cross = float((b0[:, 0] * a0[:, 1] - b0[:, 1] * a0[:, 0]).sum())
    bb, aa = float((b0**2).sum()), float((a0**2).sum())
    residual = max(aa - (dot * dot + cross * cross) / max(bb, 1e-12), 0.0)
    return float(np.sqrt(residual / len(a)) / max(np.sqrt(aa / len(a)), 1e-9))


def _angle(u, v) -> float:
    """Degrees between two vectors."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    cos = float((u * v).sum()) / max(float(np.sqrt((u * u).sum() * (v * v).sum())), 1e-12)
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _turn(u, v) -> float:
    """Signed degrees from u to v, counterclockwise positive, in (-180, 180]."""
    return float(np.degrees(np.arctan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1])))


def _wrap(a: float) -> float:
    """An angle in degrees brought into [-180, 180)."""
    return (a + 180.0) % 360.0 - 180.0


def _norm(v) -> float:
    return float(np.sqrt((np.asarray(v, float) ** 2).sum()))


def measuring_frame(lay):
    """The projection shape() measures through. A plain projected Frame is the layout's own. An unrolled frame
    (a symmetric CylinderFrame, a protein-DNA RadialFrame) bends space so that it agrees with whatever it draws;
    it is replaced by the plain Frame with its orientation (page-right, page-up, toward the viewer) at the
    nucleic acid's centre, so the distortion of unrolling is measured instead of agreed with."""
    frame = lay.frame
    if frame is None or isinstance(frame, Frame) or lay.nucleic is None or not lay.nucleic.strands:
        return frame
    centre = np.vstack([s.p for s in lay.nucleic.strands]).mean(axis=0)
    u, v, w = frame.basis(centre)
    return Frame(centre, np.asarray(u, float), np.asarray(v, float), np.asarray(w, float))


def copy_frames(lay) -> dict[str, Frame]:
    """Nucleic-acid chain -> the plain Frame its copy is measured through, on an unrolled frame that unrolls the chain
    with the one protomer holding it (CylinderFrame.held, plan 28 F2): the unrolled frame's orientation at the centre
    of that copy's nucleic acid (the chains one protomer holds), so each copy is measured as seen from its own place
    round the axis, not from the assembly's centre (3EGZ: copies 2-4 were measured seen from 90, 180 and 270 degrees
    round). Empty for any other frame."""
    frame, na = lay.frame, lay.nucleic
    held = getattr(frame, "held", None) or {}
    if na is None or not held:
        return {}
    out = {}
    for k in sorted(set(held.values())):
        chains = {c for c, x in held.items() if x == k}
        pts = [s.p for s in na.strands if s.chain in chains]
        if not pts:
            continue
        centre = np.vstack(pts).mean(axis=0)
        u, v, w = frame.basis(centre)
        f = Frame(centre, np.asarray(u, float), np.asarray(v, float), np.asarray(w, float))
        out.update({c: f for c in chains})
    return out


def shape(lay) -> dict:
    """Drawn stems against the 3D structure seen through measuring_frame(lay): the figure's own plain frame, or
    for an unrolled frame the plain projection of the same orientation; on an unrolled frame that unrolls each
    nucleic-acid copy with the protomer holding it, each copy's stems through that copy's own plain frame
    (copy_frames; "copies_measured" their number), Spearman and Procrustes taken per copy and the worst copy's kept,
    and the angles judged against half the turn between neighbouring copies ("angle_limit", unrolled_angle_limit).
    Stems are drawn along their axes as seen
    in the layout's frame, so on a plain frame angle_median, angle_max and joint_error_max are 0 by construction:
    regression guards there, not fidelity checks; on an unrolled frame they measure the unrolling.
    stems_drawn; spearman (pairwise stem-centre distances, drawn against 3D, >= 5 stems); procrustes (>= 3
    stems); angle_median and angle_max (degrees, drawn direction against the seen 3D axis, stems seen side-on
    enough); mirrored (stems seen clearly, SIDE_SEEN, whose first
    arm is drawn on the side opposite the one it is seen on); when a duplex is split into pieces, bend_3d and
    bend_drawn (net signed turn from piece to piece, degrees, counterclockwise positive) and joint_error_max (the
    largest difference of the signed drawn and seen turns at one joint, so a bend drawn the wrong way counts).
    A drawn stem that is no piece of the structure's stems raises ValueError: no stem is left out unmeasured."""
    na, plain = lay.nucleic, measuring_frame(lay)
    if na is None or plain is None:
        return {}
    pieces = {p.id: p for p in stem_pieces(na)}
    unknown = [d.id for d in lay.stems if d.id not in pieces]
    if unknown:
        raise ValueError(f"drawn stems that are no piece of the structure's stems: {unknown}")
    stems = drawn_pieces(lay)  # as drawn: a curve's pieces where it draws their nucleotides (na_curve.drawn_pieces)
    if not stems:
        return {"stems_drawn": 0}
    copies = copy_frames(lay)

    def seen_by(d):  # the frame a drawn stem is measured through: its copy's, when its chains are one copy's
        p = pieces[d.id]
        got = {id(copies.get(na.strands[s].chain)) for s in (p.a, p.b)}
        return (
            copies[na.strands[p.a].chain] if copies.get(na.strands[p.a].chain) is not None and len(got) == 1 else plain
        )

    frames = [seen_by(d) for d in stems]
    axes = [stem_axis(pieces[d.id], na.strands) for d in stems]
    c3 = np.array([c for c, _ in axes])
    drawn = np.array([[d.cx, d.cy] for d in stems])
    out: dict = {"stems_drawn": len(stems)}
    groups: dict[int, list[int]] = {}  # the stems seen through one frame
    for i, f in enumerate(frames):
        groups.setdefault(id(f), []).append(i)
    rho, fit = [], []
    for idx in groups.values():
        if len(idx) >= 5:  # fewer stems give too few distances for a rank correlation
            pairs = list(combinations(idx, 2))
            d3 = [_norm(c3[i] - c3[j]) for i, j in pairs]
            d2 = [_norm(drawn[i] - drawn[j]) for i, j in pairs]
            rho.append(round(float(spearmanr(d3, d2).statistic), 3))
        if len(idx) >= 3:
            fit.append(round(_similarity_rmsd(frames[idx[0]].project(c3[idx]) * SCALE, drawn[idx]), 3))
    if rho:
        out["spearman"] = min(rho)
    if fit:
        out["procrustes"] = max(fit)
    if len(groups) > 1:
        out["copies_measured"] = len(groups)
        out["angle_limit"] = unrolled_angle_limit(lay.frame)
    errors, mirrored = [], 0
    for d, f, (centre, axis) in zip(stems, frames, axes):
        a2 = f.direction(centre, axis)
        if _norm(a2) > 0.5:  # seen side-on enough to have a direction on the page
            errors.append(_angle(a2, d.direction))
        (a, i0, i1), (b, j0, j1) = pieces[d.id].arm5, pieces[d.id].arm3
        across = na.strands[a].c1[i0 : i1 + 1].mean(axis=0) - na.strands[b].c1[j0 : j1 + 1].mean(axis=0)
        s2 = f.direction(centre, across)
        beside = float(d.direction[0] * s2[1] - d.direction[1] * s2[0])  # Å the first arm is seen left of the second
        if abs(beside) >= SIDE_SEEN:  # the arms are seen side by side across the drawn axis
            mirrored += (1 if beside > 0 else -1) != d.first
    out["angle_median"] = round(float(np.median(errors)), 1) if errors else None
    out["angle_max"] = round(float(np.max(errors)), 1) if errors else None
    out["mirrored"] = mirrored
    joints = []  # (seen turn, drawn turn), signed, between consecutive pieces of one duplex
    for k in range(len(stems) - 1):
        d, e = stems[k], stems[k + 1]
        if d.parent and d.parent == e.parent:
            (ca, a3), (cb, b3) = axes[k], axes[k + 1]
            seen = _turn(frames[k].direction(ca, a3), frames[k + 1].direction(cb, b3))
            joints.append((seen, _turn(d.direction, e.direction)))
    if joints:
        out["bend_3d"] = round(sum(j[0] for j in joints), 1)
        out["bend_drawn"] = round(sum(j[1] for j in joints), 1)
        out["joint_error_max"] = round(max(abs(_wrap(a - b)) for a, b in joints), 1)
    return out


def unrolled_angle_limit(frame) -> float:
    """The angle limit (degrees) of a figure whose copies are measured each through its own plain frame (copy_frames):
    half the turn between neighbouring copies round the axis. A copy's own parts lie up to that far round from its
    centre, and the unrolled frame shows each part from its own place round the axis, so its stems are drawn turned
    against the copy's single view by up to that much (3EGZ, four copies: measured 10.4 at most, limit 45). Not a
    fidelity check: angles are regression guards (ROLE); this keeps the guard from failing the unrolling itself."""
    return round(180.0 / max(int(getattr(frame, "n", 1)), 1), 1)


def _gap(a, b) -> float:
    """Distance between two boxes (x0, y0, x1, y1), 0 when they touch or overlap."""
    return float(np.hypot(max(a[0] - b[2], b[0] - a[2], 0), max(a[1] - b[3], b[1] - a[3], 0)))


def contacts(lay, sses) -> dict:
    """Are the protein elements touching a stem drawn nearer to it than those that do not touch it? Per stem
    with both kinds: the AUC of the box gaps (a touching element nearer than a non-touching one; ties count half)
    and the precision of the k nearest elements, k the number touching (elements tied at the k-th gap share the
    remaining places); both averaged over the stems. stems_nearer_frac: the fraction of the stems any element
    touches whose mean gap to those elements is below their median gap to every element (the spec's criterion),
    over the stems some drawn element does not touch: where every element touches a stem the touching elements
    would be compared with themselves (6VTX), so that stem is not scored but counted in stems_nearer_not_applicable
    (plan 28 F7, ruling Q5: not applicable, never a pass; stems_nearer_frac None when no stem is left).
    ss_contacts_nearer (reported, not judged): the same for each touched single-stranded nucleotide as drawn (a
    free end, hairpin loop, run or trace; 1URN's loop held by U1A), its gap from its drawn place to the element
    boxes; None before routing or when no such nucleotide is touched."""
    na = lay.nucleic
    if na is None:
        return {
            "contact_auc": None,
            "contact_precision_at_k": None,
            "stems_nearer_frac": None,
            "stems_nearer_not_applicable": 0,
            "ss_contacts_nearer": None,
        }
    drawn = [s for s in sses if s.id in lay.placed]
    aucs, precs, nearer, everyone = [], [], [], 0
    for d, body in zip(lay.stems, piece_bodies(lay)):  # each stem's body as drawn (a curve's piece: its stretches)
        mine = set(d.axial)
        touch = {s.id for s in drawn if any(mine & na.contacts.get(r, set()) for r in range(s.start, s.end + 1))}
        rest = [s.id for s in drawn if s.id not in touch]
        if not touch:
            continue
        if not rest:  # every drawn element touches it: nothing to be nearer than (F7)
            everyone += 1
            continue
        dist = {s.id: min(_gap(lay.placed[s.id].rect, b.rect) for b in body) for s in drawn}
        nearer.append(float(np.mean([dist[t] for t in sorted(touch)])) < float(np.median(list(dist.values()))))
        wins = sum((dist[t] < dist[r]) + 0.5 * (dist[t] == dist[r]) for t in sorted(touch) for r in rest)
        aucs.append(wins / (len(touch) * len(rest)))
        kth = sorted(dist.values())[len(touch) - 1]
        inside = [k for k, v in dist.items() if v < kth]
        tied = [k for k, v in dist.items() if v == kth]
        hits = sum(k in touch for k in inside) + (len(touch) - len(inside)) * sum(k in touch for k in tied) / len(tied)
        precs.append(hits / len(touch))
    return {
        "contact_auc": round(float(np.mean(aucs)), 3) if aucs else None,
        "contact_precision_at_k": round(float(np.mean(precs)), 3) if precs else None,
        "stems_nearer_frac": round(sum(nearer) / len(nearer), 3) if nearer else None,
        "stems_nearer_not_applicable": everyone,
        "ss_contacts_nearer": _ss_contacts_nearer(lay, drawn),
    }


def ss_drawn(lay) -> dict[tuple[int, int], tuple[float, float]]:
    """Where each single-stranded nucleotide is drawn: (strand, nucleotide) -> page point, from the routed paths
    (route_loops) and the traces."""
    out = {}
    for p in lay.na_paths:
        for k, xy in zip(p.nts, p.nt_xy):
            out[(p.strand, k)] = (float(xy[0]), float(xy[1]))
    for s, pts in lay.na_traces.items():
        for k, xy in enumerate(pts):
            out[(s, k)] = (float(xy[0]), float(xy[1]))
    return out


def _ss_contacts_nearer(lay, drawn) -> float | None:
    """The fraction of touched single-stranded nucleotides drawn nearer their touching elements (mean box gap)
    than the median element (contacts)."""
    na = lay.nucleic
    where = ss_drawn(lay)
    touchers: dict[tuple[int, int], set[str]] = {}
    for s in drawn:
        for r in range(s.start, s.end + 1):
            for x in na.contacts.get(r, ()):
                if x in where:
                    touchers.setdefault(x, set()).add(s.id)
    if not touchers or len(drawn) < 2:
        return None
    nearer = []
    for x, ids in sorted(touchers.items()):
        pt = (*where[x], *where[x])
        gap = {s.id: _gap(lay.placed[s.id].rect, pt) for s in drawn}
        nearer.append(float(np.mean([gap[i] for i in sorted(ids)])) < float(np.median(list(gap.values()))))
    return round(sum(nearer) / len(nearer), 3)


def _segments(pts) -> np.ndarray:
    """(n, 4) array of a polyline's segments: x0, y0, x1, y1."""
    a = np.asarray(pts, float).reshape(-1, 2)
    return np.hstack([a[:-1], a[1:]])


def _crossings(a: np.ndarray, b: np.ndarray) -> int:
    """Segment pairs of two polylines (as _segments) that cross properly (not just touch at an end). Products
    written out elementwise, never through BLAS."""
    if not len(a) or not len(b):
        return 0
    p, q = a[:, None, 0:2], a[:, None, 2:4]
    u, v = b[None, :, 0:2], b[None, :, 2:4]

    def turn(o, e, w):
        return (e[..., 0] - o[..., 0]) * (w[..., 1] - o[..., 1]) - (e[..., 1] - o[..., 1]) * (w[..., 0] - o[..., 0])

    return int(((turn(u, v, p) * turn(u, v, q) < 0) & (turn(p, q, u) * turn(p, q, v) < 0)).sum())


def _span(seg: np.ndarray) -> tuple[float, float, float, float]:
    """The box (x0, y0, x1, y1) around a polyline's segments."""
    xs, ys = seg[:, 0::2], seg[:, 1::2]
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def _boxes_meet(a, b) -> bool:
    return a[0] <= b[2] and b[0] <= a[2] and a[1] <= b[3] and b[1] <= a[3]


def _cross_count(xs: list, ys: list) -> int:
    """Pairs (x, y) of segment sets that cross at least once, after a box test."""
    spans = [_span(y) for y in ys]
    n = 0
    for x in xs:
        sx = _span(x)
        n += sum(1 for y, sy in zip(ys, spans) if _boxes_meet(sx, sy) and _crossings(x, y))
    return n


def _through_stem(seg: np.ndarray, d) -> bool:
    """Does any segment run through the stem's true, tilted body? A first or last segment that only grazes it, no
    deeper than a bead's radius (a joint leaving its face past a corner of its own piece), does not count."""
    last = len(seg) - 1
    return any(
        crosses_body(s[:2], s[2:], d) and (0 < k < last or graze(s[:2], s[2:], d) > _BEAD / 2)
        for k, s in enumerate(seg)
    )


def _box_hits(box, seg: np.ndarray) -> bool:
    """Does any segment run into the box (x0, y0, x1, y1)?"""
    return any(clip_len(s[:2], s[2:], box) > 1e-9 for s in seg)


def _box_meets_body(box, d) -> bool:
    """Does the box overlap the stem's (or element's) tilted body: an edge of one through the other, or the box
    inside it?"""
    c, s = np.cos(d.angle), np.sin(d.angle)
    hl, hw = d.length / 2, d.width / 2
    corners = [(d.cx + x * c - y * s, d.cy + x * s + y * c) for x, y in ((-hl, -hw), (hl, -hw), (hl, hw), (-hl, hw))]
    if _box_hits(box, _segments(corners + corners[:1])):
        return True
    centre = ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2)
    return crosses_body(centre, centre, d) or in_body(centre, d)


def _whole(d):
    """What a body counts as: itself; a stretch of a curve (na_curve.Chord), its stretch of the curve na_curve._FAR
    bands long (one thing however many of its stretches a line or a label meets there, as a piece is; two passes through
    two gyres are two)."""
    n = getattr(d, "curve", -1)
    if n < 0:
        return id(d)
    from .na_curve import _FAR

    return ("curve", n, int(((d.s0 + d.s1) / 2) // max(_FAR * d.band, 1e-9)))


def _trimmed(lay) -> list[np.ndarray]:
    """Each nucleic-acid path's segments as drawn among the elements (na_paths.among): a run or joint without the
    two known rail stretches it runs along to its stems' end faces, so a pass through its own stem counts; any other
    path without its first and last stretches inside its own stems."""
    at = stems_at(lay)  # a curve's nucleotides: its stretches as drawn
    return [_segments(among(p, at)) for p in lay.na_paths]


def na_crossings(lay, loops) -> dict:
    """Crossings of the drawn nucleic-acid strands, each strand as drawn among the elements (_trimmed: a run or joint
    from its stems' end faces): na_crossings, (run, joint or gap path, hairpin arc or free end) pairs crossing at
    least once plus (such a path, stem) pairs where it runs through the stem's body, its own stems included;
    na_run_crossings, pairs of runs, joints and gaps crossing at least once, except the two strands of one joint
    (they may cross where the duplex carries on, as along a helix); na_protein_crossings, (nucleic-acid path,
    protein loop) pairs crossing at least once."""
    drawn = list(zip(lay.na_paths, _trimmed(lay)))
    moving = [(p, seg) for p, seg in drawn if p.kind in ("run", "joint", "gap")]
    fixed = [seg for p, seg in drawn if p.kind in ("hairpin", "tail")]
    bodies = stem_bodies(lay)  # as drawn: a curve's stretches, the curve counted once (_whole)
    through = sum(
        len({_whole(d) for d in bodies if _boxes_meet(_span(seg), d.rect) and _through_stem(seg, d)})
        for _, seg in moving
        if len(seg)
    )
    at = {key: d.id for d in lay.stems for key in d.axial}

    def joins(p):  # the two pieces a joint joins; None (no exemption) for anything else or a joint not placed
        if p.kind != "joint" or len(p.between) != 2 or any(k not in at for k in p.between):
            return None
        return frozenset(at[k] for k in p.between)

    runs = sum(
        1
        for k, (p, a) in enumerate(moving)
        for q, b in moving[k + 1 :]
        if len(a)
        and len(b)
        and (joins(p) is None or joins(p) != joins(q))
        and _boxes_meet(_span(a), _span(b))
        and _crossings(a, b)
    )
    return {
        "na_crossings": _cross_count([seg for _, seg in moving], fixed) + through,
        "na_run_crossings": runs,
        "na_protein_crossings": _cross_count([seg for _, seg in drawn], [_segments(l.points) for l in loops]),
    }


def render_checks(lay, loops, fig, sses, *, model=None, bb=None) -> dict:
    """Rendering defects of a drawn figure: overlaps (element or stem boxes overlapping, except two consecutive
    pieces of one bent duplex, which meet at their joint: layout.joined_pieces, and a helix drawn over a sheet with
    that sheet's strands: layout.lying_over), text_off_page,
    letter_overlaps (base letters overlapping by more than a third of the smaller glyph box), beads_missing,
    bead_colour_wrong, beads_false and beads_unmatched (verify_contacts on the figure's beads: unbeaded contacts
    counted as defects, beads not in their toucher's chain colour, beads no atoms back, beads naming a residue or
    nucleotide the model lacks), fallbacks (loop and nucleic-acid routes with no clear way), na_crossings,
    na_run_crossings and na_protein_crossings (see na_crossings), bead_overlaps (pairs of contact beads whose discs
    overlap at their drawn size) and gap_label_collisions ("… N nt …" labels overlapping another one, or meeting a
    stem's or an element's body, a base letter, a contact bead, a nucleic-acid path other than their own (outside its
    own stems) or a protein loop: one per label and thing it meets).
    The bead checks against the coordinates need the io.read_model result and the Backbone the figure was drawn
    from; without them they are None (not measured)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    from .figcheck import element_overlaps

    overlaps = len(element_overlaps(lay))
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    page = fig.bbox
    ax = fig.axes[0]
    texts = [t for t in ax.texts if t.get_text().strip()]
    extents = [t.get_window_extent(renderer) for t in texts]
    off = sum(
        1 for e in extents if e.x0 < page.x0 - 1 or e.y0 < page.y0 - 1 or e.x1 > page.x1 + 1 or e.y1 > page.y1 + 1
    )
    letters = [e for t, e in zip(texts, extents) if (t.get_gid() or "").startswith("na:nt:")]  # pixels
    bad = 0
    ink = np.array([e.extents for e in letters], float).reshape(-1, 4)
    for i, a in enumerate(letters):  # only later letters whose boxes share area with this one's can count
        rest = ink[i + 1 :]
        near = (rest[:, 0] < a.x1) & (a.x0 < rest[:, 2]) & (rest[:, 1] < a.y1) & (a.y0 < rest[:, 3])
        for j in np.flatnonzero(near).tolist():
            b = letters[i + 1 + j]
            w = min(a.x1, b.x1) - max(a.x0, b.x0)
            h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if w > 0 and h > 0 and w * h > min(a.width * a.height, b.width * b.height) / 3:
                bad += 1
    # contact beads: centres and radii on the page, in pixels
    centres, radii = [], []
    for c in ax.collections:
        if (c.get_gid() or "").startswith("na:contact:"):
            xy = ax.transData.transform(c.get_offsets())
            centres += list(xy)
            radii += [float(np.sqrt(c.get_sizes()[0])) / 2 * fig.dpi / 72] * len(xy)
    bead_overlaps = 0
    at, size = np.asarray(centres, float).reshape(-1, 2), np.asarray(radii, float)
    for i in range(len(at) - 1):  # each pair once, as combinations gives them, by the same arithmetic
        d = at[i] - at[i + 1 :]
        bead_overlaps += int(np.count_nonzero(np.hypot(d[:, 0], d[:, 1]) < size[i] + size[i + 1 :] - 1e-6))
    # gap labels: (path index, box in page units, box in pixels)
    to_page = ax.transData.inverted()
    gaps = [
        (int(t.get_gid().rsplit(":", 1)[1]), tuple(to_page.transform_bbox(e).extents), e)
        for t, e in zip(texts, extents)
        if (t.get_gid() or "").startswith("na:gap:")
    ]
    paths = _trimmed(lay)
    protein = [_segments(l.points) for l in loops]
    bodies = [*stem_bodies(lay), *lay.placed.values(), *lay.ghosts]
    collisions = sum(1 for (_, a, _), (_, b, _) in combinations(gaps, 2) if _boxes_meet(a, b))
    for n, box, px in gaps:
        collisions += len({_whole(d) for d in bodies if _boxes_meet(box, d.rect) and _box_meets_body(box, d)})
        others = [seg for k, seg in enumerate(paths) if k != n] + protein
        collisions += sum(1 for seg in others if len(seg) and _boxes_meet(box, _span(seg)) and _box_hits(box, seg))
        collisions += sum(
            1 for e in letters if min(e.x1, px.x1) > max(e.x0, px.x0) and min(e.y1, px.y1) > max(e.y0, px.y0)
        )
        collisions += sum(
            1
            for c, r in zip(centres, radii)
            if np.hypot(c[0] - np.clip(c[0], px.x0, px.x1), c[1] - np.clip(c[1], px.y0, px.y1)) < r
        )
    beads: dict = dict.fromkeys(("beads_missing", "bead_colour_wrong", "beads_false", "beads_unmatched"))
    if model is not None and bb is not None:
        report = verify_contacts(model, getattr(fig, "foldmap_beads", None), lay, bb)
        beads = {
            "beads_missing": len(report.missing),
            "bead_colour_wrong": len(report.wrong_colour),
            "beads_false": len(report.false),
            "beads_unmatched": len(report.unmatched),
        }
    return {
        "overlaps": overlaps,
        "text_off_page": off,
        "letter_overlaps": bad,
        **beads,
        "fallbacks": sum(p.fallback for p in lay.na_paths) + sum(l.fallback for l in loops),
        **na_crossings(lay, loops),
        "bead_overlaps": bead_overlaps,
        "gap_label_collisions": collisions,
    }


def judge(m: dict) -> list[str]:
    """The metrics past their THRESHOLDS, as "name value (limit limit)"; unmeasured (None) ones pass, NaN fails. A
    figure whose copies are measured each in its own frame (shape's angle_limit, plan 28 F2) has its angles judged
    against that limit."""
    failed = []
    for key, limit in THRESHOLDS.items():
        v = m.get(key)
        if key in ("angle_median", "angle_max") and m.get("angle_limit") is not None:
            limit = m["angle_limit"]
        if limit is None or v is None:
            continue
        if v != v or ((v < limit) if key in _HIGH_IS_GOOD else (v > limit)):
            failed.append(f"{key} {v} (limit {limit})")
    return failed
