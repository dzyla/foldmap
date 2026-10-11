# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Checks of a nucleic-acid figure against the deposited coordinates.

verify_contacts: every contact bead a figure draws (a nucleotide touched by protein, in that chain's colour)
must be backed by real heavy atoms.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .io import CONTACT


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
    reason: str  # "missing" (a defect), or why none is expected: "focus", "no backbone", "first toucher", "cropped"

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
    layout.focus_chains), "no backbone" (a residue the figure does not hold) and "first toucher" (the nucleotide's
    bead names an earlier residue) are expected. Touched nucleotides the figure does not draw are in `cropped`.
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
