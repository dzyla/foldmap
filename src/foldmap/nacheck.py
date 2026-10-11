# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Checks of a nucleic-acid figure against the deposited coordinates.

verify_contacts: every contact bead a figure draws (a nucleotide touched by protein, in that chain's colour)
must be backed by real heavy atoms. The check shares no code with io.load_nucleic: it asks a gemmi
NeighborSearch over the model the figure was built from, atom by atom, and reads residues from the model itself.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import gemmi

from .io import _CONTACT as CONTACT

_NO_CELL = gemmi.UnitCell()  # NeighborSearch wants a cell; an empty one means none


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
    """One bead as drawn, and whether the toucher it names is within the cutoff of the nucleotide."""

    key: tuple[int, int]  # (strand, nucleotide) in the layout
    nucleotide: str  # "C:DG7"
    chain: str  # protein chain whose colour the bead carries
    colour: str
    residue: str  # the protein residue named as the toucher, "A:ARG42"
    pair: Pair | None  # closest atom pair of that residue and the nucleotide (None: no atoms in the model)
    verified: bool


@dataclass(frozen=True)
class Missed:
    """A drawn nucleotide a protein residue touches in the coordinates, with no bead."""

    nucleotide: str
    pair: Pair
    reason: str  # "missing" (a defect), or why none is expected: "focus", "no backbone", "first toucher"

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

    @property
    def false(self) -> list[Bead]:
        return [b for b in self.beads if not b.verified]

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


def _touch(ns, model, chain_name: str, res, held) -> dict[tuple, Pair]:
    """{(chain, seq, icode): closest Pair} of the protein residues with a heavy atom within CONTACT of a heavy
    atom of the nucleotide `res`. Protein: an amino acid by name, or a residue the figure's backbone holds."""
    found: dict[tuple, Pair] = {}
    nt = _tag(chain_name, res.name, res.seqid.num, res.seqid.icode)
    for atom in res:
        if atom.element.name == "H":
            continue
        for mark in ns.find_atoms(atom.pos, "\0", radius=CONTACT):
            cra = mark.to_cra(model)
            other = cra.residue
            key = (cra.chain.name, other.seqid.num, other.seqid.icode.strip())
            if (
                cra.atom.element.name == "H"
                or key == (chain_name, res.seqid.num, res.seqid.icode.strip())
                or not (key in held or _is_protein(other.name))
            ):
                continue
            d = atom.pos.dist(cra.atom.pos)
            if d <= CONTACT and (key not in found or d < found[key].dist):
                found[key] = Pair(_tag(*key[:1], other.name, *key[1:]), cra.atom.name, nt, atom.name, d)
    return found


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


def verify_contacts(model, layout, bb, *, path=None) -> ContactReport:
    """Check the contact beads of a drawn figure (draw() has filled layout.na_beads) against `model`, the
    gemmi model (or io.read_model result) the figure was built from; `bb` is the Backbone it used.
    A bead is verified when the protein residue it names, in the bead's chain, has a heavy atom within CONTACT
    of a heavy atom of its nucleotide. Every other protein residue touching a drawn nucleotide is listed in
    `missed`; reason "missing" is a defect, while "focus" (a chain outside the drawn subunit), "no backbone"
    (a residue the figure does not hold) and "first toucher" (the bead names an earlier residue) are expected.
    Touched nucleotides the figure does not draw (cropped away, or outside the focus) are listed under `cropped`.
    path: also read the file's deposited hydrogen bonds into the report."""
    if isinstance(model, tuple):
        model = model[0]
    report = ContactReport(hydrog=deposited_hydrog(path) if path is not None else [])
    na = layout.nucleic
    if na is None or not na.strands:
        return report
    ns = gemmi.NeighborSearch(model, _NO_CELL, 5).populate()
    drawn = {
        (lab.chain, lab.seq, lab.icode): (s, i) for s, st in enumerate(na.strands) for i, lab in enumerate(st.labels)
    }
    held = {(l.chain, l.seq, l.icode) for l in bb.labels}
    beads = {key: (chain, r, colour) for chain, items in layout.na_beads.items() for key, r, colour in items}
    focus = layout.focus_chains
    for chain in model:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_nucleic_acid():
                continue
            ident = (chain.name, res.seqid.num, res.seqid.icode.strip())
            near = _touch(ns, model, chain.name, res, held)
            nt = _tag(*ident[:1], res.name, *ident[1:])
            if ident not in drawn:
                report.cropped += [
                    Missed(nt, p, "focus" if focus and k[0] not in focus else "cropped")
                    for k, p in sorted(near.items())
                ]
                continue
            bead = beads.get(drawn[ident])
            named = None
            if bead is not None:
                lab = bb.labels[bead[1]]
                named = (lab.chain, lab.seq, lab.icode)
                pair = near.get(named)
                report.beads.append(
                    Bead(
                        drawn[ident],
                        nt,
                        bead[0],
                        bead[2],
                        _tag(lab.chain, lab.name, lab.seq, lab.icode),
                        pair,
                        pair is not None and lab.chain == bead[0],
                    )
                )
            for key, pair in sorted(near.items()):
                if key == named:
                    continue
                if key not in held:
                    reason = "no backbone"
                elif focus and key[0] not in focus:
                    reason = "focus"
                else:
                    reason = "first toucher" if bead is not None else "missing"
                report.missed.append(Missed(nt, pair, reason))
    return report
