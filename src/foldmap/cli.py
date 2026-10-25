# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Command line: `foldmap summary FILE` prints secondary structure and sheet ordering."""

from __future__ import annotations

import argparse
import dataclasses
import functools
import math
import re
import sys
from collections import Counter
from collections.abc import Callable
from fractions import Fraction
from pathlib import Path

import numpy as np

from . import layoutfile
from .criteria import LIMITS
from .dssp import assign_dssp
from .features import helices_on_sheets, sse_contacts
from .frame import (
    CylinderFrame,
    Frame,
    RadialFrame,
    RosetteFrame,
    dna_radial_frame,
    dna_view_frame,
    protomer_maps,
    rosette_frame,
    symmetric_frame,
    view_frame,
)
from .io import load_backbone, load_links, load_nucleic, read_model
from .layout import (
    RES_RANGE,
    SCALE,
    AmbiguousName,
    accessory_lines,
    build_layout,
    na_torn_by_unrolling,
    resolve,
    snap,
    template_rows,
)
from .membrane import SOLUBLE_BEYOND, crosses, crosses_joined, lies_flat_here
from .membrane import joins as membrane_joins
from .model import AUTO_NT, crop_nucleic, homolog_chains, nucleic_of_chains, nucleic_touching, res_ref
from .names import all_names, apply_names, check_spec
from .nucleic import StemView, stem_views, straight_duplex
from .palette import PALETTES
from .render import LOOP_STYLES, STYLES, draw, save
from .route import route_loops
from .science import science_lines, sheet_checks, sheet_judged
from .scope import drawn_scope, scope_note
from .sheets import build_sheets
from .spec import (
    FigureSpec,
    Problem,
    SavedDomains,
    check_references,
    data_look,
    draw_spec,
    figure_spec,
    layout_document,
    parse_domains,
    parse_protomers,
    parse_vector,
)
from .ss import build_sses
from .style import (
    BETWEEN_FORMS,
    CHOICES,
    INTERACTION_FORMS,
    LAYOUT_KEYS,
    SCALES,
    THEME_KEYS,
    THEMES,
    Style,
    describe,
    resolve_style,
)
from .symmetry import TOLERANCE, detect_symmetry

_vector, _protomers, _domain_spec = parse_vector, parse_protomers, parse_domains  # older names, used by the app


class UsageError(ValueError):
    """Options that do not go together (the command line exits 2, as for an unknown option)."""


class BadReference(UsageError):
    """An element reference that cannot be drawn (spec.check_references), said before layout: 'WHERE: TEXT' (the
    command line: 'foldmap: error: WHERE: TEXT', exit 2)."""

    def __init__(self, problem: Problem):
        super().__init__(str(problem))
        self.problem = problem


def _level_accessories(path: str) -> set[str]:
    """The accessory strands the figure `foldmap plot FILE` draws level beside their partner (Layout.accessory less
    Layout.unshown): `summary` words each accessory as that figure shows it (plan 12 fix round 4). Laid out only,
    not routed or drawn; none when it cannot be laid out."""
    spec = figure_spec(path)
    try:
        layout, _, _ = make_layout(path, spec.mode, spec.rotate, spec.flip_v, **spec.layout_kwargs())
    except ValueError:
        return set()
    away = {a for a, b in layout.unshown if layout.accessory.get(a) == b}
    return {a for a in layout.accessory if a not in away}


def summarize(path: str) -> str:
    model = read_model(path)  # once, for both loaders
    bb = load_backbone(path, model=model)
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, short_helices=True, bridges=dssp.bridges)
    sheets = build_sheets(sses, dssp.bridges)
    na = load_nucleic(path, bb, model=model)
    if len(bb) or not na.strands:
        lines = [f"{path}: {len(bb)} residues, {len(set(l.chain for l in bb.labels))} chain(s)"]
    else:  # a nucleic acid alone: its nucleotides, not 0 residues
        nts = sum(len(st) for st in na.strands)
        lines = [f"{path}: {nts} nucleotides, {len({l.chain for st in na.strands for l in st.labels})} chain(s)"]
    lines.append(
        f"helices: {sum(s.kind == 'H' for s in sses)}  strands: {sum(s.kind == 'E' for s in sses)}"
        f"  3₁₀ helices: {sum(s.kind == 'G' for s in sses)}"
    )
    for s in sses:
        if s.kind == "H":
            lines.append(f"  helix  {s.id}")
    level = _level_accessories(path) if any(sh.accessory for sh in sheets) else None
    for n, sh in enumerate(sheets, 1):
        flags = ("closed " if sh.closed else "") + ("ambiguous" if sh.ambiguous else "")
        lines.append(f"sheet {n}: {len(sh.strands)} strand(s) {flags}".rstrip())
        parts = []
        for k, (s, d) in enumerate(zip(sh.strands, sh.directions)):
            parts.append(f"{s.id}({'+' if d > 0 else '-'})")
            if k < len(sh.pair_kinds):
                parts.append({"P": "=par=", "A": "=anti=", "?": "=?="}[sh.pair_kinds[k]])
        lines.append("  " + " ".join(parts))
        lines += [f"  {line}" for line in accessory_lines(sh, level=level)]  # off the run, or not side by side
        lines += [f"  pairs with {b.id}, in another sheet: {a.id}" for a, b in sh.joined]
    if sheets:  # the science checks of the sheets as built (figcheck.science_checks reads them as drawn)
        lines.append("science checks (the sheets as built):")
        found = sheet_checks(sheets, sses, dssp.bridges)
        lines += [f"  {line}" for line in science_lines(found, sheet_judged(sheets))]
    try:
        sym = detect_symmetry(bb, sses)
    except ValueError:
        sym = None
    if sym is not None:
        extra = (
            f"; twist {np.degrees(abs(sym.twist)) * (-1 if sym.left_handed else 1):.1f}°, rise {abs(sym.rise):.1f} Å"
            + sym.hand
            if sym.kind == "H"
            else ""
        )
        lines.append(
            f"symmetry: {sym.label} — protomers {' | '.join(','.join(p) for p in sym.protomers)}"
            f"; CA RMSD {sym.rmsd:.2f} Å{extra}" + (f"; unpaired chains {','.join(sym.others)}" if sym.others else "")
        )
    if na.strands:
        kind = (
            "RNA" if all(st.rna for st in na.strands) else "DNA" if not any(st.rna for st in na.strands) else "DNA/RNA"
        )
        lines.append(
            f"{kind}: {len(na.strands)} strand(s) {', '.join(st.id for st in na.strands)}; "
            f"{len(na.pairs)} base pairs; {len(na.contacts)} protein residues make contacts"
        )
        for st in na.structure.stems:
            (a, i0, i1), (b, j0, j1) = st.arm5, st.arm3
            la, lb = na.strands[a].labels, na.strands[b].labels
            lines.append(
                f"  {st.name:6s} {len(st):3d} bp  {la[i0].chain}:{la[i0].tag}-{la[i1].tag} · "
                f"{lb[j0].chain}:{lb[j0].tag}-{lb[j1].tag}"
            )
    return "\n".join(lines)


CONSERVATION_LEGEND = "conservation (entropy-based, ConSurf colour scale)"
CONSURF_LEGEND = "conservation (ConSurf grades)"


def _any_centre(bb, nucleic) -> np.ndarray:
    if len(bb):
        return bb.ca.mean(axis=0)
    if nucleic is not None and nucleic.strands:
        return np.vstack([s.p for s in nucleic.strands]).mean(axis=0)
    return np.zeros(3)  # nothing to draw (a nucleic acid with --no-dna)


NA_HOLD_MIN = 3  # residues of one protomer touching a strand for that protomer to hold it


def _na_across_protomers(nucleic, sym, bb) -> bool:
    """Whether a nucleic acid with stems has a strand held by two or more protomers of `sym` (_na_shared): the
    precondition for drawing it in a projected frame (unrolling the assembly one protomer beside the next may tear
    that strand apart; layout.na_torn_by_unrolling decides)."""
    if nucleic is None or not nucleic.structure.stems:
        return False
    return _na_shared(nucleic, sym, bb)


def _na_shared(nucleic, sym, bb) -> bool:
    """Whether a strand of the nucleic acid is held by two or more protomers of `sym`, each touching it with at
    least NA_HOLD_MIN residues (nucleosome DNA, a duplex bound by both halves of a dimer)."""
    if nucleic is None or sym is None:
        return False
    sector = {c: k for k, chains in enumerate(sym.protomers) for c in chains}
    holders: dict[int, dict[int, set[int]]] = {}  # strand -> protomer -> residues touching it
    for r, touched in sorted(nucleic.contacts.items()):
        k = sector.get(bb.labels[r].chain)
        if k is not None:
            for s, _ in sorted(touched):
                holders.setdefault(s, {}).setdefault(k, set()).add(r)
    return any(sum(len(rs) >= NA_HOLD_MIN for rs in by.values()) >= 2 for by in holders.values())


DYAD_IN_PAGE = np.radians(10)  # a two-fold axis within this of the page plane is drawn as a mirror line


def _mirror_frame(frame, sym) -> None:
    """A dimer whose dyad lies in the page (perpendicular to the duplex of a RadialFrame, along page-up of a
    Frame) is drawn as one protomer and its mirror image: the frame learns the protomers (sector, classes, as
    symmetric_frame gives them) and the page x of the mirror line, where the dyad projects (snapped).

    The mirror is exact for the default view of a dimer of matching copies (1LMB, to 6e-16). Elsewhere it is best
    effort: a dyad up to DYAD_IN_PAGE off the page (as an --up or --view can leave it) is mirrored as if it lay
    in it, copies that do not match one-to-one are tied by vote, and a sheet shared between the two protomers is
    one item that no mirror line can pass through."""
    if isinstance(frame, RadialFrame):
        in_page = abs(float(sym.axis @ frame.axis)) < np.sin(DYAD_IN_PAGE)
    else:
        in_page = abs(float(sym.axis @ frame.v)) > np.cos(DYAD_IN_PAGE)
    if not in_page:
        return
    frame.sector, frame.classes = protomer_maps(sym)
    frame.mirror_x = float(snap(frame.project(sym.centre)[0, 0] * SCALE))


def _mirrors(two_fold: str, sym, bb) -> bool:
    """Whether Style.two_fold asks for the dyad view (_dyad_frame) of a two-fold assembly: mirror always, unrolled
    never, auto only where every protomer is one protein chain (TrmD, 1UAM; nucleic-acid chains do not count). With
    two or more chains a protomer (haemoglobin, 2HHB: α1β1 | α2β2), interfaces between different chains cross the
    dyad (α1β2, the molecule's main contact) and no mirror line keeps them together: the copies are unrolled."""
    if two_fold != "auto":
        return two_fold == "mirror"
    protein = {l.chain for l in bb.labels}
    return all(sum(c in protein for c in p) == 1 for p in sym.protomers)


ROSETTE_ORDERS = (4, 7)  # cyclic=rosette: Cn assemblies of these orders are drawn seen down their axis. Plan 21's
# look test kept these two: C3, C5, C6 and C8 rosettes (and a D2's two mirrored rings) read worse than their unrolled
# figures on most entries tried, with more check failures, so those stay unrolled with a note (_rosette_declined)


def _rosette_declined(ring, focus_chains=None, focus: str = "auto") -> str:
    """The note when cyclic=rosette is asked of a cyclic or dihedral assembly (`ring`, before focus) that is not drawn
    as a rosette for what it is or how much of it is drawn ('' otherwise; a two-fold has Style.two_fold): only C4 and
    C7 are offered, and a rosette needs the whole assembly, not one subunit (focus auto or protomer) or the chains
    named (`focus_chains`, from --focus A,B)."""
    if ring is None or ring.kind not in ("C", "D") or (ring.kind == "C" and ring.n == 2):
        return ""
    if ring.kind == "C" and ring.n in ROSETTE_ORDERS:
        if not focus_chains:
            return ""
        drawn = (
            "one subunit is drawn"
            if focus in ("auto", "protomer")
            else f"only chain(s) {', '.join(sorted(focus_chains))} are drawn"
        )
        return (
            f"cyclic=rosette: the {ring.label} assembly not drawn as a rosette, {drawn} and the rosette needs the "
            "whole assembly (--focus none)"
        )
    return (
        f"cyclic=rosette: the {ring.label} assembly not drawn as a rosette, the rosette is offered for four- and "
        "seven-fold assemblies only"
    )


def _no_rosette(sses, sym, nucleic, mode: str, mem, domains=None) -> str:
    """Why a Cn assembly (n in ROSETTE_ORDERS) asked for as a rosette (Style.cyclic) is drawn unrolled ('' when it is
    not): what the rosette does not draw yet. A rosette turns protomer 0's drawing for each copy, so everything
    drawn must belong to a protomer: no nucleic acid (its stems are not turned yet), no element of a chain outside
    every protomer (an unpaired chain), and each protomer must hold elements; a stack row and a membrane view are
    not seen down the axis; domains (`domains`, --domains or a layout's) are boxed on page axes, which do not follow
    a turned copy, so their panels would cut across its neighbours."""
    if mode != "projected":
        return "a rosette is drawn in projected mode only"
    if mem is not None:
        return "a membrane view keeps the bilayer across the page"
    if domains:
        return "domain panels are not yet drawn in a rosette (drop --domains to draw the rosette)"
    if nucleic is not None and nucleic.strands:
        return "its nucleic acid is not yet drawn in a rosette (draw it with --no-dna for the protein alone)"
    inside = {c for p in sym.protomers for c in p}
    if any(s.chain not in inside for s in sses):
        return "a chain outside every protomer is not yet drawn in a rosette"
    if any(not any(s.chain in p for s in sses) for p in sym.protomers):
        return "a protomer has no helix or strand"
    return ""


def _within_protomer(bridges, sym, bb) -> list:
    """The bridges between residues of one protomer: a rosette's sheets (each copy's own; a sheet two copies share
    is drawn as each copy's part, the pairings across copies left out)."""
    sector = {c: k for k, p in enumerate(sym.protomers) for c in p}
    return [b for b in bridges if sector.get(bb.labels[b.i].chain, -1) == sector.get(bb.labels[b.j].chain, -2)]


def _pairings_across(bridges, sym, bb, sses) -> list[tuple[str, str]]:
    """The strand pairs (element ids, in sequence order) that bridge two protomers: a rosette draws each copy's own
    sheets, so these pairings are not drawn (the notes name them)."""
    sector = {c: k for k, p in enumerate(sym.protomers) for c in p}
    of = {r: s.id for s in sses if s.kind == "E" for r in range(s.start, s.end + 1)}
    out = set()
    for b in bridges:
        ka, kb = sector.get(bb.labels[b.i].chain), sector.get(bb.labels[b.j].chain)
        if ka is not None and kb is not None and ka != kb and b.i in of and b.j in of:
            out.add((of[b.i], of[b.j]))
    order = {s.id: n for n, s in enumerate(sses)}
    return sorted(out, key=lambda ab: (order[ab[0]], order[ab[1]]))


def _dyad_frame(sses, sym) -> Frame | None:
    """A two-fold assembly seen from the side (plan 18): the dyad up the page, protomer 0 on the left and protomer 1
    to its right, so that _mirror_frame draws one as the other's mirror image. Page-right runs from protomer 0's
    element centroids toward protomer 1's, perpendicular to the dyad; the view is their cross product (a rotation,
    never a mirror), turned about the dyad as _dyad_turn chooses (plan 19), so that elements close in 3D stay close.
    Both axes are snapped before they choose the page axes (the dyad and the centroids differ in the last bit between
    BLAS kernels). None when a protomer has no element, or the two sit on the dyad."""
    if sym.kind != "C" or sym.n != 2 or len(sym.protomers) != 2:
        return None
    cents = [[s.centroid for s in sses if s.chain in p] for p in sym.protomers]
    if not all(cents):
        return None
    axis = sym.axis / np.linalg.norm(sym.axis)
    right = np.mean(cents[1], axis=0) - np.mean(cents[0], axis=0)
    right = right - (right @ axis) * axis
    if np.linalg.norm(right) < 1.0:  # Å: the two halves on the dyad, no side to see them from
        return None
    up, right = snap(axis), snap(right / np.linalg.norm(right))
    face = np.cross(right, up)  # the view at turn 0: protomer 1 straight to protomer 0's right
    turn = _dyad_turn(np.array(cents[0] + cents[1]), up, right, face)
    if turn:  # turned about the dyad: protomer 1 still to the right, a little behind or in front of protomer 0
        t = np.radians(turn)  # (cos and sin may differ in the last bit between libms: snapped)
        face = snap(np.cos(t) * face + np.sin(t) * right)
    return view_frame(sses, up=up, view=face)


DYAD_TURNS = (0, 15, -15, 30, -30, 45, -45, 60, -60, 75, -75)  # degrees about the dyad, in order of preference
DYAD_TURN_GAIN = 0.01  # turns this close to the best rank correlation count as tied: the earliest is taken


def _dyad_turn(cents, up, right, face) -> int:
    """The turn about the dyad (degrees, one of DYAD_TURNS) a two-fold assembly is seen at (plan 19): the one whose
    projection keeps the order of the elements' distances best, by the rank correlation (Spearman) of their centroids'
    distances on the page and in 3D, over every pair. Every view perpendicular to the dyad draws an exact mirror image,
    but seen square-on to the line between the protomers (turn 0) an interface lying along the dyad can fall on the
    mirror line, where the layout pulls its overlapping elements apart (TrmD, 1UAM: rank correlation 0.40 projected,
    0.35 laid out; at 75°, 0.96 and 0.84, and cross-chain contacts drawn at their 3D distance instead of 4.6 times it).
    The projection predicts the laid-out figure (plan 19 Task 1's table), so one layout is enough. Turns stay within
    ±75°: protomer 1's centre stays to the right of protomer 0's (at ±90° one is behind the other). Ties, and turns
    within DYAD_TURN_GAIN of the best, go to the earliest in DYAD_TURNS (the smallest turn). Distances are snapped and
    tied ranks averaged, so the choice is the same on any CPU."""
    from scipy.spatial.distance import pdist
    from scipy.stats import rankdata

    def ranks(d):
        r = rankdata(snap(d))
        return r - r.mean()

    if len(cents) < 3:
        return DYAD_TURNS[0]
    rel = cents - cents.mean(axis=0)
    in3d = ranks(pdist(rel))
    y = (rel * up).sum(axis=1)
    score = {}
    for turn in DYAD_TURNS:
        t = np.radians(turn)
        x = (rel * (np.cos(t) * right - np.sin(t) * face)).sum(axis=1)  # page-right at this turn
        page = ranks(pdist(snap(np.stack([x, y], axis=1))))
        norm = np.sqrt((page * page).sum() * (in3d * in3d).sum())
        score[turn] = round(float((page * in3d).sum() / norm), 9) if norm > 0 else 0.0
    best = max(score.values())
    return next(t for t in DYAD_TURNS if score[t] >= best - DYAD_TURN_GAIN)


def _side_by_side(frame, sym, nucleic, bb) -> None:
    """Copies of an assembly on a straight duplex frame (no dyad in the page: _mirror_frame) drawn alike, side by
    side, as symmetric_frame draws an assembly: the frame learns the protomers (sector, classes) and no mirror line.
    Not when a strand is held by two or more protomers (_na_shared): unrolled side by side, the copies would pull
    the duplex they share apart (random-100 held-out A15: 5K7Z's two dimers, the duplex on A and B alone, were
    drawn unalike with no tie)."""
    if not _na_shared(nucleic, sym, bb):
        frame.sector, frame.classes = protomer_maps(sym)


def _unalike(lay, bb) -> float:
    """The largest deviation between copies drawn (figcheck.copies_alike, A15's measure; inf where too few match), 0
    with none judged: how unalike a layout draws its copies. Letting a split copy block go (layout._one_per_mate)
    draws 2BBJ's and 5K7Z's copies alike, but on 9GTP's whole file it moved untied sheets further from their copies'
    (M~Q 25.8 -> 52.3, c~a 1.1 -> 9.9): make_layout lays out both ways and keeps the more alike."""
    from types import SimpleNamespace

    from .figcheck import copies_alike

    return max(copies_alike(SimpleNamespace(layout=lay, bb=bb)).values(), default=0.0)


def _na_extent(na, extent: str):
    """(the nucleic acid to draw, whether `auto` cropped a large one to the stems the protein touches)."""
    if extent == "contacts":
        return crop_nucleic(na), False
    if extent == "auto" and na.size > AUTO_NT and na.contacts:
        return crop_nucleic(na, by="stems"), True
    return na, False


def _left_out(values, links, drawn: set[str], data, crosslinks) -> list:
    """Notes (spec.Problem warnings) for the values and crosslinks read on chains the figure leaves out (--chains,
    --focus): read, not drawn."""
    out = []
    if values is not None:
        missed: dict[str, int] = {}
        for chain, _ in values.values:
            if chain not in drawn:
                missed[chain] = missed.get(chain, 0) + 1
        for chain, n in missed.items():
            text = f"chain {chain}'s values ({n} residue{'s' if n > 1 else ''}) are not drawn: the figure leaves "
            out.append(Problem(f"--data {data}", text + f"chain {chain} out (--chains, --focus)", False))
    if links:
        away = [x for x in links if x[0][0] not in drawn or x[1][0] not in drawn]
        if away:
            chains = sorted({c for (a, _), (b, _), _ in away for c in (a, b) if c not in drawn})
            n = len(away)
            text = f"{n} crosslink{'s' if n > 1 else ''} reach{'' if n > 1 else 'es'} a chain the figure leaves out"
            out.append(Problem(f"--crosslinks {crosslinks}", f"{text} ({', '.join(chains)}; --chains, --focus)", False))
    return out


def _read_user_data(data, column, crosslinks, bb, na_chains):
    """userdata.read_user_data for make_layout: the warnings (each row not used) printed to standard error, then the
    first fatal problem raised as a BadReference (exit 2, before layout), else one line saying what was read from
    each file."""
    from .userdata import read_user_data, summary

    values, links, problems = read_user_data(data, column, crosslinks, bb, na_chains=na_chains)
    for problem in problems:  # userdata says the first SAID_ROWS rows of a file, then how many more
        if not problem.fatal:
            print(f"foldmap: warning: {problem}", file=sys.stderr)
    for problem in problems:
        if problem.fatal:
            raise BadReference(problem)
    for line in summary(values, links, data, crosslinks):
        print(f"foldmap: {line}", file=sys.stderr)
    return values, links, problems


def make_layout(
    path,
    mode: str = "projected",
    rotate: float = 0.0,
    flip_v: bool = False,
    dna: bool = True,
    look: Style | None = None,
    symmetry: str = "auto",
    protomers: list[list[str]] | None = None,
    symmetry_tol: float = TOLERANCE,
    assembly: str = "auto",
    rename: dict | None = None,
    swap: list | None = None,
    move: list | None = None,
    up=None,
    view=None,
    domains=None,
    msa=None,
    msa_reference: str | None = None,
    membrane: str = "off",
    chains: list[str] | None = None,
    focus: str = "auto",
    uniprot: str | None = None,
    model=None,
    names: list[str] | None = None,
    like: dict | None = None,
    like_out: dict | None = None,
    data: str | None = None,
    data_column: str | None = None,
    crosslinks: str | None = None,
):
    """`names` (--names): conventional names, each `deposited`, `beta` or a YAML file (names.apply_names), applied
    over the generic labels in order; `rename` wins over all of them. `like` (--like, plan 13a): drawn like a
    reference structure, {structure, chain, match, layout} (_drawn_like): one protein chain of this structure, the
    one `chains` names (two named is a UsageError), else its first (ruling 7: the scope line says it was chosen,
    Layout.like_chain). A frame that cannot be carried (_uncarried: a stack, a straight duplex) draws it on its own,
    the note saying why. `like_out` (a dict, foldmap compare, foldmap check's H1-H5): given what the reference gave,
    as far as it got: reference ((layout, sses, bb, model, chain) of the reference laid out, or None), alignment
    (align.ChainAlignment or None), mates (align.Mates or None), like (the like.Like this layout took, or None:
    drawn on its own) and refused (True when drawn on its own as not the reference's fold, the note saying so: H1).
    `data`, `data_column`, `crosslinks` (--data, --data-column, --crosslinks; plan 13d): the researcher's tables,
    read on every protein chain (userdata) after the references are checked: a file that cannot be read is a
    BadReference before layout, each row not used a warning on standard error (the first SAID_ROWS of a file, then
    how many more); values and crosslinks on chains the figure leaves out are noted; what was read is kept on the
    layout (Layout.residue_data, crosslinks, data_problems: their defaults, none, when no file is given)."""
    for spec in names or []:
        check_spec(spec)
    model = model or read_model(path, assembly)  # once, for every loader below (`model`: already read, same assembly)
    bb = load_backbone(path, assembly, model=model)
    read_bb = bb  # every protein chain, before --chains: what the references are checked against
    present = list(dict.fromkeys(l.chain for l in bb.labels))
    asked = [c for c in chains or [] if c in present]  # the protein chains asked for (Layout.chains_asked)
    chain_b, chosen = None, False  # the chain drawn like the reference; chosen by --like, not asked for
    if like:  # one protein chain drawn like one of the reference's (plan 13a ruling 6); nucleic-acid chains as asked
        chain_b = one_chain(chains, present, "--chains", str(path))
        chosen = not asked and len(present) > 1  # several, none named: the first, said as chosen (not asked for)
        chains = [chain_b, *(c for c in chains or [] if c not in present)]
    if chains and any(c in present for c in chains):  # protein chains named: only those (else every one)
        bb = bb.only([c for c in chains if c in present])
    if (look or Style()).color_by == "conservation" and not msa:
        raise ValueError("conservation colouring needs an alignment: add --msa ALIGNMENT")
    mem = None
    flip = membrane == "flip"
    if membrane not in ("off", None):
        from .membrane import find_membrane

        mem = find_membrane(path, read_bb)  # on every protein chain: --chains narrows the drawing, not the bilayer
        if mem is not None and flip:  # the other face is the outside: a 180° turn, never a mirror
            mem = dataclasses.replace(mem, normal=-mem.normal, centre=-mem.centre)
        if mem is None:
            print(f"foldmap: warning: no membrane found in {path}; drawing it without one", file=sys.stderr)
    dssp = assign_dssp(bb)
    matched: dict = {}  # helix pieces cut or joined so copies match (ss._copy_breaks): told by scope.drawn_scope
    in_view = mem is not None and mode == "projected"  # helix runs one or two residues apart joined (ruling T5-1)
    sses = build_sses(
        bb,
        dssp.ss,
        short_helices=(look or Style()).helices_310,
        bridges=dssp.bridges,
        matched=matched,
        join=functools.partial(membrane_joins, bb, mem) if in_view else None,
        whole_if=functools.partial(crosses_joined, bb, mem) if in_view else None,
    )
    sheets = build_sheets(sses, dssp.bridges)
    nucleic = load_nucleic(path, bb, assembly, model=model) if dna else None
    na_chains = list(dict.fromkeys(st.chain for st in nucleic.strands)) if nucleic is not None else []
    if chains:
        held = list(dict.fromkeys(st.chain for st in nucleic.strands)) if nucleic is not None else []
        missing = [c for c in chains if c not in present and c not in held]
        if missing:
            have = f"protein chains {', '.join(present) or 'none'}" + (
                f", nucleic-acid chains {', '.join(held)}" if held else ""
            )
            raise ValueError(f"no chain {', '.join(missing)}; the structure has {have}")
        if any(c in held for c in chains):  # nucleic-acid chains named: only those (else every strand)
            nucleic = nucleic_of_chains(nucleic, chains)
    sym = detect_symmetry(bb, sses, symmetry, protomers, symmetry_tol) if sses else None
    ring = sym  # the assembly as read, before focus may leave one subunit (_rosette_declined)
    shared = _homolog_labels(bb, sses)  # copies of one protein share labels, however much of each is modelled
    copies = homolog_chains(bb)
    if focus == "auto" and (look or Style()).highlight != "none":
        focus = "none"  # highlighting the asymmetric unit or a protomer asks for the whole assembly
    focus_chains = _focus_chains(focus, sym, bb, sses, nucleic)
    strands = nucleic.strands if nucleic is not None else []
    read = list(dict.fromkeys([l.chain for l in bb.labels] + [s.chain for s in strands]))
    subunits = None
    if focus_chains and sym is not None:  # what scope.drawn_scope says was left out
        subunits = (_subunits(sym, focus_chains), "helical" if sym.kind == "H" else sym.label)
    if focus_chains and nucleic is not None and nucleic.strands:  # one subunit: the nucleic acid it holds
        nucleic = nucleic_touching(nucleic, focus_chains)
    na_cropped, na_read = False, nucleic.size if nucleic is not None else 0
    if nucleic is not None and nucleic.strands:  # after focus: auto crops to the stems the subunit touches
        nucleic, na_cropped = _na_extent(nucleic, (look or Style()).na_extent)
    partial, context, unjoined = set(), None, set()
    whole, primes, relabel = sses, {}, {}  # every element read; labels shared by copies and primes of lent elements
    if focus_chains:  # a large assembly: one subunit in full, plus what its neighbours lend to its fold
        from .layout import provisional

        labels = {k: p.label for k, p in provisional(sses).placed.items()}
        labels.update(shared)
        rank = {s.id: k for k, s in enumerate(sses)}
        sses, partial = _focus_subset(sses, focus_chains, dssp.bridges)
        unjoined = {(a.id, b.id) for a, b in zip(sses, sses[1:]) if a.chain == b.chain and rank[b.id] != rank[a.id] + 1}
        sheets = build_sheets(sses, dssp.bridges)
        order = sorted(partial, key=lambda c: min(s.start for s in sses if s.chain == c))
        primes = {c: "′″‴"[min(n, 2)] for n, c in enumerate(order)}
        for s in sses:
            if s.chain in partial:
                relabel[res_ref(s.chain, s.first, s.last)] = labels[s.id] + primes[s.chain]
        context = _context(sym, focus_chains, bb, bool(partial))
        sym = None  # one subunit: no copies to draw alike
    if shared:  # applied to the elements actually drawn (focus mode may have left some out)
        for s in sses:
            if s.id in shared and s.chain not in partial:
                relabel[res_ref(s.chain, s.first, s.last)] = shared[s.id]
    refs = FigureSpec(  # the element references asked for, checked before anything is laid out (Q-A19)
        Path(path),
        look=look or Style(),
        chains=chains,
        focus=focus,
        domains=domains,
        rename=dict(rename or {}),
        swap=list(swap or []),
        move=list(move or []),
    )
    kept = set(focus_chains) | set(partial) if focus_chains else None  # the chains focus draws, lent ones included
    for problem in check_references(refs, read_bb, drawn=kept):
        if problem.fatal:  # the warnings were said where the spec was read (_spec_of, the app, the MCP server)
            raise BadReference(problem)
    user = _read_user_data(data, data_column, crosslinks, read_bb, na_chains) if data or crosslinks else None
    up_from_membrane = False  # `up` set by the membrane (its normal), not asked for: the cyclic note says so
    if mem is not None:
        if up is None and view is None:  # the bilayer horizontal, outside on top (unless the symmetry axis
            along = sym is not None and abs(float(np.dot(sym.axis, mem.normal))) > 0.9  # already stands upright)
            if not along:
                up, up_from_membrane = mem.normal, True
            elif (mem.sided and float(np.dot(sym.axis, mem.normal)) < 0) or (flip and not mem.sided):  # turn it over
                sym = dataclasses.replace(
                    sym, axis=-sym.axis, angles=[-a for a in sym.angles], rise=-sym.rise, twist=-sym.twist
                )
    duplex = nucleic is not None and straight_duplex(nucleic.structure, nucleic.strands)
    dyad = None  # a two-fold assembly seen from the side, drawn as one copy and its mirror image (plan 18)
    offered = False  # a C4 or C7 assembly drawn whole that a rosette can draw: Style.cyclic chooses how
    rosette_asked = (  # Style.cyclic rosette: a cyclic assembly of ROSETTE_ORDERS copies seen down its axis
        (look or Style()).cyclic == "rosette" and sym is not None and sym.kind == "C" and sym.n in ROSETTE_ORDERS
    )
    cyclic_note = _rosette_declined(ring, focus_chains, focus) if (look or Style()).cyclic == "rosette" else ""
    if rosette_asked and (duplex or rotate or flip_v or up is not None or view is not None):
        if duplex:
            why = "it sits on a straight duplex"
        elif up_from_membrane and not rotate and not flip_v:  # the membrane set `up`; nothing else was asked
            why = "a membrane view keeps the bilayer across the page"
        else:
            why = "a view, a turn or a flip was asked for"
        cyclic_note = f"cyclic=rosette: the {sym.label} assembly not drawn as a rosette, {why}"
    if duplex:
        frame = (  # protein on a straight duplex: the duplex across the page, the protein below it
            dna_radial_frame(sses, nucleic)
            if not rotate and not flip_v
            else dna_view_frame(sses, nucleic, rotate, flip_v)
        )
    elif (  # an assembly: every protomer drawn alike, side by side (unless that tears a nucleic acid apart)
        sym is not None
        and not rotate
        and not flip_v
        and up is None
        and view is None
        and not (_na_across_protomers(nucleic, sym, bb) and na_torn_by_unrolling(nucleic, sses, sym))
    ):
        # two_fold=unrolled (or auto with a protomer of two or more protein chains): a two-fold assembly's copies
        # side by side, as any assembly's (never the dyad view)
        dyad = _dyad_frame(sses, sym) if _mirrors((look or Style()).two_fold, sym, bb) else None
        frame = dyad or symmetric_frame(sses, sym)
        if sym.kind == "C" and sym.n in ROSETTE_ORDERS:  # what a rosette cannot draw yet ('' when it can)
            why = _no_rosette(sses, sym, nucleic, mode, mem, domains)
            offered = not why  # the app shows the choice
        if rosette_asked:  # cyclic=rosette (plan 21): seen down the axis, each copy protomer 0 turned, unless
            if why:
                cyclic_note = f"cyclic=rosette: the {sym.label} assembly drawn unrolled, {why}"
            else:
                frame = rosette_frame(sses, sym)
                sheets = build_sheets(sses, _within_protomer(dssp.bridges, sym, bb))
                cyclic_note = (
                    f"{sym.label} assembly seen down its {sym.n}-fold axis (cyclic=rosette): each copy drawn as the "
                    f"first turned by a multiple of {360 // sym.n if 360 % sym.n == 0 else round(360 / sym.n, 1)}°"
                )
                across = _pairings_across(dssp.bridges, sym, bb, sses)
                if across:
                    cyclic_note += (
                        f"; {len(across)} strand pairing(s) between copies not drawn, each copy's sheets its own: "
                        + ", ".join(f"{a}~{b}" for a, b in across)
                    )
    else:  # stem axes differ in the last bit between BLAS kernels: snapped before they choose the page axes
        if rosette_asked and not cyclic_note:
            cyclic_note = (
                f"cyclic=rosette: the {sym.label} assembly not drawn as a rosette, its nucleic acid is held by two or "
                "more copies"
            )
        stems = stem_views(nucleic.structure, nucleic.strands) if nucleic is not None else []
        stems = [StemView(snap(v.centroid), snap(v.axis), v.n) for v in stems]
        frame = view_frame(sses + stems, rotate, flip_v, up, view)
        if up_from_membrane and mode == "projected" and not rotate and not flip_v:  # turned about the normal (T5-3)
            frame = _in_sequence(frame, sses, bb, mem, copies)
    if sym is not None and sym.kind == "C" and sym.n == 2 and isinstance(frame, (RadialFrame, Frame)):
        _mirror_frame(frame, sym)
    if duplex and sym is not None and not getattr(frame, "sector", None):
        _side_by_side(frame, sym, nucleic, bb)
    uncarried = _uncarried(frame, dyad, mode, duplex, sym) if like else None
    links = load_links(path, bb, assembly, model=model)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, links).items():
        contacts[key] = contacts.get(key, 0) + n
    ann = None
    if uniprot not in (None, "off", "none"):
        from .uniprot import annotate

        ann = annotate(path, bb, None if uniprot == "auto" else uniprot)
    source, found_note = _domain_source(domains), ""
    if domains == "cath":  # CATH's published chopping (PDBe), or a note that there is none
        from .cath import domains as cath_domains

        domains, found_note = cath_domains(path, bb, sses)
    if domains == "uniprot":
        if ann is None:
            raise ValueError("--domains uniprot needs UniProt annotation: add --uniprot auto (or an accession)")
        domains = _uniprot_domains(ann, sses)
    rows = (  # βα rows are rigid units: automatic domains never split one
        template_rows(sheets, sses, copies) if mode == "projected" and (look or Style()).fold_template == "auto" else []
    )
    conventions = functools.partial(all_names, path, whole=whole)  # dep: and beta: references (layout.resolve)
    tm = _transmembrane(sses, bb, mem) if mem is not None and mode == "projected" else set()
    in_band = _band_ids(sses, bb, mem, tm) if tm else frozenset()  # auto domains never split these (R17)
    named = _domains(domains, sses, contacts, rows, copies, partial, bb, conventions, band=in_band)
    absorbed = (  # the band joined every domain contacts found (5JDN's two halves): the note says so (review-3 M1)
        bool(in_band)
        and source == "inferred"
        and not named
        and bool(_domains(domains, sses, contacts, rows, copies, partial, bb, conventions))
    )
    named_note = domain_note(source, named, found_note, band=bool(in_band), absorbed=absorbed)
    for _, ids in named:  # a domain holds together: its elements pull on each other
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                key = frozenset((ids[a], ids[b]))
                contacts[key] = contacts.get(key, 0) + _DOMAIN_PULL
    over = (  # helices packed on a sheet face, drawn over it; not in an assembly drawn alike by its symmetry, nor
        # beside nucleic acid, whose stems and free ends are placed by the protein elements they touch (plan 11 T14)
        helices_on_sheets(bb, sses, sheets)
        if (look or Style()).helix_on_sheet == "auto" and sym is None and (nucleic is None or not nucleic.strands)
        else {}
    )

    springs = _bridge_springs(bb, sses, links)
    made = None  # drawn like a reference: what this layout takes from it (_drawn_like)
    if like:
        made = _drawn_like(
            like, bb, chain_b, sses, sheets, look or Style(), mode, rotate, flip_v, up, view, dna, assembly, uncarried
        )
        if made.like is not None:  # the reference's frame carried over: the page is the reference's
            frame, dyad = made.frame, None
        if like_out is not None:
            like_out.update(
                reference=made.reference,
                alignment=made.alignment,
                mates=made.mates,
                like=made.like,
                refused=made.refused,
            )

    def lay_out_in(frame):  # build_layout through a frame (the dyad view's, or the unrolled one's), as a partial
        return functools.partial(
            build_layout,
            sses,
            sheets,
            frame,
            mode,
            contacts=contacts,
            nucleic=nucleic,
            style=look,
            bridges=springs,
            partial=partial,
            domains=named,
            rename=rename,
            swap=swap,
            move=move,
            keep_depth=mem is not None and mode == "stack",
            na_cropped=na_cropped,
            template=(look or Style()).fold_template,
            copies=copies,
            relabel=relabel,
            names=functools.partial(
                apply_names, list(names), path, whole=whole, chains=set(present), copies=copies, primes=primes
            )
            if names
            else None,
            over=over,
            conventions=conventions,
            alike_first=mem is None,  # a membrane view is judged by its bands (C17), not its copies (A15)
            **in_membrane(frame),
            **({"like": made.like} if made is not None and made.like is not None else {}),
        )

    def faces_of(frame):  # the membrane's two faces as projected through `frame` (None: no band to keep to)
        return _faces_on_page(frame, sses, bb, mem, tm) if tm else None

    def in_membrane(frame) -> dict:
        """A projected membrane view (plan 26 Task 2): build_layout's membrane arguments through `frame`
        (_membrane_args, the faces as projected through it)."""
        if mem is None or mode != "projected":
            return {}
        return _membrane_args(sses, bb, mem, tm, faces_of(frame), named, in_band)

    def finish(lay, frame, faces=None):  # what make_layout sets on a layout after build_layout through `frame`
        lay.helix_matched = {"cut": matched.get("cut", 0), "joined": matched.get("joined", 0)}
        lay.helix_matched_ids = list(matched.get("ids", []))
        lay.membrane_joined_ids = list(matched.get("membrane", []))
        lay.domain_source, lay.domain_note = (source if lay.domains else ""), named_note
        if lay.domains:  # copies of one domain share a tone
            from .features import domain_tones

            lay.domain_tones = domain_tones(lay.domains, sses, copies)
        lay.focus = _focus((look or Style()).highlight, bb, sym)
        if focus_chains:
            lay.focus, lay.focus_chains, lay.context = set(focus_chains), set(focus_chains), context
            lay.unjoined, lay.subunits = unjoined, subunits
        lay.chains_read = read
        lay.chains_asked = asked
        lay.like_chain = chain_b if chosen else ""
        lay.uncarried = _uncarried(frame, dyad if frame is dyad else None, mode, duplex, sym)
        lay.na_read = na_read if na_cropped else 0
        lay.links = links
        lay.bridges = list(dssp.bridges)
        lay.res_chain = [l.chain for l in bb.labels]
        lay.res_name = [l.name for l in bb.labels]
        lay.res_labels = list(bb.labels)
        lay.res_linked = sorted(bb.linked)
        lay.res_b = [float(x) for x in bb.b] if bb.b is not None else []
        with_elements = {s.chain for s in sses}
        lay.coil_chains = [(c, n) for c, n in Counter(lay.res_chain).items() if c not in with_elements]
        lay.uniprot = ann
        lay.cyclic_note = cyclic_note
        lay.symmetry, lay.bilayer = ring, mem  # what foldmap data reports (datafile); not drawn
        if made is not None:
            lay.like_colour = dict(made.like.colour) if made.like is not None else {}
            lay.like_extras = [e for e in made.extras if e in lay.placed]
            lay.like_note = made.note(lay)
        first = [s.centroid for s in sses if s.chain == sses[0].chain] if sses else [_any_centre(bb, nucleic)]
        centre = np.mean(first, axis=0)
        lay.page_axes = [list(map(float, a)) for a in frame.basis(centre)]
        lay.page_centre = list(map(float, bb.ca.mean(axis=0) if len(bb) else centre))
        if mem is not None:
            lay.membrane = _membrane_band(lay, sses, bb, mem, faces)
        return lay

    def laid_out(frame, mirror_only=False):
        """The layout through `frame`, as make_layout keeps it; None with `mirror_only` where build_layout drops the
        frame's mirror line (a βα row across it), before laying anything out."""
        lay_out = lay_out_in(frame)
        lay = lay_out(mirror_only=True) if mirror_only else lay_out()
        if lay is None:
            return None
        faces = faces_of(frame)  # a projected membrane view: the band is the faces as projected (C17, plan 26)
        lay = finish(lay, frame, faces)
        if lay.let_go:  # copy items let go (layout._one_per_mate): kept only when the copies come out more alike (A15)
            held = finish(lay_out(let_go=False), frame, faces)
            if _unalike(held, bb) <= _unalike(lay, bb):
                lay, lay_out = held, functools.partial(lay_out, let_go=False)
        return lay

    def off_band(lay) -> int:  # transmembrane helices drawn off their band (A18; 0 off a membrane view)
        if mem is None or mode != "projected" or lay.membrane is None:
            return 0
        from .figcheck import membrane_view

        return len(membrane_view(lay, sses, bb, mem)["tm_short"])

    lay = laid_out(frame, mirror_only=dyad is not None)
    if lay is None:  # no mirror image to draw in the dyad view (a βα row across the line): the copies unrolled
        lay = laid_out(symmetric_frame(sses, sym))
    elif dyad is not None:  # the dyad view is for the mirror image of copies drawn alike
        mirrored = any(m for _, m in lay.copy_ties.values())
        gap, off = _unalike(lay, bb), off_band(lay)
        if not mirrored or gap > LIMITS["copies_alike"] or off:
            # where no copy is tied as a mirror image, the copies are unrolled side by side as for any assembly;
            # where the mirror image leaves a transmembrane helix off its band
            # (A18) or draws the copies unalike (A15), they are unrolled when that leaves fewer helices off the
            # band, or as many and the copies more alike (bands first: alikeness never bought with a helix off). A
            # mirror image A15 passes is never traded for a smaller deviation: the choice stays off float noise.
            unrolled = laid_out(symmetric_frame(sses, sym))
            off_unrolled = off_band(unrolled)
            alike = gap <= LIMITS["copies_alike"]
            if (
                not mirrored
                or off_unrolled < off
                or (off_unrolled == off and not alike and _unalike(unrolled, bb) < gap)
            ):
                lay = unrolled
    lay.two_fold_assembly = (
        sym is not None and sym.kind == "C" and sym.n == 2 and len(sym.protomers) == 2 and not duplex
    )
    lay.cyclic_assembly = offered
    if user is not None:  # the researcher's data (plan 13d): Layout's defaults (none) when not given
        from .userdata import crosslink_lengths, crosslink_problems

        values, links, problems = user
        drawn = {l.chain for l in bb.labels}
        if focus_chains:
            drawn &= set(focus_chains) | set(partial)
        notes = _left_out(values, links, drawn, data, crosslinks)
        lengths = crosslink_lengths(links or [], read_bb)  # Cα–Cα on the structure: the span check, the drawing
        notes += crosslink_problems(links or [], lengths, crosslinks, drawn)
        for note in notes:
            print(f"foldmap: warning: {note}", file=sys.stderr)
        lay.residue_data, lay.crosslinks, lay.data_problems = values, list(links or []), [*problems, *notes]
        lay.crosslink_ca = lengths
    if msa:
        from .sequence import grades, is_consurf, read_alignment, read_consurf, residue_conservation

        if is_consurf(msa):  # ConSurf's own grades by residue, where the structure has the table's residue there
            rows: dict[str, dict[tuple[int, str], tuple[str, int]]] = {}
            for (c, n, icode), row in read_consurf(msa).items():
                rows.setdefault(c, {})[(n, icode)] = row
            names: dict[str, dict[tuple[int, str], str]] = {}
            for l in bb.labels:
                names.setdefault(l.chain, {})[(l.seq, l.icode)] = l.name

            def agreement(c: str, chain: str) -> tuple[int, int]:  # numbers in both, of them with the same residue
                shared = rows[c].keys() & names[chain].keys()
                return len(shared), sum(rows[c][k][0] == names[chain][k] for k in shared)

            def matches(c: str, chain: str, borrowed: bool) -> bool:
                """>= 90 % of table chain c's rows placed on `chain` name its residue; a chain borrowing another
                chain's grades must also share >= 20 residues or half the shorter of the two."""
                shared, same = agreement(c, chain)
                enough = not borrowed or shared >= min(20, 0.5 * min(len(rows[c]), len(names[chain])))
                return shared > 0 and same >= 0.9 * shared and enough

            def source(chain: str) -> str | None:
                if chain in rows:
                    if matches(chain, chain, False):
                        return chain
                    shared, same = agreement(chain, chain)
                    print(
                        f"foldmap: warning: {msa} does not match chain {chain} ({same} of {shared} residue names "
                        "agree); chain left without grades",
                        file=sys.stderr,
                    )
                    return None
                # assembly copies (A1, A-1) and chains identical to the one a ConSurf-DB table covers
                stripped = re.sub(r"-?\d+$", "", chain)
                candidates = ([stripped] if stripped in rows else []) + sorted(rows)
                return next((c for c in candidates if matches(c, chain, True)), None)

            sources = {c: source(c) for c in names}

            def grade(l) -> int:
                row = rows[sources[l.chain]].get((l.seq, l.icode)) if sources[l.chain] else None
                return row[1] if row and row[0] == l.name else 0

            lay.res_grade = [grade(l) for l in bb.labels]
            if not any(lay.res_grade):
                print(f"foldmap: warning: no residue of {path} matches {msa}; nothing is graded", file=sys.stderr)
            lay.cons_label = CONSURF_LEGEND
        else:
            scores = residue_conservation(bb, read_alignment(msa), msa_reference)
            lay.res_cons = [scores.get(k, float("nan")) for k in range(len(bb))]
            lay.res_grade = grades(lay.res_cons)
            lay.cons_label = CONSERVATION_LEGEND
    return lay, sses, bb


def one_chain(asked, present: list[str], flag: str, structure: str) -> str:
    """The protein chain of a structure drawn like another, or like which another is drawn (--like): the one of
    `present` (its protein chains) that `asked` names, else its first; two named is a UsageError naming `flag`."""
    named = [c for c in asked or [] if c in present]
    if len(named) > 1:
        raise UsageError(
            f"--like draws one protein chain of each structure; {flag} names {len(named)} of {structure}'s "
            f"({', '.join(named)}): name one"
        )
    if named:
        return named[0]
    if not present:
        raise ValueError(f"{structure} has no protein chain to draw like another")
    return present[0]


def _drawn_like(
    like: dict, bb, chain_b, sses, sheets, look: Style, mode, rotate, flip_v, up, view, dna, assembly, uncarried=None
):
    """What a figure drawn like a reference (make_layout's `like`: {structure, chain, match, layout}) takes from it:
    the reference laid out first, one protein chain (the one like['chain'], --like-chain, or its layout file names,
    else its first: the note says so), with the options and edits of
    like['layout'] (--like-layout, a saved layout of it) when given, else this figure's view and look; its chain
    aligned to this structure's `chain_b` (like['match']: auto, sequence, structure or a two-row alignment file;
    align.align_chains), the elements mated (align.element_mates) and like.like_from's Like built from the
    reference's finished layout, with its frame carried over (frame.frame_like). A namespace: like (None: drawn on
    its own), frame, extras (this structure's elements the reference lacks) and note(layout) (the notes' line,
    naming the extras as drawn), and what compare reads (make_layout's like_out): reference, alignment and mates,
    each None where it was not reached, and refused (not the reference's fold). Drawn on its own, the note says
    why: this figure's frame cannot be carried (`uncarried`, _uncarried: an assembly drawn by its symmetry, a stack,
    a straight duplex), nor the reference's, or it is not the reference's fold (Mates.same_fold; but not with
    like['same_protein']: two states of one protein by A19's rule, foldmap diff, whose domains may move apart past
    the fold's TM-score)."""
    from types import SimpleNamespace

    from .align import align_chains, element_mates
    from .frame import frame_like
    from .like import like_from, like_note, reference_name, refused_note, uncarried_note

    given = str(like["structure"])
    name = like.get("name") or reference_name(given)  # the app names an upload as uploaded

    def alone(note: str, reference=None, alignment=None, mates=None, refused=False) -> SimpleNamespace:
        return SimpleNamespace(
            like=None,
            frame=None,
            extras=(),
            note=lambda lay: note,
            reference=reference,
            alignment=alignment,
            mates=mates,
            refused=refused,
        )

    ref = resolve_structure(given)
    if like.get("layout"):
        spec = figure_spec(ref, layout_file=like["layout"])
        args, kw = (spec.mode, spec.rotate, spec.flip_v), spec.layout_kwargs()
    else:  # this figure's look, less its marks: they name this structure's elements, so the reference's references
        # check (Q-A19) leaves them out; only render reads a mark, so the reference is laid out as before
        mine = (  # and less its own colours and widths (elstyle): they name this structure's parts too
            dataclasses.replace(look, mark="", element_colors="", element_widths="")
            if look is not None and (look.mark or look.element_colors or look.element_widths)
            else look
        )
        args, kw = (mode, rotate, flip_v), {"dna": dna, "look": mine, "assembly": assembly, "up": up, "view": view}
    model = like.get("model")  # already read as `assembly` reads it (compare's diff page), else read here
    if model is None or like.get("layout"):
        model = read_model(ref, kw["assembly"])
    present = list(dict.fromkeys(l.chain for l in load_backbone(ref, kw["assembly"], model=model).labels))
    asked = [c.strip() for c in str(like.get("chain") or "").split(",") if c.strip()]
    missing = [c for c in asked if c not in present]
    if missing:
        raise ValueError(f"no protein chain {', '.join(missing)} in {given}; it has {', '.join(present) or 'none'}")
    named = asked or kw.get("chains")
    chain_a = one_chain(named, present, "--like-chain", given)
    first = not any(c in present for c in named or []) and len(present) > 1  # chosen, not named: the note says so
    if uncarried is not None:  # this figure's own frame: drawn on its own
        return alone(uncarried_note(name, uncarried, mine=True))
    kw.update(chains=[chain_a], like=None, data=None, data_column=None, crosslinks=None)  # the data are this figure's
    lay_a, sses_a, bb_a = make_layout(ref, *args, model=model, **kw)
    reference = (lay_a, sses_a, bb_a, model, chain_a)
    if lay_a.uncarried is not None:  # the reference's frame (a stack, a straight duplex)
        return alone(uncarried_note(name, lay_a.uncarried, mine=False), reference)
    match = like.get("match") or "auto"
    aln = align_chains(bb_a, chain_a, bb, chain_b, match)
    mates = element_mates(sses_a, sses, aln)
    if not mates.same_fold and not like.get("same_protein"):  # two states of one protein (foldmap diff): alike
        return alone(refused_note(name, aln, mates), reference, aln, mates, refused=True)
    try:
        frame = frame_like(lay_a.frame, aln.rotation, aln.translation)
    except TypeError:  # _uncarried names every frame frame_like refuses: not reached
        return alone(uncarried_note(name, ("copies", type(lay_a.frame).__name__), mine=False), reference, aln, mates)
    made = like_from(lay_a, sses_a, sses, mates, look, sheets)

    def note(lay) -> str:
        extras = [lay.placed[e].label or e for e in mates.extras if e in lay.placed]
        crossed = [
            (lay.placed[s.id].label, made.crossed[s.id]) for s in sses if s.id in made.crossed and s.id in lay.placed
        ]
        chain = f"{chain_a} (its first protein chain)" if first else chain_a
        return like_note(name, chain, aln, mates, extras, match, crossed)

    return SimpleNamespace(
        like=made,
        frame=frame,
        extras=tuple(mates.extras),
        note=note,
        reference=reference,
        alignment=aln,
        mates=mates,
        refused=False,
    )


def _uncarried(frame, dyad, mode: str, duplex: bool, sym) -> tuple[str, str] | None:
    """Why a figure's page cannot be carried to another structure drawn like it (--like; frame.frame_like carries a
    plain Frame only), as (kind, what): ('copies', ...) an assembly drawn by its symmetry (a rosette, unrolled, the
    dyad view's mirror image), ('stack', ...) --mode stack, ('duplex', ...) protein on a straight duplex; None for
    a plain frame (a membrane view, --up, --view, --rotate and --flip-v are part of one)."""
    label = sym.label if sym is not None else "two-fold"
    if isinstance(frame, RosetteFrame):
        return "copies", f"a {label} rosette"
    if isinstance(frame, CylinderFrame):
        return "copies", f"an unrolled {label} assembly"
    if dyad is not None or getattr(frame, "mirror_x", None) is not None or getattr(frame, "sector", None):
        return "copies", f"a {label} mirror image"
    if mode != "projected":
        return "stack", "--mode stack"
    if duplex or isinstance(frame, RadialFrame):
        return "duplex", "a straight duplex"
    return None


def _membrane_band(lay, sses, bb, mem, faces=None) -> dict | None:
    """The band on the page: where the transmembrane elements pass the two membrane faces (medians). A
    depth-keeping stack wrapped into rows has one band per row ("bands"), each from its own row's elements; a
    row with no transmembrane element (a soluble partner chain) has none. `faces` (a projected view whose
    transmembrane elements kept their height, _faces_on_page): the faces themselves, as projected (C17)."""
    rows = [set(ids) for ids in lay.stack_rows] or [set(lay.placed)]
    if faces is not None:
        bands = [faces]
    else:
        bands = [_face_medians(lay, [s for s in sses if s.id in ids], bb, mem) for ids in rows]
    bands = [b for b in bands if b is not None]
    if not bands:
        return None
    y_out, y_in = bands[0]
    return {
        "y": (min(y_out, y_in), max(y_out, y_in)),
        "bands": [(min(o, i), max(o, i)) for o, i in bands],
        "inside": "below" if y_in < y_out else "above",
        "source": mem.source,
        "sides": mem.sided,
    }


def _membrane_places(sses, bb, mem, tm: set[str], faces) -> tuple[set[str], dict[str, int], dict[str, int]]:
    """Where a projected membrane view draws each element that crosses no membrane (plan 26 Task 2 rulings), from its
    CA run's depths: `flat`, the helices (and 3-10 helices) that lie along the membrane (membrane.lies_flat_here:
    tilted FLAT_TILT or more from the normal, FACE_FLAT_TILT in a face's interfacial zone, ruling T5-2); `on_face`,
    those of them whose CA centre lies in the slab or within SOLUBLE_BEYOND beyond a face (amphipathic and
    interfacial helices: 32KS's NS1, A:200-219), each drawn lying on the face its centre is nearer (+1 the face higher
    on the page, -1 the lower: faces = (outer, inner) page y), and (0) every other element whose CA centre lies in
    the slab, kept at its projected height (32KS A:328-339); and `sides`, every element whose CA centre lies more than
    SOLUBLE_BEYOND beyond a face (soluble, as A18 counts it: ruling R4), and every other one not lying flat whose CA
    run lies wholly beyond a face (ruling T5-5: 6ZKT's M:50-51, 18-19 Å out, was drawn 10 Å into the slab), kept off
    the band on that face's side (+1 above, -1 below). With no band (faces None) only `flat` is read. Elements of
    `tm`, and any other crossing the membrane (a subunit in a barrel lumen), are none of these."""
    depth = mem.depth(bb.ca)
    up = 1 if faces is None or faces[0] >= faces[1] else -1  # the page side of the outside
    flat, on_face, sides = set(), {}, {}
    for s in sses:
        seg = depth[s.start : s.end + 1]
        if s.id in tm or crosses(seg, mem.half):
            continue
        centre = float(seg.mean())
        side = up if centre >= 0 else -up
        lying = s.kind in ("H", "G") and lies_flat_here(s.axis, seg, mem)
        if lying:
            flat.add(s.id)
        if faces is None:
            continue
        beyond = float(seg.min()) >= mem.half or float(seg.max()) <= -mem.half  # wholly beyond a face
        if abs(centre) > mem.half + SOLUBLE_BEYOND or (beyond and not lying):
            sides[s.id] = side
        elif lying:
            on_face[s.id] = side
        elif abs(centre) <= mem.half and s.kind != "G":  # a 3-10 turn in the slab: placed as any other
            on_face[s.id] = 0
    return flat, on_face, sides


def _band_ids(sses, bb, mem, tm: set[str]) -> set[str]:
    """What a projected membrane view's band holds (plan 26 ruling R17; review-3 M3, M4): the transmembrane elements
    `tm`, every other element crossing the membrane (membrane.crosses: 9FWZ's FimA subunit in the usher's lumen, which
    does not define it) and the elements _membrane_places sets on a face or holds in the slab. Contract: the ids
    _membrane_places puts in `on_face` do not depend on the faces' page order, only their signs do, so any faces
    serve; both orders are read and the contract asserted. Auto domains never split these (features.auto_domains
    `band`); a domain holding any is not soluble (_domain_sides)."""
    depth = mem.depth(bb.ca)
    crossing = {s.id for s in sses if crosses(depth[s.start : s.end + 1], mem.half)}
    up, down = (set(_membrane_places(sses, bb, mem, tm, faces)[1]) for faces in ((1.0, -1.0), (-1.0, 1.0)))
    assert up == down, "_membrane_places' face and slab ids depend on the faces' order"
    return set(tm) | crossing | up


def _membrane_args(sses, bb, mem, tm: set[str], faces, named, in_band: set[str]) -> dict:
    """build_layout's membrane arguments for a projected view through a frame whose faces project to `faces` ((outer,
    inner) page y; None: no band): every transmembrane element keeps its projected height across the band (`level`,
    `band`), standing upright and reaching the face its run comes nearer (`nearer`, _nearer_faces, ruling T5-4); each
    other helix stands or lies by its tilt (membrane.lies_flat_here, rulings R7 and T5-2); a flat one at a face or in
    the slab lies on that face, in rows (`on_face`, ruling R8), and a soluble element stays off the band on its own
    side (`sides`), as _membrane_places reads them. With domains (`named`; plan 26 Task 3), each
    domain holding none of the band's elements (`in_band`, _band_ids) is soluble: it and every element of it stay off
    the band on the side its centroid lies (`domain_sides`, _domain_sides; ruling T3-c)."""
    flat, on_face, sides = _membrane_places(sses, bb, mem, tm, faces)
    if faces is None:  # no band: no height to keep, only the helices turned
        return {"upright": set(), "flat": flat}
    out = {
        "upright": tm,
        "flat": flat,
        "level": tm,
        "band": tuple(sorted(faces)),
        "on_face": on_face,
        "sides": sides,
        "nearer": _nearer_faces(sses, bb, mem, tm, faces),
    }
    if named:  # each soluble domain on the side its centroid lies, every element of it with it (ruling T3-c)
        out["domain_sides"] = _domain_sides(named, sses, bb, mem, set(in_band), faces)
        for name, side in out["domain_sides"].items():
            sides.update({k: side for k in dict(named)[name]})
    return out


def _nearer_faces(sses, bb, mem, tm: set[str], faces) -> dict[str, int]:
    """Ruling T5-4: for each transmembrane element, the face its CA run comes nearer in 3D, +1 the face higher on the
    page, -1 the lower (faces = (outer, inner) page y), as build_layout's `nearer` takes it: one drawn wholly between
    the faces is moved to that face (layout._to_face). None for a run as near both (on the grid)."""
    depth = mem.depth(bb.ca)
    up = 1 if faces[0] >= faces[1] else -1  # the page side of the outside
    out = {}
    for s in sses:
        if s.id not in tm:
            continue
        seg = depth[s.start : s.end + 1]
        gap_out, gap_in = float(snap(mem.half - seg.max())), float(snap(seg.min() + mem.half))
        if gap_out != gap_in:
            out[s.id] = up if gap_out < gap_in else -up
    return out


def _domain_sides(named, sses, bb, mem, held: set[str], faces) -> dict[str, int]:
    """The soluble domains of a projected membrane view (plan 26 ruling T3-c), each by the side of the band its
    centroid lies on: +1 above, -1 below, as _membrane_places' `sides` (faces = (outer, inner) page y). A domain is
    soluble when it holds no element of `held` (the band's, _band_ids: a domain holding any keeps its place in the
    band)."""
    depth = mem.depth(bb.ca)
    up = 1 if faces[0] >= faces[1] else -1  # the page side of the outside
    by_id = {s.id: s for s in sses}
    out = {}
    for name, ids in named:
        mine = [by_id[k] for k in ids if k in by_id]
        if not mine or held.intersection(ids):
            continue
        at = np.concatenate([np.arange(s.start, s.end + 1) for s in mine])
        out[name] = up if float(depth[at].mean()) >= 0 else -up
    return out


def _transmembrane(sses, bb, mem) -> set[str]:
    """The elements that cross the membrane (membrane.crosses: their CA run crosses its mid-plane and spans half its
    thickness, as _face_medians and A18 count them) and define it (Membrane.defining): in a projected membrane view
    they keep their projected height (build_layout `level`), so each crosses the band where the membrane is
    (random-100 C17)."""
    depth = mem.depth(bb.ca)
    out = set()
    for s in sses:
        if not crosses(depth[s.start : s.end + 1], mem.half):
            continue
        if mem.defining is not None and mem.defining.isdisjoint(range(s.start, s.end + 1)):
            continue
        out.add(s.id)
    return out


AZIMUTH_STEP = 5  # degrees: the turns about the membrane normal a membrane view's page-right is chosen from (T5-3)
# A turn is taken only when it raises the worst copy's tau by 1/MIN_TAU_GAIN_INV (0.1) or more (rulings T5-3, R23:
# today's turn is kept where the gain is under 0.1). Scored by the tau summed over the chains, every turn worth less on
# the membrane campaign grew its page past the campaign's +15% gate: 6JO6 (tau +0.014, crossings +51 % against
# d132d4c), 6ZKT with domains (+0.092, +19 %) and 30BU with domains (+0.066, +32 %).
MIN_TAU_GAIN_INV = 10


def _sequence_score(frame, chains: list[list[np.ndarray]]) -> int:
    """Kendall's concordant less discordant pairs of each chain's band helices, their page x (snapped) against their
    order in sequence, summed over the chains (`chains`: each chain's helix centroids in sequence order)."""
    score = 0
    for cents in chains:
        x = [float(snap(c @ frame.u)) for c in cents]
        score += sum(int(np.sign(x[j] - x[i])) for i in range(len(x)) for j in range(i + 1, len(x)))
    return score


def _worst_tau(frame, proteins: list[list[list[np.ndarray]]]) -> Fraction:
    """Ruling R23: the band's tau as a turn is judged by it, exactly: each chain's Kendall tau (_sequence_score over
    its pairs), the least over the copies of each protein, then the least over the proteins (`proteins`: each
    protein's copies, each copy its band helix centroids in sequence order, two or more)."""
    return min(
        min(Fraction(_sequence_score(frame, [cents]), len(cents) * (len(cents) - 1) // 2) for cents in copies)
        for copies in proteins
    )


def _in_sequence(frame, sses, bb, mem, copies: dict[str, tuple[str, int]] | None = None):
    """Ruling T5-3: a membrane view's page-right turned about the membrane normal (page-up), so that the band helices
    (the transmembrane helices, _transmembrane) stand left to right as nearly in sequence order as a turn can put them:
    of the turns in AZIMUTH_STEP steps, the one whose worst copy's Kendall tau of page x against sequence is highest
    (ruling R23, _worst_tau: each protein's least tau over its copies, `copies` as model.homolog_chains gives them,
    each chain its own protein without; then the least over the proteins); on a tie the frame as view_frame chose it,
    else the least turn (anticlockwise seen from the outside first); the frame as chosen when the best gains less than
    0.1 (MIN_TAU_GAIN_INV), or when it lowers the tau pooled over every band helix of the view (ruling R27: 6ZKT, a
    photosystem of many chains, turned 105° for one chain's tau, -1.0 to -0.6, its pooled tau fell). Still a
    projection: the same frame turned, never a reordering. 32KS chain A's TM1 stood in the middle of its band, tau
    0.47; turned 35°, 0.78. In a two-fold dimer whose axis is the normal a turn that orders one copy reverses the
    other: scored by the pair's summed tau, 32KS turned 55° and drew its protomers through each other; by its worse
    copy it turns 15°, both copies at 0.156 or better."""
    tm = _transmembrane(sses, bb, mem)
    chains: dict[str, list[np.ndarray]] = {}
    for s in sses:
        if s.id in tm and s.kind == "H":
            chains.setdefault(s.chain, []).append(np.asarray(s.centroid, float))
    groups: dict[str, list[list[np.ndarray]]] = {}
    for c, cents in chains.items():
        if len(cents) > 1:
            groups.setdefault((copies or {}).get(c, (c, 0))[0], []).append(cents)
    if not groups:
        return frame
    proteins = list(groups.values())
    u0, v = np.asarray(frame.u, float), np.asarray(frame.v, float)
    side = np.cross(v, u0)  # u0 turned a quarter turn about v
    turns = sorted(range(-180 + AZIMUTH_STEP, 181, AZIMUTH_STEP), key=lambda a: (abs(a), a < 0))
    best, best_frame = None, frame
    for a in turns:
        t = np.radians(a)
        u = u0 if a == 0 else snap(np.cos(t) * u0 + np.sin(t) * side)
        turned = frame if a == 0 else Frame(frame.origin, u, frame.v, np.cross(u, v))
        score = _worst_tau(turned, proteins)
        if best is None or score > best:
            best, best_frame = score, turned
    if (best - _worst_tau(frame, proteins)) * MIN_TAU_GAIN_INV < 1:  # tau gains less than 1/10
        return frame
    every = [cents for copies in proteins for cents in copies]
    if _sequence_score(best_frame, every) < _sequence_score(frame, every):  # the pooled tau falls (R27)
        return frame
    return best_frame


def _faces_on_page(frame, sses, bb, mem, tm: set[str]) -> tuple[float, float]:
    """Page y (layout units, snapped) of the membrane's outer and inner face, projected through the frame at the
    transmembrane elements' centre (build_layout `band`, _membrane_band `faces`)."""
    first = next(s for s in sses if s.id in tm)
    at = np.mean([bb.ca[s.start : s.end + 1].mean(axis=0) for s in sses if s.id in tm], axis=0)
    mid = at - mem.normal * float(np.ravel(mem.depth(at))[0])
    y_out, y_in = (
        float(frame.project(mid + mem.normal * h, [first.chain])[0, 1]) * SCALE for h in (mem.half, -mem.half)
    )
    return float(snap(y_out)), float(snap(y_in))


def _face_medians(lay, sses, bb, mem) -> tuple[float, float] | None:
    """Where these elements pass the outer and the inner membrane face on the page (medians), if any crosses the
    membrane (membrane.crosses) and defines it (Membrane.defining: a beta-barrel's strands, not a subunit lying in its
    lumen). An element passes a face where its CA run does; one whose run stops short of a face (32KS: TM helices of
    15-19 residues in a 30 Å slab; 6MJ2's S1-S4, 1-2 Å short of the inner face) passes it where the run's straight
    line through its depths would, carried on along the element as drawn."""
    depth = mem.depth(bb.ca)
    outs, ins = [], []
    for s in sses:
        p = lay.placed.get(s.id)
        seg = depth[s.start : s.end + 1]
        if p is None or not crosses(seg, mem.half):
            continue
        if mem.defining is not None and mem.defining.isdisjoint(range(s.start, s.end + 1)):
            continue
        along = np.linspace(0.0, 1.0, len(seg))  # the line through its first and last turn (4 residues each)
        t0, t1, d0, d1 = along[:4].mean(), along[-4:].mean(), seg[:4].mean(), seg[-4:].mean()
        slope = (d1 - d0) / (t1 - t0) if t1 > t0 else 0.0
        at0 = d0 - slope * t0

        def page_y(face: float, seg=seg, p=p, slope=slope, at0=at0) -> float | None:
            k = np.flatnonzero(np.diff(np.sign(seg - face)) != 0)
            if k.size:
                k = int(k[0])
                t = (k + (face - seg[k]) / (seg[k + 1] - seg[k])) / (len(seg) - 1)
            elif abs(slope) > 1e-9:
                t = (face - at0) / slope  # beyond the run's end: its line carried on
            else:
                return None
            return float(p.n_port[1] + t * (p.c_port[1] - p.n_port[1]))

        o, i = page_y(mem.half), page_y(-mem.half)
        if o is not None and i is not None and abs(o - i) > 1.0:
            outs.append(o)
            ins.append(i)
    if not outs:
        return None
    return float(np.median(outs)), float(np.median(ins))


FOCUS_CHAINS, FOCUS_ELEMENTS = 6, 160  # symmetric assemblies larger than this are drawn one subunit at a time


def _focus_chains(focus: str, sym, bb, sses, nucleic=None) -> set | None:
    """The chains drawn in full: none ('none', or a modest assembly under 'auto'), one protomer (helical
    filaments, cubic cages, or assemblies over FOCUS_CHAINS chains / FOCUS_ELEMENTS elements), or the chains
    named. Helical filaments use their middle subunit, so it has neighbours on both sides. 'auto' never takes one
    subunit of two (a C2 histone octamer is drawn whole), nor one of a ring or other point-group assembly whose
    protomers share a nucleic acid (_na_shared: nucleosome DNA); filaments and cages keep it."""
    present = list(dict.fromkeys(l.chain for l in bb.labels))
    if focus in (None, "none"):
        return None
    if focus not in ("auto", "protomer"):
        named = [c.strip() for c in focus.split(",") if c.strip()]
        missing = [c for c in named if c not in present]
        if missing:
            raise ValueError(f"--focus: no chain {', '.join(missing)}; the structure has {', '.join(present)}")
        return set(named)
    if sym is None:
        if focus == "protomer":
            raise ValueError("--focus protomer needs a symmetric assembly; name chains instead (--focus A)")
        return None
    big = sym.kind in ("H", "T", "O", "I") or len(present) > FOCUS_CHAINS or len(sses) > FOCUS_ELEMENTS
    if focus == "auto" and not big:
        return None
    k = len(sym.protomers) // 2 if sym.kind == "H" else 0
    one = _one_of_each(sym.protomers[k], bb)
    if focus == "auto" and (
        _subunits(sym, one) <= 2 or (sym.kind in ("C", "D") and _na_shared(nucleic, sym, bb))
    ):  # one complex: drawn whole
        return None
    return one


def _one_of_each(chains: list[str], bb) -> set:
    """The smallest repeating unit inside a protomer: one chain of each distinct sequence, those nearest the
    first (a homo-oligomer's protomer gives one chain; a heterodimer's, both partners)."""
    from difflib import SequenceMatcher

    seq = {c: [l.name for l in bb.labels if l.chain == c] for c in chains}
    centre = {c: bb.ca[[k for k, l in enumerate(bb.labels) if l.chain == c]].mean(axis=0) for c in chains}
    groups: list[list[str]] = []
    for c in chains:
        for g in groups:
            if SequenceMatcher(None, seq[g[0]], seq[c], autojunk=False).ratio() > 0.9:
                g.append(c)
                break
        else:
            groups.append([c])
    first = groups[0][0]
    return {min(g, key=lambda c: float(np.linalg.norm(centre[c] - centre[first]))) for g in groups}


def _focus_subset(sses, focus: set, bridges) -> tuple[list, set]:
    """The focus chains' elements plus other chains' elements that complete them: strands paired with theirs by
    backbone H-bond bridges (donor strands, domain swaps, beta-augmentation) and helices packed along theirs
    (features.helix_partners: its neighbours in a bundle; a ring of helices lends two, not all eleven, random-100
    C25)."""
    from .features import helix_partners

    owner = {r: s for s in sses for r in range(s.start, s.end + 1)}
    own = {s.id for s in sses if s.chain in focus}
    guests = set()
    for i, j, _ in bridges:
        a, b = owner.get(i), owner.get(j)
        if a is None or b is None:
            continue
        if a.id in own and b.chain not in focus:
            guests.add(b.id)
        if b.id in own and a.chain not in focus:
            guests.add(a.id)
    for k, partners in sorted(helix_partners(sses).items()):
        if k in own:
            guests |= {m for m in partners if m not in own}
    keep = [s for s in sses if s.id in own or s.id in guests]
    return keep, {s.chain for s in keep if s.chain not in focus}


def _homolog_labels(bb, sses) -> dict[str, str]:
    """Labels shared by copies of a protein, however much of each is modelled (model.shared_labels): one element,
    one label, in every copy."""
    from .layout import provisional
    from .model import homolog_chains, shared_labels

    labels = {k: p.label for k, p in provisional(sses).placed.items()}
    return shared_labels(sses, homolog_chains(bb), labels)


def _subunits(sym, focus: set) -> int:
    """How many subunits like the `focus` chains the assembly holds: its chains over the focus's, a half rounded
    up (Python's round takes 4.5 to 4: nine chains in units of two are said as five subunits, not four)."""
    chains = sum(len(p) for p in sym.protomers)
    return int(math.floor(chains / max(len(focus), 1) + 0.5))


def _context(sym, focus: set, bb, neighbours: bool = True) -> str:
    """One line saying what the figure shows of the assembly."""
    total = len(set(l.chain for l in bb.labels))
    grey = " · grey: neighbouring subunits completing its fold" if neighbours else ""
    if sym is None:
        return f"chain {', '.join(sorted(focus))} of {total} shown in full" + grey
    copies = sum(len(p) for p in sym.protomers) / max(len(sym.protomers), 1) / max(len(focus), 1)
    n = _subunits(sym, focus)
    if sym.kind == "H":
        step = "subunit" if copies <= 1 else f"repeat of {int(round(copies * len(focus)))} chains"
        twist = f"{np.degrees(abs(sym.twist)) * (-1 if sym.left_handed else 1):.1f}".replace("-", "−")
        kind = f"helical, {twist}° twist{sym.hand}, {abs(sym.rise):.1f} Å rise per {step}"
    else:
        kind = f"{sym.label} symmetry"
    return f"1 of {n} subunits shown · {kind}" + grey


_DOMAIN_PULL = 6  # CA-pair-equivalents of attraction between any two elements of one domain


UNIPROT_DOMAINS = ("Domain", "Repeat", "Zinc finger", "DNA binding")
DOMAIN_SOURCES = {"auto": "inferred", "cath": "CATH", "uniprot": "UniProt"}  # --domains -> Layout.domain_source


def _domain_source(spec) -> str:
    """Layout.domain_source of a domains spec: inferred (auto), CATH, UniProt, a layout file's saved source, or ''
    for domains given by name and ranges."""
    if isinstance(spec, str):
        return DOMAIN_SOURCES.get(spec, "")
    return getattr(spec, "source", "") or ""


def domain_note(source: str, named, found: str = "", band: bool = False, absorbed: bool = False) -> str:
    """What the notes say of the domains (Layout.domain_note): automatic ones are inferred, not published (plan 12
    R2); CATH's come with their ids (`found`, cath.domains' note), or the note says why there are none. `band` (a
    projected membrane view, plan 26 ruling T3-a): the note says the domains holding the band's elements are one;
    `absorbed`: residue contacts split some chain but the band joined every domain, so none is drawn (5JDN's two
    halves; ruling R19, review-3 M1)."""
    if source == "CATH":
        return found
    if source != "inferred":
        return ""
    stems = list(dict.fromkeys(name.split(" · ")[0] for name, _ in named))
    joined = "; what crosses the membrane, lies on its faces or sits in it is one domain per chain" if band else ""
    if not stems and absorbed:
        return (
            "no domains drawn: residue contacts split a chain into domains, but each holds what crosses the membrane, "
            "lies on its faces or sits in it, which is one domain per chain"
        )
    if not stems:
        return "no domains inferred: each chain reads as one domain from its residue contacts" + joined
    return f"domains {', '.join(stems)} inferred from residue contacts, not a published domain assignment" + joined


def _uniprot_domains(ann, sses) -> list[tuple[str, list[str]]]:
    """Domain panels from UniProt: each domain takes the elements with most of their residues inside it."""
    chains = {c for _, c, _ in ann.regions}
    out = []
    for f, chain, ks in ann.regions:
        if f.type not in UNIPROT_DOMAINS:
            continue
        inside = set(ks)
        ids = [
            s.id
            for s in sses
            if s.chain == chain and sum(r in inside for r in range(s.start, s.end + 1)) > (s.end - s.start + 1) / 2
        ]
        if ids:
            name = f.description or f.type
            out.append((name if len(chains) == 1 else f"{name} · {chain}", ids))
    return [(name, ids) for name, ids in out]


def _element_ids(spec) -> bool:
    """Are the domains given as element ids (A:3-11: UniProt domains, a layout file) rather than references?"""
    return bool(spec) and all(
        isinstance(ids, list) and all(":" in i and "-" in i and not i.startswith("res:") for i in ids)
        for _, ids in spec
    )


def stale_domain_ids(spec, sses) -> list[tuple[str, str]]:
    """(domain, element id) for each element id the domains name that `sses` lack: none unless the domains are given
    as element ids. A layout file's (SavedDomains) are dropped with a warning; given otherwise, an error."""
    if spec in (None, "auto", "cath", "uniprot") or not _element_ids(spec):
        return []
    known = {s.id for s in sses}
    return [(name, i) for name, ids in spec for i in ids if i not in known]


def stale_domains_note(spec, sses) -> str:
    """The warning for a layout file's domains (SavedDomains) naming elements `sses` lack, '' when none: the command
    line prints it, the app shows it."""
    stale = stale_domain_ids(spec, sses) if isinstance(spec, SavedDomains) else []
    if not stale:
        return ""
    known = {s.id for s in sses}
    empty = [name for name, ids in spec if not any(i in known for i in ids)]
    return (
        "the layout file's domains name element(s) this structure does not have, "
        f"dropped: {', '.join(i for _, i in stale)}"
        + (f" (domain(s) left empty, not drawn: {', '.join(empty)})" if empty else "")
    )


def _domains(
    spec, sses, contacts, rows=(), copies=None, lent=frozenset(), bb=None, conventions=None, band=frozenset()
) -> list[tuple[str, list[str]]]:
    """[(name, element ids)] from 'auto' (chain by chain, copies of a protein (`copies`) alike, chains that only
    lend elements (`lent`) in their mates' domains, never splitting a βα row of `rows`, nor a membrane view's
    `band`: its domains holding those elements joined, features.auto_domains) or [(name, [refs])]; every
    element in at most one domain. `conventions` (names.all_names, the structure bound): the names a reference may
    give with dep: or beta: (before layout no βα row is known, so a row helix's β-scheme name is not)."""
    from .features import auto_domains
    from .layout import provisional, with_conventions

    if not spec:
        return []
    if spec == "auto":
        return auto_domains(sses, contacts, rows, copies, lent, bb=bb, band=band)
    typed = [] if isinstance(spec, SavedDomains) else stale_domain_ids(spec, sses)  # typed ids naming no element
    for name, i in typed:  # A:1-60 is read as the range res:A:1-60 (Q-A19); a range that will not parse is an error
        if RES_RANGE.fullmatch(i.split(":", 1)[1]) is None:
            raise ValueError(
                f"domain {name!r}: {i!r} is no element; write a residue range as res:CHAIN:FIRST-LAST "
                "(e.g. res:A:1-60), or name elements by label (α1, A:G)"
            )
    if _element_ids(spec) and not typed:
        known = {s.id for s in sses}
        stale = stale_domain_ids(spec, sses)
        seen: set[str] = set()  # each element in one domain: the first that names it
        kept = []
        for name, ids in spec:
            kept.append((name, [i for i in dict.fromkeys(ids) if i in known and i not in seen]))
            seen.update(kept[-1][1])
        if stale and not getattr(spec, "quiet", False):
            print(f"foldmap: warning: {stale_domains_note(spec, sses)}", file=sys.stderr)
        return [(name, ids) for name, ids in kept if ids]  # already element ids (UniProt domains, a layout file)
    lay, out, seen = with_conventions(provisional(sses), sses, conventions), [], set()
    for name, refs in spec:
        ids = []
        for ref in refs:
            ids += [k for k in resolve(lay, sses, ref) if k not in seen and k not in ids]
        if not ids:
            raise ValueError(
                f"domain {name!r} has no elements: {', '.join(refs)} already belong to an earlier domain; "
                "give each domain elements of its own"
            )
        seen |= set(ids)
        out.append((name, ids))
    return out


def _cys_place(bb, sses, r: int):
    """(element id, fraction along it) for a cysteine: its own element, or the nearest element end."""
    chain = bb.labels[r].chain
    mine = [s for s in sses if s.chain == chain]
    if not mine:
        return None
    s = min(mine, key=lambda s: 0 if s.start <= r <= s.end else min(abs(r - s.start), abs(r - s.end)))
    if s.start <= r <= s.end:
        return s.id, (r - s.start) / max(s.end - s.start, 1)
    return s.id, (1.0 if r > s.end else 0.0)


def _bridge_springs(bb, sses, links) -> list[tuple[str, float, str, float]]:
    out = []
    for i, j in links.disulfides:
        a, b = _cys_place(bb, sses, i), _cys_place(bb, sses, j)
        if a and b and a[0] != b[0]:
            out.append((a[0], a[1], b[0], b[1]))
    return out


def _bridge_contacts(bb, sses, links) -> dict[frozenset[str], int]:
    """A disulfide ties two elements together like a full contact: the cysteine's own element, or the element
    nearest it along its chain when it sits in a loop."""

    def element(r: int):
        chain = bb.labels[r].chain
        mine = [s for s in sses if s.chain == chain]
        if not mine:
            return None
        return min(mine, key=lambda s: 0 if s.start <= r <= s.end else min(abs(r - s.start), abs(r - s.end))).id

    out: dict[frozenset[str], int] = {}
    for i, j in links.disulfides:
        a, b = element(i), element(j)
        if a and b and a != b:
            out[frozenset((a, b))] = out.get(frozenset((a, b)), 0) + 20  # as strong as the contact pull gets
    return out


def _focus(highlight: str, bb, sym) -> set[str] | None:
    """Chains to keep in colour: the deposited ones, the first protomer, or the user's list."""
    chains = {l.chain for l in bb.labels}
    if highlight == "none":
        return None
    if highlight == "asu":
        if bb.asu_chains and bb.asu_chains != chains:
            return set(bb.asu_chains)
        return set(sym.protomers[0]) if sym else None  # everything was deposited: the first protomer stands in
    if highlight == "protomer":
        return set(sym.protomers[0]) if sym else None
    wanted = set(highlight.split(","))
    if not wanted <= chains:
        raise ValueError(f"highlight names chains that are not here: {', '.join(sorted(wanted - chains))}")
    return wanted


def make_figure_and_loops(
    path,
    mode: str = "projected",
    rotate: float = 0.0,
    flip_v: bool = False,
    title: str | None = None,
    *,
    palette: str | None = None,
    style: str | None = None,
    loops: str | None = None,
    dna: bool = True,
    look: Style | None = None,
    symmetry: str = "auto",
    protomers: list[list[str]] | None = None,
    symmetry_tol: float = TOLERANCE,
    assembly: str = "auto",
    **adjust,
):
    model = read_model(path, assembly)  # once: the layout and the hydrogen-bond overlay
    layout, sses, bb = make_layout(
        path, mode, rotate, flip_v, dna, look, symmetry, protomers, symmetry_tol, assembly, model=model, **adjust
    )
    routed = route_loops(layout, sses, bb)
    fig = draw(layout, routed, sses, title, model=model, look=look, palette=palette, style=style, loop_style=loops)
    return fig, routed


def make_figure(path, *args, **kwargs):
    return make_figure_and_loops(path, *args, **kwargs)[0]


def adjustments(args) -> dict:
    rename, swap, move = {}, [], []
    for item in args.rename:
        ref, sep, name = item.partition("=")
        if not sep or not name:
            raise ValueError(f"--rename needs REF=NAME, got {item!r}")
        rename[ref] = name
    for item in args.swap:
        parts = [p for p in item.split(",") if p]
        if len(parts) != 2:
            raise ValueError(f"--swap needs two elements separated by a comma, got {item!r}")
        swap.append((parts[0], parts[1]))
    for item in args.move:
        ref, sep, delta = item.partition("=")
        d = [float(x) for x in delta.split(",")] if sep else []
        if len(d) != 2:
            raise ValueError(f"--move needs REF=DX,DY, got {item!r}")
        move.append((ref, (d[0], d[1])))
    return {
        "rename": rename,
        "swap": swap,
        "move": move,
        "up": _vector(args.up, "up"),
        "view": _vector(args.view, "view"),
    }


def _layout_data(layout_file) -> bool:
    """Whether a layout file (None: none given) names a data table, which --data-column may name a column of."""
    if not layout_file:
        return False
    from . import layoutfile

    lay = layoutfile.load(layout_file).get("layout")
    return isinstance(lay, dict) and lay.get("data") is not None


def _spec_of(args) -> FigureSpec:
    """The figure the command line asks for: a layout file (if given) overlaid by the options given."""
    shorthand = [
        f"{k}={v}" for k, v in (("palette", args.palette), ("fill", args.style), ("loops", args.loops)) if v is not None
    ]
    edits = adjustments(args)
    try:
        domains = _domain_spec(args.domain)
    except ValueError as err:  # --domain N= : a usage error, as a reference that cannot be drawn is
        raise UsageError(f"error: {err}") from None
    if args.data_column and not args.data and not _layout_data(args.layout_file):
        raise UsageError("error: --data-column needs --data FILE.csv (the table it is a column of)")
    spec = figure_spec(
        args.structure,
        theme=args.theme,
        style=shorthand + args.set,
        style_file=args.style_file,
        layout_file=args.layout_file,
        title=args.title,
        mode=args.mode,
        rotate=args.rotate,
        flip_v=args.flip_v or None,
        dna=False if args.no_dna else None,
        symmetry=args.symmetry,
        symmetry_tol=args.symmetry_tol,
        protomers=args.protomers,
        assembly=args.assembly,
        up=edits["up"],
        view=edits["view"],
        domains=domains or args.domains,
        msa=args.msa,
        msa_reference=args.msa_reference,
        membrane=args.membrane,
        uniprot=args.uniprot,
        chains=args.chains,
        focus=args.focus,
        rename=edits["rename"],
        swap=edits["swap"],
        move=edits["move"],
        colors=getattr(args, "color", None) or [],
        widths=getattr(args, "width", None) or [],
        names=args.names,
        like=_like_of(args),
        data=args.data,
        data_column=args.data_column,
        crosslinks=args.crosslinks,
    )
    for problem in check_references(spec):  # pass 1, before the structure is read (Q-A19); make_layout does pass 2
        if problem.fatal:
            raise BadReference(problem)
        print(f"foldmap: warning: {problem}", file=sys.stderr)
    return spec


def _like_of(args) -> dict | None:
    """--like and its options as FigureSpec.like (cli.make_layout's like); None without --like."""
    given = {"chain": args.like_chain, "match": args.like_match, "layout": args.like_layout}
    if not args.like:
        stray = [f"--like-{k}" for k, v in given.items() if v]
        if stray:
            raise UsageError(f"{', '.join(stray)} needs --like (the reference structure)")
        return None
    return {
        "structure": args.like,
        "chain": given["chain"],
        "match": given["match"] or "auto",
        "layout": given["layout"],
    }


def experimental_method(path) -> str | None:
    """The experimental method a file declares (X-RAY DIFFRACTION, ELECTRON MICROSCOPY...); None for models."""
    import gemmi

    return dict(gemmi.read_structure(str(path)).info).get("_exptl.method") or None


def _has(module: str) -> bool:
    """Whether an optional dependency is installed."""
    from importlib.util import find_spec

    return find_spec(module) is not None


def styles_help() -> str:
    def kind(key: str) -> str:
        if key == "cyclic":
            return "one of unrolled, rosette (rosette: a C4 or C7 assembly drawn whole, seen down its axis)"
        if key in CHOICES:
            return "one of " + ", ".join(CHOICES[key])
        if key in SCALES:
            return "positive number, 1 = default size"
        if key == "na_backbone":
            return "modality, chain, ink or a colour"
        if key in ("dna_color", "rna_color", "background", "ink"):
            return "a colour (#rrggbb or a name)"
        if key == "hbonds":
            return "off, na, interchain or chains:A+B,C+D (interactions=hbond with interactions_between=...)"
        if key == "interactions":
            return INTERACTION_FORMS
        if key == "interactions_between":
            return BETWEEN_FORMS + " (auto: na when the figure holds nucleic acid, else interchain)"
        if key == "interaction_max":
            return "whole number, at most this many lines of each kind (H-bonds: hbond_max)"
        if key == "hbond_color":
            return "a colour (#rrggbb or a name)"
        if key == "hbond_max":
            return "whole number, at most this many lines (the shortest on the page kept, not the shortest in Å)"
        if key == "backbone_colors":
            return "KEY=COLOUR,... (KEY: a chain id, protein, dna or rna)"
        if key == "mark":
            return "REF[=COLOUR],... (elements in a colour of their own)"
        if key == "element_colors":
            return "REF=COLOUR,... (--color; REF an element, loop:REF, na:CHAIN or lig:CODE)"
        if key == "element_widths":
            return "REF=FACTOR,... (--width; 0.25 to 4, the same references)"
        return "true/false"

    lines = ["Themes (--theme NAME) change only how things are drawn:"]
    for name, look in THEMES.items():
        changed = {k: v for k, v in vars(look).items() if v != getattr(Style(), k)}
        lines.append(f"  {name:13s}" + (", ".join(f"{k}={v}" for k, v in changed.items()) or "the defaults"))
    lines.append("\nTheme keys (--set KEY=VALUE, or KEY: VALUE in a --style-file):")
    lines += [f"  {k:17s} {kind(k)}" for k in THEME_KEYS]
    lines.append("\nLayout keys (same syntax; these move elements):")
    lines += [f"  {k:17s} {kind(k)}" for k in LAYOUT_KEYS]
    lines.append("\nDefaults:\n" + "\n".join("  " + l for l in describe(Style()).splitlines()))
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    from . import __version__

    parser = argparse.ArgumentParser(prog="foldmap", description="Publication-grade protein topology diagrams.")
    parser.add_argument("--version", action="version", version=f"foldmap {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("summary", help="print secondary structure elements and sheet ordering")
    p.add_argument("structure")
    sub.add_parser("styles", help="list themes and every style key")
    sub.add_parser("app", help="open the Streamlit app in your browser (needs: pip install streamlit)")
    up = sub.add_parser("uniprot", help="UniProt annotation: features, domains, PDB structures and the AlphaFold model")
    up.add_argument("source", help="a UniProt accession (P04637), a PDB ID or a structure file")
    sub.add_parser(
        "mcp", help="run the MCP server (stdio) so AI agents can use foldmap (needs: pip install foldmap[mcp])"
    )
    sq = sub.add_parser("sequence", help="the sequence (or an alignment) with secondary structure drawn on top")
    sq.add_argument("structure")
    sq.add_argument("-o", "--output", action="append", required=True, help=".svg, .pdf or .png (repeatable)")
    sq.add_argument("--chain", help="protein chain to show (default: the first)")
    sq.add_argument("--msa", metavar="ALIGNMENT", help="multiple sequence alignment (FASTA, Clustal or Stockholm)")
    sq.add_argument("--reference", metavar="NAME", help="alignment row that is the structure (default: best match)")
    sq.add_argument("--columns", type=int, default=60, help="residues (alignment columns) per row (default 60)")
    sq.add_argument("--full-sequence", action="store_true", help="include unmodelled residues at both ends")
    sq.add_argument(
        "--ss-colour",
        choices=["figure", "black"],
        default="figure",
        help="colour elements as in the topology figure (default) or plain black",
    )
    sq.add_argument(
        "--similarity",
        type=float,
        default=0.7,
        metavar="F",
        help="fraction of sequences that must agree for a column to count as similar (default 0.7)",
    )
    sq.add_argument("--no-conservation-bar", action="store_true", help="leave out the conservation bars")
    sq.add_argument("--theme", "--preset", dest="theme", choices=sorted(THEMES))
    sq.add_argument("--style-file", metavar="YAML")
    sq.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    sq.add_argument("--title")
    sq.add_argument("--uniprot", metavar="auto|ACCESSION", help="add a UniProt feature track (domains, sites)")
    sq.add_argument(
        "--data",
        metavar="FILE.csv",
        help="your own per-residue values (as plot's --data), drawn on the track as the topology figure draws them",
    )
    sq.add_argument("--data-column", metavar="NAME", help="the --data column to read (default: its first of numbers)")
    q = sub.add_parser("plot", help="draw the topology figure")
    iq = sub.add_parser("interactive", help="an HTML page linking the topology, the contact map and the 3D model")
    iq.add_argument("--viewer", choices=["molstar", "3dmol"], default="molstar", help="3D viewer (default Mol*)")
    cq = sub.add_parser(
        "contacts",
        help="check the contact beads of the figure against the coordinates (exit 1 on a false, unmatched or "
        "wrong-colour bead)",
        description="Draw the figure as `plot` would and check every protein-nucleic acid contact bead against the "
        "deposited atoms (an independent atom search). Prints a summary and one line per problem bead; exits 1 "
        "on a false, unmatched or wrong-colour bead, 0 otherwise.",
    )
    kq = sub.add_parser(
        "check",
        help="check the figure: its science against the coordinates (sheet neighbours, barrel closure and strand "
        "count, strands drawn beside their partner) and its page (criteria A1-A20; exit 1 on a problem)",
        description="Draw the figure as `plot` would and check what it says about the structure against the "
        "coordinates, apart from how the page looks: S1, two strands drawn side by side in a sheet (adjacent "
        "columns, level) are joined by a bridge, and two in one column (collinear) never overlap; S2, a barrel's "
        "closure is joined; S3, a barrel's strand count is the cycle its pairing makes, and no sheet drawn open "
        "holds a ring of 5 or more; S6, a strand off a sheet's run of paired strands stands beside its partner, "
        "level with it, or its pairing is named in the notes. S4, a linked residue's star stands at its "
        "residue (linked-residue stars only: not ligand tethers, glycans, disulfide ends or UniProt marks); S5 "
        "(automatic domains) is not measured here. Then the figure criteria A1-A20 the random-100 campaign judges "
        "every figure by, with the same limits: route fallbacks, loops over labels, overlapping elements, text "
        "collisions, text off the page, repeated labels, chains left out without a note, N->C and sheet order, "
        "copies drawn alike, domains (with --domains), membrane (with --membrane), ligand tethers; a --membrane or "
        "--domains figure is judged as the campaign's membrane or domains view. Those that need "
        "more than one drawing or a tool foldmap does not use (process exit and time, the app, mkdssp, a second "
        "seed, a saved-layout round trip, RNApolis, a two-state partner) are named 'not measured here' with the "
        "reason. With --like REFERENCE, then H1-H5 of the figure drawn like the reference: H1, drawn alike exactly "
        "when one fold; H2, mates in the same places (median and 90th percentile after the best rigid 2D motion); "
        "H3, the map reads the same way (page dx, dy signs of mate pairs; each mated sheet's order); H4, mates carry "
        "one name; H5, it passes every S- and A-check the structure drawn alone passes. H2 and H3 also give the "
        "two drawn independently, for comparison; their limits are frozen on a homologue tuning set. Prints one line "
        "per check and one per problem; exits 1 on a problem, 0 otherwise.",
    )
    for cmd in (q, iq, cq, kq):
        q = cmd
        q.add_argument("structure")
        q.add_argument(
            "-o",
            "--output",
            action="append",
            required=cmd not in (cq, kq),
            default=[],
            help=".svg, .pdf or .png (repeatable)"
            + ("; optional here: also save the figure checked" if cmd in (cq, kq) else ""),
        )
        q.add_argument("--mode", choices=["projected", "stack"], help="projected (default) or stack")
        q.add_argument("--rotate", type=float, help="turn the picture counter-clockwise (degrees)")
        q.add_argument("--layout-file", metavar="YAML", help="start from a saved layout (theme, style, edits)")
        q.add_argument("--save-layout", metavar="YAML", help="write everything that defines this figure to a file")
        q.add_argument("--flip-v", action="store_true", help="mirror top/bottom (changes handedness)")
        q.add_argument("--title")
        q.add_argument(
            "--theme",
            "--preset",
            dest="theme",
            choices=list(THEMES),
            help="a ready-made look (default: publication); see `foldmap styles`",
        )
        q.add_argument("--style-file", help="YAML file of style keys (may name a theme); see `foldmap styles`")
        q.add_argument(
            "--set",
            action="append",
            default=[],
            metavar="KEY=VALUE",
            help="one style key, e.g. helix_scale=0.8 or loops=curved (repeatable; wins over the rest)",
        )
        q.add_argument("--palette", choices=list(PALETTES), help="chain colours (same as --set palette=...)")
        q.add_argument("--style", choices=STYLES, help="element fill (same as --set fill=...)")
        q.add_argument("--loops", choices=LOOP_STYLES, help="right-angled or smooth loops (same as --set loops=...)")
        q.add_argument("--no-dna", action="store_true", help="leave out DNA/RNA even when the file has it")
        q.add_argument(
            "--symmetry",
            metavar="auto|off|Cn|Dn|helical",
            help="auto finds cyclic, dihedral or helical symmetry and draws every protomer alike, side by side; "
            "name one (C3, D2, helical) to insist on it; off to ignore it",
        )
        q.add_argument(
            "--protomers",
            metavar="A,B;C,D",
            help="your own protomers, in order around the axis (chains separated by commas, protomers by ;)",
        )
        q.add_argument(
            "--assembly",
            metavar="auto|asu|ID",
            help="auto builds the file's first biological assembly when it adds copies (up to 24 chains), or "
            "when it is one of several complexes deposited together; "
            "asu keeps the file as deposited; or give an assembly id from the file",
        )
        q.add_argument(
            "--names",
            action="append",
            default=[],
            metavar="SPEC",
            help="conventional names: deposited (the file's own helix ids, e.g. globin A-H), beta (strands "
            "β1…βn, and αn after βn in a βα row) or a YAML file of element and stem names (examples in docs/names); "
            "repeatable, later ones win, --rename wins over all",
        )
        q.add_argument(
            "--rename",
            action="append",
            default=[],
            metavar="REF=NAME",
            help="your own element name, e.g. res:167-182=Gd or A:B=B' (repeatable); here and in --swap, --move, "
            "--domain and mark an element is a label, chain:label, #k, res:CHAIN:FIRST-LAST, or a name --names "
            "would give it: dep:NAME (deposited, e.g. dep:E) or beta:NAME (β-scheme, e.g. beta:α4)",
        )
        q.add_argument(
            "--swap",
            action="append",
            default=[],
            metavar="REF,REF",
            help="exchange the places of two elements, e.g. α1,α3 or A:C,A:D (repeatable)",
        )
        q.add_argument(
            "--move",
            action="append",
            default=[],
            metavar="REF=DX,DY",
            help="nudge an element (page units, one strand spacing ~1.1), e.g. α2=1,-2 (repeatable)",
        )
        q.add_argument(
            "--color",
            "--colour",
            dest="color",
            action="append",
            default=[],
            metavar="REF=COLOUR",
            help="an element's own colour, winning over the colouring, e.g. α1=red or res:A:20-25=#1f77b4; REF may "
            "also be loop:REF (the loop after REF), na:CHAIN (a nucleic-acid backbone) or lig:CODE (repeatable; "
            "style key element_colors)",
        )
        q.add_argument(
            "--width",
            action="append",
            default=[],
            metavar="REF=FACTOR",
            help="an element's drawn width times FACTOR (0.25 to 4): a helix's diameter, an arrow's width, a loop's "
            "stroke, a backbone's weight, e.g. α1=2 or loop:β2=1.5 (repeatable; style key element_widths)",
        )
        q.add_argument(
            "--domain",
            action="append",
            default=[],
            metavar="NAME=REF[,REF...]",
            help="a named domain panel, e.g. ZPN=res:A:331-440 (repeatable)",
        )
        q.add_argument(
            "--domains",
            choices=["auto", "cath", "uniprot"],
            help="auto: inferred from residue contacts (D1, D2...); cath: CATH's published domains (PDBe, needs "
            "network); uniprot: UniProt's domains (with --uniprot)",
        )
        q.add_argument(
            "--uniprot",
            metavar="auto|ACCESSION",
            help="annotate from UniProt: auto finds each chain's entry (SIFTS, or the AlphaFold accession)",
        )
        q.add_argument(
            "--membrane",
            choices=["off", "auto", "flip"],
            help="draw the lipid bilayer: auto finds it (OPM dummy atoms, or estimated); flip: the same with the "
            "two sides swapped; default off",
        )
        q.add_argument(
            "--focus",
            metavar="auto|none|protomer|A,B",
            help="large assemblies: draw one subunit in full plus what its neighbours add to its fold "
            "(auto: helical filaments, cages and assemblies over 6 chains / 160 elements, never one of two "
            "subunits nor a ring whose subunits share a nucleic acid)",
        )
        q.add_argument(
            "--chains",
            metavar="A,B",
            help="draw only these chains (e.g. one subunit of a big complex); protein and nucleic-acid chains are "
            "each restricted only if some are named",
        )
        q.add_argument(
            "--like",
            metavar="REFERENCE",
            help="draw this structure like a reference (a file or PDB id): each element matched to one of the "
            "reference's takes its place, angle, name and colour; one protein chain of each (--chains, --like-chain)",
        )
        q.add_argument(
            "--like-chain",
            metavar="CHAIN",
            help="the reference's chain to draw like (default: its first protein chain)",
        )
        q.add_argument(
            "--like-match",
            metavar="auto|sequence|structure|ALN",
            help="how the chains are aligned: auto (sequence at 40%% identity, else structure), sequence, structure, "
            "or a two-row alignment file (the reference's chain first)",
        )
        q.add_argument(
            "--like-layout", metavar="YAML", help="a saved layout of the reference: its options and edits are used"
        )
        q.add_argument(
            "--msa",
            metavar="ALIGNMENT",
            help="alignment (or a ConSurf grades table) for conservation colouring (--theme conservation)",
        )
        q.add_argument("--msa-reference", metavar="NAME", help="alignment row that is the structure (default: best)")
        q.add_argument(
            "--data",
            metavar="FILE.csv",
            help="your own per-residue values: a table with chain, residue (45, 100A) and value columns; every row "
            "that names no residue of the structure is said with its row number",
        )
        q.add_argument(
            "--data-column",
            metavar="NAME",
            help="the --data column to read (default: its first column of numbers); a column that is not all "
            "numbers is read as categories",
        )
        q.add_argument(
            "--crosslinks",
            metavar="FILE.csv",
            help="crosslinks: a table with chain1, residue1, chain2, residue2 and an optional distance column",
        )
        q.add_argument(
            "--up", metavar="X,Y,Z", help="3D direction to put at the top of the page (e.g. a membrane normal)"
        )
        q.add_argument("--view", metavar="X,Y,Z", help="3D direction pointing at the viewer")
        q.add_argument(
            "--symmetry-tol",
            type=float,
            metavar="Å",
            help=f"how far copies may differ, CA RMSD (default {TOLERANCE:g}; raise it for low-resolution X-ray)",
        )
    pq = sub.add_parser(
        "compare",
        help="two structures of one fold on one page: A on the left, B drawn like A on the right",
        description="Draw one protein chain of A, and one of B drawn like it (as plot --like draws it), side by side "
        "at one page scale, with one legend and the note under both: how many of B's elements were matched, by "
        "which alignment, or why B was drawn on its own (not the same fold). Each panel's ids start a/ and b/.",
    )
    pq.add_argument("structure_a", metavar="A", help="the reference: a file or PDB id")
    pq.add_argument("structure_b", metavar="B", help="the structure drawn like A: a file or PDB id")
    pq.add_argument(
        "-o", "--output", action="append", required=True, default=[], help=".svg, .pdf or .png (repeatable)"
    )
    pq.add_argument("--chain-a", metavar="CHAIN", help="A's protein chain (default: its first)")
    pq.add_argument("--chain-b", metavar="CHAIN", help="B's protein chain (default: its first)")
    pq.add_argument(
        "--match",
        metavar="auto|sequence|structure|ALN",
        default="auto",
        help="how the chains are aligned: auto (sequence at 40%% identity, else structure), sequence, structure, "
        "or a two-row alignment file (A's chain first)",
    )
    pq.add_argument(
        "--assembly",
        metavar="auto|asu|ID",
        default="auto",
        help="how both files are read: auto builds each file's first biological assembly when it adds copies (its "
        "chains then renamed A1, A2, ...); asu keeps the files as deposited, under their deposited chain names; or "
        "give an assembly id",
    )
    pq.add_argument(
        "--json",
        metavar="PATH",
        help="also write what was matched: the alignment (method, TM-score, RMSD, identity, residue pairs), the "
        "mates, the pieces of a helix A has whole, B's extras and A's elements B lacks",
    )
    pq.add_argument("--theme", "--preset", dest="theme", choices=list(THEMES), help="a ready-made look")
    pq.add_argument("--style-file", help="YAML file of style keys (may name a theme)")
    pq.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="one style key (repeatable)")
    _diff_parser(sub)
    _data_parser(sub)
    return parser


def _diff_parser(sub) -> None:
    dq = sub.add_parser(
        "diff",
        help="what changed between two structures of one protein (two states), drawn on their compare page",
        description="Draw one protein chain of A and the same chain of B drawn like it, as compare draws them, with "
        "what changed: an element only one structure has as a dashed outline in the other's panel, B coloured by "
        "each residue's CA displacement after superposition on the matched core, and the rest (secondary-structure "
        "changes, sheet pairing gained, lost or shifted, domain turns, interface and contact changes) named in the "
        "note under both. The chains are aligned by sequence when they are one protein, else by structure. Each "
        "panel's ids start a/ and b/.",
    )
    dq.add_argument("structure_a", metavar="A", help="the reference state: a file or PDB id")
    dq.add_argument("structure_b", metavar="B", help="the other state: a file or PDB id")
    dq.add_argument(
        "-o", "--output", action="append", required=True, default=[], help=".svg, .pdf or .png (repeatable)"
    )
    dq.add_argument("--chain-a", metavar="CHAIN", help="A's protein chain (default: its first)")
    dq.add_argument("--chain-b", metavar="CHAIN", help="B's protein chain (default: its first)")
    dq.add_argument(
        "--match",
        metavar="auto|sequence|structure|ALN",
        default="auto",
        help="how the chains are aligned: auto (sequence when they are one protein, else structure), sequence, "
        "structure, or a two-row alignment file (A's chain first)",
    )
    dq.add_argument(
        "--assembly",
        metavar="auto|asu|ID",
        default="auto",
        help="how both files are read: auto (each file's first biological assembly when it adds copies), asu (as "
        "deposited) or an assembly id",
    )
    dq.add_argument(
        "--domain",
        action="append",
        default=[],
        metavar="NAME=RANGES",
        help="a domain of A's chain whose rigid turn the note gives, e.g. LID=122-159 or CORE=1-29,60-121 (repeatable)",
    )
    dq.add_argument(
        "--json",
        metavar="PATH",
        help="also write what changed: compare's matching (alignment, mates, pieces, extras, only_a), then "
        "ss_changes, pairing, displacement, core (how many pairs the core fit kept), domain_turns, domain_angles (two "
        "domains or more: each against the first, each fitted on its own), contacts and partners (B's partner "
        "chains, each named by A's chain in its place)",
    )
    dq.add_argument(
        "--check",
        action="store_true",
        help="also judge each panel as a figure of its own (H5: it passes every check its chain drawn alone passes); "
        "exit 1 when one does not",
    )
    dq.add_argument("--theme", "--preset", dest="theme", choices=list(THEMES), help="a ready-made look")
    dq.add_argument("--style-file", help="YAML file of style keys (may name a theme)")
    dq.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="one style key (repeatable)")


def diff_domains(given: list[str]) -> dict[str, list[str]]:
    """--domain NAME=RANGES (each a comma-separated list of ranges) as diff_structures' domains, in the order
    given; a UsageError for one without a name or ranges, a name given twice, or a name with "/" in it (the --json
    domain_angles key is FIRST/NAME)."""
    out: dict[str, list[str]] = {}
    for text in given:
        name, _, ranges = text.partition("=")
        name = name.strip()
        parts = [r.strip() for r in ranges.split(",") if r.strip()]
        if not name or not parts:
            raise UsageError(f"--domain {text!r}: give NAME=RANGES, e.g. LID=122-159 or CORE=1-29,60-121")
        if "/" in name:
            raise UsageError(f"--domain {name}: a domain name cannot hold / (domain_angles names FIRST/NAME)")
        if name in out:
            raise UsageError(f"--domain {name} is given twice")
        out[name] = parts
    return out


# plot's options `foldmap data` does not take, as plot's defaults (_spec_of reads them)
_DATA_DEFAULTS = {
    "theme": None,
    "style_file": None,
    "palette": None,
    "style": None,
    "loops": None,
    "title": None,
    "mode": None,
    "rotate": None,
    "flip_v": False,
    "no_dna": False,
    "symmetry_tol": None,
    "protomers": None,
    "swap": [],
    "move": [],
    "up": None,
    "view": None,
    "msa": None,
    "msa_reference": None,
    "like": None,
    "like_chain": None,
    "like_match": None,
    "like_layout": None,
    "data": None,
    "data_column": None,
    "crosslinks": None,
}


def _data_parser(sub) -> None:
    from .datafile import SECTIONS, TABLES

    dq = sub.add_parser(
        "data",
        help="what the figure says about the structure, as JSON or CSV tables: elements and their labels, sheets "
        "(order, register, barrel n and shear S), disulfides, symmetry, membrane, ligand contacts, interfaces",
        description="Lay the structure out as `plot` would (no drawing) and write what the figure says about it as "
        "JSON (schema version 1; --schema prints it): the elements with their residues, the figure's labels and "
        "page places; the sheets with their strands in the order drawn, directions, each pairing's kind and residue "
        "register from the DSSP bridges, accessory strands, and a closed sheet's strand number n and shear number "
        "S; disulfides; symmetry and stoichiometry; with --membrane, the transmembrane segments and sides; the "
        "ligands, metal ions, glycans and fused residues (a chromophore) with the residues and elements they touch; "
        "the chain pairs in contact with their elements in contact; with --sections ...,interactions, every "
        "interaction between two chains (H-bonds, salt bridges, stacking, cation-pi, hydrophobic); the nucleic acid's "
        "strands and base pairs (those the figure draws) with their stems and pseudoknot bracket levels. Numbers are "
        "rounded to 0.001; the same structure and options give the same bytes. `-o out.csv --table NAME` writes one "
        "table as CSV; `-o DIR/` writes every table of the sections as DIR/NAME.csv; `-o out.dbn` (or .bpseq, .ct; "
        "--format for standard output) writes the nucleic acid as dot-bracket, BPSEQ or CT.",
    )
    dq.add_argument("structure", nargs="?", help="a file, PDB ID or UniProt accession (not needed with --schema)")
    dq.add_argument(
        "-o",
        "--output",
        metavar="JSON|CSV|DIR/|DBN|BPSEQ|CT",
        help="where to write it (default: standard output); a directory (ending in /, or one that exists) gets one "
        "CSV per table",
    )
    dq.add_argument("--schema", action="store_true", help="print the JSON Schema of the output and exit")
    dq.add_argument(
        "--format",
        choices=["json", "dbn", "bpseq", "ct"],
        help="json (default), or the nucleic acid as dot-bracket (dbn), BPSEQ or CT; else from -o's suffix",
    )
    dq.add_argument(
        "--sections",
        metavar="A,B",
        help=f"only these sections (default all but interactions): {', '.join(SECTIONS)}",
    )
    dq.add_argument(
        "--table",
        metavar="NAME",
        help=f"write this table as CSV (its section is computed unless --sections says otherwise): {', '.join(TABLES)}",
    )
    dq.add_argument("--layout-file", metavar="YAML", help="a saved layout: its options and edits (names, domains)")
    dq.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help="a layout style key, e.g. helices_310=false"
    )
    dq.add_argument("--chains", metavar="A,B", help="only these chains, as plot --chains")
    dq.add_argument("--assembly", metavar="auto|asu|ID", help="as plot --assembly")
    dq.add_argument("--symmetry", metavar="auto|off|Cn|Dn|helical", help="as plot --symmetry")
    dq.add_argument("--focus", metavar="auto|none|protomer|A,B", help="as plot --focus")
    dq.add_argument("--names", action="append", default=[], metavar="SPEC", help="as plot --names (repeatable)")
    dq.add_argument("--rename", action="append", default=[], metavar="REF=NAME", help="as plot --rename (repeatable)")
    dq.add_argument("--domain", action="append", default=[], metavar="NAME=REF[,REF...]", help="as plot --domain")
    dq.add_argument("--domains", choices=["auto", "cath", "uniprot"], help="as plot --domains")
    dq.add_argument("--uniprot", metavar="auto|ACCESSION", help="as plot --uniprot (for --domains uniprot)")
    dq.add_argument("--membrane", choices=["off", "auto", "flip"], help="as plot --membrane")
    dq.set_defaults(**_DATA_DEFAULTS)


def resolve_structure(text: str) -> str:
    if Path(text).expanduser().is_file():
        return text
    from .fetch import fetch

    return str(fetch(text))  # a PDB ID or UniProt accession: downloaded once, cached


def _cmd_summary(args: argparse.Namespace) -> int:
    print(summarize(args.structure))
    return 0


def _cmd_styles(args: argparse.Namespace) -> int:
    print(styles_help())
    return 0


def _cmd_app(args: argparse.Namespace) -> int:
    if not _has("streamlit"):
        print("foldmap: the app needs Streamlit: pip install 'foldmap[app]'", file=sys.stderr)
        return 1
    from .app import run

    run()
    return 0


def _cmd_uniprot(args: argparse.Namespace) -> int:
    from .uniprot import report

    print(report(args.source))
    return 0


def _cmd_mcp(args: argparse.Namespace) -> int:
    if not _has("mcp"):
        print("foldmap: the MCP server needs the mcp package: pip install 'foldmap[mcp]'", file=sys.stderr)
        return 1
    from .mcp_server import main as serve

    serve()
    return 0


def _cmd_sequence(args: argparse.Namespace) -> int:
    from .seqplot import draw_sequence
    from .sequence import read_alignment

    if args.data_column and not args.data:
        raise UsageError("error: --data-column needs --data FILE.csv (the table it is a column of)")
    theme = args.theme or "publication"
    look = data_look(resolve_style(theme, args.style_file, args.set), args.data, theme, args.style_file, args.set)
    fig = draw_sequence(
        args.structure,
        chain=args.chain,
        alignment=read_alignment(args.msa) if args.msa else None,
        reference=args.reference,
        columns=args.columns,
        look=look,
        full_sequence=args.full_sequence,
        ss_colour=args.ss_colour,
        threshold=args.similarity,
        conservation_bar=not args.no_conservation_bar,
        title=args.title,
        uniprot=args.uniprot,
        msa=args.msa,
        data=args.data,
        data_column=args.data_column,
    )
    for out in args.output:
        print(f"wrote {save(fig, out)}")
    return 0


def _cmd_figure(args: argparse.Namespace) -> int:
    """plot and interactive: the same figure spec, drawn or written as a page."""
    spec = _spec_of(args)
    if args.command == "interactive":
        from .interactive import write_page

        for n, out in enumerate(args.output):
            made: list[dict] = []
            page = write_page(spec.structure, out, viewer=args.viewer, provenance_out=made, **spec.page_kwargs())
            if n == 0:  # the same figure for every output
                print(provenance_block(made), file=sys.stderr, end="")
            print(f"wrote {page}")
    else:
        drawn = draw_spec(spec)
        scope = drawn_scope(drawn.model, drawn.layout, spec.assembly, spec.focus)
        note = scope_note({**scope, "additives": []})  # the additives hidden are in the list below
        if note:  # an automatic choice drew part of the file
            print(f"foldmap: {note}", file=sys.stderr)
        print(provenance_block(drawn.provenance), file=sys.stderr, end="")
        for out in args.output:
            print(f"wrote {save(drawn.figure, out)}")
    method = experimental_method(spec.structure) if spec.look.color_by == "plddt" else None
    if method:
        print(
            f"foldmap: warning: {spec.structure} is an experimental structure ({method.lower()}); "
            "its B-factor column holds B-factors, not pLDDT confidence",
            file=sys.stderr,
        )
    if args.save_layout:
        laid = drawn if args.command != "interactive" else None  # the figure's own layout, when it has no edits
        print(f"wrote {layoutfile.save(args.save_layout, layout_document(spec, drawn=laid))}")
    return 0


def _cmd_compare(args: argparse.Namespace) -> int:
    """compare: A and B drawn like A on one page (compare.compare); the page written to every output, what was
    matched to --json."""
    from .compare import compare, write_json

    look = resolve_style(args.theme or "publication", args.style_file, args.set)
    made = compare(args.structure_a, args.structure_b, args.chain_a, args.chain_b, args.match, look, args.assembly)
    print(provenance_block(made.provenance), file=sys.stderr, end="")
    for out in args.output:
        print(f"wrote {save(made.figure, out)}")
    if args.json:
        print(f"wrote {write_json(made.data, args.json)}")
    return 0


def _cmd_diff(args: argparse.Namespace) -> int:
    """diff: what changed from A to B drawn on their compare page (compare.diff_page); the page written to every
    output, what changed to --json; --check: H5 for each panel (compare.panel_checks), exit 1 when one fails."""
    from .compare import diff_page, panel_checks, write_json
    from .criteria import LIKE_CHECKS, NOT_JUDGED

    look = resolve_style(args.theme or "publication", args.style_file, args.set)
    domains = diff_domains(args.domain) or None
    made = diff_page(
        args.structure_a, args.structure_b, args.chain_a, args.chain_b, args.match, look, args.assembly, domains
    )
    print(provenance_block(made.provenance), file=sys.stderr, end="")
    for out in args.output:
        print(f"wrote {save(made.figure, out)}")
    if args.json:
        print(f"wrote {write_json(made.data, args.json)}")
    if not args.check:
        return 0
    failed = False
    for side, got in panel_checks(made, args.assembly).items():
        said = NOT_JUDGED if got["pass"] is None else "pass" if got["pass"] else "fail"
        print(f"H5 {LIKE_CHECKS['H5']}, {'left' if side == 'a' else 'right'} panel: {said}")
        if got.get("detail"):
            print(f"  {got['detail']}")
        failed = failed or got["pass"] is False
    return 1 if failed else 0


def _cmd_data(args: argparse.Namespace) -> int:
    """data: the figure `plot` would draw, laid out only, as datafile.structure_data's JSON (or --schema), one CSV
    table, every table, or the nucleic acid as dot-bracket, BPSEQ or CT."""
    import json

    from .datafile import (
        DEFAULT,
        NA_FORMATS,
        SECTIONS,
        TABLES,
        csv_text,
        dumps,
        schema,
        spec_data,
        tables,
        write,
        write_csv,
        write_tables,
    )

    if args.schema:
        print(json.dumps(schema(), indent=1, sort_keys=True, ensure_ascii=False))
        return 0
    if not args.structure:
        raise UsageError("data needs a structure (or --schema)")
    table = args.table
    if table is not None and table not in TABLES:
        raise UsageError(f"error: no data table {table}; the tables are {', '.join(TABLES)}")
    out = args.output
    to_dir = bool(out) and (out.endswith(("/", "\\")) or Path(out).is_dir())
    suffix = Path(out).suffix.lower().lstrip(".") if out and not to_dir else ""
    form = args.format or ({"db": "dbn"}.get(suffix, suffix) if suffix in ("dbn", "db", "bpseq", "ct") else "json")
    if form != "json" and (table is not None or to_dir):
        raise UsageError(f"error: --format {form} writes the nucleic acid as one file: not with --table or -o DIR/")
    if args.sections:
        sections = tuple(s.strip() for s in args.sections.split(",") if s.strip())
    elif form != "json":
        sections = ("structure", "nucleic_acid")
    else:
        sections = (TABLES[table][0],) if table else DEFAULT
    unknown = [s for s in sections if s not in SECTIONS]
    if unknown:
        raise UsageError(f"error: no data section {', '.join(unknown)}; the sections are {', '.join(SECTIONS)}")
    if table is not None and TABLES[table][0] not in sections:
        raise UsageError(f"error: table {table} comes from section {TABLES[table][0]}, which --sections leaves out")
    if form != "json" and "nucleic_acid" not in sections:
        raise UsageError(f"error: --format {form} comes from section nucleic_acid, which --sections leaves out")
    if out and not to_dir and table is None and form == "json" and suffix == "csv":
        raise UsageError(f"error: a .csv holds one table: name it with --table ({', '.join(TABLES)}), or give -o DIR/")
    data = spec_data(_spec_of(args), sections)
    if form != "json":
        text = NA_FORMATS[form](data)
        if out:
            path = Path(out)
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
            print(f"wrote {path}", file=sys.stderr)
        else:
            sys.stdout.write(text)
    elif to_dir:
        written = write_tables(data, out, None if table is None else (table,))
        print(f"wrote {len(written)} table{'s' if len(written) != 1 else ''} to {out}", file=sys.stderr)
    elif table is not None:
        columns, rows = tables(data, (table,))[table]
        if out:
            print(f"wrote {write_csv(columns, rows, out)}", file=sys.stderr)
        else:
            sys.stdout.write(csv_text(columns, rows))
    elif out:
        print(f"wrote {write(data, out)}", file=sys.stderr)
    else:
        sys.stdout.write(dumps(data))
    return 0


def provenance_block(items: list[dict]) -> str:
    """The command line's "How this figure was made:" block (provenance): the heading, then one indented line per
    item, those the command line has always warned of starting "warning: "; '' when the list is empty."""
    from .provenance import HEADING, lines

    return "".join(f"{line}\n" for line in [HEADING, *(f"  - {t}" for t in lines(items, warn=True))]) if items else ""


def _cmd_contacts(args: argparse.Namespace) -> int:
    """Draw the figure `plot` would, check every contact bead against the coordinates (nacheck.verify_contacts) and
    print the counts and one line per problem bead. Exit 1 when a bead is false, unmatched or the wrong colour."""
    drawn = draw_spec(_spec_of(args))
    report = drawn.verify_contacts()
    for out in args.output:
        print(f"wrote {save(drawn.figure, out)}")
    if report is None:
        print("contact beads drawn: 0 (no nucleic acid touched by protein in this figure)")
        return 0
    counts = (
        ("beads drawn", len(report.beads)),
        ("verified", len(report.verified)),
        ("false", len(report.false)),
        ("missing", len(report.missing)),
        ("unmatched", len(report.unmatched)),
        ("wrong colour", len(report.wrong_colour)),
        ("not first", len(report.not_first)),
    )
    print("contact beads: " + ", ".join(f"{name} {n}" for name, n in counts))
    for line in contact_problems(report):
        print(f"  {line}")
    return 1 if report.false or report.unmatched or report.wrong_colour else 0


def _cmd_check(args: argparse.Namespace) -> int:
    """Draw the figure `plot` would and print its science checks (figcheck.science_checks), one line per check and
    one per problem, then its figure criteria A1-A20 (criteria.check_figure), one line per criterion and its detail.
    Drawn like a reference (--like), then H1-H5 (criteria.like_checks), B also drawn alone for H5 and the baseline.
    Exit 1 when any check found a problem or any criterion judged here fails."""
    from .criteria import check_figure, criteria_lines, like_lines
    from .science import science_lines

    spec = _spec_of(args)
    got: dict | None = {} if spec.like else None
    drawn = draw_spec(spec, like_out=got)
    print(provenance_block(drawn.provenance), file=sys.stderr, end="")
    for out in args.output:
        print(f"wrote {save(drawn.figure, out)}")
    found = check_figure(drawn, spec, args.structure)
    for line in science_lines(found["science"], found["science_judged"]) + criteria_lines(found["criteria"]):
        print(line)
    failed = any(c["pass"] is False for c in found["criteria"].values())
    if got is not None:
        like = _like_checks(spec, drawn, found, got)
        for line in like_lines(like):
            print(line)
        failed = failed or any(c["pass"] is False for c in like.values())
    return 1 if any(found["science"].values()) or failed else 0


def _like_checks(spec, drawn, found: dict, got: dict) -> dict:
    """H1-H5 of the figure `drawn` like a reference (criteria.like_checks; `got`: make_layout's like_out, `found`:
    its check_figure): its chain also drawn alone, with the nucleic-acid chains it was asked with, as `foldmap check`
    draws it without --like."""
    from dataclasses import replace

    from .criteria import check_figure, like_checks

    lay = drawn.layout
    chain = drawn.bb.labels[0].chain if len(drawn.bb) else None
    chains = [chain, *(c for c in spec.chains or [] if c != chain)] if chain else spec.chains
    alone_spec = replace(spec, like=None, chains=chains)
    alone = draw_spec(alone_spec)
    reference = got.get("reference")
    return like_checks(
        reference[0] if reference else None,
        lay,
        got.get("alignment"),
        got.get("mates"),
        lay.like_note,
        got.get("like") is not None,
        alone.layout,
        found,
        check_figure(alone, alone_spec, spec.structure),
        refused=bool(got.get("refused")),
    )


def contact_problems(report) -> list[str]:
    """One line per problem of a nacheck.ContactReport: each bead with a problem, then each missing contact."""
    lines = []
    for bead in report.beads:
        if bead.problem:
            where = f" ({bead.pair})" if bead.pair else ""
            lines.append(f"{bead.problem}: {bead.nucleotide} <- {bead.residue}{where}")
    return lines + [f"missing: {miss.nucleotide} <- {miss.pair}" for miss in report.missing]


_COMMANDS: dict[str, Callable[[argparse.Namespace], int]] = {
    "summary": _cmd_summary,
    "styles": _cmd_styles,
    "app": _cmd_app,
    "uniprot": _cmd_uniprot,
    "mcp": _cmd_mcp,
    "sequence": _cmd_sequence,
    "plot": _cmd_figure,
    "interactive": _cmd_figure,
    "contacts": _cmd_contacts,
    "check": _cmd_check,
    "compare": _cmd_compare,
    "diff": _cmd_diff,
    "data": _cmd_data,
}


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if getattr(args, "structure", None):
            args.structure = resolve_structure(args.structure)
        return _COMMANDS[args.command](args)
    except UsageError as err:  # options that do not go together: exit 2, as argparse does
        print(f"foldmap: {'error: ' if isinstance(err, BadReference) else ''}{err}", file=sys.stderr)
        return 2
    except AmbiguousName as err:  # a name meaning different elements in two conventions: a usage error too
        print(f"foldmap: error: {err}", file=sys.stderr)
        return 2
    except (FileNotFoundError, ValueError, RuntimeError) as err:  # RuntimeError: gemmi cannot read the file
        print(f"foldmap: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
