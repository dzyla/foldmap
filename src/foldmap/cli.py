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
from pathlib import Path

import numpy as np

from . import layoutfile
from .criteria import LIMITS
from .dssp import assign_dssp
from .features import helices_on_sheets, sse_contacts
from .frame import (
    Frame,
    RadialFrame,
    dihedral_rings,
    dihedral_rosette_frame,
    dna_radial_frame,
    dna_view_frame,
    protomer_maps,
    rosette_frame,
    symmetric_frame,
    view_frame,
)
from .io import load_backbone, load_links, load_nucleic, read_model
from .layout import SCALE, _snap, accessory_lines, build_layout, na_torn_by_unrolling, resolve, template_rows
from .membrane import FACE_REACH
from .model import AUTO_NT, crop_nucleic, homolog_chains, nucleic_of_chains, nucleic_touching, res_ref
from .names import apply_names, check_spec
from .nucleic import StemView, stem_views, straight_duplex
from .palette import PALETTES
from .render import LOOP_STYLES, STYLES, draw, save
from .route import route_loops
from .science import science_lines, sheet_checks, sheet_judged
from .scope import drawn_scope, scope_note
from .sheets import build_sheets
from .spec import (
    FigureSpec,
    SavedDomains,
    draw_spec,
    figure_spec,
    layout_document,
    parse_domains,
    parse_protomers,
    parse_vector,
)
from .ss import build_sses
from .style import CHOICES, LAYOUT_KEYS, SCALES, THEME_KEYS, THEMES, Style, describe, resolve_style
from .symmetry import TOLERANCE, detect_symmetry, protomer_dyad

_vector, _protomers, _domain_spec = parse_vector, parse_protomers, parse_domains  # older names, used by the app


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
    frame.mirror_x = float(_snap(frame.project(sym.centre)[0, 0] * SCALE))


def _mirrors(two_fold: str, sym, bb) -> bool:
    """Whether Style.two_fold asks for the dyad view (_dyad_frame) of a two-fold assembly: mirror always, unrolled
    never, auto only where every protomer is one protein chain (TrmD, 1UAM; nucleic-acid chains do not count). With
    two or more chains a protomer (haemoglobin, 2HHB: α1β1 | α2β2), interfaces between different chains cross the
    dyad (α1β2, the molecule's main contact) and no mirror line keeps them together: the copies are unrolled."""
    if two_fold != "auto":
        return two_fold == "mirror"
    protein = {l.chain for l in bb.labels}
    return all(sum(c in protein for c in p) == 1 for p in sym.protomers)


ROSETTE_ORDERS = range(3, 9)  # cyclic=rosette: Cn assemblies of these orders are drawn seen down their axis


def _no_rosette(sses, sym, nucleic, mode: str, mem) -> str:
    """Why a Cn assembly (n in ROSETTE_ORDERS) asked for as a rosette (Style.cyclic) is drawn unrolled ('' when it is
    not): what the rosette does not draw yet. A rosette turns protomer 0's drawing for each copy, so everything
    drawn must belong to a protomer: no nucleic acid (its stems are not turned yet), no element of a chain outside
    every protomer (an unpaired chain), and each protomer must hold elements; a stack row and a membrane view are
    not seen down the axis."""
    if mode != "projected":
        return "a rosette is drawn in projected mode only"
    if mem is not None:
        return "a membrane view keeps the bilayer across the page"
    if nucleic is not None and nucleic.strands:
        return "its nucleic acid is not yet drawn in a rosette (draw it with --no-dna for the protein alone)"
    inside = {c for p in sym.protomers for c in p}
    if any(s.chain not in inside for s in sses):
        return "a chain outside every protomer is not yet drawn in a rosette"
    if any(not any(s.chain in p for s in sses) for p in sym.protomers):
        return "a protomer has no helix or strand"
    return ""


def _within_protomer(bridges, sector: dict, bb) -> list:
    """The bridges between residues of one copy (`sector`: chain -> copy, a rosette frame's): a rosette's sheets
    (each copy's own; a sheet two copies share is drawn as each copy's part, the pairings across copies left out)."""
    return [b for b in bridges if sector.get(bb.labels[b.i].chain, -1) == sector.get(bb.labels[b.j].chain, -2)]


def _pairings_across(bridges, sector: dict, bb, sses) -> list[tuple[str, str]]:
    """The strand pairs (element ids, in sequence order) that bridge two copies (`sector`: chain -> copy): a rosette
    draws each copy's own sheets, so these pairings are not drawn (the notes name them)."""
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
    up, right = _snap(axis), _snap(right / np.linalg.norm(right))
    face = np.cross(right, up)  # the view at turn 0: protomer 1 straight to protomer 0's right
    turn = _dyad_turn(np.array(cents[0] + cents[1]), up, right, face)
    if turn:  # turned about the dyad: protomer 1 still to the right, a little behind or in front of protomer 0
        t = np.radians(turn)  # (cos and sin may differ in the last bit between libms: snapped)
        face = _snap(np.cos(t) * face + np.sin(t) * right)
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
        r = rankdata(_snap(d))
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
        page = ranks(pdist(_snap(np.stack([x, y], axis=1))))
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
):
    """`names` (--names): conventional names, each `deposited`, `beta` or a YAML file (names.apply_names), applied
    over the generic labels in order; `rename` wins over all of them."""
    for spec in names or []:
        check_spec(spec)
    model = model or read_model(path, assembly)  # once, for every loader below (`model`: already read, same assembly)
    bb = load_backbone(path, assembly, model=model)
    present = list(dict.fromkeys(l.chain for l in bb.labels))
    if chains and any(c in present for c in chains):  # protein chains named: only those (else every one)
        bb = bb.only([c for c in chains if c in present])
    if (look or Style()).color_by == "conservation" and not msa:
        raise ValueError("conservation colouring needs an alignment: add --msa ALIGNMENT")
    dssp = assign_dssp(bb)
    matched: dict = {}  # helix pieces cut or joined so copies match (ss._copy_breaks): told by scope.drawn_scope
    sses = build_sses(bb, dssp.ss, short_helices=(look or Style()).helices_310, bridges=dssp.bridges, matched=matched)
    sheets = build_sheets(sses, dssp.bridges)
    nucleic = load_nucleic(path, bb, assembly, model=model) if dna else None
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
    mem = None
    if membrane not in ("off", None):
        from .membrane import find_membrane

        mem = find_membrane(path, bb)
        flip = membrane == "flip"
        if mem is not None and flip:  # the other face is the outside: a 180° turn, never a mirror
            mem = dataclasses.replace(mem, normal=-mem.normal, centre=-mem.centre)
        if mem is None:
            print(f"foldmap: warning: no membrane found in {path}; drawing it without one", file=sys.stderr)
        elif up is None and view is None:  # the bilayer horizontal, outside on top (unless the symmetry axis
            along = sym is not None and abs(float(np.dot(sym.axis, mem.normal))) > 0.9  # already stands upright)
            if not along:
                up = mem.normal
            elif (mem.sided and float(np.dot(sym.axis, mem.normal)) < 0) or (flip and not mem.sided):  # turn it over
                sym = dataclasses.replace(
                    sym, axis=-sym.axis, angles=[-a for a in sym.angles], rise=-sym.rise, twist=-sym.twist
                )
    duplex = nucleic is not None and straight_duplex(nucleic.structure, nucleic.strands)
    dyad = None  # a two-fold assembly seen from the side, drawn as one copy and its mirror image (plan 18)
    rosette_asked = (  # Style.cyclic rosette: a cyclic assembly of ROSETTE_ORDERS copies seen down its axis, a D2 as
        (look or Style()).cyclic == "rosette"  # two such rings mirrored in a dyad (plan 21 Task 3)
        and sym is not None
        and ((sym.kind == "C" and sym.n in ROSETTE_ORDERS) or (sym.kind == "D" and sym.n == 2))
    )
    cyclic_note = ""
    if (look or Style()).cyclic == "rosette" and sym is not None and sym.kind == "D" and sym.n > 2:
        cyclic_note = (
            f"cyclic=rosette: the {sym.label} assembly drawn unrolled, only a D2 is drawn as a rosette (seen down its "
            "main axis, the rings' interface lies in depth)"
        )
    if rosette_asked and (duplex or rotate or flip_v or up is not None or view is not None):
        cyclic_note = f"cyclic=rosette: the {sym.label} assembly not drawn as a rosette, " + (
            "it sits on a straight duplex" if duplex else "a view, a turn or a flip was asked for"
        )
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
        if rosette_asked:  # cyclic=rosette (plan 21): seen down the axis, each copy protomer 0 turned, unless
            why = _no_rosette(sses, sym, nucleic, mode, mem)  # what it cannot draw yet
            rings = None
            if not why and sym.kind == "D":  # a D2: the two-fold that takes one half of each protomer onto the other
                found = protomer_dyad(bb, sym, symmetry_tol)
                rings = dihedral_rings(sym, found[1]) if found is not None else None
                if rings is None:
                    why = "no two-fold across the main axis takes half of each protomer onto the other half"
            if why:
                cyclic_note = f"cyclic=rosette: the {sym.label} assembly drawn unrolled, {why}"
            else:
                if rings is not None:
                    frame = dihedral_rosette_frame(sses, sym, found[0], rings)
                    ring_a, ring_b = (", ".join(c for copy in ring for c in copy) for ring in rings)
                    cyclic_note = (
                        f"D2 assembly seen down its main two-fold axis (cyclic=rosette): ring A (chains {ring_a}) "
                        f"drawn as one copy and that copy turned by 180°, ring B (chains {ring_b}) below the dashed "
                        "dyad line as ring A's mirror image (the half-turn about the dyad, seen down the same axis)"
                    )
                else:
                    frame = rosette_frame(sses, sym)
                    step = 360 // sym.n if 360 % sym.n == 0 else round(360 / sym.n, 1)
                    cyclic_note = (
                        f"{sym.label} assembly seen down its {sym.n}-fold axis (cyclic=rosette): each copy drawn as "
                        f"the first turned by a multiple of {step}°"
                    )
                sheets = build_sheets(sses, _within_protomer(dssp.bridges, frame.sector, bb))
                across = _pairings_across(dssp.bridges, frame.sector, bb, sses)
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
        stems = [StemView(_snap(v.centroid), _snap(v.axis), v.n) for v in stems]
        frame = view_frame(sses + stems, rotate, flip_v, up, view)
    if sym is not None and sym.kind == "C" and sym.n == 2 and isinstance(frame, (RadialFrame, Frame)):
        _mirror_frame(frame, sym)
    if duplex and sym is not None and not getattr(frame, "sector", None):
        _side_by_side(frame, sym, nucleic, bb)
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
    named = _domains(domains, sses, contacts, rows, copies, partial, bb)
    named_note = domain_note(source, named, found_note)
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
            alike_first=mem is None,  # a membrane view is judged by its bands (C17), not its copies (A15)
        )

    def finish(lay, frame, faces=None):  # what make_layout sets on a layout after build_layout through `frame`
        lay.helix_matched = {"cut": matched.get("cut", 0), "joined": matched.get("joined", 0)}
        lay.helix_matched_ids = list(matched.get("ids", []))
        lay.domain_source, lay.domain_note = (source if lay.domains else ""), named_note
        if lay.domains:  # copies of one domain share a tone
            from .features import domain_tones

            lay.domain_tones = domain_tones(lay.domains, sses, copies)
        lay.focus = _focus((look or Style()).highlight, bb, sym)
        if focus_chains:
            lay.focus, lay.focus_chains, lay.context = set(focus_chains), set(focus_chains), context
            lay.unjoined, lay.subunits = unjoined, subunits
        lay.chains_read = read
        lay.chains_asked = [c for c in chains or [] if c in present]
        lay.na_read = na_read if na_cropped else 0
        lay.links = links
        lay.res_chain = [l.chain for l in bb.labels]
        lay.res_name = [l.name for l in bb.labels]
        lay.res_labels = list(bb.labels)
        lay.res_linked = sorted(bb.linked)
        lay.res_b = [float(x) for x in bb.b] if bb.b is not None else []
        with_elements = {s.chain for s in sses}
        lay.coil_chains = [(c, n) for c, n in Counter(lay.res_chain).items() if c not in with_elements]
        lay.uniprot = ann
        lay.cyclic_note = cyclic_note
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
        lay = finish(lay, frame)
        if lay.let_go:  # copy items let go (layout._one_per_mate): kept only when the copies come out more alike (A15)
            held = finish(lay_out(let_go=False), frame)
            if _unalike(held, bb) <= _unalike(lay, bb):
                lay, lay_out = held, functools.partial(lay_out, let_go=False)
        if mem is not None and mode == "projected" and lay.membrane is not None:
            from .figcheck import membrane_view

            if membrane_view(lay, sses, bb, mem)["tm_short"]:  # a transmembrane helix drawn off its band: laid out
                before = _membrane_score(lay, sses, bb, mem)  # again with every one at its projected height (C17)
                level = _transmembrane(sses, bb, mem)
                faces = _faces_on_page(frame, sses, bb, mem, level)
                again = finish(lay_out(level=level, band=tuple(sorted(faces))), frame, faces)
                after = _membrane_score(again, sses, bb, mem)
                # kept only when it helps: fewer helices off the band, no more loops in the slab, no more fallbacks
                if after[0] < before[0] and after[1] <= before[1] and after[2] <= before[2]:
                    lay = again
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
    lay.cyclic_assembly = (
        sym is not None
        and ((sym.kind == "C" and sym.n in ROSETTE_ORDERS) or (sym.kind == "D" and sym.n == 2))
        and not duplex
    )
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


def _membrane_score(lay, sses, bb, mem) -> tuple[int, int, int]:
    """A membrane view judged as A18 judges it (figcheck.membrane_view, on a routed copy of the layout): transmembrane
    helices off their band, loops running in the slab, and route fallbacks; make_layout keeps a second layout only
    when it does better (C17, fix round 1)."""
    import copy

    from .figcheck import membrane_view

    trial = copy.deepcopy(lay)
    loops = route_loops(trial, sses, bb)
    view = membrane_view(trial, sses, bb, mem, loops)
    return len(view["tm_short"]), len(view["loops_in_slab"]), sum(lp.fallback for lp in loops)


def _transmembrane(sses, bb, mem) -> set[str]:
    """The elements that cross the membrane (their CA run reaching both faces, as _face_medians counts them) and
    define it (Membrane.defining): in a projected membrane view they keep their projected height (build_layout
    `level`), so each crosses the band where the membrane is (random-100 C17)."""
    depth = mem.depth(bb.ca)
    out = set()
    for s in sses:
        seg = depth[s.start : s.end + 1]
        if len(seg) < 2 or seg.min() > -mem.half + FACE_REACH or seg.max() < mem.half - FACE_REACH:
            continue
        if mem.defining is not None and mem.defining.isdisjoint(range(s.start, s.end + 1)):
            continue
        out.add(s.id)
    return out


def _faces_on_page(frame, sses, bb, mem, tm: set[str]) -> tuple[float, float]:
    """Page y (layout units, snapped) of the membrane's outer and inner face, projected through the frame at the
    transmembrane elements' centre (build_layout `band`, _membrane_band `faces`)."""
    first = next(s for s in sses if s.id in tm)
    at = np.mean([bb.ca[s.start : s.end + 1].mean(axis=0) for s in sses if s.id in tm], axis=0)
    mid = at - mem.normal * float(np.ravel(mem.depth(at))[0])
    y_out, y_in = (
        float(frame.project(mid + mem.normal * h, [first.chain])[0, 1]) * SCALE for h in (mem.half, -mem.half)
    )
    return float(_snap(y_out)), float(_snap(y_in))


def _face_medians(lay, sses, bb, mem, reach: bool = False) -> tuple[float, float] | None:
    """Where these elements pass the outer and the inner membrane face on the page (medians), if any is TM and
    defines the membrane (Membrane.defining: a beta-barrel's strands, not a subunit lying in its lumen). When none
    passes both faces, the TM elements that only reach a face (within membrane.FACE_REACH, as they are counted TM)
    give it at their CA nearest it (`reach`: random-100 held-out 6MJ2, whose S1-S4 helices stop 1-2 Å short of the
    inner face)."""
    depth = mem.depth(bb.ca)
    outs, ins = [], []
    for s in sses:
        p = lay.placed.get(s.id)
        seg = depth[s.start : s.end + 1]
        if p is None or len(seg) < 2 or seg.min() > -mem.half + FACE_REACH or seg.max() < mem.half - FACE_REACH:
            continue
        if mem.defining is not None and mem.defining.isdisjoint(range(s.start, s.end + 1)):
            continue

        def page_y(face: float, seg=seg, p=p) -> float | None:
            k = np.flatnonzero(np.diff(np.sign(seg - face)) != 0)
            if not k.size:
                if not reach:
                    return None
                k = int(np.argmin(np.abs(seg - face)))  # the run's CA nearest the face it reaches
                return float(p.n_port[1] + k / (len(seg) - 1) * (p.c_port[1] - p.n_port[1]))
            k = int(k[0])
            t = (k + (face - seg[k]) / (seg[k + 1] - seg[k])) / (len(seg) - 1)
            return float(p.n_port[1] + t * (p.c_port[1] - p.n_port[1]))

        o, i = page_y(mem.half), page_y(-mem.half)
        if o is not None and i is not None and abs(o - i) > 1.0:
            outs.append(o)
            ins.append(i)
    if not outs:
        return None if reach else _face_medians(lay, sses, bb, mem, reach=True)
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


def domain_note(source: str, named, found: str = "") -> str:
    """What the notes say of the domains (Layout.domain_note): automatic ones are inferred, not published (plan 12
    R2); CATH's come with their ids (`found`, cath.domains' note), or the note says why there are none."""
    if source == "CATH":
        return found
    if source != "inferred":
        return ""
    stems = list(dict.fromkeys(name.split(" · ")[0] for name, _ in named))
    if not stems:
        return "no domains inferred: each chain reads as one domain from its residue contacts"
    return f"domains {', '.join(stems)} inferred from residue contacts, not a published domain assignment"


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


def _domains(spec, sses, contacts, rows=(), copies=None, lent=frozenset(), bb=None) -> list[tuple[str, list[str]]]:
    """[(name, element ids)] from 'auto' (chain by chain, copies of a protein (`copies`) alike, chains that only
    lend elements (`lent`) in their mates' domains, never splitting a βα row of `rows`) or [(name, [refs])]; every
    element in at most one domain."""
    from .features import auto_domains
    from .layout import provisional

    if not spec:
        return []
    if spec == "auto":
        return auto_domains(sses, contacts, rows, copies, lent, bb=bb)
    if _element_ids(spec):
        known = {s.id for s in sses}
        stale = stale_domain_ids(spec, sses)
        if stale and not isinstance(spec, SavedDomains):
            name, i = stale[0]
            raise ValueError(
                f"domain {name!r}: {i!r} is no element; write a residue range as res:{i}, "
                "or name elements by label (α1, A:G)"
            )
        seen: set[str] = set()  # each element in one domain: the first that names it
        kept = []
        for name, ids in spec:
            kept.append((name, [i for i in dict.fromkeys(ids) if i in known and i not in seen]))
            seen.update(kept[-1][1])
        if stale and not getattr(spec, "quiet", False):
            print(f"foldmap: warning: {stale_domains_note(spec, sses)}", file=sys.stderr)
        return [(name, ids) for name, ids in kept if ids]  # already element ids (UniProt domains, a layout file)
    lay, out, seen = provisional(sses), [], set()
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


def _adjustments(args) -> dict:
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


def _spec_of(args) -> FigureSpec:
    """The figure the command line asks for: a layout file (if given) overlaid by the options given."""
    shorthand = [
        f"{k}={v}" for k, v in (("palette", args.palette), ("fill", args.style), ("loops", args.loops)) if v is not None
    ]
    edits = _adjustments(args)
    return figure_spec(
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
        domains=_domain_spec(args.domain) or args.domains,
        msa=args.msa,
        msa_reference=args.msa_reference,
        membrane=args.membrane,
        uniprot=args.uniprot,
        chains=args.chains,
        focus=args.focus,
        rename=edits["rename"],
        swap=edits["swap"],
        move=edits["move"],
        names=args.names,
    )


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
        if key in CHOICES:
            return "one of " + ", ".join(CHOICES[key])
        if key in SCALES:
            return "positive number, 1 = default size"
        if key == "na_backbone":
            return "modality, chain, ink or a colour"
        if key in ("dna_color", "rna_color", "background", "ink"):
            return "a colour (#rrggbb or a name)"
        if key == "hbonds":
            return "off, na, interchain or chains:A+B,C+D"
        if key == "hbond_color":
            return "a colour (#rrggbb or a name)"
        if key == "hbond_max":
            return "whole number, at most this many lines (the shortest on the page kept, not the shortest in Å)"
        if key == "backbone_colors":
            return "KEY=COLOUR,... (KEY: a chain id, protein, dna or rna)"
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
        "reason. Prints one line per check and one per problem; exits 1 on a problem, 0 otherwise.",
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
            help="your own element name, e.g. res:167-182=Gd or A:B=B' (repeatable)",
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
            "--msa",
            metavar="ALIGNMENT",
            help="alignment (or a ConSurf grades table) for conservation colouring (--theme conservation)",
        )
        q.add_argument("--msa-reference", metavar="NAME", help="alignment row that is the structure (default: best)")
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
    return parser


def _resolve_structure(text: str) -> str:
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

    fig = draw_sequence(
        args.structure,
        chain=args.chain,
        alignment=read_alignment(args.msa) if args.msa else None,
        reference=args.reference,
        columns=args.columns,
        look=resolve_style(args.theme or "publication", args.style_file, args.set),
        full_sequence=args.full_sequence,
        ss_colour=args.ss_colour,
        threshold=args.similarity,
        conservation_bar=not args.no_conservation_bar,
        title=args.title,
        uniprot=args.uniprot,
        msa=args.msa,
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
    Exit 1 when any check found a problem or any criterion judged here fails."""
    from .criteria import check_figure, criteria_lines
    from .science import science_lines

    spec = _spec_of(args)
    drawn = draw_spec(spec)
    print(provenance_block(drawn.provenance), file=sys.stderr, end="")
    for out in args.output:
        print(f"wrote {save(drawn.figure, out)}")
    found = check_figure(drawn, spec, args.structure)
    for line in science_lines(found["science"], found["science_judged"]) + criteria_lines(found["criteria"]):
        print(line)
    failed = any(c["pass"] is False for c in found["criteria"].values())
    return 1 if any(found["science"].values()) or failed else 0


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
}


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if getattr(args, "structure", None):
            args.structure = _resolve_structure(args.structure)
        return _COMMANDS[args.command](args)
    except (FileNotFoundError, ValueError, RuntimeError) as err:  # RuntimeError: gemmi cannot read the file
        print(f"foldmap: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
