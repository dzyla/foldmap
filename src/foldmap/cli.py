# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Command line: `foldmap summary FILE` prints secondary structure and sheet ordering."""

from __future__ import annotations

import argparse
import dataclasses
import functools
import re
import sys
from collections.abc import Callable
from pathlib import Path

import numpy as np

from . import layoutfile
from .dssp import assign_dssp
from .features import helices_on_sheets, sse_contacts
from .frame import Frame, RadialFrame, dna_radial_frame, dna_view_frame, protomer_maps, symmetric_frame, view_frame
from .io import load_backbone, load_links, load_nucleic, read_model
from .layout import SCALE, _snap, build_layout, na_torn_by_unrolling, resolve, template_rows
from .model import AUTO_NT, crop_nucleic, homolog_chains, nucleic_of_chains, nucleic_touching, res_ref
from .names import apply_names, check_spec
from .nucleic import StemView, stem_views, straight_duplex
from .palette import PALETTES
from .render import LOOP_STYLES, STYLES, draw, save
from .route import route_loops
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
from .symmetry import TOLERANCE, detect_symmetry

_vector, _protomers, _domain_spec = parse_vector, parse_protomers, parse_domains  # older names, used by the app


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
    for n, sh in enumerate(sheets, 1):
        flags = ("closed " if sh.closed else "") + ("ambiguous" if sh.ambiguous else "")
        lines.append(f"sheet {n}: {len(sh.strands)} strand(s) {flags}".rstrip())
        parts = []
        for k, (s, d) in enumerate(zip(sh.strands, sh.directions)):
            parts.append(f"{s.id}({'+' if d > 0 else '-'})")
            if k < len(sh.pair_kinds):
                parts.append({"P": "=par=", "A": "=anti=", "?": "=?="}[sh.pair_kinds[k]])
        lines.append("  " + " ".join(parts))
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
    """Whether a nucleic acid with stems has a strand held by two or more protomers of `sym`, each touching it
    with at least NA_HOLD_MIN residues: the precondition for drawing it in a projected frame (unrolling the assembly
    one protomer beside the next may tear that strand apart; layout.na_torn_by_unrolling decides)."""
    if nucleic is None or sym is None or not nucleic.structure.stems:
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
    sses = build_sses(bb, dssp.ss, short_helices=(look or Style()).helices_310, bridges=dssp.bridges)
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
    focus_chains = _focus_chains(focus, sym, bb, sses)
    if focus_chains and nucleic is not None and nucleic.strands:  # one subunit: the nucleic acid it holds
        nucleic = nucleic_touching(nucleic, focus_chains)
    na_cropped = False
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
    if nucleic is not None and straight_duplex(nucleic.structure, nucleic.strands):
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
        frame = symmetric_frame(sses, sym)
    else:  # stem axes differ in the last bit between BLAS kernels: snapped before they choose the page axes
        stems = stem_views(nucleic.structure, nucleic.strands) if nucleic is not None else []
        stems = [StemView(_snap(v.centroid), _snap(v.axis), v.n) for v in stems]
        frame = view_frame(sses + stems, rotate, flip_v, up, view)
    if sym is not None and sym.kind == "C" and sym.n == 2 and isinstance(frame, (RadialFrame, Frame)):
        _mirror_frame(frame, sym)
    links = load_links(path, bb, assembly, model=model)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, links).items():
        contacts[key] = contacts.get(key, 0) + n
    ann = None
    if uniprot not in (None, "off", "none"):
        from .uniprot import annotate

        ann = annotate(path, bb, None if uniprot == "auto" else uniprot)
    if domains == "uniprot":
        if ann is None:
            raise ValueError("--domains uniprot needs UniProt annotation: add --uniprot auto (or an accession)")
        domains = _uniprot_domains(ann, sses)
    rows = (  # βα rows are rigid units: automatic domains never split one
        template_rows(sheets, sses, copies) if mode == "projected" and (look or Style()).fold_template == "auto" else []
    )
    named = _domains(domains, sses, contacts, rows, copies, partial)
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
    lay = build_layout(
        sses,
        sheets,
        frame,
        mode,
        contacts=contacts,
        nucleic=nucleic,
        style=look,
        bridges=_bridge_springs(bb, sses, links),
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
    )
    if lay.domains:  # copies of one domain share a tone
        from .features import domain_tones

        lay.domain_tones = domain_tones(lay.domains, sses, copies)
    lay.focus = _focus((look or Style()).highlight, bb, sym)
    if focus_chains:
        lay.focus, lay.focus_chains, lay.context = set(focus_chains), set(focus_chains), context
        lay.unjoined = unjoined
    lay.links = links
    lay.res_chain = [l.chain for l in bb.labels]
    lay.res_name = [l.name for l in bb.labels]
    lay.res_labels = list(bb.labels)
    lay.res_linked = sorted(bb.linked)
    lay.res_b = [float(x) for x in bb.b] if bb.b is not None else []
    lay.uniprot = ann
    first = [s.centroid for s in sses if s.chain == sses[0].chain] if sses else [_any_centre(bb, nucleic)]
    centre = np.mean(first, axis=0)
    lay.page_axes = [list(map(float, a)) for a in frame.basis(centre)]
    lay.page_centre = list(map(float, bb.ca.mean(axis=0) if len(bb) else centre))
    if mem is not None:
        lay.membrane = _membrane_band(lay, sses, bb, mem)
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


def _membrane_band(lay, sses, bb, mem) -> dict | None:
    """The band on the page: where the transmembrane elements pass the two membrane faces (medians). A
    depth-keeping stack wrapped into rows has one band per row ("bands"), each from its own row's elements; a
    row with no transmembrane element (a soluble partner chain) has none."""
    rows = [set(ids) for ids in lay.stack_rows] or [set(lay.placed)]
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


def _face_medians(lay, sses, bb, mem) -> tuple[float, float] | None:
    """Where these elements pass the outer and the inner membrane face on the page (medians), if any is TM."""
    depth = mem.depth(bb.ca)
    outs, ins = [], []
    for s in sses:
        p = lay.placed.get(s.id)
        seg = depth[s.start : s.end + 1]
        if p is None or len(seg) < 2 or seg.min() > -mem.half + 2 or seg.max() < mem.half - 2:
            continue

        def page_y(face: float, seg=seg, p=p) -> float | None:
            k = np.flatnonzero(np.diff(np.sign(seg - face)) != 0)
            if not k.size:
                return None
            k = int(k[0])
            t = (k + (face - seg[k]) / (seg[k + 1] - seg[k])) / (len(seg) - 1)
            return float(p.n_port[1] + t * (p.c_port[1] - p.n_port[1]))

        o, i = page_y(mem.half), page_y(-mem.half)
        if o is not None and i is not None and abs(o - i) > 1.0:
            outs.append(o)
            ins.append(i)
    if not outs:
        return None
    return float(np.median(outs)), float(np.median(ins))


FOCUS_CHAINS, FOCUS_ELEMENTS = 6, 160  # symmetric assemblies larger than this are drawn one subunit at a time


def _focus_chains(focus: str, sym, bb, sses) -> set | None:
    """The chains drawn in full: none ('none', or a modest assembly under 'auto'), one protomer (helical
    filaments, cubic cages, or assemblies over FOCUS_CHAINS chains / FOCUS_ELEMENTS elements), or the chains
    named. Helical filaments use their middle subunit, so it has neighbours on both sides."""
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
    return _one_of_each(sym.protomers[k], bb)


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
    backbone H-bond bridges (donor strands, domain swaps, beta-augmentation) and helices bundled with theirs."""
    from .features import helix_bundles

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
    for group in helix_bundles(sses):
        if own & set(group):
            guests |= {k for k in group if k not in own}
    keep = [s for s in sses if s.id in own or s.id in guests]
    return keep, {s.chain for s in keep if s.chain not in focus}


def _homolog_labels(bb, sses) -> dict[str, str]:
    """Labels shared by copies of a protein, however much of each is modelled (model.shared_labels): one element,
    one label, in every copy."""
    from .layout import provisional
    from .model import homolog_chains, shared_labels

    labels = {k: p.label for k, p in provisional(sses).placed.items()}
    return shared_labels(sses, homolog_chains(bb), labels)


def _context(sym, focus: set, bb, neighbours: bool = True) -> str:
    """One line saying what the figure shows of the assembly."""
    total = len(set(l.chain for l in bb.labels))
    grey = " · grey: neighbouring subunits completing its fold" if neighbours else ""
    if sym is None:
        return f"chain {', '.join(sorted(focus))} of {total} shown in full" + grey
    copies = sum(len(p) for p in sym.protomers) / max(len(sym.protomers), 1) / max(len(focus), 1)
    n = int(round(len(sym.protomers) * copies))
    if sym.kind == "H":
        step = "subunit" if copies <= 1 else f"repeat of {int(round(copies * len(focus)))} chains"
        twist = f"{np.degrees(abs(sym.twist)) * (-1 if sym.left_handed else 1):.1f}".replace("-", "−")
        kind = f"helical, {twist}° twist{sym.hand}, {abs(sym.rise):.1f} Å rise per {step}"
    else:
        kind = f"{sym.label} symmetry"
    return f"1 of {n} subunits shown · {kind}" + grey


_DOMAIN_PULL = 6  # CA-pair-equivalents of attraction between any two elements of one domain


UNIPROT_DOMAINS = ("Domain", "Repeat", "Zinc finger", "DNA binding")


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
    if spec in (None, "auto", "uniprot") or not _element_ids(spec):
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


def _domains(spec, sses, contacts, rows=(), copies=None, lent=frozenset()) -> list[tuple[str, list[str]]]:
    """[(name, element ids)] from 'auto' (chain by chain, copies of a protein (`copies`) alike, chains that only
    lend elements (`lent`) in their mates' domains, never splitting a βα row of `rows`) or [(name, [refs])]; every
    element in at most one domain."""
    from .features import auto_domains
    from .layout import provisional

    if not spec:
        return []
    if spec == "auto":
        return auto_domains(sses, contacts, rows, copies, lent)
    if _element_ids(spec):
        known = {s.id for s in sses}
        stale = stale_domain_ids(spec, sses)
        if stale and not isinstance(spec, SavedDomains):
            name, i = stale[0]
            raise ValueError(
                f"domain {name!r}: {i!r} is no element; write a residue range as res:{i}, "
                "or name elements by label (α1, A:G)"
            )
        kept = [(name, [i for i in ids if i in known]) for name, ids in spec]
        if stale:
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
    for cmd in (q, iq, cq):
        q = cmd
        q.add_argument("structure")
        q.add_argument(
            "-o",
            "--output",
            action="append",
            required=cmd is not cq,
            default=[],
            help=".svg, .pdf or .png (repeatable)"
            + ("; optional here: also save the figure checked" if cmd is cq else ""),
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
            choices=["auto", "uniprot"],
            help="auto: from element contacts (D1, D2...); uniprot: UniProt's domains (with --uniprot)",
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
            "(auto: helical filaments, cages and assemblies over 6 chains / 160 elements)",
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

        for out in args.output:
            print(f"wrote {write_page(spec.structure, out, viewer=args.viewer, **spec.page_kwargs())}")
    else:
        drawn = draw_spec(spec)
        if drawn.stuck:
            print(
                f"foldmap: warning: no clear route for {len(drawn.stuck)} loop(s) or strand(s), drawn as plain "
                f"curves that may cross elements: {', '.join(drawn.stuck)}",
                file=sys.stderr,
            )
        if drawn.over_labels:
            print(
                f"foldmap: warning: {len(drawn.over_labels)} loop(s) drawn over a label or chain end, the only route "
                f"found: {', '.join(drawn.over_labels)}",
                file=sys.stderr,
            )
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
        print(f"wrote {layoutfile.save(args.save_layout, layout_document(spec))}")
    return 0


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
