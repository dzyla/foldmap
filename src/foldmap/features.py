# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Facts derived from the structure that the layout and annotation layers need."""

from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .model import SSE, Backbone, Sheet, copy_classes, group_by_chain

CONTACT_CUTOFF = 8.0  # Å between CA atoms
MIN_CONTACT_PAIRS = 3  # CA pairs needed to call two elements neighbours


def sse_contacts(bb: Backbone, sses: list[SSE], cutoff: float = CONTACT_CUTOFF) -> dict[frozenset[str], int]:
    """Number of CA-CA pairs within `cutoff` between each pair of different elements."""
    owner = np.full(len(bb), -1)
    for k, s in enumerate(sses):
        owner[s.start : s.end + 1] = k
    idx = np.flatnonzero(owner >= 0)
    if len(idx) < 2:
        return {}
    pairs = cKDTree(bb.ca[idx]).query_pairs(cutoff, output_type="ndarray")
    counts: Counter[frozenset[str]] = Counter()
    for a, b in pairs:
        ka, kb = owner[idx[a]], owner[idx[b]]
        if ka != kb:
            counts[frozenset((sses[ka].id, sses[kb].id))] += 1
    return {k: v for k, v in counts.items() if v >= MIN_CONTACT_PAIRS}


BUNDLE_MIN = 10  # residues: shorter helices are turns at an interface, not a coiled-coil stalk
BUNDLE_COS = 0.75  # |cos| between axes: parallel or antiparallel
BUNDLE_DIST = 12.0  # Å between axis segments
BUNDLE_COVER = 0.8  # fraction of the shorter helix that must run that close
KINK_GAP = 2  # residues: helix pieces this close in one chain, roughly collinear, count as one helix
KINK_COS = 0.5  # a kink bends the axis by up to 60 degrees
KINK_SHORT = 7  # residues: a piece this short has no reliable axis of its own, so it continues its neighbour


def _segment(s: SSE, n: int = 20) -> np.ndarray:
    half = (len(s) - 1) * 0.75  # 1.5 Å rise per residue
    return np.linspace(s.centroid - s.axis * half, s.centroid + s.axis * half, n)


def _shape(pieces: list[SSE]) -> tuple[np.ndarray, int, np.ndarray]:
    """Points along a helix (one piece, or the kinked pieces of one), its length and its axis (longest piece).
    Kinked pieces get a point per residue, so each weighs by its length."""
    if len(pieces) == 1:
        return _segment(pieces[0]), len(pieces[0]), pieces[0].axis
    pts = np.vstack([_segment(s, max(2, len(s))) for s in pieces])
    return pts, sum(len(s) for s in pieces), max(pieces, key=len).axis


def helix_bundles(sses: list[SSE]) -> list[list[str]]:
    """Helices of different chains packed along their length (coiled-coil stalks, zippers, interface pairs):
    groups of element ids, each spanning at least two chains, sorted. Kinked pieces of one helix stay together."""
    helices = [s for s in sses if s.kind == "H"]
    group = {s.id: s.id for s in helices}

    def root(x: str) -> str:
        while group[x] != x:
            x = group[x]
        return x

    def join(a: str, b: str) -> None:
        group[root(a)] = root(b)

    kinked: dict[str, list[SSE]] = {}
    nexts = list(zip(helices, helices[1:]))
    cos = [float(a.axis @ b.axis) for a, b in nexts]
    kink = [  # pieces k and k+1 are one helix
        a.chain == b.chain and b.start - a.end - 1 <= KINK_GAP and (x > KINK_COS or min(len(a), len(b)) < KINK_SHORT)
        for (a, b), x in zip(nexts, cos)
    ]
    for k in range(1, len(helices) - 1):  # a short piece joins the pieces either side only when they are collinear
        a, b, c = helices[k - 1 : k + 2]
        if len(b) < KINK_SHORT and kink[k - 1] and kink[k] and float(a.axis @ c.axis) <= KINK_COS:
            weaker = k - 1 if cos[k - 1] < cos[k] else k
            if cos[weaker] <= KINK_COS:  # it continues the neighbour it is closer in direction to
                kink[weaker] = False
    for k in np.flatnonzero(kink):
        join(helices[k].id, helices[k + 1].id)
    for s in helices:
        kinked.setdefault(root(s.id), []).append(s)
    shapes = [[s] for s in helices] + [p for p in kinked.values() if len(p) > 1]  # pieces, and whole kinked helices
    long_ = [(p, *_shape(p)) for p in shapes if sum(len(s) for s in p) >= BUNDLE_MIN]
    linked = set()
    for i, (a, pa, la, xa) in enumerate(long_):
        for b, pb, lb, xb in long_[i + 1 :]:
            if a[0].chain == b[0].chain or abs(float(xa @ xb)) < BUNDLE_COS:
                continue
            d = np.linalg.norm(pa[:, None] - pb[None], axis=2)
            near = d.min(1) if la <= lb else d.min(0)
            if (near < BUNDLE_DIST).mean() >= BUNDLE_COVER:
                join(a[0].id, b[0].id)
                linked |= {a[0].id, b[0].id}
    out: dict[str, list[str]] = {}
    for s in helices:
        out.setdefault(root(s.id), []).append(s.id)
    chain = {s.id: s.chain for s in helices}
    return sorted(sorted(g) for g in out.values() if len({chain[k] for k in g}) > 1 and linked & set(g))


ON_SHEET_MIN = 10  # residues: a shorter helix is a turn at the sheet's edge, not one packed on its face
ON_SHEET_STRANDS = 3  # strands: a smaller sheet has no face to pack on
ON_SHEET_REACH = 10.0  # Å between a helix CA and the nearest sheet CA, counting it as touching the sheet
ON_SHEET_TOUCH = 0.3  # fraction of the helix's CAs that must touch the sheet
ON_SHEET_FLAT = 0.5  # |cos| between the helix axis and the sheet normal: the helix lies along the face
ON_SHEET_HEIGHT = (6.0, 14.0)  # Å from the sheet plane to the helix centroid: one helix packed on the face
ON_SHEET_EDGE = 2.0  # Å: the centroid, projected on the sheet plane, lies within the sheet's extent plus this


def helices_on_sheets(bb: Backbone, sses: list[SSE], sheets: list[Sheet]) -> dict[str, int]:
    """Helices packed on a sheet face: helix id -> index of the sheet in `sheets`. A helix of >= ON_SHEET_MIN residues
    packs on a sheet of >= ON_SHEET_STRANDS strands, one of them of the helix's chain, when >= ON_SHEET_TOUCH of its
    CAs have a sheet CA within ON_SHEET_REACH, its axis lies along the face (|axis . n| < ON_SHEET_FLAT, n the sheet
    normal: the smallest singular vector of the sheet's CAs), its centroid sits ON_SHEET_HEIGHT from the sheet plane
    and, projected on that plane, within the sheet CAs' in-plane extent + ON_SHEET_EDGE. A helix on several sheets
    goes to the one most of its CAs touch (the first on a tie). Coordinates, axes and the plane are snapped before
    any threshold sees them (layout._snap). Plan 11 measured: 1UBQ α1, 2ACY A:22-32 and A:55-67, 2LZM A:39-50."""
    from .layout import _snap as _grid  # layout imports this module: imported when called

    ca = _grid(bb.ca)
    faces = []
    for k, sh in enumerate(sheets):
        if len(sh.strands) < ON_SHEET_STRANDS:
            continue
        pts = np.vstack([ca[s.start : s.end + 1] for s in sh.strands])
        c = _grid(pts.mean(axis=0))
        _, _, vt = np.linalg.svd(pts - c)
        e1, e2, n = _grid(vt[0]), _grid(vt[1]), _grid(vt[2])
        flat = (pts - c) @ np.array([e1, e2]).T
        faces.append((k, {s.chain for s in sh.strands}, cKDTree(pts), c, n, e1, e2, flat.min(0), flat.max(0)))
    out = {}
    for h in sses:
        if h.kind != "H" or len(h) < ON_SHEET_MIN:
            continue
        own, axis, centre = ca[h.start : h.end + 1], _grid(h.axis), _grid(h.centroid)
        best = None
        for k, chains, tree, c, n, e1, e2, lo, hi in faces:
            if h.chain not in chains:
                continue
            touch = sum(bool(hits) for hits in tree.query_ball_point(own, ON_SHEET_REACH))
            at = _grid([(centre - c) @ e1, (centre - c) @ e2])
            if (
                touch >= ON_SHEET_TOUCH * len(h)
                and abs(float(axis @ n)) < ON_SHEET_FLAT
                and ON_SHEET_HEIGHT[0] < abs(float(_grid((centre - c) @ n))) < ON_SHEET_HEIGHT[1]
                and np.all(at >= lo - ON_SHEET_EDGE)
                and np.all(at <= hi + ON_SHEET_EDGE)
                and (best is None or touch > best[1])
            ):
                best = (k, touch)
        if best is not None:
            out[h.id] = best[0]
    return out


DOMAIN_MIN_SSE_RESIDUES = 30  # residues in helices/strands an automatic domain needs


def auto_domains(
    sses: list[SSE],
    contacts: dict[frozenset[str], int],
    together: list[list[str]] = (),
    copies: dict[str, tuple[str, int]] | None = None,
    lent: set[str] = frozenset(),
) -> list[tuple[str, list[str]]]:
    """Domains, chain by chain: each chain's split (_chain_domains), named D1, D2... in sequence order (with
    " · chain" when several chains are split). Copies of one protein (`copies`, model.homolog_chains) share their
    reference copy's split (model.copy_classes): each element joins the domain of its mate there, an element the
    reference lacks the domain of the element before it (else after it) in its chain. Chains lending a few
    elements to the drawn subunit (`lent`, focus mode) add them to their mates' domains instead of having their own.
    [] when no chain reads as more than one domain."""
    classes = copy_classes(sses, copies or {})
    mine = group_by_chain(sses)
    root = {c: (copies or {}).get(c, (c, 0))[0] for c in mine}
    owner: dict[str, tuple[str, int]] = {}  # element -> (chain, its domain's number there)
    for chain, els in mine.items():
        if classes.reference.get(root[chain], chain) == chain:
            for d, ids in enumerate(_chain_domains(els, contacts, together)):
                owner.update({k: (chain, d) for k in ids})
    for chain, els in mine.items():
        if classes.reference.get(root[chain], chain) == chain:
            continue
        got = [owner.get(classes.mate[s.id]) for s in els]
        if not any(got):
            continue
        for k, s in enumerate(els):
            near = got[k] or next((g for g in reversed(got[:k]) if g), None) or next(g for g in got[k:] if g)
            owner[s.id] = near if chain in lent else (chain, near[1])
    groups: dict[tuple[str, int], list[str]] = {}
    for s in sses:
        if s.id in owner:
            groups.setdefault(owner[s.id], []).append(s.id)
    split = list(dict.fromkeys(chain for chain, _ in groups))
    out = []
    for chain in split:
        for d in sorted(d for c, d in groups if c == chain):  # a copy's domains are numbered as its reference's
            out.append((f"D{d + 1}" if len(split) == 1 else f"D{d + 1} · {chain}", groups[(chain, d)]))
    return out


def domain_tones(
    domains: list[tuple[str, list[str]]], sses: list[SSE], copies: dict[str, tuple[str, int]]
) -> list[int]:
    """The tone (index into the domain tones) of each domain: copies of one domain in copies of a protein share one,
    distinct domains take the next, in domain order. A domain is a copy of the one holding most of its elements'
    mates (model.copy_classes); a domain whose mates lie in no domain, or in itself, is its own."""
    mate = copy_classes(sses, copies).mate
    owner = {k: d for d, (_, ids) in enumerate(domains) for k in ids}
    rep = []
    for d, (_, ids) in enumerate(domains):
        votes = Counter(owner[m] for m in (mate.get(k, k) for k in ids) if m in owner)
        rep.append(max(votes, key=lambda e: (votes[e], e == d, -e)) if votes else d)

    def root(d: int) -> int:
        seen = set()
        while rep[d] != d and d not in seen:
            seen.add(d)
            d = rep[d]
        return d

    tone: dict[int, int] = {}
    return [tone.setdefault(root(d), len(tone)) for d in range(len(domains))]


def _chain_domains(sses: list[SSE], contacts: dict[frozenset[str], int], together=()) -> list[list[str]]:
    """One chain's domains, in sequence order: communities of its element contact graph (plus chain connectivity),
    each with at least DOMAIN_MIN_SSE_RESIDUES residues of elements (smaller ones fold into their best-connected
    neighbour), when they split the chain clearly (modularity at least 0.3). A βα row of `together`
    (layout.template_rows, drawn as one rigid unit) is never split: the communities holding its strands merge and
    its helices follow its strands, so the row lies in one domain's panel (1LDM's Rossmann sheet: the contacts alone
    cut it in two and put its α2 in the catalytic domain); what that leaves too small folds again (1TIM: a
    10-residue domain of two helices). [] when the chain reads as a single domain."""
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(s.id for s in sses)
    for key, n in contacts.items():
        a, b = sorted(key)
        if a in graph and b in graph:
            graph.add_edge(a, b, weight=float(n))
    for a, b in zip(sses, sses[1:]):
        w = graph.get_edge_data(a.id, b.id, {}).get("weight", 0.0)
        graph.add_edge(a.id, b.id, weight=w + 1.0)
    if graph.number_of_edges() == 0:
        return []
    size = {s.id: len(s) for s in sses}

    def fold(groups: list[set[str]]) -> list[set[str]]:  # a domain needs enough structure: the smallest folds
        groups = list(groups)  # into its best-connected neighbour
        while len(groups) > 1:
            small = min(groups, key=lambda g: sum(size[k] for k in g))
            if sum(size[k] for k in small) >= DOMAIN_MIN_SSE_RESIDUES:
                break
            groups.remove(small)
            link = [sum(graph.get_edge_data(a, b, {}).get("weight", 0.0) for a in small for b in g) for g in groups]
            groups[int(np.argmax(link))] |= small
        return groups

    groups = fold(set(g) for g in nx.community.greedy_modularity_communities(graph, weight="weight"))
    if len(groups) < 2 or nx.community.modularity(graph, [sorted(g) for g in groups], weight="weight") < 0.3:
        return []
    kind = {s.id: s.kind for s in sses}
    for keep in together:  # a βα row stays whole: its strands' communities merge, its helices follow its strands
        strands = {k for k in keep if kind.get(k) == "E"}
        if not strands:
            continue
        mine = [g for g in groups if g & strands]
        rest = [g - set(keep) for g in groups if not g & strands]
        groups = [g for g in [*rest, set().union(*mine, set(keep) & set(kind))] if g]
    groups = fold(groups)
    if len(groups) < 2:
        return []
    start = {s.id: s.start for s in sses}
    return sorted((sorted(g) for g in groups), key=lambda g: min(start[k] for k in g))
