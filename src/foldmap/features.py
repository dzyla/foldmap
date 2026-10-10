# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Facts derived from the structure that the layout and annotation layers need."""

from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .model import SSE, Backbone

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


DOMAIN_MIN_SSE_RESIDUES = 30  # residues in helices/strands an automatic domain needs


def auto_domains(sses: list[SSE], contacts: dict[frozenset[str], int]) -> list[tuple[str, list[str]]]:
    """Domains found as communities of the element contact graph (plus chain connectivity), named D1, D2...
    in sequence order; [] when the protein reads as a single domain."""
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(s.id for s in sses)
    for key, n in contacts.items():
        a, b = sorted(key)
        if a in graph and b in graph:
            graph.add_edge(a, b, weight=float(n))
    for a, b in zip(sses, sses[1:]):
        if a.chain == b.chain:
            w = graph.get_edge_data(a.id, b.id, {}).get("weight", 0.0)
            graph.add_edge(a.id, b.id, weight=w + 1.0)
    if graph.number_of_edges() == 0:
        return []
    groups = [set(g) for g in nx.community.greedy_modularity_communities(graph, weight="weight")]
    size = {s.id: len(s) for s in sses}
    while len(groups) > 1:  # a domain needs enough structure: fold the smallest into its best-connected neighbour
        small = min(groups, key=lambda g: sum(size[k] for k in g))
        if sum(size[k] for k in small) >= DOMAIN_MIN_SSE_RESIDUES:
            break
        groups.remove(small)
        link = [sum(graph.get_edge_data(a, b, {}).get("weight", 0.0) for a in small for b in g) for g in groups]
        groups[int(np.argmax(link))] |= small
    groups = [sorted(g) for g in groups]
    if len(groups) < 2:
        return []
    if nx.community.modularity(graph, groups, weight="weight") < 0.3:
        return []
    start = {s.id: (s.chain, s.start) for s in sses}
    groups.sort(key=lambda g: min(start[k] for k in g))
    return [(f"D{k}", g) for k, g in enumerate(groups, 1)]
