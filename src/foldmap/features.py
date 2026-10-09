"""Facts derived from the structure that the layout and annotation layers need."""

from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .model import Backbone, SSE

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


def _segment(s: SSE, n: int = 20) -> np.ndarray:
    half = (len(s) - 1) * 0.75  # 1.5 Å rise per residue
    return np.linspace(s.centroid - s.axis * half, s.centroid + s.axis * half, n)


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

    pieces = {}
    for a, b in zip(helices, helices[1:]):  # kinks
        if a.chain == b.chain and b.start - a.end - 1 <= KINK_GAP and float(a.axis @ b.axis) > KINK_COS:
            pieces.setdefault(a.id, a)
            pieces.setdefault(b.id, b)
            join(a.id, b.id)
    linked = set()
    long_ = [s for s in helices if len(s) >= BUNDLE_MIN]
    for i, a in enumerate(long_):
        pa = _segment(a)
        for b in long_[i + 1 :]:
            if a.chain == b.chain or abs(float(a.axis @ b.axis)) < BUNDLE_COS:
                continue
            d = np.linalg.norm(pa[:, None] - _segment(b)[None], axis=2)
            near = d.min(1) if len(a) <= len(b) else d.min(0)
            if (near < BUNDLE_DIST).mean() >= BUNDLE_COVER:
                join(a.id, b.id)
                linked |= {a.id, b.id}
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
