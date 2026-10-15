# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Group strands into sheets from bridge pairing and order the strands side by side."""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations

import networkx as nx

from .model import SSE, Accessory, Bridge, Sheet

_SPLIT_GAP = 4  # strands of one chain this close end to end are pieces of one long strand
_BARREL_MIN = 5  # strands a sheet needs before a lone end-to-end bridge may close it into a barrel
_MIN_PAIRS = 2  # bridges needed between two strands to count them as neighbours


def _edge_table(strands: list[SSE], bridges: list[Bridge], minimum: int = _MIN_PAIRS):
    owner: dict[int, int] = {}
    for k, s in enumerate(strands):
        for r in range(s.start, s.end + 1):
            owner[r] = k
    table: dict[tuple[int, int], Counter] = defaultdict(Counter)
    for i, j, kind in bridges:
        a, b = owner.get(i), owner.get(j)
        if a is not None and b is not None and a != b:
            table[(min(a, b), max(a, b))][kind] += 1
    return {e: c for e, c in table.items() if sum(c.values()) >= minimum}


def _majority(c: Counter) -> str:
    return "P" if c["P"] > c["A"] else "A"


def _flank_one_partner(tree: nx.Graph, strands: list[SSE]) -> bool:
    """True when two end-to-end strand pieces pair with the same strand: they sit on the same side
    of it in 3D, but a side-by-side order would put them on opposite sides."""
    for k in tree:
        for a, b in combinations(tree.neighbors(k), 2):
            x, y = sorted((strands[a], strands[b]), key=lambda s: s.start)
            if x.chain == y.chain and y.start - x.end - 1 <= _SPLIT_GAP:
                return True
    return False


def _hanging(tree: nx.Graph, track: list[int]) -> list[tuple[int, int]]:
    """Every tree node off the track, with the neighbour it hangs from (the one nearer the track), parents before
    children: what pruning the tree's off-track leaves to a fixed point removes, a branch of several strands whole."""
    seen, out, frontier = set(track), [], list(track)
    while frontier:
        nxt = []
        for p in frontier:
            for q in sorted(tree.neighbors(p)):
                if q not in seen:
                    seen.add(q)
                    out.append((q, p))
                    nxt.append(q)
        frontier = nxt
    return out


def _order_tree(tree: nx.Graph, first: int) -> tuple[list[int], list[tuple[int, int]]]:
    """The tree's longest path, in order (the sheet's track), and every node off it with the node it hangs from
    (_hanging). Plan 12 (the review's R1): an off-path strand was inserted right after its parent, between the parent
    and the parent's real neighbour (9FWZ's D:662-666 between D:645-659 and D:139-153, which pair); it is now kept
    off the track and drawn beside its parent (Sheet.accessory)."""
    degrees = dict(tree.degree)
    if tree.number_of_nodes() == 1:
        return [first], []
    if max(degrees.values()) <= 2:
        end = min((n for n, d in degrees.items() if d == 1))
        return list(nx.dfs_preorder_nodes(tree, end)), []
    far = lambda src: max(nx.single_source_shortest_path_length(tree, src).items(), key=lambda kv: (kv[1], -kv[0]))[0]  # noqa: E731
    a = far(min(tree.nodes))
    path = nx.shortest_path(tree, a, far(a))
    return path, _hanging(tree, path)


def _ring(tree: nx.Graph, graph: nx.Graph) -> list[int] | None:
    """A branched closed sheet's barrel: the longest cycle its pairing closes over the tree (for each pairing left
    out of the maximum spanning tree, the tree path between its two strands; the first of the longest, by strand
    index), as that path: its ends pair (the closure). None when that cycle has fewer than _BARREL_MIN strands (a
    bifurcation, not a barrel)."""
    best: list[int] = []
    for u, v in sorted(tuple(sorted(e)) for e in graph.edges if not tree.has_edge(*e)):
        path = nx.shortest_path(tree, u, v)
        if len(path) > len(best):
            best = path
    return best if len(best) >= _BARREL_MIN else None


def _branch_columns(track: list[int], off: list[tuple[int, int]], strands: list[SSE]) -> dict[int, int]:
    """The column each off-track strand is drawn in (Accessory.column): beside the strand it hangs from, on the side
    it lies on in 3D, sharing the column of the track strand there (stacked end to end with it, layout._column_heights;
    as C2a stacked two strands along one edge of a partner) or past the track's end. The side: the
    strand's offset from its parent against the track's local across direction (from the track strand before the
    parent's column to the one after it), both seen across the parent's axis; snapped."""
    from .layout import _snap, _unit

    n = len(track)
    column = {node: k for k, node in enumerate(track)}
    for node, parent in off:
        c = column[parent]
        at = min(max(c, 0), n - 1)
        across = _snap(strands[track[min(at + 1, n - 1)]].centroid) - _snap(strands[track[max(at - 1, 0)]].centroid)
        axis = _unit(strands[parent].axis)
        v = _snap(strands[node].centroid) - _snap(strands[parent].centroid)
        side = float((v - (v @ axis) * axis) @ (across - (across @ axis) * axis))
        column[node] = c + (1 if side >= 0 else -1)
    return {node: column[node] for node, _ in off}


def _other_rings(graph: nx.Graph, nodes: set[int], ring: list[int]) -> list[list[int]]:
    """The parts of a component, its ring taken out, that close a barrel of their own (a cycle of at least
    _BARREL_MIN strands, _ring): each a sheet of its own (1EYG: subunit B's OB-fold barrel, joined to A's by one
    pairing, A:1005-1021 with B:2005-2021). The rest stays with the ring, its strands accessory."""
    rest = graph.subgraph(nodes - set(ring))
    out = []
    for part in sorted(nx.connected_components(rest), key=min):
        sub = rest.subgraph(part)
        tree = nx.maximum_spanning_tree(sub, weight="weight")
        if tree.number_of_edges() < sub.number_of_edges() and _ring(tree, sub) is not None:
            out.append(sorted(part))
    return out


def build_sheets(sses: list[SSE], bridges: list[Bridge]) -> list[Sheet]:
    """The sheets of these strands: strands joined by at least _MIN_PAIRS bridges, per connected component its
    maximum spanning tree. The track (Sheet.strands): the tree's longest path; a closed sheet's (its pairing has a
    cycle): a branched tree's longest cycle (_ring), a path tree's path. A part off that ring closing a barrel of its
    own is a sheet of its own (_other_rings). Every other strand is an accessory strand (Sheet.accessory), drawn
    beside the strand it hangs from (_branch_columns), never between two strands that pair. A sheet whose pairing
    has a cycle but whose track's ends do not pair (or whose longest cycle is under _BARREL_MIN strands) is drawn
    open, flagged `unclosed`."""
    strands = [s for s in sses if s.kind == "E"]
    table = _edge_table(strands, bridges)
    weak = _edge_table(strands, bridges, minimum=1)  # single bridges, used only to close barrels
    graph = nx.Graph()
    graph.add_nodes_from(range(len(strands)))
    for (a, b), c in table.items():
        graph.add_edge(a, b, weight=sum(c.values()), kind=_majority(c))

    sheets: list[Sheet] = []
    pending = [sorted(comp) for comp in sorted(nx.connected_components(graph), key=min)]
    comp_of = {k: comp for comp in pending for k in comp}
    while pending:
        comp = pending.pop(0)
        sub = graph.subgraph(comp)
        tree = nx.maximum_spanning_tree(sub, weight="weight")
        cyclic = tree.number_of_edges() < sub.number_of_edges()  # a cycle was cut at its weakest pairing
        order, off = _order_tree(tree, min(comp))
        ring = _ring(tree, sub) if cyclic and off else None  # a branched barrel: its ring, every other strand off it
        if ring is not None:
            order, off = ring, _hanging(tree, ring)
            split = _other_rings(graph, set(comp), ring)
            gone = {k for part in split for k in part}
            off = [(k, p) for k, p in off if k not in gone]
            pending[:0] = split  # next, in order
        ends = (min(order[0], order[-1]), max(order[0], order[-1]))
        closed = cyclic and (ring is not None or not off) and ends in weak
        unclosed = cyclic and not closed  # drawn open: its closure is not supported by pairing
        ambiguous = bool(off) or _flank_one_partner(tree, strands)
        direction = {order[0]: 1}
        for parent, child in nx.bfs_edges(tree, order[0]):
            flip = -1 if tree.edges[parent, child]["kind"] == "A" else 1
            direction[child] = direction[parent] * flip
        kinds = [tree.edges[a, b]["kind"] for a, b in zip(order, order[1:])]
        if not cyclic and not ambiguous and len(order) >= _BARREL_MIN and ends in weak:
            closed = True  # a barrel whose last pairing is a lone bridge (DSSP counts it as a ladder too)
        columns = _branch_columns(order, off, strands)
        accessory = [Accessory(strands[k], strands[p], columns[k], direction[k]) for k, p in off]
        mine = {*order, *(k for k, _ in off)}
        apart = set(comp_of[comp[0]]) - mine  # strands of its component drawn in another sheet
        joined = [(strands[a], strands[b]) for a in sorted(mine) for b in sorted(graph.neighbors(a)) if b in apart]
        sheets.append(
            Sheet(
                [strands[k] for k in order],
                [direction[k] for k in order],
                kinds,
                closed,
                ambiguous,
                accessory=accessory,
                unclosed=unclosed,
                joined=joined,
            )
        )
    return sheets
