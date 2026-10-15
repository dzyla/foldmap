# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Group strands into sheets from bridge pairing and order the strands side by side."""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations

import networkx as nx
import numpy as np

from .model import SSE, Bridge, Sheet

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


def _order_tree(tree: nx.Graph, first: int) -> tuple[list[int], bool]:
    """Order tree nodes along its longest path; branches are inserted beside their parent."""
    degrees = dict(tree.degree)
    if tree.number_of_nodes() == 1:
        return [first], False
    if max(degrees.values()) <= 2:
        end = min((n for n, d in degrees.items() if d == 1))
        return list(nx.dfs_preorder_nodes(tree, end)), False
    far = lambda src: max(nx.single_source_shortest_path_length(tree, src).items(), key=lambda kv: (kv[1], -kv[0]))[0]  # noqa: E731
    a = far(min(tree.nodes))
    b = far(a)
    path = nx.shortest_path(tree, a, b)
    order = list(path)
    placed = set(order)
    for node in nx.bfs_tree(tree, a):
        if node in placed:
            continue
        parent = next(p for p in tree.neighbors(node) if p in placed)
        order.insert(order.index(parent) + 1, node)
        placed.add(node)
    return order, True


_STRAND_RISE = 3.3  # Å per residue along a strand
_STACK_GAP = 5.0  # Å: two strands on one side of a partner whose spans along it are at most this far apart, and
_STACK_OVERLAP = 1.0  # ... overlap by no more than this, continue one another along it (a beta bulge, end to end)


def _span(s: SSE, at: np.ndarray, axis: np.ndarray) -> tuple[float, float]:
    """Where strand s lies along a partner's axis (Å from the partner's centroid `at`)."""
    c = float((s.centroid - at) @ axis)
    half = (len(s) - 1) * _STRAND_RISE / 2 * abs(float(s.axis @ axis))
    return c - half, c + half


def _columns(tree: nx.Graph, order: list[int], strands: list[SSE]) -> list[int]:
    """Column of each strand of `order`: a branch strand b set beside its partner p (_order_tree) that lies on the
    same side of p as the strand q next to it in the order, which also pairs with p, end to end with q along p (a
    bulged continuation: 5KQJ's F after B along E), shares q's column; empty when none does. Directions in 3D from
    the strands' centroids and axes, snapped."""
    from .layout import _snap, _unit

    share: dict[int, int] = {}
    for i, b in enumerate(order):
        for j in (i - 1, i + 1):
            if not 0 <= j < len(order) or b in share or order[j] in share:
                continue
            q = order[j]
            for p in sorted(set(tree.neighbors(b)) & set(tree.neighbors(q))):
                at, axis = _snap(strands[p].centroid), _snap(_unit(strands[p].axis))
                sb, sq = _snap(strands[b].centroid) - at, _snap(strands[q].centroid) - at
                if float((sb - (sb @ axis) * axis) @ (sq - (sq @ axis) * axis)) <= 0:
                    continue  # on opposite sides of p
                (b0, b1), (q0, q1) = _span(strands[b], at, axis), _span(strands[q], at, axis)
                overlap = min(b1, q1) - max(b0, q0)
                if -_STACK_GAP <= overlap <= _STACK_OVERLAP:
                    share[b] = q
                    break
    if not share:
        return []
    cols, k = [], -1
    for i, x in enumerate(order):
        if not (x in share and i > 0 and share[x] == order[i - 1]) and not (i > 0 and share.get(order[i - 1]) == x):
            k += 1
        cols.append(k)
    return cols


def build_sheets(sses: list[SSE], bridges: list[Bridge]) -> list[Sheet]:
    strands = [s for s in sses if s.kind == "E"]
    table = _edge_table(strands, bridges)
    weak = _edge_table(strands, bridges, minimum=1)  # single bridges, used only to close barrels
    graph = nx.Graph()
    graph.add_nodes_from(range(len(strands)))
    for (a, b), c in table.items():
        graph.add_edge(a, b, weight=sum(c.values()), kind=_majority(c))

    sheets: list[Sheet] = []
    for comp in sorted(nx.connected_components(graph), key=min):
        sub = graph.subgraph(comp)
        tree = nx.maximum_spanning_tree(sub, weight="weight")
        closed = tree.number_of_edges() < sub.number_of_edges()
        order, ambiguous = _order_tree(tree, min(comp))
        ambiguous = ambiguous or _flank_one_partner(tree, strands)
        direction = {order[0]: 1}
        for parent, child in nx.bfs_edges(tree, order[0]):
            flip = -1 if tree.edges[parent, child]["kind"] == "A" else 1
            direction[child] = direction[parent] * flip
        kinds = [tree.edges[a, b]["kind"] if tree.has_edge(a, b) else "?" for a, b in zip(order, order[1:])]
        ends = (min(order[0], order[-1]), max(order[0], order[-1]))
        if not closed and not ambiguous and len(order) >= _BARREL_MIN and ends in weak:
            closed = True  # a barrel whose last pairing is a lone bridge (DSSP counts it as a ladder too)
        columns = _columns(tree, order, strands) if ambiguous and not closed else []
        sheets.append(
            Sheet([strands[k] for k in order], [direction[k] for k in order], kinds, closed, ambiguous, columns)
        )
    return sheets
