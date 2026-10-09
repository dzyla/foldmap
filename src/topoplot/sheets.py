"""Group strands into sheets from bridge pairing and order the strands side by side."""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations

import networkx as nx

from .model import Bridge, SSE, Sheet

_SPLIT_GAP = 4  # strands of one chain this close end to end are pieces of one long strand
_MIN_PAIRS = 2  # bridges needed between two strands to count them as neighbours


def _edge_table(strands: list[SSE], bridges: list[Bridge]):
    owner: dict[int, int] = {}
    for k, s in enumerate(strands):
        for r in range(s.start, s.end + 1):
            owner[r] = k
    table: dict[tuple[int, int], Counter] = defaultdict(Counter)
    for i, j, kind in bridges:
        a, b = owner.get(i), owner.get(j)
        if a is not None and b is not None and a != b:
            table[(min(a, b), max(a, b))][kind] += 1
    return {e: c for e, c in table.items() if sum(c.values()) >= _MIN_PAIRS}


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


def build_sheets(sses: list[SSE], bridges: list[Bridge]) -> list[Sheet]:
    strands = [s for s in sses if s.kind == "E"]
    table = _edge_table(strands, bridges)
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
        kinds = [
            tree.edges[a, b]["kind"] if tree.has_edge(a, b) else "?" for a, b in zip(order, order[1:])
        ]
        sheets.append(
            Sheet([strands[k] for k in order], [direction[k] for k in order], kinds, closed, ambiguous)
        )
    return sheets
