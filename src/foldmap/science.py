# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Science checks (S-checks): is what a figure says about the structure true of its coordinates? They are judged
apart from the page checks (figcheck A1-A20: overlaps, text, order as drawn), since a clean page can still draw a
wrong sheet (plan 12, the external review's FimD 9FWZ: a strand wedged between two strands that pair).

- S1: every two strands drawn as plain neighbours in a sheet (adjacent columns of one strand each) are joined by at
  least one DSSP bridge; every two adjacent columns are; and a strand drawn in a shared column (an accessory strand
  end to end with its partner) pairs with a strand of its column or of a column next to it.
- S2: a closed barrel's closure (its ghost, repeating the first strand beside the last) is joined by a bridge.
- S3: a closed barrel's strand count (its accessory strands left out) equals the longest cycle of the pairing graph
  (strands joined by at least sheets._MIN_PAIRS bridges, and the drawn closure by one) among its strands.
- S4, S5: residue marks at their residue, automatic domains against a reference (plan 12 Tasks 4 and 3): hooks, not
  yet measured (None).
- S6: an accessory strand (Sheet.accessory, Layout.accessory: off its sheet's run of paired strands) is drawn at its
  partner: in its partner's own column, next to it (end to end), or in a column beside the partner's that holds no
  other strand, their spans overlapping; never in another strand's column, which reads as pairing with that strand
  (fix round 1: 9FWZ's D:662-666 stood in D:634-641's column, beside its partner D:645-659).

Two readings, one rule: sheet_checks reads the sheets as built (sheets.build_sheets, for `foldmap summary`), and
figcheck.science_checks reads the drawn figure (its sheet blocks, columns by where the strands stand)."""

from __future__ import annotations

from collections.abc import Callable, Iterable

import networkx as nx

from .model import SSE, Bridge, Sheet

SCIENCE = {
    "S1": "sheet neighbours pair",
    "S2": "barrel closure pairs",
    "S3": "barrel strand count is its pairing cycle",
    "S4": "residue marks at their residue",
    "S5": "automatic domains agree with a reference",
    "S6": "accessory drawn at its partner",
}
NOT_YET = ("S4", "S5")  # hooks: measured by later tasks; None until then
_MAX_CYCLES = 20_000  # simple cycles looked at before S3 gives up on a sheet (sheet graphs are sparse)

# A drawn sheet: its columns left to right, each the (strand id, is a ghost) standing in it.
Columns = list[list[tuple[str, bool]]]
# Where a drawn sheet's strands stand: (strand id, column, span (bottom, top) or None when not known, is a ghost).
Spot = tuple[str, int, tuple[float, float] | None, bool]


def partner_findings(spots: list[Spot], partner: dict[str, str], noted=frozenset()) -> list[str]:
    """S6 of one sheet: each accessory strand (`partner`: accessory id -> its partner's id) whose partner stands in
    the sheet too must stand in the partner's column, next to it when the spans are known (no strand between them),
    or in the column beside the partner's, alone there, their spans overlapping when known; or, when its partner
    had no end free for it (fix round 2 ruling: 7ARD), its pairing must be named in the notes (`noted`: (accessory,
    partner) pairs a note says "also pairs with")."""
    out = []
    real = {i: (c, span) for i, c, span, ghost in spots if not ghost}
    for a, p in sorted(partner.items()):
        if a not in real or p not in real or (a, p) in noted:
            continue
        (ca, sa), (cp, sp) = real[a], real[p]
        others = [(i, span) for i, c, span, _ in spots if c == ca and i != a]
        if ca == cp:
            if sa is None or sp is None:
                continue
            lo, hi = sorted((sa, sp))
            between = [i for i, span in others if i != p and span and lo[1] <= (span[0] + span[1]) / 2 <= hi[0]]
            if between:
                out.append(f"{a} drawn in its partner {p}'s column but not next to it ({', '.join(between)} between)")
            continue
        if abs(ca - cp) == 1 and not others and (sa is None or sp is None or min(sa[1], sp[1]) > max(sa[0], sp[0])):
            continue
        held = [i for i, _ in others]
        where = f"in {', '.join(held)}'s column" if held else "away from it"
        out.append(f"{a} drawn {where}, not at its partner {p}")
    return out


def bridge_counts(sses: Iterable[SSE], bridges: list[Bridge]) -> Callable[[str, str], int]:
    """The number of DSSP bridges between two strands, by id (0 for any pair not bridged)."""
    from .sheets import _edge_table

    strands = [s for s in sses if s.kind == "E"]
    table = _edge_table(strands, bridges, minimum=1)
    by_pair = {frozenset((strands[a].id, strands[b].id)): sum(c.values()) for (a, b), c in table.items()}
    return lambda a, b: by_pair.get(frozenset((a, b)), 0)


def neighbour_findings(columns: Columns, count: Callable[[str, str], int]) -> tuple[list[str], list[str]]:
    """S1 and S2 findings of one drawn sheet. Adjacent columns must share a bridge (two plain neighbours: the two
    strands); a strand in a shared column must pair with a strand of its column (an accessory end to end with its
    partner, fix round 1) or of a column beside it. A finding naming a ghost is the closure's (S2)."""
    s1: list[str] = []
    s2: list[str] = []

    def put(text: str, ghost: bool) -> None:
        (s2 if ghost else s1).append(text)

    def name(entry: tuple[str, bool]) -> str:
        return f"{entry[0]} (ghost)" if entry[1] else entry[0]

    for a, b in zip(columns, columns[1:]):
        if not any(count(x, y) for x, _ in a for y, _ in b):
            ghost = any(g for _, g in a + b)
            if len(a) == 1 and len(b) == 1:
                put(f"{name(a[0])} | {name(b[0])} drawn side by side with no bridge", ghost)
            else:
                put(f"columns {'/'.join(map(name, a))} | {'/'.join(map(name, b))} share no bridge", ghost)
    for k, col in enumerate(columns):
        if len(col) < 2:
            continue
        beside = [y for j in (k - 1, k, k + 1) if 0 <= j < len(columns) for y, _ in columns[j]]
        for entry in col:
            if not any(count(entry[0], y) for y in beside if y != entry[0]):
                put(f"{name(entry)} drawn in a shared column beside none of its partners", entry[1])
    return s1, s2


def adjacent_pairs(columns: Columns) -> list[tuple[str, str]]:
    """The strand pairs a drawn sheet stands side by side: each strand (or ghost) with each in the next column."""
    return [(x, y) for a, b in zip(columns, columns[1:]) for x, _ in a for y, _ in b]


def pairing_cycle(ids: Iterable[str], count: Callable[[str, str], int], closure=()) -> int | None:
    """The longest simple cycle among these strands in the pairing graph: strands joined by at least
    sheets._MIN_PAIRS bridges, and each pair in `closure` (drawn side by side: adjacent_pairs) by one, as a barrel's
    closing ladder may be (1TIM's beta7-beta8). 0 when there is no cycle; None when there are too many cycles to
    look at."""
    from .sheets import _MIN_PAIRS

    ids = sorted(set(ids))
    graph = nx.Graph()
    graph.add_nodes_from(ids)
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            n = count(a, b)
            if n >= _MIN_PAIRS or (n and frozenset((a, b)) in {frozenset(p) for p in closure}):
                graph.add_edge(a, b)
    core = nx.k_core(graph, 2)  # strands on no cycle cannot be in a barrel's ring
    best = 0
    for seen, cycle in enumerate(nx.simple_cycles(core)):
        if seen >= _MAX_CYCLES:
            return None
        best = max(best, len(cycle))
    return best


def count_finding(label: str, ring: list[str], others: Iterable[str], count, closure) -> str | None:
    """S3 of one closed barrel: `ring` its drawn barrel strands, `others` its accessory strands, `closure` the pairs
    it draws side by side (adjacent_pairs)."""
    cycle = pairing_cycle([*ring, *others], count, closure)
    if cycle is None:
        return f"{label}: too many pairing cycles to count"
    if cycle != len(ring):
        return f"{label}: {len(ring)} barrel strands drawn, the pairing cycle has {cycle}"
    return None


def sheet_columns(sheet: Sheet) -> Columns:
    """A built sheet's drawn columns (layout._sheet_item draws it so): its strands in order, a closed sheet's ghost of
    the first after the last, and each accessory strand in its column (Sheet.accessory)."""
    n = len(sheet.strands)
    cols: dict[int, list[tuple[str, bool]]] = {}
    for k, s in enumerate(sheet.strands):
        cols.setdefault(k, []).append((s.id, False))
    if sheet.closed:
        cols.setdefault(n, []).append((sheet.strands[0].id, True))
    for acc in sheet.accessory:
        cols.setdefault(acc.column, []).append((acc.strand.id, False))
    return [cols[k] for k in sorted(cols)]


def sheet_checks(sheets: list[Sheet], sses: Iterable[SSE], bridges: list[Bridge]) -> dict[str, list[str] | None]:
    """S1-S3 and S6 of the sheets as built (sheets.build_sheets), sheet k named "sheet k" from 1 as `foldmap summary`
    numbers them (S6 by column only: where in a column layout stacks them is the figure's); S4 and S5 None (not yet
    measured)."""
    count = bridge_counts(sses, bridges)
    out: dict[str, list[str] | None] = {"S1": [], "S2": [], "S3": [], **dict.fromkeys(NOT_YET), "S6": []}
    for n, sheet in enumerate(sheets, 1):
        columns = sheet_columns(sheet)
        spots = [(i, c, None, ghost) for c, col in enumerate(columns) for i, ghost in col]
        partner = {a.strand.id: a.parent.id for a in sheet.accessory}
        out["S6"] += [f"sheet {n}: {t}" for t in partner_findings(spots, partner)]
        s1, s2 = neighbour_findings(columns, count)
        out["S1"] += [f"sheet {n}: {t}" for t in s1]
        out["S2"] += [f"sheet {n}: {t}" for t in s2]
        if sheet.closed and len(sheet.strands) > 2:
            ring = [s.id for s in sheet.strands]
            found = count_finding(
                f"sheet {n}", ring, [a.strand.id for a in sheet.accessory], count, adjacent_pairs(columns)
            )
            if found:
                out["S3"].append(found)
    return out


def science_lines(found: dict[str, list[str] | None]) -> list[str]:
    """One line per check (pass, its count of findings, or not yet measured), then each finding indented."""
    lines = []
    for key, what in SCIENCE.items():
        got = found.get(key)
        if got is None:
            lines.append(f"{key} {what}: not yet measured")
            continue
        lines.append(f"{key} {what}: {'pass' if not got else f'{len(got)} problem(s)'}")
        lines += [f"  {t}" for t in got]
    return lines
