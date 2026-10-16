# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Science checks (S-checks): is what a figure says about the structure true of its coordinates? They are judged
apart from the page checks (figcheck A1-A20: overlaps, text, order as drawn), since a clean page can still draw a
wrong sheet (plan 12, the external review's FimD 9FWZ: a strand wedged between two strands that pair).

- S1: any two strands drawn side by side in a sheet (adjacent columns, their drawn spans overlapping) are joined by
  at least one DSSP bridge; two strands in one column (collinear, end to end: no pairing) never overlap in height
  (fix round 3 ruling: 8UTF's a:266-273 stood beside the ghost of c:440-442, 0 bridges).
- S2: the same for a pair naming a ghost: a closed barrel's closure (its ghost, repeating the first strand beside the
  last) is joined by a bridge.
- S3: a closed barrel's strand count (its accessory strands left out) equals the longest cycle of the pairing graph
  (strands joined by at least sheets._MIN_PAIRS bridges, and the drawn closure by one) among its strands; and a sheet
  drawn open holds no such cycle of _BARREL_MIN (5) strands or more.
- S4, S5: residue marks at their residue, automatic domains against a reference (plan 12 Tasks 4 and 3): hooks, not
  yet measured (None).
- S6: an accessory strand (Sheet.accessory, Layout.accessory: off its sheet's run of paired strands) is drawn at its
  partner: in a column next to the partner's, their drawn spans overlapping; or its pairing is named in the notes
  (Layout.unshown: no height beside the partner was free). Fix round 3 ruling: a shared column means collinear, so an
  accessory stands in the column on its 3D side of the partner, end to end with the strand there (9FWZ's D:662-666
  beside D:645-659, on D:634-641's side).

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


def _meet(a: tuple[float, float] | None, b: tuple[float, float] | None) -> bool:
    """Do two drawn spans (bottom, top) overlap? Unknown spans (None) count as overlapping."""
    return a is None or b is None or min(a[1], b[1]) - max(a[0], b[0]) > _TOUCH


_TOUCH = 1e-6  # page units two spans must share to overlap


def partner_findings(spots: list[Spot], partner: dict[str, str], unshown=frozenset()) -> list[str]:
    """S6 of one sheet: each accessory strand (`partner`: accessory id -> its partner's id) whose partner stands in
    the sheet too stands in a column next to the partner's, their spans overlapping (unknown spans: by column only);
    else its pairing is among `unshown` (Layout.unshown: (accessory, partner) pairs the notes name)."""
    out = []
    real = {i: (c, span) for i, c, span, ghost in spots if not ghost}
    for a, p in sorted(partner.items()):
        if a not in real or p not in real or (a, p) in unshown:
            continue
        (ca, sa), (cp, sp) = real[a], real[p]
        if abs(ca - cp) == 1 and _meet(sa, sp):
            continue
        where = (
            "in its partner's column"
            if ca == cp
            else ("beside it, not level with it" if abs(ca - cp) == 1 else "away from it")
        )
        out.append(f"{a} drawn {where}, not beside its partner {p}")
    return out


def bridge_counts(sses: Iterable[SSE], bridges: list[Bridge]) -> Callable[[str, str], int]:
    """The number of DSSP bridges between two strands, by id (0 for any pair not bridged)."""
    from .sheets import _edge_table

    strands = [s for s in sses if s.kind == "E"]
    table = _edge_table(strands, bridges, minimum=1)
    by_pair = {frozenset((strands[a].id, strands[b].id)): sum(c.values()) for (a, b), c in table.items()}
    return lambda a, b: by_pair.get(frozenset((a, b)), 0)


def neighbour_findings(
    spots: list[Spot], count: Callable[[str, str], int], closure: Iterable[str] | None = None
) -> tuple[list[str], list[str]]:
    """S1 and S2 findings of one drawn sheet (`spots`): two strands in adjacent columns whose spans overlap must
    share a bridge; two strands in one column must not overlap. A finding naming a ghost of the barrel's own closure
    (`closure`: those ghosts' ids; None: every ghost) is S2; a ghost lent by a symmetric copy is S1's."""
    closure = None if closure is None else set(closure)

    def closing(i: str, ghost: bool) -> bool:
        return ghost and (closure is None or i in closure)

    s1: list[str] = []
    s2: list[str] = []

    def name(i: str, ghost: bool) -> str:
        return f"{i} (ghost)" if ghost else i

    ordered = sorted(spots, key=lambda sp: (sp[1], sp[0], sp[3]))
    for k, (a, ca, sa, ga) in enumerate(ordered):
        for b, cb, sb, gb in ordered[k + 1 :]:
            if cb > ca + 1:
                break
            if sa is None or sb is None or not _meet(sa, sb) or a == b:
                continue  # an accessory's height is the figure's: the sheets as built do not judge it here
            if cb == ca:
                if sa is not None and sb is not None:
                    put = s2 if closing(a, ga) or closing(b, gb) else s1
                    put.append(f"{name(a, ga)} and {name(b, gb)} overlap in one column")
            elif not count(a, b):
                put = s2 if closing(a, ga) or closing(b, gb) else s1
                put.append(f"{name(a, ga)} | {name(b, gb)} drawn side by side with no bridge")
    return s1, s2


def adjacent_pairs(spots: list[Spot]) -> list[tuple[str, str]]:
    """The strand pairs a drawn sheet stands side by side: adjacent columns, spans overlapping."""
    return [(a, b) for a, ca, sa, _ in spots for b, cb, sb, _ in spots if cb == ca + 1 and _meet(sa, sb) and a != b]


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


def open_finding(label: str, ids: Iterable[str], count) -> str | None:
    """S3 of a sheet drawn open: a cycle of its pairing of _BARREL_MIN (5) strands or more is a barrel drawn open."""
    from .sheets import _BARREL_MIN

    cycle = pairing_cycle(ids, count)
    if cycle is not None and cycle >= _BARREL_MIN:
        return f"{label}: drawn open, its pairing closes a ring of {cycle} strands"
    return None


def count_finding(label: str, ring: list[str], others: Iterable[str], count, closure) -> str | None:
    """S3 of one closed barrel: `ring` its drawn barrel strands, `others` its accessory strands, `closure` the pairs
    it draws side by side (adjacent_pairs)."""
    cycle = pairing_cycle([*ring, *others], count, closure)
    if cycle is None:
        return f"{label}: too many pairing cycles to count"
    if cycle != len(ring):
        return f"{label}: {len(ring)} barrel strands drawn, the pairing cycle has {cycle}"
    return None


def sheet_spots(sheet: Sheet) -> list[Spot]:
    """A built sheet as layout._sheet_item draws it, by column: its strands in order (on the row: spans (-1, 1)), a
    closed sheet's ghost of the first after the last, and each accessory strand in its column (Sheet.accessory; its
    height is the figure's: None)."""
    n = len(sheet.strands)
    out: list[Spot] = [(s.id, k, (-1.0, 1.0), False) for k, s in enumerate(sheet.strands)]
    if sheet.closed:
        out.append((sheet.strands[0].id, n, (-1.0, 1.0), True))
    out += [(acc.strand.id, acc.column, None, False) for acc in sheet.accessory]
    return out


def sheet_checks(sheets: list[Sheet], sses: Iterable[SSE], bridges: list[Bridge]) -> dict[str, list[str] | None]:
    """S1-S3 and S6 of the sheets as built (sheets.build_sheets), sheet k named "sheet k" from 1 as `foldmap summary`
    numbers them (an accessory's height is the figure's: judged by column only); S4 and S5 None (not yet
    measured)."""
    count = bridge_counts(sses, bridges)
    out: dict[str, list[str] | None] = {"S1": [], "S2": [], "S3": [], **dict.fromkeys(NOT_YET), "S6": []}
    for n, sheet in enumerate(sheets, 1):
        spots = sheet_spots(sheet)
        partner = {a.strand.id: a.parent.id for a in sheet.accessory}
        out["S6"] += [f"sheet {n}: {t}" for t in partner_findings(spots, partner)]
        s1, s2 = neighbour_findings(spots, count)
        out["S1"] += [f"sheet {n}: {t}" for t in s1]
        out["S2"] += [f"sheet {n}: {t}" for t in s2]
        ids = [s.id for s in sheet.strands]
        if sheet.closed and len(sheet.strands) > 2:
            found = count_finding(f"sheet {n}", ids, partner, count, adjacent_pairs(spots))
        else:
            found = open_finding(f"sheet {n}", [*ids, *partner], count)
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
