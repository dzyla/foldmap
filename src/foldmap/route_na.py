# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Nucleic-acid strands on the page: route_nucleic draws every single-stranded stretch (runs and joints, na_runs;
free ends, traces and cropped gaps, na_ends) and places the text beside them (na_labels). The records and the
geometry they share are in na_paths; the names other modules use are re-exported here."""

from __future__ import annotations

import logging
import math

import numpy as np

from .layout import (
    _GLYPH_W,
    NT_ROW,
    NT_STEP,
    Layout,
    _run_band,
    glyph,
    label_boxes,
    lines_meeting,
    segments_of,
    stem_keep_box,
)
from .na_ends import SHORT_GAP, _face_label, _gap, _label_room, _tail, _trace_clear
from .na_labels import _boxes_of, _in_sight, _label_gaps, _place_letters, _ties, _unbunch, _widen
from .na_paths import (
    NAPath,
    NATie,
    _face,
    _near,
    _Strands,
    _xy,
    among,
    crosses_body,
    graze,
    in_body,
    midway,
    outside_own,
    path_name,
)
from .na_runs import _strand_between
from .route import _element_obstacles, _Grid, _Router

__all__ = [
    "SHORT_GAP",
    "NAPath",
    "NATie",
    "among",
    "crosses_body",
    "graze",
    "in_body",
    "midway",
    "outside_own",
    "path_name",
    "route_nucleic",
]

log = logging.getLogger(__name__)


def _use(grid: _Grid, pts, beside: np.ndarray | None = None) -> None:
    """Mark the cells under a drawn strand as used, so later routes avoid crossing it (at the reuse cost); and in
    `beside` the cells beside them, which a run with nucleotides on it pays to pass through as well (there its
    nucleotides would touch the strand's, and crossing the strand costs more than crossing a cell)."""
    under = {c for a, b in zip(pts, pts[1:]) for c in grid.cells_along(a, b)}
    for c in sorted(under):
        grid.used[c] += 1
    if beside is not None:
        near = {(x + dx, y + dy) for x, y in under for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))} - under
        for c in sorted(near):
            if 0 <= c[0] < grid.nx and 0 <= c[1] < grid.ny:
                beside[c] += 1


_END_TURNS = (0.0, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2, 3 * np.pi / 4, -3 * np.pi / 4, np.pi)
_END_FINE = tuple(  # every 5° between those, the least turn first, when none of them fits
    x for k in range(1, 36) if k % 9 for x in (k * np.pi / 36, -k * np.pi / 36)
)
_END_REACH = (1.0, 1.4, 1.8, 2.4, 3.0)  # times the label's distance from its end: the same first, then further out
# (2.4 and 3, only where the line to its end crosses no element but the one the end lies under: na_labels._in_sight)


def _own_ends(router, layout: Layout, fixed: dict, at: dict, ends: list, gapped: set, lines=(), loops=()) -> list:
    """The strands' end labels once every free end is drawn (fix round 1, review Q6: 9CJI's 3' end, drawn after the
    5' label was placed, started beside that label, so "5′" read as naming it): a label of a free end that lies
    nearer another drawn nucleotide (a free end's, a hairpin loop's, a stem's rail) than its own end nucleotide moves
    round that nucleotide, at the same distance, to the first turn (_END_TURNS, then the turns between them, _END_FINE)
    where its own end is the nearest and it has room (_label_room), else at the next distance out (_END_REACH; plan 12:
    30QM's 3' end, laid out anew once its sheets' branch strands moved, between nucleotides 103, 104 and 119, has room
    only about 100° round); else it stays. Run once the free ends and hairpin loops are drawn (`fixed`), and again
    with the runs (30QM's 3' label lay 0.61 from a run's nucleotide, 1.1 from its own end), when a label another
    strand's drawn line (`lines`: (strand, points) of the runs) runs through, as figcheck's A6 reads a line by a text,
    moves the same way (random-100 C2b-6: 7Y7C's 3' end of X:37-49 under run A:652-653), as does one a protein loop
    (`loops`, polylines) runs through, the label of a strand ending in a stem too (random-100 held-out 7NZ2's 5' end
    of K2:4-80 on A4:1270-1279>A4:1281-1292, where _face_label found no spot clear and took the first). A
    trace's end labels (strands drawn at their projected place, laid over the protein that holds them) move too, and
    also when the label stands on an element or stem (held-out 2FD8: the 5' end of B:501-503 on a strand). A spot
    further out than _END_REACH[2] must be in sight of its end (na_labels._in_sight): 2FD8's 5' end, under strand E,
    has no such spot with room, and its label stays where it was."""
    pts = [
        (key, np.asarray(xy, float)) for p in fixed.values() for key, xy in zip(((p.strand, k) for k in p.nts), p.nt_xy)
    ]
    pts += [((key[0], key[1]), np.asarray(d.point(key), float)) for key, d in sorted(at.items())]
    tails = {(p.strand, "5′" if p.nts[0] == 0 else "3′"): p for p in fixed.values() if p.kind == "tail"}
    traces = [p for p in fixed.values() if p.kind == "trace" and p.nts]
    tails.update({(p.strand, e): p for p in traces for e in ("5′", "3′")})

    def nearest(xy):
        return min(pts, key=lambda kp: (float(np.hypot(*(kp[1] - xy))), kp[0]))[0] if pts else None

    gh = glyph(layout)[1] / 2 + 0.05  # an end label's box, as routes keep off it
    band = _run_band(layout)
    every = segments_of([q for _, q in lines] + list(loops))
    segs = {s: segments_of([q for t, q in lines if t != s] + list(loops)) for s in sorted({t for t, _ in lines})}
    bodies = [*layout.placed.values(), *layout.stems]
    rects = _boxes_of(bodies)

    def crossed(s: int, xy) -> bool:  # another strand's line, or a protein loop, through the label at xy
        other = segs.get(s, every)
        box = (xy[0] - gh, xy[1] - gh, xy[0] + gh, xy[1] + gh)
        return bool(len(other)) and bool(lines_meeting(box, other, band).any())

    def on_body(s: int, xy) -> bool:  # a trace's label on an element or stem (_label_room's test)
        return s in layout.na_traces and any(_near(xy, e, gh) for e in bodies)

    out = []
    for s, e, xy in ends:
        path = tails.get((s, e))
        last = len(layout.nucleic.strands[s]) - 1
        face = (s, 0 if e == "5′" else last)  # a strand ending in a stem: its label is beyond the rail's end
        if (s, e) in gapped or (path is None and face not in at) or (path is not None and not path.nt_xy):
            out.append((s, e, xy))
            continue
        if path is not None:
            own = (s, path.nts[0] if e == "5′" else path.nts[-1])
            end_xy = np.asarray(path.nt_xy[0] if e == "5′" else path.nt_xy[-1], float)
            stays = nearest(np.asarray(xy, float)) == own and not crossed(s, xy) and not on_body(s, xy)
        else:  # placed by _face_label, clear of the protein loops' cells where it found room: moved off a line
            own, end_xy = face, np.asarray(at[face].point(face), float)
            stays = not crossed(s, xy)
        if stays:
            out.append((s, e, xy))
            continue
        v = np.subtract(xy, end_xy)
        r = float(np.hypot(*v))
        others = [q for t, f, q in ends if (t, f) != (s, e) and (t, f) not in gapped]
        moved = xy
        for reach, turn in ((k, t) for k in _END_REACH for t in (*_END_TURNS, *_END_FINE)):
            c, n = math.cos(turn), math.sin(turn)
            spot = end_xy + reach * np.array([v[0] * c - v[1] * n, v[0] * n + v[1] * c])
            if r > 0 and nearest(spot) == own and _label_room(router, layout, spot, others) and not crossed(s, spot):
                if reach > _END_REACH[2] and not _in_sight(spot, end_xy, bodies, rects):
                    continue  # further out: only in sight of its own end
                moved = (float(spot[0]), float(spot[1]))
                break
        out.append((s, e, moved))
    return out


def route_nucleic(layout: Layout, router: _Router | None = None, loops=()) -> list[NAPath]:
    """Every single-stranded stretch as drawn: runs routed around the elements like loops (as are the joints
    between pieces of a bent duplex), hairpin loops as arcs, free ends straight on, strands without pairs along
    their projection, cropped gaps dashed. Every drawing is checked by _clear; one that fails is a fallback.
    `loops`: the protein loops already drawn (polylines), whose crossings a run weighs against crossing strands.
    Also fills layout.na_ends (where each strand's 5' and 3' labels go) and layout.na_letters (where each
    single-stranded nucleotide's letter goes), and widens layout.bounds to hold it all."""
    na = layout.nucleic
    if router is None:
        keep_out = [b for _, b in label_boxes(layout, [], route=True)]
        router = _Router(_Grid(layout, keep_out), _element_obstacles(layout), keep_out)
    at = {key: d for d in layout.stems for key in d.axial}
    kind_of = {(r.strand, r.lo): r.kind for r in na.structure.runs}
    ends: list[tuple[int, str, tuple[float, float]]] = []
    exits: dict[tuple[int, int], tuple] = {}  # strand end -> (port, outward direction, owner) a gap leaves from

    gapped = {(g.before, "3′") for g in na.gaps} | {(g.after, "5′") for g in na.gaps}
    labels: list[tuple[float, float]] = []  # the end labels placed so far (a cropped end into a gap has none)

    def room(s: int, e: str):  # the labels an end's own label must keep clear of; None: it has no label
        return None if (s, e) in gapped else labels

    def end(s: int, e: str, xy) -> None:
        ends.append((s, e, xy))
        if (s, e) not in gapped:
            labels.append(xy)

    strands = _Strands(router.grid)  # the strands drawn so far: hard obstacles for a run's first routing pass
    crossed = _Strands(None)  # the protein loops, which a run keeps off the strands only if it crosses no more
    for pts in loops:
        crossed.add(pts)
    fixed: dict[tuple[int, int], NAPath] = {}  # (strand, nucleotide before it) -> a path of fixed shape

    def fix(key, path: NAPath) -> None:
        fixed[key] = path
        strands.add(among(path, at), fixed=path.kind in ("hairpin", "tail"))

    between: list[tuple] = []  # strands from stem to stem still to route, once the fixed shapes are known
    circles = []  # every hairpin loop's circle, with its stem: no free end is drawn into another stem's (_tail)
    for s, strand in enumerate(na.strands):
        held = [i for i in range(len(strand)) if (s, i) in at]
        for i, j in zip(held, held[1:]):
            loop = at[(s, i)]
            if loop is at[(s, j)] and loop.side[(s, i)] != loop.side[(s, j)]:
                circles.append((*loop.hairpin_circle((s, i), (s, j), j - i - 1), loop))
    goes_on = {}  # strand whose 3' end runs into a gap -> the stem point where the molecule goes on past it (_tail)
    for g in na.gaps:
        nxt = [i for i in range(len(na.strands[g.after])) if (g.after, i) in at]
        if nxt:
            goes_on[g.before] = tuple(float(c) for c in at[(g.after, nxt[0])].point((g.after, nxt[0])))
    for s, strand in enumerate(na.strands):
        n = len(strand)
        if s in layout.na_traces:
            pts = _xy(layout.na_traces[s])
            trace = NAPath("trace", s, list(range(n)), pts, pts, fallback=not _trace_clear(layout, s, pts))
            fix((s, -1), trace)
            for label, a, b in (("5′", pts[0], pts[1]), ("3′", pts[-1], pts[-2])):
                v = np.subtract(a, b)
                v = v / max(float(np.hypot(*v)), 1e-9)
                end(s, label, tuple(float(c) for c in np.add(a, v * (NT_ROW + glyph(layout)[0]))))
                exits[(s, 0 if label == "5′" else n - 1)] = (np.add(a, v * NT_STEP / 2), v, None)
            continue
        held = [i for i in range(n) if (s, i) in at]
        if not held:
            continue
        first, last = held[0], held[-1]
        d0 = at[(s, first)]
        if first > 0:
            path, tip, v = _tail(
                router, layout, d0, s, (s, first), list(range(first)), True, room(s, "5′"), strands, circles
            )
            fix((s, -1), path)
            exits[(s, 0)] = (np.add(path.nt_xy[0], v * NT_STEP / 2), v, d0)
        else:
            tip = _face_label(router, layout, d0, (s, 0), False, room(s, "5′"))
            exits[(s, 0)] = (*_face(d0, (s, 0), False), d0)
        end(s, "5′", tip)
        for i, j in zip(held, held[1:]):
            da, db = at[(s, i)], at[(s, j)]
            if da is db and da.side[(s, i)] == da.side[(s, j)]:
                continue  # along one rail of a stem (a bulge is drawn on the rail)
            if da is db:  # a hairpin loop: an arc over the tip
                fix((s, i), _strand_between(router, s, i, j, da, db, None))
            else:
                between.append((s, i, j, da, db, kind_of.get((s, i + 1))))
        dn = at[(s, last)]
        if last < n - 1:
            path, tip, v = _tail(
                router,
                layout,
                dn,
                s,
                (s, last),
                list(range(last + 1, n)),
                False,
                room(s, "3′"),
                strands,
                circles,
                goes_on.get(s),
            )
            fix((s, last), path)
            exits[(s, n - 1)] = (np.add(path.nt_xy[-1], v * NT_STEP / 2), v, dn)
        else:
            tip = _face_label(router, layout, dn, (s, last), True, room(s, "3′"))
            exits[(s, n - 1)] = (*_face(dn, (s, last), True), dn)
        end(s, "3′", tip)
    ends = _own_ends(router, layout, fixed, at, ends, gapped, loops=loops)
    labels[:] = [xy for s, e, xy in ends if (s, e) not in gapped]
    # routes keep off the strands' end labels (as a protein loop keeps off its chain's N and C)
    g = glyph(layout)[1] / 2 + 0.05  # sized from the text drawn (C4)
    router.grid.keep_off([(x - g, y - g, x + g, y + g) for s, e, (x, y) in ends if (s, e) not in gapped])
    # and pay to cross a hairpin loop or a free end, as a loop pays to share an earlier loop's cells, or to run
    # beside one; beside a protein loop costs the same, so crossing one costs as much as crossing a strand
    beside, beside_loops = np.zeros_like(router.grid.used), np.zeros_like(router.grid.used)
    for p in fixed.values():
        if p.kind in ("hairpin", "tail"):
            _use(router.grid, p.points, beside)
    for pts in loops:
        _use(router.grid, pts, beside_loops)
    names = [stem_keep_box(d, glyph(layout)[0] / _GLYPH_W) for d in layout.stems]  # each name with its margin
    # and each domain's name with its room (label_boxes): a run's bow, not routed on the grid, kept off it only so
    # (random-100 C2b-4: 7Y7C's domains view, A:1256-1257 bowed along the top of "D1 · B")
    names += [box for name, box in label_boxes(layout, []) if name.startswith("domain-label:")]
    routed = {}
    for s, i, j, da, db, kind in between:
        # a run with nucleotides may pay for the cells beside strands and protein loops too, one without
        # nucleotides for those beside protein loops
        extra = beside + beside_loops if j - i > 1 else beside_loops
        p = routed[(s, i)] = _strand_between(router, s, i, j, da, db, kind, names, strands, crossed, extra)
        strands.add(among(p, at), frozenset((da.id, db.id)) if p.kind == "joint" else None)
        _use(router.grid, among(p, at), beside)
    out = [p for _, p in sorted({**fixed, **routed}.items())]
    lines = [(p.strand, p.points) for _, p in sorted(routed.items())]
    ends = _own_ends(router, layout, {**fixed, **routed}, at, ends, gapped, lines, loops)  # again, the runs drawn too
    for g in na.gaps:  # cropped away: the 3' end of one piece to the 5' end of the next, dashed
        a, b = exits.get((g.before, len(na.strands[g.before]) - 1)), exits.get((g.after, 0))
        if a is not None and b is not None:
            out.append(_gap(router, g, a, b))
        elif a is not None or b is not None:  # one end drawn only (a piece too short to draw): a short break mark
            p, e, _ = a if a is not None else b
            mark = _xy([p, np.add(p, np.multiply(e, NT_STEP))])
            label = f"… {g.count} nt …" if g.count >= SHORT_GAP else ""
            out.append(NAPath("gap", g.before, [], mark, [], label))
    _unbunch(layout, out)
    layout.na_ends = [e for e in ends if (e[0], e[1]) not in gapped]
    layout.na_letters = _place_letters(layout, out)  # blind to the loops: _letters_off_loops
    layout.na_ties = _ties(layout, out, loops)
    _label_gaps(layout, out, loops)
    _widen(layout, out)
    if log.isEnabledFor(logging.DEBUG):
        log.debug("%d nucleic-acid paths: %d fallbacks", len(out), sum(p.fallback for p in out))
    return out
