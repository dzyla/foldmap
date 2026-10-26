# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Named fallbacks, counted per figure (plan 28 Task 9.1, ruling Q23).

A fallback is a path the code takes when its usual way gives no result, and whose result is worse than the usual
one: a loop no route was found for drawn as a curve, a name with no clear spot put at the least bad one, a count
with no clear spot left out. Each has a name in FALLBACKS; the place it is taken calls note(name). Inside
counting() (spec.draw_spec opens one for each figure) every note is counted; outside one, a note does nothing. The
counts go to Layout.fallbacks and from there to the notes the command line prints, the MCP result and the
random-100 report, so a campaign says how often each fallback fired. They never reach the drawn figure (its SVG,
PDF or PNG): counting changes no figure.

A count is how many times the path was taken while the figure was made, including in drawings tried and set aside
(a candidate layout make_layout weighs, a reference laid out for --like). Figures made outside spec.draw_spec count
nothing: compare and diff panels (compare.py), interactive pages (interactive.build_page, app_jobs.page_job) and
cli.make_figure_and_loops."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from contextvars import ContextVar

# name -> what was done instead, and where (module.function); scripts/inventory.py lists every fallback.note call and
# docs/validation/inventory/inventory.md says which check, if any, reports each
FALLBACKS: dict[str, str] = {
    # routing (route.py)
    "loop_curve": "a loop no orthogonal route was found for, drawn as a smooth curve (route.route_loops; counted over "
    "the loops as drawn, a rosette's turned copies included)",
    "loop_over_text": "a loop whose only route runs over a label or a chain end (route._Router.route, its keep-outs "
    "opened; counted over the loops as drawn, once the names have moved off lines)",
    "loop_wall_detour": "a membrane view's fallback curve that would cross its loop's wall, drawn instead as a detour "
    "out of its ports and round the page (route._walled_fallback, ruling P35)",
    "loop_wall_unkept": "... and no detour keeps to the wall either: the last one tried is drawn, logged "
    "(route._walled_fallback)",
    "route_unsnapped": "a route whose drawing snapped to its ports crossed an element, drawn unsnapped "
    "(route._Router.drawn; every drawing counted, those of searches set aside included)",
    "route_second_search_refused": "a wider second route search not run because the page's node budget was spent "
    "(route._Router.second, _MAX_NODES_PAGE)",
    "helix_name_moved": "a helix name a line runs through or near, moved to the spot lines run nearest least "
    "(route._clear_helix_names, the fallback for lines no route holds)",
    # nucleic acid (counted in route_na.route_nucleic over the paths as drawn)
    "na_tail_straight": "a free nucleic-acid end with no clear drawing, drawn straight on (na_ends._tail)",
    "na_gap_unclear": "a cropped gap's dashed line with no drawing clear of the elements (na_ends._gap)",
    "na_hairpin_unclear": "a hairpin loop's arc not clear of the elements (na_runs._strand_between)",
    "na_run_unclear": "a single strand or joint between stems with no clear route, drawn as a curve (na_runs._drawing)",
    "na_trace_unclear": "an unpaired strand's trace not clear of the elements, a nucleotide no nudge clears kept "
    "(route_na.route_nucleic, na_ends.off_unheld)",
    # layout (layout.py)
    "stem_name_least_bad": "a stem's name with no clear slot nearer its own stem with room, put in the slot with the "
    "most room (layout._choose_label_sides, _pick_slot; logged)",
    "domain_name_least_bad": "a domain's name with no clear place, put where it leaves the loops' ways out the most "
    "room (layout._place_domain_names)",
    "glyph_width_guessed": "a text measured with a character the built-in metrics lack, its width guessed "
    "(layout._text_width and _ink_line: _FALLBACK_WIDTH; render._clear_of_lines: half an em); every measuring counted",
    # drawing (render.py)
    "mark_least_bad_spot": "a ligand mark with no clear spot, drawn over other ink at the least bad one "
    "(render._least_bad_spot; counted over the marks drawn)",
    "counts_hidden": "a residue count of a long loop or chain end with no clear spot, left out (render.tail_counts, "
    "loop_counts)",
    "interactions_dropped": "an interaction line the independent re-check (interactions.verify_interactions) found "
    "breaking its rule, not drawn",
    "interactions_hidden": "an interaction line that would run wholly under text, not drawn "
    "(render._draw_interactions)",
    "crosslinks_unplaced": "a crosslink with an end not drawn or both ends at one point, not drawn (render)",
    "page_text_small": "a page held at its largest size, its text under 6 pt (layout.page_fit, render._canvas)",
    # checks (science.py): counted only when the check runs inside counting()
    "s3_cycles_capped": "S3 stopped counting a sheet's pairing cycles at _MAX_CYCLES: the barrel check says there are "
    "too many to count (science.pairing_cycle)",
}

_FIRED: ContextVar[Counter | None] = ContextVar("foldmap_fallbacks", default=None)


def note(name: str, n: int = 1) -> None:
    """Fallback `name` (a FALLBACKS key) was taken `n` times: counted in the figure being made, if one is."""
    fired = _FIRED.get()
    if fired is not None and n:
        fired[name] += n


@contextmanager
def counting() -> Iterator[Counter]:
    """Count the fallbacks taken inside: yields the Counter they are counted in. Nested, the inner counts also go
    to the outer one when the inner one closes."""
    outer = _FIRED.get()
    fired: Counter = Counter()
    token = _FIRED.set(fired)
    try:
        yield fired
    finally:
        _FIRED.reset(token)
        if outer is not None:
            outer.update(fired)


def counts(fired) -> dict[str, int]:
    """The counts as a plain dict, by name, without zeros."""
    return {k: int(fired[k]) for k in sorted(fired) if fired[k]}


def line(fired) -> str:
    """The notes' line ("fallbacks taken: loop_curve 3, ..."), "" when none was taken."""
    got = counts(fired or {})
    return "fallbacks taken: " + ", ".join(f"{k} {n}" for k, n in got.items()) if got else ""
