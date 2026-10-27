# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The readability floor (plan 28 Task 4, the review's 3.3, rulings Q39-Q41): a figure that cannot be read as a
publication figure says so, and names the options that exist today which draw less.

readability_floor(drawn) gives the reasons a drawn figure (spec.Drawn) is below the floor, each named (REASONS), with
what was measured and the limit; [] when it is above it. Two are judged:

- text_small: the figure's text is under layout.MIN_TEXT_PT (6 pt) on its page. Below that size the page grows
  instead (layout.page_scale, random-100 C4) until it reaches layout.MAX_PAGE_IN; past that the whole figure, its text
  with it, is scaled to fit (layout.page_fit's note; page_text_small among the named fallbacks). The size judged is the
  figure's own text (element names, chain ends) at Style.font_scale 1, as page_fit judges it: a theme's smaller font
  (the journal theme's 0.9, --set font_scale=0.8) prints smaller text than the size judged. The smallest text drawn is
  given with it but never judged: residue counts, base letters and ligand names are drawn at 0.6-0.9 of the figure's
  text by design, so on any page they are under 6 pt whenever the figure's text is under 10 pt
  (docs/validation/readability-floor: 17 of the 18 figures measured).
- route_budget: the page's budget for second, wider route searches (route._MAX_NODES_PAGE) was spent, so searches it
  would have run were not (counted as the named fallbacks route_second_search_refused): the loops left are drawn as
  curves or over text. A figure's route fallbacks of other kinds (one stuck loop) are A3's, not the floor's: they
  happen on small pages too, while a spent budget is the page's size.

The time a figure took is said, never judged (ruling Q40): time_notes(measures) names a figure that took longer than
the hosted app allows (criteria.LIMITS["app_s"], 90 s) to make, measured on the drawing in hand (Drawn.seconds,
spec.draw_spec's wall time), which is the machine's and its load's as much as the figure's. It never fails R1, never
changes `foldmap check`'s exit code or --strict-floor, and is never in the check's --json report (that report holds
no wall time, so it is the same on every run).

suggestions(drawn) names the existing options that draw less, those that apply: --focus protomer (one subunit of a
symmetric assembly drawn whole, in full), --focus CHAIN (its first chain in full, with what its neighbours add to its
fold) and --chains CHAIN (that chain only). Measured on 9FZL, 8XJV and the largest entries of the three random samples
(docs/validation/readability-floor, scripts/readability_floor.py): each brought every figure it was offered for above
the floor; the notes still say "may", as nothing is measured before it is drawn. Each surface names only the options
it has (_phrase, _SAID): the app has no control for one subunit, so it is never told to the app's user.

Where it is said (said(), guarded: a failure to judge is said, never raised): the command line's notes (`foldmap
plot`; with --strict-floor no figure below the floor is written and it exits 1), `foldmap check` (R1, a readability
check, in the text and the --json report; a fail exits 1), the MCP result (its note and `floor`) and the app (a
warning under the figure; never on a compare or diff page). The note never reaches the drawn figure or its files.
The judgment does in one case: focus auto draws an assembly of two subunits whole, and draws it again as one subunit
when the whole is below the floor (spec.draw_spec, cli._focus_chains; plan 29 Task 5, 31FD)."""

from __future__ import annotations

import math

from .criteria import LIMITS
from .layout import MIN_TEXT_PT

CODE = "R1"  # the floor's check id in `foldmap check` (deep.report: group "floor", kind readability)
TITLE = "readability floor: text at 6 pt or more on its page, route budget not spent"
TEXT_SMALL, ROUTE_BUDGET, TIME = "text_small", "route_budget", "time"
REASONS = {  # judged: a figure with one is below the floor
    TEXT_SMALL: "the figure's text is under 6 pt on its page (layout.MIN_TEXT_PT): its page is at its largest",
    ROUTE_BUDGET: "the page's budget for wider route searches was spent (route._MAX_NODES_PAGE)",
}
NOTES = {  # said, never judged (ruling Q40)
    TIME: "it took longer to make than the hosted app allows (criteria.LIMITS['app_s']); the machine's as much as "
    "the figure's",
}
TIME_LIMIT_S = LIMITS["app_s"]  # the hosted app's limit (app_jobs.HOSTED)


def text_sizes(fig, walk: bool = True) -> tuple[float | None, float | None]:
    """(the figure's text size, the smallest text drawn), points: render's foldmap_text_pt (layout.font_pt of its
    page's scale, at Style.font_scale 1) and the smallest font size of the visible texts on it; None where there is
    none (a figure render did not make), and for the smallest without `walk` (no text on the figure is walked)."""
    from matplotlib.text import Text

    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()] if walk else []
    own = getattr(fig, "foldmap_text_pt", None)
    return (None if own is None else float(own)), (float(min(sizes)) if sizes else None)


def measures(drawn, walk: bool = True) -> dict:
    """What the floor is judged on: the page's size (inches), the figure's text and its smallest text (points; None
    without `walk`: text_sizes) and the named fallbacks taken (Drawn.fallbacks); and the seconds it took (None when
    not timed), for time_notes only."""
    w, h = (float(x) for x in drawn.figure.get_size_inches())
    text, smallest = text_sizes(drawn.figure, walk)
    seconds = getattr(drawn, "seconds", None)
    return {
        "page_in": [round(w, 2), round(h, 2)],
        "text_pt": None if text is None else round(text, 3),
        "smallest_pt": None if smallest is None else round(smallest, 3),
        "fallbacks": dict(drawn.fallbacks),
        "seconds": None if seconds is None else round(float(seconds), 2),
    }


def reasons_of(got: dict) -> list[dict]:
    """The floor's judged reasons for measures `got` (measures()), each {reason, measured, limit, text}; [] above
    it. Never time (time_notes)."""
    out = []
    text = got.get("text_pt")
    if text is not None and text < MIN_TEXT_PT - 1e-6:
        w, h = got["page_in"]
        said = f"its text is {_down(text)} pt, under {MIN_TEXT_PT:g} pt, on its {w:.0f} x {h:.0f} in page"
        if got.get("smallest_pt") is not None:
            said += f" (its smallest text {_down(got['smallest_pt'])} pt)"
        out.append({"reason": TEXT_SMALL, "measured": text, "limit": MIN_TEXT_PT, "text": said})
    fired = got.get("fallbacks") or {}
    if refused := fired.get("route_second_search_refused", 0):
        said = f"its route budget was spent: {refused} wider route search{'' if refused == 1 else 'es'} not run"
        if curves := fired.get("loop_curve", 0):
            said += f", {curves} loop{'' if curves == 1 else 's'} drawn as curves"
        out.append({"reason": ROUTE_BUDGET, "measured": refused, "limit": 0, "text": said})
    return out


def time_notes(got: dict) -> list[dict]:
    """The time note for measures `got`, [{reason: time, measured, limit, text}] when it took longer than the hosted
    app allows, else []: said, never judged (ruling Q40)."""
    seconds = got.get("seconds")
    if seconds is None or seconds <= TIME_LIMIT_S:
        return []
    said = f"it took {seconds:.0f} s to make, over the hosted app's {TIME_LIMIT_S:g} s (said, not judged)"
    return [{"reason": TIME, "measured": seconds, "limit": TIME_LIMIT_S, "text": said}]


def _down(pt: float) -> str:
    """A size in points rounded down to a tenth, as page_fit says it (5.97 pt is 5.9, under 6)."""
    return f"{math.floor(pt * 10) / 10:.1f}"


def readability_floor(drawn) -> list[dict]:
    """The named reasons the figure `drawn` (spec.Drawn) is below the readability floor (module docstring), each
    {reason, measured, limit, text}; [] when it is above it."""
    return reasons_of(measures(drawn))


def below(drawn) -> bool:
    """Whether the figure `drawn` is below the floor (readability_floor gives a reason), judged on the same measures
    without its smallest text (said, never judged), so no text on the figure is walked (spec.draw_spec asks it of
    every assembly of two that focus auto draws whole)."""
    return bool(reasons_of(measures(drawn, walk=False)))


def suggestions(drawn) -> list[dict]:
    """The existing options that draw less of the figure `drawn`, each {options: {name: value}, what}: one subunit
    in full when it draws a symmetric assembly whole (--focus protomer), then its first chain in full and that chain
    alone (--focus CHAIN, --chains CHAIN), while it draws more than that chain; [] when none applies (one chain
    drawn)."""
    lay = drawn.layout
    focus = sorted(lay.focus_chains or [])
    sector = getattr(lay.frame, "sector", None) or {}
    strands = lay.nucleic.strands if lay.nucleic is not None else []
    protein = list(dict.fromkeys(label.chain for label in drawn.bb.labels))
    drawn_chains = set(protein) | {s.chain for s in strands}
    out = []
    if len(set(sector.values())) > 1 and not focus:
        out.append({"options": {"focus": "protomer"}, "what": "one subunit in full"})
    first = protein[0] if protein else (strands[0].chain if strands else None)
    if first is None or drawn_chains == {first}:
        return out
    if protein and focus != [first]:
        out.append(
            {"options": {"focus": first}, "what": f"chain {first} in full, with what its neighbours add to its fold"}
        )
    out.append({"options": {"chains": first}, "what": f"chain {first} only"})
    return out


def _phrase(options: dict, how: str) -> str | None:
    """Options in the words of the surface showing the note (_SAID): the command line (how='cli'), the MCP tool
    (how='mcp') or the app (how='app'); None when that surface has no control for one of them, so the note never
    tells a user to do what they cannot (the app has no subunit control: Large assemblies offers auto, none and
    chains only)."""
    said = _SAID[how]
    parts = []
    for k, v in options.items():
        word = said.get(f"{k}:{v}") if k == "focus" and v in _FOCUS_WORDS else said.get(k)  # a word: its own entry
        if word is None:
            return None
        parts.append(word.format(v=v))
    return " ".join(parts)


_FOCUS_WORDS = ("protomer",)  # focus values that are words, not chains: a surface names each only by its own entry
_SAID = {  # each surface's words for the options it has (suggestions' options); an option missing here is not named
    "cli": {"focus": "--focus {v}", "focus:protomer": "--focus {v}", "chains": "--chains {v}"},
    "mcp": {"focus": 'focus="{v}"', "focus:protomer": 'focus="{v}"', "chains": 'chains=["{v}"]'},
    "app": {  # the app's controls (app_controls.chains_section; Large assemblies is app_controls.FOCUS): no subunit
        "focus": "Chains & assembly, Large assemblies: these chains in full, {v}",
        "chains": "Chains & assembly, Chains: {v}",
    },
}


def ways(offered: list[dict], how: str = "cli") -> str:
    """The options that draw less (suggestions()) as one clause, in `how`'s words (cli, mcp or app, as _phrase); ''
    when none applies."""
    said = [f"{p} ({o['what']})" for o in offered if (p := _phrase(o["options"], how))]
    return "drawing less may bring it above: " + "; ".join(said) if said else ""


def note(reasons: list[dict], offered: list[dict], how: str = "cli", timed: list[dict] = ()) -> str:
    """One line for the notes: below the readability floor, why, and the options that draw less and may bring it
    above (ways); then the time note (`timed`, time_notes), said with or without a reason; '' when neither."""
    parts = []
    if reasons:
        parts.append("below the readability floor: " + "; ".join(r["text"] for r in reasons))
        if more := ways(offered, how):
            parts.append(more)
    parts += [t["text"] for t in timed]
    return "; ".join(parts)


def said(drawn, how: str = "cli", strict: bool = False) -> tuple[list[dict], str]:
    """(the judged reasons, the notes' line: note() with the time note) for the figure `drawn`, in `how`'s words;
    guarded: a failure to judge it is said ("readability floor not judged: ..."), never raised, so a figure is
    always written whatever the floor's code does. `strict` (`--strict-floor`): the failure is raised instead, so a
    figure the floor could not judge is not written as if it were above it (plan 28 final review F6)."""
    try:
        got = measures(drawn)
        reasons = reasons_of(got)
        return reasons, note(reasons, suggestions(drawn) if reasons else [], how, time_notes(got))
    except Exception as err:  # noqa: BLE001 (said in the notes: the floor never stops a figure, unless strict)
        if strict:
            raise
        return [], f"readability floor not judged: {type(err).__name__}"


def check(drawn) -> dict:
    """R1 for `foldmap check`: {status, detail, reasons, suggestions, measures} for deep.report (status pass or fail
    on the judged reasons; no wall time anywhere in it), and `timed` (time_notes) for the text only."""
    got = measures(drawn)
    reasons = reasons_of(got)
    offered = suggestions(drawn) if reasons else []
    return {
        "status": "fail" if reasons else "pass",
        "detail": note(reasons, offered) if reasons else "above the readability floor",
        "reasons": reasons,
        "suggestions": offered,
        "measures": {k: v for k, v in got.items() if k != "seconds"},
        "timed": time_notes(got),
    }


def lines(got: dict) -> list[str]:
    """`foldmap check`'s lines for R1 (check()): the verdict, then each reason and the options, then the time note,
    indented."""
    out = [f"{CODE} {TITLE}: {got['status']}"]
    if got["status"] == "fail":
        out += [f"  {r['text']}" for r in got["reasons"]]
        if more := ways(got["suggestions"]):
            out.append(f"  {more}")
    out += [f"  {t['text']}" for t in got.get("timed", [])]
    return out
