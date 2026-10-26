# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The readability floor (plan 28 Task 4, the review's 3.3, ruling Q39): a figure that cannot be read as a
publication figure says so, and names the options that exist today which draw less.

readability_floor(drawn) gives the reasons a drawn figure (spec.Drawn) is below the floor, each named (REASONS), with
what was measured and the limit; [] when it is above it:

- text_small: the figure's text is under layout.MIN_TEXT_PT (6 pt) on its page. Below that size the page grows
  instead (layout.page_scale, random-100 C4) until it reaches layout.MAX_PAGE_IN; past that the whole figure, its text
  with it, is scaled to fit (layout.page_fit's note; page_text_small among the named fallbacks). The size judged is the
  figure's own text (element names, chain ends) at Style.font_scale 1, as page_fit judges it: a theme's smaller font
  is the theme's choice, not the page's. The smallest text drawn is given with it but never judged: residue counts,
  base letters and ligand names are drawn at 0.6-0.9 of the figure's text by design, so on any page they are under
  6 pt whenever the figure's text is under 10 pt (docs/validation/readability-floor: on every figure measured).
- route_budget: the page's budget for second, wider route searches (route._MAX_NODES_PAGE) was spent, so searches it
  would have run were not (counted as the named fallbacks route_second_search_refused): the loops left are drawn as
  curves or over text. A figure's route fallbacks of other kinds (one stuck loop) are A3's, not the floor's: they
  happen on small pages too, while a spent budget is the page's size.
- time: making the figure took longer than the hosted app allows (criteria.LIMITS["app_s"], 90 s). Measured on the
  drawing in hand (Drawn.seconds, spec.draw_spec's wall time), so it depends on the machine and its load; a Drawn
  with no time (None: a compare page, a figure not made by draw_spec) is not judged on it.

suggestions(drawn) names the existing options that draw less, those that apply: --focus protomer (one subunit of a
symmetric assembly drawn whole, in full), --focus CHAIN (its first chain in full, with what its neighbours add to its
fold) and --chains CHAIN (that chain only). Measured on 9FZL, 8XJV and the largest entries of the three random samples
(docs/validation/readability-floor, scripts/readability_floor.py): each brought every figure it was offered for above
the floor; the notes still say "may", as nothing is measured before it is drawn. --reduced-detail (REDUCED_DETAIL, an
option preset of existing style keys that leaves out low-priority annotations: base letters, ligand marks, glycans,
3-10 helix boxes; elements, chains and their interfaces kept) removes most of the text drawn under 6 pt but measured
to bring no figure above the floor (the page's size is set by its elements), so the notes do not offer it for that.

Where it is said: the command line's notes (`foldmap plot`; with --strict-floor no figure is written below the floor
and it exits 1), `foldmap check` (R1, a readability check, in the text and the --json report), the MCP result (its
note and `floor`) and the app (a warning under the figure). It never reaches the drawn figure or its files: the
default figure is unchanged (bench --exact)."""

from __future__ import annotations

import math

from .criteria import LIMITS
from .layout import MIN_TEXT_PT

CODE = "R1"  # the floor's check id in `foldmap check` (deep.report: group "floor", kind readability)
TITLE = "readability floor: text at 6 pt or more, route budget not spent, made within the hosted time"
TEXT_SMALL, ROUTE_BUDGET, TIME = "text_small", "route_budget", "time"
REASONS = {
    TEXT_SMALL: "the figure's text is under 6 pt on its page (layout.MIN_TEXT_PT): its page is at its largest",
    ROUTE_BUDGET: "the page's budget for wider route searches was spent (route._MAX_NODES_PAGE)",
    TIME: "it took longer to make than the hosted app allows (criteria.LIMITS['app_s'])",
}
TIME_LIMIT_S = LIMITS["app_s"]  # the hosted app's limit (app_jobs.HOSTED)
# the reduced-detail preset (--reduced-detail, the MCP's reduced_detail): existing style keys only, low-priority
# annotations off; given before the user's own style keys, which win
REDUCED_DETAIL = ("nucleotide_labels=false", "ligands=none", "glycans=false", "helices_310=false")


def text_sizes(fig) -> tuple[float | None, float | None]:
    """(the figure's text size, the smallest text drawn), points: render's foldmap_text_pt (layout.font_pt of its
    page's scale, at Style.font_scale 1) and the smallest font size of the visible texts on it; None where there is
    none (a figure render did not make)."""
    from matplotlib.text import Text

    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    own = getattr(fig, "foldmap_text_pt", None)
    return (None if own is None else float(own)), (float(min(sizes)) if sizes else None)


def measures(drawn) -> dict:
    """What the floor is judged on: the page's size (inches), the figure's text and its smallest text (points), the
    named fallbacks taken (Drawn.fallbacks) and the seconds it took (None when not timed)."""
    w, h = (float(x) for x in drawn.figure.get_size_inches())
    text, smallest = text_sizes(drawn.figure)
    seconds = getattr(drawn, "seconds", None)
    return {
        "page_in": [round(w, 2), round(h, 2)],
        "text_pt": None if text is None else round(text, 3),
        "smallest_pt": None if smallest is None else round(smallest, 3),
        "fallbacks": dict(drawn.fallbacks),
        "seconds": None if seconds is None else round(float(seconds), 2),
    }


def reasons_of(got: dict) -> list[dict]:
    """The floor's reasons for measures `got` (measures()), each {reason, measured, limit, text}; [] above it."""
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
    seconds = got.get("seconds")
    if seconds is not None and seconds > TIME_LIMIT_S:
        said = f"it took {seconds:.0f} s to make, over the hosted app's {TIME_LIMIT_S:g} s"
        out.append({"reason": TIME, "measured": seconds, "limit": TIME_LIMIT_S, "text": said})
    return out


def _down(pt: float) -> str:
    """A size in points rounded down to a tenth, as page_fit says it (5.97 pt is 5.9, under 6)."""
    return f"{math.floor(pt * 10) / 10:.1f}"


def readability_floor(drawn) -> list[dict]:
    """The named reasons the figure `drawn` (spec.Drawn) is below the readability floor (module docstring), each
    {reason, measured, limit, text}; [] when it is above it."""
    return reasons_of(measures(drawn))


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
    """Options as the command line (how='cli'), the MCP tool (how='mcp') or the app (how='app') takes them; None
    where the app has no control for them (--focus protomer)."""
    parts = []
    for k, v in options.items():
        if how == "mcp":
            parts.append(f'chains=["{v}"]' if k == "chains" else f'{k}="{v}"')
        elif how == "app":
            if k == "focus" and v == "protomer":
                return None
            parts.append(_APP[k].format(v=v))
        else:
            parts.append(f"--{k} {v}")
    return " ".join(parts)


_APP = {  # the app's controls for each option (app_controls.chains_section)
    "focus": "Chains & assembly, Large assemblies: these chains in full, {v}",
    "chains": "Chains & assembly, Chains: {v}",
}


def ways(offered: list[dict], how: str = "cli") -> str:
    """The options that draw less (suggestions()) as one clause, in `how`'s words (cli, mcp or app, as _phrase); ''
    when none applies."""
    said = [f"{p} ({o['what']})" for o in offered if (p := _phrase(o["options"], how))]
    return "drawing less may bring it above: " + "; ".join(said) if said else ""


def note(reasons: list[dict], offered: list[dict], how: str = "cli") -> str:
    """One line for the notes: below the readability floor, why, and the options that draw less and may bring it
    above (ways); '' above it."""
    if not reasons:
        return ""
    said = "below the readability floor: " + "; ".join(r["text"] for r in reasons)
    more = ways(offered, how)
    return f"{said}; {more}" if more else said


def check(drawn) -> dict:
    """R1 for `foldmap check` (deep.report's floor entry): {status, detail, reasons, suggestions, measures}; status
    pass or fail, the reasons named."""
    got = measures(drawn)
    reasons = reasons_of(got)
    offered = suggestions(drawn) if reasons else []
    detail = note(reasons, offered) if reasons else "above the readability floor"
    return {
        "status": "fail" if reasons else "pass",
        "detail": detail,
        "reasons": reasons,
        "suggestions": offered,
        "measures": got,
    }


def lines(got: dict) -> list[str]:
    """`foldmap check`'s lines for R1 (check()): the verdict, then each reason and the options, indented."""
    out = [f"{CODE} {TITLE}: {got['status']}"]
    if got["status"] == "fail":
        out += [f"  {r['text']}" for r in got["reasons"]]
        if said := ways(got["suggestions"]):
            out.append(f"  {said}")
    return out
