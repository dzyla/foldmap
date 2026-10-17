# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Figure criteria A1-A20 (plan 10, "Pass criteria"): is a drawn figure fit to publish as it is? They judge how the
page reads and what it leaves out, apart from the science checks S1-S6 (science.py, figcheck.science_checks), which
judge what it says about the structure. The random-100 campaign (scripts/random100.py) judges every entry by them over
all its views, each drawn in a process of its own; `foldmap check` judges the one figure it draws (check_figure).

- measure_figure: what one drawn figure gives to be judged on (figcheck.figure_checks, route fallbacks, loops over
  labels, scope, tethers, stem-name limit, membrane, science). The campaign adds what it measures apart (layout
  signature and SVG hashes, mkdssp's codes, RNApolis base pairs) and passes it in.
- criteria: A1-A20 of one entry from its views' results ({"pass": True, False or None, "detail"}, and for None its
  "status": NOT_JUDGED, the criterion applies but the entry gave it nothing to judge, or NOT_APPLICABLE);
  figure_criteria: those judged on the drawn figures alone; pair_check: A19 of a two-state pair; outcome: one result
  as a word. An entry passes when no criterion fails: nothing to judge is never counted as a pass (final review B,
  Important 3).
- check_figure, criteria_lines: one figure's criteria as `foldmap check` prints them; the criteria one figure in one
  process cannot give are named not measured here, with the reason (NOT_HERE), never passed."""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import numpy as np

LIMITS = {
    "view_timeout_s": 900,  # a view's process is killed after this many seconds (A1 then fails)
    "default_s": 60,  # A2: wall seconds of the default view
    "default_s_many_chains": 180,  # ... when more than many_chains polymer chains are read (Q-R3)
    "many_chains": 20,
    "app_s": 90,  # A2b: the hosted app's limits (app_jobs.HOSTED, Q-R7)
    "app_memory_mb": 1024,
    "app_max_residues": 30_000,
    "dssp_agree": 0.99,  # A10: fraction of residues whose SS agrees with mkdssp 3.1.2
    "copies_alike": 0.5,  # A15: page units (Q-R12)
    "domain_overlap": 0.02,  # A16: fraction of the smaller panel
    # A17: page units from a stem name to its own letters, at most layout._INK_MARGIN em in the drawn figure's page
    # units plus 1e-3 (figcheck.stem_gap_tolerance; plan 10 Task 4 fix round 2 ruling), measured per figure
    "stem_name_gap": "figcheck.stem_gap_tolerance(layout._INK_MARGIN em)",
    "tether_length": 6.5,  # A20: page units, a ligand tether from its mark to another chain than its home
}
DOCUMENTED_WARNINGS: tuple[str, ...] = ()  # Python warnings a view may print (regexes); none at the start
FIGURE_VIEWS = ("default", "deposited", "beta", "focus", "domains", "membrane")
LABEL_VIEWS = ("default", "deposited", "beta", "focus")
CRITERIA = {
    "A1": "no crash: exit 0, files written, no traceback, no Python warning",
    "A1c": "a CA-only or P-only entry exits non-zero with foldmap's one-line message",
    "A2": "default view time",
    "A2b": "the app draws it within the hosted limits",
    "A3": "no route fallbacks",
    "A4": "no loop over a label",
    "A5": "no element or stem boxes overlapping",
    "A6": "no text collisions",
    "A7": "text and legend inside the page",
    "A8": "labels unique per chain, stem names per figure",
    "A9": "nothing left out without a note",
    "A10": "secondary structure agrees with mkdssp",
    "A11": "N->C and sheet order",
    "A12": "same figure for any PYTHONHASHSEED",
    "A13": "a saved layout redraws the same figure",
    "A14": "interactive page and sequence plot written",
    "A15": "copies drawn alike",
    "A16": "domains: one chain, alike, not tiny, panels apart, names clear",
    "A17": "nucleic-acid metrics and stem names by their letters",
    "A18": "membrane: TM helices in the band, loops out of the slab",
    "A19": "two states: same names in the same order",
    "A20": "ligand tethers local",
}


def jsonable(x):
    """x with tuple keys joined by '~', sets and tuples as lists, numpy numbers as Python ones, NaN as None (a numpy
    NaN too: final review B, Minor 8)."""
    if isinstance(x, dict):
        return {("~".join(map(str, k)) if isinstance(k, tuple) else str(k)): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set, frozenset)):
        items = sorted(x, key=str) if isinstance(x, (set, frozenset)) else x
        return [jsonable(v) for v in items]
    if hasattr(x, "item") and not isinstance(x, (str, bytes)):
        x = x.item()
    if isinstance(x, float) and x != x:
        return None
    return x


# ------------------------------------------------------------------------------------- one figure's measures


def tether_lengths(drawn, look) -> list[tuple[str, float, str, str]]:
    """(gid, length in page units, home chain, target chain) of every ligand tether drawn: from its ligand's mark
    to its target (render's ligand_marks; the drawn tether shows only where no body or text is, so its ink is
    shorter). The home chain is the one the ligand is drawn with: render._local_targets lists its targets first."""
    from .render import ligand_marks

    lay, out = drawn.layout, []
    for m in ligand_marks(lay, drawn.loops, drawn.sses, look):
        g, (cx, cy) = m["ligand"], m["centre"]
        if not m["residues"]:  # held by nucleotides alone (C33): tethered to them, its own chain's, reported only
            chains = [g.chain] * len(m["targets"])
        else:
            chains = [lay.res_chain[r] for r in m["residues"]]
        home = chains[0]
        for k, ((tx, ty), chain) in enumerate(zip(m["targets"], chains)):
            gid, length = f"ligand-tether:{g.name}:{g.chain}{g.seq}:{k}", round(float(np.hypot(tx - cx, ty - cy)), 3)
            out.append((gid, length, home, chain))
    return out


def membrane_checks(drawn, path) -> dict:
    """A18 on a --membrane auto figure (figcheck.membrane_view, the same test the membrane view lays itself out by):
    each TM helix (its CA run reaches both membrane faces) enters its band through a face and reaches past its middle
    (C27 ruling); loops run sideways inside a band (figcheck.SLAB_MARGIN from its faces) only on the first segment
    out of their own port inside it, up to clearing its element (membrane slab ruling)."""
    from . import figcheck
    from .membrane import find_membrane

    lay, bb = drawn.layout, drawn.bb
    mem = find_membrane(path, bb)
    if lay.membrane is None or mem is None:
        return {"band": lay.membrane is not None, "tm_helices": 0, "tm_short": [], "loops_in_slab": []}
    return {"band": True, **figcheck.membrane_view(lay, drawn.sses, bb, mem, drawn.loops)}


_BLOCK_NAME_EM = 0.9  # render._draw_block draws a stem's name on its block at this many ems of the figure's text


def stem_gap_limit(fig) -> float | None:
    """A17's limit for this figure: the stem-name margin (layout._INK_MARGIN em) in its page units, plus 1e-3
    (figcheck.stem_gap_tolerance). The em is read from a drawn stem name (render: 0.8 em) as points -> pixels ->
    page units, as tests/test_na_layout.py reads it from a base letter; None when no stem name is drawn."""
    from . import figcheck
    from .layout import _INK_MARGIN, _NAME_EM

    if not fig.axes:
        return None
    ax = fig.axes[0]
    name = next((t for t in ax.texts if (t.get_gid() or "").startswith("na:label:")), None)
    if name is None:
        return None
    block = any(re.fullmatch(r"na:stem:\d+:block(:\d+)?", a.get_gid() or "") for a in ax.patches)
    px_per_unit = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    size = _BLOCK_NAME_EM if block else _NAME_EM  # render draws a block's name at 0.9 em, a ladder's at 0.8
    em = name.get_fontsize() / size * fig.dpi / 72 / px_per_unit
    return figcheck.stem_gap_tolerance(_INK_MARGIN * em)


def dssp_agreement(path: Path, rows: list[tuple[str, int, str, str]] | str) -> dict:
    """foldmap's DSSP on the model as deposited against mkdssp's codes (rows; a string: why there are none, as the
    campaign's mkdssp_rows gives it), residue by residue, as tests/test_dssp_reference.disagreements compares them: a
    residue standing in for a fused chromophore, and the four on either side of it, are not compared. agree: the
    fraction of the residues both read that agree (None, with `unmeasured` saying why, when mkdssp cannot read the
    file)."""
    from .dssp import assign_dssp
    from .io import load_backbone

    if isinstance(rows, str):
        return {"agree": None, "unmeasured": rows}
    ref = {(c, n, i): code for c, n, i, code in rows}
    bb = load_backbone(Path(path), "asu")
    ss = assign_dssp(bb).ss
    keys = [(l.chain, l.seq, l.icode) for l in bb.labels]
    stood_in = {keys[k] for k in bb.linked}
    near = {keys[j] for k in bb.linked for j in range(k - 4, k + 5) if 0 <= j < len(keys)} - stood_in
    compared = [(k, key) for k, key in enumerate(keys) if key in ref and key not in near and key not in stood_in]
    differ = [
        f"{c}{n}{i} mkdssp {ref[(c, n, i)]} foldmap {ss[k]}" for k, (c, n, i) in compared if ss[k] != ref[(c, n, i)]
    ]
    if not compared:
        return {"agree": None, "unmeasured": "no residue read by both"}
    return {
        "agree": round(1 - len(differ) / len(compared), 4),
        "residues": len(compared),
        "foldmap_only": len(set(keys) - stood_in - set(ref)),
        "mkdssp_only": len(set(ref) - set(keys)),
        "differ": differ[:20],
    }


def measure_figure(view: str, drawn, spec, path: Path) -> dict:
    """Everything one drawn figure (`view`: the campaign's name for it, VIEWS) is judged on, that the figure itself
    gives: what the campaign adds (layout signature and SVG hashes for A12 and A13, na_validate's nucleic-acid
    measures for A17) it measures in its own process."""
    from . import figcheck
    from .scope import drawn_scope

    lay = drawn.layout
    checks = figcheck.figure_checks(drawn, path, assembly=spec.assembly, focus=spec.focus, chains=spec.chains)
    out = {
        "checks": jsonable(checks),
        "stuck": drawn.stuck,
        "over_labels": drawn.over_labels,
        "crowded_marks": drawn.crowded_marks,  # reported: ligand marks with no clear spot (render._least_bad_spot)
        "scope": jsonable(drawn_scope(drawn.model, lay, spec.assembly, spec.focus)),
        "chains_read": list(lay.chains_read),
        "protein_chains_drawn": list(dict.fromkeys(p.sse.chain for p in lay.placed.values())),
        "first_protein_chain": drawn.bb.labels[0].chain if len(drawn.bb) else None,
        "na_chains_drawn": sorted({s.chain for s in lay.nucleic.strands}) if lay.nucleic is not None else [],
        "elements": [
            {
                "id": p.sse.id,
                "chain": p.sse.chain,
                "kind": p.sse.kind,
                "first": p.sse.first.seq,
                "last": p.sse.last.seq,
                "label": p.label,
            }
            for p in lay.placed.values()
        ],
        "area": round(lay.width * lay.height, 2),
        "context": lay.context,
        "copy_ties": len(lay.copy_ties),
        "tethers": tether_lengths(drawn, spec.look),
        "stem_gap_limit": stem_gap_limit(drawn.figure),
    }
    if view == "beta" and lay.placed:  # A19 pairs this chain's elements with its partner's by sequence
        chain = next(iter(lay.placed.values())).sse.chain
        out["chain_residues"] = {
            "chain": chain,
            "residues": [[x.seq, x.name] for x in drawn.bb.labels if x.chain == chain],
        }
    if view == "membrane":
        out["membrane"] = membrane_checks(drawn, path)
    if view in FIGURE_VIEWS:  # S1-S4, S6 in every drawn view, judged apart from A1-A20 (batch C2b-S)
        out["science"] = jsonable(figcheck.science_checks(drawn))
        out["science_judged"] = figcheck.science_judged(drawn)  # what each had to judge: none is no pass
    return out


# ------------------------------------------------------------------------------------------------ judging


NOT_JUDGED = "not judged"  # the criterion applies to the entry, but the entry gave it nothing to judge
NOT_APPLICABLE = "does not apply"  # the criterion has no meaning for the entry (A17 with no nucleic acid)
OUTCOMES = ("pass", "fail", NOT_JUDGED, NOT_APPLICABLE)


def _ok(passed: bool | None, detail: str = "", status: str = NOT_APPLICABLE) -> dict:
    """A criterion's result; `status` says which of the two None is (NOT_JUDGED or NOT_APPLICABLE)."""
    if passed is None:
        return {"pass": None, "detail": detail, "status": status}
    return {"pass": passed, "detail": detail}


def _not_judged(detail: str) -> dict:
    return _ok(None, detail, NOT_JUDGED)


def outcome(got: dict | None) -> str:
    """One result as a word of OUTCOMES: pass, fail, not judged or does not apply (a None with no status, as reports
    from before the final review B fix wrote it, does not apply)."""
    if not got or got.get("pass") is None:
        return (got or {}).get("status") or NOT_APPLICABLE
    return "pass" if got["pass"] else "fail"


def _views(e: dict, names) -> dict[str, dict]:
    """The measured results of these views, those run."""
    return {v: (e["views"][v].get("result") or {}) for v in names if v in e["views"]}


def _count(results: dict[str, dict], key: str, limit: float = 0) -> dict:
    """A check counted in each view's figcheck results: pass when every view's count is at most limit; a count
    that could not be measured (None) or is not a finite number fails."""
    bad = []
    for v, r in results.items():
        got = (r.get("checks") or {}).get(key, None) if "checks" in r else None
        if got is not None and not (isinstance(got, (int, float)) and np.isfinite(got)) and got != float("inf"):
            bad.append(f"{v} {key} is not a number ({got!r})")
        elif got is None:
            why = ((r.get("checks") or {}).get("skipped") or {}).get(key) or r.get("measure_error", "not measured")
            bad.append(f"{v} not measured ({str(why).strip().splitlines()[-1][:120]})")
        elif got > limit:
            details = ((r["checks"].get("details") or {}).get(key)) or []
            bad.append(f"{v} {got}" + (f" ({'; '.join(map(str, details[:3]))})" if isinstance(details, list) else ""))
    return _ok(None if not results else not bad, "; ".join(bad))


def criteria(e: dict) -> dict:
    """A1-A20 of one entry (plan 10, "Pass criteria"): {"pass": True, False or None, "detail"}, a None with its
    "status" (not judged, or does not apply: _ok)."""
    out = {k: _ok(None) for k in CRITERIA}
    if e.get("error"):
        out["A1"] = _ok(False, f"cannot read the entry: {e['error']}")
        return out
    views = e["views"]
    if e.get("refusal"):
        d = views.get("default", {})
        said = d.get("foldmap_said", [])
        ok = d.get("exit", 0) != 0 and not d.get("traceback") and len(said) == 1
        out["A1c"] = _ok(ok, "" if ok else f"exit {d.get('exit')}, said {said[:2]}, traceback {d.get('traceback')}")
        return out
    problems = []
    for v, r in views.items():
        res = r.get("result") or {}
        if r["timed_out"]:
            problems.append(f"{v}: killed after {LIMITS['view_timeout_s']} s")
        elif r["exit"] != 0:
            problems.append(f"{v}: exit {r['exit']} ({' | '.join(r['stderr_tail'][-2:])})")
        if r["traceback"]:
            problems.append(f"{v}: traceback")
        warns = res.get("warnings", []) + r.get("python_warnings", [])
        bad = [w for w in warns if not any(re.search(p, w) for p in DOCUMENTED_WARNINGS)]
        if bad:
            problems.append(f"{v}: {bad[0]}")
        empty = [n for n, size in r.get("written", {}).items() if not size]
        if empty and r["exit"] == 0:
            problems.append(f"{v}: not written: {', '.join(empty)}")
        if r["exit"] == 0 and r.get("result") is None:
            problems.append(f"{v}: no result")
        if res.get("measure_error"):
            problems.append(f"{v}: not measured: {res['measure_error'].strip().splitlines()[-1][:160]}")
    out["A1"] = _ok(not problems, "; ".join(problems))
    d = _views(e, ["default"]).get("default", {})
    if "default" in views:
        many = len(d.get("chains_read", [])) > LIMITS["many_chains"]
        limit = LIMITS["default_s_many_chains" if many else "default_s"]
        took = views["default"]["seconds"]
        out["A2"] = _ok(took <= limit, f"{took} s (limit {limit} s)")
    app = (_views(e, ["app"]).get("app") or {}).get("app")
    if "app" in views:
        if app is None:
            out["A2b"] = _ok(False, "the app view gave no answer")
        elif app["refused"] is None:
            out["A2b"] = _ok(True)
        elif app["refused"] == "size" and app["residues"] > LIMITS["app_max_residues"]:  # refused as it should be
            out["A2b"] = _ok(None, f"{app['residues']:,} residues, over the hosted cap: refused with the size message")
        else:
            out["A2b"] = _ok(False, f"{app['refused']}: {app['message']}")
    out.update(figure_criteria(views))
    dssp = e.get("dssp")
    if dssp is not None:
        if dssp["agree"] is None:  # not a fail, listed (plan): never a pass
            out["A10"] = _not_judged(f"unmeasured: {dssp['unmeasured']}")
        else:
            out["A10"] = _ok(dssp["agree"] >= LIMITS["dssp_agree"], f"{dssp['agree']:.4f} of {dssp['residues']}")
    for key, other, fields in (
        ("A12", "seed1", ("signature_sha256", "loops_sha256", "svg_sha256")),
        ("A13", "roundtrip", ("signature_sha256", "svg_sha256")),
    ):
        r = _views(e, [other]).get(other)
        if r is not None and "default" in views:
            differ = [k.split("_")[0] for k in fields if d.get(k) is None or r.get(k) != d.get(k)]
            out[key] = _ok(not differ, f"{other} differs from default: {', '.join(differ)}" if differ else "")
    extra = {v: views[v] for v in ("interactive", "sequence") if v in views}
    if extra:
        bad = [v for v, r in extra.items() if r["exit"] != 0 or not all(r.get("written", {}).values())]
        out["A14"] = _ok(not bad, ", ".join(bad))
    if e.get("pair_check") is not None:
        out["A19"] = e["pair_check"]
    facts = e.get("facts") or {}
    if facts.get("na") and out["A17"]["pass"] is None:  # nucleic acid in the file, measured on no view
        drew = [v for v, r in _views(e, FIGURE_VIEWS).items() if r.get("na_chains_drawn")]
        out["A17"] = (
            _ok(False, f"nucleic acid drawn in the {', '.join(drew)} view(s), measured in none")
            if drew
            else _not_judged("no view draws the nucleic acid")
        )
    return out


def figure_criteria(views: dict[str, dict]) -> dict:
    """The criteria judged on the drawn figures alone (A3-A9, A11, A15-A18, A20), from the measured results of the
    views given ({view: {"result": measure_figure's dict, ...}}): {criterion: {"pass", "detail"}} for each one that
    the views can judge. criteria and `foldmap check` (check_figure) both judge by it."""
    e = {"views": views}
    out = {}
    d = _views(e, ["default"]).get("default", {})
    figures = _views(e, FIGURE_VIEWS)
    stuck = {v: r.get("stuck", []) for v, r in figures.items() if r.get("stuck")}
    over = {v: r.get("over_labels", []) for v, r in figures.items() if r.get("over_labels")}
    if figures:
        out["A3"] = _ok(not stuck, "; ".join(f"{v}: {', '.join(s[:4])}" for v, s in stuck.items()))
        out["A4"] = _ok(not over, "; ".join(f"{v}: {', '.join(s[:4])}" for v, s in over.items()))
    out["A5"] = _count(figures, "element_overlaps")
    out["A6"] = _count(figures, "text_collisions")
    out["A7"] = _count(figures, "off_page")
    out["A8"] = _count(_views(e, LABEL_VIEWS), "label_repeats")
    out["A9"] = _count(_views(e, ["default", "domains"]), "scope_silent")
    out["A11"] = _count(_views(e, ["default"]), "order")
    out["A15"] = _copies(_views(e, ["default", "focus"]))
    dom = _views(e, ["domains"])
    if dom:
        parts = [
            _count(dom, "domain_rules"),
            _count(dom, "domain_overlap", LIMITS["domain_overlap"]),
            _count(dom, "domain_name_collisions"),
        ]
        out["A16"] = _ok(all(p["pass"] for p in parts), "; ".join(p["detail"] for p in parts if p["detail"]))
        drew = [v for v, r in dom.items() if _details(r).get("domains") != {}]
        if out["A16"]["pass"] and not drew:  # domain_checks gives {} for a figure that drew no domain
            out["A16"] = _not_judged("the domains view drew no domain")
    na = next(((v, r) for v, r in (("default", d), *_views(e, ["focus"]).items()) if "na" in r), None)
    if na is not None:  # the default view, else the focus view when only it draws the nucleic acid (8QEM)
        out["A17"] = _nucleic(*na)
    mem = _views(e, ["membrane"]).get("membrane")
    if mem is not None:
        m = mem.get("membrane")
        if m is None:
            out["A18"] = _ok(False, "not measured")
        else:
            bad = ([] if m["band"] else ["no band drawn"]) + [f"TM helix {x}" for x in m["tm_short"]]
            bad += [f"loop in the slab {x}" for x in m["loops_in_slab"]]
            out["A18"] = _ok(not bad, "; ".join(bad))
    if "default" in views and d:
        # only a tether to another chain is bounded (ruling); a long one to the home chain is reported
        tethers = d.get("tethers", [])
        long = [(g, n, h, t) for g, n, h, t in tethers if n > LIMITS["tether_length"]]
        other = [f"{g} {n} ({h} to {t})" for g, n, h, t in long if t != h]
        same = [f"{g} {n} (home {h})" for g, n, h, t in long if t == h]
        more = lambda items: [f"(+{len(items) - 4} more)"] if len(items) > 4 else []  # noqa: E731
        shown = other[:4] + more(other)
        shown += [f"reported: {'; '.join(same[:4])}" + "".join(f" {m}" for m in more(same))] if same else []
        out["A20"] = _ok(not other, "; ".join(shown))
        if not other and not any(t != h for _, _, h, t in tethers):  # nothing A20 bounds
            why = "no tether to another chain" + (f"; {out['A20']['detail']}" if out["A20"]["detail"] else "")
            out["A20"] = _not_judged(why) if tethers else _ok(None, "no ligand tether")
    return out


def _details(r: dict) -> dict:
    return (r.get("checks") or {}).get("details") or {}


def _copies(results: dict[str, dict]) -> dict:
    """A15 over the views given: copies_alike counted (a count not measured fails), then by the ties' details
    (figcheck.figure_checks: deviation None for a tie too small to judge) does not apply when no view has copies to
    compare, not judged when none of its ties could be judged."""
    got = _count(results, "copies_alike", LIMITS["copies_alike"])
    if not got["pass"] or any("copies_alike" not in _details(r) for r in results.values()):
        return got
    ties = {f"{v} {tie}": t for v, r in results.items() for tie, t in sorted(_details(r)["copies_alike"].items())}
    if not ties:
        return _ok(None, "no copies to compare")
    if all(t.get("deviation") is None for t in ties.values()):
        return _not_judged("; ".join(f"{k} {t.get('judged')}" for k, t in ties.items()))
    return got


def _nucleic(view: str, r: dict) -> dict:
    """A17 on one view's measures: na_validate's failures, then each stem name within its limit of its own letters
    and nearer them than any other stem's (named by the view when it is not the default)."""
    bad = list(r["na"].get("failed", []))
    gaps = _details(r).get("stem_name_gaps")
    limit = r.get("stem_gap_limit")
    if gaps is None or (gaps and limit is None):
        bad.append("stem-name gaps not measured")
    else:
        for k, (own, other) in sorted(gaps.items()):
            if own > limit or own >= other:
                bad.append(f"stem name {k} {own:.3f} from its letters, {other:.3f} from another stem")
    detail = "; ".join(bad)
    if view != "default":
        detail = f"measured on the {view} view, as the default draws no nucleic acid" + (f": {detail}" if bad else "")
    return _ok(not bad, detail)


def _to_first_numbering(a: dict, b: dict) -> Callable[[int, int], tuple[int, int] | None] | None:
    """Takes a residue range of b's first chain to a's numbering by their sequences (figcheck.aligned_positions):
    the first and last of a's residues aligned to b's in the range. None when either view lacks its residues (a
    report from before C18), and the author numbers are read as they are."""
    from .figcheck import aligned_positions

    ra, rb = ((x["views"].get("beta", {}).get("result") or {}).get("chain_residues") for x in (a, b))
    if not ra or not rb:
        return None
    na, nb = ([n for n, _ in r["residues"]] for r in (ra, rb))
    pairs = aligned_positions([x for _, x in ra["residues"]], [x for _, x in rb["residues"]])

    def taken(first: int, last: int) -> tuple[int, int] | None:
        got = [na[i] for i, j in pairs if first <= nb[j] <= last]
        return (min(got), max(got)) if got else None

    return taken


def pair_check(a: dict, b: dict) -> dict:
    """A19: the beta views of the two entries of a pair, the first protein chain of each. Elements pair one-to-one
    by residue range as copies of one protein do (model.pair_places: the same kind, sharing at least half of the
    shorter, the largest overlaps first), b's ranges taken to a's numbering by sequence alignment (C18 ruling:
    2HKI numbers GAPDH ahead of 1RM5). The pairs must keep one order, and pair at least half of the larger
    chain's elements (and 3, when both have that many). One entry may hold an element the other lacks (listed, not
    failed); the generic numbers after it then shift, so every pair's names are compared modulo that shift
    (same_name; final review B, Important 1: comparing only up to the first unpaired element compared none in
    9YKJ/10UJ and 6YPX/6YQD). Not judged when no names could be compared."""
    from .model import pair_places

    ra = (a["views"].get("beta", {}).get("result") or {}).get("elements")
    rb = (b["views"].get("beta", {}).get("result") or {}).get("elements")
    if not ra or not rb:
        return _ok(False, "beta view not measured for both")

    def first_chain(els):
        chain = els[0]["chain"]
        return sorted((x for x in els if x["chain"] == chain), key=lambda x: (x["first"], x["last"]))

    ea, eb = first_chain(ra), first_chain(rb)
    taken = _to_first_numbering(a, b) or (lambda first, last: (first, last))
    places_b = [(x["kind"], *(taken(x["first"], x["last"]) or (-2, -3))) for x in eb]
    found = pair_places([(x["kind"], x["first"], x["last"]) for x in ea], places_b)
    pairs = sorted(found.items())
    bad, notes = [], [f"{x['id']} unpaired" for i, x in enumerate(ea) if i not in found]
    notes += [f"{y['id']} unpaired" for k, y in enumerate(eb) if k not in found.values()]
    # pairs sorted by the first entry's elements must take the second's in their order too: two chains' elements do
    # not overlap and both are in residue order, so this holds for residue-range pairs unless a pairing is wrong
    if [k for _, k in pairs] != sorted(k for _, k in pairs):
        bad.append("elements in another order")
    if len(pairs) < max(min(3, len(ea), len(eb)), max(len(ea), len(eb)) / 2):
        bad.append(f"only {len(pairs)} of {len(ea)} and {len(eb)} elements paired")
    paired_b, differ = set(found.values()), []
    for i, k in pairs:  # every pair's names, each number less the unpaired elements of its kind before it
        before_a = [x["label"] for j, x in enumerate(ea[:i]) if j not in found]
        before_b = [y["label"] for j, y in enumerate(eb[:k]) if j not in paired_b]
        if not same_name(ea[i]["label"], eb[k]["label"], before_a, before_b):
            differ.append(f"{ea[i]['label']}~{eb[k]['label']} ({ea[i]['id']}~{eb[k]['id']})")
    if not pairs and not bad:
        return _not_judged("no names to compare")
    bad += differ
    shown = bad + notes
    detail = f"names compared on all {len(pairs)} pairs" + (f", {len(differ)} differ" if differ else "")
    detail += "; " + "; ".join(shown[:8]) + (f" (+{len(shown) - 8} more)" if len(shown) > 8 else "") if shown else ""
    return _ok(not bad, detail)


_NUMBERED = re.compile(r"(\D*?)(\d+)(\D*)")  # a generic name: its kind (α, β, η ...), number, any mark after it


def same_name(a: str, b: str, before_a: list[str], before_b: list[str]) -> bool:
    """Are two paired elements' names the same, given the names of the unpaired elements before each in its own
    chain (A19)? Two numbered names of one kind and mark (α3, α3′) are the same when their numbers, each less the
    unpaired elements of that kind before it, agree (an extra helix in one entry shifts every later α by one); any
    other name (αN, a deposited letter) must be equal."""
    ma, mb = _NUMBERED.fullmatch(a), _NUMBERED.fullmatch(b)
    if not ma or not mb or (ma[1], ma[3]) != (mb[1], mb[3]):
        return a == b

    def shift(before: list[str], kind: str) -> int:
        return sum(1 for x in before if (m := _NUMBERED.fullmatch(x)) and m[1] == kind)

    return int(ma[2]) - shift(before_a, ma[1]) == int(mb[2]) - shift(before_b, mb[1])


# ------------------------------------------------------------------------------------------------ one figure

# The criteria one figure drawn in one process (`foldmap check`) cannot judge, and why: they are named not measured
# here, never passed. A10 and A17 only where they would apply (protein; nucleic acid).
NOT_HERE = {
    "A1": "judged on each view's own process in the campaign: exit code, files written, traceback, Python warnings",
    "A1c": "judged on the CA-only and P-only entries the campaign meets, which foldmap refuses",
    "A2": "the campaign times the default view's own process",
    "A2b": "the campaign draws it in the app's job process under the hosted limits",
    "A10": "needs mkdssp 3.1.2, which the campaign runs",
    "A12": "needs a second drawing under another PYTHONHASHSEED, the campaign's seed1 view",
    "A13": "needs a redrawing from the saved layout, the campaign's roundtrip view",
    "A14": "the campaign writes the interactive page and the sequence plot",
    "A17": "needs RNApolis base pairs, which the campaign's na_validate measures",
    "A19": "needs the other entry of a two-state pair",
}
INK_READ = {"A6": ("text_collisions",), "A7": ("off_page",), "A16": ("domain_name_collisions",)}  # figcheck.INK_CHECKS


def check_figure(drawn, spec, path) -> dict:
    """One drawn figure judged as the campaign judges the view of it its options make (measure_figure,
    figure_criteria): the membrane view when it asks for a membrane (or draws a band), else the domains view when it
    asks for domains (or draws them), else the default view. {"science": figcheck.science_checks' findings,
    "criteria": {criterion: {"pass", "detail"}, ...}, every criterion of CRITERIA}. A criterion it cannot judge has
    pass None and its reason under "here": NOT_HERE; one the campaign judges on other views, not on this one; a check
    the text font keeps from being measured (figcheck.ink_font_reason; for A16 only its domain names, its domain rules
    and panels are judged); a measure that reads "unmeasured" (unmeasured)."""
    lay = drawn.layout
    membrane = getattr(spec, "membrane", "off") not in (None, "off") or lay.membrane is not None
    domains = getattr(spec, "domains", None) not in (None, "", "none", "off") or bool(lay.domains)
    view = "membrane" if membrane else "domains" if domains else "default"  # as the options asked, as the campaign
    result = measure_figure(view, drawn, spec, path)
    judged = figure_criteria({view: {"result": result}})
    if view == "membrane" and domains:  # the campaign has no view with both: the domains are judged as its own
        judged["A16"] = figure_criteria({"domains": {"result": result}})["A16"]
    as_default = judged if view == "default" else figure_criteria({"default": {"result": result}})
    skipped = result["checks"].get("skipped") or {}
    here = {k: why for k, why in NOT_HERE.items() if k != "A10" or len(drawn.bb)}
    if lay.nucleic is None:
        here.pop("A17")
    out = {}
    for key in CRITERIA:
        got = judged.get(key, _ok(None))
        font = next((skipped[c] for c in INK_READ.get(key, ()) if c in skipped), None)
        if got["pass"] is not None and font:
            got = _domains_off_font(result, str(font)) if key == "A16" else {**_ok(None), "here": str(font)}
        elif got["pass"] is None and view != "default" and key in as_default:
            got = {
                **got,
                "here": f"the campaign judges it on other views, not on its {view} view: draw the figure "
                f"without --{view} to check it",
            }
        elif key in here:
            got = {**got, "here": here[key]}
        out[key] = unmeasured(got)
    return {"science": result["science"], "science_judged": result["science_judged"], "criteria": out}


def _domains_off_font(result: dict, font: str) -> dict:
    """A16 when the domain names cannot be measured (the text off figcheck.INK_FONT): its domain rules and panel
    overlaps measure no ink and are judged; a failure among them fails it, else it is not measured here."""
    dom = {"domains": result}
    parts = [_count(dom, "domain_rules"), _count(dom, "domain_overlap", LIMITS["domain_overlap"])]
    names = f"domain names not measured here ({font})"
    if all(p["pass"] for p in parts):
        return {**_ok(None), "here": f"{names}; domain rules and panel overlaps pass"}
    return _ok(False, "; ".join([p["detail"] for p in parts if p["detail"]] + [names]))


def unmeasured(got: dict) -> dict:
    """A result whose detail reads "unmeasured: <reason>" (the campaign's A10 without mkdssp: not judged; a pass in
    reports from before the final review B fix) as not measured: pass None, the reason under "here"."""
    detail = str(got.get("detail") or "")
    if got.get("pass") is not False and detail.startswith("unmeasured"):
        return {**_ok(None), "here": detail.removeprefix("unmeasured").lstrip(": ") or "unmeasured"}
    return got


def criteria_lines(found: dict[str, dict]) -> list[str]:
    """One line per criterion (pass, fail, not measured here and why, not judged and why, or does not apply), its
    detail indented."""
    lines = []
    for key, what in CRITERIA.items():
        got = unmeasured(found.get(key) or _ok(None))
        if got["pass"] is None:
            if got.get("here"):
                why = f"not measured here ({got['here']})"
            elif outcome(got) == NOT_JUDGED:
                why = NOT_JUDGED + (f" ({got['detail']})" if got.get("detail") else "")
            else:
                why = NOT_APPLICABLE
            lines.append(f"{key} {what}: {why}")
            continue
        lines.append(f"{key} {what}: {'pass' if got['pass'] else 'fail'}")
        if got["detail"]:
            lines.append(f"  {got['detail']}")
    return lines
