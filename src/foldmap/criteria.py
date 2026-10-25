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
  process cannot give are named not measured here, with the reason (NOT_HERE), never passed.
- like_checks, like_lines: H1-H5 of a structure drawn like a reference (plan 13a, `foldmap check B --like A`), with
  H2's and H3's baseline between the two drawn independently; their limits are frozen on the homologue tuning set
  (Task 13a-10)."""

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
    # H2, H3 (plan 13a, a figure drawn like a reference): the spec's starting limits (spec 4.5), frozen on the homologue
    # tuning set (Task 13a-10, docs/validation/homologues/tuning.json) before any held-out pair: 8 of its 11
    # positives pass them with room (H2 at most 0.28 / 0.36, H3 1.0); the 3 cross-family six-bladed propeller pairs
    # fail them (pairing differences, longer strands) and stay listed failures, not a reason to move them
    "like_place_median": 1.0,  # H2: median centre deviation of mates after the best rigid 2D motion, page units
    "like_place_p90": 2.5,  # H2: ... its 90th percentile
    "like_relations": 0.95,  # H3: share of mate pairs whose page dx, dy signs agree
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
    "A18": "membrane: TM helices in the band, loops and soluble elements out of the slab, none drawn in from beyond",
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
        if not m["targets"]:  # untethered (a rosette's axis cluster, render._axis_cluster): nothing to bound
            continue
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
    each TM helix (its CA run crosses the membrane's mid-plane and spans half its thickness, membrane.crosses) enters
    its band through a face and reaches past its middle (C27 ruling); loops run sideways inside a band
    (figcheck.SLAB_MARGIN from its faces) only on the first segment out of their own port inside it, up to clearing
    its element (membrane slab ruling); no soluble element is drawn in the slab (ruling R4), and no element whose CA
    run lies wholly beyond a face is drawn more than figcheck.SLAB_ENTRY into it (outside_in_slab, ruling T5-5). The
    membrane is the one
    the figure was drawn with (Layout.bilayer: found on every protein chain, whatever --chains draws); a figure with
    no band is judged as such, never by a membrane found again on the chains it draws."""
    from . import figcheck

    lay, bb = drawn.layout, drawn.bb
    mem = getattr(lay, "bilayer", None)
    if lay.membrane is None or mem is None:
        none = {"tm_helices": 0, "tm_short": [], "loops_in_slab": [], "soluble_in_slab": [], "outside_in_slab": []}
        return {"band": False, **none}
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
    name = next((t for t in ax.texts if re.fullmatch(r"na:(?:label:\d+|curve:\d+:label)", t.get_gid() or "")), None)
    if name is None:
        return None
    block = any(re.fullmatch(r"na:(?:stem|curve):\d+:block(:\d+)?", a.get_gid() or "") for a in ax.patches)
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


def _domain_found(*keys: str) -> Callable[[dict], list]:
    """A16's findings in one view's results: figure_checks keeps them under details["domains"][key]."""
    return lambda r: [x for k in keys for x in (_details(r).get("domains") or {}).get(k, [])]


def _count(results: dict[str, dict], key: str, limit: float = 0, found: Callable[[dict], list] | None = None) -> dict:
    """A check counted in each view's figcheck results: pass when every view's count is at most limit; a count
    that could not be measured (None) or is not a finite number fails. `found`: where its findings are, when not
    under details[key]."""
    bad = []
    for v, r in results.items():
        got = (r.get("checks") or {}).get(key, None) if "checks" in r else None
        if got is not None and not (isinstance(got, (int, float)) and np.isfinite(got)) and got != float("inf"):
            bad.append(f"{v} {key} is not a number ({got!r})")
        elif got is None:
            why = ((r.get("checks") or {}).get("skipped") or {}).get(key) or r.get("measure_error", "not measured")
            bad.append(f"{v} not measured ({str(why).strip().splitlines()[-1][:120]})")
        elif got > limit:
            details = (found(r) if found else (r["checks"].get("details") or {}).get(key)) or []
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
    if out["A11"]["pass"] and "elements" in d and not d["elements"]:  # no protein element: no order to judge
        out["A11"] = (
            _not_judged("no protein element drawn") if d.get("first_protein_chain") else _ok(None, "no protein")
        )
    labelled = _views(e, LABEL_VIEWS)
    if out["A8"]["pass"] and all(
        "elements" in r and not r["elements"] and not r.get("na_chains_drawn") for r in labelled.values()
    ):
        out["A8"] = _not_judged("no element or stem drawn")
    out["A15"] = _copies(_views(e, ["default", "focus"]))
    dom = _views(e, ["domains"])
    if dom:
        parts = [
            _count(dom, "domain_rules", found=_domain_found("cross_chain", "unalike", "small")),
            _count(dom, "domain_overlap", LIMITS["domain_overlap"], found=_domain_found("panel_overlaps")),
            _count(dom, "domain_name_collisions", found=_domain_found("name_collisions")),
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
            bad += [f"soluble element in the slab {x}" for x in m.get("soluble_in_slab", [])]
            bad += [f"element beyond a face drawn into the slab {x} Å" for x in m.get("outside_in_slab", [])]
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
    compare, not judged when none of its ties could be judged or when copies are drawn with no tie at all
    (figcheck.copies_untied: no symmetry drew them alike; random-100 held-out 5K7Z said "no copies to compare")."""
    got = _count(results, "copies_alike", LIMITS["copies_alike"])
    if not got["pass"] or any("copies_alike" not in _details(r) for r in results.values()):
        return got
    ties = {f"{v} {tie}": t for v, r in results.items() for tie, t in sorted(_details(r)["copies_alike"].items())}
    if not ties:
        untied = {v: _details(r).get("copies_untied") or [] for v, r in results.items()}
        if any(untied.values()):  # copies drawn with no symmetry tie: nothing promised alike, but copies there are
            return _not_judged(
                "copies drawn with no symmetry tie: " + "; ".join(f"{v} {', '.join(c)}" for v, c in untied.items() if c)
            )
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


PAIR_IDENTITY = 0.9  # A19: two entries' first chains are one protein at this identity over their aligned residues
PAIR_COVERAGE = 0.5  # ... with this fraction of the shorter chain's residues identical


def one_protein(a: dict, b: dict) -> dict | None:
    """Are the two entries' first chains (their beta views' chain_residues) one protein, by sequence: identity over
    the residues the global alignment pairs (figcheck.aligned_positions) at least PAIR_IDENTITY, and that many
    identical residues at least PAIR_COVERAGE of the shorter chain? A UniProt accession is not enough: a polyprotein
    (SARS-CoV-2 pp1ab, P0DTD1) names nsp1 and nsp5 alike. None when either view lacks its residues."""
    ra, rb = ((x["views"].get("beta", {}).get("result") or {}).get("chain_residues") for x in (a, b))
    if not ra or not rb:
        return None
    return one_protein_names([n for _, n in ra["residues"]], [n for _, n in rb["residues"]])


def one_protein_names(xs, ys) -> dict:
    """A19's rule (one_protein) on two chains' residue names, in order: {one, identity, aligned, lengths}. foldmap
    diff maps two states of one protein by sequence with it (plan 13e)."""
    from .figcheck import aligned_positions

    xs, ys = list(xs), list(ys)
    pairs = aligned_positions(xs, ys)
    same = sum(xs[i] == ys[j] for i, j in pairs)
    identity = same / len(pairs) if pairs else 0.0
    one = identity >= PAIR_IDENTITY and same >= PAIR_COVERAGE * min(len(xs), len(ys))
    return {"one": one, "identity": round(identity, 3), "aligned": len(pairs), "lengths": [len(xs), len(ys)]}


def pair_check(a: dict, b: dict) -> dict:
    """A19: the beta views of the two entries of a pair, the first protein chain of each. Elements pair one-to-one
    by residue range as copies of one protein do (model.pair_places: the same kind, sharing at least half of the
    shorter, the largest overlaps first), b's ranges taken to a's numbering by sequence alignment (C18 ruling:
    2HKI numbers GAPDH ahead of 1RM5). The pairs must keep one order, and pair at least half of the larger
    chain's elements (and 3, when both have that many). One entry may hold an element the other lacks (listed, not
    failed); the generic numbers after it then shift, so every pair's names are compared modulo that shift
    (same_name; final review B, Important 1: comparing only up to the first unpaired element compared none in
    9YKJ/10UJ and 6YPX/6YQD). Not judged when no names could be compared; does not apply when the two chains are
    not one protein by sequence (one_protein)."""
    from .model import pair_places

    ra = (a["views"].get("beta", {}).get("result") or {}).get("elements")
    rb = (b["views"].get("beta", {}).get("result") or {}).get("elements")
    if not ra or not rb:
        return _ok(False, "beta view not measured for both")

    def first_chain(els):
        chain = els[0]["chain"]
        return sorted((x for x in els if x["chain"] == chain), key=lambda x: (x["first"], x["last"]))

    seq = one_protein(a, b)
    if seq is not None and not seq["one"]:  # held-out P3: nsp5 against nsp1 of one polyprotein accession
        return _ok(
            None,
            f"not one protein: {seq['identity']:.0%} identity over {seq['aligned']} aligned residues "
            f"({seq['lengths'][0]} and {seq['lengths'][1]} residues)",
        )
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
    other name (αN, a deposited letter) must be equal. Equal names always match (2BBJ/2BBH: an η1 drawn in one entry
    only, the other's not among its drawn elements, must not make η2~η2 differ)."""
    if a == b:  # the same name is what A19 asks for; the shift only forgives different numbers (F1 review, I1)
        return True
    ma, mb = _NUMBERED.fullmatch(a), _NUMBERED.fullmatch(b)
    if not ma or not mb or (ma[1], ma[3]) != (mb[1], mb[3]):
        return a == b

    def shift(before: list[str], kind: str) -> int:
        return sum(1 for x in before if (m := _NUMBERED.fullmatch(x)) and m[1] == kind)

    return int(ma[2]) - shift(before_a, ma[1]) == int(mb[2]) - shift(before_b, mb[1])


# ---------------------------------------------------------------------------- drawn like a reference: H1-H5

LIKE_CHECKS = {  # plan 13a, spec 4.5: a structure B drawn like a reference A (--like); `foldmap check B --like A`
    "H1": "drawn alike exactly when one fold",
    "H2": "mates in the same places",
    "H3": "the map reads the same way",
    "H4": "mates carry one name",
    "H5": "still a good figure",
}
LIKE_MIN_MATES = 3  # H2, H3: judged on this many mates drawn at least (as figcheck.COPY_MIN_MATCHED: 2 always fit)
LIKE_APART = 1.0  # H3: two mates more than this many page units apart, in both figures, are compared ...
LIKE_AXIS = 0.5  # ... on each page axis that parts them by more than this in both figures


def like_pairs(lay_a, lay_b, mates) -> list[tuple[str, str, tuple[str, ...]]]:
    """(A element id, B element id, B's parts) of every mate (Mates.mate: one to one) both layouts drew, not as a
    ghost, in A's order. B's parts: the mate and the pieces of it (Mates.pieces: a helix B splits where A does not)
    drawn, in sequence: like.like_from lays them end to end along the mate and names the first A's name."""
    pa, pb = lay_a.placed, lay_b.placed

    def drawn(placed, k):
        return k in placed and not placed[k].ghost

    parts: dict[str, list[str]] = {}
    for b, a in [*mates.mate.items(), *mates.pieces.items()]:
        if drawn(pb, b):
            parts.setdefault(a, []).append(b)
    found = [
        (a, b, tuple(sorted(parts[a], key=lambda k: (pb[k].sse.start, k))))
        for b, a in mates.mate.items()
        if drawn(pa, a) and drawn(pb, b)
    ]
    return sorted(found, key=lambda x: (pa[x[0]].sse.start, x[0]))


def element_centres(lay, ids) -> np.ndarray:
    return np.array([(lay.placed[k].cx, lay.placed[k].cy) for k in ids], float).reshape(-1, 2)


def part_centres(lay, pairs) -> np.ndarray:
    """Each mate's centre on B's page (like_pairs' pairs): its own, or for a helix drawn in parts the middle of their
    run, from its first end to its last along the mate's axis (like.like_from centres the run on A's helix)."""
    out = []
    for _, b, parts in pairs:
        if len(parts) < 2:
            out.append((lay.placed[b].cx, lay.placed[b].cy))
            continue
        ends = np.array([e for k in parts for e in (lay.placed[k].n_port, lay.placed[k].c_port)], float)
        along = ends @ lay.placed[b].direction
        out.append(tuple((ends[int(np.argmin(along))] + ends[int(np.argmax(along))]) / 2))
    return np.array(out, float).reshape(-1, 2)


def _name(lay, k: str) -> str:
    return lay.placed[k].label or k


def _independently(baseline: dict | None) -> str:
    """The baseline (A and B drawn independently) as H2's or H3's detail ends with it; '' when there is none."""
    if baseline is None:
        return ""
    if "median" in baseline:
        return f"; drawn independently: median {baseline['median']:.3f}, 90th percentile {baseline['p90']:.3f}"
    return f"; drawn independently: {baseline['share']:.4f} of {baseline['pairs']} pairs"


def _too_few(pairs) -> dict | None:
    if len(pairs) < LIKE_MIN_MATES:
        return _not_judged(f"{len(pairs)} mate{'' if len(pairs) == 1 else 's'} drawn, fewer than {LIKE_MIN_MATES}")
    return None


def place_deviation(lay_a, lay_b, pairs) -> dict | None:
    """H2's measure: the median and 90th percentile (numpy's linear interpolation) of the mates' centre deviations
    after the best rigid 2D motion of B's centres (part_centres) onto A's (figcheck.fit_deviations, A15's
    Procrustes: a turn and a shift, no mirror), rounded to 1e-3 page units, and the farthest mate; None under
    LIKE_MIN_MATES mates."""
    from .figcheck import fit_deviations

    if len(pairs) < LIKE_MIN_MATES:
        return None
    dev = fit_deviations(element_centres(lay_a, [a for a, *_ in pairs]), part_centres(lay_b, pairs))
    far = int(np.argmax(dev))
    a, b, _ = pairs[far]
    return {
        "median": round(float(np.median(dev)), 3),
        "p90": round(float(np.percentile(dev, 90)), 3),
        "mates": len(pairs),
        "farthest": f"{_name(lay_b, b)}~{_name(lay_a, a)} {float(dev[far]):.3f}",
    }


def h2_places(lay_a, lay_b, pairs, baseline: dict | None = None) -> dict:
    """H2: matched elements sit in the same places (place_deviation within LIMITS' like_place_median and
    like_place_p90). `baseline`: the same measure between A and B drawn independently, reported."""
    got = place_deviation(lay_a, lay_b, pairs)
    if got is None:
        out = _too_few(pairs)
    else:
        median, p90 = LIMITS["like_place_median"], LIMITS["like_place_p90"]
        out = _ok(
            got["median"] <= median and got["p90"] <= p90,
            f"median {got['median']:.3f}, 90th percentile {got['p90']:.3f} page units over {got['mates']} mates "
            f"(limits {median}, {p90}); farthest {got['farthest']}",
        )
    return {**out, "detail": out["detail"] + _independently(baseline), "value": got, "baseline": baseline}


def relations(lay_a, lay_b, pairs) -> dict | None:
    """H3's measure: of every two mates more than LIKE_APART page units apart in both figures, the share whose page
    dx and dy have one sign in both, on each axis that parts them by more than LIKE_AXIS in both (none: the relation
    turned, a disagreement), rounded to 1e-4; and the mated sheets whose mated strands stand in another order
    (sheet_orders). None under LIKE_MIN_MATES mates or with no two mates far enough apart."""
    if len(pairs) < LIKE_MIN_MATES:
        return None
    ca, cb = element_centres(lay_a, [a for a, *_ in pairs]), part_centres(lay_b, pairs)
    judged, turned = 0, []
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            da, db = ca[j] - ca[i], cb[j] - cb[i]
            if np.hypot(*da) <= LIKE_APART or np.hypot(*db) <= LIKE_APART:
                continue
            judged += 1
            axes = [k for k in (0, 1) if abs(da[k]) > LIKE_AXIS and abs(db[k]) > LIKE_AXIS]
            if not axes or any(da[k] * db[k] <= 0 for k in axes):
                turned.append(f"{_name(lay_b, pairs[i][1])}-{_name(lay_b, pairs[j][1])}")
    if not judged:
        return None
    return {
        "share": round(1 - len(turned) / judged, 4),
        "pairs": judged,
        "turned": turned,
        "sheets": sheet_orders(lay_a, lay_b, pairs),
    }


def sheet_orders(lay_a, lay_b, pairs) -> list[str]:
    """The mated sheets whose mated strands stand in another order: for each sheet block of B (Layout.sheets) and the
    block of A its strands' mates are in (two or more strands), any two strands whose mates stand more than LIKE_AXIS
    apart across A's sheet (the page direction square to its strands' mean axis) and do not stand so in B, along that
    same direction. Each as "β2, β1, β3 for β1, β2, β3" (B's names in B's order, then A's in A's)."""
    mate = {b: a for a, b, _ in pairs}
    block_a = {k: n for n, ids in enumerate(lay_a.sheets) for k in ids}
    out = []
    for ids in lay_b.sheets:
        groups: dict[int, list[str]] = {}
        for b in ids:
            if b in mate and mate[b] in block_a and b in lay_b.placed and lay_b.placed[b].sse.kind == "E":
                groups.setdefault(block_a[mate[b]], []).append(b)
        for n in sorted(groups):
            strands = groups[n]
            if len(strands) < 2:
                continue
            angles = np.array([lay_a.placed[mate[b]].angle for b in strands])
            mean = 0.5 * np.arctan2(np.sin(2 * angles).sum(), np.cos(2 * angles).sum())  # an axis: modulo pi
            across = np.array([-np.sin(mean), np.cos(mean)])
            if across[0] < -1e-9 or (abs(across[0]) <= 1e-9 and across[1] < 0):  # read left to right (or upward)
                across = -across
            pa = element_centres(lay_a, [mate[b] for b in strands]) @ across
            pb = element_centres(lay_b, strands) @ across
            if all(
                (pb[j] - pb[i]) * (pa[j] - pa[i]) > 0
                for i in range(len(strands))
                for j in range(i + 1, len(strands))
                if abs(pa[j] - pa[i]) > LIKE_AXIS
            ):
                continue
            order_b = [_name(lay_b, strands[k]) for k in np.argsort(pb, kind="stable")]
            order_a = [_name(lay_a, mate[strands[k]]) for k in np.argsort(pa, kind="stable")]
            out.append(f"{', '.join(order_b)} for {', '.join(order_a)}")
    return out


def h3_relations(lay_a, lay_b, pairs, baseline: dict | None = None) -> dict:
    """H3: the map reads the same way (relations' share at least LIMITS' like_relations) and each mated
    sheet keeps one order (sheet_orders; a sheet whose pairing differs fails too: the excuse for one waits for its
    redesign, ROADMAP). `baseline`: the share between A and B drawn independently, reported."""
    got = relations(lay_a, lay_b, pairs)
    if got is None:
        out = _too_few(pairs) or _not_judged(f"no two mates more than {LIKE_APART} page unit apart in both figures")
    else:
        limit = LIMITS["like_relations"]
        bad = [f"{len(got['turned'])} turned: {', '.join(got['turned'][:6])}"] if got["share"] < limit else []
        sheets = [f"sheet order {s}" for s in got["sheets"]]
        share = f"{got['share']:.4f} of {got['pairs']} pairs agree (limit {limit})"
        out = _ok(not bad and not got["sheets"], "; ".join([share, *bad, *sheets]))
    return {**out, "detail": out["detail"] + _independently(baseline), "value": got, "baseline": baseline}


def h4_names(lay_a, lay_b, mates, note: str = "") -> dict:
    """H4: mates carry one name, by A19's rule (same_name: each number less the elements of its kind before it that
    have no mate, Mates.only_a in A and the extras in B). A helix B draws in parts carries its mate's name on its
    first part in sequence, the next ones primed (like.like_from): that first part is compared. A mate named by its
    own helix type (Like.crossed: an α-helix whose mate is a 3-10 helix) passes when the note names it so ("α1 for
    η1", like.like_note)."""
    pairs = like_pairs(lay_a, lay_b, mates)
    if not pairs:
        return _not_judged("no mates drawn")
    chain_a, chain_b = lay_a.placed[pairs[0][0]].sse.chain, lay_b.placed[pairs[0][1]].sse.chain

    def order(lay, chain):
        els = [(p.sse.start, k) for k, p in lay.placed.items() if not p.ghost and p.sse.chain == chain]
        return [k for _, k in sorted(els)]

    order_a, order_b = order(lay_a, chain_a), order(lay_b, chain_b)
    mated_a, mated_b = set(mates.mate.values()), set(mates.mate) | set(mates.pieces)
    differ, crossed = [], []
    for a, _, (b, *_) in pairs:  # the first part in sequence carries the name
        before_a = [_name(lay_a, k) for k in order_a[: order_a.index(a)] if k not in mated_a]
        before_b = [_name(lay_b, k) for k in order_b[: order_b.index(b)] if k not in mated_b]
        na, nb = lay_a.placed[a].label, lay_b.placed[b].label
        if same_name(na, nb, before_a, before_b):
            continue
        (crossed if f"{nb} for {na}" in note else differ).append(f"{nb}~{na}")
    detail = f"names compared on all {len(pairs)} mates"
    if differ:
        detail += f", {len(differ)} differ: {', '.join(differ[:8])}"
    if crossed:
        detail += f"; named by their own helix type, as the note says: {', '.join(crossed)}"
    return {**_ok(not differ, detail), "value": {"mates": len(pairs), "differ": differ}}


def h1_matched(aln, mates, alike: bool, refused: bool, kind: str | None = None, note: str = "") -> dict:
    """H1: the pair is drawn alike exactly when it is one fold (Mates.same_fold: align.SAME_FOLD_TM and
    SAME_FOLD_FRACTION). `kind` (a tuning or held-out pair's): "positive" passes when one fold and drawn alike;
    "negative" when not one fold and refused with the note (like.refused_note), so a forced negative (drawn alike)
    fails; None (`foldmap check`, no kind known): the drawing agrees with same_fold either way. `refused`: drawn on
    its own because not one fold, with the note saying so (cli._drawn_like, make_layout's like_out). Not judged when
    the chains were not aligned (a frame that cannot be carried: like.uncarried_note, `note`)."""
    from .align import SAME_FOLD_FRACTION, SAME_FOLD_TM

    if aln is None or mates is None:
        return _not_judged(f"not aligned: {note}" if note else "not aligned")
    one = mates.same_fold
    refused = refused and not alike
    facts = (
        f"TM-score {aln.tm_score:.2f} (one fold at {SAME_FOLD_TM}), {len(mates.mate)} mates of the "
        f"{mates.smaller} elements of the chain with fewer (at {SAME_FOLD_FRACTION})"
    )
    said = "drawn alike" if alike else "refused with the note" if refused else "drawn on its own, not refused"
    if kind == "negative":
        ok = not one and refused
        what = "a negative pair " + ("drawn alike: forced" if alike else said)
    elif kind == "positive":
        ok = one and alike
        what = "a positive pair " + said
    else:
        ok = (one and alike) or (not one and refused)
        what = ("one fold, " if one else "not one fold, ") + said
    value = {"tm_score": round(float(aln.tm_score), 4), "same_fold": one, "alike": alike, "refused": refused}
    return {**_ok(ok, f"{what}: {facts}"), "value": value}


def h5_good_figure(judged_like: dict, judged_alone: dict, who: str = "B") -> dict:
    """H5: likeness is not bought with a broken page: B drawn like A passes every S- and A-check that B drawn alone
    passes (each a check_figure result: {"science", "science_judged", "criteria"}; a check passes when it found
    nothing and had something to judge, a criterion when its pass is True). Not judged when B alone passes none.
    `who`: the structure the detail names as drawn alone (B; a compare or diff page's panel names its own,
    compare.panel_checks)."""

    def passed(found: dict) -> list[str]:
        judged = found.get("science_judged") or {}
        science = [k for k, v in found["science"].items() if v is not None and not v and judged.get(k) != 0]
        return science + [k for k, v in found["criteria"].items() if v.get("pass") is True]

    alone, like = passed(judged_alone), set(passed(judged_like))
    if not alone:
        return _not_judged(f"{who} drawn alone passes no S- or A-check")
    lost = [k for k in alone if k not in like]

    def why(k: str) -> str:
        got = judged_like["criteria"].get(k)
        if got is None:
            return f"{k} ({len(judged_like['science'].get(k) or [])} problem(s))"
        return f"{k} ({outcome(unmeasured(got))}{': ' + got['detail'][:120] if got.get('detail') else ''})"

    detail = f"{len(alone) - len(lost)} of the {len(alone)} checks {who} alone passes"
    detail += f"; lost: {'; '.join(why(k) for k in lost)}" if lost else ""
    return {**_ok(not lost, detail), "value": {"alone": alone, "lost": lost}}


def like_checks(
    lay_a,
    lay_b,
    aln,
    mates,
    note: str,
    alike: bool,
    lay_alone=None,
    judged_like: dict | None = None,
    judged_alone: dict | None = None,
    kind: str | None = None,
    refused: bool = False,
) -> dict[str, dict]:
    """H1-H5 of B drawn like A (spec 4.5): lay_a the reference's layout (None when it was not laid out), lay_b B's
    drawn like it, aln and mates the alignment and element mates (None where not reached), note B's like note
    (Layout.like_note), alike whether B took A's places (a Like was given to build_layout). `lay_alone`: B drawn on
    its own (the same chain), for H2's and H3's baseline (drawn independently, the manuscript's sentence);
    judged_like and judged_alone: check_figure's results of the two, for H5; `kind` and `refused` (B refused as not
    the reference's fold: make_layout's like_out): h1_matched's. A
    check with nothing to judge is not judged, never a pass; H2-H4 are not judged where B was not drawn alike."""
    out = {"H1": h1_matched(aln, mates, alike, refused, kind, note)}
    measured = lay_a is not None and mates is not None
    pairs = like_pairs(lay_a, lay_b, mates) if measured else []
    alone = like_pairs(lay_a, lay_alone, mates) if measured and lay_alone is not None else []
    place, rel = place_deviation(lay_a, lay_alone, alone), relations(lay_a, lay_alone, alone)
    if measured and alike:
        out["H2"] = h2_places(lay_a, lay_b, pairs, place)
        out["H3"] = h3_relations(lay_a, lay_b, pairs, rel)
        out["H4"] = h4_names(lay_a, lay_b, mates, note)
    else:
        why = note or "not drawn alike"  # the note says why: not the same fold, a frame that cannot be carried
        out["H2"] = {**_not_judged(why + _independently(place)), "value": None, "baseline": place}
        out["H3"] = {**_not_judged(why + _independently(rel)), "value": None, "baseline": rel}
        out["H4"] = _not_judged(why)
    out["H5"] = (
        h5_good_figure(judged_like, judged_alone)
        if judged_like is not None and judged_alone is not None
        else _not_judged("B was not drawn alone")
    )
    return out


def like_lines(found: dict[str, dict]) -> list[str]:
    """One line per check H1-H5 (pass, fail, not judged and why), its detail indented, as criteria_lines prints
    A1-A20."""
    lines = []
    for key, what in LIKE_CHECKS.items():
        got = found.get(key) or _not_judged("")
        if got["pass"] is None:
            lines.append(f"{key} {what}: {NOT_JUDGED}" + (f" ({got['detail']})" if got.get("detail") else ""))
            continue
        lines.append(f"{key} {what}: {'pass' if got['pass'] else 'fail'}")
        if got["detail"]:
            lines.append(f"  {got['detail']}")
    return lines


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
