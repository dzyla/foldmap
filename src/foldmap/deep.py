# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Deep checks (`foldmap check --deep`, plan 28 Task 7, ruling Q34) and the machine-readable check report
(`foldmap check --json`).

The deep checks are what `foldmap check` adds when asked, on top of the science checks S1-S6 and the page criteria
A1-A20, each with a status of pass, fail, unavailable (with the reason) or not applicable, and the kind of evidence it
is (science.KIND's layers):

- D1 (consistency): the figure drawn a second time in the same process is the same figure: its layout signature
  (signature) and its SVG bytes hash the same. Not A12 (that needs another PYTHONHASHSEED, a second process).
- D2 (external): foldmap's DSSP agrees with mkdssp 3.1.2 residue by residue (criteria.dssp_agreement, A10's measure)
  at criteria.LIMITS["dssp_agree"]; run only when the MKDSSP environment variable names the program.
- D3 (external): foldmap's canonical and wobble base pairs on the model as deposited against RNApolis' (Saenger
  classes XIX, XX and XXVIII), F1 at nacheck.THRESHOLDS["pair_f1"]; run only when NA_RNAPOLIS_PYTHON names the
  Python of an environment that has rnapolis.

An unavailable check is never counted as passed: the report's verdict is "fail" when any check fails, else
"incomplete" when any check is unavailable, else "pass" (so `foldmap check --json` without --deep, whose A1, A2, A10,
A12 ... are not measured in one process, always reads "incomplete" at best). The exit code is unchanged by an
unavailable check: 1 on a failure, else 0.

The report (report(); SCHEMA, VERSION): {"schema": "foldmap-check", "version": 1, "structure", "foldmap": the
version, "deep": whether D1-D3 were run, "verdict", "counts": {status: n}, "kinds": {kind: what it is evidence of},
"checks": [{"code", "group" (science, criteria, floor, like, deep), "title", "kind", "status", "detail", "findings",
"reason", "measures"}, ...]}. A field is added in a later version only with VERSION raised. The floor group (plan 28
Task 4) holds R1, the readability floor (floor.check), after the criteria: a check id added, the shape unchanged;
its measures hold no wall time (the time a figure took is said in the text only, ruling Q41), so the report is the
same on every run. S5's entry holds "measures" (plan 28 Task 5, rulings Q42, Q43) when the figure's domains are
inferred and split: {"domain_splits": each split with its confidence, "confidence": what that is, not a
probability} (domain_confidence); its status is S5's as before, and its reason, when unmeasured, says what the
measures are (domain_confidence.S5_REASON). No top-level key is added."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

from .science import CONSISTENCY, EXTERNAL, KINDS, NOTHING, READABILITY, SCIENCE, check_kinds

SCHEMA = "foldmap-check"
VERSION = 1
PASS, FAIL, UNAVAILABLE, NOT_JUDGED, NOT_APPLICABLE = "pass", "fail", "unavailable", "not_judged", "not_applicable"
STATUSES = (PASS, FAIL, UNAVAILABLE, NOT_JUDGED, NOT_APPLICABLE)
DEEP = {
    "D1": "same figure when drawn again (layout signature and SVG)",
    "D2": "secondary structure agrees with mkdssp 3.1.2 (A10's measure)",
    "D3": "base pairs agree with RNApolis",
}
DEEP_KIND = {"D1": CONSISTENCY, "D2": EXTERNAL, "D3": EXTERNAL}
# H1-H5 (a figure drawn like a reference, criteria.LIKE_CHECKS) are none of science.KIND's layers: named apart
LIKE = "like"
FLOOR = "floor"  # the group of R1, the readability floor (floor.py, plan 28 Task 4): kind readability, as A1-A20
LIKE_KIND = "drawn like a reference (H1-H5): the figure against the reference's figure and its own drawn alone"
MKDSSP_VERSION = "mkdssp version 3.1.2"
TIMEOUT = 600  # s for an external program
_TEXT = {"capture_output": True, "encoding": "utf-8", "errors": "replace"}  # a tool's output never stops the check

# RNApolis' canonical and wobble pairs of a file, run in the Python that has it (scripts/na_validate.py's reference)
_RNAPOLIS = """
import gzip, json, sys, tempfile
from pathlib import Path
from rnapolis.annotator import find_pairs
from rnapolis.parser import read_3d_structure
from rnapolis.util import handle_input_file
src = Path(sys.argv[1])
with tempfile.TemporaryDirectory() as tmp:
    name = src.name[:-3] if src.suffix == ".gz" else src.name
    path = Path(tmp) / ("model.pdb" if name.lower().endswith((".pdb", ".ent")) else "model.cif")
    path.write_bytes(gzip.decompress(src.read_bytes()) if src.suffix == ".gz" else src.read_bytes())
    found = find_pairs(read_3d_structure(handle_input_file(str(path)), 1), 1)[0]
out = []
for bp in found:
    if bp.saenger is not None and bp.saenger.value in ("XIX", "XX", "XXVIII"):
        a, b = bp.nt1.auth, bp.nt2.auth
        out.append(sorted([[a.chain, a.number, (a.icode or "").strip()], [b.chain, b.number, (b.icode or "").strip()]]))
print(json.dumps(out))
"""


def result(status: str, detail: str = "", **measures) -> dict:
    """One deep check's result: its status, a line saying what was found or why it was not measured, its
    measures."""
    assert status in STATUSES, status
    return {"status": status, "detail": detail, **measures}


# --------------------------------------------------------------------------------------------- D1 repeat run


def signature(lay, loops) -> str:
    """The drawn layout as JSON text: every element's and ghost's box, the stems, the bounds, the sheets, the
    domains, the nucleic-acid paths and each loop's points. Text, so -0.0 and 0.0 differ."""

    def box(p):
        return [p.cx, p.cy, p.angle, p.length, p.width]

    sig = {
        "placed": {k: box(p) for k, p in sorted(lay.placed.items())},
        "ghosts": [[g.label, *box(g)] for g in lay.ghosts],
        "stems": [[d.id, *box(d)] for d in lay.stems],
        "bounds": list(lay.bounds),
        "sheets": [list(x) for x in lay.sheets],
        "domains": [[n, list(ids)] for n, ids in lay.domains],
        "na": [[p.kind, p.strand, [list(q) for q in p.points]] for p in lay.na_paths],
        "loops": [[lp.a_id, lp.b_id, [list(map(float, q)) for q in lp.points]] for lp in loops],
    }
    return json.dumps(sig, sort_keys=True, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _svg(figure) -> bytes:
    from .render import save

    with tempfile.TemporaryDirectory() as tmp:
        return save(figure, Path(tmp) / "figure.svg").read_bytes()


def repeat_run(spec, drawn, draw=None) -> dict:
    """D1: the figure `drawn` (of `spec`) drawn again by `draw` (spec.draw_spec): its signature and SVG hash the
    same."""
    import matplotlib.pyplot as plt

    from .spec import draw_spec

    again = (draw or draw_spec)(spec)
    try:
        first = (_sha256(signature(drawn.layout, drawn.loops).encode()), _sha256(_svg(drawn.figure)))
        second = (_sha256(signature(again.layout, again.loops).encode()), _sha256(_svg(again.figure)))
    finally:
        plt.close(again.figure)
    hashes = {"signature_sha256": [first[0], second[0]], "svg_sha256": [first[1], second[1]]}
    differ = [k for k, (a, b) in hashes.items() if a != b]
    if differ:
        return result(FAIL, f"drawn again, its {' and its '.join(d.split('_')[0] for d in differ)} differ", **hashes)
    return result(PASS, "drawn twice: the same layout signature and SVG", **hashes)


# --------------------------------------------------------------------------------------------- D2 mkdssp


def parse_dssp(text: str) -> list[tuple[str, int, str, str]]:
    """(auth chain, number, insertion code, code) per residue of an mkdssp output, in file order (a bend S and a
    blank are '-', as foldmap does not assign bends)."""
    out, started = [], False
    for line in text.splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 17 or line[13] == "!":
            continue
        tail = line[137:].split()
        chain = tail[1] if len(tail) > 1 else line[11]
        code = line[16] if line[16] not in " S" else "-"
        out.append((chain, int(line[5:10]), line[10].strip(), code))
    return out


def mkdssp_rows(path, exe: str | None = None) -> list[tuple[str, int, str, str]] | str:
    """mkdssp 3.1.2's codes for the file (parse_dssp), or why there are none. `exe`: the program (default: the
    MKDSSP environment variable; none set: unavailable, the program is never looked for elsewhere)."""
    exe = exe if exe is not None else os.environ.get("MKDSSP", "")
    if not exe:
        return "MKDSSP is not set (the path of mkdssp 3.1.2)"
    try:
        done = subprocess.run([exe, "--version"], timeout=60, **_TEXT)
    except (OSError, subprocess.SubprocessError) as err:
        return f"cannot run mkdssp ({exe}): {err}"
    said = f"{done.stdout}\n{done.stderr}"  # some builds print their version to stderr
    if MKDSSP_VERSION not in said:
        return f"mkdssp 3.1.2 is required (found: {said.strip()[:60] or 'no version output'})"
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(path)
        if src.suffix == ".gz":
            plain = Path(tmp) / src.name[:-3]
            plain.write_bytes(gzip.decompress(src.read_bytes()))
            src = plain
        out = Path(tmp) / "out.dssp"
        try:
            done = subprocess.run([exe, str(src), str(out)], timeout=TIMEOUT, **_TEXT)
        except subprocess.TimeoutExpired:
            return f"mkdssp took longer than {TIMEOUT} s"
        if done.returncode or not out.is_file():
            tail = (done.stderr.strip().splitlines() or ["failed"])[-1]
            return f"mkdssp cannot read it: {tail[:160]}"
        return parse_dssp(out.read_text(encoding="utf-8", errors="replace"))


def dssp_check(path, rows: list | str | None = None, protein: bool = True) -> dict:
    """D2: criteria.dssp_agreement against mkdssp's codes (`rows`; None: run mkdssp, mkdssp_rows) at
    LIMITS["dssp_agree"]."""
    from .criteria import LIMITS, dssp_agreement

    if not protein:
        return result(NOT_APPLICABLE, "no protein residues")
    rows = mkdssp_rows(path) if rows is None else rows
    got = dssp_agreement(Path(path), rows)
    if got["agree"] is None:
        return result(UNAVAILABLE, got["unmeasured"])
    ok = got["agree"] >= LIMITS["dssp_agree"]
    detail = f"{got['agree']:.4f} of {got['residues']} residues agree (limit {LIMITS['dssp_agree']})"
    return result(PASS if ok else FAIL, detail, **{k: v for k, v in got.items() if k != "agree"}, agree=got["agree"])


# --------------------------------------------------------------------------------------------- D3 RNApolis


def rnapolis_pairs(path, python: str | None = None) -> set[frozenset] | str:
    """RNApolis' canonical and wobble pairs of the file, as frozensets of (chain, number, icode), or why there are
    none. `python`: the Python that has rnapolis (default: NA_RNAPOLIS_PYTHON; none set: unavailable)."""
    python = python if python is not None else os.environ.get("NA_RNAPOLIS_PYTHON", "")
    if not python:
        return "NA_RNAPOLIS_PYTHON is not set (the Python of an environment with rnapolis)"
    try:
        done = subprocess.run([python, "-I", "-c", _RNAPOLIS, str(path)], timeout=TIMEOUT, **_TEXT)
    except (OSError, subprocess.SubprocessError) as err:
        return f"cannot run RNApolis ({python}): {err}"
    if done.returncode:
        if "No module named 'rnapolis" in done.stderr:
            return f"rnapolis is not installed in {python}"
        tail = (done.stderr.strip().splitlines() or ["failed"])[-1]
        return f"RNApolis could not read it: {tail[:160]}"
    try:
        return {frozenset(tuple(x) for x in pair) for pair in json.loads(done.stdout)}
    except (ValueError, TypeError) as err:
        return f"RNApolis gave no pairs: {err}"


def deposited_pairs(path) -> set[frozenset]:
    """foldmap's canonical and wobble pairs (nucleic.pair_nucleotides) on the model as deposited, as RNApolis
    names them."""
    from .io import _nucleic_strands, read_model
    from .nucleic import pair_nucleotides

    model, _ = read_model(Path(path), "asu")
    found = _nucleic_strands(model)
    pairs = pair_nucleotides([run for _, run in found])

    def key(s, i):
        lab = found[s][0].labels[i]
        return (lab.chain, lab.seq, lab.icode)

    return {frozenset((key(p.a, p.i), key(p.b, p.j))) for p in pairs}


def pairs_check(path, theirs: set[frozenset] | str | None = None, nucleic: bool = True) -> dict:
    """D3: foldmap's pairs on the deposited model against RNApolis' (`theirs`; None: run it, rnapolis_pairs):
    precision, recall and F1, passed at nacheck.THRESHOLDS["pair_f1"]."""
    if not nucleic:
        return result(NOT_APPLICABLE, "no nucleic acid")
    theirs = rnapolis_pairs(path) if theirs is None else theirs
    if isinstance(theirs, str):
        return result(UNAVAILABLE, theirs)
    from .nacheck import THRESHOLDS

    ours = deposited_pairs(path)
    tp = len(ours & theirs)
    precision = tp / len(ours) if ours else 1.0
    recall = tp / len(theirs) if theirs else 1.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    limit = THRESHOLDS["pair_f1"]
    measures = {
        "pairs": len(ours),
        "reference_pairs": len(theirs),
        "both": tp,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
    }
    detail = f"F1 {f1:.3f} ({tp} of {len(ours)} pairs, {len(theirs)} in RNApolis; limit {limit})"
    return result(PASS if f1 >= limit else FAIL, detail, **measures)


def deep_checks(spec, drawn, path, draw=None) -> dict[str, dict]:
    """D1-D3 of the figure `drawn` of `spec` (the structure file `path`)."""
    return {
        "D1": repeat_run(spec, drawn, draw),
        "D2": dssp_check(path, protein=len(drawn.bb) > 0),
        "D3": pairs_check(path, nucleic=drawn.layout.nucleic is not None),
    }


def deep_lines(deep: dict[str, dict]) -> list[str]:
    """One line per deep check (pass, fail, unavailable and why, not applicable), its detail indented."""
    lines = [
        "deep check kinds: "
        + "; ".join(
            f"{code} {'structural consistency' if k == CONSISTENCY else 'external validation'}"
            for code, k in DEEP_KIND.items()
        )
    ]
    for code, title in DEEP.items():
        got = deep[code]
        status = got["status"]
        if status == UNAVAILABLE:
            lines.append(f"{code} {title}: unavailable ({got['detail']})")
        elif status == NOT_APPLICABLE:
            lines.append(f"{code} {title}: does not apply ({got['detail']})")
        else:
            lines.append(f"{code} {title}: {status}")
            if got["detail"]:
                lines.append(f"  {got['detail']}")
    return lines


# --------------------------------------------------------------------------------------------- the report


def _science_entry(code: str, found: list | None, judged: dict | None) -> dict:
    if found is None:
        return {"status": UNAVAILABLE, "reason": "not yet measured", "findings": []}
    if not found and judged is not None and judged.get(code) == 0:
        return {"status": NOT_JUDGED, "reason": NOTHING[code], "findings": []}
    return {"status": FAIL if found else PASS, "findings": list(found)}


def _criterion_entry(got: dict) -> dict:
    from .criteria import NOT_JUDGED as CRIT_NOT_JUDGED
    from .criteria import outcome, unmeasured

    got = unmeasured(got or {"pass": None})
    detail = str(got.get("detail") or "")
    if got["pass"] is True:
        return {"status": PASS, "detail": detail}
    if got["pass"] is False:
        return {"status": FAIL, "detail": detail}
    if got.get("here"):
        return {"status": UNAVAILABLE, "reason": str(got["here"])}
    if outcome(got) == CRIT_NOT_JUDGED:
        return {"status": NOT_JUDGED, "reason": detail}
    return {"status": NOT_APPLICABLE, "detail": detail}


def report(
    structure,
    found: dict,
    like: dict | None = None,
    deep: dict | None = None,
    floor: dict | None = None,
    splits: list | None = None,
) -> dict:
    """The machine-readable check report (module docstring) of `foldmap check`: `found` criteria.check_figure's,
    `like` the H checks (criteria.like_checks) when drawn like a reference, `deep` D1-D3 (deep_checks) when run,
    `floor` R1 (floor.check: the readability floor) when judged, `splits` the inferred domain splits with their
    confidence (domain_confidence.splits_of), S5's measures when there are any (said, never judged: its status is
    unchanged)."""
    from . import __version__
    from . import floor as floor_mod
    from .criteria import CRITERIA, LIKE_CHECKS
    from .domain_confidence import S5_REASON
    from .domain_confidence import measures as split_measures

    kinds = {**check_kinds(CRITERIA), **DEEP_KIND}
    checks = []
    for code, title in SCIENCE.items():
        entry = _science_entry(code, found["science"].get(code), found.get("science_judged"))
        if code == "S5" and splits:  # plan 28 Task 5: an optional key of the v1 shape, the status unchanged
            entry["measures"] = split_measures(splits)
            if entry["status"] == UNAVAILABLE:  # what was not measured, and what the measures are (ruling Q43 M1)
                entry["reason"] = S5_REASON
        checks.append({"code": code, "group": "science", "title": title, "kind": kinds[code], **entry})
    for code, title in CRITERIA.items():
        entry = _criterion_entry(found["criteria"].get(code))
        checks.append({"code": code, "group": "criteria", "title": title, "kind": kinds[code], **entry})
    if floor is not None:  # R1, readability like A1-A20; its reasons and the options that draw less as measures
        measures = {k: floor[k] for k in ("reasons", "suggestions", "measures")}
        entry = {"status": floor["status"], "detail": floor["detail"], "measures": measures}
        head = {"code": floor_mod.CODE, "group": FLOOR, "title": floor_mod.TITLE, "kind": READABILITY}
        checks.append({**head, **entry})
    for code, title in LIKE_CHECKS.items() if like is not None else ():
        got = like.get(code) or {"pass": None, "detail": ""}
        detail = str(got.get("detail") or "")
        status = {True: PASS, False: FAIL}.get(got["pass"], NOT_JUDGED)
        entry = {"status": status, ("reason" if status == NOT_JUDGED else "detail"): detail}
        checks.append({"code": code, "group": "like", "title": title, "kind": LIKE, **entry})
    for code, got in (deep or {}).items():
        extra = {k: v for k, v in got.items() if k not in ("status", "detail")}
        entry = {"status": got["status"], "detail": got["detail"], **({"measures": extra} if extra else {})}
        if got["status"] == UNAVAILABLE:
            entry = {"status": UNAVAILABLE, "reason": got["detail"]}
        checks.append({"code": code, "group": "deep", "title": DEEP[code], "kind": DEEP_KIND[code], **entry})
    counts = {s: sum(c["status"] == s for c in checks) for s in STATUSES}
    verdict = FAIL if counts[FAIL] else "incomplete" if counts[UNAVAILABLE] else PASS
    return {
        "schema": SCHEMA,
        "version": VERSION,
        "structure": str(structure),
        "foldmap": __version__,
        "deep": deep is not None,
        "verdict": verdict,
        "counts": counts,
        "kinds": {**KINDS, LIKE: LIKE_KIND},
        "checks": checks,
    }
