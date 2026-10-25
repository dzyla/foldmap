# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Shared paths, git access and formatting for the paper pipeline (scripts/paper; plan 14 M1, moved into the
repository by plan 28 Task 1).

Everything the pipeline reads comes out of git at one commit (`git show <commit>:<path>`), never from the working
tree, so a number in numbers.json always names the commit whose files it was read from. Only extract.py talks to git;
the other steps read what it wrote under build/ (a copy synced to another machine has no .git).

Paths are relative to this repository. The manuscript folder (its .tex text, its inputs under data/ and figures/data/,
and every output under build/) stays outside version control: by default docs/superpowers/work/manuscript/ of this
checkout (gitignored), or the folder FOLDMAP_MANUSCRIPT names.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import subprocess
from pathlib import Path

HARNESS = Path(__file__).resolve().parent  # scripts/paper
REPO = HARNESS.parents[1]
MANUSCRIPT = Path(os.environ.get("FOLDMAP_MANUSCRIPT") or REPO / "docs" / "superpowers" / "work" / "manuscript")
BUILD = MANUSCRIPT / "build"
SOURCES = BUILD / "sources"  # every source file as it is at the harness commit (extract.py)
TREE = BUILD / "tree"  # src/ and tests/data/domains/ at the harness commit (for the steps that run foldmap)
TABLES = BUILD / "tables"
DEFAULT_COMMIT = "838c579"  # the commit the current numbers were read at; a release passes its own (--commit)

# The campaign outputs the numbers come from: name -> (path in the repository, commit to read it at; None = the
# harness commit). The held-out sample's first run is read at the commit that wrote it (controller ruling
# 2026-10-17: quote the count judged with the corrected checks, 1ef0706).
CAMPAIGN = {
    "r100_tuned_report": ("docs/validation/random100/report.json", None),
    "r100_tuned_sample": ("docs/validation/random100/sample.json", None),
    "r100_held_report": ("docs/validation/random100-heldout/report.json.gz", None),
    "r100_held_sample": ("docs/validation/random100-heldout/sample.json", None),
    "r100_held_first_report": ("docs/validation/random100-heldout/report.json.gz", "1ef0706"),
    "r100_third_report": ("docs/validation/random100-third/report.json.gz", None),
    "r100_third_sample": ("docs/validation/random100-third/sample.json", None),
    "na_report": ("docs/validation/na/report.json", None),
    "na_review": ("docs/validation/na/review.md", None),
    # plan 13a's held-out homologue pairs (Task 13a-11), measured once; present from the harness commit after it
    "hom_summary": ("docs/validation/homologues/heldout-summary.json", None),
}
# Code read for constants (parsed with ast, never imported locally) and the domain fixtures (domains.py).
TREE_PATHS = (
    "src",
    "tests/data/domains",
    "tests/test_domain_score.py",
    "scripts/random100.py",
    "scripts/na_validate.py",
    "pyproject.toml",
)


def git(*args: str, binary: bool = False):
    out = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True).stdout
    return out if binary else out.decode("utf-8").strip()


def full_sha(commit: str) -> str:
    return git("rev-parse", "--verify", f"{commit}^{{commit}}")


def source_path(name: str) -> Path:
    path, at = CAMPAIGN[name]
    return SOURCES / (at or "HEAD") / path


def load_json(path: Path):
    data = path.read_bytes()
    if path.suffix == ".gz":
        data = gzip.decompress(data)
    return json.loads(data)


def load_source(name: str):
    return load_json(source_path(name))


def na_rows(report) -> list[dict]:
    """The entry rows of a nucleic-acid campaign report.json: the bare list written before plan 28 Task 1, or the
    "entries" of the stamped report written since."""
    return report["entries"] if isinstance(report, dict) else report


def stamped_commit(report) -> str:
    """What a report's own source stamp says about the code that measured it, as words for a source label: the commit
    (as git read it, or as the runner asserted it) and the fingerprint; never a commit from anywhere else (Q11)."""
    st = report.get("source_stamp") if isinstance(report, dict) else None
    if not st:
        return "no source stamp (written before plan 28 Task 1)"
    git = st.get("git") or {}
    fp = f"code fingerprint {st.get('code_fingerprint')}"
    if git.get("commit"):
        return f"at {git['commit'][:7]} ({git.get('source')}, dirty {git.get('dirty')}; {fp})"
    return f"{fp}, no commit recorded"


def stamp() -> dict:
    """What extract.py recorded: the harness commit and every source's provenance."""
    path = BUILD / "provenance.json"
    if not path.exists():
        raise SystemExit("build/provenance.json missing: run extract.py first")
    return json.loads(path.read_text(encoding="utf-8"))


def repo_module(name: str):
    """A module of this checkout's src/foldmap loaded on its own, outside any foldmap package: app_jobs (for
    code_fingerprint, the one fingerprint algorithm) and source_stamp (RUNTIME) import only the standard library at
    module level. Used where the foldmap imported is another tree's (build/tree) or none."""
    import importlib.util
    import sys

    key = f"_paper_repo_{name}"
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, REPO / "src" / "foldmap" / f"{name}.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[key] = mod
        spec.loader.exec_module(mod)
    return sys.modules[key]


def run_stamp(script: str) -> dict:
    """The source identity of a step that runs foldmap from build/tree (domains.py, dssp_mkdssp.py, membrane_opm.py),
    taken when the step starts and written into its output (plan 28 Task 1): the harness commit build/tree was
    exported at (exact: the tree is `git archive` of it; scripts/paper/verify_stamp.py checks the fingerprint against
    it), and foldmap.source_stamp of the foldmap imported (the tree's, which the step asserts). A tree from before
    foldmap.source_stamp gets the same fields this checkout's code computes for it: the tree's fingerprint (the same
    algorithm, app_jobs.code_fingerprint), Python, the platform and the runtime dependencies' versions."""
    import platform

    import foldmap

    out: dict = {"harness_commit": None}
    prov = BUILD / "provenance.json"
    if prov.exists():
        out["harness_commit"] = json.loads(prov.read_text(encoding="utf-8")).get("harness_commit")
    try:
        from foldmap.source_stamp import source_stamp
    except ImportError:
        package = Path(foldmap.__file__).resolve().parent
        runtime = repo_module("source_stamp")
        out["source_stamp"] = {
            "foldmap": getattr(foldmap, "__version__", None),
            "code_fingerprint": repo_module("app_jobs").code_fingerprint(package),
            "package": "tree (it predates foldmap.source_stamp; fields computed by this checkout's code)",
            "git": None,
            "script": runtime._script(Path(script).resolve(), REPO),
            "python": {"version": platform.python_version(), "implementation": platform.python_implementation()},
            "platform": {"platform": platform.platform(), "machine": platform.machine()},
            "packages": {name: runtime._version(name) for name in runtime.RUNTIME},
        }
        return out
    out["source_stamp"] = source_stamp(script)
    return out


def quotable(name: str, report) -> None:
    """Stop when a campaign report may not be quoted: its entries were measured by different code (mixed_code), or it
    was run outside its paper protocol (not_for_paper; plan 28 Task 1, rulings Q13, Q14)."""
    if not isinstance(report, dict):
        return
    if report.get("mixed_code"):
        raise SystemExit(f"{name}: its entries were measured by different code ({report.get('entry_code')}); "
                         "a mixed report is not quoted: rerun it whole with one code")  # fmt: skip
    if report.get("not_for_paper"):
        raise SystemExit(f"{name}: marked not for the paper ({report['not_for_paper']})")


# ---- statistics ----------------------------------------------------------------------------------------------------


def percentile(xs: list[float], q: float) -> float:
    """Linear interpolation between closest ranks (numpy.percentile's default)."""
    xs = sorted(xs)
    k = (len(xs) - 1) * q / 100
    lo, hi = math.floor(k), math.ceil(k)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def median(xs: list[float]) -> float:
    return percentile(xs, 50)


POOLED = (  # the note of a Wilson interval over units that are not independent (common.wilson)
    "pooled over residues, pairs or segments, which are not independent within an entry: the interval is narrower "
    "than one over entries would be"
)
Z95 = 1.959963984540054  # the standard normal's 97.5th percentile: a two-sided 95% interval


def wilson(k: int, n: int, z: float = Z95) -> tuple[float, float] | None:
    """The Wilson score interval of a proportion k of n (95% by default; plan 28 Task 2, ruling Q19): (low, high), or
    None when n is 0. It stays inside [0, 1] and is not degenerate at 0 or n, unlike the normal approximation."""
    if n <= 0:
        return None
    if not 0 <= k <= n:
        raise ValueError(f"wilson: {k} of {n}")
    p, z2 = k / n, z * z
    centre = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    # exact at the ends: the interval reaches 0 at k = 0 and 1 at k = n (in floating point it can stop a hair short)
    return 0.0 if k == 0 else max(0.0, centre - half), 1.0 if k == n else min(1.0, centre + half)


def quartiles(xs: list[float]) -> dict:
    """n, min, first quartile, median, third quartile and max of xs (percentile's interpolation); {"n": 0} for none."""
    if not xs:
        return {"n": 0}
    return {
        "n": len(xs),
        "min": min(xs),
        "q1": percentile(xs, 25),
        "median": percentile(xs, 50),
        "q3": percentile(xs, 75),
        "max": max(xs),
    }


# stems_nearer_frac's histogram: [0, .25), [.25, .5), [.5, .75), [.75, 1) and 1 exactly, the strict pass's own bin
NEARER_BINS = (0.0, 0.25, 0.5, 0.75, 1.0)


def nearer_applies(row: dict) -> bool | None:
    """Whether an NA report row's stems_nearer_frac is a measure (ruling Q5): not where every drawn element touches
    each stem it scores, where the touching elements are compared with themselves (6VTX). nacheck.contacts computes
    contact_auc on exactly the scored stems that leave an element untouched, so a row shows that case as contact_auc
    None with stems_nearer_frac measured. None when stems_nearer_frac is not measured at all. (A row mixing such stems
    with others keeps them in its fraction: the rows hold no per-stem split; plan 28 F7 changes the measure.)"""
    if row.get("stems_nearer_frac") is None:
        return None
    return row.get("contact_auc") is not None


def nearer_summary(rows: list[dict], limit: float = 1.0) -> dict:
    """stems_nearer_frac over NA report rows, next to its strict pass (at least `limit`, nacheck.THRESHOLDS): how
    many rows measure it, those where it is not applicable (nearer_applies; reported as such, never as a pass), and
    over the rest the strict passes, the quartiles and the histogram on NEARER_BINS."""
    measured = [r for r in rows if nearer_applies(r) is not None]
    na = [r["entry"] for r in measured if not nearer_applies(r)]
    vals = [float(r["stems_nearer_frac"]) for r in measured if nearer_applies(r)]
    edges = NEARER_BINS
    bins = [{"from": lo, "to": hi, "entries": sum(lo <= v < hi for v in vals)} for lo, hi in zip(edges, edges[1:])] + [
        {"from": edges[-1], "to": edges[-1], "entries": sum(v >= edges[-1] for v in vals)}
    ]
    return {
        "measured": len(measured),
        "not_applicable": na,
        "applicable": len(vals),
        "strict_limit": limit,
        "strict_pass": sum(v >= limit for v in vals),
        "distribution": quartiles(vals),
        "bins": bins,
    }


# The NA measures by what they are evidence of (nacheck.ROLE; plan 28 Task 2, ruling Q19)
ZERO_BY_CONSTRUCTION = "zero by construction (regression guard)"  # = nacheck.ZERO_BY_CONSTRUCTION
NA_GUARDS = ("angle_median", "angle_max", "joint_error_max")
NA_FIDELITY = ("spearman", "procrustes", "contact_auc", "contact_precision_at_k", "ss_contacts_nearer")


def na_measures(rows: list[dict], nearer_limit: float = 1.0) -> dict:
    """NA report rows summarised by role. regression_guards: angle_median, angle_max and joint_error_max, zero by
    construction on an ordinary projection (stems are drawn along their axes as seen): per measure the rows measuring
    it, how many are 0 and each other row's value (an unrolled frame's, which measures the unrolling). fidelity: the
    quartiles of Spearman, Procrustes and the contact measures; mirrored stems (rows measured, stems in all, rows with
    any); stems_nearer_frac's distribution next to its strict pass (nearer_summary)."""
    guards = {}
    for k in NA_GUARDS:
        got = {r["entry"]: r[k] for r in rows if r.get(k) is not None}
        guards[k] = {
            "tag": ZERO_BY_CONSTRUCTION,
            "measured": len(got),
            "zero": sum(v == 0 for v in got.values()),
            "nonzero": {e: v for e, v in got.items() if v != 0},
        }
    fidelity: dict = {k: quartiles([float(r[k]) for r in rows if r.get(k) is not None]) for k in NA_FIDELITY}
    mirrored = {r["entry"]: r["mirrored"] for r in rows if r.get("mirrored") is not None}
    fidelity["mirrored"] = {
        "measured": len(mirrored),
        "stems": sum(mirrored.values()),
        "entries_with_mirrored": sorted(e for e, v in mirrored.items() if v),
    }
    fidelity["stems_nearer_frac"] = nearer_summary(rows, nearer_limit)
    return {"regression_guards": guards, "fidelity": fidelity}


# ---- formatting for LaTeX -------------------------------------------------------------------------------------------


def tex_int(n: int) -> str:
    """12345 -> 12{,}345 (LaTeX thousands separator without extra space)."""
    s = f"{n:,}"
    return s.replace(",", "{,}")


def tex_frac(x: float, places: int = 3) -> str:
    return f"{x:.{places}f}"


def tex_escape(s: str) -> str:
    return s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_").replace("#", r"\#")
