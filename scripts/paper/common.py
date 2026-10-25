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


def stamp() -> dict:
    """What extract.py recorded: the harness commit and every source's provenance."""
    path = BUILD / "provenance.json"
    if not path.exists():
        raise SystemExit("build/provenance.json missing: run extract.py first")
    return json.loads(path.read_text(encoding="utf-8"))


def run_stamp(script: str) -> dict:
    """The source identity of a step that runs foldmap (domains.py, dssp_mkdssp.py, membrane_opm.py), written into its
    output as it is written (plan 28 Task 1): the harness commit build/tree was exported at (exact: the tree is
    `git archive` of it), and foldmap.source_stamp of the foldmap imported (the tree's, which the step asserts) when
    that code has it; a tree from before plan 28 Task 1 has not, and says so."""
    out: dict = {"harness_commit": None}
    prov = BUILD / "provenance.json"
    if prov.exists():
        out["harness_commit"] = json.loads(prov.read_text(encoding="utf-8")).get("harness_commit")
    try:
        from foldmap.source_stamp import source_stamp
    except ImportError:
        import platform

        out["source_stamp"] = None
        out["note"] = "the foldmap under measurement predates foldmap.source_stamp"
        out["python"] = platform.python_version()
        out["platform"] = platform.platform()
        return out
    out["source_stamp"] = source_stamp(script)
    return out


# ---- statistics ----------------------------------------------------------------------------------------------------


def percentile(xs: list[float], q: float) -> float:
    """Linear interpolation between closest ranks (numpy.percentile's default)."""
    import math

    xs = sorted(xs)
    k = (len(xs) - 1) * q / 100
    lo, hi = math.floor(k), math.ceil(k)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def median(xs: list[float]) -> float:
    return percentile(xs, 50)


# ---- formatting for LaTeX -------------------------------------------------------------------------------------------


def tex_int(n: int) -> str:
    """12345 -> 12{,}345 (LaTeX thousands separator without extra space)."""
    s = f"{n:,}"
    return s.replace(",", "{,}")


def tex_frac(x: float, places: int = 3) -> str:
    return f"{x:.{places}f}"


def tex_escape(s: str) -> str:
    return s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_").replace("#", r"\#")
