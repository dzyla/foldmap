# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Which code made a measurement (plan 28 Task 1, review 1.1): the stamp every campaign script writes into its report
when it writes it, never afterwards, so each reported number can be traced to the code that produced it.

source_stamp() returns, as plain JSON:

    foldmap            the package version (foldmap.__version__)
    code_fingerprint   app_jobs.code_fingerprint(): a hash of every .py source of the package on disk, taken with or
                       without git (a checkout synced to another machine has no .git; this still names its code)
    git                where the package sits in a git checkout: HEAD's commit, and whether src/, scripts/ or
                       pyproject.toml differ from it (dirty, with the paths, at most DIRTY_LISTED); None without .git
    script             the campaign script that wrote the report (given as `script`): its path from the repository
                       root and the SHA-256 of its text (line endings normalised), as scripts/ is not in the package
    python, platform   platform.python_version() and implementation; platform.platform() and machine
    packages           the installed versions of foldmap's runtime dependencies (RUNTIME, as pyproject.toml lists
                       them; None for one not installed)

Never raises on a missing git, an unreadable file or an absent package: the field says so instead."""

from __future__ import annotations

import hashlib
import platform
import subprocess
from pathlib import Path

RUNTIME = ("gemmi", "numpy", "scipy", "networkx", "matplotlib", "PyYAML")  # pyproject.toml [project] dependencies
CODE = ("src", "scripts", "pyproject.toml")  # what a dirty checkout can differ in that changes a measurement
DIRTY_LISTED = 50


def _repo(package: Path) -> Path | None:
    """The checkout holding the package (src/foldmap -> its root), when it is one with .git."""
    root = package.parent.parent
    return root if (root / ".git").exists() else None


def _git(repo: Path) -> dict:
    def run(*args: str) -> str:
        done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, timeout=60)
        return done.stdout

    try:
        commit = run("rev-parse", "HEAD").strip()
        changed = [line[3:] for line in run("status", "--porcelain", "--", *CODE).splitlines() if line.strip()]
    except (OSError, subprocess.SubprocessError) as err:
        return {"commit": None, "error": f"{type(err).__name__}: {err}"[:300]}
    return {"commit": commit, "dirty": bool(changed), "dirty_paths": changed[:DIRTY_LISTED]}


def _script(script: str | Path, repo_root: Path) -> dict:
    path = Path(script).resolve()
    try:
        name = path.relative_to(repo_root).as_posix()
    except ValueError:
        name = path.name
    try:
        digest = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    except OSError as err:
        return {"name": name, "sha256": None, "error": type(err).__name__}
    return {"name": name, "sha256": digest}


def _version(dist: str) -> str | None:
    from importlib import metadata

    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return None


def source_stamp(script: str | Path | None = None, package: Path | None = None) -> dict:
    """The stamp described above. `script`: the campaign script writing the report (its __file__). `package`: the
    foldmap package directory to name (default: this one; tests point it elsewhere)."""
    from . import __version__
    from .app_jobs import code_fingerprint

    package = Path(__file__).resolve().parent if package is None else Path(package).resolve()
    repo = _repo(package)
    out: dict = {
        "foldmap": __version__,
        "code_fingerprint": code_fingerprint(package),
        "git": _git(repo) if repo is not None else None,
    }
    if script is not None:
        out["script"] = _script(script, package.parent.parent)
    out["python"] = {"version": platform.python_version(), "implementation": platform.python_implementation()}
    out["platform"] = {"platform": platform.platform(), "machine": platform.machine()}
    out["packages"] = {name: _version(name) for name in RUNTIME}
    return out


def stamp_line(stamp: dict) -> str:
    """One line for a Markdown report: the commit (or that there is no git), the fingerprint, Python and packages."""
    git = stamp.get("git") or {}
    if git.get("commit"):
        dirty = git.get("dirty_paths", [])
        where = f"git {git['commit'][:12]}" + (f" (dirty: {', '.join(dirty[:5])})" if git.get("dirty") else "")
    else:
        where = "no git checkout"
    script = stamp.get("script") or {}
    packages = ", ".join(f"{k} {v}" for k, v in (stamp.get("packages") or {}).items())
    return (
        f"foldmap {stamp.get('foldmap')}, {where}, code fingerprint {stamp.get('code_fingerprint')}"
        + (f", {script['name']} sha256 {str(script.get('sha256'))[:12]}" if script else "")
        + f"; Python {(stamp.get('python') or {}).get('version')} on {(stamp.get('platform') or {}).get('platform')}"
        + f"; {packages}"
    )
