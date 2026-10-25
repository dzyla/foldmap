# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Which code made a measurement (plan 28 Task 1, review 1.1): the stamp every campaign script writes into its report,
taken when the script starts and written with the report, never added afterwards, so each reported number can be
traced to the code that produced it.

source_stamp(script) returns, as plain JSON (`script` is the writing script's __file__):

    foldmap            the package version (foldmap.__version__)
    code_fingerprint   app_jobs.code_fingerprint() of the foldmap package imported, the code that ran (installed or
                       not): a hash of its .py sources, taken with or without git
    package            "checkout" when that package is the script's checkout's src/foldmap, else "installed"
    checkout_fingerprint, installed_matches_checkout
                       the same hash of the checkout's src/foldmap, and whether the two agree (None without one): an
                       installed package that differs from the checkout is never taken for the checkout's commit
    git                the script's checkout (found from the script's own path, never from where foldmap is
                       imported): HEAD's commit and whether src/, scripts/, pyproject.toml or tests/data differ from it
                       (dirty, with the paths, at most DIRTY_LISTED), "source": "git". With no .git (a copy synced to a
                       runner), the commit and dirty state the runner asserts in FOLDMAP_GIT_COMMIT / FOLDMAP_GIT_DIRTY,
                       "source": "asserted by runner" (scripts/paper/verify_stamp.py checks such a claim against the
                       commit's tree). None when neither is there
    script             the script's path from its checkout and the SHA-256 of its text (line endings normalised), as
                       scripts/ is not in the package
    python, platform   platform.python_version() and implementation; platform.platform() and machine
    packages           the installed versions of foldmap's runtime dependencies (RUNTIME, as pyproject.toml lists
                       them; None for one not installed)

code_mix(stamps) summarises the stamps of a report's entries: a report whose entries were measured by different code
says so ("mixed_code": true, with each code and its count), and the paper pipeline refuses to quote it.

Never raises on a missing git, an unreadable file or an absent package: the field says so instead."""

from __future__ import annotations

import hashlib
import os
import platform
import subprocess
from collections import Counter
from pathlib import Path

RUNTIME = ("gemmi", "numpy", "scipy", "networkx", "matplotlib", "PyYAML")  # pyproject.toml [project] dependencies
CODE = ("src", "scripts", "pyproject.toml", "tests/data")  # what a dirty checkout can differ in that changes a result
DIRTY_LISTED = 50
ENV_COMMIT, ENV_DIRTY = "FOLDMAP_GIT_COMMIT", "FOLDMAP_GIT_DIRTY"  # set by a runner that syncs a tree without .git
ASSERTED = "asserted by runner"


def checkout_of(script: str | Path) -> Path | None:
    """The checkout a script belongs to: the nearest folder above it holding .git or pyproject.toml."""
    here = Path(script).resolve().parent
    for d in (here, *here.parents):
        if (d / ".git").exists() or (d / "pyproject.toml").is_file():
            return d
    return None


def _git(repo: Path) -> dict:
    def run(*args: str) -> str:
        done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, timeout=60)
        return done.stdout

    try:
        commit = run("rev-parse", "HEAD").strip()
        changed = [line[3:] for line in run("status", "--porcelain", "--", *CODE).splitlines() if line.strip()]
    except (OSError, subprocess.SubprocessError) as err:
        return {"commit": None, "source": "git", "error": f"{type(err).__name__}: {err}"[:300]}
    return {"commit": commit, "dirty": bool(changed), "dirty_paths": changed[:DIRTY_LISTED], "source": "git"}


def _asserted() -> dict | None:
    commit = os.environ.get(ENV_COMMIT, "").strip()
    if not commit:
        return None
    raw = os.environ.get(ENV_DIRTY, "").strip().lower()
    dirty = {"0": False, "false": False, "no": False, "1": True, "true": True, "yes": True}.get(raw)
    return {"commit": commit, "dirty": dirty, "source": ASSERTED}


def _script(path: Path, checkout: Path | None) -> dict:
    try:
        name = path.relative_to(checkout).as_posix() if checkout is not None else path.name
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


def source_stamp(script: str | Path, package: Path | None = None) -> dict:
    """The stamp described above. `script`: the campaign script writing the report (its __file__). `package`: the
    foldmap package directory that ran (default: the one imported; tests point it elsewhere)."""
    from . import __version__
    from .app_jobs import code_fingerprint

    package = Path(__file__).resolve().parent if package is None else Path(package).resolve()
    path = Path(script).resolve()
    checkout = checkout_of(path)
    fingerprint = code_fingerprint(package)
    tree = checkout / "src" / "foldmap" if checkout is not None else None
    tree_fingerprint = code_fingerprint(tree.resolve()) if tree is not None and tree.is_dir() else None
    out: dict = {
        "foldmap": __version__,
        "code_fingerprint": fingerprint,
        "package": "checkout" if tree is not None and tree.resolve() == package else "installed",
        "checkout_fingerprint": tree_fingerprint,
        "installed_matches_checkout": None if tree_fingerprint is None else tree_fingerprint == fingerprint,
    }
    if checkout is not None and (checkout / ".git").exists():
        out["git"] = _git(checkout)
    else:
        out["git"] = _asserted()
    out["script"] = _script(path, checkout)
    out["python"] = {"version": platform.python_version(), "implementation": platform.python_implementation()}
    out["platform"] = {"platform": platform.platform(), "machine": platform.machine()}
    out["packages"] = {name: _version(name) for name in RUNTIME}
    return out


def code_key(stamp: dict | None) -> str:
    """The code a stamp names, as one key: the package fingerprint and the script's hash ("unstamped" without one)."""
    if not stamp:
        return "unstamped"
    return f"{stamp.get('code_fingerprint')}+{str((stamp.get('script') or {}).get('sha256'))[:16]}"


def code_mix(stamps) -> dict:
    """The codes a report's entries were measured with: {"entry_code": {code_key: entries}, "mixed_code": bool}."""
    counts = Counter(code_key(s) for s in stamps)
    return {"entry_code": dict(sorted(counts.items())), "mixed_code": len(counts) > 1}


def stamp_line(stamp: dict) -> str:
    """One line for a Markdown report: the commit (or why there is none), the fingerprint, Python and packages."""
    git = stamp.get("git") or {}
    if git.get("commit"):
        dirty = git.get("dirty_paths", [])
        where = f"git {git['commit'][:12]}"
        if git.get("source") == ASSERTED:
            where += f" ({ASSERTED}, dirty {git.get('dirty')})"
        elif git.get("dirty"):
            where += f" (dirty: {', '.join(dirty[:5])})"
    elif git.get("error"):
        where = f"git failed: {git['error']}"
    else:
        where = "no git checkout"
    if stamp.get("installed_matches_checkout") is False:
        where += " (the foldmap that ran differs from the checkout's src/foldmap)"
    script = stamp.get("script") or {}
    packages = ", ".join(f"{k} {v}" for k, v in (stamp.get("packages") or {}).items())
    return (
        f"foldmap {stamp.get('foldmap')}, {where}, code fingerprint {stamp.get('code_fingerprint')}"
        + (f", {script['name']} sha256 {str(script.get('sha256'))[:12]}" if script else "")
        + f"; Python {(stamp.get('python') or {}).get('version')} on {(stamp.get('platform') or {}).get('platform')}"
        + f"; {packages}"
    )
