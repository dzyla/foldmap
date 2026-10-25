# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Check a report's source stamps against the commit they name (plan 28 Task 1, ruling Q11).

    python scripts/paper/verify_stamp.py REPORT [--commit SHA]

A run on a machine without .git records only the code fingerprint and the commit the runner asserts
(FOLDMAP_GIT_COMMIT). This makes the claim checkable: for every source stamp in REPORT (.json or .json.gz; the
report's own, each entry's, and those nested in a paper step's run stamp), it exports src/foldmap at the commit the
stamp names (`git archive` into a temporary folder; --commit, else a run stamp's harness_commit, comes first) and
computes its fingerprint with the same function the stamp used (app_jobs.code_fingerprint, from this checkout). The
stamp is verified when the fingerprints agree and, where the stamp names its script, the script's text at that commit
hashes to the stamp's SHA-256. Needs this checkout's .git. Exit 0 only when every stamp is verified.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO, repo_module  # noqa: E402


def stamps(doc, commit: str | None = None) -> list[tuple[dict, str | None]]:
    """Every source stamp in a report, each with the commit it is to be checked against: the stamp's own git commit,
    or the harness commit of the run stamp holding it (the code was build/tree, exported at that commit)."""
    out: list[tuple[dict, str | None]] = []

    def walk(node, harness):
        if isinstance(node, dict):
            harness = node.get("harness_commit") or harness
            st = node.get("source_stamp")
            if isinstance(st, dict) and "code_fingerprint" in st:
                out.append((st, commit or harness or (st.get("git") or {}).get("commit")))
            for v in node.values():
                walk(v, harness)
        elif isinstance(node, list):
            for v in node:
                walk(v, harness)

    walk(doc, None)
    return out


def tree_fingerprint(commit: str, repo: Path = REPO) -> str:
    """The fingerprint of src/foldmap as committed at `commit` (git archive), by app_jobs.code_fingerprint."""
    tar = subprocess.run(
        ["git", "-C", str(repo), "archive", "--format=tar", commit, "src/foldmap"], capture_output=True, check=True
    ).stdout
    with tempfile.TemporaryDirectory() as tmp, tarfile.open(fileobj=io.BytesIO(tar)) as t:
        t.extractall(tmp, filter="data")
        return repo_module("app_jobs").code_fingerprint(Path(tmp) / "src" / "foldmap")


def script_sha256(commit: str, name: str, repo: Path = REPO) -> str | None:
    done = subprocess.run(["git", "-C", str(repo), "show", f"{commit}:{name}"], capture_output=True)
    return hashlib.sha256(done.stdout.replace(b"\r\n", b"\n")).hexdigest() if done.returncode == 0 else None


def verify(stamp: dict, commit: str | None, repo: Path = REPO, _seen: dict | None = None) -> dict:
    """{commit, verified, why}: does the stamp's code (and script) match `commit`'s?"""
    if not commit:
        return {"commit": None, "verified": False, "why": "the stamp names no commit (no git, nothing asserted)"}
    try:
        full = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--verify", f"{commit}^{{commit}}"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()  # fmt: skip
    except subprocess.CalledProcessError:
        return {"commit": commit, "verified": False, "why": f"{commit} is not a commit of this repository"}
    seen = _seen if _seen is not None else {}
    if full not in seen:
        seen[full] = tree_fingerprint(full, repo)
    want = seen[full]
    if stamp.get("code_fingerprint") != want:
        why = f"code fingerprint {stamp.get('code_fingerprint')} is not {full[:7]}'s src/foldmap ({want})"
        return {"commit": full, "verified": False, "why": why}
    script = stamp.get("script") or {}
    is_run_tree = (stamp.get("git") or {}).get(
        "commit"
    ) != commit  # a run stamp: the script is the checkout's, not the tree's
    if script.get("name") and script.get("sha256") and not is_run_tree:
        got = script_sha256(full, script["name"], repo)
        if got != script["sha256"]:
            return {"commit": full, "verified": False, "why": f"{script['name']} differs from {full[:7]}'s"}
    return {"commit": full, "verified": True, "why": ""}


def verify_report(doc, commit: str | None = None, repo: Path = REPO) -> list[dict]:
    """Every distinct (fingerprint, script, commit) of a report, verified once."""
    seen: dict = {}
    done: dict[tuple, dict] = {}
    for st, c in stamps(doc, commit):
        key = (st.get("code_fingerprint"), (st.get("script") or {}).get("sha256"), c)
        if key not in done:
            done[key] = verify(st, c, repo, seen) | {"code_fingerprint": st.get("code_fingerprint")}
    return list(done.values())


def load(path: Path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == ".gz" else data)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("report", type=Path)
    ap.add_argument("--commit", help="check every stamp against this commit instead of the one it names")
    args = ap.parse_args(argv)
    results = verify_report(load(args.report), args.commit)
    if not results:
        print(f"{args.report}: no source stamp")
        return 1
    for r in results:
        mark = "verified" if r["verified"] else f"NOT VERIFIED: {r['why']}"
        print(f"{r['code_fingerprint']} at {str(r['commit'])[:12]}: {mark}")
    return 0 if all(r["verified"] for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
