# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Step 1 (local, light): copy every source out of git at one commit and check what each source was measured at.

    python scripts/paper/extract.py [--commit 838c579]

Writes
  build/sources/HEAD/<path>      each campaign output as it is at the harness commit
  build/sources/<sha>/<path>     a source read at an earlier commit (the held-out sample's first run, 1ef0706)
  build/tree/                    src/, tests/data/domains/ and pyproject.toml at the harness commit (domains.py)
  build/provenance.json          the harness commit and, per source: the commit that last wrote the file, the commit
                                 the measurement inside it was made at (as the file records it), whether that commit
                                 is an ancestor of the harness commit, and every code commit (src/, scripts/) between
                                 the two, with its subject.

The check of M1 item 3 lives here: a source measured at an earlier commit than the harness commit is not hidden; its
code changes in between are listed in provenance.json and copied into numbers.json's metadata. It exits 1 when a
source records a commit that is not an ancestor of the harness commit (a number from another line of history).
"""

from __future__ import annotations

import argparse
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
from collections import Counter

from common import BUILD, CAMPAIGN, DEFAULT_COMMIT, REPO, SOURCES, TREE, TREE_PATHS, full_sha, git, load_json

CODE = ("src", "scripts")  # what can change a measurement


def is_ancestor(a: str, b: str) -> bool:
    return subprocess.run(["git", "-C", str(REPO), "merge-base", "--is-ancestor", a, b]).returncode == 0


def code_between(older: str, newer: str) -> list[dict]:
    """Commits after `older` up to `newer` that touch src/ or scripts/, oldest first, with their changed files."""
    if older == newer:
        return []
    out = []
    log = git("log", "--reverse", "--format=%H%x09%s", f"{older}..{newer}", "--", *CODE)
    for line in filter(None, log.splitlines()):
        sha, subject = line.split("\t", 1)
        files = git("diff-tree", "--no-commit-id", "--name-only", "-r", sha, "--", *CODE).splitlines()
        out.append({"commit": sha[:7], "subject": subject, "files": files})
    return out


def recorded_commits(name: str, path, data_path) -> tuple[dict, str]:
    """The commit(s) the measurement in a source was made at, as the source records it, and where that comes from."""
    if name.endswith("_report") and name.startswith("r100"):
        d = load_json(data_path)
        per = Counter(e.get("foldmap_commit") or "unknown" for e in d["entries"])
        return {
            "report": d.get("foldmap_commit"),
            "entries": dict(per),
        }, "report.json foldmap_commit (report and per entry)"
    if name.endswith("_sample"):
        d = load_json(data_path)
        return {"sample": d.get("foldmap_commit")}, "sample.json foldmap_commit (the code that drew the sample)"
    if name == "na_report":
        d = load_json(data_path)
        if isinstance(d, dict):  # stamped as written (plan 28 Task 1): its git commit, or the code fingerprint
            st = d.get("source_stamp") or {}
            commit = (st.get("git") or {}).get("commit")
            if commit:
                return {"report": commit}, "report.json source_stamp.git.commit (stamped as written)"
            return {"report": None, "code_fingerprint": st.get("code_fingerprint")}, (
                "report.json source_stamp: no git where it ran; the code fingerprint names the code"
            )
        # an older report.json carries no commit; review.md's last "## Rerun ... (<date>, <sha>; ...)" heading names it.
        review = (SOURCES / "HEAD" / CAMPAIGN["na_review"][0]).read_text(encoding="utf-8")
        heads = re.findall(r"^## .*?\((\d{4}-\d\d-\d\d), ([0-9a-f]{7,40})[;)]", review, flags=re.M)
        if heads:
            return {
                "report": heads[-1][1]
            }, f"not stamped in report.json; review.md's last rerun heading ({heads[-1][0]})"
        return {"report": None}, "not stamped in report.json and no rerun heading in review.md"
    return {}, ""


def after_the_run(review_text: str) -> list[str]:
    """Lines of a review that say something changed after the run (quoted into the provenance, not interpreted)."""
    keep = []
    for para in re.split(r"\n\s*\n", review_text):
        flat = " ".join(para.split())
        if re.search(r"after the run|stamped|_MAX_NODES_WINDOWS|e5b4bae", flat, flags=re.I):
            keep.append(flat[:900])
    return keep


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--commit", default=DEFAULT_COMMIT)
    args = ap.parse_args()
    head = full_sha(args.commit)

    for d in (SOURCES, TREE):
        shutil.rmtree(d, ignore_errors=True)
    BUILD.mkdir(parents=True, exist_ok=True)

    # the code and fixtures at the harness commit
    tar = git("archive", "--format=tar", head, *TREE_PATHS, binary=True)
    TREE.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(tar)) as t:
        t.extractall(TREE, filter="data")

    sources, problems = {}, []
    for name, (path, at) in CAMPAIGN.items():
        read_at = full_sha(at) if at else head
        if not is_ancestor(read_at, head):
            problems.append(f"{name}: read at {read_at[:7]}, not an ancestor of {head[:7]}")
        out = SOURCES / (at or "HEAD") / path
        out.parent.mkdir(parents=True, exist_ok=True)
        if name == "hom_summary" and not git("ls-tree", "--name-only", read_at, "--", path):
            continue  # a harness commit before plan 13a's held-out run: make_numbers leaves the hom_ keys unset
        out.write_bytes(git("show", f"{read_at}:{path}", binary=True))
        written = git("log", "-1", "--format=%H%x09%ad%x09%s", "--date=short", read_at, "--", path).split("\t")
        sources[name] = {
            "path": path,
            "read_at": read_at[:7],
            "file_written_by": {"commit": written[0][:7], "date": written[1], "subject": written[2]},
        }

    for name, rec in sources.items():
        recorded, how = recorded_commits(name, rec["path"], SOURCES / (CAMPAIGN[name][1] or "HEAD") / rec["path"])
        if not recorded:
            continue
        rec["measured_at"] = recorded
        rec["measured_at_from"] = how
        if "sample" in recorded:  # a sample is a draw, not a measurement: code after it changes nothing in it
            continue
        main_commit = recorded.get("report")
        commits = {c for c in [main_commit, *recorded.get("entries", {})] if c and c != "unknown"}
        rec["code_changes_to_harness_commit"] = {}
        for c in sorted(commits):
            try:
                sha = full_sha(c)
            except subprocess.CalledProcessError:
                problems.append(f"{name}: records commit {c}, unknown to git")
                continue
            if not is_ancestor(sha, head):
                problems.append(f"{name}: records {c[:7]}, not an ancestor of {head[:7]}")
                continue
            changes = code_between(sha, head)
            rec["code_changes_to_harness_commit"][c[:7]] = changes
        rec["code_changed_since_measured"] = any(rec["code_changes_to_harness_commit"].values())

    # what the reviews say happened after the runs (quoted, so the metadata does not hide it)
    reviews = {}
    for sample in ("random100", "random100-heldout", "random100-third"):
        p = f"docs/validation/{sample}/review.md"
        reviews[p] = after_the_run(git("show", f"{head}:{p}"))

    doc = {
        "harness_commit": head,
        "harness_commit_short": head[:7],
        "harness_commit_date": git("show", "-s", "--format=%ad", "--date=iso", head),
        "harness_commit_subject": git("show", "-s", "--format=%s", head),
        "tags_at_commit": git("tag", "--points-at", head).split(),
        "sources": sources,
        "reviews_after_the_run": reviews,
        "problems": problems,
    }
    (BUILD / "provenance.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")

    print(f"harness commit {head[:7]} ({doc['harness_commit_subject'][:70]})")
    for name, rec in sources.items():
        m = rec.get("measured_at", {})
        main_c = m.get("report") or m.get("sample") or "-"
        changes = rec.get("code_changes_to_harness_commit", {})
        n = sum(len(v) for v in changes.values())
        entries = m.get("entries")
        mixed = f" entries {entries}" if entries and len(entries) > 1 else ""
        print(
            f"  {name:24s} read at {rec['read_at']}  written by {rec['file_written_by']['commit']}  "
            f"measured at {str(main_c)[:7]}{mixed}  code commits since: {n}"
        )
        for c, lst in changes.items():
            for ch in lst[:5]:
                print(f"      {c[:7]}..: {ch['commit']} {ch['subject'][:90]}")
            if len(lst) > 5:
                print(f"      {c[:7]}..: ... {len(lst) - 5} more (provenance.json)")
    for p in problems:
        print("PROBLEM:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
