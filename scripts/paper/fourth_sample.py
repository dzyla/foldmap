# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The fourth random-100 sample (plan 28 Task 1, review 3.9): drawn and run once, at the frozen release, and its count
is the one the paper quotes. The protocol is docs/validation/fourth-sample-protocol.md; every choice it fixes is a
constant here, set before any draw (2026-10-25).

    python scripts/paper/fourth_sample.py plan                    the fixed choices and the excluded ids; no network,
                                                                  writes nothing
    python scripts/paper/fourth_sample.py draw --release TAG      draw the sample (network: the RCSB search service and
                                                                  downloads) into OUT; refuses when OUT exists

The draw is scripts/random100.py's sampler unchanged (its strata, quotas, size bins, search pools and LOOK) with:
  SEED      20261025, fixed here, never re-drawn
  DATE      2026-10-13, the release-date cutoff of samples 1-3 (the same population, less their ids)
  EXCLUDED  every id of samples 1-3 (their entries and their CA-only / P-only replacements: random100's --exclude),
            read from the three committed sample.json files, which must still hash to SAMPLE_SHA256
sample.json then carries a "protocol" record: the seed, cutoff and rule, each excluded file with its SHA-256 and id
count, the release tag, and foldmap.source_stamp of the code that drew. Its run, at the same release, with every view:

    python scripts/random100.py run --views all --jobs N --sample docs/validation/random100-fourth/sample.json \\
        --out docs/validation/random100-fourth
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261025
DATE = "2026-10-13"
OUT = ROOT / "docs" / "validation" / "random100-fourth" / "sample.json"
PROTOCOL = "docs/validation/fourth-sample-protocol.md"
# samples 1-3 (tuned, held-out, third) as committed when the protocol was written: SHA-256 of the text, LF line ends
SAMPLE_SHA256 = {
    "docs/validation/random100/sample.json": "075deb5b4aa2ce50782171d1e3fecd96f1b5ceca1599664772b8a3a13709022d",
    "docs/validation/random100-heldout/sample.json": "683be66f3c0630e707960b359fde5da792fb546066938d0b14bbe0773c0f4cbd",
    "docs/validation/random100-third/sample.json": "f535ce177f42383d2cc443542df0718358d1716779aed34d5b0f6744c12b4eb2",
}
RULE = (
    "scripts/random100.py's stratified sampler, unchanged at the release (strata, quotas, size bins, search pools, "
    f"LOOK), seed {SEED}, release-date cutoff {DATE}; no id of samples 1-3 (entries and CA-only / P-only "
    "replacements) may be drawn; drawn once, at the frozen release, before any foldmap result on it is seen"
)


def random100():
    """scripts/random100.py as a module (it puts this checkout's src/ first only when run as a script)."""
    if str(ROOT / "src") not in sys.path and (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    spec = importlib.util.spec_from_file_location("random100", ROOT / "scripts" / "random100.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def excluded(r100, root: Path = ROOT) -> tuple[list[str], list[dict]]:
    """The ids samples 1-3 hold (sorted), and per sample file its path, SHA-256 and id count. Raises SystemExit when a
    file is missing or no longer hashes to SAMPLE_SHA256: the exclusion must be of the samples as committed."""
    ids: set[str] = set()
    files = []
    for rel, want in SAMPLE_SHA256.items():
        path = root / rel
        if not path.is_file():
            raise SystemExit(f"{rel} is missing: the fourth sample excludes samples 1-3 as committed")
        got = sha256_text(path)
        if got != want:
            raise SystemExit(f"{rel} hashes to {got}, not {want}: a sample of 1-3 changed after the protocol")
        mine = {e["id"] for e in r100._sample_entries(path)}
        ids |= mine
        files.append({"path": rel, "sha256": got, "ids": len(mine)})
    return sorted(ids), files


def cmd_plan(args) -> int:
    r100 = random100()
    ids, files = excluded(r100)
    print(f"seed {SEED}, release-date cutoff {DATE}, output {OUT.relative_to(ROOT).as_posix()}")
    for f in files:
        print(f"  excludes {f['path']}: {f['ids']} ids (sha256 {f['sha256'][:16]})")
    print(f"  {len(ids)} ids excluded in all")
    print(f"rule: {RULE}")
    print("drawn:", "yes, " + str(OUT) if OUT.exists() else "not yet (it is drawn once, at the frozen release)")
    return 0


def cmd_draw(args) -> int:
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f"{out} exists: the fourth sample is drawn once; it is never re-drawn")
    r100 = random100()
    ids, files = excluded(r100)
    doc = r100.sample(SEED, DATE, ids)
    drawn = {e["id"] for e in doc["entries"]}
    if drawn & set(ids):
        raise SystemExit(f"the draw holds ids of samples 1-3: {sorted(drawn & set(ids))}")
    from foldmap.source_stamp import source_stamp

    doc["protocol"] = {
        "name": "fourth random-100 sample",
        "document": PROTOCOL,
        "seed": SEED,
        "query_date": DATE,
        "rule": RULE,
        "excluded": files,
        "excluded_ids": len(ids),
        "release": args.release,
        "source_stamp": source_stamp(__file__),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(r100.dump(doc), encoding="utf-8")
    print(r100.quota_table(doc))
    short = [k for k, q in doc["quotas"].items() if q["reached"] < q["target"]]
    print(f"{len(doc['entries'])} entries written to {out}" + (f"; short: {', '.join(short)}" if short else ""))
    return 1 if short or len(doc["entries"]) != 100 else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan", help="print the fixed choices and the excluded ids (no network, writes nothing)")
    p.set_defaults(func=cmd_plan)
    d = sub.add_parser("draw", help="draw the sample once, at the frozen release (network)")
    d.add_argument("--release", required=True, help="the frozen release's tag, recorded in the sample")
    d.add_argument("--out", default=str(OUT))
    d.set_defaults(func=cmd_draw)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
