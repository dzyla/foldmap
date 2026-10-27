"""Shrink docs/validation/inventory/freeze-baseline.json after a special case was removed (plan 28 Task 9.5, Q55).

    python scripts/freeze_baseline.py

Each count becomes the smaller of the baseline's and the table's, and an entry that reaches zero is dropped. The
baseline can only shrink here: a row the table has beyond the baseline's count (a new special case) is never
admitted, so this is not a way to make tests/test_freeze.py pass. See CONTRIBUTING.md."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
BASELINE = ROOT / "docs" / "validation" / "inventory" / "freeze-baseline.json"
ANCHOR = re.compile(r"@[0-9a-f]{8}(?:#\d+)?$")
SPECIAL = ("named-entry", "measured", "unknown")


def counts(ids) -> dict[str, int]:
    out: dict[str, int] = {}
    for i in ids:
        k = ANCHOR.sub("", i)
        out[k] = out.get(k, 0) + 1
    return out


def shrink(old: dict[str, int], now: dict[str, int]) -> dict[str, int]:
    """old with every count lowered to now's where now's is smaller; never adds a key or raises a count."""
    return {k: min(n, now.get(k, 0)) for k, n in sorted(old.items()) if min(n, now.get(k, 0)) > 0}


def shrunk_baseline(doc: dict, rows: dict[str, dict], origin_kind) -> dict:
    special = counts(i for i, r in rows.items() if origin_kind(r["origin"]) in SPECIAL)
    new = {**doc, "special": shrink(doc["special"], special), "all": shrink(doc["all"], counts(rows))}
    new["count"] = {"special": sum(new["special"].values()), "all": sum(new["all"].values())}
    return new


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args(argv)
    import inventory as inv

    from foldmap.source_stamp import source_stamp

    rows = inv.read_table((BASELINE.parent / "inventory.md").read_text(encoding="utf-8"))
    doc = json.loads(BASELINE.read_text(encoding="utf-8"))
    new = shrunk_baseline(doc, rows, inv.origin_kind)
    new["source_stamp"] = source_stamp(__file__)
    BASELINE.write_text(json.dumps(new, indent=1) + "\n", encoding="utf-8", newline="\n")
    before, after = doc["count"], new["count"]
    print(f"special {before['special']} -> {after['special']}, all {before['all']} -> {after['all']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
