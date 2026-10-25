# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Step 5 (local, light): main.tex and supplement.tex with every \\NUM{key} that numbers.json fills replaced by its
value.

    python scripts/paper/fill.py

Writes build/main-filled.tex, build/supplement-filled.tex and build/keys.json:
  filled   keys replaced (with how often)
  pending  keys numbers.json lists as null (left as \\NUM{key}, which still prints [NUM:key] in blue)
  missing  keys the text uses that numbers.json lacks entirely (also left as \\NUM{key})
\\NUM in %-comments is left alone. graphicx is added to the preamble when a filled value needs \\includegraphics.
Exits 1 when a key is missing.
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter

from common import BUILD, MANUSCRIPT

NUM = re.compile(r"\\NUM\{([^}]*)\}")
COMMENT = re.compile(r"(?<!\\)%")


def split_comment(line: str) -> tuple[str, str]:
    m = COMMENT.search(line)
    return (line, "") if not m else (line[: m.start()], line[m.start() :])


def fill(text: str, keys: dict, used: Counter, filled: Counter) -> str:
    out = []
    for line in text.splitlines(keepends=True):
        body, comment = split_comment(line)

        def sub(m):
            key = m.group(1)
            used[key] += 1
            v = keys.get(key)
            if v is None or v.get("value") is None:
                return m.group(0)
            filled[key] += 1
            return v["text"]

        out.append(NUM.sub(sub, body) + comment)
    text = "".join(out)
    if "\\includegraphics" in text and "{graphicx}" not in text:
        text = text.replace("\\usepackage{xcolor}", "\\usepackage{xcolor}\n\\usepackage{graphicx}", 1)
    return text


def main() -> int:
    doc = json.loads((BUILD / "numbers.json").read_text(encoding="utf-8"))
    keys = doc["keys"]
    report = {"harness_commit": doc["_meta"]["harness_commit_short"], "files": {}}
    all_used: Counter = Counter()
    for name in ("main", "supplement"):
        used, filled = Counter(), Counter()
        src = (MANUSCRIPT / f"{name}.tex").read_text(encoding="utf-8")
        header = (
            f"% FILLED by scripts/paper/fill.py from build/numbers.json (harness commit "
            f"{doc['_meta']['harness_commit_short']}); edit {name}.tex, not this file.\n"
        )
        (BUILD / f"{name}-filled.tex").write_text(header + fill(src, keys, used, filled), encoding="utf-8")
        all_used.update(used)
        report["files"][name] = {
            "used": len(used),
            "filled": sorted(filled),
            "pending": sorted(k for k in used if k in keys and keys[k].get("value") is None),
            "missing": sorted(k for k in used if k not in keys),
        }
    report["used"] = len(all_used)
    report["filled"] = sorted(k for k in all_used if k in keys and keys[k].get("value") is not None)
    report["pending"] = {
        k: keys[k].get("pending") for k in sorted(all_used) if k in keys and keys[k].get("value") is None
    }
    report["missing"] = sorted(k for k in all_used if k not in keys)
    report["in_numbers_not_in_text"] = sorted(k for k in keys if k not in all_used)
    (BUILD / "keys.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(
        f"keys used in the text: {report['used']}; filled {len(report['filled'])}, pending {len(report['pending'])}, "
        f"missing {len(report['missing'])}"
    )
    for k, stage in report["pending"].items():
        print(f"  pending {stage:8s} {k}")
    for k in report["missing"]:
        print(f"  MISSING {k}")
    if report["in_numbers_not_in_text"]:
        print("  in numbers.json but not in the text:", ", ".join(report["in_numbers_not_in_text"]))
    return 1 if report["missing"] else 0


if __name__ == "__main__":
    sys.exit(main())
