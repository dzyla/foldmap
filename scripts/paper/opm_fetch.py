# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""M2 membrane vs OPM, step 1 (light network, run here): cache OPM's published orientation for the selected entries.

    python scripts/paper/opm_fetch.py

Selection rule (stated before any result was looked at, plan 14 M2): every entry of stratum S6 ("membrane proteins")
in the three random samples (docs/validation/random100*/sample.json at the harness commit), in sample order, for
which OPM publishes an oriented coordinate file with DUM atoms. OPM holds an entry as a primary structure or as a
secondary representation of one (same protein, another PDB entry); both have their own oriented file, whose DUM
atoms give that entry's planes. An entry OPM does not hold is listed with the reason and left out.
(First wording, 2026-10-18: "primary structures only, OPM publishes files only for those". The search showed that
secondary representations have files too, so the rule was widened before any foldmap result was computed.)

For each entry it writes under ../data/opm/:
  <id>.search.json   the primary_structures search result (opm-back.cc.lehigh.edu/opm-backend)
  <id>.json          the primary structure record: thickness, tilt, classtype, membrane sides, subunit segments
  <id>.pdb.gz        OPM's oriented coordinates (opm-assets.storage.googleapis.com/pdb/<id>.pdb), DUM atoms included
and ../data/opm/selection.json: the rule, the candidates, which were taken and why the others were not.
Files already cached are not fetched again.
"""

from __future__ import annotations

import gzip
import json
import sys
import urllib.request
from datetime import date

from common import MANUSCRIPT, SOURCES

API = "https://opm-back.cc.lehigh.edu/opm-backend"
ASSETS = "https://opm-assets.storage.googleapis.com/pdb"
OUT = MANUSCRIPT / "data" / "opm"
SAMPLES = ("random100", "random100-heldout", "random100-third")
RULE = (
    "every stratum-S6 entry of the three random samples (sample.json order) for which OPM publishes an oriented "
    "coordinate file (primary structure or secondary representation); chosen before any foldmap result was computed"
)


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "foldmap-manuscript-harness (plan 14 M2)"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def cached_json(path, url):
    if not path.exists():
        path.write_bytes(get(url))
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cands = []
    for s in SAMPLES:
        sample = json.loads((SOURCES / "HEAD" / "docs" / "validation" / s / "sample.json").read_text(encoding="utf-8"))
        cands += [(s, e["id"]) for e in sample["entries"] if e["stratum"] == "S6"]
    taken, left = [], []
    for s, pid in cands:
        low = pid.lower()
        found = cached_json(OUT / f"{pid}.search.json", f"{API}/primary_structures?search={low}&sort=&pageSize=20")
        hit = [o for o in found.get("objects", []) if o["pdbid"].lower() == low]
        if not hit:
            sec = cached_json(
                OUT / f"{pid}.secondary.json", f"{API}/secondary_representations?search={low}&sort=&pageSize=20"
            )
            secs = [o for o in sec.get("objects", []) if o.get("pdbid", "").lower() == low]
            primary = secs[0].get("primary_structure", {}).get("pdbid") if secs else None
            why = f"OPM holds it only as a secondary representation of {primary}" if secs else "not in OPM"
            left.append({"sample": s, "id": pid, "why": why})
            continue
        rec = cached_json(OUT / f"{pid}.json", f"{API}/primary_structures/{hit[0]['id']}")
        pdb = OUT / f"{pid}.pdb.gz"
        if not pdb.exists():
            pdb.write_bytes(gzip.compress(get(f"{ASSETS}/{low}.pdb")))
        taken.append(
            {
                "sample": s,
                "id": pid,
                "opm_id": rec["id"],
                "opm_primary": rec["pdbid"].upper(),
                "secondary": rec["pdbid"].lower() != low,
                "classtype": rec["family"]["superfamily"]["classtype"]["name"],
                "type": rec["family"]["superfamily"]["classtype"]["type"]["name"],
            }
        )
        print(pid, s, taken[-1]["classtype"], taken[-1]["type"], flush=True)
    for x in left:
        print("left out:", x["id"], x["sample"], x["why"])
    sel = {
        "rule": RULE,
        "fetched": date.today().isoformat(),
        "api": API,
        "assets": ASSETS,
        "candidates": [{"sample": s, "id": i} for s, i in cands],
        "taken": taken,
        "left_out": left,
    }
    (OUT / "selection.json").write_text(json.dumps(sel, indent=1) + "\n", encoding="utf-8")
    print(f"{len(taken)} taken of {len(cands)} candidates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
