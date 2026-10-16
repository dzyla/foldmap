# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Fill the domain calibration cache (needs network once): the structures of tests/data/domains/cath_reference.json
and cath_fresh.json, as RCSB mmCIF files, in $FOLDMAP_DOMAIN_CACHE (default ~/.cache/foldmap-domains), outside the
repository (46 files, about 50 MB: too big for fixtures). tests/test_domain_score.py's end-to-end check runs only
when they are all there; its CI checks run on the committed inputs (scripts/make_domain_fixture.py).

    python scripts/make_domain_cache.py               # download what is missing
    python scripts/make_domain_cache.py --from DIR    # copy DIR/<ID>.cif instead (e.g. ~/.cache/foldmap)
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

DOMAINS = Path(__file__).resolve().parent.parent / "tests" / "data" / "domains"
REFERENCES = (DOMAINS / "cath_reference.json", DOMAINS / "cath_fresh.json")


def cache_dir() -> Path:
    return Path(os.environ.get("FOLDMAP_DOMAIN_CACHE") or Path.home() / ".cache" / "foldmap-domains")


def ids() -> list[str]:
    docs = [json.loads(path.read_text(encoding="utf-8")) for path in REFERENCES]
    return sorted({e["pdb"] for doc in docs for k in ("chains", "single", "literature") for e in doc.get(k, [])})


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from", dest="src", type=Path, help="copy DIR/<ID>.cif instead of downloading")
    args = ap.parse_args(argv)
    out = cache_dir()
    out.mkdir(parents=True, exist_ok=True)
    for pdb in ids():
        dest = out / f"{pdb}.cif"
        if dest.is_file():
            continue
        if args.src is not None:
            shutil.copyfile(args.src / f"{pdb}.cif", dest)
        else:
            from foldmap.fetch import fetch

            shutil.copyfile(fetch(pdb), dest)
        print(dest)
    print(f"{len(ids())} structures in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
