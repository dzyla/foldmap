# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""M2 DSSP vs mkdssp, step 1 (zyla-lab): foldmap's secondary structure and bridge partners against mkdssp 3.1.2's.

    python scripts/paper/dssp_mkdssp.py   (M2_RERUN=1 scripts/paper/run_all.sh)

Set: every entry of the three random samples whose report measured A10 (dssp.agree not None: the protein entries
mkdssp 3.1.2 reads), the same entries and the same comparison as criterion A10 (criteria.dssp_agreement): the model
as deposited (assembly asu), residues keyed by (auth chain, number, insertion code), mkdssp's bends and blanks read
as '-', a residue standing in for a fused chromophore and the four on either side of it not compared. The A10 counts
in the reports are therefore reproduced exactly (make_numbers.py asserts it), and this run adds what A10 does not
keep: every residue's two codes and both programs' bridge partners (mkdssp's BP1/BP2 columns; foldmap's
DsspResult.bridges).

Writes under ../data/m2/dssp/:
  raw/<id>.dssp.gz   mkdssp's output as written
  <id>.json.gz       per compared residue [chain, number, icode, mkdssp code, foldmap code, mkdssp partners,
                     foldmap partners] (partners as "chain:number:icode"), plus counts and the samples it is in
  run.json           mkdssp version, foldmap location, entries, failures
Code: foldmap from build/tree/src (the harness commit; asserted); mkdssp from $MKDSSP.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import MANUSCRIPT, SOURCES, TREE, load_json, run_stamp  # noqa: E402

OUT = MANUSCRIPT / "data" / "m2" / "dssp"
REPORTS = {
    "tuned": "docs/validation/random100/report.json",
    "held": "docs/validation/random100-heldout/report.json.gz",
    "third": "docs/validation/random100-third/report.json.gz",
}


def parse(text: str):
    """mkdssp rows: (serial, chain, number, icode, code, bp1, bp2) per residue line ('!' breaks skipped), with the
    chain and code read as scripts/make_dssp_reference.parse_dssp reads them."""
    out, started = [], False
    for line in text.splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 17 or line[13] == "!":
            continue
        tail = line[137:].split()
        chain = tail[1] if len(tail) > 1 else line[11]
        code = line[16] if line[16] not in " S" else "-"
        if len(tail) >= 6:  # mmCIF input: full-width NUMBER, BP1, BP2 (the fixed columns overflow past 9999)
            serial, bp1, bp2 = int(tail[2]), int(tail[4]), int(tail[5])
        else:
            serial, bp1, bp2 = int(line[0:5]), int(line[25:29]), int(line[29:33])
        out.append((serial, chain, int(line[5:10]), line[10].strip(), code, bp1, bp2))
    return out


def af_accession(entry_id: str):
    return entry_id[5:-2] if entry_id.startswith("AF_AF") and entry_id.endswith("F1") else None


def one(job):
    pid, samples = job
    from foldmap.dssp import assign_dssp
    from foldmap.fetch import fetch
    from foldmap.io import load_backbone

    exe = os.environ["MKDSSP"]
    path = Path(fetch(af_accession(pid) or pid))
    with tempfile.TemporaryDirectory() as tmp:
        src = path
        if src.suffix == ".gz":
            src = Path(tmp) / path.name[:-3]
            src.write_bytes(gzip.decompress(path.read_bytes()))
        out = Path(tmp) / "out.dssp"
        done = subprocess.run([exe, str(src), str(out)], capture_output=True, text=True, timeout=1200)
        if done.returncode or not out.is_file():
            return {"id": pid, "error": (done.stderr.strip().splitlines() or ["failed"])[-1][:200]}
        text = out.read_text(encoding="utf-8", errors="replace")
    (OUT / "raw" / f"{pid}.dssp.gz").write_bytes(gzip.compress(text.encode("utf-8"), mtime=0))
    rows = parse(text)
    serial = {r[0]: f"{r[1]}:{r[2]}:{r[3]}" for r in rows}
    ref = {(r[1], r[2], r[3]): (r[4], sorted(serial[b] for b in (r[5], r[6]) if b and b in serial)) for r in rows}

    bb = load_backbone(path, "asu")
    res = assign_dssp(bb)
    keys = [(l.chain, l.seq, l.icode) for l in bb.labels]
    name = [f"{c}:{n}:{i}" for c, n, i in keys]
    partners = [set() for _ in keys]
    for i, j, _ in res.bridges:
        partners[i].add(name[j])
        partners[j].add(name[i])
    stood_in = {keys[k] for k in bb.linked}
    near = {keys[j] for k in bb.linked for j in range(k - 4, k + 5) if 0 <= j < len(keys)} - stood_in
    table = []
    for k, key in enumerate(keys):
        if key in ref and key not in near and key not in stood_in:
            mk, mkp = ref[key]
            table.append([*key, mk, res.ss[k], mkp, sorted(partners[k])])
    doc = {
        "id": pid,
        "samples": samples,
        "file": path.name,
        "residues_foldmap": len(keys),
        "residues_mkdssp": len(ref),
        "foldmap_only": len(set(keys) - stood_in - set(ref)),
        "mkdssp_only": len(set(ref) - set(keys)),
        "excluded_chromophore": len(stood_in) + len(near),
        "compared": len(table),
        "rows": table,
    }
    (OUT / f"{pid}.json.gz").write_bytes(gzip.compress(json.dumps(doc, separators=(",", ":")).encode(), mtime=0))
    return {"id": pid, "compared": len(table)}


def main() -> int:
    sys.path.insert(0, str(TREE / "src"))
    import foldmap

    where = Path(foldmap.__file__).resolve()
    assert str(where).startswith(str(TREE.resolve())), f"foldmap imported from {where}, not build/tree"
    version = subprocess.run([os.environ["MKDSSP"], "--version"], capture_output=True, text=True).stdout.strip()
    assert "mkdssp version 3.1.2" in version, version
    jobs: dict[str, list[str]] = {}
    for s, rel in REPORTS.items():
        rep = load_json(SOURCES / "HEAD" / rel)
        for e in rep["entries"]:
            if (e.get("dssp") or {}).get("agree") is not None:
                jobs.setdefault(e["id"], []).append(s)
    (OUT / "raw").mkdir(parents=True, exist_ok=True)
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(24) as ex:
        done = list(ex.map(one, sorted(jobs.items()), chunksize=1))
    bad = [d for d in done if "error" in d]
    run = {
        "generator": "scripts/paper/dssp_mkdssp.py",
        "source_stamp": run_stamp(__file__),
        "foldmap_from": str(where),
        "mkdssp": version,
        "entries": len(jobs),
        "written": len(done) - len(bad),
        "failures": bad,
        "samples": {s: sorted(i for i, ss in jobs.items() if s in ss) for s in REPORTS},
    }
    (OUT / "run.json").write_text(json.dumps(run, indent=1) + "\n", encoding="utf-8")
    print(f"{len(done) - len(bad)} of {len(jobs)} entries written; failures: {bad}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
