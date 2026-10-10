"""Write mkdssp's per-residue secondary structure for every test structure to tests/data/dssp/<ID>.ss.

Needs mkdssp (3.1.2 was used) on PATH; CI does not run this, it only reads the files it writes.
Each line covers residues with consecutive author numbers and no insertion code:
    <auth chain> <first number> <codes>
A residue with an insertion code gets a line of its own, its number followed by the code letter (100A).
Codes are DSSP's, with bend (S) and blank written as '-'.

    python scripts/make_dssp_reference.py [tests/data/1UBQ.cif ...]
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
OUT = DATA / "dssp"


def parse_dssp(text: str) -> list[tuple[str, int, str, str]]:
    """(auth chain, number, insertion code, code) per residue, in file order."""
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
        out.append((chain, int(line[5:10]), line[10].strip(), code))
    return out


def encode(rows: list[tuple[str, int, str, str]]) -> str:
    lines: list[list] = []
    for chain, num, icode, code in rows:
        last = lines[-1] if lines else None
        if icode:
            lines.append([chain, f"{num}{icode}", code])
        elif last and last[0] == chain and isinstance(last[1], int) and last[1] + len(last[2]) == num:
            last[2] += code
        else:
            lines.append([chain, num, code])
    return "".join(f"{c} {n} {s}\n" for c, n, s in lines)


def main(paths: list[str]) -> None:
    OUT.mkdir(exist_ok=True)
    for path in map(Path, paths or sorted(DATA.glob("*.cif"))):
        with tempfile.TemporaryDirectory() as tmp:
            dssp = Path(tmp) / "out.dssp"
            subprocess.run(["mkdssp", str(path), str(dssp)], check=True)
            text = dssp.read_text(encoding="utf-8")
        (OUT / f"{path.stem}.ss").write_text(encode(parse_dssp(text)), encoding="utf-8")
        print(path.stem)


if __name__ == "__main__":
    main(sys.argv[1:])
