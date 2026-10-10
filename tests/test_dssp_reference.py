"""foldmap's DSSP against mkdssp 3.1.2, residue by residue (reference files from scripts/make_dssp_reference.py)."""

import re
from pathlib import Path

import pytest

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone

DATA = Path(__file__).parent / "data"
IDS = sorted(p.stem for p in (DATA / "dssp").glob("*.ss"))

# Floor on per-residue agreement (%) while the DSSP rules are brought in line with mkdssp; raised task by task.
MIN_Q8 = {
    "1BL8": 100.0, "1C3W": 96.3, "1EMA": 100.0, "1LMB": 97.7, "1TIM": 100.0, "1UBQ": 100.0, "1ZAA": 100.0,
    "25IN": 96.9, "2HHB": 95.4, "2JTY": 92.3, "2LZM": 95.1, "4CH9": 99.4, "5NKT": 82.1, "6R74": 84.0,
    "6Y7S": 94.9, "6ZS5": 98.0, "6ZYA": 98.0, "8UTF": 96.8, "8UUP": 96.1, "9TX7": 98.5, "AF-P04637": 98.7,
}  # fmt: skip


def reference(path: Path) -> dict[tuple[str, int, str], str]:
    """(auth chain, number, insertion code) -> DSSP code, from one reference file."""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        chain, first, codes = line.split()
        num, icode = re.fullmatch(r"(-?\d+)([A-Za-z]?)", first).groups()
        for k, code in enumerate(codes):
            out[(chain, int(num) + k, icode)] = code
    return out


def disagreements(name: str) -> tuple[int, list[str]]:
    ref = reference(DATA / "dssp" / f"{name}.ss")
    bb = load_backbone(DATA / f"{name}.cif", "asu")
    ss = assign_dssp(bb).ss
    keys = [(l.chain, l.seq, l.icode) for l in bb.labels]
    assert set(keys) == set(ref), "foldmap and mkdssp read different residues"
    return len(keys), [
        f"{c}{n}{i} mkdssp {ref[(c, n, i)]} foldmap {ss[k]}"
        for k, (c, n, i) in enumerate(keys)
        if ss[k] != ref[(c, n, i)]
    ]


def test_every_structure_has_a_reference():
    assert IDS == sorted(p.stem for p in DATA.glob("*.cif"))


def test_reference_format(tmp_path):
    path = tmp_path / "x.ss"
    path.write_text("A -2 -EE\nA 100A H\nB 7 T\n", encoding="utf-8")
    assert reference(path) == {
        ("A", -2, ""): "-",
        ("A", -1, ""): "E",
        ("A", 0, ""): "E",
        ("A", 100, "A"): "H",
        ("B", 7, ""): "T",
    }


@pytest.mark.parametrize("name", IDS)
def test_agrees_with_mkdssp(name):
    n, bad = disagreements(name)
    q8 = 100 * (n - len(bad)) / n
    assert q8 >= MIN_Q8[name], f"Q8 {q8:.2f}% < {MIN_Q8[name]}%: " + "; ".join(bad[:10])
