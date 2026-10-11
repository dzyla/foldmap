"""foldmap's DSSP against mkdssp 3.1.2, residue by residue (reference files from scripts/make_dssp_reference.py)."""

import re
from pathlib import Path

import pytest

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone

DATA = Path(__file__).parent / "data"
IDS = sorted(p.stem for p in (DATA / "dssp").glob("*.ss"))


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
    # a chromophore foldmap keeps in the chain (CRO66) mkdssp skips, so it sees the chain broken there: compare
    # neither it nor the residues whose turns and bridges reach across it (within four)
    stood_in = {keys[k] for k in bb.linked}
    near = {keys[j] for k in bb.linked for j in range(k - 4, k + 5) if 0 <= j < len(keys)} - stood_in
    assert set(keys) - stood_in == set(ref), "foldmap and mkdssp read different residues"
    return len(keys), [
        f"{c}{n}{i} mkdssp {ref[(c, n, i)]} foldmap {ss[k]}"
        for k, (c, n, i) in enumerate(keys)
        if (c, n, i) in ref and (c, n, i) not in near and ss[k] != ref[(c, n, i)]
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
    assert not bad, f"{len(bad)} of {n} residues differ from mkdssp: " + "; ".join(bad[:10])
