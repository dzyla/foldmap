"""foldmap diff's engine (plan 13e-1): what changed between two structures of one protein chain."""

import dataclasses
import os
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.diff import CONTACT_KEYS, Diff, _side, counts, diff_sides, diff_structures

DATA = Path(__file__).parent / "data"
UBQ, HHB = DATA / "1UBQ.cif", DATA / "2HHB.cif"


def _turn(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _write(st: gemmi.Structure, path: Path, turn=None) -> Path:
    """st written as PDB to path, every atom turned by `turn` (a 3x3 rotation) and moved, when given."""
    if turn is not None:
        for model in st:
            for chain in model:
                for res in chain:
                    for atom in res:
                        p = turn @ np.array(atom.pos.tolist()) + [12.5, -7.25, 30.0]
                        atom.pos = gemmi.Position(*p)
    st.write_pdb(str(path))
    return path


@pytest.fixture(scope="module")
def self_diff() -> Diff:
    return diff_structures(UBQ, "A", UBQ, "A", domains={"N": "1-40", "C": "res:A:41-76"})


def test_a_structure_against_itself_has_no_differences(self_diff):
    d = self_diff
    assert d.alignment.method == "sequence"
    assert d.only_a == () and d.only_b == ()
    assert d.ss_changes == () and d.pairing == ()
    assert set(d.contacts) == set(CONTACT_KEYS) and all(v == () for v in d.contacts.values())
    assert len(d.displacement) == 76 and set(d.displacement.values()) == {0.0}
    assert d.domain_turns == {"N": 0.0, "C": 0.0}
    assert all(v == 0 for v in counts(d).values())
    assert len(d.mates.mate) == len(_side(UBQ, "A", "A", "auto", True).sses)


def test_a_turned_copy_has_no_differences_either(self_diff, tmp_path):
    path = _write(gemmi.read_structure(str(UBQ)), tmp_path / "turned.pdb", _turn([1.0, 2.0, -0.5], 1.1))
    d = diff_structures(UBQ, "A", path, "A", domains={"N": "1-40", "C": "res:A:41-76"})
    assert counts(d) == counts(self_diff)
    assert max(d.displacement.values()) < 1e-2  # the file's coordinates are written to 0.001 Å
    assert max(d.domain_turns.values()) < 0.1
    assert d.core == self_diff.core


def _no_helix(tmp_path) -> tuple[Path, list]:
    """1UBQ less its alpha helix (the longest H element), written to tmp_path; that helix's residue numbers."""
    side = _side(UBQ, "A", "A", "auto", True)
    helix = max((s for s in side.sses if s.kind == "H"), key=len)
    gone = [side.bb.labels[k] for k in range(helix.start, helix.end + 1)]
    st = gemmi.read_structure(str(UBQ))
    chain = st[0]["A"]
    for k in reversed(range(len(chain))):
        if (chain[k].seqid.num, chain[k].seqid.icode.strip()) in {(l.seq, l.icode) for l in gone}:
            del chain[k]
    return _write(st, tmp_path / "ubq_no_helix.pdb"), [helix.id, *gone]


@pytest.mark.parametrize("match", ["auto", "structure"])
def test_a_deleted_helix_is_exactly_that_difference(tmp_path, match):
    path, (helix, *gone) = _no_helix(tmp_path)
    d = diff_structures(UBQ, "A", path, "A", match=match)
    assert d.alignment.method == ("sequence" if match == "auto" else "structure")
    assert d.only_a == (helix,) and d.only_b == ()
    assert d.ss_changes == tuple(("A", l.tag, "H", "-") for l in gone)
    assert d.pairing == () and all(v == () for v in d.contacts.values())
    assert len(d.displacement) == 76 - len(gone) and set(d.displacement.values()) == {0.0}


_PRINT = """
import sys
from foldmap.diff import diff_structures

d = diff_structures(sys.argv[1], "A", sys.argv[2], "A", domains={"N": "1-40"})
for f in ("only_a", "only_b", "ss_changes", "pairing", "core", "domain_turns", "contacts"):
    print(f, repr(getattr(d, f)))
print(repr(sorted(d.displacement.items())))
"""


def test_deterministic_across_hash_seeds(tmp_path):
    path, _ = _no_helix(tmp_path)
    turned = _write(gemmi.read_structure(str(path)), tmp_path / "turned_no_helix.pdb", _turn([0.3, -1.0, 2.0], 2.2))
    runs = set()
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _PRINT, str(UBQ), str(turned)],
            env=dict(os.environ, PYTHONHASHSEED=seed),
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1
    assert "only_a ('A:" in runs.pop()


def test_pairing_lost_and_register_moved():
    """B's bridges edited: one strand pairing's bridges dropped (lost), another's partner residues one further on
    (register +1, antiparallel: i + j grows by one)."""
    a = _side(UBQ, "A", "A", "auto", True)
    strands = [s for s in a.sses if s.kind == "E"]
    from foldmap.diff import _pairings

    pairs = _pairings(a)  # N to C
    assert len(pairs) >= 2
    (x1, y1), (x2, y2) = list(pairs)[:2]
    by_id = {s.id: s for s in strands}

    def between(br, x, y):  # bridges have i < j: i on the earlier strand x
        sx, sy = by_id[x], by_id[y]
        return sx.start <= br.i <= sx.end and sy.start <= br.j <= sy.end

    bridges = []
    for br in a.bridges:
        if between(br, x1, y1):
            continue  # pairing 1 lost
        if between(br, x2, y2):  # pairing 2: the later strand's residue one on (j - i or j + i: one more)
            br = br._replace(j=br.j + 1)
        bridges.append(br)
    d = diff_sides(a, dataclasses.replace(a, bridges=bridges))
    assert d.pairing == ((x1, y1, "lost"), (x2, y2, "register +1"))
    assert d.only_a == d.only_b == () and d.ss_changes == ()


def test_domains_are_checked(self_diff):
    with pytest.raises(ValueError, match="domain X: 'B:1-40' is not a residue range of A's chain A"):
        diff_structures(UBQ, "A", UBQ, "A", domains={"X": "B:1-40"})
    with pytest.raises(ValueError, match="domain X: '500-600' names no residue of chain A"):
        diff_structures(UBQ, "A", UBQ, "A", domains={"X": "500-600"})
    with pytest.raises(ValueError, match="no chain 'Q' in B"):
        diff_structures(UBQ, "A", UBQ, "Q")


def test_interface_and_ligand_changes(tmp_path):
    """2HHB's chain A against a copy without chain B and without its haems: the interface with B lost, the haem
    lost; the interfaces with C and D, and chain A itself, unchanged."""
    st = gemmi.read_structure(str(HHB))
    st.setup_entities()
    st[0].remove_chain("B")
    for chain in st[0]:
        for k in reversed(range(len(chain))):
            if chain[k].name == "HEM":
                del chain[k]
    path = _write(st, tmp_path / "hhb_less.pdb")
    same = diff_structures(HHB, "A", HHB, "A", assembly="asu")
    assert all(v == () for v in same.contacts.values())
    d = diff_structures(HHB, "A", path, "A", assembly="asu")
    assert d.only_a == d.only_b == () and d.ss_changes == () and d.pairing == ()
    assert d.contacts["interfaces"] == (("B", "ca", "lost"),)
    assert d.contacts["ligands"] == (("HEM", "lost"),)
    assert d.contacts["elements"] == d.contacts["interface_elements"] == d.contacts["ligand_elements"] == ()


def test_two_proteins_take_the_structure_route():
    """1UBQ against T4 lysozyme: not one protein by A19's rule, so mapped by structure; every element of each that
    has no mate is named, and nothing fails on a small core."""
    d = diff_structures(UBQ, "A", DATA / "2LZM.cif", "A")
    assert d.alignment.method == "structure"
    assert len(d.only_a) + len(d.mates.mate) + len(d.mates.pieces) == len(_side(UBQ, "A", "A", "auto", True).sses)
    assert len(d.core) >= 3 and len(d.displacement) == len(d.alignment.pairs)


def test_a_turned_domain_shows_its_turn(tmp_path):
    """1UBQ with residues 41-76 turned 60 degrees about residue 40's CA (a hinge): the core fit, started from the
    TM-score's superposition, keeps residues 1-40 (the turned ones move past 2·d0 and are dropped), so N turns by
    nothing and C by 60 degrees."""
    st = gemmi.read_structure(str(UBQ))
    moved = [a for r in st[0]["A"] if 41 <= r.seqid.num <= 76 for a in r]
    centre = np.array(st[0]["A"]["40"][0]["CA"][0].pos.tolist())
    turn = _turn([0.0, 1.0, 1.0], np.radians(60.0))
    for a in moved:
        a.pos = gemmi.Position(*(turn @ (np.array(a.pos.tolist()) - centre) + centre))
    d = diff_structures(UBQ, "A", _write(st, tmp_path / "hinge.pdb"), "A", domains={"N": "1-40", "C": "41-76"})
    assert d.domain_turns["N"] < 0.1
    assert d.domain_turns["C"] == pytest.approx(60.0, abs=0.1)
    assert all(i < 40 for i, _ in d.core)
