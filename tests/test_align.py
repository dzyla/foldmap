"""The chain aligner (plan 13a): the sequence, structure and file routes, TM-score, the superposition."""

import os
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.align import SAME_FOLD_TM, ChainAlignment, Mates, align_chains, d0, element_mates, tm_score
from foldmap.io import load_backbone
from foldmap.model import Backbone, ResLabel
from foldmap.sequence import _one_letter

DATA = Path(__file__).parent / "data"


def _turn(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


TURN = _turn([1.0, 2.0, -0.5], 1.1)  # a fixed proper rotation
SHIFT = np.array([12.5, -7.25, 30.0])


def _moved(bb: Backbone, rot=TURN, shift=SHIFT, keep=None) -> Backbone:
    """bb's chain A turned by rot and moved by shift; keep: the residue indices kept (all by default)."""
    keep = list(range(len(bb))) if keep is None else list(keep)
    xyz = bb.xyz[keep] @ rot.T + shift
    return Backbone([bb.labels[k] for k in keep], xyz, bb.prev[keep].copy(), set(bb.asu_chains))


@pytest.fixture(scope="module")
def ubq_a():
    return load_backbone(DATA / "1UBQ.cif").only(["A"])


def test_d0():
    assert d0(10) == 0.5 and d0(15) == 0.5 and d0(19) == 0.5
    assert d0(76) == pytest.approx(1.24 * 61 ** (1 / 3) - 1.8)


def test_self_alignment_is_identity(ubq_a):
    b = _moved(ubq_a)
    out = align_chains(ubq_a, "A", b, "A")
    assert isinstance(out, ChainAlignment)
    assert out.pairs == tuple((i, i) for i in range(len(ubq_a)))
    assert out.tm_score == pytest.approx(1.0, abs=1e-9)
    assert out.rmsd < 1e-6
    assert out.identity == 1.0
    assert out.method == "sequence"


def test_rotation_recovered(ubq_a):
    b = _moved(ubq_a)
    out = align_chains(ubq_a, "A", b, "A")
    assert np.linalg.det(out.rotation) == pytest.approx(1.0, abs=1e-9)
    assert np.abs(b.ca @ out.rotation.T + out.translation - ubq_a.ca).max() < 1e-6
    # the same transform tm_score finds directly
    tm, r, t, rmsd = tm_score(ubq_a.ca, b.ca, out.pairs, len(ubq_a))
    assert tm == out.tm_score and np.array_equal(r, out.rotation) and np.array_equal(t, out.translation)


def _fasta(path, rows):
    path.write_text("".join(f">{name} some description\n{row}\n" for name, row in rows))
    return path


def test_file_route_reads_two_rows(ubq_a, tmp_path):
    seq = _one_letter([l.name for l in ubq_a.labels])
    gone = 30
    b = _moved(ubq_a, keep=[k for k in range(len(ubq_a)) if k != gone])
    row_b = seq[:gone] + "-" + seq[gone + 1 :]
    out = align_chains(ubq_a, "A", b, "A", str(_fasta(tmp_path / "two.fasta", [("ref", seq), ("other", row_b)])))
    assert out.method == "file"
    assert out.pairs == tuple((i, i) for i in range(gone)) + tuple((i, i - 1) for i in range(gone + 1, len(seq)))
    assert out.identity == 1.0 and out.tm_score == pytest.approx((len(seq) - 1) / len(seq), abs=1e-9)

    # the file is obeyed even where the sequence route would pair differently: B shifted one residue along A
    same = _moved(ubq_a)
    shifted = _fasta(tmp_path / "shifted.fasta", [("ref", "-" + seq), ("other", seq + "-")])
    out = align_chains(ubq_a, "A", same, "A", str(shifted))
    assert out.method == "file" and out.pairs == tuple((i, i + 1) for i in range(len(seq) - 1))
    assert out.identity < 0.2

    three = _fasta(tmp_path / "three.fasta", [("ref", seq), ("other", seq), ("extra", seq)])
    with pytest.raises(ValueError, match="extra"):
        align_chains(ubq_a, "A", same, "A", str(three))

    wrong = seq[:10] + ("A" if seq[10] != "A" else "G") + seq[11:]
    bad = _fasta(tmp_path / "bad.fasta", [("ref", seq), ("mutant", wrong)])
    with pytest.raises(ValueError, match="mutant"):
        align_chains(ubq_a, "A", same, "A", str(bad))
    bad = _fasta(tmp_path / "bad_ref.fasta", [("first", wrong), ("other", seq)])
    with pytest.raises(ValueError, match="first"):
        align_chains(ubq_a, "A", same, "A", str(bad))


def test_unknown_chain_or_match(ubq_a, tmp_path):
    with pytest.raises(ValueError, match="chain 'Z'"):
        align_chains(ubq_a, "A", _moved(ubq_a), "Z")
    with pytest.raises(ValueError, match="nothere"):
        align_chains(ubq_a, "A", _moved(ubq_a), "A", str(tmp_path / "nothere.fasta"))


def test_never_a_mirror(ubq_a):
    mirror = Backbone(list(ubq_a.labels), ubq_a.xyz * np.array([-1.0, 1.0, 1.0]), ubq_a.prev.copy())
    out = align_chains(ubq_a, "A", _moved(mirror), "A", "sequence")
    assert out.pairs == tuple((i, i) for i in range(len(ubq_a)))
    assert np.linalg.det(out.rotation) == pytest.approx(1.0, abs=1e-9)
    assert out.tm_score < 0.5
    auto = align_chains(ubq_a, "A", _moved(mirror), "A")  # its sequence pairs do not superpose: the structure route
    assert auto.method == "structure" and auto.tm_score < 0.5
    assert np.linalg.det(auto.rotation) == pytest.approx(1.0, abs=1e-9)


_PRINT = """
import sys
import numpy as np
from foldmap.align import align_chains
from foldmap.io import load_backbone
from foldmap.model import Backbone, ResLabel

a = load_backbone(sys.argv[1]).only(["A"])
keep = [k for k in range(len(a)) if k not in (5, 6, 40, 70)]
c, s = np.cos(0.7), np.sin(0.7)
rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
b = Backbone([a.labels[k] for k in keep], a.xyz[keep] @ rot.T + [3.0, -4.0, 5.0], a.prev[keep])
out = align_chains(a, "A", b, "A")
print(repr(out.pairs))
print(repr(out.tm_score))
"""


def test_deterministic_across_hash_seeds():
    runs = set()
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _PRINT, str(DATA / "1UBQ.cif")],
            env=dict(os.environ, PYTHONHASHSEED=seed),
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1
    pairs, tm = runs.pop().splitlines()
    assert pairs.startswith("((0, 0), ") and float(tm) > 0.9


def _run_twice(script, *args):
    """stdout of `script` under PYTHONHASHSEED 0 and 1, which must agree."""
    runs = set()
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", script, *map(str, args)],
            env=dict(os.environ, PYTHONHASHSEED=seed),
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1
    return runs.pop()


# the structure route (13a-2)

_NAMES = sorted("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def _renamed(bb: Backbone) -> Backbone:
    """bb with every residue renamed to the next standard amino acid (alphabetically): no identity left."""
    swap = {x: _NAMES[(k + 1) % len(_NAMES)] for k, x in enumerate(_NAMES)}
    return Backbone([l._replace(name=swap.get(l.name, "ALA")) for l in bb.labels], bb.xyz.copy(), bb.prev.copy())


def test_helix_deleted_keeps_the_rest(ubq_a, tmp_path):
    st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    chain = st[0]["A"]
    for k in reversed(range(len(chain))):
        if 23 <= chain[k].seqid.num <= 34 and chain[k].name != "HOH":
            del chain[k]
    path = tmp_path / "ubq_no_helix.pdb"
    st.write_pdb(str(path))
    b = _moved(load_backbone(path).only(["A"]))
    assert len(b) == len(ubq_a) - 12
    out = align_chains(ubq_a, "A", b, "A", "structure")
    assert out.method == "structure"
    kept = [i for i in range(len(ubq_a)) if not 22 <= i <= 33]
    assert out.pairs == tuple((i, j) for j, i in enumerate(kept))
    assert out.tm_score > 0.8
    assert np.abs(b.ca @ out.rotation.T + out.translation - ubq_a.ca[kept]).max() < 1e-6


def test_auto_switches_on_identity(ubq_a):
    other = _moved(_renamed(ubq_a))
    out = align_chains(ubq_a, "A", other, "A")
    assert out.identity < 0.4 and out.method == "structure"
    assert out.pairs == tuple((i, i) for i in range(len(ubq_a)))  # the fold found without the sequence
    assert out.tm_score == pytest.approx(1.0, abs=1e-9)
    assert align_chains(ubq_a, "A", _moved(ubq_a), "A").method == "sequence"
    assert align_chains(ubq_a, "A", other, "A", "sequence").method == "sequence"  # asked for: obeyed


def test_auto_takes_the_sequence_only_when_its_pairs_superpose(ubq_a):
    # the same residue names on a chain run backwards in space: identical by sequence, not one fold by it
    back = Backbone(list(ubq_a.labels), ubq_a.xyz[::-1].copy() @ TURN.T + SHIFT, ubq_a.prev.copy())
    seq = align_chains(ubq_a, "A", back, "A", "sequence")
    assert seq.identity == 1.0 and seq.tm_score < SAME_FOLD_TM
    out = align_chains(ubq_a, "A", back, "A")
    assert out.method == "structure" and out.tm_score >= seq.tm_score


def test_permutation_shows_as_low_coverage(ubq_a):
    order = list(range(39, 76)) + list(range(39))  # residues 40-76, then 1-39, as one chain
    prev = ubq_a.prev[order].copy()
    prev[0] = prev[37] = False
    permuted = Backbone([ubq_a.labels[k] for k in order], ubq_a.xyz[order] @ TURN.T + SHIFT, prev)
    out = align_chains(ubq_a, "A", permuted, "A", "structure")
    assert out.method == "structure"
    assert len(out.pairs) < 0.7 * 76
    assert all(order[j] == i for i, j in out.pairs)  # what is paired is the right residue


def test_never_a_mirror_on_the_structure_route(ubq_a):
    mirror = Backbone(list(ubq_a.labels), ubq_a.xyz * np.array([-1.0, 1.0, 1.0]), ubq_a.prev.copy())
    out = align_chains(ubq_a, "A", _moved(mirror), "A", "structure")
    assert np.linalg.det(out.rotation) == pytest.approx(1.0, abs=1e-9)
    assert out.tm_score < 0.5


_PRINT_STRUCTURE = """
import sys
import numpy as np
from foldmap.align import align_chains
from foldmap.io import load_backbone
from foldmap.model import Backbone

a = load_backbone(sys.argv[1]).only(["A"])
keep = [k for k in range(len(a)) if k not in (5, 6, 40, 70)]
c, s = np.cos(0.7), np.sin(0.7)
rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
wobble = 0.4 * np.sin(np.arange(len(keep)))[:, None, None]
b = Backbone([a.labels[k] for k in keep], a.xyz[keep] @ rot.T + [3.0, -4.0, 5.0] + wobble, a.prev[keep])
out = align_chains(a, "A", b, "A", "structure")
print(repr(out.pairs))
print(repr(out.tm_score))
"""


def test_structure_route_deterministic():
    pairs, tm = _run_twice(_PRINT_STRUCTURE, DATA / "1UBQ.cif").splitlines()
    keep = [k for k in range(76) if k not in (5, 6, 40, 70)]
    assert pairs == repr(tuple((i, j) for j, i in enumerate(keep)))
    assert float(tm) > 0.8


def _place(a, b, c, bond, angle, torsion):
    """The point d with |cd| = bond, angle bcd and torsion abcd (degrees): NeRF."""
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    ang, tor = np.radians(angle), np.radians(torsion)
    return c - bond * np.cos(ang) * bc + bond * np.sin(ang) * (np.cos(tor) * m + np.sin(tor) * n)


def _helix(count, phi=-57.0, psi=-47.0):
    """An ideal alpha helix's backbone (count, 4, 3): N, CA, C, O."""
    n, ca = np.zeros(3), np.array([1.458, 0.0, 0.0])
    c = _place(np.array([0.0, 1.0, 0.0]), n, ca, 1.525, 111.2, -60.0)
    out = []
    for k in range(count):
        if k:
            n = _place(out[-1][0], out[-1][1], out[-1][2], 1.329, 116.2, psi)
            ca = _place(out[-1][1], out[-1][2], n, 1.458, 121.7, 180.0)
            c = _place(out[-1][2], n, ca, 1.525, 111.2, phi)
        out.append([n, ca, c, _place(n, ca, c, 1.231, 120.5, psi + 180.0)])
    return np.array(out)


def _bundle(helices=18, size=24, loop=4):
    """A synthetic helix bundle of helices x size residues joined by straight loops of `loop` residues (500 with
    the defaults): antiparallel neighbours on a 10 Å grid, 6 a row, a snake through the rows."""
    h = _helix(size)
    centre = h[:, 1].mean(axis=0)
    _, _, vt = np.linalg.svd(h[:, 1] - centre)
    if np.linalg.det(vt) < 0:
        vt[2] *= -1  # a proper rotation: the helix stays right-handed
    base = (h - centre) @ vt.T  # the helix axis along x
    turn = np.diag([-1.0, -1.0, 1.0])  # a half turn about z: the next helix runs back
    parts = []
    for k in range(helices):
        row, col = divmod(k, 6)
        col = col if row % 2 == 0 else 5 - col
        piece = (base @ turn.T if k % 2 else base) + [0.0, 10.0 * col, 10.0 * row]
        if parts:
            last, first = parts[-1][-1], piece[0]
            parts.append(np.array([last + (first - last) * s / (loop + 1) for s in range(1, loop + 1)]))
        parts.append(piece)
    xyz = np.concatenate(parts)
    labels = [ResLabel("A", k + 1, "", "ALA") for k in range(len(xyz))]
    prev = np.ones(len(xyz), bool)
    prev[0] = False
    return Backbone(labels, xyz, prev, {"A"})


@pytest.mark.slow
def test_budget():
    a = _bundle()
    assert len(a) == 500
    keep = [k for k in range(len(a)) if k not in (100, 101, 300)]
    k = np.arange(len(a))[:, None]
    wobble = 0.6 * np.hstack([np.sin(0.37 * k), np.cos(0.53 * k), np.sin(0.71 * k + 1.0)])
    b = _moved(Backbone(a.labels, a.xyz + wobble[:, None, :], a.prev.copy(), {"A"}), keep=keep)
    t0 = time.perf_counter()
    out = align_chains(a, "A", b, "A", "structure")
    took = time.perf_counter() - t0
    assert out.method == "structure"
    right = sum(keep[j] == i for i, j in out.pairs)
    assert right >= 0.95 * len(keep) and out.tm_score > 0.8
    assert took < 3.0, f"{took:.2f} s"


# element mates (13a-4)


def _elements(bb: Backbone):
    from foldmap.dssp import assign_dssp
    from foldmap.ss import build_sses

    d = assign_dssp(bb)
    return build_sses(bb, d.ss, short_helices=True, bridges=d.bridges)


def _no_helix(tmp_path) -> Backbone:
    """1UBQ's chain A without its helix (residues 23-34), turned and moved."""
    st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    chain = st[0]["A"]
    for k in reversed(range(len(chain))):
        if 23 <= chain[k].seqid.num <= 34 and chain[k].name != "HOH":
            del chain[k]
    path = tmp_path / "ubq_no_helix.pdb"
    st.write_pdb(str(path))
    return _moved(load_backbone(path).only(["A"]))


def _swap(sses, old: str, new):
    """The element list with element `old` replaced by the elements `new`."""
    at = [s.id for s in sses].index(old)
    return [*sses[:at], *new, *sses[at + 1 :]]


def test_self_every_element_its_own_mate(ubq_a):
    sses = _elements(ubq_a)
    out = element_mates(sses, sses, align_chains(ubq_a, "A", _moved(ubq_a), "A"))
    assert isinstance(out, Mates)
    assert out.mate == {s.id: s.id for s in sses}
    assert out.pieces == {} and out.extras == () and out.only_a == ()
    assert out.same_fold


def test_helix_deleted_leaves_it_only_in_a(ubq_a, tmp_path):
    b = _no_helix(tmp_path)
    sses_a, sses_b = _elements(ubq_a), _elements(b)
    out = element_mates(sses_a, sses_b, align_chains(ubq_a, "A", b, "A"))
    assert out.only_a == ("A:23-34",)
    assert out.mate == {s.id: s.id for s in sses_a if s.id != "A:23-34"}
    assert out.extras == () and out.pieces == {}
    assert out.same_fold


def test_split_helix_is_one_mate_and_one_piece(ubq_a):
    from foldmap.ss import _sse

    sses = _elements(ubq_a)
    split = _swap(sses, "A:23-34", [_sse("H", ubq_a, 22, 28), _sse("H", ubq_a, 29, 33)])  # 23-29, 30-34
    out = element_mates(sses, split, align_chains(ubq_a, "A", _moved(ubq_a), "A"))
    assert out.mate["A:23-29"] == "A:23-34" and "A:30-34" not in out.mate
    assert out.pieces == {"A:30-34": "A:23-34"}
    assert out.extras == () and out.only_a == ()


def test_strand_never_mates_a_helix(ubq_a):
    from foldmap.ss import _sse

    sses = _elements(ubq_a)
    aln = align_chains(ubq_a, "A", _moved(ubq_a), "A")
    other = _swap(sses, "A:23-34", [_sse("E", ubq_a, 22, 33)])  # a strand where A has its helix
    out = element_mates(sses, other, aln)
    assert "A:23-34" not in out.mate and out.extras == ("A:23-34",) and out.only_a == ("A:23-34",)
    assert out.pieces == {}  # a strand is no piece of a helix
    # a 3-10 helix (G) and an alpha helix (H) may mate
    other = _swap(sses, "A:57-59", [_sse("H", ubq_a, 56, 58)])
    assert element_mates(sses, other, aln).mate["A:57-59"] == "A:57-59"


def test_same_fold_needs_the_tm_score_and_the_elements(ubq_a):
    sses = _elements(ubq_a)
    aln = align_chains(ubq_a, "A", _moved(ubq_a), "A")
    low = replace(aln, tm_score=0.3)  # a negative pair's TM-score
    out = element_mates(sses, sses, low)
    assert len(out.mate) == len(sses) and out.same_fold is False
    assert element_mates(sses, sses, replace(aln, tm_score=SAME_FOLD_TM)).same_fold
    few = replace(aln, pairs=aln.pairs[:20])  # only the first two strands' residues paired
    out = element_mates(sses, sses, few)
    assert set(out.mate) == {"A:2-7", "A:12-16"} and out.same_fold is False
    assert len(out.extras) == len(out.only_a) == len(sses) - 2


def test_elements_of_other_chains_are_left_out(ubq_a):
    """Element indices count the whole structure's residues; the alignment's count its chain's: chain B of a
    two-chain B (after chain A's 76 residues) still finds its mates, and chain A's elements of B are left out."""
    far = _moved(ubq_a, rot=np.eye(3), shift=np.array([100.0, 0.0, 0.0]))
    two = Backbone(
        [*ubq_a.labels, *(l._replace(chain="B") for l in far.labels)],
        np.concatenate([ubq_a.xyz, far.xyz]),
        np.concatenate([ubq_a.prev, far.prev]),
        {"A", "B"},
    )
    sses_a, sses_two = _elements(ubq_a), _elements(two)
    assert min(s.start for s in sses_two if s.chain == "B") >= len(ubq_a)
    out = element_mates(sses_a, sses_two, align_chains(ubq_a, "A", two, "B"))
    assert out.mate == {"B" + s.id[1:]: s.id for s in sses_a}
    assert out.extras == () and out.only_a == () and out.same_fold
