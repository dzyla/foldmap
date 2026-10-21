"""The chain aligner (plan 13a): the sequence and file routes, TM-score, the superposition."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from foldmap.align import ChainAlignment, align_chains, d0, tm_score
from foldmap.io import load_backbone
from foldmap.model import Backbone
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


def test_structure_route_not_yet(ubq_a):
    with pytest.raises(NotImplementedError):
        align_chains(ubq_a, "A", _moved(ubq_a), "A", "structure")


def test_unknown_chain_or_match(ubq_a, tmp_path):
    with pytest.raises(ValueError, match="chain 'Z'"):
        align_chains(ubq_a, "A", _moved(ubq_a), "Z")
    with pytest.raises(ValueError, match="nothere"):
        align_chains(ubq_a, "A", _moved(ubq_a), "A", str(tmp_path / "nothere.fasta"))


def test_never_a_mirror(ubq_a):
    mirror = Backbone(list(ubq_a.labels), ubq_a.xyz * np.array([-1.0, 1.0, 1.0]), ubq_a.prev.copy())
    out = align_chains(ubq_a, "A", _moved(mirror), "A")
    assert out.pairs == tuple((i, i) for i in range(len(ubq_a)))
    assert np.linalg.det(out.rotation) == pytest.approx(1.0, abs=1e-9)
    assert out.tm_score < 0.5


_PRINT = """
import sys
import numpy as np
from foldmap.align import align_chains
from foldmap.io import load_backbone
from foldmap.model import Backbone

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
