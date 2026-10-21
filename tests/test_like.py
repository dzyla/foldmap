"""A layout drawn like a reference (plan 13a-5): frame_like, Like and the homes build_layout takes from it."""

import functools
import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from foldmap import cli, figcheck
from foldmap.align import Mates, align_chains, element_mates
from foldmap.dssp import assign_dssp
from foldmap.frame import CylinderFrame, Frame, frame_like
from foldmap.layout import _K_HOME, _K_LIKE
from foldmap.like import Like, like_from
from foldmap.sheets import build_sheets
from foldmap.style import Style

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).resolve().parent.parent


def _bench():
    spec = importlib.util.spec_from_file_location("foldmap_bench", ROOT / "scripts" / "bench.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BENCH = _bench()
signature = BENCH.signature
OPTIONS = BENCH.OPTIONS  # every make_layout option the bench figures depend on


def _turn(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def laid_out_like(path, like: Like | None, frame, **kw):
    """cli.make_layout of `path` with build_layout given `like` and the page frame `frame` (13a-6 gives make_layout
    these itself)."""
    with (
        mock.patch.object(cli, "build_layout", functools.partial(cli.build_layout, like=like)),
        mock.patch.object(cli, "view_frame", lambda *a, **k: frame),
    ):
        return cli.make_layout(path, **{**OPTIONS, **kw})


def drawn_like(path_a, chain_a, path_b, chain_b, look=None, mates=None):
    """B (one chain) drawn like A (one chain): A laid out, the chains aligned, the elements mated (or `mates`), the
    Like built from A's finished layout, B laid out through A's frame carried over. Returns (layout B, sses B, bb B,
    layout A, Like)."""
    look = look or Style()
    lay_a, sses_a, bb_a = cli.make_layout(path_a, look=look, chains=[chain_a], **OPTIONS)
    _, sses_b, bb_b = cli.make_layout(path_b, look=look, chains=[chain_b], **OPTIONS)
    aln = align_chains(bb_a, chain_a, bb_b, chain_b)
    mates = mates or element_mates(sses_a, sses_b, aln)
    like = like_from(lay_a, sses_a, sses_b, mates, look, build_sheets(sses_b, assign_dssp(bb_b).bridges))
    frame = frame_like(lay_a.frame, aln.rotation, aln.translation)
    lay_b, sses_b, bb_b = laid_out_like(path_b, like, frame, look=look, chains=[chain_b])
    assert lay_b.frame is frame
    return lay_b, sses_b, bb_b, lay_a, like


def test_frame_like_carries_known_rotation():
    rot_a = _turn([0.3, -1.0, 0.6], 0.8)  # A's page axes: the rows of a proper rotation
    frame_a = Frame(np.array([4.0, -2.5, 11.0]), rot_a[0], rot_a[1], rot_a[2])
    r, t = _turn([1.0, 2.0, -0.5], 1.1), np.array([12.5, -7.25, 30.0])
    b = np.array([[np.sin(3.0 * k), np.cos(2.0 * k) * 7.0, 0.5 * k - 9.0] for k in range(40)])  # B's points
    a = b @ r.T + t  # their places in A: A ≈ R·B + t
    frame_b = frame_like(frame_a, r, t)
    assert np.abs(frame_b.project(b) - frame_a.project(a)).max() < 1e-9
    assert np.abs(frame_b.depth(b) - frame_a.depth(a)).max() < 1e-9
    assert np.linalg.det(np.array([frame_b.u, frame_b.v, frame_b.w])) == pytest.approx(1.0)  # never a mirror
    mirrored = Frame(frame_a.origin, frame_a.u, frame_a.v, frame_a.w, sector={"A": 0, "B": 1}, mirror_x=0.0)
    unrolled = CylinderFrame(np.zeros(3), np.array([0, 0, 1.0]), np.array([1.0, 0, 0]), 3, 10.0, {"A": 0})
    for other in (mirrored, unrolled):
        with pytest.raises(TypeError, match="plain Frame"):
            frame_like(other, r, t)


@pytest.mark.parametrize("name", ["1UBQ", "2LZM"])
def test_like_itself_reproduces_the_layout(name):
    """The invariant that tests the design (spec 4.3): a structure drawn like itself is its own layout, exactly."""
    path = DATA / f"{name}.cif"
    plain, sses, _ = cli.make_layout(path, chains=["A"], **OPTIONS)
    lay, _, _, _, like = drawn_like(path, "A", path, "A")
    assert set(like.target) == {s.id for s in sses}  # every element held at its own place
    assert like.weight == _K_LIKE / _K_HOME == 8
    assert signature(lay) == signature(plain)
    # the places are the reference's, not the projection's: the reference moved, the layout moves with it
    moved = replace(
        like,
        target={k: (x + 3.0, y - 2.0) for k, (x, y) in like.target.items()},
        unit={k: (x + 3.0, y - 2.0, ox, oy) for k, (x, y, ox, oy) in like.unit.items()},
    )
    shifted, _, _ = laid_out_like(path, moved, lay.frame, chains=["A"])
    for k, p in plain.placed.items():
        assert (shifted.placed[k].cx, shifted.placed[k].cy) == pytest.approx((p.cx + 3.0, p.cy - 2.0), abs=1e-9)


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "1TIM"])
def test_like_none_changes_nothing(name):
    path = DATA / f"{name}.cif"
    plain, _, _ = cli.make_layout(path, **OPTIONS)
    with mock.patch.object(cli, "build_layout", functools.partial(cli.build_layout, like=None)):
        given, _, _ = cli.make_layout(path, **OPTIONS)
    assert signature(given) == signature(plain)


def _strand_order(lay, ids):
    """The strands of a sheet block left to right as drawn."""
    return [k for k in sorted(ids, key=lambda k: lay.placed[k].cx)]


def test_sheet_order_kept():
    """A B whose strands' mates stand in another order in A (here 1UBQ's strands mated in reverse) keeps the strand
    order of its own pairing: Like moves the sheet as a whole, never a strand out of its place in it."""
    path = DATA / "1UBQ.cif"
    alone, sses, bb = cli.make_layout(path, chains=["A"], **OPTIONS)
    strands = [s.id for s in sses if s.kind == "E"]
    helices = [s.id for s in sses if s.kind != "E"]
    mate = {**dict(zip(strands, reversed(strands))), **{h: h for h in helices}}
    mates = Mates(mate=mate, pieces={}, extras=(), only_a=(), same_fold=True)
    lay, sses_b, bb_b, lay_a, like = drawn_like(path, "A", path, "A", mates=mates)
    assert like.target[strands[0]] == (lay_a.placed[strands[-1]].cx, lay_a.placed[strands[-1]].cy)
    assert len(lay.sheets) == len(alone.sheets) == 1
    mine, own = _strand_order(lay, lay.sheets[0]), _strand_order(alone, alone.sheets[0])
    assert mine in (own, own[::-1])
    page = SimpleNamespace(axes=[SimpleNamespace(patches=[])])
    found = figcheck.science_checks(SimpleNamespace(layout=lay, sses=sses_b, bb=bb_b, figure=page))
    assert found["S1"] == found["S2"] == found["S3"] == [] and found["S6"] == []


def test_homologue_takes_its_mates_places_names_and_angles():
    """Haemoglobin α (2HHB chain A) drawn like β (chain B): every mated helix carries its mate's label and drawn
    angle, and sits nearer its mate's place than drawn alone does (in the same frame)."""
    path = DATA / "2HHB.cif"
    lay, sses_b, _, lay_a, like = drawn_like(path, "B", path, "A")
    assert len(like.target) >= 6
    alone, _, _ = laid_out_like(path, None, lay.frame, chains=["A"])
    near, far = [], []
    for s in sses_b:
        if s.id not in like.target:
            continue
        p = lay.placed[s.id]
        assert p.angle == like.angle[s.id]
        assert p.label == like.label[f"res:{s.chain}:{s.first.tag}-{s.last.tag}"]
        tx, ty = like.target[s.id]
        near.append(np.hypot(p.cx - tx, p.cy - ty))
        far.append(np.hypot(alone.placed[s.id].cx - tx, alone.placed[s.id].cy - ty))
    assert np.median(near) < np.median(far)


_SIGNATURE = """
import json, sys
sys.path.insert(0, sys.argv[1])
from test_like import DATA, drawn_like, signature
lay = drawn_like(DATA / "2HHB.cif", "B", DATA / "2HHB.cif", "A")[0]
print(json.dumps(signature(lay)))
"""


def test_like_figure_does_not_depend_on_the_hash_seed():
    """A12 for a figure drawn like a reference: the same layout under PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        done = subprocess.run(
            [sys.executable, "-c", _SIGNATURE, str(Path(__file__).parent)],
            env=env,
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout)
    assert json.loads(runs[0]) and runs[0] == runs[1]
