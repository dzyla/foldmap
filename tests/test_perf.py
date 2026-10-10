"""Speed-ups that must not change results: each fast path checked against the plain one it replaced."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import foldmap.layout as L
from foldmap import io
from foldmap.frame import Frame
from foldmap.layout import build_layout
from helpers import fake_sse

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).parent.parent

_SIGNATURE = """
import json, sys
from foldmap.cli import make_layout
out = {}
for f in sys.argv[1:]:
    lay, _, _ = make_layout(f)
    out[f] = {k: [p.cx, p.cy, p.angle] for k, p in sorted(lay.placed.items())}
print(json.dumps(out))
"""


def test_layout_does_not_depend_on_the_hash_seed():
    files = [str(DATA / f"{p}.cif") for p in ("2HHB", "1LMB", "9TX7")]
    runs = []
    for seed in ("0", "12345", "777"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run(
            [sys.executable, "-c", _SIGNATURE, *files], env=env, capture_output=True, text=True, check=True
        )
        runs.append(json.loads(done.stdout))
    assert runs[0] and all(runs[0][f] for f in runs[0])
    assert runs[0] == runs[1] == runs[2]


def _blas_threads():
    from threadpoolctl import threadpool_info

    return [i["num_threads"] for i in threadpool_info() if i["user_api"] == "blas"]


def _pile(n):
    """n helices projected onto one point: the push-apart must untangle them by itself."""
    f = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    return [fake_sse("H", 10 * k, 10 * k + 8, centroid=(0, 0, 0), axis=(0.2, 0.1, 1)) for k in range(n)], f


def test_layout_runs_blas_on_one_thread_and_restores_it(monkeypatch):
    before = _blas_threads()
    seen = []
    real = L._attract

    def spy(*a, **k):
        seen.append(_blas_threads())
        return real(*a, **k)

    monkeypatch.setattr(L, "_attract", spy)
    sses, f = _pile(4)
    build_layout(sses, [], f, "projected")
    assert seen and all(t == 1 for t in seen[0])
    assert _blas_threads() == before


def test_overlapping_layouts_in_threads_keep_the_limit_and_restore_it():
    import threading

    before = _blas_threads()
    inside, first_out = threading.Barrier(2), threading.Event()
    seen = []

    def first():
        with L._one_blas_thread():
            inside.wait()
        first_out.set()

    def second():
        with L._one_blas_thread():
            inside.wait()
            first_out.wait()
            seen.append(_blas_threads())  # the other layout has finished: this one still runs on one thread

    threads = [threading.Thread(target=first), threading.Thread(target=second)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert seen and all(n == 1 for n in seen[0])
    assert _blas_threads() == before


def test_a_figure_reads_its_file_once(monkeypatch):
    from foldmap import cli

    calls = []
    real = io.read_model

    def counting(*a, **k):
        calls.append(a[0])
        return real(*a, **k)

    monkeypatch.setattr(io, "read_model", counting)
    monkeypatch.setattr(cli, "read_model", counting, raising=False)
    cli.make_layout(DATA / "1LMB.cif")  # protein, DNA and links all come from the one read
    assert len(calls) == 1
    cli.summarize(str(DATA / "1LMB.cif"))
    assert len(calls) == 2


def _tiled(pdb, copies):
    """Model 1 of a test entry, repeated `copies` times 80 Å apart along x (new chain names)."""
    import gemmi

    st = gemmi.read_structure(str(DATA / f"{pdb}.cif"))
    st.setup_entities()
    model = st[0]
    base = [c.clone() for c in model]
    for c in list(model):
        model.remove_chain(c.name)
    names = (a + b for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" for b in "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    for k in range(copies):
        for c in base:
            new = c.clone()
            new.name = next(names)
            for r in new:
                for a in r:
                    a.pos = gemmi.Position(a.pos.x + 80 * k, a.pos.y, a.pos.z)
            model.add_chain(new)
    return model


@pytest.mark.parametrize("pdb, copies", [("1C3W", 1), ("2HHB", 1), ("1BL8", 1), ("25IN", 1), ("25IN", 6)])
def test_bonded_ligand_parts_equal_the_all_against_all_search(pdb, copies):
    parts = io._ligand_parts(_tiled(pdb, copies))
    assert parts
    brute = {
        (a, b)
        for a in range(len(parts))
        for b in range(a + 1, len(parts))
        if not parts[a][3]
        and not parts[b][3]
        and any(x.pos.dist(y.pos) < io._BOND for x in parts[a][2] for y in parts[b][2])
    }
    assert io._bonded_parts(parts) == brute


def test_no_ligands_need_no_search():
    assert io._bonded_parts([]) == set()
