# foldmap Step 1: Structure Parsing, Secondary Structure, Sheet Ordering

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `foldmap summary FILE` command that reads a PDB/mmCIF, assigns helices and strands, groups strands into ordered sheets (including sheets completed by strands from another chain), and prints the result. This is build step 1 of 4 in the spec; steps 2-4 (layout/routing/rendering, annotation layers, membrane mode) get their own plans.

**Architecture:** Pipeline of small modules with plain dataclasses between them: `io.load_backbone` -> `dssp.assign_dssp` -> `ss.build_sses` -> `sheets.build_sheets`, driven by `cli`. Secondary structure is always computed by a built-in numpy Kabsch-Sander DSSP, because the sheet-pairing graph needs backbone H-bonds anyway and the deposited sheet records are not reliable (in 6ZS5 the same strand appears in two sheets). The deposited annotation is used only as a test cross-check.

**Tech Stack:** Python >= 3.11, gemmi (parsing), numpy, scipy (`cKDTree`), networkx (sheet graph), pytest.

**Spec:** `docs/superpowers/specs/2026-10-08-protein-topology-plotter-design.md` (sections 3, 6 step 1, 7)

## Global Constraints

- Package name `foldmap`, source in `src/foldmap/`, tests in `tests/`.
- Dependencies: gemmi, numpy, scipy, networkx, matplotlib, PyYAML (matplotlib/PyYAML are declared now for later steps); no `mkdssp` dependency.
- Only model 1 is read; only amino-acid residues with a complete N, CA, C, O backbone are kept; the first altloc is used (`res.find_atom(name, "*")`). *This replaces the spec's "highest occupancy" wording; the spec has been updated.*
- Interfaces between stages are the dataclasses in `model.py`; stages never pass gemmi objects.
- The project directory is not a git repository and no `git init` is run, so there are no commit steps; each task ends with the full test suite passing.
- Test fixtures (already in `tests/data/`): `1UBQ.cif` (ubiquitin: one helix 23-34, mixed 5-strand sheet, known strand order beta2-1-5-3-4), `6ZS5.cif` (uromodulin filament core, chains A and D plus glycan chains B, C), `6ZYA.cif`.

## Review Focus

Inputs the spec implies but a happy-path test would miss, most likely first:

1. A PDB-format file (not mmCIF) must give identical residues and coordinates. Pinned in Task 1 (`test_pdb_format_matches_mmcif`).
2. A multi-model file (NMR) must use model 1 only, not concatenate models. Pinned in Task 1 (`test_only_first_model_used`).
3. A residue with a missing backbone atom must be skipped and must break the chain there, so no peptide bond or H-bond is invented across the gap. Pinned in Task 1 (`test_residue_missing_backbone_atom_is_skipped_and_breaks_chain`).
4. A fragment with no helices or strands (a tripeptide) must print zero counts and no sheet, not crash. Pinned in Task 5 (`test_fragment_without_secondary_structure`).
5. Sheets whose pairing graph is a cycle (barrel) or is branched must be flagged `closed` / `ambiguous` rather than silently mis-ordered, and weak single-bridge contacts must not merge sheets. Pinned in Task 4 (`test_cycle_is_closed_and_cut`, `test_branch_is_ambiguous`, `test_weak_pairing_and_isolated_strand_make_separate_sheets`).

---

### Task 1: Scaffold, data model, structure loading

**Files:**
- Create: `pyproject.toml`, `src/foldmap/__init__.py`, `src/foldmap/model.py`, `src/foldmap/io.py`
- Create: `tests/conftest.py`, `tests/test_io.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `model.ResLabel(chain, seq, icode, name)`, `model.Bridge(i, j, kind)`, `model.Backbone(labels, xyz, prev)` with `.nxt`, `.ca`, `len()`, `model.SSE(kind, chain, start, end, first, last, centroid, axis)` with `.id`, `model.Sheet(strands, directions, pair_kinds, closed, ambiguous)`; `io.load_backbone(path) -> Backbone` raising `FileNotFoundError` for a missing file and `ValueError("... no protein residues ...")` when nothing usable is found.

- [ ] **Step 1: Check fixtures and create the scaffold**

Run: `ls tests/data` (expect `1UBQ.cif 6ZS5.cif 6ZYA.cif`), then create these files.

**`pyproject.toml`**

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "foldmap"
version = "0.1.0"
description = "Publication-grade protein topology diagrams"
requires-python = ">=3.11"
dependencies = ["gemmi>=0.7", "numpy", "scipy", "networkx", "matplotlib", "PyYAML"]

[project.optional-dependencies]
test = ["pytest"]

[project.scripts]
foldmap = "foldmap.cli:main"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

**`src/foldmap/__init__.py`**

```python
"""foldmap: publication-grade protein topology diagrams."""

__version__ = "0.1.0"
```

**`src/foldmap/model.py`**

```python
"""Plain data containers shared by all pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np


class ResLabel(NamedTuple):
    chain: str
    seq: int
    icode: str
    name: str


class Bridge(NamedTuple):
    """Backbone H-bond bridge between residues i < j. kind: 'P' parallel, 'A' antiparallel."""

    i: int
    j: int
    kind: str


@dataclass(eq=False)
class Backbone:
    """Protein backbone of model 1. xyz[i] is rows N, CA, C, O of residue i.

    prev[i] is True when residue i is peptide-bonded to residue i-1.
    """

    labels: list[ResLabel]
    xyz: np.ndarray  # (N, 4, 3)
    prev: np.ndarray  # (N,) bool

    def __len__(self) -> int:
        return len(self.labels)

    @property
    def nxt(self) -> np.ndarray:
        out = np.zeros(len(self), bool)
        out[:-1] = self.prev[1:]
        return out

    @property
    def ca(self) -> np.ndarray:
        return self.xyz[:, 1, :]


@dataclass(eq=False)
class SSE:
    """A helix ('H') or strand ('E'). start/end are inclusive residue indices."""

    kind: str
    chain: str
    start: int
    end: int
    first: ResLabel
    last: ResLabel
    centroid: np.ndarray
    axis: np.ndarray  # unit vector, points N-terminal -> C-terminal end

    @property
    def id(self) -> str:
        return f"{self.chain}:{self.first.seq}-{self.last.seq}"

    def __len__(self) -> int:
        return self.end - self.start + 1


@dataclass(eq=False)
class Sheet:
    """Strands in sheet order. directions[k] is +1 if strands[k] runs the same way as
    strands[0], else -1. pair_kinds[k] relates strands[k] and strands[k+1] ('P', 'A' or '?')."""

    strands: list[SSE]
    directions: list[int]
    pair_kinds: list[str]
    closed: bool = False  # barrel: a cycle was cut at its weakest pairing
    ambiguous: bool = False  # branched: ordering is a heuristic
```

**`tests/conftest.py`**

```python
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def ubq():
    return DATA / "1UBQ.cif"


@pytest.fixture(scope="session")
def zs5():
    return DATA / "6ZS5.cif"


@pytest.fixture(scope="session")
def zya():
    return DATA / "6ZYA.cif"
```


- [ ] **Step 2: Create the environment**

Run: `python3 -m venv --system-site-packages .venv && . .venv/bin/activate && pip install -q -e '.[test]'`
Expected: installs without error; `python -c "import foldmap"` succeeds.

- [ ] **Step 3: Write the failing tests**

**`tests/test_io.py`**

```python
import gemmi
import numpy as np
import pytest

from foldmap.io import load_backbone


def test_ubiquitin_backbone(ubq):
    bb = load_backbone(ubq)
    assert len(bb) == 76
    assert bb.xyz.shape == (76, 4, 3)
    assert bb.prev.sum() == 75 and not bb.prev[0]
    assert bb.labels[0].name == "MET" and bb.labels[0].seq == 1
    assert np.isfinite(bb.xyz).all()


def test_multichain_break_and_glycans_skipped(zs5):
    bb = load_backbone(zs5)
    assert len(bb) == 257  # glycan chains B, C contribute nothing
    assert {l.chain for l in bb.labels} == {"A", "D"}
    first_d = next(i for i, l in enumerate(bb.labels) if l.chain == "D")
    assert not bb.prev[first_d]


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_backbone(tmp_path / "nope.cif")


def test_no_protein_raises(tmp_path):
    f = tmp_path / "water.pdb"
    f.write_text("HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nEND\n")
    with pytest.raises(ValueError, match="no protein residues"):
        load_backbone(f)


def test_pdb_format_matches_mmcif(ubq, tmp_path):
    out = tmp_path / "u.pdb"
    gemmi.read_structure(str(ubq)).write_pdb(str(out))
    a, b = load_backbone(ubq), load_backbone(out)
    assert a.labels == b.labels
    assert np.allclose(a.xyz, b.xyz, atol=1e-3)


def test_only_first_model_used(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    st.add_model(st[0])  # second model; reading both would give 152 residues
    out = tmp_path / "two_models.cif"
    st.make_mmcif_document().write_file(str(out))
    assert len(load_backbone(out)) == 76


def test_residue_missing_backbone_atom_is_skipped_and_breaks_chain(ubq, tmp_path):
    st = gemmi.read_structure(str(ubq))
    res = next(r for r in st[0]["A"] if r.seqid.num == 30)
    res.remove_atom("O", "*")
    out = tmp_path / "gap.cif"
    st.make_mmcif_document().write_file(str(out))
    bb = load_backbone(out)
    assert len(bb) == 75
    i31 = next(i for i, l in enumerate(bb.labels) if l.seq == 31)
    assert not bb.prev[i31]  # 29 and 31 are no longer bonded
```

- [ ] **Step 4: Run to verify they fail**

Run: `python -m pytest tests/test_io.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'foldmap.io'`.

- [ ] **Step 5: Implement**

**`src/foldmap/io.py`**

```python
"""Read PDB/mmCIF into a Backbone (first model, first altloc, amino acids only)."""

from __future__ import annotations

from pathlib import Path

import gemmi
import numpy as np

from .model import Backbone, ResLabel

_BB_ATOMS = ("N", "CA", "C", "O")
_PEPTIDE_BOND_MAX = 2.0  # Å, C(i-1)-N(i)


def load_backbone(path: str | Path) -> Backbone:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        raise ValueError(f"{path}: no models")
    labels: list[ResLabel] = []
    coords: list[list[list[float]]] = []
    for chain in st[0]:
        for res in chain:
            if res.het_flag != "A":
                continue
            atoms = [res.find_atom(a, "*") for a in _BB_ATOMS]
            if any(a is None for a in atoms):
                continue
            labels.append(ResLabel(chain.name, res.seqid.num, res.seqid.icode.strip(), res.name))
            coords.append([a.pos.tolist() for a in atoms])
    if not labels:
        raise ValueError(f"{path}: no protein residues with a complete backbone")
    xyz = np.asarray(coords, float)
    prev = np.zeros(len(labels), bool)
    for i in range(1, len(labels)):
        if labels[i].chain == labels[i - 1].chain:
            prev[i] = np.linalg.norm(xyz[i - 1, 2] - xyz[i, 0]) < _PEPTIDE_BOND_MAX
    return Backbone(labels, xyz, prev)
```

- [ ] **Step 6: Run to verify they pass**

Run: `python -m pytest -q`
Expected: 7 passed.

---

### Task 2: Kabsch-Sander DSSP

**Files:**
- Create: `src/foldmap/dssp.py`, `tests/test_dssp.py`

**Interfaces:**
- Consumes: `io.load_backbone`, `model.Backbone`, `model.Bridge`.
- Produces: `dssp.assign_dssp(bb: Backbone) -> DsspResult(ss: str, bridges: list[Bridge])`. `ss` has one code per residue: `H G I E B -`. `bridges` holds every residue pair (i < j, j >= i+3, within 9 Å CA) with a parallel (`P`) or antiparallel (`A`) bridge H-bond pattern, including pairs across chains.

- [ ] **Step 1: Write the failing tests**

**`tests/test_dssp.py`**

```python
import gemmi
import numpy as np

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone


def test_ubiquitin_helix_and_strands(ubq):
    bb = load_backbone(ubq)
    res = assign_dssp(bb)
    assert len(res.ss) == len(bb)
    assert set(res.ss[22:34]) == {"H"}  # helix 23-34
    assert set(res.ss[2:7]) == {"E"}  # beta1 3-7
    assert set(res.ss[66:71]) == {"E"}  # beta5 67-71
    assert res.ss.count("H") == 12


def _deposited_strand_residues(path):
    st = gemmi.read_structure(str(path))
    out = set()
    for sheet in st.sheets:
        for s in sheet.strands:
            for k in range(s.start.res_id.seqid.num, s.end.res_id.seqid.num + 1):
                out.add((s.start.chain_name, k))
    return out


def test_agrees_with_deposited_strands(zs5):
    bb = load_backbone(zs5)
    res = assign_dssp(bb)
    mine = {(bb.labels[i].chain, bb.labels[i].seq) for i, c in enumerate(res.ss) if c == "E"}
    dep = _deposited_strand_residues(zs5)
    assert len(mine & dep) / len(mine) >= 0.95  # what we call a strand, the file agrees with
    assert len(mine & dep) / len(dep) >= 0.85  # and we recover most of the file's strands


def test_cross_chain_bridges_found(zs5):
    bb = load_backbone(zs5)
    res = assign_dssp(bb)
    cross = [b for b in res.bridges if bb.labels[b.i].chain != bb.labels[b.j].chain]
    assert len(cross) >= 10  # linker strands complete the neighbouring sheet


def test_bridge_invariants(ubq):
    res = assign_dssp(load_backbone(ubq))
    assert all(b.j >= b.i + 3 and b.kind in "PA" for b in res.bridges)
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_dssp.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'foldmap.dssp'`.

- [ ] **Step 3: Implement**

**`src/foldmap/dssp.py`**

```python
"""Kabsch-Sander DSSP: H-bond energies, helices, bridges.

Not implemented: bulge merging and the turn/bend (T, S) classes. Output codes:
H alpha, G 3-10, I pi, E strand (ladder of >= 2 bridges), B isolated bridge, '-' other.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.spatial import cKDTree

from .model import Backbone, Bridge

_Q1Q2F = 0.084 * 332.0  # kcal/mol * Å
_HBOND_E = -0.5
_CA_CUTOFF = 9.0


class DsspResult(NamedTuple):
    ss: str
    bridges: list[Bridge]


def _hydrogen_positions(bb: Backbone) -> np.ndarray:
    n = len(bb)
    h = np.full((n, 3), np.nan)
    for i in range(1, n):
        if bb.prev[i] and bb.labels[i].name != "PRO":
            v = bb.xyz[i - 1, 2] - bb.xyz[i - 1, 3]  # C(i-1) - O(i-1)
            h[i] = bb.xyz[i, 0] + v / np.linalg.norm(v)
    return h


def _hbond_set(bb: Backbone) -> set[tuple[int, int]]:
    """Set of (acceptor, donor): CO of acceptor bonded to NH of donor."""
    xyz, n = bb.xyz, len(bb)
    h = _hydrogen_positions(bb)
    close = np.array(sorted(cKDTree(bb.ca).query_pairs(_CA_CUTOFF)))
    if close.size == 0:
        return set()
    both = np.vstack([close, close[:, ::-1]])
    don, acc = both[:, 0], both[:, 1]
    keep = ~np.isnan(h[don, 0]) & ~((acc == don - 1) & bb.prev[don])
    don, acc = don[keep], acc[keep]
    hd, nd = h[don], xyz[don, 0]
    oa, ca = xyz[acc, 3], xyz[acc, 2]
    d = lambda p, q: np.linalg.norm(p - q, axis=1)  # noqa: E731
    e = _Q1Q2F * (1 / d(oa, nd) + 1 / d(ca, hd) - 1 / d(oa, hd) - 1 / d(ca, nd))
    order = np.argsort(e)
    used_d, used_a, bonds = {}, {}, set()
    for k in order:
        if e[k] >= _HBOND_E:
            break
        di, ai = int(don[k]), int(acc[k])
        if used_d.get(di, 0) < 2 and used_a.get(ai, 0) < 2:
            used_d[di] = used_d.get(di, 0) + 1
            used_a[ai] = used_a.get(ai, 0) + 1
            bonds.add((ai, di))
    return bonds


def assign_dssp(bb: Backbone) -> DsspResult:
    n, prev, nxt = len(bb), bb.prev, bb.nxt
    bonds = _hbond_set(bb)

    def hb(a: int | None, b: int | None) -> bool:
        return a is not None and b is not None and (a, b) in bonds

    def before(i: int) -> int | None:
        return i - 1 if prev[i] else None

    def after(i: int) -> int | None:
        return i + 1 if nxt[i] else None

    ss = ["-"] * n
    for k, code in ((4, "H"), (3, "G"), (5, "I")):
        turn = [
            i + k < n and (i, i + k) in bonds and all(nxt[i + m] for m in range(k)) for i in range(n)
        ]
        for i in range(1, n - k):
            if turn[i - 1] and turn[i]:
                for m in range(k):
                    if ss[i + m] == "-" or (code == "H" and ss[i + m] in "GI"):
                        ss[i + m] = code

    bridges: list[Bridge] = []
    close = cKDTree(bb.ca).query_pairs(_CA_CUTOFF)
    for i, j in sorted(close):
        if j < i + 3:
            continue
        im, ip, jm, jp = before(i), after(i), before(j), after(j)
        if (hb(im, j) and hb(j, ip)) or (hb(jm, i) and hb(i, jp)):
            bridges.append(Bridge(i, j, "P"))
        if (hb(i, j) and hb(j, i)) or (hb(im, jp) and hb(jm, ip)):
            bridges.append(Bridge(i, j, "A"))

    known = set(bridges)
    in_ladder, in_bridge = set(), set()
    for i, j, kind in bridges:
        in_bridge.update((i, j))
        if i + 1 < n and prev[i + 1]:
            nj = j + 1 if kind == "P" else j - 1
            if 0 <= nj < n and Bridge(i + 1, nj, kind) in known:
                in_ladder.update((i, j, i + 1, nj))
    for i in range(n):
        if ss[i] == "-":
            if i in in_ladder:
                ss[i] = "E"
            elif i in in_bridge:
                ss[i] = "B"
    return DsspResult("".join(ss), bridges)
```

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest -q`
Expected: 11 passed.

---

### Task 3: SSE builder

**Files:**
- Create: `src/foldmap/ss.py`, `tests/test_ss.py`

**Interfaces:**
- Consumes: `Backbone`, the `ss` string from `assign_dssp`.
- Produces: `ss.build_sses(bb: Backbone, ss: str) -> list[SSE]`, sorted by `start`. Strands: runs of `E`, length >= 2. Helices: runs of `H/G/I` with >= 4 `H` residues or length >= 6. Runs never cross a chain break. `SSE.axis` is a unit vector pointing N-terminal to C-terminal end; `SSE.centroid` is the mean CA position.

- [ ] **Step 1: Write the failing tests**

**`tests/test_ss.py`**

```python
import numpy as np

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.ss import build_sses


def test_ubiquitin_sses(ubq):
    bb = load_backbone(ubq)
    sses = build_sses(bb, assign_dssp(bb).ss)
    assert [(s.kind, s.id) for s in sses] == [
        ("E", "A:2-7"),
        ("E", "A:12-16"),
        ("H", "A:23-34"),
        ("E", "A:41-45"),
        ("E", "A:48-49"),
        ("E", "A:66-71"),
    ]  # the 3-residue 3-10 turn at 38-40 is not drawn


def test_axis_and_centroid(ubq):
    bb = load_backbone(ubq)
    for s in build_sses(bb, assign_dssp(bb).ss):
        assert abs(np.linalg.norm(s.axis) - 1) < 1e-9
        span = bb.ca[s.end] - bb.ca[s.start]
        assert np.dot(span, s.axis) > 0  # points N -> C
        assert np.allclose(s.centroid, bb.ca[s.start : s.end + 1].mean(axis=0))


def test_sses_never_span_chain_breaks(zs5):
    bb = load_backbone(zs5)
    for ss in ("E" * len(bb), "H" * len(bb)):
        for s in build_sses(bb, ss):
            assert all(bb.labels[i].chain == s.chain for i in range(s.start, s.end + 1))
            assert all(bb.prev[i] for i in range(s.start + 1, s.end + 1))


def test_short_runs_filtered(ubq):
    bb = load_backbone(ubq)
    ss = list("-" * len(bb))
    ss[5] = "E"  # lone strand residue
    ss[10:13] = "GGG"  # 3-10 turn too short
    ss[20:23] = "HHH"  # alpha run too short
    assert build_sses(bb, "".join(ss)) == []
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ss.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'foldmap.ss'`.

- [ ] **Step 3: Implement**

**`src/foldmap/ss.py`**

```python
"""Turn per-residue DSSP codes into helix/strand SSE objects with an axis and centroid."""

from __future__ import annotations

import numpy as np

from .model import Backbone, SSE

_MIN_STRAND = 2
_MIN_HELIX_H = 4  # an alpha run needs this many H residues ...
_MIN_HELIX_ANY = 6  # ... or this many helix-class residues in total (3-10 / pi)


def _runs(ss: str, classes: str, prev: np.ndarray):
    """Maximal chain-contiguous runs of residues whose code is in `classes`."""
    start = None
    for i, c in enumerate(ss):
        if c in classes:
            if start is None:
                start = i
            elif not prev[i]:
                yield start, i - 1
                start = i
        elif start is not None:
            yield start, i - 1
            start = None
    if start is not None:
        yield start, len(ss) - 1


def _fit_axis(ca: np.ndarray) -> np.ndarray:
    centred = ca - ca.mean(axis=0)
    axis = np.linalg.svd(centred, full_matrices=False)[2][0]
    head, tail = ca[: min(2, len(ca))].mean(axis=0), ca[-min(2, len(ca)) :].mean(axis=0)
    if np.dot(tail - head, axis) < 0:
        axis = -axis
    return axis / np.linalg.norm(axis)


def build_sses(bb: Backbone, ss: str) -> list[SSE]:
    out: list[SSE] = []
    for kind, classes in (("E", "E"), ("H", "HGI")):
        for a, b in _runs(ss, classes, bb.prev):
            n = b - a + 1
            if kind == "E" and n < _MIN_STRAND:
                continue
            if kind == "H" and not (ss[a : b + 1].count("H") >= _MIN_HELIX_H or n >= _MIN_HELIX_ANY):
                continue
            ca = bb.ca[a : b + 1]
            out.append(
                SSE(kind, bb.labels[a].chain, a, b, bb.labels[a], bb.labels[b], ca.mean(axis=0), _fit_axis(ca))
            )
    out.sort(key=lambda s: s.start)
    return out
```

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest -q`
Expected: 15 passed.

---

### Task 4: Sheet grouping and ordering

**Files:**
- Create: `src/foldmap/sheets.py`, `tests/test_sheets.py`

**Interfaces:**
- Consumes: `list[SSE]` from `build_sses`, `list[Bridge]` from `assign_dssp`.
- Produces: `sheets.build_sheets(sses, bridges) -> list[Sheet]`. Two strands are neighbours when at least 2 bridges join them; the pair kind is the majority of its bridges (ties give `A`). A sheet is a connected component, ordered side by side: a path is walked end to end; a cycle is cut at its weakest pairing (`closed=True`); a branched tree is ordered along its longest path with side branches inserted beside their parent (`ambiguous=True`). `directions[k]` is +1/-1 relative to `strands[0]`; `pair_kinds[k]` is `P`, `A`, or `?` when consecutive strands in the order are not directly paired. Sheets are returned sorted by their lowest strand index; an unpaired strand is a one-strand sheet. Helices are ignored.

- [ ] **Step 1: Write the failing tests**

**`tests/test_sheets.py`**

```python
import numpy as np

from foldmap.dssp import assign_dssp
from foldmap.io import load_backbone
from foldmap.model import Bridge, ResLabel, SSE
from foldmap.sheets import build_sheets
from foldmap.ss import build_sses


def fake(start, end):
    lab = ResLabel("A", start, "", "ALA")
    return SSE("E", "A", start, end, lab, lab, np.zeros(3), np.array([1.0, 0, 0]))


def ladder(a, b, kind, n=3):
    """n bridges between strands a and b (SSE objects)."""
    return [Bridge(a.start + k, b.start + k, kind) for k in range(n)]


def test_path_order_and_directions():
    s1, s2, s3 = fake(0, 9), fake(20, 29), fake(40, 49)
    sheets = build_sheets([s1, s2, s3], ladder(s1, s2, "A") + ladder(s2, s3, "P"))
    (sh,) = sheets
    assert [s.start for s in sh.strands] == [0, 20, 40]
    assert sh.directions == [1, -1, -1]
    assert sh.pair_kinds == ["A", "P"]
    assert not sh.closed and not sh.ambiguous


def test_cycle_is_closed_and_cut():
    s = [fake(0, 9), fake(20, 29), fake(40, 49)]
    br = ladder(s[0], s[1], "A", 4) + ladder(s[1], s[2], "A", 4) + ladder(s[0], s[2], "A", 2)
    (sh,) = build_sheets(s, br)
    assert sh.closed and len(sh.strands) == 3
    assert [x.start for x in sh.strands] in ([0, 20, 40], [40, 20, 0])  # weakest edge (0-2) cut


def test_branch_is_ambiguous():
    c, a, b, d = fake(0, 9), fake(20, 29), fake(40, 49), fake(60, 69)
    br = ladder(c, a, "A") + ladder(c, b, "A") + ladder(c, d, "A")
    (sh,) = build_sheets([c, a, b, d], br)
    assert sh.ambiguous and len(sh.strands) == 4
    assert sorted(s.start for s in sh.strands) == [0, 20, 40, 60]


def test_weak_pairing_and_isolated_strand_make_separate_sheets():
    s1, s2, s3 = fake(0, 9), fake(20, 29), fake(40, 49)
    sheets = build_sheets([s1, s2, s3], ladder(s1, s2, "A", 1))  # 1 bridge < minimum
    assert [len(sh.strands) for sh in sheets] == [1, 1, 1]
    assert all(sh.pair_kinds == [] and sh.directions == [1] for sh in sheets)


def test_helices_ignored():
    h = fake(0, 9)
    h.kind = "H"
    assert build_sheets([h], []) == []


def _pipeline(path):
    bb = load_backbone(path)
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss)
    return bb, sses, build_sheets(sses, d.bridges)


def test_ubiquitin_mixed_sheet(ubq):
    _, _, (sh,) = _pipeline(ubq)
    ids = [s.id for s in sh.strands]
    known = ["A:12-16", "A:2-7", "A:66-71", "A:41-45", "A:48-49"]  # beta2-1-5-3-4
    assert ids in (known, known[::-1])
    kinds = sh.pair_kinds if ids == known else sh.pair_kinds[::-1]
    assert kinds == ["A", "P", "A", "A"]  # only the 1-5 pairing is parallel


def test_uromodulin_sheets(zs5):
    bb, sses, sheets = _pipeline(zs5)
    strands = [s for s in sses if s.kind == "E"]
    assert sum(len(sh.strands) for sh in sheets) == len(strands)  # every strand in exactly one sheet
    assert len({s.id for sh in sheets for s in sh.strands}) == len(strands)
    assert any({s.chain for s in sh.strands} == {"A", "D"} for sh in sheets)  # complementation
    first = [s.id for s in sheets[0].strands]  # matches deposited sheet AA1 in the mmCIF
    assert first == ["A:331-333", "A:339-345", "A:377-383", "A:369-371"]
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_sheets.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'foldmap.sheets'`.

- [ ] **Step 3: Implement**

**`src/foldmap/sheets.py`**

```python
"""Group strands into sheets from bridge pairing and order the strands side by side."""

from __future__ import annotations

from collections import Counter, defaultdict

import networkx as nx

from .model import Bridge, SSE, Sheet

_MIN_PAIRS = 2  # bridges needed between two strands to count them as neighbours


def _edge_table(strands: list[SSE], bridges: list[Bridge]):
    owner: dict[int, int] = {}
    for k, s in enumerate(strands):
        for r in range(s.start, s.end + 1):
            owner[r] = k
    table: dict[tuple[int, int], Counter] = defaultdict(Counter)
    for i, j, kind in bridges:
        a, b = owner.get(i), owner.get(j)
        if a is not None and b is not None and a != b:
            table[(min(a, b), max(a, b))][kind] += 1
    return {e: c for e, c in table.items() if sum(c.values()) >= _MIN_PAIRS}


def _majority(c: Counter) -> str:
    return "P" if c["P"] > c["A"] else "A"


def _order_tree(tree: nx.Graph, first: int) -> tuple[list[int], bool]:
    """Order tree nodes along its longest path; branches are inserted beside their parent."""
    degrees = dict(tree.degree)
    if tree.number_of_nodes() == 1:
        return [first], False
    if max(degrees.values()) <= 2:
        end = min((n for n, d in degrees.items() if d == 1))
        return list(nx.dfs_preorder_nodes(tree, end)), False
    far = lambda src: max(nx.single_source_shortest_path_length(tree, src).items(), key=lambda kv: (kv[1], -kv[0]))[0]  # noqa: E731
    a = far(min(tree.nodes))
    b = far(a)
    path = nx.shortest_path(tree, a, b)
    order = list(path)
    placed = set(order)
    for node in nx.bfs_tree(tree, a):
        if node in placed:
            continue
        parent = next(p for p in tree.neighbors(node) if p in placed)
        order.insert(order.index(parent) + 1, node)
        placed.add(node)
    return order, True


def build_sheets(sses: list[SSE], bridges: list[Bridge]) -> list[Sheet]:
    strands = [s for s in sses if s.kind == "E"]
    table = _edge_table(strands, bridges)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(strands)))
    for (a, b), c in table.items():
        graph.add_edge(a, b, weight=sum(c.values()), kind=_majority(c))

    sheets: list[Sheet] = []
    for comp in sorted(nx.connected_components(graph), key=min):
        sub = graph.subgraph(comp)
        tree = nx.maximum_spanning_tree(sub, weight="weight")
        closed = tree.number_of_edges() < sub.number_of_edges()
        order, ambiguous = _order_tree(tree, min(comp))
        direction = {order[0]: 1}
        for parent, child in nx.bfs_edges(tree, order[0]):
            flip = -1 if tree.edges[parent, child]["kind"] == "A" else 1
            direction[child] = direction[parent] * flip
        kinds = [
            tree.edges[a, b]["kind"] if tree.has_edge(a, b) else "?" for a, b in zip(order, order[1:])
        ]
        sheets.append(
            Sheet([strands[k] for k in order], [direction[k] for k in order], kinds, closed, ambiguous)
        )
    return sheets
```

- [ ] **Step 4: Run to verify they pass**

Run: `python -m pytest -q`
Expected: 22 passed.

---

### Task 5: `foldmap summary` command

**Files:**
- Create: `src/foldmap/cli.py`, `src/foldmap/__main__.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: `load_backbone`, `assign_dssp`, `build_sses`, `build_sheets`.
- Produces: `cli.summarize(path: str) -> str` and `cli.main(argv) -> int` (0 on success; 1 with `foldmap: <message>` on stderr for a missing or unusable file). Installed as the `foldmap` console script and runnable as `python -m foldmap`.

- [ ] **Step 1: Write the failing tests**

**`tests/test_cli.py`**

```python
from foldmap.cli import main


def test_summary_output(ubq, capsys):
    assert main(["summary", str(ubq)]) == 0
    out = capsys.readouterr().out
    assert "76 residues, 1 chain(s)" in out
    assert "helix  A:23-34" in out
    assert "sheet 1: 5 strand(s)" in out


def test_missing_file_exit_code(tmp_path, capsys):
    assert main(["summary", str(tmp_path / "nope.cif")]) == 1
    assert "foldmap:" in capsys.readouterr().err


def test_fragment_without_secondary_structure(ubq, tmp_path, capsys):
    import gemmi

    full = tmp_path / "full.pdb"
    gemmi.read_structure(str(ubq)).write_pdb(str(full))
    keep = [l for l in full.read_text().splitlines() if l.startswith("ATOM") and int(l[22:26]) <= 3]
    tri = tmp_path / "tri.pdb"
    tri.write_text("\n".join(keep) + "\nEND\n")
    assert main(["summary", str(tri)]) == 0
    out = capsys.readouterr().out
    assert "3 residues" in out and "helices: 0  strands: 0" in out and "sheet" not in out
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_cli.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'foldmap.cli'`.

- [ ] **Step 3: Implement**

**`src/foldmap/cli.py`**

```python
"""Command line: `foldmap summary FILE` prints secondary structure and sheet ordering."""

from __future__ import annotations

import argparse
import sys

from .dssp import assign_dssp
from .io import load_backbone
from .sheets import build_sheets
from .ss import build_sses


def summarize(path: str) -> str:
    bb = load_backbone(path)
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss)
    sheets = build_sheets(sses, dssp.bridges)
    lines = [f"{path}: {len(bb)} residues, {len(set(l.chain for l in bb.labels))} chain(s)"]
    lines.append(f"helices: {sum(s.kind == 'H' for s in sses)}  strands: {sum(s.kind == 'E' for s in sses)}")
    for s in sses:
        if s.kind == "H":
            lines.append(f"  helix  {s.id}")
    for n, sh in enumerate(sheets, 1):
        flags = ("closed " if sh.closed else "") + ("ambiguous" if sh.ambiguous else "")
        lines.append(f"sheet {n}: {len(sh.strands)} strand(s) {flags}".rstrip())
        parts = []
        for k, (s, d) in enumerate(zip(sh.strands, sh.directions)):
            parts.append(f"{s.id}({'+' if d > 0 else '-'})")
            if k < len(sh.pair_kinds):
                parts.append({"P": "=par=", "A": "=anti=", "?": "=?="}[sh.pair_kinds[k]])
        lines.append("  " + " ".join(parts))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="foldmap")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("summary", help="print secondary structure elements and sheet ordering")
    p.add_argument("structure")
    args = parser.parse_args(argv)
    try:
        print(summarize(args.structure))
    except (FileNotFoundError, ValueError) as err:
        print(f"foldmap: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

**`src/foldmap/__main__.py`**

```python
from .cli import main

raise SystemExit(main())
```

- [ ] **Step 4: Run the full suite**

Run: `python -m pytest -q`
Expected: 25 passed.

- [ ] **Step 5: Check the real structures by eye**

Run: `foldmap summary tests/data/1UBQ.cif && foldmap summary tests/data/6ZS5.cif`
Expected for ubiquitin: `helix  A:23-34` and one 5-strand sheet `A:12-16(+) =anti= A:2-7(-) =par= A:66-71(-) =anti= A:41-45(+) =anti= A:48-49(-)` (or its reverse). Expected for 6ZS5: 21 strands in 5 sheets; sheet 1 is `A:331-333, A:339-345, A:377-383, A:369-371`; at least one sheet mixes chains A and D (the linker completing the neighbouring sheet).

---

## Self-review

- **Spec coverage (step 1):** `io` (Task 1), `ss` with DSSP (Tasks 2-3), `sheets` with pairing graph, cycle cut and branch fallback (Task 4), CLI printing SSEs and sheet ordering (Task 5), unit tests on 6ZS5/6ZYA fixtures (Tasks 2-5). Deviations from the spec text: deposited annotation is a test cross-check rather than the primary source (reason in Architecture); first altloc instead of highest occupancy.
- **Placeholder scan:** none; every step has literal code or a literal command with expected output.
- **Type consistency:** `Backbone`, `Bridge`, `SSE`, `Sheet`, `DsspResult` names and fields are identical in every task; `build_sses(bb, ss)` and `build_sheets(sses, bridges)` signatures match their callers in `cli.py`.
- **Known limits carried to later plans:** DSSP has no bulge merging, so a strand interrupted by a beta-bulge can split into two strands; ordering inside ambiguous (branched) sheets is a heuristic and shows `?` between non-paired neighbours; a 20,000-residue complex has not been timed.
