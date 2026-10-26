"""Metamorphic tests (plan 28 Task 7, review 3.6; ruling Q34): rewrite a committed fixture under a transform that
should not change what foldmap reads from it, and hold foldmap to the level each invariant should hold at.

- Exact transforms (a rigid motion that is exact once written, chain renaming that keeps the chains' name order,
  the atoms of each residue shuffled, mmCIF rewritten as PDB): the structure facts (DSSP codes, bridges, elements,
  sheet order, directions, pair kinds, closure and accessory strands, base pairs, the transmembrane set, automatic
  domains) and `foldmap check`'s verdicts are equal, and the figure's geometry is equal to the last bit. The page frame
  is derived from the structure (frame.view_frame: eigh and svd of the elements), so no 2D symmetry is needed: the
  measured gap is 0.0 on every fixture here.
- A general rotation (37 degrees about (1, 2, 3)) is not exact once written (coordinates rounded to 0.001 Å): the
  structure facts are asserted equal (measured equal), the geometry within GENERAL_GAP.
- A renaming that reverses the chains' name order: the facts and verdicts are equal; the figure of a symmetric
  assembly of several deposited chains is not (finding F1, pinned xfail(strict): the C2 dimers 1TIM and 1LMB, the C4
  tetramer 1BL8).
- Coordinate noise (sigma 0.05 and 0.1 Å on every coordinate, three seeds): measured as a change rate and asserted
  bounded (NOISE), not equal.

Every baseline is the fixture rewritten by the same writer with no transform (metamorphic.rewrite), so what the writer
leaves out of a file (categories gemmi does not keep) is not taken for the transform's effect; the round trip itself
is held to the original file's facts (test_round_trip_keeps_the_facts).

The findings and the measurements behind the limits here are recorded in docs/validation/metamorphic/report.md
(scripts/metamorphic_record.py)."""

from __future__ import annotations

import functools

import gemmi
import numpy as np
import pytest

import metamorphic as M

# fixture -> what is read from it: nucleic (base pairs), membrane (transmembrane set; the figure drawn --membrane auto)
FIXTURES = {
    "1UBQ.cif": {},
    "1TIM.cif": {},
    "1C3W.cif": {"membrane": True},
    "1LMB.cif": {"nucleic": True},
    "na/1EHZ.cif.gz": {"nucleic": True},
}
EXACT = ("motion", "rename", "shuffle", "pdb")
GENERAL_GAP = 0.01  # page units (and cos, sin of an angle): a general rotation, written to 0.001 Å (measured 6e-4)
# noise: (sigma Å) -> the most each may change (DSSP codes: the fraction of residues; bridges, base pairs and
# transmembrane residues: |x △ y| / (|x| + |y|), so one bridge of 1C3W's 6 moving is not the whole limit); the
# measured maxima over the fixtures and 3 seeds are in docs/validation/metamorphic/report.md
NOISE = {
    0.05: {"ss": 0.05, "bridges": 0.25, "base_pairs": 0.10, "tm_residues": 0.05},
    0.1: {"ss": 0.08, "bridges": 0.25, "base_pairs": 0.10, "tm_residues": 0.05},
}
F1 = (
    "T7 finding F1: protomer order follows chain names (symmetry._orbits iterates sorted(perm), symmetry.py:234; "
    "symmetry._protomers sorts orbits by first chain name, symmetry.py:255, and breaks its ties by name, "
    "symmetry.py:259); protomer 0 is laid out first and drawn "
    "on the left (cli._dyad_frame, cli.py:359) or is the copy the other mirrors (cli._mirror_frame, cli.py:257)"
)
XF1 = pytest.mark.xfail(strict=True, raises=AssertionError, reason=F1)


def _chains(name: str) -> list[str]:
    return sorted({c.name for c in gemmi.read_structure(str(M.DATA / name))[0]})


def _mapping(name: str, order: str) -> dict[str, str]:
    """A renaming of the fixture's chains: "keep" keeps their name order (A, B -> K, L), "reverse" reverses it."""
    chains = _chains(name)
    if order == "keep":
        return {c: chr(ord("K") + k) for k, c in enumerate(chains)}
    return {c: chr(ord("Z") - k) for k, c in enumerate(chains)}


@pytest.fixture(scope="module")
def work(tmp_path_factory):
    return tmp_path_factory.mktemp("metamorphic")


def _file(work, name: str, transform: str):
    """The fixture rewritten under `transform`, and the renaming to undo (new name -> old)."""
    src = M.DATA / name
    stem = name.replace("/", "_").split(".")[0]
    out = work / f"{stem}-{transform}.{'pdb' if transform == 'pdb' else 'cif'}"
    back: dict[str, str] = {}
    kw: dict = {}
    if transform == "motion":
        kw = {"rot": M.CYCLE, "shift": M.SHIFT}
    elif transform == "turn":
        kw = {"rot": M.turn([1, 2, 3], 37.0), "shift": (5.5, -3.25, 8.0)}
    elif transform in ("rename", "reverse"):
        mapping = _mapping(name, "keep" if transform == "rename" else "reverse")
        kw = {"rename": mapping}
        back = {v: k for k, v in mapping.items()}
    elif transform == "shuffle":
        kw = {"shuffle": 7}
    if not out.is_file():
        M.rewrite(src, out, **kw)
    return out, back


@functools.cache
def _read(work, name: str, transform: str, figure: bool) -> dict:
    """What foldmap reads (facts) and, with `figure`, draws (geometry) and judges (verdicts), renamed back."""
    from foldmap.spec import draw_spec, figure_spec

    path, back = _file(work, name, transform)
    kw = FIXTURES.get(name, {})
    out = {"facts": M.facts(path, **kw)}
    if figure:
        spec = figure_spec(str(path), membrane="auto" if kw.get("membrane") else None)
        drawn = draw_spec(spec)
        out["geometry"] = M.geometry(drawn)
        out["verdicts"] = M.verdicts(drawn, spec, path)
        import matplotlib.pyplot as plt

        plt.close(drawn.figure)
    return M.renamed(out, back) if back else out


def _gap(a: dict, b: dict) -> dict[str, float]:
    return {k: M.geometry_gap({k: a[k]}, {k: b[k]}) for k in a}


@pytest.mark.parametrize("name", FIXTURES)
def test_round_trip_keeps_the_facts(work, name):
    """The writer's round trip (the baseline of every transform) reads as the original file does."""
    assert _read(work, name, "base", False)["facts"] == M.facts(M.DATA / name, **FIXTURES[name])


@pytest.mark.parametrize("name", FIXTURES)
def test_the_shuffle_reorders_atoms(work, name):
    """The shuffled file really holds some residue's atoms in another order than the baseline (else "shuffle" would
    test nothing), and the same atoms."""
    a, b = (gemmi.read_structure(str(_file(work, name, t)[0])) for t in ("base", "shuffle"))
    names = [([x.name for x in r], [y.name for y in s]) for ca, cb in zip(a[0], b[0]) for r, s in zip(ca, cb)]
    assert all(sorted(x) == sorted(y) for x, y in names)
    several = [(x, y) for x, y in names if len(x) >= 3]  # a water's one atom has one order
    assert sum(x != y for x, y in several) > len(several) // 2


@pytest.mark.parametrize("transform", EXACT)
@pytest.mark.parametrize("name", FIXTURES)
def test_structure_facts_are_exact(work, name, transform):
    assert _read(work, name, transform, False)["facts"] == _read(work, name, "base", False)["facts"]


@pytest.mark.parametrize("transform", EXACT)
@pytest.mark.parametrize("name", FIXTURES)
def test_check_verdicts_and_figure_are_exact(work, name, transform):
    """`foldmap check`'s verdicts equal, and the figure equal to the last bit (an assembly read from PDB aside:
    test_pdb_assembly_operators_are_rounded)."""
    if transform == "pdb" and name == "1C3W.cif":
        pytest.skip("its assembly is built from PDB BIOMT operators written to 6 decimals: measured apart")
    got, base = _read(work, name, transform, True), _read(work, name, "base", True)
    assert got["verdicts"] == base["verdicts"]
    gaps = _gap(base["geometry"], got["geometry"])
    assert max(gaps.values()) == 0.0, gaps


def test_pdb_assembly_operators_are_rounded(work):
    """1C3W's trimer, built from the PDB file's REMARK 350 operators (written to 6 decimals), lies 2e-5 Å from the
    mmCIF one: its figure moves by up to 0.014 page units (0.0024 with --membrane auto); drawn as deposited
    (--assembly asu), it is equal to the last bit."""
    from foldmap.io import load_backbone
    from foldmap.spec import draw_spec, figure_spec

    cif, _ = _file(work, "1C3W.cif", "base")
    pdb, _ = _file(work, "1C3W.cif", "pdb")
    a, b = load_backbone(cif, "auto"), load_backbone(pdb, "auto")
    assert 0 < np.abs(a.xyz - b.xyz).max() < 1e-4
    for assembly, limit in (("asu", 0.0), ("auto", 0.05)):
        g = [M.geometry(draw_spec(figure_spec(str(p), membrane="auto", assembly=assembly))) for p in (cif, pdb)]
        assert M.geometry_gap(*g) <= limit, (assembly, _gap(*g))


@pytest.mark.parametrize("name", FIXTURES)
def test_general_rotation(work, name):
    got, base = _read(work, name, "turn", True), _read(work, name, "base", True)
    assert got["facts"] == base["facts"]
    assert got["verdicts"] == base["verdicts"]
    gaps = _gap(base["geometry"], got["geometry"])
    assert max(gaps.values()) <= GENERAL_GAP, gaps


@pytest.mark.parametrize(
    "name",
    [
        pytest.param("1TIM.cif", marks=XF1),  # sides swap: 28.9
        pytest.param("1LMB.cif", marks=XF1),  # a loop by 0.75
        pytest.param("1BL8.cif", marks=XF1),  # C4: 15.8
    ],
)
def test_reversed_chain_names_draw_the_same_figure(work, name):
    got, base = _read(work, name, "reverse", True), _read(work, name, "base", True)
    assert got["facts"] == base["facts"] and got["verdicts"] == base["verdicts"]  # these hold: F1 is the figure's
    gaps = _gap(base["geometry"], got["geometry"])
    assert max(gaps.values()) == 0.0, gaps


@pytest.mark.parametrize("name", ["1TIM.cif", "1LMB.cif", "1BL8.cif"])
def test_reversed_chain_names_keep_facts_and_verdicts(work, name):
    """What F1 leaves standing: the structure facts and `foldmap check`'s verdicts of the renamed entry."""
    got, base = _read(work, name, "reverse", True), _read(work, name, "base", True)
    assert got["facts"] == base["facts"]
    assert got["verdicts"] == base["verdicts"]


@pytest.mark.parametrize("transform", EXACT)
def test_automatic_domains_are_exact(work, transform):
    """T4 lysozyme's two automatic domains (features.auto_domains) under each exact transform."""
    from foldmap.cli import make_layout

    path, back = _file(work, "2LZM.cif", transform)
    lay, _, _ = make_layout(path, domains="auto")
    domains = M.renamed([[n, list(ids)] for n, ids in lay.domains], back)
    assert [n for n, _ in domains] == ["D1", "D2"]
    assert domains == _domains_of(work)


@functools.cache
def _domains_of(work):
    from foldmap.cli import make_layout

    lay, _, _ = make_layout(_file(work, "2LZM.cif", "base")[0], domains="auto")
    return [[n, list(ids)] for n, ids in lay.domains]


def _rate(a: list, b: list) -> float:
    x, y = {tuple(v) if isinstance(v, list) else v for v in a}, {tuple(v) if isinstance(v, list) else v for v in b}
    return len(x ^ y) / max(len(x) + len(y), 1)


@pytest.mark.parametrize("sigma", sorted(NOISE))
@pytest.mark.parametrize("name", ["1TIM.cif", "1C3W.cif", "na/1EHZ.cif.gz"])
def test_noise_changes_are_bounded(work, name, sigma):
    """Gaussian noise on every coordinate: the fraction of DSSP codes, bridges, base pairs and transmembrane
    residues that change, over three seeds, stays within NOISE (measured, not equal: a residue on a definition's
    edge may flip)."""
    base = _read(work, name, "base", False)["facts"]
    limits = NOISE[sigma]
    for seed in range(3):
        path = M.rewrite(M.DATA / name, work / f"noise-{seed}.cif", noise=(sigma, seed))
        got = M.facts(path, **FIXTURES[name])
        assert got["residues"] == base["residues"]
        rates = {
            "ss": float(np.mean([a != b for a, b in zip(base["ss"], got["ss"])])) if base["ss"] else 0.0,
            "bridges": _rate(base["bridges"], got["bridges"]),
        }
        if "base_pairs" in base:
            rates["base_pairs"] = _rate(base["base_pairs"], got["base_pairs"])
        if "transmembrane" in base:  # by residue: a helix end moving by one renames its element
            rates["tm_residues"] = _rate(base["transmembrane_residues"], got["transmembrane_residues"])
        assert all(rates[k] <= limits[k] for k in rates), (seed, rates)
