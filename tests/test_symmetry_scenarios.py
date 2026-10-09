"""Symmetry detection on assemblies built with known operators from real chains, then on real structures,
then with the user's own settings."""

from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.symmetry import detect_symmetry
from helpers import pipeline

DATA = Path(__file__).parent / "data"
NAMES = "ABCDEFGHIJKLMNOP"


def rot(axis, deg):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    t = np.radians(deg)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * k @ k


def build(
    tmp_path,
    ops,
    src="1UBQ",
    offset=(22.0, 4.0, 3.0),
    names=None,
    renumber=None,
    drop=None,
    noise=0.0,
    extra=None,
    fname="asm.pdb",
):
    """One copy of `src` chain A per operator (R, t), after moving it `offset` away from the origin."""
    st = gemmi.read_structure(str(DATA / f"{src}.cif"))
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    base = st[0][0]
    xyz0 = np.array([a.pos.tolist() for r in base for a in r])
    shift = np.asarray(offset) - xyz0.mean(axis=0)
    out = gemmi.Structure()
    out.add_model(gemmi.Model("1"))
    rng = np.random.default_rng(7)
    for k, (r, t) in enumerate(ops):
        ch = gemmi.Chain((names or NAMES)[k])
        for i, res in enumerate(base):
            if drop and k in drop and i in drop[k]:
                continue
            new = res.clone()
            if renumber and k in renumber:
                new.seqid.num += renumber[k]
            for a in new:
                p = r @ (np.array(a.pos.tolist()) + shift) + np.asarray(t) + rng.normal(0, noise, 3)
                a.pos = gemmi.Position(*p)
            ch.add_residue(new)
        out[0].add_chain(ch)
    if extra is not None:
        other = gemmi.read_structure(str(DATA / f"{extra}.cif"))[0][0]
        ch = gemmi.Chain("X")
        for res in other:
            new = res.clone()
            for a in new:
                a.pos = gemmi.Position(*(np.array(a.pos.tolist()) + [0, 0, 80]))
            ch.add_residue(new)
        out[0].add_chain(ch)
    out.setup_entities()
    path = tmp_path / fname
    out.write_pdb(str(path))
    return path


def cyclic(n, axis=(0, 0, 1)):
    return [(rot(axis, 360 * k / n), np.zeros(3)) for k in range(n)]


def found(path, **kw):
    bb, sses, _ = pipeline(path)
    return detect_symmetry(bb, sses, **kw)


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_exact_cyclic_assemblies(tmp_path, n):
    s = found(build(tmp_path, cyclic(n)))
    assert s.kind == "C" and s.n == n and s.rmsd < 0.01
    assert abs(abs(s.axis @ [0, 0, 1]) - 1) < 1e-6
    assert sorted(c for p in s.protomers for c in p) == list(NAMES[:n])


def test_axis_need_not_be_z_or_through_the_origin(tmp_path):
    ops = [(r, np.array([5.0, -7.0, 11.0]) - r @ [5.0, -7.0, 11.0]) for r, _ in cyclic(3, axis=(1, 2, 2))]
    s = found(build(tmp_path, ops))
    assert s.n == 3 and abs(abs(s.axis @ (np.array([1, 2, 2]) / 3)) - 1) < 1e-6


def test_xray_like_noise_is_tolerated(tmp_path):
    s = found(build(tmp_path, cyclic(4), noise=0.6))
    assert s.n == 4 and 0.3 < s.rmsd < 3.0


def test_renumbered_and_trimmed_copies_still_match(tmp_path):
    s = found(build(tmp_path, cyclic(3), renumber={1: 1000, 2: -40}, drop={2: range(30, 42)}))
    assert s.n == 3


def test_chain_names_do_not_matter(tmp_path):
    s = found(build(tmp_path, cyclic(4), names="QZBK"))
    assert s.n == 4 and sorted(c for p in s.protomers for c in p) == sorted("QZBK")


@pytest.mark.parametrize("n", [2, 3])
def test_dihedral_assemblies(tmp_path, n):
    ops = [(rot((1, 0, 0), 180 * f) @ rot((0, 0, 1), 360 * k / n), np.zeros(3)) for f in (0, 1) for k in range(n)]
    s = found(build(tmp_path, ops, offset=(20.0, 6.0, 9.0)))
    assert s.kind == "D" and s.n == n
    assert len(s.protomers) == n and all(len(p) == 2 for p in s.protomers)  # each protomer holds one copy of each ring


def test_helical_filament(tmp_path):
    ops = [(rot((0, 0, 1), 33 * k), np.array([0, 0, 9.0 * k])) for k in range(6)]
    s = found(build(tmp_path, ops, offset=(15.0, 0.0, 0.0)))
    assert s.kind == "H" and len(s.protomers) == 6
    assert s.twist == pytest.approx(np.radians(33), abs=np.radians(2)) and s.rise == pytest.approx(9.0, abs=0.5)
    assert [p[0] for p in s.protomers] == list(NAMES[:6])  # in order along the filament, as built


def test_translated_copies_are_not_symmetry(tmp_path):
    ops = [(np.eye(3), np.array([60.0 * k, 0, 0])) for k in range(3)]  # lattice-like neighbours, no rotation
    assert found(build(tmp_path, ops)) is None


def test_copies_that_differ_too_much_are_not_symmetry_unless_the_tolerance_allows(tmp_path):
    path = build(tmp_path, cyclic(3), noise=1.5)  # ~3.7 Å RMSD between copies
    assert found(path) is None
    assert found(path, tolerance=6.0).n == 3


def test_an_extra_single_copy_chain_is_left_out(tmp_path):
    s = found(build(tmp_path, cyclic(3), extra="2LZM"))
    assert s.n == 3 and "X" in s.others and all("X" not in p for p in s.protomers)


def test_monomer_has_no_symmetry(tmp_path):
    assert found(build(tmp_path, cyclic(1))) is None


# --- real structures --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "name, kind, n, per",
    [("1BL8", "C", 4, 1), ("8UUP", "C", 3, 2), ("1TIM", "C", 2, 1), ("2HHB", "C", 2, 2), ("1LMB", "C", 2, 1)],
)
def test_real_assemblies(name, kind, n, per):
    s = found(DATA / f"{name}.cif")
    assert (s.kind, s.n) == (kind, n) and all(len(p) == per for p in s.protomers) and s.rmsd < 3.0


@pytest.mark.parametrize("name", ["1UBQ", "5NKT", "1ZAA", "6ZS5"])
def test_real_asymmetric(name):
    assert found(DATA / f"{name}.cif") is None


# --- the user's say -------------------------------------------------------------------------------------------
def test_user_can_ask_for_a_subgroup():
    s = found(DATA / "1BL8.cif", request="C2")  # a C4 tetramer also has a 2-fold: two protomers of two chains
    assert (s.kind, s.n) == ("C", 2) and all(len(p) == 2 for p in s.protomers)


def test_asking_for_a_symmetry_that_is_not_there_is_an_error():
    with pytest.raises(ValueError, match="C3"):
        found(DATA / "1TIM.cif", request="C3")


def test_user_protomers_override_the_grouping():
    s = found(DATA / "2HHB.cif", protomers=[["A", "D"], ["C", "B"]])
    assert s.n == 2 and [sorted(p) for p in s.protomers] == [["A", "D"], ["B", "C"]]


def test_user_protomers_must_be_real_chains():
    with pytest.raises(ValueError, match="Q"):
        found(DATA / "2HHB.cif", protomers=[["A", "Q"], ["C", "D"]])


def test_symmetry_off():
    assert found(DATA / "1BL8.cif", request="off") is None


def test_cli_symmetry_options(tmp_path, capsys):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(DATA / "1BL8.cif"), "-o", str(out), "--symmetry", "C2"]) == 0
    assert main(["plot", str(DATA / "2HHB.cif"), "-o", str(out), "--protomers", "A,B;C,D", "--symmetry-tol", "2"]) == 0
    assert main(["plot", str(DATA / "1TIM.cif"), "-o", str(out), "--symmetry", "C5"]) == 1
    assert "C5" in capsys.readouterr().err
    assert main(["summary", str(DATA / "8UUP.cif")]) == 0
    assert "symmetry: C3" in capsys.readouterr().out


def test_symmetric_layout_works_for_synthetic_c4_and_helix(tmp_path):
    from foldmap.cli import make_layout

    for fname, ops in (
        ("c4.pdb", cyclic(4)),
        ("h.pdb", [(rot((0, 0, 1), 33 * k), np.array([0, 0, 9.0 * k])) for k in range(5)]),
    ):
        path = build(tmp_path, ops, offset=(18.0, 0.0, 0.0), fname=fname)
        lay, sses, _ = make_layout(path)
        first = {s.chain: s for s in sses if s.kind == "H"}
        ys = [lay.placed[first[c].id].cy for c in sorted(first)]
        xs = [lay.placed[first[c].id].cx for c in sorted(first)]
        assert len(set(np.round(np.diff(xs), 3))) == 1  # each copy one step further along
        if fname == "c4.pdb":
            assert np.ptp(ys) < 1e-6  # a ring: all copies level


def test_a_large_ring_is_found_quickly(tmp_path):
    import time

    path = build(tmp_path, cyclic(12), offset=(45.0, 0.0, 0.0))
    t = time.perf_counter()
    s = found(path)
    assert s.n == 12 and time.perf_counter() - t < 10
