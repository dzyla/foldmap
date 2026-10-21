"""Salt bridges and cation-π from the coordinates (foldmap.interactions): known contacts found, each one's atoms
meeting its rule, the limits honoured at their edges, selections as for H-bonds, output sorted and deterministic."""

import math
import os
import subprocess
import sys
from dataclasses import FrozenInstanceError
from functools import lru_cache
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.interactions import (
    CATION_PI_MAX,
    CATION_PI_OFFSET,
    HIS_NOTE,
    SALT_MAX,
    Interaction,
    find_cation_pi,
    find_salt_bridges,
)
from foldmap.io import read_model

DATA = Path(__file__).parent / "data"
LAMBDA = DATA / "1LMB.cif"

_CATIONIC = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2"}}
_ANIONIC = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}


@lru_cache(maxsize=None)
def _lambda():
    return read_model(LAMBDA)


def _atoms(model):
    """(chain, residue tag, atom) -> (residue name, gemmi.Position), the first residue of a number and the first
    conformer of each atom, as the finders read them."""
    out = {}
    for chain in model:
        for res in chain:
            for atom in res:
                key = (chain.name, f"{res.seqid.num}{res.seqid.icode.strip()}", atom.name)
                out.setdefault(key, (res.name, atom.pos))
    return out


def _sort_key(x):
    return (x.a, x.b, x.kind)


# λ repressor (chains 3 and 4) on its operator (DNA chains 1 and 2): lysines on the backbone phosphates, in both
# half-sites, as measured in the file (cation, anion, distance in Å).
LAMBDA_PHOSPHATE_SALT_BRIDGES = {
    (("3", "19", "NZ"), ("1", "4", "OP1"), 2.945),
    (("4", "19", "NZ"), ("2", "24", "OP1"), 2.640),
    (("3", "26", "NZ"), ("1", "3", "OP1"), 3.580),
    (("4", "26", "NZ"), ("2", "23", "OP1"), 3.304),
}


def test_lambda_repressor_lysines_bridge_the_dna_phosphates():
    found = find_salt_bridges(_lambda(), "na")
    assert {(x.a, x.b, x.distance) for x in found} >= LAMBDA_PHOSPHATE_SALT_BRIDGES


def test_every_lambda_salt_bridge_meets_the_rule():
    model, _ = _lambda()
    atoms = _atoms(model)
    for select in ("na", "interchain"):
        found = find_salt_bridges(_lambda(), select)
        assert found
        for x in found:
            assert x.kind == "salt" and x.angle is None
            (cname, cpos), (aname, apos) = atoms[x.a], atoms[x.b]
            assert x.a[2] in _CATIONIC.get(cname, set()) or x.a[2] == "N"
            assert x.b[2] in _ANIONIC.get(aname, set()) | {"OP1", "OP2", "O", "OXT"}
            assert x.distance == round(cpos.dist(apos), 3) <= SALT_MAX
            assert x.note == (HIS_NOTE if cname == "HIS" and x.a[2] != "N" else "")
            assert x.a[:2] != x.b[:2]
            if select == "na":
                assert x.b[2] in ("OP1", "OP2")


def test_na_keeps_protein_with_nucleic_acid_and_interchain_adds_protein_pairs():
    na = find_salt_bridges(_lambda(), "na")
    interchain = find_salt_bridges(_lambda(), "interchain")
    assert set(na) < set(interchain)
    # the dimer interface: Lys67 of each monomer on Glu89 of the other
    assert {(x.a, x.b) for x in interchain} >= {
        (("3", "67", "NZ"), ("4", "89", "OE1")),
        (("4", "67", "NZ"), ("3", "89", "OE1")),
    }


def test_chain_pairs_select_only_their_pairs():
    found = find_salt_bridges(_lambda(), "chains:2+4")
    assert found and all({x.a[0], x.b[0]} == {"2", "4"} for x in found)
    with pytest.raises(ValueError, match="no chain Z"):
        find_salt_bridges(_lambda(), "chains:3+Z")


def test_off_finds_nothing():
    assert find_salt_bridges(_lambda(), "off") == []
    assert find_cation_pi(_lambda(), "off") == []


def _model_of(*residues):
    """A model from (chain, seq, name, {atom: (x, y, z)}), each a polymer residue, in that order."""
    model = gemmi.Model("1")
    for chain, seq, name, atoms in residues:
        if chain not in [c.name for c in model]:
            model.add_chain(gemmi.Chain(chain))
        res = gemmi.Residue()
        res.name, res.seqid, res.entity_type = name, gemmi.SeqId(seq, " "), gemmi.EntityType.Polymer
        for atom_name, xyz in atoms.items():
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = atom_name, gemmi.Element(atom_name[0]), gemmi.Position(*xyz)
            res.add_atom(atom)
        model[chain].add_residue(res)
    return model


@pytest.mark.parametrize("dist, count", [(3.9, 1), (4.0, 1), (4.1, 0)])
def test_a_salt_bridge_is_within_4_angstroms(dist, count):
    model = _model_of(("A", 1, "LYS", {"NZ": (0, 0, 0)}), ("B", 1, "ASP", {"OD1": (dist, 0, 0)}))
    found = find_salt_bridges(model, "interchain")
    assert len(found) == count
    if count:
        assert found == [Interaction("salt", ("A", "1", "NZ"), ("B", "1", "OD1"), dist, None)]


def test_the_cation_is_a_whatever_the_model_order():
    model = _model_of(("A", 7, "GLU", {"OE2": (0, 0, 0)}), ("B", 3, "ARG", {"NH1": (3.0, 0, 0)}))
    assert find_salt_bridges(model, "interchain") == [
        Interaction("salt", ("B", "3", "NH1"), ("A", "7", "OE2"), 3.0, None)
    ]


def test_a_histidine_bridge_carries_its_note():
    model = _model_of(("A", 1, "HIS", {"NE2": (0, 0, 0)}), ("B", 1, "GLU", {"OE1": (3.0, 0, 0)}))
    [found] = find_salt_bridges(model, "interchain")
    assert found.note == HIS_NOTE == "His if protonated"


def test_two_like_charges_and_neutral_atoms_make_no_bridge():
    lys = ("A", 1, "LYS", {"NZ": (0, 0, 0)})
    assert find_salt_bridges(_model_of(lys, ("B", 1, "ARG", {"NH1": (3.0, 0, 0)})), "interchain") == []
    assert find_salt_bridges(_model_of(lys, ("B", 1, "ASN", {"OD1": (3.0, 0, 0)})), "interchain") == []
    asp = ("A", 1, "ASP", {"OD1": (0, 0, 0)})
    assert find_salt_bridges(_model_of(asp, ("B", 1, "GLU", {"OE1": (3.0, 0, 0)})), "interchain") == []


def test_chain_termini_are_charged_but_a_gap_is_not():
    """Chain A runs 1-3 with a gap after 1 (no peptide bond to 2): only residue 1's N and residue 3's O count."""
    backbone = {"N": (0, 0, 0), "CA": (1.46, 0, 0), "C": (2.0, 1.4, 0), "O": (1.4, 2.4, 0)}

    def shifted(dx):
        return {k: (x + dx, y, z) for k, (x, y, z) in backbone.items()}

    chain_a = [("A", 1, "GLY", shifted(0)), ("A", 2, "GLY", shifted(20)), ("A", 3, "GLY", shifted(40))]
    partners = [
        ("B", 1, "ASP", {"OD1": (-3.0, 0, 0)}),  # by A1 N (the N-terminus)
        ("B", 2, "LYS", {"NZ": (1.4, 5.4, 0)}),  # by A1 O (before the gap)
        ("B", 3, "ASP", {"OD1": (17.0, 0, 0)}),  # by A2 N (after the gap)
        ("B", 4, "LYS", {"NZ": (41.4, 5.4, 0)}),  # by A3 O (the C-terminus)
    ]
    found = find_salt_bridges(_model_of(*chain_a, *partners), "interchain")
    assert [(x.a, x.b) for x in found] == [
        (("A", "1", "N"), ("B", "1", "OD1")),
        (("B", "4", "NZ"), ("A", "3", "O")),
    ]


def test_an_oxt_marks_a_c_terminus():
    model = _model_of(
        (
            "A",
            1,
            "GLY",
            {"N": (0, 0, 0), "CA": (1.46, 0, 0), "C": (2.0, 1.4, 0), "O": (1.4, 2.4, 0), "OXT": (3.2, 1.4, 0)},
        ),
        ("A", 2, "GLY", {"N": (30, 0, 0), "CA": (31.46, 0, 0)}),
        ("B", 1, "LYS", {"NZ": (6.0, 1.4, 0)}),
    )
    [found] = find_salt_bridges(model, "interchain")
    assert (found.a, found.b) == (("B", "1", "NZ"), ("A", "1", "OXT"))


def _benzene(radius=1.39):
    """Phe ring atoms, a regular hexagon in the z = 0 plane about the origin."""
    names = ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")
    return {
        n: (radius * math.cos(k * math.pi / 3), radius * math.sin(k * math.pi / 3), 0.0) for k, n in enumerate(names)
    }


@pytest.mark.parametrize("offset, count", [(1.9, 1), (2.1, 0)])
def test_cation_pi_needs_the_cation_over_the_ring(offset, count):
    model = _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (offset, 0.0, 3.6)}))
    found = find_cation_pi(model, "interchain")
    assert len(found) == count
    if count:
        [x] = found
        assert (x.kind, x.a, x.b) == ("cation-pi", ("B", "1", "NZ"), ("A", "1", "ring:6"))
        assert x.distance == pytest.approx(math.hypot(offset, 3.6), abs=1e-3)
        assert x.angle == pytest.approx(math.degrees(math.atan2(offset, 3.6)), abs=1e-3)


@pytest.mark.parametrize("height, count", [(5.9, 1), (6.1, 0)])
def test_cation_pi_is_within_6_angstroms(height, count):
    model = _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (0.0, 0.0, -height)}))
    found = find_cation_pi(model, "interchain")
    assert len(found) == count
    if count:
        assert found[0].distance == pytest.approx(height, abs=1e-3) and found[0].angle == pytest.approx(0.0, abs=1e-3)


def test_the_limits_are_the_stated_ones():
    assert (SALT_MAX, CATION_PI_MAX, CATION_PI_OFFSET) == (4.0, 6.0, 2.0)


def test_the_arginine_cation_is_its_guanidinium_centre():
    """NE, NH1 and NH2 round a point 3.8 Å over the ring, 1.0 Å off its axis."""
    centre = np.array([1.0, 0.0, 3.8])
    arg = {
        n: tuple(centre + 1.3 * np.array([math.cos(t), 0.0, math.sin(t)]))
        for n, t in (("NE", 0.0), ("NH1", 2.0944), ("NH2", 4.1888))
    }
    model = _model_of(("A", 1, "TYR", _benzene()), ("B", 1, "ARG", arg))
    [x] = find_cation_pi(model, "interchain")
    assert (x.a, x.b) == (("B", "1", "centre"), ("A", "1", "ring:6"))
    assert x.distance == pytest.approx(math.hypot(1.0, 3.8), abs=2e-3)


def test_tryptophan_gives_both_rings_and_a_ring_needs_all_its_atoms():
    trp = {
        "CG": (-1.1, 0.7, 0),
        "CD1": (-2.0, -0.3, 0),
        "NE1": (-1.3, -1.4, 0),
        "CE2": (0.0, -1.0, 0),
        "CD2": (0.1, 0.4, 0),
        "CZ2": (1.1, -1.9, 0),
        "CH2": (2.4, -1.4, 0),
        "CZ3": (2.6, 0.0, 0),
        "CE3": (1.4, 0.9, 0),
    }
    lys = ("B", 1, "LYS", {"NZ": (0.3, -0.3, 3.5)})
    found = find_cation_pi(_model_of(("A", 1, "TRP", trp), lys), "interchain")
    assert [x.b for x in found] == [("A", "1", "ring:5"), ("A", "1", "ring:6")]
    partial = {k: v for k, v in trp.items() if k != "CZ2"}
    found = find_cation_pi(_model_of(("A", 1, "TRP", partial), lys), "interchain")
    assert [x.b for x in found] == [("A", "1", "ring:5")]


def test_the_lambda_cation_pi_contacts_meet_the_rule():
    model, _ = _lambda()
    atoms = _atoms(model)
    found = find_cation_pi(_lambda(), "interchain")
    assert found
    for x in found:
        assert x.kind == "cation-pi" and x.a[2] in ("NZ", "centre") and x.b[2] in ("ring:5", "ring:6")
        assert x.distance <= CATION_PI_MAX and 0.0 <= x.angle <= 90.0
        offset = x.distance * math.sin(math.radians(x.angle))
        assert offset <= CATION_PI_OFFSET + 1e-3
        if x.a[2] == "NZ":
            assert atoms[x.a][0] == "LYS"


def test_output_is_sorted_by_a_then_b_then_kind():
    for finder in (find_salt_bridges, find_cation_pi):
        for select in ("na", "interchain"):
            found = finder(_lambda(), select)
            assert found == sorted(found, key=_sort_key)
            assert len(set(found)) == len(found)


def test_an_interaction_is_frozen_and_names_its_kind():
    x = Interaction("salt", ("A", "1", "NZ"), ("B", "1", "OD1"), 3.0, None)
    with pytest.raises(FrozenInstanceError):
        x.distance = 2.0
    with pytest.raises(ValueError, match="cation-pi"):
        Interaction("ionic", ("A", "1", "NZ"), ("B", "1", "OD1"), 3.0, None)


_PRINT = """
import sys
from foldmap.interactions import find_cation_pi, find_salt_bridges
from foldmap.io import read_model
model = read_model(sys.argv[1])
print(repr(find_salt_bridges(model, sys.argv[2])), repr(find_cation_pi(model, sys.argv[2])))
"""


def test_output_does_not_depend_on_the_hash_seed():
    runs = set()
    for seed in ("0", "12345", "777"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run(
            [sys.executable, "-c", _PRINT, str(LAMBDA), "interchain"],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1 and "Interaction(kind='salt'" in (out := runs.pop()) and "kind='cation-pi'" in out
