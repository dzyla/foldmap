"""Salt bridges, cation-π, π-stacking and hydrophobic contacts from the coordinates (foldmap.interactions): known
contacts found, each one's atoms meeting its rule, the limits honoured at their edges, selections as for H-bonds,
output sorted and deterministic; the independent re-check (verify_interactions) passes every found line and flags a
moved or false one."""

import dataclasses
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

from foldmap.hbonds import find_hbonds
from foldmap.interactions import (
    CATION_PI_MAX,
    CATION_PI_OFFSET,
    HIS_NOTE,
    HYDROPHOBIC_MAX,
    SALT_MAX,
    STACK_MAX,
    STACK_OFFSET,
    STACK_PARALLEL,
    STACK_T_MIN,
    Interaction,
    find_cation_pi,
    find_hydrophobic,
    find_pi_stacking,
    find_salt_bridges,
    sort_key,
    verify_interactions,
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
    """a, then b, each (chain, residue number as a number, insertion code, atom), then kind."""

    def partner(p):
        digits = p[1].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
        return (p[0], int(digits), p[1][len(digits) :], p[2])

    return (partner(x.a), partner(x.b), x.kind)


_UP, _DOWN = math.nextafter(4.0, 5.0), math.nextafter(4.0, 3.0)
_FINDERS = (find_salt_bridges, find_cation_pi, find_pi_stacking, find_hydrophobic)


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
    for finder in _FINDERS:
        assert finder(_lambda(), "off") == []


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


@pytest.mark.parametrize(
    "dist, count",
    [
        (3.9, 1),
        (4.0, 1),
        (math.nextafter(4.0, 5.0), 1),
        (math.nextafter(4.0, 3.0), 1),
        (4.0004, 1),
        (4.0006, 0),
        (4.1, 0),
    ],
)
def test_a_salt_bridge_is_within_4_angstroms_on_the_grid(dist, count):
    """Distances are snapped to 1e-3 Å before the limit sees them: 1 ulp either side of 4.0 is 4.0, and so is
    4.0004; 4.0006 is 4.001."""
    model = _model_of(("A", 1, "LYS", {"NZ": (0, 0, 0)}), ("B", 1, "ASP", {"OD1": (dist, 0, 0)}))
    found = find_salt_bridges(model, "interchain")
    assert len(found) == count
    if count:
        assert found == [Interaction("salt", ("A", "1", "NZ"), ("B", "1", "OD1"), round(dist, 3), None)]


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


_EDGE_OFFSETS = [
    (1.9, 1),
    (2.0, 1),
    (math.nextafter(2.0, 3.0), 1),
    (math.nextafter(2.0, 1.0), 1),
    (2.0004, 1),
    (2.0006, 0),
    (2.1, 0),
]


@pytest.mark.parametrize("offset, count", _EDGE_OFFSETS)
def test_cation_pi_needs_the_cation_over_the_ring(offset, count):
    model = _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (offset, 0.0, 3.6)}))
    found = find_cation_pi(model, "interchain")
    assert len(found) == count
    if count:
        [x] = found
        assert (x.kind, x.a, x.b) == ("cation-pi", ("B", "1", "NZ"), ("A", "1", "ring:6"))
        assert x.distance == pytest.approx(math.hypot(offset, 3.6), abs=1e-3)
        assert x.angle == pytest.approx(math.degrees(math.atan2(offset, 3.6)), abs=1e-3)


_EDGE_HEIGHTS = [
    (5.9, 1),
    (6.0, 1),
    (math.nextafter(6.0, 7.0), 1),
    (math.nextafter(6.0, 5.0), 1),
    (6.0004, 1),
    (6.0006, 0),
    (6.1, 0),
]


@pytest.mark.parametrize("height, count", _EDGE_HEIGHTS)
def test_cation_pi_is_within_6_angstroms(height, count):
    model = _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (0.0, 0.0, -height)}))
    found = find_cation_pi(model, "interchain")
    assert len(found) == count
    if count:
        assert found[0].distance == round(height, 3) and found[0].angle == 0.0


def test_a_contact_at_the_limits_reports_the_limits():
    """Offset exactly 2.0 Å with the cation 6.0 Å from the centroid: both on their limits, both kept, and the
    reported distance and angle are on the 1e-3 grid."""
    height = math.sqrt(6.0**2 - 2.0**2)
    model = _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (2.0, 0.0, height)}))
    [x] = find_cation_pi(model, "interchain")
    assert x.distance == 6.0
    assert x.angle == round(math.degrees(math.asin(2.0 / 6.0)), 3)


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


def test_residue_numbers_sort_as_numbers():
    model = _model_of(
        ("A", 19, "LYS", {"NZ": (0, 0, 0)}),
        ("A", 5, "LYS", {"NZ": (20, 0, 0)}),
        ("B", 1, "ASP", {"OD1": (3.0, 0, 0)}),
        ("B", 2, "ASP", {"OD1": (23.0, 0, 0)}),
    )
    found = find_salt_bridges(model, "interchain")
    assert [x.a[1] for x in found] == ["5", "19"]
    tags = ["100A", "-7", "19", "100", "5", "101"]
    keys = sorted(tags, key=lambda t: sort_key(Interaction("salt", ("A", t, "NZ"), ("B", "1", "OD1"), 3.0, None)))
    assert keys == ["-7", "5", "19", "100", "100A", "101"]


def test_output_is_sorted_by_a_then_b_then_kind():
    for finder in _FINDERS:
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
from foldmap.interactions import find_cation_pi, find_hydrophobic, find_pi_stacking, find_salt_bridges
from foldmap.io import read_model
model = read_model(sys.argv[1])
for finder in (find_salt_bridges, find_cation_pi, find_pi_stacking, find_hydrophobic):
    print(repr(finder(model, sys.argv[2])))
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
    assert len(runs) == 1
    out = runs.pop()
    assert all(f"kind={k!r}" in out for k in ("salt", "cation-pi", "hydrophobic"))


# π-stacking


def _ring_at(angle, shift, names=("CG", "CD1", "CE1", "CZ", "CE2", "CD2")):
    """A regular hexagon (the Phe ring, or a base's six-membered ring under `names`) turned `angle` degrees about
    the y axis, then moved by `shift`: its normal makes `angle` with the z axis."""
    t = math.radians(angle)
    out = {}
    for n, (x, y, z) in zip(names, _benzene().values()):
        out[n] = (
            x * math.cos(t) + z * math.sin(t) + shift[0],
            y + shift[1],
            -x * math.sin(t) + z * math.cos(t) + shift[2],
        )
    return out


def _two_phe(angle, shift):
    return _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "PHE", _ring_at(angle, shift)))


def test_two_phenylalanines_stack_face_to_face():
    [x] = find_pi_stacking(_two_phe(0.0, (1.0, 0.0, 3.6)), "interchain")
    assert x == Interaction(
        "stack", ("A", "1", "ring:6"), ("B", "1", "ring:6"), round(math.hypot(1.0, 3.6), 3), 0.0, "parallel"
    )


def test_two_phenylalanines_stack_edge_to_face():
    """The second ring stands on the first, its centroid 5.0 Å over the first one's: the offset is the smaller of
    the two projections (0 here, over the face ring), as PLIP takes it."""
    [x] = find_pi_stacking(_two_phe(90.0, (0.0, 0.0, 5.0)), "interchain")
    assert x == Interaction("stack", ("A", "1", "ring:6"), ("B", "1", "ring:6"), 5.0, 90.0, "T-shaped")


def test_a_t_shaped_pair_off_both_axes_is_not_a_stack():
    """Centroid 3 Å off the face ring's axis and 4 Å off the edge ring's: both offsets over 2 Å."""
    assert find_pi_stacking(_two_phe(90.0, (3.0, 0.0, 4.0)), "interchain") == []
    [x] = find_pi_stacking(_two_phe(90.0, (1.5, 0.0, 4.0)), "interchain")
    assert x.note == "T-shaped" and x.distance == round(math.hypot(1.5, 4.0), 3)


@pytest.mark.parametrize(
    "angle, note",
    [
        (0.0, "parallel"),
        (29.9, "parallel"),
        (30.0, "parallel"),
        (30.0004, "parallel"),
        (30.0006, None),
        (45.0, None),
        (59.9994, None),
        (59.9996, "T-shaped"),
        (60.0, "T-shaped"),
        (75.0, "T-shaped"),
        (90.0, "T-shaped"),
    ],
)
def test_stacking_is_parallel_within_30_degrees_or_t_shaped_from_60(angle, note):
    found = find_pi_stacking(_two_phe(angle, (0.0, 0.0, 4.5)), "interchain")
    assert [x.note for x in found] == ([note] if note else [])
    if note:
        assert found[0].angle == round(angle, 3) and found[0].distance == 4.5


@pytest.mark.parametrize(
    "height, count",
    [
        (5.4, 1),
        (5.5, 1),
        (math.nextafter(5.5, 6.0), 1),
        (math.nextafter(5.5, 5.0), 1),
        (5.5004, 1),
        (5.5006, 0),
        (5.6, 0),
    ],
)
def test_ring_centroids_stack_within_5_5_angstroms(height, count):
    found = find_pi_stacking(_two_phe(0.0, (0.0, 0.0, height)), "interchain")
    assert len(found) == count
    if count:
        assert found[0].distance == round(height, 3)


@pytest.mark.parametrize("offset, count", _EDGE_OFFSETS)
def test_stacked_rings_are_at_most_2_angstroms_offset(offset, count):
    assert len(find_pi_stacking(_two_phe(0.0, (offset, 0.0, 3.6)), "interchain")) == count


def test_the_stacking_limits_are_the_stated_ones():
    assert (STACK_MAX, STACK_PARALLEL, STACK_T_MIN, STACK_OFFSET, HYDROPHOBIC_MAX) == (5.5, 30.0, 60.0, 2.0, 4.0)


_PYRIMIDINE = ("N1", "C2", "N3", "C4", "C5", "C6")


def test_the_stacked_bases_of_a_duplex_are_left_out():
    """Two bases stacked 3.4 Å apart, one in each strand: base-base stacking is the ladder's, so never a line, under
    any selection; a tyrosine over one of them stacks."""
    strands = (
        ("A", 1, "DC", _ring_at(0.0, (0, 0, 0), _PYRIMIDINE)),
        ("B", 1, "DC", _ring_at(0.0, (0.5, 0, 3.4), _PYRIMIDINE)),
    )
    for select in ("interchain", "chains:A+B", "na"):
        assert find_pi_stacking(_model_of(*strands), select) == []
    tyr = ("C", 1, "TYR", _ring_at(0.0, (0, 0.3, -3.5)))
    for select in ("interchain", "na"):
        [x] = find_pi_stacking(_model_of(*strands, tyr), select)
        assert (x.a, x.b, x.note) == (("A", "1", "ring:6"), ("C", "1", "ring:6"), "parallel")


def test_no_lambda_stack_pairs_two_bases_and_each_meets_the_rule():
    model, _ = _lambda()
    names = {(c, t): n for (c, t, _), (n, _) in _atoms(model).items()}
    bases = {"DA", "DC", "DG", "DT"}
    for select in ("na", "interchain"):
        found = find_pi_stacking(_lambda(), select)
        for x in found:
            assert x.kind == "stack" and x.distance <= STACK_MAX and x.a[:2] != x.b[:2]
            assert x.a[2] in ("ring:5", "ring:6") and x.b[2] in ("ring:5", "ring:6")
            assert not (names[x.a[:2]] in bases and names[x.b[:2]] in bases)
            parallel = x.angle <= STACK_PARALLEL and x.note == "parallel"
            assert parallel or (x.angle >= STACK_T_MIN and x.note == "T-shaped")
        assert verify_interactions(_lambda(), found) == []


def test_a_purine_gives_both_rings_to_a_protein_ring():
    purine = {
        "N9": (-1.1, 0.7, 0),
        "C8": (-2.0, -0.3, 0),
        "N7": (-1.3, -1.4, 0),
        "C5": (0.0, -1.0, 0),
        "C4": (0.1, 0.4, 0),
        "N3": (1.4, 0.9, 0),
        "C2": (2.6, 0.0, 0),
        "N1": (2.4, -1.4, 0),
        "C6": (1.1, -1.9, 0),
    }
    phe = ("B", 1, "PHE", _ring_at(0.0, (0.0, -0.3, 3.5)))
    found = find_pi_stacking(_model_of(("A", 1, "DG", purine), phe), "na")
    assert [x.a for x in found] == [("A", "1", "ring:5"), ("A", "1", "ring:6")]
    assert all(x.b == ("B", "1", "ring:6") for x in found)


# hydrophobic contacts

_LEU = {
    "N": (-0.5, 1.4, 0.0),
    "CA": (0.0, 0.0, 0.0),
    "C": (-0.5, -1.4, 0.0),
    "O": (-1.7, -1.6, 0.0),
    "CB": (1.5, 0.0, 0.0),
    "CG": (2.0, 1.4, 0.0),
    "CD1": (3.5, 1.4, 0.0),
    "CD2": (1.5, 2.2, -1.2),
}
# over the Leu: CG1 3.8 Å from Leu CG, CG2 3.905 Å from Leu CD1; the backbone well away
_VAL = {
    "N": (3.2, 2.8, 6.7),
    "CA": (3.2, 1.4, 6.2),
    "C": (3.2, 0.0, 6.7),
    "O": (3.2, -0.4, 7.9),
    "CB": (3.2, 1.4, 4.7),
    "CG1": (2.0, 1.4, 3.8),
    "CG2": (4.4, 1.4, 3.8),
}


def test_a_leucine_and_a_valine_touch_once_at_their_closest_carbons():
    found = find_hydrophobic(_model_of(("A", 1, "LEU", _LEU), ("B", 1, "VAL", _VAL)), "interchain")
    assert found == [Interaction("hydrophobic", ("A", "1", "CG"), ("B", "1", "CG1"), 3.8, None)]
    # the residue that comes first is a, whatever the model order
    found = find_hydrophobic(_model_of(("B", 1, "VAL", _VAL), ("A", 1, "LEU", _LEU)), "interchain")
    assert found == [Interaction("hydrophobic", ("A", "1", "CG"), ("B", "1", "CG1"), 3.8, None)]


def test_only_carbons_bonded_only_to_carbon_are_hydrophobic():
    """Ser CB (bonded to OG), Met CE (to SD), Leu CA (to N) and a backbone C (to O) never touch; Met CB and Ala CB
    do."""
    ala = ("A", 1, "ALA", {"CB": (0.0, 0.0, 0.0)})
    for name, atoms in (
        ("SER", {"CB": (3.6, 0.0, 0.0), "OG": (5.0, 0.5, 0.0)}),
        ("MET", {"CE": (3.6, 0.0, 0.0), "SD": (5.4, 0.3, 0.0)}),
        ("LEU", {"N": (5.0, 1.0, 0.0), "CA": (3.6, 0.0, 0.0)}),
        ("GLY", {"C": (3.6, 0.0, 0.0), "O": (4.4, 0.9, 0.0)}),
    ):
        assert find_hydrophobic(_model_of(ala, ("B", 1, name, atoms)), "interchain") == [], name
    met = ("B", 1, "MET", {"CB": (3.6, 0.0, 0.0), "CG": (5.1, 0.2, 0.0), "SD": (5.9, 1.8, 0.0)})
    [x] = find_hydrophobic(_model_of(ala, met), "interchain")
    assert (x.a, x.b) == (("A", "1", "CB"), ("B", "1", "CB"))


@pytest.mark.parametrize("dist, count", [(3.9, 1), (4.0, 1), (_UP, 1), (_DOWN, 1), (4.0004, 1), (4.0006, 0), (4.1, 0)])
def test_a_hydrophobic_contact_is_within_4_angstroms(dist, count):
    model = _model_of(("A", 1, "ALA", {"CB": (0, 0, 0)}), ("B", 1, "ALA", {"CB": (dist, 0, 0)}))
    found = find_hydrophobic(model, "interchain")
    assert len(found) == count
    if count:
        assert found[0].distance == round(dist, 3) and found[0].angle is None


def test_each_lambda_residue_pair_touches_once_and_meets_the_rule():
    for select in ("na", "interchain"):
        found = find_hydrophobic(_lambda(), select)
        assert found
        pairs = [(x.a[:2], x.b[:2]) for x in found]
        assert len(set(pairs)) == len(pairs)
        assert all(x.a[2].startswith("C") and x.b[2].startswith("C") and x.distance <= 4.0 for x in found)
    # the dimer interface is hydrophobic in part
    assert any({x.a[0], x.b[0]} == {"3", "4"} for x in find_hydrophobic(_lambda(), "interchain"))


# the independent re-check


def _all_lambda(select):
    found = [x for finder in _FINDERS for x in finder(_lambda(), select)]
    found += [
        Interaction(
            "hbond",
            (h.donor.res.chain, h.donor.res.tag, h.donor.name),
            (h.acceptor.res.chain, h.acceptor.res.tag, h.acceptor.name),
            h.distance,
            None,
        )
        for h in find_hbonds(_lambda(), select)
    ]
    return found


def test_every_line_the_finders_give_passes_the_re_check():
    for select in ("na", "interchain"):
        found = _all_lambda(select)
        assert {x.kind for x in found} >= {"salt", "cation-pi", "hydrophobic", "hbond"}
        assert verify_interactions(_lambda(), found) == []


def test_lines_at_the_limits_pass_the_re_check():
    lines = []
    for model in (
        _two_phe(0.0, (0.0, 0.0, 5.5004)),
        _two_phe(0.0, (2.0004, 0.0, 3.6)),
        _two_phe(30.0004, (0.0, 0.0, 4.5)),
        _two_phe(59.9996, (0.0, 0.0, 4.5)),
        _model_of(("A", 1, "ALA", {"CB": (0, 0, 0)}), ("B", 1, "ALA", {"CB": (4.0004, 0, 0)})),
        _model_of(("A", 1, "LYS", {"NZ": (0, 0, 0)}), ("B", 1, "ASP", {"OD1": (4.0004, 0, 0)})),
        _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (2.0004, 0.0, 3.6)})),
    ):
        found = [x for finder in _FINDERS for x in finder(model, "interchain")]
        assert found  # each holds its limit case (a Phe ring's carbons are hydrophobic too)
        assert verify_interactions(model, found) == []
        lines.extend(found)
    assert {x.kind for x in lines} == {"stack", "hydrophobic", "salt", "cation-pi"}
    assert {x.distance for x in lines if x.kind in ("salt", "hydrophobic") and x.a[2] in ("NZ", "CB")} == {4.0}


def test_the_re_check_flags_a_line_moved_onto_the_wrong_residue():
    [bridge] = [x for x in find_salt_bridges(_lambda(), "na") if (x.a, x.b) == (("3", "19", "NZ"), ("1", "4", "OP1"))]
    assert verify_interactions(_lambda(), [bridge]) == []
    moved = [
        dataclasses.replace(bridge, b=("1", "5", "OP1")),  # the next phosphate along the strand
        dataclasses.replace(bridge, a=("3", "20", "NZ")),  # a residue with no NZ
        dataclasses.replace(bridge, a=("4", "19", "NZ")),  # the same lysine of the other monomer
    ]
    flagged = verify_interactions(_lambda(), [bridge, *moved])
    assert len(flagged) == 3
    assert flagged[0].startswith("salt 3:19:NZ - 1:5:OP1: ")
    assert flagged[1].startswith("salt 3:20:NZ - 1:4:OP1: ")
    assert flagged[2].startswith("salt 4:19:NZ - 1:4:OP1: ")


def test_the_re_check_flags_each_broken_rule():
    lys_phe = _model_of(("A", 1, "PHE", _benzene()), ("B", 1, "LYS", {"NZ": (1.0, 0.0, 3.6)}))
    [cation_pi] = find_cation_pi(lys_phe, "interchain")
    duplex = _model_of(
        ("A", 1, "DC", _ring_at(0.0, (0, 0, 0), _PYRIMIDINE)), ("B", 1, "DC", _ring_at(0.0, (0, 0, 3.4), _PYRIMIDINE))
    )
    ser = _model_of(("A", 1, "ALA", {"CB": (0, 0, 0)}), ("B", 1, "SER", {"CB": (3.6, 0, 0), "OG": (5.0, 0.5, 0)}))
    cases = [
        (lys_phe, dataclasses.replace(cation_pi, distance=cation_pi.distance + 0.01), "distance"),
        (lys_phe, dataclasses.replace(cation_pi, b=("B", "1", "ring:6")), "ring"),
        (lys_phe, dataclasses.replace(cation_pi, kind="salt", b=("A", "1", "CG")), "anion"),
        (duplex, Interaction("stack", ("A", "1", "ring:6"), ("B", "1", "ring:6"), 3.4, 0.0, "parallel"), "nucleic"),
        (
            _two_phe(45.0, (0, 0, 4.5)),
            Interaction("stack", ("A", "1", "ring:6"), ("B", "1", "ring:6"), 4.5, 45.0, "parallel"),
            "angle",
        ),
        (
            _two_phe(0.0, (0, 0, 4.5)),
            Interaction("stack", ("A", "1", "ring:6"), ("B", "1", "ring:6"), 4.5, 0.0, "T-shaped"),
            "T-shaped",
        ),
        (ser, Interaction("hydrophobic", ("A", "1", "CB"), ("B", "1", "CB"), 3.6, None), "bonded"),
        (ser, Interaction("hydrophobic", ("A", "1", "CB"), ("A", "1", "CB"), 0.0, None), "itself"),
        (ser, Interaction("hydrophobic", ("A", "1", "CB"), ("C", "1", "CB"), 3.6, None), "no residue"),
        (ser, Interaction("hbond", ("A", "1", "CB"), ("B", "1", "OG"), 5.025, None), "H-bond"),
    ]
    for model, line, word in cases:
        flagged = verify_interactions(model, [line])
        assert len(flagged) == 1 and word in flagged[0], (line, flagged)
