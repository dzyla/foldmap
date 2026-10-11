"""Hydrogen bonds from the coordinates (foldmap.hbonds): known interactions found, selections honoured, waters
ignored, output deterministic and fast."""

import os
import subprocess
import sys
import time
from functools import lru_cache
from pathlib import Path

import gemmi
import pytest

from foldmap.hbonds import HBond, Selection, find_hbonds, parse_selection
from foldmap.io import read_model

DATA = Path(__file__).parent / "data"


def _path(name):
    return next(p for p in (DATA / "na" / f"{name}.cif.gz", DATA / f"{name}.cif") if p.is_file())


@lru_cache(maxsize=None)
def _model(name):
    return read_model(_path(name))


@lru_cache(maxsize=None)
def _found(name, select):
    return tuple(find_hbonds(_model(name), select))


def _tag(atom):
    return f"{atom.res.chain}:{atom.res.name}{atom.res.tag} {atom.name}"


def _protein(atom) -> bool:
    info = gemmi.find_tabulated_residue(atom.res.name)
    return info is not None and info.is_amino_acid()


_BACKBONE = {"N", "CA", "C", "O", "OXT"}
_NOT_BASE = {"P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}

# (donor, acceptor, distance in Å as measured in the file)
LAMBDA_SIDE_CHAIN_TO_BASE = {
    ("1:DA4 N6", "3:GLN44 OE1", 3.12),
    ("3:GLN44 NE2", "1:DA4 N7", 3.32),
    ("2:DA24 N6", "4:GLN44 OE1", 2.94),
    ("4:GLN44 NE2", "2:DA24 N7", 2.91),
    ("3:SER45 OG", "2:DG36 O6", 2.60),
    ("3:SER45 OG", "2:DG37 O6", 2.99),
    ("1:DC5 N4", "3:SER45 OG", 3.04),
    ("4:SER45 OG", "1:DG16 N7", 3.18),
    ("3:ASN55 ND2", "2:DG34 N7", 3.06),
    ("4:ASN55 ND2", "1:DG14 N7", 2.91),
    ("4:LYS4 NZ", "1:DG13 O6", 2.80),
    ("4:LYS4 NZ", "1:DG13 N7", 3.33),
    ("4:LYS4 NZ", "1:DG14 O6", 2.90),
    ("4:LYS4 NZ", "1:DG14 N7", 3.42),
    ("4:LYS3 NZ", "1:DG11 O6", 2.57),
    ("4:LYS3 NZ", "2:DG30 O6", 2.99),
}


def test_lambda_repressor_side_chains_bond_the_major_groove_bases():
    """1LMB, lambda repressor (chains 3 and 4, residues numbered 1-92 as in the file) on its operator (DNA chains
    1 and 2). The file's primary citation is Beamer & Pabo, J. Mol. Biol. 227:177-196 (1992); its citation 1
    (Clarke et al., Science 254:267-270, 1991) is on the N-terminal arm. The interactions usually described for
    this complex, as they appear in this file (each checked against the coordinates, distances measured
    heavy atom to heavy atom):

    - Gln44 to the adenine of base pair 2: Gln44 OE1 accepts from DA N6 and Gln44 NE2 donates to DA N7, in both
      half-sites (3:GLN44 with 1:DA4; 4:GLN44 with 2:DA24).
    - Ser45 OG to guanine O6/N7 (3:SER45 with 2:DG36, 2:DG37 and, as acceptor, 1:DC5 N4; 4:SER45 with 1:DG16).
      4:SER45 OG is also 2.88 Å from 1:DG16 O6, but at a CB-OG-O6 angle of 86.6°, under the 90° rule, so of that
      pair only N7 counts.
    - Asn55 ND2 to guanine N7 (3:ASN55 with 2:DG34; 4:ASN55 with 1:DG14).
    - Lys4 NZ of the N-terminal arm to guanine O6/N7 (4:LYS4 with 1:DG13 and 1:DG14; the arm of chain 3 starts
      at residue 6 and is not ordered), and Lys3 NZ to guanine O6 (4:LYS3 with 1:DG11 and 2:DG30).
    - Gln33 reaches the phosphate backbone, not a base: its main-chain N and side-chain NE2 donate to OP2
      (3:GLN33 with 1:DT3 and 1:DA4; 4:GLN33 with 2:DT23 and 2:DA24).

    With that one O6 pair left out, these are every side-chain N/O/S to base N/O pair of the file within
    2.5-3.5 Å (each a donor with an acceptor), so the side-chain to base bonds found must be exactly this set."""
    found = _found("1LMB", "na")
    side_to_base = set()
    for b in found:
        prot, nt = (b.donor, b.acceptor) if _protein(b.donor) else (b.acceptor, b.donor)
        if prot.name not in _BACKBONE and nt.name not in _NOT_BASE:
            side_to_base.add((_tag(b.donor), _tag(b.acceptor), round(b.distance, 2)))
    assert side_to_base == LAMBDA_SIDE_CHAIN_TO_BASE
    tags = {(_tag(b.donor), _tag(b.acceptor)) for b in found}
    assert {("3:GLN33 N", "1:DT3 OP2"), ("3:GLN33 NE2", "1:DA4 OP2")} <= tags
    assert {("4:GLN33 N", "2:DT23 OP2"), ("4:GLN33 NE2", "2:DA24 OP2")} <= tags


def test_na_bonds_join_protein_and_nucleic_acid():
    found = _found("1LMB", "na")
    assert found and all(_protein(b.donor) != _protein(b.acceptor) for b in found)
    assert all(2.5 <= b.distance <= 3.5 for b in found)


def test_haemoglobin_alpha1_beta1_interface():
    found = _found("2HHB", "chains:A+B")
    assert len(found) > 0
    for b in found:
        assert {b.donor.res.chain, b.acceptor.res.chain} == {"A", "B"}
        assert 2.5 <= b.distance <= 3.5
    assert found == _found("2HHB", "chains:B+A")
    assert set(found) <= set(_found("2HHB", "interchain"))


def test_interchain_covers_every_chain_pair_and_no_chain_with_itself():
    found = _found("2HHB", "interchain")
    pairs = {frozenset((b.donor.res.chain, b.acceptor.res.chain)) for b in found}
    assert all(len(p) == 2 for p in pairs)
    assert {frozenset("AB"), frozenset("CD")} <= pairs


def test_waters_are_ignored():
    model, _ = _model("1LMB")
    waters = [r for c in model for r in c if r.is_water()]
    assert waters  # the file has waters, some within bonding distance of the protein
    found = _found("1LMB", "interchain")
    assert found and all("HOH" not in (b.donor.res.name, b.acceptor.res.name) for b in found)


def test_each_atom_pair_once_and_sorted():
    found = _found("2HHB", "interchain")
    assert list(found) == sorted(found)
    pairs = [frozenset((b.donor, b.acceptor)) for b in found]
    assert len(pairs) == len(set(pairs))


def test_off_finds_nothing():
    assert find_hbonds(_model("1LMB"), "off") == []


_PRINT = """
import sys
from foldmap.hbonds import find_hbonds
from foldmap.io import read_model
print(repr(find_hbonds(read_model(sys.argv[1]), sys.argv[2])))
"""


def test_output_does_not_depend_on_the_hash_seed():
    runs = set()
    for seed in ("0", "12345", "777"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run(
            [sys.executable, "-c", _PRINT, str(_path("1LMB")), "interchain"],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1 and "HBond(" in runs.pop()


@pytest.mark.parametrize("select", ["na", "interchain"])
def test_ribosome_fixture_is_fast(select):
    model = _model("4YBB")
    start = time.perf_counter()
    found = find_hbonds(model, select)
    assert time.perf_counter() - start < 2.0
    assert found


@pytest.mark.parametrize(
    "text, canonical",
    [
        ("off", "off"),
        ("na", "na"),
        (" interchain ", "interchain"),
        ("chains:A+B", "chains:A+B"),
        ("chains: A+B , C+D", "chains:A+B,C+D"),
        ("chains:B+A", "chains:B+A"),
    ],
)
def test_selection_parses_and_prints_back(text, canonical):
    sel = parse_selection(text)
    assert isinstance(sel, Selection) and str(sel) == canonical
    assert parse_selection(str(sel)) == sel


@pytest.mark.parametrize(
    "text, words",
    [
        ("ligands", "reserved"),
        ("all", "off, na, interchain or chains:"),
        ("chains:", "chain pair"),
        ("chains:A", "chain pair"),
        ("chains:A+B,", "chain pair"),
        ("chains:A+B+C", "chain pair"),
        ("chains:A+A", "two different chains"),
    ],
)
def test_bad_selections_are_refused_clearly(text, words):
    with pytest.raises(ValueError, match="hbonds") as err:
        parse_selection(text)
    assert words in str(err.value)


def test_a_chain_the_model_lacks_is_named():
    with pytest.raises(ValueError, match="Z"):
        find_hbonds(_model("2HHB"), "chains:A+Z")


def test_hbond_carries_donor_acceptor_distance():
    b = _found("1LMB", "na")[0]
    assert isinstance(b, HBond) and b == HBond(b.donor, b.acceptor, b.distance)
