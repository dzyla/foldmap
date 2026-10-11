"""Contact beads against the deposited coordinates: every bead backed by atoms, no touched nucleotide unmarked."""

from functools import lru_cache
from pathlib import Path

import pytest

from foldmap.cli import make_layout
from foldmap.io import read_model
from foldmap.nacheck import CONTACT, verify_contacts
from foldmap.render import draw
from foldmap.route import route_loops
from foldmap.style import Style

DATA = Path(__file__).parent / "data"
WITH_PROTEIN = ["30QM", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB", "1LMB", "1ZAA"]
WITHOUT = ["1EHZ", "1ATP", "4YBB-psi"]


def _path(name):
    return next(p for p in (DATA / "na" / f"{name}.cif.gz", DATA / f"{name}.cif") if p.is_file())


def _build(name, extent="auto", plant=None):
    """(path, layout, backbone) of the figure, drawn. plant(layout, bb) may alter the contacts before drawing."""
    path, look = _path(name), Style(na_extent=extent)
    lay, sses, bb = make_layout(path, look=look)
    if plant:
        plant(lay, bb)
    draw(lay, route_loops(lay, sses, bb), sses, look=look)
    return path, lay, bb


@lru_cache(maxsize=None)
def _report(name, extent="auto"):
    path, lay, bb = _build(name, extent)
    return lay, verify_contacts(read_model(path), lay, bb, path=path)


@pytest.mark.parametrize("name", WITH_PROTEIN)
def test_every_bead_is_backed_by_atoms_and_no_contact_is_missed(name):
    lay, report = _report(name)
    assert len(report.beads) == sum(len(v) for v in lay.na_beads.values()) > 0
    assert [str(b.pair) for b in report.false] == []
    assert [(m.nucleotide, str(m.pair)) for m in report.missing] == []
    for b in report.beads:
        assert b.pair.dist <= CONTACT and b.pair.residue.split(":")[0] == b.chain


@pytest.mark.parametrize("name", WITHOUT)
def test_no_protein_no_beads(name):
    _, report = _report(name)
    assert report.beads == [] and report.missing == []


def test_closest_pair_is_recorded_per_bead():
    _, report = _report("1ZAA")
    bead = min(report.beads, key=lambda b: b.pair.dist)
    assert str(bead.pair) == "C:ARG14 NE – A:DG7 OP1 2.5 Å"
    assert (bead.residue, bead.nucleotide) == ("C:ARG14", "A:DG7")


def test_a_false_bead_is_reported():
    planted = {}

    def plant(lay, bb):  # a residue that touches nothing is made the toucher of the first nucleotide
        na = lay.nucleic
        far = next(r for r in range(len(bb)) if r not in na.contacts)
        key = (0, 0)
        for v in na.contacts.values():
            v.discard(key)
        na.contacts = {r: v for r, v in na.contacts.items() if v}
        na.contacts[far] = {key}
        na.contact_chain[far] = bb.labels[far].chain
        planted["residue"] = f"{bb.labels[far].chain}:{bb.labels[far].name}{bb.labels[far].seq}"
        planted["nucleotide"] = f"{na.strands[0].chain}:{na.strands[0].labels[0].name}{na.strands[0].labels[0].seq}"

    path, lay, bb = _build("1ZAA", plant=plant)
    report = verify_contacts(read_model(path), lay, bb)
    assert [(b.nucleotide, b.residue) for b in report.false] == [(planted["nucleotide"], planted["residue"])]
    assert len(report.verified) == len(report.beads) - 1


def test_a_missing_bead_is_reported_and_a_focus_exclusion_is_not():
    path, lay, bb = _build("3CRX")
    chains = sorted({l.chain for l in bb.labels})
    assert len(chains) > 1
    key, r, _ = next(iter(lay.na_beads.values()))[0]
    chain = bb.labels[r].chain
    for c in lay.na_beads:  # take that nucleotide's bead away
        lay.na_beads[c] = [b for b in lay.na_beads[c] if b[0] != key]
    model = read_model(path)
    assert chain in {m.pair.residue.split(":")[0] for m in verify_contacts(model, lay, bb).missing}
    lay.focus_chains = {c for c in chains if c != chain}  # that chain is outside the drawn subunit
    report = verify_contacts(model, lay, bb)
    assert chain not in {m.pair.residue.split(":")[0] for m in report.missing}
    assert "focus" in {m.reason for m in report.excluded if m.pair.residue.startswith(chain + ":")}


def test_cropped_nucleotides_are_listed_not_counted():
    _, report = _report("4OO8", "contacts")
    assert report.cropped and all(not m.counted for m in report.cropped) and report.missing == []
    assert {m.reason for m in report.cropped} == {"cropped"}


def test_deposited_hydrogen_bonds_are_read():
    _, report = _report("3CRX")
    assert report.hydrog and all(" " in h.a and ":" in h.a for h in report.hydrog)


def test_interactive_page_names_the_toucher_of_each_bead():
    import json
    import re

    from foldmap.interactive import build_page

    page = build_page(_path("1ZAA"))
    data = json.loads(re.search(r'id="topo-data">(.*?)</script>', page, re.S).group(1))
    beads = data["beads"]["C"]
    assert "DG7 – contacted by C:Arg14 (2.5 Å)" in {b["text"] for b in beads}
    group = re.search(r'<g id="na:contact:C">(.*?)</g>\s*<g id=', page, re.S).group(1)
    assert group.count("<use") == len(beads)  # one hover note per drawn bead, in drawing order
    assert all(0 <= b["k"] < len(data["residues"]) for b in beads)
