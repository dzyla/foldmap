"""Contact beads against the deposited coordinates: every bead backed by atoms, no touched nucleotide unmarked."""

from functools import lru_cache
from pathlib import Path

import gemmi
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


def _build(name, extent="auto", plant=None, **look_kw):
    """(path, layout, backbone, figure, model) of the figure, drawn. plant(layout, bb) may alter the contacts
    before drawing."""
    path, look = _path(name), Style(na_extent=extent, **look_kw)
    model = read_model(path)
    lay, sses, bb = make_layout(path, look=look, model=model)
    if plant:
        plant(lay, bb)
    fig = draw(lay, route_loops(lay, sses, bb), sses, look=look)
    return path, lay, bb, fig, model


def _check(built):
    path, lay, bb, fig, model = built
    return verify_contacts(model, getattr(fig, "foldmap_beads", None), lay, bb, path=path)


@lru_cache(maxsize=None)
def _report(name, extent="auto"):
    built = _build(name, extent)
    return getattr(built[3], "foldmap_beads", None), _check(built)


@pytest.mark.parametrize("name", WITH_PROTEIN)
def test_every_bead_is_backed_by_atoms_and_no_contact_is_missed(name):
    drawn, report = _report(name)
    assert len(report.beads) == sum(len(v) for v in drawn.marks.values()) > 0
    assert [(b.nucleotide, b.problem) for b in report.beads if b.problem] == []
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

    report = _check(_build("1ZAA", plant=plant))
    assert [(b.nucleotide, b.residue) for b in report.false] == [(planted["nucleotide"], planted["residue"])]
    assert len(report.verified) == len(report.beads) - 1


def test_a_wrong_colour_is_reported():
    built = _build("1ZAA")
    marks = built[3].foldmap_beads.marks
    key, r, colour = marks["C"][0]
    marks["C"][0] = (key, r, "#123456")
    report = _check(built)
    assert [(b.key, b.colour) for b in report.wrong_colour] == [(key, "#123456")]
    assert report.false == [] and len(report.verified) == len(report.beads) - 1


def test_a_bead_that_does_not_name_the_first_toucher_is_reported():
    built = _build("3CRX")
    drawn, report = built[3].foldmap_beads, _check(built)
    later = next(m for m in report.missed if m.reason == "first toucher")
    bb = built[2]
    r = next(k for k, l in enumerate(bb.labels) if f"{l.chain}:{l.name}{l.seq}{l.icode}" == later.pair.residue)
    chain = next(c for c, v in drawn.marks.items() for b in v if _nt(built, b[0]) == later.nucleotide)
    drawn.marks[chain] = [(k, r if _nt(built, k) == later.nucleotide else x, c) for k, x, c in drawn.marks[chain]]
    report = _check(built)
    assert [b.nucleotide for b in report.not_first] == [later.nucleotide]
    assert report.false == [] and report.wrong_colour == []


def _nt(built, key):
    lab = built[1].nucleic.strands[key[0]].labels[key[1]]
    return f"{lab.chain}:{lab.name}{lab.tag}"


def _rename(built, key, name):
    labels = built[1].nucleic.strands[key[0]].labels
    labels[key[1]] = labels[key[1]]._replace(name=name)


def _first_bead(built):
    marks = built[3].foldmap_beads.marks
    chain = sorted(marks)[0]
    return marks[chain][0][0]


def test_a_nucleotide_the_model_cannot_supply_is_unmatched_not_dropped():
    built = _build("1ZAA")
    key = _first_bead(built)
    _rename(built, key, "QQQ")  # the figure's label no longer matches the model
    report = _check(built)
    assert [b.key for b in report.unmatched] == [key] and len(report.beads) == sum(
        len(v) for v in built[3].foldmap_beads.marks.values()
    )


def test_a_bead_on_an_untabulated_nucleotide_is_verified():
    built = _build("1ZAA")
    key = _first_bead(built)
    lab = built[1].nucleic.strands[key[0]].labels[key[1]]
    for res in built[4][0][lab.chain]:
        if res.seqid.num == lab.seq and res.name == lab.name:
            res.name = "ZZZ"  # a code gemmi does not know, as in a modified nucleotide
    _rename(built, key, "ZZZ")
    report = _check(built)
    bead = next(b for b in report.beads if b.key == key)
    assert not gemmi.find_tabulated_residue("ZZZ").is_nucleic_acid()  # gemmi does not know the code
    assert bead.verified and bead.nucleotide.startswith("A:ZZZ")


def test_a_missing_bead_is_reported_and_a_focus_exclusion_is_not():
    built = _build("3CRX")
    path, lay, bb, fig, model = built
    chains = sorted({l.chain for l in bb.labels})
    assert len(chains) > 1
    marks = fig.foldmap_beads.marks
    key, r, _ = next(iter(marks.values()))[0]
    chain = bb.labels[r].chain
    for c in marks:  # take that nucleotide's bead away
        marks[c] = [b for b in marks[c] if b[0] != key]
    assert chain in {m.pair.residue.split(":")[0] for m in _check(built).missing}
    lay.focus_chains = {c for c in chains if c != chain}  # that chain is outside the drawn subunit
    report = _check(built)
    assert chain not in {m.pair.residue.split(":")[0] for m in report.missing}
    assert "focus" in {m.reason for m in report.excluded if m.pair.residue.startswith(chain + ":")}


def test_cropped_nucleotides_are_listed_not_counted():
    _, report = _report("4OO8", "contacts")
    assert report.cropped and all(not m.counted for m in report.cropped) and report.missing == []
    assert {m.reason for m in report.cropped} == {"cropped"}


def test_deposited_hydrogen_bonds_are_read():
    _, report = _report("3CRX")
    assert report.hydrog and all(" " in h.a and ":" in h.a for h in report.hydrog)


def _page_data(**kw):
    import json
    import re

    from foldmap.interactive import build_page

    page = build_page(_path("1ZAA"), **kw)
    return page, json.loads(re.search(r'id="topo-data">(.*?)</script>', page, re.S).group(1))


def test_interactive_page_names_the_toucher_of_each_bead():
    import re

    page, data = _page_data()
    beads = data["beads"]["chains"]["C"]
    assert "DG7 – contacted by C:Arg14 (2.5 Å)" in {b["text"] for b in beads}
    group = re.search(r'<g id="na:contact:C">(.*?)</g>\s*<g id=', page, re.S).group(1)
    assert group.count("<use") == len(beads)  # one hover note per drawn bead, in drawing order
    assert all(0 <= b["k"] < len(data["residues"]) for b in beads)
    assert data["beads"]["marks"]["A:7"]["text"] == "DG7 – contacted by C:Arg14 (2.5 Å)"


def test_letters_style_beads_have_the_hover_too():
    page, data = _page_data(look=Style(na_style="letters"))
    assert data["beads"]["chains"]["C"] and 'id="na:contact:' not in page  # no beads, only bold letters
    for mark in data["beads"]["marks"]:
        assert f'id="na:nt:{mark}"' in page
    assert "na:nt:" in page.split("itemAt")[1]  # the page script maps those letters to the note
