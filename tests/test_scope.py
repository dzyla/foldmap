"""What an automatic choice left out of the figure (scope.drawn_scope): assembly 1 drawn as one of several copies of
a complex, or one subunit of a large symmetric assembly drawn with its neighbours as ghosts."""

from pathlib import Path

import pytest

from foldmap.cli import make_layout
from foldmap.io import read_model
from foldmap.scope import drawn_scope, scope_note

DATA = Path(__file__).parent / "data"
NUCLEOSOME = DATA / "assemblies" / "3UTB.cif.gz"  # 8 histones + 2 DNA strands, C2: drawn as one half by focus auto


def scope_of(path, assembly="auto", focus="auto", **kw):
    read = read_model(path, assembly)
    lay, _, _ = make_layout(path, model=read, assembly=assembly, focus=focus, **kw)
    return drawn_scope(read, lay, assembly=assembly, focus=focus)


def test_kelch_assembly_one_is_one_of_two_copies():
    """4CH9: two Kelch domains with their WNK4 peptides deposited; assembly 1 (B, D) is one copy of the complex."""
    s = scope_of(DATA / "4CH9.cif")
    assert (s["drawn_chains"], s["file_chains"]) == (2, 4)
    assert s["reason"] == "assembly 1, one of 2 copies in the file"
    assert s["all_options"] == {"assembly": "asu", "focus": "none"}
    assert s["residues"] == sum(len(c.get_polymer()) for c in read_model(DATA / "4CH9.cif", "asu")[0])


def test_apoferritin_is_one_subunit_of_its_cage():
    s = scope_of(DATA / "assemblies" / "7A4M.cif")
    assert (s["drawn_chains"], s["file_chains"]) == (1, 24)
    assert s["reason"] == "one of 24 subunits of an O cage"
    assert s["all_options"] == {"focus": "none"}
    assert s["residues"] == 24 * len(read_model(DATA / "assemblies" / "7A4M.cif", "asu")[0]["A"].get_polymer())


def test_pilus_rod_is_one_subunit_of_its_filament():
    s = scope_of(DATA / "6Y7S.cif")
    assert (s["drawn_chains"], s["file_chains"]) == (1, 6)
    assert s["reason"] == "one of 6 subunits of a helical filament"
    assert s["all_options"] == {"focus": "none"}


def test_nucleosome_is_drawn_whole_with_nothing_to_report():
    """The author's report: a nucleosome (3UTB) was drawn as one copy of each histone (focus auto took one half of
    the C2 octamer, as 8 protein chains are more than FOCUS_CHAINS). focus auto never reduces a two-subunit
    assembly, nor splits one whose DNA both halves hold: all 8 histones and both DNA strands, nothing reported."""
    read = read_model(NUCLEOSOME)
    lay, sses, bb = make_layout(NUCLEOSOME, model=read)
    assert lay.focus_chains is None and not lay.partial
    assert {s.chain for s in sses} == set("ABCDEFGH") and {s.chain for s in lay.nucleic.strands} == {"I", "J"}
    assert not any(p.ghost for p in lay.placed.values())
    s = drawn_scope(read, lay)
    assert (s["drawn_chains"], s["file_chains"], s["reason"], s["all_options"]) == (10, 10, "", {})


@pytest.mark.parametrize("path", [DATA / "na" / "4OO8.cif.gz", DATA / "na" / "1KX5.cif.gz", DATA / "1UBQ.cif"])
def test_nothing_left_out_reports_nothing(path):
    s = scope_of(path)
    assert s["drawn_chains"] == s["file_chains"] and s["reason"] == "" and s["all_options"] == {}
    assert scope_note(s) == ""


def test_choices_the_user_made_are_not_reported():
    """Only automatic choices are reported: an assembly or a focus asked for by name draws what was asked."""
    assert scope_of(DATA / "4CH9.cif", assembly="1")["reason"] == ""
    assert scope_of(DATA / "4CH9.cif", assembly="asu")["file_chains"] == 4
    assert scope_of(DATA / "assemblies" / "7A4M.cif", focus="A1")["reason"] == ""
    assert scope_of(DATA / "6Y7S.cif", focus="none")["reason"] == ""


def test_scope_notes_name_the_settings_that_draw_everything():
    kelch = scope_of(DATA / "4CH9.cif")
    assert scope_note(kelch) == (
        "drew assembly 1, one of 2 copies in the file: 2 of 4 chains; --assembly asu --focus none draws all of them"
    )
    assert scope_note(kelch, how="mcp") == (
        'drew assembly 1, one of 2 copies in the file: 2 of 4 chains; assembly="asu", focus="none" draws all of them'
    )
    apo = scope_of(DATA / "assemblies" / "7A4M.cif")
    assert scope_note(apo) == "drew one of 24 subunits of an O cage: 1 of 24 chains; --focus none draws all of them"
