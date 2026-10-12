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


def test_nucleosome_is_drawn_as_one_half_and_says_so():
    """The author's report: a nucleosome (3UTB) drawn as one copy of each histone. focus auto finds the C2 axis, and
    eight protein chains are more than FOCUS_CHAINS (6): one protomer (H3, H4, H2A, H2B) in full, the other four
    as grey neighbours, both DNA strands. 6 of the 10 chains are drawn."""
    s = scope_of(NUCLEOSOME)
    assert (s["drawn_chains"], s["file_chains"]) == (6, 10)
    assert s["reason"] == "one of 2 subunits (C2 symmetry)"
    assert s["all_options"] == {"focus": "none"}
    whole = scope_of(NUCLEOSOME, focus="none")
    assert (whole["drawn_chains"], whole["file_chains"], whole["reason"], whole["all_options"]) == (10, 10, "", {})


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
    nucleosome = scope_of(NUCLEOSOME)
    assert (
        scope_note(nucleosome) == "drew one of 2 subunits (C2 symmetry): 6 of 10 chains; --focus none draws all of them"
    )
