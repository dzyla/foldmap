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


def test_na_crop_is_announced(capsys, tmp_path):
    """na_extent auto crops a large nucleic acid to the stems the protein touches (4YBB: hundreds of nucleotides
    left out). That is an automatic choice like any other: drawn_scope says it, with the setting that draws all of
    it, the command line prints it, and figcheck.scope_honest finds nothing left out silently (plan 10 Task 8)."""
    from foldmap import figcheck
    from foldmap.cli import main
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "na" / "4YBB.cif.gz"
    d = draw_spec(figure_spec(path))
    assert d.layout.na_cropped and d.layout.na_read > d.layout.nucleic.size
    s = drawn_scope(d.model, d.layout)
    size, read = d.layout.nucleic.size, d.layout.na_read
    want = f"the RNA cropped to the stems the protein touches ({size:,} of {read:,} nucleotides)"
    assert s["reason"] == want and s["all_options"] == {"set": "na_extent=all"}
    assert scope_note(s) == f"drew {want}; --set na_extent=all draws all of them"
    assert scope_note(s, how="mcp").endswith('; style={"na_extent": "all"} draws all of them')
    assert figcheck.scope_honest(d, path) == []
    assert main(["plot", str(path), "-o", str(tmp_path / "4YBB.svg")]) == 0
    assert f"foldmap: drew {want}; --set na_extent=all draws all of them" in capsys.readouterr().err
    # the crop covers its ends and the gaps between the pieces it kept, not a run missing inside a kept piece
    stem = max(d.layout.stems, key=lambda st: len(st.axial))
    keys = sorted(stem.axial)
    for key in keys[:4]:
        del stem.axial[key]
    got = figcheck.scope_honest(d, path)
    assert len(got) == 1 and got[0].startswith("nt "), got


def test_a_crop_and_a_left_out_chain_are_said_apart():
    """When a chain count and a crop are both said, the count goes with what was drawn of the chains, never after
    the crop's clause (fix round 1)."""
    head, crop = (
        "one of 6 subunits of a helical filament",
        "the RNA cropped to the stems the protein touches (8 of 90 nucleotides)",
    )
    s = {
        "file_chains": 6,
        "drawn_chains": 2,
        "reason": f"{head}; {crop}",
        "chains_reason": head,
        "crop": crop,
        "all_options": {"focus": "none", "set": "na_extent=all"},
        "residues": 900,
    }
    assert scope_note(s) == (
        "drew one of 6 subunits of a helical filament: 2 of 6 chains; the RNA cropped to the stems the protein "
        "touches (8 of 90 nucleotides); --focus none --set na_extent=all draws all of them"
    )


def test_chains_asked_with_a_crop_count_the_residues_of_those_chains():
    """--chains with na_extent auto's crop: "draw all" draws the asked-for chains whole, so the residues counted are
    theirs, not the file's (fix round 1)."""
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "na" / "4YBB.cif.gz"
    d = draw_spec(figure_spec(path))
    protein = sorted({p.sse.chain for p in d.layout.placed.values()})[:1]
    rna = sorted({st.chain for st in d.layout.nucleic.strands})
    d = draw_spec(figure_spec(path, chains=protein + rna))
    s = drawn_scope(d.model, d.layout)
    assert s["crop"] and s["chains_reason"].startswith("chains ")
    want = sum(len(c.get_polymer()) for c in d.model[0] if c.name in d.layout.chains_read)
    assert s["residues"] == want < sum(len(c.get_polymer()) for c in d.model[0])


def test_explicit_chains_note_says_what_was_asked():
    """--chains A draws what was asked: no "left out" note, the note says "chains A (asked for)" (F1 review nit);
    4CH9 --chains B is not "one of 2 copies" any more."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    s = scope_of(DATA / "4CH9.cif", chains=["B"])
    assert s["reason"] == "chains B (asked for)" and s["all_options"] == {}
    assert scope_note(s) == "drew chains B (asked for)"  # no "one of 2 copies", nothing "left out"
    d = draw_spec(figure_spec(DATA / "2HHB.cif", chains=["A", "C"]))
    s = drawn_scope(d.model, d.layout)
    assert s["reason"] == "chains A, C (asked for)" and "draws all" not in scope_note(s)
    assert figcheck.scope_honest(d, DATA / "2HHB.cif", chains=["A", "C"]) == []


@pytest.mark.parametrize("path, n", [(DATA / "assemblies" / "7A4M.cif", 24), (DATA / "6Y7S.cif", 6)])
def test_subunits_are_counted_as_drawn(path, n):
    """The number of subunits a focus says ("one of N"): 7A4M's cage 24, 6Y7S's rod 6."""
    lay, _, _ = make_layout(path)
    assert lay.subunits[0] == n


def test_half_a_subunit_rounds_up():
    """Nine chains in focus units of two: four and a half units, said as five (round() took it to four)."""
    from types import SimpleNamespace

    from foldmap.cli import _subunits

    sym = SimpleNamespace(protomers=[["A", "B", "C"], ["D", "E", "F"], ["G", "H", "I"]])
    assert _subunits(sym, {"A", "B"}) == 5
    assert _subunits(SimpleNamespace(protomers=[["A", "B", "C"]] * 5), {"A", "B"}) == 8  # 7.5
    assert _subunits(SimpleNamespace(protomers=[["A"]] * 24), {"A"}) == 24


def test_a_chain_with_no_helix_or_strand_is_said_not_drawn():
    """C15 (random-100): a protein chain with no helix or strand has no element, so it was drawn as nothing and
    said by nothing (6B5W, a 19-residue lasso peptide: an empty page; 7PIU chain P, 7MW8 chain K). The scope says
    it is not drawn and why, the note prints it, and figcheck.scope_honest counts it as said."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    path = DATA / "r100" / "6B5W.cif.gz"
    s = scope_of(path)
    assert s["undrawn_chains"] == ["A"] and s["not_drawn"] == "chain A (19 residues, no helix or strand)"
    assert s["reason"] == "not drawn: chain A (19 residues, no helix or strand)"
    assert (s["drawn_chains"], s["file_chains"]) == (0, 1) and s["all_options"] == {}
    assert scope_note(s) == "not drawn: chain A (19 residues, no helix or strand)"
    assert figcheck.scope_honest(draw_spec(figure_spec(path)), path) == []
    assert scope_of(DATA / "1UBQ.cif")["not_drawn"] == ""


def test_a_nucleic_acid_chain_too_short_to_draw_is_said(tmp_path):
    """C15 (random-100): 7MW8's chain K is one nucleotide, read but drawn by no strand; it is said too."""
    import gemmi

    st = gemmi.read_structure(str(DATA / "na" / "1URN.cif.gz"))
    rna = st[0]["R"].clone()
    for i in reversed(range(1, len(rna))):
        del rna[i]
    rna.name = "K"
    ubq = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    ubq.remove_ligands_and_waters()
    ubq[0].add_chain(rna)
    ubq.setup_entities()
    path = tmp_path / "1UBQ-K.cif"
    ubq.make_mmcif_document().write_file(str(path))
    s = scope_of(path)
    assert s["undrawn_chains"] == ["K"] and s["not_drawn"] == "chain K (1 nucleotide, too short to draw)"
