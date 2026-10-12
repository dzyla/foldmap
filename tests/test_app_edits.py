"""The app's structured controls write exactly the text the command line takes, and read it back."""

import re
from pathlib import Path

import pytest

from foldmap.app_edits import (
    COMPOSITES,
    DERIVED,
    EDITORS,
    apply_table,
    compose_into,
    control_for,
    element_choices,
    hex_colour,
    notes_text,
    palette_used,
    row_source,
    rows_of,
    sync_composites,
    text_of,
)

DATA = Path(__file__).parent / "data"

ROUND_TRIP = {
    "chains": ["", "A", "A,B"],
    "highlight": ["none", "asu", "protomer", "A,B"],
    "focus": ["auto", "none", "B", "B,C"],
    "ligands": ["auto", "all", "none", "HEM,ZN"],
    "na_backbone": ["modality", "chain", "ink", "#123456"],
    "hbonds": ["off", "na", "interchain", "chains:A+B,C+D"],
    "uniprot": ["off", "auto", "P0CG48"],
    "uniprot_sites": ["auto", "all", "none", "active,binding"],
}


@pytest.mark.parametrize("key", sorted(ROUND_TRIP))
def test_composites_round_trip_every_form(key):
    _, parts, text = COMPOSITES[key]
    for value in ROUND_TRIP[key]:
        assert text(*parts(value)) == value, (key, value)


def test_named_colours_become_hex_for_the_pickers():
    """st.color_picker takes #rrggbb only: a named colour (a layout file's teal) is converted, anything else is the
    picker's default."""
    _, parts, _ = COMPOSITES["na_backbone"]
    assert parts("teal") == ("colour", "#008080") and parts("#ABCDEF") == ("colour", "#abcdef")
    assert parts("nonsense") == ("colour", "#8c564b")
    assert hex_colour("teal", "#000000") == "#008080" and hex_colour("#fff", "#000000") == "#ffffff"
    assert hex_colour("", "#d020f8") == "#d020f8" and hex_colour(None, "#d020f8") == "#d020f8"


def test_composites_normalise_what_people_type():
    assert COMPOSITES["ligands"][2]("codes", "hem  zn,") == "HEM,ZN"
    assert COMPOSITES["hbonds"][2]("chains", "A + B, C+D") == "chains:A+B,C+D"
    assert COMPOSITES["highlight"][2]("chains", []) == "none"
    assert COMPOSITES["uniprot"][2]("accession", "p0cg48") == "P0CG48"


def test_sync_keeps_a_half_made_choice():
    ss = {"ligands": "auto"}
    sync_composites(ss)
    assert ss["ligands__mode"] == "auto" and ss["ligands__codes"] == ""
    ss["ligands__mode"] = "codes"  # chosen, nothing typed yet
    compose_into(ss, "ligands")
    assert ss["ligands"] == "auto" and ss["last_changed"] == "ligands"
    sync_composites(ss)
    assert ss["ligands__mode"] == "codes"  # not snapped back
    ss["ligands__codes"] = "hem, zn"
    compose_into(ss, "ligands")
    assert ss["ligands"] == "HEM,ZN"
    ss["ligands"] = "all"  # a layout file or a theme reset
    sync_composites(ss)
    assert ss["ligands__mode"] == "all"


EDITOR_TEXT = {
    "swap": ["α1,α3", "α1,α3; res:A:60-79,β2"],
    "move": ["α2=1,-2", "α2=1,-2; res:A:3-11=0.5,0"],
    "rename": ["res:167-182=Gd", "α1=A; res:A:60-79=C"],
    "mark": ["G", "G, A:G=#2ca02c, #3, res:150-159=red"],
    "domains": ["ZPN=res:A:331-440", "N=res:A:1-30\nC=res:A:31-76,α1"],
}


@pytest.mark.parametrize("kind", sorted(EDITOR_TEXT))
def test_editor_rows_round_trip_today_s_syntax(kind):
    for text in EDITOR_TEXT[kind]:
        rows = rows_of(kind, text)
        assert all(set(r) == set(EDITORS[kind]) for r in rows)
        assert text_of(kind, rows) == text, (kind, text)


def test_editor_text_takes_newlines_and_drops_empty_rows():
    assert rows_of("swap", "α1,α3\nα2,α4") == [{"a": "α1", "b": "α3"}, {"a": "α2", "b": "α4"}]
    assert text_of("rename", [{"element": "", "name": "x"}, {"element": "α1", "name": "Z"}]) == "α1=Z"
    assert text_of("move", [{"element": "α2", "dx": 1.0, "dy": -2.5}]) == "α2=1,-2.5"


def test_domain_ranges_become_residue_references():
    rows = [{"name": "N", "ranges": "A:1-80, 200-240, α3"}]
    assert text_of("domains", rows) == "N=res:A:1-80,res:200-240,α3"


def test_table_delta_is_applied_once():
    ss = {"swap": "α1,α3", "swap__ver": 0}
    ss["swap__table_0"] = {
        "edited_rows": {0: {"b": "α4"}},
        "added_rows": [{"a": "res:A:3-11", "b": "res:A:60-79"}],
        "deleted_rows": [],
    }
    apply_table(ss, "swap")
    assert ss["swap"] == "α1,α4; res:A:3-11,res:A:60-79" and ss["swap__ver"] == 1
    assert ss["last_changed"] == "swap"
    ss["swap__table_1"] = {"edited_rows": {}, "added_rows": [], "deleted_rows": [0]}
    apply_table(ss, "swap")
    assert ss["swap"] == "res:A:3-11,res:A:60-79" and ss["swap__ver"] == 2


def test_a_row_added_empty_stays_until_it_is_filled_in():
    """The table's add-row button sends an empty row at once; the text cannot hold it, so it waits in
    `<kind>__blank` (shown after the text's rows) until its elements are picked."""
    from foldmap.app_edits import table_rows

    ss = {"swap": "α1,α3", "swap__ver": 0, "swap__table_0": {"edited_rows": {}, "added_rows": [{}], "deleted_rows": []}}
    apply_table(ss, "swap")
    assert ss["swap"] == "α1,α3" and table_rows(ss, "swap") == [{"a": "α1", "b": "α3"}, {"a": "", "b": ""}]
    ss["swap__table_1"] = {"edited_rows": {1: {"a": "res:A:39-50"}}, "added_rows": [], "deleted_rows": []}
    apply_table(ss, "swap")  # one element of two: not written yet (it would be an error)
    assert ss["swap"] == "α1,α3" and ss["swap__blank"] == [{"a": "res:A:39-50", "b": ""}]
    ss["swap__table_2"] = {"edited_rows": {1: {"b": "res:A:60-79"}}, "added_rows": [{"b": "α2"}], "deleted_rows": []}
    apply_table(ss, "swap")
    assert ss["swap"] == "α1,α3; res:A:39-50,res:A:60-79" and ss["swap__blank"] == [{"a": "", "b": "α2"}]
    ss["swap__table_3"] = {"edited_rows": {}, "added_rows": [], "deleted_rows": [2]}
    apply_table(ss, "swap")
    assert table_rows(ss, "swap") == rows_of("swap", "α1,α3; res:A:39-50,res:A:60-79")

    ss = {"rename": "", "rename__ver": 0}  # a rename with no name yet would be `--rename needs REF=NAME`
    ss["rename__table_0"] = {"edited_rows": {}, "added_rows": [{"element": "res:A:3-11"}], "deleted_rows": []}
    apply_table(ss, "rename")
    assert ss["rename"] == "" and ss["rename__blank"] == [{"element": "res:A:3-11", "name": ""}]
    ss["rename__table_1"] = {"edited_rows": {0: {"name": "Gd"}}, "added_rows": [], "deleted_rows": []}
    apply_table(ss, "rename")
    assert ss["rename"] == "res:A:3-11=Gd" and ss["rename__blank"] == []


def test_text_typed_while_rows_wait_keeps_both():
    """Typing "as text" replaces the text's rows; the rows still waiting stay after them, and a table edit then
    indexes that same list."""
    from foldmap.app_edits import table_rows

    ss = {"swap": "α1,α3", "swap__ver": 0, "swap__blank": [{"a": "res:A:39-50", "b": ""}]}
    ss["swap"] = "α2,α4; α5,α6"  # typed
    assert table_rows(ss, "swap") == [*rows_of("swap", "α2,α4; α5,α6"), {"a": "res:A:39-50", "b": ""}]
    ss["swap__table_0"] = {"edited_rows": {2: {"b": "α7"}, 0: {"b": "α8"}}, "added_rows": [], "deleted_rows": []}
    apply_table(ss, "swap")
    assert ss["swap"] == "α2,α8; α5,α6; res:A:39-50,α7" and ss["swap__blank"] == []
    assert table_rows(ss, "swap") == rows_of("swap", ss["swap"])


@pytest.mark.parametrize(
    "kind, typed, rows",
    [
        ("swap", "α1,α3", [{"a": "res:A:3-11", "b": "res:A:60-79"}]),
        ("rename", "α1=Gd", [{"element": "res:A:3-11", "name": "Gd"}]),
        ("move", "α2=1,-2", [{"element": "res:A:39-50", "dx": 1.0, "dy": -2.0}]),
        ("mark", "α3=#2ca02c", [{"element": "res:A:60-79", "colour": "#2ca02c"}]),
    ],
)
def test_rows_from_pickers_draw_the_same_figure_as_typed_labels(kind, typed, rows):
    """Spec: table -> string -> layout_options -> the same figure as the text people typed (2LZM: α1 A:3-11,
    α2 A:39-50, α3 A:60-79)."""
    import json

    from foldmap import app_jobs

    def svg(text):
        return app_jobs.draw_job(DATA / "2LZM.cif", json.dumps({"theme": "publication", kind: text}))["svg"]

    assert svg(text_of(kind, rows)) == svg(typed)


def test_element_choices_name_the_drawn_element():
    elements = [
        {"ref": "res:A:23-34", "label": "α1", "chain": "A", "kind": "H", "first": "23", "last": "34"},
        {"ref": "res:B:23-34", "label": "α1", "chain": "B", "kind": "H", "first": "23", "last": "34"},
    ]
    assert element_choices(elements) == [("res:A:23-34", "α1 · A 23–34"), ("res:B:23-34", "α1 · B 23–34")]


def test_chain_choices_are_the_file_s_chains_unless_an_assembly_adds_copies():
    from foldmap.app_edits import chain_choices

    assert chain_choices(["1", "2", "3", "4"], ["3", "4", "1", "2"]) == ["1", "2", "3", "4"]  # file order
    assert chain_choices(["A", "B", "C", "D"], ["B", "D"]) == ["A", "B", "C", "D"]  # the rest can be asked for
    assert chain_choices(["A"], ["A1", "A2", "A3"]) == ["A1", "A2", "A3"]  # copies: named as drawn


def test_palette_is_hidden_only_where_it_colours_nothing():
    assert palette_used("chain", "black") and palette_used("shade", "black") and palette_used("domain", "black")
    assert not palette_used("plddt", "black") and not palette_used("sequence", "residue")
    assert palette_used("bfactor", "chain")


def test_the_failing_control_is_named():
    assert control_for("--swap needs two elements separated by a comma, got 'α1'", None) == "Swap (Edit layout)"
    assert control_for("no chain Z; the structure has protein chains A", None) == "Chains (Chains & assembly)"
    assert control_for("no element 'α9': ...", "rename") == "Rename (Edit layout)"
    assert control_for("something else", None) is None


def test_route_notes_read_as_one_line():
    assert notes_text([], []) == ""
    assert notes_text(["a>b"], ["c>d"]) == (
        "No clear route for a>b (drawn as plain curves). Drawn over a label or chain end, the only route found: c>d."
    )


def test_a_table_row_picks_its_structure():
    choices = [("AF-P07911-F1", "AlphaFold model"), ("4WRN", "X-ray")]
    assert row_source([1], choices) == "4WRN" and row_source([], choices) is None and row_source([5], choices) is None


# ------------------------------------------------------------- beyond the plan's list: display, symmetry, controls
def test_a_ref_not_drawn_is_shown_as_itself():
    from foldmap.app_edits import display

    choices = [("res:A:23-34", "α1 · A 23–34")]
    assert display("res:A:23-34", choices) == "α1 · A 23–34" and display("res:Q:1-9", choices) == "res:Q:1-9"


def test_symmetry_choices_are_the_ones_the_command_line_takes():
    from foldmap.app_edits import SYMMETRIES

    assert SYMMETRIES[:3] == ("auto", "off", "helical")
    assert SYMMETRIES[3:] == (*(f"C{n}" for n in range(2, 13)), *(f"D{n}" for n in range(2, 9)))


def test_every_setting_has_a_control_entry():
    from foldmap.app import HIDDEN_LAYOUT, LAYOUT_DEFAULTS, STYLE_FIELDS
    from foldmap.app_edits import CONTROLS, HIDDEN

    assert HIDDEN == HIDDEN_LAYOUT  # hard-coded in app_edits (importing app there would be a cycle)

    assert set(CONTROLS) == {*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme", "names"}
    for key, (label, section, widgets) in CONTROLS.items():
        assert (label == "" and section == "" and widgets == ()) == (key in HIDDEN_LAYOUT), key
    for key, (derived, _, _) in COMPOSITES.items():
        assert set(derived) <= set(CONTROLS[key][2]), key
    for key, (derived, _) in DERIVED.items():
        assert set(derived) <= set(CONTROLS[key][2]), key
    assert CONTROLS["swap"] == ("Swap", "Edit layout", ("swap__as_text", "swap"))
    assert CONTROLS["msa_path"] == ("Alignment", "Conservation & UniProt", ("msa_upload",))
    assert CONTROLS["names"] == ("Element names", "Edit layout", ("names__kw", "names_upload"))


def test_a_failing_layout_option_names_its_control():
    assert control_for("--move needs ref=dx,dy, got 'α2'", None) == "Move (Edit layout)"
    assert control_for("domain 'N': 'A:1-60' is no element; ...", None) == "Domains (Edit layout)"
    assert control_for("something else", "ligands") == "Ligands and ions (Show)"
    assert control_for("something else", "not_a_setting") is None


def _visible(text):
    """What a reader sees of Markdown-escaped text: the escapes resolved, the zero-width spaces gone."""
    return re.sub(r"\\(.)", r"\1", text).replace("\u200b", "")


def test_shown_text_names_no_path_carries_no_markdown_and_is_short():
    from foldmap.app_edits import shown

    kept = "/home/someone/.cache/foldmap/uploads/upload-0123456789ab.cif"
    text = shown(f"{kept}:1:0(0): expected block header", path=kept, name="bad")
    assert _visible(text) == "bad:1:0(0): expected block header"
    text = shown("cannot open /srv/app/.cache/foldmap/x.cif or ~/.cache/y.cif (chain(s) A)")
    assert _visible(text) == "cannot open x.cif or y.cif (chain(s) A)"
    crafted = shown("![x](https://example.invalid/p.png) [log in again](https://example.invalid) *b* `c` $x$ #h")
    assert "](" not in crafted and "![" not in crafted and "*b*" not in crafted and "$x$" not in crafted
    assert (
        _visible(shown("see https://github.com/dzyla/foldmap/issues."))
        == "see https://github.com/dzyla/foldmap/issues."
    )
    assert shown("x" * 1000, limit=50) == "x" * 49 + "…"
    assert _visible(shown("a/b and 1/2 and res:A:1-10 and α1")) == "a/b and 1/2 and res:A:1-10 and α1"
    assert _visible(shown("the alignment /nowhere/x.fasta", paths=False)) == "the alignment /nowhere/x.fasta"


def test_shown_text_cannot_autolink_show_icons_or_headings_and_names_no_file_url():
    """Bare URLs, www. and e-mail addresses are autolinked by GitHub-flavoured Markdown; a zero-width space breaks
    each. ':' is escaped (no :material/…: icons, no :red[…]), '=' and '-' too (no setext headings). file:// URLs and UNC
    paths are server paths like any other."""
    from foldmap.app_edits import shown

    for bare in ("http://evil.example/login", "www.evil.com", "a@b.com"):
        out = shown(f"open {bare} now")
        assert bare not in out and "\u200b" in out and _visible(out) == f"open {bare} now"
    assert ":material/home:" not in shown("go :material/home: there")
    assert "\n===" not in shown("Title\n===") and "\n---" not in shown("Title\n---")
    assert _visible(shown("read file:///home/dzyla/x.cif and file:///C:/Users/me/y.cif")) == "read x.cif and y.cif"
    assert _visible(shown("at \\\\server\\share\\z.cif")) == "at z.cif"


# ---------------------------------------------------------------------------------------------- names, domains
def test_names_keyword_is_derived_and_chosen_before_the_files():
    from foldmap import names
    from foldmap.app_edits import NAME_KEYWORDS, choose_names

    assert NAME_KEYWORDS == names.KEYWORDS  # hard-coded in app_edits (the pipeline stays out of it)
    ss = {"names": ["/x/a.yaml", "deposited"]}
    sync_composites(ss)
    assert ss["names__kw"] == "deposited"
    ss["names__kw"] = "beta"
    choose_names(ss)
    assert ss["names"] == ["beta", "/x/a.yaml"] and ss["last_changed"] == "names"
    ss["names__kw"] = "generic"
    choose_names(ss)
    assert ss["names"] == ["/x/a.yaml"]
    ss["names"] = []  # a layout file or a new structure
    sync_composites(ss)
    assert ss["names__kw"] == "generic"


def test_automatic_domains_keep_the_text_typed_before():
    from foldmap.app_edits import switch_domains

    ss = {"domains": "N=res:A:1-60"}
    sync_composites(ss)
    assert ss["domains__auto"] is False
    ss["domains__auto"] = True
    switch_domains(ss)
    assert ss["domains"] == "auto" and ss["domains__manual"] == "N=res:A:1-60" and ss["last_changed"] == "domains"
    ss["domains__auto"] = False
    switch_domains(ss)
    assert ss["domains"] == "N=res:A:1-60"
    ss["domains"] = "auto"  # a layout file
    sync_composites(ss)
    assert ss["domains__auto"] is True


def test_a_domain_row_waits_for_its_name_and_its_ranges():
    ss = {"domains": "", "domains__table_0": {"added_rows": [{"name": "N"}]}}
    apply_table(ss, "domains")
    assert ss["domains"] == "" and ss["domains__blank"] == [{"name": "N", "ranges": ""}]
    ss["domains__table_1"] = {"edited_rows": {0: {"ranges": "A:1-60"}}}
    apply_table(ss, "domains")
    assert ss["domains"] == "N=res:A:1-60" and ss["domains__blank"] == []


def test_domain_ranges_become_residue_ranges_only_where_they_map_exactly():
    rows = [{"name": "N", "ranges": "A:1-80, 200-240, 1:5-20, A1:3-7, α3, res:B:1-9"}]
    assert text_of("domains", rows) == "N=res:A:1-80,res:200-240,res:1:5-20,A1:3-7,α3,res:B:1-9"  # A1: a copy


def test_an_upload_name_cannot_carry_a_line_break():
    """A newline (or a Unicode line or paragraph separator, or a control character) in an uploaded file's name would end
    the YAML comment a hosted layout file starts with and write live keys below it; the kept name drops them."""
    from foldmap.app_edits import upload_name

    assert upload_name("x\nfoo: 1\r\n.yaml") == "xfoo: 1.yaml"
    assert upload_name("a b c\x00d​e.yaml") == "abcde.yaml"
    assert upload_name("ü names #1: a.yaml") == "ü names #1: a.yaml"
    assert upload_name("\n\n") == "upload" and len(upload_name("n" * 500)) == 120
