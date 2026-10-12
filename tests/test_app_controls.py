"""The sidebar: a short essentials list, then sections for everything else."""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")

from app_helpers import DATA, as_text, count_calls, fresh_caches, opened, svg_of, visible, widget_keys  # noqa: E402


def test_first_view_shows_at_most_twelve_controls():
    at = opened(DATA / "1UBQ.cif")
    side, main = widget_keys(at.sidebar), widget_keys(at.main)
    assert side == [
        "theme", "theme_browse", "color_by", "palette", "labels", "residue_numbers", "legend",
        "font_scale", "title", "reset_theme",
    ]  # fmt: skip
    assert main == ["pane", "new_structure"] and len(side) + len(main) <= 12  # the view switch is `pane` (R1)


def test_palette_gives_way_where_it_colours_nothing():
    at = opened(DATA / "1UBQ.cif")
    at.selectbox(key="color_by").select("sequence").run()
    keys = widget_keys(at.sidebar)
    assert "sequence_map" in keys and "palette" not in keys


def test_theme_dialog_applies_a_theme_and_draws_nothing_new(monkeypatch):
    import foldmap.spec as spec
    from foldmap.app import PREVIEW_THEMES

    at = opened(DATA / "1LMB.cif")
    calls = count_calls(monkeypatch, spec, "draw_spec")
    at.button(key="theme_browse").click().run()
    assert len(at.image) == len(PREVIEW_THEMES) and calls == []  # thumbnails are static files
    assert any(c.value.startswith("Shown on ubiquitin, the reference structure") for c in at.caption)
    at.button(key="pick_trace").click().run()
    assert not at.exception and at.selectbox(key="theme").value == "trace"
    assert at.session_state["loop_arrows"] is True and not at.image and "<svg" in svg_of(at)


def test_reset_restores_the_theme():
    from foldmap.style import THEMES

    at = opened(DATA / "1UBQ.cif")
    at.selectbox(key="palette").select("tol-bright").run()
    at.toggle(key="legend").set_value(False).run()
    at.button(key="reset_theme").click().run()
    assert at.session_state["palette"] == THEMES["publication"].palette and at.session_state["legend"] is True


def test_controls_change_the_figure():
    """Was in test_app.py (its residue-numbers half; the mode half is in Task 9)."""
    at = opened(DATA / "1UBQ.cif")
    assert "resnum:" not in svg_of(at)
    at.toggle(key="residue_numbers").set_value(True).run()
    assert "resnum:" in svg_of(at)


def test_reset_keeps_what_names_this_structure():
    """Reset to theme restores the look only: chains, assembly, edits, domains, marks, highlight, H-bonds and
    per-chain colours name this structure's parts and stay."""
    from foldmap.style import THEMES

    at = opened(DATA / "2HHB.cif")
    kept = {
        "chains": "A,B",
        "assembly": "asu",
        "rename": "res:A:1-10=Zz",
        "domains": "N=res:A:1-60",
        "mark": "res:A:1-10=#2ca02c",
        "highlight": "A",
        "backbone_colors": "A=#1f77b4",
    }
    for key in ("rename", "mark"):
        as_text(at, key, kept[key])
    as_text(at, "domains", kept["domains"])
    at.text_input(key="backbone_colors").set_value(kept["backbone_colors"])
    at.multiselect(key="chains__pick").set_value(["A", "B"])
    at.selectbox(key="assembly").select("asu")
    at.segmented_control(key="highlight__mode").set_value("chains").run()
    at.multiselect(key="highlight__chains").set_value(["A"])
    at.segmented_control(key="hbonds__mode").set_value("interchain").run()
    at.toggle(key="legend").set_value(False).run()
    assert not at.exception
    at.button(key="reset_theme").click().run()
    ss = at.session_state
    assert not at.exception and ss["legend"] is THEMES["publication"].legend and ss["last_changed"] == "theme"
    assert {k: ss[k] for k in kept} == kept and ss["hbonds"] == "interchain"
    assert at.segmented_control(key="hbonds__mode").value == "interchain"


def test_style_choices_are_segmented_with_the_style_values():
    from foldmap.style import CHOICES

    at = opened(DATA / "1UBQ.cif")
    for key in ("fill", "helix_shading", "helix_angle", "loops"):
        assert list(at.segmented_control(key=key).options) == list(CHOICES[key]), key
    at.segmented_control(key="loops").set_value("curved").run()
    assert at.session_state["loops"] == "curved" and not at.exception


def test_nucleic_acid_block_only_with_nucleic_acid():
    assert "na_style" not in widget_keys(opened(DATA / "1UBQ.cif").sidebar, skip=())
    keys = widget_keys(opened(DATA / "1LMB.cif").sidebar, skip=())
    assert {"na_style", "na_extent", "na_backbone__mode", "dna_color", "rna_color", "check_contacts"} <= set(keys)


def test_backbone_colour_and_hbonds_reach_the_figure():
    at = opened(DATA / "1LMB.cif")
    at.segmented_control(key="na_backbone__mode").set_value("colour").run()
    at.color_picker(key="na_backbone__colour").pick("#123456").run()
    assert at.session_state["na_backbone"] == "#123456"
    at.segmented_control(key="hbonds__mode").set_value("na").run()
    assert at.session_state["hbonds"] == "na" and 'id="hbond:' in svg_of(at)


def test_named_colours_reach_the_pickers_as_hex():
    """A layout file may name colours (teal); the pickers take #rrggbb."""
    at = opened(DATA / "1LMB.cif")
    text = b"foldmap: 1\ntheme: publication\nstyle: {na_backbone: teal, background: ivory}\n"
    at.file_uploader(key="layout_upload").set_value(("layout.yaml", text, "application/x-yaml")).run()
    assert not at.exception and at.segmented_control(key="na_backbone__mode").value == "colour"
    assert at.color_picker(key="na_backbone__colour").value == "#008080" and at.session_state["na_backbone"] == "teal"
    assert at.color_picker(key="background").value == "#fffff0" and "background:#fffff0" in svg_of(at)


def test_ligand_codes_reach_the_figure():
    at = opened(DATA / "2HHB.cif")
    at.segmented_control(key="ligands__mode").set_value("codes").run()
    assert at.session_state["ligands"] == "auto" and at.segmented_control(key="ligands__mode").value == "codes"
    at.text_input(key="ligands__codes").input("hem").run()
    assert at.session_state["ligands"] == "HEM" and not at.exception


def test_contacts_check_on_request():
    at = opened(DATA / "1LMB.cif")
    assert not any(c.value.startswith("Contacts checked") for c in at.caption)
    at.toggle(key="check_contacts").set_value(True).run()
    assert any(c.value == "Contacts checked against the coordinates: 21 of 21 confirmed." for c in at.caption)


def test_contacts_check_lists_at_most_ten_problems(monkeypatch):
    from foldmap import app_jobs

    report = {"drawn": 30, "verified": 18, "problems": [f"missing: DA{i} <- LYS{i}" for i in range(12)]}
    monkeypatch.setattr(app_jobs, "contacts_job", lambda path, blob: report)
    fresh_caches()
    try:
        at = opened(DATA / "1LMB.cif")
        at.toggle(key="check_contacts").set_value(True).run()
    finally:
        fresh_caches()  # the made-up report stays out of the other tests
    texts = [c.value for c in at.caption]
    assert "Contacts checked against the coordinates: 18 of 30 confirmed." in texts
    assert sum(t.startswith("missing\\:") for t in texts) == 10 and "… and 2 more" in texts


def test_contacts_check_that_fails_keeps_the_figure(monkeypatch):
    from foldmap import app_jobs

    def broken(path, blob):
        raise ValueError("no beads to read")

    monkeypatch.setattr(app_jobs, "contacts_job", broken)
    fresh_caches()
    try:
        at = opened(DATA / "1LMB.cif")
        at.toggle(key="check_contacts").set_value(True).run()
    finally:
        fresh_caches()
    assert not at.exception and "<svg" in svg_of(at)
    assert [w.value for w in at.warning] == ["Contacts check: no beads to read"]


def test_chains_are_picked_not_typed():
    at = opened(DATA / "2HHB.cif")
    pick = at.multiselect(key="chains__pick")
    assert pick.options == ["A", "B", "C", "D"] and pick.value == []
    pick.set_value(["B", "D"]).run()
    assert not at.exception and at.session_state["chains"] == "B,D"
    at.toggle(key="chains__as_text").set_value(True).run()
    assert at.text_input(key="chains").value == "B,D"


def test_chain_picker_keeps_every_chain_while_some_are_drawn():
    """The figure of chains B, D names only B and D; the picker still offers all four."""
    at = opened(DATA / "2HHB.cif")
    at.multiselect(key="chains__pick").set_value(["B", "D"]).run()
    assert at.multiselect(key="chains__pick").options == ["A", "B", "C", "D"]
    assert at.session_state["protein_chains"] == ["B", "D"]  # what the highlight and Sequence pickers offer


def test_chain_picker_offers_the_structure_s_chains():
    """From the file (the size check), not from a figure: nucleic-acid strands too, in the file's order."""
    at = opened(DATA / "1LMB.cif")
    assert at.multiselect(key="chains__pick").options == ["1", "2", "3", "4"]
    at.selectbox(key="assembly").select("asu").run()
    assert at.multiselect(key="chains__pick").options == ["1", "2", "3", "4"]


def test_top_bar_says_how_many_chains_are_drawn():
    at = opened(DATA / "2HHB.cif")
    summary = lambda: at.session_state["info"]  # noqa: E731
    assert summary().startswith("4 chains · ")
    at.multiselect(key="chains__pick").set_value(["B", "D"]).run()
    assert not at.exception and summary().startswith("2 of 4 chains · ")
    assert any(c.value.startswith("2 of 4 chains · ") for c in at.caption)
    at.multiselect(key="chains__pick").set_value([]).run()
    assert summary().startswith("4 chains · ")


def test_chain_count_of_the_summary_follows_the_assembly():
    """ "n of N chains" counts N under the assembly the figure is drawn in; until a figure of every chain has been
    drawn in a new assembly, the summary gives the count drawn only."""
    at = opened(DATA / "2HHB.cif")
    at.multiselect(key="chains__pick").set_value(["B", "D"]).run()
    assert at.session_state["info"].startswith("2 of 4 chains · ")
    at.selectbox(key="assembly").select("asu").run()
    assert not at.exception and at.session_state["info"].startswith("2 chains · ")  # N not counted in asu yet
    at.multiselect(key="chains__pick").set_value([]).run()
    assert at.session_state["chains_total"] == {"assembly": "asu", "count": 4}
    at.multiselect(key="chains__pick").set_value(["B"]).run()
    assert at.session_state["info"].startswith("1 of 4 chains · ")


def test_long_chain_and_assembly_ids_are_cut_in_the_pickers():
    from foldmap.app_controls import capped

    assert capped("A") == "A" and capped("x" * 24) == "x" * 24 and capped("x" * 40) == "x" * 23 + "…"
    at = opened(DATA / "2HHB.cif")
    up = ("layout.yaml", f"foldmap: 1\nlayout: {{assembly: {'x' * 40}}}\n".encode(), "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    assert not at.exception and at.session_state["assembly"] == "x" * 40
    assert at.selectbox(key="assembly").options[-1] == "x" * 23 + "…"


def test_highlight_symmetry_assembly_and_focus():
    from foldmap.app_edits import SYMMETRIES

    at = opened(DATA / "2HHB.cif")
    at.segmented_control(key="highlight__mode").set_value("chains").run()
    at.multiselect(key="highlight__chains").set_value(["A"]).run()
    assert at.session_state["highlight"] == "A" and not at.exception
    assert list(at.selectbox(key="symmetry").options) == list(SYMMETRIES)
    assert at.selectbox(key="assembly").options[:2] == ["auto", "asu"]
    at.segmented_control(key="focus__mode").set_value("none").run()
    assert at.session_state["focus"] == "none"


def test_highlighted_chains_grey_the_rest_and_read_back_as_text():
    from foldmap.render import MATE_GREY

    at = opened(DATA / "2HHB.cif")
    assert MATE_GREY not in svg_of(at)
    at.segmented_control(key="highlight__mode").set_value("chains").run()
    at.multiselect(key="highlight__chains").set_value(["A", "C"]).run()
    assert not at.exception and MATE_GREY in svg_of(at)
    at.toggle(key="highlight__as_text").set_value(True).run()
    assert at.text_input(key="highlight").value == "A,C"
    at.text_input(key="highlight").set_value("protomer").run()
    at.toggle(key="highlight__as_text").set_value(False).run()
    assert not at.exception and at.segmented_control(key="highlight__mode").value == "protomer"


def test_contact_lines():
    from foldmap.app_controls import contact_lines

    assert contact_lines(None) == ["No contact beads in this figure to check."]
    report = {"drawn": 3, "verified": 2, "problems": ["missing: DA1 <- LYS2"]}
    assert contact_lines(report) == [
        "Contacts checked against the coordinates: 2 of 3 confirmed.",
        "missing: DA1 <- LYS2",
    ]


def test_contacts_check_reruns_only_for_what_changes_the_contacts(monkeypatch):
    """The contacts check is cached on the settings that change the contacts drawn (chains, assembly, what DNA is
    shown, ...): a palette change redraws the figure but does not check the contacts again."""
    from foldmap import app_jobs

    fresh_caches()
    try:
        at = opened(DATA / "1LMB.cif")
        calls = count_calls(monkeypatch, app_jobs, "contacts_job")
        at.toggle(key="check_contacts").set_value(True).run()
        assert calls == [1] and at.session_state["last_changed"] == "check_contacts"
        before = svg_of(at)
        at.selectbox(key="palette").select("tol-bright").run()
        assert not at.exception and svg_of(at) != before and calls == [1]
        at.segmented_control(key="na_extent").set_value("all").run()
        assert not at.exception and calls == [1, 1]
    finally:
        fresh_caches()


def test_protein_dna_hbonds_offered_only_with_nucleic_acid():
    from foldmap.app_controls import HBONDS

    shown = opened(DATA / "2HHB.cif").segmented_control(key="hbonds__mode").options  # the labels
    assert list(shown) == [HBONDS[k] for k in ("off", "interchain", "chains")]
    assert list(opened(DATA / "1LMB.cif").segmented_control(key="hbonds__mode").options) == list(HBONDS.values())


# ---------------------------------------------------------------------------------------------- Edit layout
def test_layout_mode_and_sizes():
    at = opened(DATA / "1UBQ.cif")
    at.segmented_control(key="mode").set_value("stack").run()
    assert not at.exception and "<svg" in svg_of(at)
    # AppTest lists formatted labels; values stay the Style values
    assert list(at.segmented_control(key="fold_template").options) == ["βα row", "as projected"]
    assert list(at.segmented_control(key="helix_on_sheet").options) == ["beside the sheet", "over the sheet"]
    assert at.segmented_control(key="fold_template").value == "auto"
    at.segmented_control(key="helix_on_sheet").set_value("auto").run()
    assert at.session_state["helix_on_sheet"] == "auto" and not at.exception


def test_helices_stay_beside_the_sheet_unless_asked():
    at = opened(DATA / "1UBQ.cif")
    assert at.session_state["helix_on_sheet"] == "off" and at.segmented_control(key="helix_on_sheet").value == "off"


def test_a_layout_file_drops_the_rows_waiting_in_the_tables():
    at = opened(DATA / "2LZM.cif")
    at.session_state["swap__blank"] = [{"a": "res:A:3-11", "b": ""}]  # a swap row half made in the table
    up = ("layout.yaml", b"foldmap: 1\nedits: {swap: [[res:A:39-50, res:A:60-79]]}\n", "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    assert not at.exception and "swap__blank" not in at.session_state
    assert at.session_state["swap"] == "res:A:39-50,res:A:60-79"


def test_layout_section_holds_the_editors_and_dna_size_only_with_nucleic_acid():
    keys = widget_keys(opened(DATA / "1UBQ.cif").sidebar, skip=())
    assert {"mode", "rotate", "helix_scale", "strand_scale", "fold_template", "helix_on_sheet"} <= set(keys)
    assert {f"{kind}__as_text" for kind in ("swap", "move", "rename", "mark")} <= set(keys)
    assert "dna_scale" not in keys and "dna_scale" in widget_keys(opened(DATA / "1LMB.cif").sidebar, skip=())


def test_editor_text_and_tables_write_the_same_string():
    at = opened(DATA / "2LZM.cif")
    as_text(at, "swap", "res:A:3-11,res:A:60-79")
    assert at.session_state["swap"] == "res:A:3-11,res:A:60-79"
    assert any(visible(c.value) == "Writes: res:A:3-11,res:A:60-79" for c in at.caption)  # shown(): escaped
    at.toggle(key="swap__as_text").set_value(False).run()  # back to the table: same value
    assert at.session_state["swap"] == "res:A:3-11,res:A:60-79" and not at.exception


def test_element_pickers_list_the_drawn_elements():
    at = opened(DATA / "2LZM.cif")
    els = at.session_state["elements"]
    assert els[0]["ref"] == "res:A:3-11" and els[0]["label"] == "α1"


def test_editor_tables_show_the_rows_of_the_text():
    from foldmap.app_controls import editor_frame, element_options

    frame = editor_frame("move", [{"element": "α2", "dx": 1.0, "dy": ""}, {"element": "", "dx": None, "dy": None}])
    assert list(frame.columns) == ["element", "dx", "dy"] and frame["dx"].dtype == float
    assert frame["dx"].tolist()[0] == 1.0 and frame["dy"].isna().all()  # a cell that is not a number: empty
    assert list(editor_frame("swap", []).columns) == ["a", "b"]
    choices = [("res:A:3-11", "α1 · A 3–11"), ("res:A:60-79", "α3 · A 60–79")]
    rows = [{"a": "α1", "b": "res:A:60-79"}]
    assert element_options("swap", rows, choices) == ["res:A:3-11", "res:A:60-79", "α1"]  # typed labels kept


# ---------------------------------------------------------------------------------------------- names and domains
GLOBIN = Path(__file__).parents[1] / "docs" / "names" / "2HHB-globin.yaml"


def test_conventional_names_are_a_choice():
    at = opened(DATA / "2HHB.cif")
    generic = svg_of(at)
    at.segmented_control(key="names__kw").set_value("deposited").run()
    assert at.session_state["names"] == ["deposited"] and svg_of(at) != generic and not at.exception


def test_a_names_file_can_be_uploaded():
    from pathlib import Path

    yaml_file = Path(__file__).parents[1] / "docs" / "names" / "2HHB-globin.yaml"
    at = opened(DATA / "2HHB.cif")
    generic = svg_of(at)
    at.file_uploader(key="names_upload").set_value(("2HHB-globin.yaml", yaml_file.read_bytes(), "text/yaml")).run()
    assert not at.exception and at.session_state["names"][-1].endswith(".yaml") and svg_of(at) != generic


def test_domains_automatic_and_by_table():
    at = opened(DATA / "2LZM.cif")
    at.toggle(key="domains__auto").set_value(True).run()
    assert at.session_state["domains"] == "auto" and not at.exception
    at.toggle(key="domains__auto").set_value(False).run()
    as_text(at, "domains", "N=res:A:1-60\nC=res:A:61-164")
    assert at.session_state["domains"] == "N=res:A:1-60\nC=res:A:61-164"
    assert 'id="domain-bar:' in svg_of(at)
    at.segmented_control(key="domain_bar").set_value("off").run()
    assert 'id="domain-bar:' not in svg_of(at)


def test_a_names_file_keeps_the_keyword_shows_its_name_and_can_be_removed():
    at = opened(DATA / "2HHB.cif")
    at.segmented_control(key="names__kw").set_value("beta").run()
    at.file_uploader(key="names_upload").set_value(("my globin.yaml", GLOBIN.read_bytes(), "text/yaml")).run()
    names = at.session_state["names"]
    assert names[0] == "beta" and Path(names[1]).is_file() and len(names) == 2 and not at.exception
    assert "names_upload" not in widget_keys(at.sidebar, skip=())  # a file is set: its name and Remove instead
    assert any(visible(c.value) == "Names file: my globin.yaml" for c in at.caption)
    at.button(key="names_clear").click().run()
    assert at.session_state["names"] == ["beta"] and "names_upload" in widget_keys(at.sidebar, skip=())
    at.segmented_control(key="names__kw").set_value("generic").run()
    assert at.session_state["names"] == [] and not at.exception


@pytest.mark.parametrize(
    "data, said",
    [
        (b"- A\n- B\n", "a names file is a mapping with the keys elements and stems"),
        (b"elements: [A, B\n", "not valid YAML"),
        (b"helices:\n  res:A:3-18: A\n", "unknown key 'helices'"),
        (b"elements:\n  res:A:3-18: [A, B]\n", "'res:A:3-18' is not a name"),
        (
            b"elements:\n  r1: &a [x, x, x]\n  r2: &b [*a, *a, *a]\n  r3: [*b, *b, *b]\n",
            "aliases are not allowed in a names file",
        ),
        (b"#" * 200_001, "at most 100 kB"),
        (b"\xff\xfe", "UTF-8"),
        (b"elements: " + b"[" * 2000 + b"]" * 2000 + b"\n", "nested too deeply"),  # 4 kB, no RecursionError
        (b"elements:\n  res:A:3-18: 2001-01-01\n", "reads as a date; give a name as text"),
    ],
)
def test_an_upload_that_is_not_a_names_file_is_an_error_and_changes_nothing(data, said):
    at = opened(DATA / "2HHB.cif")
    at.segmented_control(key="names__kw").set_value("deposited").run()
    before = svg_of(at)
    at.file_uploader(key="names_upload").set_value(("bad.yaml", data, "text/yaml")).run()
    assert not at.exception and at.session_state["names"] == ["deposited"] and svg_of(at) == before
    assert any(said in visible(e.value) for e in at.error), [e.value for e in at.error]


def test_a_names_file_s_keys_and_file_name_are_shown_as_plain_text(monkeypatch, tmp_path):
    from app_helpers import rendered_text, seed_cache

    monkeypatch.setenv("FOLDMAP_APP_LOCAL_FILES", "0")
    seed_cache(monkeypatch, tmp_path, "2HHB")
    at = opened("2HHB")
    up = b"elements:\n  '[x](https://evil.example) ![i](https://evil.example/i.png)': A\n"
    at.file_uploader(key="names_upload").set_value(("[n](https://evil.example).yaml", up, "text/yaml")).run()
    text = rendered_text(at)
    assert not at.exception and at.error and "](" not in text and "uploads" not in text


def test_two_elements_of_a_chain_named_alike_is_an_error_naming_the_file():
    """Duplicate names within a chain are the names module's error, shown with the uploaded file's own name."""
    at = opened(DATA / "2HHB.cif")
    before = svg_of(at)
    dup = b"elements:\n  res:A:3-18: X\n  res:A:20-35: X\n"
    at.file_uploader(key="names_upload").set_value(("dup.yaml", dup, "text/yaml")).run()
    assert not at.exception and svg_of(at) == before
    text = visible(at.error[0].value)
    assert text.startswith("dup.yaml:") and "one name, one element" in text and "upload-" not in text
    assert text.endswith("Check Element names (Edit layout).")


def test_names_and_domains_controls_follow_a_layout_file_and_the_structure():
    at = opened(DATA / "2LZM.cif")
    at.segmented_control(key="names__kw").set_value("deposited").run()
    at.toggle(key="domains__auto").set_value(True).run()
    up = ("layout.yaml", b"foldmap: 1\nnames: [beta]\ndomains: {N: [res:A:1-60]}\n", "application/x-yaml")
    at.file_uploader(key="layout_upload").set_value(up).run()
    ss = at.session_state
    assert not at.exception and ss["names"] == ["beta"] and ss["names__kw"] == "beta"
    assert ss["domains"] == "N=res:A:1-60" and ss["domains__auto"] is False
    at.toggle(key="domains__auto").set_value(True).run()
    at.button(key="new_structure").click().run()
    at.text_input(key="source").input(str(DATA / "1UBQ.cif"))
    at.button(key="load").click().run()
    ss = at.session_state
    assert not at.exception and ss["names"] == [] and ss["names__kw"] == "generic"
    assert ss["domains"] == "" and ss["domains__auto"] is False and "domains__manual" not in ss


def test_domain_colours_and_bar_only_with_domains():
    at = opened(DATA / "2LZM.cif")
    keys = widget_keys(at.sidebar, skip=())
    assert {"names__kw", "names_upload", "domains__auto", "domains__as_text"} <= set(keys)
    assert "domain_bar" not in keys and "domain_colours" not in keys  # nothing to colour or bar without domains
    at.toggle(key="domains__auto").set_value(True).run()
    keys = widget_keys(at.sidebar, skip=())
    assert {"domain_bar", "domain_colours"} <= set(keys) and "domains__as_text" not in keys


def test_every_setting_is_reachable():
    """Every figure setting has a control in the sidebar (expanders render their contents): on 1LMB with every
    "as text" switch off and every composite in its widest mode, then on 1UBQ with an alignment."""
    from foldmap.app import HIDDEN_LAYOUT, LAYOUT_DEFAULTS, STYLE_FIELDS
    from foldmap.app_edits import CONTROLS

    def shown(at):  # expanders render their contents, so skip nothing
        return set(widget_keys(at.sidebar, skip=()))

    at = opened(DATA / "1LMB.cif")  # nucleic acid: the DNA/RNA block is shown
    seen = shown(at)
    for key, value in (("hbonds__mode", "chains"), ("na_backbone__mode", "modality"), ("ligands__mode", "codes")):
        at.segmented_control(key=key).set_value(value).run()
        seen |= shown(at)
    at.selectbox(key="color_by").select("sequence").run()  # sequence_map instead of palette
    seen |= shown(at)
    at.toggle(key="domains__auto").set_value(True).run()  # with domains: colour by domain and the domain bar
    seen |= shown(at)
    at = opened(DATA / "1UBQ.cif")  # an alignment: the reference row picker is shown
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    seen |= shown(at)
    wanted = {*STYLE_FIELDS, *LAYOUT_DEFAULTS, "theme"} - set(HIDDEN_LAYOUT)
    missing = sorted(k for k in wanted if not seen & set(CONTROLS[k][2]))
    assert missing == []


def test_alignment_reference_is_picked_from_the_rows():
    at = opened(DATA / "1UBQ.cif")
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    ref = at.selectbox(key="msa_reference")
    assert ref.options[0] == "best match" and "yeast_UBI4" in ref.options  # formatted labels
    ref.select("yeast_UBI4").run()
    assert at.session_state["msa_reference"] == "yeast_UBI4" and not at.exception


def test_uniprot_controls_write_the_option(monkeypatch):
    import json

    from foldmap import uniprot

    e = uniprot.parse_entry(json.loads((DATA / "uniprot" / "P07911.json").read_text()))
    monkeypatch.setattr(uniprot, "entry", lambda acc, cache_dir=None: e)
    at = opened(DATA / "1UBQ.cif")
    at.segmented_control(key="uniprot__mode").set_value("accession").run()
    at.text_input(key="uniprot__acc").input("p07911").run()
    assert at.session_state["uniprot"] == "P07911" and not at.exception
    at.segmented_control(key="uniprot_sites__mode").set_value("types").run()
    at.multiselect(key="uniprot_sites__types").set_value(["active", "binding"]).run()
    assert at.session_state["uniprot_sites"] == "active,binding"


def test_the_sidebar_is_the_essentials_then_six_sections():
    from foldmap import CREDIT

    at = opened(DATA / "1UBQ.cif")
    labels = [e.label for e in at.sidebar.expander]
    assert labels == ["Style", "Show", "Chains & assembly", "Edit layout", "Conservation & UniProt", "Layout file"]
    assert at.sidebar.caption[-1].value == CREDIT  # the credit closes it


def test_the_editors_are_tabs_in_edit_layout():
    """One editor in view at a time (the all-open sidebar stays short); every editor is still built, so each keeps
    its table and "as text" switch."""
    at = opened(DATA / "2LZM.cif")
    assert [t.label for t in at.sidebar.tabs] == ["Swap", "Move", "Rename", "Mark", "Domains"] and not at.main.tabs
    assert all(at.toggle(key=f"{kind}__as_text") for kind in ("swap", "move", "rename", "mark", "domains"))


def test_an_editor_tab_shows_no_second_heading_but_keeps_its_help_and_switch():
    from foldmap.app_controls import EDITOR_HELP, EDITOR_TABS

    at = opened(DATA / "2LZM.cif")
    assert not [m.value for m in at.sidebar.markdown if m.value in {f"**{label}**" for label in EDITOR_TABS}]
    for kind in ("swap", "move", "rename", "mark", "domains"):
        assert at.toggle(key=f"{kind}__as_text").proto.help == EDITOR_HELP[kind]


def test_option_labels_and_help_say_everything_in_full():
    from foldmap.app_controls import ALIGNMENT_TYPES, CONSERVATION_NOTE, HBONDS

    at = opened(DATA / "1UBQ.cif")
    assert HBONDS["chains"] == "chain pairs" and "ConSurf grades file" in CONSERVATION_NOTE
    assert all(f".{t}" in at.file_uploader(key="msa_upload").proto.help for t in ALIGNMENT_TYPES)
    assert list(at.file_uploader(key="msa_upload").proto.type) == [f".{t}" for t in ALIGNMENT_TYPES]


def test_credit_stays_visible_in_the_workspace():
    from foldmap import CREDIT

    at = opened(DATA / "1UBQ.cif")
    assert any(c.value == CREDIT for c in at.sidebar.caption)


def test_a_consurf_table_has_no_reference_row_and_clears_the_old_one():
    at = opened(DATA / "1UBQ.cif")
    at.file_uploader(key="msa_upload").set_value(
        ("ubq.fasta", (DATA / "ubq_family.fasta").read_bytes(), "text/plain")
    ).run()
    at.selectbox(key="msa_reference").select("yeast_UBI4").run()
    table = DATA / "consurf" / "1TIM_A_grades.txt"
    at.file_uploader(key="msa_upload").set_value(("1TIM_A_grades.txt", table.read_bytes(), "text/plain")).run()
    ss = at.session_state
    assert not at.exception and ss["msa_path"].endswith(".txt") and ss["msa_reference"] == ""
    assert "msa_reference" not in widget_keys(at.sidebar, skip=())
