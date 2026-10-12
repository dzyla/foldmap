"""The sidebar: a short essentials list, then sections for everything else."""

import pytest

pytest.importorskip("streamlit")

from app_helpers import DATA, count_calls, fresh_caches, opened, svg_of, widget_keys  # noqa: E402


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
    for key in ("rename", "domains", "mark", "backbone_colors"):
        (at.text_area(key=key) if key == "domains" else at.text_input(key=key)).set_value(kept[key])
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
