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
    for key, value in kept.items():
        (at.text_area(key=key) if key == "domains" else at.text_input(key=key)).set_value(value)
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


def test_contact_lines():
    from foldmap.app_controls import contact_lines

    assert contact_lines(None) == ["No contact beads in this figure to check."]
    report = {"drawn": 3, "verified": 2, "problems": ["missing: DA1 <- LYS2"]}
    assert contact_lines(report) == [
        "Contacts checked against the coordinates: 2 of 3 confirmed.",
        "missing: DA1 <- LYS2",
    ]
