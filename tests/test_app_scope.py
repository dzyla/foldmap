"""The bar above the figure that says when an automatic choice drew part of the file, and draws all of it on request."""

import pytest

pytest.importorskip("streamlit")

from app_helpers import DATA, load, opened, rendered_text, visible, widget_keys  # noqa: E402

KELCH = DATA / "4CH9.cif"  # assembly 1 (B, D) is one of two copies of the complex in the file
APO = DATA / "assemblies" / "7A4M.cif"  # one subunit of a 24-subunit cage drawn by focus auto


def bar(at) -> list[str]:
    return [k for k in widget_keys(at.main) if k in ("scope_all", "scope_keep")]


def test_the_bar_says_what_was_drawn_and_show_all_draws_every_chain():
    at = opened(KELCH)
    assert bar(at) == ["scope_all", "scope_keep"]
    assert "Drawn: assembly 1, one of 2 copies in the file — 2 of 4 chains." in visible(rendered_text(at))
    assert at.session_state["info"].startswith("2 of 4 chains · ")  # the top bar counts them too
    at.button(key="scope_all").click().run()
    assert not at.exception
    assert at.session_state["assembly"] == "asu" and at.session_state["focus"] == "none"
    assert len(at.session_state["protein_chains"]) == 4 and bar(at) == []
    assert at.session_state["info"].startswith("4 chains · ")
    assert at.selectbox(key="assembly").value == "asu"  # the sidebar shows what drew it


def test_show_all_redraws_when_auto_redraw_is_off():
    at = opened(KELCH)
    at.session_state["auto_redraw"] = False
    at.run()
    at.button(key="scope_all").click().run()
    assert not at.exception and len(at.session_state["protein_chains"]) == 4 and bar(at) == []


def test_keep_this_hides_the_bar_until_the_next_structure():
    at = opened(APO)
    assert bar(at) == ["scope_all", "scope_keep"]
    assert "Drawn: one of 24 subunits of an O cage — 1 of 24 chains." in visible(rendered_text(at))
    at.button(key="scope_keep").click().run()
    assert not at.exception and bar(at) == [] and at.session_state["focus"] == "auto"
    assert "Drawn:" not in rendered_text(at)
    at.button(key="new_structure").click().run()
    load(at, KELCH)
    assert not at.exception and bar(at) == ["scope_all", "scope_keep"]


@pytest.mark.parametrize("name", ["na/1KX5.cif.gz", "assemblies/3UTB.cif.gz"])
def test_nothing_left_out_shows_no_bar(name):
    at = opened(DATA / name)
    assert bar(at) == [] and "Drawn:" not in rendered_text(at)


def test_everything_too_large_for_the_app_is_said_not_offered(monkeypatch):
    """Drawing all 24 apoferritin chains is over the residue cap: the bar says so and how to draw it elsewhere."""
    monkeypatch.setenv("FOLDMAP_APP_MAX_RESIDUES", "1000")  # one subunit (~170) loads, 24 do not fit
    at = opened(APO)
    assert bar(at) == ["scope_keep"]
    text = visible(rendered_text(at))
    assert "All 24 chains: too large for this app — draw with foldmap plot 7A4M --focus none" in text


def test_the_bar_reads_in_both_themes():
    """The bar is the theme's tint with the text colour (the stylesheet, light-dark pairs): text 7:1, the too-large
    caption 4.5:1, the button edges 3:1, in both themes."""
    from foldmap.app import APP_THEME, MUTED, STYLESHEET, TINT
    from foldmap.render import contrast_ratio

    assert f"light-dark({TINT['light']}, {TINT['dark']})" in STYLESHEET and ".st-key-scope_bar" in STYLESHEET
    for mode in ("light", "dark"):
        t = APP_THEME[f"theme.{mode}"]
        assert contrast_ratio(t["textColor"], TINT[mode]) >= 7
        assert contrast_ratio(MUTED[mode], TINT[mode]) >= 4.5
        assert contrast_ratio(t["borderColor"], TINT[mode]) >= 3


def test_the_bar_shows_a_file_s_assembly_name_as_plain_text(tmp_path):
    """The assembly name comes from the file: it can neither link nor style anything (app_edits.shown)."""
    import gemmi

    st = gemmi.read_structure(str(KELCH))
    st.assemblies[0].name = "[x](https://e.example)*1*"
    path = tmp_path / "kelch.cif"
    st.make_mmcif_document().write_file(str(path))
    at = opened(path)
    line = next(m.value for m in at.markdown if visible(m.value).startswith("Drawn:"))
    assert "](" not in line and "\\*1\\*" in line
    assert visible(line) == "Drawn: assembly [x](https://e.example)*1*, one of 2 copies in the file — 2 of 4 chains."


def test_scoped_summary_counts_the_chains_drawn_of_the_file():
    from foldmap.app_jobs import scoped_summary

    assert scoped_summary("8 chains · 757 residues", 6, 10) == "6 of 10 chains · 757 residues"
    assert scoped_summary("1 chain · 76 residues", 1, 24) == "1 of 24 chains · 76 residues"
    assert scoped_summary("40 nucleotides", 1, 2) == "1 of 2 chains · 40 nucleotides"


def test_the_summary_counts_the_chains_the_bar_counts():
    """Protein and nucleic-acid chains, in the bar and in the top bar alike: a chain picked on 1LMB (protein 3 with
    both DNA strands) is 3 of its 4 chains; 7A4M's one subunit is 1 of 24 in both."""
    at = opened(DATA / "1LMB.cif")
    assert at.session_state["info"].startswith("4 chains · ")
    at.multiselect(key="chains__pick").set_value(["3"]).run()
    assert not at.exception and at.session_state["info"].startswith("3 of 4 chains · ")
    at = opened(APO)
    assert at.session_state["info"].startswith("1 of 24 chains · ")
    assert "— 1 of 24 chains." in visible(rendered_text(at))
