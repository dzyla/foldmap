"""figcheck: page checks for any figure (plan 10)."""

import dataclasses
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath
from matplotlib.transforms import Bbox

from foldmap import figcheck
from foldmap.spec import draw_spec, figure_spec
from helpers import FID, honest_crossings, loop_crossings

DATA = Path(__file__).parent / "data"
R100 = DATA / "r100"


def drawn(path, **kw):
    return draw_spec(figure_spec(path, **kw))


def test_crossing_counters_moved_but_unchanged():
    assert loop_crossings is figcheck.loop_crossings and honest_crossings is figcheck.honest_crossings


@pytest.mark.parametrize("font", ["this machine's", "another"])
def test_ubiquitin_is_clean(font, monkeypatch):
    """Clean on every check. Where text is not drawn in figcheck.INK_FONT (macOS, Windows: Arial or DejaVu; forced
    here by asking for another font), the ink checks are not run: None, with the reason under "skipped"."""
    if font == "another":
        monkeypatch.setattr(figcheck, "INK_FONT", "No Such Sans")
    d = drawn(DATA / "1UBQ.cif")
    got = figcheck.figure_checks(d, DATA / "1UBQ.cif", assembly="auto", focus="auto")
    reason = figcheck.ink_font_reason(d.figure)
    assert got["element_overlaps"] == 0 and got["label_repeats"] == 0 and got["order"] == 0, got
    assert got["scope_silent"] == 0
    if reason:
        assert got["skipped"] == dict.fromkeys(figcheck.INK_CHECKS, reason)
        assert got["text_collisions"] is None and got["off_page"] is None
    else:
        assert got["skipped"] == {} and got["text_collisions"] == got["off_page"] == 0, got
    assert bool(reason) == (font == "another" or figcheck.ink_font_reason(d.figure) != "")


def test_tim_copies_are_alike_and_mirror_ties_count_as_alike():
    """1TIM and 1LDM tie their copies unmirrored (1LDM: A3 to A1, A2 to A4, in a cylinder frame); 1LMB's repressor
    dimer has its dyad in the page and draws chain 4 as chain 3's mirror image."""
    assert max(figcheck.copies_alike(drawn(DATA / "1TIM.cif")).values()) <= 1e-6
    ldm = drawn(FID / "1LDM.cif.gz")
    assert ldm.layout.copy_ties == {"A2": ("A4", False), "A3": ("A1", False)}
    assert max(figcheck.copies_alike(ldm).values()) <= 1e-6
    lmb = drawn(DATA / "1LMB.cif")
    assert lmb.layout.copy_ties == {"4": ("3", True)}
    assert max(figcheck.copies_alike(lmb).values()) <= 1e-6


def test_a_tied_copy_with_nothing_matched_is_a_violation():
    d = drawn(DATA / "1TIM.cif")
    copy, (ref, _) = sorted(d.layout.copy_ties.items())[0]
    for k in [k for k, p in d.layout.placed.items() if p.sse.chain == copy]:
        del d.layout.placed[k]
    assert figcheck.copies_alike(d)[(copy, ref)] == math.inf


def test_best_fit_deviation_turns_and_mirrors_only_when_asked():
    a = np.array([[0.0, 0.0], [4.0, 0.0], [4.0, 1.0], [0.0, 3.0]])
    turn = np.array([[0.0, -1.0], [1.0, 0.0]])
    assert figcheck.fit_deviation(a, a @ turn.T + [10.0, -2.0], mirrored=False) < 1e-9
    mirror = a * [-1.0, 1.0] + [7.0, 0.0]
    assert figcheck.fit_deviation(a, mirror, mirrored=True) < 1e-9
    assert figcheck.fit_deviation(a, mirror, mirrored=False) > 0.5


def _page(*artists_and_texts):
    """A 10 x 10 unit page at 100 dpi (one unit = 50 px), its artists added, as a Drawn-like object."""
    fig = Figure(figsize=(5, 5), dpi=100)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    for make in artists_and_texts:
        make(ax)

    class D:
        figure = fig

    return D


def _text(gid, xy=(5, 5), z=4):
    def make(ax):
        t = ax.text(*xy, "WWW", ha="center", va="center", fontsize=20, zorder=z)
        if gid:
            t.set_gid(gid)

    return make


def _line(gid, a, b, z=1.0):
    def make(ax):
        ax.add_patch(PathPatch(MPath([a, b]), fc="none", ec="k", lw=2, zorder=z)).set_gid(gid)

    return make


def _body(gid, x0, y0, x1, y1, z=3):
    def make(ax):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="b", ec="none", zorder=z)).set_gid(gid)

    return make


def test_text_collisions_on_a_synthetic_page():
    """A line through a text counts whatever its z-order, unless a body drawn above it covers it where it passes
    the text (a loop under its strand's letter); a body under a text counts unless the text is its own."""
    d = _page(_text("label:A:1-5"), _line("disulfide:1-9", (0, 5), (10, 5), 3.7))
    assert figcheck.text_collisions(d) == ["label:A:1-5 × disulfide:1-9"]
    d = _page(_text("label:A:1-5"), _body("strand:A:1-5", 2, 4, 8, 6), _line("loop:A:1-5>A:9-12", (0, 5), (10, 5)))
    assert figcheck.text_collisions(d) == []  # its own body, over the loop
    d = _page(_text("label:A:1-5"), _body("strand:B:1-5", 2, 4, 8, 6))
    assert figcheck.text_collisions(d) == ["label:A:1-5 × strand:B:1-5"]
    d = _page(_text("label:A:1-5"), _text("label:B:1-5", (5.3, 5.1)), _line("stub:N:A", (0, 5), (10, 5)))
    assert figcheck.text_collisions(d) == ["label:A:1-5 × label:B:1-5"]  # a stub is no line the check reads
    assert figcheck.off_page(_page(_text("label:A:1-5", (9.9, 5)))) == ["label:A:1-5"]


def test_a_near_touch_is_no_text_collision():
    """R13: a text or a body counts against a text only reaching more than TEXT_DEPTH_PX into its ink (1 px of
    anti-aliasing or of matplotlib versions never flips a gate); deeper, it does."""
    ink = dict((t.get_gid(), b) for t, b in figcheck._ink_boxes(_page(_text("label:A:1-5")).figure))["label:A:1-5"]
    px = 10 / 500  # page units per pixel on _page
    top = ink.y1 * px
    for depth, hit in ((1.5, False), (3.0, True)):
        d = _page(_text("label:A:1-5"), _body("strand:B:1-5", 2, top - depth * px, 8, 9))
        assert bool(figcheck.text_collisions(d)) is hit, depth
        d = _page(_text("label:A:1-5"), _text("label:B:1-5", (5, 5 + (ink.height - depth) * px)))
        assert bool(figcheck.text_collisions(d)) is hit, depth


def test_overlap_of_two_outlines_is_measured_on_their_area():
    square = MPath([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)], closed=True)
    corner = MPath([(1, 1), (5, 1), (5, 5), (1, 5), (1, 1)], closed=True)
    l_shape = MPath([(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4), (0, 0)], closed=True)
    assert figcheck.outline_overlap(square, corner) == pytest.approx(0.25, abs=0.03)  # 1 of the square's 4
    inside = MPath([(2, 2), (3.5, 2), (3.5, 3.5), (2, 3.5), (2, 2)], closed=True)  # in the L's box, not on the L
    assert figcheck.outline_overlap(l_shape, inside) == 0.0


def test_a_panel_in_another_panels_hole_does_not_overlap_it():
    """C21 (random-100): a panel with a hole (another domain's element interleaved in it: the contour's inner ring
    runs the other way) is filled around the hole, as Agg fills by the nonzero rule; the overlap measure counted
    the hole as inside (Path.contains_points tests each ring alone). A panel in the hole is clear of it; one across
    the ring still overlaps."""
    ring = MPath.make_compound_path(
        MPath([(0, 0), (6, 0), (6, 6), (0, 6), (0, 0)], closed=True),
        MPath([(2, 2), (2, 4), (4, 4), (4, 2), (2, 2)], closed=True),  # clockwise: the hole
    )
    island = MPath([(2.2, 2.2), (3.8, 2.2), (3.8, 3.8), (2.2, 3.8), (2.2, 2.2)], closed=True)
    assert figcheck.outline_overlap(ring, island) == 0.0
    across = MPath([(3, 3), (5, 3), (5, 3.8), (3, 3.8), (3, 3)], closed=True)  # half in the hole, half on the ring
    assert figcheck.outline_overlap(ring, across) == pytest.approx(0.5, abs=0.05)


def test_interleaved_copies_panels_are_measured_with_their_holes():
    """C21: 3L9G's domains view (four copies, shared sheets): the panel fields share no cell, and neither do the
    drawn outlines (D1 · A3 has a hole where D3 · A4's strand sits; measured 0.18 overlap before). Fix round 1: the
    test first asserts the case it measures (the copies interleaved on one page, panels with holes), so it cannot
    pass on an exploded page as it did on be5c237."""
    from foldmap.layout import domain_fields

    d = drawn(DATA / "r100" / "3L9G-A.cif.gz", domains="auto")
    plain = drawn(DATA / "r100" / "3L9G-A.cif.gz").layout
    assert d.layout.width < 1.5 * plain.width  # the copies interleaved, not parted apart (C31): else nothing to test
    to_page = d.figure.axes[0].transData.inverted()
    rings = [
        ring
        for a in d.figure.axes[0].patches
        if (a.get_gid() or "").startswith("domain:")
        for ring in a.get_path().transformed(a.get_transform() + to_page).to_polygons(closed_only=True)
    ]
    turns = [np.sum(r[:-1, 0] * r[1:, 1] - r[1:, 0] * r[:-1, 1]) for r in rings]
    # the case measured: some panel has a hole (an interleaved domain's element), or (plan 12 Task 2: 3L9G's barrels
    # are 16-strand rings, their branch strands stacked beside them) some domain's panel is cut round another's
    assert any(t < 0 for t in turns) or len(rings) > len(d.layout.domains)
    fields = domain_fields(d.layout)
    inside = {name: v > 0 for name, v in fields.values.items()}
    assert not any((inside[a] & inside[b]).any() for a, b in combinations(sorted(inside), 2))
    assert figcheck.domain_checks(d)["panel_overlaps"] == []


def test_domain_panels_drawn_on_one_another_are_found():
    """domain_checks reads the drawn panels (domain:<name>): a copy of a panel drawn as another domain's overlaps it
    whole. test_domains.py::test_domain_panels_do_not_overlap gates real figures."""
    d = drawn(
        DATA / "2LZM.cif", assembly="asu", domains=[("N", ["res:A:13-59"]), ("C", ["res:A:1-12", "res:A:60-164"])]
    )
    assert figcheck.domain_checks(d)["panel_overlaps"] == []
    ax = d.figure.axes[0]
    (panel,) = [p for p in ax.patches if p.get_gid() == "domain:N"]
    ax.add_patch(PathPatch(panel.get_path(), transform=panel.get_transform())).set_gid("domain:copy")
    assert figcheck.domain_checks(d)["panel_overlaps"] == [("N", "copy", 1.0)]


def test_stem_names_are_measured_to_their_own_letters():
    d = drawn(DATA / "na" / "1EHZ.cif.gz")
    got = figcheck.stem_name_gaps(d)
    assert got and set(got) == {f"na:label:{k}" for k, s in enumerate(d.layout.stems) if not s.ghost}
    assert all(0 <= own < other for own, other in got.values())


def test_scope_silence_is_found():
    """A figure that drops a chain without a note is caught: fake it by removing a chain's elements after layout."""
    d = drawn(DATA / "1LMB.cif", assembly="asu")
    victim = sorted({p.sse.chain for p in d.layout.placed.values()})[-1]
    assert f"chain {victim}" not in figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "auto")
    for k in [k for k, p in d.layout.placed.items() if p.sse.chain == victim]:
        del d.layout.placed[k]
    assert f"chain {victim}" in figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "auto")


def test_a_stretch_of_nucleotides_left_out_is_found():
    """1EHZ draws every nucleotide; without its D-loop hairpin (A:14-21) the figure leaves 8 out unsaid."""
    path = DATA / "na" / "1EHZ.cif.gz"
    d = drawn(path)
    assert figcheck.scope_honest(d, path) == []
    k = next(n for n, p in enumerate(d.layout.na_paths) if p.kind == "hairpin" and len(p.nts) == 8)
    del d.layout.na_paths[k]
    assert figcheck.scope_honest(d, path) == ["nt A:14-21"]


def test_label_repeats_name_chains_and_stems():
    from types import SimpleNamespace as NS

    el = lambda chain, label: NS(sse=NS(chain=chain), label=label)  # noqa: E731
    stem = lambda name, ghost=False: NS(name=name, ghost=ghost)  # noqa: E731
    lay = NS(
        placed={1: el("A", "α1"), 2: el("A", "α1"), 3: el("B", "α1")},
        stems=[stem("P1"), stem("P1"), stem("P2"), stem("P2", ghost=True)],
    )
    assert figcheck.label_repeats(lay) == {"A": ["α1"], "<na>": ["P1"]}


def test_chains_asked_for_are_a_choice():
    d = drawn(DATA / "1LMB.cif", assembly="asu", chains=["3"])
    assert figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "auto", chains=["3"]) == []


@pytest.mark.parametrize("bad", ["swap_strands", "skip_element"])
def test_order_checks_catch_a_swapped_or_skipped_element(bad):
    d = drawn(DATA / "1UBQ.cif")
    assert figcheck.order_checks(d) == []
    if bad == "swap_strands":  # two strands of the sheet trade places on the page
        ids = next(ids for ids in d.layout.sheets if len(ids) >= 3)
        a, b = d.layout.placed[ids[0]], d.layout.placed[ids[-1]]
        (a.cx, a.cy), (b.cx, b.cy) = (b.cx, b.cy), (a.cx, a.cy)
    else:  # one loop runs from an element past the next one, whose two loops are gone
        k = next(k for k in range(len(d.loops) - 1) if d.loops[k].b_id == d.loops[k + 1].a_id)
        first, second = d.loops[k], d.loops[k + 1]
        first.b_id = second.b_id
        del d.loops[k + 1]
    got = figcheck.order_checks(d)
    assert len(got) == 1, got


# ------------------------------------------------------------------------------------- fix round 1 (review)


def _liberation(d):
    """Ink measures are calibrated on Liberation Sans (the font CI installs): elsewhere a gated check is skipped."""
    reason = figcheck.ink_font_reason(d.figure)
    if reason:
        pytest.skip(reason)


def test_ink_is_measured_under_a_fixed_hinting_and_the_font_is_recorded():
    """The ink boxes do not follow the caller's hinting settings (matplotlib 3.10 and 3.11 default differently)."""
    from matplotlib import rc_context

    d = drawn(DATA / "na" / "1EHZ.cif.gz")
    plain = figcheck.stem_name_gaps(d)
    with rc_context({"text.hinting": "force_autohint"}):
        assert figcheck.stem_name_gaps(d) == plain
    got = figcheck.figure_checks(d, DATA / "na" / "1EHZ.cif.gz")
    assert got["font"] and all(Path(f).suffix.lower() in (".ttf", ".otf") for f in got["font"])
    assert got["skipped"] == {} or set(got["skipped"]) == set(figcheck.INK_CHECKS)


def _names_its_fonts(reason: str, fig) -> bool:
    """The reason says which font the ink checks need and names every other font file the texts resolve to."""
    from matplotlib.ft2font import FT2Font

    other = [f for f in figcheck.text_fonts(fig) if FT2Font(f).family_name != figcheck.INK_FONT]
    return reason.startswith(f"ink checks are calibrated on {figcheck.INK_FONT}") and all(f in reason for f in other)


def test_ink_checks_name_the_font_they_need():
    """A text in DejaVu Sans (matplotlib's own font, on every machine) is not the calibrated font; the figure's own
    texts resolve to Liberation Sans on Linux CI, and to whatever this machine has elsewhere (Arial, DejaVu)."""
    d = _page(_text("label:A:1-5"))
    d.figure.axes[0].texts[0].set_fontfamily(["DejaVu Sans"])
    reason = figcheck.ink_font_reason(d.figure)
    assert _names_its_fonts(reason, d.figure) and "DejaVuSans" in reason
    ubq = drawn(DATA / "1UBQ.cif").figure
    reason = figcheck.ink_font_reason(ubq)
    assert reason == "" or _names_its_fonts(reason, ubq)


def test_gate_tolerances_for_later_tasks():
    assert figcheck.stem_gap_tolerance(0.02) == figcheck.STEM_GAP_FLOOR == 0.08
    assert figcheck.stem_gap_tolerance(0.2) == pytest.approx(0.201)
    assert figcheck.LINE_TEXT_CLEARANCE_PX >= 2.0
    near = _page(_text("label:A:1-5"), _line("disulfide:1-9", (0, 4.80), (10, 4.80), 3.7))  # 2.8 px under the ink
    assert figcheck.text_collisions(near) == []
    assert figcheck.text_collisions(near, clearance_px=figcheck.LINE_TEXT_CLEARANCE_PX) == [
        "label:A:1-5 × disulfide:1-9"
    ]


def test_every_text_must_reach_the_glyph_recorder():
    d = _page(_text("label:A:1-5"))
    d.figure.axes[0].texts[0].draw = lambda renderer: None  # a text drawing nothing would vanish from every check
    with pytest.raises(RuntimeError, match="label:A:1-5"):
        figcheck.text_collisions(d)


def test_a_legend_mark_off_the_page_is_found():
    def legend(ax):
        ax.add_patch(Rectangle((9.8, 1), 0.6, 0.6, fc="r")).set_gid("legend-chain:A")

    assert figcheck.off_page(_page(legend)) == ["legend-chain:A"]


def test_element_overlaps_on_synthetic_boxes():
    from types import SimpleNamespace as NS

    el = lambda k, rect: NS(sse=NS(id=k), rect=rect)  # noqa: E731
    lay = NS(
        placed={"a": el("a", (0, 0, 2, 2)), "b": el("b", (1, 1, 3, 3)), "c": el("c", (2.5, 0, 4, 0.5))},
        stems=[],
        nucleic=None,
        on_sheet={},
    )
    assert figcheck.element_overlaps(lay) == [("a", "b")]
    lay.on_sheet = {"a": ["b"]}  # a helix drawn over its sheet's strand: by design
    assert figcheck.element_overlaps(lay) == []


def _stem_page(*more):
    """A page holding one stem (P1, 4 x 2 units at the centre) and the given artists."""
    from foldmap.layout import PlacedStem

    stem = PlacedStem("s1", "P1", [0, 0], [], 5.0, 5.0, 4.0, 0.0, 2.0, {}, {})
    d = _page(*more)
    d.layout = type("L", (), {"stems": [stem], "nucleic": None, "on_sheet": {}, "na_paths": []})()
    return d


def test_text_on_a_stem_body_is_a_collision_unless_it_is_the_stems_own():
    assert figcheck.text_collisions(_stem_page(_text("label:A:1-5"))) == ["label:A:1-5 × na:stem:0"]
    assert figcheck.text_collisions(_stem_page(_text("na:label:0"))) == []


def test_a_helix_body_counts_once_and_texts_without_gid_are_told_apart():
    d = _page(_text("label:A:1-5"), _body("helix-back:B:1-9", 2, 4, 8, 6, 2.9), _body("helix:B:1-9", 2, 4, 8, 6))
    assert figcheck.text_collisions(d) == ["label:A:1-5 × helix:B:1-9"]
    d = _page(_text(None), _text(None, (5.3, 5.1)))
    assert len(figcheck.text_collisions(d)) == 1  # the same words, two texts


def test_a_loop_length_is_not_crossed_by_its_own_loop():
    own = _page(_text("loop-length:A:1-5>A:9-12"), _line("loop:A:1-5>A:9-12", (0, 5), (10, 5)))
    assert figcheck.text_collisions(own) == []
    seg = _page(_text("loop-length:A:1-5>A:9-12"), _line("loop-seg:A:1-5>A:9-12:0", (0, 5), (10, 5)))
    assert figcheck.text_collisions(seg) == []
    other = _page(_text("loop-length:A:1-5>A:9-12"), _line("loop:B:1-5>B:9-12", (0, 5), (10, 5)))
    assert figcheck.text_collisions(other) == ["loop-length:A:1-5>A:9-12 × loop:B:1-5>B:9-12"]


def test_texts_turned_alike_are_compared_in_their_turned_frame():
    """Two texts turned 30 degrees alike, one beside the other across their baseline: their page boxes overlap (a
    turned text's corners), their ink in the turned frame does not, so they are no collision (Page.pair_boxes)."""

    def turned(gid, xy):
        def make(ax):
            ax.text(*xy, "WWWWWW", ha="center", va="center", fontsize=20, rotation=30).set_gid(gid)

        return make

    across = np.array([-np.sin(np.radians(30)), np.cos(np.radians(30))]) * 0.75  # 37 px: a line and a bit
    d = _page(turned("a", (5, 5)), turned("b", tuple(np.add((5, 5), across))))
    page = figcheck.measure(d)
    (ta, _, a), (tb, _, b) = page.texts
    assert a.overlaps(b) and figcheck._texts_meet(a, b)  # by the page boxes, a collision
    assert not Bbox.overlaps(*page.pair_boxes(ta, a, tb, b))
    assert figcheck.text_collisions(d) == []


def test_a_stem_name_is_measured_across_its_row():
    """1EHZ named: the D arm has four pairs, so its name centred along the arm sits between two letters of its row;
    its gap is measured across the stem from the row (the margin), not to the nearer letter's corner (plan 10 Task 4
    fix round 1)."""
    from foldmap.layout import _INK_MARGIN

    d = drawn(DATA / "na" / "1EHZ.cif.gz", names=[str(Path(__file__).parents[1] / "docs" / "names" / "1EHZ-trna.yaml")])
    _liberation(d)
    k = next(k for k, s in enumerate(d.layout.stems) if s.name == "D")
    ax = d.figure.axes[0]
    em = next(t for t in ax.texts if t.get_gid().startswith("na:nt:")).get_fontsize() / 0.75 * d.figure.dpi / 72
    em /= ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    own, _ = figcheck.stem_name_gaps(d)[f"na:label:{k}"]
    assert own == pytest.approx(_INK_MARGIN * em, abs=2e-3)
    to_page = ax.transData.inverted()
    ink = {n: to_page.transform_bbox(b).extents for _, n, b in figcheck.measure(d).texts}
    na, s = d.layout.nucleic, d.layout.stems[k]
    tags = [f"na:nt:{na.strands[a].labels[i].chain}:{na.strands[a].labels[i].tag}" for a, i in s.axial]
    assert min(figcheck._box_gap(ink[f"na:label:{k}"], ink[t]) for t in tags) > own + 0.05  # its nearest letter


def test_a_line_is_measured_to_a_text_box_along_its_segments():
    """figcheck.stem_name_line_gaps measures a line segment by segment, not at its vertices: a long straight run
    passing a name has no vertex near it (plan 10 Task 4 fix round 2)."""
    box = (0.0, 0.0, 1.0, 1.0)
    assert figcheck._segment_box_gap((-5, 1.5), (5, 1.5), box) == pytest.approx(0.5)  # passes over it, ends far
    assert figcheck._segment_box_gap((-5, 0.5), (5, 0.5), box) == 0.0  # through it
    assert figcheck._segment_box_gap((2, 3), (3, 2), box) == pytest.approx(3 / math.sqrt(2))  # a corner to its middle
    assert figcheck._segment_box_gap((3, 0.5), (3, 0.5), box) == pytest.approx(2.0)  # a point


def test_stem_names_are_measured_to_letters_not_to_the_stem_box():
    """1EHZ named: the D arm's name sits inside its turned stem's box (page-upright) but clear of its letters (plan 10
    Task 4 put it 0.5 em beyond them; before, it sat near the box and far from the letters)."""
    d = drawn(DATA / "na" / "1EHZ.cif.gz", names=[str(Path(__file__).parents[1] / "docs" / "names" / "1EHZ-trna.yaml")])
    _liberation(d)
    got = figcheck.stem_name_gaps(d)
    k = next(k for k, s in enumerate(d.layout.stems) if s.name == "D")
    to_page = d.figure.axes[0].transData.inverted()
    name = next(t for t in d.figure.axes[0].texts if t.get_gid() == f"na:label:{k}")
    box = to_page.transform_bbox(name.get_window_extent(figcheck._renderer(d.figure))).extents
    rect_gap = figcheck._box_gap(box, d.layout.stems[k].rect)
    assert rect_gap == 0 and got[f"na:label:{k}"][0] > 0.05  # measured to the box, the gap would be 0


def test_a_strand_drawn_the_wrong_way_is_an_order_violation():
    d = drawn(DATA / "1UBQ.cif")
    ids = next(ids for ids in d.layout.sheets if len(ids) >= 3)
    d.layout.placed[ids[0]].angle += math.pi  # the sheet's edge strand turned round
    got = figcheck.order_checks(d)
    assert len(got) == 1 and ids[0] in got[0] and "drawn" in got[0], got


def test_a_domain_across_two_chains_is_found():
    d = drawn(DATA / "1TIM.cif", domains=[("both", ["res:A:1-60", "res:B:1-60"]), ("rest", ["res:A:61-250"])])
    assert figcheck.domain_checks(d)["cross_chain"] == ["both"]


def test_copies_split_differently_are_found():
    alike = [("N", ["res:A:1-100", "res:B:1-100"]), ("C", ["res:A:101-250", "res:B:101-250"])]
    unalike = [("N", ["res:A:1-100", "res:B:1-50"]), ("C", ["res:A:101-250", "res:B:51-250"])]
    assert figcheck.domain_checks(drawn(DATA / "1TIM.cif", domains=alike))["unalike"] == []
    assert figcheck.domain_checks(drawn(DATA / "1TIM.cif", domains=unalike))["unalike"]


def test_a_small_automatic_domain_is_found():
    """A domain named as automatic ones are (D1, D2 ...) holding fewer than DOMAIN_MIN_SSE_RESIDUES residues of
    elements. test_domains.py::test_domains_auto_rules_on_oligomers gates --domains auto (8UUP's D5 · F at 396e4a1)."""
    d = drawn(DATA / "1TIM.cif", domains=[("D1", ["res:A:1-20"]), ("D2", ["res:A:21-250"])])
    assert figcheck.domain_checks(d)["small"] == ["D1"]


def test_a_domain_name_on_a_body_or_under_a_line_is_found():
    """name_collisions: a domain name moved onto a helix, or crossed by a loop, is a collision; on its own plate it
    is not. test_domains.py::test_domain_panels_do_not_overlap gates real figures (Cas9 had 4 at 396e4a1)."""
    d = drawn(
        DATA / "2LZM.cif", assembly="asu", domains=[("N", ["res:A:13-59"]), ("C", ["res:A:1-12", "res:A:60-164"])]
    )
    _liberation(d)
    assert figcheck.domain_checks(d)["name_collisions"] == []
    ax = d.figure.axes[0]
    (name,) = [t for t in ax.texts if t.get_gid() == "domain-label:N"]
    helix = next(p for p in d.layout.placed.values() if p.sse.kind == "H")
    name.set_position((helix.cx, helix.cy))
    got = figcheck.domain_checks(d)["name_collisions"]
    assert f"domain-label:N × helix:{helix.sse.id}" in got, got
    loop = max(d.loops, key=lambda lp: len(lp.points))  # on the middle of a loop's longest run, off any body
    a, b = max(zip(loop.points, loop.points[1:]), key=lambda ab: math.dist(*ab))
    name.set_position(((a[0] + b[0]) / 2, (a[1] + b[1]) / 2))
    got = figcheck.domain_checks(d)["name_collisions"]
    assert f"domain-label:N × loop:{loop.a_id}>{loop.b_id}" in got, got


def _unmatched(p) -> None:
    """Make a placed element match no copy's: the other kind at the same residues."""
    p.sse = dataclasses.replace(p.sse, kind="E" if p.sse.kind == "H" else "H")


def test_copies_need_three_matches_and_half_the_smaller_copy():
    d = drawn(DATA / "1TIM.cif")
    keep = [k for k, p in d.layout.placed.items() if p.sse.chain == "B"]
    for k in keep[4:]:  # B keeps its elements, but only 4 of them still match A's
        _unmatched(d.layout.placed[k])
    assert figcheck.copy_matches(d)[("B", "A")] == (4, 20)
    assert figcheck.copies_alike(d)[("B", "A")] == math.inf  # 4 < 10, half of the smaller copy
    d = drawn(DATA / "1TIM.cif")
    for k in [k for k, p in d.layout.placed.items() if p.sse.chain == "B"][2:]:
        del d.layout.placed[k]
    assert figcheck.copy_matches(d)[("B", "A")] == (2, 2)
    assert ("B", "A") not in figcheck.copies_alike(d)  # a copy of two elements is reported, not judged (C22)


def test_copies_of_fewer_than_three_elements_are_reported_not_judged():
    """C22 (random-100): a copy with fewer than COPY_MIN_MATCHED elements (7AOD's two-helix chain D~P, 6JLN's
    single-TM subunits) can never match three; it is reported as not judged (ruling), not failed as inf. So is a
    tie left with fewer than three after its shared sheet is left out (C10; 6YQD). A copy that draws nothing at all
    stays a violation."""
    d = drawn(DATA / "1TIM.cif")
    for chain in "AB":
        for k in [k for k, p in d.layout.placed.items() if p.sse.chain == chain][2:]:
            del d.layout.placed[k]
    assert figcheck.copies_alike(d) == {} and figcheck.copy_matches(d) == {("B", "A"): (2, 2)}
    assert figcheck.copies_not_judged(d) == {("B", "A"): "n/a: 2 elements"}
    d = drawn(R100 / "7IMM-A.cif.gz")
    other = [k for k, p in d.layout.placed.items() if p.sse.chain == "A2" and p.sse.kind != "E"]
    for k in other[3:]:
        del d.layout.placed[k]  # A2 keeps its shared strands and three other elements, two of them matching A1's
    _unmatched(d.layout.placed[other[2]])
    assert figcheck.copy_matches(d)[("A2", "A1")] == (2, 3)
    assert figcheck.copies_not_judged(d) == {("A2", "A1"): "n/a: 2 matched"}  # two fit nothing (6YQD, round 2)
    assert ("A2", "A1") not in figcheck.copies_alike(d)
    _unmatched(d.layout.placed[other[1]])  # one matched or none: a real inf, not "not judged"
    assert figcheck.copies_not_judged(d) == {}
    assert figcheck.copies_alike(d)[("A2", "A1")] == math.inf
    shared = figcheck.shared_sheet_elements(d.layout, "A2", "A1")
    for k in [k for k, p in d.layout.placed.items() if p.sse.chain in ("A1", "A2") and k not in shared]:
        del d.layout.placed[k]  # the shared sheet took every element: nothing left to judge
    assert figcheck.copies_not_judged(d) == {("A2", "A1"): "n/a: shared sheet"}


def test_strands_of_a_sheet_shared_between_copies_are_left_out():
    """C10 (random-100): 7IMM's two copies make one continuous A B B A sheet: one item holding both, placed by the
    sheet, so its strands cannot sit at their copy's tied place (5.63 page units from it); A15 leaves them out
    (ruling) and judges the rest. A sheet with a third copy's strand, tied round the ring, stays in (8UTF)."""
    d = drawn(R100 / "7IMM-A.cif.gz")
    assert d.layout.copy_ties == {"A2": ("A1", False)}
    assert figcheck.shared_sheet_elements(d.layout, "A2", "A1") == {
        "A1:319-325",
        "A1:328-338",
        "A2:328-338",
        "A2:319-325",
    }
    assert figcheck.copy_matches(d)[("A2", "A1")] == (8, 8)
    assert figcheck.copies_alike(d)[("A2", "A1")] <= 0.5
    moved = next(k for k, p in d.layout.placed.items() if p.sse.chain == "A2" and k not in {"A2:328-338", "A2:319-325"})
    d.layout.placed[moved].cx += 5.0  # any element outside the shared sheet still counts
    assert figcheck.copies_alike(d)[("A2", "A1")] > 1.0
    utf = drawn(DATA / "8UTF.cif")  # tied ring sheets stay measured (fix round 1): nothing left out, all matched
    assert figcheck.shared_sheet_elements(utf.layout, "b", "a") == set()
    assert figcheck.copy_matches(utf)[("b", "a")] == (25, 25) and figcheck.copies_alike(utf)[("b", "a")] < 1e-6
    uup = drawn(DATA / "8UUP.cif")
    assert figcheck.shared_sheet_elements(uup.layout, "D", "B") == set()
    assert figcheck.copy_matches(uup)[("D", "B")] == (33, 33)
    utf.layout.placed[
        "b:"
        + next(k for k in utf.layout.placed if k.startswith("b:") and utf.layout.placed[k].sse.kind == "E").split(
            ":", 1
        )[1]
    ].cx += 3.0
    assert figcheck.copies_alike(utf)[("b", "a")] > 1.0  # drift in a tied ring sheet is seen


def test_copies_are_matched_by_sequence_not_numbers(tmp_path):
    """C18 (random-100): copies pair their elements by sequence, not by author numbers (ruling). 1TIM with chain B
    renumbered from residue 101 on (+1000): by numbers only 8 of its 20 elements paired A's and the tie read inf;
    aligned, all 20 pair and the drawing is judged."""
    import gemmi

    st = gemmi.read_structure(str(DATA / "1TIM.cif"))
    st.setup_entities()
    for r in st[0]["B"]:
        if r.seqid.num > 100:
            r.seqid.num += 1000
    path = tmp_path / "1TIM-renumbered.cif"
    st.make_mmcif_document().write_file(str(path))
    d = drawn(path)
    assert d.layout.copy_ties == {"B": ("A", False)}
    assert figcheck.copy_matches(d) == {("B", "A"): (20, 20)}
    assert math.isfinite(figcheck.copies_alike(d)[("B", "A")])


def test_untied_copies_of_a_symmetric_assembly_are_checked_too():
    """C18 (random-100): 2MJW ties A~C but not B~D (DSSP split D's elements), so B and D were drawn unlike with no
    check. A chain no tie holds, in a figure that ties copies, is checked against its copy found by sequence
    (25IN with D's tie dropped: D is matched to C by sequence and judged)."""
    d = drawn(DATA / "25IN.cif")
    del d.layout.copy_ties["D"]
    assert set(figcheck.copies_alike(d)) == {("B", "A"), ("D", "C")}
    assert figcheck.copies_alike(d)[("D", "C")] <= 1e-6
    moved = next(p for p in d.layout.placed.values() if p.sse.chain == "D")
    moved.cx += 3.0
    assert figcheck.copies_alike(d)[("D", "C")] > 0.5
    assert figcheck.aligned_positions(["ALA", "GLY", "SER", "TRP", "LYS"], ["GLY", "SER", "LYS"]) == [
        (1, 0),
        (2, 1),
        (4, 2),
    ]


def test_an_explicit_focus_still_checks_the_nucleotides_of_chains_it_draws():
    d = drawn(DATA / "1LMB.cif", assembly="asu", focus="3")
    assert figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "3") == []
    d.layout.stems.clear()  # the duplex gone, its chains' 5' tails still drawn
    got = figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "3")
    assert [s.split(":")[0] for s in got] == ["nt 1", "nt 2"], got


def test_same_protein_needs_a_length_ratio():
    """A 30-residue chain aligned inside a 300-residue one is no copy of it (C30 fix round 1)."""
    long = ["ALA", "GLY", "LEU", "SER", "VAL", "LYS", "GLU", "ASP", "THR", "PRO"] * 30
    short = long[100:130]
    assert not figcheck.same_protein(long, short, figcheck.aligned_positions(long, short))
    assert figcheck.same_protein(long, long[:260], figcheck.aligned_positions(long, long[:260]))


@pytest.mark.parametrize("name, named, of", [("2OMF-A", 3, 3), ("3L9G-A", 8, 16)])
def test_accessories_named_not_drawn_level_are_reported(name, named, of):
    """Plan 12 ruling (fix round 4): S6 passes an accessory drawn level beside its partner or named in the notes, so
    the count named rather than drawn level is the measure of placement. figure_checks reports it (not a gate):
    OmpF's three A:2-6 are named, 8 of 3L9G's 16 accessories."""
    path = R100 / f"{name}.cif.gz"
    d = drawn(path)
    lines = figcheck.accessories_named(d.layout)
    assert len(lines) == named and all(" also pairs with " in line for line in lines)
    got = figcheck.figure_checks(d, path)
    assert (got["accessories_named"], got["accessories"]) == (named, of)
    assert got["details"]["accessories_named"] == lines
