"""figcheck: page checks for any figure (plan 10)."""

import math
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


def drawn(path, **kw):
    return draw_spec(figure_spec(path, **kw))


def test_crossing_counters_moved_but_unchanged():
    assert loop_crossings is figcheck.loop_crossings and honest_crossings is figcheck.honest_crossings


def test_ubiquitin_is_clean():
    d = drawn(DATA / "1UBQ.cif")
    got = figcheck.figure_checks(d, DATA / "1UBQ.cif", assembly="auto", focus="auto")
    assert got["element_overlaps"] == got["text_collisions"] == got["off_page"] == 0, got
    assert got["label_repeats"] == 0 and got["order"] == 0 and got["scope_silent"] == 0


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


def test_ink_checks_name_the_font_they_need():
    d = _page(_text("label:A:1-5"))
    d.figure.axes[0].texts[0].set_fontfamily(["DejaVu Sans"])
    assert "Liberation Sans" in figcheck.ink_font_reason(d.figure)
    assert figcheck.ink_font_reason(drawn(DATA / "1UBQ.cif").figure) in ("", figcheck.ink_font_reason(d.figure))


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


def test_copies_need_three_matches_and_half_the_smaller_copy():
    d = drawn(DATA / "1TIM.cif")
    keep = [k for k, p in d.layout.placed.items() if p.sse.chain == "B"]
    for k in keep[4:]:  # B keeps its elements, but only 4 of them still match A's
        d.layout.placed[k + "~"] = d.layout.placed.pop(k)
    assert figcheck.copy_matches(d)[("B", "A")] == (4, 20)
    assert figcheck.copies_alike(d)[("B", "A")] == math.inf  # 4 < 10, half of the smaller copy
    d = drawn(DATA / "1TIM.cif")
    for k in [k for k, p in d.layout.placed.items() if p.sse.chain == "B"][2:]:
        del d.layout.placed[k]
    assert figcheck.copy_matches(d)[("B", "A")] == (2, 2)
    assert figcheck.copies_alike(d)[("B", "A")] == math.inf  # two points always fit


def test_an_explicit_focus_still_checks_the_nucleotides_of_chains_it_draws():
    d = drawn(DATA / "1LMB.cif", assembly="asu", focus="3")
    assert figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "3") == []
    d.layout.stems.clear()  # the duplex gone, its chains' 5' tails still drawn
    got = figcheck.scope_honest(d, DATA / "1LMB.cif", "asu", "3")
    assert [s.split(":")[0] for s in got] == ["nt 1", "nt 2"], got


def test_na_crop_is_not_yet_counted_as_silent():
    """Pins a known gap: 4YBB's automatic na_extent crop leaves nucleotides out and no note says so yet, and
    scope_honest does not count it. Task 8 announces the crop (scope.drawn_scope) and turns this test into
    test_scope.py::test_na_crop_is_announced."""

    path = DATA / "na" / "4YBB.cif.gz"
    d = drawn(path)
    assert d.layout.na_cropped
    shown = sum(len(s.axial) for s in d.layout.stems) + sum(len(p.nts) for p in d.layout.na_paths)
    held = sum(len(c.get_polymer()) for c in d.model[0] if c.get_polymer().check_polymer_type().name in ("Rna", "Dna"))
    assert held - shown > 100  # the crop leaves out hundreds of nucleotides ...
    assert not [s for s in figcheck.scope_honest(d, path) if s.startswith("nt ")]  # ... and nothing says so
