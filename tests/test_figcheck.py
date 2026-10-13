"""figcheck: page checks for any figure (plan 10)."""

import math
from pathlib import Path

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath

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


def test_measles_f_copies_measured_unalike_today():
    """Measured at bcdec7a: a~c 33.0 page units. Task 2 turns this into the gate (<= 0.5). The ties are the ones
    asked for, kept although _separate dropped them, and chain b is matched through its copy classes."""
    d = drawn(DATA / "8UTF.cif")
    assert {c: ref for c, (ref, _) in d.layout.copy_ties.items()} == {"B": "A", "C": "A", "b": "a", "c": "a"}
    got = figcheck.copies_alike(d)
    assert set(got) == {("B", "A"), ("C", "A"), ("b", "a"), ("c", "a")} and all(math.isfinite(v) for v in got.values())
    assert max(got.values()) > 0.5


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


def test_disulfide_bar_over_a_label_is_a_text_collision():
    d = drawn(DATA / "r100" / "7FAB-HL.cif.gz")  # Fab, chains H and L (Task 1 fixture)
    hits = [h for h in figcheck.text_collisions(d) if "disulfide:" in h]
    assert hits  # Task 6 turns this into == []


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
        ax.text(*xy, "WWW", ha="center", va="center", fontsize=20, zorder=z).set_gid(gid)

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


def test_overlap_of_two_outlines_is_measured_on_their_area():
    square = MPath([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)], closed=True)
    corner = MPath([(1, 1), (5, 1), (5, 5), (1, 5), (1, 1)], closed=True)
    l_shape = MPath([(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4), (0, 0)], closed=True)
    assert figcheck.outline_overlap(square, corner) == pytest.approx(0.25, abs=0.03)  # 1 of the square's 4
    inside = MPath([(2, 2), (3.5, 2), (3.5, 3.5), (2, 3.5), (2, 2)], closed=True)  # in the L's box, not on the L
    assert figcheck.outline_overlap(l_shape, inside) == 0.0


def test_domain_panel_overlap_is_measured_on_1ldm():
    got = figcheck.domain_checks(drawn(FID / "1LDM.cif.gz", domains="auto"))
    assert max(f for *_, f in got["panel_overlaps"]) > 0.5  # D1·A1 / D1·A4, outline areas at bcdec7a; Task 3 gates it
    assert got["cross_chain"] == got["small"] == []


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
