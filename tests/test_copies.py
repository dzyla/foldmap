"""Copies of one protein tied by symmetry are drawn alike (figcheck.copies_alike), the measles F trimer too."""

import math
from pathlib import Path

import numpy as np
import pytest

from foldmap import figcheck
from foldmap.spec import draw_spec, figure_spec
from helpers import FID

DATA = Path(__file__).parent / "data"
ALIKE = 0.5  # page units: the largest element-centre deviation between two tied copies (A15, Q-R12)


def drawn(path, **kw):
    return draw_spec(figure_spec(path, **kw))


@pytest.mark.parametrize(
    "name, ties",
    [
        ("8UTF", {"B": "A", "C": "A", "b": "a", "c": "a"}),
        ("8UUP", {"C": "A", "E": "A", "D": "B", "F": "B"}),
        ("25IN", {"B": "A", "D": "C"}),
    ],
)
def test_measles_f_trimer_copies_drawn_alike(name, ties):
    """Every element two tied copies share sits alike, matched through the copy classes. Before: 8UTF c~a 33.0 and
    b~a 14.6 (the sheets that hold a neighbour's strand were never tied, and the 6-helix bundle was one item);
    8UUP D~B 25.1; 25IN B~A 20.4 (its swapped helix pairs were one bundle each)."""
    d = drawn(DATA / f"{name}.cif")
    assert {c: ref for c, (ref, _) in d.layout.copy_ties.items()} == ties
    got = figcheck.copies_alike(d)
    assert set(got) == set(ties.items())
    assert max(got.values()) <= ALIKE, {f"{a}~{b}": round(v, 3) for (a, b), v in got.items()}


@pytest.mark.parametrize("path", [DATA / "1TIM.cif", DATA / "6ZYA.cif", FID / "1LDM.cif.gz"], ids=lambda p: p.name)
def test_symmetric_controls_drawn_alike(path):
    got = figcheck.copies_alike(drawn(path))
    assert all(math.isfinite(v) and v <= 1e-6 for v in got.values()), got


def test_a_strand_lent_across_the_cut_leaves_a_ghost_where_it_pairs():
    """8UTF: each protomer lends its S and T strands (440-442, 445-446) to the next protomer's sheets. Unrolled into a
    row a, b, c, chain c's strands would pair in a's sheets at the far end; they are drawn past c instead, where its
    copies' strands sit, and a's sheets keep a ghost of each. The order check reads a sheet with its ghosts."""
    from foldmap.layout import block_ghosts

    d = drawn(DATA / "8UTF.cif")
    lay = d.layout
    hosts = {g.sse.id: lay.sheets[k] for k in range(len(lay.sheets)) for g in block_ghosts(lay, k)}
    hosts = {i: block for i, block in hosts.items() if i.startswith("c:")}  # their partners' ghosts aside
    assert sorted(hosts) == ["c:440-442", "c:445-446"]
    assert "a:325-329" in hosts["c:440-442"] and hosts["c:445-446"] == ["A:61-62"]
    right = max(p.rect[2] for k, p in lay.placed.items() if p.sse.chain in "Cc" and k not in hosts)
    assert all(lay.placed[k].cx > right for k in hosts)  # past the rest of copy c
    assert figcheck.order_checks(d) == []


def test_a_sheet_panel_covers_its_own_ghosts_only():
    """8UUP: A's and B's strands lent to F's sheet (protomer 2 holds most of it) across the cut stand before the
    row's start as a block of their own; its panel covers them and their partners' ghosts, not their own ghosts in
    F's sheet at the other end of the page."""
    from foldmap.layout import block_ghosts

    d = drawn(DATA / "8UUP.cif")
    lay = d.layout
    ghosted = {g.sse.id for g in lay.ghosts}
    lent = [k for k, ids in enumerate(lay.sheets) if set(ids) <= ghosted and len(ids) > 1]
    host = [k for k in range(len(lay.sheets)) if set(lay.sheets[lent[0]]) & {g.sse.id for g in block_ghosts(lay, k)}]
    assert len(lent) == 1 and len(host) == 1
    assert sorted(lay.sheets[lent[0]]) == ["A:98-100", "B:116-119"]
    panels = {a.get_gid(): a for a in d.figure.axes[0].patches if (a.get_gid() or "").startswith("sheet-panel:")}
    for k in lent + host:
        rects = np.array([lay.placed[i].rect for i in lay.sheets[k]] + [g.rect for g in block_ghosts(lay, k)])
        panel = panels[f"sheet-panel:{k}"]
        assert panel.get_width() <= rects[:, 2].max() - rects[:, 0].min() + 0.61  # its strands and ghosts, 0.3 pad


def test_a_ring_stalk_is_drawn_with_each_copy():
    """8UUP's C-terminal stalk (α9-α11, residues 456-485 of B, D and F) is one bundle of the three protomers. It is
    cut by protomer, so each copy's part is drawn with that copy, alike; the kinked pieces of each chain's helix
    stay one column, joined straight."""
    d = drawn(DATA / "8UUP.cif")
    lay = d.layout
    loops = {(lp.a_id, lp.b_id): lp for lp in d.loops}
    for c in "BDF":
        pieces = sorted(
            (p for p in lay.placed.values() if p.sse.chain == c and p.sse.first.seq >= 456), key=lambda p: p.sse.start
        )
        assert len(pieces) == 3
        assert all(abs(abs(p.direction @ pieces[0].direction) - 1) < 1e-9 for p in pieces)
        side = np.array([-pieces[0].direction[1], pieces[0].direction[0]])
        assert len({round(float(np.array([p.cx, p.cy]) @ side), 3) for p in pieces}) == 1  # one column
        for a, b in zip(pieces, pieces[1:]):
            lp = loops[(a.sse.id, b.sse.id)]
            assert not lp.fallback and len(lp.points) == 2


def test_a_cage_drawn_whole_has_its_copies_alike():
    """7A4M (an O cage of 24 chains) drawn whole: every tied copy alike (before the ring's bundles and lent elements
    were tied: 4.1 to 19.6 page units)."""
    got = figcheck.copies_alike(drawn(DATA / "assemblies" / "7A4M.cif", focus="none"))
    assert got and max(got.values()) <= ALIKE, {f"{a}~{b}": round(v, 3) for (a, b), v in got.items()}


def test_a_filament_drawn_whole_has_its_copies_alike():
    """6Y7S (FimA rod, six subunits) drawn whole: each subunit's sheet holds the next subunit's donor strand, and an
    item holding two subunits' elements belongs to the one holding most of it along a filament too, as round a ring
    (plan 10 Task 2 follow-up): every tied copy alike (5.2 to 6.7 page units before), the end subunit, which no
    neighbour lends a strand, within 1.2."""
    got = figcheck.copies_alike(drawn(DATA / "6Y7S.cif", focus="none"))
    shown = {f"{a}~{b}": round(v, 3) for (a, b), v in got.items()}
    assert len(got) == 5 and max(got.values()) <= 1.2, shown
    assert sorted(got.values())[-2] <= ALIKE, shown


def test_copies_of_one_domain_get_alike_panels():
    """A domain's panel counts a ghost only where it repeats a strand of its own block (a barrel's); a strand lent
    across the cut leaves its ghost in another copy's sheet, which is no part of this copy's domain. Taken on the
    default figure: with --domain, _part_domains still moves copies' domains apart one by one (Task 3)."""
    from foldmap.layout import domain_panel

    lay = drawn(DATA / "8UTF.cif").layout
    span = {"a": [], "c": []}
    for k, p in lay.placed.items():
        if p.sse.chain in span and 374 <= p.sse.first.seq and p.sse.last.seq <= 456:
            span[p.sse.chain].append(k)
    lay.domains = [("DIIa", span["a"]), ("DIIc", span["c"])]
    widths = {n: domain_panel(lay, n)[2] - domain_panel(lay, n)[0] for n in ("DIIa", "DIIc")}
    assert abs(widths["DIIa"] - widths["DIIc"]) < 1e-6, widths


@pytest.mark.parametrize("name", ["8UTF", "8UUP"])
def test_a_lent_strand_past_the_row_end_shows_the_strands_it_pairs_with(name):
    """A strand drawn past the row's end (it pairs in protomer 0's sheet, across the cut) stands beside dashed ghosts
    of its partners, one strand pitch across, their spans meeting, as it stands beside the real ones in every other
    copy (plan 12 Task 2: 8UTF's 440-442, off its sheet's track, stands in the next column, above the partner's
    centre)."""
    from foldmap.layout import PITCH, block_ghosts

    d = drawn(DATA / f"{name}.cif")
    lay = d.layout
    ghosted = {g.sse.id for g in lay.ghosts}
    lent = [k for k, ids in enumerate(lay.sheets) if set(ids) <= ghosted]
    assert lent
    for k in lent:
        partners = [g for g in block_ghosts(lay, k) if g.sse.id not in lay.sheets[k]]
        assert partners, lay.sheets[k]
        for i in lay.sheets[k]:
            p = lay.placed[i]
            near = min(
                abs(g.cx - p.cx) if abs(g.cy - p.cy) < (g.length + p.length) / 2 else np.inf
                for g in [*partners, *(lay.placed[j] for j in lay.sheets[k] if j != i)]
            )
            assert near <= PITCH + 1e-6, (i, near)
    assert figcheck.order_checks(d) == []


def test_an_item_of_two_protomers_belongs_to_the_one_holding_most_of_it():
    from foldmap.layout import Placed, _host
    from helpers import fake_sse

    def m(chain, a, b):
        return Placed(fake_sse("E", a, b, chain=chain), 0, 0, 1.0, 0.0, 1.0)

    sector, classes = {"A": 0, "B": 0, "E": 2, "F": 2}, {"A": 0, "B": 1, "E": 0, "F": 1}
    # 8UUP: A's and B's strands (7 residues) with F's three (8 residues) lent by protomer 2
    held = [m("A", 98, 100), m("B", 116, 119), m("F", 429, 431), m("F", 422, 424), m("F", 396, 397)]
    assert _host(held, sector, classes) == 2
    # a tie on residues: the first element by chain class, then residue number
    assert _host([m("B", 61, 62), m("E", 445, 446)], sector, classes) == 2
    assert _host([m("A", 61, 62), m("F", 445, 446)], sector, classes) == 0
    # two protomers' copies of one element first: symmetric in itself
    assert _host([m("A", 10, 20), m("E", 10, 20)], sector, classes) is None


def test_copies_numbered_apart_are_tied():
    """2MJW (random-100, revealed by C18): S100P chain D is numbered 95-188 where its copy B is 1-94. Copies were
    tied by shared residue numbers, so B and D shared none, were never tied, and were drawn unlike (A15 B~D 3.34).
    Numbers are taken to the first copy's by homolog_chains' offset before matching: D is tied to B, drawn alike."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(Path(__file__).parent / "data" / "r100" / "2MJW-BD.cif.gz"))
    assert d.layout.copy_ties == {"D": ("B", False)}
    got = figcheck.copies_alike(d)
    assert got[("D", "B")] <= 0.5, got
