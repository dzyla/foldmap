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
def test_symmetric_controls_unchanged(path):
    got = figcheck.copies_alike(drawn(path))
    assert all(math.isfinite(v) and v <= 1e-6 for v in got.values()), got


def test_a_strand_lent_across_the_cut_leaves_a_ghost_where_it_pairs():
    """8UTF: each protomer lends its S and T strands (440-442, 445-446) to the next protomer's sheets. Unrolled into a
    row a, b, c, chain c's strands would pair in a's sheets at the far end; they are drawn past c instead, where its
    copies' strands sit, and a's sheets keep a ghost of each. The order check reads a sheet with its ghosts."""
    from foldmap.layout import block_ghosts

    d = drawn(DATA / "8UTF.cif")
    lay = d.layout
    assert sorted(g.sse.id for g in lay.ghosts) == ["c:440-442", "c:445-446"]
    hosts = {g.sse.id: lay.sheets[k] for k in range(len(lay.sheets)) for g in block_ghosts(lay, k)}
    assert "a:325-329" in hosts["c:440-442"] and hosts["c:445-446"] == ["A:61-62"]
    right = max(p.rect[2] for k, p in lay.placed.items() if p.sse.chain in "Cc" and k not in hosts)
    assert all(lay.placed[k].cx > right for k in hosts)  # past the rest of copy c
    assert figcheck.order_checks(d) == []


def test_a_sheet_panel_covers_its_own_ghosts_only():
    """8UUP: chain F's three strands lent across the cut stand past the row's end as a block of their own; its panel
    covers them, not their ghosts in chain B's sheet at the other end of the page."""
    from foldmap.layout import block_ghosts

    d = drawn(DATA / "8UUP.cif")
    lay = d.layout
    ghosted = {g.sse.id for g in lay.ghosts}
    lent = [k for k, ids in enumerate(lay.sheets) if set(ids) <= ghosted]
    host = [k for k in range(len(lay.sheets)) if block_ghosts(lay, k)]
    assert len(lent) == 1 and len(host) == 1 and not block_ghosts(lay, lent[0])
    assert sorted(lay.sheets[lent[0]]) == ["F:396-397", "F:422-424", "F:429-431"]
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
