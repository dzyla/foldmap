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
    copy. Never a ghost in its own column: one column means collinear, no pairing (fix round 4, the re-review: 8UTF's
    c:440-442, an accessory end to end with a:340-343 beside its partner A:35-57, took a:340-343's ghost, which pairs
    with nothing in its block)."""
    from foldmap.layout import PITCH, block_ghosts

    d = drawn(DATA / f"{name}.cif")
    lay = d.layout
    ghosted = {g.sse.id for g in lay.ghosts}
    lent = [k for k, ids in enumerate(lay.sheets) if set(ids) <= ghosted]
    assert lent
    for k in lent:
        partners = [g for g in block_ghosts(lay, k) if g.sse.id not in lay.sheets[k]]
        assert partners, lay.sheets[k]
        mine = [lay.placed[j] for j in lay.sheets[k]]
        for g in partners:  # each ghost stands beside a strand of the block, not end to end with one
            assert any(abs(g.cx - p.cx) > 1e-6 and abs(g.cy - p.cy) < (g.length + p.length) / 2 for p in mine), g.sse.id
        for i in lay.sheets[k]:
            p = lay.placed[i]
            near = min(  # across, spans meeting
                abs(g.cx - p.cx) if abs(g.cy - p.cy) < (g.length + p.length) / 2 and abs(g.cx - p.cx) > 1e-6 else np.inf
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
    Numbers are taken to the first copy's by homolog_chains' offset before matching: D is tied to B, drawn alike (as
    its mirror image: a two-fold assembly, plan 18)."""
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    d = draw_spec(figure_spec(Path(__file__).parent / "data" / "r100" / "2MJW-BD.cif.gz"))
    assert d.layout.copy_ties == {"D": ("B", True)}
    got = figcheck.copies_alike(d)
    assert got[("D", "B")] <= 0.5, got


@pytest.mark.parametrize("turned", [False, True])
def test_a_dimer_sheet_shared_by_two_copies_keeps_them_alike(monkeypatch, turned):
    """1SLV (plan 12 Task 2 fix round 2, pinned in fix round 4 on a cut fixture): assembly 1 is the A-B heterodimer
    twice, with one continuous sheet across the two copies, an item no tie holds. It sat between A1's eta2 and its
    copy, the tied pushes never cleared and _separate drew the copies one by one, unalike (A15). The shared items now
    yield to the tied copies: A2~A1 and B2~B1 are drawn exactly alike. The dyad is a half-turn whose axis may point
    either way (symmetry._rotation signs it): the copies are alike drawn either way up, and turned over against
    that sign the fixture still needs the shared items to yield (without that, b029264's parent: no `shared`, it is
    drawn unalike), so it holds the case whatever the sign. Drawn unrolled side by side, as before plan 18: in the
    dyad view (one copy and its mirror image) the shared sheet no longer stands between the copies."""
    import dataclasses

    import foldmap.cli as cli
    import foldmap.layout as layout

    monkeypatch.setattr(cli, "_dyad_frame", lambda *a: None)

    if turned:
        real_sym = cli.detect_symmetry
        monkeypatch.setattr(
            cli, "detect_symmetry", lambda *a, **k: (lambda s: dataclasses.replace(s, axis=-s.axis))(real_sym(*a, **k))
        )
    path = DATA / "r100" / "1SLV.cif.gz"
    got = figcheck.copies_alike(drawn(path))
    assert set(got) == {("A2", "A1"), ("B2", "B1")} and max(got.values()) <= 1e-6, got
    if not turned:
        return
    real = layout._separate

    def without_shared(*args, **kw):
        kw.pop("shared", None)
        return real(*args, **kw)

    monkeypatch.setattr(layout, "_separate", without_shared)
    monkeypatch.setattr(layout, "_tied_again", lambda *a, **k: False)  # nor the second chance (C2b-3), which yields
    before = figcheck.copies_alike(drawn(path))
    assert max(before.values()) > ALIKE, before  # the fixture still exercises the shared-sheet case


def _placed(path, monkeypatch=None, nudge=None):
    """Every placed element's centre and angle (cli.make_layout, as draw_spec lays it out); `nudge` (seed): each
    element axis moved by -1, 0 or +1 ulp per component first, as another CPU's BLAS kernel computes them."""
    import foldmap.cli as cli

    if nudge is not None:
        import dataclasses

        real = cli.build_sses
        rng = np.random.default_rng(nudge)

        def nudged(*args, **kw):
            out = real(*args, **kw)
            return [
                dataclasses.replace(s, axis=s.axis + rng.integers(-1, 2, 3) * np.spacing(np.abs(s.axis))) for s in out
            ]

        monkeypatch.setattr(cli, "build_sses", nudged)
    spec = figure_spec(path)
    lay = cli.make_layout(spec.structure, spec.mode, spec.rotate, spec.flip_v, **spec.layout_kwargs())[0]
    return {k: (p.cx, p.cy, p.angle) for k, p in sorted(lay.placed.items())}


@pytest.mark.parametrize("seed", range(4))
def test_a_sheet_seen_across_the_page_is_drawn_the_same_on_any_cpu(monkeypatch, seed):
    """1SLV's sheet shared by the two copies runs across the page: by symmetry its strands' mean direction is
    perpendicular to the dyad (page-up) and its end strands, dyad mates, stand level, so its height and their step
    across it are 0 up to the last bits of the element axes, which differ between BLAS kernels
    (OPENBLAS_CORETYPE=Nehalem, Windows, numpy 1.24's OpenBLAS). Their signs chose which way the strands point and
    in which order (layout._sheet_order), and with them the whole figure (on Windows the copies were drawn 1.49
    apart). The figure no longer depends on those bits."""
    path = DATA / "r100" / "1SLV.cif.gz"
    want = _placed(path)
    assert _placed(path, monkeypatch, nudge=seed) == want


def test_a_dyad_points_the_same_way_on_any_cpu(monkeypatch):
    """1SLV's dyad is a crystallographic half-turn: the copies' rotation has no antisymmetric part, and
    symmetry._rotation took the axis from np.linalg.eigh, whose sign LAPACK leaves free (OPENBLAS_CORETYPE=Nehalem,
    Windows and numpy 1.24's OpenBLAS turned it over, and the figure upside down). It is now signed by a rule of its
    own, whatever sign eigh returns."""
    import foldmap.cli as cli

    path = DATA / "r100" / "1SLV.cif.gz"
    found = []
    real = cli.detect_symmetry

    def kept(*args, **kw):
        found.append(real(*args, **kw))
        return found[-1]

    monkeypatch.setattr(cli, "detect_symmetry", kept)
    want = _placed(path)
    eigh = np.linalg.eigh
    monkeypatch.setattr(np.linalg, "eigh", lambda a, *args, **kw: (lambda w, v: (w, -v))(*eigh(a, *args, **kw)))
    assert _placed(path) == want
    one, other = found
    assert one.n == 2 and np.array_equal(one.axis, other.axis)


@pytest.mark.parametrize("name, layouts", [("1SLV", 1), ("2ILU", 2)])
def test_a_shared_sheet_whose_way_up_is_a_symmetry_tie_is_laid_out_both_ways(monkeypatch, name, layouts):
    """A sheet shared by two tied copies that runs straight across the page (_sheet_order) may point either way, and
    drawn one way the copies can take too much room to stay alike: 2ILU's two 20-strand sheets, strands[0] up, need
    1.07 times the one-by-one room (_SEP_SYM_ROOM 1.05) and its copies were drawn 3.1 apart. Such a layout is made
    again with those sheets turned down, and the turned one is kept when only it draws the copies alike (CI ruling,
    plan 12). 1SLV is alike the first way and is laid out once."""
    import foldmap.layout as layout

    real, runs = layout._attract, []

    def counted(*args, **kw):
        runs.append(1)
        return real(*args, **kw)

    monkeypatch.setattr(layout, "_attract", counted)
    got = figcheck.copies_alike(drawn(DATA / "r100" / f"{name}.cif.gz"))
    assert len(got) == 2 and max(got.values()) <= 1e-6, got
    assert len(runs) == layouts


@pytest.mark.parametrize("seed", [None, 0, 1])
def test_a_sheet_laid_out_both_ways_is_drawn_the_same_on_any_cpu(monkeypatch, seed):
    """2ILU's choice between its two layouts depends on no float noise: the element axes nudged by an ulp (seed) or
    eigh's eigenvectors negated (None) give the identical layout."""
    path = DATA / "r100" / "2ILU.cif.gz"
    want = _placed(path)
    if seed is None:
        eigh = np.linalg.eigh
        monkeypatch.setattr(np.linalg, "eigh", lambda a, *args, **kw: (lambda w, v: (w, -v))(*eigh(a, *args, **kw)))
        assert _placed(path) == want
    else:
        assert _placed(path, monkeypatch, nudge=seed) == want


def _laid_out(path, **kw):
    """The figure's layout and backbone (cli.make_layout, as draw_spec makes them; no loops routed or drawn): enough
    for figcheck.copies_alike, which measures element centres."""
    from types import SimpleNamespace

    from foldmap.cli import make_layout

    spec = figure_spec(path, **kw)
    lay, _, bb = make_layout(spec.structure, spec.mode, spec.rotate, spec.flip_v, **spec.layout_kwargs())
    return SimpleNamespace(layout=lay, bb=bb)


def test_crowded_copies_are_drawn_alike_in_more_room():
    """Random-100 C2b-3 A15, cause A (8QEM focus 26.2): the twelve copies of a C12 ring drawn whole cleared their
    overlaps together, but in 1.61 times the one-by-one room (_SEP_SYM_ROOM 1.05), so every item was moved alone and
    the copies drifted 2.3 to 26.2 page units apart. The tied pushes are kept, the copy-less items yielding to them,
    within _SEP_TIE_ROOM."""
    from helpers import r100_entry

    got = figcheck.copies_alike(_laid_out(r100_entry("8QEM"), focus="none"))
    assert len(got) == 22 and max(got.values()) <= ALIKE, {
        f"{a}~{b}": round(v, 3) for (a, b), v in got.items() if v > ALIKE
    }


@pytest.mark.parametrize(
    "entry, view, ties",
    [("9P1X", {"focus": "none"}, 20), ("7QUP", {"focus": "none"}, 60), ("6JLN", {}, 19)],
    ids=["9P1X-focus", "7QUP-focus", "6JLN-default"],
)
def test_copies_whose_elements_are_grouped_apart_are_drawn_alike(entry, view, ties):
    """Random-100 C2b-3 A15, cause B: DSSP joins and cuts the copies' sheets differently, so a copy's elements sit in
    items that differ from its mate's, and were tied in part or not at all. 9P1X: J's four strands and E:106-109 are
    a sheet of their own in protomer 0, one sheet with O, P and R in protomer 1 (Y~E 16.9, d~J 12.3); C's two strands
    share too few residues with V's for an item tie (C~V 10.0). 7QUP: 7D's sheet holds a strand 1D's lacks (7D~1D
    32.5). 6JLN: B:336-340 and B:430-433 are a sheet of their own, b's are in b's sheet (b~B 3.8); E's helices are
    bundled with R, a chain with no copy (e~E 18.3). A copy's matched elements are grouped as its mate's are
    (_regroup_copies), each group tied to the mate's item, so every copy is drawn alike."""
    from helpers import r100_entry

    d = _laid_out(r100_entry(entry), **view)
    got = figcheck.copies_alike(d)
    bad = {f"{a}~{b}": round(v, 3) for (a, b), v in got.items() if v > ALIKE}
    assert len(got) + len(figcheck.copies_not_judged(d)) == ties and not bad, bad


def test_what_a_regrouped_copy_holds_alone_stays_beside_it():
    """Random-100 C2b-3 VIEW: 7QUP's 7D:244-246, a strand 1D lacks, was left untied where its sheet had been
    projected, 42 page units from its copy drawn alike with 1D. Its sheet (which keeps its own shape, R2) is tied
    beside its mates' copy: the strand stands 4 page units from 7D's median element."""
    from helpers import r100_entry

    lay = _laid_out(r100_entry("7QUP"), focus="none").layout
    mine = np.array([(p.cx, p.cy) for p in lay.placed.values() if p.sse.chain == "7D"])
    lone = lay.placed["7D:244-246"]
    assert np.hypot(*(np.subtract((lone.cx, lone.cy), np.median(mine, axis=0)))) < 6.0  # 4.0 drawn


def test_a_copys_sheet_is_drawn_by_its_own_pairing():
    """C2b-3 fix round 1, ruling R2 (science outranks A15): _regroup_copies split 6JLN b's one DSSP sheet as B's two
    to draw b like B, which put b:353-356 beside b:398-402 with no bridge (S1) and drew panels no sheet of the
    bridges (A11). A sheet is never dissolved, split, merged or reshaped: S1-S3, S6 and A11 hold, the copies' differing
    sheets are left out of A15 and reported (figcheck.copies_sheets_differ), the rest drawn alike."""
    from helpers import r100_entry

    d = drawn(r100_entry("6JLN"))
    assert {k: v for k, v in figcheck.science_checks(d).items() if v} == {}
    assert figcheck.order_checks(d) == []
    assert "b" in {c for c, _ in figcheck.copies_sheets_differ(d)}
    got = figcheck.copies_alike(d)
    assert got and max(got.values()) <= ALIKE, {f"{a}~{b}": round(v, 3) for (a, b), v in got.items() if v > ALIKE}


def test_a_copys_barrel_keeps_its_ghost():
    """C2b-3 fix round 1, ruling R2: 9P1X G's sheet (G, H and I strands) closes into a barrel, a's does not;
    regrouping G's strands as a's dropped its closure (4787c55). The copy keeps its own: G:29-35 repeated as the
    ghost of the very block that holds it, as at 20bb525 (which never lost it)."""
    from helpers import r100_entry

    lay = _laid_out(r100_entry("9P1X"), focus="none").layout
    (n,) = [n for n, ids in enumerate(lay.sheets) if "G:29-35" in ids]
    assert lay.sheet_ghosts[n] == ["G:29-35"]
    assert {lay.placed[k].sse.chain for k in lay.sheets[n]} == {"G", "H", "I"}


def test_a_helix_piece_let_go_stays_with_its_helix():
    """C2b-3 fix round 1 (concern 3): 7QUP 9C:240-243, cut from 9C:224-238 where 1C:224-243 is one helix, was let go
    by _regroup_copies and drifted 59 to 62 page units from its copy. It joins its helix end to end (_on_end)."""
    from helpers import r100_entry

    lay = _laid_out(r100_entry("7QUP"), focus="none").layout
    piece, helix = lay.placed["9C:240-243"], lay.placed["9C:224-238"]
    assert np.hypot(piece.cx - helix.cx, piece.cy - helix.cy) < (piece.length + helix.length) / 2 + 1.0


@pytest.fixture(scope="module")
def p1x_focus_science():
    """9P1X --focus none (the whole file: 8 D3 protomers, 20 copy ties) laid out as draw_spec lays it out, and the
    sheet science checks of that layout (S1-S3, S6; no marks drawn, so S4 has none to judge): random-100 C2b-S."""
    from types import SimpleNamespace

    from foldmap.cli import make_layout
    from helpers import r100_entry

    spec = figure_spec(r100_entry("9P1X"), focus="none")
    lay, sses, bb = make_layout(spec.structure, spec.mode, spec.rotate, spec.flip_v, **spec.layout_kwargs())
    page = SimpleNamespace(axes=[SimpleNamespace(patches=[])])
    found = figcheck.science_checks(SimpleNamespace(layout=lay, sses=sses, bb=bb, figure=page))
    return lay, found


def _findings_naming(found: dict, *ids: str) -> list[str]:
    return [f"{k}: {t}" for k, got in found.items() for t in got or [] if any(i in t for i in ids)]


def test_an_accessory_never_stands_beside_a_strand_it_does_not_pair_with(p1x_focus_science):
    """Random-100 C2b-S, 9P1X focus: _column_heights stood an accessory where its own column was free, but beside
    strands of the next columns it does not pair with. a:106-109 (partner c:85-91, no room level with it) stood
    below c:2-8, beside c:50-61 (0 bridges); I:94-97, level with its partner H:85-97, stood beside I:2-8 and I:29-30
    (0 bridges); Q:106-109 is a:106-109's copy. An accessory's height must also keep it off every strand it does not
    pair with in the columns beside its own: where it cannot stand level with its partner so, it stands as near as
    it can and the notes name the pairing (S6, Layout.unshown)."""
    _, found = p1x_focus_science
    assert _findings_naming(found, "a:106-109", "Q:106-109", "I:94-97") == []


def test_a_copy_takes_its_mates_places_only_where_its_pairing_holds(p1x_focus_science):
    """Random-100 C2b-S, 9P1X focus: a copy tied in part (_symmetry_ties' else branch, DSSP grouped the copies' strands
    differently) took its mate's places for the strands they share (_take_shared), the rest moving by the mean step.
    L's sheet (N:16-17, L:110-111, L:74-79, L:142-147, L:172-176) took V's three strands in V's order, which runs the
    other way, so L:110-111 stood beside L:172-176 (0 bridges); O's barrel took E's places for R:29-35 and R:52-61
    with its ghost of O:29-36 between them, R:52-61 drawn away from its partner R:29-35 though the notes said beside.
    A copy takes its mate's places only when every strand it then draws side by side pairs, and every accessory
    drawn beside its partner stays so; else it keeps its own (ruling R2: S1-S3, S6 outrank A15)."""
    _, found = p1x_focus_science
    assert _findings_naming(found, "L:110-111", "R:52-61") == []


def test_a_block_another_chain_finds_grouped_apart_differs_for_every_chain_in_it():
    """Random-100 C2b-S, ruling C2b-3/R2 read per copy, not per chain: 9P1X protomer 1's O/P/R block is a 9-strand
    barrel (R:29-35 in its ring), protomer 0's E/F/J block an 8-strand one with J's strands beside it and J's others
    in a block of their own. O~E and R~J find the blocks grouped differently; P's strands lie in the same block, so
    they are placed by its own pairing too (the copy keeps its own shape: its mate's places drew R:52-61 away from
    R:29-35) and are left out of A15 and reported, not judged against F's (4.2 page units)."""
    from helpers import r100_entry

    d = _laid_out(r100_entry("9P1X"), focus="none")
    notes = figcheck.copies_sheets_differ(d)
    for copy, ref in (("P", "F"), ("Z", "F")):
        assert (copy, ref) in notes, sorted(notes)
    got = figcheck.copies_alike(d)
    assert max(got.values()) <= ALIKE, {f"{a}~{b}": round(v, 3) for (a, b), v in got.items() if v > ALIKE}


def test_a_filaments_copies_keep_their_own_strands_beside_their_partners():
    """Random-100 C2b-S: 6Y7S drawn whole. Each copy's sheet takes its mate's places for the strands they share,
    which run the other way across the page; the mean step left B:122-123 (which A's sheet lacks) beside B:29-38,
    0 bridges (and C's, D's, E's likewise). It turns with the shared strands, beside B:68-74, which it pairs with."""
    found = figcheck.science_checks(drawn(DATA / "6Y7S.cif", focus="none"))
    assert {k: v for k, v in found.items() if v} == {}


def test_a_copy_given_its_mates_place_may_name_a_pairing_its_mate_names():
    """Random-100 C2b-S: 4XG1's B barrel takes A's places one to one; B:174-175 then stands where A:174-175 does, not
    level with its partner B:145-152, and the notes name the pairing, as A's do (S6's ruled case). That is no reason
    to keep B's own shape (2.36 page units from A's). Only a strand that merely moves with the taken ones must keep
    every partner it stood beside (9P1X R:2-8)."""
    from helpers import r100_entry

    d = _laid_out(r100_entry("4XG1"))
    got = figcheck.copies_alike(d)
    assert got and max(got.values()) <= ALIKE, {f"{a}~{b}": round(v, 3) for (a, b), v in got.items()}


def _strand_item(chain: str, strands: list[tuple[int, int, float]]):
    """A sheet item of upright strands (start, end, x), 2 long and 1 wide, all at height 0."""
    from foldmap.layout import Placed, _Item
    from helpers import fake_sse

    members = [Placed(fake_sse("E", a, b, chain=chain), 0.0, 0.0, 2.0, np.pi / 2, 1.0) for a, b, _ in strands]
    return _Item(members, [(x, 0.0) for _, _, x in strands], 0.0, 0.0, (0.0, 0.0), strands[0][0], chain)


def test_a_copy_names_a_pairing_by_taking_its_mates_place_only_where_its_mate_names_it():
    """C2b-S fix round 1 (review minor 1): a strand given its mate's place may leave a partner the notes then name
    only when the mate's own drawing names the matching pairing (4XG1 B:174-175, as A:174-175). Here B:20-22 pairs
    with B:1-7 (its partner) and B:10-16, A:20-22 with A:10-16 alone, drawn level beside it: A's places would leave
    B:20-22 away from B:1-7 where A names nothing, so B keeps its own drawing, B:20-22 beside B:1-7."""
    from types import SimpleNamespace

    from foldmap.layout import _symmetry_ties

    mate = _strand_item("A", [(1, 7, 0.0), (10, 16, 1.1), (20, 22, 2.2)])
    copy = _strand_item("B", [(1, 7, 0.0), (10, 16, 1.1), (20, 22, -1.1)])
    frame = SimpleNamespace(sector={"A": 0, "B": 1}, classes={"A": "p", "B": "p"}, mirror_x=None, rise=0.0)
    pairs = {frozenset(p) for p in [("A:1-7", "A:10-16"), ("A:10-16", "A:20-22")]}
    pairs |= {frozenset(p) for p in [("B:1-7", "B:10-16"), ("B:1-7", "B:20-22"), ("B:10-16", "B:20-22")]}
    partner = {"A:20-22": "A:10-16", "B:20-22": "B:1-7"}
    own: set[int] = set()
    ties = _symmetry_ties([mate, copy], frame, False, None, (pairs, partner), own)
    assert ties == {1: (0, 1)}
    assert copy.offsets == [(0.0, 0.0), (1.1, 0.0), (-1.1, 0.0)] and id(copy) in own


def test_an_element_the_mate_lacks_is_slid_off_the_shared_ones():
    """Random-100 held-out 2NAD (A3, A5): chain B's βα row holds a helix, B:267-270, between two strands where chain
    A's row has none. Taking A's places for the shared strands (_take_shared) stood them side by side and the helix,
    moved by the mean step, between them, on both: its two loops had no route. It is slid off them, clear by MARGIN,
    and the shared strands keep the mate's places."""
    from foldmap.layout import HELIX_W, PITCH, Placed, _clear_of, _Item, _take_shared
    from helpers import fake_sse

    def member(kind, start, chain):
        width = PITCH if kind == "E" else HELIX_W
        return Placed(fake_sse(kind, start, start + 5, chain=chain), 0.0, 0.0, 3.0, np.pi / 2, width)

    mate = _Item([member("E", 0, "A"), member("E", 20, "A")], [(0.0, 0.0), (1.4, 0.0)], 0.0, 0.0, (0.0, 0.0), 0, "A")
    copy = _Item(
        [member("E", 0, "B"), member("H", 10, "B"), member("E", 20, "B")],
        [(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)],
        0.0,
        0.0,
        (0.0, 0.0),
        0,
        "B",
    )
    _take_shared(copy, mate, [(0, 0), (2, 1)])
    assert copy.offsets[0] == (0.0, 0.0) and copy.offsets[2] == (1.4, 0.0)
    others = [(m, o) for k, (m, o) in enumerate(zip(copy.members, copy.offsets)) if k != 1]
    assert _clear_of(copy.members[1], copy.offsets[1], others)


@pytest.mark.parametrize("view", [{}, {"domains": "auto"}], ids=["default", "domains"])
def test_2nad_copy_helix_stands_off_its_rows_strands(view):
    """Random-100 held-out 2NAD: no element on another, every loop routed."""
    from helpers import r100_entry

    d = draw_spec(figure_spec(r100_entry("2NAD"), **view))
    assert d.stuck == []
    assert figcheck.element_overlaps(d.layout) == []
