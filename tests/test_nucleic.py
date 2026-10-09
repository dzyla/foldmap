from pathlib import Path

import numpy as np
import pytest

from foldmap.io import load_nucleic
from helpers import pipeline

DATA = Path(__file__).parent / "data"


def test_protein_only_file_has_no_nucleic_acid(ubq):
    bb, _, _ = pipeline(ubq)
    na = load_nucleic(ubq, bb)
    assert na.strands == [] and na.pairs == [] and na.contacts == {}


@pytest.mark.parametrize("pdb, n_nt, min_pairs", [("1ZAA", 11, 10), ("1LMB", 20, 18)])
def test_duplex_strands_and_watson_crick_pairs(ubq, pdb, n_nt, min_pairs):
    path = ubq.parent / f"{pdb}.cif"
    bb, _, _ = pipeline(path)
    na = load_nucleic(path, bb)
    assert len(na.strands) == 2
    assert all(len(s) == n_nt and not s.rna for s in na.strands)
    assert all(s.p.shape == (len(s), 3) and s.c1.shape == (len(s), 3) for s in na.strands)
    assert len(na.pairs) >= min_pairs
    sums = {i + j for a, i, b, j in na.pairs}
    assert all((a, b) == (0, 1) for a, i, b, j in na.pairs)
    assert max(sums) - min(sums) <= 2  # antiparallel: i on one strand meets n-1-i on the other (± overhang)
    used = [(a, i) for a, i, _, _ in na.pairs] + [(b, j) for _, _, b, j in na.pairs]
    assert len(used) == len(set(used))  # every nucleotide pairs at most once


def test_each_zinc_finger_recognition_helix_touches_the_dna(ubq):
    path = ubq.parent / "1ZAA.cif"
    bb, sses, _ = pipeline(path)
    na = load_nucleic(path, bb)
    helices = [s for s in sses if s.kind == "H"]
    assert len(helices) == 3
    for h in helices:
        assert any(r in na.contacts for r in range(h.start, h.end + 1)), h.id
    for r, nts in na.contacts.items():
        assert 0 <= r < len(bb) and all(0 <= s < 2 and 0 <= k < len(na.strands[s]) for s, k in nts)


def test_strand_ids_and_ends(ubq):
    path = ubq.parent / "1ZAA.cif"
    bb, _, _ = pipeline(path)
    na = load_nucleic(path, bb)
    assert [s.id for s in na.strands] == ["A:1-11", "B:12-22"] or all(":" in s.id for s in na.strands)
    for s in na.strands:
        gaps = np.linalg.norm(np.diff(s.p, axis=0), axis=1)
        assert gaps.max() < 8.0  # one continuous strand, 5' -> 3'


def _dna_layout(path, mode, with_dna=True):
    from foldmap.features import sse_contacts
    from foldmap.frame import view_frame
    from foldmap.layout import build_layout

    bb, sses, sheets = pipeline(path)
    na = load_nucleic(path, bb)
    from foldmap.frame import dna_radial_frame

    f = dna_radial_frame(sses, na) if with_dna and na.strands else view_frame(sses)
    lay = build_layout(sses, sheets, f, mode, contacts=sse_contacts(bb, sses), nucleic=na if with_dna else None)
    return bb, sses, na, f, lay


def _overlap(a, b):
    return min(a[2], b[2]) - max(a[0], b[0]) > 1e-6 and min(a[3], b[3]) - max(a[1], b[1]) > 1e-6


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
@pytest.mark.parametrize("mode", ["projected", "stack"])
def test_duplex_is_placed_true_to_scale_and_clear_of_the_protein(ubq, pdb, mode):
    from foldmap.layout import SCALE

    _, _, na, f, lay = _dna_layout(ubq.parent / f"{pdb}.cif", mode)
    (d,) = lay.dna
    assert sorted(d.strands) == [0, 1]
    n_bp = max(len(s) for s in na.strands)
    assert d.length == pytest.approx(n_bp * 3.4 * SCALE, abs=1.5)  # B-DNA rises 3.4 Å per base pair
    assert d.width == pytest.approx(20.0 * SCALE)  # and is 20 Å across
    for p in [*lay.placed.values(), *lay.ghosts]:
        assert not _overlap(d.rect, p.rect), p.sse.id
    x0, y0, x1, y1 = lay.bounds
    assert x0 <= d.rect[0] and d.rect[2] <= x1 and y0 <= d.rect[1] and d.rect[3] <= y1


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_duplex_runs_along_its_projected_axis_with_strands_antiparallel(ubq, pdb):
    _, _, na, f, lay = _dna_layout(ubq.parent / f"{pdb}.cif", "projected")
    (d,) = lay.dna
    pts = np.vstack([s.p for s in na.strands])
    axis = np.linalg.svd(pts - pts.mean(axis=0))[2][0]
    a2 = f.project(f.origin + axis)[0]
    if np.linalg.norm(a2) > 0.3:
        assert abs(a2 @ d.direction) / np.linalg.norm(a2) > 0.99
    s0 = [d.axial[(0, i)] for i in range(len(na.strands[0]))]
    s1 = [d.axial[(1, i)] for i in range(len(na.strands[1]))]
    assert np.all(np.diff(s0) > 0) and np.all(np.diff(s1) < 0)  # 5' -> 3' run opposite ways
    for a, i, b, j in na.pairs:
        assert abs(d.axial[(a, i)] - d.axial[(b, j)]) < 0.6  # partners sit across from each other


def _gaps(path):
    """(chain, touches the DNA, box gap to the duplex) per element."""
    bb, sses, na, _, lay = _dna_layout(path, "projected")
    (d,) = lay.dna

    def gap(a, b):  # between boxes; 0 when they touch
        return float(np.hypot(max(a[0] - b[2], b[0] - a[2], 0), max(a[1] - b[3], b[1] - a[3], 0)))

    return [
        (s.chain, any(r in na.contacts for r in range(s.start, s.end + 1)), gap(lay.placed[s.id].rect, d.rect))
        for s in sses
    ]


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_every_protein_chain_reaches_the_dna(ubq, pdb):
    from foldmap.layout import HELIX_W, MARGIN

    rows = _gaps(ubq.parent / f"{pdb}.cif")
    for chain in {c for c, _, _ in rows}:
        assert min(g for c, t, g in rows if c == chain and t) <= MARGIN + HELIX_W, chain  # within a helix width


def test_stack_mode_puts_the_dna_in_its_own_row_on_top(ubq):
    _, _, _, _, lay = _dna_layout(ubq.parent / "1LMB.cif", "stack")
    (d,) = lay.dna
    assert abs(np.sin(d.angle)) < 1e-9
    assert d.rect[1] > max(p.rect[3] for p in lay.placed.values())


def test_protein_only_layout_has_no_dna(ubq):
    _, _, _, _, lay = _dna_layout(ubq, "projected")
    assert lay.dna == []


def _figure(path, **kw):
    from foldmap.cli import make_figure_and_loops

    return make_figure_and_loops(path, **kw)


def _gid(fig, prefix):
    return [
        a
        for a in [*fig.axes[0].patches, *fig.axes[0].texts, *fig.axes[0].lines, *fig.axes[0].collections]
        if (a.get_gid() or "").startswith(prefix)
    ]


def test_duplex_is_drawn_as_a_double_helix_from_one_repeated_turn(ubq):
    fig, _ = _figure(ubq.parent / "1LMB.cif")
    (front,) = _gid(fig, "dna-front:")
    (back,) = _gid(fig, "dna-back:")
    polys = [p for p in front.get_path().to_polygons(closed_only=False) if len(p) > 2]
    assert len(polys) >= 4  # two backbones, ~2 turns for 20 bp, each turn its own pair of front faces
    sizes = {len(p) for p in polys}
    assert len(sizes) <= 3  # tiles of one fragment (full turns, plus clipped end pieces)
    assert len(_gid(fig, "dna-rung:")) == 1  # base pairs: one compound path


def test_dna_ends_are_labelled_five_and_three_prime(ubq):
    fig, _ = _figure(ubq.parent / "1ZAA.cif")
    ends = sorted(t.get_text() for t in _gid(fig, "dna-end:"))
    assert ends == ["3′", "3′", "5′", "5′"]


def test_contacts_are_marked_in_the_protein_chain_colour(ubq):
    from matplotlib.colors import to_hex

    from foldmap.palette import chain_colors

    fig, _ = _figure(ubq.parent / "1LMB.cif")
    marks = _gid(fig, "dna-contact:")
    assert marks
    colours = {to_hex(c) for m in marks for c in m.get_facecolors()}
    assert colours <= {c.lower() for c in chain_colors(["3", "4"]).values()}  # protein chains 3 and 4
    assert _gid(fig, "dna-tether:")


def test_loops_route_around_the_duplex(ubq):
    from foldmap.cli import make_layout
    from foldmap.route import route_loops
    from test_route import segments_clear

    for pdb in ("1ZAA", "1LMB"):
        lay, sses, bb = make_layout(ubq.parent / f"{pdb}.cif")
        (d,) = lay.dna
        for l in route_loops(lay, sses, bb):
            if not l.fallback:
                assert segments_clear(l.points, [d.rect]), (pdb, l.a_id, l.b_id)


def test_no_dna_option_hides_it(ubq, tmp_path):
    from foldmap.cli import main

    out = tmp_path / "x.svg"
    assert main(["plot", str(ubq.parent / "1ZAA.cif"), "-o", str(out), "--no-dna"]) == 0
    assert "dna-front:" not in out.read_text()
    assert main(["plot", str(ubq.parent / "1ZAA.cif"), "-o", str(out)]) == 0
    assert "dna-front:" in out.read_text()


def test_summary_lists_the_nucleic_acid(ubq, capsys):
    from foldmap.cli import main

    assert main(["summary", str(ubq.parent / "1LMB.cif")]) == 0
    out = capsys.readouterr().out
    assert "DNA" in out and "19 base pairs" in out and "contacts" in out


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_tethers_join_an_element_to_its_nearest_contact_mark(ubq, pdb):
    from foldmap.cli import make_layout

    fig, _ = _figure(ubq.parent / f"{pdb}.cif")
    lay, _, _ = make_layout(ubq.parent / f"{pdb}.cif")
    marks = np.vstack([m.get_offsets() for m in _gid(fig, "dna-contact:")])
    for line in _gid(fig, "dna-tether:"):
        p = lay.placed[line.get_gid().split(":", 1)[1]]
        (x0, x1), (y0, y1) = line.get_xdata(), line.get_ydata()
        a, b = np.array(p.n_port), np.array(p.c_port)
        t = np.clip((np.array([x0, y0]) - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1)
        assert np.hypot(*(a + t * (b - a) - [x0, y0])) < 1e-6  # starts on the element's axis
        assert np.min(np.hypot(*(marks - [x1, y1]).T)) < 1e-6  # ends on a contact mark
        own = _contact_points(lay, p.sse)
        ab = np.linspace(a, b, 400)
        nearest = min(np.min(np.hypot(*(ab - q).T)) for q in own)
        assert np.hypot(x1 - x0, y1 - y0) == pytest.approx(nearest, abs=0.02)  # the shortest way to its contacts


def _contact_points(lay, sse):
    from foldmap.layout import dna_to_page, dna_y

    (d,) = lay.dna
    na = lay.nucleic
    slot = {k: n % 2 for n, k in enumerate(d.strands)}
    hit = {h for r in range(sse.start, sse.end + 1) for h in na.contacts.get(r, ()) if h[0] in slot}
    return [dna_to_page(d, d.axial[h], dna_y(slot[h[0]], d.axial[h])) for h in hit]


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
@pytest.mark.parametrize("mode", ["projected", "stack"])
def test_tethers_never_cross_other_elements(ubq, pdb, mode):
    from foldmap.cli import make_layout
    from test_route import segments_clear

    fig, _ = _figure(ubq.parent / f"{pdb}.cif", mode=mode)
    lay, _, _ = make_layout(ubq.parent / f"{pdb}.cif", mode)
    for line in _gid(fig, "dna-tether:"):
        own = line.get_gid().split(":", 1)[1]
        others = [p.rect for k, p in lay.placed.items() if k != own] + [g.rect for g in lay.ghosts]
        assert segments_clear(list(zip(line.get_xdata(), line.get_ydata())), others), (mode, own)


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_dna_complexes_are_viewed_with_the_dna_across_and_the_protein_below(ubq, pdb):
    from foldmap.cli import make_layout

    lay, _, _ = make_layout(ubq.parent / f"{pdb}.cif")
    (d,) = lay.dna
    assert abs(np.sin(d.angle)) < 0.05  # the duplex runs across the page
    ys = [p.cy for p in lay.placed.values()]
    assert np.mean(ys) < d.cy  # the protein hangs below it


def test_dna_view_is_a_rotation_not_a_mirror(ubq):
    from foldmap.frame import dna_view_frame

    path = ubq.parent / "1LMB.cif"
    bb, sses, _ = pipeline(path)
    na = load_nucleic(path, bb)
    f = dna_view_frame(sses, na)
    m = np.array([f.u, f.v, f.w])
    assert np.allclose(m @ m.T, np.eye(3), atol=1e-9) and np.linalg.det(m) == pytest.approx(1.0)


def test_every_nucleotide_is_labelled_with_its_base(ubq):
    fig, _ = _figure(ubq.parent / "1ZAA.cif")
    bb, _, _ = pipeline(ubq.parent / "1ZAA.cif")
    na = load_nucleic(ubq.parent / "1ZAA.cif", bb)
    labels = _gid(fig, "nt:")
    assert len(labels) == sum(len(s) for s in na.strands)
    want = sorted(l.name.strip()[-1] for s in na.strands for l in s.labels)  # DA -> A
    assert sorted(t.get_text() for t in labels) == want


def test_nucleotide_labels_can_be_switched_off(ubq):
    from foldmap.style import Style

    fig, _ = _figure(ubq.parent / "1ZAA.cif", look=Style(nucleotide_labels=False))
    assert not _gid(fig, "nt:")


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_every_protein_element_hangs_below_the_dna(ubq, pdb):
    from foldmap.cli import make_layout

    lay, _, _ = make_layout(ubq.parent / f"{pdb}.cif")
    (d,) = lay.dna
    assert all(p.rect[3] <= d.rect[1] + 1e-6 for p in lay.placed.values())  # nothing above or beside it


def test_zinc_fingers_are_compact_under_their_sites(ubq):
    from foldmap.cli import make_layout

    lay, sses, _ = make_layout(ubq.parent / "1ZAA.cif")
    (d,) = lay.dna
    depth = d.rect[1] - min(p.rect[1] for p in lay.placed.values())
    assert depth <= 3 * d.width  # three ββα fingers fit in a band about three duplex widths deep
    x0, x1 = d.rect[0], d.rect[2]
    assert all(x0 - 2 <= p.cx <= x1 + 2 for p in lay.placed.values())  # under the DNA, not off to the side


def test_radial_dna_frame_is_a_rotation_everywhere(ubq):
    from foldmap.frame import dna_radial_frame

    path = ubq.parent / "1ZAA.cif"
    bb, sses, _ = pipeline(path)
    na = load_nucleic(path, bb)
    f = dna_radial_frame(sses, na)
    for s in sses:
        m = np.array(f.basis(s.centroid))
        assert np.allclose(m @ m.T, np.eye(3), atol=1e-9) and np.linalg.det(m) == pytest.approx(1.0)


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_stack_mode_keeps_the_protein_close_and_the_dna_over_its_binders(ubq, pdb):
    from foldmap.cli import make_layout
    from foldmap.layout import MARGIN, label_boxes

    lay, sses, bb = make_layout(ubq.parent / f"{pdb}.cif", "stack")
    na = load_nucleic(ubq.parent / f"{pdb}.cif", bb)
    (d,) = lay.dna
    protein = [p.rect for p in lay.placed.values()] + [b for n, b in label_boxes(lay, sses) if not n.startswith("dna")]
    top = max(r[3] for r in protein)  # the protein row's top, its labels included
    assert d.rect[1] - top <= MARGIN + 0.6  # the next row down, not a row and a half
    binders = [lay.placed[s.id].cx for s in sses if any(r in na.contacts for r in range(s.start, s.end + 1))]
    assert abs(d.cx - np.mean(binders)) < 1.0  # centred over what it touches


@pytest.mark.parametrize("pdb", ["1ZAA", "1LMB"])
def test_contact_beads_leave_the_base_letters_readable(ubq, pdb):
    fig, _ = _figure(ubq.parent / f"{pdb}.cif")
    ax = fig.axes[0]
    unit_px = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    px_per_pt = fig.dpi / 72
    letters = np.array([t.get_position() for t in _gid(fig, "nt:")])
    glyph = 0.5 * max(t.get_fontsize() for t in _gid(fig, "nt:")) * px_per_pt / unit_px  # half a letter
    for beads in _gid(fig, "dna-contact:"):
        r = np.sqrt(beads.get_sizes()[0]) / 2 * px_per_pt / unit_px  # marker radius, page units
        for b in beads.get_offsets():
            assert np.min(np.hypot(*(letters - b).T)) >= r + glyph * 0.8


def test_empty_turn_list_gives_an_empty_path():
    from foldmap.render import _compound

    assert len(_compound([]).vertices) == 0 and len(_compound([], closed=False).vertices) == 0


def _touched(na):
    out = {}
    for v in na.contacts.values():
        for k, i in v:
            out.setdefault(k, set()).add(i)
    return out


def test_dna_cropped_to_the_contacted_stretch():
    from foldmap.io import load_backbone, load_nucleic
    from foldmap.model import crop_nucleic

    path = DATA / "1LMB.cif"
    bb = load_backbone(path, "asu")
    na = load_nucleic(path, bb, "asu")
    cut = crop_nucleic(na, flank=0)
    assert len(cut.strands) == len(na.strands)
    for k, s in enumerate(cut.strands):
        assert len(s) <= len(na.strands[k])
        assert s.labels[0].seq >= na.strands[k].labels[0].seq
    for a, i, b, j in cut.pairs:  # indices stay valid and still name the same nucleotides
        assert 0 <= i < len(cut.strands[a]) and 0 <= j < len(cut.strands[b])
    old = {(na.strands[a].labels[i], na.strands[b].labels[j]) for a, i, b, j in na.pairs}
    assert {(cut.strands[a].labels[i], cut.strands[b].labels[j]) for a, i, b, j in cut.pairs} <= old
    assert sum(len(v) for v in _touched(cut).values()) == sum(len(v) for v in _touched(na).values())
    assert any(s.cut5 or s.cut3 for s in cut.strands)


def test_untouched_duplexes_are_dropped():
    from foldmap.io import load_backbone, load_nucleic
    from foldmap.model import crop_nucleic

    path = DATA / "1LMB.cif"
    bb = load_backbone(path, "asu")
    na = load_nucleic(path, bb, "asu")
    na.contacts = {}
    assert crop_nucleic(na).strands == []


def test_cropped_ends_are_marked_and_full_dna_on_request():
    from foldmap.cli import make_figure

    texts = lambda fig: [t.get_text() for t in fig.axes[0].texts]  # noqa: E731
    from foldmap.style import Style

    cropped = make_figure(DATA / "1LMB.cif", look=Style(dna_extent="contacts"))
    full = make_figure(DATA / "1LMB.cif", look=Style(dna_extent="all"))
    assert "5′" in texts(full) and texts(full).count("5′") >= 2
    assert len(texts(cropped)) <= len(texts(full))
