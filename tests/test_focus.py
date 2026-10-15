"""Large assemblies: one subunit in full, with the parts of its neighbours that complete its fold."""

from pathlib import Path

import pytest

from foldmap.cli import main, make_figure, make_layout
from foldmap.layout import termini

DATA = Path(__file__).parent / "data"
PILUS = DATA / "6Y7S.cif"  # type 1 pilus rod (FimA), helical, six subunits deposited


def _chains(sses):
    return {s.chain for s in sses}


def test_helical_filament_shows_its_central_subunit():
    lay, sses, _ = make_layout(PILUS)
    focus = lay.focus_chains
    assert len(focus) == 1 and next(iter(focus)) in "CD"  # the middle of A..F: neighbours on both sides
    own = [s for s in sses if s.chain in focus]
    full, _, _ = make_layout(PILUS, focus="none")
    assert len(own) == sum(1 for p in full.placed.values() if p.sse.chain in focus)  # every own element
    assert len(sses) - len(own) <= 6  # only the neighbours' fold-completing elements


def test_neighbours_donor_strand_completes_the_fold():
    """FimA's Ig fold lacks a strand; the next subunit's N-terminal donor strand fills the groove."""
    lay, sses, _ = make_layout(PILUS)
    focus = next(iter(lay.focus_chains))
    blocks = [ids for ids in lay.sheets if any(k.startswith(focus + ":") for k in ids)]
    assert any(any(not k.startswith(focus + ":") for k in ids) for ids in blocks)  # a sheet with a guest strand
    guests = [p for p in lay.placed.values() if p.sse.chain != focus]
    assert guests and all(p.label.endswith(("′", "″")) for p in guests)


def test_fragments_of_neighbours_have_no_termini():
    lay, sses, _ = make_layout(PILUS)
    assert {chain for _, chain, _, _ in termini(lay, sses)} == lay.focus_chains


def test_figure_says_what_it_shows():
    fig = make_figure(PILUS)
    texts = " ".join(t.get_text() for t in fig.axes[0].texts)
    assert "1 of 6 subunits" in texts and "helical" in texts and "115" in texts
    assert "neighbour" in texts


def test_everything_on_request():
    lay, sses, _ = make_layout(PILUS, focus="none")
    assert lay.focus_chains is None and _chains(sses) == set("ABCDEF")


def test_chosen_subunit():
    lay, sses, _ = make_layout(PILUS, focus="B")
    assert lay.focus_chains == {"B"}


@pytest.mark.parametrize("name", ["8UTF", "1BL8", "2HHB", "1UBQ"])
def test_modest_assemblies_are_drawn_whole(name):
    lay, _, _ = make_layout(DATA / f"{name}.cif")
    assert lay.focus_chains is None


def test_cubic_cage_focuses_on_one_protomer():
    lay, sses, _ = make_layout(DATA / "assemblies" / "7A4M.cif")
    assert lay.focus_chains is not None and len(_chains(sses)) < 24


def test_cli_and_layout_file(tmp_path):
    import yaml

    out, lf = tmp_path / "p.svg", tmp_path / "p.yaml"
    assert main(["plot", str(PILUS), "-o", str(out), "--focus", "none", "--save-layout", str(lf)]) == 0
    assert yaml.safe_load(lf.read_text())["layout"]["focus"] == "none"


def test_only_the_focus_subunit_is_traced():
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(PILUS)
    loops = route_loops(lay, sses, bb)
    assert loops and all(lay.placed[l.a_id].sse.chain in lay.focus_chains for l in loops)


def test_neighbour_fragments_join_only_where_the_chain_really_does():
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(PILUS)
    full, fsses, _ = make_layout(PILUS, focus="none")
    order = {s.id: k for k, s in enumerate(fsses)}
    for loop in route_loops(lay, sses, bb):
        assert order[loop.b_id] == order[loop.a_id] + 1, (loop.a_id, loop.b_id)  # no element skipped


def test_links_of_neighbour_fragments_are_left_out():
    fig = make_figure(PILUS)
    lay, _, _ = make_layout(PILUS)
    for a in fig.axes[0].patches:
        gid = a.get_gid() or ""
        if gid.startswith(("disulfide:", "glycan:")):
            i = int(gid.split(":")[1].split("-")[0])
            assert lay.res_chain[i] in lay.focus_chains, gid


def test_one_copy_of_each_distinct_chain():
    """Apoferritin is 24 identical chains: the unit is one chain, and the figure says '1 of 24'."""
    lay, sses, _ = make_layout(DATA / "assemblies" / "7A4M.cif")
    assert len(lay.focus_chains) == 1
    assert "1 of 24" in lay.context


def test_heteromeric_unit_keeps_one_of_each_partner():
    lay, _, _ = make_layout(DATA / "8UUP.cif", focus="protomer")  # C3 of A,B-like heterodimers
    assert len(lay.focus_chains) == 2


def test_ligands_only_when_they_touch_the_subunit():
    path = DATA / "assemblies" / "7A4M.cif"
    fig = make_figure(path)
    lay, _, _ = make_layout(path)
    markers = [a for a in fig.axes[0].patches if (a.get_gid() or "").startswith("ligand:")]
    touching = [
        g for g in lay.links.ligands if not g.additive and any(lay.res_chain[r] in lay.focus_chains for r in g.contacts)
    ]
    assert len(markers) == len(touching)


CAPSID = DATA / "9TX7.cif"  # MPMV capsid 3-fold interface: full CA chains A, C, E; CTD-only chains B, D, F


def test_copies_of_one_protein_share_element_labels():
    lay, sses, _ = make_layout(CAPSID, focus="none")
    label = {(p.sse.chain, p.sse.first.seq): p.label for p in lay.placed.values()}
    assert label[("B", 156)] == label[("A", 156)] == "α7"  # B starts at the CTD: its first helix is still α7
    assert label[("D", 206)] == label[("A", 206)]


def test_no_neighbour_legend_when_no_neighbour_is_drawn():
    lay, _, _ = make_layout(CAPSID, focus="A")
    if not lay.partial:
        fig = make_figure(CAPSID, focus="A")
        texts = " ".join(t.get_text() for t in fig.axes[0].texts)
        assert "neighbour" not in texts and "1 of" in texts


NUCLEOSOME = DATA / "assemblies" / "3UTB.cif.gz"  # C2 histone octamer wrapped by one DNA duplex


def test_focus_auto_never_reduces_a_two_subunit_assembly():
    """The histone octamer without its DNA: C2, 8 protein chains (more than FOCUS_CHAINS), two subunits. focus
    auto draws both halves; asked for, the protomer is still one half."""
    lay, sses, _ = make_layout(NUCLEOSOME, dna=False)
    assert lay.focus_chains is None and _chains(sses) == set("ABCDEFGH")
    lay, _, _ = make_layout(NUCLEOSOME, dna=False, focus="protomer")
    assert lay.focus_chains is not None and len(lay.focus_chains) == 4


def test_focus_auto_does_not_split_a_nucleic_acid_its_protomers_share():
    """A strand held by two or more protomers (NA_HOLD_MIN residues each) of a point-group assembly ties them into
    one complex; one with a strand per protomer, or none, does not. Filaments and cages keep focus."""
    from types import SimpleNamespace

    from foldmap.cli import _na_shared
    from foldmap.model import ResLabel

    bb = SimpleNamespace(labels=[ResLabel(c, k, "", "ALA") for k, c in enumerate("AAABBBCCC")])
    sym = SimpleNamespace(protomers=[["A"], ["B"], ["C"]])
    shared = SimpleNamespace(contacts={r: {(0, r)} for r in range(6)})
    each = SimpleNamespace(contacts={r: {(r // 3, r)} for r in range(9)})
    glancing = SimpleNamespace(contacts={0: {(0, 1)}, 1: {(0, 2)}, 2: {(0, 3)}, 3: {(0, 9)}})
    assert _na_shared(shared, sym, bb)  # no stems needed: a single strand two protomers hold is shared too
    assert not _na_shared(each, sym, bb) and not _na_shared(glancing, sym, bb)
    assert not _na_shared(None, sym, bb) and not _na_shared(shared, None, bb)


def test_4v4f_shared_rna_ring_draws_clear(monkeypatch):
    """F1's shared-nucleic-acid ring rule, end to end (K24): TRAP's 11-subunit ring (4V4F assembly 1, backbone only)
    with the RNA wound round it is one complex. focus auto draws the whole ring and every RNA piece, says nothing was
    left out, and routes every loop and strand clear (the full entry's two fallbacks, A2:102-103 and AZ:102-103, went
    with plan 10 Task 2's ring layout). Without the rule the same figure would be one subunit."""
    from foldmap import cli
    from foldmap.scope import drawn_scope
    from foldmap.spec import draw_spec, figure_spec

    ring = DATA / "r100" / "4V4F-ring.cif.gz"
    d = draw_spec(figure_spec(ring))
    assert d.layout.focus_chains is None and len(_chains(d.sses)) == 11
    assert len(d.layout.nucleic.strands) == 11 and len(d.layout.na_paths) == 11
    assert drawn_scope(d.model, d.layout, "auto", "auto")["reason"] == ""
    assert d.stuck == []
    monkeypatch.setattr(cli, "_na_shared", lambda *args: False)
    lay, sses, _ = make_layout(ring)
    assert lay.focus_chains is not None and len(lay.focus_chains) == 1


def test_shared_na_cage_stays_focused():
    """A cage keeps focus auto on one subunit even when a nucleic acid touches two or more of its protomers (F1
    review): only a ring or other point group (C, D) whose protomers share a nucleic acid is drawn whole (the 4V4F
    ring test). 7A4M's O cage with a strand held by four of its subunits."""
    from types import SimpleNamespace

    from foldmap.cli import _focus_chains, _na_shared
    from foldmap.io import load_backbone
    from foldmap.symmetry import detect_symmetry

    cage = DATA / "assemblies" / "7A4M.cif"
    _, sses, _ = make_layout(cage, focus="none")
    bb = load_backbone(cage)
    sym = detect_symmetry(bb, sses)
    assert sym.kind == "O"
    held = [k for c in sorted({l.chain for l in bb.labels})[:4] for k in range(len(bb)) if bb.labels[k].chain == c]
    shared = SimpleNamespace(contacts={r: {(0, 1)} for r in held}, strands=[SimpleNamespace(chain="X")])
    assert _na_shared(shared, sym, bb)  # the strand is held by four subunits
    focus = _focus_chains("auto", sym, bb, sses, shared)
    assert focus is not None and len(focus) == 1


def test_a_ring_of_helices_lends_a_subunit_its_two_neighbours():
    """C25 (random-100): a subunit drawn in focus borrowed every helix of a bundle one of its helices is in, so
    8QEM's ring of twelve C-terminal helices put eleven neighbours' copies on the page as primed ghosts in a row
    (15 lent elements; 6 now). It borrows the helices packed along its own (features.helix_partners): in a ring,
    the two beside it; in a three-helix stalk, both others."""
    import numpy as np

    from foldmap.cli import _focus_subset
    from helpers import fake_sse

    radius = 10.0 / (2 * np.sin(np.pi / 12))  # neighbours 10 Å apart round the ring
    ring = [
        fake_sse("H", 438, 451, (radius * np.cos(t), radius * np.sin(t), 0.0), (0, 0, 1), chain=c)
        for c, t in zip("ABCDEFGHIJKL", np.linspace(0, 2 * np.pi, 12, endpoint=False))
    ]
    kept, partial = _focus_subset(ring, {"A"}, [])
    assert sorted(s.chain for s in kept) == ["A", "B", "L"] and partial == {"B", "L"}
    stalk = [
        fake_sse("H", 438, 451, (6 * np.cos(t), 6 * np.sin(t), 0.0), (0, 0, 1), chain=c)
        for c, t in zip("ABC", (0, 2.1, 4.2))
    ]
    assert _focus_subset(stalk, {"A"}, [])[1] == {"B", "C"}
