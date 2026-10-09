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
