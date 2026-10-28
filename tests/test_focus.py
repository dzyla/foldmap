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
    auto lays out both halves (and says it drew an assembly of two whole: spec.draw_spec judges that figure's floor);
    asked for, the protomer is still one half."""
    lay, sses, _ = make_layout(NUCLEOSOME, dna=False)
    assert lay.focus_chains is None and _chains(sses) == set("ABCDEFGH") and lay.whole_of_two
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


PSII = DATA / "app29" / "31FD-O-small.cif.gz"  # a crop of 31FD (plan 29's tune half), the C2 photosystem II dimer


@pytest.fixture
def page_held_at_16_in(monkeypatch):
    """Pages held at 16 in on their longer side (layout.MAX_PAGE_IN, 48 in): the crop of 31FD drawn whole (22-23 in at
    6 pt) is below the readability floor, as the whole entry is at 48 in, and one subunit (11-12 in) is above it."""
    from foldmap import layout

    monkeypatch.setattr(layout, "MAX_PAGE_IN", 16.0)


def test_a_two_subunit_complex_below_the_floor_draws_one_subunit(page_held_at_16_in):
    """Plan 29 Task 5 (ruling P29-7, 31FD): focus auto drew this C2 photosystem II dimer whole although that figure is
    below the readability floor (the app's warning offered only chain A). Now the whole is judged and, below the floor,
    drawn again as one subunit, said as any other subunit focus auto draws. The fixture is 31FD's two copies of PsbO
    (O, o) and the ten single-helix subunits (E F H I K L M T U X and their copies), backbone (N, CA, C, O) only,
    rewritten by gemmi (B 50, occupancy 1): the whole entry so trimmed is 298 kB, so the page is held small instead
    (the whole entry, at 48 in, is test_31fd_is_drawn_as_one_subunit)."""
    from foldmap.scope import drawn_scope
    from foldmap.spec import draw_spec, figure_spec

    whole = draw_spec(figure_spec(PSII, focus="none"))
    assert whole.floor and whole.floor[0]["reason"] == "text_small"  # the fault: the whole is below the floor
    d = draw_spec(figure_spec(PSII))
    assert d.floor == []
    assert d.layout.subunits == (2, "C2")
    assert len({c.upper() for c in d.layout.focus_chains}) == len(d.layout.focus_chains) == 11  # one of each chain
    assert d.layout.context.startswith("1 of 2 subunits shown · C2 symmetry")
    scope = drawn_scope(d.model, d.layout, "auto", "auto")
    assert scope["reason"].startswith("one of 2 subunits (C2 symmetry)") and scope["all_options"] == {"focus": "none"}
    lay, _, _ = make_layout(PSII)  # laid out alone, nothing judged: whole, and said to be an assembly of two
    assert lay.focus_chains is None and lay.whole_of_two


def test_a_two_subunit_complex_above_the_floor_stays_whole():
    """The control (Task 3 review): the C2 histone octamer without its DNA (3UTB), above the floor drawn whole, is
    drawn whole, both halves and every chain; drawn once (no second figure)."""
    from foldmap import spec as spec_mod

    drawn = []
    once = spec_mod._draw_once

    def counted(*args, **kw):
        drawn.append(kw.get("below_floor", False))
        return once(*args, **kw)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(spec_mod, "_draw_once", counted)
        d = spec_mod.draw_spec(spec_mod.figure_spec(NUCLEOSOME, dna=False))
    assert d.floor == [] and d.layout.whole_of_two and d.layout.focus_chains is None
    assert _chains(d.sses) == set("ABCDEFGH") and drawn == [False]


def test_a_shared_dna_keeps_the_nucleosome_whole_below_the_floor(monkeypatch):
    """The nucleosome's DNA ties its two halves (_na_shared): below the floor too, focus auto draws it whole (no
    subunit of two is taken where its protomers share a nucleic acid)."""
    from foldmap import layout
    from foldmap.spec import draw_spec, figure_spec

    monkeypatch.setattr(layout, "MAX_PAGE_IN", 6.0)
    d = draw_spec(figure_spec(NUCLEOSOME))
    assert d.floor and not d.layout.whole_of_two and d.layout.focus_chains is None
    assert _chains(d.sses) == set("ABCDEFGH")


def test_a_reference_to_the_other_subunit_keeps_the_whole(page_held_at_16_in):
    """An element the figure names on a chain one subunit leaves out (a rename) cannot be drawn by that subunit
    (BadReference): the whole figure is kept, below the floor, as before."""
    from foldmap.spec import draw_spec, figure_spec

    kept = draw_spec(figure_spec(PSII)).layout.focus_chains
    lay, sses, _ = make_layout(PSII)
    other = next(s for s in sses if s.chain not in kept)
    d = draw_spec(figure_spec(PSII, rename={f"res:{other.chain}:{other.first.seq}-{other.last.seq}": "X"}))
    assert d.floor and d.layout.focus_chains is None and other.chain in _chains(d.sses)


def test_the_app_draws_one_subunit_below_the_floor(page_held_at_16_in):
    """The app's path (app_jobs.draw_settings, the app's default settings) reaches the same decision, with no floor
    note and the scope said in focus auto's words."""
    pytest.importorskip("streamlit")  # the app extra: absent in the sdist and minimum-deps jobs
    from foldmap import app, app_jobs
    from foldmap.floor import said
    from foldmap.scope import drawn_scope

    settings = {**app.LAYOUT_DEFAULTS, "theme": "publication"}
    d, _ = app_jobs.draw_settings(PSII, settings)
    assert d.layout.subunits == (2, "C2") and said(d, "app") == ([], "")
    assert drawn_scope(d.model, d.layout, "auto", "auto")["reason"].startswith("one of 2 subunits")


def test_31fd_is_drawn_as_one_subunit():
    """The same, end to end on the whole 31FD (the download; skipped where it is not): drawn whole, its text was
    5.0 pt on a 48 x 9 in page; one subunit is above the floor."""
    from foldmap.spec import draw_spec, figure_spec
    from helpers import r100_entry

    d = draw_spec(figure_spec(r100_entry("31FD")))
    assert d.floor == [] and d.layout.subunits == (2, "C2") and len(d.layout.focus_chains) == 15


def _page_chains(html: str) -> set[str]:
    """The chains of the elements an interactive page holds (its topo-data JSON)."""
    import json
    import re

    blob = re.search(r'<script type="application/json" id="topo-data">(.*?)</script>', html, re.S).group(1)
    return {e["chain"] for e in json.loads(blob)["elements"]}


def test_every_route_takes_the_figures_subunit(page_held_at_16_in, tmp_path):
    """Plan 29 Task 5 fix round 1 (ruling P29-12): the data, the interactive page and the layout file are of the figure
    draw_spec drew, one subunit of the 31FD crop below the floor, by the decision carried from the figure
    (spec.below_floor_of): given (the figure, a callable, a bool) or drawn to decide when none is at hand."""
    from foldmap import layoutfile
    from foldmap.datafile import spec_data
    from foldmap.interactive import build_page
    from foldmap.spec import draw_spec, figure_spec, layout_document

    spec = figure_spec(PSII)
    d = draw_spec(spec)
    assert d.took_subunit and d.layout.subunits == (2, "C2")
    drawn_chains = {s.chain for s in d.sses}
    for given in (None, d, lambda: d, True):
        assert spec_data(spec, drawn=given)["structure"]["focus_chains"] == sorted(d.layout.focus_chains)
    assert spec_data(spec, drawn=False)["structure"]["focus_chains"] is None  # the decision is taken as given
    assert _page_chains(build_page(PSII)) == _page_chains(build_page(PSII, drawn=d)) == drawn_chains
    label = next(p.label for p in d.layout.placed.values() if p.sse.chain in d.layout.focus_chains)
    renamed = figure_spec(PSII, rename={label: "X"})  # a label both copies share, saved as the subunit's residues
    figure = draw_spec(renamed)
    for given in (True, None, figure):
        out = layoutfile.save(tmp_path / "psii.yaml", layout_document(renamed, drawn=given))
        back = draw_spec(figure_spec(PSII, layout_file=out))
        assert back.layout.focus_chains == figure.layout.focus_chains and back.floor == []
    whole = layoutfile.save(tmp_path / "whole.yaml", layout_document(renamed, drawn=False))
    assert draw_spec(figure_spec(PSII, layout_file=whole)).layout.focus_chains is None  # both copies named: kept whole


def test_a_figure_drawn_twice_says_its_warnings_once(page_held_at_16_in, tmp_path, capsys):
    """M2 (fix round 1): what the whole figure's drawing writes to standard error (a data row not used, what was read)
    is held and dropped when the figure is drawn again as one subunit: said once, of the figure drawn."""
    from foldmap.spec import draw_spec, figure_spec

    table = tmp_path / "values.csv"
    table.write_text("chain,residue,value\nO,120,1.5\nO,9999,2.0\n", encoding="utf-8")
    d = draw_spec(figure_spec(PSII, data=str(table)))
    err = capsys.readouterr().err
    assert d.took_subunit and err.count("9999") == 1 and err.count("foldmap: read ") <= 1, err


def test_the_apps_downloads_are_the_figures(page_held_at_16_in, tmp_path):
    """The app's Data JSON, its 3D view's page and its layout file (with a rename) are of the figure it shows,
    one subunit of the 31FD crop: the data and page jobs take the decision from the kept figure, the layout file
    from the draw job's took_subunit, drawing nothing."""
    pytest.importorskip("streamlit")  # the app extra: absent in the sdist and minimum-deps jobs
    import json

    import yaml

    from foldmap import app, app_jobs
    from foldmap.app_views import layout_yaml
    from foldmap.spec import draw_spec, figure_spec

    settings = {**app.LAYOUT_DEFAULTS, "theme": "publication"}
    blob = json.dumps(settings, sort_keys=True)
    got = app_jobs.draw_job(str(PSII), blob)
    d, _ = app_jobs.draw_settings(PSII, settings)
    assert got["took_subunit"] and d.took_subunit
    data = json.loads(app_jobs.data_job(str(PSII), blob, "json"))
    assert data["structure"]["focus_chains"] == sorted(d.layout.focus_chains)
    assert _page_chains(app_jobs.page_job(str(PSII), blob)) == {s.chain for s in d.sses}
    label = next(p.label for p in d.layout.placed.values() if p.sse.chain in d.layout.focus_chains)
    renamed = json.dumps({**settings, "rename": f"{label}=X"}, sort_keys=True)  # a label both copies share
    figure = draw_spec(app.spec_of(PSII, json.loads(renamed)))
    text = layout_yaml(str(PSII), renamed, None, app_jobs.draw_job(str(PSII), renamed)["took_subunit"])
    saved = tmp_path / "psii.yaml"
    saved.write_bytes(text)
    saved_as = set(yaml.safe_load(text.decode())["edits"]["rename"])  # res:CHAIN:a-b, the subunit's chains only
    assert saved_as and {r.split(":")[1] for r in saved_as} <= figure.layout.focus_chains and figure.took_subunit
    back = draw_spec(figure_spec(PSII, layout_file=saved))
    assert back.layout.focus_chains == figure.layout.focus_chains


def _captured(name: str):
    """A captured input of cli._one_of_each (tests/data/app29/NAME.json.gz): each chain's residue numbers, names and CA
    as io.load_backbone read them from the whole entry. The other backbone atoms are not kept (stood in by the CA):
    _one_of_each, _focus_plan and _context read no other. Returns the backbone and the document."""
    import gzip
    import json

    import numpy as np

    from foldmap.model import Backbone, ResLabel

    with gzip.open(DATA / "app29" / f"{name}.json.gz", "rt", encoding="utf-8") as fh:
        doc = json.load(fh)
    labels, ca = [], []
    for c, d in doc["chains"].items():
        icodes = d["icode"] or [""] * len(d["seq"])
        labels += [ResLabel(c, s, i, n) for s, i, n in zip(d["seq"], icodes, d["name"].split())]
        ca += d["ca"]
    ca = np.asarray(ca, float).reshape(-1, 3)
    prev = np.array([k > 0 and labels[k - 1].chain == l.chain for k, l in enumerate(labels)])
    return Backbone(labels, np.repeat(ca[:, None, :], 4, axis=1), prev), doc


def test_a_copy_modelled_in_part_is_no_protein_of_its_own():
    """Plan 29 Task 8b (ruling P29-14; docs/validation/app29/diagnosis-t8.md): 26EL is twenty copies of one protein
    (residues 1-682), C2. Its protomer holds F, modelled 1-291, and M, 292-682: against A's whole sequence their
    ratios are 0.59 and 0.73, so _one_of_each took them as two more proteins and focus auto drew A, F and M, 132 to
    301 Å apart, said as "1 of 7". model.homolog_chains, the rule the labels share by, calls them A's copies."""
    from types import SimpleNamespace

    from foldmap.cli import _context, _focus_plan

    bb, doc = _captured("26EL-one-of-each")  # 76 kB; and the C2 symmetry detect_symmetry found: two protomers of ten
    sym = SimpleNamespace(**doc["sym"])
    focus, whole_of_two = _focus_plan("auto", sym, bb, [None] * doc["n_sses"])
    assert focus == {"A"} and not whole_of_two
    assert _context(sym, focus, bb, neighbours=False) == "1 of 20 subunits shown · C2 symmetry"


def test_26el_is_drawn_as_one_chain():
    """The same, end to end on the whole 26EL (the download; skipped where it is not)."""
    from foldmap.spec import draw_spec, figure_spec
    from helpers import r100_entry

    d = draw_spec(figure_spec(r100_entry("26EL")))
    assert d.layout.focus_chains == {"A"} and d.layout.context.startswith("1 of 20 subunits shown · C2 symmetry")


def test_one_protein_by_its_numbers_is_today_one_protein_even_as_two_entities():
    """The control the fix's risk asks for (ruling P29-14), pinned as it is: two chains that share their residue
    numbers and are alike at them are one protein to model.homolog_chains, and so now to _one_of_each, even when the
    file holds them as two entities. No near-identical paralogs sharing numbering exist among the 710 structures at
    hand (tests/data, plan 29's tune half and the download cache less its check half): every pair of chains of two
    entities alike at 70 % or more of 20 or more shared numbers is 97 % alike or more, and homolog_chains joins each.
    Three of them it joins where the whole-sequence ratio did not: HIV-1 reverse transcriptase's p66 and p51 (8VBG and
    1HYS; one gene, p51 its first 440 residues, folded otherwise) and 6V6S's unassigned poly-UNK chains of different
    entities. None is symmetric, so focus auto never reaches _one_of_each on them, and no figure changes. The risk,
    named: a symmetric assembly of p66/p51 heterodimers (or of poly-UNK models) would draw p66 alone as its subunit.
    8VBG's p66 (A, 1-542) and p51 (B, 6-428), captured (app29/8VBG-p66-p51.json.gz, 8 kB)."""
    from difflib import SequenceMatcher

    from foldmap.cli import _one_of_each
    from foldmap.model import homolog_chains

    bb, _ = _captured("8VBG-p66-p51")
    seq = {c: [l.name for l in bb.labels if l.chain == c] for c in "AB"}
    assert SequenceMatcher(None, seq["A"], seq["B"], autojunk=False).ratio() < 0.9  # two proteins, by the ratio
    assert homolog_chains(bb)["B"] == ("A", 0)  # one protein, by the numbers
    assert _one_of_each(["A", "B"], bb) == {"A"}


def test_near_identical_paralogs_below_the_copy_rule_stay_two_proteins():
    """The rule's own edge, on a synthetic pair: two chains numbered alike, 1-100, differing at 15 of the 100 (85 %
    alike: under homolog_chains' 90 % and under the ratio's 0.9) stay two proteins; at 5 they are one."""
    import numpy as np

    from foldmap.cli import _one_of_each
    from foldmap.model import Backbone, ResLabel

    names = ["ALA", "GLY", "LEU", "SER", "VAL", "THR", "ASP", "GLU", "LYS", "ARG"]

    def pair(differ: int) -> Backbone:
        a = [names[k % 10] for k in range(100)]
        b = [("TRP" if k % 7 == 0 and k // 7 < differ else n) for k, n in enumerate(a)]
        labels = [ResLabel("A", k + 1, "", n) for k, n in enumerate(a)] + [
            ResLabel("B", k + 1, "", n) for k, n in enumerate(b)
        ]
        ca = np.array([(k * 3.8, 0.0, 0.0) for k in range(100)] + [(k * 3.8, 10.0, 0.0) for k in range(100)])
        prev = np.array([k % 100 > 0 for k in range(200)])
        return Backbone(labels, np.repeat(ca[:, None, :], 4, axis=1), prev)

    assert _one_of_each(["A", "B"], pair(15)) == {"A", "B"}
    assert _one_of_each(["A", "B"], pair(5)) == {"A"}
