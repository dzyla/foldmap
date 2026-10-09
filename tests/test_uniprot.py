"""UniProt annotation: entries, chain mapping (SIFTS or alignment), features on residues. Offline: fixtures."""

import shutil
from pathlib import Path

import pytest

from foldmap import uniprot
from foldmap.io import load_backbone

DATA = Path(__file__).parent / "data"
FIX = DATA / "uniprot"


@pytest.fixture
def cache(tmp_path, monkeypatch):
    """A cache holding the fixtures; any download attempt fails the test."""
    for f in FIX.iterdir():
        shutil.copy(f, tmp_path / f.name)
    monkeypatch.setattr(uniprot, "CACHE", tmp_path.parent)
    (tmp_path.parent / "uniprot").mkdir(exist_ok=True)
    for f in FIX.iterdir():
        shutil.copy(f, tmp_path.parent / "uniprot" / f.name)

    def offline(*a, **k):
        raise AssertionError("network used")

    monkeypatch.setattr(uniprot.urllib.request, "urlopen", offline)
    return tmp_path.parent / "uniprot"


def test_entry_parses_names_features_and_structures(cache):
    e = uniprot.entry("P07911")
    assert e.name == "Uromodulin" and e.gene == "UMOD" and e.organism == "Homo sapiens" and len(e.sequence) == 640
    assert ("Domain", "ZP") in {(f.type, f.description) for f in e.regions()}
    assert any(f.type == "Glycosylation" for f in e.sites())
    assert e.structures and e.alphafold == "P07911"
    assert e.structures[0].coverage >= e.structures[-1].coverage  # widest coverage first


def test_binding_sites_name_their_ligand(cache):
    sites = uniprot.entry("P69905").sites()
    heme = next(f for f in sites if f.type == "Binding site" and f.start == 88)
    assert heme.ligand == "heme b"


def test_bad_accession_is_explained(cache):
    with pytest.raises(ValueError, match="UniProt accession"):
        uniprot.entry("nonsense")


def test_sifts_maps_each_chain_of_a_pdb_entry(cache):
    maps = uniprot.chain_maps(DATA / "2HHB.cif")
    acc = {m.chain: m.accession for m in maps}
    assert acc == {"A": "P69905", "C": "P69905", "B": "P68871", "D": "P68871"}
    a = next(m for m in maps if m.chain == "A")
    assert a.to_label[2] == 1  # UniProt numbering counts the initiator Met, the structure does not


def test_features_land_on_the_right_residues(cache):
    bb = load_backbone(DATA / "2HHB.cif", "asu")
    ann = uniprot.annotate(DATA / "2HHB.cif", bb)
    heme = [(f, c, ks) for f, c, ks in ann.sites if f.type == "Binding site" and f.start == 88]
    assert {c for _, c, _ in heme} == {"A", "C"}
    for _, _c, ks in heme:  # UniProt His88 of alpha globin is the proximal His87 of the structure
        assert [(bb.labels[k].name, bb.labels[k].seq) for k in ks] == [("HIS", 87)]
    globin = [r for r in ann.regions if r[0].description == "Globin"]
    assert len(globin) == 4


def test_alignment_maps_a_chain_when_the_accession_is_given(cache):
    bb = load_backbone(DATA / "6ZS5.cif", "asu")
    maps = uniprot.chain_maps(DATA / "6ZS5.cif", accession="P07911")
    assert {m.chain for m in maps} == {"A", "D"}
    ann = uniprot.annotate(DATA / "6ZS5.cif", bb, accession="P07911")
    zp = [ks for f, c, ks in ann.regions if f.description == "ZP-N" and c == "A"]
    assert zp and bb.labels[zp[0][0]].seq >= 331


def test_describe_lists_structures_and_alphafold(cache):
    text = uniprot.describe(uniprot.entry("P07911"))
    assert "Uromodulin" in text and "ZP" in text and "AF-P07911-F1" in text and "Experimental structures" in text


def _gids(fig, prefix):
    return [a for a in fig.axes[0].patches if (a.get_gid() or "").startswith(prefix)]


def test_figure_marks_binding_sites_and_names_them_in_the_legend(cache):
    from foldmap.cli import make_figure

    fig = make_figure(DATA / "2HHB.cif", assembly="asu", uniprot="auto")
    assert len(_gids(fig, "feature:")) >= 4  # heme-binding His of each of four chains
    assert any("binding site (UniProt)" == t.get_text() for t in fig.axes[0].texts)
    plain = make_figure(DATA / "2HHB.cif", assembly="asu")
    assert not _gids(plain, "feature:")


def test_uniprot_domains_become_panels(cache):
    from foldmap.cli import make_layout

    lay, _, _ = make_layout(DATA / "2HHB.cif", assembly="asu", uniprot="auto", domains="uniprot")
    assert sorted(name for name, _ in lay.domains) == ["Globin · A", "Globin · B", "Globin · C", "Globin · D"]


def test_domains_uniprot_needs_the_annotation():
    from foldmap.cli import make_layout

    with pytest.raises(ValueError, match="--uniprot"):
        make_layout(DATA / "2HHB.cif", domains="uniprot")


def test_cli_uniprot_for_an_accession_and_a_structure(cache, capsys):
    from foldmap.cli import main

    assert main(["uniprot", "P07911"]) == 0
    assert "Uromodulin" in capsys.readouterr().out
    assert main(["uniprot", str(DATA / "2HHB.cif")]) == 0
    out = capsys.readouterr().out
    assert "chain A → P69905" in out and "chain B → P68871" in out and "Hemoglobin subunit beta" in out


def test_mcp_uniprot_tool(cache):
    import asyncio

    pytest.importorskip("mcp")
    from foldmap.mcp_server import server

    text = asyncio.run(server.call_tool("uniprot_info", {"source": "P69905"})).content[0].text
    assert "Hemoglobin subunit alpha" in text and "AF-P69905-F1" in text


def test_sequence_view_shows_the_uniprot_track(cache):
    from foldmap.seqplot import draw_sequence

    fig = draw_sequence(DATA / "6ZS5.cif", chain="A", uniprot="P07911", columns=60)
    regions = {
        a.get_gid().split(":")[1] for a in fig.axes[0].patches if (a.get_gid() or "").startswith("uniprot-region:")
    }
    assert len(regions) >= 2  # ZP and ZP-N at least, in separate lanes
    plain = draw_sequence(DATA / "6ZS5.cif", chain="A", columns=60)
    assert not any((a.get_gid() or "").startswith("uniprot-") for a in plain.axes[0].patches)


def test_interactive_page_lists_uniprot_features(cache):
    import json
    import re

    from foldmap.interactive import build_page

    page = build_page(DATA / "2HHB.cif", assembly="asu", uniprot="auto")
    data = json.loads(re.search(r'id="topo-data">(.*?)</script>', page, re.S).group(1))
    feats = data["links"]["features"]
    assert feats and all('id="feature:' + f["key"] + '"' in page for f in feats)
    assert any("heme" in f["text"] or "binding" in f["text"] for f in feats)
