from foldmap.features import MIN_CONTACT_PAIRS, helix_bundles, sse_contacts
from helpers import fake_sse, pipeline


def test_paired_strands_touch_and_distant_elements_do_not(ubq):
    bb, sses, sheets = pipeline(ubq)
    c = sse_contacts(bb, sses)
    ids = {s.id: s for s in sses}
    (sheet,) = sheets
    for a, b in zip(sheet.strands, sheet.strands[1:]):
        assert c.get(frozenset((a.id, b.id)), 0) >= MIN_CONTACT_PAIRS  # neighbours in a sheet
    assert all(v >= MIN_CONTACT_PAIRS for v in c.values())
    assert all(len(k) == 2 and all(i in ids for i in k) for k in c)


def test_contacts_span_chains(zs5):
    bb, sses, _ = pipeline(zs5)
    c = sse_contacts(bb, sses)
    chain = {s.id: s.chain for s in sses}
    assert any(len({chain[i] for i in k}) == 2 for k in c)  # the linker strands touch the other chain


def test_no_elements_no_contacts(ubq):
    bb, _, _ = pipeline(ubq)
    assert sse_contacts(bb, []) == {}


def _bundles(path, assembly="asu"):
    from foldmap.features import helix_bundles

    _, sses, _ = pipeline(path, assembly)
    return [set(g) for g in helix_bundles(sses)]


def test_8uup_stalk_and_heterodimer_pairs_are_bundles(ubq):
    groups = _bundles(ubq.parent / "8UUP.cif")
    stalk = {"B:456-472", "B:474-478", "B:480-485", "D:456-472", "D:474-477", "D:480-485"}
    assert stalk | {"F:456-472", "F:474-477", "F:480-483"} <= next(g for g in groups if "D:456-472" in g)
    for pair in ({"A:70-95", "B:195-217"}, {"C:70-95", "D:195-217"}, {"E:70-95", "F:195-217"}):
        assert any(pair <= g for g in groups), pair


def test_lambda_repressor_dimerisation_helices_are_a_bundle(ubq):
    assert _bundles(ubq.parent / "1LMB.cif") == [{"3:78-88", "4:78-88"}]


def test_no_bundles_without_a_long_parallel_interface(ubq):
    for name in ("1UBQ", "2HHB", "1TIM", "6ZS5", "2LZM"):
        assert _bundles(ubq.parent / f"{name}.cif") == [], name


def test_bacteriorhodopsin_trimer_interface_helices_are_not_a_bundle(ubq):
    # helix B of one protomer runs within 12 Å of helix D of the next over exactly 3/4 of its length
    assert _bundles(ubq.parent / "1C3W.cif", assembly="auto") == []


def test_a_short_kinked_piece_stays_with_its_helix_in_a_bundle():
    # chain A: a 16-residue helix, then after one residue a 4-residue piece whose fitted axis is meaningless
    a1 = fake_sse("H", 10, 25, centroid=(0, 0, 0), axis=(0, 0, 1), chain="A")
    a2 = fake_sse("H", 27, 30, centroid=(0, 0, 14), axis=(1, 0, 0), chain="A")
    b = fake_sse("H", 10, 35, centroid=(10, 0, 6), axis=(0, 0, 1), chain="B")  # packed alongside, 10 Å away
    assert helix_bundles([a1, a2, b]) == [["A:10-25", "A:27-30", "B:10-35"]]


def test_a_short_piece_does_not_join_antiparallel_helices():
    # chain A: a helix, a 4-residue turn piece, then a helix running back down beside the first (2LZM A:126-155);
    # chain B packs against the returning helix only, so only that helix is in the bundle
    a1 = fake_sse("H", 10, 25, centroid=(0, 0, 0), axis=(0, 0, 1), chain="A")
    a2 = fake_sse("H", 27, 30, centroid=(5, 0, 14), axis=(1, 0, 0.3), chain="A")
    a3 = fake_sse("H", 32, 47, centroid=(10, 0, 0), axis=(0, 0, -1), chain="A")
    b = fake_sse("H", 10, 25, centroid=(20, 0, 0), axis=(0, 0, 1), chain="B")
    assert helix_bundles([a1, a2, a3, b]) == [["A:32-47", "B:10-25"]]


def test_a_small_domain_folds_though_another_protein_has_one_domain_left():
    """features._fold_small_copies: a protein whose only domain is small cannot fold it, and that no longer stops
    the small domains of the other proteins from folding (it broke out of the loop)."""
    from foldmap.features import _fold_small_copies

    sses = [
        fake_sse("H", 1, 10, chain="A"),  # protein A: one small domain
        fake_sse("H", 1, 10, chain="B"),  # protein B: a small domain ...
        fake_sse("H", 20, 60, chain="B"),  # ... and a large one
        fake_sse("H", 70, 110, chain="B"),  # ... and another
    ]
    owner = {sses[0].id: ("A", 0), sses[1].id: ("B", 0), sses[2].id: ("B", 1), sses[3].id: ("B", 2)}
    _fold_small_copies(owner, {"A": "A", "B": "B"}, sses, {})
    assert sses[0].id not in owner  # A is not split
    assert owner[sses[1].id] == owner[sses[2].id] == ("B", 0) and owner[sses[3].id] == ("B", 1)
