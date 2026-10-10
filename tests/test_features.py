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


def _bundles(path):
    from foldmap.features import helix_bundles

    _, sses, _ = pipeline(path)
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


def test_a_short_kinked_piece_stays_with_its_helix_in_a_bundle():
    # chain A: a 16-residue helix, then after one residue a 4-residue piece whose fitted axis is meaningless
    a1 = fake_sse("H", 10, 25, centroid=(0, 0, 0), axis=(0, 0, 1), chain="A")
    a2 = fake_sse("H", 27, 30, centroid=(0, 0, 14), axis=(1, 0, 0), chain="A")
    b = fake_sse("H", 10, 35, centroid=(10, 0, 6), axis=(0, 0, 1), chain="B")  # packed alongside, 10 Å away
    assert helix_bundles([a1, a2, b]) == [["A:10-25", "A:27-30", "B:10-35"]]
