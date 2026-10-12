"""Protein fixtures for the published-figure comparisons (plan 11) and the loop-crossing counter."""

import gemmi

from helpers import FID, loop_crossings

IDS = ["1TEN", "1LDM", "2ACY", "1MBN", "7TIM", "1PV7", "3CHY"]


def test_fixtures_are_small_and_complete():
    sizes = {i: (FID / f"{i}.cif.gz").stat().st_size for i in IDS}
    assert all(s < 100_000 for s in sizes.values()), sizes
    assert sum(sizes.values()) < 250_000
    for i in IDS:
        st = gemmi.read_structure(str(FID / f"{i}.cif.gz"))
        assert st.assemblies, i  # assembly records kept (1LDM's tetramer is built from them)
        assert not any(r.is_water() for c in st[0] for r in c), i


def test_1ldm_builds_its_tetramer():
    from foldmap.io import load_backbone

    assert len({l.chain for l in load_backbone(FID / "1LDM.cif.gz", "auto").labels}) == 4


def test_fixtures_stay_out_of_the_bench():
    assert not list(FID.glob("*.cif"))  # scripts/bench.py globs tests/data/*.cif only


def test_loop_crossings_counts_proper_crossings_only():
    class L:  # duck-typed route.Loop
        def __init__(self, pts):
            self.points = pts

    plus = [L([(0, -1), (0, 1)]), L([(-1, 0), (1, 0)])]
    touch = [L([(0, 0), (0, 1)]), L([(0, 0), (1, 0)])]  # sharing an end point is not a crossing
    assert loop_crossings(plus) == 1 and loop_crossings(touch) == 0
    zigzag = [L([(0, -1), (0, 1), (2, 1), (2, -1)]), L([(-1, 0), (3, 0)])]  # one segment pair per leg it cuts
    assert loop_crossings(zigzag) == 2
    collinear = [L([(0, 0), (2, 0)]), L([(1, 0), (3, 0)])]  # overlapping on one line is not a proper crossing
    assert loop_crossings(collinear) == 0
