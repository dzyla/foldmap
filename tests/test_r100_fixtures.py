"""Campaign fixtures stay small (plan 10 budget, 1 MB since plan 12) and out of the bench."""

from pathlib import Path

R100 = Path(__file__).parent / "data" / "r100"


def test_fixtures_are_small():
    sizes = {p.name: p.stat().st_size for p in R100.glob("*.cif.gz")}
    assert sizes and all(s < 100_000 for s in sizes.values()), sizes
    assert sum(sizes.values()) < 1_000_000, sum(sizes.values())  # raised from 900 KB for 1SLV (plan 12)


def test_fixtures_stay_out_of_the_bench():
    assert not list(R100.glob("*.cif"))
