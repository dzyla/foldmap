"""Campaign fixtures stay small (plan 10 budget) and out of the bench."""

from pathlib import Path

R100 = Path(__file__).parent / "data" / "r100"


def test_fixtures_are_small():
    sizes = {p.name: p.stat().st_size for p in R100.glob("*.cif.gz")}
    assert sizes and all(s < 100_000 for s in sizes.values()), sizes
    assert sum(sizes.values()) < 800_000, sum(sizes.values())


def test_fixtures_stay_out_of_the_bench():
    assert not list(R100.glob("*.cif"))
