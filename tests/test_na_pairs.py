"""Base pairing against an independent annotator (RNApolis, stored in tests/data/na/reference.json)."""

import json
from pathlib import Path

NA = Path(__file__).parent / "data" / "na"
FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB"]
REFERENCE = json.loads((NA / "reference.json").read_text(encoding="utf-8"))


def test_every_fixture_has_reference_pairs():
    assert sorted(REFERENCE) == sorted(FIXTURES)
    for pdb in FIXTURES:
        assert (NA / f"{pdb}.cif.gz").is_file()
        assert REFERENCE[pdb], pdb
