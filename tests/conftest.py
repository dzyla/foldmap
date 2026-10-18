from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def ubq():
    return DATA / "1UBQ.cif"


@pytest.fixture(scope="session")
def zs5():
    return DATA / "6ZS5.cif"


@pytest.fixture(scope="session")
def zya():
    return DATA / "6ZYA.cif"


@pytest.fixture(autouse=True)
def _temp_cache(tmp_path, monkeypatch):
    """No test reads, writes or prunes the real ~/.cache/foldmap (uploads, downloads, UniProt entries, CATH
    choppings)."""
    from foldmap import cath, fetch, uniprot

    monkeypatch.setattr(fetch, "CACHE", tmp_path / "cache")
    monkeypatch.setattr(uniprot, "CACHE", tmp_path / "cache")
    monkeypatch.setattr(cath, "CACHE", tmp_path / "cache")


def pytest_addoption(parser):
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        help="rewrite tests/data/snapshots/*.json from the current figures (review the diff before committing)",
    )


@pytest.fixture(scope="session")
def update_snapshots(request):
    return request.config.getoption("--update-snapshots")
