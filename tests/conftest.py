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
