"""Release metadata agrees with itself: one version, declared floors, matching citation files."""

import json
import re
import tomllib
from pathlib import Path

import foldmap

ROOT = Path(__file__).resolve().parents[1]


def _pyproject() -> dict:
    return tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def test_pyproject_takes_version_from_the_package():
    meta = _pyproject()
    assert "version" not in meta["project"] and "version" in meta["project"]["dynamic"]
    assert meta["tool"]["setuptools"]["dynamic"]["version"] == {"attr": "foldmap.__version__"}


def test_citation_version_matches_package():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    m = re.search(r"^version:\s*['\"]?([^'\"\s]+)", text, re.M)
    assert m and m.group(1) == foldmap.__version__, (
        f"CITATION.cff says {m and m.group(1)!r}; bump it (and date-released) to {foldmap.__version__!r}"
    )


def test_zenodo_version_matches_package():
    meta = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    assert meta.get("version") == foldmap.__version__, f'.zenodo.json needs "version": "{foldmap.__version__}"'


FLOORS = {
    "gemmi": "0.7",
    "numpy": "1.24",
    "scipy": "1.10",
    "networkx": "3.0",
    "matplotlib": "3.10",
    "PyYAML": "6.0.1",
}


def test_every_dependency_has_the_tested_floor():
    deps = dict(re.match(r"([A-Za-z0-9_.-]+)>=(\S+)", d).groups() for d in _pyproject()["project"]["dependencies"])
    assert deps == FLOORS, "keep pyproject floors and the minimum-deps CI job in step with FLOORS"
