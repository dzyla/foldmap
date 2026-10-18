"""Release metadata agrees with itself: one version, declared floors, matching citation files."""

import json
import re
import shlex
import tomllib
from pathlib import Path

import pytest
import yaml

import foldmap

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "tests.yml"
needs_workflow = pytest.mark.skipif(not WORKFLOW.exists(), reason="the CI workflow is not part of the sdist")


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
    "networkx": "3.1",  # science.pairing_cycle: nx.simple_cycles on an undirected graph (new in 3.1)
    "matplotlib": "3.10",
    "PyYAML": "6.0.1",
}


def test_every_dependency_has_the_tested_floor():
    deps = dict(re.match(r"([A-Za-z0-9_.-]+)>=(\S+)", d).groups() for d in _pyproject()["project"]["dependencies"])
    assert deps == FLOORS, "keep pyproject floors and the minimum-deps CI job in step with FLOORS"


TESTED = {"3.11", "3.12", "3.13", "3.14"}


def test_classifiers_name_the_tested_pythons():
    meta = _pyproject()["project"]
    listed = {
        c.rsplit(" :: ", 1)[1] for c in meta["classifiers"] if re.match(r"Programming Language :: Python :: 3\.\d+$", c)
    }
    assert listed == TESTED and meta["requires-python"] == ">=3.11"


@needs_workflow
def test_ci_matrix_covers_the_tested_pythons():
    tests_yml = WORKFLOW.read_text(encoding="utf-8")
    assert all(f'"{v}"' in tests_yml for v in TESTED)


@needs_workflow
def test_linux_jobs_that_run_pytest_install_the_snapshot_font():
    """The snapshot gate needs Liberation Sans; without it the fallback font moves text by several points."""
    jobs = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))["jobs"]
    checked = 0
    for name, job in jobs.items():
        steps = job["steps"]
        if not any("pytest" in s.get("run", "") for s in steps):
            continue
        if not str(job["runs-on"]).startswith("ubuntu") and "matrix.os" not in str(job["runs-on"]):
            continue
        installs = [s for s in steps if "fonts-liberation" in s.get("run", "")]
        installs += [s for x in steps if x.get("uses") == FONT_ACTION for s in _font_action_steps()]
        assert installs, f"job {name} runs pytest on Linux without installing fonts-liberation"
        if "matrix.os" in str(job["runs-on"]):
            assert all(s.get("if") == "runner.os == 'Linux'" for s in installs), f"job {name}: guard the install"
        checked += 1
    assert checked >= 3


FONT_ACTION = "./.github/actions/liberation-sans"  # Liberation Sans 2.1.5 on every OS (tests.yml, determinism.yml)


def _font_action_steps() -> list[dict]:
    """The composite action's steps that install fonts-liberation (apt: Linux only)."""
    action = yaml.safe_load((ROOT / FONT_ACTION / "action.yml").read_text(encoding="utf-8"))
    return [s for s in action["runs"]["steps"] if "fonts-liberation" in s.get("run", "")]


@needs_workflow
def test_minimum_deps_job_pins_equal_the_floors():
    """The floor "1.24" is pinned as "==1.24.*" and "0.7" as "==0.7.0"; nothing else is accepted.

    The minimum-deps pins must stay on a single `pip install` line; this test reads that one line.
    """
    workflow = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    runs = [s["run"] for s in workflow["jobs"]["minimum-deps"]["steps"] if "run" in s]
    (install,) = [r for r in runs if r.startswith("pip install ") and "==" in r]
    pins = dict(re.fullmatch(r"([A-Za-z0-9_.-]+)==(\S+)", t).groups() for t in shlex.split(install) if "==" in t)
    assert set(pins) == set(FLOORS)
    for name, floor in FLOORS.items():
        assert pins[name] in (f"{floor}.*", f"{floor}.0", floor), (
            f"minimum-deps pins {name}=={pins[name]} but the floor is {floor}"
        )


HEADER = (
    "# SPDX-License-Identifier: AGPL-3.0-or-later\n"
    "# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.\n"
)


def test_notice_names_the_attribution_the_code_shows():
    from foldmap import ATTRIBUTION, CREDIT

    notice = (ROOT / "NOTICE").read_text(encoding="utf-8")
    assert ATTRIBUTION == "foldmap, by Dawid Zyla"
    assert f'"{ATTRIBUTION}"' in notice and ATTRIBUTION in CREDIT
    assert "AGPL-3.0" in CREDIT and "cite" in CREDIT


def test_output_exception_is_in_notice_and_the_credit_points_to_it():
    from foldmap import CREDIT

    notice = (ROOT / "NOTICE").read_text(encoding="utf-8")
    assert "Output exception" in notice
    assert "output exception: see NOTICE" in CREDIT


def test_citation_file_ships_in_the_wheel_metadata():
    assert "CITATION.cff" in _pyproject()["project"]["license-files"]


def test_every_source_file_names_its_license():
    files = sorted((ROOT / "src" / "foldmap").glob("*.py")) + sorted((ROOT / "scripts").glob("*.py"))
    missing = [f.name for f in files if not f.read_text(encoding="utf-8").startswith(HEADER)]
    assert not missing, f"add the SPDX header to: {', '.join(missing)}"
    js = (ROOT / "site" / "app.js").read_text(encoding="utf-8")
    assert js.startswith("// SPDX-License-Identifier: AGPL-3.0-or-later\n")
