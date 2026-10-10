"""The release guard: a tag must name the version the package carries."""

import importlib.util
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / ".github" / "workflows" / "release.yml"
needs_workflows = pytest.mark.skipif(not RELEASE.exists(), reason="the workflows are not part of the sdist")


def _check():
    spec = importlib.util.spec_from_file_location("check_release", ROOT / "scripts" / "check_release.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.check


def test_tag_with_or_without_v_matches():
    check = _check()
    assert check("v1.2.3", "1.2.3") is None
    assert check("1.2.3", "1.2.3") is None
    assert check("refs/tags/v1.2.3", "1.2.3") is None


def test_mismatched_tag_is_refused():
    problem = _check()("v1.2.4", "1.2.3")
    assert problem is not None and "1.2.4" in problem and "1.2.3" in problem


def test_prerelease_must_match_exactly():
    assert _check()("v1.2.3rc1", "1.2.3") is not None


@needs_workflows
def test_release_runs_the_tests_and_checks_the_tag_before_building():
    jobs = yaml.safe_load(RELEASE.read_text(encoding="utf-8"))["jobs"]
    assert jobs["tests"]["uses"] == "./.github/workflows/tests.yml"
    assert jobs["build"]["needs"] == "tests" and jobs["publish"]["needs"] == "build"
    assert any("check_release.py" in s.get("run", "") for s in jobs["build"]["steps"])
