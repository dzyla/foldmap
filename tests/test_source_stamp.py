"""Every campaign report says which code produced it (plan 28 Task 1, review 1.1): foldmap.source_stamp names the
commit and a dirty flag in a git checkout, the code fingerprint with or without git, the writing script, Python, the
platform and the runtime dependencies; and each campaign script's report writer stamps it as it writes."""

import ast
import hashlib
import importlib.util
import platform
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

from foldmap import __version__, app_jobs
from foldmap.source_stamp import RUNTIME, source_stamp, stamp_line

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="needs git")


def test_the_stamp_names_the_code_the_interpreter_and_the_dependencies():
    s = source_stamp()
    assert s["foldmap"] == __version__
    assert s["code_fingerprint"] == app_jobs.code_fingerprint()
    assert s["python"]["version"] == platform.python_version()
    assert s["platform"]["platform"] == platform.platform()
    assert set(s["packages"]) == set(RUNTIME) and s["packages"]["numpy"]
    if (ROOT / ".git").exists() and (ROOT / "src" / "foldmap").resolve() == Path(app_jobs.__file__).resolve().parent:
        assert len(s["git"]["commit"]) == 40 and isinstance(s["git"]["dirty"], bool)


def test_the_runtime_list_is_pyprojects():
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    names = {d.split(">")[0].split("=")[0].split("<")[0].strip().lower() for d in deps}
    assert names == {n.lower() for n in RUNTIME}


def _package(root: Path) -> Path:
    pkg = root / "src" / "foldmap"
    pkg.mkdir(parents=True)
    (pkg / "a.py").write_text("X = 1\n", encoding="utf-8")
    return pkg


def test_without_git_the_fingerprint_names_the_code(tmp_path):
    pkg = _package(tmp_path)
    s = source_stamp(package=pkg)
    assert s["git"] is None
    assert s["code_fingerprint"] == app_jobs.code_fingerprint(pkg.resolve())
    assert "no git checkout" in stamp_line(s)


@needs_git
def test_in_a_checkout_the_commit_and_whether_the_code_differs_from_it(tmp_path):
    pkg = _package(tmp_path)

    def git(*args):
        return subprocess.run(
            ["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t", *args],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()

    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "one")
    head = git("rev-parse", "HEAD")
    clean = source_stamp(package=pkg)["git"]
    assert clean == {"commit": head, "dirty": False, "dirty_paths": []}
    (tmp_path / "notes.txt").write_text("outside the code\n", encoding="utf-8")
    assert source_stamp(package=pkg)["git"]["dirty"] is False  # only src/, scripts/ and pyproject.toml count
    (pkg / "a.py").write_text("X = 2\n", encoding="utf-8")
    dirty = source_stamp(package=pkg)["git"]
    assert dirty["commit"] == head and dirty["dirty"] is True and dirty["dirty_paths"] == ["src/foldmap/a.py"]
    assert f"git {head[:12]} (dirty: src/foldmap/a.py)" in stamp_line(source_stamp(package=pkg))


def test_the_script_that_wrote_the_report_is_named_and_hashed(tmp_path):
    pkg = _package(tmp_path)
    script = tmp_path / "scripts" / "campaign.py"
    script.parent.mkdir()
    script.write_bytes(b"print(1)\r\n")
    s = source_stamp(script, package=pkg)["script"]
    assert s["name"] == "scripts/campaign.py"
    assert s["sha256"] == hashlib.sha256(b"print(1)\n").hexdigest()
    assert s == source_stamp(script, package=pkg)["script"]
    script.write_bytes(b"print(1)\n")  # line endings do not change the hash
    assert source_stamp(script, package=pkg)["script"] == s


# ---------------------------------------------------------------------------------------------- the writers

# Each campaign script's functions that write a report (or a record joined into one): each must call the stamp.
WRITERS = {
    "bench.py": ["_meta"],
    "random100.py": ["sample", "run_entry", "report"],
    "na_validate.py": ["main"],
    "homologues.py": ["write_usalign_report", "write_report", "run_pair", "diff_pair"],
    "sweep.py": ["main"],
    "paper/domains.py": ["main"],
    "paper/dssp_mkdssp.py": ["main"],
    "paper/membrane_opm.py": ["main"],
    "paper/make_figures.py": ["provenance"],
    "paper/fourth_sample.py": ["cmd_draw"],
}
STAMPS = {"source_stamp", "_stamp", "run_stamp"}
# files a campaign writes as its report: a function writing one of these must stamp it
REPORT_FILES = {"report.json", "results.json", "summary.json", "usalign.json", "domains.json", "run.json"}


def _functions(path: Path) -> dict[str, ast.FunctionDef]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}


def _stamps(fn: ast.AST) -> bool:
    for node in ast.walk(fn):
        if isinstance(node, ast.Call):
            f = node.func
            name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else None
            if name in STAMPS:
                return True
    return False


def _writes_a_report(fn: ast.AST) -> bool:
    for node in ast.walk(fn):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in ("write_text", "write_bytes"):
                names = {c.value for c in ast.walk(node.func.value) if isinstance(c, ast.Constant)}
                if names & REPORT_FILES:
                    return True
    return False


@pytest.mark.parametrize("script", sorted(WRITERS))
def test_each_campaign_writer_stamps_its_report(script):
    found = _functions(SCRIPTS / script)
    missing = [f for f in WRITERS[script] if f not in found]
    assert not missing, f"{script}: no function {missing}"
    unstamped = [f for f in WRITERS[script] if not _stamps(found[f])]
    assert not unstamped, f"{script}: {unstamped} write a report without foldmap.source_stamp"


def test_every_function_writing_a_campaign_report_is_listed_and_stamps():
    """A new campaign writing report.json (or another report file) without the stamp fails here."""
    unlisted = []
    for path in sorted(SCRIPTS.rglob("*.py")):
        rel = path.relative_to(SCRIPTS).as_posix()
        for name, fn in _functions(path).items():
            if _writes_a_report(fn) and (name not in WRITERS.get(rel, []) or not _stamps(fn)):
                unlisted.append(f"{rel}:{name}")
    assert not unlisted, f"report writers not stamping (or not in WRITERS): {unlisted}"


def test_bench_baselines_record_the_stamp():
    spec = importlib.util.spec_from_file_location("bench_for_stamp", SCRIPTS / "bench.py")
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    meta = bench._meta()
    assert meta["source_stamp"]["code_fingerprint"] == app_jobs.code_fingerprint()
    assert meta["source_stamp"]["script"]["name"] == "scripts/bench.py"
