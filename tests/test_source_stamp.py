"""Every campaign report says which code produced it (plan 28 Task 1, review 1.1; fix round 1, rulings Q10-Q13).

foldmap.source_stamp finds the checkout from the writing script's own path (never from where foldmap is imported, so a
non-editable install as CI's still names it), records the fingerprint of the foldmap that ran and of the checkout's
src/foldmap and whether they agree, git's commit and dirty state (or the commit a runner asserts, marked so), the
script, Python, the platform and the runtime dependencies. Every JSON or Markdown writer in scripts/ takes or uses
the stamp, found by call shape, except a short list of non-report writers."""

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
from foldmap.source_stamp import ASSERTED, RUNTIME, code_mix, source_stamp, stamp_line

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="needs git")


@pytest.fixture(autouse=True)
def _no_runner_claims(monkeypatch):
    monkeypatch.delenv("FOLDMAP_GIT_COMMIT", raising=False)
    monkeypatch.delenv("FOLDMAP_GIT_DIRTY", raising=False)


def test_the_stamp_names_the_script_from_its_checkout_however_foldmap_is_installed():
    """Review C1: CI installs foldmap non-editable; the script is still scripts/bench.py of this checkout."""
    s = source_stamp(SCRIPTS / "bench.py")
    assert s["foldmap"] == __version__
    assert s["code_fingerprint"] == app_jobs.code_fingerprint()  # the foldmap that ran, wherever it is
    assert s["script"]["name"] == "scripts/bench.py"
    here = app_jobs.code_fingerprint((ROOT / "src" / "foldmap").resolve())
    assert s["checkout_fingerprint"] == here
    assert s["installed_matches_checkout"] == (here == s["code_fingerprint"])
    editable = Path(app_jobs.__file__).resolve().parent == (ROOT / "src" / "foldmap").resolve()
    assert s["package"] == ("checkout" if editable else "installed")
    assert s["python"]["version"] == platform.python_version()
    assert s["platform"]["platform"] == platform.platform()
    assert set(s["packages"]) == set(RUNTIME) and s["packages"]["numpy"]
    if (ROOT / ".git").exists():
        assert s["git"]["source"] == "git" and len(s["git"]["commit"]) == 40 and isinstance(s["git"]["dirty"], bool)
    else:
        assert s["git"] is None


def test_the_runtime_list_is_pyprojects():
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    names = {d.split(">")[0].split("=")[0].split("<")[0].strip().lower() for d in deps}
    assert names == {n.lower() for n in RUNTIME}


def _checkout(root: Path, text: str = "X = 1\n") -> tuple[Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    (root / "pyproject.toml").write_text("[project]\nname = 'x'\n", encoding="utf-8")
    pkg = root / "src" / "foldmap"
    pkg.mkdir(parents=True)
    (pkg / "a.py").write_text(text, encoding="utf-8")
    script = root / "scripts" / "campaign.py"
    script.parent.mkdir()
    script.write_bytes(b"print(1)\r\n")
    return pkg, script


def test_without_git_the_fingerprint_names_the_code(tmp_path):
    pkg, script = _checkout(tmp_path)
    s = source_stamp(script, package=pkg)
    assert s["git"] is None and s["package"] == "checkout" and s["installed_matches_checkout"] is True
    assert s["code_fingerprint"] == app_jobs.code_fingerprint(pkg.resolve())
    assert "no git checkout" in stamp_line(s)


def test_a_runners_commit_is_recorded_as_asserted(tmp_path, monkeypatch):
    """Ruling Q11: a synced copy has no .git; the runner exports the commit it synced, marked as its claim."""
    pkg, script = _checkout(tmp_path)
    monkeypatch.setenv("FOLDMAP_GIT_COMMIT", "a" * 40)
    monkeypatch.setenv("FOLDMAP_GIT_DIRTY", "0")
    s = source_stamp(script, package=pkg)
    assert s["git"] == {"commit": "a" * 40, "dirty": False, "source": ASSERTED}
    assert "asserted by runner" in stamp_line(s)
    monkeypatch.setenv("FOLDMAP_GIT_DIRTY", "1")
    assert source_stamp(script, package=pkg)["git"]["dirty"] is True


def test_an_installed_copy_is_compared_with_the_checkout(tmp_path):
    """Ruling Q10: the package that ran may be installed elsewhere; the stamp says whether it is the checkout's code."""
    pkg, script = _checkout(tmp_path / "repo")
    site = tmp_path / "site" / "foldmap"
    shutil.copytree(pkg, site)
    s = source_stamp(script, package=site)
    assert s["package"] == "installed" and s["installed_matches_checkout"] is True
    assert s["script"]["name"] == "scripts/campaign.py"
    changed = tmp_path / "site2" / "foldmap"  # another path: code_fingerprint is taken once per path
    shutil.copytree(pkg, changed)
    (changed / "a.py").write_text("X = 2\n", encoding="utf-8")
    other = source_stamp(script, package=changed)
    assert other["installed_matches_checkout"] is False
    assert "differs from the checkout" in stamp_line(other)


@needs_git
def test_in_a_checkout_the_commit_and_whether_the_code_differs_from_it(tmp_path):
    pkg, script = _checkout(tmp_path)

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
    clean = source_stamp(script, package=pkg)["git"]
    assert clean == {"commit": head, "dirty": False, "dirty_paths": [], "source": "git"}
    (tmp_path / "notes.txt").write_text("outside the code\n", encoding="utf-8")
    assert source_stamp(script, package=pkg)["git"]["dirty"] is False  # only the code and tests/data count
    (tmp_path / "tests" / "data").mkdir(parents=True)
    (tmp_path / "tests" / "data" / "x.cif").write_text("data_x\n", encoding="utf-8")
    assert source_stamp(script, package=pkg)["git"]["dirty_paths"] == ["tests/data/"]
    (tmp_path / "tests" / "data" / "x.cif").unlink()
    (pkg / "a.py").write_text("X = 2\n", encoding="utf-8")
    dirty = source_stamp(script, package=pkg)["git"]
    assert dirty["commit"] == head and dirty["dirty"] is True and dirty["dirty_paths"] == ["src/foldmap/a.py"]
    assert f"git {head[:12]} (dirty: src/foldmap/a.py)" in stamp_line(source_stamp(script, package=pkg))


def test_a_failed_git_is_said(tmp_path):
    line = stamp_line({"git": {"commit": None, "source": "git", "error": "OSError: no git"}, "code_fingerprint": "f"})
    assert "git failed: OSError: no git" in line


def test_the_script_is_named_and_hashed_whatever_its_line_endings(tmp_path):
    pkg, script = _checkout(tmp_path)
    s = source_stamp(script, package=pkg)["script"]
    assert s == {"name": "scripts/campaign.py", "sha256": hashlib.sha256(b"print(1)\n").hexdigest()}
    script.write_bytes(b"print(1)\n")
    assert source_stamp(script, package=pkg)["script"] == s


def test_a_report_says_when_its_entries_were_measured_by_different_code(tmp_path):
    """Ruling Q13."""
    pkg, script = _checkout(tmp_path)
    one = source_stamp(script, package=pkg)
    assert code_mix([one, one]) == {"entry_code": {f"{one['code_fingerprint']}+{one['script']['sha256'][:16]}": 2},
                                    "mixed_code": False}  # fmt: skip
    other = one | {"code_fingerprint": "0" * 16}
    assert code_mix([one, other])["mixed_code"] is True
    assert code_mix([one, None])["mixed_code"] is True  # an unstamped entry beside a stamped one
    assert code_mix([None, None]) == {"entry_code": {"unstamped": 2}, "mixed_code": False}


def test_bench_baselines_record_the_stamp():
    spec = importlib.util.spec_from_file_location("bench_for_stamp", SCRIPTS / "bench.py")
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    meta = bench._meta()
    assert meta["source_stamp"]["code_fingerprint"] == app_jobs.code_fingerprint()
    assert meta["source_stamp"]["script"]["name"] == "scripts/bench.py"


# ---------------------------------------------------------------------------------------------- the writers

# References that count as taking or using the stamp in a writing function.
STAMP_CALLS = {"source_stamp", "_stamp", "run_stamp"}
STAMP_NAMES = {"SOURCE_STAMP", "_STAMP"}
REPORT_SUFFIXES = (".json", ".json.gz", ".md")
# Writers that write no measurement report, each with why (path from scripts/, function).
NOT_REPORTS = {
    ("bench.py", "_run"): "writes the baseline with _meta(), which carries SOURCE_STAMP",
    ("homologues.py", "cmd_diff"): "a pair's spec handed to its child process (an input); the child stamps its result",
    ("random100.py", "measure_view"): "a view's measures inside an entry's process; the entry record is stamped",
    ("homologues.py", "cmd_sample"): "a drawn pair set (the input of a run, not a measurement); its runs are stamped",
    ("make_na_fixtures.py", "main"): "committed test fixtures",
    ("random100.py", "cmd_sample"): "writes the sample sample() built, which carries the stamp",
    ("random100.py", "cmd_view"): "a view's measures inside an entry's process; the entry record is stamped",
    ("random100.py", "finish"): "re-judges stamped entries and rewrites them with their measuring stamps",
    ("random100.py", "cmd_entry"): "writes the entry run_entry() built, which carries the stamp",
    ("paper/extract.py", "main"): "provenance.json: the record of the sources' stamps",
    ("paper/make_numbers.py", "main"): "numbers derived from stamped sources (build/provenance.json)",
    ("paper/fill.py", "main"): "the filled manuscript and its key list",
    ("paper/crosscheck.py", "main"): "the prose checked against the sources",
    ("paper/table_s1.py", "main"): "Table S1, read from the code at the harness commit (named in the file)",
    ("paper/runtime.py", "main"): "a plot of the stamped reports' timings",
    ("paper/opm_fetch.py", "main"): "OPM's records as downloaded (an input)",
    ("paper/make_figures.py", "f5"): "the data F5 plots, from stamped sources",
    ("paper/make_figures.py", "main"): "figures.json and checks.json carry provenance() (the stamp) per figure",
}


def _uses_stamp(fn: ast.AST) -> bool:
    for node in ast.walk(fn):
        if isinstance(node, ast.Call):
            f = node.func
            name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else None
            if name in STAMP_CALLS:
                return True
        if isinstance(node, ast.Name) and node.id in STAMP_NAMES:
            return True
        if isinstance(node, ast.Constant) and node.value == "source_stamp":
            return True
    return False


def _report_path(node: ast.AST) -> bool:
    for n in ast.walk(node):
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value.endswith(REPORT_SUFFIXES):
            return True
    return False


def _dumps_json(node: ast.AST) -> bool:
    for n in ast.walk(node):
        if isinstance(n, ast.Call) and isinstance(n.func, (ast.Attribute, ast.Name)):
            name = n.func.attr if isinstance(n.func, ast.Attribute) else n.func.id
            if name in ("dumps", "dump"):
                return True
    return False


def _writes_a_report(call: ast.Call) -> bool:
    """json.dump(...), or write_text / write_bytes of JSON, or of a path ending .json / .json.gz / .md."""
    f = call.func
    if isinstance(f, ast.Attribute) and f.attr == "dump" and isinstance(f.value, ast.Name) and f.value.id == "json":
        return True
    if isinstance(f, ast.Attribute) and f.attr in ("write_text", "write_bytes"):
        return _report_path(f.value) or any(_dumps_json(a) for a in call.args)
    return False


def report_writers(source: str) -> list[tuple[str, bool]]:
    """(function, uses the stamp) for every function (or <module>) of a script that writes a report by call shape.
    A nested function counts as its outermost function, so a closure that writes what its function stamped passes."""
    tree = ast.parse(source)
    out: dict[str, bool] = {}
    top = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    for fn in top:
        if any(_writes_a_report(c) for c in ast.walk(fn) if isinstance(c, ast.Call)):
            out[fn.name] = _uses_stamp(fn)
    rest = ast.Module(body=[n for n in tree.body if n not in top], type_ignores=[])
    if any(_writes_a_report(c) for c in ast.walk(rest) if isinstance(c, ast.Call)):
        out["<module>"] = _uses_stamp(rest)
    return sorted(out.items())


def test_every_report_writer_in_scripts_takes_the_stamp():
    """Ruling Q12: a new campaign writing a report without the stamp fails here, whatever it names its file."""
    unstamped, listed_but_gone = [], set(NOT_REPORTS)
    for path in sorted(SCRIPTS.rglob("*.py")):
        rel = path.relative_to(SCRIPTS).as_posix()
        for fn, stamped in report_writers(path.read_text(encoding="utf-8")):
            listed_but_gone.discard((rel, fn))
            if not stamped and (rel, fn) not in NOT_REPORTS:
                unstamped.append(f"{rel}:{fn}")
    assert not unstamped, f"writers without foldmap.source_stamp (stamp them, or list why in NOT_REPORTS): {unstamped}"
    assert not listed_but_gone, f"NOT_REPORTS lists writers that no longer write: {sorted(listed_but_gone)}"


@pytest.mark.parametrize(
    "planted",
    [
        'import json\ndef run(out):\n    json.dump({"x": 1}, open(out, "w"))\n',
        'import json\ndef run(out):\n    p = out / "results.json"\n    p.write_text(json.dumps({"x": 1}))\n',
        'def run(out, s):\n    (out / f"{s}.md").write_text("# report\\n")\n',
        "import json\ndef run(out):\n    def inner():\n        out.write_text(json.dumps([1]))\n    inner()\n",
        'import json\nfrom pathlib import Path\nPath("x.json").write_text(json.dumps({}))\n',
    ],
)
def test_the_net_catches_a_planted_unstamped_writer(planted):
    found = report_writers(planted)
    assert found and not any(stamped for _, stamped in found)


def test_the_net_passes_a_stamped_writer():
    ok = (
        "import json\nfrom foldmap.source_stamp import source_stamp\nSOURCE_STAMP = source_stamp(__file__)\n"
        'def run(out):\n    doc = {"source_stamp": SOURCE_STAMP}\n'
        '    (out / "report.json").write_text(json.dumps(doc))\n'
    )
    assert report_writers(ok) == [("run", True)]
