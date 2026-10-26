"""The sensitivity of the tuned constants (plan 28 Task 9.2, ruling Q29): scripts/perturb.py reaches every reader of
a perturbed constant (a planted package, then foldmap in random100.py's own process), and nothing changes without
it; scripts/sensitivity.py's perturbations, stage A subset, campaign keys and classes."""

import importlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
VALIDATION = ROOT / "docs" / "validation"
needs_records = pytest.mark.skipif(
    not (VALIDATION / "inventory" / "inventory.md").exists()
    or not (VALIDATION / "random100-third" / "sample.json").exists(),
    reason="the validation records are not part of the sdist",
)
needs_report = pytest.mark.skipif(
    not (VALIDATION / "sensitivity" / "report.json").exists(), reason="the sensitivity report is not part of the sdist"
)


def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"foldmap_{name}_t", SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def perturb():
    return _script("perturb")


@pytest.fixture(scope="module")
def sens():
    return _script("sensitivity")


# ------------------------------------------------------------------------------------------------ the patch

PLANTED = {
    "__init__.py": "",
    "a.py": "LIMIT = 3.0  # a planted constant\nDOUBLE = LIMIT * 2\n\n\ndef f(x=LIMIT):\n    return x\n\n\n"
    "def g():\n    return LIMIT\n",
    "b.py": "from .a import LIMIT\nfrom .a import LIMIT as BOUND\n\n\ndef h():\n    from .a import LIMIT as late\n\n"
    "    return late\n",
    "c.py": "from .b import BOUND\n",
    "d.py": "from . import a\n\n\ndef k():\n    return a.LIMIT\n",
}


def _forget(package: str) -> None:
    for m in [m for m in sys.modules if m == package or m.startswith(f"{package}.")]:
        del sys.modules[m]


@pytest.fixture
def planted(tmp_path, monkeypatch):
    pkg = tmp_path / "plantpkg"
    pkg.mkdir()
    for name, text in PLANTED.items():
        (pkg / name).write_text(text, encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))
    _forget("plantpkg")
    yield pkg
    _forget("plantpkg")


def _load(*names):
    return [importlib.import_module(f"plantpkg.{n}") for n in names]


def test_the_patch_reaches_every_reader_of_a_planted_constant(perturb, planted):
    """Ruling Q29: the module's attribute, a from-import binder (under its own name and an alias), a re-export of it,
    a function-level import, a reader through the module, a default and a constant derived from it at import."""
    finder = perturb.install("a:LIMIT=4.5", package="plantpkg")
    try:
        a, b, c, d = _load("a", "b", "c", "d")
    finally:
        perturb.uninstall(finder)
    assert finder.loader.bound == 1
    assert (a.LIMIT, a.g(), a.f(), a.DOUBLE) == (4.5, 4.5, 4.5, 9.0)
    assert (b.LIMIT, b.BOUND, b.h(), c.BOUND, d.k()) == (4.5, 4.5, 4.5, 4.5, 4.5)
    assert not list((planted / "__pycache__").glob("a.*.pyc"))  # the rewritten module is never cached
    found = perturb.binders(planted, "a", "LIMIT", package="plantpkg")
    assert sorted((x["module"], x["as"], x["at"]) for x in found) == [
        ("b", "BOUND", "module"),
        ("b", "LIMIT", "module"),
        ("b", "late", "function"),
        ("c", "BOUND", "module"),
    ]


def test_without_the_variable_nothing_is_installed_and_nothing_changes(perturb, planted, monkeypatch):
    monkeypatch.delenv(perturb.ENV, raising=False)
    assert perturb.install(package="plantpkg") is None and perturb.install("  ", package="plantpkg") is None
    a, b = _load("a", "b")
    assert (a.LIMIT, a.f(), a.DOUBLE, b.BOUND) == (3.0, 3.0, 6.0, 3.0)


def test_reach_says_what_holds_the_value(perturb, planted):
    got = perturb.reach("a:LIMIT=4.5", src=planted, package="plantpkg")
    assert got["reached"] and got["module"] == 4.5 and got["statements"] == 1
    assert got["binders"] == {"b.LIMIT": 4.5, "b.BOUND": 4.5, "b.late": "read when called", "c.BOUND": 4.5}
    assert not any(isinstance(f, perturb._Finder) for f in sys.meta_path)  # reach uninstalls its finder


def test_install_refuses_what_it_cannot_perturb(perturb, planted):
    for bad in ("a:LIMIT", "LIMIT=3", "a:LIMIT=(1, 2)", "a:LIMIT=True", "a:LIMIT='x'", "a:1X=2"):
        with pytest.raises(ValueError):
            perturb.parse(bad)
    finder = perturb.install("c:NOT_THERE=1", package="plantpkg")
    try:
        with pytest.raises(ImportError, match="binds NOT_THERE"):
            _load("c")
    finally:
        perturb.uninstall(finder)
    _load("a")
    with pytest.raises(RuntimeError, match="imported already"):
        perturb.install("a:LIMIT=4.5", package="plantpkg")


def test_the_reading_lines_are_function_lines_naming_the_constant(perturb, planted):
    got = perturb.reading_lines(planted, ["a:LIMIT"], package="plantpkg")
    lines = sorted((Path(f).name, n) for f, n in got)
    # a.g's return, b.h's use of its function-level binder, d.k's attribute; not module lines, defaults, imports
    assert lines == [("a.py", 10), ("b.py", 8), ("d.py", 5)]
    assert all(v == {"a:LIMIT"} for v in got.values())


@pytest.mark.skipif(sys.version_info < (3, 12), reason="sys.monitoring")
def test_the_trace_says_which_constants_a_process_read(perturb, tmp_path):
    """FOLDMAP_TRACE_READS in a child: a constant read in a function that ran is recorded, one only defined or used as
    a default is not, and the record is written at exit."""
    pkg = tmp_path / "pkgs" / "tracepkg"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("", encoding="utf-8")
    (pkg / "m.py").write_text(
        "USED = 2.5\nIDLE = 4.0\nDEFAULT = 1.5\n\n\ndef f(x=DEFAULT):\n    return x * USED\n\n\ndef g():\n"
        "    return IDLE\n",
        encoding="utf-8",
    )
    folder = tmp_path / "reads"
    folder.mkdir()
    (folder / "names.json").write_text(json.dumps(["m:USED", "m:IDLE", "m:DEFAULT"]), encoding="utf-8")
    code = (
        "import importlib.util, sys\n"
        f"sys.path.insert(0, {str(pkg.parent)!r})\n"
        f"spec = importlib.util.spec_from_file_location('p', {str(SCRIPTS / 'perturb.py')!r})\n"
        "p = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(p)\n"
        "p.install_trace(package='tracepkg')\n"
        "from tracepkg import m\n"
        "print(m.f())\n"
    )
    env = {**os.environ, "FOLDMAP_TRACE_READS": str(folder)}
    done = subprocess.run([sys.executable, "-I", "-c", code], env=env, capture_output=True, text=True, timeout=120)
    assert done.returncode == 0, done.stderr
    records = [json.loads(p.read_text(encoding="utf-8")) for p in folder.glob("*.json") if p.name != "names.json"]
    assert len(records) == 1 and records[0]["reads"] == ["m:USED"]


def test_reads_are_counted_per_entry(sens, tmp_path):
    report = {"entries": [{"id": "1ABC", "name": "1ABC"}, {"id": "AF-X", "name": "AF_X"}]}
    for pid, argv, reads in (
        (1, ["r.py", "view", "--name", "AF_X", "--work", "w", "default", "f.cif"], ["m:A"]),
        (2, ["r.py", "entry", "1ABC", "--out", "o"], ["m:A", "m:B"]),
        (3, ["r.py", "run", "--ids", "1ABC"], ["m:B"]),
    ):
        (tmp_path / f"{pid}.json").write_text(json.dumps({"argv": argv, "reads": reads}), encoding="utf-8")
    (tmp_path / "names.json").write_text("[]", encoding="utf-8")
    assert sens.read_traces(tmp_path, report) == {"m:A": ["1ABC", "AF-X"], "m:B": ["(judging)", "1ABC"]}


def _child(code: str, perturbed: str | None) -> list:
    env = {k: v for k, v in os.environ.items() if k != "FOLDMAP_PERTURB"}
    if perturbed:
        env["FOLDMAP_PERTURB"] = perturbed
    done = subprocess.run([sys.executable, "-I", "-c", code], env=env, capture_output=True, text=True, timeout=300)
    assert done.returncode == 0, done.stderr[-2000:]
    return json.loads(done.stdout.strip().splitlines()[-1])


def test_random100_perturbs_foldmap_in_its_own_process_before_importing_it():
    """random100.py installs the hook before its first foldmap import, so every process of a run (the run, each
    entry, each view: all random100.py) draws with the value: layout._K_LIKE, the LIKE_WEIGHT derived from it (and
    like.py's from-import of that), and _attract's default held_k all move together."""
    code = (
        "import importlib.util, inspect, json\n"
        f"spec = importlib.util.spec_from_file_location('r100', {str(SCRIPTS / 'random100.py')!r})\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(m)\n"
        "from foldmap import layout, like\n"
        "held = inspect.signature(layout._attract).parameters['held_k'].default\n"
        "print(json.dumps([layout._K_LIKE, layout.LIKE_WEIGHT * layout._K_HOME, like.LIKE_WEIGHT * layout._K_HOME,"
        " held]))\n"
    )
    plain = _child(code, None)
    assert plain == pytest.approx([2.56] * 4)
    assert _child(code, "layout:_K_LIKE=1.92") == pytest.approx([1.92] * 4)


def test_perturb_reach_on_foldmap_in_a_child():
    env = {**os.environ, "FOLDMAP_PERTURB": "align:SAME_FOLD_TM=0.375"}
    done = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / "perturb.py"), "reach"], env=env, capture_output=True, text=True
    )
    got = json.loads(done.stdout)
    assert done.returncode == 0 and got["reached"] and got["module"] == 0.375
    assert got["binders"] == {"criteria.SAME_FOLD_TM": "read when called"}  # criteria imports it in a function


def test_a_perturbed_random100_run_needs_not_for_paper(tmp_path, monkeypatch):
    r100 = _script("random100")  # loaded unperturbed
    monkeypatch.setenv("FOLDMAP_PERTURB", "layout:_STACK_PAD=0.975")
    out = tmp_path / "o"
    with pytest.raises(SystemExit, match="not-for-paper"):
        r100.main(["run", "--ids", "1ABC", "--out", str(out), "--sample", str(tmp_path / "none.json")])
    assert not out.exists()


# ------------------------------------------------------------------------------------------------ the campaign


def test_perturbations_round_ints_and_skip_what_does_not_change(sens):
    assert sens.perturbed(0.6, 0.75) == (0.45, None) and sens.perturbed(0.6, 1.25) == (0.75, None)
    assert sens.perturbed(1.1, 1.25) == (1.375, None)  # no float noise
    assert sens.perturbed(6, 0.75) == (5, None) and sens.perturbed(6, 1.25) == (8, None)  # 4.5 -> 5, 7.5 -> 8
    assert sens.perturbed(2, 0.75)[0] is None and sens.perturbed(2, 1.25) == (3, None)
    assert sens.perturbed(1, 0.75)[0] is None and sens.perturbed(0.0, 1.25)[0] is None
    assert sens.perturbed((1.0, 2.0), 0.75)[1].startswith("not perturbable") and sens.perturbed(True, 1.25)[0] is None


@needs_records
def test_the_plan_is_the_inventorys_measured_and_named_entry_constants(sens):
    got = sens.plan()
    assert got["constants"] == 64  # ruling Q28: 67 less the three demoted git rows
    assert {s["id"] for s in got["not_perturbable"]} == {
        "cli:ROSETTE_ORDERS",
        "na_ends:_OFF_TURNS",
        "na_labels:_GAP_OUT",
        "na_labels:_GAP_TURNS",
        "na_runs:_JOINT_BENDS",
        "na_runs:_JOINT_SWINGS",
        "route:_NAME_ALONG",
    }
    assert {(s["id"], s["factor"]) for s in got["skipped"]} >= {
        ("align:SCREEN", 0.75),
        ("align:SCREEN", 1.25),
        ("membrane:JOIN_GAP", 0.75),
    }
    labels = [p["label"] for p in got["perturbations"]]
    assert len(labels) == len(set(labels)) == 2 * 57 - len(got["skipped"])


def _pairs(*ids):
    a, b = ids
    return [
        {"id": a, "stratum": "S11", "pair": a, "role": "seed", "partner": b},
        {"id": b, "stratum": "S11", "pair": a, "role": "partner", "partner": a},
    ]


def test_the_subset_is_stratified_fixed_keeps_pairs_and_leaves_out_the_slow(sens):
    entries = [{"id": f"S1-{k}", "stratum": "S1"} for k in range(6)] + [
        {"id": f"S2-{k}", "stratum": "S2"} for k in range(3)
    ]
    sample = {"entries": entries + _pairs("P-a", "P-b") + _pairs("Q-a", "Q-b"), "replaced": []}
    seconds = {e["id"]: 10.0 for e in sample["entries"]} | {"S1-0": 500.0, "Q-b": 300.0}
    got = sens.subset(sample, seconds, size=6, seed=1)
    assert got["left_out"] == {"S1-0": 500.0, "Q-a": 10.0, "Q-b": 300.0}
    assert got["quota"] == {"S1": 3, "S2": 2, "S11": 1}  # 6 x 5/10, 3/10, 2/10: 3, 1.8, 1.2 by largest remainder
    ids = got["ids"]
    assert len(ids) == 7 and {"P-a", "P-b"} <= set(ids) and not {"S1-0", "Q-a", "Q-b"} & set(ids)
    assert got == sens.subset(sample, seconds, size=6, seed=1) and got != sens.subset(sample, seconds, size=6, seed=2)


def _entry(i, ok=True, failing=(), sci=(), fb=None, area=10.0):
    return {
        "id": i,
        "pass": ok,
        "criteria": {k: {"pass": k not in failing} for k in ("A1", "A2")},
        "science": {"S1": {"pass": "S1" not in sci}, "S5": {"pass": None}},
        "views": {"default": {"result": {"fallbacks": fb or {}, "area": area, "svg_sha256": "x"}}},
        "process": {"seconds": 3.0},
    }


def test_the_classes_follow_the_ruling(sens):
    def keys(*entries):
        return sens.keys_of({"entries": list(entries)})

    base = keys(_entry("1"), _entry("2"), _entry("3", ok=False, failing=("A1",)))
    assert base["3"]["failing"] == ["A1"] and base["1"]["failing"] == []
    same = sens.judge(base, base)
    assert same["class"] == "stable" and not same["moved"]
    area = sens.judge(base, keys(_entry("1", area=11.0), _entry("2"), _entry("3", ok=False, failing=("A1",))))
    assert area["class"] == "sensitive" and area["area_moved"] == ["1"] and area["moved"] == ["1: area 10.0 -> 11.0"]
    fb = sens.judge(base, keys(_entry("1", fb={"loop_curve": 2}), _entry("2"), _entry("3", ok=False, failing=("A1",))))
    assert fb["class"] == "sensitive" and fb["fallbacks"] == {"loop_curve": [0, 2]}
    one = sens.judge(
        base, keys(_entry("1", ok=False, failing=("A1",)), _entry("2"), _entry("3", ok=False, failing=("A1",)))
    )
    assert one["class"] == "sensitive" and one["pass_delta"] == -1 and one["criteria"] == {"A1": [1, 2]}
    two = sens.judge(
        base,
        keys(
            _entry("1", ok=False, failing=("A1",)),
            _entry("2", ok=False, failing=("A1",)),
            _entry("3", ok=False, failing=("A1",)),
        ),
    )
    assert two["class"] == "fragile" and two["pass_delta"] == -2 and not two["new_failing"]
    new = sens.judge(base, keys(_entry("1"), _entry("2", sci=("S1",)), _entry("3", ok=False, failing=("A1",))))
    assert new["class"] == "fragile" and new["new_failing"] == ["S1"] and new["pass_delta"] == 0
    assert sens.worst(["stable", "fragile", "sensitive"]) == "fragile" and sens.worst([]) == "stable"


def test_sensitivity_never_writes_into_the_third_sample():
    sens = _script("sensitivity")
    with pytest.raises(SystemExit, match="refused"):
        sens.main(["run", "--stage", "A", "--out", str(VALIDATION / "random100-third" / "runs")])


@needs_report
def test_the_committed_report_is_stamped_and_its_classes_add_up():
    doc = json.loads((VALIDATION / "sensitivity" / "report.json").read_text(encoding="utf-8"))
    assert doc["source_stamp"]["script"]["name"] == "scripts/sensitivity.py"
    assert sum(doc["counts"].values()) == len(doc["constants"])
    assert doc["stage_a"]["determinism"]["identical"]
    assert doc["limitations"] == [c["id"] for c in doc["constants"] if c["class"] == "fragile"]
