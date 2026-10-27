"""Pruning on evidence (plan 28 Task 9.3, ruling Q46): every candidate is planned once, each disabling value is one
scripts/perturb.py's hook sets on a name its module binds, the verdicts follow the prune condition, and the committed
record is stamped and adds up."""

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
    or not (VALIDATION / "sensitivity" / "report.json").exists(),
    reason="the validation records are not part of the sdist",
)
needs_report = pytest.mark.skipif(
    not (VALIDATION / "prune" / "report.json").exists(), reason="the prune record is not part of the sdist"
)


def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"foldmap_{name}_prune_t", SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def prune():
    return _script("prune")


@pytest.fixture(scope="module")
def perturb():
    return _script("perturb")


def test_the_tables_are_disjoint_and_a_shared_value_is_tried(prune):
    tables = [set(prune.DISABLE), set(prune.SHARED), set(prune.NOT_TRIED)]
    assert all(not (a & b) for i, a in enumerate(tables) for b in tables[i + 1 :])
    assert set(prune.SHARED.values()) <= set(prune.DISABLE)
    assert set(prune.PRUNED) <= set(prune.DISABLE) and set(prune.HELD) <= set(prune.DISABLE)
    assert not set(prune.PRUNED) & set(prune.HELD)


def test_what_was_pruned_is_gone_from_the_code(prune, perturb):
    """Q46: a pruned special case is removed by a code commit: its module no longer binds the name."""
    import ast

    for cid in prune.PRUNED:
        module, name = cid.split(":", 1)
        tree = ast.parse((ROOT / "src" / "foldmap" / f"{module}.py").read_text(encoding="utf-8"))
        assert perturb.rewrite(tree, name, 0) == 0, cid


def test_every_disabling_value_is_a_number_the_hook_sets_on_a_name_its_module_binds(prune, perturb):
    """No silent no-op: each value parses (an int or a float, ruling Q29) and its module binds the name at top level
    (perturb.rewrite inserts the value after that binding)."""
    import ast

    for cid in set(prune.DISABLE) - set(prune.PRUNED):
        module, name, value = perturb.parse(prune.perturb_text(cid))
        assert value == prune.DISABLE[cid][0] and isinstance(value, (int, float)) and not isinstance(value, bool)
        tree = ast.parse((ROOT / "src" / "foldmap" / f"{module}.py").read_text(encoding="utf-8"))
        assert perturb.rewrite(tree, name, value) >= 1, cid
        assert prune.label(cid) == prune._script("sensitivity").safe_name(prune.label(cid))


@needs_records
def test_every_candidate_is_planned_once_and_no_plan_is_stale(prune):
    """Q46's candidates: the inventory's named-entry rows and the sensitivity study's not exercised constants."""
    items = prune.candidates()
    ids = [c["id"] for c in items]
    assert len(ids) == len(set(ids))
    report = json.loads((VALIDATION / "sensitivity" / "report.json").read_text(encoding="utf-8"))
    assert set(report["not_exercised"]) <= set(ids)
    assert all(c["origin"].startswith("named-entry") or c["not_exercised"] for c in items)
    p = prune.plan(items)
    assert p["unplanned"] == []
    assert {*prune.DISABLE, *prune.SHARED, *prune.NOT_TRIED} == set(ids)
    assert len(p["tried"]) + len(p["shared"]) + len(p["not_tried"]) == len(items)


@needs_records
def test_run_refuses_the_third_and_fourth_samples(prune):
    for where in (VALIDATION / "random100-third" / "x", VALIDATION / "random100-fourth" / "x"):
        args = type("A", (), {"out": str(where), "measure": "bench", "labels": None, "jobs": 1, "parallel": 1})
        with pytest.raises(SystemExit, match="refused"):
            prune.cmd_run(args)


def _entry(passed=True, failing=(), fallbacks=None, area=10.0):
    return {"pass": passed, "failing": list(failing), "fallbacks": fallbacks or {}, "area": area}


def test_a_campaign_verdict_orders_only_what_a_judge_orders(prune):
    base = {"A": _entry(), "B": _entry(False, ["A17"], {"loop_curve": 2})}
    assert prune.campaign_verdict(base, dict(base))["verdict"] == "identical"
    assert prune.campaign_verdict(base, {**base, "B": _entry(True, [], {"loop_curve": 1})})["verdict"] == "better"
    worse = prune.campaign_verdict(base, {**base, "A": _entry(False, ["A6"])})
    assert worse["verdict"] == "worse" and "A: PASS -> FAIL" in worse["worse"] and "A: A6 now fails" in worse["worse"]
    assert prune.campaign_verdict(base, {**base, "A": _entry(area=11.0)})["verdict"] == "moved"
    mixed = prune.campaign_verdict(base, {"A": _entry(area=11.0), "B": _entry(False, ["A17"], {"loop_curve": 3})})
    assert mixed["verdict"] == "worse"  # worse outranks moved
    assert prune.campaign_verdict(base, {"A": _entry()})["verdict"] == "worse"  # an entry lost


def test_the_other_verdicts(prune):
    s = {"sets": {"dom_cal": {"1ABC_A": {"overlap": 0.9, "f1": 1.0, "recall": 1.0, "meets": True}}}}
    lower = {"sets": {"dom_cal": {"1ABC_A": {"overlap": 0.8, "f1": 1.0, "recall": 1.0, "meets": True}}}}
    assert prune.domsets_verdict(s, s)["verdict"] == "identical"
    assert prune.domsets_verdict(s, lower)["verdict"] == "worse"
    assert prune.domsets_verdict(lower, s)["verdict"] == "better"
    bench = {"exit": 0, "changed": []}
    assert prune.bench_verdict(bench, bench)["verdict"] == "identical"
    assert prune.bench_verdict(bench, {"exit": 1, "changed": ["1UBQ"]})["verdict"] == "moved"
    suite = {"tests": 3, "failures": []}
    assert prune.suite_verdict(suite, suite)["verdict"] == "identical"
    assert prune.suite_verdict(suite, {"tests": 3, "failures": ["t::x"]})["worse"] == ["t::x"]
    assert prune.suite_verdict(suite, {"tests": 2, "failures": []})["verdict"] == "worse"  # tests not run


def test_prune_only_when_every_measure_is_identical_or_better(prune):
    ok, bad = {"verdict": "identical"}, {"verdict": "moved", "moved": ["X: area 1 -> 2"]}
    wanted = ["r100", "suite"]
    assert prune.decide({"r100": ok, "suite": {"verdict": "better", "better": ["t"]}}, wanted)["decision"] == "prune"
    kept = prune.decide({"r100": bad, "suite": ok}, wanted)
    assert kept["decision"] == "keep" and kept["kept_by"][0]["measure"] == "r100"
    assert prune.decide({"r100": ok, "suite": None}, wanted)["decision"] == "incomplete"
    assert (
        prune.decide({"r100": {"verdict": "incomplete", "why": "no report"}, "suite": ok}, wanted)["decision"] == "keep"
    )


def test_junit_names_the_tests_that_failed_or_erred(prune, tmp_path):
    xml = tmp_path / "j.xml"
    xml.write_text(
        '<testsuites><testsuite><testcase classname="tests.test_a" name="test_x"/>'
        '<testcase classname="tests.test_a" name="test_y"><failure message="m"/></testcase>'
        '<testcase classname="tests.test_b" name="test_z"><error message="e"/></testcase>'
        '<testcase classname="tests.test_b" name="test_s"><skipped/></testcase></testsuite></testsuites>',
        encoding="utf-8",
    )
    got = prune.junit(xml)
    assert got == {"tests": 4, "failed": 2, "skipped": 1, "failures": ["tests.test_a::test_y", "tests.test_b::test_z"]}
    assert prune.junit(tmp_path / "none.xml") == {"error": "no junit xml"}


def test_the_suite_plugin_installs_the_value_before_foldmap_is_imported(prune, tmp_path):
    (tmp_path / "foldmap_prune_plugin.py").write_text(
        prune.PLUGIN.format(perturb=str(SCRIPTS / "perturb.py"), env=prune.PRUNE_ENV), encoding="utf-8"
    )
    env = {k: v for k, v in os.environ.items() if k not in ("FOLDMAP_PERTURB", prune.PRUNE_ENV)}
    env |= {"PYTHONPATH": str(tmp_path), prune.PRUNE_ENV: "symmetry:MIN_HELIX_COPIES=3"}
    code = "import foldmap_prune_plugin\nimport foldmap.symmetry as s\nprint(s.MIN_HELIX_COPIES)"
    done = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "3"
    env.pop(prune.PRUNE_ENV)
    done = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, cwd=tmp_path)
    assert done.stdout.strip() == "4"  # without the variable the plugin installs nothing


def test_prune_refuses_to_run_perturbed_itself(prune, monkeypatch):
    monkeypatch.setenv("FOLDMAP_PERTURB", "layout:_STACK_PAD=1.0")
    with pytest.raises(SystemExit):
        prune.main(["plan"])


@needs_report
def test_the_committed_record_is_stamped_and_adds_up(prune):
    doc = json.loads((VALIDATION / "prune" / "report.json").read_text(encoding="utf-8"))
    stamp = doc["source_stamp"]
    assert stamp["git"]["commit"] and not stamp["git"].get("dirty")
    assert doc["code_fingerprint"]
    c = doc["counts"]
    rows = doc["candidates"] + doc["shared"]
    assert c["candidates"] == len(rows) + len(doc["not_tried"])
    assert c["tried"] == len(rows) and c["not_tried"] == len(doc["not_tried"])
    assert c["kept"] + c["pruned"] + c["held"] + c["incomplete"] == c["tried"]
    assert c["incomplete"] == 0
    assert sorted(doc["pruned"]) == sorted(r["id"] for r in doc["candidates"] if r["decision"] == "prune")
    for r in doc["candidates"]:
        if r["decision"] == "keep":
            assert r["kept_by"], r["id"]
        if r["decision"] in ("prune", "held"):
            assert all(v and v["verdict"] in ("identical", "better") for v in r["measures"].values()), r["id"]
    assert doc["baselines"]["r100_deterministic"]["identical"]
    assert doc["removed"] == sorted(doc["pruned"])  # every special case that met the condition and was not held
    assert set(doc["held"]) == {r["id"] for r in doc["candidates"] if r["decision"] == "held"}
    # the md is the JSON's (render), and the record is the code's tables (review finding 4)
    assert prune.render(doc) == (VALIDATION / "prune" / "report.md").read_text(encoding="utf-8")
    every = doc["candidates"] + doc["shared"] + doc["not_tried"]
    assert {r["id"] for r in every} == {*prune.DISABLE, *prune.SHARED, *prune.NOT_TRIED}
    assert all(r["off"] == prune.DISABLE[r["id"]][0] for r in doc["candidates"])
    assert {s["id"]: s["with"] for s in doc["shared"]} == prune.SHARED
    assert set(doc["held"]) == set(prune.HELD) and doc["removed"] == sorted(prune.PRUNED)
    assert {k: c[k] for k in ("candidates", "tried", "kept", "pruned", "held", "not_tried")} == {
        "candidates": 45,
        "tried": 32,
        "kept": 28,
        "pruned": 2,
        "held": 2,
        "not_tried": 13,
    }
    # the pruned code draws what was measured: the post-prune baselines byte-identical to each "-off" run
    eq = doc["equivalence"]
    assert eq["runs_code"] == doc["code_fingerprint"] and eq["after_code"] != eq["runs_code"]
    assert set(eq["measures"]) == {"r100", "domains", "na"}
    for g in eq["measures"].values():
        assert g["off"] and all(o["identical"] == g["entries"] and not o["differ"] for o in g["off"].values())
    name_open = next(r for r in doc["candidates"] if r["id"] == "layout:_NAME_OPEN")
    assert name_open["measures"]["domains"]["svg"] == ["5NF8", "7F5S", "9FZL"]  # live logic, keys identical (Q49)
