"""`foldmap check --deep` and `--json` (plan 28 Task 7, ruling Q34(3)): D1 repeat run, D2 mkdssp, D3 RNApolis, each
pass, fail or unavailable with its kind; an unavailable check is never counted as passed; the default `check` output
is unchanged. The external programs are not needed: where they are not set up (CI) D2 and D3 are unavailable, and
their comparisons are tested on given references; the tests that run them for real are skipped there."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

from foldmap import deep
from foldmap.cli import main
from foldmap.science import CONSISTENCY, EXTERNAL, READABILITY
from foldmap.spec import draw_spec, figure_spec

DATA = Path(__file__).parent / "data"
UBQ = DATA / "1UBQ.cif"
TRNA = DATA / "na" / "1EHZ.cif.gz"


@pytest.fixture
def no_tools(monkeypatch):
    monkeypatch.delenv("MKDSSP", raising=False)
    monkeypatch.delenv("NA_RNAPOLIS_PYTHON", raising=False)


def _run(capsys, *argv) -> tuple[int, str]:
    code = main(["check", *map(str, argv)])
    return code, capsys.readouterr().out


def test_deep_appends_to_the_default_text(capsys, no_tools):
    """`check` prints what it printed before Task 7; --deep prints the same, then its kinds line and D1-D3."""
    code0, plain = _run(capsys, UBQ)
    code1, deeper = _run(capsys, UBQ, "--deep")
    assert code0 == code1 == 0
    assert "D1 " not in plain and "deep check kinds" not in plain
    assert deeper.startswith(plain)
    added = deeper[len(plain) :].splitlines()
    assert added[0].startswith("deep check kinds: D1 structural consistency; D2 external validation")
    assert added[1] == f"D1 {deep.DEEP['D1']}: pass"
    assert f"D2 {deep.DEEP['D2']}: unavailable (MKDSSP is not set" in deeper
    assert f"D3 {deep.DEEP['D3']}: does not apply (no nucleic acid)" in deeper


def test_json_report(capsys, no_tools):
    code, out = _run(capsys, UBQ, "--deep", "--json", "-")
    assert code == 0
    doc = json.loads(out)  # nothing but the report on standard output
    assert doc["schema"] == deep.SCHEMA == "foldmap-check" and doc["version"] == deep.VERSION == 1
    assert {"structure", "foldmap", "deep", "verdict", "counts", "kinds", "checks"} <= set(doc)
    assert doc["deep"] is True
    by = {c["code"]: c for c in doc["checks"]}
    assert [c["group"] for c in doc["checks"]].count("deep") == 3
    for c in doc["checks"]:
        assert {"code", "group", "title", "kind", "status"} <= set(c)
        assert c["status"] in deep.STATUSES
        assert c["kind"] in doc["kinds"]
        if c["status"] == deep.UNAVAILABLE:
            assert c["reason"], c  # an unavailable check says why
    assert by["S1"]["kind"] == CONSISTENCY and by["S5"]["kind"] == EXTERNAL and by["A3"]["kind"] == READABILITY
    assert by["D1"]["kind"] == CONSISTENCY and by["D2"]["kind"] == by["D3"]["kind"] == EXTERNAL
    assert by["D1"]["status"] == deep.PASS
    a, b = by["D1"]["measures"]["svg_sha256"]
    assert a == b and len(a) == 64
    assert by["D2"]["status"] == deep.UNAVAILABLE and "MKDSSP" in by["D2"]["reason"]
    assert by["D3"]["status"] == deep.NOT_APPLICABLE
    assert by["S5"]["status"] == by["A10"]["status"] == by["A12"]["status"] == deep.UNAVAILABLE
    assert sum(doc["counts"].values()) == len(doc["checks"])
    assert doc["counts"][deep.PASS] == sum(c["status"] == deep.PASS for c in doc["checks"])
    assert doc["verdict"] == "incomplete"  # never pass while a check is unavailable


def test_json_to_a_file_keeps_the_text(capsys, tmp_path, no_tools):
    path = tmp_path / "report.json"
    _, plain = _run(capsys, UBQ)
    code, out = _run(capsys, UBQ, "--json", path)
    assert code == 0 and out == plain + f"wrote {path}\n"
    doc = json.loads(path.read_text(encoding="utf-8"))
    assert doc["deep"] is False and not [c for c in doc["checks"] if c["group"] == "deep"]


def _found(science=None, criteria=None):
    sci = {"S1": [], "S2": [], "S3": [], "S4": [], "S5": [], "S6": []} | (science or {})
    return {
        "science": sci,
        "science_judged": dict.fromkeys(sci, 1),
        "criteria": {"A3": {"pass": True, "detail": ""}} | (criteria or {}),
    }


def test_unavailable_is_never_counted_as_passed():
    assert deep.report("x", _found())["verdict"] == deep.PASS
    doc = deep.report("x", _found(science={"S5": None}))
    assert doc["verdict"] == "incomplete" and doc["counts"][deep.UNAVAILABLE] == 1
    doc = deep.report("x", _found(), deep={"D2": deep.result(deep.UNAVAILABLE, "MKDSSP is not set")})
    assert doc["verdict"] == "incomplete" and doc["counts"][deep.PASS] == 7  # S1-S6, A3: the D2 is not among them
    doc = deep.report("x", _found(criteria={"A10": {"pass": None, "detail": "", "here": "needs mkdssp"}}))
    assert doc["verdict"] == "incomplete"
    assert deep.report("x", _found(science={"S1": ["a | b drawn side by side with no bridge"]}))["verdict"] == "fail"
    nothing = _found()
    nothing["science_judged"]["S2"] = 0
    doc = deep.report("x", nothing)
    assert {c["code"]: c["status"] for c in doc["checks"]}["S2"] == deep.NOT_JUDGED and doc["verdict"] == deep.PASS


def test_repeat_run():
    spec = figure_spec(str(UBQ))
    drawn = draw_spec(spec)
    assert deep.repeat_run(spec, drawn)["status"] == deep.PASS
    other = deep.repeat_run(spec, drawn, draw=lambda _: draw_spec(figure_spec(str(DATA / "1EMA.cif"))))
    assert other["status"] == deep.FAIL and "signature and its svg differ" in other["detail"]


def test_a_failed_deep_check_fails_the_command(capsys, monkeypatch, no_tools):
    monkeypatch.setattr(deep, "repeat_run", lambda *a, **k: deep.result(deep.FAIL, "drawn again, it differs"))
    code, out = _run(capsys, UBQ, "--deep")
    assert code == 1 and f"D1 {deep.DEEP['D1']}: fail" in out


def _reference_rows(name: str) -> list[tuple[str, int, str, str]]:
    sys.path.insert(0, str(Path(__file__).parent))
    from test_dssp_reference import reference

    return [(c, n, i, code) for (c, n, i), code in reference(DATA / "dssp" / f"{name}.ss").items()]


def test_dssp_check_on_given_rows():
    rows = _reference_rows("1UBQ")
    got = deep.dssp_check(UBQ, rows)
    assert got["status"] == deep.PASS and got["agree"] == 1.0
    wrong = [(c, n, i, "H") for c, n, i, _ in rows]
    assert deep.dssp_check(UBQ, wrong)["status"] == deep.FAIL
    assert deep.dssp_check(UBQ, "MKDSSP is not set")["status"] == deep.UNAVAILABLE
    assert deep.dssp_check(TRNA, protein=False)["status"] == deep.NOT_APPLICABLE


def test_mkdssp_unavailable_says_why(tmp_path, no_tools):
    assert deep.mkdssp_rows(UBQ).startswith("MKDSSP is not set")
    assert deep.mkdssp_rows(UBQ, str(tmp_path / "no-such-mkdssp")).startswith("cannot run mkdssp")
    assert deep.mkdssp_rows(UBQ, sys.executable).startswith("mkdssp 3.1.2 is required")


def test_parse_dssp():
    text = (
        "header\n  #  RESIDUE AA STRUCTURE BP1 BP2  ACC\n"
        "    1    1 A M              0   0  \n"
        "    2    2 A Q  E           0   0  \n"
        "    3        !              0   0  \n"
        "    4   10AB K  S           0   0  \n"
    )
    assert deep.parse_dssp(text) == [("A", 1, "", "-"), ("A", 2, "", "E"), ("B", 10, "A", "-")]


def test_pairs_check_on_given_pairs(no_tools):
    ours = deep.deposited_pairs(TRNA)
    assert len(ours) == 21  # tRNA-Phe's canonical and wobble pairs (tests/test_independent_kernels.py)
    got = deep.pairs_check(TRNA, ours)
    assert got["status"] == deep.PASS and got["f1"] == 1.0
    half = set(sorted(ours, key=sorted)[: len(ours) // 2])
    assert deep.pairs_check(TRNA, half)["status"] == deep.FAIL
    assert deep.pairs_check(TRNA)["status"] == deep.UNAVAILABLE
    assert deep.pairs_check(UBQ, nucleic=False)["status"] == deep.NOT_APPLICABLE


def test_rnapolis_unavailable_says_why(no_tools):
    assert deep.rnapolis_pairs(TRNA).startswith("NA_RNAPOLIS_PYTHON is not set")
    got = deep.rnapolis_pairs(TRNA, sys.executable)
    if not isinstance(got, str):
        pytest.skip("this Python has rnapolis")
    assert got.startswith("RNApolis could not read it")


@pytest.mark.skipif(not os.environ.get("MKDSSP"), reason="MKDSSP not set (mkdssp 3.1.2): D2 unavailable here")
def test_d2_with_mkdssp():
    got = deep.dssp_check(UBQ)
    assert got["status"] == deep.PASS and got["agree"] == 1.0, got


@pytest.mark.skipif(not os.environ.get("NA_RNAPOLIS_PYTHON"), reason="NA_RNAPOLIS_PYTHON not set: D3 unavailable")
def test_d3_with_rnapolis(tmp_path):
    """RNApolis reads the trimmed mmCIF fixtures no further than their entity sequences (IndexError in its
    get_one_letter_name: D3 is then unavailable, and says so); the same model written as PDB it reads."""
    import gemmi

    assert deep.pairs_check(TRNA)["status"] in (deep.PASS, deep.UNAVAILABLE)
    pdb = tmp_path / "1EHZ.pdb"
    gemmi.read_structure(str(TRNA)).write_pdb(str(pdb))
    got = deep.pairs_check(pdb)
    assert got["status"] == deep.PASS, got
