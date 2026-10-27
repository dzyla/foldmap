"""The confidence of each inferred domain split (plan 28 Task 5, the review's 3.4, rulings Q42, Q43):
residue_domains records its cuts without changing the domains, auto_domains keeps the splits the drawn domains keep,
each with its confidence (1 - normalized cut / DOMAIN_NCUT), and the notes, `foldmap check` (text and S5's measures in
--json), the MCP result and the app say them; the figure is unchanged and a failure to say them never stops it."""

from __future__ import annotations

import json
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from foldmap import deep, domain_confidence, features
from foldmap import domain_score as ds
from foldmap.cli import main

DATA = Path(__file__).parent / "data"
T4L = DATA / "2LZM.cif"  # T4 lysozyme: two domains, the C-terminal one holding residues 1-12 (test_domain_score)
FIXTURE = DATA / "domains" / "cath_inputs.npz"
FRESH_FIXTURE = DATA / "domains" / "fresh_inputs.npz"
RECORD = Path(__file__).resolve().parents[1] / "docs" / "validation" / "domain-uncertainty" / "report.json"


def _blob(n: int, centre) -> np.ndarray:
    """n points 3.8 Å apart on a snake through a cube about `centre` (test_domain_score's)."""
    side = int(np.ceil(n ** (1 / 3)))
    pts = [
        (x if (y + z * side) % 2 == 0 else side - 1 - x, y if z % 2 == 0 else side - 1 - y, z)
        for z in range(side)
        for y in range(side)
        for x in range(side)
    ]
    return np.array(pts[:n], float) * 3.8 + np.asarray(centre, float)


def test_the_record_of_cuts_leaves_the_domains_as_they_are():
    ca = np.vstack([_blob(80, (0, 0, 0)), _blob(80, (60, 0, 0))])
    cuts: list = []
    assert features.residue_domains(ca, cuts=cuts).tolist() == features.residue_domains(ca).tolist()
    made = [c for c in cuts if c["made"]]
    (cut,) = made
    assert cut["ncut"] < features.DOMAIN_NCUT and cut["across"] == 0 and cut["density_ratio"] == 0
    assert sorted([len(cut["a"]), len(cut["b"])]) == [80, 80] and len(cut["part"]) == 160
    assert len(cuts) == 3 and not any(c["made"] for c in cuts if c is not cut)  # each side weighed, left whole


def test_confidence_is_the_margin_below_the_limit():
    assert features.split_confidence(0.0) == 1.0
    assert features.split_confidence(features.DOMAIN_NCUT) == 0.0
    assert features.split_confidence(features.DOMAIN_NCUT / 2) == 0.5
    assert features.split_confidence(1.0) == 0.0  # never negative


@pytest.fixture(scope="module")
def t4l():
    """2LZM drawn with --domains auto, once for the tests that only read it."""
    from foldmap.spec import draw_spec, figure_spec

    return draw_spec(figure_spec(T4L, domains="auto"))


def _plain(inputs, **kw):
    bb = SimpleNamespace(labels=inputs.labels, ca=inputs.ca)
    return features.auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb, **kw)


def test_auto_domains_say_their_splits_and_stay_the_same():
    """On every chain of both domain fixtures (44), the domains are the same with the record of splits as without
    it, and each split said is one between the drawn domains."""
    chains = {**ds.load_fixture(FIXTURE), **ds.load_fixture(FRESH_FIXTURE)}
    assert len(chains) == 44
    split_chains = 0
    for key, inputs in chains.items():
        splits: list = []
        plain = _plain(inputs)
        assert _plain(inputs, splits=splits) == plain, key
        assert len(splits) <= max(len(plain) - 1, 0), key  # one split at most per boundary the drawn domains keep
        split_chains += bool(splits)
        names = {name for name, _ in plain}
        for s in splits:
            a, b = set(s["between"][0]), set(s["between"][1])
            assert a | b <= names and not a & b, key
            assert 0 <= s["confidence"] <= 1 and s["confidence"] == features.split_confidence(s["ncut"]), key
            assert s["residues"][0] > 0 and s["residues"][1] > 0 and s["sides"][0] and s["sides"][1], key
        json.dumps(splits)  # plain JSON: the MCP result and check's --json carry it
    assert split_chains >= 30


def test_which_cuts_the_drawn_domains_keep():
    """A cut is kept when the domains holding most of their residues in its part lie on both of its sides (each on
    the side holding more of it: an element the rules moved across keeps it); undone by a join, it is none; one
    parting the same domains again is said once."""
    labels = [SimpleNamespace(seq=k + 1, icode="") for k in range(100)]
    sses = [SimpleNamespace(id=f"E{k}", chain="A", start=10 * k, end=10 * k + 7) for k in range(10)]
    cut = {"made": True, "a": np.arange(50), "b": np.arange(50, 100), "ncut": 0.1, "across": 3, "density_ratio": 0.1}
    apart = [("D1", [f"E{k}" for k in range(5)]), ("D2", [f"E{k}" for k in range(5, 10)])]
    (got,) = features.domain_splits({"A": [cut]}, apart, sses, labels)
    assert got["between"] == [["D1"], ["D2"]] and got["sides"] == ["1-50", "51-100"] and got["confidence"] == 0.73
    moved = [("D1", [f"E{k}" for k in range(6)]), ("D2", [f"E{k}" for k in range(6, 10)])]
    assert [s["between"] for s in features.domain_splits({"A": [cut]}, moved, sses, labels)] == [[["D1"], ["D2"]]]
    joined = [("D1", [f"E{k}" for k in range(10)])]
    assert features.domain_splits({"A": [cut]}, joined, sses, labels) == []
    inner = {**cut, "a": np.arange(30), "b": np.arange(30, 50), "ncut": 0.2}  # parts D1's half: D2 is not in it
    assert features.domain_splits({"A": [inner]}, apart, sses, labels) == []
    again = {**cut, "a": np.arange(40), "b": np.arange(40, 100), "ncut": 0.3}
    assert len(features.domain_splits({"A": [cut, again]}, apart, sses, labels)) == 1
    assert features.domain_splits({"A": [{**cut, "made": False}]}, apart, sses, labels) == []
    assert features.domain_splits({"A": [{"error": "ValueError"}]}, apart, sses, labels) == [
        {"chain": "A", "error": "ValueError"}
    ]


def test_a_copys_domains_are_not_in_its_reference_chains_split():
    """Several chains split: names carry " · chain", and a cut of chain A parts A's domains only, though its copy
    B's domains share their numbers (auto_domains names them D1 · A, D1 · B ...)."""
    labels = [SimpleNamespace(seq=k + 1, icode="") for k in range(100)] * 2
    sses = [SimpleNamespace(id=f"A{k}", chain="A", start=10 * k, end=10 * k + 7) for k in range(10)]
    sses += [SimpleNamespace(id=f"B{k}", chain="B", start=100 + 10 * k, end=100 + 10 * k + 7) for k in range(10)]
    cut = {"made": True, "a": np.arange(50), "b": np.arange(50, 100), "ncut": 0.1, "across": 3, "density_ratio": 0.1}
    named = [
        ("D1 · A", [f"A{k}" for k in range(5)]),
        ("D2 · A", [f"A{k}" for k in range(5, 10)]),
        ("D1 · B", [f"B{k}" for k in range(5)]),
        ("D2 · B", [f"B{k}" for k in range(5, 10)]),
    ]
    (got,) = features.domain_splits({"A": [cut]}, named, sses, labels)
    assert got["chain"] == "A" and got["between"] == [["D1 · A"], ["D2 · A"]]
    assert domain_confidence.describe(got).startswith("D1 · A | D2 · A 0.73")


SPLIT = {
    "chain": "A",
    "between": [["D1"], ["D2", "D3"]],
    "ncut": 0.185,
    "confidence": 0.5,
    "across": 120,
    "density_ratio": 0.21,
    "residues": [80, 160],
    "sides": ["1-80", "81-240"],
}


def test_the_words():
    """Each interface says the limit's value, that it is no probability and that it moves with the limit."""
    said = domain_confidence.definition()
    assert f"1 - normalized cut / {features.DOMAIN_NCUT:g} (DOMAIN_NCUT" in said
    assert "not a probability" in said and "moves with the limit" in said
    got = domain_confidence.note([SPLIT])
    assert got.startswith(f"inferred domain splits, with their confidence ({said}; said, not judged): ")
    assert "D1 | D2, D3 0.50 (normalized cut 0.18" in got and "120 contacts across, 0.21 of the contact density" in got
    assert domain_confidence.note([]) == "" and domain_confidence.lines([]) == []
    assert domain_confidence.note([{"error": "KeyError"}]) == "domain split confidence not judged: KeyError"
    head, line = domain_confidence.lines([SPLIT])
    assert head.startswith("  inferred domain splits") and line.startswith("    D1 | D2, D3 0.50")


def test_the_limit_followed(monkeypatch):
    monkeypatch.setattr(features, "DOMAIN_NCUT", 0.3)
    assert "1 - normalized cut / 0.3 (DOMAIN_NCUT" in domain_confidence.definition()


def test_one_failed_chain_hides_no_split():
    mixed = [SPLIT, {"chain": "B", "error": "KeyError"}]
    got = domain_confidence.note(mixed)
    assert "D1 | D2, D3 0.50" in got and got.endswith("; domain split confidence not judged for chain B: KeyError")
    said = domain_confidence.lines(mixed)
    assert (
        said[1].startswith("    D1 | D2, D3")
        and said[-1] == "  domain split confidence not judged for chain B: KeyError"
    )


def test_said_is_guarded():
    assert domain_confidence.said(object()) == (
        [{"error": "AttributeError"}],
        "domain split confidence not judged: AttributeError",
    )
    assert domain_confidence.said(SimpleNamespace(layout=SimpleNamespace(domain_splits=[]))) == ([], "")


def test_the_split_lines_go_under_s5_after_its_findings():
    science = [
        "S4 x: pass",
        "S5 automatic domains agree with a reference: 1 problem(s)",
        "  A: overlap 0.5",
        "S6 y: pass",
    ]
    got = domain_confidence.under_s5(science, [SPLIT])
    assert got[:3] == science[:3] and got[3].startswith("  inferred domain splits") and got[-1] == "S6 y: pass"
    unmeasured = ["S5 automatic domains agree with a reference: not yet measured", "S6 y: pass"]
    got = domain_confidence.under_s5(unmeasured, [SPLIT])
    assert got[0] == f"S5 automatic domains agree with a reference: {domain_confidence.S5_LINE}"
    assert got[1].startswith("  inferred domain splits") and got[-1] == "S6 y: pass"
    assert domain_confidence.under_s5(unmeasured, []) == unmeasured
    failed = domain_confidence.under_s5(unmeasured, [{"error": "RuntimeError"}])
    assert failed[1] == "  domain split confidence not judged: RuntimeError"


def test_the_check_report_holds_them_as_s5_measures():
    """Within schema v1: S5's optional "measures", its reason saying what they are; no top-level key, no status
    change."""
    found = {"science": {}, "criteria": {}, "science_judged": {}}
    split = {"chain": "A", "between": [["D1"], ["D2"]], "confidence": 0.4}
    doc = deep.report("x.cif", found, splits=[split])
    (s5,) = [c for c in doc["checks"] if c["code"] == "S5"]
    bare_doc = deep.report("x.cif", found)
    (bare,) = [c for c in bare_doc["checks"] if c["code"] == "S5"]
    assert s5["measures"] == {"domain_splits": [split], "confidence": domain_confidence.definition()}
    assert s5["status"] == bare["status"] == deep.UNAVAILABLE and s5["reason"] == domain_confidence.S5_REASON
    assert "measures" not in bare and bare["reason"] == "not yet measured"
    assert set(doc) == set(bare_doc) and doc["version"] == 1


def test_plot_says_them_and_the_file_does_not(tmp_path, capsys):
    out = tmp_path / "f.svg"
    assert main(["plot", str(T4L), "--domains", "auto", "-o", str(out)]) == 0
    err = capsys.readouterr().err
    assert "foldmap: inferred domain splits, with their confidence" in err and "D1 | D2 " in err
    assert "not a probability" in err
    svg = out.read_text(encoding="utf-8")
    assert "confidence" not in svg and "normalized cut" not in svg
    plain = tmp_path / "g.svg"
    assert main(["plot", str(T4L), "-o", str(plain)]) == 0
    assert "inferred domain splits" not in capsys.readouterr().err  # no inferred domains, nothing said


def test_the_figure_is_the_same_with_or_without_the_record(t4l, monkeypatch):
    """The record is taken beside the domains: with it emptied of every split the layout signature is the same."""
    from foldmap.spec import draw_spec, figure_spec

    assert t4l.layout.domain_splits
    monkeypatch.setattr(features, "domain_splits", lambda *a, **k: [])
    again = draw_spec(figure_spec(T4L, domains="auto"))
    assert again.layout.domain_splits == [] and again.layout.domains == t4l.layout.domains
    assert deep.signature(again.layout, again.loops) == deep.signature(t4l.layout, t4l.loops)


def test_said_of_a_drawn_figure(t4l):
    splits, line = domain_confidence.said(t4l)
    (split,) = splits
    assert split["between"] == [["D1"], ["D2"]] and line == domain_confidence.note(splits)


def test_a_failing_record_never_stops_a_figure(monkeypatch, tmp_path, capsys):
    def broken(*a, **k):
        raise RuntimeError("no")

    monkeypatch.setattr(features, "domain_splits", broken)
    out = tmp_path / "f.svg"
    assert main(["plot", str(T4L), "--domains", "auto", "-o", str(out)]) == 0 and out.is_file()
    assert "foldmap: domain split confidence not judged: RuntimeError" in capsys.readouterr().err


def test_check_says_them(tmp_path, capsys):
    report = tmp_path / "check.json"
    main(["check", str(T4L), "--domains", "auto", "--json", str(report)])
    lines = capsys.readouterr().out.splitlines()
    at = next(k for k, line in enumerate(lines) if line.startswith("S5 "))
    assert lines[at].endswith(domain_confidence.S5_LINE)
    assert lines[at + 1].startswith("  inferred domain splits, with their confidence")
    doc = json.loads(report.read_text(encoding="utf-8"))
    (s5,) = [c for c in doc["checks"] if c["code"] == "S5"]
    assert doc["version"] == 1 and s5["status"] == deep.UNAVAILABLE  # said, never judged: S5 as before
    assert s5["measures"]["domain_splits"][0]["between"] == [["D1"], ["D2"]]


def test_check_says_a_failure_to_say_them(monkeypatch, tmp_path, capsys):
    def broken(drawn):
        raise RuntimeError("no")

    monkeypatch.setattr(domain_confidence, "splits_of", broken)
    report = tmp_path / "check.json"
    main(["check", str(T4L), "--domains", "auto", "--json", str(report)])
    assert "  domain split confidence not judged: RuntimeError" in capsys.readouterr().out.splitlines()
    (s5,) = [c for c in json.loads(report.read_text(encoding="utf-8"))["checks"] if c["code"] == "S5"]
    assert s5["measures"]["domain_splits"] == [{"error": "RuntimeError"}]


def test_the_mcp_result_says_them(monkeypatch):
    pytest.importorskip("mcp")
    monkeypatch.delenv("FOLDMAP_OUTPUT_DIR", raising=False)
    import asyncio

    from foldmap.mcp_server import server

    r = asyncio.run(server.call_tool("draw_topology", {"source": str(T4L), "domains_from": "auto"}))
    (split,) = r.structured_content["domain_splits"]
    assert split["between"] == [["D1"], ["D2"]] and "inferred domain splits" in r.content[0].text


def test_the_app_says_them():
    from foldmap import app, app_jobs

    blob = json.dumps({**app.LAYOUT_DEFAULTS, "theme": "publication", "domains": "auto"}, sort_keys=True)
    out = app_jobs.draw_job(T4L, blob)
    assert out["splits_note"].startswith("inferred domain splits, with their confidence")


@pytest.mark.skipif(not RECORD.is_file(), reason="docs/validation is not part of the sdist")
def test_the_record_says_what_the_code_says():
    """The committed record (scripts/domain_uncertainty.py): each split's confidence is split_confidence of its
    normalized cut, each failure has one class and the class counts add up, the domains were compared with b72f65d's
    on every chain and none differs, and each bin's interval is the Wilson interval of its counts."""
    import importlib.util

    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    for s in doc["splits"]:
        assert s["confidence"] == features.split_confidence(s["ncut"]), s
    failed = doc["diagnosis"]["failures"]
    assert failed and all(f["class"] in doc["diagnosis"]["classes"] for f in failed)
    assert sum(doc["diagnosis"]["counts"]["class"].values()) == len(failed)
    assert doc["against"]["against"] == "b72f65d" and doc["against"]["differ"] == 0
    assert doc["against"]["chains"] == sum(doc["set_sizes"].values())
    spec = importlib.util.spec_from_file_location("_t5_common", RECORD.parents[3] / "scripts" / "paper" / "common.py")
    common = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(common)
    for name, got in doc["calibration"].items():
        for b in got.get("bins", []):
            if b["n"]:
                lo, hi = common.wilson(b["right"], b["n"])
                assert b["wilson95"] == [round(lo, 2), round(hi, 2)], name


# ---- the typed copies of the limit and its sensitivity (plan 28 final review F8) --------------------------------

ROOT = Path(__file__).resolve().parents[1]
SENSITIVITY = ROOT / "docs" / "validation" / "sensitivity" / "report.json"
LIMIT = re.compile(r"normalized cut(?: at|\s*/)\s*(\d+\.\d+)")
MOVES = re.compile(r"S5 failures (?:from )?(\d+) (?:to|→) (\d+)")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text).replace("`", "")


def _texts() -> dict[str, str]:
    """Every text that types DOMAIN_NCUT's value or its sensitivity: the check help, the MCP tool's docstring (its
    source: the mcp package may be absent), this module's docstring, the record script's text, and the README and
    the model page where they ship."""
    from foldmap.cli import build_parser

    sub = next(a for a in build_parser()._actions if a.dest == "command")
    out = {
        "check help": sub.choices["check"].description,
        "mcp_server.py": (ROOT / "src" / "foldmap" / "mcp_server.py").read_text(encoding="utf-8"),
        "domain_confidence": domain_confidence.__doc__,
    }
    for rel in ("README.md", "docs/model/domains.md", "scripts/domain_uncertainty.py"):
        if (ROOT / rel).is_file():
            out[rel] = (ROOT / rel).read_text(encoding="utf-8")
    return {k: _flat(v) for k, v in out.items()}


def test_the_typed_limit_is_the_constant():
    texts = _texts()
    for name in ("check help", "mcp_server.py") + (("README.md",) if "README.md" in texts else ()):
        assert LIMIT.search(texts[name]), f"{name} no longer says the limit"
    for name, text in texts.items():
        for got in LIMIT.findall(text):
            assert got == f"{features.DOMAIN_NCUT:g}", f"{name} types the limit as {got}"


@pytest.mark.skipif(not SENSITIVITY.is_file(), reason="docs/validation is not part of the sdist")
def test_the_typed_sensitivity_is_the_records():
    """ "x1.25 moves the third sample's S5 failures from 12 to 19" is stage B's S5 count for DOMAIN_NCUT x1.25 in the
    sensitivity record; a rerun of that record that moves it fails here until every copy is changed."""
    doc = json.loads(SENSITIVITY.read_text(encoding="utf-8"))
    row = next(c for c in doc["constants"] if c["id"] == "features:DOMAIN_NCUT")
    want = tuple(str(n) for n in row["runs"]["x1.25"][row["stage"]]["criteria"]["S5"])
    texts = _texts()
    for name in ("mcp_server.py", "domain_confidence") + tuple(n for n in ("README.md",) if n in texts):
        assert MOVES.search(texts[name]), f"{name} no longer says the sensitivity"
    for name, text in texts.items():
        for got in MOVES.findall(text):
            assert got == want, f"{name} says S5 failures {got[0]} to {got[1]}; the record {want[0]} to {want[1]}"
