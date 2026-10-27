"""The confidence of each inferred domain split (plan 28 Task 5, the review's 3.4, ruling Q42): residue_domains records
its cuts without changing the domains, auto_domains keeps the splits the drawn domains keep, each with its confidence
(1 - normalized cut / DOMAIN_NCUT), and the notes, `foldmap check` (text and S5's measures in --json), the MCP result
and the app say them; the figure is unchanged and a failure to say them never stops it."""

from __future__ import annotations

import json
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


def test_auto_domains_say_their_splits_and_stay_the_same():
    inputs = ds.load_fixture(FIXTURE)["1LDM_A"]
    bb = SimpleNamespace(labels=inputs.labels, ca=inputs.ca)
    plain = features.auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb)
    splits: list = []
    said = features.auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb, splits=splits)
    assert said == plain and len(plain) >= 2
    assert len(splits) == len(plain) - 1  # one split per boundary the drawn domains keep
    names = {name for name, _ in plain}
    for s in splits:
        assert set(s["between"][0]) | set(s["between"][1]) <= names and not set(s["between"][0]) & set(s["between"][1])
        assert 0 < s["confidence"] <= 1 and s["confidence"] == features.split_confidence(s["ncut"])
        assert s["residues"][0] > 0 and s["residues"][1] > 0 and s["sides"][0] and s["sides"][1]
    json.dumps(splits)  # plain JSON: the MCP result and check's --json carry it


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


def test_the_words():
    split = {
        "chain": "A",
        "between": [["D1"], ["D2", "D3"]],
        "ncut": 0.185,
        "confidence": 0.5,
        "across": 120,
        "density_ratio": 0.21,
        "residues": [80, 160],
        "sides": ["1-80", "81-240"],
    }
    got = domain_confidence.note([split])
    assert got.startswith("inferred domain splits, with their confidence (1 - (the split's normalized cut")
    assert "D1 | D2, D3 0.50 (normalized cut 0.18" in got and "120 contacts across, 0.21 of the contact density" in got
    assert domain_confidence.note([]) == "" and domain_confidence.lines([]) == []
    assert domain_confidence.note([{"error": "KeyError"}]) == "domain split confidence not judged: KeyError"
    head, line = domain_confidence.lines([split])
    assert head.startswith("  inferred domain splits") and line.startswith("    D1 | D2, D3 0.50")


def test_said_is_guarded():
    assert domain_confidence.said(object()) == ([], "domain split confidence not judged: AttributeError")
    assert domain_confidence.said(SimpleNamespace(layout=SimpleNamespace(domain_splits=[]))) == ([], "")


def test_the_check_report_holds_them_as_s5_measures():
    found = {"science": {}, "criteria": {}, "science_judged": {}}
    split = {"chain": "A", "between": [["D1"], ["D2"]], "confidence": 0.4}
    (s5,) = [c for c in deep.report("x.cif", found, splits=[split])["checks"] if c["code"] == "S5"]
    (bare,) = [c for c in deep.report("x.cif", found)["checks"] if c["code"] == "S5"]
    assert s5["measures"] == {"domain_splits": [split], "confidence": domain_confidence.CONFIDENCE}
    assert "measures" not in bare and {k: v for k, v in s5.items() if k != "measures"} == bare


def test_plot_says_them_and_the_file_does_not(tmp_path, capsys):
    out = tmp_path / "f.svg"
    assert main(["plot", str(T4L), "--domains", "auto", "-o", str(out)]) == 0
    err = capsys.readouterr().err
    assert "foldmap: inferred domain splits, with their confidence" in err and "D1 | D2 " in err
    svg = out.read_text(encoding="utf-8")
    assert "confidence" not in svg and "normalized cut" not in svg
    plain = tmp_path / "g.svg"
    assert main(["plot", str(T4L), "-o", str(plain)]) == 0
    assert "inferred domain splits" not in capsys.readouterr().err  # no inferred domains, nothing said


def test_the_figure_is_the_same_with_or_without_the_record(monkeypatch):
    """The record is taken beside the domains: with it emptied of every split the layout signature is the same."""
    from foldmap.spec import draw_spec, figure_spec

    spec = figure_spec(T4L, domains="auto")
    drawn = draw_spec(spec)
    assert drawn.layout.domain_splits
    monkeypatch.setattr(features, "domain_splits", lambda *a, **k: [])
    again = draw_spec(spec)
    assert again.layout.domain_splits == [] and again.layout.domains == drawn.layout.domains
    assert deep.signature(again.layout, again.loops) == deep.signature(drawn.layout, drawn.loops)


def test_a_failing_record_never_stops_a_figure(monkeypatch, tmp_path, capsys):
    def broken(*a, **k):
        raise RuntimeError("no")

    monkeypatch.setattr(features, "domain_splits", broken)
    out = tmp_path / "f.svg"
    assert main(["plot", str(T4L), "--domains", "auto", "-o", str(out)]) == 0 and out.is_file()
    assert "foldmap: domain split confidence not judged: RuntimeError" in capsys.readouterr().err


def test_check_says_them(capsys):
    main(["check", str(T4L), "--domains", "auto"])
    out = capsys.readouterr().out
    lines = out.splitlines()
    at = next(k for k, line in enumerate(lines) if line.startswith("S5 "))
    assert lines[at + 1].startswith("  inferred domain splits, with their confidence")
    main(["check", str(T4L), "--domains", "auto", "--json", "-"])
    doc = json.loads(capsys.readouterr().out)
    (s5,) = [c for c in doc["checks"] if c["code"] == "S5"]
    assert doc["version"] == 1 and s5["status"] == deep.UNAVAILABLE  # said, never judged: S5 as before
    assert s5["measures"]["domain_splits"][0]["between"] == [["D1"], ["D2"]]


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
    normalized cut, each failure has one class, and the class counts add up."""
    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    for s in doc["splits"]:
        assert s["confidence"] == features.split_confidence(s["ncut"]), s
    failed = doc["diagnosis"]["failures"]
    assert failed and all(f["class"] in doc["diagnosis"]["classes"] for f in failed)
    assert sum(doc["diagnosis"]["counts"]["class"].values()) == len(failed)
