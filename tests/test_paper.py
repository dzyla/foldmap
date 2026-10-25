"""The paper pipeline in scripts/paper (plan 28 Task 1): Table S1 is read from the code's constants (every value
found, by name), the nucleic-acid report is read in both its formats, and the fourth random sample's sampler keeps
its protocol: a fixed seed and cutoff, no id of samples 1-3, drawn once."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "scripts" / "paper"
SAMPLES = ROOT / "docs" / "validation"
needs_samples = pytest.mark.skipif(
    not (SAMPLES / "random100-third" / "sample.json").exists(), reason="the samples are not part of the sdist"
)


def _load(name: str):
    if str(PAPER) not in sys.path:
        sys.path.insert(0, str(PAPER))  # the scripts import common from their own folder
    spec = importlib.util.spec_from_file_location(f"paper_{name}", PAPER / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # a dataclass resolves its annotations through its module
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------------------------- Table S1


@pytest.fixture(scope="module")
def table_s1():
    return _load("table_s1")


def test_table_s1_finds_every_value_in_this_code(table_s1):
    """A renamed or reworded constant stops the table instead of dropping out of it."""
    rows, values = table_s1.build(table_s1.Source(ROOT / "src" / "foldmap"))
    assert len(rows) == len(table_s1.ROWS) and all(r.endswith(r" \\") for r in rows)
    assert len(values) == sum(len(refs) for *_, refs in table_s1.ROWS)
    by_row = {r.split(" & ")[1]: r.split(" & ")[2].removesuffix(r" \\") for r in rows}
    assert by_row["node limits"] == "200,000; 1,600,000; 8,000,000"
    assert by_row[r"energy constant; bond cut; floor; C$\alpha$ cut-off; closest atoms"].startswith(
        r"$0.084\times332$; $-0.5$; $-9.9$"
    )
    assert by_row[r"dyad in the page; dyad turns; turn gain"] == r"10$^\circ$; 0 to $\pm 75^\circ$ by 15$^\circ$; 0.01"
    assert by_row["grid"] == "$10^{-4}$ page units"
    assert by_row[r"$\beta\alpha$ row: strands; share of units"] == "4--10; 3 of 5"


def test_table_s1_reads_constants_defaults_and_inline_numbers(table_s1, tmp_path):
    (tmp_path / "m.py").write_text(
        "import numpy as np\n"
        "A = 2.5\n"
        "B, C = 3, (1.0, 2.0)\n"
        "D = 2 * A\n"
        "E = np.radians(10)\n"
        "F = {'x': 0.8}\n"
        "def f(x, n=300, *, tol=1e-5):\n"
        "    return len(x) >= 16\n",
        encoding="utf-8",
    )
    src = table_s1.Source(tmp_path)
    C, D, Inline = table_s1.C, table_s1.D, table_s1.Inline
    assert src.read(C("m", "A"))[0] == 2.5
    assert src.read(C("m", "C"))[0] == (1.0, 2.0)
    assert src.read(C("m", "D"))[0] == 5.0
    assert table_s1.fmt(src.read(C("m", "E"))[0], "deg") == "10"
    assert src.read(C("m", "F", key="x"))[0] == 0.8
    assert src.read(D("m", "f", "n"))[0] == 300 and src.read(D("m", "f", "tol"))[0] == 1e-5
    assert src.read(Inline("m", "f", r"len\(x\) >= (\d+)"))[0] == 16
    with pytest.raises(LookupError):
        src.read(C("m", "GONE"))
    with pytest.raises(LookupError):
        src.read(Inline("m", "f", r"len\(y\) >= (\d+)"))


def test_table_s1_compares_with_a_typed_table(table_s1, tmp_path):
    rows, _ = table_s1.build(table_s1.Source(ROOT / "src" / "foldmap"))
    tex = table_s1.longtable(rows, "abc1234")
    assert table_s1.typed_rows(tex) == [(r.split(" & ")[0], r.split(" & ")[2].removesuffix(r" \\")) for r in rows]
    typed = tex.replace("200,000; 1,600,000", "100,000; 1,600,000")
    (tmp_path / "supplement.tex").write_text(typed, encoding="utf-8")
    code = table_s1.main(["--src", str(ROOT / "src" / "foldmap"), "--out", str(tmp_path / "t"), "--check-against",
                          str(tmp_path / "supplement.tex")])  # fmt: skip
    assert code == 1
    assert json.loads((tmp_path / "t" / "params.json").read_text(encoding="utf-8"))["values"]


# ---------------------------------------------------------------------------------------------- the NA report


def test_the_na_report_is_read_in_both_formats():
    common = _load("common")
    rows = [{"entry": "1EHZ", "failed": []}]
    assert common.na_rows(rows) == rows
    assert common.na_rows({"source_stamp": {}, "entries": rows}) == rows


# ---------------------------------------------------------------------------------------------- the fourth sample


@pytest.fixture(scope="module")
def fourth():
    return _load("fourth_sample")


def test_the_fourth_samples_choices_are_fixed(fourth):
    assert (fourth.SEED, fourth.DATE) == (20261025, "2026-10-13")
    assert fourth.OUT == ROOT / "docs" / "validation" / "random100-fourth" / "sample.json"
    assert len(fourth.SAMPLE_SHA256) == 3 and (ROOT / fourth.PROTOCOL).name == "fourth-sample-protocol.md"


@needs_samples
def test_samples_one_to_three_are_as_the_protocol_pinned_them(fourth):
    for rel, digest in fourth.SAMPLE_SHA256.items():
        assert fourth.sha256_text(ROOT / rel) == digest, rel
    seeds = [json.loads((ROOT / rel).read_text(encoding="utf-8"))["seed"] for rel in fourth.SAMPLE_SHA256]
    assert fourth.SEED not in seeds


@needs_samples
def test_the_draw_excludes_samples_one_to_three_and_records_its_protocol(fourth, tmp_path, monkeypatch):
    r100 = fourth.random100()
    want, _ = fourth.excluded(r100)
    assert len(want) >= 300
    asked = {}

    def fake_sample(seed, date, exclude=()):
        asked.update(seed=seed, date=date, exclude=list(exclude))
        entries = [{"id": f"9Z{k:02d}"} for k in range(100)]
        return {"seed": seed, "query_date": date, "entries": entries, "replaced": [], "quotas": {}}

    monkeypatch.setattr(r100, "sample", fake_sample)
    monkeypatch.setattr(fourth, "random100", lambda: r100)
    out = tmp_path / "random100-fourth" / "sample.json"
    assert fourth.main(["draw", "--release", "v9.9.9", "--out", str(out)]) == 0
    assert asked == {"seed": 20261025, "date": "2026-10-13", "exclude": want}
    doc = json.loads(out.read_text(encoding="utf-8"))
    p = doc["protocol"]
    assert (p["seed"], p["query_date"], p["release"], p["excluded_ids"]) == (
        20261025,
        "2026-10-13",
        "v9.9.9",
        len(want),
    )
    assert [f["path"] for f in p["excluded"]] == list(fourth.SAMPLE_SHA256)
    assert p["source_stamp"]["code_fingerprint"] and p["document"] == "docs/validation/fourth-sample-protocol.md"
    with pytest.raises(SystemExit, match="drawn once"):
        fourth.main(["draw", "--release", "v9.9.9", "--out", str(out)])


@needs_samples
def test_a_draw_holding_an_excluded_id_is_refused(fourth, tmp_path, monkeypatch):
    r100 = fourth.random100()
    want, _ = fourth.excluded(r100)
    monkeypatch.setattr(r100, "sample", lambda seed, date, exclude=(): {"entries": [{"id": want[0]}], "quotas": {}})
    monkeypatch.setattr(fourth, "random100", lambda: r100)
    with pytest.raises(SystemExit, match="samples 1-3"):
        fourth.main(["draw", "--release", "v9.9.9", "--out", str(tmp_path / "s.json")])
    assert not (tmp_path / "s.json").exists()


def test_a_changed_sample_stops_the_draw(fourth, tmp_path):
    for rel in fourth.SAMPLE_SHA256:
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text('{"entries": [], "replaced": []}\n', encoding="utf-8")
    with pytest.raises(SystemExit, match="changed after the protocol"):
        fourth.excluded(None, root=tmp_path)
