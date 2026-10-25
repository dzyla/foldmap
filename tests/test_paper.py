"""The paper pipeline in scripts/paper (plan 28 Task 1): Table S1 is read from the code's constants (every value
found, by name), the nucleic-acid report is read in both its formats, and the fourth random sample's sampler keeps
its protocol: a fixed seed and cutoff, no id of samples 1-3, drawn once."""

import importlib.util
import json
import shutil
import subprocess
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

    out = tmp_path / "random100-fourth" / "sample.json"
    monkeypatch.setattr(r100, "sample", fake_sample)
    monkeypatch.setattr(fourth, "random100", lambda: r100)
    monkeypatch.setattr(fourth, "OUT", out)  # its one path, here in a temporary folder
    monkeypatch.setattr(fourth, "release_problem", lambda tag, stamp: None)  # as at the release
    assert fourth.main(["draw", "--release", "v9.9.9"]) == 0
    assert asked == {"seed": 20261025, "date": "2026-10-13", "exclude": want}
    p = json.loads(out.read_text(encoding="utf-8"))["protocol"]
    assert (p["seed"], p["query_date"], p["release"], p["excluded_ids"]) == (20261025, "2026-10-13", "v9.9.9",
                                                                             len(want))  # fmt: skip
    assert [f["path"] for f in p["excluded"]] == list(fourth.SAMPLE_SHA256)
    assert p["source_stamp"]["code_fingerprint"] and p["document"] == "docs/validation/fourth-sample-protocol.md"
    with pytest.raises(SystemExit, match="drawn once"):
        fourth.main(["draw", "--release", "v9.9.9"])
    with pytest.raises(SystemExit):  # no --out: it writes only to its fixed path (ruling Q14)
        fourth.main(["draw", "--release", "v9.9.9", "--out", str(tmp_path / "elsewhere.json")])


@needs_samples
def test_a_draw_holding_an_excluded_id_is_refused(fourth, tmp_path, monkeypatch):
    """Review M10: what random100 run will draw, the entries and the A1c refusals, is checked."""
    r100 = fourth.random100()
    want, _ = fourth.excluded(r100)
    refused = {"entries": [{"id": "9ZZZ"}], "replaced": [{"id": want[0], "replacement": None}], "quotas": {}}
    monkeypatch.setattr(r100, "sample", lambda seed, date, exclude=(): refused)
    monkeypatch.setattr(fourth, "random100", lambda: r100)
    monkeypatch.setattr(fourth, "OUT", tmp_path / "s.json")
    monkeypatch.setattr(fourth, "release_problem", lambda tag, stamp: None)
    with pytest.raises(SystemExit, match="samples 1-3"):
        fourth.main(["draw", "--release", "v9.9.9"])
    assert not (tmp_path / "s.json").exists()


def test_the_draw_runs_only_in_a_clean_checkout_at_the_release_tag(fourth, tmp_path):
    """Ruling Q14: HEAD must be the tag's commit, the tree clean, the foldmap the checkout's; no .git: refused."""
    if shutil.which("git") is None:
        pytest.skip("needs git")

    def git(*args):
        return subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t", *args],
                              capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip

    git("init", "-q")
    (tmp_path / "a.txt").write_text("1\n", encoding="utf-8")
    git("add", ".")
    git("commit", "-q", "-m", "one")
    git("tag", "v1.0.0")
    head = git("rev-parse", "HEAD")
    clean = {"git": {"commit": head, "dirty": False, "source": "git"}, "installed_matches_checkout": True}
    assert fourth.release_problem("v1.0.0", clean, tmp_path) is None
    assert "not v1.0.0" in fourth.release_problem("v1.0.0", clean | {"git": clean["git"] | {"commit": "0" * 40}},
                                                   tmp_path)  # fmt: skip
    assert "no-such-tag" in fourth.release_problem("no-such-tag", clean, tmp_path)
    dirty = clean | {"git": clean["git"] | {"dirty": True, "dirty_paths": ["src/x.py"]}}
    assert "clean tree" in fourth.release_problem("v1.0.0", dirty, tmp_path)
    asserted = clean | {"git": clean["git"] | {"source": "asserted by runner"}}
    assert "no git checkout" in fourth.release_problem("v1.0.0", asserted, tmp_path)
    assert "no git checkout" in fourth.release_problem("v1.0.0", {"git": None}, tmp_path)
    other = clean | {"installed_matches_checkout": False}
    assert "not the checkout" in fourth.release_problem("v1.0.0", other, tmp_path)


def test_a_changed_sample_stops_the_draw(fourth, tmp_path):
    for rel in fourth.SAMPLE_SHA256:
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text('{"entries": [], "replaced": []}\n', encoding="utf-8")
    with pytest.raises(SystemExit, match="changed after the protocol"):
        fourth.excluded(None, root=tmp_path)


# ---------------------------------------------------------------------------------------------- stamps, checked


def _repo_with_foldmap(root: Path):
    """A git checkout holding src/foldmap/a.py and scripts/c.py, committed; returns (git, package, script, commit)."""

    def git(*args):
        return subprocess.run(["git", "-C", str(root), "-c", "user.name=t", "-c", "user.email=t@t", *args],
                              capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip

    (root / "pyproject.toml").write_text("[project]\nname = 'x'\n", encoding="utf-8")
    pkg = root / "src" / "foldmap"
    pkg.mkdir(parents=True)
    (pkg / "a.py").write_text("X = 1\n", encoding="utf-8")
    script = root / "scripts" / "c.py"
    script.parent.mkdir()
    script.write_text("print(1)\n", encoding="utf-8")
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "one")
    return git, pkg, script, git("rev-parse", "HEAD")


needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="needs git")


@needs_git
def test_a_stamp_is_verified_against_the_commit_it_names(tmp_path, monkeypatch):
    """Ruling Q11: the fingerprint of the commit's src/foldmap (git archive), by the same function, must match."""
    from foldmap.source_stamp import source_stamp

    verify = _load("verify_stamp")
    git, pkg, script, head = _repo_with_foldmap(tmp_path)
    st = source_stamp(script, package=pkg)
    assert st["git"]["commit"] == head
    good = verify.verify_report({"source_stamp": st, "entries": [{"source_stamp": st}]}, repo=tmp_path)
    assert [r["verified"] for r in good] == [True]  # one distinct code, checked once
    forged = st | {"code_fingerprint": "0" * 16}
    assert "is not" in verify.verify_report({"source_stamp": forged}, repo=tmp_path)[0]["why"]
    other_script = st | {"script": {"name": "scripts/c.py", "sha256": "0" * 64}}
    assert "differs" in verify.verify_report({"source_stamp": other_script}, repo=tmp_path)[0]["why"]
    asserted = st | {"git": {"commit": "f" * 40, "dirty": False, "source": "asserted by runner"}}
    assert "not a commit" in verify.verify_report({"source_stamp": asserted}, repo=tmp_path)[0]["why"]
    nothing = st | {"git": None}
    assert "names no commit" in verify.verify_report({"source_stamp": nothing}, repo=tmp_path)[0]["why"]
    # ruling Q16: a dirty tree, or one of unknown state, never verifies, even when src/foldmap and the script match
    sibling = tmp_path / "scripts" / "na_validate.py"  # a script the campaign loads, outside the code key
    sibling.write_text("print(2)\n", encoding="utf-8")
    git("add", ".")
    git("commit", "-q", "-m", "sibling")
    head2 = git("rev-parse", "HEAD")
    sibling.write_text("print(3)\n", encoding="utf-8")  # an uncommitted edit
    dirty = source_stamp(script, package=pkg)
    assert dirty["git"]["dirty"] is True and dirty["git"]["dirty_paths"] == ["scripts/na_validate.py"]
    got = verify.verify_report({"source_stamp": dirty}, repo=tmp_path)[0]
    assert got["verified"] is False and "dirty" in got["why"] and got["commit"] == head2
    unknown = st | {"git": {"commit": head, "dirty": None, "source": "asserted by runner"}}
    got = verify.verify_report({"source_stamp": unknown}, repo=tmp_path)[0]
    assert got["verified"] is False and "unknown" in got["why"]
    asserted_clean = st | {"git": {"commit": head, "dirty": False, "source": "asserted by runner"}}
    assert verify.verify_report({"source_stamp": asserted_clean}, repo=tmp_path)[0]["verified"] is True
    # a paper step's run stamp: checked against the harness commit its tree was exported at
    run = {"source_stamp": {"harness_commit": head, "source_stamp": st | {"git": None}}}
    assert verify.verify_report(run, repo=tmp_path)[0]["verified"] is True
    out = tmp_path / "report.json"
    out.write_text(json.dumps({"source_stamp": st}), encoding="utf-8")
    monkeypatch.setattr(verify, "REPO", tmp_path)
    assert verify.verify_report(verify.load(out), repo=tmp_path)[0]["verified"] is True


@needs_git
def test_extract_reads_commits_from_stamps_only_and_refuses_mixed_or_unquotable_reports(tmp_path):
    """Review I1 (5): both NA formats, legacy files, verified stamps, mixed code and not-for-paper."""
    from foldmap.source_stamp import code_mix, source_stamp

    extract = _load("extract")
    rec, how, stop = extract.recorded_commits("na_report", [{"entry": "1EHZ", "failed": []}])
    assert rec == {"report": None} and "no source stamp" in how and not stop
    rec, how, stop = extract.recorded_commits("r100_third_report", {"foldmap_commit": "7696947", "entries": []})
    assert rec == {"report": None, "legacy_foldmap_commit": "7696947"} and not stop  # asserted after the run: unused
    git, pkg, script, head = _repo_with_foldmap(tmp_path)
    st = source_stamp(script, package=pkg)
    doc = {"source_stamp": st, **code_mix([st]), "entries": [{"entry": "1EHZ", "source_stamp": st}]}
    rec, _, stop = extract.recorded_commits("na_report", doc, repo=tmp_path)
    assert rec["report"] == head and rec["verified"][0]["verified"] and not stop
    mixed = doc | code_mix([st, st | {"code_fingerprint": "0" * 16}])
    assert any("mixed_code" in p for p in extract.recorded_commits("na_report", mixed, repo=tmp_path)[2])
    marked = doc | {"not_for_paper": "run with --not-for-paper (asked)"}
    assert any("not for the paper" in p for p in extract.recorded_commits("na_report", marked, repo=tmp_path)[2])
    forged = doc | {"source_stamp": st | {"code_fingerprint": "0" * 16}}
    assert any("not verified" in p for p in extract.recorded_commits("na_report", forged, repo=tmp_path)[2])
    assert rec["dirty"] is False and rec["source"] == "git" and rec["commits"] == [head]
    # ruling Q17: verified stamps at two commits (scripts/ is outside the code key) are a mixed report
    git("commit", "-q", "--allow-empty", "-m", "two")
    later = st | {"git": st["git"] | {"commit": git("rev-parse", "HEAD")}}
    two = doc | {"entries": [{"entry": "1EHZ", "source_stamp": st}, {"entry": "1FFK", "source_stamp": later}]}
    rec2, _, stop2 = extract.recorded_commits("na_report", two, repo=tmp_path)
    assert len(rec2["commits"]) == 2 and any("2 commits" in p for p in stop2)


def test_make_numbers_refuses_a_mixed_or_unquotable_report():
    common = _load("common")
    common.quotable("ok", {"mixed_code": False})
    common.quotable("legacy", [{"entry": "1EHZ"}])
    with pytest.raises(SystemExit, match="different code"):
        common.quotable("x", {"mixed_code": True, "entry_code": {"a": 1, "b": 1}})
    with pytest.raises(SystemExit, match="not for the paper"):
        common.quotable("x", {"not_for_paper": "asked"})
    assert common.stamped_commit([]) == "no source stamp (written before plan 28 Task 1)"
    st = {"code_fingerprint": "f", "git": {"commit": "a" * 40, "dirty": False, "source": "asserted by runner"}}
    assert common.stamped_commit({"source_stamp": st}).startswith("at aaaaaaa (asserted by runner, dirty False")


def test_a_tree_older_than_the_stamp_is_still_fingerprinted():
    """Review M9: run_stamp's fallback computes the tree's fingerprint (the same algorithm) and the versions."""
    from foldmap import app_jobs

    common = _load("common")
    aj = common.repo_module("app_jobs")
    here = (ROOT / "src" / "foldmap").resolve()
    assert aj.code_fingerprint(here) == app_jobs.code_fingerprint(here)
    assert common.repo_module("source_stamp").RUNTIME[0] == "gemmi"
