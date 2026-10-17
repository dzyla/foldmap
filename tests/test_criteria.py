"""The figure criteria A1-A20 in the package (foldmap.criteria, plan 10 Task 13a): the random-100 script judges by
them, and `foldmap check` prints them under the science checks for the one figure it draws."""

import importlib.util
from pathlib import Path

import pytest

from foldmap import criteria, figcheck
from foldmap.cli import main
from foldmap.spec import draw_spec, figure_spec

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"
FIGURE_JUDGED = ("A3", "A4", "A5", "A6", "A7", "A8", "A9", "A11", "A15", "A20")


@pytest.fixture(scope="module")
def r100():
    spec = importlib.util.spec_from_file_location("random100", ROOT / "scripts" / "random100.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_check_prints_the_figure_criteria_under_the_science_checks(capsys):
    assert main(["check", str(DATA / "1UBQ.cif")]) == 0
    out = capsys.readouterr().out.splitlines()
    assert "S1 sheet neighbours pair: pass" in out
    assert "A5 no element or stem boxes overlapping: pass" in out and "A11 N->C and sheet order: pass" in out
    assert out.index("S6 accessory drawn at its partner: pass") < out.index(
        "A1 no crash: exit 0, files written, "
        "no traceback, no Python warning: "
        f"not measured here ({criteria.NOT_HERE['A1']})"
    )
    assert f"A10 secondary structure agrees with mkdssp: not measured here ({criteria.NOT_HERE['A10']})" in out
    assert "A17 nucleic-acid metrics and stem names by their letters: does not apply" in out  # no nucleic acid
    assert "A18 membrane: TM helices in the band, loops out of the slab: does not apply" in out  # no --membrane
    assert [line.split()[0] for line in out if line[:1] in "SA"] == [
        "S1",
        "S2",
        "S3",
        "S4",
        "S5",
        "S6",
        *criteria.CRITERIA,
    ]


def test_check_exits_1_when_a_figure_criterion_fails(capsys, monkeypatch):
    monkeypatch.setattr(figcheck, "element_overlaps", lambda lay: [("A:2-7", "A:12-17")])
    assert main(["check", str(DATA / "1UBQ.cif")]) == 1
    out = capsys.readouterr().out.splitlines()
    assert "A5 no element or stem boxes overlapping: fail" in out
    assert "  default 1 (['A:2-7', 'A:12-17'])" in out  # as the campaign's report details it
    assert "S1 sheet neighbours pair: pass" in out  # the science checks are still printed, and pass


def test_a_check_the_font_keeps_from_measuring_is_not_measured_here_not_failed(monkeypatch):
    """figcheck.figure_checks does not measure the ink checks off its calibration font (counts None, "skipped"); the
    campaign fails such a view closed, `foldmap check` names it not measured here, with the reason."""
    real = figcheck.figure_checks

    def off_font(*args, **kw):
        got = real(*args, **kw)
        why = "the text resolved to DejaVu Sans, not Liberation Sans"
        return {**got, "text_collisions": None, "off_page": None, "skipped": dict.fromkeys(figcheck.INK_CHECKS, why)}

    monkeypatch.setattr(figcheck, "figure_checks", off_font)
    spec = figure_spec(DATA / "1UBQ.cif")
    got = criteria.check_figure(draw_spec(spec), spec, DATA / "1UBQ.cif")["criteria"]
    for key in ("A6", "A7"):
        assert got[key]["pass"] is None and got[key]["here"] == "the text resolved to DejaVu Sans, not Liberation Sans"
    assert got["A5"]["pass"] is True


@pytest.mark.parametrize("name", ["1UBQ.cif", "2HHB.cif"])
def test_the_script_and_the_package_judge_one_figure_alike(r100, tmp_path, name):
    """The script's criteria are the package's, and the figure criteria `foldmap check` judges one figure by are those
    the campaign gives its default view, value and detail."""
    path = DATA / name
    spec = figure_spec(path)
    drawn = draw_spec(spec)
    result = r100.measure_view("default", drawn, spec, path, name, tmp_path)
    view = {"exit": 0, "seconds": 1.0, "timed_out": False, "traceback": False, "written": {}, "result": result}
    e = {"views": {"default": view}}
    assert r100.criteria is criteria.criteria and r100.criteria(e) == criteria.criteria(e)
    script = r100.criteria(e)
    package = criteria.check_figure(drawn, spec, path)["criteria"]
    assert {k: script[k] for k in FIGURE_JUDGED} == {
        k: {"pass": package[k]["pass"], "detail": package[k]["detail"]} for k in FIGURE_JUDGED
    }
    assert all(script[k]["pass"] is not None for k in FIGURE_JUDGED)
