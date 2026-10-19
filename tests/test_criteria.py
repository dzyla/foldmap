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
    assert "S2 barrel closure pairs: not judged (no closed barrel drawn)" in out  # nothing to judge: never a pass
    assert out.index("S6 accessory drawn at its partner: not judged (no accessory strand)") < out.index(
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
    off = _font_off(result)
    assert bool(off) == bool(figcheck.ink_font_reason(drawn.figure))  # none where text is in Liberation Sans
    alike = [k for k in FIGURE_JUDGED if k not in off]
    assert {k: script[k] for k in alike} == {k: {x: v for x, v in package[k].items() if x != "here"} for k in alike}
    for k in set(FIGURE_JUDGED) & set(off):  # the campaign fails it closed; `foldmap check` names it not measured
        assert script[k]["pass"] is False and script[k]["detail"].startswith("default not measured ("), k
        assert package[k]["pass"] is None and package[k]["here"] == off[k], k
    assert all(script[k]["pass"] is not None for k in FIGURE_JUDGED if k not in ("A15", "A20"))
    # nothing to judge is no pass (final review B, Important 3): ubiquitin has no copies and no ligand
    want = {"1UBQ.cif": ("does not apply", "does not apply"), "2HHB.cif": ("pass", "not judged")}[name]
    assert (criteria.outcome(script["A15"]), criteria.outcome(script["A20"])) == want


def _font_off(result: dict) -> dict[str, str]:
    """The criteria whose ink checks the text font kept from being measured in this drawing (criteria.INK_READ), and
    the reason figcheck gave. The ink checks are calibrated on Liberation Sans (figcheck.INK_FONT); where text resolves
    to another font (Arial on macOS and Windows: render._FONTS lists it first) the campaign fails these closed and
    `foldmap check` names them not measured here, so there the two say different things by design. None on Linux."""
    skipped = result["checks"].get("skipped") or {}
    return {k: str(skipped[c]) for k, reads in criteria.INK_READ.items() for c in reads if c in skipped}


OFF_FONT = "the text resolved to DejaVu Sans, not Liberation Sans"


def _off_font(real):
    def off_font(*args, **kw):
        got = real(*args, **kw)
        return {
            **got,
            "text_collisions": None,
            "off_page": None,
            "domain_name_collisions": None,
            "skipped": dict.fromkeys(figcheck.INK_CHECKS, OFF_FONT),
        }

    return off_font


def test_off_font_a16_still_judges_the_domain_rules_and_panels(monkeypatch):
    """T13a review Important 1: off the calibration font only A16's domain names go unmeasured; its domain rules and
    panel overlaps measure no ink and are judged (a cross-chain domain fails), and a figure without domains reports A16
    as not applying, not as unmeasured."""
    real_checks, real_domains = figcheck.figure_checks, figcheck.domain_checks
    monkeypatch.setattr(figcheck, "figure_checks", _off_font(real_checks))
    spec = figure_spec(DATA / "1UBQ.cif")
    got = criteria.check_figure(draw_spec(spec), spec, DATA / "1UBQ.cif")["criteria"]
    assert got["A16"] == {"pass": None, "detail": "", "status": "does not apply"}  # no domains

    spec = figure_spec(DATA / "2LZM.cif", domains="auto")  # two automatic domains
    drawn = draw_spec(spec)
    got = criteria.check_figure(drawn, spec, DATA / "2LZM.cif")["criteria"]
    assert got["A16"]["pass"] is None and OFF_FONT in got["A16"]["here"] and "names" in got["A16"]["here"]

    def crossing(d, collisions=None):
        return {**real_domains(d, collisions), "cross_chain": ["D1 holds A:10-20 and B:30-40"]}

    monkeypatch.setattr(figcheck, "domain_checks", crossing)
    got = criteria.check_figure(drawn, spec, DATA / "2LZM.cif")["criteria"]
    assert got["A16"]["pass"] is False and "domains 1" in got["A16"]["detail"]
    assert OFF_FONT in got["A16"]["detail"]  # the names part is said to be unmeasured


@pytest.mark.parametrize(
    "name, options, view",
    [("1C3W.cif", {"membrane": "auto", "assembly": "asu"}, "membrane"), ("2HHB.cif", {"domains": "auto"}, "domains")],
)
def test_a_membrane_or_domains_figure_is_judged_as_the_campaign_judges_that_view(name, options, view):
    """T13a review Important 2: a --membrane or --domains figure is judged by the criteria the campaign applies to its
    membrane or domains view (figure_criteria under that view's name), value and detail; the criteria the campaign
    judges on the default view only (A8, A11, A15 ...) are not measured here, with the reason."""
    path = DATA / name
    spec = figure_spec(path, **options)
    drawn = draw_spec(spec)
    result = criteria.measure_figure(view, drawn, spec, path)
    campaign = criteria.figure_criteria({view: {"result": result}})
    got = criteria.check_figure(drawn, spec, path)["criteria"]
    judged = {k for k, v in campaign.items() if v["pass"] is not None}
    assert {"A3", "A5", "A6", f"A1{6 if view == 'domains' else 8}"} <= judged
    off = _font_off(result)  # none where text is in Liberation Sans
    assert bool(off) == bool(figcheck.ink_font_reason(drawn.figure))
    for k in judged - set(off):
        assert {"pass": got[k]["pass"], "detail": got[k]["detail"]} == campaign[k], k
    for k in judged & set(off):  # the campaign fails it closed; `foldmap check` says the font kept it unmeasured
        assert campaign[k]["pass"] is False and f"{view} not measured (" in campaign[k]["detail"], k
        assert got[k]["pass"] is not True and off[k] in (got[k].get("here") or got[k]["detail"]), (k, got[k])
        assert k == "A16" or got[k]["pass"] is None  # A16 still judges its domain rules and panels (above)
    for k in ("A8", "A11", "A15", "A20"):
        assert got[k]["pass"] is None and f"{view} view" in got[k]["here"], (k, got[k])


def test_an_unmeasured_criterion_never_prints_as_a_pass():
    """T13a review Minor 6 ruling: in `foldmap check` a criterion that was not measured prints as not measured, never
    as a pass (the campaign's A10 reads pass with "unmeasured: <reason>")."""
    lines = criteria.criteria_lines({"A10": {"pass": True, "detail": "unmeasured: no mkdssp on PATH"}})
    assert "A10 secondary structure agrees with mkdssp: not measured here (no mkdssp on PATH)" in lines
    assert not any(line.endswith(": pass") for line in lines if line.startswith("A10"))


def test_a_criterion_with_nothing_to_judge_prints_not_judged():
    """Final review B, Important 3: `foldmap check` prints a criterion the figure gave nothing to judge as not judged,
    with why, never as a pass."""
    lines = criteria.criteria_lines(
        {"A15": {"pass": None, "detail": "default B~A n/a: 2 elements", "status": "not judged"}}
    )
    assert "A15 copies drawn alike: not judged (default B~A n/a: 2 elements)" in lines
    assert "A16 domains: one chain, alike, not tiny, panels apart, names clear: does not apply" in lines
