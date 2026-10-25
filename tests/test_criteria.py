"""The figure criteria A1-A20 in the package (foldmap.criteria, plan 10 Task 13a): the random-100 script judges by
them, and `foldmap check` prints them under the science checks for the one figure it draws."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from foldmap import criteria, figcheck
from foldmap.align import Mates
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
    a18 = (
        "A18 membrane: only the band set in the band (TM helices across it, slab helices at their depths, crossing "
        "strands), no element or loop nearer than its residues or on the wrong side"
    )
    assert f"{a18}: does not apply" in out  # no --membrane
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


# ------------------------------------------------------------------------- drawn like a reference: H1-H5 (13a-9)


def _stub(points: dict, sheets=()):
    """A layout stub as the H checks read it: placed {id: (cx, cy, label[, kind, angle])}, residue order by insertion,
    and its sheet blocks."""
    placed = {}
    for n, (k, (x, y, label, *rest)) in enumerate(points.items()):
        kind, angle = [*rest, "H", 0.0][:2]
        sse = SimpleNamespace(start=10 * n, chain="A", kind=kind)
        placed[k] = SimpleNamespace(cx=x, cy=y, label=label, ghost=False, sse=sse, angle=angle)
    return SimpleNamespace(placed=placed, sheets=[list(s) for s in sheets])


def _mates(pairs: dict, same_fold=True, extras=(), only_a=()):
    smaller = len(pairs) + min(len(only_a), len(extras))
    return Mates(
        mate=dict(pairs), pieces={}, extras=tuple(extras), only_a=tuple(only_a), same_fold=same_fold, smaller=smaller
    )


def test_check_like_itself_passes_every_h_check(capsys):
    """1UBQ drawn like itself: H1 one fold and drawn alike, H2 0 (median and 90th percentile), H3 every pair agrees,
    H4 every name the same, H5 nothing lost; drawn independently it is the same figure, so the baseline is 0 too."""
    assert main(["check", str(DATA / "1UBQ.cif"), "--like", str(DATA / "1UBQ.cif")]) == 0
    out = capsys.readouterr().out.splitlines()
    assert [line.split()[0] for line in out if line[:1] == "H"] == list(criteria.LIKE_CHECKS)
    for key, what in criteria.LIKE_CHECKS.items():
        assert f"{key} {what}: pass" in out, key
    detail = {out[k].split()[0]: out[k + 1] for k, line in enumerate(out) if line[:1] == "H"}
    assert detail["H1"].startswith("  one fold, drawn alike: TM-score 1.00")
    assert detail["H2"].startswith("  median 0.000, 90th percentile 0.000 page units over ")
    assert "drawn independently: median 0.000, 90th percentile 0.000" in detail["H2"]
    assert detail["H3"].startswith("  1.0000 of ") and "provisional" not in detail["H3"]
    assert detail["H4"].startswith("  names compared on all ") and "differ" not in detail["H4"]
    assert " of the " in detail["H5"] and "lost" not in detail["H5"]


def test_check_without_like_prints_no_h_check(capsys):
    assert main(["check", str(DATA / "1UBQ.cif")]) == 0
    assert not any(line[:1] == "H" for line in capsys.readouterr().out.splitlines())


def test_check_like_a_refused_pair(capsys):
    """2LZM is not ubiquitin's fold: H1 passes (refused with the note), H2-H4 are not judged, never passed; H5 still
    compares the figure with the note against the figure alone."""
    main(["check", str(DATA / "2LZM.cif"), "--like", str(DATA / "1UBQ.cif")])
    out = capsys.readouterr().out.splitlines()
    assert "H1 drawn alike exactly when one fold: pass" in out
    assert out[out.index("H1 drawn alike exactly when one fold: pass") + 1].startswith(
        "  not one fold, refused with the note: TM-score 0.33"
    )
    for key in ("H2", "H3", "H4"):
        line = next(x for x in out if x.startswith(key))
        assert f": {criteria.NOT_JUDGED} (not drawn like 1UBQ.cif: not the same fold" in line, line
    assert any(x.startswith("H5 still a good figure: ") for x in out)


def test_h2_is_zero_after_a_known_shift_and_turn():
    pts = {"a1": (0.0, 0.0, "α1"), "a2": (3.0, 1.0, "β1"), "a3": (-2.0, 4.0, "β2"), "a4": (5.0, -3.0, "α2")}
    turn = np.array([[np.cos(0.7), -np.sin(0.7)], [np.sin(0.7), np.cos(0.7)]])
    moved = {f"b{k[1]}": (*(turn @ np.array(p[:2]) + [12.5, -7.25]), p[2]) for k, p in pts.items()}
    lay_a, lay_b = _stub(pts), _stub(moved)
    pairs = criteria.like_pairs(lay_a, lay_b, _mates({f"b{k[1]}": k for k in pts}))
    got = criteria.h2_places(lay_a, lay_b, pairs)
    assert got["pass"] is True and got["value"]["median"] == 0.0 and got["value"]["p90"] == 0.0
    assert got["value"]["mates"] == 4
    shifted = _stub({k: (p[0] + 0.5, p[1], p[2]) if k == "b2" else p for k, p in moved.items()})
    assert criteria.h2_places(lay_a, shifted, pairs)["value"]["p90"] > 0


def test_h2_h3_not_judged_under_three_mates_h4_under_one():
    pts = {"a1": (0.0, 0.0, "α1"), "a2": (3.0, 1.0, "β1")}
    lay = _stub(pts)
    two = criteria.like_checks(lay, lay, SimpleNamespace(tm_score=1.0), _mates({"a1": "a1", "a2": "a2"}), "", True)
    for key in ("H2", "H3"):
        assert criteria.outcome(two[key]) == criteria.NOT_JUDGED and "2 mates drawn" in two[key]["detail"], key
    assert two["H4"]["pass"] is True
    none = criteria.like_checks(lay, lay, SimpleNamespace(tm_score=1.0), _mates({}), "", True)
    assert criteria.outcome(none["H4"]) == criteria.NOT_JUDGED
    assert criteria.outcome(none["H5"]) == criteria.NOT_JUDGED  # B not drawn alone: nothing to compare


def test_h1_a_forced_negative_fails():
    aln = SimpleNamespace(tm_score=0.31)
    no = _mates({"b1": "a1"}, same_fold=False, extras=("b2", "b3"), only_a=("a2", "a3"))
    assert criteria.h1_matched(aln, no, False, True, "negative")["pass"] is True
    assert criteria.h1_matched(aln, no, False, True)["pass"] is True  # check: refused as it should be
    assert criteria.h1_matched(aln, no, False, False, "negative")["pass"] is False  # on its own, but not refused
    forced = criteria.h1_matched(aln, no, True, False, "negative")
    assert forced["pass"] is False and "forced" in forced["detail"]
    assert criteria.h1_matched(aln, no, True, True, "negative")["pass"] is False  # drawn alike is never refused
    yes = _mates({"b1": "a1", "b2": "a2"}, extras=("b3",), only_a=("a3",))
    assert criteria.h1_matched(SimpleNamespace(tm_score=0.6), yes, True, False, "positive")["pass"]
    assert criteria.h1_matched(aln, no, False, True, "positive")["pass"] is False  # a positive refused
    assert "1 mates of the 3 elements" in criteria.h1_matched(aln, no, False, True)["detail"]  # Mates.smaller
    unaligned = criteria.h1_matched(None, None, False, False, note="not drawn like A.cif: ...")
    assert criteria.outcome(unaligned) == criteria.NOT_JUDGED and "not drawn like A.cif" in unaligned["detail"]


def test_h3_a_turned_map_and_a_sheet_in_another_order():
    pts = {
        "a1": (0.0, 0.0, "β1", "E", np.pi / 2),
        "a2": (1.5, 0.0, "β2", "E", -np.pi / 2),
        "a3": (3.0, 0.0, "β3", "E", np.pi / 2),
        "a4": (1.5, 5.0, "α1", "H", 0.0),
    }
    mates = _mates({k: k for k in pts})
    lay_a = _stub(pts, sheets=[["a1", "a2", "a3"]])
    pairs = criteria.like_pairs(lay_a, lay_a, mates)
    same = criteria.h3_relations(lay_a, lay_a, pairs)
    assert same["pass"] is True and same["value"]["share"] == 1.0 and same["value"]["sheets"] == []
    mirrored = _stub({k: (-x, y, *rest) for k, (x, y, *rest) in pts.items()}, sheets=[["a1", "a2", "a3"]])
    got = criteria.h3_relations(lay_a, mirrored, pairs)
    assert got["pass"] is False and got["value"]["share"] < 0.95
    swapped = dict(pts, a1=(1.5, 0.0, "β1", "E", np.pi / 2), a2=(0.0, 0.0, "β2", "E", -np.pi / 2))
    lay_b = _stub(swapped, sheets=[["a2", "a1", "a3"]])
    got = criteria.h3_relations(lay_a, lay_b, pairs)
    assert got["pass"] is False and "sheet order β2, β1, β3 for β1, β2, β3" in got["detail"]


def test_h4_names_by_a19s_rule():
    lay_a = _stub({"a1": (0, 0, "β1"), "a2": (2, 0, "α1"), "a3": (4, 0, "β2")})
    alike = _stub({"b1": (0, 0, "β1"), "bx": (1, 0, "β2"), "b2": (2, 0, "α1"), "b3": (4, 0, "β3")})
    mates = _mates({"b1": "a1", "b2": "a2", "b3": "a3"}, extras=("bx",))
    assert criteria.h4_names(lay_a, alike, mates)["pass"] is True  # β3 less the extra strand before it is β2
    other = _stub({"b1": (0, 0, "β1"), "b2": (2, 0, "α2"), "b3": (4, 0, "β2")})
    got = criteria.h4_names(lay_a, other, _mates({"b1": "a1", "b2": "a2", "b3": "a3"}))
    assert got["pass"] is False and "α2~α1" in got["detail"]
    crossed = _stub({"b1": (0, 0, "β1"), "b2": (2, 0, "α1"), "b3": (4, 0, "β2")})
    lay_eta = _stub({"a1": (0, 0, "β1"), "a2": (2, 0, "η1"), "a3": (4, 0, "β2")})
    note = "drawn like A.cif chain A: 3 of 3 elements matched; named by their own helix type: α1 for η1"
    assert criteria.h4_names(lay_eta, crossed, _mates({"b1": "a1", "b2": "a2", "b3": "a3"}), note)["pass"] is True
    assert criteria.h4_names(lay_eta, crossed, _mates({"b1": "a1", "b2": "a2", "b3": "a3"}))["pass"] is False


def test_h5_fails_when_the_like_figure_fails_a_check_b_alone_passes(monkeypatch):
    spec = figure_spec(DATA / "1UBQ.cif")
    alone = criteria.check_figure(draw_spec(spec), spec, DATA / "1UBQ.cif")
    assert criteria.h5_good_figure(alone, alone)["pass"] is True
    assert "checks B alone passes" in criteria.h5_good_figure(alone, alone)["detail"]
    assert "checks 1UBQ.cif chain A alone passes" in criteria.h5_good_figure(alone, alone, "1UBQ.cif chain A")["detail"]
    monkeypatch.setattr(figcheck, "element_overlaps", lambda lay: [("A:2-7", "A:12-17")])
    broken = criteria.check_figure(draw_spec(spec), spec, DATA / "1UBQ.cif")
    got = criteria.h5_good_figure(broken, alone)
    assert got["pass"] is False and got["value"]["lost"] == ["A5"] and "A5 (fail" in got["detail"]
    assert criteria.h5_good_figure(alone, broken)["pass"] is True  # only what B alone passes must hold


def test_a_helix_drawn_in_parts_is_judged_as_its_run():
    """A helix B splits where A does not (Mates.pieces) is laid end to end along its mate, the first part in
    sequence named as A's helix and the next primed (like.like_from): H2 and H3 take the middle of the run, H4 the
    first part's name."""
    from foldmap.layout import Placed

    def placed(k, x, y, label, start, length=1.0):
        sse = SimpleNamespace(start=start, chain="A", kind="H")
        return k, Placed(sse=sse, cx=x, cy=y, length=length, angle=0.0, width=0.5, label=label)

    a = dict(
        [placed("a1", 0, 5, "α1", 0), placed("a2", 5, 5, "α2", 10), placed("a3", 5, -4, "α3", 20)]
        + [placed("a6", 0, 0, "α6", 30, 4.0)]
    )
    b = dict(
        [placed("b1", 0, 5, "α1", 0), placed("b2", 5, 5, "α2", 10), placed("b3", 5, -4, "α3", 20)]
        + [placed("p", -1.25, 0, "α6", 30, 1.5), placed("m", 1.25, 0, "α6′", 40, 1.5)]
    )
    lay_a, lay_b = SimpleNamespace(placed=a, sheets=[]), SimpleNamespace(placed=b, sheets=[])
    mates = Mates(
        mate={"b1": "a1", "b2": "a2", "b3": "a3", "m": "a6"},
        pieces={"p": "a6"},
        extras=(),
        only_a=(),
        same_fold=True,
        smaller=4,
    )
    pairs = criteria.like_pairs(lay_a, lay_b, mates)
    assert pairs[-1] == ("a6", "m", ("p", "m"))
    assert criteria.place_deviation(lay_a, lay_b, pairs)["p90"] == 0.0
    assert criteria.relations(lay_a, lay_b, pairs)["share"] == 1.0
    assert criteria.h4_names(lay_a, lay_b, mates)["pass"] is True
