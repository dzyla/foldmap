"""The readability floor (plan 28 Task 4, the review's 3.3, rulings Q39-Q41): floor.readability_floor names why a
figure is below it (text under 6 pt on its page, the route budget spent), the notes say it with the existing options
that draw less, `--strict-floor` writes nothing below it, `foldmap check` reports R1 (a fail exits 1); the time a
figure took is said, never judged; the default figure is unchanged. The floor's logic is tested on synthetic
measures; a figure below it is made from a committed fixture by holding its page small (layout.MAX_PAGE_IN), and a
slow one by its measures, never by timing anything."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from foldmap import deep, floor, layout
from foldmap.cli import main
from foldmap.science import READABILITY
from foldmap.spec import draw_spec, figure_spec

DATA = Path(__file__).parent / "data"
UBQ = DATA / "1UBQ.cif"
TETRAMER = DATA / "1BL8.cif"  # a C4 potassium channel: four chains, one per subunit
RECORD = Path(__file__).resolve().parents[1] / "docs" / "validation" / "readability-floor" / "report.json"


def _measures(**given) -> dict:
    return {"page_in": [7.0, 5.0], "text_pt": 7.37, "smallest_pt": 4.4, "fallbacks": {}, "seconds": 1.0} | given


@pytest.fixture
def small_page(monkeypatch):
    """Pages held at 3 in on their longer side: 1UBQ's text falls under 6 pt (page_scale scales the figure)."""
    monkeypatch.setattr(layout, "MAX_PAGE_IN", 3.0)


def test_each_reason_from_its_measure():
    assert floor.reasons_of(_measures()) == []
    assert floor.reasons_of(_measures(text_pt=6.0, smallest_pt=3.6)) == []  # the floor is the figure's text
    (got,) = floor.reasons_of(_measures(page_in=[48.0, 9.8], text_pt=2.049, smallest_pt=1.23))
    assert got["reason"] == floor.TEXT_SMALL and got["limit"] == layout.MIN_TEXT_PT == 6.0
    assert got["text"] == "its text is 2.0 pt, under 6 pt, on its 48 x 10 in page (its smallest text 1.2 pt)"
    (got,) = floor.reasons_of(_measures(text_pt=5.97))  # rounded down, as page_fit says it: never "6.0, under 6"
    assert "5.9 pt, under 6 pt" in got["text"]
    fired = {"route_second_search_refused": 2, "loop_curve": 3, "helix_name_moved": 9}
    (got,) = floor.reasons_of(_measures(fallbacks=fired))
    assert got["reason"] == floor.ROUTE_BUDGET and got["measured"] == 2
    assert got["text"] == "its route budget was spent: 2 wider route searches not run, 3 loops drawn as curves"
    assert floor.reasons_of(_measures(fallbacks={"loop_curve": 3})) == []  # a stuck loop is A3's, not the floor's
    assert set(floor.REASONS) == {floor.TEXT_SMALL, floor.ROUTE_BUDGET}  # time is said, never judged (Q40)


def test_below_judges_the_reasons_without_the_smallest_text():
    """floor.below (spec.draw_spec's judgment, plan 29 Task 5) gives a reason exactly when reasons_of does: the text
    size the page gives, the route budget spent; the smallest text and the time are never judged."""

    def drawn(pt, fired=None):
        fig = SimpleNamespace(get_size_inches=lambda: (48.0, 9.0), foldmap_text_pt=pt)
        return SimpleNamespace(figure=fig, fallbacks=fired or {})

    assert floor.below(drawn(5.026)) and not floor.below(drawn(6.0))
    assert floor.below(drawn(7.0, {"route_second_search_refused": 1}))
    assert not floor.below(drawn(7.0, {"loop_curve": 2}))


def test_time_is_a_note_never_a_reason():
    assert floor.reasons_of(_measures(seconds=107.9)) == []
    (got,) = floor.time_notes(_measures(seconds=107.9))
    assert got["reason"] == floor.TIME and got["limit"] == floor.TIME_LIMIT_S == 90
    assert got["text"] == "it took 108 s to make, over the hosted app's 90 s (said, not judged)"
    assert floor.time_notes(_measures(seconds=90.0)) == [] and floor.time_notes(_measures(seconds=None)) == []
    both = _measures(text_pt=5.48, seconds=107.9)
    assert [r["reason"] for r in floor.reasons_of(both)] == [floor.TEXT_SMALL]
    line = floor.note(floor.reasons_of(both), [], timed=floor.time_notes(both))
    assert line.startswith("below the readability floor: its text is 5.4 pt") and line.endswith("(said, not judged)")
    assert floor.note([], [], timed=floor.time_notes(both)).startswith("it took 108 s")


def test_the_options_in_each_interfaces_words():
    offered = [
        {"options": {"focus": "protomer"}, "what": "one subunit in full"},
        {"options": {"focus": "2M"}, "what": "chain 2M in full"},
        {"options": {"chains": "2M"}, "what": "chain 2M only"},
    ]
    reasons = floor.reasons_of(_measures(text_pt=2.0))
    cli = floor.note(reasons, offered)
    assert cli.startswith("below the readability floor: its text is 2.0 pt")
    assert cli.endswith(
        "drawing less may bring it above: --focus protomer (one subunit in full); --focus 2M (chain 2M in full); "
        "--chains 2M (chain 2M only)"
    )
    assert 'focus="protomer"' in floor.note(reasons, offered, "mcp") and 'chains=["2M"]' in floor.note(
        reasons, offered, "mcp"
    )
    app = floor.note(reasons, offered, "app")
    assert "protomer" not in app  # the app has no control for it
    assert "Chains & assembly, Large assemblies: these chains in full, 2M" in app
    assert floor.note([], offered) == "" and floor.note(reasons, []).endswith("pt)")


def test_each_surface_names_only_the_options_it_has():
    """Plan 29 Task 5 (ruling P29-7 b): the note never tells a user to do what their surface cannot. Every option
    suggestions() offers is named by the command line and the MCP tool; the app names the two it has controls for
    (Chains & assembly: Large assemblies' "these chains in full", Chains) and drops one subunit, which it has no
    control for, so a note offering only that offers nothing in the app."""
    one = {"focus": "protomer"}
    for options in (one, {"focus": "A"}, {"chains": "A"}):
        assert floor._phrase(options, "cli") and floor._phrase(options, "mcp")
    assert floor._phrase(one, "app") is None and floor.ways([{"options": one, "what": "x"}], "app") == ""
    assert floor._phrase({"focus": "A"}, "app") == "Chains & assembly, Large assemblies: these chains in full, A"
    assert floor._phrase({"chains": "A"}, "app") == "Chains & assembly, Chains: A"


def test_the_apps_words_are_its_controls():
    """The app's words (floor._SAID) are the controls app_controls.chains_section draws: Large assemblies offers
    auto, none and "these chains in full", no subunit."""
    pytest.importorskip("streamlit")  # the app extra: absent in the sdist and minimum-deps jobs
    import inspect

    from foldmap import app_controls

    assert "protomer" not in app_controls.FOCUS and "protomer" not in app_controls.FOCUS.values()
    assert floor._SAID["app"]["focus"].endswith(f"Large assemblies: {app_controls.FOCUS['chains']}, {{v}}")
    source = inspect.getsource(app_controls.chains_section)
    assert '"Chains & assembly"' in source and '"Large assemblies"' in source and 'heading(st, "Chains"' in source


def test_an_ordinary_figure_is_above_the_floor():
    drawn = draw_spec(figure_spec(UBQ))
    got = floor.measures(drawn)
    assert got["text_pt"] >= 6.0 and got["smallest_pt"] < got["text_pt"]
    assert drawn.seconds is not None and got["seconds"] == round(drawn.seconds, 2)
    assert drawn.floor == [] and floor.suggestions(drawn) == []  # one chain: nothing to draw less of
    assert not floor.below(drawn)


def test_a_page_held_small_is_below_the_floor(small_page):
    drawn = draw_spec(figure_spec(UBQ))
    (got,) = drawn.floor
    assert got["reason"] == floor.TEXT_SMALL and got["measured"] < 6.0
    assert floor.below(drawn)  # the judgment spec.draw_spec makes, without walking the figure's text
    assert drawn.fallbacks.get("page_text_small") == 1  # the same condition page_fit says and counts
    assert drawn.page_notes and "under 6 pt" in drawn.page_notes[0]


def test_suggestions_draw_less(tmp_path):
    drawn = draw_spec(figure_spec(TETRAMER))
    offered = floor.suggestions(drawn)
    first = drawn.bb.labels[0].chain
    assert [o["options"] for o in offered][-2:] == [{"focus": first}, {"chains": first}]
    alone = draw_spec(figure_spec(TETRAMER, chains=[first]))
    assert {p.sse.chain for p in alone.layout.placed.values()} == {first}


def _stub(sector=None, focus=None, chains=("A", "B")):
    """A drawn figure as suggestions() reads it: its frame's protomers, its focus, its protein chains."""
    lay = SimpleNamespace(frame=SimpleNamespace(sector=sector), focus_chains=focus, nucleic=None)
    return SimpleNamespace(layout=lay, bb=SimpleNamespace(labels=[SimpleNamespace(chain=c) for c in chains]))


def test_which_options_are_offered():
    def offered(drawn):
        return [o["options"] for o in floor.suggestions(drawn)]

    assert offered(_stub({"A": 0, "B": 1})) == [{"focus": "protomer"}, {"focus": "A"}, {"chains": "A"}]
    assert offered(_stub({"A": 0, "B": 0})) == [{"focus": "A"}, {"chains": "A"}]  # one protomer: no subunit to pick
    assert offered(_stub({"A": 0, "B": 1}, focus={"A"})) == [{"chains": "A"}]  # already one subunit, chain A
    assert offered(_stub(None, chains=("A",))) == []  # one chain: nothing to draw less of


def test_plot_says_it_and_strict_writes_nothing(small_page, tmp_path, capsys):
    out = tmp_path / "f.svg"
    assert main(["plot", str(UBQ), "-o", str(out)]) == 0
    err = capsys.readouterr().err
    assert "foldmap: warning: below the readability floor: its text is" in err and out.is_file()
    strict = tmp_path / "g.svg"
    assert main(["plot", str(UBQ), "-o", str(strict), "--strict-floor"]) == 1
    err = capsys.readouterr().err
    assert "no figure written: it is below the readability floor (--strict-floor)" in err
    assert not strict.exists()


def test_strict_writes_a_figure_above_the_floor(tmp_path, capsys):
    out = tmp_path / "f.svg"
    assert main(["plot", str(UBQ), "-o", str(out), "--strict-floor"]) == 0 and out.is_file()
    assert "readability floor" not in capsys.readouterr().err


def test_the_floor_is_never_in_the_file(small_page, tmp_path):
    """The floor is said in the notes, never drawn or written into the figure's file (its <desc> is the
    provenance's, which says the page, as before)."""
    out = tmp_path / "f.svg"
    assert main(["plot", str(UBQ), "-o", str(out)]) == 0
    svg = out.read_text(encoding="utf-8")
    assert "readability floor" not in svg and "route budget" not in svg


def test_check_reports_r1(small_page, capsys):
    code = main(["check", str(UBQ)])
    out = capsys.readouterr().out
    assert code == 1
    assert f"R1 {floor.TITLE}: fail" in out and "  its text is" in out
    code = main(["check", str(UBQ), "--json", "-"])
    doc = json.loads(capsys.readouterr().out)
    (r1,) = [c for c in doc["checks"] if c["code"] == floor.CODE]
    assert r1["group"] == deep.FLOOR and r1["kind"] == READABILITY and r1["status"] == deep.FAIL
    assert [r["reason"] for r in r1["measures"]["reasons"]] == [floor.TEXT_SMALL]
    assert doc["verdict"] == "fail" and doc["version"] == 1


def test_check_passes_r1_above_the_floor(capsys):
    assert main(["check", str(UBQ)]) == 0
    assert f"R1 {floor.TITLE}: pass" in capsys.readouterr().out


@pytest.fixture
def slow(monkeypatch):
    """Every figure measured as if it took 107.9 s (8XJV's time in the record), its text as drawn: no timing."""
    measured = floor.measures

    monkeypatch.setattr(floor, "measures", lambda drawn: {**measured(drawn), "seconds": 107.9})


def test_time_alone_passes_check(slow, capsys):
    """Q40: the time a figure took never fails R1 nor changes the exit code, and is not in the --json report."""
    assert main(["check", str(UBQ)]) == 0
    out = capsys.readouterr().out
    assert f"R1 {floor.TITLE}: pass" in out and "  it took 108 s to make, over the hosted app's 90 s" in out
    assert main(["check", str(UBQ), "--json", "-"]) == 0
    text = capsys.readouterr().out
    (r1,) = [c for c in json.loads(text)["checks"] if c["code"] == floor.CODE]
    assert r1["status"] == deep.PASS and "seconds" not in r1["measures"]["measures"]
    assert "107.9" not in text and "took" not in json.dumps(r1)  # no wall time in the report: the same every run


def test_strict_writes_a_figure_only_time_is_said_of(slow, tmp_path, capsys):
    out = tmp_path / "f.svg"
    assert main(["plot", str(UBQ), "-o", str(out), "--strict-floor"]) == 0 and out.is_file()
    err = capsys.readouterr().err
    assert "foldmap: warning: it took 108 s to make" in err and "below the readability floor" not in err


def test_a_floor_that_cannot_judge_never_stops_a_figure(small_page, monkeypatch, tmp_path, capsys):
    def broken(drawn):
        raise RuntimeError("no")

    monkeypatch.setattr(floor, "suggestions", broken)
    out = tmp_path / "f.svg"
    assert main(["plot", str(UBQ), "-o", str(out)]) == 0 and out.is_file()
    assert "foldmap: warning: readability floor not judged: RuntimeError" in capsys.readouterr().err


@pytest.mark.parametrize("where", ["measures", "suggestions"])
def test_strict_writes_nothing_the_floor_could_not_judge(small_page, monkeypatch, tmp_path, capsys, where):
    """Final review F6: --strict-floor fails closed: a floor that cannot judge the figure is no pass."""

    def broken(*args, **kwargs):
        raise RuntimeError("no")

    monkeypatch.setattr(floor, where, broken)
    out = tmp_path / "f.svg"
    assert main(["plot", str(UBQ), "-o", str(out), "--strict-floor"]) == 1
    err = capsys.readouterr().err
    assert "no figure written: the readability floor could not be judged (RuntimeError) (--strict-floor)" in err
    assert not out.exists()


def test_the_mcp_result_says_it(small_page, monkeypatch):
    pytest.importorskip("mcp")
    monkeypatch.delenv("FOLDMAP_OUTPUT_DIR", raising=False)
    import asyncio

    from foldmap.mcp_server import server

    r = asyncio.run(server.call_tool("draw_topology", {"source": str(UBQ)}))
    assert [x["reason"] for x in r.structured_content["floor"]] == [floor.TEXT_SMALL]
    assert "warning: below the readability floor" in r.content[0].text


def test_the_app_says_it(small_page):
    from foldmap import app, app_jobs

    blob = json.dumps({**app.LAYOUT_DEFAULTS, "theme": "publication"}, sort_keys=True)
    out = app_jobs.draw_job(UBQ, blob)
    assert out["floor_note"].startswith("below the readability floor: its text is")


def test_the_app_says_nothing_of_a_pair_page(small_page):
    """A compare or diff page is two panels on one page, not the figure the floor judges (B's layout under the
    page's text): the app says no floor of it."""
    from foldmap import app, app_jobs

    like = app.like_entry(str(UBQ), "1UBQ", b="1UBQ")
    like.update(page="compare")
    blob = json.dumps({"theme": "publication", "like": like}, sort_keys=True, default=str)
    out = app_jobs.draw_job(UBQ, blob)
    assert app.pair_kind(json.loads(blob)) == "compare" and out["floor_note"] == ""


@pytest.mark.skipif(not RECORD.is_file(), reason="docs/validation is not part of the sdist")
def test_the_record_is_what_the_floor_says_of_its_measures():
    """The committed record (scripts/readability_floor.py): each figure's reasons are reasons_of its measures, and
    every option the notes offered for a figure below the floor brought it above."""
    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    rows = doc["entries"]
    assert rows and {"9FZL", "8XJV"} <= {r["id"] for r in rows}
    for r in rows:
        for fig in [r["default"], *r["tried"]]:
            if "error" in fig:
                continue
            got = floor.reasons_of(fig["measures"])
            assert [x["reason"] for x in got] == fig["reasons"], (r["id"], fig.get("options"))
        if "error" not in r["default"]:
            timed = [x["reason"] for x in floor.time_notes(r["default"]["measures"])]
            assert timed == r["default"]["notes"], r["id"]  # said, never among the reasons (Q40)
        offered = [t for t in r["tried"] if t.get("offered")]
        if r["default"]["reasons"]:
            assert offered, r["id"]
        for t in offered:
            assert t["reasons"] == [], (r["id"], t["options"])
