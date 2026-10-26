"""The inventory of thresholds, special cases and fallbacks (plan 28 Task 9.1, ruling Q23): scripts/inventory.py's
scan against docs/validation/inventory/inventory.md, and the named fallbacks counted per figure (foldmap.fallback):
in Layout.fallbacks, the command line's notes and the random-100 report, never in the drawn figure."""

import importlib.util
import json
import re
from pathlib import Path

import pytest

from foldmap import fallback, route
from foldmap.cli import main
from foldmap.cli import make_layout as cli_make_layout
from foldmap.io import read_model
from foldmap.render import draw, save_svg
from foldmap.route import route_loops
from foldmap.spec import draw_spec, figure_spec

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "foldmap"
INVENTORY = ROOT / "docs" / "validation" / "inventory"
needs_inventory = pytest.mark.skipif(
    not (INVENTORY / "inventory.md").exists(), reason="the inventory is not part of the sdist"
)


def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"foldmap_{name}", ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def inv():
    return _script("inventory")


@pytest.fixture(scope="module")
def scanned(inv):
    return inv.scan(SRC)


# ------------------------------------------------------------------------------------------------ the scan

PLANTED = """
LIMIT = 3.5  # Å: a planted threshold
TRIVIAL_ONE = 1
A, B = 0.25, "x"


def f(x, tol=0.35, n=2):
    margin = 4.5
    if x > 7.25 and x != 0:
        return max(x * 0.8, 12)
    y = x[3] + x[:12]
    fallback.note("planted_curve")
    return round(y, 3)  # falls back to rounding
"""


def test_the_scan_finds_a_planted_constant_comparison_argument_and_fallback(inv):
    items = {it["id"]: it for it in inv.scan_source("planted", PLANTED)}
    kinds = {i: it["kind"] for i, it in items.items()}
    assert kinds["planted:LIMIT"] == "constant" and "planted threshold" in items["planted:LIMIT"]["comment"]
    assert kinds["planted:TRIVIAL_ONE"] == "constant"  # every module constant, trivial or not
    assert kinds["planted:A"] == "constant" and "planted:B" not in kinds  # a tuple target: each name its value
    assert kinds["planted:f:default:tol"] == "default" and "planted:f:default:n" not in kinds  # 2 is trivial
    assert kinds["planted:f:local:margin"] == "local"
    comparisons = [i for i, k in kinds.items() if k == "comparison"]
    assert len(comparisons) == 1 and comparisons[0].startswith("planted:f:comparison:7.25@")  # 0 is trivial
    arguments = sorted(it["value"] for it in items.values() if it["kind"] == "argument")
    assert arguments == ["0.8", "12", "3"]  # x[3] and x[:12] are indices, not listed
    assert items["planted:f:fallback:planted_curve"]["name"] == "planted_curve"
    assert kinds["planted:f:fallback-marker"] == "fallback-marker"


def test_an_inline_id_follows_its_code_not_its_line_number_comment_or_spacing(inv):
    """Ruling Q25 (a): the anchor is the line's code, its comment cut off, its whitespace dropped."""
    one = {it["id"] for it in inv.scan_source("m", PLANTED)}
    moved = {it["id"] for it in inv.scan_source("m", "\n\n# a comment\n" + PLANTED)}
    spaced = {it["id"] for it in inv.scan_source("m", PLANTED.replace("x > 7.25 and", "x>7.25  and"))}
    reworded = inv.scan_source("m", PLANTED.replace("# falls back to rounding", "# digits kept"))
    assert moved == one and spaced == one
    assert {it["id"] for it in reworded if it["kind"] != "fallback-marker"} == {
        i for i in one if not i.endswith("fallback-marker")
    }
    commented = inv.scan_source("m", PLANTED.replace("return max(x * 0.8, 12)", "return max(x * 0.8, 12)  # a note"))
    assert {it["id"] for it in commented} == one
    changed = {it["id"] for it in inv.scan_source("m", PLANTED.replace("x > 7.25", "x >= 7.25"))}
    assert changed != one and any(i.startswith("m:f:comparison:7.25@") for i in changed - one)


def test_the_merge_carries_a_row_to_its_edited_line_and_gives_new_rows_their_role(inv):
    """Ruling Q25 (b), (c): an edited line's row moves to its new id when the two alone share
    module:function:kind:value; a new row takes the origin its number's role fixes, else unknown."""
    before = inv.scan_source("m", PLANTED)
    rows = {
        it["id"]: {"id": it["id"], "origin": "named-entry:1ABC", "evidence": "by hand", "check": "-"} for it in before
    }
    after = inv.scan_source("m", PLANTED.replace("x > 7.25", "x >= 7.25") + "\n\ndef g(v):\n    return v < 1e-9\n")
    merged, diff = inv.merge(after, rows)
    old = next(i for i in rows if i.startswith("m:f:comparison:7.25@"))
    new = next(it["id"] for it in after if it["id"].startswith("m:f:comparison:7.25@"))
    assert diff["carried"] == [(old, new)] and old not in merged
    assert merged[new]["origin"] == "named-entry:1ABC" and merged[new]["evidence"] == "by hand"
    tol = next(it["id"] for it in after if it["function"] == "g")
    assert diff["added"] == [tol] and merged[tol]["origin"] == "convention"  # a float tolerance, by role
    assert not diff["dropped"]
    twice = PLANTED.replace("return max(x * 0.8, 12)", "return max(x * 0.8, 12) + max(x * 0.8, 12)")
    two = inv.scan_source("m", twice.replace("x > 7.25", "x >= 7.25"))
    _, diff = inv.merge(two, rows)  # two 0.8s and 12s on the edited line: no unique match, nothing carried there
    assert not any(o.startswith("m:f:argument:0.8@") for o, _ in diff["carried"])


def test_the_scan_reads_with_items_and_except_clauses(inv, scanned):
    """Ruling Q26 M1: a with statement's call (cath's urlopen timeout) is scanned as any other."""
    src = "def f(p):\n    with open(p, buffering=4096) as fh:\n        return fh\n"
    assert [it["value"] for it in inv.scan_source("m", src) if it["kind"] == "argument"] == ["4096"]
    assert any(it["module"] == "cath" and it["kind"] == "argument" and it.get("keyword") == "timeout" for it in scanned)
    assert any(it["module"] == "app" and it.get("keyword") == "width" and it["value"] == "150" for it in scanned)


def test_tests_naming_an_item_skip_short_and_common_names(inv):
    index = {"test_a": {"at", "compare", "draw", "layout", "_text_width"}, "test_b": {"_text_width"}}
    assert inv._naming("at", "compare", index) == [] and inv._naming("draw", "layout", index) == []
    assert inv._naming("_text_width", "layout", index) == ["test_a"]  # test_b names the function, not its module


def test_a_fallback_note_with_no_literal_name_is_reported(inv):
    """Ruling Q26 M6: the scan cannot read a computed name, so none may exist."""
    s = inv._scanned("m", "def f(name):\n    fallback.note(name)\n    fallback.note('loop_curve')\n")
    assert s.dynamic == ["m:2"] and [it["name"] for it in s.items if it["kind"] == "fallback"] == ["loop_curve"]
    assert inv.dynamic_notes(SRC) == []


@needs_inventory
def test_the_table_lists_every_item_the_scan_finds_and_no_other(inv, scanned):
    """Ruling Q23: a new threshold or fallback the table lacks fails here, and so does a row the code no longer has.
    Run python scripts/inventory.py, then complete the new rows (origin, evidence, check)."""
    rows = inv.read_table((INVENTORY / "inventory.md").read_text(encoding="utf-8"))
    got = inv.compare(scanned, rows)
    assert not got["missing"], f"items the table lacks (run scripts/inventory.py): {got['missing'][:20]}"
    assert not got["stale"], f"rows the scan no longer finds (run scripts/inventory.py): {got['stale'][:20]}"
    assert not got["bad_origin"], f"rows whose origin is none of {inv.ORIGINS}: {got['bad_origin'][:20]}"


@needs_inventory
def test_the_staleness_check_catches_both_ways(inv, scanned):
    rows = inv.read_table((INVENTORY / "inventory.md").read_text(encoding="utf-8"))
    planted = inv.scan_source("planted", PLANTED)
    got = inv.compare(scanned + planted, rows)
    assert set(got["missing"]) == {it["id"] for it in planted} and not got["stale"]
    extra = {**rows, "layout:NO_SUCH_THRESHOLD": {**next(iter(rows.values())), "id": "layout:NO_SUCH_THRESHOLD"}}
    got = inv.compare(scanned, extra)
    assert got["stale"] == ["layout:NO_SUCH_THRESHOLD"] and not got["missing"]
    bad = {**rows, scanned[0]["id"]: {**rows[scanned[0]["id"]], "origin": "a hunch"}}
    assert inv.compare(scanned, bad)["bad_origin"] == [scanned[0]["id"]]


@needs_inventory
def test_the_table_reads_back_what_it_writes(inv, scanned):
    rows = inv.read_table((INVENTORY / "inventory.md").read_text(encoding="utf-8"))
    again = inv.read_table(inv.render_table(scanned, rows, "# header\n"))
    assert {i: (r["origin"], r["evidence"], r["check"]) for i, r in again.items()} == {
        i: (r["origin"], r["evidence"], r["check"]) for i, r in rows.items()
    }


@needs_inventory
def test_the_committed_json_is_the_scan(scanned):
    doc = json.loads((INVENTORY / "inventory.json").read_text(encoding="utf-8"))
    assert [(it["id"], it["kind"], it["value"]) for it in doc["items"]] == [
        (it["id"], it["kind"], it["value"]) for it in scanned
    ], "inventory.json is stale: run python scripts/inventory.py"
    assert doc["source_stamp"]["script"]["name"] == "scripts/inventory.py"


# ------------------------------------------------------------------------------------------------ the fallbacks


def test_every_fallback_note_names_a_registered_fallback_and_every_one_is_noted(scanned):
    noted = {it["name"] for it in scanned if it["kind"] == "fallback"}
    assert noted == set(fallback.FALLBACKS)


def test_every_nucleic_acid_path_kind_is_counted():
    """route_na.route_nucleic counts the fallbacks of the paths as drawn by kind: a new kind would go uncounted."""
    kinds = set()
    for path in SRC.glob("*.py"):
        kinds |= set(re.findall(r'NAPath\(\s*"(\w+)"', path.read_text(encoding="utf-8")))
    kinds |= {"run", "joint"}  # na_runs._drawing's kind_drawn
    counted = re.findall(r'fell\.count\("(\w+)"\)', (SRC / "route_na.py").read_text(encoding="utf-8"))
    assert kinds == set(counted) == {"tail", "gap", "hairpin", "run", "joint", "trace"}


def test_counting_is_per_context_and_nested_counts_reach_the_outer():
    fallback.note("loop_curve")  # outside counting(): nothing to count it in, nothing happens
    with fallback.counting() as outer:
        fallback.note("loop_curve", 2)
        with fallback.counting() as inner:
            fallback.note("helix_name_moved")
            fallback.note("loop_curve", 0)
        assert fallback.counts(inner) == {"helix_name_moved": 1}
    assert fallback.counts(outer) == {"helix_name_moved": 1, "loop_curve": 2}
    assert fallback.line(outer) == "fallbacks taken: helix_name_moved 1, loop_curve 2" and fallback.line({}) == ""


@pytest.fixture
def no_routes(monkeypatch):
    """Every loop falls back: the router finds no route (a fallback made to fire)."""
    monkeypatch.setattr(route._Router, "route", lambda self, *a, **k: (None, False))


def _drawn_outside_counting(spec):
    """The figure draw_spec draws, made the same way outside fallback.counting()."""
    model = read_model(spec.structure, spec.assembly)
    layout, sses, bb = cli_make_layout(
        spec.structure, spec.mode, spec.rotate, spec.flip_v, model=model, like_out=None, **spec.layout_kwargs()
    )
    loops = route_loops(layout, sses, bb)
    return draw(layout, loops, sses, spec.title, model=model, look=spec.look), layout


def test_a_fallback_taken_is_counted_in_the_layout_and_never_drawn(ubq, no_routes):
    spec = figure_spec(ubq)
    drawn = draw_spec(spec)
    stuck = drawn.stuck
    assert stuck and drawn.fallbacks["loop_curve"] == len(stuck) == drawn.layout.fallbacks["loop_curve"]
    svg = save_svg(drawn.figure)
    assert "fallbacks taken" not in svg and "loop_curve" not in svg
    plain, layout = _drawn_outside_counting(spec)
    assert layout.fallbacks == {}  # made outside counting(): nothing counted
    assert save_svg(plain) == svg  # counting changes no figure


def test_the_command_line_notes_say_how_often_each_fallback_fired(ubq, tmp_path, capsys, no_routes):
    out = tmp_path / "u.svg"
    assert main(["plot", str(ubq), "-o", str(out)]) == 0
    err = capsys.readouterr().err
    said = re.search(r"foldmap: fallbacks taken: (.*)", err)
    assert said and re.search(r"\bloop_curve [1-9]\d*\b", said.group(1))
    assert "fallbacks taken" not in out.read_text(encoding="utf-8")


def test_the_notes_line_comes_only_with_a_fallback_taken(ubq, tmp_path, capsys):
    drawn = draw_spec(figure_spec(ubq))
    assert set(drawn.fallbacks) <= set(fallback.FALLBACKS)
    assert main(["plot", str(ubq), "-o", str(tmp_path / "u.svg")]) == 0
    err = capsys.readouterr().err
    assert ("fallbacks taken" in err) == bool(drawn.fallbacks)


def test_the_random100_report_counts_each_fallback_per_entry_and_in_all():
    r100 = _script("random100")
    entries = [
        {
            "id": "1ABC",
            "views": {
                "default": {"result": {"fallbacks": {"loop_curve": 2}}},
                "beta": {"result": {"fallbacks": {"loop_curve": 1, "helix_name_moved": 3}}},
            },
        },
        {"id": "2DEF", "views": {"default": {"result": {"fallbacks": {}}}}},
    ]
    got = r100.fallback_totals(entries)
    assert got["names"]["loop_curve"] == {"entries": 1, "default_view": 2, "every_view": 3}
    assert got["names"]["helix_name_moved"] == {"entries": 0, "default_view": 0, "every_view": 3}
    assert set(got["names"]) == set(fallback.FALLBACKS)
    assert got["per_entry"] == {"1ABC": {"loop_curve": 2}}
    md = "\n".join(r100.fallbacks_section(entries))
    assert "| loop_curve | 1 | 2 | 3 |" in md and "| 1ABC | loop_curve 2 |" in md
