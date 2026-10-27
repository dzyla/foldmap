"""The freeze (plan 28 Task 9.5, ruling Q53): no new special case and no new tuned number without a named mechanism
and a measurement on an unseen set (CONTRIBUTING.md).

The inventory (docs/validation/inventory/inventory.md) gives every constant, threshold and fallback an origin. A row
whose origin is named-entry (set for a named entry), measured (a tuned number) or unknown (not yet classified) is a
special case or a tuned number. freeze-baseline.json lists those that existed when the freeze began; a row outside it
fails here unless its evidence names a committed record under docs/validation/, the measurement CONTRIBUTING.md asks
for. docs/validation is not in the sdist, so this module is skipped where it is absent."""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "docs" / "validation" / "inventory"

if not (INVENTORY / "inventory.md").is_file():  # the sdist's signature
    pytest.skip("docs/validation is not part of the sdist", allow_module_level=True)

SPECIAL = ("named-entry", "measured", "unknown")  # the origins that mean a special case or a tuned number
RECORD = re.compile(r"docs/validation/[A-Za-z0-9_./-]*[A-Za-z0-9_]")
ANCHOR = re.compile(r"@[0-9a-f]{8}(?:#\d+)?$")  # a comparison's or argument's line anchor and repeat number
HOW = (
    "see CONTRIBUTING.md: a named mechanism and a measurement on an unseen set, recorded under docs/validation/, "
    "named in the row's evidence"
)


def _inventory():
    spec = importlib.util.spec_from_file_location("foldmap_inventory_freeze", ROOT / "scripts" / "inventory.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def key(ident: str) -> str:
    """An id less its line anchor (scripts/inventory.py's _key): reformatting the line keeps the key, a changed
    value, function or kind does not."""
    return ANCHOR.sub("", ident)


def keys(ids) -> dict[str, int]:
    """The keys of ids with how many ids share each."""
    out: dict[str, int] = {}
    for i in ids:
        out[key(i)] = out.get(key(i), 0) + 1
    return out


def row_names(ident: str) -> set[str]:
    """What a record must say to name the row `ident`: its `module:NAME`, or its name (the constant's or the
    function's, each dotted part of it too)."""
    module, name = key(ident).split(":")[:2]
    return {f"{module}:{name}", name, *name.split(".")} - {""}


def names_record(evidence: str, root: Path = ROOT, ident: str | None = None) -> bool:
    """True when the evidence names a docs/validation file that exists (a folder is no record: plan 28 final review
    F13) and, given the row's id, the file names the row (row_names)."""
    for m in RECORD.finditer(evidence):
        path = root / m.group(0)
        if not path.is_file():
            continue
        if ident is None or any(n in path.read_text(encoding="utf-8", errors="replace") for n in row_names(ident)):
            return True
    return False


def _beyond(ids: list[str], known: dict[str, int]) -> list[str]:
    """The ids past the count `known` allows for their key (the later ones by id)."""
    by: dict[str, list[str]] = {}
    for i in sorted(ids):
        by.setdefault(key(i), []).append(i)
    return sorted(i for k, v in by.items() for i in v[known.get(k, 0) :])


def new_special(rows: dict[str, dict], baseline: dict, origin_kind, root: Path = ROOT) -> list[str]:
    """Special-case rows (named-entry, measured, unknown) past the baseline's count for their key (the count first,
    then the record: a twin of a baseline row needs its own record), whose evidence names no existing record."""
    special = [i for i, r in rows.items() if origin_kind(r["origin"]) in SPECIAL]
    beyond = _beyond(special, baseline["special"])
    return [i for i in beyond if not names_record(rows[i].get("evidence", ""), root, i)]


def unevidenced(rows: dict[str, dict], baseline: dict, origin_kind) -> list[str]:
    """New rows (past the baseline's count of all rows for their key) given biology, geometry or convention by hand
    with no evidence. The scan's own role_origin gives its rows an evidence text, so '-' means set by hand."""
    new = _beyond(list(rows), baseline["all"])
    return [
        i
        for i in new
        if origin_kind(rows[i]["origin"]) not in SPECIAL and rows[i].get("evidence", "-").strip() in ("-", "")
    ]


def _freeze_script():
    spec = importlib.util.spec_from_file_location("foldmap_freeze_baseline", ROOT / "scripts" / "freeze_baseline.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def failure(bad: list[str], unevid: list[str] = ()) -> str:
    """The message of the failing test."""
    out = []
    if bad:
        out.append(
            f"{len(bad)} inventory rows are named-entry, measured or unknown and not in "
            f"docs/validation/inventory/freeze-baseline.json: {bad[:20]}. A new special case or tuned number needs "
            f"its evidence to name a docs/validation/... record"
        )
    if unevid:
        out.append(f"{len(unevid)} new rows have a hand-set origin and no evidence: {list(unevid)[:20]}")
    return "; ".join(out) + f"; {HOW}." if out else ""


@pytest.fixture(scope="module")
def inv():
    return _inventory()


@pytest.fixture(scope="module")
def rows(inv):
    return inv.read_table((INVENTORY / "inventory.md").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def baseline():
    return json.loads((INVENTORY / "freeze-baseline.json").read_text(encoding="utf-8"))


def test_no_new_special_case_without_a_record(inv, rows, baseline):
    bad, unevid = new_special(rows, baseline, inv.origin_kind), unevidenced(rows, baseline, inv.origin_kind)
    assert not (bad or unevid), failure(bad, unevid)


def test_the_baseline_is_not_padded(rows, baseline):
    """Every key in the baseline is a row of the table, with at most the rows' count: a special case removed from
    the code leaves the baseline too (regenerate it), so its place cannot be reused for a new one."""
    have = keys(rows)
    for part in ("special", "all"):
        over = {k: n for k, n in baseline[part].items() if have.get(k, 0) < n}
        assert not over, (
            f"baseline {part} keys the table no longer has (shrink it: python scripts/freeze_baseline.py): "
            f"{list(over)[:20]}"
        )


ROW = {"origin": "named-entry:1ABC", "evidence": "-", "check": "-"}


def _plant(rows, baseline, ident, **row):
    base = {"special": dict(baseline["special"]), "all": dict(baseline["all"])}
    return {**rows, ident: {**ROW, **row}}, base


def test_the_freeze_catches_a_new_row_and_accepts_one_with_a_record(inv, rows, baseline, tmp_path):
    new, base = _plant(rows, baseline, "layout:_PLANTED")
    assert new_special(new, base, inv.origin_kind) == ["layout:_PLANTED"]
    for origin in ("measured", "unknown"):
        new, base = _plant(rows, baseline, "layout:_P", origin=origin)
        assert new_special(new, base, inv.origin_kind) == ["layout:_P"]
    # a row that is not a special case is not this gate's concern
    for origin in ("biology", "geometry", "convention"):
        new, base = _plant(rows, baseline, "layout:_P", origin=origin, evidence="a quoted comment")
        assert not new_special(new, base, inv.origin_kind)
    # evidence naming a record that exists and names the row passes; one that does not exist fails
    made = tmp_path / "docs" / "validation" / "mine"
    made.mkdir(parents=True)
    (made / "report.md").write_text("`layout:_P` measured: x", encoding="utf-8")
    new, base = _plant(rows, baseline, "layout:_P", evidence="docs/validation/mine/report.md")
    assert not new_special(new, base, inv.origin_kind, tmp_path)
    new, base = _plant(rows, baseline, "layout:_P", evidence="docs/validation/not-there/report.md")
    assert new_special(new, base, inv.origin_kind, tmp_path) == ["layout:_P"]
    # F13: a folder is no record, and a record that never names the row is not its measurement
    new, base = _plant(rows, baseline, "layout:_P", evidence="docs/validation/mine")
    assert new_special(new, base, inv.origin_kind, tmp_path) == ["layout:_P"]
    new, base = _plant(rows, baseline, "layout:_Q", evidence="docs/validation/mine/report.md")
    assert new_special(new, base, inv.origin_kind, tmp_path) == ["layout:_Q"]
    # a function's row is named by its function
    (made / "fn.md").write_text("na_ends._reach: the 2 * grid.res reach, measured", encoding="utf-8")
    new, base = _plant(rows, baseline, "na_ends:_reach:comparison:2@deadbeef", evidence="docs/validation/mine/fn.md")
    assert not new_special(new, base, inv.origin_kind, tmp_path)


def test_a_reformatted_line_keeps_its_key_and_a_new_one_does_not(inv, rows, baseline):
    """Q54: the baseline is keyed on the id less its line anchor, with multiplicity."""
    ident = next(i for i, r in rows.items() if ANCHOR.search(i) and inv.origin_kind(r["origin"]) in SPECIAL)
    # the same literal on a reformatted line: the id's anchor changes, the key does not
    moved = {i: r for i, r in rows.items() if i != ident}
    moved[ANCHOR.sub("@deadbeef", ident)] = rows[ident]
    assert not new_special(moved, baseline, inv.origin_kind)
    # a second literal of the same key in the same function is new
    twin = {**rows, ANCHOR.sub("@deadbeef", ident): rows[ident]}
    got = new_special(twin, baseline, inv.origin_kind)
    assert len(got) == 1 and key(got[0]) == key(ident)
    # another value is a different key: new
    other = {i: r for i, r in rows.items() if i != ident}
    other[ANCHOR.sub("", ident).rsplit(":", 1)[0] + ":99.5@deadbeef"] = rows[ident]
    assert len(new_special(other, baseline, inv.origin_kind)) == 1


def test_a_hand_set_origin_needs_evidence(inv, rows, baseline):
    new, base = _plant(rows, baseline, "layout:_P", origin="convention")
    assert unevidenced(new, base, inv.origin_kind) == ["layout:_P"]
    new, base = _plant(rows, baseline, "layout:_P", origin="convention", evidence="what the unit is")
    assert not unevidenced(new, base, inv.origin_kind)
    assert not unevidenced(rows, baseline, inv.origin_kind)  # the table as committed


def test_the_failure_points_to_contributing(inv, rows, baseline):
    assert (ROOT / "CONTRIBUTING.md").is_file()
    assert "CONTRIBUTING.md" in failure(["layout:_P"]) and "CONTRIBUTING.md" in failure([], ["layout:_P"])


def test_a_baseline_with_a_record_does_not_excuse_a_twin_without_one(inv, rows, baseline, tmp_path):
    """A new row whose key matches a baseline row, with no record, still fails even when the baseline row itself
    names a record (the count is applied before the record filter)."""
    ident = next(i for i, r in rows.items() if inv.origin_kind(r["origin"]) in SPECIAL and ANCHOR.search(i))
    made = tmp_path / "docs" / "validation" / "mine"
    made.mkdir(parents=True)
    (made / "report.md").write_text(f"`{key(ident)}`", encoding="utf-8")  # it names the row
    twin = ANCHOR.sub("@ffffffff", ident)  # sorts after any real anchor: the later id is the one past the count
    assert ANCHOR.search(ident).group(0) < ANCHOR.search(twin).group(0)
    planted = {**rows, ident: {**rows[ident], "evidence": "docs/validation/mine/report.md"}, twin: dict(rows[ident])}
    assert new_special(planted, baseline, inv.origin_kind, tmp_path) == [twin]


def test_regenerating_the_baseline_only_shrinks_it(inv, rows, baseline):
    """Q55: scripts/freeze_baseline.py lowers counts and drops entries; a new special case is never admitted."""
    fb = _freeze_script()
    assert fb.shrunk_baseline(baseline, rows, inv.origin_kind)["special"] == baseline["special"]  # as committed
    ident = next(i for i, r in rows.items() if inv.origin_kind(r["origin"]) in SPECIAL and ANCHOR.search(i))
    # a special case removed: the baseline shrinks by one
    gone = {i: r for i, r in rows.items() if i != ident}
    new = fb.shrunk_baseline(baseline, gone, inv.origin_kind)
    assert sum(new["special"].values()) == sum(baseline["special"].values()) - 1
    assert all(n <= baseline["special"].get(k, 0) for k, n in new["special"].items())
    # a new special case and a new key: not admitted
    planted = {**rows, ANCHOR.sub("@deadbeef", ident): dict(rows[ident]), "layout:_PLANTED": dict(ROW)}
    new = fb.shrunk_baseline(baseline, planted, inv.origin_kind)
    assert new["special"] == baseline["special"] and new["all"] == baseline["all"]
