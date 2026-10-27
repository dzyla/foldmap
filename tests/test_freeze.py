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
HOW = "see CONTRIBUTING.md: a named mechanism and a measurement on an unseen set, recorded under docs/validation/"


def _inventory():
    spec = importlib.util.spec_from_file_location("foldmap_inventory_freeze", ROOT / "scripts" / "inventory.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def special_ids(rows: dict[str, dict], origin_kind) -> set[str]:
    """The ids whose origin is a special case or a tuned number."""
    return {i for i, r in rows.items() if origin_kind(r["origin"]) in SPECIAL}


def names_record(evidence: str, root: Path = ROOT) -> bool:
    """True when the evidence names a docs/validation path that exists (a directory or a file)."""
    return any((root / m.group(0)).exists() for m in RECORD.finditer(evidence))


def new_special(rows: dict[str, dict], baseline: set[str], origin_kind, root: Path = ROOT) -> list[str]:
    """Special-case rows outside the baseline whose evidence names no existing record."""
    return sorted(
        i for i in special_ids(rows, origin_kind) - baseline if not names_record(rows[i].get("evidence", ""), root)
    )


@pytest.fixture(scope="module")
def inv():
    return _inventory()


@pytest.fixture(scope="module")
def rows(inv):
    return inv.read_table((INVENTORY / "inventory.md").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def baseline():
    doc = json.loads((INVENTORY / "freeze-baseline.json").read_text(encoding="utf-8"))
    return set(doc["ids"])


def test_no_new_special_case_without_a_record(inv, rows, baseline):
    bad = new_special(rows, baseline, inv.origin_kind)
    assert not bad, (
        f"{len(bad)} inventory rows are named-entry, measured or unknown and not in "
        f"docs/validation/inventory/freeze-baseline.json: {bad[:20]}. A new special case or tuned number needs "
        f"its evidence to name a docs/validation/... record; {HOW}."
    )


def test_the_baseline_is_ids_of_the_table(inv, rows, baseline):
    """The baseline is not padded: every id in it is a row of the table (a row removed from the code leaves the
    baseline too, so it cannot be reused for a new special case)."""
    assert not baseline - set(rows), f"baseline ids the table no longer has: {sorted(baseline - set(rows))[:20]}"


def test_the_freeze_catches_a_new_row_and_accepts_one_with_a_record(inv, rows, baseline, tmp_path):
    row = {"origin": "named-entry:1ABC", "evidence": "-", "check": "-"}
    new = {**rows, "layout:_PLANTED": row}
    assert new_special(new, baseline, inv.origin_kind) == ["layout:_PLANTED"]
    for origin in ("measured", "unknown"):
        assert new_special({**rows, "layout:_P": {**row, "origin": origin}}, baseline, inv.origin_kind) == ["layout:_P"]
    # a row that is not a special case is not the freeze's concern
    for origin in ("biology", "geometry", "convention"):
        assert not new_special({**rows, "layout:_P": {**row, "origin": origin}}, baseline, inv.origin_kind)
    # evidence naming a record that exists passes; one that does not exist fails
    made = tmp_path / "docs" / "validation" / "mine"
    made.mkdir(parents=True)
    (made / "report.md").write_text("x", encoding="utf-8")
    ok = {**row, "evidence": "measured on random100-fourth, docs/validation/mine/report.md"}
    assert not new_special({**rows, "layout:_P": ok}, baseline, inv.origin_kind, tmp_path)
    gone = {**row, "evidence": "docs/validation/not-there/report.md"}
    assert new_special({**rows, "layout:_P": gone}, baseline, inv.origin_kind, tmp_path) == ["layout:_P"]
    # an entry of the baseline is accepted without one
    assert not new_special(new, baseline | {"layout:_PLANTED"}, inv.origin_kind)


def test_the_failure_points_to_contributing(inv, rows, baseline):
    assert (ROOT / "CONTRIBUTING.md").is_file()
    with pytest.raises(AssertionError, match=r"CONTRIBUTING\.md"):
        new = {**rows, "layout:_PLANTED": {"origin": "unknown", "evidence": "-", "check": "-"}}
        bad = new_special(new, baseline, inv.origin_kind)
        assert not bad, f"{bad}: {HOW}"
