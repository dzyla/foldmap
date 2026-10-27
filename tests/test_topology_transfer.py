"""The record of the topology --like carries (scripts/topology_transfer.py, plan 28 Task 6, ruling Q44): the reasons
it gives where Like.order says nothing are like._order's, its summary counts what its rows hold, and the committed
record (skipped in the sdist, which has no docs/validation) is what its rows say."""

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from foldmap.align import Mates
from foldmap.like import _order

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "topology_transfer.py"
RECORD = ROOT / "docs" / "validation" / "homologues" / "topology-transfer" / "report.json"


@pytest.fixture(scope="module")
def tt():
    spec = importlib.util.spec_from_file_location("topology_transfer", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _case(xs, ups, directions, blocks=None):
    """A reference layout, a B sheet and mates: B's strands b0, b1, ... mated to a0, a1, ... drawn at x = xs, up (+1)
    or down (-1), in A's sheet blocks `blocks` (one block: all of them); test_like's _order_of."""
    n = len(xs)
    placed = {f"a{k}": SimpleNamespace(cx=x, angle=np.pi / 2 * u) for k, (x, u) in enumerate(zip(xs, ups))}
    lay_a = SimpleNamespace(placed=placed, sheets=blocks or [[f"a{k}" for k in range(n)]])
    sheet = SimpleNamespace(strands=[SimpleNamespace(id=f"b{k}") for k in range(n)], directions=directions)
    mates = Mates(mate={f"b{k}": f"a{k}" for k in range(n)}, pieces={}, extras=(), only_a=(), same_fold=True, smaller=n)
    return lay_a, sheet, mates


def test_each_reason_is_one_of_like_order(tt):
    """Where like._order says nothing, _why_none names its condition, in like._order's order."""
    cases = {
        "fewer than two mated": ([0.0], [1], [1], None),
        "mates in two of A's sheets": ([0.0, 1.1, 2.2], [1, -1, 1], [1, -1, 1], [["a0", "a1"], ["a2"]]),
        "mates point otherwise than B's pairing": ([0.0, 1.1, 2.2], [1, 1, 1], [1, -1, 1], None),
        "mates not in one order across A's sheet": ([0.0, 2.2, 1.1], [1, -1, 1], [1, -1, 1], None),
    }
    for why, args in cases.items():
        case = _case(*args)
        assert _order(*case) is None, why
        assert tt._why_none(*case) == why
    assert _order(*_case([0.0, 1.1, 2.2], [1, -1, 1], [1, -1, 1])) is True  # the reasons are only for None


def test_a_sheet_row_says_whether_b_reads_and_points_as_a(tt):
    """sheet_rows: a sheet Like.order sets, drawn in A's order and way up, reads and points as A's; drawn turned half
    round (order and every strand reversed), it does neither."""
    lay_a, sheet, mates = _case([0.0, 1.1, 2.2], [1, -1, 1], [1, -1, 1])
    sheet.closed, sheet.accessory = False, []
    like = SimpleNamespace(order=[True])
    same = SimpleNamespace(placed={f"b{k}": lay_a.placed[f"a{k}"] for k in range(3)})
    turned = SimpleNamespace(
        placed={f"b{k}": SimpleNamespace(cx=-p.cx, angle=-p.angle) for k, p in enumerate(lay_a.placed.values())}
    )
    (row,) = tt.sheet_rows(lay_a, same, [sheet], like, mates)
    assert row["order"] is True and row["reads"] and row["points"] and row["reads_all"] and "why" not in row
    (row,) = tt.sheet_rows(lay_a, turned, [sheet], like, mates)
    assert not row["reads"] and not row["points"]
    (row,) = tt.sheet_rows(lay_a, same, [sheet], SimpleNamespace(order=[None]), mates, withheld=True)
    assert row["order"] is None and row["why"] == "withheld" and row["reads"]  # like._order not asked


CHECKS = ("H1", "H2", "H3", "H4", "H5")


def _row(name, kind, carried, withheld, orders=((), ())):
    """A record row: each way's verdicts (H1-H5 in order) and its sheets in another order."""

    def way(verdicts, found):
        return {**dict(zip(CHECKS, verdicts)), "h2": None, "h3": None, "sheet_orders": list(found), "sheets": []}

    return {"pair": name, "kind": kind, "carried": way(carried, orders[0]), "withheld": way(withheld, orders[1])}


def test_the_summary_counts_its_rows_and_names_a_lost_check(tt):
    """summary: positives passing all five and each check, negatives refused, per way; the pairs drawn otherwise and
    those losing a check carried."""
    rows = [
        _row("fixed", "positive", [True] * 5, [True, False, False, True, True], ((), ("D, A for A, D",))),
        _row("lost", "positive", [True, True, False, True, True], [True] * 5),
        _row("same", "positive", [True] * 5, [True] * 5),
        _row("negative", "negative", [True, None, None, None, True], [True, None, None, None, True]),
    ]
    got = tt.summary(rows)
    assert got["carried"]["all_pass"] == 2 and got["withheld"]["all_pass"] == 2
    assert got["carried"]["pass"]["H3"] == 2 and got["withheld"]["pass"]["H3"] == 2
    assert got["carried"]["refused"] == 1 and got["carried"]["negatives"] == 1
    assert got["withheld"]["pairs_with_sheet_orders"] == ["fixed"] and got["carried"]["pairs_with_sheet_orders"] == []
    assert [c["pair"] for c in got["changed"]] == ["fixed", "lost"]
    assert got["lost_any"] == ["lost"]


@pytest.mark.skipif(not RECORD.is_file(), reason="docs/validation is not part of the sdist")
def test_the_record_is_what_its_rows_say(tt):
    """The committed record: its summary is summary() of its rows; 3WZ0-2AV5 passes H1-H5 with A's order carried
    and draws A's mirror image without it (the held-out run's failure); no pair loses an H check carried; every
    sheet Like.order sets reads and points as A's; the sets carry sets.json's tags."""
    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    for name, rows in doc["sets"].items():
        assert not [r for r in rows if "error" in r], name
        assert doc["summary"][name] == json.loads(json.dumps(tt.summary(rows))), name
        assert doc["summary"][name]["lost_any"] == []
        carried, withheld = doc["summary"][name]["carried"], doc["summary"][name]["withheld"]
        assert carried["reads_and_points"] == carried["order_set"] > 0 and withheld["order_set"] == 0
        assert list(withheld["none_why"]) == ["withheld"] and "withheld" not in carried["none_why"]
    (reg,) = [r for r in doc["sets"]["hom_heldout"] if r["pair"] == tt.REGRESSION]
    assert all(reg["carried"][k] is True for k in tt.CHECKS) and reg["carried"]["sheet_orders"] == []
    assert reg["withheld"]["H3"] is False and reg["withheld"]["sheet_orders"] == ["D, A, C, B for B, C, A, D"]
    assert doc["tags"] == tt.tags()
    assert len(doc["sets"]["hom_tuning"]) == 15 and len(doc["sets"]["hom_heldout"]) == 40
