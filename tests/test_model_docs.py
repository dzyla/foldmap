"""The model of the system (plan 28 Task 9.4, ruling Q51): docs/model holds one page per subsystem, each saying what
it assumes, guarantees, how it fails and which checks see the failure. These tests keep the pages from drifting from
the code and the records: every code name a page quotes exists, every check id is defined, every record path and
test it names exists, and every page has the six sections.

What a page quotes, in backticks (docs/model/index.md says the same):

- `module.symbol` (`layout._open_at_ends`, `like.Like.order`): a name of live code in foldmap, resolved by importing
  the module and walking the attributes (a dataclass field counts);
- `module:NAME` (`layout:_NAME_OPEN`): a record id, as the inventory, sensitivity and prune records write them, for a
  name the code may no longer have (a pruned constant); it must appear in one of those records;
- a check id (`S1`, `A17`, `A1c`, `R1`, `D2`, `H3`): defined by science.SCIENCE, criteria.CRITERIA,
  criteria.LIKE_CHECKS, floor.CODE or deep.DEEP;
- a path under docs/, tests/, scripts/ or src/ (`docs/validation/prune/report.md`): it exists; `tests/x.py::name`
  names a test function defined in that file.

Relative links between pages and to the records must resolve too. docs/model and docs/validation are not in the
sdist (MANIFEST.in grafts only docs/examples, docs/names and docs/domains), so this module is skipped where they are
absent."""

from __future__ import annotations

import ast
import dataclasses
import importlib
import pkgutil
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "docs" / "model"
VALIDATION = ROOT / "docs" / "validation"

if not (MODEL / "index.md").exists() or not VALIDATION.is_dir():
    pytest.skip("docs/model and docs/validation are not part of the sdist", allow_module_level=True)

import foldmap  # noqa: E402
from foldmap import criteria, deep, floor, science  # noqa: E402

SECTIONS = ("Scope", "Assumes", "Guarantees", "Fails", "Seen by", "Evidence")
PAGES = sorted(p for p in MODEL.glob("*.md") if p.name != "index.md")
MODULES = {m.name for m in pkgutil.iter_modules(foldmap.__path__)}
RECORDS = (
    VALIDATION / "inventory" / "inventory.md",
    VALIDATION / "sensitivity" / "report.md",
    VALIDATION / "prune" / "report.md",
)
FILE_SUFFIXES = (".md", ".json", ".py", ".cif", ".gz", ".svg", ".png", ".txt", ".toml", ".yml", ".yaml", ".sh", ".in")

TICKED = re.compile(r"`([^`\n]+)`")
CODE = re.compile(r"[a-z_][a-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)+")
RECORD_ID = re.compile(r"([a-z_][a-z0-9_]*):([A-Za-z_][A-Za-z0-9_]*)")
CHECK = re.compile(r"(?:S\d+|A\d+[a-z]?|R\d+|D\d+|H\d+)")
PATH = re.compile(r"(?:docs|tests|scripts|src)/[^\s:`]+")
TEST = re.compile(r"(tests/[\w/]+\.py)::(\w+)(?:\[[^\]]*\])?")
LINK = re.compile(r"\]\(([^)\s]+)\)")


def _checks() -> set[str]:
    return {*science.SCIENCE, *criteria.CRITERIA, *criteria.LIKE_CHECKS, floor.CODE, *deep.DEEP}


def _resolves(name: str) -> bool:
    """Does `module.attr[.attr...]` name live code in foldmap?"""
    module, *parts = name.split(".")
    obj = importlib.import_module(f"foldmap.{module}")
    for part in parts:
        if hasattr(obj, part):
            obj = getattr(obj, part)
            continue
        fields = {f.name for f in dataclasses.fields(obj)} if dataclasses.is_dataclass(obj) else set()
        return part in fields or part in getattr(obj, "__annotations__", {})
    return True


def _defines(path: Path, name: str) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return any(isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name for n in ast.walk(tree))


def problems(page: Path) -> list[str]:
    """What on a page does not resolve: code names, record ids, check ids, paths, tests and relative links."""
    text = page.read_text(encoding="utf-8")
    records = "\n".join(p.read_text(encoding="utf-8") for p in RECORDS)
    checks = _checks()
    out = []
    for token in TICKED.findall(text):
        token = token.strip()
        if m := TEST.fullmatch(token):
            path = ROOT / m.group(1)
            if not path.exists() or not _defines(path, m.group(2)):
                out.append(f"test {token}: not defined")
        elif PATH.fullmatch(token):
            if not (ROOT / token.rstrip("/")).exists():
                out.append(f"path {token}: missing")
        elif CHECK.fullmatch(token):
            if token not in checks:
                out.append(f"check {token}: not defined in the code")
        elif m := RECORD_ID.fullmatch(token):
            if m.group(1) not in MODULES:
                out.append(f"record id {token}: no module {m.group(1)}")
            elif f"`{token}`" not in records:
                out.append(f"record id {token}: in no inventory, sensitivity or prune record")
        elif CODE.fullmatch(token) and not token.endswith(FILE_SUFFIXES):
            if token.split(".")[0] not in MODULES:
                out.append(f"code name {token}: no module {token.split('.')[0]}")
            elif not _resolves(token):
                out.append(f"code name {token}: not in the code")
    for target in LINK.findall(text):
        if target.startswith(("http://", "https://", "#")):
            continue
        if not (page.parent / target.split("#")[0]).exists():
            out.append(f"link {target}: missing")
    return out


def test_every_page_is_listed_in_the_index_and_every_listed_page_exists():
    index = (MODEL / "index.md").read_text(encoding="utf-8")
    listed = {Path(t).name for t in LINK.findall(index) if not t.startswith(("http", "#", "../"))}
    assert listed == {p.name for p in PAGES}
    assert not problems(MODEL / "index.md")


def test_the_pages_are_the_ten_subsystems():
    assert len(PAGES) == 10


def test_no_file_name_has_a_colon():
    assert not [p.name for p in MODEL.rglob("*") if ":" in p.name]


@pytest.mark.parametrize("page", PAGES, ids=lambda p: p.stem)
def test_each_page_has_the_six_sections_in_order(page):
    heads = [h.strip() for h in re.findall(r"^## (.+)$", page.read_text(encoding="utf-8"), flags=re.M)]
    assert heads == list(SECTIONS)


@pytest.mark.parametrize("page", PAGES, ids=lambda p: p.stem)
def test_every_name_check_record_and_test_a_page_quotes_exists(page):
    assert problems(page) == []


@pytest.mark.parametrize("page", PAGES, ids=lambda p: p.stem)
def test_every_page_quotes_code_and_a_record(page):
    """A page with no code name or no record path would be prose the anti-drift checks cannot hold."""
    ticked = [t.strip() for t in TICKED.findall(page.read_text(encoding="utf-8"))]
    assert any(CODE.fullmatch(t) and t.split(".")[0] in MODULES for t in ticked)
    assert any(t.startswith("docs/validation/") for t in ticked)


def test_a_drifted_name_is_caught(tmp_path):
    page = tmp_path / "x.md"
    page.write_text(
        "`layout._no_such_name` `layout:_NO_SUCH` `A99` `docs/validation/no-such` `tests/test_layout.py::test_none` "
        "[x](no-such.md) `layout._open_at_ends` `like.Like.order` `S1` `layout:_NAME_OPEN`",
        encoding="utf-8",
    )
    found = problems(page)
    assert len(found) == 6
    assert all("_open_at_ends" not in p and "Like.order" not in p and "_NAME_OPEN" not in p for p in found)
