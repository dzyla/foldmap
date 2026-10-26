# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The inventory of foldmap's thresholds, special cases and fallbacks (plan 28 Task 9.1, ruling Q23).

python scripts/inventory.py          # scan src/foldmap, write docs/validation/inventory/inventory.json and merge
                                     # inventory.md: rows kept by id with their hand columns, new rows added as
                                     # origin unknown, rows the scan no longer finds dropped (each printed)
python scripts/inventory.py --check  # exit 1 when inventory.md and the scan disagree (tests/test_inventory.py)

The scan is an AST scan of every module of src/foldmap. It lists, each item under a stable id:

  constant         a module-level assignment whose value holds a number (int or float, never bool), or a class's
                   UPPER_CASE one: the whole value, one item; id module:NAME (module:Class.NAME)
  default          a function parameter's numeric default holding a non-trivial number; id module:function:default:param
  local            a name bound in a function to a value made of numbers alone (x = 0.35, x = (3, 4.5)) holding a
                   non-trivial one; id module:function:local:name
  comparison       a non-trivial number inside a comparison (d < 3.5, 0.2 <= x <= y)
  argument         a non-trivial number inside a call's argument (f(x, 3.5), g(tol=0.2), h(x * 0.8))
  fallback         a fallback.note(NAME) call: a named fallback path, counted per figure (foldmap.fallback);
                   id module:function:fallback:NAME
  fallback-marker  a function whose text (code, comments, docstring) speaks of falling back (fallback, falls back,
                   fell back, last resort, gives up, degraded; not the counting of fallbacks, COUNTING): the table
                   says which named fallback it holds, or why it holds none; id module:function:fallback-marker

Comparison and argument ids are module:function:kind:value@anchor, the anchor the first 8 hex digits of the SHA-1 of
the item's source line (whitespace collapsed), so an id survives lines moving and changes when its line does; the
k-th repeat of an id in one function takes "#k". Trivial numbers (TRIVIAL): 0, 1, 2, -1 and 0.5, as int or float.
Not listed: numbers in arithmetic outside those places (a bare x * 0.35, a return value), indices and slices (a[3],
a[:12]), numbers in f-strings, and dataclass field defaults (Style's keys, documented in style.py); a module-level
statement other than an assignment is scanned as function "<module>".

Each item carries the comments at it (the trailing comments of its statement's lines, a compound statement's header
only; the comment block directly above it; for a constant, the block carried on under it when a blank line ends it)
and the test files naming it (a constant by its name, a function's item by the function's name, a fallback
by its quoted name): grep evidence, not proof that a test pins the value. The hand columns of inventory.md (origin,
evidence, check) are kept across merges by id."""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import json
import re
import sys
import tokenize
from collections import Counter
from pathlib import Path

from foldmap.source_stamp import source_stamp

SOURCE_STAMP = source_stamp(__file__)  # the code this inventory was taken of (a code inventory, not a measurement)
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "foldmap"
TESTS = ROOT / "tests"
OUT = ROOT / "docs" / "validation" / "inventory"
TRIVIAL = (0, 1, 2, -1, 0.5)
KINDS = ("constant", "default", "local", "comparison", "argument", "fallback", "fallback-marker")
ORIGINS = ("biology", "geometry", "measured", "named-entry", "convention", "unknown")
MARKERS_SKIPPED = ("fallback",)  # the registry itself: every line of it speaks of fallbacks
# the counting of fallbacks (a note, the registry, the notes' line, Layout.fallbacks), not a fallback: no marker
COUNTING = re.compile(
    r"fallback\.(?:note|counting|counts|line|FALLBACKS)\b|import fallback\b|\bfallback_line\b|\.fallbacks\b"
    r"|\bFALLBACKS\b|named fallbacks?\b|from \.fallback import"
)
MARKER = re.compile(r"(?i)fall(?:s|ing)?[ -]?back(?!s)|fell back|last resort|gives? up|giving up|degraded")
GENERATED = "<!-- inventory:generated (everything below is rewritten by scripts/inventory.py) -->"
COLUMNS = ("id", "kind", "value", "at", "origin", "evidence", "tests naming it", "check")
HAND = ("origin", "evidence", "check")


def _number(node: ast.AST) -> int | float | None:
    """The number a literal is (an int or float constant, or one under a unary minus or plus), else None."""
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        got = _number(node.operand)
        return None if got is None else (-got if isinstance(node.op, ast.USub) else got)
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return node.value
    return None


def _trivial(value) -> bool:
    return any(value == t for t in TRIVIAL)


def _literals(node: ast.AST) -> list[ast.AST]:
    """The number literals under `node` (a negative number once, as its unary minus), outside f-strings."""
    out = []

    def walk(n):
        if isinstance(n, ast.JoinedStr):
            return
        if _number(n) is not None:
            out.append(n)
            return
        for c in ast.iter_child_nodes(n):
            walk(c)

    walk(node)
    return out


def _numbers_only(node: ast.AST) -> bool:
    """A value made of numbers alone: literals joined by operators, tuples and lists."""
    if _number(node) is not None:
        return True
    if isinstance(node, ast.BinOp):
        return _numbers_only(node.left) and _numbers_only(node.right)
    if isinstance(node, (ast.Tuple, ast.List)):
        return bool(node.elts) and all(_numbers_only(e) for e in node.elts)
    return False


def _squash(text: str, most: int = 0) -> str:
    text = " ".join(text.split())
    return text if not most or len(text) <= most else text[: most - 1] + "…"


def _anchor(line: str) -> str:
    return hashlib.sha1(_squash(line).encode("utf-8")).hexdigest()[:8]


def _comments(source: str) -> tuple[dict[int, str], set[int]]:
    """Each line's comment, and the lines holding a comment alone."""
    said, alone = {}, set()
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, SyntaxError):
        return said, alone
    code_lines = set()
    for t in tokens:
        if t.type == tokenize.COMMENT:
            said[t.start[0]] = t.string.lstrip("#").strip()
        elif t.type not in (tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT, tokenize.ENDMARKER):
            code_lines.update(range(t.start[0], t.end[0] + 1))
    alone = {n for n in said if n not in code_lines}
    return said, alone


class _Scan(ast.NodeVisitor):
    """One module's items (see the module docstring)."""

    def __init__(self, module: str, source: str):
        self.module, self.source = module, source
        self.lines = source.splitlines()
        self.said, self.alone = _comments(source)
        self.items: list[dict] = []
        self.stack: list[str] = []
        self.seen: Counter = Counter()
        self.statement: list[ast.stmt] = []

    # ----------------------------------------------------------------------------------------------- helpers
    def qual(self) -> str:
        return ".".join(self.stack) or "<module>"

    def comment_at(self, first: int, last: int, below: bool = False) -> str:
        """The comment block directly above line `first`, then the comments on lines first..last; with `below`, then
        the comment block directly under line `last` when a blank line or the end ends it (a constant's comment
        carried on under it, not the next one's)."""
        above, n = [], first - 1
        while n >= 1 and n in self.alone:
            above.append(self.said[n])
            n -= 1
        on = [self.said[k] for k in range(first, last + 1) if k in self.said and k not in self.alone]
        under, n = [], last + 1
        while below and n in self.alone:
            under.append(self.said[n])
            n += 1
        if under and n <= len(self.lines) and self.lines[n - 1].strip():
            under = []
        return _squash(" ".join([*reversed(above), *on, *under]), 600)

    def add(
        self, kind: str, ident: str, node: ast.AST, value: str, stmt: ast.AST | None = None, below=False, **extra
    ) -> None:
        self.seen[ident] += 1
        if self.seen[ident] > 1:
            ident = f"{ident}#{self.seen[ident]}"
        stmt = stmt or node
        line = node.lineno
        last = getattr(stmt, "end_lineno", None) or stmt.lineno
        body = getattr(stmt, "body", None)
        if isinstance(body, list) and body:  # a compound statement (if, for, def ...): its header lines only
            last = max(stmt.lineno, body[0].lineno - 1)
        self.items.append(
            {
                "id": ident,
                "kind": kind,
                "module": self.module,
                "function": self.qual(),
                "line": line,
                "value": value,
                "at": _squash(self.lines[line - 1], 160),
                "comment": self.comment_at(stmt.lineno, last, below),
                **extra,
            }
        )

    def segment(self, node: ast.AST, most: int = 120) -> str:
        return _squash(ast.get_source_segment(self.source, node) or ast.unparse(node), most)

    # ----------------------------------------------------------------------------------------------- module level
    def scan_module(self, tree: ast.Module) -> None:
        self.parents = {c: p for p in ast.walk(tree) for c in ast.iter_child_nodes(p)}
        self.constant_nodes: set[int] = set()
        self.body(tree.body, cls=None)
        self.visit_rest(tree)
        if self.module not in MARKERS_SKIPPED:
            self.markers(tree)

    def body(self, stmts, cls: str | None) -> None:
        """Module-level (or class-level) constants, through top-level if / try blocks."""
        for s in stmts:
            if isinstance(s, (ast.If, ast.Try)):
                self.body(s.body + s.orelse + getattr(s, "finalbody", []), cls)
                for h in getattr(s, "handlers", []):
                    self.body(h.body, cls)
                continue
            if isinstance(s, ast.ClassDef) and cls is None:
                self.body(s.body, cls=s.name)
                continue
            if isinstance(s, ast.Assign):
                targets, value = s.targets, s.value
            elif isinstance(s, ast.AnnAssign) and s.value is not None:
                targets, value = [s.target], s.value
            else:
                continue
            pairs = []  # (name, its value): A, B = 1, 2 gives each name its own
            for t in targets:
                if isinstance(t, ast.Name):
                    pairs.append((t.id, value))
                elif isinstance(t, ast.Tuple) and all(isinstance(e, ast.Name) for e in t.elts):
                    same = isinstance(value, ast.Tuple) and len(value.elts) == len(t.elts)
                    pairs += [(e.id, value.elts[k] if same else value) for k, e in enumerate(t.elts)]
            pairs = [(n, v) for n, v in pairs if not n.startswith("__")]
            if cls is not None:  # a class body: its UPPER_CASE constants only (dataclass fields are not listed)
                pairs = [(n, v) for n, v in pairs if n.upper() == n and n.strip("_").isupper()]
            for name, v in pairs:
                if not [n for n in _literals(v) if not self._indexed(n, v)]:
                    continue
                for n in ast.walk(value):
                    self.constant_nodes.add(id(n))
                ident = f"{self.module}:{cls + '.' if cls else ''}{name}"
                self.add("constant", ident, s, self.segment(v), stmt=s, below=True)

    def _indexed(self, lit: ast.AST, top: ast.AST) -> bool:
        """`lit` is an index or slice bound (a[3], a[:12]), up to `top`."""
        node = lit
        while node is not top and node in self.parents:
            p = self.parents[node]
            if isinstance(p, ast.Slice) or (isinstance(p, ast.Subscript) and node is p.slice):
                return True
            node = p
        return False

    # ----------------------------------------------------------------------------------------------- functions
    def visit_rest(self, tree: ast.Module) -> None:
        self.visit(tree)

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self.stack.append(node.name)
        self.generic_visit(node)
        self.stack.pop()

    def visit_FunctionDef(self, node) -> None:
        self.stack.append(node.name)
        self.inline(node)  # its decorators' numbers
        args = node.args
        pos = args.posonlyargs + args.args
        for a, d in zip(pos[len(pos) - len(args.defaults) :], args.defaults):
            self.default(a, d, node)
        for a, d in zip(args.kwonlyargs, args.kw_defaults):
            if d is not None:
                self.default(a, d, node)
        for s in node.body:
            self.visit(s)
        self.stack.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def default(self, arg: ast.arg, value: ast.AST, fn) -> None:
        if _numbers_only(value) and any(not _trivial(_number(n)) for n in _literals(value)):
            self.add("default", f"{self.module}:{self.qual()}:default:{arg.arg}", value, self.segment(value), stmt=fn)

    def generic_visit(self, node: ast.AST) -> None:
        if id(node) in self.constant_nodes:
            return
        if isinstance(node, ast.stmt):
            self.statement.append(node)
            self.inline(node)
        super().generic_visit(node)
        if isinstance(node, ast.stmt):
            self.statement.pop()

    def visit_Assign(self, node: ast.Assign) -> None:
        if id(node.value) in self.constant_nodes:
            return
        if self.stack and _numbers_only(node.value) and any(not _trivial(_number(n)) for n in _literals(node.value)):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    self.add("local", f"{self.module}:{self.qual()}:local:{t.id}", node, self.segment(node.value))
            return
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        f = node.func
        if (
            isinstance(f, ast.Attribute)
            and f.attr == "note"
            and isinstance(f.value, ast.Name)
            and f.value.id == "fallback"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            name = node.args[0].value
            stmt = self.statement[-1] if self.statement else node
            self.add("fallback", f"{self.module}:{self.qual()}:fallback:{name}", node, name, stmt=stmt, name=name)
        self.generic_visit(node)

    def inline(self, stmt: ast.stmt) -> None:
        """The comparisons' and arguments' numbers of one statement's own expressions (not its nested statements)."""
        for _, value in ast.iter_fields(stmt):
            nodes = value if isinstance(value, list) else [value]
            for n in nodes:
                if isinstance(n, ast.expr):
                    for lit in _literals(n):
                        self.place(lit, stmt)

    def place(self, lit: ast.AST, stmt: ast.stmt) -> None:
        value = _number(lit)
        if _trivial(value):
            return
        node, kw = lit, None
        while node is not stmt:
            p = self.parents[node]
            if isinstance(p, ast.Slice) or (isinstance(p, ast.Subscript) and node is p.slice):
                return
            if isinstance(p, (ast.JoinedStr, ast.FormattedValue)):
                return
            if isinstance(p, ast.keyword):
                kw = p.arg
            if isinstance(p, ast.Compare):
                kind, callee = "comparison", None
                break
            if isinstance(p, ast.Call) and node is not p.func:
                kind, callee = "argument", self.segment(p.func, 60)
                break
            if isinstance(p, ast.Lambda) and node is not p.body:
                return
            node = p
        else:
            return
        text = self.segment(lit, 30)
        ident = f"{self.module}:{self.qual()}:{kind}:{text}@{_anchor(self.lines[lit.lineno - 1])}"
        extra = {"callee": callee, "keyword": kw} if kind == "argument" else {}
        self.add(kind, ident, lit, text, stmt=stmt, **extra)

    # ----------------------------------------------------------------------------------------------- markers
    def markers(self, tree: ast.Module) -> None:
        spans = []  # (first, last, qualname) of every function, innermost found last when sorted by first line

        def walk(n, stack):
            for c in ast.iter_child_nodes(n):
                if isinstance(c, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    q = [*stack, c.name]
                    if not isinstance(c, ast.ClassDef):
                        first = min([c.lineno, *(d.lineno for d in c.decorator_list)])
                        spans.append((first, c.end_lineno, ".".join(q)))
                    walk(c, q)
                else:
                    walk(c, stack)

        walk(tree, [])
        hits: dict[str, list[int]] = {}
        for k, line in enumerate(self.lines, 1):
            if MARKER.search(COUNTING.sub("", line)):
                inner = [s for s in spans if s[0] <= k <= s[1]]
                q = max(inner, key=lambda s: s[0])[2] if inner else "<module>"
                hits.setdefault(q, []).append(k)
        for q, lines in hits.items():
            first = lines[0]
            self.items.append(
                {
                    "id": f"{self.module}:{q}:fallback-marker",
                    "kind": "fallback-marker",
                    "module": self.module,
                    "function": q,
                    "line": first,
                    "value": str(len(lines)),
                    "at": _squash(self.lines[first - 1], 160),
                    "comment": _squash(" / ".join(self.lines[k - 1] for k in lines), 400),
                    "lines": lines,
                }
            )


def scan_source(module: str, source: str) -> list[dict]:
    """The items of one module's source."""
    s = _Scan(module, source)
    s.scan_module(ast.parse(source))
    return s.items


def scan(src: Path = SRC, tests: Path | None = TESTS) -> list[dict]:
    """Every item of every module under `src`, in module then line order, each with the test files naming it."""
    items = []
    for path in sorted(src.glob("*.py")):
        items += scan_source(path.stem, path.read_text(encoding="utf-8"))
    if tests is not None:
        index = {
            p.stem: set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", p.read_text(encoding="utf-8")))
            for p in sorted(tests.glob("test_*.py"))
        }
        quoted = {
            p.stem: set(re.findall(r"[\"']([a-z_]+)[\"']", p.read_text(encoding="utf-8")))
            for p in sorted(tests.glob("test_*.py"))
        }
        for it in items:
            if it["kind"] == "fallback":
                it["tests"] = [t for t, words in quoted.items() if it["name"] in words]
            else:
                word = (
                    it["id"].split(":")[1].split(".")[-1] if it["kind"] == "constant" else it["function"].split(".")[-1]
                )
                it["tests"] = [t for t, words in index.items() if word in words] if word != "<module>" else []
    return items


# ------------------------------------------------------------------------------------------------- the table


def _cells(row: str) -> list[str]:
    parts = re.split(r"(?<!\\)\|", row.strip())
    return [p.strip().replace("\\|", "|") for p in parts[1:-1]]


def _cell(text) -> str:
    return _squash(str(text)).replace("|", "\\|")


def read_table(text: str) -> dict[str, dict]:
    """inventory.md's rows by id: {column: text}, the id unquoted."""
    rows = {}
    for line in text.splitlines():
        if not line.startswith("| `"):
            continue
        cells = _cells(line)
        if len(cells) != len(COLUMNS):
            continue
        row = dict(zip(COLUMNS, cells))
        row["id"] = row["id"].strip("`")
        rows[row["id"]] = row
    return rows


def origin_kind(origin: str) -> str:
    """The origin's class (named-entry:7VVK,1EYG is named-entry); "" when it is none of ORIGINS."""
    head = origin.split(":", 1)[0].strip()
    return head if head in ORIGINS and (head != "named-entry" or ":" in origin) else ""


def compare(items: list[dict], rows: dict[str, dict]) -> dict[str, list[str]]:
    """What the scan and the table disagree on: ids the table lacks (missing), ids the scan no longer finds
    (stale), and rows whose origin is none of ORIGINS (bad_origin)."""
    ids = [it["id"] for it in items]
    have = set(ids)
    return {
        "missing": [i for i in ids if i not in rows],
        "stale": sorted(i for i in rows if i not in have),
        "bad_origin": sorted(i for i, r in rows.items() if i in have and not origin_kind(r["origin"])),
    }


def counts(items: list[dict], rows: dict[str, dict]) -> dict:
    """Items by kind, and by kind and origin class."""
    by = Counter()
    for it in items:
        o = origin_kind((rows.get(it["id"]) or {}).get("origin", "unknown")) or "unknown"
        by[(it["kind"], o)] += 1
    return {
        "by_kind": dict(Counter(it["kind"] for it in items)),
        "by_kind_origin": {k: {o: by[(k, o)] for o in ORIGINS if by[(k, o)]} for k in KINDS},
        "by_origin": {o: sum(by[(k, o)] for k in KINDS) for o in ORIGINS},
    }


def _tests(names: list[str]) -> str:
    names = [n.removeprefix("test_") for n in names]
    return ", ".join(names[:4]) + (f" +{len(names) - 4}" if len(names) > 4 else "") if names else "-"


def render_table(items: list[dict], rows: dict[str, dict], header: str) -> str:
    """inventory.md: `header` (the hand-written part above GENERATED) then the counts, the fallback registry and one
    table per module, each row's hand columns taken from `rows` (origin unknown, evidence and check '-' for a new
    one)."""
    from foldmap.fallback import FALLBACKS

    c = counts(items, rows)
    out = [header.rstrip(), "", GENERATED, "", "## Counts", ""]
    out += ["| kind | items | " + " | ".join(ORIGINS) + " |", "|---|---|" + "---|" * len(ORIGINS)]
    for k in KINDS:
        got = c["by_kind_origin"][k]
        out.append(f"| {k} | {c['by_kind'].get(k, 0)} | " + " | ".join(str(got.get(o, 0)) for o in ORIGINS) + " |")
    out.append(f"| all | {len(items)} | " + " | ".join(str(c["by_origin"][o]) for o in ORIGINS) + " |")
    sites = Counter(it["name"] for it in items if it["kind"] == "fallback")
    out += ["", "## Named fallbacks (foldmap.fallback.FALLBACKS)", "", "| name | sites | what it is |", "|---|---|---|"]
    out += [f"| `{n}` | {sites[n]} | {_cell(what)} |" for n, what in FALLBACKS.items()]
    module = None
    for it in items:
        if it["module"] != module:
            module = it["module"]
            out += ["", f"## {module}", "", "| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
        r = rows.get(it["id"]) or {}
        hand = [r.get("origin") or "unknown", r.get("evidence") or "-", r.get("check") or "-"]
        out.append(
            f"| `{it['id']}` | {it['kind']} | {_cell(it['value'])} | {_cell(it['at'])} | "
            + " | ".join(_cell(h) for h in hand[:2])
            + f" | {_tests(it.get('tests', []))} | {_cell(hand[2])} |"
        )
    return "\n".join(out) + "\n"


def write(items: list[dict], out: Path = OUT) -> dict[str, list[str]]:
    """inventory.json, then inventory.md merged by id. -> what the merge added and dropped."""
    out.mkdir(parents=True, exist_ok=True)
    md = out / "inventory.md"
    text = md.read_text(encoding="utf-8") if md.is_file() else "# foldmap inventory\n"
    rows = read_table(text)
    header = text.split(GENERATED)[0]
    diff = compare(items, rows)
    doc = {
        "what": "foldmap's thresholds, special cases and fallbacks, scanned from src/foldmap (scripts/inventory.py)",
        "source_stamp": SOURCE_STAMP,
        "trivial": list(TRIVIAL),
        "counts": counts(items, rows),
        "items": items,
    }
    (out / "inventory.json").write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    md.write_text(render_table(items, rows, header), encoding="utf-8")
    return diff


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--check", action="store_true", help="only compare inventory.md with the scan; exit 1 if they differ"
    )
    args = ap.parse_args(argv)
    items = scan()
    if args.check:
        md = OUT / "inventory.md"
        got = compare(items, read_table(md.read_text(encoding="utf-8") if md.is_file() else ""))
        for k, ids in got.items():
            for i in ids:
                print(f"{k}: {i}")
        return 1 if any(got.values()) else 0
    diff = write(items)
    for i in diff["missing"]:
        print(f"added (origin unknown): {i}")
    for i in diff["stale"]:
        print(f"dropped (no longer found): {i}")
    print(f"{len(items)} items: " + ", ".join(f"{k} {n}" for k, n in Counter(it["kind"] for it in items).items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
