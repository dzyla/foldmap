# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""One module constant of foldmap perturbed in this process (plan 28 Task 9.2, ruling Q29): the hook
scripts/sensitivity.py runs scripts/random100.py under, in each of its processes (the run, every entry, every view).

    FOLDMAP_PERTURB="layout:_STACK_PAD=0.975" python -I scripts/random100.py run --not-for-paper ...
    FOLDMAP_PERTURB="layout:_STACK_PAD=0.975" python -I scripts/perturb.py reach   # what holds the value: JSON

install() reads FOLDMAP_PERTURB (module:NAME=VALUE, VALUE an int or float literal) and puts a finder first on
sys.meta_path that loads foldmap.<module> from its source with `NAME = VALUE` inserted after every module-level
statement that binds NAME (one inside a top-level if or try included), compiled from that tree and never written as
bytecode. The module's attribute is VALUE from that line on, so everything that reads it sees VALUE: the module's own
uses, what the module derives from it while it is imported (another constant, a function's default, a class
attribute: layout.LIKE_WEIGHT = _K_LIKE / _K_HOME), every module that binds it by `from .module import NAME` (they
import it after), and every reader through the module. That is more than setting the attribute after import, which
leaves the derived values and the defaults at the old value. install() must run before foldmap.<module> is imported,
and raises when it is imported already or no statement binds NAME (no silent no-op). Without FOLDMAP_PERTURB it does
nothing: random100.py and every other process run as before.

reach (this script's command) imports the perturbed module and every module that binds the name by `from ... import`
at module level, transitively through re-exports (binders, found by an AST scan of src/foldmap), and prints what each
holds; imports inside a function read the module's attribute when called and are listed as such.

FOLDMAP_TRACE_READS=FOLDER (install_trace, Python 3.12+: sys.monitoring) does not perturb: it records which of the
constants FOLDER/names.json lists are read while the process runs, so a constant no figure reads is told from one
that reads it and does not move. A read is the run of a line, inside a function body of src/foldmap, that names the
constant (its module's own name for it, a from-import binder's, or an attribute NAME); a module-level line, a default
or a constant derived at import is no read, and neither is a read of what was derived from it. Each process writes
FOLDER/<pid>.json at exit: its argv and the constants read. Each line location is reported once, then disabled, so
the drawing is not slowed much and not changed."""

from __future__ import annotations

import ast
import functools
import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import json
import os
import sys
from pathlib import Path

ENV = "FOLDMAP_PERTURB"
PACKAGE = "foldmap"
ROOT = Path(__file__).resolve().parents[1]
_TRY = tuple({ast.Try, getattr(ast, "TryStar", ast.Try)})


def parse(text: str) -> tuple[str, str, int | float]:
    """'module:NAME=VALUE' -> (module, NAME, VALUE); ValueError unless VALUE is an int or float literal."""
    target, eq, literal = text.partition("=")
    module, colon, name = target.strip().partition(":")
    if not (eq and colon and module and name.isidentifier() and literal.strip()):
        raise ValueError(f"{ENV}={text!r}: expected module:NAME=VALUE")
    value = ast.literal_eval(literal.strip())
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{ENV}={text!r}: VALUE must be an int or a float")
    return module, name, value


def _binds(stmt: ast.stmt, name: str) -> bool:
    if isinstance(stmt, ast.Assign):
        targets = stmt.targets
    elif isinstance(stmt, (ast.AnnAssign, ast.AugAssign)) and getattr(stmt, "value", None) is not None:
        targets = [stmt.target]
    else:
        return False
    return any(isinstance(n, ast.Name) and n.id == name for t in targets for n in ast.walk(t))


def rewrite(tree: ast.Module, name: str, value: int | float) -> int:
    """`name = value` inserted after every module-level statement of `tree` binding `name` (through top-level if and
    try blocks). -> how many were inserted."""
    count = 0

    def body(stmts: list[ast.stmt]) -> list[ast.stmt]:
        nonlocal count
        out = []
        for s in stmts:
            if isinstance(s, ast.If):
                s.body, s.orelse = body(s.body), body(s.orelse)
            elif isinstance(s, _TRY):
                s.body, s.orelse, s.finalbody = body(s.body), body(s.orelse), body(s.finalbody)
                for h in s.handlers:
                    h.body = body(h.body)
            out.append(s)
            if _binds(s, name):
                new = ast.Assign(targets=[ast.Name(id=name, ctx=ast.Store())], value=ast.Constant(value=value))
                out.append(ast.copy_location(new, s))
                count += 1
        return out

    tree.body = body(tree.body)
    ast.fix_missing_locations(tree)
    return count


class _Loader(importlib.machinery.SourceFileLoader):
    """A module's source loader that compiles the rewritten tree (rewrite) and never reads or writes bytecode."""

    def __init__(self, fullname: str, path: str, name: str, value: int | float):
        super().__init__(fullname, path)
        self.name, self.value, self.bound = name, value, 0

    def get_code(self, fullname):
        tree = ast.parse(self.get_data(self.path), self.path)
        self.bound = rewrite(tree, self.name, self.value)
        if not self.bound:
            raise ImportError(f"{ENV}: no module-level statement of {fullname} binds {self.name}")
        return compile(tree, self.path, "exec", dont_inherit=True)


class _Finder(importlib.abc.MetaPathFinder):
    def __init__(self, fullname: str, name: str, value: int | float):
        self.fullname, self.name, self.value = fullname, name, value
        self.loader: _Loader | None = None

    def find_spec(self, fullname, path, target=None):
        if fullname != self.fullname:
            return None
        found = importlib.machinery.PathFinder.find_spec(fullname, path)
        if found is None or not str(found.origin).endswith(".py"):
            raise ImportError(f"{ENV}: {fullname} has no source to perturb")
        self.loader = _Loader(fullname, found.origin, self.name, self.value)
        return importlib.util.spec_from_file_location(fullname, found.origin, loader=self.loader)


def install(text: str | None = None, package: str = PACKAGE) -> _Finder | None:
    """Perturb what `text` (default: FOLDMAP_PERTURB) names in this process, before its module is imported. None when
    there is nothing to perturb; RuntimeError when the module is imported already."""
    text = os.environ.get(ENV, "") if text is None else text
    if not text.strip():
        return None
    module, name, value = parse(text)
    fullname = f"{package}.{module}"
    if fullname in sys.modules:
        raise RuntimeError(f"{ENV}: {fullname} is imported already; the perturbation must be installed first")
    finder = _Finder(fullname, name, value)
    sys.meta_path.insert(0, finder)
    return finder


def uninstall(finder: _Finder | None) -> None:
    if finder is not None and finder in sys.meta_path:
        sys.meta_path.remove(finder)


# ------------------------------------------------------------------------------------------------- reach


def _module_level(tree: ast.Module) -> set[int]:
    """ids of the statements run when the module is imported (top level, through if and try blocks)."""
    out: set[int] = set()

    def body(stmts):
        for s in stmts:
            out.add(id(s))
            if isinstance(s, ast.If):
                body(s.body + s.orelse)
            elif isinstance(s, _TRY):
                body(s.body + s.orelse + s.finalbody)
                for h in s.handlers:
                    body(h.body)

    body(tree.body)
    return out


@functools.cache
def _trees(src: str) -> dict[str, ast.Module]:
    """Every module of the package at `src`, parsed once per process."""
    return {p.stem: ast.parse(p.read_text(encoding="utf-8")) for p in sorted(Path(src).glob("*.py"))}


@functools.cache
def _imports(src: str, package: str) -> tuple[tuple[str, str, str, str, str], ...]:
    """Every `from ... import` of the package at `src`, once per process: (importer, source module, name, bound as,
    "module" | "function")."""
    out = []
    for importer, tree in sorted(_trees(src).items()):
        top = _module_level(tree)
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom):
                continue
            if node.level == 1 and node.module:
                source = node.module
            elif node.level == 0 and (node.module or "").startswith(f"{package}."):
                source = node.module.split(".", 1)[1]
            else:
                continue
            at = "module" if id(node) in top else "function"
            out += [(importer, source, a.name, a.asname or a.name, at) for a in node.names]
    return tuple(out)


def binders(src: Path, module: str, name: str, package: str = PACKAGE) -> list[dict]:
    """Every `from ... import` of `name` from `module` in the package at `src`, transitively through module-level
    re-exports: [{"module", "as", "at": "module" | "function"}]."""
    imports = _imports(str(src), package)
    out, todo, seen = [], [(module, name)], set()
    while todo:
        source, bound = todo.pop()
        if (source, bound) in seen:
            continue
        seen.add((source, bound))
        for importer, src_module, imported, alias, at in imports:
            if src_module == source and imported == bound:
                out.append({"module": importer, "as": alias, "at": at})
                if at == "module":
                    todo.append((importer, alias))
    return out


def reach(text: str | None = None, src: Path | None = None, package: str = PACKAGE) -> dict:
    """Install the perturbation, import its module and every module-level binder, and say what each holds."""
    finder = install(text, package)
    if finder is None:
        raise SystemExit(f"{ENV} is not set")
    try:
        module = finder.fullname.split(".", 1)[1]
        if src is None:
            src = Path(importlib.util.find_spec(package).submodule_search_locations[0])
        got = {"perturbed": f"{module}:{finder.name}", "value": finder.value}
        mod = importlib.import_module(finder.fullname)
        got["module"] = getattr(mod, finder.name)
        got["statements"] = finder.loader.bound if finder.loader else 0
        held = {}
        for b in binders(src, module, finder.name, package):
            key = f"{b['module']}.{b['as']}"
            if b["at"] == "module":
                held[key] = getattr(importlib.import_module(f"{package}.{b['module']}"), b["as"])
            else:
                held[key] = "read when called"
    finally:
        uninstall(finder)  # what is imported keeps the value; nothing imported later is rewritten
    got["binders"] = held
    got["reached"] = got["module"] == finder.value and all(
        v == finder.value for v in held.values() if v != "read when called"
    )
    return got


# ------------------------------------------------------------------------------------------------- reads


TRACE_ENV = "FOLDMAP_TRACE_READS"


def reading_lines(src: Path, wanted: list[str], package: str = PACKAGE) -> dict[tuple[str, int], set[str]]:
    """(the module file's real path, line) -> the constants ("module:NAME" of `wanted`) a line inside a function body
    of the package at `src` names: as the defining module's name or a binder's (binders), or as an attribute NAME."""
    names: dict[str, dict[str, set[str]]] = {}  # module -> the name it reads a constant by -> constants
    attrs: dict[str, set[str]] = {}  # NAME -> constants read as an attribute x.NAME
    for cid in wanted:
        module, name = cid.split(":", 1)
        names.setdefault(module, {}).setdefault(name, set()).add(cid)
        attrs.setdefault(name, set()).add(cid)
        for b in binders(src, module, name, package):
            names.setdefault(b["module"], {}).setdefault(b["as"], set()).add(cid)
    out: dict[tuple[str, int], set[str]] = {}
    for stem, tree in _trees(str(src)).items():
        own = names.get(stem, {})
        real = os.path.realpath(src / f"{stem}.py")
        for fn in ast.walk(tree):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for stmt in fn.body:
                for n in ast.walk(stmt):
                    got = set()
                    if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load) and n.id in own:
                        got = own[n.id]
                    elif isinstance(n, ast.Attribute) and isinstance(n.ctx, ast.Load) and n.attr in attrs:
                        got = attrs[n.attr]
                    if got:
                        out.setdefault((real, n.lineno), set()).update(got)
    return out


def install_trace(folder: str | None = None, package: str = PACKAGE) -> dict | None:
    """Record, until this process exits, which constants of FOLDER/names.json are read (the module docstring);
    written to FOLDER/<pid>.json at exit. None without FOLDMAP_TRACE_READS. -> the live record {"reads": set}."""
    import atexit

    folder = os.environ.get(TRACE_ENV, "") if folder is None else folder
    if not folder.strip():
        return None
    mon = getattr(sys, "monitoring", None)
    if mon is None:
        raise RuntimeError(f"{TRACE_ENV} needs Python 3.12 or later (sys.monitoring)")
    wanted = json.loads((Path(folder) / "names.json").read_text(encoding="utf-8"))
    src = Path(importlib.util.find_spec(package).submodule_search_locations[0])
    targets = reading_lines(src, wanted, package)
    record = {"argv": sys.argv, "reads": set()}
    real: dict[str, str] = {}

    def line(code, lineno):
        f = real.get(code.co_filename)
        if f is None:
            f = real[code.co_filename] = os.path.realpath(code.co_filename)
        got = targets.get((f, lineno))
        if got:
            record["reads"].update(got)
        return mon.DISABLE  # each line location once

    def write():
        mon.set_events(tool, 0)
        doc = {"argv": record["argv"], "reads": sorted(record["reads"])}
        (Path(folder) / f"{os.getpid()}.json").write_text(json.dumps(doc) + "\n", encoding="utf-8")

    tool = next(t for t in range(6) if mon.get_tool(t) is None)
    mon.use_tool_id(tool, "foldmap-trace-reads")
    mon.register_callback(tool, mon.events.LINE, line)
    mon.set_events(tool, mon.events.LINE)
    atexit.register(write)
    return record


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv != ["reach"]:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    if (ROOT / "src" / PACKAGE).is_dir():  # the checkout's foldmap, as random100.py (-I: no PYTHONPATH)
        sys.path.insert(0, str(ROOT / "src"))
    got = reach()
    print(json.dumps(got, sort_keys=True))
    return 0 if got["reached"] else 1


if __name__ == "__main__":
    sys.exit(main())
