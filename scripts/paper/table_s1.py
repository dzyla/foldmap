# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Table S1 (the parameters) from the code's constants at the release, never typed in (plan 28 Task 1).

    python scripts/paper/table_s1.py [--src DIR] [--check-against supplement.tex]

Reads each value from the foldmap sources at the harness commit (build/tree/src/foldmap, exported by extract.py;
--src names another package directory, such as src/foldmap of this checkout) with ast, never by importing them:

  C(module, name)              a module-level constant (an expression of literals, other constants of the module and
                               np.radians / math.radians), `key` for one entry of a dict
  C(module, name, operands)    the literal operands of a product such as dssp._Q1Q2F = 0.084 * 332.0
  D(module, function, arg)     a function's default argument (layout._lbfgs's maxiter)
  Inline(module, function, re) a number written inline in a function: the regex's one group, which must match
                               exactly once in that function's source

A value that cannot be found fails the run, naming it, so a renamed or reworded constant cannot drop out of the
table unseen. Writes build/tables/params.tex (the longtable, caption and label included: \\input it in place of the
typed table) and build/tables/params.json (each value with its module, name and line). --check-against compares the
generated values with a typed table (the rows of tab:params in that .tex, in order) and lists every difference.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import TABLES, TREE  # noqa: E402


@dataclass(frozen=True)
class C:
    module: str
    name: str
    key: str | None = None
    operands: bool = False


@dataclass(frozen=True)
class D:
    module: str
    function: str
    arg: str


@dataclass(frozen=True)
class Inline:
    module: str
    function: str
    pattern: str


# (module, quantity, value template, refs). The template's #name or #name:spec is replaced by that ref's value:
# spec a format spec (",", ".2f", "g") or deg (radians in degrees), pct (a fraction as a percentage), range (a pair as
# a--b), list (a tuple, comma-separated), of (a pair as "a of b"), times (a product's operands), pow10 (10^n),
# inv10 (1/x as 10^n), turns (signed turns as "0 to ±max by step"). No spec: str(value).
# fmt: off
ROWS: list[tuple[str, str, str, dict]] = [
    ("io", r"peptide bond C--N; carbonyl C=O rebuilt; farthest O counted", r"#a\angst; #b\angst; #c\angst",
     {"a": C("io", "_PEPTIDE_BOND_MAX"), "b": C("io", "_CO_BOND"), "c": C("io", "_CO_MAX")}),
    ("io", r"\texttt{auto} assembly, at most", "#a chains", {"a": C("io", "AUTO_MAX")}),
    ("io", r"duplicate copy: atoms within distance", r"#a:pct, #b\angst",
     {"a": C("io", "COINCIDENT_FRAC"), "b": C("io", "COINCIDENT_DIST")}),
    ("io", r"one protein: identity, length", "#a, #b", {"a": C("io", "SAME_PROTEIN"), "b": C("io", "SAME_LENGTH")}),
    ("io", r"protein--nucleotide contact (heavy atoms)", r"#a\angst", {"a": C("io", "CONTACT")}),
    ("io", r"disulfide S$\gamma$--S$\gamma$; glycan link; metal reach; ligand reach; covalent",
     r"#a; #b; #c; #d; #e\angst",
     {"a": C("io", "_SS_MAX"), "b": C("io", "_GLYCAN_LINK_MAX"), "c": C("io", "_METAL_REACH"),
      "d": C("io", "_LIGAND_REACH"), "e": C("io", "_BOND")}),
    ("dssp", r"energy constant; bond cut; floor; C$\alpha$ cut-off; closest atoms",
     r"#a:times; $#b$; $#c$ kcal\,mol$^{-1}$; #d\angst; #e\angst",
     {"a": C("dssp", "_Q1Q2F", operands=True), "b": C("dssp", "_HBOND_E"), "c": C("dssp", "_HBOND_E_MIN"),
      "d": C("dssp", "_CA_CUTOFF"), "e": C("dssp", "_MIN_DIST")}),
    ("ss", r"shortest strand; helix (H residues or all)", "#a; #b or #c residues",
     {"a": C("ss", "_MIN_STRAND"), "b": C("ss", "_MIN_HELIX_H"), "c": C("ss", "_MIN_HELIX_ANY")}),
    ("ss", r"kink split; shortest arm", r"#a:g$^\circ$; #b residues",
     {"a": C("ss", "KINK_SPLIT"), "b": C("ss", "_MIN_KINK_ARM")}),
    ("ss", r"bulge join: partner gap; copy break", "#a; #b residues",
     {"a": C("ss", "_PARTNER_GAP"), "b": C("ss", "_BREAK_MAX")}),
    ("sheets", r"bridges for neighbours; barrel minimum; split gap", "#a; #b strands; #c residues",
     {"a": C("sheets", "_MIN_PAIRS"), "b": C("sheets", "_BARREL_MIN"), "c": C("sheets", "_SPLIT_GAP")}),
    ("features", r"element contact: C$\alpha$ distance, pairs", r"#a\angst, #b",
     {"a": C("features", "CONTACT_CUTOFF"), "b": C("features", "MIN_CONTACT_PAIRS")}),
    ("features", r"bundle: length, $|\cos|$, distance, cover", r"#a residues, #b, #c\angst, #d",
     {"a": C("features", "BUNDLE_MIN"), "b": C("features", "BUNDLE_COS"), "c": C("features", "BUNDLE_DIST"),
      "d": C("features", "BUNDLE_COVER")}),
    ("features", r"helix on sheet (opt-in): length, strands, reach, touch, $|\cos|$, height, edge",
     r"#a, #b, #c\angst, #d, #e, #f:range\angst, #g\angst",
     {"a": C("features", "ON_SHEET_MIN"), "b": C("features", "ON_SHEET_STRANDS"), "c": C("features", "ON_SHEET_REACH"),
      "d": C("features", "ON_SHEET_TOUCH"), "e": C("features", "ON_SHEET_FLAT"), "f": C("features", "ON_SHEET_HEIGHT"),
      "g": C("features", "ON_SHEET_EDGE")}),
    ("features", r"domains: contact, sequence separation", r"#a\angst, #b",
     {"a": C("features", "DOMAIN_CONTACT"), "b": C("features", "DOMAIN_SEPARATION")}),
    ("features", r"domains: normalized cut, domain size, piece size, blocks", "#a, #b, #c residues, #d",
     {"a": C("features", "DOMAIN_NCUT"), "b": C("features", "DOMAIN_MIN_RESIDUES"),
      "c": C("features", "DOMAIN_MIN_SEGMENT"), "d": C("features", "DOMAIN_MAX_POINTS")}),
    ("features", r"domains: stray run, element residues per domain", "#a, #b",
     {"a": C("features", "DOMAIN_STRAY"), "b": C("features", "DOMAIN_MIN_SSE_RESIDUES")}),
    ("symmetry", r"RMSD tolerance; identity; shortest chain", r"#a\angst; #b; #c residues",
     {"a": C("symmetry", "TOLERANCE"), "b": C("symmetry", "SAME"), "c": C("symmetry", "MIN_RES")}),
    ("symmetry", r"least rotation; highest order; candidates per subunit", r"#a:deg$^\circ$; #b; #c",
     {"a": C("symmetry", "MIN_TURN"), "b": C("symmetry", "MAX_ORDER"), "c": C("symmetry", "_NEAREST")}),
    ("symmetry", r"helical: rise, twist, copies (asked); pseudo-$2_1$",
     r"#a\angst, #b:deg$^\circ$, #c (#d); #e:g$^\circ$",
     {"a": C("symmetry", "MIN_RISE"), "b": C("symmetry", "MIN_TWIST"), "c": C("symmetry", "MIN_HELIX_COPIES"),
      "d": C("symmetry", "MIN_HELIX_COPIES_ASKED"), "e": C("symmetry", "_PSEUDO_21")}),
    ("cli", r"focus: chains, elements", "#a, #b", {"a": C("cli", "FOCUS_CHAINS"), "b": C("cli", "FOCUS_ELEMENTS")}),
    ("cli", r"dyad in the page; dyad turns; turn gain", r"#a:deg$^\circ$; #b:turns; #c",
     {"a": C("cli", "DYAD_IN_PAGE"), "b": C("cli", "DYAD_TURNS"), "c": C("cli", "DYAD_TURN_GAIN")}),
    ("cli", r"rosette orders (\texttt{cyclic=rosette})", "#a:list", {"a": C("cli", "ROSETTE_ORDERS")}),
    ("elstyle", r"width factor (\texttt{--width})", "#a:g--#b:g",
     {"a": C("elstyle", "WIDTH_MIN"), "b": C("elstyle", "WIDTH_MAX")}),
    ("layout", r"page units per \AA; strand pitch; margin", "#a; #b; #c",
     {"a": C("layout", "SCALE"), "b": C("layout", "PITCH"), "c": C("layout", "MARGIN")}),
    ("layout", r"rise: strand, $\alpha$, $3_{10}$, DNA (least)",
     r"#a, #b, #c, #d (#e)\angst{} per residue or base pair",
     {"a": C("layout", "STRAND_RISE"), "b": C("layout", "HELIX_RISE"), "c": C("layout", "TEN_RISE"),
      "d": C("layout", "DNA_RISE"), "e": C("layout", "_MIN_RISE")}),
    ("layout", r"home $k$; loop weight, slack; disulfide weight, slack", "#a; #b, #c; #d, #e",
     {"a": C("layout", "_K_HOME"), "b": C("layout", "_W_SEQ"), "c": C("layout", "_LOOP_SLACK"),
      "d": C("layout", "_W_SS"), "e": C("layout", "_SS_SLACK")}),
    ("layout", r"contact weight, cap; duplex joint weight", "#a, #b pairs; #c",
     {"a": C("layout", "_W_CONTACT"), "b": C("layout", "_CONTACT_CAP"), "c": C("layout", "_W_JOINT")}),
    ("layout", r"overlap $\lambda$ stages", "#a:list", {"a": C("layout", "_LAMBDAS")}),
    ("layout", r"L-BFGS history, iterations, gradient tolerance; Wolfe $c_1$, $c_2$", "#a, #b, #c:pow10; #d:pow10, #e",
     {"a": D("layout", "_lbfgs", "m"), "b": D("layout", "_lbfgs", "maxiter"), "c": D("layout", "_lbfgs", "pgtol"),
      "d": D("layout", "_wolfe", "c1"), "e": D("layout", "_wolfe", "c2")}),
    ("layout", r"grid", "#a:inv10 page units", {"a": C("layout", "_GRID")}),
    ("layout", r"push sweeps (items, orbits); extra push; slide step", "#a:,, #b:,; #c; #d",
     {"a": C("layout", "_SEP_SWEEPS"), "b": C("layout", "_SEP_ORBIT_SWEEPS"), "c": C("layout", "_SEP_PUSH"),
      "d": C("layout", "_SLIDE_STEP")}),
    ("layout", r"area kept for copies alike (translates, mirror images)", "#a:.2f, #b:.2f",
     {"a": C("layout", "_SEP_SYM_ROOM"), "b": C("layout", "_SEP_MIRROR_ROOM")}),
    ("layout", r"$\beta\alpha$ row: strands; share of units", "#a--#b; #c:of",
     {"a": Inline("layout", "can_be_row", r"return (\d+) <= len\(sheet\.strands\) <= _BA_MAX"),
      "b": C("layout", "_BA_MAX"),
      "c": C("layout", "_BA_SHARE")}),
    ("layout", r"room for stems (no protein)", "#a times their area", {"a": C("layout", "_STEM_ROOM")}),
    ("layout", r"page: width, text floor, longest side; PNG resolution", "#a in, #b:g pt, #c:g in; #d dpi",
     {"a": C("layout", "MAX_WIDTH_IN"), "b": C("layout", "MIN_TEXT_PT"), "c": C("layout", "MAX_PAGE_IN"),
      "d": C("render", "PNG_DPI")}),
    ("route", r"grid; box growth", "#a; #b page units", {"a": C("route", "RES"), "b": C("route", "_INFLATE")}),
    ("route", r"costs: step, turn, reuse, slab cell, slab run", "#a, #b:g, #c:g, #d:g, #e:g",
     {"a": Inline("route", "_astar", r"cost = g \+ (\d+) \+ \(_TURN"), "b": C("route", "_TURN"),
      "c": C("route", "_REUSE"),
      "d": C("route", "_SLAB_COST"), "e": C("route", "_SLAB_RUN")}),
    ("route", r"node limits", "#a:,; #b:,; #c:,",
     {"a": C("route", "_MAX_NODES"), "b": C("route", "_MAX_NODES_WIDE"), "c": C("route", "_MAX_NODES_PAGE")}),
    ("nucleic", r"H-bond; stretched; C1$'$--C1$'$; planes", r"#a\angst; #b\angst; #c:range\angst; #d:g$^\circ$",
     {"a": C("nucleic", "HBOND_MAX"), "b": C("nucleic", "HBOND_STRETCH"), "c": C("nucleic", "C1_PAIR"),
      "d": C("nucleic", "PLANE_MAX")}),
    ("nucleic", r"hairpin loop; axis window; segment turn; shortest piece", r"#a nt; #b bp; #c:g$^\circ$; #d bp",
     {"a": C("nucleic", "MIN_LOOP"), "b": C("nucleic", "AXIS_WINDOW"), "c": C("nucleic", "SEGMENT_TURN"),
      "d": C("nucleic", "PIECE_MIN")}),
    ("nucleic", r"G-tetrad H-bond", r"#a\angst", {"a": C("nucleic", "TETRAD_HBOND")}),
    ("model", r"crop large nucleic acids above", "#a nt", {"a": C("model", "AUTO_NT")}),
    ("membrane", r"helical slab half-thickness; barrel slab half-thickness", r"#a\angst; #b\angst",
     {"a": C("membrane", "HALF"), "b": C("membrane", "BARREL_HALF")}),
    ("membrane", r"helical slab: hydrophobic, polar (many crossing)",
     r"$\ge #a$, $\le #b$ ($\le #c:.2f$ with $\ge #d$)",
     {"a": C("membrane", "_MIN_HYDROPHOBIC"), "b": C("membrane", "_MAX_POLAR"),
      "c": C("membrane", "_MAX_POLAR_CROSSED"), "d": C("membrane", "_MANY_CROSSING")}),
    ("membrane", r"barrel: strands, length; outward hydrophobic, charged; girdle",
     r"#a, #b residues; $\ge #c$, $\le #d$; #e\angst",
     {"a": C("membrane", "_BARREL_STRANDS"), "b": C("membrane", "_BARREL_STRAND_LEN"),
      "c": C("membrane", "_MIN_OUTWARD_HYDROPHOBIC"), "d": C("membrane", "_MAX_OUTWARD_CHARGED"),
      "e": C("membrane", "_GIRDLE")}),
    ("membrane", r"long elements: helix, strand; crossing reach", r"#a, #b residues; #c\angst",
     {"a": Inline("membrane", "_long", r's\.kind == "H" and len\(s\) >= (\d+)'),
      "b": Inline("membrane", "_long", r's\.kind == "E" and len\(s\) >= (\d+)'), "c": C("membrane", "_MIN_CROSSING")}),
    ("membrane", r"face reach for a spanning element", r"#a\angst", {"a": C("membrane", "FACE_REACH")}),
    ("membrane", r"sides: ratio, difference (loops for barrels)", "#a, #b (#c)",
     {"a": C("membrane", "_SIDED_RATIO"), "b": C("membrane", "_SIDED_DIFF"), "c": C("membrane", "_SIDED_LOOPS")}),
    ("hbonds", r"donor--acceptor; angle (rotating H)", r"#a--#b\angst; #c:g$^\circ$ (#d:g$^\circ$)",
     {"a": C("hbonds", "MIN_DIST"), "b": C("hbonds", "MAX_DIST"), "c": C("hbonds", "SP2_ANGLE"),
      "d": C("hbonds", "ROTATING_ANGLE")}),
    ("domain\\_score", r"boundary tolerance; S5 targets (overlap, F1, recall)", "#a residues; #b:.2f, #c:.2f, #d:.2f",
     {"a": C("domain_score", "TOLERANCE"), "b": C("domain_score", "TARGETS", key="overlap"),
      "c": C("domain_score", "TARGETS", key="f1"), "d": C("domain_score", "TARGETS", key="recall")}),
]
# fmt: on

CAPTION = r"Thresholds and constants, as defined in the code at the release described (\NUM{release_tag})."
PLACEHOLDER = re.compile(r"#([a-z])(?::(deg|pct|range|list|of|times|pow10|inv10|turns|,|\.\d+f|g))?")


class Source:
    """The foldmap package's modules, parsed (never imported)."""

    def __init__(self, package: Path):
        self.package = package
        self._parsed: dict[str, tuple[str, ast.Module]] = {}

    def module(self, name: str) -> tuple[str, ast.Module]:
        if name not in self._parsed:
            text = (self.package / f"{name}.py").read_text(encoding="utf-8")
            self._parsed[name] = (text, ast.parse(text))
        return self._parsed[name]

    def assignment(self, module: str, name: str) -> tuple[ast.expr, int]:
        """The expression assigned to a module-level name (also one of a tuple assignment `A, B = 1, 2`)."""
        _, tree = self.module(module)
        for node in tree.body:
            if isinstance(node, ast.Assign):
                targets, value = node.targets, node.value
            elif isinstance(node, ast.AnnAssign) and node.value is not None:
                targets, value = [node.target], node.value
            else:
                continue
            for t in targets:
                if isinstance(t, ast.Name) and t.id == name:
                    return value, node.lineno
                if isinstance(t, ast.Tuple) and isinstance(value, ast.Tuple):
                    for el, v in zip(t.elts, value.elts):
                        if isinstance(el, ast.Name) and el.id == name:
                            return v, node.lineno
        raise LookupError(f"{module}.{name}: no module-level assignment")

    def evaluate(self, module: str, node: ast.expr):
        """Literals, tuples, lists and dicts of them, + - * / and unary -, names of the module's other constants,
        and np.radians / math.radians."""
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, (ast.Tuple, ast.List)):
            return tuple(self.evaluate(module, e) for e in node.elts)
        if isinstance(node, ast.Dict):
            return {self.evaluate(module, k): self.evaluate(module, v) for k, v in zip(node.keys, node.values)}
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -self.evaluate(module, node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            a, b = self.evaluate(module, node.left), self.evaluate(module, node.right)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            return a / b
        if isinstance(node, ast.Name):
            return self.evaluate(module, self.assignment(module, node.id)[0])
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in ("radians", "deg2rad")
            and len(node.args) == 1
        ):
            return math.radians(self.evaluate(module, node.args[0]))
        raise ValueError(f"{module}: cannot evaluate {ast.unparse(node)!r} without running the code")

    def read(self, ref) -> tuple[object, str, int]:
        """A ref's value, what it names and the line it is on."""
        if isinstance(ref, C):
            node, line = self.assignment(ref.module, ref.name)
            where = f"{ref.module}.{ref.name}"
            if ref.operands:
                if not isinstance(node, ast.BinOp):
                    raise ValueError(f"{where}: not a product of literals")
                return (self.evaluate(ref.module, node.left), self.evaluate(ref.module, node.right)), where, line
            value = self.evaluate(ref.module, node)
            if ref.key is not None:
                return value[ref.key], f"{where}[{ref.key!r}]", line
            return value, where, line
        if isinstance(ref, D):
            fn = self._function(ref.module, ref.function)
            args = fn.args.posonlyargs + fn.args.args
            defaults = dict(zip([a.arg for a in args][len(args) - len(fn.args.defaults) :], fn.args.defaults))
            defaults |= {a.arg: d for a, d in zip(fn.args.kwonlyargs, fn.args.kw_defaults) if d is not None}
            if ref.arg not in defaults:
                raise LookupError(f"{ref.module}.{ref.function}: no default for {ref.arg}")
            return self.evaluate(ref.module, defaults[ref.arg]), f"{ref.module}.{ref.function}({ref.arg}=)", fn.lineno
        if isinstance(ref, Inline):
            fn = self._function(ref.module, ref.function)
            text = ast.get_source_segment(self.module(ref.module)[0], fn) or ""
            found = re.findall(ref.pattern, text)
            if len(found) != 1:
                raise LookupError(f"{ref.module}.{ref.function}: {ref.pattern!r} matches {len(found)} times, not once")
            number = found[0]
            return (float(number) if "." in number else int(number)), f"{ref.module}.{ref.function} (inline)", fn.lineno
        raise TypeError(ref)

    def _function(self, module: str, name: str) -> ast.FunctionDef:
        _, tree = self.module(module)
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == name:
                return node
        raise LookupError(f"{module}.{name}: no such function")


def _g(x) -> str:
    return f"{x:g}" if isinstance(x, float) else str(x)


def _pow10(x: float) -> str:
    n = round(math.log10(x))
    if not math.isclose(x, 10.0**n):
        raise ValueError(f"{x} is not a power of ten")
    return f"$10^{{{n}}}$"


def fmt(value, spec: str | None) -> str:
    if spec is None:
        return str(value)
    if spec == "deg":
        return _g(round(math.degrees(value), 6))
    if spec == "pct":
        return _g(round(value * 100, 6)) + r"\%"
    if spec == "range":
        return f"{_g(value[0])}--{_g(value[1])}"
    if spec == "list":
        return ", ".join(_g(v) for v in value)
    if spec == "of":
        return f"{value[0]} of {value[1]}"
    if spec == "times":
        return rf"${_g(value[0])}\times{_g(value[1])}$"
    if spec == "pow10":
        return _pow10(value)
    if spec == "inv10":
        return _pow10(1 / value)
    if spec == "turns":
        steps = sorted({abs(v) for v in value if v})
        if sorted(value, key=abs)[0] != 0 or any(b - a != steps[0] for a, b in zip([0, *steps], steps)):
            raise ValueError(f"turns {value}: not 0 to ±max in equal steps")
        if sorted(value) != sorted([0, *steps, *(-s for s in steps)]):
            raise ValueError(f"turns {value}: not symmetric")
        return rf"0 to $\pm {_g(steps[-1])}^\circ$ by {_g(steps[0])}$^\circ$"
    return format(value, spec)


def build(source: Source) -> tuple[list[str], list[dict]]:
    """The table's rows as LaTeX, and every value with what it names."""
    rows, values = [], []
    for module, quantity, template, refs in ROWS:
        read = {k: source.read(r) for k, r in refs.items()}
        used = {m.group(1) for m in PLACEHOLDER.finditer(template)}
        if used != set(refs):
            raise ValueError(f"{module} {quantity}: refs {sorted(set(refs) ^ used)} unused or undefined")
        cell = PLACEHOLDER.sub(lambda m, read=read: fmt(read[m.group(1)][0], m.group(2)), template)
        rows.append(f"{module} & {quantity} & {cell} \\\\")
        for k in refs:
            value, where, line = read[k]
            values.append({"row": f"{module}: {quantity}", "names": where, "line": line, "value": value})
    return rows, values


def longtable(rows: list[str], commit: str | None) -> str:
    head = [
        f"% Table S1, generated by scripts/paper/table_s1.py from the code at {commit or 'the given sources'}; "
        "do not edit by hand",
        r"\begin{longtable}{@{}p{3.2cm}p{8.4cm}p{3.6cm}@{}}",
        rf"\caption{{{CAPTION}}}\label{{tab:params}}\\",
        r"\toprule",
        r"Module & Quantity & Value \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule Module & Quantity & Value \\ \midrule \endhead",
        r"\bottomrule \endfoot",
    ]
    return "\n".join([*head, *rows, r"\end{longtable}"]) + "\n"


def typed_rows(tex: str) -> list[tuple[str, str]]:
    """(module, value) of each row of tab:params in a .tex file, in order."""
    body = tex.split(r"\label{tab:params}", 1)[1].split(r"\end{longtable}", 1)[0]
    body = body.split(r"\bottomrule \endfoot", 1)[1]
    out = []
    for line in body.strip().splitlines():
        cells = [c.strip() for c in line.rstrip().removesuffix("\\\\").split(" & ")]
        if len(cells) == 3:
            out.append((cells[0], cells[2]))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", type=Path, default=TREE / "src" / "foldmap", help="the foldmap package to read")
    ap.add_argument("--out", type=Path, default=TABLES, help="where params.tex and params.json go")
    ap.add_argument("--check-against", type=Path, help="a .tex holding a typed tab:params to compare with")
    args = ap.parse_args(argv)
    commit = None
    prov = TREE.parent / "provenance.json"
    if args.src.resolve() == (TREE / "src" / "foldmap").resolve() and prov.exists():
        commit = json.loads(prov.read_text(encoding="utf-8")).get("harness_commit_short")
    rows, values = build(Source(args.src))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "params.tex").write_text(longtable(rows, commit), encoding="utf-8")
    doc = {"generator": "scripts/paper/table_s1.py", "commit": commit, "package": str(args.src), "values": values}
    (args.out / "params.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8")
    print(f"{len(rows)} rows, {len(values)} values -> {args.out / 'params.tex'}")
    if args.check_against:
        typed = typed_rows(args.check_against.read_text(encoding="utf-8"))
        made = [(r.split(" & ")[0], r.split(" & ")[2].removesuffix(" \\\\")) for r in rows]
        diffs = [(i, t, m) for i, (t, m) in enumerate(zip(typed, made)) if t != m]
        if len(typed) != len(made):
            print(f"typed table: {len(typed)} rows, generated: {len(made)}")
        for i, t, m in diffs:
            print(f"row {i + 1} {m[0]}:\n  typed     {t[1]}\n  generated {m[1]}")
        print(f"{len(diffs)} of {min(len(typed), len(made))} rows differ")
        return 1 if diffs or len(typed) != len(made) else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
