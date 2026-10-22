# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The researcher's own data, read against the structure (plan 13d): per-residue values (--data FILE.csv) and
crosslinks (--crosslinks FILE.csv).

A data file is a table (comma, tab or semicolon separated, the first line not starting # its header) with a
chain column, a residue column (the number as the structure file writes it, with its insertion code: 45, -7, 100A;
or the code in a column of its own, icode) and value columns:

    chain,residue,uptake,site
    A,45,0.82,core
    A,100A,0.15,surface

A column is numbers when every non-empty cell is a number (empty, NA, N/A and NaN cells are empty), else
categories. Without --data-column the first column of numbers is read. An optional residue-name column (resname,
aa: ALA or A) is checked against the structure, so a file numbered another way (label_seq_id, a construct's own
numbering) is caught where the names differ instead of colouring the wrong residue.

A crosslinks file has chain1, residue1, chain2, residue2 (or chain_a ... residue_b) and an optional distance: the
crosslinker's span, the most Cα–Cα distance (Å) it bridges (CROSSLINK_SPAN, 30 Å, when not given). Each crosslink is
measured on the structure (crosslink_lengths); one longer than its span is a warning, and is drawn all the same.

Every row is checked (spec.Problem, as check_references' pass 2 checks a reference): a chain the structure lacks, a
residue its chain lacks, a residue named otherwise than in the structure, a cell that is no residue number, a
residue given twice, each named with its row (the file's line number, the header being row 1); these are warnings,
the row is not used. A file whose header lacks what is needed, or none of whose rows names a residue of the
structure, is an error. A chain the structure holds only as assembly copies (A1, A2: io.read_model names copies so)
is read for every copy, and a warning says so. Deterministic: the same file and structure give the same result."""

from __future__ import annotations

import csv
import io
import math
import re
from dataclasses import dataclass
from pathlib import Path

from .spec import Problem

CHAIN = ("chain", "chainid", "chainname", "authasymid")
RESIDUE = ("residue", "residuenumber", "resnum", "resno", "resi", "resid", "res", "resseq", "authseqid", "number")
ICODE = ("icode", "insertioncode", "inscode", "pdbxpdbinscode")
RESNAME = ("resname", "residuename", "aa", "aminoacid", "authcompid")
DISTANCE = ("distance", "dist", "length", "span", "maxdistance", "cutoff")  # the crosslinker's span (Å)
# The Cα–Cα distance a crosslink is checked against when its row gives none: 30 Å, the usual limit for the
# lysine-reactive NHS esters (DSS, BS3; their 11.4 Å spacer, two lysine side chains and backbone motion: Merkley et
# al. 2014, Protein Sci 23:747; Leitner et al. 2016, Trends Biochem Sci 41:20). A row's own distance (the
# crosslinker's span, Å) replaces it, e.g. 20 Å for a zero-length (EDC) or short crosslinker.
CROSSLINK_SPAN = 30.0
EMPTY = {"", "na", "n/a", "nan", "none", "null"}  # cells that hold no value (any case)
_NUMBER = re.compile(r"(-?\d+)\s*([A-Za-z]?)")


@dataclass(frozen=True)
class ResidueData:
    """Per-residue values read from the researcher's file: the column read, its kind ("number" or "category") and
    the value of each residue, keyed by (chain, residue tag with its insertion code) as the structure names them
    (an assembly copy under its own name, A1)."""

    column: str
    kind: str  # "number" | "category"
    values: dict[tuple[str, str], float | str]


Crosslink = tuple[tuple[str, str], tuple[str, str], float | None]


def _key(name: str) -> str:
    """A column name as matched: lower case, without spaces, _, - and dots."""
    return re.sub(r"[\s_.\-]", "", name).lower()


def _find(header: list[str], names, suffix: str = "") -> int | None:
    """The first column of `header` named one of `names` (+ suffix), or None."""
    wanted = {n + suffix for n in names}
    return next((k for k, h in enumerate(header) if _key(h) in wanted), None)


def _rows(path, where: str) -> tuple[list[str], list[tuple[int, list[str]]], list[Problem]]:
    """(header, [(row number, cells)], problems) of a table: leading lines starting # and blank lines skipped,
    the delimiter (comma, tab or semicolon) taken from the header line. A row's number is the file's line number."""
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except FileNotFoundError:
        return [], [], [Problem(where, "no such file", True)]
    except UnicodeDecodeError:
        return [], [], [Problem(where, "not a text table (UTF-8): save it as CSV", True)]
    except OSError as err:  # a directory, no permission: a message, never a traceback
        return [], [], [Problem(where, f"cannot be read ({type(err).__name__})", True)]
    try:
        return _parse(text, where)
    except csv.Error as err:  # an unbalanced quote, a field past the csv module's limit
        return [], [], [Problem(where, f"not a readable table ({err})", True)]


MAX_COLUMNS = 64  # a data or crosslink table has a few columns; more is refused before any row is kept
MAX_ROWS = 100_000  # rows kept per table (a protein has far fewer residues; a crafted file must not fill memory)


def _parse(text: str, where: str) -> tuple[list[str], list[tuple[int, list[str]]], list[Problem]]:
    """_rows' parse of the text, its size bounded: header width and row count capped before rows are kept, so memory
    stays proportional to the cells actually written (a short row is padded only up to MAX_COLUMNS)."""
    lines = text.splitlines()
    head = next((k for k, line in enumerate(lines) if line.strip() and not line.lstrip().startswith("#")), None)
    if head is None:
        return [], [], [Problem(where, "the file is empty; its first line is a header (chain,residue,...)", True)]
    first = lines[head]
    delimiter = "\t" if "\t" in first and "," not in first else ";" if ";" in first and "," not in first else ","
    reader = csv.reader(io.StringIO("\n".join(lines[head:])), delimiter=delimiter)
    header = [h.strip() for h in next(reader)]
    if len(header) > MAX_COLUMNS:
        return [], [], [Problem(where, f"has {len(header)} columns; a table has at most {MAX_COLUMNS}", True)]
    rows = []
    for cells in reader:
        number = head + reader.line_num
        if not any(c.strip() for c in cells) or cells[0].lstrip().startswith("#"):
            continue
        if len(rows) == MAX_ROWS:
            return [], [], [Problem(where, f"has more than {MAX_ROWS:,} rows; split it or keep one chain", True)]
        cells = [c.strip() for c in cells[:MAX_COLUMNS]]  # a row longer than the header keeps its cells, capped
        rows.append((number, cells + [""] * (len(header) - len(cells))))
    return header, rows, []


class _Residues:
    """The protein residues of a backbone by chain, for matching rows: numbers, names, and the chains each chain
    named in a file stands for (itself, or its assembly copies)."""

    def __init__(self, bb, na_chains=()):
        self.names: dict[str, dict[str, str]] = {}  # chain -> residue tag -> residue name
        self.span: dict[str, tuple[int, int]] = {}
        for label in bb.labels:
            self.names.setdefault(label.chain, {})[label.tag] = label.name
            lo, hi = self.span.get(label.chain, (label.seq, label.seq))
            self.span[label.chain] = (min(lo, label.seq), max(hi, label.seq))
        self.na_chains = set(na_chains) - set(self.names)

    def copies(self, chain: str, named: set[str]) -> list[str]:
        """The structure's chains a file's `chain` stands for: itself when the structure has it, else the assembly
        copies named after it (A -> A1, A2; io.read_model names copies so), each copy going to the longest chain of
        the file (`named`) its name is followed by digits after (figcheck._file_chain's rule)."""
        if chain in self.names:
            return [chain]

        def deposited(c: str) -> str | None:
            found = [n for n in named if c.startswith(n) and c[len(n) :].isdigit()]
            return max(found, key=len) if found else None

        return [c for c in self.names if c not in named and deposited(c) == chain]

    def no_chain(self, chain: str) -> str:
        if chain in self.na_chains:
            return f"chain {chain} is a nucleic acid; only protein residues are read"
        return f"no protein chain {chain} in the structure; its protein chains are {', '.join(self.names) or 'none'}"

    def no_residue(self, chain: str, tag: str) -> str:
        lo, hi = self.span[chain]
        return f"chain {chain} has no residue {tag} (the structure file numbers its residues {lo} to {hi})"


def _tag(cell: str, icode: str = "") -> str | None:
    """A residue tag (45, -7, 100A) from a number cell and an optional insertion-code cell (empty, . or ?: none);
    None when it is none."""
    m = _NUMBER.fullmatch(cell.strip())
    icode = "" if icode.strip() in (".", "?") else icode.strip()
    if m is None or (m[2] and icode) or len(icode) > 1:
        return None
    return f"{int(m[1])}{m[2] or icode}"


def _same_residue(given: str, name: str) -> bool:
    """A residue name cell (ALA or A, any case) names the structure's residue `name`."""
    from .sequence import one_letter

    given = given.strip().upper()
    return given == name.upper() if len(given) != 1 else given == one_letter([name])


class _Reader:
    """Rows matched to residues: each end (chain cell, residue cell, ...) to the structure's (chain, tag) pairs, the
    problems said once per row; the copies read for each chain said once. `named`: every chain the file names;
    `single`: the chain read when the file has no chain column (the structure's one protein chain)."""

    def __init__(self, where: str, residues: _Residues, named: set[str]):
        self.where, self.res, self.named = where, residues, named
        self.single = next(iter(residues.names)) if len(residues.names) == 1 else None
        self.problems: list[Problem] = []
        self.copied: dict[str, list[str]] = {}
        self.suffix: dict[str, str] = {}  # assembly copy -> the digits after its chain's name (A2 -> 2)

    def row(self, number: int, text: str) -> None:
        self.problems.append(Problem(f"{self.where} row {number}", text + "; the row is not used", False))

    def end(self, number: int, chain: str | None, cell: str, icode: str = "", resname: str = ""):
        """[(chain, tag)] the row's residue is in the structure (one per copy), or None (said for the row)."""
        if chain is None:
            chain = self.single
        if not chain:
            self.row(number, "no chain given")
            return None
        tag = _tag(cell, icode)
        if tag is None:
            self.row(number, f"{cell!r} is no residue number (45, -7 or 100A)")
            return None
        chains = self.res.copies(chain, self.named)
        if not chains:
            self.row(number, self.res.no_chain(chain))
            return None
        if chains != [chain]:
            self.copied.setdefault(chain, chains)
            self.suffix.update({c: c[len(chain) :] for c in chains})
        found = []
        for c in chains:
            name = self.res.names[c].get(tag)
            if name is None:
                self.row(number, self.res.no_residue(c, tag))
                return None
            if resname and not _same_residue(resname, name):
                self.row(number, f"{c} {tag} is {name} in the structure, the row says {resname}")
                return None
            found.append((c, tag))
        return found

    def notes(self) -> list[Problem]:
        """One warning per chain read for its assembly copies."""
        return [
            Problem(
                self.where,
                f"chain {chain} is drawn as its assembly copies {', '.join(copies)}; each copy takes chain {chain}'s "
                "rows",
                False,
            )
            for chain, copies in sorted(self.copied.items())
        ]


def _header_problem(where: str, header: list[str], text: str) -> Problem:
    return Problem(where, f"{text}; its columns are {', '.join(header) or 'none'}", True)


def _number(cell: str) -> float | None:
    """The number a cell holds, None for an empty one; ValueError when it holds something else."""
    if cell.strip().lower() in EMPTY:
        return None
    return float(cell)


def _is_number(cell: str) -> bool:
    try:
        _number(cell)
    except ValueError:
        return False
    return True


def read_residue_data(path, column: str | None, bb, *, na_chains=()) -> tuple[ResidueData, list[Problem]]:
    """The values of one column of the researcher's per-residue table (see the module text) on the residues of
    `bb` (the protein backbone as read, every chain, before --chains). `column`: the column to read (exact, else
    any case); None reads the first column of numbers, a warning naming the others. `na_chains`: the nucleic-acid
    chains read, named as such when a row names one. The problems: one per row not used (with its row number), the
    assembly copies read, and fatal ones (nothing read) when the header lacks a chain (with several protein chains),
    residue or value column, or no row gives a value on a residue of the structure."""
    where = f"--data {path}"
    header, rows, problems = _rows(path, where)
    if problems:
        return ResidueData(column or "", "number", {}), problems
    chain_k, residue_k = _find(header, CHAIN), _find(header, RESIDUE)
    icode_k, resname_k = _find(header, ICODE), _find(header, RESNAME)
    residues = _Residues(bb, na_chains)
    empty = ResidueData(column or "", "number", {})
    if residue_k is None:
        return empty, [_header_problem(where, header, "no residue column (residue, resnum or resi)")]
    if chain_k is None and len(residues.names) != 1:
        text = f"no chain column, and the structure has protein chains {', '.join(residues.names) or 'none'}"
        return empty, [_header_problem(where, header, text)]
    ids = {chain_k, residue_k, icode_k, resname_k}
    candidates = [k for k in range(len(header)) if k not in ids]
    numeric = [k for k in candidates if all(_is_number(cells[k]) for _, cells in rows)]
    if column is not None:
        value_k = next((k for k in candidates if header[k] == column), None)
        if value_k is None:
            value_k = next((k for k in candidates if header[k].lower() == column.lower()), None)
        if value_k is None:
            return empty, [_header_problem(where, header, f"no value column {column!r}")]
    elif not numeric:
        text = "no column of numbers; name the column to read with --data-column NAME"
        return empty, [_header_problem(where, header, text)]
    else:
        value_k = numeric[0]
        if len(numeric) > 1:
            others = ", ".join(header[k] for k in numeric[1:])
            text = (
                f"read column {header[value_k]}, the first of numbers; the file also has {others} (--data-column NAME)"
            )
            problems.append(Problem(where, text, False))
    kind = "number" if value_k in numeric else "category"
    reader = _Reader(where, residues, {cells[chain_k] for _, cells in rows if chain_k is not None and cells[chain_k]})
    values: dict[tuple[str, str], float | str] = {}
    first_row: dict[tuple[str, str], int] = {}
    on_residue = False
    for number, cells in rows:
        found = reader.end(
            number,
            cells[chain_k] if chain_k is not None else None,
            cells[residue_k],
            cells[icode_k] if icode_k is not None else "",
            cells[resname_k] if resname_k is not None else "",
        )
        if found is None:
            continue
        on_residue = True
        twice = [first_row[r] for r in found if r in first_row]
        if twice:
            reader.row(number, f"{found[0][0]} {found[0][1]} is given by row {twice[0]} already")
            continue
        cell = cells[value_k]
        if kind == "number":
            value = _number(cell)
            if value is not None and not math.isfinite(value):
                reader.row(number, f"{header[value_k]} {cell!r} is not a finite number")
                continue
        else:
            value = None if cell.lower() in EMPTY else cell
        for r in found:
            first_row[r] = number
            if value is not None:
                values[r] = value
    problems += reader.problems + reader.notes()
    if not on_residue:
        text = "no row names a residue of the structure (residue numbers as the structure file writes them)"
        problems.append(Problem(where, text, True))
    elif not values:
        problems.append(Problem(where, f"no row gives a value in column {header[value_k]}", True))
    return ResidueData(header[value_k], kind, dict(sorted(values.items(), key=_order(bb)))), problems


def _order(bb):
    """Sort key for (chain, tag) keys: the backbone's residue order."""
    rank = {(label.chain, label.tag): k for k, label in enumerate(bb.labels)}
    return lambda item: rank[item[0]]


def read_crosslinks(path, bb, *, na_chains=()) -> tuple[list[Crosslink], list[Problem]]:
    """The crosslinks of the researcher's table (columns chain1, residue1, chain2, residue2, or _a/_b; optional
    insertion-code and residue-name columns per end, and distance) on the residues of `bb`, each ((chain, tag),
    (chain, tag), distance or None), in file order. A chain read as its assembly copies (A -> A1, A2): a crosslink
    between two copied chains joins the copies of one number (A1-B1, A2-B2), one with a chain the structure has as
    named joins it to every copy; a warning says which chains were read so. The problems as read_residue_data's."""
    where = f"--crosslinks {path}"
    header, rows, problems = _rows(path, where)
    if problems:
        return [], problems
    residues = _Residues(bb, na_chains)
    ends = []
    for one, other in (("1", "a"), ("2", "b")):
        cols = []
        for names in (CHAIN, RESIDUE, ICODE, RESNAME):
            k = _find(header, names, one)
            cols.append(k if k is not None else _find(header, names, other))
        ends.append(cols)
    if ends[0][1] is None or ends[1][1] is None:
        return [], [_header_problem(where, header, "no residue1 and residue2 columns")]
    if (ends[0][0] is None or ends[1][0] is None) and len(residues.names) != 1:
        text = f"no chain1 and chain2 columns, and the structure has protein chains {', '.join(residues.names)}"
        return [], [_header_problem(where, header, text)]
    distance_k = _find(header, DISTANCE)
    named = {cells[e[0]] for _, cells in rows for e in ends if e[0] is not None and cells[e[0]]}
    reader = _Reader(where, residues, named)
    out: list[Crosslink] = []
    on_residue = False
    for number, cells in rows:
        found = []
        for chain_k, residue_k, icode_k, resname_k in ends:
            found.append(
                reader.end(
                    number,
                    cells[chain_k] if chain_k is not None else None,
                    cells[residue_k],
                    cells[icode_k] if icode_k is not None else "",
                    cells[resname_k] if resname_k is not None else "",
                )
            )
            if found[-1] is None:
                break
        if None in found:
            continue
        on_residue = True
        distance = None
        if distance_k is not None:
            try:
                distance = _number(cells[distance_k])
            except ValueError:
                distance = math.nan
            if distance is not None and not math.isfinite(distance):
                reader.problems.append(
                    Problem(
                        f"{where} row {number}",
                        f"distance {cells[distance_k]!r} is no number; the crosslink is read without one",
                        False,
                    )
                )
                distance = None
        pairs = _pairs(*found, reader.suffix)
        if not pairs:
            reader.row(number, "the assembly copies of its two chains share no number (A1-B1)")
        out += [(a, b, distance) for a, b in pairs]
    problems += reader.problems + reader.notes()
    if not on_residue:
        text = "no row names two residues of the structure (residue numbers as the structure file writes them)"
        problems.append(Problem(where, text, True))
    return out, problems


def _pairs(a: list[tuple[str, str]], b: list[tuple[str, str]], suffix: dict[str, str]) -> list:
    """The residue pairs a crosslink joins, its ends a and b each one residue or one per assembly copy (`suffix`:
    copy -> the digits after its chain's name): both ends copies, the copies of one number (A1-B1, A2-B2); else each
    residue of one end to each of the other's."""
    if all(c in suffix for c, _ in a + b):
        by = {suffix[c]: (c, t) for c, t in b}
        return [(x, by[suffix[x[0]]]) for x in a if suffix[x[0]] in by]
    return [(x, y) for x in a for y in b]


def read_user_data(data, column: str | None, crosslinks, bb, *, na_chains=()):
    """(ResidueData or None, crosslinks or None, problems) for a figure's --data FILE [--data-column NAME] and
    --crosslinks FILE (each None when not given), read on `bb` (every protein chain read, before --chains). The
    problems of both files, in order; a fatal one means that file could not be read (cli.make_layout refuses the
    figure)."""
    values, links, problems = None, None, []
    if data is not None:
        values, found = read_residue_data(data, column, bb, na_chains=na_chains)
        problems += found
    if crosslinks is not None:
        links, found = read_crosslinks(crosslinks, bb, na_chains=na_chains)
        problems += found
    return values, links, problems


def summary(values: ResidueData | None, links: list | None, data=None, crosslinks=None) -> list[str]:
    """What was read, one line per file: '--data FILE: column C (numbers) on N residues', '--crosslinks FILE: N
    crosslinks'."""
    out = []
    if values is not None:
        kind = "numbers" if values.kind == "number" else "categories"
        out.append(f"--data {data}: column {values.column} ({kind}) on {len(values.values)} residues")
    if links is not None:
        out.append(f"--crosslinks {crosslinks}: {len(links)} crosslinks")
    return out


def crosslink_lengths(links: list, bb) -> list[float]:
    """The Cα–Cα distance (Å, rounded to 0.001) of each crosslink, measured on the structure as read (`bb`, the
    backbone the files were read on), in the order of `links`."""
    import numpy as np

    where = {(label.chain, label.tag): k for k, label in enumerate(bb.labels)}
    out = []
    for a, b, _ in links:
        d = float(np.linalg.norm(bb.xyz[where[a], 1] - bb.xyz[where[b], 1]))
        out.append(round(d, 3))
    return out


def crosslink_name(link) -> str:
    """A crosslink as the notes name it: A:45 – B:120."""
    (ca, ta), (cb, tb), _ = link
    return f"{ca}:{ta} – {cb}:{tb}"


def over_span(link, length: float) -> bool:
    """Whether a crosslink's Cα–Cα distance is longer than its span (its row's distance, else CROSSLINK_SPAN)."""
    return length > (link[2] if link[2] is not None else CROSSLINK_SPAN)


def crosslink_problems(links: list, lengths: list[float], crosslinks, drawn=None) -> list:
    """A warning for each crosslink whose ends lie further apart (Cα–Cα, in the structure) than its crosslinker
    spans: drawn all the same, dashed (when both its chains are among `drawn`, None: all), so the figure shows what
    the structure cannot explain."""
    out = []
    for link, length in zip(links, lengths):
        if over_span(link, length):
            span = link[2] if link[2] is not None else CROSSLINK_SPAN
            text = f"crosslink {crosslink_name(link)}: Cα–Cα {length:.1f} Å, longer than its {span:g} Å span"
            if drawn is None or (link[0][0] in drawn and link[1][0] in drawn):
                text += " (drawn dashed)"
            out.append(Problem(f"--crosslinks {crosslinks}", text, False))
    return out
