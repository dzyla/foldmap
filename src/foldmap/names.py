# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Conventional names (--names): the helix ids the file itself gives (`deposited`: globin A-H), strands numbered
β1…βn with each βα-row helix named after its strand (`beta`), or a YAML file naming elements and stems:

    elements:            # element ref -> name (refs as in --rename: label, chain:label, #k, res:a-b, res:chain:a-b)
      res:A:3-18: A
    stems:               # a nucleotide (chain:number) held by the stem -> name
      A:1: acceptor

Names are given after the layout is placed and labelled, before the user's own renames (layout.build_layout)."""

from __future__ import annotations

import datetime
import gzip
import re
from collections import Counter
from pathlib import Path

from .layout import AmbiguousName, Layout, provisional, resolve
from .model import SSE, Nucleic, copy_classes, group_by_chain, pair_places, res_ref
from .yamlfile import read_yaml

KEYWORDS = ("deposited", "beta")
_GENERIC = re.compile(r"(HELX_P)?\d+|H\d+")  # numbered helices name nothing
_PRIMES = ("", "′", "″", "‴")
_LETTERS = "abcdefghijklmnopqrstuvwxyz"  # beta names of the helices after a strand past what _PRIMES can mark
_LETTER = re.compile(r"[A-Za-z]")  # a deposited name starting with a letter, as strand labels do


def check_spec(spec: str) -> None:
    """A names spec is a keyword or an existing YAML file."""
    if spec not in KEYWORDS and not Path(spec).is_file():
        raise ValueError(f"--names {spec!r}: use deposited, beta or a YAML file (no such file)")


def load_names(
    spec: str,
    path,
    layout: Layout,
    sses: list[SSE],
    nucleic: Nucleic | None,
    *,
    whole: list[SSE] | None = None,
    chains: set[str] | None = None,
) -> tuple[dict[str, str], dict[int, str]]:
    """(element renames, ref -> name as --rename takes them; stem renames, index in layout.stems -> name) for one
    spec: `deposited`, `beta` or a YAML file. `path`: the structure file. `whole` (the elements of the chains read,
    before focus left some out) and `chains` (every protein chain of the file) let a file name elements this figure
    leaves out: those refs are skipped, any other ref that names nothing is an error."""
    drawn = [s for s in sses if s.id in layout.placed]
    if spec == "deposited":
        return deposited_names(path, drawn, whole), {}
    if spec == "beta":
        return beta_names(layout, drawn, whole), {}
    return _file_names(Path(spec), layout, drawn, nucleic, whole, chains)


def all_names(path, layout: Layout, sses: list[SSE], whole: list[SSE] | None = None) -> dict[str, dict[str, str]]:
    """Every drawn element's names in the other conventions (Q-A12), computed and never applied: {"dep": {element id:
    deposited name}, "beta": {element id: β-scheme name}}, as --names deposited and --names beta would name them
    (before primes tell repeats apart). layout.resolve reads them for dep: and beta: references (Layout.other_names).
    `path`: the structure file (None: no deposited names); `whole` as for load_names."""
    drawn = [s for s in sses if s.id in layout.placed]

    def ids(named: dict[str, str]) -> dict[str, str]:  # res:A:3-18 -> A:3-18, the element's id
        return {ref.removeprefix("res:"): name for ref, name in named.items()}

    return {
        "dep": ids(deposited_names(path, drawn, whole)) if path is not None else {},
        "beta": ids(beta_names(layout, drawn, whole)),
    }


# ---------------------------------------------------------------------------------------------- deposited


def deposited_names(path, sses: list[SSE], whole: list[SSE] | None = None) -> dict[str, str]:
    """Each helix (H or G) named after the deposited helix record (mmCIF _struct_conf, PDB HELIX) at its place, when
    that record's id is a name (helix_names): records and helices of a chain pair one-to-one as copies' elements do
    (model.pair_places). In a chain holding strands (in `whole`, the chain as read, when given) a name starting with
    a letter takes an α (CheY's A-E: αA-αE), as the strands keep their letters (the sheet records are numbered in
    every file checked)."""
    rows = _helix_rows(path)
    names = helix_names([(hid, chain) for hid, chain, _, _ in rows])
    named = [(chain, *sorted((lo, hi)), name) for (_, chain, lo, hi), name in zip(rows, names) if name]
    files = {chain for chain, _, _, _ in named}
    stranded = {s.chain for s in whole or sses if s.kind == "E"}
    out = {}
    for chain, mine in group_by_chain(s for s in sses if s.kind in ("H", "G")).items():
        records = [(lo, hi, name) for c, lo, hi, name in named if c == _deposited_chain(chain, files)]
        helices = [("H", *sorted((s.first.seq, s.last.seq))) for s in mine]
        for i, j in pair_places(helices, [("H", lo, hi) for lo, hi, _ in records]).items():
            name = records[j][2]
            out[res_ref(chain, mine[i].first, mine[i].last)] = (
                f"α{name}" if chain in stranded and _LETTER.match(name) else name
            )
    return out


def helix_names(rows: list[tuple[str, str]]) -> list[str | None]:
    """The name each helix record (id, chain) gives, or None. A chain's ids lose the chain id when every one of them
    starts with it and something is left (2HHB AA, AB -> A, B; 7TIM AH1 -> H1). Numbered ids (1, HELX_P1, H1; A1
    once the chain is gone), the wwPDB serial ids (AA1, AA2 … AA9, AB1 … across the whole file) and ids used twice
    in one chain name nothing."""
    if [hid for hid, _ in rows] == [_serial(k) for k in range(len(rows))]:
        return [None] * len(rows)
    ids: dict[str, list[str]] = {}
    for hid, chain in rows:
        ids.setdefault(chain, []).append(hid)
    strip = {c: all(h.startswith(c) and len(h) > len(c) for h in hs) for c, hs in ids.items()}
    out = []
    for hid, chain in rows:
        name = hid[len(chain) :] if strip[chain] else hid
        generic = _GENERIC.fullmatch(hid) or name.isdigit() or ids[chain].count(hid) > 1
        out.append(None if generic else name)
    return out


def _serial(k: int) -> str:
    """The k-th wwPDB serial helix id: AA1 … AA9, AB1 … AZ9, BA1 …"""
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return f"{letters[k // 234 % 26]}{letters[k // 9 % 26]}{k % 9 + 1}"


def _deposited_chain(chain: str, files: set[str]) -> str:
    """The deposited chain a drawn chain is (an assembly copy is named after it plus a number: A -> A1, A2)."""
    if chain in files:
        return chain
    found = [c for c in files if chain.startswith(c) and chain[len(c) :].isdigit()]
    return max(found, key=len) if found else chain


def _helix_rows(path) -> list[tuple[str, str, int, int]]:
    """(id, chain, first, last residue number) of each deposited helix, in file order."""
    import gemmi

    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except (RuntimeError, ValueError):  # not mmCIF: a PDB file's HELIX records
        return _pdb_helix_rows(path)
    table = block.find(
        "_struct_conf.",
        ["pdbx_PDB_helix_id", "beg_auth_asym_id", "beg_auth_seq_id", "end_auth_seq_id", "?conf_type_id"],
    )
    out = []
    for row in table:
        if row.has(4) and not row[4].startswith("HELX"):
            continue  # a turn record
        try:
            out.append((*map(gemmi.cif.as_string, (row[0], row[1])), int(row[2]), int(row[3])))
        except ValueError:
            continue
    return out


def _pdb_helix_rows(path) -> list[tuple[str, str, int, int]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = [line for line in f if line.startswith("HELIX ")]
    out = []
    for line in lines:
        try:
            out.append((line[11:14].strip(), line[19].strip(), int(line[21:25]), int(line[33:37])))
        except (ValueError, IndexError):
            continue
    return out


# ---------------------------------------------------------------------------------------------- beta


def beta_names(layout: Layout, sses: list[SSE], whole: list[SSE] | None = None) -> dict[str, str]:
    """The strands of each chain become β1…βn in sequence (numbered over `whole`, the chain as read, when focus
    left some out). In a chain with a βα row (layout.ba_rows), the row helix after βn becomes αn; its other helices
    are named after the strand before them, in sequence: αn while that name is free, then αn′, αn″ … (αN, αN′ …
    before the first strand). Without a row the helices keep their labels; 3-10 helices always do."""
    rows = {k for row in layout.ba_rows for k in row}
    drawn = {s.id for s in sses}
    out: dict[str, str] = {}
    for _, mine in group_by_chain(whole or sses).items():
        mine = sorted(mine, key=lambda s: s.start)
        strands = [s for s in mine if s.kind == "E"]
        number = {s.id: n for n, s in enumerate(strands, 1)}
        names = {s.id: f"β{number[s.id]}" for s in strands}
        if any(s.id in rows for s in strands):
            names.update(_row_helices(mine, strands, number, rows))
        out.update({res_ref(s.chain, s.first, s.last): names[s.id] for s in mine if s.id in names and s.id in drawn})
    return out


def _row_helices(mine: list[SSE], strands: list[SSE], number: dict[str, int], rows: set[str]) -> dict[str, str]:
    """The helices of a chain with a βα row: each row helix after βn is αn; the others, in sequence, αn, αn′ …
    after the last strand before them, the first name not taken. When a strand is followed by more of them than
    αn, αn′, αn″, αn‴ can name, they are lettered αna, αnb … instead (9YKJ: nine helices after β10 were α10′ to
    α10′′′′′′′′′, names whose width the loops could not route around)."""

    def after(h: SSE, pool: list[SSE]) -> int:  # the number of the pool's last strand before the helix (0: none)
        return max((number[s.id] for s in pool if s.start < h.start), default=0)

    helices = [s for s in mine if s.kind == "H"]
    in_row = [s for s in strands if s.id in rows]

    def name(n: int) -> str:  # before the first strand: αN
        return f"α{n or 'N'}"

    out = {h.id: name(after(h, in_row)) for h in helices if h.id in rows}
    taken = set(out.values())
    rest = [h for h in helices if h.id not in rows]
    many = Counter(after(h, strands) for h in rest)
    for h in rest:
        n, k = after(h, strands), 0
        if many[n] + (name(n) in taken) > len(_PRIMES):  # more than αn‴ can mark: lettered αna, αnb …
            while _lettered(name(n), k) in taken:
                k += 1
            out[h.id] = _lettered(name(n), k)
        else:
            while _primed(name(n), k) in taken:
                k += 1
            out[h.id] = _primed(name(n), k)
        taken.add(out[h.id])
    return out


def _lettered(name: str, k: int) -> str:
    """αna, αnb … αnz, then αnz′, αnz″ … past the alphabet (never seen: 26 helices after one strand)."""
    return name + _LETTERS[k] if k < len(_LETTERS) else _primed(name + _LETTERS[-1], k - len(_LETTERS) + 1)


def _primed(name: str, k: int) -> str:
    return name + (_PRIMES[k] if k < len(_PRIMES) else "′" * k)


# ---------------------------------------------------------------------------------------------- YAML files


_SCALAR = (str, int, float)  # a reference or a name: YAML may read `1` or `1.5` as a number


def parse_names(text: str, where: str) -> dict:
    """A names file's text as its mapping ({} when empty), checked: the keys elements and stems only, each a mapping
    of a reference to a name (a single value, never a list or a mapping; read with no YAML anchors or aliases). `where`
    names the file in the messages. The references are checked against the figure later (load_names)."""
    doc = read_yaml(text, "names", where) or {}
    if not isinstance(doc, dict):
        raise ValueError(f"{where}: a names file is a mapping with the keys elements and stems")
    for key in doc:
        if key not in ("elements", "stems"):
            raise ValueError(f"{where}: unknown key {key!r}; a names file has elements and stems")
    for key in ("elements", "stems"):
        part = doc.get(key) or {}
        if not isinstance(part, dict):
            raise ValueError(f"{where}: {key} is a mapping of reference to name, e.g. res:A:3-18: A")
        for ref, name in part.items():
            if isinstance(name, datetime.date):  # YAML reads 2001-01-01 as a date
                raise ValueError(
                    f"{where}: the name of {str(ref)[:80]!r} reads as a date; give a name as text, e.g. quote it: "
                    f"'{name.isoformat()}'"
                )
            if not isinstance(ref, _SCALAR) or not isinstance(name, _SCALAR):
                raise ValueError(f"{where}: {str(ref)[:80]!r} is not a name (give each reference one name)")
    return doc


def _file_names(path: Path, layout: Layout, sses, nucleic, whole, chains) -> tuple[dict[str, str], dict[int, str]]:
    doc = parse_names(path.read_text(encoding="utf-8"), str(path))
    chain_of = {s.id: s.chain for s in sses}
    elements: dict[str, str] = {}
    given: dict[tuple[str, str], tuple[str, str]] = {}  # (chain, name) -> (ref, element) that took it
    for ref, name in (doc.get("elements") or {}).items():
        ref, name = str(ref), str(name)
        try:
            ids = _one_piece(_resolves(layout, sses, ref), sses, ref)
        except AmbiguousName as e:
            raise ValueError(f"{path}: {e}") from None
        if not ids and not _left_out(ref, sses, whole, chains):
            raise ValueError(
                f"{path}: no element {ref!r} in the figure; use a label such as A or α1, chain:label, #k or res:a-b"
            )
        for sid in ids:
            _one_each(path, given, (chain_of[sid], name), (ref, sid), "element")
        if ids:
            elements[ref] = name
    stems: dict[int, str] = {}
    taken: dict[tuple[str, str], tuple[str, str]] = {}
    for key, name in (doc.get("stems") or {}).items():
        key, name = str(key), str(name)
        named = _stem_names(path, key, name, layout, nucleic)
        for k in named:
            d = layout.stems[k]
            for chain in sorted({nucleic.strands[s].chain for s in d.strands}):
                _one_each(path, taken, (chain, name), (key, d.parent or d.id), "stem")
        stems.update(named)
    return elements, stems


_PIECE_GAP = 2  # residues: helices of one chain this close are pieces of one helix (ss: a kink, a copy's break)


def _one_piece(ids: list[str], sses: list[SSE], ref: str) -> list[str]:
    """A res: range that spans the pieces of one helix (consecutive helices of a chain no more than _PIECE_GAP
    residues apart: a kink split, or a break made alike in every copy, ss._copy_breaks) names the piece it overlaps
    most, the first on a tie, not every piece (fix round 1: 7TIM's res:B:178-203 named B:178-196 and B:198-203 once
    B was cut like A). Other references are returned as they are."""
    m = re.fullmatch(r"res:(?:[^:]+:)?(-?\d+)[A-Za-z]?(?:-(-?\d+)[A-Za-z]?)?", ref.strip())
    if m is None or len(ids) < 2:
        return ids
    lo, hi = sorted((int(m[1]), int(m[2] if m[2] is not None else m[1])))
    by_id = {s.id: s for s in sses}
    els = sorted((by_id[k] for k in ids), key=lambda s: s.start)

    def overlap(s: SSE) -> int:
        return max(0, min(hi, s.last.seq) - max(lo, s.first.seq) + 1)

    out, run = [], [els[0]]
    for s in els[1:]:
        last = run[-1]
        if s.kind == last.kind == "H" and s.chain == last.chain and s.start - last.end - 1 <= _PIECE_GAP:
            run.append(s)
            continue
        out.append(max(run, key=lambda x: (overlap(x), -x.start)))
        run = [s]
    out.append(max(run, key=lambda x: (overlap(x), -x.start)))
    return [k for k in ids if k in {s.id for s in out}]


def _one_each(path: Path, given: dict, at: tuple[str, str], by: tuple[str, str], what: str) -> None:
    """One name, one element (or stem) per chain: `at` (chain, name) is taken by `by` (ref, which element), unless
    another one already took it. The same name in two chains is fine (H3 α1 and H4 α1)."""
    ref, thing = by
    first = given.setdefault(at, by)
    if first[1] != thing:
        chain, name = at
        also = f"{first[0]!r} and {ref!r} both name" if first[0] != ref else f"{ref!r} names two"
        raise ValueError(
            f"{path}: {also} {what}s of chain {chain} {name!r} ({first[1]}, {thing}); one name, one {what}"
        )


def _resolves(layout: Layout, sses: list[SSE], ref: str) -> list[str]:
    try:
        return resolve(layout, sses, ref)
    except AmbiguousName:
        raise
    except ValueError:
        return []


def _left_out(ref: str, sses: list[SSE], whole: list[SSE] | None, chains: set[str] | None) -> bool:
    """Does the ref name an element of the structure this figure does not draw (another chain, or a part focus
    left out)?"""
    body = ref[4:] if ref.startswith("res:") else ref
    chain = body.split(":", 1)[0] if body.count(":") == 1 and not ref.startswith("#") else None
    if chain is not None and chain not in {s.chain for s in sses}:
        return chains is None or chain in chains
    return bool(whole) and ref.startswith("res:") and bool(_resolves(provisional(whole), whole, ref))


def _stem_names(path: Path, key: str, name: str, layout: Layout, nucleic: Nucleic | None) -> dict[int, str]:
    """The stems (pieces of one drawn stem) holding the nucleotide `key` (chain:number) -> their names: the name, or
    the name plus the piece's suffix when the stem is drawn in straight pieces (nucleic.piece_name)."""
    from .nucleic import piece_name

    chain = key.rsplit(":", 1)[0] if ":" in key else None
    drawn = {st.chain for st in nucleic.strands} if nucleic is not None else set()
    if chain not in drawn:
        if drawn and chain is None:
            raise ValueError(f"{path}: stem key {key!r} must name a nucleotide as chain:number, e.g. A:1")
        return {}  # no nucleic acid drawn, or not this chain: nothing to name

    def holds(d) -> bool:
        return any(
            f"{lab.chain}:{lab.tag}" == key
            for a, i, b, j in d.pairs
            for lab in (nucleic.strands[a].labels[i], nucleic.strands[b].labels[j])
        )

    hit = next((d for d in layout.stems if holds(d)), None)
    if hit is None:
        raise ValueError(f"{path}: no stem holds nucleotide {key!r}; a stem key is a paired nucleotide (chain:number)")
    pieces = [k for k, d in enumerate(layout.stems) if d is hit or (hit.parent and d.parent == hit.parent)]
    if len(pieces) == 1:
        return {pieces[0]: name}
    return {k: piece_name(name, n) for n, k in enumerate(pieces)}


# ---------------------------------------------------------------------------------------------- applying them


def apply_names(
    specs: list[str],
    path,
    layout: Layout,
    sses: list[SSE],
    *,
    whole: list[SSE] | None = None,
    chains: set[str] | None = None,
    copies: dict[str, tuple[str, int]] | None = None,
    primes: dict[str, str] | None = None,
) -> None:
    """Name the placed elements and stems of `layout` in place, spec after spec (a later one wins; a file's label
    refs read the names given before it). Then copies of one protein (`copies`, model.homolog_chains) share names:
    an element no spec named takes the name of the overlapping element of a named copy, so mirror images and other
    copies read alike. Elements lent by neighbouring chains (focus mode) keep their chain's prime (`primes`)."""
    named: dict[str, str] = {}
    for spec in specs:
        elements, stems = load_names(spec, path, layout, sses, layout.nucleic, whole=whole, chains=chains)
        given = set()
        for ref, name in elements.items():
            for sid in _one_piece(resolve(layout, sses, ref), sses, ref):
                layout.placed[sid].label = named[sid] = name
                given.add(sid)
        for k, name in stems.items():
            layout.stems[k].name = name
        _once_per_chain(spec, layout, sses, given)
    shared = _shared(layout, sses, named, copies or {}, whole or sses)
    for sid, name in shared.items():
        layout.placed[sid].label = name
    for sid in sorted({*named, *shared}):
        layout.placed[sid].label += (primes or {}).get(layout.placed[sid].sse.chain, "")


def _once_per_chain(spec: str, layout: Layout, sses: list[SSE], given: set[str]) -> None:
    """After a names spec, no two drawn elements of a chain carry one label. A file naming an element (`given`, the
    ids this spec named) with a label another element of its chain carries is an error; a generated name
    (deposited, beta) is told apart by primes: in sequence, the first element keeps it, each later one takes the
    first prime its chain leaves free."""
    for chain, mine in group_by_chain(
        sorted((s for s in sses if s.id in layout.placed), key=lambda s: s.start)
    ).items():
        holders: dict[str, list[SSE]] = {}
        for s in mine:
            holders.setdefault(layout.placed[s.id].label, []).append(s)
        used = set(holders)
        for name, held in holders.items():
            if len(held) < 2 or not any(s.id in given for s in held):
                continue
            if spec not in KEYWORDS:
                ours = next(s for s in held if s.id in given)
                other = next(s for s in held if s is not ours)
                raise ValueError(
                    f"{spec}: {name!r} names {ours.id} but {other.id} carries it in chain {chain} too; one name, one "
                    "element (name that one too, or choose another name)"
                )
            for s in held[1:]:
                k = 1
                while _primed(name, k) in used:
                    k += 1
                layout.placed[s.id].label = _primed(name, k)
                used.add(_primed(name, k))


def _shared(
    layout: Layout, sses: list[SSE], named: dict[str, str], copies: dict[str, tuple[str, int]], whole: list[SSE]
):
    """sid -> name for each unnamed drawn element of a copy whose mate (model.copy_classes over `whole`, the chains
    as read) in another copy was named: the reference copy's name, else the first named copy's, unless its chain
    already uses that name."""
    drawn = [s for s in sses if s.id in layout.placed]
    mate = copy_classes(whole, copies).mate
    name: dict[str, str] = {}  # mate -> name
    for s in sorted((s for s in drawn if s.id in named and s.id in mate), key=lambda s: s.id != mate[s.id]):
        name.setdefault(mate[s.id], named[s.id])
    used: dict[str, set[str]] = {}
    for s in drawn:
        used.setdefault(s.chain, set()).add(layout.placed[s.id].label)
    out = {}
    for s in drawn:
        given = name.get(mate.get(s.id, ""))
        if s.id in named or given is None or given in used[s.chain]:
            continue
        out[s.id] = given
        used[s.chain].add(given)
    return out
