# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""What changed between two structures of one protein chain (plan 13e, `foldmap diff`): A the reference, B the other
state (apo and holo, open and closed) or a homologue.

The residues are mapped as 13a maps them (align.align_chains): by sequence when the two chains are one protein by
A19's rule (criteria.one_protein_names: the rule criteria.pair_check uses on the random-100 pairs), else by
structure; `match` may also name the route or a two-row alignment file. The elements are the figure's
(ss.build_sses on the chain alone, as cli.make_layout builds them for one drawn chain), paired by
align.element_mates. From those:

- only_a, only_b: the elements of one structure with no mate in the other (Mates.only_a, Mates.extras);
- ss_changes: each residue whose element kind differs (H, G 3-10, E strand, C none; "-" a residue the other
  structure lacks): an aligned pair whose kinds differ, a residue of A in an element with no partner in B, a
  residue of B in an element with no partner in A;
- pairing: each two strands paired (sheets._edge_table: at least sheets._MIN_PAIRS DSSP bridges) in one structure
  and not in the other, and each pairing in both whose register moved (datafile._pair's bridges, the shift of the
  median_low of j - s·i as datafile._shear reads a pairing: in B the residue across from each of the first strand's
  residues lies k further along the second, counted in A's residues). Strands are named by A's ids (B's through
  their mates); a pairing whose kind turned (parallel / antiparallel) is lost and gained;
- displacement: each aligned B residue's Cα distance to its A residue after the core superposition (_core_fit);
- domain_turns: with domains (name -> residue ranges of A's chain, res:A:1-60, A:1-60 or 1-60), each domain's rigid
  rotation, in degrees, left after the core superposition (the Kabsch fit of B's residues of it onto A's);
- contacts: element and interface changes measured as foldmap data's tables measure them (_contacts).

The core superposition (_core_fit), a fixed rule: the core is the aligned pairs inside mated elements (mates and
the pieces of a split helix; every aligned pair when fewer than 3 are). It starts from the superposition the
alignment's TM-score found (align.tm_score: TM-align's search, which settles on the largest part that superposes,
where one fit on every pair would split the difference between two domains that moved): the core pairs within
CORE_CUT·d0(L) Å of each other under it (L: A's chain length, align.d0; the whole core when fewer than 3 are) are
fitted by Kabsch, the core pairs farther apart than that under the new fit dropped and the rest refitted, until the
kept pairs repeat, at most align.ROUNDS fits, never fewer than 3 pairs.

Determinism (A12): distances are snapped to 1/align.GRID Å (align._dist) before the core's cut, turns are rounded to
0.001 degree, every list is in a fixed order (A's residues and elements N to C, then B's; partners and ligand
names sorted), and every set compared is of names, never of floats.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from statistics import median_low
from types import SimpleNamespace

import numpy as np

from .align import ROUNDS, ChainAlignment, Mates, _dist, align_chains, d0, element_mates
from .model import SSE, Backbone

CORE_CUT = 2.0  # core pairs farther apart than CORE_CUT·d0 after a fit are dropped from the next (_core_fit)
NONE, ABSENT = "C", "-"  # a residue in no element; a residue the other structure lacks
CONTACT_KEYS = ("elements", "interfaces", "interface_elements", "ligands", "ligand_elements")


@dataclass(frozen=True)
class Diff:
    """What changed from A (the reference) to B. Residues are written (chain, number with insertion code)."""

    alignment: ChainAlignment
    mates: Mates
    only_a: tuple[str, ...]  # A's element ids with no mate in B, N to C
    only_b: tuple[str, ...]  # B's element ids with no mate in A, N to C
    # (chain, residue, kind in A, kind in B): A's chain and number, B's for a residue only B has (kind in A "-")
    ss_changes: tuple[tuple[str, str, str, str], ...]
    pairing: tuple[tuple[str, str, str], ...]  # (A strand id, A strand id, "gained" | "lost" | "register ±k")
    displacement: dict[tuple[str, str], float]  # B residue -> Cα distance to its A residue after the core fit, Å
    domain_turns: dict[str, float]  # domain -> its rigid rotation after the core fit, degrees (domains given only)
    contacts: dict[str, tuple]  # CONTACT_KEYS -> changes (_contacts)
    core: tuple[tuple[int, int], ...] = ()  # the aligned pairs the core fit kept (indices as alignment.pairs)
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3), compare=False)  # A ≈ R·B + t, the core fit
    translation: np.ndarray = field(default_factory=lambda: np.zeros(3), compare=False)


@dataclass(frozen=True)
class _Side:
    """One structure as diff reads it: the chain alone (its backbone, elements and DSSP bridges, as a figure of that
    chain builds them), and the whole structure (every protein chain, the model, the ligands) for the contacts."""

    name: str
    chain: str
    bb: Backbone
    sses: list
    bridges: list
    full: Backbone
    model: object
    links: object


def _side(source, chain: str | None, side: str, assembly: str, helices_310: bool) -> _Side:
    from .cli import resolve_structure
    from .dssp import assign_dssp
    from .io import load_backbone, load_links, read_model
    from .ss import build_sses

    path = resolve_structure(str(source))
    model = read_model(path, assembly)
    full = load_backbone(path, assembly, model=model)
    present = list(dict.fromkeys(l.chain for l in full.labels))
    if chain is None:
        if not present:
            raise ValueError(f"{side} ({source}) has no protein chain")
        chain = present[0]
    elif chain not in present:
        raise ValueError(f"no chain {chain!r} in {side} ({source}); its protein chains: {', '.join(present)}")
    bb = full.only([chain])
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, short_helices=helices_310, bridges=dssp.bridges)
    links = load_links(path, full, assembly, model=model)
    return _Side(str(source), chain, bb, sses, list(dssp.bridges), full, model, links)


def diff_structures(
    a,
    chain_a: str | None,
    b,
    chain_b: str | None,
    match: str = "auto",
    domains: dict | None = None,
    *,
    assembly: str = "auto",
    look=None,
) -> Diff:
    """What changed from chain_a of `a` to chain_b of `b` (files or PDB ids; a chain None: the structure's first
    protein chain). match: auto (sequence when the chains are one protein by A19's rule, else structure),
    sequence, structure, or a two-row alignment file (A's chain first). domains: {name: residue range or ranges of
    A's chain} for domain_turns. Both read as `assembly` reads them; `look` (a Style) decides whether short 3-10
    helices are elements (Style.helices_310), as the figure does."""
    from .style import Style

    helices_310 = (look or Style()).helices_310
    sa = _side(a, chain_a, "A", assembly, helices_310)
    sb = _side(b, chain_b, "B", assembly, helices_310)
    return diff_sides(sa, sb, match, domains)


def _route(sa: _Side, sb: _Side, match: str) -> str:
    """auto: the sequence route when the two chains are one protein (A19's rule), else the structure route."""
    if match != "auto":
        return match
    from .criteria import one_protein_names

    one = one_protein_names([l.name for l in sa.bb.labels], [l.name for l in sb.bb.labels])
    return "sequence" if one["one"] else "structure"


def diff_sides(sa: _Side, sb: _Side, match: str = "auto", domains: dict | None = None) -> Diff:
    """diff_structures on two structures already read (_side)."""
    aln = align_chains(sa.bb, sa.chain, sb.bb, sb.chain, _route(sa, sb, match))
    mates = element_mates(sa.sses, sb.sses, aln)
    to_a = {**mates.mate, **mates.pieces}  # B element id -> the A element it is (whole or a piece of)
    core, r, t, dist = _core_fit(sa, sb, aln, to_a)
    displacement = {(sb.chain, aln.residues_b[j].tag): float(d) for (_, j), d in zip(aln.pairs, dist)}
    return Diff(
        alignment=aln,
        mates=mates,
        only_a=mates.only_a,
        only_b=mates.extras,
        ss_changes=_ss_changes(sa, sb, aln),
        pairing=_pairing(sa, sb, aln, mates),
        displacement=displacement,
        domain_turns=_domain_turns(sa, sb, aln, r, t, domains) if domains else {},
        contacts=_contacts(sa, sb, to_a),
        core=core,
        rotation=r,
        translation=t,
    )


def _owner(sses: list[SSE], n: int) -> list[SSE | None]:
    """The element each of the chain's n residues is in (None: a loop)."""
    out: list[SSE | None] = [None] * n
    for s in sses:
        out[s.start : s.end + 1] = [s] * (s.end + 1 - s.start)
    return out


def _core_fit(sa: _Side, sb: _Side, aln: ChainAlignment, to_a: dict[str, str]):
    """(core pairs kept, R, t, every aligned pair's snapped distance under them): the fixed rule of the module
    docstring."""
    from .symmetry import _kabsch

    if not aln.pairs:
        return (), np.eye(3), np.zeros(3), np.zeros(0)
    ia, ib = (np.array(k, dtype=int) for k in zip(*aln.pairs))
    pa, pb = np.asarray(sa.bb.ca, float)[ia], np.asarray(sb.bb.ca, float)[ib]
    own_a, own_b = _owner(sa.sses, len(sa.bb)), _owner(sb.sses, len(sb.bb))
    core = np.array(
        [own_a[i] is not None and own_b[j] is not None and to_a.get(own_b[j].id) == own_a[i].id for i, j in aln.pairs]
    )
    if core.sum() < 3:
        core = np.ones(len(aln.pairs), bool)
    if core.sum() < 3:  # one or two pairs: no rotation to fit
        r, t = np.eye(3), (pa - pb).mean(axis=0)
        return tuple(aln.pairs), r, t, _dist(pa, pb, r, t)
    cut = CORE_CUT * d0(len(sa.bb))
    keep = core & (_dist(pa, pb, aln.rotation, aln.translation) <= cut)  # from the TM-score's superposition
    if keep.sum() < 3:
        keep = core.copy()
    for _ in range(ROUNDS):
        r, t, _ = _kabsch(pb[keep], pa[keep])
        near = core & (_dist(pa, pb, r, t) <= cut)
        if near.sum() < 3 or np.array_equal(near, keep):
            break
        keep = near
    r, t, _ = _kabsch(pb[keep], pa[keep])
    kept = tuple(p for p, k in zip(aln.pairs, keep) if k)
    return kept, r, t, _dist(pa, pb, r, t)


def _kinds(side: _Side) -> list[str]:
    return [s.kind if s is not None else NONE for s in _owner(side.sses, len(side.bb))]


def _ss_changes(sa: _Side, sb: _Side, aln: ChainAlignment) -> tuple[tuple[str, str, str, str], ...]:
    """Along the alignment: the pairs whose kinds differ, then before each pair the residues of A it skips that are
    in an element (kind in B "-"), then those of B (kind in A "-"); the ends after the last pair likewise."""
    ka, kb = _kinds(sa), _kinds(sb)
    la, lb = sa.bb.labels, sb.bb.labels
    out: list[tuple[str, str, str, str]] = []
    i0 = j0 = 0
    for i, j in [*aln.pairs, (len(la), len(lb))]:
        out += [(sa.chain, la[k].tag, ka[k], ABSENT) for k in range(i0, i) if ka[k] != NONE]
        out += [(sb.chain, lb[k].tag, ABSENT, kb[k]) for k in range(j0, j) if kb[k] != NONE]
        if i < len(la) and ka[i] != kb[j]:
            out.append((sa.chain, la[i].tag, ka[i], kb[j]))
        i0, j0 = i + 1, j + 1
    return tuple(out)


def _pairings(side: _Side) -> dict[tuple[str, str], tuple[str, list]]:
    """Each two strands of the chain paired by at least sheets._MIN_PAIRS bridges (sheets._edge_table): (id, id)
    N to C -> (kind P or A, the bridges as (residue of the first, residue of the second), datafile._pair's)."""
    from .datafile import _pair
    from .sheets import _edge_table

    strands = [s for s in side.sses if s.kind == "E"]
    out = {}
    for x, y in sorted(_edge_table(strands, side.bridges)):
        p = _pair(strands[x], strands[y], side.bridges, side.bb.labels)
        out[(strands[x].id, strands[y].id)] = (p["kind"], [(i, j) for i, j, _ in p["_reg"]])
    return out


def _register(kind: str, reg: list[tuple[int, int]]) -> int | None:
    """A pairing's register c (datafile._shear: residue i of the first strand faces c + i parallel, c - i
    antiparallel), the median_low over its bridges; None without one."""
    s = 1 if kind == "P" else -1
    return median_low([j - s * i for i, j in reg]) if reg else None


def _pairing(sa: _Side, sb: _Side, aln: ChainAlignment, mates: Mates) -> tuple[tuple[str, str, str], ...]:
    rank = {s.id: k for k, s in enumerate(sa.sses)}
    to_a_res = {j: i for i, j in aln.pairs}
    in_a = _pairings(sa)
    in_b: dict[tuple[str, str], tuple[str, list]] = {}
    for (x, y), (kind, reg) in _pairings(sb).items():
        if x not in mates.mate or y not in mates.mate:
            continue  # an extra strand's pairings: the strand itself is in only_b
        mx, my = mates.mate[x], mates.mate[y]
        moved = [(to_a_res[i], to_a_res[j]) for i, j in reg if i in to_a_res and j in to_a_res]
        if rank[mx] > rank[my]:
            mx, my, moved = my, mx, [(j, i) for i, j in moved]
        in_b[(mx, my)] = (kind, sorted(moved))
    matched = set(mates.mate.values())
    out = []
    for key in sorted({*in_a, *in_b}, key=lambda k: (rank[k[0]], rank[k[1]])):
        if key[0] not in matched or key[1] not in matched:
            continue  # a strand only A has: in only_a
        a, b = in_a.get(key), in_b.get(key)
        if a is None:
            out.append((*key, "gained"))
        elif b is None:
            out.append((*key, "lost"))
        elif a[0] != b[0]:  # parallel <-> antiparallel: A's pairing lost, another gained
            out += [(*key, "lost"), (*key, "gained")]
        else:
            ca, cb = _register(a[0], a[1]), _register(b[0], b[1])
            if ca is not None and cb is not None and ca != cb:
                out.append((*key, f"register {cb - ca:+d}"))
    return tuple(out)


def _domain_residues(name: str, refs, side: _Side) -> set[int]:
    """The indices of A's chain residues a domain's ranges name (res:A:1-60, A:1-60, 1-60; insertion codes as in
    the file); ValueError naming the domain for a range that is malformed, on another chain, or names nothing."""
    from .layout import RES_RANGE

    got: set[int] = set()
    for ref in [refs] if isinstance(refs, str) else list(refs):
        body = ref.strip().removeprefix("res:")
        chain = None
        if body.count(":") == 1:
            chain, body = body.split(":")
        m = RES_RANGE.fullmatch(body)
        if m is None or chain not in (None, side.chain):
            raise ValueError(
                f"domain {name}: {ref!r} is not a residue range of A's chain {side.chain} "
                f"(res:{side.chain}:1-60, {side.chain}:1-60 or 1-60)"
            )
        lo = (int(m[1]), m[2])
        hi = (int(m[3]), m[4]) if m[3] is not None else lo
        lo, hi = min(lo, hi), max(lo, hi)
        mine = {k for k, l in enumerate(side.bb.labels) if lo <= (l.seq, l.icode) <= hi}
        if not mine:
            raise ValueError(f"domain {name}: {ref!r} names no residue of chain {side.chain} in A ({side.name})")
        got |= mine
    return got


def _domain_turns(sa: _Side, sb: _Side, aln: ChainAlignment, r, t, domains: dict) -> dict[str, float]:
    """Each domain's rotation (degrees, 0.001) of B's aligned residues of it, after the core fit, onto A's."""
    from .symmetry import _kabsch, _rotation

    ia, ib = (np.array(k, dtype=int) for k in zip(*aln.pairs)) if aln.pairs else (np.zeros(0, int),) * 2
    pa = np.asarray(sa.bb.ca, float)[ia]
    pb = np.asarray(sb.bb.ca, float)[ib] @ np.asarray(r).T + t
    out = {}
    for name, refs in domains.items():
        mine = _domain_residues(name, refs, sa)
        sel = np.array([i in mine for i in ia], bool)
        if sel.sum() < 3:
            raise ValueError(
                f"domain {name}: {int(sel.sum())} of its residues aligned with B; its turn needs at least 3"
            )
        rd, _, _ = _kabsch(pb[sel], pa[sel])
        out[name] = round(float(np.degrees(_rotation(rd)[0])), 3) + 0.0
    return out


def _changes(a: set, b: set) -> list[tuple]:
    """What b has that a lacks ("gained") and a has that b lacks ("lost"), each a tuple of the item's parts."""
    return [(*k, "gained") for k in b - a] + [(*k, "lost") for k in a - b]


def _contacts(sa: _Side, sb: _Side, to_a: dict[str, str]) -> dict[str, tuple]:
    """Changes of what the chain touches, as foldmap data measures it (datafile; elements named by A's ids, B's
    through their mates or as pieces; an element only one structure has is in only_a / only_b, not here):

    - elements: two matched elements in contact (features.sse_contacts: at least MIN_CONTACT_PAIRS Cα pairs within
      CONTACT_CUTOFF), (id, id, gained | lost);
    - interfaces: the chains the diffed chain has an interface with (datafile's interfaces table: protein chains by
      Cα pairs, nucleic-acid chains by atom contacts), named by chain id, (partner, measure, gained | lost);
    - interface_elements: for each protein partner of both, the chain's elements touching it (the same
      sse_contacts rule, the partner whole), (id, partner, gained | lost);
    - ligands: the ligands touching the chain (io.load_links, the figure's ligand rule; additives left out, as the
      figure hides them), by name, (name, gained | lost);
    - ligand_elements: for each ligand name of both, the chain's elements touching it, (id, name, gained | lost).
    Partners and ligands are named, not mapped: two files that name one partner differently show it lost and
    gained."""
    rank_a = {s.id: k for k, s in enumerate(sa.sses)}
    matched = set(to_a.values())

    def name_a(side: _Side, i: str) -> str | None:
        return (i if i in matched else None) if side is sa else to_a.get(i)

    def element_pairs(side: _Side) -> set:
        from .features import sse_contacts

        out = set()
        for key in sse_contacts(side.bb, side.sses):
            x, y = (name_a(side, i) for i in key)
            if x is not None and y is not None and x != y:
                out.add(tuple(sorted((x, y), key=rank_a.__getitem__)))
        return out

    def faces(side: _Side) -> dict[str, str]:
        from .datafile import _interfaces

        layout = SimpleNamespace(nucleic=None, chains_read=())
        model = side.model[0] if isinstance(side.model, tuple) else side.model
        out = {}
        for row in _interfaces(layout, [], side.full, None, model):
            if side.chain in row["chains"] and row["contacts"]:
                partner = next(c for c in row["chains"] if c != side.chain)
                out[partner] = row["measure"]
        return out

    def chain_index(side: _Side) -> dict[int, int]:
        """The whole structure's residue index -> the chain's."""
        mine = [k for k, l in enumerate(side.full.labels) if l.chain == side.chain]
        return {k: n for n, k in enumerate(mine)}

    def touching(side: _Side, partners: set[str]) -> set:
        """(element as named in A, partner) for the chain's elements touching each partner protein chain."""
        from .features import sse_contacts

        index = chain_index(side)
        back = {n: k for k, n in index.items()}
        spans = [SimpleNamespace(id=s.id, start=back[s.start], end=back[s.end]) for s in side.sses]
        labels = side.full.labels
        k = 0
        while k < len(labels):  # each partner chain's runs of residues, as one element named "\0chain"
            e = k
            while e + 1 < len(labels) and labels[e + 1].chain == labels[k].chain:
                e += 1
            if labels[k].chain in partners:
                spans.append(SimpleNamespace(id="\0" + labels[k].chain, start=k, end=e))
            k = e + 1
        out = set()
        for key in sse_contacts(side.full, spans):
            ids = sorted(key)
            if len(ids) == 2 and ids[0].startswith("\0") != ids[1].startswith("\0"):
                partner, own = sorted(ids)  # "\0..." sorts first
                x = name_a(side, own)
                if x is not None:
                    out.add((x, partner[1:]))
        return out

    def ligands(side: _Side) -> dict[str, set]:
        """Ligand name -> the chain's elements (as named in A, None: a loop) touching it."""
        index = chain_index(side)
        own = _owner(side.sses, len(side.bb))
        out: dict[str, set] = {}
        for lig in side.links.ligands if side.links is not None else []:
            mine = [index[k] for k in lig.contacts if k in index]
            if lig.additive or not mine:
                continue
            got = out.setdefault(lig.name, set())
            got |= {name_a(side, own[k].id) for k in mine if own[k] is not None}
        return {name: {x for x in got if x is not None} for name, got in out.items()}

    fa, fb = faces(sa), faces(sb)
    both = sorted(c for c in set(fa) & set(fb) if fa[c] == fb[c] == "ca")
    la, lb = ligands(sa), ligands(sb)
    common = sorted(set(la) & set(lb))

    def by_rank(rows):
        return tuple(sorted(rows, key=lambda r: (rank_a[r[0]], *r[1:])))

    return {
        "elements": tuple(
            sorted(_changes(element_pairs(sa), element_pairs(sb)), key=lambda r: (rank_a[r[0]], rank_a[r[1]], r[2]))
        ),
        "interfaces": tuple(sorted(_changes(set(fa.items()), set(fb.items())))),
        "interface_elements": by_rank(_changes(touching(sa, set(both)), touching(sb, set(both)))),
        "ligands": tuple(sorted(_changes({(n,) for n in la}, {(n,) for n in lb}))),
        "ligand_elements": by_rank(
            _changes({(x, n) for n in common for x in la[n]}, {(x, n) for n in common for x in lb[n]})
        ),
    }


def counts(diff: Diff) -> Counter:
    """How many changes of each kind (only_a, only_b, ss_changes, pairing, and each contacts key): zero for every
    kind when B is A."""
    out = Counter({k: len(getattr(diff, k)) for k in ("only_a", "only_b", "ss_changes", "pairing")})
    out.update({f"contacts.{k}": len(v) for k, v in diff.contacts.items()})
    return out
