# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Facts derived from the structure that the layout and annotation layers need."""

from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .model import SSE, Backbone, Sheet, copy_classes, group_by_chain

CONTACT_CUTOFF = 8.0  # Å between CA atoms
MIN_CONTACT_PAIRS = 3  # CA pairs needed to call two elements neighbours


def sse_contacts(bb: Backbone, sses: list[SSE], cutoff: float = CONTACT_CUTOFF) -> dict[frozenset[str], int]:
    """Number of CA-CA pairs within `cutoff` between each pair of different elements."""
    owner = np.full(len(bb), -1)
    for k, s in enumerate(sses):
        owner[s.start : s.end + 1] = k
    idx = np.flatnonzero(owner >= 0)
    if len(idx) < 2:
        return {}
    pairs = cKDTree(bb.ca[idx]).query_pairs(cutoff, output_type="ndarray")
    counts: Counter[frozenset[str]] = Counter()
    for a, b in pairs:
        ka, kb = owner[idx[a]], owner[idx[b]]
        if ka != kb:
            counts[frozenset((sses[ka].id, sses[kb].id))] += 1
    return {k: v for k, v in counts.items() if v >= MIN_CONTACT_PAIRS}


BUNDLE_MIN = 10  # residues: shorter helices are turns at an interface, not a coiled-coil stalk
BUNDLE_COS = 0.75  # |cos| between axes: parallel or antiparallel
BUNDLE_DIST = 12.0  # Å between axis segments
BUNDLE_COVER = 0.8  # fraction of the shorter helix that must run that close
KINK_GAP = 2  # residues: helix pieces this close in one chain, roughly collinear, count as one helix
KINK_COS = 0.5  # a kink bends the axis by up to 60 degrees
KINK_SHORT = 7  # residues: a piece this short has no reliable axis of its own, so it continues its neighbour


def _segment(s: SSE, n: int = 20) -> np.ndarray:
    half = (len(s) - 1) * 0.75  # 1.5 Å rise per residue
    return np.linspace(s.centroid - s.axis * half, s.centroid + s.axis * half, n)


def _shape(pieces: list[SSE]) -> tuple[np.ndarray, int, np.ndarray]:
    """Points along a helix (one piece, or the kinked pieces of one), its length and its axis (longest piece).
    Kinked pieces get a point per residue, so each weighs by its length."""
    if len(pieces) == 1:
        return _segment(pieces[0]), len(pieces[0]), pieces[0].axis
    pts = np.vstack([_segment(s, max(2, len(s))) for s in pieces])
    return pts, sum(len(s) for s in pieces), max(pieces, key=len).axis


def helix_bundles(sses: list[SSE]) -> list[list[str]]:
    """Helices of different chains packed along their length (coiled-coil stalks, zippers, interface pairs):
    groups of element ids, each spanning at least two chains, sorted. Kinked pieces of one helix stay together."""
    return _helix_packing(sses)[0]


def helix_partners(sses: list[SSE]) -> dict[str, set[str]]:
    """Each helix (element id) -> the helices of other chains packed along it directly (every piece of a kinked
    one): its neighbours in a bundle, not the whole bundle (a ring of twelve helices lends a subunit two, random-100
    C25)."""
    return _helix_packing(sses)[1]


def _helix_packing(sses: list[SSE]) -> tuple[list[list[str]], dict[str, set[str]]]:
    """helix_bundles' groups and helix_partners' direct packing, from one pass."""
    helices = [s for s in sses if s.kind == "H"]
    group = {s.id: s.id for s in helices}

    def root(x: str) -> str:
        while group[x] != x:
            x = group[x]
        return x

    def join(a: str, b: str) -> None:
        group[root(a)] = root(b)

    kinked: dict[str, list[SSE]] = {}
    nexts = list(zip(helices, helices[1:]))
    cos = [float(a.axis @ b.axis) for a, b in nexts]
    kink = [  # pieces k and k+1 are one helix
        a.chain == b.chain and b.start - a.end - 1 <= KINK_GAP and (x > KINK_COS or min(len(a), len(b)) < KINK_SHORT)
        for (a, b), x in zip(nexts, cos)
    ]
    for k in range(1, len(helices) - 1):  # a short piece joins the pieces either side only when they are collinear
        a, b, c = helices[k - 1 : k + 2]
        if len(b) < KINK_SHORT and kink[k - 1] and kink[k] and float(a.axis @ c.axis) <= KINK_COS:
            weaker = k - 1 if cos[k - 1] < cos[k] else k
            if cos[weaker] <= KINK_COS:  # it continues the neighbour it is closer in direction to
                kink[weaker] = False
    for k in np.flatnonzero(kink):
        join(helices[k].id, helices[k + 1].id)
    for s in helices:
        kinked.setdefault(root(s.id), []).append(s)
    shapes = [[s] for s in helices] + [p for p in kinked.values() if len(p) > 1]  # pieces, and whole kinked helices
    long_ = [(p, *_shape(p)) for p in shapes if sum(len(s) for s in p) >= BUNDLE_MIN]
    linked = set()
    partners: dict[str, set[str]] = {}
    for i, (a, pa, la, xa) in enumerate(long_):
        for b, pb, lb, xb in long_[i + 1 :]:
            if a[0].chain == b[0].chain or abs(float(xa @ xb)) < BUNDLE_COS:
                continue
            d = np.linalg.norm(pa[:, None] - pb[None], axis=2)
            near = d.min(1) if la <= lb else d.min(0)
            if (near < BUNDLE_DIST).mean() >= BUNDLE_COVER:
                join(a[0].id, b[0].id)
                linked |= {a[0].id, b[0].id}
                for x, y in ((a, b), (b, a)):
                    for piece in x:
                        partners.setdefault(piece.id, set()).update(q.id for q in y)
    out: dict[str, list[str]] = {}
    for s in helices:
        out.setdefault(root(s.id), []).append(s.id)
    chain = {s.id: s.chain for s in helices}
    groups = sorted(sorted(g) for g in out.values() if len({chain[k] for k in g}) > 1 and linked & set(g))
    return groups, partners


ON_SHEET_MIN = 10  # residues: a shorter helix is a turn at the sheet's edge, not one packed on its face
ON_SHEET_STRANDS = 3  # strands: a smaller sheet has no face to pack on
ON_SHEET_REACH = 10.0  # Å between a helix CA and the nearest sheet CA, counting it as touching the sheet
ON_SHEET_TOUCH = 0.3  # fraction of the helix's CAs that must touch the sheet
ON_SHEET_FLAT = 0.5  # |cos| between the helix axis and the sheet normal: the helix lies along the face
ON_SHEET_HEIGHT = (6.0, 14.0)  # Å from the sheet plane to the helix centroid: one helix packed on the face
ON_SHEET_EDGE = 2.0  # Å: the centroid, projected on the sheet plane, lies within the sheet's extent plus this


def helices_on_sheets(bb: Backbone, sses: list[SSE], sheets: list[Sheet]) -> dict[str, int]:
    """Helices packed on a sheet face: helix id -> index of the sheet in `sheets`. A helix of >= ON_SHEET_MIN residues
    packs on a sheet of >= ON_SHEET_STRANDS strands, one of them of the helix's chain, when >= ON_SHEET_TOUCH of its
    CAs have a sheet CA within ON_SHEET_REACH, its axis lies along the face (|axis . n| < ON_SHEET_FLAT, n the sheet
    normal: the smallest singular vector of the sheet's CAs), its centroid sits ON_SHEET_HEIGHT from the sheet plane
    and, projected on that plane, within the sheet CAs' in-plane extent + ON_SHEET_EDGE. A helix on several sheets
    goes to the one most of its CAs touch (the first on a tie). Coordinates, axes and the plane are snapped before
    any threshold sees them (layout._snap). Plan 11 measured: 1UBQ α1, 2ACY A:22-32 and A:55-67, 2LZM A:39-50."""
    from .layout import _snap as _grid  # layout imports this module: imported when called

    ca = _grid(bb.ca)
    faces = []
    for k, sh in enumerate(sheets):
        if len(sh.strands) < ON_SHEET_STRANDS:
            continue
        pts = np.vstack([ca[s.start : s.end + 1] for s in sh.strands])
        c = _grid(pts.mean(axis=0))
        _, _, vt = np.linalg.svd(pts - c)
        e1, e2, n = _grid(vt[0]), _grid(vt[1]), _grid(vt[2])
        flat = (pts - c) @ np.array([e1, e2]).T
        faces.append((k, {s.chain for s in sh.strands}, cKDTree(pts), c, n, e1, e2, flat.min(0), flat.max(0)))
    out = {}
    for h in sses:
        if h.kind != "H" or len(h) < ON_SHEET_MIN:
            continue
        own, axis, centre = ca[h.start : h.end + 1], _grid(h.axis), _grid(h.centroid)
        best = None
        for k, chains, tree, c, n, e1, e2, lo, hi in faces:
            if h.chain not in chains:
                continue
            touch = sum(bool(hits) for hits in tree.query_ball_point(own, ON_SHEET_REACH))
            at = _grid([(centre - c) @ e1, (centre - c) @ e2])
            if (
                touch >= ON_SHEET_TOUCH * len(h)
                and abs(float(axis @ n)) < ON_SHEET_FLAT
                and ON_SHEET_HEIGHT[0] < abs(float(_grid((centre - c) @ n))) < ON_SHEET_HEIGHT[1]
                and np.all(at >= lo - ON_SHEET_EDGE)
                and np.all(at <= hi + ON_SHEET_EDGE)
                and (best is None or touch > best[1])
            ):
                best = (k, touch)
        if best is not None:
            out[h.id] = best[0]
    return out


DOMAIN_MIN_SSE_RESIDUES = 30  # residues in helices/strands an automatic domain needs
# The residue-level split (residue_domains), calibrated on 22 of the 33 CATH chains of
# tests/data/domains/cath_reference.json (plan 12 Task 3; tests/test_domain_score.py holds the targets):
DOMAIN_CONTACT = 12.0  # Å between CA atoms: two residues in contact, counted for the split
DOMAIN_SEPARATION = 6  # residues apart in sequence at least: nearer pairs are the chain, not the fold
DOMAIN_NCUT = 0.37  # a part splits when its best cut's normalized cut is below this
DOMAIN_MIN_RESIDUES = 40  # residues a domain needs
DOMAIN_MIN_SEGMENT = 20  # residues each piece of a cut needs (a domain may be two pieces: 5JDN's alternate)
DOMAIN_MAX_POINTS = 800  # longer chains are split on blocks of consecutive residues, this many blocks at most


def residue_domains(ca: np.ndarray, keep=()) -> np.ndarray:
    """One chain's domains, residue by residue (CA coordinates in chain order): the domain number of each residue,
    numbered by first residue. The chain is cut in two again and again, as long as a part splits clearly: each cut
    takes one or two runs of the part's residues (in sequence) as one side and the rest as the other, so a domain
    may be discontinuous (5JDN's two domains alternate along its chain, as CATH chops it). The best cut has the lowest
    normalized cut (contacts across / contacts of either side, summed: Shi and Malik 2000) over residue contacts
    (DOMAIN_CONTACT, DOMAIN_SEPARATION); it is made when that is below DOMAIN_NCUT and both sides keep
    DOMAIN_MIN_RESIDUES, each piece DOMAIN_MIN_SEGMENT (_best_cut). No cut parts the residues of a group of
    `keep` (residue indices: a βα row's elements, drawn as one rigid unit)."""
    n = len(ca)
    block = max(1, -(-n // DOMAIN_MAX_POINTS))  # residues per point: a long chain is split on blocks of residues
    m = -(-n // block)
    pairs = cKDTree(ca).query_pairs(DOMAIN_CONTACT, output_type="ndarray") if n > 1 else np.zeros((0, 2), int)
    pairs = pairs[np.abs(pairs[:, 0] - pairs[:, 1]) >= DOMAIN_SEPARATION]
    graph = np.zeros((m, m))
    np.add.at(graph, (pairs[:, 0] // block, pairs[:, 1] // block), 1.0)
    graph += graph.T
    least, piece = -(-DOMAIN_MIN_RESIDUES // block), -(-DOMAIN_MIN_SEGMENT // block)
    groups = [np.unique(np.asarray(g, dtype=int) // block) for g in keep if len(g)]
    parts, done = [np.arange(m)], []
    while parts:
        part = parts.pop()
        whole = [np.isin(part, g) for g in groups]
        side = _best_cut(graph[np.ix_(part, part)], least, piece, [w for w in whole if w.any()])
        if side is None:
            done.append(part)
        else:
            parts += [part[~side], part[side]]
    out = np.zeros(m, dtype=int)
    for d, part in enumerate(sorted(done, key=lambda p: int(p.min()))):
        out[part] = d
    return np.repeat(out, block)[:n]


def _best_cut(graph: np.ndarray, least: int, piece: int, whole=()) -> np.ndarray | None:
    """The side of the part's best cut (residue_domains) as a mask over the part's points, None when the part does
    not split: the run [a, b) of lowest normalized cut, or that run with one more run of the rest when that cuts
    lower (5JDN's domains alternate four times along the chain). Each side keeps `least` points and each piece
    `piece`; no side holds only some of a group of `whole` (masks over the part's points). Every run is weighed at
    once from the cumulative sums of the contact matrix."""
    m = len(graph)
    if m < 2 * least:
        return None
    cum = np.zeros((m + 1, m + 1))
    cum[1:, 1:] = graph.cumsum(0).cumsum(1)
    degree = np.concatenate([[0.0], graph.sum(1).cumsum()])
    total = degree[m]
    i, j = np.meshgrid(np.arange(m + 1), np.arange(m + 1), indexing="ij")
    inner = cum[j, j] - cum[i, j] - cum[j, i] + cum[i, i]  # contacts within the run [i, j) (each counted twice)
    vol = degree[j] - degree[i]  # contacts of the run's points
    size = j - i
    counts = [np.concatenate([[0], np.cumsum(mask)]) for mask in whole]

    def ncut(within, volume, ok):  # alike on any CPU: the sums are exact integer counts, each division correctly
        # rounded; the rounding to 1e-9 only keeps the printed values short, it guards no tie
        with np.errstate(divide="ignore", invalid="ignore"):
            across = volume - within
            value = across / volume + across / (total - volume)
        return np.where(ok & (volume > 0) & (total - volume > 0), np.round(value, 9), np.inf)

    ok = (size >= least) & (m - size >= least) & (size >= piece) & ~((i == 0) & (j == m))
    ok &= ((i == 0) | (i >= piece)) & ((j == m) | (m - j >= piece))
    for count in counts:
        inside = count[j] - count[i]
        ok &= (inside == 0) | (inside == count[-1])
    one = ncut(inner, vol, ok)
    a, b = np.unravel_index(int(np.argmin(one)), one.shape)
    if not np.isfinite(one[a, b]):
        return None
    side = np.zeros(m, dtype=bool)
    side[a:b] = True
    best = one[a, b]
    with_a = np.concatenate([[0.0], graph[:, a:b].sum(1).cumsum()])  # each point's contacts with [a, b), summed
    ok = (size >= piece) & (m - (b - a) - size >= least)  # a second run [i, j) beside [a, b), the rest kept whole
    ok &= ((j <= a) & ((i == 0) | (i >= piece)) & ((j == a) | (a - j >= piece))) | (
        (i >= b) & ((i == b) | (i - b >= piece)) & ((j == m) | (m - j >= piece))
    )
    for count in counts:
        inside = count[b] - count[a] + count[j] - count[i]
        ok &= (inside == 0) | (inside == count[-1])
    two = ncut(inner[a, b] + 2 * (with_a[j] - with_a[i]) + inner, vol[a, b] + vol, ok)
    c, d = np.unravel_index(int(np.argmin(two)), two.shape)
    if two[c, d] < best:
        side[c:d], best = True, two[c, d]
    return side if best < DOMAIN_NCUT else None


def auto_domains(
    sses: list[SSE],
    contacts: dict[frozenset[str], int],
    together: list[list[str]] = (),
    copies: dict[str, tuple[str, int]] | None = None,
    lent: set[str] = frozenset(),
    *,
    bb: Backbone,
) -> list[tuple[str, list[str]]]:
    """Domains, chain by chain: each chain's split (_chain_domains, from its residues in the backbone `bb`), named
    D1, D2... in sequence order (with " · chain" when several chains are split). Copies of one protein (`copies`,
    model.homolog_chains) share their reference copy's split (model.copy_classes): each element joins the domain
    of its mate there, an element the reference lacks the domain of the element before it (else after it) in its
    chain. Chains lending a few elements to the drawn subunit (`lent`, focus mode) add them to the domain they
    complete instead of having their own: the drawn domain they have most contacts with (1LDM's N-terminal strand
    of a neighbour pairs in the catalytic domain's sheet, though its mate is in the Rossmann domain), else their
    mates'. [] when no chain reads as more than one domain."""
    classes = copy_classes(sses, copies or {})
    mine = group_by_chain(sses)
    root = {c: (copies or {}).get(c, (c, 0))[0] for c in mine}
    owner: dict[str, tuple[str, int]] = {}  # element -> (chain, its domain's number there)
    for chain, els in mine.items():
        if classes.reference.get(root[chain], chain) == chain:
            index = [k for k, l in enumerate(bb.labels) if l.chain == chain]
            for d, ids in enumerate(_chain_domains(els, contacts, together, (bb.ca[index], index))):
                owner.update({k: (chain, d) for k in ids})
    drawn: dict[str, tuple[str, int]] = {}  # the drawn chains' elements and their domains, once all are known
    for chain, els in sorted(mine.items(), key=lambda kv: kv[0] in lent):  # the lending chains last
        if classes.reference.get(root[chain], chain) == chain:
            continue
        if chain in lent and not drawn:
            drawn = {k: key for k, key in owner.items() if key[0] not in lent}
        got = [owner.get(classes.mate[s.id]) for s in els]
        if not any(got):
            continue
        for k, s in enumerate(els):
            near = got[k] or next((g for g in reversed(got[:k]) if g), None) or next(g for g in got[k:] if g)
            owner[s.id] = (_host(s.id, drawn, contacts) or near) if chain in lent else (chain, near[1])
    ref = {c: classes.reference.get(root[c], c) for c in mine}
    _fold_small_copies(owner, ref, sses, contacts)
    groups: dict[tuple[str, int], list[str]] = {}
    for s in sses:
        if s.id in owner:
            groups.setdefault(owner[s.id], []).append(s.id)
    split = list(dict.fromkeys(chain for chain, _ in groups))
    out = []
    for chain in split:
        for d in sorted(d for c, d in groups if c == chain):  # a copy's domains are numbered as its reference's
            out.append((f"D{d + 1}" if len(split) == 1 else f"D{d + 1} · {chain}", groups[(chain, d)]))
    return out


def _host(k: str, drawn: dict[str, tuple[str, int]], contacts: dict[frozenset[str], int]) -> tuple | None:
    """The domain of the `drawn` elements (element -> its domain) element `k` has most contacts with, None when it
    touches none; a tie goes to the earlier domain."""
    near: Counter = Counter()
    for m, key in drawn.items():
        near[key] += contacts.get(frozenset((k, m)), 0)
    best = [key for key in sorted(near) if near[key] > 0]
    return max(best, key=lambda key: near[key]) if best else None


def _fold_small_copies(
    owner: dict[str, tuple[str, int]], ref: dict[str, str], sses: list[SSE], contacts: dict[frozenset[str], int]
) -> None:
    """In place: a domain that some copy of the protein holds with fewer than DOMAIN_MIN_SSE_RESIDUES residues of
    elements (a copy missing elements its reference has: 8UUP chain F's D5) folds, in every copy alike, into the
    domain of the reference it is best connected to (contacts, plus one per neighbour in sequence); the domains are
    then numbered again from 0 in each reference, and a protein left with one domain is not split."""
    size = {s.id: len(s) for s in sses}
    order = [s.id for s in sses]

    def root_of(key: tuple[str, int]) -> tuple[str, int]:  # a domain's place in its reference
        return ref.get(key[0], key[0]), key[1]

    alone: set[tuple[str, int]] = set()  # small domains with nothing to fold into: their protein is not split
    while True:
        total: Counter = Counter()
        for k, key in owner.items():
            total[key] += size[k]
        small = sorted({root_of(key) for key, n in total.items() if n < DOMAIN_MIN_SSE_RESIDUES} - alone)
        if not small:
            break
        r, d = small[0]
        mine = {k for k, key in owner.items() if key == (r, d)}
        link: Counter = Counter()
        for k, key in owner.items():
            if key[0] == r and key[1] != d:
                link[key[1]] += sum(contacts.get(frozenset((k, m)), 0) for m in sorted(mine))
        for a, b in zip(order, order[1:]):
            for x, y in ((a, b), (b, a)):
                if x in mine and owner.get(y, (None,))[0] == r and owner[y][1] != d:
                    link[owner[y][1]] += 1
        others = sorted({key[1] for key in owner.values() if key[0] == r and key[1] != d})
        if not others:
            alone.add((r, d))
            continue
        into = max(others, key=lambda e: (link[e], -e))
        for k, key in list(owner.items()):
            if root_of(key) == (r, d):
                owner[k] = (key[0], into)
    for r in sorted(set(ref.values())):
        numbers = sorted({key[1] for key in owner.values() if root_of(key)[0] == r})
        for k, key in list(owner.items()):
            if root_of(key)[0] == r:
                if len(numbers) < 2:
                    del owner[k]
                else:
                    owner[k] = (key[0], numbers.index(key[1]))


def domain_tones(
    domains: list[tuple[str, list[str]]], sses: list[SSE], copies: dict[str, tuple[str, int]]
) -> list[int]:
    """The tone (index into the domain tones) of each domain: copies of one domain in copies of a protein share one,
    distinct domains take the next, in domain order. A domain is a copy of the one holding most of its elements'
    mates (model.copy_classes); a domain whose mates lie in no domain, or in itself, is its own."""
    mate = copy_classes(sses, copies).mate
    owner = {k: d for d, (_, ids) in enumerate(domains) for k in ids}
    rep = []
    for d, (_, ids) in enumerate(domains):
        votes = Counter(owner[m] for m in (mate.get(k, k) for k in ids) if m in owner)
        rep.append(max(votes, key=lambda e: (votes[e], e == d, -e)) if votes else d)

    def root(d: int) -> int:
        seen = set()
        while rep[d] != d and d not in seen:
            seen.add(d)
            d = rep[d]
        return d

    tone: dict[int, int] = {}
    return [tone.setdefault(root(d), len(tone)) for d in range(len(domains))]


DOMAIN_STRAY = 20  # residues of elements: a run of one domain's elements this small, between two runs of another


def _no_strays(domain: list[int], size: list[int]) -> list[int]:
    """Each element's domain (in chain order) with the strays taken home: a run of one domain's elements holding
    fewer than DOMAIN_STRAY element residues, with one other domain's elements on both sides, joins that domain (a
    helix a cut left in the domain across: 1LDM's 247-261 in the Rossmann domain, its panel then reaching over the
    catalytic one). A run at either end of the chain stays (a terminal helix packed on the other domain). One run
    moves at a time, the first from the N-terminus, and the runs are read again after each move."""
    out = list(domain)
    while True:
        runs = []
        for k, d in enumerate(out):
            if runs and runs[-1][0] == d:
                runs[-1][2] = k + 1
            else:
                runs.append([d, k, k + 1])
        stray = next(
            (
                (a, b, before)
                for (before, _, _), (d, a, b), (after, _, _) in zip(runs, runs[1:], runs[2:])
                if before == after != d and sum(size[a:b]) < DOMAIN_STRAY
            ),
            None,
        )
        if stray is None:
            return out
        a, b, home = stray
        out[a:b] = [home] * (b - a)


def _chain_domains(
    sses: list[SSE], contacts: dict[frozenset[str], int], together, residues: tuple[np.ndarray, list[int]]
) -> list[list[str]]:
    """One chain's domains, in sequence order: its residue domains (residue_domains; `residues`: the CA coordinates
    of the chain's residues and their indices in the backbone), each element in the one holding most of its
    residues, a stray element taken home (_no_strays). Each holds at least DOMAIN_MIN_SSE_RESIDUES residues of
    elements (smaller ones fold into their best-connected neighbour: element contacts, plus one per neighbour in
    sequence). A βα row of `together` (layout.template_rows, drawn as one rigid unit) is never split: no cut
    parts its residues, and should its elements still lie in several domains, the domains holding its strands
    merge and its helices follow its strands, so the row lies in one domain's panel (1LDM's Rossmann sheet); what
    that leaves too small folds again (1TIM: a 10-residue domain of two helices). [] when the chain reads as a
    single domain."""
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(s.id for s in sses)
    for key, n in contacts.items():
        a, b = sorted(key)
        if a in graph and b in graph:
            graph.add_edge(a, b, weight=float(n))
    for a, b in zip(sses, sses[1:]):
        w = graph.get_edge_data(a.id, b.id, {}).get("weight", 0.0)
        graph.add_edge(a.id, b.id, weight=w + 1.0)
    size = {s.id: len(s) for s in sses}

    def fold(groups: list[set[str]]) -> list[set[str]]:  # a domain needs enough structure: the smallest folds
        groups = list(groups)  # into its best-connected neighbour
        while len(groups) > 1:
            small = min(groups, key=lambda g: sum(size[k] for k in g))
            if sum(size[k] for k in small) >= DOMAIN_MIN_SSE_RESIDUES:
                break
            groups.remove(small)
            link = [sum(graph.get_edge_data(a, b, {}).get("weight", 0.0) for a in small for b in g) for g in groups]
            groups[int(np.argmax(link))] |= small
        return groups

    ca, index = residues
    at = {k: n for n, k in enumerate(index)}
    span = {s.id: range(s.start, s.end + 1) for s in sses}
    keep = [[at[k] for i in row if i in span for k in span[i] if k in at] for row in together]
    label = dict(zip(index, residue_domains(ca, keep).tolist()))
    mine = []
    for s in sses:
        votes = Counter(label[k] for k in span[s.id] if k in label)
        mine.append(min(votes, key=lambda d: (-votes[d], d)) if votes else 0)
    found: dict[int, set[str]] = {}
    for s, d in zip(sses, _no_strays(mine, [len(s) for s in sses])):
        found.setdefault(d, set()).add(s.id)
    groups = fold(found.values())
    if len(groups) < 2:
        return []
    kind = {s.id: s.kind for s in sses}
    for keep in together:  # a βα row stays whole: its strands' communities merge, its helices follow its strands
        strands = {k for k in keep if kind.get(k) == "E"}
        if not strands:
            continue
        mine = [g for g in groups if g & strands]
        rest = [g - set(keep) for g in groups if not g & strands]
        groups = [g for g in [*rest, set().union(*mine, set(keep) & set(kind))] if g]
    groups = fold(groups)
    if len(groups) < 2:
        return []
    start = {s.id: s.start for s in sses}
    return sorted((sorted(g) for g in groups), key=lambda g: min(start[k] for k in g))
