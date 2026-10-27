"""Independent kernels (plan 28 Task 7, review 3.8): small implementations written here from the published
definitions, not from foldmap's code, compared with what foldmap computes on committed fixtures. A disagreement is a
finding, recorded in docs/validation/metamorphic/report.md (scripts/metamorphic_record.py); each is pinned with
xfail(strict) naming it, the disagreements themselves are pinned exactly (FINDINGS), and the agreement measured is
asserted not to fall.

- H-bond energy: Kabsch & Sander (1983) Biopolymers 22:2577, the electrostatic model of the backbone H-bond,
  E = q1 q2 (1/r(ON) + 1/r(CH) - 1/r(OH) - 1/r(CN)) f, q1 = 0.42 e, q2 = 0.20 e, f = 332 (E in kcal/mol, r in Å),
  the amide H placed 1.0 Å from N opposite the preceding residue's C=O; a bond where E < -0.5 kcal/mol.
- Bridges: the same paper, parallel bridge(i, j) = [Hbond(i-1, j) and Hbond(j, i+1)] or [Hbond(j-1, i) and
  Hbond(i, j+1)]; antiparallel = [Hbond(i, j) and Hbond(j, i)] or [Hbond(i-1, j+1) and Hbond(j-1, i+1)]
  (Hbond(a, b): CO of a to NH of b), i-1..i+1 and j-1..j+1 not overlapping (j >= i + 3), each residue's neighbours
  peptide-bonded to it (C-N <= 2.5 Å, the paper's chain-break rule). Ladder: consecutive bridges of one type; sheet:
  ladders joined by a shared residue.
- Canonical base pairs: Watson-Crick pairs and the G-U wobble (Saenger 1984, Principles of Nucleic Acid Structure,
  ch. 6; Leontis & Westhof 2001, RNA 7:499, the cis Watson-Crick/Watson-Crick family, both glycosidic bonds on one
  side of the pair), with the geometry of Olson et al. 2001 J Mol Biol 313:229 (C1'-C1' about 10.5 Å, each
  glycosidic bond at about 55 degrees to the C1'-C1' line): the Watson-Crick-edge donor-acceptor pairs (G-C: N1-N3,
  N2-O2, O6-N4; A-U/T: N1-N3, N6-O4; G-U/T: N1-O2, O6-N3), at least two within 3.5 Å, C1'-C1' within 9.5-11.5 Å,
  both glycosidic angles within 40-70 degrees (30-80 for the wobble, whose two angles are unequal: Masquida &
  Westhof 2000, RNA 6:9), cis; partners on one chain at least HAIRPIN_MIN + 1 apart (a hairpin loop holds three
  nucleotides or more: Tinoco et al. 1973 Nature New Biol 246:40). Standard nucleotides only: a modified one is
  outside the kernel's scope (its disagreements are named "modified", not findings).

The bonds, bridges and sheets read the coordinates foldmap reads (io.load_backbone's N, CA, C, O per residue), so what
is compared is the definitions, not the file reading; the base-pair kernel reads the file with gemmi itself.
"""

from __future__ import annotations

import functools
from pathlib import Path
from typing import NamedTuple

import gemmi
import numpy as np
import pytest

from foldmap.dssp import _hbond_set, assign_dssp
from foldmap.io import load_backbone, load_nucleic
from foldmap.sheets import build_sheets
from foldmap.ss import build_sses

DATA = Path(__file__).parent / "data"
NA = DATA / "na"
PROTEINS = ["1UBQ", "1TIM", "1EMA", "2LZM", "1C3W"]

# ------------------------------------------------------------------------------------------- Kabsch-Sander 1983

KS_F = 0.42 * 0.20 * 332.0  # q1 q2 f, kcal/mol Å
KS_CUT = -0.5  # kcal/mol
PEPTIDE = 2.5  # Å, C(i-1)-N(i): longer is a chain break


def peptide_bonded(bb) -> np.ndarray:
    """linked[i]: residue i is peptide-bonded to residue i-1 (same chain, C(i-1)-N(i) <= PEPTIDE)."""
    n = len(bb.labels)
    out = np.zeros(n, bool)
    for i in range(1, n):
        if bb.labels[i].chain == bb.labels[i - 1].chain:
            out[i] = np.linalg.norm(bb.xyz[i - 1, 2] - bb.xyz[i, 0]) <= PEPTIDE
    return out


def ks_energies(bb) -> np.ndarray:
    """E[acceptor, donor] (kcal/mol) for every pair, NaN where the donor has no amide H (a chain start, proline)."""
    n_, c_, o_ = bb.xyz[:, 0], bb.xyz[:, 2], bb.xyz[:, 3]
    linked = peptide_bonded(bb)
    h = np.full_like(n_, np.nan)
    for i in range(1, len(n_)):
        if linked[i] and bb.labels[i].name != "PRO":
            co = c_[i - 1] - o_[i - 1]
            h[i] = n_[i] + co / np.linalg.norm(co)  # 1.0 Å from N, along O->C of the preceding residue

    def dist(a, b):  # a[acceptor], b[donor] -> matrix
        return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)

    with np.errstate(divide="ignore", invalid="ignore"):
        e = KS_F * (1 / dist(o_, n_) + 1 / dist(c_, h) - 1 / dist(o_, h) - 1 / dist(c_, n_))
    np.fill_diagonal(e, np.nan)
    return e


def ks_bonds(bb) -> set[tuple[int, int]]:
    """(acceptor, donor) with E < -0.5 kcal/mol."""
    e = ks_energies(bb)
    acc, don = np.nonzero(np.nan_to_num(e, nan=0.0) < KS_CUT)
    return set(zip(acc.tolist(), don.tolist()))


def ks_bridges(bb, bonds: set[tuple[int, int]]) -> set[tuple[int, int, str]]:
    """(i, j, 'P' or 'A'), i < j, by the paper's patterns."""
    linked = peptide_bonded(bb)
    n = len(bb.labels)

    def hb(a, b):
        return (a, b) in bonds

    out = set()
    for i in range(1, n - 1):
        if not (linked[i] and linked[i + 1]):
            continue
        for j in range(i + 3, n - 1):
            if not (linked[j] and linked[j + 1]):
                continue
            if (hb(i - 1, j) and hb(j, i + 1)) or (hb(j - 1, i) and hb(i, j + 1)):
                out.add((i, j, "P"))
            elif (hb(i, j) and hb(j, i)) or (hb(i - 1, j + 1) and hb(j - 1, i + 1)):
                out.add((i, j, "A"))
    return out


def ks_ladders(bridges: set[tuple[int, int, str]]) -> list[set[tuple[int, int, str]]]:
    """Ladders: bridges of one type, each the next of the one before (parallel i+1, j+1; antiparallel i+1, j-1)."""
    left = set(bridges)
    out = []
    while left:
        start = min(left)
        ladder, frontier = {start}, [start]
        left.discard(start)
        while frontier:
            i, j, t = frontier.pop()
            step = 1 if t == "P" else -1
            for nxt in ((i + 1, j + step, t), (i - 1, j - step, t)):
                if nxt in left:
                    left.discard(nxt)
                    ladder.add(nxt)
                    frontier.append(nxt)
        out.append(ladder)
    return out


def ks_strand_graph(sses, bridges) -> dict[frozenset, int]:
    """Bridges between two foldmap strands (their residue spans), by strand-id pair: the strands a sheet makes
    adjacent by the paper's definition."""
    strands = [s for s in sses if s.kind == "E"]
    where = {r: s.id for s in strands for r in range(s.start, s.end + 1)}
    out: dict[frozenset, int] = {}
    for i, j, _ in bridges:
        a, b = where.get(i), where.get(j)
        if a and b and a != b:
            out[frozenset((a, b))] = out.get(frozenset((a, b)), 0) + 1
    return out


def ks_sheets(sses, bridges) -> list[set[str]]:
    """The paper's sheets (ladders joined by shared residues) as the sets of foldmap strands they hold."""
    import networkx as nx

    g = nx.Graph()
    g.add_nodes_from(s.id for s in sses if s.kind == "E")
    for pair in ks_strand_graph(sses, bridges):
        g.add_edge(*pair)
    return sorted((set(c) for c in nx.connected_components(g) if len(c) > 1), key=sorted)


# ------------------------------------------------------------------------------------------- canonical base pairs

HB_MAX = 3.5  # Å, a Watson-Crick-edge donor-acceptor distance
HB_NEED = 2  # of the edge's donor-acceptor pairs within HB_MAX
C1_RANGE = (9.5, 11.5)  # Å, C1'-C1' (10.5 in the standard pair, Olson et al. 2001)
LAMBDA = {"WC": (40.0, 70.0), "GU": (30.0, 80.0)}  # degrees; the wobble's two angles are unequal (about 40 and 70)
HAIRPIN_MIN = 3  # unpaired nucleotides a hairpin loop needs at least (Tinoco et al. 1973; Zuker & Stiegler 1981)
EDGE = {  # (purine, pyrimidine) -> (Watson-Crick-edge donor/acceptor pairs (purine atom, pyrimidine atom), kind)
    ("G", "C"): ((("N1", "N3"), ("N2", "O2"), ("O6", "N4")), "WC"),
    ("A", "U"): ((("N1", "N3"), ("N6", "O4")), "WC"),
    ("A", "T"): ((("N1", "N3"), ("N6", "O4")), "WC"),
    ("G", "U"): ((("N1", "O2"), ("O6", "N3")), "GU"),
    ("G", "T"): ((("N1", "O2"), ("O6", "N3")), "GU"),
}
GLYCOSIDIC = {"A": "N9", "G": "N9", "C": "N1", "U": "N1", "T": "N1"}
RING = ("N1", "C2", "N3", "C4", "C5", "C6")


class Nt(NamedTuple):
    key: str  # CHAIN:NUMBER
    base: str  # A C G U T: the standard base (modified nucleotides are not read)
    res: object  # gemmi residue
    c1: np.ndarray
    n: np.ndarray  # N9 or N1
    chain: str
    k: int  # its place among its chain's nucleotides


def _parent(res) -> str | None:
    info = gemmi.find_tabulated_residue(res.name)
    if info is None or not info.is_nucleic_acid() or not info.is_standard():
        return None
    code = info.one_letter_code.upper()
    return code if code in GLYCOSIDIC else None


def _xyz(res, name):
    a = res.find_atom(name, "*")
    return None if a is None else np.array(a.pos.tolist())


@functools.cache
def nucleotides(path) -> dict[str, Nt]:
    """The standard nucleotides of the deposited model 1 with a C1' and the N bonded to it, by key."""
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    out = {}
    for chain in st[0]:
        k = 0
        for res in chain:
            base = _parent(res)
            c1, n = _xyz(res, "C1'"), _xyz(res, GLYCOSIDIC.get(base or "", ""))
            if base and c1 is not None and n is not None:
                key = f"{chain.name}:{res.seqid.num}{res.seqid.icode.strip()}"
                out[key] = Nt(key, base, res, c1, n, chain.name, k)
                k += 1
    return out


def measure(a: Nt, b: Nt) -> dict | None:
    """The pair's geometry by the definitions above: kind, donor-acceptor distances, C1'-C1', the two glycosidic
    angles, cis; None when the two bases are not a canonical combination."""
    pur, pyr = (a, b) if a.base in "AG" else (b, a)
    rule = EDGE.get((pur.base, pyr.base))
    if rule is None:
        return None
    edge, kind = rule
    span = b.c1 - a.c1
    d = float(np.linalg.norm(span))
    u = span / d
    hbonds = []
    for p, q in edge:
        pp, qq = _xyz(pur.res, p), _xyz(pyr.res, q)
        if pp is not None and qq is not None:
            hbonds.append(float(np.linalg.norm(pp - qq)))
    ga, gb = a.n - a.c1, b.n - b.c1  # C1' -> N9/N1
    lam = (
        float(np.degrees(np.arccos(np.clip(ga @ u / np.linalg.norm(ga), -1, 1)))),
        float(np.degrees(np.arccos(np.clip(gb @ -u / np.linalg.norm(gb), -1, 1)))),
    )
    cis = float(np.cross(u, ga) @ np.cross(u, gb)) > 0  # both glycosidic bonds on one side of the C1'-C1' line
    return {"kind": kind, "hbonds": hbonds, "c1": d, "lambda": lam, "cis": cis}


def failed(m: dict) -> list[str]:
    """The definitions a measured pair fails ([]: a canonical pair)."""
    out = []
    if sum(h <= HB_MAX for h in m["hbonds"]) < HB_NEED:
        out.append("hbond")
    if not C1_RANGE[0] <= m["c1"] <= C1_RANGE[1]:
        out.append("c1")
    lo, hi = LAMBDA[m["kind"]]
    if not all(lo <= x <= hi for x in m["lambda"]):
        out.append("lambda")
    if not m["cis"]:
        out.append("trans")
    return out


@functools.cache
def canonical_pairs(path) -> dict[frozenset, str]:
    """Canonical pairs of the deposited model 1: {frozenset of the two keys: "WC" or "GU"}. The candidates are the
    nucleotide pairs whose C1' atoms lie within C1_RANGE's upper end (a k-d tree: a distance prefilter only)."""
    from scipy.spatial import cKDTree

    nts = list(nucleotides(path).values())
    out = {}
    if len(nts) < 2:
        return out
    for x, y in sorted(cKDTree(np.array([n.c1 for n in nts])).query_pairs(C1_RANGE[1])):
        a, b = nts[x], nts[y]
        if a.chain == b.chain and abs(a.k - b.k) <= HAIRPIN_MIN:
            continue
        m = measure(a, b)
        if m is not None and not failed(m):
            out[frozenset((a.key, b.key))] = m["kind"]
    return out


def plane_angle(a: Nt, b: Nt) -> float:
    """Degrees between the two bases' six-membered rings (least-squares planes)."""

    def normal(res):
        pts = np.array([_xyz(res, x) for x in RING if _xyz(res, x) is not None])
        return np.linalg.svd(pts - pts.mean(axis=0))[2][2]

    return float(np.degrees(np.arccos(min(1.0, abs(float(normal(a.res) @ normal(b.res)))))))


@functools.cache
def foldmap_pairs(path) -> dict[frozenset, str]:
    bb = load_backbone(path, "asu")
    na = load_nucleic(path, bb, "asu")

    def key(s, i):
        lab = na.strands[s].labels[i]
        return f"{lab.chain}:{lab.seq}{lab.icode}"

    return {frozenset((key(p.a, p.i), key(p.b, p.j))): p.kind for p in na.base_pairs}


def disagreements(path) -> dict[str, list[str]]:
    """Each pair foldmap and the kernel disagree on, by why: what the kernel's definitions fail for a pair only
    foldmap calls ("modified": a nucleotide the kernel does not read), "plane" for a pair only the kernel calls whose
    bases are more than 45 degrees apart (foldmap's PLANE_MAX), "kind" where both call it but differently,
    "unexplained" otherwise."""
    ours, theirs = foldmap_pairs(path), canonical_pairs(path)
    nts = nucleotides(path)
    out: dict[str, list[str]] = {}

    def add(why, pair):
        out.setdefault(why, []).append("/".join(sorted(pair)))

    for pair in sorted(ours.keys() - theirs.keys(), key=sorted):
        a, b = sorted(pair)
        if a not in nts or b not in nts:
            add("modified", pair)
            continue
        m = measure(nts[a], nts[b])
        for why in failed(m) if m is not None else ["unexplained"]:
            add(why, pair)
    for pair in sorted(theirs.keys() - ours.keys(), key=sorted):
        a, b = sorted(pair)
        add("plane" if plane_angle(nts[a], nts[b]) > 45.0 else "unexplained", pair)
    for pair in sorted(ours.keys() & theirs.keys(), key=sorted):
        if ours[pair] != theirs[pair]:
            add("kind", pair)
    return out


# ------------------------------------------------------------------------------------------------------- tests

PROTEIN_FIXTURES = [*PROTEINS, "6ZS5"]
NA_FIXTURES = sorted(p.name for p in NA.glob("*.cif.gz")) + ["../1LMB.cif"]
MF2 = (
    "T7 finding MF2: foldmap calls a base pair whose Watson-Crick-edge H-bonds reach 3.8 Å (nucleic.HBOND_MAX, "
    "nucleic.py:16; one stretched to 4.0, HBOND_STRETCH, nucleic.py:17); the kernel holds them to 3.5 Å"
)
MF3 = "T7 finding MF3: foldmap's C1'-C1' window is 9-12 Å (nucleic.C1_PAIR, nucleic.py:18); the kernel's 9.5-11.5"
MF4 = (
    "T7 finding MF4: foldmap drops a pair whose base planes are more than 45 degrees apart (nucleic.PLANE_MAX, "
    "nucleic.py:19, applied in pair_nucleotides, nucleic.py:156); the kernel has no planarity bound"
)


def _protein(name):
    return load_backbone(DATA / f"{name}.cif", "asu")


@pytest.mark.parametrize("name", PROTEIN_FIXTURES)
def test_hbonds_are_the_papers(name):
    """foldmap's H-bonds (dssp._hbond_set) are exactly the paper's E < -0.5 kcal/mol bonds but for the peptide
    neighbour's (CO(i-1) to NH(i)): the formula gives every one of those below -0.5, and the DSSP program skips them,
    as foldmap does. Neither mkdssp's two-best-acceptors rule nor its 9 Å CA cutoff drops a bond here."""
    bb = _protein(name)
    ks, ours = ks_bonds(bb), _hbond_set(bb)
    linked = peptide_bonded(bb)
    neighbours = {(i - 1, i) for i in range(1, len(linked)) if linked[i]}
    assert ours == ks - neighbours
    assert neighbours - {(i - 1, i) for i in range(len(linked)) if bb.labels[i].name == "PRO"} <= ks


@pytest.mark.parametrize("name", PROTEIN_FIXTURES)
def test_bridges_ladders_and_sheets_are_the_papers(name):
    """Bridges from the kernel's H-bonds are foldmap's; the paper's sheets (ladders joined by a shared residue) hold
    the strands foldmap's sheets hold (track and accessory strands), and every two strands next to each other on a
    sheet's track (and a closed sheet's last and first) share a bridge."""
    bb = _protein(name)
    d = assign_dssp(bb)
    ks = ks_bridges(bb, ks_bonds(bb))
    assert {(b.i, b.j, b.kind) for b in d.bridges} == ks
    for ladder in ks_ladders(ks):  # a ladder is one type throughout and steps one residue at a time
        assert len({t for _, _, t in ladder}) == 1
    sses = build_sses(bb, d.ss, bridges=d.bridges)
    sheets = build_sheets(sses, d.bridges)
    ours = sorted(({s.id for s in sh.strands} | {a.strand.id for a in sh.accessory} for sh in sheets), key=sorted)
    assert ours == ks_sheets(sses, ks)
    graph = ks_strand_graph(sses, ks)
    for sh in sheets:
        ids = [s.id for s in sh.strands]
        steps = list(zip(ids, ids[1:])) + ([(ids[-1], ids[0])] if sh.closed and len(ids) > 2 else [])
        assert all(graph.get(frozenset(p), 0) >= 1 for p in steps), (ids, graph)


@functools.cache
def _disagreements(name: str) -> dict[str, list[str]]:
    return disagreements(NA / name)


def test_base_pairs_agree_with_the_kernel():
    """Over every nucleic-acid fixture (and 1LMB's operator DNA): 820 pairs both call, 54 only foldmap, 5 only the
    kernel (docs/validation/metamorphic/report.md); held at no less than 0.92 of the union."""
    both = union = 0
    for name in NA_FIXTURES:
        ours, theirs = foldmap_pairs(NA / name).keys(), canonical_pairs(NA / name).keys()
        both += len(ours & theirs)
        union += len(ours | theirs)
    assert both >= 800 and both / union >= 0.92, (both, union)


@pytest.mark.parametrize("name", NA_FIXTURES)
def test_every_base_pair_disagreement_is_explained(name):
    """Each disagreement is one the findings name (MF2 H-bond window, MF3 C1'-C1' window, MF4 planarity; glycosidic
    angle only together with one of them) or the kernel's scope (a modified nucleotide it does not read); never a
    pair of another kind, nor one both call with different kinds."""
    found = _disagreements(name)
    assert not set(found) - {"hbond", "c1", "plane", "lambda", "modified"}, found
    for pair in found.get("lambda", []):  # never the only reason
        assert pair in found.get("hbond", []) + found.get("c1", []), pair


def _any(why: str) -> list[str]:
    return [f"{name}: {p}" for name in NA_FIXTURES for p in _disagreements(name).get(why, [])]


# The disagreements exactly as recorded (docs/validation/metamorphic/report.md): MF2 by fixture (count), MF3, MF4, the
# glycosidic angle (with MF2) and the modified nucleotides by pair. A new one, or one gone, fails
# test_the_findings_are_exactly_these, not only the xfails below.
FINDINGS = {
    "hbond": {"1AOI.cif.gz": 16, "1HYS.cif.gz": 2, "1KX5.cif.gz": 8, "1URN.cif.gz": 1, "3CRX.cif.gz": 9,
              "4OO8.cif.gz": 1, "4YBB.cif.gz": 7},
    "c1": ["1AOI.cif.gz: I:11/J:282", "1AOI.cif.gz: I:28/J:265", "1HYS.cif.gz: E:867/F:887", "4LD0.cif.gz: D:12/F:22",
           "4YBB.cif.gz: AA:112/AA:330", "4YBB.cif.gz: AA:582/AA:760"],
    "plane": ["1IHF.cif.gz: C:-37/E:37", "4YBB.cif.gz: AA:456/AA:476", "4YBB.cif.gz: AA:457/AA:475",
              "4YBB.cif.gz: AA:80/AA:89", "4YBB.cif.gz: AA:82/AA:87"],
    "lambda": ["1AOI.cif.gz: I:25/J:268"],
    "modified": ["1EHZ.cif.gz: A:10/A:25", "1EHZ.cif.gz: A:30/A:40", "1EHZ.cif.gz: A:49/A:65", "2OEU.cif.gz: A:1/B:20",
                 "2OEU.cif.gz: A:10/B:11", "4YBB-psi.cif.gz: DA:2590/DA:2604", "4YBB.cif.gz: AA:522/AA:527"],
}  # fmt: skip


def test_the_findings_are_exactly_these():
    by_fixture = {name: len(_disagreements(name).get("hbond", [])) for name in NA_FIXTURES}
    assert {k: v for k, v in by_fixture.items() if v} == FINDINGS["hbond"]
    assert sum(FINDINGS["hbond"].values()) == 44
    for why in ("c1", "plane", "lambda", "modified"):
        assert sorted(_any(why)) == sorted(FINDINGS[why]), why


def mf2_pairs(name: str) -> list[str]:
    """MF2's pairs on one fixture ("CHAIN:NUMBER/CHAIN:NUMBER")."""
    return list(_disagreements(name).get("hbond", []))


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=MF2)
def test_foldmap_pairs_hold_their_hbonds_within_3_5():
    assert not _any("hbond")


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=MF3)
def test_foldmap_pairs_span_the_kernels_c1_window():
    assert not _any("c1")


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=MF4)
def test_foldmap_finds_every_kernel_pair():
    assert not _any("plane")


RECORD = Path(__file__).resolve().parents[1] / "docs" / "validation" / "metamorphic" / "report.json"


@pytest.mark.skipif(not RECORD.is_file(), reason="the metamorphic record is not part of the sdist")
def test_the_record_names_these_findings():
    """The committed record (scripts/metamorphic_record.py) holds the disagreements the kernels find now."""
    import json

    record = json.loads(RECORD.read_text(encoding="utf-8"))["base_pairs"]["disagreements"]
    for why in ("hbond", "c1", "plane", "lambda", "modified"):
        assert sorted(record.get(why, [])) == sorted(_any(why)), why
