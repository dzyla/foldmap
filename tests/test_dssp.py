from collections import Counter
from pathlib import Path

import gemmi

from foldmap.dssp import _hbond_set, assign_dssp
from foldmap.io import load_backbone


def test_ubiquitin_helix_and_strands(ubq):
    bb = load_backbone(ubq)
    res = assign_dssp(bb)
    assert len(res.ss) == len(bb)
    assert set(res.ss[22:34]) == {"H"}  # helix 23-34
    assert set(res.ss[2:7]) == {"E"}  # beta1 3-7
    assert set(res.ss[66:71]) == {"E"}  # beta5 67-71
    assert res.ss.count("H") == 12


def _deposited_strand_residues(path):
    st = gemmi.read_structure(str(path))
    out = set()
    for sheet in st.sheets:
        for s in sheet.strands:
            for k in range(s.start.res_id.seqid.num, s.end.res_id.seqid.num + 1):
                out.add((s.start.chain_name, k))
    return out


def test_agrees_with_deposited_strands(zs5):
    bb = load_backbone(zs5)
    res = assign_dssp(bb)
    mine = {(bb.labels[i].chain, bb.labels[i].seq) for i, c in enumerate(res.ss) if c == "E"}
    dep = _deposited_strand_residues(zs5)
    assert len(mine & dep) / len(mine) >= 0.95  # what we call a strand, the file agrees with
    assert len(mine & dep) / len(dep) >= 0.85  # and we recover most of the file's strands


def test_cross_chain_bridges_found(zs5):
    bb = load_backbone(zs5)
    res = assign_dssp(bb)
    cross = [b for b in res.bridges if bb.labels[b.i].chain != bb.labels[b.j].chain]
    assert len(cross) >= 10  # linker strands complete the neighbouring sheet


def test_bridge_invariants(ubq):
    res = assign_dssp(load_backbone(ubq))
    assert all(b.j >= b.i + 3 and b.kind in "PA" for b in res.bridges)


def test_turns_marked_like_dssp(ubq):
    ss = assign_dssp(load_backbone(ubq)).ss
    assert "T" in ss
    assert ss.replace("T", "-").count("H") == ss.count("H")  # turns never displace other codes


DATA = Path(__file__).parent / "data"


def test_insertion_codes_do_not_change_the_assignment(tmp_path):
    st = gemmi.read_structure(str(DATA / "1UBQ.cif"))
    for k, res in enumerate(list(st[0][0])[:17]):  # residues 1-17 become 1A, 1B, ... 1Q (antibody-style)
        res.seqid = gemmi.SeqId(1, chr(ord("A") + k))
    path = tmp_path / "icodes.cif"
    st.make_mmcif_document().write_file(str(path))
    plain = assign_dssp(load_backbone(DATA / "1UBQ.cif", "asu")).ss
    assert assign_dssp(load_backbone(path, "asu")).ss == plain


def _index(bb, chain, seq):
    return next(k for k, l in enumerate(bb.labels) if l.chain == chain and l.seq == seq)


def test_an_acceptor_may_take_more_than_two_hbonds():
    # mkdssp keeps the two best acceptors per NH and never limits a CO; 1LMB O(37) -> N(41) is one it keeps
    bb = load_backbone(DATA / "1LMB.cif", "asu")
    bonds = _hbond_set(bb)
    assert max(Counter(a for a, _ in bonds).values()) == 3
    assert max(Counter(d for _, d in bonds).values()) <= 2
    assert (_index(bb, "3", 37), _index(bb, "3", 41)) in bonds


def _starts(n, **turns):
    """Turn starts for _assign_codes: _starts(12, k4=[0, 1]) has 4-turns starting at residues 0 and 1."""
    return {k: [i in set(turns.get(f"k{k}", ())) for i in range(n)] for k in (3, 4, 5)}


def test_alpha_replaces_strand_and_bridge():
    from foldmap.dssp import Ladder, _assign_codes

    ss = "".join(_assign_codes(12, [Ladder("A", [2, 3], [8, 9]), Ladder("P", [5], [10])], _starts(12, k4=[1, 2])))
    assert ss[2:6] == "HHHH" and ss[8:11] == "EEB"


def test_bridges_need_unbroken_neighbours():
    bb = load_backbone(DATA / "6Y7S.cif", "asu")
    for i, j, _ in assign_dssp(bb).bridges:
        assert bb.prev[i] and bb.nxt[i] and bb.prev[j] and bb.nxt[j], (bb.labels[i], bb.labels[j])


def test_peptide_too_short_for_any_turn(tmp_path):
    from helpers import tripeptide

    res = assign_dssp(load_backbone(tripeptide(DATA / "1UBQ.cif", tmp_path)))
    assert res.ss == "---" and res.bridges == []
