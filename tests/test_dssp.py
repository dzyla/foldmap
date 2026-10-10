from collections import Counter
from pathlib import Path

import gemmi

from foldmap.dssp import _hbond_set, _resolve, assign_dssp
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


def test_strand_outranks_310_and_pi_but_not_alpha():
    # DSSP priority: H > E > B > G > I
    codes = _resolve(8, alpha={0}, g={0, 1, 2}, pi={3, 4}, ladder={0, 2, 3}, bridge={1, 3, 5})
    assert codes[0] == "H"  # alpha beats strand
    assert codes[1] == "B"  # bridge beats 3-10
    assert codes[2] == "E"  # strand beats 3-10
    assert codes[3] == "E"  # strand beats pi
    assert codes[4] == "I"
    assert codes[5] == "B"
    assert codes[6] == "-" and codes[7] == "-"


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
