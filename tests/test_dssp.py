import gemmi

from foldmap.dssp import _resolve, assign_dssp
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
