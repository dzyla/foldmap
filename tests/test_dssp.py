from collections import Counter
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.dssp import _hbond_set, assign_dssp
from foldmap.io import load_backbone
from foldmap.model import Bridge


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


def test_310_needs_every_residue_free():
    from foldmap.dssp import Ladder, _assign_codes

    # 3-turns at 0 and 1 would make residues 1-3 a 3-10 helix, but 2 and 3 are strand: no G at all
    ss = "".join(_assign_codes(12, [Ladder("A", [2, 3], [8, 9])], _starts(12, k3=[0, 1])))
    assert ss[1:4] == "TEE"


def test_pi_helix_takes_priority_over_alpha():
    from foldmap.dssp import _assign_codes

    starts = _starts(12, k4=[0, 1, 2, 3], k5=[0, 1])
    assert "".join(_assign_codes(12, [], starts))[1:7] == "IIIIIH"  # mkdssp >= 2.1
    assert "".join(_assign_codes(12, [], starts, prefer_pi=False))[1:7] == "HHHHHH"  # DSSP 1983


def test_pi_never_overwrites_a_strand():
    from foldmap.dssp import Ladder, _assign_codes

    ss = "".join(_assign_codes(12, [Ladder("A", [5, 6], [9, 10])], _starts(12, k4=[0, 1], k5=[0, 1])))
    assert ss[1:7] == "HHHHEE"  # the pi segment 1-5 would cover strand residue 5


def test_310_turns_do_not_join_alpha_helices():
    from foldmap.ss import build_sses

    bb = load_backbone(DATA / "8UUP.cif", "asu")
    ids = {s.id for s in build_sses(bb, assign_dssp(bb).ss)}
    assert {"B:456-472", "B:474-478"} <= ids and "B:456-478" not in ids  # mkdssp: T at 473


def test_ladders_join_across_a_bulge():
    from foldmap.dssp import _assign_codes, _ladders

    # antiparallel 10-12 / 40-38, then 14-15 / 36-35: one extra residue on each strand (13 and 37)
    bridges = [Bridge(10, 40, "A"), Bridge(11, 39, "A"), Bridge(12, 38, "A"), Bridge(14, 36, "A"), Bridge(15, 35, "A")]
    (lad,) = _ladders(bridges, np.ones(60, bool))
    assert lad.i == [10, 11, 12, 14, 15] and lad.j == [35, 36, 38, 39, 40]
    ss = "".join(_assign_codes(60, [lad], _starts(60)))
    assert ss[10:16] == "EEEEEE" and ss[35:41] == "EEEEEE"


def test_bulge_join_stops_at_chain_breaks_and_wide_gaps():
    from foldmap.dssp import _ladders

    bridges = [Bridge(10, 40, "A"), Bridge(11, 39, "A"), Bridge(12, 38, "A"), Bridge(14, 36, "A"), Bridge(15, 35, "A")]
    prev = np.ones(60, bool)
    prev[13] = False
    assert len(_ladders(bridges, prev)) == 2
    far = [Bridge(10, 40, "A"), Bridge(11, 39, "A"), Bridge(20, 30, "A"), Bridge(21, 29, "A")]
    assert len(_ladders(far, np.ones(60, bool))) == 2


def test_lone_bridge_joined_across_a_bulge_is_strand():
    from foldmap.dssp import _assign_codes, _ladders

    lads = _ladders([Bridge(10, 40, "A"), Bridge(12, 38, "A"), Bridge(13, 37, "A")], np.ones(60, bool))
    assert len(lads) == 1
    ss = "".join(_assign_codes(60, lads, _starts(60)))
    assert ss[10:14] == "EEEE" and ss[37:41] == "EEEE"


def test_bulged_strands_stay_whole():
    from foldmap.ss import build_sses

    bb = load_backbone(DATA / "5NKT.cif", "asu")
    assert sum(s.kind == "E" for s in build_sses(bb, assign_dssp(bb).ss)) == 12  # mkdssp: 12 strands


def test_hbond_energy_is_floored_like_dssp():
    from foldmap.dssp import _hbond_energy

    e = _hbond_energy([0.3, 3.0, 5.0], [2.0, 2.1, 5.0], [1.9, 2.0, 0.6], [3.0, 3.0, 0.6])
    assert e[0] == pytest.approx(-9.9)  # an atom pair closer than 0.5 Å
    assert e[1] == pytest.approx(-0.664, abs=1e-3)  # an ordinary H-bond, unchanged
    assert e[2] == pytest.approx(-9.9)  # -81.8 kcal/mol, floored


def test_coincident_atoms_give_the_floor_without_warnings():
    import warnings

    from foldmap.dssp import _hbond_energy

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _hbond_energy([0.0], [0.0], [0.0], [0.0])[0] == pytest.approx(-9.9)


def _energy_distances(e):
    """O-N, C-H, O-H and C-N distances whose energy is e: C-H and C-N cancel, O-H is 2 Å."""
    from foldmap.dssp import _Q1Q2F

    return [1 / (e / _Q1Q2F + 0.5)], [3.0], [2.0], [3.0]


def test_hbond_energy_is_rounded_to_a_thousandth_like_mkdssp():
    from foldmap.dssp import _HBOND_E, _hbond_energy

    assert _hbond_energy(*_energy_distances(-0.5003))[0] == -0.5  # mkdssp: not a bond (the test is < -0.5)
    assert not _hbond_energy(*_energy_distances(-0.5003))[0] < _HBOND_E
    assert _hbond_energy(*_energy_distances(-0.5006))[0] == -0.501  # a bond
    assert _hbond_energy(*_energy_distances(-0.5006))[0] < _HBOND_E


def test_equal_energies_keep_the_lower_acceptors(monkeypatch):
    import foldmap.dssp as dssp
    from foldmap.model import Backbone, ResLabel

    # eight residues packed within 9 Å, every CO-NH pair at the same energy: each NH keeps its two
    # lowest-numbered acceptors, as mkdssp does (strict < while it walks acceptors in order)
    rng = np.random.default_rng(0)
    n = 8
    xyz = rng.uniform(0, 4, (n, 4, 3))
    bb = Backbone([ResLabel("A", k + 1, "", "ALA") for k in range(n)], xyz, np.r_[False, np.ones(n - 1, bool)])
    monkeypatch.setattr(dssp, "_hbond_energy", lambda d_on, *_: np.full(len(d_on), -1.0))
    bonds = dssp._hbond_set(bb)
    for donor in range(1, n):
        acceptors = [a for a in range(n) if a not in (donor, donor - 1)][:2]
        assert {a for a, d in bonds if d == donor} == set(acceptors), donor
