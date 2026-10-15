"""Biological assemblies built from the file's operators, cubic point groups, and the highlight mode."""

from pathlib import Path

import numpy as np
import pytest

from foldmap.io import load_backbone

DATA = Path(__file__).parent / "data"
APO = DATA / "assemblies" / "7A4M.cif"  # mouse apoferritin: one subunit deposited, 24 by the O operators


def test_auto_builds_the_24mer_and_remembers_the_deposited_subunit():
    bb = load_backbone(APO)
    assert len({l.chain for l in bb.labels}) == 24
    assert len(bb.asu_chains) == 1 and bb.asu_chains <= {l.chain for l in bb.labels}


def test_asu_keeps_the_file_as_deposited():
    bb = load_backbone(APO, assembly="asu")
    assert {l.chain for l in bb.labels} == {"A"} and bb.asu_chains == {"A"}


def test_unknown_assembly_is_an_error():
    with pytest.raises(ValueError, match="assembly"):
        load_backbone(APO, assembly="7")


def test_single_copy_files_are_unchanged(ubq):
    bb = load_backbone(ubq)
    assert {l.chain for l in bb.labels} == {"A"} and bb.asu_chains == {"A"}


def test_apoferritin_is_octahedral():
    from foldmap.symmetry import detect_symmetry

    bb = load_backbone(APO)
    s = detect_symmetry(bb)
    assert s.kind == "O" and s.label == "O" and s.n == 4
    assert len(s.protomers) == 4 and all(len(p) == 6 for p in s.protomers)


def test_tetrahedral_group_is_recognised(tmp_path):
    from itertools import product

    from foldmap.symmetry import detect_symmetry
    from test_symmetry_scenarios import build, rot

    gens = [rot((1, 1, 1), 120), rot((0, 0, 1), 180)]
    group = [np.eye(3)]
    while True:  # close the group under the generators
        new = [g @ h for g, h in product(gens, group)]
        add = [m for m in new if not any(np.allclose(m, x) for x in group)]
        if not add:
            break
        group += add
    assert len(group) == 12
    path = build(tmp_path, [(m, np.zeros(3)) for m in group], offset=(16.0, 7.0, 2.0))
    s = detect_symmetry(load_backbone(path))
    assert s.kind == "T" and s.n == 3


def test_highlight_keeps_the_asu_coloured_and_greys_the_mates(tmp_path):
    from matplotlib.colors import to_hex

    from foldmap.cli import make_figure, make_layout
    from foldmap.render import MATE_GREY
    from foldmap.style import Style

    lay, sses, bb = make_layout(APO)
    fig = make_figure(APO, look=Style(highlight="asu"))
    ax = fig.axes[0]
    fills = {
        a.get_gid().split(":", 1)[1]: to_hex(a.get_facecolor())
        for a in ax.patches
        if (a.get_gid() or "").startswith("helix:")
    }
    asu = [s.id for s in sses if s.chain in bb.asu_chains]
    mates = [s.id for s in sses if s.chain not in bb.asu_chains]
    assert asu and all(fills[k] != MATE_GREY for k in asu)
    assert mates and all(fills[k] == MATE_GREY for k in mates)


def test_highlight_by_chain_list():
    from matplotlib.colors import to_hex

    from foldmap.cli import make_figure
    from foldmap.render import MATE_GREY
    from foldmap.style import Style

    fig = make_figure(DATA / "2HHB.cif", look=Style(highlight="A,B"))
    for a in fig.axes[0].patches:
        gid = a.get_gid() or ""
        if gid.startswith("helix:"):
            chain = gid.split(":")[1]
            assert (to_hex(a.get_facecolor()) == MATE_GREY) == (chain in "CD"), gid


def test_bad_highlight_value_is_rejected():
    from foldmap.style import resolve_style

    with pytest.raises(ValueError, match="highlight"):
        resolve_style(overrides=["highlight=A;B"])


def test_cli_assembly_and_highlight(tmp_path, capsys):
    from foldmap.cli import main

    out = tmp_path / "apo.svg"
    assert main(["plot", str(APO), "-o", str(out), "--assembly", "asu"]) == 0
    assert main(["summary", str(APO)]) == 0
    text = capsys.readouterr().out
    assert "24 chain(s)" in text and "symmetry: O" in text


def _model(chains):
    """A gemmi Model of one-atom residues: {chain name: [xyz, ...]}."""
    import gemmi

    model = gemmi.Model("1")
    for name, xyz in chains.items():
        chain = gemmi.Chain(name)
        for k, p in enumerate(xyz):
            res = gemmi.Residue()
            res.name, res.seqid = "ALA", gemmi.SeqId(k + 1, " ")
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = "CA", gemmi.Element("C"), gemmi.Position(*p)
            res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    return model


def test_a_copy_lying_on_a_chain_already_present_is_dropped():
    """A copy with >= 90 % of its atoms within 1 Å of one chain already present is the same molecule seen twice
    (an operator mapping a chain on a special position onto itself or its partner); a copy that only partly
    overlaps, or lies on the union of two chains, is kept."""
    from foldmap.io import _drop_coincident

    rng = np.random.default_rng(0)
    a = rng.uniform(0, 40, (20, 3))
    b = a + 20.0
    on_a = a + 0.4  # every atom 0.69 Å from its twin
    half = np.vstack([a[:10] + 0.4, a[10:] + 6.0])  # half of it coincides: a real neighbour, kept
    split = np.vstack([a[:9] + 0.4, b[9:] + 0.4])  # 45 % on A and 55 % on B: on no single chain
    model = _model({"A1": a, "B1": b, "A2": on_a, "C2": half, "D2": split})
    assert _drop_coincident(model) == ["A2"]
    assert [c.name for c in model] == ["A1", "B1", "C2", "D2"]


def test_holliday_junction_assembly_keeps_each_dna_strand_once():
    """3CRX: the two-fold operator lays a copy of every DNA strand on a deposited one (C2 on E1, D2 on F1, ...).
    The assembly keeps the four Cre and the four strands, which pair into the junction's four arms."""
    from scipy.spatial import cKDTree

    from foldmap.io import load_nucleic, read_model

    path = DATA / "na" / "3CRX.cif.gz"
    model = read_model(path)
    bb = load_backbone(path, model=model)
    assert len({l.chain for l in bb.labels}) == 4  # four Cre
    na = load_nucleic(path, bb, model=model)
    assert len(na.strands) == 4
    for i, s in enumerate(na.strands):
        for t in na.strands[:i]:
            d, _ = cKDTree(t.p).query(s.p)
            assert (d < 1.5).mean() < 0.5, (s.chain, t.chain)
    assert len(na.structure.stems) == 4


def test_an_assembly_whose_every_added_copy_is_dropped_keeps_the_deposited_names(tmp_path):
    """auto: when every copy the assembly adds lies exactly on a chain already present (C5 drops them all), the
    result is the molecule as deposited, under its deposited chain names (A, not A1)."""
    import gemmi

    from foldmap.io import read_model

    st = gemmi.read_structure(str(DATA / "na" / "1EHZ.cif.gz"))
    gen = gemmi.Assembly.Gen()
    gen.chains = [c.name for c in st[0]]
    gen.subchains = sorted({r.subchain for c in st[0] for r in c})
    for name in ("1", "2"):  # the second operator maps the tRNA onto itself
        op = gemmi.Assembly.Operator()
        op.name, op.transform = name, gemmi.Transform()
        gen.operators.append(op)
    assembly = gemmi.Assembly("1")
    assembly.generators.append(gen)
    st.assemblies.clear()
    st.assemblies.append(assembly)
    path = tmp_path / "twice.cif"
    st.make_mmcif_document().write_file(str(path))
    model, deposited = read_model(path)
    assert [c.name for c in model] == ["A"] and deposited == {"A"}
    bb = load_backbone(path, model=(model, deposited))
    assert bb.asu_chains <= {"A"}


def test_a_deposited_chain_is_never_dropped_and_one_operator_never_drops_its_own_copies():
    """Only copies made by another operator are the same molecule seen twice. Chains of one operator that lie on
    each other were deposited that way (alternative conformations: 1C7Y's DNA chains B/D, C/E ...); the
    identity copy of a deposited chain is kept whatever it lies on."""
    from foldmap.io import _drop_coincident

    rng = np.random.default_rng(1)
    a = rng.uniform(0, 40, (20, 3))
    model = _model({"A1": a, "B1": a + 0.3, "A2": a + 100.0, "B2": a + 100.3, "C2": a + 0.2})
    assert _drop_coincident(model, keep={"A1", "B1"}, tags=[1, 1, 2, 2, 2]) == ["C2"]
    assert [c.name for c in model] == ["A1", "B1", "A2", "B2"]
    model = _model({"A1": a, "B1": a + 0.3})
    assert _drop_coincident(model, keep={"A1", "B1"}) == []  # without operator tags: the identity copies stay


def _with_assembly(st, ops, tmp_path, name="asm.cif"):
    """Write st with one assembly applying the (rotation, translation) operators to every chain."""
    import gemmi

    gen = gemmi.Assembly.Gen()
    gen.subchains = sorted({r.subchain for c in st[0] for r in c})
    for k, (rot, shift) in enumerate(ops):
        op = gemmi.Assembly.Operator()
        op.name = str(k + 1)
        op.transform = gemmi.Transform(gemmi.Mat33(rot), gemmi.Vec3(*shift))
        gen.operators.append(op)
    assembly = gemmi.Assembly("1")
    assembly.generators.append(gen)
    st.assemblies.clear()
    st.assemblies.append(assembly)
    path = tmp_path / name
    st.make_mmcif_document().write_file(str(path))
    return path


def test_chains_deposited_on_each_other_are_kept_in_every_copy(tmp_path):
    """1C7Y-like: two DNA chains deposited on top of each other (B a copy of A), an assembly of two copies far
    apart. All four chains stay (C5 dropped the deposited B1, and B2 with it)."""
    import gemmi

    from foldmap.io import read_model

    st = gemmi.read_structure(str(DATA / "na" / "1EHZ.cif.gz"))
    st.setup_entities()
    twin = st[0]["A"].clone()
    twin.name = "B"
    for r in twin:
        r.subchain = "Z"
    st[0].add_chain(twin)
    eye = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    path = _with_assembly(st, [(eye, (0, 0, 0)), (eye, (150.0, 0, 0))], tmp_path)
    model, asu = read_model(path)
    assert [c.name for c in model] == ["A1", "B1", "A2", "B2"]
    assert asu == {"A1", "B1"}


def test_holliday_junction_drops_exactly_the_copies_of_the_second_operator():
    from foldmap.io import read_model

    model, asu = read_model(DATA / "na" / "3CRX.cif.gz")
    assert sorted(c.name for c in model) == ["A1", "A2", "B1", "B2", "C1", "D1", "E1", "F1"]
    assert asu == {"A1", "B1", "C1", "D1", "E1", "F1"}


def _copies_in_asu(tmp_path, n, *, keep=None, name="ncs.cif"):
    """1URN-like: the U1A–RNA complex of the 1URN fixture (RNA R, protein C) deposited `n` times in the asymmetric
    unit, each copy turned 120° about z and raised 25 Å (a screw, as crystal packing lays them), one assembly per
    copy built by the identity operator (assembly 1 = the first copy). keep: chain names of the fixture kept."""
    import gemmi

    st = gemmi.read_structure(str(DATA / "na" / "1URN.cif.gz"))
    st.setup_entities()
    model = st[0]
    for c in [c.name for c in model if keep is not None and c.name not in keep]:
        model.remove_chain(c)
    base = [c.clone() for c in model]
    entity_of = {sub: e for e in st.entities for sub in e.subchains}
    names = {"R": "RSTU", "C": "CDEF"}
    assemblies = []
    for k in range(n):
        turn = np.radians(120.0 * k)
        rot = [[np.cos(turn), -np.sin(turn), 0.0], [np.sin(turn), np.cos(turn), 0.0], [0.0, 0.0, 1.0]]
        subchains = []
        for c in base:
            copy = c.clone()
            copy.name = names[c.name][k]
            for r in copy:
                old = r.subchain
                r.subchain = f"{old}{k}" if k else old
                if k and r.subchain not in entity_of[old].subchains:
                    entity_of[old].subchains = [*entity_of[old].subchains, r.subchain]  # a copy is no new entity
                subchains.append(r.subchain)
                for a in r:
                    p = np.array(a.pos.tolist())
                    a.pos = gemmi.Position(*(np.array(rot) @ p + (0.0, 0.0, 25.0 * k)))
            if k:
                model.add_chain(copy)
            else:
                model.remove_chain(c.name)
                model.add_chain(copy)
        gen = gemmi.Assembly.Gen()
        gen.subchains = sorted(set(subchains))
        op = gemmi.Assembly.Operator()
        op.name, op.transform = "1", gemmi.Transform()
        gen.operators.append(op)
        assembly = gemmi.Assembly(str(k + 1))
        assembly.generators.append(gen)
        assemblies.append(assembly)
    st.assemblies.clear()
    for a in assemblies:
        st.assemblies.append(a)
    path = tmp_path / name
    st.make_mmcif_document().write_file(str(path))
    return path


def test_auto_draws_one_complex_when_the_asymmetric_unit_holds_several(tmp_path):
    """4OO8, 1URN: the asymmetric unit holds two or three copies of a protein–nucleic-acid complex and assembly 1
    is one of them (a strict subset, the identity operator). auto draws that biological assembly, under the
    deposited chain names; asu still draws every copy."""
    from foldmap.io import read_model

    path = _copies_in_asu(tmp_path, 3)
    model, deposited = read_model(path)
    assert [c.name for c in model] == ["R", "C"] and deposited == {"R", "C"}
    model, deposited = read_model(path, "asu")
    assert sorted(c.name for c in model) == ["C", "D", "E", "R", "S", "T"]
    bb = load_backbone(path)
    assert {l.chain for l in bb.labels} == {"C"} and bb.asu_chains == {"C"}


def test_crystal_packing_of_copies_in_the_asymmetric_unit_is_not_reported_as_symmetry(tmp_path):
    """1URN reported its three U1A–RNA copies, related by crystal packing (a 120° screw), as a helical assembly.
    Symmetry is looked for among the copies the chosen assembly holds: under auto, one complex, no symmetry."""
    from foldmap.cli import make_layout
    from foldmap.symmetry import detect_symmetry

    path = _copies_in_asu(tmp_path, 3)
    asu = load_backbone(path, assembly="asu")
    assert detect_symmetry(asu) is None  # three copies on a screw: crystal packing, not a filament
    assert detect_symmetry(asu, request="helical").kind == "H"  # still found when asked for
    assert detect_symmetry(load_backbone(path)) is None
    lay, sses, bb = make_layout(path)
    assert {s.chain for s in lay.nucleic.strands} == {"R"} and {s.chain for s in sses} == {"C"}


def test_one_protein_complex_of_several_draws_assembly_one(tmp_path):
    """Protein-only now as for nucleic acid (COMPARE #4): two copies of one complex deposited together, assembly 1
    the first of them (the identity operator, a strict subset): auto draws that one; asu still draws both."""
    from foldmap.io import read_model

    path = _copies_in_asu(tmp_path, 2, keep={"C"})
    model, deposited = read_model(path)
    assert [c.name for c in model] == ["C"] and deposited == {"C"}
    model, _ = read_model(path, "asu")
    assert sorted(c.name for c in model) == ["C", "D"]


def test_one_complex_of_several_counts_copies_of_a_protein_by_sequence():
    """C30 (random-100): 4UM3 deposits eight pentamers (40 chains); assembly 1 is one of them, P-T. Its chains P and
    S are entities of their own (a deletion and a point mutant of the same 229-residue protein), each deposited
    once, so counting by entity kept the 40-chain asymmetric unit (301 s, past the app's limit). Entities whose
    sequences align at least 90 % identical count as one protein: auto draws the pentamer and says so."""
    from foldmap.io import read_model

    path = DATA / "r100" / "4UM3-two.cif.gz"
    model = read_model(path)
    assert sorted(c.name for c in model[0]) == list("PQRST")
    assert model.several == {"assembly": "1", "copies": 2, "chains": 10, "residues": model.several["residues"]}
    assert sorted(c.name for c in read_model(path, "asu")[0]) == list("KLMNOPQRST")


def test_kelch_dimer_asu_draws_its_biological_assembly():
    from foldmap.cli import make_layout
    from foldmap.io import read_model
    from foldmap.symmetry import detect_symmetry

    model, _ = read_model(DATA / "4CH9.cif")
    assert [c.name for c in model] == ["B", "D"]  # assembly 1: one Kelch domain and its WNK4 peptide
    lay, sses, bb = make_layout(DATA / "4CH9.cif")
    assert {l.chain for l in bb.labels} == {"B", "D"} and detect_symmetry(bb) is None


def test_nucleic_acid_outside_assembly_one_does_not_choose_it(tmp_path):
    """Assembly 1 is chosen only when it holds every kind of polymer the asymmetric unit holds: two protein copies
    whose assembly 1 is one of them, beside an RNA that no assembly holds, keep the asymmetric unit (the RNA is
    never silently dropped)."""
    import gemmi

    from foldmap.io import read_model

    path = _copies_in_asu(tmp_path, 2, keep={"C"})
    st = gemmi.read_structure(str(path))
    rna = gemmi.read_structure(str(DATA / "na" / "1URN.cif.gz"))[0]["R"].clone()
    st[0].add_chain(rna)
    out = tmp_path / "with_rna.cif"
    st.make_mmcif_document().write_file(str(out))
    model, deposited = read_model(out)
    assert sorted(c.name for c in model) == ["C", "D", "R"] and deposited == {"C", "D", "R"}


def test_a_chromophore_is_a_protein_residue_for_the_kind_gate():
    """GFP's CRO66 has no N/CA/C names and no C1': it is neither a nucleotide nor a new polymer kind."""
    from foldmap.io import _is_nucleotide, _polymer_kinds, read_model

    model, _ = read_model(DATA / "1EMA.cif", "asu")
    chain = model["A"]
    k = next(i for i, r in enumerate(chain) if r.name == "CRO")
    assert not _is_nucleotide(chain, k)
    assert _polymer_kinds(model) == {"protein"}


def test_dna_outside_assembly_one_keeps_the_asymmetric_unit(tmp_path):
    """The same gate for a protein–DNA entry: assembly 1 holds one protein copy and no DNA, so auto keeps the
    asymmetric unit, DNA included."""
    import gemmi

    from foldmap.io import read_model

    path = _copies_in_asu(tmp_path, 2, keep={"C"})
    st = gemmi.read_structure(str(path))
    dna = gemmi.read_structure(str(DATA / "na" / "1KX5.cif.gz"))[0]["I"].clone()
    st[0].add_chain(dna)
    st.setup_entities()
    out = tmp_path / "with_dna.cif"
    st.make_mmcif_document().write_file(str(out))
    model, deposited = read_model(out)
    assert sorted(c.name for c in model) == ["C", "D", "I"] and deposited == {"C", "D", "I"}
    model, _ = read_model(out, "asu")
    assert sorted(c.name for c in model) == ["C", "D", "I"]


def test_a_ligand_outside_assembly_one_is_not_a_polymer_kind(tmp_path):
    """A nucleotide ligand (ATP) beside two protein copies is no DNA or RNA: assembly 1 still holds every kind of
    polymer the asymmetric unit holds, and auto draws it."""
    import gemmi

    from foldmap.io import _polymer_kinds, read_model

    path = _copies_in_asu(tmp_path, 2, keep={"C"})
    st = gemmi.read_structure(str(path))
    atp = gemmi.read_structure(str(DATA / "na" / "1ATP.cif.gz"))
    res = next(r for c in atp[0] for r in c if r.name == "ATP")
    ligand = gemmi.Chain("L")
    ligand.add_residue(res.clone())
    st[0].add_chain(ligand)
    st.setup_entities()
    out = tmp_path / "with_atp.cif"
    st.make_mmcif_document().write_file(str(out))
    model, _ = read_model(out, "asu")
    assert "L" in [c.name for c in model] and _polymer_kinds(model) == {"protein"}
    model, deposited = read_model(out)
    assert [c.name for c in model] == ["C"] and deposited == {"C"}


def test_a_peptide_is_not_a_copy_of_its_receptor():
    """C30 fix round 1: a global alignment scatters a short sequence over a long one, so a 30-residue peptide
    counted as the same protein as a 300-residue receptor (and the bound peptide of an asymmetric unit was merged
    with it): the shorter must also be at least 0.8 as long as the longer. A mutant 90 % as long still merges."""
    import gemmi

    from foldmap.io import _same_proteins

    seq = [
        ("ALA", "GLY", "LEU", "SER", "VAL", "LYS", "GLU", "ASP", "THR", "PRO")[(7 * i + i // 3) % 10]
        for i in range(300)
    ]

    def entity(name, residues):
        e = gemmi.Entity(name)
        e.entity_type = gemmi.EntityType.Polymer
        e.polymer_type = gemmi.PolymerType.PeptideL
        e.full_sequence = list(residues)
        return e

    st = gemmi.Structure()
    for e in (entity("1", seq), entity("2", seq[100:130]), entity("3", seq[:270])):
        st.entities.append(e)
    assert _same_proteins(st) == {"1": "1", "2": "2", "3": "1"}


def test_chains_that_pair_are_one_object_and_auto_keeps_both(tmp_path):
    """3OWW ruling: the two RNAs of a glycine riboswitch dimer are identical entities deposited once each, and
    assembly 1 is chain A alone (assembly 3 the dimer); counting copies by sequence drew A alone. Chains held together
    by base pairs or β-strands between them are one object, so auto keeps A and B."""
    from foldmap.io import read_model

    model = read_model(DATA / "r100" / "3OWW-AB.cif.gz")
    assert sorted(c.name for c in model[0]) == ["A", "B"] and model.several is None


def _chain(name, resname, atoms, n, shift):
    """A synthetic chain: n residues `resname`, each holding atoms (name, dx, dy), 4 A apart along x from `shift`."""
    import gemmi

    ch = gemmi.Chain(name)
    for i in range(n):
        r = gemmi.Residue()
        r.name, r.seqid = resname, gemmi.SeqId(i + 1, " ")
        for an, dx, dy in atoms:
            a = gemmi.Atom()
            a.name, a.element = an, gemmi.Element(an[0])
            a.pos = gemmi.Position(4.0 * i + dx + shift, dy, 0.0)
            r.add_atom(a)
        ch.add_residue(r)
    return ch


def test_protein_copies_in_contact_are_not_paired_but_base_pairs_are():
    """Fix round 2: 7OKA's ten trimers touch in up to 10 residue pairs of backbone N and O (crystal and assembly
    contacts), which merged them into one object (and drew chain A alone). Only nucleotides' base pairs count: two
    protein chains with 20 backbone N...O contacts are not paired; two RNAs with 4 base-atom contacts (3OWW's measured
    count; PAIRED_MIN is 2) are, and one contact pair is not."""
    from foldmap.io import PAIRED_MIN, _paired

    prot = [("N", 0.0, 0.0), ("O", 1.0, 0.0)]
    assert not _paired([_chain("A", "ALA", prot, 20, 0.0)], [_chain("B", "ALA", prot, 20, 0.2)])
    base = [("N1", 0.0, 0.0), ("O6", 1.0, 0.0)]
    rna = lambda name, shift, n: _chain(name, "G", base, n, shift)  # noqa: E731
    assert PAIRED_MIN == 2 and _paired([rna("A", 0.0, 4)], [rna("B", 0.2, 4)])
    assert not _paired([rna("A", 0.0, 1)], [rna("B", 0.2, 1)]) and PAIRED_MIN > 1
