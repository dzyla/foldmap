# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Read PDB/mmCIF into a Backbone (first model, first altloc, amino acids only)."""

from __future__ import annotations

from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .model import Backbone, Ligand, Links, NAStrand, Nucleic, ResLabel, nucleic_from_pairs
from .nucleic import NAGap, g_tetrads, pair_nucleotides, parent_base

_BB_ATOMS = ("N", "CA", "C", "O")
_PEPTIDE_BOND_MAX = 2.5  # Å, C(i-1)-N(i), as in DSSP
_CO_BOND = 1.231  # Å, carbonyl C=O
_CO_MAX = 1.5  # Å, the farthest an oxygen still counts as the carbonyl one


def _substitute_backbone(res, before, after):
    """The N, CA, C, O stand-ins of a polymer residue with no N/CA/C of its own (GFP's chromophore CRO66, which
    fuses three residues into one ring system): its N bonded to the C of the residue before, its C bonded to the N
    of the one after, its CA the other carbon nearest that N, its O the oxygen on that C (None: rebuilt later).
    None when the residue is not between two residues it is peptide-bonded to."""
    if before is None or after is None:
        return None
    prev_c, next_n = before.find_atom("C", "*"), after.find_atom("N", "*")
    if prev_c is None or next_n is None:
        return None
    elements = {el: [a for a in res if a.element.name == el] for el in ("N", "C", "O")}
    n_in = min(elements["N"], key=lambda a: a.pos.dist(prev_c.pos), default=None)
    c_out = min(elements["C"], key=lambda a: a.pos.dist(next_n.pos), default=None)
    if n_in is None or c_out is None:
        return None
    if n_in.pos.dist(prev_c.pos) >= _PEPTIDE_BOND_MAX or c_out.pos.dist(next_n.pos) >= _PEPTIDE_BOND_MAX:
        return None
    ca = min((a for a in elements["C"] if a is not c_out), key=lambda a: a.pos.dist(n_in.pos), default=None)
    if ca is None:
        return None
    o = min(elements["O"], key=lambda a: a.pos.dist(c_out.pos), default=None)
    if o is not None and o.pos.dist(c_out.pos) > _CO_MAX:
        o = None
    return [n_in, ca, c_out, o]


def _chain_residues(chain):
    """(label, coords, is_hetero, CA B-factor, linked) per residue with N, CA and C; one residue per seqid
    (microheterogeneity). coords[3] (O) is None when the file has no carbonyl oxygen. A polymer residue of
    another kind (not a nucleotide) between two residues it is peptide-bonded to takes substitute backbone atoms
    (_substitute_backbone) and is `linked`."""
    seen: set[tuple[int, str]] = set()
    out = []
    residues = list(chain)
    for k, res in enumerate(residues):
        key = (res.seqid.num, res.seqid.icode.strip())
        atoms = [res.find_atom(a, "*") for a in _BB_ATOMS]
        linked = False
        if (
            all(a is None for a in atoms[:3])  # none of N, CA, C: a damaged residue is skipped, not stood in for
            and res.entity_type == gemmi.EntityType.Polymer
            and res.find_atom("C1'", "*") is None
        ):
            sub = _substitute_backbone(
                res, residues[k - 1] if k else None, residues[k + 1] if k + 1 < len(residues) else None
            )
            if sub is not None:
                atoms, linked = sub, True
        if key in seen or any(a is None for a in atoms[:3]):
            continue
        seen.add(key)
        label = ResLabel(chain.name, res.seqid.num, res.seqid.icode.strip(), res.name)
        coords = [a.pos.tolist() if a is not None else None for a in atoms]
        out.append((label, coords, res.het_flag == "H", atoms[1].b_iso, linked))
    return out


def _rebuild_oxygens(found):
    """Place a missing carbonyl O in the peptide plane, opposite the bisector of C->CA and C->N(next), 1.231 Å
    from C. A residue missing O with no peptide-bonded successor (a chain end or break) is dropped."""
    keep = []
    for k, item in enumerate(found):
        n, ca, c, o = (None if p is None else np.asarray(p, float) for p in item[1])
        if o is None:
            nxt = np.asarray(found[k + 1][1][0], float) if k + 1 < len(found) else None
            if nxt is None or np.linalg.norm(nxt - c) >= _PEPTIDE_BOND_MAX:
                continue
            u = (c - ca) / np.linalg.norm(c - ca) + (c - nxt) / np.linalg.norm(c - nxt)
            o = c + _CO_BOND * u / np.linalg.norm(u)
            item = (item[0], [*item[1][:3], o.tolist()], *item[2:])
        keep.append(item)
    return keep


def _drop_free_hetero(found):
    """Keep HETATM residues (MSE, SEP, ...) only when peptide-bonded to a neighbouring residue."""

    def bonded(a, b) -> bool:
        return np.linalg.norm(np.subtract(a[1][2], b[1][0])) < _PEPTIDE_BOND_MAX  # C(a) - N(b)

    keep = []
    for k, item in enumerate(found):
        if not item[2]:
            keep.append(item)
            continue
        before = k > 0 and bonded(found[k - 1], item)
        after = k + 1 < len(found) and bonded(item, found[k + 1])
        if before or after:
            keep.append(item)
    return keep


AUTO_MAX = 24  # `auto` builds assemblies up to this many chains (a capsid is better asked for by name)


def _copies(assembly, model) -> int:
    """Chains an assembly would build. PDB files name the chains of each generator; mmCIF files name subchains
    (label_asym_id), which are mapped back to the chains that hold them."""
    total = 0
    for g in assembly.generators:
        subchains = set(g.subchains)
        chains = set(g.chains) or {c.name for c in model if any(r.subchain in subchains for r in c)}
        total += len(g.operators) * len(chains)
    return total


COINCIDENT_FRAC = 0.90  # of a copy's atoms ...
COINCIDENT_DIST = 1.0  # ... within this many Å of one chain already present: the same molecule seen twice


def _drop_coincident(model, keep=frozenset(), tags=None) -> list[str]:
    """Remove, in place, every chain with COINCIDENT_FRAC of its atoms within COINCIDENT_DIST Å of one chain kept
    before it: an assembly operator that maps a molecule on a special position (DNA on a crystallographic
    two-fold) onto itself or its partner copies it exactly on top. `keep`: chains never dropped (the identity
    copy of each deposited chain), compared against first. `tags`: the operator that made each chain, in model
    order; chains of one operator are never compared (lying on each other, they were deposited so: alternative
    conformations). Returns the names dropped, in model order."""
    chains = list(model)
    tags = list(tags) if tags is not None else list(range(len(chains)))
    kept: list[tuple[np.ndarray, np.ndarray, cKDTree, object]] = []  # (low, high corner, tree, tag) per chain
    dropped = []
    order = sorted(range(len(chains)), key=lambda k: chains[k].name not in keep)  # stable: kept chains first
    for k in order:
        chain = chains[k]
        xyz = np.array([a.pos.tolist() for r in chain for a in r], float).reshape(-1, 3)
        if not len(xyz):
            continue
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        on_top = False
        for klo, khi, tree, tag in [] if chain.name in keep else kept:
            if tag == tags[k] or (lo > khi + COINCIDENT_DIST).any() or (hi < klo - COINCIDENT_DIST).any():
                continue  # one operator's copies, or boxes apart: no atom can be that close
            d, _ = tree.query(xyz, distance_upper_bound=COINCIDENT_DIST)
            if (d <= COINCIDENT_DIST).mean() >= COINCIDENT_FRAC:
                on_top = True
                break
        if on_top:
            dropped.append(chain.name)
        else:
            kept.append((lo, hi, cKDTree(xyz), tags[k]))
    dropped.sort(key=[c.name for c in chains].index)
    for name in dropped:
        model.remove_chain(name)
    return dropped


def _operator_tags(assembly, model, built) -> list[int] | None:
    """The operator (serial number over the generators) that made each chain of `built` =
    gemmi.make_assembly(assembly, model): it builds generator by generator, operator by operator, the model's
    chains in order. Checked atom by atom on each chain's first atom; None when the order is not that."""
    tags, expect = [], []
    serial = 0
    for g in assembly.generators:
        subchains, names = set(g.subchains), set(g.chains)
        members = []
        for c in model:
            first = next((r for r in c if r.subchain in subchains or c.name in names), None)
            if first is not None and len(first):
                members.append(first[0].pos)
        for op in g.operators:
            tags += [serial] * len(members)
            expect += [gemmi.Position(op.transform.apply(p)) for p in members]
            serial += 1
    chains = list(built)
    if len(tags) != len(chains):
        return None
    for c, p in zip(chains, expect):
        if not len(c) or not len(c[0]) or c[0][0].pos.dist(p) > 1e-3:
            return None
    return tags


def _polymer_entities(st, chains) -> list[str]:
    """The polymer entity of each chain's polymer, in order ('' where the file names none), entities of one
    protein (_same_proteins) under the first one's name."""
    one = _same_proteins(st)
    out = []
    for c in chains:
        polymer = c.get_polymer()
        if len(polymer):
            entity = st.get_entity_of(polymer)
            out.append(one.get(entity.name, entity.name) if entity is not None else "")
    return out


SAME_PROTEIN = 0.9  # two polymer entities are one protein when their full sequences align this identical (of the
# shorter) and the shorter is at least SAME_LENGTH of the longer (a bound peptide is not its receptor)
SAME_LENGTH = 0.8


def _same_proteins(st) -> dict[str, str]:
    """Polymer entity -> the first entity (file order) of its kind whose full sequence aligns with it at least
    SAME_PROTEIN identical over the shorter one, which is at least SAME_LENGTH as long as the longer (4UM3: a
    deletion and a point mutant of the protein deposited 38 times are entities of their own; a peptide bound to its
    receptor is not one with it). Entities without a full sequence stand alone."""
    seqs = [
        (e.name, e.polymer_type, [gemmi.Entity.first_mon(m) for m in e.full_sequence])
        for e in st.entities
        if e.entity_type == gemmi.EntityType.Polymer
    ]
    out: dict[str, str] = {}
    firsts: list[tuple[str, object, list[str]]] = []
    for name, kind, seq in seqs:
        out[name] = name
        for first, k, ref in firsts:
            if (
                k != kind
                or not seq
                or not ref
                or min(len(seq), len(ref)) < max(20, SAME_LENGTH * max(len(seq), len(ref)))
            ):
                continue
            same = gemmi.align_string_sequences(ref, seq, [], gemmi.AlignmentScoring()).match_count
            if same >= SAME_PROTEIN * min(len(seq), len(ref)):
                out[name] = first
                break
        else:
            firsts.append((name, kind, seq))
    return out


_POLYMER_KINDS = {
    gemmi.PolymerType.PeptideL: {"protein"},
    gemmi.PolymerType.PeptideD: {"protein"},
    gemmi.PolymerType.Dna: {"DNA"},
    gemmi.PolymerType.Rna: {"RNA"},
    gemmi.PolymerType.DnaRnaHybrid: {"DNA", "RNA"},
}


def _polymer_kinds(chains) -> set[str]:
    """The kinds of polymer the chains hold: protein, DNA, RNA. From each chain's polymer type, and from its
    nucleotides (_is_nucleotide; RNA where the sugar has an O2') for chains whose file gives no entity."""
    out = set()
    for c in chains:
        polymer = c.get_polymer()
        if len(polymer):
            out |= _POLYMER_KINDS.get(polymer.check_polymer_type(), set())
        for k in range(len(c)):
            if _is_nucleotide(c, k):
                out.add("RNA" if c[k].find_atom("O2'", "*") is not None else "DNA")
    return out


PAIRED_MIN = 3  # residue pairs of two chains' nucleotides with polar base atoms this near: held together (_paired)
_PAIR_DIST = 3.4  # angstrom
_BASE_POLAR = frozenset({"N1", "N3", "N2", "N4", "N6", "N7", "O2", "O4", "O6"})
# 3OWW's two RNAs pair in 4 residue pairs at _PAIR_DIST (A41-B84, A84-B41, A64-B7, A65-B65); 4PLX's three copies touch
# in the crystal by one swapped pair of bases (A28-B71, A71-B28: 2), which is no pairing (it drew the ASU). Chains
# of a protein only never count (their backbone contacts in a crystal or an assembly are no pairing: 7OKA's ten
# trimers touch in up to 10 residue pairs).


def _pair_atoms(chain) -> list[tuple[tuple[float, float, float], int]]:
    """(position, residue index) of the base N and O atoms of a chain's nucleotides (not amino acids)."""
    out = []
    for k, r in enumerate(chain):
        info = gemmi.find_tabulated_residue(r.name)
        if info is not None and info.is_amino_acid():
            continue
        out += [(a.pos.tolist(), k) for a in r if a.name in _BASE_POLAR]
    return out


def _paired(members, others) -> bool:
    """Are chains of `members` held to chains of `others` by base pairs between their nucleic acids (at least
    PAIRED_MIN residue pairs with polar base atoms within _PAIR_DIST)? Such chains are one object (3OWW: the two
    chains of a glycine riboswitch dimer pair their RNAs), however alike their sequences (3OWW ruling; protein
    β-strand pairing is not used)."""
    from scipy.spatial import cKDTree

    for a in members:
        pa = _pair_atoms(a)
        for b in others:
            pb = _pair_atoms(b)
            if not pa or not pb:
                continue
            tree = cKDTree([x for x, _ in pb])
            found = {(k, pb[j][1]) for x, k in pa for j in tree.query_ball_point(x, _PAIR_DIST)}
            if len(found) >= PAIRED_MIN:
                return True
    return False


def _one_complex_of_several(st, assembly, model) -> bool:
    """`assembly` is one of several copies of a complex deposited in the asymmetric unit (4OO8: two Cas9
    complexes, assembly 1 = A, B, C; 1URN: three U1A–RNA complexes; 4CH9: two Kelch domains with their WNK4
    peptides, assembly 1 = B, D; 4UM3: eight pentamers): every operator the identity, a strict subset of the
    deposited chains, each of its polymer entities (entities of one protein counted as one: _same_proteins)
    deposited at least twice as often as it holds it, and every kind of polymer the asymmetric
    unit holds (protein, DNA, RNA) among its own chains: a protein–DNA entry whose assembly 1 leaves the DNA out
    keeps the asymmetric unit, never silently losing the DNA. The other copies lie next to it by crystal packing,
    not by biology."""
    if not all(op.transform.is_identity() for g in assembly.generators for op in g.operators):
        return False
    subchains = {s for g in assembly.generators for s in g.subchains}
    names = {n for g in assembly.generators for n in g.chains}
    members = [c for c in model if c.name in names or any(r.subchain in subchains for r in c)]
    if not members or len(members) >= len(model):
        return False
    held, deposited = _polymer_entities(st, members), _polymer_entities(st, model)
    if not held or "" in held + deposited:
        return False
    if _polymer_kinds(members) != _polymer_kinds(model):
        return False
    if not all(deposited.count(e) >= 2 * held.count(e) for e in set(held)):
        return False
    return not _paired(members, [c for c in model if c not in members])  # copies paired to it are one object


class ModelRead(tuple):
    """read_model's result: (model, deposited chain names), unpacked as a pair. `several`: when 'auto' drew assembly
    1 as one of several copies of a complex the asymmetric unit holds, what the asymmetric unit holds (assembly
    name, copies, polymer chains, polymer residues); else None."""

    several: dict | None

    def __new__(cls, model, names, several: dict | None = None):
        out = super().__new__(cls, (model, names))
        out.several = several
        return out


def _polymers(model) -> tuple[int, int]:
    """(chains, residues) of the model's polymers: protein and nucleic acid."""
    sizes = [len(c.get_polymer()) for c in model]
    return sum(n > 0 for n in sizes), sum(sizes)


def read_model(path: str | Path, assembly: str = "auto"):
    """(model, deposited chain names). assembly: 'asu' as deposited; an assembly id from the file; or 'auto':
    the first assembly when it adds copies and stays within AUTO_MAX chains, or when it is one of several copies
    of a complex the asymmetric unit holds (_one_complex_of_several; drawn under the deposited chain
    names, and only its copies are looked at for symmetry). A copy the assembly lays exactly on a chain already
    present is dropped (_drop_coincident); under 'auto', when that drops every copy it adds, the model is returned
    as deposited, with its deposited chain names. A ModelRead: a (model, names) pair that also says when 'auto' left
    copies of the complex out (`several`)."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        raise ValueError(f"{path}: no models")
    model = st[0]
    deposited = {c.name for c in model}
    if assembly == "asu" or not st.assemblies:
        if assembly not in ("asu", "auto"):
            raise ValueError(f"{path}: no assembly {assembly!r} in the file (it has none; use asu)")
        return ModelRead(model, deposited)
    if assembly == "auto":
        chosen = st.assemblies[0]
        copies = _copies(chosen, model)
        if copies < len(deposited) and _one_complex_of_several(st, chosen, model):
            one = gemmi.make_assembly(chosen, model, gemmi.HowToNameCopiedChain.Dup)  # identity: names kept
            (chains, residues), held = _polymers(model), _polymers(one)[0]
            several = {"assembly": chosen.name, "copies": round(chains / max(held, 1)), "chains": chains}
            return ModelRead(one, {c.name for c in one}, {**several, "residues": residues})
        if copies <= len(deposited) or copies > AUTO_MAX:
            return ModelRead(model, deposited)
    else:
        found = [a for a in st.assemblies if a.name == assembly]
        if not found:
            raise ValueError(
                f"{path}: no assembly {assembly!r}; the file has {', '.join(a.name for a in st.assemblies)}"
            )
        chosen = found[0]
    built = gemmi.make_assembly(chosen, model, gemmi.HowToNameCopiedChain.AddNumber)
    tags = _operator_tags(chosen, model, built)
    original = [(c.name, np.array([a.pos.tolist() for r in c for a in r][:20])) for c in model]
    asu, claimed = set(), set()
    for c in built:  # the copy made by the identity operator sits exactly on a deposited chain (one per chain)
        xyz = np.array([a.pos.tolist() for r in c for a in r][:20])
        for name, o in original:
            if name not in claimed and len(o) == len(xyz) and np.allclose(o, xyz, atol=1e-3):
                asu.add(c.name)
                claimed.add(name)
                break
    _drop_coincident(built, keep=asu, tags=tags)
    if assembly == "auto" and len(asu) == len(built) == len(model):  # every added copy dropped: as deposited
        return ModelRead(model, deposited)
    return ModelRead(built, asu)


def load_backbone(
    path: str | Path, assembly: str = "auto", *, model: tuple[gemmi.Model, set[str]] | None = None
) -> Backbone:
    """`model`: the result of read_model(path, assembly), with this same assembly, already in hand (a figure reads
    the file once for every loader)."""
    model, asu = model or read_model(path, assembly)
    labels: list[ResLabel] = []
    coords: list[list[list[float]]] = []
    bfac: list[float] = []
    linked: set[int] = set()
    for chain in model:
        found = _chain_residues(chain)
        found = _rebuild_oxygens(_drop_free_hetero(found))
        linked.update(len(labels) + k for k, f in enumerate(found) if f[4])
        labels.extend(f[0] for f in found)
        coords.extend(f[1] for f in found)
        bfac.extend(f[3] for f in found)
    if not labels:
        if any(
            gemmi.find_tabulated_residue(r.name) is not None
            and gemmi.find_tabulated_residue(r.name).is_amino_acid()
            and r.find_atom("CA", "*") is not None
            and r.find_atom("N", "*") is None
            for c in model
            for r in c
        ):
            raise ValueError(
                f"{path}: only Cα atoms; foldmap needs the N, CA and C atoms of each residue to assign secondary "
                "structure (rebuild the backbone first, e.g. with PULCHRA or modeller)"
            )
        if any(_is_nucleotide(chain, k) for chain in model for k in range(len(chain))):  # nucleic acid only:
            return Backbone([], np.zeros((0, 4, 3)), np.zeros(0, bool), set(), np.zeros(0))  # drawn on its own
        raise ValueError(f"{path}: no protein residues with a complete backbone")
    xyz = np.asarray(coords, float)
    prev = np.zeros(len(labels), bool)
    for i in range(1, len(labels)):
        if labels[i].chain == labels[i - 1].chain:
            prev[i] = np.linalg.norm(xyz[i - 1, 2] - xyz[i, 0]) < _PEPTIDE_BOND_MAX
    chains = {l.chain for l in labels}
    return Backbone(labels, xyz, prev, (asu & chains) or chains, np.asarray(bfac, float), frozenset(linked))


_NA_LINK_MAX = 2.0  # Å, O3'(i-1)-P(i)
CONTACT = 4.0  # Å, heavy atoms of a protein residue and a nucleotide


def _linked(a, b) -> bool:
    """The O3' of nucleotide a is bonded to the P of b."""
    o3, p = a.find_atom("O3'", "*"), b.find_atom("P", "*")
    return o3 is not None and p is not None and o3.pos.dist(p.pos) <= _NA_LINK_MAX


def _is_nucleotide(chain, k: int) -> bool:
    """chain[k] is a nucleotide: tabulated as one (modified ones included), or, for a code gemmi does not know
    (5MC, YYG), a sugar with its C1', C4' and O4' plus a phosphate or O3' that is part of a polymer entity or,
    where the file sets no entity types, O3'-P linked to a neighbouring nucleotide. So a nucleotide ligand (ATP,
    GTP, SAM) is never a strand."""
    res = chain[k]
    if res.find_atom("C1'", "*") is None:
        return False
    info = gemmi.find_tabulated_residue(res.name)
    if info is not None and info.is_nucleic_acid():
        return True
    names = {a.name for a in res}
    if not ({"C4'", "O4'"} <= names and ("P" in names or "O3'" in names)):
        return False
    if res.entity_type != gemmi.EntityType.Unknown:
        return res.entity_type == gemmi.EntityType.Polymer
    before = k > 0 and chain[k - 1].find_atom("C1'", "*") is not None and _linked(chain[k - 1], res)
    after = k + 1 < len(chain) and chain[k + 1].find_atom("C1'", "*") is not None and _linked(res, chain[k + 1])
    return before or after


def _nucleic_strands(model):
    """[(NAStrand, [gemmi residues])], split where the O3'-P link is missing."""
    out = []
    for chain in model:
        run: list = []
        for k, res in enumerate(chain):
            if not _is_nucleotide(chain, k):
                continue
            if run and not _linked(run[-1], res):
                out.append((chain.name, run))
                run = []
            run.append(res)
        if run:
            out.append((chain.name, run))
    strands = []
    for name, run in out:

        def point(r):
            a = r.find_atom("P", "*") or r.find_atom("C4'", "*") or r.find_atom("C1'", "*")
            return a.pos.tolist()

        labels = [ResLabel(name, r.seqid.num, r.seqid.icode.strip(), r.name) for r in run]
        rna = any(r.find_atom("O2'", "*") is not None for r in run)
        strands.append(
            (
                NAStrand(
                    name,
                    labels,
                    np.array([point(r) for r in run]),
                    np.array([r.find_atom("C1'", "*").pos.tolist() for r in run]),
                    rna,
                    bases="".join(parent_base(r) or "N" for r in run),
                ),
                run,
            )
        )
    return strands


def load_nucleic(
    path: str | Path, bb: Backbone, assembly: str = "auto", *, model: tuple[gemmi.Model, set[str]] | None = None
) -> Nucleic:
    """Nucleic-acid strands of model 1, their base pairs, and which protein residues touch which nucleotides.
    `model`: as for load_backbone (same assembly)."""
    model, _ = model or read_model(path, assembly)
    st = gemmi.Structure()  # NeighborSearch wants a cell; an empty one means none

    found = _nucleic_strands(model)
    if not found:
        return Nucleic()
    ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    where = {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}
    contacts: dict[int, set[tuple[int, int]]] = {}
    for s, (_, run) in enumerate(found):
        for i, r in enumerate(run):
            for atom in r:
                if atom.element.name == "H":
                    continue
                for mark in ns.find_atoms(atom.pos, "\0", radius=CONTACT):
                    cra = mark.to_cra(model)
                    if cra.atom.element.name == "H":
                        continue
                    k = where.get((cra.chain.name, cra.residue.seqid.num, cra.residue.seqid.icode.strip()))
                    if k is not None:
                        contacts.setdefault(k, set()).add((s, i))
    strands = [s for s, _ in found]
    runs = [run for _, run in found]
    chains = {k: bb.labels[k].chain for k in contacts}
    return nucleic_from_pairs(strands, pair_nucleotides(runs), contacts, chains, _breaks(found), g_tetrads(runs))


def _breaks(found) -> list[NAGap]:
    """Where a chain is broken inside itself (unmodelled nucleotides, or an O3'-P link too long): a gap from one
    strand piece's 3' end to the next one's 5' end, counting the nucleotides left out by sequence position
    (label_seq), else by residue number (4P9R is numbered 22 -> 216 with one nucleotide missing). The pieces' ends
    there are marked cut: the molecule goes on beyond them."""
    gaps = []
    for k, ((a, ra), (b, rb)) in enumerate(zip(found, found[1:])):
        if a.chain != b.chain:
            continue
        x, y = ra[-1], rb[0]
        if x.label_seq is not None and y.label_seq is not None:
            count = y.label_seq - x.label_seq - 1
        else:
            count = y.seqid.num - x.seqid.num - 1
        gaps.append(NAGap(k, k + 1, max(count, 0)))
        a.cut3 = b.cut5 = True
    return gaps


SUGARS = {
    "NAG",
    "NDG",
    "A2G",
    "NGA",
    "MAN",
    "BMA",
    "GAL",
    "GLA",
    "GLC",
    "BGC",
    "FUC",
    "FUL",
    "SIA",
    "XYS",
    "XYP",
    "GCS",
    "GLB",
    "SLB",
    "NGC",
    "KDN",
}
ANOMERIC = {  # sugars whose glycosidic bond leaves from C2 (the sialic acids: Neu5Ac, Neu5Gc, Kdn); others: C1
    "SIA": "C2",
    "SLB": "C2",
    "NGC": "C2",
    "KDN": "C2",
}


def anomeric(res):
    """The sugar's anomeric carbon: the atom of its glycosidic bond to the protein or the next sugar in."""
    return res.find_atom(ANOMERIC.get(res.name, "C1"), "*")


_SS_MAX = 2.5  # Å between cysteine SG atoms
_GLYCAN_LINK_MAX = 2.0  # Å between a sugar's anomeric carbon and the atom it is bonded to


def load_links(
    path: str | Path, bb: Backbone, assembly: str = "auto", *, model: tuple[gemmi.Model, set[str]] | None = None
) -> Links:
    """Disulfide bridges (SG-SG) and glycans (a sugar's anomeric carbon on Asn ND2, Ser OG or Thr OG1, then
    sugar-to-sugar). Sugars in no glycan (bound saccharides, receptor analogues) are ligands.
    `model`: as for load_backbone (same assembly)."""
    model, _ = model or read_model(path, assembly)
    where = {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}

    def key(chain, res):
        return (chain.name, res.seqid.num, res.seqid.icode.strip())

    sg = [(key(c, r), r.find_atom("SG", "*")) for c in model for r in c if r.name == "CYS" and r.find_atom("SG", "*")]
    disulfides = []
    for a in range(len(sg)):
        for b in range(a + 1, len(sg)):
            if sg[a][1].pos.dist(sg[b][1].pos) < _SS_MAX and sg[a][0] in where and sg[b][0] in where:
                disulfides.append(tuple(sorted((where[sg[a][0]], where[sg[b][0]]))))
    sugars = [(c, r) for c in model for r in c if r.name in SUGARS and anomeric(r) is not None]
    anchors = [
        (key(c, r), r.find_atom(n, "*"))
        for c in model
        for r in c
        for n in {"ASN": ("ND2",), "SER": ("OG",), "THR": ("OG1",)}.get(r.name, ())
        if r.find_atom(n, "*")
    ]
    glycans, members = [], []
    used = set()
    for k_res, atom in anchors:  # every tree is grown, so a left-out chain's glycan is not taken for a ligand
        root = next(
            ((c, r) for c, r in sugars if id(r) not in used and anomeric(r).pos.dist(atom.pos) < _GLYCAN_LINK_MAX),
            None,
        )
        if root is None:
            continue
        tree, frontier = [root], [root]
        used.add(id(root[1]))
        while frontier:  # sugars whose anomeric carbon bonds to an oxygen of a sugar already in the tree
            nxt = []
            for _, parent in frontier:
                oxygens = [a for a in parent if a.element.name == "O"]
                for c, r in sugars:
                    if id(r) in used:
                        continue
                    link = anomeric(r)
                    if any(link.pos.dist(o.pos) < _GLYCAN_LINK_MAX for o in oxygens):
                        used.add(id(r))
                        tree.append((c, r))
                        nxt.append((c, r))
            frontier = nxt
        if k_res in where:
            glycans.append((where[k_res], [r.name for _, r in tree]))
            members.append([(c.name, r.seqid.num) for c, r in tree])
    in_glycans = {key(c, r) for c, r in sugars if id(r) in used}
    return Links(sorted(set(disulfides)), glycans, _ligands(model, where, in_glycans), members)


ADDITIVES = {  # crystallisation and purification leftovers: hidden from figures unless asked for
    "SO4",
    "PO4",
    "GOL",
    "EDO",
    "PEG",
    "PG4",
    "PGE",
    "1PE",
    "P6G",
    "ACT",
    "ACY",
    "FMT",
    "DMS",
    "MPD",
    "TRS",
    "EPE",
    "MES",
    "BME",
    "IMD",
    "CIT",
    "FLC",
    "TAR",
    "MLI",
    "NO3",
    "SCN",
    "AZI",
    "IPA",
    "EOH",
    "MOH",
    "BU3",
    "CL",
    "BR",
    "IOD",
    "NA",
    "UNX",
    "UNL",
    "BOG",
    "LDA",
    "LMT",
    "DDM",
    "OLC",
    "OLA",
    "PLM",
    "LI1",
    "SQU",
    "CPS",
    "HTG",
    "C8E",
    "UMQ",
    "PC1",
    "POV",
    "CDL",
    "LHG",
}
_METAL_REACH = 3.1  # Å, metal to a coordinating atom (K+ to carbonyl O is ~2.8-3.0)
_LIGAND_REACH = 4.0  # Å, ligand heavy atom to a residue atom


SHORT_NAMES = {  # sugar residues by their usual short names (SNFG); other ligands keep their CCD code
    "NAG": "GlcNAc",
    "NDG": "GlcNAc",
    "GN1": "GlcNAc-P",
    "A1E9H": "GlcNAc-P",
    "RAM": "Rha",
    "GZL": "Galf",
    "GAL": "Gal",
    "GLA": "Gal",
    "GLC": "Glc",
    "BGC": "Glc",
    "MAN": "Man",
    "BMA": "Man",
    "FUC": "Fuc",
    "XYL": "Xyl",
    "SIA": "Neu5Ac",
    "SLB": "Neu5Ac",
    "NGC": "Neu5Gc",
    "KDN": "Kdn",
    "A2G": "GalNAc",
    "NGA": "GalNAc",
}
_BOND = 1.9  # Å between heavy atoms of two residues that are covalently joined
_OPM_DUMMY = "DUM"  # OPM/PPM pseudo-atoms marking the membrane boundaries


def _ligand_label(codes: list[str]) -> str:
    """A1E8P·GlcNAc-P·Rha·Galf×3: residues in bond order, short sugar names, runs compressed."""
    names = [SHORT_NAMES.get(c, c) for c in codes]
    out: list[str] = []
    for n in names:
        if out and out[-1].split("×")[0] == n:
            k = int(out[-1].split("×")[1]) if "×" in out[-1] else 1
            out[-1] = f"{n}×{k + 1}"
        else:
            out.append(n)
    return "·".join(out)


def _ligands(model, where, skip: set | None = None) -> list[Ligand]:
    """Ligands touching the protein, with the residues that hold them. Non-polymer residues bonded to each other
    (a lipid-linked oligosaccharide, a cofactor built of parts) are one ligand; metals stay single ions. skip:
    (chain, number, icode) of residues drawn otherwise (the sugars of glycans)."""
    skip = skip or set()
    search = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate()
    parts = _ligand_parts(model, skip)
    group = list(range(len(parts)))

    def root(x: int) -> int:
        while group[x] != x:
            group[x] = group[group[x]]
            x = group[x]
        return x

    bonded = _bonded_parts(parts)
    for a, b in sorted(bonded):  # the order the pairwise loop met them in
        group[root(a)] = root(b)
    members: dict[int, list[int]] = {}
    for k in range(len(parts)):
        members.setdefault(root(k), []).append(k)

    out = []
    for ks in members.values():
        ks = _bond_order(ks, bonded)
        touching = set()
        for k in ks:
            _, res, heavy, metal = parts[k]
            reach = _METAL_REACH if metal else _LIGAND_REACH
            for atom in heavy:
                for mark in search.find_atoms(atom.pos, radius=reach):
                    cra = mark.to_cra(model)
                    key = (cra.chain.name, cra.residue.seqid.num, cra.residue.seqid.icode.strip())
                    if key in where and (not metal or cra.atom.element.name in ("O", "N", "S")):
                        touching.add(where[key])
        if not touching:
            continue
        chain, res, heavy, metal = parts[ks[0]]
        codes = [parts[k][1].name for k in ks]
        symbol = heavy[0].element.name.capitalize() if metal else _ligand_label(codes)
        name = codes[0] if len(codes) == 1 else "+".join(dict.fromkeys(codes))
        out.append(
            Ligand(
                name,
                chain,
                res.seqid.num,
                metal,
                symbol,
                sorted(touching),
                all(c in ADDITIVES for c in codes),
                [(parts[k][0], parts[k][1].seqid.num) for k in ks],
                codes,
            )
        )
    return out


def _ligand_parts(model, skip: set | None = None) -> list[tuple[str, gemmi.Residue, list[gemmi.Atom], bool]]:
    """(chain name, residue, heavy atoms, is a metal ion) for every residue outside the polymers, waters,
    skipped glycan sugars and OPM's DUM membrane-boundary markers (read by membrane._opm, not molecules): the
    pieces _ligands joins into ligands."""
    skip = skip or set()
    parts = []
    for chain in model:
        polymer = {(r.seqid.num, r.seqid.icode) for r in chain.get_polymer()}
        for res in chain:
            if (
                (res.seqid.num, res.seqid.icode) in polymer
                or res.is_water()
                or res.name == _OPM_DUMMY
                or (chain.name, res.seqid.num, res.seqid.icode.strip()) in skip
            ):
                continue
            heavy = [a for a in res if a.element.name != "H"]
            if heavy:
                parts.append((chain.name, res, heavy, len(heavy) == 1 and heavy[0].element.is_metal))
    return parts


def _bonded_parts(parts) -> set[tuple[int, int]]:
    """(a, b), a < b: ligand parts with heavy atoms closer than _BOND. Metal ions bond to nothing."""
    atoms = [(k, a) for k, (_, _, heavy, metal) in enumerate(parts) if not metal for a in heavy]
    if not atoms:
        return set()
    xyz = np.array([a.pos.tolist() for _, a in atoms])
    out = set()
    for i, j in sorted(cKDTree(xyz).query_pairs(_BOND)):
        (a, x), (b, y) = atoms[i], atoms[j]
        if a != b and x.pos.dist(y.pos) < _BOND:  # query_pairs includes the boundary; the cutoff does not
            out.add((min(a, b), max(a, b)))
    return out


def _bond_order(ks: list[int], bonded: set[tuple[int, int]]) -> list[int]:
    """The residues of one ligand from an end (a residue with one bonded partner), walking along the bonds."""
    if len(ks) == 1:
        return ks
    near = {k: [m for m in ks if m != k and (min(k, m), max(k, m)) in bonded] for k in ks}
    start = min((k for k in ks if len(near[k]) <= 1), default=ks[0], key=lambda k: ks.index(k))
    order, seen = [start], {start}
    while len(order) < len(ks):
        nxt = [m for m in near[order[-1]] if m not in seen] or [m for m in ks if m not in seen]
        order.append(nxt[0])
        seen.add(nxt[0])
    return order
