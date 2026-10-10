"""Read PDB/mmCIF into a Backbone (first model, first altloc, amino acids only)."""

from __future__ import annotations

from pathlib import Path

import gemmi
import numpy as np

from .model import Backbone, Ligand, Links, NAStrand, Nucleic, ResLabel

_BB_ATOMS = ("N", "CA", "C", "O")
_PEPTIDE_BOND_MAX = 2.5  # Å, C(i-1)-N(i), as in DSSP
_CO_BOND = 1.231  # Å, carbonyl C=O


def _chain_residues(chain):
    """(label, coords, is_hetero, CA B-factor) per residue with N, CA and C; one residue per seqid
    (microheterogeneity). coords[3] (O) is None when the file has no carbonyl oxygen."""
    seen: set[tuple[int, str]] = set()
    out = []
    for res in chain:
        key = (res.seqid.num, res.seqid.icode.strip())
        atoms = [res.find_atom(a, "*") for a in _BB_ATOMS]
        if key in seen or any(a is None for a in atoms[:3]):
            continue
        seen.add(key)
        label = ResLabel(chain.name, res.seqid.num, res.seqid.icode.strip(), res.name)
        coords = [a.pos.tolist() if a is not None else None for a in atoms]
        out.append((label, coords, res.het_flag == "H", atoms[1].b_iso))
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


def read_model(path: str | Path, assembly: str = "auto"):
    """(model, deposited chain names). assembly: 'asu' as deposited; an assembly id from the file; or 'auto':
    the first assembly when it adds copies and stays within AUTO_MAX chains."""
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
        return model, deposited
    if assembly == "auto":
        chosen = st.assemblies[0]
        copies = _copies(chosen, model)
        if copies <= len(deposited) or copies > AUTO_MAX:
            return model, deposited
    else:
        found = [a for a in st.assemblies if a.name == assembly]
        if not found:
            raise ValueError(
                f"{path}: no assembly {assembly!r}; the file has {', '.join(a.name for a in st.assemblies)}"
            )
        chosen = found[0]
    built = gemmi.make_assembly(chosen, model, gemmi.HowToNameCopiedChain.AddNumber)
    original = {c.name: np.array([a.pos.tolist() for r in c for a in r][:20]) for c in model}
    asu = set()
    for c in built:  # the copy made by the identity operator sits exactly on a deposited chain
        xyz = np.array([a.pos.tolist() for r in c for a in r][:20])
        if any(len(o) == len(xyz) and np.allclose(o, xyz, atol=1e-3) for o in original.values()):
            asu.add(c.name)
    return built, asu


def load_backbone(path: str | Path, assembly: str = "auto") -> Backbone:
    model, asu = read_model(path, assembly)
    labels: list[ResLabel] = []
    coords: list[list[list[float]]] = []
    bfac: list[float] = []
    for chain in model:
        found = _chain_residues(chain)
        found = _rebuild_oxygens(_drop_free_hetero(found))
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
        raise ValueError(f"{path}: no protein residues with a complete backbone")
    xyz = np.asarray(coords, float)
    prev = np.zeros(len(labels), bool)
    for i in range(1, len(labels)):
        if labels[i].chain == labels[i - 1].chain:
            prev[i] = np.linalg.norm(xyz[i - 1, 2] - xyz[i, 0]) < _PEPTIDE_BOND_MAX
    chains = {l.chain for l in labels}
    return Backbone(labels, xyz, prev, (asu & chains) or chains, np.asarray(bfac, float))


_PURINES = {"A", "G", "DA", "DG", "I", "DI"}
_PYRIMIDINES = {"C", "U", "T", "DC", "DT", "DU"}
_NA_LINK_MAX = 2.0  # Å, O3'(i-1)-P(i)
_WC_MAX = 3.4  # Å, purine N1 - pyrimidine N3 in a Watson-Crick pair
_CONTACT = 4.0  # Å, heavy atoms of a protein residue and a nucleotide


def _is_nucleotide(res) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_nucleic_acid()) and res.find_atom("C1'", "*") is not None


def _nucleic_strands(model):
    """[(NAStrand, [gemmi residues])], split where the O3'-P link is missing."""
    out = []
    for chain in model:
        run: list = []
        for res in chain:
            if not _is_nucleotide(res):
                continue
            if run:
                o3, p = run[-1].find_atom("O3'", "*"), res.find_atom("P", "*")
                if o3 is None or p is None or o3.pos.dist(p.pos) > _NA_LINK_MAX:
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
                ),
                run,
            )
        )
    return strands


def _watson_crick(strands) -> list[tuple[int, int, int, int]]:
    """Each nucleotide's closest partner with purine N1 - pyrimidine N3 under _WC_MAX, on another strand or
    far along the same one (a hairpin), mutual best only."""
    atoms = []  # (strand, index, is_purine, position)
    for a, (_, run) in enumerate(strands):
        for i, r in enumerate(run):
            name = r.name.strip()
            pur = name in _PURINES
            if not pur and name not in _PYRIMIDINES:
                continue
            at = r.find_atom("N1" if pur else "N3", "*")
            if at is not None:
                atoms.append((a, i, pur, np.array(at.pos.tolist())))
    best: dict[tuple[int, int], tuple[float, tuple[int, int]]] = {}
    for x in atoms:
        for y in atoms:
            if x[2] == y[2] or (x[0] == y[0] and abs(x[1] - y[1]) < 3):
                continue
            d = float(np.linalg.norm(x[3] - y[3]))
            if d < _WC_MAX and d < best.get((x[0], x[1]), (np.inf,))[0]:
                best[(x[0], x[1])] = (d, (y[0], y[1]))
    pairs = set()
    for k, (_, other) in best.items():
        if best.get(other, (0, None))[1] == k:
            pairs.add(min((k[0], k[1], other[0], other[1]), (other[0], other[1], k[0], k[1])))
    return sorted(pairs)


def load_nucleic(path: str | Path, bb: Backbone, assembly: str = "auto") -> Nucleic:
    """Nucleic-acid strands of model 1, their base pairs, and which protein residues touch which nucleotides."""
    model, _ = read_model(path, assembly)
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
                for mark in ns.find_atoms(atom.pos, "\0", radius=_CONTACT):
                    cra = mark.to_cra(model)
                    if cra.atom.element.name == "H":
                        continue
                    k = where.get((cra.chain.name, cra.residue.seqid.num, cra.residue.seqid.icode.strip()))
                    if k is not None:
                        contacts.setdefault(k, set()).add((s, i))
    return Nucleic([s for s, _ in found], _watson_crick(found), contacts, {k: bb.labels[k].chain for k in contacts})


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
}
ANOMERIC = {  # sugars whose glycosidic bond leaves from C2 (sialic acids, Kdn, Kdo, fructose); every other one: C1
    "SIA": "C2",
    "SLB": "C2",
    "NGC": "C2",
    "KDN": "C2",
    "KDO": "C2",
    "FRU": "C2",
}


def anomeric(res):
    """The sugar's anomeric carbon: the atom of its glycosidic bond to the protein or the next sugar in."""
    return res.find_atom(ANOMERIC.get(res.name, "C1"), "*")


_SS_MAX = 2.5  # Å between cysteine SG atoms
_GLYCAN_LINK_MAX = 2.0  # Å between a sugar's anomeric carbon and the atom it is bonded to


def load_links(path: str | Path, bb: Backbone, assembly: str = "auto") -> Links:
    """Disulfide bridges (SG-SG) and glycans (a sugar's anomeric carbon on Asn ND2, Ser OG or Thr OG1, then
    sugar-to-sugar). Sugars in no glycan (bound saccharides, receptor analogues) are ligands."""
    model, _ = read_model(path, assembly)
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
    "A2G": "GalNAc",
    "NGA": "GalNAc",
}
_BOND = 1.9  # Å between heavy atoms of two residues that are covalently joined


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
    parts = []
    for chain in model:
        polymer = {(r.seqid.num, r.seqid.icode) for r in chain.get_polymer()}
        for res in chain:
            if (
                (res.seqid.num, res.seqid.icode) in polymer
                or res.is_water()
                or (chain.name, res.seqid.num, res.seqid.icode.strip()) in skip
            ):
                continue
            heavy = [a for a in res if a.element.name != "H"]
            if heavy:
                parts.append((chain.name, res, heavy, len(heavy) == 1 and heavy[0].element.is_metal))
    group = list(range(len(parts)))

    def root(x: int) -> int:
        while group[x] != x:
            group[x] = group[group[x]]
            x = group[x]
        return x

    for a in range(len(parts)):
        for b in range(a + 1, len(parts)):
            if parts[a][3] or parts[b][3]:
                continue
            if any(x.pos.dist(y.pos) < _BOND for x in parts[a][2] for y in parts[b][2]):
                group[root(a)] = root(b)
    members: dict[int, list[int]] = {}
    for k in range(len(parts)):
        members.setdefault(root(k), []).append(k)

    out = []
    for ks in members.values():
        ks = _bond_order(ks, parts)
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


def _bond_order(ks: list[int], parts) -> list[int]:
    """The residues of one ligand from an end (a residue with one bonded partner), walking along the bonds."""
    if len(ks) == 1:
        return ks
    near = {
        k: [m for m in ks if m != k and any(x.pos.dist(y.pos) < _BOND for x in parts[k][2] for y in parts[m][2])]
        for k in ks
    }
    start = min((k for k in ks if len(near[k]) <= 1), default=ks[0], key=lambda k: ks.index(k))
    order, seen = [start], {start}
    while len(order) < len(ks):
        nxt = [m for m in near[order[-1]] if m not in seen] or [m for m in ks if m not in seen]
        order.append(nxt[0])
        seen.add(nxt[0])
    return order
