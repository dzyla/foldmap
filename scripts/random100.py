# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The random-100 campaign: a stratified random sample of 100 PDB entries and AlphaFold models, each drawn in every
view foldmap offers, measured and shown on a contact sheet (plan 10, needs network; not run in CI).

    python scripts/random100.py sample [--seed 20261013] [--date YYYY-MM-DD] [--out .../sample.json] [--check]

sample draws the sample from the RCSB search API (https://search.rcsb.org/rcsbsearch/v2/query) and writes
docs/validation/random100/sample.json. Every pool is one stratum x size-bin query, all hits, sorted; it is drawn
from in the order random.Random(f"{seed}:{stratum}:{bin}") gives it (a string seed is hashed with SHA-512, so
PYTHONHASHSEED cannot change the draw), at most LOOK candidates deep. A candidate is downloaded (foldmap.fetch, into
~/.cache/foldmap) and read with gemmi; it is skipped, with the reason logged, when it is obsolete or cannot be
downloaded, holds no protein or nucleic-acid polymer, is CA-only or P-only, or is in the sample already. The extra
quotas inside a stratum (S4 helical, big, cap; S8 ribosome) are drawn first from their own sub-pools; the pairs of S11
follow PAIR_TYPES; the cross-cutting quotas (QUOTAS) are then topped up. The same seed, date and pools give the same
file byte for byte. --check re-queries the ids and lists those now obsolete, without changing the file.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import random
import subprocess
import sys
import urllib.error
import urllib.request
from collections.abc import Callable, Mapping
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    # the foldmap of this checkout, not whichever one is installed (children run with python -I: no PYTHONPATH)
    sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

SEED = 20261013
LOOK = 200  # candidates looked at per pool, at most
OUT = ROOT / "docs" / "validation" / "random100"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

# --------------------------------------------------------------------------------------------- strata


def term(attribute: str, operator: str, value=None, *, negation: bool = False) -> dict:
    """One terminal node of an RCSB search query (text service)."""
    params: dict = {"attribute": attribute, "operator": operator}
    if value is not None:
        params["value"] = value
    if negation:
        params["negation"] = True
    return {"type": "terminal", "service": "text", "parameters": params}


def either(*nodes: dict) -> dict:
    return {"type": "group", "logical_operator": "or", "nodes": list(nodes)}


def size_range(lo: int | None, hi: int) -> dict:
    """deposited polymer residues in [lo, hi)."""
    attr = "rcsb_entry_info.deposited_polymer_monomer_count"
    if lo is None:
        return term(attr, "less", hi)
    return term(attr, "range", {"from": lo, "to": hi, "include_lower": True, "include_upper": False})


# size bins by deposited polymer residues, [lo, hi); XL-cap (20,000-30,000) is a sub-bin of XL
BINS = {"S": (None, 300), "M": (300, 1500), "L": (1500, 8000), "XL": (8000, 30001)}
CAP = (20000, 30001)

METHOD = "rcsb_entry_info.experimental_method"
TYPES = "rcsb_entry_info.selected_polymer_entity_types"
XRAY = term(METHOD, "exact_match", "X-ray")
EM = term(METHOD, "exact_match", "EM")
NMR = term(METHOD, "exact_match", "NMR")
PROTEIN_ONLY = term(TYPES, "exact_match", "Protein (only)")
PROTEIN = term("rcsb_entry_info.polymer_entity_count_protein", "greater", 0)
DNA = term("rcsb_entry_info.polymer_entity_count_DNA", "greater", 0)
NO_DNA = term("rcsb_entry_info.polymer_entity_count_DNA", "equals", 0)
RNA = term("rcsb_entry_info.polymer_entity_count_RNA", "greater", 0)
NO_RNA = term("rcsb_entry_info.polymer_entity_count_RNA", "equals", 0)
ONE_PROTEIN = term("rcsb_entry_info.polymer_entity_count_protein", "equals", 1)
UNIPROT = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers"

# The table of plan 10's "Sample design": n entries per stratum, per size bin; extra quotas drawn first from the
# stratum's pool narrowed by their nodes (their bin is the bin pool that holds the id; `bin` limits the sub-pool).
STRATA: list[dict] = [
    {
        "name": "S1",
        "title": "X-ray, 1 protein chain",
        "n": 12,
        "bins": {"S": 6, "M": 6},
        "nodes": [XRAY, PROTEIN_ONLY, term("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1)],
    },
    {
        "name": "S2",
        "title": "X-ray homo-oligomer (symmetric)",
        "n": 12,
        "bins": {"S": 2, "M": 4, "L": 4, "XL": 2},
        "nodes": [
            XRAY,
            PROTEIN_ONLY,
            ONE_PROTEIN,
            term("rcsb_assembly_info.polymer_entity_instance_count", "greater_or_equal", 2),
            term("rcsb_struct_symmetry.symbol", "exact_match", "C1", negation=True),
        ],
    },
    {
        "name": "S3",
        "title": "X-ray hetero-complex",
        "n": 10,
        "bins": {"S": 1, "M": 3, "L": 4, "XL": 2},
        "nodes": [XRAY, PROTEIN_ONLY, term("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 2)],
    },
    {
        "name": "S4",
        "title": "cryo-EM",
        "n": 15,
        "bins": {"M": 2, "L": 5, "XL": 8},
        "nodes": [EM, PROTEIN_ONLY],
        "extra": [
            {"tag": "helical", "n": 1, "nodes": [term("rcsb_struct_symmetry.symbol", "exact_match", "H")]},
            {"tag": "big", "n": 2, "nodes": [term("rcsb_assembly_info.polymer_entity_instance_count", "greater", 20)]},
            {"tag": "cap", "n": 3, "nodes": [size_range(*CAP)], "bin": "XL"},
        ],
    },
    {
        "name": "S5",
        "title": "NMR, multi-model",
        "n": 8,
        "bins": {"S": 7, "M": 1},
        "nodes": [NMR, term("rcsb_entry_info.deposited_model_count", "greater", 1)],
    },
    {
        "name": "S6",
        "title": "Membrane",
        "n": 8,
        "bins": {"M": 3, "L": 4, "XL": 1},
        "nodes": [
            either(
                term("rcsb_polymer_entity_annotation.type", "exact_match", "OPM"),
                term("rcsb_polymer_entity_annotation.type", "exact_match", "PDBTM"),
            )
        ],
    },
    {
        "name": "S7",
        "title": "Protein-DNA",
        "n": 7,
        "bins": {"S": 1, "M": 3, "L": 2, "XL": 1},
        "nodes": [DNA, PROTEIN, NO_RNA],
    },
    {
        "name": "S8",
        "title": "Protein-RNA",
        "n": 7,
        "bins": {"S": 1, "M": 2, "L": 2, "XL": 2},
        "nodes": [RNA, PROTEIN],
        "extra": [  # ribosome-scale: >= RIBOSOME_NT nucleotides of RNA, checked on the file
            {
                "tag": "ribosome",
                "n": 1,
                "nodes": [term("entity_poly.rcsb_sample_sequence_length", "greater_or_equal", 2500)],
            }
        ],
    },
    {
        "name": "S9",
        "title": "RNA only",
        "n": 3,
        "bins": {"S": 2, "M": 1},
        "nodes": [term(TYPES, "exact_match", "Nucleic acid (only)"), RNA, NO_DNA],
    },
    {
        "name": "S10",
        "title": "Glycoproteins",
        "n": 5,
        "bins": {"M": 2, "L": 2, "XL": 1},
        "nodes": [term("rcsb_entry_info.branched_entity_count", "greater", 0), PROTEIN],
    },
    {
        "name": "S11",
        "title": "One protein, two states (5 pairs)",
        "n": 10,
        "bins": {},
        # the seed pool: X-ray or EM, protein only, one protein entity, <= L, with a UniProt accession
        "nodes": [
            either(XRAY, EM),
            PROTEIN_ONLY,
            ONE_PROTEIN,
            size_range(None, BINS["L"][1]),
            term(f"{UNIPROT}.database_name", "exact_match", "UniProt"),
        ],
        "pairs": {"apo/holo": 2, "monomer/oligomer": 2, "open/closed": 1},
    },
    {
        "name": "S12",
        "title": "AlphaFold DB models",
        "n": 3,
        "bins": {"S": 1, "M": 2},
        "nodes": [term("rcsb_comp_model_provenance.source_db", "exact_match", "AlphaFoldDB")],
        "content": "computational",
    },
]
PAIR_TYPES = ("apo/holo", "monomer/oligomer", "open/closed")
PAIR_PARTNERS = 8  # partners tried per seed, at most
OPEN_CLOSED_RMSD = 3.0  # Å, CA RMSD over the common residues after superposition
QUOTAS = {"ligand": 30, "glycan": 8, "modified": 8}  # entries with the feature, at least (cross-cutting)
TOP_UP = ("S1", "S2", "S3", "S7", "S8")  # strata a short quota is topped up from, in this order
RIBOSOME_NT = 2500
BIG_CHAINS = 20  # S4 "big": more than this many polymer chains in assembly 1


def date_node(date: str) -> dict:
    return term("rcsb_accession_info.initial_release_date", "less_or_equal", date)


def stratum_nodes(stratum: dict, bin_: str | None, extra: dict | None = None) -> list[dict]:
    """The stratum's nodes, narrowed to one size bin (None: the whole stratum) and an extra quota's nodes."""
    nodes = list(stratum["nodes"])
    if bin_ is not None:
        nodes.append(size_range(*BINS[bin_]))
    elif stratum["bins"] and extra is not None:  # an extra quota over the stratum's bins
        bins = [BINS[b] for b in stratum["bins"]]
        lo = None if any(b[0] is None for b in bins) else min(b[0] for b in bins)
        nodes.append(size_range(lo, max(b[1] for b in bins)))
    return nodes + (extra["nodes"] if extra else [])


def query(nodes: list[dict], date: str, content: str = "experimental") -> dict:
    """The search request for every entry matching all nodes, released on or before date, every hit as a bare id."""
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": [*nodes, date_node(date)]},
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "results_verbosity": "compact",
            "results_content_type": [content],
        },
    }


def _post(url: str, body: dict, timeout: float = 600) -> bytes:
    req = urllib.request.Request(url, json.dumps(body).encode(), {"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as answer:
        return b"" if answer.status == 204 else answer.read()


def fetch_pool(q: dict) -> list[str]:
    """The ids a search request finds, sorted (the service answers 204 when nothing matches)."""
    raw = _post(SEARCH, q)
    if not raw:
        return []
    found = json.loads(raw)["result_set"]
    return sorted(x if isinstance(x, str) else x["identifier"] for x in found)


def pool_record(ids: list[str]) -> dict:
    """A pool as it is stored: its size and the SHA-256 of its sorted ids, one per line."""
    ids = sorted(ids)
    return {"count": len(ids), "sha256": hashlib.sha256("\n".join(ids).encode()).hexdigest()}


# ------------------------------------------------------------------------------------------------ draw


def order(pool: list[str], key: str) -> list[str]:
    """The pool, sorted, in the order Random(key) draws it: at most LOOK ids."""
    pool = sorted(pool)
    return random.Random(key).sample(pool, k=min(len(pool), LOOK))


def draw(
    pool: list[str], key: str, k: int, accept: Callable[[str], tuple[bool, str]]
) -> tuple[list[str], list[tuple[str, str]]]:
    """The first k ids of order(pool, key) that accept(id) -> (ok, reason) takes, and the ids it refused before
    them, each with its reason."""
    taken, skipped = [], []
    for i in order(pool, key):
        if len(taken) == k:
            break
        ok, reason = accept(i)
        if ok:
            taken.append(i)
        else:
            skipped.append((i, reason))
    return taken, skipped


def top_up(
    sample: list[dict],
    features: Mapping[str, set],
    rest: dict[tuple[str, str], list[str]],
    targets: dict[str, int],
    accept: Callable[[str], tuple[bool, str]],
) -> list[dict]:
    """Reach each short quota (targets: feature -> entries that must have it): within TOP_UP, in that order and
    round after round, the last drawn entry of the stratum that lacks the feature (never one drawn for an extra
    quota or a pair, never the last holder of another feature at its target) is replaced, in place, by the next id
    of its pool's seeded order (rest[(stratum, bin)], consumed) that accept takes and that has the feature.
    features: id -> feature names. Returns the replacements; a quota no stratum can reach stays short."""
    replaced = []

    def reached(f: str) -> int:
        return sum(1 for e in sample if f in features[e["id"]])

    for feature, target in targets.items():
        while reached(feature) < target:
            moved = False
            for name in TOP_UP:
                if reached(feature) >= target:
                    break
                victims = [
                    e
                    for e in sample
                    if e["stratum"] == name
                    and not e.get("quota")
                    and not e.get("pair")
                    and feature not in features[e["id"]]
                    and all(reached(f) > targets[f] or f not in features[e["id"]] for f in targets if f != feature)
                ]
                if not victims:
                    continue
                victim = max(victims, key=lambda e: e["order"])
                pool = rest.get((name, victim["bin"]), [])
                while pool:
                    cand = pool.pop(0)
                    ok, _ = accept(cand)
                    if ok and feature in features[cand]:
                        replaced.append({"id": victim["id"], "reason": f"top-up: lacks {feature}", "replacement": cand})
                        victim["id"] = cand
                        moved = True
                        break
            if not moved:
                break
    return replaced


def pairs(
    seed_pool: list[str],
    partners_of: Callable[[str], list[str]],
    classify: Callable[[str, str], set],
    wanted: dict[str, int] | None = None,
    seed: int = SEED,
) -> tuple[list[dict], list[tuple[str, str]]]:
    """S11: seeds in the seeded order of seed_pool; each seed's partners (partners_of) sorted and shuffled with
    Random(f"{seed}:S11:{seed_id}"), at most PAIR_PARTNERS of them tried; a pair is kept for the first type of
    PAIR_TYPES that classify(seed, partner) -> {types} names and that is not filled yet (wanted: type -> pairs).
    Returns the pairs ({seed, partner, type}) and the seeds that filled nothing, with the reason."""
    wanted = dict(wanted or next(s for s in STRATA if s["name"] == "S11")["pairs"])
    got: list[dict] = []
    skipped = []
    used: set[str] = set()  # an entry stands in one pair only
    for s in order(seed_pool, f"{seed}:S11"):
        open_ = [t for t in PAIR_TYPES if sum(1 for p in got if p["type"] == t) < wanted.get(t, 0)]
        if not open_:
            break
        if s in used:
            continue
        partners = sorted(set(partners_of(s)) - {s})
        random.Random(f"{seed}:S11:{s}").shuffle(partners)
        kept = None
        for p in partners[:PAIR_PARTNERS]:
            if p in used:
                continue
            types = classify(s, p)
            kind = next((t for t in open_ if t in types), None)
            if kind is not None:
                kept = {"seed": s, "partner": p, "type": kind}
                break
        if kept is None:
            skipped.append((s, "no partner of an open type" if partners else "no partner"))
        else:
            got.append(kept)
            used |= {s, kept["partner"]}
    return got, skipped


# --------------------------------------------------------------------------------------------- features

STANDARD = {
    *"ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    *"A C G U DA DC DG DT".split(),
}
NOT_MODIFIED = {"UNK", "N", "DN"}  # unknown residues are linking types, not modifications
PROTEIN_BACKBONE = {"N", "CA", "C", "O"}
NA_BACKBONE = {"P", "O5'", "C5'", "C4'", "C3'", "O3'"}
LIGAND_ATOMS = 6  # heavy atoms a ligand has at least (not water, not an ion)
TOUCH = 4.0  # Å: a ligand touches the protein


def _gemmi_structure(path):
    import gemmi

    block = gemmi.cif.read(str(path)).sole_block()
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    st.remove_hydrogens()
    return block, st


def features(path, contacts: bool = False) -> dict:
    """What the sampler and the contact sheet need to know of a structure file, read with gemmi: method, models,
    polymer chains and residues (modelled, and deposited: entity sequence x copies), protein / nucleic acid,
    RNA nucleotides, CA-only or P-only (fewer than 3 backbone atoms on more than half the residues), ligands
    (non-polymer, not water or sugar, >= LIGAND_ATOMS heavy atoms), glycans (branched entities and saccharides)
    and whether one is linked to an ASN ND2, modified polymer residues (a linking type, not one of the 20 + 8
    standard ones), polymer chains of assembly 1, UniProt accessions, and the quota features it holds (`has`).
    contacts: also whether a ligand touches the protein (within TOUCH Å)."""
    import gemmi

    block, st = _gemmi_structure(path)
    types = {
        gemmi.cif.as_string(r[0]).upper(): gemmi.cif.as_string(r[1]).upper()
        for r in block.find("_chem_comp.", ["id", "type"])
    }
    sugar = {n for n, t in types.items() if "SACCHARIDE" in t}
    poly = {}  # subchain -> polymer type
    deposited = 0
    for e in st.entities:
        if e.entity_type == gemmi.EntityType.Polymer:
            deposited += len(e.full_sequence) * len(e.subchains)
            for sc in e.subchains:
                poly[sc] = e.polymer_type
        if e.entity_type == gemmi.EntityType.Branched:
            sugar |= {gemmi.Entity.first_mon(m) for m in e.full_sequence}
    protein_types = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)
    na_types = (gemmi.PolymerType.Dna, gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid)
    model = st[0]
    residues = thin = rna = 0
    chains, kinds = set(), set()
    ligands, glycans, modified = set(), set(), set()
    ligand_res = []
    for chain in model:
        for res in chain:
            sc = res.subchain
            if sc in poly:
                kind = "protein" if poly[sc] in protein_types else "na" if poly[sc] in na_types else "other"
                kinds.add(kind)
                chains.add(chain.name)
                residues += 1
                rna += poly[sc] == gemmi.PolymerType.Rna
                names = {a.name for a in res}
                thin += len(names & (PROTEIN_BACKBONE if kind == "protein" else NA_BACKBONE)) < 3
                if res.name not in STANDARD | NOT_MODIFIED and "LINKING" in types.get(res.name, ""):
                    modified.add(res.name)
            elif res.is_water():
                continue
            elif res.name in sugar or res.entity_type == gemmi.EntityType.Branched:
                glycans.add(res.name)
            elif len(res) >= LIGAND_ATOMS:
                ligands.add(res.name)
                ligand_res.append(res)
    n_linked = any(
        ("ASN", "ND2")
        in {(c.partner1.res_id.name, c.partner1.atom_name), (c.partner2.res_id.name, c.partner2.atom_name)}
        and bool({c.partner1.res_id.name, c.partner2.res_id.name} & glycans)
        for c in st.connections
    )
    copies = len(chains)
    if st.assemblies:
        polymer_scs = set(poly)
        copies = sum(
            len(g.operators) * (sum(1 for sc in g.subchains if sc in polymer_scs) or len(set(g.chains) & chains))
            for g in st.assemblies[0].generators
        )
    accessions = sorted(
        {
            gemmi.cif.as_string(r[1])
            for r in block.find("_struct_ref.", ["db_name", "pdbx_db_accession"])
            if gemmi.cif.as_string(r[0]) == "UNP"
        }
    )
    model_name = " ".join(gemmi.cif.as_string(v) for v in block.find_values("_ma_model_list.model_group_name"))
    out = {
        "method": "; ".join(gemmi.cif.as_string(m) for m in block.find_values("_exptl.method"))
        or ("COMPUTED MODEL" if model_name else "unknown"),
        "model_name": model_name,
        "models": len(st),
        "chains": len(chains),
        "residues": residues,
        "deposited_residues": deposited,
        "protein": "protein" in kinds,
        "na": "na" in kinds,
        "rna_nt": rna,
        "ca_only": residues > 0 and thin > residues / 2,
        "ligands": sorted(ligands),
        "glycans": sorted(glycans),
        "n_linked": n_linked,
        "modified": sorted(modified),
        "assembly1_chains": copies,
        "accessions": accessions,
        "has": [k for k, v in (("ligand", ligands), ("glycan", glycans), ("modified", modified)) if v],
    }
    if contacts:
        ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
        out["ligand_touching"] = any(
            poly.get(m.to_cra(model).residue.subchain) in protein_types
            for res in ligand_res
            for atom in res
            for m in ns.find_atoms(atom.pos, "\0", radius=TOUCH)
        )
    return out


def ca_trace(path) -> dict[tuple[int, str], object]:
    """(author number, insertion code) -> CA position, of the first protein chain of the first model."""
    import gemmi

    _, st = _gemmi_structure(path)
    for chain in st[0]:
        span = chain.get_polymer()
        if len(span) and span.check_polymer_type() in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            out = {}
            for res in span:
                ca = res.find_atom("CA", "*")
                if ca is not None:
                    out[(res.seqid.num, res.seqid.icode.strip())] = ca.pos
            return out
    return {}


def ca_rmsd(path_a, path_b, least: int = 30) -> float | None:
    """CA RMSD (Å) of the first protein chains over their common residues (by author number) after superposition;
    None under `least` common residues."""
    import gemmi

    a, b = ca_trace(path_a), ca_trace(path_b)
    common = sorted(a.keys() & b.keys())
    if len(common) < least:
        return None
    sup = gemmi.superpose_positions([a[k] for k in common], [b[k] for k in common])
    return round(float(sup.rmsd), 3)


def pair_types(fa: dict, fb: dict, rmsd: float | None) -> set[str]:
    """The S11 types a pair of entries (their features, with contacts) can stand for."""
    out = set()
    if fa["ligand_touching"] != fb["ligand_touching"]:
        out.add("apo/holo")
    copies = sorted([fa["assembly1_chains"], fb["assembly1_chains"]])
    if copies[0] == 1 and copies[1] >= 2:
        out.add("monomer/oligomer")
    if fa["ligand_touching"] == fb["ligand_touching"] and rmsd is not None and rmsd >= OPEN_CLOSED_RMSD:
        out.add("open/closed")
    return out


def bin_of(residues: int) -> str:
    """The size bin of a number of deposited polymer residues ('XL+' above XL)."""
    return next((b for b, (lo, hi) in BINS.items() if (lo or 0) <= residues < hi), "XL+")


def af_accession(entry_id: str) -> str | None:
    """The UniProt accession of an RCSB computed-model id AF_AF<acc>F1 (None for a PDB id)."""
    if entry_id.startswith("AF_AF") and entry_id.endswith("F1"):
        return entry_id[5:-2]
    return None


def _commit() -> str:
    try:
        done = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT, check=True)
        return done.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _log(*parts) -> None:
    print(*parts, file=sys.stderr, flush=True)


def resolve(entry_id: str) -> Path:
    """The local file of a sampled id (downloaded once into ~/.cache/foldmap): a PDB entry, or the AlphaFold DB
    model of an AF_AF<acc>F1 id; a local path as it is."""
    from foldmap.fetch import fetch

    if Path(entry_id).is_file():
        return Path(entry_id)
    return Path(fetch(af_accession(entry_id) or entry_id))


class Sampler:
    """The network side of the draw: pools from the search service, files from foldmap.fetch, features cached."""

    def __init__(self, seed: int, date: str):
        self.seed, self.date = seed, date
        self.feats: dict[str, dict] = {}
        self.taken: set[str] = set()
        self.records: list[dict] = []  # per stratum, what sample.json says of it
        self.entries: list[dict] = []
        self.skipped: list[dict] = []
        self.replaced: list[dict] = []
        self.thin: list[str] = []  # CA-only / P-only ids met (A1c)
        self.orders: dict[tuple[str, str], list[str]] = {}  # (stratum, bin) -> its seeded order

    def pool(self, nodes: list[dict], content: str = "experimental") -> list[str]:
        return fetch_pool(query(nodes, self.date, content))

    def feature(self, i: str, contacts: bool = False) -> dict:
        f = self.feats.get(i)
        if f is None or (contacts and "ligand_touching" not in f):
            f = self.feats[i] = features(resolve(i), contacts=contacts)
        return f

    def check(self, i: str, also: Callable[[dict], str] | None = None) -> tuple[bool, str]:
        """(usable, why not): in the sample already, cannot be downloaded (obsolete), no polymer, CA/P-only, or
        `also`'s reason."""
        if i in self.taken:
            return False, "already in the sample"
        try:
            f = self.feature(i)
        except Exception as err:  # obsolete, withdrawn or unreadable: one id, not the campaign
            return False, f"download or read failed: {type(err).__name__}: {str(err)[:120]}"
        if not f["protein"] and not f["na"]:
            return False, "no protein or nucleic-acid polymer"
        if f["ca_only"]:
            if i not in self.thin:
                self.thin.append(i)
            return False, "CA-only or P-only"
        reason = also(f) if also else ""
        return (not reason), reason

    def add(self, i: str, stratum: str, bin_: str, at: int, **more) -> None:
        f = self.feats[i]
        self.taken.add(i)
        if af_accession(i):  # the AlphaFold DB file foldmap fetches for it, and the model it holds
            from foldmap.fetch import AFDB_VERSION

            more.update(source=af_accession(i), model_version=f"AlphaFold DB v{AFDB_VERSION}: {f['model_name']}")
        self.entries.append(
            {
                "id": i,
                "stratum": stratum,
                "bin": bin_,
                "order": at,
                **more,
                "method": f["method"],
                "chains": f["chains"],
                "residues": f["deposited_residues"],
                "has": f["has"],
            }
        )

    def skip(self, stratum: str, where: str, skipped: list[tuple[str, str]]) -> None:
        self.skipped += [{"id": i, "stratum": stratum, "pool": where, "reason": r} for i, r in skipped]

    def stratum(self, s: dict) -> dict[tuple[str, str], list[str]]:
        """Draw one stratum (not S11): its extra quotas first, then each bin. Returns what is left of each bin's
        seeded order (for the top-up)."""
        name, content = s["name"], s.get("content", "experimental")
        rec = {
            "name": name,
            "title": s["title"],
            "n": s["n"],
            "query": query(stratum_nodes(s, None), self.date, content),
            "bins": {},
            "extra": [],
        }
        pools = {}
        for b, k in s["bins"].items():
            pools[b] = self.pool(stratum_nodes(s, b), content)
            lo, hi = BINS[b]
            rec["bins"][b] = {"n": k, "residues": [lo or 0, hi - 1], **pool_record(pools[b])}
            _log(f"{name} {b}: pool {len(pools[b])}")
        room = dict(s["bins"])
        also = self.extra_check(s)
        for x in s.get("extra", []):
            sub = self.pool(stratum_nodes(s, x.get("bin"), x), content)
            rec["extra"].append({"tag": x["tag"], "n": x["n"], "nodes": x["nodes"], **pool_record(sub)})
            key = f"{self.seed}:{name}:{x['tag']}"
            where = {i: b for b in pools for i in pools[b]}
            check = self.extra_check(s, x["tag"])

            def accept(i, check=check, where=where):
                if i not in where or room[where[i]] <= 0:
                    return False, "not in a bin of the stratum" if i not in where else f"bin {where[i]} full"
                return self.check(i, check)

            taken, skipped = draw(sub, key, x["n"], accept)
            seen = order(sub, key)
            for i in taken:
                room[where[i]] -= 1
                self.add(i, name, where[i], seen.index(i), quota=x["tag"])
            self.skip(name, x["tag"], skipped)
            _log(f"{name} {x['tag']}: pool {len(sub)}, took {taken}")
        rest = {}
        for b, pool in pools.items():
            key = f"{self.seed}:{name}:{b}"
            taken, skipped = draw(pool, key, room[b], lambda i: self.check(i, also))
            seen = self.orders[(name, b)] = order(pool, key)
            for i in taken:
                self.add(i, name, b, seen.index(i))
            self.skip(name, b, skipped)
            last = max((seen.index(i) for i in taken), default=-1)
            last = max([last] + [seen.index(i) for i, _ in skipped])
            rest[(name, b)] = seen[last + 1 :]
            _log(f"{name} {b}: took {taken}")
        self.records.append(rec)
        return rest

    def extra_check(self, s: dict, tag: str | None = None) -> Callable[[dict], str] | None:
        """The checks on the file a stratum or extra quota adds."""
        if s["name"] == "S10":
            return lambda f: "" if f["n_linked"] else "no N-linked glycan (ASN ND2)"
        if tag == "big":
            return lambda f: (
                "" if f["assembly1_chains"] > BIG_CHAINS else f"assembly 1 has {f['assembly1_chains']} chains"
            )
        if tag == "ribosome":
            return lambda f: "" if f["rna_nt"] >= RIBOSOME_NT else f"{f['rna_nt']} nt of RNA"
        return None

    def pairs(self, s: dict) -> None:
        """S11: pairs of entries of one protein (one UniProt accession) in two states."""
        seed_pool = self.pool(s["nodes"])
        rec = {"name": s["name"], "title": s["title"], "n": s["n"], "query": query(s["nodes"], self.date)}
        rec["seed_pool"] = pool_record(seed_pool)
        rec["pairs"] = dict(s["pairs"])
        _log(f"S11 seed pool {len(seed_pool)}")
        partners: dict[str, list[str]] = {}

        def partners_of(i: str) -> list[str]:
            ok, why = self.check(i)
            if not ok:
                return []
            acc = self.feats[i]["accessions"]
            if len(acc) != 1:
                return []
            nodes = [*s["nodes"], term(f"{UNIPROT}.database_accession", "exact_match", acc[0])]
            partners[i] = self.pool(nodes)
            return partners[i]

        def classify(a: str, b: str) -> set[str]:
            ok, _ = self.check(b)
            if not ok:
                return set()
            fa, fb = self.feature(a, contacts=True), self.feature(b, contacts=True)
            if fb["accessions"] != fa["accessions"]:
                return set()
            types = pair_types(fa, fb, None)
            if fa["ligand_touching"] == fb["ligand_touching"]:
                types = pair_types(fa, fb, ca_rmsd(resolve(a), resolve(b)))
            return types

        got, skipped = pairs(seed_pool, partners_of, classify, s["pairs"], self.seed)
        seen = order(seed_pool, f"{self.seed}:S11")
        for k, p in enumerate(got, 1):
            rmsd = ca_rmsd(resolve(p["seed"]), resolve(p["partner"]))
            for role in ("seed", "partner"):
                i = p[role]
                self.add(
                    i,
                    "S11",
                    bin_of(self.feats[i]["deposited_residues"]),
                    seen.index(p["seed"]),
                    pair=f"P{k}",
                    pair_type=p["type"],
                    role=role,
                    partner=p["partner"] if role == "seed" else p["seed"],
                    ligand_touching=self.feats[i]["ligand_touching"],
                    assembly1_chains=self.feats[i]["assembly1_chains"],
                    ca_rmsd=rmsd,
                )
            _log(f"S11 pair P{k} {p['type']}: {p['seed']} {p['partner']}")
        self.skip("S11", "seeds", skipped)
        self.records.append(rec)


def sample(seed: int, date: str) -> dict:
    """Draw the whole sample (network): every stratum, the pairs, then the cross-cutting top-up."""
    s = Sampler(seed, date)
    rest: dict = {}
    for st in STRATA:
        if st["name"] == "S11":
            s.pairs(st)
        else:
            rest.update(s.stratum(st))

    class Has(dict):
        def __missing__(self, i):
            return set(s.feature(i)["has"])

    def accept(i: str) -> tuple[bool, str]:  # TOP_UP strata add no check on the file
        if any(e["id"] == i for e in s.entries):
            return False, "already in the sample"
        return s.check(i)

    s.replaced = top_up(s.entries, Has(), rest, QUOTAS, accept)
    for r in s.replaced:
        e = next(e for e in s.entries if e["id"] == r["replacement"])
        f = s.feats[e["id"]]
        e.update(
            order=s.orders[(e["stratum"], e["bin"])].index(e["id"]),
            method=f["method"],
            chains=f["chains"],
            residues=f["deposited_residues"],
            has=f["has"],
        )
        s.taken.discard(r["id"])
        s.taken.add(e["id"])
    s.replaced += [{"id": i, "reason": "CA-only or P-only (run once for A1c)", "replacement": None} for i in s.thin]
    reached = {f: sum(1 for e in s.entries if f in e["has"]) for f in QUOTAS}
    quotas = {f: {"target": t, "reached": reached[f]} for f, t in QUOTAS.items()}
    for st in STRATA:
        for x in st.get("extra", []):
            n = sum(1 for e in s.entries if e["stratum"] == st["name"] and e.get("quota") == x["tag"])
            quotas[f"{st['name']}:{x['tag']}"] = {"target": x["n"], "reached": n}
        for b, k in st["bins"].items():
            n = sum(1 for e in s.entries if e["stratum"] == st["name"] and e["bin"] == b)
            quotas[f"{st['name']}:{b}"] = {"target": k, "reached": n}
        for t, k in st.get("pairs", {}).items():
            n = sum(1 for e in s.entries if e.get("pair_type") == t and e.get("role") == "seed")
            quotas[f"S11:{t}"] = {"target": k, "reached": n}
    names = [st["name"] for st in STRATA]
    entries = sorted(s.entries, key=lambda e: (names.index(e["stratum"]), e.get("pair", ""), e["bin"], e["order"]))
    return {
        "seed": seed,
        "query_date": date,
        "foldmap_commit": _commit(),
        "look": LOOK,
        "bins": {b: [lo or 0, hi - 1] for b, (lo, hi) in BINS.items()},
        "strata": s.records,
        "entries": entries,
        "replaced": s.replaced,
        "skipped": s.skipped,
        "quotas": quotas,
    }


def dump(doc: dict) -> str:
    return json.dumps(doc, indent=1, ensure_ascii=False) + "\n"


def obsolete(doc: dict) -> list[str]:
    """The sampled ids (entries and replacements) the search service no longer finds (network)."""
    ids = sorted({e["id"] for e in doc["entries"]} | {r["replacement"] for r in doc["replaced"] if r["replacement"]})
    found: set[str] = set()
    for content in ("experimental", "computational"):
        group = [i for i in ids if (af_accession(i) is not None) == (content == "computational")]
        if group:
            nodes = [term("rcsb_entry_container_identifiers.entry_id", "in", group)]
            found |= set(fetch_pool(query(nodes, doc["query_date"], content)))
    return [i for i in ids if i not in found]


def quota_table(doc: dict) -> str:
    rows = [f"{'quota':<22} target reached"]
    for name, q in doc["quotas"].items():
        mark = "" if q["reached"] >= q["target"] else "  SHORT"
        rows.append(f"{name:<22} {q['target']:>6} {q['reached']:>7}{mark}")
    return "\n".join(rows)


def cmd_sample(args) -> int:
    out = Path(args.out)
    if args.check:
        doc = json.loads(out.read_text(encoding="utf-8"))
        gone = obsolete(doc)
        print(f"{len(gone)} obsolete id(s){': ' + ', '.join(gone) if gone else ''}")
        return 1 if gone else 0
    doc = sample(args.seed, args.date)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(dump(doc), encoding="utf-8")
    print(quota_table(doc))
    short = [k for k, q in doc["quotas"].items() if q["reached"] < q["target"]]
    print(f"{len(doc['entries'])} entries written to {out}" + (f"; short: {', '.join(short)}" if short else ""))
    return 1 if short or len(doc["entries"]) != 100 else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("sample", help="draw the stratified sample (network)")
    s.add_argument("--seed", type=int, default=SEED)
    s.add_argument("--date", default=datetime.date.today().isoformat(), help="release-date cutoff (YYYY-MM-DD)")
    s.add_argument("--out", default=str(OUT / "sample.json"))
    s.add_argument("--check", action="store_true", help="list the sampled ids now obsolete; change nothing")
    s.set_defaults(func=cmd_sample)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
