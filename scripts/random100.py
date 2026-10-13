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

    python scripts/random100.py run [--ids A,B | --failed-from report.json] [--views all|quick|v1,v2]
                                    [--out docs/validation/random100] [--jobs 1]

run draws every entry (plus each CA-only / P-only id met while sampling, which must be refused: A1c) with
`entry <id> --out DIR`, a child process (python -I, OPENBLAS_NUM_THREADS=1, PYTHONHASHSEED=0). The entry draws each
view (VIEWS) in a grandchild of its own (`view`, run with -W default), killed after LIMITS["view_timeout_s"], so a
crash or a hang in one view does not hide the others; wall time and max RSS come from os.wait4. Each view writes
what it measured (figcheck.figure_checks, route statistics from foldmap.route's DEBUG summary record, layout
signature and SVG hashes, na_validate's nucleic-acid measures, membrane and tether checks); the entry adds mkdssp
agreement and judges A1-A20 (criteria, LIMITS) into <out>/<id>.json. Then the pairs are compared (A19), a contact
sheet <out>/<id>.png is drawn per entry (the figure, the 3D CA/P traces seen through the figure's page axes in its
element colours, the same turned 90° about page-up, the partner's figure for a pair, facts and criteria), and
report.json, report.md and index.html are written. Exit 0 only if every entry passes. Each view's files are kept in
<out>/work/<id>/. Pairing against RNApolis needs NA_RNAPOLIS_PYTHON (the Python of a venv that has it); mkdssp 3.1.2
is found on PATH or as MKDSSP.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import json
import logging
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request
from collections.abc import Callable, Mapping
from pathlib import Path

import numpy as np

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


# ------------------------------------------------------------------------------------------------ runner

SCRIPT = Path(__file__).resolve()
LIMITS = {
    "view_timeout_s": 900,  # a view's process is killed after this many seconds (A1 then fails)
    "default_s": 60,  # A2: wall seconds of the default view
    "default_s_many_chains": 180,  # ... when more than many_chains polymer chains are read (Q-R3)
    "many_chains": 20,
    "app_s": 90,  # A2b: the hosted app's limits (app_jobs.HOSTED, Q-R7)
    "app_memory_mb": 1024,
    "app_max_residues": 30_000,
    "dssp_agree": 0.99,  # A10: fraction of residues whose SS agrees with mkdssp 3.1.2
    "copies_alike": 0.5,  # A15: page units (Q-R12)
    "domain_overlap": 0.02,  # A16: fraction of the smaller panel
    # A17: page units from a stem name to its own letters, at most layout._INK_MARGIN em in the drawn figure's page
    # units plus 1e-3 (figcheck.stem_gap_tolerance; plan 10 Task 4 fix round 2 ruling), measured per figure
    "stem_name_gap": "figcheck.stem_gap_tolerance(layout._INK_MARGIN em)",
    "tm_overlap": 0.75,  # A18: a TM helix overlaps the band by this fraction of the band's height
    "slab_margin": 0.3,  # A18: a loop's sideways run inside the band, this far from its faces, is in the slab
    "tether_length": 6.5,  # A20: page units, a ligand tether from its mark to the element it touches most
}
DOCUMENTED_WARNINGS: tuple[str, ...] = ()  # Python warnings a view may print (regexes); none at the start
VIEWS = (
    "default",
    "deposited",
    "beta",
    "focus",
    "domains",
    "membrane",
    "seed1",
    "roundtrip",
    "interactive",
    "sequence",
    "app",
)
QUICK = ("default", "seed1", "roundtrip")
FIGURE_VIEWS = ("default", "deposited", "beta", "focus", "domains", "membrane")
LABEL_VIEWS = ("default", "deposited", "beta", "focus")
CRITERIA = {
    "A1": "no crash: exit 0, files written, no traceback, no Python warning",
    "A1c": "a CA-only or P-only entry exits non-zero with foldmap's one-line message",
    "A2": "default view time",
    "A2b": "the app draws it within the hosted limits",
    "A3": "no route fallbacks",
    "A4": "no loop over a label",
    "A5": "no element or stem boxes overlapping",
    "A6": "no text collisions",
    "A7": "text and legend inside the page",
    "A8": "labels unique per chain, stem names per figure",
    "A9": "nothing left out without a note",
    "A10": "secondary structure agrees with mkdssp",
    "A11": "N->C and sheet order",
    "A12": "same figure for any PYTHONHASHSEED",
    "A13": "a saved layout redraws the same figure",
    "A14": "interactive page and sequence plot written",
    "A15": "copies drawn alike",
    "A16": "domains: one chain, alike, not tiny, panels apart, names clear",
    "A17": "nucleic-acid metrics and stem names by their letters",
    "A18": "membrane: TM helices in the band, loops out of the slab",
    "A19": "two states: same names in the same order",
    "A20": "ligand tethers local",
}
_WARNING = re.compile(r"^\S.*:\d+: (\w*Warning): (.*)$")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _jsonable(x):
    """x with tuple keys joined by '~', sets and tuples as lists, numpy numbers as Python ones."""
    if isinstance(x, dict):
        return {("~".join(map(str, k)) if isinstance(k, tuple) else str(k)): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set, frozenset)):
        items = sorted(x, key=str) if isinstance(x, (set, frozenset)) else x
        return [_jsonable(v) for v in items]
    if hasattr(x, "item") and not isinstance(x, (str, bytes)):
        return x.item()
    if isinstance(x, float) and x != x:
        return None
    return x


def _script_module(name: str):
    """A sibling script (na_validate, bench) as a module, without putting scripts/ on sys.path."""
    import importlib.util

    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT.parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def run_child(cmd: list[str], *, env: dict | None = None, timeout: float, out: Path, err: Path) -> dict:
    """Run cmd in its own session, stdout and stderr to files; killed (with everything it started) after timeout
    seconds. Exit code, wall seconds and the child's own max RSS (os.wait4: its rusage, the waited-for children it
    reaped included)."""
    t0 = time.perf_counter()
    with open(out, "wb") as fo, open(err, "wb") as fe:
        proc = subprocess.Popen(cmd, stdout=fo, stderr=fe, env=env, start_new_session=True)
        timed_out = False
        while True:
            pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.perf_counter() - t0 > timeout:
                timed_out = True
                os.killpg(proc.pid, signal.SIGKILL)
                pid, status, usage = os.wait4(proc.pid, 0)
                break
            time.sleep(0.05)
    proc.returncode = os.waitstatus_to_exitcode(status)
    return {
        "exit": proc.returncode,
        "seconds": round(time.perf_counter() - t0, 2),
        "max_rss_mb": round(usage.ru_maxrss / 1024, 1),
        "timed_out": timed_out,
    }


def view_argv(view: str, path: Path, name: str, work: Path, focus: str | None = None) -> list[str]:
    """The foldmap command line of a view (after `foldmap`); [] for the app view (not a command)."""
    stem = work / f"{name}-{view}"
    figure = [str(path), "-o", f"{stem}.svg", "-o", f"{stem}.png"]
    return {
        "default": ["plot", *figure, "--save-layout", str(work / f"{name}-layout.yaml")],
        "deposited": ["plot", *figure, "--names", "deposited"],
        "beta": ["plot", *figure, "--names", "beta"],
        "focus": ["plot", *figure, "--focus", focus or "none"],
        "domains": ["plot", *figure, "--domains", "auto"],
        "membrane": ["plot", *figure, "--membrane", "auto"],
        "seed1": ["plot", *figure],
        "roundtrip": ["plot", *figure, "--layout-file", str(work / f"{name}-layout.yaml")],
        "interactive": ["interactive", str(path), "-o", str(work / f"{name}.html")],
        "sequence": ["sequence", str(path), "-o", f"{stem}.svg"],
        "app": [],
    }[view]


def view_command(view: str, path: Path, name: str, work: Path, focus: str | None = None) -> tuple[list[str], dict]:
    """The command and environment of one view's process: python -I, warnings shown (-W default), one BLAS thread,
    PYTHONHASHSEED 0 (1 for seed1)."""
    cmd = [sys.executable, "-I", "-W", "default", str(SCRIPT), "view", "--name", name, "--work", str(work)]
    cmd += [view, str(path), "--", *view_argv(view, path, name, work, focus)]
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "PYTHONHASHSEED": "1" if view == "seed1" else "0"}
    return cmd, env


class _RouteLog(logging.Handler):
    """The summary records route.route_loops logs at DEBUG, summed over every routing of the figure."""

    KEYS = (  # the record's arguments in order; "inflate" is the clearance it names, not a count
        "loops",
        "fallback_curves",
        "fallback_curves_crossing_an_element",
        "kink_links_rejected",
        "tighter_than_inflate",
        "inflate",
        "rejected_as_crossing",
    )

    def __init__(self):
        super().__init__(logging.DEBUG)
        self.sums = {k: 0 for k in self.KEYS if k != "inflate"}
        self.records = 0

    def emit(self, record):
        if not str(record.msg).startswith("%d loops:") or len(record.args) != len(self.KEYS):
            return
        self.records += 1
        for k, v in zip(self.KEYS, record.args):
            if k in self.sums:
                self.sums[k] += int(v)

    def summary(self) -> dict:
        return {**self.sums, "records": self.records}


def _colour(c) -> str:
    from matplotlib.colors import to_hex

    return to_hex(c, keep_alpha=False)


def element_colours(fig) -> tuple[dict[str, str], dict[int, str]]:
    """Each element's colour (helix:, strand:, eta: bodies) and each stem's (na:stem:<k>: block, rails or front),
    read from the figure's artists by gid."""
    elements, stems = {}, {}
    if not fig.axes:
        return elements, stems
    for a in fig.axes[0].get_children():
        gid = a.get_gid() or ""
        m = re.fullmatch(r"(helix|strand|eta):(.+)", gid)
        if m:
            face = a.get_facecolor() if hasattr(a, "get_facecolor") else None
            edge = a.get_edgecolor() if hasattr(a, "get_edgecolor") else None
            c = face if face is not None and len(face) == 4 and face[3] > 0 else edge
            if c is not None:
                elements.setdefault(m.group(2), _colour(c))
            continue
        m = re.fullmatch(r"na:stem:(\d+):(block|rails|front)", gid)
        if m:
            if hasattr(a, "get_facecolor") and m.group(2) == "block":
                stems.setdefault(int(m.group(1)), _colour(a.get_facecolor()))
            elif hasattr(a, "get_color"):
                stems.setdefault(int(m.group(1)), _colour(a.get_color()))
            elif hasattr(a, "get_edgecolor"):
                stems.setdefault(int(m.group(1)), _colour(a.get_edgecolor()))
    return elements, stems


def trace_of(drawn) -> dict:
    """The 3D traces for the contact sheet: each protein chain's CA and each nucleic-acid strand's P, with each
    residue's colour as the figure draws its element or stem (None outside elements), the page axes and centre."""
    from foldmap.layout import stem_pieces

    lay, bb = drawn.layout, drawn.bb
    elements, stems = element_colours(drawn.figure)
    colour = [None] * len(bb)
    for s in drawn.sses:
        if s.id in elements:
            for k in range(s.start, s.end + 1):
                colour[k] = elements[s.id]
    chains = []
    for c in dict.fromkeys(l.chain for l in bb.labels):
        idx = [k for k, l in enumerate(bb.labels) if l.chain == c]
        chains.append(
            {
                "chain": c,
                "xyz": [[round(float(v), 2) for v in bb.ca[k]] for k in idx],
                "joined": [bool(bb.prev[k]) for k in idx],
                "colour": [colour[k] for k in idx],
            }
        )
    strands = []
    na = lay.nucleic
    if na is not None:
        nt_colour: dict = {}
        pieces = {p.id: p for p in stem_pieces(na)}
        for k, d in enumerate(lay.stems):
            piece = pieces.get(d.id)
            if piece is not None and k in stems:
                for x in piece.nucleotides():
                    nt_colour[tuple(x)] = stems[k]
        for s, strand in enumerate(na.strands):
            strands.append(
                {
                    "chain": strand.chain,
                    "xyz": [[round(float(v), 2) for v in p] for p in strand.p],
                    "colour": [nt_colour.get((s, i)) for i in range(len(strand.p))],
                }
            )
    return {"axes": lay.page_axes, "centre": lay.page_centre, "chains": chains, "strands": strands}


def tether_lengths(drawn, look) -> list[tuple[str, float]]:
    """(gid, length in page units) of every ligand tether drawn: from its ligand's mark to its target (render's
    ligand_marks; the drawn tether shows only where no body or text is, so its ink is shorter)."""
    from foldmap.render import ligand_marks

    out = []
    for m in ligand_marks(drawn.layout, drawn.loops, drawn.sses, look):
        g, (cx, cy) = m["ligand"], m["centre"]
        for k, (tx, ty) in enumerate(m["targets"]):
            out.append((f"ligand-tether:{g.name}:{g.chain}{g.seq}:{k}", round(float(np.hypot(tx - cx, ty - cy)), 3)))
    return out


def membrane_checks(drawn, path) -> dict:
    """A18 on a --membrane auto figure: each TM helix (its CA run reaches both membrane faces, as cli._face_medians
    counts it) overlaps its band by at least tm_overlap of the band's height; loops run sideways inside a band
    (slab_margin from its faces) only next to one of their own ports inside it."""
    from foldmap.membrane import find_membrane

    lay, bb = drawn.layout, drawn.bb
    if lay.membrane is None:
        return {"band": False, "tm_helices": 0, "tm_short": [], "loops_in_slab": []}
    mem = find_membrane(path, bb)
    bands = lay.membrane.get("bands") or [lay.membrane["y"]]
    depth = mem.depth(bb.ca) if mem is not None else None
    short, tm = [], 0
    for s in drawn.sses:
        p = lay.placed.get(s.id)
        if p is None or s.kind != "H" or depth is None:
            continue
        seg = depth[s.start : s.end + 1]
        if len(seg) < 2 or seg.min() > -mem.half + 2 or seg.max() < mem.half - 2:
            continue
        tm += 1
        lo, hi = sorted((p.n_port[1], p.c_port[1]))
        y0, y1 = min(bands, key=lambda b: abs((b[0] + b[1]) / 2 - (lo + hi) / 2))
        frac = (min(hi, y1) - max(lo, y0)) / max(y1 - y0, 1e-9)
        if frac < LIMITS["tm_overlap"]:
            short.append(f"{s.id} {frac:.2f}")
    margin = LIMITS["slab_margin"]

    def inside(y):
        return any(lo + margin < y < hi - margin for lo, hi in bands)

    in_slab = []
    for lp in drawn.loops:
        pts, n = lp.points, len(lp.points)
        for k, ((xa, ya), (xb, yb)) in enumerate(zip(pts, pts[1:])):
            if abs(ya - yb) < 1e-9 and abs(xa - xb) > 0.5 and inside(ya):
                if not ((k <= 1 and inside(pts[0][1])) or (k >= n - 3 and inside(pts[-1][1]))):
                    in_slab.append(f"{lp.a_id}>{lp.b_id}")
                    break
    return {"band": True, "tm_helices": tm, "tm_short": short, "loops_in_slab": in_slab}


def na_metrics(drawn, path) -> dict:
    """na_validate's measures of the figure's nucleic acid (pairing against RNApolis when NA_RNAPOLIS_PYTHON names
    its Python, shape, contacts, render checks), with the metrics past nacheck.THRESHOLDS under "failed"."""
    from foldmap.io import read_model
    from foldmap.nacheck import contacts, judge, render_checks, shape

    nav = _script_module("na_validate")
    lay = drawn.layout
    m: dict = {"pair_f1": None, "stem_recall": None}
    try:
        theirs = nav.reference_pairs(Path(path))
    except (RuntimeError, subprocess.TimeoutExpired) as err:
        theirs, m["note"] = None, str(err)
    else:
        m.update(nav.pairing(Path(path), theirs))
        asu_model, deposited = read_model(Path(path), "asu")
        model = drawn.model
        as_deposited = model[1] <= deposited and {c.name for c in model[0]} == model[1]
        m.update(nav.drawn_pairing(lay.nucleic, theirs, asu_model, as_deposited))
    m.update(shape(lay))
    m.update(contacts(lay, drawn.sses))
    m.update(render_checks(lay, drawn.loops, drawn.figure, drawn.sses, model=drawn.model, bb=drawn.bb))
    m["failed"] = judge(m)
    if m.get("unmapped_chains"):
        m["failed"].append(
            f"pairing: drawn chains not mapped to the deposited model: {', '.join(m['unmapped_chains'])}"
        )
    if m.get("pair_f1") is None:
        m["failed"].append(f"pairing not measured ({m.get('note', 'no reason given')})")
    return m


def stem_gap_limit(fig) -> float | None:
    """A17's limit for this figure: the stem-name margin (layout._INK_MARGIN em) in its page units, plus 1e-3
    (figcheck.stem_gap_tolerance). The em is read from a drawn stem name (render: 0.8 em) as points -> pixels ->
    page units, as tests/test_na_layout.py reads it from a base letter; None when no stem name is drawn."""
    from foldmap import figcheck
    from foldmap.layout import _INK_MARGIN, _NAME_EM

    if not fig.axes:
        return None
    ax = fig.axes[0]
    name = next((t for t in ax.texts if (t.get_gid() or "").startswith("na:label:")), None)
    if name is None:
        return None
    px_per_unit = ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    em = name.get_fontsize() / _NAME_EM * fig.dpi / 72 / px_per_unit
    return figcheck.stem_gap_tolerance(_INK_MARGIN * em)


def measure_view(view: str, drawn, spec, path: Path, name: str, work: Path) -> dict:
    """Everything one drawn figure view is judged on."""
    from foldmap import figcheck
    from foldmap.scope import drawn_scope

    bench = _script_module("bench")
    lay = drawn.layout
    checks = figcheck.figure_checks(drawn, path, assembly=spec.assembly, focus=spec.focus, chains=spec.chains)
    sig = json.dumps(_jsonable(bench.signature(lay)), sort_keys=True, default=str)
    loops = json.dumps([[l.a_id, l.b_id, [list(map(float, p)) for p in l.points]] for l in drawn.loops])
    svg = work / f"{name}-{view}.svg"
    out = {
        "checks": _jsonable(checks),
        "stuck": drawn.stuck,
        "over_labels": drawn.over_labels,
        "signature_sha256": _sha256(sig.encode()),
        "loops_sha256": _sha256(loops.encode()),
        "svg_sha256": _sha256(svg.read_bytes()) if svg.is_file() else None,
        "scope": _jsonable(drawn_scope(drawn.model, lay, spec.assembly, spec.focus)),
        "chains_read": list(lay.chains_read),
        "protein_chains_drawn": list(dict.fromkeys(p.sse.chain for p in lay.placed.values())),
        "first_protein_chain": drawn.bb.labels[0].chain if len(drawn.bb) else None,
        "na_chains_drawn": sorted({s.chain for s in lay.nucleic.strands}) if lay.nucleic is not None else [],
        "elements": [
            {
                "id": p.sse.id,
                "chain": p.sse.chain,
                "kind": p.sse.kind,
                "first": p.sse.first.seq,
                "last": p.sse.last.seq,
                "label": p.label,
            }
            for p in lay.placed.values()
        ],
        "area": round(lay.width * lay.height, 2),
        "context": lay.context,
        "copy_ties": len(lay.copy_ties),
        "tethers": tether_lengths(drawn, spec.look),
        "stem_gap_limit": stem_gap_limit(drawn.figure),
    }
    if view == "membrane":
        out["membrane"] = membrane_checks(drawn, path)
    if view == "default":
        if lay.nucleic is not None:
            out["na"] = _jsonable(na_metrics(drawn, path))
        (work / f"{name}-trace.json").write_text(json.dumps(trace_of(drawn)), encoding="utf-8")
    return out


def app_view(path: Path, name: str) -> dict:
    """The app's way: the size check, then app_jobs.draw_job in a job process under the hosted limits, with the
    app's default settings."""
    from foldmap import app, app_jobs

    lim = app_jobs.limits(hosted=True)
    size = app_jobs.run_job(app_jobs.size_job, str(path), limits=lim)
    if lim.max_residues and size["residues"] > lim.max_residues:
        return {"refused": "size", "residues": size["residues"], "message": app_jobs.size_message(name, size, lim)}
    settings = json.dumps({**app.LAYOUT_DEFAULTS, "theme": "publication"}, sort_keys=True)
    try:
        got = app_jobs.run_job(app_jobs.draw_job, str(path), settings, lim.slow, limits=lim)
    except app_jobs.JobLimit as err:
        return {"refused": err.kind, "residues": size["residues"], "message": app_jobs.limit_message(name, err, lim)}
    return {"refused": None, "residues": size["residues"], "svg_bytes": len(got["svg"])}


def cmd_view(args) -> int:
    """One view of one entry, in its own process (run by `entry`): the foldmap command, with what it drew kept
    (cli.draw_spec wrapped) and route's summary records caught, then measured. Writes <work>/<name>-<view>.json;
    exits with the command's exit code (an exception escaping foldmap prints its traceback here, as foldmap's)."""
    import logging as _logging
    import warnings

    from foldmap import cli

    work, path = Path(args.work), Path(args.path)
    caught: dict = {}
    real = cli.draw_spec

    def keep(spec):
        drawn = real(spec)
        caught["drawn"], caught["spec"] = drawn, spec
        return drawn

    cli.draw_spec = keep
    route_log = _RouteLog()
    logger = _logging.getLogger("foldmap.route")
    logger.addHandler(route_log)
    logger.setLevel(_logging.DEBUG)
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as said:
        warnings.simplefilter("default")
        if args.view == "app":
            result: dict = {"app": app_view(path, args.name)}
            code = 0
        else:
            result = {}
            try:
                code = cli.main(args.argv)
            except SystemExit as stop:  # argparse
                code = stop.code if isinstance(stop.code, int) else 1
    result.update(
        view=args.view,
        argv=args.argv,
        exit=code,
        seconds_inside=round(time.perf_counter() - t0, 2),
        warnings=[f"{w.category.__name__}: {w.message}" for w in said],
        route=route_log.summary(),
    )
    logger.removeHandler(route_log)
    if code == 0 and "drawn" in caught:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # the measuring's own warnings are not the figure's
            try:
                result.update(measure_view(args.view, caught["drawn"], caught["spec"], path, args.name, work))
            except Exception:  # the measure failed: the view is not measured, which fails it (fail closed)
                import traceback

                result["measure_error"] = traceback.format_exc()[-2000:]
    (work / f"{args.name}-{args.view}.json").write_text(json.dumps(_jsonable(result), indent=1), encoding="utf-8")
    return code


# --------------------------------------------------------------------------------------------- dssp


def mkdssp_rows(path: Path) -> list[tuple[str, int, str, str]] | str:
    """mkdssp 3.1.2's per-residue codes for the file (make_dssp_reference.parse_dssp), or why there are none."""
    exe = os.environ.get("MKDSSP") or shutil.which("mkdssp")
    if not exe:
        return "no mkdssp on PATH"
    try:
        version = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60).stdout
    except (OSError, subprocess.SubprocessError) as err:
        return f"cannot run mkdssp: {err}"
    if "mkdssp version 3.1.2" not in version:
        return f"mkdssp 3.1.2 is required (found: {version.strip()[:60] or 'no version output'})"
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(path)
        if src.suffix == ".gz":
            plain = Path(tmp) / src.name[:-3]
            plain.write_bytes(gzip.decompress(src.read_bytes()))
            src = plain
        dssp = Path(tmp) / "out.dssp"
        try:
            done = subprocess.run([exe, str(src), str(dssp)], capture_output=True, text=True, timeout=600)
        except subprocess.TimeoutExpired:
            return "mkdssp took longer than 600 s"
        if done.returncode or not dssp.is_file():
            tail = (done.stderr.strip().splitlines() or ["failed"])[-1]
            return f"mkdssp cannot read it: {tail[:160]}"
        text = dssp.read_text(encoding="utf-8", errors="replace")
    return _script_module("make_dssp_reference").parse_dssp(text)


def dssp_agreement(path: Path, rows: list[tuple[str, int, str, str]] | None = None) -> dict:
    """foldmap's DSSP on the model as deposited against mkdssp's codes (rows; None: run mkdssp), residue by residue,
    as tests/test_dssp_reference.disagreements compares them: a residue standing in for a fused chromophore, and the
    four on either side of it, are not compared. agree: the fraction of the residues both read that agree (None,
    with `unmeasured` saying why, when mkdssp cannot read the file)."""
    from foldmap.dssp import assign_dssp
    from foldmap.io import load_backbone

    if rows is None:
        rows = mkdssp_rows(path)
        if isinstance(rows, str):
            return {"agree": None, "unmeasured": rows}
    ref = {(c, n, i): code for c, n, i, code in rows}
    bb = load_backbone(Path(path), "asu")
    ss = assign_dssp(bb).ss
    keys = [(l.chain, l.seq, l.icode) for l in bb.labels]
    stood_in = {keys[k] for k in bb.linked}
    near = {keys[j] for k in bb.linked for j in range(k - 4, k + 5) if 0 <= j < len(keys)} - stood_in
    compared = [(k, key) for k, key in enumerate(keys) if key in ref and key not in near and key not in stood_in]
    differ = [
        f"{c}{n}{i} mkdssp {ref[(c, n, i)]} foldmap {ss[k]}" for k, (c, n, i) in compared if ss[k] != ref[(c, n, i)]
    ]
    if not compared:
        return {"agree": None, "unmeasured": "no residue read by both"}
    return {
        "agree": round(1 - len(differ) / len(compared), 4),
        "residues": len(compared),
        "foldmap_only": len(set(keys) - stood_in - set(ref)),
        "mkdssp_only": len(set(ref) - set(keys)),
        "differ": differ[:20],
    }


# --------------------------------------------------------------------------------------------- entry


def entry_name(source: str) -> str:
    """The name an entry's files carry: the id, or a local file's name without its extensions."""
    return Path(source).name.split(".")[0] if Path(source).is_file() else source


def _outputs(view: str, name: str, work: Path) -> list[Path]:
    stem = work / f"{name}-{view}"
    if view == "interactive":
        return [work / f"{name}.html"]
    if view == "sequence":
        return [stem.with_suffix(".svg")]
    if view == "app":
        return []
    return [stem.with_suffix(".svg"), stem.with_suffix(".png")]


def run_view(view: str, path: Path, name: str, work: Path, focus: str | None = None) -> dict:
    """One view in its own process: what it exited with, its time and memory, what it printed, what it measured."""
    cmd, env = view_command(view, path, name, work, focus)
    err = work / f"{name}-{view}.err"
    r = run_child(cmd, env=env, timeout=LIMITS["view_timeout_s"], out=work / f"{name}-{view}.out", err=err)
    stderr = err.read_text(encoding="utf-8", errors="replace") if err.is_file() else ""
    lines = stderr.splitlines()
    r["traceback"] = "Traceback (most recent call last)" in stderr
    r["python_warnings"] = [f"{m.group(1)}: {m.group(2)}" for m in map(_WARNING.match, lines) if m]
    r["foldmap_said"] = [line for line in lines if line.startswith("foldmap:")]
    r["stderr_tail"] = [line[:300] for line in lines[-6:]]
    r["written"] = {p.name: p.stat().st_size if p.is_file() else 0 for p in _outputs(view, name, work)}
    got = work / f"{name}-{view}.json"
    r["result"] = json.loads(got.read_text(encoding="utf-8")) if got.is_file() else None
    return r


def run_entry(source: str, out: Path, views=VIEWS, info: dict | None = None) -> dict:
    """Draw and measure one entry (a sampled id, or a local file) in every view that applies; judged by criteria.
    info: the entry's line of sample.json (stratum, bin, pair ...); refusal: a CA-only or P-only id (A1c), drawn in
    the default view only."""
    info = dict(info or {})
    name = entry_name(source)
    work = out / "work" / name
    work.mkdir(parents=True, exist_ok=True)
    keep = ("stratum", "bin", "order", "quota", "pair", "pair_type", "role", "partner", "refusal")
    e: dict = {"id": info.get("id", name), "name": name, **{k: info[k] for k in keep if k in info}, "views": {}}
    try:
        path = resolve(source)
        e["file"], e["sha256"] = str(path), _sha256(path.read_bytes())
        e["facts"] = features(path)
    except Exception as err:  # a download that fails fails the entry, not the campaign
        e["error"] = f"{type(err).__name__}: {err}"
        e["criteria"] = criteria(e)
        return _verdict(e)
    e["facts"]["resolution"] = _resolution(path)
    wanted = ("default",) if e.get("refusal") else tuple(views)

    def run(view: str, focus: str | None = None) -> None:
        if view in wanted:
            e["views"][view] = run_view(view, path, name, work, focus)

    run("default")
    d = e["views"].get("default", {}).get("result") or {}
    f = e["facts"]
    if not e.get("refusal"):
        if f["protein"]:
            run("deposited")
            run("beta")
        if f["chains"] >= 2:  # the whole file when the default drew part, else the first protein chain alone
            scoped = (d.get("scope") or {}).get("reason")
            run("focus", "none" if scoped or not d.get("first_protein_chain") else d["first_protein_chain"])
        if len(d.get("protein_chains_drawn", [])) >= 2:
            run("domains")
        if e.get("stratum") == "S6":
            run("membrane")
        for view in ("seed1", "roundtrip", "interactive"):
            run(view)
        if f["protein"]:
            run("sequence")
        run("app")
        if f["protein"] and "default" in wanted:
            e["dssp"] = dssp_agreement(path)
    e["criteria"] = criteria(e)
    return _verdict(e)


def _resolution(path: Path) -> float | None:
    import gemmi

    block = gemmi.cif.read(str(path)).sole_block()
    for tag in ("_refine.ls_d_res_high", "_em_3d_reconstruction.resolution", "_reflns.d_resolution_high"):
        v = block.find_value(tag)
        try:
            return round(float(gemmi.cif.as_string(v)), 2) if v else None
        except ValueError:
            continue
    return None


def _verdict(e: dict) -> dict:
    e["failed"] = [f"{k}: {c['detail']}" for k, c in e["criteria"].items() if c["pass"] is False]
    e["pass"] = not e["failed"]
    return e


# ------------------------------------------------------------------------------------------------ criteria


def _ok(passed: bool | None, detail: str = "") -> dict:
    return {"pass": passed, "detail": detail}


def _views(e: dict, names) -> dict[str, dict]:
    """The measured results of these views, those run."""
    return {v: (e["views"][v].get("result") or {}) for v in names if v in e["views"]}


def _count(results: dict[str, dict], key: str, limit: float = 0) -> dict:
    """A check counted in each view's figcheck results: pass when every view's count is at most limit; a count
    that could not be measured (None) fails."""
    bad = []
    for v, r in results.items():
        got = (r.get("checks") or {}).get(key, None) if "checks" in r else None
        if got is None:
            why = ((r.get("checks") or {}).get("skipped") or {}).get(key) or r.get("measure_error", "not measured")
            bad.append(f"{v} not measured ({str(why).strip().splitlines()[-1][:120]})")
        elif got > limit:
            details = ((r["checks"].get("details") or {}).get(key)) or []
            bad.append(f"{v} {got}" + (f" ({'; '.join(map(str, details[:3]))})" if isinstance(details, list) else ""))
    return _ok(None if not results else not bad, "; ".join(bad))


def criteria(e: dict) -> dict:
    """A1-A20 of one entry (plan 10, "Pass criteria"): {"pass": True, False or None (does not apply), "detail"}."""
    out = {k: _ok(None) for k in CRITERIA}
    if e.get("error"):
        out["A1"] = _ok(False, f"cannot read the entry: {e['error']}")
        return out
    views = e["views"]
    if e.get("refusal"):
        d = views.get("default", {})
        said = d.get("foldmap_said", [])
        ok = d.get("exit", 0) != 0 and not d.get("traceback") and len(said) == 1
        out["A1c"] = _ok(ok, "" if ok else f"exit {d.get('exit')}, said {said[:2]}, traceback {d.get('traceback')}")
        return out
    problems = []
    for v, r in views.items():
        res = r.get("result") or {}
        if r["timed_out"]:
            problems.append(f"{v}: killed after {LIMITS['view_timeout_s']} s")
        elif r["exit"] != 0:
            problems.append(f"{v}: exit {r['exit']} ({' | '.join(r['stderr_tail'][-2:])})")
        if r["traceback"]:
            problems.append(f"{v}: traceback")
        warns = res.get("warnings", []) + r.get("python_warnings", [])
        bad = [w for w in warns if not any(re.search(p, w) for p in DOCUMENTED_WARNINGS)]
        if bad:
            problems.append(f"{v}: {bad[0]}")
        empty = [n for n, size in r.get("written", {}).items() if not size]
        if empty and r["exit"] == 0:
            problems.append(f"{v}: not written: {', '.join(empty)}")
        if r["exit"] == 0 and r.get("result") is None:
            problems.append(f"{v}: no result")
        if res.get("measure_error"):
            problems.append(f"{v}: not measured: {res['measure_error'].strip().splitlines()[-1][:160]}")
    out["A1"] = _ok(not problems, "; ".join(problems))
    d = _views(e, ["default"]).get("default", {})
    if "default" in views:
        many = len(d.get("chains_read", [])) > LIMITS["many_chains"]
        limit = LIMITS["default_s_many_chains" if many else "default_s"]
        took = views["default"]["seconds"]
        out["A2"] = _ok(took <= limit, f"{took} s (limit {limit} s)")
    app = (_views(e, ["app"]).get("app") or {}).get("app")
    if "app" in views:
        if app is None:
            out["A2b"] = _ok(False, "the app view gave no answer")
        elif app["refused"] is None:
            out["A2b"] = _ok(True)
        else:
            fine = app["refused"] == "size" and app["residues"] > LIMITS["app_max_residues"]
            out["A2b"] = _ok(fine, "" if fine else f"{app['refused']}: {app['message']}")
    figures = _views(e, FIGURE_VIEWS)
    stuck = {v: r.get("stuck", []) for v, r in figures.items() if r.get("stuck")}
    over = {v: r.get("over_labels", []) for v, r in figures.items() if r.get("over_labels")}
    if figures:
        out["A3"] = _ok(not stuck, "; ".join(f"{v}: {', '.join(s[:4])}" for v, s in stuck.items()))
        out["A4"] = _ok(not over, "; ".join(f"{v}: {', '.join(s[:4])}" for v, s in over.items()))
    out["A5"] = _count(figures, "element_overlaps")
    out["A6"] = _count(figures, "text_collisions")
    out["A7"] = _count(figures, "off_page")
    out["A8"] = _count(_views(e, LABEL_VIEWS), "label_repeats")
    out["A9"] = _count(_views(e, ["default", "domains"]), "scope_silent")
    dssp = e.get("dssp")
    if dssp is not None:
        if dssp["agree"] is None:
            out["A10"] = _ok(True, f"unmeasured: {dssp['unmeasured']}")
        else:
            out["A10"] = _ok(dssp["agree"] >= LIMITS["dssp_agree"], f"{dssp['agree']:.4f} of {dssp['residues']}")
    out["A11"] = _count(_views(e, ["default"]), "order")
    for key, other, fields in (
        ("A12", "seed1", ("signature_sha256", "loops_sha256", "svg_sha256")),
        ("A13", "roundtrip", ("signature_sha256", "svg_sha256")),
    ):
        r = _views(e, [other]).get(other)
        if r is not None and "default" in views:
            differ = [k.split("_")[0] for k in fields if d.get(k) is None or r.get(k) != d.get(k)]
            out[key] = _ok(not differ, f"{other} differs from default: {', '.join(differ)}" if differ else "")
    extra = {v: views[v] for v in ("interactive", "sequence") if v in views}
    if extra:
        bad = [v for v, r in extra.items() if r["exit"] != 0 or not all(r.get("written", {}).values())]
        out["A14"] = _ok(not bad, ", ".join(bad))
    out["A15"] = _count(_views(e, ["default", "focus"]), "copies_alike", LIMITS["copies_alike"])
    dom = _views(e, ["domains"])
    if dom:
        parts = [
            _count(dom, "domain_rules"),
            _count(dom, "domain_overlap", LIMITS["domain_overlap"]),
            _count(dom, "domain_name_collisions"),
        ]
        out["A16"] = _ok(all(p["pass"] for p in parts), "; ".join(p["detail"] for p in parts if p["detail"]))
    if "na" in d:
        bad = list(d["na"].get("failed", []))
        gaps = ((d.get("checks") or {}).get("details") or {}).get("stem_name_gaps")
        limit = d.get("stem_gap_limit")
        if gaps is None or (gaps and limit is None):
            bad.append("stem-name gaps not measured")
        else:
            for k, (own, other) in sorted(gaps.items()):
                if own > limit or own >= other:
                    bad.append(f"stem name {k} {own:.3f} from its letters, {other:.3f} from another stem")
        out["A17"] = _ok(not bad, "; ".join(bad))
    mem = _views(e, ["membrane"]).get("membrane")
    if mem is not None:
        m = mem.get("membrane")
        if m is None:
            out["A18"] = _ok(False, "not measured")
        else:
            bad = ([] if m["band"] else ["no band drawn"]) + [f"TM helix {x}" for x in m["tm_short"]]
            bad += [f"loop in the slab {x}" for x in m["loops_in_slab"]]
            out["A18"] = _ok(not bad, "; ".join(bad))
    if e.get("pair_check") is not None:
        out["A19"] = e["pair_check"]
    if "default" in views and d:
        long = [f"{g} {n}" for g, n in d.get("tethers", []) if n > LIMITS["tether_length"]]
        out["A20"] = _ok(not long, "; ".join(long[:4]))
    return out


def pair_check(a: dict, b: dict) -> dict:
    """A19: the beta views of the two entries of a pair, the first protein chain of each. Elements pair one-to-one
    by author residue range as copies of one protein do (model.pair_places: the same kind, sharing at least half of
    the shorter, the largest overlaps first). The pairs must keep one order, and pair at least half of the larger
    chain's elements (and 3, when both have that many). One entry may hold an element the other lacks (listed, not
    failed); the generic names after it then differ, so names are compared only up to the first unpaired element."""
    from foldmap.model import pair_places

    ra = (a["views"].get("beta", {}).get("result") or {}).get("elements")
    rb = (b["views"].get("beta", {}).get("result") or {}).get("elements")
    if not ra or not rb:
        return _ok(False, "beta view not measured for both")

    def first_chain(els):
        chain = els[0]["chain"]
        return sorted((x for x in els if x["chain"] == chain), key=lambda x: (x["first"], x["last"]))

    ea, eb = first_chain(ra), first_chain(rb)
    found = pair_places(*([(x["kind"], x["first"], x["last"]) for x in els] for els in (ea, eb)))
    pairs = sorted(found.items())
    bad, notes = [], [f"{x['id']} unpaired" for i, x in enumerate(ea) if i not in found]
    notes += [f"{y['id']} unpaired" for k, y in enumerate(eb) if k not in found.values()]
    if [k for _, k in pairs] != sorted(k for _, k in pairs):
        bad.append("elements in another order")
    if len(pairs) < max(min(3, len(ea), len(eb)), max(len(ea), len(eb)) / 2):
        bad.append(f"only {len(pairs)} of {len(ea)} and {len(eb)} elements paired")
    for i, k in pairs:  # names, while every element before is paired (generic numbering shifts after an extra one)
        if any(j not in found for j in range(i)) or any(j not in found.values() for j in range(k)):
            break
        if ea[i]["label"] != eb[k]["label"]:
            bad.append(f"{ea[i]['label']}~{eb[k]['label']} ({ea[i]['id']}~{eb[k]['id']})")
    shown = bad + notes
    detail = "; ".join(shown[:8]) + (f" (+{len(shown) - 8} more)" if len(shown) > 8 else "")
    return _ok(not bad, detail)


# ------------------------------------------------------------------------------------------- contact sheet


def _project(xyz, axes, centre, side: bool) -> np.ndarray:
    """Page (x, y) of 3D points seen through the figure's page axes; side: turned 90° about page-up (what faced
    the viewer now points left)."""
    q = (np.asarray(xyz, float).reshape(-1, 3) - np.asarray(centre, float)) @ np.asarray(axes, float).T
    return np.column_stack([-q[:, 2], q[:, 1]]) if side else q[:, :2]


def _traces(ax, trace: dict, side: bool) -> list[np.ndarray]:
    """Draw every chain's CA trace and strand's P trace, residues in their element's or stem's colour (grey
    outside), chain ends marked; returns the points drawn."""
    from matplotlib.collections import LineCollection

    drawn = []
    axes, centre = trace["axes"], trace["centre"]
    for kind, items in (("protein", trace["chains"]), ("na", trace["strands"])):
        for c in items:
            q = _project(c["xyz"], axes, centre, side)
            if not len(q):
                continue
            drawn.append(q)
            if kind == "protein":
                joined = c["joined"][1:]
            else:  # a strand's P atoms more than 10 Å apart are not one step
                joined = list(np.linalg.norm(np.diff(np.asarray(c["xyz"], float), axis=0), axis=1) < 10.0)
            segs, colours, widths = [], [], []
            for k in range(len(q) - 1):
                if not joined[k]:
                    continue
                a, b = c["colour"][k], c["colour"][k + 1]
                coloured = a is not None and a == b
                segs.append(q[k : k + 2])
                colours.append(a if coloured else ("#555555" if kind == "na" else "#b0b0b0"))
                widths.append(1.8 if coloured else 0.7)
            if segs:
                ax.add_collection(LineCollection(segs, colors=colours, linewidths=widths, capstyle="round"))
            first, last = ("N", "C") if kind == "protein" else ("5'", "3'")
            for (x, y), text in ((q[0], f"{c['chain']} {first}"), (q[-1], f"{c['chain']} {last}")):
                ax.plot([x], [y], "o", ms=3, color="k")
                ax.annotate(text, (x, y), fontsize=6, xytext=(2, 2), textcoords="offset points")
    return drawn


def _facts(e: dict) -> list[str]:
    f, d = e.get("facts") or {}, (e.get("views", {}).get("default", {}).get("result") or {})
    scope = d.get("scope") or {}
    said = [s for s in e.get("views", {}).get("default", {}).get("foldmap_said", []) if "warning" in s]
    lines = [
        f"{e['id']}  {e.get('stratum', '')} {e.get('bin', '')}"
        + (f"  pair {e['pair']} {e['pair_type']} ({e['role']}, partner {e['partner']})" if e.get("pair") else ""),
        f"method {f.get('method')}{' (' + f['model_name'] + ')' if f.get('model_name') else ''}; "
        f"resolution {f.get('resolution')}; models {f.get('models')}; "
        f"{f.get('chains')} polymer chains, {f.get('deposited_residues')} residues deposited",
        f"read {len(d.get('chains_read', []))} chains; drew protein {', '.join(d.get('protein_chains_drawn', [])[:12])}"
        f" NA {', '.join(d.get('na_chains_drawn', [])[:12])}; assembly 1: {f.get('assembly1_chains')} chains",
        f"symmetry: {d.get('context') or ('copies tied: ' + str(d.get('copy_ties', 0)))}; "
        f"scope: {scope.get('reason') or 'everything drawn'}",
        f"ligands {', '.join(f.get('ligands', [])[:10]) or '-'}; glycans {', '.join(f.get('glycans', [])[:6]) or '-'}"
        f"{' (N-linked)' if f.get('n_linked') else ''}; modified {', '.join(f.get('modified', [])[:10]) or '-'}",
    ]
    return lines + [s[:200] for s in said[:3]]


def _metric_items(e: dict) -> list[str]:
    items = []
    for v, r in e.get("views", {}).items():
        items.append(f"{v}={r['seconds']}s/{r['max_rss_mb']:.0f}MB")
    d = e.get("views", {}).get("default", {}).get("result") or {}
    ch = d.get("checks") or {}
    for k in ("honest_crossings", "loop_crossings", "copies_alike", "stem_name_gap"):
        if ch.get(k) is not None:
            items.append(f"{k}={ch[k]}")
    if d.get("area") is not None:
        items.append(f"area={d['area']}")
    for k, v in (d.get("route") or {}).items():
        if v is not None and k != "records":
            items.append(f"route.{k}={v}")
    if e.get("dssp"):
        items.append(f"dssp={e['dssp'].get('agree')}")
    for k, v in (d.get("na") or {}).items():
        if k not in ("failed", "note", "missed_pairs", "extra_pairs") and v is not None:
            items.append(f"{k}={v}")
    return items


def _thumbnail(path: Path, most: int = 2400) -> np.ndarray:
    """A PNG as an RGB(A) array, at most `most` pixels on its longer side (a 300 dpi figure of a large entry would
    not fit in memory as floats)."""
    from PIL import Image

    with Image.open(path) as img:
        img.thumbnail((most, most))
        return np.asarray(img)


def sheet(e: dict, out: Path, partner: dict | None = None) -> dict:
    """The contact sheet <out>/<name>.png, 110 dpi: the default figure | the 3D traces seen through its page axes,
    in its elements' colours | the same turned 90° about page-up | (pairs) the partner's figure; facts and every
    criterion below, failed ones in red. Returns what it drew."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nav = _script_module("na_validate")
    name = e["name"]
    work = out / "work" / name
    figure_png = work / f"{name}-default.png"
    trace_path = work / f"{name}-trace.json"
    trace = json.loads(trace_path.read_text(encoding="utf-8")) if trace_path.is_file() else None
    images = [_thumbnail(figure_png)] if figure_png.is_file() else [None]
    titles = [f"{e['id']}: foldmap (default)"]
    if trace is not None and trace.get("axes") and any(c["xyz"] for c in trace["chains"] + trace["strands"]):
        seen = [c["xyz"] for c in trace["chains"] + trace["strands"] if c["xyz"]]
        boxes = [nav._bounds([_project(xyz, trace["axes"], trace["centre"], s) for xyz in seen]) for s in (False, True)]
    else:
        boxes = []
    pw = out / "work" / partner["name"] / f"{partner['name']}-default.png" if partner else None
    panel_h, gap = 6.0, 0.25
    aspects = [min(max(img.shape[1] / img.shape[0], 0.4), 2.5) if img is not None else 1.0 for img in images]
    aspects += [min(max((b[1] - b[0]) / max(b[3] - b[2], 1e-9), 0.4), 2.5) for b in boxes]
    titles += ["3D, the figure's page axes (CA and P; element colours)", "3D turned 90° about page-up"][: len(boxes)]
    if pw is not None and pw.is_file():
        partner_img = _thumbnail(pw)
        aspects.append(min(max(partner_img.shape[1] / partner_img.shape[0], 0.4), 2.5))
        titles.append(f"partner {partner['id']} ({e.get('pair_type')})")
    else:
        partner_img = None
    widths = [panel_h * a for a in aspects]
    total_w = max(sum(widths) + gap * (len(widths) + 1), 12.0)
    size = 8
    chars = int((total_w - 0.4) / (size * 0.62 / 72))
    facts = _facts(e)
    crit = e.get("criteria") or {}
    marks = "  ".join(f"{k} {'-' if c['pass'] is None else 'ok' if c['pass'] else 'FAIL'}" for k, c in crit.items())
    failed = [f"{k}: {c['detail']}"[: chars * 2] for k, c in crit.items() if c["pass"] is False]
    lines = facts + nav.wrap_metrics(marks.split("  "), chars) + nav.wrap_metrics(_metric_items(e), chars)
    text_h = 0.4 + (len(lines) + sum(1 + len(f) // chars for f in failed)) * size * 1.3 / 72
    top = 0.5
    high = panel_h + text_h + top
    fig = plt.figure(figsize=(total_w, high))
    axes, x = [], gap
    for w in widths:
        axes.append(fig.add_axes([x / total_w, text_h / high, w / total_w, panel_h / high]))
        x += w + gap
    for a, t in zip(axes, titles):
        a.set_title(nav.wrap_title(t, a.get_position().width * total_w, 9), loc="left", fontsize=9)
        a.set_axis_off()
        a.set_anchor("N")
    if images[0] is not None:
        axes[0].imshow(images[0])
    else:
        axes[0].text(0.5, 0.5, "no figure", ha="center", va="center", color="red", transform=axes[0].transAxes)
    drawn = {"panels": len(axes), "protein_traces": 0, "na_traces": 0}
    for k, side in enumerate((False, True)[: len(boxes)]):
        a = axes[1 + k]
        _traces(a, trace, side)
        a.set_xlim(*boxes[k][:2])
        a.set_ylim(*boxes[k][2:])
        a.set_aspect("equal", adjustable="box")
    if boxes:
        drawn["protein_traces"] = len(trace["chains"])
        drawn["na_traces"] = len(trace["strands"])
    if partner_img is not None:
        axes[-1].imshow(partner_img)
    y = text_h - 0.2
    for line in failed:
        for part in [line[i : i + chars] for i in range(0, len(line), chars)]:
            fig.text(0.15 / total_w, y / high, part, fontsize=size, color="red", family="monospace", va="top")
            y -= size * 1.3 / 72
    for line in lines:
        fig.text(0.15 / total_w, y / high, line[: chars + 20], fontsize=size, family="monospace", va="top")
        y -= size * 1.3 / 72
    path = out / f"{name}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    drawn["path"] = str(path)
    return drawn


# ------------------------------------------------------------------------------------------------ report


def _versions() -> dict:
    import gemmi
    import matplotlib

    out = {"gemmi": gemmi.__version__, "matplotlib": matplotlib.__version__, "numpy": np.__version__}
    rows = mkdssp_version()
    out["mkdssp"] = rows
    py = os.environ.get("NA_RNAPOLIS_PYTHON", "")
    out["rnapolis"] = None
    if py:
        try:
            code = "import importlib.metadata as m; print(m.version('rnapolis'))"
            out["rnapolis"] = subprocess.run(
                [py, "-I", "-c", code], capture_output=True, text=True, timeout=60
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    return out


def mkdssp_version() -> str | None:
    exe = os.environ.get("MKDSSP") or shutil.which("mkdssp")
    if not exe:
        return None
    try:
        return (
            subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60)
            .stdout.strip()
            .splitlines()[0]
        )
    except (OSError, subprocess.SubprocessError, IndexError):
        return None


def report(entries: list[dict], out: Path, sample_path: Path | None = None, total: int | None = None) -> dict:
    """report.json (every value, LIMITS, commit, versions, sample.json's SHA-256, date), report.md (the counts,
    a table per stratum and per entry with the A1-A20 marks and what failed, the reported-only statistics, the
    limits) and index.html (a grid of the contact sheets, linking to each)."""
    total = total or len(entries)
    passed = sum(1 for e in entries if e.get("pass"))
    doc = {
        "date": datetime.date.today().isoformat(),
        "foldmap_commit": _commit(),
        "versions": _versions(),
        "sample_sha256": _sha256(sample_path.read_bytes()) if sample_path and sample_path.is_file() else None,
        "limits": LIMITS,
        "documented_warnings": list(DOCUMENTED_WARNINGS),
        "criteria": CRITERIA,
        "passed": passed,
        "entries": entries,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(_jsonable(doc), indent=1) + "\n", encoding="utf-8")
    keys = list(CRITERIA)
    strata: dict[str, list] = {}
    for e in entries:
        strata.setdefault(e.get("stratum") or "-", []).append(e)
    lines = [
        "# foldmap random-100 campaign",
        "",
        f"{passed} of {total} pass.",
        "",
        "| stratum | entries | pass |",
        "|---|---|---|",
    ]
    lines += [f"| {s} | {len(es)} | {sum(1 for e in es if e.get('pass'))} |" for s, es in strata.items()]
    lines += ["", "| id | stratum | bin | chains | seconds | " + " | ".join(keys) + " | failed |"]
    lines.append("|" + "---|" * (len(keys) + 6))
    for e in entries:
        crit = e.get("criteria") or {}
        marks = ["-" if crit.get(k, {}).get("pass") is None else "✓" if crit[k]["pass"] else "✗" for k in keys]
        took = (e.get("views", {}).get("default") or {}).get("seconds", "")
        detail = "; ".join(e.get("failed", [])).replace("|", "/")[:600]
        chains = (e.get("facts") or {}).get("chains", "")
        lines.append(
            f"| {e['id']} | {e.get('stratum', '')} | {e.get('bin', '')} | {chains} | {took} | "
            + " | ".join(marks)
            + f" | {detail} |"
        )
    lines += [
        "",
        "## Reported, not gated",
        "",
        "| id | " + " | ".join(REPORTED) + " |",
        "|" + "---|" * (len(REPORTED) + 1),
    ]
    for e in entries:
        lines.append(f"| {e['id']} | " + " | ".join(str(v) for v in reported(e).values()) + " |")
    lines += ["", "## Criteria", ""] + [f"- **{k}** {v}" for k, v in CRITERIA.items()]
    lines += ["", "## Limits", ""] + [f"- {k}: {v}" for k, v in LIMITS.items()]
    lines += [f"- documented warnings: {', '.join(DOCUMENTED_WARNINGS) or 'none'}", ""]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    cells = "\n".join(
        f'<a href="{e["name"]}.png" title="{e["id"]}: {"pass" if e.get("pass") else "FAIL"}">'
        f'<img src="{e["name"]}.png" alt="{e["id"]}" class="{"pass" if e.get("pass") else "fail"}"></a>'
        for e in entries
    )
    (out / "index.html").write_text(
        "<!doctype html><meta charset=utf-8><title>random-100</title><style>"
        "body{font-family:sans-serif}div{display:grid;grid-template-columns:repeat(10,1fr);gap:4px}"
        "img{width:100%;border:3px solid #2a2}img.fail{border-color:#c22}</style>"
        f"<h1>foldmap random-100: {passed} of {total} pass</h1><div>{cells}</div>\n",
        encoding="utf-8",
    )
    return doc


REPORTED = (
    "honest_crossings",
    "loop_crossings",
    "na_run_crossings",
    "fallback_curves",
    "fallback_curves_crossing_an_element",
    "kink_links_rejected",
    "tighter_than_inflate",
    "rejected_as_crossing",
    "area",
    "max_rss_mb",
)


def reported(e: dict) -> dict:
    """The statistics reported, not gated (default view; max RSS over every view)."""
    d = (e.get("views", {}).get("default") or {}).get("result") or {}
    ch, route = d.get("checks") or {}, d.get("route") or {}
    return {
        "honest_crossings": ch.get("honest_crossings"),
        "loop_crossings": ch.get("loop_crossings"),
        "na_run_crossings": (d.get("na") or {}).get("na_run_crossings"),
        "fallback_curves": route.get("fallback_curves"),
        "fallback_curves_crossing_an_element": route.get("fallback_curves_crossing_an_element"),
        "kink_links_rejected": route.get("kink_links_rejected"),
        "tighter_than_inflate": route.get("tighter_than_inflate"),
        "rejected_as_crossing": route.get("rejected_as_crossing"),
        "area": d.get("area"),
        "max_rss_mb": max((r["max_rss_mb"] for r in e.get("views", {}).values()), default=None),
    }


# ------------------------------------------------------------------------------------------------ commands


def _sample_entries(sample_path: Path) -> list[dict]:
    """sample.json's entries, then its CA-only / P-only replacements (A1c: drawn once, expected to be refused)."""
    doc = json.loads(sample_path.read_text(encoding="utf-8"))
    refusals = [{"id": r["id"], "stratum": "A1c", "refusal": True} for r in doc["replaced"] if r["replacement"] is None]
    return doc["entries"] + refusals


def cmd_entry(args) -> int:
    out = Path(args.out)
    info = {"id": entry_name(args.id)}
    sample_path = Path(args.sample)
    if sample_path.is_file():
        info = next((x for x in _sample_entries(sample_path) if x["id"] == args.id), info)
    e = run_entry(args.id, out, _view_list(args.views), info)
    (out / f"{e['name']}.json").write_text(json.dumps(_jsonable(e), indent=1) + "\n", encoding="utf-8")
    print(f"{e['id']}: {'PASS' if e['pass'] else 'FAIL ' + '; '.join(e['failed'])[:400]}", flush=True)
    return 0 if e["pass"] else 1


def _view_list(text: str) -> tuple[str, ...]:
    if text == "quick":
        return QUICK
    if text in ("all", "default", ""):
        return VIEWS if text != "default" else ("default",)
    views = tuple(v.strip() for v in text.split(","))
    unknown = [v for v in views if v not in VIEWS]
    if unknown:
        raise SystemExit(f"unknown view(s): {', '.join(unknown)}; choose from {', '.join(VIEWS)}, quick or all")
    return views


def _load_entry(out: Path, x: dict) -> dict | None:
    p = out / f"{entry_name(x['id'])}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None


def finish(out: Path, sample_path: Path, ids: list[str] | None = None) -> dict:
    """A19 for the pairs, the contact sheets (of `ids`, and their partners; None: all) and the report, from the
    entries' JSON files in out."""
    listed = _sample_entries(sample_path) if sample_path.is_file() else [{"id": i} for i in ids or []]
    entries = {x["id"]: e for x in listed if (e := _load_entry(out, x)) is not None}
    for e in entries.values():
        partner = entries.get(e.get("partner"))
        if e.get("pair") and partner is not None and e.get("criteria"):
            e["pair_check"] = pair_check(e, partner)
            e["criteria"] = criteria(e)
            _verdict(e)
            (out / f"{e['name']}.json").write_text(json.dumps(_jsonable(e), indent=1) + "\n", encoding="utf-8")
    redo = set(ids) if ids is not None else set(entries)
    redo |= {e["partner"] for e in entries.values() if e.get("partner") and e["id"] in redo}
    for i in sorted(redo):
        if i in entries:
            sheet(entries[i], out, entries.get(entries[i].get("partner")))
    return report(list(entries.values()), out, sample_path, total=len(listed))


def cmd_run(args) -> int:
    from concurrent.futures import ThreadPoolExecutor

    out, sample_path = Path(args.out), Path(args.sample)
    listed = _sample_entries(sample_path) if sample_path.is_file() else []
    if args.ids:
        ids = args.ids.split(",")
    elif args.failed_from:
        old = json.loads(Path(args.failed_from).read_text(encoding="utf-8"))
        ids = [e["id"] for e in old["entries"] if not e.get("pass")]
    else:
        ids = [x["id"] for x in listed]
    views = _view_list(args.views)
    out.mkdir(parents=True, exist_ok=True)
    log = out / "run.log"

    def one(i: str) -> str:
        cmd = [sys.executable, "-I", str(SCRIPT), "entry", i, "--out", str(out), "--sample", str(sample_path)]
        cmd += ["--views", args.views]
        env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "PYTHONHASHSEED": "0"}
        work = out / "work" / entry_name(i)
        work.mkdir(parents=True, exist_ok=True)
        timeout = LIMITS["view_timeout_s"] * len(views) + 600
        r = run_child(cmd, env=env, timeout=timeout, out=work / "entry.out", err=work / "entry.err")
        e = _load_entry(out, {"id": i})
        if e is None or r["timed_out"]:  # the entry process died before writing: a failure, said so
            tail = (work / "entry.err").read_text(encoding="utf-8", errors="replace").strip().splitlines()[-3:]
            e = {"id": i, "name": entry_name(i), "views": {}, "error": f"entry process: exit {r['exit']}"}
            e["error"] += f", killed after {timeout} s" if r["timed_out"] else f" ({' | '.join(tail)})"
            info = next((x for x in listed if x["id"] == i), {})
            e.update({k: info[k] for k in ("stratum", "bin", "pair", "pair_type", "role", "partner") if k in info})
            e["criteria"] = criteria(e)
            _verdict(e)
        e["process"] = r
        (out / f"{e['name']}.json").write_text(json.dumps(_jsonable(e), indent=1) + "\n", encoding="utf-8")
        line = f"{i}: {'PASS' if e['pass'] else 'FAIL ' + '; '.join(e['failed'])[:300]} ({r['seconds']} s)"
        with open(log, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
        print(line, flush=True)
        return line

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        list(pool.map(one, ids))
    doc = finish(out, sample_path, ids if args.ids or args.failed_from else None)
    print(f"{doc['passed']} of {len(listed) or len(ids)} pass; report in {out / 'report.md'}")
    return 0 if doc["passed"] == (len(listed) or len(ids)) else 1


def cmd_sheet(args) -> int:
    out = Path(args.out)
    entries = {i: _load_entry(out, {"id": i}) for i in args.ids}
    for i, e in entries.items():
        if e is None:
            print(f"{i}: no {entry_name(i)}.json in {out}", file=sys.stderr)
            return 1
        partner = _load_entry(out, {"id": e["partner"]}) if e.get("partner") else None
        print(sheet(e, out, partner)["path"])
    return 0


def cmd_report(args) -> int:
    doc = finish(Path(args.out), Path(args.sample), [])
    print(f"{doc['passed']} pass; report in {Path(args.out) / 'report.md'}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("sample", help="draw the stratified sample (network)")
    s.add_argument("--seed", type=int, default=SEED)
    s.add_argument("--date", default=datetime.date.today().isoformat(), help="release-date cutoff (YYYY-MM-DD)")
    s.add_argument("--out", default=str(OUT / "sample.json"))
    s.add_argument("--check", action="store_true", help="list the sampled ids now obsolete; change nothing")
    s.set_defaults(func=cmd_sample)
    views = "all, quick (default, seed1, roundtrip) or a comma-separated list of " + ", ".join(VIEWS)
    r = sub.add_parser("run", help="draw, measure and report every sampled entry (network for downloads)")
    r.add_argument("--ids", help="comma-separated ids (default: every entry of the sample)")
    r.add_argument("--failed-from", metavar="REPORT", help="the entries a report.json says failed")
    r.add_argument("--views", default="all", help=views)
    r.add_argument("--out", default=str(OUT))
    r.add_argument("--sample", default=str(OUT / "sample.json"))
    r.add_argument("--jobs", type=int, default=1, help="entries drawn at once")
    r.set_defaults(func=cmd_run)
    e = sub.add_parser("entry", help="draw and measure one entry in every view that applies (run's child)")
    e.add_argument("id", help="a sampled id, or a local structure file")
    e.add_argument("--out", default=str(OUT))
    e.add_argument("--views", default="all", help=views)
    e.add_argument("--sample", default=str(OUT / "sample.json"))
    e.set_defaults(func=cmd_entry)
    v = sub.add_parser("view", help="one view of one entry, measured (entry's child)")
    v.add_argument("view", choices=VIEWS)
    v.add_argument("path")
    v.add_argument("--name", required=True)
    v.add_argument("--work", required=True)
    v.add_argument("argv", nargs=argparse.REMAINDER, help="-- then the foldmap command line")
    v.set_defaults(func=cmd_view)
    c = sub.add_parser("sheet", help="draw the contact sheets of measured entries again")
    c.add_argument("ids", nargs="+")
    c.add_argument("--out", default=str(OUT))
    c.set_defaults(func=cmd_sheet)
    p = sub.add_parser("report", help="write report.md, report.json and index.html from the measured entries")
    p.add_argument("--out", default=str(OUT))
    p.add_argument("--sample", default=str(OUT / "sample.json"))
    p.set_defaults(func=cmd_report)
    args = ap.parse_args(argv)
    if getattr(args, "argv", None) and args.argv[0] == "--":
        args.argv = args.argv[1:]
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
