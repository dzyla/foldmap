# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The random-100 campaign: a stratified random sample of 100 PDB entries and AlphaFold models, each drawn in every
view foldmap offers, measured and shown on a contact sheet (plan 10, needs network; not run in CI).

    python scripts/random100.py sample [--seed 20261013] [--date YYYY-MM-DD] [--out .../sample.json] [--check]
                                       [--exclude other/sample.json]

sample draws the sample from the RCSB search API (https://search.rcsb.org/rcsbsearch/v2/query) and writes
docs/validation/random100/sample.json. Every pool is one stratum x size-bin query, all hits, sorted; it is drawn
from in the order random.Random(f"{seed}:{stratum}:{bin}") gives it (a string seed is hashed with SHA-512, so
PYTHONHASHSEED cannot change the draw), at most LOOK candidates deep. A candidate is downloaded (foldmap.fetch, into
~/.cache/foldmap) and read with gemmi; it is skipped, with the reason logged, when it is obsolete or cannot be
downloaded, holds no protein or nucleic-acid polymer, is CA-only or P-only, or is in the sample already. The extra
quotas inside a stratum (S4 helical, big, cap; S8 ribosome) are drawn first from their own sub-pools; the pairs of S11
follow PAIR_TYPES; the cross-cutting quotas (QUOTAS) are then topped up. --exclude takes another sample's ids out of
every pool (the held-out sample: final review A, Important 2). The same seed, date, pools and files give the same
draw; the file also records the foldmap commit, the gemmi version and the SHA-256 of every candidate file read, so a
changed file is seen. --check re-queries the ids and lists those now obsolete, without changing the file.

    python scripts/random100.py run [--ids A,B | --failed-from report.json] [--views all|quick|v1,v2]
                                    [--out docs/validation/random100] [--jobs 1]

run draws every entry (plus each CA-only / P-only id met while sampling, which must be refused: A1c) with `entry <id>
--out DIR`, a child process (python -I, OPENBLAS_NUM_THREADS=1, PYTHONHASHSEED=0). The entry draws each view (VIEWS) in
a grandchild of its own (`view`, run with -W default), killed after LIMITS["view_timeout_s"], so a crash or a hang in
one view does not hide the others; wall time and max RSS come from os.wait4 (where there is one: run_child). Each view
writes what it measured (figcheck.figure_checks, route statistics from foldmap.route's DEBUG summary record, layout
signature and SVG hashes, na_validate's nucleic-acid measures, membrane and tether checks); the entry adds mkdssp
agreement and judges A1-A20 (foldmap.criteria: criteria, LIMITS; `foldmap check` judges one figure by them too)
into <out>/<id>.json. Then the pairs are compared (A19), a contact sheet <out>/<id>.png is drawn per entry (the
figure, the 3D CA/P traces seen through the figure's page axes in its element colours, the same turned 90° about
page-up, the partner's figure for a pair, facts and criteria), and report.json, report.md and index.html are written.
Exit 0 only if every entry passes. Each view's files are kept in <out>/work/<id>/. Pairing against RNApolis needs
NA_RNAPOLIS_PYTHON (the Python of a venv that has it); mkdssp 3.1.2 is found on PATH or as MKDSSP.
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
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    # the foldmap of this checkout, not whichever one is installed (children run with python -I: no PYTHONPATH)
    sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

# A1-A20 live in the package (plan 10 Task 13a: `foldmap check` judges by them too), imported once the path is set
from foldmap import criteria as criteria_module  # noqa: E402
from foldmap.criteria import (  # noqa: E402
    CRITERIA,
    DOCUMENTED_WARNINGS,
    FIGURE_VIEWS,
    LIMITS,
    NOT_APPLICABLE,
    NOT_JUDGED,
    OUTCOMES,
    _ok,
    criteria,
    measure_figure,
    outcome,
    pair_check,
)
from foldmap.criteria import jsonable as _jsonable  # noqa: E402

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
    contacts: also whether a ligand touches the protein (within TOUCH Å). The additives foldmap hides from its
    figures (io.ADDITIVES: glycerol, sulfate ...) are no ligand here either (final review B, Important 4: 6YPX was
    "holo" by glycerol alone)."""
    import gemmi

    from foldmap.io import ADDITIVES

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
            elif len(res) >= LIGAND_ATOMS and res.name not in ADDITIVES:
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


def ca_trace(path) -> list[tuple[str, object]]:
    """(residue name, CA position) of each residue with a CA of the first protein chain of the first model, in order."""
    import gemmi

    _, st = _gemmi_structure(path)
    for chain in st[0]:
        span = chain.get_polymer()
        if len(span) and span.check_polymer_type() in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            return [(res.name, ca.pos) for res in span if (ca := res.find_atom("CA", "*")) is not None]
    return []


def ca_rmsd(path_a, path_b, least: int = 30) -> float | None:
    """CA RMSD (Å) of the first protein chains over their residues aligned by sequence (figcheck.aligned_positions,
    as A19 pairs them: one protein numbered differently in two entries, C18) after superposition; None under `least`
    aligned residues."""
    import gemmi

    from foldmap.figcheck import aligned_positions

    a, b = ca_trace(path_a), ca_trace(path_b)
    common = aligned_positions([n for n, _ in a], [n for n, _ in b])
    if len(common) < least:
        return None
    sup = gemmi.superpose_positions([a[i][1] for i, _ in common], [b[j][1] for _, j in common])
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
    """HEAD of this checkout, or FOLDMAP_COMMIT where there is no .git (a copy synced to another machine)."""
    if os.environ.get("FOLDMAP_COMMIT"):
        return os.environ["FOLDMAP_COMMIT"]
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

    def __init__(self, seed: int, date: str, exclude: Iterable[str] = ()):
        self.seed, self.date = seed, date
        self.exclude = set(exclude)  # another sample's ids, never drawn (the held-out sample)
        self.feats: dict[str, dict] = {}
        self.files: dict[str, str] = {}  # id -> SHA-256 of the file its features were read from
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
            path = resolve(i)
            f = self.feats[i] = features(path, contacts=contacts)
            self.files[i] = _sha256(path.read_bytes())
        return f

    def check(self, i: str, also: Callable[[dict], str] | None = None) -> tuple[bool, str]:
        """(usable, why not): in the sample already, cannot be downloaded (obsolete), no polymer, CA/P-only, or
        `also`'s reason."""
        if i in self.taken:
            return False, "already in the sample"
        if i in self.exclude:
            return False, "in the excluded sample"
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


def sample(seed: int, date: str, exclude: Iterable[str] = ()) -> dict:
    """Draw the whole sample (network): every stratum, the pairs, then the cross-cutting top-up; `exclude`: ids never
    drawn (another sample's)."""
    import gemmi

    s = Sampler(seed, date, exclude)
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
        if i in s.exclude:
            return False, "in the excluded sample"
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
    out = {
        "seed": seed,
        "query_date": date,
        "foldmap_commit": _commit(),
        "gemmi": gemmi.__version__,
        "look": LOOK,
        "bins": {b: [lo or 0, hi - 1] for b, (lo, hi) in BINS.items()},
        "strata": s.records,
        "entries": entries,
        "replaced": s.replaced,
        "skipped": s.skipped,
        "quotas": quotas,
    }
    if s.exclude:
        out["excluded"] = {"ids": len(s.exclude), "sha256": _sha256("\n".join(sorted(s.exclude)).encode())}
    out["files"] = dict(sorted(s.files.items()))  # every candidate file read, by id: its SHA-256
    return out


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
    exclude = [e["id"] for e in _sample_entries(Path(args.exclude))] if args.exclude else []
    doc = sample(args.seed, args.date, exclude)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(dump(doc), encoding="utf-8")
    print(quota_table(doc))
    short = [k for k, q in doc["quotas"].items() if q["reached"] < q["target"]]
    print(f"{len(doc['entries'])} entries written to {out}" + (f"; short: {', '.join(short)}" if short else ""))
    return 1 if short or len(doc["entries"]) != 100 else 0


# ------------------------------------------------------------------------------------------------ runner

SCRIPT = Path(__file__).resolve()
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
_WARNING = re.compile(r"^\S.*:\d+: (\w*Warning): (.*)$")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


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


def _rss_mb(ru_maxrss: float) -> float:
    """Max RSS in MB from getrusage's ru_maxrss (kilobytes on Linux, bytes on macOS)."""
    return round(ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)


def run_child(cmd: list[str], *, env: dict | None = None, timeout: float, out: Path, err: Path) -> dict:
    """Run cmd in its own session, stdout and stderr to files; killed (with everything it started) after timeout
    seconds. Exit code, wall seconds and the child's own max RSS (os.wait4: its rusage, the waited-for children it
    reaped included). Without os.wait4 (Windows) the child is waited for by subprocess and its max RSS is the
    largest of this process's reaped children (resource.RUSAGE_CHILDREN), or None where there is no resource."""
    if not hasattr(os, "wait4"):
        return _run_child_portable(cmd, env=env, timeout=timeout, out=out, err=err)
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
        "max_rss_mb": _rss_mb(usage.ru_maxrss),
        "timed_out": timed_out,
    }


def _run_child_portable(cmd: list[str], *, env: dict | None, timeout: float, out: Path, err: Path) -> dict:
    """run_child where os.wait4 is missing."""
    try:
        import resource
    except ImportError:
        resource = None
    t0 = time.perf_counter()
    with open(out, "wb") as fo, open(err, "wb") as fe:
        proc = subprocess.Popen(cmd, stdout=fo, stderr=fe, env=env, start_new_session=os.name == "posix")
        timed_out = False
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            if hasattr(os, "killpg"):
                os.killpg(proc.pid, signal.SIGKILL)
            else:
                proc.kill()
            proc.wait()
    rss = _rss_mb(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) if resource is not None else None
    return {
        "exit": proc.returncode,
        "seconds": round(time.perf_counter() - t0, 2),
        "max_rss_mb": rss,
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


def view_command(
    view: str, path: Path, name: str, work: Path, focus: str | None = None, na: bool = False
) -> tuple[list[str], dict]:
    """The command and environment of one view's process: python -I, warnings shown (-W default), one BLAS thread,
    PYTHONHASHSEED 0 (1 for seed1). na: measure its nucleic acid (A17) though it is not the default view."""
    cmd = [sys.executable, "-I", "-W", "default", str(SCRIPT), "view", "--name", name, "--work", str(work)]
    cmd += ["--na"] if na else []
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
        self.sums["names_moved"] = 0  # helix names route._clear_helix_names moved off lines (random-100 C6)
        # route's second searches (C2b-2 fix round 1): nodes spent of the page's budget, searches run, refused
        self.sums.update(second_search_nodes=0, second_searches=0, second_searches_refused=0)
        self.records = 0

    def emit(self, record):
        if str(record.msg).startswith("%d helix names moved"):
            self.sums["names_moved"] += int(record.args[0])
            return
        if str(record.msg).startswith("%d second-search nodes spent"):
            spent, _, runs, refused = record.args
            self.sums["second_search_nodes"] += int(spent)
            self.sums["second_searches"] += int(runs)
            self.sums["second_searches_refused"] += int(refused)
            return
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


def measure_view(view: str, drawn, spec, path: Path, name: str, work: Path, na: bool = False) -> dict:
    """Everything one drawn figure view is judged on: what the figure gives (criteria.measure_figure), then what the
    campaign measures apart: its layout signature, loops and SVG hashes (A12, A13), na_validate's nucleic-acid
    measures (A17: the default view's, or with `na` this view's, when the default draws no nucleic acid: 8QEM's
    focus view, final review B, Important 5), and its 3D trace for the contact sheet."""
    bench = _script_module("bench")
    lay = drawn.layout
    out = measure_figure(view, drawn, spec, path)
    sig = json.dumps(_jsonable(bench.signature(lay)), sort_keys=True, default=str)
    loops = json.dumps([[l.a_id, l.b_id, [list(map(float, p)) for p in l.points]] for l in drawn.loops])
    svg = work / f"{name}-{view}.svg"
    out["signature_sha256"] = _sha256(sig.encode())
    out["loops_sha256"] = _sha256(loops.encode())
    out["svg_sha256"] = _sha256(svg.read_bytes()) if svg.is_file() else None
    if (view == "default" or na) and lay.nucleic is not None:
        out["na"] = _jsonable(na_metrics(drawn, path))
    if view == "default":
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
    out = {"refused": None, "residues": size["residues"], "svg_bytes": len(got["svg"])}
    try:
        out["science"], out["science_judged"] = app_science(path, got["svg"])
    except Exception:  # not measured, which fails the view (fail closed, as measure_view's errors)
        import traceback

        out["measure_error"] = traceback.format_exc()[-2000:]
    return out


def app_science(path: Path, svg: str) -> tuple[dict, dict]:
    """The science checks of the figure the app drew: the job process hands back its SVG (`svg`), not the figure, so
    the same spec (app.spec_of, the app's default settings) is drawn again here, outside the hosted limits, and
    measured, once its SVG is found equal to the app's (else ValueError: not that figure). The job's own warnings
    are the app view's; this measuring's are not. With what each check had to judge (figcheck.science_judged)."""
    import warnings

    from foldmap import app, figcheck, render
    from foldmap.spec import draw_spec

    settings = {**app.LAYOUT_DEFAULTS, "theme": "publication"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        drawn = draw_spec(app.spec_of(path, settings))
        again = render.save_svg(drawn.figure)
        if again != svg:
            raise ValueError(
                f"the redraw differs from the figure the app drew (SVG sha256 {_sha256(again.encode())[:12]} vs "
                f"{_sha256(svg.encode())[:12]})"
            )
        return _jsonable(figcheck.science_checks(drawn)), figcheck.science_judged(drawn)


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
            for key in ("science", "science_judged", "measure_error"):  # judged with the other views' (science, A1)
                if key in result["app"]:
                    result[key] = result["app"].pop(key)
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
                drawn, spec = caught["drawn"], caught["spec"]
                result.update(measure_view(args.view, drawn, spec, path, args.name, work, na=args.na))
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
    """criteria.dssp_agreement against mkdssp's codes (rows; None: run mkdssp, mkdssp_rows)."""
    return criteria_module.dssp_agreement(path, mkdssp_rows(path) if rows is None else rows)


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


def run_view(view: str, path: Path, name: str, work: Path, focus: str | None = None, na: bool = False) -> dict:
    """One view in its own process: what it exited with, its time and memory, what it printed, what it measured."""
    cmd, env = view_command(view, path, name, work, focus, na)
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


def focus_target(d: dict) -> str:
    """The focus view's --focus: 'none' (the whole file) when the default view drew part of it (a scope the chains
    or the nucleic-acid crop left out), else the protein chain with the most elements drawn alone (the first of
    them in drawing order; C25: chain "0" of 7JU4, a two-helix coiled coil, was taken by sort order), or the first
    protein chain when no element is recorded. A chain "not drawn" because it has nothing to draw (C15) is no
    scope: 7ARD and 7PIU keep their single-chain focus view."""
    scope = d.get("scope") or {}
    scoped = scope.get("chains_reason") or scope.get("crop")
    if scoped or not d.get("first_protein_chain"):
        return "none"
    count = Counter(e["chain"] for e in d.get("elements") or [])
    order = list(dict.fromkeys(e["chain"] for e in d.get("elements") or []))
    return max(order, key=lambda c: (count[c], -order.index(c))) if order else d["first_protein_chain"]


def run_entry(source: str, out: Path, views=VIEWS, info: dict | None = None) -> dict:
    """Draw and measure one entry (a sampled id, or a local file) in every view that applies; judged by criteria.
    info: the entry's line of sample.json (stratum, bin, pair ...); refusal: a CA-only or P-only id (A1c), drawn in
    the default view only."""
    info = dict(info or {})
    name = entry_name(source)
    work = out / "work" / name
    clear_work(work)
    keep = ("stratum", "bin", "order", "quota", "pair", "pair_type", "role", "partner", "refusal")
    e: dict = {"id": info.get("id", name), "name": name, **{k: info[k] for k in keep if k in info}, "views": {}}
    e["foldmap_commit"] = _commit()  # this entry's own: a report may join runs (final review B, Minor 6)
    try:
        path = resolve(source)
        e["file"], e["sha256"] = str(path), _sha256(path.read_bytes())
        e["facts"] = features(path, contacts=bool(e.get("pair")))  # a pair's type is re-derived (pair_retyped)
    except Exception as err:  # a download that fails fails the entry, not the campaign
        e["error"] = f"{type(err).__name__}: {err}"
        e["criteria"] = criteria(e)
        return _verdict(e)
    e["facts"]["resolution"] = _resolution(path)
    wanted = ("default",) if e.get("refusal") else tuple(views)

    def run(view: str, focus: str | None = None, na: bool = False) -> None:
        if view in wanted:
            e["views"][view] = run_view(view, path, name, work, focus, na)

    run("default")
    d = e["views"].get("default", {}).get("result") or {}
    f = e["facts"]
    if not e.get("refusal"):
        if f["protein"]:
            run("deposited")
            run("beta")
        if f["chains"] >= 2:  # the whole file when the default drew part, else the first protein chain alone
            run("focus", focus_target(d), na=bool(f["na"]) and not d.get("na_chains_drawn"))
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
            e["domains_cath"] = domain_reference(path, e["id"])  # S5, where the CATH reference holds the entry
    e["criteria"] = criteria(e)
    e["science"] = science(e)
    return _verdict(e)


KEEP_IN_WORK = ("entry.out", "entry.err")  # the entry process's own output, written by `run` before it starts


def clear_work(work: Path) -> None:
    """An entry's work folder, emptied of an earlier run's files (a view that writes nothing must not be judged on
    them: final review B, Minor 5), the entry process's own output kept."""
    work.mkdir(parents=True, exist_ok=True)
    for f in sorted(work.iterdir()):
        if f.name in KEEP_IN_WORK:
            continue
        if f.is_dir():
            shutil.rmtree(f)
        else:
            f.unlink()


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
        rss = "?" if r.get("max_rss_mb") is None else f"{r['max_rss_mb']:.0f}"
        items.append(f"{v}={r['seconds']}s/{rss}MB")
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
    items += [f"accessories_named.{v}={n}/{of}" for v, (n, of) in accessories_named(e).items() if of]
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


MARKS = {"pass": "✓", "fail": "✗", NOT_JUDGED: "○", NOT_APPLICABLE: "-"}


def criteria_counts(entries: list[dict]) -> list[str]:
    """report.md's per-criterion outcomes: pass, fail, not judged, does not apply (an entry passes when no criterion
    fails, so a criterion's pass count is the entries it judged and passed, never one with nothing to judge: final
    review B, Important 3), then every entry a criterion did not judge, with why (the plan lists A10's unmeasured)."""
    lines = ["", "## Outcomes per criterion", ""]
    lines += [
        f"✓ pass, ✗ fail, ○ {NOT_JUDGED} (the criterion applies but the entry gave it nothing to judge: no tie, "
        f"no tether to another chain, no domain drawn, mkdssp unmeasured ...), - {NOT_APPLICABLE}.",
        "",
        "| criterion | " + " | ".join(OUTCOMES) + " |",
        "|---|---|---|---|---|",
    ]
    for key, what in CRITERIA.items():
        got = Counter(outcome((e.get("criteria") or {}).get(key)) for e in entries)
        lines.append(f"| {key} {what} | " + " | ".join(str(got[o]) for o in OUTCOMES) + " |")
    lines += ["", "### Not judged, per criterion", ""]
    for key in CRITERIA:
        left = [
            f"{e['id']} ({((e.get('criteria') or {})[key].get('detail') or 'no detail')[:160]})".replace("|", "/")
            for e in entries
            if outcome((e.get("criteria") or {}).get(key)) == NOT_JUDGED
        ]
        if left:
            lines.append(f"- **{key}** ({len(left)}): " + "; ".join(left))
    if not any(line.startswith("- **") for line in lines):
        lines.append("None.")
    return lines


def pairs_section(entries: list[dict]) -> list[str]:
    """report.md's two-state pairs (S11): each pair's type as sampled and as its files read now (pair_retyped, with
    the additives foldmap hides left out of the ligands, final review B, Important 4; CA RMSD over residues aligned
    by sequence, Minor 9), and A19 on both entries."""
    seeds = [e for e in entries if e.get("role") == "seed"]
    if not seeds:
        return []
    lines = ["", "## Two-state pairs (S11)", ""]
    lines += [
        "| pair | seed | partner | sampled as | types now | ligand touching | assembly 1 chains | CA RMSD Å | A19 |"
    ]
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for e in seeds:
        r = e.get("pair_retyped") or {}
        a19 = (e.get("criteria") or {}).get("A19")
        lines.append(
            f"| {e.get('pair')} | {e['id']} | {e.get('partner')} | {e.get('pair_type')} | "
            f"{r.get('note') if r.get('types') is None else ', '.join(r['types']) or 'none'} | "
            f"{r.get('ligand_touching')} | {r.get('assembly1_chains')} | "
            f"{r.get('ca_rmsd')} | {MARKS[outcome(a19)]} {(a19 or {}).get('detail', '')[:300].replace('|', '/')} |"
        )
    return lines


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
    lines += criteria_counts(entries)
    lines += [
        "",
        "## Entries",
        "",
        "process s: the default view's whole process (drawing, measuring and, for nucleic acid, RNApolis), which A2 "
        "gates; drawing s: the drawing inside it.",
        "",
        "| id | stratum | bin | chains | process s | drawing s | " + " | ".join(keys) + " | failed |",
    ]
    lines.append("|" + "---|" * (len(keys) + 7))
    for e in entries:
        crit = e.get("criteria") or {}
        marks = [MARKS[outcome(crit.get(k))] for k in keys]
        default = e.get("views", {}).get("default") or {}
        took, inside = default.get("seconds", ""), (default.get("result") or {}).get("seconds_inside", "")
        detail = "; ".join(e.get("failed", [])).replace("|", "/")[:600]
        chains = (e.get("facts") or {}).get("chains", "")
        lines.append(
            f"| {e['id']} | {e.get('stratum', '')} | {e.get('bin', '')} | {chains} | {took} | {inside} | "
            + " | ".join(marks)
            + f" | {detail} |"
        )
    lines += pairs_section(entries)
    lines += [
        "",
        "## Reported, not gated",
        "",
        "| id | " + " | ".join(REPORTED) + " |",
        "|" + "---|" * (len(REPORTED) + 1),
    ]
    for e in entries:
        lines.append(f"| {e['id']} | " + " | ".join(str(v) for v in reported(e).values()) + " |")
    lines += science_section(entries)
    lines += domains_section(entries)
    lines += accessories_section(entries)
    lines += sheets_differ_section(entries)
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


DOMAIN_REFERENCE = ROOT / "tests" / "data" / "domains" / "cath_reference.json"  # S5's choppings (CATH, CC BY 4.0)


def domain_reference(path: Path, pdb: str) -> list[dict] | None:
    """S5's measures: each reference chain of entry `pdb` (DOMAIN_REFERENCE's multi-domain CATH chains) scored with
    its automatic domains (domain_score.chain_score: overlap, boundary recall, precision and F1); None when the
    reference does not hold the entry."""
    from foldmap import domain_score

    chains = [c for c in domain_score.load_reference(DOMAIN_REFERENCE) if c["pdb"] == pdb.upper()]
    return [domain_score.chain_score(path, c) for c in chains] or None


SCIENCE_VIEWS = (*FIGURE_VIEWS, "app")  # the views whose figure the science checks judge (batch C2b-S)


def _view_science(e: dict) -> dict[str, dict]:
    """Each measured view's science findings ({check: [findings] or None}), in SCIENCE_VIEWS order."""
    views = e.get("views") or {}
    out = {}
    for v in SCIENCE_VIEWS:
        found = ((views.get(v) or {}).get("result") or {}).get("science")
        if found is not None:
            out[v] = found
    return out


def _view_judged(e: dict) -> dict[str, dict]:
    """Each measured view's science_judged (what each check had to judge there: figcheck.science_judged)."""
    views = e.get("views") or {}
    return {v: (views[v].get("result") or {}).get("science_judged") for v in _view_science(e)}


def science(e: dict) -> dict:
    """The science checks (foldmap.science.SCIENCE) of every drawn view (SCIENCE_VIEWS), apart from A1-A20 (never in
    "failed" or "pass"): {"S1": {"pass": True, False or None, "detail", and for None its "status"}, ...}. A check
    fails when any view has a finding, each named by its view ("focus: sheet 27: ..."); it passes when some view had
    something to judge (science_judged: a sheet, a barrel, a mark, an accessory) and none has a finding; else it
    does not apply to an entry with no protein, and is not judged for one with protein (final review B, Important 3:
    nothing to judge is no pass). S5 is the entry's, once (chain-level, against the CATH reference)."""
    from foldmap.science import NOTHING, SCIENCE, domain_findings

    per_view, judged = _view_science(e), _view_judged(e)
    protein = bool((e.get("facts") or {}).get("protein"))
    out = {}
    for key in SCIENCE:
        if key == "S5":
            got = domain_findings(e["domains_cath"]) if per_view and e.get("domains_cath") else None
            if got is not None:
                out[key] = _ok(not got, "; ".join(got))
            else:
                out[key] = _ok(None, "no CATH reference chopping", NOT_JUDGED if protein else NOT_APPLICABLE)
            continue
        measured = [(v, f[key]) for v, f in per_view.items() if f.get(key) is not None]
        found = [f"{v}: {t}" for v, got in measured for t in got]
        some = any((judged.get(v) or {}).get(key, 1) for v, _ in measured)  # a view from before: judged
        if found or (measured and some):
            out[key] = _ok(not found, "; ".join(found))
        elif not measured:
            out[key] = _ok(None, "not measured", NOT_JUDGED if protein else NOT_APPLICABLE)
        else:
            out[key] = _ok(None, f"{NOTHING[key]} in any view", NOT_JUDGED if protein else NOT_APPLICABLE)
    return out


def science_views(e: dict) -> dict[str, dict[str, int]]:
    """Per measured view, its count of findings per check (checks with none left out); S5 is the entry's."""
    return {v: {k: len(t) for k, t in f.items() if k != "S5" and t} for v, f in _view_science(e).items()}


def science_section(entries: list[dict]) -> list[str]:
    """report.md's science checks: per check how many entries pass, fail or are not measured (an entry fails a check
    when any of its views has a finding), then per view how many entries it was measured in and how many of them it
    fails per check, then each entry with its marks and findings (each named by its view)."""
    from foldmap.science import SCIENCE

    rows = [(e, science(e)) for e in entries]
    lines = ["", "## Science checks (S1-S6), apart from A1-A20", ""]
    lines += [f"Judged in every drawn view ({', '.join(SCIENCE_VIEWS)}); S5 once per entry.", ""]
    lines += [
        "A check passes only where it had something to judge; an entry with protein and nothing to judge (no sheet "
        "drawn, no barrel, no residue mark, no accessory strand, no CATH reference) is not judged, one with no "
        "protein does not apply.",
        "",
    ]
    lines += ["| check | " + " | ".join(OUTCOMES) + " |", "|---|---|---|---|---|"]
    for key, what in SCIENCE.items():
        got = Counter(outcome(s[key]) for _, s in rows)
        lines.append(f"| {key} {what} | " + " | ".join(str(got[o]) for o in OUTCOMES) + " |")
    keys = [k for k in SCIENCE if k != "S5"]
    lines += [
        "",
        "### Science by view",
        "",
        "Entries measured in each view, and how many of them fail each check there.",
        "",
    ]
    lines += ["| view | measured | " + " | ".join(keys) + " |", "|" + "---|" * (len(keys) + 2)]
    counts = [science_views(e) for e, _ in rows]
    for v in SCIENCE_VIEWS:
        here = [c[v] for c in counts if v in c]
        if here:
            lines.append(
                f"| {v} | {len(here)} | " + " | ".join(str(sum(1 for c in here if c.get(k))) for k in keys) + " |"
            )
    lines += ["", "| id | " + " | ".join(SCIENCE) + " | findings |", "|" + "---|" * (len(SCIENCE) + 2)]
    for (e, s), per in zip(rows, (_view_science(e) for e, _ in rows)):
        marks = [MARKS[outcome(s[k])] for k in SCIENCE]
        found = [f"{v}: {k} {t}" for v, f in per.items() for k in keys for t in (f.get(k) or [])]
        if s["S5"]["pass"] is False:
            found.append(f"S5 {s['S5']['detail']}")
        lines.append(f"| {e['id']} | " + " | ".join(marks) + f" | {'; '.join(found).replace('|', '/')[:600]} |")
    return lines


def domains_section(entries: list[dict]) -> list[str]:
    """report.md's S5 measures: each entry the CATH reference holds, per chain its reference and automatic domain
    counts, overlap, boundary F1, recall and precision (domain_reference), then their means. Reported, not gated;
    in-sample (the reference chains are the calibration set)."""
    rows = [(e["id"], s) for e in entries for s in e.get("domains_cath") or []]
    lines = ["", "## S5: automatic domains against CATH (reported, not gated)", ""]
    lines += [
        "In-sample: these chains are among the 33 the automatic domains were calibrated on "
        "(tests/data/domains/cath_reference.json); the never-seen check is the fresh set (cath_fresh.json).",
        "",
    ]
    if not rows:
        return [*lines, "No entry is in the CATH reference."]
    lines += ["| id | chain | CATH domains | automatic | overlap | boundary F1 | recall | precision |"]
    lines.append("|---|---|---|---|---|---|---|---|")
    for name, s in rows:
        lines.append(
            f"| {name} | {s['chain']} | {s['n_ref']} | {s['n_auto']} | {s['overlap']:.2f} | {s['f1']:.2f} | "
            f"{s['recall']:.2f} | {s['precision']:.2f} |"
        )
    mean = {k: sum(s[k] for _, s in rows) / len(rows) for k in ("overlap", "f1", "recall", "precision")}
    lines.append(
        f"| mean of {len(rows)} | | | | {mean['overlap']:.2f} | {mean['f1']:.2f} | {mean['recall']:.2f} | "
        f"{mean['precision']:.2f} |"
    )
    return lines


def accessories_named(e: dict) -> dict[str, tuple[int, int]]:
    """Per figure view, (accessory strands named rather than drawn level beside their partner, accessory strands):
    figcheck.figure_checks' "accessories_named" of "accessories" (plan 12 ruling: S6 passes a named one, so this
    count measures placement). Reported, never in "failed" or "pass"; views without a measure are left out."""
    out = {}
    for v, r in e.get("views", {}).items():
        ch = (r.get("result") or {}).get("checks") or {}
        if ch.get("accessories") is not None:
            out[v] = (ch.get("accessories_named") or 0, ch["accessories"])
    return out


def sheets_differ_section(entries: list[dict]) -> list[str]:
    """report.md's copies whose sheets DSSP groups differently (figcheck._sheets_differ, ruling C2b-3/R2): for each
    entry and A15 view, every tie whose differing sheets' strands A15 left out, in the ruling's own words, with
    what was still judged."""
    lines = ["", "## Copies' sheets that differ, left out of A15 (ruling C2b-3/R2: reported, not judged)", ""]
    lines += ["| id | view | tie | sheets left out | judged |", "|---|---|---|---|---|"]
    for e in entries:
        for v in ("default", "focus"):
            details = ((((e.get("views") or {}).get(v) or {}).get("result") or {}).get("checks") or {}).get(
                "details"
            ) or {}
            for tie, got in sorted((details.get("copies_alike") or {}).items()):
                if got.get("sheets_differ"):
                    judged = got.get("judged")
                    judged = f"yes, {got.get('matched')} matched" if judged == "yes" else judged
                    lines.append(f"| {e['id']} | {v} | {tie} | {got['sheets_differ']} | {judged} |")
    return lines


def accessories_section(entries: list[dict]) -> list[str]:
    """report.md's accessory placement: each entry with accessory strands, per view how many are named rather than
    drawn level beside their partner (accessories_named)."""
    lines = ["", "## Accessories named, not drawn level beside their partner (reported, not gated)", ""]
    lines += ["| id | named of accessory strands, per view |", "|---|---|"]
    for e in entries:
        got = accessories_named(e)
        if any(of for _, of in got.values()):
            lines.append(f"| {e['id']} | " + "; ".join(f"{v} {n} of {of}" for v, (n, of) in got.items()) + " |")
    return lines


REPORTED = (
    "honest_crossings",
    "loop_crossings",
    "na_run_crossings",
    "fallback_curves",
    "fallback_curves_crossing_an_element",
    "kink_links_rejected",
    "tighter_than_inflate",
    "rejected_as_crossing",
    "names_moved",
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
        "names_moved": route.get("names_moved"),  # by the post-route mover (C6): lines no route holds
        "area": d.get("area"),
        "max_rss_mb": max(
            (r["max_rss_mb"] for r in e.get("views", {}).values() if r.get("max_rss_mb") is not None), default=None
        ),
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
        if e.get("role") == "seed" and partner is not None:
            e["pair_retyped"] = pair_retyped(e, partner)
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


def pair_retyped(seed: dict, partner: dict) -> dict:
    """A pair's S11 types as its files read now (pair_types: features with contacts, the additives left out, and CA
    RMSD over sequence-aligned residues): sample.json keeps the type it was drawn as, the report says both."""
    fa, fb = seed.get("facts") or {}, partner.get("facts") or {}
    if "ligand_touching" not in fa or "ligand_touching" not in fb:
        return {"types": None, "note": "no contacts measured"}
    rmsd = ca_rmsd(seed["file"], partner["file"]) if seed.get("file") and partner.get("file") else None
    return {
        "types": sorted(pair_types(fa, fb, rmsd)),
        "ligand_touching": [fa["ligand_touching"], fb["ligand_touching"]],
        "ligands": [fa.get("ligands"), fb.get("ligands")],
        "assembly1_chains": [fa.get("assembly1_chains"), fb.get("assembly1_chains")],
        "ca_rmsd": rmsd,
    }


def read_report(path: Path) -> dict:
    """A report.json, or the report.json.gz a large report is committed as."""
    path = Path(path)
    if not path.is_file() and Path(f"{path}.gz").is_file():
        path = Path(f"{path}.gz")
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == ".gz" else data)


def cmd_run(args) -> int:
    from concurrent.futures import ThreadPoolExecutor

    out, sample_path = Path(args.out), Path(args.sample)
    listed = _sample_entries(sample_path) if sample_path.is_file() else []
    if args.ids:
        ids = args.ids.split(",")
    elif args.failed_from:
        old = read_report(Path(args.failed_from))
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
    s.add_argument("--exclude", metavar="SAMPLE", help="another sample.json whose ids are never drawn (held-out)")
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
    v.add_argument("--na", action="store_true", help="measure the nucleic acid (A17) in this view too")
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
