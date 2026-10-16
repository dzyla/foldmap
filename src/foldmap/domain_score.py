# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Automatic domains against a reference chopping (plan 12 Task 3, the external review's R2; science check S5).

A chain's domains are compared residue by residue, over the chain's residues in the file (its Backbone labels):
- the reference: each residue in a reference segment takes that domain; residues no segment covers are not judged;
- the automatic domains: each residue of a domain's element takes that domain, a loop residue the domain of the
  nearest element residue (the earlier one on a tie), so the whole chain is assigned.

Two measures, per chain:
- overlap: the reference and automatic domains matched one to one (Hungarian, most shared residues), the fraction of
  judged residues in their matched domain. Splitting a reference domain costs what the extra pieces hold;
- boundaries: a boundary is where the domain changes between consecutive residues (judged residues only, for the
  reference); a reference boundary is found when an automatic one lies within TOLERANCE residues of it, an
  automatic one is right when a reference one lies that near. Recall, precision and their F1 (a chain with no
  boundary on either side scores 1 on both).

The calibration set (tests/data/domains/cath_reference.json, CATH via the PDBe API, CC BY 4.0) holds 33
multi-domain chains (plan 12 ruling: mean overlap at least TARGETS["overlap"], mean boundary F1 at least
TARGETS["f1"], mean recall at least TARGETS["recall"])."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

TOLERANCE = 20  # residues: a boundary this near the other side's counts as found
TARGETS = {"overlap": 0.80, "f1": 0.70, "recall": 0.85}  # plan 12 ruling R2, means over the CATH set


def reference_labels(seqs, domains) -> np.ndarray:
    """Per residue (author numbers `seqs`, in chain order) the index of the reference domain whose segments
    (`domains`: per domain [[first, last], ...], author numbers) hold it, -1 where none does."""
    seqs = np.asarray(seqs)
    out = np.full(len(seqs), -1)
    for d, segments in enumerate(domains):
        for a, b in segments:
            out[(seqs >= a) & (seqs <= b)] = d
    return out


def element_labels(n: int, spans) -> np.ndarray:
    """Per residue 0..n-1 the domain of the element holding it (`spans`: (first index, last index, domain)), a
    loop residue the domain of the nearest element residue (the earlier on a tie); all 0 with no element."""
    out = np.full(n, -1)
    for a, b, d in spans:
        out[max(a, 0) : b + 1] = d
    known = np.flatnonzero(out >= 0)
    if not len(known):
        return np.zeros(n, dtype=int)
    for k in np.flatnonzero(out < 0):
        out[k] = out[known[np.argmin(np.abs(known - k))]]
    return out


def overlap(ref: np.ndarray, auto: np.ndarray) -> float:
    """The fraction of judged residues (ref >= 0) in the automatic domain matched one to one to theirs."""
    from scipy.optimize import linear_sum_assignment

    judged = ref >= 0
    if not judged.any():
        return 1.0
    r, a = ref[judged], np.unique(auto, return_inverse=True)[1].reshape(-1)[judged]
    k = max(r.max(), a.max()) + 1
    shared = np.zeros((k, k))
    np.add.at(shared, (r, a), 1)
    rows, cols = linear_sum_assignment(-shared)
    return float(shared[rows, cols].sum() / judged.sum())


def boundaries(labels: np.ndarray, judged: np.ndarray | None = None) -> list[float]:
    """Where the domain changes along the chain: the midpoint of each pair of consecutive judged residues of
    different domains (residue indices)."""
    idx = np.flatnonzero(judged if judged is not None else np.ones(len(labels), bool))
    return [(a + b) / 2 for a, b in zip(idx, idx[1:]) if labels[a] != labels[b]]


def boundary_scores(ref: np.ndarray, auto: np.ndarray, tolerance: float = TOLERANCE) -> tuple[float, float, float]:
    """(recall, precision, F1) of the automatic boundaries against the reference ones, within `tolerance`."""
    rb, ab = boundaries(ref, ref >= 0), boundaries(auto)
    found = sum(any(abs(b - x) <= tolerance for x in ab) for b in rb)
    right = sum(any(abs(b - x) <= tolerance for b in rb) for x in ab)
    recall = found / len(rb) if rb else 1.0
    precision = right / len(ab) if ab else (1.0 if not rb else 0.0)
    f1 = 2 * recall * precision / (recall + precision) if recall + precision else 0.0
    return recall, precision, f1


def score(ref: np.ndarray, auto: np.ndarray) -> dict:
    """overlap, recall, precision, f1 and the domain counts of one chain."""
    recall, precision, f1 = boundary_scores(ref, auto)
    return {
        "overlap": overlap(ref, auto),
        "recall": recall,
        "precision": precision,
        "f1": f1,
        "n_ref": int(len(set(ref[ref >= 0].tolist()))),
        "n_auto": int(len(set(auto.tolist()))),
    }


def summary(scores: list[dict]) -> dict:
    """Means over chains of overlap, recall, precision and F1 (the ruled measures), with the chain count."""
    keys = ("overlap", "recall", "precision", "f1")
    return {"n": len(scores), **{k: float(np.mean([s[k] for s in scores])) if scores else 0.0 for k in keys}}


def layout_labels(domains, sses, labels, chain: str) -> np.ndarray:
    """A chain's residues (its `labels`, Backbone order) labelled by the domains (Layout.domains: (name, element
    ids)) through their elements (element_labels): domains numbered in their order."""
    idx = [k for k, l in enumerate(labels) if l.chain == chain]
    pos = {k: n for n, k in enumerate(idx)}
    owner = {i: d for d, (_, ids) in enumerate(domains) for i in ids}
    spans = []
    for s in sses:
        if s.chain == chain and s.id in owner:
            inside = [pos[k] for k in range(s.start, s.end + 1) if k in pos]
            if inside:
                spans.append((inside[0], inside[-1], owner[s.id]))
    return element_labels(len(idx), spans)


@dataclass
class ChainInputs:
    """What one chain's automatic domains are found from (features.auto_domains): its residues (labels and CA
    coordinates), its elements, their contacts (with disulfide ties) and the βα rows a template draws whole."""

    labels: list
    ca: np.ndarray
    sses: list
    contacts: dict
    rows: list


def chain_inputs(path, chain: str) -> ChainInputs:
    """A chain's ChainInputs as `--domains auto --chains <chain> --assembly asu` reads them (cli.make_layout)."""
    from .cli import _bridge_contacts
    from .dssp import assign_dssp
    from .features import sse_contacts
    from .io import load_backbone, load_links, read_model
    from .layout import template_rows
    from .model import homolog_chains
    from .sheets import build_sheets
    from .ss import build_sses

    model = read_model(path, "asu")
    bb = load_backbone(path, "asu", model=model).only([chain])
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, short_helices=True, bridges=dssp.bridges, matched={})
    sheets = build_sheets(sses, dssp.bridges)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, load_links(path, bb, "asu", model=model)).items():
        contacts[key] = contacts.get(key, 0) + n
    rows = template_rows(sheets, sses, homolog_chains(bb))
    return ChainInputs(list(bb.labels), np.asarray(bb.ca, dtype=float), sses, contacts, rows)


def inputs_domains(inputs: ChainInputs) -> list:
    """The automatic domains of one chain's inputs (features.auto_domains, as cli._domains calls it)."""
    from types import SimpleNamespace

    from .features import auto_domains

    bb = SimpleNamespace(labels=inputs.labels, ca=inputs.ca)
    return auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb)


def auto_chain(path, chain: str) -> tuple[list, np.ndarray, list]:
    """One chain's automatic domains, without the figure: (its residue labels, per residue its domain
    (layout_labels), the domains)."""
    inputs = chain_inputs(path, chain)
    domains = inputs_domains(inputs)
    return inputs.labels, layout_labels(domains, inputs.sses, inputs.labels, chain), domains


def inputs_score(inputs: ChainInputs, entry: dict) -> dict:
    """One reference entry scored from its chain's inputs (a fixture's, or chain_inputs')."""
    chain = entry["chain"]
    auto = layout_labels(inputs_domains(inputs), inputs.sses, inputs.labels, chain)
    ref = reference_labels([l.seq for l in inputs.labels if l.chain == chain], entry["domains"])
    return {"pdb": entry["pdb"], "chain": chain, **score(ref, auto)}


# A fixture of ChainInputs (tests/data/domains/*.npz, scripts/make_domain_fixture.py): per chain "<PDB>_<chain>",
# CA coordinates in thousandths of an Å (as the mmCIF file writes them), residue numbers and insertion codes,
# element spans (indices into the residues) and kinds, element contacts (element index pairs and counts) and the
# βα rows (element indices, -1 between rows).
def fixture_arrays(inputs: ChainInputs) -> dict[str, np.ndarray]:
    """One chain's inputs as the fixture's arrays."""
    at = {s.id: k for k, s in enumerate(inputs.sses)}
    contacts = sorted(
        (min(at[a], at[b]), max(at[a], at[b]), n)
        for (a, b), n in ((sorted(k), n) for k, n in inputs.contacts.items())
        if a in at and b in at
    )
    rows = [k for row in inputs.rows for k in [*(at[i] for i in row if i in at), -1]]
    return {
        "ca": np.round(np.asarray(inputs.ca) * 1000).astype(np.int32),
        "seq": np.array([l.seq for l in inputs.labels], dtype=np.int32),
        "icode": np.array([l.icode for l in inputs.labels], dtype="<U2"),
        "chain": np.array([inputs.labels[0].chain if inputs.labels else ""], dtype="<U4"),
        "sse": np.array([[s.start, s.end] for s in inputs.sses], dtype=np.int32).reshape(-1, 2),
        "kind": np.array([s.kind for s in inputs.sses], dtype="<U1"),
        "contacts": np.array(contacts, dtype=np.int32).reshape(-1, 3),
        "rows": np.array(rows, dtype=np.int32),
    }


def save_fixture(path, chains: dict[str, ChainInputs]) -> None:
    """Write a fixture of ChainInputs, keyed "<PDB>_<chain>" (np.savez_compressed)."""
    arrays = {f"{key}__{name}": a for key, inputs in chains.items() for name, a in fixture_arrays(inputs).items()}
    np.savez_compressed(path, **arrays)


def load_fixture(path) -> dict[str, ChainInputs]:
    """A fixture's ChainInputs by "<PDB>_<chain>": the elements rebuilt from their spans (no centroid or axis:
    auto_domains does not read them)."""
    from .model import SSE, ResLabel

    data = np.load(path)
    keys = sorted({name.split("__")[0] for name in data.files})
    out = {}
    for key in keys:
        a = {name.split("__")[1]: data[name] for name in data.files if name.split("__")[0] == key}
        chain = str(a["chain"][0])
        labels = [ResLabel(chain, int(n), str(i), "") for n, i in zip(a["seq"], a["icode"])]
        zero = np.zeros(3)
        sses = [
            SSE(str(kind), chain, int(lo), int(hi), labels[lo], labels[hi], zero, zero)
            for (lo, hi), kind in zip(a["sse"], a["kind"])
        ]
        contacts = {frozenset((sses[i].id, sses[j].id)): int(n) for i, j, n in a["contacts"]}
        rows, row = [], []
        for k in a["rows"]:
            if k < 0:
                rows.append(row)
                row = []
            else:
                row.append(sses[k].id)
        out[key] = ChainInputs(labels, a["ca"].astype(float) / 1000, sses, contacts, rows)
    return out


def load_reference(path) -> list[dict]:
    """The reference choppings: [{pdb, chain, tag, source, domains: [[[first, last], ...], ...]}]."""
    return json.loads(Path(path).read_text(encoding="utf-8"))["chains"]


def chain_score(path, entry: dict) -> dict:
    """One reference entry scored from its structure file: its automatic domains against entry["domains"]."""
    return inputs_score(chain_inputs(path, entry["chain"]), entry)


HELD_OUT_SEED = "20261016:domains"  # random.Random seed of the held-out third


def held_out(entries: list[dict]) -> set[tuple[str, str]]:
    """The held-out third of the calibration chains, (pdb, chain), drawn by homology group (each entry's "group":
    chains sharing a CATH superfamily) so that no held-out chain has a homologue among the rest: the groups in the
    order HELD_OUT_SEED shuffles them, taken whole until a third of the chains is held out.

    Fix round 1 (plan 12 Task 3 review): the first split drew chains, not groups, and leaked near-duplicates of
    training chains into the held-out third (6FY9/6FYF, 2YBE/3PGK, the hexokinases). The parameters were chosen on
    that split's training two thirds, and that split was looked at twice, so this one is not untouched either: it
    is reported, and the never-seen check is the fresh set (cath_fresh.json)."""
    import random

    groups: dict[str, list[tuple[str, str]]] = {}
    for e in sorted(entries, key=lambda e: (e["pdb"], e["chain"])):
        groups.setdefault(e.get("group") or e["pdb"], []).append((e["pdb"], e["chain"]))
    order = sorted(groups)
    random.Random(HELD_OUT_SEED).shuffle(order)
    held: set[tuple[str, str]] = set()
    for g in order:
        if len(held) >= len(entries) // 3:
            break
        held.update(groups[g])
    return held
