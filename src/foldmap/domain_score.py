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


def auto_chain(path, chain: str) -> tuple[list, np.ndarray, list]:
    """One chain's automatic domains as `--domains auto --chains <chain> --assembly asu` finds them, without the
    figure: (its residue labels, per residue its domain (layout_labels), the domains). The same inputs as
    cli.make_layout: the elements, their contacts and disulfide ties, and the βα rows a template draws whole."""
    from .cli import _bridge_contacts, _domains
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
    copies = homolog_chains(bb)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, load_links(path, bb, "asu", model=model)).items():
        contacts[key] = contacts.get(key, 0) + n
    domains = _domains("auto", sses, contacts, template_rows(sheets, sses, copies), copies)
    return list(bb.labels), layout_labels(domains, sses, bb.labels, chain), domains


def load_reference(path) -> list[dict]:
    """The reference choppings: [{pdb, chain, tag, source, domains: [[[first, last], ...], ...]}]."""
    return json.loads(Path(path).read_text(encoding="utf-8"))["chains"]


def chain_score(path, entry: dict) -> dict:
    """One reference entry scored: its auto_chain domains against entry["domains"]."""
    labels, auto, domains = auto_chain(path, entry["chain"])
    ref = reference_labels([l.seq for l in labels if l.chain == entry["chain"]], entry["domains"])
    return {"pdb": entry["pdb"], "chain": entry["chain"], **score(ref, auto)}


HELD_OUT_SEED = "20261016:domains"  # random.Random seed of the held-out third (plan 12 Task 3: tuned on the rest)


def held_out(entries: list[dict]) -> set[tuple[str, str]]:
    """The held-out third of the calibration chains, (pdb, chain): drawn once with HELD_OUT_SEED, never tuned on."""
    import random

    keys = sorted((e["pdb"], e["chain"]) for e in entries)
    return set(random.Random(HELD_OUT_SEED).sample(keys, len(keys) // 3))
