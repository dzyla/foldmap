# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The record of the automatic domains' uncertainty (plan 28 Task 5, the review's 3.4, rulings Q42, Q43): why they
fail against CATH, what each split's confidence is worth, and what keeping a weakly supported region undivided would
do.

    python scripts/domain_uncertainty.py [--out docs/validation/domain-uncertainty] [--jobs 12]
                                         [--against REV [--against-source FILE]]

Writes report.json and report.md (the summaries) and report-chains.json (each chain's cuts, trace and scores at every
t), all stamped with the code that measured (foldmap.source_stamp). Reads the structures of the third random sample
from ~/.cache/foldmap (fetched when missing) and CATH's choppings through foldmap.cath (cached in
~/.cache/foldmap/cath); the domain sets from their committed fixtures. Not run in CI.

--against REV compares the domains (features.auto_domains) of every chain with those of features.py at commit REV,
loaded beside this one (`git show REV:src/foldmap/features.py`, or --against-source FILE where the checkout has no
git, as on a runner: FILE is that text, and its SHA-256 is recorded): {against, source_sha256, chains, differ}.

Sets, with their tags read from scripts/paper/sets.json: dom_cal (the 22 calibration chains: cath_reference.json
minus the third held out chain by chain, seed domain_score.HELD_OUT_SEED, as the parameters were chosen on them and
as scripts/paper/domains.py splits them), dom_held (that third), dom_fresh (cath_fresh.json), all from their
committed fixtures; r100_third (the third random sample: every chain CATH chops, scored as
random100.domain_reference scores it; looked at again here). Every chain is scored by domain_score (overlap, boundary
recall, precision and F1 against CATH, TOLERANCE residues) and meets the targets when its overlap, F1 and recall are
each at least domain_score.TARGETS (S5, science.domain_findings).

1. diagnosis (read-only): each chain of the third sample with two CATH domains or more that fails S5, classed by
   what went wrong (exclusive: fewer automatic domains than CATH's, "boundary missed"; as many, "boundary moved"; more,
   "extra split") and tagged by what the chain is like ("discontinuous domain": a CATH domain in two runs or more
   along the chain; "tightly packed": a CATH domain whose normalized cut against the rest of the chain, on the
   contacts residue_domains weighs, is at DOMAIN_NCUT or above, so the criterion cannot see it; "large": three CATH
   domains or more, or more residues than DOMAIN_MAX_POINTS), with the mechanism in features.py named from what the
   cut did (its record: features.residue_domains' `cuts`). A missed boundary whose residue cut was made is traced
   through _chain_domains' steps (trace: the element vote, _no_strays, fold, the βα row merge, the second fold) to
   the first that brings the domain count below CATH's, with its line in features.py.
2. splits: every inferred split the drawn domains keep (features.domain_splits) on every set, with its confidence,
   called right when each place its two sides meet lies within TOLERANCE residues of a CATH boundary; per set, the
   fraction right in each confidence bin with its Wilson 95% interval (scripts/paper/common.wilson), and the AUC of
   the confidence between right and wrong splits.
3. undivided: the reading measured is the literal one, a region whose cut has confidence t or less is not cut and
   never weighed again (its stronger sub-cuts go with it): the same as making a cut only below DOMAIN_NCUT * (1 - t),
   as the confidence is the margin below that limit (up to its rounding to 0.01); measured by setting
   features.DOMAIN_NCUT so, for t on GRID, on every set. Merging only across a weak drawn split is another reading,
   not measured. t is chosen on dom_cal (tuning) only: the largest t at which no dom_cal chain's overlap, F1 or recall
   is lower than at t = 0. It is then measured on dom_held, dom_fresh and the third sample; ruling Q42 commits the
   change only when no set is worse (chains meeting the targets, mean overlap, F1 and recall; the third sample's
   entries failing S5 too)."""

from __future__ import annotations

import argparse
import functools
import hashlib
import importlib.util
import json
import random
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__" and (ROOT / "src" / "foldmap").is_dir():
    sys.path.insert(0, str(ROOT / "src"))  # the foldmap of this checkout (child processes inherit sys.path)

from foldmap import domain_score as ds  # noqa: E402
from foldmap import features  # noqa: E402
from foldmap.source_stamp import source_stamp, stamp_line  # noqa: E402

SOURCE_STAMP = source_stamp(__file__)  # the code as this process starts
OUT = ROOT / "docs" / "validation" / "domain-uncertainty"
DOMAINS = ROOT / "tests" / "data" / "domains"
THIRD = ROOT / "docs" / "validation" / "random100-third"
SETS_FILE = ROOT / "scripts" / "paper" / "sets.json"
GRID = tuple(round(0.05 * k, 2) for k in range(11))  # t: 0 (as drawn) to 0.5
BASE_NCUT = features.DOMAIN_NCUT
CLASSES = ("boundary missed", "boundary moved", "extra split")
TAGS = ("discontinuous domain", "tightly packed", "large")
BINS = ((0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.01))
SET_NAMES = ("dom_cal", "dom_held", "dom_fresh", "r100_third")
LOOKED_AGAIN = {"r100_third": "looked at again here"}
# _chain_domains' steps (trace), each found in features.py by the text of its line
STEPS = (
    ("the element vote", "mine.append(min(votes"),
    ("_no_strays", "_no_strays(mine,"),
    ("fold (DOMAIN_MIN_SSE_RESIDUES)", "groups = fold(found.values())"),
    ("the βα row merge", "for keep in together:"),
    ("the second fold", "    groups = fold(groups)"),
)


def _wilson():
    """scripts/paper/common.wilson, loaded by path under a private name (as na_validate loads it)."""
    name = "_domain_uncertainty_paper_common"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "paper" / "common.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[name].wilson


def set_tags() -> dict[str, str]:
    """Each set's tag as scripts/paper/sets.json gives it, with what this record adds."""
    sets = json.loads(SETS_FILE.read_text(encoding="utf-8"))["sets"]
    return {k: sets[k]["tag"] + (f"; {LOOKED_AGAIN[k]}" if k in LOOKED_AGAIN else "") for k in SET_NAMES}


# ------------------------------------------------------------------------------------------------ one chain


def _meets(s: dict) -> bool:
    return all(s[k] >= ds.TARGETS[k] for k in ("overlap", "f1", "recall"))


def _score(ref: np.ndarray, auto: np.ndarray) -> dict:
    got = ds.score(ref, auto)
    return {**{k: round(v, 4) if isinstance(v, float) else v for k, v in got.items()}, "meets": _meets(got)}


def _domains(inputs, splits=None, module=features):
    bb = SimpleNamespace(labels=inputs.labels, ca=inputs.ca)
    kw = {} if splits is None else {"splits": splits}
    return module.auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb, **kw)


def _meet_points(a, b) -> list[float]:
    """Where a split's two sides (residue indices) meet along the chain: the midpoints between consecutive residues
    of the part on different sides."""
    side = {int(k): 0 for k in a} | {int(k): 1 for k in b}
    order = sorted(side)
    return [(x + y) / 2 for x, y in zip(order, order[1:]) if side[x] != side[y]]


def _cath_ncut(ca: np.ndarray, ref: np.ndarray) -> list[float | None]:
    """Each CATH domain's normalized cut against the rest of the chain's judged residues, on the contacts
    residue_domains weighs (CA pairs within DOMAIN_CONTACT, DOMAIN_SEPARATION apart or more, residue by residue)."""
    from scipy.spatial import cKDTree

    pairs = cKDTree(ca).query_pairs(features.DOMAIN_CONTACT, output_type="ndarray") if len(ca) > 1 else None
    if pairs is None or not len(pairs):
        return [None] * (int(ref.max()) + 1 if (ref >= 0).any() else 0)
    pairs = pairs[np.abs(pairs[:, 0] - pairs[:, 1]) >= features.DOMAIN_SEPARATION]
    pairs = pairs[(ref[pairs[:, 0]] >= 0) & (ref[pairs[:, 1]] >= 0)]
    la, lb = ref[pairs[:, 0]], ref[pairs[:, 1]]
    out = []
    for d in range(int(ref.max()) + 1):
        ina, inb = la == d, lb == d
        across = float(np.sum(ina != inb))
        vol = float(np.sum(ina) + np.sum(inb))
        rest = 2.0 * len(pairs) - vol
        out.append(round(across / vol + across / rest, 4) if vol > 0 and rest > 0 else None)
    return out


def _runs(labels: np.ndarray) -> dict[int, int]:
    """Each CATH domain's number of runs along the judged residues."""
    judged = labels[labels >= 0]
    out: dict[int, int] = {}
    for k, d in enumerate(judged):
        if k == 0 or judged[k - 1] != d:
            out[int(d)] = out.get(int(d), 0) + 1
    return out


def _boundary_status(ref: np.ndarray, auto: np.ndarray) -> dict:
    """Each CATH boundary found (an automatic one within TOLERANCE), moved (none within TOLERANCE, but an automatic
    boundary no CATH one claims lies between its two neighbouring CATH boundaries) or missed; the automatic
    boundaries no CATH one claims left over are extra."""
    rb, ab = ds.boundaries(ref, ref >= 0), ds.boundaries(auto)
    free = [x for x in ab if not any(abs(b - x) <= ds.TOLERANCE for b in rb)]
    found = moved = missed = 0
    for k, b in enumerate(rb):
        if any(abs(b - x) <= ds.TOLERANCE for x in ab):
            found += 1
            continue
        lo = rb[k - 1] if k else -1.0
        hi = rb[k + 1] if k + 1 < len(rb) else float(len(ref))
        near = [x for x in free if lo < x < hi]
        if near:
            free.remove(min(near, key=lambda x: abs(x - b)))
            moved += 1
        else:
            missed += 1
    return {"found": found, "moved": moved, "missed": missed, "extra": len(free)}


@functools.cache
def step_lines() -> dict[str, int]:
    """Each _chain_domains step's line in features.py (STEPS: the first line after `def _chain_domains` holding its
    text)."""
    text = Path(features.__file__).read_text(encoding="utf-8").splitlines()
    start = next(k for k, line in enumerate(text) if line.startswith("def _chain_domains("))
    out = {}
    for step, marker in STEPS:
        out[step] = next(k + 1 for k in range(start, len(text)) if marker in text[k])
    return out


def trace(inputs) -> list[dict]:
    """_chain_domains' steps for a one-chain `inputs`, the code of features.py (its fold and βα row merge copied, its
    other steps called), with the domain count after each: [{step, count, line}]; the domains it ends with are
    checked against _chain_domains'."""
    import networkx as nx

    sses, contacts, together = inputs.sses, inputs.contacts, inputs.rows
    graph = nx.Graph()
    graph.add_nodes_from(s.id for s in sses)
    for key, n in contacts.items():
        a, b = sorted(key)
        if a in graph and b in graph:
            graph.add_edge(a, b, weight=float(n))
    for a, b in zip(sses, sses[1:]):
        w = graph.get_edge_data(a.id, b.id, {}).get("weight", 0.0)
        graph.add_edge(a.id, b.id, weight=w + 1.0)
    size = {s.id: len(s) for s in sses}

    def fold(groups):  # features._chain_domains' fold
        groups = list(groups)
        while len(groups) > 1:
            small = min(groups, key=lambda g: sum(size[k] for k in g))
            if sum(size[k] for k in small) >= features.DOMAIN_MIN_SSE_RESIDUES:
                break
            groups.remove(small)
            link = [sum(graph.get_edge_data(a, b, {}).get("weight", 0.0) for a in small for b in g) for g in groups]
            groups[int(np.argmax(link))] |= small
        return groups

    index = list(range(len(inputs.labels)))
    span = {s.id: range(s.start, s.end + 1) for s in sses}
    keep = [[k for i in row if i in span for k in span[i]] for row in together]
    label = dict(zip(index, features.residue_domains(np.asarray(inputs.ca), keep).tolist()))
    lines = step_lines()
    out = [{"step": "residue_domains", "count": len(set(label.values())), "line": None}]
    mine = []
    for s in sses:
        votes = Counter(label[k] for k in span[s.id] if k in label)
        mine.append(min(votes, key=lambda d: (-votes[d], d)) if votes else 0)
    out.append({"step": STEPS[0][0], "count": len(set(mine)), "line": lines[STEPS[0][0]]})
    strays = features._no_strays(mine, [len(s) for s in sses])
    out.append({"step": STEPS[1][0], "count": len(set(strays)), "line": lines[STEPS[1][0]]})
    found: dict[int, set[str]] = {}
    for s, d in zip(sses, strays):
        found.setdefault(d, set()).add(s.id)
    groups = fold(found.values())
    out.append({"step": STEPS[2][0], "count": len(groups), "line": lines[STEPS[2][0]]})
    if len(groups) >= 2:
        kind = {s.id: s.kind for s in sses}
        for row in together:
            strands = {k for k in row if kind.get(k) == "E"}
            if not strands:
                continue
            mine_g = [g for g in groups if g & strands]
            rest = [g - set(row) for g in groups if not g & strands]
            groups = [g for g in [*rest, set().union(*mine_g, set(row) & set(kind))] if g]
        out.append({"step": STEPS[3][0], "count": len(groups), "line": lines[STEPS[3][0]]})
        groups = fold(groups)
        out.append({"step": STEPS[4][0], "count": len(groups), "line": lines[STEPS[4][0]]})
    made = features._chain_domains(sses, contacts, together, (np.asarray(inputs.ca), index))
    ended = sorted(sorted(g) for g in groups) if len(groups) >= 2 else []
    if ended != sorted(made):
        raise AssertionError("the trace of _chain_domains does not end where _chain_domains does")
    return out


def _undone_by(steps: list[dict], n_ref: int) -> dict | None:
    """The first step of a trace that brings the domain count below n_ref, when the residue cut gave enough."""
    if steps[0]["count"] < n_ref:
        return None
    for s in steps[1:]:
        if s["count"] < n_ref:
            return {"step": s["step"], "line": f"src/foldmap/features.py:{s['line']}"}
    return {"step": "after _chain_domains (auto_domains, or the scoring's labels)", "line": None}


def _mechanism(row: dict) -> str:
    """The mechanism in features.py a failure's class points to, from what the cut did."""
    lim = f"DOMAIN_NCUT {BASE_NCUT}"
    cath = [x for x in row["cath_ncut"] if x is not None]
    if row["class"] == "boundary missed":
        if row["undone_by"]:
            u = row["undone_by"]
            where = f" ({u['line']})" if u["line"] else ""
            return (
                f"the residue cut was made ({row['trace'][0]['count']} residue domains); {u['step']}{where} brought "
                f"the count below CATH's {row['n_ref']}"
            )
        if row["tags"].count("tightly packed"):
            return f"not cut: a CATH domain's own normalized cut is {max(cath):.2f}, at or above {lim}"
        if any(c["ncut"] is None for c in row["leaves"]):
            return "not cut: no cut keeps the sizes (DOMAIN_MIN_RESIDUES, DOMAIN_MIN_SEGMENT)"
        if max(row["runs"].values(), default=1) > 2:
            return "not cut: CATH's domain takes more runs than one cut can (two at most, _lowest_cut)"
        return f"not cut: the part's lowest cut was at or above {lim}"
    if row["class"] == "extra split":
        made = [c for c in row["splits"] if not c["right"]]
        said = ", ".join(f"{c['confidence']:.2f}" for c in made) or "none kept"
        return f"a cut below {lim} inside a CATH domain (its confidence: {said})"
    if row["residue_level"]["meets"]:
        return (
            "the residue cut agrees with CATH; the element rules moved it (each element in the domain holding most "
            "of its residues, _no_strays, the folds of _chain_domains)"
        )
    if "discontinuous domain" in row["tags"]:
        return "a discontinuous CATH domain: one cut takes two runs at most (_lowest_cut), so its pieces part"
    found = [c["ncut"] for c in row["splits"]]
    best = f"{min(found):.2f}" if found else "-"
    return (
        f"the lowest cut lies elsewhere: foldmap's cut at normalized cut {best}, CATH's domains at "
        f"{', '.join(f'{x:.2f}' for x in cath) or '-'} (the criterion prefers another boundary)"
    )


def analyse(inputs, entry: dict, deep: bool, against=None) -> dict:
    """One chain: its scores at each t of GRID (features.DOMAIN_NCUT = BASE_NCUT * (1 - t)), its splits with their
    confidence and whether CATH agrees (at t = 0), whether its domains are those of the `against` module (features.py
    at another commit), and, `deep`, what the diagnosis reads."""
    chain = entry["chain"]
    ref = ds.reference_labels([lab.seq for lab in inputs.labels if lab.chain == chain], entry["domains"])
    out = {"pdb": entry["pdb"], "chain": chain, "residues": len(inputs.labels), "grid": {}}
    try:
        for t in GRID:
            features.DOMAIN_NCUT = BASE_NCUT * (1 - t)
            auto = ds.layout_labels(_domains(inputs), inputs.sses, inputs.labels, chain)
            out["grid"][str(t)] = _score(ref, auto)
    finally:
        features.DOMAIN_NCUT = BASE_NCUT
    splits: list = []
    domains = _domains(inputs, splits)
    if _domains(inputs) != domains:
        raise AssertionError(f"{entry['pdb']} {chain}: the domains differ with the record of splits")
    if against is not None:
        out["same_as_against"] = _domains(inputs, module=against) == domains
    auto = ds.layout_labels(domains, inputs.sses, inputs.labels, chain)
    rb = ds.boundaries(ref, ref >= 0)
    index = list(range(len(inputs.labels)))  # the inputs hold one chain: backbone indices are chain positions
    cuts: list = []
    features._chain_domains(inputs.sses, inputs.contacts, inputs.rows, (inputs.ca[index], index), cuts)
    out["splits"] = []
    said = set()
    for c in cuts:  # each kept split beside its cut's sides (features.domain_splits, one cut at a time, said once)
        if "error" in c or not c["made"]:
            continue
        for s in features.domain_splits({chain: [c]}, domains, inputs.sses, inputs.labels):
            if (key := json.dumps(s["between"])) in said:
                continue
            said.add(key)
            meet = _meet_points(c["a"], c["b"])
            right = bool(meet) and all(any(abs(m - b) <= ds.TOLERANCE for b in rb) for m in meet)
            out["splits"].append({**s, "meets_at": meet, "right": right})
    if [{k: v for k, v in s.items() if k not in ("meets_at", "right")} for s in out["splits"]] != splits:
        raise AssertionError(f"{entry['pdb']} {chain}: the splits auto_domains says are not its cuts'")
    out.update(_score(ref, auto))
    if not deep:
        return out
    leaves = [c for c in cuts if not c.get("made") and "error" not in c]
    residue = np.zeros(len(index), dtype=int)
    for d, c in enumerate(sorted(leaves, key=lambda c: int(np.min(c["part"])))):
        residue[np.asarray(c["part"], dtype=int)] = d
    cath = _cath_ncut(np.asarray(inputs.ca, float), ref)
    runs = _runs(ref)
    tags = []
    if any(n >= 2 for n in runs.values()):
        tags.append("discontinuous domain")
    if any(x is not None and x >= BASE_NCUT for x in cath):
        tags.append("tightly packed")
    if out["n_ref"] >= 3 or len(index) > features.DOMAIN_MAX_POINTS:
        tags.append("large")
    steps = trace(inputs)
    out.update(
        cath_ncut=cath,
        runs=runs,
        tags=tags,
        boundaries=_boundary_status(ref, auto),
        residue_level=_score(ref, residue),
        trace=steps,
        cuts=[
            {"made": c["made"], "ncut": None if c["ncut"] is None else round(c["ncut"], 4), "part": len(c["part"])}
            for c in cuts
            if "error" not in c
        ],
        leaves=[{"ncut": None if c["ncut"] is None else round(c["ncut"], 4), "part": len(c["part"])} for c in leaves],
    )
    n_auto, n_ref = out["n_auto"], out["n_ref"]
    out["class"] = CLASSES[0] if n_auto < n_ref else CLASSES[2] if n_auto > n_ref else CLASSES[1]
    out["undone_by"] = _undone_by(steps, n_ref) if out["class"] == CLASSES[0] else None
    out["mechanism"] = _mechanism(out)
    return out


# ------------------------------------------------------------------------------------------------ the sets


def fixture_sets() -> dict[str, list[tuple]]:
    """dom_cal, dom_held and dom_fresh: (inputs, entry) of each chain, from the committed fixtures; dom_held the third
    held out chain by chain (sets.json; scripts/paper/domains.py), dom_cal the rest."""
    entries = ds.load_reference(DOMAINS / "cath_reference.json")
    keys = sorted((e["pdb"], e["chain"]) for e in entries)
    held = set(random.Random(ds.HELD_OUT_SEED).sample(keys, len(keys) // 3))
    reference = ds.load_fixture(DOMAINS / "cath_inputs.npz")
    fresh_entries = json.loads((DOMAINS / "cath_fresh.json").read_text(encoding="utf-8"))["chains"]
    fresh = ds.load_fixture(DOMAINS / "fresh_inputs.npz")

    def key(e):
        return f"{e['pdb']}_{e['chain']}"

    return {
        "dom_cal": [(reference[key(e)], e) for e in entries if (e["pdb"], e["chain"]) not in held],
        "dom_held": [(reference[key(e)], e) for e in entries if (e["pdb"], e["chain"]) in held],
        "dom_fresh": [(fresh[key(e)], e) for e in fresh_entries],
    }


@functools.cache
def _against_module(path: str | None):
    """features.py of another commit (its text at `path`), loaded as foldmap._features_against; None without one."""
    if path is None:
        return None
    spec = importlib.util.spec_from_file_location("foldmap._features_against", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def third_entry(pid: str, against: str | None = None) -> dict:
    """Every chain of the third sample's entry `pid` that CATH chops, analysed (random100.cath_chains); chains that
    could not be scored with the reason, as random100.domain_reference counts them."""
    from random100 import NO_RESIDUE, cath_chains

    from foldmap.fetch import fetch

    out = {"id": pid, "chains": [], "unscored": []}
    try:
        chopped = cath_chains(pid)
        path = fetch(pid) if chopped else None
    except Exception as err:  # noqa: BLE001 (PDBe or RCSB unreachable: said, not scored)
        out["error"] = f"{type(err).__name__}: {str(err)[:160]}"
        return out
    for c in sorted(chopped, key=lambda c: c["chain"]):
        try:
            inputs = ds.chain_inputs(path, c["chain"])
            got = analyse(inputs, c, deep=True, against=_against_module(against))
        except Exception as err:  # noqa: BLE001
            out["unscored"].append({"chain": c["chain"], "reason": type(err).__name__, "detail": str(err)[:160]})
            continue
        if not got["n_ref"]:
            out["unscored"].append({"chain": c["chain"], "reason": NO_RESIDUE, "detail": ""})
            continue
        out["chains"].append({**got, "cath": c["cath"]})
    return out


def third_ids() -> list[str]:
    from random100 import PDB_ID

    sample = json.loads((THIRD / "sample.json").read_text(encoding="utf-8"))
    return [e["id"] for e in sample["entries"] if PDB_ID.fullmatch(e["id"])]


# ------------------------------------------------------------------------------------------------ the record


def _summary(rows: list[dict], t) -> dict:
    got = [r["grid"][str(t)] for r in rows]
    means = ds.summary(got)
    return {
        "chains": len(got),
        "meet": sum(g["meets"] for g in got),
        **{k: round(means[k], 4) for k in ("overlap", "f1", "recall", "precision")},
    }


def _entries_failing(entries: list[dict], t) -> int:
    return sum(any(not c["grid"][str(t)]["meets"] for c in e["chains"]) for e in entries if e["chains"])


def choose(rows: list[dict]) -> float | None:
    """t on dom_cal: the largest of GRID above 0 at which no chain's overlap, F1 or recall is lower than at t = 0;
    None when there is none."""
    ok = []
    for t in GRID[1:]:
        if all(r["grid"][str(t)][k] >= r["grid"]["0.0"][k] - 1e-9 for r in rows for k in ("overlap", "f1", "recall")):
            ok.append(t)
    return max(ok) if ok else None


def choose_most(rows: list[dict]) -> float:
    """t on dom_cal read most favourably: the t of GRID with most chains meeting the targets, the smallest of those
    (0 when none beats it). Measured beside choose()'s, so the record shows what even this reading does elsewhere."""
    return min(GRID, key=lambda t: (-sum(r["grid"][str(t)]["meets"] for r in rows), t))


def _no_worse(now: dict, then: dict) -> bool:
    return now["meet"] >= then["meet"] and all(now[k] >= then[k] - 1e-9 for k in ("overlap", "f1", "recall"))


def _fell(rows: list[dict], now: dict, then: dict, t) -> list[str]:
    """What fell at t against t = 0 on one set, in words: the chains meeting the targets, each mean that fell, and the
    chains whose overlap, F1 or recall fell."""
    out = []
    if now["meet"] < then["meet"]:
        out.append(f"{then['meet']} -> {now['meet']} meet")
    names = {"overlap": "overlap", "f1": "F1", "recall": "recall"}
    out += [f"mean {names[k]} {then[k]:.3f} -> {now[k]:.3f}" for k in names if now[k] < then[k] - 1e-9]
    lowered = [
        f"{r['pdb']} {r['chain']}"
        for r in rows
        if any(r["grid"][str(t)][k] < r["grid"]["0.0"][k] - 1e-9 for k in ("overlap", "f1", "recall"))
    ]
    if lowered:
        out.append(f"lowered: {', '.join(lowered)}")
    return out


def _auc(right: list[float], wrong: list[float]) -> float | None:
    """The chance a right split's confidence is above a wrong one's (ties half): the area under the ROC curve."""
    if not right or not wrong:
        return None
    wins = sum((r > w) + 0.5 * (r == w) for r in right for w in wrong)
    return round(wins / (len(right) * len(wrong)), 3)


def calibration(splits: list[dict]) -> dict:
    """Per group of sets: in each confidence bin, the splits right of all with the Wilson 95% interval of the
    fraction right; the AUC; the medians."""
    wilson = _wilson()
    third = [s for s in splits if s["set"] == "r100_third"]
    first: dict[str, str] = {}  # the third sample's first chain of each entry (copies of one protein said once)
    for s in third:
        first.setdefault(s["pdb"], s["chain"])
    groups = {
        "dom_cal (tuning: DOMAIN_NCUT was chosen on it)": [s for s in splits if s["set"] == "dom_cal"],
        "dom_held + dom_fresh (regression)": [s for s in splits if s["set"] in ("dom_held", "dom_fresh")],
        "r100_third, every chain": third,
        "r100_third, its first chain per entry": [s for s in third if first[s["pdb"]] == s["chain"]],
    }
    out = {}
    for name, group in groups.items():
        bins = []
        for lo, hi in BINS:
            got = [s for s in group if lo <= s["confidence"] < hi]
            k = sum(s["right"] for s in got)
            ci = wilson(k, len(got))
            bins.append(
                {
                    "bin": f"{lo:.2f}-{min(hi, 1.0):.2f}",
                    "right": k,
                    "n": len(got),
                    "wilson95": None if ci is None else [round(ci[0], 2), round(ci[1], 2)],
                }
            )
        right = [s["confidence"] for s in group if s["right"]]
        wrong = [s["confidence"] for s in group if not s["right"]]
        out[name] = {
            "splits": len(group),
            "right": len(right),
            "bins": bins,
            "auc": _auc(right, wrong),
            "median_right": None if not right else round(float(np.median(right)), 2),
            "median_wrong": None if not wrong else round(float(np.median(wrong)), 2),
        }
    low_wrong = [s for s in third if s["confidence"] < BINS[0][1] and not s["right"]]
    top = Counter(s["pdb"] for s in low_wrong).most_common(1)
    out["_third_low_wrong"] = {
        "n": len(low_wrong),
        "top_entry": top[0][0] if top else None,
        "top": top[0][1] if top else 0,
    }
    return out


def build(jobs: int, against: dict | None) -> tuple[dict, dict]:
    """(the record, its per-chain detail)."""
    path = against["path"] if against else None
    module = _against_module(path)
    fixtures = fixture_sets()
    rows = {name: [analyse(i, e, deep=False, against=module) for i, e in pairs] for name, pairs in fixtures.items()}
    ids = third_ids()
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        third = list(pool.map(functools.partial(third_entry, against=path), ids))
    rows["r100_third"] = [c for e in third for c in e["chains"]]
    multi = [c for c in rows["r100_third"] if c["n_ref"] >= 2]
    failures = [c for c in multi if not c["meets"]]
    undone = [f for f in failures if f["undone_by"]]
    counts = {
        "class": {k: sum(f["class"] == k for f in failures) for k in CLASSES},
        "class_entries": {k: len({f["pdb"] for f in failures if f["class"] == k}) for k in CLASSES},
        "tag": {k: sum(k in f["tags"] for f in failures) for k in TAGS},
        "untagged": sum(not f["tags"] for f in failures),
        "boundaries": {k: sum(f["boundaries"][k] for f in failures) for k in ("found", "moved", "missed", "extra")},
        "residue_level_meets": sum(f["residue_level"]["meets"] for f in failures),
        "undone_by": dict(Counter(f"{f['undone_by']['step']} ({f['undone_by']['line']})" for f in undone)),
        "undone_entries": len({f["pdb"] for f in undone}),
        "undone_residue_level_meets": sum(f["residue_level"]["meets"] for f in undone),
    }
    single_fail = [c for c in rows["r100_third"] if c["n_ref"] < 2 and not c["meets"]]
    keep = ("pdb", "chain", "residues", "n_ref", "n_auto", "overlap", "f1", "recall", "precision", "meets")
    lean = (
        "class",
        "tags",
        "cath",
        "cath_ncut",
        "runs",
        "boundaries",
        "residue_level",
        "trace",
        "undone_by",
        "mechanism",
    )
    lean_failures = [{k: f[k] for k in (*keep, *lean) if k in f} for f in failures]
    for lf, f in zip(lean_failures, failures):
        lf["splits"] = [{k: s[k] for k in ("between", "ncut", "confidence", "right")} for s in f["splits"]]
    splits = [
        {"set": name, "pdb": r["pdb"], "chain": r["chain"], **s}
        for name, rs in rows.items()
        for r in rs
        for s in r["splits"]
        if "error" not in s
    ]
    grid = {name: {str(t): _summary(rs, t) for t in GRID} for name, rs in rows.items()}
    for t in GRID:
        grid["r100_third"][str(t)]["entries_failing_s5"] = _entries_failing(third, t)
    chosen, most = choose(rows["dom_cal"]), choose_most(rows["dom_cal"])

    def decide(t: float) -> dict:
        """What t does on every set against t = 0, and whether ruling Q42 would commit it."""
        per = {name: _no_worse(grid[name][str(t)], grid[name]["0.0"]) for name in rows}
        third_ok = grid["r100_third"][str(t)]["entries_failing_s5"] <= grid["r100_third"]["0.0"]["entries_failing_s5"]
        changed = {
            name: [
                {"pdb": r["pdb"], "chain": r["chain"], "at_0": r["grid"]["0.0"], "at_t": r["grid"][str(t)]}
                for r in rs
                if r["grid"]["0.0"] != r["grid"][str(t)]
            ]
            for name, rs in rows.items()
        }
        fell = {name: _fell(rs, grid[name][str(t)], grid[name]["0.0"], t) for name, rs in rows.items()}
        return {
            "t": t,
            "no_worse": per,
            "third_entries_no_worse": third_ok,
            "commit": all(per.values()) and third_ok,
            "fell": fell,
            "changed": changed,
        }

    # the t chosen, the most favourable reading and the smallest step: what fell on each set is said at each
    decisions = {str(t): decide(t) for t in dict.fromkeys(x for x in (chosen, most, GRID[1]) if x)}
    best_third = min(GRID, key=lambda t: (-grid["r100_third"][str(t)]["meet"], t))
    gained = [r for r in rows["r100_third"] if r["grid"][str(best_third)]["meets"] and not r["grid"]["0.0"]["meets"]]
    top = Counter(r["pdb"] for r in gained).most_common(1)
    third_gain = {
        "t": best_third,
        "meet": grid["r100_third"][str(best_third)]["meet"],
        "meet_at_0": grid["r100_third"]["0.0"]["meet"],
        "gained": len(gained),
        "top_entry": top[0][0] if top else None,
        "top": top[0][1] if top else 0,
        "ts_gaining": [t for t in GRID if grid["r100_third"][str(t)]["meet"] > grid["r100_third"]["0.0"]["meet"]],
    }
    compared = None
    if against:
        every = [r for rs in rows.values() for r in rs]
        compared = {
            "against": against["rev"],
            "source_sha256": against["sha256"],
            "chains": len(every),
            "differ": sum(not r["same_as_against"] for r in every),
            "differing": [f"{r['pdb']} {r['chain']}" for r in every if not r["same_as_against"]],
        }
    doc = {
        "sets": set_tags(),
        "set_sizes": {name: len(rs) for name, rs in rows.items()},
        "targets": ds.TARGETS,
        "tolerance": ds.TOLERANCE,
        "domain_ncut": BASE_NCUT,
        "grid": list(GRID),
        "against": compared,
        "third": {
            "entries": len(ids),
            "entries_scored": sum(bool(e["chains"]) for e in third),
            "chains": len(rows["r100_third"]),
            "multi": len(multi),
            "multi_meet": len(multi) - len(failures),
            "single_failing": [f"{c['pdb']} {c['chain']}" for c in single_fail],
            "unscored": [{"id": e["id"], **u} for e in third for u in e["unscored"]],
            "errors": [{"id": e["id"], "error": e["error"]} for e in third if "error" in e],
        },
        "diagnosis": {"classes": list(CLASSES), "tags": list(TAGS), "counts": counts, "failures": lean_failures},
        "splits": splits,
        "calibration": calibration(splits),
        "undivided": {"grid": grid, "chosen": chosen, "most": most, "decisions": decisions, "third_gain": third_gain},
    }
    detail = {name: [{k: v for k, v in r.items() if k != "splits"} for r in rs] for name, rs in rows.items()}
    return doc, detail


def _count(dg: dict, cls: str) -> str:
    """A class's failing chains, and the entries they are in."""
    got = [f for f in dg["failures"] if f["class"] == cls]
    return f"{len(got)} ({len({f['pdb'] for f in got})})"


def _ci(b: dict) -> str:
    if not b["n"]:
        return "-"
    lo, hi = b["wilson95"]
    return f"{b['right']}/{b['n']} [{lo:.2f}, {hi:.2f}]"


def render(doc: dict) -> str:
    th, dg, un, cal = doc["third"], doc["diagnosis"], doc["undivided"], doc["calibration"]
    c = dg["counts"]
    lines = [
        "# The automatic domains with their uncertainty",
        "",
        f"Plan 28 Task 5 (the review's 3.4, rulings Q42, Q43). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "Sets, tagged as scripts/paper/sets.json tags them: "
        + "; ".join(f"{k}, {doc['set_sizes'][k]} chains ({v})" for k, v in doc["sets"].items())
        + ". dom_cal and dom_held are the reference chains split chain by chain (seed 20261016:domains), as the "
        "parameters were chosen. A chain meets the targets (S5) when its overlap, boundary F1 and recall are each at "
        f"least {doc['targets']['overlap']}, {doc['targets']['f1']} and {doc['targets']['recall']} against CATH "
        f"(domain_score, boundaries within {doc['tolerance']} residues). DOMAIN_NCUT is {doc['domain_ncut']}. "
        "report-chains.json holds each chain's cuts, trace and scores at every t.",
    ]
    if doc["against"]:
        a = doc["against"]
        lines += [
            "",
            f"Domains against features.py at {a['against']} (sha256 {a['source_sha256'][:12]}): {a['chains']} chains "
            f"compared (every chain of every set), {a['differ']} differ"
            + (f" ({', '.join(a['differing'])})" if a["differing"] else "")
            + ".",
        ]
    lines += [
        "",
        "## 1. Why the multi-domain chains of the third sample fail (read-only)",
        "",
        f"{th['entries']} PDB entries in the sample; {th['entries_scored']} with a chain CATH chops scored; "
        f"{th['chains']} chains scored, {th['multi']} of them with two CATH domains or more: {th['multi_meet']} meet "
        f"the targets, {len(dg['failures'])} fail, in {len({f['pdb'] for f in dg['failures']})} entries (copies of "
        "one protein in an entry count once each: 1GQ2 has 16 chains). Single-domain chains failing (an extra "
        f"split): {len(th['single_failing'])}"
        + (f" ({', '.join(th['single_failing'])})" if th["single_failing"] else "")
        + ".",
        "",
        "What went wrong (exclusive, by the automatic and CATH domain counts):",
        "",
        "| class | chains (entries) | what it means | the mechanism in features.py |",
        "|---|---|---|---|",
        f"| boundary missed | {_count(dg, 'boundary missed')} | fewer automatic domains than CATH's | "
        "the part's lowest cut at DOMAIN_NCUT or above (residue_domains, _best_cut), or no cut keeps "
        "DOMAIN_MIN_RESIDUES / DOMAIN_MIN_SEGMENT, or a cut made and undone by a step of _chain_domains (traced "
        "below) |",
        f"| boundary moved | {_count(dg, 'boundary moved')} | as many, with a boundary elsewhere | "
        "the lowest normalized cut is not CATH's boundary (the criterion prefers another), a discontinuous CATH "
        "domain one cut cannot take (two runs at most, _lowest_cut), or the element rules moved a right residue cut |",
        f"| extra split | {_count(dg, 'extra split')} | more automatic domains than CATH's | a cut below "
        "DOMAIN_NCUT inside one CATH domain (a lobe or subdomain CATH keeps whole) |",
        "",
        f"Missed boundaries whose residue cut was made: {sum(c['undone_by'].values())} chains in "
        f"{c['undone_entries']} entries, {c['undone_residue_level_meets']} of them meeting the targets at the residue "
        "level; the first step of _chain_domains that brings the count below CATH's: "
        + "; ".join(f"{k} {v}" for k, v in c["undone_by"].items())
        + ".",
        "",
        "What the failing chains are like (tags, not exclusive): "
        + ", ".join(f"{k} {v}" for k, v in c["tag"].items())
        + f", none of them {c['untagged']}. Discontinuous: a CATH domain in two runs or more along the "
        "chain; tightly packed: a CATH domain whose normalized cut against the rest of the chain (the contacts "
        "residue_domains weighs) is at DOMAIN_NCUT or above; large: three CATH domains or more, or more residues "
        "than DOMAIN_MAX_POINTS. Their CATH boundaries: "
        + ", ".join(f"{k} {v}" for k, v in c["boundaries"].items())
        + " (moved: none within the tolerance, but an automatic boundary no CATH one claims lies between its two "
        "neighbours; extra: automatic boundaries no CATH one claims, left over). "
        f"{c['residue_level_meets']} of the failures meet the targets at the residue level, before the "
        "element rules.",
        "",
        "| id | chain | residues | CATH | auto | overlap | F1 | recall | class | tags | boundaries f/m/x/e | "
        "CATH domains' ncut | mechanism |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for f in sorted(dg["failures"], key=lambda f: (CLASSES.index(f["class"]), f["pdb"], f["chain"])):
        b = f["boundaries"]
        cath = ", ".join("-" if x is None else f"{x:.2f}" for x in f["cath_ncut"])
        lines.append(
            f"| {f['pdb']} | {f['chain']} | {f['residues']} | {f['n_ref']} | {f['n_auto']} | {f['overlap']:.2f} | "
            f"{f['f1']:.2f} | {f['recall']:.2f} | {f['class']} | {', '.join(f['tags']) or '-'} | "
            f"{b['found']}/{b['moved']}/{b['missed']}/{b['extra']} | {cath} | {f['mechanism']} |"
        )
    names = [k for k in cal if not k.startswith("_")]
    low = cal["_third_low_wrong"]
    lines += [
        "",
        "## 2. The confidence of each split",
        "",
        "A split's confidence is 1 - (its normalized cut / DOMAIN_NCUT): 0 at the limit, where the part would have "
        "stayed whole, 1 for a part with no contact across (features.split_confidence, said by domain_confidence in "
        "the notes, `foldmap check`, its --json report, the MCP result and the app; never drawn). It is the margin "
        "below a tuned limit, not a probability, and moves with that limit. A split is right here when each place "
        f"its two sides meet lies within {doc['tolerance']} residues of a CATH boundary. Per group: right of all in "
        "each bin, with the Wilson 95% interval of the fraction right (scripts/paper/common.wilson); the AUC is the "
        "chance a right split's confidence is above a wrong one's.",
        "",
        "| confidence | " + " | ".join(names) + " |",
        "|---|" + "---|" * len(names),
    ]
    for k, (lo, hi) in enumerate(BINS):
        lines.append(f"| {lo:.2f}-{min(hi, 1.0):.2f} | " + " | ".join(_ci(cal[n]["bins"][k]) for n in names) + " |")
    lines.append(
        "| AUC (right vs wrong) | "
        + " | ".join("-" if cal[n]["auc"] is None else f"{cal[n]['auc']:.2f}" for n in names)
        + " |"
    )
    lines.append(
        "| splits (right); median confidence right / wrong | "
        + " | ".join(
            f"{cal[n]['splits']} ({cal[n]['right']}); {cal[n]['median_right']} / {cal[n]['median_wrong']}"
            for n in names
        )
        + " |"
    )
    lines += [
        "",
        f"On the third sample, {low['top']} of the {low['n']} wrong splits under 0.25 are {low['top_entry']}'s "
        "(copies of one protein); counted once per entry (its first chain) the bins hold a few splits each. dom_cal "
        "is where DOMAIN_NCUT was chosen, so its column is in-sample; the regression column and the per-entry column "
        "are the independent views.",
    ]
    lines += [
        "",
        "## 3. Keeping a weakly supported region undivided",
        "",
        "The reading measured: a region whose cut has confidence t or less is not cut, and never weighed again, so "
        'its stronger sub-cuts go with it (the literal "kept undivided"). That is a cut made only below DOMAIN_NCUT '
        "x (1 - t), as the confidence is the margin below that limit (up to its rounding to 0.01), measured by "
        "setting features.DOMAIN_NCUT so. Merging only across a weak drawn split, keeping the cuts inside it, is "
        "another reading, not measured. For each t on every set (chains meeting the targets of all chains; mean "
        "overlap, F1, recall):",
        "",
        "| t | DOMAIN_NCUT x (1 - t) | "
        + " | ".join(f"{name} meet | {name} overlap / F1 / recall" for name in doc["sets"])
        + " | third: entries failing S5 |",
        "|---|---|" + "---|---|" * len(doc["sets"]) + "---|",
    ]
    for t in doc["grid"]:
        cells = []
        for name in doc["sets"]:
            g = un["grid"][name][str(t)]
            cells += [f"{g['meet']}/{g['chains']}", f"{g['overlap']:.3f} / {g['f1']:.3f} / {g['recall']:.3f}"]
        lines.append(
            f"| {t:.2f} | {doc['domain_ncut'] * (1 - t):.4f} | "
            + " | ".join(cells)
            + f" | {un['grid']['r100_third'][str(t)]['entries_failing_s5']} |"
        )
    lines.append("")
    rule = "the largest t at which no dom_cal chain's overlap, F1 or recall drops"
    if un["chosen"] is None:
        lines.append(f"Chosen on dom_cal ({rule}): none; every t above 0 lowers some dom_cal chain.")
    else:
        lines.append(f"Chosen on dom_cal ({rule}): t = {un['chosen']:.2f}.")
    if un["most"]:
        lines.append(
            f"Read most favourably (the t with most dom_cal chains meeting the targets, the smallest of them): "
            f"t = {un['most']:.2f}, measured on every set too, to show what even that reading does elsewhere."
        )
    else:
        lines.append(
            "Read most favourably (the t with most dom_cal chains meeting the targets): t = 0, as drawn; no t above "
            "0 has as many."
        )
    lines.append(
        f"What each set loses is said below at the smallest step, t = {GRID[1]:.2f}, and at any t a reading chose."
    )
    for key, dec in un["decisions"].items():
        verdict = "no set is worse" if dec["commit"] else "some set is worse"
        lines += ["", f'At t = {float(key):.2f}: {verdict}. What fell, per set ("worse" is any of these):', ""]
        for name, fell in dec["fell"].items():
            lines.append(f"- {name}: {'; '.join(fell) if fell else 'nothing fell'}")
        if not dec["third_entries_no_worse"]:
            lines.append("- r100_third: more entries fail S5")
        lines += ["", "Chains whose scores change:", ""]
        for name, got in dec["changed"].items():
            for ch in got:
                a, b = ch["at_0"], ch["at_t"]
                lines.append(
                    f"- {name} {ch['pdb']} {ch['chain']}: domains {a['n_auto']} -> {b['n_auto']}, overlap "
                    f"{a['overlap']:.2f} -> {b['overlap']:.2f}, F1 {a['f1']:.2f} -> {b['f1']:.2f}, recall "
                    f"{a['recall']:.2f} -> {b['recall']:.2f}, meets {a['meets']} -> {b['meets']}"
                )
    g = un["third_gain"]
    if g["ts_gaining"]:
        lines += [
            "",
            f"The third sample alone would gain: more of its chains meet the targets at t = "
            f"{', '.join(f'{t:.2f}' for t in g['ts_gaining'])}, most at t = {g['t']:.2f} ({g['meet_at_0']} -> "
            f"{g['meet']}), where {g['top']} of its {g['gained']} new passes are {g['top_entry']}'s copies. It cannot "
            "choose t: it is no tuning set (sets.json), and one entry's copies decide its count.",
        ]
    commit = any(d["commit"] for d in un["decisions"].values() if float(d["t"]) == un["chosen"])
    lines += [
        "",
        "Ruling Q42 commits the change only when the t chosen on dom_cal leaves no set worse: "
        + ("it does, and the change is committed." if commit else "it is not committed; the confidence ships alone."),
    ]
    if th["unscored"] or th["errors"]:
        lines += [
            "",
            "Not scored: "
            + "; ".join(
                f"{u['id']} {u.get('chain', '')} ({u.get('reason', u.get('error'))})"
                for u in [*th["unscored"], *th["errors"]]
            ),
        ]
    lines += [
        "",
        "## Limits",
        "",
        "- The third sample has been run and looked at since plan 15; its failures are described here, never tuned "
        "on. dom_cal is a tuning set, dom_held and dom_fresh regression sets (sets.json).",
        "- CATH is one chopping; another (SCOP, ECOD) may draw some boundaries elsewhere. A wrong split here is one "
        "CATH does not draw.",
        "- The classes are read from domain counts and boundaries; a chain's mechanism is the first that its cut "
        "record explains, not every one that may contribute.",
        "- The confidence measures how far a cut is from not being made; it says nothing of where else a boundary "
        "could lie. It moves with DOMAIN_NCUT, a sensitive constant (docs/validation/sensitivity: x1.25 moves the "
        "third sample's S5 failures from 12 to 19).",
        "",
    ]
    return "\n".join(lines)


def _against(rev: str | None, source: Path | None) -> dict | None:
    """--against: the revision, the file holding its features.py (from git show, or --against-source) and its
    SHA-256."""
    if rev is None:
        return None
    if source is None:
        text = subprocess.run(
            ["git", "-C", str(ROOT), "show", f"{rev}:src/foldmap/features.py"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        source = Path(tempfile.mkdtemp()) / f"features_{rev}.py"  # removed after the run
        source.write_text(text, encoding="utf-8")
    data = Path(source).read_bytes()
    return {"rev": rev, "path": str(Path(source).resolve()), "sha256": hashlib.sha256(data).hexdigest()}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--against", help="a commit whose features.py the domains are compared with")
    parser.add_argument("--against-source", type=Path, help="that commit's features.py, where there is no git")
    args = parser.parse_args(argv)
    against = _against(args.against, args.against_source)
    record, detail = build(args.jobs, against)
    doc = {"source_stamp": SOURCE_STAMP, **record}  # the stamp taken as this process started
    chains = {"source_stamp": SOURCE_STAMP, "chains": detail}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8", newline="\n")
    (args.out / "report-chains.json").write_text(json.dumps(chains, default=str) + "\n", encoding="utf-8", newline="\n")
    (args.out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    if args.against and args.against_source is None:
        Path(against["path"]).unlink(missing_ok=True)
    print(f"wrote {args.out / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
