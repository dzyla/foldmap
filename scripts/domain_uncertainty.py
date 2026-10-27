# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The record of the automatic domains' uncertainty (plan 28 Task 5, the review's 3.4, ruling Q42): why they fail
against CATH, what each split's confidence is worth, and what keeping a weakly supported region undivided would do.

    python scripts/domain_uncertainty.py [--out docs/validation/domain-uncertainty] [--jobs 12]

Writes report.json and report.md, stamped with the code that measured (foldmap.source_stamp). Reads the structures of
the third random sample from ~/.cache/foldmap (fetched when missing) and CATH's choppings through foldmap.cath (cached
in ~/.cache/foldmap/cath); the domain sets from their committed fixtures. Not run in CI.

Sets (scripts/paper/sets.json): dom_cal (tuning: the 22 calibration chains), dom_held (regression: the 11 held out by
homology group), dom_fresh (regression: 11 chains of superfamilies none of the 33 has), all from tests/data/domains;
r100_third (the third random sample: every chain CATH chops, scored as random100.domain_reference scores it; looked
at, not unseen). Every chain is scored by domain_score (overlap, boundary recall, precision and F1 against CATH,
TOLERANCE residues) and meets the targets when its overlap, F1 and recall are each at least domain_score.TARGETS (S5,
science.domain_findings).

1. diagnosis (read-only): each chain of the third sample with two CATH domains or more that fails S5, classed by
   what went wrong (exclusive: fewer automatic domains than CATH's, "boundary missed"; as many, "boundary moved"; more,
   "extra split") and tagged by what the chain is like ("discontinuous domain": a CATH domain in two runs or more
   along the chain; "tightly packed": a CATH domain whose normalized cut against the rest of the chain, on the
   contacts residue_domains weighs, is at DOMAIN_NCUT or above, so the criterion cannot see it; "large": three CATH
   domains or more, or more residues than DOMAIN_MAX_POINTS), with the mechanism in features.py named from what the
   cut did (its record: features.residue_domains' `cuts`), each boundary found / moved / missed / extra, and whether
   the residue-level domains (before the element rules) met the targets.
2. splits: every inferred split the drawn domains keep (features.domain_splits) on every set, with its confidence,
   called right when each place its two sides meet lies within TOLERANCE residues of a CATH boundary.
3. undivided: keeping a region undivided when its split's confidence is below t is the same as making a cut only
   below DOMAIN_NCUT * (1 - t), as the confidence is the margin below that limit; measured by setting
   features.DOMAIN_NCUT so, for t on GRID, on every set. t is chosen on dom_cal (tuning) only: the largest t at which
   no dom_cal chain's overlap, F1 or recall is lower than at t = 0. It is then measured on dom_held, dom_fresh and the
   third sample; ruling Q42 commits the change only when no set is worse (chains meeting the targets, mean overlap,
   F1 and recall; the third sample's entries failing S5 too)."""

from __future__ import annotations

import argparse
import json
import sys
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
GRID = tuple(round(0.05 * k, 2) for k in range(11))  # t: 0 (as drawn) to 0.5
BASE_NCUT = features.DOMAIN_NCUT
CLASSES = ("boundary missed", "boundary moved", "extra split")
TAGS = ("discontinuous domain", "tightly packed", "large")
BINS = ((0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.01))
SETS = {
    "dom_cal": "tuning",
    "dom_held": "regression",
    "dom_fresh": "regression",
    "r100_third": "looked at (third random sample; not unseen)",
}


# ------------------------------------------------------------------------------------------------ one chain


def _meets(s: dict) -> bool:
    return all(s[k] >= ds.TARGETS[k] for k in ("overlap", "f1", "recall"))


def _score(ref: np.ndarray, auto: np.ndarray) -> dict:
    got = ds.score(ref, auto)
    return {**{k: round(v, 4) if isinstance(v, float) else v for k, v in got.items()}, "meets": _meets(got)}


def _domains(inputs, splits=None):
    bb = SimpleNamespace(labels=inputs.labels, ca=inputs.ca)
    return features.auto_domains(inputs.sses, inputs.contacts, inputs.rows, {}, frozenset(), bb=bb, splits=splits)


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


def _mechanism(row: dict) -> str:
    """The mechanism in features.py a failure's class points to, from what the cut did."""
    lim = f"DOMAIN_NCUT {BASE_NCUT}"
    cath = [x for x in row["cath_ncut"] if x is not None]
    if row["class"] == "boundary missed":
        if row["residue_level"]["n_auto"] >= row["n_ref"]:
            return (
                "the residue cut was made, and the element rules undid it (_chain_domains: a domain under "
                "DOMAIN_MIN_SSE_RESIDUES folded, a βα row kept whole, _no_strays)"
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


def analyse(inputs, entry: dict, deep: bool) -> dict:
    """One chain: its scores at each t of GRID (features.DOMAIN_NCUT = BASE_NCUT * (1 - t)), its splits with their
    confidence and whether CATH agrees (at t = 0), and, `deep`, what the diagnosis reads."""
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
    out.update(
        cath_ncut=cath,
        runs=runs,
        tags=tags,
        boundaries=_boundary_status(ref, auto),
        residue_level=_score(ref, residue),
        cuts=[
            {"made": c["made"], "ncut": None if c["ncut"] is None else round(c["ncut"], 4), "part": len(c["part"])}
            for c in cuts
            if "error" not in c
        ],
        leaves=[{"ncut": None if c["ncut"] is None else round(c["ncut"], 4), "part": len(c["part"])} for c in leaves],
    )
    n_auto, n_ref = out["n_auto"], out["n_ref"]
    out["class"] = CLASSES[0] if n_auto < n_ref else CLASSES[2] if n_auto > n_ref else CLASSES[1]
    out["mechanism"] = _mechanism(out)
    return out


# ------------------------------------------------------------------------------------------------ the sets


def fixture_sets() -> dict[str, list[tuple]]:
    """dom_cal, dom_held and dom_fresh: (inputs, entry) of each chain, from the committed fixtures."""
    entries = ds.load_reference(DOMAINS / "cath_reference.json")
    held = ds.held_out(entries)
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


def third_entry(pid: str) -> dict:
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
            got = analyse(inputs, c, deep=True)
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


def build(jobs: int) -> dict:
    fixtures = fixture_sets()
    rows = {name: [analyse(i, e, deep=False) for i, e in pairs] for name, pairs in fixtures.items()}
    ids = third_ids()
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        third = list(pool.map(third_entry, ids))
    rows["r100_third"] = [c for e in third for c in e["chains"]]
    multi = [c for c in rows["r100_third"] if c["n_ref"] >= 2]
    failures = [c for c in multi if not c["meets"]]
    counts = {
        "class": {k: sum(f["class"] == k for f in failures) for k in CLASSES},
        "class_entries": {k: len({f["pdb"] for f in failures if f["class"] == k}) for k in CLASSES},
        "tag": {k: sum(k in f["tags"] for f in failures) for k in TAGS},
        "untagged": sum(not f["tags"] for f in failures),
        "boundaries": {k: sum(f["boundaries"][k] for f in failures) for k in ("found", "moved", "missed", "extra")},
        "residue_level_meets": sum(f["residue_level"]["meets"] for f in failures),
    }
    single_fail = [c for c in rows["r100_third"] if c["n_ref"] < 2 and not c["meets"]]
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
        return {
            "t": t,
            "no_worse": per,
            "third_entries_no_worse": third_ok,
            "commit": all(per.values()) and third_ok,
            "changed": changed,
        }

    decisions = {str(t): decide(t) for t in dict.fromkeys(x for x in (chosen, most) if x)}
    return {
        "sets": SETS,
        "targets": ds.TARGETS,
        "tolerance": ds.TOLERANCE,
        "domain_ncut": BASE_NCUT,
        "grid": list(GRID),
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
        "diagnosis": {"classes": list(CLASSES), "tags": list(TAGS), "counts": counts, "failures": failures},
        "splits": splits,
        "undivided": {"grid": grid, "chosen": chosen, "most": most, "decisions": decisions},
        "chains": {name: [{k: v for k, v in r.items() if k != "splits"} for r in rs] for name, rs in rows.items()},
    }


def _count(dg: dict, cls: str) -> str:
    """A class's failing chains, and the entries they are in."""
    got = [f for f in dg["failures"] if f["class"] == cls]
    return f"{len(got)} ({len({f['pdb'] for f in got})})"


def _bins(splits: list[dict]) -> list[tuple[str, int, int]]:
    out = []
    for lo, hi in BINS:
        got = [s for s in splits if lo <= s["confidence"] < hi]
        out.append((f"{lo:.2f}-{min(hi, 1.0):.2f}", sum(s["right"] for s in got), sum(not s["right"] for s in got)))
    return out


def render(doc: dict) -> str:
    th, dg, un = doc["third"], doc["diagnosis"], doc["undivided"]
    lines = [
        "# The automatic domains with their uncertainty",
        "",
        f"Plan 28 Task 5 (the review's 3.4, ruling Q42). Written by: {stamp_line(doc['source_stamp'])}",
        "",
        "Sets (scripts/paper/sets.json): "
        + "; ".join(f"{k} ({v})" for k, v in doc["sets"].items())
        + f". A chain meets the targets (S5) when its overlap, boundary F1 and recall are each at least "
        f"{doc['targets']['overlap']}, {doc['targets']['f1']} and {doc['targets']['recall']} against CATH "
        f"(domain_score, boundaries within {doc['tolerance']} residues). DOMAIN_NCUT is {doc['domain_ncut']}.",
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
        "DOMAIN_MIN_RESIDUES / DOMAIN_MIN_SEGMENT, or a cut made and undone by the element rules (_chain_domains: "
        "DOMAIN_MIN_SSE_RESIDUES, βα rows, _no_strays) |",
        f"| boundary moved | {_count(dg, 'boundary moved')} | as many, with a boundary elsewhere | "
        "the lowest normalized cut is not CATH's boundary (the criterion prefers another), a discontinuous CATH "
        "domain one cut cannot take (two runs at most, _lowest_cut), or the element rules moved a right residue cut |",
        f"| extra split | {_count(dg, 'extra split')} | more automatic domains than CATH's | a cut below "
        "DOMAIN_NCUT inside one CATH domain (a lobe or subdomain CATH keeps whole) |",
        "",
        "What the failing chains are like (tags, not exclusive): "
        + ", ".join(f"{k} {v}" for k, v in dg["counts"]["tag"].items())
        + f", none of them {dg['counts']['untagged']}. Discontinuous: a CATH domain in two runs or more along the "
        "chain; tightly packed: a CATH domain whose normalized cut against the rest of the chain (the contacts "
        "residue_domains weighs) is at DOMAIN_NCUT or above; large: three CATH domains or more, or more residues "
        "than DOMAIN_MAX_POINTS. Their CATH boundaries: "
        + ", ".join(f"{k} {v}" for k, v in dg["counts"]["boundaries"].items())
        + " (moved: none within the tolerance, but an automatic boundary no CATH one claims lies between its two "
        "neighbours; extra: automatic boundaries no CATH one claims, left over). "
        f"{dg['counts']['residue_level_meets']} of the failures meet the targets at the residue level, before the "
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
    splits = doc["splits"]
    sets = [s for s in splits if s["set"] != "r100_third"]
    third = [s for s in splits if s["set"] == "r100_third"]
    first = {}  # the third sample's first chain of each entry (copies of one protein said once)
    for s in third:
        first.setdefault(s["pdb"], s["chain"])
    once = [s for s in third if first[s["pdb"]] == s["chain"]]
    lines += [
        "",
        "## 2. The confidence of each split",
        "",
        "A split's confidence is 1 - (its normalized cut / DOMAIN_NCUT): 0 at the limit, where the part would have "
        "stayed whole, 1 for a part with no contact across (features.split_confidence, said by domain_confidence in "
        "the notes, `foldmap check`, its --json report, the MCP result and the app; never drawn). A split is right "
        f"here when each place its two sides meet lies within {doc['tolerance']} residues of a CATH boundary.",
        "",
        "| confidence | domain sets: right | wrong | third sample: right | wrong | its first chain per entry: right "
        "| wrong |",
        "|---|---|---|---|---|---|---|",
    ]
    for (rng, r1, w1), (_, r2, w2), (_, r3, w3) in zip(_bins(sets), _bins(third), _bins(once)):
        lines.append(f"| {rng} | {r1} | {w1} | {r2} | {w2} | {r3} | {w3} |")
    for name, group in (
        ("domain sets", sets),
        ("third sample", third),
        ("third sample, its first chain per entry", once),
    ):
        right = [s["confidence"] for s in group if s["right"]]
        wrong = [s["confidence"] for s in group if not s["right"]]
        if group:
            lines.append(
                f"\n{name}: {len(group)} splits, {len(right)} right; median confidence of the right ones "
                f"{np.median(right) if right else float('nan'):.2f}, of the wrong ones "
                f"{np.median(wrong) if wrong else float('nan'):.2f}."
            )
    lines += [
        "",
        "## 3. Keeping a weakly supported region undivided",
        "",
        "A region whose split's confidence is below t kept whole is a cut made only below DOMAIN_NCUT x (1 - t). "
        "Measured for each t on every set (chains meeting the targets of all chains; mean overlap, F1, recall):",
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
    if not un["decisions"]:
        lines.append("No t above 0 is chosen by either reading: the change is not made.")
    for key, dec in un["decisions"].items():
        worse = [n for n, ok in dec["no_worse"].items() if not ok]
        if not dec["third_entries_no_worse"]:
            worse.append("r100_third's entries failing S5")
        verdict = "no set is worse" if dec["commit"] else f"worse on {', '.join(worse)}"
        lines += ["", f"At t = {float(key):.2f}: {verdict}. Chains whose scores change:", ""]
        for name, got in dec["changed"].items():
            for c in got:
                a, b = c["at_0"], c["at_t"]
                lines.append(
                    f"- {name} {c['pdb']} {c['chain']}: domains {a['n_auto']} -> {b['n_auto']}, overlap "
                    f"{a['overlap']:.2f} -> {b['overlap']:.2f}, F1 {a['f1']:.2f} -> {b['f1']:.2f}, recall "
                    f"{a['recall']:.2f} -> {b['recall']:.2f}, meets {a['meets']} -> {b['meets']}"
                )
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
        "on. The domain sets are tuning (dom_cal) and regression (dom_held, dom_fresh) sets.",
        "- CATH is one chopping; another (SCOP, ECOD) may draw some boundaries elsewhere. A wrong split here is one "
        "CATH does not draw.",
        "- The classes are read from domain counts and boundaries; a chain's mechanism is the first that its cut "
        "record explains, not every one that may contribute.",
        "- The confidence measures how far a cut is from not being made; it says nothing of where else a boundary "
        "could lie.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--jobs", type=int, default=12)
    args = parser.parse_args(argv)
    doc = {"source_stamp": SOURCE_STAMP, **build(args.jobs)}  # the stamp taken as this process started
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "report.json").write_text(json.dumps(doc, indent=1, default=str) + "\n", encoding="utf-8", newline="\n")
    (args.out / "report.md").write_text(render(doc), encoding="utf-8", newline="\n")
    print(f"wrote {args.out / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
