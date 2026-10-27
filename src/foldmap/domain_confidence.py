# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The confidence of each inferred domain split (plan 28 Task 5, the review's 3.4, ruling Q42): the automatic domains
(--domains auto) are said with how strongly the residue contacts support each split; the figure is unchanged.

A split is a cut the normalized cut made (features.residue_domains: the part's lowest cut, made when below
DOMAIN_NCUT) that the drawn domains keep (features.domain_splits: no domain on both of its sides; a cut the later
rules undid, a domain too small to keep, a βα row kept whole or a membrane band joined, is none). Its confidence
(features.split_confidence) is

    1 - (its normalized cut / DOMAIN_NCUT)

0 at the limit, where the part would have stayed whole, 1 for a part with no contact across. It is said with what the
cut weighed: the CA contacts across it, the density of contacts across over that within its two sides, and the sides'
sizes and residues. It says how far a split is from not being made, nothing more: it is not a probability that CATH
or any other chopping draws the same boundary. What it is worth against CATH is measured, not assumed, in
docs/validation/domain-uncertainty (scripts/domain_uncertainty.py): on the domain sets and the third random sample,
most splits under 0.25 are boundaries CATH does not draw and most at 0.75 or more are ones it draws, but not all: a
split of high confidence can still be a subdomain CATH keeps whole. No threshold is drawn from it: keeping a split
under some confidence undivided made the tuning set worse at every level tried (ruling Q42; the record).

Where it is said (said(), guarded: a failure to say it is said, never raised, so a figure is always written): the
command line's notes (`foldmap plot`), `foldmap check` (a line under S5, and S5's measures in the --json report: an
added field, the v1 shape unchanged), the MCP result (its note and `domain_splits`) and the app (a line under the
figure). Never in the drawn figure or its files: the figure is the same with or without it (bench --exact)."""

from __future__ import annotations

CONFIDENCE = "1 - (the split's normalized cut / DOMAIN_NCUT): 0 at the limit, 1 with no contact across"


def splits_of(drawn) -> list[dict]:
    """The inferred splits of the figure `drawn` (spec.Drawn), each with its confidence (Layout.domain_splits,
    features.domain_splits); [] when its domains are not inferred or none is split."""
    return list(getattr(drawn.layout, "domain_splits", None) or [])


def _side(names: list[str]) -> str:
    return ", ".join(names)


def describe(split: dict) -> str:
    """One split in words: its domains on either side, its confidence and what it weighed."""
    ratio = split.get("density_ratio")
    dense = f", {ratio:.2f} of the contact density within" if ratio is not None else ""
    a, b = split["residues"]
    return (
        f"{_side(split['between'][0])} | {_side(split['between'][1])} {split['confidence']:.2f} "
        f"(normalized cut {split['ncut']:.2f}; {split['across']} contacts across{dense}; {a} and {b} residues)"
    )


def note(splits: list[dict]) -> str:
    """One line for the notes: each inferred split with its confidence; '' when there is none."""
    if not splits:
        return ""
    if errors := sorted({s["error"] for s in splits if "error" in s}):
        return f"domain split confidence not judged: {', '.join(errors)}"
    said = "; ".join(describe(s) for s in splits)
    return f"inferred domain splits, with their confidence ({CONFIDENCE}; said, not judged): {said}"


def said(drawn) -> tuple[list[dict], str]:
    """(the inferred splits, the notes' line: note()) for the figure `drawn`; guarded: a failure is said ("domain
    split confidence not judged: ..."), never raised."""
    try:
        splits = splits_of(drawn)
        return splits, note(splits)
    except Exception as err:  # noqa: BLE001 (said in the notes: the confidence never stops a figure)
        return [], f"domain split confidence not judged: {type(err).__name__}"


def measures(splits: list[dict]) -> dict:
    """S5's measures in `foldmap check --json` (deep.report): the splits and what their confidence is."""
    return {"domain_splits": splits, "confidence": CONFIDENCE}


def lines(splits: list[dict]) -> list[str]:
    """`foldmap check`'s lines under S5: each inferred split with its confidence, said, not judged; [] without."""
    if not splits:
        return []
    if errors := sorted({s["error"] for s in splits if "error" in s}):
        return [f"  domain split confidence not judged: {', '.join(errors)}"]
    head = f"  inferred domain splits, with their confidence ({CONFIDENCE}; said, not judged):"
    return [head, *(f"    {describe(s)}" for s in splits)]
