# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The confidence of each inferred domain split (plan 28 Task 5, the review's 3.4, rulings Q42, Q43): the automatic
domains (--domains auto) are said with how strongly the residue contacts support each split; the figure is unchanged.

A split is a cut the normalized cut made (features.residue_domains: the part's lowest cut, made when below
DOMAIN_NCUT) that the drawn domains keep (features.domain_splits: no domain on both of its sides; a cut the later
rules undid, a domain too small to keep, a βα row kept whole or a membrane band joined, is none). Its confidence
(features.split_confidence) is

    1 - (its normalized cut / DOMAIN_NCUT)

0 at the limit, where the part would have stayed whole, 1 for a part with no contact across. It is said with what the
cut weighed: the CA contacts across it, the density of contacts across over that within its two sides, and the sides'
sizes and residues. It is the margin below a tuned limit, nothing more: not a probability that CATH or any other
chopping draws the same boundary, and it moves with the limit (DOMAIN_NCUT is sensitive: docs/validation/sensitivity,
x1.25 moves the third random sample's S5 failures from 12 to 19). Every interface says so (definition()).

What it is worth against CATH is measured in docs/validation/domain-uncertainty (scripts/domain_uncertainty.py):
informative, not decisive. On the 44 splits of the two regression domain sets (dom_held, dom_fresh), 3 of 11 under
0.25 fall on a CATH boundary and 5 of 5 at 0.75 or more do (AUC 0.75); on dom_cal, where DOMAIN_NCUT was chosen, the AUC
is 0.61; n is small (Wilson intervals in the record). No threshold is drawn from it: keeping splits under a confidence
undivided lowers some calibration chain at every level tried (ruling Q42; the record).

Where it is said (said(), guarded: a failure to say it is said, never raised, so a figure is always written): the
command line's notes (`foldmap plot`), `foldmap check` (lines under S5, and S5's measures in the --json report: an
optional key of the v1 shape, S5's status unchanged), the MCP result (its note and `domain_splits`) and the app (a
line under the figure). Never in the drawn figure or its files: the figure is the same with or without it (bench
--exact)."""

from __future__ import annotations

NOT_JUDGED = "domain split confidence not judged"


def definition() -> str:
    """What the confidence is, as every interface says it: the limit's value read from features.DOMAIN_NCUT (so it
    follows the constant), and that it is no probability and moves with the limit."""
    from .features import DOMAIN_NCUT

    return (
        f"1 - normalized cut / {DOMAIN_NCUT:g} (DOMAIN_NCUT, the tuned split limit): the margin below that limit, "
        "not a probability; it moves with the limit"
    )


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


def _failed(splits: list[dict]) -> str:
    """The records that failed, in words ('' when none): per chain when the chain is known."""
    errors = [s for s in splits if "error" in s]
    if not errors:
        return ""
    said = ", ".join(f"chain {s['chain']}: {s['error']}" if s.get("chain") else s["error"] for s in errors)
    return f"{NOT_JUDGED} for {said}" if any(s.get("chain") for s in errors) else f"{NOT_JUDGED}: {said}"


def note(splits: list[dict]) -> str:
    """One line for the notes: each inferred split with its confidence, then any record that failed (a failure never
    hides the splits that were said); '' when there is none."""
    good = [s for s in splits if "error" not in s]
    parts = []
    if good:
        said = "; ".join(describe(s) for s in good)
        parts.append(f"inferred domain splits, with their confidence ({definition()}; said, not judged): {said}")
    if failed := _failed(splits):
        parts.append(failed)
    return "; ".join(parts)


def said(drawn) -> tuple[list[dict], str]:
    """(the inferred splits, the notes' line: note()) for the figure `drawn`; guarded: a failure is said ("domain
    split confidence not judged: ...", and the splits are [{"error": its type}]), never raised."""
    try:
        splits = splits_of(drawn)
        return splits, note(splits)
    except Exception as err:  # noqa: BLE001 (said in the notes: the confidence never stops a figure)
        return [{"error": type(err).__name__}], f"{NOT_JUDGED}: {type(err).__name__}"


S5_REASON = "not measured here (no reference chopping); measures: each inferred split's confidence, said, never judged"
S5_LINE = "not yet measured (no reference chopping); inferred splits said below, not judged"


def measures(splits: list[dict]) -> dict:
    """S5's measures in `foldmap check --json` (deep.report): the splits (a failed record as {"error"}) and what their
    confidence is."""
    return {"domain_splits": splits, "confidence": definition()}


def lines(splits: list[dict]) -> list[str]:
    """`foldmap check`'s lines under S5: each inferred split with its confidence, said, not judged, then any record
    that failed; [] without."""
    good = [s for s in splits if "error" not in s]
    out = []
    if good:
        out.append(f"  inferred domain splits, with their confidence ({definition()}; said, not judged):")
        out += [f"    {describe(s)}" for s in good]
    if failed := _failed(splits):
        out.append(f"  {failed}")
    return out


def under_s5(science: list[str], splits: list[dict]) -> list[str]:
    """`foldmap check`'s science lines with the split lines (lines()) under S5, after S5's own findings; S5's line says
    the splits follow (S5_LINE) when they do and S5 was not measured. Unchanged without splits."""
    said = lines(splits)
    if not said:
        return list(science)
    out = list(science)
    at = next((k for k, line in enumerate(out) if line.startswith("S5 ")), None)
    if at is None:
        return out + said
    if out[at].endswith(": not yet measured"):
        out[at] = out[at].removesuffix("not yet measured") + S5_LINE
    end = at + 1
    while end < len(out) and out[end].startswith("  "):
        end += 1
    return out[:end] + said + out[end:]
