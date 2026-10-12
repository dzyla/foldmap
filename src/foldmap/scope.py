# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""What an automatic choice left out of a figure, so the app, the command line and the MCP server can say so:
`--assembly auto` drawing assembly 1 when the file holds several copies of a complex (io.read_model), and
`--focus auto` drawing one subunit of a cage, a filament or a large symmetric assembly (cli._focus_chains)."""

from __future__ import annotations

from .io import _polymers

CAGES = ("T", "O", "I")


def _subunits_reason(n: int, kind: str) -> str:
    if kind == "helical":
        return f"one of {n} subunits of a helical filament"
    if kind in CAGES:
        return f"one of {n} subunits of {'an' if kind in 'OI' else 'a'} {kind} cage"
    return f"one of {n} subunits ({kind} symmetry)"


def drawn_scope(model, layout, assembly: str = "auto", focus: str = "auto") -> dict:
    """What the figure drew of the structure, given the io.read_model result it was built from (`model`), its
    layout and the assembly and focus settings it was drawn with. Only automatic choices count (an assembly or a
    focus asked for by name draws what was asked):

    file_chains: chains (protein and nucleic acid) of the structure as automatic choices did not cut it (the
        asymmetric unit when assembly 1 is one of its copies; else the assembly read); drawn_chains: those drawn
        in full; reason: what was drawn, '' when nothing was left out; all_options: the settings that draw every
        chain ({} when nothing was left out); residues: polymer residues those settings would draw."""
    several = getattr(model, "several", None) if assembly == "auto" else None
    read = list(getattr(layout, "chains_read", None) or [])
    reasons, options = [], {}
    total, drawn, residues = len(read), len(read), None
    if several:
        reasons.append(f"assembly {several['assembly']}, one of {several['copies']} copies in the file")
        options["assembly"] = "asu"
        if focus == "auto":  # the whole asymmetric unit, never one subunit of it
            options["focus"] = "none"
        total, residues = several["chains"], several["residues"]
    subunits = getattr(layout, "subunits", None)
    if focus == "auto" and layout.focus_chains and subunits:
        strands = layout.nucleic.strands if layout.nucleic is not None else []
        drawn = len(set(layout.focus_chains) | {s.chain for s in strands})
        reasons.append(_subunits_reason(*subunits))
        options["focus"] = "none"
    if not reasons or drawn >= total:
        return {"file_chains": total, "drawn_chains": drawn, "reason": "", "all_options": {}, "residues": 0}
    if residues is None:
        residues = _polymers(model[0])[1]
    return {
        "file_chains": total,
        "drawn_chains": drawn,
        "reason": "; ".join(reasons),
        "all_options": options,
        "residues": residues,
    }


def scope_note(scope: dict, how: str = "cli") -> str:
    """One line saying what was left out and how to draw everything ('' when nothing was): with command-line flags
    (how='cli') or MCP tool arguments (how='mcp')."""
    if not scope.get("reason"):
        return ""
    options = scope["all_options"]
    if how == "mcp":
        draw = ", ".join(f'{k}="{v}"' for k, v in options.items())
    else:
        draw = " ".join(f"--{k} {v}" for k, v in options.items())
    return f"drew {scope['reason']}: {scope['drawn_chains']} of {scope['file_chains']} chains; {draw} draws all of them"
