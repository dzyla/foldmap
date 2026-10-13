# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""What an automatic choice left out of a figure, so the app, the command line and the MCP server can say so:
`--assembly auto` drawing assembly 1 when the file holds several copies of a complex (io.read_model), `--focus
auto` drawing one subunit of a cage, a filament or a large symmetric assembly (cli._focus_chains), and `na_extent
auto` cropping a large nucleic acid to the stems the protein touches (cli._na_extent). Chains asked for with
--chains are said as asked for, not as left out."""

from __future__ import annotations

from .io import _polymers

CAGES = ("T", "O", "I")


def _subunits_reason(n: int, kind: str) -> str:
    if kind == "helical":
        return f"one of {n} subunits of a helical filament"
    if kind in CAGES:
        return f"one of {n} subunits of {'an' if kind in 'OI' else 'a'} {kind} cage"
    return f"one of {n} subunits ({kind} symmetry)"


def _na_kind(layout) -> str:
    strands = layout.nucleic.strands if layout.nucleic is not None else []
    kinds = {s.rna for s in strands}
    return "RNA" if kinds == {True} else "DNA" if kinds == {False} else "nucleic acid"


def _undrawn_chains(model, layout) -> list[tuple[str, str]]:
    """(chain, why) for each polymer chain read for the figure (and asked for, and in the focus) that it has nothing
    of: a protein chain with no helix or strand (Layout.coil_chains: a bound peptide, a lasso peptide), or a
    nucleic-acid chain no strand of the figure holds (7MW8 chain K: a single nucleotide)."""
    focus, asked = layout.focus_chains, set(getattr(layout, "chains_asked", None) or [])
    coil = dict(getattr(layout, "coil_chains", None) or [])
    protein = set(getattr(layout, "res_chain", None) or [])
    strands = {s.chain for s in layout.nucleic.strands} if layout.nucleic is not None else set()
    out = []
    for c in getattr(layout, "chains_read", None) or []:
        if c in layout.partial or (focus and c not in focus) or (asked and c not in asked):
            continue
        if c in coil:
            out.append((c, f"{coil[c]} residues, no helix or strand"))
        elif c not in protein and c not in strands and not getattr(layout, "na_cropped", False):
            n = next((len(ch.get_polymer()) for ch in model[0] if ch.name == c), 0)
            out.append((c, f"{n} nucleotide{'' if n == 1 else 's'}, too short to draw"))
    return out


def _hidden_additives(layout) -> list[str]:
    """The additives (io.ADDITIVES: crystallisation and purification leftovers) touching the drawn chains that the
    figure hides (Style.ligands auto), by code: the note names them (random-100 C20 ruling)."""
    links = getattr(layout, "links", None)
    if getattr(layout, "ligands", "auto") != "auto" or links is None:
        return []
    focus = layout.focus_chains
    return sorted(
        {
            g.name
            for g in getattr(links, "ligands", None) or []
            if g.additive and (not focus or any(layout.res_chain[r] in focus for r in g.contacts))
        }
    )


def drawn_scope(model, layout, assembly: str = "auto", focus: str = "auto") -> dict:
    """What the figure drew of the structure, given the io.read_model result it was built from (`model`), its
    layout and the assembly and focus settings it was drawn with. Only automatic choices count (an assembly or a
    focus asked for by name draws what was asked; chains asked for, Layout.chains_asked, are said as asked for):

    file_chains: chains (protein and nucleic acid) of the structure as automatic choices did not cut it (the
        asymmetric unit when assembly 1 is one of its copies; else the assembly read); drawn_chains: those drawn
        in full; reason: what was drawn, '' when nothing was left out and no chain was asked for, made of
        chains_reason (what was drawn of the chains) and crop (what of the nucleic acid, '' when none was
        cropped); all_options: the settings that draw everything ({} when nothing was left out: `set` is a style
        key=value, as --set takes it); residues: polymer residues those settings would draw (of the chains asked
        for, when some were); undrawn_chains: chains the figure has nothing of (_undrawn_chains), not_drawn: them
        with why ('' when none), part of the reason (random-100 C15); additives: the additives it hides on purpose
        (_hidden_additives), named by scope_note but not part of the reason (C20 ruling)."""
    several = getattr(model, "several", None) if assembly == "auto" else None
    read = list(getattr(layout, "chains_read", None) or [])
    asked = list(getattr(layout, "chains_asked", None) or [])
    reasons, options = [], {}
    total, drawn, residues = len(read), len(read), None
    if asked:  # the visitor chose the chains: nothing an automatic choice left out of them
        several = None
    if several:
        reasons.append(f"assembly {several['assembly']}, one of {several['copies']} copies in the file")
        options["assembly"] = "asu"
        if focus == "auto":  # the whole asymmetric unit, never one subunit of it
            options["focus"] = "none"
        total, residues = several["chains"], several["residues"]
    subunits = getattr(layout, "subunits", None)
    if focus == "auto" and layout.focus_chains and subunits and not asked:
        strands = layout.nucleic.strands if layout.nucleic is not None else []
        drawn = len(set(layout.focus_chains) | {s.chain for s in strands})
        reasons.append(_subunits_reason(*subunits))
        options["focus"] = "none"
    if reasons and drawn >= total:  # every chain drawn after all
        reasons, options = [], {}
    crop = ""
    if getattr(layout, "na_cropped", False) and layout.nucleic is not None:
        crop = (
            f"the {_na_kind(layout)} cropped to the stems the protein touches ({layout.nucleic.size:,} of "
            f"{layout.na_read:,} nucleotides)"
        )
        options["set"] = "na_extent=all"
    if asked:
        reasons.insert(0, f"chains {', '.join(asked)} (asked for)")
    undrawn = _undrawn_chains(model, layout)
    additives = _hidden_additives(layout)
    not_drawn = "; ".join(f"chain {c} ({why})" for c, why in undrawn)
    drawn -= len(undrawn)
    unsaid = {"undrawn_chains": [c for c, _ in undrawn], "additives": additives, "not_drawn": not_drawn}
    if not reasons and not crop and not not_drawn:
        return {
            "file_chains": total,
            "drawn_chains": drawn,
            "reason": "",
            "chains_reason": "",
            "crop": "",
            "all_options": {},
            "residues": 0,
            **unsaid,
        }
    if residues is None and options:
        residues = sum(len(c.get_polymer()) for c in model[0] if c.name in read) if asked else _polymers(model[0])[1]
    return {
        "file_chains": total,
        "drawn_chains": drawn,
        "reason": "; ".join([*reasons, *([crop] if crop else []), *([f"not drawn: {not_drawn}"] if not_drawn else [])]),
        "chains_reason": "; ".join(reasons),
        "crop": crop,
        "all_options": options,
        "residues": residues or 0,
        **unsaid,
    }


def scope_note(scope: dict, how: str = "cli") -> str:
    """One line saying what was drawn of the structure and how to draw everything ('' when nothing was left out
    and no chain asked for), then what it has not drawn (chains with nothing to draw, additives hidden): with
    command-line flags (how='cli') or MCP tool arguments (how='mcp'); the chain count when chains were left out."""
    if not scope.get("reason") and not scope.get("additives"):
        return ""
    options = scope["all_options"]
    if how == "mcp":
        draw = ", ".join(
            f'style={{"{v.split("=", 1)[0]}": "{v.split("=", 1)[1]}"}}' if k == "set" else f'{k}="{v}"'
            for k, v in options.items()
        )
    else:
        draw = " ".join(f"--{k} {v}" for k, v in options.items())
    what = scope_what(scope)
    parts = [f"drew {what}" + (f"; {draw} draws all of them" if options else "")] if what else []
    unsaid = [scope["not_drawn"]] if scope.get("not_drawn") else []
    if scope.get("additives"):  # hidden on purpose, but named (C20 ruling)
        shows = 'style={"ligands": "all"}' if how == "mcp" else "--set ligands=all"
        unsaid.append(f"additives {', '.join(scope['additives'])} ({shows} draws them)")
    if unsaid:
        parts.append(f"not drawn: {'; '.join(unsaid)}")
    return "; ".join(parts)


def scope_what(scope: dict, count: str = ": {drawn} of {total} chains") -> str:
    """What was drawn, its clauses apart: what of the chains with their count (`count`, when chains were left out),
    then what of the nucleic acid, so the chain count never reads as part of the crop."""
    head = scope.get("chains_reason", scope["reason"])
    if head and scope["drawn_chains"] < scope["file_chains"]:
        head += count.format(drawn=scope["drawn_chains"], total=scope["file_chains"])
    return "; ".join(part for part in (head, scope.get("crop", "")) if part)
