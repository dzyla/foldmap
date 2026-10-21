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


def _nucleic_kinds(model) -> dict[str, str]:
    """Chain -> DNA, RNA or nucleic acid, for each nucleic-acid chain of the model read."""
    import gemmi

    kinds = {gemmi.PolymerType.Dna: "DNA", gemmi.PolymerType.Rna: "RNA", gemmi.PolymerType.DnaRnaHybrid: "nucleic acid"}
    out = {}
    for c in model[0]:
        polymer = c.get_polymer()
        kind = kinds.get(polymer.check_polymer_type()) if len(polymer) else None
        if kind:
            out[c.name] = kind
    return out


def _na_outside(model, layout) -> list[str]:
    """The nucleic-acid chains read that a subunit focus leaves out: in no focus chain and no strand drawn (8QEM: a
    C12 ring's 63-bp duplex, held by no protomer; F1 review, Important 3)."""
    strands = {s.chain for s in layout.nucleic.strands} if layout.nucleic is not None else set()
    focus = set(layout.focus_chains or [])
    kinds = _nucleic_kinds(model)
    return [c for c in getattr(layout, "chains_read", None) or [] if c in kinds and c not in strands | focus]


def na_outside_text(model, chains: list[str]) -> str:
    """'DNA Y, Z' (or RNA, or nucleic acid when mixed) for the nucleic-acid chains a focus leaves out."""
    kinds = {_nucleic_kinds(model).get(c, "nucleic acid") for c in chains}
    kind = kinds.pop() if len(kinds) == 1 else "nucleic acid"
    shown = chains[:6] + ([f"+{len(chains) - 6} more"] if len(chains) > 6 else [])
    return f"{kind} {', '.join(shown)}"


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
    asked = set(getattr(layout, "chains_asked", None) or [])
    asked |= {getattr(layout, "like_chain", "")} - {""}  # --like's chosen chain: the others as if not asked for
    focus = layout.focus_chains
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


def _na_ions(layout) -> dict[str, int]:
    """Metal ions held by nucleotides alone (Links.na_ions), by name: never marked (a ribosome holds hundreds), so
    the note names them (random-100 C33)."""
    links = getattr(layout, "links", None)
    return dict(sorted((getattr(links, "na_ions", None) or {}).items())) if layout.nucleic is not None else {}


def _tetrads(layout) -> int:
    """G-tetrads of the drawn nucleic acid: foldmap draws none, their guanines unpaired (C29 ruling: said)."""
    return len(layout.nucleic.tetrads) if layout.nucleic is not None else 0


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
        (_hidden_additives), named by scope_note but not part of the reason (C20 ruling); tetrads: the G-tetrads of
        the drawn nucleic acid, drawn unpaired, named by scope_note but not part of the reason (C29 ruling); na_ions:
        metal ions held by nucleotides alone, by name, never marked, named by scope_note likewise (C33);
        helix_matched: how many helix runs were cut or joined so copies of a protein match (ss._copy_breaks: a fact,
        not part of the reason or the note; fix round 1 ruling)."""
    several = getattr(model, "several", None) if assembly == "auto" else None
    left_out = list(getattr(model, "left_out", None) or []) if assembly == "auto" else []
    read = list(getattr(layout, "chains_read", None) or [])
    asked = list(getattr(layout, "chains_asked", None) or [])
    chosen = getattr(layout, "like_chain", "") or ""  # one chain of several, chosen to draw like a reference
    reasons, options, na_outside = [], {}, []
    total, drawn, residues = len(read), len(read), None
    if asked or chosen:  # the visitor chose the chains (or --like one): nothing an automatic choice left out of them
        several, left_out = None, []
    if left_out:  # assembly 1 holds no copy of some chains of the file (1RM5's chain O: final review B, Important 2)
        name = getattr(model, "assembly", None) or "1"
        held = f"chain{'s' if len(left_out) > 1 else ''} {', '.join(left_out)}"
        reasons.append(f"assembly {name}, which holds no copy of {held} of the file")
        options["assembly"] = "asu"
        if focus == "auto":
            options["focus"] = "none"
        # counted in the file's chains, as --assembly asu draws them (F1 review, Minor 1: not copies and file chains)
        total = getattr(model, "asu_chains", None) or total + len(left_out)
        drawn = total - len(left_out)
        residues = getattr(model, "asu_residues", None)
    if several:
        reasons.append(f"assembly {several['assembly']}, one of {several['copies']} copies in the file")
        options["assembly"] = "asu"
        if focus == "auto":  # the whole asymmetric unit, never one subunit of it
            options["focus"] = "none"
        total, residues = several["chains"], several["residues"]
    subunits = getattr(layout, "subunits", None)
    if focus == "auto" and layout.focus_chains and subunits and not asked and not chosen:
        strands = layout.nucleic.strands if layout.nucleic is not None else []
        drawn = len(set(layout.focus_chains) | {s.chain for s in strands})
        na_outside = _na_outside(model, layout)
        reasons.append(
            _subunits_reason(*subunits) + (f", without the {na_outside_text(model, na_outside)}" if na_outside else "")
        )
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
    elif chosen:  # never said as asked for (13a-7): --like drew its first protein chain
        reasons.insert(0, f"chain {chosen} (the first protein chain: one chain of each is drawn like a reference)")
    undrawn = _undrawn_chains(model, layout)
    additives = _hidden_additives(layout)
    not_drawn = "; ".join(f"chain {c} ({why})" for c, why in undrawn)
    drawn -= len(undrawn)
    unsaid = {
        "na_outside": na_outside,
        "undrawn_chains": [c for c, _ in undrawn],
        "additives": additives,
        "not_drawn": not_drawn,
        "tetrads": _tetrads(layout),
        "na_ions": _na_ions(layout),
        "helix_matched": {"cut": 0, "joined": 0, **(getattr(layout, "helix_matched", None) or {})},
    }
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
        residues = (
            sum(len(c.get_polymer()) for c in model[0] if c.name in read) if asked or chosen else _polymers(model[0])[1]
        )
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
    and no chain asked for), then what it has not drawn (chains with nothing to draw, additives hidden, metal ions
    held by the nucleic acid alone): with
    command-line flags (how='cli') or MCP tool arguments (how='mcp'); the chain count when chains were left out;
    G-tetrads, which no setting draws."""
    if not scope.get("reason") and not scope.get("additives") and not scope.get("tetrads") and not scope.get("na_ions"):
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
    unsaid = not_drawn(scope, how)
    if unsaid:
        parts.append(f"not drawn: {'; '.join(unsaid)}")
    return "; ".join(parts)


def not_drawn(scope: dict, how: str = "cli") -> list[str]:
    """What the figure has not drawn, each as the notes say it: chains with nothing to draw, additives hidden (with
    the command-line flag, how='cli', or MCP tool argument, how='mcp', that draws them; none named for how='app',
    whose Ligands control says it hides them), metal ions held by the nucleic acid alone, G-tetrads. The command
    line, MCP and the app say these in the same words (final review C-I2)."""
    unsaid = [scope["not_drawn"]] if scope.get("not_drawn") else []
    if scope.get("additives") and how != "app":  # hidden on purpose, but named (C20 ruling)
        shows = 'style={"ligands": "all"}' if how == "mcp" else "--set ligands=all"
        unsaid.append(f"additives {', '.join(scope['additives'])} ({shows} draws them)")
    if scope.get("na_ions"):  # read, never marked (C33)
        unsaid.append(na_ions_text(scope["na_ions"]))
    if scope.get("tetrads"):  # read, never drawn (C29 ruling)
        n = scope["tetrads"]
        unsaid.append(f"{n} G-tetrad{'' if n == 1 else 's'} of a G-quadruplex (their Gs drawn unpaired)")
    return unsaid


def na_ions_text(ions: dict[str, int]) -> str:
    """The metal ions held by the nucleic acid alone (drawn_scope na_ions), as the notes say them."""
    return (
        "metal ions held by the nucleic acid alone ("
        + ", ".join(f"{name} ×{n}" if n > 1 else name for name, n in ions.items())
        + ")"
    )


def scope_what(scope: dict, count: str = ": {drawn} of {total} chains") -> str:
    """What was drawn, its clauses apart: what of the chains with their count (`count`, when chains were left out),
    then what of the nucleic acid, so the chain count never reads as part of the crop."""
    head = scope.get("chains_reason", scope["reason"])
    if head and scope["drawn_chains"] < scope["file_chains"]:
        head += count.format(drawn=scope["drawn_chains"], total=scope["file_chains"])
    return "; ".join(part for part in (head, scope.get("crop", "")) if part)
