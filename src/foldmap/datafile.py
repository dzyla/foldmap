# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""`foldmap data` (plan 13d): what a figure says about a structure, as data. A laid-out structure (cli.make_layout:
no routing, no drawing) becomes one versioned dict, its names the figure's names:

- structure: the file's name, the chains read and drawn, focus;
- elements: each helix, 3-10 helix and strand, its residues (insertion codes kept), the label the figure gives it
  and where the figure puts it (page units, y up);
- sheets: each sheet's strands in the order the figure draws them, their directions, each pairing's kind and
  residue register (the DSSP bridges), accessory strands, and for a closed sheet the barrel's strand number n and
  shear number S (the bridges followed once round the barrel);
- disulfides; symmetry and stoichiometry; with --membrane, the elements that span the drawn slab (A18's rule) or
  cross its mid-plane, and the sides of elements and termini.

Every list has a fixed order, every number is rounded to 0.001 (±0 written 0.0, angles in (-180, 180]), so the
same structure and options give the same bytes (A12). `schema()` is the format, as a JSON Schema; SCHEMA_VERSION
is its version, the dict's "foldmap_data" key."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from statistics import median_low

import numpy as np

SCHEMA_VERSION = 1
SECTIONS = ("structure", "elements", "sheets", "disulfides", "symmetry", "membrane")
ALL = SECTIONS
KINDS = {"H": "helix", "E": "strand", "G": "helix_3_10"}


def _num(x) -> float | None:
    """A number as written: rounded to 0.001, -0.0 as 0.0; None for nan."""
    v = float(x)
    if v != v:
        return None
    return round(v, 3) + 0.0


def _angle(radians: float) -> float:
    """An angle in degrees, rounded, in (-180, 180] (-180 and 180 are one direction: written 180)."""
    v = _num(np.degrees(radians))
    while v <= -180.0:
        v = _num(v + 360.0)
    while v > 180.0:
        v = _num(v - 360.0)
    return v


def _vec(v) -> list[float]:
    return [_num(x) for x in np.asarray(v, float)]


def _res(label) -> dict:
    """A residue (model.ResLabel): its chain, number and insertion code."""
    return {"chain": label.chain, "seq": int(label.seq), "icode": label.icode}


def structure_data(layout, sses, bb, model=None, sections: tuple[str, ...] = ALL, *, name: str | None = None) -> dict:
    """The data of a laid-out structure (cli.make_layout's layout, sses and backbone; `model`: io.read_model's, not
    needed by these sections), the sections asked for (SECTIONS; ValueError naming them for another). `name`: the
    structure as the user named it (written as its file name)."""
    unknown = [s for s in sections if s not in SECTIONS]
    if unknown:
        raise ValueError(f"no data section {', '.join(unknown)}; the sections are {', '.join(SECTIONS)}")
    from . import __version__

    out: dict = {"foldmap_data": SCHEMA_VERSION, "foldmap": __version__}
    sheets = _sheets(layout, sses, bb) if {"elements", "sheets"} & set(sections) else []
    for section in SECTIONS:
        if section not in sections:
            continue
        if section == "structure":
            out[section] = _structure(layout, sses, bb, name)
        elif section == "elements":
            out[section] = _elements(layout, sses, sheets)
        elif section == "sheets":
            out[section] = sheets
        elif section == "disulfides":
            out[section] = _disulfides(layout, sses, bb)
        elif section == "symmetry":
            out[section] = _symmetry(layout, bb)
        else:
            out[section] = _membrane(layout, sses, bb)
    return out


def dumps(data: dict) -> str:
    """The dict as written: JSON, keys sorted, one-space indent, LF line ends, a final newline."""
    return json.dumps(data, indent=1, sort_keys=True, ensure_ascii=False) + "\n"


def write(data: dict, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(dumps(data))
    return path


# ---------------------------------------------------------------------------------------------- sections


def _structure(layout, sses, bb, name) -> dict:
    protein = list(dict.fromkeys(l.chain for l in bb.labels))
    na = layout.nucleic.strands if layout.nucleic is not None else []
    return {
        "name": Path(str(name)).name if name else None,
        "residues": len(bb),
        "protein_chains": protein,
        "nucleic_acid_chains": list(dict.fromkeys(st.chain for st in na)),
        "chains_read": list(layout.chains_read),
        "focus_chains": sorted(layout.focus_chains) if layout.focus_chains else None,
        "lent_chains": sorted(layout.partial),
        "context": layout.context,
    }


def _elements(layout, sses, sheets: list[dict]) -> list[dict]:
    domain = {i: name for name, ids in layout.domains for i in ids}
    sheet_of = {}
    for sh in sheets:
        for i in sh["strands"]:
            sheet_of[i] = sh["id"]
        for a in sh["accessory"]:
            sheet_of[a["strand"]] = sh["id"]
    out = []
    for s in sses:
        p = layout.placed.get(s.id)
        page = None
        if p is not None:
            page = {
                "x": _num(p.cx),
                "y": _num(p.cy),
                "angle": _angle(p.angle),
                "length": _num(p.length),
                "width": _num(p.width),
            }
        out.append(
            {
                "id": s.id,
                "kind": KINDS.get(s.kind, s.kind),
                "chain": s.chain,
                "first": _res(s.first),
                "last": _res(s.last),
                "residues": len(s),
                "label": p.label if p is not None else None,
                "page": page,
                "sheet": sheet_of.get(s.id),
                "domain": domain.get(s.id),
                "lent": s.chain in layout.partial,
            }
        )
    return out


def _bridges(layout, bb) -> list:
    """The DSSP bridges the figure's sheets are built from (cli.make_layout): a rosette's within each copy."""
    from .dssp import assign_dssp

    bridges = assign_dssp(bb).bridges
    if getattr(layout, "rosette", 0):
        sector, labels = layout.frame.sector, bb.labels
        bridges = [b for b in bridges if sector.get(labels[b.i].chain, -1) == sector.get(labels[b.j].chain, -2)]
    return bridges


def _drawn_order(layout, ids: list[str]) -> list[str] | None:
    """These strands in the order the figure draws them side by side, left to right (bottom to top when they lie
    across the page): across their mean direction. None when one is not drawn or they stand in two sheet blocks."""
    from .figcheck import local_layout

    lay = local_layout(layout)
    if any(i not in lay.placed for i in ids):
        return None
    blocks = {k for k, block in enumerate(lay.sheets) for i in block if i in ids}
    if len(blocks) != 1:
        return None
    at = [lay.placed[i] for i in ids]
    ref = at[0].direction
    mean = sum(p.direction * (1.0 if p.direction @ ref >= 0 else -1.0) for p in at)
    across = np.array([mean[1], -mean[0]])
    if round(float(across[0]), 9) < 0 or (round(float(across[0]), 9) == 0 and across[1] < 0):
        across = -across
    return [i for _, i in sorted((round(float(np.array([p.cx, p.cy]) @ across), 6), i) for i, p in zip(ids, at))]


def _same_cycle(got: list[str], want: list[str], closed: bool) -> bool:
    turns = [want[k:] + want[:k] for k in range(len(want))] if closed else [want]
    return any(got in (t, t[::-1]) for t in turns)


def _pair(a, b, bridges, labels) -> dict:
    """Strands a and b as a pairing: its kind (the bridges' majority, P or A, as sheets._majority) and register
    (each bridge, a's residue then b's, in a's order)."""
    reg = []
    for i, j, kind in bridges:
        if a.start <= i <= a.end and b.start <= j <= b.end:
            reg.append((i, j, kind))
        elif a.start <= j <= a.end and b.start <= i <= b.end:
            reg.append((j, i, kind))
    reg.sort()
    kinds = Counter(k for _, _, k in reg)
    return {
        "a": a.id,
        "b": b.id,
        "kind": ("P" if kinds["P"] > kinds["A"] else "A") if reg else None,
        "bridges": len(reg),
        "register": [[_res(labels[i]), _res(labels[j])] for i, j, _ in reg],
        "_reg": reg,
    }


def _shear(pairs: list[dict]) -> int | None:
    """The shear number S of a barrel (its pairings in order round it, the last joining the last strand to the
    first): a residue followed across each pairing to its partner (x -> x + c parallel, x -> c - x antiparallel, c
    the median over the pairing's bridges) comes back to its own strand S residues along. None when a pairing has
    no bridge."""
    sign, shift = 1, 0
    for pair in pairs:
        reg = pair["_reg"]
        if not reg:
            return None
        s = 1 if pair["kind"] == "P" else -1
        c = median_low([j - s * i for i, j, _ in reg])
        sign, shift = s * sign, s * shift + c
    return abs(shift) if sign == 1 else None


def _sheets(layout, sses, bb) -> list[dict]:
    from .sheets import build_sheets

    bridges = _bridges(layout, bb)
    out = []
    for n, sheet in enumerate(build_sheets(sses, bridges), 1):
        track = [s.id for s in sheet.strands]
        by_id = {s.id: s for s in sheet.strands}
        direction = dict(zip(track, sheet.directions))
        drawn = _drawn_order(layout, track)
        as_drawn = drawn is not None and _same_cycle(drawn, track, sheet.closed)
        order = drawn if as_drawn else track
        sense = direction[order[0]]
        steps = list(zip(order, order[1:])) + ([(order[-1], order[0])] if sheet.closed and len(order) > 2 else [])
        pairs = [_pair(by_id[a], by_id[b], bridges, bb.labels) for a, b in steps]
        barrel = {"n": len(order), "shear": _shear(pairs)} if sheet.closed else None
        for p in pairs:
            del p["_reg"]
        out.append(
            {
                "id": f"S{n}",
                "strands": order,
                "directions": [direction[i] * sense for i in order],
                "as_drawn": as_drawn,
                "pairs": pairs,
                "closed": sheet.closed,
                "barrel": barrel,
                "accessory": [
                    {
                        "strand": a.strand.id,
                        "parent": a.parent.id,
                        "direction": a.direction * sense,
                        "pairs_with": [o.id for o in a.others],
                    }
                    for a in sheet.accessory
                ],
                "chords": [[a.id, b.id] for a, b in sheet.chords],
                "joined": [[a.id, b.id] for a, b in sheet.joined],
            }
        )
    return out


def _element_at(sses):
    def at(r: int) -> str | None:
        return next((s.id for s in sses if s.start <= r <= s.end), None)

    return at


def _disulfides(layout, sses, bb) -> list[dict]:
    links = layout.links
    if links is None:
        return []
    at = _element_at(sses)
    return [
        {"a": _res(bb.labels[i]), "b": _res(bb.labels[j]), "elements": [at(i), at(j)]}
        for i, j in sorted(tuple(sorted(p)) for p in links.disulfides)
    ]


def _letters(k: int) -> str:
    out = ""
    k += 1
    while k:
        k, r = divmod(k - 1, 26)
        out = chr(65 + r) + out
    return out


def _stoichiometry(bb) -> dict:
    """The protein chains read, grouped as copies of one protein (model.homolog_chains): each group's chains, and
    the formula, groups lettered from the most copies (A2B2; a lone chain A)."""
    from .model import homolog_chains

    groups: dict[str, list[str]] = {}
    for chain, (first, _) in homolog_chains(bb).items():
        groups.setdefault(first, []).append(chain)
    ranked = sorted(groups.values(), key=lambda g: -len(g))  # stable: ties in order of first appearance
    entities = [{"letter": _letters(k), "chains": g, "copies": len(g)} for k, g in enumerate(ranked)]
    formula = "".join(e["letter"] + (str(e["copies"]) if e["copies"] > 1 else "") for e in entities)
    return {"formula": formula, "entities": entities}


def _symmetry(layout, bb) -> dict:
    sym = getattr(layout, "symmetry", None)
    found = None
    if sym is not None:
        found = {
            "label": sym.label,
            "kind": sym.kind,
            "n": int(sym.n),
            "protomers": [list(p) for p in sym.protomers],
            "unpaired_chains": list(sym.others),
            "rmsd": _num(sym.rmsd),
            "axis": _vec(sym.axis),
            "centre": _vec(sym.centre),
            "twist": _num(np.degrees(sym.twist)) if sym.kind == "H" else None,
            "rise": _num(sym.rise) if sym.kind == "H" else None,
            "left_handed": bool(sym.left_handed) if sym.kind == "H" else None,
        }
    if layout.rosette:
        drawn = "rosette"
    elif any(m for _, m in layout.copy_ties.values()):
        drawn = "mirror"
    elif layout.focus_chains:
        drawn = "one subunit"
    elif sym is not None:
        drawn = "side by side"
    else:
        drawn = None
    return {"found": found, "drawn": drawn, "stoichiometry": _stoichiometry(bb)}


def _side(depth: float, half: float) -> str:
    return "out" if depth >= half else "in" if depth <= -half else "slab"


def _membrane(layout, sses, bb) -> dict | None:
    """The membrane the figure draws (--membrane): None without one. Two measures of an element crossing it, neither
    stated as more than it measures: `spans_slab`, its CA run within membrane.FACE_REACH (2 Å) of both faces of the
    drawn slab (the rule A18 counts a transmembrane helix by: figcheck.membrane_view); `crosses_midplane`, its CA run
    crosses the slab's centre plane and spans at least half the slab's thickness (the usual reading of a
    transmembrane segment, lenient to a slab estimated too thick or tilted). Its side: "across" when it crosses the
    mid-plane, else out, in or slab by its mean CA depth; `from` and `to`, the sides of its first and last CA. Sides
    are named by the membrane's normal, which points out (the non-cytoplasmic side; `sided` says whether the
    positive-inside rule settled which side that is)."""
    mem = getattr(layout, "bilayer", None)
    if mem is None:
        return None
    from .membrane import FACE_REACH

    depth = mem.depth(bb.ca)

    def sign(r: int) -> str:
        return "out" if depth[r] >= 0 else "in"

    elements, spans, crosses = [], [], []
    for s in sses:
        seg = depth[s.start : s.end + 1]
        lo, hi = float(seg.min()), float(seg.max())
        span = len(seg) >= 2 and lo <= -mem.half + FACE_REACH and hi >= mem.half - FACE_REACH  # A18's rule
        cross = len(seg) >= 2 and lo < 0 < hi and hi - lo >= mem.half  # mid-plane crossed, half the thickness spanned
        if span:
            spans.append(s.id)
        if cross:
            crosses.append(s.id)
        elements.append(
            {
                "element": s.id,
                "side": "across" if cross else _side(float(np.mean(seg)), mem.half),
                "from": sign(s.start),
                "to": sign(s.end),
                "depth": [_num(lo), _num(hi)],
                "spans_slab": bool(span),
                "crosses_midplane": bool(cross),
            }
        )
    termini = []
    for chain in dict.fromkeys(l.chain for l in bb.labels):
        idx = [k for k, l in enumerate(bb.labels) if l.chain == chain]
        termini.append(
            {
                "chain": chain,
                "n": _side(float(depth[idx[0]]), mem.half),
                "c": _side(float(depth[idx[-1]]), mem.half),
            }
        )
    band = layout.membrane
    return {
        "source": mem.source,
        "sided": bool(mem.sided),
        "thickness": _num(2 * mem.half),
        "normal": _vec(mem.normal),
        "spans_slab": spans,
        "crosses_midplane": crosses,
        "elements": elements,
        "termini": termini,
        "page": {"y": [_num(band["y"][0]), _num(band["y"][1])], "inside": band["inside"]} if band else None,
    }


# ---------------------------------------------------------------------------------------------- schema


def _obj(props: dict, required: list[str] | None = None, nullable: bool = False) -> dict:
    out = {
        "type": ["object", "null"] if nullable else "object",
        "properties": props,
        "required": list(props) if required is None else required,
        "additionalProperties": False,
    }
    return out


def _t(kind, description: str = "", nullable: bool = False, **extra) -> dict:
    out = {"type": [kind, "null"] if nullable else kind, **extra}
    if description:
        out["description"] = description
    return out


def _arr(items: dict, description: str = "") -> dict:
    return _t("array", description, items=items)


RESIDUE = _obj(
    {
        "chain": _t("string"),
        "seq": _t("integer", "residue number (auth_seq_id)"),
        "icode": _t("string", "insertion code, '' when none"),
    }
)


def schema() -> dict:
    """The format of structure_data's dict, as a JSON Schema (draft 2020-12). Every section is optional (asked for
    by `sections`); within a section every key is present."""
    ids = _arr(_t("string"))
    structure = _obj(
        {
            "name": _t("string", "the structure's file name", nullable=True),
            "residues": _t("integer", "protein residues read (after --chains)"),
            "protein_chains": ids,
            "nucleic_acid_chains": ids,
            "chains_read": _arr(_t("string"), "every chain read and asked for, protein then nucleic acid"),
            "focus_chains": _t(
                "array", "the subunit drawn in full (--focus); null: everything", True, items=_t("string")
            ),
            "lent_chains": _arr(_t("string"), "chains drawn only as the elements they lend to the subunit's fold"),
            "context": _t("string", "what the figure shows of the assembly", nullable=True),
        }
    )
    page = _obj(
        {
            "x": _t("number"),
            "y": _t("number"),
            "angle": _t("number", "direction N->C, degrees counter-clockwise from page right, (-180, 180]"),
            "length": _t("number"),
            "width": _t("number"),
        },
        nullable=True,
    )
    element = _obj(
        {
            "id": _t("string", "CHAIN:FIRST-LAST, residue numbers with insertion codes"),
            "kind": _t("string", enum=sorted(KINDS.values())),
            "chain": _t("string"),
            "first": RESIDUE,
            "last": RESIDUE,
            "residues": _t("integer"),
            "label": _t("string", "the name the figure gives it", nullable=True),
            "page": page,
            "sheet": _t("string", "the sheet id it belongs to", nullable=True),
            "domain": _t("string", "the domain panel it is drawn in", nullable=True),
            "lent": _t("boolean", "a neighbour's element lent to the subunit's fold (--focus)"),
        }
    )
    pair = _obj(
        {
            "a": _t("string"),
            "b": _t("string"),
            "kind": _t("string", "P parallel, A antiparallel (the bridges' majority)", True, enum=["P", "A", None]),
            "bridges": _t("integer"),
            "register": _arr(_arr(RESIDUE), "each DSSP bridge: a's residue, then b's, in a's order"),
        }
    )
    sheet = _obj(
        {
            "id": _t("string", "S1, S2, ..."),
            "strands": _arr(_t("string"), "its strands side by side, left to right as drawn (as_drawn)"),
            "directions": _arr(_t("integer", enum=[1, -1]), "+1 running as the first strand, -1 against it"),
            "as_drawn": _t("boolean", "false: the figure does not draw the order (the bridges' order is given)"),
            "pairs": _arr(pair, "each pair of neighbours; a closed sheet's last pairing joins last and first"),
            "closed": _t("boolean"),
            "barrel": _obj(
                {
                    "n": _t("integer", "strands round the barrel"),
                    "shear": _t("integer", "shear number S, residues", nullable=True),
                },
                nullable=True,
            ),
            "accessory": _arr(
                _obj(
                    {
                        "strand": _t("string"),
                        "parent": _t("string"),
                        "direction": _t("integer", enum=[1, -1]),
                        "pairs_with": ids,
                    }
                ),
                "strands off the sheet's run of paired strands, each beside the strand it pairs with",
            ),
            "chords": _arr(ids, "strands of the run that pair but are not neighbours"),
            "joined": _arr(ids, "a strand of this sheet and one of another it pairs with"),
        }
    )
    disulfide = _obj(
        {"a": RESIDUE, "b": RESIDUE, "elements": _arr(_t("string", nullable=True), "the elements a and b are in")}
    )
    symmetry = _obj(
        {
            "found": _obj(
                {
                    "label": _t("string", "C2, D3, T, O, I or helical"),
                    "kind": _t("string", enum=["C", "D", "T", "O", "I", "H"]),
                    "n": _t("integer"),
                    "protomers": _arr(ids),
                    "unpaired_chains": ids,
                    "rmsd": _t("number", "CA RMSD between copies, Å"),
                    "axis": _arr(_t("number")),
                    "centre": _arr(_t("number")),
                    "twist": _t("number", "helical: degrees per copy", nullable=True),
                    "rise": _t("number", "helical: Å per copy", nullable=True),
                    "left_handed": _t("boolean", nullable=True),
                },
                nullable=True,
            ),
            "drawn": _t(
                "string",
                "how the figure draws the copies",
                True,
                enum=["rosette", "mirror", "one subunit", "side by side", None],
            ),
            "stoichiometry": _obj(
                {
                    "formula": _t("string", "protein chains lettered from the most copies, e.g. A2B2"),
                    "entities": _arr(
                        _obj({"letter": _t("string"), "chains": ids, "copies": _t("integer")}),
                    ),
                }
            ),
        }
    )
    sides = ["out", "in", "slab"]
    spans = "its CA run within 2 Å of both faces of the drawn slab (the rule A18 counts a transmembrane helix by)"
    crosses = "its CA run crosses the slab's centre plane and spans at least half the slab's thickness"
    membrane = _obj(
        {
            "source": _t("string", "opm (dummy atoms in the file) or estimate"),
            "sided": _t("boolean", "the positive-inside rule settled which side is out (non-cytoplasmic)"),
            "thickness": _t("number", "the drawn slab's, Å"),
            "normal": _arr(_t("number"), "unit vector pointing out"),
            "spans_slab": _arr(_t("string"), "the elements whose " + spans),
            "crosses_midplane": _arr(_t("string"), "the elements whose " + crosses),
            "elements": _arr(
                _obj(
                    {
                        "element": _t("string"),
                        "side": _t(
                            "string",
                            "across: it crosses the mid-plane (crosses_midplane); else out, in or slab by its mean "
                            "CA depth",
                            enum=["across", *sides],
                        ),
                        "from": _t("string", "the side its first CA is on", enum=["out", "in"]),
                        "to": _t("string", "the side its last CA is on", enum=["out", "in"]),
                        "depth": _arr(_t("number"), "its CAs' least and greatest depth, Å from the mid-plane, + out"),
                        "spans_slab": _t("boolean", spans),
                        "crosses_midplane": _t("boolean", crosses),
                    },
                )
            ),
            "termini": _arr(
                _obj({"chain": _t("string"), "n": _t("string", enum=sides), "c": _t("string", enum=sides)})
            ),
            "page": _obj(
                {"y": _arr(_t("number")), "inside": _t("string", enum=["above", "below"])},
                nullable=True,
            ),
        },
        nullable=True,
    )
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "foldmap data",
        "description": "What a foldmap figure says about a structure (foldmap data). Numbers rounded to 0.001; page "
        "places in page units, y up.",
        "type": "object",
        "properties": {
            "foldmap_data": _t("integer", "the format's version", const=SCHEMA_VERSION),
            "foldmap": _t("string", "the foldmap version that wrote it"),
            "structure": structure,
            "elements": _arr(element, "in sequence order"),
            "sheets": _arr(sheet),
            "disulfides": _arr(disulfide),
            "symmetry": symmetry,
            "membrane": membrane,
        },
        "required": ["foldmap_data", "foldmap"],
        "additionalProperties": False,
    }
