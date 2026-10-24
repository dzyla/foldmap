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
  cross its mid-plane, and the sides of elements and termini;
- ligands: each ligand, metal ion, glycan and residue fused into the chain (a chromophore), the residues it touches
  (io's ligand rule) and the elements they are in;
- interfaces: each pair of chains in contact (protein chains by CA pairs, a protein and a nucleic-acid chain by atom
  contacts, two nucleic-acid chains by base pairs), with their elements in contact (features.sse_contacts) and,
  with the interactions section, the interactions between them by element pair;
- nucleic_acid: the strands the figure draws, their sequences and base pairs (nucleic.py's, the pairs the figure
  draws) with the stem each is in and its bracket level (nucleic._knots' rule, extended to every pair), written as
  dot-bracket, BPSEQ and CT by `dot_bracket_text`, `bpseq_text` and `ct_text`;
- interactions (only when asked for): every interaction 13c's finders find between two chains, re-checked as the
  figure's lines are (interactions.failing_lines), with the element at each end.

Every list has a fixed order, every number is rounded to 0.001 (±0 written 0.0, angles in (-180, 180]), so the
same structure and options give the same bytes (A12). `schema()` is the format, as a JSON Schema; SCHEMA_VERSION
is its version, the dict's "foldmap_data" key. `tables()` gives the list-shaped parts as flat tables (TABLES: one
row per element, sheet strand, pairing, bridge, disulfide, membrane element, ligand, ligand contact, interface,
interface element pair, interaction, nucleotide, base pair), written as CSV by `csv_text` and `write_tables`.
`spec_data` is the whole route from a figure spec (cli.make_layout on io.read_model's model), the one the command
line, the MCP server and the app share."""

from __future__ import annotations

import csv
import io
import json
import re
from collections import Counter
from pathlib import Path
from statistics import median_low

import numpy as np

SCHEMA_VERSION = 1
SECTIONS = (
    "structure",
    "elements",
    "sheets",
    "disulfides",
    "symmetry",
    "membrane",
    "ligands",
    "interfaces",
    "nucleic_acid",
    "interactions",
)
ALL = SECTIONS
DEFAULT = tuple(s for s in SECTIONS if s != "interactions")  # interactions run 13c's finders: asked for by name
INTERACTION_KINDS = ("hbond", "salt", "stack", "cation-pi", "hydrophobic")  # interactions.KINDS (a test checks)
BETWEEN = "interchain"  # the interactions section's partners (hbonds.parse_selection): residues of two chains
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


def structure_data(
    layout, sses, bb, model=None, sections: tuple[str, ...] = DEFAULT, *, name: str | None = None
) -> dict:
    """The data of a laid-out structure (cli.make_layout's layout, sses and backbone; `model`: io.read_model's, the
    one the layout was built from, needed by the ligands and interactions sections), the sections asked for
    (SECTIONS, by default DEFAULT: all but interactions; ValueError naming them for another). `name`: the structure
    as the user named it (written as its file name)."""
    unknown = [s for s in sections if s not in SECTIONS]
    if unknown:
        raise ValueError(f"no data section {', '.join(unknown)}; the sections are {', '.join(SECTIONS)}")
    needs = [s for s in ("ligands", "interactions") if s in sections]
    if model is None and needs:
        raise ValueError(f"the {' and '.join(needs)} section{'s' if len(needs) > 1 else ''} need the model (model=)")
    if isinstance(model, tuple):  # io.read_model's ModelRead: (gemmi.Model, chains)
        model = model[0]
    from . import __version__

    out: dict = {"foldmap_data": SCHEMA_VERSION, "foldmap": __version__}
    sheets = _sheets(layout, sses, bb) if {"elements", "sheets"} & set(sections) else []
    found = _interactions(layout, sses, bb, model) if "interactions" in sections else None
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
        elif section == "membrane":
            out[section] = _membrane(layout, sses, bb)
        elif section == "ligands":
            out[section] = _ligands(layout, sses, bb, model)
        elif section == "interfaces":
            out[section] = chain_interfaces(
                bb, model, sses=sses, nucleic=layout.nucleic, chains_read=layout.chains_read, found=found
            )
        elif section == "nucleic_acid":
            na = layout.nucleic
            drawn = na is not None and na.strands
            out[section] = nucleic_acid(na, cropped=bool(layout.na_cropped)) if drawn else None
        else:
            out[section] = found
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


def spec_data(spec, sections: tuple[str, ...] = DEFAULT, name: str | None = None) -> dict:
    """The data of the figure a spec.FigureSpec describes: its structure read (io.read_model, its assembly) and laid
    out as `plot` lays it out (cli.make_layout: no routing, no drawing), then structure_data. The one route the
    command line, the MCP server and the app take, so they give the same bytes. `name`: the structure's name in the
    data (default its file's; the app names an uploaded structure as uploaded, not by the copy it keeps)."""
    from .cli import make_layout
    from .io import read_model

    model = read_model(spec.structure, spec.assembly)
    layout, sses, bb = make_layout(
        spec.structure, spec.mode, spec.rotate, spec.flip_v, model=model, **spec.layout_kwargs()
    )
    return structure_data(layout, sses, bb, model, sections, name=name or str(spec.structure))


def zip_bytes(data: dict) -> bytes:
    """Every table of `data` as NAME.csv, and with a nucleic acid its dot-bracket, BPSEQ and CT (nucleic_acid.dbn,
    .bpseq, .ct), in one zip: the same bytes for the same data (fixed time stamps, TABLES order)."""
    import zipfile

    files = [(f"{name}.csv", csv_text(columns, rows)) for name, (columns, rows) in tables(data).items()]
    if data.get("nucleic_acid") is not None:
        files += [(f"nucleic_acid.{ext}", text_of(data)) for ext, text_of in NA_FORMATS.items()]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, text in files:
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, text.encode("utf-8"))
    return buf.getvalue()


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


def strand_pairing(a, b, bridges, labels) -> dict:
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
        pairs = [strand_pairing(by_id[a], by_id[b], bridges, bb.labels) for a, b in steps]
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
    drawn slab (A18's rule before plan 26); `crosses_midplane`, its CA run crosses the slab's centre plane and spans
    at least half the slab's thickness (membrane.crosses: the usual reading of a transmembrane segment, lenient to a
    slab estimated too thick or tilted, and the rule the membrane view and A18 count a transmembrane element by).
    Its side: "across" when it crosses the mid-plane, else out, in or slab by its mean CA depth; `from` and `to`, the
    sides of its first and last CA. Sides are named by the membrane's normal, which points out (the non-cytoplasmic
    side; `sided` says whether the positive-inside rule settled which side that is)."""
    mem = getattr(layout, "bilayer", None)
    if mem is None:
        return None
    from .membrane import FACE_REACH
    from .membrane import crosses as crossing

    depth = mem.depth(bb.ca)

    def sign(r: int) -> str:
        return "out" if depth[r] >= 0 else "in"

    elements, spans, crosses = [], [], []
    for s in sses:
        seg = depth[s.start : s.end + 1]
        lo, hi = float(seg.min()), float(seg.max())
        span = len(seg) >= 2 and lo <= -mem.half + FACE_REACH and hi >= mem.half - FACE_REACH
        cross = crossing(seg, mem.half)  # mid-plane crossed, half the thickness spanned: the view's and A18's rule
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


def _owners(bb, sses) -> list[str | None]:
    """The element each backbone residue is in (None: a loop)."""
    owner: list[str | None] = [None] * len(bb)
    for s in sses:
        owner[s.start : s.end + 1] = [s.id] * (s.end + 1 - s.start)
    return owner


def _where(bb) -> dict[tuple[str, int, str], int]:
    return {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}


def _heavy(res) -> list:
    return [a for a in res if a.element.name not in ("H", "D")]


def _bonded(search, model, heavy, where: dict, own: set) -> set[int]:
    """The backbone residues covalently joined to these atoms: a heavy atom within io._BOND (the rule io joins the
    parts of a ligand by), the group's own residues `own` ((chain, number, icode)) left out."""
    from .io import _BOND

    out = set()
    for atom in heavy:
        for mark in search.find_atoms(atom.pos, radius=_BOND):
            cra = mark.to_cra(model)
            key = (cra.chain.name, cra.residue.seqid.num, cra.residue.seqid.icode.strip())
            if key in where and key not in own and cra.atom.element.name not in ("H", "D"):
                if cra.atom.pos.dist(atom.pos) < _BOND:
                    out.add(where[key])
    return out


def _ligands(layout, sses, bb, model) -> list[dict]:
    """What is bound to the protein (layout.links, as io.load_links found it for the figure): each ligand and metal
    ion (io's Ligand, its contacts as io found them: any atom within 4.0 Å, a metal's O, N or S within 3.1 Å), each
    glycan (its sugars, by the same ligand rule through io.touched) and each residue fused into the chain (a
    chromophore: Backbone.linked; the ligand rule, itself left out). For each: the residues it touches, the element
    each is in, the elements in sequence order, and the residues covalently joined to it (attached_to)."""
    links = layout.links
    if links is None:
        return []
    import gemmi

    from .io import _ligand_label, touched

    owner, where = _owners(bb, sses), _where(bb)
    order = {s.id: n for n, s in enumerate(sses)}
    search = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate()
    index: dict[tuple[str, int, str], object] = {}
    for chain in model:
        for res in chain:
            index.setdefault((chain.name, res.seqid.num, res.name), (chain.name, res))

    def entry(kind, name, codes, symbol, parts, at, contacts, attached, additive=False, nucleotides=()) -> dict:
        contacts = sorted(contacts)
        elements = sorted({owner[k] for k in contacts if owner[k] is not None}, key=order.__getitem__)
        return {
            "id": f"{kind}:{at['chain']}:{at['seq']}{at['icode']}",
            "kind": kind,
            "name": name,
            "codes": list(codes),
            "symbol": symbol,
            "at": at,
            "residues": [{"chain": c, "seq": int(r.seqid.num), "icode": r.seqid.icode.strip()} for c, r in parts],
            "additive": bool(additive),
            "attached_to": [_res(bb.labels[k]) for k in sorted(attached)],
            "contacts": [{"residue": _res(bb.labels[k]), "element": owner[k]} for k in contacts],
            "elements": elements,
            "nucleotides": [{"chain": c, "seq": int(n), "icode": i} for c, n, i in sorted(nucleotides)],
        }

    def group(members, codes) -> tuple[list, set, set]:
        """A group's residues in the model, the residues its atoms touch and those joined to it."""
        parts = [index[(c, n, code)] for (c, n), code in zip(members, codes) if (c, n, code) in index]
        own = {(c, r.seqid.num, r.seqid.icode.strip()) for c, r in parts}
        contacts, attached = set(), set()
        for _, res in parts:
            contacts |= touched(search, model, _heavy(res), False, where)[0]
            attached |= _bonded(search, model, _heavy(res), where, own)
        return parts, contacts, attached

    out = []
    for g in links.ligands:
        codes = g.codes or [g.name]
        parts = [index[(c, n, code)] for (c, n), code in zip(g.members, codes) if (c, n, code) in index]
        own = {(c, r.seqid.num, r.seqid.icode.strip()) for c, r in parts}
        attached = set()
        if not g.metal:  # io: a metal ion bonds to nothing
            for _, res in parts:
                attached |= _bonded(search, model, _heavy(res), where, own)
        icode = parts[0][1].seqid.icode.strip() if parts else ""
        at = {"chain": g.chain, "seq": int(g.seq), "icode": icode}
        kind = "metal" if g.metal else "ligand"
        out.append(entry(kind, g.name, codes, g.symbol, parts, at, g.contacts, attached, g.additive, g.nucleotides))
    for (_, sugars), members in zip(links.glycans, links.glycan_members):
        parts, contacts, attached = group(members, sugars)
        if not parts:
            continue
        c, root = parts[0]
        at = {"chain": c, "seq": int(root.seqid.num), "icode": root.seqid.icode.strip()}
        name = sugars[0] if len(sugars) == 1 else "+".join(dict.fromkeys(sugars))
        out.append(entry("glycan", name, sugars, _ligand_label(sugars), parts, at, contacts, attached))
    for k in sorted(bb.linked):
        lab = bb.labels[k]
        res = next(
            (
                r
                for chain in model
                if chain.name == lab.chain
                for r in chain
                if r.seqid.num == lab.seq and r.seqid.icode.strip() == lab.icode and r.name == lab.name
            ),
            None,
        )
        if res is None:
            continue
        own = {(lab.chain, lab.seq, lab.icode)}
        contacts = touched(search, model, _heavy(res), False, where)[0] - {k}
        attached = _bonded(search, model, _heavy(res), where, own)
        out.append(entry("linked", lab.name, [lab.name], lab.name, [(lab.chain, res)], _res(lab), contacts, attached))
    return out


def _counts(kinds: Counter | None) -> dict | None:
    return None if kinds is None else {k: int(kinds[k]) for k in INTERACTION_KINDS}


def _na_faces(na, chains_read, bb, model) -> tuple[dict, dict]:
    """The nucleic acid's chain pairs in contact: (protein chain, nucleic-acid chain) -> (residue-nucleotide pairs
    with a heavy atom of each within io.CONTACT, 4 Å: the rule the figure's contacts are read by, io.na_contacts;
    the protein residues and the nucleotides in them), measured on the structure (`model`: every nucleotide of the
    chains read, cropped or not; without it, the nucleic acid the figure draws); and (nucleic-acid chain,
    nucleic-acid chain) -> (base pairs between them, the nucleotides of each in them), the base pairs the figure
    draws (nucleic.py's). `na`: the nucleic acid the figure draws (None: none); `chains_read`: the chains read
    (none: every one)."""
    protein: dict[tuple[str, str], list] = {}
    if model is not None:
        import gemmi

        from .io import _nucleic_strands, na_contacts

        read = set(chains_read)
        found = [(st, run) for st, run in _nucleic_strands(model) if not read or st.chain in read]
        strand_chain = [st.chain for st, _ in found]
        contacts = {}
        if found and len(bb):
            ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate()
            contacts = na_contacts(ns, model, [run for _, run in found], _where(bb))
    else:
        strand_chain = [st.chain for st in na.strands] if na is not None else []
        contacts = na.contacts if na is not None else {}
    for k, touched in contacts.items():
        for s, i in touched:
            face = protein.setdefault((bb.labels[k].chain, strand_chain[s]), [0, set(), set()])
            face[0] += 1
            face[1].add(k)
            face[2].add((s, i))
    paired: dict[tuple[str, str], list] = {}
    for p in na.base_pairs if na is not None else []:
        ca, cb = na.strands[p.a].chain, na.strands[p.b].chain
        if ca != cb:
            face = paired.setdefault((ca, cb), [0, set(), set()])
            face[0] += 1
            face[1].add((p.a, p.i))
            face[2].add((p.b, p.j))
    return protein, paired


def chain_interfaces(bb, model=None, *, sses=(), nucleic=None, chains_read=(), found: dict | None = None) -> list[dict]:
    """Each pair of chains in contact, by its measure: two protein chains of `bb` with CAs within
    features.CONTACT_CUTOFF (8 Å), "ca"; a protein and a nucleic-acid chain with a residue and a nucleotide in
    contact (_na_faces: heavy atoms within 4 Å, on the structure `model`, else on `nucleic`), "atoms"; two
    nucleic-acid chains base-paired (`nucleic`'s pairs: the figure's), "base_pairs"; and with the interactions
    section any two chains with an interaction between them. For each, the pairs in contact by its measure and the
    residues (nucleotides) of each chain in them, the CA pairs (protein chains), its element pairs among `sses`
    (features.sse_contacts: CA pairs within 8 Å, listed from MIN_CONTACT_PAIRS), and with `found` (the interactions
    section) the interactions by kind, in all and by element pair (an end outside the elements: null). Chains are
    ranked as `chains_read` (the figure's) lists them, then as `bb`, then the nucleic acid's."""
    from scipy.spatial import cKDTree

    from .features import CONTACT_CUTOFF, sse_contacts

    chains = [l.chain for l in bb.labels]
    protein_chains = set(chains)
    na = nucleic
    na_chains = [st.chain for st in na.strands] if na is not None else []
    rank = {c: n for n, c in enumerate(dict.fromkeys([*chains_read, *chains, *na_chains]))}
    order = {s.id: n for n, s in enumerate(sses)}
    chain_of = {s.id: s.chain for s in sses}
    pairs: dict[tuple[str, str], dict] = {}

    def measure(a: str, b: str) -> str:
        n = (a in protein_chains) + (b in protein_chains)
        return "ca" if n == 2 else "atoms" if n == 1 else "base_pairs"

    def at(a: str, b: str) -> tuple[dict, bool]:
        """The interface of chains a and b, and whether a is its second chain."""
        swap = rank.get(a, len(rank)) > rank.get(b, len(rank)) or (a > b and rank.get(a) == rank.get(b))
        key = (b, a) if swap else (a, b)
        if key not in pairs:
            ca = 0 if measure(*key) == "ca" else None
            pairs[key] = {"contacts": 0, "ca_pairs": ca, "residues": [0, 0], "elements": {}, "kinds": Counter()}
        return pairs[key], swap

    if len(bb) >= 2:
        code = np.array([rank[c] for c in chains], int)
        ij = cKDTree(bb.ca).query_pairs(CONTACT_CUTOFF, output_type="ndarray")
        if len(ij):
            i, j = ij[:, 0], ij[:, 1]
            cross = code[i] != code[j]
            i, j = i[cross], j[cross]
            low = code[i] < code[j]
            i, j = np.where(low, i, j), np.where(low, j, i)
            for ci, cj in sorted({(int(a), int(b)) for a, b in zip(code[i], code[j])}):
                mask = (code[i] == ci) & (code[j] == cj)
                face, _ = at(chains[int(i[mask][0])], chains[int(j[mask][0])])
                face["ca_pairs"] = face["contacts"] = int(mask.sum())
                face["residues"] = [len(np.unique(i[mask])), len(np.unique(j[mask]))]
    for faces in _na_faces(nucleic, chains_read, bb, model):
        for (a, b), (n, mine, theirs) in sorted(faces.items()):
            face, swap = at(a, b)
            face["contacts"] = n
            face["residues"] = [len(theirs), len(mine)] if swap else [len(mine), len(theirs)]
    for key, n in sse_contacts(bb, sses).items():
        a, b = sorted(key, key=order.__getitem__)
        if chain_of[a] == chain_of[b]:
            continue
        face, swap = at(chain_of[a], chain_of[b])
        ends = (b, a) if swap else (a, b)
        face["elements"].setdefault(ends, {"ca_pairs": None, "kinds": Counter()})["ca_pairs"] = int(n)
    for x in (found or {}).get("found", []):
        face, swap = at(x["a"]["chain"], x["b"]["chain"])
        ends = tuple(x["elements"][::-1] if swap else x["elements"])
        face["kinds"][x["kind"]] += 1
        face["elements"].setdefault(ends, {"ca_pairs": None, "kinds": Counter()})["kinds"][x["kind"]] += 1
    big = len(order)

    def place(e):
        return (e is None, order.get(e, big))

    out = []
    for (a, b), face in sorted(pairs.items(), key=lambda kv: (rank.get(kv[0][0], big), rank.get(kv[0][1], big))):
        out.append(
            {
                "chains": [a, b],
                "measure": measure(a, b),
                "contacts": face["contacts"],
                "ca_pairs": face["ca_pairs"],
                "residues": face["residues"],
                "elements": [
                    {
                        "a": ea,
                        "b": eb,
                        "ca_pairs": pair["ca_pairs"],
                        "interactions": _counts(pair["kinds"] if found is not None else None),
                    }
                    for (ea, eb), pair in sorted(
                        face["elements"].items(), key=lambda kv: (*place(kv[0][0]), *place(kv[0][1]))
                    )
                ],
                "interactions": _counts(face["kinds"] if found is not None else None),
            }
        )
    return out


def _interactions(layout, sses, bb, model) -> dict:
    """Every interaction 13c's finders find between residues of two chains (hbonds.find_hbonds and the
    interactions finders, BETWEEN), among the chains read; each re-checked by interactions.failing_lines as the
    figure's lines are, a flagged one left out and named. Each with its partners' residues and atoms (an H-bond:
    donor then acceptor; else as interactions.Interaction), and the element each end is in (null: a loop or a
    nucleotide). In interactions.sort_key order."""
    from .hbonds import find_hbonds
    from .interactions import KINDS, _partner_key, failing_lines, sort_key
    from .render import INTERACTION_FINDERS, _as_interaction

    assert KINDS == INTERACTION_KINDS
    found = [_as_interaction(b) for b in find_hbonds(model, BETWEEN)]
    for kind in KINDS:
        if kind != "hbond":
            found += INTERACTION_FINDERS[kind](model, BETWEEN)
    read = set(layout.chains_read)
    found = [x for x in found if x.a[0] in read and x.b[0] in read]
    failing = failing_lines(model, found)
    flagged = {n for n, _ in failing}
    kept = sorted((x for n, x in enumerate(found) if n not in flagged), key=sort_key)
    owner, where = _owners(bb, sses), _where(bb)

    def end(p) -> tuple[dict, str | None]:
        chain, seq, icode, atom = _partner_key(p)
        k = where.get((chain, seq, icode))
        return {"chain": chain, "seq": int(seq), "icode": icode, "atom": atom}, (owner[k] if k is not None else None)

    rows = []
    for x in kept:
        (a, ea), (b, eb) = end(x.a), end(x.b)
        rows.append(
            {
                "kind": x.kind,
                "a": a,
                "b": b,
                "distance": _num(x.distance),
                "angle": _num(x.angle) if x.angle is not None else None,
                "note": x.note,
                "elements": [ea, eb],
            }
        )
    return {
        "between": BETWEEN,
        "kinds": list(KINDS),
        "found": rows,
        "failed_recheck": [text for _, text in failing],
    }


# ---------------------------------------------------------------------------------------------- nucleic acid

# bracket pairs by level: level 0 the nested pairs, 1 the first pseudoknot level, then upper- and lower-case
# letters (the extended dot-bracket RNApolis and others read)
BRACKETS = (("(", ")"), ("[", "]"), ("{", "}"), ("<", ">"), *((chr(c), chr(c + 32)) for c in range(65, 91)))


def _crosses(x: tuple[int, int], placed: list[tuple[int, int]]) -> bool:
    x0, x1 = x
    return any(x0 < y0 < x1 < y1 or y0 < x0 < y1 < x1 for y0, y1 in placed)


def pair_levels(na) -> dict[tuple[int, int, int, int], int]:
    """The bracket level of each base pair of a model.Nucleic ((a, i, b, j) -> 0, 1, ...), nucleotides in strand
    order (the strands one after another, as they are listed). foldmap's pseudoknot rule (nucleic._knots, the one
    that names a stem PK) decides level 0: every stem in one chain it keeps nested. The rest, in the order that rule
    takes stems (longest first, of equal ones the 5' one; then the isolated pairs, 5' first), each go whole to the
    lowest level where none of their pairs crosses a pair already there: a PK stem to 1 ([]), one crossing it too
    to 2 ({}), a stem between two chains (a duplex, a codon and anticodon) to 0 unless it crosses a pair there."""
    from .nucleic import _knots

    start, pos = 0, {}
    for s, st in enumerate(na.strands):
        pos[s] = start
        start += len(st)
    stems, strands = na.structure.stems, na.strands
    knots = _knots(stems, strands)

    def spans(pairs) -> list[tuple[int, int]]:
        return [(pos[p.a] + p.i, pos[p.b] + p.j) for p in pairs]

    levels: list[list[tuple[int, int]]] = [[]]
    out: dict[tuple[int, int, int, int], int] = {}
    rest = []
    for k, st in enumerate(stems):
        if k not in knots and strands[st.a].chain == strands[st.b].chain:  # nested by the rule: never crossing
            levels[0] += spans(st.pairs)
            out.update({tuple(p[:4]): 0 for p in st.pairs})
        else:
            rest.append(st.pairs)
    rest = sorted(rest, key=lambda pairs: -len(pairs)) + [[p] for p in sorted(na.structure.ties)]
    for pairs in rest:
        mine = spans(pairs)
        level = next((n for n, placed in enumerate(levels) if not any(_crosses(x, placed) for x in mine)), None)
        if level is None:
            level = len(levels)
            levels.append([])
        if level >= len(BRACKETS):
            raise ValueError(f"the base pairs need more than {len(BRACKETS)} bracket levels")
        levels[level] += mine
        out.update({tuple(p[:4]): level for p in pairs})
    return out


def nucleic_acid(na, cropped: bool = False) -> dict:
    """The nucleic-acid section of a model.Nucleic (the figure's: layout.nucleic): its strands (each a strand piece
    as drawn, 5'->3', with its residues, sequence and dot-bracket), its base pairs (nucleic.py's, numbered 1.. in
    strand order, with their stem and bracket level: pair_levels), stems, gaps and G-tetrads. `cropped`: na_extent
    cropped it to what the protein touches."""
    from .nucleic import _knots

    levels = pair_levels(na)
    knots = _knots(na.structure.stems, na.strands)
    stem_of: dict[tuple, dict] = {}
    stems = []
    for k, st in enumerate(na.structure.stems):
        entry = {
            "id": st.id,
            "name": st.name,
            "pairs": len(st),
            "level": levels[tuple(st.pairs[0][:4])],
            "pseudoknot": k in knots,
        }
        stems.append(entry)
        stem_of.update({tuple(p[:4]): entry for p in st.pairs})
    start, pos = 1, {}
    for s, st in enumerate(na.strands):
        pos[s] = start
        start += len(st)
    marks = ["."] * (start - 1)
    pairs = []
    for p in na.base_pairs:
        key = tuple(p[:4])
        i, j = pos[p.a] + p.i, pos[p.b] + p.j
        opening, closing = BRACKETS[levels[key]]
        marks[i - 1], marks[j - 1] = opening, closing
        stem = stem_of.get(key)
        pairs.append(
            {
                "i": i,
                "j": j,
                "a": _res(na.strands[p.a].labels[p.i]),
                "b": _res(na.strands[p.b].labels[p.j]),
                "kind": p.kind,
                "stem": stem["name"] if stem else None,
                "level": levels[key],
            }
        )
    strands = []
    for s, st in enumerate(na.strands):
        strands.append(
            {
                "id": st.id,
                "chain": st.chain,
                "rna": bool(st.rna),
                "first": pos[s],
                "residues": [_res(label) for label in st.labels],
                "sequence": st.bases,
                "dot_bracket": "".join(marks[pos[s] - 1 : pos[s] - 1 + len(st)]),
                "cut5": bool(st.cut5),
                "cut3": bool(st.cut3),
            }
        )
    return {
        "nucleotides": start - 1,
        "cropped": bool(cropped),
        "strands": strands,
        "pairs": pairs,
        "levels": 1 + max(levels.values(), default=0),
        "stems": stems,
        "gaps": [
            {"before": na.strands[g.before].id, "after": na.strands[g.after].id, "nucleotides": int(g.count)}
            for g in na.gaps
        ],
        "tetrads": [[_res(na.strands[s].labels[i]) for s, i in t] for t in na.tetrads],
    }


def _section(data: dict) -> dict:
    na = data.get("nucleic_acid")
    if na is None:
        raise ValueError("no nucleic acid is drawn: dot-bracket, BPSEQ and CT need one (or the nucleic_acid section)")
    return na


def _partners(na: dict) -> list[int]:
    partner = [0] * (na["nucleotides"] + 1)
    for p in na["pairs"]:
        partner[p["i"]], partner[p["j"]] = p["j"], p["i"]
    return partner


def dot_bracket_text(data: dict) -> str:
    """The nucleic acid as dot-bracket: a '>' line naming each strand (CHAIN:FIRST-LAST), its sequence, its
    structure; a pair between strands bracketed across them (the strands read one after another), pseudoknot levels
    as [] {} <> then Aa Bb ... (BRACKETS). The multi-strand form RNApolis reads."""
    na = _section(data)
    lines = []
    for st in na["strands"]:
        lines += [f">{st['id']}", st["sequence"], st["dot_bracket"]]
    return "\n".join(lines) + "\n"


def bpseq_text(data: dict) -> str:
    """The nucleic acid as BPSEQ: one line per nucleotide, its number (1.. through every strand in order), base and
    partner's number (0 unpaired)."""
    na = _section(data)
    partner = _partners(na)
    seq = "".join(st["sequence"] for st in na["strands"])
    return "".join(f"{n} {seq[n - 1]} {partner[n]}\n" for n in range(1, na["nucleotides"] + 1))


def ct_text(data: dict) -> str:
    """The nucleic acid as CT (connect table): a header with the nucleotide count and the structure's name, then one
    line per nucleotide: its number, base, the numbers before and after it (0 at a strand's end: a new strand, a
    chain break or a cropped stretch), its partner's (0 unpaired) and its residue number in the file (the insertion
    code is in the JSON)."""
    na = _section(data)
    partner = _partners(na)
    name = (data.get("structure") or {}).get("name") or "foldmap"
    lines = [f"{na['nucleotides']} {name}"]
    for st in na["strands"]:
        first, last = st["first"], st["first"] + len(st["residues"]) - 1
        for n, (base, res) in enumerate(zip(st["sequence"], st["residues"]), first):
            before = n - 1 if n > first else 0
            after = n + 1 if n < last else 0
            lines.append(f"{n} {base} {before} {after} {partner[n]} {res['seq']}")
    return "\n".join(lines) + "\n"


NA_FORMATS = {"dbn": dot_bracket_text, "bpseq": bpseq_text, "ct": ct_text}  # foldmap data --format, by suffix


# ---------------------------------------------------------------------------------------------- tables


def _ref(r: dict) -> str:
    return f"{r['chain']}:{r['seq']}{r['icode']}"


def _kind_cells(counts: dict | None) -> list:
    return [counts[k] if counts is not None else None for k in INTERACTION_KINDS]


def _t_elements(d: dict) -> list[list]:
    rows = []
    for e in d["elements"]:
        page = e["page"] or {}
        rows.append(
            [
                e["id"],
                e["kind"],
                e["chain"],
                e["first"]["seq"],
                e["first"]["icode"],
                e["last"]["seq"],
                e["last"]["icode"],
                e["residues"],
                e["label"],
                *(page.get(k) for k in ("x", "y", "angle", "length", "width")),
                e["sheet"],
                e["domain"],
                e["lent"],
            ]
        )
    return rows


def _t_sheets(d: dict) -> list[list]:
    rows = []
    for sh in d["sheets"]:
        barrel = sh["barrel"] or {}
        tail = [sh["as_drawn"], sh["closed"], barrel.get("n"), barrel.get("shear")]
        for n, (strand, sense) in enumerate(zip(sh["strands"], sh["directions"]), 1):
            rows.append([sh["id"], strand, n, "run", sense, None, *tail])
        for a in sh["accessory"]:
            rows.append([sh["id"], a["strand"], None, "accessory", a["direction"], a["parent"], *tail])
    return rows


def _t_sheet_pairs(d: dict) -> list[list]:
    return [[sh["id"], p["a"], p["b"], p["kind"], p["bridges"]] for sh in d["sheets"] for p in sh["pairs"]]


def _t_sheet_register(d: dict) -> list[list]:
    return [
        [sh["id"], p["a"], p["b"], ra["chain"], ra["seq"], ra["icode"], rb["chain"], rb["seq"], rb["icode"]]
        for sh in d["sheets"]
        for p in sh["pairs"]
        for ra, rb in p["register"]
    ]


def _t_disulfides(d: dict) -> list[list]:
    return [
        [
            s["a"]["chain"],
            s["a"]["seq"],
            s["a"]["icode"],
            s["b"]["chain"],
            s["b"]["seq"],
            s["b"]["icode"],
            *s["elements"],
        ]
        for s in d["disulfides"]
    ]


def _t_membrane(d: dict) -> list[list]:
    mem = d["membrane"]
    if mem is None:
        return []
    return [
        [e["element"], e["side"], e["from"], e["to"], *e["depth"], e["spans_slab"], e["crosses_midplane"]]
        for e in mem["elements"]
    ]


def _t_ligands(d: dict) -> list[list]:
    return [
        [
            g["id"],
            g["kind"],
            g["name"],
            " ".join(g["codes"]),
            g["symbol"],
            g["at"]["chain"],
            g["at"]["seq"],
            g["at"]["icode"],
            g["additive"],
            " ".join(_ref(r) for r in g["attached_to"]),
            len(g["contacts"]),
            " ".join(g["elements"]),
            len(g["nucleotides"]),
        ]
        for g in d["ligands"]
    ]


def _t_ligand_contacts(d: dict) -> list[list]:
    return [
        [g["id"], g["kind"], g["name"], c["residue"]["chain"], c["residue"]["seq"], c["residue"]["icode"], c["element"]]
        for g in d["ligands"]
        for c in g["contacts"]
    ]


def _t_interfaces(d: dict) -> list[list]:
    return [
        [
            *f["chains"],
            f["measure"],
            f["contacts"],
            f["ca_pairs"],
            *f["residues"],
            len(f["elements"]),
            *_kind_cells(f["interactions"]),
        ]
        for f in d["interfaces"]
    ]


def _t_interface_elements(d: dict) -> list[list]:
    return [
        [*f["chains"], e["a"], e["b"], e["ca_pairs"], *_kind_cells(e["interactions"])]
        for f in d["interfaces"]
        for e in f["elements"]
    ]


def _t_interactions(d: dict) -> list[list]:
    return [
        [
            x["kind"],
            *(x["a"][k] for k in ("chain", "seq", "icode", "atom")),
            *(x["b"][k] for k in ("chain", "seq", "icode", "atom")),
            x["distance"],
            x["angle"],
            x["note"],
            *x["elements"],
        ]
        for x in d["interactions"]["found"]
    ]


def _t_nucleotides(d: dict) -> list[list]:
    na = d["nucleic_acid"]
    if na is None:
        return []
    partner = _partners(na)
    level = {n: p["level"] for p in na["pairs"] for n in (p["i"], p["j"])}
    stem = {n: p["stem"] for p in na["pairs"] for n in (p["i"], p["j"])}
    rows = []
    for st in na["strands"]:
        for n, (base, mark, r) in enumerate(zip(st["sequence"], st["dot_bracket"], st["residues"]), st["first"]):
            rows.append(
                [
                    n,
                    st["id"],
                    r["chain"],
                    r["seq"],
                    r["icode"],
                    base,
                    partner[n] or None,
                    mark,
                    level.get(n),
                    stem.get(n),
                ]
            )
    return rows


def _t_base_pairs(d: dict) -> list[list]:
    na = d["nucleic_acid"]
    if na is None:
        return []
    return [
        [p["i"], p["j"], *(p["a"][k] for k in _RES), *(p["b"][k] for k in _RES), p["kind"], p["stem"], p["level"]]
        for p in na["pairs"]
    ]


_RES = ("chain", "seq", "icode")
TABLES: dict[str, tuple[str, tuple[str, ...], object]] = {  # name -> (its section, its columns, its rows)
    "elements": (
        "elements",
        (
            "id",
            "kind",
            "chain",
            "first_seq",
            "first_icode",
            "last_seq",
            "last_icode",
            "residues",
            "label",
            "x",
            "y",
            "angle",
            "length",
            "width",
            "sheet",
            "domain",
            "lent",
        ),
        _t_elements,
    ),
    "sheets": (
        "sheets",
        ("sheet", "strand", "position", "role", "direction", "parent", "as_drawn", "closed", "barrel_n", "shear"),
        _t_sheets,
    ),
    "sheet_pairs": ("sheets", ("sheet", "a", "b", "kind", "bridges"), _t_sheet_pairs),
    "sheet_register": (
        "sheets",
        ("sheet", "a", "b", *(f"a_{k}" for k in _RES), *(f"b_{k}" for k in _RES)),
        _t_sheet_register,
    ),
    "disulfides": (
        "disulfides",
        (*(f"a_{k}" for k in _RES), *(f"b_{k}" for k in _RES), "a_element", "b_element"),
        _t_disulfides,
    ),
    "membrane": (
        "membrane",
        ("element", "side", "from", "to", "depth_min", "depth_max", "spans_slab", "crosses_midplane"),
        _t_membrane,
    ),
    "ligands": (
        "ligands",
        (
            "id",
            "kind",
            "name",
            "codes",
            "symbol",
            "chain",
            "seq",
            "icode",
            "additive",
            "attached_to",
            "contacts",
            "elements",
            "nucleotides",
        ),
        _t_ligands,
    ),
    "ligand_contacts": ("ligands", ("ligand", "kind", "name", *_RES, "element"), _t_ligand_contacts),
    "interfaces": (
        "interfaces",
        (
            "chain_a",
            "chain_b",
            "measure",
            "contacts",
            "ca_pairs",
            "residues_a",
            "residues_b",
            "element_pairs",
            *INTERACTION_KINDS,
        ),
        _t_interfaces,
    ),
    "interface_elements": (
        "interfaces",
        ("chain_a", "chain_b", "element_a", "element_b", "ca_pairs", *INTERACTION_KINDS),
        _t_interface_elements,
    ),
    "nucleotides": (
        "nucleic_acid",
        ("index", "strand", *_RES, "base", "partner", "bracket", "level", "stem"),
        _t_nucleotides,
    ),
    "base_pairs": (
        "nucleic_acid",
        ("i", "j", *(f"a_{k}" for k in _RES), *(f"b_{k}" for k in _RES), "kind", "stem", "level"),
        _t_base_pairs,
    ),
    "interactions": (
        "interactions",
        (
            "kind",
            *(f"a_{k}" for k in (*_RES, "atom")),
            *(f"b_{k}" for k in (*_RES, "atom")),
            "distance",
            "angle",
            "note",
            "a_element",
            "b_element",
        ),
        _t_interactions,
    ),
}


def tables(data: dict, names=None) -> dict[str, tuple[tuple[str, ...], list[list]]]:
    """The tables of the sections `data` holds (or those of `names` among them), in TABLES order: name -> (columns,
    rows), rows in the order of the lists they come from."""
    out = {}
    for name, (section, columns, rows) in TABLES.items():
        if section in data and (names is None or name in names):
            out[name] = (columns, rows(data))
    return out


_FORMULA = ("=", "+", "-", "@", "\t", "\r")  # a spreadsheet reads a cell starting so as a formula
_NUMERIC = re.compile(r"[-+]?\d+(\.\d+)?([eE][-+]?\d+)?[A-Za-z]?")  # a number, or a residue number (-7, -7A)


def _cell(v) -> str:
    """A value as a CSV cell: null empty, booleans true/false (as JSON), numbers as written in the JSON. Text read
    from the structure (chain ids, residue and ligand names: a crafted mmCIF may hold anything) that a spreadsheet
    would run as a formula (=, +, -, @ first) is written with a ' before it, unless it is a number or a residue
    number (-7, -7A); the JSON keeps every value as it is."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    text = str(v)
    if isinstance(v, str) and text.startswith(_FORMULA) and not _NUMERIC.fullmatch(text):
        return "'" + text
    return text


def csv_text(columns, rows) -> str:
    """One table as CSV: a header row, then the rows; LF line ends, quoted only where needed."""
    buf = io.StringIO()
    out = csv.writer(buf, lineterminator="\n")
    out.writerow(columns)
    out.writerows([_cell(v) for v in row] for row in rows)
    return buf.getvalue()


def write_csv(columns, rows, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(csv_text(columns, rows))
    return path


def write_tables(data: dict, directory, names=None) -> list[Path]:
    """Every table of `data` (or those `names` names) as DIRECTORY/<name>.csv, in TABLES order."""
    directory = Path(directory)
    return [
        write_csv(columns, rows, directory / f"{name}.csv") for name, (columns, rows) in tables(data, names).items()
    ]


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
    spans = "its CA run within 2 Å of both faces of the drawn slab"
    crosses = (
        "its CA run crosses the slab's centre plane and spans at least half the slab's thickness (the rule the "
        "membrane view and A18 count a transmembrane element by)"
    )
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
    element_or_null = _t("string", "an element id; null: a loop (or a nucleotide)", nullable=True)
    contact = _obj({"residue": RESIDUE, "element": element_or_null})
    ligand = _obj(
        {
            "id": _t("string", "KIND:CHAIN:NUMBER, where it sits"),
            "kind": _t(
                "string",
                "ligand (a small molecule, its parts joined), metal (an ion), glycan (sugars on Asn, Ser or Thr), "
                "linked (a non-standard residue fused into the chain, such as a chromophore)",
                enum=["ligand", "metal", "glycan", "linked"],
            ),
            "name": _t("string", "its residue code (codes joined by + when several)"),
            "codes": _arr(_t("string"), "the residue codes of its parts, in bond order"),
            "symbol": _t("string", "a metal's element; else its short name (sugars by their SNFG names)"),
            "at": RESIDUE,
            "residues": _arr(RESIDUE, "its own residues"),
            "additive": _t("boolean", "a buffer, cryoprotectant, detergent or counter-ion the figure hides unasked"),
            "attached_to": _arr(RESIDUE, "protein residues covalently joined to it: a heavy atom within 1.9 Å"),
            "contacts": _arr(
                contact,
                "the protein residues it touches, in backbone order: a heavy atom within 4.0 Å of one of its heavy "
                "atoms (a metal: an O, N or S within 3.1 Å); its own residue left out",
            ),
            "elements": _arr(_t("string"), "the elements of its contacts, in sequence order"),
            "nucleotides": _arr(RESIDUE, "for one touching no protein residue, the nucleotides holding it"),
        }
    )
    kinds = _obj(
        {k: _t("integer") for k in INTERACTION_KINDS},
        nullable=True,
    )
    kinds["description"] = "interactions by kind (the interactions section); null when it is not asked for"
    interface = _obj(
        {
            "chains": _arr(_t("string"), "the two chains, in the order read"),
            "measure": _t(
                "string",
                "how contact is measured: ca, two protein chains (CA pairs within 8 Å); atoms, a protein and a "
                "nucleic-acid chain (residue-nucleotide pairs with heavy atoms within 4 Å, on the structure); "
                "base_pairs, two nucleic-acid chains (the base pairs the figure draws)",
                enum=["ca", "atoms", "base_pairs"],
            ),
            "contacts": _t("integer", "the pairs in contact by its measure (0: the chains meet only by interactions)"),
            "ca_pairs": _t(
                "integer", "CA pairs within 8 Å, one CA of each chain; null unless both are protein", nullable=True
            ),
            "residues": _arr(_t("integer"), "of each chain, the residues (nucleotides) in those pairs"),
            "elements": _arr(
                _obj(
                    {
                        "a": element_or_null,
                        "b": element_or_null,
                        "ca_pairs": _t(
                            "integer",
                            "CA pairs within 8 Å between the two elements (features.sse_contacts, the contacts the "
                            "layout weighs); null under its 3 pairs",
                            nullable=True,
                        ),
                        "interactions": kinds,
                    }
                ),
                "element pairs in contact, a in the first chain, in sequence order (null ends last)",
            ),
            "interactions": kinds,
        }
    )
    nucleic = _obj(
        {
            "nucleotides": _t("integer", "nucleotides drawn, numbered 1.. through the strands in order"),
            "cropped": _t("boolean", "na_extent cropped a large nucleic acid to the stems the protein touches"),
            "strands": _arr(
                _obj(
                    {
                        "id": _t("string", "CHAIN:FIRST-LAST"),
                        "chain": _t("string"),
                        "rna": _t("boolean"),
                        "first": _t("integer", "the number of its first nucleotide"),
                        "residues": _arr(RESIDUE, "its nucleotides, 5'->3'"),
                        "sequence": _t(
                            "string", "A C G U T, I (inosine), N (no base read); a modified one as its parent"
                        ),
                        "dot_bracket": _t("string", "its part of the dot-bracket"),
                        "cut5": _t("boolean", "the molecule goes on beyond its 5' end (a break, or cropped)"),
                        "cut3": _t("boolean", "the molecule goes on beyond its 3' end"),
                    }
                ),
                "strand pieces as drawn, in the order read",
            ),
            "pairs": _arr(
                _obj(
                    {
                        "i": _t("integer"),
                        "j": _t("integer", "i < j"),
                        "a": RESIDUE,
                        "b": RESIDUE,
                        "kind": _t("string", "WC (Watson-Crick) or GU (G·U wobble)", enum=["WC", "GU"]),
                        "stem": _t("string", "the stem's name (P1, PK, DNA 1); null: an isolated pair", nullable=True),
                        "level": _t("integer", "bracket level: 0 (), 1 [], 2 {}, 3 <>, then Aa, Bb ..."),
                    }
                ),
                "the base pairs the figure draws (nucleic.py: Watson-Crick and G·U wobble), in order of i",
            ),
            "levels": _t("integer", "bracket levels used"),
            "stems": _arr(
                _obj(
                    {
                        "id": _t("string"),
                        "name": _t("string"),
                        "pairs": _t("integer"),
                        "level": _t("integer"),
                        "pseudoknot": _t("boolean", "foldmap's pseudoknot rule calls it a knot (named PK)"),
                    }
                ),
                "stems, 5'->3'",
            ),
            "gaps": _arr(
                _obj({"before": _t("string"), "after": _t("string"), "nucleotides": _t("integer")}),
                "nucleotides left out between two pieces of one molecule (a break, or cropping)",
            ),
            "tetrads": _arr(_arr(RESIDUE), "G-tetrads (read, not drawn)"),
        },
        nullable=True,
    )
    end = _obj({**RESIDUE["properties"], "atom": _t("string", "an atom name, ring:5 / ring:6, or centre (Arg)")})
    interactions = _obj(
        {
            "between": _t("string", "the partners: interchain, residues of two different chains"),
            "kinds": _arr(_t("string")),
            "found": _arr(
                _obj(
                    {
                        "kind": _t("string", enum=list(INTERACTION_KINDS)),
                        "a": end,
                        "b": end,
                        "distance": _t("number", "Å"),
                        "angle": _t("number", "degrees (stacking, cation-pi)", nullable=True),
                        "note": _t("string", "parallel / T-shaped (stacking), His if protonated (salt bridges)"),
                        "elements": _arr(element_or_null, "the elements of a and b"),
                    }
                ),
                "each interaction (an H-bond: a donor, b acceptor; a salt bridge or cation-pi: a the cation), in "
                "order of a, then b, then kind",
            ),
            "failed_recheck": _arr(
                _t("string"), "interactions found but flagged by the independent re-check, left out of found"
            ),
        }
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
            "ligands": _arr(ligand, "ligands and metal ions (as io found them), then glycans, then linked residues"),
            "interfaces": _arr(interface, "chain pairs in contact, in the order the chains were read"),
            "nucleic_acid": nucleic,
            "interactions": interactions,
        },
        "required": ["foldmap_data", "foldmap"],
        "additionalProperties": False,
    }
