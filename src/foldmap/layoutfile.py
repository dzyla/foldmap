# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Layout files: one YAML document holding everything that defines a figure, so a hand-tuned figure can be
re-rendered later, after the model changes, or shared.

    foldmap: 1                       # file format version
    structure: model.cif              # for the record; the structure is always given on the command line
    theme: trace
    style: {residue_numbers: true}    # only what differs from the theme
    layout: {mode: projected, rotate: 0, symmetry: auto, assembly: auto,
             like: {structure: 2HHB, chain: B, match: auto},   # --like, the reference as given
             data: hdx.csv, data_column: uptake, crosslinks: xl.csv}   # the researcher's files (userdata), as given
    names: [deposited, my-names.yaml]  # --names, as given
    edits:
      rename: {res:A:167-182: Gd}
      swap: [[res:A:3-11, res:A:60-79]]
      move: {res:A:39-50: [1.0, -2.0]}

Element references are saved as residue ranges (res:chain:first-last): labels such as G or α3 can shift when
a model gains or loses an element, residue ranges keep pointing at the same stretch of chain."""

from __future__ import annotations

import math
import reprlib
from dataclasses import asdict
from pathlib import Path

import yaml

from .layout import Layout, panel_members, resolve
from .model import SSE, res_ref
from .style import _KEY_ALIASES, THEMES, Style
from .yamlfile import read_yaml

VERSION = 1
DOMAIN_SOURCES = ("inferred", "CATH", "UniProt")  # Layout.domain_source, as a layout file keeps it
TOP = ("foldmap", "structure", "theme", "style", "layout", "domains", "domain_source", "names", "edits")
LAYOUT = (
    "mode",
    "rotate",
    "flip_v",
    "dna",
    "symmetry",
    "symmetry_tol",
    "assembly",
    "protomers",
    "up",
    "view",
    "title",
    "msa",
    "msa_reference",
    "membrane",
    "chains",
    "focus",
    "uniprot",
    "like",
    "data",
    "data_column",
    "crosslinks",
)
EDITS = ("rename", "swap", "move")


def canonical(layout: Layout, sses: list[SSE], ref: str) -> list[str]:
    by_id = {s.id: s for s in sses}
    return [_range(by_id[k]) for k in resolve(layout, sses, ref)]


def document(
    structure, theme: str, look: Style, options: dict, edits: dict, layout: Layout, sses: list[SSE], names=()
) -> dict:
    """The layout document for a figure: theme, style changes, layout options, names (--names, as given), edits in
    residue form."""
    base = THEMES[theme]
    style = {k: v for k, v in asdict(look).items() if v != getattr(base, k)}
    if style.get("mark"):
        marks = []
        for entry in filter(None, (e.strip() for e in style["mark"].split(","))):
            ref, sep, colour = entry.partition("=")
            marks += [f"{c}{sep}{colour}" for c in canonical(layout, sses, ref)]
        style["mark"] = ",".join(marks)
    rename = {c: name for ref, name in (edits.get("rename") or {}).items() for c in canonical(layout, sses, ref)}
    swap = [[canonical(layout, sses, a)[0], canonical(layout, sses, b)[0]] for a, b in edits.get("swap") or []]
    move = {}
    for ref, (dx, dy) in edits.get("move") or []:
        for c in canonical(layout, sses, ref):
            move[c] = [float(dx), float(dy)]
    doc = {
        "foldmap": VERSION,
        "structure": str(structure),
        "theme": theme,
        "style": style,
        "layout": {k: v for k, v in options.items() if k in LAYOUT and v is not None},
    }
    if layout.domains:  # saved as resolved element ranges, so 'auto' domains are frozen as found
        by_id = {s.id: s for s in sses}
        doc["domains"] = {  # with the elements neighbours lend (panel_members): read back, they are lent again
            name: [_range(by_id[k]) for k in panel_members(layout, name)] for name, _ in layout.domains
        }
        if layout.domain_source:  # inferred, CATH or UniProt: said again when the file is read back
            doc["domain_source"] = layout.domain_source
    if names:
        doc["names"] = [str(n) for n in names]
    doc["edits"] = {"rename": rename, "swap": swap, "move": move}
    return doc


def _range(s: SSE) -> str:
    return res_ref(s.chain, s.first, s.last)


def save(path, doc: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    head = "# foldmap layout: re-render with  foldmap plot STRUCTURE --layout-file " + path.name + "\n"
    path.write_text(head + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return path


def load(path) -> dict:
    """The document, checked: unknown keys and other versions are errors that say what is wrong, and so is a
    section of the wrong kind (style, layout and edits are mappings, names a list of text). The values are the
    reader's to check (spec.figure_spec, app.settings_from_document)."""
    doc = read_yaml(Path(path).read_text(encoding="utf-8"), "layout", str(path))
    if not isinstance(doc, dict):
        raise ValueError(f"{path}: a layout file is a mapping of keys (foldmap, theme, style, layout, edits)")
    if doc.get("foldmap", VERSION) != VERSION:
        raise ValueError(
            f"{path}: layout file version {_brief(doc.get('foldmap'))} is not supported (this is version {VERSION})"
        )
    for key in doc:
        if key not in TOP:
            raise ValueError(f"{path}: unknown key {_brief(key)}; a layout file has {', '.join(TOP)}")
    theme = doc.get("theme")
    if theme is not None and not (isinstance(theme, str) and theme in THEMES):
        raise ValueError(f"{path}: unknown theme {_brief(theme)}")
    for section, what in (("style", "style key: value"), ("layout", "option: value"), ("edits", "rename, swap, move")):
        if not isinstance(doc.get(section) or {}, dict):
            raise ValueError(f"{path}: {section} is a mapping ({what}); got {_brief(doc[section])}")
    for key in doc.get("style") or {}:
        if not (isinstance(key, str) and hasattr(Style(), _KEY_ALIASES.get(key, key))):  # earlier names are read too
            raise ValueError(f"{path}: unknown style key {_brief(key)}")
    for key in doc.get("layout") or {}:
        if key not in LAYOUT:
            raise ValueError(f"{path}: unknown layout key {_brief(key)}; use {', '.join(LAYOUT)}")
    source = doc.get("domain_source")
    if source is not None and source not in DOMAIN_SOURCES:
        raise ValueError(f"{path}: domain_source is one of {', '.join(DOMAIN_SOURCES)}; got {_brief(source)}")
    names = doc.get("names")
    if names is not None and not (isinstance(names, list) and all(isinstance(n, str) for n in names)):
        raise ValueError(f"{path}: names is a list, e.g. [deposited] or [beta, my-names.yaml]")
    for key in doc.get("edits") or {}:
        if key not in EDITS:
            raise ValueError(f"{path}: unknown edit {_brief(key)}; use {', '.join(EDITS)}")
    return doc


NAME_MAX = 80  # characters of an element's name (edits: rename)
REF_MAX = 200  # characters of an element reference (res:A:3-18, α3, A:G)
MOVE_MAX = 100.0  # page units an edit may move an element by, each way (one strand spacing is about 1.1)


def _brief(value) -> str:
    """A value as a message shows it: its repr, shortened (and never deeper than a few levels)."""
    return reprlib.repr(value)


def _word(edit: str, value, limit: int) -> str:
    """A reference or a name in the edits: text (YAML may read 1 or 1.5 as a number), at most `limit` characters."""
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ValueError(f"edits: {edit} takes element references and names as text; got {_brief(value)}")
    text = str(value)
    if len(text) > limit:
        raise ValueError(f"edits: {edit} has a value of {len(text)} characters; at most {limit}")
    return text


def _offset(ref: str, d) -> float:
    if isinstance(d, bool) or not isinstance(d, (int, float)) or abs(d) > MOVE_MAX or not math.isfinite(d):
        raise ValueError(
            f"edits: move of {ref} must be two numbers from -{MOVE_MAX:g} to {MOVE_MAX:g}; got {_brief(d)}"
        )
    return float(d)


def edits_of(doc: dict) -> dict:
    """The document's edits in the form the layout takes (rename dict, swap pairs, move list), checked: a value of
    the wrong kind is a ValueError naming the edit."""
    e = doc.get("edits") or {}
    if not isinstance(e, dict):
        raise ValueError(f"edits is a mapping (rename, swap, move); got {_brief(e)}")
    rename, swap, move = e.get("rename") or {}, e.get("swap") or [], e.get("move") or {}
    if not isinstance(rename, dict):
        raise ValueError(f"edits: rename is a mapping of element to name, e.g. res:A:167-182: Gd; got {_brief(rename)}")
    if not (isinstance(swap, list) and all(isinstance(p, list) and len(p) == 2 for p in swap)):
        raise ValueError(
            f"edits: swap is a list of element pairs, e.g. [[res:A:3-11, res:A:60-79]]; got {_brief(swap)}"
        )
    if not (isinstance(move, dict) and all(isinstance(d, list) and len(d) == 2 for d in move.values())):
        raise ValueError(
            f"edits: move is a mapping of element to [dx, dy], e.g. res:A:39-50: [1.0, -2.0]; got {_brief(move)}"
        )
    moves = []
    for ref, d in move.items():
        ref = _word("move", ref, REF_MAX)
        moves.append((ref, (_offset(ref, d[0]), _offset(ref, d[1]))))
    return {
        "rename": {_word("rename", k, REF_MAX): _word("rename", v, NAME_MAX) for k, v in rename.items()},
        "swap": [(_word("swap", a, REF_MAX), _word("swap", b, REF_MAX)) for a, b in swap],
        "move": moves,
    }
