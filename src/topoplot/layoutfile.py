"""Layout files: one YAML document holding everything that defines a figure, so a hand-tuned figure can be
re-rendered later, after the model changes, or shared.

    topoplot: 1                       # file format version
    structure: model.cif              # for the record; the structure is always given on the command line
    theme: trace
    style: {residue_numbers: true}    # only what differs from the theme
    layout: {mode: projected, rotate: 0, symmetry: auto, assembly: auto}
    edits:
      rename: {res:A:167-182: Gd}
      swap: [[res:A:3-11, res:A:60-79]]
      move: {res:A:39-50: [1.0, -2.0]}

Element references are saved as residue ranges (res:chain:first-last): labels such as G or α3 can shift when
a model gains or loses an element, residue ranges keep pointing at the same stretch of chain."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import yaml

from .layout import Layout, resolve
from .model import SSE
from .style import THEMES, Style

VERSION = 1
TOP = ("topoplot", "structure", "theme", "style", "layout", "domains", "edits")
LAYOUT = ("mode", "rotate", "flip_v", "symmetry", "symmetry_tol", "assembly", "protomers", "up", "view", "title")
EDITS = ("rename", "swap", "move")


def canonical(layout: Layout, sses: list[SSE], ref: str) -> list[str]:
    by_id = {s.id: s for s in sses}
    return [_range(by_id[k]) for k in resolve(layout, sses, ref)]


def document(structure, theme: str, look: Style, options: dict, edits: dict, layout: Layout,
             sses: list[SSE]) -> dict:
    """The layout document for a figure: theme, style changes, layout options, edits in residue form."""
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
    doc = {"topoplot": VERSION, "structure": str(structure), "theme": theme, "style": style,
           "layout": {k: v for k, v in options.items() if k in LAYOUT and v is not None}}
    if layout.domains:  # saved as resolved element ranges, so 'auto' domains are frozen as found
        by_id = {s.id: s for s in sses}
        doc["domains"] = {name: [_range(by_id[k]) for k in ids] for name, ids in layout.domains}
    doc["edits"] = {"rename": rename, "swap": swap, "move": move}
    return doc


def _range(s: SSE) -> str:
    return f"res:{s.chain}:{s.first.seq}-{s.last.seq}"


def save(path, doc: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    head = "# topoplot layout: re-render with  topoplot plot STRUCTURE --layout-file " + path.name + "\n"
    path.write_text(head + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True))
    return path


def load(path) -> dict:
    """The document, checked: unknown keys and other versions are errors that say what is wrong."""
    doc = yaml.safe_load(Path(path).read_text())
    if not isinstance(doc, dict):
        raise ValueError(f"{path}: a layout file is a mapping of keys (topoplot, theme, style, layout, edits)")
    if doc.get("topoplot", VERSION) != VERSION:
        raise ValueError(f"{path}: layout file version {doc.get('topoplot')} is not supported (this is version {VERSION})")
    for key in doc:
        if key not in TOP:
            raise ValueError(f"{path}: unknown key {key!r}; a layout file has {', '.join(TOP)}")
    if doc.get("theme") is not None and doc["theme"] not in THEMES:
        raise ValueError(f"{path}: unknown theme {doc['theme']!r}")
    for key in doc.get("style") or {}:
        if not hasattr(Style(), key):
            raise ValueError(f"{path}: unknown style key {key!r}")
    for key in doc.get("layout") or {}:
        if key not in LAYOUT:
            raise ValueError(f"{path}: unknown layout key {key!r}; use {', '.join(LAYOUT)}")
    for key in doc.get("edits") or {}:
        if key not in EDITS:
            raise ValueError(f"{path}: unknown edit {key!r}; use {', '.join(EDITS)}")
    return doc


def edits_of(doc: dict) -> dict:
    """The document's edits in the form the layout takes (rename dict, swap pairs, move list)."""
    e = doc.get("edits") or {}
    return {"rename": dict(e.get("rename") or {}),
            "swap": [tuple(p) for p in e.get("swap") or []],
            "move": [(ref, (float(d[0]), float(d[1]))) for ref, d in (e.get("move") or {}).items()]}
