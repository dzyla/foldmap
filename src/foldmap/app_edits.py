# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""The app's structured controls, translated to and from the text foldmap's command line takes.

Every figure setting keeps one canonical string in the session (what `app.spec_of`, layout files and the cache key
read). Pickers, segmented choices and tables live on derived keys (`<key>__mode`, `<key>__pick`, ...) and write the
canonical string; `sync_composites` derives them back when the canonical value changes elsewhere (a layout file, a
theme reset). Pure: no Streamlit, no pandas; any mutable mapping stands in for the session state."""

from __future__ import annotations

import math
import re
import unicodedata
from collections.abc import Callable, MutableMapping

from foldmap.fetch import ALIGNMENT_SUFFIXES, STRUCTURE_SUFFIXES, upload_suffix

# ------------------------------------------------------------------------------------------ composites
NA_BACKBONE_COLOUR = "#8c564b"  # the picker's colour until one is chosen


def _list(text: str) -> list[str]:
    return [c.strip() for c in (text or "").split(",") if c.strip()]


def _chains_parts(text: str) -> tuple:
    return (_list(text),)


def _chains_text(pick) -> str:
    return ",".join(pick or [])


def _highlight_parts(text: str) -> tuple:
    text = (text or "").strip()
    return (text or "none", []) if text in ("", "none", "asu", "protomer") else ("chains", _list(text))


def _highlight_text(mode: str, chains) -> str:
    return (",".join(chains or []) or "none") if mode == "chains" else mode


def _focus_parts(text: str) -> tuple:
    text = (text or "").strip()
    return (text or "auto", []) if text in ("", "auto", "none") else ("chains", _list(text))


def _focus_text(mode: str, chains) -> str:
    return (",".join(chains or []) or "auto") if mode == "chains" else mode


def _ligands_parts(text: str) -> tuple:
    text = (text or "").strip()
    return (text or "auto", "") if text in ("", "auto", "all", "none") else ("codes", ", ".join(_list(text)))


def _ligands_text(mode: str, codes: str) -> str:
    if mode != "codes":
        return mode
    return ",".join(c.upper() for c in re.split(r"[,\s]+", codes or "") if c) or "auto"


def hex_colour(value, default: str) -> str:
    """A colour as st.color_picker takes it (#rrggbb): a named colour (teal) or #rgb converted, anything that is
    not a colour gives `default`."""
    from matplotlib.colors import is_color_like, to_hex

    value = str(value or "").strip()
    return to_hex(value) if value and is_color_like(value) else default


def _na_backbone_parts(text: str) -> tuple:
    text = (text or "").strip()
    if text in ("", "modality", "chain", "ink"):
        return (text or "modality", NA_BACKBONE_COLOUR)
    return ("colour", hex_colour(text, NA_BACKBONE_COLOUR))


def _na_backbone_text(mode: str, colour: str) -> str:
    return (colour or NA_BACKBONE_COLOUR) if mode == "colour" else mode


def _hbonds_parts(text: str) -> tuple:
    text = (text or "").strip()
    if text.startswith("chains:"):
        return ("chains", ", ".join(_list(text[len("chains:") :])))
    return (text or "off", "")


def _hbonds_text(mode: str, pairs: str) -> str:
    if mode != "chains":
        return mode
    pairs = [p.replace(" ", "") for p in _list(pairs)]
    return "chains:" + ",".join(pairs) if pairs else "off"


def _uniprot_parts(text: str) -> tuple:
    text = (text or "").strip()
    return (text or "off", "") if text in ("", "off", "auto") else ("accession", text)


def _uniprot_text(mode: str, acc: str) -> str:
    return ((acc or "").strip().upper() or "off") if mode == "accession" else mode


def _sites_parts(text: str) -> tuple:
    text = (text or "").strip()
    return (text or "auto", []) if text in ("", "auto", "all", "none") else ("types", _list(text))


def _sites_text(mode: str, types) -> str:
    from foldmap.render import SITE_GROUPS

    if mode != "types":
        return mode
    return ",".join(t for t in SITE_GROUPS if t in set(types or [])) or "auto"


COMPOSITES: dict[str, tuple[tuple[str, ...], Callable, Callable]] = {
    "chains": (("chains__pick",), _chains_parts, _chains_text),
    "highlight": (("highlight__mode", "highlight__chains"), _highlight_parts, _highlight_text),
    "focus": (("focus__mode", "focus__chains"), _focus_parts, _focus_text),
    "ligands": (("ligands__mode", "ligands__codes"), _ligands_parts, _ligands_text),
    "na_backbone": (("na_backbone__mode", "na_backbone__colour"), _na_backbone_parts, _na_backbone_text),
    "hbonds": (("hbonds__mode", "hbonds__pairs"), _hbonds_parts, _hbonds_text),
    "uniprot": (("uniprot__mode", "uniprot__acc"), _uniprot_parts, _uniprot_text),
    "uniprot_sites": (("uniprot_sites__mode", "uniprot_sites__types"), _sites_parts, _sites_text),
}


NAME_KEYWORDS = ("deposited", "beta")  # names.KEYWORDS (the pipeline stays out of this module)


def _names_parts(names) -> tuple:
    keywords = [n for n in names or [] if n in NAME_KEYWORDS]
    return (keywords[-1] if keywords else "generic",)  # a later spec wins


def _domains_parts(text) -> tuple:
    return (str(text or "").strip() == "auto",)


DERIVED: dict[str, tuple[tuple[str, ...], Callable]] = {  # widget keys derived like a composite's, set by their
    "names": (("names__kw",), _names_parts),  # own callbacks (choose_names, switch_domains)
    "domains": (("domains__auto",), _domains_parts),
}


def _noted(ss: MutableMapping, key: str) -> None:
    ss["_synced"] = {**(ss.get("_synced") or {}), key: ss[key]}
    ss["last_changed"] = key


def choose_names(ss: MutableMapping) -> None:
    """The Element names choice (`names__kw`: generic, deposited or beta): that keyword (none for generic), then the
    names files, which win over it."""
    files = [n for n in ss.get("names") or [] if n not in NAME_KEYWORDS]
    kw = ss.get("names__kw")
    ss["names"] = ([kw] if kw in NAME_KEYWORDS else []) + files
    _noted(ss, "names")


def switch_domains(ss: MutableMapping) -> None:
    """The Automatic domains switch (`domains__auto`): on, `auto` (the domains typed so far kept in
    `domains__manual`); off, those domains back."""
    if ss.get("domains__auto"):
        if not _domains_parts(ss.get("domains"))[0]:
            ss["domains__manual"] = ss.get("domains") or ""
        ss["domains"] = "auto"
    else:
        ss["domains"] = ss.get("domains__manual", "")
    _noted(ss, "domains")


def sync_composites(ss: MutableMapping) -> None:
    """Derive each composite's widget keys (and DERIVED's) from its canonical value, but only when that value changed
    since it was last derived or composed: a half-made choice (a mode picked, nothing typed yet) is never snapped
    back."""
    synced = dict(ss.get("_synced") or {})
    every = {**{key: (derived, parts) for key, (derived, parts, _) in COMPOSITES.items()}, **DERIVED}
    for key, (derived, parts) in every.items():
        if key not in ss or (key in synced and synced[key] == ss[key]):
            continue
        for k, v in zip(derived, parts(ss[key])):
            ss[k] = v
        synced[key] = ss[key]
    ss["_synced"] = synced


def compose_into(ss: MutableMapping, key: str) -> None:
    """Write a composite's canonical value from its widget keys (the widgets' callback)."""
    derived, parts, text = COMPOSITES[key]
    current = parts(ss.get(key, ""))
    ss[key] = text(*(ss.get(k, d) for k, d in zip(derived, current)))
    ss["_synced"] = {**(ss.get("_synced") or {}), key: ss[key]}
    ss["last_changed"] = key


# --------------------------------------------------------------------------------------------- editors
EDITORS = {
    "swap": ("a", "b"),
    "move": ("element", "dx", "dy"),
    "rename": ("element", "name"),
    "mark": ("element", "colour"),
    "domains": ("name", "ranges"),
}
NEEDED = {"swap": ("a", "b"), "rename": ("element", "name"), "domains": ("name", "ranges")}  # a table row is
# written once these are filled (else its first column): a swap of one element, a rename with no name or a domain
# with no ranges would be an error
_RANGE_IN_CHAIN = re.compile(r"(?![A-Za-z]+\d+:)[A-Za-z0-9]+:-?\d+[A-Z]?(--?\d+[A-Z]?)?")  # not a copy's (A1:)
_RANGE = re.compile(r"-?\d+[A-Z]?(--?\d+[A-Z]?)?")


def _cell(v) -> str:
    """A table cell as text: an empty cell (None, NaN) is ''."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return str(v).strip()


def _number(v) -> str:
    """A move offset as written: numbers with :g, an empty cell as 0, anything else as typed."""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int | float):
        return "0" if isinstance(v, float) and math.isnan(v) else f"{v:g}"
    return _cell(v) or "0"


def _float_or_text(text: str):
    try:
        return float(text)
    except ValueError:
        return text.strip()


def _residues(ref: str) -> str:
    """A domain range cell item as a reference: A:1-80 -> res:A:1-80, 1-80 -> res:1-80. Anything else stays as
    written: labels, residue references, and the element ids of an assembly's copies (A1:3-7: a chain id with a
    copy number), which stay element ids rather than become a guess at a residue range."""
    ref = ref.strip()
    return f"res:{ref}" if _RANGE_IN_CHAIN.fullmatch(ref) or _RANGE.fullmatch(ref) else ref


def rows_of(kind: str, text: str) -> list[dict]:
    """Today's text of an editor as table rows (columns `EDITORS[kind]`)."""
    text = text or ""
    if kind == "mark":
        items = _list(text)
    else:
        items = [p.strip() for p in re.split(r"[;\n]", text) if p.strip()]
    rows = []
    for item in items:
        if kind == "swap":
            a, _, b = item.partition(",")
            rows.append({"a": a.strip(), "b": b.strip()})
        elif kind == "move":
            ref, _, offset = item.partition("=")
            dx, _, dy = offset.partition(",")
            rows.append({"element": ref.strip(), "dx": _float_or_text(dx), "dy": _float_or_text(dy)})
        elif kind == "rename":
            ref, _, name = item.partition("=")
            rows.append({"element": ref.strip(), "name": name.strip()})
        elif kind == "mark":
            ref, _, colour = item.partition("=")
            rows.append({"element": ref.strip(), "colour": colour.strip()})
        elif kind == "domains":
            if item == "auto":  # automatic domains are a switch, not a row
                continue
            name, _, refs = item.partition("=")
            rows.append({"name": name.strip(), "ranges": ", ".join(_list(refs))})
        else:
            raise KeyError(kind)
    return rows


def text_of(kind: str, rows: list[dict]) -> str:
    """Table rows as the text the command line takes; rows whose first column is empty are dropped."""
    first = EDITORS[kind][0]
    items = []
    for row in rows:
        if not _cell(row.get(first)):
            continue
        if kind == "swap":
            a, b = _cell(row.get("a")), _cell(row.get("b"))
            items.append(f"{a},{b}" if b else a)
        elif kind == "move":
            items.append(f"{_cell(row.get('element'))}={_number(row.get('dx'))},{_number(row.get('dy'))}")
        elif kind == "rename":
            items.append(f"{_cell(row.get('element'))}={_cell(row.get('name'))}")
        elif kind == "mark":
            ref, colour = _cell(row.get("element")), _cell(row.get("colour"))
            items.append(f"{ref}={colour}" if colour else ref)
        elif kind == "domains":
            refs = ",".join(_residues(r) for r in _list(_cell(row.get("ranges"))))
            items.append(f"{_cell(row.get('name'))}={refs}")
    return {"mark": ", ", "domains": "\n"}.get(kind, "; ").join(items)


def table_rows(ss: MutableMapping, kind: str) -> list[dict]:
    """The rows an editor's table shows: its text's rows, then the rows not filled in yet (NEEDED cells empty: a
    row just added, a swap with one element picked), which are not written to the text: they wait in
    `<kind>__blank`."""
    return rows_of(kind, ss.get(kind, "")) + [dict(row) for row in ss.get(f"{kind}__blank") or []]


def apply_table(ss: MutableMapping, kind: str) -> None:
    """Apply a data editor's edits (`{"edited_rows", "added_rows", "deleted_rows"}`, indices into the rows it was
    given: table_rows) to the editor's text; rows not filled in yet stay in `<kind>__blank` (the add-row button
    sends an empty row at once, and a half-made swap would be an error). The version moves on, so the next run
    builds a new editor key and the same delta is never applied twice."""
    ver = ss.get(f"{kind}__ver", 0)
    delta = ss.get(f"{kind}__table_{ver}")
    if delta is None:
        return
    rows = table_rows(ss, kind)
    for i, change in (delta.get("edited_rows") or {}).items():
        if 0 <= int(i) < len(rows):
            rows[int(i)].update(change)
    gone = {int(i) for i in delta.get("deleted_rows") or []}
    rows = [r for i, r in enumerate(rows) if i not in gone]
    rows += [{c: added.get(c, "") for c in EDITORS[kind]} for added in delta.get("added_rows") or []]
    needed = NEEDED.get(kind, EDITORS[kind][:1])
    filled = [all(_cell(row.get(c)) for c in needed) for row in rows]
    ss[f"{kind}__blank"] = [row for row, done in zip(rows, filled) if not done]
    ss[kind] = text_of(kind, [row for row, done in zip(rows, filled) if done])
    ss[f"{kind}__ver"] = ver + 1
    ss["last_changed"] = kind


def element_choices(elements: list[dict]) -> list[tuple[str, str]]:
    """(reference, what the picker shows) for each drawn element: ('res:A:23-34', 'α1 · A 23–34')."""
    return [(e["ref"], f"{e['label']} · {e['chain']} {e['first']}–{e['last']}") for e in elements]


def display(ref: str, choices: list[tuple[str, str]]) -> str:
    """What a picker shows for a reference: the element's name, or the reference itself when it is not drawn."""
    return dict(choices).get(ref, ref)


# ------------------------------------------------------------------------------------------- choices
def chain_choices(deposited: list[str], drawn: list[str]) -> list[str]:
    """The chains the Chains picker offers: the file's own, in its order (a chain the first figure left out can
    still be asked for), unless that figure drew an assembly's copies, whose names the file does not hold (A1, A2):
    then the chains drawn."""
    return list(deposited) if set(drawn) <= set(deposited) else list(dict.fromkeys(drawn))


_UNPALETTED = ("plddt", "bfactor", "hydropathy", "conservation", "sequence")  # colourings with their own scale


def palette_used(color_by: str, loop_color: str) -> bool:
    """Whether the palette colours anything: not under a colouring with its own scale, unless loops take their
    chain's colour."""
    return color_by not in _UNPALETTED or loop_color == "chain"


SYMMETRIES = (
    "auto",
    "off",
    "helical",
    *(f"C{n}" for n in range(2, 13)),
    *(f"D{n}" for n in range(2, 9)),
)

# ------------------------------------------------------------------------------------ controls, errors
_TOP, _STYLE, _SHOW = "top of the sidebar", "Style", "Show"
_CHAINS, _LAYOUT, _CONSERVATION = "Chains & assembly", "Edit layout", "Conservation & UniProt"
HIDDEN = ("flip_v", "dna", "symmetry_tol", "protomers", "up", "view")  # app.HIDDEN_LAYOUT: layout files only

CONTROLS: dict[str, tuple[str, str, tuple[str, ...]]] = {
    # the essentials
    "theme": ("Theme", _TOP, ("theme", "theme_browse")),
    "color_by": ("Colour by", _TOP, ("color_by",)),
    "palette": ("Palette", _TOP, ("palette",)),
    "sequence_map": ("N → C colours", _TOP, ("sequence_map",)),
    "labels": ("Labels", _TOP, ("labels",)),
    "residue_numbers": ("Residue numbers", _TOP, ("residue_numbers",)),
    "legend": ("Legend", _TOP, ("legend",)),
    "font_scale": ("Text size", _TOP, ("font_scale",)),
    "title": ("Title", _TOP, ("title",)),
    # Style
    "fill": ("Fill", _STYLE, ("fill",)),
    "helix_shading": ("Helix shading", _STYLE, ("helix_shading",)),
    "helix_angle": ("Helix angle", _STYLE, ("helix_angle",)),
    "loops": ("Loop shape", _STYLE, ("loops",)),
    "loop_color": ("Loop colour", _STYLE, ("loop_color",)),
    "loop_arrows": ("Loop arrows", _STYLE, ("loop_arrows",)),
    "loop_width": ("Loop weight", _STYLE, ("loop_width",)),
    "background": ("Background", _STYLE, ("background",)),
    "ink": ("Ink", _STYLE, ("ink",)),
    # Show
    "disulfides": ("Disulfides", _SHOW, ("disulfides",)),
    "glycans": ("Glycans", _SHOW, ("glycans",)),
    "sheet_panels": ("Sheet panels", _SHOW, ("sheet_panels",)),
    "helices_310": ("3₁₀ helices as boxes", _SHOW, ("helices_310",)),
    "ligands": ("Ligands and ions", _SHOW, ("ligands__mode", "ligands__codes")),
    "membrane": ("Membrane", _SHOW, ("membrane",)),
    "hbonds": ("Hydrogen bonds", _SHOW, ("hbonds__mode", "hbonds__pairs")),
    "hbond_color": ("H-bond colour", _SHOW, ("hbond_color",)),
    "hbond_max": ("Most H-bonds drawn", _SHOW, ("hbond_max",)),
    "nucleotide_labels": ("Base letters", _SHOW, ("nucleotide_labels",)),
    "na_style": ("Stems", _SHOW, ("na_style",)),
    "na_extent": ("Shown", _SHOW, ("na_extent",)),
    "na_backbone": ("Backbone", _SHOW, ("na_backbone__mode", "na_backbone__colour")),
    "dna_color": ("DNA colour", _SHOW, ("dna_color",)),
    "rna_color": ("RNA colour", _SHOW, ("rna_color",)),
    "backbone_colors": ("Per-chain colours", _SHOW, ("backbone_colors",)),
    # Chains & assembly
    "chains": ("Chains", _CHAINS, ("chains__pick", "chains__as_text", "chains")),
    "highlight": ("Highlight", _CHAINS, ("highlight__mode", "highlight__chains", "highlight__as_text", "highlight")),
    "assembly": ("Assembly", _CHAINS, ("assembly",)),
    "symmetry": ("Symmetry", _CHAINS, ("symmetry",)),
    "focus": ("Large assemblies", _CHAINS, ("focus__mode", "focus__chains")),
    # Edit layout
    "mode": ("Arrangement", _LAYOUT, ("mode",)),
    "rotate": ("Rotate", _LAYOUT, ("rotate",)),
    "helix_scale": ("Helix size", _LAYOUT, ("helix_scale",)),
    "strand_scale": ("Strand size", _LAYOUT, ("strand_scale",)),
    "dna_scale": ("DNA size", _LAYOUT, ("dna_scale",)),
    "fold_template": ("Parallel α/β sheets", _LAYOUT, ("fold_template",)),
    "helix_on_sheet": ("Helices packed on a sheet", _LAYOUT, ("helix_on_sheet",)),
    "names": ("Element names", _LAYOUT, ("names__kw", "names_upload")),
    "swap": ("Swap", _LAYOUT, ("swap__as_text", "swap")),
    "move": ("Move", _LAYOUT, ("move__as_text", "move")),
    "rename": ("Rename", _LAYOUT, ("rename__as_text", "rename")),
    "mark": ("Mark", _LAYOUT, ("mark__as_text", "mark")),
    "domains": ("Domains", _LAYOUT, ("domains__auto", "domains__as_text", "domains")),
    "domains_from_file": ("Domains", _LAYOUT, ("domains__auto", "domains__as_text", "domains")),
    "domain_colours": ("Colour by domain", _LAYOUT, ("domain_colours",)),
    "domain_bar": ("Domain bar", _LAYOUT, ("domain_bar",)),
    # Conservation & UniProt
    "msa_path": ("Alignment", _CONSERVATION, ("msa_upload",)),
    "msa_reference": ("Reference row", _CONSERVATION, ("msa_reference",)),
    "uniprot": ("UniProt", _CONSERVATION, ("uniprot__mode", "uniprot__acc")),
    "uniprot_domains": ("UniProt domains as panels", _CONSERVATION, ("uniprot_domains",)),
    "uniprot_sites": ("UniProt sites", _CONSERVATION, ("uniprot_sites__mode", "uniprot_sites__types")),
    **{key: ("", "", ()) for key in HIDDEN},
}

_NEEDLES = (  # (in the lower-cased error message, the setting it points at); the first match wins
    ("swap", "swap"),
    ("--move", "move"),
    ("rename", "rename"),
    ("mark", "mark"),
    ("domain", "domains"),
    ("highlight", "highlight"),
    ("no chain", "chains"),
    ("assembly", "assembly"),
    ("symmetry", "symmetry"),
    ("focus", "focus"),
    ("ligands", "ligands"),
    ("hbonds", "hbonds"),
    ("backbone_colors", "backbone_colors"),
    ("na_backbone", "na_backbone"),
    ("uniprot", "uniprot"),
    ("needs an alignment", "msa_path"),
    ("--names", "names"),
)


def control_for(message: str, last_changed: str | None) -> str | None:
    """The control to check when a figure fails: 'Swap (Edit layout)', from the message, else the setting changed
    last; None when neither names one."""
    text = message.lower()
    key = next((k for needle, k in _NEEDLES if needle in text), last_changed)
    label, section, _ = CONTROLS.get(key or "", ("", "", ()))
    return f"{label} ({section})" if label else None


def notes_text(fallback: list[str], over_labels: list[str], crowded: list[str] = (), page: list[str] = ()) -> str:
    """The routing and placement notes under a figure, as one line ('' when there are none): loops with no clear
    route, loops over a label, ligand marks with no clear spot (Drawn.crowded_marks), and what is said of its page
    (Drawn.page_notes: held at its largest with text under 6 pt; a PNG under 300 dpi)."""
    notes = [f"{note[0].upper()}{note[1:]}." for note in page]
    if fallback:
        notes.append(f"No clear route for {', '.join(fallback)} (drawn as plain curves).")
    if over_labels:
        notes.append(f"Drawn over a label or chain end, the only route found: {', '.join(over_labels)}.")
    if crowded:
        notes.append(f"{len(crowded)} ligand mark(s) with no clear spot, drawn over other ink: {', '.join(crowded)}.")
    return " ".join(notes)


def row_source(rows: list[int], choices: list[tuple[str, str]]) -> str | None:
    """The id of the structure in the selected table row; None when nothing (or nothing listed) is selected."""
    if not rows or not 0 <= rows[0] < len(choices):
        return None
    return choices[rows[0]][0]


# --------------------------------------------------------------------------------------------- shown text
_PATH = re.compile(r"(?<![\w:/.~-])(?:~?/|[A-Za-z]:\\|\\\\)[^\s'\"`<>|/\\]+(?:[/\\][^\s'\"`<>|/\\]*)*")
_FILE_URL = re.compile(r"\bfile:/+", re.IGNORECASE)
_MARKDOWN = re.compile(r"([\\`*_\[\]<>#|~$:=-])")
_AUTOLINK = re.compile(
    r"(?<=www)(?=\.)|(?=@)|(?<=:)(?=//)", re.IGNORECASE
)  # where a zero-width space keeps GFM from linking


UPLOADS = {  # kind -> (suffixes, the MIME type also named, what it is called): Streamlit checks no suffix of an
    # uploader that names a MIME type (on Windows its check split a name at a colon and crashed the page), so the
    # app checks it (wrong_upload)
    "structure": (STRUCTURE_SUFFIXES, "chemical/x-cif", "a structure file"),
    "layout": ((".yaml", ".yml"), "application/yaml", "a layout file"),
    "names": ((".yaml", ".yml"), "application/yaml", "a names file"),
    "alignment": (ALIGNMENT_SUFFIXES, "text/x-fasta", "an alignment"),
}


def upload_types(kind: str) -> list[str]:
    """An uploader's `type`: its suffixes and a MIME type."""
    suffixes, mime, _ = UPLOADS[kind]
    return [s.lstrip(".") for s in suffixes] + [mime]


def wrong_upload(name, kind: str) -> str:
    """'' when the upload's suffix is one of its kind's, else what to tell the visitor."""
    suffixes, _, what = UPLOADS[kind]
    if upload_suffix(name) in suffixes:
        return ""
    listed = ", ".join(suffixes[:-1]) + f" or {suffixes[-1]}"
    return f"{upload_name(name)} is not {what}: upload a {listed} file."


def upload_name(name, limit: int = 120) -> str:
    """An uploaded file's name as the app keeps and writes it: control, format and line/paragraph-separator characters
    dropped (a newline in a name would end a YAML comment and start live keys), at most `limit` characters."""
    kept = "".join(c for c in str(name) if unicodedata.category(c) not in {"Cc", "Cf", "Zl", "Zp"})
    return kept[:limit] or "upload"


def shown(text, limit: int = 400, path=None, name: str = "", paths: bool = True) -> str:
    """Text that came from a visitor or from an error (an id, a file name, a YAML key, an exception's message), made
    safe to show in any Streamlit text element: the loaded file's `path` given as its `name`, any other absolute path
    as its last part (paths=False keeps them: a path the visitor named on their own machine), at most `limit`
    characters, and Markdown escaped, so it can neither link, embed nor style anything."""
    text = str(text)
    if path:
        text = text.replace(str(path), name)
    if paths:
        text = _PATH.sub(lambda m: re.split(r"[/\\]", m.group().rstrip("/\\"))[-1], _FILE_URL.sub("/", text))
    if len(text) > limit:
        text = text[: limit - 1] + "…"
    return _AUTOLINK.sub("​", _MARKDOWN.sub(r"\\\1", text).replace("](", "]\\("))
