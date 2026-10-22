# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""One description of a figure, shared by the command line, the app and the MCP server: the structure, how it is
laid out (layout options and edits) and how it is drawn (theme and style). A layout file supplies a starting
point; values given explicitly win; defaults fill the rest."""

from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
from pathlib import Path
from typing import Any

from . import layoutfile
from .layout import CONVENTIONS, RES_RANGE
from .render import draw
from .route import route_loops
from .style import THEMES, Style, hbonds_superseded, resolve_style
from .symmetry import TOLERANCE

MODES = ("projected", "stack")
MEMBRANES = ("off", "auto", "flip")


@dataclass
class FigureSpec:
    structure: Path
    theme: str = "publication"
    look: Style = field(default_factory=Style)
    title: str | None = None
    mode: str = "projected"
    rotate: float = 0.0
    flip_v: bool = False
    dna: bool = True
    symmetry: str = "auto"
    symmetry_tol: float = TOLERANCE
    protomers: list[list[str]] | None = None
    assembly: str = "auto"
    up: list[float] | None = None
    view: list[float] | None = None
    domains: str | list[tuple[str, list[str]]] | None = None
    msa: str | None = None
    msa_reference: str | None = None
    membrane: str = "off"
    uniprot: str | None = None
    chains: list[str] | None = None
    focus: str = "auto"
    rename: dict[str, str] = field(default_factory=dict)
    swap: list[tuple[str, str]] = field(default_factory=list)
    move: list[tuple[str, tuple[float, float]]] = field(default_factory=list)
    names: list[str] = field(default_factory=list)  # --names: deposited, beta or YAML files, later ones winning
    # --like: drawn like a reference, {structure, chain, match, layout} (cli.make_layout's like); None: on its own
    like: dict | None = None
    # the researcher's own data (plan 13d, userdata): a per-residue table, the column read (None: its first of
    # numbers) and a crosslinks table; paths as given
    data: str | None = None
    data_column: str | None = None
    crosslinks: str | None = None

    def layout_kwargs(self, edits: bool = True) -> dict:
        """Keyword arguments for cli.make_layout after (path, mode, rotate, flip_v)."""
        kw = {
            "dna": self.dna,
            "look": self.look,
            "symmetry": self.symmetry,
            "protomers": self.protomers,
            "symmetry_tol": self.symmetry_tol,
            "assembly": self.assembly,
            "up": self.up,
            "view": self.view,
            "domains": self.domains,
            "msa": self.msa,
            "msa_reference": self.msa_reference,
            "membrane": self.membrane,
            "uniprot": self.uniprot,
            "chains": self.chains,
            "focus": self.focus,
            "names": list(self.names),  # not an edit: the user's renames are written for the names
            "like": dict(self.like) if self.like else None,
            "data": self.data,
            "data_column": self.data_column,
            "crosslinks": self.crosslinks,
        }
        if edits:
            kw.update(self.edits())
        return kw

    def page_kwargs(self) -> dict:
        """Keyword arguments for interactive.build_page / write_page after the path."""
        return {
            "mode": self.mode,
            "title": self.title,
            "rotate": self.rotate,
            "flip_v": self.flip_v,
            **self.layout_kwargs(),
        }

    def options(self) -> dict:
        """The layout options as a layout file stores them (like: the keys given, the reference as given)."""
        out = {k: getattr(self, k) for k in layoutfile.LAYOUT}
        if self.like:
            out["like"] = {k: v for k, v in self.like.items() if v is not None}
        return out

    def edits(self) -> dict:
        return {"rename": dict(self.rename), "swap": list(self.swap), "move": list(self.move)}


_GIVEN = {f.name for f in fields(FigureSpec)} - {"structure", "theme", "look"}


def figure_spec(
    structure,
    *,
    theme: str | None = None,
    style: list[str] | tuple[str, ...] = (),
    palette: str | None = None,
    fill: str | None = None,
    loops: str | None = None,
    style_file: str | Path | None = None,
    layout_file: str | Path | None = None,
    **given: Any,
) -> FigureSpec:
    """The figure for `structure`: the layout file (if any), overlaid by every value in `given` that is not None,
    then the defaults. style: KEY=VALUE overrides, applied over the layout file's style (palette, fill and loops
    are shorthands for the style keys of those names, applied before them). Edits (rename, swap,
    move) from the layout file and from `given` are combined, the given ones last; so are names (a spec given twice
    keeps its last place)."""
    unknown = sorted(set(given) - _GIVEN)
    if unknown:
        raise ValueError(f"unknown figure option(s): {', '.join(unknown)}")
    doc = layoutfile.load(layout_file) if layout_file else {}
    theme = theme or doc.get("theme") or "publication"
    if theme not in THEMES:
        raise ValueError(f"unknown theme {theme!r}; choose from {', '.join(THEMES)}")
    file_style = [f"{k}={v}" for k, v in (doc.get("style") or {}).items()]
    shorthand = [f"{k}={v}" for k, v in (("palette", palette), ("fill", fill), ("loops", loops)) if v is not None]
    look = resolve_style(theme, style_file, file_style + shorthand + list(style))
    lay = {k: v for k, v in (doc.get("layout") or {}).items() if v is not None}
    defaults = FigureSpec(Path(structure))

    def pick(name: str):
        """The given value, else the file's, else the default; None and '' count as not given."""
        if name == "title" and given.get("title") == "":
            return None  # an explicit empty title means no title, whatever the file says
        for source in (given, lay):
            value = source.get(name)
            if value is not None and value != "":
                return "off" if value is False and name in ("membrane", "symmetry") else value
        return getattr(defaults, name)

    mode = pick("mode")
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r}; choose from {', '.join(MODES)}")
    membrane = str(pick("membrane"))
    if membrane not in MEMBRANES:
        raise ValueError(f"unknown membrane {membrane!r}; choose from {', '.join(MEMBRANES)}")
    saved = doc.get("domains")
    domains = given.get("domains") or saved or None
    if isinstance(domains, dict):  # {name: [refs]}; auto, uniprot and [(name, [refs])] are taken as they are
        pairs = [(name, list(refs)) for name, refs in domains.items()]
        domains = SavedDomains(pairs, source=str(doc.get("domain_source") or "")) if domains is saved else pairs
    edits = layoutfile.edits_of(doc)
    edits["rename"].update(given.get("rename") or {})
    edits["swap"] += [_pair(p) for p in given.get("swap") or []]
    moves = given.get("move") or []
    for ref, d in moves.items() if isinstance(moves, dict) else moves:
        if len(d) != 2:
            raise ValueError(f"a move needs two numbers (dx, dy), got {d!r} for {ref!r}")
        edits["move"].append((ref, (float(d[0]), float(d[1]))))
    protomers, chains = pick("protomers"), pick("chains")
    names = [*(doc.get("names") or []), *(given.get("names") or [])]
    names = [n for k, n in enumerate(names) if n not in names[k + 1 :]]
    user = _user_data(pick("data"), pick("data_column"), pick("crosslinks"))
    look = data_look(look, user["data"], theme, style_file, file_style + shorthand + list(style))
    return FigureSpec(
        structure=Path(structure),
        theme=theme,
        look=look,
        title=pick("title"),
        mode=mode,
        rotate=float(pick("rotate")),
        flip_v=parse_flag(pick("flip_v"), "flip_v"),
        dna=parse_flag(pick("dna"), "dna"),
        symmetry=str(pick("symmetry")),
        symmetry_tol=float(pick("symmetry_tol")),
        protomers=parse_protomers(protomers) if isinstance(protomers, str) else protomers,
        assembly=str(pick("assembly")),
        up=parse_vector(pick("up"), "up"),
        view=parse_vector(pick("view"), "view"),
        domains=domains,
        msa=pick("msa"),
        msa_reference=pick("msa_reference"),
        membrane=membrane,
        uniprot=pick("uniprot"),
        chains=[c.strip() for c in chains.split(",") if c.strip()] if isinstance(chains, str) else chains,
        focus=str(pick("focus")),
        names=names,
        like=parse_like(pick("like")),
        **user,
        **edits,
    )


LIKE_KEYS = ("structure", "chain", "match", "layout")  # FigureSpec.like, as --like and its options give it


def parse_like(value) -> dict | None:
    """FigureSpec.like from --like's options or a layout file's `like` (a mapping of LIKE_KEYS, structure required),
    checked and with every key: {structure, chain, match, layout}, match 'auto' when not given. None stays None."""
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError(
            f"like is a mapping (structure, chain, match), e.g. {{structure: 2HHB, chain: B}}; got {value!r}"
        )
    unknown = [k for k in value if k not in LIKE_KEYS]
    if unknown:
        raise ValueError(f"like: unknown key(s) {', '.join(map(str, unknown))}; use {', '.join(LIKE_KEYS)}")
    out = {}
    for key in LIKE_KEYS:
        v = value.get(key)
        if v is not None and (isinstance(v, bool) or not isinstance(v, (str, int, float))):
            raise ValueError(f"like: {key} is text; got {v!r}")
        out[key] = str(v) if v is not None and str(v) != "" else None
    if out["structure"] is None:
        raise ValueError("like: structure (the reference, a file or PDB id) is required")
    out["match"] = out["match"] or "auto"
    return out


def colouring_asked(style_file, overrides) -> bool:
    """Whether a colouring is asked for by name: color_by among the KEY=VALUE overrides or in the style file."""
    if any(item.partition("=")[0].strip() == "color_by" for item in overrides):
        return True
    if style_file is None:
        return False
    from .yamlfile import read_yaml

    loaded = read_yaml(Path(style_file).read_text(encoding="utf-8"), "style", str(style_file)) or {}
    return isinstance(loaded, dict) and "color_by" in loaded


def data_look(look: Style, data, theme: str, style_file=None, overrides=()) -> Style:
    """The look of a figure given the researcher's values (`data`, a path or None): drawn by them (color_by=data)
    unless a colouring is asked for (by the theme, an override or the style file); color_by=data with no data is a
    ValueError."""
    if data is not None and look.color_by == THEMES[theme].color_by == Style().color_by:
        if not colouring_asked(style_file, overrides):
            look = replace(look, color_by="data")
    if look.color_by == "data" and data is None:
        raise ValueError("color_by=data draws your own per-residue values: give them with --data FILE.csv")
    return look


def _user_data(data, column, crosslinks) -> dict:
    """FigureSpec's data, data_column and crosslinks (paths and a column name as text; a number YAML read stays a
    name), checked: a column needs its data file."""
    out = {}
    for key, value in (("data", data), ("data_column", column), ("crosslinks", crosslinks)):
        if value is not None and (isinstance(value, bool) or not isinstance(value, (str, int, float))):
            raise ValueError(f"{key} is text (a file path or a column name); got {layoutfile._brief(value)}")
        out[key] = str(value) if value is not None else None
    if out["data_column"] is not None and out["data"] is None:
        raise ValueError(f"data_column {out['data_column']!r} needs data (--data FILE.csv): the file it is a column of")
    return out


def parse_flag(value, what: str) -> bool:
    """A true/false value from a real boolean or yes/no/true/false text."""
    if isinstance(value, bool):
        return value
    word = value.strip().lower() if isinstance(value, str) else None
    if word in ("true", "yes"):
        return True
    if word in ("false", "no"):
        return False
    raise ValueError(f"{what} must be true or false, got {value!r}")


def _pair(item) -> tuple[str, str]:
    """Two element references from a pair or from 'REF,REF' text."""
    parts = [p.strip() for p in item.split(",") if p.strip()] if isinstance(item, str) else list(item)
    if len(parts) != 2:
        raise ValueError(f"a swap needs two elements, like A,B; got {item!r}")
    return parts[0], parts[1]


def parse_vector(value, what: str) -> list[float] | None:
    """Three numbers from 'X,Y,Z' text or a sequence; None stays None."""
    if value is None:
        return None
    if isinstance(value, str):
        try:
            v = [float(x) for x in value.split(",")]
        except ValueError:
            v = []
    else:
        v = [float(x) for x in value]
    if len(v) != 3:
        raise ValueError(f"--{what} needs three numbers like 0,0,1; got {value!r}")
    return v


def parse_protomers(text: str | None) -> list[list[str]] | None:
    if not text:
        return None
    return [[c.strip() for c in group.split(",") if c.strip()] for group in text.split(";") if group.strip()]


class SavedDomains(list):
    """Domains read from a layout file, [(name, [refs])]. A file outlives the structure it was saved from, so
    element ids it names that the structure no longer has are dropped with a warning (cli._domains); given on the
    command line, they are an error. quiet: the caller shows the warning itself (cli.stale_domains_note); none is
    printed. source: where the saved domains came from (the file's domain_source)."""

    def __init__(self, items=(), quiet: bool = False, source: str = ""):
        super().__init__(items)
        self.quiet = quiet
        self.source = source  # the saved figure's Layout.domain_source: "inferred" domains stay said to be


def parse_domains(items: list[str]) -> list[tuple[str, list[str]]]:
    out = []
    for item in items:
        name, sep, refs = item.partition("=")
        if not sep or not name.strip() or not refs.strip():
            raise ValueError(f"--domain needs NAME=REF[,REF...], got {item!r}")
        out.append((name.strip(), [r.strip() for r in refs.split(",") if r.strip()]))
    return out


@dataclass(frozen=True)
class Problem:
    """What is wrong with one element reference of a figure (check_references), said before it is laid out."""

    where: str  # the option and the item: "--domain N=A:1-60", "--rename α1=Gd", "--swap α1,α3", "mark α3"
    text: str  # what is wrong and the accepted forms
    fatal: bool  # False: a warning (A:1-60 read as res:A:1-60), the figure is drawn

    def __str__(self) -> str:
        return f"{self.where}: {self.text}"


FORMS = "a label (α1, A), chain:label (A:G), #k (the k-th element) or a residue range (res:A:1-60)"
RANGE_FORMS = "res:150-159, res:A:150-159, res:A:-7-4 or res:H:100A-100K"


def _range(text: str) -> tuple[int, int] | None:
    """(first, last) residue numbers of a range written without res: (1-60, -7-4, 100A-100K); None for anything
    else, a single number included (it may be a label)."""
    m = RES_RANGE.fullmatch(text)
    return (int(m[1]), int(m[3])) if m is not None and m[3] is not None else None


def _reference(ref: str, where: str, domain: bool, elements) -> tuple[list[Problem], tuple | None]:
    """Pass 1 for one reference (layout.resolve's forms): its problems, and what pass 2 checks against the structure,
    (chain or None, (first, last) residue numbers or None); None when it names no chain and no residues. A range
    written without res: (A:1-60, 1-60) that is no element id (`elements`) is read as one, with a warning. A hyphen
    in what can only be a label is an error in a domain (generic labels have none); elsewhere a name given by
    --names or --rename may hold one."""
    ref = ref.strip()
    if not ref:
        return [Problem(where, f"an empty reference; name an element as {FORMS}", True)], None
    if ref.startswith("res:"):
        chain, body = ref[4:].split(":") if ref[4:].count(":") == 1 else (None, ref[4:])
        m = RES_RANGE.fullmatch(body) if chain != "" else None
        if m is None:
            return [Problem(where, f"{ref} is no residue range; write {RANGE_FORMS}", True)], None
        return [], (chain, (int(m[1]), int(m[3] if m[3] is not None else m[1])))
    if ref.startswith("#"):
        if not ref[1:].isdigit() or int(ref[1:]) < 1:
            return [Problem(where, f"{ref}: #k counts elements from the N terminus, #1 the first", True)], None
        return [], None
    prefix = ref.split(":", 1)[0]
    if prefix in CONVENTIONS and ":" in ref:  # dep:helix A, beta:A:α4 (Q-A12): a name, its chain checked in pass 2
        chain, name = ref[len(prefix) + 1 :].split(":", 1) if ref.count(":") > 1 else (None, ref[len(prefix) + 1 :])
        if not name.strip() or (chain is not None and not chain.strip()):
            return [Problem(where, f"{ref} names no element; write {prefix}:NAME or {prefix}:CHAIN:NAME", True)], None
        return [], ((chain.strip(), None) if chain is not None else None)
    chain, label = ref.split(":", 1) if ":" in ref else (None, ref)
    if chain is not None and (not chain.strip() or not label.strip()):
        return [Problem(where, f"{ref} names no element; name one as {FORMS}", True)], None
    span = _range(label)
    if span is not None:
        if ref in elements:  # an element id (a layout file's, an assembly copy's: A1:3-7): the element itself
            return [], (chain, span)
        return [Problem(where, f"read {ref} as the residue range res:{ref}", False)], (chain, span)
    if domain and "-" in label:
        return [
            Problem(
                where,
                f"{ref} is neither a residue range nor an element; write a range as res:{chain or 'A'}:FIRST-LAST "
                f"(e.g. res:{chain or 'A'}:1-60), an element as {FORMS}",
                True,
            )
        ], None
    return [], ((chain, None) if chain is not None else None)


def _references(spec: FigureSpec) -> tuple[list[Problem], list[tuple[str, str, bool]]]:
    """Every element reference the figure names, as (where, reference, in a domain), with the problems of the
    domains themselves (one with no name or no reference). A layout file's domains (SavedDomains: element ids it
    outlived are dropped with a warning) and automatic ones are not references."""
    problems, refs = [], []
    domains = spec.domains
    if isinstance(domains, list) and not isinstance(domains, SavedDomains):
        for name, items in domains:
            items = [items] if isinstance(items, str) else list(items or [])
            if not str(name).strip():
                where = f"--domain ={','.join(map(str, items))}"
                problems.append(Problem(where, "a domain needs a name: NAME=REF[,REF...]", True))
            if not items:
                text = f"a domain needs elements: NAME=REF[,REF...], each REF {FORMS}"
                problems.append(Problem(f"--domain {name}=", text, True))
            refs += [(f"--domain {name}={ref}", str(ref), True) for ref in items]
    refs += [(f"--rename {ref}={new}", str(ref), False) for ref, new in (spec.rename or {}).items()]
    refs += [(f"--swap {a},{b}", str(r), False) for a, b in spec.swap or [] for r in (a, b)]
    refs += [(f"--move {ref}", str(ref), False) for ref, _ in spec.move or []]
    marks = (getattr(spec.look, "mark", "") or "").split(",")
    refs += [(f"mark {entry.strip()}", entry.partition("=")[0], False) for entry in marks if entry.strip()]
    return problems, refs


def check_references(spec: FigureSpec, bb=None, *, drawn: set[str] | None = None, elements=()) -> list[Problem]:
    """What is wrong with the element references the figure names (domains given, rename, swap, move, mark), said
    before it is laid out (Q-A19). Pass 1 (no structure): each reference's form; a range written without res:
    (A:1-60, 1-60) is read as res:A:1-60 with a warning, unless it is one of `elements` (element ids: an assembly
    copy's A1:3-7). Pass 2 (`bb`, the protein backbone as read, before --chains): the chain a reference names
    exists, is drawn (--chains; `drawn`: the chains --focus keeps, None when it keeps all) and its residue range
    meets residues of it. A reference that passes may still name no element (a loop's residues): the layout says
    so. Pass 1 also warns when the style gives both hbonds and interactions (style.hbonds_superseded)."""
    problems, refs = _references(spec)
    look = spec.look or Style()
    superseded = hbonds_superseded(look.hbonds, look.interactions)
    if superseded:  # a style key the other overrides (plan 13c): said with the references, a warning
        problems.append(Problem(f"hbonds={look.hbonds}", superseded, False))
    elements = set(elements)
    checked = []
    for where, ref, domain in refs:
        found, target = _reference(ref, where, domain, elements)
        problems += found
        if target is not None and not any(p.fatal for p in found):
            checked.append((where, ref.strip(), *target))
    if bb is None:
        return problems
    numbers: dict[str, list[int]] = {}
    for label in bb.labels:
        numbers.setdefault(label.chain, []).append(label.seq)
    present = list(numbers)
    asked = [c for c in spec.chains or [] if c in numbers]  # protein chains named: only those are drawn
    kept = [c for c in present if (not asked or c in asked) and (drawn is None or c in drawn)]
    for where, ref, chain, span in checked:
        if chain is not None and chain not in numbers:
            text = f"no protein chain {chain} in the structure; its protein chains are {', '.join(present) or 'none'}"
            problems.append(Problem(where, text, True))
            continue
        if chain is not None and chain not in kept:
            if asked and chain not in asked:
                text = f"chain {chain} is not drawn: --chains draws {', '.join(asked)} only"
            else:
                text = (
                    f"chain {chain} is not drawn: --focus {spec.focus} draws chain(s) {', '.join(kept)} "
                    "(add --focus none to draw every chain)"
                )
            problems.append(Problem(where, text, True))
            continue
        if span is None:
            continue
        lo, hi = min(span), max(span)
        chains = [chain] if chain is not None else kept
        if not any(min(numbers[c]) <= hi and max(numbers[c]) >= lo for c in chains):
            if chain is not None:
                said = f"chain {chain}: its residues are numbered {min(numbers[chain])} to {max(numbers[chain])}"
            else:
                held = ", ".join(f"{c} {min(numbers[c])} to {max(numbers[c])}" for c in chains) or "none"
                said = f"any chain drawn (numbered {held})"
            problems.append(Problem(where, f"{ref} meets no residue of {said}", True))
    return problems


@dataclass
class Drawn:
    """A drawn figure and what it was drawn from."""

    figure: Any  # matplotlib Figure
    loops: list
    layout: Any  # layout.Layout
    sses: list
    bb: Any  # model.Backbone
    model: Any = None  # io.read_model result the figure was built from

    @property
    def beads(self):
        """render.ContactBeads: the contact marks drawn (for nacheck.verify_contacts), None without any."""
        return getattr(self.figure, "foldmap_beads", None)

    def verify_contacts(self):
        """nacheck.ContactReport of the drawn contact beads against the coordinates."""
        from .nacheck import verify_contacts

        return verify_contacts(self.model, self.beads, self.layout, self.bb) if self.beads else None

    @property
    def provenance(self) -> list[dict]:
        """How the figure was made (provenance.figure_provenance): what in it is inferred, ambiguous or moved for
        readability, each {kind, text, count, elements}; [] when nothing beyond the standard method."""
        return [dict(i) for i in getattr(self.figure, "foldmap_provenance", [])]

    @property
    def stuck(self) -> list[str]:
        """Loops with no clear route, drawn as plain curves, then nucleic-acid strands drawn with no clear way."""
        from .provenance import stuck

        return stuck(self.loops, self.layout)

    @property
    def crowded_marks(self) -> list[str]:
        """Ligand marks (name:chain+number) with no clear spot anywhere near, drawn over an element or text at the
        least bad spot (render._least_bad_spot)."""
        return list(getattr(self.figure, "foldmap_crowded_marks", []))

    @property
    def page_notes(self) -> list[str]:
        """What is said of the figure's page (render.page_notes): held at layout.MAX_PAGE_IN with its text under
        6 pt; its PNG at less than 300 dpi."""
        from .render import page_notes

        return page_notes(self.figure)

    @property
    def sheet_notes(self) -> list[str]:
        """What is said of the sheets as drawn (layout.sheet_notes): barrels drawn apart that pair; the strands off
        a sheet's run of paired strands, each with its partner; the pairings not drawn side by side."""
        return list(getattr(self.layout, "sheet_notes", []))

    @property
    def over_labels(self) -> list[str]:
        """Loops routed over a label or a chain-end stub, the only route found (route.Loop.crosses_labels)."""
        from .provenance import over_labels

        return over_labels(self.loops)


def draw_spec(spec: FigureSpec, like_out: dict | None = None) -> Drawn:
    """Lay out, route and draw the figure once. `like_out`: cli.make_layout's (a figure drawn like a reference: what
    the reference gave, for the H checks)."""
    from .cli import make_layout
    from .io import read_model

    model = read_model(spec.structure, spec.assembly)
    layout, sses, bb = make_layout(
        spec.structure, spec.mode, spec.rotate, spec.flip_v, model=model, like_out=like_out, **spec.layout_kwargs()
    )
    loops = route_loops(layout, sses, bb)
    return Drawn(draw(layout, loops, sses, spec.title, model=model, look=spec.look), loops, layout, sses, bb, model)


def layout_document(spec: FigureSpec, structure_name: str | None = None, drawn: Drawn | None = None) -> dict:
    """The layout file document for the figure. Element references in the edits are resolved against the
    layout without the edits (the labels they were written for) and saved as residue ranges. With no edits that is
    the layout `drawn` (the figure of this spec, draw_spec) when given, which is not laid out again (random-100
    C33: 7JU4's figure was laid out twice, 27 of its 183 s)."""
    from .cli import make_layout

    if drawn is not None and not any(spec.edits().values()):
        layout, sses = drawn.layout, drawn.sses
    else:
        kw = {**spec.layout_kwargs(edits=False), "data": None, "data_column": None, "crosslinks": None}  # read once
        layout, sses, _ = make_layout(spec.structure, spec.mode, spec.rotate, spec.flip_v, **kw)
    return layoutfile.document(
        structure_name or str(spec.structure),
        spec.theme,
        spec.look,
        spec.options(),
        spec.edits(),
        layout,
        sses,
        names=spec.names,
    )
