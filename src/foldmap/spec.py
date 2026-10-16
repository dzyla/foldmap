# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""One description of a figure, shared by the command line, the app and the MCP server: the structure, how it is
laid out (layout options and edits) and how it is drawn (theme and style). A layout file supplies a starting
point; values given explicitly win; defaults fill the rest."""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

from . import layoutfile
from .render import draw
from .route import route_loops
from .style import THEMES, Style, resolve_style
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
        """The layout options as a layout file stores them."""
        return {k: getattr(self, k) for k in layoutfile.LAYOUT}

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
        **edits,
    )


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


def draw_spec(spec: FigureSpec) -> Drawn:
    """Lay out, route and draw the figure once."""
    from .cli import make_layout
    from .io import read_model

    model = read_model(spec.structure, spec.assembly)
    layout, sses, bb = make_layout(
        spec.structure, spec.mode, spec.rotate, spec.flip_v, model=model, **spec.layout_kwargs()
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
        layout, sses, _ = make_layout(
            spec.structure, spec.mode, spec.rotate, spec.flip_v, **spec.layout_kwargs(edits=False)
        )
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
