"""Command line: `foldmap summary FILE` prints secondary structure and sheet ordering."""

from __future__ import annotations

import argparse
import sys

import numpy as np

from . import layoutfile
from .dssp import assign_dssp
from .features import sse_contacts
from .frame import dna_radial_frame, dna_view_frame, symmetric_frame, view_frame
from .io import load_backbone, load_links, load_nucleic
from .layout import build_layout, resolve
from .palette import PALETTES
from .render import LOOP_STYLES, STYLES, draw, save
from .route import route_loops
from .sheets import build_sheets
from .ss import build_sses
from .style import CHOICES, LAYOUT_KEYS, SCALES, THEME_KEYS, THEMES, Style, describe, resolve_style
from .symmetry import TOLERANCE, detect_symmetry


def summarize(path: str) -> str:
    bb = load_backbone(path)
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, short_helices=True)
    sheets = build_sheets(sses, dssp.bridges)
    lines = [f"{path}: {len(bb)} residues, {len(set(l.chain for l in bb.labels))} chain(s)"]
    lines.append(
        f"helices: {sum(s.kind == 'H' for s in sses)}  strands: {sum(s.kind == 'E' for s in sses)}"
        f"  3₁₀ helices: {sum(s.kind == 'G' for s in sses)}"
    )
    for s in sses:
        if s.kind == "H":
            lines.append(f"  helix  {s.id}")
    for n, sh in enumerate(sheets, 1):
        flags = ("closed " if sh.closed else "") + ("ambiguous" if sh.ambiguous else "")
        lines.append(f"sheet {n}: {len(sh.strands)} strand(s) {flags}".rstrip())
        parts = []
        for k, (s, d) in enumerate(zip(sh.strands, sh.directions)):
            parts.append(f"{s.id}({'+' if d > 0 else '-'})")
            if k < len(sh.pair_kinds):
                parts.append({"P": "=par=", "A": "=anti=", "?": "=?="}[sh.pair_kinds[k]])
        lines.append("  " + " ".join(parts))
    try:
        sym = detect_symmetry(bb, sses)
    except ValueError:
        sym = None
    if sym is not None:
        extra = f"; twist {np.degrees(sym.twist):.1f}°, rise {sym.rise:.1f} Å" if sym.kind == "H" else ""
        lines.append(
            f"symmetry: {sym.label} — protomers {' | '.join(','.join(p) for p in sym.protomers)}"
            f"; CA RMSD {sym.rmsd:.2f} Å{extra}" + (f"; unpaired chains {','.join(sym.others)}" if sym.others else "")
        )
    na = load_nucleic(path, bb)
    if na.strands:
        kind = (
            "RNA" if all(st.rna for st in na.strands) else "DNA" if not any(st.rna for st in na.strands) else "DNA/RNA"
        )
        lines.append(
            f"{kind}: {len(na.strands)} strand(s) {', '.join(st.id for st in na.strands)}; "
            f"{len(na.pairs)} base pairs; {len(na.contacts)} protein residues make contacts"
        )
    return "\n".join(lines)


def make_layout(
    path,
    mode: str = "projected",
    rotate: float = 0.0,
    flip_v: bool = False,
    dna: bool = True,
    look: Style | None = None,
    symmetry: str = "auto",
    protomers: list[list[str]] | None = None,
    symmetry_tol: float = TOLERANCE,
    assembly: str = "auto",
    rename: dict | None = None,
    swap: list | None = None,
    move: list | None = None,
    up=None,
    view=None,
    domains=None,
):
    bb = load_backbone(path, assembly)
    dssp = assign_dssp(bb)
    sses = build_sses(bb, dssp.ss, short_helices=(look or Style()).helices_310)
    sheets = build_sheets(sses, dssp.bridges)
    nucleic = load_nucleic(path, bb, assembly) if dna else None
    sym = detect_symmetry(bb, sses, symmetry, protomers, symmetry_tol) if sses else None
    if nucleic is not None and nucleic.strands:
        frame = (
            dna_radial_frame(sses, nucleic)
            if not rotate and not flip_v
            else dna_view_frame(sses, nucleic, rotate, flip_v)
        )
    elif (
        sym is not None and not rotate and not flip_v and up is None and view is None
    ):  # an assembly: every protomer drawn alike, side by side
        frame = symmetric_frame(sses, sym)
    else:
        frame = view_frame(sses, rotate, flip_v, up, view)
    links = load_links(path, bb, assembly)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, links).items():
        contacts[key] = contacts.get(key, 0) + n
    named = _domains(domains, sses, contacts)
    for _, ids in named:  # a domain holds together: its elements pull on each other
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                key = frozenset((ids[a], ids[b]))
                contacts[key] = contacts.get(key, 0) + _DOMAIN_PULL
    lay = build_layout(
        sses,
        sheets,
        frame,
        mode,
        contacts=contacts,
        nucleic=nucleic,
        style=look,
        bridges=_bridge_springs(bb, sses, links),
        domains=named,
        rename=rename,
        swap=swap,
        move=move,
    )
    lay.focus = _focus((look or Style()).highlight, bb, sym)
    lay.links = links
    lay.res_chain = [l.chain for l in bb.labels]
    lay.res_name = [l.name for l in bb.labels]
    lay.res_b = [float(x) for x in bb.b] if bb.b is not None else []
    return lay, sses, bb


_DOMAIN_PULL = 6  # CA-pair-equivalents of attraction between any two elements of one domain


def _domains(spec, sses, contacts) -> list[tuple[str, list[str]]]:
    """[(name, element ids)] from 'auto' or [(name, [refs])]; every element in at most one domain."""
    from .features import auto_domains
    from .layout import provisional

    if not spec:
        return []
    if spec == "auto":
        return auto_domains(sses, contacts)
    lay, out, seen = provisional(sses), [], set()
    for name, refs in spec:
        ids = []
        for ref in refs:
            ids += [k for k in resolve(lay, sses, ref) if k not in seen and k not in ids]
        seen |= set(ids)
        out.append((name, ids))
    return out


def _cys_place(bb, sses, r: int):
    """(element id, fraction along it) for a cysteine: its own element, or the nearest element end."""
    chain = bb.labels[r].chain
    mine = [s for s in sses if s.chain == chain]
    if not mine:
        return None
    s = min(mine, key=lambda s: 0 if s.start <= r <= s.end else min(abs(r - s.start), abs(r - s.end)))
    if s.start <= r <= s.end:
        return s.id, (r - s.start) / max(s.end - s.start, 1)
    return s.id, (1.0 if r > s.end else 0.0)


def _bridge_springs(bb, sses, links) -> list[tuple[str, float, str, float]]:
    out = []
    for i, j in links.disulfides:
        a, b = _cys_place(bb, sses, i), _cys_place(bb, sses, j)
        if a and b and a[0] != b[0]:
            out.append((a[0], a[1], b[0], b[1]))
    return out


def _bridge_contacts(bb, sses, links) -> dict[frozenset[str], int]:
    """A disulfide ties two elements together like a full contact: the cysteine's own element, or the element
    nearest it along its chain when it sits in a loop."""

    def element(r: int):
        chain = bb.labels[r].chain
        mine = [s for s in sses if s.chain == chain]
        if not mine:
            return None
        return min(mine, key=lambda s: 0 if s.start <= r <= s.end else min(abs(r - s.start), abs(r - s.end))).id

    out: dict[frozenset[str], int] = {}
    for i, j in links.disulfides:
        a, b = element(i), element(j)
        if a and b and a != b:
            out[frozenset((a, b))] = out.get(frozenset((a, b)), 0) + 20  # as strong as the contact pull gets
    return out


def _focus(highlight: str, bb, sym) -> set[str] | None:
    """Chains to keep in colour: the deposited ones, the first protomer, or the user's list."""
    chains = {l.chain for l in bb.labels}
    if highlight == "none":
        return None
    if highlight == "asu":
        if bb.asu_chains and bb.asu_chains != chains:
            return set(bb.asu_chains)
        return set(sym.protomers[0]) if sym else None  # everything was deposited: the first protomer stands in
    if highlight == "protomer":
        return set(sym.protomers[0]) if sym else None
    wanted = set(highlight.split(","))
    if not wanted <= chains:
        raise ValueError(f"highlight names chains that are not here: {', '.join(sorted(wanted - chains))}")
    return wanted


def make_figure_and_loops(
    path,
    mode: str = "projected",
    rotate: float = 0.0,
    flip_v: bool = False,
    title: str | None = None,
    *,
    palette: str | None = None,
    style: str | None = None,
    loops: str | None = None,
    dna: bool = True,
    look: Style | None = None,
    symmetry: str = "auto",
    protomers: list[list[str]] | None = None,
    symmetry_tol: float = TOLERANCE,
    assembly: str = "auto",
    **adjust,
):
    layout, sses, bb = make_layout(
        path, mode, rotate, flip_v, dna, look, symmetry, protomers, symmetry_tol, assembly, **adjust
    )
    routed = route_loops(layout, sses, bb)
    return draw(layout, routed, sses, title, look=look, palette=palette, style=style, loop_style=loops), routed


def make_figure(path, *args, **kwargs):
    return make_figure_and_loops(path, *args, **kwargs)[0]


def _vector(text: str | None, what: str):
    if text is None:
        return None
    try:
        v = [float(x) for x in text.split(",")]
    except ValueError:
        v = []
    if len(v) != 3:
        raise ValueError(f"--{what} needs three numbers like 0,0,1; got {text!r}")
    return v


def _adjustments(args) -> dict:
    rename, swap, move = {}, [], []
    for item in args.rename:
        ref, sep, name = item.partition("=")
        if not sep or not name:
            raise ValueError(f"--rename needs REF=NAME, got {item!r}")
        rename[ref] = name
    for item in args.swap:
        parts = [p for p in item.split(",") if p]
        if len(parts) != 2:
            raise ValueError(f"--swap needs two elements separated by a comma, got {item!r}")
        swap.append((parts[0], parts[1]))
    for item in args.move:
        ref, sep, delta = item.partition("=")
        d = [float(x) for x in delta.split(",")] if sep else []
        if len(d) != 2:
            raise ValueError(f"--move needs REF=DX,DY, got {item!r}")
        move.append((ref, (d[0], d[1])))
    return {
        "rename": rename,
        "swap": swap,
        "move": move,
        "up": _vector(args.up, "up"),
        "view": _vector(args.view, "view"),
    }


def _figure_spec(args):
    """(theme, look, options, edits): a layout file (if given) overlaid by the command line, then defaults."""
    doc = layoutfile.load(args.layout_file) if args.layout_file else {}
    theme = args.theme or doc.get("theme") or "publication"
    file_style = [f"{k}={v}" for k, v in (doc.get("style") or {}).items()]
    shorthand = [
        f"{k}={v}" for k, v in (("palette", args.palette), ("fill", args.style), ("loops", args.loops)) if v is not None
    ]
    look = resolve_style(theme, args.style_file, file_style + shorthand + args.set)
    lay = doc.get("layout") or {}

    def pick(name, default):
        value = getattr(args, name, None)
        return value if value is not None else lay.get(name, default)

    cli = _adjustments(args)
    vec = lambda v, what: _vector(v, what) if isinstance(v, str) else v  # noqa: E731
    protomers = args.protomers or lay.get("protomers")
    opts = {
        "mode": pick("mode", "projected"),
        "rotate": float(pick("rotate", 0.0)),
        "flip_v": bool(args.flip_v or lay.get("flip_v", False)),
        "symmetry": pick("symmetry", "auto"),
        "symmetry_tol": float(pick("symmetry_tol", TOLERANCE)),
        "assembly": str(pick("assembly", "auto")),
        "protomers": _protomers(protomers) if isinstance(protomers, str) else protomers,
        "up": cli["up"] or vec(lay.get("up"), "up"),
        "view": cli["view"] or vec(lay.get("view"), "view"),
        "title": args.title if args.title is not None else lay.get("title"),
    }
    domains = _domain_spec(args.domain) or args.domains or doc.get("domains") or None
    if isinstance(domains, dict):
        domains = [(name, list(refs)) for name, refs in domains.items()]
    opts["domains"] = domains
    edits = layoutfile.edits_of(doc)
    edits["rename"].update(cli["rename"])
    edits["swap"] += cli["swap"]
    edits["move"] += cli["move"]
    return theme, look, opts, edits


def _domain_spec(items: list[str]) -> list[tuple[str, list[str]]]:
    out = []
    for item in items:
        name, sep, refs = item.partition("=")
        if not sep or not name.strip() or not refs.strip():
            raise ValueError(f"--domain needs NAME=REF[,REF...], got {item!r}")
        out.append((name.strip(), [r.strip() for r in refs.split(",") if r.strip()]))
    return out


def _protomers(text: str | None) -> list[list[str]] | None:
    if not text:
        return None
    return [[c.strip() for c in group.split(",") if c.strip()] for group in text.split(";") if group.strip()]


def styles_help() -> str:
    def kind(key: str) -> str:
        if key in CHOICES:
            return "one of " + ", ".join(CHOICES[key])
        if key in SCALES:
            return "positive number, 1 = default size"
        return "true/false"

    lines = ["Themes (--theme NAME) change only how things are drawn:"]
    for name, look in THEMES.items():
        changed = {k: v for k, v in vars(look).items() if v != getattr(Style(), k)}
        lines.append(f"  {name:13s}" + (", ".join(f"{k}={v}" for k, v in changed.items()) or "the defaults"))
    lines.append("\nTheme keys (--set KEY=VALUE, or KEY: VALUE in a --style-file):")
    lines += [f"  {k:17s} {kind(k)}" for k in THEME_KEYS]
    lines.append("\nLayout keys (same syntax; these move elements):")
    lines += [f"  {k:17s} {kind(k)}" for k in LAYOUT_KEYS]
    lines.append("\nDefaults:\n" + "\n".join("  " + l for l in describe(Style()).splitlines()))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="foldmap")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("summary", help="print secondary structure elements and sheet ordering")
    p.add_argument("structure")
    sub.add_parser("styles", help="list themes and every style key")
    sub.add_parser("app", help="open the Streamlit app in your browser (needs: pip install streamlit)")
    q = sub.add_parser("plot", help="draw the topology figure")
    iq = sub.add_parser("interactive", help="an HTML page linking the topology, the contact map and the 3D model")
    for cmd in (q, iq):
        q = cmd
        q.add_argument("structure")
        q.add_argument("-o", "--output", action="append", required=True, help=".svg, .pdf or .png (repeatable)")
        q.add_argument("--mode", choices=["projected", "stack"], help="projected (default) or stack")
        q.add_argument("--rotate", type=float, help="turn the picture counter-clockwise (degrees)")
        q.add_argument("--layout-file", metavar="YAML", help="start from a saved layout (theme, style, edits)")
        q.add_argument("--save-layout", metavar="YAML", help="write everything that defines this figure to a file")
        q.add_argument("--flip-v", action="store_true", help="mirror top/bottom (changes handedness)")
        q.add_argument("--title")
        q.add_argument(
            "--theme",
            "--preset",
            dest="theme",
            choices=list(THEMES),
            help="a ready-made look (default: publication); see `foldmap styles`",
        )
        q.add_argument("--style-file", help="YAML file of style keys (may name a theme); see `foldmap styles`")
        q.add_argument(
            "--set",
            action="append",
            default=[],
            metavar="KEY=VALUE",
            help="one style key, e.g. helix_scale=0.8 or loops=curved (repeatable; wins over the rest)",
        )
        q.add_argument("--palette", choices=list(PALETTES), help="chain colours (same as --set palette=...)")
        q.add_argument("--style", choices=STYLES, help="element fill (same as --set fill=...)")
        q.add_argument("--loops", choices=LOOP_STYLES, help="right-angled or smooth loops (same as --set loops=...)")
        q.add_argument("--no-dna", action="store_true", help="leave out DNA/RNA even when the file has it")
        q.add_argument(
            "--symmetry",
            metavar="auto|off|Cn|Dn|helical",
            help="auto finds cyclic, dihedral or helical symmetry and draws every protomer alike, side by side; "
            "name one (C3, D2, helical) to insist on it; off to ignore it",
        )
        q.add_argument(
            "--protomers",
            metavar="A,B;C,D",
            help="your own protomers, in order around the axis (chains separated by commas, protomers by ;)",
        )
        q.add_argument(
            "--assembly",
            metavar="auto|asu|ID",
            help="auto builds the file's first biological assembly when it adds copies (up to 24 chains); "
            "asu keeps the file as deposited; or give an assembly id from the file",
        )
        q.add_argument(
            "--rename",
            action="append",
            default=[],
            metavar="REF=NAME",
            help="your own element name, e.g. res:167-182=Gd or A:B=B' (repeatable)",
        )
        q.add_argument(
            "--swap",
            action="append",
            default=[],
            metavar="REF,REF",
            help="exchange the places of two elements, e.g. α1,α3 or A:C,A:D (repeatable)",
        )
        q.add_argument(
            "--move",
            action="append",
            default=[],
            metavar="REF=DX,DY",
            help="nudge an element (page units, one strand spacing ~1.1), e.g. α2=1,-2 (repeatable)",
        )
        q.add_argument(
            "--domain",
            action="append",
            default=[],
            metavar="NAME=REF[,REF...]",
            help="a named domain panel, e.g. ZPN=res:A:331-440 (repeatable)",
        )
        q.add_argument("--domains", choices=["auto"], help="find domains from element contacts (D1, D2...)")
        q.add_argument(
            "--up", metavar="X,Y,Z", help="3D direction to put at the top of the page (e.g. a membrane normal)"
        )
        q.add_argument("--view", metavar="X,Y,Z", help="3D direction pointing at the viewer")
        q.add_argument(
            "--symmetry-tol",
            type=float,
            metavar="Å",
            help=f"how far copies may differ, CA RMSD (default {TOLERANCE:g}; raise it for low-resolution X-ray)",
        )
    args = parser.parse_args(argv)
    try:
        if args.command == "app":
            from .app import run

            run()
        if args.command == "styles":
            print(styles_help())
            return 0
        if args.command in ("plot", "interactive"):
            theme, look, opts, edits = _figure_spec(args)
            kw = dict(
                dna=not args.no_dna,
                look=look,
                symmetry=opts["symmetry"],
                protomers=opts["protomers"],
                symmetry_tol=opts["symmetry_tol"],
                assembly=opts["assembly"],
                up=opts["up"],
                view=opts["view"],
                domains=opts["domains"],
                **edits,
            )
            if args.command == "interactive":
                from .interactive import write_page

                page = dict(mode=opts["mode"], title=opts["title"], rotate=opts["rotate"], flip_v=opts["flip_v"])
                for out in args.output:
                    print(f"wrote {write_page(args.structure, out, **page, **kw)}")
            else:
                fig, routed = make_figure_and_loops(
                    args.structure, opts["mode"], opts["rotate"], opts["flip_v"], opts["title"], **kw
                )
                stuck = [f"{l.a_id}>{l.b_id}" for l in routed if l.fallback]
                if stuck:
                    print(
                        f"foldmap: warning: no clear route for {len(stuck)} loop(s), drawn as plain curves "
                        f"that may cross elements: {', '.join(stuck)}",
                        file=sys.stderr,
                    )
                for out in args.output:
                    print(f"wrote {save(fig, out)}")
            if args.save_layout:
                plain = {k: v for k, v in kw.items() if k not in ("rename", "swap", "move")}
                layout, sses, _ = make_layout(args.structure, opts["mode"], opts["rotate"], opts["flip_v"], **plain)
                doc = layoutfile.document(args.structure, theme, look, opts, edits, layout, sses)
                print(f"wrote {layoutfile.save(args.save_layout, doc)}")
            return 0
        print(summarize(args.structure))
    except (FileNotFoundError, ValueError) as err:
        print(f"foldmap: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
