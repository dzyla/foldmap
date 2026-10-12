# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Everything about how a figure looks, in one place: named themes, a YAML style file, key=value overrides.

Theme keys change only how things are drawn; layout keys (LAYOUT_KEYS) change where they go. Themes set
theme keys only, so switching theme never moves an element."""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path

import yaml

from .hbonds import parse_selection
from .palette import PALETTES

CHOICES = {
    "palette": tuple(PALETTES),
    "fill": ("bold", "pale", "outline"),  # solid colour / pale fills / white fills with coloured edges
    "na_extent": ("auto", "contacts", "all"),
    "na_style": ("ladder", "ribbon", "block", "letters"),
    "color_by": (
        "chain",
        "sequence",
        "shade",
        "sstype",
        "bfactor",
        "hydropathy",
        "plddt",
        "conservation",
    ),  # chain / N->C ramp / chain hue
    # light N -> dark C / helix vs strand (Richardson) / mean B-factor, rigid -> flexible / Kyte-Doolittle
    "loop_color": (
        "black",
        "chain",
        "element",
        "residue",
    ),  # loops in black / their chain's colour / the element they leave / each residue's colour
    "sequence_map": ("turbo", "viridis", "plasma", "cividis"),
    "loops": ("orthogonal", "curved"),
    "helix_shading": ("depth", "none", "gloss"),  # back face darker / no shading / darker back + glint
    "helix_angle": ("snap", "upright", "tilted"),  # nearest of vertical/horizontal / always vertical / as in 3D
    "fold_template": ("auto", "off"),  # all-parallel α/β sheets as one βα row (TIM, Rossmann) / as projected
    "helix_on_sheet": ("auto", "off"),  # a helix packed on a sheet face drawn over the sheet / beside it, as projected
}
_ALIASES = {"fill": {"classic": "bold", "flat": "pale"}}
_KEY_ALIASES = {"dna_extent": "na_extent"}  # earlier names of style keys
LAYOUT_KEYS = (
    "helix_angle",
    "helix_scale",
    "strand_scale",
    "dna_scale",
    "helices_310",
    "na_extent",
    "fold_template",
    "helix_on_sheet",
)  # these move elements
NA_BACKBONES = ("modality", "chain", "ink")
MODALITIES = ("protein", "dna", "rna")


def backbone_overrides(text: str) -> dict[str, str]:
    """`backbone_colors` as {key: colour}: keys are chain ids or protein / dna / rna; ValueError names the key."""
    from matplotlib.colors import is_color_like

    out: dict[str, str] = {}
    for entry in filter(None, (e.strip() for e in text.split(","))):
        key, eq, colour = (part.strip() for part in entry.partition("="))
        if not eq or not re.fullmatch(r"[A-Za-z0-9]+", key):
            raise ValueError(
                f"backbone_colors entries are KEY=COLOUR, KEY a chain id or protein/dna/rna; got {entry!r}"
            )
        if not is_color_like(colour):
            raise ValueError(f"backbone_colors: {key} has {colour!r}, which is not a colour (use a name or #rrggbb)")
        out[key.lower() if key.lower() in MODALITIES else key] = colour
    return out


SCALES = ("helix_scale", "strand_scale", "dna_scale", "loop_width", "font_scale")


@dataclass(frozen=True)
class Style:
    palette: str = "okabe-ito"
    fill: str = "bold"
    loops: str = "orthogonal"
    helix_shading: str = "depth"
    helix_angle: str = "snap"
    fold_template: str = "auto"  # auto: an all-parallel α/β sheet (TIM barrel, Rossmann fold) is drawn as one row,
    # each crossover helix between the two strands it joins; off: helices where the projection puts them
    helix_on_sheet: str = "auto"  # auto: a helix packed on a sheet face (ubiquitin's α1) is drawn over the sheet,
    # semi-transparent, as it lies on it; off: beside the sheet, where the projection puts it
    helix_scale: float = 1.0  # helix ribbon diameter
    strand_scale: float = 1.0  # arrow width and the spacing of strands in a sheet
    dna_scale: float = 1.0  # duplex diameter
    loop_width: float = 1.0  # line weight of loops
    font_scale: float = 1.0
    labels: bool = True  # strand letters, helix names
    legend: bool = True
    color_by: str = "chain"
    sequence_map: str = "turbo"  # blue (N) -> red (C), like PyMOL's spectrum
    nucleotide_labels: bool = True  # base letters on DNA/RNA
    na_extent: str = "auto"  # DNA/RNA drawn: auto (all up to 300 nt, else the stems the protein touches),
    # contacts (the stretch the protein touches, +6 nt) or all
    na_style: str = "ladder"  # stems as ladder (rails and rungs), ribbon (double helix), block (bars) or letters
    helices_310: bool = True  # short 3-10 helices as small boxes (η1, η2...)
    sheet_panels: bool = True  # a light panel behind each β-sheet, so sandwiches read as two sheets
    loop_color: str = "black"
    loop_arrows: bool = False  # a small chevron on each loop pointing N -> C
    residue_numbers: bool = False  # first and last residue number at each element's ends
    highlight: str = "none"  # none / asu (as deposited) / protomer (the first) / chains "A,B": the rest in grey
    disulfides: bool = True  # yellow bars between bridged cysteines
    glycans: bool = True  # SNFG symbols on glycosylated residues
    background: str = "#ffffff"  # page colour
    ink: str = "#000000"  # text and plain loops
    uniprot_sites: str = "auto"  # with --uniprot: auto (active + binding sites), all, none, or types "active,modified"
    ligands: str = "auto"  # bound ligands/metal ions: auto (skip buffer additives), all, none, or codes "HEM,ZN"
    na_backbone: str = "modality"  # DNA/RNA backbones: modality (dna_color / rna_color), chain (one colour per chain),
    # ink (as loops) or one colour (#rrggbb or a name) for all nucleic acid
    dna_color: str = "#8c564b"  # DNA backbone with na_backbone=modality
    rna_color: str = "#6a3d9a"  # RNA backbone with na_backbone=modality
    backbone_colors: str = ""  # overrides, "A=#1f77b4,dna=#555555,I=crimson": a chain id or protein / dna / rna
    hbonds: str = (
        "off"  # dashed lines for hydrogen bonds: off, na (protein with nucleic acid), interchain or chains:A+B,C+D
    )
    hbond_color: str = "#d020f8"  # the dashed lines (a violet-magenta clear of every theme palette; see THEMES)
    hbond_max: int = 200  # at most this many lines, the shortest kept
    mark: str = ""  # elements in a colour of their own: "G", "A:G=#2ca02c", "#3", "res:150-159=red", comma-separated

    def validate(self) -> Style:
        for key, allowed in CHOICES.items():
            if getattr(self, key) not in allowed:
                raise ValueError(f"{key} must be one of {', '.join(allowed)}; got {getattr(self, key)!r}")
        from matplotlib.colors import is_color_like

        for key in ("background", "ink"):
            if not is_color_like(getattr(self, key)):
                raise ValueError(f"{key} must be a colour such as #0d2a4a or white; got {getattr(self, key)!r}")
        for key in ("dna_color", "rna_color"):
            if not is_color_like(getattr(self, key)):
                raise ValueError(f"{key} must be a colour such as #8c564b or brown; got {getattr(self, key)!r}")
        if self.na_backbone not in NA_BACKBONES and not is_color_like(self.na_backbone):
            raise ValueError(
                f"na_backbone must be {', '.join(NA_BACKBONES)} or a colour such as #8c564b; got {self.na_backbone!r}"
            )
        backbone_overrides(self.backbone_colors)
        parse_selection(self.hbonds)
        if not is_color_like(self.hbond_color):
            raise ValueError(f"hbond_color must be a colour such as #d020f8 or teal; got {self.hbond_color!r}")
        if isinstance(self.hbond_max, bool) or not isinstance(self.hbond_max, int) or self.hbond_max < 1:
            raise ValueError(f"hbond_max must be a whole number of 1 or more; got {self.hbond_max!r}")
        if not re.fullmatch(r"auto|all|none|[a-z-]+(,[a-z-]+)*", self.uniprot_sites):
            raise ValueError(
                f"uniprot_sites must be auto, all, none or types like active,binding; got {self.uniprot_sites!r}"
            )
        if not re.fullmatch(r"auto|all|none|[A-Za-z0-9]{1,5}(,[A-Za-z0-9]{1,5})*", self.ligands):
            raise ValueError(f"ligands must be auto, all, none or residue codes like HEM,ZN; got {self.ligands!r}")
        if not re.fullmatch(r"none|asu|protomer|[A-Za-z0-9]+(,[A-Za-z0-9]+)*", self.highlight):
            raise ValueError(f"highlight must be none, asu, protomer or chains like A,B; got {self.highlight!r}")
        for key in SCALES:
            if not getattr(self, key) > 0:
                raise ValueError(f"{key} must be a positive number; got {getattr(self, key)!r}")
        return self


THEMES = {
    "publication": Style(),
    "minimal": Style(fill="outline", helix_shading="none", loop_width=0.8, na_style="block"),
    "print": Style(
        palette="greys",
        fill="pale",
        helix_shading="none",
        na_style="block",
        dna_color="#4d4d4d",
        rna_color="#808080",
        hbond_color="#b3b3b3",
    ),
    "presentation": Style(loop_width=1.6, font_scale=1.4, helix_shading="gloss", loops="curved", na_style="ribbon"),
    "cartoon": Style(loops="curved", helix_shading="gloss", palette="tol-bright", na_style="ribbon"),
    "rainbow": Style(color_by="sequence", helix_shading="gloss"),
    "shaded": Style(color_by="shade"),
    "richardson": Style(color_by="sstype", helix_shading="depth", loops="curved"),
    "flexibility": Style(color_by="bfactor"),
    "hydropathy": Style(color_by="hydropathy"),
    "blueprint": Style(
        background="#0d2a4a",
        ink="#e8f1ff",
        palette="blueprint",
        fill="outline",
        helix_shading="none",
        loop_width=0.9,
        dna_color="#e8743b",
        rna_color="#a06cd5",
        hbond_color="#989890",
    ),
    "conservation": Style(color_by="conservation", loop_color="residue"),
    "alphafold": Style(color_by="plddt", loop_color="residue", loop_width=1.3),
    "goodsell": Style(palette="tol-bright", fill="pale", helix_shading="none", loop_width=1.3, sheet_panels=False),
    "journal": Style(palette="tol-muted", loop_width=0.8, font_scale=0.9, rna_color="#0070a0"),
    "trace": Style(color_by="shade", loop_color="element", loop_arrows=True),
}
PRESETS = THEMES  # earlier name
_TYPES = {f.name: f.type for f in fields(Style)}
THEME_KEYS = tuple(k for k in _TYPES if k not in LAYOUT_KEYS)


def _coerce(key: str, value):
    key = _KEY_ALIASES.get(key, key)
    if key not in _TYPES:
        raise ValueError(f"unknown style key {key!r}; valid keys: {', '.join(_TYPES)}")
    kind = _TYPES[key]
    value = _ALIASES.get(key, {}).get(value, value)
    if kind == "bool":
        if isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if text in ("true", "yes", "on", "1"):
            return True
        if text in ("false", "no", "off", "0"):
            return False
        raise ValueError(f"{key} must be true or false; got {value!r}")
    if kind == "int":
        try:
            number = float(value)
            if number != int(number):
                raise ValueError
        except (TypeError, ValueError, OverflowError):
            raise ValueError(f"{key} must be a whole number; got {value!r}") from None
        return int(number)
    if kind == "float":
        try:
            return float(value)
        except (TypeError, ValueError):
            raise ValueError(f"{key} must be a number; got {value!r}") from None
    return str(value)


def resolve_style(
    theme: str | None = None, style_file: str | Path | None = None, overrides: list[str] | tuple[str, ...] = ()
) -> Style:
    """Theme (or the file's `theme:`), then the file's keys, then key=value overrides; validated."""
    settings: dict = {}
    if style_file is not None:
        loaded = yaml.safe_load(Path(style_file).read_text(encoding="utf-8")) or {}
        if not isinstance(loaded, dict):
            raise ValueError(f"{style_file}: a style file is a mapping of key: value")
        named = loaded.pop("theme", None) or loaded.pop("preset", None)
        loaded.pop("preset", None)
        theme = theme or named
        settings.update({_KEY_ALIASES.get(k, k): _coerce(k, v) for k, v in loaded.items()})
    name = theme or "publication"
    if name not in THEMES:
        raise ValueError(f"unknown theme {name!r}; choose from {', '.join(THEMES)}")
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"style override {item!r} must be key=value")
        key, value = item.split("=", 1)
        settings[_KEY_ALIASES.get(key.strip(), key.strip())] = _coerce(key.strip(), value.strip())
    return replace(THEMES[name], **settings).validate()


def describe(style: Style) -> str:
    return "\n".join(f"{k}: {v}" for k, v in asdict(style).items())
