"""Everything about how a figure looks, in one place: named themes, a YAML style file, key=value overrides.

Theme keys change only how things are drawn; layout keys (LAYOUT_KEYS) change where they go. Themes set
theme keys only, so switching theme never moves an element."""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path

import yaml

from .palette import PALETTES

CHOICES = {
    "palette": tuple(PALETTES),
    "fill": ("bold", "pale", "outline"),  # solid colour / pale fills / white fills with coloured edges
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
}
_ALIASES = {"fill": {"classic": "bold", "flat": "pale"}}
LAYOUT_KEYS = ("helix_angle", "helix_scale", "strand_scale", "dna_scale", "helices_310")  # these move elements
SCALES = ("helix_scale", "strand_scale", "dna_scale", "loop_width", "font_scale")


@dataclass(frozen=True)
class Style:
    palette: str = "okabe-ito"
    fill: str = "bold"
    loops: str = "orthogonal"
    helix_shading: str = "depth"
    helix_angle: str = "snap"
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
    helices_310: bool = True  # short 3-10 helices as small boxes (η1, η2...)
    sheet_panels: bool = True  # a light panel behind each β-sheet, so sandwiches read as two sheets
    loop_color: str = "black"
    loop_arrows: bool = False  # a small chevron on each loop pointing N -> C
    residue_numbers: bool = False  # first and last residue number at each element's ends
    highlight: str = "none"  # none / asu (as deposited) / protomer (the first) / chains "A,B": the rest in grey
    disulfides: bool = True  # yellow bars between bridged cysteines
    glycans: bool = True  # SNFG symbols on glycosylated residues
    mark: str = ""  # elements in a colour of their own: "G", "A:G=#2ca02c", "#3", "res:150-159=red", comma-separated

    def validate(self) -> Style:
        for key, allowed in CHOICES.items():
            if getattr(self, key) not in allowed:
                raise ValueError(f"{key} must be one of {', '.join(allowed)}; got {getattr(self, key)!r}")
        if not re.fullmatch(r"none|asu|protomer|[A-Za-z0-9]+(,[A-Za-z0-9]+)*", self.highlight):
            raise ValueError(f"highlight must be none, asu, protomer or chains like A,B; got {self.highlight!r}")
        for key in SCALES:
            if not getattr(self, key) > 0:
                raise ValueError(f"{key} must be a positive number; got {getattr(self, key)!r}")
        return self


THEMES = {
    "publication": Style(),
    "minimal": Style(fill="outline", helix_shading="none", loop_width=0.8),
    "print": Style(palette="greys", fill="pale", helix_shading="none"),
    "presentation": Style(loop_width=1.6, font_scale=1.4, helix_shading="gloss", loops="curved"),
    "cartoon": Style(loops="curved", helix_shading="gloss", palette="tol-bright"),
    "rainbow": Style(color_by="sequence", helix_shading="gloss"),
    "shaded": Style(color_by="shade"),
    "richardson": Style(color_by="sstype", helix_shading="depth", loops="curved"),
    "flexibility": Style(color_by="bfactor"),
    "hydropathy": Style(color_by="hydropathy"),
    "conservation": Style(color_by="conservation", loop_color="residue"),
    "alphafold": Style(color_by="plddt", loop_color="residue", loop_width=1.3),
    "goodsell": Style(palette="tol-bright", fill="pale", helix_shading="none", loop_width=1.3, sheet_panels=False),
    "journal": Style(palette="tol-muted", loop_width=0.8, font_scale=0.9),
    "trace": Style(color_by="shade", loop_color="element", loop_arrows=True),
}
PRESETS = THEMES  # earlier name
_TYPES = {f.name: f.type for f in fields(Style)}
THEME_KEYS = tuple(k for k in _TYPES if k not in LAYOUT_KEYS)


def _coerce(key: str, value):
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
        loaded = yaml.safe_load(Path(style_file).read_text()) or {}
        if not isinstance(loaded, dict):
            raise ValueError(f"{style_file}: a style file is a mapping of key: value")
        named = loaded.pop("theme", None) or loaded.pop("preset", None)
        loaded.pop("preset", None)
        theme = theme or named
        settings.update({k: _coerce(k, v) for k, v in loaded.items()})
    name = theme or "publication"
    if name not in THEMES:
        raise ValueError(f"unknown theme {name!r}; choose from {', '.join(THEMES)}")
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"style override {item!r} must be key=value")
        key, value = item.split("=", 1)
        settings[key.strip()] = _coerce(key.strip(), value.strip())
    return replace(THEMES[name], **settings).validate()


def describe(style: Style) -> str:
    return "\n".join(f"{k}: {v}" for k, v in asdict(style).items())
