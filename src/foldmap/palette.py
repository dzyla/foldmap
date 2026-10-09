"""Chain colours: named colour-blind-safe sets (Okabe-Ito, Paul Tol), each in a fixed order."""

from __future__ import annotations

PALETTES = {
    "okabe-ito": ["#0072B2", "#E69F00", "#CC79A7", "#009E73", "#56B4E9", "#D55E00", "#F0E442", "#999999"],
    "tol-bright": ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"],
    "tol-muted": ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"],
    "blueprint": ["#9fd3ff", "#ffd27f", "#b6f0b1", "#ffb3cf", "#d7c4ff", "#8ff0e6", "#ffffff"],  # pale, for dark pages
    "greys": ["#404040", "#8C8C8C", "#BFBFBF", "#262626", "#A6A6A6", "#666666"],
}


def chain_colors(chains: list[str], palette: str = "okabe-ito") -> dict[str, str]:
    """One colour per chain, assigned in the order given (cycles when the palette runs out)."""
    if palette not in PALETTES:
        raise ValueError(f"unknown palette {palette!r}; choose from {', '.join(PALETTES)}")
    cols = PALETTES[palette]
    return {c: cols[i % len(cols)] for i, c in enumerate(chains)}


def _rgb(hex_color: str) -> tuple[float, float, float]:
    return tuple(int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5))  # type: ignore[return-value]


def darken(hex_color: str, factor: float = 0.6) -> str:
    return "#" + "".join(f"{round(c * factor * 255):02x}" for c in _rgb(hex_color))


def tint(hex_color: str, amount: float = 0.6) -> str:
    """Mix with white; amount=1 gives white."""
    return "#" + "".join(f"{round((c + (1 - c) * amount) * 255):02x}" for c in _rgb(hex_color))


def text_color_on(hex_color: str) -> str:
    r, g, b = _rgb(hex_color)
    return "black" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.6 else "white"
