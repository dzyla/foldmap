# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
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


def _lab(colour: str) -> tuple[float, float, float]:
    """CIE L*a*b* (D65) of a colour matplotlib can read."""
    from matplotlib.colors import to_rgb

    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(colour)]
    x = (0.4124564 * lin[0] + 0.3575761 * lin[1] + 0.1804375 * lin[2]) / 0.95047
    y = 0.2126729 * lin[0] + 0.7151522 * lin[1] + 0.0721750 * lin[2]
    z = (0.0193339 * lin[0] + 0.1191920 * lin[1] + 0.9503041 * lin[2]) / 1.08883

    def f(t: float) -> float:
        return t ** (1 / 3) if t > 216 / 24389 else (24389 / 27 * t + 16) / 116

    fx, fy, fz = f(x), f(y), f(z)
    return 116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)


def delta_e(a: str, b: str) -> float:
    """CIEDE2000 colour difference between two colours (Sharma, Wu, Dalal 2005)."""
    return delta_e_lab(_lab(a), _lab(b))


def delta_e_lab(lab1: tuple[float, float, float], lab2: tuple[float, float, float]) -> float:
    """CIEDE2000 between two L*a*b* colours."""
    import math

    l1, a1, b1 = lab1
    l2, a2, b2 = lab2
    c1, c2 = math.hypot(a1, b1), math.hypot(a2, b2)
    g = 0.5 * (1 - math.sqrt(((c1 + c2) / 2) ** 7 / (((c1 + c2) / 2) ** 7 + 25**7)))
    a1p, a2p = (1 + g) * a1, (1 + g) * a2
    c1p, c2p = math.hypot(a1p, b1), math.hypot(a2p, b2)
    h1p = math.degrees(math.atan2(b1, a1p)) % 360 if c1p else 0.0
    h2p = math.degrees(math.atan2(b2, a2p)) % 360 if c2p else 0.0
    dl, dc = l2 - l1, c2p - c1p
    if not c1p * c2p:
        dh = 0.0
    elif abs(h2p - h1p) <= 180:
        dh = h2p - h1p
    else:
        dh = h2p - h1p + (360 if h2p <= h1p else -360)
    dhh = 2 * math.sqrt(c1p * c2p) * math.sin(math.radians(dh) / 2)
    lm, cm = (l1 + l2) / 2, (c1p + c2p) / 2
    if not c1p * c2p:
        hm = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hm = (h1p + h2p) / 2
    else:
        hm = (h1p + h2p + (360 if h1p + h2p < 360 else -360)) / 2
    t = (
        1
        - 0.17 * math.cos(math.radians(hm - 30))
        + 0.24 * math.cos(math.radians(2 * hm))
        + 0.32 * math.cos(math.radians(3 * hm + 6))
        - 0.20 * math.cos(math.radians(4 * hm - 63))
    )
    sl = 1 + 0.015 * (lm - 50) ** 2 / math.sqrt(20 + (lm - 50) ** 2)
    sc, sh = 1 + 0.045 * cm, 1 + 0.015 * cm * t
    rc = 2 * math.sqrt(cm**7 / (cm**7 + 25**7))
    rt = -math.sin(math.radians(60 * math.exp(-(((hm - 275) / 25) ** 2)))) * rc
    return math.sqrt((dl / sl) ** 2 + (dc / sc) ** 2 + (dhh / sh) ** 2 + rt * (dc / sc) * (dhh / sh))
