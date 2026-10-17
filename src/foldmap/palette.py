# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Chain colours: named colour-blind-safe sets (Okabe-Ito, Paul Tol), each in a fixed order."""

from __future__ import annotations

import threading

PALETTES = {
    "okabe-ito": ["#0072B2", "#E69F00", "#CC79A7", "#009E73", "#56B4E9", "#D55E00", "#F0E442", "#999999"],
    "tol-bright": ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"],
    "tol-muted": ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"],
    "blueprint": ["#9fd3ff", "#ffd27f", "#b6f0b1", "#ffb3cf", "#d7c4ff", "#8ff0e6", "#ffffff"],  # pale, for dark pages
    "greys": ["#404040", "#8C8C8C", "#BFBFBF", "#262626", "#A6A6A6", "#666666"],
}


BACK_SHADE = 0.72  # a helix's back face is its colour darkened this much (render)
CHAIN_DE = 10.0  # CIEDE2000: a chain's colour past the palette's first pass stands this far from every other
_PASS_STEPS = (("tint", 0.45), ("darken", 0.75), ("tint", 0.7), ("darken", 0.55), ("tint", 0.25), ("darken", 0.85))
_LAB_STEPS = tuple(range(-60, 61, 3))  # L* steps tried when none of _PASS_STEPS stands CHAIN_DE apart


def chain_colors(chains: list[str], palette: str = "okabe-ito") -> dict[str, str]:
    """One colour per chain, assigned in the order given: the palette's colours, then (random-100 C35: chains A and I
    drew one #0072B2) each further pass through it a lightness step of the colour it repeats: the pass's own step of
    _PASS_STEPS first (tint 0.45, then darken 0.75, ...), the first in that order that stands CHAIN_DE (CIEDE2000)
    from every colour already assigned, from the white page and from the helix back faces (BACK_SHADE) of the colours
    it repeats; else the step of its L* (_LAB_STEPS) that stands furthest from them. A palette with few hues (greys)
    cannot keep CHAIN_DE past its first passes; its colours still differ. Deterministic: the order given, ties by
    the colour's hex. A chain's colour depends only on its place in the order, so each palette's colours are worked
    out once in this process and kept (_KEPT; final review C: 3.7 s for 240 chains, asked two or three times a
    figure)."""
    if palette not in PALETTES:
        raise ValueError(f"unknown palette {palette!r}; choose from {', '.join(PALETTES)}")
    with _KEPT_LOCK:
        out = _KEPT.setdefault(palette, [])
        _more_colours(out, PALETTES[palette], len(chains))
        return dict(zip(chains, out))


_KEPT: dict[str, list[str]] = {}  # palette -> its chain colours in order, as many as asked for so far
_KEPT_LOCK = threading.Lock()  # figures are drawn in threads (the app's downloads)


def _more_colours(out: list[str], cols: list[str], count: int) -> None:
    """Extend `out`, the first chain colours of the palette `cols` (chain_colors), to `count` colours."""
    n = len(cols)
    labs: dict[str, tuple[float, float, float]] = {}

    def lab(colour: str) -> tuple[float, float, float]:  # _lab, each colour once
        if colour not in labs:
            labs[colour] = _lab(colour)
        return labs[colour]

    for i in range(len(out), count):
        base = cols[i % n]
        if i < n:
            out.append(base)
            continue
        family = out[i % n :: n]  # the colours of the earlier passes over this palette colour
        refs = [*out, "#ffffff", *(darken(c, BACK_SHADE) for c in family)]

        def apart(colour: str, refs=refs) -> float:  # delta_e, its colours' L*a*b* kept
            return min(delta_e_lab(lab(colour), lab(r)) for r in refs)

        k = (i // n - 1) % len(_PASS_STEPS)
        steps = _PASS_STEPS[k:] + _PASS_STEPS[:k]
        made = [tint(base, a) if how == "tint" else darken(base, a) for how, a in steps]
        pick = next((c for c in made if apart(c) >= CHAIN_DE), None)
        if pick is None:
            at = _lab(base)
            ladder = [_from_lab(at[0] + d, at[1], at[2]) for d in _LAB_STEPS if 8 <= at[0] + d <= 92]
            pick = max(ladder, key=lambda c: (apart(c), c))
        out.append(pick)


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


def _from_lab(L: float, a: float, b: float) -> str:  # noqa: N803 (CIE's own letter)
    """The sRGB hex of a CIE L*a*b* (D65) colour, clipped into the gamut (the inverse of _lab)."""
    fy = (L + 16) / 116
    fx, fz = fy + a / 500, fy - b / 200

    def finv(t: float) -> float:
        return t**3 if t**3 > 216 / 24389 else (116 * t - 16) / (24389 / 27)

    x, y, z = finv(fx) * 0.95047, finv(fy), finv(fz) * 1.08883
    rgb = (
        3.2404542 * x - 1.5371385 * y - 0.4985314 * z,
        -0.9692660 * x + 1.8760108 * y + 0.0415560 * z,
        0.0556434 * x - 0.2040259 * y + 1.0572252 * z,
    )

    def gamma(c: float) -> float:
        c = min(max(c, 0.0), 1.0)
        return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055

    return "#" + "".join(f"{round(gamma(c) * 255):02x}" for c in rgb)


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
