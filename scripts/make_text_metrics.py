# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Write src/foldmap/text_metrics.py: the advance widths, pair kerning and glyph ink outlines (convex hulls) of
Liberation Sans Bold and Regular (metrically equal to Arial, the figures' first font), so the layout measures its
labels the same way on every machine, whatever fonts are installed.

    python scripts/make_text_metrics.py [/path/to/LiberationSans-Bold.ttf]

The Regular font is read from the same folder (LiberationSans-Regular.ttf).
"""

from __future__ import annotations

import sys
from pathlib import Path

from fontTools.pens.basePen import BasePen
from fontTools.ttLib import TTFont

FONT = Path("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf")
OUT = Path(__file__).resolve().parent.parent / "src" / "foldmap" / "text_metrics.py"
HEADER = (
    "# SPDX-License-Identifier: AGPL-3.0-or-later\n"
    "# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.\n"
)
GREEK = "".join(chr(c) for c in range(0x3B1, 0x3CA)) + "".join(chr(c) for c in range(0x391, 0x3AA) if c != 0x3A2)
CHARS = "".join(chr(c) for c in range(0x20, 0x7F)) + "′″‴·–—‘’“”×±°" + GREEK  # what names are written with
CURVE_STEPS = 16  # points each outline curve is flattened into
HULL_TOLERANCE = 8  # font units a hull vertex may be dropped by: the hull only shrinks (an inner bound on the ink)


def _quote(s: str) -> str:
    """A string literal as ruff format writes it: double quotes unless the text holds one and no single quote."""
    if '"' in s and "'" not in s:
        return repr(s)
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


class _Flatten(BasePen):
    """A glyph's outline as points: its line ends and its curves flattened (points on the curves, never their
    control points, so the hull of the points lies inside the drawn ink)."""

    def __init__(self, glyphs):
        super().__init__(glyphs)
        self.points: list[tuple[float, float]] = []

    def _moveTo(self, p):
        self.points.append(p)

    def _lineTo(self, p):
        self.points.append(p)

    def _qCurveToOne(self, p1, p2):
        p0 = self._getCurrentPoint()
        for k in range(1, CURVE_STEPS + 1):
            t = k / CURVE_STEPS
            u = 1 - t
            self.points.append(tuple(u * u * a + 2 * u * t * b + t * t * c for a, b, c in zip(p0, p1, p2)))

    def _curveToOne(self, p1, p2, p3):
        p0 = self._getCurrentPoint()
        for k in range(1, CURVE_STEPS + 1):
            t = k / CURVE_STEPS
            u = 1 - t
            self.points.append(
                tuple(u**3 * a + 3 * u * u * t * b + 3 * u * t * t * c + t**3 * d for a, b, c, d in zip(p0, p1, p2, p3))
            )


def _hull(points) -> list[tuple[float, float]]:
    """Convex hull, counter-clockwise (monotone chain)."""
    pts = sorted(set(points))
    if len(pts) < 3:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list = []
    upper: list = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _thin(hull, tol: float) -> list[tuple[float, float]]:
    """Drop the hull vertices nearest the chord of their neighbours while within `tol`: a vertex dropped from a
    convex polygon only cuts a corner off, so the thinned hull stays inside the ink."""
    hull = list(hull)
    while len(hull) > 3:
        best, k = None, None
        for i in range(len(hull)):
            a, p, b = hull[i - 1], hull[i], hull[(i + 1) % len(hull)]
            dx, dy = b[0] - a[0], b[1] - a[1]
            dist = abs((p[0] - a[0]) * dy - (p[1] - a[1]) * dx) / max((dx * dx + dy * dy) ** 0.5, 1e-9)
            if best is None or dist < best:
                best, k = dist, i
        if best > tol:
            break
        del hull[k]
    return hull


def _metrics(path: Path):
    font = TTFont(path)
    cmap, hmtx, glyphs = font.getBestCmap(), font["hmtx"], font.getGlyphSet()
    missing = [c for c in CHARS if ord(c) not in cmap]
    if missing:  # never silently: the layout would measure them at the fallback width while render draws another font
        raise SystemExit(f"{path.name} has no glyph for {' '.join(missing)}")
    glyph = {c: cmap[ord(c)] for c in CHARS}
    advance = {c: hmtx[g][0] for c, g in glyph.items()}
    pairs = {}
    for table in font["kern"].kernTables if "kern" in font else []:
        for (a, b), v in sorted(table.kernTable.items()):
            for ca in (c for c, g in glyph.items() if g == a):
                for cb in (c for c, g in glyph.items() if g == b):
                    if v:
                        pairs[ca + cb] = v
    hulls = {}
    for c, g in glyph.items():
        pen = _Flatten(glyphs)
        glyphs[g].draw(pen)
        if pen.points:
            hulls[c] = " ".join(f"{round(x)} {round(y)}" for x, y in _thin(_hull(pen.points), HULL_TOLERANCE))
    os2 = font["OS/2"]
    return font, advance, pairs, hulls, (os2.sTypoAscender, -os2.sTypoDescender)


def _hull_lines(hulls: dict[str, str]) -> list[str]:
    """Each hull as `"c": "x y x y ...",`, split over parenthesised string pieces when longer than a line."""
    out = []
    for c, h in hulls.items():
        line = f"    {_quote(c)}: {_quote(h)},"
        if len(line) <= 120:
            out.append(line)
            continue
        words, pieces, piece = h.split(), [], ""
        for i in range(0, len(words), 2):
            pair = f"{words[i]} {words[i + 1]} "  # each piece ends in the space before the next one's first number
            if piece and len(piece) + len(pair) > 110 < len(h):  # pieces fill their lines: ruff format keeps them
                pieces.append(piece)
                piece = ""
            piece += pair
        pieces.append(piece.rstrip())
        out += [f"    {_quote(c)}: (", *(f"        {_quote(p)}" for p in pieces), "    ),"]
    return out


def main(path: Path) -> None:
    bold, advance, pairs, hulls, (ascent, descent) = _metrics(path)
    regular_path = path.with_name(path.name.replace("Bold", "Regular"))
    _, advance_n, pairs_n, hulls_n, line_n = _metrics(regular_path)
    if line_n != (ascent, descent):
        raise SystemExit(f"{path.name} and {regular_path.name} differ in their typographic ascender and descender")
    units = bold["head"].unitsPerEm
    version = bold["name"].getDebugName(5)
    lines = [
        *HEADER.splitlines(),
        '"""Text metrics of Liberation Sans Bold and Regular (metrically Arial; SIL Open Font License 1.1),',
        f"generated by scripts/make_text_metrics.py from {path.name} and {regular_path.name}",
        f'({version}). Do not edit by hand."""',
        "",
        f"UNITS = {units}  # font units per em",
        f"ASCENT = {ascent}  # OS/2 typographic ascender and descender, both weights: the line a text is centred on",
        f"DESCENT = {descent}",
        "ADVANCE = {  # character: advance width, font units (Bold)",
        *(f"    {_quote(c)}: {v!r}," for c, v in advance.items()),
        "}",
        "KERN = {  # character pair: kerning added between them, font units (Bold)",
        *(f"    {_quote(k)}: {v}," for k, v in pairs.items()),
        "}",
        "HULL = {  # character: its ink's convex hull, x y pairs in font units from its pen origin (Bold)",
        *_hull_lines(hulls),
        "}",
        "ADVANCE_NORMAL = {  # Regular",
        *(f"    {_quote(c)}: {v!r}," for c, v in advance_n.items()),
        "}",
        "KERN_NORMAL = {  # Regular",
        *(f"    {_quote(k)}: {v}," for k, v in pairs_n.items()),
        "}",
        "HULL_NORMAL = {  # Regular",
        *_hull_lines(hulls_n),
        "}",
        "",
    ]
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT}: {len(advance)} characters, {len(pairs)} and {len(pairs_n)} kerning pairs")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else FONT)
