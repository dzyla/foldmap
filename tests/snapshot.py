"""Structural fingerprints of SVG figures: every named item (gid) with its text and rounded extent, plus counts
of unnamed artists. Stable across runs; tolerant of sub-point numeric noise; blind to pixels and fonts.
Colours (fill and stroke, own and descendants) are compared exactly. `transform` attributes are not applied
(matplotlib's SVG output does not use them for artists)."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

SVG = "{http://www.w3.org/2000/svg}"
TOL = 0.75  # points
_AUTO = re.compile(r"^(figure|axes|patch|line2d|text|PathCollection|PolyCollection|LineCollection|image|legend)_\d+$")
# matplotlib names clip paths, markers and hatches p/m/h + a salted hash of their full-precision geometry,
# so a 1-ulp change renames them; they carry no layout information of their own.
_HASHED = re.compile(r"^[pmh][0-9a-f]{10}$")
_NUM = re.compile(r"-?\d+(?:\.\d+)?(?:e-?\d+)?")


def _bbox(el) -> list[float] | None:
    xs, ys = [], []
    for sub in el.iter():
        if sub.tag == SVG + "path" and sub.get("d"):
            nums = [float(v) for v in _NUM.findall(sub.get("d"))]
            xs += nums[0::2]
            ys += nums[1::2]
        for a, b in (("x", "y"), ("cx", "cy")):
            if sub.get(a) is not None and sub.get(b) is not None:
                xs.append(float(sub.get(a)))
                ys.append(float(sub.get(b)))
    if not xs:
        return None
    return [round(min(xs), 1), round(min(ys), 1), round(max(xs), 1), round(max(ys), 1)]


def _hex(value: str) -> str:
    v = value.strip().lower()
    if v == "none":
        return v
    if re.fullmatch(r"#[0-9a-f]{3}", v):
        return "#" + "".join(c * 2 for c in v[1:])
    m = re.fullmatch(r"rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", v)
    if m:
        return "#%02x%02x%02x" % tuple(int(g) for g in m.groups())
    if v.startswith("#"):
        return v
    try:
        from matplotlib.colors import to_hex

        return to_hex(v)
    except ValueError:
        return v


def _colours(el) -> list[str]:
    found = set()
    for sub in el.iter():
        decls = dict(d.split(":", 1) for d in (sub.get("style") or "").split(";") if ":" in d)
        for prop in ("fill", "stroke"):
            val = decls.get(prop, decls.get(" " + prop, sub.get(prop)))
            if val is not None and val.strip() and not val.strip().startswith("url("):
                found.add(f"{prop}:{_hex(val)}")
    return sorted(found)


def fingerprint(svg_text: str) -> dict:
    root = ET.fromstring(svg_text)
    size = [round(float(root.get(k, "0").removesuffix("pt")), 1) for k in ("width", "height")]
    items: dict[str, dict] = {}
    unnamed: Counter[str] = Counter()
    for el in root.iter():
        gid = el.get("id")
        if not gid:
            continue
        if _HASHED.match(gid):
            continue
        auto = _AUTO.match(gid)
        if auto:
            if auto.group(1) not in ("figure", "axes"):
                unnamed[auto.group(1)] += 1
            continue
        key, n = gid, 1
        while key in items:
            n += 1
            key = f"{gid}#{n}"
        items[key] = {
            "bbox": _bbox(el),
            "text": ["".join(t.itertext()) for t in el.iter(SVG + "text")],
            "colours": _colours(el),
        }
    return {"size": size, "items": items, "unnamed": dict(sorted(unnamed.items()))}


def compare(expected: dict, actual: dict, tol: float = TOL) -> list[str]:
    """Human-readable differences; empty when the figures match."""
    out = []
    if any(abs(a - b) > tol for a, b in zip(expected["size"], actual["size"])):
        out.append(f"page size {expected['size']} -> {actual['size']}")
    if expected["unnamed"] != actual["unnamed"]:
        out.append(f"unnamed artists {expected['unnamed']} -> {actual['unnamed']}")
    old, new = expected["items"], actual["items"]
    out += [f"missing {k}" for k in sorted(set(old) - set(new))]
    out += [f"new {k}" for k in sorted(set(new) - set(old))]
    for k in sorted(set(old) & set(new)):
        e, a = old[k], new[k]
        if e["text"] != a["text"]:
            out.append(f"{k}: text {e['text']} -> {a['text']}")
        if e["colours"] != a["colours"]:
            out.append(f"{k}: colours {e['colours']} -> {a['colours']}")
        if (e["bbox"] is None) != (a["bbox"] is None):
            out.append(f"{k}: extent {e['bbox']} -> {a['bbox']}")
        elif e["bbox"] and max(abs(x - y) for x, y in zip(e["bbox"], a["bbox"])) > tol:
            out.append(f"{k}: moved {e['bbox']} -> {a['bbox']}")
    return out


def dump(path, fp: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(fp, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def load(path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))
