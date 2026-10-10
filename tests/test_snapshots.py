"""Figures compared with recorded structural fingerprints (ids, text, rounded extents), not pixels.

Refactors of the layout and drawing code must keep these unchanged. After an intended visual change:
    python -m pytest tests/test_snapshots.py --update-snapshots
and review the JSON diff before committing it.
"""

import sys
from pathlib import Path

import pytest
from matplotlib import font_manager

from foldmap import render
from foldmap.cli import main
from snapshot import compare, dump, fingerprint, load

DATA = Path(__file__).parent / "data"
SNAPSHOTS = DATA / "snapshots"

SAMPLE = """<svg xmlns="http://www.w3.org/2000/svg" width="100pt" height="50.25pt">
<g id="figure_1"><g id="patch_1"><path d="M 0 0 L 100 50 z"/></g>
<g id="helix:A:2-9"><path d="M 10 20 L 30.04 40 z" style="fill: #ABC; stroke: rgb(0, 0, 255)"/></g>
<g id="label:A:2-9"><text x="12.5" y="8">α1</text></g>
<g id="loop:A:2-9>A:12-15"><path d="M 1 1 L 2 2"/></g>
<g id="loop:A:2-9>A:12-15"><path d="M 3 3 L 4 4"/></g>
<g id="text_7"><text x="1" y="1">Chain A</text></g></g></svg>"""


def test_fingerprint_keeps_named_items_text_and_extent():
    fp = fingerprint(SAMPLE)
    assert fp["size"] == [100.0, 50.2]
    assert "figure_1" not in fp["items"] and "patch_1" not in fp["items"]
    assert fp["items"]["helix:A:2-9"] == {
        "bbox": [10.0, 20.0, 30.0, 40.0],
        "text": [],
        "colours": ["fill:#aabbcc", "stroke:#0000ff"],
    }
    assert fp["items"]["label:A:2-9"]["text"] == ["α1"]
    assert {"loop:A:2-9>A:12-15", "loop:A:2-9>A:12-15#2"} <= set(fp["items"])
    assert fp["unnamed"] == {"patch": 1, "text": 1}


def test_compare_reports_moves_beyond_tolerance_only():
    a = fingerprint(SAMPLE)
    b = fingerprint(SAMPLE.replace("30.04 40", "30.5 40"))
    c = fingerprint(SAMPLE.replace("30.04 40", "32 40"))
    assert compare(a, b) == []
    assert any(p.startswith("helix:A:2-9: moved") for p in compare(a, c))
    d = fingerprint(SAMPLE.replace(">α1<", ">α2<"))
    assert any("label:A:2-9: text" in p for p in compare(a, d))
    e = fingerprint(SAMPLE.replace('id="helix:A:2-9"', 'id="helix:A:2-10"'))
    assert any("missing helix:A:2-9" in p for p in compare(a, e))
    assert any("new helix:A:2-10" in p for p in compare(a, e))


HASHED = """<svg xmlns="http://www.w3.org/2000/svg" width="10pt" height="10pt">
<defs><clipPath id="p1a2b3c4d5e"><rect x="0" y="0" width="5" height="5"/></clipPath>
<path id="m0123456789" d="M 0 0 L 1 1"/><pattern id="hfedcba9876"/></defs>
<g id="helix:A:1-5" clip-path="url(#p1a2b3c4d5e)"><path d="M 1 1 L 2 2" style="fill: url(#hfedcba9876)"/></g></svg>"""


def test_hashed_matplotlib_ids_are_ignored():
    fp = fingerprint(HASHED)
    assert list(fp["items"]) == ["helix:A:1-5"]
    renamed = HASHED.replace("p1a2b3c4d5e", "p00000000ff").replace("m0123456789", "maaaaaaaaaa")
    assert compare(fp, fingerprint(renamed)) == []
    assert fp["unnamed"] == {}


def test_compare_detects_a_colour_only_change():
    a = fingerprint(SAMPLE)
    b = fingerprint(SAMPLE.replace("#ABC", "#ABD"))
    assert compare(a, fingerprint(SAMPLE.replace("#ABC", "#aabbcc"))) == []
    assert any(p.startswith("helix:A:2-9: colours") for p in compare(a, b))
    assert fingerprint(SAMPLE.replace("#ABC", "none"))["items"]["helix:A:2-9"]["colours"][0] == "fill:none"


def test_snapshot_roundtrip_keeps_unicode_and_lf(tmp_path):
    fp = fingerprint(SAMPLE.replace(">α1<", ">η1 A′<"))
    path = tmp_path / "s.json"
    dump(path, fp)
    raw = path.read_bytes()
    assert "η1 A′".encode() in raw and b"\r\n" not in raw
    assert load(path) == fp


CASES = {
    "1UBQ-publication": ("1UBQ.cif", []),  # sheet panel, helix, terminus
    "1TIM-symmetry": ("1TIM.cif", []),  # symmetric assembly: protomer ghosts
    "1LMB-dna": ("1LMB.cif", []),  # DNA ladder, contacts, protein-DNA
    "1C3W-membrane": ("1C3W.cif", ["--membrane", "auto"]),  # membrane slab and heads, ligands
    "6ZS5-disulfide-glycan": ("6ZS5.cif", []),  # disulfide lane brackets, glycans, ss-dots
    "2HHB-journal": ("2HHB.cif", ["--theme", "journal"]),  # journal colours and fonts
    "2HHB-shaded": ("2HHB.cif", ["--theme", "shaded"]),  # shade branch of the colouring
    "1TIM-shaded": ("1TIM.cif", ["--theme", "shaded"]),  # shade branch with symmetry copies
    "2HHB-trace": ("2HHB.cif", ["--theme", "trace"]),  # trace theme
    "2HHB-highlight": ("2HHB.cif", ["--set", "highlight=A"]),  # dimmed MATE_GREY branch, other chains grey
    "2LZM-domains-rainbow": ("2LZM.cif", ["--domains", "auto", "--theme", "rainbow", "--set", "residue_numbers=true"]),
}


def _resolved_font_family() -> str:
    """Family of the font matplotlib picks for foldmap's font list (what the figures are measured in)."""
    try:
        path = font_manager.findfont(font_manager.FontProperties(family=render._FONTS))
        return font_manager.FontProperties(fname=path).get_name()
    except Exception:  # noqa: BLE001 - an unresolvable font simply disables the snapshots
        return ""


_FONT = _resolved_font_family() if sys.platform == "linux" else ""
_SNAPSHOT_SKIP = (
    "snapshots need Linux with Liberation Sans installed (CI: apt-get install fonts-liberation): "
    f"fonts differ by machine, so text metrics and glyph paths do; matplotlib resolved {_FONT!r}"
)


@pytest.mark.skipif(sys.platform != "linux" or _FONT != "Liberation Sans", reason=_SNAPSHOT_SKIP)
@pytest.mark.parametrize("name", sorted(CASES))
def test_figure_matches_snapshot(name, tmp_path, update_snapshots):
    file, flags = CASES[name]
    out = tmp_path / "f.svg"
    assert main(["plot", str(DATA / file), "-o", str(out), *flags]) == 0
    actual = fingerprint(out.read_text(encoding="utf-8"))
    path = SNAPSHOTS / f"{name}.json"
    if update_snapshots:
        dump(path, actual)
        pytest.skip(f"wrote {path.name}")
    if not path.exists():
        pytest.fail(f"no snapshot {path.name}: run  python -m pytest tests/test_snapshots.py --update-snapshots")
    problems = compare(load(path), actual)
    assert not problems, (
        f"{name} differs from its snapshot:\n  "
        + "\n  ".join(problems[:40])
        + "\nIf this change is intended: python -m pytest tests/test_snapshots.py --update-snapshots"
    )
