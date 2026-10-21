"""foldmap compare (plan 13a-8): A on the left, B drawn like A on the right, at one scale, with one legend and the note
under both; what was matched as JSON; the MCP tool."""

import asyncio
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from matplotlib.collections import Collection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text

from foldmap import cli
from foldmap.compare import compare, write_json

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).resolve().parent.parent
UBQ, LZM, HHB = DATA / "1UBQ.cif", DATA / "2LZM.cif", DATA / "2HHB.cif"


def _bench():
    spec = importlib.util.spec_from_file_location("foldmap_bench", ROOT / "scripts" / "bench.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BENCH = _bench()


def panel_axes(fig):
    """The two panels' axes (A's, B's): the page's own axes (legend, note) comes first."""
    base, a, b = fig.axes
    return a, b


def ink(ax, prefix: str) -> list:
    """Every artist of a panel with a gid: (gid less the panel's prefix, type, where it is in the panel's own page
    units), sorted. Every gid there carries the prefix."""
    out = []
    for artist in ax.findobj(lambda o: o is not ax and o.get_gid() is not None):
        gid = artist.get_gid()
        assert gid.startswith(prefix), gid
        if isinstance(artist, Text):
            where = (artist.get_text(), tuple(np.round(artist.get_position(), 9)), artist.get_fontsize())
        elif isinstance(artist, Patch):
            where = np.round(artist.get_patch_transform().transform(artist.get_path().vertices), 9).tobytes()
        elif isinstance(artist, Line2D):
            where = np.round(artist.get_xydata(), 9).tobytes()
        elif isinstance(artist, Collection):
            where = (np.round(artist.get_offsets(), 9).tobytes(), len(artist.get_paths()))
        else:
            where = None
        out.append((gid[len(prefix) :], type(artist).__name__, where))
    return sorted(out, key=repr)


def test_itself_gives_two_identical_panels(tmp_path):
    """1UBQ compared with itself: the panels carry the same gids after their prefixes, every one at the same place
    in its panel, and on the page B's sits exactly the panel offset right of A's; the page holds the note."""
    made = compare(UBQ, UBQ)
    a, b = panel_axes(made.figure)
    ink_a, ink_b = ink(a, "a/"), ink(b, "b/")
    assert len(ink_a) > 20 and ink_a == ink_b
    pa, pb = a.get_position(), b.get_position()
    assert pa.width == pytest.approx(pb.width) and pa.height == pytest.approx(pb.height) and pa.y0 == pb.y0
    width = made.figure.get_size_inches()[0] * made.figure.dpi
    offset = (pb.x0 - pa.x0) * width
    assert offset > 0
    texts_a = [t for t in a.texts if t.get_gid()]
    texts_b = {t.get_gid()[2:]: t for t in b.texts if t.get_gid()}
    for t in texts_a:  # on the page: B's text is A's moved right by the panel offset
        u = texts_b[t.get_gid()[2:]]
        xa = t.get_transform().transform(t.get_position())
        xb = u.get_transform().transform(u.get_position())
        assert np.allclose(xb - xa, (offset, 0.0), atol=1e-6)
    assert made.note == (
        "drawn like 1UBQ.cif chain A: 8 of 8 elements matched (sequence alignment, TM-score 1.00 over 76 residues)"
    )
    (note,) = [t for t in made.figure.axes[0].texts if t.get_gid() == "compare-note"]
    assert " ".join(note.get_text().split("\n")) == made.note
    out = tmp_path / "pair.svg"
    cli.save(made.figure, out)
    ids = re.findall(r'id="([ab])/([^"]+)"', out.read_text(encoding="utf-8"))
    assert ids and sorted(g for p, g in ids if p == "a") == sorted(g for p, g in ids if p == "b")
    assert {"strand", "helix", "label"} <= {g.split(":")[0] for _, g in ids}


def test_one_scale_and_one_legend():
    """The panels share the larger page's scale (a page unit as long in both) and the page has one legend: β like α
    coloured by type has the type entries once, under both panels, none in a panel."""
    from foldmap.style import Style

    made = compare(HHB, HHB, "A", "B", look=Style(color_by="sstype"))
    a, b = panel_axes(made.figure)

    def per_unit(ax):  # inches per page unit across the panel
        x0, x1 = ax.get_xlim()
        return ax.get_position().width * made.figure.get_size_inches()[0] / (x1 - x0)

    assert per_unit(a) == pytest.approx(per_unit(b), rel=1e-9)
    base = made.figure.axes[0]
    words = [t.get_text() for t in base.texts if t.get_gid() != "compare-note"]
    assert "helix" in words and "3₁₀ helix" in words and len(words) == len(set(words))
    assert not [t for ax in (a, b) for t in ax.texts if t.get_text() in ("helix", "strand", "3₁₀ helix")]


def test_json_schema_and_order(tmp_path):
    """--json: {alignment: {method, tm_score, rmsd, identity, pairs}, mates, pieces, extras, only_a}, numbers to 1e-3,
    lists in sequence order, keys sorted in the file; β (2HHB B) like α (A): its mates take α's names, its two
    extras are named as drawn."""
    made = compare(HHB, HHB, "A", "B")
    data = made.data
    assert set(data) == {"alignment", "mates", "pieces", "extras", "only_a"}
    aln = data["alignment"]
    assert set(aln) == {"method", "tm_score", "rmsd", "identity", "pairs"}
    assert aln["method"] == "sequence" and round(aln["tm_score"], 2) == 0.89
    for key in ("tm_score", "rmsd", "identity"):
        assert aln[key] == round(aln[key], 3)
    assert len(aln["pairs"]) == 139
    firsts, seconds = ([int(re.match(r"-?\d+", p[k]).group()) for p in aln["pairs"]] for k in (0, 1))
    assert firsts == sorted(firsts) and seconds == sorted(seconds)

    def start(ref):  # res:B:3-17 -> 3
        return int(ref.split(":")[2].split("-")[0])

    assert len(data["mates"]) + len(data["pieces"]) == 7
    for m in data["mates"]:
        assert set(m) == {"a", "a_residues", "b", "b_residues"} and m["a"] == m["b"]
        assert m["a_residues"].startswith("res:A:") and m["b_residues"].startswith("res:B:")
    assert [start(m["b_residues"]) for m in data["mates"]] == sorted(start(m["b_residues"]) for m in data["mates"])
    assert [e["b"] for e in data["extras"]] == ["η1", "α3a"]
    lay_b = made.panels[1][0]
    assert [lay_b.placed[k].label for k in lay_b.like_extras] == ["η1", "α3a"]
    out = write_json(data, tmp_path / "mates.json")

    def keys_sorted(pairs):
        keys = [k for k, _ in pairs]
        assert keys == sorted(keys)
        return dict(pairs)

    assert json.loads(out.read_text(encoding="utf-8"), object_pairs_hook=keys_sorted) == data


def test_negative_pair_draws_both_panels_on_their_own():
    """T4 lysozyme (2LZM) is not ubiquitin's fold: both panels are each structure drawn on its own (bench
    signature), the note says why under both, and the JSON still says what was aligned."""
    made = compare(UBQ, LZM)
    (lay_a, _, _), (lay_b, _, _) = made.panels
    alone_a, _, _ = cli.make_layout(UBQ, chains=["A"], **BENCH.OPTIONS)
    alone_b, _, _ = cli.make_layout(LZM, chains=["A"], **BENCH.OPTIONS)
    assert BENCH.signature(lay_a) == BENCH.signature(alone_a)
    assert BENCH.signature(lay_b) == BENCH.signature(alone_b)
    assert made.note == "not drawn like 1UBQ.cif: not the same fold (TM-score 0.14; 4 of 13 elements matched)"
    (note,) = [t for t in made.figure.axes[0].texts if t.get_gid() == "compare-note"]
    assert " ".join(note.get_text().split("\n")) == made.note
    assert round(made.data["alignment"]["tm_score"], 2) == 0.14
    assert len(made.data["mates"]) + len(made.data["pieces"]) == 4
    assert any(i["kind"] == "like" and i["text"].startswith("right: not drawn like") for i in made.provenance)


def test_command_line(tmp_path, capsys):
    """foldmap compare A B -o pair.svg --json mates.json: both written, the note on stderr; a chain A lacks is an
    error naming it."""
    svg, js = tmp_path / "pair.svg", tmp_path / "mates.json"
    assert cli.main(["compare", str(UBQ), str(UBQ), "-o", str(svg), "--json", str(js)]) == 0
    err = capsys.readouterr().err
    assert "right: drawn like 1UBQ.cif chain A: 8 of 8 elements matched" in err
    assert 'id="a/' in svg.read_text(encoding="utf-8") and json.loads(js.read_text(encoding="utf-8"))["mates"]
    assert cli.main(["compare", str(UBQ), str(UBQ), "--chain-a", "Z", "-o", str(svg)]) == 1
    assert "no protein chain Z" in capsys.readouterr().err


_PAGE = """
import hashlib, sys
from foldmap import cli
svg, js = sys.argv[2], sys.argv[3]
assert cli.main(["compare", sys.argv[1], sys.argv[1], "--chain-a", "A", "--chain-b", "B", "-o", svg, "--json", js]) == 0
print(hashlib.sha256(open(svg, "rb").read() + open(js, "rb").read()).hexdigest())
"""


def test_page_does_not_depend_on_the_hash_seed(tmp_path):
    """A12: the compare page and its JSON have the same bytes under PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _PAGE, str(HHB), str(tmp_path / f"{seed}.svg"), str(tmp_path / f"{seed}.json")],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout.split()[-1])
    assert len(runs[0]) == 64 and runs[0] == runs[1]


def test_mcp_returns_the_preview_and_the_json(tmp_path, monkeypatch):
    """compare_structures writes the page and the JSON and returns the note, a PNG preview and what was matched."""
    pytest.importorskip("mcp")
    from foldmap.mcp_server import server

    monkeypatch.delenv("FOLDMAP_OUTPUT_DIR", raising=False)
    out, js = tmp_path / "pair.png", tmp_path / "mates.json"
    r = asyncio.run(
        server.call_tool(
            "compare_structures",
            {"source_a": str(UBQ), "source_b": str(UBQ), "output_path": str(out), "json_path": str(js)},
        )
    )
    assert [c.type for c in r.content] == ["text", "image"]
    assert out.is_file() and str(out) in r.content[0].text and "drawn like 1UBQ.cif chain A" in r.content[0].text
    data = r.structured_content
    assert data["alignment"]["method"] == "sequence" and len(data["mates"]) == 8
    assert any(i["kind"] == "like" for i in data["provenance"])
    saved = json.loads(js.read_text(encoding="utf-8"))
    assert saved == {k: v for k, v in data.items() if k != "provenance"}
