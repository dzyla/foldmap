"""The linked interactive page: topology SVG + contact map + 3D model, one self-contained HTML file."""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from foldmap.interactive import THREEDMOL, build_page

DATA = Path(__file__).parent / "data"


def _data(page: str) -> dict:
    blob = re.search(r'<script type="application/json" id="topo-data">(.*?)</script>', page, re.S).group(1)
    return json.loads(blob)


def _script(page: str, sid: str) -> str:
    return re.search(rf'<script id="{sid}">(.*?)</script>', page, re.S).group(1)


@pytest.mark.parametrize("name", ["1UBQ", "1LMB", "5NKT"])
def test_page_carries_svg_data_and_model(name):
    page = build_page(DATA / f"{name}.cif")
    data = _data(page)
    assert "<svg" in page and THREEDMOL in page
    assert len(data["residues"]) == len(data["ca"]) > 0
    ids = set(re.findall(r'id="([^"]+)"', page))
    for e in data["elements"]:
        assert e["id"] in {i.split(":", 1)[1] for i in ids if i.startswith(("strand:", "helix:", "eta:"))}
        assert 0 <= e["start"] <= e["end"] < len(data["residues"]) and re.fullmatch(r"#[0-9a-f]{6}", e["colour"])
    assert "_atom_site" in data["model"]
    assert all(r["e"] is None or r["e"] in {e["id"] for e in data["elements"]} for r in data["residues"])


def test_residues_point_back_to_their_elements():
    data = _data(build_page(DATA / "1UBQ.cif"))
    for e in data["elements"]:
        assert all(data["residues"][k]["e"] == e["id"] for k in range(e["start"], e["end"] + 1))


def test_page_honours_the_theme():
    from foldmap.style import Style

    plain = _data(build_page(DATA / "1UBQ.cif"))
    rainbow = _data(build_page(DATA / "1UBQ.cif", look=Style(color_by="sequence")))
    assert [e["colour"] for e in plain["elements"]] != [e["colour"] for e in rainbow["elements"]]


node = shutil.which("node")


@pytest.mark.skipif(node is None, reason="node not installed")
def test_page_scripts_parse(tmp_path):
    page = build_page(DATA / "1UBQ.cif")
    for sid in ("topo-lib", "topo-app"):
        f = tmp_path / f"{sid}.js"
        f.write_text(_script(page, sid))
        assert subprocess.run([node, "--check", str(f)], capture_output=True).returncode == 0, sid


@pytest.mark.skipif(node is None, reason="node not installed")
def test_library_logic_in_node(tmp_path):
    page = build_page(DATA / "1UBQ.cif")
    data = _data(page)
    probe = tmp_path / "probe.js"
    probe.write_text(_script(page, "topo-lib") + """
const data = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const lib = globalThis.TopoLib;
const idx = lib.index(data);
const out = {
  bins: lib.binning(data.ca.length, 50),
  dist: lib.distance(data.ca, 0, 1),
  el: data.elements.map(e => lib.elementAt(data, e.start)),
  key: idx.byKey[data.residues[5].c + ':' + data.residues[5].n],
  span: lib.span(data, data.elements[0].id),
  img: Array.from(lib.matrix(data.ca, 20).values.slice(0, 3)).map(v => Math.round(v * 10) / 10),
};
console.log(JSON.stringify(out));
""")
    blob = tmp_path / "data.json"
    blob.write_text(json.dumps(data))
    res = subprocess.run([node, str(probe), str(blob)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    out = json.loads(res.stdout)
    assert out["bins"] == {"size": 2, "count": 38}  # 76 residues into at most 50 cells
    assert abs(out["dist"] - 3.8) < 0.3  # neighbouring CA atoms
    assert out["el"] == [e["id"] for e in data["elements"]]
    assert out["key"] == 5
    assert out["span"] == [data["elements"][0]["start"], data["elements"][0]["end"]]
    assert out["img"][0] == 0.0  # the diagonal


def test_cli_interactive_writes_the_page(tmp_path):
    from foldmap.cli import main

    out = tmp_path / "view.html"
    assert main(["interactive", str(DATA / "1LMB.cif"), "-o", str(out), "--theme", "trace"]) == 0
    assert out.stat().st_size > 10_000 and "topo-data" in out.read_text()
