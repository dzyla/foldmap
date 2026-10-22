"""foldmap data (plan 13d-1): the versioned dict of a laid-out structure. Its schema, its labels and sheet orders
against the figure drawn from the same spec, a TIM barrel's n and S, and the same bytes for any hash seed (A12)."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main
from foldmap.datafile import ALL, SCHEMA_VERSION, SECTIONS, dumps, schema, structure_data
from foldmap.spec import draw_spec, figure_spec

DATA = Path(__file__).parent / "data"


def _conforms(value, sch: dict, where: str = "$") -> list[str]:
    """Where `value` breaks `sch` (the JSON Schema subset schema() uses: type, properties, required,
    additionalProperties false, items, enum, const)."""
    kinds = sch.get("type")
    if kinds is not None:
        kinds = [kinds] if isinstance(kinds, str) else kinds
        ok = {
            "object": isinstance(value, dict),
            "array": isinstance(value, list),
            "string": isinstance(value, str),
            "integer": isinstance(value, int) and not isinstance(value, bool),
            "number": isinstance(value, (int, float)) and not isinstance(value, bool),
            "boolean": isinstance(value, bool),
            "null": value is None,
        }
        if not any(ok[k] for k in kinds):
            return [f"{where}: {value!r} is not {kinds}"]
    if "enum" in sch and value not in sch["enum"]:
        return [f"{where}: {value!r} not in {sch['enum']}"]
    if "const" in sch and value != sch["const"]:
        return [f"{where}: {value!r} != {sch['const']!r}"]
    out = []
    if isinstance(value, dict):
        props = sch.get("properties", {})
        out += [f"{where}: missing {k}" for k in sch.get("required", []) if k not in value]
        if sch.get("additionalProperties") is False:
            out += [f"{where}: unexpected {k}" for k in value if k not in props]
        for k, v in value.items():
            if k in props:
                out += _conforms(v, props[k], f"{where}.{k}")
    if isinstance(value, list) and "items" in sch:
        for n, v in enumerate(value):
            out += _conforms(v, sch["items"], f"{where}[{n}]")
    return out


def _floats(value):
    if isinstance(value, float):
        yield value
    elif isinstance(value, dict):
        for v in value.values():
            yield from _floats(v)
    elif isinstance(value, list):
        for v in value:
            yield from _floats(v)


@pytest.fixture(scope="module")
def drawn():
    """1TIM and 1EMA drawn as `plot` draws them (the figure the data must agree with), and their data."""
    out = {}
    for name in ("1TIM", "1EMA"):
        d = draw_spec(figure_spec(DATA / f"{name}.cif"))
        out[name] = (d, structure_data(d.layout, d.sses, d.bb, d.model, name=str(DATA / f"{name}.cif")))
    return out


def test_schema_keys_and_version(drawn):
    sch = schema()
    assert SCHEMA_VERSION == 1 and ALL == SECTIONS
    assert sch["properties"]["foldmap_data"]["const"] == 1
    assert set(sch["properties"]) == {"foldmap_data", "foldmap", *SECTIONS}
    assert set(SECTIONS) == {"structure", "elements", "sheets", "disulfides", "symmetry", "membrane"}
    for name, (_, data) in drawn.items():
        assert data["foldmap_data"] == 1 and set(data) == {"foldmap_data", "foldmap", *SECTIONS}
        assert _conforms(data, sch) == [], name
        assert data["structure"]["name"] == f"{name}.cif"
        assert data["membrane"] is None  # none asked for
    json.dumps(sch)  # the schema itself is JSON


def test_numbers_are_rounded_and_never_negative_zero(drawn):
    for _, data in drawn.values():
        for x in _floats(data):
            assert x == round(x, 3) and not (x == 0 and np.signbit(x)), x
        for e in data["elements"]:
            if e["page"]:
                assert -180 < e["page"]["angle"] <= 180


@pytest.mark.parametrize("name", ["1TIM", "1EMA"])
def test_labels_and_sheet_orders_are_the_figures(drawn, name):
    """Each element's label is the text the figure writes on it (gid label:<id>); each sheet's strands stand left to
    right as the figure draws them (the strand arrows, gid strand:<id>, by their centre's page x)."""
    d, data = drawn[name]
    texts = {t.get_gid()[len("label:") :]: t.get_text() for t in d.figure.findobj() if _gid(t).startswith("label:")}
    elements = {e["id"]: e for e in data["elements"]}
    assert texts and set(texts) <= set(elements)
    for i, text in texts.items():
        assert elements[i]["label"] == text, i
    centre = {
        _gid(p)[len("strand:") :]: float(np.mean(p.get_xy()[:, 0]))
        for p in d.figure.findobj()
        if _gid(p).startswith("strand:")
    }
    assert data["sheets"]
    for sheet in data["sheets"]:
        assert sheet["as_drawn"]
        assert sheet["strands"] == sorted(sheet["strands"], key=lambda i: centre[i]), sheet["id"]
        for i in sheet["strands"]:
            assert elements[i]["sheet"] == sheet["id"] and elements[i]["kind"] == "strand"
        assert sheet["directions"][0] == 1 and len(sheet["directions"]) == len(sheet["strands"])
        n = len(sheet["strands"])
        assert len(sheet["pairs"]) == (n if sheet["closed"] else n - 1)
        for k, pair in enumerate(sheet["pairs"]):
            assert (pair["a"], pair["b"]) == (sheet["strands"][k], sheet["strands"][(k + 1) % n])
            assert pair["bridges"] == len(pair["register"]) > 0
            same = sheet["directions"][k] == sheet["directions"][(k + 1) % n]
            assert pair["kind"] == ("P" if same else "A")
            for a, b in pair["register"]:  # a's residue on strand a, b's on strand b
                for res, sid in ((a, pair["a"]), (b, pair["b"])):
                    e = elements[sid]
                    assert res["chain"] == e["chain"] and e["first"]["seq"] <= res["seq"] <= e["last"]["seq"]


def _gid(artist) -> str:
    return artist.get_gid() or ""


def test_tim_barrel_n_and_shear(drawn):
    """1TIM: each subunit's (βα)8 barrel, n = 8 strands, shear S = 8."""
    _, data = drawn["1TIM"]
    barrels = [s["barrel"] for s in data["sheets"] if s["closed"]]
    assert barrels == [{"n": 8, "shear": 8}, {"n": 8, "shear": 8}]
    assert all(set(s["directions"]) == {1} for s in data["sheets"])  # all parallel


def test_gfp_barrel(drawn):
    """1EMA: GFP's 11-stranded barrel, one parallel pairing (the strands A:12-22 and A:118-128)."""
    _, data = drawn["1EMA"]
    (sheet,) = data["sheets"]
    assert sheet["closed"] and sheet["barrel"]["n"] == 11 and sheet["barrel"]["shear"] is not None
    assert [p["kind"] for p in sheet["pairs"]].count("P") == 1
    labels = {e["id"]: e["label"] for e in data["elements"]}
    assert [labels[i] for i in sheet["strands"]][:3] == ["H", "I", "D"]


def test_symmetry_and_stoichiometry(drawn):
    _, data = drawn["1TIM"]
    sym = data["symmetry"]
    assert sym["found"]["label"] == "C2" and sym["found"]["protomers"] == [["A"], ["B"]]
    assert sym["stoichiometry"] == {"formula": "A2", "entities": [{"letter": "A", "chains": ["A", "B"], "copies": 2}]}
    _, gfp = drawn["1EMA"]
    assert gfp["symmetry"] == {
        "found": None,
        "drawn": None,
        "stoichiometry": {"formula": "A", "entities": [{"letter": "A", "chains": ["A"], "copies": 1}]},
    }


def test_disulfides_name_their_cysteines():
    path = DATA / "6ZS5.cif"
    d = draw_spec(figure_spec(path))
    data = structure_data(d.layout, d.sses, d.bb, d.model, ("disulfides", "elements"))
    assert set(data) == {"foldmap_data", "foldmap", "disulfides", "elements"}
    names = {(l.chain, l.seq, l.icode): l.name for l in d.bb.labels}
    assert len(data["disulfides"]) == len(d.layout.links.disulfides) > 0
    ids = {e["id"] for e in data["elements"]}
    for ss in data["disulfides"]:
        for end in (ss["a"], ss["b"]):
            assert names[(end["chain"], end["seq"], end["icode"])] == "CYS"
        assert all(e is None or e in ids for e in ss["elements"])


def test_membrane_segments_and_sides(capsys):
    """1C3W (bacteriorhodopsin) with --membrane auto: the helices whose CA run reaches both faces of the slab are its
    transmembrane segments, as A18 counts them (figcheck.membrane_view), each entering from one side and leaving to
    the other; every other element is out, in or in the slab by its depth (helices C and E stop short of a face of
    the estimated 30 Å slab: in it)."""
    from foldmap.cli import make_layout
    from foldmap.figcheck import membrane_view

    path = DATA / "1C3W.cif"
    assert main(["data", str(path), "--membrane", "auto", "--assembly", "asu"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert _conforms(data, schema()) == []
    mem = data["membrane"]
    assert mem is not None and mem["page"] is not None and mem["termini"] == [{"chain": "A", "n": "out", "c": "in"}]
    tm = [e for e in mem["elements"] if e["side"] == "tm"]
    lay, sses, bb = make_layout(path, membrane="auto", assembly="asu")
    assert [e["element"] for e in tm] == mem["transmembrane"]
    assert len(tm) == membrane_view(lay, sses, bb, lay.bilayer)["tm_helices"] >= 5
    assert [e["from"] for e in tm] == ["out", "in", "in", "in", "out"]
    assert all(e["from"] != e["to"] for e in tm)
    half = mem["thickness"] / 2
    for e in mem["elements"]:
        lo, hi = e["depth"]
        assert lo <= hi
        if e["side"] == "tm":
            assert lo <= -half + 2 and hi >= half - 2  # membrane.FACE_REACH of both faces
        elif e["side"] == "out":
            assert hi > half
        elif e["side"] == "in":
            assert lo < -half


def test_cli_writes_json_and_prints_schema(tmp_path, capsys):
    out = tmp_path / "d" / "ubq.json"
    assert main(["data", str(DATA / "1UBQ.cif"), "-o", str(out), "--sections", "elements,sheets"]) == 0
    data = json.loads(out.read_text())
    assert set(data) == {"foldmap_data", "foldmap", "elements", "sheets"}
    assert out.read_text() == dumps(data)
    assert main(["data", "--schema"]) == 0
    assert json.loads(capsys.readouterr().out) == schema()
    assert main(["data", str(DATA / "1UBQ.cif"), "--sections", "elements,colours"]) == 2
    assert "no data section colours" in capsys.readouterr().err
    assert main(["data"]) == 2


def test_names_and_renames_are_the_figures(capsys):
    """--rename reaches the data's labels, as it reaches the figure's."""
    assert main(["data", str(DATA / "1UBQ.cif"), "--rename", "res:A:23-34=Hx", "--sections", "elements"]) == 0
    labels = {e["id"]: e["label"] for e in json.loads(capsys.readouterr().out)["elements"]}
    assert labels["A:23-34"] == "Hx"


_RUN = """
import hashlib, sys
from foldmap import cli
assert cli.main(["data", sys.argv[1], "-o", sys.argv[2]]) == 0
print(hashlib.sha256(open(sys.argv[2], "rb").read()).hexdigest())
"""


def test_same_bytes_for_any_hash_seed(tmp_path):
    """A12: `foldmap data` writes the same bytes under PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _RUN, str(DATA / "1EMA.cif"), str(tmp_path / f"{seed}.json")],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout.split()[-1])
    assert len(runs[0]) == 64 and runs[0] == runs[1]
    assert hashlib.sha256((tmp_path / "0.json").read_bytes()).hexdigest() == runs[0]
