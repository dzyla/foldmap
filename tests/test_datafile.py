"""foldmap data (plan 13d-1, 13d-2, 13d-3): the versioned dict of a laid-out structure. Its schema, its labels and
sheet orders against the figure drawn from the same spec, a TIM barrel's n and S, ligand contacts (GFP's chromophore,
hemoglobin's hemes) and interfaces (hemoglobin's chain pairs, a repressor's protein-DNA pairs) against the
coordinates, the CSV tables, the nucleic acid as dot-bracket, BPSEQ and CT (read back, and by RNApolis where it is
installed), and the same bytes for any hash seed (A12)."""

import csv
import hashlib
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from foldmap.cli import main
from foldmap.datafile import (
    ALL,
    BRACKETS,
    DEFAULT,
    INTERACTION_KINDS,
    SCHEMA_VERSION,
    SECTIONS,
    TABLES,
    bpseq_text,
    csv_text,
    ct_text,
    dot_bracket_text,
    dumps,
    nucleic_acid,
    pair_levels,
    schema,
    structure_data,
    tables,
    zip_bytes,
)
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


@pytest.fixture(scope="module")
def hemoglobin():
    """2HHB (deoxy hemoglobin, α2β2) as `plot` draws it, and its data with every section, interactions included."""
    d = draw_spec(figure_spec(DATA / "2HHB.cif"))
    return d, structure_data(d.layout, d.sses, d.bb, d.model, SECTIONS, name=str(DATA / "2HHB.cif"))


def test_schema_keys_and_version(drawn, hemoglobin):
    from foldmap.interactions import KINDS

    sch = schema()
    assert SCHEMA_VERSION == 1 and ALL == SECTIONS
    assert sch["properties"]["foldmap_data"]["const"] == 1
    assert set(sch["properties"]) == {"foldmap_data", "foldmap", *SECTIONS}
    assert set(SECTIONS) == {
        "structure",
        "elements",
        "sheets",
        "disulfides",
        "symmetry",
        "membrane",
        "ligands",
        "interfaces",
        "nucleic_acid",
        "interactions",
    }
    assert DEFAULT == tuple(s for s in SECTIONS if s != "interactions")  # the finders run only when asked for
    assert INTERACTION_KINDS == KINDS
    _, hb = hemoglobin
    assert set(hb) == {"foldmap_data", "foldmap", *SECTIONS} and _conforms(hb, sch) == []
    for name, (_, data) in drawn.items():
        assert data["foldmap_data"] == 1 and set(data) == {"foldmap_data", "foldmap", *DEFAULT}
        assert _conforms(data, sch) == [], name
        assert data["structure"]["name"] == f"{name}.cif"
        assert data["membrane"] is None  # none asked for
    json.dumps(sch)  # the schema itself is JSON


def test_numbers_are_rounded_and_never_negative_zero(drawn, hemoglobin):
    for _, data in [*drawn.values(), hemoglobin]:
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


@pytest.fixture(scope="module")
def glycoprotein():
    """6ZS5, with disulfides and glycans."""
    d = draw_spec(figure_spec(DATA / "6ZS5.cif"))
    return d, structure_data(d.layout, d.sses, d.bb, d.model, ("disulfides", "elements", "ligands"))


def test_disulfides_name_their_cysteines(glycoprotein):
    d, data = glycoprotein
    assert set(data) == {"foldmap_data", "foldmap", "disulfides", "elements", "ligands"}
    names = {(l.chain, l.seq, l.icode): l.name for l in d.bb.labels}
    assert len(data["disulfides"]) == len(d.layout.links.disulfides) > 0
    ids = {e["id"] for e in data["elements"]}
    for ss in data["disulfides"]:
        for end in (ss["a"], ss["b"]):
            assert names[(end["chain"], end["seq"], end["icode"])] == "CYS"
        assert all(e is None or e in ids for e in ss["elements"])


def test_glycans_hang_on_their_residue(glycoprotein):
    """6ZS5: one glycan entry per glycan the figure draws (layout.links.glycans), in that order: its sugars as
    codes, sitting at its root sugar, covalently attached to the one Asn, Ser or Thr it hangs on (also a contact),
    its contacts those a brute-force reading of the coordinates gives (4.0 Å, heavy atoms, over all its sugars)."""
    d, data = glycoprotein
    glycans = [g for g in data["ligands"] if g["kind"] == "glycan"]
    links = d.layout.links
    assert len(glycans) == len(links.glycans) > 0
    protein = {(l.chain, l.seq) for l in d.bb.labels}
    for g, (site, sugars), members in zip(glycans, links.glycans, links.glycan_members):
        anchor = d.bb.labels[site]
        assert g["codes"] == sugars and [(r["chain"], r["seq"]) for r in g["residues"]] == members
        assert (g["at"]["chain"], g["at"]["seq"]) == members[0] and g["id"] == f"glycan:{_ref(g['at'])}"
        assert [_ref(r) for r in g["attached_to"]] == [f"{anchor.chain}:{anchor.tag}"]
        assert anchor.name in ("ASN", "SER", "THR")
        seen = {(c["residue"]["chain"], c["residue"]["seq"]) for c in g["contacts"]}
        assert (anchor.chain, anchor.seq) in seen
        assert seen == _touching(d.model, set(members), protein, 4.0)
        for c in g["contacts"]:
            assert c["element"] == _element_of(data, c["residue"]["chain"], c["residue"]["seq"])


def test_membrane_segments_and_sides(capsys):
    """1C3W (bacteriorhodopsin) with --membrane auto. spans_slab (CA run within 2 Å of both faces of the drawn slab,
    A18's rule before plan 26): 5 helices, C and E stopping short of a face of the estimated 30 Å slab.
    crosses_midplane (the CA run crosses the centre plane and spans half the thickness: membrane.crosses, the rule
    figcheck.membrane_view counts by) reads the usual 7 helices, alternating in and out; every other element is out,
    in or in the slab by its depth."""
    from foldmap.cli import make_layout
    from foldmap.figcheck import membrane_view

    path = DATA / "1C3W.cif"
    assert main(["data", str(path), "--membrane", "auto", "--assembly", "asu"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert _conforms(data, schema()) == []
    mem = data["membrane"]
    assert mem is not None and mem["page"] is not None and mem["termini"] == [{"chain": "A", "n": "out", "c": "in"}]
    lay, sses, bb = make_layout(path, membrane="auto", assembly="asu")
    spans = [e["element"] for e in mem["elements"] if e["spans_slab"]]
    assert spans == mem["spans_slab"] and len(spans) == 5
    across = [e for e in mem["elements"] if e["crosses_midplane"]]
    assert [e["element"] for e in across] == mem["crosses_midplane"] and len(across) == 7
    assert membrane_view(lay, sses, bb, lay.bilayer)["tm_helices"] == 7
    assert set(spans) <= set(mem["crosses_midplane"])
    assert all(e["side"] == "across" for e in across)
    for a, b in zip(across, across[1:]):  # each helix leaves on the side the next one enters from
        assert a["from"] != a["to"] and a["to"] == b["from"]
    assert across[0]["from"] == "out" and across[-1]["to"] == "in"
    half = mem["thickness"] / 2
    for e in mem["elements"]:
        lo, hi = e["depth"]
        assert lo <= hi
        if e["spans_slab"]:
            assert lo <= -half + 2 and hi >= half - 2  # membrane.FACE_REACH of both faces
        if e["crosses_midplane"]:
            assert lo < 0 < hi and hi - lo >= half
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


@pytest.mark.parametrize("name", ["1EMA", "1LMB"])
def test_same_bytes_for_any_hash_seed(tmp_path, name):
    """A12: `foldmap data` writes the same bytes under PYTHONHASHSEED 0 and 1 (1LMB: a protein-DNA complex, its
    interfaces and nucleic acid)."""
    runs = []
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _RUN, str(DATA / f"{name}.cif"), str(tmp_path / f"{seed}.json")],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout.split()[-1])
    assert len(runs[0]) == 64 and runs[0] == runs[1]
    assert hashlib.sha256((tmp_path / "0.json").read_bytes()).hexdigest() == runs[0]


# ---------------------------------------------------------------------------------------------- 13d-2


def _ref(r: dict) -> str:
    return f"{r['chain']}:{r['seq']}{r['icode']}"


def _element_of(data: dict, chain: str, seq: int) -> str | None:
    """The element a residue is in, read from the data's own elements (first..last), None in a loop."""
    for e in data["elements"]:
        if e["chain"] == chain and e["first"]["seq"] <= seq <= e["last"]["seq"]:
            return e["id"]
    return None


def _touching(model, own: set, wanted: set, reach: float) -> set:
    """(chain, number) of the residues in `wanted` with a heavy atom within `reach` of a heavy atom of the residues
    `own` ((chain, number)): read here from the coordinates, by brute force."""
    if isinstance(model, tuple):
        model = model[0]
    atoms = {}
    for chain in model:
        for res in chain:
            key = (chain.name, res.seqid.num)
            if key in own or key in wanted:
                atoms.setdefault(key, []).extend(a.pos.tolist() for a in res if a.element.name not in ("H", "D"))
    mine = np.array([p for k in own for p in atoms[k]], float)
    out = set()
    for key in wanted - own:
        theirs = np.array(atoms.get(key, []), float).reshape(-1, 3)
        if len(theirs) and np.min(np.linalg.norm(theirs[:, None] - mine[None], axis=2)) <= reach:
            out.add(key)
    return out


def test_chromophore_contacts_on_its_elements(drawn):
    """1EMA: GFP's chromophore CRO66, fused into the chain, is listed as a linked residue. Its contacts are the
    residues a brute-force reading of the coordinates puts within 4.0 Å of it (heavy atoms), each on the element it
    is in, among them the pocket's Gln94, Arg96, His148, Thr203 and Glu222 on strands D, D, G, J and K; it is joined
    to Phe64 and Val68 by peptide bonds."""
    d, data = drawn["1EMA"]
    (cro,) = data["ligands"]
    assert cro["id"] == "linked:A:66" and cro["kind"] == "linked" and cro["name"] == "CRO"
    assert [_ref(r) for r in cro["attached_to"]] == ["A:64", "A:68"]
    protein = {(l.chain, l.seq) for l in d.bb.labels}
    seen = [(c["residue"]["chain"], c["residue"]["seq"]) for c in cro["contacts"]]
    assert set(seen) == _touching(d.model, {("A", 66)}, protein, 4.0) and len(seen) == len(set(seen))
    assert seen == sorted(seen, key=lambda r: r[1])  # backbone order
    labels = {e["id"]: e["label"] for e in data["elements"]}
    kinds = {e["id"]: e["kind"] for e in data["elements"]}
    for c in cro["contacts"]:
        assert c["element"] == _element_of(data, c["residue"]["chain"], c["residue"]["seq"])
    pocket = {c["residue"]["seq"]: labels.get(c["element"]) for c in cro["contacts"]}
    assert {k: pocket.get(k) for k in (94, 96, 148, 203, 222)} == {94: "D", 96: "D", 148: "G", 203: "J", 222: "K"}
    order = [e["id"] for e in data["elements"]]
    assert cro["elements"] == sorted({c["element"] for c in cro["contacts"] if c["element"]}, key=order.index)
    assert sum(kinds[e] == "strand" for e in cro["elements"]) >= 6  # it sits inside the barrel


def test_hemes_touch_their_helices(hemoglobin):
    """2HHB: each heme touches its own chain, the proximal His (α His87, β His92) on the F helix (the figure's α6);
    the phosphate is marked an additive. The contacts are those of io's ligand rule (4.0 Å, heavy atoms)."""
    d, data = hemoglobin
    labels = {e["id"]: e["label"] for e in data["elements"]}
    protein = {(l.chain, l.seq) for l in d.bb.labels}
    hemes = [g for g in data["ligands"] if g["name"] == "HEM"]
    assert [g["id"] for g in hemes] == ["ligand:A:142", "ligand:B:148", "ligand:C:142", "ligand:D:148"]
    for g in hemes:
        chain = g["at"]["chain"]
        seen = {(c["residue"]["chain"], c["residue"]["seq"]) for c in g["contacts"]}
        assert seen == _touching(d.model, {(chain, g["at"]["seq"])}, protein, 4.0)
        assert {c for c, _ in seen} == {chain} and g["attached_to"] == []
        his = 87 if chain in "AC" else 92
        (proximal,) = [c for c in g["contacts"] if c["residue"]["seq"] == his]
        assert labels[proximal["element"]] == "α6"
    (po4,) = [g for g in data["ligands"] if g["name"] == "PO4"]
    assert po4["additive"] and po4["kind"] == "ligand"


def test_hemoglobin_interfaces_by_chain_pair(hemoglobin):
    """2HHB: the chain pairs in contact are α1β1 and α2β2 (A-B, C-D), α1β2 and α2β1 (A-D, B-C), and α1α2 (A-C); the
    β chains do not touch in deoxy hemoglobin (2,3-BPG's cavity). CA pairs are counted within 8 Å; each
    element pair is features.sse_contacts' and the like interfaces pair the same helices. The interactions (13c's
    finders, between chains) add up by element pair, and include the T-state salt bridges: Lys40α to the C-terminal
    carboxylate of His146β across α1β2, and Arg141α to Asp126α across α1α2."""
    from foldmap.features import sse_contacts

    d, data = hemoglobin
    faces = {tuple(f["chains"]): f for f in data["interfaces"]}
    assert list(faces) == [("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("C", "D")]
    chains = np.array([l.chain for l in d.bb.labels])
    dist = np.linalg.norm(d.bb.ca[:, None] - d.bb.ca[None], axis=2) <= 8.0
    for (a, b), f in faces.items():
        near = dist[np.ix_(chains == a, chains == b)]
        assert f["ca_pairs"] == int(near.sum()) and f["residues"] == [int(near.any(1).sum()), int(near.any(0).sum())]
    chain_of = {s.id: s.chain for s in d.sses}
    want = {k: v for k, v in sse_contacts(d.bb, d.sses).items() if len({chain_of[i] for i in k}) == 2}
    got = {frozenset((e["a"], e["b"])): e["ca_pairs"] for f in faces.values() for e in f["elements"] if e["ca_pairs"]}
    assert got == want
    labels = {e["id"]: e["label"] for e in data["elements"]}

    def helices(pair):
        return sorted((labels[e["a"]], labels[e["b"]]) for e in faces[pair]["elements"] if e["ca_pairs"])

    assert helices(("A", "B")) == helices(("C", "D")) == [("α2", "α8"), ("α7", "α7"), ("α8", "α2")]
    assert helices(("A", "D")) == [("α6", "α3")] and helices(("B", "C")) == [("α3", "α6")]
    for (a, b), f in faces.items():
        for e in f["elements"]:
            assert e["a"] is None or chain_of[e["a"]] == a
            assert e["b"] is None or chain_of[e["b"]] == b
    found = data["interactions"]["found"]
    assert data["interactions"]["between"] == "interchain" and data["interactions"]["failed_recheck"] == []
    for (a, b), f in faces.items():
        mine = [x for x in found if {x["a"]["chain"], x["b"]["chain"]} == {a, b}]
        assert f["interactions"] == {k: sum(x["kind"] == k for x in mine) for k in INTERACTION_KINDS}
        for k in INTERACTION_KINDS:
            assert sum(e["interactions"][k] for e in f["elements"]) == f["interactions"][k]
    assert not [x for x in found if {x["a"]["chain"], x["b"]["chain"]} == {"B", "D"}]
    for x in found:
        assert x["a"]["chain"] != x["b"]["chain"]
        for end, element in zip((x["a"], x["b"]), x["elements"]):
            assert element == _element_of(data, end["chain"], end["seq"])
    salt = {(_ref(x["a"]), _ref(x["b"])) for x in found if x["kind"] == "salt"}
    assert {("A:40", "D:146"), ("C:40", "B:146"), ("A:141", "C:126"), ("C:141", "A:126")} <= salt


def _read_csv(text: str) -> list[list[str]]:
    return list(csv.reader(io.StringIO(text)))


def test_csv_tables_headers_and_row_order(tmp_path, capsys):
    """`-o out.csv --table NAME` writes one table: its header is TABLES' columns and its rows are the JSON's list in
    the JSON's order; with no -o the CSV goes to standard output. A .csv without --table, an unknown table, and a
    table whose section --sections leaves out are refused (exit 2)."""
    path = str(DATA / "1TIM.cif")
    assert main(["data", path, "--sections", "elements,sheets"]) == 0
    data = json.loads(capsys.readouterr().out)
    out = tmp_path / "e.csv"
    assert main(["data", path, "-o", str(out), "--table", "elements"]) == 0
    text = out.read_text(encoding="utf-8")
    assert "\r" not in text and text.endswith("\n")
    rows = _read_csv(text)
    assert tuple(rows[0]) == TABLES["elements"][1]
    assert [r[0] for r in rows[1:]] == [e["id"] for e in data["elements"]]
    col = {c: k for k, c in enumerate(rows[0])}
    for row, e in zip(rows[1:], data["elements"]):
        assert row[col["label"]] == e["label"] and float(row[col["x"]]) == e["page"]["x"]
        assert row[col["first_seq"]] == str(e["first"]["seq"]) and row[col["lent"]] == "false"
        assert row[col["domain"]] == ""  # null is empty
    capsys.readouterr()
    assert main(["data", path, "--table", "sheet_pairs"]) == 0
    rows = _read_csv(capsys.readouterr().out)
    assert tuple(rows[0]) == ("sheet", "a", "b", "kind", "bridges")
    assert rows[1:] == [
        [s["id"], p["a"], p["b"], p["kind"], str(p["bridges"])] for s in data["sheets"] for p in s["pairs"]
    ]
    assert main(["data", path, "--table", "sheets"]) == 0
    rows = _read_csv(capsys.readouterr().out)
    run = [r for r in rows[1:] if r[3] == "run"]
    assert [(r[0], r[1], int(r[2]), int(r[4])) for r in run] == [
        (s["id"], i, n, sense)
        for s in data["sheets"]
        for n, (i, sense) in enumerate(zip(s["strands"], s["directions"]), 1)
    ]
    assert {(r[8], r[9]) for r in run} == {("8", "8")}  # barrel n and shear on every strand of a TIM barrel
    assert main(["data", path, "-o", str(tmp_path / "x.csv")]) == 2
    assert "name it with --table" in capsys.readouterr().err
    assert main(["data", path, "--table", "colours"]) == 2
    assert "no data table colours" in capsys.readouterr().err
    assert main(["data", path, "--table", "interactions", "--sections", "elements"]) == 2
    assert "leaves out" in capsys.readouterr().err
    assert csv_text(("a", "b", "c"), [[None, True, "x,y"]]) == 'a,b,c\n,true,"x,y"\n'


def test_dir_writes_every_table(tmp_path, hemoglobin):
    """`-o DIR/` writes one CSV per table of the sections asked for, each the table of the data structure_data gives
    for the same structure (every section asked for: every table)."""
    _, data = hemoglobin
    out = tmp_path / "tables"
    assert main(["data", str(DATA / "2HHB.cif"), "-o", f"{out}/", "--sections", ",".join(SECTIONS)]) == 0
    assert sorted(p.name for p in out.iterdir()) == sorted(f"{name}.csv" for name in TABLES)
    want = tables(data)
    assert list(want) == list(TABLES)
    for name, (columns, rows) in want.items():
        text = (out / f"{name}.csv").read_text(encoding="utf-8")
        assert text == csv_text(columns, rows), name
        assert tuple(_read_csv(text)[0]) == TABLES[name][1]
    assert len(want["interactions"][1]) == len(data["interactions"]["found"]) > 0
    assert len(want["ligand_contacts"][1]) == sum(len(g["contacts"]) for g in data["ligands"])
    assert len(want["interface_elements"][1]) == sum(len(f["elements"]) for f in data["interfaces"])
    default = tmp_path / "default"
    default.mkdir()
    assert main(["data", str(DATA / "1UBQ.cif"), "-o", str(default)]) == 0  # an existing directory, no slash
    assert sorted(p.stem for p in default.iterdir()) == sorted(n for n, t in TABLES.items() if t[0] in DEFAULT)


_RUN_TABLES = """
import hashlib, pathlib, sys
from foldmap import cli
from foldmap.datafile import SECTIONS
assert cli.main(["data", sys.argv[1], "-o", sys.argv[2] + "/", "--sections", ",".join(SECTIONS)]) == 0
h = hashlib.sha256()
for p in sorted(pathlib.Path(sys.argv[2]).iterdir()):
    h.update(p.name.encode() + b"\\0" + p.read_bytes())
print(h.hexdigest())
"""


def test_tables_same_bytes_for_any_hash_seed(tmp_path):
    """A12: every CSV table of 2HHB (ligands, interfaces and interactions included) has the same bytes under
    PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _RUN_TABLES, str(DATA / "2HHB.cif"), str(tmp_path / seed)],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout.split()[-1])
    assert len(runs[0]) == 64 and runs[0] == runs[1]


# ---------------------------------------------------------------------------------------------- nucleic acid (13d-3)

NA = DATA / "na"
NA_FIXTURES = ["30QM", "1EHZ", "4OO8", "1KX5", "1IHF", "3CRX", "4YBB", "1RVA", "4LD0", "1HYS", "1AOI", "2GIS"]


def _nucleic(pdb):
    from foldmap.io import load_backbone, load_nucleic, read_model

    path = NA / f"{pdb}.cif.gz"
    model = read_model(path, "asu")
    return load_nucleic(path, load_backbone(path, "asu", model=model), "asu", model=model)


def _read_dot_bracket(text: str) -> tuple[list[str], str, set[tuple[int, int]]]:
    """(strand headers, sequence, pairs as 1-based (i, j)) of a multi-strand dot-bracket: '>' lines, then a sequence
    and a structure line each; brackets matched across the strands read one after another."""
    lines = text.splitlines()
    assert len(lines) % 3 == 0
    heads, seq, marks = lines[0::3], "".join(lines[1::3]), "".join(lines[2::3])
    assert all(h.startswith(">") for h in heads) and len(seq) == len(marks)
    for s_, m_ in zip(lines[1::3], lines[2::3]):
        assert len(s_) == len(m_)
    opening = {o: c for o, c in BRACKETS}
    closing = {c: o for o, c in BRACKETS}
    open_at: dict[str, list[int]] = {o: [] for o in opening}
    pairs = set()
    for n, c in enumerate(marks, 1):
        if c in opening:
            open_at[c].append(n)
        elif c in closing:
            pairs.add((open_at[closing[c]].pop(), n))
        else:
            assert c == "."
    assert not any(open_at.values())
    return heads, seq, pairs


def _read_bpseq(text: str) -> tuple[str, set[tuple[int, int]]]:
    rows = [line.split() for line in text.splitlines()]
    assert all(len(r) == 3 for r in rows)
    assert [int(r[0]) for r in rows] == list(range(1, len(rows) + 1))
    partner = {int(r[0]): int(r[2]) for r in rows}
    assert all(partner[j] == i for i, j in partner.items() if j)  # each pair written from both ends
    return "".join(r[1] for r in rows), {(i, j) for i, j in partner.items() if 0 < i < j}


def _read_ct(text: str) -> tuple[str, set[tuple[int, int]], list[list[str]]]:
    lines = text.splitlines()
    n = int(lines[0].split()[0])
    rows = [line.split() for line in lines[1:]]
    assert len(rows) == n and all(len(r) == 6 for r in rows)
    assert [int(r[0]) for r in rows] == list(range(1, n + 1))
    partner = {int(r[0]): int(r[4]) for r in rows}
    assert all(partner[j] == i for i, j in partner.items() if j)
    return "".join(r[1] for r in rows), {(i, j) for i, j in partner.items() if 0 < i < j}, rows


def _positions(na) -> dict[tuple[int, int], int]:
    out, n = {}, 0
    for s, st in enumerate(na.strands):
        for i in range(len(st)):
            n += 1
            out[(s, i)] = n
    return out


@pytest.mark.parametrize("pdb", NA_FIXTURES)
def test_dot_bracket_bpseq_and_ct_give_back_the_pairs_they_were_written_from(pdb):
    """Every NA fixture: the base pairs read back from its dot-bracket, BPSEQ and CT are the pairs written
    (nucleic.py's, the figure's), on the same sequence; a pair between strands is bracketed across them."""
    na = _nucleic(pdb)
    data = {"nucleic_acid": nucleic_acid(na)}
    where = _positions(na)
    want = {(where[(p.a, p.i)], where[(p.b, p.j)]) for p in na.base_pairs}
    sequence = "".join(st.bases for st in na.strands)
    heads, seq, got = _read_dot_bracket(dot_bracket_text(data))
    assert heads == [f">{st.id}" for st in na.strands] and seq == sequence and got == want
    assert _read_bpseq(bpseq_text(data)) == (sequence, want)
    seq, got, _ = _read_ct(ct_text(data))
    assert (seq, got) == (sequence, want)
    assert {(p["i"], p["j"]) for p in data["nucleic_acid"]["pairs"]} == want
    assert data["nucleic_acid"]["nucleotides"] == len(sequence) == na.size


@pytest.mark.parametrize("pdb", ["1EHZ", "4LD0", "1KX5", "1AOI"])
def test_dot_bracket_pairs_agree_with_rnapolis(pdb):
    """The pairs read back from the dot-bracket, named by residue, against RNApolis' canonical and wobble pairs
    (tests/data/na/reference.json, the NA campaign's): the same agreement test_na_pairs asks of the pairing."""
    reference = json.loads((NA / "reference.json").read_text(encoding="utf-8"))
    na = _nucleic(pdb)
    data = {"nucleic_acid": nucleic_acid(na)}
    residues = [tuple(r.values()) for st in data["nucleic_acid"]["strands"] for r in st["residues"]]
    _, _, pairs = _read_dot_bracket(dot_bracket_text(data))
    ours = {frozenset((residues[i - 1], residues[j - 1])) for i, j in pairs}
    theirs = {frozenset((tuple(a), tuple(b))) for a, b, _ in reference[pdb]}
    both = ours & theirs
    assert len(both) / len(ours) >= 0.95 and len(both) / len(theirs) >= 0.90, pdb


_RNAPOLIS = """
import json, sys
from pathlib import Path
from rnapolis.common import BpSeq, MultiStrandDotBracket
db = MultiStrandDotBracket.from_multiline_string(Path(sys.argv[1]).read_text())
bp = BpSeq.from_string(Path(sys.argv[2]).read_text())
print(json.dumps({
    "dbn": sorted([i + 1, j + 1] for i, j in db.pairs),
    "dbn_sequence": db.sequence,
    "strands": len(db.strands),
    "bpseq": sorted([i, j] for i, _, j in bp.entries if 0 < i < j),
}))
"""


@pytest.mark.skipif(not os.environ.get("NA_RNAPOLIS_PYTHON"), reason="RNApolis' venv (NA_RNAPOLIS_PYTHON): zyla-office")
@pytest.mark.parametrize("pdb", ["1EHZ", "2GIS", "1KX5", "4YBB"])
def test_rnapolis_reads_our_dot_bracket_and_bpseq(pdb, tmp_path):
    """RNApolis' own readers (the NA campaign's venv: MultiStrandDotBracket, BpSeq) take our dot-bracket and BPSEQ
    and find the pairs, strands and sequence we wrote: a pseudoknot (2GIS), strands paired across (1KX5) and a
    ribosome's many strands and levels (4YBB)."""
    na = _nucleic(pdb)
    data = {"nucleic_acid": nucleic_acid(na)}
    (tmp_path / "x.dbn").write_text(dot_bracket_text(data))
    (tmp_path / "x.bpseq").write_text(bpseq_text(data))
    files = [str(tmp_path / f"x.{ext}") for ext in ("dbn", "bpseq")]
    done = subprocess.run(
        [os.environ["NA_RNAPOLIS_PYTHON"], "-I", "-c", _RNAPOLIS, *files], capture_output=True, text=True, timeout=600
    )
    assert done.returncode == 0, done.stderr[-2000:]
    got = json.loads(done.stdout.splitlines()[-1])
    pairs = sorted([p["i"], p["j"]] for p in data["nucleic_acid"]["pairs"])
    assert got["dbn"] == pairs and got["bpseq"] == pairs
    assert got["dbn_sequence"] == "".join(st.bases for st in na.strands)
    assert got["strands"] == len(na.strands)


def test_a_pseudoknot_is_bracketed_by_foldmaps_rule():
    """2GIS (SAM-I riboswitch): the stem foldmap's rule calls the pseudoknot (PK, A:25-28/A:65-68, the one the figure
    names) is written with [], every other stem with (); the levels are those of pair_levels."""
    na = _nucleic("2GIS")
    sec = nucleic_acid(na)
    (strand,) = sec["strands"]
    marks = dict(zip((r["seq"] for r in strand["residues"]), strand["dot_bracket"]))
    (pk,) = [st for st in na.structure.stems if st.name == "PK"]
    assert pk.id == "A:25-28/A:65-68"
    for bp in pk.pairs:
        assert marks[na.strands[bp.a].labels[bp.i].seq] == "[" and marks[na.strands[bp.b].labels[bp.j].seq] == "]"
    assert set(strand["dot_bracket"]) <= {".", "(", ")", "[", "]"}
    assert strand["dot_bracket"].count("[") == len(pk)
    for st in sec["stems"]:
        assert st["pseudoknot"] == (st["name"] == "PK") and st["level"] == (1 if st["name"] == "PK" else 0)
    assert sec["levels"] == 2
    levels = pair_levels(na)
    assert [p["level"] for p in sec["pairs"]] == [levels[tuple(p[:4])] for p in na.base_pairs]
    for p in sec["pairs"]:
        assert (p["stem"] == "PK") == (p["level"] == 1)


def test_crossing_pseudoknots_climb_the_bracket_levels():
    """Synthetic: a hairpin, a knot crossing it, a shorter knot crossing both, and an isolated pair crossing none
    of them: (), [], {} and () again, in the order foldmap's rule takes stems (longest first)."""
    from foldmap.model import NAStrand, ResLabel, nucleic_from_pairs
    from foldmap.nucleic import BasePair

    n = 60
    strand = NAStrand("A", [ResLabel("A", k + 1, "", "A") for k in range(n)], np.zeros((n, 3)), np.zeros((n, 3)), True)
    pairs = [(k, 20 - k) for k in range(4)]  # the hairpin 0-3/17-20
    pairs += [(10 + k, 30 - k) for k in range(3)]  # crosses the hairpin
    pairs += [(15 + k, 40 - k) for k in range(2)]  # crosses both
    pairs += [(50, 58)]  # an isolated pair inside nothing
    na = nucleic_from_pairs([strand], [BasePair(0, i, 0, j, "WC") for i, j in pairs], {}, {})
    sec = nucleic_acid(na)
    db = sec["strands"][0]["dot_bracket"]
    assert db[0:4] == "((((" and db[17:21] == "))))"
    assert db[10:13] == "[[[" and db[28:31] == "]]]"
    assert db[15:17] == "{{" and db[39:41] == "}}"
    assert db[50] == "(" and db[58] == ")"
    assert sec["levels"] == 3
    _, _, got = _read_dot_bracket(dot_bracket_text({"nucleic_acid": sec}))
    assert got == {(i + 1, j + 1) for i, j in pairs}


def test_cli_writes_dot_bracket_bpseq_and_ct(tmp_path, capsys):
    """`foldmap data X -o F.dbn | .bpseq | .ct` (or --format to standard output) writes the nucleic acid of the
    figure `plot` draws: the base pairs of its section, the columns each format has (BPSEQ: number, base, partner;
    CT: number, base, the one before, the one after, partner, residue number), the same for every route. A
    structure with no nucleic acid is an error; a format with --table is refused."""
    path = str(NA / "2GIS.cif.gz")
    assert main(["data", path, "--sections", "structure,nucleic_acid"]) == 0
    data = json.loads(capsys.readouterr().out)
    sec = data["nucleic_acid"]
    pairs = {(p["i"], p["j"]) for p in sec["pairs"]}
    sequence = "".join(st["sequence"] for st in sec["strands"])
    for ext in ("dbn", "bpseq", "ct"):
        assert main(["data", path, "-o", str(tmp_path / f"x.{ext}")]) == 0
        assert main(["data", path, "--format", ext]) == 0
        assert capsys.readouterr().out == (tmp_path / f"x.{ext}").read_text(encoding="utf-8")
    text = (tmp_path / "x.dbn").read_text(encoding="utf-8")
    assert text == dot_bracket_text(data)
    assert _read_dot_bracket(text)[1:] == (sequence, pairs)
    assert _read_bpseq((tmp_path / "x.bpseq").read_text(encoding="utf-8")) == (sequence, pairs)
    ct = (tmp_path / "x.ct").read_text(encoding="utf-8")
    assert ct.splitlines()[0] == f"{len(sequence)} 2GIS.cif.gz"
    seq, got, rows = _read_ct(ct)
    assert (seq, got) == (sequence, pairs)
    residues = [r for st in sec["strands"] for r in st["residues"]]
    assert [int(r[5]) for r in rows] == [r["seq"] for r in residues]
    for st in sec["strands"]:
        first, last = st["first"], st["first"] + len(st["residues"]) - 1
        for n in range(first, last + 1):
            assert int(rows[n - 1][2]) == (n - 1 if n > first else 0)
            assert int(rows[n - 1][3]) == (n + 1 if n < last else 0)
    assert main(["data", str(DATA / "1UBQ.cif"), "--format", "dbn"]) == 1
    assert "no nucleic acid" in capsys.readouterr().err
    assert main(["data", path, "--format", "ct", "--table", "elements"]) == 2
    assert main(["data", path, "-o", str(tmp_path / "y.ct"), "--sections", "elements"]) == 2
    assert "leaves out" in capsys.readouterr().err


def test_multistrand_dna_breaks_the_ct_chain_between_strands():
    """1KX5 (nucleosome): the two DNA strands are two dot-bracket blocks whose brackets pair across them; in CT the
    one-after column is 0 at the end of each strand and the one-before 0 at the start of the next."""
    na = _nucleic("1KX5")
    sec = nucleic_acid(na)
    assert len(sec["strands"]) >= 2 and {st["chain"] for st in sec["strands"]} == {"I", "J"}
    across = [p for p in sec["pairs"] if p["a"]["chain"] != p["b"]["chain"]]
    assert len(across) > 100 and all(p["level"] == 0 for p in across)
    _, _, rows = _read_ct(ct_text({"nucleic_acid": sec}))
    for st in sec["strands"][1:]:
        n = st["first"]
        assert rows[n - 2][3] == "0" and rows[n - 1][2] == "0"


def test_protein_nucleic_acid_interfaces_without_the_interactions_section():
    """1LMB (λ repressor on DNA): the interfaces list each protein chain with each DNA strand it touches, measured
    on the structure as residue-nucleotide pairs with heavy atoms within 4 Å (io.CONTACT), without the interactions
    section; the two DNA strands are one interface by the base pairs the figure draws. Counts against a brute-force
    reading of the coordinates."""
    import gemmi

    from foldmap.io import CONTACT

    d = draw_spec(figure_spec(DATA / "1LMB.cif"))
    data = structure_data(d.layout, d.sses, d.bb, d.model, DEFAULT)
    faces = {tuple(f["chains"]): f for f in data["interfaces"]}
    na = d.layout.nucleic
    dna = sorted({st.chain for st in na.strands})
    protein = sorted({l.chain for l in d.bb.labels})
    model = d.model[0] if isinstance(d.model, tuple) else d.model

    def heavy(res):
        return np.array([a.pos.tolist() for a in res if a.element.name not in ("H", "D")])

    residues = {c.name: [r for r in c] for c in model}
    nucleotides = {c: [r for r in residues[c] if r.find_atom("C1'", "*") is not None] for c in dna}
    where = {(l.chain, l.seq, l.icode) for l in d.bb.labels}
    amino = {c: [r for r in residues[c] if (c, r.seqid.num, r.seqid.icode.strip()) in where] for c in protein}
    seen = 0
    for p in protein:
        for q in dna:
            n, mine, theirs = 0, set(), set()
            for r in amino[p]:
                x = heavy(r)
                for t in nucleotides[q]:
                    if (np.linalg.norm(x[:, None] - heavy(t)[None], axis=2) <= CONTACT).any():
                        n += 1
                        mine.add(r.seqid.num)
                        theirs.add(t.seqid.num)
            if n:
                seen += 1
                f = faces[(p, q)]
                assert f["measure"] == "atoms" and f["ca_pairs"] is None and f["interactions"] is None
                assert f["contacts"] == n and f["residues"] == [len(mine), len(theirs)], (p, q)
            else:
                assert (p, q) not in faces
    assert seen >= 2
    pair = tuple(dna)
    across = sum(1 for bp in na.base_pairs if na.strands[bp.a].chain != na.strands[bp.b].chain)
    assert faces[pair]["measure"] == "base_pairs" and faces[pair]["contacts"] == across > 10
    for f in data["interfaces"]:
        if f["measure"] == "ca":
            assert f["contacts"] == f["ca_pairs"]
    sch = schema()
    assert _conforms(data, sch) == []
    sec = data["nucleic_acid"]
    assert {(p["i"], p["j"]) for p in sec["pairs"]} == {
        (_positions(na)[(p.a, p.i)], _positions(na)[(p.b, p.j)]) for p in na.base_pairs
    }
    assert gemmi  # the model is gemmi's


def test_nucleic_tables_and_zip(tmp_path):
    """The nucleotides and base_pairs tables follow the section (one row per nucleotide, per pair); the zip holds
    every table and the three nucleic-acid files, the same bytes each time."""
    import zipfile

    na = _nucleic("1EHZ")
    data = {"structure": {"name": "1EHZ.cif.gz"}, "nucleic_acid": nucleic_acid(na)}
    t = tables(data)
    assert list(t) == ["nucleotides", "base_pairs"]
    columns, rows = t["nucleotides"]
    assert len(rows) == na.size and [r[0] for r in rows] == list(range(1, na.size + 1))
    col = {c: k for k, c in enumerate(columns)}
    partner = {r[0]: r[col["partner"]] for r in rows}
    assert {(i, j) for i, j in partner.items() if j and i < j} == {
        (p["i"], p["j"]) for p in data["nucleic_acid"]["pairs"]
    }
    assert "".join(r[col["bracket"]] for r in rows) == data["nucleic_acid"]["strands"][0]["dot_bracket"]
    assert len(t["base_pairs"][1]) == len(na.base_pairs)
    blob = zip_bytes(data)
    assert blob == zip_bytes(data)
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        assert zf.namelist() == [
            "nucleotides.csv",
            "base_pairs.csv",
            "nucleic_acid.dbn",
            "nucleic_acid.bpseq",
            "nucleic_acid.ct",
        ]
        assert zf.read("nucleic_acid.ct").decode() == ct_text(data)
        assert zf.read("nucleotides.csv").decode() == csv_text(*t["nucleotides"])


def test_csv_cells_from_the_structure_never_start_a_spreadsheet_formula():
    """Section review 13b-13d, Minor 20: a crafted mmCIF's chain id or residue name starting =, +, - or @ is written
    with a ' before it in the CSV tables (a spreadsheet would run it); numbers and residue numbers (-7, -7A) are not,
    and the JSON keeps every value as it is."""
    rows = [["=HYPERLINK(1)", "+x", 5, "-7A"], ["@a", "-1", -3, "A"], ["\tz", "-1e5", 2.5, "-"]]
    assert csv_text(("chain", "name", "n", "tag"), rows) == (
        "chain,name,n,tag\n'=HYPERLINK(1),'+x,5,-7A\n'@a,-1,-3,A\n'\tz,-1e5,2.5,'-\n"
    )
