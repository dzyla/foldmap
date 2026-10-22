"""foldmap diff's engine (plan 13e-1): what changed between two structures of one protein chain."""

import dataclasses
import os
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
import pytest

from foldmap.diff import CONTACT_KEYS, Diff, _side, counts, diff_sides, diff_structures

DATA = Path(__file__).parent / "data"
UBQ, HHB = DATA / "1UBQ.cif", DATA / "2HHB.cif"


def _turn(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _write(st: gemmi.Structure, path: Path, turn=None) -> Path:
    """st written as PDB to path, every atom turned by `turn` (a 3x3 rotation) and moved, when given."""
    if turn is not None:
        for model in st:
            for chain in model:
                for res in chain:
                    for atom in res:
                        p = turn @ np.array(atom.pos.tolist()) + [12.5, -7.25, 30.0]
                        atom.pos = gemmi.Position(*p)
    st.write_pdb(str(path))
    return path


@pytest.fixture(scope="module")
def self_diff() -> Diff:
    return diff_structures(UBQ, "A", UBQ, "A", domains={"N": "1-40", "C": "res:A:41-76"})


def test_a_structure_against_itself_has_no_differences(self_diff):
    d = self_diff
    assert d.alignment.method == "sequence"
    assert d.only_a == () and d.only_b == ()
    assert d.ss_changes == () and d.pairing == ()
    assert set(d.contacts) == set(CONTACT_KEYS) and all(v == () for v in d.contacts.values())
    assert len(d.displacement) == 76 and set(d.displacement.values()) == {0.0}
    assert d.domain_turns == {"N": 0.0, "C": 0.0}
    assert all(v == 0 for v in counts(d).values())
    assert len(d.mates.mate) == len(_side(UBQ, "A", "A", "auto", True).sses)


def test_a_turned_copy_has_no_differences_either(self_diff, tmp_path):
    path = _write(gemmi.read_structure(str(UBQ)), tmp_path / "turned.pdb", _turn([1.0, 2.0, -0.5], 1.1))
    d = diff_structures(UBQ, "A", path, "A", domains={"N": "1-40", "C": "res:A:41-76"})
    assert counts(d) == counts(self_diff)
    assert max(d.displacement.values()) < 1e-2  # the file's coordinates are written to 0.001 Å
    assert max(d.domain_turns.values()) < 0.1
    assert d.core == self_diff.core


def _no_helix(tmp_path) -> tuple[Path, list]:
    """1UBQ less its alpha helix (the longest H element), written to tmp_path; that helix's residue numbers."""
    side = _side(UBQ, "A", "A", "auto", True)
    helix = max((s for s in side.sses if s.kind == "H"), key=len)
    gone = [side.bb.labels[k] for k in range(helix.start, helix.end + 1)]
    st = gemmi.read_structure(str(UBQ))
    chain = st[0]["A"]
    for k in reversed(range(len(chain))):
        if (chain[k].seqid.num, chain[k].seqid.icode.strip()) in {(l.seq, l.icode) for l in gone}:
            del chain[k]
    return _write(st, tmp_path / "ubq_no_helix.pdb"), [helix.id, *gone]


@pytest.mark.parametrize("match", ["auto", "structure"])
def test_a_deleted_helix_is_exactly_that_difference(tmp_path, match):
    path, (helix, *gone) = _no_helix(tmp_path)
    d = diff_structures(UBQ, "A", path, "A", match=match)
    assert d.alignment.method == ("sequence" if match == "auto" else "structure")
    assert d.only_a == (helix,) and d.only_b == ()
    assert d.ss_changes == tuple(("A", l.tag, "H", "-") for l in gone)
    assert d.pairing == () and all(v == () for v in d.contacts.values())
    assert len(d.displacement) == 76 - len(gone) and set(d.displacement.values()) == {0.0}


_PRINT = """
import sys
from foldmap.diff import diff_structures

d = diff_structures(sys.argv[1], "A", sys.argv[2], "A", domains={"N": "1-40"})
for f in ("only_a", "only_b", "ss_changes", "pairing", "core", "domain_turns", "contacts"):
    print(f, repr(getattr(d, f)))
print(repr(sorted(d.displacement.items())))
"""


def test_deterministic_across_hash_seeds(tmp_path):
    path, _ = _no_helix(tmp_path)
    turned = _write(gemmi.read_structure(str(path)), tmp_path / "turned_no_helix.pdb", _turn([0.3, -1.0, 2.0], 2.2))
    runs = set()
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _PRINT, str(UBQ), str(turned)],
            env=dict(os.environ, PYTHONHASHSEED=seed),
            capture_output=True,
            text=True,
            check=True,
        )
        runs.add(done.stdout)
    assert len(runs) == 1
    assert "only_a ('A:" in runs.pop()


def test_pairing_lost_and_register_moved():
    """B's bridges edited: one strand pairing's bridges dropped (lost), another's partner residues one further on
    (register +1, antiparallel: i + j grows by one)."""
    a = _side(UBQ, "A", "A", "auto", True)
    strands = [s for s in a.sses if s.kind == "E"]
    from foldmap.diff import _pairings

    pairs = _pairings(a)  # N to C
    assert len(pairs) >= 2
    (x1, y1), (x2, y2) = list(pairs)[:2]
    by_id = {s.id: s for s in strands}

    def between(br, x, y):  # bridges have i < j: i on the earlier strand x
        sx, sy = by_id[x], by_id[y]
        return sx.start <= br.i <= sx.end and sy.start <= br.j <= sy.end

    bridges = []
    for br in a.bridges:
        if between(br, x1, y1):
            continue  # pairing 1 lost
        if between(br, x2, y2):  # pairing 2: the later strand's residue one on (j - i or j + i: one more)
            br = br._replace(j=br.j + 1)
        bridges.append(br)
    d = diff_sides(a, dataclasses.replace(a, bridges=bridges))
    assert d.pairing == ((x1, y1, "lost"), (x2, y2, "register +1"))
    assert d.only_a == d.only_b == () and d.ss_changes == ()


def test_domains_are_checked(self_diff):
    with pytest.raises(ValueError, match="domain X: 'B:1-40' is not a residue range of A's chain A"):
        diff_structures(UBQ, "A", UBQ, "A", domains={"X": "B:1-40"})
    with pytest.raises(ValueError, match="domain X: '500-600' names no residue of chain A"):
        diff_structures(UBQ, "A", UBQ, "A", domains={"X": "500-600"})
    with pytest.raises(ValueError, match="no chain 'Q' in B"):
        diff_structures(UBQ, "A", UBQ, "Q")


def test_interface_and_ligand_changes(tmp_path):
    """2HHB's chain A against a copy without chain B and without its haems: the interface with B lost, the haem
    lost; the interfaces with C and D, and chain A itself, unchanged."""
    st = gemmi.read_structure(str(HHB))
    st.setup_entities()
    st[0].remove_chain("B")
    for chain in st[0]:
        for k in reversed(range(len(chain))):
            if chain[k].name == "HEM":
                del chain[k]
    path = _write(st, tmp_path / "hhb_less.pdb")
    same = diff_structures(HHB, "A", HHB, "A", assembly="asu")
    assert all(v == () for v in same.contacts.values())
    d = diff_structures(HHB, "A", path, "A", assembly="asu")
    assert d.only_a == d.only_b == () and d.ss_changes == () and d.pairing == ()
    assert d.contacts["interfaces"] == (("B", "ca", "lost"),)
    assert d.contacts["ligands"] == (("HEM", "lost"),)
    assert d.contacts["elements"] == d.contacts["interface_elements"] == d.contacts["ligand_elements"] == ()


def test_two_proteins_take_the_structure_route():
    """1UBQ against T4 lysozyme: not one protein by A19's rule, so mapped by structure; every element of each that
    has no mate is named, and nothing fails on a small core."""
    d = diff_structures(UBQ, "A", DATA / "2LZM.cif", "A")
    assert d.alignment.method == "structure"
    assert len(d.only_a) + len(d.mates.mate) + len(d.mates.pieces) == len(_side(UBQ, "A", "A", "auto", True).sses)
    assert len(d.core) >= 3 and len(d.displacement) == len(d.alignment.pairs)


@pytest.mark.parametrize("angle", [45.0, 60.0])
def test_a_turned_domain_shows_its_turn(tmp_path, angle):
    """1UBQ with residues 41-76 turned about residue 40's CA (a hinge): the core fit, started from the TM-score's
    superposition, keeps residues 1-40 (the turned ones move past d0 and are dropped), so N turns by nothing and C by
    the whole angle. At 45 degrees the cut first planned, 2·d0, kept 4 of C's residues and split the turn 7.7 / 41.1
    (plan 13e-3)."""
    st = gemmi.read_structure(str(UBQ))
    moved = [a for r in st[0]["A"] if 41 <= r.seqid.num <= 76 for a in r]
    centre = np.array(st[0]["A"]["40"][0]["CA"][0].pos.tolist())
    turn = _turn([0.0, 1.0, 1.0], np.radians(angle))
    for a in moved:
        a.pos = gemmi.Position(*(turn @ (np.array(a.pos.tolist()) - centre) + centre))
    d = diff_structures(UBQ, "A", _write(st, tmp_path / "hinge.pdb"), "A", domains={"N": "1-40", "C": "41-76"})
    assert d.domain_turns["N"] < 0.1
    assert d.domain_turns["C"] == pytest.approx(angle, abs=0.1)
    assert all(i < 40 for i, _ in d.core)


# ------------------------------------------------------------------------------- the diff page (plan 13e-2)

PAIRS = DATA / "diff"  # two states of one protein each: adenylate kinase open (4AKE, apo) and closed (1AKE, AP5A);
AKE = (PAIRS / "4AKE.cif.gz", PAIRS / "1AKE.cif.gz")  # calmodulin with calcium (1CLL) and without (1CFD, NMR)
CAM = (PAIRS / "1CLL.cif.gz", PAIRS / "1CFD.cif.gz")


def _gids(fig, prefix: str) -> list[str]:
    return sorted(
        a.get_gid()[len(prefix) :]
        for ax in fig.axes
        for a in ax.findobj(lambda o, ax=ax: o is not ax and (o.get_gid() or "").startswith(prefix))
    )


def test_a_self_diff_page_is_its_compare_page(tmp_path):
    """1UBQ against itself: nothing changed, so nothing is drawn or said beyond the compare page: the same SVG,
    byte for byte, and its JSON says nothing changed."""
    from foldmap.cli import save
    from foldmap.compare import DIFF_KEYS, compare, diff_page

    made, plain = diff_page(UBQ, UBQ), compare(UBQ, UBQ)
    assert made.note == plain.note and made.provenance == plain.provenance
    svg, svg_plain = tmp_path / "diff.svg", tmp_path / "compare.svg"
    save(made.figure, svg)
    save(plain.figure, svg_plain)
    assert svg.read_bytes() == svg_plain.read_bytes()
    assert not [g for p in ("a/absent:", "b/absent:") for g in _gids(made.figure, p)]
    data = made.data
    assert {k: data[k] for k in plain.data} == plain.data
    assert data["ss_changes"] == data["pairing"] == [] and data["domain_turns"] == {}
    assert all(v == [] for v in data["contacts"].values()) and {r[2] for r in data["displacement"]} == {0.0}
    assert set(data) == set(plain.data) | set(DIFF_KEYS)


def test_the_deleted_helix_is_drawn_dashed_in_bs_panel(tmp_path):
    """1UBQ less its helix: the helix outlined dashed in B's panel, where A's panel draws it (the panels are alike:
    the mates' rigid fit carries it), with the legend's entry; A's panel has no outline (B has nothing A lacks);
    nothing moved, so B keeps its colours; the note names the helix."""
    from foldmap.compare import diff_page
    from foldmap.figcheck import fit_deviations
    from foldmap.render import ABSENT_LEGEND, absent_outline

    path, (helix, *_) = _no_helix(tmp_path)
    made = diff_page(UBQ, path)
    assert made.diff.only_a == (helix,)
    (lay_a, _, _), (lay_b, _, _) = made.panels
    assert _gids(made.figure, "b/absent:") == [helix] and _gids(made.figure, "a/absent:") == []
    (out,) = lay_b.absent
    p = lay_a.placed[helix]
    assert (out.length, out.width, out.sse.kind) == (p.length, p.width, "H")
    mates = [(a, b) for b, a in made.diff.mates.mate.items() if a in lay_a.placed and b in lay_b.placed]
    at_a = [(lay_a.placed[a].cx, lay_a.placed[a].cy) for a, _ in mates] + [(p.cx, p.cy)]
    at_b = [(lay_b.placed[b].cx, lay_b.placed[b].cy) for _, b in mates] + [(out.cx, out.cy)]
    assert fit_deviations(at_a, at_b)[-1] < 1.0  # where the mates' fit puts A's helix on B's page
    (line,) = [a for a in made.figure.axes[2].lines if a.get_gid() == f"b/absent:{helix}"]
    assert line.get_linestyle() == "--" and len(absent_outline(out)) > 8
    legend = [t.get_text() for t in made.figure.axes[0].texts]
    assert ABSENT_LEGEND in legend and not any(t.startswith("Cα displacement") for t in legend)
    assert lay_b.residue_data is None
    assert f"only in 1UBQ.cif: {p.label} (dashed outlines on the right)" in made.note
    assert f"{len(made.diff.ss_changes)} residues change secondary structure" in made.note


def test_diff_json_schema(tmp_path):
    """--json: compare's matching, then ss_changes, pairing, displacement, core, domain_turns and contacts; numbers to
    1e-3; written with keys sorted, short rows on one line, and read back the same."""
    import json

    from foldmap.compare import DIFF_KEYS, diff_page, write_json

    path, (helix, *gone) = _no_helix(tmp_path)
    data = diff_page(UBQ, path, domains={"N": ["1-40"], "C": ["41-76"]}).data
    assert set(data) == {"alignment", "mates", "pieces", "extras", "only_a", *DIFF_KEYS}
    assert [e["a_residues"] for e in data["only_a"]] == [f"res:A:{gone[0].tag}-{gone[-1].tag}"]
    assert data["ss_changes"] == [{"chain": "A", "residue": l.tag, "a": "H", "b": "-"} for l in gone]
    assert data["pairing"] == [] and set(data["domain_turns"]) == {"N", "C"}
    assert [r[1] for r in data["displacement"]] == [p[1] for p in data["alignment"]["pairs"]]
    assert all(len(r) == 3 and r[0] == "A" and r[2] == round(r[2], 3) for r in data["displacement"])
    assert isinstance(data["core"], int) and data["core"] >= 3
    assert set(data["contacts"]) == set(CONTACT_KEYS)
    out = write_json(data, tmp_path / "diff.json")
    text = out.read_text(encoding="utf-8")
    assert json.loads(text) == data
    assert '[\n  "A",\n  "1",\n' not in text and '"displacement": [\n  ["A", "1", 0.0],' in text


@pytest.fixture(scope="module", params=["ake", "cam"])
def two_states(request):
    from foldmap.compare import diff_page

    a, b = AKE if request.param == "ake" else CAM
    domains = {"N": ["5-75"], "C": ["82-146"]} if request.param == "cam" else {"LID": ["122-159"], "NMP": ["30-59"]}
    return request.param, diff_page(a, b, domains=domains)


def test_two_states_are_drawn_alike_with_their_changes(two_states):
    """Adenylate kinase open and closed, calmodulin with and without calcium: one protein each (sequence route),
    B drawn like A even where the domains moved past the fold's TM-score (calmodulin, 0.38); each element only one
    state has outlined in the other's panel; B coloured by displacement, with its colour bar; the lid and the
    C-terminal lobe turn."""
    name, made = two_states
    d = made.diff
    assert d.alignment.method == "sequence" and made.note.startswith("drawn like ")
    for side, ids in (("a/absent:", d.only_b), ("b/absent:", d.only_a)):
        assert _gids(made.figure, side) == sorted(ids)
    assert d.only_a if name == "ake" else d.only_b
    legend = [t.get_text() for t in made.figure.axes[0].texts]
    assert any(t.startswith("Cα displacement (Å): ") for t in legend)
    assert _gids(made.figure, "b/data-stripe:") and not _gids(made.figure, "a/data-stripe:")
    turns = d.domain_turns
    assert (turns["LID"] > 30 and turns["NMP"] > 20) if name == "ake" else (turns["C"] > 90 and turns["N"] < 20)
    assert "the right panel coloured by Cα displacement" in made.note


def test_each_panel_stays_a_good_figure(two_states):
    """H5 for each panel: drawn on its own with its marks and colouring, it passes every S- and A-check its chain
    drawn alone passes, A6 (no text crossed, the outlines read as lines) among them."""
    from foldmap.compare import panel_checks

    _, made = two_states
    got = panel_checks(made)
    assert set(got) == {"a", "b"}
    for side in got.values():
        assert side["pass"] is True, side["detail"]
        assert "A6" in side["value"]["alone"]


_DIFF_PAGE = """
import hashlib, sys
from foldmap import cli
svg, js = sys.argv[3], sys.argv[4]
assert cli.main(["diff", sys.argv[1], sys.argv[2], "-o", svg, "--json", js, "--domain", "LID=122-159"]) == 0
print(hashlib.sha256(open(svg, "rb").read() + open(js, "rb").read()).hexdigest())
"""


def test_the_diff_page_does_not_depend_on_the_hash_seed(tmp_path):
    """A12: the diff page (outlines and displacement colouring) and its JSON have the same bytes under
    PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        done = subprocess.run(
            [
                sys.executable,
                "-c",
                _DIFF_PAGE,
                *map(str, AKE),
                str(tmp_path / f"{seed}.svg"),
                str(tmp_path / f"{seed}.json"),
            ],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout.split()[-1])
    assert len(runs[0]) == 64 and runs[0] == runs[1]


def test_command_line(tmp_path, capsys):
    """foldmap diff A B -o diff.svg --json diff.json --check: the page and the JSON written, the note on stderr, H5
    of each panel printed; a malformed --domain is a usage error (exit 2)."""
    import json

    from foldmap import cli

    svg, js = tmp_path / "diff.svg", tmp_path / "diff.json"
    args = ["diff", *map(str, AKE), "-o", str(svg), "--json", str(js), "--domain", "LID=122-159", "--check"]
    assert cli.main(args) == 0
    out, err = capsys.readouterr()
    assert "changes from 4AKE.cif.gz to 1AKE.cif.gz: only in 4AKE.cif.gz:" in err and "LID " in err
    assert "H5 still a good figure, left panel: pass" in out and "H5 still a good figure, right panel: pass" in out
    assert 'id="b/absent:' in svg.read_text(encoding="utf-8")
    assert json.loads(js.read_text(encoding="utf-8"))["domain_turns"]["LID"] > 30
    assert cli.main(["diff", str(UBQ), str(UBQ), "-o", str(svg), "--domain", "LID"]) == 2
    assert "give NAME=RANGES" in capsys.readouterr().err
