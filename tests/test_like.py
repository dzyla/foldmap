"""A layout drawn like a reference (plan 13a-5): frame_like, Like and the homes build_layout takes from it."""

import functools
import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from foldmap import cli, figcheck
from foldmap.align import Mates, align_chains, element_mates
from foldmap.dssp import assign_dssp
from foldmap.frame import CylinderFrame, Frame, frame_like
from foldmap.layout import LIKE_WEIGHT
from foldmap.like import Like, like_from
from foldmap.sheets import build_sheets
from foldmap.style import Style

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).resolve().parent.parent


def _bench():
    spec = importlib.util.spec_from_file_location("foldmap_bench", ROOT / "scripts" / "bench.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BENCH = _bench()
signature = BENCH.signature
OPTIONS = BENCH.OPTIONS  # every make_layout option the bench figures depend on


def _turn(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def laid_out_like(path, like: Like | None, frame, **kw):
    """cli.make_layout of `path` with build_layout given `like` and the page frame `frame` (13a-6 gives make_layout
    these itself)."""
    with (
        mock.patch.object(cli, "build_layout", functools.partial(cli.build_layout, like=like)),
        mock.patch.object(cli, "view_frame", lambda *a, **k: frame),
    ):
        return cli.make_layout(path, **{**OPTIONS, **kw})


def drawn_like(path_a, chain_a, path_b, chain_b, look=None, mates=None):
    """B (one chain) drawn like A (one chain): A laid out, the chains aligned, the elements mated (or `mates`), the
    Like built from A's finished layout, B laid out through A's frame carried over. Returns (layout B, sses B, bb B,
    layout A, Like)."""
    look = look or Style()
    lay_a, sses_a, bb_a = cli.make_layout(path_a, look=look, chains=[chain_a], **OPTIONS)
    _, sses_b, bb_b = cli.make_layout(path_b, look=look, chains=[chain_b], **OPTIONS)
    aln = align_chains(bb_a, chain_a, bb_b, chain_b)
    mates = mates or element_mates(sses_a, sses_b, aln)
    like = like_from(lay_a, sses_a, sses_b, mates, look, build_sheets(sses_b, assign_dssp(bb_b).bridges))
    frame = frame_like(lay_a.frame, aln.rotation, aln.translation)
    lay_b, sses_b, bb_b = laid_out_like(path_b, like, frame, look=look, chains=[chain_b])
    assert lay_b.frame is frame
    return lay_b, sses_b, bb_b, lay_a, like


def test_frame_like_carries_known_rotation():
    rot_a = _turn([0.3, -1.0, 0.6], 0.8)  # A's page axes: the rows of a proper rotation
    frame_a = Frame(np.array([4.0, -2.5, 11.0]), rot_a[0], rot_a[1], rot_a[2])
    r, t = _turn([1.0, 2.0, -0.5], 1.1), np.array([12.5, -7.25, 30.0])
    b = np.array([[np.sin(3.0 * k), np.cos(2.0 * k) * 7.0, 0.5 * k - 9.0] for k in range(40)])  # B's points
    a = b @ r.T + t  # their places in A: A ≈ R·B + t
    frame_b = frame_like(frame_a, r, t)
    assert np.abs(frame_b.project(b) - frame_a.project(a)).max() < 1e-9
    assert np.abs(frame_b.depth(b) - frame_a.depth(a)).max() < 1e-9
    assert np.linalg.det(np.array([frame_b.u, frame_b.v, frame_b.w])) == pytest.approx(1.0)  # never a mirror
    mirrored = Frame(frame_a.origin, frame_a.u, frame_a.v, frame_a.w, sector={"A": 0, "B": 1}, mirror_x=0.0)
    unrolled = CylinderFrame(np.zeros(3), np.array([0, 0, 1.0]), np.array([1.0, 0, 0]), 3, 10.0, {"A": 0})
    for other in (mirrored, unrolled):
        with pytest.raises(TypeError, match="plain Frame"):
            frame_like(other, r, t)


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "1TIM"])  # 1TIM: βα rows (the reference's row way up)
def test_like_itself_reproduces_the_layout(name):
    """The invariant that tests the design (spec 4.3): a structure drawn like itself is its own layout, exactly."""
    path = DATA / f"{name}.cif"
    plain, sses, _ = cli.make_layout(path, chains=["A"], **OPTIONS)
    lay, _, _, _, like = drawn_like(path, "A", path, "A")
    assert set(like.target) == {s.id for s in sses}  # every element held at its own place
    assert like.weight == LIKE_WEIGHT == 8
    if name == "1TIM":
        assert plain.ba_rows and like.template and any(like.template)
    assert signature(lay) == signature(plain)
    # the places are the reference's, not the projection's: the reference moved (its targets, not their units, which
    # then no longer agree and are not read), the layout moves with it, to the grid
    moved = replace(like, target={k: (x + 3.0, y - 2.0) for k, (x, y) in like.target.items()})
    shifted, _, _ = laid_out_like(path, moved, lay.frame, chains=["A"])
    for k, p in plain.placed.items():
        assert (shifted.placed[k].cx, shifted.placed[k].cy) == pytest.approx((p.cx + 3.0, p.cy - 2.0), abs=1e-4)


BEFORE = json.loads((DATA / "like" / "signatures-before.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("name", ["1UBQ", "2LZM", "1TIM"])
def test_like_none_changes_nothing(name):
    """Off means identical: with like=None (the default) the three figures are the ones recorded before 13a-5
    (bench signatures under bench's options, tests/data/like/signatures-before.json)."""
    assert BEFORE["_options"] == OPTIONS
    path = DATA / f"{name}.cif"
    with mock.patch.object(cli, "build_layout", functools.partial(cli.build_layout, like=None)):
        given, _, _ = cli.make_layout(path, **OPTIONS)
    assert json.loads(json.dumps(signature(given))) == BEFORE[name]


def _strand_order(lay, ids):
    """The strands of a sheet block left to right as drawn."""
    return [k for k in sorted(ids, key=lambda k: lay.placed[k].cx)]


def test_sheet_order_kept():
    """A B whose strands' mates stand in another order in A (here 1UBQ's strands mated in reverse) keeps the strand
    order of its own pairing: Like moves the sheet as a whole, never a strand out of its place in it."""
    path = DATA / "1UBQ.cif"
    alone, sses, bb = cli.make_layout(path, chains=["A"], **OPTIONS)
    strands = [s.id for s in sses if s.kind == "E"]
    helices = [s.id for s in sses if s.kind != "E"]
    mate = {**dict(zip(strands, reversed(strands))), **{h: h for h in helices}}
    mates = Mates(mate=mate, pieces={}, extras=(), only_a=(), same_fold=True)
    lay, sses_b, bb_b, lay_a, like = drawn_like(path, "A", path, "A", mates=mates)
    assert like.target[strands[0]] == (lay_a.placed[strands[-1]].cx, lay_a.placed[strands[-1]].cy)
    assert len(lay.sheets) == len(alone.sheets) == 1
    mine, own = _strand_order(lay, lay.sheets[0]), _strand_order(alone, alone.sheets[0])
    assert mine in (own, own[::-1])
    # turned half round in the page (order reversed and every strand turned), never seen in a mirror
    turned = mine != own
    for k in own:
        up_alone, up_here = np.sin(alone.placed[k].angle) > 0, np.sin(lay.placed[k].angle) > 0
        assert up_here == (not up_alone if turned else up_alone), k
    page = SimpleNamespace(axes=[SimpleNamespace(patches=[])])
    found = figcheck.science_checks(SimpleNamespace(layout=lay, sses=sses_b, bb=bb_b, figure=page))
    assert found["S1"] == found["S2"] == found["S3"] == [] and found["S6"] == []


def test_homologue_takes_its_mates_places_names_and_angles():
    """Haemoglobin α (2HHB chain A) drawn like β (chain B): every mated helix carries its mate's label and drawn
    angle, and sits nearer its mate's place than drawn alone does (in the same frame)."""
    path = DATA / "2HHB.cif"
    lay, sses_b, _, lay_a, like = drawn_like(path, "B", path, "A")
    assert len(like.target) >= 6
    alone, _, _ = laid_out_like(path, None, lay.frame, chains=["A"])
    near, far = [], []
    for s in sses_b:
        if s.id not in like.target:
            continue
        p = lay.placed[s.id]
        assert p.angle == like.angle[s.id]
        assert p.label == like.label[f"res:{s.chain}:{s.first.tag}-{s.last.tag}"]
        tx, ty = like.target[s.id]
        near.append(np.hypot(p.cx - tx, p.cy - ty))
        far.append(np.hypot(alone.placed[s.id].cx - tx, alone.placed[s.id].cy - ty))
    assert np.median(near) < np.median(far)


def test_split_helix_parts_and_relations_only_b_has():
    """2HHB α drawn like β: α's helix split in two is laid end to end on its mate's place, the first part in
    sequence named as the mate and the next primed; neighbours in α whose mates in β are not neighbours (β's D helix
    lies between) are no relation the reference settled (Like.settled, which a disulfide only α had would be judged
    by); in α drawn like itself every neighbour pair is settled. And α sits at β's places: the H2 starting limits
    (spec 4.5: median <= 1.0, 90th percentile <= 2.5 page units from the targets)."""
    path = DATA / "2HHB.cif"
    lay, sses_b, _, lay_a, like = drawn_like(path, "B", path, "A")
    _, sses_a, bb_a = cli.make_layout(path, chains=["B"], **OPTIONS)
    _, _, bb_b = cli.make_layout(path, chains=["A"], **OPTIONS)
    mates = element_mates(sses_a, sses_b, align_chains(bb_a, "B", bb_b, "A"))
    assert mates.pieces
    by_id = {s.id: s for s in sses_b}
    for a in sorted(set(mates.pieces.values())):
        parts = sorted(
            [b for b, x in [*mates.mate.items(), *mates.pieces.items()] if x == a], key=lambda b: by_id[b].start
        )
        names = [lay.placed[b].label for b in parts]
        assert names[0] == lay_a.placed[a].label and names[1] == names[0] + "′"
        for b in parts:  # its places on the grid (cos and sin of the mate's angle stay off the optimiser)
            assert like.target[b] == tuple(np.round(np.array(like.target[b]) * 1e4) / 1e4)
    seq = sorted(sses_b, key=lambda s: s.start)
    neighbours = {frozenset((x.id, y.id)) for x, y in zip(seq, seq[1:]) if x.id in like.target and y.id in like.target}
    assert neighbours - like.settled  # a relation only α has
    itself = drawn_like(path, "A", path, "A")[4]
    seq_a = sorted(cli.make_layout(path, chains=["A"], **OPTIONS)[1], key=lambda s: s.start)
    assert {frozenset((x.id, y.id)) for x, y in zip(seq_a, seq_a[1:])} <= itself.settled
    off = [np.hypot(lay.placed[k].cx - x, lay.placed[k].cy - y) for k, (x, y) in like.target.items()]
    assert np.median(off) <= 1.0 and np.percentile(off, 90) <= 2.5


def _two_held(bridged: bool, settled):
    """Two helices of one chain, consecutive in sequence (a loop), held 11.5 page units apart at their places
    (_attract's `held`), optionally bridged by a disulfide: their centres' distance after _attract."""
    from foldmap.layout import _attract, _helix_item
    from helpers import fake_sse

    frame = Frame(np.zeros(3), np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0]))
    a = fake_sse("H", 0, 9, centroid=(0, 0, 0), axis=(0, 1, 0))
    b = fake_sse("H", 14, 23, centroid=(50, 0, 0), axis=(0, 1, 0))
    items = [_helix_item(a, frame), _helix_item(b, frame)]
    held = {k: (it.cx, it.cy) for k, it in enumerate(items)}
    bridges = [(a.id, 0.5, b.id, 0.5)] if bridged else []
    _attract(items, [a, b], {}, bridges=bridges, held=held, settled=settled(a, b))
    return float(np.hypot(items[1].cx - items[0].cx, items[1].cy - items[0].cy))


def test_only_a_disulfide_this_structure_has_pulls_held_items():
    """Between two items held at the reference's places a loop never pulls (the reference places them, a loop only
    this structure has included), nor a disulfide the reference has too; a disulfide only this structure has does
    (a bond drawn across the page is wrong whatever the reference)."""
    apart = 50 * 0.23  # layout.SCALE
    none = lambda a, b: frozenset()  # noqa: E731 - nothing settled: the loop and any bridge are this structure's own
    both = lambda a, b: frozenset({frozenset((a.id, b.id))})  # noqa: E731 - the reference bridges their mates
    assert _two_held(False, none) == pytest.approx(apart, abs=1e-12)  # the loop: no pull
    assert _two_held(True, both) == pytest.approx(apart, abs=1e-12)  # the reference's own bridge: no pull
    assert _two_held(True, none) < apart - 1.0  # a bridge only this structure has: pulled together


def test_targets_a_bit_off_give_the_same_figure():
    """Every target moved by one unit in the last place (as sin and cos may differ between CPUs) draws the same
    figure: units read where they agree to 1e-9, targets on the grid otherwise (A12 on any CPU)."""
    path = DATA / "2HHB.cif"
    lay, _, _, _, like = drawn_like(path, "B", path, "A")
    off = replace(
        like, target={k: (np.nextafter(x, np.inf), np.nextafter(y, -np.inf)) for k, (x, y) in like.target.items()}
    )
    again, _, _ = laid_out_like(path, off, lay.frame, chains=["A"])
    assert signature(again) == signature(lay)


_SIGNATURE = """
import json, sys
sys.path.insert(0, sys.argv[1])
from test_like import DATA, drawn_like, signature
lay = drawn_like(DATA / "2HHB.cif", "B", DATA / "2HHB.cif", "A")[0]
print(json.dumps(signature(lay)))
"""


def test_like_figure_does_not_depend_on_the_hash_seed():
    """A12 for a figure drawn like a reference: the same layout under PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        done = subprocess.run(
            [sys.executable, "-c", _SIGNATURE, str(Path(__file__).parent)],
            env=env,
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout)
    assert json.loads(runs[0]) and runs[0] == runs[1]


# --like on the command line, names, colours and the note (plan 13a-6) ---------------------------------------------


def like_of(path, chain=None, match="auto"):
    """make_layout's `like` (FigureSpec.like, --like): the reference `path`, its chain `chain`."""
    return {"structure": str(path), "chain": chain, "match": match, "layout": None}


def test_note_on_itself():
    """1UBQ drawn like itself: the note, every element matched, and the figure its own (through make_layout)."""
    path = DATA / "1UBQ.cif"
    lay, _, _ = cli.make_layout(path, like=like_of(path), **OPTIONS)
    plain, _, _ = cli.make_layout(path, chains=["A"], **OPTIONS)
    assert lay.like_note == (
        "drawn like 1UBQ.cif chain A: 8 of 8 elements matched (sequence alignment, TM-score 1.00 over 76 residues)"
    )
    assert lay.like_extras == []
    assert signature(lay) == signature(plain)


def test_note_names_the_extras_as_drawn():
    """Haemoglobin β (2HHB chain B) drawn like α (chain A): of β's two elements α lacks, the α-helix is named after
    α's helix before it, lettered (α3a), and the 3-10 helix, with no 3-10 helix of α before it, keeps its own name
    (η1: never α3b, which would call a 3-10 helix an α-helix; 13a-7), in the figure and in the note."""
    path = DATA / "2HHB.cif"
    lay, _, _ = cli.make_layout(path, chains=["B"], like=like_of(path, "A"), **OPTIONS)
    assert lay.like_note == (
        "drawn like 2HHB.cif chain A: 7 of 9 elements matched (sequence alignment, TM-score 0.89 over 139 residues); "
        "only in this structure: η1, α3a"
    )
    assert [lay.placed[k].label for k in lay.like_extras] == ["η1", "α3a"]


def test_not_the_same_fold_is_drawn_on_its_own():
    """T4 lysozyme (2LZM) is not drawn like ubiquitin: the note says why, and the figure is 2LZM's own."""
    lay, _, _ = cli.make_layout(DATA / "2LZM.cif", like=like_of(DATA / "1UBQ.cif"), **OPTIONS)
    alone, _, _ = cli.make_layout(DATA / "2LZM.cif", chains=["A"], **OPTIONS)
    assert lay.like_note == "not drawn like 1UBQ.cif: not the same fold (TM-score 0.14; 4 of 13 elements matched)"
    assert lay.like_colour == {} and lay.like_extras == []
    assert signature(lay) == signature(alone)


def test_note_words():
    """The note's alignment words, and the low-coverage clause: how much of the shorter chain aligned, rounded down,
    with no guess at why (a hinge between domains, or a circular permutation; 13a-7)."""
    from foldmap.align import ChainAlignment
    from foldmap.like import like_note, refused_note

    aln = ChainAlignment(
        pairs=tuple((k, k) for k in range(50)),
        tm_score=0.414,
        rmsd=1.0,
        identity=0.2,
        method="structure",
        rotation=np.eye(3),
        translation=np.zeros(3),
        residues_a=(None,) * 100,
        residues_b=(None,) * 80,
    )
    mates = Mates(mate={"b1": "a1"}, pieces={"b2": "a1"}, extras=("b3",), only_a=(), same_fold=False)
    assert refused_note("A.cif", aln, mates) == (
        "not drawn like A.cif: not the same fold (TM-score 0.41; 2 of 3 elements matched); only 62% of the shorter "
        "chain aligned"
    )
    filed = replace(aln, method="file", pairs=tuple((k, k) for k in range(80)))
    assert like_note("A.cif", "B", filed, mates, ["α1a"], "x/pair.fasta", [("η3", "α3")]) == (
        "drawn like A.cif chain B: 2 of 3 elements matched (alignment pair.fasta, TM-score 0.41 over 80 residues); "
        "only in this structure: α1a; named by their own helix type: η3 for α3"
    )


def test_extras_named_after_the_element_before_them():
    """An extra takes the name of the matched element of its type before it: a 3-10 helix after η1 lettered (η1a,
    never after the α-helix between them: 13a-7), a strand named by a capital letter primed, as Ig's C′ (never
    'Ba')."""
    lay_a, sses, _ = cli.make_layout(DATA / "1UBQ.cif", chains=["A"], **OPTIONS)
    seq = sorted(sses, key=lambda s: s.start)
    helices = [s for s in seq if s.kind == "G"]  # η1, η2
    strands = [s for s in seq if s.kind == "E"]
    helix, strand = helices[1], strands[2]  # each with a matched element of its type before it
    extras = (helix.id, strand.id)
    mates = Mates(
        mate={s.id: s.id for s in sses if s.id not in extras},
        pieces={},
        extras=tuple(s.id for s in seq if s.id in extras),
        only_a=tuple(s.id for s in seq if s.id in extras),
        same_fold=True,
    )
    like = like_from(lay_a, sses, sses, mates, Style())
    named = {k: like.label[f"res:A:{s.first.tag}-{s.last.tag}"] for k, s in (("helix", helix), ("strand", strand))}
    assert named["helix"] == lay_a.placed[helices[0].id].label + "a"
    assert named["strand"] == lay_a.placed[strands[1].id].label + "′"
    assert named["strand"][0].isupper() and named["strand"][1:] == "′"


def _haemoglobins(colour: str):
    """β (2HHB B) drawn like α (A) and α drawn alone, coloured by `colour`: (fill of β as drawn, of α, of β alone,
    the mates)."""
    from foldmap.render import element_colours

    path, look = DATA / "2HHB.cif", Style(color_by=colour)
    lay_b, sses_b, bb_b = cli.make_layout(path, chains=["B"], look=look, like=like_of(path, "A"), **OPTIONS)
    lay_a, sses_a, bb_a = cli.make_layout(path, chains=["A"], look=look, **OPTIONS)
    alone, sses_alone, _ = cli.make_layout(path, chains=["B"], look=look, **OPTIONS)
    mates = element_mates(sses_a, sses_b, align_chains(bb_a, "A", bb_b, "B"))
    return (
        element_colours(lay_b, sses_b, look)[0],
        element_colours(lay_a, sses_a, look)[0],
        element_colours(alone, sses_alone, look)[0],
        mates,
    )


@pytest.mark.parametrize("colour", ["chain", "sequence"])
def test_identity_colouring_takes_the_reference_colours(colour):
    """Coloured by identity, a mate is drawn in its mate's colour and an extra in a tint of its neighbour's."""
    fill_b, fill_a, _, mates = _haemoglobins(colour)
    assert mates.mate and mates.extras
    for b, a in mates.mate.items():
        assert fill_b[b] == fill_a[a], b
    for e in mates.extras:
        assert fill_b[e] not in fill_a.values(), e


def test_data_colouring_keeps_its_own_values():
    """Coloured by data (B-factors), each structure keeps its own values: β drawn like α is coloured as β alone."""
    fill_b, fill_a, alone, mates = _haemoglobins("bfactor")
    assert fill_b == alone
    assert any(fill_b[b] != fill_a[a] for b, a in mates.mate.items())


def test_like_on_the_command_line(tmp_path, capsys):
    """`plot B --chains B --like A --like-chain A`: the figure is drawn and its notes (printed, and in the SVG's
    <desc>) say what it was drawn like."""
    path, out = DATA / "2HHB.cif", tmp_path / "beta.svg"
    argv = ["plot", str(path), "--chains", "B", "--like", str(path), "--like-chain", "A", "-o", str(out)]
    assert cli.main(argv) == 0
    assert "drawn like 2HHB.cif chain A: 7 of 9 elements matched" in capsys.readouterr().err
    assert "drawn like 2HHB.cif chain A: 7 of 9 elements matched" in out.read_text(encoding="utf-8")


def test_like_with_two_chains_is_a_usage_error(tmp_path, capsys):
    """One protein chain of each structure: --like with two chains named is a usage error (exit 2) naming
    --chains; so is a --like option without --like."""
    path, out = DATA / "2HHB.cif", str(tmp_path / "x.svg")
    assert cli.main(["plot", str(path), "--chains", "A,B", "--like", str(path), "-o", out]) == 2
    assert "--chains" in capsys.readouterr().err
    assert cli.main(["plot", str(path), "--like-chain", "A", "-o", out]) == 2
    assert "--like-chain needs --like" in capsys.readouterr().err


# the layout file, MCP and frames that cannot be carried (plan 13a-7) -----------------------------------------------


def test_a_mate_of_the_other_helix_type_keeps_its_type():
    """A mate of the other helix type takes its own type's sign with its mate's number (a 3-10 helix mated to α3:
    η3), never the mate's sign; where the reference names another element so (its α1), it is named as an extra is
    (lettered: α1a); Like.crossed records both for the note."""
    from foldmap.like import _labels

    _, sses, _ = cli.make_layout(DATA / "1UBQ.cif", **OPTIONS)
    by_label = dict(zip(["A", "B", "α1", "η1", "C", "D", "η2", "E"], sses))
    placed = {  # the reference: an α-helix α1, then α2 and α3, its own η1 elsewhere
        "x1": SimpleNamespace(label="α2"),
        "x2": SimpleNamespace(label="η1"),
        "x3": SimpleNamespace(label="α3"),
        "x4": SimpleNamespace(label="α1"),
    }
    mates = Mates(
        mate={by_label["α1"].id: "x1", by_label["η1"].id: "x2", by_label["η2"].id: "x3"},
        pieces={},
        extras=(),
        only_a=("x4",),
        same_fold=True,
    )
    ref = {k: f"res:A:{s.first.tag}-{s.last.tag}" for k, s in by_label.items()}
    labels, crossed = _labels(placed, sses, mates)
    assert labels[ref["α1"]] == "α2"  # the same type: its mate's name
    assert labels[ref["η1"]] == "η1"  # the same type
    assert labels[ref["η2"]] == "η3"  # a 3-10 helix mated to α3: η3, free
    assert crossed == {by_label["η2"].id: "α3"}
    swapped = replace(mates, mate={by_label["α1"].id: "x2", by_label["η1"].id: "x1"})
    labels, crossed = _labels(placed, sses, swapped)
    assert labels[ref["α1"]] == "α1a"  # an α-helix mated to η1: α1 is the reference's other helix, so lettered
    assert labels[ref["η1"]] == "η2"  # a 3-10 helix mated to α2: η2, free in the reference
    assert crossed == {by_label["α1"].id: "η1", by_label["η1"].id: "α2"}


UNCARRIED = {  # one chain whose frame still cannot be carried: drawn on its own, with the note
    "duplex": (
        "1LMB.cif",
        {},
        "not drawn like 1UBQ.cif: this structure's frame (a straight duplex) cannot be carried; draw it with --no-dna",
    ),
    "stack": (
        "1UBQ.cif",
        {"mode": "stack"},
        "not drawn like 1UBQ.cif: this structure's frame (--mode stack) cannot be carried; draw it with --mode "
        "projected",
    ),
}


@pytest.mark.parametrize("case", list(UNCARRIED))
def test_a_frame_that_cannot_be_carried_is_drawn_on_its_own(case):
    """One chain (its first protein chain, none named) whose frame is not a plain Frame (protein on a straight
    duplex, a stack) is not drawn like the reference: it is drawn exactly as that chain on its own, and the note says
    why and what would draw it like the reference."""
    name, options, note = UNCARRIED[case]
    kw = {**OPTIONS, **options}
    lay, _, bb = cli.make_layout(DATA / name, like=like_of(DATA / "1UBQ.cif"), **kw)
    alone, _, _ = cli.make_layout(DATA / name, chains=[bb.labels[0].chain], **kw)
    assert lay.like_note == note
    assert lay.like_colour == {} and lay.like_extras == []
    assert signature(lay) == signature(alone)


@pytest.mark.parametrize(
    "name, options",
    [("2HHB.cif", {}), ("6LTO.cif.gz", {"focus": "none", "look": Style(cyclic="rosette")})],
    ids=["unrolled", "rosette"],
)
def test_an_assembly_is_drawn_by_its_first_chain(name, options):
    """Ruling 7: an assembly (2HHB unrolled, 6LTO as a rosette) with no chain named on either side is drawn by its
    first protein chain like the reference's first: the note and the scope line say the chains were chosen, never
    asked for, and the figure is that chain drawn like the reference."""
    from foldmap.io import read_model
    from foldmap.scope import drawn_scope

    path, kw = DATA / name, {**OPTIONS, **options}
    lay, _, _ = cli.make_layout(path, like=like_of(path), **kw)
    named, _, _ = cli.make_layout(path, chains=["A"], like=like_of(path, "A"), **kw)
    assert lay.like_note.startswith(f"drawn like {name} chain A (its first protein chain): ")
    assert lay.like_chain == "A" and lay.chains_asked == [] and lay.uncarried is None
    reason = drawn_scope(read_model(path, "auto"), lay)["reason"]
    assert reason == "chain A (the first protein chain: one chain of each is drawn like a reference)"
    assert signature(lay) == signature(named)


def test_a_reference_frame_that_cannot_be_carried(tmp_path):
    """A reference laid out from its layout file as a stack: this figure is drawn on its own, the note naming the
    reference's frame."""
    path, saved = DATA / "1UBQ.cif", tmp_path / "ubq-stack.yaml"
    argv = ["plot", str(path), "--mode", "stack", "-o", str(tmp_path / "s.svg"), "--save-layout", str(saved)]
    assert cli.main(argv) == 0
    lay, _, _ = cli.make_layout(path, like={**like_of(path), "layout": str(saved)}, **OPTIONS)
    assert lay.like_note == (
        "not drawn like 1UBQ.cif: its frame (--mode stack) cannot be carried; draw it with --mode projected"
    )


def test_a_reference_drawn_as_an_assembly_takes_its_first_chain(tmp_path, capsys):
    """8UTF, drawn as an assembly by its symmetry, as the reference with no chain named: its first protein chain, the
    note saying so (ruling 7); with --like-chain A the note names the chain as asked."""
    path, out = DATA / "8UTF.cif", tmp_path / "x.svg"
    assert cli.main(["plot", str(path), "--chains", "B", "--like", str(path), "-o", str(out)]) == 0
    assert "drawn like 8UTF.cif chain A (its first protein chain): 2 of 2 elements matched" in capsys.readouterr().err
    assert cli.main(["plot", str(path), "--chains", "B", "--like", str(path), "--like-chain", "A", "-o", str(out)]) == 0
    assert "drawn like 8UTF.cif chain A: 2 of 2 elements matched" in capsys.readouterr().err


def test_a_chain_like_chose_is_not_said_as_asked_for():
    """Of several protein chains with none named, --like draws the first: the scope line says so, not "asked for";
    a structure with one protein chain says nothing of chains."""
    from foldmap.io import read_model
    from foldmap.scope import drawn_scope

    path = DATA / "2HHB.cif"
    lay, _, _ = cli.make_layout(path, like=like_of(path, "B"), **OPTIONS)
    assert lay.like_chain == "A" and lay.chains_asked == []
    assert lay.like_note.startswith("drawn like 2HHB.cif chain B: ")
    scope = drawn_scope(read_model(path, "auto"), lay)
    assert scope["chains_reason"] == "chain A (the first protein chain: one chain of each is drawn like a reference)"
    assert "asked for" not in scope["reason"] and "nucleotide" not in scope["reason"]
    ubq = DATA / "1UBQ.cif"
    one, _, _ = cli.make_layout(ubq, like=like_of(ubq), **OPTIONS)
    assert one.like_chain == "" and drawn_scope(read_model(ubq, "auto"), one)["reason"] == ""


_SVG = """
import hashlib, sys
from foldmap import cli
out = sys.argv[2]
assert cli.main(["plot", sys.argv[1], "--chains", "B", "--like", sys.argv[1], "--like-chain", "A", "-o", out]) == 0
print(hashlib.sha256(open(out, "rb").read()).hexdigest())
"""


def test_like_svg_does_not_depend_on_the_hash_seed(tmp_path):
    """A12 for the drawn figure: the SVG of β drawn like α has the same hash under PYTHONHASHSEED 0 and 1."""
    runs = []
    for seed in ("0", "1"):
        done = subprocess.run(
            [sys.executable, "-c", _SVG, str(DATA / "2HHB.cif"), str(tmp_path / f"{seed}.svg")],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        runs.append(done.stdout.split()[-1])
    assert len(runs[0]) == 64 and runs[0] == runs[1]
