"""The sequence figure: secondary structure drawn over the sequence (or an alignment), in rows of N columns."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap.cli import main
from foldmap.seqplot import draw_sequence
from foldmap.sequence import chain_track, column_classes, map_alignment, read_alignment
from foldmap.style import Style

DATA = Path(__file__).parent / "data"
UBQ = DATA / "1UBQ.cif"
FAMILY = DATA / "ubq_family.fasta"


def _artists(fig, prefix):
    ax = fig.axes[0]
    return [a for a in [*ax.patches, *ax.texts, *ax.lines] if (a.get_gid() or "").startswith(prefix)]


def _texts(fig):
    return [t.get_text() for t in fig.axes[0].texts]


def test_rows_of_n_columns():
    for columns, blocks in ((60, 2), (40, 2), (25, 4), (80, 1)):
        fig = draw_sequence(UBQ, columns=columns)
        assert len(_artists(fig, "block:")) == blocks


def test_every_residue_once_with_numbers_every_ten():
    fig = draw_sequence(UBQ, columns=30)
    letters = sorted(_artists(fig, "res:0:"), key=lambda t: int(t.get_gid().split(":")[2]))
    assert "".join(t.get_text() for t in letters) == chain_track(UBQ).letters
    assert {"10", "20", "30", "40", "50", "60", "70"} <= set(_texts(fig))


def test_elements_drawn_over_their_residues_with_figure_labels():
    track = chain_track(UBQ)
    fig = draw_sequence(UBQ, columns=80)
    for e in track.elements:
        pieces = _artists(fig, f"ss:{e.id}:")
        assert pieces, e.label
        x0 = min(p.get_path().vertices[:, 0].min() for p in pieces)
        x1 = max(p.get_path().vertices[:, 0].max() for p in pieces)
        assert x0 == pytest.approx(e.start - 0.5, abs=0.05) and x1 == pytest.approx(e.end + 0.5, abs=0.05)
    assert {e.label for e in track.elements} <= set(_texts(fig))
    assert "TT" in _texts(fig)


def test_elements_split_across_rows_keep_one_label():
    track = chain_track(UBQ)
    a1 = next(e for e in track.elements if e.label == "α1")
    split = a1.start + 2  # break the row inside α1
    fig = draw_sequence(UBQ, columns=split)
    assert len(_artists(fig, f"ss:{a1.id}:")) >= 2 and _texts(fig).count("α1") == 1


def test_glyphs_take_the_topology_colours():
    from foldmap.cli import make_layout
    from foldmap.render import element_colours

    look = Style(color_by="sequence")
    lay, sses, _ = make_layout(UBQ, look=look)
    colours, _ = element_colours(lay, sses, look)
    fig = draw_sequence(UBQ, look=look)
    for e in chain_track(UBQ).elements:
        piece = _artists(fig, f"ss:{e.id}:")[0]
        c = piece.get_facecolor() if e.kind == "E" else piece.get_edgecolor()
        assert to_hex(c) == to_hex(colours[e.id])
    plain = draw_sequence(UBQ, ss_colour="black")
    assert to_hex(_artists(plain, "ss:")[0].get_edgecolor()) == "#000000"


def test_unmodelled_residues_are_grey_lowercase():
    fig = draw_sequence(DATA / "6ZS5.cif", chain="A", full_sequence=True, columns=100)
    first = next(t for t in _artists(fig, "res:0:") if t.get_gid() == "res:0:0")
    assert first.get_text().islower() and to_hex(first.get_color()) != "#000000"


def test_alignment_shows_every_row_and_conservation():
    aln = read_alignment(FAMILY)
    fig = draw_sequence(UBQ, alignment=aln, columns=40)
    assert {n for n in aln.names} <= set(_texts(fig))
    for r in range(len(aln.rows)):
        assert len(_artists(fig, f"res:{r}:")) == aln.width
    strict, similar = column_classes(aln)
    assert len(_artists(fig, "strict:")) == int(strict.sum())
    red = {int(t.get_gid().split(":")[2]) for t in _artists(fig, "res:") if to_hex(t.get_color()) == "#c81e1e"}
    assert red == set(np.flatnonzero(similar))
    assert _artists(fig, "frame:") and len(_artists(fig, "cons:")) == aln.width


def test_alignment_structure_follows_the_reference_row():
    aln = read_alignment(FAMILY)
    track = chain_track(UBQ)
    m = map_alignment(aln, track)
    fig = draw_sequence(UBQ, alignment=aln, columns=100)
    col = {p: c for c, p in enumerate(m.position) if p is not None}
    for e in track.elements:
        pieces = _artists(fig, f"ss:{e.id}:")
        x1 = max(p.get_path().vertices[:, 0].max() for p in pieces)
        assert x1 == pytest.approx(col[e.end] + 0.5, abs=0.05)


def test_cli_sequence(tmp_path):
    out = tmp_path / "seq.svg"
    assert main(["sequence", str(UBQ), "-o", str(out), "--columns", "50"]) == 0
    assert "ss:" in out.read_text()
    out2 = tmp_path / "aln.png"
    assert main(["sequence", str(UBQ), "-o", str(out2), "--msa", str(FAMILY), "--reference", "1UBQ"]) == 0
    assert out2.stat().st_size > 1000
    assert main(["sequence", str(UBQ), "-o", str(out), "--chain", "Q"]) == 1
