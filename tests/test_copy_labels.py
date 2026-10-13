"""Copies of one protein share element labels (model.copy_classes): one element, one label, in every chain."""

from pathlib import Path

import pytest

from foldmap.cli import make_layout
from foldmap.figcheck import label_repeats as _repeats
from helpers import FID

DATA = Path(__file__).parent / "data"
FIXTURES = sorted([*DATA.glob("*.cif"), *FID.glob("*.cif.gz"), *(DATA / "na").glob("*.cif*")])


VIEWS = [  # (focus, names, domains): every focus with every naming, and the automatic domains (K6, K19)
    *((focus, names, None) for focus in ("auto", "none", "first") for names in ("none", "deposited", "beta")),
    ("auto", "none", "auto"),
]


@pytest.mark.parametrize("focus, names, domains", VIEWS, ids=lambda v: str(v))
@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.name.split(".")[0])
def test_no_chain_draws_one_label_twice(path, focus, names, domains):
    """Every fixture in its default view, whole, and focused on its first protein chain (a structure without
    protein: whole), with its generic, deposited and beta names, and with automatic domains: no chain carries one
    label on two elements, and no two stems one name (figcheck.label_repeats, the campaign's A8)."""
    if focus == "first":
        _, sses, _ = make_layout(path, focus="none")
        focus = sses[0].chain if sses else "none"
    named = None if names == "none" else [names]
    lay, _, _ = make_layout(path, focus=focus, names=named, domains=domains)
    assert _repeats(lay) == {}


def _strands(lay, sses, chain):
    return [lay.placed[s.id].label for s in sses if s.chain == chain and s.kind == "E" and s.id in lay.placed]


def test_pilus_subunit_reads_a_to_i_and_its_neighbours_lend_their_mates_labels():
    """6Y7S (focus D): FimA strands A-I once each (C and I were drawn twice: chain A, modelled with fewer strands, was
    the reference). Each lent element carries the label it has in the whole rod, plus its chain's prime."""
    lay, sses, _ = make_layout(DATA / "6Y7S.cif")
    assert lay.focus_chains == {"D"}
    assert _strands(lay, sses, "D") == list("ABCDEFGHI")
    whole, _, _ = make_layout(DATA / "6Y7S.cif", focus="none")
    lent = [p for p in lay.placed.values() if p.sse.chain in lay.partial]
    primes = {p.sse.chain: p.label[-1] for p in lent}
    assert lent and sorted(primes.values()) == sorted(set(primes.values())) and set(primes.values()) <= set("′″‴")
    for p in lent:
        assert p.label == whole.placed[p.sse.id].label + primes[p.sse.chain]


def test_whole_pilus_gives_each_copy_one_label_per_strand():
    """6Y7S drawn whole: chain A (strands 22-23 … 146-158, no 3-17 nor 122-123) reads like its mates, and chain F's two
    strands no other copy has (77-79, 82-86) take labels no copy uses (J, K) instead of E and F."""
    lay, sses, _ = make_layout(DATA / "6Y7S.cif", focus="none")
    label = {(p.sse.chain, p.sse.first.seq): p.label for p in lay.placed.values()}
    assert _strands(lay, sses, "A") == list("BCDEFHI")
    assert _strands(lay, sses, "B") == list("ABCDEFGHI")
    assert (label[("F", 77)], label[("F", 82)]) == ("J", "K")
    assert label[("F", 99)] == label[("B", 99)] == "F"


def test_helicase_neighbours_keep_their_primes():
    """6ZYA (focus A1): A2, A3, A4 lend elements primed ′, ″, ‴ as before; A2's 554-555, which A1 does not have,
    no longer reads as A1's Q (520-523) but takes a label of its own."""
    lay, sses, _ = make_layout(DATA / "6ZYA.cif")
    label = {(p.sse.chain, p.sse.first.seq): p.label for p in lay.placed.values()}
    assert {p.sse.chain: p.label[-1] for p in lay.placed.values() if p.sse.chain in lay.partial} == {
        "A2": "′",
        "A3": "″",
        "A4": "‴",
    }
    assert [label[("A3", n)] for n in (418, 430, 438)] == ["I″", "J″", "K″"]
    assert [label[("A4", n)] for n in (445, 476, 563)] == ["L‴", "N‴", "T‴"]
    assert label[("A1", 520)] == "Q" and label[("A2", 554)] not in {f"{c}′" for c in "ABCDEFGHIJKLMNOPQRST"}


def test_a_name_given_to_one_copy_reaches_one_element_of_each_other():
    """7TIM: B:178-203 is one helix, A holds it as 178-196 + 198-203. Named in B only, the name reaches A's
    178-196 (the larger overlap) and not 198-203 too."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        names = Path(tmp) / "b.yaml"
        names.write_text("elements:\n  res:B:178-203: H8\n", encoding="utf-8")
        lay, _, _ = make_layout(FID / "7TIM.cif.gz", names=[str(names)])
    assert lay.placed["B:178-203"].label == lay.placed["A:178-196"].label == "H8"
    assert lay.placed["A:198-203"].label != "H8"
    assert _repeats(lay) == {}


def test_tim_copies_name_each_helix_once():
    """7TIM: B's 178-203 is A's 178-196 (α8); A's 198-203, which B does not have, keeps its α9: no chain has a
    second α8."""
    lay, sses, _ = make_layout(FID / "7TIM.cif.gz")
    labels = {c: [lay.placed[s.id].label for s in sses if s.chain == c] for c in "AB"}
    assert labels["A"].count("α8") == labels["B"].count("α8") == 1
    assert lay.placed["A:198-203"].label == "α9" and lay.placed["B:178-203"].label == "α8"
