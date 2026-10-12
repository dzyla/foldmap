"""--names: conventional element and stem names (the file's own helix ids, β1…βn numbering, YAML name files)."""

from pathlib import Path

import pytest

from foldmap.cli import make_layout
from helpers import FID

DATA = Path(__file__).parent / "data"
DOCS = Path(__file__).parent.parent / "docs" / "names"


def _labels(lay, sses, chain):
    return {s.first.seq: lay.placed[s.id].label for s in sses if s.chain == chain and s.id in lay.placed}


def test_myoglobin_helices_take_their_deposited_names():
    """1MBN struct_conf: A 3-18, B 20-35, C 36-42, D 51-57, E 58-77, F 86-94, G 100-118, H 125-148."""
    lay, sses, _ = make_layout(FID / "1MBN.cif.gz", names=["deposited"])
    label = {s.id: lay.placed[s.id].label for s in sses}
    assert [label[k] for k in ("A:4-17", "A:21-35", "A:59-76", "A:83-97", "A:101-118", "A:125-149")] == list("ABEFGH")
    assert label["A:44-47"] == "α4"  # the CD turn overlaps no deposited helix: it keeps its own label


def test_haemoglobin_alpha_a_and_b_are_named_after_the_kink_split():
    lay, sses, _ = make_layout(DATA / "2HHB.cif", names=["deposited"])
    labels = _labels(lay, sses, "A")
    assert labels[4] == "A" and labels[21] == "B"  # Task 2 split; deposited AA 3-18, AB 20-35


def test_beta_names_strands_in_sequence():
    lay, sses, _ = make_layout(DATA / "1UBQ.cif", names=["beta"])
    assert [lay.placed[s.id].label for s in sses if s.kind == "E"] == ["β1", "β2", "β3", "β4", "β5"]
    assert lay.placed["A:23-34"].label == "α1"


def test_beta_names_number_template_helices_after_their_strand():
    lay, sses, _ = make_layout(DATA / "1TIM.cif", chains=["A"], names=["beta"])
    row = [lay.placed[k] for k in next(b for b in lay.sheets if len(b) == 8)]
    names = {p.label for p in lay.placed.values() if p.sse.kind == "H"}
    assert {f"α{n}" for n in range(1, 9)} <= names and len(row) == 8


@pytest.mark.parametrize(
    "file, structure, kw",
    [
        ("2HHB-globin.yaml", DATA / "2HHB.cif", {}),
        ("1TEN-fn3.yaml", FID / "1TEN.cif.gz", {}),
        ("1KX5-histones.yaml", DATA / "na" / "1KX5.cif.gz", {"focus": "A,B"}),
        ("1EHZ-trna.yaml", DATA / "na" / "1EHZ.cif.gz", {}),
    ],
)
def test_shipped_name_files_apply_cleanly(file, structure, kw):
    lay, sses, _ = make_layout(structure, names=[str(DOCS / file)], **kw)
    assert lay.placed or lay.stems


def test_trna_arms_by_name():
    lay, _, _ = make_layout(DATA / "na" / "1EHZ.cif.gz", names=[str(DOCS / "1EHZ-trna.yaml")])
    assert sorted(d.name for d in lay.stems) == ["D", "T", "acceptor", "anticodon"]


def test_unknown_reference_in_a_names_file_is_an_error(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("elements:\n  res:A:900-910: Z\n", encoding="utf-8")
    with pytest.raises(ValueError, match="bad.yaml.*res:A:900-910"):
        make_layout(DATA / "1UBQ.cif", names=[str(bad)])


def test_names_survive_a_layout_file_round_trip(tmp_path):
    from foldmap import layoutfile

    doc = {"foldmap": 1, "names": ["beta"]}
    path = layoutfile.save(tmp_path / "f.yaml", doc)
    assert layoutfile.load(path)["names"] == ["beta"]


# ---------------------------------------------------------------------------------------------- beyond the brief


def test_deposited_ids_lose_the_chain_prefix_and_generic_ids_are_ignored():
    """2HHB AA/AB -> A/B and CA/CB -> A/B; 7TIM AD1 -> D1; numeric, H<n> and the wwPDB serial ids (AA1, AA2 ... AB1,
    numbered across the whole file, as in 8UTF) name nothing; an id used twice in a chain names nothing."""
    from foldmap.names import helix_names

    rows = [("AA", "A"), ("AB", "A"), ("CA", "C"), ("CB", "C"), ("1", "D"), ("H2", "D"), ("X", "E"), ("X", "E")]
    assert helix_names(rows) == ["A", "B", "A", "B", None, None, None, None]
    assert helix_names([("AD1", "A"), ("AD2", "A"), ("AH1", "A")]) == ["D1", "D2", "H1"]
    assert helix_names([("A1F", "A"), ("AH", "A")]) == ["1F", "H"]
    assert helix_names([("A1", "A"), ("A2", "A")]) == [None, None]  # A + a number: generic once the chain is gone
    serial = [("AA1", "A"), ("AA2", "A"), ("AA3", "B")] + [(f"AA{k}", "a") for k in range(4, 10)] + [("AB1", "a")]
    assert helix_names(serial) == [None] * len(serial)
    assert helix_names([("A", "A"), ("B", "A")]) == ["A", "B"]  # 1MBN: nothing to strip (nothing would be left)


def test_copies_of_haemoglobin_take_the_same_deposited_names():
    lay, sses, _ = make_layout(DATA / "2HHB.cif", names=["deposited"])
    assert _labels(lay, sses, "C") == _labels(lay, sses, "A")
    beta = _labels(lay, sses, "B")
    assert [beta[k] for k in (5, 20, 36, 51, 58, 78, 101, 124)] == list("ABCDEFGH") and beta[43] == "η1"


def test_beta_names_extra_helices_with_a_prime_after_their_strand():
    """1TIM A: the row holds the longest helix of each βα unit. β4 is followed by 96-102 and 106-118 (in the row:
    α4), so 96-102 is α4′; β8's helix 239-244 closes the barrel outside the row (opened at β8|β1): α8, the name
    being free. 3-10 helices keep their η names."""
    lay, sses, _ = make_layout(DATA / "1TIM.cif", chains=["A"], names=["beta"])
    label = {s.id: lay.placed[s.id].label for s in sses}
    assert label["A:106-118"] == "α4" and label["A:96-102"] == "α4′"
    assert label["A:131-136"] == "α5′" and label["A:139-153"] == "α5"
    assert label["A:239-244"] == "α8" and label["A:167-169"] == "η1"
    assert [label[s.id] for s in sses if s.kind == "E"] == [f"β{n}" for n in range(1, 9)]


def test_later_names_win_and_an_explicit_rename_wins_over_all(tmp_path):
    names = tmp_path / "n.yaml"
    names.write_text("elements:\n  β2: second\n  res:A:23-34: helix\n", encoding="utf-8")
    lay, _, _ = make_layout(DATA / "1UBQ.cif", names=["beta", str(names)], rename={"res:A:23-34": "mine"})
    assert lay.placed["A:12-16"].label == "second"  # the file's label refs read the names given before it
    assert lay.placed["A:23-34"].label == "mine"


def test_mirrored_copies_take_their_mates_names(tmp_path):
    """1LMB: the two repressor protomers are drawn as mirror images (Task 8). Names given for one copy only are
    shared with the other, so both read alike."""
    names = tmp_path / "n.yaml"
    names.write_text("elements:\n  res:3:9-29: h1\n  res:3:44-51: h3\n", encoding="utf-8")
    lay, sses, _ = make_layout(DATA / "1LMB.cif", names=[str(names)])
    assert lay.frame.mirror_x is not None
    assert lay.placed["4:9-29"].label == "h1" and lay.placed["4:44-51"].label == "h3"
    assert _labels(lay, sses, "4") == _labels(lay, sses, "3")


def test_neighbouring_fragments_keep_their_prime_under_names():
    """1KX5 drawn as H3 + H4 (focus A,B): elements lent by neighbours are primed; a name keeps the prime."""
    lay, sses, _ = make_layout(DATA / "na" / "1KX5.cif.gz", focus="A,B", names=[str(DOCS / "1KX5-histones.yaml")])
    assert lay.placed["A:86-113"].label == "α2" and lay.placed["B:50-75"].label == "α2"
    lent = [lay.placed[s.id].label for s in sses if s.chain in lay.partial and s.kind == "H"]
    assert lent and all(name[-1] in "′″‴" for name in lent)
    assert any(name.startswith("α") and name[1] in "123NC" for name in lent)


def test_names_for_chains_left_out_of_the_figure_are_skipped():
    lay, sses, _ = make_layout(DATA / "2HHB.cif", chains=["A"], names=[str(DOCS / "2HHB-globin.yaml")])
    assert _labels(lay, sses, "A")[4] == "A"


def test_stem_names_are_measured_as_named():
    """Stem name sides are chosen once the names are final: choosing again changes nothing."""
    from foldmap.layout import _choose_label_sides
    from foldmap.style import Style

    lay, sses, _ = make_layout(DATA / "na" / "1EHZ.cif.gz", names=[str(DOCS / "1EHZ-trna.yaml")])
    slots = [(d.label_side, d.label_out) for d in lay.stems]
    _choose_label_sides(lay, sses, Style())
    assert [(d.label_side, d.label_out) for d in lay.stems] == slots


def test_names_reach_the_layout_file_and_the_command_line(tmp_path):
    from foldmap import layoutfile
    from foldmap.cli import _spec_of, build_parser
    from foldmap.spec import figure_spec, layout_document

    args = build_parser().parse_args(["plot", str(DATA / "1UBQ.cif"), "--names", "beta", "-o", "x.svg"])
    spec = _spec_of(args)
    assert spec.names == ["beta"] and spec.layout_kwargs(edits=False)["names"] == ["beta"]
    doc = layout_document(spec)
    assert doc["names"] == ["beta"]
    path = layoutfile.save(tmp_path / "f.yaml", doc)
    assert figure_spec(DATA / "1UBQ.cif", layout_file=path).names == ["beta"]
    both = figure_spec(DATA / "1UBQ.cif", layout_file=path, names=["deposited", "beta"])
    assert both.names == ["deposited", "beta"]  # the file's, then the ones given (a repeat keeps its last place)


def test_a_names_spec_that_is_neither_a_keyword_nor_a_file_is_an_error():
    with pytest.raises(ValueError, match="deposited, beta or a YAML file"):
        make_layout(DATA / "1UBQ.cif", names=["betta"])


def test_a_names_file_holds_elements_and_stems_only(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("helices:\n  α1: A\n", encoding="utf-8")
    with pytest.raises(ValueError, match="bad.yaml.*helices"):
        make_layout(DATA / "1UBQ.cif", names=[str(bad)])


def test_a_stem_key_must_name_a_paired_nucleotide(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("stems:\n  A:8: hinge\n", encoding="utf-8")  # U8 pairs with nothing
    with pytest.raises(ValueError, match="bad.yaml.*A:8"):
        make_layout(DATA / "na" / "1EHZ.cif.gz", names=[str(bad)])


def test_layout_text_measures_every_prime_a_name_may_hold():
    from foldmap import text_metrics

    assert all(c in text_metrics.ADVANCE for c in "αβη′″‴")
