"""Plan 22: a helix's name keeps its gap to its bar at any Style.font_scale. Its centre stands LABEL_GAP times the
text's glyph_scale beyond the bar (helix_label_pos; further beside an upright bar, so its keep-out box reaches no
further into the bar than at font_scale 1), and its keep-out box with it; at a fixed LABEL_GAP, names drawn
at font_scale 2 or more sat on their own bars (Fig. 2: 2HHB --chains A,B and 1C3W one subunit at 2.3-2.8), and A6
never read it: a text over its own body was exempt. A helix's or 3-10 helix's bar under its own name now counts
(figcheck.BESIDE_BODIES); a strand's letter in its arrow, a residue number beside its end, stay the element's own."""

import math
from pathlib import Path

import pytest

from foldmap import criteria, figcheck
from foldmap.layout import _GLYPH_K, _GLYPH_W, _NAME_PAD, HELIX_W, LABEL_GAP, Placed, helix_label_pos
from foldmap.spec import draw_spec, figure_spec
from helpers import fake_sse
from test_figcheck import _body, _page, _text

DATA = Path(__file__).parent / "data"


def test_a_helix_name_stands_out_in_proportion_to_its_text():
    p = Placed(fake_sse("H", 1, 10), 0.0, 0.0, 4.0, 0.0, HELIX_W, "α1")
    assert _GLYPH_K.get() == 1.0 and helix_label_pos(p) == helix_label_pos(p, 1.0) == (0.0, HELIX_W / 2 + LABEL_GAP)
    assert helix_label_pos(p, 2.5) == pytest.approx((0.0, HELIX_W / 2 + 2.5 * LABEL_GAP))
    p.label_side = -1
    assert helix_label_pos(p, 2.5) == pytest.approx((0.0, -(HELIX_W / 2 + 2.5 * LABEL_GAP)))
    # an upright bar: the name's width lies across its side; its box (three glyphs) reaches no further into the bar
    # at 2.5 than at 1
    p.angle, p.label_side = math.pi / 2, 1
    reach = _GLYPH_W * 1.5 + _NAME_PAD
    for k in (1.0, 2.5):
        x, y = helix_label_pos(p, k)
        assert y == pytest.approx(0.0) and -x - reach * k == pytest.approx(HELIX_W / 2 + LABEL_GAP - reach)


def test_a_helix_name_on_its_own_bar_is_a_text_collision():
    """A helix's or 3-10 helix's name stands beside its bar: on it, A6 reads it, once per element; a strand's letter
    in its arrow and a residue number on its own end are the element's own."""
    for bar in ("helix:A:1-5", "helix-back:A:1-5", "eta:A:1-5"):
        d = _page(_text("label:A:1-5"), _body(bar, 2, 4, 8, 6))
        assert figcheck.text_collisions(d) == ["label:A:1-5 × " + bar.replace("helix-back:", "helix:")]
    d = _page(_text("label:A:1-5"), _body("helix-back:A:1-5", 2, 4, 8, 6), _body("helix:A:1-5", 2, 4, 8, 6, z=3.1))
    assert figcheck.text_collisions(d) == ["label:A:1-5 × helix:A:1-5"]  # its two faces: one element
    assert figcheck.text_collisions(_page(_text("label:A:1-5"), _body("strand:A:1-5", 2, 4, 8, 6))) == []
    assert figcheck.text_collisions(_page(_text("resnum:A:1-5:N"), _body("helix:A:1-5", 2, 4, 8, 6))) == []
    assert figcheck.text_collisions(_page(_text("label:A:1-5"), _body("helix:A:1-5", 2, 7, 8, 9))) == []  # beside


CASES = {
    "1UBQ": ("1UBQ.cif", {}),
    "2HHB A,B": ("2HHB.cif", {"chains": ["A", "B"]}),
    "1C3W one subunit": ("1C3W.cif", {"membrane": "auto", "focus": "protomer", "style": ["ligands=none"]}),
}


@pytest.mark.parametrize("font_scale", [1.0, 2.0, 3.0])
@pytest.mark.parametrize("case", list(CASES))
def test_names_keep_off_every_body_at_any_text_size(case, font_scale):
    """No helix name's ink box lies on any element's body (its own bar included; a strand's letter sits in its own
    arrow), and A4 and A6 pass, at text 1, 2 and 3 times its size (Fig. 2's structures; at a fixed LABEL_GAP 2HHB A,B
    and 1C3W had names on their bars from font_scale 2)."""
    name, options = CASES[case]
    options = dict(options)
    style = options.pop("style", []) + [f"font_scale={font_scale}"]
    spec = figure_spec(DATA / name, style=style, **options)
    drawn = draw_spec(spec)
    page = figcheck.measure(drawn)
    bodies = [(g, path) for g, path, _, _ in page.bodies if g.startswith(figcheck.ELEMENT_BODIES)]
    names = [(t.get_gid(), box) for t, _, box in page.texts if (t.get_gid() or "").startswith("label:")]
    assert names
    on = [
        (gid, g)
        for gid, box in names
        for g, path in bodies
        if not (figcheck._owner(g) == figcheck._owner(gid) and g.startswith(("strand:", "ghost:")))
        and figcheck._overlap(box, path.get_extents()) > 0
        and figcheck._area_on(box, path) > 0
    ]
    assert on == []
    got = criteria.check_figure(drawn, spec, DATA / name)["criteria"]
    if reason := figcheck.ink_font_reason(drawn.figure):  # A6's ink is read in Liberation Sans only: not measured here
        assert got["A4"]["pass"] and got["A6"]["pass"] is None and got["A6"]["here"] == reason, (got["A4"], got["A6"])
    else:
        assert got["A4"]["pass"] and got["A6"]["pass"], (got["A4"], got["A6"])
