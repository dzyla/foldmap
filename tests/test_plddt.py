"""AlphaFold confidence colouring: pLDDT (stored in the B-factor column) in AlphaFold's four fixed bands."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.colors import to_hex

from foldmap.cli import make_figure, make_layout
from foldmap.style import THEMES, Style

DATA = Path(__file__).parent / "data"
P53 = DATA / "AF-P04637.cif"
BANDS = {"very high": "#0053d6", "confident": "#65cbf3", "low": "#ffdb13", "very low": "#ff7d45"}


def _band(v):
    return "very high" if v >= 90 else "confident" if v >= 70 else "low" if v >= 50 else "very low"


def test_elements_take_the_band_of_their_mean_plddt():
    from foldmap.render import element_colours

    lay, sses, _ = make_layout(P53)
    colours, residue = element_colours(lay, sses, Style(color_by="plddt"))
    seen = set()
    for s in sses:
        band = _band(np.mean(lay.res_b[s.start : s.end + 1]))
        assert to_hex(colours[s.id]) == BANDS[band], s.id
        seen.add(band)
    assert len(seen) >= 2
    tail = int(np.argmin(lay.res_b))  # p53's disordered tails hold no elements: their residues carry the band
    assert lay.res_b[tail] < 50 and to_hex(residue(tail, lay.res_chain[tail])) == BANDS["very low"]


def test_bands_are_absolute_not_rescaled():
    from foldmap.render import element_colours

    lay, sses, _ = make_layout(DATA / "AF-P04637.cif")
    lay.res_b = [95.0] * len(lay.res_b)  # all very high: no stretching to fill the scale
    colours, residue = element_colours(lay, sses, Style(color_by="plddt"))
    assert {to_hex(c) for c in colours.values()} == {BANDS["very high"]}
    assert to_hex(residue(0, sses[0].chain)) == BANDS["very high"]


def test_legend_shows_the_four_bands():
    fig = make_figure(P53, look=Style(color_by="plddt"))
    texts = {t.get_text() for t in fig.axes[0].texts}
    assert {"very high (> 90)", "confident (70–90)", "low (50–70)", "very low (< 50)"} <= texts


def test_alphafold_theme():
    assert THEMES["alphafold"].color_by == "plddt"


def test_cli_warns_when_b_factors_do_not_look_like_plddt(tmp_path, capsys):
    from foldmap.cli import main

    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(tmp_path / "x.svg"), "--theme", "alphafold"]) == 0
    assert "pLDDT" in capsys.readouterr().err
    assert main(["plot", str(P53), "-o", str(tmp_path / "y.svg"), "--theme", "alphafold"]) == 0
    assert "pLDDT" not in capsys.readouterr().err


@pytest.mark.parametrize(
    "text, url",
    [
        ("AF-P04637-F1", "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.cif"),
        ("P04637", "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.cif"),
        ("af-q9y6k9", "https://alphafold.ebi.ac.uk/files/AF-Q9Y6K9-F1-model_v6.cif"),
    ],
)
def test_app_fetches_alphafold_models(text, url, monkeypatch, tmp_path):
    import foldmap.app as app

    asked = []

    class Reply:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self):
            return P53.read_bytes()

    def fake(request, timeout=0):
        asked.append(getattr(request, "full_url", request))
        return Reply()

    monkeypatch.setattr(app, "CACHE", tmp_path)
    monkeypatch.setattr(app.urllib.request, "urlopen", fake)
    path = app.fetch(text)
    assert asked == [url] and path.read_bytes() == P53.read_bytes()
    assert app.fetch(text) == path and len(asked) == 1  # cached


def _segs(fig):
    out = {}
    for a in fig.axes[0].patches:
        gid = a.get_gid() or ""
        if gid.startswith("loop-seg:"):
            name, k = gid[len("loop-seg:") :].rsplit(":", 1)
            out.setdefault(name, []).append((int(k), to_hex(a.get_edgecolor())))
    return {name: [c for _, c in sorted(v)] for name, v in out.items()}


def test_residue_loops_follow_the_residues_they_span():
    from foldmap.route import route_loops

    lay, sses, bb = make_layout(P53)
    fig = make_figure(P53, look=Style(color_by="plddt", loop_color="residue"))
    segs = _segs(fig)
    by_id = {s.id: s for s in sses}
    assert segs
    for loop in route_loops(lay, sses, bb):
        a, b = by_id[loop.a_id], by_id[loop.b_id]
        inside = lay.res_b[a.end + 1 : b.start] or [lay.res_b[a.end], lay.res_b[b.start]]
        colours = segs[f"{loop.a_id}>{loop.b_id}"]
        assert colours[0] == BANDS[_band(inside[0])] and colours[-1] == BANDS[_band(inside[-1])]
        assert set(colours) <= {BANDS[_band(v)] for v in inside}


def test_disordered_tails_show_on_the_termini():
    lay, sses, _ = make_layout(P53)
    fig = make_figure(P53, look=THEMES["alphafold"])
    stubs = {
        a.get_gid(): to_hex(a.get_edgecolor()) for a in fig.axes[0].patches if (a.get_gid() or "").startswith("stub:")
    }
    n_tail = lay.res_b[: min(s.start for s in sses)]
    assert stubs["stub:N:A"] == BANDS[_band(np.mean(n_tail))] == BANDS["very low"]


def test_residue_loops_work_for_every_colouring():
    for by in ("chain", "sequence", "bfactor", "hydropathy", "sstype"):
        fig = make_figure(DATA / "1UBQ.cif", look=Style(color_by=by, loop_color="residue"))
        assert _segs(fig), by
