"""Layout files: save everything that defines a figure, reload it, get the same figure."""

from pathlib import Path

import pytest
import yaml

from foldmap import layoutfile
from foldmap.cli import main

DATA = Path(__file__).parent / "data"


def _svg(tmp_path, name, *args):
    out = tmp_path / name
    assert main(["plot", *args, "-o", str(out)]) == 0
    return out.read_bytes()


def test_round_trip_gives_the_same_figure(tmp_path):
    lf = tmp_path / "fig.yaml"
    args = [
        str(DATA / "2LZM.cif"),
        "--theme",
        "trace",
        "--set",
        "residue_numbers=true",
        "--swap",
        "α1,α3",
        "--move",
        "α2=1,-2",
        "--rename",
        "α5=HelixE",
        "--set",
        "mark=α4=#2ca02c",
    ]
    first = _svg(tmp_path, "a.svg", *args, "--save-layout", str(lf))
    again = _svg(tmp_path, "b.svg", str(DATA / "2LZM.cif"), "--layout-file", str(lf))
    assert first == again


def test_saved_references_are_residue_ranges(tmp_path):
    lf = tmp_path / "fig.yaml"
    _svg(
        tmp_path,
        "a.svg",
        str(DATA / "2LZM.cif"),
        "--swap",
        "α1,α3",
        "--move",
        "α2=1,-2",
        "--rename",
        "α5=E",
        "--set",
        "mark=α4",
        "--save-layout",
        str(lf),
    )
    doc = yaml.safe_load(lf.read_text())
    refs = [*doc["edits"]["swap"][0], *doc["edits"]["move"], *doc["edits"]["rename"]]
    assert all(r.startswith("res:A:") for r in refs)
    assert doc["style"]["mark"].startswith("res:A:")
    assert doc["foldmap"] == 1 and doc["theme"] == "publication"


def test_command_line_wins_over_the_file(tmp_path):
    lf = tmp_path / "fig.yaml"
    _svg(tmp_path, "a.svg", str(DATA / "1UBQ.cif"), "--theme", "trace", "--save-layout", str(lf))
    resaved = tmp_path / "again.yaml"
    _svg(
        tmp_path,
        "b.svg",
        str(DATA / "1UBQ.cif"),
        "--layout-file",
        str(lf),
        "--theme",
        "print",
        "--mode",
        "stack",
        "--save-layout",
        str(resaved),
    )
    doc = yaml.safe_load(resaved.read_text())
    assert doc["theme"] == "print" and doc["layout"]["mode"] == "stack"


def test_edits_from_file_and_command_line_add_up(tmp_path):
    lf = tmp_path / "fig.yaml"
    _svg(tmp_path, "a.svg", str(DATA / "2LZM.cif"), "--move", "α2=1,0", "--save-layout", str(lf))
    both = tmp_path / "both.yaml"
    _svg(
        tmp_path,
        "b.svg",
        str(DATA / "2LZM.cif"),
        "--layout-file",
        str(lf),
        "--move",
        "α3=0,1",
        "--save-layout",
        str(both),
    )
    assert len(yaml.safe_load(both.read_text())["edits"]["move"]) == 2


@pytest.mark.parametrize(
    "text, msg",
    [
        ("foldmap: 1\ncolour: red\n", "colour"),
        ("foldmap: 9\n", "version"),
        ("- a\n- b\n", "mapping"),
        ("names: " + "[" * 3000 + "]" * 3000 + "\n", "nested too deeply"),  # not a RecursionError traceback
    ],
    ids=["unknown-key", "version", "list", "deep"],  # short ids (PYTEST_CURRENT_TEST has a size limit on Windows)
)
def test_bad_layout_files_are_explained(tmp_path, text, msg, capsys):
    lf = tmp_path / "bad.yaml"
    lf.write_text(text)
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(tmp_path / "x.svg"), "--layout-file", str(lf)]) == 1
    assert msg in capsys.readouterr().err


def test_layout_file_drives_the_interactive_page_too(tmp_path):
    lf = tmp_path / "fig.yaml"
    _svg(tmp_path, "a.svg", str(DATA / "1UBQ.cif"), "--theme", "rainbow", "--save-layout", str(lf))
    out = tmp_path / "v.html"
    assert main(["interactive", str(DATA / "1UBQ.cif"), "-o", str(out), "--layout-file", str(lf)]) == 0


def test_app_settings_round_trip_through_a_layout_file(tmp_path):
    from foldmap.app import layout_document, settings_from_document

    settings = {
        "theme": "trace",
        "residue_numbers": True,
        "mode": "stack",
        "swap": "α1,α3",
        "move": "α2=1,-2",
        "rename": "",
        "rotate": 30,
        "symmetry": "auto",
        "assembly": "auto",
        "title": "T4L",
    }
    doc = layout_document(DATA / "2LZM.cif", settings)
    back = settings_from_document(doc)
    assert back["theme"] == "trace" and back["residue_numbers"] is True and back["mode"] == "stack"
    assert back["rotate"] == 30 and back["title"] == "T4L" and back["swap"].startswith("res:A:")


def test_version_flag(capsys):
    from foldmap import __version__

    with pytest.raises(SystemExit):
        main(["--version"])
    assert capsys.readouterr().out.strip() == f"foldmap {__version__}"


def test_interactive_page_credits_foldmap():
    from foldmap.interactive import build_page

    page = build_page(DATA / "1UBQ.cif")
    assert 'id="credit"' in page and "foldmap" in page and "AGPL-3.0" in page
    assert 'see <a href="https://github.com/dzyla/foldmap/blob/main/NOTICE">NOTICE</a>' in page


def test_app_credits_foldmap():
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(str(Path(__file__).parents[1] / "src" / "foldmap" / "app.py")).run(timeout=60)
    assert any("AGPL-3.0" in c.value and "cite" in c.value for c in at.caption)


def test_invalid_yaml_is_a_readable_error(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("layout: {rotate: 30\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"bad\.yaml: not valid YAML"):
        layoutfile.load(bad)


def test_na_extent_round_trips_and_dna_extent_is_read_as_its_old_name(tmp_path):
    from foldmap.app import layout_document, settings_from_document
    from foldmap.spec import figure_spec

    doc = layout_document(DATA / "1UBQ.cif", {"theme": "publication", "na_extent": "contacts"})
    assert doc["style"]["na_extent"] == "contacts"
    assert settings_from_document(doc)["na_extent"] == "contacts"
    lf = tmp_path / "old.yaml"
    lf.write_text("foldmap: 1\nstyle:\n  dna_extent: all\n")
    assert figure_spec(DATA / "1UBQ.cif", layout_file=lf).look.na_extent == "all"
    assert settings_from_document(layoutfile.load(lf))["na_extent"] == "all"


def test_two_fold_round_trips(tmp_path):
    """Plan 18 Task 2: the two-fold choice is saved with the figure and read back (app and command line alike)."""
    from foldmap.app import layout_document, settings_from_document
    from foldmap.spec import figure_spec

    doc = layout_document(DATA / "1UBQ.cif", {"theme": "publication", "two_fold": "unrolled"})
    assert doc["style"]["two_fold"] == "unrolled"
    assert settings_from_document(doc)["two_fold"] == "unrolled"
    lf = tmp_path / "fig.yaml"
    lf.write_text(yaml.safe_dump(doc), encoding="utf-8")
    assert figure_spec(DATA / "1UBQ.cif", layout_file=lf).look.two_fold == "unrolled"
    assert "two_fold" not in layout_document(DATA / "1UBQ.cif", {"theme": "publication"})["style"]  # the default


def test_copy_ghosts_round_trips(tmp_path):
    """Plan 19 Task 2: the copy-ghosts choice is saved with the figure and read back (app and command line alike)."""
    from foldmap.app import layout_document, settings_from_document
    from foldmap.spec import figure_spec

    doc = layout_document(DATA / "1UBQ.cif", {"theme": "publication", "copy_ghosts": "off"})
    assert doc["style"]["copy_ghosts"] == "off"
    assert settings_from_document(doc)["copy_ghosts"] == "off"
    lf = tmp_path / "fig.yaml"
    lf.write_text(yaml.safe_dump(doc), encoding="utf-8")
    assert figure_spec(DATA / "1UBQ.cif", layout_file=lf).look.copy_ghosts == "off"
    assert "copy_ghosts" not in layout_document(DATA / "1UBQ.cif", {"theme": "publication"})["style"]  # the default


def test_an_old_layout_file_with_two_fold_mirror_still_mirrors(tmp_path):
    """Plan 20 Task 1: two_fold gained auto as its default; a layout file written with mirror (plan 18's default,
    written out by name) still draws haemoglobin as one copy and its mirror image, which auto no longer does."""
    from foldmap.cli import make_layout
    from foldmap.frame import Frame
    from foldmap.spec import figure_spec

    lf = tmp_path / "fig.yaml"
    lf.write_text("foldmap: 1\nstyle:\n  two_fold: mirror\n")
    look = figure_spec(DATA / "2HHB.cif", layout_file=lf).look
    assert look.two_fold == "mirror"
    lay, _, _ = make_layout(DATA / "2HHB.cif", look=look)
    assert type(lay.frame) is Frame and lay.frame.mirror_x is not None


def test_a_figure_drawn_like_a_reference_round_trips(tmp_path):
    """A13 with --like: the layout file records the reference as given, its chain and the match, and the figure
    read back from it is the same (the SVG, and bench's signature of the layout)."""
    import importlib.util

    from foldmap.cli import make_layout
    from foldmap.spec import figure_spec

    script = Path(__file__).parent.parent / "scripts" / "bench.py"
    bench = importlib.util.spec_from_file_location("foldmap_bench", script)
    mod = importlib.util.module_from_spec(bench)
    bench.loader.exec_module(mod)
    path, lf = DATA / "2HHB.cif", tmp_path / "beta.yaml"
    args = [str(path), "--chains", "B", "--like", str(path), "--like-chain", "A", "--set", "color_by=sequence"]
    first = _svg(tmp_path, "a.svg", *args, "--save-layout", str(lf))
    assert yaml.safe_load(lf.read_text(encoding="utf-8"))["layout"]["like"] == {
        "structure": str(path),
        "chain": "A",
        "match": "auto",
    }
    again = _svg(tmp_path, "b.svg", str(path), "--layout-file", str(lf))
    assert first == again
    spec = figure_spec(path, layout_file=lf)
    assert spec.like == {"structure": str(path), "chain": "A", "match": "auto", "layout": None}
    read_back, _, _ = make_layout(path, spec.mode, spec.rotate, spec.flip_v, **spec.layout_kwargs())
    direct, _, _ = make_layout(
        path, chains=["B"], look=spec.look, like={"structure": str(path), "chain": "A", "match": "auto", "layout": None}
    )
    assert read_back.like_note.startswith("drawn like 2HHB.cif chain A: ")
    assert mod.signature(read_back) == mod.signature(direct)


@pytest.mark.parametrize(
    "like, msg",
    [("2HHB", "like is a mapping"), ({"chain": "A"}, "structure"), ({"structure": "2HHB", "colour": 1}, "colour")],
    ids=["text", "no-structure", "unknown"],
)
def test_bad_like_in_a_layout_file_is_explained(tmp_path, like, msg, capsys):
    lf = tmp_path / "bad.yaml"
    lf.write_text(yaml.safe_dump({"foldmap": 1, "layout": {"like": like}}))
    assert main(["plot", str(DATA / "1UBQ.cif"), "-o", str(tmp_path / "x.svg"), "--layout-file", str(lf)]) == 1
    assert msg in capsys.readouterr().err
