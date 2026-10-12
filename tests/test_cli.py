from foldmap.cli import main


def test_summary_output(ubq, capsys):
    assert main(["summary", str(ubq)]) == 0
    out = capsys.readouterr().out
    assert "76 residues, 1 chain(s)" in out
    assert "helix  A:23-34" in out
    assert "sheet 1: 5 strand(s)" in out


def test_missing_file_exit_code(tmp_path, capsys):
    assert main(["summary", str(tmp_path / "nope.cif")]) == 1
    assert "foldmap:" in capsys.readouterr().err


def test_fragment_without_secondary_structure(ubq, tmp_path, capsys):
    import gemmi

    full = tmp_path / "full.pdb"
    gemmi.read_structure(str(ubq)).write_pdb(str(full))
    keep = [l for l in full.read_text().splitlines() if l.startswith("ATOM") and int(l[22:26]) <= 3]
    tri = tmp_path / "tri.pdb"
    tri.write_text("\n".join(keep) + "\nEND\n")
    assert main(["summary", str(tri)]) == 0
    out = capsys.readouterr().out
    assert "3 residues" in out and "helices: 0  strands: 0" in out and "sheet" not in out


def test_plot_writes_every_requested_file(ubq, tmp_path, capsys):
    svg, png = tmp_path / "out" / "u.svg", tmp_path / "out" / "u.png"
    assert main(["plot", str(ubq), "-o", str(svg), "-o", str(png), "--mode", "stack", "--rotate", "90"]) == 0
    assert svg.stat().st_size > 1000 and png.read_bytes()[:4] == b"\x89PNG"
    assert str(svg) in capsys.readouterr().out


def test_plot_missing_file_exit_code(tmp_path, capsys):
    assert main(["plot", str(tmp_path / "nope.cif"), "-o", str(tmp_path / "x.svg")]) == 1
    assert "foldmap:" in capsys.readouterr().err


def test_plot_bad_extension_exit_code(ubq, tmp_path, capsys):
    assert main(["plot", str(ubq), "-o", str(tmp_path / "x.tiff")]) == 1
    assert "extension" in capsys.readouterr().err


def test_plot_layout_uses_3d_contacts(zs5):
    from foldmap.cli import _bridge_contacts, _bridge_springs, make_layout
    from foldmap.features import sse_contacts
    from foldmap.frame import view_frame
    from foldmap.io import load_links
    from foldmap.layout import build_layout
    from helpers import pipeline

    bb, sses, sheets = pipeline(zs5)
    links = load_links(zs5, bb)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, links).items():  # disulfides count as contacts and as springs
        contacts[key] = contacts.get(key, 0) + n
    want = build_layout(
        sses, sheets, view_frame(sses), "projected", contacts=contacts, bridges=_bridge_springs(bb, sses, links)
    )
    got, _, _ = make_layout(zs5)
    assert {k: (p.cx, p.cy) for k, p in got.placed.items()} == {k: (p.cx, p.cy) for k, p in want.placed.items()}


def test_plot_accepts_palette_style_and_loop_options(ubq, tmp_path):
    out = tmp_path / "s.svg"
    args = ["plot", str(ubq), "-o", str(out), "--palette", "tol-muted", "--style", "pale", "--loops", "curved"]
    assert main(args) == 0 and out.stat().st_size > 1000
    import subprocess
    import sys

    bad = subprocess.run(
        [sys.executable, "-m", "foldmap", "plot", str(ubq), "-o", str(out), "--palette", "rainbow"],
        capture_output=True,
        text=True,
    )
    assert bad.returncode == 2 and "invalid choice" in bad.stderr


def test_plot_warns_about_loops_drawn_as_fallback_curves(ubq, tmp_path, capsys, monkeypatch):
    import foldmap.cli as cli
    import foldmap.spec as spec

    real = cli.route_loops

    def one_stuck(*a, **k):  # the test structures all route cleanly; mark one loop as a fallback
        loops = real(*a, **k)
        loops[1].fallback = True
        return loops

    monkeypatch.setattr(spec, "route_loops", one_stuck)
    assert main(["plot", str(ubq), "-o", str(tmp_path / "u.svg")]) == 0
    err = capsys.readouterr().err
    assert "A:12-16>A:23-34" in err and "no clear route" in err


def test_plot_is_quiet_when_every_loop_is_routed(ubq, tmp_path, capsys):
    assert main(["plot", str(ubq), "-o", str(tmp_path / "u.svg")]) == 0
    assert capsys.readouterr().err == ""


def test_missing_extras_are_explained(monkeypatch, capsys):
    import foldmap.cli as cli

    monkeypatch.setattr(cli, "_has", lambda name: False, raising=False)
    assert cli.main(["mcp"]) == 1
    assert "pip install 'foldmap[mcp]'" in capsys.readouterr().err
    assert cli.main(["app"]) == 1
    assert "pip install 'foldmap[app]'" in capsys.readouterr().err


def test_empty_title_option_clears_the_layout_files_title(ubq, tmp_path):
    lf = tmp_path / "fig.yaml"
    lf.write_text("foldmap: 1\nlayout: {title: Layout title}\n", encoding="utf-8")
    kept, cleared = tmp_path / "kept.svg", tmp_path / "cleared.svg"
    assert main(["plot", str(ubq), "-o", str(kept), "--layout-file", str(lf)]) == 0
    assert main(["plot", str(ubq), "-o", str(cleared), "--layout-file", str(lf), "--title", ""]) == 0
    assert "Layout title" in kept.read_text() and "Layout title" not in cleared.read_text()


def test_every_subcommand_has_a_handler():
    from foldmap.cli import _COMMANDS, build_parser

    sub = next(a for a in build_parser()._actions if a.dest == "command")
    assert set(sub.choices) == set(_COMMANDS)


def test_loops_drawn_over_text_are_listed_and_warned(tmp_path, capsys):
    """A loop routed over a label or chain-end stub as a last resort (Loop.crosses_labels) is listed like a loop with
    no clear route (Drawn.over_labels, beside Drawn.stuck) and warned about. 2ACY with helix_on_sheet=auto: the loop
    into α2, whose port lies on β1 by β1's letter."""
    from foldmap.spec import draw_spec, figure_spec
    from helpers import FID

    drawn = draw_spec(figure_spec(FID / "2ACY.cif.gz", style=["helix_on_sheet=auto"]))
    assert drawn.over_labels == ["A:47-54>A:55-67"] and not drawn.stuck
    out = tmp_path / "a.svg"
    assert main(["plot", str(FID / "2ACY.cif.gz"), "-o", str(out), "--set", "helix_on_sheet=auto"]) == 0
    err = capsys.readouterr().err
    assert "warning" in err and "over a label" in err and "A:47-54>A:55-67" in err
    assert main(["plot", str(FID / "2ACY.cif.gz"), "-o", str(out)]) == 0
    assert "over a label" not in capsys.readouterr().err


def test_plot_says_when_an_automatic_choice_drew_part_of_the_file(tmp_path, capsys):
    """One stderr note naming what was drawn and the flags that draw everything; none when everything is drawn."""
    from pathlib import Path

    apo = Path(__file__).parent / "data" / "assemblies" / "7A4M.cif"
    assert main(["plot", str(apo), "-o", str(tmp_path / "a.svg")]) == 0
    err = capsys.readouterr().err
    note = "foldmap: drew one of 24 subunits of an O cage: 1 of 24 chains; --focus none draws all of them"
    assert err.count("foldmap: drew ") == 1 and note in err
    assert main(["plot", str(apo), "-o", str(tmp_path / "b.svg"), "--focus", "none"]) == 0
    assert "foldmap: drew " not in capsys.readouterr().err
    kelch = Path(__file__).parent / "data" / "4CH9.cif"
    assert main(["plot", str(kelch), "-o", str(tmp_path / "k.svg")]) == 0
    assert (
        "foldmap: drew assembly 1, one of 2 copies in the file: 2 of 4 chains; --assembly asu --focus none draws all "
        "of them" in capsys.readouterr().err
    )
