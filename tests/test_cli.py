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
    from helpers import pipeline

    from foldmap.cli import make_layout
    from foldmap.features import sse_contacts
    from foldmap.frame import view_frame
    from foldmap.layout import build_layout

    from foldmap.cli import _bridge_contacts, _bridge_springs
    from foldmap.io import load_links

    bb, sses, sheets = pipeline(zs5)
    links = load_links(zs5, bb)
    contacts = sse_contacts(bb, sses)
    for key, n in _bridge_contacts(bb, sses, links).items():  # disulfides count as contacts and as springs
        contacts[key] = contacts.get(key, 0) + n
    want = build_layout(sses, sheets, view_frame(sses), "projected", contacts=contacts,
                        bridges=_bridge_springs(bb, sses, links))
    got, _, _ = make_layout(zs5)
    assert {k: (p.cx, p.cy) for k, p in got.placed.items()} == {k: (p.cx, p.cy) for k, p in want.placed.items()}


def test_plot_accepts_palette_style_and_loop_options(ubq, tmp_path):
    out = tmp_path / "s.svg"
    args = ["plot", str(ubq), "-o", str(out), "--palette", "tol-muted", "--style", "pale", "--loops", "curved"]
    assert main(args) == 0 and out.stat().st_size > 1000
    import subprocess
    import sys

    bad = subprocess.run([sys.executable, "-m", "foldmap", "plot", str(ubq), "-o", str(out), "--palette", "rainbow"],
                         capture_output=True, text=True)
    assert bad.returncode == 2 and "invalid choice" in bad.stderr


def test_plot_warns_about_loops_drawn_as_fallback_curves(ubq, tmp_path, capsys, monkeypatch):
    import foldmap.cli as cli

    real = cli.route_loops

    def one_stuck(*a, **k):  # the test structures all route cleanly; mark one loop as a fallback
        loops = real(*a, **k)
        loops[1].fallback = True
        return loops

    monkeypatch.setattr(cli, "route_loops", one_stuck)
    assert main(["plot", str(ubq), "-o", str(tmp_path / "u.svg")]) == 0
    err = capsys.readouterr().err
    assert "A:12-16>A:23-34" in err and "no clear route" in err


def test_plot_is_quiet_when_every_loop_is_routed(ubq, tmp_path, capsys):
    assert main(["plot", str(ubq), "-o", str(tmp_path / "u.svg")]) == 0
    assert capsys.readouterr().err == ""
