"""The cross-platform determinism run (scripts/paper/determinism.py): its comparison fails on any difference, a signed
zero in a layout signature or one byte of an SVG included, and a run has the shape CI uploads and the paper reads."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "paper" / "determinism.py"
DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def det():
    spec = importlib.util.spec_from_file_location("paper_determinism", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(det, label, system="Linux", angle=0.0, svg=b"<svg>figure</svg>", drawn=True):
    sig = det.signature_text({"H1": [12.5, -3.0, angle], "_bounds": [0.0, 0.0, 100.0, 80.0]})
    rows = [
        det.figure_row("1UBQ.cif", 0, sig, svg),
        det.figure_row("2LZM.cif", 0, sig, b"<svg>other</svg>") if drawn else det.figure_row("2LZM.cif", 1, None, None),
    ]
    return {
        "label": label,
        "system": system,
        "python": "3.12.13",
        "numpy": "2.5.4",
        "matplotlib": "3.11.2",
        "gemmi": "0.7.5",
        "font": {"font_family": "Liberation Sans"},
        "figures": rows,
    }


def test_identical_runs_pass(det):
    runs = [_run(det, "ubuntu-latest"), _run(det, "macos-latest", "Darwin"), _run(det, "windows-latest", "Windows")]
    assert det.compare(runs) == []
    assert det.compare(runs, systems=3) == []


def test_a_signed_zero_in_the_signature_fails(det):
    """0.0 == -0.0, so parsed floats would agree: the signature is compared as text (and hash)."""
    runs = [_run(det, "seed0"), _run(det, "seed0-noavx", angle=-0.0)]
    assert runs[0]["figures"][0]["svg_sha256"] == runs[1]["figures"][0]["svg_sha256"]
    problems = det.compare(runs)
    assert any(p.startswith("1UBQ.cif: signature differs: [seed0] vs [seed0-noavx]") for p in problems), problems
    assert any("-0.0" in p for p in problems), problems  # the first difference is shown
    assert any(p.startswith("1UBQ.cif: signature_sha256 differs") for p in problems), problems
    assert not any(p.startswith("2LZM.cif: svg") for p in problems)


def test_one_svg_byte_fails(det):
    runs = [_run(det, "ubuntu-latest"), _run(det, "windows-latest", "Windows", svg=b"<svg>figurE</svg>")]
    assert det.compare(runs) == ["1UBQ.cif: svg_sha256 differs: [ubuntu-latest] vs [windows-latest]"]


def test_missing_figures_platforms_and_version_mismatches_fail(det):
    assert det.compare([_run(det, "a"), _run(det, "b", drawn=False)]) == ["2LZM.cif: not drawn in b"]
    assert det.compare([_run(det, "a"), _run(det, "b")], systems=3) == ["1 platform(s) (Linux), expected 3"]
    other = _run(det, "b") | {"matplotlib": "3.10.9"}
    assert det.compare([_run(det, "a"), other]) == ["matplotlib differs: a 3.11.2, b 3.10.9"]
    assert det.compare([_run(det, "a")]) == ["1 run(s): nothing to compare"]


def test_compare_exits_1_on_a_difference(det, tmp_path, capsys):
    same = [tmp_path / "a.json", tmp_path / "b.json"]
    for p, sys_ in zip(same, ("Linux", "Windows")):
        p.write_text(json.dumps(_run(det, p.stem, sys_)), encoding="utf-8")
    assert det.main(["--compare", "--systems", "2", *map(str, same)]) == 0
    assert "2 runs on Linux, Windows: 2 figures, identical" in capsys.readouterr().out
    (tmp_path / "c.json").write_text(json.dumps(_run(det, "c", angle=-0.0)), encoding="utf-8")
    assert det.main(["--compare", str(same[0]), str(tmp_path / "c.json")]) == 1
    assert "figures, differ" in capsys.readouterr().out


def test_a_run_has_the_shape_ci_uploads(tmp_path):
    """One small figure drawn by the script (any font here): what tests.yml uploads and the paper's harness reads."""
    out = tmp_path / "run.json"
    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(out), "--label", "here", "--any-font", "--svg-dir", str(tmp_path)]
        + [str(DATA / "1UBQ.cif")],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 0, done.stderr[-2000:]
    run = json.loads(out.read_text(encoding="utf-8"))
    assert {"label", "system", "machine", "python", "numpy", "matplotlib", "gemmi", "env", "font"} <= set(run)
    assert run["label"] == "here" and {"font_file", "font_family"} <= set(run["font"])
    (fig,) = run["figures"]
    assert fig["structure"] == "1UBQ.cif" and fig["exit"] == 0
    assert {"signature", "signature_sha256", "svg_sha256"} <= set(fig)
    assert (tmp_path / "1UBQ.svg").is_file()
