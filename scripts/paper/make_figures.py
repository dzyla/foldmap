# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Manuscript figures (plan 14, M4): draw every foldmap panel with the installed CLI, then assemble each figure.

    python scripts/paper/make_figures.py TARGET ...

    draw [PANEL ...]          foldmap panels -> out/panels/<id>.svg/.pdf/.png + out/panels/panels.json
    f1 | f2 | f3 | f4 | f5    one figure -> out/F<n>.pdf/.svg/.png (300 dpi); several may follow
    fs                        the supplementary gallery -> out/FSgallery.pdf/.svg/.png
    gallery                   draw Fig. 2's and the supplementary gallery's panels, assemble both
    all                       draw, then every figure, then the rule checks (out/checks.json)
    check                     the rule checks only, on what is in out/

Run through figures.sh (here, or on a runner). It calls `foldmap` (the installed CLI of the paper environment) and
never edits repository code. Outputs go to build/figures/ of the manuscript folder (common.MANUSCRIPT; the commands
above write out/ there). Assembly: tectonic (XeLaTeX; $TECTONIC or on PATH) places each panel's PDF at its
final printed size and writes the panel letters; matplotlib draws the schematic (F1) and the charts (F5). SVG and
PNG come from the assembled PDF (pdftocairo, pdftoppm). Every size here is the printed size: single column 86 mm,
double column 178 mm (Bioinformatics).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent  # scripts/paper
sys.path.insert(0, str(HERE))
from common import BUILD, REPO, na_rows  # noqa: E402
from common import MANUSCRIPT as MS  # noqa: E402

OUT = BUILD / "figures"  # every output (gitignored, with the rest of the manuscript folder)
PANELS = OUT / "panels"
DATA = MS / "figures" / "data"  # inputs: the CI determinism artifacts F3d reads
TECTONIC = Path(os.environ.get("TECTONIC") or shutil.which("tectonic") or "tectonic")

MM = 72 / 25.4  # pt per mm
SINGLE_MM, DOUBLE_MM = 86.0, 178.0
MIN_PT = 6.0  # the figure rules' text floor, at the printed size
LETTER_PT = 9.0  # panel letters
FONT = "Liberation Sans"  # foldmap's own (render._FONTS; the only one of them on zyla-lab)

OKABE_ITO = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]

# ---------------------------------------------------------------------------------------------------------- panels
# Each figure is a list of boxes, in mm on the printed page: (letter, panel id, x, y, w, h), y from the top. A box
# keeps LETTER_ROOM mm at its top for the letter; the panel is scaled to fit the rest and centred in it.
LETTER_ROOM = 4.0
LAYOUTS: dict[str, list[tuple]] = {
    # Fig. 2, four large labelled panels (the author, 2026-10-21: "for gallery large is fine")
    # a, b on top, c, d below: every panel at font_scale <= 2.8, where `foldmap check` passes A4-A7 (NOTES.md,
    # 2026-10-21: gallery, second pass)
    "F2": [
        ("a", "F2a", 0, 0, 90, 74),
        ("b", "F2b", 92, 0, 86, 74),
        ("c", "F2c", 0, 74, 64, 104),
        ("d", "F2d", 66, 74, 112, 104),
    ],
    # Supplementary gallery: the other four structures, each panel half or full page width
    "FSgallery": [
        ("a", "S1a", 0, 0, 90, 86),
        ("b", "S1b", 92, 0, 86, 86),
        ("c", "S1c", 0, 86, 178, 56),
        ("d", "S1d", 0, 142, 100, 54),
        ("e", "S1e", 102, 142, 76, 54),
    ],
    "F3": [
        ("a", "F3a", 0, 0, 96, 105),
        ("b", "F3b", 98, 0, 80, 105),
        ("c", "F3c", 0, 105, 88, 66),
        ("d", "F3d", 90, 105, 88, 42),
    ],
    "F4": [
        ("a", "F4a", 0, 0, 82, 92),
        ("b", "F4b", 84, 0, 94, 142),
        ("c", "F4c", 0, 92, 82, 50),
        ("d", "F4d", 0, 142, 178, 92),
    ],
    "F1": [("", "F1x", 0.5, 30.0, 56.6, 34.0)],
}
# Panel ids follow the letters: F2a..F2d, S1a..S1e (NOTES.md, 2026-10-21: gallery)
# A panel drawn only to judge a choice, fitted in another panel's box (not placed in any figure): none now (1KX5 was
# tried in F2d's box, NOTES.md 2026-10-21: gallery)
BOX_AS: dict[str, str] = {}

# id -> foldmap arguments after `foldmap plot` (structure first). Defaults unless the caption asks for an option.
PANEL_ARGS: dict[str, list[str]] = {
    "F1x": ["1UBQ"],  # the small example in the schematic
    # Fig. 2, gallery: four large panels, labels on
    "F2a": ["1UBQ"],
    "F2b": ["1C3W", "--membrane", "auto", "--focus", "protomer", "--set", "ligands=none"],  # one subunit, retinal
    # not drawn: its name (2/3 of the label size) cannot print at 6 pt with the labels clear of the helices
    "F2c": ["1LMB"],
    "F2d": ["2HHB", "--chains", "A,B"],  # one alpha-beta protomer: the tetramer cannot print 6 pt at font_scale <= 2
    # Supplementary gallery
    "S1a": ["9FWZ"],
    "S1b": ["1EMA"],
    "S1c": ["1TIM"],
    "S1d": ["2HHB"],  # the tetramer, all four chains (author, 2026-10-22); unlabelled: its labels need a 134 x 104 mm
    # box at font_scale 3.5 (3.95 pt in 178 x 65 mm), which with the other panels does not fit the page
    "S1e": ["1KX5"],
    # Fig. 3, method details
    "F3a": ["9FWZ", "--chains", "D"],
    "F3b": ["1EMA"],
    "F3c": ["2LZM"],
    "F3d": ["tests/data/1TIM.cif"],  # the file the determinism run draws (scripts/paper/determinism.py)
    # Fig. 4, nucleic acids
    "F4a": ["1EHZ", "--names", "docs/names/1EHZ-trna.yaml"],
    "F4b": ["2GIS"],
    "F4c": ["4OO8", "--domains", "auto"],
    "F4d": ["1LMB", "--set", "hbonds=na"],
}
# The text floor. A panel is first drawn as is; if its smallest text prints under MIN_PT in its box, it is drawn again
# with a larger font_scale (foldmap's own style key: the layout makes room for the larger text), up to FS_MAX. Panels
# whose labels are the point (KEEP_LABELS) stop there; the others, and any that still fail, are drawn bare: no element,
# nucleotide or residue-number labels and no legend (BARE), the remaining text fitted the same way; a panel that
# still fails is drawn SILENT: the bare figure with every text hidden (foldmap_silent.py), so none prints under 6 pt.
FS_MAX = 3.0
KEEP_LABELS = {"F3c", "F4a", "F4b", "F4d"}
BARE = [
    "--set",
    "labels=false",
    "--set",
    "nucleotide_labels=false",
    "--set",
    "legend=false",
    "--set",
    "residue_numbers=false",
]
# F3d's picture is fitted like the others; the hashes are of the default figure, as the determinism run draws it
# (scripts/paper/determinism.py: Liberation Sans forced, PYTHONHASHSEED=0), run here on the remote as a fourth system
HASH_ONLY: dict = {}
FS_MAX_KEEP = 3.5  # panels whose labels are the point may go a little further
FS_MAX_LABELS = 2.5  # other panels, labels on: past this, labels outgrow narrow strands (9FWZ chain D at 3)
# Fig. 2's four panels keep their labels, font_scale as high as it takes (no bare or silent fallback); a panel that
# cannot reach the floor is kept at its best try and flagged under_floor
GALLERY_LABELS = {"F2a", "F2b", "F2c", "F2d"}
FS_MAX_GALLERY = 2.8  # past this the labels collide with elements (foldmap check A4/A6 fail)
FS_MAX_PANEL = {"F2b": 2.4}  # 1C3W: A4 and A6 fail from 2.6 (foldmap check)
# the supplementary gallery's labelled panels may go to 3.5 before they are drawn bare
FS_MAX_SUPP = {"S1a": 3.5, "S1b": 3.5, "S1c": 3.5, "S1d": 3.5, "S1e": 3.5}
# where the fit starts (labels mode): font_scale read off a sweep of each panel at 1..4 for its box (the layouts
# change with the text, not smoothly, so a start near the answer avoids a jump to a worse layout)
FS_START = {"F2a": 1.3, "F2b": 2.3, "F2c": 2.8, "F2d": 2.8, "S1b": 2.6, "S1c": 2.8}


def box_of(pid: str) -> tuple[float, float]:
    pid = BOX_AS.get(pid, pid)
    for boxes in LAYOUTS.values():
        for _, p, _x, _y, w, h in boxes:
            if p == pid:
                return w, h - LETTER_ROOM - 0.5
    raise KeyError(pid)


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


_FONT_RE = re.compile(r"font(?:-size)?\s*:\s*(?:[a-z\- ]*\s)?([0-9.]+)px")


def svg_facts(svg: Path) -> dict:
    """Page size (pt), the text sizes (pt, as drawn: matplotlib's SVG user unit is 1 pt) and the fill colours."""
    t = svg.read_text(encoding="utf-8")
    m = re.search(r'<svg[^>]*width="([0-9.]+)pt"[^>]*height="([0-9.]+)pt"', t)
    w, h = float(m.group(1)), float(m.group(2))
    sizes = [float(s) for s in _FONT_RE.findall(t)]
    fills = Counter(c.upper() for c in re.findall(r"fill:\s*(#[0-9a-fA-F]{6})", t))
    return {
        "width_pt": w,
        "height_pt": h,
        "texts": len(sizes),
        "min_font_pt": min(sizes) if sizes else None,
        "font_sizes": sorted(set(round(s, 2) for s in sizes)),
        "fills": dict(fills.most_common()),
        "svg_sha256": hashlib.sha256(svg.read_bytes()).hexdigest(),
    }


def _draw_once(pid: str, args: list[str], tag: str, silent: bool = False) -> dict:
    base = PANELS / "try" / f"{pid}-{tag}"
    base.parent.mkdir(parents=True, exist_ok=True)
    head = [sys.executable, str(HERE / "foldmap_silent.py")] if silent else ["foldmap", "plot"]
    cmd = [*head, *args, "-o", f"{base}.svg", "-o", f"{base}.pdf"]
    p = sh(cmd, cwd=REPO)
    r = {"args": args, "exit": p.returncode, "stderr": p.stderr[-2000:], "stdout": p.stdout[-2000:], "base": str(base)}
    if Path(f"{base}.svg").exists():
        r.update(svg_facts(Path(f"{base}.svg")))
    return r


def printed_min(r: dict, box: tuple[float, float]) -> tuple[float, float | None]:
    s = min(box[0] * MM / r["width_pt"], box[1] * MM / r["height_pt"])
    return s, (r["min_font_pt"] * s if r.get("min_font_pt") else None)


def fit_panel(pid: str) -> dict:
    """Draw pid so that its smallest text prints at MIN_PT or more in its box (see FS_MAX, KEEP_LABELS, BARE)."""
    if pid in HASH_ONLY:
        r = _draw_once(pid, HASH_ONLY[pid], "default")
        r["tries"] = []
        return r
    box = box_of(pid)
    tries = []
    keep = pid in KEEP_LABELS or pid in GALLERY_LABELS
    modes = [("labels", []), ("bare", BARE)]
    for mode, extra in modes:
        fs = FS_START.get(pid, 1.0) if mode == "labels" else 1.0
        fs_max = (
            FS_MAX_PANEL[pid]
            if pid in FS_MAX_PANEL
            else FS_MAX_GALLERY
            if pid in GALLERY_LABELS
            else FS_MAX_KEEP
            if pid in KEEP_LABELS
            else FS_MAX_SUPP.get(pid, FS_MAX_LABELS)
            if mode == "labels"
            else FS_MAX
        )
        for k in range(8 if pid in GALLERY_LABELS else 6):
            args = PANEL_ARGS[pid] + extra + (["--set", f"font_scale={fs:.2f}"] if fs != 1.0 else [])
            r = _draw_once(pid, args, f"{mode}-{k}")
            s, pm = printed_min(r, box) if r.get("width_pt") else (None, None)
            r.update(mode=mode, font_scale=round(fs, 2), scale=s, min_printed_pt=pm)
            tries.append(
                {
                    k2: r.get(k2)
                    for k2 in (
                        "args",
                        "exit",
                        "mode",
                        "font_scale",
                        "scale",
                        "min_font_pt",
                        "min_printed_pt",
                        "width_pt",
                        "height_pt",
                        "texts",
                    )
                }
            )
            if r["exit"] != 0 and not r.get("width_pt"):
                break
            if pm is None or pm >= MIN_PT - 1e-6:
                r["tries"] = tries
                return r
            if k == 0 and pm * fs_max / fs < MIN_PT and not keep:
                break  # the page grows with the text, so not even FS_MAX would reach the floor
            new = fs * MIN_PT / pm * 1.03
            if new > fs_max + 1e-9:
                if fs >= fs_max - 1e-9:
                    break
                new = fs_max
            fs = new
        if keep and mode == "labels":
            # labels are the point of this panel: keep the best labelled try, flagged
            best = max((t for t in tries if t["min_printed_pt"]), key=lambda t: t["min_printed_pt"])
            k = tries.index(best)
            r = _draw_once(pid, best["args"], "kept")
            s, pm = printed_min(r, box)
            r.update(
                mode="labels", font_scale=best["font_scale"], scale=s, min_printed_pt=pm, tries=tries, under_floor=True
            )
            return r
    # no text can print at 6 pt in this box: the bare figure, every text hidden
    r = _draw_once(pid, PANEL_ARGS[pid] + BARE, "silent", silent=True)
    s, pm = printed_min(r, box)
    r.update(mode="silent", font_scale=1.0, scale=s, min_printed_pt=pm, tries=tries)
    if r.get("texts"):
        r["under_floor"] = True
    return r


def determinism_here() -> None:
    """The CI's determinism draw (21 test structures, default view) on this machine, for F3d's fourth hash."""
    out = PANELS / "determinism-remote.json"
    env = dict(os.environ, PYTHONHASHSEED="0", OPENBLAS_NUM_THREADS="1")
    p = sh(
        [sys.executable, "scripts/paper/determinism.py", "--label", os.uname().nodename, "--out", str(out)],
        cwd=REPO,
        env=env,
    )
    print("determinism draw:", p.returncode, p.stdout[-300:], p.stderr[-600:], flush=True)


_STAMP: dict | None = None


def _stamp() -> dict:
    """foldmap.source_stamp of the code that draws (the installed foldmap the CLI runs), taken once as main() starts
    (plan 28 Task 1, review M8) and written with every panel and figure record."""
    global _STAMP
    if _STAMP is None:
        from foldmap.source_stamp import source_stamp

        _STAMP = source_stamp(__file__)
    return _STAMP


def tool_versions() -> dict:
    """The external tools the figures pass through (review M2): tectonic, pdftocairo and pdftoppm (poppler)."""
    out = {}
    for name, cmd in (("tectonic", [str(TECTONIC), "--version"]), ("pdftocairo", ["pdftocairo", "-v"]),
                      ("pdftoppm", ["pdftoppm", "-v"])):  # fmt: skip
        try:
            done = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            out[name] = ((done.stdout or done.stderr).strip().splitlines() or [None])[0]
        except (OSError, subprocess.SubprocessError) as err:
            out[name] = f"unavailable: {type(err).__name__}"
    return out


def draw(ids: list[str]) -> None:
    from concurrent.futures import ThreadPoolExecutor

    PANELS.mkdir(parents=True, exist_ok=True)
    rec_path = PANELS / "panels.json"
    rec = json.loads(rec_path.read_text(encoding="utf-8")) if rec_path.exists() else {}
    with ThreadPoolExecutor(max_workers=min(16, len(ids))) as ex:
        for pid, r in zip(ids, ex.map(fit_panel, ids)):
            base = Path(r.pop("base"))
            for ext in (".svg", ".pdf"):
                if base.with_suffix(ext).exists():
                    shutil.copy(base.with_suffix(ext), PANELS / f"{pid}{ext}")
            how = "python foldmap_silent.py " if r.get("mode") == "silent" else "foldmap plot "
            r["cmd"] = how + " ".join(r["args"]) + f" -o {pid}.svg"
            r["source_stamp"] = _stamp()  # the code that drew this panel
            rec[pid] = r
            print(
                pid,
                r["exit"],
                r.get("mode"),
                "font_scale",
                r.get("font_scale"),
                "printed min",
                r.get("min_printed_pt") and round(r["min_printed_pt"], 2),
                "UNDER" if r.get("under_floor") else "",
                flush=True,
            )
    shutil.rmtree(PANELS / "try", ignore_errors=True)
    for stale in [k for k in rec if k not in PANEL_ARGS]:  # panels no figure uses any more (old gallery ids)
        del rec[stale]
        for ext in (".svg", ".pdf", ".png"):
            (PANELS / f"{stale}{ext}").unlink(missing_ok=True)
    rec_path.write_text(json.dumps(rec, indent=1, sort_keys=True), encoding="utf-8")


def panels() -> dict:
    return json.loads((PANELS / "panels.json").read_text(encoding="utf-8"))


# -------------------------------------------------------------------------------------------------------- assembly
def _tex_escape(s: str) -> str:
    return (
        s.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def assemble(
    name: str, width_mm: float, height_mm: float, items: list[dict], letters: list[tuple], texts: list[tuple] = ()
) -> dict:
    """One figure page, width_mm x height_mm, built by tectonic. items: {pdf, x, y, w | h} in mm from the top left
    (the panel is scaled to width w, or to height h); letters: (letter, x, y) in mm, top left of the letter;
    texts: (text, x, y, size_pt, anchor 'l'|'c'|'r', mono). Returns the scale of each item."""
    build = OUT / "tex" / name
    shutil.rmtree(build, ignore_errors=True)
    build.mkdir(parents=True)
    body, scales = [], {}
    for k, it in enumerate(items):
        src = Path(it["pdf"])
        dst = build / f"p{k}.pdf"
        shutil.copy(src, dst)
        size = it.get("size_pt")  # (w, h) of the panel's page in pt
        if "w" in it:
            opt, scale = f"width={it['w']}mm", it["w"] * MM / size[0]
        else:
            opt, scale = f"height={it['h']}mm", it["h"] * MM / size[1]
        scales[it.get("id", str(k))] = scale
        y_bottom = height_mm - it["y"] - size[1] * scale / MM
        body.append(rf"\put({it['x']:.2f},{y_bottom:.2f}){{\includegraphics[{opt}]{{p{k}.pdf}}}}")
    for letter, x, y in letters:
        body.append(
            rf"\put({x:.2f},{height_mm - y:.2f}){{\makebox(0,0)[lt]{{\fontsize{{{LETTER_PT}}}{{{LETTER_PT}}}"
            rf"\selectfont\bfseries {letter}}}}}"
        )
    for text, x, y, size, anchor, mono in texts:
        fam = r"\ttfamily" if mono else ""
        body.append(
            rf"\put({x:.2f},{height_mm - y:.2f}){{\makebox(0,0)[{anchor}t]{{\fontsize{{{size}}}{{{size * 1.2:.1f}}}"
            rf"\selectfont{fam} {_tex_escape(text)}}}}}"
        )
    tex = rf"""\documentclass{{article}}
\usepackage[paperwidth={width_mm}mm,paperheight={height_mm}mm,margin=0mm]{{geometry}}
\usepackage{{graphicx,fontspec}}
\setmainfont{{{FONT}}}\setmonofont{{Liberation Mono}}
\pagestyle{{empty}}\setlength{{\unitlength}}{{1mm}}\setlength{{\parindent}}{{0pt}}
\begin{{document}}\noindent\begin{{picture}}({width_mm},{height_mm})
{chr(10).join(body)}
\end{{picture}}\end{{document}}
"""
    (build / f"{name}.tex").write_text(tex, encoding="utf-8")
    p = sh([str(TECTONIC), "-X", "compile", "--keep-logs", f"{name}.tex"], cwd=build)
    if p.returncode:
        sys.exit(f"tectonic failed for {name}:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
    shutil.copy(build / f"{name}.pdf", OUT / f"{name}.pdf")
    derive(name)
    return scales


def derive(name: str) -> None:
    """SVG and 300-dpi PNG from the assembled PDF."""
    pdf = OUT / f"{name}.pdf"
    subprocess.run(["pdftocairo", "-svg", str(pdf), str(OUT / f"{name}.svg")], check=True)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile", str(pdf), str(OUT / name)], check=True)


def item(pid: str, x: float, y: float, w: float) -> dict:
    f = panels()[pid]
    return {"id": pid, "pdf": PANELS / f"{pid}.pdf", "x": x, "y": y, "w": w, "size_pt": (f["width_pt"], f["height_pt"])}


def fit(pid: str, w: float, h: float) -> tuple[float, float, float]:
    """The panel scaled to fit a w x h mm box: (w, h, scale) of the drawn panel."""
    f = panels()[pid]
    s = min(w * MM / f["width_pt"], h * MM / f["height_pt"])
    return f["width_pt"] * s / MM, f["height_pt"] * s / MM, s


def box_items(name: str, letters_from: dict | None = None) -> tuple[list[dict], list[tuple]]:
    """The panels of LAYOUTS[name], each centred in its box under the letter room, and the letters."""
    items, letters = [], []
    for letter, pid, x, y, w, h in LAYOUTS[name]:
        bw, bh = w, h - LETTER_ROOM - 0.5
        pw, ph, _ = fit(pid, bw, bh)
        items.append(item(pid, x + (bw - pw) / 2, y + LETTER_ROOM + (bh - ph) / 2, w=pw))
        if letter:  # at the panel's left edge, not the box's
            letters.append(((letters_from or {}).get(pid, letter), max(x, x + (bw - pw) / 2 - 1.0), y + 0.5))
    return items, letters


def box_figure(name: str, texts: list[tuple] = (), letters_from: dict | None = None, extra_h: float = 0.0) -> dict:
    items, letters = box_items(name, letters_from)
    H = round(max(b[3] + b[5] for b in LAYOUTS[name]) + extra_h, 1)
    scales = assemble(name, DOUBLE_MM, H, items, letters, list(texts))
    return {"width_mm": DOUBLE_MM, "height_mm": H, "scales": scales}


# ----------------------------------------------------------------------------------------------- matplotlib setup
def mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [FONT],
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "lines.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    return plt


def mpl_text_floor(fig) -> float:
    from matplotlib.text import Text

    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_text().strip() and t.get_visible()]
    return min(sizes)


def mpl_save(fig, name: str) -> float:
    floor = mpl_text_floor(fig)
    assert floor >= MIN_PT, f"{name}: text at {floor} pt"
    fig.savefig(OUT / f"{name}.pdf")
    derive(name)
    return floor


# ------------------------------------------------------------------------------------------------------- F1
def f1() -> dict:
    """Pipeline schematic: the stages as boxes, the self-check of the drawn figure, its outputs, and a small drawn
    example (1UBQ, F1x)."""
    plt = mpl()
    from matplotlib.patches import FancyBboxPatch

    W, H = DOUBLE_MM, 64.0
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    blue, green, verm = "#0072B2", "#009E73", "#D55E00"

    def box(x, y, w, h, title, body, colour):
        ax.add_patch(
            FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", lw=0.7, ec=colour, fc=colour + "1c")
        )
        ax.text(x + w / 2, y + 1.8, title, ha="center", va="top", fontsize=7, weight="bold", linespacing=1.1)
        ax.text(
            x + w / 2, y + (8.6 if "\n" in title else 5.6), body, ha="center", va="top", fontsize=6, linespacing=1.25
        )

    def path(pts, colour="black", ls="-", head=True):
        xs, ys = zip(*pts)
        ax.plot(xs, ys, color=colour, lw=0.7, ls=ls, solid_capstyle="butt")
        if head:
            (x0, y0), (x1, y1) = pts[-2], pts[-1]
            ax.annotate(
                "",
                (x1, y1),
                (x0 + (x1 - x0) * 0.6, y0 + (y1 - y0) * 0.6),
                arrowprops=dict(arrowstyle="-|>", color=colour, lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0),
            )

    stages = [
        ("Read", "PDB / mmCIF file\nor identifier;\nassembly chosen"),
        ("Secondary\nstructure", "DSSP; elements,\nsheets, barrels,\nbundles"),
        ("Symmetry", "copies found;\nevery copy\ndrawn alike"),
        ("Project", "3D order onto\nthe page; energy\nminimised; overlaps\nremoved"),
        ("Route", "loops orthogonal,\naround elements\nand labels"),
        ("Annotate\nand draw", "ligands, nucleic\nacids, contacts,\nmembrane, domains"),
    ]
    bw, gap, top, bh = 26.4, 3.8, 1.0, 24.0
    xs = [0.5 + i * (bw + gap) for i in range(len(stages))]
    for i, (x, (t, b)) in enumerate(zip(xs, stages)):
        box(x, top, bw, bh, t, b, blue)
        if i:
            path([(xs[i - 1] + bw, top + bh / 2), (x, top + bh / 2)])
    # the self-check, under Symmetry .. Route
    cx, cy, cw, ch = xs[2], 34.0, xs[4] + bw - xs[2], 17.0
    box(
        cx,
        cy,
        cw,
        ch,
        "Self-check of the drawn figure",
        "science S1–S6 against the coordinates:\npairing, barrel closure, residue marks, domains\n"
        "page criteria A1–A20: overlaps, crossings, text size, …",
        verm,
    )
    xa = xs[5] + bw / 2
    path([(xa, top + bh), (xa, 29.5), (xs[4] + bw / 2, 29.5), (xs[4] + bw / 2, cy)])
    # problems are reported with the figure; the author's options redraw it
    xr = xs[3] + bw / 2 - 6
    path([(xr, cy), (xr, top + bh)], colour=verm, ls=(0, (2.5, 1.5)))
    ax.text(
        xr - 1.2,
        29.5,
        "problems reported;\noptions redraw",
        ha="right",
        va="center",
        fontsize=6,
        color=verm,
        linespacing=1.15,
    )
    # outputs, right of the self-check
    ox, ow = xs[5], bw
    box(
        ox,
        cy,
        ow,
        ch + 12.0,
        "Outputs",
        "SVG / PDF / PNG;\nprovenance list:\nwhat was inferred\n(assembly, membrane,\ndomains) or moved;\ncheck report",
        green,
    )
    path([(cx + cw, cy + ch / 2), (ox, cy + ch / 2)], colour=green)
    floor = mpl_text_floor(fig)
    assert floor >= MIN_PT, floor
    fig.savefig(OUT / "tex_F1_base.pdf")
    plt.close(fig)
    # the example: 1UBQ as drawn (F1x), under Read and Secondary structure
    items, _ = box_items("F1")
    base = {"id": "base", "pdf": OUT / "tex_F1_base.pdf", "x": 0, "y": 0, "w": W, "size_pt": (W * MM, H * MM)}
    _, _, x, y, w, h = LAYOUTS["F1"][0]
    scales = assemble(
        "F1",
        W,
        H,
        [base, *items],
        [],
        texts=[
            (
                "Example: ubiquitin (1UBQ)"
                + ("" if panels()["F1x"]["mode"] == "labels" else ", labels omitted at this size"),
                x + w / 2,
                y + 0.6,
                6,
                "c",
                False,
            )
        ],
    )
    (OUT / "tex_F1_base.pdf").unlink()
    return {"width_mm": W, "height_mm": H, "scales": scales, "mpl_floor_pt": floor}


# ------------------------------------------------------------------------------------------------------- F2
def f2() -> dict:
    return box_figure("F2")


def fsgallery() -> dict:
    """Supplementary gallery: 9FWZ, 1EMA, 1TIM and 1KX5, the structures Fig. 2 no longer shows."""
    return box_figure("FSgallery")


# ------------------------------------------------------------------------------------------------------- F3
def f3() -> dict:
    """Panel d: 1TIM (tests/data/1TIM.cif, fitted for the text floor) and, under it, the SHA-256 of its default SVG
    as the CI's determinism run drew it on three systems at cb23732, and as drawn here (F3d-default)."""
    hashes = {}
    for os_ in ("ubuntu-latest", "macos-latest", "windows-latest"):
        d = json.loads((DATA / f"determinism-cb23732/determinism-{os_}/{os_}.json").read_text(encoding="utf-8"))
        row = next(r for r in d["figures"] if r["structure"] == "1TIM.cif")
        hashes[f"{d['system']} {d['machine']} (CI)"] = row["svg_sha256"]
    here = json.loads((PANELS / "determinism-remote.json").read_text(encoding="utf-8"))
    row = next(r for r in here["figures"] if r["structure"] == "1TIM.cif")
    ours = row["svg_sha256"]
    hashes[f"{here['system']} {here['machine']} ({here['label']})"] = ours
    _, _, x, y, w, h = next(b for b in LAYOUTS["F3"] if b[1] == "F3d")
    ty = y + h + 0.5
    texts = [("SHA-256 of the default SVG (foldmap plot 1TIM.cif -o 1TIM.svg) at cb23732:", x, ty, 6, "l", False)]
    for i, (k, v) in enumerate(hashes.items()):
        texts.append((k, x, ty + 3.2 + i * 5.4, 6, "l", False))
        texts.append((v, x, ty + 5.6 + i * 5.4, 6, "l", True))
    info = box_figure("F3", texts, extra_h=0.0)
    info["hashes"] = hashes
    info["identical"] = len(set(hashes.values())) == 1
    return info


# ------------------------------------------------------------------------------------------------------- F4
def f4() -> dict:
    return box_figure("F4")


# ------------------------------------------------------------------------------------------------------- F5
OUTCOMES = [("pass", "#009E73"), ("fail", "#D55E00"), ("not judged", "#E69F00"), ("does not apply", "#BBBBBB")]


def f5_data() -> dict:
    """Everything F5 plots, from the harness's outputs (build/, data/m2/); written to out/F5-data.json."""
    # (a) mkdssp agreement: confusion of the 8 codes over every compared residue
    # rows: chain, number, icode, mkdssp code, foldmap code, partners (dssp_mkdssp.py); as A10 and the
    # pipeline compare them, mkdssp's S and blank read as '-'
    conf = Counter()
    n_entries = 0
    for f in sorted((MS / "data/m2/dssp").glob("*.json.gz")):
        d = json.loads(gzip.decompress(f.read_bytes()))
        n_entries += 1
        for row in d["rows"]:
            mk, fm = (row[3] or "-"), (row[4] or "-")
            mk = "-" if mk in (" ", "S") else mk
            fm = "-" if fm in (" ", "S") else fm
            conf[(fm, mk)] += 1
    seen = {c for k in conf for c in k}
    codes = [c for c in "HGIPEBT-" if c in seen] + sorted(seen - set("HGIPEBT-"))
    # (b) base pairs per entry (as drawn)
    na = na_rows(json.loads((MS / "data/m2/na/report.json").read_text(encoding="utf-8")))
    bp = [
        {
            "entry": e["entry"],
            "precision": e["pair_precision"],
            "recall": e["pair_recall"],
            "pairs": e["pairs_reference"],
        }
        for e in na
        if e.get("pairs_reference")
    ]
    # (c) domains
    dom = json.loads((BUILD / "domains.json").read_text(encoding="utf-8"))
    sets = defaultdict(list)
    for r in dom["reference"]:
        sets["tuning" if r["split"] == "tuning" else "held-out"].append(r)
    sets["fresh"] = dom["fresh"]
    doms = {
        k: [{"pdb": r["pdb"], "chain": r["chain"], "overlap": r["overlap"], "f1": r["f1"]} for r in v]
        for k, v in sets.items()
    }
    # (d) random samples: outcome per criterion
    with open(BUILD / "tables/r100-entries.csv", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    crit = [c for c in rows[0] if re.fullmatch(r"[AS]\d+[a-z]?", c)]
    out = {c: Counter(r[c] for r in rows) for c in crit}
    samples = Counter(r["sample"] for r in rows)
    return {
        "dssp": {"codes": codes, "confusion": {f"{a}{b}": n for (a, b), n in conf.items()}, "entries": n_entries},
        "bp": bp,
        "domains": doms,
        "targets": dom["targets"],
        "r100": {"criteria": crit, "outcomes": {c: dict(v) for c, v in out.items()}, "samples": dict(samples)},
    }


def f5() -> dict:
    d = f5_data()
    (OUT / "F5-data.json").write_text(json.dumps(d, indent=1), encoding="utf-8")
    plt = mpl()
    import numpy as np

    W, H = DOUBLE_MM, 118.0
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))

    # top row: a, b, c; bottom: d over the full width (positions in figure fractions)
    def fx(mm):
        return mm / W

    def fy(mm):
        return 1 - mm / H

    ax_a = fig.add_axes([fx(12), fy(52), fx(42), 46 / H])
    ax_b = fig.add_axes([fx(80), fy(52), fx(38), 46 / H])
    ax_c = fig.add_axes([fx(132), fy(52), fx(44), 46 / H])
    ax_d = fig.add_axes([fx(12), fy(108), fx(164), 42 / H])
    # (a) confusion heat map (log counts)
    codes = d["dssp"]["codes"]
    m = np.array([[d["dssp"]["confusion"].get(f"{a}{b}", 0) for b in codes] for a in codes], float)
    total, agree = m.sum(), np.trace(m)
    from matplotlib.colors import LogNorm

    im = ax_a.imshow(m + 1, cmap="cividis", aspect="equal", norm=LogNorm(vmin=1, vmax=m.max() + 1))
    n = len(codes)
    for i in range(n):
        for j in range(n):
            if m[i, j] and i != j:
                ax_a.text(j, i, f"{int(m[i, j])}", ha="center", va="center", fontsize=6, color="white")
    ax_a.set_xticks(range(n), [c if c != "-" else "–" for c in codes])
    ax_a.set_yticks(range(n), [c if c != "-" else "–" for c in codes])
    ax_a.set_xlabel("mkdssp 3.1.2")
    ax_a.set_ylabel("foldmap")
    ax_a.tick_params(length=0)
    for s in ax_a.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.03)
    cb.set_label("residues", fontsize=6, labelpad=1)
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

    cb.ax.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1_000, 10_000, 100_000]))
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    cb.ax.yaxis.set_minor_formatter(NullFormatter())
    cb.ax.tick_params(labelsize=6, length=2)
    ax_a.set_title(f"{int(agree):,} of {int(total):,} residues agree", fontsize=6.5, loc="left")
    # (b) base pairs per entry
    bp = d["bp"]
    pr = Counter((round(e["precision"], 3), round(e["recall"], 3)) for e in bp)
    xs = [k[1] for k in pr]
    ys = [k[0] for k in pr]
    ss = [8 + 10 * n for n in pr.values()]
    ax_b.scatter(xs, ys, s=ss, color="#0072B2", alpha=0.75, lw=0)
    lo = min(min(xs), min(ys))
    lo = np.floor(lo * 20) / 20
    ax_b.set_xlim(lo - 0.01, 1.01)
    ax_b.set_ylim(lo - 0.01, 1.01)
    ax_b.set_aspect("equal")
    ax_b.set_xlabel("pair recall")
    ax_b.set_ylabel("pair precision")
    n11 = pr.get((1.0, 1.0), 0)
    ax_b.annotate(
        f"{n11} entries at\nprecision = recall = 1",
        (1.0, 1.0),
        xytext=(lo + 0.02, lo + 0.035),
        textcoords="data",
        ha="left",
        va="center",
        fontsize=6,
        arrowprops=dict(arrowstyle="-", color="#555555", lw=0.5, shrinkA=1, shrinkB=6),
    )
    ax_b.text(lo + 0.02, lo + 0.012, "marker area: entries at that point", fontsize=6, color="#555555")
    ax_b.set_title(f"{len(bp)} entries vs RNApolis", fontsize=6.5, loc="left")
    # (c) domains: overlap and boundary F1 per chain, by set
    sets = ["tuning", "held-out", "fresh"]
    rng = np.random.default_rng(20261020)
    for k, metric in enumerate(["overlap", "f1"]):
        for j, s in enumerate(sets):
            v = [r[metric] for r in d["domains"][s]]
            x0 = j * 2.6 + k
            ax_c.scatter(
                x0 + rng.uniform(-0.22, 0.22, len(v)), v, s=6, lw=0, color=["#0072B2", "#E69F00"][k], alpha=0.8
            )
            ax_c.plot([x0 - 0.32, x0 + 0.32], [np.mean(v)] * 2, color="black", lw=0.9)
    ax_c.axhline(d["targets"]["overlap"], color="#0072B2", lw=0.6, ls="--")
    ax_c.axhline(d["targets"]["f1"], color="#E69F00", lw=0.6, ls="--")
    ax_c.set_xticks([j * 2.6 + 0.5 for j in range(3)], [f"{s}\n(n = {len(d['domains'][s])})" for s in sets])
    ax_c.tick_params(axis="x", length=0)
    ax_c.set_ylim(-0.02, 1.05)
    ax_c.set_ylabel("score per chain")
    from matplotlib.lines import Line2D

    ax_c.legend(
        handles=[
            Line2D([], [], marker="o", lw=0, ms=3, color="#0072B2", label="overlap"),
            Line2D([], [], marker="o", lw=0, ms=3, color="#E69F00", label="boundary F1"),
            Line2D([], [], color="black", lw=0.9, label="mean"),
            Line2D([], [], color="grey", lw=0.6, ls="--", label="target"),
        ],
        frameon=False,
        loc="lower left",
        ncol=2,
        handlelength=1.4,
        columnspacing=0.8,
        borderaxespad=0.1,
    )
    ax_c.set_title("inferred domains vs CATH", fontsize=6.5, loc="left")
    # (d) random samples: stacked outcome counts per criterion
    crit = d["r100"]["criteria"]
    bottom = np.zeros(len(crit))
    n_rows = sum(d["r100"]["samples"].values())
    for label, colour in OUTCOMES:
        v = np.array([d["r100"]["outcomes"][c].get(label, 0) for c in crit], float)
        ax_d.bar(range(len(crit)), v, bottom=bottom, color=colour, width=0.75, label=label, lw=0)
        if label == "fail":
            for i, n in enumerate(v):
                if n:
                    ax_d.text(i, n_rows + 3, f"{int(n)}", ha="center", va="bottom", fontsize=6, color="#D55E00")
        bottom += v
    ax_d.set_xticks(range(len(crit)), crit, rotation=0)
    ax_d.tick_params(axis="x", length=0)
    ax_d.set_xlim(-0.6, len(crit) - 0.4)
    ax_d.set_ylim(0, n_rows * 1.12)
    ax_d.set_ylabel("entries")
    ax_d.legend(
        frameon=False,
        ncol=4,
        loc="lower left",
        bbox_to_anchor=(0, 1.04),
        handlelength=1.0,
        columnspacing=1.0,
        borderaxespad=0,
    )
    samples = d["r100"]["samples"]
    names = {"tuned": "tuned", "held": "held-out", "third": "third"}
    ax_d.text(
        1.0,
        1.08,
        f"{n_rows} entries: "
        + ", ".join(f"{names.get(k, k)} {v}" for k, v in samples.items())
        + "; number of failures above each bar",
        transform=ax_d.transAxes,
        ha="right",
        va="bottom",
        fontsize=6,
    )
    for letter, (x, y) in zip("abcd", [(1, 2), (62, 2), (120, 2), (1, 60)]):
        fig.text(fx(x), fy(y), letter, fontsize=LETTER_PT, weight="bold", va="top", ha="left")
    floor = mpl_save(fig, "F5")
    plt.close(fig)
    return {
        "width_mm": W,
        "height_mm": H,
        "mpl_floor_pt": floor,
        "dssp_agree": [int(agree), int(total)],
        "bp_entries": len(bp),
        "bp_at_1_1": n11,
        "r100_rows": n_rows,
    }


# -------------------------------------------------------------------------------------------------- rule checks
def pdf_text_floor(pdf: Path) -> dict:
    """Independent check on the assembled PDF: pdftotext -bbox word heights. A word box is the font's ascent plus
    descent (Liberation Sans 0.905 + 0.212 = 1.117 em), so size = height / 1.117."""
    p = sh(["pdftotext", "-bbox", str(pdf), "-"])
    hs = []
    for m in re.finditer(
        r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', p.stdout
    ):
        x0, y0, x1, y1 = map(float, m.groups()[:4])
        hs.append((min(y1 - y0, x1 - x0) if False else (y1 - y0), m.group(5)))
    hs.sort()
    return {
        "words": len(hs),
        "min_box_pt": hs[0][0] if hs else None,
        "min_font_pt_est": round(hs[0][0] / 1.117, 2) if hs else None,
        "smallest_words": [w for _, w in hs[:8]],
    }


def check(results: dict) -> dict:
    rec = panels()
    out = {}
    for name, info in results.items():
        c = {"size_mm": [info["width_mm"], info["height_mm"]]}
        floors = {}
        for pid, s in info.get("scales", {}).items():
            if pid in rec:
                f = rec[pid]
                floors[pid] = {
                    "scale": round(s, 3),
                    "min_drawn_pt": f["min_font_pt"],
                    "min_printed_pt": round(f["min_font_pt"] * s, 2) if f["min_font_pt"] else None,
                }
        c["panels"] = floors
        if "mpl_floor_pt" in info:
            c["matplotlib_min_pt"] = info["mpl_floor_pt"]
        c["pdf_words"] = pdf_text_floor(OUT / f"{name}.pdf")
        printed = [v["min_printed_pt"] for v in floors.values() if v["min_printed_pt"]]
        c["min_printed_pt"] = min(printed + [info.get("mpl_floor_pt", 99), LETTER_PT])
        c["text_floor_ok"] = c["min_printed_pt"] >= MIN_PT - 1e-6
        out[name] = c
    return out


TARGETS = {
    "f1": ("F1", f1),
    "f2": ("F2", f2),
    "f3": ("F3", f3),
    "f4": ("F4", f4),
    "f5": ("F5", f5),
    "fs": ("FSgallery", fsgallery),
}
GALLERY = ["F2a", "F2b", "F2c", "F2d", "S1a", "S1b", "S1c", "S1d", "S1e"]


def provenance() -> dict:
    """The code the panels were drawn with: the repository state figures.sh read before the run (HEAD, and whether
    src/ or scripts/ differ from it: build/figures/provenance-local.json; a runner syncs the working tree, not a
    commit), foldmap.source_stamp of the foldmap that drew (taken as the process started) and the versions of the
    external tools."""
    p = OUT / "provenance-local.json"
    out = (
        json.loads(p.read_text(encoding="utf-8"))
        if p.exists()
        else {"head": None, "note": "not run through figures.sh"}
    )
    out["source_stamp"] = _stamp()
    out["tools"] = tool_versions()
    return out


def main(argv: list[str]) -> None:
    _stamp()  # the code as this process starts
    OUT.mkdir(parents=True, exist_ok=True)
    what = argv[0] if argv else "all"
    if what == "draw":
        draw(argv[1:] or list(PANEL_ARGS))
        if not argv[1:]:
            determinism_here()
        return
    if what == "gallery":  # Fig. 2 and the supplementary gallery: their panels, then both figures
        draw(GALLERY)
        names = ["f2", "fs"]
    elif what in ("all", "figures"):
        names = list(TARGETS)
    else:
        names = argv
    if what == "all":
        draw(list(PANEL_ARGS))
        determinism_here()
    res_path = OUT / "figures.json"
    res = json.loads(res_path.read_text(encoding="utf-8")) if res_path.exists() else {}
    if what != "check":
        res["_provenance"] = provenance()
    for n in names:
        if n == "check":
            continue
        print("assembling", n, flush=True)
        fig_name, make = TARGETS[n]
        res[fig_name] = make()
        res[fig_name]["provenance"] = res["_provenance"]
    res_path.write_text(json.dumps(res, indent=1, sort_keys=True, default=str), encoding="utf-8")
    checks = check({k: v for k, v in res.items() if not k.startswith("_")})
    (OUT / "checks.json").write_text(json.dumps(checks, indent=1, sort_keys=True), encoding="utf-8")
    for n, c in checks.items():
        print(
            n,
            "min printed text",
            c["min_printed_pt"],
            "pt; ok" if c["text_floor_ok"] else "pt; UNDER 6 pt",
            "| pdf words min ~",
            c["pdf_words"]["min_font_pt_est"],
            "pt",
        )
    if res["_provenance"].get("dirty"):
        print(
            "WARNING: drawn from a working tree that differs from HEAD",
            res["_provenance"]["head"],
            "in src/ or scripts/:",
            " ".join(res["_provenance"].get("changed", [])),
            flush=True,
        )


if __name__ == "__main__":
    main(sys.argv[1:])
