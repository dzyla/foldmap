# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Draw the app's mark (src/foldmap/assets/icon.png and icon.svg): a three-turn helix ribbon in Okabe-Ito blue
beside a strand arrow in Okabe-Ito orange, on a transparent 64 × 64 canvas. The page icon of `foldmap app` and the
mark on its Start screen. Reruns write the same bytes (no date, a fixed SVG hash salt). Re-run after changing it:

    python scripts/make_app_icon.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

ASSETS = Path(__file__).resolve().parents[1] / "src" / "foldmap" / "assets"
BLUE, ORANGE = "#0072B2", "#E69F00"
SIZE = 64  # px, and the drawing's units


def mark():
    """The figure: the helix (back half faded, front half on top) and the strand arrow, in pixel units."""
    fig = plt.figure(figsize=(1, 1), dpi=SIZE)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, SIZE)
    ax.set_ylim(0, SIZE)
    ax.set_axis_off()
    t = np.linspace(0, 3 * 2 * np.pi, 541)  # three turns, seen from the side: each turn an ellipse's halves
    x, y = 20 + 11 * np.sin(t), 10 + 44 * t / t[-1] + 5 * np.cos(t)
    front = np.cos(t) >= 0
    edges = np.flatnonzero(np.diff(front.astype(int))) + 1
    for faded in (True, False):  # the back half of each turn first, the front half over it
        for run in np.split(np.arange(len(t)), edges):
            if front[run[0]] != faded:
                ends = np.r_[run, run[-1] + 1] if run[-1] + 1 < len(t) else run
                ax.plot(x[ends], y[ends], color=BLUE, alpha=0.4 if faded else 1.0, lw=6, solid_capstyle="round")
    arrow = [(42, 6), (52, 6), (52, 42), (60, 42), (47, 59), (34, 42), (42, 42)]
    ax.add_patch(Polygon(arrow, closed=True, facecolor=ORANGE, edgecolor="none"))
    return fig


def main() -> None:
    plt.rcParams["svg.hashsalt"] = "foldmap-icon"
    fig = mark()
    fig.savefig(ASSETS / "icon.png", dpi=SIZE, transparent=True, metadata={"Software": None})
    fig.savefig(ASSETS / "icon.svg", transparent=True, metadata={"Date": None})
    plt.close(fig)
    print("wrote", ASSETS / "icon.png", "and", ASSETS / "icon.svg")


if __name__ == "__main__":
    main()
