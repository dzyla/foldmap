# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""`foldmap plot ARGS...` with every text of the figure hidden before it is saved: for a panel too small to print
any of its text at 6 pt (make_figures.fit_panel's last step). It runs the installed CLI in this process and wraps
the `save` it calls; repository code is not changed.

    python foldmap_silent.py STRUCTURE [OPTIONS] -o x.svg -o x.pdf
"""

import sys

from matplotlib.text import Text

from foldmap import cli

_save = cli.save


def save(fig, path, *a, **k):
    for t in fig.findobj(Text):
        t.set_visible(False)
    return _save(fig, path, *a, **k)


cli.save = save
sys.exit(cli.main(["plot", *sys.argv[1:]]))
