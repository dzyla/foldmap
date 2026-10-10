# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Release guard: the git tag must name the version the package carries (v0.2.0 or 0.2.0 for 0.2.0).

Usage: python scripts/check_release.py "$GITHUB_REF_NAME"
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from foldmap import __version__  # noqa: E402


def check(tag: str, version: str = __version__) -> str | None:
    """None when the tag matches the version, otherwise why it does not."""
    bare = tag.removeprefix("refs/tags/").removeprefix("v")
    return None if bare == version else f"tag {tag!r} does not match foldmap.__version__ {version!r}"


if __name__ == "__main__":
    problem = check(sys.argv[1])
    if problem:
        sys.exit(problem)
    print(f"tag {sys.argv[1]} matches foldmap {__version__}")
