"""Speed-ups that must not change results: each fast path checked against the plain one it replaced."""

import json
import os
import subprocess
import sys
from pathlib import Path

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).parent.parent

_SIGNATURE = """
import json, sys
from foldmap.cli import make_layout
out = {}
for f in sys.argv[1:]:
    lay, _, _ = make_layout(f)
    out[f] = {k: [p.cx, p.cy, p.angle] for k, p in sorted(lay.placed.items())}
print(json.dumps(out))
"""


def test_layout_does_not_depend_on_the_hash_seed():
    files = [str(DATA / f"{p}.cif") for p in ("2HHB", "1LMB", "9TX7")]
    runs = []
    for seed in ("0", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        done = subprocess.run(
            [sys.executable, "-c", _SIGNATURE, *files], env=env, capture_output=True, text=True, check=True
        )
        runs.append(json.loads(done.stdout))
    assert runs[0] == runs[1]
