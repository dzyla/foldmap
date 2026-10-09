"""The project website: built from site/ by scripts/build_site.py, every figure rendered by foldmap."""

import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _builder():
    spec = importlib.util.spec_from_file_location("build_site", ROOT / "scripts" / "build_site.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_quick_build_has_every_referenced_file(tmp_path):
    out = tmp_path / "site"
    data = _builder().build(out, quick=True)
    index = (out / "index.html").read_text()
    assert "{{" not in index and (out / "app.js").is_file() and (out / ".nojekyll").is_file()
    embedded = json.loads(re.search(r'id="site-data">(.*?)</script>', index, re.S).group(1))
    assert embedded == json.loads(json.dumps(data))
    paths = [f for s in data["structures"] for f in s["figures"].values()]
    paths += [p["src"] for m in data["modes"] for p in m["panels"]] + [x["src"] for x in data["explorers"]]
    assert paths and all((out / p).is_file() for p in paths)
    assert all((out / p).read_text().lstrip().startswith(("<?xml", "<svg", "<!doctype")) for p in paths)


def test_every_theme_has_a_note():
    from foldmap.style import THEMES

    assert set(THEMES) <= set(_builder().THEME_NOTES)
