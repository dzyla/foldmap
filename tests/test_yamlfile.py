"""User YAML (layout, names and style files) is read with no anchors or aliases: a few hundred bytes of nested
aliases ("billion laughs") would otherwise expand to gigabytes in the process that reads them."""

import time
import tracemalloc

import pytest

from foldmap import layoutfile
from foldmap.names import parse_names
from foldmap.style import resolve_style
from foldmap.yamlfile import read_yaml


def laughs(levels: int = 8) -> str:
    """Nested aliases: 10**levels ones from a few hundred bytes (the final review's laughs8.yaml, 415 bytes)."""
    nest = "&a0 [1,1,1,1,1,1,1,1,1,1]"
    for i in range(1, levels):
        nest = f"&a{i} [{nest}, " + ", ".join([f"*a{i - 1}"] * 9) + "]"
    return nest


def refused_quickly(read, *args) -> str:
    """read(*args) raises ValueError at once, holding next to no memory while it does; its message."""
    tracemalloc.start()
    t0 = time.perf_counter()
    try:
        with pytest.raises(ValueError) as info:
            read(*args)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert time.perf_counter() - t0 < 1.0
    assert peak < 2_000_000, peak  # bytes: the expansion would be 10**8 list items
    return str(info.value)


def test_the_laughs_layout_file_is_refused_at_once_without_growing(tmp_path):
    f = tmp_path / "laughs8.yaml"
    f.write_text(f"foldmap: 1\nlayout:\n  protomers: {laughs()}\n")
    assert len(f.read_bytes()) < 450
    said = refused_quickly(layoutfile.load, f)
    assert said == f"{f}: aliases are not allowed in a layout file (no &anchor or *alias)"


def test_names_and_style_files_refuse_aliases_too(tmp_path):
    said = refused_quickly(parse_names, f"elements: {{res:A:1-5: {laughs()}}}\n", "n.yaml")
    assert said.startswith("n.yaml: aliases are not allowed in a names file")
    style = tmp_path / "s.yaml"
    style.write_text(f"mark: {laughs()}\n")
    assert "aliases are not allowed in a style file" in refused_quickly(resolve_style, None, style)


@pytest.mark.parametrize(
    "text",
    [
        "a: &x 1\n",  # an anchor alone
        "a: &x 1\nb: *x\n",  # and its alias
        "base: &b {fill: pale}\nlayout:\n  <<: *b\n",  # a merge key needs an alias
    ],
)
def test_any_anchor_or_alias_is_refused(text):
    with pytest.raises(ValueError, match="w.yaml: aliases are not allowed in a layout file"):
        read_yaml(text, "layout", "w.yaml")


def test_plain_yaml_reads_as_before():
    assert read_yaml("", "names", "w.yaml") is None
    assert read_yaml("a: [1, 'x&y', '*z']\nb: {c: off}\n", "layout", "w.yaml") == {
        "a": [1, "x&y", "*z"],
        "b": {"c": False},
    }
    with pytest.raises(ValueError, match="w.yaml: not valid YAML"):
        read_yaml("a: [1\n", "layout", "w.yaml")
    with pytest.raises(ValueError, match="w.yaml: nested too deeply"):
        read_yaml("a: " + "[" * 3000 + "]" * 3000 + "\n", "layout", "w.yaml")
