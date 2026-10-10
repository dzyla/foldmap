"""Downloads and the cache: only whole files are kept, and a broken cache entry is fetched again. Offline."""

import json
from pathlib import Path

import pytest

import foldmap.fetch as fetch_mod
from foldmap import uniprot
from foldmap.fetch import atomic_write, fetch

DATA = Path(__file__).parent / "data"
UBQ = (DATA / "1UBQ.cif").read_bytes()


class _Reply:
    def __init__(self, payload: bytes):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self):
        return self.payload


def _serve(monkeypatch, cache: Path, payload: bytes) -> list[str]:
    """Point the cache at `cache` and answer every download with `payload`; returns the URLs asked for."""
    asked = []

    def fake(request, timeout=0):
        asked.append(getattr(request, "full_url", request))
        return _Reply(payload)

    monkeypatch.setattr(fetch_mod, "CACHE", cache)
    monkeypatch.setattr(fetch_mod.urllib.request, "urlopen", fake)
    return asked


def test_atomic_write_writes_text_as_utf8_and_leaves_nothing_else(tmp_path):
    target = tmp_path / "sub" / "x.json"
    atomic_write(target, '{"name": "α-helix"}')
    assert json.loads(target.read_bytes().decode("utf-8")) == {"name": "α-helix"}
    assert [p.name for p in target.parent.iterdir()] == ["x.json"]


def test_failed_atomic_write_keeps_the_old_file(tmp_path, monkeypatch):
    target = tmp_path / "x.cif"
    target.write_bytes(b"old")

    def disk_full(src, dst):
        raise OSError("No space left on device")

    monkeypatch.setattr(fetch_mod.os, "replace", disk_full)
    with pytest.raises(OSError):
        atomic_write(target, b"new")
    assert target.read_bytes() == b"old"
    assert [p.name for p in tmp_path.iterdir()] == ["x.cif"]


def test_pdb_entry_is_downloaded_once_then_cached(monkeypatch, tmp_path):
    asked = _serve(monkeypatch, tmp_path, UBQ)
    path = fetch("1ubq")
    assert path == tmp_path / "1UBQ.cif" and path.read_bytes() == UBQ
    assert fetch("1UBQ") == path
    assert asked == ["https://files.rcsb.org/download/1UBQ.cif"]


def test_cut_off_cache_file_is_downloaded_again(monkeypatch, tmp_path):
    (tmp_path / "1UBQ.cif").write_bytes(UBQ[:3000])  # left by an interrupted download
    asked = _serve(monkeypatch, tmp_path, UBQ)
    assert fetch("1UBQ").read_bytes() == UBQ
    assert len(asked) == 1


def test_incomplete_download_is_refused_and_not_cached(monkeypatch, tmp_path):
    _serve(monkeypatch, tmp_path, UBQ[:3000])
    with pytest.raises(ValueError, match="incomplete"):
        fetch("1UBQ")
    assert not (tmp_path / "1UBQ.cif").exists()


def test_error_page_served_as_a_structure_is_refused(monkeypatch, tmp_path):
    _serve(monkeypatch, tmp_path, b"<html><body>Service temporarily unavailable</body></html>")
    with pytest.raises(ValueError, match="incomplete"):
        fetch("1UBQ")
    assert list(tmp_path.iterdir()) == []


def test_broken_uniprot_cache_entry_is_fetched_again(monkeypatch, tmp_path):
    good = (DATA / "uniprot" / "P69905.json").read_bytes()
    folder = tmp_path / "uniprot"
    folder.mkdir()
    (folder / "P69905.json").write_bytes(good[:500])
    asked = []

    def fake(request, timeout=0):
        asked.append(request.full_url)
        return _Reply(good)

    monkeypatch.setattr(uniprot, "CACHE", tmp_path)
    monkeypatch.setattr(uniprot.urllib.request, "urlopen", fake)
    assert uniprot.entry("P69905").sequence
    assert len(asked) == 1
    assert json.loads((folder / "P69905.json").read_bytes()) == json.loads(good)
