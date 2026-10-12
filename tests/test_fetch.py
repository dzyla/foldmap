"""Downloads and the cache: only whole files are kept, and a broken cache entry is fetched again. Offline."""

import http.client
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


class _Dropped:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self):
        raise http.client.IncompleteRead(b"abc", 100)


def test_connection_cut_mid_download_gives_a_friendly_error_and_caches_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch_mod, "CACHE", tmp_path)
    monkeypatch.setattr(fetch_mod.urllib.request, "urlopen", lambda request, timeout=0: _Dropped())
    with pytest.raises(ValueError, match="could not download"):
        fetch("1UBQ")
    assert list(tmp_path.iterdir()) == []


def test_connection_cut_mid_uniprot_request_gives_a_friendly_error(monkeypatch, tmp_path):
    monkeypatch.setattr(uniprot, "CACHE", tmp_path)
    monkeypatch.setattr(uniprot.urllib.request, "urlopen", lambda request, timeout=0: _Dropped())
    with pytest.raises(ValueError, match="could not reach"):
        uniprot.entry("P69905")


def test_cut_off_cache_file_is_replaced_without_deleting_it_first(monkeypatch, tmp_path):
    (tmp_path / "1UBQ.cif").write_bytes(UBQ[:3000])
    _serve(monkeypatch, tmp_path, UBQ)

    def no_unlink(self, *a, **k):
        raise AssertionError("the cut-off file must be replaced atomically, not unlinked")

    monkeypatch.setattr(Path, "unlink", no_unlink)
    assert fetch("1UBQ").read_bytes() == UBQ


def test_uploads_are_kept_by_content_and_old_ones_removed(monkeypatch, tmp_path):
    import os
    import time

    monkeypatch.setattr(fetch_mod, "CACHE", tmp_path)
    (tmp_path / "1UBQ.cif").write_bytes(UBQ)  # a downloaded structure, not an upload
    a = fetch("", ("x.cif", UBQ))
    assert a == fetch("", ("renamed.cif", UBQ)) and a.parent == tmp_path / "uploads"
    odd = fetch("", ("../../evil.sh", b"data_x\n"))
    assert odd.parent == tmp_path / "uploads" and odd.suffix == ".cif"
    stranger = tmp_path / "uploads" / "notes.txt"  # not made by foldmap
    stranger.write_text("mine")
    outside = tmp_path / "outside.txt"
    outside.write_text("keep me")
    link = tmp_path / "uploads" / "upload-000000000000.cif"
    link.symlink_to(outside)
    old = time.time() - fetch_mod.UPLOAD_MAX_AGE - 60
    for f in (a, stranger, tmp_path / "1UBQ.cif"):
        os.utime(f, (old, old))
    os.utime(link, (old, old), follow_symlinks=False)
    os.utime(outside, (old, old))
    fetch("", ("y.pdb", b"other"))
    assert not a.exists() and odd.exists()
    assert stranger.exists() and (tmp_path / "1UBQ.cif").exists() and outside.exists()
    assert link.is_symlink()  # not made by foldmap as a file: left alone


def test_alignment_uploads_keep_only_allowed_suffixes(monkeypatch, tmp_path):
    monkeypatch.setattr(fetch_mod, "CACHE", tmp_path)
    ok = fetch_mod.keep_upload("a.STO", b">x\nAA\n", fetch_mod.ALIGNMENT_SUFFIXES)
    odd = fetch_mod.keep_upload("a.py", b">x\nAA\n", fetch_mod.ALIGNMENT_SUFFIXES)
    assert ok.suffix == ".sto" and odd.suffix == ".fasta"


def test_paths_can_be_switched_off_for_a_hosted_app():
    path = DATA / "1UBQ.cif"
    assert fetch(str(path)) == path
    with pytest.raises(ValueError, match="is not a PDB ID"):
        fetch(str(path), local=False)


P53 = DATA / "AF-P04637.cif"


def _offline(monkeypatch, tmp_path, asked: list):
    class Reply:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self):
            return P53.read_bytes()

    def fake(request, timeout=0):
        asked.append(getattr(request, "full_url", request))
        return Reply()

    monkeypatch.setattr(fetch_mod, "CACHE", tmp_path)
    monkeypatch.setattr(fetch_mod.urllib.request, "urlopen", fake)


def test_alphafold_models_are_cached_per_release(monkeypatch, tmp_path):
    asked = []
    _offline(monkeypatch, tmp_path, asked)
    first = fetch("P04637")
    assert first.name == "AF-P04637-F1.cif" and first.parent.name == "afdb-v6"
    monkeypatch.setattr(fetch_mod, "AFDB_VERSION", 7)
    second = fetch("P04637")
    assert second != first and len(asked) == 2 and asked[-1].endswith("-model_v7.cif")


def test_share_ids_are_ids_only():
    from foldmap.fetch import source_id

    assert source_id(" 1lmb ") == "1LMB" and source_id("p04637") == "P04637"
    assert source_id("AF-P04637-F1") == "AF-P04637-F1"
    assert source_id("/tmp/x.cif") is None and source_id("hello") is None and source_id("") is None
