# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Structures by name: PDB IDs from RCSB, UniProt accessions / AlphaFold DB ids from the AlphaFold DB, local
paths, or uploaded bytes. Downloads are cached in ~/.cache/foldmap: at most DOWNLOAD_MAX bytes each, kept a week
after their last use and DOWNLOAD_CACHE_MAX bytes in all."""

from __future__ import annotations

import hashlib
import http.client
import os
import re
import stat
import tempfile
import threading
import time
import urllib.request
from pathlib import Path

CACHE = Path.home() / ".cache" / "foldmap"
AFDB_VERSION = 6  # AlphaFold DB model version to fetch
AGENT = "foldmap (https://pypi.org/project/foldmap)"  # the AlphaFold DB refuses requests without a user agent
STRUCTURE_SUFFIXES = (".cif", ".mmcif", ".pdb", ".ent")
ALIGNMENT_SUFFIXES = (".fasta", ".fa", ".aln", ".sto", ".txt")
UPLOAD_MAX_AGE = 24 * 3600  # seconds an uploaded file is kept in the cache
MB = 1024 * 1024
DOWNLOAD_MAX = 100 * MB  # bytes a download may take (the largest PDB entry, 3J3Q, is 242 MB as mmCIF)
CHUNK = MB  # bytes read at a time: a download is never held whole in memory
DOWNLOAD_MAX_AGE = 7 * 24 * 3600  # seconds a downloaded structure is kept unused in the cache
DOWNLOAD_CACHE_MAX = 2 * 1024 * MB  # bytes of downloaded structures kept; the oldest go first
PRUNE_EVERY = 3600  # seconds: the cache is pruned at most this often, from the download path
_UPLOAD_NAME = re.compile(r"upload-[0-9a-f]{12}\.[a-z]+")
_AF = re.compile(r"(?:AF-)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-F1)?")
_PDB = re.compile(r"[0-9][A-Za-z0-9]{3}")
_DOWNLOADED = re.compile(r"[0-9][A-Z0-9]{3}\.cif|AF-[A-Z0-9]+-F1\.cif")  # the cache's names for downloads


def atomic_write(path: Path, data: bytes | str) -> None:
    """Write a file so no reader ever sees it half-written: a temporary file beside it, then a rename.
    Text is written as UTF-8."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = data.encode("utf-8") if isinstance(data, str) else data
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".part")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(raw)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _complete(data: bytes) -> bool:
    """Whether a downloaded mmCIF file is whole: RCSB and AlphaFold DB files end with a '#' line."""
    return data.rstrip().endswith(b"#")


def _download(url: str, path: Path, what: str) -> Path:
    """`url` into the cache at `path`, unless a whole copy is already there (a cut-off one is fetched again). The
    download streams to a file beside `path`, stopped past DOWNLOAD_MAX bytes; only a whole file takes its place.
    The cache is pruned first, at most once an hour (prune_downloads)."""
    maybe_prune()
    if path.is_file():
        with path.open("rb") as f:
            f.seek(max(0, path.stat().st_size - 64))
            if _complete(f.read()):
                os.utime(path)  # used now: pruned last
                return path
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, part = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".part")
    try:
        with os.fdopen(fd, "wb") as out:
            tail = _stream(url, out, what)
        if not _complete(tail):
            raise ValueError(f"could not download {what} (the file arrived incomplete; try again)")
        os.replace(part, path)
    except BaseException:
        Path(part).unlink(missing_ok=True)
        raise
    return path


def _stream(url: str, out, what: str) -> bytes:
    """The download of `url` written to `out`, a chunk at a time; its last bytes. ValueError when it fails or
    passes DOWNLOAD_MAX bytes (told by its Content-Length at once, else counted as it arrives)."""
    too_big = (
        f"could not download {what}: it is larger than {DOWNLOAD_MAX // MB} MB, more than foldmap downloads. "
        "Download the file yourself and draw it on your own computer: foldmap plot FILE"
    )
    tail, size = b"", 0
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": AGENT}), timeout=30) as r:
            length = (getattr(r, "headers", None) or {}).get("Content-Length")
            if length and length.isdigit() and int(length) > DOWNLOAD_MAX:
                raise ValueError(too_big)
            while chunk := r.read(CHUNK):
                size += len(chunk)
                if size > DOWNLOAD_MAX:
                    raise ValueError(too_big)
                out.write(chunk)
                tail = (tail + chunk)[-64:]
    except (OSError, http.client.HTTPException) as err:
        raise ValueError(f"could not download {what} ({err})") from None
    return tail


def _downloads() -> list[Path]:
    """The structures foldmap downloaded into the cache (PDB entries, AlphaFold models): regular files only."""
    found = [*CACHE.glob("*.cif"), *CACHE.glob("afdb-v*/AF-*-F1.cif")]
    return [p for p in found if _DOWNLOADED.fullmatch(p.name) and p.is_file() and not p.is_symlink()]


def prune_downloads(now: float | None = None) -> list[Path]:
    """Remove downloads older than DOWNLOAD_MAX_AGE (since last used), then the oldest until the rest take at
    most DOWNLOAD_CACHE_MAX bytes; what was removed. Uploads have their own pruning (keep_upload)."""
    now = time.time() if now is None else now
    files = []
    for p in _downloads():
        try:
            info = p.stat()
        except OSError:
            continue
        files.append((info.st_mtime, info.st_size, p))
    files.sort()  # oldest first
    total = sum(size for _, size, _ in files)
    removed = []
    for mtime, size, p in files:
        if now - mtime <= DOWNLOAD_MAX_AGE and total <= DOWNLOAD_CACHE_MAX:
            break  # this and every newer file stay
        p.unlink(missing_ok=True)
        total -= size
        removed.append(p)
    return removed


_PRUNED = [0.0]  # when the download cache was last pruned (time.monotonic), process-wide
_PRUNE_LOCK = threading.Lock()


def maybe_prune() -> None:
    """prune_downloads, unless this process pruned in the last PRUNE_EVERY seconds."""
    with _PRUNE_LOCK:
        if _PRUNED[0] and time.monotonic() - _PRUNED[0] < PRUNE_EVERY:
            return
        _PRUNED[0] = time.monotonic()
    try:
        prune_downloads()
    except OSError:  # a cache that cannot be listed: the download goes ahead
        pass


def _prune_uploads(folder: Path) -> None:
    """Remove uploads older than UPLOAD_MAX_AGE. Only regular files named like foldmap's own uploads go; symlinks,
    other files and everything outside the folder are left alone."""
    if folder.is_symlink():
        return
    now = time.time()
    for old in folder.iterdir():
        if not _UPLOAD_NAME.fullmatch(old.name):
            continue
        try:
            info = old.lstat()
            if stat.S_ISREG(info.st_mode) and now - info.st_mtime > UPLOAD_MAX_AGE:
                old.unlink(missing_ok=True)
        except OSError:
            continue


def keep_upload(name: str, data: bytes, suffixes: tuple[str, ...] = STRUCTURE_SUFFIXES) -> Path:
    """An uploaded file, saved under the cache's uploads/ folder by content hash (so the same upload is the same
    path). Old uploads are removed first. The name only lends its suffix, when allowed."""
    folder = CACHE / "uploads"
    folder.mkdir(parents=True, exist_ok=True)
    _prune_uploads(folder)
    suffix = Path(name).suffix.lower()
    target = folder / f"upload-{hashlib.sha1(data).hexdigest()[:12]}{suffix if suffix in suffixes else suffixes[0]}"
    atomic_write(target, data)
    return target


def fetch(source: str, upload: tuple[str, bytes] | None = None, *, local: bool = True) -> Path:
    """A local structure file from an upload, a path, a PDB ID (from RCSB) or a UniProt accession / AlphaFold DB
    id (the AlphaFold model); downloads are cached. local=False never reads a path given as source (for a hosted
    app, whose visitors must not read the server's files)."""
    if upload is not None:
        return keep_upload(*upload)
    text = (source or "").strip()
    if local and text and Path(text).expanduser().is_file():
        return Path(text).expanduser()
    af = _AF.fullmatch(text.upper())
    if af:  # a UniProt accession or AlphaFold DB id: the AlphaFold model
        acc = af.group(1)
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{AFDB_VERSION}.cif"
        folder = CACHE / f"afdb-v{AFDB_VERSION}"  # a new AlphaFold DB release is a new download
        return _download(url, folder / f"AF-{acc}-F1.cif", f"AF-{acc} from the AlphaFold DB")
    if _PDB.fullmatch(text):
        pdb = text.upper()
        return _download(f"https://files.rcsb.org/download/{pdb}.cif", CACHE / f"{pdb}.cif", f"{pdb} from the PDB")
    what = "a structure file, a PDB ID" if local else "a PDB ID"
    raise ValueError(f"{text!r} is not {what} (e.g. 1LMB) or a UniProt accession (e.g. P04637)")


def source_id(text: str) -> str | None:
    """The id a share link may carry for this source: a PDB ID, UniProt accession or AlphaFold DB id (upper-cased,
    as typed otherwise); None for anything else (a path, a file name, an upload)."""
    text = (text or "").strip().upper()
    return text if _AF.fullmatch(text) or _PDB.fullmatch(text) else None
