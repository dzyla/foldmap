# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Structures by name: PDB IDs from RCSB, UniProt accessions / AlphaFold DB ids from the AlphaFold DB, local
paths, or uploaded bytes. Downloads are cached in ~/.cache/foldmap."""

from __future__ import annotations

import hashlib
import http.client
import os
import re
import stat
import tempfile
import time
import urllib.request
from pathlib import Path

CACHE = Path.home() / ".cache" / "foldmap"
AFDB_VERSION = 6  # AlphaFold DB model version to fetch
AGENT = "foldmap (https://pypi.org/project/foldmap)"  # the AlphaFold DB refuses requests without a user agent
STRUCTURE_SUFFIXES = (".cif", ".mmcif", ".pdb", ".ent")
ALIGNMENT_SUFFIXES = (".fasta", ".fa", ".aln", ".sto", ".txt")
UPLOAD_MAX_AGE = 24 * 3600  # seconds an uploaded file is kept in the cache
_UPLOAD_NAME = re.compile(r"upload-[0-9a-f]{12}\.[a-z]+")


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
    """`url` into the cache at `path`, unless a whole copy is already there (a cut-off one is fetched again)."""
    if path.is_file():
        with path.open("rb") as f:
            f.seek(max(0, path.stat().st_size - 64))
            if _complete(f.read()):
                return path
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": AGENT}), timeout=30) as r:
            data = r.read()
    except (OSError, http.client.HTTPException) as err:
        raise ValueError(f"could not download {what} ({err})") from None
    if not _complete(data):
        raise ValueError(f"could not download {what} (the file arrived incomplete; try again)")
    atomic_write(path, data)
    return path


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
    af = re.fullmatch(
        r"(?:AF-)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-F1)?", text.upper()
    )
    if af:  # a UniProt accession or AlphaFold DB id: the AlphaFold model
        acc = af.group(1)
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{AFDB_VERSION}.cif"
        folder = CACHE / f"afdb-v{AFDB_VERSION}"  # a new AlphaFold DB release is a new download
        return _download(url, folder / f"AF-{acc}-F1.cif", f"AF-{acc} from the AlphaFold DB")
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", text):
        pdb = text.upper()
        return _download(f"https://files.rcsb.org/download/{pdb}.cif", CACHE / f"{pdb}.cif", f"{pdb} from the PDB")
    what = "a structure file, a PDB ID" if local else "a PDB ID"
    raise ValueError(f"{text!r} is not {what} (e.g. 1LMB) or a UniProt accession (e.g. P04637)")
