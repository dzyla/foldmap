"""Structures by name: PDB IDs from RCSB, UniProt accessions / AlphaFold DB ids from the AlphaFold DB, local
paths, or uploaded bytes. Downloads are cached in ~/.cache/foldmap."""

from __future__ import annotations

import hashlib
import os
import re
import tempfile
import urllib.request
from pathlib import Path

CACHE = Path.home() / ".cache" / "foldmap"
AFDB_VERSION = 6  # AlphaFold DB model version to fetch
AGENT = "foldmap (https://pypi.org/project/foldmap)"  # the AlphaFold DB refuses requests without a user agent


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
        path.unlink()
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": AGENT}), timeout=30) as r:
            data = r.read()
    except OSError as err:
        raise ValueError(f"could not download {what} ({err})") from None
    if not _complete(data):
        raise ValueError(f"could not download {what} (the file arrived incomplete; try again)")
    atomic_write(path, data)
    return path


def fetch(source: str, upload: tuple[str, bytes] | None = None) -> Path:
    """A local structure file from an upload, a path, a PDB ID (from RCSB) or a UniProt accession / AlphaFold DB
    id (the AlphaFold model); downloads are cached."""
    if upload is not None:
        name, data = upload
        path = CACHE / f"upload-{hashlib.sha1(data).hexdigest()[:12]}{Path(name).suffix or '.cif'}"
        atomic_write(path, data)
        return path
    text = (source or "").strip()
    if text and Path(text).expanduser().is_file():
        return Path(text).expanduser()
    af = re.fullmatch(
        r"(?:AF-)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-F1)?", text.upper()
    )
    if af:  # a UniProt accession or AlphaFold DB id: the AlphaFold model
        acc = af.group(1)
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{AFDB_VERSION}.cif"
        return _download(url, CACHE / f"AF-{acc}-F1.cif", f"AF-{acc} from the AlphaFold DB")
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", text):
        pdb = text.upper()
        return _download(f"https://files.rcsb.org/download/{pdb}.cif", CACHE / f"{pdb}.cif", f"{pdb} from the PDB")
    raise ValueError(f"{text!r} is not a structure file, a PDB ID (e.g. 1LMB) or a UniProt accession (e.g. P04637)")
