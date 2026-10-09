"""Structures by name: PDB IDs from RCSB, UniProt accessions / AlphaFold DB ids from the AlphaFold DB, local
paths, or uploaded bytes. Downloads are cached in ~/.cache/foldmap."""

from __future__ import annotations

import hashlib
import re
import urllib.request
from pathlib import Path

CACHE = Path.home() / ".cache" / "foldmap"
AFDB_VERSION = 6  # AlphaFold DB model version to fetch
AGENT = "foldmap (https://pypi.org/project/foldmap)"  # the AlphaFold DB refuses requests without a user agent


def fetch(source: str, upload: tuple[str, bytes] | None = None) -> Path:
    """A local structure file from an upload, a path, a PDB ID (from RCSB) or a UniProt accession / AlphaFold DB
    id (the AlphaFold model); downloads are cached."""
    if upload is not None:
        name, data = upload
        CACHE.mkdir(parents=True, exist_ok=True)
        path = CACHE / f"upload-{hashlib.sha1(data).hexdigest()[:12]}{Path(name).suffix or '.cif'}"
        path.write_bytes(data)
        return path
    text = (source or "").strip()
    if text and Path(text).expanduser().is_file():
        return Path(text).expanduser()
    af = re.fullmatch(
        r"(?:AF-)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-F1)?", text.upper()
    )
    if af:  # a UniProt accession or AlphaFold DB id: the AlphaFold model
        acc = af.group(1)
        CACHE.mkdir(parents=True, exist_ok=True)
        path = CACHE / f"AF-{acc}-F1.cif"
        if not path.is_file():
            url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{AFDB_VERSION}.cif"
            try:
                with urllib.request.urlopen(
                    urllib.request.Request(url, headers={"User-Agent": AGENT}), timeout=30
                ) as r:
                    path.write_bytes(r.read())
            except OSError as err:
                raise ValueError(f"could not download AF-{acc} from the AlphaFold DB ({err})") from None
        return path
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", text):
        CACHE.mkdir(parents=True, exist_ok=True)
        path = CACHE / f"{text.upper()}.cif"
        if not path.is_file():
            try:
                with urllib.request.urlopen(f"https://files.rcsb.org/download/{text.upper()}.cif", timeout=30) as r:
                    path.write_bytes(r.read())
            except OSError as err:
                raise ValueError(f"could not download {text.upper()} from the PDB ({err})") from None
        return path
    raise ValueError(f"{text!r} is not a structure file, a PDB ID (e.g. 1LMB) or a UniProt accession (e.g. P04637)")
