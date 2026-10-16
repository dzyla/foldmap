# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""Published domains: CATH's chopping of a PDB entry, through the PDBe API (`--domains cath`, plan 12 Task 3).

Each CATH domain is one or more segments of a chain in author residue numbers; a panel takes the elements with most
of their residues inside its segments, named by its CATH domain id (1ldmA01). Assembly copies of a chain (A1, A-2)
take its domains. Responses are downloaded once and cached in ~/.cache/foldmap/cath; an entry PDBe has no CATH
chopping for (a 404) is cached as such. CATH data: CC BY 4.0 (https://www.cathdb.info)."""

from __future__ import annotations

import http.client
import json
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .fetch import AGENT, CACHE, atomic_write

CATH_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/cath/{pdb}"
SOURCE = "CATH"  # Layout.domain_source of these domains
NO_CHOPPING_MAX_AGE = 30 * 24 * 3600  # seconds a cached "no CATH chopping" (404) is trusted: CATH grows


@dataclass
class CathDomain:
    domain: str  # CATH domain id, e.g. 1ldmA01
    code: str  # CATH superfamily, e.g. 3.40.50.720
    topology: str  # e.g. Rossmann fold
    chain: str  # author chain id
    segments: list[tuple[tuple[int, str], tuple[int, str]]]  # ((first number, icode), (last number, icode))


def chopping(pdb_id: str, cache_dir: Path | None = None) -> list[CathDomain]:
    """CATH's domains of a PDB entry, in chain and domain order; [] when PDBe has no CATH chopping for it."""
    pdb = pdb_id.lower()
    folder = Path(cache_dir) if cache_dir else CACHE / "cath"
    cache = folder / f"cath_{pdb}.json"
    data = None
    if cache.is_file():
        try:
            data = json.loads(cache.read_text(encoding="utf-8"))
        except (ValueError, UnicodeDecodeError):  # cut off: fetch it again
            cache.unlink(missing_ok=True)
        if data == {} and time.time() - cache.stat().st_mtime > NO_CHOPPING_MAX_AGE:
            data = None  # no chopping a month ago: ask again
    if data is None:
        url = CATH_URL.format(pdb=pdb)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                body = r.read()
            try:
                data = json.loads(body.decode())
            except (ValueError, UnicodeDecodeError):
                raise ValueError(
                    f"PDBe's CATH answer for {pdb} is not JSON (the service may be down); try again later"
                ) from None
        except urllib.error.HTTPError as err:
            if err.code != 404:
                raise ValueError(f"could not reach www.ebi.ac.uk ({err})") from None
            data = {}  # PDBe: "Requested endpoint does not contain any data"
        except (OSError, http.client.HTTPException) as err:
            raise ValueError(f"could not reach www.ebi.ac.uk ({err})") from None
        atomic_write(cache, json.dumps(data))
    found: dict[str, CathDomain] = {}
    for code, family in (data.get(pdb, {}).get("CATH") or {}).items():
        for m in sorted(family.get("mappings", []), key=lambda m: (m.get("domain", ""), m.get("segment_id", 0))):
            start, end = m["start"], m["end"]
            d = found.setdefault(
                m["domain"], CathDomain(m["domain"], code, family.get("topology", ""), m["chain_id"], [])
            )
            d.segments.append(
                (
                    (int(start["author_residue_number"]), start.get("author_insertion_code") or ""),
                    (int(end["author_residue_number"]), end.get("author_insertion_code") or ""),
                )
            )
    return sorted(found.values(), key=lambda d: (d.chain, d.domain))


def domains(path, bb, sses, cache_dir: Path | None = None) -> tuple[list[tuple[str, list[str]]], str]:
    """(domain panels [(CATH domain id, element ids)], a note for the user) for a structure file: its PDB id's CATH
    chopping (chopping) laid on the chains of `bb`, each element in the domain holding most of its residues (none
    when no domain holds half of them). With several chains split, a panel's name ends in " · chain"."""
    from .uniprot import pdb_id_of

    pdb = pdb_id_of(path)
    if pdb is None:
        return [], "no CATH domains: the structure is not a PDB entry (CATH chops PDB entries); drawn without domains"
    found = chopping(pdb, cache_dir)
    if not found:
        return [], f"no CATH chopping for {pdb} at PDBe; drawn without domains (--domains auto infers them)"
    order = {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}
    of_chain: dict[str, list[int]] = {}
    for k, l in enumerate(bb.labels):
        of_chain.setdefault(l.chain, []).append(k)
    present = list(of_chain)
    chopped = {d.chain for d in found}  # an author chain CATH chops itself is nobody's copy (an "A1" beside "A")
    panels = []
    for d in found:
        copies = [
            c for c in present if c == d.chain or (c not in chopped and re.fullmatch(re.escape(d.chain) + r"-?\d+", c))
        ]
        for copy in copies:
            inside: set[int] = set()
            for (a, ai), (b, bi) in d.segments:
                lo, hi = order.get((copy, a, ai)), order.get((copy, b, bi))
                if lo is not None and hi is not None:  # in file order: insertion codes fall where they stand
                    inside.update(range(lo, hi + 1))
                else:  # an end residue not modelled: by number
                    inside.update(k for k in of_chain[copy] if a <= bb.labels[k].seq <= b)
            ids = [
                s.id
                for s in sses
                if s.chain == copy and sum(k in inside for k in range(s.start, s.end + 1)) > (s.end - s.start + 1) / 2
            ]
            if ids:
                panels.append((d.domain, copy, ids))
    split = list(dict.fromkeys(c for _, c, _ in panels))
    named = [(name if len(split) == 1 else f"{name} · {c}", ids) for name, c, ids in panels]
    listed = ", ".join(f"{d.domain} ({d.code}{', ' + d.topology if d.topology else ''})" for d in found)
    if not named:
        return [], f"CATH domains of {pdb} ({listed}) hold no drawn element; drawn without domains"
    return named, f"domains from CATH ({pdb}, via PDBe): {listed}"
