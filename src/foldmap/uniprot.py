"""UniProt annotation: which UniProt entry each chain is, its features (domains, sites, modifications...) mapped
onto the structure's residues, and the structures known for a protein (PDB entries and the AlphaFold model).

Chains are matched through PDBe's SIFTS residue mapping when the structure is a PDB entry, otherwise by aligning
the UniProt sequence to the chain. Entries and mappings are downloaded once and cached in ~/.cache/foldmap."""

from __future__ import annotations

import json
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from .fetch import AGENT, CACHE

ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"
SEARCH_URL = (
    "https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,protein_name,organism_name&size=5"
)
ACCESSION = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")

REGIONS = (
    "Domain",
    "Region",
    "Repeat",
    "Zinc finger",
    "DNA binding",
    "Coiled coil",
    "Motif",
    "Compositional bias",
    "Transmembrane",
    "Intramembrane",
    "Topological domain",
    "Signal",
    "Propeptide",
    "Transit peptide",
)
SITES = (
    "Active site",
    "Binding site",
    "Site",
    "Modified residue",
    "Lipidation",
    "Cross-link",
    "Glycosylation",
    "Disulfide bond",
)
SITE_SYMBOL = {
    "Active site": "★",
    "Binding site": "◆",
    "Site": "●",
    "Modified residue": "P",
    "Lipidation": "L",
    "Cross-link": "X",
    "Glycosylation": "G",
    "Disulfide bond": "S",
}
MIN_IDENTITY = 0.5  # a chain matched by alignment must be at least this identical to the UniProt sequence


@dataclass
class Feature:
    type: str
    description: str
    start: int  # UniProt numbering, 1-based, inclusive
    end: int
    ligand: str = ""

    @property
    def label(self) -> str:
        text = self.description or self.ligand or self.type
        return text if len(text) <= 40 else text[:38] + "…"


@dataclass
class StructureRef:
    id: str
    method: str
    resolution: float | None
    chains: str  # as UniProt gives it: "A/B=295-610"
    start: int | None = None  # UniProt range covered
    end: int | None = None

    @property
    def coverage(self) -> int:
        return (self.end - self.start + 1) if self.start and self.end else 0


@dataclass
class Entry:
    accession: str
    entry_name: str
    name: str
    gene: str
    organism: str
    sequence: str
    features: list[Feature] = field(default_factory=list)
    structures: list[StructureRef] = field(default_factory=list)
    alphafold: str | None = None

    def regions(self) -> list[Feature]:
        return [f for f in self.features if f.type in REGIONS]

    def sites(self) -> list[Feature]:
        return [f for f in self.features if f.type in SITES]


def _get(url: str, cache: Path) -> dict:
    if cache.is_file():
        return json.loads(cache.read_text())
    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": AGENT, "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = json.loads(r.read().decode())
    except OSError as err:
        raise ValueError(f"could not reach {url.split('/')[2]} ({err})") from None
    cache.write_text(json.dumps(data))
    return data


def parse_entry(data: dict) -> Entry:
    desc = data.get("proteinDescription", {})
    name = (
        (desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
    )
    genes = data.get("genes") or [{}]
    features = []
    for f in data.get("features", []):
        loc = f.get("location", {})
        start, end = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        features.append(
            Feature(
                f["type"], f.get("description", "") or "", int(start), int(end), (f.get("ligand") or {}).get("name", "")
            )
        )
    structures, alphafold = [], None
    for ref in data.get("uniProtKBCrossReferences", []):
        props = {p["key"]: p["value"] for p in ref.get("properties", [])}
        if ref["database"] == "AlphaFoldDB":
            alphafold = ref["id"]
        elif ref["database"] == "PDB":
            res = re.match(r"([\d.]+)", props.get("Resolution", "") or "")
            span = re.search(r"=(\d+)-(\d+)", props.get("Chains", ""))
            structures.append(
                StructureRef(
                    ref["id"],
                    props.get("Method", ""),
                    float(res.group(1)) if res else None,
                    props.get("Chains", ""),
                    *(int(x) for x in span.groups()) if span else (),
                )
            )
    structures.sort(key=lambda s: (-s.coverage, s.resolution if s.resolution is not None else 99.0, s.id))
    return Entry(
        data["primaryAccession"],
        data.get("uniProtkbId", ""),
        name,
        (genes[0].get("geneName") or {}).get("value", ""),
        data.get("organism", {}).get("scientificName", ""),
        data.get("sequence", {}).get("value", ""),
        features,
        structures,
        alphafold,
    )


def entry(accession: str, cache_dir: Path | None = None) -> Entry:
    acc = accession.strip().upper()
    if not ACCESSION.fullmatch(acc):
        raise ValueError(f"{accession!r} is not a UniProt accession (e.g. P04637)")
    folder = Path(cache_dir) if cache_dir else CACHE / "uniprot"
    return parse_entry(_get(ENTRY_URL.format(acc=acc), folder / f"{acc}.json"))


def sifts(pdb_id: str, cache_dir: Path | None = None) -> dict[str, list[dict]]:
    """{accession: [mapping segments]} for a PDB entry, from PDBe SIFTS."""
    pdb = pdb_id.lower()
    folder = Path(cache_dir) if cache_dir else CACHE / "uniprot"
    data = _get(SIFTS_URL.format(pdb=pdb), folder / f"sifts_{pdb}.json")
    return {acc: v.get("mappings", []) for acc, v in data.get(pdb, {}).get("UniProt", {}).items()}


def pdb_id_of(path) -> str | None:
    """The PDB ID a structure file belongs to (its _entry.id), or None for models and other files."""
    import gemmi

    try:
        name = gemmi.read_structure(str(path)).name.strip()
    except (RuntimeError, ValueError):
        return None
    return name.upper() if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", name) else None


@dataclass
class ChainMap:
    chain: str
    accession: str
    to_label: dict[int, int]  # UniProt position -> label_seq_id of the chain


def chain_maps(path, accession: str | None = None, cache_dir: Path | None = None) -> list[ChainMap]:
    """UniProt entry and residue mapping for each protein chain of the deposited structure. With an accession,
    every chain that matches it (by alignment); otherwise SIFTS for PDB entries, else the accession named in the
    file (AlphaFold models), each chain aligned to it."""
    from .sequence import _global, chain_track

    path = Path(path)
    out: list[ChainMap] = []
    pdb = pdb_id_of(path)
    if accession is None and pdb is not None:
        for acc, segments in sifts(pdb, cache_dir).items():
            for seg in segments:
                start, ustart = seg["start"]["residue_number"], seg["unp_start"]
                n = seg["unp_end"] - ustart + 1
                m = next((c for c in out if c.chain == seg["chain_id"] and c.accession == acc), None)
                if m is None:
                    m = ChainMap(seg["chain_id"], acc, {})
                    out.append(m)
                for k in range(n):
                    m.to_label[ustart + k] = start + k
        return out
    acc = accession or _accession_in_file(path)
    if acc is None:
        raise ValueError(f"{path.name}: not a PDB entry and no UniProt accession given; use --uniprot ACCESSION")
    e = entry(acc, cache_dir)
    from .io import load_backbone

    for chain in dict.fromkeys(l.chain for l in load_backbone(path, "asu").labels):
        track = chain_track(path, chain, full_sequence=True)
        pairs = _global(e.sequence, track.letters)
        same = sum(e.sequence[i] == track.letters[j] for i, j in pairs)
        if same >= MIN_IDENTITY * min(len(track.letters), len(e.sequence)):
            out.append(ChainMap(chain, e.accession, {i + 1: j + 1 for i, j in pairs}))
    return out


def _accession_in_file(path: Path) -> str | None:
    """A UniProt accession recorded in the file (AlphaFold models: AF-P04637-F1, _ma_target_ref_db_details)."""
    m = re.search(r"AF-(" + ACCESSION.pattern + r")-F", path.name.upper())
    if m:
        return m.group(1)
    try:
        import gemmi

        block = gemmi.cif.read(str(path)).sole_block()
        for tag in ("_ma_target_ref_db_details.db_accession", "_struct_ref.pdbx_db_accession"):
            values = list(block.find_values(tag))
            if values and ACCESSION.fullmatch(gemmi.cif.as_string(values[0])):
                return gemmi.cif.as_string(values[0])
    except (RuntimeError, ValueError, IndexError):
        pass
    return None


@dataclass
class Annotation:
    """Features mapped onto a Backbone: regions as residue ranges, sites as residues, per chain."""

    entries: dict[str, Entry]  # accession -> entry
    chains: dict[str, str]  # backbone chain -> accession
    regions: list[tuple[Feature, str, list[int]]] = field(default_factory=list)  # (feature, chain, residue indices)
    sites: list[tuple[Feature, str, list[int]]] = field(default_factory=list)


def annotate(path, bb, accession: str | None = None, cache_dir: Path | None = None) -> Annotation:
    """Map every chain's UniProt features onto bb's residues (assembly copies of a chain take its features)."""
    from .sequence import chain_track

    maps = chain_maps(path, accession, cache_dir)
    entries = {m.accession: entry(m.accession, cache_dir) for m in maps}
    index = {(l.chain, l.seq, l.icode): k for k, l in enumerate(bb.labels)}
    present = list(dict.fromkeys(l.chain for l in bb.labels))
    ann = Annotation(entries, {})
    from .io import load_backbone

    asu = load_backbone(path, "asu")
    for m in maps:
        track = chain_track(path, m.chain, full_sequence=True)
        # UniProt position -> (author number, insertion code) via the chain's label_seq positions
        label_to_auth = {}
        for p, r in enumerate(track.residue):
            if r is not None:
                label_to_auth[p + 1] = (asu.labels[r].seq, asu.labels[r].icode)
        copies = [c for c in present if c == m.chain or re.fullmatch(re.escape(m.chain) + r"-?\d+", c)]
        for copy in copies:
            ann.chains[copy] = m.accession

            def residues(start: int, end: int, copy=copy, m=m, label_to_auth=label_to_auth) -> list[int]:
                out = []
                for u in range(start, end + 1):
                    auth = label_to_auth.get(m.to_label.get(u, -1))
                    if auth is not None and (copy, *auth) in index:
                        out.append(index[(copy, *auth)])
                return out

            e = entries[m.accession]
            for f in e.regions():
                ks = residues(f.start, f.end)
                if ks:
                    ann.regions.append((f, copy, ks))
            for f in e.sites():
                ends = (f.start, f.end) if f.type == "Disulfide bond" else range(f.start, f.end + 1)
                ks = [k for u in ends for k in residues(u, u)]
                if ks:
                    ann.sites.append((f, copy, ks))
    return ann


def report(source: str) -> str:
    """Text for an accession, or for a structure: which entry each chain is, then each entry described."""
    if ACCESSION.fullmatch(source.strip().upper()):
        return describe(entry(source))
    from .fetch import fetch

    path = fetch(source)
    maps = chain_maps(path)
    if not maps:
        return f"{source}: no chain matched a UniProt entry"
    lines = [f"chain {m.chain} → {m.accession}" for m in maps]
    for acc in dict.fromkeys(m.accession for m in maps):
        lines += ["", describe(entry(acc))]
    return "\n".join(lines)


def describe(e: Entry, limit: int = 12) -> str:
    """A text summary: names, regions, sites, structures and the AlphaFold model."""
    lines = [f"{e.accession} ({e.entry_name}) · {e.name} · {e.gene or '-'} · {e.organism} · {len(e.sequence)} aa"]
    if e.regions():
        lines.append("Domains and regions:")
        lines += [f"  {f.type:<18} {f.start:>5}-{f.end:<5} {f.label}" for f in e.regions()]
    sites = e.sites()
    if sites:
        lines.append(f"Sites and modifications ({len(sites)}):")
        lines += [
            f"  {f.type:<18} {f.start:>5}{'-' + str(f.end) if f.end != f.start else '':<6} {f.label}"
            for f in sites[:limit]
        ]
        if len(sites) > limit:
            lines.append(f"  … {len(sites) - limit} more")
    if e.structures:
        lines.append(f"Experimental structures ({len(e.structures)}), widest coverage first:")
        for s in e.structures[:limit]:
            res = f"{s.resolution:.2f} Å" if s.resolution is not None else "-"
            lines.append(f"  {s.id}  {s.method:<6} {res:>8}  {s.chains}")
        if len(e.structures) > limit:
            lines.append(f"  … {len(e.structures) - limit} more")
    lines.append(f"AlphaFold model: {'AF-' + e.alphafold + '-F1' if e.alphafold else 'none'}")
    return "\n".join(lines)
