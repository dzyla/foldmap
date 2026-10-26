"""Metamorphic test helpers (plan 28 Task 7, review 3.6): rewrite a committed fixture under a transform that should not
change what foldmap reads from it, and read back what foldmap says, in a form two runs can be compared in.

Transforms (`rewrite`): a rigid motion (a proper rotation and a translation; the assembly operators carried with it),
chain renaming, the atoms of each residue shuffled, the file written as mmCIF or PDB (by gemmi), and Gaussian noise
on every coordinate. The exact motion is a cyclic permutation of the axes and a translation by whole Ångströms, so
every coordinate is written back with the digits it was read with: an exact transform (structure facts asserted
equal). A general rotation is not exact once written (coordinates rounded to 0.001 Å): measured, not asserted.

What is compared (`facts`): the structure facts, chains named as strings "CHAIN:..." throughout so a renaming maps
back by its prefix (`renamed`); `geometry`: the drawn layout's element and stem boxes, bounds and loop points."""

from __future__ import annotations

from pathlib import Path

import gemmi
import numpy as np

DATA = Path(__file__).parent / "data"
CYCLE = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)  # (x, y, z) -> (z, x, y): det +1, exact in decimal
SHIFT = (17.0, -23.0, 41.0)  # Å, whole numbers: exact once written with three decimals


def turn(axis, degrees: float) -> np.ndarray:
    """A proper rotation about `axis` (Rodrigues)."""
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    a = np.radians(degrees)
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(a) * kx + (1 - np.cos(a)) * kx @ kx


def _transform(rot: np.ndarray, shift) -> gemmi.Transform:
    tr = gemmi.Transform()
    tr.mat.fromlist(np.asarray(rot, float).tolist())
    tr.vec.fromlist([float(x) for x in shift])
    return tr


def _move(st: gemmi.Structure, rot: np.ndarray, shift) -> None:
    """Every atom moved by x -> rot x + shift; each assembly operator conjugated (T op T^-1) so the assembly it
    builds is the moved assembly."""
    rot = np.asarray(rot, float)
    shift = np.asarray(shift, float)
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    p = rot @ np.array([atom.pos.x, atom.pos.y, atom.pos.z]) + shift
                    atom.pos = gemmi.Position(*p.tolist())
    tr = _transform(rot, shift)
    inv = tr.inverse()
    for assembly in st.assemblies:
        for gen in assembly.generators:
            for op in gen.operators:
                op.transform = tr.combine(op.transform).combine(inv)


def _shuffle(st: gemmi.Structure, seed: int) -> None:
    """The atoms of every residue in a random order (residues and chains stay in order)."""
    rng = np.random.default_rng(seed)
    for model in st:
        for chain in model:
            for res in chain:
                atoms = [res[k].clone() for k in range(len(res))]
                for k in range(len(res) - 1, -1, -1):
                    del res[k]
                for k in rng.permutation(len(atoms)):
                    res.add_atom(atoms[int(k)])


def _noise(st: gemmi.Structure, sigma: float, seed: int) -> None:
    rng = np.random.default_rng(seed)
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    d = rng.normal(0.0, sigma, 3)
                    atom.pos = gemmi.Position(atom.pos.x + d[0], atom.pos.y + d[1], atom.pos.z + d[2])


def rewrite(
    src: Path,
    out: Path,
    *,
    rot=None,
    shift=(0.0, 0.0, 0.0),
    rename: dict[str, str] | None = None,
    shuffle: int | None = None,
    noise: tuple[float, int] | None = None,
) -> Path:
    """`src` read by gemmi, transformed, written to `out` (.cif: mmCIF, .pdb: PDB). With no transform it is the
    writer's round trip, the baseline each transform is compared with."""
    st = gemmi.read_structure(str(src))
    if rot is not None:
        _move(st, rot, shift)
    for old in rename or {}:  # through a temporary name, so a swap (A->B, B->A) works
        st.rename_chain(old, f"~{old}")
    for old, new in (rename or {}).items():
        st.rename_chain(f"~{old}", new)
    if shuffle is not None:
        _shuffle(st, shuffle)
    if noise is not None:
        _noise(st, *noise)
    out = Path(out)
    if out.suffix == ".pdb":
        st.write_pdb(str(out))
    else:
        st.make_mmcif_document().write_file(str(out))
    return out


# ------------------------------------------------------------------------------------------------ what is compared


def _key(label) -> str:
    return f"{label.chain}:{label.seq}{label.icode}"


def facts(path: Path, nucleic: bool = False, membrane: bool = False) -> dict:
    """What foldmap reads from the structure as deposited (assembly asu): residues, DSSP codes, bridges, elements,
    sheets (order, directions, pair kinds, closure, accessory strands), base pairs (nucleic), the transmembrane set
    (membrane). Chains appear only as the prefix of "CHAIN:..." strings."""
    from foldmap.dssp import assign_dssp
    from foldmap.io import load_backbone, load_nucleic
    from foldmap.sheets import build_sheets
    from foldmap.ss import build_sses

    bb = load_backbone(path, "asu")
    keys = [_key(l) for l in bb.labels]
    d = assign_dssp(bb)
    sses = build_sses(bb, d.ss, bridges=d.bridges)
    out: dict = {
        "residues": keys,
        "ss": d.ss,
        "bridges": sorted([keys[b.i], keys[b.j], b.kind] for b in d.bridges),
        "elements": [[s.kind, s.id] for s in sses],
        "sheets": [
            [
                [s.id for s in sh.strands],
                list(sh.directions),
                list(sh.pair_kinds),
                sh.closed,
                sorted([a.strand.id, a.parent.id, a.column] for a in sh.accessory),
            ]
            for sh in build_sheets(sses, d.bridges)
        ],
    }
    if nucleic:
        na = load_nucleic(path, bb, "asu")
        out["base_pairs"] = sorted(
            [_key(na.strands[p.a].labels[p.i]), _key(na.strands[p.b].labels[p.j]), p.kind] for p in na.base_pairs
        )
    if membrane:
        from foldmap.cli import _transmembrane
        from foldmap.membrane import find_membrane

        mem = find_membrane(path, bb)
        tm = set() if mem is None else _transmembrane(sses, bb, mem)
        out["transmembrane"] = None if mem is None else sorted(tm)
        out["transmembrane_residues"] = [keys[r] for s in sses if s.id in tm for r in range(s.start, s.end + 1)]
    return out


def _chain(name: str, mapping: dict[str, str]) -> str:
    """A chain name mapped; an assembly copy's name (the deposited name plus the operator's number, A1) by its
    deposited part."""
    if name in mapping:
        return mapping[name]
    stem = name.rstrip("0123456789")
    return mapping[stem] + name[len(stem) :] if stem and stem in mapping else name


def renamed(x, mapping: dict[str, str]):
    """x with the chain of every "CHAIN:..." string mapped (each part of a "CHAIN:.../CHAIN:..." stem id), and the
    lists that are sorted by name ("bridges", "base_pairs") sorted again, so two runs compare equal under a
    renaming exactly when they agree."""
    if isinstance(x, str):
        if "/" in x:
            return "/".join(renamed(part, mapping) for part in x.split("/"))
        head, sep, tail = x.partition(":")
        return f"{_chain(head, mapping)}{sep}{tail}" if sep else x
    if isinstance(x, dict):
        out = {renamed(k, mapping): renamed(v, mapping) for k, v in x.items()}
        for key in ("bridges", "base_pairs"):
            if isinstance(out.get(key), list):
                out[key] = sorted(out[key])
        return out
    if isinstance(x, list):
        return [renamed(v, mapping) for v in x]
    if isinstance(x, tuple):
        return tuple(renamed(v, mapping) for v in x)
    return x


def verdicts(drawn, spec, path) -> dict:
    """`foldmap check`'s verdicts on a drawn figure: each science check's findings (None: not measured) and each
    criterion's pass (True, False or None)."""
    from foldmap.criteria import check_figure

    found = check_figure(drawn, spec, str(path))
    return {
        "science": {k: (None if v is None else sorted(v)) for k, v in found["science"].items()},
        "criteria": {k: v["pass"] for k, v in found["criteria"].items()},
    }


def _box(p) -> list:
    """(cx, cy, cos, sin of the angle, length, width): an angle of pi and one of -pi are one direction."""
    return [p.cx, p.cy, float(np.cos(p.angle)), float(np.sin(p.angle)), p.length, p.width]


def geometry(drawn) -> dict:
    """The drawn layout's geometry: per element id its box (_box), ghosts, stems, bounds, sheets, and each loop's
    points (by its two element ids)."""
    lay = drawn.layout
    out = {
        "placed": {k: _box(p) for k, p in sorted(lay.placed.items())},
        "ghosts": [[g.label, *_box(g)] for g in lay.ghosts],
        "stems": [[d.id, *_box(d)] for d in lay.stems],
        "bounds": list(lay.bounds),
        "sheets": [list(x) for x in lay.sheets],
        "loops": [[lp.a_id, lp.b_id, [list(map(float, q)) for q in lp.points]] for lp in drawn.loops],
    }
    return out


def _numbers(x) -> list[float]:
    if isinstance(x, (int, float, np.floating, np.integer)) and not isinstance(x, bool):
        return [float(x)]
    if isinstance(x, (list, tuple)):
        return [v for y in x for v in _numbers(y)]
    if isinstance(x, dict):
        return [v for k in sorted(x) for v in _numbers(x[k])]
    return []


def _shape(x):
    """x with every number replaced by 0: what must agree exactly (ids, labels, counts)."""
    if isinstance(x, (int, float, np.floating, np.integer)) and not isinstance(x, bool):
        return 0
    if isinstance(x, (list, tuple)):
        return [_shape(y) for y in x]
    if isinstance(x, dict):
        return {k: _shape(v) for k, v in sorted(x.items())}
    return x


def geometry_gap(a: dict, b: dict) -> float:
    """The largest difference between two geometries' numbers (page units, degrees); inf when their ids, labels or
    counts differ."""
    if _shape(a) != _shape(b):
        return float("inf")
    x, y = np.array(_numbers(a)), np.array(_numbers(b))
    return float(np.abs(x - y).max()) if len(x) else 0.0
