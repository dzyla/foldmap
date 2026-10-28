# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dawid Zyla. Additional terms (AGPL-3.0 section 7(b)): see NOTICE.
"""A plan 29 fix's short record under ruling P29-13 (the lighter loop, Tasks 6 on): a fix that must change no figure,
measured on the tune half and the bench only, plus its own measure, written as report.json and report.md, stamped
with the code that wrote them (foldmap.source_stamp).

    python scripts/app29_light.py helical --old FILE --out FILE [--jobs N] [--max-copies 64]
    python scripts/app29_light.py science --old FILE --out FILE [--jobs N]
    python scripts/app29_light.py campaign --before DIR --after DIR --label LABEL --out FILE [--na]
    python scripts/app29_light.py navalidate --before DIR --after DIR --out FILE
    python scripts/app29_light.py report --name NAME --before RUN --after RUN --bench FILE --measure FILE... --out DIR

helical (Task 6) runs symmetry._helical as it was (--old: that symmetry.py, e.g. `git show <commit>:src/foldmap/
symmetry.py`) and as it is on the same _Assembly of every structure: every .cif/.cif.gz under tests/data, plan 29's
tune half and every cached download (~/.cache/foldmap/*.cif) less plan 29's check half, which it never names. It
compares the two returns exactly (floats bit for bit), counts _Assembly.pair calls and times each, with both
MIN_HELIX_COPIES and MIN_HELIX_COPIES_ASKED. A structure whose largest cluster holds more than --max-copies copies is
skipped for time (the old search is n^4). Exit 1 when any return differs.

science (Task 7) draws plan 29's tune half and every .cif/.cif.gz under tests/data as `foldmap check` does
(spec.draw_spec(spec.figure_spec(path)): the default view) and judges each drawing with figcheck.science_checks as it
was (--old: that figcheck.py) and as it is: S1-S6 before and after, per entry, and the entries whose findings
changed, and whether any finding is new.

campaign (Task 8b) reads two scripts/random100.py run folders (--views default; report.json) of the same entries,
before and after a fix that may move figures: each entry's campaign keys (sensitivity.keys_of: pass, the failing
criteria, the named fallbacks, the figure's area) judged as scripts/prune.py judges them (campaign_verdict: worse,
better or moved), the entries whose SVG moved, and with --na each entry's A17 measures (na_fixes_record._row, every
entry of the run, not only that script's list) judged as na_fixes_record.compare judges them. LABEL names the set
(the third sample, the NA controls).

navalidate (Task 8b) reads two scripts/na_validate.py reports (the NA campaign, docs/validation/na) and compares each
entry's judged measures (na_fixes_record.JUDGED) the same way, with each entry's pass.

report reads two scripts/app29.py run records of the tune half (load_tune: report.json alone; RUN, a folder with it),
the bench's --check --exact output (bench.txt, then "bench exit N") and the fix's own measure (--measure, a JSON
file this script's other commands wrote; it is copied into the record), and says: every tune entry whose figure (the
app's SVG and PNG, sha256) or criteria changed, the draw seconds of the entries the fix names (FIXES' named), the
bench's verdict, and the measure's summary. The verdict is P29-13's: kept only when every figure is byte-identical,
the bench is identical on every entry, and the measure changed only what the fix may change (FIXES' may_change),
with no new finding. A fix that may move figures (FIXES' explained: entry -> why, Task 8b on) is kept when nothing
got worse (a tune entry's criteria, the bench's exit, a campaign key or an NA measure) and every figure that moved
(tune, bench, campaign) is one the fix explains. --measure takes several files."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "docs" / "validation" / "app29" / "sample.json"

# each fix: what it is, the entries whose draw seconds the record names, and what it was expected to change
FIXES: dict[str, dict] = {
    "fix-screw": {
        "title": "Plan 29 Task 6: the screw search without n^4 pair fits (triage proposal 3, ruling P29-7)",
        "mechanism": "symmetry._helical tried every pair of copies as the screw and tested each on every copy against "
        "every member of its cluster: n^4/2 pair fits (23FR's 42 copies: 1.26 million _Assembly.pair calls). A fit's "
        "RMSD is at least the distance between the moved and the target centroids of its matched atoms, so a member "
        "whose centroid lands farther than the tolerance (plus 1e-6 Å for rounding) is passed over unfitted; the "
        "first member within the tolerance, in cluster order, is the one the full search took.",
        "named": ["23FR"],
        "expected": "23FR's time only; every figure byte-identical",
        "may_change": [],
    },
    "fix-ghost": {
        "title": "Plan 29 Task 7: a lent ghost is no barrel strand (triage proposal 4, ruling P29-7)",
        "mechanism": "figcheck.science_checks counted every ghost lent to a barrel block (layout._unwrap: a "
        "neighbouring copy's strand drawn across an unrolled ring's cut leaves its ghost where it pairs) as a ring "
        "strand. 9X4L's sheet 7 (chain 3's 5-strand barrel) holds the ghost of 7:7-8, lent beside the accessory "
        "3:208-213 and paired with it alone: ring 6 against science.pairing_cycle's 5. science.barrel_finding now "
        "counts a lent ghost in the ring only when it is joined to a ring strand in the pairing graph (pairing_cycle's "
        "own edge rule, science._joined), directly or through another such ghost (science.lent_ring); the others "
        "are counted as accessories are. A check fix: no figure may change.",
        "named": ["9X4L"],
        "expected": "9X4L's S3 only (foldmap check; the app does not show S3); every figure byte-identical",
        "may_change": ["9X4L"],
        "control": "No barrel block in tests/data or the tune half holds a ring strand lent across a cut: the "
        "measles F sheets (8UTF) and the other blocks with lent ghosts (8UUP, 1SLV, 4V4F-ring, 14CU) are drawn "
        "open, judged by science.open_finding, which this fix leaves as it was. The guard is a synthetic ring "
        "(tests/test_science.py::test_a_ring_strand_lent_across_the_cut_still_counts); the barrels whose S3 passes "
        "(1EMA, 1TIM, 1BXW-A: tests/test_science.py::test_a_control_barrel_figure_passes_and_a_moved_strand_fails) "
        "and 8UTF's S3 are pinned.",
    },
    "fix-copies": {
        "title": "Plan 29 Task 8b: copies modelled in part are one protein to focus auto (26EL, ruling P29-14)",
        "mechanism": "cli._one_of_each, which picks focus auto's chains (one of each distinct protein in a protomer), "
        "grouped chains by the whole modelled sequence (difflib ratio > 0.9). A copy modelled in part falls under it: "
        "26EL's protomer holds A (1-682), F (1-291) and M (292-682), ratios 0.59 and 0.73, so A, F and M were drawn, "
        "132-301 Å apart on a diagonal, said as '1 of 7 subunits'. Two chains are now also one protein when "
        "model.homolog_chains (the copy rule the labels already share by: at least 20 residues at the same numbers, "
        "90 % alike, the numbering offset allowed) gives them one root. No new constant.",
        "named": ["26EL"],
        "expected": "26EL's figure (A alone, '1 of 20 subunits shown'); elsewhere only symmetric entries whose "
        "protomer holds a copy modelled in part",
        "explained": {
            "26EL": "the named entry: focus auto draws A alone (58 elements) instead of A, F and M on a diagonal; the "
            "note reads '1 of 20 subunits shown · C2 symmetry' (was '1 of 7'); the page 60.0 x 75.3 -> 27.9 x 32.9 "
            "units (diagnosis-t8.md); no criterion moved",
            "7SPI": "a C13 type IV secretion core whose protomer holds two copies of TraV, A1 (150-204) and B1 "
            "(91-203), each modelled in part (ratio 0.90, not above it): both were drawn in full, with the same labels "
            "B-E twice. They are one protein by their numbers (54 residues at the same numbers, all alike), so A1 is "
            "drawn once and B1 only where it completes A1's sheet (D' and E', grey); the note says 1 of 26 subunits "
            "(78 chains over A1, C1, E1; C1/D1 and E1/F1 were already one protein each, by the ratio), was 1 of 20. "
            "Area 426.01 -> 338.03 (no judge orders it); pass, criteria and fallbacks identical",
        },
        "control": "No near-identical paralogs sharing residue numbering exist among the 710 structures at hand "
        "(tests/data, the tune half, the zl download cache less the check half): every pair of chains of two "
        "polymer entities alike at 70 % or more of 20 or more shared numbers is 96.7 % alike or more, and "
        "homolog_chains joins each (so no real entry shows it keeping such paralogs apart). Three of those pairs "
        "the fix newly joins (whole-sequence ratio 0.9 or less): HIV-1 reverse transcriptase's p66 and p51 in 8VBG "
        "and 1HYS (one gene, p51 its first 440 residues, folded otherwise) and 6V6S's poly-UNK chains of different "
        "entities (unassigned poly-alanine models, their residues all UNK). None is symmetric, so focus auto never "
        "reaches _one_of_each on them, and their figures do not change (8VBG and 1HYS are NA controls: identical). "
        "Pinned as it is: tests/test_focus.py::test_one_protein_by_its_numbers_is_today_one_protein_even_as_two_"
        "entities (8VBG captured); the rule's edge on a synthetic pair (85 % alike: two proteins): "
        "test_near_identical_paralogs_below_the_copy_rule_stay_two_proteins. The risk, named: a symmetric assembly "
        "of p66/p51 heterodimers, or of poly-UNK models of different proteins, would draw one chain where it drew "
        "two. Reach: of the 13 structures at hand where the fix joins a new pair, 4 are symmetric and reach it: "
        "26EL and 7SPI (both measured above), 7SPB (7SPI's sibling: the same TraV pair) and 9GTP (C6; its '1l' and "
        "'2l' are copies of one baseplate protein: 8 chains -> 7, '1 of 8' -> '1 of 9'); 9TX7 (tests/data, C3) "
        "changes _one_of_each but its focus auto stays whole (6 chains).",
        "notes": "The third sample's committed report (7696947, views all) against this fix's run: 83 -> 82 pass, 38 "
        "entries differ; 37 of them differ identically against the BASE run (cdd315f), so they come from the code "
        "since 7696947, not from this fix; the 38th is 7SPI (area only). Read-only: measures/A-third-committed.json.",
    },
    "fix-letters": {
        "title": "Plan 29 Task 8b: a base letter sits nearest its own strand (22MZ, ruling P29-14)",
        "mechanism": "na_labels._place_letters took a single-stranded nucleotide's first spot (beside its strand, "
        "further out, or a little along it) clear of text and lines, whichever strand it was nearest. 22MZ's "
        "pseudoknot routes runs 2 and 6 side by side 0.65 NT_ROW apart: U57's spot beside run 6 lay 0.35 NT_ROW from "
        "run 2's line, on run 2's far side, and 7 of its 26 unpaired letters were nearer another strand than their "
        "own. A spot beside or along the strand is now clear only where no other strand's drawn line is nearer it "
        "than its own strand's (na_labels._own_nearest; a tie is its own). The later passes (by the drawn outlines, "
        "the rings) are unchanged. No new constant.",
        "named": ["22MZ"],
        "expected": "22MZ's letters (7 nearer another strand -> 3); elsewhere only letters whose first clear spot was "
        "nearer another strand",
        "explained": {
            **dict.fromkeys(
                (
                    "1EHZ 1GID 1HQ1 1MSW 1QTQ 1RLG 1Y26 1Z43 2GIS 2NZ4 2OEU 30QM 3D0U 3EGZ 3IGI 4AS1 4LVV 4OO8 4P9R "
                    "4PJO 4YBB 6QIK 7F5S 7XSP 8XJV 9K4F 9OWM"
                ).split(),
                "base letters moved off spots nearer another strand (the SVG only); no key moved",
            ),
            "22MZ": "the named entry: letters nearer another strand 7 -> 3 (0:50 G, 0:56 C, 0:67 G stay: no spot "
            "beside or along the strand both clear and nearest their own); no criterion moved",
            "9YPK": "letters nearer another strand 408 -> 130 of 1196 (262 moved; 2 newly nearer another strand); "
            "moved letters' median distance to their nucleotide 1.0 -> 1.24 NT_ROW; no criterion moved",
            "9YPN": "letters nearer another strand 347 -> 103 of 1144 (219 moved; 2 newly nearer another strand); no "
            "criterion moved",
            "3G6E": "better: A6 now passes (FAIL -> PASS)",
            "6IZP": "area 85.29 -> 77.33 (moved); letters moved",
            "9IWF": "area 232.25 -> 229.37 (moved); letters moved",
            "5F9R": "better: na_protein_crossings 7 -> 6",
            "2GDI": "WORSE: A6 now fails, a letter on a tie (na:nt:X:45 x na:tie:20-45): the letter moved to a spot "
            "nearest its own strand, where the tie's label is placed after it",
            "9WB3": "WORSE: A6 now fails, a letter on a tie mark (na:nt:R:75 x na:tie:43-75:mark), as 2GDI",
            "1URN": "WORSE: na_protein_crossings 1 -> 2 (a protein loop routed across a moved letter's strand)",
        },
        "control": "",
        "notes": "Rejected and reverted (43e7637); the code is kept as reverted.patch beside this record. Measured at "
        "87952c1 against Fix A's runs (749648c's code), on zyla-office: the NA set (na_fixes_record's cases and "
        "controls plus the NA campaign's entries, 69, default view) 54 -> 52 pass, 2GDI and 9WB3 failing A6, 1URN "
        "one more NA-protein crossing; so the keys are not identical or better (ruling P29-14). Also against it: "
        "tests/test_snapshots.py's 1EHZ-trna (na:nt:A:15 moves about 40 px) and 30QM-ncrna (7 nucleotides move) "
        "fail, and the moves do not read clearly better; _place_letters on the whole 4YBB SSU RNA takes 1.30 -> "
        "2.40 s on one host (gpu-desktop, 6 runs each), over test_na_layout's 1.5 s limit there (_own_nearest "
        "measures every strand segment for every spot; a cell hash would bound it, but not the A6 failures). The "
        "outline-pass variant (3 -> 0 on 22MZ, letters 3.8 rows out) was not taken. 22MZ's letters stay as they "
        "were (7 of 26 nearer another strand), for the author.",
    },
}


def _stamp() -> dict:
    if str(ROOT / "src") not in sys.path and (ROOT / "src" / "foldmap").is_dir():
        sys.path.insert(0, str(ROOT / "src"))
    from foldmap.source_stamp import source_stamp

    return source_stamp(__file__)


def _halves() -> tuple[list[str], set[str]]:
    split = json.loads(SAMPLE.read_text(encoding="utf-8"))["split"]
    return list(split["tune"]), set(split["check"])


def _module(path: Path, name: str):
    """A foldmap module as it was (its file), imported inside the foldmap package so its relative imports resolve."""
    import foldmap  # noqa: F401

    spec = importlib.util.spec_from_file_location(f"foldmap.{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _script(name: str):
    """A script of scripts/ imported as a module (random100, sensitivity, prune, na_fixes_record)."""
    if str(ROOT / "scripts") not in sys.path:
        sys.path.insert(0, str(ROOT / "scripts"))
    return importlib.import_module(name)


def structures() -> dict[str, str]:
    """path -> where it comes from: tests/data, the tune half (fetched), the download cache (less the check half)."""
    from foldmap.fetch import fetch

    tune, check = _halves()
    paths: dict[str, str] = {}
    for p in sorted((ROOT / "tests" / "data").rglob("*")):
        if p.is_file() and p.name.endswith((".cif", ".cif.gz")):
            paths[str(p)] = "tests/data"
    for i in tune:
        paths.setdefault(str(fetch(i)), "tune")
    for p in sorted((Path.home() / ".cache" / "foldmap").glob("*.cif")):
        if p.stem.upper() not in check:
            paths.setdefault(str(p), "cache")
    return paths


# ------------------------------------------------------------------------------------------------ helical


def _same_helical(a, b) -> bool:
    import numpy as np

    if a is None or b is None:
        return a is None and b is None
    (sa, ta, axa, ca, ra, cha, rma), (sb, tb, axb, cb, rb, chb, rmb) = a, b
    return (sa, ta, ra, cha, rma) == (sb, tb, rb, chb, rmb) and np.array_equal(axa, axb) and np.array_equal(ca, cb)


def _counting(asm) -> list[int]:
    calls, real = [0], asm.pair

    def pair(a, b):
        calls[0] += 1
        return real(a, b)

    asm.pair = pair
    return calls


def helical_one(path: str, old_file: str, max_copies: int) -> dict:
    from foldmap import symmetry as new
    from foldmap.io import load_backbone, read_model

    old = _module(Path(old_file), "_symmetry_before")
    name = Path(path).name.split(".")[0]
    try:
        bb = load_backbone(path, "auto", model=read_model(path, "auto"))
    except Exception as e:  # not a structure foldmap reads: said, not measured
        return {"name": name, "skipped": f"{type(e).__name__}: {e}"[:200]}
    asm_new, asm_old = new._Assembly(bb), new._Assembly(bb)
    largest = max((len(c) for c in asm_new.clusters), default=0)
    out: dict = {"name": name, "chains": len(asm_new.ca), "largest_cluster": largest}
    if largest > max_copies:
        out["skipped"] = f"largest cluster {largest} > {max_copies}"
        return out
    for mc in (new.MIN_HELIX_COPIES, new.MIN_HELIX_COPIES_ASKED):
        n_calls, o_calls = _counting(asm_new), _counting(asm_old)
        t0 = time.perf_counter()
        b_new = new._helical(asm_new, new.TOLERANCE, mc)
        t1 = time.perf_counter()
        b_old = old._helical(asm_old, new.TOLERANCE, mc)
        t2 = time.perf_counter()
        out[f"min{mc}"] = {
            "same": _same_helical(b_new, b_old),
            "copies": len(b_new[5][0]) if b_new is not None else None,
            "pair_calls_before": o_calls[0],
            "pair_calls_after": n_calls[0],
            "seconds_before": round(t2 - t1, 3),
            "seconds_after": round(t1 - t0, 3),
        }
        del asm_new.pair, asm_old.pair  # the class's own pair again for the next count
    return out


def cmd_helical(args) -> int:
    paths = structures()
    rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = [(pool.submit(helical_one, p, args.old, args.max_copies), src) for p, src in paths.items()]
        for f, src in futs:
            rows.append({**f.result(), "source": src})
    runs = [(r, r[k]) for r in rows for k in ("min4", "min3") if k in r]
    summary = {
        "structures": len(rows),
        "measured": sum(1 for r in rows if "min4" in r),
        "skipped": sorted(f"{r['name']}: {r['skipped']}" for r in rows if "skipped" in r),
        "differ": sorted({r["name"] for r, m in runs if not m["same"]}),
        "helical_found": {k: sum(1 for r in rows if r.get(k, {}).get("copies")) for k in ("min4", "min3")},
        "pair_calls": [sum(m["pair_calls_before"] for _, m in runs), sum(m["pair_calls_after"] for _, m in runs)],
        "seconds": [
            round(sum(m["seconds_before"] for _, m in runs), 2),
            round(sum(m["seconds_after"] for _, m in runs), 2),
        ],
        "more_pair_calls": sorted({r["name"] for r, m in runs if m["pair_calls_after"] > m["pair_calls_before"]}),
        "slowest_by": max(
            ((round(m["seconds_after"] - m["seconds_before"], 3), r["name"]) for r, m in runs), default=None
        ),
    }
    doc = {"kind": "helical", "source_stamp": _stamp(), "summary": summary, "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 1 if summary["differ"] else 0


def _helical_lines(m: dict, named: list[str]) -> list[str]:
    s = m["summary"]
    lines = [
        f"- {s['measured']} structures measured of {s['structures']} (tests/data, the tune half, the download cache "
        f"less the check half), {len(s['skipped'])} skipped: the return of _helical before and after differs on "
        f"{len(s['differ'])}{': ' + ', '.join(s['differ']) if s['differ'] else ''} (MIN_HELIX_COPIES and "
        "MIN_HELIX_COPIES_ASKED; every float compared bit for bit)",
        f"- a screw found on {s['helical_found']['min4']} (4 copies or more) and {s['helical_found']['min3']} (3 or "
        "more)",
        f"- _Assembly.pair calls, all structures: {s['pair_calls'][0]:,} before, {s['pair_calls'][1]:,} after; "
        f"_helical's seconds: {s['seconds'][0]} before, {s['seconds'][1]} after (run in parallel: indicative)",
        f"- more pair calls after on: {', '.join(s['more_pair_calls']) or 'none'} (a screw passes, so every copy's "
        "images are gathered at once)",
    ]
    rows = {r["name"]: r for r in m["rows"] if r["source"] in ("tune", "tests/data")}
    top = sorted((r for r in rows.values() if "min4" in r), key=lambda r: r["min4"]["pair_calls_before"], reverse=True)[
        :6
    ]
    lines += [
        "",
        "| entry | where | largest cluster | screw (copies) | pair calls | _helical s |",
        "|---|---|---|---|---|---|",
    ]
    for r in sorted(top, key=lambda r: r["name"] not in named):
        x = r["min4"]
        lines.append(
            f"| {r['name']} | {r['source']} | {r['largest_cluster']} | {x['copies'] or '-'} | "
            f"{x['pair_calls_before']:,} -> {x['pair_calls_after']:,} | {x['seconds_before']} -> {x['seconds_after']} |"
        )
    return lines


# ------------------------------------------------------------------------------------------------ science


def science_one(path: str, old_file: str) -> dict:
    from foldmap import figcheck
    from foldmap.spec import draw_spec, figure_spec

    old = _module(Path(old_file), "_figcheck_before")
    name = Path(path).name.split(".")[0]
    try:
        drawn = draw_spec(figure_spec(path))
    except Exception as e:  # not drawn: said, not judged
        return {"name": name, "error": f"{type(e).__name__}: {e}"[:200]}
    before, after = old.science_checks(drawn), figcheck.science_checks(drawn)
    return {"name": name, "before": before, "after": after, "changed": before != after}


def cmd_science(args) -> int:
    tune, _ = _halves()
    from foldmap.fetch import fetch

    paths = {str(fetch(i)): "tune" for i in tune}
    for p in sorted((ROOT / "tests" / "data").rglob("*")):
        if p.is_file() and p.name.endswith((".cif", ".cif.gz")):
            paths.setdefault(str(p), "tests/data")
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = [(pool.submit(science_one, p, args.old), src) for p, src in paths.items()]
        rows = [{**f.result(), "source": src} for f, src in futs]
    changed = [r for r in rows if r.get("changed")]

    def new_findings(r) -> list[str]:
        return [
            f"{k}: {t}" for k, found in r["after"].items() for t in found or [] if t not in (r["before"].get(k) or [])
        ]

    summary = {
        "entries": len(rows),
        "judged": sum(1 for r in rows if "after" in r),
        "errors": sorted(f"{r['name']}: {r['error']}" for r in rows if "error" in r),
        "differ": sorted(r["name"] for r in changed),
        "new_findings": {r["name"]: new_findings(r) for r in changed if new_findings(r)},
        "changes": {
            r["name"]: {k: [r["before"][k], r["after"][k]] for k in r["after"] if r["before"].get(k) != r["after"][k]}
            for r in changed
        },
        "S3_before": {r["name"]: r["before"]["S3"] for r in rows if r.get("before", {}).get("S3")},
        "S3_after": {r["name"]: r["after"]["S3"] for r in rows if r.get("after", {}).get("S3")},
    }
    doc = {"kind": "science", "source_stamp": _stamp(), "summary": summary, "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0


def _science_lines(m: dict, fix: dict) -> list[str]:
    s = m["summary"]
    lines = [
        f"- {s['judged']} entries judged of {s['entries']} (the tune half and tests/data, drawn as foldmap check "
        f"draws), {len(s['errors'])} not drawn{': ' + '; '.join(s['errors']) if s['errors'] else ''}",
        f"- science findings changed on: {', '.join(s['differ']) or 'none'}; new findings: "
        f"{', '.join(s['new_findings']) or 'none'}",
    ]
    for name, ch in s["changes"].items():
        for k, (a, b) in ch.items():
            lines.append(f"- {name} {k}: {a or 'none'} -> {b or 'none'}")
    lines.append(f"- S3 findings left after: {', '.join(f'{n}: {v}' for n, v in s['S3_after'].items()) or 'none'}")
    if fix.get("control"):
        lines += ["", "### The control", "", fix["control"]]
    return lines


# ------------------------------------------------------------------------------------------------ campaign


def _na_changes(before: dict[str, dict], after: dict[str, dict]) -> dict:
    """na_fixes_record.compare over every entry of both runs: per entry the kept A17 measures that changed and the
    judged ones that got worse (the named fallbacks' total and an error among them), and those that got better."""
    na = _script("na_fixes_record")
    out = {}
    for i in sorted(set(before) | set(after)):
        a, b = before.get(i, {}), after.get(i, {})
        changed = {
            k: [a.get(k), b.get(k)] for k in na.KEPT if a.get(k) != b.get(k) and (a.get(k), b.get(k)) != (None, 0)
        }
        worse = [k for k in na.JUDGED if na.worse(k, a.get(k), b.get(k))]
        better = [k for k in na.JUDGED if na.worse(k, b.get(k), a.get(k))]
        fa, fb = sum((a.get("named_fallbacks") or {}).values()), sum((b.get("named_fallbacks") or {}).values())
        if fa != fb:
            changed["named_fallbacks"] = [fa, fb]
            (worse if fb > fa else better).append("named_fallbacks")
        if a.get("error") != b.get("error"):
            changed["error"] = [a.get("error"), b.get("error")]
            (worse if b.get("error") else better).append("error")
        if changed:
            out[i] = {"changed": changed, "worse": worse, "better": better}
    return out


def cmd_campaign(args) -> int:
    r100, sens, prune = _script("random100"), _script("sensitivity"), _script("prune")
    reps = [r100.read_report(Path(d) / "report.json") for d in (args.before, args.after)]
    keys = [sens.keys_of(r) for r in reps]
    verdict = prune.campaign_verdict(*keys)
    svg = sorted(
        i
        for i in set(keys[0]) | set(keys[1])
        if keys[0].get(i, {}).get("svg_sha256") != keys[1].get(i, {}).get("svg_sha256")
    )
    summary = {
        "label": args.label,
        "entries": [len(keys[0]), len(keys[1])],
        "passed": [sens.passed(keys[0]), sens.passed(keys[1])],
        "verdict": verdict["verdict"],
        "worse": verdict["worse"],
        "better": verdict["better"],
        "moved": verdict["moved"],
        "svg_moved": svg,
        "differ": sorted({line.split(":")[0] for k in ("worse", "better", "moved") for line in verdict[k]} | set(svg)),
        "stamps": [(r.get("source_stamp") or {}).get("git") for r in reps],
        "reference": bool(args.reference),
    }
    rows = {"keys": {i: [keys[0].get(i), keys[1].get(i)] for i in summary["differ"]}}
    if args.na:
        na = _script("na_fixes_record")
        nr = [{e["id"]: na._row(e) for e in r["entries"]} for r in reps]
        ch = _na_changes(*nr)
        summary["na"] = {
            "entries": [len(nr[0]), len(nr[1])],
            "failing_a17": [sorted(i for i, r in x.items() if r.get("na_failed")) for x in nr],
            "changed": sorted(ch),
            "worse": {i: d["worse"] for i, d in ch.items() if d["worse"]},
            "better": {i: d["better"] for i, d in ch.items() if d["better"]},
        }
        rows["na"] = ch
    doc = {"kind": "campaign", "source_stamp": _stamp(), "summary": summary, "rows": rows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0


def cmd_navalidate(args) -> int:
    na = _script("na_fixes_record")
    reps = [json.loads((Path(d) / "report.json").read_text(encoding="utf-8")) for d in (args.before, args.after)]
    rows = [{e["entry"]: e for e in r["entries"]} for r in reps]
    ch = {}
    for i in sorted(set(rows[0]) | set(rows[1])):
        a, b = rows[0].get(i, {}), rows[1].get(i, {})
        changed = {k: [a.get(k), b.get(k)] for k in na.KEPT if a.get(k) != b.get(k)}
        worse = [k for k in na.JUDGED if na.worse(k, a.get(k), b.get(k))]
        better = [k for k in na.JUDGED if na.worse(k, b.get(k), a.get(k))]
        if (a.get("failed") or []) != (b.get("failed") or []):
            changed["failed"] = [a.get("failed"), b.get("failed")]
        if changed:
            ch[i] = {"changed": changed, "worse": worse, "better": better}
    summary = {
        "label": "the NA campaign (scripts/na_validate.py, docs/validation/na)",
        "entries": [len(rows[0]), len(rows[1])],
        "passed": [sum(not r.get("failed") for r in x.values()) for x in rows],
        "differ": sorted(ch),
        "worse": {i: d["worse"] for i, d in ch.items() if d["worse"]},
        "better": {i: d["better"] for i, d in ch.items() if d["better"]},
        "stamps": [(r.get("source_stamp") or {}).get("git") for r in reps],
    }
    doc = {"kind": "navalidate", "source_stamp": _stamp(), "summary": summary, "rows": ch}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0


def _campaign_lines(m: dict, fix: dict) -> list[str]:
    s, why = m["summary"], fix.get("explained", {})
    st = [(x or {}).get("commit", "?")[:7] for x in s["stamps"]]
    lines = [
        f"### {s['label']}",
        "",
        f"- before git {st[0]}, after git {st[1]}; entries {s['entries'][0]} -> {s['entries'][1]}; pass "
        f"{s['passed'][0]} -> {s['passed'][1]}; keys: {s['verdict']}",
        f"- worse: {'; '.join(s['worse']) or 'none'}",
        f"- better: {'; '.join(s['better']) or 'none'}",
        f"- moved (no judge orders it): {'; '.join(s['moved']) or 'none'}",
        f"- the SVG moved on: {', '.join(s['svg_moved']) or 'none'}",
    ]
    if "na" in s:
        n = s["na"]
        lines += [
            f"- A17 measures ({n['entries'][1]} entries): changed on {', '.join(n['changed']) or 'none'}; worse: "
            f"{'; '.join(f'{i} {v}' for i, v in n['worse'].items()) or 'none'}; better: "
            f"{'; '.join(f'{i} {v}' for i, v in n['better'].items()) or 'none'}",
            f"- failing A17 before: {', '.join(n['failing_a17'][0]) or 'none'}; after: "
            f"{', '.join(n['failing_a17'][1]) or 'none'}",
        ]
        for i, d in sorted(m["rows"].get("na", {}).items()):
            lines.append(f"  - {i}: " + "; ".join(f"{k} {a} -> {b}" for k, (a, b) in sorted(d["changed"].items())))
    if s.get("reference"):
        lines.insert(2, "For reference, not judged: the before is older code, so its differences are not this fix's.")
        lines.insert(3, "")
        return lines + [""]
    for i in s["differ"]:
        lines.append(f"- **{i}**: {why.get(i, 'NOT EXPLAINED')}")
    return lines + [""]


def _navalidate_lines(m: dict, fix: dict) -> list[str]:
    s, why = m["summary"], fix.get("explained", {})
    st = [(x or {}).get("commit", "?")[:7] for x in s["stamps"]]
    lines = [
        f"### {s['label']}",
        "",
        f"- before git {st[0]}, after git {st[1]}; entries {s['entries'][0]} -> {s['entries'][1]}; pass "
        f"{s['passed'][0]} -> {s['passed'][1]}",
        f"- measures changed on: {', '.join(s['differ']) or 'none'}; worse: "
        f"{'; '.join(f'{i} {v}' for i, v in s['worse'].items()) or 'none'}; better: "
        f"{'; '.join(f'{i} {v}' for i, v in s['better'].items()) or 'none'}",
    ]
    for i, d in sorted(m["rows"].items()):
        lines.append(f"  - {i}: " + "; ".join(f"{k} {a} -> {b}" for k, (a, b) in sorted(d["changed"].items())))
        lines.append(f"    {why.get(i, 'NOT EXPLAINED')}")
    return lines + [""]


def _measure_ok(m: dict, explained: set) -> bool:
    s = m["summary"]
    if m["kind"] == "campaign" and s.get("reference"):  # older code against the fix: said, not judged
        return True
    if m["kind"] == "campaign":
        return not s["worse"] and not (s.get("na") or {}).get("worse") and set(s["differ"]) <= explained
    if m["kind"] == "navalidate":
        return not s["worse"] and set(s["differ"]) <= explained
    return False


# ------------------------------------------------------------------------------------------------ report


def bench(path: Path) -> dict:
    """bench.txt: each entry's note ("identical" or what changed) and the exit."""
    notes, code = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("bench exit "):
            code = int(line.split()[-1])
            continue
        parts = line.split()
        if parts and ("identical" in line or "changed" in line or "MISSING" in line):
            notes[parts[0]] = "identical" if line.rstrip().endswith("identical") else line.split("  ")[-1].strip()
    return {
        "exit": code,
        "entries": len(notes),
        "identical": sum(v == "identical" for v in notes.values()),
        "notes": notes,
    }


def tune(before: dict, after: dict, named: list[str]) -> dict:
    a = {e["id"]: e for e in before["entries"] if e.get("half") == "tune"}
    b = {e["id"]: e for e in after["entries"] if e.get("half") == "tune"}
    common = sorted(a.keys() & b.keys())

    def moved(i, f):
        return a[i].get("files", {}).get(f) != b[i].get("files", {}).get(f)

    return {
        "entries": [len(a), len(b)],
        "missing": sorted(a.keys() ^ b.keys()),
        "svg_changed": [i for i in common if moved(i, "app.svg")],
        "png_changed": [i for i in common if moved(i, "app.png")],
        "criteria_changed": [
            i
            for i in common
            if (a[i]["criteria_failed"], a[i]["floor_ok"], a[i]["error"], a[i]["refused"])
            != (b[i]["criteria_failed"], b[i]["floor_ok"], b[i]["error"], b[i]["refused"])
        ],
        "criteria_worse": [
            i
            for i in common
            if set(b[i]["criteria_failed"] or []) - set(a[i]["criteria_failed"] or [])
            or (a[i]["floor_ok"] and not b[i]["floor_ok"])
            or (b[i]["error"] and not a[i]["error"])
            or (b[i]["refused"] and not a[i]["refused"])
        ],
        "named_seconds": {
            i: {"draw": [a[i]["draw_seconds"], b[i]["draw_seconds"]], "entry": [a[i]["seconds"], b[i]["seconds"]]}
            for i in named
            if i in a and i in b
        },
        "draw_seconds_total": [
            round(sum(a[i]["draw_seconds"] or 0 for i in common), 2),
            round(sum(b[i]["draw_seconds"] or 0 for i in common), 2),
        ],
        "stamps": [before.get("source_stamp", {}).get("git"), after.get("source_stamp", {}).get("git")],
    }


def cmd_report(args) -> int:
    fix = FIXES[args.name]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    measures = []
    for f in args.measure:
        m = json.loads(Path(f).read_text(encoding="utf-8"))
        name = Path(f).name if len(args.measure) > 1 else f"{m['kind']}.json"
        if Path(f).resolve() != (out / name).resolve():
            shutil.copyfile(f, out / name)
        measures.append((name, m))
    load = lambda d: json.loads((Path(d) / "report.json").read_text(encoding="utf-8"))  # noqa: E731
    t = tune(load(args.before), load(args.after), fix["named"])
    b = bench(Path(args.bench))
    if "explained" in fix:  # a fix that may move figures (Task 8b on): nothing worse, every move explained
        why = set(fix["explained"])
        moved = set(t["svg_changed"]) | set(t["png_changed"])
        bench_moved = {n for n, v in b["notes"].items() if v != "identical"}
        ok = (
            not t["missing"]
            and not t["criteria_worse"]
            and moved <= why
            and b["exit"] is not None
            and (b["exit"] == 0 or bench_moved)
            and bench_moved <= why
            and all(_measure_ok(m, why) for _, m in measures)
        )
    else:
        identical = not (t["svg_changed"] or t["png_changed"] or t["missing"])
        bench_ok = b["exit"] == 0 and b["identical"] == b["entries"]
        may = set(fix.get("may_change", []))
        ok = identical and bench_ok
        for _, m in measures:
            ok = ok and set(m["summary"].get("differ", [])) <= may and not m["summary"].get("new_findings")
    verdict = "kept" if ok else "rejected"
    rep = {
        "name": args.name,
        "source_stamp": _stamp(),
        "fix": fix,
        "before": str(args.before),
        "after": str(args.after),
        "tune": t,
        "bench": b,
        "measures": [{"file": n, "kind": m["kind"], "summary": m["summary"]} for n, m in measures],
        "verdict": verdict,
    }
    if len(measures) == 1:  # the short records of Tasks 6 and 7 as they were written
        rep["measure"] = rep.pop("measures")[0]
    (out / "report.json").write_text(json.dumps(rep, indent=1) + "\n", encoding="utf-8")
    (out / "report.md").write_text(render(rep, measures), encoding="utf-8")
    print(f"{args.name}: {verdict}")
    return 0


def render(rep: dict, measures: list[tuple[str, dict]]) -> str:
    t, b, fix = rep["tune"], rep["bench"], rep["fix"]
    st = rep["source_stamp"].get("git", {})
    why = fix.get("explained", {})
    scope = (
        "the tune half and the bench, the check half, the third sample and the NA controls left to Task N+1"
        if "explained" not in fix
        else "the tune half, the bench and the measures ruling P29-14 names; the check half left to Task N+1"
    )
    lines = [
        f"# {fix['title']}",
        "",
        f"Written by scripts/app29_light.py report at {st.get('commit', '?')[:7]} "
        f"({'dirty' if st.get('dirty') else 'clean'}); ruling P29-13's light loop: {scope}.",
        "",
        f"**Verdict: {rep['verdict']}.** Expected to change: {fix['expected']}.",
        "",
        "## Mechanism",
        "",
        fix["mechanism"],
        "",
        "## The tune half (scripts/app29.py run, the app's way)",
        "",
        f"- before: `{rep['before']}` (git {(t['stamps'][0] or {}).get('commit', '?')[:7]}); after: `{rep['after']}` "
        f"(git {(t['stamps'][1] or {}).get('commit', '?')[:7]})",
        f"- entries {t['entries'][0]} -> {t['entries'][1]}; missing: {', '.join(t['missing']) or 'none'}",
        f"- the app's SVG changed on: {', '.join(t['svg_changed']) or 'none'}; its PNG on: "
        f"{', '.join(t['png_changed']) or 'none'}",
        f"- criteria, floor, error or refusal changed on: {', '.join(t['criteria_changed']) or 'none'}; worse on: "
        f"{', '.join(t.get('criteria_worse', [])) or 'none'}",
        f"- draw seconds, all tune entries: {t['draw_seconds_total'][0]} -> {t['draw_seconds_total'][1]} (8 entries "
        "drawn at once: indicative)",
    ]
    for i, s in t["named_seconds"].items():
        lines.append(
            f"- {i}: draw {s['draw'][0]} -> {s['draw'][1]} s; the entry's process (draw, PNG, checks) "
            f"{s['entry'][0]} -> {s['entry'][1]} s"
        )
    if why:
        for i in sorted(set(t["svg_changed"]) | set(t["png_changed"])):
            lines.append(f"- **{i}**: {why.get(i, 'NOT EXPLAINED')}")
    lines += [
        "",
        "## bench --check --exact",
        "",
        f"- {b['identical']} of {b['entries']} identical, exit {b['exit']}",
    ]
    lines += [
        f"- {n}: {v}{' -- ' + why[n] if n in why else ''}" for n, v in sorted(b["notes"].items()) if v != "identical"
    ]
    for name, m in measures:
        lines += ["", f"## The fix's own measure ({name})", ""]
        if m["kind"] == "helical":
            lines += _helical_lines(m, fix["named"])
        elif m["kind"] == "science":
            lines += _science_lines(m, fix)
        elif m["kind"] == "campaign":
            lines += _campaign_lines(m, fix)
        elif m["kind"] == "navalidate":
            lines += _navalidate_lines(m, fix)
    for k in ("control", "notes"):
        if fix.get(k) and "explained" in fix:
            lines += ["", f"## {'The control' if k == 'control' else 'Notes'}", "", fix[k]]
    lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("helical", help="symmetry._helical before and after on every structure at hand")
    h.add_argument("--old", required=True, help="symmetry.py as it was")
    h.add_argument("--out", required=True)
    h.add_argument("--jobs", type=int, default=24)
    h.add_argument("--max-copies", type=int, default=64)
    c = sub.add_parser("science", help="figcheck.science_checks before and after on the tune half and tests/data")
    c.add_argument("--old", required=True, help="figcheck.py as it was")
    c.add_argument("--out", required=True)
    c.add_argument("--jobs", type=int, default=24)
    g = sub.add_parser("campaign", help="two random100 runs' keys (and A17 measures) before and after")
    g.add_argument("--before", required=True)
    g.add_argument("--after", required=True)
    g.add_argument("--label", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--na", action="store_true", help="also each entry's A17 measures (na_fixes_record._row)")
    g.add_argument("--reference", action="store_true", help="for reference only: said, never judged (older code)")
    v = sub.add_parser("navalidate", help="two na_validate reports before and after")
    v.add_argument("--before", required=True)
    v.add_argument("--after", required=True)
    v.add_argument("--out", required=True)
    r = sub.add_parser("report", help="the short record")
    r.add_argument("--name", required=True, choices=sorted(FIXES))
    r.add_argument("--before", required=True)
    r.add_argument("--after", required=True)
    r.add_argument("--bench", required=True)
    r.add_argument("--measure", required=True, nargs="+")
    r.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    cmds = {"helical": cmd_helical, "science": cmd_science, "campaign": cmd_campaign, "navalidate": cmd_navalidate}
    return {**cmds, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
