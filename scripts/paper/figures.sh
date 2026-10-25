#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Every manuscript figure (make_figures.py), into build/figures/ of the manuscript folder.
#
#   scripts/paper/figures.sh            draw every foldmap panel, assemble F1-F5 and the supplementary gallery, check
#   scripts/paper/figures.sh f2         one target (draw [IDS] | gallery | f1 | f2 | f3 | f4 | f5 | fs | check) after
#                                       a draw; several figure targets may follow each other (figures.sh f2 fs)
#   PAPER_ZL=/path/to/zl scripts/paper/figures.sh ...   on a runner with zl's interface (see run_all.sh); exit 75
#                                       (busy) or 23 (rsync raced another writer) is retried
#
# Inputs: the manuscript folder's figures/data/ (the CI determinism artifacts F3d reads), and what F5 reads
# (build/domains.json, build/tables/, data/m2/dssp/, data/m2/na/). Needs foldmap on PATH (the paper environment),
# tectonic ($TECTONIC or on PATH), pdftocairo and pdftoppm.
# Provenance: a runner syncs the working tree, not a commit. Before each run this records HEAD and whether src/ or
# scripts/ differ from it (build/figures/provenance-local.json; make_figures.py copies it into figures.json beside
# foldmap.source_stamp of the code that drew) and warns when they do.
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$H/../.." && pwd)
M=${FOLDMAP_MANUSCRIPT:-$REPO/docs/superpowers/work/manuscript}
out=$M/build/figures
[[ $# -gt 0 ]] || set -- all
mkdir -p "$out"
for _ in $(seq 1 20); do
  head=$(git -C "$REPO" rev-parse HEAD)
  changed=$(git -C "$REPO" status --porcelain -- src scripts)
  dirty=false
  if [[ -n $changed ]]; then
    dirty=true
    echo "WARNING: src/ or scripts/ differ from HEAD $head; the figures are drawn from the working tree as it is:" >&2
    echo "$changed" | sed 's/^/  /' >&2
  fi
  {
    printf '{\n "head": "%s",\n "dirty": %s,\n "changed": [' "$head" "$dirty"
    sep=""
    while IFS= read -r l; do [[ -n $l ]] && printf '%s"%s"' "$sep" "$l" && sep=", "; done <<< "$changed"
    printf '],\n "recorded": "%s"\n}\n' "$(date -Is)"
  } > "$out/provenance-local.json"
  if [[ -n ${PAPER_ZL:-} ]]; then
    "$PAPER_ZL" --push "$M/figures/data" --push "$out" --push "$M/build/domains.json" --push "$M/build/tables" \
      --push "$M/data/m2/dssp" --push "$M/data/m2/na" \
      --pull "$out" -- env FOLDMAP_MANUSCRIPT="$M" ${TECTONIC:+TECTONIC="$TECTONIC"} python "$H/make_figures.py" "$@"
    rc=$?
  else
    (cd "$REPO" && python "$H/make_figures.py" "$@")
    rc=$?
  fi
  [[ $rc -ne 75 && $rc -ne 23 ]] && exit $rc
  echo "runner busy or raced (exit $rc); retrying in 60 s" >&2
  sleep 60
done
exit 1
