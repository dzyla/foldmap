#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Regenerate the paper's numbers and tables (build/) from one commit: a release passes its tag.
#
#   scripts/paper/run_all.sh [COMMIT]               every step on this machine (COMMIT default 838c579)
#   M2_RERUN=1 scripts/paper/run_all.sh [COMMIT]    the M2 measurements first (membrane vs OPM, DSSP vs mkdssp, the NA
#                                                   rerun with pseudoknots, the determinism Linux leg), then the rest
#   PAPER_ZL=/path/to/zl scripts/paper/run_all.sh   the steps that compute or draw through a runner wrapper with zl's
#                                                   interface, `zl [--push P]... [--pull P]... -- CMD` (sync this
#                                                   checkout to the same path elsewhere, run CMD there, copy each --pull
#                                                   back); exit 75 = busy, retried
#
# The manuscript folder (its text, inputs and build/) is docs/superpowers/work/manuscript of this checkout (gitignored)
# or $FOLDMAP_MANUSCRIPT. `python` is the paper environment's (scripts/paper/requirements.lock, foldmap installed).
set -euo pipefail
COMMIT=${1:-838c579}
H=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$H/../.." && pwd)
M=${FOLDMAP_MANUSCRIPT:-$REPO/docs/superpowers/work/manuscript}
export FOLDMAP_MANUSCRIPT=$M
D=$M/data
RUN=${TMPDIR:-/tmp}/foldmap-paper

# A step that computes or draws: through $PAPER_ZL with its --push/--pull options, or here with them dropped.
heavy() {
  if [[ -n ${PAPER_ZL:-} ]]; then
    local rc
    for _ in 1 2 3 4 5 6 7 8 9 10; do
      set +e; "$PAPER_ZL" "$@"; rc=$?; set -e
      [[ $rc -ne 75 ]] && return $rc
      echo "runner busy; waiting 60 s" >&2; sleep 60
    done
    return 75
  fi
  while [[ $1 != -- ]]; do shift 2; done
  shift
  (cd "$REPO" && "$@")
}

python "$H/extract.py" --commit "$COMMIT"                                    # build/sources, tree, provenance.json
heavy --push "$M/build" --pull "$M/build/domains.json" -- env FOLDMAP_MANUSCRIPT="$M" python "$H/domains.py"
heavy --push "$M/build" --pull "$M/build" -- env FOLDMAP_MANUSCRIPT="$M" python "$H/runtime.py"
if [[ "${M2_RERUN:-0}" == 1 ]]; then  # the M2 measurements (raw outputs under data/m2/; make_numbers reads what is there)
  python "$H/opm_fetch.py"                                                   # data/opm (cached; network)
  heavy --push "$M/build" --push "$D/opm" --pull "$D/m2/membrane" -- \
    env FOLDMAP_MANUSCRIPT="$M" python "$H/membrane_opm.py"
  heavy --push "$M/build" --pull "$D/m2/dssp" -- env FOLDMAP_MANUSCRIPT="$M" python "$H/dssp_mkdssp.py"
  heavy --pull "$RUN/na" -- bash -c "rm -rf '$RUN/na' && python scripts/na_validate.py --out '$RUN/na' \
    --require-rnapolis; mkdir -p '$RUN/na'"
  mkdir -p "$D/m2/na" && cp "$RUN/na/report.json" "$RUN/na/report.md" "$D/m2/na/" \
    && git -C "$REPO" rev-parse HEAD > "$D/m2/na/COMMIT"
  # the Linux leg: determinism.py's three configurations, compared (exit 1 on a difference)
  heavy --pull "$RUN/determinism" -- bash -c "rm -rf '$RUN/determinism' && \
    python scripts/paper/determinism.py --configurations '$RUN/determinism'; rc=\$?; mkdir -p '$RUN/determinism'; exit \$rc" \
    || echo "determinism: the Linux leg's runs differ (listed above); make_numbers reports which" >&2
  mkdir -p "$D/m2/determinism" && cp "$RUN"/determinism/linux-*.json "$D/m2/determinism/" \
    && git -C "$REPO" rev-parse HEAD > "$D/m2/determinism/COMMIT"
fi
python "$H/make_numbers.py"                                                  # build/numbers.json, tables/
python "$H/table_s1.py" || echo "table_s1: a value is not in the code at $COMMIT (named above)" >&2  # Table S1
python "$H/fill.py"                                                          # build/*-filled.tex, keys.json
python "$H/crosscheck.py" || echo "crosscheck: disagreements listed in build/crosscheck.json"
