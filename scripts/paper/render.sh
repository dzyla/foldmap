#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Typeset build/main-filled.tex and build/supplement-filled.tex to PDF with Tectonic, into build/pdf/.
#
#   scripts/paper/render.sh                         here (run fill.py first so the filled files are current)
#   PAPER_ZL=/path/to/zl scripts/paper/render.sh    on a runner with zl's interface (see run_all.sh)
#
# Tectonic: $TECTONIC, else `tectonic` on PATH, else 0.15.0 installed on first use into
# ${XDG_CACHE_HOME:-~/.cache}/foldmap-paper/tectonic. The figures come from build/figures/ (make_figures.py) and are
# placed where the text includes them (figures/out/F<n>.pdf). With any TeX installation instead:
# cd build/pdf-src && latexmk -pdf main-filled.tex (after this script made pdf-src).
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$H/../.." && pwd)
M=${FOLDMAP_MANUSCRIPT:-$REPO/docs/superpowers/work/manuscript}
src=$M/build/pdf-src out=$M/build/pdf
rm -rf "$src" && mkdir -p "$src" "$out"
cp "$M"/build/main-filled.tex "$M"/build/supplement-filled.tex "$M"/build/runtime.pdf "$M"/refs.bib "$src"/
cp -r "$M"/build/tables "$src"/tables
if compgen -G "$M/build/figures/F*.pdf" > /dev/null; then
  mkdir -p "$src"/figures/out && cp "$M"/build/figures/F*.pdf "$src"/figures/out/
fi
# shellcheck disable=SC2016  # expanded by the shell that typesets
typeset_cmd='
set -e
t=${TECTONIC:-$(command -v tectonic || true)}
if [[ -z $t ]]; then
  d=${XDG_CACHE_HOME:-$HOME/.cache}/foldmap-paper/tectonic
  t=$d/tectonic
  if [[ ! -x $t ]]; then
    mkdir -p "$d"
    curl -fsSL https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.15.0/tectonic-0.15.0-x86_64-unknown-linux-musl.tar.gz | tar -xz -C "$d"
  fi
fi
mkdir -p "$OUT"
cd "$SRC"
for f in main-filled supplement-filled; do
  "$t" -X compile --keep-logs --outdir "$OUT" $f.tex 2>&1 | grep -v "^note: downloading" | tail -15
done
ls -la "$OUT"
'
if [[ -n ${PAPER_ZL:-} ]]; then
  "$PAPER_ZL" --push "$src" --pull "$out" -- env SRC="$src" OUT="$out" ${TECTONIC:+TECTONIC="$TECTONIC"} bash -c \
    "$typeset_cmd"
else
  SRC=$src OUT=$out bash -c "$typeset_cmd"
fi
