#!/bin/sh
# Compile the report only; never runs analyses. Tectonic lives in _local/toolchain (not in git).
set -eu
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(dirname "$(dirname "$DIR")")
TECTONIC="$ROOT/_local/toolchain/tectonic"
[ -x "$TECTONIC" ] || TECTONIC=$(command -v tectonic)
cd "$DIR"
"$TECTONIC" --keep-logs --outdir "$DIR/output" main.tex
echo "PDF: $DIR/output/main.pdf"
