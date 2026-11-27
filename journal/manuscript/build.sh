#!/bin/sh
# Build the journal manuscript (and SI). Tectonic (XeTeX) for drafts; the source stays pdfLaTeX-compatible for submission.
set -eu
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TECTONIC=/Users/alejandrotreny/Proyectos/spain-holiday-roots/_local/toolchain/tectonic
[ -x "$TECTONIC" ] || TECTONIC=$(command -v tectonic)
cd "$DIR"
for f in main supplementary; do
  [ -f "$f.tex" ] && "$TECTONIC" --keep-logs --outdir "$DIR/output" "$f.tex"
done
echo "PDFs in $DIR/output"
