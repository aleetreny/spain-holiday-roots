#!/usr/bin/env bash
# Build the arXiv source bundle from the journal manuscript: no referee spacing or line numbers,
# figures beside the source, the supplementary PDF appended, and the .bbl included.
# Output: _local/arxiv/ (folder, compiled check PDF and arxiv_source.zip).
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
MS="$ROOT/journal/manuscript"
OUT="$ROOT/_local/arxiv"
SRC="$OUT/source"
TECTONIC=${TECTONIC:-$HOME/Proyectos/spain-holiday-roots/_local/toolchain/tectonic}
[ -x "$TECTONIC" ] || TECTONIC=$(command -v tectonic)

sh "$MS/build.sh" >/dev/null
rm -rf "$OUT" && mkdir -p "$SRC" "$OUT/check"
cp "$MS/sn-jnl.cls" "$MS/sn-mathphys-num.bst" "$MS/references.bib" "$SRC/"
cp "$ROOT"/journal/figures/out/fig{1,2,3,4}.pdf "$SRC/"
cp "$MS/output/supplementary.pdf" "$SRC/supplementary.pdf"

python3 - "$MS/main.tex" "$SRC/ms.tex" <<'PY'
import sys
s = open(sys.argv[1]).read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old
    s = s.replace(old, new)
rep("[pdflatex,sn-mathphys-num,referee]{sn-jnl}", "[pdflatex,sn-mathphys-num]{sn-jnl}")
rep("\\usepackage[mathlines]{lineno}\n\\renewcommand{\\linenumberfont}{\\normalfont\\tiny}\n", "\\usepackage{pdfpages}\n")
rep("\\maketitle\n\\linenumbers\n", "\\maketitle\n")
rep("\\graphicspath{{../figures/out/}}", "\\graphicspath{{./}}")
rep("\\end{document}", "\\clearpage\n\\includepdf[pages=-]{supplementary.pdf}\n\\end{document}")
open(sys.argv[2], "w").write(s)
PY

"$TECTONIC" --keep-intermediates --outdir "$OUT/check" "$SRC/ms.tex" >"$OUT/check/tectonic.txt" 2>&1
cp "$OUT/check/ms.bbl" "$SRC/ms.bbl"
rm "$SRC/references.bib" "$SRC/sn-mathphys-num.bst"
cp "$OUT/check/ms.pdf" "$OUT/arxiv_preview.pdf"
(cd "$SRC" && zip -qX "$OUT/arxiv_source.zip" *)
echo "Source: $OUT/arxiv_source.zip ($(du -h "$OUT/arxiv_source.zip" | cut -f1))"
echo "Preview: $OUT/arxiv_preview.pdf"
