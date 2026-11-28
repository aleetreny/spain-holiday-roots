#!/usr/bin/env bash
# Build the Zenodo replication package from a committed state of this repository.
# Usage: journal/replication/build_zenodo.sh [git-ref] (default HEAD). Output: _local/zenodo/.
set -euo pipefail

ROOT=$(git rev-parse --show-toplevel)
REF=${1:-HEAD}
VERSION=${VERSION:-1.0.0}
NAME="spain-holiday-roots-replication-v${VERSION}"
OUT="$ROOT/_local/zenodo"
STAGE="$OUT/$NAME"

rm -rf "$STAGE" "$OUT/$NAME.zip"
mkdir -p "$STAGE"

git -C "$ROOT" archive "$REF" requirements.txt explore pipeline tests results \
  journal/analysis journal/figures journal/pilots/design | tar -x -C "$STAGE"

# Working notes and run logs are not part of the package.
rm -f "$STAGE/journal/figures/FIGURES_SPEC.md" "$STAGE/journal/figures/REVIEW_AESTHETICS.md" \
  "$STAGE/journal/figures/REPORT_FIGURES.md"
find "$STAGE/results" \( -name '*.log' -o -name '*.md' \) -delete

for f in README.md LICENSE.md CITATION.cff; do
  git -C "$ROOT" show "$REF:journal/replication/$f" > "$STAGE/$f"
done
printf 'repository https://github.com/aleetreny/spain-holiday-roots\nbranch journal\ncommit %s\nversion %s\n' \
  "$(git -C "$ROOT" rev-parse "$REF")" "$VERSION" > "$STAGE/VERSION.txt"

(cd "$OUT" && rm -f "$NAME.zip" && zip -qrX "$NAME.zip" "$NAME")
echo "Folder: $STAGE"
echo "Zip:    $OUT/$NAME.zip ($(du -h "$OUT/$NAME.zip" | cut -f1))"
echo "Files:  $(find "$STAGE" -type f | wc -l | tr -d ' ')"
