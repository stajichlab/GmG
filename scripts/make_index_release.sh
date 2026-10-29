#!/bin/bash
# usage: make_index_release.sh INDEX.gmgi [INDEX.gmgi ...]   (run in the directory that will be uploaded)
# Writes SHA256SUMS and MANIFEST.txt (the `gmg info` of each index) for an index release (docs/index-distribution.md).
set -euo pipefail
GMG=${GMG:-gmg}
: > SHA256SUMS; : > MANIFEST.txt
for f in "$@"; do
  if command -v sha256sum >/dev/null; then sha256sum "$f" >> SHA256SUMS; else shasum -a 256 "$f" >> SHA256SUMS; fi
  { echo "## $(basename "$f")"; $GMG info "$f"; echo; } >> MANIFEST.txt
done
cat SHA256SUMS
