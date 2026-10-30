#!/bin/bash
# usage: fetch_examples.sh OUT_DIR [ID ...]   (IDs from examples/examples.tsv; default: all)
# Downloads the example genomes from NCBI so `gmg run OUT_DIR --index ...` can be tried without your own data.
set -euo pipefail
out=${1:?output directory}; shift || true
tsv=$(cd "$(dirname "$0")/.." && pwd)/examples/examples.tsv
mkdir -p "$out"
tail -n +2 "$tsv" | while IFS=$'\t' read -r id acc name expected note url; do
  if [ $# -gt 0 ]; then case " $* " in *" $id "*) ;; *) continue ;; esac; fi
  f="$out/${acc}.fna.gz"
  [ -s "$f" ] && { echo "have $acc ($name)"; continue; }
  echo "fetching $acc ($name)"
  curl -fsSL --retry 3 -o "$f.part" "$url" && mv "$f.part" "$f"
done
