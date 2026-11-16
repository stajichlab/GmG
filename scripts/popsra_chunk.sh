#!/bin/bash
# Which manifest rows a chunk job takes. Sourced by popsra_fetch_chunks.slurm and tested by test_popsra_chunks.sh.
#
#   chunk_runs MANIFEST CHUNK I   prints the run_accession values of rows I*CHUNK .. I*CHUNK+CHUNK-1 (0-based, header
#                                 skipped), one per line; nothing if the chunk is past the end.
#   chunk_count MANIFEST CHUNK    prints the number of chunks (the array size is that minus 1).
# The manifest is a CSV or TSV with a run_accession column, as popseq fetch reads it.

# The columns are read with Python's csv module (popsra_manifest.py), not cut: a quoted field with a comma must not shift them.
_manifest_runs() { # MANIFEST: the run_accession of every data row
  python3 "$(dirname "${BASH_SOURCE[0]}")/popsra_manifest.py" "$1"
}

# The runs are read into a variable first so that a manifest error (no run_accession column) is returned, not lost
# on the left side of a pipe: a chunk job on a bad manifest must fail, not quietly process nothing.
chunk_runs() {
  local m=$1 chunk=$2 i=$3 runs
  runs=$(_manifest_runs "$m") || return $?
  [ -z "$runs" ] || printf "%s\n" "$runs" | sed -n "$((i * chunk + 1)),$((i * chunk + chunk))p"
}

chunk_count() {
  local m=$1 chunk=$2 runs
  runs=$(_manifest_runs "$m") || return $?
  [ -z "$runs" ] && { echo 0; return 0; }
  echo $(((($(printf "%s\n" "$runs" | wc -l) + chunk - 1)) / chunk))
}
