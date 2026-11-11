#!/bin/bash
# Which manifest rows a chunk job takes. Sourced by popsra_fetch_chunks.slurm and tested by test_popsra_chunks.sh.
#
#   chunk_runs MANIFEST CHUNK I   prints the run_accession values of rows I*CHUNK .. I*CHUNK+CHUNK-1 (0-based, header
#                                 skipped), one per line; nothing if the chunk is past the end.
#   chunk_count MANIFEST CHUNK    prints the number of chunks (the array size is that minus 1).
# The manifest is a CSV or TSV with a run_accession column, as popseq fetch reads it.

_manifest_runs() { # MANIFEST: the run_accession of every data row
  local m=$1 delim=","
  head -1 "$m" | grep -q $'\t' && delim=$'\t'
  local col
  col=$(head -1 "$m" | tr "$delim" '\n' | grep -n -x run_accession | cut -d: -f1)
  [ -n "$col" ] || { echo "no run_accession column in $m" >&2; return 2; }
  # grep exits 1 when there is no data row; an empty manifest is zero chunks, not an error
  tail -n +2 "$m" | cut -d"$delim" -f"$col" | { grep -v '^[[:space:]]*$' || true; }
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
