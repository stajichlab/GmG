#!/bin/bash
# Sketch hook for `popseq fetch --exec` (PopSeq issue #13): turns one verified run into a sketch so the reads
# can be deleted. `popseq fetch` calls:  CMD RUN FASTQ_1 [FASTQ_2]   (one path for a single-end run).
#
#   popseq fetch SRR123 --tier coverage --coverage 25 --out "$SCRATCH/x" \
#       --exec scripts/popsra_sketch_hook.sh
#
# Environment:
#   SKETCH_DIR   required. The folder of sketches (shared by all samples of one set).
#   POPSKETCH    the popsketch binary (default: popsketch on PATH).
#   SKETCH_ARGS  extra arguments for `popsketch sketch` (for example "--scaled 5" in tests).
#   HOOK_LOG     optional. One line per run is appended: run, seconds, bytes of the sketch files.
#   CHECK_LOG    optional. The `popsketch check` line of each run is appended (header first).
#   POPSKETCH_CHECK  `off` skips the sampling check.
# `popseq fetch` also sets (when it knows them): POPSEQ_GENOME_MB (the genome size) and POPSEQ_COVERAGE (the x the kept
# reads should give if sampled evenly), with POPSEQ_SAMPLING, POPSEQ_TIER and POPSEQ_ROUTE for information.
#
# Sampling check: when POPSEQ_GENOME_MB is set, `popsketch check` compares the sketch with the genome. A sample that
# misses part of the genome or piles up on part of it (contiguous windows of a position-sorted run) fails with 75: the
# hook exits 75, which `popseq fetch` takes as "resample this run by stride". 65 means the check cannot judge the sampling
# (too little k-mer coverage, or not this genome): resampling would not help, so popseq fetch reports hook_failed. 2 is an
# error of the check itself. In every failing case the sketch files of the run are removed, so that only samples that
# passed the check are in the folder that `popsketch run --from-sketches` reads.
#
# --max-bases 0: `popseq fetch` has already limited the bases (for example to 25x of the genome), so the sketch
# must not cap them again (the default cap of 1 Gb would truncate a genome larger than 40 Mb).
# A non-zero exit becomes exit 14 (`hook_failed`) in popseq fetch; the download itself is still `ok`.
set -euo pipefail
: "${SKETCH_DIR:?set SKETCH_DIR to the folder of sketches}"
run=${1:?usage: popsra_sketch_hook.sh RUN FASTQ_1 [FASTQ_2]}
shift
[ "$#" -ge 1 ] || { echo "no FASTQ files given for $run" >&2; exit 2; }
for f in "$@"; do [ -s "$f" ] || { echo "missing or empty: $f" >&2; exit 2; }; done

start=$SECONDS
mkdir -p "$SKETCH_DIR"
# shellcheck disable=SC2086
"${POPSKETCH:-popsketch}" sketch --id "$run" -o "$SKETCH_DIR" --max-bases 0 ${SKETCH_ARGS:-} "$@" </dev/null

# The sketch must be complete: the whole sketch, both split halves and the metadata.
for part in "$run.json" "$run.gsk" "$run.half_a.gsk" "$run.half_b.gsk"; do
  [ -s "$SKETCH_DIR/$part" ] || { echo "popsketch left no $part for $run" >&2; exit 3; }
done
if [ -n "${POPSEQ_GENOME_MB:-}" ] && [ "${POPSKETCH_CHECK:-on}" != off ]; then
  rc=0
  # --header always: the header line is written to CHECK_LOG only when this hook creates the file (noclobber: exactly one
  # of several hooks that start together does), and the data line is appended whole
  # shellcheck disable=SC2086
  out=$("${POPSKETCH:-popsketch}" check --id "$run" -s "$SKETCH_DIR" --genome-mb "$POPSEQ_GENOME_MB" \
        ${POPSEQ_COVERAGE:+--coverage $POPSEQ_COVERAGE} --header </dev/null) || rc=$?
  header=$(printf "%s\n" "$out" | sed -n 1p); line=$(printf "%s\n" "$out" | sed -n 2p)
  if [ -n "${CHECK_LOG:-}" ] && [ -n "$line" ]; then
    ( set -C; printf "%s\n" "$header" > "$CHECK_LOG" ) 2>/dev/null || true
    printf "%s\n" "$line" >> "$CHECK_LOG"
  fi
  if [ "$rc" -ne 0 ]; then
    echo "$run: sampling check failed (exit $rc): $(printf "%s" "$line" | cut -f9)" >&2
    rm -f "$SKETCH_DIR/$run".json "$SKETCH_DIR/$run".gsk "$SKETCH_DIR/$run".half_a.gsk "$SKETCH_DIR/$run".half_b.gsk
    exit "$rc"
  fi
fi
if [ -n "${HOOK_LOG:-}" ]; then
  bytes=$(cat "$SKETCH_DIR/$run".json "$SKETCH_DIR/$run".gsk "$SKETCH_DIR/$run".half_a.gsk "$SKETCH_DIR/$run".half_b.gsk | wc -c)
  printf "%s\t%s\t%s\n" "$run" "$((SECONDS - start))" "$bytes" >> "$HOOK_LOG"
fi
