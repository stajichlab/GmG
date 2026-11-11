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
if [ -n "${HOOK_LOG:-}" ]; then
  bytes=$(cat "$SKETCH_DIR/$run".json "$SKETCH_DIR/$run".gsk "$SKETCH_DIR/$run".half_a.gsk "$SKETCH_DIR/$run".half_b.gsk | wc -c)
  printf "%s\t%s\t%s\n" "$run" "$((SECONDS - start))" "$bytes" >> "$HOOK_LOG"
fi
