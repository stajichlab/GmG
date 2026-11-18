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
#   SPECIES_LOG  the species check (PopSeq#52) runs only when this is set: the file to which the verdict row of each run is appended
#                (header first, once). SPECIES_MATCH: `refuse` (default) refuses a run whose verdict is excluded by default, `warn`
#                keeps it and warns, `off` skips the check. SPECIES_INDEX: the reference index (use the small bfd_s10000.gmgi: 0.6 to
#                1.5 s and 1.4 GB per run, against 6 to 8 s and 13.7 GB for bfd_s1000). SPECIES_FILED: a table `id<TAB>filed species
#                name` (scripts/popsra_filed_table.py makes it from the manifest). A check that cannot run (no readable index, no
#                filed name, `popsketch species-match` failing) keeps the run, records an `unchecked` row with the reason and warns.
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
# Species check: a verdict that is excluded by default (mismatch, host_dominated, no_fungal_call_suspect; the `exclude_default`
# column of `popsketch species-match`) makes the hook exit 66 under SPECIES_MATCH=refuse. The sketch is MOVED to $SKETCH_DIR/refused/
# (not deleted): it stays for a re-check against a newer index, and `popsketch run --from-sketches` does not read the subfolder, so
# the folder still holds only samples that passed. The verdict row is in SPECIES_LOG either way.
# A non-zero exit becomes exit 14 (`hook_failed`) in popseq fetch; the download itself is still `ok`.
set -euo pipefail
: "${SKETCH_DIR:?set SKETCH_DIR to the folder of sketches}"
run=${1:?usage: popsra_sketch_hook.sh RUN FASTQ_1 [FASTQ_2]}
shift
[ "$#" -ge 1 ] || { echo "no FASTQ files given for $run" >&2; exit 2; }
for f in "$@"; do [ -s "$f" ] || { echo "missing or empty: $f" >&2; exit 2; }; done

case "${SPECIES_MATCH:-refuse}" in refuse|warn|off) ;; *) echo "SPECIES_MATCH must be refuse, warn or off, got '${SPECIES_MATCH}'" >&2; exit 2 ;; esac

# One `unchecked` row (the check could not run), with the reason, to SPECIES_LOG and a warning to stderr.
species_unchecked() { # RUN FILED REASON
  ( set -C; printf "id\tfiled_name\tstatus\texclude_default\tdetected_genus\tdetected_species\tbest_shared_hashes\tbest_share\tpurity\treason\n" > "$SPECIES_LOG" ) 2>/dev/null || true
  printf "%s\t%s\tunchecked\tno\t\t\t\t\t\t%s\n" "$1" "$2" "$(printf "%s" "$3" | tr '\t\n' '  ')" >> "$SPECIES_LOG"
  echo "$1: WARNING species check not run: $3" >&2
}

# The species check of one run whose sketch is complete. Returns 66 when the run is refused, else 0.
species_check() { # RUN
  local run=$1 filed idx rc=0 tmp header line status excl
  idx=${SPECIES_INDEX:-}
  filed=$(awk -F'\t' -v r="$run" '$1 == r {print $2; exit}' "${SPECIES_FILED:-/dev/null}" 2>/dev/null || true)
  if [ -z "$idx" ] || [ ! -r "$idx" ]; then species_unchecked "$run" "$filed" "no readable index (SPECIES_INDEX='$idx')"; return 0; fi
  if [ -z "$filed" ]; then species_unchecked "$run" "" "no filed species name for the run (SPECIES_FILED='${SPECIES_FILED:-}')"; return 0; fi
  tmp=$(mktemp -d "${SCRATCH:-${TMPDIR:-/tmp}}/speciesmatch.XXXXXX")
  printf "id\tfiled\n%s\t%s\n" "$run" "$filed" > "$tmp/filed.tsv"
  "${POPSKETCH:-popsketch}" species-match --index "$idx" -s "$SKETCH_DIR" --filed "$tmp/filed.tsv" -o "$tmp/out.tsv" </dev/null >/dev/null 2>"$tmp/err" || rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$tmp/out.tsv" ]; then
    species_unchecked "$run" "$filed" "popsketch species-match failed (exit $rc): $(head -c 200 "$tmp/err")"
    rm -rf "$tmp"; return 0
  fi
  header=$(sed -n 1p "$tmp/out.tsv"); line=$(sed -n 2p "$tmp/out.tsv")
  ( set -C; printf "%s\n" "$header" > "$SPECIES_LOG" ) 2>/dev/null || true
  printf "%s\n" "$line" >> "$SPECIES_LOG"
  ( set -C; cp "$tmp/out.tsv.meta.json" "${SPECIES_LOG%.tsv}.meta.json" ) 2>/dev/null || true   # the index and cutoffs, once per task
  status=$(printf "%s" "$line" | cut -f3); excl=$(printf "%s" "$line" | cut -f4)
  rm -rf "$tmp"
  if [ "$excl" = yes ]; then
    if [ "${SPECIES_MATCH:-refuse}" = refuse ]; then
      mkdir -p "$SKETCH_DIR/refused"
      mv "$SKETCH_DIR/$run".json "$SKETCH_DIR/$run".gsk "$SKETCH_DIR/$run".half_a.gsk "$SKETCH_DIR/$run".half_b.gsk "$SKETCH_DIR/refused/"
      echo "$run: species check refused the run ($status): $(printf "%s" "$line" | cut -f10)" >&2
      return 66
    fi
    echo "$run: WARNING species check says $status; the run is kept (SPECIES_MATCH=${SPECIES_MATCH})" >&2
  fi
  return 0
}

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
if [ -n "${SPECIES_LOG:-}" ] && [ "${SPECIES_MATCH:-refuse}" != off ]; then
  species_check "$run" || exit $?
fi
if [ -n "${HOOK_LOG:-}" ]; then
  bytes=$(cat "$SKETCH_DIR/$run".json "$SKETCH_DIR/$run".gsk "$SKETCH_DIR/$run".half_a.gsk "$SKETCH_DIR/$run".half_b.gsk | wc -c)
  printf "%s\t%s\t%s\n" "$run" "$((SECONDS - start))" "$bytes" >> "$HOOK_LOG"
fi
