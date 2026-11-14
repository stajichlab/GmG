#!/bin/bash
# Phase 1 (one species) and later phases of docs/plan_popsra_scaleup.md: fetch and sketch a chunk of a manifest, a few fetches
# at a time INSIDE one job, reads never kept. Array task i takes rows i*CHUNK .. i*CHUNK+CHUNK-1 of the manifest.
#
#   export POPSEQ_ROOT=<PopSeq checkout> GMG_ROOT=<GmG checkout> POPSKETCH=<popsketch binary> WORK=<work dir> \
#          MANIFEST=<manifest> CHUNK=32 PAR=4 CONCURRENCY=<PAR x the array limit>
#   N=$(bash -c '. scripts/popsra_chunk.sh; chunk_count "$MANIFEST" "$CHUNK"')
#   sbatch --export=ALL --array=0-$((N-1))%5 -o "$WORK/logs/%x.%A_%a.log" scripts/popsra_fetch_chunks.slurm
#
# Why chunks with PAR fetches inside each job: PopSeq's concurrency sweep (issue #27, spec 002) found the speed per fetch
# depends on how many fetches share a NODE (about 17 MB/s at 2 to 5 per node, about 5 MB/s at 20 or more), and the array
# limit %N does not control that. One job holds one node's worth of fetches: PAR (default 4). The cluster-wide number of
# fetches is PAR x the %N limit, and that product is CONCURRENCY, the value recorded in downloads.csv; fetch does not enforce it.
# Sizing: CHUNK = 32 runs at PAR = 4 is 8 rounds of about 4 to 8 minutes for a 1 Gb run taken by stride (phase 1 pilot, 3
# runs), so about 30 to 70 minutes of work per job, in line with AGENTS section 5 (1 to 1.5 h per unit). The memory (a
# 20 GB request for 4 fetches, from a peak of 3.1 GB per fetch at stride, n = 3) and the time (4 h, a ceiling for ENA-routed
# runs that take 10 to 30 minutes each) are requests, not measurements at PAR = 4: read MaxRSS and Elapsed from sacct.
#
# The job exits 0 when some fetches fail (their outcomes are results), and 1 when every fetch of the chunk exited 126 or 127
# (the command is not there: nothing was tried) or a module did not load (exit 3). A bad manifest exits 2.
# Each task writes its TSV files (results, profile, exits, hook_times, sampling_check) to $WORK/tasks/<i>/, so that tasks on
# different nodes never append to one file. After the run: `scripts/popsra_merge.py $WORK` joins them into $WORK/merged/, and
#   scripts/popsra_retry_list.py MANIFEST $WORK/merged/results.tsv -o retry.csv --exit-status
# turns the results into a retry manifest by the exit-code policy (exit 10 while runs remain). profile.tsv needs the PopSeq
# checkout at POPSEQ_ROOT to include PopSeq #45.
#SBATCH -J popsrachunk
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem 20G
#SBATCH -t 04:00:00
set -uo pipefail
: "${POPSEQ_ROOT:?}" "${GMG_ROOT:?}" "${POPSKETCH:?}" "${WORK:?}" "${MANIFEST:?}" "${SCRATCH:?run under SLURM}"
: "${CONCURRENCY:?set CONCURRENCY to PAR x the array limit (what popseq fetch records)}"
CHUNK=${CHUNK:-32}; PAR=${PAR:-4}; COVERAGE=${COVERAGE:-25}; SAMPLING=${SAMPLING:-stride}
i=${SLURM_ARRAY_TASK_ID:?run as an array task}
THREADS=$(( ${SLURM_CPUS_PER_TASK:-$PAR} / PAR )); [ "$THREADS" -ge 1 ] || THREADS=1

# shellcheck source=popsra_chunk.sh
. "$GMG_ROOT/scripts/popsra_chunk.sh"
runs=$(chunk_runs "$MANIFEST" "$CHUNK" "$i") || { echo "chunk $i: cannot read the manifest $MANIFEST" >&2; exit 2; }
[ -n "$runs" ] || { echo "chunk $i is past the end of the manifest: nothing to do"; exit 0; }

if [ "${SKIP_MODULES:-0}" != 1 ]; then
  # a failed module load would show up later as 127 on every fetch: stop here, with the reason
  module load sratoolkit/3.2.0 || { echo "chunk $i: module load sratoolkit/3.2.0 failed" >&2; exit 3; }
  module load parallel-fastq-dump/0.6.7 || { echo "chunk $i: module load parallel-fastq-dump/0.6.7 failed" >&2; exit 3; }
fi
echo "chunk $i: $(printf '%s\n' "$runs" | wc -l | tr -d ' ') runs, $PAR at a time, $THREADS threads each, coverage $COVERAGE ($SAMPLING), declared concurrency $CONCURRENCY, on $(hostname)"

# The TSV files of this task are its own ($WORK/tasks/<i>/): array tasks run on different nodes, and appends from several
# nodes to one file on a shared file system can interleave or lose lines. popsra_merge.sh joins them afterwards.
TASKDIR="$WORK/tasks/$i"
mkdir -p "$WORK/sketches" "$WORK/reads" "$WORK/logs/runs" "$TASKDIR"
export SKETCH_DIR="$WORK/sketches" HOOK_LOG="$TASKDIR/hook_times.tsv" CHECK_LOG="$TASKDIR/sampling_check.tsv" POPSKETCH
export WORK TASKDIR GMG_ROOT COVERAGE SAMPLING THREADS CONCURRENCY
FETCH_CMD=${FETCH_CMD:-"$HOME/.local/bin/uv run --project $POPSEQ_ROOT popseq fetch"}
export FETCH_CMD

one() { # RUN: one fetch, its own logs, its exit code on one line
  local run=$1 rc
  # shellcheck disable=SC2086
  $FETCH_CMD "$run" --tier coverage --coverage "$COVERAGE" --sampling "$SAMPLING" --out "$WORK/reads" \
      --exec "$GMG_ROOT/scripts/popsra_sketch_hook.sh" --threads "$THREADS" --require-pairs \
      --concurrency "$CONCURRENCY" --results "$TASKDIR/results.tsv" --profile "$TASKDIR/profile.tsv" </dev/null > "$WORK/logs/runs/$run.out" 2> "$WORK/logs/runs/$run.err"
  rc=$?
  printf "%s\t%s\t%s\n" "$run" "$rc" "$(date -u +%FT%TZ)" >> "$TASKDIR/exits.tsv"
  return 0
}
export -f one
# lines of exits.tsv at the start: an earlier attempt of this array task may have written some
exits_before=$(wc -l < "$TASKDIR/exits.tsv" 2>/dev/null || echo 0)
# the run accession is an argument of the function, never text pasted into shell code
printf '%s\n' "$runs" | xargs -P "$PAR" -I{} bash -c 'one "$1"' _ {}
echo "chunk $i: finished; $(wc -l < "$TASKDIR/exits.tsv" 2>/dev/null || echo 0) lines in $TASKDIR/exits.tsv"
# Failing fetches are results (exit codes 10 to 14, recorded), and the job still exits 0. But if EVERY fetch of the chunk
# exited 126 or 127 (the command was not found or could not run: no uv, a broken checkout), nothing was tried at all:
# that is a failed job, not a set of failed runs.
n_runs=$(printf '%s\n' "$runs" | wc -l | tr -d ' ')
n_notfound=$(tail -n +"$((exits_before + 1))" "$TASKDIR/exits.tsv" 2>/dev/null | awk -F'\t' '$2 == 126 || $2 == 127 {n++} END {print n + 0}')
if [ "$n_notfound" -ge "$n_runs" ]; then
  echo "chunk $i: all $n_runs fetches exited 126 or 127 (command not found: $FETCH_CMD): see $WORK/logs/runs/*.err" >&2
  exit 1
fi
exit 0
