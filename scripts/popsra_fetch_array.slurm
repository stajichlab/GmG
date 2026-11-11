#!/bin/bash
# Phase 1 of docs/plan_popsra_scaleup.md: fetch one run at a coverage target and sketch it, one run per array
# task, reads never kept. Task i takes row i of the manifest (a CSV or TSV with a run_accession column).
#
#   export POPSEQ_ROOT=<PopSeq checkout> GMG_ROOT=<GmG checkout> POPSKETCH=<popsketch binary> \
#          WORK=<work dir> MANIFEST=<manifest> COVERAGE=25 CONCURRENCY=5
#   sbatch --export=ALL --array=0-9%5 -o "$WORK/logs/%x.%A_%a.log" scripts/popsra_fetch_array.slurm
#
# --require-pairs is always on: the PopSeq panel is selected for paired isolates, and a run that NCBI holds as
# single-end (spots_with_mates 0) can come back from fastq-dump as mate files that cannot be verified. With the
# flag such a run is taken from ENA as verified pairs instead (found in the phase 1 pilot: ERR2299868).
#
# CONCURRENCY must equal the array's %N throttle: it is the value recorded in benchmarks/downloads.csv, and
# `popseq fetch` does not enforce it. Outputs in $WORK: sketches/ (one .gsk, two halves and a .json per run),
# results.tsv (popseq fetch result lines), hook_times.tsv (sketch seconds and bytes), reads/ (sidecars only).
# Resources: 4 cores for parallel-fastq-dump; 8 GB is a request, not a measurement (peak memory of the fetch
# plus sketch step is UNMEASURED; read MaxRSS from sacct after the run). Time: short partition, 2 h.
#SBATCH -J popsra
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 02:00:00
set -uo pipefail
: "${POPSEQ_ROOT:?}" "${GMG_ROOT:?}" "${POPSKETCH:?}" "${WORK:?}" "${MANIFEST:?}" "${SCRATCH:?run under SLURM}"
COVERAGE=${COVERAGE:-25}; CONCURRENCY=${CONCURRENCY:-1}
i=${SLURM_ARRAY_TASK_ID:?run as an array task}

module load sratoolkit/3.2.0
module load parallel-fastq-dump/0.6.7

# Row i of the manifest (header skipped); the delimiter is a tab or a comma, as in popseq fetch.
delim=","; head -1 "$MANIFEST" | grep -q $'\t' && delim=$'\t'
col=$(head -1 "$MANIFEST" | tr "$delim" '\n' | grep -n -x run_accession | cut -d: -f1)
RUN=$(tail -n +2 "$MANIFEST" | sed -n "$((i + 1))p" | cut -d"$delim" -f"$col")
[ -n "$RUN" ] || { echo "task $i: no manifest row" >&2; exit 2; }
echo "task $i: $RUN on $(hostname), coverage $COVERAGE, declared concurrency $CONCURRENCY"

mkdir -p "$WORK/sketches" "$WORK/reads"
export SKETCH_DIR="$WORK/sketches" HOOK_LOG="$WORK/hook_times.tsv" POPSKETCH
"$HOME/.local/bin/uv" run --project "$POPSEQ_ROOT" popseq fetch "$RUN" --tier coverage --coverage "$COVERAGE" \
    --out "$WORK/reads" --exec "$GMG_ROOT/scripts/popsra_sketch_hook.sh" --threads "$SLURM_CPUS_PER_TASK" \
    --require-pairs --concurrency "$CONCURRENCY" --results "$WORK/results.tsv" </dev/null
rc=$?
echo "task $i: $RUN popseq fetch exit $rc"
exit $rc
