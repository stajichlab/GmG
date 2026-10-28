#!/bin/bash -l
#SBATCH -J gmg_reps
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:30:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/06_reps.%j.log
# How accurate and how small is a representatives-only reference?
#   reps1 = best genome per species (BFD rule), reps3 = best 3 per species.
#   Queries = genomes NOT in the reference set (leave-one-out semantics, all of them).
set -uo pipefail
source /bigdata/stajichlab/shared/projects/GmG_runs/src/GmG/scripts/hpcc/env.sh
cd $GR/index
for S in 20000 50000; do
  [ -s bfd_s${S}.gmgi ] || $GMG --threads 4 build --meta meta.tsv --sketches $GR/sketch --scaled $S -o bfd_s${S}.gmgi
done
mkdir -p $GR/bench/reps
for K in 1 3; do
  for S in 10000 20000 50000; do
    # size of the standalone representatives-only index
    $GMG --threads 4 build --meta reps.meta_reps${K}.tsv --sketches $GR/sketch --scaled $S -o reps${K}_s${S}.gmgi 2>&1 | tail -1
    $GMG convert reps${K}_s${S}.gmgi -o reps${K}_s${S}.v1.gmgi 2>&1 | tail -1
    # accuracy of the same reference set, queried with the genomes left out of it
    $GMG --threads 4 bench --index bfd_s${S}.gmgi --sketches $GR/sketch --ref-list reps.reps${K}.txt --queries-outside \
        --mode loo -o $GR/bench/reps/bench_loo_reps${K}_s${S}.tsv
  done
done
cd $GR/bench/reps
for K in 1 3; do for S in 10000 20000 50000; do
  python3 $GR/src/GmG/scripts/analyze_bench.py summary_reps${K}_s${S} bench_loo_reps${K}_s${S}.tsv > /dev/null
done; done
ls -la $GR/index/reps*.gmgi
