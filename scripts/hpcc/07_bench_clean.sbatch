#!/bin/bash -l
#SBATCH -J gmg_bclean
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/07_bench_clean.%j.log
# leave-out benchmarks with the suspected mislabeled genomes (reports/label_audit/exclude_ids.txt) removed as references;
# the same ids are dropped from the queries in the analysis (GMG_EXCLUDE).
set -uo pipefail
source /bigdata/stajichlab/shared/projects/GmG_runs/src/GmG/scripts/hpcc/env.sh
mkdir -p $GR/bench/clean
cd $GR/bench/clean
for S in 10000 2000; do
  N=0; [ $S = 2000 ] && N=6000
  args=""; [ "$N" != "0" ] && args="--sample $N"
  for M in loo los log; do
    out=bench_${M}_s${S}.tsv
    [ -s $out ] || $GMG --threads 4 bench --index $GR/index/bfd_s${S}.gmgi --sketches $GR/sketch --ref-list $GR/index/clean_refs.txt \
        --mode $M $args --seed 1 -o $out
  done
  GMG_EXCLUDE=$GR/index/exclude_ids.txt python3 $GR/src/GmG/scripts/analyze_bench.py summary_s${S} bench_*_s${S}.tsv > /dev/null
done
ls -la
