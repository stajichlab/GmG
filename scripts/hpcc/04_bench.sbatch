#!/bin/bash -l
#SBATCH -J gmg_bench
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/04_bench.%j.log
# leave-one-out (loo), leave-species-out (los), leave-genus-out (log) on BFD genomes.
# scaled=10000: all genomes; scaled=2000 and 1000: a seeded sample (cost grows with index size).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
mkdir -p $GR/bench
run() { # scaled sample(0=all)
  local S=$1 N=$2 args=""
  [ "$N" != "0" ] && args="--sample $N"
  for M in loo los log; do
    out=$GR/bench/bench_${M}_s${S}.tsv
    [ -s $out ] && continue
    /usr/bin/time -f "bench $M s=$S wall=%e s maxrss=%M KB" $GMG --threads 4 bench --index $GR/index/bfd_s${S}.gmgi --sketches $GR/sketch \
       --mode $M $args --seed 1 -o $out
  done
}
run 10000 0
run 2000 6000
run 1000 3000
cd $GR/bench
for S in 10000 2000 1000; do
  ls bench_*_s${S}.tsv >/dev/null 2>&1 && python3 $GMG_REPO/scripts/analyze_bench.py summary_s${S} bench_*_s${S}.tsv > /dev/null
done
ls -la $GR/bench
