#!/bin/bash -l
#SBATCH -J gmg_absent
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/14_absent.%j.log
# Reads calibration (#32): index WITHOUT Coccidioides immitis, then classify the Coccidioides read sets.
# immitis reads must not get an immitis call (absent species); the closest relative is C. posadasii, so false species calls are visible.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/reads/cocci_absent; mkdir -p $R
grep -v "Coccidioides immitis" $GR/index/meta.tsv > $R/meta_no_immitis.tsv
echo "meta: $(wc -l < $GR/index/meta.tsv) -> $(wc -l < $R/meta_no_immitis.tsv)"
[ -s $R/no_immitis.gmgi ] || $GMG --threads 4 build --meta $R/meta_no_immitis.tsv --sketches $GR/sketch --scaled 10000 -o $R/no_immitis.gmgi
for N in 500000 2000000; do
  [ -s $R/n$N/results.tsv ] && continue
  $GMG --threads ${SLURM_CPUS_PER_TASK} run $GR/reads/cocci_validation/farm --index $R/no_immitis.gmgi -o $R/n$N --no-report --max-reads $N
done
