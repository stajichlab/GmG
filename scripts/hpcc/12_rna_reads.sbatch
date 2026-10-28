#!/bin/bash -l
#SBATCH -J gmg_rna
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/12_rna.%j.log
# BFD RNA-seq read sets (species-labelled): first 100k / 500k / 2M reads, classified against the DNA genome index.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/reads/rna_validation; FARM=$R/farm
mkdir -p $FARM
while read -r p; do ln -sf "$p" $FARM/$(basename "$p"); done < $R/rna.list
echo "farm: $(ls $FARM | wc -l) read files"
for N in 100000 500000 2000000; do
  [ -s $R/n$N/results.tsv ] && continue
  /usr/bin/time -f "n=$N wall=%e s maxrss=%M KB" $GMG --threads ${SLURM_CPUS_PER_TASK} run $FARM --index $GR/index/bfd_s10000.v1.gmgi -o $R/n$N --no-report --max-reads $N
done
