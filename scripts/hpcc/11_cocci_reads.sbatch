#!/bin/bash -l
#SBATCH -J gmg_cocci
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/11_cocci.%j.log
# Real Illumina reads (Coccidioides, species-labelled): classify R1 of every read set from the first 100k / 500k / 2M reads.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/reads/cocci_validation; FARM=$R/farm
mkdir -p $FARM
while read -r p; do ln -sf "$p" $FARM/$(basename "$p"); done < $GMG_REPO/data/validation/cocci_reads.r1.list
echo "farm: $(ls $FARM | wc -l) read files"
for N in 100000 500000 2000000; do
  [ -s $R/n$N/results.tsv ] && continue
  /usr/bin/time -f "n=$N wall=%e s maxrss=%M KB" $GMG --threads ${SLURM_CPUS_PER_TASK} run $FARM --index $GR/index/bfd_s10000.v1.gmgi -o $R/n$N --no-report --max-reads $N
done
