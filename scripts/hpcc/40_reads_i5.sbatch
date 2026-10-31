#!/bin/bash -l
#SBATCH -J gmg_reads_i5
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/reads_i5.%A_%a.log
# Cocci and RNA-seq read sets (as 11_cocci_reads / 12_rna_reads) against the 2026.10.5 index. usage: sbatch --array=0-5 reads_i5.sbatch
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
NS=(100000 500000 2000000); N=${NS[$((SLURM_ARRAY_TASK_ID % 3))]}
if [ $SLURM_ARRAY_TASK_ID -lt 3 ]; then R=$GR/reads/cocci_validation; else R=$GR/reads/rna_validation; fi
$GMG --threads 8 -q run $R/farm --index $GR/release/2026.10.5/bfd-2026.10.5-s10000.gmgi -o $R/n$N.i5 --no-report --max-reads $N
