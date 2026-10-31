#!/bin/bash -l
#SBATCH -J gmg_bench_i5
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/bench_i5.%A_%a.log
# leave-out bench on the 2026.10.5 reference set (23,597 fungal genomes; no decoys), same method as 07_bench_clean.sbatch
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
mkdir -p $GR/bench/i5; cd $GR/bench/i5
MODES=(loo los log); M=${MODES[$SLURM_ARRAY_TASK_ID]}
$GMG --threads 4 bench --index $GR/index/bfd_s10000.gmgi --sketches $GR/sketch --ref-list $GR/index/refs_2026.10.5.txt --mode $M --seed 1 -o bench_${M}_s10000.tsv
