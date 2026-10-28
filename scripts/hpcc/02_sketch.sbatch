#!/bin/bash -l
#SBATCH -J gmg_sketch
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 01:30:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/02_sketch.%A_%a.log
# usage: sbatch --array=0-119%50 02_sketch.sbatch   (chunks are $GR/index/chunks/chunk_NNN; re-runnable, existing sketches are kept)
set -euo pipefail
source /bigdata/stajichlab/shared/projects/GmG_runs/src/GmG/scripts/hpcc/env.sh
CH=$(printf "%s/index/chunks/chunk_%03d" $GR ${SLURM_ARRAY_TASK_ID})
mkdir -p $GR/sketch
$GMG --threads ${SLURM_CPUS_PER_TASK} sketch -k 21 -s 1000 -o $GR/sketch --list $CH
