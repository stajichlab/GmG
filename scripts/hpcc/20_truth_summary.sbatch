#!/bin/bash -l
#SBATCH -J gmg_truth_sum
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/20_truth_summary.%j.log
# #39: list the NCBI RefSeq fungal assemblies (type material and reference flags) for the trusted truth set.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
T=$GR/truthset; mkdir -p $T; cd $T
datasets summary genome taxon fungi --assembly-source RefSeq --as-json-lines > refseq_fungi.jsonl
wc -l refseq_fungi.jsonl
head -1 refseq_fungi.jsonl | python3 -m json.tool | head -80
