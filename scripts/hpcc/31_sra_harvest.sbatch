#!/bin/bash -l
#SBATCH -J gmg_sraget
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/31_sra_harvest.%j.log
# Real SRA/ENA reads for #42, step 2: stream the first 500,000 read pairs of every run in the accession list (never a whole run; <= 24 GB in total).
# Raw reads go to a scratch directory under GmG_runs and are deleted after scoring (32_sra_titration.sbatch).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${READS_SCRIPTS:-$GMG_REPO/scripts}
R=$GR/reads_sra; mkdir -p $R/raw; cd $R
python3 $S/sra_harvest.py sra_runs.tsv raw ${PAIRS:-500000} 3 24e9
cat raw/bytes.tsv | awk -F'\t' 'NR>1{b+=$2+$3; n++} END{printf "%d runs, %.2f GB of compressed bytes read\n", n, b/1e9}'
