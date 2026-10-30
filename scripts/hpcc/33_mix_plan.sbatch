#!/bin/bash -l
#SBATCH -J mix_plan
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/33_mix_plan.%j.log
# #45: choose the pairs and triples of real read sets to mix (scripts/mix_plan.py; counting read pairs with zcat is the slow part).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${MIX_SCRIPTS:-$GR/reads/mix/scripts}
R=$GR/reads/mix; mkdir -p $R/plan
python3 $S/mix_plan.py $GR/reads/titration $GMG_REPO/data/validation/cocci_reads.truth.tsv \
  /bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/input $R/plan
cut -f1,2,3 $R/plan/plan.tsv | column -t | head -50
cut -f1,2 $R/plan/triples.tsv | column -t
