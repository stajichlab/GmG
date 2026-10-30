#!/bin/bash -l
#SBATCH -J gmg_extset
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/20_extset.%j.log
# External validation set (#40): classify $GR/external/$SET/in (symlinks made by scripts/prep_external.py) with the decoy index,
# score against the provider names in labels.tsv when there are any. Per-genome output stays on the HPCC.
# usage: SET=kasson sbatch 20_external_set.sbatch      (SET: kasson antarctica ccfee_aaftf ccfee_shovill zygolife)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
: "${SET:?SET}"
INDEX=${INDEX:-$GR/release/2026.10.2/bfd-2026.10.2-s10000.gmgi}
R=$GR/external/$SET; cd $R
echo "set $SET: $(ls in | wc -l) inputs; index $INDEX; $($GMG --version)"
rm -rf out
$GMG --threads ${SLURM_CPUS_PER_TASK} run in --index $INDEX -o out --no-report
cp out/results.tsv results.tsv
if [ "$(awk -F'\t' 'NR>1 && $2!=""' labels.tsv | wc -l)" -gt 0 ]; then
  python3 $GR/external/scripts/score_by_name.py results.tsv labels.tsv $GR/taxdump_20261029 score
fi
