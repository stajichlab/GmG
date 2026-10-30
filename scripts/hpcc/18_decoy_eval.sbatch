#!/bin/bash -l
#SBATCH -J gmg_decoy_eval
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/18_decoy_eval.%j.log
# #79 evaluation: (1) held-out non-fungal genomes (not decoys) against the index with decoys, and against the one without;
# (2) the decoy genomes themselves. usage: sbatch 18_decoy_eval.sbatch [VER]   (index: $GR/release/VER/bfd-VER-s10000.gmgi)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
VER=${1:-2026.10.2}
NEW=$GR/release/$VER/bfd-$VER-s10000.gmgi
OLD=$GR/release/2026.10.1/bfd-2026.10.1-s10000.gmgi
R=$GR/decoys/eval; mkdir -p $R
for w in old:$OLD new:$NEW; do
  n=${w%%:*}; ix=${w#*:}
  $GMG --threads ${SLURM_CPUS_PER_TASK} run $GR/decoys/held --index $ix -o $R/held_$n --no-report -q
done
# the decoys against the index that contains them: sanity check (must be called as themselves and flagged)
$GMG --threads ${SLURM_CPUS_PER_TASK} run $GR/decoys/dl --index $NEW -o $R/decoys_self --no-report -q
ls -la $R
