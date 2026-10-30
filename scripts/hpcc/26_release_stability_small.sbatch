#!/bin/bash -l
#SBATCH -J gmg_relstab2
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/26_relstab_small.%j.log
# #48: the IndigoAg assemblies and the held-out non-fungal genomes against every release (same binary).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
VERS=${VERS:-"2026.10 2026.10.1 2026.10.2 2026.10.3"}
R=$GR/relstab; mkdir -p $R
for v in $VERS; do
  I=$GR/release/$v/bfd-$v-s10000.gmgi
  [ -s $I ] || { echo "no index for $v"; continue; }
  $GMG --threads ${SLURM_CPUS_PER_TASK} -q run $GR/external/indigo/in --index $I -o $R/indigo_$v --no-report
  $GMG --threads ${SLURM_CPUS_PER_TASK} -q run $GR/decoys/held --index $I -o $R/held_$v --no-report
done
