#!/bin/bash -l
#SBATCH -J gmg_margin_ext
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/36_margin_ext.%j.log
# #58: the IndigoAg assemblies (index 2026.10.3) with the margin thresholds (THF=thresholds.json), next to relstab/indigo_2026.10.3
# (the default thresholds before the margin). Same binary; only the thresholds differ.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
THF=${THF:?}
V=${V:-2026.10.3}
$GMG --threads ${SLURM_CPUS_PER_TASK} -q run $GR/external/indigo/in --index $GR/release/$V/bfd-$V-s10000.gmgi -o $GR/relstab/indigo_$V.m58 --no-report --thresholds $THF
