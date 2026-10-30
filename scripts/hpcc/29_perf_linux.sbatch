#!/bin/bash -l
#SBATCH -J gmg_perf
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/29_perf.%j.log
# #49: the performance table on a cluster node (compare with the macOS numbers). The golden panel must be cached by 28_parity_linux.sbatch.
# usage: sbatch 29_perf_linux.sbatch SCRIPT_DIR [VER]
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
SD=$1; VER=${2:-2026.10.2}
echo "node: $(hostname), $(nproc) cpus; binary: $(readlink -f $GMG)"
bash $SD/perf_native.sh $GMG $GR/release/$VER/bfd-$VER-s10000.gmgi $GR/release/$VER/bfd-reps-$VER-s20000.gmgi $GR/parity/cache
