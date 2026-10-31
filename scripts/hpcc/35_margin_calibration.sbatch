#!/bin/bash -l
#SBATCH -J gmg_margin
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_margin.%j.log
# #58: derive the genus/family safety margin with the clade_calibration protocol (scripts/margin_calibration.py); needs pandas.
# usage: REPO=<dir with scripts/> SHIPPED=<known.json before the margin> sbatch 35_margin_calibration.sbatch
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REPO=${REPO:-$GR/margin/repo}
mkdir -p $GR/margin
python3 $REPO/scripts/margin_calibration.py $GR/bench/clean $GR/index/exclude_ids.txt $GR/margin/margin_cal 10 0.95 0.02 \
  --shipped-json ${SHIPPED:-$GR/margin/known_before.json}
