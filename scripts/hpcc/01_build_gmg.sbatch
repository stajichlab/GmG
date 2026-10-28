#!/bin/bash -l
#SBATCH -J gmg_build
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/01_build_gmg.%j.log
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
cd $GMG_REPO
pixi run test 2>&1 | grep 'test result'
pixi run build
bash scripts/hpcc/deploy_gmg.sh
$GMG --version
