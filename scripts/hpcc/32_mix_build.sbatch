#!/bin/bash -l
#SBATCH -J mix_build
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/32_mix_build.%j.log
# #45: build current main privately (GmG_runs/reads/mix/bin/gmg) so the mixed-read runs have every flag of main (weak_phylum, ...)
# without swapping the shared GmG_runs/bin/gmg that other jobs are using.
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/reads/mix; mkdir -p $R/src $R/bin; cd $R/src
[ -d repo ] || git clone -q https://github.com/stajichlab/GmG.git repo
cd repo
git fetch -q origin
git reset -q --hard origin/main
pixi run build
h=$(git rev-parse --short HEAD)
cp target/release/gmg $R/bin/gmg-$h
ln -sfn gmg-$h $R/bin/gmg
$R/bin/gmg --version
