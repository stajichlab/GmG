#!/bin/bash -l
#SBATCH -J mixdet_build
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_mixdet_build.%j.log
# #31: build the mixture-detector branch privately (GmG_runs/mixdet/bin/gmg); the shared GmG_runs/bin/gmg is not touched.
# The source tree is copied to GmG_runs/mixdet/src first (rsync from the working copy, or a git checkout of the branch).
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/mixdet; cd $R/src
export CARGO_TARGET_DIR=$R/target
cargo build --release -p gmg-cli
cp $CARGO_TARGET_DIR/release/gmg $R/bin/gmg-${TAG:-dev}
ln -sfn gmg-${TAG:-dev} $R/bin/gmg
$R/bin/gmg --version
