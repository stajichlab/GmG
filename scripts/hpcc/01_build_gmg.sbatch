#!/bin/bash -l
#SBATCH -J gmg_build
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/01_build_gmg.%j.log
set -euo pipefail
source /bigdata/stajichlab/shared/projects/GmG_runs/src/GmG/scripts/hpcc/env.sh
cd $GR/src/GmG
cargo --version
cargo test --release -p gmg-core --lib 2>&1 | tail -3
cargo build --release -p gmg-cli
mkdir -p $GR/bin && cp target/release/gmg $GR/bin/gmg
$GMG --version
