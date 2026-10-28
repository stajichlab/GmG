#!/bin/bash -l
#SBATCH -J gmg_test
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/00_test.%j.log
# Test the checkout in ~/projects/GmG on a compute node (the login node kills the heavy test build): unit + oracle tests,
# the wasm/native parity test and the release build. Use after `git pull`.
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
cd $GMG_REPO
git log -1 --oneline
export CARGO_BUILD_JOBS=$SLURM_CPUS_PER_TASK
# keep the output (tail, not grep) so failures stay visible; pipefail makes a failing pixi fail the job
pixi run test 2>&1 | tail -40
pixi run wasm-test 2>&1 | tail -3
pixi run build 2>&1 | tail -2
ls -la target/release/gmg
