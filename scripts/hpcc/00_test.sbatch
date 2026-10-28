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
git log --oneline | head -1
export CARGO_BUILD_JOBS=$SLURM_CPUS_PER_TASK
pixi run test 2>&1 | grep -E "test result|FAILED|error"
pixi run wasm-test 2>&1 | tail -1
pixi run build 2>&1 | tail -1
