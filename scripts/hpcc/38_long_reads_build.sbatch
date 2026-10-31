#!/bin/bash -l
#SBATCH -J lr_build
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/38_lr_build.%j.log
# Long reads (#44): a PRIVATE build of the branch with the `high_error_reads` flag, in $GR/long_reads/bin/gmg (the shared $GR/bin/gmg is not touched).
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
BRANCH=${BRANCH:-long-reads}
R=$GR/long_reads/repo
[ -d $R ] || git clone -q --branch $BRANCH https://github.com/stajichlab/GmG.git $R
cd $R && git fetch -q origin $BRANCH && git checkout -q -B $BRANCH origin/$BRANCH
export PATH=$HOME/.cargo/bin:$PATH
cargo build --release -q -p gmg-cli
mkdir -p $GR/long_reads/bin
cp target/release/gmg $GR/long_reads/bin/gmg.new && mv $GR/long_reads/bin/gmg.new $GR/long_reads/bin/gmg
$GR/long_reads/bin/gmg --version; git log --oneline | head -1
