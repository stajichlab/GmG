#!/bin/bash
# Run on the HPCC after `cargo build --release -p gmg-cli`. Never overwrite a binary that a job may be running
# (that crashed running jobs with SIGSEGV/SIGBUS on 2026-10-28): install a versioned copy and swap the symlink.
set -eu
GR=/bigdata/stajichlab/shared/projects/GmG_runs
src=$GR/src/GmG/target/release/gmg
v=$(date +%Y%m%d-%H%M%S)
cp $src $GR/bin/gmg-$v
ln -sfn gmg-$v $GR/bin/gmg.new && mv -Tf $GR/bin/gmg.new $GR/bin/gmg
echo "deployed gmg-$v"
